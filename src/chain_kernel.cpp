#include "panvar/chain_kernel.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>

namespace panvar {

namespace {

// The per-edge shift that keeps exp(log psi) representable. Constant across the edge, so it cancels
// in the per-block normalisation; accumulated into log_scale so the partition stays recoverable.
double linkage_shift(const ChainEdgeLinkage& e) {
    double m = -std::numeric_limits<double>::infinity();
    for (double v : e.log_psi) m = std::max(m, v);
    return std::isfinite(m) ? m : 0.0;
}

inline double psi_weight(const ChainEdgeLinkage& e, double shift,
                         std::size_t i, std::size_t j, std::size_t i2, std::size_t j2) {
    // Homologue 1 takes (allele of i at the LEFT block, allele of i2 at the RIGHT).
    // Homologue 2 takes (allele of j at the LEFT block, allele of j2 at the RIGHT).
    // In that order, with no implicit swap.
    const std::size_t a1 = e.allele_a[i],  b1 = e.allele_b[i2];
    const std::size_t a2 = e.allele_a[j],  b2 = e.allele_b[j2];
    // RANGE-CHECKED AT THE POINT OF USE. build_allele_mapping() already refuses anything out of
    // range, but this is where a bad index would corrupt memory rather than produce a wrong number,
    // so the check is repeated here and a violation degrades to a neutral 1.0 instead of reading
    // past the table.
    if (a1 >= e.n_a || a2 >= e.n_a || b1 >= e.n_b || b2 >= e.n_b) return 1.0;
    const std::size_t idx = ((a1 * e.n_b + b1) * e.n_a + a2) * e.n_b + b2;
    if (idx >= e.log_psi.size()) return 1.0;
    return std::exp(e.log_psi[idx] - shift);
}

}  // namespace

namespace {

// One grouped correction pass. `src` is the incoming message, `dst` the factorised baseline already
// computed; the corrections are ADDED to dst. `forward` selects which endpoint's alleles index the
// source and which the target: forward goes A -> B, backward B -> A, and getting that backwards
// would silently apply the transpose.
void apply_grouped_corrections(const SparseEdgeLinkage& e, std::size_t nh, double r,
                               const std::vector<double>& src, std::vector<double>& dst,
                               bool forward, ChainKernelStats& st) {
    const double one_r = 1.0 - r, rn = r / static_cast<double>(nh);
    const std::vector<std::uint32_t>& sa = forward ? e.allele_a : e.allele_b;
    const std::vector<std::vector<std::uint32_t>>& gs = forward ? e.group_a : e.group_b;
    const std::vector<std::vector<std::uint32_t>>& gt = forward ? e.group_b : e.group_a;
    const std::size_t n_src = gs.size();
    // GROUP ROW AND COLUMN SUMS over the SOURCE allele partition. Each haplotype belongs to exactly
    // one group, so all of these together cost O(n_h^2) -- the same order as the factorised step --
    // rather than being recomputed per class.
    std::vector<std::vector<double>> rs(n_src, std::vector<double>(nh, 0.0));
    std::vector<std::vector<double>> cs(n_src, std::vector<double>(nh, 0.0));
    for (std::size_t i = 0; i < nh; ++i) {
        for (std::size_t j = 0; j < nh; ++j) {
            const double v = src[i * nh + j];
            if (v == 0.0) continue;
            rs[sa[j]][i] += v;   // sum over j in group a2, indexed by i
            cs[sa[i]][j] += v;   // sum over i in group a1, indexed by j
        }
    }
    const auto total = [&](std::uint32_t a1, std::uint32_t a2) {
        double t = 0.0;
        for (std::uint32_t i : gs[a1]) t += rs[a2][i];
        return t;
    };
    const auto add = [&](std::uint32_t a1, std::uint32_t b1, std::uint32_t a2, std::uint32_t b2,
                         double cm1) {
        if (cm1 == 0.0) return;
        if (a1 >= n_src || a2 >= n_src || b1 >= gt.size() || b2 >= gt.size()) return;
        const double tot = total(a1, a2);
        for (std::uint32_t i2 : gt[b1]) {
            const bool i2_in_a1 = sa[i2] == a1;
            for (std::uint32_t j2 : gt[b2]) {
                // The four terms of t(i->i2) t(j->j2) with t = (1-r)I + (r/n) 11^T.
                double w = rn * rn * tot;
                if (i2_in_a1) w += one_r * rn * rs[a2][i2];
                if (sa[j2] == a2) w += rn * one_r * cs[a1][j2];
                if (i2_in_a1 && sa[j2] == a2) w += one_r * one_r * src[i2 * nh + j2];
                dst[i2 * nh + j2] += cm1 * w;
                ++st.corrections_applied;
            }
        }
    };
    for (const SparsePhaseClass& c : e.classes) {
        // Both homologue orders of each biological phase, applied exactly once each.
        add(c.amin, c.bmin, c.amax, c.bmax, c.straight_m1);
        add(c.amax, c.bmax, c.amin, c.bmin, c.straight_m1);
        add(c.amin, c.bmax, c.amax, c.bmin, c.crossed_m1);
        add(c.amax, c.bmin, c.amin, c.bmax, c.crossed_m1);
    }
    // A correction can nearly cancel the baseline where psi approaches zero. Tiny negatives are
    // roundoff and are clamped with a count; a materially negative weight is an internal error.
    for (double& x : dst) {
        if (x >= 0.0) continue;
        st.worst_negative = std::min(st.worst_negative, x);
        ++st.clamped_negatives;
        x = 0.0;
    }
}

}  // namespace

void chain_forward_backward(
    std::size_t nh, std::size_t nb, double r,
    const std::function<void(std::size_t, std::vector<double>&)>& block_emissions,
    const std::vector<ChainEdgeLinkage>* edges,
    std::vector<std::vector<double>>& fwd,
    std::vector<std::vector<double>>& bwd,
    ChainKernelStats* stats,
    const std::vector<SparseEdgeLinkage>* sparse) {
    if (nh == 0 || nb == 0) return;
    const std::size_t ns = nh * nh;
    fwd.assign(nb, std::vector<double>(ns, 0.0));
    bwd.assign(nb, std::vector<double>(ns, 0.0));
    ChainKernelStats st;

    const auto sparse_at = [&](std::size_t b) -> const SparseEdgeLinkage* {
        if (sparse == nullptr || sparse->size() != nb) return nullptr;
        const SparseEdgeLinkage& e = (*sparse)[b];
        return e.has_corrections() ? &e : nullptr;
    };
    const auto edge_at = [&](std::size_t b) -> const ChainEdgeLinkage* {
        if (edges == nullptr || edges->size() != nb) return nullptr;
        const ChainEdgeLinkage& e = (*edges)[b];
        return e.active ? &e : nullptr;
    };
    auto normalize = [&](std::vector<double>& v) {
        const double s = std::accumulate(v.begin(), v.end(), 0.0);
        if (s > 0.0) {
            for (double& x : v) x /= s;
            st.log_weight_sum += std::log(s);
        }
    };

    std::vector<double> e(ns), beta(ns), rowsum(nh), colsum(nh);
    block_emissions(0, e);
    fwd[0] = e;
    normalize(fwd[0]);
    for (std::size_t bi = 1; bi < nb; ++bi) {
        const std::vector<double>& prev = fwd[bi - 1];
        const SparseEdgeLinkage* sp = sparse_at(bi);
        const ChainEdgeLinkage* lk = sp ? nullptr : edge_at(bi);
        block_emissions(bi, e);
        if (sp != nullptr) {
            // GROUPED PATH: the factorised baseline, then sparse corrections. The baseline is the
            // SAME arithmetic the linkage-free branch runs, so an edge with no classes reduces to
            // it exactly rather than to a correction loop that happens to add zero.
            ++st.grouped_edges;
            std::fill(rowsum.begin(), rowsum.end(), 0.0);
            for (std::size_t i = 0; i < nh; ++i)
                for (std::size_t j = 0; j < nh; ++j) rowsum[i] += prev[i * nh + j];
            for (std::size_t i = 0; i < nh; ++i)
                for (std::size_t j = 0; j < nh; ++j)
                    beta[i * nh + j] = (1.0 - r) * prev[i * nh + j] + (r / nh) * rowsum[i];
            std::fill(colsum.begin(), colsum.end(), 0.0);
            for (std::size_t i = 0; i < nh; ++i)
                for (std::size_t j = 0; j < nh; ++j) colsum[j] += beta[i * nh + j];
            std::vector<double> acc(ns, 0.0);
            for (std::size_t i = 0; i < nh; ++i)
                for (std::size_t j = 0; j < nh; ++j)
                    acc[i * nh + j] = (1.0 - r) * beta[i * nh + j] + (r / nh) * colsum[j];
            apply_grouped_corrections(*sp, nh, r, prev, acc, true, st);
            for (std::size_t k = 0; k < ns; ++k) fwd[bi][k] = e[k] * acc[k];
        } else if (lk == nullptr) {
            // FACTORISED PATH, preserved exactly as the legacy recurrence wrote it.
            ++st.factorised_edges;
            std::fill(rowsum.begin(), rowsum.end(), 0.0);
            for (std::size_t i = 0; i < nh; ++i)
                for (std::size_t j = 0; j < nh; ++j) rowsum[i] += prev[i * nh + j];
            for (std::size_t i = 0; i < nh; ++i)
                for (std::size_t j = 0; j < nh; ++j)
                    beta[i * nh + j] = (1.0 - r) * prev[i * nh + j] + (r / nh) * rowsum[i];
            std::fill(colsum.begin(), colsum.end(), 0.0);
            for (std::size_t i = 0; i < nh; ++i)
                for (std::size_t j = 0; j < nh; ++j) colsum[j] += beta[i * nh + j];
            for (std::size_t i = 0; i < nh; ++i)
                for (std::size_t j = 0; j < nh; ++j)
                    fwd[bi][i * nh + j] =
                        e[i * nh + j] * ((1.0 - r) * beta[i * nh + j] + (r / nh) * colsum[j]);
        } else {
            // LINKED PATH. psi is not separable, so the transition cannot be factorised and this is
            // O(n_h^4). The potential is used UNNORMALISED.
            ++st.linked_edges;
            const double shift = linkage_shift(*lk);
            st.log_weight_sum += shift;
            for (std::size_t i2 = 0; i2 < nh; ++i2)
                for (std::size_t j2 = 0; j2 < nh; ++j2) {
                    double acc = 0.0;
                    for (std::size_t i = 0; i < nh; ++i) {
                        const double ti = (i == i2 ? (1.0 - r) : 0.0) + r / static_cast<double>(nh);
                        for (std::size_t j = 0; j < nh; ++j) {
                            const double p = prev[i * nh + j];
                            if (p == 0.0) continue;
                            const double tj =
                                (j == j2 ? (1.0 - r) : 0.0) + r / static_cast<double>(nh);
                            acc += p * ti * tj * psi_weight(*lk, shift, i, j, i2, j2);
                        }
                    }
                    fwd[bi][i2 * nh + j2] = e[i2 * nh + j2] * acc;
                }
        }
        normalize(fwd[bi]);
    }

    std::fill(bwd[nb - 1].begin(), bwd[nb - 1].end(), 1.0);
    {
        const double s = static_cast<double>(ns);
        for (double& x : bwd[nb - 1]) x /= s;
    }
    for (std::size_t bi = nb - 1; bi-- > 0;) {
        block_emissions(bi + 1, e);
        const SparseEdgeLinkage* sp = sparse_at(bi + 1);
        const ChainEdgeLinkage* lk = sp ? nullptr : edge_at(bi + 1);
        std::vector<double> t(ns);
        for (std::size_t k = 0; k < ns; ++k) t[k] = bwd[bi + 1][k] * e[k];
        if (sp != nullptr) {
            // THE SAME DECOMPOSITION, independently, for the backward message. Reusing the forward
            // pass with the roles unswapped would apply the transpose of the intended correction.
            std::fill(rowsum.begin(), rowsum.end(), 0.0);
            for (std::size_t i = 0; i < nh; ++i)
                for (std::size_t j = 0; j < nh; ++j) rowsum[i] += t[i * nh + j];
            for (std::size_t i = 0; i < nh; ++i)
                for (std::size_t j = 0; j < nh; ++j)
                    beta[i * nh + j] = (1.0 - r) * t[i * nh + j] + (r / nh) * rowsum[i];
            std::fill(colsum.begin(), colsum.end(), 0.0);
            for (std::size_t i = 0; i < nh; ++i)
                for (std::size_t j = 0; j < nh; ++j) colsum[j] += beta[i * nh + j];
            std::vector<double> acc(ns, 0.0);
            for (std::size_t i = 0; i < nh; ++i)
                for (std::size_t j = 0; j < nh; ++j)
                    acc[i * nh + j] = (1.0 - r) * beta[i * nh + j] + (r / nh) * colsum[j];
            apply_grouped_corrections(*sp, nh, r, t, acc, false, st);
            bwd[bi] = acc;
        } else if (lk == nullptr) {
            std::fill(rowsum.begin(), rowsum.end(), 0.0);
            for (std::size_t i = 0; i < nh; ++i)
                for (std::size_t j = 0; j < nh; ++j) rowsum[i] += t[i * nh + j];
            for (std::size_t i = 0; i < nh; ++i)
                for (std::size_t j = 0; j < nh; ++j)
                    beta[i * nh + j] = (1.0 - r) * t[i * nh + j] + (r / nh) * rowsum[i];
            std::fill(colsum.begin(), colsum.end(), 0.0);
            for (std::size_t i = 0; i < nh; ++i)
                for (std::size_t j = 0; j < nh; ++j) colsum[j] += beta[i * nh + j];
            for (std::size_t i = 0; i < nh; ++i)
                for (std::size_t j = 0; j < nh; ++j)
                    bwd[bi][i * nh + j] = (1.0 - r) * beta[i * nh + j] + (r / nh) * colsum[j];
        } else {
            const double shift = linkage_shift(*lk);
            for (std::size_t i = 0; i < nh; ++i)
                for (std::size_t j = 0; j < nh; ++j) {
                    double acc = 0.0;
                    for (std::size_t i2 = 0; i2 < nh; ++i2) {
                        const double ti = (i == i2 ? (1.0 - r) : 0.0) + r / static_cast<double>(nh);
                        for (std::size_t j2 = 0; j2 < nh; ++j2) {
                            const double v = t[i2 * nh + j2];
                            if (v == 0.0) continue;
                            const double tj =
                                (j == j2 ? (1.0 - r) : 0.0) + r / static_cast<double>(nh);
                            acc += v * ti * tj * psi_weight(*lk, shift, i, j, i2, j2);
                        }
                    }
                    bwd[bi][i * nh + j] = acc;
                }
        }
        const double s = std::accumulate(bwd[bi].begin(), bwd[bi].end(), 0.0);
        if (s > 0.0) for (double& x : bwd[bi]) x /= s;
    }
    if (stats != nullptr) *stats = st;
}

}  // namespace panvar
