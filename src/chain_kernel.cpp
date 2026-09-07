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

void chain_forward_backward(
    std::size_t nh, std::size_t nb, double r,
    const std::function<void(std::size_t, std::vector<double>&)>& block_emissions,
    const std::vector<ChainEdgeLinkage>* edges,
    std::vector<std::vector<double>>& fwd,
    std::vector<std::vector<double>>& bwd,
    ChainKernelStats* stats) {
    if (nh == 0 || nb == 0) return;
    const std::size_t ns = nh * nh;
    fwd.assign(nb, std::vector<double>(ns, 0.0));
    bwd.assign(nb, std::vector<double>(ns, 0.0));
    ChainKernelStats st;

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
        const ChainEdgeLinkage* lk = edge_at(bi);
        block_emissions(bi, e);
        if (lk == nullptr) {
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
        const ChainEdgeLinkage* lk = edge_at(bi + 1);
        std::vector<double> t(ns);
        for (std::size_t k = 0; k < ns; ++k) t[k] = bwd[bi + 1][k] * e[k];
        if (lk == nullptr) {
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
