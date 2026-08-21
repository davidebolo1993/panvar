#include "panvar/pangenie_model.hpp"

#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <unordered_set>

namespace panvar {
namespace {

// src/probabilitytable.cpp:7 -- the CN0 geometric parameter, a step function of coverage.
double error_param(double kmer_coverage) {
    if (kmer_coverage < 10.0) return 0.99;
    if (kmer_coverage < 20.0) return 0.95;
    if (kmer_coverage < 40.0) return 0.90;
    return 0.80;
}

// src/probabilitytable.cpp:75 -- computed in log space there too, for the same overflow reason.
long double poisson_pmf(long double mean, unsigned int value) {
    long double log_fact = 0.0L;
    for (unsigned int i = 1; i <= value; ++i) log_fact += std::log(static_cast<long double>(i));
    if (mean <= 0.0L) return value == 0 ? 1.0L : 0.0L;
    return std::exp(-mean + static_cast<long double>(value) * std::log(mean) - log_fact);
}

// src/probabilitytable.cpp:83
long double geometric_pmf(long double p, unsigned int value) {
    return std::pow(1.0L - p, static_cast<long double>(value)) * p;
}

// src/copynumber.cpp:22 -- the regularised constructor normalises, so the three copy-number
// probabilities sum to one and no single k-mer can drive the product to zero.
struct CopyNumber {
    long double p[3] = {0.0L, 0.0L, 0.0L};
};

CopyNumber copy_number_probs(double coverage, unsigned int count, long double reg) {
    const long double p0 = geometric_pmf(error_param(coverage), count);
    const long double p1 = poisson_pmf(static_cast<long double>(coverage) / 2.0L, count);
    const long double p2 = poisson_pmf(static_cast<long double>(coverage), count);
    CopyNumber cn;
    if (reg > 0.0L) {
        const long double sum = p0 + p1 + p2 + 3.0L * reg;
        cn.p[0] = (p0 + reg) / sum;
        cn.p[1] = (p1 + reg) / sum;
        cn.p[2] = 1.0L - cn.p[0] - cn.p[1];
    } else {
        cn.p[0] = p0; cn.p[1] = p1; cn.p[2] = p2;
    }
    return cn;
}

// src/transitionprobabilitycomputer.cpp:14 -- three switch classes, distance-scaled.
struct Transition {
    long double p[3] = {1.0L, 1.0L, 1.0L};
};

Transition transition_for(std::size_t from_pos, std::size_t to_pos, double recombrate,
                          std::size_t nr_paths, double effective_N) {
    const long double distance = static_cast<long double>(to_pos - from_pos) * 0.000004L *
                                 static_cast<long double>(recombrate) *
                                 static_cast<long double>(effective_N);
    const long double n = static_cast<long double>(std::max<std::size_t>(1, nr_paths));
    const long double recomb = (1.0L - std::exp(-distance / n)) * (1.0L / n);
    const long double no_recomb = std::exp(-distance / n) + recomb;
    Transition t;
    t.p[0] = no_recomb * no_recomb;
    t.p[1] = no_recomb * recomb;
    t.p[2] = recomb * recomb;
    return t;
}

} // namespace

std::vector<PanGenieCall> genotype_pangenie(
    const std::vector<Block>& chain,
    const std::vector<BlockAlleles>& blocks,
    const ReadPanel& panel,
    const ReadCounts& counts,
    const std::vector<BlockDepth>& depth,
    const std::vector<std::string>& haplotype_names,
    const PanGenieOptions& options) {

    const std::size_t nb = chain.size();
    const std::size_t nh = haplotype_names.size();
    std::vector<PanGenieCall> calls(nb);
    if (nb == 0 || nh == 0) return calls;

    // haplotype -> allele per block, as their ColumnIndexer provides.
    std::vector<std::vector<int>> allele_of(nb, std::vector<int>(nh, -1));
    for (std::size_t bi = 0; bi < nb; ++bi) {
        for (std::size_t hi = 0; hi < nh; ++hi) {
            const auto it = blocks[bi].allele_of.find(haplotype_names[hi]);
            if (it != blocks[bi].allele_of.end()) allele_of[bi][hi] = static_cast<int>(it->second);
        }
    }

    // Their marker rule. A k-mer must occur EXACTLY ONCE in EXACTLY ONE allele; presence is stored,
    // not count. Everything a tandem array contributes is discarded here by construction, because a
    // unit k-mer occurs N times in every allele that carries the array.
    std::vector<std::vector<std::uint32_t>> block_markers(nb);
    std::vector<std::vector<std::vector<char>>> on_allele(nb);   // [block][marker][allele] presence
    for (std::size_t bi = 0; bi < nb; ++bi) {
        const auto& per_allele = panel.by_block[bi];
        if (per_allele.empty()) continue;
        std::unordered_map<std::uint32_t, std::pair<std::uint32_t, std::uint32_t>> seen;  // slot -> (alleles, mult)
        for (std::size_t ai = 0; ai < per_allele.size(); ++ai) {
            for (const auto& [slot, mult] : per_allele[ai].nodes) {
                auto& e = seen[slot];
                e.first += 1;
                e.second = std::max(e.second, mult);
            }
        }
        for (const auto& [slot, e] : seen) {
            if (e.first == 1 && e.second == 1) block_markers[bi].push_back(slot);
        }
        std::sort(block_markers[bi].begin(), block_markers[bi].end());
        on_allele[bi].assign(block_markers[bi].size(), std::vector<char>(per_allele.size(), 0));
        for (std::size_t mi = 0; mi < block_markers[bi].size(); ++mi) {
            for (std::size_t ai = 0; ai < per_allele.size(); ++ai) {
                for (const auto& [slot, mult] : per_allele[ai].nodes) {
                    (void)mult;
                    if (slot == block_markers[bi][mi]) { on_allele[bi][mi][ai] = 1; break; }
                }
            }
        }
        calls[bi].n_markers = block_markers[bi].size();
    }

    // src/emissionprobabilitycomputer.cpp:37 -- a product over k-mers of the copy-number probability
    // at the expected count, which is presence(a1) + presence(a2) and so lies in {0,1,2}. With no
    // markers every pair scores 1.0, which their code calls "all_zeros -> uniform".
    std::vector<std::vector<long double>> emis(nb);
    for (std::size_t bi = 0; bi < nb; ++bi) {
        const std::size_t na = blocks[bi].n_alleles;
        emis[bi].assign(na * na, 1.0L);
        if (block_markers[bi].empty()) { calls[bi].undefined = true; continue; }
        const double cov = depth[bi].usable ? depth[bi].lambda_hap * 2.0 : 1.0;   // their "kmer_coverage"
        std::vector<CopyNumber> cn_of(block_markers[bi].size());
        for (std::size_t mi = 0; mi < block_markers[bi].size(); ++mi) {
            cn_of[mi] = copy_number_probs(cov, counts.node[block_markers[bi][mi]],
                                          static_cast<long double>(options.regularization));
        }
        bool all_zero = true;
        for (std::size_t a1 = 0; a1 < na; ++a1) {
            for (std::size_t a2 = 0; a2 < na; ++a2) {
                long double r = 1.0L;
                for (std::size_t mi = 0; mi < block_markers[bi].size(); ++mi) {
                    const int expected = on_allele[bi][mi][a1] + on_allele[bi][mi][a2];
                    r *= cn_of[mi].p[expected];
                }
                emis[bi][a1 * na + a2] = r;
                if (r > 0.0L) all_zero = false;
            }
        }
        if (all_zero) {
            std::fill(emis[bi].begin(), emis[bi].end(), 1.0L);
            calls[bi].undefined = true;
        }
    }

    // src/hmm.cpp:196 -- forward, then backward, over ordered path pairs, with the same factorised
    // helper sums and per-column normalisation.
    std::vector<std::vector<long double>> fwd(nb, std::vector<long double>(nh * nh, 0.0L));
    std::vector<std::vector<long double>> bwd(nb, std::vector<long double>(nh * nh, 0.0L));
    auto emission_at = [&](std::size_t bi, std::size_t h1, std::size_t h2) -> long double {
        const int a1 = allele_of[bi][h1];
        const int a2 = allele_of[bi][h2];
        const std::size_t na = blocks[bi].n_alleles;
        if (a1 < 0 || a2 < 0 || na == 0) return 1.0L;
        return emis[bi][static_cast<std::size_t>(a1) * na + static_cast<std::size_t>(a2)];
    };
    // Base-pair position, not block index. PanGenie scales its transition by the distance between
    // consecutive variants (transitionprobabilitycomputer.cpp:14), so standing the block index in for
    // it made every interval one unit apart and removed the distance term the port exists to reproduce.
    std::vector<double> block_bp(nb, 0.0);
    {
        double cum = 0.0;
        for (std::size_t b = 0; b < nb; ++b) {
            block_bp[b] = cum;
            std::vector<std::size_t> lens = blocks[b].allele_bp;
            if (!lens.empty()) {
                std::sort(lens.begin(), lens.end());
                cum += static_cast<double>(lens[lens.size() / 2]);
            }
        }
    }
    auto pos_of = [&](std::size_t bi) { return block_bp[bi]; };

    for (std::size_t bi = 0; bi < nb; ++bi) {
        std::vector<long double> hi(nh, 0.0L), hj(nh, 0.0L);
        long double hij = 0.0L;
        Transition tr;
        if (bi > 0) {
            for (std::size_t a = 0; a < nh; ++a) {
                for (std::size_t b = 0; b < nh; ++b) {
                    const long double v = fwd[bi - 1][a * nh + b];
                    hi[a] += v; hj[b] += v; hij += v;
                }
            }
            tr = transition_for(pos_of(bi - 1), pos_of(bi), options.recombrate, nh, options.effective_N);
        }
        long double norm = 0.0L;
        for (std::size_t a = 0; a < nh; ++a) {
            for (std::size_t b = 0; b < nh; ++b) {
                long double prev = 1.0L;
                if (bi > 0) {
                    const long double self = fwd[bi - 1][a * nh + b];
                    prev = tr.p[0] * self
                         + tr.p[1] * (hi[a] + hj[b] - 2.0L * self)
                         + tr.p[2] * (hij - hi[a] - hj[b] + self);
                }
                const long double cur = prev * emission_at(bi, a, b);
                fwd[bi][a * nh + b] = cur;
                norm += cur;
            }
        }
        if (norm > 0.0L) for (long double& x : fwd[bi]) x /= norm;
    }

    for (std::size_t bi = nb; bi-- > 0;) {
        if (bi + 1 == nb) { std::fill(bwd[bi].begin(), bwd[bi].end(), 1.0L); continue; }
        std::vector<long double> t(nh * nh, 0.0L);
        for (std::size_t a = 0; a < nh; ++a) {
            for (std::size_t b = 0; b < nh; ++b) {
                t[a * nh + b] = bwd[bi + 1][a * nh + b] * emission_at(bi + 1, a, b);
            }
        }
        std::vector<long double> hi(nh, 0.0L), hj(nh, 0.0L);
        long double hij = 0.0L;
        for (std::size_t a = 0; a < nh; ++a) {
            for (std::size_t b = 0; b < nh; ++b) {
                const long double v = t[a * nh + b];
                hi[a] += v; hj[b] += v; hij += v;
            }
        }
        const Transition tr =
            transition_for(pos_of(bi), pos_of(bi + 1), options.recombrate, nh, options.effective_N);
        long double norm = 0.0L;
        for (std::size_t a = 0; a < nh; ++a) {
            for (std::size_t b = 0; b < nh; ++b) {
                const long double self = t[a * nh + b];
                const long double v = tr.p[0] * self
                                    + tr.p[1] * (hi[a] + hj[b] - 2.0L * self)
                                    + tr.p[2] * (hij - hi[a] - hj[b] + self);
                bwd[bi][a * nh + b] = v;
                norm += v;
            }
        }
        if (norm > 0.0L) for (long double& x : bwd[bi]) x /= norm;
    }

    // src/hmm.cpp:368 -- posterior mass is accumulated onto the ALLELE pair the path pair induces.
    for (std::size_t bi = 0; bi < nb; ++bi) {
        calls[bi].block_index = chain[bi].index;
        std::unordered_map<std::uint64_t, long double> ap;
        long double total = 0.0L;
        for (std::size_t a = 0; a < nh; ++a) {
            for (std::size_t b = 0; b < nh; ++b) {
                const long double p = fwd[bi][a * nh + b] * bwd[bi][a * nh + b];
                if (p <= 0.0L) continue;
                const int x = allele_of[bi][a];
                const int y = allele_of[bi][b];
                if (x < 0 || y < 0) continue;
                const std::uint32_t lo = static_cast<std::uint32_t>(std::min(x, y));
                const std::uint32_t hi2 = static_cast<std::uint32_t>(std::max(x, y));
                ap[(static_cast<std::uint64_t>(lo) << 32) | hi2] += p;
                total += p;
            }
        }
        long double best = 0.0L;
        std::uint64_t key = 0;
        for (const auto& [k, p] : ap) if (p > best) { best = p; key = k; }
        calls[bi].allele1 = static_cast<std::size_t>(key >> 32);
        calls[bi].allele2 = static_cast<std::size_t>(key & 0xffffffffu);
        const long double frac = total > 0.0L ? best / total : 0.0L;
        calls[bi].gq = frac <= 0.0L ? 0.0
                                    : std::min(99.0, -10.0 * std::log10(std::max(1e-10,
                                                 1.0 - static_cast<double>(frac))));
    }
    return calls;
}

} // namespace panvar
