#include "panvar/genotype.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <unordered_map>

namespace panvar {
namespace {

// Negative binomial rather than Poisson: real short-read depth is overdispersed (GC, mappability,
// duplicates), and a Poisson emission would make the model overconfident -- the failure mode that
// matters most here, since a sizeable share of calls involve alleles the panel cannot supply.
double log_nb(double x, double mean, double phi) {
    if (mean <= 0.0) mean = 1e-9;
    if (phi <= 0.0 || !std::isfinite(phi)) {          // degenerate dispersion -> Poisson
        return -mean + x * std::log(mean) - std::lgamma(x + 1.0);
    }
    return std::lgamma(x + phi) - std::lgamma(phi) - std::lgamma(x + 1.0)
         + phi * std::log(phi / (phi + mean)) + x * std::log(mean / (phi + mean));
}

// Dispersion from the invariant depth anchors: they are carried by every allele at one copy, so
// their spread is coverage noise alone rather than genotype signal.
double estimate_phi(const ReadPanel& panel, const ReadCounts& counts) {
    std::vector<double> v;
    for (const auto& slots : panel.anchor_slots) {
        for (const std::uint32_t s : slots) v.push_back(static_cast<double>(counts.node[s]));
    }
    if (v.size() < 20) return 0.0;
    const double mean = std::accumulate(v.begin(), v.end(), 0.0) / static_cast<double>(v.size());
    double var = 0.0;
    for (const double x : v) var += (x - mean) * (x - mean);
    var /= static_cast<double>(v.size() - 1);
    if (var <= mean || mean <= 0.0) return 0.0;       // not overdispersed -> Poisson is fine
    return mean * mean / (var - mean);
}

} // namespace

std::vector<BlockCall> genotype_sample(
    const std::vector<Block>& chain,
    const std::vector<BlockAlleles>& blocks,
    const ReadPanel& panel,
    const ReadCounts& counts,
    const std::vector<BlockDepth>& depth,
    const std::vector<std::string>& haplotype_names,
    const GenotypeOptions& options,
    GenotypeSummary* summary,
    const std::vector<int>* truth_allele1,
    const std::vector<int>* truth_allele2) {

    const std::size_t nb = chain.size();
    const std::size_t nh = haplotype_names.size();
    std::vector<BlockCall> calls(nb);
    if (nh == 0 || nb == 0) return calls;

    const double phi = options.overdispersion > 0.0 ? options.overdispersion
                                                    : estimate_phi(panel, counts);
    double mu = options.error_background;
    if (mu <= 0.0) {
        double lam = 0.0;
        std::size_t n = 0;
        for (const BlockDepth& d : depth) { if (d.usable) { lam += d.lambda_hap; ++n; } }
        lam = n ? lam / static_cast<double>(n) : 1.0;
        mu = std::max(0.01, 0.02 * lam);
    }

    // haplotype -> allele index, per block (-1 when the haplotype does not traverse it).
    std::vector<std::vector<int>> allele_of(nb, std::vector<int>(nh, -1));
    for (std::size_t bi = 0; bi < nb; ++bi) {
        for (std::size_t hi = 0; hi < nh; ++hi) {
            const auto it = blocks[bi].allele_of.find(haplotype_names[hi]);
            if (it != blocks[bi].allele_of.end()) allele_of[bi][hi] = static_cast<int>(it->second);
        }
    }

    // ---- per-block emission over allele pairs ----
    std::vector<std::vector<double>> emis(nb);          // [block][x*kn + y], log scale
    std::vector<std::vector<std::uint32_t>> kept(nb);
    std::vector<std::vector<double>> explained(nb);
    std::vector<std::vector<double>> detected(nb);
    std::vector<double> baseline_of(nb, 0.0);
    std::vector<double> obs_universe_of(nb, 0.0);

    for (std::size_t bi = 0; bi < nb; ++bi) {
        const std::size_t na = panel.by_block[bi].size();
        if (na == 0) continue;
        const double lambda = depth[bi].usable ? depth[bi].lambda_hap : 1.0;

        // Prune candidates by how much of each allele's own marker set was actually detected: an
        // allele the sample does not carry leaves most of its discriminative markers unobserved.
        std::vector<std::pair<double, std::uint32_t>> score(na);
        for (std::size_t ai = 0; ai < na; ++ai) {
            const auto& ms = panel.by_block[bi][ai].nodes;
            std::size_t seen = 0;
            for (const auto& [slot, mult] : ms) { (void)mult; if (counts.node[slot] > 0) ++seen; }
            score[ai] = {ms.empty() ? 0.0 : static_cast<double>(seen) / static_cast<double>(ms.size()),
                         static_cast<std::uint32_t>(ai)};
        }
        std::sort(score.begin(), score.end(), [](const auto& x, const auto& y) {
            if (x.first != y.first) return x.first > y.first;
            return x.second < y.second;                 // total order: keep it reproducible
        });
        // Half the budget by that score, half spread evenly over the LENGTH distribution.
        //
        // Detected-marker fraction cannot rank alleles that differ mainly in copy number: they share
        // the repeat unit's syncmers, so nearly every allele scores alike and the top-K is decided by
        // the index tie-break. Measured on lpa's KIV-2 block, 457 alleles: the allele nearest the
        // truth in length was pruned before pairing. Stratifying by length guarantees the candidate
        // set spans the range whatever the score does, which needs no knowledge of the locus -- a
        // block whose alleles are all the same length simply gets its top-scoring alleles back, since
        // the stratified picks collapse onto them.
        const std::size_t keep_n = std::min(na, std::max<std::size_t>(2, options.max_alleles_per_block));
        kept[bi].reserve(keep_n);
        std::vector<char> taken(na, 0);
        const std::size_t by_score = keep_n - keep_n / 2;
        for (std::size_t i = 0; i < by_score && i < score.size(); ++i) {
            kept[bi].push_back(score[i].second);
            taken[score[i].second] = 1;
        }
        if (kept[bi].size() < keep_n) {
            std::vector<std::uint32_t> by_len(na);
            for (std::size_t i = 0; i < na; ++i) by_len[i] = static_cast<std::uint32_t>(i);
            std::sort(by_len.begin(), by_len.end(), [&](std::uint32_t x, std::uint32_t y) {
                const std::size_t bx = x < blocks[bi].allele_bp.size() ? blocks[bi].allele_bp[x] : 0;
                const std::size_t by = y < blocks[bi].allele_bp.size() ? blocks[bi].allele_bp[y] : 0;
                if (bx != by) return bx < by;
                return x < y;                           // total order: keep it reproducible
            });
            const std::size_t want = keep_n - kept[bi].size();
            for (std::size_t k = 0; k < want && kept[bi].size() < keep_n; ++k) {
                // Evenly spaced ranks across the sorted-by-length list, then the next untaken allele.
                std::size_t start = want > 1 ? (k * (na - 1)) / (want - 1) : na / 2;
                for (std::size_t step = 0; step < na; ++step) {
                    const std::size_t idx = (start + step) % na;
                    if (!taken[by_len[idx]]) {
                        kept[bi].push_back(by_len[idx]);
                        taken[by_len[idx]] = 1;
                        break;
                    }
                }
            }
        }
        std::sort(kept[bi].begin(), kept[bi].end());

        // Every candidate pair must be scored over the SAME marker set. Summing only over the markers
        // a pair carries would reward pairs with smaller support, since each extra term is another
        // negative log-probability -- which showed up as a spurious preference for homozygous calls.
        // Written as a per-block baseline over the full universe plus a correction on the union:
        //   E(a1,a2) = SUM_all logNB(obs; mu) + SUM_union [logNB(obs; lambda*(m1+m2)+mu) - logNB(obs; mu)]
        std::vector<std::uint32_t> universe;
        for (const auto& mset : panel.by_block[bi]) {
            for (const auto& [slot, mult] : mset.nodes) { (void)mult; universe.push_back(slot); }
        }
        std::sort(universe.begin(), universe.end());
        universe.erase(std::unique(universe.begin(), universe.end()), universe.end());
        // How many of this block's alleles carry each marker -- the basis of the specificity weight.
        std::unordered_map<std::uint32_t, std::uint32_t> carriers;
        for (const auto& mset : panel.by_block[bi]) {
            for (const auto& [slot, mult] : mset.nodes) { (void)mult; ++carriers[slot]; }
        }
        std::unordered_map<std::uint32_t, double> wt;
        if (options.carrier_weight > 0.0 && !universe.empty()) {
            const double na = static_cast<double>(std::max<std::size_t>(1, panel.by_block[bi].size()));
            double sum = 0.0;
            for (const std::uint32_t slot : universe) {
                const double c = static_cast<double>(std::max<std::uint32_t>(1, carriers[slot]));
                const double w = std::pow(na / c, options.carrier_weight);
                wt[slot] = w;
                sum += w;
            }
            const double scale = static_cast<double>(universe.size()) / std::max(1e-12, sum);
            for (auto& [slot, w] : wt) { (void)slot; w *= scale; }
        }
        auto weight_of = [&](std::uint32_t slot) {
            if (wt.empty()) return 1.0;
            const auto it = wt.find(slot);
            return it == wt.end() ? 1.0 : it->second;
        };

        double baseline = 0.0;
        double obs_universe = 0.0;
        for (const std::uint32_t slot : universe) {
            baseline += log_nb(static_cast<double>(counts.node[slot]), mu, phi);
            obs_universe += counts.node[slot];
        }
        baseline_of[bi] = baseline;
        obs_universe_of[bi] = obs_universe;

        // Effective-sample-size discount. Adjacent syncmers sit a few bp apart, so ~25 of them share
        // one 150 bp read; treating each as independent overstates the evidence and yields
        // overconfident GQ. The independent count is closer to span/fragment_len than to n_markers.
        double rho = 1.0;
        if (options.fragment_len > 0.0 && !universe.empty()) {
            std::vector<std::size_t> lens = blocks[bi].allele_bp;
            std::sort(lens.begin(), lens.end());
            const double span = lens.empty() ? 0.0 : static_cast<double>(lens[lens.size() / 2]);
            rho = std::min(1.0, (span / options.fragment_len) / static_cast<double>(universe.size()));
            if (!(rho > 0.0)) rho = 1e-6;
        }

        const std::size_t kn = kept[bi].size();
        emis[bi].assign(kn * kn, 0.0);
        explained[bi].assign(kn * kn, 0.0);
        detected[bi].assign(kn * kn, 0.0);
        for (std::size_t x = 0; x < kn; ++x) {
            for (std::size_t y = 0; y < kn; ++y) {
                // Expected count per marker is lambda * (copies in allele1 + copies in allele2), with
                // integer multiplicity and no cap -- which is what lets a tandem array with 20 copies
                // be modelled at all.
                std::unordered_map<std::uint32_t, std::uint32_t> tot;
                for (const auto& [slot, mult] : panel.by_block[bi][kept[bi][x]].nodes) tot[slot] += mult;
                for (const auto& [slot, mult] : panel.by_block[bi][kept[bi][y]].nodes) tot[slot] += mult;
                double ll = 0.0;
                double obs_in = 0.0;
                double pred_in = 0.0;
                for (const auto& [slot, m] : tot) {
                    const double o = static_cast<double>(counts.node[slot]);
                    ll += weight_of(slot) * (log_nb(o, lambda * m + mu, phi) - log_nb(o, mu, phi));
                    obs_in += o;
                    pred_in += lambda * m + mu;
                }
                emis[bi][x * kn + y] = baseline + rho * ll;
                // Share of the block's observed marker mass this pair accounts for. The denominator is
                // the UNION of the block's markers -- summing per allele would count every shared
                // marker once per carrier and make `explained` scale like 1/n_alleles.
                explained[bi][x * kn + y] = obs_universe > 0.0 ? obs_in / obs_universe : 0.0;
                // Raw ratio, NOT capped. Capping at 1 made this blind in one direction: a call that
                // UNDER-states multiplicity predicts fewer counts than the reads carry, so obs/pred
                // rises above 1 -- and the cap turned that into a clean 1.0. Measured on a synthetic
                // tandem array, 7 of 8 samples under-called the total by 9-54 kb while reporting
                // detected = 1.0. A fit ratio has to be two-sided to be a fit ratio.
                detected[bi][x * kn + y] = pred_in > 0.0 ? obs_in / pred_in : 1.0;
            }
        }
    }

    // Where does the truth's allele pair rank on emission alone? Separates "the emission is wrong"
    // from "linkage overrode a correct emission".
    if (truth_allele1 != nullptr && truth_allele2 != nullptr) {
        for (std::size_t bi = 0; bi < nb; ++bi) {
            const int t1 = (*truth_allele1)[bi];
            const int t2 = (*truth_allele2)[bi];
            calls[bi].truth_allele1 = t1;
            calls[bi].truth_allele2 = t2;
            if (t1 < 0 || t2 < 0 || kept[bi].empty()) continue;
            const auto xi = std::find(kept[bi].begin(), kept[bi].end(), static_cast<std::uint32_t>(t1));
            const auto yi = std::find(kept[bi].begin(), kept[bi].end(), static_cast<std::uint32_t>(t2));
            if (xi == kept[bi].end() || yi == kept[bi].end()) { calls[bi].truth_emission_rank = -2; continue; }
            const std::size_t kn = kept[bi].size();
            const double tv = emis[bi][static_cast<std::size_t>(xi - kept[bi].begin()) * kn +
                                       static_cast<std::size_t>(yi - kept[bi].begin())];
            double best = -std::numeric_limits<double>::infinity();
            int rank = 1;
            for (std::size_t x = 0; x < kn; ++x) {
                for (std::size_t y = x; y < kn; ++y) {
                    const double v = emis[bi][x * kn + y];
                    best = std::max(best, v);
                    if (v > tv + 1e-9) ++rank;
                }
            }
            calls[bi].truth_emission_rank = rank;
            calls[bi].truth_emission_delta = tv - best;
        }
    }

    // ---- diploid Li-Stephens forward-backward over haplotype pairs ----
    // Transitions factorize (each haplotype switches independently), so a step is O(n^2) via row and
    // column sums rather than O(n^4).
    auto emission_for = [&](std::size_t bi, std::size_t i, std::size_t j) -> double {
        const int a = allele_of[bi][i];
        const int b = allele_of[bi][j];
        if (a < 0 || b < 0 || kept[bi].empty()) return baseline_of[bi];
        const auto xi = std::lower_bound(kept[bi].begin(), kept[bi].end(), static_cast<std::uint32_t>(a));
        const auto yi = std::lower_bound(kept[bi].begin(), kept[bi].end(), static_cast<std::uint32_t>(b));
        // A pruned allele falls back to the all-background likelihood rather than -inf. Zeroing it
        // would remove that haplotype from the forward recursion for the rest of the chain, so one
        // weak block could silently veto the correct haplotype everywhere.
        if (xi == kept[bi].end() || *xi != static_cast<std::uint32_t>(a)) return baseline_of[bi];
        if (yi == kept[bi].end() || *yi != static_cast<std::uint32_t>(b)) return baseline_of[bi];
        const std::size_t kn = kept[bi].size();
        return emis[bi][static_cast<std::size_t>(xi - kept[bi].begin()) * kn +
                        static_cast<std::size_t>(yi - kept[bi].begin())];
    };

    const double r = std::min(0.5, options.recomb_rate / static_cast<double>(std::max<std::size_t>(1, nb)));
    std::vector<std::vector<double>> fwd(nb, std::vector<double>(nh * nh, 0.0));
    std::vector<std::vector<double>> bwd(nb, std::vector<double>(nh * nh, 0.0));

    auto normalize = [](std::vector<double>& v) {
        const double s = std::accumulate(v.begin(), v.end(), 0.0);
        if (s > 0.0) for (double& x : v) x /= s;
    };
    auto compute_emissions = [&](std::size_t bi, std::vector<double>& e) {
        e.assign(nh * nh, 0.0);
        double best = -std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i < nh; ++i)
            for (std::size_t j = 0; j < nh; ++j)
                best = std::max(best, emission_for(bi, i, j));
        for (std::size_t i = 0; i < nh; ++i)
            for (std::size_t j = 0; j < nh; ++j)
                e[i * nh + j] = std::exp(emission_for(bi, i, j) - best);
    };
    // Provenance re-runs the recursion once per block with that block's evidence removed. Neutralizing
    // means a flat emission -- the block still exists and still costs a recombination step, it simply
    // says nothing -- and the cache keeps the expensive part (the emissions) from being recomputed
    // n_blocks times.
    std::size_t neutral_block = std::numeric_limits<std::size_t>::max();
    std::vector<std::vector<double>> e_cache;
    auto block_emissions = [&](std::size_t bi, std::vector<double>& e) {
        if (bi == neutral_block) { e.assign(nh * nh, 1.0); return; }
        if (!e_cache.empty()) { e = e_cache[bi]; return; }
        compute_emissions(bi, e);
    };

    std::vector<double> e(nh * nh);
    std::vector<double> beta(nh * nh), rowsum(nh), colsum(nh);
    auto run_forward = [&]() {
    block_emissions(0, e);
    fwd[0] = e;
    normalize(fwd[0]);
    for (std::size_t bi = 1; bi < nb; ++bi) {
        const std::vector<double>& prev = fwd[bi - 1];
        std::fill(rowsum.begin(), rowsum.end(), 0.0);
        for (std::size_t i = 0; i < nh; ++i)
            for (std::size_t j = 0; j < nh; ++j) rowsum[i] += prev[i * nh + j];
        for (std::size_t i = 0; i < nh; ++i)
            for (std::size_t j = 0; j < nh; ++j)
                beta[i * nh + j] = (1.0 - r) * prev[i * nh + j] + (r / nh) * rowsum[i];
        std::fill(colsum.begin(), colsum.end(), 0.0);
        for (std::size_t i = 0; i < nh; ++i)
            for (std::size_t j = 0; j < nh; ++j) colsum[j] += beta[i * nh + j];
        block_emissions(bi, e);
        for (std::size_t i = 0; i < nh; ++i)
            for (std::size_t j = 0; j < nh; ++j)
                fwd[bi][i * nh + j] = e[i * nh + j] * ((1.0 - r) * beta[i * nh + j] + (r / nh) * colsum[j]);
        normalize(fwd[bi]);
    }
    };

    auto run_backward = [&]() {
    std::fill(bwd[nb - 1].begin(), bwd[nb - 1].end(), 1.0);
    normalize(bwd[nb - 1]);
    for (std::size_t bi = nb - 1; bi-- > 0;) {
        block_emissions(bi + 1, e);
        std::vector<double> t(nh * nh);
        for (std::size_t i = 0; i < nh * nh; ++i) t[i] = bwd[bi + 1][i] * e[i];
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
        normalize(bwd[bi]);
    }
    };

    run_forward();
    run_backward();

    // Most probable allele pair per block, from whatever fwd/bwd currently hold. Used both for the
    // reported call and, with a block neutralized, to see which blocks that call actually depended on.
    auto map_allele_pairs = [&]() {
        std::vector<std::pair<std::size_t, std::size_t>> out(nb, {0, 0});
        std::unordered_map<std::uint64_t, double> ap;
        for (std::size_t bi = 0; bi < nb; ++bi) {
            ap.clear();
            for (std::size_t i = 0; i < nh; ++i) {
                for (std::size_t j = 0; j < nh; ++j) {
                    const double p = fwd[bi][i * nh + j] * bwd[bi][i * nh + j];
                    if (p <= 0.0) continue;
                    const int a = allele_of[bi][i];
                    const int b = allele_of[bi][j];
                    if (a < 0 || b < 0) continue;
                    const std::uint32_t lo = static_cast<std::uint32_t>(std::min(a, b));
                    const std::uint32_t hi = static_cast<std::uint32_t>(std::max(a, b));
                    ap[(static_cast<std::uint64_t>(lo) << 32) | hi] += p;
                }
            }
            double bestp = -1.0;
            for (const auto& [key, p] : ap) {
                if (p > bestp) {
                    bestp = p;
                    out[bi] = {static_cast<std::size_t>(key >> 32),
                               static_cast<std::size_t>(key & 0xffffffffu)};
                }
            }
        }
        return out;
    };

    // ---- posteriors -> per-block allele pair ----
    double gq_sum = 0.0;
    std::size_t called = 0;
    std::size_t nocall = 0;
    std::size_t offpanel = 0;
    for (std::size_t bi = 0; bi < nb; ++bi) {
        BlockCall& c = calls[bi];
        c.block_index = chain[bi].index;
        c.is_bubble = chain[bi].kind == BlockKind::Bubble;
        c.bubble_id = chain[bi].bubble_id;
        {
            std::vector<std::uint32_t> u;
            for (const auto& mset : panel.by_block[bi]) {
                for (const auto& [slot, mult] : mset.nodes) { (void)mult; u.push_back(slot); }
            }
            std::sort(u.begin(), u.end());
            u.erase(std::unique(u.begin(), u.end()), u.end());
            c.n_markers = u.size();
        }

        std::vector<double> post(nh * nh);
        double tot = 0.0;
        for (std::size_t i = 0; i < nh * nh; ++i) { post[i] = fwd[bi][i] * bwd[bi][i]; tot += post[i]; }
        if (tot <= 0.0) { c.filter = "NOCALL"; ++nocall; continue; }
        for (double& x : post) x /= tot;

        // Collapse haplotype pairs onto the allele pair they induce -- the reportable unit. Many
        // haplotype pairs give the same alleles, and their posterior mass belongs together.
        std::unordered_map<std::uint64_t, double> allele_post;
        double best_hp = 0.0;
        for (std::size_t i = 0; i < nh; ++i) {
            for (std::size_t j = 0; j < nh; ++j) {
                const double p = post[i * nh + j];
                if (p <= 0.0) continue;
                if (p > best_hp) { best_hp = p; c.hap1 = i; c.hap2 = j; }
                const int a = allele_of[bi][i];
                const int b = allele_of[bi][j];
                if (a < 0 || b < 0) continue;
                const std::uint32_t lo = static_cast<std::uint32_t>(std::min(a, b));
                const std::uint32_t hi = static_cast<std::uint32_t>(std::max(a, b));
                allele_post[(static_cast<std::uint64_t>(lo) << 32) | hi] += p;
            }
        }
        c.hap_posterior = best_hp;
        if (allele_post.empty()) { c.filter = "NOCALL"; ++nocall; continue; }

        double p1 = 0.0;
        double p2 = 0.0;
        std::uint64_t best_key = 0;
        for (const auto& [key, p] : allele_post) {
            if (p > p1) { p2 = p1; p1 = p; best_key = key; }
            else if (p > p2) { p2 = p; }
        }
        c.allele1 = static_cast<std::size_t>(best_key >> 32);
        c.allele2 = static_cast<std::size_t>(best_key & 0xffffffffu);
        c.gq = (p1 <= 0.0) ? 0.0 : std::min(99.0, -10.0 * std::log10(std::max(1e-10, 1.0 - p1)));

        const auto xi = std::find(kept[bi].begin(), kept[bi].end(), static_cast<std::uint32_t>(c.allele1));
        const auto yi = std::find(kept[bi].begin(), kept[bi].end(), static_cast<std::uint32_t>(c.allele2));
        if (xi != kept[bi].end() && yi != kept[bi].end()) {
            const std::size_t kn = kept[bi].size();
            const std::size_t off = static_cast<std::size_t>(xi - kept[bi].begin()) * kn +
                                    static_cast<std::size_t>(yi - kept[bi].begin());
            c.explained = explained[bi][off];
            c.detected = detected[bi][off];
        }

        // A block with no markers of its own is NOT evidence-free: the forward-backward carries the
        // haplotype assignment in from its neighbours, which is the whole point of chaining blocks.
        // On ankrd36c only the two flank blocks retain markers after region uniqueness, yet all 23
        // blocks are called correctly from linkage alone -- so treating "no local markers" as a
        // no-call discarded 21 correct calls. Such blocks are marked LINKED and still face the same
        // GQ gate, so the user can tell locally-supported calls from inherited ones.
        // `explained` is undefined without local markers, so the off-panel test is skipped there
        // rather than being read as "the panel explains none of this".
        const bool has_markers = c.n_markers > 0;
        const bool measurable = has_markers && obs_universe_of[bi] >= 1.0;
        if (measurable && c.explained < options.min_explained) { c.filter = "OFFPANEL"; ++offpanel; }
        // The reads do not account for what the call predicts, so the block is under-covered and the
        // likelihood has been minimising predicted counts rather than matching them. GQ cannot see
        // this: it is a comparison between alternatives, and when every alternative fits badly the
        // least-bad one still wins by a wide margin.
        // Either direction is a misfit: below the bound the reads do not account for what the call
        // predicts (a coverage dropout), above its reciprocal the call does not account for what the
        // reads carry (multiplicity under-stated, which is how a copy-number call collapses).
        else if (has_markers && options.min_detected > 0.0 &&
                 (c.detected < options.min_detected || c.detected > 1.0 / options.min_detected)) {
            c.filter = "LOWCOV"; ++nocall;
        }
        else if (c.gq < options.min_gq) { c.filter = "LOWGQ"; ++nocall; }
        else { c.filter = "PASS"; ++called; gq_sum += c.gq; }
        // Provenance is orthogonal to quality: a block with no markers of its own is not a worse call,
        // it is a differently sourced one, and a reader must be able to filter on each separately.
        c.evidence = has_markers ? "local" : "linked";
    }

    // ---- provenance: which blocks did each call actually depend on? ----
    // Neutralize one block's evidence at a time and re-run the recursion; a block whose removal
    // changes another block's call is what determined it. Model-exact rather than a heuristic, and
    // it costs only the recursion (the emissions are cached), not the emissions themselves.
    if (options.provenance) {
        const std::vector<std::pair<std::size_t, std::size_t>> base = map_allele_pairs();
        // Posterior mass on a FIXED allele pair per block, so the perturbed run is compared against
        // the same hypothesis. Asking only whether the MAP flipped would report "nothing influenced
        // this" whenever two blocks redundantly support the same call -- which is the common case, and
        // exactly the case worth describing.
        auto pair_posterior = [&](const std::vector<std::pair<std::size_t, std::size_t>>& target) {
            std::vector<double> out(nb, 0.0);
            for (std::size_t bi = 0; bi < nb; ++bi) {
                double hit = 0.0;
                double tot = 0.0;
                for (std::size_t i = 0; i < nh; ++i) {
                    for (std::size_t j = 0; j < nh; ++j) {
                        const double p = fwd[bi][i * nh + j] * bwd[bi][i * nh + j];
                        if (p <= 0.0) continue;
                        const int a = allele_of[bi][i];
                        const int b = allele_of[bi][j];
                        if (a < 0 || b < 0) continue;
                        tot += p;
                        const std::size_t lo = static_cast<std::size_t>(std::min(a, b));
                        const std::size_t hi = static_cast<std::size_t>(std::max(a, b));
                        if (lo == target[bi].first && hi == target[bi].second) hit += p;
                    }
                }
                out[bi] = tot > 0.0 ? hit / tot : 0.0;
            }
            return out;
        };
        const std::vector<double> base_post = pair_posterior(base);
        e_cache.resize(nb);
        for (std::size_t bi = 0; bi < nb; ++bi) compute_emissions(bi, e_cache[bi]);
        std::vector<std::vector<std::pair<double, std::size_t>>> influence(nb);
        for (std::size_t j = 0; j < nb; ++j) {
            neutral_block = j;
            run_forward();
            run_backward();
            const std::vector<double> alt = pair_posterior(base);
            for (std::size_t b = 0; b < nb; ++b) {
                const double drop = base_post[b] - alt[b];
                if (drop > 0.01) influence[b].emplace_back(drop, j);
            }
        }
        neutral_block = std::numeric_limits<std::size_t>::max();
        e_cache.clear();
        e_cache.shrink_to_fit();
        for (std::size_t b = 0; b < nb; ++b) {
            std::sort(influence[b].begin(), influence[b].end(),
                      [](const auto& x, const auto& y) {
                          return x.first != y.first ? x.first > y.first : x.second < y.second;
                      });
            for (const auto& [drop, j] : influence[b]) {
                (void)drop;
                if (calls[b].influencers.size() >= 5) break;
                calls[b].influencers.push_back(j);
            }
            // Named for the block that contributes most, so the label answers "where did this call
            // come from" rather than "what could have overturned it".
            if (influence[b].empty()) calls[b].provenance = "none";
            else {
                const std::size_t top = influence[b].front().second;
                calls[b].provenance = top == b ? "self"
                                    : (top + 1 == b || b + 1 == top) ? "neighbours" : "distant";
            }
            calls[b].evidence = calls[b].provenance;   // the finer answer, when we paid for it
        }
    }

    if (summary != nullptr) {
        summary->blocks = nb;
        summary->called = called;
        summary->no_calls = nocall;
        summary->off_panel = offpanel;
        summary->mean_gq = called ? gq_sum / static_cast<double>(called) : 0.0;
        summary->overdispersion = phi;
        summary->error_background = mu;
        double lam = 0.0;
        std::size_t n = 0;
        for (const BlockDepth& d : depth) { if (d.usable) { lam += d.lambda_hap; ++n; } }
        summary->lambda_hap = n ? lam / static_cast<double>(n) : 0.0;
    }
    return calls;
}

void write_genotypes(
    const std::string& out_prefix,
    const std::vector<Block>& chain,
    const std::vector<BlockAlleles>& blocks,
    const std::vector<BlockCall>& calls,
    const std::vector<std::string>& haplotype_names) {

    // One row per block, in chain order, carrying everything needed to read a call without joining
    // another file: where the block is (source -> sink), what was called, which panel haplotypes it
    // came from, how confident, WHY, and whether it passed quality.
    //
    // `evidence` and `filter` answer different questions and used to be conflated in one column.
    // `evidence` is provenance -- did this block's own syncmers decide it (`local`), or was it carried
    // by the rest of the chain (`linked`)? With --provenance that refines to which blocks. `filter` is
    // quality alone. A LINKED call is not lower quality, it is differently sourced, and a reader needs
    // to be able to filter on each independently.
    const std::string gpath = out_prefix + ".genotypes.tsv";
    std::ofstream g(gpath);
    if (!g) throw std::runtime_error("genotype: cannot write " + gpath);
    g << "block_index\tblock_kind\tbubble_id\tsource\tsink\tn_alleles\tn_markers"
         "\tallele1\tallele2\thaplotype1\thaplotype2\thap_posterior\tgq\texplained\tdetected"
         "\tevidence\tfilter\ttruth1\ttruth2\ttruth_rank\tinfluencers\n";
    for (std::size_t bi = 0; bi < calls.size(); ++bi) {
        const BlockCall& c = calls[bi];
        g << c.block_index << '\t'
          << (chain[bi].kind == BlockKind::Bubble ? "bubble"
              : chain[bi].kind == BlockKind::Flank ? "flank" : "backbone")
          << '\t' << c.bubble_id << '\t'
          << (chain[bi].source.empty() ? "." : chain[bi].source) << '\t'
          << (chain[bi].sink.empty() ? "." : chain[bi].sink) << '\t'
          << blocks[bi].n_alleles << '\t' << c.n_markers << '\t'
          << c.allele1 << '\t' << c.allele2 << '\t'
          << (c.hap1 < haplotype_names.size() ? haplotype_names[c.hap1] : ".") << '\t'
          << (c.hap2 < haplotype_names.size() ? haplotype_names[c.hap2] : ".") << '\t'
          << c.hap_posterior << '\t' << c.gq << '\t' << c.explained << '\t' << c.detected << '\t'
          << c.evidence << '\t' << c.filter << '\t'
          << c.truth_allele1 << '\t' << c.truth_allele2 << '\t' << c.truth_emission_rank << '\t';
        for (std::size_t k = 0; k < c.influencers.size(); ++k) {
            if (k) g << ',';
            g << c.influencers[k];
        }
        if (c.influencers.empty()) g << '.';
        g << '\n';
    }
}

} // namespace panvar
