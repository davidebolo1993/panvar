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
    const std::vector<int>* truth_allele2,
    const CoverageEvidence* coverage) {

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
    std::vector<double> universe_size(nb, 0.0);
    std::vector<std::size_t> max_copies_of(nb, 1);

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

        // A marker the candidate pair does not carry is predicted at the error background, so observing
        // real counts there is nearly impossible and the marker becomes a veto. Under leave-one-out
        // that is guaranteed to happen: the sample's own allele is gone, so EVERY candidate lacks some
        // of the sample's sequence, and the one that lacks least wins whatever else it gets wrong.
        //
        // Measured at lpa's KIV-2 block, comparing the called pair against the most identical available
        // pair: 57 markers carried by one and not the other, holding 797 of 451,000 observed counts,
        // contributed +2704 log units, while the 951 markers both carry at differing copy number --
        // holding 404,138 counts, 90% of the data -- contributed +453. A five-thousandth of the reads
        // outvoted nine tenths of them, six to one. And the bias has a direction: a longer tandem array
        // carries more distinct unit variants, so it is more likely to contain any given marker, which
        // is why the call came out two repeat units too long.
        //
        // So every marker gets a mixture: mostly the allele's own prediction, and with probability
        // `marker_outlier` a broad component standing for sequence the candidate does not model. That
        // bounds what one marker can say without changing what many markers say together.
        // The outlier component is FLAT, not another Poisson. Unmodelled sequence is present at an
        // unknown copy number, so marginalising the count over that unknown gives a near-uniform
        // density across the range a marker could reach -- about 1/(lambda * max copies). A Poisson at
        // the block's mean count was tried first and is too narrow: it still charges ~20 log units for
        // a marker the candidate lacks, which halved the problem instead of removing it (the identity
        // oracle moved from emission rank 29 to 14 and the call did not change).
        const double eps = std::min(0.5, std::max(0.0, options.marker_outlier));
        double max_mult = 1.0;
        for (const auto& mset : panel.by_block[bi]) {
            for (const auto& [slot, m] : mset.nodes) { (void)slot; max_mult = std::max(max_mult, (double)m); }
        }
        max_copies_of[bi] = static_cast<std::size_t>(max_mult);
        const double broad_ll = -std::log(std::max(2.0, lambda * max_mult + 1.0));
        auto mix = [&](double o, double mean) {
            if (eps <= 0.0) return log_nb(o, mean, phi);
            const double a = std::log1p(-eps) + log_nb(o, mean, phi);
            const double b = std::log(eps) + broad_ll;
            const double hi = std::max(a, b);
            return hi + std::log(std::exp(a - hi) + std::exp(b - hi));
        };

        double baseline = 0.0;
        double obs_universe = 0.0;
        for (const std::uint32_t slot : universe) {
            baseline += mix(static_cast<double>(counts.node[slot]), mu);
            obs_universe += counts.node[slot];
        }
        baseline_of[bi] = baseline;
        obs_universe_of[bi] = obs_universe;
        universe_size[bi] = static_cast<double>(universe.size());

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

        // ---- total-mass term ----
        //
        // The per-marker product above is discounted by rho because adjacent syncmers share reads, so
        // thousands of markers do not carry thousands of independent observations. That discount is
        // right for COMPOSITION -- which markers are present, at what relative multiplicity -- but it
        // is wrong for the TOTAL, and the two are different statistics with different effective sample
        // sizes. One fragment contributes to many markers at once, so summing the block's counts adds
        // fragments rather than markers: the total's relative error is 1/sqrt(number of fragments
        // crossing the block), independent of how many markers it holds.
        //
        // Scored jointly under rho, the total is shrunk by the same factor as the composition and stops
        // deciding anything. Measured at lpa's KIV-2 array: the pair whose total length is correct
        // ranks 4270th of 457x457 by emission alone, while a pair one repeat unit short ranks first at
        // GQ 99 -- the likelihood was trading a uniform 3.7% shortfall in predicted counts, which costs
        // it almost nothing per marker, against unit-variant composition, which costs it a lot. The
        // total is exactly the statistic that distinguishes them, so it gets its own weight.
        //
        // Worked in multiplicity rather than in bases: informative-marker density is not a property of
        // sequence -- an allele with many confined markers and one with few can be the same length --
        // so a bases-to-multiplicity conversion would be a second assumption. The emission already
        // needs each pair's total multiplicity, so this costs nothing to compute.
        double mass_target = 0.0;
        double mass_sd = 0.0;
        if (options.mass_weight > 0.0 && lambda > 0.0 && !universe.empty()) {
            mass_target = std::max(0.0, obs_universe - static_cast<double>(universe.size()) * mu) / lambda;
            std::vector<std::size_t> lens = blocks[bi].allele_bp;
            std::sort(lens.begin(), lens.end());
            const double span = lens.empty() ? 0.0 : static_cast<double>(lens[lens.size() / 2]);
            // Fragments over both haplotypes of the block. lambda under-states the base depth (it is
            // per-syncmer, so it already carries the (readlen-k+1)/readlen factor), which makes this
            // an under-count and the resulting sd conservative.
            const double frags = options.fragment_len > 0.0 ? 2.0 * lambda * span / options.fragment_len : 0.0;
            mass_sd = frags > 1.0 && mass_target > 0.0 ? mass_target / std::sqrt(frags) : 0.0;
        }

        // Coverage emission for this block, if it was asked for. Same candidate set, same chain, same
        // filters -- only the evidence changes, so a difference in the result is attributable.
        const bool use_cov = coverage != nullptr && bi < coverage->use_block.size() &&
                             coverage->use_block[bi] != 0 &&
                             bi < coverage->block_allele_nodes.size() &&
                             !coverage->block_allele_nodes[bi].empty();
        if (use_cov) {
            const auto& av = coverage->block_allele_nodes[bi];
            const double lam = coverage->lambda > 0.0 ? coverage->lambda : 1.0;
            const double cmu = std::max(0.01, 0.02 * lam);
            // Only the nodes some allele here traverses. Anything else belongs to another block and
            // would add the same constant to every candidate while diluting the comparison.
            std::vector<std::uint32_t> cnodes;
            for (const auto& v : av) {
                for (std::size_t i = 0; i < v.size(); ++i) {
                    if (v[i] > 0) cnodes.push_back(static_cast<std::uint32_t>(i));
                }
            }
            std::sort(cnodes.begin(), cnodes.end());
            cnodes.erase(std::unique(cnodes.begin(), cnodes.end()), cnodes.end());
            const std::size_t kn2 = kept[bi].size();
            emis[bi].assign(kn2 * kn2, 0.0);
            explained[bi].assign(kn2 * kn2, 1.0);
            detected[bi].assign(kn2 * kn2, 1.0);
            for (std::size_t x = 0; x < kn2; ++x) {
                for (std::size_t y = 0; y < kn2; ++y) {
                    const std::size_t ax = kept[bi][x];
                    const std::size_t ay = kept[bi][y];
                    if (ax >= av.size() || ay >= av.size()) continue;
                    double ll = 0.0, obs_in = 0.0, pred_in = 0.0;
                    for (const std::uint32_t n : cnodes) {
                        const double m = static_cast<double>(av[ax][n] + av[ay][n]);
                        const double o = n < coverage->node.size() ? coverage->node[n] : 0.0;
                        ll += log_nb(o, lam * m + cmu, phi);
                        obs_in += o;
                        pred_in += lam * m + cmu;
                    }
                    emis[bi][x * kn2 + y] = ll;
                    detected[bi][x * kn2 + y] = pred_in > 0.0 ? obs_in / pred_in : 1.0;
                }
            }
            // The same effective-sample-size discount the marker emission gets, and for the same
            // reason: one fragment covers many adjacent nodes, so a sum over thousands of nodes is not
            // thousands of independent observations. Without it the coverage blocks are enormously
            // more decisive than the marker blocks they are chained to -- a scale mismatch, not a
            // difference in evidence -- and they drag the haplotype assignment across the whole locus.
            // Measured on lpa leave-one-out before this: routing four blocks to coverage took the pairs
            // from 15/13/18/11 exact blocks down to 12/11/15/7.
            double rho_cov = 1.0;
            if (options.fragment_len > 0.0 && !cnodes.empty()) {
                std::vector<std::size_t> lens = blocks[bi].allele_bp;
                std::sort(lens.begin(), lens.end());
                const double span = lens.empty() ? 0.0 : static_cast<double>(lens[lens.size() / 2]);
                rho_cov = std::min(1.0, (span / options.fragment_len) / static_cast<double>(cnodes.size()));
                if (!(rho_cov > 0.0)) rho_cov = 1e-6;
            }
            const double base_cov = *std::min_element(emis[bi].begin(), emis[bi].end());
            for (double& v : emis[bi]) v = base_cov + rho_cov * (v - base_cov);
            baseline_of[bi] = base_cov;
            continue;
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
                double obs_in = 0.0;
                double pred_in = 0.0;
                double mult_in = 0.0;
                for (const auto& [slot, m] : tot) {
                    obs_in += static_cast<double>(counts.node[slot]);
                    pred_in += lambda * m + mu;
                    mult_in += m;
                }
                // Shape and scale are scored separately, because the reads measure them at very
                // different precisions and folding both into one discounted product gets the second
                // one wrong.
                //
                // `scale` profiles the pair's overall depth out of the per-marker term: every predicted
                // mean is multiplied by the factor that makes the pair's total match the observed
                // total, so what remains is the SHAPE of the count vector -- which markers, in what
                // proportion -- and that is what the ESS discount is the right correction for.
                // `mass_ll` then puts the scale back as its own term at its own precision: the block's
                // total count gains one observation per fragment crossing it, not one per marker, so
                // its relative error is 1/sqrt(fragments) and does not shrink with rho.
                //
                // Scored the old way -- scale buried inside the discounted product -- lpa's KIV-2 block
                // preferred a pair two repeat units too long by 208 log units, because a uniform 4.4%
                // excess in predicted counts costs almost nothing per marker after rho, while the
                // unit-variant composition it buys is worth a great deal. The total, which says
                // unambiguously that the pair is 4.4% too big, was being shrunk 22-fold alongside it.
                const double scale = (options.mass_weight > 0.0 && pred_in > 0.0 && obs_in > 0.0)
                                         ? obs_in / pred_in : 1.0;
                double ll = 0.0;
                for (const auto& [slot, m] : tot) {
                    const double o = static_cast<double>(counts.node[slot]);
                    ll += weight_of(slot) * (mix(o, scale * (lambda * m + mu)) - mix(o, mu));
                }
                double mass_ll = 0.0;
                if (mass_sd > 0.0) {
                    const double z = (mult_in - mass_target) / mass_sd;
                    mass_ll = -0.5 * options.mass_weight * z * z;
                }
                emis[bi][x * kn + y] = baseline + rho * ll + mass_ll;
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

        // How much sequence the reads say is here, measured without reference to the panel's allele
        // list. Every marker in the block is expected at lambda*(copies in hap1 + copies in hap2), so
        // the observed mass over the whole marker universe, less the error background, divided by
        // lambda, IS the total multiplicity the sample carries. The called pair's own
        // multiplicity-per-base converts that to bases -- a ratio that is stable inside a tandem array,
        // where every allele is the same unit repeated.
        //
        // This matters because an allele call is quantised to the lengths the panel happens to hold.
        // When the sample's own array is not among them the call lands on the nearest whole array,
        // which at lpa's KIV-2 is a whole repeat unit away; the mass estimate is continuous and does
        // not have to round to a panel member.
        {
            const double bp1 = c.allele1 < blocks[bi].allele_bp.size()
                                   ? static_cast<double>(blocks[bi].allele_bp[c.allele1]) : 0.0;
            const double bp2 = c.allele2 < blocks[bi].allele_bp.size()
                                   ? static_cast<double>(blocks[bi].allele_bp[c.allele2]) : 0.0;
            c.called_bp = bp1 + bp2;
            double m_called = 0.0;
            if (c.allele1 < panel.by_block[bi].size()) {
                for (const auto& [slot, m] : panel.by_block[bi][c.allele1].nodes) { (void)slot; m_called += m; }
            }
            if (c.allele2 < panel.by_block[bi].size()) {
                for (const auto& [slot, m] : panel.by_block[bi][c.allele2].nodes) { (void)slot; m_called += m; }
            }
            c.max_copies = max_copies_of[bi];
            // Three copies of a marker is already past anything ordinary variation produces: a
            // substitution or a small indel changes a marker's presence, not how many times it recurs.
            // Recurrence at that level means the block contains a repeat, and it is the repeat that
            // makes the allele call and the mass estimate answer different questions.
            c.is_array = c.max_copies >= 3;
            const double lambda = bi < depth.size() ? depth[bi].lambda_hap : 0.0;
            const double density = c.called_bp > 0.0 ? m_called / c.called_bp : 0.0;
            if (lambda > 0.0 && density > 0.0) {
                const double m_hat = std::max(0.0, obs_universe_of[bi] - universe_size[bi] * mu) / lambda;
                c.mass_bp = m_hat / density;
                // The mass is a sum over fragments, not over markers: one fragment contributes to many
                // markers at once, so its relative error is set by how many fragments cross the block,
                // not by how many markers it holds. That is also why the per-marker ESS discount must
                // not be read as the precision of this number -- they measure different things.
                const double frags = options.fragment_len > 0.0
                                         ? 2.0 * lambda * c.called_bp / (2.0 * options.fragment_len) : 0.0;
                c.mass_bp_sd = frags > 1.0 ? c.mass_bp / std::sqrt(frags) : 0.0;
            }
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
         "\tcalled_bp\tmass_bp\tmass_bp_sd\tmax_copies\tblock_class"
         "\tevidence\tfilter\ttruth1\ttruth2\ttruth_rank\ttruth_delta\tinfluencers\n";
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
          << c.called_bp << '\t' << c.mass_bp << '\t' << c.mass_bp_sd << '\t'
          << c.max_copies << '\t' << (c.is_array ? "array" : "simple") << '\t'
          << c.evidence << '\t' << c.filter << '\t'
          << c.truth_allele1 << '\t' << c.truth_allele2 << '\t' << c.truth_emission_rank << '\t'
          << c.truth_emission_delta << '\t';
        for (std::size_t k = 0; k < c.influencers.size(); ++k) {
            if (k) g << ',';
            g << c.influencers[k];
        }
        if (c.influencers.empty()) g << '.';
        g << '\n';
    }
}

} // namespace panvar
