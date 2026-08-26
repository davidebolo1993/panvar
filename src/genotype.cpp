#include "panvar/genotype.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <unordered_map>
#include <cstdio>
#include <cstdlib>

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
        // the index tie-break, so the allele nearest the truth in length can be pruned before
        // pairing. Stratifying by length guarantees the candidate
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

        // A marker the candidate pair does not carry is predicted at the error background, so real
        // counts there are nearly impossible and the marker becomes a veto. Under leave-one-out that is
        // guaranteed: the sample's own allele is gone, so EVERY candidate lacks some of its sequence,
        // and the one that lacks least wins whatever else it gets wrong. The bias has a direction --
        // a longer tandem array carries more distinct unit variants, so it is likelier to contain any
        // given marker -- which is why such calls come out too long.
        //
        // So every marker gets a mixture: mostly the allele's own prediction, and with probability
        // `marker_outlier` a broad component standing for sequence the candidate does not model. That
        // bounds what one marker can say without changing what many markers say together.
        //
        // The outlier component is FLAT, not another Poisson. Unmodelled sequence is present at an
        // unknown copy number, so marginalising the count over that unknown gives a near-uniform
        // density across the range a marker could reach, about 1/(lambda * max copies). A Poisson at the
        // block's mean is too narrow: it still charges tens of log units for an absent marker, which
        // halves the problem rather than removing it.
        const double eps = std::min(0.5, std::max(0.0, options.marker_outlier));
        double max_mult = 1.0;
        for (const auto& mset : panel.by_block[bi]) {
            for (const auto& [slot, m] : mset.nodes) { (void)slot; max_mult = std::max(max_mult, (double)m); }
        }
        max_copies_of[bi] = static_cast<std::size_t>(max_mult);
        const double broad_ll = -std::log(std::max(2.0, lambda * max_mult + 1.0));
        // log_nb(o, mean) = G(o) + H(o, mean), where
        //   G(o) = lgamma(o+phi) - lgamma(phi) - lgamma(o+1)   depends ONLY on the observation
        //   H(o, mean) = phi*log(phi/(phi+mean)) + o*log(mean/(phi+mean))
        // Every emission term is a DIFFERENCE of two log_nb at the same o, so G cancels exactly.
        // Computing it was six lgamma calls per marker-pair summing to zero, and lgamma dominated
        // the inner loop: at 457 candidates one block is ~100M evaluations.
        auto nb_h = [&](double o, double mean) {
            if (mean <= 0.0) mean = 1e-9;
            if (phi <= 0.0 || !std::isfinite(phi)) return -mean + o * std::log(mean);
            return phi * std::log(phi / (phi + mean)) + o * std::log(mean / (phi + mean));
        };
        auto mix = [&](double o, double mean) {
            if (eps <= 0.0) return log_nb(o, mean, phi);
            const double a = std::log1p(-eps) + log_nb(o, mean, phi);
            const double b = std::log(eps) + broad_ll;
            const double hi = std::max(a, b);
            return hi + std::log(std::exp(a - hi) + std::exp(b - hi));
        };

        // The difference of two mixed likelihoods at one observation. Without the outlier component
        // the G terms cancel and only H is needed; with it the mixture is not a plain difference and
        // the full form is used, which is the price of that option rather than of every run.
        auto mix_diff = [&](double o, double mean_a, double mean_b) {
            if (eps <= 0.0) return nb_h(o, mean_a) - nb_h(o, mean_b);
            return mix(o, mean_a) - mix(o, mean_b);
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
            // Independent observations = the fragment-length WINDOWS the markers actually occupy, not
            // the block's span. Using the span assumes the markers are spread across it; a deletion
            // junction puts all of them inside a single fragment, where they share reads and move
            // together, and treating five such markers as four independent votes is what turned a
            // coverage fluctuation into a confident wrong homozygous call.
            const double clumps = bi < panel.marker_clumps.size() && panel.marker_clumps[bi] > 0.0
                ? panel.marker_clumps[bi] : span / options.fragment_len;
            rho = std::min(1.0, clumps / static_cast<double>(universe.size()));
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
        // deciding anything. At a large array that is decisive: the pair with the correct total length
        // can rank thousands of places below one a repeat unit short, because the likelihood trades a
        // uniform shortfall in predicted counts -- almost free per marker -- against unit-variant
        // composition, which is not. The
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
            // Left unscaled, routing even a few blocks to coverage measurably degrades the rest of the
            // chain.
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

        // Independent observations behind the block's TOTAL count: fragments crossing it, not markers.
        double scale_neff = 0.0;
        if (options.fragment_len > 0.0) {
            std::vector<std::size_t> lens2 = blocks[bi].allele_bp;
            std::sort(lens2.begin(), lens2.end());
            const double span2 = lens2.empty() ? 0.0 : static_cast<double>(lens2[lens2.size() / 2]);
            scale_neff = 2.0 * lambda * span2 / options.fragment_len;
        }
        // Use the block's TOTAL only where SHAPE cannot decide, and decide that from the PANEL rather
        // than from the sample, so the rule cannot be tuned by the data it judges.
        //
        // Shape can only separate alleles whose marker PROPORTIONS differ. When an allele contributes
        // nothing -- a bypass, a haplotype that deletes the block -- every genotype containing it has
        // the same proportions as the homozygote and only magnitude separates them, so there the total
        // is the entire signal. Everywhere else shape is available and the total is the more fragile
        // half: it is what a mis-estimated lambda corrupts, what a paralogue inflates, and what
        // divergence from a folded consensus deflates.
        //
        // A data-driven alternative -- also use the total wherever alleles differ in magnitude -- was
        // measured and bought nothing while costing the folded-consensus case. The simpler rule stands.
        //
        // The cost to state: at a tandem array with no bypass allele, copy number IS the scale and this
        // rule discards it, so `called_bp` there is less accurate. `mass_bp` remains the copy-number
        // answer at such a block, which is what block_class=array already says.
        if (blocks[bi].bypass_allele < 0) scale_neff = 0.0;
        double cov_target = 0.0, cov_sd = 0.0;
        if (coverage != nullptr && bi < coverage->target_bp.size() && coverage->target_bp[bi] > 0.0 &&
            bi < coverage->use_block.size() && coverage->use_block[bi] != 0) {
            cov_target = coverage->target_bp[bi];
            cov_sd = coverage->target_sd[bi];
        }

        // Target total multiplicity implied by the observed marker mass, for the array window.
        double win_target = 0.0;
        // Not where an allele contributes nothing: there the observed mass legitimately falls short of
        // any full traversal, so a window built from it excludes the right answer. Same gate the scale
        // term uses, and for the same reason.
        if (options.mass_window > 0.0 && max_copies_of[bi] >= 3 && lambda > 0.0 &&
            blocks[bi].bypass_allele < 0) {
            win_target = std::max(0.0, obs_universe - static_cast<double>(universe.size()) * mu) / lambda;

        }

        // The floor a PRUNED allele falls back to has to sit on the same scale as the alleles that
        // survived pruning, or the chain compares a multinomial against a negative binomial and the
        // pruned pair wins or loses for arithmetic reasons. In the default branch `ll` is a CORRECTION
        // over `baseline`, so a candidate carrying nothing scores exactly `baseline` and the existing
        // fallback is already right. The compositional and robust branches are absolute, so their floor
        // has to be computed the same way they compute everything else: by scoring the null candidate,
        // one that carries no marker at all. Pruning fires at every block above --max-alleles (64),
        // which is exactly the large array blocks those two modes were built for.
        if (options.compositional) {
            const double pred_tot = static_cast<double>(universe.size()) * mu;
            double cll = 0.0;
            if (pred_tot > 0.0) {
                for (const std::uint32_t slot : universe) {
                    const double o = static_cast<double>(counts.node[slot]);
                    if (o <= 0.0) continue;
                    cll += o * std::log(std::max(1e-300, mu / pred_tot));
                }
            }
            double sll = 0.0;
            if (scale_neff > 0.0 && pred_tot > 0.0 && obs_universe > 0.0) {
                const double c = std::max(1.0, obs_universe / scale_neff);
                const double oo = obs_universe / c;
                const double pp = pred_tot / c;
                sll = -pp + oo * std::log(std::max(1e-300, pp)) - std::lgamma(oo + 1.0);
            }
            baseline_of[bi] = rho * cll + options.scale_weight * sll;
        } else if (options.robust_c > 0.0) {
            const double c = options.robust_c;
            double ll = 0.0;
            for (const std::uint32_t slot : universe) {
                const double o = static_cast<double>(counts.node[slot]);
                const double z = (o - mu) / std::sqrt(mu + 1.0);
                const double a = std::fabs(z);
                ll -= weight_of(slot) * (a <= c ? 0.5 * z * z : c * (a - 0.5 * c));
            }
            baseline_of[bi] = rho * ll;
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
                // With the scale buried inside the discounted product, a large array prefers a pair
                // several repeat units too long: a uniform excess in predicted counts costs almost
                // nothing per marker after rho, while the unit-variant composition it buys is worth a
                // great deal. The total, which says
                // unambiguously that the pair is 4.4% too big, was being shrunk 22-fold alongside it.
                const double scale = (options.mass_weight > 0.0 && pred_in > 0.0 && obs_in > 0.0)
                                         ? obs_in / pred_in : 1.0;
                double ll = 0.0;
                if (options.robust_c > 0.0) {
                    // Huber on the Pearson residual, negated so larger is better and the scale stays
                    // comparable with the likelihood path. Every marker of the block's universe is
                    // scored, so a candidate is charged for markers it lacks -- but only up to the cap.
                    const double c = options.robust_c;
                    for (const std::uint32_t slot : universe) {
                        const auto it = tot.find(slot);
                        const double m = it == tot.end() ? 0.0 : static_cast<double>(it->second);
                        const double pred = lambda * m + mu;
                        const double o = static_cast<double>(counts.node[slot]);
                        const double z = (o - pred) / std::sqrt(pred + 1.0);
                        const double a = std::fabs(z);
                        ll -= weight_of(slot) * (a <= c ? 0.5 * z * z : c * (a - 0.5 * c));
                    }
                } else {
                    for (const auto& [slot, m] : tot) {
                        const double o = static_cast<double>(counts.node[slot]);
                        ll += weight_of(slot) * mix_diff(o, scale * (lambda * m + mu), mu);
                    }
                }
                // Outside the window the candidate is not scored: the reads say how much sequence is
                // present, and a pair that contradicts that by more than the estimate's own precision
                // is not a candidate whatever its composition.
                if (win_target > 0.0 &&
                    std::fabs(mult_in - win_target) > options.mass_window * win_target) {
                    emis[bi][x * kn + y] = -1e18;
                    explained[bi][x * kn + y] = 0.0;
                    detected[bi][x * kn + y] = 1.0;
                    continue;
                }
                double mass_ll = 0.0;
                if (mass_sd > 0.0) {
                    const double z = (mult_in - mass_target) / mass_sd;
                    mass_ll = -0.5 * options.mass_weight * z * z;
                }
                // Coverage says how much sequence is here; the markers say which alleles it is. The
                // two are good at different halves: the marker emission places one haplotype almost exactly
                // and puts the whole error in the other, while the coverage emission gets the TOTAL
                // nearly right and splits the error across both. So the markers keep the choice and
                // coverage only constrains the sum.
                //
                // This is not the earlier mass term, which failed. That took its target from marker
                // mass -- the same evidence it was correcting -- and enforced it hard enough to pick
                // badly among pairs of the right length. This target is alignment-derived, independent,
                // and accurate to about 1%, so it can be applied at its own precision.
                double cov_ll = 0.0;
                if (cov_target > 0.0 && cov_sd > 0.0) {
                    const double bp_pair =
                        static_cast<double>(kept[bi][x] < blocks[bi].allele_bp.size()
                                                ? blocks[bi].allele_bp[kept[bi][x]] : 0) +
                        static_cast<double>(kept[bi][y] < blocks[bi].allele_bp.size()
                                                ? blocks[bi].allele_bp[kept[bi][y]] : 0);
                    const double z = (bp_pair - cov_target) / cov_sd;
                    cov_ll = -0.5 * z * z;
                }
                if (options.compositional) {
                    // Multinomial over the block's marker universe: log L = SUM o_j log p_j, with
                    // p_j the pair's predicted share of the block's total predicted mass. The
                    // multinomial coefficient depends only on the observations, so it is identical for
                    // every candidate and cancels in the comparison.
                    double pred_tot = 0.0;
                    for (const std::uint32_t slot : universe) {
                        const auto it = tot.find(slot);
                        pred_tot += lambda * (it == tot.end() ? 0.0 : it->second) + mu;
                    }
                    double cll = 0.0;
                    if (pred_tot > 0.0) {
                        for (const std::uint32_t slot : universe) {
                            const double o = static_cast<double>(counts.node[slot]);
                            if (o <= 0.0) continue;
                            const auto it = tot.find(slot);
                            const double pj = (lambda * (it == tot.end() ? 0.0 : it->second) + mu) / pred_tot;
                            cll += o * std::log(std::max(1e-300, pj));
                        }
                    }
                    // A Poisson likelihood factorises exactly into a total and a composition:
                    //   L(o; mean) = Poisson(SUM o; SUM mean) x Multinomial(o; mean / SUM mean)
                    // and the two halves deserve different effective sample sizes. Composition is
                    // spread over many correlated markers, so it takes the ESS discount. The total adds
                    // one observation per fragment crossing the block, so it is far sharper and must
                    // NOT be discounted alongside it -- and it is the only thing that separates a
                    // deletion from a homozygote, where both alleles predict the same proportions and
                    // only the magnitude differs.
                    double sll = 0.0;
                    if (scale_neff > 0.0 && pred_tot > 0.0 && obs_universe > 0.0) {
                        const double c = std::max(1.0, obs_universe / scale_neff);
                        const double oo = obs_universe / c;
                        const double pp = pred_tot / c;
                        sll = -pp + oo * std::log(std::max(1e-300, pp)) - std::lgamma(oo + 1.0);
                    }
                    emis[bi][x * kn + y] = rho * cll + options.scale_weight * sll;
                } else {
                    // Adjacency term, same NB against the same depth, over the pair's summed edge
                    // multiplicities. Discounted by the SAME rho as the nodes: the two are counted
                    // from one set of reads, so they do not carry independent evidence and must not
                    // be added at full weight -- measured strictly worse at low depth when they were.
                    double ell = 0.0;
                    if (options.edge_weight > 0.0) {
                        std::unordered_map<std::uint32_t, std::uint32_t> etot;
                        for (const auto& [slot, mult] : panel.by_block[bi][kept[bi][x]].edges)
                            etot[slot] += mult;
                        for (const auto& [slot, mult] : panel.by_block[bi][kept[bi][y]].edges)
                            etot[slot] += mult;
                        for (const auto& [slot, m] : etot) {
                            const double o = slot < counts.edge.size()
                                                 ? static_cast<double>(counts.edge[slot]) : 0.0;
                            ell += log_nb(o, lambda * m + mu, phi) - log_nb(o, mu, phi);
                        }
                    }
                    emis[bi][x * kn + y] = (options.robust_c > 0.0 ? 0.0 : baseline)
                                           + rho * ll + mass_ll + cov_ll
                                           + options.edge_weight * rho * ell;
                }
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
            int ties = 0;
            for (std::size_t x = 0; x < kn; ++x) {
                for (std::size_t y = x; y < kn; ++y) {
                    if (emis[bi][x * kn + y] > best - 1e-9) ++ties;
                }
            }
            calls[bi].truth_emission_rank = rank;
            calls[bi].truth_emission_delta = tv - best;
            calls[bi].truth_emission_ties = ties;
            calls[bi].n_scored_alleles = static_cast<int>(kn);
        }
    }

    // ---- linkage constraint -------------------------------------------------------------------
    // A floor per block, below which a state is unreachable by the chain. Computed from the
    // UNCONSTRAINED emission and applied only inside emission_for, so every diagnostic above still
    // describes the emission as scored -- otherwise a constrained run could not be compared with an
    // unconstrained one.
    std::vector<double> emis_floor(nb, -std::numeric_limits<double>::infinity());
    if (std::isfinite(options.max_linkage_emission_loss)) {
        for (std::size_t bi = 0; bi < nb; ++bi) {
            if (emis[bi].empty()) continue;
            const double best = *std::max_element(emis[bi].begin(), emis[bi].end());
            emis_floor[bi] = best - options.max_linkage_emission_loss;
        }
    }

    // ---- diploid Li-Stephens forward-backward over haplotype pairs ----
    // Transitions factorize (each haplotype switches independently), so a step is O(n^2) via row and
    // column sums rather than O(n^4).
    auto emission_for = [&](std::size_t bi, std::size_t i, std::size_t j) -> double {
        // EVERY exit goes through the floor, the background fallbacks included. A pruned pair takes
        // the all-background likelihood, and if that is left unguarded the constraint has a hole
        // exactly where the block is largest: at a block above the allele cap, linkage could still
        // reach a pruned pair while every scored pair was being held to within tau. At tau = inf the
        // floor is -inf and this is the identity, so default behaviour is untouched.
        // -1e18 rather than -inf, matching the mass window: the recursion multiplies these and a
        // true -inf would poison sums that still have to normalize.
        const auto guard = [&](double v) { return v < emis_floor[bi] ? -1e18 : v; };
        // At finite tau an UNSCORED pair is excluded outright. Its background likelihood is a
        // constant, not an emission, so "within tau of the block optimum" is not a statement about
        // it -- and if the constant happens to sit above the floor, linkage could reach a pair the
        // emission never ranked while every scored pair was being held to tau.
        const bool restrict_to_scored = std::isfinite(options.max_linkage_emission_loss);
        const int a = allele_of[bi][i];
        const int b = allele_of[bi][j];
        if (a < 0 || b < 0 || kept[bi].empty())
            return restrict_to_scored ? -1e18 : guard(baseline_of[bi]);
        const auto xi = std::lower_bound(kept[bi].begin(), kept[bi].end(), static_cast<std::uint32_t>(a));
        const auto yi = std::lower_bound(kept[bi].begin(), kept[bi].end(), static_cast<std::uint32_t>(b));
        // A pruned allele falls back to the all-background likelihood rather than -inf. Zeroing it
        // would remove that haplotype from the forward recursion for the rest of the chain, so one
        // weak block could silently veto the correct haplotype everywhere.
        if (xi == kept[bi].end() || *xi != static_cast<std::uint32_t>(a))
            return restrict_to_scored ? -1e18 : guard(baseline_of[bi]);
        if (yi == kept[bi].end() || *yi != static_cast<std::uint32_t>(b))
            return restrict_to_scored ? -1e18 : guard(baseline_of[bi]);
        const std::size_t kn = kept[bi].size();
        const double v = emis[bi][static_cast<std::size_t>(xi - kept[bi].begin()) * kn +
                                  static_cast<std::size_t>(yi - kept[bi].begin())];
        return guard(v);
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

        // The block's emission optimum exists whether or not the CALLED pair survived pruning, so
        // it is computed first and unconditionally. Nesting it inside the called-pair lookup left it
        // unreported at exactly the blocks where the call is a pruned pair -- the diagnostic case.
        if (!kept[bi].empty()) {
            const std::size_t kn = kept[bi].size();
            double best = -std::numeric_limits<double>::infinity();
            for (std::size_t x = 0; x < kn; ++x)
                for (std::size_t y = x; y < kn; ++y) best = std::max(best, emis[bi][x * kn + y]);
            int ties = 0, ba = -1, bb = -1;
            for (std::size_t x = 0; x < kn; ++x) {
                for (std::size_t y = x; y < kn; ++y) {
                    if (emis[bi][x * kn + y] > best - 1e-9) {
                        ++ties;
                        if (ba < 0) { ba = static_cast<int>(kept[bi][x]); bb = static_cast<int>(kept[bi][y]); }
                    }
                }
            }
            c.local_best_a = ba; c.local_best_b = bb; c.local_best_ties = ties;

            const auto xi = std::find(kept[bi].begin(), kept[bi].end(), static_cast<std::uint32_t>(c.allele1));
            const auto yi = std::find(kept[bi].begin(), kept[bi].end(), static_cast<std::uint32_t>(c.allele2));
            if (xi != kept[bi].end() && yi != kept[bi].end()) {
                const std::size_t off = static_cast<std::size_t>(xi - kept[bi].begin()) * kn +
                                        static_cast<std::size_t>(yi - kept[bi].begin());
                c.explained = explained[bi][off];
                c.detected = detected[bi][off];
                const double cv = emis[bi][off];
                int rank = 1;
                for (std::size_t x = 0; x < kn; ++x)
                    for (std::size_t y = x; y < kn; ++y)
                        if (emis[bi][x * kn + y] > cv + 1e-9) ++rank;
                c.called_emission_rank = rank;
                c.called_emission_delta = cv - best;
            }
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
        // which at a tandem array can be a whole repeat unit away; the mass estimate is continuous
        // and does not have to round to a panel member.
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
            // The BYPASS allele is excluded: a haplotype that does not traverse the block carries no
            // sequence, so having no markers is the correct description of it rather than a failure of
            // evidence, and the machinery for that case (predicted count zero, `detected`, mass) is
            // already in place. What matters here is a real allele with real sequence that confinement
            // has stripped to nothing -- then the model genuinely cannot see it.
            for (std::size_t ai = 0; ai < panel.by_block[bi].size(); ++ai) {
                if (blocks[bi].bypass_allele >= 0 &&
                    ai == static_cast<std::size_t>(blocks[bi].bypass_allele)) continue;
                if (panel.by_block[bi][ai].nodes.empty()) ++c.alleles_without_markers;
            }
            c.scale_only = c.alleles_without_markers > 0 && c.allele1 == c.allele2 &&
                           c.allele1 < panel.by_block[bi].size() &&
                           !panel.by_block[bi][c.allele1].nodes.empty();
            if (coverage != nullptr && bi < coverage->target_bp.size()) c.cov_bp = coverage->target_bp[bi];
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
        // At a locus where region uniqueness leaves only the flanks with markers, the interior blocks
        // are still called correctly from linkage alone, so treating "no local markers" as a no-call
        // throws away correct calls. Such blocks are marked LINKED and still face the same
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
        // A homozygous call at a block holding an invisible allele is decided by absolute depth alone,
        // and measurement says that is not safe here: on the paralogous fixture the same five markers
        // read 1.78 x lambda for a true heterozygote and 2.33 x lambda for a true homozygote, against
        // predictions of 1.0 and 2.0 -- the classes are half as far apart as the model believes and the
        // heterozygote falls nearer the homozygous prediction. Reported rather than guessed.
        else if (c.scale_only) { c.filter = "SCALEONLY"; ++nocall; }
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
    const std::vector<std::string>& haplotype_names,
    bool probe_target) {

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
         "\tcalled_bp\tmass_bp\tmass_bp_sd\tcov_bp\tno_marker_alleles\tmax_copies\tblock_class"
         "\tevidence\tfilter\ttruth1\ttruth2\ttruth_rank\ttruth_delta\ttruth_ties"
          "\tn_scored_alleles\tcalled_emission_rank\tcalled_emission_delta"
          "\tlocal_best_a\tlocal_best_b\tlocal_best_ties\trank_target\tinfluencers\n";
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
          << c.cov_bp << '\t' << c.alleles_without_markers << '\t' << c.max_copies << '\t' << (c.is_array ? "array" : "simple") << '\t'
          << c.evidence << '\t' << c.filter << '\t'
          << c.truth_allele1 << '\t' << c.truth_allele2 << '\t' << c.truth_emission_rank << '\t'
          << c.truth_emission_delta << '\t' << c.truth_emission_ties << '\t'
          << c.n_scored_alleles << '\t' << c.called_emission_rank << '\t'
          << c.called_emission_delta << '\t' << c.local_best_a << '\t' << c.local_best_b
          << '\t' << c.local_best_ties << '\t' << (probe_target ? "probe" : "truth") << '\t';
        for (std::size_t k = 0; k < c.influencers.size(); ++k) {
            if (k) g << ',';
            g << c.influencers[k];
        }
        if (c.influencers.empty()) g << '.';
        g << '\n';
    }
}

} // namespace panvar
