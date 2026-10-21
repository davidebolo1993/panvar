// M12-A (DIAGNOSTIC): repeated-unit placement-multiplicity and exposure consistency fixtures.
//
// Haplotypes H_k = flankL + unit^k + flankR with identical tandem units and explicit boundaries. Fragments are generated
// deterministically at known starts. Each gate compares the production scorer (score_diplotypes) and the exhaustive oracle
// against an INDEPENDENT calculation written here: physical origins counted by brute-force string comparison, the
// per-origin mass of an exact in-band placement (orientation prior x insert probability; band 0 gives emission 1), the
// analytic exposure sum_j P(j) (L - j + 1)^+, and the Poisson-process score
//   sum_f log((1-eta) lambda (M_a,f + M_b,f) + eta P_bg(f)) - (1-eta) lambda (E_a + E_b).
// This is a consistency test of the declared objective, not a new objective: legitimate copy multiplicity is kept.
//
// G1  origin multiplicity on H_k (k = 0..5) for interior, junction-spanning and boundary fragments = brute force = analytic
// G2  scorer masses = oracle masses = independent n * m0; exposure = analytic; E(H_{k+1}) - E(H_k) = U exactly; the
//     scorer's fragment and exposure terms = the independent formula, for homozygous and heterozygous candidates
// G3  noise-free tiling at the declared intensity: the score peaks at the true total copy number; allocations with the same
//     total copies tie exactly (identical units make allocation unidentifiable, which is correct); the independent score
//     reproduces every candidate; the rate sensitivity of the argmax is reported (not gated)
// G4  equivalent representations: different allele paths spelling the same bytes intern to one haplotype, score identically
// G5  gapped model: an indel-ambiguous mate is ONE origin with several latent alignments (summed, not multiplied); the
//     same mate on two tandem copies is two origins; scorer masses = oracle masses on gapped junction/interior fragments

#include "panvar/fragment_likelihood.hpp"
#include "panvar/genotype_model.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr double kNegInf = -std::numeric_limits<double>::infinity();
std::size_t g_checks = 0, g_failed = 0;

void check(bool ok, const std::string& what) {
    ++g_checks;
    if (!ok) {
        ++g_failed;
        std::cout << "FAIL " << what << '\n';
    }
}

bool close(double a, double b, double tol) {
    if (a == kNegInf || b == kNegInf) return a == b;
    return std::fabs(a - b) <= tol * std::max(1.0, std::fabs(b));
}

std::string revcomp(const std::string& s) {
    std::string out(s.size(), 'N');
    for (std::size_t i = 0; i < s.size(); ++i) {
        const char b = s[s.size() - 1 - i];
        out[i] = b == 'A' ? 'T' : b == 'C' ? 'G' : b == 'G' ? 'C' : 'A';
    }
    return out;
}

std::string random_sequence(std::mt19937& g, std::size_t n) {
    std::string s(n, 'A');
    for (char& c : s) c = "ACGT"[g() & 3];
    return s;
}

panvar::SequencingRead rd(const std::string& name, const std::string& bases) { return {name, bases, std::string(bases.size(), 'I')}; }

constexpr std::size_t R = 20, FLANK = 60, U = 80, INSERT = 40;

panvar::LibraryModel library(std::size_t max_gap = 0) {
    panvar::LibraryModel m;
    m.insert_mean_bp = 40.0;
    m.insert_sd_bp = 3.0;
    m.support_sigmas = 3.0;
    m.substitution_error = 0.01;
    m.outlier_weight = 0.001;
    m.fragment_start_intensity = 1.0 / (1.0 - m.outlier_weight);   // (1 - eta) lambda = one template per base per copy
    m.max_gap_bp = max_gap;
    m.provenance = "fixture: m12a";
    return m;
}

struct Frag {
    panvar::ReadPair pair;
    std::size_t start = 0;
    bool first_forward = true;
    std::string cls;   // origin class on the generating haplotype: boundary | junction | interior
};

// A template [s, s + INSERT) on H in one library orientation.
Frag make_fragment(const std::string& h, std::size_t s, bool first_forward, std::size_t copies) {
    const std::string left = h.substr(s, R), right = revcomp(h.substr(s + INSERT - R, R));
    Frag f;
    f.start = s;
    f.first_forward = first_forward;
    f.pair = first_forward ? panvar::ReadPair{rd("f/1", left), rd("f/2", right)} : panvar::ReadPair{rd("f/1", right), rd("f/2", left)};
    const std::size_t array_begin = FLANK, array_end = FLANK + copies * U, e = s + INSERT;
    if (s < array_begin || e > array_end) f.cls = "boundary";
    else if ((s - array_begin) / U != (e - 1 - array_begin) / U) f.cls = "junction";
    else f.cls = "interior";
    return f;
}

// Independent physical-origin count: brute force over every (left endpoint, insert) in the support and both orientations.
std::size_t brute_origins(const std::string& h, const panvar::ReadPair& p, const panvar::InsertSupport& sup) {
    std::size_t n = 0;
    for (std::size_t x = 0; x + R <= h.size(); ++x) {
        for (std::size_t j = std::max(sup.min_bp, R); j <= sup.max_bp && x + j <= h.size(); ++j) {
            const std::string l = h.substr(x, R), r = revcomp(h.substr(x + j - R, R));
            if (l == p.first.bases && r == p.second.bases) ++n;    // first mate forward at the left endpoint
            if (l == p.second.bases && r == p.first.bases) ++n;    // second mate forward at the left endpoint
        }
    }
    return n;
}

std::size_t analytic_origins(const std::string& cls, std::size_t k) {
    if (cls == "interior") return k;
    if (cls == "junction") return k >= 2 ? k - 1 : 0;
    return k >= 1 ? 1 : 0;   // boundary: the template carries flank sequence, so one place at most (fragments < U)
}

double analytic_exposure(std::size_t length, const panvar::DiscreteInsertDistribution& ins) {
    double e = 0.0;
    for (std::size_t o = 0; o < ins.log_probability_by_offset.size(); ++o) {
        const std::size_t j = ins.support.min_bp + o;
        if (length + 1 > j) e += std::exp(ins.log_probability_by_offset[o]) * static_cast<double>(length + 1 - j);
    }
    return e;
}

panvar::FragmentScoringOptions options(std::size_t band) {
    panvar::FragmentScoringOptions o;
    o.mismatch_band = band;
    o.minimum_piece_bp = 1;
    o.threads = 1;
    return o;
}

struct World {
    std::string flank_l, unit, flank_r;
    std::string hap(std::size_t k) const {
        std::string h = flank_l;
        for (std::size_t i = 0; i < k; ++i) h += unit;
        return h + flank_r;
    }
};

}  // namespace

int main() {
    std::mt19937 g(20261021);
    World w{random_sequence(g, FLANK), random_sequence(g, U), random_sequence(g, FLANK)};
    const panvar::LibraryModel lib = library(0);
    const panvar::DiscreteInsertDistribution ins = panvar::build_library_insert_distribution(lib);
    const double p_insert = std::exp(ins.log_probability_by_offset[INSERT - ins.support.min_bp]);
    const double log_bg_independent = -2.0 * static_cast<double>(R) * std::log(4.0);
    const double rate = (1.0 - lib.outlier_weight) * lib.fragment_start_intensity;
    std::cout << std::setprecision(10) << "fixture: R=" << R << " flank=" << FLANK << " unit=" << U << " insert=" << INSERT
              << " support=[" << ins.support.min_bp << "," << ins.support.max_bp << "] (1-eta)lambda=" << rate << '\n';

    // ---- G1: origin multiplicity --------------------------------------------------------------------------------------
    const std::size_t kt = 3;
    const std::string truth = w.hap(kt);
    std::vector<Frag> frags;
    for (std::size_t s = 0; s + INSERT <= truth.size(); ++s) frags.push_back(make_fragment(truth, s, s % 2 == 0, kt));
    std::map<std::string, std::size_t> by_class;
    for (const Frag& f : frags) ++by_class[f.cls];
    std::size_t g1 = 0;
    for (std::size_t k = 0; k <= 5; ++k) {
        const std::string h = w.hap(k);
        for (const Frag& f : frags) {
            const panvar::InsertSupport sup = panvar::fragment_insert_support(f.pair, lib);
            const std::size_t oracle = panvar::enumerate_placements_exhaustive(h, f.pair, sup, lib, 0).size();
            const std::size_t brute = brute_origins(h, f.pair, sup);
            const std::size_t analytic = analytic_origins(f.cls, k);
            check(oracle == brute && brute == analytic, "G1 origins k=" + std::to_string(k) + " start=" + std::to_string(f.start) + " " + f.cls +
                                                            " oracle=" + std::to_string(oracle) + " brute=" + std::to_string(brute) +
                                                            " analytic=" + std::to_string(analytic));
            ++g1;
        }
    }
    std::cout << "G1 origins: " << g1 << " (fragment, haplotype) comparisons; classes interior=" << by_class["interior"]
              << " junction=" << by_class["junction"] << " boundary=" << by_class["boundary"] << '\n';

    // ---- G2 / G3: scorer vs oracle vs independent; copy-number balance under a noise-free tiling --------------------
    // The observed data: one template per start on EACH of the two truth copies (homozygous H_3), i.e. exactly the
    // declared intensity, all at insert INSERT. Candidates: every unordered pair (a, b), 0 <= a <= b <= 5.
    std::vector<panvar::ReadPair> reads;
    std::vector<const Frag*> origin_of;
    for (int copy = 0; copy < 2; ++copy) {
        for (const Frag& f : frags) {
            reads.push_back(f.pair);
            origin_of.push_back(&f);
        }
    }
    std::vector<panvar::MosaicPath> paths;
    std::vector<std::string> seqs;
    for (std::size_t k = 0; k <= 5; ++k) {
        paths.push_back(panvar::MosaicPath{{k}});
        seqs.push_back(w.hap(k));
    }
    const panvar::InternedHaplotypes interned = panvar::intern_haplotypes(paths, seqs);
    std::vector<panvar::HaplotypeCandidate> cands;
    for (std::size_t a = 0; a <= 5; ++a) for (std::size_t b = a; b <= 5; ++b) cands.push_back({interned.index_of[a], interned.index_of[b]});
    const panvar::FragmentScoringResult scored = panvar::score_diplotypes(interned, cands, reads, lib, options(0));
    const std::size_t H = scored.unique_haplotypes;

    // independent per-haplotype exposure and per-fragment mass
    std::vector<double> e_ind(6);
    std::size_t g2 = 0;
    for (std::size_t k = 0; k <= 5; ++k) {
        const std::size_t u = interned.index_of[k];
        e_ind[k] = analytic_exposure(seqs[k].size(), ins);
        check(close(scored.exposure[u], e_ind[k], 1e-12), "G2 exposure k=" + std::to_string(k));
        if (k >= 1 && k <= 5 && seqs[k].size() > ins.support.max_bp) {
            check(close(e_ind[k] - e_ind[k - 1], static_cast<double>(U), 1e-12) || seqs[k - 1].size() <= ins.support.max_bp,
                  "G2 one unit adds exactly U to the exposure, k=" + std::to_string(k));
        }
        for (std::size_t f = 0; f < frags.size(); ++f) {
            const panvar::ReadPair& p = frags[f].pair;
            const panvar::InsertSupport sup = panvar::fragment_insert_support(p, lib);
            const auto placements = panvar::enumerate_placements_exhaustive(seqs[k], p, sup, lib, 0);
            const double oracle = panvar::log_placement_mass(placements, p, ins, lib);
            const std::size_t n = brute_origins(seqs[k], p, sup);
            const double independent = n == 0 ? kNegInf : std::log(static_cast<double>(n)) + std::log(0.5) + std::log(p_insert);
            check(close(scored.log_mass[f * H + u], oracle, 1e-9) && close(oracle, independent, 1e-9),
                  "G2 mass k=" + std::to_string(k) + " start=" + std::to_string(frags[f].start));
            ++g2;
        }
    }
    check(close(panvar::log_background_probability(reads[0]), log_bg_independent, 1e-12), "G2 background = (1/4)^(2R)");
    // independent score of every candidate
    std::map<std::pair<std::size_t, std::size_t>, double> score_of, frag_of, expo_of;
    std::size_t ci = 0;
    for (std::size_t a = 0; a <= 5; ++a) {
        for (std::size_t b = a; b <= 5; ++b, ++ci) {
            double fragment = 0.0;
            for (std::size_t f = 0; f < reads.size(); ++f) {
                const Frag& fr = *origin_of[f];
                const double m = 0.5 * p_insert * static_cast<double>(analytic_origins(fr.cls, a) + analytic_origins(fr.cls, b));
                fragment += std::log(rate * m + lib.outlier_weight * std::exp(log_bg_independent));
            }
            const double exposure = -rate * (e_ind[a] + e_ind[b]);
            const auto& sc = scored.scores[ci];
            check(close(sc.fragment_term, fragment, 1e-9), "G2 fragment term (" + std::to_string(a) + "," + std::to_string(b) + ")");
            check(close(sc.exposure_term, exposure, 1e-12), "G2 exposure term (" + std::to_string(a) + "," + std::to_string(b) + ")");
            score_of[{a, b}] = sc.fragment_term + sc.exposure_term;
            frag_of[{a, b}] = sc.fragment_term;
            expo_of[{a, b}] = sc.exposure_term;
        }
    }
    std::cout << "G2 masses: " << g2 << " (fragment, haplotype) comparisons; terms checked for " << ci << " candidates\n";

    // G3: the balance. Best total copy number must be the truth's (2 * kt = 6); equal-total allocations tie exactly.
    std::map<std::size_t, double> best_by_total;
    for (const auto& [ab, s] : score_of) {
        const std::size_t t = ab.first + ab.second;
        if (!best_by_total.count(t) || s > best_by_total[t]) best_by_total[t] = s;
    }
    std::size_t argmax_total = 0;
    for (const auto& [t, s] : best_by_total) if (s > best_by_total[argmax_total]) argmax_total = t;
    check(argmax_total == 2 * kt, "G3 the score peaks at the true total copy number (" + std::to_string(argmax_total) + ")");
    check(close(score_of[{3, 3}], score_of[{2, 4}], 1e-12) && close(score_of[{3, 3}], score_of[{1, 5}], 1e-12),
          "G3 identical units: equal-total allocations tie exactly");
    std::cout << "G3 balance (fragment + exposure, relative to (3,3)); total copies: ";
    for (const auto& [t, s] : best_by_total) std::cout << t << ":" << std::fixed << std::setprecision(3) << s - score_of[{3, 3}] << ' ';
    std::cout << std::setprecision(10) << '\n';
    std::cout << "    one copy more (3,4) vs (3,3): fragment " << frag_of[{3, 4}] - frag_of[{3, 3}] << ", exposure "
              << expo_of[{3, 4}] - expo_of[{3, 3}] << "; one copy fewer (2,3): fragment " << frag_of[{2, 3}] - frag_of[{3, 3}]
              << ", exposure " << expo_of[{2, 3}] - expo_of[{3, 3}] << '\n';
    // rate sensitivity (reported, not gated): the argmax under a mis-declared intensity, same data
    for (double scale : {0.8, 0.9, 1.1, 1.25}) {
        std::size_t best_t = 0;
        double best_s = kNegInf;
        for (std::size_t t = 0; t <= 10; ++t) {
            double s = 0.0;
            const std::size_t a = t / 2, b = t - t / 2;
            if (b > 5) continue;
            for (std::size_t f = 0; f < reads.size(); ++f) {
                const Frag& fr = *origin_of[f];
                const double m = 0.5 * p_insert * static_cast<double>(analytic_origins(fr.cls, a) + analytic_origins(fr.cls, b));
                s += std::log(scale * rate * m + lib.outlier_weight * std::exp(log_bg_independent));
            }
            s -= scale * rate * (e_ind[a] + e_ind[b]);
            if (s > best_s) { best_s = s; best_t = t; }
        }
        std::cout << "    rate x" << scale << ": independent argmax total copies " << best_t << " (truth " << 2 * kt << ")\n";
    }

    // ---- G4: equivalent representations --------------------------------------------------------------------------------
    {
        const std::vector<panvar::MosaicPath> two = {panvar::MosaicPath{{7, 1}}, panvar::MosaicPath{{9, 2, 4}}};
        const panvar::InternedHaplotypes same = panvar::intern_haplotypes(two, {w.hap(3), w.hap(3)});
        check(same.sequences.size() == 1 && same.index_of[0] == same.index_of[1], "G4 two paths spelling the same bytes are one haplotype");
        const auto r1 = panvar::score_diplotypes(same, {{same.index_of[0], same.index_of[1]}}, reads, lib, options(0));
        check(close(r1.scores[0].fragment_term + r1.scores[0].exposure_term, score_of[{3, 3}], 1e-12),
              "G4 the duplicate-path homozygote scores as (3,3)");
    }

    // ---- G5: gapped model -- duplicate latent alignments vs distinct physical starts ----------------------------------
    {
        const panvar::LibraryModel gl = library(3);
        // a haplotype with a homopolymer run inside the unit; the forward mate carries a one-base deletion within the run
        std::string unit = w.unit;
        unit.replace(30, 8, "AAAAAAAA");
        const std::string h1 = w.flank_l + unit + w.flank_r, h2 = w.flank_l + unit + unit + w.flank_r;
        const std::size_t s = FLANK + 22;   // the mate starts outside the run and reads into it
        std::string left = h1.substr(s, R + 1);
        left.erase(12, 1);                  // delete one A inside the run (run occupies offsets 8..15 of the mate)
        const std::string right = revcomp(h1.substr(s + INSERT - R, R));
        const panvar::ReadPair p{rd("g/1", left), rd("g/2", right)};
        const panvar::InsertSupport sup = panvar::fragment_insert_support(p, gl);
        const auto on1 = panvar::enumerate_placements_exhaustive(h1, p, sup, gl, 2);
        const auto on2 = panvar::enumerate_placements_exhaustive(h2, p, sup, gl, 2);
        std::set<std::pair<std::size_t, std::size_t>> origins1, origins2;
        std::size_t latent = 0;
        for (const auto& o : on1) { origins1.insert({o.forward_start, o.insert_bp}); latent = std::max(latent, o.latent_alignments); }
        for (const auto& o : on2) origins2.insert({o.forward_start, o.insert_bp});
        // the true origin must be present on one copy, and present once per tandem copy on two
        const bool true_on1 = origins1.count({s, INSERT}) == 1;
        std::size_t true_like2 = 0;
        for (const auto& o : origins2) if (o.second == INSERT && (o.first == s || o.first == s + U)) ++true_like2;
        check(true_on1 && latent > 1, "G5 an indel-ambiguous mate is one origin with " + std::to_string(latent) + " latent alignments");
        check(true_like2 == 2, "G5 the same mate on two tandem copies is two origins");
        check(on2.size() == 2 * on1.size(), "G5 two tandem copies double the origin set (" + std::to_string(on1.size()) + " -> " +
                                                std::to_string(on2.size()) + ")");
        // scorer vs oracle masses under the gapped model, interior and junction fragments on H_1..H_3 (one substitution each)
        std::vector<panvar::ReadPair> greads;
        std::vector<std::string> gseq = {w.hap(1), w.hap(2), w.hap(3)};
        for (std::size_t st = FLANK - 10; st + INSERT <= FLANK + 2 * U + 10; st += 7) {
            Frag f = make_fragment(w.hap(2), st, st % 2 == 0, 2);
            std::string b = f.pair.first.bases;
            b[5] = b[5] == 'A' ? 'C' : 'A';
            f.pair.first.bases = b;
            greads.push_back(f.pair);
        }
        const panvar::InternedHaplotypes gi = panvar::intern_haplotypes({panvar::MosaicPath{{1}}, panvar::MosaicPath{{2}}, panvar::MosaicPath{{3}}}, gseq);
        const auto gs = panvar::score_diplotypes(gi, {{0, 1}, {1, 2}}, greads, gl, options(2));
        std::size_t compared = 0;
        for (std::size_t f = 0; f < greads.size(); ++f) {
            for (std::size_t k = 0; k < 3; ++k) {
                const auto pl = panvar::enumerate_placements_exhaustive(gseq[k], greads[f], panvar::fragment_insert_support(greads[f], gl), gl, 2);
                const double oracle = panvar::log_placement_mass(pl, greads[f], panvar::build_library_insert_distribution(gl), gl);
                check(close(gs.log_mass[f * gs.unique_haplotypes + gi.index_of[k]], oracle, 1e-9), "G5 gapped scorer mass = oracle mass");
                ++compared;
            }
        }
        std::cout << "G5 gapped: one-copy origins " << on1.size() << " (max latent alignments " << latent << "), two-copy origins "
                  << on2.size() << "; " << compared << " gapped mass comparisons\n";
    }

    std::cout << (g_failed == 0 ? "ALL GATES PASS" : "GATES FAILED") << ": " << g_checks - g_failed << " of " << g_checks << " checks\n";
    return g_failed == 0 ? 0 : 1;
}
