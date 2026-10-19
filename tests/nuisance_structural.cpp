// M6-C: nuisance-marginalised structural likelihood, SYNTHETIC FEASIBILITY ONLY (diagnostic; not a caller, not
// sequence reconstruction, never real-locus fitting).
//
//   P(reads | G) = sum_theta P(theta | G) P(reads | G, theta)
//
// G: a synthetic structural scaffold -- haplotype A = F_A + nA unit copies + G_A, haplotype B = F_B + nB unit copies +
// G_B, each copy carrying a declared base allele at two predeclared biallelic unit sites. theta: ONE assignment of the
// allele at every (copy, site) shared by the whole read set (2^(2(nA+nB)) states, at most 1,024; enumerated completely).
// P(theta | G): independent per (copy, site) flips from the scaffold's base allele with declared probability p (normalised:
// sums to one). P(reads | G, theta): the EXISTING paired-fragment likelihood (score_diplotypes: placement alternatives
// summed, exposure, outlier/background mixture) of the complete diplotype. theta is summed (long double log-sum-exp),
// never maximised; the count score is not used; no copying prior is applied to structures.
//
//   nuisance_structural gates
//   nuisance_structural run <fixture> <seed> <out prefix>
#include "panvar/fragment_likelihood.hpp"
#include "panvar/genotype_model.hpp"

#include <sys/resource.h>

#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr double kNuisanceFlip = 0.05;   // declared nuisance prior: per (copy, site) flip probability
constexpr int kReadLength = 100;

panvar::LibraryModel library() {
    panvar::LibraryModel lib;
    lib.insert_mean_bp = 350;
    lib.insert_sd_bp = 35;
    lib.support_sigmas = 4;
    lib.substitution_error = 0.002;
    lib.outlier_weight = 0.001;
    lib.fragment_start_intensity = 20.0 / (4.0 * kReadLength);   // diploid depth 20, lambda = depth / (4 L)
    lib.provenance = "M6-C synthetic declared";
    return lib;
}

panvar::FragmentScoringOptions scoring(std::size_t threads) {
    panvar::FragmentScoringOptions o;
    o.threads = threads;
    o.mismatch_band = 6;
    o.minimum_piece_bp = 12;
    return o;
}

std::string rnd(std::mt19937_64& rng, std::size_t n) {
    std::string s(n, 'A');
    for (char& c : s) c = "ACGT"[rng() % 4];
    return s;
}
std::string revcomp(const std::string& x) {
    std::string r(x.rbegin(), x.rend());
    for (char& c : r) c = c == 'A' ? 'T' : c == 'C' ? 'G' : c == 'G' ? 'C' : 'A';
    return r;
}
char alt_of(char c) { return c == 'A' ? 'C' : c == 'C' ? 'G' : c == 'G' ? 'T' : 'A'; }

struct World {
    std::string FA, GA, FB, GB, unit;
    std::array<std::size_t, 2> site{};
    std::array<char, 2> ref{}, alt{};
};

World make_world(std::uint64_t seed, std::size_t unit_bp, bool synonymous_site2 = false) {
    std::mt19937_64 rng(seed * 7919 + unit_bp);
    World w;
    w.FA = rnd(rng, 1200);
    w.GA = rnd(rng, 1200);
    w.FB = rnd(rng, 1200);
    w.GB = rnd(rng, 1200);
    w.unit = rnd(rng, unit_bp);
    w.site = {unit_bp / 3, (2 * unit_bp) / 3};   // predeclared positions (never chosen from data)
    for (int s = 0; s < 2; ++s) {
        w.ref[s] = w.unit[w.site[s]];
        w.alt[s] = alt_of(w.ref[s]);
    }
    if (synonymous_site2) w.alt[1] = w.ref[1];   // a label-only site: both alleles spell identical bytes (gate)
    return w;
}

struct Scaffold {
    int nA = 0, nB = 0;
    std::vector<std::array<int, 2>> base;   // per copy (A copies first), allele index 0 ref / 1 alt
    std::string label() const { return std::to_string(nA) + "+" + std::to_string(nB); }
};

std::string spell_hap(const World& w, const Scaffold& g, const std::vector<std::array<int, 2>>& alleles, bool hapA) {
    std::string s = hapA ? w.FA : w.FB;
    const int from = hapA ? 0 : g.nA, to = hapA ? g.nA : g.nA + g.nB;
    for (int c = from; c < to; ++c) {
        std::string u = w.unit;
        for (int k = 0; k < 2; ++k) u[w.site[k]] = alleles[c][k] ? w.alt[k] : w.ref[k];
        s += u;
    }
    return s + (hapA ? w.GA : w.GB);
}

// Reads: fragment starts per haplotype base ~ Poisson(lambda); insert ~ N(350, 35) within +-4 sd; both orientations;
// substitution errors at 0.002. Deterministic per seed.
void simulate(const std::string& hap, double lambda, std::mt19937_64& rng, std::vector<panvar::ReadPair>& out, const std::string& tag) {
    std::poisson_distribution<int> starts(lambda);
    std::normal_distribution<double> insert(350.0, 35.0);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    const auto noisy = [&](std::string r) {
        for (char& c : r) if (u(rng) < 0.002) c = alt_of(c);
        return r;
    };
    for (std::size_t p = 0; p < hap.size(); ++p) {
        for (int k = starts(rng); k > 0; --k) {
            const int len = static_cast<int>(std::lround(std::max(210.0, std::min(490.0, insert(rng)))));
            if (p + static_cast<std::size_t>(len) > hap.size()) continue;
            const std::string frag = hap.substr(p, len);
            std::string m1 = frag.substr(0, kReadLength), m2 = revcomp(frag.substr(len - kReadLength));
            if (u(rng) < 0.5) std::swap(m1, m2);
            const std::string name = tag + std::to_string(out.size());
            out.push_back({{name + "/1", noisy(m1), std::string(kReadLength, 'I')}, {name + "/2", noisy(m2), std::string(kReadLength, 'I')}});
        }
    }
}

long double logsumexp(const std::vector<long double>& v) {
    long double m = -INFINITY;
    for (auto x : v) m = std::max(m, x);
    if (!std::isfinite(static_cast<double>(m))) return m;
    long double s = 0;
    for (auto x : v) s += std::exp(x - m);
    return m + std::log(s);
}

struct Evaluation {
    std::string label;
    std::size_t states = 0, unique_haplotypes = 0;
    long double marginal = 0, maximum = 0, fixed = 0;   // log P(reads|G): summed theta, best theta (negative control), base theta
    long double prior_mass_check = 0;                 // sum of P(theta|G) (must be 1)
    double fixed_fragment = 0, fixed_exposure = 0, best_fragment = 0, best_exposure = 0, best_log_prior = 0;
    double top_state_posterior = 0, effective_states = 0;
    int best_flips = 0;
};

// The batched, interned evaluation (the implementation under test).
Evaluation evaluate(const World& w, const Scaffold& g, const std::vector<panvar::ReadPair>& reads, double p, std::size_t threads) {
    const int bits = 2 * (g.nA + g.nB);
    if (bits > 10) throw std::runtime_error("more than 1,024 nuisance states");
    const std::size_t states = std::size_t(1) << bits;
    std::vector<std::string> seqs;
    std::vector<panvar::MosaicPath> paths;
    std::vector<long double> log_prior(states);
    std::vector<int> flips(states);
    std::vector<std::pair<std::size_t, std::size_t>> pair_of(states);
    long double mass = 0;
    for (std::size_t t = 0; t < states; ++t) {
        std::vector<std::array<int, 2>> al = g.base;
        int f = 0;
        long double lp = 0;
        for (int b = 0; b < bits; ++b) {
            const bool flip = (t >> b) & 1u;
            if (flip) al[b / 2][b % 2] ^= 1;
            f += flip;
            lp += flip ? std::log(static_cast<long double>(p)) : std::log1p(static_cast<long double>(-p));
        }
        if (p == 0.0 && f > 0) lp = -INFINITY;
        log_prior[t] = lp;
        flips[t] = f;
        mass += std::exp(lp);
        seqs.push_back(spell_hap(w, g, al, true));
        seqs.push_back(spell_hap(w, g, al, false));
        paths.push_back(panvar::MosaicPath{{2 * t}});
        paths.push_back(panvar::MosaicPath{{2 * t + 1}});
        pair_of[t] = {2 * t, 2 * t + 1};
    }
    const panvar::InternedHaplotypes in = panvar::intern_haplotypes(paths, seqs);
    std::vector<panvar::HaplotypeCandidate> cand;
    std::map<std::pair<std::size_t, std::size_t>, std::size_t> slot;
    std::vector<std::size_t> slot_of(states);
    for (std::size_t t = 0; t < states; ++t) {
        std::size_t a = in.index_of[pair_of[t].first], b = in.index_of[pair_of[t].second];
        const auto key = std::minmax(a, b);
        auto it = slot.find({key.first, key.second});
        if (it == slot.end()) {
            it = slot.emplace(std::make_pair(key.first, key.second), cand.size()).first;
            cand.push_back({key.first, key.second});
        }
        slot_of[t] = it->second;
    }
    const panvar::FragmentScoringResult r = panvar::score_diplotypes(in, cand, reads, library(), scoring(threads));
    Evaluation e;
    e.label = g.label();
    e.states = states;
    e.unique_haplotypes = in.sequences.size();
    e.prior_mass_check = mass;
    std::vector<long double> terms(states);
    std::size_t best = 0;
    for (std::size_t t = 0; t < states; ++t) {
        terms[t] = log_prior[t] + r.scores[slot_of[t]].log_likelihood;
        if (terms[t] > terms[best]) best = t;
    }
    e.marginal = logsumexp(terms);
    std::size_t argmax_l = 0;
    for (std::size_t t = 0; t < states; ++t) if (r.scores[slot_of[t]].log_likelihood > r.scores[slot_of[argmax_l]].log_likelihood) argmax_l = t;
    e.maximum = r.scores[slot_of[argmax_l]].log_likelihood;   // negative control: best-fitting state, no prior, no sum
    e.fixed = r.scores[slot_of[0]].log_likelihood;            // theta = base content (state 0)
    e.fixed_fragment = r.scores[slot_of[0]].fragment_term;
    e.fixed_exposure = r.scores[slot_of[0]].exposure_term;
    e.best_fragment = r.scores[slot_of[best]].fragment_term;
    e.best_exposure = r.scores[slot_of[best]].exposure_term;
    e.best_log_prior = static_cast<double>(log_prior[best]);
    e.best_flips = flips[best];
    long double h = 0;
    for (std::size_t t = 0; t < states; ++t) {
        const long double q = std::exp(terms[t] - e.marginal);
        h += q * q;
    }
    e.top_state_posterior = static_cast<double>(std::exp(terms[best] - e.marginal));
    e.effective_states = static_cast<double>(1.0L / h);
    return e;
}

// Independent exhaustive oracle: every theta scored ALONE (its own interning, its own call), no deduplication, summed.
long double oracle(const World& w, const Scaffold& g, const std::vector<panvar::ReadPair>& reads, double p) {
    const int bits = 2 * (g.nA + g.nB);
    std::vector<long double> terms;
    for (std::size_t t = 0; t < (std::size_t(1) << bits); ++t) {
        std::vector<std::array<int, 2>> al = g.base;
        long double lp = 0;
        for (int b = 0; b < bits; ++b) {
            const bool flip = (t >> b) & 1u;
            if (flip) al[b / 2][b % 2] ^= 1;
            lp += flip ? std::log(static_cast<long double>(p)) : std::log1p(static_cast<long double>(-p));
        }
        const std::string a = spell_hap(w, g, al, true), b2 = spell_hap(w, g, al, false);
        const panvar::InternedHaplotypes in = panvar::intern_haplotypes({panvar::MosaicPath{{0}}, panvar::MosaicPath{{1}}}, {a, b2});
        const auto r = panvar::score_diplotypes(in, {{in.index_of[0], in.index_of[1]}}, reads, library(), scoring(1));
        terms.push_back(lp + r.scores[0].log_likelihood);
    }
    return logsumexp(terms);
}

Scaffold scaffold(int nA, int nB, int base_allele) {
    Scaffold g;
    g.nA = nA;
    g.nB = nB;
    g.base.assign(nA + nB, {base_allele, base_allele});
    return g;
}
std::vector<Scaffold> catalogue(int base_balanced, int base_unbalanced) {
    std::vector<Scaffold> c;
    for (int total = 3; total <= 5; ++total)
        for (int a = 1; a < total; ++a) c.push_back(scaffold(a, total - a, a == total - a ? base_balanced : base_unbalanced));
    return c;
}

std::size_t g_failed = 0;
std::ostringstream g_rows;
void row(const std::string& t, const std::string& item, double v, const std::string& expect, bool pass) {
    g_rows << t << '\t' << item << '\t' << std::setprecision(15) << v << '\t' << expect << '\t' << (pass ? "PASS" : "FAIL") << '\n';
    g_failed += pass ? 0 : 1;
}
std::vector<panvar::ReadPair> reads_of(const World& w, const Scaffold& truth, std::uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::vector<panvar::ReadPair> r;
    simulate(spell_hap(w, truth, truth.base, true), library().fragment_start_intensity, rng, r, "a");
    simulate(spell_hap(w, truth, truth.base, false), library().fragment_start_intensity, rng, r, "b");
    return r;
}

int gates() {
    const World w = make_world(5, 200);
    Scaffold truth = scaffold(1, 1, 0);
    truth.base[1] = {1, 0};
    const auto reads = reads_of(w, truth, 99);
    const Scaffold g11 = scaffold(1, 1, 0), g21 = scaffold(2, 1, 0);
    // GA: nuisance probabilities sum to one
    for (const Scaffold& g : {g11, g21}) {
        const Evaluation e = evaluate(w, g, reads, kNuisanceFlip, 4);
        row("GA_prior_mass", g.label() + ": sum_theta P(theta|G)", static_cast<double>(e.prior_mass_check), "1 (|.-1| <= 1e-15)", std::fabs(static_cast<double>(e.prior_mass_check) - 1.0) <= 1e-15);
    }
    // GB: nuisance disabled (p = 0) reproduces fixed-sequence scoring
    {
        const Evaluation e = evaluate(w, g11, reads, 0.0, 4);
        row("GB_disabled", "p = 0: marginal - fixed", static_cast<double>(e.marginal - e.fixed), "0", e.marginal == e.fixed);
    }
    // GC: independent exhaustive oracle (each theta scored alone, no dedup)
    for (const Scaffold& g : {g11, g21}) {
        const Evaluation e = evaluate(w, g, reads, kNuisanceFlip, 4);
        const long double o = oracle(w, g, reads, kNuisanceFlip);
        row("GC_oracle", g.label() + ": batched marginal - oracle", static_cast<double>(e.marginal - o), "|.| <= 1e-6", std::fabs(static_cast<double>(e.marginal - o)) <= 1e-6);
    }
    // GD: negative control -- maximisation is detected (differs from the marginal; the oracle sides with the sum)
    {
        const Evaluation e = evaluate(w, g11, reads, kNuisanceFlip, 4);
        const long double o = oracle(w, g11, reads, kNuisanceFlip);
        row("GD_negative_control", "max-theta log L - marginal", static_cast<double>(e.maximum - e.marginal), "!= 0 (substituting max for sum is detectable)",
            std::fabs(static_cast<double>(e.maximum - e.marginal)) > 1e-6 && std::fabs(static_cast<double>(e.maximum - o)) > 1e-6);
    }
    // GE: equivalent sequence states keep their summed prior mass (site 2 synonymous: theta labels collapse to one sequence)
    {
        const World ws = make_world(5, 200, true);
        const auto rs = reads_of(ws, truth, 98);
        const Evaluation e = evaluate(ws, g11, rs, kNuisanceFlip, 4);
        const long double o = oracle(ws, g11, rs, kNuisanceFlip);
        row("GE_equivalent_states", "batched (interned) - oracle with a synonymous site", static_cast<double>(e.marginal - o), "|.| <= 1e-6",
            std::fabs(static_cast<double>(e.marginal - o)) <= 1e-6);
        // the 16 states spell 4 distinct diplotypes: unique haplotypes 2 per side
        row("GE_equivalent_states", "unique haplotypes (synonymous site collapses labels)", static_cast<double>(e.unique_haplotypes), "4", e.unique_haplotypes == 4);
    }
    // GF: whole-haplotype exchange and harmless relabelling are invariant
    {
        World sw = w;
        std::swap(sw.FA, sw.FB);
        std::swap(sw.GA, sw.GB);
        Scaffold g21s = scaffold(1, 2, 0);   // (2,1) with the haplotypes exchanged
        const Evaluation e1 = evaluate(w, g21, reads, kNuisanceFlip, 4), e2 = evaluate(sw, g21s, reads, kNuisanceFlip, 4);
        row("GF_invariance", "whole-haplotype exchange: marginal difference", static_cast<double>(e1.marginal - e2.marginal), "|.| <= 1e-9",
            std::fabs(static_cast<double>(e1.marginal - e2.marginal)) <= 1e-9);
        World rl = w;
        std::swap(rl.site[0], rl.site[1]);
        std::swap(rl.ref[0], rl.ref[1]);
        std::swap(rl.alt[0], rl.alt[1]);
        Scaffold t2 = truth;
        std::swap(t2.base[1][0], t2.base[1][1]);
        const Evaluation e3 = evaluate(rl, g21, reads, kNuisanceFlip, 4);
        row("GF_invariance", "site relabelling: marginal difference", static_cast<double>(e1.marginal - e3.marginal), "|.| <= 1e-9",
            std::fabs(static_cast<double>(e1.marginal - e3.marginal)) <= 1e-9);
    }
    // GG: identical complete observation models tie; pooled-count equality alone is not identity
    {
        const World wl = make_world(6, 1500);   // long unit: no fragment spans a copy
        const auto rl = reads_of(wl, scaffold(3, 1, 0), 77);
        const Evaluation a = evaluate(wl, scaffold(3, 1, 0), rl, kNuisanceFlip, 4), b = evaluate(wl, scaffold(2, 2, 0), rl, kNuisanceFlip, 4);
        row("GG_ties", "long unit (no spanning fragment): 3+1 vs 2+2 marginal", static_cast<double>(a.marginal - b.marginal), "|.| <= 1e-6 (tied)",
            std::fabs(static_cast<double>(a.marginal - b.marginal)) <= 1e-6);
        const World ws = make_world(6, 200);    // short unit: same pooled 31-mer content, spanning fragments exist
        const auto rs = reads_of(ws, scaffold(3, 1, 0), 77);
        const Evaluation c = evaluate(ws, scaffold(3, 1, 0), rs, kNuisanceFlip, 4), d = evaluate(ws, scaffold(2, 2, 0), rs, kNuisanceFlip, 4);
        row("GG_ties", "short unit (same pooled counts, informative adjacency): 3+1 - 2+2", static_cast<double>(c.marginal - d.marginal), "!= 0 (not tied)",
            std::fabs(static_cast<double>(c.marginal - d.marginal)) > 1e-3);
    }
    std::cout << "test\titem\tvalue\texpectation\tresult\n" << g_rows.str() << "nuisance_gates\t" << (g_failed ? "FAIL" : "PASS") << '\n';
    return g_failed ? 1 : 0;
}

// The frozen fixtures (DESIGN.md): unit length, truth scaffold (allocation + true copy alleles), catalogue base contents.
struct Fixture {
    std::string name;
    std::size_t unit_bp;
    Scaffold truth;
    int base_balanced, base_unbalanced;   // catalogue copy alleles for balanced / unbalanced scaffolds
    bool contamination;
};
std::vector<Fixture> fixtures() {
    std::vector<Fixture> f;
    f.push_back({"F1_control_2+2", 600, scaffold(2, 2, 0), 0, 0, false});
    f.push_back({"F2a_missing_variants_invisible_allocation", 600, scaffold(2, 2, 1), 0, 1, false});
    f.push_back({"F2b_missing_variants_informative_allocation", 200, scaffold(2, 2, 1), 0, 1, false});
    f.push_back({"F3_matched_content_true_3+1", 1500, scaffold(3, 1, 0), 0, 0, false});
    f.push_back({"F3_matched_content_true_2+2", 1500, scaffold(2, 2, 0), 0, 0, false});
    f.push_back({"F4_adjacency_true_3+1", 200, scaffold(3, 1, 0), 0, 0, false});
    f.push_back({"F4_adjacency_true_2+2", 200, scaffold(2, 2, 0), 0, 0, false});
    f.push_back({"F6_dosage_true_2+3", 600, scaffold(2, 3, 0), 0, 0, false});
    f.push_back({"F7_contamination_2+2", 200, scaffold(2, 2, 0), 0, 0, true});
    return f;
}

double peak_rss_gib() {
    struct rusage u {};
    getrusage(RUSAGE_SELF, &u);
    return static_cast<double>(u.ru_maxrss) / (1024.0 * 1024.0 * 1024.0);   // bytes on macOS
}

int run(const std::string& name, std::uint64_t seed, const std::string& out) {
    const auto t0 = std::chrono::steady_clock::now();
    Fixture fx;
    bool found = false;
    for (const auto& f : fixtures()) if (f.name == name) { fx = f; found = true; }
    if (!found) throw std::runtime_error("unknown fixture " + name);
    const World w = make_world(seed, fx.unit_bp);
    auto reads = reads_of(w, fx.truth, seed * 1000 + 1);
    std::size_t contaminants = 0;
    if (fx.contamination) {   // 5% fragments from a diverged paralog array, ~2% from unrelated sequence
        std::mt19937_64 rng(seed * 1000 + 2);
        std::string para = rnd(rng, 1200);
        std::string pu = w.unit;
        for (char& c : pu) if (rng() % 100 < 4) c = alt_of(c);
        for (int k = 0; k < 3; ++k) para += pu;
        para += rnd(rng, 1200);
        const std::size_t before = reads.size();
        const double n_true = static_cast<double>(reads.size());
        const std::string off = rnd(rng, 3000);
        simulate(para, 0.05 * n_true / static_cast<double>(para.size()), rng, reads, "p");   // expected 5% of the true fragments
        simulate(off, 0.02 * n_true / static_cast<double>(off.size()), rng, reads, "o");     // expected 2%
        contaminants = reads.size() - before;
    }
    std::ostringstream t;
    t << std::setprecision(12) << "fixture\tseed\tscaffold\tnuisance_states\tunique_haplotypes\tlog_marginal\tlog_fixed\tlog_max_negative_control\t"
         "marginal_minus_best_marginal\tfixed_minus_best_fixed\tfixed_fragment\tfixed_exposure\tbest_state_fragment\tbest_state_exposure\t"
         "best_state_log_prior\tbest_state_flips\ttop_state_posterior\teffective_states\tposterior_structure\n";
    std::vector<Evaluation> ev;
    for (const Scaffold& g : catalogue(fx.base_balanced, fx.base_unbalanced)) ev.push_back(evaluate(w, g, reads, kNuisanceFlip, 4));
    long double bm = -INFINITY, bf = -INFINITY;
    std::vector<long double> marg;
    for (const auto& e : ev) { bm = std::max(bm, e.marginal); bf = std::max(bf, e.fixed); marg.push_back(e.marginal); }
    const long double z = logsumexp(marg);   // uniform over the declared catalogue: structural posterior with NO structural prior
    for (const auto& e : ev) {
        t << name << '\t' << seed << '\t' << e.label << '\t' << e.states << '\t' << e.unique_haplotypes << '\t' << static_cast<double>(e.marginal) << '\t'
          << static_cast<double>(e.fixed) << '\t' << static_cast<double>(e.maximum) << '\t' << static_cast<double>(e.marginal - bm) << '\t'
          << static_cast<double>(e.fixed - bf) << '\t' << e.fixed_fragment << '\t' << e.fixed_exposure << '\t' << e.best_fragment << '\t' << e.best_exposure << '\t'
          << e.best_log_prior << '\t' << e.best_flips << '\t' << e.top_state_posterior << '\t' << e.effective_states << '\t'
          << static_cast<double>(std::exp(e.marginal - z)) << '\n';
    }
    const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::ofstream(out + ".scaffolds.tsv") << t.str();
    std::ofstream(out + ".run.tsv") << "key\tvalue\nfixture\t" << name << "\nseed\t" << seed << "\ntruth\t" << fx.truth.label() << "\ntruth_copy_alleles\t"
                                    << fx.truth.base[0][0] << "\nunit_bp\t" << fx.unit_bp << "\nfragments\t" << reads.size() << "\ncontaminant_fragments\t" << contaminants
                                    << "\nnuisance_flip_prior\t" << kNuisanceFlip << "\nseconds\t" << secs << "\npeak_rss_gib\t" << peak_rss_gib() << '\n';
    std::cout << name << " seed " << seed << ": " << reads.size() << " fragments, " << secs << " s\n";
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const std::string mode = argc > 1 ? argv[1] : "";
        if (mode == "gates") return gates();
        if (mode == "run" && argc == 5) return run(argv[2], std::stoull(argv[3]), argv[4]);
        std::cerr << "usage: nuisance_structural gates | run <fixture> <seed> <out prefix>\n";
        return 2;
    } catch (const std::exception& e) {
        std::cerr << "nuisance_structural: ERROR " << e.what() << '\n';
        return 2;
    }
}
