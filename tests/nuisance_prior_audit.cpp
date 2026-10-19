// M6-D: audit of candidate-centred nuisance priors, SYNTHETIC DIAGNOSTIC ONLY (not a caller, not a proposed biological
// prior, never real-locus fitting).
//
// Each scaffold's absolute sequence states (one allele per (copy, site), 2^(2(nA+nB)) <= 1,024) are enumerated and scored
// ONCE with the existing paired-fragment likelihood (score_diplotypes, as in M6-C); two declared nuisance priors are then
// applied to the same likelihood table:
//   original        independent flips with p = 0.05 relative to the scaffold's base alleles (M6-C)
//   common_content  independent absolute alleles, P(World.alt) = 0.05, P(World.ref) = 0.95 at every (copy, site),
//                   identical for every scaffold; the base annotation does not enter (diagnostic control only)
// Per prior: log marginal = posterior-mean log likelihood - KL(posterior over nuisance states || nuisance prior), reported
// as an exact identity (the best state is not substituted for it).
//
// The world, scaffolds, read simulator and the four fixtures are copied verbatim from tests/nuisance_structural.cpp
// (M6-C), so the reads are M6-C's deterministic reads; the original-prior marginals must reproduce M6-C's recorded ones.
//
//   nuisance_prior_audit gates
//   nuisance_prior_audit run <fixture> <seed> <out prefix>
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

constexpr double kNuisanceFlip = 0.05;   // original prior: per (copy, site) flip probability from the base allele
constexpr double kAltProbability = 0.05; // common-content control: per (copy, site) P(World.alt)
constexpr int kReadLength = 100;

// ---- verbatim from M6-C (tests/nuisance_structural.cpp) -------------------------------------------------------------
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
    w.site = {unit_bp / 3, (2 * unit_bp) / 3};
    for (int s = 0; s < 2; ++s) {
        w.ref[s] = w.unit[w.site[s]];
        w.alt[s] = alt_of(w.ref[s]);
    }
    if (synonymous_site2) w.alt[1] = w.ref[1];
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
std::vector<panvar::ReadPair> reads_of(const World& w, const Scaffold& truth, std::uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::vector<panvar::ReadPair> r;
    simulate(spell_hap(w, truth, truth.base, true), library().fragment_start_intensity, rng, r, "a");
    simulate(spell_hap(w, truth, truth.base, false), library().fragment_start_intensity, rng, r, "b");
    return r;
}
// ---- end verbatim ----------------------------------------------------------------------------------------------------

long double logsumexp(const std::vector<long double>& v) {
    long double m = -INFINITY;
    for (auto x : v) m = std::max(m, x);
    if (!std::isfinite(static_cast<double>(m))) return m;
    long double s = 0;
    for (auto x : v) s += std::exp(x - m);
    return m + std::log(s);
}

// One scaffold's absolute-state likelihood table. States are enumerated as flips from the scaffold's base (M6-C's order),
// but every state records its ABSOLUTE alleles; priors are functions of the absolute alleles (and, for the original
// prior only, of the base annotation).
struct Table {
    Scaffold g;
    std::size_t states = 0, unique_haplotypes = 0, unique_diplotypes = 0;
    std::vector<std::vector<std::array<int, 2>>> alleles;   // absolute allele per (copy, site)
    std::vector<std::size_t> slot;                          // unique-diplotype slot (equivalent states share one)
    std::vector<double> log_likelihood, fragment, exposure;
};

Table build(const World& w, const Scaffold& g, const std::vector<panvar::ReadPair>& reads, std::size_t threads) {
    const int bits = 2 * (g.nA + g.nB);
    if (bits > 10) throw std::runtime_error("more than 1,024 nuisance states");
    Table T;
    T.g = g;
    T.states = std::size_t(1) << bits;
    std::vector<std::string> seqs;
    std::vector<panvar::MosaicPath> paths;
    for (std::size_t t = 0; t < T.states; ++t) {
        std::vector<std::array<int, 2>> al = g.base;
        for (int b = 0; b < bits; ++b) if ((t >> b) & 1u) al[b / 2][b % 2] ^= 1;
        T.alleles.push_back(al);
        seqs.push_back(spell_hap(w, g, al, true));
        seqs.push_back(spell_hap(w, g, al, false));
        paths.push_back(panvar::MosaicPath{{2 * t}});
        paths.push_back(panvar::MosaicPath{{2 * t + 1}});
    }
    const panvar::InternedHaplotypes in = panvar::intern_haplotypes(paths, seqs);
    std::vector<panvar::HaplotypeCandidate> cand;
    std::map<std::pair<std::size_t, std::size_t>, std::size_t> slot;
    for (std::size_t t = 0; t < T.states; ++t) {
        const auto key = std::minmax(in.index_of[2 * t], in.index_of[2 * t + 1]);
        auto it = slot.find({key.first, key.second});
        if (it == slot.end()) {
            it = slot.emplace(std::make_pair(key.first, key.second), cand.size()).first;
            cand.push_back({key.first, key.second});
        }
        T.slot.push_back(it->second);
    }
    const panvar::FragmentScoringResult r = panvar::score_diplotypes(in, cand, reads, library(), scoring(threads));
    T.unique_haplotypes = in.sequences.size();
    T.unique_diplotypes = cand.size();
    for (std::size_t t = 0; t < T.states; ++t) {
        T.log_likelihood.push_back(r.scores[T.slot[t]].log_likelihood);
        T.fragment.push_back(r.scores[T.slot[t]].fragment_term);
        T.exposure.push_back(r.scores[T.slot[t]].exposure_term);
    }
    return T;
}

// The two priors. q[s] = P(allele index 1 at site s) for the common-content control (renaming gate uses q -> 1 - q).
long double log_prior_original(const std::vector<std::array<int, 2>>& al, const Scaffold& g, double p) {
    long double lp = 0;
    for (std::size_t c = 0; c < al.size(); ++c)
        for (int s = 0; s < 2; ++s) lp += al[c][s] != g.base[c][s] ? std::log(static_cast<long double>(p)) : std::log1p(static_cast<long double>(-p));
    return lp;
}
long double log_prior_common(const std::vector<std::array<int, 2>>& al, const std::array<double, 2>& q) {
    long double lp = 0;
    for (const auto& c : al)
        for (int s = 0; s < 2; ++s) lp += c[s] ? std::log(static_cast<long double>(q[s])) : std::log1p(static_cast<long double>(-q[s]));
    return lp;
}

struct Marginal {
    long double log_marginal = 0, prior_mass = 0, posterior_mean_log_likelihood = 0, kl = 0, identity_residual = 0;
    double top_state_posterior = 0, effective_states = 0;
};

Marginal marginalise(const std::vector<double>& logl, const std::vector<long double>& lp) {
    Marginal m;
    std::vector<long double> terms(logl.size());
    std::size_t best = 0;
    for (std::size_t t = 0; t < logl.size(); ++t) {
        m.prior_mass += std::exp(lp[t]);
        terms[t] = lp[t] + logl[t];
        if (terms[t] > terms[best]) best = t;
    }
    m.log_marginal = logsumexp(terms);
    long double h = 0;
    for (std::size_t t = 0; t < logl.size(); ++t) {
        if (!std::isfinite(static_cast<double>(terms[t]))) continue;
        const long double q = std::exp(terms[t] - m.log_marginal);
        m.posterior_mean_log_likelihood += q * logl[t];
        m.kl += q * (terms[t] - m.log_marginal - lp[t]);   // q log(q / prior), with log q = lp + logL - log m
        h += q * q;
    }
    m.identity_residual = m.log_marginal - (m.posterior_mean_log_likelihood - m.kl);
    m.top_state_posterior = static_cast<double>(std::exp(terms[best] - m.log_marginal));
    m.effective_states = static_cast<double>(1.0L / h);
    return m;
}

std::vector<long double> priors_original(const Table& T, double p) {
    std::vector<long double> v;
    for (const auto& al : T.alleles) v.push_back(log_prior_original(al, T.g, p));
    return v;
}
std::vector<long double> priors_common(const Table& T, const std::array<double, 2>& q) {
    std::vector<long double> v;
    for (const auto& al : T.alleles) v.push_back(log_prior_common(al, q));
    return v;
}

// Independent exhaustive oracle: each absolute state scored ALONE (own interning, own call, one thread), its prior
// computed in-loop from the absolute alleles, summed; no deduplication, no table.
std::array<long double, 2> oracle(const World& w, const Scaffold& g, const std::vector<panvar::ReadPair>& reads, double p, double q) {
    const int n = g.nA + g.nB, bits = 2 * n;
    std::vector<long double> orig, comm;
    for (std::size_t t = 0; t < (std::size_t(1) << bits); ++t) {
        std::vector<std::array<int, 2>> al(n);
        long double lo = 0, lc = 0;
        for (int c = 0; c < n; ++c)
            for (int s = 0; s < 2; ++s) {
                al[c][s] = static_cast<int>((t >> (2 * c + s)) & 1u);   // absolute enumeration (independent of the base)
                lo += al[c][s] == g.base[c][s] ? std::log1p(static_cast<long double>(-p)) : std::log(static_cast<long double>(p));
                lc += al[c][s] == 1 ? std::log(static_cast<long double>(q)) : std::log1p(static_cast<long double>(-q));
            }
        const std::string a = spell_hap(w, g, al, true), b = spell_hap(w, g, al, false);
        const panvar::InternedHaplotypes in = panvar::intern_haplotypes({panvar::MosaicPath{{0}}, panvar::MosaicPath{{1}}}, {a, b});
        const double l = panvar::score_diplotypes(in, {{in.index_of[0], in.index_of[1]}}, reads, library(), scoring(1)).scores[0].log_likelihood;
        orig.push_back(lo + l);
        comm.push_back(lc + l);
    }
    return {logsumexp(orig), logsumexp(comm)};
}

std::size_t g_failed = 0;
std::ostringstream g_rows;
void row(const std::string& t, const std::string& item, double v, const std::string& expect, bool pass) {
    g_rows << t << '\t' << item << '\t' << std::setprecision(15) << v << '\t' << expect << '\t' << (pass ? "PASS" : "FAIL") << '\n';
    g_failed += pass ? 0 : 1;
}
const std::array<double, 2> kQ = {kAltProbability, kAltProbability};

int gates() {
    const World w = make_world(5, 200);   // M6-C's gate world, truth and reads
    Scaffold truth = scaffold(1, 1, 0);
    truth.base[1] = {1, 0};
    const auto reads = reads_of(w, truth, 99);

    // HA: both priors normalise within every scaffold (all nine allocations, both base annotations)
    {
        double worst = 0;
        for (int base = 0; base <= 1; ++base)
            for (const Scaffold& g : catalogue(base, base)) {
                const int bits = 2 * (g.nA + g.nB);
                long double mo = 0, mc = 0;
                for (std::size_t t = 0; t < (std::size_t(1) << bits); ++t) {
                    std::vector<std::array<int, 2>> al = g.base;
                    for (int b = 0; b < bits; ++b) if ((t >> b) & 1u) al[b / 2][b % 2] ^= 1;
                    mo += std::exp(log_prior_original(al, g, kNuisanceFlip));
                    mc += std::exp(log_prior_common(al, kQ));
                }
                worst = std::max({worst, std::fabs(static_cast<double>(mo) - 1.0), std::fabs(static_cast<double>(mc) - 1.0)});
            }
        row("HA_normalisation", "max |sum P(theta|G) - 1| over 9 allocations x 2 base annotations x 2 priors", worst, "<= 1e-15", worst <= 1e-15);
    }
    // HB: under the common-content control, changing ONLY the base annotation leaves the marginal unchanged; the
    // original prior changes (contrast). Separate tables, different state order.
    {
        const Table r0 = build(w, scaffold(2, 1, 0), reads, 4), r1 = build(w, scaffold(2, 1, 1), reads, 4);
        const auto c0 = marginalise(r0.log_likelihood, priors_common(r0, kQ)), c1 = marginalise(r1.log_likelihood, priors_common(r1, kQ));
        const auto o0 = marginalise(r0.log_likelihood, priors_original(r0, kNuisanceFlip)), o1 = marginalise(r1.log_likelihood, priors_original(r1, kNuisanceFlip));
        row("HB_base_annotation", "common-content: marginal(base ref) - marginal(base alt), 2+1", static_cast<double>(c0.log_marginal - c1.log_marginal), "|.| <= 1e-9",
            std::fabs(static_cast<double>(c0.log_marginal - c1.log_marginal)) <= 1e-9);
        row("HB_base_annotation", "original (contrast): marginal(base ref) - marginal(base alt), 2+1", static_cast<double>(o0.log_marginal - o1.log_marginal),
            "!= 0 (the original prior depends on the annotation)", std::fabs(static_cast<double>(o0.log_marginal - o1.log_marginal)) > 1e-6);
    }
    // HC: renaming ref/alt with consistently transformed probabilities leaves results unchanged
    {
        World rn = w;
        std::swap(rn.ref, rn.alt);   // allele index 0 now spells the old alt
        Scaffold g = scaffold(2, 1, 0), gr = scaffold(2, 1, 1);   // complemented base: same spelled base sequences
        const Table a = build(w, g, reads, 4), b = build(rn, gr, reads, 4);
        const auto ca = marginalise(a.log_likelihood, priors_common(a, kQ));
        const auto cb = marginalise(b.log_likelihood, priors_common(b, {1.0 - kAltProbability, 1.0 - kAltProbability}));
        const auto oa = marginalise(a.log_likelihood, priors_original(a, kNuisanceFlip)), ob = marginalise(b.log_likelihood, priors_original(b, kNuisanceFlip));
        row("HC_renaming", "common-content: original names - renamed (q -> 1-q)", static_cast<double>(ca.log_marginal - cb.log_marginal), "|.| <= 1e-9",
            std::fabs(static_cast<double>(ca.log_marginal - cb.log_marginal)) <= 1e-9);
        row("HC_renaming", "original: original names - renamed (base complemented)", static_cast<double>(oa.log_marginal - ob.log_marginal), "|.| <= 1e-9",
            std::fabs(static_cast<double>(oa.log_marginal - ob.log_marginal)) <= 1e-9);
    }
    // HD: splitting or merging synonymous state labels preserves total prior mass and the marginal
    {
        const World ws = make_world(5, 200, true);
        const auto rs = reads_of(ws, truth, 98);
        const Table T = build(ws, scaffold(1, 1, 0), rs, 4);
        for (int which = 0; which < 2; ++which) {
            const std::vector<long double> lp = which == 0 ? priors_original(T, kNuisanceFlip) : priors_common(T, kQ);
            const std::string name = which == 0 ? "original" : "common-content";
            const Marginal lab = marginalise(T.log_likelihood, lp);
            // merge: one label per unique diplotype, prior mass summed
            std::vector<long double> mass(T.unique_diplotypes, 0);
            std::vector<double> ll(T.unique_diplotypes, 0);
            for (std::size_t t = 0; t < T.states; ++t) { mass[T.slot[t]] += std::exp(lp[t]); ll[T.slot[t]] = T.log_likelihood[t]; }
            std::vector<long double> lm;
            long double tot = 0;
            for (auto x : mass) { lm.push_back(std::log(x)); tot += x; }
            const Marginal mer = marginalise(ll, lm);
            // split: every label into two labels of half its mass
            std::vector<double> l2;
            std::vector<long double> p2;
            for (std::size_t t = 0; t < T.states; ++t)
                for (int k = 0; k < 2; ++k) { l2.push_back(T.log_likelihood[t]); p2.push_back(lp[t] - std::log(2.0L)); }
            const Marginal spl = marginalise(l2, p2);
            row("HD_synonymous_labels", name + ": " + std::to_string(T.states) + " labels -> " + std::to_string(T.unique_diplotypes) + " merged: marginal difference",
                static_cast<double>(lab.log_marginal - mer.log_marginal), "|.| <= 1e-12", std::fabs(static_cast<double>(lab.log_marginal - mer.log_marginal)) <= 1e-12);
            row("HD_synonymous_labels", name + ": merged total prior mass", static_cast<double>(tot), "1 (|.-1| <= 1e-15)", std::fabs(static_cast<double>(tot) - 1.0) <= 1e-15);
            row("HD_synonymous_labels", name + ": split (2x labels) marginal difference", static_cast<double>(lab.log_marginal - spl.log_marginal), "|.| <= 1e-12",
                std::fabs(static_cast<double>(lab.log_marginal - spl.log_marginal)) <= 1e-12);
            row("HD_synonymous_labels", name + ": KL(labels) - KL(merged)", static_cast<double>(lab.kl - mer.kl), "|.| <= 1e-12",
                std::fabs(static_cast<double>(lab.kl - mer.kl)) <= 1e-12);
        }
        row("HD_synonymous_labels", "unique diplotypes (synonymous site collapses labels)", static_cast<double>(T.unique_diplotypes), "4", T.unique_diplotypes == 4);
    }
    // HE: a likelihood constant across nuisance states is returned exactly, with no state-count penalty (1,024 states)
    {
        Table T;
        T.g = scaffold(2, 3, 0);
        T.states = 1024;
        for (std::size_t t = 0; t < T.states; ++t) {
            std::vector<std::array<int, 2>> al = T.g.base;
            for (int b = 0; b < 10; ++b) if ((t >> b) & 1u) al[b / 2][b % 2] ^= 1;
            T.alleles.push_back(al);
        }
        const std::vector<double> flat(T.states, -1234.5);
        const auto mo = marginalise(flat, priors_original(T, kNuisanceFlip)), mc = marginalise(flat, priors_common(T, kQ));
        row("HE_constant_likelihood", "original: marginal - constant (1,024 states)", static_cast<double>(mo.log_marginal + 1234.5L), "|.| <= 1e-12 (no -log 1024)",
            std::fabs(static_cast<double>(mo.log_marginal + 1234.5L)) <= 1e-12);
        row("HE_constant_likelihood", "common-content: marginal - constant (1,024 states)", static_cast<double>(mc.log_marginal + 1234.5L), "|.| <= 1e-12 (no -log 1024)",
            std::fabs(static_cast<double>(mc.log_marginal + 1234.5L)) <= 1e-12);
    }
    // HF: table marginals match an independent exhaustive calculation, both priors; HG: the KL identity holds exactly
    for (const Scaffold& g : {scaffold(1, 1, 0), scaffold(2, 1, 1)}) {
        const Table T = build(w, g, reads, 4);
        const auto mo = marginalise(T.log_likelihood, priors_original(T, kNuisanceFlip)), mc = marginalise(T.log_likelihood, priors_common(T, kQ));
        const auto o = oracle(w, g, reads, kNuisanceFlip, kAltProbability);
        const std::string tag = g.label() + (g.base[0][0] ? " base alt" : " base ref");
        row("HF_oracle", tag + ": original table - oracle", static_cast<double>(mo.log_marginal - o[0]), "|.| <= 1e-6", std::fabs(static_cast<double>(mo.log_marginal - o[0])) <= 1e-6);
        row("HF_oracle", tag + ": common-content table - oracle", static_cast<double>(mc.log_marginal - o[1]), "|.| <= 1e-6", std::fabs(static_cast<double>(mc.log_marginal - o[1])) <= 1e-6);
        row("HG_kl_identity", tag + ": original: log m - (E_q log L - KL)", static_cast<double>(mo.identity_residual), "|.| <= 1e-9", std::fabs(static_cast<double>(mo.identity_residual)) <= 1e-9);
        row("HG_kl_identity", tag + ": common-content: log m - (E_q log L - KL)", static_cast<double>(mc.identity_residual), "|.| <= 1e-9", std::fabs(static_cast<double>(mc.identity_residual)) <= 1e-9);
    }
    std::cout << "test\titem\tvalue\texpectation\tresult\n" << g_rows.str() << "prior_audit_gates\t" << (g_failed ? "FAIL" : "PASS") << '\n';
    return g_failed ? 1 : 0;
}

// The four existing fixtures (M6-C DESIGN.md), copied verbatim; nothing else is run.
struct Fixture {
    std::string name;
    std::size_t unit_bp;
    Scaffold truth;
    int base_balanced, base_unbalanced;
};
std::vector<Fixture> fixtures() {
    return {{"F2a_missing_variants_invisible_allocation", 600, scaffold(2, 2, 1), 0, 1},
            {"F2b_missing_variants_informative_allocation", 200, scaffold(2, 2, 1), 0, 1},
            {"F4_adjacency_true_3+1", 200, scaffold(3, 1, 0), 0, 0},
            {"F4_adjacency_true_2+2", 200, scaffold(2, 2, 0), 0, 0}};
}

double peak_rss_gib() {
    struct rusage u {};
    getrusage(RUSAGE_SELF, &u);
    return static_cast<double>(u.ru_maxrss) / (1024.0 * 1024.0 * 1024.0);   // bytes on macOS
}

std::string alleles_string(const std::vector<std::array<int, 2>>& al, int nA) {
    std::string s;
    for (std::size_t c = 0; c < al.size(); ++c) {
        if (static_cast<int>(c) == nA) s += '|';
        s += static_cast<char>('0' + al[c][0]);
        s += static_cast<char>('0' + al[c][1]);
    }
    return s;
}

int run(const std::string& name, std::uint64_t seed, const std::string& out) {
    const auto t0 = std::chrono::steady_clock::now();
    Fixture fx;
    bool found = false;
    for (const auto& f : fixtures()) if (f.name == name) { fx = f; found = true; }
    if (!found) throw std::runtime_error("not one of the four M6-D fixtures: " + name);
    const World w = make_world(seed, fx.unit_bp);
    const auto reads = reads_of(w, fx.truth, seed * 1000 + 1);
    std::ostringstream st, sc;
    st << std::setprecision(17) << "scaffold\tstate\talleles_absolute\talt_count\tflips_from_base\tdiplotype_slot\tlog_likelihood\tfragment_term\texposure_term\t"
                                   "log_prior_original\tlog_prior_common_content\n";
    sc << std::setprecision(15) << "fixture\tseed\tscaffold\tbase_allele\tnuisance_states\tunique_haplotypes\tunique_diplotypes\tprior\tlog_marginal\tprior_mass\t"
                                   "posterior_mean_log_likelihood\tkl_posterior_prior\tidentity_residual\ttop_state_posterior\teffective_states\n";
    for (const Scaffold& g : catalogue(fx.base_balanced, fx.base_unbalanced)) {
        const Table T = build(w, g, reads, 4);
        const auto po = priors_original(T, kNuisanceFlip), pc = priors_common(T, kQ);
        for (std::size_t t = 0; t < T.states; ++t) {
            int alt = 0, flips = 0;
            for (std::size_t c = 0; c < T.alleles[t].size(); ++c)
                for (int s = 0; s < 2; ++s) { alt += T.alleles[t][c][s]; flips += T.alleles[t][c][s] != g.base[c][s]; }
            st << g.label() << '\t' << t << '\t' << alleles_string(T.alleles[t], g.nA) << '\t' << alt << '\t' << flips << '\t' << T.slot[t] << '\t'
               << T.log_likelihood[t] << '\t' << T.fragment[t] << '\t' << T.exposure[t] << '\t' << static_cast<double>(po[t]) << '\t' << static_cast<double>(pc[t]) << '\n';
        }
        for (int which = 0; which < 2; ++which) {
            const Marginal m = marginalise(T.log_likelihood, which == 0 ? po : pc);
            sc << name << '\t' << seed << '\t' << g.label() << '\t' << g.base[0][0] << '\t' << T.states << '\t' << T.unique_haplotypes << '\t' << T.unique_diplotypes << '\t'
               << (which == 0 ? "original" : "common_content") << '\t' << static_cast<double>(m.log_marginal) << '\t' << static_cast<double>(m.prior_mass) << '\t'
               << static_cast<double>(m.posterior_mean_log_likelihood) << '\t' << static_cast<double>(m.kl) << '\t' << static_cast<double>(m.identity_residual) << '\t'
               << m.top_state_posterior << '\t' << m.effective_states << '\n';
        }
    }
    const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::ofstream(out + ".states.tsv") << st.str();
    std::ofstream(out + ".scaffolds.tsv") << sc.str();
    std::ofstream(out + ".run.tsv") << "key\tvalue\nfixture\t" << name << "\nseed\t" << seed << "\ntruth\t" << fx.truth.label() << "\ntruth_copy_alleles\t"
                                    << fx.truth.base[0][0] << "\nunit_bp\t" << fx.unit_bp << "\nfragments\t" << reads.size() << "\noriginal_flip_prior\t" << kNuisanceFlip
                                    << "\ncommon_content_alt_probability\t" << kAltProbability << "\nseconds\t" << secs << "\npeak_rss_gib\t" << peak_rss_gib() << '\n';
    std::cout << name << " seed " << seed << ": " << reads.size() << " fragments, " << secs << " s\n";
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const std::string mode = argc > 1 ? argv[1] : "";
        if (mode == "gates") return gates();
        if (mode == "run" && argc == 5) return run(argv[2], std::stoull(argv[3]), argv[4]);
        std::cerr << "usage: nuisance_prior_audit gates | run <fixture> <seed> <out prefix>\n";
        return 2;
    } catch (const std::exception& e) {
        std::cerr << "nuisance_prior_audit: ERROR " << e.what() << '\n';
        return 2;
    }
}
