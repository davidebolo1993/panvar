// M6-G: a tiny EXACT prototype bridging the C1 sequence prior to read-based structural inference (diagnostic only; not a
// caller, not reconstruction, never a real-sample run).
//
// Read model, layer by layer:
//   G      structural arrangement: homologue A = F_A + nA copies + G_A, homologue B = F_B + nB copies + G_B (synthetic
//          unique flanks); catalogue 1+1, 1+2, 2+1, 1+3, 2+2, 3+1 with uniform structural weights (explicit; the same
//          for every compared method)
//   X_j    each copy's biological sequence, ONE state shared by the whole read set (never redrawn per read):
//          X_j in S, the admitted set; P(X_j) = C1 conditioned on S, renormalised:
//            P_C1(X) = sum_k w_k prod_x e_res(X_x | b_k(x)),  e_res = 1 - eps_res on match, eps_res / 3 otherwise
//          method RES:   S = {library window k} x {ancestral / declared novel allele at 2 admitted sites}  (8 per copy)
//          method PANEL: S = {library window k}                                                       (2 per copy)
//          copies independent given G (no coupling); the joint prior is the product
//   reads  the existing paired-fragment likelihood of the complete diplotype (score_diplotypes: placements summed,
//          sequencing substitution error 0.002, exposure, outlier/background mixture)
// Marginal: log P(reads | G) = log sum_X P(X | G) P(reads | G, X), exhaustively over every distinct admitted state.
// Negative control (per-read independent marginalisation, never used for inference):
//   NC(G) = sum_f log E_X[exp(term_f(X))] + log E_X[exp(-expected signal(X))]
//
//   c1_structural gates <library.tsv>
//   c1_structural run <library.tsv> <fixture> <out prefix>
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

constexpr int kReadLength = 100;

panvar::LibraryModel make_library() {
    panvar::LibraryModel lib;
    lib.insert_mean_bp = 350;
    lib.insert_sd_bp = 35;
    lib.support_sigmas = 4;
    lib.substitution_error = 0.002;   // sequencing error: distinct from the residual divergence rate eps_res
    lib.outlier_weight = 0.001;
    lib.fragment_start_intensity = 20.0 / (4.0 * kReadLength);
    lib.provenance = "M6-G synthetic declared";
    return lib;
}
const panvar::LibraryModel& library() {
    static const panvar::LibraryModel lib = make_library();
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
long double lse(const std::vector<long double>& v) {
    long double m = -INFINITY;
    for (auto x : v) m = std::max(m, x);
    if (!std::isfinite(static_cast<double>(m))) return m;
    long double s = 0;
    for (auto x : v) s += std::exp(x - m);
    return m + std::log(s);
}
struct Sum {
    long double s = 0, c = 0;
    void add(long double x) {
        const long double t = s + x;
        c += std::fabs(s) >= std::fabs(x) ? (s - t) + x : (x - t) + s;
        s = t;
    }
    long double value() const { return s + c; }
};

// ---- the frozen library -------------------------------------------------------------------------------------------
struct Library {
    double eps_res = 0;
    std::vector<std::string> names;
    std::vector<double> weights;
    std::map<int, std::vector<std::string>> window;   // L -> window sequence per library entry
    std::map<int, std::array<int, 2>> admitted;
    std::map<int, std::array<char, 2>> novel;
};
Library read_library(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot read " + path);
    std::map<std::string, std::string> kv;
    for (std::string line; std::getline(in, line);) {
        const auto t = line.find('\t');
        if (t != std::string::npos) kv[line.substr(0, t)] = line.substr(t + 1);
    }
    const auto split = [](const std::string& s) {
        std::vector<std::string> v;
        std::stringstream ss(s);
        for (std::string x; std::getline(ss, x, ',');) v.push_back(x);
        return v;
    };
    Library L;
    L.eps_res = std::stod(kv.at("eps_res"));
    L.names = split(kv.at("components"));
    for (const auto& w : split(kv.at("weights"))) L.weights.push_back(std::stod(w));
    for (int len : {200, 600}) {
        for (const auto& n : L.names) L.window[len].push_back(kv.at("window_" + std::to_string(len) + "_seq_" + n));
        const auto a = split(kv.at("window_" + std::to_string(len) + "_admitted_sites"));
        const auto nv = split(kv.at("window_" + std::to_string(len) + "_novel_alleles"));
        L.admitted[len] = {std::stoi(a[0]), std::stoi(a[1])};
        L.novel[len] = {nv[0][0], nv[1][0]};
    }
    return L;
}

// The admitted per-copy state space and its exact conditioned C1 prior. Entries may alias (the same window twice with
// split weight): states are DISTINCT sequences, and P(X) sums over every entry.
struct CopySpace {
    std::vector<std::string> seq;        // distinct admitted sequences
    std::vector<int> backbone, pattern;  // backbone = index of the distinct window; pattern bit s = novel allele at site s
    std::vector<long double> log_prior;  // normalised over the admitted set
    long double prior_mass = 0;
};
CopySpace copy_space(const std::vector<std::string>& windows, const std::vector<double>& weights, const std::array<int, 2>& sites,
                     const std::array<char, 2>& novel, double eps, bool residual) {
    std::vector<std::string> distinct;
    for (const auto& w : windows) {
        bool seen = false;
        for (const auto& d : distinct) seen |= d == w;
        if (!seen) distinct.push_back(w);
    }
    CopySpace cs;
    for (std::size_t b = 0; b < distinct.size(); ++b)
        for (int p = 0; p < (residual ? 4 : 1); ++p) {
            std::string x = distinct[b];
            for (int s = 0; s < 2; ++s) if ((p >> s) & 1) x[sites[s]] = novel[s];
            cs.seq.push_back(x);
            cs.backbone.push_back(static_cast<int>(b));
            cs.pattern.push_back(p);
        }
    std::vector<long double> raw;
    for (const auto& x : cs.seq) {
        std::vector<long double> t;
        for (std::size_t e = 0; e < windows.size(); ++e) {
            std::size_t mis = 0;
            for (std::size_t i = 0; i < x.size(); ++i) mis += x[i] != windows[e][i];
            t.push_back(std::log(static_cast<long double>(weights[e])) + (x.size() - mis) * std::log1p(-static_cast<long double>(eps)) +
                        mis * std::log(static_cast<long double>(eps) / 3.0L));
        }
        raw.push_back(lse(t));
    }
    const long double z = lse(raw);
    Sum m;
    for (auto r : raw) { cs.log_prior.push_back(r - z); m.add(std::exp(r - z)); }
    cs.prior_mass = m.value();
    return cs;
}

struct Structure {
    int nA, nB;
    std::string label() const { return std::to_string(nA) + "+" + std::to_string(nB); }
};
const std::vector<Structure> kCatalogue = {{1, 1}, {1, 2}, {2, 1}, {1, 3}, {2, 2}, {3, 1}};

struct World {
    std::string FA, GA, FB, GB;
};
World make_world(std::uint64_t seed) {
    std::mt19937_64 rng(seed * 7919 + 13);
    return {rnd(rng, 1200), rnd(rng, 1200), rnd(rng, 1200), rnd(rng, 1200)};
}
std::string spell(const World& w, const Structure& g, const std::vector<std::string>& copy_seq, bool hapA) {
    std::string s = hapA ? w.FA : w.FB;
    const int from = hapA ? 0 : g.nA, to = hapA ? g.nA : g.nA + g.nB;
    for (int c = from; c < to; ++c) s += copy_seq[c];
    return s + (hapA ? w.GA : w.GB);
}

// Every distinct admitted state of every structure, scored in ONE engine call; per-fragment terms reconstructed.
struct Scored {
    std::vector<std::vector<std::vector<int>>> states;   // [structure][state] -> per-copy index into the copy space
    std::vector<std::vector<double>> logl;               // [structure][state]
    std::vector<std::vector<std::size_t>> cand;          // [structure][state] -> candidate index
    std::vector<std::pair<std::size_t, std::size_t>> cand_haps;
    panvar::FragmentScoringResult result;
    std::vector<double> background;
    std::size_t unique_haplotypes = 0;
    double max_reconstruction_error = 0;
};
std::vector<std::vector<int>> enumerate(const Structure& g, std::size_t space, const std::vector<int>& base_offset) {
    const int c = g.nA + g.nB;
    std::size_t total = 1;
    for (int i = 0; i < c; ++i) total *= space;
    if (total > 4096) throw std::runtime_error("more than 4,096 distinct admitted states for " + g.label());
    std::vector<std::vector<int>> out;
    for (std::size_t t = 0; t < total; ++t) {
        std::vector<int> d(c);
        std::size_t x = t;
        for (int i = 0; i < c; ++i) { d[i] = static_cast<int>((x % space + base_offset[i]) % space); x /= space; }
        out.push_back(d);
    }
    return out;
}
Scored score_all(const World& w, const CopySpace& cs, const std::vector<Structure>& cat, const std::vector<panvar::ReadPair>& reads,
                 const std::vector<std::vector<int>>& base_offsets, std::size_t threads) {
    Scored S;
    std::vector<std::string> seqs;
    std::vector<panvar::MosaicPath> paths;
    std::vector<std::pair<std::size_t, std::size_t>> raw;
    for (std::size_t gi = 0; gi < cat.size(); ++gi) {
        S.states.push_back(enumerate(cat[gi], cs.seq.size(), base_offsets[gi]));
        S.cand.emplace_back();
        for (const auto& st : S.states.back()) {
            std::vector<std::string> cp;
            for (int d : st) cp.push_back(cs.seq[d]);
            seqs.push_back(spell(w, cat[gi], cp, true));
            seqs.push_back(spell(w, cat[gi], cp, false));
            paths.push_back(panvar::MosaicPath{{seqs.size() - 2}});
            paths.push_back(panvar::MosaicPath{{seqs.size() - 1}});
            raw.push_back({seqs.size() - 2, seqs.size() - 1});
        }
    }
    const panvar::InternedHaplotypes in = panvar::intern_haplotypes(paths, seqs);
    std::map<std::pair<std::size_t, std::size_t>, std::size_t> slot;
    std::vector<panvar::HaplotypeCandidate> cands;
    std::size_t r = 0;
    for (std::size_t gi = 0; gi < cat.size(); ++gi)
        for (std::size_t t = 0; t < S.states[gi].size(); ++t, ++r) {
            const std::pair<std::size_t, std::size_t> key{in.index_of[raw[r].first], in.index_of[raw[r].second]};   // ordered: A, B
            auto it = slot.find(key);
            if (it == slot.end()) {
                it = slot.emplace(key, cands.size()).first;
                cands.push_back({key.first, key.second});
                S.cand_haps.push_back(key);
            }
            S.cand[gi].push_back(it->second);
        }
    panvar::FragmentScoringOptions opt = scoring(threads);
    S.result = panvar::score_diplotypes(in, cands, reads, library(), opt);
    S.unique_haplotypes = in.sequences.size();
    for (const auto& f : reads) S.background.push_back(panvar::log_background_probability(f));
    for (std::size_t gi = 0; gi < cat.size(); ++gi) {
        S.logl.emplace_back();
        for (std::size_t t = 0; t < S.states[gi].size(); ++t) S.logl[gi].push_back(S.result.scores[S.cand[gi][t]].log_likelihood);
    }
    return S;
}
// Per-fragment terms of one candidate and its exposure term, from the engine's own masses (gate GF certifies the sum).
void fragment_terms(const Scored& S, std::size_t c, const std::vector<panvar::ReadPair>& reads, std::vector<double>& terms, double& exposure_term) {
    const std::size_t U = S.result.unique_haplotypes, a = S.cand_haps[c].first, b = S.cand_haps[c].second;
    terms.resize(reads.size());
    for (std::size_t f = 0; f < reads.size(); ++f) {
        const double ma = S.result.log_mass[f * U + a], mb = S.result.log_mass[f * U + b];
        const double comb = std::isinf(ma) && ma < 0 ? mb : std::isinf(mb) && mb < 0 ? ma : std::max(ma, mb) + std::log1p(std::exp(-std::fabs(ma - mb)));
        terms[f] = panvar::log_fragment_contribution(comb, S.background[f], library());
    }
    exposure_term = -panvar::expected_signal_fragments(S.result.exposure[a] + S.result.exposure[b], library());
}

struct Marginal {
    long double log_m = 0, elogl = 0, kl = 0, residual = 0, mass = 0, novel_mean = 0, any_novel = 0;
    double top = 0, eff = 0;
};
Marginal marginalise(const std::vector<double>& logl, const std::vector<long double>& lp, const std::vector<int>& novel_count) {
    Marginal m;
    std::vector<long double> t(logl.size());
    std::size_t best = 0;
    Sum mass;
    for (std::size_t i = 0; i < logl.size(); ++i) {
        mass.add(std::exp(lp[i]));
        t[i] = lp[i] + logl[i];
        if (t[i] > t[best]) best = i;
    }
    m.mass = mass.value();
    const long double top = t[best];
    std::vector<long double> w(t.size());
    Sum ws;
    for (std::size_t i = 0; i < t.size(); ++i) { w[i] = std::isfinite(static_cast<double>(t[i])) ? std::exp(t[i] - top) : 0; ws.add(w[i]); }
    const long double W = ws.value(), lW = std::log(W);
    m.log_m = top + lW;
    Sum el, kl, h, nm, an;
    for (std::size_t i = 0; i < t.size(); ++i) {
        if (w[i] == 0) continue;
        const long double q = w[i] / W;
        el.add(q * logl[i]);
        kl.add(q * (((t[i] - top) - lW) - lp[i]));
        h.add(q * q);
        nm.add(q * novel_count[i]);
        an.add(novel_count[i] > 0 ? q : 0);
    }
    m.elogl = el.value();
    m.kl = kl.value();
    m.residual = m.log_m - (m.elogl - m.kl);
    m.top = static_cast<double>(w[best] / W);
    m.eff = static_cast<double>(1.0L / h.value());
    m.novel_mean = nm.value();
    m.any_novel = an.value();
    return m;
}
struct Method {
    std::string name;
    bool residual;
};
// prior of every state of a structure under a method (states outside the method's admitted set get -inf)
std::vector<long double> state_prior(const CopySpace& full, const CopySpace& sub, const std::vector<std::vector<int>>& states, bool residual) {
    std::vector<long double> lp;
    for (const auto& st : states) {
        long double s = 0;
        for (int d : st) {
            if (!residual) {
                if (full.pattern[d] != 0) { s = -INFINITY; break; }
                s += sub.log_prior[full.backbone[d]];
            } else {
                s += full.log_prior[d];
            }
        }
        lp.push_back(s);
    }
    return lp;
}
std::vector<int> novel_counts(const CopySpace& cs, const std::vector<std::vector<int>>& states) {
    std::vector<int> v;
    for (const auto& st : states) {
        int n = 0;
        for (int d : st) n += __builtin_popcount(static_cast<unsigned>(cs.pattern[d]));
        v.push_back(n);
    }
    return v;
}
long double negative_control(const Scored& S, std::size_t gi, const std::vector<long double>& lp, const std::vector<panvar::ReadPair>& reads) {
    const std::size_t n = S.states[gi].size(), F = reads.size();
    std::vector<std::vector<long double>> per_f(F);
    std::vector<long double> expo;
    std::vector<double> terms;
    double e = 0;
    for (std::size_t t = 0; t < n; ++t) {
        if (!std::isfinite(static_cast<double>(lp[t]))) continue;
        fragment_terms(S, S.cand[gi][t], reads, terms, e);
        for (std::size_t f = 0; f < F; ++f) per_f[f].push_back(lp[t] + terms[f]);
        expo.push_back(lp[t] + e);
    }
    Sum s;
    for (std::size_t f = 0; f < F; ++f) s.add(lse(per_f[f]));
    s.add(lse(expo));
    return s.value();
}

double peak_rss_gib() {
    struct rusage u {};
    getrusage(RUSAGE_SELF, &u);
    return static_cast<double>(u.ru_maxrss) / (1024.0 * 1024.0 * 1024.0);
}

// ---- fixtures (frozen) -----------------------------------------------------------------------------------------------
struct Fixture {
    std::string name;
    int L;
    Structure truth;
    std::vector<std::pair<int, int>> copies;   // (library entry, novel pattern) per copy, A copies first
    std::uint64_t seed;
};
std::vector<Fixture> fixtures() {
    return {{"A_missing_variation", 600, {1, 2}, {{0, 3}, {0, 3}, {1, 3}}, 101},
            {"B_unbalanced_allocation", 200, {3, 1}, {{0, 0}, {1, 0}, {0, 0}, {1, 0}}, 102},
            {"C_uninformative_allocation", 600, {2, 2}, {{0, 0}, {0, 0}, {0, 0}, {0, 0}}, 103}};
}
std::vector<panvar::ReadPair> fixture_reads(const World& w, const Library& lib, const Fixture& fx) {
    std::vector<std::string> cp;
    for (auto [e, p] : fx.copies) {
        std::string x = lib.window.at(fx.L)[e];
        for (int s = 0; s < 2; ++s) if ((p >> s) & 1) x[lib.admitted.at(fx.L)[s]] = lib.novel.at(fx.L)[s];
        cp.push_back(x);
    }
    std::mt19937_64 rng(fx.seed * 1000 + 1);
    std::vector<panvar::ReadPair> r;
    simulate(spell(w, fx.truth, cp, true), library().fragment_start_intensity, rng, r, "a");
    simulate(spell(w, fx.truth, cp, false), library().fragment_start_intensity, rng, r, "b");
    return r;
}
std::vector<std::vector<int>> zero_offsets(const std::vector<Structure>& cat) {
    std::vector<std::vector<int>> o;
    for (const auto& g : cat) o.push_back(std::vector<int>(g.nA + g.nB, 0));
    return o;
}

std::size_t g_failed = 0;
std::ostringstream g_rows;
void row(const std::string& t, const std::string& item, double v, const std::string& expect, bool pass) {
    g_rows << t << '\t' << item << '\t' << std::setprecision(15) << v << '\t' << expect << '\t' << (pass ? "PASS" : "FAIL") << '\n';
    g_failed += pass ? 0 : 1;
}

int gates(const std::string& libpath) {
    const Library lib = read_library(libpath);
    const Fixture fx = fixtures()[1];   // the 200-bp fixture's world and reads (cheap oracle)
    const World w = make_world(fx.seed);
    const auto reads = fixture_reads(w, lib, fx);
    const CopySpace res = copy_space(lib.window.at(200), lib.weights, lib.admitted.at(200), lib.novel.at(200), lib.eps_res, true);
    const CopySpace pan = copy_space(lib.window.at(200), lib.weights, lib.admitted.at(200), lib.novel.at(200), lib.eps_res, false);
    // GA normalisation (per copy and joint)
    row("GA_normalisation", "per-copy prior mass, RES (8 states)", static_cast<double>(res.prior_mass), "1 (|.-1| <= 1e-12)", std::fabs(static_cast<double>(res.prior_mass) - 1) <= 1e-12);
    row("GA_normalisation", "per-copy prior mass, PANEL (2 states)", static_cast<double>(pan.prior_mass), "1 (|.-1| <= 1e-12)", std::fabs(static_cast<double>(pan.prior_mass) - 1) <= 1e-12);
    const std::vector<Structure> small = {{1, 1}, {1, 2}};
    const Scored S = score_all(w, res, small, reads, zero_offsets(small), 4);
    double jm = 0;
    for (std::size_t gi = 0; gi < small.size(); ++gi)
        for (bool r : {true, false}) {
            const auto m = marginalise(S.logl[gi], state_prior(res, pan, S.states[gi], r), novel_counts(res, S.states[gi]));
            jm = std::max(jm, std::fabs(static_cast<double>(m.mass) - 1));
        }
    row("GA_normalisation", "joint prior mass over every enumerated state (1+1, 1+2; RES and PANEL)", jm, "<= 1e-12", jm <= 1e-12);
    // GF per-fragment reconstruction equals the engine's candidate totals
    double worst = 0;
    std::vector<double> terms;
    double e = 0;
    for (std::size_t c = 0; c < S.cand_haps.size(); ++c) {
        fragment_terms(S, c, reads, terms, e);
        Sum s;
        for (double t : terms) s.add(t);
        s.add(e);
        worst = std::max(worst, std::fabs(static_cast<double>(s.value()) - S.result.scores[c].log_likelihood));
    }
    row("GF_fragment_terms", "sum of reconstructed per-fragment terms + exposure - engine total (all candidates)", worst, "<= 1e-6", worst <= 1e-6);
    // GB independent exhaustive oracle: each state scored alone (own interning, own call, one thread), prior in-loop
    for (std::size_t gi = 0; gi < small.size(); ++gi)
        for (bool r : {true, false}) {
            const auto lp = state_prior(res, pan, S.states[gi], r);
            const auto m = marginalise(S.logl[gi], lp, novel_counts(res, S.states[gi]));
            std::vector<long double> ot;
            for (std::size_t t = 0; t < S.states[gi].size(); ++t) {
                if (!std::isfinite(static_cast<double>(lp[t]))) continue;
                std::vector<std::string> cp;
                long double pr = 0;
                for (int d : S.states[gi][t]) {
                    cp.push_back(res.seq[d]);
                    // independent prior: C1 over the library entries, conditioned on the method's admitted set
                    std::vector<long double> num;
                    for (std::size_t k = 0; k < lib.window.at(200).size(); ++k) {
                        std::size_t mis = 0;
                        for (std::size_t i = 0; i < res.seq[d].size(); ++i) mis += res.seq[d][i] != lib.window.at(200)[k][i];
                        num.push_back(std::log(static_cast<long double>(lib.weights[k])) + (res.seq[d].size() - mis) * std::log1p(-static_cast<long double>(lib.eps_res)) +
                                      mis * std::log(static_cast<long double>(lib.eps_res) / 3.0L));
                    }
                    const CopySpace& sp = r ? res : pan;
                    std::vector<long double> den;
                    for (const auto& y : sp.seq) {
                        std::vector<long double> dd;
                        for (std::size_t k = 0; k < lib.window.at(200).size(); ++k) {
                            std::size_t mis = 0;
                            for (std::size_t i = 0; i < y.size(); ++i) mis += y[i] != lib.window.at(200)[k][i];
                            dd.push_back(std::log(static_cast<long double>(lib.weights[k])) + (y.size() - mis) * std::log1p(-static_cast<long double>(lib.eps_res)) +
                                         mis * std::log(static_cast<long double>(lib.eps_res) / 3.0L));
                        }
                        den.push_back(lse(dd));
                    }
                    pr += lse(num) - lse(den);
                }
                const std::string a = spell(w, small[gi], cp, true), b = spell(w, small[gi], cp, false);
                const panvar::InternedHaplotypes in = panvar::intern_haplotypes({panvar::MosaicPath{{0}}, panvar::MosaicPath{{1}}}, {a, b});
                ot.push_back(pr + panvar::score_diplotypes(in, {{in.index_of[0], in.index_of[1]}}, reads, library(), scoring(1)).scores[0].log_likelihood);
            }
            const double d = static_cast<double>(m.log_m - lse(ot));
            row("GB_oracle", small[gi].label() + (r ? " RES" : " PANEL") + ": prototype - independent exhaustive marginal", d, "|.| <= 1e-6", std::fabs(d) <= 1e-6);
        }
    // GC equivalent labels: library entry 0 split into two aliases of half weight
    {
        std::vector<std::string> wins = lib.window.at(200);
        std::vector<double> ws = lib.weights;
        wins.push_back(wins[0]);
        ws[0] /= 2;
        ws.push_back(ws[0]);
        const CopySpace ra = copy_space(wins, ws, lib.admitted.at(200), lib.novel.at(200), lib.eps_res, true);
        const CopySpace pa = copy_space(wins, ws, lib.admitted.at(200), lib.novel.at(200), lib.eps_res, false);
        double d = 0;
        for (std::size_t i = 0; i < res.seq.size(); ++i) d = std::max(d, std::fabs(static_cast<double>(res.log_prior[i] - ra.log_prior[i])));
        for (std::size_t i = 0; i < pan.seq.size(); ++i) d = std::max(d, std::fabs(static_cast<double>(pan.log_prior[i] - pa.log_prior[i])));
        row("GC_equivalent_labels", "alias split: per-sequence prior unchanged (distinct states " + std::to_string(ra.seq.size()) + ")", d, "<= 1e-12 and 8 states",
            d <= 1e-12 && ra.seq.size() == res.seq.size());
    }
    // GD base annotation: the same states enumerated from a different base annotation give the same marginal
    {
        std::vector<std::vector<int>> off = {{3, 5}, {7, 2, 6}};
        const Scored T = score_all(w, res, small, reads, off, 4);
        double d = 0;
        for (std::size_t gi = 0; gi < small.size(); ++gi)
            for (bool r : {true, false}) {
                const auto a = marginalise(S.logl[gi], state_prior(res, pan, S.states[gi], r), novel_counts(res, S.states[gi]));
                const auto b = marginalise(T.logl[gi], state_prior(res, pan, T.states[gi], r), novel_counts(res, T.states[gi]));
                d = std::max(d, std::fabs(static_cast<double>(a.log_m - b.log_m)));
            }
        row("GD_base_annotation", "re-annotated enumeration: marginal difference (1+1, 1+2; RES, PANEL)", d, "<= 1e-9", d <= 1e-9);
    }
    // GE negative control detects per-read independent marginalisation; GG KL identity; GH constant likelihood
    {
        const auto lp = state_prior(res, pan, S.states[1], true);
        const auto m = marginalise(S.logl[1], lp, novel_counts(res, S.states[1]));
        const double nc = static_cast<double>(negative_control(S, 1, lp, reads) - m.log_m);
        row("GE_negative_control", "1+2 RES: per-read marginalisation - correct marginal", nc, "|.| > 1e-3 (detected)", std::fabs(nc) > 1e-3);
        row("GG_kl_identity", "1+2 RES: log m - (E_q log L - KL)", static_cast<double>(m.residual), "|.| <= 1e-9", std::fabs(static_cast<double>(m.residual)) <= 1e-9);
        const std::vector<double> flat(S.states[1].size(), -2500.0);
        const auto c = marginalise(flat, lp, novel_counts(res, S.states[1]));
        row("GH_constant_likelihood", "1+2 RES (512 states): marginal - constant", static_cast<double>(c.log_m + 2500.0L), "|.| <= 1e-9 (no state-count penalty)",
            std::fabs(static_cast<double>(c.log_m + 2500.0L)) <= 1e-9);
    }
    std::cout << "test\titem\tvalue\texpectation\tresult\n" << g_rows.str() << "m6g_gates\t" << (g_failed ? "FAIL" : "PASS") << '\n';
    return g_failed ? 1 : 0;
}

int run(const std::string& libpath, const std::string& name, const std::string& out) {
    const auto t0 = std::chrono::steady_clock::now();
    const Library lib = read_library(libpath);
    Fixture fx;
    bool found = false;
    for (const auto& f : fixtures()) if (f.name == name) { fx = f; found = true; }
    if (!found) throw std::runtime_error("unknown fixture " + name);
    const World w = make_world(fx.seed);
    const auto reads = fixture_reads(w, lib, fx);
    const CopySpace res = copy_space(lib.window.at(fx.L), lib.weights, lib.admitted.at(fx.L), lib.novel.at(fx.L), lib.eps_res, true);
    const CopySpace pan = copy_space(lib.window.at(fx.L), lib.weights, lib.admitted.at(fx.L), lib.novel.at(fx.L), lib.eps_res, false);
    const Scored S = score_all(w, res, kCatalogue, reads, zero_offsets(kCatalogue), 4);
    std::ostringstream t;
    t << std::setprecision(12) << "fixture\tstructure\ttotal\tmethod\tstates_admitted\tlog_marginal\tposterior_mean_log_likelihood\tkl\tidentity_residual\tprior_mass\t"
                                  "top_state_posterior\teffective_states\tposterior_mean_novel_site_copies\tposterior_any_novel\tbest_state_log_likelihood_diagnostic\t"
                                  "negative_control_per_read\n";
    for (std::size_t gi = 0; gi < kCatalogue.size(); ++gi)
        for (const Method& me : {Method{"RES", true}, Method{"PANEL", false}}) {
            const auto lp = state_prior(res, pan, S.states[gi], me.residual);
            const auto m = marginalise(S.logl[gi], lp, novel_counts(res, S.states[gi]));
            std::size_t admitted = 0;
            double best = -INFINITY;
            for (std::size_t i = 0; i < lp.size(); ++i)
                if (std::isfinite(static_cast<double>(lp[i]))) { ++admitted; best = std::max(best, S.logl[gi][i]); }
            const long double nc = negative_control(S, gi, lp, reads);
            t << fx.name << '\t' << kCatalogue[gi].label() << '\t' << kCatalogue[gi].nA + kCatalogue[gi].nB << '\t' << me.name << '\t' << admitted << '\t'
              << static_cast<double>(m.log_m) << '\t' << static_cast<double>(m.elogl) << '\t' << static_cast<double>(m.kl) << '\t' << static_cast<double>(m.residual) << '\t'
              << static_cast<double>(m.mass) << '\t' << m.top << '\t' << m.eff << '\t' << static_cast<double>(m.novel_mean) << '\t' << static_cast<double>(m.any_novel) << '\t'
              << best << '\t' << static_cast<double>(nc) << '\n';
        }
    // truth state's prior under RES: cumulative residual cost relative to the same backbones without the novel alleles
    long double truth_lp = 0, backbone_lp = 0;
    int truth_novel = 0;
    for (auto [e, p] : fx.copies) {
        const int idx = e * 4 + p, idx0 = e * 4;
        truth_lp += res.log_prior[idx];
        backbone_lp += res.log_prior[idx0];
        truth_novel += __builtin_popcount(static_cast<unsigned>(p));
    }
    const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::ofstream(out + ".structures.tsv") << t.str();
    std::ofstream(out + ".run.tsv") << std::setprecision(12) << "key\tvalue\nfixture\t" << fx.name << "\nunit_bp\t" << fx.L << "\ntruth\t" << fx.truth.label()
                                    << "\nfragments\t" << reads.size() << "\nunique_haplotypes\t" << S.unique_haplotypes << "\ncandidates_scored\t" << S.cand_haps.size()
                                    << "\ndistinct_states_per_copy_RES\t" << res.seq.size() << "\ndistinct_states_per_copy_PANEL\t" << pan.seq.size()
                                    << "\ntruth_novel_site_copies\t" << truth_novel << "\ntruth_state_log_prior_RES\t" << static_cast<double>(truth_lp)
                                    << "\nsame_backbones_without_novel_log_prior_RES\t" << static_cast<double>(backbone_lp)
                                    << "\ncumulative_residual_prior_cost\t" << static_cast<double>(truth_lp - backbone_lp)
                                    << "\nper_copy_novel_pattern3_vs_0\t" << static_cast<double>(res.log_prior[3] - res.log_prior[0])
                                    << "\nseconds\t" << secs << "\npeak_rss_gib\t" << peak_rss_gib() << '\n';
    std::cout << fx.name << ": " << reads.size() << " fragments, " << S.cand_haps.size() << " candidates, " << secs << " s\n";
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const std::string mode = argc > 1 ? argv[1] : "";
        if (mode == "gates" && argc == 3) return gates(argv[2]);
        if (mode == "run" && argc == 5) return run(argv[2], argv[3], argv[4]);
        std::cerr << "usage: c1_structural gates <library.tsv> | run <library.tsv> <fixture> <out prefix>\n";
        return 2;
    } catch (const std::exception& e) {
        std::cerr << "c1_structural: ERROR " << e.what() << '\n';
        return 2;
    }
}
