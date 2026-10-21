// M13-A (DIAGNOSTIC): gates for the variation-tolerant truncated marginal (tests/vt_marginal.hpp) on small synthetic
// worlds in M6-G's style (synthetic flanks, tandem units, 100-bp reads at diploid depth 20, sequencing error 0.002).
// The exact oracle is M6-G's GB construction: every joint admitted state scored ALONE by score_diplotypes, the full C1
// mixture prior recomputed per state (never factorised), summed exhaustively.
//
//   V1 fixed content: no admitted site -> base log likelihood = an independent score_diplotypes call; log m_res = log m_panel
//   V2 constant likelihood (no fragments): log m_res - base = the analytic prior mass of the admitted set (no state-count bonus)
//   V3 equivalent labels: splitting a component into two half-weight aliases leaves every value unchanged
//   V4 homologue exchange: (A,B) and (B,A) give identical values
//   V5 shared state: the per-read independent marginalisation (negative control) differs from the shared-state marginal
//   V6 exactness: engine = brute-force exact enumeration, including sites coupled across identical tandem copies
//   V7 informative unbalanced arrangement (truth 3+1 with variation absent from the panel): the marginal prefers 3+1 over 2+2
//   V8 uninformative arrangement (identical flanks, units longer than inserts, equal component weights): 2+2 and 1+3 tie
#include "vt_marginal.hpp"

#include <chrono>
#include <iomanip>
#include <iostream>
#include <random>

namespace {

std::size_t g_checks = 0, g_failed = 0;
void check(bool ok, const std::string& what) {
    ++g_checks;
    std::cout << (ok ? "PASS " : "FAIL ") << what << '\n';
    if (!ok) ++g_failed;
}

panvar::LibraryModel library() {
    panvar::LibraryModel lib;
    lib.insert_mean_bp = 350;
    lib.insert_sd_bp = 35;
    lib.support_sigmas = 4;
    lib.substitution_error = 0.002;
    lib.outlier_weight = 0.001;
    lib.fragment_start_intensity = 20.0 / (4.0 * 100.0);
    lib.provenance = "M13-A synthetic declared";
    return lib;
}
panvar::FragmentScoringOptions options() {
    panvar::FragmentScoringOptions o;
    o.threads = 4;
    o.mismatch_band = 6;
    o.minimum_piece_bp = 12;
    return o;
}
std::string rnd(std::mt19937_64& g, std::size_t n) {
    std::string s(n, 'A');
    for (char& c : s) c = "ACGT"[g() % 4];
    return s;
}
char alt_of(char c) { return c == 'A' ? 'C' : c == 'C' ? 'G' : c == 'G' ? 'T' : 'A'; }
void simulate(const std::string& hap, std::mt19937_64& rng, std::vector<panvar::ReadPair>& out, const std::string& tag) {
    std::poisson_distribution<int> starts(20.0 / (4.0 * 100.0));
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
            std::string m1 = frag.substr(0, 100), m2 = vt::revcomp(frag.substr(len - 100));
            if (u(rng) < 0.5) std::swap(m1, m2);
            const std::string name = tag + std::to_string(out.size());
            out.push_back({{name + "/1", noisy(m1), std::string(100, 'I')}, {name + "/2", noisy(m2), std::string(100, 'I')}});
        }
    }
}

struct World {
    std::string fa_l, fa_r, fb_l, fb_r;
};
vt::Candidate make(const std::string& label, const World& w, const std::string& block_a, const std::string& block_b,
                   const std::vector<vt::Component>& comps) {
    vt::Candidate c;
    c.label = label;
    c.hap = {w.fa_l + block_a + w.fa_r, w.fb_l + block_b + w.fb_r};
    c.block = {panvar::SequenceInterval{w.fa_l.size(), w.fa_l.size() + block_a.size()},
               panvar::SequenceInterval{w.fb_l.size(), w.fb_l.size() + block_b.size()}};
    c.comps = {comps, comps};
    return c;
}
std::string units(const std::vector<std::string>& v) {
    std::string s;
    for (const auto& x : v) s += x;
    return s;
}

// Exact oracle: every joint state over the candidate's admitted sites, each diplotype scored alone, full mixture prior.
double brute_force(const vt::Candidate& c, const std::vector<vt::Site>& sites, const std::vector<panvar::ReadPair>& frags, double eps) {
    std::vector<std::pair<std::string, std::string>> pairs;
    std::vector<std::array<std::string, 2>> blocks;
    for (unsigned m = 0; m < (1U << sites.size()); ++m) {
        std::array<std::string, 2> h = c.hap;
        for (std::size_t j = 0; j < sites.size(); ++j) if ((m >> j) & 1U) h[sites[j].hom][sites[j].pos] = sites[j].alt;
        pairs.push_back({h[0], h[1]});
        blocks.push_back({h[0].substr(c.block[0].begin, c.block[0].size()), h[1].substr(c.block[1].begin, c.block[1].size())});
    }
    std::vector<double> terms;
    for (std::size_t i = 0; i < pairs.size(); ++i) {
        const auto sc = vt::score_pairs({pairs[i]}, frags, library(), options());   // each state scored alone
        terms.push_back(sc[0].first + sc[0].second + vt::log_c1(blocks[i][0], c.comps[0], eps) + vt::log_c1(blocks[i][1], c.comps[1], eps));
    }
    return vt::lse(terms);
}

}  // namespace

int main() {
    const auto t0 = std::chrono::steady_clock::now();
    std::cout << std::setprecision(12);
    std::mt19937_64 g(20261021);
    World w{rnd(g, 800), rnd(g, 800), rnd(g, 800), rnd(g, 800)};
    const std::string U = rnd(g, 500);
    std::string Uv = U;   // the same unit carrying one substitution absent from every panel component
    Uv[250] = alt_of(Uv[250]);
    const vt::Settings S;
    const double eps = S.eps_res;
    const std::vector<vt::Component> panel = {{units({U}), std::log(1.0 / 3)}, {units({U, U}), std::log(1.0 / 3)}, {units({U, U, U}), std::log(1.0 / 3)}};

    // ---- world 1: truth A = 3 copies (copy 2 carries the novel substitutions), B = 1 copy; distinct flanks ------------
    std::vector<panvar::ReadPair> reads;
    std::mt19937_64 rng(7);
    simulate(w.fa_l + units({U, Uv, U}) + w.fa_r, rng, reads, "a");
    simulate(w.fb_l + units({U}) + w.fb_r, rng, reads, "b");
    const vt::Candidate c31 = make("3+1", w, units({U, U, U}), units({U}), panel);
    const vt::Candidate c22 = make("2+2", w, units({U, U}), units({U, U}), panel);
    const vt::Candidate c13 = make("1+3", w, units({U}), units({U, U, U}), panel);
    const auto scan = vt::base_scan({c31, c22, c13}, reads, library(), options());
    const auto r = vt::evaluate({c31, c22, c13}, reads, library(), options(), S, scan);
    for (const auto& x : r)
        std::cout << "  world1 " << x.label << ": status " << x.status << " sites " << x.sites << " groups " << x.groups << " max_group " << x.max_group
                  << " base " << x.base_loglik << " m_panel " << x.log_m_panel << " m_res " << x.log_m_res << " fit_gain " << x.eq_fit_gain
                  << " kl " << x.kl << " interaction " << x.interaction_error << '\n';

    // V1 fixed content
    {
        vt::Settings f = S;
        f.admit = false;
        const auto rf = vt::evaluate({c31}, reads, library(), options(), f, scan);
        const auto direct = vt::score_pairs({{c31.hap[1], c31.hap[0]}}, reads, library(), options());
        check(std::fabs(rf[0].base_loglik - (direct[0].first + direct[0].second)) <= 1e-9 && std::fabs(rf[0].log_m_res - rf[0].log_m_panel) <= 1e-9,
              "V1 fixed content: base log likelihood = an independent scorer call; m_res = m_panel");
    }
    // V2 constant likelihood: forced sites, no fragments
    {
        vt::Settings f = S;
        f.forced_hom = {0, 0, 1};
        f.forced_pos = {c31.block[0].begin + 120, c31.block[0].begin + 620, c31.block[1].begin + 250};
        f.forced_alt = {alt_of(c31.hap[0][f.forced_pos[0]]), alt_of(c31.hap[0][f.forced_pos[1]]), alt_of(c31.hap[1][f.forced_pos[2]])};
        const std::vector<panvar::ReadPair> none;
        const auto sc0 = vt::base_scan({c31}, none, library(), options());
        const auto rc = vt::evaluate({c31}, none, library(), options(), f, sc0);
        // analytic: sum over states of the full mixture prior, by direct enumeration of the mixture (8 states)
        std::vector<double> t;
        for (unsigned m = 0; m < 8; ++m) {
            std::array<std::string, 2> h = c31.hap;
            for (int j = 0; j < 3; ++j) if ((m >> j) & 1U) h[f.forced_hom[j]][f.forced_pos[j]] = f.forced_alt[j];
            t.push_back(vt::log_c1(h[0].substr(c31.block[0].begin, c31.block[0].size()), panel, eps) +
                        vt::log_c1(h[1].substr(c31.block[1].begin, c31.block[1].size()), panel, eps));
        }
        const double mass = vt::lse(t);
        check(std::fabs((rc[0].log_m_res - rc[0].base_loglik) - mass) <= 1e-9 && mass <= 0.0 && rc[0].groups == 3,
              "V2 constant likelihood: marginal - base = analytic admitted prior mass (" + std::to_string(mass) + " <= 0; no state-count bonus)");
    }
    // V3 alias split
    {
        std::vector<vt::Component> split = panel;
        split[2].log_w = std::log(1.0 / 6);
        split.push_back({units({U, U, U}), std::log(1.0 / 6)});
        vt::Candidate a = c31;
        a.comps = {split, split};
        const auto ra = vt::evaluate({a}, reads, library(), options(), S, scan);
        check(std::fabs(ra[0].log_m_res - r[0].log_m_res) <= 1e-9 && std::fabs(ra[0].log_m_panel - r[0].log_m_panel) <= 1e-9,
              "V3 equivalent labels: splitting a component into two half-weight aliases changes nothing");
    }
    // V4 homologue exchange
    {
        vt::Candidate x = c31;
        std::swap(x.hap[0], x.hap[1]);
        std::swap(x.block[0], x.block[1]);
        std::swap(x.comps[0], x.comps[1]);
        const auto sx = vt::base_scan({x}, reads, library(), options());
        const auto rx = vt::evaluate({x}, reads, library(), options(), S, sx);
        check(std::fabs(rx[0].log_m_res - r[0].log_m_res) <= 1e-6 && std::fabs(rx[0].base_loglik - r[0].base_loglik) <= 1e-9,
              "V4 homologue exchange: (A,B) and (B,A) identical");
    }
    // V6 exactness against brute force (world 1, 3+1 candidate: sites across tandem copies)
    {
        const double exact = brute_force(c31, r[0].site_list, reads, eps);
        std::cout << "  V6: engine " << r[0].log_m_res << " brute force " << exact << " (sites " << r[0].sites << ", groups " << r[0].groups << ")\n";
        check(r[0].sites >= 1 && r[0].sites <= 10 && std::fabs(r[0].log_m_res - exact) <= 1e-6,
              "V6 exactness: engine = brute-force enumeration over every joint state (" + std::to_string(r[0].sites) + " sites)");
    }
    // V5 shared state vs per-read independent marginalisation (negative control, never used for inference)
    {
        const auto& sites = r[0].site_list;
        std::vector<std::string> h0s, h1s;
        std::vector<double> lp;
        for (unsigned m = 0; m < (1U << sites.size()); ++m) {
            std::array<std::string, 2> h = c31.hap;
            for (std::size_t j = 0; j < sites.size(); ++j) if ((m >> j) & 1U) h[sites[j].hom][sites[j].pos] = sites[j].alt;
            h0s.push_back(h[0]);
            h1s.push_back(h[1]);
            lp.push_back(vt::log_c1(h[0].substr(c31.block[0].begin, c31.block[0].size()), panel, eps) +
                         vt::log_c1(h[1].substr(c31.block[1].begin, c31.block[1].size()), panel, eps));
        }
        const double lpz = vt::lse(lp);
        std::vector<panvar::MosaicPath> paths;
        std::vector<std::string> seqs;
        for (std::size_t i = 0; i < h0s.size(); ++i) {
            paths.push_back(panvar::MosaicPath{{999999, 2 * i}});
            seqs.push_back(h0s[i]);
            paths.push_back(panvar::MosaicPath{{999999, 2 * i + 1}});
            seqs.push_back(h1s[i]);
        }
        const auto in = panvar::intern_haplotypes(paths, seqs);
        std::vector<panvar::HaplotypeCandidate> cands;
        for (std::size_t i = 0; i < h0s.size(); ++i) cands.push_back({in.index_of[2 * i], in.index_of[2 * i + 1]});
        const auto sc = panvar::score_diplotypes(in, cands, reads, library(), options());
        const std::size_t H = sc.unique_haplotypes;
        const panvar::LibraryModel L = library();
        double nc = 0;
        for (std::size_t f = 0; f < reads.size(); ++f) {
            std::vector<double> v;
            const double bg = panvar::log_background_probability(reads[f]);
            for (std::size_t i = 0; i < h0s.size(); ++i) {
                const double m = panvar::log_add_exp(sc.log_mass[f * H + in.index_of[2 * i]], sc.log_mass[f * H + in.index_of[2 * i + 1]]);
                v.push_back(lp[i] - lpz + panvar::log_fragment_contribution(m, bg, L));
            }
            nc += vt::lse(v);
        }
        std::vector<double> ex;
        for (std::size_t i = 0; i < h0s.size(); ++i) ex.push_back(lp[i] - lpz + sc.scores[i].exposure_term);
        nc += vt::lse(ex);
        const double shared = r[0].log_m_res - lpz;   // the shared-state marginal under the same normalised prior
        std::cout << "  V5: shared-state " << shared << " per-read independent " << nc << '\n';
        check(std::fabs(nc - shared) > 1e-3, "V5 shared state: per-read independent marginalisation differs from the shared-state marginal");
    }
    // V7 informative unbalanced arrangement
    {
        std::cout << "  V7: m_res 3+1 " << r[0].log_m_res << " vs 2+2 " << r[1].log_m_res << "; fixed-content (panel) 3+1 " << r[0].log_m_panel << " vs 2+2 "
                  << r[1].log_m_panel << '\n';
        check(r[0].status == "OK" && r[1].status == "OK" && r[0].log_m_res > r[1].log_m_res && r[0].log_m_res > r[2].log_m_res,
              "V7 informative unbalanced arrangement: the marginal prefers the true 3+1");
        check(std::fabs(r[0].interaction_error) <= 1e-6 && std::fabs(r[1].interaction_error) <= 1e-6, "V7 joint MAP state is additive across groups");
    }
    // ---- world 2: identical flanks on both homologues, identical units; truth 2+2 --------------------------------------
    {
        World w2{w.fa_l, w.fa_r, w.fa_l, w.fa_r};
        std::vector<panvar::ReadPair> r2;
        std::mt19937_64 rng2(11);
        simulate(w2.fa_l + units({U, U}) + w2.fa_r, rng2, r2, "a");
        simulate(w2.fb_l + units({U, U}) + w2.fb_r, rng2, r2, "b");
        const vt::Candidate a22 = make("2+2", w2, units({U, U}), units({U, U}), panel);
        const vt::Candidate a13 = make("1+3", w2, units({U}), units({U, U, U}), panel);
        const auto s2 = vt::base_scan({a22, a13}, r2, library(), options());
        const auto q = vt::evaluate({a22, a13}, r2, library(), options(), S, s2);
        std::cout << "  V8: m_res 2+2 " << q[0].log_m_res << " 1+3 " << q[1].log_m_res << " (sites " << q[0].sites << "/" << q[1].sites << ")\n";
        check(q[0].status == "OK" && q[1].status == "OK" && std::fabs(q[0].log_m_res - q[1].log_m_res) <= 1e-6,
              "V8 uninformative arrangement remains ambiguous (2+2 and 1+3 tie)");
    }
    const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::cout << (g_failed == 0 ? "ALL GATES PASS" : "GATES FAILED") << ": " << g_checks - g_failed << " of " << g_checks << " checks, " << secs << " s\n";
    return g_failed == 0 ? 0 : 1;
}
