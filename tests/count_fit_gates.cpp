// M6-A gates for the count-fit comparator core (tests/count_fit.hpp), on synthetic sequences with exact observed counts.
// Every check is a TSV row; any FAIL makes the exit status 1. No gate assumes that every structural difference is
// observable through 31-mers: G3 asserts the opposite for pooled-content-identical allocations.
#include "count_fit.hpp"

#include <iomanip>
#include <iostream>
#include <random>
#include <set>
#include <sstream>

namespace {
using namespace countfit;
std::size_t g_failed = 0;
std::ostringstream g_rows;
void row(const std::string& t, const std::string& item, double v, const std::string& expect, bool pass) {
    g_rows << t << '\t' << item << '\t' << std::setprecision(15) << v << '\t' << expect << '\t' << (pass ? "PASS" : "FAIL") << '\n';
    g_failed += pass ? 0 : 1;
}
std::mt19937_64 rng(6060);
std::string rnd(std::size_t n) {
    std::string s(n, 'A');
    for (char& c : s) c = "ACGT"[rng() % 4];
    return s;
}
std::string rep(const std::string& u, int n) {
    std::string s;
    for (int i = 0; i < n; ++i) s += u;
    return s;
}
Dictionary dict_of(const std::vector<std::string>& seqs) {
    Dictionary d;
    for (const auto& s : seqs) for_each_kmer(s, [&](std::uint64_t k, std::size_t) { d.add(k); });
    return d;
}
std::vector<double> exact_y(const Dictionary& d, const std::string& h1, const std::string& h2) {
    std::vector<double> y(d.kmer.size(), 0.0);
    for (const std::string* h : {&h1, &h2}) for_each_kmer(*h, [&](std::uint64_t k, std::size_t) { const long i = d.find(k); if (i >= 0) y[i] += 1.0; });
    return y;
}
struct Fixture {
    Context ctx;
    std::vector<std::string> alleles;
    Dictionary d;
    std::vector<Sparse> m1, m2;
    Sparse b1, b2;
    void build(const std::vector<std::string>& extra_panel = {}) {
        std::vector<std::string> panel = extra_panel;
        for (const auto& a : alleles) for (int h = 0; h < 2; ++h) panel.push_back(ctx.left[h] + a + ctx.right[h]);
        d = dict_of(panel);
        b1 = count_string(ctx.left[0], d); { Sparse r = count_string(ctx.right[0], d); merge(b1, r); }
        b2 = count_string(ctx.left[1], d); { Sparse r = count_string(ctx.right[1], d); merge(b2, r); }
        m1.clear(); m2.clear();
        for (const auto& a : alleles) {
            m1.push_back(count_string(middle_window(ctx.left[0], a, ctx.right[0]), d));
            m2.push_back(count_string(middle_window(ctx.left[1], a, ctx.right[1]), d));
        }
    }
    static void merge(Sparse& a, const Sparse& b) {
        std::map<std::uint32_t, std::uint32_t> m(a.v.begin(), a.v.end());
        for (const auto& [i, c] : b.v) m[i] += c;
        a.v.assign(m.begin(), m.end());
        a.hash += b.hash;
        a.outside += b.outside;
    }
    std::string hap(int h, std::size_t a) const { return ctx.left[h] + alleles[a] + ctx.right[h]; }
};
const Candidate& find(const std::vector<Candidate>& c, std::size_t a, std::size_t b) {
    for (const auto& x : c) if (x.a == a && x.b == b) return x;
    throw std::runtime_error("candidate not found");
}

void g1_multiplicities() {
    const std::string U = rnd(40);
    Fixture f;
    f.ctx.left[0] = f.ctx.left[1] = rnd(300);
    f.ctx.right[0] = rnd(300);
    f.ctx.right[1] = f.ctx.right[0];
    f.ctx.right[1][5] = f.ctx.right[1][5] == 'A' ? 'C' : 'A';   // homologue flanks differ near the junction
    f.alleles = {"", "ACGTA", rep(U, 3), rnd(90)};
    f.build();
    std::size_t bad = 0;
    for (int h = 0; h < 2; ++h) {
        for (std::size_t a = 0; a < f.alleles.size(); ++a) {
            Sparse whole = count_string(f.hap(h, a), f.d), parts = h == 0 ? f.b1 : f.b2;
            Fixture::merge(parts, h == 0 ? f.m1[a] : f.m2[a]);
            bad += whole.v != parts.v ? 1 : 0;
            bad += whole.hash != parts.hash ? 1 : 0;
        }
    }
    row("G1_multiplicity", "base + middle == whole-spelling counts (empty, 5-bp, tandem U^3, 90-bp alleles; both homologues)", static_cast<double>(bad), "0", bad == 0);
    // the tandem copy: a U k-mer that wraps a copy boundary occurs twice in U^3 (three copies, two boundaries) plus flank effects
    const Sparse t = count_string(rep(U, 3), f.d);
    std::uint32_t mx = 0;
    for (const auto& [i, c] : t.v) mx = std::max(mx, c);
    row("G1_multiplicity", "repeated copies counted (max multiplicity inside U^3)", mx, "3", mx == 3);
}

void g2_content_change() {
    const std::string U = rnd(45);
    Fixture f;
    f.ctx.left[0] = f.ctx.left[1] = rnd(250);
    f.ctx.right[0] = f.ctx.right[1] = rnd(250);
    f.alleles = {rep(U, 2), rep(U, 3)};
    f.build();
    const auto y = exact_y(f.d, f.hap(0, 0), f.hap(1, 0));   // truth (2, 2)
    const auto c = score_all(y, f.b1, f.b2, f.m1, f.m2);
    const auto& t = find(c, 0, 0), &alt = find(c, 1, 0);
    row("G2_content", "truth (2,2) loss", t.loss, "0", std::fabs(t.loss) < 1e-9);
    row("G2_content", "(3,2) loss > truth", alt.loss, "> 0 (total copy content change is detected)", alt.loss > 1.0 && alt.profile != t.profile);
}

void g3_pooled_allocation() {
    const std::string U = rnd(50);
    Fixture f;
    f.ctx.left[0] = f.ctx.left[1] = rnd(250);
    f.ctx.right[0] = f.ctx.right[1] = rnd(250);
    f.alleles = {rep(U, 1), rep(U, 2), rep(U, 3)};
    f.build();
    const auto y = exact_y(f.d, f.hap(0, 1), f.hap(1, 1));   // truth 2+2
    const auto c = score_all(y, f.b1, f.b2, f.m1, f.m2);
    const auto& two = find(c, 1, 1), &tri = find(c, 2, 0), &tri2 = find(c, 0, 2);
    row("G3_allocation", "2+2 and 3+1 same profile hash", static_cast<double>(two.profile == tri.profile && tri.profile == tri2.profile), "1", two.profile == tri.profile && tri.profile == tri2.profile);
    row("G3_allocation", "|loss(3+1) - loss(2+2)|", std::fabs(tri.loss - two.loss), "<= 1e-9 (identical count scores)", std::fabs(tri.loss - two.loss) <= 1e-9);
    const auto groups = group_profiles(c, f.m1, f.m2);
    bool together = false;
    for (const auto& g : groups) {
        std::set<std::pair<std::size_t, std::size_t>> m;
        for (std::size_t i : g) m.insert({c[i].a, c[i].b});
        if (m.count({1, 1}) && m.count({2, 0}) && m.count({0, 2})) together = g.size() == 3 && g == groups.front();
    }
    row("G3_allocation", "best group holds 2+2 and both 3+1 assignments: allocation UNRESOLVED", static_cast<double>(together), "1", together);
}

void g4_junction() {
    const std::string U = rnd(50);
    std::string V = U;
    V[25] = V[25] == 'A' ? 'G' : 'A';
    Fixture f;
    f.ctx.left[0] = f.ctx.left[1] = rnd(250);
    f.ctx.right[0] = f.ctx.right[1] = rnd(250);
    f.alleles = {U + V, V + U, U + U};   // same unit multiset {U,V} vs {V,U}: only junction content differs
    f.build();
    const auto y = exact_y(f.d, f.hap(0, 0), f.hap(1, 2));   // truth (UV, UU)
    const auto c = score_all(y, f.b1, f.b2, f.m1, f.m2);
    const auto& t = find(c, 0, 2), &alt = find(c, 1, 2);
    row("G4_junction", "(VU,UU) profile differs from truth (UV,UU)", static_cast<double>(alt.profile != t.profile), "1", alt.profile != t.profile);
    row("G4_junction", "truth loss 0 < arrangement-swap loss", alt.loss, "> 0 (junction content discriminates)", std::fabs(t.loss) < 1e-9 && alt.loss > 1.0);
}

void g5_alias() {
    Fixture f, g;
    f.ctx.left[0] = f.ctx.left[1] = rnd(200);
    f.ctx.right[0] = f.ctx.right[1] = rnd(200);
    const std::string a = rnd(80), b = rnd(80);
    f.alleles = {a, b};
    g.ctx = f.ctx;
    g.alleles = {a, b, a};   // index 2 is a label alias of index 0
    f.build();
    g.build();
    const auto y = exact_y(f.d, f.hap(0, 0), f.hap(1, 1));
    const auto yg = exact_y(g.d, g.hap(0, 0), g.hap(1, 1));
    const auto cf = score_all(y, f.b1, f.b2, f.m1, f.m2), cg = score_all(yg, g.b1, g.b2, g.m1, g.m2);
    row("G5_alias", "dictionary unchanged by an alias", static_cast<double>(f.d.kmer.size() == g.d.kmer.size()), "1", f.d.kmer.size() == g.d.kmer.size());
    row("G5_alias", "alias candidate has the original's loss and profile (evidence not multiplied)", std::fabs(find(cg, 2, 1).loss - find(cf, 0, 1).loss),
        "0", find(cg, 2, 1).loss == find(cg, 0, 1).loss && find(cg, 2, 1).profile == find(cg, 0, 1).profile && std::fabs(find(cg, 0, 1).loss - find(cf, 0, 1).loss) < 1e-9);
}

void g6_oracle_and_order() {
    const std::string U = rnd(40);
    Fixture f;
    f.ctx.left[0] = rnd(300);
    f.ctx.left[1] = f.ctx.left[0];
    f.ctx.left[1][290] = f.ctx.left[1][290] == 'C' ? 'T' : 'C';
    f.ctx.right[0] = f.ctx.right[1] = rnd(300);
    f.alleles = {"", rep(U, 2), rep(U, 3), rnd(70), rnd(35) + U};
    f.build({rnd(400)});   // an extra panel spelling: dictionary k-mers no candidate carries
    std::vector<double> y = exact_y(f.d, f.hap(0, 2), f.hap(1, 3));
    std::normal_distribution<double> noise(0.0, 0.3);
    for (double& x : y) x = std::max(0.0, x + noise(rng));   // observed counts with noise
    const auto c = score_all(y, f.b1, f.b2, f.m1, f.m2);
    double worst = 0.0;
    for (const auto& x : c) worst = std::max(worst, std::fabs(x.loss - naive_loss(y, f.hap(0, x.a), f.hap(1, x.b), f.d)) / std::max(1.0, x.loss));
    row("G6_oracle", "incremental loss vs independent scalar loss (max relative diff, all candidates)", worst, "<= 1e-12", worst <= 1e-12);
    // order invariance: permute the alleles; every (sequence on h1, sequence on h2) keeps its loss
    std::vector<std::size_t> perm{3, 0, 4, 2, 1};
    Fixture p = f;
    p.alleles.clear();
    for (std::size_t i : perm) p.alleles.push_back(f.alleles[i]);
    p.build({});
    p.d = f.d;   // same dictionary
    p.b1 = count_string(p.ctx.left[0], p.d); { Sparse r = count_string(p.ctx.right[0], p.d); Fixture::merge(p.b1, r); }
    p.b2 = count_string(p.ctx.left[1], p.d); { Sparse r = count_string(p.ctx.right[1], p.d); Fixture::merge(p.b2, r); }
    p.m1.clear(); p.m2.clear();
    for (const auto& a : p.alleles) { p.m1.push_back(count_string(middle_window(p.ctx.left[0], a, p.ctx.right[0]), p.d)); p.m2.push_back(count_string(middle_window(p.ctx.left[1], a, p.ctx.right[1]), p.d)); }
    const auto cp = score_all(y, p.b1, p.b2, p.m1, p.m2);
    double od = 0.0;
    for (const auto& x : cp) od = std::max(od, std::fabs(x.loss - find(c, perm[x.a], perm[x.b]).loss));
    row("G6_order", "candidate-order invariance (max abs diff)", od, "0", od == 0.0);
}
// M6-B: exact grouping under FORCED hash collisions, and complete-haplotype swap invariance.
void g7_exact_grouping() {
    const std::string U = rnd(50);
    std::string V = U;
    V[20] = V[20] == 'A' ? 'G' : 'A';
    Fixture f;
    f.ctx.left[0] = f.ctx.left[1] = rnd(250);
    f.ctx.right[0] = f.ctx.right[1] = rnd(250);
    f.alleles = {rep(U, 1), rep(U, 2), rep(U, 3), U + V, V + U};
    f.build();
    const auto y = exact_y(f.d, f.hap(0, 1), f.hap(1, 1));
    const auto c = score_all(y, f.b1, f.b2, f.m1, f.m2);
    std::size_t col_default = 0, col_forced = 0;
    const auto normal = group_profiles(c, f.m1, f.m2, ~std::uint64_t(0), &col_default);
    const auto forced = group_profiles(c, f.m1, f.m2, 0, &col_forced);   // every candidate shares one hash bucket
    // ground truth: distinct complete profiles by direct comparison
    std::vector<Sparse> profiles;
    for (const auto& x : c) {
        const Sparse p = sum_sparse(f.m1[x.a], f.m2[x.b]);
        if (std::none_of(profiles.begin(), profiles.end(), [&](const Sparse& q) { return q.v == p.v; })) profiles.push_back(p);
    }
    const auto as_sets = [&](const std::vector<std::vector<std::size_t>>& g) {
        std::set<std::set<std::size_t>> s;
        for (const auto& v : g) s.insert(std::set<std::size_t>(v.begin(), v.end()));
        return s;
    };
    row("G7_exact_grouping", "forced collisions: groups == distinct complete profiles", static_cast<double>(forced.size()), std::to_string(profiles.size()),
        forced.size() == profiles.size());
    row("G7_exact_grouping", "forced collisions: same partition as the default hash", static_cast<double>(as_sets(forced) == as_sets(normal)), "1", as_sets(forced) == as_sets(normal));
    row("G7_exact_grouping", "forced bucket reports the collision; default reports none", static_cast<double>(col_forced * 10 + col_default), "10", col_forced == 1 && col_default == 0);
    bool equal_grouped = false;
    for (const auto& g : forced) {
        std::set<std::pair<std::size_t, std::size_t>> m;
        for (std::size_t i : g) m.insert({c[i].a, c[i].b});
        equal_grouped = equal_grouped || (m.count({1, 1}) && m.count({2, 0}) && m.count({0, 2}));
    }
    row("G7_exact_grouping", "equal profiles (2+2, 3+1, 1+3) still grouped under forced collisions", static_cast<double>(equal_grouped), "1", equal_grouped);
}
void g8_haplotype_swap() {
    Fixture f, g;
    f.ctx.left[0] = rnd(300);
    f.ctx.left[1] = rnd(300);
    f.ctx.right[0] = rnd(300);
    f.ctx.right[1] = rnd(300);
    f.alleles = {rnd(60), rnd(80), ""};
    g.ctx.left[0] = f.ctx.left[1];
    g.ctx.left[1] = f.ctx.left[0];
    g.ctx.right[0] = f.ctx.right[1];
    g.ctx.right[1] = f.ctx.right[0];
    g.alleles = f.alleles;
    f.build();
    g.build();
    const auto y = exact_y(f.d, f.hap(0, 0), f.hap(1, 1));
    const auto cf = score_all(y, f.b1, f.b2, f.m1, f.m2);
    std::vector<double> yg(g.d.kmer.size(), 0.0);
    for (std::size_t i = 0; i < f.d.kmer.size(); ++i) yg[static_cast<std::size_t>(g.d.find(f.d.kmer[i]))] = y[i];
    const auto cg = score_all(yg, g.b1, g.b2, g.m1, g.m2);
    double worst = 0.0;
    bool profiles = true;
    for (const auto& x : cf) {
        const auto& z = find(cg, x.b, x.a);   // the complete haplotypes swapped: (ctx1+a, ctx2+b) == (ctx2+b, ctx1+a)
        worst = std::max(worst, std::fabs(x.loss - z.loss));
        profiles = profiles && x.profile == z.profile;
    }
    row("G8_haplotype_swap", "swapping the two complete haplotypes leaves loss and profile unchanged", worst, "0", worst == 0.0 && profiles);
    // target-allele swap with fixed, different flanks is NOT invariant in general (junctions differ)
    row("G8_haplotype_swap", "target swap (a,b) vs (b,a) under different flanks changes the profile", static_cast<double>(find(cf, 0, 1).profile != find(cf, 1, 0).profile),
        "1 (junction counts differ)", find(cf, 0, 1).profile != find(cf, 1, 0).profile);
}
}  // namespace

int main() {
    g1_multiplicities();
    g2_content_change();
    g3_pooled_allocation();
    g4_junction();
    g5_alias();
    g6_oracle_and_order();
    g7_exact_grouping();
    g8_haplotype_swap();
    std::cout << "test\titem\tvalue\texpectation\tresult\n" << g_rows.str() << "count_fit_gates\t" << (g_failed ? "FAIL" : "PASS") << '\n';
    return g_failed ? 1 : 0;
}
