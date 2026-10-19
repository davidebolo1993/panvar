// M4-A: invariance gates for the proposed distance-aware copying-prior contract (diagnostic; not a production test).
//
// The contract under test (tests/analysis/m4a/CONTRACT.md) lives ONLY here, as a prototype:
//   * geometry: each block's copying state sits at the midpoint of its polymorphic core on the reference axis -- the
//     reference allele's interval along the canonical chain with the maximal prefix/suffix shared by every allele of the
//     frame removed (a monomorphic block has an empty core at its interval's midpoint); d_b = core_mid(b) - core_mid(b-1);
//   * transition: T_b = (1 - r_b) I + r_b U with r_b = 1 - exp(-d_b / lambda), lambda a FIXED copying length (bp) and U the
//     declared redraw measure (uniform over the n template rows, the current one included);
//   * the existing library recursions (log_haploid_mosaic_prior / log_haploid_sequence_prior) consume the per-boundary r.
// The independent oracle is tests/genotype_mosaic_core.cpp's brute_haploid_prior (every template history enumerated),
// included unchanged as M3-A did. Geometry, lambda and U are frozen across every invariance comparison. Each check is a TSV
// row; any FAIL makes the exit status 1. "CONTROL" rows document a dependence the contract removes (legacy constant r,
// M3-A panel-mean span, untrimmed block midpoints); they assert the dependence exists, so they PASS when it does.
#define main genotype_mosaic_core_main
#include "genotype_mosaic_core.cpp"
#undef main

#include <iomanip>
#include <random>
#include <set>
#include <sstream>

namespace {

std::size_t g_failed = 0;
std::ostringstream g_rows;
constexpr double kLambda = 40.0;   // fixed copying length for every gate (bp); never re-matched after an edit

void row(const std::string& test, const std::string& item, double value, const std::string& expectation, bool pass) {
    g_rows << test << '\t' << item << '\t' << std::setprecision(15) << value << '\t' << expectation << '\t' << (pass ? "PASS" : "FAIL") << '\n';
    g_failed += pass ? 0 : 1;
}

// ---- the prototype contract ------------------------------------------------------------------------------------------
struct Geometry {
    std::vector<double> position;   // per block: core midpoint on the reference axis (bp)
    std::vector<double> distance;   // per boundary
};

// `ref` = the reference allele id per block (canonical chain order); `frame` = per block, every allele sequence of the
// FRAME (the alleles the core is trimmed against; by default the model's own alleles). `core` false = untrimmed block
// midpoints (a control).
Geometry reference_geometry(const panvar::PreparedGenotypeModel& m, const std::vector<std::size_t>& ref,
                            const std::vector<std::vector<std::string>>* frame = nullptr, bool core = true) {
    Geometry g;
    double x = static_cast<double>(m.invariant_segments()[0].size());
    for (std::size_t b = 0; b < m.blocks().size(); ++b) {
        std::vector<std::string> al;
        if (frame) al = (*frame)[b];
        else for (const auto& a : m.blocks()[b].alleles) al.push_back(a.sequence);
        const std::string& rs = m.blocks()[b].alleles[ref[b]].sequence;
        const double L = static_cast<double>(rs.size());
        std::size_t pre = rs.size(), suf = rs.size();
        if (al.size() > 1) {
            for (const auto& a : al) {
                std::size_t p = 0, q = 0;
                while (p < rs.size() && p < a.size() && rs[p] == a[p]) ++p;
                while (q < rs.size() && q < a.size() && rs[rs.size() - 1 - q] == a[a.size() - 1 - q]) ++q;
                pre = std::min(pre, p);
                suf = std::min(suf, q);
            }
        }
        double lo = x + (core ? std::min<double>(pre, L) : 0.0), hi = x + L - (core ? std::min<double>(suf, L) : 0.0);
        if (lo > hi) lo = hi = 0.5 * (lo + hi);
        g.position.push_back(0.5 * (lo + hi));
        x += L + static_cast<double>(m.invariant_segments()[b + 1].size());
    }
    for (std::size_t b = 1; b < g.position.size(); ++b) g.distance.push_back(g.position[b] - g.position[b - 1]);
    return g;
}

std::vector<double> hazard(const std::vector<double>& d, double lambda = kLambda) {
    std::vector<double> r;
    for (double x : d) r.push_back(x > 0.0 ? -std::expm1(-x / lambda) : 0.0);
    return r;
}

// M3-A's panel-mean span (start of b-1 to start of b, averaged over template rows): a control only.
std::vector<double> panel_mean(const panvar::PreparedGenotypeModel& m) {
    std::vector<double> d;
    for (std::size_t b = 1; b < m.blocks().size(); ++b) {
        double mean = 0.0;
        for (const auto& t : m.panel().template_alleles) mean += static_cast<double>(m.blocks()[b - 1].alleles[t[b - 1]].sequence.size());
        d.push_back(mean / m.panel().template_count() + static_cast<double>(m.invariant_segments()[b].size()));
    }
    return d;
}

// ---- oracles ---------------------------------------------------------------------------------------------------------
std::vector<std::vector<std::size_t>> all_paths(const panvar::PreparedGenotypeModel& m) {
    std::vector<std::vector<std::size_t>> paths{{}};
    for (const auto& block : m.blocks()) {
        std::vector<std::vector<std::size_t>> next;
        for (const auto& p : paths) for (std::size_t a = 0; a < block.alleles.size(); ++a) { auto q = p; q.push_back(a); next.push_back(q); }
        paths = next;
    }
    return paths;
}

// sequence-level oracle: the sum over EVERY allele path spelling the bytes of the brute-force history sum
double oracle_sequence(const panvar::PreparedGenotypeModel& m, const std::string& s, const std::vector<double>& r) {
    long double total = 0.0L;
    for (const auto& p : all_paths(m)) {
        if (panvar::spell_mosaic(m, panvar::MosaicPath{p}).bases == s) total += brute_haploid_prior(m.panel(), panvar::MosaicPath{p}, r);
    }
    return total > 0 ? std::log(static_cast<double>(total)) : -std::numeric_limits<double>::infinity();
}

std::set<std::string> spellable(const panvar::PreparedGenotypeModel& m) {
    std::set<std::string> out;
    for (const auto& p : all_paths(m)) out.insert(panvar::spell_mosaic(m, panvar::MosaicPath{p}).bases);
    return out;
}

double diff(double a, double b) { return (std::isinf(a) && std::isinf(b)) ? 0.0 : a - b; }

std::vector<std::size_t> ref_ids(const panvar::PreparedGenotypeModel& m, const std::string& name) {
    for (std::size_t t = 0; t < m.panel().template_count(); ++t) {
        if (m.panel().template_names[t] == name) return std::vector<std::size_t>(m.panel().template_alleles[t].begin(), m.panel().template_alleles[t].end());
    }
    throw std::runtime_error("no reference template " + name);
}

panvar::PanelTemplateInput tpl(const std::string& n, std::vector<std::string> ids) {
    return {n, panvar::TemplateFrameStatus::Complete, std::move(ids), panvar::ChainOrientation::Forward, ""};
}

// Compare every spellable sequence's sequence-level prior between two models (same templates), given their vectors.
void same_priors(const std::string& test, const std::string& what, const panvar::PreparedGenotypeModel& a, const std::vector<double>& ra,
                 const panvar::PreparedGenotypeModel& b, const std::vector<double>& rb) {
    const auto sa = spellable(a), sb = spellable(b);
    row(test, what + ":same_spellable_set", static_cast<double>(sa == sb), "1", sa == sb);
    double worst = 0.0;
    for (const auto& s : sa) worst = std::max(worst, std::fabs(diff(panvar::log_haploid_sequence_prior(a, s, ra), panvar::log_haploid_sequence_prior(b, s, rb))));
    row(test, what + ":max_abs_prior_difference", worst, "<= 1e-12", worst <= 1e-12);
}

double max_prior_change(const panvar::PreparedGenotypeModel& a, const std::vector<double>& ra, const panvar::PreparedGenotypeModel& b,
                        const std::vector<double>& rb) {
    double worst = 0.0;
    for (const auto& s : spellable(a)) worst = std::max(worst, std::fabs(diff(panvar::log_haploid_sequence_prior(a, s, ra), panvar::log_haploid_sequence_prior(b, s, rb))));
    return worst;
}

// A 3-block base fixture: polymorphic, polymorphic, polymorphic, with invariant sequence between.
panvar::StructuralLocusInput base_fixture() {
    panvar::StructuralLocusInput in;
    in.locus_id = "base";
    in.blocks = {{"b0", {{"a", "TTACGA"}, {"c", "TTCCGA"}}},       // shared prefix TT, shared suffix GA
                 {"b1", {{"g", "GGG"}, {"t", "GTG"}, {"d", "G"}}},   // includes a deletion-like short allele
                 {"b2", {{"x", "AAAC"}, {"y", "AAAT"}}}};
    in.invariant_segments = {"CC", "AAAAAAAAAAAA", "TTTTTTTT", "GG"};
    in.templates = {tpl("ref", {"a", "g", "x"}), tpl("h1", {"c", "t", "y"}), tpl("h2", {"a", "d", "y"}), tpl("h3", {"c", "g", "x"})};
    return in;
}

// G1: valid distances and transitions.
void g1_validity() {
    std::mt19937_64 rng(4040);
    std::size_t bad = 0, checked = 0;
    for (int trial = 0; trial < 200; ++trial) {
        auto in = base_fixture();
        in.invariant_segments[1] = std::string(rng() % 30, 'A');
        in.invariant_segments[2] = std::string(rng() % 30, 'T');
        const auto m = panvar::prepare_genotype_model(in);
        const auto g = reference_geometry(m, ref_ids(m, "ref"));
        const double lambda = 1.0 + static_cast<double>(rng() % 500);
        for (double r : hazard(g.distance, lambda)) {
            const std::size_t n = m.panel().template_count();
            for (std::size_t i = 0; i < n; ++i) {
                double sum = 0.0;
                for (std::size_t j = 0; j < n; ++j) {
                    const double t = (i == j ? 1.0 - r : 0.0) + r / static_cast<double>(n);
                    bad += t < 0.0 ? 1 : 0;
                    sum += t;
                }
                bad += std::fabs(sum - 1.0) > 1e-15 ? 1 : 0;
            }
            bad += (r < 0.0 || r >= 1.0) ? 1 : 0;
            ++checked;
        }
        for (double d : g.distance) bad += d < 0.0 ? 1 : 0;
    }
    row("G1_validity", "violations(d<0, r outside [0,1), negative entry, row sum != 1) over 400 boundaries", static_cast<double>(bad), "0", bad == 0 && checked == 400);
}

// G2: zero distance is the identity transition: only same-donor histories survive.
void g2_zero_distance() {
    panvar::StructuralLocusInput in;
    in.locus_id = "zero";
    in.blocks = {{"b0", {{"a", "A"}, {"c", "C"}}}, {"b1", {{"g", "G"}, {"t", "T"}}}};   // adjacent, both single-base cores
    in.invariant_segments = {"", "", ""};
    in.templates = {tpl("ref", {"a", "g"}), tpl("h1", {"c", "t"})};
    const auto m = panvar::prepare_genotype_model(in);
    // force zero distance explicitly (two zero-extent states at the same point), as an insertion pair would give
    const std::vector<double> r = hazard({0.0});
    row("G2_zero_distance", "r_at_d=0", r[0], "0 (T = I)", r[0] == 0.0);
    const double parental = panvar::log_haploid_sequence_prior(m, "AG", r), recomb = panvar::log_haploid_sequence_prior(m, "AT", r);
    row("G2_zero_distance", "parental_AG_minus_oracle", diff(parental, oracle_sequence(m, "AG", r)), "|.| <= 1e-13", std::fabs(diff(parental, oracle_sequence(m, "AG", r))) <= 1e-13);
    row("G2_zero_distance", "recombinant_AT_is_-inf", std::isinf(recomb) ? 1.0 : 0.0, "1 (no redraw at zero distance)", std::isinf(recomb) && recomb < 0);
    const auto g = reference_geometry(m, ref_ids(m, "ref"));
    row("G2_zero_distance", "adjacent single-base cores: d", g.distance[0], "1 (core midpoints 0.5 and 1.5: adjacency alone is not zero distance)", g.distance[0] == 1.0);
}

// G3: neutral subdivision -- a monomorphic block inserted inside an invariant segment, at every split point.
void g3_neutral_insertion() {
    const auto in = base_fixture();
    const auto m = panvar::prepare_genotype_model(in);
    const auto g = reference_geometry(m, ref_ids(m, "ref"));
    const std::string seg = in.invariant_segments[1];
    for (std::size_t cut : {std::size_t(0), std::size_t(5), seg.size() - 4}) {
        auto sub = in;
        sub.blocks.insert(sub.blocks.begin() + 1, {"m", {{"m", seg.substr(cut, 4)}}});
        sub.invariant_segments = {in.invariant_segments[0], seg.substr(0, cut), seg.substr(cut + 4), in.invariant_segments[2], in.invariant_segments[3]};
        for (auto& t : sub.templates) t.allele_source_ids.insert(t.allele_source_ids.begin() + 1, "m");
        const auto ms = panvar::prepare_genotype_model(sub);
        const auto gs = reference_geometry(ms, ref_ids(ms, "ref"));
        const double composed = gs.distance[0] + gs.distance[1];
        row("G3_neutral_insertion", "cut=" + std::to_string(cut) + ":d01+d12_minus_d01", composed - g.distance[0], "0 (distances add)", std::fabs(composed - g.distance[0]) <= 1e-12);
        same_priors("G3_neutral_insertion", "cut=" + std::to_string(cut) + "(distance-aware, lambda fixed)", m, hazard(g.distance), ms, hazard(gs.distance));
        const double legacy = max_prior_change(m, std::vector<double>(2, 0.01), ms, std::vector<double>(3, 0.01));
        row("G3_neutral_insertion", "cut=" + std::to_string(cut) + ":CONTROL_constant_r_changes_prior", legacy, "> 1e-6 (legacy is not subdivision-invariant)", legacy > 1e-6);
    }
}

// G4: shared sequence written inside a polymorphic block vs as invariant sequence vs as a monomorphic block.
void g4_shared_affix() {
    const auto in = base_fixture();
    const auto m = panvar::prepare_genotype_model(in);
    // (a) the shared prefix "TT" of b0 moved into the preceding invariant segment
    auto moved = in;
    moved.blocks[0].alleles = {{"a", "ACGA"}, {"c", "CCGA"}};
    moved.invariant_segments[0] = in.invariant_segments[0] + "TT";
    // (b) the same "TT" carved into its own monomorphic block
    auto carved = in;
    carved.blocks.insert(carved.blocks.begin(), {"m", {{"m", "TT"}}});
    carved.blocks[1].alleles = {{"a", "ACGA"}, {"c", "CCGA"}};
    carved.invariant_segments.insert(carved.invariant_segments.begin() + 1, "");
    for (auto& t : carved.templates) t.allele_source_ids.insert(t.allele_source_ids.begin(), "m");
    const auto mm = panvar::prepare_genotype_model(moved), mc = panvar::prepare_genotype_model(carved);
    const auto g = reference_geometry(m, ref_ids(m, "ref")), gm = reference_geometry(mm, ref_ids(mm, "ref")), gc = reference_geometry(mc, ref_ids(mc, "ref"));
    row("G4_shared_affix", "moved_into_invariant:max|d - d0|", std::max(std::fabs(gm.distance[0] - g.distance[0]), std::fabs(gm.distance[1] - g.distance[1])), "0", gm.distance == g.distance);
    same_priors("G4_shared_affix", "moved_into_invariant", m, hazard(g.distance), mm, hazard(gm.distance));
    row("G4_shared_affix", "carved_mono_block:d(m,b0)+d(b0,b1)_minus_d(b0,b1)... (b0 position unchanged)", gc.position[1] - g.position[0], "0", std::fabs(gc.position[1] - g.position[0]) <= 1e-12);
    std::vector<double> carved_r = hazard(gc.distance);
    same_priors("G4_shared_affix", "carved_mono_block", m, hazard(g.distance), mc, carved_r);
    // controls: untrimmed block midpoints and M3-A panel-mean both move with the affix
    const auto bm = reference_geometry(m, ref_ids(m, "ref"), nullptr, false), bmm = reference_geometry(mm, ref_ids(mm, "ref"), nullptr, false);
    row("G4_shared_affix", "CONTROL_untrimmed_block_midpoint_changes_d01", std::fabs(bmm.distance[0] - bm.distance[0]), "> 0 (block midpoints depend on where shared bytes are written)", std::fabs(bmm.distance[0] - bm.distance[0]) > 0);
    const auto pm = panel_mean(m), pmm = panel_mean(mm);
    row("G4_shared_affix", "CONTROL_panel_mean_changes_d01", std::fabs(pmm[0] - pm[0]), "> 0 (M3-A span depends on it too)", std::fabs(pmm[0] - pm[0]) > 0);
}

std::string rc(std::string s) {
    std::reverse(s.begin(), s.end());
    for (char& c : s) c = c == 'A' ? 'T' : c == 'C' ? 'G' : c == 'G' ? 'C' : c == 'T' ? 'A' : c;
    return s;
}

// G5: orientation -- the same locus written in the reverse-complement chain orientation.
void g5_orientation() {
    const auto in = base_fixture();
    panvar::StructuralLocusInput rev;
    rev.locus_id = "rev";
    rev.chain_orientation = panvar::ChainOrientation::Reverse;
    for (auto it = in.blocks.rbegin(); it != in.blocks.rend(); ++it) {
        panvar::GenotypeBlockInput b{it->block_id, {}, it->kind};
        for (const auto& a : it->alleles) b.alleles.push_back({a.source_id, rc(a.sequence)});
        rev.blocks.push_back(b);
    }
    for (auto it = in.invariant_segments.rbegin(); it != in.invariant_segments.rend(); ++it) rev.invariant_segments.push_back(rc(*it));
    for (const auto& t : in.templates) {
        auto ids = t.allele_source_ids;
        std::reverse(ids.begin(), ids.end());
        rev.templates.push_back(tpl(t.name, ids));
    }
    const auto m = panvar::prepare_genotype_model(in), mr = panvar::prepare_genotype_model(rev);
    const auto g = reference_geometry(m, ref_ids(m, "ref")), gr = reference_geometry(mr, ref_ids(mr, "ref"));
    auto back = gr.distance;
    std::reverse(back.begin(), back.end());
    double worst = 0.0;
    for (std::size_t i = 0; i < back.size(); ++i) worst = std::max(worst, std::fabs(back[i] - g.distance[i]));
    row("G5_orientation", "reverse_chain_distances_mirror_forward(max abs diff)", worst, "0", worst <= 1e-12);
    double pw = 0.0;
    for (const auto& s : spellable(m)) pw = std::max(pw, std::fabs(diff(panvar::log_haploid_sequence_prior(m, s, hazard(g.distance)), panvar::log_haploid_sequence_prior(mr, rc(s), hazard(gr.distance)))));
    row("G5_orientation", "reverse_chain_sequence_priors_equal(max abs diff)", pw, "<= 1e-12 (uniform start, symmetric T: reversible)", pw <= 1e-12);
}

// G6: polymorphic subdivision is NOT neutral.
void g6_polymorphic_split() {
    panvar::StructuralLocusInput whole, split;
    whole.locus_id = "whole";
    split.locus_id = "split";
    whole.blocks = {{"b0", {{"p", "AXXXXC"}, {"q", "GXXXXT"}}}, {"b1", {{"g", "G"}, {"t", "T"}}}};
    whole.invariant_segments = {"N", "NNNN", "N"};
    whole.templates = {tpl("ref", {"p", "g"}), tpl("h1", {"q", "t"})};
    split.blocks = {{"b0a", {{"a", "A"}, {"g", "G"}}}, {"b0b", {{"c", "C"}, {"t", "T"}}}, {"b1", {{"g", "G"}, {"t", "T"}}}};
    split.invariant_segments = {"N", "XXXX", "NNNN", "N"};
    split.templates = {tpl("ref", {"a", "c", "g"}), tpl("h1", {"g", "t", "t"})};
    const auto mw = panvar::prepare_genotype_model(whole), ms = panvar::prepare_genotype_model(split);
    const auto sw = spellable(mw), ss = spellable(ms);
    std::size_t new_admissible = 0;
    const auto gs = reference_geometry(ms, ref_ids(ms, "ref"));
    for (const auto& s : ss) if (!sw.count(s) && std::isfinite(panvar::log_haploid_sequence_prior(ms, s, hazard(gs.distance)))) ++new_admissible;
    row("G6_polymorphic_split", "new_admissible_recombinant_sequences", static_cast<double>(new_admissible), "> 0 => classified NON_NEUTRAL (never an invariance case)", new_admissible > 0);
    const auto gw = reference_geometry(mw, ref_ids(mw, "ref"));
    row("G6_polymorphic_split", "core_extent_whole_vs_split_positions", gw.position[0] - 0.5 * (gs.position[0] + gs.position[1]), "info: the split places two states inside the old core", true);
}

// G7: sequence-level summation over equivalent spellings, with distance-aware r.
void g7_equivalent_spellings() {
    panvar::StructuralLocusInput in;
    in.locus_id = "spell";
    in.blocks = {{"b0", {{"a", "A"}, {"ac", "AC"}}}, {"b1", {{"cg", "CG"}, {"g", "G"}}}, {"b2", {{"t", "T"}, {"a", "A"}}}};
    in.invariant_segments = {"", "", "TTTT", ""};
    in.templates = {tpl("ref", {"a", "cg", "t"}), tpl("h1", {"ac", "g", "a"}), tpl("h2", {"a", "cg", "a"}), tpl("h3", {"ac", "g", "t"})};
    const auto m = panvar::prepare_genotype_model(in);
    const auto g = reference_geometry(m, ref_ids(m, "ref"));
    const auto r = hazard(g.distance);
    std::size_t multi = 0;
    double worst = 0.0;
    for (const auto& s : spellable(m)) {
        std::size_t spellings = 0;
        for (const auto& p : all_paths(m)) spellings += panvar::spell_mosaic(m, panvar::MosaicPath{p}).bases == s ? 1 : 0;
        multi += spellings > 1 ? 1 : 0;
        worst = std::max(worst, std::fabs(diff(panvar::log_haploid_sequence_prior(m, s, r), oracle_sequence(m, s, r))));
    }
    row("G7_equivalent_spellings", "sequences_with_more_than_one_spelling", static_cast<double>(multi), ">= 1 (the fixture exercises aggregation)", multi >= 1);
    row("G7_equivalent_spellings", "library_sequence_prior_minus_summed_oracle(max abs)", worst, "<= 1e-12", worst <= 1e-12);
}

// G8: donor multiplicity and redraw semantics stay explicit.
void g8_donors() {
    auto in = base_fixture();
    const auto m = panvar::prepare_genotype_model(in);
    const auto g = reference_geometry(m, ref_ids(m, "ref"));
    const auto r = hazard(g.distance);
    auto dup = in;
    dup.templates.push_back(tpl("h1_dup", in.templates[1].allele_source_ids));
    const auto md = panvar::prepare_genotype_model(dup);
    const auto gd = reference_geometry(md, ref_ids(md, "ref"));
    row("G8_donors", "duplicated_row_leaves_geometry_unchanged", static_cast<double>(gd.distance == g.distance), "1", gd.distance == g.distance);
    double worst = 0.0, moved = 0.0;
    for (const auto& s : spellable(md)) {
        worst = std::max(worst, std::fabs(diff(panvar::log_haploid_sequence_prior(md, s, r), oracle_sequence(md, s, r))));
        moved = std::max(moved, std::fabs(diff(panvar::log_haploid_sequence_prior(md, s, r), panvar::log_haploid_sequence_prior(m, s, r))));
    }
    row("G8_donors", "duplicated_row:library_minus_oracle(U counts rows)", worst, "<= 1e-12", worst <= 1e-12);
    row("G8_donors", "duplicated_row_reweights_prior(max abs change)", moved, "> 0 (multiplicity is a declared property of U, not hidden)", moved > 0);
    auto perm = in;
    std::reverse(perm.templates.begin(), perm.templates.end());
    const auto mp = panvar::prepare_genotype_model(perm);
    same_priors("G8_donors", "template_order", m, r, mp, hazard(reference_geometry(mp, ref_ids(mp, "ref")).distance));
    const double n = static_cast<double>(m.panel().template_count());
    row("G8_donors", "actual_donor_change_prob_b0b1 = r(n-1)/n", r[0] * (n - 1) / n, "info (redraw includes the current donor)", true);
}

// G9: leaving out a donor. Cores trimmed against the FRAME are unchanged; cores trimmed against the callable alleles move
// when the removed donor carried the only divergent allele; M3-A's panel-mean moves whenever a removed row's allele length
// differs from the mean.
void g9_donor_removal() {
    const auto in = base_fixture();
    const auto m = panvar::prepare_genotype_model(in);
    std::vector<std::vector<std::string>> frame;
    for (const auto& b : in.blocks) { std::vector<std::string> v; for (const auto& a : b.alleles) v.push_back(a.sequence); frame.push_back(v); }
    auto loo = in;   // remove h2, the only carrier of the short b1 allele "G"
    loo.templates.erase(loo.templates.begin() + 2);
    loo.blocks[1].alleles.pop_back();
    const auto ml = panvar::prepare_genotype_model(loo);
    const auto gf = reference_geometry(m, ref_ids(m, "ref"), &frame), glf = reference_geometry(ml, ref_ids(ml, "ref"), &frame);
    row("G9_donor_removal", "frame_cores:max|d_loo - d|", std::max(std::fabs(glf.distance[0] - gf.distance[0]), std::fabs(glf.distance[1] - gf.distance[1])), "0 (geometry fixed by the frame)", glf.distance == gf.distance);
    const auto gc = reference_geometry(m, ref_ids(m, "ref")), glc = reference_geometry(ml, ref_ids(ml, "ref"));
    row("G9_donor_removal", "CONTROL_callable_cores:max|d_loo - d|", std::max(std::fabs(glc.distance[0] - gc.distance[0]), std::fabs(glc.distance[1] - gc.distance[1])), "> 0 (callable-allele cores are panel-dependent)", glc.distance != gc.distance);
    const auto pm = panel_mean(m), pml = panel_mean(ml);
    row("G9_donor_removal", "CONTROL_panel_mean:max|d_loo - d|", std::max(std::fabs(pml[0] - pm[0]), std::fabs(pml[1] - pm[1])), "> 0 (M3-A span is panel-dependent)", pml != pm);
}

// G10: the scale is a length, not a per-boundary count: refining the representation does not change total hazard, while
// the hazard-matched control (lambda re-derived to keep (B-1) * -log(0.99)) does -- which is why it is a control only.
void g10_scale() {
    const auto in = base_fixture();
    const auto m = panvar::prepare_genotype_model(in);
    const auto g = reference_geometry(m, ref_ids(m, "ref"));
    auto sub = in;
    sub.blocks.insert(sub.blocks.begin() + 2, {"m", {{"m", "TTTT"}}});
    sub.invariant_segments = {in.invariant_segments[0], in.invariant_segments[1], "", "TTTT", in.invariant_segments[3]};
    for (auto& t : sub.templates) t.allele_source_ids.insert(t.allele_source_ids.begin() + 2, "m");
    const auto ms = panvar::prepare_genotype_model(sub);
    const auto gs = reference_geometry(ms, ref_ids(ms, "ref"));
    double h0 = 0, h1 = 0;
    for (double d : g.distance) h0 += d / kLambda;
    for (double d : gs.distance) h1 += d / kLambda;
    row("G10_scale", "fixed_lambda_total_hazard_change_after_neutral_subdivision", h1 - h0, "0", std::fabs(h1 - h0) <= 1e-12);
    const double matched0 = std::accumulate(g.distance.begin(), g.distance.end(), 0.0) / (static_cast<double>(g.distance.size()) * -std::log1p(-0.01));
    const double matched1 = std::accumulate(gs.distance.begin(), gs.distance.end(), 0.0) / (static_cast<double>(gs.distance.size()) * -std::log1p(-0.01));
    row("G10_scale", "CONTROL_hazard_matched_lambda_changes_with_block_count", matched1 - matched0, "!= 0 (block-count dependent: diagnostic control only)", std::fabs(matched1 - matched0) > 0);
}

}  // namespace

int main() {
    try {
        test_li_stephens_prior();   // the oracle's own regression fixtures first
        row("G0_oracle", "mosaic_core brute-force prior fixtures", 1.0, "pass", true);
        g1_validity();
        g2_zero_distance();
        g3_neutral_insertion();
        g4_shared_affix();
        g5_orientation();
        g6_polymorphic_split();
        g7_equivalent_spellings();
        g8_donors();
        g9_donor_removal();
        g10_scale();
    } catch (const std::exception& e) {
        std::cerr << "m4a distance prior: ERROR " << e.what() << '\n';
        return 2;
    }
    std::cout << "test\titem\tvalue\texpectation\tresult\n" << g_rows.str() << "m4a_distance_prior\t" << (g_failed ? "FAIL" : "PASS") << '\n';
    return g_failed ? 1 : 0;
}
