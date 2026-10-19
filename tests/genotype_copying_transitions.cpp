// M4-B: gates for the SHARED copying-transition implementation -- the adapter's geometry rule (compute_locus_geometry,
// the function build_locus_input calls) and the one builder (build_copying_transitions) the caller and score-set use --
// against the brute-force history oracle of tests/genotype_mosaic_core.cpp (included unchanged, as M3-A and M4-A did).
// The M4-A invariance gates G1-G10 are re-run here against the shared code (not the M4-A prototype), plus N1-N6:
// parameter validation and explicit refusal, constant-mode bit identity, saturation, zero-distance exclusion, digests
// and the hazard-matched control. Geometry, lambda and U are frozen across every invariance comparison.
#define main genotype_mosaic_core_main
#include "genotype_mosaic_core.cpp"
#undef main

#include "panvar/copying_transitions.hpp"
#include "panvar/genotype_locus.hpp"

#include <iomanip>
#include <random>
#include <set>
#include <sstream>

namespace {

std::size_t g_failed = 0;
std::ostringstream g_rows;
constexpr double kLambda = 40.0;

void row(const std::string& test, const std::string& item, double value, const std::string& expectation, bool pass) {
    g_rows << test << '\t' << item << '\t' << std::setprecision(15) << value << '\t' << expectation << '\t' << (pass ? "PASS" : "FAIL") << '\n';
    g_failed += pass ? 0 : 1;
}

panvar::PanelTemplateInput tpl(const std::string& n, std::vector<std::string> ids) {
    return {n, panvar::TemplateFrameStatus::Complete, std::move(ids), panvar::ChainOrientation::Forward, ""};
}

// The shared geometry rule on a fixture: the frame = every allele the input lists (pre-curation), the reference = the
// named template's alleles.
panvar::LocusGeometry geometry(const panvar::StructuralLocusInput& in, const std::string& ref,
                               const std::vector<std::vector<std::string>>* frame_override = nullptr) {
    std::vector<std::string> ids;
    std::vector<std::vector<std::string>> frame;
    for (const auto& b : in.blocks) {
        ids.push_back(b.block_id);
        std::vector<std::string> v;
        for (const auto& a : b.alleles) v.push_back(a.sequence);
        frame.push_back(v);
    }
    std::vector<std::string> ref_alleles;
    for (const auto& t : in.templates) {
        if (t.name != ref) continue;
        for (std::size_t b = 0; b < in.blocks.size(); ++b) {
            for (const auto& a : in.blocks[b].alleles) if (a.source_id == t.allele_source_ids[b]) ref_alleles.push_back(a.sequence);
        }
    }
    return panvar::compute_locus_geometry(ids, in.invariant_segments, frame_override ? *frame_override : frame, &ref_alleles, ref);
}

panvar::CopyingTransitions transitions(const panvar::StructuralLocusInput& in, const std::string& ref = "ref", double lambda = kLambda,
                                       const std::vector<std::vector<std::string>>* frame = nullptr) {
    const panvar::LocusGeometry g = geometry(in, ref, frame);
    panvar::CopyingTransitionSpec spec;
    spec.mode = "distance_ref_core_v1_experimental";
    spec.copying_length_bp = lambda;
    return panvar::build_copying_transitions(in.blocks.size(), &g, spec);
}

std::vector<std::vector<std::size_t>> all_paths(const panvar::PreparedGenotypeModel& m) {
    std::vector<std::vector<std::size_t>> paths{{}};
    for (const auto& block : m.blocks()) {
        std::vector<std::vector<std::size_t>> next;
        for (const auto& p : paths) for (std::size_t a = 0; a < block.alleles.size(); ++a) { auto q = p; q.push_back(a); next.push_back(q); }
        paths = next;
    }
    return paths;
}
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
double diff(double a, double b) { return (std::isinf(a) && std::isinf(b) && (a < 0) == (b < 0)) ? 0.0 : a - b; }
double max_change(const panvar::PreparedGenotypeModel& a, const std::vector<double>& ra, const panvar::PreparedGenotypeModel& b,
                  const std::vector<double>& rb) {
    double worst = 0.0;
    for (const auto& s : spellable(a)) worst = std::max(worst, std::fabs(diff(panvar::log_haploid_sequence_prior(a, s, ra), panvar::log_haploid_sequence_prior(b, s, rb))));
    return worst;
}
void same_priors(const std::string& test, const std::string& what, const panvar::PreparedGenotypeModel& a, const std::vector<double>& ra,
                 const panvar::PreparedGenotypeModel& b, const std::vector<double>& rb) {
    row(test, what + ":same_spellable_set", static_cast<double>(spellable(a) == spellable(b)), "1", spellable(a) == spellable(b));
    const double w = max_change(a, ra, b, rb);
    row(test, what + ":max_abs_prior_difference", w, "<= 1e-12", w <= 1e-12);
}
double oracle_gap(const panvar::PreparedGenotypeModel& m, const std::vector<double>& r) {
    double worst = 0.0;
    for (const auto& s : spellable(m)) worst = std::max(worst, std::fabs(diff(panvar::log_haploid_sequence_prior(m, s, r), oracle_sequence(m, s, r))));
    return worst;
}

panvar::StructuralLocusInput base_fixture() {
    panvar::StructuralLocusInput in;
    in.locus_id = "base";
    in.blocks = {{"b0", {{"a", "TTACGA"}, {"c", "TTCCGA"}}}, {"b1", {{"g", "GGG"}, {"t", "GTG"}, {"d", "G"}}}, {"b2", {{"x", "AAAC"}, {"y", "AAAT"}}}};
    in.invariant_segments = {"CC", "AAAAAAAAAAAA", "TTTTTTTT", "GG"};
    in.templates = {tpl("ref", {"a", "g", "x"}), tpl("h1", {"c", "t", "y"}), tpl("h2", {"a", "d", "y"}), tpl("h3", {"c", "g", "x"})};
    return in;
}
std::string rc(std::string s) {
    std::reverse(s.begin(), s.end());
    for (char& c : s) c = c == 'A' ? 'T' : c == 'C' ? 'G' : c == 'G' ? 'C' : c == 'T' ? 'A' : c;
    return s;
}

void g1_validity() {
    std::mt19937_64 rng(4040);
    std::size_t bad = 0, checked = 0;
    for (int trial = 0; trial < 200; ++trial) {
        auto in = base_fixture();
        in.invariant_segments[1] = std::string(rng() % 30, 'A');
        in.invariant_segments[2] = std::string(rng() % 30, 'T');
        const auto t = transitions(in, "ref", 1.0 + static_cast<double>(rng() % 500));
        for (std::size_t b = 0; b < t.r.size(); ++b) {
            bad += (t.distance_bp[b] < 0.0 || t.r[b] < 0.0 || t.r[b] > 1.0) ? 1 : 0;
            const double n = 4.0;
            bad += std::fabs(((1.0 - t.r[b]) + t.r[b] / n) + 3.0 * t.r[b] / n - 1.0) > 1e-15 ? 1 : 0;   // a T row sums to 1
            ++checked;
        }
    }
    row("G1_validity", "violations over 400 boundaries (shared builder)", static_cast<double>(bad), "0", bad == 0 && checked == 400);
}

void g2_zero_distance_and_n4_exclusion() {
    // Two informative states at ONE reference position: adjacent zero-reference-span insertion blocks.
    panvar::StructuralLocusInput in;
    in.locus_id = "zero";
    in.blocks = {{"i0", {{"e0", ""}, {"c", "C"}}}, {"i1", {{"e1", ""}, {"g", "G"}}}};
    in.invariant_segments = {"AAAA", "", "TTTT"};
    in.templates = {tpl("ref", {"e0", "e1"}), tpl("h1", {"c", "g"}), tpl("h2", {"c", "e1"})};
    const auto t = transitions(in);
    row("G2_zero_distance", "d(i0,i1) for two insertions at one reference point", t.distance_bp[0], "0", t.distance_bp[0] == 0.0);
    row("G2_zero_distance", "r at d = 0", t.r[0], "0 (T = I)", t.r[0] == 0.0);
    row("N4_zero_distance_inventory", "zero_distance_boundaries / informative", static_cast<double>(t.zero_distance_boundaries * 10 + t.zero_distance_informative_boundaries),
        "11 (one boundary, both states informative)", t.zero_distance_boundaries == 1 && t.zero_distance_informative_boundaries == 1);
    const auto m = panvar::prepare_genotype_model(in);
    const double recomb = panvar::log_haploid_sequence_prior(m, "AAAAGTTTT", t.r);   // allele "" at i0 with "G" at i1: no template carries it
    const double recomb_const = panvar::log_haploid_sequence_prior(m, "AAAAGTTTT", {0.01});
    row("N4_zero_distance_exclusion", "recombinant (empty, G) under distance mode", std::isinf(recomb) && recomb < 0 ? 1.0 : 0.0,
        "1 (-inf: excluded at zero distance; no probability floor)", std::isinf(recomb) && recomb < 0);
    row("N4_zero_distance_exclusion", "same recombinant under constant r = 0.01 is finite", std::isfinite(recomb_const) ? 1.0 : 0.0, "1", std::isfinite(recomb_const));
    row("G2_zero_distance", "library_minus_oracle (all sequences)", oracle_gap(m, t.r), "<= 1e-13", oracle_gap(m, t.r) <= 1e-13);
}

void g3_neutral_insertion() {
    const auto in = base_fixture();
    const auto m = panvar::prepare_genotype_model(in);
    const auto t = transitions(in);
    const std::string seg = in.invariant_segments[1];
    for (std::size_t cut : {std::size_t(0), std::size_t(5), seg.size() - 4}) {
        auto sub = in;
        sub.blocks.insert(sub.blocks.begin() + 1, {"m", {{"m", seg.substr(cut, 4)}}});
        sub.invariant_segments = {in.invariant_segments[0], seg.substr(0, cut), seg.substr(cut + 4), in.invariant_segments[2], in.invariant_segments[3]};
        for (auto& tp : sub.templates) tp.allele_source_ids.insert(tp.allele_source_ids.begin() + 1, "m");
        const auto ms = panvar::prepare_genotype_model(sub);
        const auto ts = transitions(sub);
        row("G3_neutral_insertion", "cut=" + std::to_string(cut) + ":distances_add", ts.distance_bp[0] + ts.distance_bp[1] - t.distance_bp[0], "0",
            ts.distance_bp[0] + ts.distance_bp[1] == t.distance_bp[0]);
        same_priors("G3_neutral_insertion", "cut=" + std::to_string(cut), m, t.r, ms, ts.r);
        const double legacy = max_change(m, std::vector<double>(2, 0.01), ms, std::vector<double>(3, 0.01));
        row("G3_neutral_insertion", "cut=" + std::to_string(cut) + ":CONTROL_constant_r_changes_prior", legacy, "> 1e-6", legacy > 1e-6);
    }
}

void g4_shared_affix() {
    const auto in = base_fixture();
    const auto m = panvar::prepare_genotype_model(in);
    auto moved = in;
    moved.blocks[0].alleles = {{"a", "ACGA"}, {"c", "CCGA"}};
    moved.invariant_segments[0] = in.invariant_segments[0] + "TT";
    auto carved = in;
    carved.blocks.insert(carved.blocks.begin(), {"m", {{"m", "TT"}}});
    carved.blocks[1].alleles = {{"a", "ACGA"}, {"c", "CCGA"}};
    carved.invariant_segments.insert(carved.invariant_segments.begin() + 1, "");
    for (auto& tp : carved.templates) tp.allele_source_ids.insert(tp.allele_source_ids.begin(), "m");
    const auto mm = panvar::prepare_genotype_model(moved), mc = panvar::prepare_genotype_model(carved);
    const auto t = transitions(in), tm = transitions(moved), tc = transitions(carved);
    row("G4_shared_affix", "moved_into_invariant:distances_equal", static_cast<double>(tm.distance_bp == t.distance_bp), "1", tm.distance_bp == t.distance_bp);
    same_priors("G4_shared_affix", "moved_into_invariant", m, t.r, mm, tm.r);
    // Amendment 1: the carved neutral block precedes b0, so the ORIGINAL blocks' distances are unchanged (the M4-B freeze
    // asserted tc[1] == t[0] - tc[0], a test defect); the neutral state's own boundary does not enter any prior
    row("G4_shared_affix", "A1:carved_mono_block:original_distances_unchanged", (tc.distance_bp[1] - t.distance_bp[0]) + (tc.distance_bp[2] - t.distance_bp[1]),
        "0", tc.distance_bp[1] == t.distance_bp[0] && tc.distance_bp[2] == t.distance_bp[1]);
    same_priors("G4_shared_affix", "carved_mono_block", m, t.r, mc, tc.r);
}

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
    for (const auto& tp : in.templates) {
        auto ids = tp.allele_source_ids;
        std::reverse(ids.begin(), ids.end());
        rev.templates.push_back(tpl(tp.name, ids));
    }
    const auto m = panvar::prepare_genotype_model(in), mr = panvar::prepare_genotype_model(rev);
    const auto t = transitions(in), tr = transitions(rev);
    auto back = tr.distance_bp;
    std::reverse(back.begin(), back.end());
    row("G5_orientation", "reverse_chain_distances_mirror", static_cast<double>(back == t.distance_bp), "1", back == t.distance_bp);
    double w = 0.0;
    for (const auto& s : spellable(m)) w = std::max(w, std::fabs(diff(panvar::log_haploid_sequence_prior(m, s, t.r), panvar::log_haploid_sequence_prior(mr, rc(s), tr.r))));
    row("G5_orientation", "reverse_chain_priors_equal", w, "<= 1e-12", w <= 1e-12);
}

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
    const auto ts = transitions(split);
    std::size_t fresh = 0;
    const auto sw = spellable(mw);
    for (const auto& s : spellable(ms)) if (!sw.count(s) && std::isfinite(panvar::log_haploid_sequence_prior(ms, s, ts.r))) ++fresh;
    row("G6_polymorphic_split", "new_admissible_recombinants", static_cast<double>(fresh), "> 0 => NON_NEUTRAL", fresh > 0);
}

void g7_equivalent_spellings() {
    panvar::StructuralLocusInput in;
    in.locus_id = "spell";
    in.blocks = {{"b0", {{"a", "A"}, {"ac", "AC"}}}, {"b1", {{"cg", "CG"}, {"g", "G"}}}, {"b2", {{"t", "T"}, {"a", "A"}}}};
    in.invariant_segments = {"", "", "TTTT", ""};
    in.templates = {tpl("ref", {"a", "cg", "t"}), tpl("h1", {"ac", "g", "a"}), tpl("h2", {"a", "cg", "a"}), tpl("h3", {"ac", "g", "t"})};
    const auto m = panvar::prepare_genotype_model(in);
    const auto t = transitions(in);
    row("G7_equivalent_spellings", "library_sequence_prior_minus_summed_oracle", oracle_gap(m, t.r), "<= 1e-12", oracle_gap(m, t.r) <= 1e-12);
}

void g8_donors() {
    auto in = base_fixture();
    const auto m = panvar::prepare_genotype_model(in);
    const auto t = transitions(in);
    auto dup = in;
    dup.templates.push_back(tpl("h1_dup", in.templates[1].allele_source_ids));
    const auto md = panvar::prepare_genotype_model(dup);
    const auto td = transitions(dup);
    row("G8_donors", "duplicated_row_leaves_geometry_unchanged", static_cast<double>(td.geometry_digest == t.geometry_digest), "1", td.geometry_digest == t.geometry_digest);
    row("G8_donors", "duplicated_row:library_minus_oracle (U counts rows)", oracle_gap(md, td.r), "<= 1e-12", oracle_gap(md, td.r) <= 1e-12);
    row("G8_donors", "duplicated_row_reweights_prior", max_change(md, td.r, m, t.r), "> 0", max_change(md, td.r, m, t.r) > 0);
    auto perm = in;
    std::reverse(perm.templates.begin(), perm.templates.end());
    same_priors("G8_donors", "template_order", m, t.r, panvar::prepare_genotype_model(perm), transitions(perm).r);
}

void g9_curation_and_donor_removal() {
    // the frame is the pre-curation allele set: removing a donor (and its allele) from the CALLABLE model does not change
    // the geometry when the geometry is computed against the frame, as the adapter does
    auto with_e = base_fixture();
    with_e.blocks[0].alleles.push_back({"e", "TTACGT"});
    with_e.templates[2].allele_source_ids[0] = "e";
    std::vector<std::vector<std::string>> frame;
    for (const auto& b : with_e.blocks) { std::vector<std::string> v; for (const auto& a : b.alleles) v.push_back(a.sequence); frame.push_back(v); }
    auto curated = with_e;
    curated.templates.erase(curated.templates.begin() + 2);
    curated.blocks[0].alleles.pop_back();
    curated.blocks[1].alleles.pop_back();
    const auto t = transitions(with_e), tc_frame = transitions(curated, "ref", kLambda, &frame), tc_callable = transitions(curated);
    row("G9_curation", "frame_geometry_digest_unchanged_by_curation", static_cast<double>(tc_frame.geometry_digest == t.geometry_digest), "1",
        tc_frame.geometry_digest == t.geometry_digest);
    row("G9_curation", "CONTROL_callable_allele_geometry_moves", static_cast<double>(tc_callable.distance_bp != t.distance_bp), "1",
        tc_callable.distance_bp != t.distance_bp);
}

void g10_scale() {
    const auto in = base_fixture();
    auto sub = in;
    sub.blocks.insert(sub.blocks.begin() + 2, {"m", {{"m", "TTTT"}}});
    sub.invariant_segments = {in.invariant_segments[0], in.invariant_segments[1], "", "TTTT", in.invariant_segments[3]};
    for (auto& tp : sub.templates) tp.allele_source_ids.insert(tp.allele_source_ids.begin() + 2, "m");
    const auto t = transitions(in), ts = transitions(sub);
    double h0 = 0, h1 = 0;
    for (double d : t.distance_bp) h0 += d;
    for (double d : ts.distance_bp) h1 += d;
    row("G10_scale", "total_distance_unchanged_by_neutral_subdivision", h1 - h0, "0", h1 == h0);
    const double l0 = panvar::hazard_matched_copying_length(geometry(in, "ref"), 0.01), l1 = panvar::hazard_matched_copying_length(geometry(sub, "ref"), 0.01);
    row("G10_scale", "CONTROL_hazard_matched_length_depends_on_block_count", l1 - l0, "!= 0", l1 != l0);
    const double expect = h0 / (static_cast<double>(t.distance_bp.size()) * -std::log1p(-0.01));
    row("N6_hazard_matched", "formula: sum d / ((B-1) * -log(1-r))", l0 - expect, "0", l0 == expect);
}

void n1_validation() {
    const auto in = base_fixture();
    const panvar::LocusGeometry g = geometry(in, "ref");
    std::size_t refused = 0, cases = 0;
    const auto expect_refusal = [&](const std::string& what, auto&& f) {
        ++cases;
        try { f(); row("N1_validation", what, 0.0, "refused", false); }
        catch (const panvar::UnsupportedCopyingTransition&) { ++refused; row("N1_validation", what, 1.0, "refused (UnsupportedCopyingTransition)", true); }
    };
    panvar::CopyingTransitionSpec d;
    d.mode = "distance_ref_core_v1_experimental";
    expect_refusal("distance: lambda 0", [&] { d.copying_length_bp = 0.0; panvar::build_copying_transitions(3, &g, d); });
    expect_refusal("distance: lambda NaN", [&] { d.copying_length_bp = std::nan(""); panvar::build_copying_transitions(3, &g, d); });
    expect_refusal("distance: lambda inf", [&] { d.copying_length_bp = INFINITY; panvar::build_copying_transitions(3, &g, d); });
    d.copying_length_bp = 10.0;
    expect_refusal("distance: no geometry", [&] { panvar::build_copying_transitions(3, nullptr, d); });
    expect_refusal("distance: block count mismatch", [&] { panvar::build_copying_transitions(4, &g, d); });
    const panvar::LocusGeometry missing = panvar::compute_locus_geometry(g.block_ids, in.invariant_segments, {{"A"}, {"C"}, {"G"}}, nullptr, "ref",
                                                                         "the reference path's frame is partial");
    row("N1_validation", "unavailable geometry status is explicit", missing.available() ? 0.0 : 1.0, missing.status, !missing.available() &&
        missing.status == "unavailable: the reference path's frame is partial");
    expect_refusal("distance: unavailable geometry", [&] { panvar::build_copying_transitions(3, &missing, d); });
    panvar::CopyingTransitionSpec c;
    expect_refusal("constant: r = 0", [&] { c.constant_r = 0.0; panvar::build_copying_transitions(3, nullptr, c); });
    expect_refusal("constant: r = 0.6", [&] { c.constant_r = 0.6; panvar::build_copying_transitions(3, nullptr, c); });
    panvar::CopyingTransitionSpec u;
    u.mode = "something_else";
    expect_refusal("unknown mode", [&] { panvar::build_copying_transitions(3, nullptr, u); });
    row("N1_validation", "all refusals explicit", static_cast<double>(refused), std::to_string(cases), refused == cases);
}

void n2_constant_identity_and_n5_digests() {
    panvar::CopyingTransitionSpec c;
    c.constant_r = 0.01;
    const auto t = panvar::build_copying_transitions(11, nullptr, c);
    const std::vector<double> legacy(10, 0.01);
    row("N2_constant_identity", "constant_v1 vector equals the legacy vector bit for bit", static_cast<double>(t.r == legacy), "1", t.r == legacy);
    row("N2_constant_identity", "constant_v1 needs no geometry", static_cast<double>(t.geometry_digest == "none"), "1", t.geometry_digest == "none");
    const auto in = base_fixture();
    const auto a = transitions(in, "ref", 40.0), b = transitions(in, "ref", 40.0), x = transitions(in, "ref", 41.0);
    row("N5_digest", "same spec, same digest", static_cast<double>(a.vector_digest == b.vector_digest), "1", a.vector_digest == b.vector_digest);
    row("N5_digest", "different lambda, different digest", static_cast<double>(a.vector_digest != x.vector_digest), "1", a.vector_digest != x.vector_digest);
    row("N5_digest", "constant and distance digests differ", static_cast<double>(a.vector_digest != t.vector_digest), "1", a.vector_digest != t.vector_digest);
}

void n3_saturation() {
    const auto in = base_fixture();
    const auto m = panvar::prepare_genotype_model(in);
    const auto t = transitions(in, "ref", 1e-6);   // d >> lambda: r saturates to 1 (T -> U)
    bool all_one = true;
    for (double r : t.r) all_one = all_one && r == 1.0;
    row("N3_saturation", "r == 1 exactly when d >> lambda", all_one ? 1.0 : 0.0, "1 (T = U)", all_one);
    row("N3_saturation", "saturated priors equal the oracle", oracle_gap(m, t.r), "<= 1e-12", oracle_gap(m, t.r) <= 1e-12);
    bool finite = true;
    for (const auto& s : spellable(m)) finite = finite && std::isfinite(panvar::log_haploid_sequence_prior(m, s, t.r));
    row("N3_saturation", "every spellable sequence keeps a finite prior", finite ? 1.0 : 0.0, "1", finite);
    const auto tiny = transitions(in, "ref", 1e12);   // d << lambda: r ~ d / lambda, computed by expm1 without cancellation
    const double expect = t.distance_bp[0] / 1e12;
    row("N3_small_distance", "r = -expm1(-d/lambda) relative error at d/lambda ~ 1e-11", std::fabs(tiny.r[0] - expect) / expect, "<= 1e-9",
        std::fabs(tiny.r[0] - expect) / expect <= 1e-9);
}

}  // namespace

int main() {
    try {
        test_li_stephens_prior();
        row("G0_oracle", "mosaic_core brute-force prior fixtures", 1.0, "pass", true);
        g1_validity();
        g2_zero_distance_and_n4_exclusion();
        g3_neutral_insertion();
        g4_shared_affix();
        g5_orientation();
        g6_polymorphic_split();
        g7_equivalent_spellings();
        g8_donors();
        g9_curation_and_donor_removal();
        g10_scale();
        n1_validation();
        n2_constant_identity_and_n5_digests();
        n3_saturation();
    } catch (const std::exception& e) {
        std::cerr << "genotype copying transitions: ERROR " << e.what() << '\n';
        return 2;
    }
    std::cout << "test\titem\tvalue\texpectation\tresult\n" << g_rows.str() << "genotype_copying_transitions\t" << (g_failed ? "FAIL" : "PASS") << '\n';
    return g_failed ? 1 : 0;
}
