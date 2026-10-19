// M3-A: representation tests of the copying prior (diagnostic; not a production test).
//
// The independent oracles are REUSED, not rewritten: tests/genotype_mosaic_core.cpp is included unchanged (its main is
// renamed), which brings brute_haploid_prior / brute_best_haploid_history (brute force over every template history) and
// its existing fixtures. Every number here is printed as a row of a TSV; a check that fails makes the exit status 1.
//
// Prior models compared: the production constant-r Li-Stephens prior (T = (1-r)I + (r/n)11^T at every block boundary),
// and a distance model r_b = 1 - exp(-d_b / L) with d_b = panel-mean span (bp) from the start of block b-1 to the start
// of block b and a FIXED scale L (the same in every subdivision test).
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

void row(const std::string& test, const std::string& item, double value, const std::string& expectation, bool pass) {
    g_rows << test << '\t' << item << '\t' << std::setprecision(15) << value << '\t' << expectation << '\t' << (pass ? "PASS" : "FAIL") << '\n';
    g_failed += pass ? 0 : 1;
}

void note(const std::string& test, const std::string& item, double value, const std::string& meaning) {
    g_rows << test << '\t' << item << '\t' << std::setprecision(15) << value << '\t' << meaning << "\tINFO\n";
}

std::vector<double> distance_switches(const panvar::PreparedGenotypeModel& model, double scale) {
    const auto& panel = model.panel();
    std::vector<double> r;
    for (std::size_t b = 1; b < model.blocks().size(); ++b) {
        double mean = 0.0;
        for (const auto& t : panel.template_alleles) mean += static_cast<double>(model.blocks()[b - 1].alleles[t[b - 1]].sequence.size());
        const double d = mean / static_cast<double>(panel.template_count()) + static_cast<double>(model.invariant_segments()[b].size());
        r.push_back(d > 0.0 ? -std::expm1(-d / scale) : 0.0);
    }
    return r;
}

// Every sequence some allele path spells (brute force over paths), with its spelled bytes.
std::vector<std::string> spellable(const panvar::PreparedGenotypeModel& model) {
    std::vector<std::vector<std::size_t>> paths{{}};
    for (const auto& block : model.blocks()) {
        std::vector<std::vector<std::size_t>> next;
        for (const auto& p : paths) for (std::size_t a = 0; a < block.alleles.size(); ++a) { auto q = p; q.push_back(a); next.push_back(q); }
        paths = next;
    }
    std::set<std::string> out;
    for (const auto& p : paths) out.insert(panvar::spell_mosaic(model, panvar::MosaicPath{p}).bases);
    return {out.begin(), out.end()};
}

// R1: neutral subdivision (a monomorphic block inserted between two blocks) composes transitions.
void r1_neutral_subdivision() {
    panvar::PanelAlleleMatrix two, three;
    two.template_names = three.template_names = {"h0", "h1", "h2", "h3"};
    two.template_alleles = {{0, 0}, {1, 1}, {0, 1}, {1, 0}};
    three.template_alleles = {{0, 5, 0}, {1, 5, 1}, {0, 5, 1}, {1, 5, 0}};
    two.template_alleles[3] = {1, 1};   // keep an unbalanced panel so stay/redraw differ visibly
    three.template_alleles[3] = {1, 5, 1};
    for (const auto& rr : std::vector<std::pair<double, double>>{{0.01, 0.01}, {0.001, 0.05}, {0.2, 0.3}}) {
        const double composed = 1.0 - (1.0 - rr.first) * (1.0 - rr.second);
        for (const auto& p : {panvar::MosaicPath{{0, 1}}, panvar::MosaicPath{{0, 0}}, panvar::MosaicPath{{1, 1}}}) {
            const panvar::MosaicPath p3{{p.alleles[0], 5, p.alleles[1]}};
            const double sub = panvar::log_haploid_mosaic_prior(three, p3, {rr.first, rr.second});
            const double und = panvar::log_haploid_mosaic_prior(two, p, {composed});
            const double brute = std::log(static_cast<double>(brute_haploid_prior(three, p3, {rr.first, rr.second})));
            std::ostringstream item;
            item << "r1=" << rr.first << ",r2=" << rr.second << ",path=" << p.alleles[0] << p.alleles[1];
            row("R1_composition", item.str() + ":subdivided_minus_composed", sub - und, "|.| <= 1e-13 (r_combined = 1-(1-r1)(1-r2))", std::fabs(sub - und) <= 1e-13);
            row("R1_composition", item.str() + ":production_minus_brute", sub - brute, "|.| <= 1e-12 (oracle)", std::fabs(sub - brute) <= 1e-12);
            if (rr.first == 0.01 && rr.second == 0.01) {
                const double plain = panvar::log_haploid_mosaic_prior(two, p, {0.01});
                note("R1_composition", item.str() + ":subdivided_minus_undivided_at_0.01", sub - plain,
                     "constant r is not invariant: a neutral boundary changes the prior by this many nats");
            }
        }
    }
    // the same at sequence level: a shared segment written as invariant context vs as a one-allele block
    panvar::StructuralLocusInput ctx, blk;
    ctx.locus_id = "ctx";
    blk.locus_id = "blk";
    ctx.blocks = {{"b0", {{"a", "A"}, {"c", "C"}}}, {"b1", {{"g", "G"}, {"t", "T"}}}};
    ctx.invariant_segments = {"", "AAAA", ""};
    blk.blocks = {{"b0", {{"a", "A"}, {"c", "C"}}}, {"bm", {{"m", "AAAA"}}}, {"b1", {{"g", "G"}, {"t", "T"}}}};
    blk.invariant_segments = {"", "", "", ""};
    const std::vector<std::pair<std::string, std::vector<std::string>>> rows{{"h0", {"a", "g"}}, {"h1", {"c", "t"}}, {"h2", {"a", "t"}}};
    for (const auto& r : rows) {
        ctx.templates.push_back({r.first, panvar::TemplateFrameStatus::Complete, r.second});
        blk.templates.push_back({r.first, panvar::TemplateFrameStatus::Complete, {r.second[0], "m", r.second[1]}});
    }
    const auto mc = panvar::prepare_genotype_model(ctx), mb = panvar::prepare_genotype_model(blk);
    const auto sc = spellable(mc), sb = spellable(mb);
    row("R1_sequence_level", "same_spellable_sequence_set", static_cast<double>(sc == sb), "1 (the edit is neutral: no new admissible sequence)", sc == sb);
    const double scale = 50.0;
    const auto dc = distance_switches(mc, scale), db = distance_switches(mb, scale);
    for (const std::string& s : sc) {
        const double c01 = panvar::log_haploid_sequence_prior(mc, s, {0.01});
        const double b01 = panvar::log_haploid_sequence_prior(mb, s, {0.01, 0.01});
        const double c_comp = panvar::log_haploid_sequence_prior(mc, s, {1.0 - 0.99 * 0.99});
        row("R1_sequence_level", s + ":block_form_minus_composed_context_form(const)", b01 - c_comp, "|.| <= 1e-13", std::fabs(b01 - c_comp) <= 1e-13);
        note("R1_sequence_level", s + ":block_form_minus_context_form(const 0.01)", b01 - c01, "representation dependence of the constant prior");
        const double cd = panvar::log_haploid_sequence_prior(mc, s, dc), bd = panvar::log_haploid_sequence_prior(mb, s, db);
        row("R1_sequence_level", s + ":block_form_minus_context_form(distance, L fixed)", bd - cd, "|.| <= 1e-13 (subdivision-consistent)", std::fabs(bd - cd) <= 1e-13);
    }
}

// R2: a subdivision that is NOT neutral -- splitting a polymorphic block admits new recombinant sequences.
void r2_non_neutral_split() {
    panvar::StructuralLocusInput whole, split;
    whole.locus_id = "whole";
    split.locus_id = "split";
    whole.blocks = {{"b0", {{"ac", "AC"}, {"gt", "GT"}}}};
    whole.invariant_segments = {"N", "N"};
    split.blocks = {{"b0a", {{"a", "A"}, {"g", "G"}}}, {"b0b", {{"c", "C"}, {"t", "T"}}}};
    split.invariant_segments = {"N", "", "N"};
    whole.templates = {{"h0", panvar::TemplateFrameStatus::Complete, {"ac"}}, {"h1", panvar::TemplateFrameStatus::Complete, {"gt"}}};
    split.templates = {{"h0", panvar::TemplateFrameStatus::Complete, {"a", "c"}}, {"h1", panvar::TemplateFrameStatus::Complete, {"g", "t"}}};
    const auto mw = panvar::prepare_genotype_model(whole), ms = panvar::prepare_genotype_model(split);
    std::size_t finite_w = 0, finite_s = 0;
    for (const std::string& s : {std::string("NACN"), std::string("NGTN"), std::string("NATN"), std::string("NGCN")}) {
        finite_w += std::isfinite(panvar::log_haploid_sequence_prior(mw, s, {})) ? 1 : 0;
        finite_s += std::isfinite(panvar::log_haploid_sequence_prior(ms, s, {0.01})) ? 1 : 0;
    }
    row("R2_non_neutral", "admissible_sequences_before", static_cast<double>(finite_w), "2", finite_w == 2);
    row("R2_non_neutral", "admissible_sequences_after_split", static_cast<double>(finite_s), "4 (NATN, NGCN become admissible: NOT neutral)", finite_s == 4);
    note("R2_non_neutral", "log_prior_new_recombinant_NATN(r=0.01)", panvar::log_haploid_sequence_prior(ms, "NATN", {0.01}), "was -inf before the split");
}

// R3: summed-history vs best-history ranking (production = sum; oracle both), incl. at the production r = 0.01.
void r3_sum_vs_best() {
    std::mt19937_64 rng(31031);
    std::size_t comparisons = 0, disagreements = 0, production_mismatch = 0;
    for (int trial = 0; trial < 400; ++trial) {
        panvar::PanelAlleleMatrix panel;
        const std::size_t t = 3 + rng() % 4, b = 3 + rng() % 3;
        for (std::size_t h = 0; h < t; ++h) {
            panel.template_names.push_back("t" + std::to_string(h));
            std::vector<std::size_t> row;
            for (std::size_t k = 0; k < b; ++k) row.push_back(rng() % 2);
            panel.template_alleles.push_back(row);
        }
        for (std::size_t k = 0; k < b; ++k) { panel.template_alleles[0][k] = 0; panel.template_alleles[1][k] = 1; }   // both alleles present
        std::vector<panvar::MosaicPath> paths;
        for (std::size_t m = 0; m < (1u << b); ++m) {
            panvar::MosaicPath p;
            for (std::size_t k = 0; k < b; ++k) p.alleles.push_back((m >> k) & 1u);
            paths.push_back(p);
        }
        const std::vector<double> r(b - 1, 0.01);
        std::vector<long double> sum, best;
        for (const auto& p : paths) {
            sum.push_back(brute_haploid_prior(panel, p, r));
            best.push_back(brute_best_haploid_history(panel, p, r));
            const double prod = panvar::log_haploid_mosaic_prior(panel, p, r);
            if (sum.back() > 0 && std::fabs(prod - std::log(static_cast<double>(sum.back()))) > 1e-12) ++production_mismatch;
        }
        for (std::size_t i = 0; i < paths.size(); ++i) for (std::size_t j = i + 1; j < paths.size(); ++j) {
            if (sum[i] == 0 || sum[j] == 0) continue;
            ++comparisons;
            if ((sum[i] > sum[j]) != (best[i] > best[j]) && sum[i] != sum[j] && best[i] != best[j]) ++disagreements;
        }
    }
    row("R3_sum_vs_best", "production_equals_summed_oracle_mismatches", static_cast<double>(production_mismatch), "0", production_mismatch == 0);
    note("R3_sum_vs_best", "random_panels_pairs_compared(r=0.01)", static_cast<double>(comparisons), "400 random panels, 3-6 templates, 3-5 biallelic blocks");
    note("R3_sum_vs_best", "pairs_ranked_differently_sum_vs_best(r=0.01)", static_cast<double>(disagreements), "order flips between summed and best history");
}

// R5: metadata aliases vs genuine panel frequency.
void r5_aliases_and_multiplicity() {
    const auto build = [](bool alias, bool duplicate) {
        panvar::StructuralLocusInput in;
        in.locus_id = "alias";
        in.blocks = {{"b0", {{"a", "A"}, {"c", "C"}}}, {"b1", {{"g", "G"}, {"t", "T"}}}};
        if (alias) in.blocks[0].alleles.push_back({"a_alias", "A"});   // a second source id with identical bytes
        in.invariant_segments = {"", "", ""};
        in.templates = {{"h0", panvar::TemplateFrameStatus::Complete, {"a", "g"}},
                        {"h1", panvar::TemplateFrameStatus::Complete, {"c", "t"}},
                        {"h2", panvar::TemplateFrameStatus::Complete, {alias ? "a_alias" : "a", "t"}}};
        if (duplicate) in.templates.push_back({"h0_dup", panvar::TemplateFrameStatus::Complete, {"a", "g"}});
        return panvar::prepare_genotype_model(in);
    };
    const auto base = build(false, false), alias = build(true, false), dup = build(false, true);
    row("R5_alias", "alias_block_alleles_after_merge", static_cast<double>(alias.blocks()[0].alleles.size()), "2 (identical bytes merged)",
        alias.blocks()[0].alleles.size() == 2);
    for (const std::string& s : {std::string("AG"), std::string("AT"), std::string("CG"), std::string("CT")}) {
        const double b = panvar::log_haploid_sequence_prior(base, s, {0.01}), a = panvar::log_haploid_sequence_prior(alias, s, {0.01});
        row("R5_alias", s + ":alias_minus_base", (std::isinf(a) && std::isinf(b)) ? 0.0 : a - b, "0 (a source-id alias changes nothing)",
            (std::isinf(a) && std::isinf(b)) || std::fabs(a - b) <= 1e-14);
        const double d = panvar::log_haploid_sequence_prior(dup, s, {0.01});
        note("R5_multiplicity", s + ":duplicated_template_minus_base", (std::isinf(d) && std::isinf(b)) ? 0.0 : d - b,
             "a duplicated template row reweights the prior exactly as a genuine extra carrier would (indistinguishable to the model)");
    }
}

// R6: donor redraw vs actual donor change vs allele change (constant model), checked against the brute oracle.
void r6_redraw_semantics() {
    panvar::PanelAlleleMatrix panel;
    panel.template_names = {"h0", "h1", "h2", "h3"};
    panel.template_alleles = {{0, 0}, {0, 1}, {1, 1}, {1, 1}};
    const double r = 0.01, n = 4.0;
    const double p00 = std::exp(panvar::log_haploid_mosaic_prior(panel, panvar::MosaicPath{{0, 0}}, {r}));
    const double p01 = std::exp(panvar::log_haploid_mosaic_prior(panel, panvar::MosaicPath{{0, 1}}, {r}));
    // start on allele 0 (h0 or h1, each 1/4). P(next allele 1 | allele 0 now) = [ (1/4)(1-r)*[h1 carries 1] + (1/2)(r/n)*3 ] / (1/2)
    const double expected = ((0.25 * (1.0 - r)) + 0.5 * (r / n) * 3.0) / 0.5;
    row("R6_semantics", "P(allele 1 at block 1 | allele 0 at block 0)", p01 / (p00 + p01), "matches (1-r)*stay-carriers + r*carrier-share", std::fabs(p01 / (p00 + p01) - expected) <= 1e-13);
    note("R6_semantics", "redraw_probability_r", r, "the parameter: a donor REDRAW, uniform over all n templates (the current one included)");
    note("R6_semantics", "actual_donor_change_probability", r * (n - 1.0) / n, "r(n-1)/n");
    note("R6_semantics", "stay_on_same_donor", 1.0 - r + r / n, "(1-r) + r/n");
}

}  // namespace

int main() {
    try {
        // the existing fixtures first (regression): brute-force prior, sequence-equivalent aggregation, conceded counterexamples
        test_li_stephens_prior();
        test_sequence_equivalent_prior_aggregation();
        test_conceded_counterexamples();
        row("R0_existing_oracles", "mosaic_core prior fixtures", 1.0, "pass", true);
        r1_neutral_subdivision();
        r2_non_neutral_split();
        r3_sum_vs_best();
        r5_aliases_and_multiplicity();
        r6_redraw_semantics();
    } catch (const std::exception& e) {
        std::cerr << "m3a prior representation: ERROR " << e.what() << '\n';
        return 2;
    }
    std::cout << "test\titem\tvalue\texpectation\tresult\n" << g_rows.str() << "m3a_prior_representation\t" << (g_failed ? "FAIL" : "PASS") << '\n';
    return g_failed ? 1 : 0;
}
