// Regression fixtures for the 2026-10-14 short-read code review (docs/reports/genotype-shortread-code-review-20261014.md).
//
// Each case evaluates the CORRECT property on a minimal adversarial input and prints `<case> holds=<0|1>`. The test
// passes when every observed value equals tests/genotype_review_regressions.expected, which records the current
// behaviour: a defect case is expected to read holds=0 until it is fixed, and its control must read holds=1 (so the
// property check itself is shown to work). Fixing a defect therefore changes this test's expectation file in the same
// commit, which makes the behaviour change explicit and reviewable. No case asserts which historical calls a defect
// affected.
//
//   C1 marker multiplicity: every retained anchor's inventory multiplicity equals its multiplicity in the complete
//      spelling of every haplotype the model admits (panel templates and, for these small models, every mosaic).
//      C1_panel_two_block: the review's reproduction (33 A + an A/C allele); C1_multi_short_block: a junction crossing
//      several short blocks, admitted only by mosaics; C1_control: the same shapes with a separating base.
//   C2 projection-aware annotation: when a lower-ranked sequence group has an admissible projection carrying a
//      block genotype different from the called one, the called block cannot have GQ = +inf (it must carry error mass).
//      C2_lower_group_projection: the review's reproduction; C2_control: the alternative genotype is the lower group's
//      representative, so the existing code accounts for it.
//   C3 dosage route invariance: two models with the same fingerprint (same canonical alleles and panel), encoded as a
//      repeated node versus two sequence-identical nodes, get the same census route.
//      C3_dosage_route_representation: the review's reproduction; C3_control: two identical encodings.
//
// Usage: genotype_review_regressions <expected file>
#include "panvar/genotype_dosage.hpp"
#include "panvar/genotype_locus.hpp"
#include "panvar/genotype_model.hpp"
#include "panvar/genotype_repair.hpp"
#include "panvar/marker_model.hpp"
#include "panvar/mosaic_spelling.hpp"
#include "panvar/syncmer.hpp"

#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {

using Blocks = std::vector<std::vector<std::string>>;   // per block: allele sequences

panvar::PreparedGenotypeModel model_of(const Blocks& blocks, const std::vector<std::vector<std::size_t>>& rows) {
    panvar::StructuralLocusInput input;
    input.locus_id = "review-regression";
    for (std::size_t b = 0; b < blocks.size(); ++b) {
        panvar::GenotypeBlockInput block;
        block.block_id = (blocks[b].size() == 1 ? "ctx:" : "sv:") + std::to_string(b);
        block.kind = blocks[b].size() == 1 ? panvar::GenotypeBlockKind::Context : panvar::GenotypeBlockKind::Variant;
        for (std::size_t a = 0; a < blocks[b].size(); ++a) {
            panvar::BlockAlleleInput allele;
            allele.source_id = "b" + std::to_string(b) + "a" + std::to_string(a);
            allele.sequence = blocks[b][a];
            block.alleles.push_back(allele);
        }
        input.blocks.push_back(block);
    }
    input.invariant_segments.assign(blocks.size() + 1, "");
    for (std::size_t t = 0; t < rows.size(); ++t) {
        panvar::PanelTemplateInput row;
        row.name = "t" + std::to_string(t);
        row.frame_status = panvar::TemplateFrameStatus::Complete;
        for (std::size_t b = 0; b < blocks.size(); ++b) {
            row.allele_source_ids.push_back("b" + std::to_string(b) + "a" + std::to_string(rows[t][b]));
        }
        input.templates.push_back(row);
    }
    return panvar::prepare_genotype_model(input);
}

// C1: true when every anchor's inventory multiplicity equals its count in every admitted haplotype's spelling.
bool anchors_consistent(const panvar::PreparedGenotypeModel& model, bool include_mosaics) {
    const auto markers = panvar::build_marker_panel(model, 31, 15);
    std::vector<std::string> haplotypes;
    for (std::size_t t = 0; t < model.panel().template_count(); ++t) {
        haplotypes.push_back(panvar::spell_panel_template(model, t).bases);
    }
    if (include_mosaics) {
        const std::size_t nb = model.blocks().size();
        std::vector<std::size_t> path(nb, 0);
        std::function<void(std::size_t)> walk = [&](std::size_t b) {
            if (b == nb) {
                panvar::MosaicPath m;
                m.alleles = path;
                haplotypes.push_back(panvar::spell_mosaic(model, m).bases);
                return;
            }
            for (std::size_t a = 0; a < model.blocks()[b].alleles.size(); ++a) {
                path[b] = a;
                walk(b + 1);
            }
        };
        walk(0);
    }
    for (std::size_t j = 0; j < markers.anchor_slots.size(); ++j) {
        const auto code = markers.codes[markers.anchor_slots[j]];
        for (const std::string& h : haplotypes) {
            std::size_t count = 0;
            for (const auto& k : panvar::collect_syncmers(h, 31, 15)) count += k.code == code;
            if (count != markers.anchor_haploid_multiplicity[j]) return false;
        }
    }
    return true;
}

panvar::PoolCandidate candidate(std::vector<std::size_t> first, std::vector<std::size_t> second, const std::string& s1,
                                const std::string& s2, double prior, double path_prior) {
    panvar::PoolCandidate c;
    c.first = std::move(first);
    c.second = std::move(second);
    c.first_sequence = panvar::sequence_digest(s1);
    c.second_sequence = panvar::sequence_digest(s2);
    c.score.log_likelihood = 0;
    c.score.fragment_term = 0;
    c.score.exposure_term = 0;
    c.sequence_prior = prior;
    c.path_prior = path_prior;
    return c;
}

panvar::PreparedGenotypeModel repeat_model(bool unrolled) {
    panvar::StructuralLocusInput x;
    x.locus_id = "same-sequence-repeat-representation";
    x.blocks = {{"sv:0", {{"walk:x+", std::string(600, 'A')}, {unrolled ? "walk:x+,y+" : "walk:x+,x+", std::string(1200, 'A')}},
                 panvar::GenotypeBlockKind::Variant},
                {"ctx:1", {{"walk:z+", "CCC"}}, panvar::GenotypeBlockKind::Context}};
    x.invariant_segments = {"", "", ""};
    for (const auto& [name, first] : {std::pair<std::string, std::string>{"short", "walk:x+"},
                                      {"long", unrolled ? "walk:x+,y+" : "walk:x+,x+"}}) {
        panvar::PanelTemplateInput row;
        row.name = name;
        row.frame_status = panvar::TemplateFrameStatus::Complete;
        row.allele_source_ids = {first, "walk:z+"};
        x.templates.push_back(row);
    }
    return panvar::prepare_genotype_model(x);
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "usage: genotype_review_regressions <expected file>\n";
        return 2;
    }
    std::map<std::string, int> observed;

    // ---- C1
    observed["C1_panel_two_block"] = anchors_consistent(model_of({{std::string(33, 'A')}, {"A", "C"}}, {{0, 0}, {0, 1}}), false);
    // a G-delimited 40 A run (the all-A 31-mer is an anchor there), then several short blocks whose A alleles join
    // with A-runs of the flanking contexts into extra all-A 31-mers; no panel template takes all A alleles, a mosaic can
    const Blocks multi = {{"G" + std::string(40, 'A') + "G"}, {std::string(10, 'A')}, {"AAAAAAAA", "C"}, {"AAAAAAAA", "T"},
                          {std::string(10, 'A') + "G"}};
    const std::vector<std::vector<std::size_t>> multi_rows = {{0, 0, 0, 1, 0}, {0, 0, 1, 0, 0}, {0, 0, 1, 1, 0}};
    observed["C1_multi_short_block_panel_only"] = anchors_consistent(model_of(multi, multi_rows), false);
    observed["C1_multi_short_block"] = anchors_consistent(model_of(multi, multi_rows), true);
    observed["C1_control"] = anchors_consistent(
        model_of({{"G" + std::string(40, 'A') + "G"}, {"CA", "CC"}, {"TTAG"}}, {{0, 0, 0}, {0, 1, 0}}), true);

    // ---- C2 (block 0 alleles A / AA; block 1 alleles C / AC). Winner AC + AAAC = (A,C)/(AA,AC).
    {
        const auto best = candidate({0, 0}, {1, 1}, "AC", "AAAC", 0, 0);
        auto runner = candidate({0, 1}, {1, 0}, "AAC", "AAC", -1, -1);     // projection (A,AC)/(AA,C): winner's genotypes
        auto alternative = runner;                                        // projection (AA,C)/(AA,C): another genotype
        alternative.first = alternative.second = std::vector<std::size_t>{1, 0};
        alternative.path_prior = -2;
        const auto d = panvar::decide_pooled_genotype({best, runner, alternative}, 2, 1, 0);
        observed["C2_lower_group_projection"] = !(d.has_other_genotype && std::isinf(d.block_genotype_quality[0]));
        // control: the alternative projection is the lower group's representative (greater path prior)
        alternative.path_prior = 0;
        runner.path_prior = -2;
        const auto c = panvar::decide_pooled_genotype({best, alternative, runner}, 2, 1, 0);
        observed["C2_control"] = !std::isinf(c.block_genotype_quality[0]);
    }

    // ---- C3
    {
        const panvar::dosage::NodeSequences nodes = {{"x", std::string(600, 'A')}, {"y", std::string(600, 'A')}, {"z", "CCC"}};
        const panvar::dosage::DosageSettings ds;
        const auto repeated = repeat_model(false), unrolled = repeat_model(true);
        const bool same_fp = panvar::model_fingerprint(repeated) == panvar::model_fingerprint(unrolled);
        const auto r1 = panvar::dosage::run_census(repeated, nodes, 500, ds);
        const auto r2 = panvar::dosage::run_census(unrolled, nodes, 500, ds);
        observed["C3_dosage_route_representation"] = same_fp && r1[0].route == r2[0].route;
        const auto r3 = panvar::dosage::run_census(repeat_model(false), nodes, 500, ds);
        observed["C3_control"] = r1[0].route == r3[0].route;
        std::cout << "# C3 same_fingerprint=" << same_fp << " routes " << panvar::dosage::census_route_name(r1[0].route)
                  << " / " << panvar::dosage::census_route_name(r2[0].route) << '\n';
    }

    std::map<std::string, int> expected;
    std::ifstream in(argv[1]);
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        std::string name;
        int value = -1;
        ss >> name >> value;
        expected[name] = value;
    }
    int failures = 0;
    for (const auto& [name, value] : observed) {
        const auto e = expected.find(name);
        const bool ok = e != expected.end() && e->second == value;
        failures += !ok;
        std::cout << name << " holds=" << value << " expected=" << (e == expected.end() ? -1 : e->second)
                  << (ok ? "" : "  MISMATCH") << '\n';
    }
    for (const auto& [name, value] : expected) {
        if (!observed.count(name)) {
            ++failures;
            std::cout << name << " expected but not evaluated  MISMATCH\n";
        }
    }
    std::cout << failures << " mismatch(es)\n";
    return failures == 0 ? 0 : 1;
}
