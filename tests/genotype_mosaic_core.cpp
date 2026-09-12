#include "panvar/genotype_model.hpp"
#include "panvar/ls_hmm.hpp"
#include "panvar/mosaic_spelling.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

void require_close(double observed, double expected, double tolerance, const std::string& message) {
    if (!std::isfinite(observed) || !std::isfinite(expected) ||
        std::abs(observed - expected) > tolerance) {
        throw std::runtime_error(message + ": observed=" + std::to_string(observed) +
                                 " expected=" + std::to_string(expected));
    }
}

void require_throws(const std::function<void()>& operation, const std::string& message) {
    try {
        operation();
    } catch (const std::invalid_argument&) {
        return;
    }
    throw std::runtime_error(message);
}

long double brute_haploid_prior(
    const panvar::PanelAlleleMatrix& panel,
    const panvar::MosaicPath& mosaic,
    const std::vector<double>& switch_probability) {
    const std::size_t n = panel.template_count();
    long double total = 0.0L;

    std::function<void(std::size_t, std::size_t, long double)> extend;
    extend = [&](std::size_t block, std::size_t previous, long double probability) {
        if (block == panel.block_count()) {
            total += probability;
            return;
        }
        for (std::size_t h = 0; h < n; ++h) {
            if (panel.template_alleles[h][block] != mosaic.alleles[block]) continue;
            const long double r = switch_probability[block - 1];
            const long double transition =
                (h == previous ? 1.0L - r : 0.0L) + r / static_cast<long double>(n);
            extend(block + 1, h, probability * transition);
        }
    };

    for (std::size_t h = 0; h < n; ++h) {
        if (panel.template_alleles[h][0] != mosaic.alleles[0]) continue;
        extend(1, h, 1.0L / static_cast<long double>(n));
    }
    return total;
}

long double brute_best_haploid_history(
    const panvar::PanelAlleleMatrix& panel,
    const panvar::MosaicPath& mosaic,
    const std::vector<double>& switch_probability) {
    const std::size_t n = panel.template_count();
    long double best = 0.0L;
    std::function<void(std::size_t, std::size_t, long double)> extend;
    extend = [&](std::size_t block, std::size_t previous, long double probability) {
        if (block == panel.block_count()) {
            best = std::max(best, probability);
            return;
        }
        for (std::size_t h = 0; h < n; ++h) {
            if (panel.template_alleles[h][block] != mosaic.alleles[block]) continue;
            const long double r = switch_probability[block - 1];
            const long double transition =
                (h == previous ? 1.0L - r : 0.0L) + r / static_cast<long double>(n);
            extend(block + 1, h, probability * transition);
        }
    };
    for (std::size_t h = 0; h < n; ++h) {
        if (panel.template_alleles[h][0] != mosaic.alleles[0]) continue;
        extend(1, h, 1.0L / static_cast<long double>(n));
    }
    return best;
}

panvar::SequencingRead make_read(const std::string& name, std::size_t length) {
    return {name, std::string(length, 'A'), std::string(length, 'I')};
}

void test_library_model() {
    panvar::LibraryModel model;
    model.insert_mean_bp = 350.0;
    model.insert_sd_bp = 50.0;
    model.support_sigmas = 6.0;
    model.provenance = "fixture: declared";

    const panvar::ReadPair overlapping{make_read("r/1", 150), make_read("r/2", 150)};
    const panvar::InsertSupport support = panvar::fragment_insert_support(overlapping, model);
    require(support.min_bp == 150 && support.max_bp == 650,
            "overlapping mates must use max(read lengths), giving support [150,650]");

    const panvar::DiscreteInsertDistribution distribution =
        panvar::build_fragment_insert_distribution(overlapping, model);
    double probability_sum = 0.0;
    for (std::size_t length = support.min_bp; length <= support.max_bp; ++length) {
        probability_sum += std::exp(distribution.log_probability(length));
    }
    require_close(probability_sum, 1.0, 2e-14,
                  "discrete insert probabilities must sum to one on declared support");
    require(std::isinf(distribution.log_probability(149)) &&
                distribution.log_probability(149) < 0.0,
            "insert below physical support must have zero probability");

    const panvar::ReadPair variable{make_read("v/1", 100), make_read("v/2", 250)};
    const panvar::InsertSupport variable_support = panvar::fragment_insert_support(variable, model);
    require(variable_support.min_bp == 250 && variable_support.max_bp == 650,
            "insert floor must be fragment-specific for variable read lengths");

    const panvar::ReadPair unsupported{make_read("u/1", 700), make_read("u/2", 150)};
    require(panvar::fragment_insert_support(unsupported, model).empty(),
            "a read longer than declared insert support must be an explicit unsupported state");
    require(panvar::build_fragment_insert_distribution(unsupported, model)
                .log_probability_by_offset.empty(),
            "unsupported insert distribution must not expose a partial table");

    panvar::ReadPair bad_quality = overlapping;
    bad_quality.first.qualities.pop_back();
    require_throws([&] { panvar::validate_read_pair(bad_quality); },
                   "base/quality length mismatch must refuse");

    panvar::LibraryModel no_provenance = model;
    no_provenance.provenance.clear();
    require_throws([&] { panvar::validate_library_model(no_provenance); },
                   "an unprovenanced library model must refuse");
}

panvar::StructuralLocusInput structural_fixture(bool reordered) {
    panvar::StructuralLocusInput input;
    input.locus_id = "fixture-locus";
    input.chain_orientation = panvar::ChainOrientation::Reverse;
    input.invariant_segments = {"LEFT", "MIDDLE", "RIGHT"};
    input.blocks = {
        {"block-0", {{"walk-alt", "TT"}, {"walk-ref-2", "AC"}, {"walk-ref-1", "AC"}}},
        {"block-1", {{"walk-bypass", ""}, {"walk-one", "CC"}, {"walk-two", "GG"}}},
    };
    input.templates = {
        {"hap-z", panvar::TemplateFrameStatus::Complete, {"walk-alt", "walk-two"}},
        {"hap-a", panvar::TemplateFrameStatus::Complete, {"walk-ref-1", "walk-one"}},
        {"hap-m", panvar::TemplateFrameStatus::Complete, {"walk-ref-2", "walk-bypass"}},
    };
    input.templates[2].source_orientation = panvar::ChainOrientation::Reverse;
    if (reordered) {
        std::reverse(input.blocks[0].alleles.begin(), input.blocks[0].alleles.end());
        std::reverse(input.blocks[1].alleles.begin(), input.blocks[1].alleles.end());
        std::reverse(input.templates.begin(), input.templates.end());
    }
    return input;
}

void test_prepared_genotype_model() {
    static_assert(!std::is_default_constructible<panvar::PreparedGenotypeModel>::value,
                  "a prepared model must only exist after validated construction");

    const panvar::StructuralLocusInput input = structural_fixture(false);
    const panvar::PreparedGenotypeModel model = panvar::prepare_genotype_model(input);
    const panvar::PreparedGenotypeModel reordered =
        panvar::prepare_genotype_model(structural_fixture(true));

    require(model.locus_id() == "fixture-locus", "prepared model must retain its locus ID");
    require(model.chain_orientation() == panvar::ChainOrientation::Reverse,
            "prepared model must retain chain orientation");
    require(model.frame_status() == panvar::TemplateFrameStatus::Complete,
            "a constructed model must certify a complete panel frame");
    require(model.invariant_segments() == std::vector<std::string>({"LEFT", "MIDDLE", "RIGHT"}),
            "prepared model must retain every ordered invariant segment exactly once");
    require(model.blocks().size() == 2 && model.panel().block_count() == 2,
            "prepared structural and panel block counts must agree");
    require(model.panel().template_names == std::vector<std::string>({"hap-a", "hap-m", "hap-z"}),
            "template rows must have deterministic name order");

    // Sequence-identical walks are one statistical allele, but their graph aliases and panel
    // carrier multiplicity survive canonicalization. Empty sequence is a valid bypass allele.
    require(model.blocks()[0].alleles.size() == 2,
            "sequence-identical source walks must collapse to one canonical allele");
    const std::size_t ref_id = model.blocks()[0].allele_id_for_source("walk-ref-1");
    require(ref_id == model.blocks()[0].allele_id_for_source("walk-ref-2"),
            "sequence-identical source aliases must resolve to the same allele ID");
    require(model.allele(0, ref_id).source_ids ==
                std::vector<std::string>({"walk-ref-1", "walk-ref-2"}),
            "canonical allele must retain sorted source aliases");
    require(model.allele(0, ref_id).panel_carrier_count == 2,
            "canonical allele must retain panel carrier multiplicity");
    require(model.allele(1, model.blocks()[1].allele_id_for_source("walk-bypass")).sequence.empty(),
            "empty bypass sequence must remain an explicit allele rather than missing data");

    // Every input template/source mapping round-trips through the canonical block-local allele ID.
    std::map<std::string, const panvar::PanelTemplateInput*> input_template_by_name;
    for (const auto& input_template : input.templates) {
        input_template_by_name.emplace(input_template.name, &input_template);
    }
    for (std::size_t h = 0; h < model.panel().template_count(); ++h) {
        const auto found = input_template_by_name.find(model.panel().template_names[h]);
        require(found != input_template_by_name.end(), "prepared template must come from panel input");
        for (std::size_t b = 0; b < model.panel().block_count(); ++b) {
            const std::size_t expected =
                model.blocks()[b].allele_id_for_source(found->second->allele_source_ids[b]);
            require(model.panel().template_alleles[h][b] == expected,
                    "template/block/allele mapping must round-trip exactly");
            require(model.allele(b, expected).id == expected,
                    "canonical allele index and stored stable ID index must agree");
        }
    }

    // Reordering source records and template records cannot alter canonical IDs or the panel matrix.
    require(model.panel().template_names == reordered.panel().template_names &&
                model.panel().template_alleles == reordered.panel().template_alleles &&
                model.template_source_orientations() ==
                    reordered.template_source_orientations(),
            "panel mapping must be stable under input record reordering");
    require(model.blocks().size() == reordered.blocks().size(),
            "reordering must preserve the number of prepared blocks");
    for (std::size_t b = 0; b < model.blocks().size(); ++b) {
        require(model.blocks()[b].alleles.size() == reordered.blocks()[b].alleles.size(),
                "reordering must preserve canonical allele count");
        for (std::size_t a = 0; a < model.blocks()[b].alleles.size(); ++a) {
            const auto& lhs = model.blocks()[b].alleles[a];
            const auto& rhs = reordered.blocks()[b].alleles[a];
            require(lhs.id == rhs.id && lhs.stable_id == rhs.stable_id &&
                        lhs.sequence == rhs.sequence && lhs.source_ids == rhs.source_ids &&
                        lhs.panel_carrier_count == rhs.panel_carrier_count,
                    "canonical allele identity must be stable under input record reordering");
        }
    }

    panvar::StructuralLocusInput duplicate_template = input;
    duplicate_template.templates.push_back(duplicate_template.templates.front());
    require_throws([&] { panvar::prepare_genotype_model(duplicate_template); },
                   "duplicate template names must refuse");

    panvar::StructuralLocusInput missing_mapping = input;
    missing_mapping.templates.front().allele_source_ids[0] = "not-in-block";
    require_throws([&] { panvar::prepare_genotype_model(missing_mapping); },
                   "a template mapping to an absent block allele must refuse");

    panvar::StructuralLocusInput inconsistent_blocks = input;
    inconsistent_blocks.templates.front().allele_source_ids.pop_back();
    require_throws([&] { panvar::prepare_genotype_model(inconsistent_blocks); },
                   "a template missing one block mapping must refuse");

    panvar::StructuralLocusInput partial_frame = input;
    partial_frame.templates.front().frame_status = panvar::TemplateFrameStatus::Partial;
    require_throws([&] { panvar::prepare_genotype_model(partial_frame); },
                   "a partial panel frame must refuse explicitly");

    panvar::StructuralLocusInput duplicate_block = input;
    duplicate_block.blocks[1].block_id = duplicate_block.blocks[0].block_id;
    require_throws([&] { panvar::prepare_genotype_model(duplicate_block); },
                   "duplicate block names must refuse");

    panvar::StructuralLocusInput duplicate_source = input;
    duplicate_source.blocks[0].alleles[1].source_id =
        duplicate_source.blocks[0].alleles[0].source_id;
    require_throws([&] { panvar::prepare_genotype_model(duplicate_source); },
                   "duplicate source allele names within one block must refuse");

    panvar::StructuralLocusInput missing_context = input;
    missing_context.invariant_segments.pop_back();
    require_throws([&] { panvar::prepare_genotype_model(missing_context); },
                   "an incomplete invariant-context frame must refuse");
}

std::string independent_reverse_complement(const std::string& sequence) {
    std::string result;
    result.reserve(sequence.size());
    for (auto it = sequence.rbegin(); it != sequence.rend(); ++it) {
        switch (*it) {
            case 'A': result.push_back('T'); break;
            case 'C': result.push_back('G'); break;
            case 'G': result.push_back('C'); break;
            case 'T': result.push_back('A'); break;
            default: result.push_back('N'); break;
        }
    }
    return result;
}

std::size_t template_id(const panvar::PreparedGenotypeModel& model, const std::string& name) {
    const auto found = std::find(
        model.panel().template_names.begin(), model.panel().template_names.end(), name);
    if (found == model.panel().template_names.end()) {
        throw std::runtime_error("fixture template not found: " + name);
    }
    return static_cast<std::size_t>(found - model.panel().template_names.begin());
}

void test_mosaic_spelling() {
    const panvar::StructuralLocusInput input = structural_fixture(false);
    const panvar::PreparedGenotypeModel model = panvar::prepare_genotype_model(input);

    const panvar::HaplotypeSequence hap_a =
        panvar::spell_panel_template(model, template_id(model, "hap-a"));
    const panvar::HaplotypeSequence hap_m =
        panvar::spell_panel_template(model, template_id(model, "hap-m"));
    const panvar::HaplotypeSequence hap_z =
        panvar::spell_panel_template(model, template_id(model, "hap-z"));
    require(hap_a.bases == "LEFTACMIDDLECCRIGHT" &&
                hap_m.bases == "LEFTACMIDDLERIGHT" &&
                hap_z.bases == "LEFTTTMIDDLEGGRIGHT",
            "every panel template must spell in canonical chain orientation");

    const std::vector<panvar::AuthoritativePanelSequence> authoritative = {
        {"hap-z", "LEFTTTMIDDLEGGRIGHT"},
        {"hap-m", independent_reverse_complement("LEFTACMIDDLERIGHT")},
        {"hap-a", "LEFTACMIDDLECCRIGHT"},
    };
    panvar::verify_complete_panel_round_trip(model, authoritative);

    // This allele combination is carried by no panel template. It must still spell directly from
    // the block alphabet, with every invariant segment present once and the bypass allele explicit.
    const panvar::MosaicPath recombinant{{
        model.blocks()[0].allele_id_for_source("walk-alt"),
        model.blocks()[1].allele_id_for_source("walk-bypass"),
    }};
    const panvar::HaplotypeSequence novel = panvar::spell_mosaic(model, recombinant);
    require(novel.bases == "LEFTTTMIDDLERIGHT",
            "off-panel recombinant must spell the hand-computed locus sequence");
    require(novel.invariant_intervals ==
                std::vector<panvar::SequenceInterval>({{0, 4}, {6, 12}, {12, 17}}),
            "invariant contexts must appear exactly once at their expected coordinates");
    require(novel.block_intervals ==
                std::vector<panvar::SequenceInterval>({{4, 6}, {12, 12}}),
            "block intervals must cover variable lengths and a zero-width bypass exactly");

    // An antiparallel source path has mirrored source coordinates. Comparing this exact result to
    // the raw chain coordinate makes the historical raw-walk-coordinate mutation non-vacuous.
    const panvar::SequenceInterval chain_block = hap_m.block_intervals[0];
    const panvar::SequenceInterval source_block = panvar::interval_in_source_frame(
        chain_block, hap_m.bases.size(), panvar::ChainOrientation::Reverse);
    require(chain_block == panvar::SequenceInterval{4, 6},
            "fixture must pin the forward chain block interval");
    require(source_block == panvar::SequenceInterval{11, 13} && source_block != chain_block,
            "reverse source frame must mirror rather than reuse raw chain coordinates");
    require(panvar::sequence_in_source_frame(
                hap_m.bases, panvar::ChainOrientation::Reverse) == authoritative[1].source_frame_bases,
            "reverse source frame must be the exact reverse complement of chain spelling");

    // Source record names and graph-route multiplicity cannot affect the spelled statistical state.
    panvar::StructuralLocusInput alternative_topology = input;
    alternative_topology.blocks[0].alleles.push_back({"other-route-ref", "AC"});
    alternative_topology.blocks[1].alleles.push_back({"other-route-one", "CC"});
    for (auto& panel_template : alternative_topology.templates) {
        if (panel_template.name == "hap-a") {
            panel_template.allele_source_ids = {"other-route-ref", "other-route-one"};
        }
    }
    const panvar::PreparedGenotypeModel alternative =
        panvar::prepare_genotype_model(alternative_topology);
    require(model.panel().template_names == alternative.panel().template_names &&
                model.panel().template_alleles == alternative.panel().template_alleles,
            "sequence-equivalent graph route aliases must not change the HMM state matrix");
    for (const std::string& name : model.panel().template_names) {
        require(panvar::spell_panel_template(model, template_id(model, name)).bases ==
                    panvar::spell_panel_template(alternative, template_id(alternative, name)).bases,
                "sequence-equivalent graph route aliases must not change panel spelling");
    }
    require(panvar::spell_mosaic(alternative, recombinant).bases == novel.bases,
            "sequence-equivalent graph route aliases must not change recombinant spelling");

    // The round-trip, rather than construction alone, is what detects a context duplicated or lost
    // during graph decomposition.
    panvar::StructuralLocusInput omitted_context_input = input;
    omitted_context_input.invariant_segments[1].clear();
    const panvar::PreparedGenotypeModel omitted_context =
        panvar::prepare_genotype_model(omitted_context_input);
    require_throws(
        [&] { panvar::verify_complete_panel_round_trip(omitted_context, authoritative); },
        "omitting an invariant context must fail authoritative panel round-trip");

    panvar::StructuralLocusInput duplicated_context_input = input;
    duplicated_context_input.invariant_segments[1] += input.invariant_segments[1];
    const panvar::PreparedGenotypeModel duplicated_context =
        panvar::prepare_genotype_model(duplicated_context_input);
    require_throws(
        [&] { panvar::verify_complete_panel_round_trip(duplicated_context, authoritative); },
        "duplicating an invariant context must fail authoritative panel round-trip");

    require_throws(
        [&] { panvar::spell_mosaic(model, panvar::MosaicPath{{0}}); },
        "mosaic missing a block allele must refuse");
    require_throws(
        [&] { panvar::spell_mosaic(model, panvar::MosaicPath{{0, 99}}); },
        "mosaic with an out-of-range allele must refuse");
}

void test_li_stephens_prior() {
    panvar::PanelAlleleMatrix panel;
    panel.template_names = {"h0", "h1", "h2"};
    panel.template_alleles = {
        {0, 0, 1, 0},
        {0, 1, 1, 1},
        {1, 0, 0, 1},
    };
    panvar::MosaicPath recombinant{{0, 0, 1, 1}};

    for (const std::vector<double>& recombination :
         std::vector<std::vector<double>>{{0.0, 0.0, 0.0},
                                          {1.0, 1.0, 1.0},
                                          {0.05, 0.3, 0.9}}) {
        const long double brute = brute_haploid_prior(panel, recombinant, recombination);
        const double exact =
            panvar::log_haploid_mosaic_prior(panel, recombinant, recombination);
        if (brute == 0.0L) {
            require(std::isinf(exact) && exact < 0.0,
                    "an unreachable exhaustive mosaic must have exactly zero prior");
            continue;
        }
        require_close(exact, std::log(static_cast<double>(brute)), 2e-14,
                      "O(blocks*templates) mosaic prior must equal exhaustive template histories");
    }

    const std::vector<double> intermediate{0.05, 0.3, 0.9};
    const long double summed_history = brute_haploid_prior(panel, recombinant, intermediate);
    const long double best_history = brute_best_haploid_history(panel, recombinant, intermediate);
    require(summed_history > best_history * 1.01L,
            "fixture must distinguish exact marginalization from the best template history");
    const panvar::MosaicPath second{{1, 0, 0, 1}};
    const double first_prior =
        panvar::log_haploid_mosaic_prior(panel, recombinant, intermediate);
    const double second_prior =
        panvar::log_haploid_mosaic_prior(panel, second, intermediate);

    panvar::MosaicDiplotype heterozygous{recombinant, second};
    require_close(
        panvar::log_unordered_diplotype_prior(panel, heterozygous, intermediate),
        first_prior + second_prior + std::log(2.0),
        2e-14,
        "distinct unordered mosaics must aggregate both homologue orders exactly once");

    panvar::MosaicDiplotype homozygous{recombinant, recombinant};
    require_close(
        panvar::log_unordered_diplotype_prior(panel, homozygous, intermediate),
        2.0 * first_prior,
        2e-14,
        "identical mosaics must not receive a spurious log(2) multiplicity");

    panvar::MosaicPath absent{{9, 9, 9, 9}};
    const double absent_prior = panvar::log_haploid_mosaic_prior(panel, absent, intermediate);
    require(std::isinf(absent_prior) && absent_prior < 0.0,
            "a mosaic with an unavailable allele must have zero prior");

    require_throws(
        [&] { panvar::log_haploid_mosaic_prior(panel, recombinant, {0.2, 0.2}); },
        "missing boundary probability must refuse");
    require_throws(
        [&] { panvar::log_haploid_mosaic_prior(panel, recombinant, {0.2, -0.1, 0.2}); },
        "invalid boundary probability must refuse");
}

} // namespace

int main() {
    try {
        test_library_model();
        test_prepared_genotype_model();
        test_mosaic_spelling();
        test_li_stephens_prior();
        std::cout << "genotype mosaic core: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "genotype mosaic core: FAIL: " << error.what() << '\n';
        return 1;
    }
}
