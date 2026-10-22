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
#include <utility>
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
    model.substitution_error = 0.01;
    model.outlier_weight = 0.001;
    model.fragment_start_intensity = 0.2;
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

    // Every scoring parameter is part of the same declared policy, so each must refuse on its own.
    for (const auto& mutation : std::vector<std::pair<std::string, panvar::LibraryModel>>{
             {"zero substitution error", [&] { auto m = model; m.substitution_error = 0.0; return m; }()},
             {"uninformative substitution error",
              [&] { auto m = model; m.substitution_error = 0.75; return m; }()},
             {"outlier weight of one", [&] { auto m = model; m.outlier_weight = 1.0; return m; }()},
             {"negative outlier weight", [&] { auto m = model; m.outlier_weight = -0.1; return m; }()},
             {"zero fragment-start intensity",
              [&] { auto m = model; m.fragment_start_intensity = 0.0; return m; }()},
         }) {
        require_throws([&] { panvar::validate_library_model(mutation.second); },
                       "library model must refuse a " + mutation.first);
    }

    // The background term depends only on the fragment's own read lengths, so it is the same for
    // every candidate and can never reorder them.
    const panvar::ReadPair shorter{make_read("s/1", 100), make_read("s/2", 100)};
    require_close(panvar::log_background_probability(overlapping), 300.0 * std::log(0.25), 1e-12,
                  "background probability must be uniform over both mates' bases");
    require(panvar::log_background_probability(shorter) >
                panvar::log_background_probability(overlapping),
            "a shorter fragment must have a larger background probability");

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


// Two different allele paths can spell identical complete sequences, and a read-based method cannot
// tell them apart. The prior of a sequence must aggregate every such path, and the unordered pair's
// log(2) must be decided on bytes. The oracle here is brute force over every allele path and every
// ordered path pair, grouped by spelled bytes -- no forward recursion involved.
void test_sequence_equivalent_prior_aggregation() {
    const double kNegInf = -std::numeric_limits<double>::infinity();
    panvar::StructuralLocusInput input;
    input.locus_id = "sequence-equivalence";
    input.blocks = {
        {"b0", {{"a", "A"}, {"at", "AT"}}},
        {"b1", {{"tg", "TG"}, {"g", "G"}}},
    };
    input.invariant_segments = {"C", "", "C"};
    input.templates = {
        {"h0", panvar::TemplateFrameStatus::Complete, {"a", "tg"}},    // CATGC
        {"h1", panvar::TemplateFrameStatus::Complete, {"at", "g"}},    // CATGC, a different path
        {"h2", panvar::TemplateFrameStatus::Complete, {"a", "g"}},     // CAGC
        {"h3", panvar::TemplateFrameStatus::Complete, {"at", "tg"}},   // CATTGC
    };
    const panvar::PreparedGenotypeModel model = panvar::prepare_genotype_model(input);

    for (double r : {0.0, 0.3, 1.0}) {
        const std::vector<double> switches = {r};
        std::vector<panvar::MosaicPath> paths;
        for (std::size_t a0 = 0; a0 < 2; ++a0) {
            for (std::size_t a1 = 0; a1 < 2; ++a1) paths.push_back(panvar::MosaicPath{{a0, a1}});
        }
        std::map<std::string, double> by_sequence;
        std::vector<std::string> spelled;
        std::vector<double> path_prior;
        for (const panvar::MosaicPath& path : paths) {
            spelled.push_back(panvar::spell_mosaic(model, path).bases);
            path_prior.push_back(panvar::log_haploid_mosaic_prior(model.panel(), path, switches));
            double& slot = by_sequence.emplace(spelled.back(), kNegInf).first->second;
            slot = panvar::log_add_exp(slot, path_prior.back());
        }
        require(by_sequence.count("CATGC") == 1 && spelled.size() == 4 && by_sequence.size() == 3,
                "fixture: exactly two distinct allele paths must spell CATGC");

        for (const auto& item : by_sequence) {
            const double dp = panvar::log_haploid_sequence_prior(model, item.first, switches);
            if (item.second == kNegInf) {
                require(dp == kNegInf, "an unreachable sequence must have zero prior");
            } else {
                require_close(dp, item.second, 1e-12,
                              "sequence prior must equal the sum over every path spelling it");
            }
        }
        require(panvar::log_haploid_sequence_prior(model, "CAAC", switches) == kNegInf,
                "a sequence no allele path spells must have zero prior");

        // Ordered path pairs, grouped by the unordered pair of spelled sequences.
        std::map<std::pair<std::string, std::string>, double> by_pair;
        for (std::size_t i = 0; i < paths.size(); ++i) {
            for (std::size_t j = 0; j < paths.size(); ++j) {
                std::pair<std::string, std::string> key = {spelled[i], spelled[j]};
                if (key.second < key.first) std::swap(key.first, key.second);
                double& slot = by_pair.emplace(key, kNegInf).first->second;
                slot = panvar::log_add_exp(slot, path_prior[i] + path_prior[j]);
            }
        }
        for (const auto& item : by_pair) {
            const double expected = item.second;
            const double actual = panvar::log_unordered_sequence_pair_prior(
                model, item.first.first, item.first.second, switches);
            if (expected == kNegInf) {
                require(actual == kNegInf, "an unreachable pair must have zero prior");
            } else {
                require_close(actual, expected, 1e-12,
                              "unordered sequence-pair prior must equal the brute-force sum");
            }
        }

        // Mutation: deciding log(2) on allele paths. The two CATGC paths form a SEQUENCE homozygote,
        // which the path-based prior scores as a heterozygote with log(2) and without aggregation.
        if (r > 0.0 && r < 1.0) {
            const double path_based = panvar::log_unordered_diplotype_prior(
                model.panel(), panvar::MosaicDiplotype{paths[1], paths[2]}, switches);
            const double sequence_based =
                panvar::log_unordered_sequence_pair_prior(model, "CATGC", "CATGC", switches);
            require(spelled[1] == "CATGC" && spelled[2] == "CATGC",
                    "fixture: paths 1 and 2 must be the two CATGC representations");
            require(std::abs(path_based - sequence_based) > 1e-6,
                    "the path-based prior must disagree on a sequence homozygote, or this gate is empty");
        }
    }
}

} // namespace

// The conceded review points (Astra, 2026-10-15), as mathematical counterexamples on the production prior and the
// existing brute-force history helpers. Each has a negative control that a wrong implementation or a wrong premise
// would fail. Verified quantities are named in each message.
void test_conceded_counterexamples() {
    // (a) Distant evidence does not factorize along an adjacent chain. One fragment with additive placement mass
    // m1(a1) + m3(a3) on blocks 1 and 3 (block 2 monomorphic between them) contributes c*(m1 + m3) + beta. With
    // m1 = (1, 2), m3 = (2, 3), c = 1, beta = 0 that factor is [[3,4],[4,5]]: rank 2, so no g(a1)*h(a3) through
    // the single block-2 state equals it. Control: a distant component that is invariant across the pool (m3 equal
    // for both alleles) leaves a rank-1 factor, which does factorize -- invariant contributions only shift beta.
    {
        const auto factor = [](const double m1[2], const double m3[2], double c, double beta) {
            std::vector<std::vector<double>> f(2, std::vector<double>(2));
            for (int a = 0; a < 2; ++a) for (int b = 0; b < 2; ++b) f[a][b] = c * (m1[a] + m3[b]) + beta;
            return f;
        };
        const double m1[2] = {1.0, 2.0}, varying[2] = {2.0, 3.0}, invariant[2] = {2.5, 2.5};
        const auto f = factor(m1, varying, 1.0, 0.0);
        require(f[0][0] == 3.0 && f[0][1] == 4.0 && f[1][0] == 4.0 && f[1][1] == 5.0, "(a) fixture is [[3,4],[4,5]]");
        const double det = f[0][0] * f[1][1] - f[0][1] * f[1][0];
        require(det != 0.0, "(a) a varying distant component must give a rank-2 (non-factorizable) evidence factor");
        // the best attempt through one chain state: g(a1) = f[a1][0], h(a3) = f[0][a3] / f[0][0] -- fails at (1,1)
        require(std::abs(f[1][0] * f[0][1] / f[0][0] - f[1][1]) > 0.3, "(a) chain reconstruction must fail");
        const auto g = factor(m1, invariant, 1.0, 0.0);
        require(std::abs(g[0][0] * g[1][1] - g[0][1] * g[1][0]) < 1e-12,
                "(a) control: an invariant distant component must leave a rank-1 factor");
    }

    // (b) Allele states are not Markov under Li-Stephens copying. Donors 000 and 101: after a1 = 1, a2 = 0 the copied
    // donor is most likely 101, so a3 = 1 is likely; after a1 = 0, a2 = 0 it is most likely 000. P(a3 | a2, a1) then
    // depends on a1. Control: r = 1 copies every block independently, and the dependence vanishes.
    {
        panvar::PanelAlleleMatrix panel;
        panel.template_names = {"d000", "d101"};
        panel.template_alleles = {{0, 0, 0}, {1, 0, 1}};
        const auto conditional = [&](double r, std::size_t a1) {
            const std::vector<double> rr{r, r};
            const double p1 = std::exp(panvar::log_haploid_mosaic_prior(panel, panvar::MosaicPath{{a1, 0, 1}}, rr));
            const double p0 = std::exp(panvar::log_haploid_mosaic_prior(panel, panvar::MosaicPath{{a1, 0, 0}}, rr));
            return p1 / (p0 + p1);   // P(a3 = 1 | a1, a2 = 0)
        };
        const double after1 = conditional(0.1, 1), after0 = conditional(0.1, 0);
        require(after1 > 0.9 && after0 < 0.1,
                "(b) P(a3=1 | a2=0) must depend on a1 under copying (allele states are not a Markov chain)");
        require(std::abs(conditional(1.0, 1) - conditional(1.0, 0)) < 1e-12,
                "(b) control: at r = 1 the allele process is independent and the dependence must vanish");
    }

    // (c) Summed history probability and best history rank sequences differently. Panel (5 templates, 4 blocks,
    // r = 0.9): 0101 has the larger summed prior but 0001 (a duplicated donor) has the larger single history. The
    // production prior must rank as the sum (it is the sum, brute-checked), so a best-history implementation is
    // detectably wrong on this fixture. Control: at r = 0 every sequence is one donor copy and the two orders agree.
    {
        panvar::PanelAlleleMatrix panel;
        panel.template_names = {"t0", "t1", "t2", "t3", "t4"};
        panel.template_alleles = {{1, 1, 1, 0}, {0, 0, 0, 1}, {0, 1, 1, 1}, {0, 0, 0, 1}, {1, 1, 0, 1}};
        const std::vector<double> r(3, 0.9);
        const panvar::MosaicPath summed_winner{{0, 1, 0, 1}}, best_history_winner{{0, 0, 0, 1}};
        const long double sum_a = brute_haploid_prior(panel, summed_winner, r);
        const long double sum_b = brute_haploid_prior(panel, best_history_winner, r);
        const long double best_a = brute_best_haploid_history(panel, summed_winner, r);
        const long double best_b = brute_best_haploid_history(panel, best_history_winner, r);
        require(sum_a > sum_b * 1.1L && best_b > best_a * 1.5L, "(c) fixture must flip the sum and max orders");
        require(panvar::log_haploid_mosaic_prior(panel, summed_winner, r) >
                    panvar::log_haploid_mosaic_prior(panel, best_history_winner, r),
                "(c) the production prior must rank by summed histories, not by the best history");
        const std::vector<double> r0(3, 0.0);
        const panvar::MosaicPath donor{{0, 0, 0, 1}}, other{{1, 1, 1, 0}};
        require((brute_haploid_prior(panel, donor, r0) > brute_haploid_prior(panel, other, r0)) ==
                    (brute_best_haploid_history(panel, donor, r0) >= brute_best_haploid_history(panel, other, r0)),
                "(c) control: at r = 0 the summed and best-history orders must agree");
    }

    // (d) Neutral subdivision is not neutral. Inserting a monomorphic block between two blocks composes two copying
    // steps of r into one of 1 - (1 - r)^2: 0.01 becomes 0.0199. Control: r = 1 is a fixed point of the composition.
    {
        panvar::PanelAlleleMatrix two, three;
        two.template_names = three.template_names = {"h0", "h1", "h2"};
        two.template_alleles = {{0, 0}, {1, 1}, {0, 1}};
        three.template_alleles = {{0, 7, 0}, {1, 7, 1}, {0, 7, 1}};
        for (const panvar::MosaicPath& p : {panvar::MosaicPath{{1, 0}}, panvar::MosaicPath{{0, 0}}}) {
            const panvar::MosaicPath p3{{p.alleles[0], 7, p.alleles[1]}};
            const double subdivided = panvar::log_haploid_mosaic_prior(three, p3, {0.01, 0.01});
            require_close(subdivided, panvar::log_haploid_mosaic_prior(two, p, {0.0199}), 1e-14,
                          "(d) a monomorphic subdivision must compose r = 0.01 into 0.0199");
            require(std::abs(subdivided - panvar::log_haploid_mosaic_prior(two, p, {0.01})) > 1e-4,
                    "(d) the subdivided prior must differ from the undivided one at r = 0.01");
            require_close(panvar::log_haploid_mosaic_prior(three, p3, {1.0, 1.0}),
                          panvar::log_haploid_mosaic_prior(two, p, {1.0}), 1e-14,
                          "(d) control: r = 1 is invariant under subdivision");
        }
    }
}

int main() {
    try {
        test_library_model();
        test_prepared_genotype_model();
        test_mosaic_spelling();
        test_li_stephens_prior();
        test_sequence_equivalent_prior_aggregation();
        test_conceded_counterexamples();
        std::cout << "genotype mosaic core: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "genotype mosaic core: FAIL: " << error.what() << '\n';
        return 1;
    }
}
