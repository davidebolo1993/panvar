// Phase 3 and Phase 4 both defer an acceptance gate to the graph-facing adapter. This fixture is
// those gates: it runs the real C4 panel through the adapter and proves the prepared model respells
// every panel template byte-for-byte, that the block census is what the model contract assumes, and
// that a representation change which preserves sequence does not change the model fingerprint.

#include "panvar/bubbles.hpp"
#include "panvar/genotype_locus.hpp"
#include "panvar/genotype_model.hpp"
#include "panvar/gfa.hpp"
#include "panvar/gfa_io.hpp"
#include "panvar/graph_sort.hpp"
#include "panvar/graph_utils.hpp"
#include "panvar/integrated_snarls.hpp"
#include "panvar/marker_model.hpp"
#include "panvar/mosaic_spelling.hpp"

#include <algorithm>
#include <cstdio>
#include <random>
#include <filesystem>
#include <fstream>
#include <array>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

void require_throws(const std::function<void()>& operation, const std::string& message) {
    try {
        operation();
    } catch (const std::invalid_argument&) {
        return;
    }
    throw std::runtime_error(message);
}

// ---------------------------------------------------------------------------------------------
// Synthetic gates: adapter behaviour that must hold with no fixture file present.
// ---------------------------------------------------------------------------------------------

panvar::StructuralLocusInput tiny_locus() {
    panvar::StructuralLocusInput input;
    input.locus_id = "tiny";
    input.blocks = {
        {"ctx:0", {{"walk:a", "GGG"}}, panvar::GenotypeBlockKind::Context},
        {"sv:1", {{"walk:ref", "AAAA"}, {"walk:alt", "A"}, {"bypass", ""}},
         panvar::GenotypeBlockKind::Variant},
        {"ctx:2", {{"walk:b", "TTT"}, {"walk:c", "TTA"}}, panvar::GenotypeBlockKind::Context},
    };
    input.invariant_segments = {"", "", "", ""};
    input.templates = {
        {"t0", panvar::TemplateFrameStatus::Complete, {"walk:a", "walk:ref", "walk:b"}},
        {"t1", panvar::TemplateFrameStatus::Complete, {"walk:a", "walk:alt", "walk:c"}},
        {"t2", panvar::TemplateFrameStatus::Complete, {"walk:a", "bypass", "walk:b"},
         panvar::ChainOrientation::Reverse},
    };
    return input;
}

void test_fingerprint_invariance() {
    const panvar::StructuralLocusInput input = tiny_locus();
    const std::string reference =
        panvar::model_fingerprint(panvar::prepare_genotype_model(input));

    // Renaming source walks models a node split or merge: the bytes are untouched, so the
    // statistical model, and therefore the fingerprint, must not move.
    panvar::StructuralLocusInput renamed = input;
    for (auto& block : renamed.blocks) {
        for (auto& allele : block.alleles) {
            if (allele.source_id != "bypass") allele.source_id += ":split";
        }
    }
    for (auto& tpl : renamed.templates) {
        for (auto& source : tpl.allele_source_ids) {
            if (source != "bypass") source += ":split";
        }
    }
    require(panvar::model_fingerprint(panvar::prepare_genotype_model(renamed)) == reference,
            "renaming source walks must not change the model fingerprint");

    // Reordering allele and template records must not move it either.
    panvar::StructuralLocusInput reordered = input;
    std::reverse(reordered.templates.begin(), reordered.templates.end());
    for (auto& block : reordered.blocks) {
        std::reverse(block.alleles.begin(), block.alleles.end());
    }
    require(panvar::model_fingerprint(panvar::prepare_genotype_model(reordered)) == reference,
            "record order must not change the model fingerprint");

    // Adding a sequence-identical route is also representation, not model.
    panvar::StructuralLocusInput aliased = input;
    aliased.blocks[1].alleles.push_back({"walk:ref-alias", "AAAA"});
    require(panvar::model_fingerprint(panvar::prepare_genotype_model(aliased)) == reference,
            "a sequence-identical alias must not change the model fingerprint");

    // Changing a block boundary IS a model change and must move the fingerprint.
    panvar::StructuralLocusInput reblocked = input;
    reblocked.blocks[2].alleles[1].sequence = "TTG";
    for (auto& tpl : reblocked.templates) {
        if (tpl.name == "t1") tpl.allele_source_ids[2] = "walk:c";
    }
    require(panvar::model_fingerprint(panvar::prepare_genotype_model(reblocked)) != reference,
            "changing an allele sequence must change the model fingerprint");

    // A block kind is part of the model: it decides what is reported as a call.
    panvar::StructuralLocusInput rekinded = input;
    rekinded.blocks[0].kind = panvar::GenotypeBlockKind::Variant;
    require(panvar::model_fingerprint(panvar::prepare_genotype_model(rekinded)) != reference,
            "changing a block kind must change the model fingerprint");
}

void test_partial_frame_refusal() {
    panvar::StructuralLocusInput truncated = tiny_locus();
    truncated.templates[1].frame_status = panvar::TemplateFrameStatus::Partial;
    require_throws([&] { panvar::prepare_genotype_model(truncated); },
                   "a partial template frame must refuse, not be silently completed");

    panvar::StructuralLocusInput short_frame = tiny_locus();
    short_frame.templates[1].allele_source_ids.pop_back();
    require_throws([&] { panvar::prepare_genotype_model(short_frame); },
                   "a template missing an allele at a block must refuse");
}

void test_bypass_is_a_state() {
    const panvar::PreparedGenotypeModel model =
        panvar::prepare_genotype_model(tiny_locus());
    const panvar::PreparedGenotypeBlock& sv = model.blocks()[1];
    const std::size_t bypass = sv.allele_id_for_source("bypass");
    require(sv.alleles[bypass].sequence.empty(),
            "a bypass allele is an explicit empty-sequence state");
    require(sv.alleles.size() == 3,
            "a bypass state must not collapse into a traversing allele");

    // It must spell, and it must spell nothing, so a deletion predicts zero marker observations
    // rather than being unrepresentable.
    const panvar::HaplotypeSequence spelled = panvar::spell_mosaic(
        model, panvar::MosaicPath{{0, bypass, 0}});
    require(spelled.block_intervals[1].size() == 0, "a bypass block interval must be zero-width");

    const panvar::PreparedMarkerPanel panel = panvar::build_marker_panel(model, 3, 3);
    require(panel.blocks[1].allele_markers[bypass].empty(),
            "a bypass allele must carry no markers");
}

// ---------------------------------------------------------------------------------------------
// Real C4 gates.
// ---------------------------------------------------------------------------------------------

struct PreparedLocus {
    panvar::LocusInputResult built;
    panvar::PreparedGenotypeModel model;
    panvar::Graph graph;
    std::vector<panvar::Bubble> bubbles;
};

panvar::Graph sorted_graph(panvar::GfaModel& model, const std::string& reference,
                           std::string& resolved, const std::string& stage_path) {
    panvar::GraphSortOptions sort_options;
    sort_options.reference_path = reference;
    sort_options.flip = true;
    resolved = panvar::sort_graph_reference(model, sort_options).resolved_reference;
    panvar::write_gfa_model(stage_path, model);
    panvar::ParseGfaOptions parse_options;
    return panvar::parse_gfa(stage_path, parse_options);
}

std::vector<panvar::Bubble> call_sv_bubbles(
    const panvar::Graph& graph,
    const panvar::GfaModel& model,
    const std::string& reference) {
    panvar::BubbleCallOptions options;
    options.reference_path = reference;
    options.quiet = true;
    options.min_variant_bp = 50;   // the production SV threshold: panvar genotypes structural variants
    options.snarl_pairs_override =
        panvar::find_top_level_snarls_cactus(panvar::snarl_input_from_model(model));
    options.snarl_source_supplied = true;
    return panvar::call_bubbles(graph, options);
}

PreparedLocus prepare_from_gfa(const std::string& gfa_path, const std::string& reference,
                               const std::string& stage_path,
                               const std::function<void(panvar::GfaModel&)>& rewrite) {
    panvar::GfaModel model = panvar::read_gfa_model(gfa_path);
    if (rewrite) rewrite(model);
    std::string resolved;
    panvar::Graph graph = sorted_graph(model, reference, resolved, stage_path);
    std::vector<panvar::Bubble> bubbles = call_sv_bubbles(graph, model, resolved);

    panvar::LocusInputOptions options;
    options.locus_id = "c4";
    options.reference_path_name = resolved;
    panvar::LocusInputResult built = panvar::build_locus_input(graph, bubbles, options);
    panvar::PreparedGenotypeModel prepared = panvar::prepare_genotype_model(built.input);
    return PreparedLocus{std::move(built), std::move(prepared), std::move(graph),
                         std::move(bubbles)};
}

// Splitting every node in half is the sharpest representation change that preserves sequence: every
// node identifier, every walk signature and every edge changes, while the bytes each path spells do
// not. The prepared model must be identical.
void split_every_node(panvar::GfaModel& model) {
    std::vector<std::string> order;
    std::unordered_map<std::string, std::string> seq;
    std::map<std::string, std::pair<std::string, std::string>> halves;
    for (const std::string& id : model.node_order) {
        const std::string& bases = model.seq.at(id);
        if (bases.size() < 2) {
            order.push_back(id);
            seq.emplace(id, bases);
            continue;
        }
        const std::string left = id + "a";
        const std::string right = id + "b";
        const std::size_t cut = bases.size() / 2;
        order.push_back(left);
        order.push_back(right);
        seq.emplace(left, bases.substr(0, cut));
        seq.emplace(right, bases.substr(cut));
        halves.emplace(id, std::make_pair(left, right));
    }
    std::vector<panvar::GfaEdge> edges;
    for (const auto& item : halves) {
        // Interior link, forward orientation: the two halves are read left then right.
        panvar::GfaEdge edge;
        edge.from = item.second.first;
        edge.to = item.second.second;
        edges.push_back(edge);
    }
    const auto outgoing = [&](const std::string& id, char orient) {
        const auto found = halves.find(id);
        if (found == halves.end()) return id;
        return orient == '+' ? found->second.second : found->second.first;
    };
    const auto incoming = [&](const std::string& id, char orient) {
        const auto found = halves.find(id);
        if (found == halves.end()) return id;
        return orient == '+' ? found->second.first : found->second.second;
    };
    for (const panvar::GfaEdge& edge : model.edges) {
        panvar::GfaEdge out = edge;
        out.from = outgoing(edge.from, edge.from_orient);
        out.to = incoming(edge.to, edge.to_orient);
        edges.push_back(out);
    }
    for (panvar::GfaPath& path : model.paths) {
        std::vector<panvar::PathStep> steps;
        for (const panvar::PathStep& step : path.steps) {
            const auto found = halves.find(step.node_id);
            if (found == halves.end()) {
                steps.push_back(step);
                continue;
            }
            panvar::PathStep first = step;
            panvar::PathStep second = step;
            // A reverse step reads the right half first, then the left.
            first.node_id = step.reverse ? found->second.second : found->second.first;
            second.node_id = step.reverse ? found->second.first : found->second.second;
            steps.push_back(first);
            steps.push_back(second);
        }
        path.steps = std::move(steps);
    }
    model.node_order = std::move(order);
    model.seq = std::move(seq);
    model.edges = std::move(edges);
}

void report_census(const panvar::LocusInputReport& report) {
    std::cout << "  blocks: " << report.variant_blocks << " SV + " << report.context_blocks
              << " context = " << (report.variant_blocks + report.context_blocks) << '\n'
              << "  alleles: " << report.total_source_alleles << " source walks -> "
              << report.total_canonical_alleles << " canonical sequences\n"
              << "  templates: " << report.forward_templates << " forward, "
              << report.reverse_templates << " reverse, " << report.refused_templates.size()
              << " refused\n"
              << "  contained bubbles dropped: " << report.contained_bubbles_dropped
              << ", abutting boundaries trimmed: " << report.abutting_boundaries_trimmed << '\n'
              << "  boundary bases moved into context blocks: "
              << report.boundary_bases_moved_to_context << '\n';
    std::cout << "  " << "block                      kind      alleles   min bp   max bp  bypass\n";
    for (const panvar::LocusBlockReport& block : report.blocks) {
        std::printf("  %-26s %-8s %8zu %8zu %8zu %7zu\n", block.block_id.c_str(),
                    block.kind == panvar::GenotypeBlockKind::Variant ? "SV" : "context",
                    block.canonical_alleles, block.min_bp, block.max_bp, block.bypass_templates);
    }
}

void test_c4(const std::string& gfa_path, const std::string& stage_prefix) {
    const std::string reference = "HG00096#1#haplotype1-0000024:31848049-32074186";
    const PreparedLocus locus =
        prepare_from_gfa(gfa_path, reference, stage_prefix + ".sorted.gfa", nullptr);
    report_census(locus.built.report);

    // Gate 1: every complete panel template respells byte-for-byte in its own frame.
    panvar::verify_complete_panel_round_trip(locus.model, locus.built.authoritative);
    require(locus.built.authoritative.size() == 131,
            "all 131 C4 panel paths must reach the round-trip verifier");

    // Gate 3: the census the model contract was written against.
    const panvar::LocusInputReport& report = locus.built.report;
    require(report.refused_templates.empty(), "no C4 panel template may be refused");
    require(report.forward_templates + report.reverse_templates == 131,
            "every C4 panel path must receive a chain orientation");
    require(report.reverse_templates > 0,
            "C4 carries antiparallel templates; reverse-frame spelling must be exercised");
    require(report.variant_blocks >= 1 && report.context_blocks == report.variant_blocks + 1,
            "the chain must tile: one context block between and beyond every SV block");
    for (const panvar::LocusBlockReport& block : report.blocks) {
        require(block.canonical_alleles >= 1, "every block must carry at least one allele");
    }

    // Blocks must tile the locus exactly: spelling a template from the chain reproduces its length,
    // which the round-trip above already proves byte-wise, but assert the interval frame too.
    const panvar::HaplotypeSequence spelled = panvar::spell_panel_template(locus.model, 0);
    std::size_t covered = 0;
    for (const panvar::SequenceInterval& interval : spelled.block_intervals) {
        covered += interval.size();
    }
    require(covered == spelled.bases.size(),
            "block intervals must tile the spelled haplotype with no gap or overlap");

    // Phase 4 gate 3: every canonical allele is scored, with no cap.
    const panvar::PreparedMarkerPanel panel = panvar::build_marker_panel(locus.model, 31, 15);
    require(panel.total_alleles == panel.scored_alleles,
            "every C4 block allele must be scored with no fallback");
    require(panel.total_alleles == report.total_canonical_alleles,
            "the marker panel must score exactly the canonical allele set");
    // `scored_alleles == total_alleles` holds by construction because version one has no allele cap,
    // so on its own that equality proves nothing. The load-bearing check is the one above — the panel
    // covers exactly the model's canonical alleles — plus the coverage census below, which shows the
    // markers are actually distributed over those alleles rather than concentrated in a few.
    std::size_t alleles_with_markers = 0;
    std::size_t empty_alleles = 0;
    for (std::size_t b = 0; b < panel.blocks.size(); ++b) {
        require(panel.blocks[b].allele_markers.size() == locus.model.blocks()[b].alleles.size(),
                "every block allele must have a marker record, even an empty one");
        for (std::size_t a = 0; a < panel.blocks[b].allele_markers.size(); ++a) {
            if (!panel.blocks[b].allele_markers[a].empty()) ++alleles_with_markers;
            if (locus.model.blocks()[b].alleles[a].sequence.empty()) ++empty_alleles;
        }
    }
    std::cout << "  alleles carrying >=1 discriminating marker: " << alleles_with_markers << " / "
              << panel.total_alleles << " (" << empty_alleles << " are empty deletion alleles)\n";
    require(alleles_with_markers + empty_alleles >= panel.total_alleles - panel.blocks.size(),
            "all but a handful of non-empty alleles must carry a discriminating marker");
    // Real-C4 depth-fit convergence: the actual anchor multiplicities of the C4 panel, with counts from
    // three regimes -- an interior optimum, an optimum on the background floor from noisy Poisson counts,
    // and an exact negative-intercept fit that forces the boundary. Every one must converge and say how.
    {
        std::mt19937_64 rng(2026);
        const auto fit = [&](const std::function<std::uint64_t(std::uint32_t)>& draw, const std::string& label) {
            panvar::MarkerObservations observations;
            observations.counts.assign(panel.codes.size(), 0);
            for (std::size_t k = 0; k < panel.anchor_slots.size(); ++k) {
                observations.counts[panel.anchor_slots[k]] = draw(panel.anchor_haploid_multiplicity[k]);
            }
            const panvar::MarkerDepthEstimate estimate =
                panvar::estimate_marker_depth(panel, observations, std::nullopt);
            require(estimate.converged,
                    "real-C4 depth fit must converge (" + label + "), stopped by " + estimate.termination);
            std::cout << "  depth fit on C4 anchors, " << label << ": rate " << estimate.rate_per_copy
                      << ", background " << estimate.background << ", " << estimate.termination << " in "
                      << estimate.iterations << " iterations\n";
            return estimate;
        };
        const auto interior = fit([&](std::uint32_t m) {
            return static_cast<std::uint64_t>(std::poisson_distribution<long>(1.5 + 9.32 * 2.0 * m)(rng));
        }, "interior background 1.5");
        require(std::abs(interior.rate_per_copy - 9.32) < 0.02 * 9.32 && interior.termination == "interior_stationary",
                "the interior real-C4 fit must recover the per-copy rate in the interior");
        const auto noisy_floor = fit([&](std::uint32_t m) {
            return static_cast<std::uint64_t>(std::poisson_distribution<long>(9.77 * 2.0 * m)(rng));
        }, "no background, Poisson counts");
        require(std::abs(noisy_floor.rate_per_copy - 9.77) < 0.02 * 9.77,
                "the no-background real-C4 fit must recover the per-copy rate");
        const auto boundary = fit([](std::uint32_t m) { return static_cast<std::uint64_t>(22 * m - 4); },
                                  "negative intercept");
        require(boundary.termination == "boundary_kkt" && boundary.background_at_floor &&
                    boundary.gradient_background < 0.0,
                "an exact negative-intercept fit on C4 anchors must stop on the boundary by KKT");
    }
    std::cout << "  markers: universe " << panel.codes.size() << ", anchors "
              << panel.anchor_slots.size() << ", dropped multi-block "
              << panel.dropped_multiblock_markers << '\n';
    std::size_t blocks_without_markers = 0;
    for (const panvar::BlockMarkerInventory& block : panel.blocks) {
        if (block.marker_slots.empty()) ++blocks_without_markers;
    }
    std::cout << "  blocks left with no discriminating marker: " << blocks_without_markers << " / "
              << panel.blocks.size() << '\n';

    // Gate 2: a node-split equivalent graph must produce an identical prepared model.
    const std::string reference_fingerprint = panvar::model_fingerprint(locus.model);
    const PreparedLocus split = prepare_from_gfa(
        gfa_path, reference, stage_prefix + ".split.gfa", split_every_node);
    std::cout << "  --- node-split graph ---\n";
    report_census(split.built.report);
    std::cout << "  fingerprint: " << reference_fingerprint << " (original) vs "
              << panvar::model_fingerprint(split.model) << " (node-split)\n";
    for (std::size_t b = 0; b < std::min(report.blocks.size(), split.built.report.blocks.size());
         ++b) {
        const auto& lhs = locus.model.blocks()[b];
        const auto& rhs = split.model.blocks()[b];
        if (lhs.alleles.size() != rhs.alleles.size()) {
            std::cout << "  block " << b << " allele count " << lhs.alleles.size() << " vs "
                      << rhs.alleles.size() << '\n';
            continue;
        }
        for (std::size_t a = 0; a < lhs.alleles.size(); ++a) {
            if (lhs.alleles[a].sequence != rhs.alleles[a].sequence) {
                std::cout << "  block " << b << " allele " << a << " differs: "
                          << lhs.alleles[a].sequence.size() << " vs "
                          << rhs.alleles[a].sequence.size() << " bp\n";
                break;
            }
        }
    }
    require(panvar::model_fingerprint(split.model) == reference_fingerprint,
            "splitting every graph node must not change the prepared model fingerprint");
    panvar::verify_complete_panel_round_trip(split.model, split.built.authoritative);

    // The evaluation projection frame is the model frame: every panel template projected from the
    // graph must reproduce its own model row, block by block, across both chain orientations.
    std::size_t projected = 0;
    for (std::size_t t = 0; t < locus.model.panel().template_count(); ++t) {
        const std::string& name = locus.model.panel().template_names[t];
        const auto projection = panvar::project_path_onto_locus(
            locus.graph, locus.bubbles, locus.built.report, name);
        require(projection.has_value(), "every C4 panel template must project onto the chain");
        require(projection->orientation == locus.model.template_source_orientations()[t],
                "projection must recover the template's chain orientation");
        for (std::size_t b = 0; b < locus.model.blocks().size(); ++b) {
            const std::string& allele =
                locus.model.blocks()[b].alleles[locus.model.panel().template_alleles[t][b]].sequence;
            require(projection->blocks[b].has_value() && *projection->blocks[b] == allele,
                    "projection of " + name + " must equal its model row at " +
                        locus.model.blocks()[b].block_id);
        }
        ++projected;
    }
    require(projected == 131, "all 131 C4 templates must be projected and checked");

    // Leave-one-out frame: NA19240 is removed from every proposal source before the model is built,
    // then projected from the full graph onto that model's blocks. Sorting depends only on the
    // reference path, so node identifiers name the same sequence in both graphs.
    const PreparedLocus loo = prepare_from_gfa(
        gfa_path, reference, stage_prefix + ".loo.gfa", [](panvar::GfaModel& gfa) {
            gfa.paths.erase(std::remove_if(gfa.paths.begin(), gfa.paths.end(),
                                           [](const panvar::GfaPath& path) {
                                               return panvar::gfa_path_name(path).rfind("NA19240#", 0) == 0;
                                           }),
                            gfa.paths.end());
        });
    require(loo.model.panel().template_count() == 129,
            "leave-one-out must remove exactly NA19240's two haplotypes");
    for (const panvar::PathRecord& path : locus.graph.paths) {
        if (path.name.rfind("NA19240#", 0) != 0) continue;
        const auto projection =
            panvar::project_path_onto_locus(locus.graph, loo.bubbles, loo.built.report, path.name);
        require(projection.has_value(), "the held-out donor must span the leave-one-out chain");
        std::size_t defined = 0;
        std::size_t available = 0;
        for (std::size_t b = 0; b < loo.model.blocks().size(); ++b) {
            if (!projection->blocks[b].has_value()) continue;
            ++defined;
            for (const auto& allele : loo.model.blocks()[b].alleles) {
                if (allele.sequence == *projection->blocks[b]) {
                    ++available;
                    break;
                }
            }
        }
        std::cout << "  leave-one-out " << path.name << ": " << defined << "/"
                  << loo.model.blocks().size() << " blocks defined, " << available
                  << " exactly available in the retained panel\n";
        require(defined == loo.model.blocks().size(),
                "every block of the held-out donor must have a positional frame at C4");
    }
}

// Fail closed on representation dependence. Template C deletes bubble b1 together with both
// neighbouring backbones, joining b0's sink directly to b2's source. Its spelling still round-trips,
// so the round-trip gate does not catch this: the traversing alleles of b1 share boundary-node bytes
// C does not have, the divergence trim cannot move them, and where the graph cut those boundary nodes
// would decide what sits in b1 and both its neighbours.
void test_bypass_representation_dependence(const std::string& stage_prefix) {
    std::mt19937_64 rng(88);
    const auto bases = [&](std::size_t length) {
        std::string out(length, 'A');
        for (char& base : out) base = "ACGT"[rng() % 4];
        return out;
    };
    // Numeric identifiers in the reference template's order, so node order is reference order.
    const std::vector<std::pair<std::string, std::string>> nodes = {
        {"1", bases(300)}, {"2", bases(40)}, {"3", bases(120)}, {"4", bases(70)},
        {"5", bases(40)}, {"6", bases(200)}, {"7", bases(40)}, {"8", bases(130)},
        {"9", bases(60)}, {"10", bases(40)}, {"11", bases(200)}, {"12", bases(40)},
        {"13", bases(110)}, {"14", bases(55)}, {"15", bases(40)}, {"16", bases(300)},
    };
    const std::vector<std::pair<std::string, std::vector<std::string>>> all_templates = {
        {"A#1#x", {"1", "2", "3", "5", "6", "7", "8", "10", "11", "12", "13", "15", "16"}},
        {"B#1#x", {"1", "2", "4", "5", "6", "7", "9", "10", "11", "12", "14", "15", "16"}},
        {"C#1#x", {"1", "2", "3", "5", "12", "13", "15", "16"}},
    };
    const auto write_gfa = [&](const std::string& file, bool with_bypass) {
        std::ofstream out(file);
        out << "H\tVN:Z:1.0\n";
        for (const auto& node : nodes) out << "S\t" << node.first << '\t' << node.second << '\n';
        std::set<std::pair<std::string, std::string>> edges;
        for (const auto& item : all_templates) {
            if (!with_bypass && item.first == "C#1#x") continue;
            for (std::size_t i = 1; i < item.second.size(); ++i) {
                edges.emplace(item.second[i - 1], item.second[i]);
            }
        }
        for (const auto& edge : edges) {
            out << "L\t" << edge.first << "\t+\t" << edge.second << "\t+\t0M\n";
        }
        for (const auto& item : all_templates) {
            if (!with_bypass && item.first == "C#1#x") continue;
            out << "P\t" << item.first << '\t';
            for (std::size_t i = 0; i < item.second.size(); ++i) {
                out << (i ? "," : "") << item.second[i] << '+';
            }
            out << "\t*\n";
        }
    };
    std::vector<panvar::Bubble> bubbles(3);
    const std::vector<std::array<std::string, 4>> spec = {
        {"2", "5", "3", "4"}, {"7", "10", "8", "9"}, {"12", "15", "13", "14"}};
    for (std::size_t i = 0; i < 3; ++i) {
        bubbles[i].id = i + 1;
        bubbles[i].source = spec[i][0];
        bubbles[i].sink = spec[i][1];
        bubbles[i].inside = {spec[i][2], spec[i][3]};
    }
    const auto build = [&](bool with_bypass) {
        const std::string file = stage_prefix + (with_bypass ? ".bypass.gfa" : ".nobypass.gfa");
        write_gfa(file, with_bypass);
        panvar::ParseGfaOptions parse_options;
        const panvar::Graph graph = panvar::parse_gfa(file, parse_options);
        panvar::LocusInputOptions options;
        options.locus_id = "bypass-fixture";
        options.reference_path_name = "A#1#x";
        return panvar::build_locus_input(graph, bubbles, options);
    };

    const panvar::LocusInputResult dependent = build(true);
    const std::vector<std::string> expected = {"ctx:2", "sv:3", "ctx:4"};
    require(dependent.report.representation_dependent_blocks == expected,
            "a bypassed variant block and both its neighbours must be flagged representation-dependent");
    require(dependent.report.blocks[3].bypass_templates == 1,
            "fixture: exactly one template must bypass the middle variant block");
    // The round trip passes: this is the case the byte-for-byte gate cannot see.
    const panvar::PreparedGenotypeModel model = panvar::prepare_genotype_model(dependent.input);
    panvar::verify_complete_panel_round_trip(model, dependent.authoritative);
    bool refused = false;
    try {
        panvar::require_representation_independent(dependent.report);
    } catch (const std::invalid_argument& error) {
        refused = std::string(error.what()).find("sv:3") != std::string::npos;
    }
    require(refused, "a representation-dependent model must be refused, naming its blocks");

    const panvar::LocusInputResult control = build(false);
    require(control.report.representation_dependent_blocks.empty(),
            "without the bypassing template no block may be flagged");
    panvar::require_representation_independent(control.report);
    require(control.report.blocks[3].trimmed_prefix_bp > 0 && control.report.blocks[3].trimmed_suffix_bp > 0,
            "without the bypass the middle block must be divergence-trimmed on both sides");
}

// Locus coverage. A template that does not reach a terminal anchor never observed the blocks out
// there: it has not deleted them. Recording that as a bypass invents a zero-length deletion allele
// and silently drops the bases the path does carry, which is what cost cyp2d6 1,978 bases. These
// gates fix the distinction in both directions -- uncovered must refuse, bracketed must not.
void test_locus_coverage_partial() {
    std::mt19937_64 rng(4242);
    const auto bases = [&](std::size_t length) {
        std::string out(length, 'A');
        for (char& base : out) base = "ACGT"[rng() % 4];
        return out;
    };
    const std::vector<std::pair<std::string, std::string>> nodes = {
        {"1", bases(300)}, {"2", bases(40)}, {"3", bases(120)}, {"4", bases(70)},
        {"5", bases(40)}, {"6", bases(200)}, {"7", bases(40)}, {"8", bases(130)},
        {"9", bases(60)}, {"10", bases(40)}, {"11", bases(200)}, {"12", bases(40)},
        {"13", bases(110)}, {"14", bases(55)}, {"15", bases(40)}, {"16", bases(300)},
    };
    // Chain: ctx:0 | sv:1 (2..5) | ctx:2 | sv:3 (7..10) | ctx:4 | sv:5 (12..15) | ctx:6
    const std::vector<std::string> full_a =
        {"1", "2", "3", "5", "6", "7", "8", "10", "11", "12", "13", "15", "16"};
    const std::vector<std::string> full_b =
        {"1", "2", "4", "5", "6", "7", "9", "10", "11", "12", "14", "15", "16"};
    const std::vector<std::string> bracketed_deletion = {"1", "2", "3", "5", "12", "13", "15", "16"};
    const std::vector<std::string> empty_interior =
        {"1", "2", "3", "5", "7", "8", "10", "11", "12", "13", "15", "16"};
    const std::vector<std::string> no_left_terminal =
        {"6", "7", "8", "10", "11", "12", "13", "15", "16"};
    const std::vector<std::string> no_right_terminal = {"1", "2", "3", "5", "6", "7", "8", "10", "11"};

    struct FixtureTemplate {
        std::string name;
        std::vector<std::string> nodes;
        bool reverse = false;
    };
    const auto build = [&](const std::vector<FixtureTemplate>& templates, const std::string& tag,
                           const std::vector<std::string>& curated = {}) {
        const std::string file =
            (std::filesystem::temp_directory_path() / ("genotype_coverage_" + tag + ".gfa")).string();
        {
            std::ofstream out(file);
            out << "H\tVN:Z:1.0\n";
            for (const auto& node : nodes) out << "S\t" << node.first << '\t' << node.second << '\n';
            std::set<std::pair<std::string, std::string>> edges;
            for (const auto& item : templates) {
                for (std::size_t i = 1; i < item.nodes.size(); ++i) {
                    edges.emplace(item.nodes[i - 1], item.nodes[i]);
                }
            }
            for (const auto& edge : edges) {
                out << "L\t" << edge.first << "\t+\t" << edge.second << "\t+\t0M\n";
            }
            for (const auto& item : templates) {
                out << "P\t" << item.name << '\t';
                for (std::size_t i = 0; i < item.nodes.size(); ++i) {
                    const std::string& node =
                        item.reverse ? item.nodes[item.nodes.size() - 1 - i] : item.nodes[i];
                    out << (i ? "," : "") << node << (item.reverse ? '-' : '+');
                }
                out << "\t*\n";
            }
        }
        std::vector<panvar::Bubble> bubbles(3);
        const std::vector<std::array<std::string, 4>> spec = {
            {"2", "5", "3", "4"}, {"7", "10", "8", "9"}, {"12", "15", "13", "14"}};
        for (std::size_t i = 0; i < 3; ++i) {
            bubbles[i].id = i + 1;
            bubbles[i].source = spec[i][0];
            bubbles[i].sink = spec[i][1];
            bubbles[i].inside = {spec[i][2], spec[i][3]};
        }
        panvar::ParseGfaOptions parse_options;
        const panvar::Graph graph = panvar::parse_gfa(file, parse_options);
        panvar::LocusInputOptions options;
        options.locus_id = "coverage-fixture";
        options.reference_path_name = "A#1#x";
        options.excluded_templates = curated;
        panvar::LocusInputResult built = panvar::build_locus_input(graph, bubbles, options);
        std::filesystem::remove(file);
        return built;
    };
    const auto refusal_message = [](const panvar::StructuralLocusInput& input) {
        try {
            panvar::prepare_genotype_model(input);
        } catch (const std::invalid_argument& error) {
            return std::string(error.what());
        }
        return std::string();
    };
    const auto allele_census = [](const panvar::LocusInputResult& built) {
        std::string census;
        for (const auto& block : built.input.blocks) {
            census += block.block_id + "{";
            std::set<std::string> alleles;
            for (const auto& allele : block.alleles) alleles.insert(allele.source_id + "=" + allele.sequence);
            for (const auto& allele : alleles) census += allele + ";";
            census += "}";
        }
        return census;
    };
    const auto uncovered_of = [](const panvar::LocusInputResult& built, const std::string& name) {
        for (const auto& partial : built.report.partial_templates) {
            if (partial.name == name) return partial;
        }
        return panvar::PartialTemplateReport{};
    };

    // 1. A forward template that does not reach the chain's left terminal anchor.
    const panvar::LocusInputResult left_gap = build(
        {{"A#1#x", full_a, false}, {"B#1#x", full_b, false}, {"E#1#x", no_left_terminal, false}},
        "left");
    const panvar::PartialTemplateReport left = uncovered_of(left_gap, "E#1#x");
    require(left.name == "E#1#x", "a template missing the left terminal must be reported partial");
    require(left.missing_terminal_anchor == "left", "the missing terminal anchor must be named");
    require(left.uncovered_blocks == std::vector<std::string>({"ctx:0", "sv:1", "ctx:2"}),
            "the uncovered blocks must be named exactly, in chain order");
    require(refusal_message(left_gap.input).find("locus_coverage_partial") != std::string::npos,
            "the refusal must be locus_coverage_partial, not a generic partial frame");
    for (const auto& block : left_gap.report.blocks) {
        require(block.bypass_templates == 0,
                "an uncovered template must contribute no bypass count at any block");
    }
    require(left_gap.report.representation_dependent_blocks.empty(),
            "an uncovered template must not make any block representation-dependent");

    // 2. The same gap on a reverse-oriented template: cyp2d6's NA18989 in miniature.
    const panvar::LocusInputResult reverse_gap = build(
        {{"A#1#x", full_a, false}, {"B#1#x", full_b, false}, {"E#1#x", no_left_terminal, true}},
        "reverse");
    require(reverse_gap.report.reverse_templates == 1, "fixture: the gapped template must be reverse");
    require(uncovered_of(reverse_gap, "E#1#x").missing_terminal_anchor == "left",
            "a reverse template's missing terminal anchor is still the chain's left end");
    require(uncovered_of(reverse_gap, "E#1#x").uncovered_blocks ==
                std::vector<std::string>({"ctx:0", "sv:1", "ctx:2"}),
            "orientation must not change which blocks are uncovered");

    // 3. A template missing the trailing end.
    const panvar::LocusInputResult right_gap = build(
        {{"A#1#x", full_a, false}, {"B#1#x", full_b, false}, {"F#1#x", no_right_terminal, false}},
        "right");
    require(uncovered_of(right_gap, "F#1#x").missing_terminal_anchor == "right",
            "a template stopping before the last bubble is uncovered on the right");
    require(uncovered_of(right_gap, "F#1#x").uncovered_blocks ==
                std::vector<std::string>({"ctx:4", "sv:5", "ctx:6"}),
            "the uncovered trailing blocks must be named exactly");

    // 3b. Curating the uncovered template away. Round 13 shipped this path with external cyp2d6 as its
    //     only gate -- a whole pipeline away from the adapter that implements it, and unavailable to
    //     anyone without that graph. The contract belongs here, where the code is.
    const panvar::LocusInputResult curated = build(
        {{"A#1#x", full_a, false}, {"B#1#x", full_b, false}, {"F#1#x", no_right_terminal, false}},
        "curated", {"F#1#x"});

    // It disappears from the panel.
    require(curated.input.templates.size() == 2, "a curated panel keeps exactly the other templates");
    for (const auto& item : curated.input.templates) {
        require(item.name != "F#1#x", "a curated template must not remain a panel template");
    }

    // It is recorded once, as curated and NOT as refused. Both headings would say the panel lost two
    // templates, and worse, that a partial frame reached a model that never saw it.
    require(curated.report.curated_exclusions.size() == 1 &&
                curated.report.curated_exclusions.front().template_name == "F#1#x",
            "the curated template must be recorded exactly once, by name");
    require(curated.report.refused_templates.empty(),
            "a curated template must not also be counted refused: one template, one heading");
    require(uncovered_of(curated, "F#1#x").name.empty(),
            "a curated template must leave no partial-template report behind either");

    // Its reason survives, so the record can be audited instead of trusted.
    const std::string reason = curated.report.curated_exclusions.front().reason;
    require(reason.find("locus_coverage_partial") != std::string::npos,
            "the curated record must keep the locus_coverage_partial reason");
    require(reason.find("right") != std::string::npos,
            "the curated record must keep which terminal anchor the template failed to reach");
    require(!curated.report.curation_digest.empty(), "a curated panel must carry a curation digest");

    // The remaining model builds, where the same panel uncurated refuses.
    require(refusal_message(right_gap.input).find("locus_coverage_partial") != std::string::npos,
            "fixture: uncurated, this same panel must refuse");
    require(refusal_message(curated.input).empty(),
            "the curated panel must build where the uncurated one refused");

    // No allele came from it. Checked as equivalence rather than by name -- allele source ids are
    // block-local, not template names, so searching for the name would pass whatever curation did.
    // A curated panel must be indistinguishable from one the template never entered.
    const panvar::LocusInputResult absent =
        build({{"A#1#x", full_a, false}, {"B#1#x", full_b, false}}, "absent");
    require(allele_census(curated) == allele_census(absent),
            "a curated panel must be identical to one where the template was never present");
    require(curated.report.blocks.size() == absent.report.blocks.size(),
            "curation must not change the block chain");

    // A manifest that has drifted from its graph is refused, never quietly ignored.
    bool drifted = false;
    try {
        build({{"A#1#x", full_a, false}, {"B#1#x", full_b, false}}, "drift", {"Z#1#x"});
    } catch (const std::invalid_argument& error) {
        drifted = std::string(error.what()).find("Z#1#x") != std::string::npos;
    }
    require(drifted, "curation naming a template the graph does not have must refuse, and name it");

    // 4. A genuine internal deletion, bracketed on both sides, stays a valid bypass. No real locus
    //    of the six exercises this, and the branch must survive anyway: future loci will.
    const panvar::LocusInputResult deletion = build(
        {{"A#1#x", full_a, false}, {"B#1#x", full_b, false}, {"C#1#x", bracketed_deletion, false}},
        "deletion");
    require(deletion.report.partial_templates.empty(),
            "a deletion bracketed on both sides must NOT be refused as uncovered");
    require(deletion.report.blocks[3].bypass_templates == 1,
            "a bracketed deletion must still be recorded as a bypass");
    require(deletion.report.representation_dependent_blocks ==
                std::vector<std::string>({"ctx:2", "sv:3", "ctx:4"}),
            "a bracketed bypass must still flag its block and both neighbours");
    const panvar::PreparedGenotypeModel deletion_model =
        panvar::prepare_genotype_model(deletion.input);
    panvar::verify_complete_panel_round_trip(deletion_model, deletion.authoritative);

    // 5. A traversing empty allele is not a bypass: the template walks both boundaries and finds
    //    nothing between them. gstm1 has 179 of these and myom2 286; they must stay callable.
    const panvar::LocusInputResult empty = build(
        {{"A#1#x", full_a, false}, {"B#1#x", full_b, false}, {"G#1#x", empty_interior, false}},
        "empty");
    require(empty.report.partial_templates.empty(),
            "a traversing empty allele must not be refused");
    require(empty.report.blocks[2].bypass_templates == 0,
            "a traversing empty allele must not be counted as a bypass");
    // The state is what must survive, not the byte count: the neighbouring variant block's
    // divergence trim moves its common flanks into this one, so an empty interior can legitimately
    // end up carrying trimmed bases. What must never happen is its collapsing into a bypass.
    bool has_empty_allele = false;
    bool has_bypass_allele = false;
    for (const auto& allele : empty.input.blocks[2].alleles) {
        if (allele.source_id == "empty") has_empty_allele = true;
        if (allele.source_id == "bypass") has_bypass_allele = true;
    }
    require(has_empty_allele, "an empty interior must be its own explicit allele state");
    require(!has_bypass_allele, "an empty interior must not be recorded as a bypass");
    panvar::prepare_genotype_model(empty.input);

    // 6. The mutation gate. If uncovered were folded back into bypass, the uncovered template would
    //    contribute an allele and widen three blocks. Adding it must change nothing at all.
    const panvar::LocusInputResult without =
        build({{"A#1#x", full_a, false}, {"B#1#x", full_b, false}}, "without");
    require(allele_census(left_gap) == allele_census(without),
            "an uncovered template must not add or widen a single allele anywhere");
    for (const auto& tpl : left_gap.input.templates) {
        if (tpl.name != "E#1#x") continue;
        require(tpl.frame_status == panvar::TemplateFrameStatus::Partial,
                "an uncovered template's frame must be partial");
        require(tpl.allele_source_ids.empty(), "an uncovered template must emit no allele record");
    }
}

} // namespace

int main(int argc, char** argv) {
    try {
        test_fingerprint_invariance();
        test_partial_frame_refusal();
        test_bypass_is_a_state();
        test_locus_coverage_partial();
        if (argc > 2) {
            std::cout << "genotype locus adapter: C4 gates\n";
            test_c4(argv[1], argv[2]);
            test_bypass_representation_dependence(argv[2]);
        } else {
            std::cout << "genotype locus adapter: C4 gates SKIPPED (no fixture path given)\n";
        }
        std::cout << "genotype locus adapter: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "genotype locus adapter: FAIL: " << error.what() << '\n';
        return 1;
    }
}
