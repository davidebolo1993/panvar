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

} // namespace

int main(int argc, char** argv) {
    try {
        test_fingerprint_invariance();
        test_partial_frame_refusal();
        test_bypass_is_a_state();
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
