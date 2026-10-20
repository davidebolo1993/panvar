// M10-B diagnostic (no science, no reads, no truth): the CYP2D6 pilot's block-allele dump plus each allele's oriented
// node walks and the sorted graph's node sequences -- the model information a structural projection can reuse.
//   additionally writes <out>.walks.tsv (block_index, allele_index, walk_index, walk as node+/- steps) and
//   <out>.nodes.tsv (node, length, sequence) for every node any walk visits
// M10-B.1 (additive): <out>.blocks.tsv (the locus report's per-block divergence trims, bypass and representation flags),
//   <out>.sources.tsv (per allele: counts of walk / empty / bypass / other source records) and <out>.segments.tsv
//   (invariant segment lengths), so the boundary transformation can be audited instead of assumed
//
// Builds the model exactly as genotype_probe's load_model does (reference sort with flipping, cactus snarls, bubbles,
// locus input with the panel curation, prepare_genotype_model) and prints its fingerprint, so a caller can require it
// to equal the frozen model's. Writes:
//   <out>.block_alleles.tsv   block_index, block_id, kind, allele_index, length, sequence
//   <out>.panel.tsv           template, comma-separated allele index per block
//   <out>.summary.tsv         fingerprint, blocks, templates, canonical alleles
// Usage: m10b_walk_dump <gfa> <reference path> <curation file | none> <out prefix>
#include "panvar/bubbles.hpp"
#include "panvar/genotype_locus.hpp"
#include "panvar/genotype_model.hpp"
#include "panvar/gfa.hpp"
#include "panvar/gfa_io.hpp"
#include "panvar/graph_sort.hpp"
#include "panvar/graph_utils.hpp"
#include "panvar/integrated_snarls.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

std::vector<std::string> read_curation(const std::string& path) {
    std::vector<std::string> names;
    if (path == "none") return names;
    std::ifstream file(path);
    if (!file) throw std::runtime_error("cannot open curation " + path);
    std::string line;
    while (std::getline(file, line)) {
        const std::size_t first = line.find_first_not_of(" \t\r\n");
        if (first == std::string::npos || line[first] == '#') continue;
        const std::size_t last = line.find_last_not_of(" \t\r\n");
        names.push_back(line.substr(first, last - first + 1));
    }
    return names;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 5) {
        std::cerr << "usage: m10b_walk_dump <gfa> <reference path> <curation|none> <out prefix>\n";
        return 2;
    }
    try {
        const std::string gfa_path = argv[1], reference = argv[2], out = argv[4];
        const std::string stage = out + ".sorted.gfa";
        panvar::GfaModel gfa = panvar::read_gfa_model(gfa_path);
        panvar::GraphSortOptions sort_options;
        sort_options.reference_path = reference;
        sort_options.flip = true;
        const panvar::GraphSortResult sorted = panvar::sort_graph_reference(gfa, sort_options);
        panvar::write_gfa_model(stage, gfa);
        panvar::ParseGfaOptions parse_options;
        const panvar::Graph graph = panvar::parse_gfa(stage, parse_options);
        panvar::validate_graph_paths(graph, "cyp2d6_block_alleles", true, true);
        panvar::BubbleCallOptions bubble_options;
        bubble_options.reference_path = sorted.resolved_reference;
        bubble_options.quiet = true;
        bubble_options.min_variant_bp = 50;
        bubble_options.snarl_pairs_override =
            panvar::find_top_level_snarls_cactus(panvar::snarl_input_from_model(gfa));
        bubble_options.snarl_source_supplied = true;
        const std::vector<panvar::Bubble> bubbles = panvar::call_bubbles(graph, bubble_options);
        panvar::LocusInputOptions locus_options;
        locus_options.locus_id = sorted.resolved_reference;
        locus_options.reference_path_name = sorted.resolved_reference;
        locus_options.excluded_templates = read_curation(argv[3]);
        const panvar::LocusInputResult built = panvar::build_locus_input(graph, bubbles, locus_options);
        const panvar::PreparedGenotypeModel model = panvar::prepare_genotype_model(built.input);
        std::filesystem::remove(stage);
        {
            std::ofstream walks(out + ".walks.tsv");
            walks << "block_index\tallele_index\twalk_index\twalk\n";
            std::set<std::string> used;
            for (std::size_t b = 0; b < model.blocks().size(); ++b) {
                const panvar::PreparedGenotypeBlock& block = model.blocks()[b];
                for (std::size_t a = 0; a < block.alleles.size(); ++a) {
                    std::size_t w = 0;
                    for (const std::string& source : block.alleles[a].source_ids) {
                        const std::size_t at = source.find("walk:");
                        if (at == std::string::npos) continue;
                        const std::string steps = source.substr(at + 5);
                        walks << b << '\t' << a << '\t' << w++ << '\t' << steps << '\n';
                        std::stringstream ss(steps);
                        for (std::string step; std::getline(ss, step, ',');) used.insert(step.substr(0, step.size() - 1));
                    }
                }
            }
            std::ofstream sources(out + ".sources.tsv");
            sources << "block_index\tallele_index\twalk_sources\tempty_sources\tbypass_sources\tother_sources\n";
            for (std::size_t b = 0; b < model.blocks().size(); ++b) {
                const panvar::PreparedGenotypeBlock& block = model.blocks()[b];
                for (std::size_t a = 0; a < block.alleles.size(); ++a) {
                    std::size_t w = 0, e = 0, y = 0, o = 0;
                    for (const std::string& source : block.alleles[a].source_ids) {
                        if (source.find("walk:") != std::string::npos) ++w;
                        else if (source == "empty") ++e;
                        else if (source == "bypass") ++y;
                        else ++o;
                    }
                    sources << b << '\t' << a << '\t' << w << '\t' << e << '\t' << y << '\t' << o << '\n';
                }
            }
            std::ofstream blocks(out + ".blocks.tsv");
            blocks << "block_index\tmodel_block_id\treport_block_id\tkind\ttrimmed_prefix_bp\ttrimmed_suffix_bp\t"
                      "representation_dependent\tbypass_templates\tsource_node\tsink_node\n";
            const auto& report = built.report.blocks;
            for (std::size_t b = 0; b < model.blocks().size(); ++b) {
                blocks << b << '\t' << model.blocks()[b].block_id << '\t';
                if (b < report.size()) {
                    const panvar::LocusBlockReport& r = report[b];
                    blocks << r.block_id << '\t' << (r.kind == panvar::GenotypeBlockKind::Variant ? "SV" : "context") << '\t'
                           << r.trimmed_prefix_bp << '\t' << r.trimmed_suffix_bp << '\t' << (r.representation_dependent ? 1 : 0)
                           << '\t' << r.bypass_templates << '\t' << (r.source_node.empty() ? "-" : r.source_node) << '\t'
                           << (r.sink_node.empty() ? "-" : r.sink_node) << '\n';
                } else {
                    blocks << "MISSING\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n";
                }
            }
            std::ofstream segments(out + ".segments.tsv");
            segments << "segment_index\tlength\n";
            for (std::size_t i = 0; i < built.input.invariant_segments.size(); ++i) {
                segments << i << '\t' << built.input.invariant_segments[i].size() << '\n';
            }
            std::ofstream nodes(out + ".nodes.tsv");
            nodes << "node\tlength\tsequence\n";
            for (const std::string& id : used) {
                const auto it = graph.nodes.find(id);
                if (it == graph.nodes.end()) throw std::runtime_error("walk node missing from the graph: " + id);
                nodes << id << '\t' << it->second.sequence.size() << '\t' << it->second.sequence << '\n';
            }
        }

        std::ofstream alleles(out + ".block_alleles.tsv");
        alleles << "block_index\tblock_id\tkind\tallele_index\tlength\tsequence\n";
        std::size_t canonical = 0;
        for (std::size_t b = 0; b < model.blocks().size(); ++b) {
            const panvar::PreparedGenotypeBlock& block = model.blocks()[b];
            for (std::size_t a = 0; a < block.alleles.size(); ++a) {
                ++canonical;
                alleles << b << '\t' << block.block_id << '\t'
                        << (block.kind == panvar::GenotypeBlockKind::Variant ? "SV" : "context") << '\t' << a << '\t'
                        << block.alleles[a].sequence.size() << '\t' << block.alleles[a].sequence << '\n';
            }
        }
        std::ofstream panel(out + ".panel.tsv");
        panel << "template\talleles\n";
        for (std::size_t t = 0; t < model.panel().template_count(); ++t) {
            panel << model.panel().template_names[t] << '\t';
            for (std::size_t b = 0; b < model.panel().template_alleles[t].size(); ++b) {
                panel << (b ? "," : "") << model.panel().template_alleles[t][b];
            }
            panel << '\n';
        }
        std::ofstream summary(out + ".summary.tsv");
        summary << "key\tvalue\nmodel_fingerprint\t" << panvar::model_fingerprint(model) << "\nblocks\t"
                << model.blocks().size() << "\ntemplates\t" << model.panel().template_count()
                << "\ncanonical_alleles\t" << canonical << '\n';
        std::cout << "fingerprint " << panvar::model_fingerprint(model) << ", " << model.blocks().size() << " blocks, "
                  << canonical << " alleles\n";
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << '\n';
        return 1;
    }
    return 0;
}
