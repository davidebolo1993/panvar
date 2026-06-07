#include "panvar/allele_command.hpp"

#include "panvar/allele.hpp"
#include "panvar/cli_utils.hpp"
#include "panvar/gfa.hpp"

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace panvar {
namespace {

void print_allele_help() {
    std::cout
        << "Usage:\n"
        << "  panvar allele -i <graph.gfa> --bubble-prefix-in <module1-prefix> [options]\n\n"
        << "Options:\n"
        << "  -i, --gfa <path>                 Input GFA file (required)\n"
        << "  -o, --out-prefix <prefix>        Output prefix (default: allele_calls)\n"
        << "      --bubble-prefix-in <prefix>  Module-1 output prefix from 'panvar bubble'\n"
        << "                                   (auto-uses <prefix>.bubbles.csv)\n"
        << "      --clusters-csv <path>        Explicit allele-cluster CSV output path\n"
        << "      --assignments-csv <path>     Explicit per-path assignment CSV output path\n"
        << "      --clusters-json <path>       Optional predefined path->cluster-label JSON map\n"
        << "      --cluster-sequences-csv <path> Optional representative sequence export\n"
        << "                                   (one representative sequence row per bubble+cluster)\n"
        << "      --min-similarity <X>         Min similarity in [0,1] or percent (e.g. 0.9 or 90)\n"
        << "                                    (default: 0.90)\n"
        << "      --cluster-mode <mode>        similarity mode: walk|sequence (default: walk)\n"
        << "      --threads <N>                Worker threads for distance computation (0=auto)\n"
        << "      --distance-mode <mode>       auto|exact (default: auto)\n"
        << "      --max-upgma-alleles <N>      auto-switch to threshold-graph clustering above N unique\n"
        << "                                    alleles per bubble (default: 256, 0=disable)\n"
        << "      --skip-no-reference-bubbles   Skip bubbles where --reference-path has no assignment\n"
        << "                                    (useful when preparing call-ready clusters)\n"
        << "      --quiet                      Disable progress logs\n"
        << "      --similarity-out-dir <dir>   Optional per-bubble similarity diagnostics folder\n"
        << "                                    (distance stats, pairwise cluster matrix,\n"
        << "                                     clustering quality, dendrogram SVG, UPGMA Newick)\n"
        << "\n"
        << "ODGI export options:\n"
        << "      --odgi-viz-out-dir <dir>     Emit per-bubble ODGI viz inputs in this folder\n"
        << "      --reference-path <path>      Reference path used for odgi viz -r windows\n"
        << "      --flank-nodes <N>            Add N nodes on each side of the bubble window (default: 1)\n"
        << "      --max-paths-per-bubble <N>   Cap exported paths per bubble after cluster ordering\n"
        << "                                    (default: 0 = all)\n"
        << "      --odgi-input <path>          Input graph for odgi viz (default: --gfa)\n"
        << "      --odgi-bin <path>            odgi executable (default: odgi)\n"
        << "      --odgi-width <N>             odgi viz width (default: 2400)\n"
        << "      --odgi-path-height <N>       odgi path height (default: 8)\n"
        << "      --run-odgi                   Execute odgi viz scripts while exporting\n"
        << "  -h, --help                       Show this help\n";
}

ClusterMode parse_cluster_mode_arg(const std::string& value) {
    if (value == "sequence") {
        return ClusterMode::Sequence;
    }
    if (value == "sequence-fast" || value == "seq-fast" || value == "fast-sequence") {
        return ClusterMode::SequenceFast;
    }
    if (value == "walk" || value == "walk-signature" || value == "signature") {
        return ClusterMode::Walk;
    }
    throw std::runtime_error(
        "Invalid --cluster-mode: " + value +
        " (expected sequence|walk)");
}

const char* cluster_mode_label(ClusterMode mode) {
    switch (mode) {
        case ClusterMode::Sequence:
            return "sequence";
        case ClusterMode::SequenceFast:
            return "sequence";
        case ClusterMode::Walk:
            return "walk";
    }
    return "unknown";
}

} // namespace

int run_allele_command(const std::vector<std::string>& args) {
    std::string gfa_path;
    std::string out_prefix = "allele_calls";
    std::string bubble_prefix_in;
    std::string clusters_csv_path;
    std::string assignments_csv_path;

    AlleleCallOptions options;

    for (std::size_t i = 0; i < args.size(); ++i) {
        const std::string& arg = args[i];
        auto require_value = [&](const std::string& flag) -> const std::string& {
            if (i + 1 >= args.size()) {
                throw std::runtime_error("Missing value after " + flag);
            }
            return args[++i];
        };

        if (arg == "-h" || arg == "--help") {
            print_allele_help();
            return 0;
        }
        if (arg == "-i" || arg == "--gfa") {
            gfa_path = require_value(arg);
            continue;
        }
        if (arg == "-o" || arg == "--out-prefix") {
            out_prefix = require_value(arg);
            continue;
        }
        if (arg == "--bubble-prefix-in") {
            bubble_prefix_in = require_value(arg);
            continue;
        }
        if (arg == "--clusters-csv") {
            clusters_csv_path = require_value(arg);
            continue;
        }
        if (arg == "--assignments-csv") {
            assignments_csv_path = require_value(arg);
            continue;
        }
        if (arg == "--clusters-json") {
            options.predefined_clusters_json_path = require_value(arg);
            continue;
        }
        if (arg == "--cluster-sequences-csv") {
            options.cluster_sequences_csv_path = require_value(arg);
            options.write_cluster_sequences = true;
            continue;
        }
        if (arg == "--min-similarity") {
            options.min_similarity = cli::parse_similarity_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--cluster-mode") {
            options.cluster_mode = parse_cluster_mode_arg(require_value(arg));
            continue;
        }
        if (arg == "--threads") {
            options.threads = cli::parse_size_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--distance-mode") {
            const std::string mode = require_value(arg);
            if (mode == "auto") {
                options.fast_distance = true;
            } else if (mode == "exact") {
                options.fast_distance = false;
            } else {
                throw std::runtime_error("Invalid --distance-mode: " + mode + " (expected auto|exact)");
            }
            continue;
        }
        if (arg == "--max-upgma-alleles") {
            options.max_upgma_alleles = cli::parse_size_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--skip-no-reference-bubbles") {
            options.skip_bubbles_without_reference = true;
            continue;
        }
        if (arg == "--quiet") {
            options.show_progress = false;
            continue;
        }
        if (arg == "--reference-path") {
            options.reference_path = require_value(arg);
            continue;
        }
        if (arg == "--similarity-out-dir") {
            options.similarity_out_dir = require_value(arg);
            options.write_similarity_reports = true;
            continue;
        }
        if (arg == "--odgi-viz-out-dir") {
            options.odgi_viz_out_dir = require_value(arg);
            options.write_odgi_viz_inputs = true;
            continue;
        }
        if (arg == "--flank-nodes") {
            options.flank_nodes = cli::parse_size_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--max-paths-per-bubble") {
            options.max_paths_per_bubble = cli::parse_size_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--odgi-input") {
            options.odgi_input_path = require_value(arg);
            continue;
        }
        if (arg == "--odgi-bin") {
            options.odgi_bin = require_value(arg);
            continue;
        }
        if (arg == "--odgi-width") {
            options.odgi_width = cli::parse_size_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--odgi-path-height") {
            options.odgi_path_height = cli::parse_size_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--run-odgi") {
            options.run_odgi = true;
            options.write_odgi_viz_inputs = true;
            continue;
        }

        throw std::runtime_error("Unknown option: " + arg);
    }

    if (gfa_path.empty()) {
        throw std::runtime_error("Missing required input: --gfa <path>");
    }
    if (bubble_prefix_in.empty()) {
        throw std::runtime_error(
            "Module 'allele' requires module-1 bubbles prefix: use --bubble-prefix-in <prefix>.");
    }
    options.bubbles_csv_in = bubble_prefix_in + ".bubbles.csv";
    if (options.skip_bubbles_without_reference && options.reference_path.empty()) {
        throw std::runtime_error(
            "--skip-no-reference-bubbles requires --reference-path in module 'allele'.");
    }

    if (options.write_odgi_viz_inputs && options.odgi_viz_out_dir.empty()) {
        options.odgi_viz_out_dir = out_prefix + ".odgi_viz";
    }
    if (options.write_odgi_viz_inputs && options.odgi_input_path.empty()) {
        options.odgi_input_path = gfa_path;
    }
    if (options.run_odgi && options.odgi_bin.empty()) {
        throw std::runtime_error("--odgi-bin cannot be empty when --run-odgi is set");
    }
    if (options.write_cluster_sequences && options.cluster_sequences_csv_path.empty()) {
        options.cluster_sequences_csv_path = out_prefix + ".cluster_sequences.csv";
    }
    if (options.cluster_mode == ClusterMode::SequenceFast) {
        options.fast_distance = true;
    }

    if (clusters_csv_path.empty()) {
        clusters_csv_path = out_prefix + ".allele_clusters.csv";
    }
    if (assignments_csv_path.empty()) {
        assignments_csv_path = out_prefix + ".allele_assignments.csv";
    }
    cli::ensure_parent_dir_for_file(clusters_csv_path);
    cli::ensure_parent_dir_for_file(assignments_csv_path);
    if (options.write_cluster_sequences) {
        cli::ensure_parent_dir_for_file(options.cluster_sequences_csv_path);
    }

    ParseGfaOptions parse_options;
    parse_options.include_paths = true;
    parse_options.include_sequences = true;

    const Graph graph = parse_gfa(gfa_path, parse_options);
    if (graph.paths.empty()) {
        throw std::runtime_error("Input GFA has no P/W paths; allele clustering requires paths");
    }
    AlleleCallSummary summary;
    call_alleles_to_csv(graph, options, clusters_csv_path, assignments_csv_path, &summary);

    std::cout
        << "Input graph: " << gfa_path << '\n'
        << "Bubble source: " << options.bubbles_csv_in << "\n"
        << "Nodes: " << graph.nodes.size() << "\n"
        << "P/W paths loaded: " << graph.paths.size() << "\n"
        << "Min similarity: " << options.min_similarity << "\n"
        << "Cluster mode: " << cluster_mode_label(options.cluster_mode) << "\n"
        << "Distance mode: " << (options.fast_distance ? "auto" : "exact") << "\n"
        << "Skip no-reference bubbles: " << (options.skip_bubbles_without_reference ? "yes" : "no") << "\n"
        << "Predefined clusters: "
        << (options.predefined_clusters_json_path.empty() ? "off" : options.predefined_clusters_json_path)
        << "\n"
        << "Max UPGMA alleles: " << options.max_upgma_alleles << " (0=disabled)\n"
        << "Threads: " << options.threads << " (0=auto)\n"
        << "Bubbles processed: " << summary.bubbles_processed << "\n"
        << "Bubbles with haplotype assignments: " << summary.bubbles_with_assignments << "\n"
        << "Bubbles skipped (no reference assignment): " << summary.bubbles_skipped_no_reference << "\n"
        << "Unique alleles: " << summary.unique_alleles << "\n"
        << "Allele clusters: " << summary.clusters << "\n"
        << "Path assignments: " << summary.assignments << "\n"
        << "Wrote: " << clusters_csv_path << "\n"
        << "Wrote: " << assignments_csv_path << "\n";
    if (options.write_cluster_sequences) {
        std::cout
            << "Representative cluster sequences: " << options.cluster_sequences_csv_path << "\n"
            << "Cluster sequence rows: " << summary.cluster_sequences_written << "\n";
    }
    if (options.write_similarity_reports) {
        std::cout
            << "Similarity diagnostics dir: " << options.similarity_out_dir << "\n"
            << "Similarity reports (bubbles): " << summary.similarity_reports_written << "\n";
    }
    if (options.write_odgi_viz_inputs) {
        std::cout
            << "ODGI export dir: " << options.odgi_viz_out_dir << "\n"
            << "ODGI exports (bubbles): " << summary.bubbles_with_odgi_exports << "\n"
            << "ODGI runs ok/failed: " << summary.odgi_runs_ok << "/" << summary.odgi_runs_failed << "\n";
    }

    return 0;
}

} // namespace panvar
