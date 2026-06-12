#include "panvar/bubble_command.hpp"

#include "panvar/bubbles.hpp"
#include "panvar/cli_utils.hpp"
#include "panvar/gfa.hpp"
#include "panvar/output.hpp"

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace panvar {
namespace {

void print_bubble_help() {
    std::cout
        << "Usage:\n"
        << "  panvar bubble -i <graph.gfa> --snarls-in <snarls.jsonl> [options]\n\n"
        << "Options:\n"
        << "  -i, --gfa <path>                 Input GFA file (required)\n"
        << "      --snarls-in <path>           Snarl JSONL from 'vg view -R -j' (required)\n"
        << "  -o, --out-prefix <prefix>        Output prefix (default: bubble_calls)\n"
        << "      --bubbles-csv <path>         Explicit bubbles CSV output path\n"
        << "      --bandage-csv <path>         Explicit Bandage color CSV output path\n"
        << "      --snarl-debug-tsv <path>     Optional diagnostics TSV for snarl candidates\n"
        << "      --min-variant-bp <N>         Keep bubbles with at least one path carrying >= N bp\n"
        << "                                    inside the bubble (default: 50, 0=disable)\n"
        << "      --min-path-support <N>       Require at least N supporting P/W paths (default: 0)\n"
        << "      --merge-nearby-bp <N>        Merge nearby bubbles only after base filters\n"
        << "                                    (min-path/min-variant) when sink->source shortest-path\n"
        << "                                    distance is <= N bp (default: 0, disabled)\n"
        << "      --quiet                      Disable the progress bar\n"
        << "  -h, --help                       Show this help\n";
}

} // namespace

int run_bubble_command(const std::vector<std::string>& args) {
    std::string gfa_path;
    std::string out_prefix = "bubble_calls";
    std::string bubbles_csv_path;
    std::string bandage_csv_path;
    std::string snarl_debug_tsv_path;

    BubbleCallOptions options;

    for (std::size_t i = 0; i < args.size(); ++i) {
        const std::string& arg = args[i];
        auto require_value = [&](const std::string& flag) -> const std::string& {
            if (i + 1 >= args.size()) {
                throw std::runtime_error("Missing value after " + flag);
            }
            return args[++i];
        };

        if (arg == "-h" || arg == "--help") {
            print_bubble_help();
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
        if (arg == "--bubbles-csv") {
            bubbles_csv_path = require_value(arg);
            continue;
        }
        if (arg == "--bandage-csv") {
            bandage_csv_path = require_value(arg);
            continue;
        }
        if (arg == "--snarls-in") {
            options.snarls_input_path = require_value(arg);
            continue;
        }
        if (arg == "--snarl-debug-tsv") {
            snarl_debug_tsv_path = require_value(arg);
            options.collect_snarl_debug = true;
            continue;
        }
        if (arg == "--min-variant-bp") {
            options.min_variant_bp = cli::parse_size_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--min-path-support") {
            options.min_path_support = cli::parse_size_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--merge-nearby-bp") {
            options.merge_nearby_bp = cli::parse_size_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--quiet") {
            options.quiet = true;
            continue;
        }
        throw std::runtime_error("Unknown option: " + arg);
    }

    if (gfa_path.empty()) {
        throw std::runtime_error("Missing required input: --gfa <path>");
    }

    if (bubbles_csv_path.empty()) {
        bubbles_csv_path = out_prefix + ".bubbles.csv";
    }
    if (bandage_csv_path.empty()) {
        bandage_csv_path = out_prefix + ".bandage_nodes.csv";
    }

    cli::ensure_parent_dir_for_file(bubbles_csv_path);
    cli::ensure_parent_dir_for_file(bandage_csv_path);
    if (!snarl_debug_tsv_path.empty()) {
        cli::ensure_parent_dir_for_file(snarl_debug_tsv_path);
    }

    if (options.snarls_input_path.empty()) {
        throw std::runtime_error(
            "Missing required input: --snarls-in <path> (JSONL from 'vg view -R -j')");
    }

    ParseGfaOptions parse_options;
    parse_options.include_paths = true;
    parse_options.include_sequences = true;

    const Graph graph = parse_gfa(gfa_path, parse_options);
    if (graph.paths.empty()) {
        throw std::runtime_error("Input GFA has no P/W paths; snarl refinement requires path walks");
    }
    const auto report = call_bubbles_report(graph, options);
    const auto& bubbles = report.bubbles;

    write_bubbles_csv(bubbles_csv_path, bubbles);
    write_bandage_node_colors_csv(bandage_csv_path, bubbles, report.non_snp_bubbles);
    if (!snarl_debug_tsv_path.empty()) {
        write_snarl_debug_tsv(snarl_debug_tsv_path, report.snarl_debug);
    }

    std::size_t inversion_signal_count = 0;
    std::size_t long_path_positive_count = 0;
    for (const auto& bubble : bubbles) {
        if (bubble.inversion_signal) {
            ++inversion_signal_count;
        }
        if (bubble.long_path_support > 0) {
            ++long_path_positive_count;
        }
    }

    std::cout
        << "Input graph: " << gfa_path << '\n'
        << "Site mode: snarl (JSONL import)\n"
        << "Snarl input: " << options.snarls_input_path << "\n"
        << "Nodes: " << graph.nodes.size() << "\n"
        << "P/W paths loaded: " << graph.paths.size()
        << "\n"
        << "Bubbles called: " << bubbles.size() << "\n"
        << "Min variant bp: " << options.min_variant_bp << " (0 = disabled)\n"
        << "Min path support: " << options.min_path_support << " (0 = disabled)\n"
        << "Merge nearby bp: " << options.merge_nearby_bp << " (0 = disabled)\n"
        << "Bubbles with >=min bp path: " << long_path_positive_count << "\n"
        << "Bubbles kept by inversion signal: " << inversion_signal_count << "\n"
        << "Non-SNP candidate bubbles (for Bandage context): " << report.non_snp_bubbles.size() << "\n"
        << "Wrote: " << bubbles_csv_path << "\n"
        << "Wrote: " << bandage_csv_path << "\n";
    if (!snarl_debug_tsv_path.empty()) {
        std::cout << "Wrote: " << snarl_debug_tsv_path << "\n";
    }

    return 0;
}

} // namespace panvar
