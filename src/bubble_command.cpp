#include "panvar/bubble_command.hpp"

#include "panvar/bubbles.hpp"
#include "panvar/cli_utils.hpp"
#include "panvar/gfa.hpp"
#include "panvar/gfa_io.hpp"
#include "panvar/graph_sort.hpp"
#include "panvar/integrated_snarls.hpp"
#include "panvar/output.hpp"

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace panvar {
namespace {

void print_bubble_help() {
    std::cout
        << "Usage:\n"
        << "  panvar bubble -i <graph.gfa> -r <name> [-o <prefix>] [--superbubbles] [options]\n\n"
        << "By default the graph is sorted+flipped along the reference internally and snarls are\n"
        << "found with an internal cactus decomposition (the same 3-edge-connected-component method\n"
        << "as 'vg snarls'); with --superbubbles only the acyclic snarls (ultrabubbles, i.e. the\n"
        << "acyclic subset of the cactus snarls) are kept. Pass --snarls-in to use an external vg\n"
        << "snarls JSONL instead.\n\n"
        << "Options:\n"
        << "  -i, --gfa <path>                 Input GFA file (required)\n"
        << "  -r, --reference-path <name>      Reference path name or unique substring; orders the\n"
        << "                                    internal sort/flip + snarl finder (required unless\n"
        << "                                    --snarls-in is given)\n"
        << "  -s, --superbubbles               Emit only acyclic superbubbles (default: all snarls)\n"
        << "      --no-flip                    Do not reorient nodes to the reference forward strand\n"
        << "      --sorted-gfa-out <path>      Internally-sorted GFA output (default: <prefix>.sorted.gfa)\n"
        << "      --emit-snarls-jsonl <path>   Also write the internal snarls as a vg-style JSONL\n"
        << "      --snarls-in <path>           Override: snarl JSONL from 'vg snarls -A integrated | vg view -R -j' (skips\n"
        << "                                    internal sort + finding)\n"
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
        << "  -q, --quiet                      Disable the progress bar\n"
        << "  -h, --help                       Show this help\n";
}

} // namespace

int run_bubble_command(const std::vector<std::string>& args) {
    if (args.empty()) {
        print_bubble_help();
        return 0;
    }

    std::string gfa_path;
    std::string out_prefix = "bubble_calls";
    std::string bubbles_csv_path;
    std::string bandage_csv_path;
    std::string snarl_debug_tsv_path;
    std::string sorted_gfa_path;
    std::string emit_snarls_jsonl_path;
    bool no_flip = false;

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
        if (arg == "-r" || arg == "--reference-path") {
            options.reference_path = require_value(arg);
            continue;
        }
        if (arg == "-s" || arg == "--superbubbles") {
            options.superbubbles_only = true;
            continue;
        }
        if (arg == "--no-flip") {
            no_flip = true;
            continue;
        }
        if (arg == "--sorted-gfa-out") {
            sorted_gfa_path = require_value(arg);
            continue;
        }
        if (arg == "--emit-snarls-jsonl") {
            emit_snarls_jsonl_path = require_value(arg);
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
        if (arg == "-q" || arg == "--quiet") {
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

    ParseGfaOptions parse_options;
    parse_options.include_paths = true;
    parse_options.include_sequences = true;

    std::string site_mode;
    std::string effective_gfa = gfa_path;
    Graph graph;

    if (options.snarls_input_path.empty()) {
        // Default: internally sort+flip along the reference and find snarls (no vg/odgi).
        if (options.reference_path.empty()) {
            throw std::runtime_error(
                "Missing required input: --reference-path <name> (or pass --snarls-in to use vg)");
        }
        if (sorted_gfa_path.empty()) {
            sorted_gfa_path = out_prefix + ".sorted.gfa";
        }
        cli::ensure_parent_dir_for_file(sorted_gfa_path);

        GfaModel model = read_gfa_model(gfa_path);
        GraphSortOptions sort_opts;
        sort_opts.reference_path = options.reference_path;
        sort_opts.flip = !no_flip;
        sort_graph_reference(model, sort_opts);
        write_gfa_model(sorted_gfa_path, model);
        effective_gfa = sorted_gfa_path;

        graph = parse_gfa(sorted_gfa_path, parse_options);
        if (graph.paths.empty()) {
            throw std::runtime_error("Input GFA has no P/W paths; snarl finding requires path walks");
        }

        // Find boundary pairs internally with the vg-faithful cactus finder. --superbubbles
        // then keeps only the acyclic snarls (= superbubbles), filtered in call_bubbles_report.
        options.snarl_pairs_override = find_top_level_snarls_cactus(snarl_input_from_model(model));
        site_mode = options.superbubbles_only ? "superbubble (internal, cactus + acyclic)"
                                              : "snarl (internal, cactus)";
    } else {
        // Legacy override: external vg snarls on the graph as-is (no internal sort).
        graph = parse_gfa(gfa_path, parse_options);
        if (graph.paths.empty()) {
            throw std::runtime_error("Input GFA has no P/W paths; snarl refinement requires path walks");
        }
        site_mode = "snarl (JSONL import)";
    }

    const auto report = call_bubbles_report(graph, options);
    const auto& bubbles = report.bubbles;

    write_bubbles_csv(bubbles_csv_path, bubbles);
    write_bandage_node_colors_csv(bandage_csv_path, bubbles, report.non_snp_bubbles);
    if (!snarl_debug_tsv_path.empty()) {
        write_snarl_debug_tsv(snarl_debug_tsv_path, report.snarl_debug);
    }
    if (!emit_snarls_jsonl_path.empty()) {
        // Emit the internally found (cactus) snarl pairs. In --snarls-in (legacy) mode the
        // snarls already exist as the input file, so there is nothing new to emit.
        if (options.snarl_pairs_override.empty()) {
            std::cerr << "note: --emit-snarls-jsonl is a no-op with --snarls-in (the input IS the snarls)\n";
        } else {
            cli::ensure_parent_dir_for_file(emit_snarls_jsonl_path);
            std::ofstream js(emit_snarls_jsonl_path);
            for (const auto& [s, t] : options.snarl_pairs_override) {
                js << "{\"start\": {\"node_id\": \"" << s << "\"}, \"end\": {\"node_id\": \"" << t
                   << "\"}}\n";
            }
        }
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
        << "Site mode: " << site_mode << "\n"
        << (options.snarls_input_path.empty()
                ? ("Sorted graph: " + effective_gfa + "\n")
                : ("Snarl input: " + options.snarls_input_path + "\n"))
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
