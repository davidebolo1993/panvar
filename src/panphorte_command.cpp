#include "panvar/panphorte_command.hpp"

#include "panvar/cli_utils.hpp"
#include "panvar/panphorte.hpp"

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace panvar {
namespace {

void print_panphorte_help() {
    std::cout
        << "Usage:\n"
        << "  panvar panphorte -i <graph.gfa> --bubble-prefix-in <module1-prefix> -o <out_prefix> [options]\n\n"
        << "Normalizes tandem-repeat bubbles into a compact representation (repeat unit + cycle)\n"
        << "and writes a new GFA that can be re-processed by bubble/inspect.\n\n"
        << "Options:\n"
        << "  -i, --gfa <path>                 Input GFA file (required)\n"
        << "      --bubble-prefix-in <prefix>  Module-1 output prefix (auto-uses <prefix>.bubbles.csv)\n"
        << "      --bubbles-csv <path>         Module-1 bubbles CSV (required if no prefix)\n"
        << "  -o, --out-prefix <prefix>        Output prefix for <prefix>.normalized.gfa + report (required)\n"
        << "      --bubble-id <N>              Restrict to one bubble ID (repeatable)\n"
        << "      --min-unit-bp <N>            Minimum repeat-unit span to normalize (default: 50)\n"
        << "      --min-copies <N>             Minimum tandem copies to normalize (default: 2)\n"
        << "      --max-interruption-frac <f>  Tolerance for interrupting bases in an array (default: 0.25)\n"
        << "      --min-similarity <f>         Min identity to treat a block as a unit copy\n"
        << "                                   (1.0=exact; <1.0 enables approximate, lossy collapse; default: 1.0)\n"
        << "      --threads <N>                Worker threads for approximate detection (0=auto; default: 0)\n"
        << "      --reference-path <name>      If set, internally sort+flip the normalized graph along\n"
        << "                                   this reference and re-snarl (cactus), writing\n"
        << "                                   <prefix>.normalized.sorted.gfa + <prefix>.bubbles.csv so\n"
        << "                                   'call' runs with no external vg/odgi tools\n"
        << "      --no-flip                    With --reference-path, skip reorienting to the ref strand\n"
        << "      --quiet                      Disable progress logs\n"
        << "  -h, --help                       Show this help\n";
}

} // namespace

int run_panphorte_command(const std::vector<std::string>& args) {
    std::string bubble_prefix_in;
    PanphorteOptions options;

    for (std::size_t i = 0; i < args.size(); ++i) {
        const std::string& arg = args[i];
        auto require_value = [&](const std::string& flag) -> const std::string& {
            if (i + 1 >= args.size()) {
                throw std::runtime_error("Missing value after " + flag);
            }
            return args[++i];
        };

        if (arg == "-h" || arg == "--help") {
            print_panphorte_help();
            return 0;
        }
        if (arg == "-i" || arg == "--gfa") {
            options.gfa_path = require_value(arg);
            continue;
        }
        if (arg == "--bubble-prefix-in") {
            bubble_prefix_in = require_value(arg);
            continue;
        }
        if (arg == "--bubbles-csv") {
            options.bubbles_csv_in = require_value(arg);
            continue;
        }
        if (arg == "-o" || arg == "--out-prefix") {
            options.out_prefix = require_value(arg);
            continue;
        }
        if (arg == "--bubble-id") {
            const std::size_t id = cli::parse_size_arg(arg, require_value(arg));
            if (id == 0) {
                throw std::runtime_error("--bubble-id must be > 0");
            }
            options.bubble_ids.push_back(id);
            continue;
        }
        if (arg == "--min-unit-bp") {
            options.min_unit_bp = cli::parse_size_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--min-copies") {
            options.min_copies = cli::parse_size_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--max-interruption-frac") {
            options.max_interruption_frac = cli::parse_unit_fraction_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--min-similarity") {
            options.min_similarity = cli::parse_similarity_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--threads") {
            options.threads = cli::parse_size_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--reference-path") {
            options.reference_path = require_value(arg);
            continue;
        }
        if (arg == "--no-flip") {
            options.no_flip = true;
            continue;
        }
        if (arg == "--quiet") {
            options.quiet = true;
            continue;
        }
        throw std::runtime_error("Unknown option for panphorte: " + arg);
    }

    if (options.gfa_path.empty()) {
        throw std::runtime_error("panphorte requires -i/--gfa <graph.gfa>");
    }
    if (options.out_prefix.empty()) {
        throw std::runtime_error("panphorte requires -o/--out-prefix <prefix>");
    }
    if (!bubble_prefix_in.empty()) {
        const std::string derived = bubble_prefix_in + ".bubbles.csv";
        if (options.bubbles_csv_in.empty()) {
            options.bubbles_csv_in = derived;
        } else if (options.bubbles_csv_in != derived) {
            throw std::runtime_error(
                "Conflicting bubble inputs: --bubble-prefix-in resolves to '" +
                derived + "' but --bubbles-csv is '" + options.bubbles_csv_in + "'");
        }
    }
    if (options.bubbles_csv_in.empty()) {
        throw std::runtime_error("panphorte requires --bubble-prefix-in <prefix> or --bubbles-csv <path>");
    }
    if (options.min_copies < 2) {
        throw std::runtime_error("--min-copies must be >= 2");
    }

    PanphorteSummary summary;
    panphorte_normalize(options, &summary);

    std::cout
        << "Input graph: " << options.gfa_path << "\n"
        << "Bubble source: " << options.bubbles_csv_in << "\n"
        << "Output GFA: " << options.out_prefix << ".normalized.gfa\n"
        << "Report: " << options.out_prefix << ".panphorte.report.tsv\n"
        << "Bubbles seen: " << summary.bubbles_seen << "\n"
        << "Bubbles normalized: " << summary.bubbles_normalized << "\n"
        << "Paths rewritten: " << summary.paths_rewritten << "\n"
        << "Nodes removed: " << summary.nodes_removed << "\n"
        << "Nodes added: " << summary.nodes_added << "\n"
        << "Edges added: " << summary.edges_added << "\n";

    if (summary.sorted) {
        std::cout
            << "Sorted GFA: " << options.out_prefix << ".normalized.sorted.gfa\n"
            << "Bubbles CSV: " << options.out_prefix << ".bubbles.csv\n"
            << "Re-snarled bubbles: " << summary.resnarled_bubbles << "\n"
            << "Ready for: panvar call -i " << options.out_prefix
            << ".normalized.sorted.gfa --bubble-prefix-in " << options.out_prefix << "\n";
    }

    return 0;
}

} // namespace panvar
