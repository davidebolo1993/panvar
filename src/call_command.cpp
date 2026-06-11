#include "panvar/call_command.hpp"

#include "panvar/cli_utils.hpp"
#include "panvar/gfa.hpp"
#include "panvar/variant_call.hpp"

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace panvar {
namespace {

void print_call_help() {
    std::cout
        << "Usage:\n"
        << "  panvar call -i <graph.gfa> (--bubble-prefix-in <prefix> | --bubbles-csv-in <path>)\n"
        << "             --reference-path <name> -o <out_prefix> [options]\n\n"
        << "Graph-native structural variant caller. Compares each haplotype's bubble walk to a\n"
        << "reference walk, merges fragmented same-type events within a bubble and equivalent\n"
        << "events across haplotypes, and writes a multi-sample VCF per bubble plus a region VCF.\n"
        << "Input is expected to be a panphorte-normalized GFA (a tandem DUP is a REP node looped).\n\n"
        << "Options:\n"
        << "  -i, --gfa <path>                 Input GFA file (required)\n"
        << "      --bubble-prefix-in <prefix>  Module-1 output prefix (auto-uses <prefix>.bubbles.csv)\n"
        << "      --bubbles-csv-in <path>      Module-1 bubbles CSV (required if no prefix)\n"
        << "      --reference-path <name>      Reference path name used as the diff baseline (required)\n"
        << "  -o, --out-prefix <prefix>        Output prefix for VCFs (required)\n"
        << "      --min-sv-bp <N>              Minimum event size to report (default: 50)\n"
        << "      --merge-distance-bp <N>      Coalesce nearby same-type events within a bubble (default: 100)\n"
        << "      --merge-jaccard <X>          Cross-haplotype node-set Jaccard to merge events (default: 0.80)\n"
        << "      --classify-ins               Refine INS subtype NOVEL/DUP with minimap2\n"
        << "      --minimap-preset <name>      minimap2 preset for INS refinement: asm5|asm10|asm20 (default: asm20)\n"
        << "      --minimap-best-n <N>         minimap2 best_n chains (default: 8)\n"
        << "      --ins-dup-min-identity <X>   Identity for an INS to be subtyped DUP (default: 0.90)\n"
        << "      --bubble-id <N>              Restrict to one bubble ID (repeatable)\n"
        << "      --no-per-bubble-vcf          Only write the concatenated region VCF\n"
        << "      --quiet                      Disable progress logs\n"
        << "  -h, --help                       Show this help\n";
}

} // namespace

int run_call_command(const std::vector<std::string>& args) {
    std::string gfa_path;
    std::string bubble_prefix_in;
    std::string out_prefix;
    VariantCallOptions options;

    for (std::size_t i = 0; i < args.size(); ++i) {
        const std::string& arg = args[i];
        auto require_value = [&](const std::string& flag) -> const std::string& {
            if (i + 1 >= args.size()) {
                throw std::runtime_error("Missing value after " + flag);
            }
            return args[++i];
        };

        if (arg == "-h" || arg == "--help") {
            print_call_help();
            return 0;
        }
        if (arg == "-i" || arg == "--gfa") {
            gfa_path = require_value(arg);
            continue;
        }
        if (arg == "--bubble-prefix-in") {
            bubble_prefix_in = require_value(arg);
            continue;
        }
        if (arg == "--bubbles-csv-in") {
            options.bubbles_csv_in = require_value(arg);
            continue;
        }
        if (arg == "--reference-path") {
            options.reference_path = require_value(arg);
            continue;
        }
        if (arg == "-o" || arg == "--out-prefix") {
            out_prefix = require_value(arg);
            continue;
        }
        if (arg == "--min-sv-bp") {
            options.min_sv_bp = cli::parse_size_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--merge-distance-bp") {
            options.merge_distance_bp = cli::parse_size_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--merge-jaccard") {
            options.merge_jaccard = cli::parse_unit_fraction_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--classify-ins") {
            options.classify_ins = true;
            continue;
        }
        if (arg == "--minimap-preset") {
            options.minimap_preset = require_value(arg);
            continue;
        }
        if (arg == "--minimap-best-n") {
            options.minimap_best_n = cli::parse_size_arg(arg, require_value(arg));
            if (options.minimap_best_n == 0) {
                throw std::runtime_error("--minimap-best-n must be >= 1");
            }
            continue;
        }
        if (arg == "--ins-dup-min-identity") {
            options.ins_dup_min_identity = cli::parse_unit_fraction_arg(arg, require_value(arg));
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
        if (arg == "--no-per-bubble-vcf") {
            options.write_per_bubble_vcf = false;
            continue;
        }
        if (arg == "--quiet") {
            options.quiet = true;
            continue;
        }
        throw std::runtime_error("Unknown option for call: " + arg);
    }

    if (gfa_path.empty()) {
        throw std::runtime_error("Missing required input: --gfa <path>");
    }
    if (!bubble_prefix_in.empty()) {
        const std::string derived = bubble_prefix_in + ".bubbles.csv";
        if (options.bubbles_csv_in.empty()) {
            options.bubbles_csv_in = derived;
        } else if (options.bubbles_csv_in != derived) {
            throw std::runtime_error(
                "Conflicting bubble inputs: --bubble-prefix-in resolves to '" +
                derived + "' but --bubbles-csv-in is '" + options.bubbles_csv_in + "'");
        }
    }
    if (options.bubbles_csv_in.empty()) {
        throw std::runtime_error("call requires --bubble-prefix-in <prefix> or --bubbles-csv-in <path>");
    }
    if (options.reference_path.empty()) {
        throw std::runtime_error("--reference-path is required for module 'call'");
    }
    if (out_prefix.empty()) {
        throw std::runtime_error("call requires -o/--out-prefix <prefix>");
    }
    options.out_prefix = out_prefix;

    ParseGfaOptions parse_options;
    parse_options.include_paths = true;
    parse_options.include_sequences = true;
    const Graph graph = parse_gfa(gfa_path, parse_options);
    if (graph.paths.empty()) {
        throw std::runtime_error("Input GFA has no P/W paths; variant calling requires paths");
    }

    VariantCallSummary summary;
    call_variants(graph, options, &summary);

    std::cout
        << "Input graph: " << gfa_path << "\n"
        << "Bubble source: " << options.bubbles_csv_in << "\n"
        << "Reference path: " << options.reference_path << "\n"
        << "Nodes: " << graph.nodes.size() << "\n"
        << "P/W paths loaded: " << graph.paths.size() << "\n"
        << "Bubbles seen: " << summary.bubbles_seen << "\n"
        << "Bubbles with reference: " << summary.bubbles_with_reference << "\n"
        << "Bubbles with calls: " << summary.bubbles_with_calls << "\n"
        << "VCF records: " << summary.records_written << "\n"
        << "  DEL: " << summary.del << "  INS: " << summary.ins
        << "  INV: " << summary.inv << "  DUP: " << summary.dup << "\n"
        << "Region VCF: " << options.out_prefix << ".region.vcf\n";
    if (options.write_per_bubble_vcf) {
        std::cout << "Per-bubble VCFs: " << options.out_prefix << ".bubble_<id>.vcf\n";
    }

    return 0;
}

} // namespace panvar
