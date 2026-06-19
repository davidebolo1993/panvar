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
        << "  panvar call -i <graph.gfa> (-b <prefix> | -c <bubbles.csv>) -r <name> -o <prefix> [options]\n\n"
        << "Graph-native structural variant caller. Compares each haplotype's bubble walk to a\n"
        << "reference walk, merges fragmented same-type events within a bubble and equivalent\n"
        << "events across haplotypes, and writes a multi-sample VCF per bubble plus a region VCF.\n"
        << "Input is expected to be a panphorte-normalized GFA (a tandem DUP is a REP node looped).\n\n"
        << "Options:\n"
        << "  -i, --gfa <path>                 panphorte-normalized/sorted GFA, i.e.\n"
        << "                                   <panphorte_prefix>.normalized.sorted.gfa (required)\n"
        << "  -b, --bubble-prefix-in <prefix>  panphorte output prefix (auto-uses <prefix>.bubbles.csv)\n"
        << "  -c, --bubbles-csv-in <path>      panphorte bubbles CSV (required if no prefix)\n"
        << "  -r, --reference-path <name>      Reference path name used as the diff baseline (required)\n"
        << "  -o, --out-prefix <prefix>        Output prefix for VCFs (required)\n"
        << "      --min-sv-bp <N>              Minimum size of a reported (merged) event (default: 50)\n"
        << "      --merge-distance-bp <N>      Coalesce nearby same-type events within a bubble (default: 100)\n"
        << "      --merge-jaccard <X>          Cross-haplotype node-set Jaccard to merge events (default: 0.80)\n"
        << "      --merge-seq-identity <X>     Cross-haplotype event-sequence identity to merge (default: 0.80)\n"
        << "      --merge-size-ratio <X>       Length-ratio floor for the sequence merge; lower to merge\n"
        << "                                   same-locus events of different sizes, e.g. STR alleles\n"
        << "                                   (default: 0 = use --merge-seq-identity)\n"
        << "      --min-haplotypes <N>         Drop records carried by fewer than N haplotypes (default: 1)\n"
        << "      --min-maf <X>                Drop records with AF (carriers/traversing-haps) below X (default: 0=off)\n"
        << "      --multiallelic-loci          Collapse a bounded locus into ONE multiallelic record (REF+ALT1,ALT2..)\n"
        << "      --multiallelic-max-bp <N>    Skip multiallelic collapse above this allele size (default: 5000)\n"
        << "      --rescue-min-bp <N>          Floor for sub-threshold events kept for rescue (default: min-sv-bp/2)\n"
        << "      --classify-ins               Refine INS subtype NOVEL/DUP with minimap2\n"
        << "      --minimap-preset <name>      minimap2 preset for INS refinement: asm5|asm10|asm20 (default: asm20)\n"
        << "      --minimap-best-n <N>         minimap2 best_n chains (default: 8)\n"
        << "      --ins-dup-min-identity <X>   Identity for an INS to be subtyped DUP (default: 0.90)\n"
        << "      --cn-from-multiplicity       Emit DUP from peak node multiplicity for folded bubbles\n"
        << "                                   with no self-loop (e.g. GSTM1) that panphorte left intact\n"
        << "      --cn-from-coverage           Emit total-module copy number (spelled-bp/unit) on folded\n"
        << "                                   paralog clusters (reference folds >=2x, e.g. CYP2D6);\n"
        << "                                   recovers deletions+gains; precedence over --cn-from-multiplicity\n"
        << "      --bubble-id <N>              Restrict to one bubble ID (repeatable)\n"
        << "      --no-per-bubble-vcf          Only write the concatenated region VCF\n"
        << "      --no-variant-paths           Skip the variant_paths.tsv + node_track.tsv sidecars\n"
        << "      --threads <N>                Worker threads for the per-bubble loop (0 = auto)\n"
        << "  -q, --quiet                      Disable progress logs\n"
        << "  -h, --help                       Show this help\n";
}

} // namespace

int run_call_command(const std::vector<std::string>& args) {
    if (args.empty()) {
        print_call_help();
        return 0;
    }

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
        if (arg == "-b" || arg == "--bubble-prefix-in") {
            bubble_prefix_in = require_value(arg);
            continue;
        }
        if (arg == "-c" || arg == "--bubbles-csv-in") {
            options.bubbles_csv_in = require_value(arg);
            continue;
        }
        if (arg == "-r" || arg == "--reference-path") {
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
        if (arg == "--merge-seq-identity") {
            options.merge_seq_identity = cli::parse_unit_fraction_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--merge-size-ratio") {
            options.merge_size_ratio = cli::parse_unit_fraction_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--min-haplotypes") {
            options.min_haplotypes = cli::parse_size_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--min-maf") {
            options.min_maf = cli::parse_unit_fraction_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--multiallelic-loci") {
            options.multiallelic_loci = true;
            continue;
        }
        if (arg == "--multiallelic-max-bp") {
            options.multiallelic_max_bp = cli::parse_size_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--rescue-min-bp") {
            options.rescue_min_bp = cli::parse_size_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--classify-ins") {
            options.classify_ins = true;
            continue;
        }
        if (arg == "--cn-from-multiplicity") {
            options.cn_from_multiplicity = true;
            continue;
        }
        if (arg == "--cn-from-coverage") {
            options.cn_from_coverage = true;
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
        if (arg == "--no-variant-paths") {
            options.write_variant_paths = false;
            continue;
        }
        if (arg == "--threads") {
            options.threads = cli::parse_size_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "-q" || arg == "--quiet") {
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
