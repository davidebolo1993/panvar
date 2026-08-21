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
        << "      --min-alt-af <X>             Drop records whose ALT carrier frequency (carriers over\n"
        << "                                   traversing haplotypes) is below X. Not a minor-allele\n"
        << "                                   frequency: a DUP's reference-like state is not an allele\n"
        << "                                   here (default: 0=off; --min-maf is an accepted alias)\n"
        << "      --allele-vcf                 Also write <prefix>.alleles.vcf: one record per bubble carrying\n"
        << "                                   EVERY distinct allele as explicit sequence, each haplotype's GT\n"
        << "                                   indexing its own. Lossless companion to the merged region VCF,\n"
        << "                                   which describes one walk with several interpreted records\n"
        << "      --allele-vcf-max-bp <N>      Skip a bubble in the allele VCF if any allele exceeds N bp\n"
        << "                                   (0 = no limit, the default)\n"
        << "      --multiallelic-loci          Collapse a bounded locus into ONE multiallelic record (REF+ALT1,ALT2..)\n"
        << "      --multiallelic-max-bp <N>    Skip multiallelic collapse above this allele size (default: 5000)\n"
        << "      --rescue-min-bp <N>          Floor for sub-threshold events kept for rescue (default: min-sv-bp/2)\n"
        << "      --classify-ins               Refine INS subtype NOVEL/DUP with minimap2\n"
        << "      --minimap-preset <name>      minimap2 preset for INS refinement: asm5|asm10|asm20 (default: asm20)\n"
        << "      --minimap-best-n <N>         minimap2 best_n chains (default: 8)\n"
        << "      --ins-dup-min-identity <X>   Identity for an INS to be subtyped DUP (default: 0.90)\n"
        << "      --tangle-min-hubs <N>        Suppress DUP/CN calling in a low-complexity tangle: a bubble\n"
        << "                                   with >= N high-degree hub interior nodes (default: 10, 0=off)\n"
        << "      --tangle-hub-degree <N>      Distinct-neighbour degree for an interior node to count as a\n"
        << "                                   tangle hub (default: 20)\n"
        << "      --cn-unit-spacing            Take the MODULE_BP copy-number step from the panel's own\n"
        << "                                   cluster spacing rather than ref_bp/ref_fold. This makes\n"
        << "                                   the calibration panel-dependent and is therefore opt-in\n"
        << "      --max-cn-model-residual <F>  Refuse a MODULE_BP copy-number call when more than this\n"
        << "                                   FRACTION of traversers round ambiguously\n"
        << "                                   (CN_ROUND_AMBIGUOUS_FRAC, 0..1); the bubble then\n"
        << "                                   gets no CN call at all, not a fallback one. 0 = off\n"
        << "      --max-dup-region-frac <F>    Suppress a peak DUP spanning more than this fraction of the\n"
        << "                                   reference (a tangle artifact, not a real dup; default: 0.8, 0=off)\n"
        << "      --cn                         Enable inferred MODULE_BP and PEAK copy-number routes. REP\n"
        << "                                   self-loop traversal counting is automatic when that topology\n"
        << "                                   exists. Routes resolve as REP > MODULE_BP > PEAK and emit the\n"
        << "                                   total copy number of the represented module. With --gtf, a\n"
        << "                                   separate table estimates per-gene CN where paralogs are\n"
        << "                                   separable. Sequence-resolved\n"
        << "                                   DEL/INS/INV are kept alongside the CN call.\n"
        << "      --gtf <path>                 Reference-coordinate GTF: annotate variants with the genes\n"
        << "                                   they touch (INFO GENES), write <prefix>.node_genes.tsv and a\n"
        << "                                   per-gene DUP copy-number table; needs a PanSN reference path\n"
        << "      --bubble-id <N>              Restrict to one bubble ID (repeatable)\n"
        << "      --no-per-bubble-vcf          Only write the concatenated region VCF\n"
        << "      --no-variant-nodes           Skip the variant_nodes.tsv sidecar\n"
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
        if (arg == "--min-alt-af" || arg == "--min-maf") {   // --min-maf is the historical spelling
            options.min_maf = cli::parse_unit_fraction_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--allele-vcf") {
            options.allele_vcf = true;
            continue;
        }
        if (arg == "--allele-vcf-max-bp") {
            options.allele_vcf_max_bp = cli::parse_size_arg(arg, require_value(arg));
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
        if (arg == "--cn") {                 // copy-number calling: all routes, resolved by topology
            options.cn = true;
            continue;
        }
        if (arg == "--tangle-hub-degree") {
            options.tangle_hub_degree = cli::parse_size_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--tangle-min-hubs") {    // 0 disables the low-complexity tangle guard
            options.tangle_min_hubs = cli::parse_size_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--max-dup-region-frac") {  // suppress a peak DUP spanning > this fraction of the reference
            // A fraction of the reference, so 0..1. Raw std::stod accepted 2, and "nan" -- which then
            // compared false against everything and silently disabled the guard.
            options.max_dup_region_frac = cli::parse_unit_fraction_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--cn-unit-spacing") {   // take the MODULE_BP copy step from the panel, not ref_bp/ref_fold
            options.cn_unit_spacing = true;
            continue;
        }
        if (arg == "--max-cn-model-residual") {
            options.max_cn_model_residual = cli::parse_unit_fraction_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--gtf") {
            options.gtf_path = require_value(arg);
            continue;
        }
        if (arg == "--minimap-preset") {
            // Only the three the help advertises; anything else reaches minimap2 as an unknown preset.
            options.minimap_preset = require_value(arg);
            if (options.minimap_preset != "asm5" && options.minimap_preset != "asm10" &&
                options.minimap_preset != "asm20") {
                throw std::runtime_error("--minimap-preset must be asm5, asm10 or asm20 (got '" +
                                         options.minimap_preset + "')");
            }
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
        if (arg == "--no-variant-nodes") {
            options.write_variant_nodes = false;
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
    options.gfa_path = gfa_path;
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

    cli::RunLog log("call", options.quiet);
    log.info("input " + gfa_path + " (" + std::to_string(graph.nodes.size()) + " nodes, " +
             std::to_string(graph.paths.size()) + " paths); reference " + options.reference_path);

    VariantCallSummary summary;
    call_variants(graph, options, &summary);

    log.info("called " + std::to_string(summary.records_written) + " records across " +
             std::to_string(summary.bubbles_with_calls) + "/" + std::to_string(summary.bubbles_seen) +
             " bubbles (DEL " + std::to_string(summary.del) + ", INS " + std::to_string(summary.ins) +
             ", INV " + std::to_string(summary.inv) + ", DUP " + std::to_string(summary.dup) +
             ", MULTI " + std::to_string(summary.multi) + ")");
    if (summary.oversized_dups > 0) {
        log.info("suppressed " + std::to_string(summary.oversized_dups) +
                 " copy-number DUP call(s) as low-complexity-tangle or oversized (> max-dup-region-frac)" +
                 (summary.tangle_bubbles > 0
                      ? "; " + std::to_string(summary.tangle_bubbles) + " bubble(s) flagged tangle"
                      : ""));
    }
    // Both of these are guaranteed losses, and a quiet run previously left no record of either: the
    // counters reached the summary struct and nothing printed them, so a missing CN call or an
    // unaligned divergent block had no durable reason anywhere in the output.
    if (summary.skipped_large_segments > 0) {
        log.info("skipped " + std::to_string(summary.skipped_large_segments) +
                 " divergent block(s) too large to align; each is a false negative, not a quiet pass");
    }
    if (summary.declined_cn_model > 0) {
        log.info("declined " + std::to_string(summary.declined_cn_model) +
                 " module copy-number call(s) over --max-cn-model-residual, leaving " +
                 std::to_string(summary.declined_cn_model_bubbles) +
                 " bubble(s) with no CN record (sequence-resolved events are unaffected)");
    }

    std::vector<std::string> outputs;
    outputs.push_back(options.out_prefix + ".region.vcf");
    if (options.allele_vcf) {
        outputs.push_back(options.out_prefix + ".alleles.vcf");
    }
    if (options.write_per_bubble_vcf) {
        outputs.push_back(options.out_prefix + ".bubble_<id>.vcf");
    }
    if (options.write_variant_nodes) {
        outputs.push_back(options.out_prefix + ".variant_nodes.tsv");
    }
    if (!options.gtf_path.empty()) {
        outputs.push_back(options.out_prefix + ".node_genes.tsv");
        outputs.push_back(options.out_prefix + ".dup_gene_cn.tsv");
    }
    log.wrote(outputs);
    log.done();
    return 0;
}

} // namespace panvar
