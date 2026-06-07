#include "panvar/call_command.hpp"

#include "panvar/allele.hpp"
#include "panvar/call.hpp"
#include "panvar/cli_utils.hpp"
#include "panvar/gfa.hpp"

#include <algorithm>
#include <cctype>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace panvar {
namespace {

void print_call_help() {
    std::cout
        << "Usage:\n"
        << "  panvar call -i <graph.gfa> [--bubble-prefix-in <module1-prefix> | --bubbles-csv-in <module1.bubbles.csv>] [--allele-prefix-in <module2-prefix> | (--clusters-csv-in <module2.allele_clusters.csv> --assignments-csv-in <module2.allele_assignments.csv>)] [options]\n\n"
        << "Options:\n"
        << "  -i, --gfa <path>                 Input GFA file (required)\n"
        << "  -o, --out-prefix <prefix>        Output prefix (default: call_calls)\n"
        << "      --bubble-prefix-in <prefix>  Module-1 output prefix from 'panvar bubble'\n"
        << "                                   (auto-uses <prefix>.bubbles.csv)\n"
        << "      --allele-prefix-in <prefix>  Module-2 output prefix from 'panvar allele'\n"
        << "                                   (auto-uses <prefix>.allele_clusters.csv and\n"
        << "                                    <prefix>.allele_assignments.csv)\n"
        << "      --bubbles-csv-in <path>      Module-1 bubbles CSV input (required if no prefix)\n"
        << "      --clusters-csv-in <path>     Module-2 allele-cluster CSV input (required if no allele prefix)\n"
        << "      --assignments-csv-in <path>  Module-2 per-path assignment CSV input (required if no allele prefix)\n"
        << "      --quiet                      Disable progress logs\n"
        << "      --reference-path <path>      Reference path name (required)\n"
        << "      --skip-no-reference-bubbles   Skip bubbles where reference has no assignment\n"
        << "      --min-sv-bp <N>              Minimum size for INS/DEL-style calls (default: 50)\n"
        << "      --vcf-out <path>             Region-level multi-sample VCF path (default: <out>.region.vcf)\n"
        << "      --debug                      Write per-cluster debug artifacts (PAF, FASTA, dotplot, VCF)\n"
        << "      --debug-out-dir <dir>        Debug output directory (enables --debug)\n"
        << "      --dotplot-gtf <path>         Optional GTF (gene features) overlay for debug dotplots\n"
        << "      --dotplot-gene-match <expr>  Optional gene regex/term filter for dotplot highlighting\n"
        << "                                   (repeat flag to pass multiple patterns)\n"
        << "      --minimap-preset <name>      minimap2 preset: asm5|asm10|asm20 (default: asm20)\n"
        << "      --minimap-best-n <N>         minimap2 best_n chain count (default: 8)\n"
        << "      --minimap-no-secondary       Disable secondary/supplementary alignments in minimap2 output\n"
        << "      --split-ins-svlen-mode <m>   query-span|geometric (default: geometric)\n"
        << "      --classify-ins               Classify INS as NOVEL/DUP-like and estimate DUP copy numbers\n"
        << "      --pangene-bed <path>         Optional pangene BED/BED.GZ for gene copy annotations\n"
        << "      --pangene-gene-match <expr>  Optional pangene gene regex/term filter (repeatable)\n"
        << "      --pangene-tune-ins           Use pangene copy gains to set INS subtype DUP_PANGENE\n"
        << "      --pangene-copy-tsv <path>    Optional per-bubble/cluster pangene copy-count TSV\n"
        << "      --vcf-merge-window-bp <N>    Cross-cluster event merge window in bp (default: 20)\n"
        << "      --vcf-merge-mode <m>         strict|lenient (default: strict)\n"
        << "      --vcf-merge-lenient-window-bp <N>\n"
        << "                                   Wider start window used in lenient mode (default: 100)\n"
        << "      --vcf-merge-lenient-min-ref-jaccard <X>\n"
        << "                                   Min reference interval Jaccard overlap in lenient mode\n"
        << "                                   when strict merge checks fail (default: 0.30)\n"
        << "      --vcf-merge-min-seq-sim <X>  Cross-cluster min sequence similarity in [0,1] (default: 0.80)\n"
        << "      --vcf-merge-max-edit-frac <X> Max edit fraction for merge comparisons in [0,1] (default: 0.35)\n"
        << "  -h, --help                       Show this help\n";
}

std::string to_lower_ascii(std::string value) {
    std::transform(
        value.begin(),
        value.end(),
        value.begin(),
        [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return value;
}

std::string normalize_minimap_preset(const std::string& value) {
    const std::string lowered = to_lower_ascii(value);
    if (lowered == "asm5" || lowered == "asm10" || lowered == "asm20") {
        return lowered;
    }
    throw std::runtime_error(
        "Invalid value for --minimap-preset: " + value +
        " (supported: asm5, asm10, asm20)");
}

bool parse_split_ins_svlen_geometric(const std::string& value) {
    const std::string lowered = to_lower_ascii(value);
    if (lowered == "geometric" || lowered == "svim" || lowered == "deviation") {
        return true;
    }
    if (lowered == "query-span" || lowered == "query_span" || lowered == "span" || lowered == "legacy") {
        return false;
    }
    throw std::runtime_error(
        "Invalid --split-ins-svlen-mode: " + value +
        " (expected query-span|geometric)");
}

std::string normalize_vcf_merge_mode(const std::string& value) {
    const std::string lowered = to_lower_ascii(value);
    if (lowered == "strict" || lowered == "lenient") {
        return lowered;
    }
    throw std::runtime_error(
        "Invalid --vcf-merge-mode: " + value +
        " (expected strict|lenient)");
}

} // namespace

int run_call_command(const std::vector<std::string>& args) {
    std::string gfa_path;
    std::string out_prefix = "call_calls";
    std::string bubble_prefix_in;
    std::string allele_prefix_in;
    std::string clusters_csv_in_path;
    std::string assignments_csv_in_path;

    AlleleCallOptions options;
    options.write_region_vcf = true;

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
        if (arg == "-o" || arg == "--out-prefix") {
            out_prefix = require_value(arg);
            continue;
        }
        if (arg == "--bubble-prefix-in") {
            bubble_prefix_in = require_value(arg);
            continue;
        }
        if (arg == "--allele-prefix-in") {
            allele_prefix_in = require_value(arg);
            continue;
        }
        if (arg == "--bubbles-csv-in") {
            options.bubbles_csv_in = require_value(arg);
            continue;
        }
        if (arg == "--clusters-csv-in") {
            clusters_csv_in_path = require_value(arg);
            continue;
        }
        if (arg == "--assignments-csv-in") {
            assignments_csv_in_path = require_value(arg);
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
        if (arg == "--skip-no-reference-bubbles") {
            options.skip_bubbles_without_reference = true;
            continue;
        }
        if (arg == "--min-sv-bp") {
            options.min_sv_bp = cli::parse_size_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--vcf-out") {
            options.region_vcf_path = require_value(arg);
            options.write_region_vcf = true;
            continue;
        }
        if (arg == "--debug") {
            options.write_debug_reports = true;
            continue;
        }
        if (arg == "--debug-out-dir") {
            options.debug_out_dir = require_value(arg);
            options.write_debug_reports = true;
            continue;
        }
        if (arg == "--dotplot-gtf") {
            options.dotplot_gtf_path = require_value(arg);
            continue;
        }
        if (arg == "--dotplot-gene-match") {
            options.dotplot_gene_matches.push_back(require_value(arg));
            continue;
        }
        if (arg == "--minimap-preset") {
            options.minimap_preset = normalize_minimap_preset(require_value(arg));
            continue;
        }
        if (arg == "--minimap-best-n") {
            options.minimap_best_n = cli::parse_size_arg(arg, require_value(arg));
            if (options.minimap_best_n == 0) {
                throw std::runtime_error("--minimap-best-n must be >= 1");
            }
            continue;
        }
        if (arg == "--minimap-no-secondary") {
            options.minimap_emit_secondary = false;
            continue;
        }
        if (arg == "--split-ins-svlen-mode") {
            options.split_ins_use_geometric_svlen =
                parse_split_ins_svlen_geometric(require_value(arg));
            continue;
        }
        if (arg == "--classify-ins") {
            options.classify_ins = true;
            continue;
        }
        if (arg == "--pangene-bed") {
            options.pangene_bed_path = require_value(arg);
            continue;
        }
        if (arg == "--pangene-gene-match") {
            options.pangene_gene_matches.push_back(require_value(arg));
            continue;
        }
        if (arg == "--pangene-tune-ins") {
            options.pangene_tune_ins = true;
            continue;
        }
        if (arg == "--pangene-copy-tsv") {
            options.pangene_copy_tsv_path = require_value(arg);
            continue;
        }
        if (arg == "--vcf-merge-window-bp") {
            options.vcf_merge_window_bp = cli::parse_size_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--vcf-merge-mode") {
            options.vcf_merge_mode = normalize_vcf_merge_mode(require_value(arg));
            continue;
        }
        if (arg == "--vcf-merge-lenient-window-bp") {
            options.vcf_merge_lenient_window_bp = cli::parse_size_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--vcf-merge-lenient-min-ref-jaccard") {
            options.vcf_merge_lenient_min_ref_jaccard =
                cli::parse_unit_fraction_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--vcf-merge-min-seq-sim") {
            options.vcf_merge_min_seq_similarity =
                cli::parse_unit_fraction_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--vcf-merge-max-edit-frac") {
            options.vcf_merge_max_seq_edit_fraction =
                cli::parse_unit_fraction_arg(arg, require_value(arg));
            continue;
        }
        throw std::runtime_error("Unknown option: " + arg);
    }

    if (gfa_path.empty()) {
        throw std::runtime_error("Missing required input: --gfa <path>");
    }
    if (!bubble_prefix_in.empty()) {
        const std::string derived_bubbles_csv = bubble_prefix_in + ".bubbles.csv";
        if (options.bubbles_csv_in.empty()) {
            options.bubbles_csv_in = derived_bubbles_csv;
        } else if (options.bubbles_csv_in != derived_bubbles_csv) {
            throw std::runtime_error(
                "Conflicting bubble inputs: --bubble-prefix-in resolves to '" +
                derived_bubbles_csv + "' but --bubbles-csv-in is '" + options.bubbles_csv_in + "'");
        }
    }
    if (!allele_prefix_in.empty()) {
        const std::string derived_clusters_csv = allele_prefix_in + ".allele_clusters.csv";
        const std::string derived_assignments_csv = allele_prefix_in + ".allele_assignments.csv";
        if (clusters_csv_in_path.empty()) {
            clusters_csv_in_path = derived_clusters_csv;
        } else if (clusters_csv_in_path != derived_clusters_csv) {
            throw std::runtime_error(
                "Conflicting cluster inputs: --allele-prefix-in resolves to '" +
                derived_clusters_csv + "' but --clusters-csv-in is '" + clusters_csv_in_path + "'");
        }
        if (assignments_csv_in_path.empty()) {
            assignments_csv_in_path = derived_assignments_csv;
        } else if (assignments_csv_in_path != derived_assignments_csv) {
            throw std::runtime_error(
                "Conflicting assignment inputs: --allele-prefix-in resolves to '" +
                derived_assignments_csv + "' but --assignments-csv-in is '" + assignments_csv_in_path + "'");
        }
    }
    if (options.bubbles_csv_in.empty()) {
        throw std::runtime_error(
            "Module 'call' requires module-1 bubbles: use --bubble-prefix-in <prefix> "
            "or --bubbles-csv-in <path>.");
    }
    if (options.reference_path.empty()) {
        throw std::runtime_error("--reference-path is required for module 'call'");
    }
    if (clusters_csv_in_path.empty() || assignments_csv_in_path.empty()) {
        throw std::runtime_error(
            "Module 'call' requires module-2 inputs: use --allele-prefix-in <prefix> "
            "or both --clusters-csv-in and --assignments-csv-in.");
    }
    if (options.region_vcf_path.empty()) {
        options.region_vcf_path = out_prefix + ".region.vcf";
    }
    cli::ensure_parent_dir_for_file(options.region_vcf_path);
    if (options.write_debug_reports && options.debug_out_dir.empty()) {
        options.debug_out_dir = out_prefix + ".debug";
    }
    if (!options.pangene_bed_path.empty() && options.pangene_copy_tsv_path.empty()) {
        options.pangene_copy_tsv_path = out_prefix + ".pangene_copy.tsv";
    }
    if (!options.pangene_copy_tsv_path.empty()) {
        cli::ensure_parent_dir_for_file(options.pangene_copy_tsv_path);
    }
    ParseGfaOptions parse_options;
    parse_options.include_paths = true;
    parse_options.include_sequences = true;
    const Graph graph = parse_gfa(gfa_path, parse_options);

    if (graph.paths.empty()) {
        throw std::runtime_error("Input GFA has no P/W paths; variant calling requires paths");
    }

    AlleleCallSummary summary;
    call_variants_from_precomputed_alleles(
        graph,
        options,
        clusters_csv_in_path,
        assignments_csv_in_path,
        &summary);

    std::cout
        << "Input graph: " << gfa_path << '\n'
        << "Bubble source: " << options.bubbles_csv_in << "\n"
        << "Nodes: " << graph.nodes.size() << "\n"
        << "P/W paths loaded: " << graph.paths.size() << "\n"
        << "Bubbles processed: " << summary.bubbles_processed << "\n"
        << "Bubbles with haplotype assignments: " << summary.bubbles_with_assignments << "\n"
        << "Allele clusters: " << summary.clusters << "\n"
        << "Path assignments: " << summary.assignments << "\n"
        << "Skip no-reference bubbles: " << (options.skip_bubbles_without_reference ? "yes" : "no") << "\n"
        << "Bubbles skipped (no reference assignment): " << summary.bubbles_skipped_no_reference << "\n"
        << "Minimap2 preset: " << options.minimap_preset << "\n"
        << "Minimap2 best_n: " << options.minimap_best_n << "\n"
        << "Minimap2 secondary: " << (options.minimap_emit_secondary ? "on" : "off") << "\n"
        << "Split INS SVLEN mode: "
        << (options.split_ins_use_geometric_svlen ? "geometric" : "query-span") << "\n"
        << "INS classification: " << (options.classify_ins ? "on" : "off") << "\n"
        << "Pangene BED: " << (options.pangene_bed_path.empty() ? "off" : options.pangene_bed_path) << "\n"
        << "Pangene tune INS: " << (options.pangene_tune_ins ? "on" : "off") << "\n"
        << "VCF merge window bp: " << options.vcf_merge_window_bp << "\n"
        << "VCF merge mode: " << options.vcf_merge_mode << "\n"
        << "VCF merge lenient window bp: " << options.vcf_merge_lenient_window_bp << "\n"
        << "VCF merge lenient min ref jaccard: " << options.vcf_merge_lenient_min_ref_jaccard << "\n"
        << "VCF merge min seq sim: " << options.vcf_merge_min_seq_similarity << "\n"
        << "VCF merge max edit frac: " << options.vcf_merge_max_seq_edit_fraction << "\n"
        << "VCF merge effective max edit frac: "
        << std::min(options.vcf_merge_max_seq_edit_fraction, 1.0 - options.vcf_merge_min_seq_similarity) << "\n"
        << "Precomputed clusters CSV: " << clusters_csv_in_path << "\n"
        << "Precomputed assignments CSV: " << assignments_csv_in_path << "\n";
    if (!options.pangene_gene_matches.empty()) {
        std::cout
            << "Pangene gene matches: "
            << cli::join_with_comma(options.pangene_gene_matches)
            << "\n";
    }
    if (!options.pangene_copy_tsv_path.empty()) {
        std::cout << "Pangene copy TSV: " << options.pangene_copy_tsv_path << "\n";
    }
    std::cout
        << "Region-level VCF: " << options.region_vcf_path << "\n"
        << "VCF records written: " << summary.region_vcf_records << "\n";
    if (options.write_debug_reports) {
        std::cout
            << "Call debug dir: " << options.debug_out_dir << "\n"
            << "Debug reports written: " << summary.debug_reports_written << "\n";
        if (!options.dotplot_gtf_path.empty()) {
            std::cout << "Dotplot GTF: " << options.dotplot_gtf_path << "\n";
        }
        if (!options.dotplot_gene_matches.empty()) {
            std::cout
                << "Dotplot gene matches: "
                << cli::join_with_comma(options.dotplot_gene_matches)
                << "\n";
        }
    }

    return 0;
}

} // namespace panvar
