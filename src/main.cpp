#include <algorithm>
#include <cctype>
#include <exception>
#include <filesystem>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <sstream>
#include <string>
#include <vector>

#include "panvar/allele.hpp"
#include "panvar/bubbles.hpp"
#include "panvar/call.hpp"
#include "panvar/describe.hpp"
#include "panvar/gfa.hpp"
#include "panvar/output.hpp"

namespace {

void print_general_help() {
    std::cout
        << "panvar - modular pangenome graph toolkit\n\n"
        << "Usage:\n"
        << "  panvar <subcommand> [options]\n\n"
        << "Subcommands:\n"
        << "  bubble    Module 1: refine/import sites from 'vg snarls'\n"
        << "  allele    Module 2: allele extraction and clustering from module-1 sites\n"
        << "  call      Module 3: variant calling on module-2 clustered alleles\n"
        << "  describe  Module 4: per-bubble haplotype feature description\n\n"
        << "Run 'panvar <subcommand> --help' for options.\n";
}

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
        << "  -h, --help                       Show this help\n";
}

void print_allele_help() {
    std::cout
        << "Usage:\n"
        << "  panvar allele -i <graph.gfa> [--bubble-prefix-in <module1-prefix> | --bubbles-csv-in <module1.bubbles.csv>] [options]\n\n"
        << "Options:\n"
        << "  -i, --gfa <path>                 Input GFA file (required)\n"
        << "  -o, --out-prefix <prefix>        Output prefix (default: allele_calls)\n"
        << "      --bubble-prefix-in <prefix>  Module-1 output prefix from 'panvar bubble'\n"
        << "                                   (auto-uses <prefix>.bubbles.csv)\n"
        << "      --bubbles-csv-in <path>      Module-1 bubbles CSV input (required if no prefix)\n"
        << "      --clusters-csv <path>        Explicit allele-cluster CSV output path\n"
        << "      --assignments-csv <path>     Explicit per-path assignment CSV output path\n"
        << "      --clusters-json <path>       Optional predefined path->cluster-label JSON map\n"
        << "      --cluster-sequences-csv <path> Optional representative sequence export\n"
        << "                                   (one representative sequence row per bubble+cluster)\n"
        << "      --min-similarity <X>         Min similarity in [0,1] or percent (e.g. 0.9 or 90)\n"
        << "                                    (default: 0.90)\n"
        << "      --cluster-mode <mode>        similarity mode: sequence|walk (default: sequence)\n"
        << "      --threads <N>                Worker threads for distance computation (0=auto)\n"
        << "      --distance-mode <mode>       auto|exact (default: auto)\n"
        << "      --max-upgma-alleles <N>      auto-switch to threshold-graph clustering above N unique\n"
        << "                                    alleles per bubble (default: 256, 0=disable)\n"
        << "      --quiet                      Disable progress logs\n"
        << "      --similarity-out-dir <dir>   Optional per-bubble similarity diagnostics folder\n"
        << "                                    (distance stats, pairwise cluster matrix,\n"
        << "                                     dendrogram SVG, UPGMA Newick)\n"
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

void print_describe_help() {
    std::cout
        << "Usage:\n"
        << "  panvar describe --vcf-in <call.region.vcf> --out-dir <dir> [options]\n\n"
        << "Options:\n"
        << "      --vcf-in <path>              Input region-level VCF from panvar call (required)\n"
        << "      --out-dir <dir>              Output directory for per-bubble tables\n"
        << "                                   (default: describe_out)\n"
        << "      --gtf <path>                 Optional GTF/GTF.GZ for gene-region overlaps\n"
        << "      --gene-match <expr>          Optional case-insensitive gene regex filter\n"
        << "                                   (repeat flag to pass multiple patterns)\n"
        << "      --size-bins <csv>            Size bins as comma-separated bp thresholds\n"
        << "                                   (default: 100,1000)\n"
        << "      --quiet                      Disable progress logs\n"
        << "  -h, --help                       Show this help\n";
}

std::string trim_ascii(const std::string& text) {
    std::size_t lo = 0;
    while (lo < text.size() && std::isspace(static_cast<unsigned char>(text[lo]))) {
        ++lo;
    }
    std::size_t hi = text.size();
    while (hi > lo && std::isspace(static_cast<unsigned char>(text[hi - 1]))) {
        --hi;
    }
    return text.substr(lo, hi - lo);
}

std::size_t parse_size_arg(const std::string& name, const std::string& value) {
    try {
        return static_cast<std::size_t>(std::stoull(value));
    } catch (const std::exception&) {
        throw std::runtime_error("Invalid value for " + name + ": " + value);
    }
}

double parse_similarity_arg(const std::string& name, const std::string& value) {
    double parsed = 0.0;
    try {
        parsed = std::stod(value);
    } catch (const std::exception&) {
        throw std::runtime_error("Invalid value for " + name + ": " + value);
    }
    if (parsed > 1.0 && parsed <= 100.0) {
        parsed /= 100.0;
    }
    if (!(parsed > 0.0 && parsed <= 1.0)) {
        throw std::runtime_error(name + " must be in (0,1] or (0,100]");
    }
    return parsed;
}

std::string normalize_minimap_preset(const std::string& value) {
    std::string lowered = value;
    std::transform(
        lowered.begin(),
        lowered.end(),
        lowered.begin(),
        [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    if (lowered == "asm5" || lowered == "asm10" || lowered == "asm20") {
        return lowered;
    }
    throw std::runtime_error(
        "Invalid value for --minimap-preset: " + value +
        " (supported: asm5, asm10, asm20)");
}

double parse_unit_fraction_arg(const std::string& name, const std::string& value) {
    double parsed = 0.0;
    try {
        parsed = std::stod(value);
    } catch (const std::exception&) {
        throw std::runtime_error("Invalid value for " + name + ": " + value);
    }
    if (!(parsed >= 0.0 && parsed <= 1.0)) {
        throw std::runtime_error(name + " must be in [0,1]");
    }
    return parsed;
}

bool parse_split_ins_svlen_geometric(const std::string& value) {
    std::string lowered = value;
    std::transform(
        lowered.begin(),
        lowered.end(),
        lowered.begin(),
        [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
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
    std::string lowered = value;
    std::transform(
        lowered.begin(),
        lowered.end(),
        lowered.begin(),
        [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    if (lowered == "strict" || lowered == "lenient") {
        return lowered;
    }
    throw std::runtime_error(
        "Invalid --vcf-merge-mode: " + value +
        " (expected strict|lenient)");
}

std::string join_with_comma(const std::vector<std::string>& values) {
    if (values.empty()) {
        return {};
    }
    std::string out;
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (i > 0) {
            out += ", ";
        }
        out += values[i];
    }
    return out;
}

void ensure_parent_dir_for_file(const std::string& path_text) {
    const std::filesystem::path p(path_text);
    const auto parent = p.parent_path();
    if (!parent.empty()) {
        std::filesystem::create_directories(parent);
    }
}

std::vector<std::size_t> parse_size_bins_arg(const std::string& csv) {
    std::vector<std::size_t> bins;
    std::string token;
    std::istringstream iss(csv);
    while (std::getline(iss, token, ',')) {
        const std::string trimmed = trim_ascii(token);
        if (trimmed.empty()) {
            continue;
        }
        const std::size_t v = parse_size_arg("--size-bins", trimmed);
        if (v == 0) {
            throw std::runtime_error("--size-bins values must be > 0");
        }
        bins.push_back(v);
    }
    if (bins.empty()) {
        throw std::runtime_error("--size-bins requires at least one integer threshold");
    }
    std::sort(bins.begin(), bins.end());
    bins.erase(std::unique(bins.begin(), bins.end()), bins.end());
    return bins;
}

panvar::ClusterMode parse_cluster_mode_arg(const std::string& value) {
    if (value == "sequence") {
        return panvar::ClusterMode::Sequence;
    }
    if (value == "sequence-fast" || value == "seq-fast" || value == "fast-sequence") {
        return panvar::ClusterMode::SequenceFast;
    }
    if (value == "walk" || value == "walk-signature" || value == "signature") {
        return panvar::ClusterMode::Walk;
    }
    throw std::runtime_error(
        "Invalid --cluster-mode: " + value +
        " (expected sequence|walk)");
}

const char* cluster_mode_label(panvar::ClusterMode mode) {
    switch (mode) {
        case panvar::ClusterMode::Sequence:
            return "sequence";
        case panvar::ClusterMode::SequenceFast:
            return "sequence";
        case panvar::ClusterMode::Walk:
            return "walk";
    }
    return "unknown";
}

int run_bubble(const std::vector<std::string>& args) {
    std::string gfa_path;
    std::string out_prefix = "bubble_calls";
    std::string bubbles_csv_path;
    std::string bandage_csv_path;
    std::string snarl_debug_tsv_path;

    panvar::BubbleCallOptions options;

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
            options.min_variant_bp = parse_size_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--min-path-support") {
            options.min_path_support = parse_size_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--merge-nearby-bp") {
            options.merge_nearby_bp = parse_size_arg(arg, require_value(arg));
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

    ensure_parent_dir_for_file(bubbles_csv_path);
    ensure_parent_dir_for_file(bandage_csv_path);
    if (!snarl_debug_tsv_path.empty()) {
        ensure_parent_dir_for_file(snarl_debug_tsv_path);
    }

    if (options.snarls_input_path.empty()) {
        throw std::runtime_error(
            "Missing required input: --snarls-in <path> (JSONL from 'vg view -R -j')");
    }

    panvar::ParseGfaOptions parse_options;
    parse_options.include_paths = true;
    parse_options.include_sequences = true;

    const panvar::Graph graph = panvar::parse_gfa(gfa_path, parse_options);
    if (graph.paths.empty()) {
        throw std::runtime_error("Input GFA has no P/W paths; snarl refinement requires path walks");
    }
    const auto report = panvar::call_bubbles_report(graph, options);
    const auto& bubbles = report.bubbles;

    panvar::write_bubbles_csv(bubbles_csv_path, bubbles);
    panvar::write_bandage_node_colors_csv(bandage_csv_path, bubbles, report.non_snp_bubbles);
    if (!snarl_debug_tsv_path.empty()) {
        panvar::write_snarl_debug_tsv(snarl_debug_tsv_path, report.snarl_debug);
    }

    std::size_t simple = 0;
    std::size_t super = 0;
    std::size_t insertion = 0;
    std::size_t inversion_signal_count = 0;
    std::size_t long_path_positive_count = 0;
    for (const auto& bubble : bubbles) {
        if (bubble.inversion_signal) {
            ++inversion_signal_count;
        }
        if (bubble.long_path_support > 0) {
            ++long_path_positive_count;
        }
        switch (bubble.type) {
            case panvar::BubbleType::Simple:
                ++simple;
                break;
            case panvar::BubbleType::Insertion:
                ++insertion;
                break;
            case panvar::BubbleType::Super:
                ++super;
                break;
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
        << "Bubble topology counts: simple=" << simple
        << ", super=" << super
        << ", insertion=" << insertion << "\n"
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

int run_allele(const std::vector<std::string>& args) {
    std::string gfa_path;
    std::string out_prefix = "allele_calls";
    std::string bubble_prefix_in;
    std::string clusters_csv_path;
    std::string assignments_csv_path;

    panvar::AlleleCallOptions options;

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
        if (arg == "--bubbles-csv-in") {
            options.bubbles_csv_in = require_value(arg);
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
            options.min_similarity = parse_similarity_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--cluster-mode") {
            options.cluster_mode = parse_cluster_mode_arg(require_value(arg));
            continue;
        }
        if (arg == "--threads") {
            options.threads = parse_size_arg(arg, require_value(arg));
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
            options.max_upgma_alleles = parse_size_arg(arg, require_value(arg));
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
            options.flank_nodes = parse_size_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--max-paths-per-bubble") {
            options.max_paths_per_bubble = parse_size_arg(arg, require_value(arg));
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
            options.odgi_width = parse_size_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--odgi-path-height") {
            options.odgi_path_height = parse_size_arg(arg, require_value(arg));
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
    if (options.bubbles_csv_in.empty()) {
        throw std::runtime_error(
            "Module 'allele' requires module-1 bubbles: use --bubble-prefix-in <prefix> "
            "or --bubbles-csv-in <path>.");
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
    if (options.cluster_mode == panvar::ClusterMode::SequenceFast) {
        options.fast_distance = true;
    }

    if (clusters_csv_path.empty()) {
        clusters_csv_path = out_prefix + ".allele_clusters.csv";
    }
    if (assignments_csv_path.empty()) {
        assignments_csv_path = out_prefix + ".allele_assignments.csv";
    }
    ensure_parent_dir_for_file(clusters_csv_path);
    ensure_parent_dir_for_file(assignments_csv_path);
    if (options.write_cluster_sequences) {
        ensure_parent_dir_for_file(options.cluster_sequences_csv_path);
    }

    panvar::ParseGfaOptions parse_options;
    parse_options.include_paths = true;
    parse_options.include_sequences = true;

    const panvar::Graph graph = panvar::parse_gfa(gfa_path, parse_options);
    if (graph.paths.empty()) {
        throw std::runtime_error("Input GFA has no P/W paths; allele clustering requires paths");
    }
    panvar::AlleleCallSummary summary;
    panvar::call_alleles_to_csv(graph, options, clusters_csv_path, assignments_csv_path, &summary);

    std::cout
        << "Input graph: " << gfa_path << '\n'
        << "Bubble source: " << options.bubbles_csv_in << "\n"
        << "Nodes: " << graph.nodes.size() << "\n"
        << "P/W paths loaded: " << graph.paths.size() << "\n"
        << "Min similarity: " << options.min_similarity << "\n"
        << "Cluster mode: " << cluster_mode_label(options.cluster_mode) << "\n"
        << "Distance mode: " << (options.fast_distance ? "auto" : "exact") << "\n"
        << "Predefined clusters: "
        << (options.predefined_clusters_json_path.empty() ? "off" : options.predefined_clusters_json_path)
        << "\n"
        << "Max UPGMA alleles: " << options.max_upgma_alleles << " (0=disabled)\n"
        << "Threads: " << options.threads << " (0=auto)\n"
        << "Bubbles processed: " << summary.bubbles_processed << "\n"
        << "Bubbles with haplotype assignments: " << summary.bubbles_with_assignments << "\n"
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

int run_call(const std::vector<std::string>& args) {
    std::string gfa_path;
    std::string out_prefix = "call_calls";
    std::string bubble_prefix_in;
    std::string allele_prefix_in;
    std::string clusters_csv_in_path;
    std::string assignments_csv_in_path;

    panvar::AlleleCallOptions options;
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
        if (arg == "--min-sv-bp") {
            options.min_sv_bp = parse_size_arg(arg, require_value(arg));
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
            options.minimap_best_n = parse_size_arg(arg, require_value(arg));
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
            options.vcf_merge_window_bp = parse_size_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--vcf-merge-mode") {
            options.vcf_merge_mode = normalize_vcf_merge_mode(require_value(arg));
            continue;
        }
        if (arg == "--vcf-merge-lenient-window-bp") {
            options.vcf_merge_lenient_window_bp = parse_size_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--vcf-merge-lenient-min-ref-jaccard") {
            options.vcf_merge_lenient_min_ref_jaccard =
                parse_unit_fraction_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--vcf-merge-min-seq-sim") {
            options.vcf_merge_min_seq_similarity =
                parse_unit_fraction_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--vcf-merge-max-edit-frac") {
            options.vcf_merge_max_seq_edit_fraction =
                parse_unit_fraction_arg(arg, require_value(arg));
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
    ensure_parent_dir_for_file(options.region_vcf_path);
    if (options.write_debug_reports && options.debug_out_dir.empty()) {
        options.debug_out_dir = out_prefix + ".debug";
    }
    if (!options.pangene_bed_path.empty() && options.pangene_copy_tsv_path.empty()) {
        options.pangene_copy_tsv_path = out_prefix + ".pangene_copy.tsv";
    }
    if (!options.pangene_copy_tsv_path.empty()) {
        ensure_parent_dir_for_file(options.pangene_copy_tsv_path);
    }
    panvar::ParseGfaOptions parse_options;
    parse_options.include_paths = true;
    parse_options.include_sequences = true;
    const panvar::Graph graph = panvar::parse_gfa(gfa_path, parse_options);

    if (graph.paths.empty()) {
        throw std::runtime_error("Input GFA has no P/W paths; variant calling requires paths");
    }

    panvar::AlleleCallSummary summary;
    panvar::call_variants_from_precomputed_alleles(
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
            << join_with_comma(options.pangene_gene_matches)
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
                << join_with_comma(options.dotplot_gene_matches)
                << "\n";
        }
    }

    return 0;
}

int run_describe(const std::vector<std::string>& args) {
    std::string vcf_in_path;
    std::string out_dir = "describe_out";

    panvar::DescribeOptions options;

    for (std::size_t i = 0; i < args.size(); ++i) {
        const std::string& arg = args[i];
        auto require_value = [&](const std::string& flag) -> const std::string& {
            if (i + 1 >= args.size()) {
                throw std::runtime_error("Missing value after " + flag);
            }
            return args[++i];
        };

        if (arg == "-h" || arg == "--help") {
            print_describe_help();
            return 0;
        }
        if (arg == "--vcf-in") {
            vcf_in_path = require_value(arg);
            continue;
        }
        if (arg == "--out-dir") {
            out_dir = require_value(arg);
            continue;
        }
        if (arg == "--gtf") {
            options.gtf_path = require_value(arg);
            continue;
        }
        if (arg == "--gene-match") {
            options.gene_match_patterns.push_back(require_value(arg));
            continue;
        }
        if (arg == "--size-bins") {
            options.size_bins = parse_size_bins_arg(require_value(arg));
            continue;
        }
        if (arg == "--quiet") {
            options.quiet = true;
            continue;
        }

        throw std::runtime_error("Unknown option for describe: " + arg);
    }

    if (vcf_in_path.empty()) {
        throw std::runtime_error("describe requires --vcf-in");
    }

    options.vcf_in_path = vcf_in_path;
    options.out_dir = out_dir;

    panvar::DescribeSummary summary;
    panvar::describe_from_region_vcf(options, &summary);

    std::cout
        << "Input VCF: " << options.vcf_in_path << "\n"
        << "Output dir: " << options.out_dir << "\n"
        << "Bubbles described: " << summary.bubbles << "\n"
        << "Events processed: " << summary.events << "\n"
        << "Haplotype rows written: " << summary.haplotype_rows << "\n"
        << "Files written: " << summary.files_written << "\n";

    if (!options.gtf_path.empty()) {
        std::cout << "GTF: " << options.gtf_path << "\n";
    }
    if (!options.gene_match_patterns.empty()) {
        std::cout << "Gene filters: " << join_with_comma(options.gene_match_patterns) << "\n";
    }

    return 0;
}

} // namespace

int main(int argc, char** argv) {
    try {
        if (argc <= 1) {
            print_general_help();
            return 1;
        }

        const std::string subcommand = argv[1];

        std::vector<std::string> args;
        args.reserve(static_cast<std::size_t>(argc));
        for (int i = 2; i < argc; ++i) {
            args.emplace_back(argv[i]);
        }

        if (subcommand == "-h" || subcommand == "--help") {
            print_general_help();
            return 0;
        }

        if (subcommand == "bubble") {
            return run_bubble(args);
        }
        if (subcommand == "allele") {
            return run_allele(args);
        }
        if (subcommand == "call") {
            return run_call(args);
        }
        if (subcommand == "describe") {
            return run_describe(args);
        }

        throw std::runtime_error("Unknown subcommand: " + subcommand);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
        return 2;
    }
}
