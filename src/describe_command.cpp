#include "panvar/describe_command.hpp"

#include "panvar/cli_utils.hpp"
#include "panvar/describe.hpp"

#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace panvar {
namespace {

DescribeFeatureMode parse_feature_mode_arg(const std::string& value) {
    if (value == "all" || value == "kmer" || value == "kmers") {
        return DescribeFeatureMode::AllKmers;
    }
    if (value == "syncmer" || value == "syncmers") {
        return DescribeFeatureMode::Syncmer;
    }
    throw std::runtime_error("--feature-mode must be one of: all, syncmer");
}

std::string feature_mode_label(DescribeFeatureMode mode) {
    switch (mode) {
        case DescribeFeatureMode::AllKmers:
            return "all";
        case DescribeFeatureMode::Syncmer:
            return "syncmer";
    }
    return "all";
}

void print_describe_help() {
    std::cout
        << "Usage:\n"
        << "  panvar describe -i <graph.gfa> (-b <prefix> | -c <bubbles.csv>) [-o <dir>] [options]\n\n"
        << "Options:\n"
        << "  -i, --gfa <path>                 panphorte-normalized/sorted GFA, i.e.\n"
        << "                                   <panphorte_prefix>.normalized.sorted.gfa (required)\n"
        << "  -b, --bubble-prefix-in <prefix>  panphorte output prefix\n"
        << "                                   (auto-uses <prefix>.bubbles.csv)\n"
        << "  -c, --bubbles-csv <path>         panphorte bubbles CSV (required if no prefix)\n"
        << "  -o, --out-dir <dir>              Output directory (default: describe_out)\n"
        << "      --bubble-id <N>              Restrict to one bubble ID (repeatable)\n"
        << "  -k, --kmer-size <K>              K-mer size, 1..31 for 2-bit encoding (default: 31)\n"
        << "      --feature-mode <mode>        all|syncmer (default: syncmer)\n"
        << "      --syncmer-s <S>              Internal s-mer size for closed syncmer mode (default: auto)\n"
        << "      --min-paths <N>              Drop features with min(present,absent) paths <= N,\n"
        << "                                   keeping copy-number features (default: 1; 0 keeps all)\n"
        << "      --max-wide-features <N>      Skip wide matrix above N features (default: 250000; 0=no cap)\n"
        << "      --force-wide                 Write wide matrix even above safety cap\n"
        << "      --no-wide-matrix             Write only feature map + sparse JSONL counts\n"
        << "      --variant-nodes <tsv>        Restrict k-mers to call <prefix>.variant_nodes.tsv\n"
        << "                                   (only those bubbles' variant nodes contribute)\n"
        << "      --variant-flank-bp <N>       With --variant-nodes, also keep nodes within N bp of a\n"
        << "                                   variant node so junction/flanking k-mers are retained\n"
        << "                                   (default: k-1 under --variant-nodes so short variant nodes\n"
        << "                                   still yield k-mers; pass 0 for strict variant-node-only)\n"
        << "      --samples <tsv>              cosigt sample->haplotype-path table; also writes the\n"
        << "                                   sample-level bimbam_{kmers,graph}.samples.bimbam.gz (summed dosage)\n"
        << "      --variant-vcf <vcf>          call region VCF; also emit the VARIANT-level BIMBAM\n"
        << "                                   (bimbam_variant.* + feature_annot.variant.tsv.gz): one\n"
        << "                                   dosage row per SV call -- the honest GWAS unit for associate\n"
        << "      --no-bimbam                  Do not write the pooled BIMBAM dosage + feature_annot.tsv.gz\n"
        << "      --only-kmers                 Emit only the k-mer substrate\n"
        << "      --only-graph                 Emit only the node/edge graph substrate\n"
        << "      --only-variant               Emit only the variant substrate (needs --variant-vcf; no GFA needed)\n"
        << "      --scale-dosage               Rescale each feature's BIMBAM dosage to 0..2 (per-feature min-max);\n"
        << "                                   for tools that assume a 0..2 diploid dosage (e.g. GEMMA). Linear-model\n"
        << "                                   p-values are unchanged; default off (raw counts)\n"
        << "      --threads <N>                Worker threads for the per-bubble loop (0 = auto)\n"
        << "  -q, --quiet                      Disable the progress bar\n"
        << "  -h, --help                       Show this help\n";
}

} // namespace

int run_describe_command(const std::vector<std::string>& args) {
    if (args.empty()) {
        print_describe_help();
        return 0;
    }

    std::string bubble_prefix_in;
    DescribeOptions options;
    bool variant_flank_set = false;

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
        if (arg == "-i" || arg == "--gfa") {
            options.gfa_path = require_value(arg);
            continue;
        }
        if (arg == "-b" || arg == "--bubble-prefix-in") {
            bubble_prefix_in = require_value(arg);
            continue;
        }
        if (arg == "-c" || arg == "--bubbles-csv") {
            options.bubbles_csv_in = require_value(arg);
            continue;
        }
        if (arg == "-o" || arg == "--out-dir") {
            options.out_dir = require_value(arg);
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
        if (arg == "-k" || arg == "--kmer-size") {
            options.kmer_size = cli::parse_size_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--feature-mode") {
            options.feature_mode = parse_feature_mode_arg(require_value(arg));
            continue;
        }
        if (arg == "--syncmer-s") {
            options.syncmer_s = cli::parse_size_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--min-paths") {
            options.min_feature_paths = cli::parse_size_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--max-wide-features") {
            options.max_wide_features = cli::parse_size_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--force-wide") {
            options.force_wide_matrix = true;
            continue;
        }
        if (arg == "--no-wide-matrix") {
            options.write_wide_matrix = false;
            continue;
        }
        if (arg == "--variant-nodes") {
            options.variant_nodes_path = require_value(arg);
            continue;
        }
        if (arg == "--variant-flank-bp") {
            options.variant_flank_bp = cli::parse_size_arg(arg, require_value(arg));
            variant_flank_set = true;
            continue;
        }
        if (arg == "--samples") {
            options.samples_path = require_value(arg);
            continue;
        }
        if (arg == "--variant-vcf") {
            options.variant_vcf_path = require_value(arg);
            continue;
        }
        if (arg == "--scale-dosage") {
            options.scale_dosage = true;
            continue;
        }
        if (arg == "--no-bimbam") {
            options.bimbam = false;
            continue;
        }
        if (arg == "--only-kmers") {
            options.emit_graph = false; options.emit_variant = false;
            continue;
        }
        if (arg == "--only-graph") {
            options.emit_kmers = false; options.emit_variant = false;
            continue;
        }
        if (arg == "--only-variant") {
            options.emit_kmers = false; options.emit_graph = false;
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
        if (arg == "--vcf-in" || arg == "--gtf" || arg == "--gene-match" || arg == "--size-bins") {
            throw std::runtime_error(
                "The describe module is now k-mer based and no longer consumes VCF/GTF inputs; "
                "use -i/--gfa with --bubble-prefix-in or --bubbles-csv.");
        }

        throw std::runtime_error("Unknown option for describe: " + arg);
    }

    // The k-mer and graph substrates are built from the GFA; the variant substrate is built only from
    // call's VCF, so --only-variant needs neither -i nor a bubble source.
    const bool graph_substrates = options.emit_kmers || options.emit_graph;
    if (graph_substrates && options.gfa_path.empty()) {
        throw std::runtime_error("describe requires -i/--gfa <graph.gfa>");
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
    if (graph_substrates && options.bubbles_csv_in.empty()) {
        throw std::runtime_error("describe requires --bubble-prefix-in <prefix> or --bubbles-csv <path>");
    }
    if (options.emit_variant && options.variant_vcf_path.empty() && !graph_substrates) {
        throw std::runtime_error("--only-variant requires --variant-vcf <call.region.vcf>");
    }
    if (options.kmer_size == 0 || options.kmer_size > 31) {
        throw std::runtime_error("--kmer-size must be in [1,31]");
    }
    // Under --variant-nodes, k-mers are confined to the kept variant runs. With zero flank a variant
    // node shorter than k yields NO k-mers (no k-window fits), so the k-mer layer comes out empty for
    // single short-variant bubbles. Default the flank to k-1 so junction-spanning k-mers are retained;
    // an explicit --variant-flank-bp (including 0 for strict variant-only) always wins.
    if (!options.variant_nodes_path.empty() && !variant_flank_set) {
        options.variant_flank_bp = options.kmer_size - 1;
    }
    if (options.feature_mode == DescribeFeatureMode::Syncmer) {
        const std::size_t syncmer_s =
            options.syncmer_s == 0 ? std::max<std::size_t>(1, std::min<std::size_t>(11, (options.kmer_size + 2) / 3))
                                   : options.syncmer_s;
        if (syncmer_s >= options.kmer_size) {
            throw std::runtime_error("--syncmer-s must be > 0 and < --kmer-size");
        }
    }
    if (options.force_wide_matrix && !options.write_wide_matrix) {
        throw std::runtime_error("--force-wide and --no-wide-matrix cannot be used together");
    }
    if (!options.emit_kmers && !options.emit_graph && !options.emit_variant) {
        throw std::runtime_error("no substrate selected (the --only-* flags are mutually exclusive)");
    }

    cli::RunLog log("describe", options.quiet);
    log.info("input " + options.gfa_path + " (feature mode " + feature_mode_label(options.feature_mode) +
             ", k=" + std::to_string(options.kmer_size) + ")");

    DescribeSummary summary;
    if (graph_substrates) {
        describe_kmers_from_graph(options, &summary);
    }
    if (options.emit_variant && !options.variant_vcf_path.empty()) {
        describe_variant_from_vcf(options, &summary);
    }

    log.info("processed " + std::to_string(summary.bubbles_processed) + " bubbles; kept " +
             std::to_string(summary.features_written) + "/" + std::to_string(summary.features_candidates) +
             " k-mer and " + std::to_string(summary.node_edge_features_written) + "/" +
             std::to_string(summary.node_edge_candidates) + " graph features");
    log.info("wrote " + std::to_string(summary.files_written) + " files in " + options.out_dir +
             "/ (index: describe.index.tsv)");
    log.done();
    return 0;
}

} // namespace panvar
