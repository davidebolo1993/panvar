#include "panvar/panphorte_command.hpp"

#include "panvar/cli_utils.hpp"
#include "panvar/panphorte.hpp"

#include <filesystem>

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace panvar {
namespace {

void print_panphorte_help() {
    std::cout
        << "Usage:\n"
        << "  panvar panphorte -i <graph.gfa> (-b <prefix> | -c <bubbles.csv>) -o <prefix> [options]\n\n"
        << "Normalizes tandem-repeat bubbles into a compact representation (repeat unit + cycle)\n"
        << "and writes a new GFA that can be re-processed by bubble/inspect.\n\n"
        << "Options:\n"
        << "  -i, --gfa <path>                 Input GFA file (required)\n"
        << "  -b, --bubble-prefix-in <prefix>  Module-1 output prefix (auto-uses <prefix>.bubbles.csv)\n"
        << "  -c, --bubbles-csv <path>         Module-1 bubbles CSV (required if no prefix)\n"
        << "  -o, --out-prefix <prefix>        Output prefix for the normalized GFA + report (required)\n"
        << "      --bubble-id <N>              Restrict to one bubble ID (repeatable)\n"
        << "      --min-unit-bp <N>            Minimum repeat-unit span to normalize (default: 50)\n"
        << "      --allow-partial-boundary     Fold a site even when a copy cannot be bounded by any\n"
        << "                                   step range at all. Since copies keep their flanks as\n"
        << "                                   fragment nodes this no longer arises on any measured\n"
        << "                                   input; the guard and this override remain for the case\n"
        << "                                   where a haplotype would otherwise reach the site\n"
        << "                                   literally and be read CN 0 while carrying copies\n"
        << "      --resnarl-min-variant-bp <N> Interior-span filter for the re-snarled call-ready CSV\n"
        << "                                   under --reference-path (default: 50, 0 = keep all)\n"
        << "      --min-copies <N>             Minimum tandem copies to normalize (default: 2)\n"
        << "      --min-array-prevalence <f>   Min fraction of traversing haplotypes carrying a >=min-copies\n"
        << "                                   array for a bubble to be folded -- separates true population\n"
        << "                                   VNTRs (folded) from rare private gene/module dups left for\n"
        << "                                   `call` (e.g. CYP2D6x2) (default: 0.5)\n"
        << "      --max-interruption-frac <f>  Tolerance for interrupting bases in an array (default: 0.25)\n"
        << "      --min-similarity <f>         Min identity to treat a block as a unit copy\n"
        << "                                   (1.0=exact; <1.0 enables approximate, lossy collapse; default: 1.0)\n"
        << "      --threads <N>                Worker threads for approximate detection (0=auto; default: 0)\n"
        << "  -r, --reference-path <name>      If set, internally sort+flip the normalized graph along\n"
        << "                                   this reference and re-snarl (cactus), writing only\n"
        << "                                   <prefix>.normalized.sorted.gfa + <prefix>.bubbles.csv +\n"
        << "                                   <prefix>.bandage_nodes.csv so 'call' runs with no\n"
        << "                                   external vg/odgi tools (the unsorted GFA is not written)\n"
        << "      --no-flip                    With --reference-path, skip reorienting to the ref strand\n"
        << "      --gtf <path>                 Reference-coordinate GTF; after re-sorting, project genes\n"
        << "                                   onto the normalized graph's nodes and write\n"
        << "                                   <prefix>.bandage_genes.csv (needs a PanSN --reference-path)\n"
        << "  -q, --quiet                      Disable progress logs\n"
        << "  -h, --help                       Show this help\n";
}

} // namespace

int run_panphorte_command(const std::vector<std::string>& args) {
    if (args.empty()) {
        print_panphorte_help();
        return 0;
    }

    std::string bubble_prefix_in;
    std::string gtf_path;
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
        if (arg == "-b" || arg == "--bubble-prefix-in") {
            bubble_prefix_in = require_value(arg);
            continue;
        }
        if (arg == "-c" || arg == "--bubbles-csv") {
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
        if (arg == "--allow-partial-boundary") {
            options.allow_partial_boundary = true;
            continue;
        }
        if (arg == "--resnarl-min-variant-bp") {
            options.resnarl_min_variant_bp = cli::parse_size_arg(arg, require_value(arg));
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
        if (arg == "--min-array-prevalence") {
            options.min_array_prevalence = cli::parse_unit_fraction_arg(arg, require_value(arg));
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
        if (arg == "-r" || arg == "--reference-path") {
            options.reference_path = require_value(arg);
            continue;
        }
        if (arg == "--no-flip") {
            options.no_flip = true;
            continue;
        }
        if (arg == "--gtf") {
            gtf_path = require_value(arg);
            continue;
        }
        if (arg == "-q" || arg == "--quiet") {
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

    cli::RunLog log("panphorte", options.quiet);
    log.info("input " + options.gfa_path);

    // The gene annotation is written next to the other outputs, so it can name the --gtf input.
    if (!gtf_path.empty()) {
        std::error_code e1, e2;
        const auto gp = std::filesystem::weakly_canonical(gtf_path, e1);
        const auto op = std::filesystem::weakly_canonical(options.out_prefix + ".bandage_genes.csv", e2);
        if (!e1 && !e2 && !gp.empty() && gp == op)
            throw std::runtime_error("panphorte: output '" + options.out_prefix +
                                     ".bandage_genes.csv' is the same file as input '" + gtf_path + "'");
    }

    // Gene annotation on the collapsed graph: node ids differ from the bubble graph, so this Bandage
    // CSV is distinct from the one bubble emits. Needs a PanSN --reference-path. Projected inside
    // panphorte_normalize so it shares the staged transaction rather than landing after it.
    options.gtf_path = options.reference_path.empty() ? std::string() : gtf_path;
    if (!gtf_path.empty() && options.reference_path.empty()) {
        // Projecting genes needs reference coordinates, so without --reference-path there is nothing to
        // project onto. It was skipped in silence, which reads as "the GTF had no genes here".
        log.info("--gtf needs --reference-path to place genes; no gene annotation will be written");
    }

    PanphorteSummary summary;
    panphorte_normalize(options, &summary);

    log.info("normalized " + std::to_string(summary.bubbles_normalized) + "/" +
             std::to_string(summary.bubbles_seen) + " bubbles (" +
             std::to_string(summary.paths_rewritten) + " paths rewritten; +" +
             std::to_string(summary.nodes_added) + " (" +
             std::to_string(summary.fragment_nodes_added) + " fragments) −" +
             std::to_string(summary.nodes_removed) + " nodes, +" +
             std::to_string(summary.edges_added) + " −" + std::to_string(summary.edges_removed) +
             " edges)");

    std::vector<std::string> outputs;
    if (summary.sorted) {
        log.info("sorted along the reference; re-snarled " +
                 std::to_string(summary.resnarled_bubbles) + " bubbles");
        outputs.push_back(options.out_prefix + ".normalized.sorted.gfa");
        outputs.push_back(options.out_prefix + ".bubbles.csv");
        outputs.push_back(options.out_prefix + ".bandage_nodes.csv");
        if (summary.genes_written) {
            outputs.push_back(options.out_prefix + ".bandage_genes.csv");
        }
    } else {
        outputs.push_back(options.out_prefix + ".normalized.gfa");
    }
    outputs.push_back(options.out_prefix + ".panphorte.report.tsv");
    // Both were written on every run and neither was named. The provenance table in particular is the
    // join key between a REP node and the site it stands for, so a consumer has to be told it exists.
    outputs.push_back(options.out_prefix + ".panphorte.rep_provenance.tsv");
    if (options.min_similarity < 1.0) {
        outputs.push_back(options.out_prefix + ".panphorte.copies.tsv");
    }

    log.wrote(outputs);
    log.done();
    return 0;
}

} // namespace panvar
