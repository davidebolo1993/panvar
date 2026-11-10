#include "panvar/panphorte_command.hpp"

#include "panvar/cli_utils.hpp"
#include "panvar/panphorte.hpp"

#include <fstream>
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
        << "Fold tandem-repeat sites onto their repeat unit, so copy number is the number of passes through it.\n\n"
        << "Options:\n"
        << "  -i, --gfa <path>                 Sorted GFA from bubble (required)\n"
        << "  -b, --bubble-prefix <prefix>     bubble's output prefix (reads <prefix>.bubbles.csv)\n"
        << "  -c, --bubbles-csv <path>         bubble's site table, instead of -b\n"
        << "  -o, --out-prefix <prefix>        Output prefix (required)\n"
        << "  -r, --reference-path <name>      Sort the folded graph along this reference and find its sites\n"
        << "                                   again, ready for call (writes <prefix>.normalized.sorted.gfa and\n"
        << "                                   <prefix>.bubbles.csv). Without -r, <prefix>.normalized.gfa,\n"
        << "                                   for viewing: call needs -r\n"
        << "      --bubble-id <N>              Only this site (repeatable)\n"
        << "      --min-similarity <f>         Identity for a block to count as a copy of the unit (default 1.0,\n"
        << "                                   identical copies only; below 1, copies are found by alignment and\n"
        << "                                   folded onto a unit graph; a value above 1 is a percentage)\n"
        << "      --lossy-units                With --min-similarity below 1, keep only the reference's\n"
        << "                                   differences between copies; other haplotypes keep their copy\n"
        << "                                   number but not their own differences\n"
        << "      --min-unit-bp <N>            Shortest unit folded, at least 1 (default 50; with --min-similarity\n"
        << "                                   below 1, units under 16 bp are never folded)\n"
        << "      --min-copies <N>             Copies a haplotype needs for the site to be an array (default 2)\n"
        << "      --min-array-prevalence <f>   Share of the site's haplotypes that must carry an array for it to be\n"
        << "                                   folded; rarer duplications are left to call (default 0.5)\n"
        << "      --max-interruption-frac <f>  Share of an array's bases that may lie between copies (default 0.25)\n"
        << "      --refine                     Realign the flanks around each folded array (needs -r; writes\n"
        << "                                   <prefix>.panphorte.refine.tsv)\n"
        << "      --refine-max-poa-bp <N>      Longest flank realigned (default 5000)\n"
        << "      --refine-max-walks <N>       Most distinct flanks per side (default 500)\n"
        << "      --refine-max-poa-work <N>    Longest x total flank bp per site (default 0 = no limit)\n"
        << "      --resnarl-min-variant-bp <N> With -r: bubble's --min-variant-bp for the new sites (default 50)\n"
        << "      --merge-nearby-bp <N>        With -r: bubble's --merge-nearby-bp (default 0 = off, so a folded\n"
        << "                                   site is not merged with its neighbours)\n"
        << "      --no-fuse-ambiguous          With -r: as in bubble\n"
        << "      --no-flip                    With -r: do not reorient nodes to the reference strand\n"
        << "      --gtf <path>                 With -r: also write genes on reference nodes,\n"
        << "                                   <prefix>.bandage_genes.csv\n"
        << "  -t, --threads <N>                Worker threads (default 0 = all cores)\n"
        << "  -q, --quiet                      No progress logs\n"
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
        const std::string arg = cli::renamed_input_flag(args[i], "panphorte");
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
        if (arg == "-b" || arg == "--bubble-prefix") {
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
        if (arg == "--merge-nearby-bp") { options.resnarl_merge_nearby_bp = cli::parse_size_arg(arg, require_value(arg)); continue; }
        if (arg == "--no-fuse-ambiguous") { options.resnarl_fuse_ambiguous = false; continue; }
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
        if (arg == "--lossy-units") {
            options.lossless_units = false;
            continue;
        }
        if (arg == "-t" || arg == "--threads") {
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
        if (arg == "--refine") {
            options.refine = true;
            continue;
        }
        if (arg == "--refine-max-poa-bp") {
            options.refine_max_poa_bp = cli::parse_size_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--refine-max-walks") {
            options.refine_max_walks = cli::parse_size_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--refine-max-poa-work") {
            options.refine_max_poa_work = cli::parse_size_arg(arg, require_value(arg));
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
                "Conflicting bubble inputs: --bubble-prefix resolves to '" +
                derived + "' but --bubbles-csv is '" + options.bubbles_csv_in + "'");
        }
    }
    if (options.bubbles_csv_in.empty()) {
        throw std::runtime_error("panphorte requires --bubble-prefix <prefix> or --bubbles-csv <path>");
    }
    if (options.refine && options.reference_path.empty()) {
        throw std::runtime_error("--refine needs -r/--reference-path");
    }
    if (options.min_unit_bp == 0) {
        throw std::runtime_error("--min-unit-bp must be at least 1");
    }
    if (!gtf_path.empty()) {
        // Before the fold, not after it.
        cli::require_regular_input("panphorte", "--gtf", gtf_path);
        if (!std::ifstream(gtf_path)) throw std::runtime_error("panphorte: cannot read --gtf " + gtf_path);
    }
    if (options.min_copies < 2) {
        throw std::runtime_error("--min-copies must be >= 2");
    }
    if (!options.lossless_units && options.min_similarity >= 1.0) {
        throw std::runtime_error("--lossy-units needs --min-similarity below 1 (exact mode folds only identical "
                                 "copies, so there is nothing to discard)");
    }

    cli::RunLog log("panphorte", options.quiet);
    log.info("reading: " + options.gfa_path);
    if (options.reference_path.empty()) {
        // The unit-graph table names this run's node ids; sorting the output elsewhere renumbers them.
        std::cerr << cli::log_prefix("panphorte") << "warning: without -r/--reference-path the folded graph is "
                     "not sorted and its sites are not found again; it is for viewing, and call needs a run with -r\n";
    }
    if (!options.lossless_units && options.reference_path.empty()) {
        // Without a reference no path is kept exact: every haplotype walks the plain unit. Shown even with -q.
        std::cerr << cli::log_prefix("panphorte") << "warning: --lossy-units without -r/--reference-path keeps no "
                     "path's sequence: every haplotype, a path named like the reference included, walks the plain "
                     "unit once per copy\n";
    }

    // Gene annotation on the folded graph (node ids differ from the bubble graph's), projected inside
    // panphorte_normalize so it shares the staged transaction. Passed through even when unused without a
    // reference, so the input/output alias preflight still sees it as an input.
    options.gtf_path = gtf_path;
    if (!gtf_path.empty() && options.reference_path.empty()) {
        // Genes need reference coordinates to be placed; say so rather than skip in silence, which would
        // read as "the GTF had no genes here".
        std::cerr << cli::log_prefix("panphorte")
                  << "warning: --gtf needs --reference-path to place genes; no gene annotation will be written\n";
    }

    PanphorteSummary summary;
    panphorte_normalize(options, &summary);

    std::vector<std::string> outputs;
    if (summary.sorted) {
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
    // Named so a consumer knows they exist: the unit graph table is the join key between a fold's nodes
    // and the site it stands for.
    outputs.push_back(options.out_prefix + ".panphorte.unit_graph.tsv");
    if (options.min_similarity < 1.0) {
        outputs.push_back(options.out_prefix + ".panphorte.copies.tsv");
    }
    if (options.refine) {
        outputs.push_back(options.out_prefix + ".panphorte.refine.tsv");
    }

    log.wrote(outputs);
    log.done();
    return 0;
}

} // namespace panvar
