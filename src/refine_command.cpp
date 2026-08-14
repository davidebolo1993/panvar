#include "panvar/refine_command.hpp"

#include "panvar/cli_utils.hpp"
#include "panvar/refine.hpp"

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace panvar {
namespace {

void print_refine_help() {
    std::cout
        << "Usage:\n"
        << "  panvar refine -i <panphorte.normalized.sorted.gfa> --bubble-prefix-in <prefix> \\\n"
        << "                -r <name> -o <prefix> [--gtf <gtf>] [options]\n\n"
        << "Re-align the actual haplotype sequences inside selected bubbles with POA (abPOA) so pggb\n"
        << "alignment artifacts `call` reports as split records (e.g. a spurious INS+DEL pair for one\n"
        << "indel) collapse into a single clean record. Folded self-loop REP nodes (DUPs) are held fixed\n"
        << "and only the residual flanks around them are re-aligned; bubbles carrying an unfolded\n"
        << "copy-number revisit are skipped. Emits panphorte's output family (sorted GFA + bubbles CSV +\n"
        << "Bandage colours [+ genes with --gtf]) so the result is a drop-in for `call`/`describe`.\n\n"
        << "Options:\n"
        << "  -i, --gfa <path>              Input panphorte normalized+sorted GFA (required)\n"
        << "  -b, --bubble-prefix-in <pfx>  Input bubbles from <pfx>.bubbles.csv (required unless --bubbles-csv-in)\n"
        << "      --bubbles-csv-in <path>   Explicit input bubbles CSV\n"
        << "  -r, --reference-path <name>   Reference path name or unique case-insensitive substring (required)\n"
        << "  -o, --out-prefix <prefix>     Output prefix (required)\n"
        << "      --gtf <path>              Reference GTF; also write <prefix>.bandage_genes.csv\n"
        << "      --bubble-id <id[,id...]>  Refine only these bubble ids (default: auto over the whole locus)\n"
        << "      --max-poa-bp <N>          Skip a residual segment whose LONGEST sequence exceeds this, and\n"
        << "                                bound total POA cost with it (default 5000)\n"
        << "      --resnarl-min-variant-bp <N>  Interior-span filter for the re-snarled CSV (default 50)\n"
        << "      --max-walks <N>           Skip a residual segment with more than N distinct walks (default 500)\n"
        << "      --min-bubbles <N>         Only rebuild regions fusing >= N bubbles (default 1)\n"
        << "      --no-flip                 Do not reorient nodes to the reference forward strand\n"
        << "  -q, --quiet                   Disable progress logs\n"
        << "  -h, --help                    Show this help\n";
}

std::vector<std::string> split_csv(const std::string& s) {
    std::vector<std::string> out;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, ',')) {
        if (!item.empty()) out.push_back(item);
    }
    return out;
}

}  // namespace

int run_refine_command(const std::vector<std::string>& args) {
    if (args.empty()) {
        print_refine_help();
        return 0;
    }
    RefineOptions opt;
    std::string bubble_prefix;

    for (std::size_t i = 0; i < args.size(); ++i) {
        const std::string& arg = args[i];
        auto require_value = [&](const std::string& flag) -> const std::string& {
            if (i + 1 >= args.size()) throw std::runtime_error("Missing value after " + flag);
            return args[++i];
        };
        if (arg == "-h" || arg == "--help") { print_refine_help(); return 0; }
        else if (arg == "-i" || arg == "--gfa") opt.gfa_path = require_value(arg);
        else if (arg == "-b" || arg == "--bubble-prefix-in") bubble_prefix = require_value(arg);
        else if (arg == "--bubbles-csv-in") opt.bubbles_csv_in = require_value(arg);
        else if (arg == "-r" || arg == "--reference-path") opt.reference_path = require_value(arg);
        else if (arg == "-o" || arg == "--out-prefix") opt.out_prefix = require_value(arg);
        else if (arg == "--gtf") opt.gtf_path = require_value(arg);
        else if (arg == "--bubble-id") opt.only_bubble_ids = split_csv(require_value(arg));
        else if (arg == "--max-poa-bp") opt.max_poa_bp = cli::parse_size_arg(arg, require_value(arg));
        else if (arg == "--max-walks") opt.max_walks = cli::parse_size_arg(arg, require_value(arg));
        else if (arg == "--min-bubbles") opt.min_bubbles = cli::parse_size_arg(arg, require_value(arg));
        else if (arg == "--no-flip") opt.no_flip = true;
        else if (arg == "-q" || arg == "--quiet") opt.quiet = true;
        else throw std::runtime_error("Unknown option: " + arg);
    }

    if (opt.gfa_path.empty()) throw std::runtime_error("Missing required input: --gfa <path>");
    if (opt.reference_path.empty()) throw std::runtime_error("Missing required input: --reference-path <name>");
    if (opt.out_prefix.empty()) throw std::runtime_error("Missing required input: --out-prefix <prefix>");
    if (opt.bubbles_csv_in.empty()) {
        if (bubble_prefix.empty()) {
            throw std::runtime_error("Missing bubbles input: --bubble-prefix-in <prefix> or --bubbles-csv-in <path>");
        }
        opt.bubbles_csv_in = bubble_prefix + ".bubbles.csv";
    }
    cli::ensure_parent_dir_for_file(opt.out_prefix + ".normalized.sorted.gfa");

    cli::RunLog log("refine", opt.quiet);
    const RefineSummary s = refine_graph(opt);
    log.info("rebuilt " + std::to_string(s.regions_rebuilt) + " region(s), skipped " +
             std::to_string(s.regions_skipped) + " (+" + std::to_string(s.nodes_added) + " nodes, -" +
             std::to_string(s.nodes_removed) + "); " + std::to_string(s.bubbles_after) + " bubbles after");
    std::vector<std::string> outs{opt.out_prefix + ".normalized.sorted.gfa",
                                  opt.out_prefix + ".bubbles.csv",
                                  opt.out_prefix + ".bandage_nodes.csv",
                                  opt.out_prefix + ".refine.report.tsv"};
    // The gene annotation was produced and never reported, so a --gtf run understated what it wrote.
    if (s.wrote_gene_annotation) outs.push_back(opt.out_prefix + ".bandage_genes.csv");
    log.wrote(outs);
    log.done();
    return 0;
}

}  // namespace panvar
