#include "panvar/refine_command.hpp"

#include "panvar/cli_utils.hpp"
#include "panvar/refine.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <algorithm>

namespace panvar {
namespace {

void print_refine_help() {
    std::cout
        << "Usage:\n"
        << "  panvar refine -i <graph.gfa> -b <prefix> -r <name> -o <prefix> [options]\n\n"
        << "Realign the haplotypes inside each site (abPOA), so a variant the graph spelled in pieces\n"
        << "becomes one record. Folded repeats are kept as they are and only the sequence around them is\n"
        << "realigned; a site holding an unfolded copy is skipped. Writes panphorte's file set, ready for call.\n\n"
        << "Options:\n"
        << "  -i, --gfa <path>                  Sorted GFA from panphorte (required)\n"
        << "  -b, --bubble-prefix <prefix>      Sites from <prefix>.bubbles.csv (required unless -c)\n"
        << "  -c, --bubbles-csv <path>          Site table, by path\n"
        << "      --panphorte-prefix <prefix>   Unit-graph table <prefix>.panphorte.unit_graph.tsv (default:\n"
        << "                                    the site table's prefix)\n"
        << "  -r, --reference-path <name>       Reference path, exact name or unique substring (required)\n"
        << "  -o, --out-prefix <prefix>         Output prefix (required)\n"
        << "      --gtf <path>                  Also write genes on reference nodes, <prefix>.bandage_genes.csv\n"
        << "      --bubble-id <id[,id...]>      Only these sites; repeatable (default: all)\n"
        << "      --max-poa-bp <N>              Skip a site whose longest sequence to align exceeds N bp\n"
        << "                                    (default 5000)\n"
        << "      --max-poa-work <N>            Skip a site whose longest x total bp to align exceeds N\n"
        << "                                    (default 0 = no cap)\n"
        << "      --max-walks <N>               Skip a site with more than N distinct sequences to align\n"
        << "                                    (default 500)\n"
        << "      --min-dup-bp <N>              Realign a site whose repeated node visits total under N bp\n"
        << "                                    instead of skipping it; realigning turns those copies into an\n"
        << "                                    indel (default 0 = any repeat skips)\n"
        << "      --min-bubbles <N>             Only realign regions of at least N adjacent sites (default 1)\n"
        << "      --resnarl-min-variant-bp <N>  bubble's --min-variant-bp for the sites found again (default 50)\n"
        << "      --merge-nearby-bp <N>         bubble's --merge-nearby-bp for the sites found again (default 0)\n"
        << "      --no-fuse-ambiguous           As in bubble, for the sites found again\n"
        << "      --no-flip                     Do not reorient nodes to the reference strand\n"
        << "  -q, --quiet                       No progress logs\n"
        << "  -h, --help                        Show this help\n";
}

// One --bubble-id value: decimal ids separated by commas. An empty value or item is refused, not
// dropped: an empty selection means "every site", which is not what an empty argument asked for.
std::vector<std::string> split_bubble_ids(const std::string& s) {
    std::vector<std::string> out;
    std::stringstream ss(s + ",");
    std::string item;
    while (std::getline(ss, item, ',')) {
        if (item.empty() || item.find_first_not_of("0123456789") != std::string::npos)
            throw std::runtime_error("refine: --bubble-id takes bubble ids separated by commas, got '" + s + "'");
        item.erase(0, std::min(item.find_first_not_of('0'), item.size() - 1));   // 007 is bubble 7
        out.push_back(item);
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
    bool panphorte_prefix_given = false;

    for (std::size_t i = 0; i < args.size(); ++i) {
        const std::string arg = cli::renamed_input_flag(args[i], "refine");
        auto require_value = [&](const std::string& flag) -> const std::string& {
            if (i + 1 >= args.size()) throw std::runtime_error("Missing value after " + flag);
            return args[++i];
        };
        if (arg == "-h" || arg == "--help") { print_refine_help(); return 0; }
        else if (arg == "-i" || arg == "--gfa") opt.gfa_path = require_value(arg);
        else if (arg == "-b" || arg == "--bubble-prefix") bubble_prefix = require_value(arg);
        else if (arg == "-c" || arg == "--bubbles-csv") opt.bubbles_csv_in = require_value(arg);
        else if (arg == "--panphorte-prefix") { opt.panphorte_prefix = require_value(arg); panphorte_prefix_given = true; }
        else if (arg == "-r" || arg == "--reference-path") opt.reference_path = require_value(arg);
        else if (arg == "-o" || arg == "--out-prefix") opt.out_prefix = require_value(arg);
        else if (arg == "--gtf") opt.gtf_path = require_value(arg);
        else if (arg == "--bubble-id") {
            for (std::string& id : split_bubble_ids(require_value(arg))) opt.only_bubble_ids.push_back(std::move(id));
        }
        else if (arg == "--max-poa-bp") opt.max_poa_bp = cli::parse_size_arg(arg, require_value(arg));
        else if (arg == "--max-poa-work") opt.max_poa_work = cli::parse_size_arg(arg, require_value(arg));
        else if (arg == "--min-dup-bp") opt.min_dup_bp = cli::parse_size_arg(arg, require_value(arg));
        else if (arg == "--merge-nearby-bp") opt.resnarl_merge_nearby_bp = cli::parse_size_arg(arg, require_value(arg));
        else if (arg == "--no-fuse-ambiguous") opt.resnarl_fuse_ambiguous = false;
        else if (arg == "--resnarl-min-variant-bp")
            opt.resnarl_min_variant_bp = cli::parse_size_arg(arg, require_value(arg));
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
            throw std::runtime_error("Missing bubbles input: --bubble-prefix <prefix> or --bubbles-csv <path>");
        }
        opt.bubbles_csv_in = bubble_prefix + ".bubbles.csv";
    }
    // -b and --bubbles-csv naming different runs would read the sites from one and the repeat tables
    // from the other, so the tables that go with the CSV must then be named explicitly.
    if (!bubble_prefix.empty() && opt.bubbles_csv_in != bubble_prefix + ".bubbles.csv" && !panphorte_prefix_given)
        throw std::runtime_error("refine: -b " + bubble_prefix + " and --bubbles-csv " + opt.bubbles_csv_in +
                                 " name different runs; give one of them, or --panphorte-prefix for the "
                                 "repeat tables that go with the CSV");
    if (!panphorte_prefix_given) {
        const std::string tail = ".bubbles.csv";
        const std::string& csv = opt.bubbles_csv_in;
        opt.panphorte_prefix = !bubble_prefix.empty() ? bubble_prefix
                             : csv.size() > tail.size() && csv.compare(csv.size() - tail.size(), tail.size(), tail) == 0
                                 ? csv.substr(0, csv.size() - tail.size()) : std::string();
    }
    if (panphorte_prefix_given && !std::ifstream(opt.panphorte_prefix + ".panphorte.unit_graph.tsv"))
        throw std::runtime_error("refine: --panphorte-prefix " + opt.panphorte_prefix + ": no " +
                                 opt.panphorte_prefix + ".panphorte.unit_graph.tsv");
    cli::RunLog log("refine", opt.quiet);
    log.info("reading: " + opt.gfa_path);
    const RefineSummary s = refine_graph(opt);
    for (const std::string& p : s.stale_removed) log.info("removed " + p + ", left by an earlier run at this prefix");
    std::vector<std::string> outs{opt.out_prefix + ".normalized.sorted.gfa",
                                  opt.out_prefix + ".bubbles.csv",
                                  opt.out_prefix + ".bandage_nodes.csv",
                                  opt.out_prefix + ".refine.report.tsv"};
    // The gene table is listed only when --gtf produced one.
    if (s.wrote_gene_annotation) outs.push_back(opt.out_prefix + ".bandage_genes.csv");
    outs.insert(outs.end(), s.rep_tables_written.begin(), s.rep_tables_written.end());
    log.wrote(outs);
    log.done();
    return 0;
}

}  // namespace panvar
