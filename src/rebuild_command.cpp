#include "panvar/rebuild_command.hpp"

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "panvar/cli_utils.hpp"
#include "panvar/rebuild.hpp"

namespace panvar {
namespace {

void print_rebuild_help() {
    std::cout
        << "Usage:\n"
        << "  panvar rebuild -i <graph.gfa> -o <out.gfa> [options]\n\n"
        << "Re-induce a pathological (fragmented low-complexity) locus graph BEFORE bubble decomposition.\n"
        << "A cheap degree gate decides whether the graph is pathological (a cluster of high-degree hubs\n"
        << "from seqwish over-merging); if so, the locus is rebuilt by progressive graph generation using\n"
        << "minigraph's engine, with haplotypes added most-complete-first (k-mer richness). The emitted GFA\n"
        << "carries per-haplotype P lines and preserves link orientation, so inversion bubbles survive.\n"
        << "Healthy graphs pass through unchanged, so validated loci are never touched.\n\n"
        << "Options:\n"
        << "  -i, --gfa <path>            Input GFA (required)\n"
        << "  -o, --out <path>            Output GFA (required)\n"
        << "      --kmer <N>              k for the k-mer richness metric: haplotypes are ordered by\n"
        << "                              distinct k-mers, ties broken by total k-mers (default 21)\n"
        << "      --min-var <N>           Minimum variant length augmented into the graph, i.e.\n"
        << "                              minigraph -L (default 50). Lower values keep more of the\n"
        << "                              small variation, at the cost of a denser graph\n"
        << "      --min-align-len <N>     Minimum alignment length that may contribute events, i.e.\n"
        << "                              minigraph -l (default auto: half the seed haplotype, since\n"
        << "                              minigraph's own default assumes chromosome-scale input)\n"
        << "      --tmp-dir <path>        Parent dir for the per-haplotype FASTA scratch (default: beside\n"
        << "                              --out); a dedicated subfolder under it is created and removed\n"
        << "      --hub-degree <N>        Node degree that counts as a hub (default 50)\n"
        << "      --min-hubs <N>          >= this many hubs => pathological (default 10)\n"
        << "  -t, --threads <N>           Worker threads (0 = auto)\n"
        << "      --force                 Rebuild even if the gate says healthy (testing/small inputs)\n"
        << "  -q, --quiet                 Disable progress logs\n"
        << "  -h, --help                  Show this help\n";
}

} // namespace

int run_rebuild_command(const std::vector<std::string>& args) {
    if (args.empty()) {
        print_rebuild_help();
        return 0;
    }
    RebuildOptions opt;

    for (std::size_t i = 0; i < args.size(); ++i) {
        const std::string& arg = args[i];
        auto require_value = [&](const std::string& flag) -> const std::string& {
            if (i + 1 >= args.size()) throw std::runtime_error("Missing value after " + flag);
            return args[++i];
        };
        if (arg == "-h" || arg == "--help") { print_rebuild_help(); return 0; }
        else if (arg == "-i" || arg == "--gfa") opt.gfa_path = require_value(arg);
        else if (arg == "-o" || arg == "--out") opt.out_path = require_value(arg);
        else if (arg == "--kmer") opt.kmer = cli::parse_size_arg(arg, require_value(arg));
        else if (arg == "--min-var") opt.min_var = cli::parse_size_arg(arg, require_value(arg));
        else if (arg == "--min-align-len") {
            const std::string& v = require_value(arg);
            opt.min_align_len = (v == "auto") ? 0 : cli::parse_size_arg(arg, v);
        }
        else if (arg == "--tmp-dir") opt.tmp_dir = require_value(arg);
        else if (arg == "--hub-degree") opt.hub_degree = cli::parse_size_arg(arg, require_value(arg));
        else if (arg == "--min-hubs") opt.min_hubs = cli::parse_size_arg(arg, require_value(arg));
        else if (arg == "-t" || arg == "--threads") opt.threads = cli::parse_size_arg(arg, require_value(arg));
        else if (arg == "--force") opt.force = true;
        else if (arg == "-q" || arg == "--quiet") opt.quiet = true;
        else throw std::runtime_error("Unknown option: " + arg);
    }

    if (opt.gfa_path.empty()) throw std::runtime_error("Missing required input: --gfa <path>");
    if (opt.out_path.empty()) throw std::runtime_error("Missing required input: --out <path>");
    cli::ensure_parent_dir_for_file(opt.out_path);

    cli::RunLog log("rebuild", opt.quiet);
    const RebuildSummary s = rebuild_graph(opt);
    log.info((s.ran ? "rebuilt" : "passed through") + std::string(": ") + std::to_string(s.raw_nodes) +
             " -> " + std::to_string(s.out_nodes) + " nodes; seed=" + (s.seed.empty() ? "-" : s.seed));
    log.wrote({opt.out_path});
    log.done();
    return 0;
}

} // namespace panvar
