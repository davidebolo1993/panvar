#include "panvar/rebuild_command.hpp"

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "panvar/cli_utils.hpp"
#include "panvar/rebuild.hpp"

#include <filesystem>

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
        else if (arg == "--min-recovered-identity")
            opt.min_recovered_identity = cli::parse_unit_fraction_arg(arg, require_value(arg));
        else if (arg == "--min-matched-cover")
            opt.min_matched_cover = cli::parse_unit_fraction_arg(arg, require_value(arg));
        else if (arg == "-r" || arg == "--reference-path") opt.reference_path = require_value(arg);
        else if (arg == "--allow-loss") opt.allow_loss = true;
        else if (arg == "--audit") opt.audit_path = require_value(arg);
        else if (arg == "-q" || arg == "--quiet") opt.quiet = true;
        else throw std::runtime_error("Unknown option: " + arg);
    }

    if (opt.gfa_path.empty()) throw std::runtime_error("Missing required input: --gfa <path>");
    if (opt.out_path.empty()) throw std::runtime_error("Missing required input: --out <path>");
    // Writing to the input path truncated it on open, so the graph was destroyed before it could be
    // read -- a healthy pass-through turned a 2 kB file into 0 bytes. Compare the resolved paths, and
    // also `equivalent` so a symlink or a second mount point to the same file is caught.
    {
        std::error_code ec;
        const std::filesystem::path in_p = std::filesystem::weakly_canonical(opt.gfa_path, ec);
        const std::filesystem::path out_p = std::filesystem::weakly_canonical(opt.out_path, ec);
        bool same = (!ec && !in_p.empty() && in_p == out_p);
        if (!same && std::filesystem::exists(opt.out_path)) {
            std::error_code eq;
            same = std::filesystem::equivalent(opt.gfa_path, opt.out_path, eq) && !eq;
        }
        if (same)
            throw std::runtime_error("rebuild: --out is the same file as --gfa (" + opt.out_path +
                                     "); refusing to overwrite the input");
    }
    cli::ensure_parent_dir_for_file(opt.out_path);

    cli::RunLog log("rebuild", opt.quiet);
    const RebuildSummary s = rebuild_graph(opt);
    // Report what was actually WRITTEN. Saying "rebuilt: 4740 -> 31 nodes" after the contract rejected
    // the rebuild and the original was passed through describes a file that does not exist.
    if (!s.ran) {
        log.info("passed through: " + std::to_string(s.raw_nodes) + " nodes (graph is not pathological)");
    } else if (s.accepted || opt.allow_loss) {
        log.info(std::string(s.accepted ? "rebuilt" : "rebuilt (contract not met, accepted with "
                                                     "--allow-loss)") + ": " +
                 std::to_string(s.raw_nodes) + " -> " + std::to_string(s.out_nodes) +
                 " nodes; seed=" + (s.seed.empty() ? "-" : s.seed) + "; " +
                 std::to_string(s.paths_recovered) + "/" + std::to_string(s.haplotypes) +
                 " paths recovered");
    } else {
        log.info("rejected: the rebuild would have been " + std::to_string(s.raw_nodes) + " -> " +
                 std::to_string(s.out_nodes) + " nodes, but " + s.reject_reason +
                 "; the original graph was written unchanged");
    }
    std::vector<std::string> written{opt.out_path};
    if (s.audit_written)
        written.push_back(opt.audit_path.empty() ? opt.out_path + ".rebuild_audit.tsv" : opt.audit_path);
    log.wrote(written);
    log.done();
    return 0;
}

} // namespace panvar
