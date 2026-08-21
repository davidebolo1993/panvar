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
        << "A cheap degree gate decides whether the graph is pathological (a cluster of hubs with a high\n"
        << "single-END degree, from seqwish over-merging); if so, the locus is rebuilt by progressive graph\n"
        << "generation using minigraph's engine, with haplotypes added most-complete-first by CANONICAL\n"
        << "k-mer richness (so the order does not depend on which strand the graph is stored on), ties\n"
        << "broken by path name. The emitted GFA carries per-haplotype P lines and preserves link\n"
        << "orientation, so inversion bubbles survive.\n"
        << "\n"
        << "The result is ACCEPTED only if every haplotype comes back: each path recovered as a walk whose\n"
        << "consecutive steps are joined by links that exist, with matched cover and recovered-walk\n"
        << "identity above --min-matched-cover / --min-recovered-identity. minigraph augments variation\n"
        << "above --min-var, so sub-threshold differences are collapsed by construction and a walk is\n"
        << "never byte-identical -- the contract is fidelity within those bounds, not losslessness.\n"
        << "Otherwise the rebuilt graph is discarded and the ORIGINAL written unchanged. Acceptance proves\n"
        << "fidelity, NOT that the graph got less tangled; the run reports hubs and max handle degree\n"
        << "before and after so that can be judged separately.\n"
        << "Healthy graphs pass through unchanged, so validated loci are never touched.\n\n"
        << "Options:\n"
        << "  -i, --gfa <path>            Input GFA (required)\n"
        << "  -o, --out <path>            Output GFA (required)\n"
        << "      --kmer <N>              k for the k-mer richness metric: haplotypes are ordered by\n"
        << "                              distinct k-mers, ties broken by total k-mers (default 21)\n"
        << "      --min-var <N>           Minimum variant length augmented into the graph, i.e.\n"
        << "                              minigraph -L (default 50). Lower values keep more of the\n"
        << "                              small variation, at the cost of a denser graph\n"
        << "      --min-align-len <N>     Minimum alignment length that may contribute events. 0 = auto:\n"
        << "                              scaled from the MEDIAN haplotype length and capped at half the\n"
        << "                              shortest, so no haplotype is gated out by a bar it could never\n"
        << "                              clear. A larger value is honoured but warned about.\n"
        << "      --tmp-dir <path>        Parent dir for the per-haplotype FASTA scratch (default: beside\n"
        << "                              --out); a dedicated subfolder under it is created and removed\n"
        << "      --hub-degree <N>        Node degree that counts as a hub (default 50)\n"
        << "      --min-hubs <N>          >= this many hubs => pathological (default 10)\n"
        << "  -t, --threads <N>           Worker threads (0 = auto)\n"
        << "  ---- acceptance contract (the rebuild is discarded and the ORIGINAL written if unmet) ----\n"
        << "      --min-recovered-identity <X>  Per-path identity of the walk re-spelled from the REBUILT\n"
        << "                              graph against the original haplotype (default 0.98)\n"
        << "      --min-matched-cover <X> Per-path matching bases / haplotype length (default 0.95).\n"
        << "                              Distinct from the chain's outer envelope, which hides internal gaps\n"
        << "  -r, --reference-path <name>  This path must be recovered within the same bounds. Exact match\n"
        << "                              wins, else a unique substring; an ambiguous name is refused.\n"
        << "                              Seeding stays richness-driven -- this is about RECOVERY\n"
        << "      --allow-loss            Accept a rebuild that fails the contract, recording what it violated\n"
        << "      --audit <path>          Per-path audit TSV (default <out>.rebuild_audit.tsv): coverage,\n"
        << "                              identity and a status per path, plus the global verdict\n"
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
        // Acceptance proves FIDELITY -- that every haplotype came back within the bounds. It does not
        // prove the graph got less tangled, which is the reason for rebuilding in the first place, so
        // report the before/after on the same per-handle measure and let the reader judge.
        // At least one metric must strictly improve and none may worsen. `<=` alone called an entirely
        // unchanged graph "untangled", which is the one thing the word must not mean.
        const bool no_worse = s.out_hubs <= s.raw_hubs && s.out_maxdeg <= s.raw_maxdeg &&
                              s.out_selfloops <= s.raw_selfloops;
        const bool strictly_better = s.out_hubs < s.raw_hubs || s.out_maxdeg < s.raw_maxdeg ||
                                     s.out_selfloops < s.raw_selfloops;
        const bool better = no_worse && strictly_better;
        log.info(std::string("structure: hubs ") + std::to_string(s.raw_hubs) + " -> " +
                 std::to_string(s.out_hubs) + ", max handle degree " + std::to_string(s.raw_maxdeg) +
                 " -> " + std::to_string(s.out_maxdeg) + ", self-loops " +
                 std::to_string(s.raw_selfloops) + " -> " + std::to_string(s.out_selfloops) +
                 (better ? " (untangled)"
                         : (no_worse ? " (unchanged -- the rebuild preserved the haplotypes but did not "
                                       "simplify the graph)"
                                     : " (NOT untangled -- some structure measure got worse)")));
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
