#include "panvar/rebuild_command.hpp"

#include <iostream>
#include <limits>
#include <utility>
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
        << "Re-induce a fragmented locus graph before bubble decomposition.\n\n"
        << "Options:\n"
        << "  -i, --gfa <path>               Input GFA, a regular file (required)\n"
        << "  -o, --out <path>               Output GFA, ending in .gfa or .gfa.gz (gzipped) (required). Also\n"
        << "                                 writes <prefix>.rebuild.report.tsv when a rebuilt graph is evaluated\n"
        << "  -r, --reference-path <name>    Reference path: it seeds the graph and must come back exactly\n"
        << "                                 (ignoring case), or the rebuild is refused. Exact name or\n"
        << "                                 unique substring\n"
        << "      --mode <folded|minigraph>  folded: each copy-number module becomes one unit subgraph its\n"
        << "                                 copies walk (default); minigraph: the graph minigraph builds\n"
        << "      --fold-identity <X>        Identity of a repeat that makes a module, 0.80 to 1 (default 0.95)\n"
        << "      --fold-min-copy-bp <N>     Shortest copy that makes a module, 1000 or more (default 3000;\n"
        << "                                 shorter tandem units are left to panphorte)\n"
        << "      --min-var <N>              Shortest variant added to the graph (minigraph -L; default 50)\n"
        << "      --min-align-len <N>        Shortest alignment that may add variation (minigraph -l).\n"
        << "                                 0 = auto (default): minigraph's 100 kb, or half the median\n"
        << "                                 haplotype if shorter (at least 1 kb); at most half the shortest\n"
        << "      --base-align               Base-level alignment in generation (minigraph -c; up to 2x slower)\n"
        << "      --clip-to-reference        With -r: drop haplotype sequence running past the reference's ends\n"
        << "      --kmer <N>                 k-mer size for ordering haplotypes by k-mer richness, 11 to 255\n"
        << "                                 (default 21)\n"
        << "      --hub-degree <N>           Node degree that counts as a hub (default 50)\n"
        << "      --min-hubs <N>             Hubs that make a graph fragmented (default 10)\n"
        << "      --force                    Rebuild even a graph that is not fragmented\n"
        << "      --tmp-dir <path>           Parent directory for scratch files (default: beside --out)\n"
        << "  -t, --threads <N>              Worker threads (default 0 = all cores)\n"
        << "  -q, --quiet                    Disable progress logs\n"
        << "  -h, --help                     Show this help\n\n"
        << "Acceptance (when unmet, the input graph is written unchanged):\n"
        << "      --min-recovered-identity <X>  Identity of each rebuilt walk to its haplotype (default 0.98)\n"
        << "      --min-matched-cover <X>    Share of each haplotype's bases aligned and matching (default 0.95)\n"
        << "      --allow-loss               Accept a rebuild that fails these, recording why; never one whose\n"
        << "                                 reference does not come back exactly\n";
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
        else if (arg == "--fold-identity") {
            opt.fold_identity = cli::parse_unit_fraction_arg(arg, require_value(arg));
            if (opt.fold_identity < 0.80) throw std::runtime_error("--fold-identity must be 0.80 or more");
        }
        else if (arg == "--fold-min-copy-bp") {
            opt.fold_min_copy_bp = cli::parse_size_arg(arg, require_value(arg));
            if (opt.fold_min_copy_bp < 1000) throw std::runtime_error("--fold-min-copy-bp must be 1000 or more");
        }
        else if (arg == "--min-align-len") {
            const std::string& v = require_value(arg);
            opt.min_align_len = (v == "auto") ? 0 : cli::parse_size_arg(arg, v);
        }
        else if (arg == "--tmp-dir") opt.tmp_dir = require_value(arg);
        else if (arg == "--hub-degree" || arg == "--min-hubs") {
            // 0 makes every graph pathological, which forces a rebuild without --force's warning
            // that variation under --min-var is collapsed.
            const std::size_t v = cli::parse_size_arg(arg, require_value(arg));
            if (v == 0) throw std::runtime_error(arg + " must be 1 or more");
            (arg == "--hub-degree" ? opt.hub_degree : opt.min_hubs) = v;
        }
        else if (arg == "-t" || arg == "--threads") opt.threads = cli::parse_size_arg(arg, require_value(arg));
        else if (arg == "--force") opt.force = true;
        else if (arg == "--base-align") opt.base_align = true;
        else if (arg == "--clip-to-reference") opt.clip_to_reference = true;
        else if (arg == "--min-recovered-identity")
            opt.min_recovered_identity = cli::parse_unit_fraction_arg(arg, require_value(arg));
        else if (arg == "--min-matched-cover")
            opt.min_matched_cover = cli::parse_unit_fraction_arg(arg, require_value(arg));
        else if (arg == "-r" || arg == "--reference-path") {
            opt.reference_path = require_value(arg);
            if (opt.reference_path.empty()) throw std::runtime_error(arg + " is empty: name the reference path");
        }
        else if (arg == "--allow-loss") opt.allow_loss = true;
        else if (arg == "--mode") {
            const std::string m = require_value(arg);
            if (m == "folded") opt.fold_modules = true;
            else if (m == "minigraph") opt.fold_modules = false;
            else throw std::runtime_error("--mode is '" + m + "', expected folded or minigraph");
        }
        else if (arg == "--audit")
            throw std::runtime_error("--audit was removed: the report is always written beside the output, "
                                     "as <prefix>.rebuild.report.tsv");
        else if (arg == "-q" || arg == "--quiet") opt.quiet = true;
        else throw std::runtime_error("Unknown option: " + arg);
    }

    if (opt.gfa_path.empty()) throw std::runtime_error("Missing required input: --gfa <path>");
    if (opt.out_path.empty()) throw std::runtime_error("Missing required input: --out <path>");
    if (opt.min_var < 1) throw std::runtime_error("--min-var must be at least 1");
    // under 3 kb a tandem repeat is panphorte's to fold
    if (opt.fold_min_copy_bp < 3000)
        std::cerr << cli::log_prefix("rebuild") << "warning: --fold-min-copy-bp " << opt.fold_min_copy_bp
                  << ": modules of copies under 3 kb are usually tandem repeats left to panphorte\n";
    if (opt.clip_to_reference && opt.reference_path.empty())
        throw std::runtime_error("--clip-to-reference needs -r: the reference's ends define what is clipped");
    // Below 11 the k-mer space (4^k) is smaller than a locus, so distinct counts saturate and stop ranking
    // anything; 255 is only a sanity limit.
    if (opt.kmer < 11 || opt.kmer > 255)
        throw std::runtime_error("--kmer must be between 11 and 255 (got " + std::to_string(opt.kmer) + ")");
    // These reach minigraph and the thread pool as int: values that would not fit are refused.
    const std::pair<const char*, std::size_t> as_int[] = {
        {"--min-var", opt.min_var}, {"--min-align-len", opt.min_align_len}, {"--threads", opt.threads}};
    for (const auto& [flag, v] : as_int)
        if (v > static_cast<std::size_t>(std::numeric_limits<int>::max()))
            throw std::runtime_error(std::string(flag) + " is too large: " + std::to_string(v));
    // Writing to the input path would truncate it on open, before it is read. Compare the resolved paths, and
    // `equivalent` so a symlink or a second mount point to the same file is caught; the report gets the same
    // check.
    auto same_file = [](const std::string& a, const std::string& b) {
        std::error_code ec;
        const std::filesystem::path pa = std::filesystem::weakly_canonical(a, ec);
        const std::filesystem::path pb = std::filesystem::weakly_canonical(b, ec);
        bool same = (!ec && !pa.empty() && pa == pb);
        if (!same && std::filesystem::exists(b)) {
            std::error_code eq;
            same = std::filesystem::equivalent(a, b, eq) && !eq;
        }
        return same;
    };
    // The input is read twice -- parsed, and copied as it stands on a pass-through -- so it must be a regular
    // file, not a pipe.
    cli::require_regular_input("rebuild", "--gfa", opt.gfa_path);
    // The report is named from the output's prefix, so the output must say it is a GFA.
    {
        const auto ends = [&](const std::string& e) {
            return opt.out_path.size() > e.size() && opt.out_path.compare(opt.out_path.size() - e.size(), e.size(), e) == 0;
        };
        if (!ends(".gfa") && !ends(".gfa.gz"))
            throw std::runtime_error("rebuild: --out must end in .gfa or .gfa.gz (got '" + opt.out_path + "')");
    }
    if (same_file(opt.gfa_path, opt.out_path))
        throw std::runtime_error("rebuild: --out is the same file as --gfa (" + opt.out_path +
                                 "); refusing to overwrite the input");
    if (same_file(opt.gfa_path, rebuild_report_path(opt.out_path)))
        throw std::runtime_error("rebuild: the report path " + rebuild_report_path(opt.out_path) +
                                 " is the input graph; refusing to overwrite it");
    // Directories made for the output are removed again if the run fails, when nothing else is in them.
    std::vector<std::filesystem::path> made;
    for (std::filesystem::path d = std::filesystem::path(opt.out_path).parent_path();
         !d.empty() && !std::filesystem::exists(d); d = d.parent_path())
        made.push_back(d);
    cli::ensure_parent_dir_for_file(opt.out_path);

    cli::RunLog log("rebuild", opt.quiet);
    RebuildSummary s;
    try {
        s = rebuild_graph(opt);
    } catch (...) {
        std::error_code ec;
        for (const std::filesystem::path& d : made) std::filesystem::remove(d, ec);   // only if empty
        throw;
    }
    // The run logs its verdict and both graphs' shapes; here only what was written.
    std::vector<std::string> written{opt.out_path};
    if (s.report_written) written.push_back(rebuild_report_path(opt.out_path));
    log.wrote(written);
    log.done();
    return 0;
}

} // namespace panvar
