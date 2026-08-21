#include "panvar/describe_command.hpp"

#include <filesystem>
#include <system_error>
#include <unistd.h>

#include "panvar/cli_utils.hpp"
#include "panvar/describe.hpp"

#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
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
        << "                                   keeping copy-number features (default: 1; 0 keeps every\n"
        << "                                   DISCRIMINATIVE feature -- constant ones are always dropped)\n"
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

    std::string only_flag;
    bool flank_explicit = false;
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
            flank_explicit = true;
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
        // One --only-* selects one substrate. Applying each as "clear the other two" made the flags
        // order-dependent and let `--only-kmers --only-graph` disable every substrate, producing a
        // successful run that emitted nothing.
        if (arg == "--only-kmers" || arg == "--only-graph" || arg == "--only-variant") {
            if (!only_flag.empty() && only_flag != arg)
                throw std::runtime_error("describe: " + only_flag + " and " + arg +
                                         " are mutually exclusive; pass exactly one --only-* flag");
            only_flag = arg;
            options.emit_kmers = (arg == "--only-kmers");
            options.emit_graph = (arg == "--only-graph");
            options.emit_variant = (arg == "--only-variant");
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
    // Every BIMBAM, not just the pooled graph/k-mer ones. The variant substrate IS a BIMBAM and the
    // sample level is BIMBAM-only, so --no-bimbam with either of them asked for output it then
    // suppressed -- or, for --variant-vcf, silently wrote the matrix anyway.
    if (!options.bimbam && options.emit_variant && !options.variant_vcf_path.empty())
        throw std::runtime_error("--no-bimbam and --variant-vcf conflict: the variant substrate is a "
                                 "BIMBAM matrix, so there would be nothing to write");
    if (!options.bimbam && !options.samples_path.empty())
        throw std::runtime_error("--no-bimbam and --samples conflict: the sample level is emitted only "
                                 "as a BIMBAM matrix");
    // A flank with no variant restriction to widen is a silently inert argument.
    if (options.variant_flank_bp != 0 && options.variant_nodes_path.empty() && flank_explicit)
        throw std::runtime_error("--variant-flank-bp only means anything with --variant-nodes, which "
                                 "defines the scope it widens");
    // Compressed VCF input is not supported; the parsers read plain text.
    auto reject_gz = [](const std::string& path, const char* flag) {
        if (path.size() > 3 && path.compare(path.size() - 3, 3, ".gz") == 0)
            throw std::runtime_error(std::string(flag) + ": compressed VCF input is not supported; "
                                     "decompress it first");
    };
    reject_gz(options.variant_vcf_path, "--variant-vcf");
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

    // Everything describe owns is built in a sibling staging directory and moved into place only once
    // BOTH substrate passes have succeeded. Writing straight into the output directory meant a failure
    // part-way left a mixed family, and -- worse -- a rerun with different --only-*, --bubble-id or
    // --no-wide-matrix settings left the previous run's files sitting beside the new ones, looking
    // current. The transaction lives here rather than inside either pass because the two share one
    // directory: a per-pass commit would make the second treat the first's output as stale.
    const std::filesystem::path final_dir(options.out_dir);
    const std::filesystem::path stage_dir =
        final_dir.parent_path() /
        (final_dir.filename().string() + ".describe-staging." + std::to_string(::getpid()));
    std::error_code ec;
    std::filesystem::remove_all(stage_dir, ec);
    struct StageGuard {
        std::filesystem::path dir;
        bool committed = false;
        ~StageGuard() {
            if (!committed) { std::error_code e; std::filesystem::remove_all(dir, e); }
        }
    } guard{stage_dir, false};

    DescribeSummary summary;
    {
        DescribeOptions staged = options;
        staged.out_dir = stage_dir.string();
        if (graph_substrates) {
            describe_kmers_from_graph(staged, &summary);
        }
        if (staged.emit_variant && !staged.variant_vcf_path.empty()) {
            describe_variant_from_vcf(staged, &summary);
        }
    }

    // Commit. Two properties the first version did not have:
    //
    //   * OWNERSHIP IS EXACT. Anything starting with "bubble_" was treated as ours, so a user's
    //     `bubble_notes` file or directory would have been deleted. Only bubble_<digits> is generated.
    //   * NEITHER RESULT IS LOST. Deleting the old family and then moving the new one in entry by entry
    //     leaves neither complete if a move fails midway. The old entries are moved aside to a backup
    //     first and restored if anything goes wrong.
    auto is_owned = [](const std::string& n) {
        if (n == "describe.index.tsv" || n == "describe.params.json") return true;
        if (n == "haplotype" || n == "sample") return true;
        if (n.rfind("bubble_", 0) != 0) return false;
        const std::string suffix = n.substr(7);
        return !suffix.empty() &&
               suffix.find_first_not_of("0123456789") == std::string::npos;
    };
    std::filesystem::create_directories(final_dir);

    // An INPUT sitting under an owned output name would be consumed and then replaced by the commit
    // below. `reject_output_collisions` compares files; the owned set here is a set of directory
    // ENTRIES, some of them directories, so the containment test is done directly.
    {
        const auto canon = [](const std::string& p) {
            std::error_code e;
            const auto c = std::filesystem::weakly_canonical(p, e);
            return (e || c.empty()) ? std::filesystem::path(p) : c;
        };
        const std::filesystem::path canon_out = canon(final_dir.string());
        for (const std::string* in : {&options.gfa_path, &options.bubbles_csv_in,
                                      &options.variant_nodes_path, &options.samples_path,
                                      &options.variant_vcf_path}) {
            if (in->empty()) continue;
            std::filesystem::path c = canon(*in);
            for (std::filesystem::path up = c; !up.empty() && up != up.parent_path();
                 up = up.parent_path()) {
                if (up.parent_path() != canon_out) continue;
                if (!is_owned(up.filename().string())) break;
                throw std::runtime_error(
                    "describe: input '" + *in + "' lives under an output this run owns (" +
                    up.filename().string() + " in --out-dir); the commit would consume it and then "
                    "replace it");
            }
        }
    }

    const std::filesystem::path backup_dir =
        final_dir.parent_path() /
        (final_dir.filename().string() + ".describe-backup." + std::to_string(::getpid()));
    std::filesystem::remove_all(backup_dir, ec);
    std::vector<std::pair<std::filesystem::path, std::filesystem::path>> moved;   // backup <- original
    // Every destination this run installed, in order. Restoring only `moved` left an entry that had
    // NO predecessor sitting in the output directory after a rollback, so a failed run published half
    // a new family -- exactly the state the transaction exists to prevent.
    std::vector<std::filesystem::path> installed;
    std::size_t stale = 0;
    auto restore = [&]() {
        std::vector<std::string> unrestored;
        for (auto it = installed.rbegin(); it != installed.rend(); ++it) {
            std::error_code e;
            std::filesystem::remove_all(*it, e);
            if (e) unrestored.push_back("could not remove " + it->string() + ": " + e.message());
        }
        for (const auto& [bak, orig] : moved) {
            std::error_code e;
            std::filesystem::remove_all(orig, e);
            std::filesystem::rename(bak, orig, e);
            // A restore that itself fails is the one outcome the caller cannot recover from: the
            // previous output is gone and the backup is the only remaining copy. Naming it is the
            // difference between a recoverable situation and a silent one.
            if (e) unrestored.push_back("could not restore " + orig.string() + " (kept as " +
                                        bak.string() + "): " + e.message());
        }
        if (!unrestored.empty()) {
            std::string msg = "describe: rollback did not complete:";
            for (const auto& u : unrestored) msg += "\n  " + u;
            throw std::runtime_error(msg);
        }
        std::error_code e;
        std::filesystem::remove_all(backup_dir, e);
    };

    // Fault injection, test-only, mirroring the shared StagedOutputs contract: fail before installing
    // the Nth entry. Without it the rollback path is unreachable from a test.
    std::size_t fail_at = 0;
    if (const char* env = std::getenv("PANVAR_TEST_FAIL_COMMIT_AT")) {
        fail_at = static_cast<std::size_t>(std::strtoull(env, nullptr, 10));
    }
    std::size_t n_installed = 0;

    try {
        std::filesystem::create_directories(backup_dir);
        for (const auto& e : std::filesystem::directory_iterator(final_dir)) {
            if (!is_owned(e.path().filename().string())) continue;
            const std::filesystem::path bak = backup_dir / e.path().filename();
            std::filesystem::rename(e.path(), bak, ec);
            if (ec) throw std::runtime_error("describe: cannot set aside stale output " +
                                             e.path().string() + ": " + ec.message());
            moved.emplace_back(bak, e.path());
            ++stale;
        }
        for (const auto& e : std::filesystem::directory_iterator(stage_dir)) {
            const std::filesystem::path dst = final_dir / e.path().filename();
            if (fail_at != 0 && ++n_installed == fail_at) {
                throw std::runtime_error("PANVAR_TEST_FAIL_COMMIT_AT: injected failure installing " +
                                         dst.string());
            }
            std::filesystem::rename(e.path(), dst, ec);
            if (ec) {   // across filesystems rename fails; fall back to copy
                ec.clear();
                std::filesystem::copy(e.path(), dst,
                                      std::filesystem::copy_options::recursive |
                                      std::filesystem::copy_options::overwrite_existing, ec);
                if (ec) throw std::runtime_error("describe: cannot install output " + dst.string() +
                                                 ": " + ec.message());
            }
            installed.push_back(dst);
        }
    } catch (...) {
        restore();
        throw;
    }
    guard.committed = true;
    std::filesystem::remove_all(backup_dir, ec);
    std::filesystem::remove_all(stage_dir, ec);
    if (stale) log.info("replaced " + std::to_string(stale) + " stale describe output(s) from a previous run");

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
