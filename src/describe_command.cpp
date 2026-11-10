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
#include <cstdlib>
#include <iterator>

namespace panvar {
namespace {

DescribeFeatureMode parse_feature_mode_arg(const std::string& value) {
    if (value == "all") {
        return DescribeFeatureMode::AllKmers;
    }
    if (value == "syncmer") {
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
        << "Turn each haplotype's walks through the sites into association features: k-mers, node and\n"
        << "edge traversal counts, and one dosage per call record (per ALT on a record with several).\n"
        << "Each is written as a BIMBAM matrix for associate. A copy number stays a count, not a\n"
        << "presence flag.\n\n"
        << "Options:\n"
        << "  -i, --gfa <path>                  Sorted GFA the calls were made on (required unless\n"
        << "                                    --only-variant)\n"
        << "  -b, --bubble-prefix <prefix>      Sites from <prefix>.bubbles.csv (required unless -c or\n"
        << "                                    --only-variant)\n"
        << "  -c, --bubbles-csv <path>          Site table, by path\n"
        << "  -o, --out-dir <dir>               Output directory (default describe_out)\n"
        << "  -r, --reference-path <name>       Reference path, exact name or unique substring: k-mer and\n"
        << "                                    graph rows get its coordinates (records have the VCF's)\n"
        << "      --bubble-id <N>               Only this site; repeatable (default: all)\n"
        << "      --variant-vcf <vcf>           call's region VCF, uncompressed and per haplotype (not\n"
        << "                                    --diploid): also write the variant matrix\n"
        << "      --variant-nodes <tsv>         call's <prefix>.variant_nodes.tsv: k-mers and graph features\n"
        << "                                    only from its sites' variant nodes and boundaries; the\n"
        << "                                    variant rows' nodes\n"
        << "      --variant-flank-bp <N>        With --variant-nodes, also N bases along the walk on either\n"
        << "                                    side for k-mers, every node they reach for graph features\n"
        << "                                    (default k-1)\n"
        << "      --samples <tsv>               cosigt table (sample, then its haplotype paths): also write\n"
        << "                                    per-sample matrices, summed over each sample's haplotypes\n"
        << "      --only-kmers                  Only the k-mer matrix (the --only-* flags combine; default:\n"
        << "                                    all three, the variant one with --variant-vcf)\n"
        << "      --only-graph                  Only the node and edge matrix\n"
        << "      --only-variant                Only the variant matrix (needs --variant-vcf, not -i)\n"
        << "  -k, --kmer-size <K>               K-mer size, 1 to 31 (default 31)\n"
        << "      --feature-mode <mode>         syncmer (closed syncmers only) or all (default syncmer)\n"
        << "      --syncmer-s <S>               S-mer size of the closed syncmers (default (k+2)/3, at most 11)\n"
        << "      --min-paths <N>               Drop a feature present in, or absent from, N paths or fewer,\n"
        << "                                    unless its count varies (default 1; 0 drops only constant\n"
        << "                                    features)\n"
        << "      --scale-dosage                Rescale each feature to 0..2, for tools that expect a diploid\n"
        << "                                    dosage (default: raw counts)\n"
        << "      --no-bimbam                   No BIMBAM matrices, only the per-site files\n"
        << "      --no-wide-matrix              Per site, only the feature map and sparse counts\n"
        << "      --max-wide-features <N>       Skip a site's dense matrix above N features (default 250000;\n"
        << "                                    0 = no cap)\n"
        << "      --force-wide                  Write the dense matrix above that cap\n"
        << "  -t, --threads <N>                 Worker threads (default 0 = all cores)\n"
        << "  -q, --quiet                       No progress logs\n"
        << "  -h, --help                        Show this help\n";
}

} // namespace

int run_describe_command(const std::vector<std::string>& args) {
    if (args.empty()) {
        print_describe_help();
        return 0;
    }

    std::string bubble_prefix_in;
    DescribeOptions options;

    std::string only_flags;   // the --only-* flags given, for messages ("" = none: every substrate)
    bool flank_explicit = false;
    bool syncmer_explicit = false;
    for (std::size_t i = 0; i < args.size(); ++i) {
        const std::string arg = cli::renamed_input_flag(args[i], "describe");
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
        if (arg == "-b" || arg == "--bubble-prefix") {
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
            // 0 is the unset value inside the options, so an explicit 0 is refused here.
            if (options.syncmer_s == 0) throw std::runtime_error("--syncmer-s must be > 0 and < --kmer-size");
            syncmer_explicit = true;
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
        if (arg == "-r" || arg == "--reference-path") {
            options.reference_path = require_value(arg);
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
        // The --only-* flags combine: the first clears the default (every substrate) and each adds its own.
        if (arg == "--only-kmers" || arg == "--only-graph" || arg == "--only-variant") {
            if (only_flags.empty()) options.emit_kmers = options.emit_graph = options.emit_variant = false;
            if (only_flags.find(arg) == std::string::npos) only_flags += (only_flags.empty() ? "" : " ") + arg;
            if (arg == "--only-kmers") options.emit_kmers = true;
            if (arg == "--only-graph") options.emit_graph = true;
            if (arg == "--only-variant") options.emit_variant = true;
            continue;
        }
        if (arg == "-t" || arg == "--threads") {
            options.threads = cli::parse_size_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "-q" || arg == "--quiet") {
            options.quiet = true;
            continue;
        }
        throw std::runtime_error("Unknown option for describe: " + arg);
    }

    // The k-mer and graph substrates come from the GFA, the variant one from call's VCF alone.
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
                "Conflicting bubble inputs: --bubble-prefix resolves to '" +
                derived + "' but --bubbles-csv is '" + options.bubbles_csv_in + "'");
        }
    }
    if (graph_substrates && options.bubbles_csv_in.empty()) {
        throw std::runtime_error("describe requires --bubble-prefix <prefix> or --bubbles-csv <path>");
    }
    // Asked for by name, the variant substrate needs a VCF; by default a run without one emits the other
    // two and says so.
    bool variant_skipped = false;
    if (options.emit_variant && options.variant_vcf_path.empty()) {
        if (!only_flags.empty())
            throw std::runtime_error("--only-variant requires --variant-vcf <call.region.vcf>");
        options.emit_variant = false;
        variant_skipped = true;
    }
    // The variant substrate and the sample level are BIMBAM matrices only.
    if (!options.bimbam && options.emit_variant && !options.variant_vcf_path.empty())
        throw std::runtime_error("--no-bimbam and --variant-vcf conflict: the variant substrate is a "
                                 "BIMBAM matrix, so there would be nothing to write");
    if (!options.bimbam && !options.samples_path.empty())
        throw std::runtime_error("--no-bimbam and --samples conflict: the sample level is emitted only "
                                 "as a BIMBAM matrix");
    // An input or setting for output the run does not write is refused.
    if (!options.emit_variant && !options.variant_vcf_path.empty())
        throw std::runtime_error("--variant-vcf with " + only_flags + ": the variant substrate it feeds is not "
                                 "emitted; add --only-variant or drop --variant-vcf");
    if (!graph_substrates && flank_explicit)
        throw std::runtime_error("--variant-flank-bp with --only-variant: it widens the k-mer and graph "
                                 "scope, which is not written");
    if (options.variant_nodes_path.empty() && flank_explicit)
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
    // A variant node shorter than k holds no k-mer, so the flank defaults to k-1: the k-mers across its
    // junctions are kept.
    if (!options.variant_nodes_path.empty() && !flank_explicit) {
        options.variant_flank_bp = options.kmer_size - 1;
    }
    if (options.feature_mode == DescribeFeatureMode::Syncmer) {
        const std::size_t syncmer_s =
            options.syncmer_s == 0 ? std::max<std::size_t>(1, std::min<std::size_t>(11, (options.kmer_size + 2) / 3))
                                   : options.syncmer_s;
        if (options.kmer_size < 2)
            throw std::runtime_error("--feature-mode syncmer needs --kmer-size 2 or more (the s-mer is "
                                     "shorter than the k-mer); use --feature-mode all for k=1");
        if (syncmer_s >= options.kmer_size) {
            throw std::runtime_error("--syncmer-s must be > 0 and < --kmer-size");
        }
    }
    if (syncmer_explicit && options.feature_mode != DescribeFeatureMode::Syncmer)
        throw std::runtime_error("--syncmer-s only applies with --feature-mode syncmer");
    if (options.scale_dosage && !options.bimbam)
        throw std::runtime_error("--scale-dosage and --no-bimbam conflict: the scaling applies to the BIMBAM "
                                 "matrices, which would not be written");
    if (options.force_wide_matrix && !options.write_wide_matrix) {
        throw std::runtime_error("--force-wide and --no-wide-matrix cannot be used together");
    }
    if (!options.emit_kmers && !options.emit_graph && !options.emit_variant) {
        throw std::runtime_error("no substrate selected");
    }

    // A trailing slash would leave the directory no name to stage beside.
    while (options.out_dir.size() > 1 && options.out_dir.back() == '/') options.out_dir.pop_back();
    // An input both passes read must be a file they can each open: a pipe gives its content to the first.
    if (graph_substrates && options.emit_variant) {
        if (!options.samples_path.empty()) cli::require_regular_input("describe", "--samples", options.samples_path);
        if (!options.variant_nodes_path.empty())
            cli::require_regular_input("describe", "--variant-nodes", options.variant_nodes_path);
    }

    // describe's own files in the output directory, known by name: the index, and the files of haplotype/,
    // sample/ and each bubble_<N>/. A rerun replaces them and removes those it did not write, as the other
    // modules do with the files at their prefix; anything else in the directory is left alone.
    const auto is_owned_file = [](const std::filesystem::path& rel) {
        std::vector<std::string> parts;
        for (const auto& part : rel) parts.push_back(part.string());
        if (parts.size() == 1) return parts[0] == "describe.index.tsv";
        if (parts.size() == 2) {
            const std::string& d = parts[0];
            if (d.size() <= 7 || d.rfind("bubble_", 0) != 0 || d[7] == '0' ||
                d.find_first_not_of("0123456789", 7) != std::string::npos)
                return false;
            for (const char* f : {"kmer_features.tsv.gz", "kmer_counts.jsonl.gz", "graph_features.tsv.gz",
                                  "kmer_matrix.tsv.gz", "graph_matrix.tsv.gz"})
                if (parts[1] == f) return true;
            return false;
        }
        if (parts.size() == 3) {
            const std::string& sub = parts[1];
            if ((parts[0] != "haplotype" && parts[0] != "sample") ||
                (sub != "kmers" && sub != "graph" && sub != "variant"))
                return false;
            return parts[2] == "samples.txt.gz" || parts[2] == "feature_annot." + sub + ".tsv.gz" ||
                   parts[2] == "bimbam_" + sub + ".bimbam.gz";
        }
        return false;
    };

    cli::RunLog log("describe", options.quiet);
    if (graph_substrates) log.info("reading: " + options.gfa_path);
    if (variant_skipped)
        log.info("variant matrix: skipped, it needs call's VCF (--variant-vcf <call>.region.vcf); "
                 "--only-kmers --only-graph selects the other two explicitly");

    // The output directory and the parents it lacks are made through StagedOutputs, which removes them
    // again, when left empty, if the run fails, and holds them against a concurrent run's cleanup until
    // then. The files themselves are installed below, so nothing is committed through it.
    const std::filesystem::path final_dir(options.out_dir);
    cli::StagedOutputs out_dirs("describe");
    out_dirs.stage((final_dir / "describe.index.tsv").string());

    // Both passes write into a sibling staging directory, installed only once both succeed, so a failure
    // leaves the previous output whole.
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

    // Commit, file by file: each file this run wrote replaces its namesake, each one replaced set aside
    // first, so a failure puts them back and removes what was installed, leaving the previous output whole.
    // Then describe's files this run did not write are removed, and folders left empty with them.
    std::filesystem::create_directories(final_dir);
    std::vector<std::filesystem::path> written;   // relative to the output directory
    for (std::filesystem::recursive_directory_iterator it(stage_dir), end; it != end; ++it)
        if (it->is_regular_file()) written.push_back(std::filesystem::relative(it->path(), stage_dir));
    std::sort(written.begin(), written.end());

    // An input that is one of describe's files in the output directory would be replaced or removed.
    const auto canon = [](const std::filesystem::path& p) {
        std::error_code e;
        const auto c = std::filesystem::weakly_canonical(p, e);
        return (e || c.empty()) ? p : c;
    };
    const std::filesystem::path canon_out = canon(final_dir);
    std::vector<std::filesystem::path> inputs;
    for (const std::string* in : {&options.gfa_path, &options.bubbles_csv_in, &options.variant_nodes_path,
                                  &options.samples_path, &options.variant_vcf_path}) {
        if (in->empty()) continue;
        const std::filesystem::path c = canon(*in);
        inputs.push_back(c);
        const std::filesystem::path rel = c.lexically_relative(canon_out);
        if (!rel.empty() && *rel.begin() != ".." && is_owned_file(rel))
            throw std::runtime_error("describe: input '" + *in + "' is one of describe's own files in "
                                     "--out-dir, which this run would replace or remove");
    }

    const std::filesystem::path backup_dir =
        final_dir.parent_path() /
        (final_dir.filename().string() + ".describe-backup." + std::to_string(::getpid()));
    std::filesystem::remove_all(backup_dir, ec);
    std::vector<std::pair<std::filesystem::path, std::filesystem::path>> moved;   // backup <- original
    std::vector<std::filesystem::path> installed;
    auto restore = [&]() {
        std::vector<std::string> unrestored;
        for (auto it = installed.rbegin(); it != installed.rend(); ++it) {
            std::error_code e;
            std::filesystem::remove(*it, e);
            if (e) unrestored.push_back("could not remove " + it->string() + ": " + e.message());
        }
        for (const auto& [bak, orig] : moved) {
            std::error_code e;
            std::filesystem::rename(bak, orig, e);
            // A failed restore names the backup, then the only copy of the previous file.
            if (e) unrestored.push_back("could not restore " + orig.string() + " (kept as " +
                                        bak.string() + "): " + e.message());
        }
        // Folders the install made and left empty.
        for (auto it = installed.rbegin(); it != installed.rend(); ++it)
            for (std::filesystem::path d = it->parent_path(); d != final_dir && !d.empty(); d = d.parent_path()) {
                std::error_code e;
                if (!std::filesystem::is_empty(d, e) || e) break;
                std::filesystem::remove(d, e);
            }
        if (!unrestored.empty()) {
            std::string msg = "describe: rollback did not complete:";
            for (const auto& u : unrestored) msg += "\n  " + u;
            throw std::runtime_error(msg);
        }
        std::error_code e;
        std::filesystem::remove_all(backup_dir, e);
    };

    // Test-only fault injection: fail before installing the Nth file.
    std::size_t fail_at = 0;
    if (const char* env = std::getenv("PANVAR_TEST_FAIL_COMMIT_AT")) {
        fail_at = static_cast<std::size_t>(std::strtoull(env, nullptr, 10));
    }
    try {
        std::size_t n_installed = 0;
        for (const std::filesystem::path& rel : written) {
            const std::filesystem::path dst = final_dir / rel;
            if (fail_at != 0 && ++n_installed == fail_at)
                throw std::runtime_error("PANVAR_TEST_FAIL_COMMIT_AT: injected failure installing " + dst.string());
            const auto st = std::filesystem::symlink_status(dst, ec);
            if (std::filesystem::is_directory(st))
                throw std::runtime_error("describe: output '" + dst.string() + "' is an existing directory");
            if (std::filesystem::exists(st)) {
                const std::filesystem::path bak = backup_dir / rel;
                std::filesystem::create_directories(bak.parent_path());
                std::filesystem::rename(dst, bak, ec);
                if (ec) throw std::runtime_error("describe: cannot set aside " + dst.string() + ": " + ec.message());
                moved.emplace_back(bak, dst);
            }
            std::filesystem::create_directories(dst.parent_path());
            // Registered before the move, so a failed move is rolled back too.
            installed.push_back(dst);
            // Staging is a sibling of the output directory, so the rename never crosses file systems.
            std::filesystem::rename(stage_dir / rel, dst, ec);
            if (ec) throw std::runtime_error("describe: cannot install output " + dst.string() + ": " + ec.message());
        }
    } catch (...) {
        restore();
        throw;
    }
    guard.committed = true;
    std::filesystem::remove_all(backup_dir, ec);
    std::filesystem::remove_all(stage_dir, ec);

    // describe's files an earlier run left that this run did not write. Only the top level and describe's
    // own folders are read, so an output directory holding other data is not walked.
    std::size_t removed = 0;
    {
        const std::vector<std::filesystem::path> wrote_set = written;
        const auto wrote = [&](const std::filesystem::path& rel) {
            return std::binary_search(wrote_set.begin(), wrote_set.end(), rel);
        };
        const auto is_input = [&](const std::filesystem::path& p) {
            const std::filesystem::path c = canon(p);
            return std::find(inputs.begin(), inputs.end(), c) != inputs.end();
        };
        std::vector<std::filesystem::path> stale, folders;
        for (std::filesystem::directory_iterator it(final_dir, ec), end; !ec && it != end; it.increment(ec)) {
            const std::filesystem::path rel = it->path().filename();
            if (!it->is_directory()) {
                if (is_owned_file(rel) && !wrote(rel)) stale.push_back(rel);
                continue;
            }
            const std::string n = rel.string();
            if (n != "haplotype" && n != "sample" && !is_owned_file(rel / "kmer_features.tsv.gz")) continue;
            std::error_code e2;
            for (std::filesystem::recursive_directory_iterator r(it->path(), e2), rend; !e2 && r != rend; r.increment(e2)) {
                const std::filesystem::path sub = std::filesystem::relative(r->path(), final_dir);
                if (r->is_directory()) folders.push_back(sub);
                else if (is_owned_file(sub) && !wrote(sub)) stale.push_back(sub);
            }
            folders.push_back(rel);
        }
        for (const std::filesystem::path& rel : stale) {
            std::error_code e;
            if (!is_input(final_dir / rel) && std::filesystem::remove(final_dir / rel, e)) ++removed;
        }
        // Deepest first, so a folder emptied by its subfolders' removal goes too.
        std::sort(folders.begin(), folders.end(), [](const auto& x, const auto& y) {
            return std::distance(x.begin(), x.end()) > std::distance(y.begin(), y.end());
        });
        for (const std::filesystem::path& rel : folders) {
            std::error_code e;
            if (std::filesystem::is_empty(final_dir / rel, e) && !e) std::filesystem::remove(final_dir / rel, e);
        }
    }
    std::vector<std::string> top_files, top_dirs;   // for the closing log line
    std::size_t bubble_dirs = 0;
    {
        std::vector<std::string> seen;
        for (const std::filesystem::path& rel : written) {
            const std::string first = rel.begin()->string();
            if (std::find(seen.begin(), seen.end(), first) != seen.end()) continue;
            seen.push_back(first);
            if (std::distance(rel.begin(), rel.end()) == 1) top_files.push_back(first);
            else if (first.rfind("bubble_", 0) == 0) ++bubble_dirs;
            else top_dirs.push_back(first + "/");
        }
    }
    if (graph_substrates) {
        log.info("input: haplotypes: " + std::to_string(summary.haplotypes) + "; nodes: " +
                 std::to_string(summary.nodes));
        if (options.emit_kmers)
            log.info("k-mers: " + feature_mode_label(options.feature_mode) + "; k: " +
                     std::to_string(options.kmer_size));
        log.info("sites: " + std::to_string(summary.bubbles_processed));
        if (options.emit_kmers)
            log.info("k-mer features kept: " + std::to_string(summary.features_written) + " of " +
                     std::to_string(summary.features_candidates));
        if (options.emit_graph)
            log.info("node and edge features kept: " + std::to_string(summary.node_edge_features_written) +
                     " of " + std::to_string(summary.node_edge_candidates));
    }
    if (options.emit_variant && !options.variant_vcf_path.empty())
        log.info("variant records: " + std::to_string(summary.variant_rows));
    if (summary.samples) log.info("samples: " + std::to_string(summary.samples));
    if (!moved.empty()) log.info("replaced: files of a previous run: " + std::to_string(moved.size()));
    if (removed > 0)
        log.info("removed " + std::to_string(removed) + " file(s) an earlier run left in " + options.out_dir);
    std::sort(top_files.begin(), top_files.end());
    std::sort(top_dirs.begin(), top_dirs.end());
    std::vector<std::string> listed = top_files;
    listed.insert(listed.end(), top_dirs.begin(), top_dirs.end());
    if (bubble_dirs) listed.push_back(std::to_string(bubble_dirs) + " bubble_<id>/");
    log.info("wrote " + std::to_string(summary.files_written) + " files in " + options.out_dir + "/: " +
             cli::join_with_comma(listed));
    log.done();
    return 0;
}

} // namespace panvar
