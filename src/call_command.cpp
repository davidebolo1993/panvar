#include "panvar/call_command.hpp"

#include "panvar/cli_utils.hpp"
#include "panvar/gfa.hpp"
#include "panvar/graph_utils.hpp"
#include "panvar/variant_call.hpp"

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace panvar {
namespace {

void print_call_help() {
    std::cout
        << "Usage:\n"
        << "  panvar call -i <graph.gfa> (-b <prefix> | -c <bubbles.csv>) -r <name> -o <prefix> [options]\n"
        << "\n"
        << "Call structural variants: each haplotype's walk through a site is compared with the\n"
        << "reference's, pieces of one event are joined, and the same event in different haplotypes\n"
        << "becomes one VCF record. Folded repeats are called as copy numbers.\n"
        << "\n"
        << "Options:\n"
        << "  -i, --gfa <path>                  Sorted GFA from bubble, panphorte or refine (required)\n"
        << "  -b, --bubble-prefix <prefix>      Sites from <prefix>.bubbles.csv (required unless -c)\n"
        << "  -c, --bubbles-csv <path>          Site table, by path\n"
        << "      --panphorte-prefix <prefix>   Unit-graph table <prefix>.panphorte.unit_graph.tsv, to count\n"
        << "                                    folded repeats' copies (default: the site table's prefix)\n"
        << "  -r, --reference-path <name>       Reference path, exact name or unique substring (required)\n"
        << "  -o, --out-prefix <prefix>         Output prefix (required)\n"
        << "      --gtf <path>                  Genes each variant touches, and copy number per gene (needs a\n"
        << "                                    PanSN reference name)\n"
        << "      --cn                          Also count the copies of a duplicated stretch the graph did not\n"
        << "                                    fold (CN_METHOD=MODULE)\n"
        << "      --bubble-id <N>               Only this site; repeatable (default: all)\n"
        << "      --min-sv-bp <N>               Smallest event reported (default 50)\n"
        << "      --rescue-min-bp <N>           Smallest piece kept for joining into a larger event (default:\n"
        << "                                    half --min-sv-bp)\n"
        << "      --no-align-alleles            Compare with the reference node by node, not by base alignment\n"
        << "      --merge-distance-bp <N>       Merge events of different haplotypes within N bp plus the\n"
        << "                                    smaller event's size (default 100)\n"
        << "      --merge-jaccard <X>           Shared nodes needed to merge events (default 0.80)\n"
        << "      --merge-seq-identity <X>      Sequence identity needed to merge events (default 0.80)\n"
        << "      --merge-size-ratio <X>        Smaller/larger size needed to merge events (default 0 = as\n"
        << "                                    --merge-seq-identity)\n"
        << "      --min-haplotypes <N>          Drop records carried by fewer than N haplotypes (default 1)\n"
        << "      --min-alt-af <X>              Drop records whose AF (AC/AN, the haplotypes or with --diploid\n"
        << "                                    the haplotype slots that have an allele there) is below X\n"
        << "                                    (default 0)\n"
        << "      --classify-ins                Tell new insertions from duplications with minimap2\n"
        << "      --minimap-preset <name>       asm5, asm10 or asm20, for --classify-ins (default asm20)\n"
        << "      --minimap-best-n <N>          Chains minimap2 keeps, for --classify-ins (default 8)\n"
        << "      --ins-dup-min-identity <X>    Identity for an insertion to count as a duplication\n"
        << "                                    (default 0.90)\n"
        << "      --diploid                     One VCF column per sample (paths SAMPLE#HAP#...), not per path\n"
        << "      --no-per-bubble-vcf           Only the region VCF, not one per site\n"
        << "      --max-seq-bp <N>              Leave INSSEQ, DELSEQ, INVSEQ and each CN_FORMS sequence out above\n"
        << "                                    N bases, so VCF tools read the lines (default 20000; 0 = no limit)\n"
        << "      --no-variant-nodes            Do not write <prefix>.variant_nodes.tsv (the nodes describe\n"
        << "                                    reads) nor <prefix>.carrier_events.tsv\n"
        << "  -t, --threads <N>                 Worker threads (default 0 = all cores)\n"
        << "  -q, --quiet                       No progress logs\n"
        << "  -h, --help                        Show this help\n";
}

} // namespace

int run_call_command(const std::vector<std::string>& args) {
    if (args.empty()) {
        print_call_help();
        return 0;
    }

    std::string gfa_path;
    std::string bubble_prefix_in;
    std::string panphorte_prefix;
    bool panphorte_prefix_given = false;
    std::string out_prefix;
    VariantCallOptions options;

    for (std::size_t i = 0; i < args.size(); ++i) {
        const std::string arg = cli::renamed_input_flag(args[i], "call");
        auto require_value = [&](const std::string& flag) -> const std::string& {
            if (i + 1 >= args.size()) {
                throw std::runtime_error("Missing value after " + flag);
            }
            return args[++i];
        };

        if (arg == "-h" || arg == "--help") {
            print_call_help();
            return 0;
        }
        if (arg == "-i" || arg == "--gfa") {
            gfa_path = require_value(arg);
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
        if (arg == "--panphorte-prefix") {
            panphorte_prefix = require_value(arg);
            panphorte_prefix_given = true;
            continue;
        }
        if (arg == "-r" || arg == "--reference-path") {
            options.reference_path = require_value(arg);
            continue;
        }
        if (arg == "-o" || arg == "--out-prefix") {
            out_prefix = require_value(arg);
            continue;
        }
        if (arg == "--min-sv-bp") {
            options.min_sv_bp = cli::parse_size_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--max-seq-bp") {
            options.max_seq_bp = cli::parse_size_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--merge-distance-bp") {
            options.merge_distance_bp = cli::parse_size_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--merge-jaccard") {
            options.merge_jaccard = cli::parse_unit_fraction_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--merge-seq-identity") {
            options.merge_seq_identity = cli::parse_unit_fraction_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--merge-size-ratio") {
            options.merge_size_ratio = cli::parse_unit_fraction_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--min-haplotypes") {
            options.min_haplotypes = cli::parse_size_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--min-alt-af") {
            options.min_maf = cli::parse_unit_fraction_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--no-align-alleles") {
            options.align_alleles = false;
            continue;
        }
        if (arg == "--diploid") {
            options.diploid = true;
            continue;
        }
        if (arg == "--rescue-min-bp") {
            options.rescue_min_bp = cli::parse_size_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--classify-ins") {
            options.classify_ins = true;
            continue;
        }
        if (arg == "--cn") {                 // copy-number calling: the module route
            options.cn = true;
            continue;
        }
        if (arg == "--gtf") {
            options.gtf_path = require_value(arg);
            continue;
        }
        if (arg == "--minimap-preset") {
            // Only the three the help advertises; anything else reaches minimap2 as an unknown preset.
            options.minimap_preset = require_value(arg);
            if (options.minimap_preset != "asm5" && options.minimap_preset != "asm10" &&
                options.minimap_preset != "asm20") {
                throw std::runtime_error("--minimap-preset must be asm5, asm10 or asm20 (got '" +
                                         options.minimap_preset + "')");
            }
            continue;
        }
        if (arg == "--minimap-best-n") {
            options.minimap_best_n = cli::parse_size_arg(arg, require_value(arg));
            if (options.minimap_best_n == 0) {
                throw std::runtime_error("--minimap-best-n must be >= 1");
            }
            continue;
        }
        if (arg == "--ins-dup-min-identity") {
            options.ins_dup_min_identity = cli::parse_unit_fraction_arg(arg, require_value(arg));
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
        if (arg == "--no-per-bubble-vcf") {
            options.write_per_bubble_vcf = false;
            continue;
        }
        if (arg == "--no-variant-nodes") {
            options.write_variant_nodes = false;
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
        throw std::runtime_error("Unknown option for call: " + arg);
    }

    if (gfa_path.empty()) {
        throw std::runtime_error("Missing required input: --gfa <path>");
    }
    options.gfa_path = gfa_path;
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
        throw std::runtime_error("call requires --bubble-prefix <prefix> or --bubbles-csv <path>");
    }
    if (options.reference_path.empty()) {
        throw std::runtime_error("--reference-path is required for module 'call'");
    }
    if (out_prefix.empty()) {
        throw std::runtime_error("call requires -o/--out-prefix <prefix>");
    }
    options.out_prefix = out_prefix;
    // panphorte's unit graphs travel with the site table (refine forwards them under its own prefix). Named
    // explicitly, the table must exist; found by default, it is used when it is there.
    {
        const std::string tail = ".bubbles.csv";
        const std::string& csv = options.bubbles_csv_in;
        if (!panphorte_prefix_given)
            panphorte_prefix = !bubble_prefix_in.empty() ? bubble_prefix_in
                             : csv.size() > tail.size() && csv.compare(csv.size() - tail.size(), tail.size(), tail) == 0
                                 ? csv.substr(0, csv.size() - tail.size()) : std::string();
        const std::string ug = panphorte_prefix.empty() ? std::string() : panphorte_prefix + ".panphorte.unit_graph.tsv";
        std::error_code ec;
        if (!ug.empty() && std::filesystem::is_regular_file(ug, ec)) options.unit_graph_tsv = ug;
        else if (panphorte_prefix_given)
            throw std::runtime_error("call: --panphorte-prefix " + panphorte_prefix + ": no " + ug);
    }

    cli::RunLog log("call", options.quiet);
    log.info("reading: " + gfa_path);
    ParseGfaOptions parse_options;
    parse_options.include_paths = true;
    parse_options.include_sequences = true;
    Graph graph = parse_gfa(gfa_path, parse_options);
    require_graph_paths(graph.nodes.size(), graph.paths.size(), "call");
    // Uppercase, and any letter other than ACGT as N, so a soft-masked or IUPAC graph still writes valid
    // REF/ALT, does not split one allele into case variants, and merges sequences regardless of case.
    for (auto& [id, node] : graph.nodes) {
        (void)id;
        for (char& c : node.sequence) {
            const char u = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
            c = (u == 'A' || u == 'C' || u == 'G' || u == 'T') ? u : 'N';
        }
    }

    log.info("input: haplotypes: " + std::to_string(graph.paths.size()) + "; nodes: " +
             std::to_string(graph.nodes.size()));
    // The path the query resolves to (exact, else a unique substring), as the VCF's ##reference names it.
    log.info("reference: " + resolve_reference_path_name(graph, options.reference_path, "call"));
    log.info(options.unit_graph_tsv.empty() ? std::string("folded repeats: no unit-graph table")
                                            : "folded repeats: " + options.unit_graph_tsv);

    VariantCallSummary summary;
    call_variants(graph, options, &summary);

    log.info("sites: " + std::to_string(summary.bubbles_seen) + "; with records: " +
             std::to_string(summary.bubbles_with_calls));
    log.info("records: " + std::to_string(summary.records_written) + "; DEL: " + std::to_string(summary.del) +
             "; INS: " + std::to_string(summary.ins) + "; INV: " + std::to_string(summary.inv) + "; CNV: " +
             std::to_string(summary.dup) + "; SUB: " + std::to_string(summary.sub));
    // A module copy number not called because the graph is tangled there is a loss, and the remedy is
    // upstream: said even with -q.
    auto id_list = [](std::vector<std::size_t> ids) {
        std::sort(ids.begin(), ids.end());
        ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
        std::string out;
        for (std::size_t k = 0; k < ids.size(); ++k) out += (k ? "," : "") + std::to_string(ids[k]);
        return out;
    };
    if (!summary.tangle_ids.empty())
        std::cerr << cli::log_prefix("call") << "warning: no module copy number at site(s) " << id_list(summary.tangle_ids)
                  << ": ten or more nodes there have fifty neighbours or more on one end, a tangle whose node counts are "
                     "no copy number; `panvar rebuild` may untangle the graph\n";
    if (!summary.oversized_ids.empty())
        std::cerr << cli::log_prefix("call") << "warning: no module copy number at site(s) " << id_list(summary.oversized_ids)
                  << ": the module would span more than 80% of the reference, scattered revisits rather than "
                     "copies; `panvar rebuild` may untangle the graph\n";
    // Genotypes set missing (GT .) for a reason the run cannot avoid, each logged with its count so the
    // loss has a stated reason.
    if (summary.skipped_large_segments > 0) {
        std::string ids;
        for (std::size_t k = 0; k < summary.unaligned_bubbles.size(); ++k)
            ids += (k ? "," : "") + std::to_string(summary.unaligned_bubbles[k]);
        log.info("blocks too large to align: " + std::to_string(summary.skipped_large_segments) +
                 (ids.empty() ? std::string() : " (sites " + ids + ")") + "; GT . there: " +
                 std::to_string(summary.unaligned_missing));
    }
    if (summary.covered_missing > 0)
        log.info("GT . within a larger deletion of the haplotype's own: " + std::to_string(summary.covered_missing));
    if (summary.copy_carriers + summary.copy_own_event > 0)
        log.info("module copy changes: carriers " + std::to_string(summary.copy_carriers) + ", a larger event of their own " +
                 std::to_string(summary.copy_own_event) + "; carriers with no change holding a whole copy " +
                 std::to_string(summary.copy_no_form) + ", whose copies' anchor passes differ from CN " +
                 std::to_string(summary.copy_pass_off));
    if (summary.ref_gap_bubbles > 0) {
        std::string ids;
        for (std::size_t k = 0; k < summary.ref_gap_ids.size(); ++k)
            ids += (k ? "," : "") + std::to_string(summary.ref_gap_ids[k]);
        log.info("sites skipped, assembly gap in the reference: " + std::to_string(summary.ref_gap_bubbles) +
                 " (" + ids + ")");
    }
    if (summary.gap_missing > 0) {
        log.info("GT . across an assembly gap: " + std::to_string(summary.gap_missing));
    }
    // Detection only: walks and copy numbers on these bubbles are measured as anywhere else, and every
    // record there carries INFO/MULTI_CROSSING.
    if (summary.multi_crossing_bubbles > 0) {
        log.info("paths crossing a site more than once: " + std::to_string(summary.multi_crossing_paths) +
                 "; sites: " + std::to_string(summary.multi_crossing_bubbles) + " (records flagged MULTI_CROSSING)");
    }
    if (summary.diploid_records_dropped > 0)
        log.info("--diploid: records with no sample carrier, not written: " +
                 std::to_string(summary.diploid_records_dropped));
    if (summary.multi_crossing_unscored_paths > 0) {
        std::string ids;
        for (std::size_t k = 0; k < summary.multi_crossing_unscored_ids.size() && k < 10; ++k)
            ids += (k ? "," : "") + std::to_string(summary.multi_crossing_unscored_ids[k]);
        if (summary.multi_crossing_unscored_ids.size() > 10) ids += ",...";
        std::cerr << cli::log_prefix("call") << "warning: " << summary.multi_crossing_unscored_paths << " path(s) at "
                  << summary.multi_crossing_unscored_bubbles << " bubble(s) (" << ids
                  << ") cross a second time with an allele --min-sv-bp or more longer or shorter than every "
                     "reference crossing, outside the walk that is scored: an SV-sized change no record "
                     "describes (MULTI_CROSSING_UNSCORED)\n";
    }

    // Every output shares the prefix's directory. The per-site VCFs, one per site with records, are named
    // once as a pattern with their count, so the line stays short and the file count is the real one.
    const std::filesystem::path out_dir = std::filesystem::path(options.out_prefix).parent_path();
    const std::string base = std::filesystem::path(options.out_prefix).filename().string();
    std::vector<std::string> outputs;
    outputs.push_back(base + ".region.vcf");
    std::size_t n_files = 1;
    if (summary.per_bubble_vcfs > 0) {
        outputs.push_back(base + ".bubble_<id>.vcf (" + std::to_string(summary.per_bubble_vcfs) +
                          (summary.per_bubble_vcfs == 1 ? " file)" : " files)"));
        n_files += summary.per_bubble_vcfs;
    }
    if (options.write_variant_nodes) {
        outputs.push_back(base + ".variant_nodes.tsv");
        outputs.push_back(base + ".carrier_events.tsv");
        n_files += 2;
    }
    if (summary.gene_tables_written) {
        outputs.push_back(base + ".node_genes.tsv");
        outputs.push_back(base + ".dup_gene_cn.tsv");
        n_files += 2;
    }
    std::string listed;
    for (std::size_t k = 0; k < outputs.size(); ++k) listed += (k ? ", " : "") + outputs[k];
    log.info("wrote " + std::to_string(n_files) + (n_files == 1 ? " file" : " files") +
             (out_dir.empty() ? std::string() : " in " + out_dir.string() + "/") + ": " + listed);
    log.done();
    return 0;
}

} // namespace panvar
