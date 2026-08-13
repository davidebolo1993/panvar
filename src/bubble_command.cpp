#include "panvar/bubble_command.hpp"

#include "panvar/bubbles.hpp"
#include <filesystem>
#include "panvar/graph_utils.hpp"
#include "panvar/cli_utils.hpp"
#include "panvar/gfa.hpp"
#include "panvar/gfa_io.hpp"
#include "panvar/graph_sort.hpp"
#include "panvar/gtf.hpp"
#include "panvar/integrated_snarls.hpp"
#include "panvar/output.hpp"

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace panvar {
namespace {

void print_bubble_help() {
    std::cout
        << "Usage:\n"
        << "  panvar bubble -i <graph.gfa> -r <name> [-o <prefix>] [--superbubbles] [options]\n\n"
        << "By default the graph is sorted+flipped along the reference internally and snarls are\n"
        << "found with an internal cactus decomposition (the same 3-edge-connected-component method\n"
        << "as 'vg snarls'); with --superbubbles only the acyclic snarls (ultrabubbles, i.e. the\n"
        << "acyclic subset of the cactus snarls) are kept. Pass --snarls-in to use an external vg\n"
        << "snarls JSONL instead.\n\n"
        << "Options:\n"
        << "  -i, --gfa <path>                 Input GFA file (required)\n"
        << "  -r, --reference-path <name>      Reference path name or unique substring; orders the\n"
        << "                                    internal sort/flip + snarl finder (required unless\n"
        << "                                    --snarls-in is given)\n"
        << "  -s, --superbubbles               Emit only acyclic superbubbles (default: all snarls)\n"
        << "      --no-flip                    Do not reorient nodes to the reference forward strand\n"
        << "      --sorted-gfa-out <path>      Internally-sorted GFA output (default: <prefix>.sorted.gfa)\n"
        << "      --emit-snarls-jsonl <path>   Also write the internal snarls as a vg-style JSONL\n"
        << "      --snarls-in <path>           Override: snarl JSONL from 'vg snarls -A integrated | vg view -R -j' (skips\n"
        << "                                    internal sort + finding)\n"
        << "  -o, --out-prefix <prefix>        Output prefix (default: bubble_calls)\n"
        << "      --bubbles-csv <path>         Explicit bubbles CSV output path\n"
        << "      --bandage-csv <path>         Explicit Bandage color CSV output path\n"
        << "      --gtf <path>                 Reference-coordinate GTF; project genes onto reference\n"
        << "                                    nodes and write <prefix>.bandage_genes.csv (Bandage).\n"
        << "                                    Requires a PanSN reference path; skipped otherwise.\n"
        << "      --snarl-debug-tsv <path>     Optional diagnostics TSV for snarl candidates\n"
        << "      --min-variant-bp <N>         Keep bubbles with at least one path carrying >= N bp\n"
        << "      --max-variant-bp <N>         Largest variant to keep: drop bubbles whose longest path\n"
        << "                                   exceeds N bp (0 = no cap; tames hypervariable tangles)\n"
        << "                                    inside the bubble (default: 50, 0=disable)\n"
        << "      --min-path-support <N>       Require at least N supporting P/W paths (default: 0).\n"
        << "                                   This is TRAVERSAL support -- on a fully-typed panel nearly\n"
        << "                                   every haplotype crosses nearly every bubble, so it says\n"
        << "                                   little about any particular allele\n"
        << "      --min-alt-support <N>        Require the best-supported NON-REFERENCE allele to have at\n"
        << "                                   least N supporting paths. This is what a support filter is\n"
        << "                                   usually wanted for (default: 0)\n"
        << "      --merge-nearby-bp <N>        Merge nearby bubbles, after the base filters and again\n"
        << "                                   after merging, when the reference sequence STRICTLY\n"
        << "                                   BETWEEN the facing boundaries is <= N bp. Ordering and\n"
        << "                                   distance are both in reference coordinates, so bubbles\n"
        << "                                   that abut have a gap of 0 (default: 0, disabled)\n"
        << "  -q, --quiet                      Disable the progress bar\n"
        << "  -h, --help                       Show this help\n";
}

} // namespace

int run_bubble_command(const std::vector<std::string>& args) {
    if (args.empty()) {
        print_bubble_help();
        return 0;
    }

    std::string gfa_path;
    std::string out_prefix = "bubble_calls";
    std::string bubbles_csv_path;
    std::string bandage_csv_path;
    std::string snarl_debug_tsv_path;
    std::string sorted_gfa_path;
    std::string emit_snarls_jsonl_path;
    std::string gtf_path;
    bool no_flip = false;

    BubbleCallOptions options;

    for (std::size_t i = 0; i < args.size(); ++i) {
        const std::string& arg = args[i];
        auto require_value = [&](const std::string& flag) -> const std::string& {
            if (i + 1 >= args.size()) {
                throw std::runtime_error("Missing value after " + flag);
            }
            return args[++i];
        };

        if (arg == "-h" || arg == "--help") {
            print_bubble_help();
            return 0;
        }
        if (arg == "-i" || arg == "--gfa") {
            gfa_path = require_value(arg);
            continue;
        }
        if (arg == "-o" || arg == "--out-prefix") {
            out_prefix = require_value(arg);
            continue;
        }
        if (arg == "--bubbles-csv") {
            bubbles_csv_path = require_value(arg);
            continue;
        }
        if (arg == "--bandage-csv") {
            bandage_csv_path = require_value(arg);
            continue;
        }
        if (arg == "--gtf") {
            gtf_path = require_value(arg);
            continue;
        }
        if (arg == "--snarls-in") {
            options.snarls_input_path = require_value(arg);
            continue;
        }
        if (arg == "-r" || arg == "--reference-path") {
            options.reference_path = require_value(arg);
            continue;
        }
        if (arg == "-s" || arg == "--superbubbles") {
            options.superbubbles_only = true;
            continue;
        }
        if (arg == "--no-flip") {
            no_flip = true;
            continue;
        }
        if (arg == "--sorted-gfa-out") {
            sorted_gfa_path = require_value(arg);
            continue;
        }
        if (arg == "--emit-snarls-jsonl") {
            emit_snarls_jsonl_path = require_value(arg);
            continue;
        }
        if (arg == "--snarl-debug-tsv") {
            snarl_debug_tsv_path = require_value(arg);
            options.collect_snarl_debug = true;
            continue;
        }
        if (arg == "--min-variant-bp") {
            options.min_variant_bp = cli::parse_size_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--max-variant-bp") {
            options.max_variant_bp = cli::parse_size_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--min-path-support") {
            options.min_path_support = cli::parse_size_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--min-alt-support") {
            options.min_alt_support = cli::parse_size_arg(arg, require_value(arg));
            continue;
        }
        // Aliases that say what the bp filters actually measure: the interior SPAN between the
        // boundaries, not the size of the difference between alleles.
        if (arg == "--min-interior-bp") {
            options.min_variant_bp = cli::parse_size_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--max-interior-bp") {
            options.max_variant_bp = cli::parse_size_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--merge-nearby-bp") {
            options.merge_nearby_bp = cli::parse_size_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "-q" || arg == "--quiet") {
            options.quiet = true;
            continue;
        }
        throw std::runtime_error("Unknown option: " + arg);
    }

    if (gfa_path.empty()) {
        throw std::runtime_error("Missing required input: --gfa <path>");
    }

    // Every output path is resolved here, BEFORE anything is opened, so the alias check below sees the
    // defaults too -- the sorted GFA used to be named after the check had already run, which left the
    // one output written before discovery unchecked.
    if (bubbles_csv_path.empty()) {
        bubbles_csv_path = out_prefix + ".bubbles.csv";
    }
    if (bandage_csv_path.empty()) {
        bandage_csv_path = out_prefix + ".bandage_nodes.csv";
    }
    if (sorted_gfa_path.empty() && options.snarls_input_path.empty()) {
        sorted_gfa_path = out_prefix + ".sorted.gfa";
    }
    std::string bandage_genes_path;
    if (!gtf_path.empty()) {
        bandage_genes_path = out_prefix + ".bandage_genes.csv";
    }

    ParseGfaOptions parse_options;
    parse_options.include_paths = true;
    parse_options.include_sequences = true;

    // Refuse to write any output over any input. `bubble` writes the sorted GFA before discovery even
    // begins, so an aliased path destroys the graph it is about to read.
    {
        const std::string* inputs[] = {&gfa_path, &options.snarls_input_path, &gtf_path};
        const std::string* outs[] = {&bubbles_csv_path, &bandage_csv_path, &sorted_gfa_path,
                                     &bandage_genes_path, &snarl_debug_tsv_path,
                                     &emit_snarls_jsonl_path};
        for (const std::string* in : inputs) {
            if (in->empty()) continue;
            std::error_code ec;
            const std::filesystem::path in_p = std::filesystem::weakly_canonical(*in, ec);
            if (ec || in_p.empty()) continue;
            for (const std::string* out : outs) {
                if (out->empty()) continue;
                std::error_code e2;
                const std::filesystem::path o = std::filesystem::weakly_canonical(*out, e2);
                if (!e2 && in_p == o)
                    throw std::runtime_error("bubble: output '" + *out + "' is the same file as input '" +
                                             *in + "'; refusing to overwrite it");
            }
        }
    }

    // Nothing lands in its final location until the run has succeeded: a malformed input used to exit
    // non-zero having already left a complete-looking .sorted.gfa behind, which the next command in a
    // pipeline would happily consume.
    cli::StagedOutputs staged("bubble");
    const std::string bubbles_csv_staged = staged.stage(bubbles_csv_path);
    const std::string bandage_csv_staged = staged.stage(bandage_csv_path);
    const std::string sorted_gfa_staged =
        sorted_gfa_path.empty() ? std::string() : staged.stage(sorted_gfa_path);
    const std::string bandage_genes_staged =
        bandage_genes_path.empty() ? std::string() : staged.stage(bandage_genes_path);
    const std::string snarl_debug_staged =
        snarl_debug_tsv_path.empty() ? std::string() : staged.stage(snarl_debug_tsv_path);
    const std::string emit_snarls_staged =
        emit_snarls_jsonl_path.empty() ? std::string() : staged.stage(emit_snarls_jsonl_path);

    std::string site_mode;
    std::string effective_gfa = gfa_path;
    Graph graph;

    if (options.snarls_input_path.empty()) {
        // Default: internally sort+flip along the reference and find snarls (no vg/odgi).
        if (options.reference_path.empty()) {
            throw std::runtime_error(
                "Missing required input: --reference-path <name> (or pass --snarls-in to use vg)");
        }
        GfaModel model = read_gfa_model(gfa_path);
        GraphSortOptions sort_opts;
        sort_opts.reference_path = options.reference_path;
        sort_opts.flip = !no_flip;
        // The sorter resolves an exact name, a unique case-insensitive match or a unique substring;
        // everything downstream compares against the path name EXACTLY. Keeping the user's query here
        // meant `-r FULL` for a path named `full` sorted correctly and then matched nothing: bubbles
        // came back unoriented (source=3, sink=1 on a three-node graph), reference allele support read
        // 0, and merging had no reference coordinates to order by -- all silently.
        const GraphSortResult sort_result = sort_graph_reference(model, sort_opts);
        options.reference_path = sort_result.resolved_reference;
        write_gfa_model(sorted_gfa_staged, model);
        effective_gfa = sorted_gfa_path;

        graph = parse_gfa(sorted_gfa_staged, parse_options);
        // Sequences are needed because the bp filters measure interior span; overlaps must be a
        // verified 0M for the same reason -- span is summed over whole segments.
        validate_graph_paths(graph, "bubble", true, true);
        if (graph.paths.empty()) {
            throw std::runtime_error("Input GFA has no P/W paths; snarl finding requires path walks");
        }

        // Find boundary pairs internally with the vg-faithful cactus finder. --superbubbles
        // then keeps only the acyclic snarls (= superbubbles), filtered in call_bubbles_report.
        options.snarl_pairs_override = find_top_level_snarls_cactus(snarl_input_from_model(model));
    options.snarl_source_supplied = true;
        site_mode = options.superbubbles_only ? "superbubble (internal, cactus + acyclic)"
                                              : "snarl (internal, cactus)";
    } else {
        // Legacy override: external vg snarls on the graph as-is (no internal sort).
        graph = parse_gfa(gfa_path, parse_options);
        // The same graph contract as internal mode. Skipping it here meant a malformed graph -- a path
        // step with no link behind it -- was accepted through this door with exit 0 and an empty table,
        // while the identical file was refused through the other one.
        validate_graph_paths(graph, "bubble", true, true);
        if (graph.paths.empty()) {
            throw std::runtime_error("Input GFA has no P/W paths; snarl refinement requires path walks");
        }
        // Imported snarl boundaries are an unordered pair and this mode does not sort, so without a
        // reference there is nothing to orient them by: source/sink do not mean reference-left/right,
        // and every consumer that reads them as an interval is then reading a coin flip.
        if (options.reference_path.empty()) {
            std::cerr << "warning: --snarls-in without --reference-path: bubble boundaries are "
                         "unordered, so source/sink do not mean reference-left/right\n";
        }
        site_mode = "snarl (JSONL import)";
    }

    const auto report = call_bubbles_report(graph, options);
    const auto& bubbles = report.bubbles;

    write_bubbles_csv(bubbles_csv_staged, bubbles);
    write_bandage_node_colors_csv(bandage_csv_staged, bubbles, report.non_snp_bubbles);
    if (!bandage_genes_path.empty()) {
        if (!emit_gene_annotation(graph, options.reference_path, gtf_path, bandage_genes_staged)) {
            bandage_genes_path.clear();
        }
    }
    if (!snarl_debug_staged.empty()) {
        write_snarl_debug_tsv(snarl_debug_staged, report.snarl_debug);
    }
    if (!emit_snarls_jsonl_path.empty()) {
        // Emit the internally found (cactus) snarl pairs. In --snarls-in (legacy) mode the
        // snarls already exist as the input file, so there is nothing new to emit.
        if (options.snarl_pairs_override.empty()) {
            std::cerr << "note: --emit-snarls-jsonl is a no-op with --snarls-in (the input IS the snarls)\n";
        } else {
            std::ofstream js(emit_snarls_staged);
            for (const auto& [s, t] : options.snarl_pairs_override) {
                js << "{\"start\": {\"node_id\": \"" << s << "\"}, \"end\": {\"node_id\": \"" << t
                   << "\"}}\n";
            }
        }
    }

    std::size_t inversion_signal_count = 0;
    std::size_t long_path_positive_count = 0;
    for (const auto& bubble : bubbles) {
        if (bubble.inversion_signal) {
            ++inversion_signal_count;
        }
        if (bubble.long_path_support > 0) {
            ++long_path_positive_count;
        }
    }

    cli::RunLog log("bubble", options.quiet);
    log.info("input " + gfa_path + " (" + std::to_string(graph.nodes.size()) + " nodes, " +
             std::to_string(graph.paths.size()) + " paths; " + site_mode + ")");
    log.info("found " + std::to_string(bubbles.size()) + " bubbles (" +
             std::to_string(long_path_positive_count) + " ≥ min-bp, " +
             std::to_string(inversion_signal_count) + " inversion-flagged)");

    // Everything succeeded: move the staged files into place. Reported only after this, so the log
    // never names a file that is not there.
    staged.commit();

    std::vector<std::string> outputs;
    if (options.snarls_input_path.empty()) {
        outputs.push_back(effective_gfa);
    }
    outputs.push_back(bubbles_csv_path);
    outputs.push_back(bandage_csv_path);
    if (!bandage_genes_path.empty()) {
        outputs.push_back(bandage_genes_path);
    }
    if (!snarl_debug_tsv_path.empty()) {
        outputs.push_back(snarl_debug_tsv_path);
    }
    log.wrote(outputs);
    log.done();

    return 0;
}

} // namespace panvar
