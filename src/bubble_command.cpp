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
        << "  panvar bubble -i <graph.gfa> -r <name> [-o <prefix>] [options]\n\n"
        << "Find variant sites (bubbles) in a locus graph, ordered along the reference.\n\n"
        << "Options:\n"
        << "  -i, --gfa <path>                Input GFA (required)\n"
        << "  -r, --reference-path <name>     Reference path, exact name or unique substring: orders the graph\n"
        << "                                  and the sites (required unless --snarls)\n"
        << "  -o, --out-prefix <prefix>       Output prefix (default: bubble_calls)\n"
        << "      --gtf <path>                Also write genes on reference nodes, <prefix>.bandage_genes.csv;\n"
        << "                                  needs a PanSN reference name with its interval\n"
        << "                                  (sample#hap#contig:start-end), else skipped with a warning\n"
        << "      --min-variant-bp <N>        Keep a site where some path has N bp or more between its\n"
        << "                                  boundaries (default 50; 0 = off; alias --min-interior-bp)\n"
        << "      --max-variant-bp <N>        Drop a site whose longest interior exceeds N bp (default 0 = no\n"
        << "                                  cap; alias --max-interior-bp)\n"
        << "      --min-path-support <N>      Paths that must cross a site (default 0)\n"
        << "      --min-alt-support <N>       Paths that must carry its most common non-reference allele\n"
        << "                                  (default 0)\n"
        << "      --merge-nearby-bp <N>       Merge sites at most N bp apart on the reference (default 100;\n"
        << "                                  0 = off)\n"
        << "      --no-fuse-ambiguous         Keep apart sites holding one variant the graph placed in pieces\n"
        << "      --no-flip                   Do not reorient nodes to the reference strand\n"
        << "      --snarls <path>             Take the sites' boundaries from vg snarls JSONL instead of finding\n"
        << "                                  them (vg snarls | vg view -R -j); without -r, diagnostic only\n"
        << "      --emit-snarls-jsonl <path>  Also write the boundaries found, as vg snarls JSONL\n"
        << "      --sorted-gfa-out <path>     Sorted graph (default: <prefix>.sorted.gfa)\n"
        << "      --bubbles-csv <path>        Site table (default: <prefix>.bubbles.csv)\n"
        << "      --bandage-csv <path>        Bandage colours (default: <prefix>.bandage_nodes.csv)\n"
        << "      --snarl-debug-tsv <path>    Also write every candidate considered, with why it was kept\n"
        << "  -q, --quiet                     No progress logs\n"
        << "  -h, --help                      Show this help\n";
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
    // Sites this close describe one region: a short-unit repeat split into two sites leaves each with too
    // few copies to fold, and call reports compensating length changes in place of one unit-level record.
    // Only this command's default: panphorte and refine re-detect sites with the library default.
    options.merge_nearby_bp = 100;
    // Sites whose variants could be placed in each other are one site (see fuse_ambiguous_placements): the
    // graph can split one expansion over several sites, and the size floor would keep only one piece.
    options.fuse_ambiguous = true;

    for (std::size_t i = 0; i < args.size(); ++i) {
        const std::string arg = cli::renamed_input_flag(args[i], "bubble");
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
        if (arg == "--snarls") {
            options.snarls_input_path = require_value(arg);
            continue;
        }
        if (arg == "-r" || arg == "--reference-path") {
            options.reference_path = require_value(arg);
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
        // Aliases that say what the bp filters measure: the interior span between the boundaries, not the
        // size of the difference between alleles.
        if (arg == "--min-interior-bp") {
            options.min_variant_bp = cli::parse_size_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--max-interior-bp") {
            options.max_variant_bp = cli::parse_size_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "--no-fuse-ambiguous") {
            options.fuse_ambiguous = false;
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

    // Every output path is resolved here, before anything is opened, so the alias check below sees the
    // defaults too.
    if (bubbles_csv_path.empty()) {
        bubbles_csv_path = out_prefix + ".bubbles.csv";
    }
    if (bandage_csv_path.empty()) {
        bandage_csv_path = out_prefix + ".bandage_nodes.csv";
    }
    if (!options.snarls_input_path.empty()) {
        // With --snarls the graph is taken as it is: nothing is sorted, flipped or re-emitted.
        if (!sorted_gfa_path.empty())
            std::cerr << cli::log_prefix("bubble") << "note: --sorted-gfa-out is ignored with --snarls (the graph is not sorted)\n";
        if (no_flip) std::cerr << cli::log_prefix("bubble") << "note: --no-flip is ignored with --snarls (the graph is not flipped)\n";
        if (!emit_snarls_jsonl_path.empty())
            std::cerr << cli::log_prefix("bubble") << "note: --emit-snarls-jsonl is a no-op with --snarls (the input IS the snarls)\n";
        sorted_gfa_path.clear();
        emit_snarls_jsonl_path.clear();
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

    // Refuse to write any output over any input (`bubble` writes the sorted GFA before discovery, so an
    // aliased path would destroy the graph it is about to read), and refuse two outputs naming one file. Both
    // rules live in cli_utils so this and `inspect` apply the same ones.
    cli::reject_output_collisions(
        "bubble",
        {bubbles_csv_path, bandage_csv_path, sorted_gfa_path, bandage_genes_path,
         snarl_debug_tsv_path, emit_snarls_jsonl_path},
        {gfa_path, options.snarls_input_path, gtf_path});

    // Nothing lands in its final location until the run has succeeded, so a failed run never leaves a
    // complete-looking .sorted.gfa for the next command in a pipeline.
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

    cli::RunLog log("bubble", options.quiet);
    log.info("reading: " + gfa_path);
    std::string effective_gfa = gfa_path;
    Graph graph;

    if (options.snarls_input_path.empty()) {
        // Default: internally sort+flip along the reference and find snarls (no vg/odgi).
        if (options.reference_path.empty()) {
            throw std::runtime_error(
                "Missing required input: --reference-path <name> (or pass --snarls to use vg)");
        }
        GfaModel model = read_gfa_model(gfa_path);
        // Every node a walk or a link names must be defined before the sort renumbers them: afterwards an
        // undefined id can coincide with a renumbered one, and the graph would be accepted with a haplotype
        // it never had.
        require_graph_paths(model.node_order.size(), model.paths.size(), "bubble");
        for (const GfaPath& p : model.paths)
            for (const PathStep& st : p.steps)
                if (!model.has_node(st.node_id))
                    throw std::runtime_error("bubble: path " + gfa_path_name(p) +
                                             " references a node that is not in the graph: " + st.node_id);
        for (const GfaEdge& e : model.edges)
            if (!model.has_node(e.from) || !model.has_node(e.to))
                throw std::runtime_error("bubble: link " + e.from + e.from_orient + " -> " + e.to + e.to_orient +
                                         " names a node that is not in the graph");
        // Every segment needs a sequence, on a path or not: a component of '*' segments has no bases for
        // the snarl finder to work on. The path check below only sees segments some path walks.
        for (const std::string& id : model.node_order)
            if (model.seq.at(id).empty())
                throw std::runtime_error("bubble: segment " + id + " has no sequence ('*'); every segment "
                                         "needs one");
        log.info("input: haplotypes: " + std::to_string(model.paths.size()) + "; nodes: " +
                 std::to_string(model.node_order.size()) + "; edges: " + std::to_string(model.edges.size()));
        GraphSortOptions sort_opts;
        sort_opts.reference_path = options.reference_path;
        sort_opts.flip = !no_flip;
        // The sorter resolves an exact name or a unique substring (as written, else in any case), while
        // everything downstream compares against the path name exactly: the resolved name is kept, or a
        // query differing in case would sort correctly and then match nothing downstream.
        const GraphSortResult sort_result = sort_graph_reference(model, sort_opts);
        options.reference_path = sort_result.resolved_reference;
        log.info("reference: " + options.reference_path);
        log.info("sorted along the reference: nodes turned around: " + std::to_string(sort_result.flipped_nodes) +
                 "; loops placed as one block: " + std::to_string(sort_result.cycle_blocks));
        write_gfa_model(sorted_gfa_staged, model);
        effective_gfa = sorted_gfa_path;

        graph = parse_gfa(sorted_gfa_staged, parse_options);
        // Sequences are needed because the bp filters measure interior span; overlaps must be a verified 0M
        // for the same reason, as span is summed over whole segments.
        validate_graph_paths(graph, "bubble", true, true);

        // Find boundary pairs internally with the vg-faithful cactus finder.
        options.snarl_pairs_override =
            find_top_level_snarls_cactus(snarl_input_from_model(model), &options.snarl_inward_override);
        options.snarl_source_supplied = true;
        log.info("snarls: " + std::to_string(options.snarl_pairs_override.size()) + " (built-in finder)");
    } else {
        // Override: external vg snarls on the graph as is (no internal sort).
        graph = parse_gfa(gfa_path, parse_options);
        // The same graph contract as internal mode, so a malformed graph (a path step with no link behind
        // it) is refused through either door.
        validate_graph_paths(graph, "bubble", true, true);
        // The same name resolution as internal mode, so the reference is matched by its exact name.
        options.reference_path =
            resolve_reference_path_name(graph, options.reference_path, "bubble");
        // Imported boundaries are an unordered pair and this mode does not sort, so without a reference
        // there is nothing to orient them by, and every consumer reading them as an interval would read a
        // coin flip.
        if (options.reference_path.empty()) {
            std::cerr << cli::log_prefix("bubble") << "warning: --snarls without --reference-path produces DIAGNOSTIC-ONLY "
                         "output, not input for the rest of the pipeline. Imported boundaries are an "
                         "unordered pair and this mode does not sort, so source/sink do not mean "
                         "reference-left/right, boundary orientations are defaults rather than "
                         "measurements, reference-coordinate merging is skipped, and no sorted graph "
                         "is written. `call` anchors coordinates on the source, so a downstream run "
                         "would place records by a coin flip. Supply --reference-path for output any "
                         "later module should consume.\n";
        }
        log.info("input: haplotypes: " + std::to_string(graph.paths.size()) + "; nodes: " +
                 std::to_string(graph.nodes.size()));
        log.info("reference: " + (options.reference_path.empty() ? std::string("none") : options.reference_path));
        log.info("snarls: read from " + options.snarls_input_path);
    }

    const auto report = call_bubbles_report(graph, options);
    const auto& bubbles = report.bubbles;

    write_bubbles_csv(bubbles_csv_staged, bubbles);
    write_bandage_node_colors_csv(bandage_csv_staged, bubbles, report.non_snp_bubbles);
    if (!bandage_genes_path.empty()) {
        if (!emit_gene_annotation(graph, options.reference_path, gtf_path, bandage_genes_staged, "", "bubble")) {
            bandage_genes_path.clear();
        }
    }
    if (!snarl_debug_staged.empty()) {
        write_snarl_debug_tsv(snarl_debug_staged, report.snarl_debug);
    }
    if (!emit_snarls_jsonl_path.empty()) {
        // The internally found (cactus) snarls, as vg writes them: each boundary is a handle, and `backward`
        // says which side the interior hangs off (start: false = its end side; end: false = its start side),
        // since the pair alone is unordered and re-imported with --snarls the sides would be voted from
        // the walks. Written even when empty, so an earlier run's file at this path is replaced.
        std::ofstream js(emit_snarls_staged);
        if (!js) throw std::runtime_error("Failed to write snarls JSONL: " + emit_snarls_jsonl_path);
        const auto& pairs = options.snarl_pairs_override;
        const auto& inward = options.snarl_inward_override;
        for (std::size_t i = 0; i < pairs.size(); ++i) {
            const bool start_back = inward.size() == pairs.size() ? !inward[i].source_end : false;
            const bool end_back = inward.size() == pairs.size() ? inward[i].sink_end : false;
            js << "{\"start\": {\"node_id\": \"" << pairs[i].first << "\", \"backward\": "
               << (start_back ? "true" : "false") << "}, \"end\": {\"node_id\": \"" << pairs[i].second
               << "\", \"backward\": " << (end_back ? "true" : "false") << "}}\n";
        }
        cli::close_output(js, emit_snarls_jsonl_path);
    }

    // One line per step, every run, so two runs' logs compare line by line.
    const BubbleRunCounts& c = report.counts;
    log.info("candidate sites: " + std::to_string(c.candidates) + "; widened around a loop: " + std::to_string(c.widened));
    log.info("passing the filters: " + std::to_string(c.passed));
    log.info("merged nearby: " + std::to_string(c.merged) + "; fused over placements: " + std::to_string(c.fused) +
             "; dropped as overlapping: " + std::to_string(c.overlaps_dropped) + "; rescued from pieces under the floor: " +
             std::to_string(c.rescued));
    log.info("inversion sites: " + std::to_string(c.inversion_sites));
    log.info("sites: " + std::to_string(bubbles.size()));

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
    if (!emit_snarls_jsonl_path.empty()) {
        outputs.push_back(emit_snarls_jsonl_path);
    }
    // Default-named files an earlier run left at this prefix describe that run, not this one, when
    // this run wrote none: the sorted graph under --snarls, the gene table without one (no --gtf,
    // or a reference that names no region). Never an input of this run.
    // Only when the graph was read from a regular file: the identity test below cannot tell whether a
    // pipe, a fifo or a /dev/fd path (which reports the redirected file as regular) was that file.
    const auto from_device = [](const std::string& p) {
        std::error_code ec;
        const std::string canon = std::filesystem::weakly_canonical(p, ec).string();
        return p.rfind("/dev/", 0) == 0 || canon.rfind("/dev/", 0) == 0;
    };
    std::error_code reg_ec;
    if (std::filesystem::is_regular_file(gfa_path, reg_ec) && !from_device(gfa_path)) {
        std::vector<std::string> stale;
        if (!options.snarls_input_path.empty()) stale.push_back(out_prefix + ".sorted.gfa");
        if (bandage_genes_path.empty()) stale.push_back(out_prefix + ".bandage_genes.csv");
        for (const std::string& path : stale) {
            std::error_code ec;
            // Nor a file this run just wrote under that name (--bubbles-csv <prefix>.bandage_genes.csv, say).
            bool keep = false;
            for (const std::string& in : {gfa_path, options.snarls_input_path, gtf_path})
                keep = keep || (!in.empty() && std::filesystem::equivalent(path, in, ec));
            for (const std::string& out : outputs)
                keep = keep || std::filesystem::equivalent(path, out, ec);
            if (!keep && std::filesystem::remove(path, ec))
                log.info("removed " + path + ", left by an earlier run at this prefix");
        }
    }
    log.wrote(outputs);
    log.done();

    return 0;
}

} // namespace panvar
