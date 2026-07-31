#include "panvar/genotype_command.hpp"

#include "panvar/align.hpp"
#include "panvar/cli_utils.hpp"
#include "panvar/genotype_blocks.hpp"
#include "panvar/genotype_markers.hpp"
#include "panvar/genotype.hpp"
#include "panvar/genotype_index.hpp"
#include "panvar/graph_utils.hpp"
#include "panvar/genotype_reads.hpp"
#include "panvar/gfa.hpp"
#include "panvar/output.hpp"

#include "panvar/parallel.hpp"
#include "panvar/syncmer.hpp"

#include <unordered_set>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>

#include <iostream>
#include <stdexcept>

namespace panvar {
namespace {

void print_genotype_help() {
    std::cout
        << "Usage:\n"
        << "  panvar genotype -i <graph.gfa> -b <bubble-prefix> -r <ref-path> -o <out-prefix> --audit\n"
        << "\n"
        << "Genotype a short-read sample against a regional graph, per bubble.\n"
        << "\n"
        << "Only --audit is implemented so far: it builds the per-bubble syncmer graph and reports,\n"
        << "for every allele, how many discriminative markers exist in each candidate unit -- syncmer\n"
        << "nodes (present in this allele, absent from its siblings) and syncmer adjacencies (the\n"
        << "same rule on consecutive-syncmer pairs, which is what carries a deletion junction). The\n"
        << "two counts side by side decide which unit the genotyper should use.\n"
        << "\n"
        << "Options:\n"
        << "  -i, --gfa <path>            Input GFA (required)\n"
        << "  -b, --bubble-prefix-in <p>  Bubble prefix; reads <p>.bubbles.csv\n"
        << "  -c, --bubbles-csv-in <path> Bubbles CSV (alternative to --bubble-prefix-in)\n"
        << "  -r, --reference-path <name> Reference path name (required)\n"
        << "  -o, --out-prefix <path>     Output prefix (required)\n"
        << "      --audit                 Per-bubble marker feasibility audit\n"
        << "      --build-index <path>    Build the panel index and write it here, then exit. The\n"
        << "                              index depends only on the graph, so a cohort builds it once\n"
        << "      --index <path>          Genotype using a prebuilt index instead of rebuilding the\n"
        << "                              panel from the graph (much faster per sample)\n"
        << "  -R, --reads <path>          Short reads (FASTA/FASTQ, plain or gzipped); repeatable.\n"
        << "                              Projects them onto the block marker panel and reports\n"
        << "                              per-block depth. No genotypes are called yet.\n"
        << "      --min-anchors <N>       Minimum invariant markers for per-block depth (default 20)\n"
        << "      --marker-rule <r>       panvar (default) = keep markers whose multiplicity varies\n"
        << "                              across alleles; unique = carried by exactly one allele (any\n"
        << "                              copy number); pangenie = unique AND occurring once,\n"
        << "                              presence/absence only (PanGenie's actual rule);\n"
        << "                              mixed = presence/absence outside bubbles,\n"
        << "                              multiplicity inside them\n"
        << "      --fragment-len <N>      Library fragment length, used to discount correlated\n"
        << "                              markers when computing GQ (default 350; 0 disables)\n"
        << "      --carrier-weight <b>    Down-weight markers by how many of the block's alleles\n"
        << "                              carry them: weight = (n_alleles/carriers)^b, mean 1.\n"
        << "                              At blocks with hundreds of alleles the set is swamped\n"
        << "                              by markers shared across many of them, which\n"
        << "                              discriminate little but outvote the specific ones.\n"
        << "                              0 (default) disables\n"
        << "      --recomb-rate <x>       Li-Stephens switch scaling; 1.0 (default) is about one\n"
        << "                              expected haplotype switch across the locus. Raising it\n"
        << "                              makes blocks nearly independent, lowering it locks the\n"
        << "                              chain to one haplotype pair -- useful for telling\n"
        << "                              emission error apart from linkage error\n"
        << "      --provenance            Attribute each call to the blocks that determined it, by\n"
        << "                              neutralizing one block at a time and re-running the chain.\n"
        << "                              Adds provenance (self/neighbours/distant/none) and the\n"
        << "                              influencing block list to the genotypes table. Diagnostic:\n"
        << "                              caches every block's emissions, so memory grows with\n"
        << "                              n_blocks * n_haplotypes^2\n"
        << "      --max-alleles <N>       Candidate alleles kept per block before pairing (default\n"
        << "                              64). Rare alleles fall outside it and cannot be called\n"
        << "      --exclude-haplotypes <a,b> Drop these from the panel before genotyping. With\n"
        << "                              --truth-haplotypes this is the leave-one-out design: the\n"
        << "                              sample is no longer in the panel, which is the only\n"
        << "                              honest test, since identifying a haplotype that IS in the\n"
        << "                              panel is trivial and measures nothing\n"
        << "      --truth-haplotypes <a,b> Two panel haplotype names the reads came from; scores the\n"
        << "                              called allele pair against theirs, per block\n"
        << "      --audit-linkage         Block-chain + linkage identifiability audit: how far the\n"
        << "                              bubble/backbone chain narrows which panel haplotype a\n"
        << "                              sample is on. Run this on the pre-panphorte bubble graph\n"
        << "  -k, --kmer-size <N>         k for the syncmer markers (default 31, max 31)\n"
        << "      --syncmer-s <N>         s for the closed-syncmer test (0 = auto)\n"
        << "      --min-markers <N>       Per-allele threshold reported as usable (default 10)\n"
        << "      --max-multiplicity <N>  Drop markers repeated more than N times in an allele\n"
        << "                              (default 0 = no cap). Capping discards the copy-number\n"
        << "                              signal a tandem array carries: at LPA's KIV-2 block a cap\n"
        << "                              of 3 cuts median allele separation from 261 to 48\n"
        << "      --all-kmers             Use every k-mer instead of the ~17%% closed-syncmer sample\n"
        << "                              (~6x more markers, ~6x the memory and time)\n"
        << "      --no-region-unique      Skip the region-uniqueness filter (diagnostic)\n"
        << "      --max-dense-alleles <N> Above N alleles in a block, score pairwise separation with\n"
        << "                              an O(n) sparse accumulator instead of the O(n^2) dense one\n"
        << "                              (default 2048; identical results, ~40%% slower, but the\n"
        << "                              dense matrix would need 200 MB/block at 5000 alleles)\n"
        << "      --sep-top-k <N>         Approximate: score each allele only against its N most\n"
        << "                              similar siblings (MinHash), for panels where even the\n"
        << "                              sparse path is too slow. 0 = exact (default); compare\n"
        << "                              against an exact run before trusting a value. The error is\n"
        << "                              one-directional: it can only OVERSTATE separability\n"
        << "  -t, --threads <N>           Worker threads (0 = auto)\n"
        << "  -q, --quiet                 Disable progress logs\n"
        << "  -h, --help                  Show this help\n";
}

} // namespace

int run_genotype_command(const std::vector<std::string>& args) {
    if (args.empty()) {
        print_genotype_help();
        return 0;
    }

    std::string gfa_path;
    std::string bubble_prefix_in;
    std::string bubbles_csv_in;
    std::string reference_path;
    std::string out_prefix;
    bool audit = false;
    bool audit_linkage = false;
    std::vector<std::string> read_paths;
    std::size_t min_anchors = 20;
    std::string truth_haplotypes;
    std::string exclude_haplotypes;
    std::string index_out;
    std::string index_in;
    std::size_t max_alleles = 64;
    double fragment_len = 350.0;
    double recomb_rate = 1.0;
    double carrier_weight = 0.0;
    bool provenance = false;
    double uneven_tolerance = 0.35;
    bool quiet = false;
    MarkerOptions options;

    for (std::size_t i = 0; i < args.size(); ++i) {
        const std::string& arg = args[i];
        auto require_value = [&](const std::string& flag) -> const std::string& {
            if (i + 1 >= args.size()) throw std::runtime_error("Missing value after " + flag);
            return args[++i];
        };
        if (arg == "-h" || arg == "--help") { print_genotype_help(); return 0; }
        else if (arg == "-i" || arg == "--gfa") gfa_path = require_value(arg);
        else if (arg == "-b" || arg == "--bubble-prefix-in") bubble_prefix_in = require_value(arg);
        else if (arg == "-c" || arg == "--bubbles-csv-in") bubbles_csv_in = require_value(arg);
        else if (arg == "-r" || arg == "--reference-path") reference_path = require_value(arg);
        else if (arg == "-o" || arg == "--out-prefix") out_prefix = require_value(arg);
        else if (arg == "--audit") audit = true;
        else if (arg == "--audit-linkage") audit_linkage = true;
        else if (arg == "-R" || arg == "--reads") read_paths.push_back(require_value(arg));
        else if (arg == "--min-anchors") min_anchors = cli::parse_size_arg(arg, require_value(arg));
        else if (arg == "--truth-haplotypes") truth_haplotypes = require_value(arg);
        else if (arg == "--exclude-haplotypes") exclude_haplotypes = require_value(arg);
        else if (arg == "--build-index") index_out = require_value(arg);
        else if (arg == "--index") index_in = require_value(arg);
        else if (arg == "--marker-rule") {
            const std::string v = require_value(arg);
            if (v == "pangenie") options.rule = MarkerRule::PanGenie;
            else if (v == "unique") options.rule = MarkerRule::Unique;
            else if (v == "mixed") options.rule = MarkerRule::Mixed;
            else if (v != "panvar") throw std::runtime_error("genotype: --marker-rule must be panvar|unique|pangenie|mixed");
        }
        else if (arg == "--max-alleles") max_alleles = cli::parse_size_arg(arg, require_value(arg));
        else if (arg == "--fragment-len") fragment_len = std::stod(require_value(arg));
        else if (arg == "--recomb-rate") recomb_rate = std::stod(require_value(arg));
        else if (arg == "--carrier-weight") carrier_weight = std::stod(require_value(arg));
        else if (arg == "--provenance") provenance = true;
        else if (arg == "-k" || arg == "--kmer-size") options.kmer_size = cli::parse_size_arg(arg, require_value(arg));
        else if (arg == "--syncmer-s") options.syncmer_s = cli::parse_size_arg(arg, require_value(arg));
        else if (arg == "--min-markers") options.min_markers = cli::parse_size_arg(arg, require_value(arg));
        else if (arg == "--max-multiplicity") options.max_multiplicity = cli::parse_size_arg(arg, require_value(arg));
        else if (arg == "--no-region-unique") options.require_region_unique = false;
        else if (arg == "--all-kmers") options.all_kmers = true;
        else if (arg == "--max-dense-alleles") options.max_dense_alleles = cli::parse_size_arg(arg, require_value(arg));
        else if (arg == "--sep-top-k") options.sep_top_k = cli::parse_size_arg(arg, require_value(arg));
        else if (arg == "-t" || arg == "--threads") options.threads = cli::parse_size_arg(arg, require_value(arg));
        else if (arg == "-q" || arg == "--quiet") quiet = true;
        else throw std::runtime_error("Unknown option for genotype: " + arg);
    }

    if (gfa_path.empty() && index_in.empty()) {
        throw std::runtime_error("genotype requires -i/--gfa (or --index for a prebuilt panel)");
    }
    if (!index_in.empty() && read_paths.empty()) {
        throw std::runtime_error("genotype: --index is for genotyping; pass -R/--reads too");
    }
    if (out_prefix.empty()) throw std::runtime_error("genotype requires -o/--out-prefix");
    if (reference_path.empty() && index_in.empty()) {
        throw std::runtime_error("genotype requires -r/--reference-path");
    }
    if (options.kmer_size == 0 || options.kmer_size > 31) {
        throw std::runtime_error("genotype: --kmer-size must be in [1, 31]");
    }
    if (!audit && !audit_linkage && read_paths.empty() && index_out.empty()) {
        throw std::runtime_error(
            "genotype: pass --audit, --audit-linkage, --build-index, or -R/--reads");
    }
    // Reads need the block chain and the marker panel, but not the linkage/novelty audits or the
    // separation statistics -- those are diagnostics and were roughly half the runtime.
    const bool need_blocks = audit_linkage || !read_paths.empty() || !index_out.empty();
    const bool want_audit_stats = audit_linkage;
    if (!bubble_prefix_in.empty()) {
        if (!bubbles_csv_in.empty()) {
            throw std::runtime_error("genotype: use either --bubble-prefix-in or --bubbles-csv-in");
        }
        bubbles_csv_in = bubble_prefix_in + ".bubbles.csv";
    }
    if (bubbles_csv_in.empty() && index_in.empty()) {
        throw std::runtime_error("genotype requires --bubble-prefix-in or --bubbles-csv-in");
    }
    cli::ensure_parent_dir_for_file(out_prefix);

    cli::RunLog log("genotype", quiet);

    if (!index_in.empty()) {
        // Everything the panel contributes was precomputed; only the reads are new.
        const GenotypeIndex idx = read_genotype_index(index_in);
        log.info("index " + index_in + ": " + std::to_string(idx.chain.size()) + " blocks, " +
                 std::to_string(idx.haplotype_names.size()) + " haplotypes, " +
                 std::to_string(idx.panel.node_codes.size()) + " markers");
        const ReadCounts rc = count_reads(read_paths, idx.panel, options.threads);
        log.info("reads: " + std::to_string(rc.reads) + " (" + std::to_string(rc.bases / 1000) +
                 " kb); " + std::to_string(rc.syncmers) + " syncmers, " +
                 std::to_string(100 * rc.matched_syncmers / std::max<std::uint64_t>(1, rc.syncmers)) +
                 "% matched a panel marker");
        const std::vector<BlockDepth> depth =
            estimate_depth(idx.panel, rc, min_anchors, uneven_tolerance);
        GenotypeOptions gopt;
        gopt.threads = options.threads;
        gopt.max_alleles_per_block = max_alleles;
        gopt.fragment_len = fragment_len;
        GenotypeSummary gsum;
        const std::vector<BlockCall> calls = genotype_sample(idx.chain, idx.blocks, idx.panel, rc,
                                                             depth, idx.haplotype_names, gopt, &gsum);
        log.info("calls: " + std::to_string(gsum.called) + " PASS, " + std::to_string(gsum.no_calls) +
                 " no-call, " + std::to_string(gsum.off_panel) + " off-panel; mean GQ " +
                 std::to_string(gsum.mean_gq));
        write_read_audit(out_prefix, idx.chain, idx.panel, rc, depth);
        write_genotypes(out_prefix, idx.chain, idx.blocks, calls, idx.haplotype_names);
        log.wrote({out_prefix + ".reads.depth.tsv", out_prefix + ".genotypes.tsv",
                   out_prefix + ".haplotypes.tsv"});
        log.done();
        return 0;
    }
    ParseGfaOptions parse_options;
    parse_options.include_paths = true;
    parse_options.include_sequences = true;
    const Graph graph = parse_gfa(gfa_path, parse_options);
    if (graph.paths.empty()) {
        throw std::runtime_error("genotype: no paths in " + gfa_path);
    }
    // Leave-one-out: the excluded haplotypes are removed from the panel entirely, but their walks
    // are kept aside so the truth can still be resolved against the REDUCED panel's alleles. A
    // held-out haplotype whose block sequence no other haplotype carries is simply not representable
    // -- that is the mosaic ceiling, not a failure of the caller.
    Graph panel_graph = graph;
    std::vector<PathRecord> held_out;
    if (!exclude_haplotypes.empty()) {
        std::vector<std::string> names;
        std::string cur;
        for (const char ch : exclude_haplotypes) {
            if (ch == ',') { if (!cur.empty()) names.push_back(cur); cur.clear(); }
            else cur.push_back(ch);
        }
        if (!cur.empty()) names.push_back(cur);
        std::vector<PathRecord> keep;
        for (const PathRecord& p : panel_graph.paths) {
            if (std::find(names.begin(), names.end(), p.name) != names.end()) held_out.push_back(p);
            else keep.push_back(p);
        }
        if (held_out.size() != names.size()) {
            throw std::runtime_error("genotype: --exclude-haplotypes named a path not in the graph");
        }
        panel_graph.paths = std::move(keep);
    }

    const std::vector<Bubble> bubbles = read_bubbles_csv(bubbles_csv_in);
    log.info("input " + gfa_path + " (" + std::to_string(graph.nodes.size()) + " nodes, " +
             std::to_string(graph.paths.size()) + " paths); " + std::to_string(bubbles.size()) +
             " bubbles; reference " + reference_path);

    if (need_blocks) {
        const std::vector<Block> chain = build_block_chain(bubbles);
        std::vector<BubblePathIndex> path_indexes(panel_graph.paths.size());
        run_parallel(panel_graph.paths.size(), options.threads, [&](std::size_t i) {
            path_indexes[i] = build_bubble_path_index(panel_graph.paths[i]);
        });
        std::vector<BlockAlleles> blocks(chain.size());
        for (std::size_t i = 0; i < chain.size(); ++i) {
            blocks[i] = enumerate_block_alleles(panel_graph, path_indexes, bubbles, chain[i], options.threads);
        }
        std::size_t n_bb = 0;
        std::size_t n_fl = 0;
        for (const Block& b : chain) {
            if (b.kind == BlockKind::Backbone) ++n_bb;
            else if (b.kind == BlockKind::Flank) ++n_fl;
        }
        log.info("block chain: " + std::to_string(chain.size()) + " blocks (" +
                 std::to_string(chain.size() - n_bb - n_fl) + " bubble, " + std::to_string(n_bb) +
                 " backbone, " + std::to_string(n_fl) + " flank)");
        if (want_audit_stats) {
        const LinkageReport rep = measure_linkage(panel_graph, blocks);
        log.info("uniquely identified by the full chain: " + std::to_string(rep.uniquely_identified) +
                 "/" + std::to_string(rep.n_haplotypes) + " haplotypes");
        if (!rep.collapse_curve.empty()) {
            std::string curve;
            for (std::size_t i = 0; i < rep.collapse_curve.size(); ++i) {
                if (i && i % 4) continue;
                curve += (curve.empty() ? "" : " -> ") +
                         std::to_string(static_cast<long>(rep.collapse_curve[i] + 0.5));
            }
            log.info("mean compatible haplotypes: " + std::to_string(panel_graph.paths.size()) + " -> " + curve);
        }
        LinkageReport rep2 = rep;
        measure_leave_one_out(panel_graph, chain, blocks, rep2);
        log.info("leave-one-out ceiling: nearest remaining haplotype recovers " +
                 std::to_string(static_cast<long>(100.0 * rep2.loo_mean_agreement + 0.5)) +
                 "% of bubble alleles on average; " + std::to_string(rep2.loo_perfect) + "/" +
                 std::to_string(rep2.n_haplotypes) + " fully recovered");
        log.info("mosaic ceiling (any panel haplotype may supply each block): " +
                 std::to_string(static_cast<long>(100.0 * rep2.mosaic_ceiling + 0.5)) +
                 "% of bubble alleles recoverable; " +
                 std::to_string(static_cast<long>(100.0 * rep2.singleton_mass + 0.5)) +
                 "% of (haplotype, bubble) cells hold a private allele");
        write_linkage_audit(out_prefix, panel_graph, chain, blocks, rep2);
        }

        ReadPanel read_panel;
        const std::vector<BlockMarkerStats> bstats =
            build_block_marker_panel(chain, blocks, options,
                                     (read_paths.empty() && index_out.empty()) ? nullptr : &read_panel,
                                     (read_paths.empty() && index_out.empty()) ? nullptr : &panel_graph,
                                     want_audit_stats);
        if (want_audit_stats) {
            double tot = 0.0; std::size_t inf_n = 0, inf_e = 0, over = 0, maxmult = 0;
            double mass_n = 0.0, mass_e = 0.0; std::size_t hap_tot = 0;
            for (const BlockMarkerStats& b : bstats) {
                tot += b.seconds; inf_n += b.n_informative_nodes; inf_e += b.n_informative_edges;
                over += b.markers_over_cap; maxmult = std::max(maxmult, b.max_marker_multiplicity_seen);
                mass_n += b.separable_mass_nodes * static_cast<double>(b.n_haplotypes);
                mass_e += b.separable_mass_edges * static_cast<double>(b.n_haplotypes);
                hap_tot += b.n_haplotypes;
            }
            log.info("block markers: " + std::to_string(inf_n) + " informative nodes, " +
                     std::to_string(inf_e) + " informative adjacencies; separable haplotype mass at >=" +
                     std::to_string(options.min_markers) + " markers: nodes " +
                     std::to_string(static_cast<long>(100.0 * mass_n / std::max<std::size_t>(1, hap_tot) + 0.5)) +
                     "%, edges " +
                     std::to_string(static_cast<long>(100.0 * mass_e / std::max<std::size_t>(1, hap_tot) + 0.5)) + "%");
            log.info("multiplicity cap (--max-multiplicity " + std::to_string(options.max_multiplicity) +
                     ") would drop " + std::to_string(over) + " markers; highest multiplicity seen " +
                     std::to_string(maxmult) + "; marker scoring took " +
                     std::to_string(static_cast<long>(tot * 1000.0)) + " ms of CPU");
            write_block_marker_audit(out_prefix, bstats);
        }

        if (!index_out.empty()) {
            GenotypeIndex idx;
            idx.chain = chain;
            idx.blocks = blocks;
            for (BlockAlleles& b : idx.blocks) b.allele_seq.clear();   // not needed to genotype
            idx.panel = read_panel;
            idx.kmer_size = options.kmer_size;
            idx.syncmer_s = read_panel.syncmer_s;
            for (const PathRecord& p : panel_graph.paths) idx.haplotype_names.push_back(p.name);
            write_genotype_index(index_out, idx);
            log.info("index: " + std::to_string(idx.chain.size()) + " blocks, " +
                     std::to_string(idx.haplotype_names.size()) + " haplotypes, " +
                     std::to_string(idx.panel.node_codes.size()) + " markers");
            log.wrote({index_out});
            log.done();
            return 0;
        }

        if (!read_paths.empty()) {
            log.info("marker panel: " + std::to_string(read_panel.node_codes.size()) + " syncmers, " +
                     std::to_string(read_panel.edge_keys.size()) + " adjacencies; " +
                     std::to_string(read_panel.region_filtered_markers) +
                     " dropped for occurring elsewhere in the region (" +
                     std::to_string(read_panel.dropped_multi_block) + " appear in >1 block, " +
                     std::to_string(read_panel.dropped_over_expected) + " exceed what the blocks own, of " +
                     std::to_string(read_panel.informative_before_filter) + " informative)");
            log.info("confinement by context: 1-syncmer " +
                     std::to_string(100 * read_panel.confined_vary_nodes /
                                    std::max<std::size_t>(1, read_panel.vary_nodes)) +
                     "% (" + std::to_string(read_panel.confined_vary_nodes) + "/" +
                     std::to_string(read_panel.vary_nodes) + "), 2-syncmer " +
                     std::to_string(100 * read_panel.confined_vary_edges /
                                    std::max<std::size_t>(1, read_panel.vary_edges)) +
                     "% (" + std::to_string(read_panel.confined_vary_edges) + "/" +
                     std::to_string(read_panel.vary_edges) + ")");
            for (std::size_t b = 0; b < read_panel.block_overlap.size(); ++b) {
                if (read_panel.block_overlap[b].empty()) continue;
                std::string msg = "  block " + std::to_string(b) + " (" +
                                  (chain[b].kind == BlockKind::Bubble ? "bubble" :
                                   chain[b].kind == BlockKind::Flank ? "flank" : "backbone") +
                                  ") shares markers with:";
                for (const auto& [other, n] : read_panel.block_overlap[b]) {
                    msg += " " + std::to_string(other) + "(" + std::to_string(n) + ")";
                }
                log.info(msg);
            }
            log.info("multi-block markers: " + std::to_string(read_panel.dropped_adjacent_blocks) +
                     " span ADJACENT blocks only (shared boundary), " +
                     std::to_string(read_panel.dropped_distant_blocks) +
                     " span DISTANT blocks (duplication elsewhere in the region)");
            const ReadCounts rc = count_reads(read_paths, read_panel, options.threads);
            log.info("reads: " + std::to_string(rc.reads) + " (" + std::to_string(rc.bases / 1000) +
                     " kb); " + std::to_string(rc.syncmers) + " syncmers, " +
                     std::to_string(100 * rc.matched_syncmers / std::max<std::uint64_t>(1, rc.syncmers)) +
                     "% matched a panel marker; " + std::to_string(rc.novel_adjacencies) +
                     " novel adjacencies");
            const std::vector<BlockDepth> depth =
                estimate_depth(read_panel, rc, min_anchors, uneven_tolerance);
            std::size_t uneven = 0;
            std::vector<double> lam;
            for (const BlockDepth& d : depth) { if (d.uneven) ++uneven; if (d.usable) lam.push_back(d.lambda_hap); }
            std::sort(lam.begin(), lam.end());
            log.info("per-haplotype depth (lambda): median " +
                     std::to_string(lam.empty() ? 0.0 : lam[lam.size() / 2]) + " over " +
                     std::to_string(lam.size()) + " usable blocks; " + std::to_string(uneven) +
                     " flagged UNEVEN");
            write_read_audit(out_prefix, chain, read_panel, rc, depth);

            std::vector<std::string> hap_names;
            hap_names.reserve(panel_graph.paths.size());
            for (const PathRecord& p : panel_graph.paths) hap_names.push_back(p.name);
            GenotypeOptions gopt;
            gopt.threads = options.threads;
            gopt.max_alleles_per_block = max_alleles;
            gopt.fragment_len = fragment_len;
            gopt.provenance = provenance;
            gopt.recomb_rate = recomb_rate;
            gopt.carrier_weight = carrier_weight;
            GenotypeSummary gsum;
            std::vector<int> ta1;
            std::vector<int> ta2;
            std::vector<std::string> ts1;   // true spelled sequence per block, even when unrepresentable
            std::vector<std::string> ts2;
            if (!truth_haplotypes.empty()) {
                const auto comma = truth_haplotypes.find(',');
                if (comma == std::string::npos) {
                    throw std::runtime_error("genotype: --truth-haplotypes needs two comma-separated names");
                }
                const std::string a = truth_haplotypes.substr(0, comma);
                const std::string b = truth_haplotypes.substr(comma + 1);
                ta1.assign(chain.size(), -1);
                ta2.assign(chain.size(), -1);
                // Under leave-one-out the truth haplotypes are not in the panel, so their allele is
                // resolved by spelling their own walk and matching it against the reduced panel's
                // allele sequences. No match means the panel simply cannot represent that allele --
                // the mosaic ceiling -- and the block is scored as unrepresentable rather than wrong.
                auto resolve = [&](const std::string& name, std::vector<int>& out_alleles,
                                   std::vector<std::string>& out_seq) {
                    const auto direct = [&](std::size_t bi) {
                        const auto it = blocks[bi].allele_of.find(name);
                        return it == blocks[bi].allele_of.end() ? -1 : static_cast<int>(it->second);
                    };
                    const PathRecord* held = nullptr;
                    for (const PathRecord& p : held_out) if (p.name == name) held = &p;
                    if (held == nullptr) {
                        for (std::size_t bi = 0; bi < chain.size(); ++bi) {
                            out_alleles[bi] = direct(bi);
                            if (out_alleles[bi] >= 0 &&
                                static_cast<std::size_t>(out_alleles[bi]) < blocks[bi].allele_seq.size()) {
                                out_seq[bi] = blocks[bi].allele_seq[static_cast<std::size_t>(out_alleles[bi])];
                            }
                        }
                        return;
                    }
                    const BubblePathIndex idx = build_bubble_path_index(*held);
                    for (std::size_t bi = 0; bi < chain.size(); ++bi) {
                        std::optional<std::vector<PathStep>> st;
                        if (chain[bi].kind == BlockKind::Bubble) {
                            for (const Bubble& bb : bubbles) {
                                if (bb.id == chain[bi].bubble_id) { st = bubble_steps(*held, idx, bb); break; }
                            }
                        } else if (chain[bi].kind == BlockKind::Flank) {
                            const bool leading = chain[bi].source.empty();
                            st = flank_steps(*held, idx, leading ? chain[bi].sink : chain[bi].source, leading);
                        } else {
                            st = interval_interior_steps(*held, idx, chain[bi].source, chain[bi].sink);
                        }
                        if (!st.has_value() || st->empty()) continue;
                        std::string seq = spell_path_steps_sequence(graph, *st);
                        for (std::size_t ai = 0; ai < blocks[bi].allele_seq.size(); ++ai) {
                            if (blocks[bi].allele_seq[ai] == seq) { out_alleles[bi] = static_cast<int>(ai); break; }
                        }
                        // Kept even when no allele matches: an unrepresentable block still has a true
                        // sequence, and that is exactly the case the graded score exists to measure.
                        out_seq[bi] = std::move(seq);
                    }
                };
                ts1.assign(chain.size(), std::string());
                ts2.assign(chain.size(), std::string());
                resolve(a, ta1, ts1);
                resolve(b, ta2, ts2);
            }
            const std::vector<BlockCall> calls =
                genotype_sample(chain, blocks, read_panel, rc, depth, hap_names, gopt, &gsum,
                                ta1.empty() ? nullptr : &ta1, ta2.empty() ? nullptr : &ta2);
            log.info("model: lambda " + std::to_string(gsum.lambda_hap) + ", overdispersion phi " +
                     std::to_string(gsum.overdispersion) + ", error background " +
                     std::to_string(gsum.error_background));
            log.info("calls: " + std::to_string(gsum.called) + " PASS, " + std::to_string(gsum.no_calls) +
                     " no-call, " + std::to_string(gsum.off_panel) + " off-panel; mean GQ " +
                     std::to_string(gsum.mean_gq));
            write_genotypes(out_prefix, chain, blocks, calls, hap_names);

            if (!truth_haplotypes.empty()) {
                const auto comma = truth_haplotypes.find(',');
                if (comma == std::string::npos) {
                    throw std::runtime_error("genotype: --truth-haplotypes needs two comma-separated names");
                }
                const std::string t1 = truth_haplotypes.substr(0, comma);
                const std::string t2 = truth_haplotypes.substr(comma + 1);
                std::size_t ok = 0;
                std::size_t partial = 0;
                std::size_t wrong = 0;
                std::size_t scored = 0;
                std::size_t ok_bub = 0;
                std::size_t scored_bub = 0;
                std::size_t unrepresentable = 0;
                for (std::size_t bi = 0; bi < chain.size(); ++bi) {
                    if (ta1[bi] < 0 || ta2[bi] < 0) { ++unrepresentable; continue; }
                    const std::pair<const std::string, std::size_t> i1v{t1, static_cast<std::size_t>(ta1[bi])};
                    const std::pair<const std::string, std::size_t> i2v{t2, static_cast<std::size_t>(ta2[bi])};
                    const auto* i1 = &i1v;
                    const auto* i2 = &i2v;
                    ++scored;
                    const bool is_bub = chain[bi].kind == BlockKind::Bubble;
                    if (is_bub) ++scored_bub;
                    const std::size_t lo = std::min(i1->second, i2->second);
                    const std::size_t hi = std::max(i1->second, i2->second);
                    const std::size_t clo = std::min(calls[bi].allele1, calls[bi].allele2);
                    const std::size_t chi = std::max(calls[bi].allele1, calls[bi].allele2);
                    if (lo == clo && hi == chi) { ++ok; if (is_bub) ++ok_bub; }
                    else if (lo == clo || hi == chi || lo == chi || hi == clo) ++partial;
                    else ++wrong;
                }
                std::size_t em1 = 0;
                std::size_t empruned = 0;
                std::size_t emscored = 0;
                for (const BlockCall& c : calls) {
                    if (c.truth_emission_rank == -2) { ++empruned; ++emscored; }
                    else if (c.truth_emission_rank > 0) { ++emscored; if (c.truth_emission_rank == 1) ++em1; }
                }
                log.info("emission alone ranks the true pair first in " + std::to_string(em1) + "/" +
                         std::to_string(emscored) + " blocks; " + std::to_string(empruned) +
                         " had the true allele pruned before pairing");
                log.info("unrepresentable blocks (panel has no such allele): " +
                         std::to_string(unrepresentable) + "/" + std::to_string(chain.size()));
                log.info("truth check: " + std::to_string(ok) + "/" + std::to_string(scored) +
                         " blocks with the exact allele pair (" + std::to_string(partial) +
                         " one allele right, " + std::to_string(wrong) + " both wrong); bubble blocks " +
                         std::to_string(ok_bub) + "/" + std::to_string(scored_bub));

                // Graded accuracy. Exact allele identity is the strict test, but it discards every
                // block whose true allele the reduced panel cannot represent -- around 40% under
                // leave-one-out -- and it scores a near miss exactly like a wild one. Aligning the
                // called allele against the true sequence scores every block and says how wrong a
                // call is, on the same edit-distance/QV scale benchmark already uses.
                //
                // Reported alongside is the best any panel allele could have achieved. Without it the
                // metric flatters: at a block whose alleles are all but identical, any pick scores
                // high, and the number would say more about the block than about the caller.
                if (!blocks.empty() && !blocks[0].allele_seq.empty()) {
                    const std::string acc_path = out_prefix + ".accuracy.tsv";
                    std::ofstream acc(acc_path);
                    acc << "block_index\tblock_kind\tbubble_id\tn_alleles\trepresentable\texact\t"
                           "dbp\tbest_dbp\tidentity\tbest_identity\toracle_rank\tid_h1\tid_h2\trank_h1\trank_h2"
                           "\tcross_c1t2\tcross_c2t1\tqv\ttrue_bp\tcalled_bp\tfilter\n";
                    // Identity oracle: the most similar allele the panel could have offered. `best_dbp`
                    // only says some allele matched the truth's LENGTH, which among hundreds of alleles
                    // happens by coincidence -- it cannot distinguish "the panel had nothing better"
                    // from "the panel had something better and we missed it". Shortlist by syncmer
                    // Jaccard (cheap, over all alleles), then align the shortlist exactly.
                    const std::size_t kk = read_panel.kmer_size;
                    const std::size_t ss = read_panel.syncmer_s != 0 ? read_panel.syncmer_s
                                                                     : default_syncmer_s(kk);
                    auto sync_set = [&](const std::string& s) {
                        std::unordered_set<std::uint64_t> out;
                        for (const KmerOccurrence& o : collect_syncmers(s, kk, ss)) out.insert(o.code);
                        return out;
                    };
                    auto len_of = [&](std::size_t bi, std::size_t ai) {
                        return ai < blocks[bi].allele_seq.size() ? blocks[bi].allele_seq[ai].size() : 0;
                    };
                    double id_sum = 0.0, qv_sum = 0.0, id_sum_b = 0.0;
                    std::size_t graded = 0, graded_b = 0;
                    std::size_t dbp_sum = 0, best_sum = 0, dbp_sum_b = 0, best_sum_b = 0;
                    // Per-block identity averages blocks of wildly different sizes equally; the
                    // length-weighted figure is the one to compare across loci.
                    std::size_t tot_edits = 0, tot_aln = 0;
                    std::size_t dbp_rep = 0, best_rep = 0, graded_rep = 0;
                    std::size_t dbp_unrep = 0, graded_norep = 0;
                    for (std::size_t bi = 0; bi < chain.size(); ++bi) {
                        if (ts1[bi].empty() && ts2[bi].empty()) continue;
                        const std::string* tv[2] = {&ts1[bi], &ts2[bi]};
                        const std::size_t cv[2] = {calls[bi].allele1, calls[bi].allele2};
                        // Which called haplotype to compare against which true one is itself a choice,
                        // and length is a poor way to make it -- two alleles of equal length can be
                        // entirely different sequence. Align all four combinations and take the
                        // assignment with the smaller total edit distance, the same NW kernel the rest
                        // of the metric uses.
                        auto seq_of = [&](std::size_t ai, const std::string& fallback) -> const std::string& {
                            return ai < blocks[bi].allele_seq.size() ? blocks[bi].allele_seq[ai] : fallback;
                        };
                        NwAlign nwm[2][2];
                        for (int c = 0; c < 2; ++c) {
                            for (int t = 0; t < 2; ++t) {
                                nwm[c][t] = nw_edit_distance(seq_of(cv[c], *tv[t]), *tv[t]);
                            }
                        }
                        const int swap = (nwm[1][0].edits + nwm[0][1].edits) <
                                         (nwm[0][0].edits + nwm[1][1].edits) ? 1 : 0;
                        std::size_t dbp = 0, best = 0, tbp = 0, cbp = 0;
                        double idsum = 0.0, qvsum = 0.0;
                        for (int h = 0; h < 2; ++h) {
                            const std::string& truth = *tv[h];
                            const std::size_t ca = cv[h ^ swap];
                            const std::string& called = seq_of(ca, truth);
                            tbp += truth.size();
                            cbp += called.size();
                            dbp += static_cast<std::size_t>(
                                std::labs(static_cast<long>(called.size()) - static_cast<long>(truth.size())));
                            std::size_t bd = SIZE_MAX;
                            for (const std::string& s : blocks[bi].allele_seq) {
                                bd = std::min(bd, static_cast<std::size_t>(std::labs(
                                    static_cast<long>(s.size()) - static_cast<long>(truth.size()))));
                            }
                            best += bd == SIZE_MAX ? 0 : bd;
                            const NwAlign& nw = nwm[h ^ swap][h];
                            const double denom = static_cast<double>(std::max<std::size_t>(1, nw.aln_len));
                            idsum += 1.0 - static_cast<double>(nw.edits) / denom;
                            qvsum += -10.0 * std::log10(std::max(0.5, static_cast<double>(nw.edits)) / denom);
                            tot_edits += nw.edits;
                            tot_aln += nw.aln_len;
                        }
                        // Best identity any panel allele could have reached, and where our pick ranked.
                        double best_id = 0.0;
                        long oracle_rank = -1;
                        double id_h0 = 0.0, id_h1 = 0.0;
                        long rank_h0 = -1, rank_h1 = -1;
                        if (blocks[bi].allele_seq.size() <= 1024) {
                            std::vector<std::unordered_set<std::uint64_t>> asets;
                            asets.reserve(blocks[bi].allele_seq.size());
                            for (const std::string& s : blocks[bi].allele_seq) asets.push_back(sync_set(s));
                            double idsum_best = 0.0;
                            long ranksum = 0;
                            int nh_scored = 0;
                            // Per haplotype, not averaged: a mean hides the common failure where one
                            // haplotype is placed perfectly and the other is essentially arbitrary.
                            double id_h[2] = {0.0, 0.0};
                            long rank_h[2] = {-1, -1};
                            for (int h = 0; h < 2; ++h) {
                                const std::string& truth = *tv[h];
                                if (truth.empty()) continue;
                                const auto tset = sync_set(truth);
                                std::vector<std::pair<double, std::size_t>> jac;
                                jac.reserve(asets.size());
                                for (std::size_t ai = 0; ai < asets.size(); ++ai) {
                                    std::size_t inter = 0;
                                    const auto& A = asets[ai];
                                    const auto& S = A.size() <= tset.size() ? A : tset;
                                    const auto& L = A.size() <= tset.size() ? tset : A;
                                    for (const std::uint64_t c : S) if (L.count(c) != 0) ++inter;
                                    const std::size_t uni = A.size() + tset.size() - inter;
                                    jac.emplace_back(uni == 0 ? 1.0 : static_cast<double>(inter) /
                                                                          static_cast<double>(uni), ai);
                                }
                                std::sort(jac.begin(), jac.end(),
                                          [](const auto& x, const auto& y) {
                                              return x.first != y.first ? x.first > y.first
                                                                        : x.second < y.second;
                                          });
                                const std::size_t called_ai = cv[h ^ swap];
                                for (std::size_t r = 0; r < jac.size(); ++r) {
                                    if (jac[r].second == called_ai) {
                                        ranksum += static_cast<long>(r) + 1;
                                        rank_h[h] = static_cast<long>(r) + 1;
                                        break;
                                    }
                                }
                                double bid = 0.0;
                                for (std::size_t r = 0; r < std::min<std::size_t>(16, jac.size()); ++r) {
                                    const NwAlign n2 = nw_edit_distance(blocks[bi].allele_seq[jac[r].second], truth);
                                    const double d2 = static_cast<double>(std::max<std::size_t>(1, n2.aln_len));
                                    bid = std::max(bid, 1.0 - static_cast<double>(n2.edits) / d2);
                                }
                                idsum_best += bid;
                                id_h[h] = 1.0 - static_cast<double>(nwm[h ^ swap][h].edits) /
                                              static_cast<double>(std::max<std::size_t>(1, nwm[h ^ swap][h].aln_len));
                                ++nh_scored;
                            }
                            if (nh_scored > 0) {
                                best_id = idsum_best / nh_scored;
                                oracle_rank = ranksum / nh_scored;
                                id_h0 = id_h[0]; id_h1 = id_h[1];
                                rank_h0 = rank_h[0]; rank_h1 = rank_h[1];
                            }
                        }
                        const bool is_bub = chain[bi].kind == BlockKind::Bubble;
                        const bool repr = ta1[bi] >= 0 && ta2[bi] >= 0;
                        const bool exact = repr &&
                            std::min<std::size_t>(static_cast<std::size_t>(ta1[bi]), static_cast<std::size_t>(ta2[bi])) ==
                                std::min(calls[bi].allele1, calls[bi].allele2) &&
                            std::max<std::size_t>(static_cast<std::size_t>(ta1[bi]), static_cast<std::size_t>(ta2[bi])) ==
                                std::max(calls[bi].allele1, calls[bi].allele2);
                        acc << chain[bi].index << '\t'
                            << (is_bub ? "bubble" : chain[bi].kind == BlockKind::Flank ? "flank" : "backbone") << '\t'
                            << chain[bi].bubble_id << '\t' << blocks[bi].n_alleles << '\t'
                            << (repr ? 1 : 0) << '\t' << (exact ? 1 : 0) << '\t'
                            << dbp << '\t' << best << '\t' << (idsum / 2.0) << '\t' << best_id << '\t'
                            << oracle_rank << '\t' << id_h0 << '\t' << id_h1 << '\t'
                            << rank_h0 << '\t' << rank_h1 << '\t'
                            // Both called alleles matching the SAME truth haplotype is the failure a
                            // per-haplotype identity cannot show: distinct allele indices can still be
                            // near-identical sequence, leaving the other haplotype unexplained.
                            << (1.0 - static_cast<double>(nwm[0][1].edits) /
                                    static_cast<double>(std::max<std::size_t>(1, nwm[0][1].aln_len))) << '\t'
                            << (1.0 - static_cast<double>(nwm[1][0].edits) /
                                    static_cast<double>(std::max<std::size_t>(1, nwm[1][0].aln_len))) << '\t'
                            << (qvsum / 2.0) << '\t'
                            << tbp << '\t' << cbp << '\t' << calls[bi].filter << '\n';
                        id_sum += idsum / 2.0; qv_sum += qvsum / 2.0; ++graded;
                        dbp_sum += dbp; best_sum += best;
                        if (is_bub) { id_sum_b += idsum / 2.0; ++graded_b; dbp_sum_b += dbp; best_sum_b += best; }
                        // Scoring a block the caller declined to call measures the metric, not the
                        // caller. Reported separately so a refused call cannot be counted as an error.
                        const bool reported = calls[bi].filter == "PASS" || calls[bi].filter == "LINKED";
                        if (reported) { dbp_rep += dbp; best_rep += best; ++graded_rep; }
                        else { dbp_unrep += dbp; ++graded_norep; }
                    }
                    if (graded > 0) {
                        auto pct = [](double v) { return std::to_string(100.0 * v); };
                        log.info("graded accuracy over " + std::to_string(graded) + " blocks (vs " +
                                 std::to_string(scored) + " scored exactly): mean identity " +
                                 pct(id_sum / static_cast<double>(graded)) + "%, mean QV " +
                                 std::to_string(qv_sum / static_cast<double>(graded)) + ", total dbp " +
                                 std::to_string(dbp_sum) + " (best any panel allele could do: " +
                                 std::to_string(best_sum) + "); length-weighted identity " +
                                 pct(1.0 - static_cast<double>(tot_edits) /
                                               static_cast<double>(std::max<std::size_t>(1, tot_aln))) + "%");
                        log.info("graded accuracy, reported calls only (PASS/LINKED, " +
                                 std::to_string(graded_rep) + " blocks): total dbp " +
                                 std::to_string(dbp_rep) + " (best: " + std::to_string(best_rep) +
                                 "); the " + std::to_string(graded_norep) +
                                 " blocks the caller declined carry dbp " + std::to_string(dbp_unrep));
                        if (graded_b > 0) {
                            log.info("graded accuracy, bubble blocks only (" + std::to_string(graded_b) +
                                     "): mean identity " + pct(id_sum_b / static_cast<double>(graded_b)) +
                                     "%, total dbp " + std::to_string(dbp_sum_b) + " (best: " +
                                     std::to_string(best_sum_b) + ")");
                        }
                    }
                    log.wrote({acc_path});
                }
            }
        }

        if (want_audit_stats) {
        const NoveltyReport nov = measure_novelty(chain, blocks, options.kmer_size, options.syncmer_s);
        log.info("private bubble alleles: " + std::to_string(nov.private_alleles) + "; syncmer content reused elsewhere in the block " +
                 std::to_string(static_cast<long>(100.0 * nov.mean_syncmer_reuse + 0.5)) + "%, adjacencies " +
                 std::to_string(static_cast<long>(100.0 * nov.mean_adjacency_reuse + 0.5)) + "%; " +
                 std::to_string(nov.fully_reusable) + " are pure rearrangements, " +
                 std::to_string(nov.with_novel_sequence) + " carry novel sequence");
        log.wrote({out_prefix + ".audit.blocks.tsv", out_prefix + ".audit.linkage.tsv",
                   out_prefix + ".audit.blockmarkers.tsv"});
        }
        if (!read_paths.empty()) {
            log.wrote({out_prefix + ".reads.depth.tsv", out_prefix + ".genotypes.tsv",
                       out_prefix + ".haplotypes.tsv"});
        }
        if (!audit) { log.done(); return 0; }
    }

    const MarkerPanel panel = build_marker_panel(graph, bubbles, reference_path, options);

    std::size_t green_nodes = 0;
    std::size_t green_edges = 0;
    std::size_t scored = 0;
    for (const AlleleMarkers& m : panel.markers) {
        if (m.n_haplotypes < 2) continue;   // singletons are not validatable; report but do not score
        ++scored;
        if (m.min_separating_nodes >= options.min_markers) ++green_nodes;
        if (m.min_separating_edges >= options.min_markers) ++green_edges;
    }
    log.info("k=" + std::to_string(panel.kmer_size) + " s=" + std::to_string(panel.syncmer_s) +
             "; " + std::to_string(panel.markers.size()) + " alleles over " +
             std::to_string(panel.bubbles.size()) + " bubbles");
    log.info("alleles with >=2 haplotypes: " + std::to_string(scored) + "; reaching --min-markers " +
             std::to_string(options.min_markers) + ": nodes " + std::to_string(green_nodes) +
             ", edges " + std::to_string(green_edges));

    write_marker_audit(out_prefix, panel);
    log.wrote({out_prefix + ".audit.bubbles.tsv", out_prefix + ".audit.alleles.tsv"});
    log.done();
    return 0;
}

} // namespace panvar
