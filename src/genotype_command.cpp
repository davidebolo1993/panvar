#include "panvar/genotype_command.hpp"

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

#include <algorithm>

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
        << "                              presence/absence only (PanGenie's actual rule)\n"
        << "      --fragment-len <N>      Library fragment length, used to discount correlated\n"
        << "                              markers when computing GQ (default 350; 0 disables)\n"
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
            else if (v != "panvar") throw std::runtime_error("genotype: --marker-rule must be panvar|unique|pangenie");
        }
        else if (arg == "--max-alleles") max_alleles = cli::parse_size_arg(arg, require_value(arg));
        else if (arg == "--fragment-len") fragment_len = std::stod(require_value(arg));
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
            GenotypeSummary gsum;
            std::vector<int> ta1;
            std::vector<int> ta2;
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
                auto resolve = [&](const std::string& name, std::vector<int>& out_alleles) {
                    const auto direct = [&](std::size_t bi) {
                        const auto it = blocks[bi].allele_of.find(name);
                        return it == blocks[bi].allele_of.end() ? -1 : static_cast<int>(it->second);
                    };
                    const PathRecord* held = nullptr;
                    for (const PathRecord& p : held_out) if (p.name == name) held = &p;
                    if (held == nullptr) {
                        for (std::size_t bi = 0; bi < chain.size(); ++bi) out_alleles[bi] = direct(bi);
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
                            st = interval_steps(*held, idx, chain[bi].source, chain[bi].sink);
                        }
                        if (!st.has_value() || st->empty()) continue;
                        const std::string seq = spell_path_steps_sequence(graph, *st);
                        for (std::size_t ai = 0; ai < blocks[bi].allele_seq.size(); ++ai) {
                            if (blocks[bi].allele_seq[ai] == seq) { out_alleles[bi] = static_cast<int>(ai); break; }
                        }
                    }
                };
                resolve(a, ta1);
                resolve(b, ta2);
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
