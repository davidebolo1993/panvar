#include "panvar/genotype_frag_command.hpp"

#include "panvar/bubble_path.hpp"
#include "panvar/bubbles.hpp"
#include "panvar/cli_utils.hpp"
#include "panvar/genotype_blocks.hpp"
#include "panvar/genotype_fragments.hpp"
#include "panvar/gfa.hpp"
#include "panvar/output.hpp"
#include "panvar/parallel.hpp"

#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace panvar {

namespace {

void print_help() {
    std::cout
        << "Usage:\n"
        << "  panvar genotype-frag -i <graph.gfa> -b <bubble-prefix> -o <out-prefix>\n"
        << "                       -R <reads.fq.gz> [-R <reads2.fq.gz>]\n"
        << "\n"
        << "PROTOTYPE. Arm F4 of experiments/FRAGMENT_EVIDENCE_PREREGISTRATION.md, built to be run\n"
        << "beside `panvar genotype` rather than instead of it. No production default moves.\n"
        << "\n"
        << "The difference under test is the UNIT OF OBSERVATION. `genotype` reduces the reads to a\n"
        << "vector of marker counts and then reconstructs the lost read linkage through confinement,\n"
        << "over-expected filtering, clump correction and an ESS discount. This scores each PHYSICAL\n"
        << "FRAGMENT once, against the spelled sequence of each candidate allele plus its flanks:\n"
        << "\n"
        << "  log L(a,b) = SUM over fragments log( 0.5 P(frag|a) + 0.5 P(frag|b) )\n"
        << "\n"
        << "P(frag|allele) combines both mates' infix alignment edit distance and, when both mates\n"
        << "place, the implied insert length. Two consequences that the marker model cannot reach:\n"
        << "a marker shared between two blocks is localised by the block-specific marker on the SAME\n"
        << "fragment rather than deleted by a filter, and an allele shorter than k -- which carries\n"
        << "no k-mer at all and is invisible to the emission -- is still voted on, because the\n"
        << "fragment spans it and both its boundaries.\n"
        << "\n"
        << "NOT IMPLEMENTED here, by design: no block chain, no HMM, no linkage between blocks, no\n"
        << "depth/dosage channel. Each block is scored on its own fragments. Those layers are only\n"
        << "worth building if this one passes its gate.\n"
        << "\n"
        << "Options:\n"
        << "  -i, --gfa <path>            Input GFA (required)\n"
        << "  -b, --bubble-prefix-in <p>  Bubble prefix; reads <p>.bubbles.csv\n"
        << "  -c, --bubbles-csv-in <path> Bubbles CSV (alternative to --bubble-prefix-in)\n"
        << "  -o, --out-prefix <path>     Output prefix (required)\n"
        << "  -R, --reads <path>          FASTA/FASTQ, plain or gzipped; repeatable. Mates are joined\n"
        << "                              by name, so interleaved and split R1/R2 both work\n"
        << "      --haplotype-mode        Score whole panel HAPLOTYPE pairs end to end and project\n"
        << "                              the answer onto blocks, instead of scoring each block on\n"
        << "                              its own candidates. Measured reason it exists, at cyp2d6\n"
        << "                              block 5 under leave-ZERO-out: the panel's allele 7 (556 bp)\n"
        << "                              contains the truth label's allele 0 (261 bp) and both begin\n"
        << "                              at the same position in the sample's own haplotype -- the\n"
        << "                              rest of the sample's real sequence lives in the NEIGHBOURING\n"
        << "                              block. Any block-local context must guess a flank there, the\n"
        << "                              sample's flank is unknown, and a fragment spanning the\n"
        << "                              junction then prefers whichever candidate supplies more real\n"
        << "                              sequence before the guess begins. The evidence that decides\n"
        << "                              the block sits across the block boundary, so the block is\n"
        << "                              the wrong unit to score. In this mode nothing is guessed:\n"
        << "                              every haplotype carries its own flank\n"
        << "                              LIMIT: a pair of complete panel haplotypes cannot represent\n"
        << "                              a mosaic, so under leave-one-out this is bounded by the\n"
        << "                              panel's mosaic ceiling. Measuring that bound is the point\n"
        << "      --max-haplotypes <N>    Shortlist size in --haplotype-mode (default 48)\n"
        << "      --probe-haplotypes <a,b>  Report the likelihood's opinion of this named haplotype\n"
        << "                              pair: whether it was shortlisted at all, its rank, its\n"
        << "                              score delta to the best pair, and how many fragments placed\n"
        << "                              on each. Separates three failures a truncated top-N list\n"
        << "                              cannot -- never shortlisted, shortlisted and scored badly,\n"
        << "                              or scored and close. Repeatable; costs nothing, the pair\n"
        << "                              scores already exist\n"
        << "      --max-anchor-occ <N>    Occurrences PER HAPLOTYPE above which a syncmer anchors\n"
        << "                              nothing (default 8): inside a tandem array it points\n"
        << "                              everywhere. Per haplotype, not across the shortlist -- a\n"
        << "                              shared cap makes anchoring depend on shortlist size, and\n"
        << "                              it did: at cyp2d6 NA18939 the same reads put the truth at\n"
        << "                              rank 2 with 48 haplotypes and rank 1 with 96\n"
        << "      --anchor-slack <N>      Window either side of an anchored read start (default 40)\n"
        << "      --marginalise-placements  Sum a fragment's likelihood over EVERY placement it has\n"
        << "                              on a haplotype instead of taking its best one. A fragment\n"
        << "                              compatible with ten positions is evidence for a haplotype\n"
        << "                              offering ten, and a maximum scores that identically to a\n"
        << "                              haplotype offering one -- which is worst exactly where\n"
        << "                              cyp2d6 is hardest. Implies --length-normalize, because a\n"
        << "                              sum without the 1/N factor rewards a repetitive haplotype\n"
        << "                              for offering more places to land\n"
        << "      --placement-topk <N>    Implied-start bins kept per mate per haplotype (default 2).\n"
        << "                              A sum over placements is only meaningful if the placements\n"
        << "                              are enumerated, so raise it with the flag above\n"
        << "      --coverage-weight <w>   Weight on the depth channel in --haplotype-mode (default\n"
        << "                              1.0; 0 disables). Sequence compatibility and copy number are\n"
        << "                              different signals: reads from a duplicated segment align\n"
        << "                              perfectly to a single copy, so alignment alone cannot see\n"
        << "                              dosage. This scores each haplotype's own length in windows,\n"
        << "                              with the rate fitted per haplotype. It reuses the same\n"
        << "                              fragments as the sequence term, so at weight 1 the two\n"
        << "                              double count -- which is why the default is 0, not 1:\n"
        << "                              measured at cyp2d6 leave-ZERO-out, weight 1 moved NA18939's\n"
        << "                              own true pair from rank 1 to rank 2 and cost 8 of its 19\n"
        << "                              blocks. Real mechanism, not yet a safe default\n"
        << "      --coverage-window <N>   Window size for that channel (default 500)\n"
        << "      --blocks <a,b,c>        Score only these chain indices. Default: every bubble block\n"
        << "      --all-blocks            Score backbone and flank blocks too. gstm1's worst blocks\n"
        << "                              are BACKBONE blocks carrying 16 and 21 alleles, so the\n"
        << "                              default view omits them\n"
        << "      --truth-haplotypes <a,b>  Two panel haplotype names the reads came from. Adds the\n"
        << "                              truth pair's rank, tie count and delta per block\n"
        << "      --exclude-haplotypes <a,b>  Drop these from the panel first. With\n"
        << "                              --truth-haplotypes this is leave-one-out; the truth allele\n"
        << "                              is then resolved by spelling the held-out walk and matching\n"
        << "                              it against the reduced panel, exactly as `genotype` does\n"
        << "      --flank-bp <N>          Neighbouring sequence glued to each side of a candidate\n"
        << "                              (default 500). Identical for every candidate of a block, so\n"
        << "                              it buys reachability for boundary-spanning fragments and\n"
        << "                              cannot shift a comparison\n"
        << "      --max-alleles <N>       Candidates kept per block after coarse containment\n"
        << "                              shortlisting (default 64). truth_rank -2 means the truth\n"
        << "                              was representable but dropped here\n"
        << "      --min-hits <N>          Syncmer hits needed to recruit a fragment to a block\n"
        << "                              (default 2)\n"
        << "      --error-rate <p>        Per-base edit probability (default 0.01)\n"
        << "      --bg-divergence <d>     Divergence at which a fragment stops being able to prefer\n"
        << "                              anything (default 0.10). Alignment likelihood is unbounded\n"
        << "                              below, so without this one fragment from sequence no\n"
        << "                              candidate models outvotes hundreds of ordinary ones -- and\n"
        << "                              under leave-one-out every candidate lacks some of the\n"
        << "                              sample's sequence, so that asymmetry has a direction\n"
        << "      --frag-outlier <e>      Prior weight on that background component (default 1e-3)\n"
        << "      --fragment-len <N>      Library mean insert (default 350)\n"
        << "      --fragment-sd <N>       Its standard deviation (default 50)\n"
        << "      --no-insert-size        Score sequence compatibility only. The channel separation\n"
        << "                              matters: run both to say which one carries a result\n"
        << "      --no-length-normalize   Drop the 1/L factor. Without it P(fragment|allele) is not\n"
        << "                              a generative likelihood and a SUPERSET allele is never\n"
        << "                              penalised -- every read from a short allele also fits an\n"
        << "                              allele containing it plus extra, and the extra is free.\n"
        << "                              Measured: leaving it off is what made leave-zero-out fail\n"
        << "      --no-compete            Score the block in isolation instead of against the rest\n"
        << "                              of the locus. A block scored alone must explain every\n"
        << "                              fragment handed to it, so a read from a paralogous copy\n"
        << "                              elsewhere votes here -- at CYP2D6 for whichever candidate\n"
        << "                              is most paralog-like\n"
        << "      --max-divergence <d>    Alignment band as a fraction of read length (default 0.20)\n"
        << "      --top-pairs <N>         Pairs written per block to <prefix>.frag_pairs.tsv\n"
        << "                              (default 20)\n"
        << "      --debug-block <N>       Write the per-fragment, per-candidate log-likelihood and\n"
        << "                              edit distance for chain block N to\n"
        << "                              <prefix>.frag_debug.tsv. A block score is a sum over\n"
        << "                              hundreds of fragments and nothing about it is checkable\n"
        << "                              from the sum\n"
        << "  -k, --kmer-size <N>         k for recruitment syncmers (default 31)\n"
        << "      --syncmer-s <N>         s for the closed-syncmer test (0 = auto)\n"
        << "  -t, --threads <N>           Worker threads (0 = auto)\n"
        << "  -q, --quiet                 Disable progress logs\n"
        << "  -h, --help                  Show this help\n";
}

std::vector<std::string> split_commas(const std::string& text) {
    std::vector<std::string> out;
    std::string cur;
    for (const char ch : text) {
        if (ch == ',') { if (!cur.empty()) out.push_back(cur); cur.clear(); }
        else cur.push_back(ch);
    }
    if (!cur.empty()) out.push_back(cur);
    return out;
}

} // namespace

int run_genotype_frag_command(const std::vector<std::string>& args) {
    if (args.empty()) { print_help(); return 0; }

    std::string gfa_path, bubble_prefix_in, bubbles_csv_in, out_prefix;
    std::vector<std::string> read_paths;
    std::string truth_haplotypes, exclude_haplotypes, blocks_arg;
    bool all_blocks = false, quiet = false, hap_mode = false, length_normalize_set = false;
    std::size_t top_pairs = 20;
    HaplotypeScoreOptions hopt;
    FragmentScoreOptions& opt = hopt;

    const auto value = [&](std::size_t& i, const std::string& name) {
        if (i + 1 >= args.size()) throw std::runtime_error("genotype-frag: " + name + " needs a value");
        return args[++i];
    };

    for (std::size_t i = 0; i < args.size(); ++i) {
        const std::string& a = args[i];
        if (a == "-h" || a == "--help") { print_help(); return 0; }
        else if (a == "-i" || a == "--gfa") gfa_path = value(i, a);
        else if (a == "-b" || a == "--bubble-prefix-in") bubble_prefix_in = value(i, a);
        else if (a == "-c" || a == "--bubbles-csv-in") bubbles_csv_in = value(i, a);
        else if (a == "-o" || a == "--out-prefix") out_prefix = value(i, a);
        else if (a == "-R" || a == "--reads") read_paths.push_back(value(i, a));
        else if (a == "--blocks") blocks_arg = value(i, a);
        else if (a == "--all-blocks") all_blocks = true;
        else if (a == "--truth-haplotypes") truth_haplotypes = value(i, a);
        else if (a == "--exclude-haplotypes") exclude_haplotypes = value(i, a);
        else if (a == "--flank-bp") opt.flank_bp = cli::parse_size_arg(a, value(i, a));
        else if (a == "--max-alleles") opt.max_alleles = cli::parse_size_arg(a, value(i, a));
        else if (a == "--min-hits") opt.min_recruit_hits = cli::parse_size_arg(a, value(i, a));
        else if (a == "--error-rate") opt.error_rate = std::stod(value(i, a));
        else if (a == "--bg-divergence") opt.bg_divergence = std::stod(value(i, a));
        else if (a == "--frag-outlier") opt.outlier_mix = std::stod(value(i, a));
        else if (a == "--fragment-len") opt.fragment_len = std::stod(value(i, a));
        else if (a == "--fragment-sd") opt.fragment_sd = std::stod(value(i, a));
        else if (a == "--no-insert-size") opt.use_insert_size = false;
        else if (a == "--no-length-normalize") { opt.length_normalize = false; length_normalize_set = true; }
        else if (a == "--length-normalize") { opt.length_normalize = true; length_normalize_set = true; }
        else if (a == "--no-compete") opt.compete = false;
        else if (a == "--haplotype-mode") hap_mode = true;
        else if (a == "--max-haplotypes") hopt.max_haplotypes = cli::parse_size_arg(a, value(i, a));
        else if (a == "--max-anchor-occ") hopt.max_anchor_occ = cli::parse_size_arg(a, value(i, a));
        else if (a == "--anchor-slack") hopt.anchor_slack = cli::parse_size_arg(a, value(i, a));
        else if (a == "--probe-haplotypes") {
            const std::vector<std::string> two = split_commas(value(i, a));
            if (two.size() != 2) {
                throw std::runtime_error("genotype-frag: --probe-haplotypes needs exactly two names");
            }
            hopt.probe_pairs.emplace_back(two[0], two[1]);
        }
        else if (a == "--marginalise-placements") hopt.marginalise_placements = true;
        else if (a == "--placement-topk") hopt.placement_topk = cli::parse_size_arg(a, value(i, a));
        else if (a == "--coverage-weight") hopt.coverage_weight = std::stod(value(i, a));
        else if (a == "--coverage-window") hopt.coverage_window = cli::parse_size_arg(a, value(i, a));
        else if (a == "--max-divergence") opt.max_divergence = std::stod(value(i, a));
        else if (a == "--top-pairs") top_pairs = cli::parse_size_arg(a, value(i, a));
        else if (a == "--debug-block") opt.debug_block =
            static_cast<long>(cli::parse_size_arg(a, value(i, a)));
        else if (a == "-k" || a == "--kmer-size") opt.kmer_size = cli::parse_size_arg(a, value(i, a));
        else if (a == "--syncmer-s") opt.syncmer_s = cli::parse_size_arg(a, value(i, a));
        else if (a == "-t" || a == "--threads") opt.threads = cli::parse_size_arg(a, value(i, a));
        else if (a == "-q" || a == "--quiet") quiet = true;
        else throw std::runtime_error("genotype-frag: unknown option " + a);
    }

    if (gfa_path.empty() || out_prefix.empty()) {
        throw std::runtime_error("genotype-frag requires --gfa and --out-prefix");
    }
    if (read_paths.empty()) throw std::runtime_error("genotype-frag requires at least one --reads");
    if (!bubble_prefix_in.empty()) {
        if (!bubbles_csv_in.empty()) {
            throw std::runtime_error("genotype-frag: use either --bubble-prefix-in or --bubbles-csv-in");
        }
        bubbles_csv_in = bubble_prefix_in + ".bubbles.csv";
    }
    if (bubbles_csv_in.empty()) {
        throw std::runtime_error("genotype-frag requires --bubble-prefix-in or --bubbles-csv-in");
    }
    if (opt.error_rate <= 0.0 || opt.error_rate >= 1.0) {
        throw std::runtime_error("genotype-frag: --error-rate must be in (0,1)");
    }
    if (opt.bg_divergence <= 0.0 || opt.bg_divergence >= 1.0) {
        throw std::runtime_error("genotype-frag: --bg-divergence must be in (0,1)");
    }
    if (opt.outlier_mix <= 0.0 || opt.outlier_mix >= 1.0) {
        throw std::runtime_error("genotype-frag: --frag-outlier must be in (0,1)");
    }
    if (opt.fragment_sd <= 0.0) throw std::runtime_error("genotype-frag: --fragment-sd must be > 0");
    if (opt.kmer_size == 0 || opt.kmer_size > 31) {
        throw std::runtime_error("genotype-frag: --kmer-size must be in [1,31]");
    }
    if (opt.debug_block >= 0) opt.debug_path = out_prefix + ".frag_debug.tsv";
    cli::ensure_parent_dir_for_file(out_prefix);

    cli::RunLog log("genotype-frag", quiet);

    ParseGfaOptions parse_options;
    parse_options.include_paths = true;
    parse_options.include_sequences = true;
    const Graph graph = parse_gfa(gfa_path, parse_options);
    if (graph.paths.empty()) throw std::runtime_error("genotype-frag: no paths in " + gfa_path);

    // Leave-one-out, with the same contract as `genotype`: the named haplotypes leave the panel but
    // their walks are kept, so the truth allele can still be resolved against the REDUCED panel. A
    // block where no remaining haplotype spells the held-out sequence is unrepresentable -- the
    // mosaic ceiling -- and is reported as truth_rank -1 rather than counted as a wrong call.
    Graph panel_graph = graph;
    std::vector<PathRecord> held_out;
    if (!exclude_haplotypes.empty()) {
        const std::vector<std::string> names = split_commas(exclude_haplotypes);
        std::vector<PathRecord> keep;
        for (const PathRecord& p : panel_graph.paths) {
            if (std::find(names.begin(), names.end(), p.name) != names.end()) held_out.push_back(p);
            else keep.push_back(p);
        }
        if (held_out.size() != names.size()) {
            throw std::runtime_error("genotype-frag: --exclude-haplotypes named a path not in the graph");
        }
        panel_graph.paths = std::move(keep);
    }

    const std::vector<Bubble> bubbles = read_bubbles_csv(bubbles_csv_in);
    const std::vector<Block> chain = build_block_chain(bubbles);
    std::vector<BubblePathIndex> path_indexes(panel_graph.paths.size());
    run_parallel(panel_graph.paths.size(), opt.threads, [&](std::size_t i) {
        path_indexes[i] = build_bubble_path_index(panel_graph.paths[i]);
    });
    std::vector<BlockAlleles> blocks(chain.size());
    for (std::size_t i = 0; i < chain.size(); ++i) {
        blocks[i] = enumerate_block_alleles(panel_graph, path_indexes, bubbles, chain[i], opt.threads);
    }
    log.info("input " + gfa_path + " (" + std::to_string(graph.nodes.size()) + " nodes, " +
             std::to_string(graph.paths.size()) + " paths); " + std::to_string(bubbles.size()) +
             " bubbles; chain of " + std::to_string(chain.size()) + " blocks" +
             (held_out.empty() ? "" : "; " + std::to_string(held_out.size()) + " held out"));

    // ---- truth, resolved the same way `genotype` resolves it ---------------------------------
    std::vector<int> truth1, truth2;
    bool have_truth = false;
    if (!truth_haplotypes.empty()) {
        const std::vector<std::string> names = split_commas(truth_haplotypes);
        if (names.size() != 2) {
            throw std::runtime_error("genotype-frag: --truth-haplotypes needs exactly two names");
        }
        std::vector<BlockAlleles> held_blocks(chain.size());
        if (!held_out.empty()) {
            Graph held_graph = graph;
            held_graph.paths = held_out;
            std::vector<BubblePathIndex> held_idx(held_graph.paths.size());
            for (std::size_t kk = 0; kk < held_graph.paths.size(); ++kk) {
                held_idx[kk] = build_bubble_path_index(held_graph.paths[kk]);
            }
            for (std::size_t bi = 0; bi < chain.size(); ++bi) {
                held_blocks[bi] = enumerate_block_alleles(held_graph, held_idx, bubbles,
                                                          chain[bi], opt.threads);
            }
        }
        const auto resolve = [&](const std::string& name, std::vector<int>& out_alleles) {
            out_alleles.assign(chain.size(), -1);
            const bool is_held =
                std::any_of(held_out.begin(), held_out.end(),
                            [&](const PathRecord& p) { return p.name == name; });
            for (std::size_t bi = 0; bi < chain.size(); ++bi) {
                if (!is_held) {
                    const auto it = blocks[bi].allele_of.find(name);
                    if (it != blocks[bi].allele_of.end()) out_alleles[bi] = static_cast<int>(it->second);
                    continue;
                }
                const auto it = held_blocks[bi].allele_of.find(name);
                if (it == held_blocks[bi].allele_of.end()) continue;
                const std::size_t hai = it->second;
                if (held_blocks[bi].bypass_allele >= 0 &&
                    hai == static_cast<std::size_t>(held_blocks[bi].bypass_allele)) {
                    if (blocks[bi].bypass_allele >= 0) out_alleles[bi] = blocks[bi].bypass_allele;
                    continue;
                }
                if (hai >= held_blocks[bi].allele_seq.size()) continue;
                const std::string& seq = held_blocks[bi].allele_seq[hai];
                if (seq.empty()) continue;
                for (std::size_t ai = 0; ai < blocks[bi].allele_seq.size(); ++ai) {
                    if (blocks[bi].allele_seq[ai] == seq) {
                        out_alleles[bi] = static_cast<int>(ai);
                        break;
                    }
                }
            }
        };
        resolve(names[0], truth1);
        resolve(names[1], truth2);
        have_truth = true;
        std::size_t representable = 0;
        for (std::size_t bi = 0; bi < chain.size(); ++bi) {
            if (truth1[bi] >= 0 && truth2[bi] >= 0) ++representable;
        }
        log.info("truth " + names[0] + " / " + names[1] + ": both alleles representable at " +
                 std::to_string(representable) + "/" + std::to_string(chain.size()) + " blocks");
    }

    if (hap_mode) {
        // Summing over placements without dividing by the number of positions a haplotype offers
        // rewards a repetitive haplotype for offering more places to land. The two go together.
        if (!length_normalize_set) hopt.length_normalize = hopt.marginalise_placements;
        FragmentLoadStats hs;
        const std::vector<Fragment> frags = load_fragments(read_paths, &hs);
        log.info("reads: " + std::to_string(hs.reads) + " -> " + std::to_string(hs.fragments) +
                 " fragments (" + std::to_string(hs.paired) + " paired, " +
                 std::to_string(hs.singleton) + " single)");
        if (hs.over_paired > 0) {
            log.info("WARNING: " + std::to_string(hs.over_paired) + " reads shared a name with an "
                     "already-complete pair and were dropped; fragment counts are unreliable");
        }
        if (hs.paired == 0 && hopt.use_insert_size) {
            hopt.use_insert_size = false;
            log.info("no paired fragments: insert-size channel off");
        }
        std::vector<std::string> names;
        for (const PathRecord& p : panel_graph.paths) names.push_back(p.name);
        log.info("scoring " + std::to_string(names.size()) + " panel haplotypes (shortlist " +
                 std::to_string(hopt.max_haplotypes) + ") over " + std::to_string(frags.size()) +
                 " fragments");
        const HaplotypeResult hr =
            genotype_haplotype_pairs(chain, blocks, names, frags, hopt,
                                     have_truth ? &truth1 : nullptr,
                                     have_truth ? &truth2 : nullptr, top_pairs);
        if (!hr.top_pairs.empty()) {
            log.info("best pair: " + hr.shortlist[hr.top_pairs[0].hap1] + " / " +
                     hr.shortlist[hr.top_pairs[0].hap2] + " (posterior " +
                     std::to_string(hr.top_pairs[0].posterior) + ")");
        }
        log.info(std::to_string(hr.n_informative) + " of " + std::to_string(hr.n_fragments) +
                 " fragments discriminate between shortlisted haplotypes");
        if (have_truth) {
            // Whether the coarse stage kept the answer is reported, never assumed: a candidate
            // generator that drops the truth loses it outright and no score below can recover it.
            const std::vector<std::string> tn = split_commas(truth_haplotypes);
            for (const std::string& n : tn) {
                const bool kept = std::find(hr.shortlist.begin(), hr.shortlist.end(), n) != hr.shortlist.end();
                const bool in_panel =
                    std::any_of(panel_graph.paths.begin(), panel_graph.paths.end(),
                                [&](const PathRecord& p) { return p.name == n; });
                if (in_panel && !kept) {
                    log.info("WARNING: truth haplotype " + n + " is in the panel but did NOT survive "
                             "the coarse shortlist; raise --max-haplotypes");
                }
            }
            std::size_t rep = 0, ex = 0, bub_rep = 0, bub_ex = 0;
            for (const BlockProjection& p : hr.blocks) {
                if (!p.truth_representable) continue;
                ++rep; if (p.exact) ++ex;
                if (p.kind == BlockKind::Bubble) { ++bub_rep; if (p.exact) ++bub_ex; }
            }
            log.info("blocks with a representable truth: " + std::to_string(rep) + "; exact " +
                     std::to_string(ex) + " (" + std::to_string(rep ? 100 * ex / rep : 0) +
                     "%); bubble blocks " + std::to_string(bub_ex) + "/" + std::to_string(bub_rep));
        }
        write_haplotype_results(out_prefix, hr, have_truth);
        log.wrote({out_prefix + ".hap_blocks.tsv", out_prefix + ".hap_pairs.tsv",
                   out_prefix + ".hap_scores.tsv"});
        log.done();
        return 0;
    }

    // ---- which blocks ------------------------------------------------------------------------
    std::vector<std::size_t> targets;
    if (!blocks_arg.empty()) {
        for (const std::string& tok : split_commas(blocks_arg)) {
            const std::size_t bi = cli::parse_size_arg("--blocks", tok);
            if (bi >= chain.size()) {
                throw std::runtime_error("genotype-frag: --blocks " + tok + " is past the chain (" +
                                         std::to_string(chain.size()) + " blocks)");
            }
            targets.push_back(bi);
        }
    } else {
        for (std::size_t bi = 0; bi < chain.size(); ++bi) {
            if (all_blocks || chain[bi].kind == BlockKind::Bubble) targets.push_back(bi);
        }
    }
    if (targets.empty()) throw std::runtime_error("genotype-frag: no blocks selected");

    FragmentLoadStats fstats;
    const std::vector<Fragment> fragments = load_fragments(read_paths, &fstats);
    log.info("reads: " + std::to_string(fstats.reads) + " -> " + std::to_string(fstats.fragments) +
             " fragments (" + std::to_string(fstats.paired) + " paired, " +
             std::to_string(fstats.singleton) + " single)");
    if (fstats.over_paired > 0) {
        log.info("WARNING: " + std::to_string(fstats.over_paired) + " reads shared a name with an "
                 "already-complete pair and were dropped. Fragment counts are unreliable; check that "
                 "the input is one library and that read names are not duplicated across files");
    }
    if (fstats.paired == 0 && opt.use_insert_size) {
        opt.use_insert_size = false;
        log.info("no paired fragments: the insert-size channel is off, sequence compatibility only");
    }

    log.info("scoring " + std::to_string(targets.size()) + " block(s) from " +
             std::to_string(fragments.size()) + " fragments");
    const std::vector<BlockFragmentResult> results =
        genotype_fragments(chain, blocks, fragments, targets, opt,
                           have_truth ? &truth1 : nullptr, have_truth ? &truth2 : nullptr,
                           top_pairs);

    std::size_t scored = 0, exact = 0, unrep = 0, pruned = 0, no_frag = 0;
    for (const BlockFragmentResult& r : results) {
        if (r.n_fragments == 0) ++no_frag;
        if (!have_truth) continue;
        if (r.truth_rank == -1) { ++unrep; continue; }
        if (r.truth_rank == -2) { ++pruned; continue; }
        ++scored;
        const std::size_t ta = std::min<std::size_t>(r.truth_a, r.truth_b);
        const std::size_t tb = std::max<std::size_t>(r.truth_a, r.truth_b);
        if (ta == std::min(r.best_a, r.best_b) && tb == std::max(r.best_a, r.best_b)) ++exact;
    }
    if (no_frag > 0) {
        log.info(std::to_string(no_frag) + " block(s) recruited no fragment at all -- their call is "
                 "arbitrary, not merely uncertain");
    }
    if (have_truth) {
        log.info("truth scored at " + std::to_string(scored) + " blocks: " + std::to_string(exact) +
                 " exact (" +
                 std::to_string(scored == 0 ? 0 : 100 * exact / scored) + "%); " +
                 std::to_string(unrep) + " unrepresentable, " + std::to_string(pruned) +
                 " pruned before scoring");
    }

    write_fragment_results(out_prefix, results, have_truth);
    log.wrote({out_prefix + ".frag_blocks.tsv", out_prefix + ".frag_pairs.tsv"});
    log.done();
    return 0;
}

} // namespace panvar
