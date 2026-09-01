#include "panvar/genotype_frag_command.hpp"

#include "panvar/align.hpp"
#include "panvar/bubble_path.hpp"
#include "panvar/bubbles.hpp"
#include "panvar/cli_utils.hpp"
#include "panvar/genotype_blocks.hpp"
#include "panvar/genotype_fragments.hpp"
#include "panvar/gfa.hpp"
#include "panvar/graph_utils.hpp"
#include "panvar/md5.hpp"
#include "panvar/output.hpp"
#include "panvar/parallel.hpp"
#include "panvar/syncmer.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <set>
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
        << "      --haplotype-mode        DIAGNOSTIC ORACLE, not the product. Scores whole panel\n"
        << "                              HAPLOTYPE pairs end to end and projects\n"
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
        << "      --project-marginal      Call each block from the allele pair with the most\n"
        << "                              posterior mass summed over haplotype pairs, instead of from\n"
        << "                              the best pair. It sounds stricter and measures worse: at\n"
        << "                              cyp2d6 under leave-one-out the best pair reconstructs the\n"
        << "                              donor to a median 130 edits and the marginal projection of\n"
        << "                              the SAME posterior to 385, because a per-block argmax of\n"
        << "                              marginals assembles a combination no single pair realises.\n"
        << "                              Kept so that stays reproducible\n"
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
        << "      --total-depth           Score the pair's TOTAL length against how many fragments\n"
        << "                              were seen: N ~ Poisson(lambda * (L(h1)+L(h2))), with lambda\n"
        << "                              fitted once from the observed count over twice the panel's\n"
        << "                              median haplotype length -- outside any candidate, so the\n"
        << "                              term cannot be self-fulfilling. Each fragment is counted\n"
        << "                              once by construction, because the observation is a single\n"
        << "                              scalar; that is the defect the per-haplotype window channel\n"
        << "                              could not escape, where one fragment made every candidate\n"
        << "                              copy it aligned to look covered\n"
        << "      --haploid-depth <x>     Supply lambda (fragments per bp per haplotype copy) instead\n"
        << "                              of estimating it; implies --total-depth\n"
        << "      --truth-total-bp <N>    DIAGNOSTIC: rank pairs by how close their total length is to\n"
        << "                              N, the sample's true diploid total, which no caller can\n"
        << "                              know. Bounds what perfect dosage knowledge could buy -- if\n"
        << "                              the tail does not close under this, it is not a dosage\n"
        << "                              problem and no depth model will fix it\n"
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
        << "      --exact-distance <a.fa> <b.fa>\n"
        << "                              Global (Needleman-Wunsch) edit distance between two\n"
        << "                              sequences, printed and nothing else. The panel floor is\n"
        << "                              computed from minimap2, which reports alignments rather\n"
        << "                              than a distance: split or overlapping records can double\n"
        << "                              count edits, and anything it declines to align is charged\n"
        << "                              as a whole unaligned base. Use this to check the floor\n"
        << "                              exactly on the few nearest candidates -- it is affordable\n"
        << "                              there and not over a whole panel\n"
        << "      --reference-score <a.fa> <b.fa>\n"
        << "                              Score this haplotype pair against --reads with the EXACT\n"
        << "                              reference implementation of the model contract: every\n"
        << "                              fragment start on both haplotypes enumerated, the insert\n"
        << "                              prior summed over, no syncmer index, no anchor cap and no\n"
        << "                              --placement-topk. It is the oracle the fast path is tested\n"
        << "                              against, not a caller -- it is O(fragments x length x insert\n"
        << "                              range) and is meant for small synthetic haplotypes. Pass the\n"
        << "                              same file twice for a homozygous pair; exposure and the\n"
        << "                              placement set both double, as the contract requires\n"
        << "      --mosaic-floor          Three floors instead of one, needing --truth-haplotypes\n"
        << "                              and no reads: COMPLETE (one panel haplotype per homologue\n"
        << "                              across the locus), FREE (the nearest panel allele at every\n"
        << "                              block, chosen independently -- a bound, since it may switch\n"
        << "                              source at every boundary with no evidence for any switch),\n"
        << "                              and PENALISED (the same with a cost per switch). Their\n"
        << "                              difference is the most a mosaic model could ever recover,\n"
        << "                              which is the measurement that decides whether a factor\n"
        << "                              graph is worth building at all\n"
        << "      --switch-penalties <a,b,c>  Penalties to evaluate (default 0,10,100,1000)\n"
        << "      --spell-calls <gt.tsv>  Read a per-block call table (either caller's) and write the\n"
        << "                              two sequences it claims, to <out-prefix>.called.fa. Needed\n"
        << "                              because the two callers report in different shapes and\n"
        << "                              cannot otherwise be compared on the objective that matters:\n"
        << "                              production emits a per-block allele pair, which is a mosaic\n"
        << "                              and has no called haplotype to align. Needs no reads\n"
        << "      --all-blocks            Score backbone and flank blocks too. gstm1's worst blocks\n"
        << "                              are BACKBONE blocks carrying 16 and 21 alleles, so the\n"
        << "                              default view omits them\n"
        << "      --spell-pair <hap_pairs.tsv>  Spell the rank-1 pair BY PATH NAME, straight from the\n"
        << "                              GFA walks. Prefer this to --spell-calls in haplotype mode:\n"
        << "                              path names survive a change of graph, decomposition or\n"
        << "                              exclusion, and integer allele indices do not\n"
        << "      --distance-band <N>     Cap --exact-distance banding at N; a pair outside it prints\n"
        << "                              \">N\" instead of falling back to the unbanded computation.\n"
        << "                              Certifying a panel floor only needs to know what is NEARER\n"
        << "                              than the best so far\n"
        << "      --dump-scored-sequences <prefix>  Write the exact sequences whole-haplotype mode\n"
        << "                              scores, AFTER --exclude-haplotypes, as <prefix>.scored_sequences\n"
        << "                              .fa plus a .tsv of name, length and md5 beside the raw GFA path\n"
        << "                              spelling of the same name. Runs standalone with no --reads. Use\n"
        << "                              it to check the graph against the assembly it was built from:\n"
        << "                              the round-trip invariant only checks blocks against the GFA.\n"
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

namespace {
// THE target list, built once. The oracle previously built its own and ignored --blocks, so it could
// score a different chain from the one production and the incidence table use -- and then "adjacent
// target" means two different things in two places.
std::vector<std::size_t> build_targets(const std::vector<Block>& chain,
                                       const std::string& blocks_arg,
                                       bool all_blocks) {
    std::vector<std::size_t> targets;
    if (!blocks_arg.empty()) {
        for (const std::string& tok : split_commas(blocks_arg)) {
            const std::size_t bi = cli::parse_size_arg("--blocks", tok);
            if (bi >= chain.size()) {
                throw std::runtime_error("genotype-frag: --blocks " + tok + " is past the chain (" +
                                         std::to_string(chain.size()) + " blocks)");
            }
            // CHAIN ORDER IS PART OF THE CONTRACT, not a convention. classify_fragment_factors and
            // chain_span_sequence both index by RANK and derive adjacency from it, so an unsorted or
            // duplicated list silently redefines what "adjacent target" means. Rejected rather than
            // sorted: reordering what was asked for would hide the mistake instead of naming it.
            if (!targets.empty() && bi <= targets.back()) {
                throw std::runtime_error(
                    "genotype-frag: --blocks must be strictly increasing in chain order; got " +
                    std::to_string(bi) + " after " + std::to_string(targets.back()) +
                    ". Adjacency between targets is derived from this order, so an unsorted or "
                    "duplicated list would change which fragments count as boundary-spanning.");
            }
            targets.push_back(bi);
        }
    } else {
        for (std::size_t bi = 0; bi < chain.size(); ++bi) {
            if (all_blocks || chain[bi].kind == BlockKind::Bubble) targets.push_back(bi);
        }
    }
    if (targets.empty()) throw std::runtime_error("genotype-frag: no blocks selected");
    return targets;
}
} // namespace

int run_genotype_frag_command(const std::vector<std::string>& args) {
    if (args.empty()) { print_help(); return 0; }

    std::string gfa_path, bubble_prefix_in, bubbles_csv_in, out_prefix;
    std::vector<std::string> read_paths;
    std::string truth_haplotypes, exclude_haplotypes, blocks_arg, spell_calls, dump_sequences;
    bool mosaic_floor = false;
    std::vector<std::string> reference_pair;
    std::vector<std::string> ref_block, ref_block_pair;
    // Escape hatch for the equality gate, which must compare the factor oracle against the
    // whole-haplotype reference over the SAME fragment set. Not for measurement.
    std::string ref_subset;
    std::string origin_universe;
    double scope_tol = 1e-6;
    std::vector<std::string> reconcile_scope;
    std::string switch_penalties_arg = "0,10,100,1000";
    std::vector<std::string> exact_distance;
    std::size_t distance_band = 0;
    std::string spell_pair;
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
        else if (a == "--spell-calls") spell_calls = value(i, a);
        else if (a == "--spell-pair") spell_pair = value(i, a);
        else if (a == "--distance-band") distance_band = cli::parse_size_arg(a, value(i, a));
        else if (a == "--dump-scored-sequences") dump_sequences = value(i, a);
        else if (a == "--mosaic-floor") mosaic_floor = true;
        else if (a == "--reference-score") { reference_pair.push_back(value(i, a));
                                             reference_pair.push_back(value(i, a)); }
        else if (a == "--reference-all-fragments") ref_subset = "all";
        else if (a == "--origin-universe") origin_universe = value(i, a);
        else if (a == "--scope-tol") {
            const std::string sv = value(i, a);
            scope_tol = std::stod(sv);
            // Zero is meaningful -- it asks for the exact structural scope, every block whose
            // removal moves the likelihood at all. Negative or non-finite is not: a negative
            // tolerance puts every touched block in scope and then reports a bound that "holds"
            // against a threshold no residual can meet, and a NaN makes every comparison false,
            // so the joint bound silently passes for any scope whatsoever.
            if (!std::isfinite(scope_tol) || scope_tol < 0.0) {
                throw std::runtime_error(
                    "genotype-frag: --scope-tol must be finite and >= 0; got " + sv);
            }
        }
        else if (a == "--reconcile-scope") { reconcile_scope.push_back(value(i, a));
                                             reconcile_scope.push_back(value(i, a)); }
        else if (a == "--reference-block") {
            // <target_index> <allele_a> <allele_b> -- the exact LOCAL oracle for one target.
            // target_index is a RANK into the scored-block list (bubbles only unless --all-blocks),
            // not a raw block index: targets 0 and 1 are typically raw blocks 1 and 3.
            for (int q = 0; q < 3; ++q) ref_block.push_back(value(i, a));
        }
        else if (a == "--reference-block-pair") {
            // <target_index> <a1> <b1> <a2> <b2> -- the exact oracle over two ADJACENT TARGETS. Its
            // difference from the two unary scores is what a transition factor is worth.
            for (int q = 0; q < 5; ++q) ref_block_pair.push_back(value(i, a));
        }
        else if (a == "--switch-penalties") switch_penalties_arg = value(i, a);
        else if (a == "--exact-distance") { exact_distance.push_back(value(i, a));
                                            exact_distance.push_back(value(i, a)); }
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
        else if (a == "--project-marginal") hopt.project_map = false;
        else if (a == "--marginalise-placements") hopt.marginalise_placements = true;
        else if (a == "--placement-topk") hopt.placement_topk = cli::parse_size_arg(a, value(i, a));
        else if (a == "--joint-depth") hopt.joint_depth = true;
        else if (a == "--joint-marginal") { hopt.joint_depth = true; hopt.joint_marginal = true; }
        else if (a == "--rung-zero") {
            // The accelerated event enumeration with NONE of the approximations, so that a
            // difference from the exact reference can only be the enumeration itself.
            hopt.joint_depth = true; hopt.joint_marginal = true; hopt.rung_zero = true;
            hopt.max_anchor_occ = std::numeric_limits<std::size_t>::max();
            hopt.placement_topk = std::numeric_limits<std::size_t>::max();
            hopt.placement_bin = 1;
            hopt.placement_dedup = 0;      // keep every distinct (start, end) state
            hopt.joint_top_pairs = 0;
        }
        else if (a == "--placement-bin") hopt.placement_bin = cli::parse_size_arg(a, value(i, a));
        else if (a == "--hamming-emission") hopt.hamming_emission = true;
        else if (a == "--band-floor") hopt.band_floor = true;
        else if (a == "--incidence") opt.incidence_path = value(i, a);
        else if (a == "--unplaced-neutral") hopt.unplaced_neutral = true;
        else if (a == "--force-haplotypes") {
            for (const std::string& n : split_commas(value(i, a))) hopt.force_haplotypes.push_back(n);
        }
        else if (a == "--dump-fragment-mass") hopt.dump_fragment_mass = value(i, a);
        else if (a == "--dump-mass-pair") { const std::vector<std::string> two = split_commas(value(i, a));
            if (two.size() != 2) throw std::runtime_error("genotype-frag: --dump-mass-pair needs two names");
            hopt.dump_mass_pair1 = two[0]; hopt.dump_mass_pair2 = two[1]; }
        else if (a == "--placement-dedup") hopt.placement_dedup = cli::parse_size_arg(a, value(i, a));
        else if (a == "--mate-rescue") hopt.mate_rescue = true;
        else if (a == "--coordinate-join") hopt.coordinate_join = true;
        else if (a == "--no-coordinate-join") { hopt.coordinate_join = false; hopt.force_join = false; }
        else if (a == "--force-join") { hopt.coordinate_join = true; hopt.force_join = true; }
        else if (a == "--zero-seed-fallback") hopt.zero_seed_fallback = true;
        else if (a == "--zero-seed-complete") { hopt.zero_seed_fallback = true;
                                                hopt.zero_seed_complete = true; }
        else if (a == "--zero-seed-exhaustive") { hopt.zero_seed_fallback = true;
                                                  hopt.zero_seed_exhaustive = true; }
        else if (a == "--multiplicity-aware") {
            hopt.multiplicity_aware = true;
            hopt.max_anchor_occ = std::numeric_limits<std::size_t>::max();
            hopt.placement_topk = std::numeric_limits<std::size_t>::max();
        }
        else if (a == "--mass-tolerance") hopt.mass_tolerance = std::stod(value(i, a));
        else if (a == "--joint-reverse-order") hopt.joint_reverse_order = true;
        else if (a == "--equivalence-tolerance") hopt.equivalence_tolerance = std::stod(value(i, a));
        else if (a == "--joint-top-pairs") hopt.joint_top_pairs = cli::parse_size_arg(a, value(i, a));
        else if (a == "--joint-window") hopt.joint_window = cli::parse_size_arg(a, value(i, a));
        else if (a == "--joint-iterations") hopt.joint_iterations = cli::parse_size_arg(a, value(i, a));
        else if (a == "--total-depth") hopt.total_depth = true;
        else if (a == "--haploid-depth") { hopt.haploid_depth = std::stod(value(i, a)); hopt.total_depth = true; }
        else if (a == "--truth-total-bp") hopt.truth_total_bp = std::stod(value(i, a));
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

    // Flag-consistency checks run HERE, before any early-exit mode (--reference-score,
    // --exact-distance, --mosaic-floor ...) returns. Placed after them, as they first were, an
    // inconsistent combination was accepted silently in exactly the modes used for oracle
    // comparison -- where a silent mismatch does the most damage.
    if (hopt.coordinate_join && !hopt.use_insert_size) {
        // The join enumerates insert lengths over the prior's support. Without the insert term every
        // FR pair with the reverse mate downstream contributes, the support is the whole haplotype,
        // and the join would silently become a restriction rather than a re-ordering of the same sum.
        throw std::runtime_error(
            "genotype-frag: --coordinate-join requires the insert-size term. The join sweeps insert "
            "lengths over the prior's bounded support, which is exact only because combinations "
            "outside it carry zero mass; with the insert term off they do not");
    }
    if (hopt.zero_seed_fallback && !hopt.hamming_emission) {
        throw std::runtime_error(
            "genotype-frag: --zero-seed-fallback requires --hamming-emission. The fallback verifies "
            "candidate starts at FIXED positions by Hamming distance; under the banded "
            "local-alignment emission a read can slide, and the placements it returns would not be "
            "the ones the scorer then evaluates");
    }
    if (hopt.mate_rescue && !hopt.hamming_emission) {
        // The rescue matches the unanchored mate at fixed positions by Hamming distance, because the
        // interval it searches is defined by coordinates. The default emission is banded local
        // alignment, which can slide. Allowing both at once would silently mix two emission models
        // within a single fragment -- one mate scored by alignment, its partner by Hamming.
        throw std::runtime_error(
            "genotype-frag: --mate-rescue currently requires --hamming-emission. The rescued mate is "
            "scored at fixed positions by Hamming distance, so combining it with the banded "
            "local-alignment emission would mix two emission models inside one fragment. A "
            "fixed-start edit likelihood for the normal emission would lift this restriction");
    }
    if (hopt.multiplicity_aware && !hopt.joint_marginal) {
        // Grouping by likelihood is exact only where POSITION is irrelevant, which is the global
        // marginal scorer. In the windowed or hard-assignment models, equal-likelihood placements in
        // different windows cannot share one representative coordinate, and collapsing them would
        // move coverage between windows.
        throw std::runtime_error(
            "genotype-frag: --multiplicity-aware requires --joint-marginal. Grouping placements by "
            "likelihood is exact only for the global marginal model, where a placement's position "
            "does not enter the score; under the windowed model equal-likelihood placements in "
            "different windows are not interchangeable");
    }

    if (!reference_pair.empty()) {
        if (read_paths.empty()) throw std::runtime_error("genotype-frag: --reference-score needs --reads");
        const auto slurp_fa = [](const std::string& path) {
            std::ifstream in(path);
            if (!in) throw std::runtime_error("genotype-frag: cannot read " + path);
            std::string line, seq;
            while (std::getline(in, line)) {
                if (!line.empty() && line[0] == '>') continue;
                while (!line.empty() && (line.back() == '\r' || line.back() == '\n')) line.pop_back();
                seq += line;
            }
            return seq;
        };
        FragmentLoadStats st;
        const std::vector<Fragment> frags = load_fragments(read_paths, &st);
        ReferenceParams rp;
        rp.lambda = hopt.haploid_depth > 0.0 ? hopt.haploid_depth : 0.05;
        rp.eta = hopt.outlier_mix;
        rp.error_rate = hopt.error_rate;
        rp.fragment_len = hopt.fragment_len;
        rp.fragment_sd = hopt.fragment_sd;
        rp.bg_divergence = hopt.bg_divergence;
        std::vector<double> mass, contrib;
        const double v = reference_pair_loglik(slurp_fa(reference_pair[0]), slurp_fa(reference_pair[1]),
                                               frags, rp,
                                               hopt.dump_fragment_mass.empty() ? nullptr : &mass,
                                               hopt.dump_fragment_mass.empty() ? nullptr : &contrib);
        if (!hopt.dump_fragment_mass.empty()) {
            std::ofstream mf(hopt.dump_fragment_mass);
            if (!mf) throw std::runtime_error("genotype-frag: cannot write " + hopt.dump_fragment_mass);
            mf.precision(17);
            mf << "# reference\n" << "fragment\tlog_mass\tmates_seeded\tcontrib\n";
            for (std::size_t i = 0; i < mass.size() && i < frags.size(); ++i) {
                mf << frags[i].name << '\t' << mass[i] << "\tNA\t" << contrib[i] << '\n';
            }
        }
        std::printf("%.17g\n", v);
        return 0;
    }
    if (!exact_distance.empty()) {
        const auto slurp = [](const std::string& path) {
            std::ifstream in(path);
            if (!in) throw std::runtime_error("genotype-frag: cannot read " + path);
            std::string line, seq;
            while (std::getline(in, line)) {
                if (!line.empty() && line[0] == '>') continue;
                while (!line.empty() && (line.back() == '\r' || line.back() == '\n')) line.pop_back();
                seq += line;
            }
            return seq;
        };
        const std::string a = slurp(exact_distance[0]);
        const std::string b = slurp(exact_distance[1]);
        if (a.empty() || b.empty()) { std::cout << (a.size() + b.size()) << '\n'; return 0; }
        // Progressive banding. An unbanded global alignment of two 200 kb haplotypes is quadratic in
        // practice and does not finish in useful time, but the distances this is used on are small --
        // a panel floor is tens to a few thousand edits over 200 kb. edlib returns the EXACT distance
        // whenever it is within the band, so doubling from a small band is exact and fast in the
        // common case and only falls back to the full computation when the sequences really are far
        // apart. The band actually used is reported on stderr so a number can never be mistaken for
        // an unbanded one.
        // --distance-band caps the escalation. Certifying a panel floor asks 464 x 2 questions of
        // the form "is anything nearer than the best so far", and for all but a handful the answer
        // is a large number nobody needs -- computing it unbanded is most of the cost and none of
        // the information. With the cap, a pair outside the band reports ">N" and the floor is
        // certified as long as the winning distance is <= N.
        std::vector<std::size_t> bands{1024, 4096, 16384, 65536, 262144};
        if (distance_band != 0) {
            bands.clear();
            for (std::size_t b2 = 1024; b2 < distance_band; b2 *= 4) bands.push_back(b2);
            bands.push_back(distance_band);
        }
        for (const std::size_t band : bands) {
            const NwBanded r = nw_edit_distance_banded(a, b, band);
            if (r.ok) {
                std::fprintf(stderr, "band %zu sufficed\n", band);
                std::cout << r.edits << '\n';
                return 0;
            }
        }
        if (distance_band != 0) {
            std::fprintf(stderr, "exceeds band %zu\n", distance_band);
            std::cout << '>' << distance_band << '\n';
            return 0;
        }
        std::fprintf(stderr, "no band sufficed; computing unbanded\n");
        std::cout << nw_edit_distance(a, b).edits << '\n';
        return 0;
    }

    if (gfa_path.empty() || out_prefix.empty()) {
        throw std::runtime_error("genotype-frag requires --gfa and --out-prefix");
    }
    if (read_paths.empty() && spell_calls.empty() && !mosaic_floor && dump_sequences.empty() &&
        spell_pair.empty() && ref_block.empty() && ref_block_pair.empty() &&
        origin_universe.empty() && reconcile_scope.empty()) {
        throw std::runtime_error("genotype-frag requires at least one --reads");
    }
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

    // The held-out haplotypes get their own decomposition: they are gone from `blocks`, so their
    // own walks can only be spelled from a chain enumerated over them alone. Hoisted out of the
    // truth branch because it depends on the exclusion, not on whether truth names were given --
    // the sequence dump needs it with no --truth-haplotypes at all.
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

    // ---- the sequence dump ---------------------------------------------------------------------
    // What whole-haplotype mode actually scores, byte for byte, after the exclusion has been
    // applied -- beside the raw GFA spelling of the same path name.
    //
    // This exists because the LPA accuracy numbers were charging a representation difference to the
    // genotyper. The caller scores sequences spelled from the graph it is given; the truth FASTAs
    // were spelled from a DIFFERENT stage of the same pipeline, and every distance reported was the
    // sum of a genotyping error and a drift no one had measured. The existing round-trip invariant
    // could not see it: it checks the block spelling against the GFA path, never the GFA path
    // against the assembly the graph was built from. Both columns are emitted here so the two
    // comparisons are separable, and the md5 is the one that travels -- it can be checked against
    // an external FASTA with md5sum and nothing else.
    if (!dump_sequences.empty()) {
        cli::ensure_parent_dir_for_file(dump_sequences);
        const std::string tsv = dump_sequences + ".scored_sequences.tsv";
        const std::string fa = dump_sequences + ".scored_sequences.fa";
        std::ofstream tf(tsv), ff(fa);
        if (!tf) throw std::runtime_error("genotype-frag: cannot write " + tsv);
        if (!ff) throw std::runtime_error("genotype-frag: cannot write " + fa);
        tf << "group\tname\tscored_bp\tscored_md5\tblock_bp\tblock_md5\tgfa_bp\tgfa_md5"
              "\tgfa_complete\tround_trips\tframe\n";

        const auto by_name = path_records_by_name(graph);
        std::size_t n_mismatch = 0, n_incomplete = 0, n_rows = 0;

        const auto emit = [&](const char* group, const std::vector<BlockAlleles>& src,
                              const std::string& name) {
            // THE SEQUENCE ACTUALLY SCORED. Whole-haplotype mode scores the graph WALK, so dumping
            // the block concatenation here would recreate exactly the diagnostic mismatch this
            // change removed: a dump that disagrees with the thing it claims to describe. The block
            // spelling is still emitted, in its own column, because the two disagreeing is the
            // finding the dump exists to surface.
            std::string scored;
            {
                const auto wi = by_name.find(name);
                bool wok = false;
                if (wi != by_name.end() && wi->second != nullptr) {
                    scored = spell_path_steps_sequence(graph, wi->second->steps, &wok);
                }
                if (!wok) scored.clear();
            }
            const std::string block_spelling = spell_block_haplotype(src, name);
            if (scored.empty()) {
                // No fallback. A dump that quietly substituted the block spelling here would report
                // a sequence the scorer never sees -- the exact defect this dump exists to detect,
                // reintroduced inside the detector.
                throw std::runtime_error(
                    "genotype-frag: the graph cannot completely spell path '" + name +
                    "', so there is no scored sequence to dump. Refusing rather than substituting "
                    "the block concatenation, which is not what would be scored.");
            }
            std::string raw;
            bool complete = false;
            const auto it = by_name.find(name);
            if (it != by_name.end() && it->second != nullptr) {
                raw = spell_path_steps_sequence(graph, it->second->steps, &complete);
            }
            // A path the graph cannot fully spell has no raw sequence to compare against, so it is
            // reported as such rather than as a round-trip failure -- those are different faults.
            // FRAME: the block chain is reference-oriented, so a path antiparallel to the reference
            // spells the reverse complement of its walk. Same sequence, different frame -- reported
            // as a round trip, with the frame named, rather than refused.
            // round_trips / frame describe the BLOCK spelling against the walk -- that is the
            // decomposition question. `scored` is the walk and is authoritative regardless.
            const std::string rc = complete ? reverse_complement(raw) : std::string();
            const bool same_frame = complete && block_spelling == raw;
            const bool rc_frame = complete && !same_frame && block_spelling == rc;
            const bool round_trips = same_frame || rc_frame;
            const char* frame = !complete ? "." : (same_frame ? "fwd" : (rc_frame ? "rc" : "MISMATCH"));
            if (!complete) ++n_incomplete;
            else if (!round_trips) ++n_mismatch;
            ++n_rows;
            tf << group << '\t' << name << '\t' << scored.size() << '\t' << md5_hex(scored)
               << '\t' << block_spelling.size() << '\t' << md5_hex(block_spelling)
               << '\t' << (complete ? raw.size() : std::size_t{0}) << '\t'
               << (complete ? md5_hex(raw) : std::string(".")) << '\t' << (complete ? "yes" : "no")
               << '\t' << (complete ? (round_trips ? "yes" : "NO") : ".") << '\t' << frame << '\n';
            ff << '>' << name << ' ' << group << '\n';
            for (std::size_t off = 0; off < scored.size(); off += 60) {
                ff << scored.substr(off, 60) << '\n';
            }
            if (scored.empty()) ff << '\n';
        };

        for (const PathRecord& p : panel_graph.paths) emit("panel", blocks, p.name);
        for (const PathRecord& p : held_out) emit("held_out", held_blocks, p.name);

        tf.flush();
        ff.flush();
        if (!tf) throw std::runtime_error("genotype-frag: write failed for " + tsv);
        if (!ff) throw std::runtime_error("genotype-frag: write failed for " + fa);
        log.info("dumped " + std::to_string(n_rows) + " scored sequences (" +
                 std::to_string(panel_graph.paths.size()) + " panel, " +
                 std::to_string(held_out.size()) + " held out); " +
                 std::to_string(n_mismatch) + " do not round-trip against the GFA path, " +
                 std::to_string(n_incomplete) + " have no complete GFA spelling");
        log.wrote({tsv, fa});
        // Standalone when nothing else was asked for: the audit case is a graph and a panel and
        // nothing else. Combined with another mode the dump is a side effect and the run continues,
        // so --spell-calls and --mosaic-floor still do their own work.
        if (read_paths.empty() && spell_calls.empty() && !mosaic_floor) {
            log.done();
            return 0;
        }
    }

    // ---- spell a called PAIR by path name ------------------------------------------------------
    // The safe route for haplotype mode. Path names survive a change of graph, decomposition or
    // exclusion; integer allele indices do not, which is the whole reason the manifest below exists.
    // Spelling straight from the GFA walks also bypasses the block decomposition entirely, so this
    // cannot inherit a decomposition fault.
    if (!spell_pair.empty()) {
        std::ifstream pf(spell_pair);
        if (!pf) throw std::runtime_error("genotype-frag: cannot read " + spell_pair);
        std::string line;
        if (!std::getline(pf, line)) throw std::runtime_error("genotype-frag: empty " + spell_pair);
        std::vector<std::string> header;
        {
            std::size_t start = 0;
            for (std::size_t i = 0; i <= line.size(); ++i) {
                if (i == line.size() || line[i] == '\t') {
                    header.push_back(line.substr(start, i - start));
                    start = i + 1;
                }
            }
        }
        const auto col = [&](const std::string& n) {
            for (std::size_t i = 0; i < header.size(); ++i) if (header[i] == n) return static_cast<long>(i);
            throw std::runtime_error("genotype-frag: --spell-pair table has no column '" + n + "'");
        };
        const long c1 = col("hap1"), c2 = col("hap2");
        if (!std::getline(pf, line)) {
            throw std::runtime_error("genotype-frag: " + spell_pair + " has a header but no rows");
        }
        std::vector<std::string> f;
        {
            std::size_t start = 0;
            for (std::size_t i = 0; i <= line.size(); ++i) {
                if (i == line.size() || line[i] == '\t') { f.push_back(line.substr(start, i - start)); start = i + 1; }
            }
        }
        if (f.size() <= static_cast<std::size_t>(std::max(c1, c2))) {
            throw std::runtime_error("genotype-frag: " + spell_pair + " rank-1 row is truncated");
        }
        const std::string n1 = f[static_cast<std::size_t>(c1)], n2 = f[static_cast<std::size_t>(c2)];
        const auto by_name = path_records_by_name(graph);
        const auto spell_one = [&](const std::string& nm) {
            const auto it = by_name.find(nm);
            if (it == by_name.end() || it->second == nullptr) {
                throw std::runtime_error("genotype-frag: --spell-pair names path '" + nm +
                                         "' which is not in " + gfa_path);
            }
            bool complete = false;
            std::string out = spell_path_steps_sequence(graph, it->second->steps, &complete);
            if (!complete) {
                throw std::runtime_error("genotype-frag: the graph cannot fully spell path '" + nm + "'");
            }
            // Emit in the same FRAME as the block spelling. The panel, the truth and every distance
            // in the benchmarks are reference-oriented block spellings; a walk-oriented sequence for
            // an antiparallel path would compare as a whole-length mismatch. The flip is DECIDED by
            // comparison, never assumed, so a path whose two spellings genuinely disagree still
            // surfaces rather than being silently reverse-complemented into looking fine.
            if (block_spelling_frame(graph, blocks, nm) == SpellFrame::ReverseComplement) {
                out = reverse_complement(out);
            }
            return out;
        };
        const std::string s1 = spell_one(n1), s2 = spell_one(n2);
        const std::string fa = out_prefix + ".called.fa";
        std::ofstream of(fa);
        if (!of) throw std::runtime_error("genotype-frag: cannot write " + fa);
        of << '>' << n1 << "\n" << s1 << "\n>" << n2 << "\n" << s2 << '\n';
        of.flush();
        if (!of) throw std::runtime_error("genotype-frag: write failed for " + fa);
        log.info("spelled the rank-1 pair by NAME: " + n1 + " (" + std::to_string(s1.size()) +
                 " bp, md5 " + md5_hex(s1) + ") and " + n2 + " (" + std::to_string(s2.size()) +
                 " bp, md5 " + md5_hex(s2) + ")");
        log.wrote({fa});
        log.done();
        return 0;
    }

    if (!reconcile_scope.empty()) {
        // THE GATE. Scope comes from the origin universe over the WHOLE candidate set, computed
        // before the pair is looked at, so the factor topology cannot depend on what is being
        // scored. Then the same pair is scored twice: once by the whole-locus reference, once by
        // summing only in-scope origins with exposure charged once.
        const std::vector<Fragment> rfr = load_fragments(read_paths);
        ReferenceParams rp;
        rp.lambda = hopt.haploid_depth > 0.0 ? hopt.haploid_depth : 0.05;
        rp.eta = opt.outlier_mix; rp.error_rate = opt.error_rate;
        rp.fragment_len = opt.fragment_len; rp.fragment_sd = opt.fragment_sd;
        rp.bg_divergence = opt.bg_divergence;

        std::vector<CandidateFrame> frames;
        std::vector<std::string> cname;
        {
            const auto by_name = path_records_by_name(graph);
            for (const PathRecord& pr : panel_graph.paths) {
                const auto it = by_name.find(pr.name);
                if (it == by_name.end() || it->second == nullptr) continue;
                bool complete = false;
                const std::string walk =
                    spell_path_steps_sequence(graph, it->second->steps, &complete);
                if (!complete) continue;
                const CandidateFrame cf = build_candidate_frame(blocks, pr.name, walk);
                if (!cf.ok) continue;
                frames.push_back(cf); cname.push_back(pr.name);
            }
        }
        if (frames.size() != panel_graph.paths.size()) {
            throw std::runtime_error(
                "genotype-frag: --reconcile-scope has a verified walk-to-block map for only " +
                std::to_string(frames.size()) + " of " +
                std::to_string(panel_graph.paths.size()) + " panel candidates. Reconciling over the "
                "reduced set would certify a factorisation on a smaller, easier state space than a "
                "caller would score. Exclude those candidates explicitly with --exclude-haplotypes, "
                "or fix the decomposition, so the oracle's candidate set and the scored set are the "
                "same by construction. Run --origin-universe to list them.");
        }
        std::vector<std::vector<std::uint32_t>> scopes(rfr.size());
        for (std::size_t fi = 0; fi < rfr.size(); ++fi) {
            scopes[fi] = enumerate_fragment_origins(rfr[fi], frames, rp,
                                                    hopt.placement_topk, scope_tol).scope;
        }
        long ia = -1, ib = -1;
        for (std::size_t i2 = 0; i2 < cname.size(); ++i2) {
            if (cname[i2] == reconcile_scope[0]) ia = static_cast<long>(i2);
            if (cname[i2] == reconcile_scope[1]) ib = static_cast<long>(i2);
        }
        if (ia < 0 || ib < 0) {
            throw std::runtime_error("genotype-frag: --reconcile-scope named a path not in the panel");
        }
        const double whole = reference_pair_loglik(frames[ia].seq, frames[ib].seq, rfr, rp);
        const double fact = scope_restricted_pair_loglik(frames[ia], frames[ib], rfr, scopes, rp);
        // BLOCK-LEVEL: the part a block-factored model could express. Unattributable mass is
        // dropped here on purpose. `whole - fact` being zero is guaranteed by always retaining that
        // mass and so proves nothing; `fact - block` is the component that no block state can carry,
        // and two candidates with the same block alleles must agree on `block` whatever it is.
        const double block =
            scope_restricted_pair_loglik(frames[ia], frames[ib], rfr, scopes, rp, false);
        std::printf("%.17g\t%.17g\t%.17g\t%.17g\t%.17g\n",
                    whole, fact, whole - fact, block, fact - block);
        log.done();
        return 0;
    }

    if (!origin_universe.empty()) {
        // The candidate-independent origin universe, per fragment, over EVERY panel candidate. It is
        // enumerated before any candidate is selected precisely so that the scope it yields cannot
        // depend on the candidate being scored.
        const std::vector<Fragment> ofr = load_fragments(read_paths);
        ReferenceParams rp;
        rp.lambda = hopt.haploid_depth > 0.0 ? hopt.haploid_depth : 0.05;
        rp.eta = opt.outlier_mix; rp.error_rate = opt.error_rate;
        rp.fragment_len = opt.fragment_len; rp.fragment_sd = opt.fragment_sd;
        rp.bg_divergence = opt.bg_divergence;

        // AUTHORITATIVE WALKS, with a verified coordinate map. Rebuilding candidates by
        // concatenating block alleles here would reintroduce the defect the walk rule removed:
        // wrong bytes for antiparallel paths and short bytes for unprojectable ones.
        std::vector<CandidateFrame> frames;
        std::vector<std::string> cand_names;
        std::vector<std::pair<std::string, std::string>> skipped_names;
        {
            const auto by_name = path_records_by_name(graph);
            std::size_t skipped = 0, flipped = 0;
            for (const PathRecord& pr : panel_graph.paths) {
                const auto it = by_name.find(pr.name);
                if (it == by_name.end() || it->second == nullptr) {
                    ++skipped; skipped_names.emplace_back(pr.name, "no path record"); continue;
                }
                bool complete = false;
                const std::string walk =
                    spell_path_steps_sequence(graph, it->second->steps, &complete);
                if (!complete) {
                    ++skipped; skipped_names.emplace_back(pr.name, "graph cannot spell the walk");
                    continue;
                }
                const CandidateFrame cf = build_candidate_frame(blocks, pr.name, walk);
                if (!cf.ok) {
                    ++skipped;
                    skipped_names.emplace_back(pr.name, "no verified walk-to-block map");
                    continue;
                }
                if (cf.reverse_frame) ++flipped;
                frames.push_back(cf);
                cand_names.push_back(pr.name);
            }
            if (!skipped_names.empty()) {
                const std::string sp = origin_universe + ".skipped";
                std::ofstream sf(sp);
                if (!sf) throw std::runtime_error("genotype-frag: cannot write " + sp);
                sf << "candidate\treason\n";
                for (const auto& [nm, why] : skipped_names) sf << nm << '\t' << why << '\n';
                sf.flush();
                log.wrote({sp});
            }
            // ALL-SKIPPED IS A FAILURE, not a result. With no frames every fragment reports zero
            // origins, empty scope, and a joint bound of 0 -- a table of plausible-looking zeros
            // that reads exactly like "the scope is trivially small" rather than "the instrument
            // never ran". A fixture whose segment names the sorter renumbers lands here, and the
            // zeros it produced passed three thresholds before the sidecar was read.
            if (frames.empty()) {
                throw std::runtime_error(
                    "genotype-frag: --origin-universe has no usable candidate frames; all " +
                    std::to_string(skipped_names.size()) + " panel candidates were skipped (see " +
                    origin_universe + ".skipped). Every value in the table would be a zero "
                    "produced by an instrument that never ran.");
            }
            log.info("candidate frames: " + std::to_string(frames.size()) + " usable (" +
                     std::to_string(flipped) + " antiparallel, mapped by mirroring), " +
                     std::to_string(skipped) +
                     " skipped for want of a verified walk-to-block map");
        }
        std::ofstream of(origin_universe);
        if (!of) throw std::runtime_error("genotype-frag: cannot write " + origin_universe);
        of.precision(17);
        of << "fragment\tn_origins\tn_candidates_hit\texact_lse\tretained_lse\tomitted_lse"
              "\tscope_blocks\tscope\tachieved_bound\tbound_holds"
              // appended, never inserted: three assertions once changed meaning silently when
              // columns were added in the middle of this header
              "\tinitial_bound\tblocks_added\tworst_pair\tunmapped_lse\n";
        for (const Fragment& F : ofr) {
            const OriginUniverse u =
                enumerate_fragment_origins(F, frames, rp, hopt.placement_topk, scope_tol);
            std::set<std::uint32_t> hits;
            for (const FragmentOrigin& o : u.origins) hits.insert(o.hap);
            of << F.name << '\t' << u.origins.size() << '\t' << hits.size() << '\t'
               << u.exact_lse << '\t' << u.retained_lse << '\t' << u.omitted_lse << '\t'
               << u.scope.size() << '\t';
            for (std::size_t q = 0; q < u.scope.size(); ++q) of << (q ? "," : "") << u.scope[q];
            if (u.scope.empty()) of << '.';
            of << '\t' << u.achieved_bound << '\t' << (u.bound_holds ? "yes" : "NO")
               << '\t' << u.initial_bound << '\t' << u.blocks_added
               << '\t' << u.worst_pair_a << ':' << u.worst_pair_b
               << '\t' << u.unmapped_lse << '\n';
        }
        of.flush();
        if (!of) throw std::runtime_error("genotype-frag: write failed for " + origin_universe);
        log.info("origin universe for " + std::to_string(ofr.size()) + " fragments over " +
                 std::to_string(frames.size()) + " candidates");
        log.wrote({origin_universe});
        log.done();
        return 0;
    }

    if (!ref_block.empty() || !ref_block_pair.empty()) {
        if (read_paths.empty()) {
            throw std::runtime_error("genotype-frag: --reference-block needs --reads");
        }
        std::vector<Fragment> rf_all = load_fragments(read_paths);

        // INCIDENCE-SELECTED EVIDENCE. A unary oracle handed the whole read set includes boundary and
        // unrelated fragments; a transition oracle handed the whole set scores everything. Their
        // difference then measures the sequences, not the linkage. So each factor sees only the
        // fragments the incidence rule assigns to it: `local` for a unary factor, `boundary` for an
        // adjacent-target factor. Ambiguous and unrecruited fragments belong to neither and are
        // excluded from both -- they need their own factor, not an arbitrary home.
        ReferenceParams rp;
        rp.lambda = hopt.haploid_depth > 0.0 ? hopt.haploid_depth : 0.05;
        rp.eta = opt.outlier_mix;
        rp.error_rate = opt.error_rate;
        rp.fragment_len = opt.fragment_len;
        rp.fragment_sd = opt.fragment_sd;
        rp.bg_divergence = opt.bg_divergence;

        // ONE classifier, shared with the production recruiter and the incidence dump. Rebuilding
        // the rule here is how a previous version selected zero boundary fragments where the real
        // rule selects 177.
        const std::vector<std::size_t> targets = build_targets(chain, blocks_arg, all_blocks);
        const FragmentFactorIncidence inc =
            classify_fragment_factors(blocks, targets, rf_all, opt.kmer_size, opt.syncmer_s,
                                      opt.flank_bp, opt.min_recruit_hits);

        if (!ref_block.empty()) {
            const std::size_t t = static_cast<std::size_t>(std::stoul(ref_block[0]));
            if (t >= targets.size()) {
                throw std::runtime_error("genotype-frag: --reference-block takes a TARGET INDEX (a "
                    "rank into the scored-block list), not a raw block index; there are " +
                    std::to_string(targets.size()) + " targets");
            }
            const int aa = std::stoi(ref_block[1]), ab = std::stoi(ref_block[2]);
            const std::string ca = chain_span_sequence(blocks, targets, t, t, {aa}, opt.flank_bp);
            const std::string cb = chain_span_sequence(blocks, targets, t, t, {ab}, opt.flank_bp);
            std::vector<Fragment> sel;
            for (std::size_t fi = 0; fi < rf_all.size(); ++fi) {
                if (inc.factor_of[fi] == FragmentFactor::Local &&
                    inc.targets_of[fi].front() == t) sel.push_back(rf_all[fi]);
            }
            if (!ref_subset.empty()) sel = rf_all;                  // --reference-all-fragments
            std::fprintf(stderr, "[reference-block] %zu of %zu fragments are local to target %zu\n",
                         sel.size(), rf_all.size(), t);
            // Per-fragment mass, so the factor score can be DIFFERENCED against the whole-locus one
            // fragment by fragment. A total tells you the factorisation is wrong; only the
            // per-fragment split says which fragments carry it and how they were classified.
            std::vector<double> fm, fc;
            const double v = reference_pair_loglik(ca, cb, sel, rp,
                                                   opt.incidence_path.empty() ? nullptr : &fm,
                                                   opt.incidence_path.empty() ? nullptr : &fc);
            if (!opt.incidence_path.empty()) {
                std::ofstream mf(opt.incidence_path);
                mf.precision(17);
                mf << "fragment\tfactor\ttarget\tlog_mass\tcontrib\n";
                for (std::size_t q = 0; q < sel.size(); ++q) {
                    mf << sel[q].name << "\tlocal\t" << t << '\t'
                       << (q < fm.size() ? fm[q] : 0.0) << '\t'
                       << (q < fc.size() ? fc[q] : 0.0) << '\n';
                }
            }
            std::printf("%.17g\n", v);
        }
        if (!ref_block_pair.empty()) {
            const std::size_t t = static_cast<std::size_t>(std::stoul(ref_block_pair[0]));
            if (t + 1 >= targets.size()) {
                throw std::runtime_error("genotype-frag: --reference-block-pair takes a TARGET INDEX "
                    "and needs a following target; there are " + std::to_string(targets.size()) +
                    " targets");
            }
            const int a1 = std::stoi(ref_block_pair[1]), b1 = std::stoi(ref_block_pair[2]);
            const int a2 = std::stoi(ref_block_pair[3]), b2 = std::stoi(ref_block_pair[4]);
            const std::string c1 = chain_span_sequence(blocks, targets, t, t + 1, {a1, b1}, opt.flank_bp);
            const std::string c2 = chain_span_sequence(blocks, targets, t, t + 1, {a2, b2}, opt.flank_bp);
            std::vector<Fragment> sel;
            for (std::size_t fi = 0; fi < rf_all.size(); ++fi) {
                if (inc.factor_of[fi] == FragmentFactor::Boundary &&
                    inc.targets_of[fi].front() == t &&
                    inc.targets_of[fi].back() == t + 1) sel.push_back(rf_all[fi]);
            }
            if (!ref_subset.empty()) sel = rf_all;
            std::fprintf(stderr,
                         "[reference-block-pair] %zu of %zu fragments span targets %zu/%zu\n",
                         sel.size(), rf_all.size(), t, t + 1);
            std::vector<double> fm, fc;
            const double v = reference_pair_loglik(c1, c2, sel, rp,
                                                   opt.incidence_path.empty() ? nullptr : &fm,
                                                   opt.incidence_path.empty() ? nullptr : &fc);
            if (!opt.incidence_path.empty()) {
                std::ofstream mf(opt.incidence_path);
                mf.precision(17);
                mf << "fragment\tfactor\ttarget\tlog_mass\tcontrib\n";
                for (std::size_t q = 0; q < sel.size(); ++q) {
                    mf << sel[q].name << "\tboundary\t" << t << '\t'
                       << (q < fm.size() ? fm[q] : 0.0) << '\t'
                       << (q < fc.size() ? fc[q] : 0.0) << '\n';
                }
            }
            std::printf("%.17g\n", v);
        }
        log.done();
        return 0;
    }

    if (!spell_calls.empty()) {
        std::ifstream cf(spell_calls);
        if (!cf) throw std::runtime_error("genotype-frag: cannot read " + spell_calls);
        std::string line;
        if (!std::getline(cf, line)) throw std::runtime_error("genotype-frag: empty " + spell_calls);

        // PROVENANCE. An allele index is an index into one exact catalogue -- the product of this
        // graph, this bubble decomposition and this exclusion set. Spelling a table against a
        // different catalogue produces a wrong sequence and no error, which cost one full LPA pilot:
        // the completion arm read 56042 edits from truth when the true figure was 21, because the
        // spelling run omitted the --exclude-haplotypes the scoring run had used.
        const std::string mine = allele_catalogue_fingerprint(blocks);
        if (line.rfind("# panvar-allele-catalogue", 0) == 0) {
            const std::size_t tab = line.find('\t');
            const std::string theirs = tab == std::string::npos ? std::string()
                                                               : line.substr(tab + 1);
            if (theirs != mine) {
                throw std::runtime_error(
                    "genotype-frag: --spell-calls table was produced against a DIFFERENT allele "
                    "catalogue (table " + theirs + ", this run " + mine + "). Allele indices are "
                    "only meaningful for one combination of graph, bubble decomposition and "
                    "--exclude-haplotypes; spelling across a mismatch yields a wrong sequence "
                    "silently. Re-run this command with the same -i / --bubble-prefix-in / "
                    "--exclude-haplotypes the call table was produced with.");
            }
            if (!std::getline(cf, line)) {
                throw std::runtime_error("genotype-frag: " + spell_calls + " has no header after "
                                         "its provenance line");
            }
        } else {
            // Production's own table carries no manifest, and comparing the two callers is a
            // supported use. Unverifiable is not the same as wrong, so this warns rather than
            // refusing -- but it says so, because the silent case is what caused the defect.
            std::fprintf(stderr,
                "[genotype-frag] WARNING: %s carries no allele-catalogue manifest, so its "
                "provenance cannot be checked. If it came from genotype-frag, re-run the scoring "
                "with a build that emits one. Allele indices are meaningless against a different "
                "graph, bubble decomposition or --exclude-haplotypes.\n", spell_calls.c_str());
        }
        std::vector<std::string> header;
        {
            std::size_t start = 0;
            for (std::size_t i = 0; i <= line.size(); ++i) {
                if (i == line.size() || line[i] == '\t') {
                    header.push_back(line.substr(start, i - start));
                    start = i + 1;
                }
            }
        }
        const auto col = [&](const std::string& n) {
            for (std::size_t i = 0; i < header.size(); ++i) if (header[i] == n) return static_cast<long>(i);
            throw std::runtime_error("genotype-frag: --spell-calls table has no column '" + n + "'");
        };
        // Both callers' tables are accepted: production names the column block_index, the prototype
        // names it block. Requiring one of them would make the comparison depend on which caller
        // happened to write the file.
        const long c_blk = std::find(header.begin(), header.end(), std::string("block_index")) != header.end()
            ? col("block_index") : col("block");
        const long c_a1 = col("allele1");
        const long c_a2 = col("allele2");
        std::vector<int> a1(chain.size(), -1), a2(chain.size(), -1);
        while (std::getline(cf, line)) {
            std::vector<std::string> f;
            std::size_t start = 0;
            for (std::size_t i = 0; i <= line.size(); ++i) {
                if (i == line.size() || line[i] == '\t') { f.push_back(line.substr(start, i - start)); start = i + 1; }
            }
            if (f.size() <= static_cast<std::size_t>(std::max(c_blk, std::max(c_a1, c_a2)))) continue;
            const long bi = std::stol(f[static_cast<std::size_t>(c_blk)]);
            if (bi < 0 || static_cast<std::size_t>(bi) >= chain.size()) continue;
            a1[static_cast<std::size_t>(bi)] = std::stoi(f[static_cast<std::size_t>(c_a1)]);
            a2[static_cast<std::size_t>(bi)] = std::stoi(f[static_cast<std::size_t>(c_a2)]);
        }
        std::string s1, s2;
        spell_called_pair(blocks, a1, a2, s1, s2);
        const std::string fa = out_prefix + ".called.fa";
        std::ofstream of(fa);
        if (!of) throw std::runtime_error("genotype-frag: cannot write " + fa);
        of << ">called_1\n" << s1 << "\n>called_2\n" << s2 << '\n';
        of.flush();
        if (!of) throw std::runtime_error("genotype-frag: write failed for " + fa);
        log.info("spelled " + spell_calls + ": " + std::to_string(s1.size()) + " bp and " +
                 std::to_string(s2.size()) + " bp");
        log.wrote({fa});
        log.done();
        return 0;
    }

    // ---- truth, resolved the same way `genotype` resolves it ---------------------------------
    std::vector<int> truth1, truth2;
    std::vector<std::string> truth_seq1, truth_seq2;
    bool have_truth = false;
    if (!truth_haplotypes.empty()) {
        const std::vector<std::string> names = split_commas(truth_haplotypes);
        if (names.size() != 2) {
            throw std::runtime_error("genotype-frag: --truth-haplotypes needs exactly two names");
        }
        const auto resolve = [&](const std::string& name, std::vector<int>& out_alleles,
                                 std::vector<std::string>* out_seq = nullptr) {
            out_alleles.assign(chain.size(), -1);
            if (out_seq != nullptr) out_seq->assign(chain.size(), std::string());
            const bool is_held =
                std::any_of(held_out.begin(), held_out.end(),
                            [&](const PathRecord& p) { return p.name == name; });
            for (std::size_t bi = 0; bi < chain.size(); ++bi) {
                if (!is_held) {
                    const auto it = blocks[bi].allele_of.find(name);
                    if (it != blocks[bi].allele_of.end()) {
                        out_alleles[bi] = static_cast<int>(it->second);
                        if (out_seq != nullptr && it->second < blocks[bi].allele_seq.size()) {
                            (*out_seq)[bi] = blocks[bi].allele_seq[it->second];
                        }
                    }
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
                if (out_seq != nullptr) (*out_seq)[bi] = seq;
                if (seq.empty()) continue;
                for (std::size_t ai = 0; ai < blocks[bi].allele_seq.size(); ++ai) {
                    if (blocks[bi].allele_seq[ai] == seq) {
                        out_alleles[bi] = static_cast<int>(ai);
                        break;
                    }
                }
            }
        };
        resolve(names[0], truth1, &truth_seq1);
        resolve(names[1], truth2, &truth_seq2);
        have_truth = true;
        std::size_t representable = 0;
        for (std::size_t bi = 0; bi < chain.size(); ++bi) {
            if (truth1[bi] >= 0 && truth2[bi] >= 0) ++representable;
        }
        log.info("truth " + names[0] + " / " + names[1] + ": both alleles representable at " +
                 std::to_string(representable) + "/" + std::to_string(chain.size()) + " blocks");
    }

    if (mosaic_floor) {
        if (!have_truth) throw std::runtime_error("genotype-frag: --mosaic-floor needs --truth-haplotypes");
        std::vector<double> pens;
        for (const std::string& t : split_commas(switch_penalties_arg)) pens.push_back(std::stod(t));
        std::vector<std::string> names;
        for (const PathRecord& p : panel_graph.paths) names.push_back(p.name);
        const std::string fp = out_prefix + ".mosaic_floor.tsv";
        std::ofstream mf(fp);
        if (!mf) throw std::runtime_error("genotype-frag: cannot write " + fp);
        mf << "homologue\tfloor\tswitch_penalty\tedits\tswitches\tblocks_scored\tblocks\n";
        std::size_t tot_complete = 0, tot_free = 0;
        std::vector<std::size_t> tot_pen(pens.size(), 0);
        for (int side = 0; side < 2; ++side) {
            const MosaicFloors r = mosaic_floors(blocks, names, side == 0 ? truth_seq1 : truth_seq2,
                                                 pens, opt.threads);
            const int hn = side + 1;
            mf << hn << "\tcomplete\tNA\t" << r.complete << "\tNA\t" << r.blocks_scored << '\t'
               << r.blocks << '\n';
            mf << hn << "\tfree\t0\t" << r.free_mosaic << '\t' << r.switches_free << '\t'
               << r.blocks_scored << '\t' << r.blocks << '\n';
            for (std::size_t i = 0; i < r.penalised.size(); ++i) {
                mf << hn << "\tpenalised\t" << r.penalised[i].first << '\t' << r.penalised[i].second
                   << '\t' << r.penalised_switches[i] << '\t' << r.blocks_scored << '\t'
                   << r.blocks << '\n';
                tot_pen[i] += r.penalised[i].second;
            }
            tot_complete += r.complete;
            tot_free += r.free_mosaic;
        }
        mf.flush();
        if (!mf) throw std::runtime_error("genotype-frag: write failed for " + fp);
        log.info("floors, both homologues: complete " + std::to_string(tot_complete) +
                 ", free mosaic " + std::to_string(tot_free) + " (H_mosaic " +
                 std::to_string(tot_complete > tot_free ? tot_complete - tot_free : 0) + ")");
        for (std::size_t i = 0; i < pens.size(); ++i) {
            log.info("  switch penalty " + std::to_string(pens[i]) + ": " + std::to_string(tot_pen[i]));
        }
        log.wrote({fp});
        log.done();
        return 0;
    }

    if (hopt.multiplicity_aware && !(hopt.mass_tolerance > 0.0 && hopt.mass_tolerance < 1.0)) {
        throw std::runtime_error("genotype-frag: --mass-tolerance must lie strictly between 0 and 1");
    }
    if (hap_mode && hopt.joint_depth && hopt.haploid_depth <= 0.0) {
        throw std::runtime_error(
            "genotype-frag: --joint-depth requires --haploid-depth. The depth rate must come from "
            "outside the candidate set -- invariant flanks, or genome-wide depth, or a known "
            "simulation rate. Fitting it from the alignment-best pair lets the expectation follow "
            "the hypothesis it is meant to test, and that is what the first implementation did");
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
        // Refuse to score before the decomposition is known to round-trip. Without this the caller
        // silently scores sequences the panel does not contain.
        verify_block_spelling(panel_graph, blocks, names);
        log.info("scoring " + std::to_string(names.size()) + " panel haplotypes (shortlist " +
                 std::to_string(hopt.max_haplotypes) + ") over " + std::to_string(frags.size()) +
                 " fragments");
        // THE HAPLOTYPE IS THE WALK. Spell every panel path from its own P-line steps -- the same
        // thing `odgi paths -f` does -- and hand those to the scorer as the authoritative sequence.
        // Concatenating per-block alleles is a decomposition artifact that is NOT guaranteed to
        // reproduce the walk: at cyp2d6 one path's block spelling is an exact 205236 bp prefix of a
        // 207214 bp walk, the chain stopping 1978 bp early with the tail dropped silently.
        //
        // Canonicalised into the BLOCK frame, because the block chain is reference-oriented and a
        // path antiparallel to the reference spells the reverse complement of its walk. The flip is
        // decided by comparison, never assumed.
        std::unordered_map<std::string, std::string> walk_sequences;
        {
            const auto by_name = path_records_by_name(graph);
            std::size_t flipped = 0, gapped = 0;
            for (const std::string& nm : names) {
                const auto it = by_name.find(nm);
                if (it == by_name.end() || it->second == nullptr) continue;
                bool complete = false;
                std::string w = spell_path_steps_sequence(graph, it->second->steps, &complete);
                if (!complete) continue;
                const std::string blk = spell_block_haplotype(blocks, nm);
                if (blk != w) {
                    // Reported, NOT converted. The scorer is strand-symmetric -- every fragment is
                    // aligned in both orientations, and a reverse-complemented panel is asserted to
                    // score identically -- so converting the walk into the block frame would buy
                    // nothing and would make the bytes we dump differ from the bytes we score.
                    // Strand-canonical deduplication handles the rest.
                    if (blk == reverse_complement(w)) ++flipped;
                    else ++gapped;
                }
                walk_sequences.emplace(nm, std::move(w));
            }
            log.info("scoring the GRAPH WALK for " + std::to_string(walk_sequences.size()) + "/" +
                     std::to_string(names.size()) + " panel paths (" + std::to_string(flipped) +
                     " antiparallel to the block frame and scored as-is, " + std::to_string(gapped) +
                     " where the block decomposition does not reproduce the walk and the walk wins)");
        }

        const HaplotypeResult hr =
            genotype_haplotype_pairs(chain, blocks, names, frags, hopt,
                                     have_truth ? &truth1 : nullptr,
                                     have_truth ? &truth2 : nullptr, top_pairs,
                                     &walk_sequences);
        if (!hr.top_pairs.empty()) {
            log.info("best pair: " + hr.shortlist[hr.top_pairs[0].hap1] + " / " +
                     hr.shortlist[hr.top_pairs[0].hap2] + " (posterior " +
                     std::to_string(hr.top_pairs[0].posterior) + ")");
        }
        if (hr.convergence.pairs_rescored > 0 && !hopt.joint_marginal) {
            log.info("joint assignment: " + std::to_string(hr.convergence.pairs_converged) + "/" +
                     std::to_string(hr.convergence.pairs_rescored) + " pairs reached a fixed point; "
                     "max iterations used " + std::to_string(hr.convergence.max_iterations_used) +
                     "; fragments still moving in the last pass " +
                     std::to_string(hr.convergence.total_moves_last_iteration));
        }
        log.info(std::to_string(hr.n_informative) + " of " + std::to_string(hr.n_fragments) +
                 " fragments discriminate between shortlisted haplotypes");
        if (hopt.unplaced_neutral) {
            // Reported, never assumed: this arm buys neutrality by DISCARDING evidence, and how much
            // it discards is the cost that has to be weighed against what it fixes.
            log.info("--unplaced-neutral dropped " + std::to_string(hr.n_neutralised) + " of " +
                     std::to_string(hr.n_fragments) + " fragments as unseeded on some shortlisted "
                     "haplotype (search never happened there)");
        }
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
        write_haplotype_results(out_prefix, hr, have_truth, allele_catalogue_fingerprint(blocks));
        {
            const auto& c = hr.completeness;
            const double kept_pct = c.clusters_found == 0 ? 100.0
                : 100.0 * static_cast<double>(c.clusters_kept) / static_cast<double>(c.clusters_found);
            const double anch_pct = c.anchor_occurrences_seen == 0 ? 0.0
                : 100.0 * static_cast<double>(c.anchor_occurrences_dropped) /
                  static_cast<double>(c.anchor_occurrences_seen);
            const double capped_pct = c.fragments_total == 0 ? 0.0
                : 100.0 * static_cast<double>(c.fragments_truncated) /
                  static_cast<double>(c.fragments_total);
            if (hopt.multiplicity_aware) {
                char buf[160];
                std::snprintf(buf, sizeof(buf),
                              "POST-RECRUITMENT pruning omitted: max %.3e, mean %.3e per "
                              "fragment-haplotype against --mass-tolerance %.1e (bound %s). This is "
                              "mass discarded from what recruitment FOUND -- it says nothing about "
                              "states recruitment never generated, which is a separate and currently "
                              "larger loss",
                              c.omitted_mass_max, c.omitted_mass_mean, hopt.mass_tolerance,
                              c.omitted_mass_max <= hopt.mass_tolerance ? "met" : "NOT MET");
                log.info(buf);
            }
            {
                // Reported whichever path ran, and SEPARATED: the Cartesian figure is what the
                // product would cost and grows as the square of the copy number, the join figure is
                // what was actually performed, and the rescue figures are neither -- the interval
                // scan is linear in the prior's width per anchored placement. Collapsed into one
                // "cost" number, a quadratic hidden behind a linear one is invisible.
                char jb[320];
                std::snprintf(jb, sizeof(jb),
                              "work: %llu Cartesian combinations hypothetical, %llu coordinate-join "
                              "probes actual, %llu rescue interval positions examined -> %llu "
                              "rescued placements; dispatch %llu join / %llu cartesian",
                              (unsigned long long)c.cartesian_combinations,
                              (unsigned long long)c.join_operations,
                              (unsigned long long)c.rescue_positions,
                              (unsigned long long)c.rescue_placements,
                              (unsigned long long)c.join_chosen,
                              (unsigned long long)c.cartesian_chosen);
                log.info(jb);
            }
            if (hopt.zero_seed_fallback) {
                char zb[320];
                std::snprintf(zb, sizeof(zb),
                              "no-primary-placement fallback (--zero-seed-*): ran for %llu "
                              "fragment-haplotype pairs (%llu pigeonhole, "
                              "%llu exhaustive); %llu candidate starts -> %llu verified -> %llu "
                              "placements",
                              (unsigned long long)c.zs_invocations,
                              (unsigned long long)c.zs_pigeonhole,
                              (unsigned long long)c.zs_exhaustive,
                              (unsigned long long)c.zs_candidate_starts,
                              (unsigned long long)c.zs_verified_starts,
                              (unsigned long long)c.zs_placements);
                log.info(zb);
            }
            if (hopt.joint_depth) {
                char cb[220];
                std::snprintf(cb, sizeof(cb),
                              "cost: %llu anchor hits, %llu mate combinations scored, %llu placements "
                              "before grouping, %llu groups after -- compression is DOWNSTREAM, so it "
                              "restores evidence without avoiding the work",
                              (unsigned long long)c.anchor_hits,
                              (unsigned long long)c.mate_combinations,
                              (unsigned long long)c.placements_before_grouping,
                              (unsigned long long)c.groups_after_grouping);
                log.info(cb);
            }
            log.info("post-anchor cluster retention " + std::to_string(static_cast<int>(kept_pct)) +
                     "% (NOT total completeness); anchor occurrences dropped by --max-anchor-occ " +
                     std::to_string(static_cast<int>(anch_pct)) + "%; fragments capped by "
                     "--placement-topk " + std::to_string(static_cast<int>(capped_pct)) + "%");
            if (c.recall_measured) {
                const double recall = 100.0 * static_cast<double>(c.truth_recovered) /
                                      static_cast<double>(c.truth_resolvable);
                const double spur = static_cast<double>(c.spurious_placements) /
                                    static_cast<double>(c.truth_resolvable);
                log.info("true-origin placement recall " + std::to_string(static_cast<int>(recall)) +
                         "% over " + std::to_string(c.truth_resolvable) + " fragments; " +
                         std::to_string(spur).substr(0, 4) + " additional placements per fragment "
                         "elsewhere on the origin haplotype (inside a repeat these are legitimate "
                         "alternatives, not errors). Nothing may be called an information limit "
                         "while recall is low");
            } else {
                log.info("true-origin recall not measurable here (needs reads simulated from "
                         "haplotypes that are IN the panel); retention alone cannot distinguish "
                         "restored true placements from added false ones");
            }
        }
        log.info("call margin " + std::to_string(hr.equivalence.margin) + " nats; " +
                 std::to_string(hr.equivalence.size) + " pair(s) within tolerance; " +
                 std::to_string(hr.equivalence.blocks_determined) + "/" +
                 std::to_string(hr.equivalence.blocks_total) + " blocks determined by all of them");
        log.wrote({out_prefix + ".hap_blocks.tsv", out_prefix + ".hap_pairs.tsv",
                   out_prefix + ".hap_scores.tsv", out_prefix + ".equivalence.tsv"});
        log.done();
        return 0;
    }

    // ---- which blocks ------------------------------------------------------------------------
    const std::vector<std::size_t> targets = build_targets(chain, blocks_arg, all_blocks);

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
