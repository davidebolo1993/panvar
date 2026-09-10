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
#include <atomic>
#include <cstdio>
#include <sstream>
#include <map>
#include <mutex>
#include <random>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <set>
#include <string>
#include <vector>

#include <sys/resource.h>

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
        << "      --bounded-search <out.tsv>  STAGE 1 (Hamming). For every (mate, strand,\n"
           "                              haplotype), find every placement within the divergence\n"
           "                              band by pigeonhole and verify each against the exhaustive\n"
           "                              scan. No occurrence cap, no top-k. Exits nonzero if the\n"
           "                              two ever disagree. O(panel x reference) -- for fixtures.\n"
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
    std::string bounded_search;
    bool bounded_verify = false;
    std::string interval_score, interval_cands;
    double interval_tau = 0.0, interval_tol = 1.0;
    bool interval_contrib = false;
    std::string interval_ckpt;
    std::size_t interval_batch = 500;
    bool interval_lazy = false, interval_screen_only = false;
    std::string interval_rounds;
    std::string ownership_table;
    std::string exposure_probe;
    std::string linkage_potential_out;
    bool linkage_selftest = false;
    bool hybrid_oracle = false;
    bool mapping_selftest = false;
    bool support_selftest = false;
    bool coordinate_selftest = false;
    bool budget_selftest = false;
    bool grouping_selftest = false;
    bool interval_selftest = false;
    bool normalisation_selftest = false;
    bool factor_selftest = false;
    bool contribution_selftest = false;
    bool hoinfer_selftest = false;
    bool completeness_selftest = false;
    bool activation_selftest = false;
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
        else if (a == "--bounded-search") bounded_search = value(i, a);
        else if (a == "--bounded-verify") bounded_verify = true;
        else if (a == "--interval-score") interval_score = value(i, a);
        else if (a == "--interval-candidates") interval_cands = value(i, a);
        else if (a == "--interval-tau") interval_tau = std::stod(value(i, a));
        else if (a == "--interval-tol") interval_tol = std::stod(value(i, a));
        else if (a == "--interval-contrib") interval_contrib = true;
        else if (a == "--interval-checkpoint") interval_ckpt = value(i, a);
        else if (a == "--interval-batch") interval_batch = cli::parse_size_arg(a, value(i, a));
        else if (a == "--interval-lazy") interval_lazy = true;
        else if (a == "--interval-screen-only") { interval_lazy = true; interval_screen_only = true; }
        else if (a == "--interval-rounds") interval_rounds = value(i, a);
        else if (a == "--ownership-table") ownership_table = value(i, a);
        else if (a == "--exposure-probe") exposure_probe = value(i, a);
        else if (a == "--linkage-potential") linkage_potential_out = value(i, a);
        else if (a == "--linkage-selftest") linkage_selftest = true;
        else if (a == "--hybrid-oracle") hybrid_oracle = true;
        else if (a == "--mapping-selftest") mapping_selftest = true;
        else if (a == "--support-selftest") support_selftest = true;
        else if (a == "--coordinate-selftest") coordinate_selftest = true;
        else if (a == "--budget-selftest") budget_selftest = true;
        else if (a == "--grouping-selftest") grouping_selftest = true;
        else if (a == "--interval-selftest") interval_selftest = true;
        else if (a == "--normalisation-selftest") normalisation_selftest = true;
        else if (a == "--factor-selftest") factor_selftest = true;
        else if (a == "--contribution-selftest") contribution_selftest = true;
        else if (a == "--hoinfer-selftest") hoinfer_selftest = true;
        else if (a == "--completeness-selftest") completeness_selftest = true;
        else if (a == "--activation-selftest") activation_selftest = true;
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
        else if (a == "--dump-containment") hopt.dump_containment = value(i, a);
        else if (a == "--insert-sigmas") opt.insert_sigmas = std::stoi(value(i, a));
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
        // EVERY PRIOR-DEFINING FIELD, copied. Letting ReferenceParams keep its own defaults for
        // any of these is how the oracle came to integrate over [100,270] while the accelerated
        // scorer used [100,230]: the two produced different numbers for the same reads and the
        // difference was read as a modelling failure for an entire session.
        rp.insert_sigmas = hopt.insert_sigmas;
        rp.discordant_rate = hopt.discordant_rate;
        rp.allow_overlapping_pairs = hopt.allow_overlapping_pairs;
        // AND THE PRIORS MUST THEN BE IDENTICAL. Asserted, not assumed: a differential between two
        // arms that integrate different insert supports is not a measurement, so this refuses
        // rather than printing a score.
        {
            const long fl_p = fragment_insert_floor(frags, hopt.allow_overlapping_pairs);
            const long fl_r = fragment_insert_floor(frags, rp.allow_overlapping_pairs);
            const InsertPrior pp = make_insert_prior(hopt.fragment_len, hopt.fragment_sd,
                                                     hopt.discordant_rate, hopt.insert_sigmas,
                                                     fl_p);
            const InsertPrior pr2 = make_insert_prior(rp.fragment_len, rp.fragment_sd,
                                                      rp.discordant_rate, rp.insert_sigmas, fl_r);
            bool same = (fl_p == fl_r) && (pp.lo == pr2.lo) && (pp.hi == pr2.hi) &&
                        (pp.logp.size() == pr2.logp.size());
            for (std::size_t q = 0; same && q < pp.logp.size(); ++q)
                if (std::abs(pp.logp[q] - pr2.logp[q]) > 1e-12) same = false;
            if (same && (std::abs(pp.log_residual_below - pr2.log_residual_below) > 1e-12 ||
                         std::abs(pp.log_residual_above - pr2.log_residual_above) > 1e-12))
                same = false;
            if (!same) {
                throw std::runtime_error(
                    "genotype-frag: the accelerated and reference arms would use DIFFERENT insert "
                    "priors (floor " + std::to_string(fl_p) + " vs " + std::to_string(fl_r) +
                    ", support " + std::to_string(pp.lo) + "-" + std::to_string(pp.hi) + " vs " +
                    std::to_string(pr2.lo) + "-" + std::to_string(pr2.hi) +
                    "); a differential between different insert-state universes is not a "
                    "measurement, so no score is produced");
            }
        }
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
    // --linkage-selftest builds its own geometry and emissions, so it needs no reads. Requiring a
    // fixture would defeat it: the neutrality cases it covers -- an edge with no fragments, and one
    // whose emissions are identical everywhere -- are exactly what a read fixture cannot produce on
    // demand.
    if (read_paths.empty() && spell_calls.empty() && !mosaic_floor && dump_sequences.empty() &&
        spell_pair.empty() && ref_block.empty() && ref_block_pair.empty() &&
        origin_universe.empty() && reconcile_scope.empty() && !linkage_selftest &&
        !hybrid_oracle && !mapping_selftest && !completeness_selftest && !activation_selftest &&
        !support_selftest && !coordinate_selftest && !budget_selftest && !grouping_selftest &&
        !interval_selftest && !normalisation_selftest && !factor_selftest &&
        !hoinfer_selftest && !contribution_selftest) {
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
    // ---- HIGHER-ORDER INFERENCE, UNOPTIMISED, AGAINST COMPLETE BRUTE FORCE ---------------------
    // The exact model is
    //     prod_b U_b(X_b) * prod_b T(X_{b-1}, X_b) * F1(c1..c4) * F2(c3..c5)
    // with X_b an ORDERED PAIR of panel-template identities. Li-Stephens permits a different
    // template at each block, so F1 is genuinely FOURTH order: grouping makes its lookup cheap but
    // does not reduce its order, and the stay term (1-r)I depends on exact template identity, so a
    // message may not be collapsed to signature classes without a contraction proof.
    //
    // This establishes what any such contraction must REPRODUCE. It carries the full history the
    // factors need and is deliberately unoptimised; the reference enumerates every path.
    if (hoinfer_selftest) {
        std::size_t fails = 0;
        const auto ok_ = [&](bool c, const std::string& w) {
            std::printf("%s\t%s\n", c ? "ok" : "FAIL", w.c_str());
            if (!c) ++fails;
        };
        const auto sci = [](double x) {
            char b[32]; std::snprintf(b, sizeof b, "%.3e", x); return std::string(b);
        };
        std::mt19937_64 rng(20260909);
        const std::size_t NH = 3, NB = 6;              // templates, blocks
        const std::size_t NS = NH * NH;                // ordered diploid template pairs
        // Each template's allele at each block, and the factors' spans.
        // THREE alleles per block, mapped by a NON-IDENTITY function to TWO signature classes, so
        // the class map is genuinely a map and not a relabelling -- and so several templates share
        // a class while carrying different unary weights, which is the case where collapsing to
        // classes for the MESSAGE (rather than only for the factor lookup) would be wrong.
        std::vector<std::vector<std::uint32_t>> tal(NH, std::vector<std::uint32_t>(NB, 0));
        for (std::size_t t = 0; t < NH; ++t)
            for (std::size_t b = 0; b < NB; ++b)
                tal[t][b] = static_cast<std::uint32_t>(rng() % 3);   // three alleles per block
        // allele -> class: {0,2} -> 0, {1} -> 1. Non-identity, non-injective.
        const auto cls = [](std::uint32_t a) { return static_cast<std::uint32_t>(a == 1 ? 1 : 0); };
        const std::vector<std::size_t> span1 = {1, 2, 3, 4};   // F1
        const std::vector<std::size_t> span2 = {3, 4, 5};      // F2, overlapping on 3 and 4
        // Unary weights per diploid state per block, deliberately asymmetric.
        std::vector<std::vector<double>> U(NB, std::vector<double>(NS, 0.0));
        for (std::size_t b = 0; b < NB; ++b)
            for (std::size_t x = 0; x < NS; ++x)
                U[b][x] = std::uniform_real_distribution<double>(0.2, 1.8)(rng);
        // The two factors, as functions of the ORDERED allele tuples of the two homologues, made
        // invariant under GLOBAL SWAP so they are legitimate phase factors.
        const auto key_of = [&](const std::vector<std::size_t>& span,
                                const std::vector<std::size_t>& xs) {
            std::string k;
            for (std::size_t q = 0; q < span.size(); ++q) {
                const std::size_t x = xs[q];
                k.push_back(static_cast<char>('0' + cls(tal[x / NH][span[q]])));
            }
            k.push_back('|');
            for (std::size_t q = 0; q < span.size(); ++q) {
                const std::size_t x = xs[q];
                k.push_back(static_cast<char>('0' + cls(tal[x % NH][span[q]])));
            }
            return k;
        };
        std::unordered_map<std::string, double> F1v, F2v;
        const auto factor = [&](std::unordered_map<std::string, double>& M,
                                const std::vector<std::size_t>& span,
                                const std::vector<std::size_t>& xs) {
            std::string k = key_of(span, xs);
            // Global swap: exchange the two homologues everywhere.
            const std::size_t half = k.find('|');
            std::string sw = k.substr(half + 1) + "|" + k.substr(0, half);
            const std::string canon = std::min(k, sw);
            auto it = M.find(canon);
            if (it == M.end())
                it = M.emplace(canon,
                               std::uniform_real_distribution<double>(0.3, 3.0)(rng)).first;
            return it->second;
        };
        double r = 0.25;
        const auto Tw = [&](std::size_t from, std::size_t to) {
            const double stay = (from == to) ? (1.0 - r) : 0.0;
            return stay + r / static_cast<double>(NH);
        };
        const auto Tdip = [&](std::size_t xp, std::size_t x) {
            return Tw(xp / NH, x / NH) * Tw(xp % NH, x % NH);
        };
        // Every arm: r at both extremes and between, and a SYMMETRIC-emission arm in which the
        // unaries are invariant under exchanging the two homologues -- the only setting where
        // homologue-swap invariance of the posterior is a property of the KERNEL rather than of
        // the fixture.
        struct Arm { double r; bool symmetric; const char* name; };
        const Arm arms[] = {{0.0, false, "r=0"}, {0.25, false, "r=0.25"}, {1.0, false, "r=1"},
                            {0.25, true, "r=0.25 symmetric emissions"}};
        std::vector<std::vector<double>> Usave = U;
        for (const Arm& arm : arms) {
        r = arm.r;
        U = Usave;
        if (arm.symmetric) {
            for (std::size_t b = 0; b < NB; ++b)
                for (std::size_t i = 0; i < NH; ++i)
                    for (std::size_t j = 0; j < NH; ++j)
                        U[b][j * NH + i] = U[b][i * NH + j];
        }
        // ---- COMPLETE BRUTE FORCE: every path of diploid states -------------------------------
        std::vector<std::vector<double>> marg_bf(NB, std::vector<double>(NS, 0.0));
        double Z_bf = 0.0;
        {
            std::vector<std::size_t> path(NB, 0);
            std::size_t total = 1;
            for (std::size_t b = 0; b < NB; ++b) total *= NS;
            for (std::size_t code = 0; code < total; ++code) {
                std::size_t c = code;
                for (std::size_t b = 0; b < NB; ++b) { path[b] = c % NS; c /= NS; }
                double w = 1.0;
                for (std::size_t b = 0; b < NB; ++b) w *= U[b][path[b]];
                for (std::size_t b = 1; b < NB; ++b) w *= Tdip(path[b - 1], path[b]);
                std::vector<std::size_t> xs1, xs2;
                for (std::size_t q : span1) xs1.push_back(path[q]);
                for (std::size_t q : span2) xs2.push_back(path[q]);
                w *= factor(F1v, span1, xs1);
                w *= factor(F2v, span2, xs2);
                Z_bf += w;
                for (std::size_t b = 0; b < NB; ++b) marg_bf[b][path[b]] += w;
            }
        }
        ok_(Z_bf > 0.0, std::string(arm.name) + ": brute force over " + std::to_string(NS) + "^" +
                        std::to_string(NB) + " = " +
                        std::to_string(static_cast<std::size_t>(std::pow(NS, NB))) +
                        " diploid paths gives Z = " + sci(Z_bf));
        // ---- UNOPTIMISED INFERENCE: carry exactly the history the factors need ----------------
        // Eliminating left to right, the running scope is the set of block variables an
        // un-applied factor still needs. F1 needs blocks 1..4, F2 needs 3..5, so the maximum
        // carried history is three earlier variables plus the current one.
        std::vector<std::vector<double>> marg_inf(NB, std::vector<double>(NS, 0.0));
        double Z_inf = 0.0;
        {
            // State: (x1, x2, x3, x_cur) as needed. Implemented as a map from the tuple of
            // retained variables to weight, advanced block by block.
            std::map<std::vector<std::size_t>, double> msg;
            for (std::size_t x = 0; x < NS; ++x) msg[{x}] = U[0][x];
            for (std::size_t b = 1; b < NB; ++b) {
                std::map<std::vector<std::size_t>, double> nxt;
                for (const auto& kv : msg) {
                    for (std::size_t x = 0; x < NS; ++x) {
                        double w = kv.second * Tdip(kv.first.back(), x) * U[b][x];
                        if (w == 0.0) continue;
                        std::vector<std::size_t> hist = kv.first;
                        hist.push_back(x);
                        // Apply a factor as soon as its last block is reached.
                        if (b == span1.back()) {
                            std::vector<std::size_t> xs;
                            for (std::size_t q : span1) xs.push_back(hist[q]);
                            w *= factor(F1v, span1, xs);
                        }
                        if (b == span2.back()) {
                            std::vector<std::size_t> xs;
                            for (std::size_t q : span2) xs.push_back(hist[q]);
                            w *= factor(F2v, span2, xs);
                        }
                        nxt[hist] += w;
                    }
                }
                msg.swap(nxt);
            }
            for (const auto& kv : msg) {
                Z_inf += kv.second;
                for (std::size_t b = 0; b < NB; ++b) marg_inf[b][kv.first[b]] += kv.second;
            }
        }
        ok_(std::abs(Z_inf - Z_bf) / Z_bf < 1e-12,
            std::string(arm.name) +
            ": the unoptimised inference reproduces the partition weight (" + sci(Z_inf) +
            " vs " + sci(Z_bf) + ", relative " + sci(std::abs(Z_inf - Z_bf) / Z_bf) + ")");
        double worst_m = 0.0;
        for (std::size_t b = 0; b < NB; ++b)
            for (std::size_t x = 0; x < NS; ++x)
                worst_m = std::max(worst_m,
                                   std::abs(marg_inf[b][x] - marg_bf[b][x]) / Z_bf);
        ok_(worst_m < 1e-12, std::string(arm.name) + ": and every block marginal, over all " +
                             std::to_string(NB * NS) + " (worst relative " + sci(worst_m) + ")");
        // HOMOLOGUE-SWAP INVARIANCE, only meaningful on the symmetric arm: with asymmetric unaries
        // (i,j) and (j,i) carry different mass with NO factors at all, so a swap check there would
        // measure the fixture.
        if (arm.symmetric) {
            double worst_sw = 0.0;
            for (std::size_t b = 0; b < NB; ++b)
                for (std::size_t i = 0; i < NH; ++i)
                    for (std::size_t j = 0; j < NH; ++j)
                        worst_sw = std::max(worst_sw,
                                            std::abs(marg_bf[b][i * NH + j] -
                                                     marg_bf[b][j * NH + i]) / Z_bf);
            ok_(worst_sw < 1e-12,
                "homologue-swap invariance holds on SYMMETRIC emissions (worst relative " +
                sci(worst_sw) + ")");
        }
        // ---- THE STRUCTURED CONTRACTION -------------------------------------------------------
        // T = (1-r)I + (r/n)11^T. The STAY component leaves the template and the open run
        // untouched; the SWITCH component closes the run -- appending the classes that run
        // contributed, which follow from the template it held -- and redraws the template
        // uniformly. Exact template identity is preserved for as long as a run is OPEN, which is
        // what MH6 shows is required; only CLOSED runs are reduced to class history.
        //
        // THE OPERATION COUNTER counts MESSAGE-ENTRY UPDATES: one per accumulation into a state of
        // the next message, each being a multiply and an add in linear space. It is not a count of
        // floating-point instructions, and it is reported as such.
        {
            struct St { std::vector<std::uint32_t> h1, h2; std::size_t x; };
            // clamp_b < 0 means "no clamp"; otherwise block clamp_b is restricted to state clamp_x,
            // so the returned Z IS that block's unnormalised marginal. Marginals therefore come
            // from the SAME contraction rather than from a second code path.
            // Shared by both directions, so the two passes cannot disagree about state identity.
            const auto keyof = [](const St& s2) {
                std::string k;
                for (std::uint32_t v : s2.h1) k.push_back(static_cast<char>('a' + v));
                k.push_back('|');
                for (std::uint32_t v : s2.h2) k.push_back(static_cast<char>('a' + v));
                k.push_back('|');
                k += std::to_string(s2.x);
                return k;
            };
            std::uint64_t ops_tot = 0, peak_tot = 0;
            const auto contract = [&](int clamp_b, std::size_t clamp_x) {
            // Span 1 is blocks 1..4; span 2 is 3..5. Run the contraction over the WHOLE chain,
            // closing each factor when its last block is passed.
            std::map<std::string, std::pair<St, double>> msg;
            {
                for (std::size_t x = 0; x < NS; ++x) {
                    if (clamp_b == 0 && x != clamp_x) continue;
                    St s2; s2.x = x;
                    msg[keyof(s2)] = {s2, U[0][x]};
                }
            }
            std::uint64_t ops = 0, peak_entries = msg.size();
            for (std::size_t b = 1; b < NB; ++b) {
                std::map<std::string, std::pair<St, double>> nxt;
                for (const auto& kv : msg) {
                    const St& cur = kv.second.first;
                    const double w = kv.second.second;
                    const std::size_t t1 = cur.x / NH, t2 = cur.x % NH;
                    // The four stay/switch combinations of the ordered diploid edge.
                    for (int c1 = 0; c1 < 2; ++c1) {
                        for (int c2 = 0; c2 < 2; ++c2) {
                            const double w1 = (c1 == 0) ? (1.0 - r) : (r / double(NH));
                            const double w2 = (c2 == 0) ? (1.0 - r) : (r / double(NH));
                            if (w1 == 0.0 || w2 == 0.0) continue;
                            // Which templates the new state can take.
                            const std::size_t lo1 = (c1 == 0) ? t1 : 0;
                            const std::size_t hi1 = (c1 == 0) ? t1 + 1 : NH;
                            const std::size_t lo2 = (c2 == 0) ? t2 : 0;
                            const std::size_t hi2 = (c2 == 0) ? t2 + 1 : NH;
                            for (std::size_t n1 = lo1; n1 < hi1; ++n1) {
                                for (std::size_t n2 = lo2; n2 < hi2; ++n2) {
                                    St s2 = cur;
                                    // A SWITCH closes the open run: the classes it contributed at
                                    // every position from its start to b-1 are appended, and they
                                    // follow from the template that run held.
                                    if (c1 == 1)
                                        for (std::size_t q = s2.h1.size(); q < b; ++q)
                                            s2.h1.push_back(cls(tal[t1][q]));
                                    if (c2 == 1)
                                        for (std::size_t q = s2.h2.size(); q < b; ++q)
                                            s2.h2.push_back(cls(tal[t2][q]));
                                    s2.x = n1 * NH + n2;
                                    if (clamp_b == static_cast<int>(b) && s2.x != clamp_x) continue;
                                    double nw = w * w1 * w2 * U[b][s2.x];
                                    // Close a factor the moment its last block is passed: the open
                                    // run's remaining classes come from the CURRENT template.
                                    const auto close = [&](const std::vector<std::size_t>& span,
                                                           std::unordered_map<std::string, double>& M) {
                                        std::vector<std::size_t> xs;
                                        for (std::size_t q : span) {
                                            const std::uint32_t a1 = q < s2.h1.size()
                                                ? s2.h1[q] : cls(tal[n1][q]);
                                            const std::uint32_t a2 = q < s2.h2.size()
                                                ? s2.h2[q] : cls(tal[n2][q]);
                                            xs.push_back(a1 * NH + a2);   // placeholder, see below
                                            (void)a1; (void)a2;
                                        }
                                        return xs;
                                    };
                                    (void)close;
                                    if (b == span1.back() || b == span2.back()) {
                                        const auto& span = (b == span1.back()) ? span1 : span2;
                                        auto& M = (b == span1.back()) ? F1v : F2v;
                                        std::string k;
                                        for (std::size_t q : span)
                                            k.push_back(static_cast<char>('0' +
                                                (q < s2.h1.size() ? s2.h1[q] : cls(tal[n1][q]))));
                                        k.push_back('|');
                                        for (std::size_t q : span)
                                            k.push_back(static_cast<char>('0' +
                                                (q < s2.h2.size() ? s2.h2[q] : cls(tal[n2][q]))));
                                        const std::size_t half = k.find('|');
                                        const std::string sw =
                                            k.substr(half + 1) + "|" + k.substr(0, half);
                                        const std::string canon = std::min(k, sw);
                                        auto it = M.find(canon);
                                        if (it == M.end())
                                            it = M.emplace(canon,
                                                std::uniform_real_distribution<double>(0.3, 3.0)(rng)).first;
                                        nw *= it->second;
                                    }
                                    // Both factors closed: the histories are spent.
                                    if (b == span2.back()) { s2.h1.clear(); s2.h2.clear(); }
                                    ++ops;
                                    const std::string kk = keyof(s2);
                                    auto f = nxt.find(kk);
                                    if (f == nxt.end()) nxt.emplace(kk, std::make_pair(s2, nw));
                                    else f->second.second += nw;
                                }
                            }
                        }
                    }
                }
                msg.swap(nxt);
                peak_entries = std::max<std::uint64_t>(peak_entries, msg.size());
            }
            double Zc = 0.0;
            for (const auto& kv : msg) Zc += kv.second.second;
            ops_tot += ops;
            peak_tot = std::max(peak_tot, peak_entries);
            return Zc;
            };
            const double Z_con = contract(-1, 0);
            ok_(std::abs(Z_con - Z_bf) / Z_bf < 1e-12,
                std::string(arm.name) + ": the STRUCTURED CONTRACTION reproduces the partition "
                "weight (" + sci(Z_con) + " vs " + sci(Z_bf) + ", relative " +
                sci(std::abs(Z_con - Z_bf) / Z_bf) + "); peak " + std::to_string(peak_tot) +
                " message entries, " + std::to_string(ops_tot) + " message-entry updates");
            // ---- BACKWARD CONTRACTION AND THE FORWARD/BACKWARD JOIN ----------------------
            // Clamping proves correctness but is not the production marginal algorithm: it costs
            // one full contraction per (block, state). Production needs messages from both ends.
            //
            // A FACTOR STRADDLING BLOCK b IS APPLIED AT THE JOIN, not in either pass, or it would
            // be counted twice. The forward message carries the classes of positions whose run has
            // closed to the LEFT, the backward message those closed to the RIGHT, and the open run
            // on both sides is covered by the SAME X_b -- which is what makes the two histories
            // joinable into one class tuple.
            const auto contract_back = [&](int clamp_b, std::size_t clamp_x) {
                std::map<std::string, std::pair<St, double>> msg;
                for (std::size_t x = 0; x < NS; ++x) {
                    if (clamp_b == static_cast<int>(NB) - 1 && x != clamp_x) continue;
                    St s2; s2.x = x;
                    msg[keyof(s2)] = {s2, U[NB - 1][x]};
                }
                for (std::size_t bb = NB - 1; bb-- > 0;) {
                    std::map<std::string, std::pair<St, double>> nxt;
                    for (const auto& kv : msg) {
                        const St& cur = kv.second.first;
                        const double w = kv.second.second;
                        const std::size_t t1 = cur.x / NH, t2 = cur.x % NH;
                        for (int c1 = 0; c1 < 2; ++c1) {
                            for (int c2 = 0; c2 < 2; ++c2) {
                                const double w1 = (c1 == 0) ? (1.0 - r) : (r / double(NH));
                                const double w2 = (c2 == 0) ? (1.0 - r) : (r / double(NH));
                                if (w1 == 0.0 || w2 == 0.0) continue;
                                const std::size_t lo1 = (c1 == 0) ? t1 : 0;
                                const std::size_t hi1 = (c1 == 0) ? t1 + 1 : NH;
                                const std::size_t lo2 = (c2 == 0) ? t2 : 0;
                                const std::size_t hi2 = (c2 == 0) ? t2 + 1 : NH;
                                for (std::size_t n1 = lo1; n1 < hi1; ++n1) {
                                    for (std::size_t n2 = lo2; n2 < hi2; ++n2) {
                                        St s2 = cur;
                                        // Going LEFT, a switch closes the run that extended to the
                                        // right; its classes are prepended.
                                        if (c1 == 1)
                                            for (std::size_t q = NB - 1 - s2.h1.size(); q > bb; --q)
                                                s2.h1.push_back(cls(tal[t1][q]));
                                        if (c2 == 1)
                                            for (std::size_t q = NB - 1 - s2.h2.size(); q > bb; --q)
                                                s2.h2.push_back(cls(tal[t2][q]));
                                        s2.x = n1 * NH + n2;
                                        if (clamp_b == static_cast<int>(bb) && s2.x != clamp_x)
                                            continue;
                                        double nw = w * w1 * w2 * U[bb][s2.x];
                                        // A factor closes going left when its FIRST block is
                                        // reached.
                                        for (int fsel = 0; fsel < 2; ++fsel) {
                                            const auto& span = fsel == 0 ? span1 : span2;
                                            if (bb != span.front()) continue;
                                            auto& M = fsel == 0 ? F1v : F2v;
                                            std::string k;
                                            for (std::size_t q : span) {
                                                const std::size_t back = NB - 1 - q;
                                                k.push_back(static_cast<char>('0' +
                                                    (back < s2.h1.size() ? s2.h1[back]
                                                                         : cls(tal[n1][q]))));
                                            }
                                            k.push_back('|');
                                            for (std::size_t q : span) {
                                                const std::size_t back = NB - 1 - q;
                                                k.push_back(static_cast<char>('0' +
                                                    (back < s2.h2.size() ? s2.h2[back]
                                                                         : cls(tal[n2][q]))));
                                            }
                                            const std::size_t half = k.find('|');
                                            const std::string sw =
                                                k.substr(half + 1) + "|" + k.substr(0, half);
                                            const std::string canon = std::min(k, sw);
                                            auto it = M.find(canon);
                                            if (it == M.end())
                                                it = M.emplace(canon,
                                                    std::uniform_real_distribution<double>(0.3, 3.0)(rng)).first;
                                            nw *= it->second;
                                        }
                                        if (bb == span1.front()) { s2.h1.clear(); s2.h2.clear(); }
                                        const std::string kk = keyof(s2);
                                        auto f = nxt.find(kk);
                                        if (f == nxt.end()) nxt.emplace(kk, std::make_pair(s2, nw));
                                        else f->second.second += nw;
                                    }
                                }
                            }
                        }
                    }
                    msg.swap(nxt);
                }
                double Zb = 0.0;
                for (const auto& kv : msg) Zb += kv.second.second;
                return Zb;
            };
            // The backward contraction must reproduce the same partition weight, computed by an
            // independent traversal in the opposite direction.
            const double Z_back = contract_back(-1, 0);
            ok_(std::abs(Z_back - Z_bf) / Z_bf < 1e-12,
                std::string(arm.name) + ": the BACKWARD contraction independently reproduces the "
                "partition weight (" + sci(Z_back) + " vs " + sci(Z_bf) + ", relative " +
                sci(std::abs(Z_back - Z_bf) / Z_bf) + ")");
            double worst_bm = 0.0;
            for (std::size_t b = 0; b < NB; ++b)
                for (std::size_t x = 0; x < NS; ++x)
                    worst_bm = std::max(worst_bm,
                                        std::abs(contract_back(static_cast<int>(b), x) -
                                                 marg_bf[b][x]) / Z_bf);
            ok_(worst_bm < 1e-12,
                std::string(arm.name) + ": and the BACKWARD marginals match brute force over all " +
                std::to_string(NB * NS) + " (worst relative " + sci(worst_bm) + ")");

            // ---- THE ADJOINT SWEEP: marginals from the FORWARD circuit ---------------------
            // The forward contraction is an arithmetic circuit for Z, and Z is LINEAR in every
            // unary, so a reverse-mode sweep through that same circuit yields every marginal:
            //
            //     marginal(b, x) = SUM over states s at block b with X(s) = x of  msg_b(s) * bar_b(s)
            //
            // where bar_b(s) = dZ/d msg_b(s) is the weight of everything downstream. Each forward
            // update  dst += src * m  has adjoint  bar_src += bar_dst * m, so the sweep reuses the
            // forward history organisation and never builds the right-to-left history that makes an
            // independent backward contraction expensive. The independent backward implementation
            // stays as an ORACLE; it is not the production algorithm.
            {
                // Forward again, retaining every message and its update list.
                struct UpdK { std::string src, dst; double m; };
                std::vector<std::map<std::string, std::pair<St, double>>> keep;
                std::vector<std::vector<UpdK>> tape;
                std::map<std::string, std::pair<St, double>> msg;
                for (std::size_t x = 0; x < NS; ++x) { St s2; s2.x = x; msg[keyof(s2)] = {s2, U[0][x]}; }
                keep.push_back(msg);
                std::uint64_t adj_ops = 0;
                std::uint64_t fw_mult = 0, fw_lookup = 0, fw_hist = 0;
                for (std::size_t b = 1; b < NB; ++b) {
                    std::map<std::string, std::pair<St, double>> nxt;
                    std::vector<UpdK> tpk;
                    for (const auto& kv : msg) {
                        const St& cur = kv.second.first; const double w = kv.second.second;
                        const std::size_t t1 = cur.x / NH, t2 = cur.x % NH;
                        for (int c1 = 0; c1 < 2; ++c1) for (int c2 = 0; c2 < 2; ++c2) {
                            const double w1 = (c1 == 0) ? (1.0 - r) : (r / double(NH));
                            const double w2 = (c2 == 0) ? (1.0 - r) : (r / double(NH));
                            if (w1 == 0.0 || w2 == 0.0) continue;
                            const std::size_t lo1 = (c1 == 0) ? t1 : 0, hi1 = (c1 == 0) ? t1 + 1 : NH;
                            const std::size_t lo2 = (c2 == 0) ? t2 : 0, hi2 = (c2 == 0) ? t2 + 1 : NH;
                            for (std::size_t n1 = lo1; n1 < hi1; ++n1)
                            for (std::size_t n2 = lo2; n2 < hi2; ++n2) {
                                St s2 = cur;
                                if (c1 == 1) for (std::size_t q = s2.h1.size(); q < b; ++q)
                                    { s2.h1.push_back(cls(tal[t1][q])); ++fw_hist; }
                                if (c2 == 1) for (std::size_t q = s2.h2.size(); q < b; ++q)
                                    { s2.h2.push_back(cls(tal[t2][q])); ++fw_hist; }
                                s2.x = n1 * NH + n2;
                                double m = w1 * w2 * U[b][s2.x];
                                ++fw_mult;
                                if (b == span1.back() || b == span2.back()) {
                                    const auto& span = (b == span1.back()) ? span1 : span2;
                                    auto& M = (b == span1.back()) ? F1v : F2v;
                                    std::string k;
                                    for (std::size_t q : span) k.push_back(static_cast<char>('0' +
                                        (q < s2.h1.size() ? s2.h1[q] : cls(tal[n1][q]))));
                                    k.push_back('|');
                                    for (std::size_t q : span) k.push_back(static_cast<char>('0' +
                                        (q < s2.h2.size() ? s2.h2[q] : cls(tal[n2][q]))));
                                    const std::size_t half = k.find('|');
                                    const std::string sw = k.substr(half+1) + "|" + k.substr(0,half);
                                    const std::string canon = std::min(k, sw);
                                    auto it = M.find(canon);
                                    ++fw_lookup;
                                    if (it == M.end()) it = M.emplace(canon,
                                        std::uniform_real_distribution<double>(0.3,3.0)(rng)).first;
                                    m *= it->second;
                                }
                                if (b == span2.back()) { s2.h1.clear(); s2.h2.clear(); }
                                const std::string kk = keyof(s2);
                                auto f = nxt.find(kk);
                                if (f == nxt.end()) nxt.emplace(kk, std::make_pair(s2, w * m));
                                else f->second.second += w * m;
                                // The tape records the DESTINATION KEY, not an index into a map
                                // that is still growing: a position taken mid-insert would be
                                // invalidated by the next insertion.
                                tpk.push_back({kv.first, kk, m});
                                ++adj_ops;
                            }
                        }
                    }
                    tape.push_back(std::move(tpk));
                    msg.swap(nxt);
                    keep.push_back(msg);
                }
                // ---- THE TAPE-FREE ADJOINT --------------------------------------------------
                // The tape above costs one entry per forward update -- at C4 scale ~1.9e9 entries,
                // dwarfing everything else -- so production must RECOMPUTE each multiplier by
                // re-walking the same loop over the retained messages. Same message-entry update
                // count; different CPU work, which is counted separately below.
                std::uint64_t fwd_updates = adj_ops;   // everything so far was the forward pass

                // The count is PREDICTABLE before either sweep runs: each retained state expands
                // into (stay ? 1 : 0) + (switch ? NH : 0) targets per homologue, and the forward
                // pass inserts every one of them, so the adjoint -- enumerating the same loop in
                // reverse -- must visit exactly the same number. An inequality here means the two
                // enumerations have drifted apart, which no marginal check would necessarily catch
                // if the drift happened to be mass-preserving.
                const std::uint64_t per_hom =
                    ((1.0 - r) != 0.0 ? 1u : 0u) + ((r / double(NH)) != 0.0 ? NH : 0u);
                std::uint64_t pred_updates = 0;
                for (std::size_t b = 0; b + 1 < NB; ++b)
                    pred_updates += std::uint64_t(keep[b].size()) * per_hom * per_hom;
                ok_(fwd_updates == pred_updates, std::string(arm.name) +
                    ": the forward update count matches its closed-form prediction (" +
                    std::to_string(fwd_updates) + " = " + std::to_string(pred_updates) + ")");
                // PARALLEL SAFETY. The outer loop is over SOURCES and each iteration
                // accumulates only into barf[b][its own source key], reading barf[b+1] read-only.
                // So the adjoint partitions by source with no contention and no reduction buffer.
                // (The FORWARD pass is the contended one -- many sources reach one destination.)
                // The `rev` parameter exists to prove it: visiting sources in the opposite order
                // must give BITWISE identical adjoints, which it cannot if any two sources shared
                // an accumulator.
                std::uint64_t tf_updates = 0, tf_mult_recon = 0, tf_factor_lookups = 0,
                              tf_hist_ops = 0;
                auto sweep = [&](bool rev) {
                std::vector<std::map<std::string, double>> barf(NB);
                for (const auto& kv : keep[NB - 1]) barf[NB - 1][kv.first] = 1.0;
                for (std::size_t b = NB - 1; b-- > 0;) {
                    for (const auto& kv : keep[b]) barf[b][kv.first] = 0.0;
                    std::vector<const std::pair<const std::string,
                                std::pair<St, double>>*> order;
                    for (const auto& e : keep[b]) order.push_back(&e);
                    if (rev) std::reverse(order.begin(), order.end());
                    for (const auto* ep : order) { const auto& kv = *ep;
                        const St& cur = kv.second.first;
                        const std::size_t t1 = cur.x / NH, t2 = cur.x % NH;
                        for (int c1 = 0; c1 < 2; ++c1) for (int c2 = 0; c2 < 2; ++c2) {
                            const double w1 = (c1 == 0) ? (1.0 - r) : (r / double(NH));
                            const double w2 = (c2 == 0) ? (1.0 - r) : (r / double(NH));
                            if (w1 == 0.0 || w2 == 0.0) continue;
                            const std::size_t lo1 = (c1 == 0) ? t1 : 0, hi1 = (c1 == 0) ? t1 + 1 : NH;
                            const std::size_t lo2 = (c2 == 0) ? t2 : 0, hi2 = (c2 == 0) ? t2 + 1 : NH;
                            for (std::size_t n1 = lo1; n1 < hi1; ++n1)
                            for (std::size_t n2 = lo2; n2 < hi2; ++n2) {
                                St s2 = cur;
                                if (c1 == 1) { for (std::size_t q = s2.h1.size(); q < b + 1; ++q)
                                    { s2.h1.push_back(cls(tal[t1][q])); ++tf_hist_ops; } }
                                if (c2 == 1) { for (std::size_t q = s2.h2.size(); q < b + 1; ++q)
                                    { s2.h2.push_back(cls(tal[t2][q])); ++tf_hist_ops; } }
                                s2.x = n1 * NH + n2;
                                double m = w1 * w2 * U[b + 1][s2.x];
                                ++tf_mult_recon;
                                const std::size_t bb = b + 1;
                                if (bb == span1.back() || bb == span2.back()) {
                                    const auto& span = (bb == span1.back()) ? span1 : span2;
                                    auto& M = (bb == span1.back()) ? F1v : F2v;
                                    std::string k;
                                    for (std::size_t q : span) k.push_back(static_cast<char>('0' +
                                        (q < s2.h1.size() ? s2.h1[q] : cls(tal[n1][q]))));
                                    k.push_back('|');
                                    for (std::size_t q : span) k.push_back(static_cast<char>('0' +
                                        (q < s2.h2.size() ? s2.h2[q] : cls(tal[n2][q]))));
                                    const std::size_t half = k.find('|');
                                    const std::string sw = k.substr(half+1) + "|" + k.substr(0,half);
                                    auto it = M.find(std::min(k, sw));
                                    ++tf_factor_lookups;
                                    if (it != M.end()) m *= it->second;
                                }
                                if (bb == span2.back()) { s2.h1.clear(); s2.h2.clear(); }
                                const auto itb = barf[b + 1].find(keyof(s2));
                                if (itb == barf[b + 1].end()) continue;
                                barf[b][kv.first] += itb->second * m;
                                ++tf_updates;
                            }
                        }
                    }
                }
                return barf;
                };
                std::uint64_t z1 = 0, z2 = 0, z3 = 0, z4 = 0;
                const auto barf = sweep(false);
                z1 = tf_updates; z2 = tf_mult_recon; z3 = tf_factor_lookups; z4 = tf_hist_ops;
                tf_updates = tf_mult_recon = tf_factor_lookups = tf_hist_ops = 0;
                const auto barf_rev = sweep(true);
                bool bitwise = (tf_updates == z1);
                for (std::size_t b = 0; b < NB && bitwise; ++b) {
                    if (barf[b].size() != barf_rev[b].size()) { bitwise = false; break; }
                    auto i1 = barf[b].begin(); auto i2 = barf_rev[b].begin();
                    for (; i1 != barf[b].end(); ++i1, ++i2)
                        if (i1->first != i2->first ||
                            std::memcmp(&i1->second, &i2->second, sizeof(double)) != 0)
                            { bitwise = false; break; }
                }
                tf_updates = z1; tf_mult_recon = z2; tf_factor_lookups = z3; tf_hist_ops = z4;
                ok_(bitwise, std::string(arm.name) +
                    ": the adjoint PARTITIONS BY SOURCE -- reversing the source visit order gives "
                    "bitwise identical adjoints, so threads over sources need no reduction and no "
                    "atomic accumulation");
                double worst_tf = 0.0;
                for (std::size_t b = 0; b < NB; ++b) {
                    std::vector<double> mg(NS, 0.0);
                    for (const auto& kv : keep[b]) {
                        const auto it = barf[b].find(kv.first);
                        if (it == barf[b].end()) continue;
                        mg[kv.second.first.x] += kv.second.second * it->second;
                    }
                    for (std::size_t x = 0; x < NS; ++x)
                        worst_tf = std::max(worst_tf, std::abs(mg[x] - marg_bf[b][x]) / Z_bf);
                }
                ok_(tf_updates == fwd_updates && tf_updates == pred_updates,
                    std::string(arm.name) + ": forward and reverse enumerate the SAME loop (" +
                    std::to_string(fwd_updates) + " forward, " + std::to_string(tf_updates) +
                    " adjoint, " + std::to_string(pred_updates) + " predicted)");
                ok_(worst_tf < 1e-12, std::string(arm.name) +
                    ": the TAPE-FREE adjoint gives every marginal (worst relative " +
                    sci(worst_tf) + "); forward " + std::to_string(fwd_updates) + " updates, "
                    "adjoint " + std::to_string(tf_updates) + ", total " +
                    std::to_string(fwd_updates + tf_updates) + "; " +
                    std::to_string(tf_mult_recon) + " multiplier reconstructions, " +
                    std::to_string(tf_factor_lookups) + " factor lookups, " +
                    std::to_string(tf_hist_ops) + " history appends (forward, separately: " +
                    std::to_string(fw_mult) + " / " + std::to_string(fw_lookup) + " / " +
                    std::to_string(fw_hist) + ")");

                // ---- THE TAPED ADJOINT, kept only as a cross-check ---------------------------
                // bar at the last block is 1 for every state: Z is their plain sum. Then every
                // forward update  dst += src * m  contributes  bar_src += bar_dst * m, walked in
                // reverse. This is the ONLY reverse traversal production needs.
                std::vector<std::map<std::string, double>> bar(NB);
                for (const auto& kv : keep[NB - 1]) bar[NB - 1][kv.first] = 1.0;
                for (std::size_t b = NB - 1; b-- > 0;) {
                    for (const auto& kv : keep[b]) bar[b][kv.first] = 0.0;
                    for (const UpdK& u : tape[b]) {
                        const auto it = bar[b + 1].find(u.dst);
                        if (it == bar[b + 1].end()) continue;
                        bar[b][u.src] += it->second * u.m;
                        ++adj_ops;
                    }
                }
                // marginal(b, x) = sum over states at b carrying X = x of msg * bar
                double worst_a = 0.0;
                for (std::size_t b = 0; b < NB; ++b) {
                    std::vector<double> mg(NS, 0.0);
                    for (const auto& kv : keep[b]) {
                        const auto it = bar[b].find(kv.first);
                        if (it == bar[b].end()) continue;
                        mg[kv.second.first.x] += kv.second.second * it->second;
                    }
                    for (std::size_t x = 0; x < NS; ++x)
                        worst_a = std::max(worst_a, std::abs(mg[x] - marg_bf[b][x]) / Z_bf);
                }
                ok_(worst_a < 1e-12, std::string(arm.name) +
                    ": the TAPED adjoint agrees too (worst relative " + sci(worst_a) +
                    ") -- it is a cross-check on the tape-free one, not the production path");
                double worst_tt = 0.0;
                for (std::size_t b = 0; b < NB; ++b)
                    for (const auto& kv : keep[b]) {
                        const auto a1 = bar[b].find(kv.first);
                        const auto a2 = barf[b].find(kv.first);
                        if (a1 == bar[b].end() || a2 == barf[b].end()) continue;
                        worst_tt = std::max(worst_tt, std::abs(a1->second - a2->second));
                    }
                ok_(worst_tt < 1e-12, std::string(arm.name) +
                    ": taped and TAPE-FREE adjoints agree entry for entry (worst " +
                    sci(worst_tt) + ")");
            }

            // EVERY MARGINAL, by clamping each block to each state in turn.
            double worst_c = 0.0;
            for (std::size_t b = 0; b < NB; ++b)
                for (std::size_t x = 0; x < NS; ++x)
                    worst_c = std::max(worst_c,
                                       std::abs(contract(static_cast<int>(b), x) -
                                                marg_bf[b][x]) / Z_bf);
            ok_(worst_c < 1e-12,
                std::string(arm.name) + ": and every block marginal from the contraction, over "
                "all " + std::to_string(NB * NS) + " (worst relative " + sci(worst_c) + ")");
        }

        // NON-VACUITY: the factors must actually move the answer.
        {
            double Z_nf = 0.0;
            std::vector<std::size_t> path(NB, 0);
            std::size_t total = 1;
            for (std::size_t b = 0; b < NB; ++b) total *= NS;
            for (std::size_t code = 0; code < total; ++code) {
                std::size_t c = code;
                for (std::size_t b = 0; b < NB; ++b) { path[b] = c % NS; c /= NS; }
                double w = 1.0;
                for (std::size_t b = 0; b < NB; ++b) w *= U[b][path[b]];
                for (std::size_t b = 1; b < NB; ++b) w *= Tdip(path[b - 1], path[b]);
                Z_nf += w;
            }
            ok_(std::abs(Z_bf - Z_nf) / Z_nf > 0.05,
                std::string(arm.name) +
                ": the higher-order factors MATERIALLY change the partition weight (" + sci(Z_bf) +
                " with, " + sci(Z_nf) + " without)");
        }
        }   // arms

        // ==========================================================================================
        // THROUGH THE PRODUCTION FUNCTION, with REAL factor tables.
        //
        // Everything above is a fixture reimplementation. It proves the mathematics; it does not
        // prove that hybrid_higher_order does the same thing, and a gate that only exercises the
        // fixture would let the production recurrence drift away from it silently. So this arm
        // builds two OVERLAPPING IntervalFactorTables through build_interval_factor, hands them to
        // hybrid_higher_order, and requires agreement with hybrid_bruteforce -- which reads the
        // factors straight off each path and carries no message, no history and no refinement.
        {
            std::mt19937_64 rng2(20260910);
            const auto rseq = [&](std::size_t n) {
                static const char* B = "ACGT";
                std::string t(n, 'A');
                for (std::size_t i = 0; i < n; ++i) t[i] = B[rng2() & 3];
                return t;
            };
            const auto sub = [](std::string x, std::size_t at, char ch) {
                x[at] = (x[at] == ch) ? (ch == 'A' ? 'C' : 'A') : ch; return x;
            };
            // Build one real factor. `ident[j][a]` says WHICH SEQUENCE allele a of block j carries:
            // two alleles sharing an id are byte-identical, so they collapse into one signature
            // class while remaining distinct HMM states. That is the lever that lets two factors
            // induce DIFFERENT class partitions on a shared block -- without it both partitions are
            // the discrete one, either factor's history would do, and the common refinement the
            // recurrence carries would never be tested.
            const auto make_factor = [&](const std::vector<std::vector<int>>& ident,
                                         IntervalGeometry& G) {
                G.blocks.clear();
                for (std::size_t j = 0; j < ident.size(); ++j)
                    G.blocks.push_back(static_cast<std::uint32_t>(j));
                G.alleles.clear();
                for (std::size_t j = 0; j < ident.size(); ++j) {
                    const std::string base = rseq(200);
                    int nseq = 0;
                    for (int id : ident[j]) nseq = std::max(nseq, id + 1);
                    std::vector<std::string> distinct{base};
                    for (int q = 1; q < nseq; ++q)
                        distinct.push_back(sub(base, 40 + 7 * static_cast<std::size_t>(q),
                                               "CGTA"[q & 3]));
                    std::vector<std::string> alt;
                    for (int id : ident[j]) alt.push_back(distinct[static_cast<std::size_t>(id)]);
                    G.alleles.push_back(alt);
                }
                G.contexts.assign(ident.size() - 1, rseq(9));
                G.lflank.clear(); G.rflank.clear();
                G.ok = true; G.exposure_affine = true;
                InsertPrior ip; ip.lo = 200; ip.hi = 700;
                ip.logp.assign(static_cast<std::size_t>(ip.hi - ip.lo + 1),
                               -std::log(static_cast<double>(ip.hi - ip.lo + 1)));
                const double lep = std::log(0.001 / 3.0), l1m = std::log(1.0 - 0.001);
                const auto ix = build_interval_seed_index(G, 16);
                std::vector<Fragment> frags;
                std::vector<std::uint32_t> ch;
                for (std::size_t cell = 0; cell < G.cells(); cell += 3) {
                    G.cell_choice(cell, ch);
                    std::vector<const std::string*> al2, cx2;
                    for (std::size_t j = 0; j < ident.size(); ++j)
                        al2.push_back(&G.alleles[j][ch[j]]);
                    for (const std::string& cx : G.contexts) cx2.push_back(&cx);
                    VirtualWindow vw; vw.bind_chain(G.lflank, al2, cx2, G.rflank);
                    for (std::size_t stt : {std::size_t(30), std::size_t(210)}) {
                        if (stt + 400 > vw.size()) continue;
                        Fragment f;
                        f.name = "g" + std::to_string(cell) + "_" + std::to_string(stt);
                        for (std::size_t i = 0; i < 150; ++i) f.r1.push_back(vw.base_at(stt + i));
                        std::string t2;
                        for (std::size_t i = 0; i < 150; ++i) t2.push_back(vw.base_at(stt + 250 + i));
                        f.r2 = reverse_complement(t2);
                        frags.push_back(f);
                    }
                }
                std::vector<IntervalEmission> ems;
                std::vector<std::vector<std::string>> sigs;
                for (const Fragment& f : frags) {
                    IntervalEmission E = interval_emission(f, G, ip, 0.05, lep, l1m, -420.0, &ix,
                                                           nullptr, true, false);
                    sigs.push_back(E.cell_signature);
                    ems.push_back(std::move(E));
                }
                const IntervalGrouping GR = build_interval_grouping(G, sigs);
                return build_interval_factor(G, GR, ems, 0.05, std::log1p(-0.05), std::log(0.05));
            };

            const std::size_t NHp = 4, NBp = 5, NSp = NHp * NHp;
            // Allele counts per chain block. Blocks 2 and 3 sit in BOTH factors, so their history
            // must carry the common refinement of two different class partitions.
            const std::vector<std::size_t> nA = {2, 3, 3, 2, 3};
            IntervalGeometry GA, GB;
            // Chain block 2 is shared, and the two factors partition its three alleles
            // INCOMPATIBLY: FA merges alleles 0 and 1, FB merges alleles 1 and 2. Neither
            // partition refines the other, so the history must carry their common refinement --
            // all three alleles apart -- and storing either factor's own classes loses the
            // other's information.
            const IntervalFactorTable FA = make_factor({{0, 1, 2}, {0, 0, 1}, {0, 1}}, GA);
            const IntervalFactorTable FB = make_factor({{0, 1, 1}, {0, 1}, {0, 1, 2}}, GB);
            ok_(FA.ok && FB.ok, std::string("production gate: both real factors build (") +
                std::to_string(FA.classes_stored) + " and " +
                std::to_string(FB.classes_stored) + " classes)");
            // MEASURED, not asserted: the common refinement of the two partitions at the shared
            // block must be STRICTLY FINER than either. Equal class counts would not show this --
            // two different 2-class partitions and two identical ones both report 2.
            std::size_t refined_here = 0;
            if (FA.allele_class.size() == 3 && FB.allele_class.size() == 3) {
                std::set<std::pair<std::uint32_t, std::uint32_t>> pairs;
                for (std::size_t a = 0; a < nA[2]; ++a)
                    pairs.insert({FA.allele_class[1][a], FB.allele_class[0][a]});
                (void)0;
                refined_here = pairs.size();
            }
            ok_(refined_here > FA.classes_per_block[1] && refined_here > FB.classes_per_block[0],
                "production gate: the common refinement at the shared block is STRICTLY FINER "
                "than either factor's partition (" + std::to_string(refined_here) + " vs FA " +
                std::to_string(FA.classes_per_block[1]) + " and FB " +
                std::to_string(FB.classes_per_block[0]) + "), so neither refines the other");

            HybridChain ch;
            ch.n_hap = NHp; ch.n_blocks = NBp; ch.recomb = 0.2;
            ch.edges.assign(NBp, HybridEdge{});
            ch.hap_allele.assign(NHp, std::vector<std::uint32_t>(NBp, 0));
            for (std::size_t t = 0; t < NHp; ++t)
                for (std::size_t b = 0; b < NBp; ++b)
                    ch.hap_allele[t][b] = static_cast<std::uint32_t>((t + b) % nA[b]);
            std::mt19937_64 rng3(20260911);
            std::uniform_real_distribution<double> ud(-3.0, 0.0);
            ch.log_emission.assign(NBp, std::vector<double>(NSp, 0.0));
            for (std::size_t b = 0; b < NBp; ++b)
                for (std::size_t k = 0; k < NSp; ++k) ch.log_emission[b][k] = ud(rng3);
            // A RETAINED PAIRWISE EDGE, carried in the same list. Without one the pairwise
            // branch is present and never taken, and the aggregation over source templates would
            // silently destroy the left endpoint it depends on. Edge 0-1 is OUTSIDE both interval
            // factors, so it also makes block 0 enter the carried history on its own account.
            SparseEdgeLinkage PW;
            PW.active = true;
            PW.n_a = nA[0]; PW.n_b = nA[1];
            PW.allele_a.resize(NHp); PW.allele_b.resize(NHp);
            for (std::size_t t = 0; t < NHp; ++t) {
                PW.allele_a[t] = static_cast<std::uint32_t>(t % nA[0]);
                PW.allele_b[t] = static_cast<std::uint32_t>((t + 1) % nA[1]);
            }
            {
                std::mt19937_64 rp(20260912);
                std::uniform_real_distribution<double> pd(-0.9, 0.9);
                for (std::uint32_t amin = 0; amin < nA[0]; ++amin)
                for (std::uint32_t amax = amin; amax < nA[0]; ++amax)
                for (std::uint32_t bmin = 0; bmin < nA[1]; ++bmin)
                for (std::uint32_t bmax = bmin; bmax < nA[1]; ++bmax) {
                    if (amin == amax || bmin == bmax) continue;   // no phase to express
                    SparsePhaseClass k;
                    k.amin = amin; k.amax = amax; k.bmin = bmin; k.bmax = bmax;
                    const double d = pd(rp);
                    k.straight_m1 = std::expm1(d);
                    k.crossed_m1 = std::expm1(-d);
                    PW.classes.push_back(k);
                }
            }
            ch.higher.push_back(HybridHigherFactor{{0, 1}, nullptr, &PW});
            ch.higher.push_back(HybridHigherFactor{{1, 2, 3}, &FA, nullptr});
            ch.higher.push_back(HybridHigherFactor{{2, 3, 4}, &FB, nullptr});
            ok_(PW.has_corrections(),
                "production gate: the fixture carries a RETAINED PAIRWISE edge with " +
                std::to_string(PW.classes.size()) + " non-neutral phase classes, alongside the "
                "two interval factors");

            HigherOrderStats hs;
            HigherOrderTrace tr;
            const HybridPosterior got = hybrid_higher_order(ch, 0, &hs, &tr);
            ok_(got.ok && !hs.refused,
                std::string("production gate: hybrid_higher_order runs") +
                (hs.refused ? (" -- REFUSED: " + hs.refusal) : ""));
            const HybridPosterior bf2 = hybrid_bruteforce(ch);
            ok_(bf2.ok, "production gate: brute force runs over " +
                std::to_string(NSp) + "^" + std::to_string(NBp) + " paths");
            if (got.ok && bf2.ok) {
                double worst = 0.0, worstZ = std::abs(got.log_partition_unnormalised -
                                                      bf2.log_partition_unnormalised);
                for (std::size_t b = 0; b < NBp; ++b)
                    for (std::size_t x = 0; x < NSp; ++x) {
                        const double a = std::exp(got.log_marginal[b][x]);
                        const double e = std::exp(bf2.log_marginal[b][x]);
                        worst = std::max(worst, std::abs(a - e));
                    }
                ok_(worstZ < 1e-9, "production gate: the PARTITION agrees with brute force (log "
                    "difference " + sci(worstZ) + ")");
                ok_(worst < 1e-9, "production gate: EVERY block marginal agrees with brute force "
                    "(worst absolute " + sci(worst) + ")");
            }
            ok_(hs.forward_updates == hs.adjoint_updates,
                "production gate: forward and adjoint enumerate the same loop (" +
                std::to_string(hs.forward_updates) + " each, total " +
                std::to_string(hs.forward_updates + hs.adjoint_updates) + ")");
            ok_(hs.history_ops > 0 && hs.factor_lookups > 0,
                "production gate: auxiliary work is counted separately -- " +
                std::to_string(hs.multiplier_reconstructions) + " multiplier reconstructions, " +
                std::to_string(hs.factor_lookups) + " factor lookups, " +
                std::to_string(hs.history_ops) + " history appends, " +
                std::to_string(hs.grouping_ops) + " aggregate accumulations; peak message " +
                std::to_string(hs.peak_message_entries) + " entries, predicted payload " +
                std::to_string(hs.predicted_payload_bytes) + " bytes");
            // ---- THE PLAN MUST EQUAL WHAT HAPPENED ------------------------------------------
            // Predicted == actual, counter by counter. A plan that merely BOUNDS the work would
            // let --plan-only report a number nobody can act on; equality is what makes it a
            // contract. The plan is exact because the diploid reachable set is the square of the
            // haploid one, which the planner enumerates outright.
            {
                const HigherOrderPlan pl = plan_higher_order(ch, 4);
                ok_(pl.ok, pl.ok ? "production gate: the resource plan builds"
                                 : ("production gate: the plan REFUSED: " + pl.refusal));
                std::vector<std::string> off;
                const auto eq = [&](const char* nm, std::uint64_t p2, std::uint64_t a) {
                    if (p2 != a) off.push_back(std::string(nm) + " planned " +
                                               std::to_string(p2) + " actual " + std::to_string(a));
                };
                eq("forward_updates", pl.forward_updates, hs.forward_updates);
                eq("adjoint_updates", pl.adjoint_updates, hs.adjoint_updates);
                eq("multiplier_reconstructions", pl.multiplier_reconstructions,
                   hs.multiplier_reconstructions);
                eq("factor_lookups", pl.factor_lookups, hs.factor_lookups);
                eq("history_ops", pl.history_ops, hs.history_ops);
                eq("grouping_ops", pl.grouping_ops, hs.grouping_ops);
                eq("peak_message_entries", pl.peak_message_entries, hs.peak_message_entries);
                eq("total_message_entries", pl.total_message_entries, hs.total_message_entries);
                eq("dense_equivalent", pl.dense_equivalent_updates, hs.dense_equivalent_updates);
                std::string why;
                for (std::size_t i = 0; i < off.size(); ++i) why += (i ? "; " : "") + off[i];
                ok_(hs.all_initial_emissions_finite && hs.initial_states_dropped == 0,
                    "production gate: every initial state has a finite nonzero scaled emission, "
                    "which is the condition under which the plan is EXACT rather than an upper "
                    "bound");
                ok_(off.empty(), off.empty()
                    ? ("production gate: EVERY planned count equals the realised one (" +
                       std::to_string(pl.forward_updates) + " forward updates, " +
                       std::to_string(pl.peak_message_entries) + " peak entries, " +
                       std::to_string(pl.total_bytes) + " total bytes over " +
                       std::to_string(pl.threads) + " threads)")
                    : ("production gate: the plan DISAGREES with what happened -- " + why));
                // And the plan's per-block message sizes must match, not just the peak.
                bool per_block = true;
                for (std::size_t b = 0; b < tr.messages.size(); ++b)
                    if (pl.message_entries[b] != tr.messages[b].size()) per_block = false;
                ok_(per_block, "production gate: the planned message size is right at EVERY "
                    "block, not only at the peak");
                ok_(pl.total_bytes > pl.payload_bytes,
                    "production gate: the byte account separates payload " +
                    std::to_string(pl.payload_bytes) + ", container " +
                    std::to_string(pl.container_bytes) + ", temporaries " +
                    std::to_string(pl.temporary_bytes) + " and " + std::to_string(pl.threads) +
                    " x " + std::to_string(pl.per_thread_bytes) + " per-thread");
            }

            // ---- AND WHAT HAPPENS WHEN AN EMISSION IS NOT FINITE ---------------------------
            // With a -inf emission the initial message drops states the plan counted, so the plan
            // becomes an upper BOUND. Both halves are asserted: the run must say so, and the
            // counts must stay under the plan rather than merely differ from it.
            {
                HybridChain deadch = ch;
                deadch.log_emission[0][0] = -std::numeric_limits<double>::infinity();
                deadch.log_emission[0][5] = -std::numeric_limits<double>::infinity();
                HigherOrderStats hz;
                const HybridPosterior zp = hybrid_higher_order(deadch, 0, &hz);
                const HigherOrderPlan pz = plan_higher_order(deadch, 1);
                ok_(zp.ok && !hz.all_initial_emissions_finite && hz.initial_states_dropped == 2,
                    "production gate: a -inf emission is REPORTED as dropping initial states (" +
                    std::to_string(hz.initial_states_dropped) + " dropped, all-finite=" +
                    std::to_string(hz.all_initial_emissions_finite ? 1 : 0) + ")");
                ok_(pz.ok && hz.forward_updates <= pz.forward_updates &&
                    hz.peak_message_entries <= pz.peak_message_entries &&
                    hz.total_message_entries <= pz.total_message_entries,
                    "production gate: and the plan then BOUNDS the work rather than equalling it "
                    "(" + std::to_string(hz.forward_updates) + " <= " +
                    std::to_string(pz.forward_updates) + " updates, " +
                    std::to_string(hz.peak_message_entries) + " <= " +
                    std::to_string(pz.peak_message_entries) + " peak entries)");
                ok_(hz.forward_updates < pz.forward_updates,
                    "production gate: the bound is STRICT here, so the two contracts are really "
                    "distinguishable and the exact one is not passing by accident");
            }

            // ---- WHAT THE BENCHMARK IS ALLOWED TO CLAIM ------------------------------------
            // A runtime benchmark may use uniform emissions ONLY IF values never change topology.
            // So: run the same chain with uniform emissions and with deterministic nonuniform
            // ones, and require IDENTICAL message sizes and identical operation counters. If that
            // holds, the benchmark measures computation -- not numerical behaviour, and not
            // accuracy, neither of which it is entitled to say anything about.
            {
                HybridChain uni = ch;
                for (auto& row : uni.log_emission) row.assign(NSp, 0.0);
                HigherOrderStats hu;
                HigherOrderTrace tu;
                const HybridPosterior up = hybrid_higher_order(uni, 0, &hu, &tu);
                bool same_shape = up.ok && tu.messages.size() == tr.messages.size();
                for (std::size_t b = 0; b < tr.messages.size() && same_shape; ++b)
                    if (tu.messages[b].size() != tr.messages[b].size() ||
                        tu.adjoints[b].size() != tr.adjoints[b].size()) same_shape = false;
                const bool same_ops =
                    hu.forward_updates == hs.forward_updates &&
                    hu.adjoint_updates == hs.adjoint_updates &&
                    hu.multiplier_reconstructions == hs.multiplier_reconstructions &&
                    hu.factor_lookups == hs.factor_lookups &&
                    hu.history_ops == hs.history_ops &&
                    hu.grouping_ops == hs.grouping_ops &&
                    hu.peak_message_entries == hs.peak_message_entries &&
                    hu.total_message_entries == hs.total_message_entries;
                ok_(same_shape && same_ops,
                    "production gate: UNIFORM and nonuniform emissions give identical message "
                    "sizes and identical operation counters, so a runtime benchmark on uniform "
                    "emissions measures computation and claims nothing about accuracy");
            }

            // ---- AGAINST THE DENSE ORACLE, ENTRY BY ENTRY ----------------------------------
            // Marginals agreeing is necessary and not sufficient: two different messages can give
            // the same marginals. Both implementations share one HoPlan, so their keys mean the
            // same thing and every message and adjoint entry is directly comparable.
            HigherOrderStats hd;
            HigherOrderTrace td;
            const HybridPosterior dn = hybrid_higher_order_dense_oracle(ch, 0, &hd, &td);
            bool same_keys = dn.ok && td.messages.size() == tr.messages.size();
            double worst_msg = 0.0, worst_bar = 0.0;
            std::size_t entries_compared = 0;
            for (std::size_t b = 0; b < tr.messages.size() && same_keys; ++b) {
                if (td.messages[b].size() != tr.messages[b].size() ||
                    td.adjoints[b].size() != tr.adjoints[b].size()) { same_keys = false; break; }
                for (std::size_t i = 0; i < tr.messages[b].size(); ++i) {
                    if (td.messages[b][i].first != tr.messages[b][i].first) {
                        same_keys = false; break;
                    }
                    const double a = tr.messages[b][i].second, e = td.messages[b][i].second;
                    worst_msg = std::max(worst_msg,
                                         std::abs(a - e) / (std::abs(e) > 0.0 ? std::abs(e) : 1.0));
                    ++entries_compared;
                }
                for (std::size_t i = 0; i < tr.adjoints[b].size() && same_keys; ++i) {
                    if (td.adjoints[b][i].first != tr.adjoints[b][i].first) {
                        same_keys = false; break;
                    }
                    const double a = tr.adjoints[b][i].second, e = td.adjoints[b][i].second;
                    worst_bar = std::max(worst_bar,
                                         std::abs(a - e) / (std::abs(e) > 0.0 ? std::abs(e) : 1.0));
                }
            }
            ok_(same_keys, "production gate: the factorised recurrence reaches EXACTLY the same "
                "states as the dense oracle, block by block");
            ok_(same_keys && worst_msg < 1e-12 && worst_bar < 1e-12,
                "production gate: every MESSAGE and ADJOINT entry agrees with the dense oracle (" +
                std::to_string(entries_compared) + " message entries, worst relative " +
                sci(worst_msg) + " and " + sci(worst_bar) + ")");
            ok_(dn.ok && std::abs(dn.log_partition_unnormalised -
                                  got.log_partition_unnormalised) < 1e-9,
                "production gate: and the same partition (log difference " +
                sci(std::abs(dn.log_partition_unnormalised -
                             got.log_partition_unnormalised)) + ")");
            ok_(hd.forward_updates == hd.dense_equivalent_updates,
                "production gate: the DENSE oracle costs exactly its dense-equivalent (" +
                std::to_string(hd.forward_updates) + " = " +
                std::to_string(hd.dense_equivalent_updates) + "), so the ratio below is real work "
                "removed and not an accounting artefact");
            ok_(hs.forward_updates <= hs.dense_equivalent_updates,
                "production gate: the Li-Stephens FACTORISATION never costs more than the dense "
                "fan-out (" + std::to_string(hs.forward_updates) + " updates against " +
                std::to_string(hs.dense_equivalent_updates) + " dense-equivalent)");
            // WITHOUT the factors the same chain must reduce to the legacy pairwise path.
            HybridChain plain = ch;
            plain.higher.clear();
            const HybridPosterior lo = hybrid_forward_backward(plain);
            const HybridPosterior lb = hybrid_bruteforce(plain);
            double worst_leg = 0.0;
            if (lo.ok && lb.ok)
                for (std::size_t b = 0; b < NBp; ++b)
                    for (std::size_t x = 0; x < NSp; ++x)
                        worst_leg = std::max(worst_leg, std::abs(std::exp(lo.log_marginal[b][x]) -
                                                                 std::exp(lb.log_marginal[b][x])));
            ok_(lo.ok && lb.ok && worst_leg < 1e-9,
                "production gate: with NO higher factors the chain still takes the legacy pairwise "
                "path and matches brute force (worst " + sci(worst_leg) + ")");
            // AND THE RECURRENCE ITSELF MUST DEGENERATE. With no factors there is no history and
            // no refinement, so hybrid_higher_order reduces to plain Li-Stephens -- if it does not,
            // the history machinery is contributing something it should not.
            HigherOrderStats hs0;
            const HybridPosterior deg = hybrid_higher_order(plain, 0, &hs0);
            double worst_deg = 0.0;
            if (deg.ok && lb.ok)
                for (std::size_t b = 0; b < NBp; ++b)
                    for (std::size_t x = 0; x < NSp; ++x)
                        worst_deg = std::max(worst_deg, std::abs(std::exp(deg.log_marginal[b][x]) -
                                                                 std::exp(lb.log_marginal[b][x])));
            ok_(deg.ok && worst_deg < 1e-9 && hs0.history_ops == 0 && hs0.factor_lookups == 0,
                "production gate: the recurrence DEGENERATES to plain Li-Stephens with no factors "
                "(worst " + sci(worst_deg) + ", " + std::to_string(hs0.history_ops) +
                " history appends, " + std::to_string(hs0.factor_lookups) + " factor lookups)");
        }

        std::printf("higher-order inference selftest: %zu failure(s)\n", fails);
        return fails == 0 ? 0 : 1;
    }

    // =============================================================================================
    // THE CONTRIBUTION AUDIT, on fixtures that separate the four ways it can be wrong.
    //
    // What a fragment's omitted tail can do is NOT its raw mass. The caller evaluates
    //     log( (1-eta) * lambda * (M_a + M_b) + eta * P_bg )
    // so the tail matters only RELATIVE to eta * P_bg. A real C4 run has P_bg at e^-240 and tails
    // at e^-55, which is 186 nats per fragment -- and dismissing those tails on their raw size is
    // exactly the error these fixtures exist to catch.
    if (contribution_selftest) {
        std::size_t fails = 0;
        const auto ok_ = [&](bool c, const std::string& what) {
            std::printf("%s\t%s\n", c ? "ok" : "FAIL", what.c_str());
            if (!c) ++fails;
        };
        const auto sci = [](double x) {
            char b[32]; std::snprintf(b, sizeof b, "%.3e", x); return std::string(b);
        };
        const double eta = 0.05;
        const double lmix = std::log1p(-eta), lbgw = std::log(eta), llam = std::log(0.05);
        const double ninf = -std::numeric_limits<double>::infinity();
        // The width of the contribution interval for a fragment with NO in-band mass and a
        // certified tail bound `omitted`, against a background floor `log_p_bg`.
        const auto width = [&](double omitted, double log_p_bg) {
            const MassInterval m{ninf, omitted};
            const MassInterval c = fragment_contribution(m, m, false, lmix, llam, lbgw, log_p_bg);
            return c.upper - c.lower;
        };

        // ---- 1. ZERO IN-BAND STATES, BUT THE TAIL IS CANDIDATE-SPECIFIC -------------------------
        // states = 0 does not make a fragment neutral. Two candidates with DIFFERENT tail bounds
        // give different contributions, so the fragment can still order them.
        {
            const double bg = -240.0;
            const double wa = width(-55.0, bg), wb = width(-70.0, bg);
            ok_(std::abs(wa - wb) > 1.0,
                "zero in-band states still ORDER two candidates whose tails differ: widths " +
                sci(wa) + " and " + sci(wb) + " nats apart by " + sci(std::abs(wa - wb)));
        }

        // ---- 2. TINY RAW MASS, LARGE CONTRIBUTION ----------------------------------------------
        // The C4 case. e^-55 is negligible against 1e-6 and catastrophic against e^-240.
        {
            const double w_low_bg = width(-55.0, -240.0);
            const double w_high_bg = width(-55.0, -40.0);
            ok_(w_low_bg > 100.0 && w_high_bg < 1e-3,
                "the SAME raw tail of e^-55 is worth " + sci(w_low_bg) +
                " nats against a background of e^-240 and " + sci(w_high_bg) +
                " against e^-40 -- so the raw mass alone decides nothing");
        }

        // ---- 3. INDIVIDUALLY SMALL WIDTHS WHOSE SUM FAILS --------------------------------------
        // The reason the guarantee is on the AGGREGATE. Each of these would pass a per-fragment
        // test at 1e-3 nats; 3,740 of them do not.
        // The width rises monotonically with the tail bound, so a target width is reached by
        // bisecting on it. Scanning and hoping to land in the band is how the first version of
        // this fixture accidentally chose a width that passed BOTH tests and proved nothing.
        const auto omitted_for_width = [&](double target, double bg) {
            double lo = -600.0, hi = 0.0;
            for (int it = 0; it < 200; ++it) {
                const double mid = 0.5 * (lo + hi);
                if (width(mid, bg) < target) lo = mid; else hi = mid;
            }
            return 0.5 * (lo + hi);
        };
        {
            const double tol = 1e-3;
            const std::size_t n = 3740;
            // Individually a hundredth of the tolerance; together thirty-seven times it.
            const double target = tol / 100.0;
            const double one = width(omitted_for_width(target, -240.0), -240.0);
            const double summed = one * static_cast<double>(n);
            ok_(one < tol && summed > tol,
                "each of " + std::to_string(n) + " fragments is individually negligible (" +
                sci(one) + " nats < " + sci(tol) + ") yet their SUM is not (" + sci(summed) +
                ") -- a per-fragment threshold would have passed this");
        }

        // ---- 4. A GENUINELY NEGLIGIBLE AGGREGATE PASSES ----------------------------------------
        // The gate must be able to say yes, or it is not a test -- and it must do so with a
        // NONZERO width, or the pass is just underflow.
        {
            const double tol = 1e-3;
            const std::size_t n = 3740;
            const double target = tol / (10.0 * static_cast<double>(n));
            const double one = width(omitted_for_width(target, -240.0), -240.0);
            const double summed = one * static_cast<double>(n);
            ok_(one > 0.0 && summed <= tol,
                "and a genuinely negligible set PASSES with nonzero widths: " +
                std::to_string(n) + " x " + sci(one) + " = " + sci(summed) + " nats, within " +
                sci(tol));
        }

        // ---- 5. THE BOUND IS CONSERVATIVE IN THE DIRECTION CLAIMED -----------------------------
        // The audit charges both homologues the all-candidate aggregate. That must never
        // UNDERstate a per-candidate interval.
        {
            const double bg = -240.0;
            const MassInterval agg{ninf, -55.0};
            const MassInterval per{ninf, -58.0};   // one candidate's share of the same tail
            const double w_agg = fragment_contribution(agg, agg, false, lmix, llam, lbgw, bg).upper -
                                 fragment_contribution(agg, agg, false, lmix, llam, lbgw, bg).lower;
            const double w_per = fragment_contribution(per, per, false, lmix, llam, lbgw, bg).upper -
                                 fragment_contribution(per, per, false, lmix, llam, lbgw, bg).lower;
            ok_(w_agg >= w_per,
                "the all-candidate aggregate bound never understates a per-candidate one (" +
                sci(w_agg) + " >= " + sci(w_per) + "), so it is an UPPER bound on the width");
        }

        // ---- 6. OVERLAPPING MATES ARE A VALID FRAGMENT ----------------------------------------
        // The insert floor used to be r1 + r2, on the reasoning that "an insert shorter than the
        // two mates is not a state at all". That is only true if mates cannot OVERLAP. A
        // 350 +- 50 library with 150 bp mates puts 15.9% of fragments under 300 bp, and on C4
        // exactly that share -- 3,740 of 23,953 -- had no in-band placement at all. The true floor
        // is max(r1, r2): an insert cannot be shorter than its longest mate.
        {
            std::mt19937_64 rq(20260913);
            std::string ref(3000, 'A');
            static const char* B = "ACGT";
            for (char& c : ref) c = B[rq() & 3];
            std::vector<BlockAlleles> blk;   // one block spanning the whole reference
            CandidateFrame fr;
            fr.seq = ref; fr.ok = true; fr.partial = false;
            fr.offsets = {0}; fr.block_at = {0};
            fr.mapped_lo = 0; fr.mapped_hi = ref.size();
            std::vector<CandidateFrame> frames{fr};
            std::vector<char> block_var{1};

            // An OVERLAPPING pair: r1 = [400,550), r2 = revcomp of [500,650). The physical
            // fragment is 250 bp, shorter than r1 + r2 = 300.
            Fragment f;
            f.name = "overlap";
            f.r1 = ref.substr(400, 150);
            f.r2 = reverse_complement(ref.substr(500, 150));
            const long true_insert = 250;

            const double lep = std::log(0.001 / 3.0), l1m = std::log1p(-0.001);
            const auto owner_with_floor = [&](long floor_len) {
                const InsertPrior ip = make_insert_prior(350.0, 50.0, 0.0, 4, floor_len);
                return assign_fragment_owner_panel_domain(f, frames, block_var, ip, 0.05, lep, l1m, 1e-6,
                                             nullptr);
            };
            const FragmentOwner sum_floor = owner_with_floor(
                static_cast<long>(f.r1.size() + f.r2.size()));          // 300 -- excludes it
            const FragmentOwner max_floor = owner_with_floor(
                static_cast<long>(std::max(f.r1.size(), f.r2.size()))); // 150 -- admits it

            ok_(sum_floor.origins == 0 &&
                sum_floor.why == UnusableReason::NoInBandOrigins,
                "an overlapping pair (insert " + std::to_string(true_insert) +
                " bp, mates 150+150) has NO in-band origin under the r1+r2 floor -- which is the "
                "3,740-fragment C4 signature, reproduced in one fixture");
            ok_(max_floor.origins > 0 && max_floor.kind != OwnerKind::Unusable,
                "and it places normally under the max(r1,r2) floor: " +
                std::to_string(max_floor.origins) + " origin(s), kind " +
                owner_kind_name(max_floor.kind));
            // The mates must place INDIVIDUALLY either way -- it is the pair join, not the search
            // depth, that rejects them. Adaptive deepening could never have repaired this.
            ok_(sum_floor.omitted_bound > -1e300,
                "the mates place individually under both floors (out-of-band bound " +
                sci(sum_floor.omitted_bound) + "), so the rejection was the INSERT SUPPORT and "
                "not the edit band -- deepening would not have found them");
        }

        // ---- 7. BOTH INSERT BOUNDS, AND WHAT RESTORING EITHER ONE COSTS ------------------------
        // The two tails of one distribution. A 350 +- 50 library with 150 bp mates puts 15.9% of
        // fragments below a 300 floor and about 3.2e-5 above a 550 ceiling; on C4 that was 3,746
        // and 1 fragment respectively. Neither is an edit-band case: the C4 survivor matched both
        // mates with ZERO mismatches at 573 bp.
        {
            std::mt19937_64 rq(20260914);
            std::string ref(4000, 'A');
            static const char* B = "ACGT";
            for (char& c : ref) c = B[rq() & 3];
            CandidateFrame fr;
            fr.seq = ref; fr.ok = true; fr.partial = false;
            fr.offsets = {0}; fr.block_at = {0};
            fr.mapped_lo = 0; fr.mapped_hi = ref.size();
            std::vector<CandidateFrame> frames{fr};
            std::vector<char> block_var{1};
            const double lep = std::log(0.001 / 3.0), l1m = std::log1p(-0.001);

            // A pair whose insert is 573 bp -- 4.46 sd out, zero mismatches, correctly oriented.
            Fragment lng;
            lng.name = "long_insert_573";
            lng.r1 = ref.substr(1000, 150);
            lng.r2 = reverse_complement(ref.substr(1000 + 573 - 150, 150));
            const auto owner_at = [&](int sigmas, long floor_len) {
                const InsertPrior ip = make_insert_prior(350.0, 50.0, 0.0, sigmas, floor_len);
                return assign_fragment_owner_panel_domain(lng, frames, block_var, ip, 0.05, lep, l1m, 1e-6,
                                             nullptr);
            };
            const FragmentOwner at4 = owner_at(4, 150);   // support 150-550: rejects it
            const FragmentOwner at6 = owner_at(6, 150);   // support 150-650: admits it
            ok_(at4.origins == 0,
                "restoring the FOUR-sigma ceiling (550) loses a 573 bp pair with zero mismatches "
                "-- the C4 survivor's exact signature");
            ok_(at6.origins > 0 && at6.kind != OwnerKind::Unusable,
                "and the declared SIX-sigma support (650) recovers it: " +
                std::to_string(at6.origins) + " origin(s), kind " + owner_kind_name(at6.kind));

            // And the lower bound, restored, must still fail on the overlapping pair.
            Fragment ov;
            ov.name = "overlap_250";
            ov.r1 = ref.substr(2000, 150);
            ov.r2 = reverse_complement(ref.substr(2100, 150));
            const InsertPrior ip300 = make_insert_prior(350.0, 50.0, 0.0, 6, 300);
            const InsertPrior ip150 = make_insert_prior(350.0, 50.0, 0.0, 6, 150);
            const FragmentOwner lo_bad =
                assign_fragment_owner_panel_domain(ov, frames, block_var, ip300, 0.05, lep, l1m, 1e-6, nullptr);
            const FragmentOwner lo_ok =
                assign_fragment_owner_panel_domain(ov, frames, block_var, ip150, 0.05, lep, l1m, 1e-6, nullptr);
            ok_(lo_bad.origins == 0 && lo_ok.origins > 0,
                "restoring the |r1|+|r2| floor (300) loses the overlapping pair that max(|r1|,"
                "|r2|) (150) keeps -- both bounds are gated, not just the one just changed");

            // THE RESIDUAL IS MEASURED, NOT DECLARED ZERO.
            ok_(std::isfinite(ip150.log_residual_above) &&
                ip150.log_residual_above < std::log(1e-8),
                "the six-sigma support reports the tail it does NOT cover rather than calling it "
                "zero: log residual above = " + sci(ip150.log_residual_above) +
                " (below = " + sci(ip150.log_residual_below) + ")");
            const InsertPrior ip4 = make_insert_prior(350.0, 50.0, 0.0, 4, 150);
            ok_(ip4.log_residual_above > ip150.log_residual_above,
                "and a narrower support reports MORE residual (" + sci(ip4.log_residual_above) +
                " at four sigma against " + sci(ip150.log_residual_above) + " at six), so the "
                "number tracks the support rather than being a constant");
        }

        // ---- 8. IS AN ENLARGED BLOCK CONTEXT, OR A DEPENDENCY? --------------------------------
        // A span enlarged for EXPOSURE is supposed to add normalisation context, not evidence.
        // That claim is testable: build a factor over four blocks from fragments that span only
        // the last three, then vary the FIRST block's allele while holding the others fixed. If
        // the factor's value moves, the enlargement made that block a real dependency and the
        // statistical model changed -- which matters because a per-donor planner may enlarge
        // differently for different donors.
        {
            std::mt19937_64 rq(20260915);
            const auto rseq = [&](std::size_t n) {
                static const char* B = "ACGT";
                std::string t(n, 'A');
                for (std::size_t i = 0; i < n; ++i) t[i] = B[rq() & 3];
                return t;
            };
            const auto sub = [](std::string x, std::size_t at, char c) {
                x[at] = (x[at] == c) ? (c == 'A' ? 'C' : 'A') : c; return x;
            };
            IntervalGeometry G;
            G.blocks = {0, 1, 2, 3};
            const std::string b0 = rseq(300), b1 = rseq(300), b2 = rseq(300), b3 = rseq(300);
            G.alleles = {{b0, sub(b0, 40, 'G')},        // the ENLARGED block: 2 alleles
                         {b1, sub(b1, 50, 'G')},
                         {b2, sub(b2, 60, 'T')},
                         {b3, sub(b3, 70, 'C')}};
            G.contexts = {rseq(8), rseq(8), rseq(8)};
            G.lflank.clear(); G.rflank.clear(); G.ok = true; G.exposure_affine = true;
            InsertPrior ip; ip.lo = 200; ip.hi = 700;
            ip.logp.assign(static_cast<std::size_t>(ip.hi - ip.lo + 1),
                           -std::log(static_cast<double>(ip.hi - ip.lo + 1)));
            const double lep = std::log(0.001 / 3.0), l1m = std::log1p(-0.001);
            const auto ix = build_interval_seed_index(G, 16);
            // FRAGMENTS THAT NEVER TOUCH BLOCK 0: they start after it.
            std::vector<Fragment> frags;
            std::vector<std::uint32_t> ch;
            for (std::size_t cell = 0; cell < G.cells(); cell += 3) {
                G.cell_choice(cell, ch);
                std::vector<const std::string*> al2, cx2;
                for (std::size_t j = 0; j < 4; ++j) al2.push_back(&G.alleles[j][ch[j]]);
                for (const std::string& cx : G.contexts) cx2.push_back(&cx);
                VirtualWindow vw; vw.bind_chain(G.lflank, al2, cx2, G.rflank);
                const std::size_t after_b0 = 320;
                for (std::size_t st : {after_b0, after_b0 + 200}) {
                    if (st + 400 > vw.size()) continue;
                    Fragment f;
                    f.name = "e" + std::to_string(cell) + "_" + std::to_string(st);
                    for (std::size_t i = 0; i < 150; ++i) f.r1.push_back(vw.base_at(st + i));
                    std::string t2;
                    for (std::size_t i = 0; i < 150; ++i) t2.push_back(vw.base_at(st + 250 + i));
                    f.r2 = reverse_complement(t2);
                    frags.push_back(f);
                }
            }
            std::vector<IntervalEmission> ems;
            std::vector<std::vector<std::string>> sigs;
            for (const Fragment& f : frags) {
                IntervalEmission E = interval_emission(f, G, ip, 0.05, lep, l1m, -420.0, &ix,
                                                       nullptr, true, false);
                sigs.push_back(E.cell_signature);
                ems.push_back(std::move(E));
            }
            const IntervalGrouping GR = build_interval_grouping(G, sigs);
            const IntervalFactorTable T =
                build_interval_factor(G, GR, ems, 0.05, std::log1p(-0.05), std::log(0.05));
            ok_(T.ok, "enlargement gate: the four-block factor builds");
            // Vary block 0 only; hold blocks 1..3 fixed on both homologues.
            double worst = 0.0; std::size_t compared = 0;
            for (std::uint32_t a1 = 0; a1 < 2; ++a1)
            for (std::uint32_t a2 = 0; a2 < 2; ++a2)
            for (std::uint32_t a3 = 0; a3 < 2; ++a3)
            for (std::uint32_t b1 = 0; b1 < 2; ++b1)
            for (std::uint32_t b2 = 0; b2 < 2; ++b2)
            for (std::uint32_t b3 = 0; b3 < 2; ++b3) {
                const double ref = T.log_psi({0, a1, a2, a3}, {0, b1, b2, b3});
                for (std::uint32_t x = 0; x < 2; ++x)
                for (std::uint32_t y = 0; y < 2; ++y) {
                    const double v = T.log_psi({x, a1, a2, a3}, {y, b1, b2, b3});
                    worst = std::max(worst, std::abs(v - ref)); ++compared;
                }
            }
            // AND THE MECHANISM, not just the observation. Fragments that never touch the
            // enlarged block have identical signatures across its alleles, so signature grouping
            // collapses it to ONE class. With one class the two homologues always agree there, it
            // never counts toward effective_m, and the mean-one centring log(2^(m-1)) is
            // unchanged -- which is precisely why enlarging for exposure cannot alter the
            // statistical model. If that block ever carried more than one class, the enlargement
            // WOULD have added a phase dimension the evidence cannot inform.
            ok_(!GR.classes_per_block.empty() && GR.classes_per_block[0] == 1,
                "enlargement gate: the enlarged block collapses to ONE signature class (" +
                std::to_string(GR.classes_per_block.empty() ? 0 : GR.classes_per_block[0]) +
                "), so it adds no phase dimension and no change of centring");
            ok_(worst < 1e-12,
                std::string("enlargement gate: block 0 is CONTEXT, not a dependency -- varying it "
                "over ") + std::to_string(compared) + " configurations moves the factor by " +
                sci(worst) + " nats");
        }

        std::printf("contribution selftest: %zu failure(s)\n", fails);
        return fails == 0 ? 0 : 1;
    }

    // ---- HIGHER-ORDER FACTOR CONSTRUCTION SELF-TEST ---------------------------------------------
    // The production constructor against an oracle that computes the diploid formula from its OWN
    // arithmetic. It must not call mix(): a defect there would move both sides together and the
    // comparison would certify nothing.
    if (factor_selftest) {
        std::size_t fails = 0;
        const auto ok_ = [&](bool c, const std::string& w) {
            std::printf("%s\t%s\n", c ? "ok" : "FAIL", w.c_str());
            if (!c) ++fails;
        };
        const auto sci = [](double x) {
            char b[32]; std::snprintf(b, sizeof b, "%.3e", x); return std::string(b);
        };
        std::mt19937_64 rng(20260909);
        const auto rseq = [&](std::size_t n) {
            static const char* B = "ACGT";
            std::string t(n, 'A');
            for (std::size_t i = 0; i < n; ++i) t[i] = B[rng() & 3];
            return t;
        };
        const auto sub = [](std::string x, std::size_t at, char c) {
            x[at] = (x[at] == c) ? (c == 'A' ? 'C' : 'A') : c; return x;
        };
        // A COLLAPSING BLOCK: block 0 carries the SAME sequence twice, so its two alleles are one
        // signature class while remaining two distinct HMM states. Effective heterozygosity there
        // is zero even though the genotype is biologically heterozygous, and swapping them must be
        // exactly neutral -- which is only true if compression stayed an evidence lookup.
        const std::string A0 = rseq(200), B0 = rseq(180), C0 = rseq(220);
        IntervalGeometry G;
        G.blocks = {0, 1, 2};
        G.alleles = {{A0, A0},                                   // one class, two states
                     {B0, sub(B0, 90, 'G'), sub(B0, 91, 'T')},
                     {C0, sub(C0, 110, 'G')}};
        G.contexts = {rseq(10), rseq(8)};
        G.lflank = ""; G.rflank = "";
        G.ok = true; G.exposure_affine = true;
        InsertPrior ip; ip.lo = 200; ip.hi = 650;
        ip.logp.assign(static_cast<std::size_t>(ip.hi - ip.lo + 1),
                       -std::log(static_cast<double>(ip.hi - ip.lo + 1)));
        const double lep = std::log(0.001 / 3.0), l1m = std::log(1.0 - 0.001);
        const auto ix = build_interval_seed_index(G, 16);
        ok_(ix.ok && ix.complete, "the factor fixture's seed index is complete");
        // Fragments planted so that different phases are genuinely preferred.
        std::vector<Fragment> frags;
        std::vector<std::uint32_t> ch;
        for (std::size_t c : {std::size_t(2), std::size_t(7), std::size_t(9)}) {
            if (c >= G.cells()) continue;
            G.cell_choice(c, ch);
            std::vector<const std::string*> al2, cx2;
            for (std::size_t j = 0; j < 3; ++j) al2.push_back(&G.alleles[j][ch[j]]);
            for (const std::string& cx : G.contexts) cx2.push_back(&cx);
            VirtualWindow vw; vw.bind_chain(G.lflank, al2, cx2, G.rflank);
            for (std::size_t st : {std::size_t(60), std::size_t(180), std::size_t(300)}) {
                if (st + 400 > vw.size()) continue;
                Fragment f;
                f.name = "f" + std::to_string(c) + "_" + std::to_string(st);
                for (std::size_t i = 0; i < 150; ++i) f.r1.push_back(vw.base_at(st + i));
                std::string t2;
                for (std::size_t i = 0; i < 150; ++i) t2.push_back(vw.base_at(st + 250 + i));
                f.r2 = reverse_complement(t2);
                frags.push_back(f);
            }
        }
        std::vector<IntervalEmission> ems;
        std::vector<std::vector<std::string>> sigs;
        for (const Fragment& f : frags) {
            const double bgf = -420.0;
            IntervalEmission E = interval_emission(f, G, ip, 0.05, lep, l1m, bgf, &ix,
                                                   nullptr, true, false);
            sigs.push_back(E.cell_signature);
            ems.push_back(std::move(E));
        }
        ok_(!frags.empty() && ems.size() == frags.size(),
            "the fixture plants " + std::to_string(frags.size()) + " fragments");
        const IntervalGrouping GR = build_interval_grouping(G, sigs);
        ok_(GR.ok && GR.joint_equality_verified,
            "grouping validates, joint equality over all " + std::to_string(GR.cells_checked) +
            " cells");
        ok_(GR.classes_per_block[0] == 1,
            "the collapsing block's TWO alleles form ONE signature class (" +
            std::to_string(GR.classes_per_block[0]) + ")");
        const double eta = 0.05, lambda = 0.05;
        const double log_mix = std::log1p(-eta), log_bg = std::log(eta);
        const IntervalFactorTable T = build_interval_factor(G, GR, ems, lambda, log_mix, log_bg);
        // A REFUSAL MUST NAME ITSELF. "0 classes stored" is indistinguishable from a legitimately
        // neutral factor, so the reason is reported and the table is required to be empty.
        ok_(T.ok, T.ok ? ("the factor builds: " + std::to_string(T.classes_stored) +
                          " classes stored, " + std::to_string(T.phase_values_stored) +
                          " phase values, " + std::to_string(T.classes_neutral) + " neutral")
                       : ("the factor REFUSED: \"" + T.refusal + "\", table empty=" +
                          std::to_string(T.phase_value.empty() && T.class_offset.empty() ? 1 : 0)));
        ok_(T.canonical_payload_bytes == T.predicted_payload_bytes &&
            T.total_factor_bytes == T.predicted_total_bytes,
            "measured bytes equal the PREDICTION, payload " +
            std::to_string(T.canonical_payload_bytes) + " and total " +
            std::to_string(T.total_factor_bytes) + " -- compared against their matching estimates");

        // ---- THE ORACLE: its own diploid formula, never mix() ------------------------------
        const auto oracle_S = [&](const std::vector<std::uint32_t>& h1,
                                  const std::vector<std::uint32_t>& h2) {
            const std::size_t c1 = G.cell_index(h1), c2 = G.cell_index(h2);
            double acc = 0.0;
            for (const IntervalEmission& E : ems) {
                // log( (1-eta) * lambda * (M_h1 + M_h2) + eta * P_bg ), written out directly.
                const double m1 = E.mass[c1], m2 = E.mass[c2];
                double lin = 0.0;
                if (m1 != -std::numeric_limits<double>::infinity()) lin += std::exp(m1);
                if (m2 != -std::numeric_limits<double>::infinity()) lin += std::exp(m2);
                const double sig = (1.0 - eta) * lambda * lin;
                const double bg = eta * std::exp(E.log_p_bg);
                acc += std::log(sig + bg);
            }
            return acc;
        };
        // Every ordered diploid configuration over the ORIGINAL alleles.
        std::vector<std::uint32_t> h1(3, 0), h2(3, 0);
        double worst = 0.0; std::size_t checked = 0, nonzero = 0;
        for (std::size_t i1 = 0; i1 < G.cells(); ++i1) {
            G.cell_choice(i1, h1);
            for (std::size_t i2 = 0; i2 < G.cells(); ++i2) {
                G.cell_choice(i2, h2);
                // The RAW content class and its ordered configurations.
                std::vector<std::size_t> het;
                for (std::size_t j = 0; j < 3; ++j) if (h1[j] != h2[j]) het.push_back(j);
                const std::size_t mr = het.size();
                std::vector<double> Sall(std::size_t(1) << mr, 0.0);
                std::vector<std::uint32_t> a(3), b(3);
                for (std::size_t bits = 0; bits < Sall.size(); ++bits) {
                    for (std::size_t j = 0; j < 3; ++j) {
                        a[j] = std::min(h1[j], h2[j]); b[j] = std::max(h1[j], h2[j]);
                    }
                    for (std::size_t q = 0; q < mr; ++q)
                        if (bits & (std::size_t(1) << q)) std::swap(a[het[q]], b[het[q]]);
                    Sall[bits] = oracle_S(a, b);
                }
                double mx = -std::numeric_limits<double>::infinity();
                for (double x : Sall) mx = std::max(mx, x);
                double sm = 0.0;
                for (double x : Sall) sm += std::exp(x - mx);
                const double Z = mx + std::log(sm);
                // This configuration's own bit pattern.
                std::size_t bits = 0;
                for (std::size_t q = 0; q < mr; ++q)
                    if (h1[het[q]] > h2[het[q]]) bits |= (std::size_t(1) << q);
                const double psi_oracle =
                    Sall[bits] - Z + std::log(static_cast<double>(std::size_t(1) << mr));
                const double psi_prod = T.log_psi(h1, h2);
                worst = std::max(worst, std::abs(psi_oracle - psi_prod));
                if (std::abs(psi_oracle) > 1e-9) ++nonzero;
                ++checked;
            }
        }
        ok_(worst < 1e-9, "the constructor equals the independent oracle over all " +
                          std::to_string(checked) + " ordered diploid configurations (worst " +
                          sci(worst) + ")");
        ok_(nonzero > 0, "the fixture is NON-VACUOUS: " + std::to_string(nonzero) +
                         " configurations carry a non-zero psi");
        // THE COLLAPSING BLOCK: swapping its two alleles must change nothing, because they are one
        // signature class -- while remaining two distinct states the caller can still tell apart.
        {
            double worst_sw = 0.0; std::size_t tested = 0;
            for (std::size_t i1 = 0; i1 < G.cells(); ++i1) {
                G.cell_choice(i1, h1);
                for (std::size_t i2 = 0; i2 < G.cells(); ++i2) {
                    G.cell_choice(i2, h2);
                    if (h1[0] == h2[0]) continue;
                    std::vector<std::uint32_t> s1 = h1, s2 = h2;
                    std::swap(s1[0], s2[0]);
                    worst_sw = std::max(worst_sw,
                                        std::abs(T.log_psi(h1, h2) - T.log_psi(s1, s2)));
                    ++tested;
                }
            }
            ok_(tested > 0 && worst_sw < 1e-12,
                "swapping the collapsing block's two alleles is EXACTLY neutral over " +
                std::to_string(tested) + " configurations (worst " + sci(worst_sw) +
                ") -- effective heterozygosity is counted over classes, not alleles");
        }
        // GLOBAL-SWAP EQUALITY IS A THEOREM HERE, NOT A RUNTIME CONDITION, and saying so is more
        // honest than a test that cannot fail. An ordered configuration and its global-swap
        // partner map to the SAME unordered cell pair, and the diploid combiner is symmetric in
        // its two homologues, so their raw scores are identical by construction. No corruption of
        // the emissions can break that -- the first version of this check tried, and could not.
        //
        // What CAN break it is the combiner losing its symmetry, so the precondition is asserted
        // directly and the runtime guard is exercised by mutating mix() rather than by a fixture.
        {
            double worst_sym = 0.0;
            const double lm2 = std::log1p(-eta), lb2 = std::log(eta), ll2 = std::log(lambda);
            const auto mix2 = [&](double p, double q, double pb) {
                double sig = -std::numeric_limits<double>::infinity();
                if (p != -std::numeric_limits<double>::infinity()) sig = p;
                if (q != -std::numeric_limits<double>::infinity())
                    sig = (sig == -std::numeric_limits<double>::infinity())
                              ? q : std::max(sig, q) + std::log1p(std::exp(-std::abs(sig - q)));
                if (sig != -std::numeric_limits<double>::infinity()) sig += lm2 + ll2;
                const double bg = lb2 + pb;
                return (sig == -std::numeric_limits<double>::infinity())
                           ? bg : std::max(sig, bg) + std::log1p(std::exp(-std::abs(sig - bg)));
            };
            for (const IntervalEmission& E : ems)
                for (std::size_t c1 = 0; c1 < G.cells(); ++c1)
                    for (std::size_t c2 = 0; c2 < G.cells(); ++c2)
                        worst_sym = std::max(worst_sym,
                                             std::abs(mix2(E.mass[c1], E.mass[c2], E.log_p_bg) -
                                                      mix2(E.mass[c2], E.mass[c1], E.log_p_bg)));
            ok_(worst_sym == 0.0,
                "the diploid combiner is SYMMETRIC in the two homologues (worst " +
                sci(worst_sym) + ") -- which is what makes global-swap equality automatic, so the "
                "runtime guard is tested by mutating the combiner, not by corrupting a fixture");
        }
        ok_(T.swap_partner_failures == 0,
            "and the constructor's own swap check finds no disagreement (" +
            std::to_string(T.swap_partner_failures) + ")");
        ok_(T.worst_bound_slack >= -1e-9,
            "every class respects its own log psi bound (worst slack " +
            sci(T.worst_bound_slack) + ", max log psi " + sci(T.max_log_psi) + ")");
        std::printf("factor selftest: %zu failure(s)\n", fails);
        return fails == 0 ? 0 : 1;
    }

    // ---- HIGHER-ORDER MEAN-ONE NORMALISATION SELF-TEST ------------------------------------------
    // Two representations of the same object, computed independently and required to agree. The
    // ordered form is what the HMM consumes, because its state is ordered; the canonical
    // biological-phase form is the check on it.
    if (normalisation_selftest) {
        std::size_t fails = 0;
        const auto ok_ = [&](bool c, const std::string& w) {
            std::printf("%s\t%s\n", c ? "ok" : "FAIL", w.c_str());
            if (!c) ++fails;
        };
        const auto sci = [](double x) {
            char b[32]; std::snprintf(b, sizeof b, "%.3e", x); return std::string(b);
        };
        std::mt19937_64 rng(20260909);
        const double kTol = 1e-12;
        // The constants themselves, per heterozygosity.
        for (std::size_t m = 0; m <= 4; ++m) {
            const ContentClassNorm N = content_class_norm(m);
            const std::size_t exp_ord = std::size_t(1) << m;
            const std::size_t exp_bio = m == 0 ? 1 : (std::size_t(1) << (m - 1));
            ok_(N.ok && N.ordered_configs == exp_ord && N.biological_phases == exp_bio,
                "m=" + std::to_string(m) + ": " + std::to_string(N.ordered_configs) +
                " ordered configs, " + std::to_string(N.biological_phases) +
                " biological phases, bound log psi <= " + sci(N.max_log_psi_bound));
        }
        // THE TRAP, as its own assertion: the ordered centring constant is log(2^m). Using
        // log(2^(m-1)) there -- the canonical form's constant -- makes the mean psi one HALF.
        {
            const ContentClassNorm N = content_class_norm(3);
            ok_(std::abs(N.centring_ordered - N.centring_canonical - std::log(2.0)) < kTol,
                "the ordered and canonical centring constants differ by exactly log 2 (" +
                sci(N.centring_ordered) + " vs " + sci(N.centring_canonical) + ")");
        }
        // Now the two normalisations, over random raw scores, for every m.
        for (std::size_t m = 0; m <= 4; ++m) {
            const ContentClassNorm N = content_class_norm(m);
            const std::size_t no = N.ordered_configs;
            // A raw score per BIOLOGICAL phase; each ordered configuration inherits its phase's
            // score, which is what "global swap leaves the score unchanged" means.
            std::vector<double> Sb(N.biological_phases);
            for (double& x : Sb) x = std::uniform_real_distribution<double>(-8.0, 8.0)(rng);
            // Ordered index o in [0, 2^m); its biological phase is o with the top bit folded,
            // because flipping every homologue maps o to its complement.
            const auto phase_of = [&](std::size_t o) {
                if (m == 0) return std::size_t(0);
                const std::size_t comp = (~o) & (no - 1);
                return std::min(o, comp) % N.biological_phases;
            };
            std::vector<double> So(no);
            for (std::size_t o = 0; o < no; ++o) So[o] = Sb[phase_of(o)];
            // GLOBAL-SWAP PARTNERS must exist and score identically.
            bool swap_ok = true;
            for (std::size_t o = 0; o < no; ++o) {
                const std::size_t partner = (~o) & (no - 1);
                if (partner >= no) { swap_ok = false; break; }
                if (So[o] != So[partner]) { swap_ok = false; break; }
            }
            ok_(swap_ok, "m=" + std::to_string(m) +
                         ": every ordered configuration has its global-swap partner with the "
                         "identical raw score" +
                         std::string(m == 0 ? " (and is its own partner)" : ""));
            const auto lse = [](const std::vector<double>& v) {
                double mx = -std::numeric_limits<double>::infinity();
                for (double x : v) mx = std::max(mx, x);
                double acc = 0.0;
                for (double x : v) acc += std::exp(x - mx);
                return mx + std::log(acc);
            };
            const double Zo = lse(So), Zb = lse(Sb);
            std::vector<double> psi_o(no), psi_b(N.biological_phases);
            for (std::size_t o = 0; o < no; ++o) psi_o[o] = So[o] - Zo + N.centring_ordered;
            for (std::size_t b = 0; b < N.biological_phases; ++b)
                psi_b[b] = Sb[b] - Zb + N.centring_canonical;
            // (1) MEAN ONE, in the ordered representation.
            double mean_o = 0.0;
            for (double x : psi_o) mean_o += std::exp(x);
            mean_o /= static_cast<double>(no);
            ok_(std::abs(mean_o - 1.0) < kTol,
                "m=" + std::to_string(m) + ": mean psi over ORDERED configurations is one (" +
                sci(mean_o) + ")");
            // (2) THE TWO REPRESENTATIONS AGREE per ordered configuration.
            double worst = 0.0;
            for (std::size_t o = 0; o < no; ++o)
                worst = std::max(worst, std::abs(psi_o[o] - psi_b[phase_of(o)]));
            ok_(worst < kTol,
                "m=" + std::to_string(m) + ": ordered and canonical psi agree per configuration (" +
                sci(worst) + ")");
            // (3) THE PER-CLASS BOUND holds, and is the SAME in both forms.
            double mx = -std::numeric_limits<double>::infinity();
            for (double x : psi_o) mx = std::max(mx, x);
            ok_(mx <= N.max_log_psi_bound + kTol,
                "m=" + std::to_string(m) + ": max log psi " + sci(mx) + " <= bound " +
                sci(N.max_log_psi_bound));
            // (4) m <= 1 IS EXACTLY NEUTRAL: one biological phase, no preference expressible.
            if (m <= 1) {
                double worst_n = 0.0;
                for (double x : psi_o) worst_n = std::max(worst_n, std::abs(x));
                ok_(worst_n < kTol,
                    "m=" + std::to_string(m) + ": a single biological phase gives an EXACTLY "
                    "neutral factor (" + sci(worst_n) + ")");
            }
        }
        // (5) THE BOUND MUST BE TIGHT, not merely valid. A "max <= bound" check passes for any
        //     bound that is too large, so it cannot detect one. Saturate the class -- put all the
        //     mass on ONE biological phase -- and the maximum must EQUAL log(2^(m-1)) exactly.
        for (std::size_t m = 2; m <= 4; ++m) {
            const ContentClassNorm N = content_class_norm(m);
            const std::size_t no = N.ordered_configs;
            const auto phase_of = [&](std::size_t o) {
                const std::size_t comp = (~o) & (no - 1);
                return std::min(o, comp) % N.biological_phases;
            };
            std::vector<double> Sb(N.biological_phases, -400.0);
            Sb[0] = 0.0;                       // one phase carries everything
            std::vector<double> So(no);
            for (std::size_t o = 0; o < no; ++o) So[o] = Sb[phase_of(o)];
            double mxs = -std::numeric_limits<double>::infinity();
            for (double x : So) mxs = std::max(mxs, x);
            double acc = 0.0;
            for (double x : So) acc += std::exp(x - mxs);
            const double Zo = mxs + std::log(acc);
            double top = -std::numeric_limits<double>::infinity();
            for (std::size_t o = 0; o < no; ++o)
                top = std::max(top, So[o] - Zo + N.centring_ordered);
            ok_(std::abs(top - N.max_log_psi_bound) < 1e-9,
                "m=" + std::to_string(m) + ": a SATURATED class ATTAINS the bound exactly (" +
                sci(top) + " vs " + sci(N.max_log_psi_bound) + ") -- so the bound is tight, not "
                "merely valid");
        }

        // (6) ZERO EVIDENCE is neutral: with no fragments every raw score is equal.
        {
            const ContentClassNorm N = content_class_norm(3);
            std::vector<double> So(N.ordered_configs, 0.0);
            double mx = -std::numeric_limits<double>::infinity();
            double Z;
            {
                double acc = 0.0;
                for (double x : So) acc += std::exp(x);
                Z = std::log(acc);
            }
            for (std::size_t o = 0; o < So.size(); ++o)
                mx = std::max(mx, std::abs(So[o] - Z + N.centring_ordered));
            ok_(mx < kTol, "a factor owning NO fragments is exactly neutral (" + sci(mx) + ")");
        }
        std::printf("normalisation selftest: %zu failure(s)\n", fails);
        return fails == 0 ? 0 : 1;
    }

    // ---- K-BLOCK INTERVAL EMISSION SELF-TEST ----------------------------------------------------
    // The symbolic path must equal the exhaustive oracle PER CELL AND PER STATE, not in totals.
    // Totals agree whenever two errors cancel, and the specific error this guards against --
    // one state expanded into the same cell twice by overlapping symbolic coverage -- is exactly
    // the kind that leaves totals intact while corrupting multiplicity.
    if (interval_selftest) {
        std::size_t fails = 0;
        const auto ok_ = [&](bool c, const std::string& w) {
            std::printf("%s\t%s\n", c ? "ok" : "FAIL", w.c_str());
            if (!c) ++fails;
        };
        const auto sci = [](double x) {
            char b[32]; std::snprintf(b, sizeof b, "%.3e", x); return std::string(b);
        };
        std::mt19937_64 rng(20260909);
        const auto rseq = [&](std::size_t n) {
            static const char* B = "ACGT";
            std::string t(n, 'A');
            for (std::size_t i = 0; i < n; ++i) t[i] = B[rng() & 3];
            return t;
        };
        const auto sub = [](std::string x, std::size_t at, char c) {
            x[at] = (x[at] == c) ? (c == 'A' ? 'C' : 'A') : c; return x;
        };
        // THREE variable blocks with substitution alleles, long enough that no allele can sit
        // wholly inside a seed -- the completeness predicate must hold, and is asserted.
        IntervalGeometry G;
        G.blocks = {0, 1, 2};
        const std::string A0 = rseq(180), B0 = rseq(160), C0 = rseq(200);
        // THE MIDDLE BLOCK CARRIES TWO DISTINCT LENGTHS. With substitution-only alleles every
        // window is the same length, so an intermediate block's length can never affect an insert
        // and a mutation that omits it cannot be observed. A 30 bp indel allele is what makes that
        // condition exist.
        const std::string B_long = B0 + rseq(30);
        // A TANDEM allele in the last block: two identical 400 bp copies, so a fragment planted
        // inside one matches BOTH at the same relative offset. That yields two DISTINCT physical
        // origins with identical edits and identical insert -- the only shape in which collapsing
        // origins by their statistics is observable, and the reason a signature multiset alone
        // cannot certify this search.
        const std::string Q = rseq(400);
        G.alleles = {{A0, sub(A0, 90, 'G')},
                     {B0, sub(B0, 80, 'T'), B_long},
                     {C0, sub(C0, 100, 'G'), Q + Q}};
        // CONTEXTS SHORTER THAN A PIECE, so no seed can sit wholly inside one. At 40 and 35
        // against a 16 bp piece the index was incomplete and did not say so; the geometry passed
        // only because some allele piece always happened to match too.
        G.contexts = {rseq(12), rseq(10)};
        G.lflank = ""; G.rflank = "";
        G.ok = true; G.exposure_affine = true;
        const std::size_t NC = G.cells();
        // BOUNDS SET TO REAL STATES, so both comparisons can be observed. Planted inserts are
        // offset + 150, and offsets 199/200/395/396 give 349/350/545/546 -- one below the lower
        // bound, one exactly on it, one exactly on the upper, one above.
        InsertPrior ip; ip.lo = 350; ip.hi = 545;
        ip.logp.assign(static_cast<std::size_t>(ip.hi - ip.lo + 1),
                       -std::log(static_cast<double>(ip.hi - ip.lo + 1)));
        const double lep = std::log(0.001 / 3.0), l1m = std::log(1.0 - 0.001);
        const auto ix = build_interval_seed_index(G, 16);
        ok_(ix.ok && ix.complete,
            "the seed index is COMPLETE: shortest allele " +
            std::to_string(ix.shortest_allele) + " >= piece - 1, longest context " +
            std::to_string(ix.longest_context) + " < piece");
        {
            IntervalGeometry LC = G;
            LC.contexts[0] = rseq(40);   // long enough to hold a whole seed
            const auto lx = build_interval_seed_index(LC, 16);
            ok_(lx.ok && !lx.complete,
                "a context at least as long as a piece makes the index INCOMPLETE (" +
                std::to_string(lx.longest_context) + " >= 16) -- a seed could sit wholly inside it");
        }
        // ---- SEED-START COVERAGE ORACLE ---------------------------------------------------
        // The completeness predicate is a claim about EVERY seed shape, so it is validated as a
        // whole rather than by accumulating one inequality per shape I happened to think of: for
        // every haploid cell and every possible seed start, either the index carries a hit that
        // reproduces that exact window position, or complete is false. Shapes I never enumerated --
        // allele into context, context into allele, allele across a short context into the next
        // allele, and anything touching a flank -- are covered by construction.
        {
            const auto prefix_at = [&](const std::vector<std::uint32_t>& cc, std::size_t b) {
                std::size_t pfx = G.lflank.size();
                for (std::size_t j = 0; j < b; ++j)
                    pfx += G.alleles[j][cc[j]].size() + (j + 1 < G.alleles.size()
                                                             ? G.contexts[j].size() : 0);
                return pfx;
            };
            const auto coverage = [&](const IntervalGeometry& GG, const IntervalSeedIndex& XX,
                                      std::size_t* uncovered, std::size_t* checked) {
                *uncovered = 0; *checked = 0;
                std::vector<std::uint32_t> cc;
                for (std::size_t c = 0; c < GG.cells(); ++c) {
                    GG.cell_choice(c, cc);
                    std::vector<const std::string*> alle, ctxp;
                    for (std::size_t j = 0; j < GG.alleles.size(); ++j)
                        alle.push_back(&GG.alleles[j][cc[j]]);
                    for (const std::string& cx : GG.contexts) ctxp.push_back(&cx);
                    VirtualWindow vw;
                    vw.bind_chain(GG.lflank, alle, ctxp, GG.rflank);
                    const std::string W = vw.materialize();
                    for (std::size_t at = 0; at + XX.piece <= W.size(); ++at) {
                        ++(*checked);
                        std::uint64_t code = 0;
                        bool enc = true;
                        code = 0;
                        for (std::size_t q = 0; q < XX.piece && enc; ++q) {
                            const char ch = W[at + q];
                            const int v = ch == 'A' ? 0 : ch == 'C' ? 1 : ch == 'G' ? 2
                                        : ch == 'T' ? 3 : -1;
                            if (v < 0) enc = false; else code = (code << 2) | std::uint64_t(v);
                        }
                        if (!enc) continue;
                        bool cov = false;
                        for (int pass = 0; pass < 2 && !cov; ++pass) {
                            const auto& mp = pass == 0 ? XX.inside : XX.boundary;
                            const auto it = mp.find(code);
                            if (it == mp.end()) continue;
                            for (const IntervalSeedHit& h : it->second) {
                                if (cc[h.block] != h.allele) continue;
                                if (h.spans_boundary && h.block + 1 < cc.size() &&
                                    cc[h.block + 1] != h.next_allele) continue;
                                if (prefix_at(cc, h.block) + h.offset == at) { cov = true; break; }
                            }
                        }
                        if (!cov) ++(*uncovered);
                    }
                }
            };
            std::size_t unc = 0, chk = 0;
            coverage(G, ix, &unc, &chk);
            ok_(unc == 0 && chk > 0,
                "SEED-START COVERAGE: every one of " + std::to_string(chk) +
                " seed starts over every cell is represented by the index (" +
                std::to_string(unc) + " uncovered)");
            // And the predicate must be HONEST in the other direction: a geometry it calls
            // incomplete really does have uncovered starts.
            IntervalGeometry LG = G;
            LG.contexts[0] = rseq(40);
            const auto lx = build_interval_seed_index(LG, 16);
            std::size_t unc2 = 0, chk2 = 0;
            coverage(LG, lx, &unc2, &chk2);
            ok_(!lx.complete && unc2 > 0,
                "a geometry the predicate calls INCOMPLETE genuinely has " +
                std::to_string(unc2) + " uncovered seed starts -- the predicate is not merely "
                "conservative here");
        }

        // Fragments planted across every segment and boundary of a chosen cell.
        std::vector<Fragment> frags;
        std::vector<std::uint32_t> ch;
        for (std::size_t c : {std::size_t(0), std::size_t(5), std::size_t(11),
                              std::size_t(2), std::size_t(17)}) {
            if (c >= G.cells()) continue;
            G.cell_choice(c, ch);
            const bool cc_has_tandem = ch[2] == 2;
            std::vector<const std::string*> alle, ctxp;
            for (std::size_t j = 0; j < 3; ++j) alle.push_back(&G.alleles[j][ch[j]]);
            for (const std::string& cx : G.contexts) ctxp.push_back(&cx);
            VirtualWindow vw;
            vw.bind_chain(G.lflank, alle, ctxp, G.rflank);
            const std::size_t wl = vw.size();
            // READ LENGTH IS NOT FREE: the band is floor(0.05 * len) + 1 and the piece is
            // len / (band + 1), so 150 bp gives piece 16 and matches the index. At 100 bp the
            // piece is 14, the symbolic path refuses on a piece mismatch, and every comparison is
            // silently skipped -- which is how the first version of this fixture "passed" three
            // assertions while scoring nothing at all.
            for (std::size_t st : {std::size_t(5), std::size_t(120), std::size_t(190),
                                   std::size_t(240), std::size_t(260)}) {
                if (st + 350 > wl) continue;
                Fragment f;
                f.name = "c" + std::to_string(c) + "_" + std::to_string(st);
                for (std::size_t i = 0; i < 150; ++i) f.r1.push_back(vw.base_at(st + i));
                std::string tail;
                for (std::size_t i = 0; i < 150; ++i) tail.push_back(vw.base_at(st + 200 + i));
                f.r2 = reverse_complement(tail);
                frags.push_back(f);
            }
            // A fragment planted inside the tandem unit, which therefore has two origins.
            if (cc_has_tandem) {
                const std::size_t base = G.lflank.size() + G.alleles[0][0].size() +
                                         G.contexts[0].size() + G.alleles[1][0].size() +
                                         G.contexts[1].size();
                Fragment tf;
                tf.name = "tandem" + std::to_string(c);
                for (std::size_t i = 0; i < 150; ++i) tf.r1.push_back(vw.base_at(base + 5 + i));
                std::string tt;
                for (std::size_t i = 0; i < 150; ++i) tt.push_back(vw.base_at(base + 205 + i));
                tf.r2 = reverse_complement(tt);
                frags.push_back(tf);
            }
            // THE FOUR BOUNDARY OBSERVATIONS. Without a state exactly ON each bound, shifting
            // that bound changes nothing and the mutation survives on an absent condition.
            for (std::size_t off : {std::size_t(199), std::size_t(200), std::size_t(395),
                                    std::size_t(396)}) {
                if (5 + off + 150 > wl) continue;
                Fragment fb2;
                fb2.name = "ins" + std::to_string(off + 150) + "_" + std::to_string(c);
                for (std::size_t i = 0; i < 150; ++i) fb2.r1.push_back(vw.base_at(5 + i));
                std::string tb;
                for (std::size_t i = 0; i < 150; ++i) tb.push_back(vw.base_at(5 + off + i));
                fb2.r2 = reverse_complement(tb);
                frags.push_back(fb2);
            }
            // A LONG-INSERT FRAGMENT whose mates sit in the OUTER blocks, leaving the middle one
            // free: the only shape whose insert an intermediate allele length can change.
            if (wl > 560) {
                Fragment g2;
                g2.name = "span" + std::to_string(c);
                for (std::size_t i = 0; i < 150; ++i) g2.r1.push_back(vw.base_at(5 + i));
                std::string t2;
                for (std::size_t i = 0; i < 150; ++i) t2.push_back(vw.base_at(400 + i));
                g2.r2 = reverse_complement(t2);
                frags.push_back(g2);
                // THE SECOND LIBRARY ORIENTATION: r2 forward and r1 reverse-complemented. Without
                // one of these, dropping that orientation changes nothing and the mutation
                // survives on an absent condition rather than on correctness.
                Fragment g3;
                g3.name = "revorient" + std::to_string(c);
                for (std::size_t i = 0; i < 150; ++i) g3.r2.push_back(vw.base_at(5 + i));
                std::string t3;
                for (std::size_t i = 0; i < 150; ++i) t3.push_back(vw.base_at(400 + i));
                g3.r1 = reverse_complement(t3);
                frags.push_back(g3);
            }
        }
        ok_(frags.size() >= 9, "the fixture plants " + std::to_string(frags.size()) +
                               " fragments across segments and boundaries");
        std::size_t cells_differ = 0, mult_differ = 0, sig_differ = 0, org_differ = 0,
                    con_differ = 0;
        std::size_t at_lo = 0, at_hi = 0, below_lo = 0, above_hi = 0;
        std::size_t free_mid = 0, or_a = 0, or_b = 0, dup_origin_cells = 0;
        std::size_t fin_o = 0, fin_s = 0;
        double worst = 0.0;
        std::size_t tot_states_o = 0, tot_states_s = 0;
        std::size_t sym_states = 0, expansions = 0, joins = 0;
        for (const Fragment& f : frags) {
            const double bgf = -400.0;
            const IntervalEmission O = interval_emission_oracle(f, G, ip, 0.05, lep, l1m, bgf,
                                                                true, true);
            const IntervalEmission S = interval_emission(f, G, ip, 0.05, lep, l1m, bgf, &ix,
                                                         nullptr, true, true);
            // A REFUSAL IS A FAILURE HERE, not a skip: skipping is what made the comparison
            // vacuous, because the loop below never ran and every per-cell assertion passed on an
            // empty set.
            if (!O.ok || !S.ok) {
                ++cells_differ;
                std::printf("FAIL\tinterval emission refused on %s: oracle_ok=%d fast_ok=%d %s\n",
                            f.name.c_str(), O.ok ? 1 : 0, S.ok ? 1 : 0, S.refusal.c_str());
                ++fails;
                continue;
            }
            // THESE COUNT DIFFERENT THINGS. The oracle visits every cell, so its
            // verified_fr_states is a count of state-CELL writes; the symbolic path counts JOINED
            // states before expansion, each covering many cells. Comparing them directly is a
            // category error -- the comparable quantity is the total cell-state writes, which is
            // the sum of cell_states on both sides.
            for (std::size_t c = 0; c < NC; ++c) {
                tot_states_o += O.cell_states[c];
                tot_states_s += S.cell_states[c];
            }
            sym_states += S.verified_fr_states;
            free_mid += S.states_with_free_intermediate;
            or_a += S.states_orientation_a;
            or_b += S.states_orientation_b;
            // Cells where two different origins share (edits, insert): the tandem case.
            for (std::size_t c = 0; c < NC; ++c) {
                const std::string& sg = S.cell_signature[c];
                const std::string& ob = S.cell_origin[c];
                const std::size_t ns = sg.size() / 12;
                if (ns < 2 || ob.size() / 40 != ns) continue;
                bool dup = false;
                for (std::size_t x = 0; x + 1 < ns && !dup; ++x)
                    for (std::size_t y = x + 1; y < ns; ++y)
                        if (std::memcmp(&sg[x * 12], &sg[y * 12], 12) == 0 &&
                            std::memcmp(&ob[x * 40], &ob[y * 40], 40) != 0) { dup = true; break; }
                if (dup) ++dup_origin_cells;
            }
            // The four boundary observations, read off the ORIGINS rather than assumed from the
            // planted offsets: a planted fragment only counts if it actually produced a state.
            for (std::size_t c = 0; c < NC; ++c) {
                const std::string& ob = S.cell_origin[c];
                for (std::size_t q = 0; q + 40 <= ob.size(); q += 40) {
                    long ins = 0; std::memcpy(&ins, &ob[q + 32], 8);
                    if (ins == ip.lo) ++at_lo;
                    if (ins == ip.hi) ++at_hi;
                    if (ins < ip.lo) ++below_lo;
                    if (ins > ip.hi) ++above_hi;
                }
            }
            expansions += S.tuple_expansions;
            joins += S.joined_pairs;
            for (std::size_t c = 0; c < NC; ++c) {
                const bool fo = O.mass[c] != -std::numeric_limits<double>::infinity();
                const bool fsx = S.mass[c] != -std::numeric_limits<double>::infinity();
                if (fo) ++fin_o;
                if (fsx) ++fin_s;
                if (fo != fsx) {
                    ++cells_differ;
                    if (cells_differ <= 3) {
                        std::vector<std::uint32_t> dc;
                        G.cell_choice(c, dc);
                        std::printf("....\tdiff %s cell %zu (%u,%u,%u): oracle=%s fast=%s "
                                    "oracle_states=%u fast_states=%u\n",
                                    f.name.c_str(), c, dc[0], dc[1], dc[2],
                                    fo ? "finite" : "-inf", fsx ? "finite" : "-inf",
                                    O.cell_states[c], S.cell_states[c]);
                    }
                    continue;
                }
                if (O.cell_states[c] != S.cell_states[c]) ++mult_differ;
                if (O.cell_signature[c] != S.cell_signature[c]) ++sig_differ;
                if (O.cell_origin[c] != S.cell_origin[c]) ++org_differ;
                if (O.cell_contrib[c].size() != S.cell_contrib[c].size()) ++con_differ;
                else for (std::size_t q = 0; q < O.cell_contrib[c].size(); ++q)
                    if (O.cell_contrib[c][q] != S.cell_contrib[c][q]) { ++con_differ; break; }
                if (fo) worst = std::max(worst, std::abs(O.mass[c] - S.mass[c]));
            }
        }
        ok_(cells_differ == 0, "the same cells are finite in both paths (" +
                               std::to_string(cells_differ) + " differ)");
        ok_(mult_differ == 0, "state MULTIPLICITY is identical per cell (" +
                              std::to_string(mult_differ) + " differ)");
        ok_(sig_differ == 0, "the structural signature MULTISET is identical per cell (" +
                             std::to_string(sig_differ) + " differ)");
        // THE STRONGER COMPARISON: full origin identity -- both mate starts, both strands and the
        // insert. Equal signatures cannot tell one repeat origin from another with the same
        // statistics, so this is what makes the origin-collapse mutation meaningful.
        ok_(org_differ == 0, "the canonical ORIGIN multiset is identical per cell -- both mate "
                             "starts, strands and insert (" + std::to_string(org_differ) +
                             " differ)");
        // MASS CARRIES A TOLERANCE; THE DISCRETE STRUCTURES DO NOT. The origin multisets are
        // proved identical above, so the two paths sum the SAME terms -- but log_add is not
        // associative in floating point and they sum them in different orders. The exact claims
        // therefore sit where exactness is meaningful (finite cells, multiplicity, signature
        // multiset, origin multiset) and mass inherits a declared tolerance from them, rather than
        // the tolerance being chosen to make a disagreement pass.
        ok_(con_differ == 0, "the per-origin log CONTRIBUTION multiset is identical per cell (" +
                             std::to_string(con_differ) + " differ) -- so only reduction order can "
                             "differ below");
        // With the contribution computed in one shared place, the two paths reduce identical
        // term multisets and the aggregate agrees exactly too. The tolerance stays declared rather
        // than removed: reduction order could differ on another fixture without anything being
        // wrong, and the exact per-origin comparison above is what would still hold if it did.
        ok_(worst <= 1e-12, "log mass agrees per cell (worst " + sci(worst) +
                            ", tolerance 1.000e-12, and exact on this fixture)");
        ok_(fin_o > 0 && fin_s == fin_o,
            "the fixture is NON-VACUOUS: " + std::to_string(fin_o) + " finite cells in the oracle");
        ok_(tot_states_o > 0 && tot_states_o == tot_states_s,
            "total cell-state WRITES agree (" + std::to_string(tot_states_o) + ")");
        ok_(sym_states > 0 && sym_states <= tot_states_s,
            "the symbolic path forms " + std::to_string(sym_states) +
            " joined states covering " + std::to_string(tot_states_s) +
            " cell writes -- free dimensions are never enumerated at the seed");
        std::printf("....\t%zu mate joins considered, %zu intervening-length tuples\n",
                    joins, expansions);
        // ---- BASELINE NON-VACUITY OBSERVABLES ---------------------------------------------
        // Every condition a mutation perturbs, counted BEFORE any mutation runs. A mutation that
        // fails an assertion proves nothing unless its condition existed.
        ok_(at_lo > 0 && at_hi > 0,
            "states sit EXACTLY on both insert bounds: " + std::to_string(at_lo) + " at lo=" +
            std::to_string(ip.lo) + ", " + std::to_string(at_hi) + " at hi=" +
            std::to_string(ip.hi));
        ok_(below_lo == 0 && above_hi == 0,
            "and none outside them: " + std::to_string(below_lo) + " below lo, " +
            std::to_string(above_hi) + " above hi (planted 349 and 546 are REJECTED)");
        ok_(free_mid > 0,
            "states whose mates straddle a FREE intermediate block: " + std::to_string(free_mid) +
            " -- the only shape an intermediate allele length can affect");
        ok_(dup_origin_cells > 0,
            "cells carrying TWO DISTINCT origins with identical edits and insert: " +
            std::to_string(dup_origin_cells) + " -- without one, collapsing origins by their "
            "statistics changes nothing and that mutation survives on an absent condition");
        ok_(or_a > 0 && or_b > 0,
            "both library orientations are populated: " + std::to_string(or_a) + " and " +
            std::to_string(or_b));
        // ---- PERMANENT REGRESSION: a seed ending inside a context must not pin the next allele.
        // It produced one physical origin written once; pinning the untouched allele duplicated it
        // once per candidate and showed up as a log 2 mass difference over eight cells.
        {
            std::size_t ctx_only = 0, spans = 0;
            for (const auto& kv : ix.boundary)
                for (const IntervalSeedHit& h : kv.second)
                    if (h.spans_boundary) ++spans; else ++ctx_only;
            ok_(ctx_only > 0 && spans > 0,
                "the boundary index separates seeds that STOP in the context (" +
                std::to_string(ctx_only) + ") from those REACHING the next allele (" +
                std::to_string(spans) + ")");
        }

        // A SHORT ALLELE must trip the completeness predicate, not be scored incompletely.
        {
            IntervalGeometry H = G;
            H.alleles[1].push_back(rseq(8));   // shorter than piece - 1
            const auto hx = build_interval_seed_index(H, 16);
            ok_(hx.ok && !hx.complete,
                "an allele shorter than piece - 1 makes the index INCOMPLETE (shortest " +
                std::to_string(hx.shortest_allele) + ")");
            const IntervalEmission R = interval_emission(frags[0], H, ip, 0.05, lep, l1m, -400.0,
                                                         &hx, nullptr, false);
            ok_(R.work_refused && R.refusal == "interval-seed-index-incomplete",
                "and the emission REFUSES on it rather than scoring incompletely");
        }
        std::printf("interval selftest: %zu failure(s)\n", fails);
        return fails == 0 ? 0 : 1;
    }

    // ---- EXACT SIGNATURE GROUPING SELF-TEST -----------------------------------------------------
    // Equality of distinct pattern SETS is equality of vocabulary, not of the factor application. A
    // set comparison discards how many content classes realise each pattern, which alleles realise
    // it, the SIGN attached to each occurrence, and the weights of those members. So this compares
    // EVERY content class, member-weighted and signed, and then carries the comparison through psi,
    // the forward-backward weight sum and every marginal.
    if (grouping_selftest) {
        std::size_t fails = 0;
        const auto ok_ = [&](bool c, const std::string& what) {
            std::printf("%s\t%s\n", c ? "ok" : "FAIL", what.c_str());
            if (!c) ++fails;
        };
        // RESIDUALS IN SCIENTIFIC NOTATION, against a DECLARED tolerance. "0.000000" cannot
        // distinguish an exact zero from 1e-9, and that difference is the whole claim.
        const auto sci = [](double x) {
            char buf[32];
            std::snprintf(buf, sizeof buf, "%.3e", x);
            return std::string(buf);
        };
        const double kExact = 0.0;        // bit-identical required
        const double kTolMarg = 1e-15;    // marginals, after normalisation
        const double kTolZ = 1e-12;       // partition weight, psi, swap
        std::mt19937_64 rng(20260908);
        const auto rseq = [&](std::size_t n) {
            static const char* B = "ACGT";
            std::string t(n, 'A');
            for (std::size_t i = 0; i < n; ++i) t[i] = B[rng() & 3];
            return t;
        };
        // A FIXTURE THAT ACTUALLY COLLAPSES: three distinct A sequences each carried by two
        // alleles, four distinct B sequences each carried by two. If the row and column classes do
        // not come out at 3 and 4 the fixture is vacuous and every later assertion is meaningless,
        // so that is asserted first.
        // THE ALLELES MUST DIFFER BY SUBSTITUTIONS, not be independent sequences. Independent
        // sequences make every non-matching cell -inf, so a cell's mismatch count never varies
        // while its insert structure stays fixed -- and then dropping mismatch counts from the
        // signature merges nothing and the mutation cannot fail. Real C4 is the substitution case:
        // 119 of 126 B alleles matched a fragment perfectly and 7 cost exactly one mismatch.
        const auto sub = [](std::string x, std::size_t at, char c) {
            x[at] = (x[at] == c) ? (c == 'A' ? 'C' : 'A') : c;
            return x;
        };
        const std::string A0 = rseq(120), B0 = rseq(140);
        std::vector<std::string> A3 = {A0, sub(A0, 40, 'G'), sub(A0, 41, 'T')};
        std::vector<std::string> B4 = {B0, sub(B0, 60, 'G'), sub(B0, 61, 'T'),
                                       sub(sub(B0, 60, 'G'), 61, 'T')};
        LinkageGeometry g;
        g.block_a = 0; g.block_b = 1;
        // UNEQUAL CLASS SIZES on purpose: a grouped implementation that treated members uniformly
        // would pass on classes that all happen to be the same size.
        // AND INTERLEAVED, not contiguous. Laid out in blocks, class(a1) <= class(a2) holds for
        // every a1 < a2, so only ONE orientation of each class pair is ever realisable and a
        // mutation that merges straight with crossed has nothing to collide -- the test would pass
        // while measuring nothing. Real classes are scattered through the allele order, as C4's
        // are, and scattering them here is what makes the orientation gate load-bearing.
        for (std::size_t k : {0u, 1u, 2u, 1u, 2u, 2u}) g.alleles_a.push_back(A3[k]);
        for (std::size_t k : {0u, 1u, 2u, 3u, 1u, 2u, 3u, 3u}) g.alleles_b.push_back(B4[k]);
        g.context = rseq(60); g.lflank = rseq(900); g.rflank = rseq(900);
        const std::size_t NA = g.alleles_a.size(), NB = g.alleles_b.size();
        g.window_len.assign(NA * NB, 0);
        g.exposure.assign(NA * NB, 1000.0);
        g.exposure_affine = true; g.ok = true;
        InsertPrior ip; ip.lo = 200; ip.hi = 600;
        ip.logp.assign(static_cast<std::size_t>(ip.hi - ip.lo + 1),
                       -std::log(static_cast<double>(ip.hi - ip.lo + 1)));
        const auto ix = build_allele_product_index(g, 16);
        const double lep = std::log(0.001 / 3.0), l1m = std::log(1.0 - 0.001);

        // Fragments planted across several windows so the signatures actually differ.
        std::vector<Fragment> frags;
        for (std::size_t al : {std::size_t(0), std::size_t(1), std::size_t(2),
                               std::size_t(4), std::size_t(5)}) {
            for (std::size_t be : {std::size_t(0), std::size_t(3), std::size_t(6)}) {
                VirtualWindow vw;
                vw.bind_pair(g, static_cast<std::uint32_t>(al),
                             static_cast<std::uint32_t>(be));
                const std::size_t start = g.lflank.size() - 30;
                Fragment f;
                f.name = "f" + std::to_string(al) + "_" + std::to_string(be);
                for (std::size_t i = 0; i < 150; ++i) f.r1.push_back(vw.base_at(start + i));
                std::string tail;
                for (std::size_t i = 0; i < 150; ++i) tail.push_back(vw.base_at(start + 260 + i));
                f.r2 = reverse_complement(tail);
                frags.push_back(f);
            }
        }
        std::vector<LinkageEmission> ems;
        std::vector<std::vector<std::string>> sigs;
        for (const Fragment& f : frags) {
            std::vector<std::string> cs;
            const std::size_t len = f.bases();
            const double bgf = static_cast<double>(len) * 0.10 * lep +
                               static_cast<double>(len) * 0.90 * l1m;
            ems.push_back(linkage_emission_supported(f, g, ip, 0.05, lep, l1m, bgf,
                                                     nullptr, &ix, nullptr, &cs));
            sigs.push_back(std::move(cs));
        }
        // The joint signature per cell is the concatenation over fragments, length-prefixed.
        std::vector<std::string> joint(NA * NB);
        for (std::size_t k = 0; k < NA * NB; ++k) {
            for (const auto& cs : sigs) {
                const std::uint32_t n = static_cast<std::uint32_t>(cs[k].size());
                joint[k].append(reinterpret_cast<const char*>(&n), 4);
                joint[k].append(cs[k]);
            }
        }
        const SignatureMatrix M = build_signature_matrix(joint, NA, NB);
        ok_(M.ok, "the signature matrix validates: classes partition and every cell reconstructs");
        ok_(M.row_members.size() == 3 && M.col_members.size() == 3,
            "the fixture COLLAPSES as designed: " + std::to_string(M.row_members.size()) +
            " row classes of " + std::to_string(NA) + ", " +
            std::to_string(M.col_members.size()) + " column classes of " + std::to_string(NB));
        // WHY THREE COLUMN CLASSES FROM FOUR DISTINCT B SEQUENCES, and not four. Two alleles
        // carrying substitutions at DIFFERENT positions are different sequences, yet each sits
        // exactly one mismatch from the observed reads, so both yield the same
        // (m1_edits, m2_edits, insert) multiset and land in one signature class. Signature
        // equivalence is equidistance from THESE fragments, not sequence identity -- which is
        // exactly the C4 phenomenon that makes the collapse large, and is worth asserting rather
        // than discovering again.
        {
            bool distinct_share = false;
            for (const auto& mem : M.col_members)
                for (std::size_t i = 0; i + 1 < mem.size() && !distinct_share; ++i)
                    for (std::size_t j = i + 1; j < mem.size(); ++j)
                        if (g.alleles_b[mem[i]] != g.alleles_b[mem[j]]) { distinct_share = true; break; }
            ok_(distinct_share,
                "two DIFFERENT B sequences share a signature class -- equidistance, not identity");
        }

        // THE FACTOR APPLICATION, per content class. mix() is the sparse builder's own combiner,
        // reproduced here so the two cannot drift; delta is straight minus crossed, SIGNED.
        const double eta = 0.05, lambda = 0.05;
        const double log_mix = std::log1p(-eta), log_bg_weight = std::log(eta);
        const double log_lambda = std::log(lambda);
        const double NEG = -std::numeric_limits<double>::infinity();
        const auto log_add_d = [&](double x, double y) {
            if (x == NEG) return y;
            if (y == NEG) return x;
            const double hi = std::max(x, y), lo = std::min(x, y);
            return hi + std::log1p(std::exp(lo - hi));
        };
        const auto mix = [&](double p, double q, double log_p_bg) {
            double sig = NEG;
            if (p != NEG) sig = p;
            if (q != NEG) sig = (sig == NEG) ? q : log_add_d(sig, q);
            if (sig != NEG) sig += log_mix + log_lambda;
            const double bg = log_bg_weight + log_p_bg;
            return (sig == NEG) ? bg : log_add_d(sig, bg);
        };
        // Delta for one content class, computed from the emissions themselves.
        const auto delta_direct = [&](std::size_t a1, std::size_t a2, std::size_t b1,
                                      std::size_t b2, int sig_mode) {
            double d = 0.0;
            for (const LinkageEmission& m : ems) {
                const double s = mix(m.mass[a1 * NB + b1], m.mass[a2 * NB + b2], m.log_p_bg);
                const double c = mix(m.mass[a1 * NB + b2], m.mass[a2 * NB + b1], m.log_p_bg);
                d += (sig_mode == 2) ? (s + c) : (s - c);   // mode 2 merges straight with crossed
            }
            return d;
        };
        // A 16-BYTE pattern key. The earlier four 16-bit fields packed into 64 bits would silently
        // alias two patterns once a locus carried more than 65535 signatures, and a silent
        // collision in a CHECKER makes the equivalence it tests pass spuriously.
        const auto pattern_key = [&](std::uint32_t q11, std::uint32_t q22, std::uint32_t q12,
                                     std::uint32_t q21, bool merge_phase) {
            std::uint32_t sA = std::min(q11, q22), sB = std::max(q11, q22);
            std::uint32_t cA = std::min(q12, q21), cB = std::max(q12, q21);
            if (merge_phase && (std::make_pair(cA, cB) < std::make_pair(sA, sB))) {
                std::swap(sA, cA); std::swap(sB, cB);   // the mutation: sorts ACROSS the phases
            }
            std::string k(16, '\0');
            std::memcpy(&k[0], &sA, 4); std::memcpy(&k[4], &sB, 4);
            std::memcpy(&k[8], &cA, 4); std::memcpy(&k[12], &cB, 4);
            return k;
        };
        // THE COMPARISON over every content class: a pattern's delta is computed ONCE, from the
        // first class realising it, and every later member must reproduce it exactly.
        const auto run = [&](int sig_mode, bool merge_phase, bool members_matter,
                             double* worst_out, std::size_t* patterns_out) {
            std::unordered_map<std::string, double> cache;
            double worst = 0.0;
            for (std::size_t a1 = 0; a1 + 1 < NA; ++a1)
            for (std::size_t a2 = a1 + 1; a2 < NA; ++a2)
            for (std::size_t b1 = 0; b1 + 1 < NB; ++b1)
            for (std::size_t b2 = b1 + 1; b2 < NB; ++b2) {
                std::uint32_t q11 = M.signature_at(a1, b1), q22 = M.signature_at(a2, b2);
                std::uint32_t q12 = M.signature_at(a1, b2), q21 = M.signature_at(a2, b1);
                if (sig_mode == 1) {   // the mutation: signature ignores WHICH alleles, using only
                    q11 = q22 = q12 = q21 = 0;   // membership-free grouping
                }
                const std::string k = pattern_key(q11, q22, q12, q21, merge_phase);
                const double d = delta_direct(a1, a2, b1, b2, merge_phase ? 2 : 0);
                const auto it = cache.find(k);
                if (it == cache.end()) cache.emplace(k, d);
                else worst = std::max(worst, std::abs(it->second - d));
                (void)members_matter;
            }
            if (worst_out != nullptr) *worst_out = worst;
            if (patterns_out != nullptr) *patterns_out = cache.size();
            return worst;
        };
        double worst = 0.0; std::size_t npat = 0;
        run(0, false, true, &worst, &npat);
        ok_(worst == kExact,
            "every content class sharing a pattern has the IDENTICAL signed delta (worst " +
            sci(worst) + ", tolerance " + sci(kExact) + ", over " + std::to_string(npat) + " patterns)");

        // NON-VACUITY: some delta must be non-zero, or the equality above is trivially satisfied.
        double maxabs = 0.0;
        for (std::size_t a1 = 0; a1 + 1 < NA; ++a1)
        for (std::size_t a2 = a1 + 1; a2 < NA; ++a2)
        for (std::size_t b1 = 0; b1 + 1 < NB; ++b1)
        for (std::size_t b2 = b1 + 1; b2 < NB; ++b2)
            maxabs = std::max(maxabs, std::abs(delta_direct(a1, a2, b1, b2, 0)));
        ok_(maxabs > 1e-6, "the fixture carries real phase signal (max |delta| " +
                           sci(maxabs) + ")");

        // THE SIGN GATE. Swapping the two B alleles exchanges straight with crossed, so the delta
        // must be the exact negative. Merging the phases would make these equal instead.
        double worst_sign = 0.0;
        for (std::size_t a1 = 0; a1 + 1 < NA; ++a1)
        for (std::size_t a2 = a1 + 1; a2 < NA; ++a2)
        for (std::size_t b1 = 0; b1 + 1 < NB; ++b1)
        for (std::size_t b2 = b1 + 1; b2 < NB; ++b2) {
            const double d1 = delta_direct(a1, a2, b1, b2, 0);
            const double d2 = delta_direct(a1, a2, b2, b1, 0);
            worst_sign = std::max(worst_sign, std::abs(d1 + d2));
        }
        ok_(worst_sign < 1e-9,
            "reversing straight and crossed gives exactly the opposite delta (worst " +
            sci(worst_sign) + ", tolerance 1.000e-09)");

        // THE FLAT GATE. Two alleles in one row class are indistinguishable to every fragment, so
        // straight and crossed coincide and the class must be neutral.
        double worst_flat = 0.0; std::size_t n_flat = 0;
        for (const auto& mem : M.row_members) {
            if (mem.size() < 2) continue;
            for (std::size_t b1 = 0; b1 + 1 < NB; ++b1)
            for (std::size_t b2 = b1 + 1; b2 < NB; ++b2) {
                ++n_flat;
                worst_flat = std::max(worst_flat,
                                      std::abs(delta_direct(mem[0], mem[1], b1, b2, 0)));
            }
        }
        ok_(n_flat > 0 && worst_flat < 1e-12,
            "a class whose two A alleles share a signature row is exactly neutral (" +
            std::to_string(n_flat) + " classes, worst " + sci(worst_flat) + ", tolerance 1.000e-12)");

        // THE COLUMN GATE, symmetric to the row one. Two B alleles sharing a signature column are
        // indistinguishable to every fragment, so their class must be neutral too. Testing only the
        // row direction would leave an asymmetric implementation passing.
        double worst_flat_b = 0.0; std::size_t n_flat_b = 0;
        for (const auto& mem : M.col_members) {
            if (mem.size() < 2) continue;
            for (std::size_t a1 = 0; a1 + 1 < NA; ++a1)
            for (std::size_t a2 = a1 + 1; a2 < NA; ++a2) {
                ++n_flat_b;
                worst_flat_b = std::max(worst_flat_b,
                                        std::abs(delta_direct(a1, a2, mem[0], mem[1], 0)));
            }
        }
        ok_(n_flat_b > 0 && worst_flat_b < 1e-12,
            "a class whose two B alleles share a signature column is exactly neutral (" +
            std::to_string(n_flat_b) + " classes, worst " + sci(worst_flat_b) + ", tolerance 1.000e-12)");

        // ---- END TO END: the grouped delta must survive the WEIGHTED chain computation ----------
        // Member-resolved factor equivalence is not member-weighted HMM equivalence: marker
        // unaries, Li-Stephens weights and the forward/backward messages have not yet passed
        // through the grouped representation. They do here, against deliberately adversarial
        // weights -- unequal unaries among alleles that share a signature class, unequal class
        // sizes, and r at both ends of its range as well as in between.
        {
            // The sparse edge's own key shape. It is file-local to genotype_fragments.cpp, so it
            // is restated here -- and if the restatement were wrong, log_psi would find no entry,
            // every class would read as neutral, and the non-vacuity gate below would fail. The
            // test cannot silently pass on a mismatched key.
            const auto ckey = [](std::size_t amin, std::size_t amax,
                                 std::size_t bmin, std::size_t bmax) {
                return (static_cast<std::uint64_t>(amin) << 48) |
                       (static_cast<std::uint64_t>(amax) << 32) |
                       (static_cast<std::uint64_t>(bmin) << 16) |
                       static_cast<std::uint64_t>(bmax);
            };
            SparseLinkageEdge orig, grp;
            orig.block_a = 0; orig.block_b = 1; orig.n_a = NA; orig.n_b = NB;
            grp = orig;
            std::unordered_map<std::string, double> pcache;
            for (std::size_t a1 = 0; a1 + 1 < NA; ++a1)
            for (std::size_t a2 = a1 + 1; a2 < NA; ++a2)
            for (std::size_t b1 = 0; b1 + 1 < NB; ++b1)
            for (std::size_t b2 = b1 + 1; b2 < NB; ++b2) {
                const double d = delta_direct(a1, a2, b1, b2, 0);
                if (d != 0.0) orig.delta.emplace(ckey(a1, a2, b1, b2), d);
                // THE GROUPED PATH: one delta per pattern, reused by every member class.
                const std::string k = pattern_key(M.signature_at(a1, b1), M.signature_at(a2, b2),
                                                  M.signature_at(a1, b2), M.signature_at(a2, b1),
                                                  false);
                auto it = pcache.find(k);
                if (it == pcache.end()) it = pcache.emplace(k, d).first;
                if (it->second != 0.0) grp.delta.emplace(ckey(a1, a2, b1, b2), it->second);
            }
            orig.status = LinkageStatus::Ok; grp.status = LinkageStatus::Ok;
            orig.stored_classes = orig.delta.size(); grp.stored_classes = grp.delta.size();
            bool maps_equal = orig.delta.size() == grp.delta.size();
            if (maps_equal)
                for (const auto& kv : orig.delta) {
                    const auto it = grp.delta.find(kv.first);
                    if (it == grp.delta.end() || it->second != kv.second) { maps_equal = false; break; }
                }
            ok_(maps_equal && !orig.delta.empty(),
                "the grouped delta map equals the direct one entry for entry (" +
                std::to_string(orig.delta.size()) + " classes, " +
                std::to_string(pcache.size()) + " patterns)");
            // THE PRODUCTION BUILDERS, compared against each other rather than against a map
            // assembled by the test. This is what production will actually run.
            {
                GroupedBuildStats gs;
                const SparseLinkageEdge ub =
                    build_sparse_linkage_edge(ems, g, lambda, log_mix, log_bg_weight);
                const SparseLinkageEdge gb = build_sparse_linkage_edge_grouped(
                    ems, sigs, g, lambda, log_mix, log_bg_weight, {}, &gs);
                ok_(ub.usable() && gb.usable() && gb.grouped,
                    "both builders produce a usable edge, and the grouped one is grouped");
                ok_(gs.oracle_visits == 0,
                    "the grouped build NEVER enumerates a content class (oracle visits " +
                    std::to_string(gs.oracle_visits) + ")");
                double wp = 0.0;
                for (std::size_t a1 = 0; a1 < NA; ++a1)
                for (std::size_t a2 = 0; a2 < NA; ++a2)
                for (std::size_t b1 = 0; b1 < NB; ++b1)
                for (std::size_t b2 = 0; b2 < NB; ++b2)
                    wp = std::max(wp, std::abs(ub.log_psi(a1, b1, a2, b2) -
                                               gb.log_psi(a1, b1, a2, b2)));
                ok_(wp == kExact,
                    "the GROUPED BUILDER's psi equals the ungrouped builder's for every ordered "
                    "content class (worst " + sci(wp) + ", tolerance " + sci(kExact) + "); " +
                    std::to_string(ub.stored_classes) + " content classes vs " +
                    std::to_string(gb.stored_classes) + " class quadruples, " +
                    std::to_string(gs.representative_visits) + " representative visits");
            }
            // PSI for every real content class, both orientations.
            double worst_psi = 0.0;
            for (std::size_t a1 = 0; a1 < NA; ++a1)
            for (std::size_t a2 = 0; a2 < NA; ++a2)
            for (std::size_t b1 = 0; b1 < NB; ++b1)
            for (std::size_t b2 = 0; b2 < NB; ++b2)
                worst_psi = std::max(worst_psi, std::abs(orig.log_psi(a1, b1, a2, b2) -
                                                          grp.log_psi(a1, b1, a2, b2)));
            ok_(worst_psi == kExact, "log psi is identical for every ordered content class (worst " +
                                  sci(worst_psi) + ", tolerance " + sci(kExact) + ")");
            // THE WEIGHTED CHAIN. Six haplotypes, deliberately spanning signature classes so that
            // members of one class carry different unaries and different forward mass.
            const std::size_t nh = 6, nbk = 3, ns = nh * nh;
            AlleleMapping ma, mb;
            ma.n_alleles = NA; mb.n_alleles = NB;
            ma.status = MappingStatus::Ok; mb.status = MappingStatus::Ok;
            ma.allele = {0, 1, 2, 3, 4, 5};        // one haplotype per A allele, classes 1/2/3
            mb.allele = {0, 1, 2, 3, 5, 7};        // spanning all four B classes
            std::vector<std::vector<double>> em(nbk, std::vector<double>(ns, 1.0));
            for (std::size_t b = 0; b < nbk; ++b)
                for (std::size_t i = 0; i < nh; ++i)
                    for (std::size_t j = 0; j < nh; ++j)
                        em[b][i * nh + j] = std::exp(-0.37 * (b + 1) * (i + 1) -
                                                     0.11 * (j + 2) * (b + 3) - 0.05 * i * j);
            const auto emit_fn = [&](std::size_t b, std::vector<double>& ev) { ev = em[b]; };
            const SparseEdgeLinkage se_o = make_sparse_kernel_edge(orig, ma, mb);
            const SparseEdgeLinkage se_g = make_sparse_kernel_edge(grp, ma, mb);
            // SWAP INVARIANCE NEEDS SYMMETRIC EMISSIONS. The adversarial emissions above are
            // asymmetric in i and j on purpose -- which is right for the marginal and partition
            // identities and useless here: (i,j) and (j,i) then carry different mass with NO
            // linkage at all, and the check would measure the fixture instead of the kernel.
            std::vector<std::vector<double>> sym(nbk, std::vector<double>(ns, 1.0));
            for (std::size_t b = 0; b < nbk; ++b)
                for (std::size_t i = 0; i < nh; ++i)
                    for (std::size_t j = 0; j < nh; ++j) {
                        const double v = std::exp(-0.21 * (b + 1) * (i + 1) * (j + 1) -
                                                  0.05 * (i + j));
                        sym[b][i * nh + j] = v; sym[b][j * nh + i] = v;
                    }
            const auto emit_sym = [&](std::size_t b, std::vector<double>& ev) { ev = sym[b]; };
            double worst_marg = 0.0, worst_z = 0.0, worst_swap = 0.0, max_shift = 0.0;
            for (double r : {0.0, 0.3, 1.0}) {
                std::vector<SparseEdgeLinkage> so(nbk), sg(nbk), snone(nbk);
                so[2] = se_o; sg[2] = se_g;
                std::vector<std::vector<double>> f1, b1v, f2, b2v, f0, b0v;
                ChainKernelStats s1, s2, s0;
                chain_forward_backward(nh, nbk, r, emit_fn, nullptr, f1, b1v, &s1, &so);
                chain_forward_backward(nh, nbk, r, emit_fn, nullptr, f2, b2v, &s2, &sg);
                chain_forward_backward(nh, nbk, r, emit_fn, nullptr, f0, b0v, &s0, &snone);
                worst_z = std::max(worst_z, std::abs(s1.log_weight_sum - s2.log_weight_sum));
                for (std::size_t b = 0; b < nbk; ++b) {
                    double z1 = 0.0, z2 = 0.0, z0 = 0.0;
                    for (std::size_t k = 0; k < ns; ++k) {
                        z1 += f1[b][k] * b1v[b][k];
                        z2 += f2[b][k] * b2v[b][k];
                        z0 += f0[b][k] * b0v[b][k];
                    }
                    for (std::size_t i = 0; i < nh; ++i)
                    for (std::size_t j = 0; j < nh; ++j) {
                        const std::size_t k = i * nh + j, ks = j * nh + i;
                        const double p1 = z1 > 0 ? f1[b][k] * b1v[b][k] / z1 : 0.0;
                        const double p2 = z2 > 0 ? f2[b][k] * b2v[b][k] / z2 : 0.0;
                        const double p0 = z0 > 0 ? f0[b][k] * b0v[b][k] / z0 : 0.0;
                        const double p2s = z2 > 0 ? f2[b][ks] * b2v[b][ks] / z2 : 0.0;
                        worst_marg = std::max(worst_marg, std::abs(p1 - p2));
                        (void)p2s;
                        max_shift = std::max(max_shift, std::abs(p2 - p0));
                    }
                }
            }
            ok_(worst_marg <= kTolMarg, "every block marginal agrees across r in {0, 0.3, 1} (worst " +
                                    sci(worst_marg) + ", tolerance " + sci(kTolMarg) + ")");
            ok_(worst_z <= kTolZ, "the unnormalised partition weight agrees (worst " +
                                 sci(worst_z) + ", tolerance " + sci(kTolZ) + ")");
            for (double r : {0.0, 0.3, 1.0}) {
                std::vector<SparseEdgeLinkage> sg(nbk);
                sg[2] = se_g;
                std::vector<std::vector<double>> f2, b2v;
                ChainKernelStats s2;
                chain_forward_backward(nh, nbk, r, emit_sym, nullptr, f2, b2v, &s2, &sg);
                for (std::size_t b = 0; b < nbk; ++b) {
                    double z2 = 0.0;
                    for (std::size_t k = 0; k < ns; ++k) z2 += f2[b][k] * b2v[b][k];
                    for (std::size_t i = 0; i < nh; ++i)
                    for (std::size_t j = 0; j < nh; ++j) {
                        const std::size_t k = i * nh + j, ks = j * nh + i;
                        const double p = z2 > 0 ? f2[b][k] * b2v[b][k] / z2 : 0.0;
                        const double ps = z2 > 0 ? f2[b][ks] * b2v[b][ks] / z2 : 0.0;
                        worst_swap = std::max(worst_swap, std::abs(p - ps));
                    }
                }
            }
            ok_(worst_swap <= kTolZ, "homologue swap invariance holds under the grouped edge, on "
                                    "SYMMETRIC emissions (worst " +
                                    sci(worst_swap) + ", tolerance " + sci(kTolZ) + ")");
            // NON-VACUITY: the linkage must actually move a posterior, or every identity above is
            // an identity between two copies of the linkage-free answer.
            ok_(max_shift > 1e-6, "linkage MATERIALLY changes a posterior against no edge (max " +
                                  sci(max_shift) + ")");
        }

        // MUTATION 1: drop mismatch counts from the signature. Cells differing only by a mismatch
        // then share a signature, and their deltas must stop agreeing.
        {
            std::vector<std::string> jm(NA * NB);
            for (std::size_t k = 0; k < NA * NB; ++k) {
                for (const auto& cs : sigs) {
                    // keep only the insert field of each triple, discarding both edit counts
                    for (std::size_t j = 0; j + 12 <= cs[k].size(); j += 12)
                        jm[k].append(&cs[k][j + 8], 4);
                    jm[k].push_back('|');
                }
            }
            const SignatureMatrix Mm = build_signature_matrix(jm, NA, NB);
            std::unordered_map<std::string, double> cache;
            double w = 0.0;
            for (std::size_t a1 = 0; a1 + 1 < NA; ++a1)
            for (std::size_t a2 = a1 + 1; a2 < NA; ++a2)
            for (std::size_t b1 = 0; b1 + 1 < NB; ++b1)
            for (std::size_t b2 = b1 + 1; b2 < NB; ++b2) {
                const std::string k = pattern_key(Mm.signature_at(a1, b1), Mm.signature_at(a2, b2),
                                                  Mm.signature_at(a1, b2), Mm.signature_at(a2, b1),
                                                  false);
                const double d = delta_direct(a1, a2, b1, b2, 0);
                const auto it = cache.find(k);
                if (it == cache.end()) cache.emplace(k, d);
                else w = std::max(w, std::abs(it->second - d));
            }
            ok_(w > 1e-6, "MUTATION removing mismatch counts BREAKS the equality (worst " +
                          sci(w) + ", must exceed 1.000e-06)");
        }
        // MUTATION 2: merge straight with crossed in the key. The sign gate must then fail.
        {
            double w = 0.0;
            std::unordered_map<std::string, double> cache;
            for (std::size_t a1 = 0; a1 + 1 < NA; ++a1)
            for (std::size_t a2 = a1 + 1; a2 < NA; ++a2)
            for (std::size_t b1 = 0; b1 + 1 < NB; ++b1)
            for (std::size_t b2 = b1 + 1; b2 < NB; ++b2) {
                const std::string k = pattern_key(M.signature_at(a1, b1), M.signature_at(a2, b2),
                                                  M.signature_at(a1, b2), M.signature_at(a2, b1),
                                                  true);
                const double d = delta_direct(a1, a2, b1, b2, 0);
                const auto it = cache.find(k);
                if (it == cache.end()) cache.emplace(k, d);
                else w = std::max(w, std::abs(it->second - d));
            }
            ok_(w > 1e-6, "MUTATION merging straight with crossed BREAKS the equality (worst " +
                          sci(w) + ", must exceed 1.000e-06)");
        }
        // MUTATION 3: ignore membership entirely -- every cell gets the same signature id, so all
        // content classes collapse into one pattern.
        {
            double w = 0.0; std::size_t np = 0;
            run(1, false, false, &w, &np);
            ok_(w > 1e-6 && np == 1,
                "MUTATION ignoring membership BREAKS the equality (" + std::to_string(np) +
                " pattern, worst " + sci(w) + ")");
        }
        std::printf("grouping selftest: %zu failure(s)\n", fails);
        return fails == 0 ? 0 : 1;
    }

    // ---- OPERATIONAL WORK BUDGET SELF-TEST ------------------------------------------------------
    // The budget must bound work BEFORE it is done, and exhaustion must be a refusal carrying no
    // mass. A budget checked afterwards, or one that returns a partially enumerated support, is
    // worse than none: it reports a smaller model rather than an unfinished one.
    if (budget_selftest) {
        std::size_t fails = 0;
        const auto ok_ = [&](bool c, const char* what) {
            std::printf("%s\t%s\n", c ? "ok" : "FAIL", what);
            if (!c) ++fails;
        };
        std::mt19937_64 rng(20260908);
        const auto rseq = [&](std::size_t n) {
            static const char* B = "ACGT";
            std::string t(n, 'A');
            for (std::size_t i = 0; i < n; ++i) t[i] = B[rng() & 3];
            return t;
        };
        LinkageGeometry g;
        g.block_a = 0; g.block_b = 1;
        g.alleles_a = {rseq(120), rseq(120), rseq(120)};
        g.alleles_b = {rseq(140), rseq(140), rseq(140)};
        g.context = rseq(60); g.lflank = rseq(800); g.rflank = rseq(800);
        g.window_len.assign(9, 0); g.exposure.assign(9, 1000.0);
        g.exposure_affine = true; g.ok = true;
        const auto ix = build_allele_product_index(g, 16);
        InsertPrior ip; ip.lo = 200; ip.hi = 600;
        ip.logp.assign(static_cast<std::size_t>(ip.hi - ip.lo + 1),
                       -std::log(static_cast<double>(ip.hi - ip.lo + 1)));
        VirtualWindow vw; vw.bind_pair(g, 1, 2);
        const std::size_t at = g.lflank.size() - 60;
        Fragment f; f.name = "b";
        for (std::size_t i = 0; i < 150; ++i) f.r1.push_back(vw.base_at(at + i));
        std::string tail;
        for (std::size_t i = 0; i < 150; ++i) tail.push_back(vw.base_at(at + 250 + i));
        f.r2 = reverse_complement(tail);
        const double lep = std::log(0.001 / 3.0), l1m = std::log(1.0 - 0.001);

        // (1) INERT WHEN NOT BINDING: a generous budget must not change one bit of the answer.
        const auto ref = linkage_emission_supported(f, g, ip, 0.05, lep, l1m, -300.0, nullptr, &ix,
                                                    nullptr);
        HybridWorkBudget big; big.max_proposed_cells = 0; big.max_full_read_verifications = 0;
        const auto same = linkage_emission_supported(f, g, ip, 0.05, lep, l1m, -300.0, nullptr, &ix,
                                                     &big);
        ok_(!same.work_refused && same.mass == ref.mass,
            "an unlimited budget leaves the emission bit-identical");
        ok_(big.full_read_verifications > 0 && big.proposed_cells > 0 &&
            big.bases_compared_upper_bound >= big.full_read_verifications,
            "the budget counts the work it observed, and bases >= verifications");

        // (2) THE CELL LIMIT refuses, and carries NO mass.
        HybridWorkBudget c1; c1.max_proposed_cells = 1; c1.max_full_read_verifications = 0;
        const auto rc1 = linkage_emission_supported(f, g, ip, 0.05, lep, l1m, -300.0, nullptr, &ix,
                                                    &c1);
        bool no_mass = !rc1.mass.empty();
        for (double m : rc1.mass) if (m != -std::numeric_limits<double>::infinity()) no_mass = false;
        ok_(rc1.work_refused && !rc1.ok && no_mass &&
            rc1.work_refusal == "support-search-proposed-cell-limit",
            "the proposed-cell limit refuses with no mass at all");

        // (3) THE VERIFICATION LIMIT is charged BEFORE the comparison, so it stops the work.
        HybridWorkBudget v1; v1.max_proposed_cells = 0; v1.max_full_read_verifications = 2;
        const auto rv1 = linkage_emission_supported(f, g, ip, 0.05, lep, l1m, -300.0, nullptr, &ix,
                                                    &v1);
        // THE ORDERING IS THE POINT. Charged BEFORE the comparison, the number of whole-read
        // comparisons actually PERFORMED never exceeds the limit. Charged after, the limit is
        // discovered by overrunning it -- one comparison too late, every time.
        ok_(rv1.work_refused && rv1.work_refusal == "support-search-verification-limit" &&
            rv1.full_read_verifications <= v1.max_full_read_verifications,
            "no more whole-read comparisons are PERFORMED than the limit allows");

        // (4) THE FALLBACK MUST NOT BEGIN. A non-ACGT read forces the dense path; with a budget it
        //     cannot afford, no window may be enumerated at all.
        Fragment nf = f; nf.r1[10] = 'N';
        HybridWorkBudget fb; fb.max_proposed_cells = 4; fb.max_full_read_verifications = 10;
        const auto rfb = linkage_emission_supported(nf, g, ip, 0.05, lep, l1m, -300.0, nullptr, &ix,
                                                    &fb);
        bool fb_empty = !rfb.mass.empty();
        for (double m : rfb.mass) if (m != -std::numeric_limits<double>::infinity()) fb_empty = false;
        ok_(rfb.work_refused && rfb.work_refusal == "support-search-fallback-work-limit" &&
            fb_empty && fb.full_read_verifications == 0,
            "an unaffordable dense fallback refuses BEFORE enumerating any window");
        // And the same fallback with room must actually produce the dense answer.
        HybridWorkBudget fb2; fb2.max_proposed_cells = 0; fb2.max_full_read_verifications = 0;
        const auto rfb2 = linkage_emission_supported(nf, g, ip, 0.05, lep, l1m, -300.0, nullptr, &ix,
                                                     &fb2);
        ok_(!rfb2.work_refused && rfb2.ok,
            "the same fallback with budget to spare still computes the dense emission");
        std::printf("budget selftest: %zu failure(s)\n", fails);
        return fails == 0 ? 0 : 1;
    }

    // ---- COORDINATE VALIDATION SELF-TEST -------------------------------------------------------
    // The derived start is the one thing direct verification cannot recover from being wrong: a
    // misplaced coordinate does not fail loudly, it silently verifies the wrong bases and drops the
    // placement. So check the coordinate itself, three ways.
    if (coordinate_selftest) {
        std::size_t fails = 0;
        const auto ok_ = [&](bool c, const char* what) {
            std::printf("%s\t%s\n", c ? "ok" : "FAIL", what);
            if (!c) ++fails;
        };
        std::mt19937_64 rng(20260908);
        const auto rseq = [&](std::size_t n) {
            static const char* B = "ACGT";
            std::string t(n, 'A');
            for (std::size_t i = 0; i < n; ++i) t[i] = B[rng() & 3];
            return t;
        };
        const auto geom_of = [&](const std::vector<std::string>& A,
                                 const std::vector<std::string>& Bv, const std::string& ctx,
                                 const std::string& lf, const std::string& rf) {
            LinkageGeometry g;
            g.block_a = 0; g.block_b = 1;
            g.alleles_a = A; g.alleles_b = Bv; g.context = ctx; g.lflank = lf; g.rflank = rf;
            g.window_len.assign(A.size() * Bv.size(), 0);
            g.exposure.assign(A.size() * Bv.size(), 1000.0);
            g.exposure_affine = true; g.ok = true;
            return g;
        };

        // (1) THE VIEW AND THE STRING AGREE. materialize() is the oracle's window; base_at() is the
        // verifier's. If these ever diverge the two paths silently score different sequence.
        {
            const auto g = geom_of({rseq(40), rseq(400)}, {rseq(90), rseq(7)}, rseq(30),
                                   rseq(600), rseq(600));
            bool same = true;
            for (std::uint32_t al = 0; al < 2 && same; ++al) {
                for (std::uint32_t be = 0; be < 2 && same; ++be) {
                    VirtualWindow vw; vw.bind_pair(g, al, be);
                    const std::string m = vw.materialize();
                    if (m.size() != vw.size()) { same = false; break; }
                    for (std::size_t i = 0; i < m.size(); ++i)
                        if (m[i] != vw.base_at(i)) { same = false; break; }
                }
            }
            ok_(same, "base_at agrees with materialize on every position of every window");
        }

        // (2) THE LENGTH SHIFT, isolated. One B allele, two A alleles differing greatly in length,
        // and a read planted INSIDE B. The two derived starts must differ by exactly |A2|-|A1| --
        // that shift is the whole positional content of a B-side seed.
        {
            const std::string A1 = rseq(50), A2 = rseq(950);
            const auto g = geom_of({A1, A2}, {rseq(500)}, rseq(40), rseq(700), rseq(700));
            const std::size_t piece = 16;
            const auto ix = build_allele_product_index(g, piece);
            VirtualWindow w1; w1.bind_pair(g, 0, 0);
            VirtualWindow w2; w2.bind_pair(g, 1, 0);
            // A 150 bp read wholly inside B of window (alpha=0), and the SAME sequence inside
            // window (alpha=1), where it sits |A2|-|A1| further along.
            const std::size_t inb = g.lflank.size() + A1.size() + g.context.size() + 120;
            std::string read;
            for (std::size_t i = 0; i < 150; ++i) read.push_back(w1.base_at(inb + i));
            Fragment f; f.name = "c";
            f.r1 = read;
            std::string tail;
            for (std::size_t i = 0; i < 150; ++i) tail.push_back(w1.base_at(inb + 200 + i));
            f.r2 = reverse_complement(tail);
            InsertPrior ip; ip.lo = 200; ip.hi = 600;
            ip.logp.assign(static_cast<std::size_t>(ip.hi - ip.lo + 1),
                              -std::log(static_cast<double>(ip.hi - ip.lo + 1)));
            const auto sup = propose_allele_pairs(f, g, ip, 0.05, &ix);
            long s0 = -1, s1 = -1;
            for (const auto& st : sup.mate_states[0]) {
                const std::uint32_t al = static_cast<std::uint32_t>(st.first >> 32);
                if (al == 0 && w1.count_mismatches(f.r1, st.second, 0) == 0) s0 = st.second;
                if (al == 1 && w2.count_mismatches(f.r1, st.second, 0) == 0) s1 = st.second;
            }
            ok_(s0 == static_cast<long>(inb), "the B-side read verifies at its planted start");
            ok_(s1 >= 0 && s1 - s0 == static_cast<long>(A2.size()) - static_cast<long>(A1.size()),
                "the two derived B starts differ by exactly |A2| - |A1|");
        }

        // (3) THE SEED INVARIANT ON REAL PROPOSALS: at every proposed start that verifies, the
        // window content under the read equals the read. Checked through the view, on a geometry
        // with an EMPTY allele -- the case the previous six-map index could not seed at all.
        {
            const std::string A1 = rseq(60);
            const auto g = geom_of({A1, ""}, {rseq(80), rseq(3)}, rseq(25), rseq(500), rseq(500));
            const std::size_t piece = 16;
            const auto ix = build_allele_product_index(g, piece);
            ok_(ix.ok && ix.complete, "an empty allele still yields a COMPLETE index");
            VirtualWindow w; w.bind_pair(g, 1, 1);   // empty A, short B
            const std::size_t at = g.lflank.size() - 40;             // spans L, (empty A), C, B, R
            std::string read;
            for (std::size_t i = 0; i < 150; ++i) read.push_back(w.base_at(at + i));
            Fragment f; f.name = "e"; f.r1 = read;
            std::string tail;
            for (std::size_t i = 0; i < 150; ++i) tail.push_back(w.base_at(at + 200 + i));
            f.r2 = reverse_complement(tail);
            InsertPrior ip; ip.lo = 200; ip.hi = 600;
            ip.logp.assign(static_cast<std::size_t>(ip.hi - ip.lo + 1),
                              -std::log(static_cast<double>(ip.hi - ip.lo + 1)));
            const auto sup = propose_allele_pairs(f, g, ip, 0.05, &ix);
            bool found = false, all_consistent = true;
            for (const auto& st : sup.mate_states[0]) {
                const std::uint32_t al = static_cast<std::uint32_t>(st.first >> 32);
                const std::uint32_t be = static_cast<std::uint32_t>(st.first & 0xFFFFFFFFu);
                VirtualWindow v; v.bind_pair(g, al, be);
                if (v.count_mismatches(f.r1, st.second, 0) != 0) continue;
                for (std::size_t i = 0; i < f.r1.size(); ++i)
                    if (v.base_at(static_cast<std::size_t>(st.second) + i) != f.r1[i])
                        all_consistent = false;
                if (al == 1 && be == 1 && st.second == static_cast<long>(at)) found = true;
            }
            ok_(found, "a read spanning an EMPTY allele is seeded at its true start");
            ok_(all_consistent, "window content under every verified start equals the read");
        }
        std::printf("coordinate selftest: %zu failure(s)\n", fails);
        return fails == 0 ? 0 : 1;
    }

    // ---- ALLELE-PRODUCT SUPPORT SEARCH SELF-TEST ------------------------------------------------
    // The support-restricted emission must equal the DENSE ORACLE cell for cell: same finite cells,
    // same log mass, same informative classification. Not "the same phase call" -- that would pass
    // while multiplicity or an off-panel combination went missing.
    if (support_selftest) {
        // STAGE COUNTERS reported separately, so it is visible whether cost moves out of window
        // verification and into tuple construction or joining rather than disappearing.
        std::printf("case\tn_a\tn_b\tdense_pairs\tproposed\tverified\tseed_occurrences"
                    "\tfallback\tcells_differ\tworst_mass_diff\tfinite_dense\tfinite_supported"
                    "\tinformative_match\treduction\tseed_start_proposals\tunique_seed_starts"
                    "\tseed_compatible_joins\tfull_read_verifications\taccepted_mate_placements"
                    "\tverified_fr_states\tfinite_emission_cells\n");
        std::mt19937_64 rng(20260908);
        const auto rseq = [&](std::size_t n) {
            static const char* B = "ACGT";
            std::string s2(n, 'A');
            for (std::size_t i = 0; i < n; ++i) s2[i] = B[rng() & 3];
            return s2;
        };
        const auto mutate = [&](std::string x, std::size_t k) {
            for (std::size_t i = 0; i < k; ++i) {
                const std::size_t p2 = rng() % x.size();
                const char c = "ACGT"[rng() & 3];
                x[p2] = (x[p2] == c) ? "ACGT"[(rng() & 3)] : c;
            }
            return x;
        };
        const double lep = std::log(0.001 / 3.0), l1m = std::log1p(-0.001);
        const auto run = [&](const char* name, const LinkageGeometry& g, const Fragment& f,
                             const InsertPrior& ip) {
            const LinkageEmission D = linkage_emission(f, g, ip, 0.05, lep, l1m, -400.0);
            // THE INDEX IS BUILT HERE, at this fragment's own piece length. Without it the search
            // takes the exhaustive fallback and the comparison silently stops testing the search.
            const std::size_t dd = mate_band_edits(0.05, f.r1.size());
            const std::size_t pp = f.r1.size() / (dd + 1);
            const AlleleProductIndex aix = (pp >= 8 && pp <= 32)
                ? build_allele_product_index(g, pp) : AlleleProductIndex{};
            AlleleProductSupport sup;
            const LinkageEmission S2 = linkage_emission_supported(f, g, ip, 0.05, lep, l1m,
                                                                  -400.0, &sup,
                                                                  aix.ok ? &aix : nullptr);
            std::size_t differ = 0, fin_d = 0, fin_s = 0;
            double worst = 0.0;
            for (std::size_t k = 0; k < D.mass.size(); ++k) {
                const double x = D.mass[k], y = S2.mass[k];
                const bool fx = x != -std::numeric_limits<double>::infinity();
                const bool fy = y != -std::numeric_limits<double>::infinity();
                if (fx) ++fin_d;
                if (fy) ++fin_s;
                if (fx != fy) { ++differ; continue; }
                if (fx) worst = std::max(worst, std::abs(x - y));
            }
            const std::size_t ver = sup.exhaustive_fallback ? sup.dense_pairs
                                                            : sup.proposals.size();
            std::printf("%s\t%zu\t%zu\t%zu\t%zu\t%zu\t%zu\t%d\t%zu\t%.3g\t%zu\t%zu"
                        "\t%d\t%.2f\t%zu\t%zu\t%zu\t%zu\t%zu\t%zu\t%zu\n",
                        name, D.n_a, D.n_b, sup.dense_pairs, sup.proposed_states, ver,
                        sup.seed_occurrences, sup.exhaustive_fallback ? 1 : 0,
                        differ, worst, fin_d, fin_s,
                        D.informative == S2.informative ? 1 : 0,
                        ver == 0 ? 0.0 : static_cast<double>(sup.dense_pairs) / ver,
                        sup.seed_start_proposals, sup.unique_seed_starts,
                        sup.seed_compatible_joins, S2.full_read_verifications,
                        S2.accepted_mate_placements, S2.verified_fr_states,
                        S2.finite_emission_cells);
        };
        const auto geom_of = [&](const std::vector<std::string>& A,
                                 const std::vector<std::string>& Bv, const std::string& ctx,
                                 const std::string& lf, const std::string& rf) {
            LinkageGeometry g;
            g.block_a = 0; g.block_b = 1;
            g.alleles_a = A; g.alleles_b = Bv; g.context = ctx; g.lflank = lf; g.rflank = rf;
            g.window_len.assign(A.size() * Bv.size(), 0);
            g.exposure.assign(A.size() * Bv.size(), 1000.0);
            g.exposure_affine = true; g.ok = true;
            return g;
        };
        const auto frag_from = [&](const std::string& win, std::size_t start, std::size_t ins) {
            Fragment f;
            f.name = "f";
            f.r1 = win.substr(start, 100);
            const std::string tail = win.substr(start + ins - 100, 100);
            f.r2 = reverse_complement(tail);
            return f;
        };
        InsertPrior ip = make_insert_prior(300.0, 40.0, 0.01, 4, 200);

        // Three A alleles and three B alleles; the sample's truth (A2, B1) is a combination NO
        // panel haplotype need carry. This is the case that fails if support is taken from panel
        // placements rather than from seeds.
        const std::string L = rseq(400), C = rseq(60), R = rseq(400);
        std::vector<std::string> A = {rseq(500), rseq(500), rseq(500)};
        std::vector<std::string> Bv = {rseq(500), rseq(500), rseq(500)};
        const LinkageGeometry g3 = geom_of(A, Bv, C, L, R);
        {
            // OFF-PANEL COMBINATION: read pair drawn from A2 + C + B1.
            const std::string win = L + A[2] + C + Bv[1] + R;
            run("offpanel_A2_B1", g3, frag_from(win, 850, 300), ip);
            // A SEED CROSSING THE DIRECT A->B JUNCTION, context emptied.
            const LinkageGeometry g0 = geom_of(A, Bv, "", L, R);
            const std::string w0 = L + A[1] + Bv[2] + R;
            run("junction_crossing_seed", g0, frag_from(w0, 860, 280), ip);
            // A SEED SPANNING A -> short context -> B.
            const LinkageGeometry gs = geom_of(A, Bv, rseq(12), L, R);
            const std::string ws = L + A[0] + gs.context + Bv[0] + R;
            run("seed_spans_A_ctx_B", gs, frag_from(ws, 870, 260), ip);
            // A-ONLY: both mates land inside A, leaving beta unconstrained.
            run("A_only_constraint", g3, frag_from(win, 450, 260), ip);
            // B-ONLY: both mates inside B.
            run("B_only_constraint", g3, frag_from(win, 1000, 260), ip);
            // INVARIANT-ONLY: both mates inside the left flank, constraining neither allele.
            run("invariant_only_seed", g3, frag_from(win, 20, 260), ip);
            // REVERSE STRAND, as the library actually produces it: mate 1 is the reverse-
            // complement of the DOWNSTREAM end and mate 2 the forward UPSTREAM end, which is the
            // other valid-FR orientation. Swapping and complementing both mates instead gives a
            // fragment that is not FR-valid at all -- both arms then agree on ZERO cells, which
            // compares nothing.
            {
                Fragment r;
                r.name = "rev";
                r.r1 = reverse_complement(win.substr(850 + 300 - 100, 100));
                r.r2 = win.substr(850, 100);
                run("reverse_strand", g3, r, ip);
            }
            // MORE THAN EIGHT IDENTICAL ORIGINS: a repeat placed nine times in the context.
            {
                const std::string unit = rseq(120);
                std::string rep;
                for (int i = 0; i < 9; ++i) rep += unit;
                const LinkageGeometry gr = geom_of(A, Bv, rep, L, R);
                const std::string wr = L + A[0] + rep + Bv[0] + R;
                run("nine_identical_origins", gr, frag_from(wr, 950, 240), ip);
            }
            // NON-ACGT: must take the exhaustive fallback over the virtual allele product.
            {
                Fragment f = frag_from(win, 850, 300);
                f.r1[10] = 'N';
                run("non_acgt_fallback", g3, f, ip);
            }
            // ZERO-STATE: a read from unrelated sequence. Every cell must still be emitted, as -inf.
            {
                Fragment f;
                f.name = "z"; f.r1 = rseq(100); f.r2 = rseq(100);
                run("zero_state_cells", g3, f, ip);
            }
            // AT C4 SCALE, and with the DERIVED flank of zero that C4 actually has -- so an
            // invariant-only seed cannot occur and the reduction is the one production would see.
            {
                std::vector<std::string> A2v, B2v;
                for (int i = 0; i < 118; ++i) A2v.push_back(rseq(500));
                for (int i = 0; i < 119; ++i) B2v.push_back(rseq(500));
                const LinkageGeometry gl = geom_of(A2v, B2v, "", "", "");
                const std::string wl = A2v[40] + B2v[77];
                run("c4_scale_118x119", gl, frag_from(wl, 400, 300), ip);
                // A repeated allele: several A alleles sharing sequence must all be proposed.
                std::vector<std::string> A3v = A2v;
                A3v[7] = A3v[40]; A3v[91] = A3v[40];
                const LinkageGeometry gd = geom_of(A3v, B2v, "", "", "");
                run("c4_scale_dup_alleles", gd, frag_from(wl, 400, 300), ip);
            }
        }
        return 0;
    }

    // ---- TRANSACTIONAL ACTIVATION SELF-TEST -----------------------------------------------------
    // Subtracting linkage-owned fragments and activating their edges is ONE transaction. Half of it
    // makes those fragments vanish from BOTH models -- gone from the marker counts, consumed by no
    // edge -- leaving a run quietly weaker than the legacy caller it extends.
    if (activation_selftest) {
        // FOUR CONDITIONS, reported separately. Collapsing them into one "complete" is how an
        // INCOMPLETE call comes to look finished -- mapping_refused is ownership-complete and
        // factor-INcomplete, and its final status must be INCOMPLETE.
        std::printf("case\townership_complete\tfactors_buildable\thybrid_activated\tcall_status"
                    "\tactive_edges\texcluded\tconsumed\tunique\tequality\trefusal\n");
        const auto frag = [](const char* n2) { Fragment f; f.name = n2; f.r1 = "A"; f.r2 = "C";
                                               return f; };
        const auto own = [](OwnerKind k, std::uint32_t lo, std::uint32_t hi) {
            FragmentOwner o; o.kind = k; o.block_lo = lo; o.block_hi = hi;
            if (k == OwnerKind::Wide) o.panel_domain_var_scope = {lo, (lo + hi) / 2, hi};
            return o;
        };
        LinkageEdge good;
        good.n_a = 2; good.n_b = 2; good.status = LinkageStatus::Ok; good.log_psi.assign(16, 0.0);
        const AlleleMapping mok = build_allele_mapping({0, 1, 1}, 2, -1);
        const AlleleMapping mbad = build_allele_mapping({0, -1, 1}, 2, -1);   // refused
        const auto run = [&](const char* name, const std::vector<Fragment>& fr,
                             const std::vector<FragmentOwner>& ow,
                             const std::vector<EdgeStatusEntry>& es,
                             const std::map<std::pair<std::uint32_t, std::uint32_t>,
                                            LinkageEdge>& ed,
                             const std::vector<AlleleMapping>& mp) {
            const HybridActivation A = plan_hybrid_activation(fr, ow, es, ed, mp, 3);
            // THE EQUALITY: excluded == consumed by ACTIVE edges, each exactly once.
            std::set<std::string> ex(A.excluded_fragments.begin(), A.excluded_fragments.end());
            std::set<std::string> consumed;
            for (std::size_t i = 0; i < ow.size() && i < fr.size(); ++i) {
                if (ow[i].kind != OwnerKind::Linkage) continue;
                const std::uint32_t b = ow[i].block_hi;
                if (b < A.kernel_edges.size() && A.kernel_edges[b].active) consumed.insert(fr[i].name);
            }
            const bool eq = (ex == consumed);
            std::printf("%s\t%d\t%d\t%d\t%s\t%zu\t%zu\t%zu\t%zu\t%d\t%s\n", name,
                        A.ownership_complete ? 1 : 0, A.factors_buildable ? 1 : 0,
                        A.hybrid_activated ? 1 : 0, hybrid_call_status_name(A.call_status),
                        A.active_edges, A.excluded_fragments.size(), A.consumed_fragments,
                        ex.size(), eq ? 1 : 0, A.refusal.empty() ? "-" : A.refusal.c_str());
        };
        using OK = OwnerKind;  using LS = LinkageStatus;
        // No linkage evidence at all: complete, activated, nothing excluded -> legacy inference.
        run("no_linkage", {frag("u1"), frag("u2")}, {own(OK::Unary,1,1), own(OK::Unary,2,2)},
            {}, {}, {mok, mok, mok});
        // A complete hybrid: two fragments consumed by one active edge, both excluded.
        run("complete", {frag("u1"), frag("l1"), frag("l2")},
            {own(OK::Unary,1,1), own(OK::Linkage,1,2), own(OK::Linkage,1,2)},
            {{1,2,LS::Ok,2}}, {{{1,2}, good}}, {mok, mok, mok});
        // A REFUSED EDGE: nothing excluded, nothing activated. Its fragments stay in the markers.
        run("refused_edge", {frag("u1"), frag("l1")},
            {own(OK::Unary,1,1), own(OK::Linkage,1,2)},
            {{1,2,LS::TooManyConfigurations,1}}, {}, {mok, mok, mok});
        // A WIDE fragment: no pairwise consumer, so the whole activation is refused.
        run("wide_present", {frag("u1"), frag("w1")},
            {own(OK::Unary,1,1), own(OK::Wide,0,2)}, {}, {}, {mok, mok, mok});
        // Completeness passes but a MAPPING is refused: the edge cannot be built, so the
        // transaction is abandoned rather than run with a hole.
        run("mapping_refused", {frag("l1")}, {own(OK::Linkage,1,2)},
            {{1,2,LS::Ok,1}}, {{{1,2}, good}}, {mok, mbad, mok});
        // An UNUSABLE fragment: explicitly missing evidence, so no activation.
        run("unusable_present", {frag("u1"), frag("x1")},
            {own(OK::Unary,1,1), own(OK::Unusable,0,0)}, {}, {}, {mok, mok, mok});
        return 0;
    }

    // ---- HYBRID COMPLETENESS SELF-TEST ----------------------------------------------------------
    // The third situation is the one that gets lost: linkage or Wide evidence EXISTS but cannot be
    // represented. An inactive ChainEdgeLinkage cannot express it -- the kernel reads that exactly
    // like a legitimately linkage-free edge -- so a refused edge would silently become neutral and
    // the posterior would look complete while a reduced evidence model ran.
    if (completeness_selftest) {
        std::printf("case\townership_complete\towned\tinvariant\tunary\tlinkage\twide"
                    "\trefused_edge\tunusable\tn_refusals\treasons\twide_scopes\n");
        const auto mkown = [](OwnerKind k, std::uint32_t lo, std::uint32_t hi,
                              std::vector<std::uint32_t> vs = {}) {
            FragmentOwner o;
            o.kind = k; o.block_lo = lo; o.block_hi = hi; o.panel_domain_var_scope = std::move(vs);
            return o;
        };
        const auto run = [&](const char* name, const std::vector<FragmentOwner>& owners,
                             const std::vector<EdgeStatusEntry>& edges) {
            const HybridCompletenessReport R = assess_hybrid_completeness(owners, edges);
            std::string reasons;
            for (const EdgeRefusal& r : R.refusals) {
                reasons += (reasons.empty() ? "" : ";") + std::to_string(r.block_a) + "-" +
                           std::to_string(r.block_b) + ":" + linkage_status_name(r.status) +
                           "(" + std::to_string(r.n_fragments) + ")";
            }
            if (reasons.empty()) reasons = "-";
            std::string scopes;
            for (const auto& sc : R.wide_scopes) {
                std::string one;
                for (std::uint32_t b : sc) one += (one.empty() ? "" : ".") + std::to_string(b);
                scopes += (scopes.empty() ? "" : ";") + one;
            }
            if (scopes.empty()) scopes = "-";
            std::printf("%s\t%d\t%zu\t%zu\t%zu\t%zu\t%zu\t%zu\t%zu\t%zu\t%s\t%s\n",
                        name, R.ownership_complete ? 1 : 0, R.owned_total, R.invariant,
                        R.consumed_unary,
                        R.consumed_linkage, R.unconsumed_wide, R.unconsumed_refused_edge,
                        R.unconsumed_unusable, R.refusals.size(), reasons.c_str(), scopes.c_str());
        };
        using OK = OwnerKind;  using LS = LinkageStatus;
        // 1. No owned linkage at all: pure legacy shape, complete.
        run("no_linkage", {mkown(OK::Unary,1,1), mkown(OK::Unary,2,2), mkown(OK::Invariant,0,0)}, {});
        // 2. Linkage present and representable: complete.
        run("linkage_ok", {mkown(OK::Unary,1,1), mkown(OK::Linkage,1,2)},
            {{1,2,LS::Ok,1}});
        // 3-6. Each refusal reason leaves its fragments without a consumer -> INCOMPLETE.
        run("too_many_configs", {mkown(OK::Linkage,1,2)}, {{1,2,LS::TooManyConfigurations,1}});
        run("exposure_refused", {mkown(OK::Linkage,1,2)}, {{1,2,LS::ExposureDoesNotCancel,1}});
        run("invalid_emission", {mkown(OK::Linkage,1,2)}, {{1,2,LS::InvalidEmissions,1}});
        run("mapping_refused",  {mkown(OK::Linkage,1,2)}, {{1,2,LS::NotComputed,1}});
        // 7. A Wide fragment has no pairwise consumer; its VARIABLE scope is reported.
        run("wide_fragment", {mkown(OK::Unary,1,1), mkown(OK::Wide,1,3,{1,2,3})}, {});
        // 8. Several refusals, for DIFFERENT reasons: all must survive, not only the last.
        run("multi_refusal",
            {mkown(OK::Linkage,1,2), mkown(OK::Linkage,3,4), mkown(OK::Linkage,5,6)},
            {{1,2,LS::TooManyConfigurations,1}, {3,4,LS::ExposureDoesNotCancel,1},
             {5,6,LS::InvalidEmissions,1}});
        // 9. An Unusable fragment is explicitly reported missing evidence, not silence.
        run("unusable_fragment", {mkown(OK::Unary,1,1), mkown(OK::Unusable,0,0)}, {});
        // 10. Invariant fragments need no consumer and must not make the run incomplete.
        run("invariant_only", {mkown(OK::Invariant,0,0), mkown(OK::Invariant,0,0)}, {});
        // NO REFUSED EDGE TABLE REACHES THE KERNEL. make_kernel_edge() is the one construction
        // path, and it yields an INACTIVE entry with no table for any refused edge or refused
        // mapping -- because the kernel cannot tell a refused edge from a linkage-free one.
        std::printf("kernel_edge_case\tactive\tn_a\tn_b\tpsi_size\n");
        {
            LinkageEdge good;
            good.n_a = 2; good.n_b = 2; good.status = LinkageStatus::Ok;
            good.log_psi.assign(16, 0.0);
            const AlleleMapping ma = build_allele_mapping({0, 1, 1}, 2, -1);
            const AlleleMapping mb = build_allele_mapping({0, 0, 1}, 2, -1);
            const auto emit = [&](const char* nm, const ChainEdgeLinkage& k) {
                std::printf("%s\t%d\t%zu\t%zu\t%zu\n", nm, k.active ? 1 : 0, k.n_a, k.n_b,
                            k.log_psi.size());
            };
            emit("usable_edge", make_kernel_edge(good, ma, mb));
            for (const auto& bad : {LinkageStatus::TooManyConfigurations,
                                    LinkageStatus::ExposureDoesNotCancel,
                                    LinkageStatus::InvalidEmissions,
                                    LinkageStatus::NotComputed}) {
                LinkageEdge e = good;
                e.status = bad;
                emit(linkage_status_name(bad), make_kernel_edge(e, ma, mb));
            }
            emit("refused_mapping_a",
                 make_kernel_edge(good, build_allele_mapping({0, -1, 1}, 2, -1), mb));
            emit("refused_mapping_b",
                 make_kernel_edge(good, ma, build_allele_mapping({0, 5, 1}, 2, -1)));
            LinkageEdge mismatched = good;
            mismatched.n_a = 3;   // allele count disagrees with the mapping
            emit("allele_count_mismatch", make_kernel_edge(mismatched, ma, mb));
        }
        return 0;
    }

    // ---- HAPLOTYPE -> ALLELE MAPPING SELF-TEST --------------------------------------------------
    // The int -> unsigned boundary is where this breaks silently: BlockAlleles reports -1 for "no
    // allele here", and a bare cast makes that 4294967295, which then indexes log_psi far out of
    // bounds. Every case below is checked at the boundary rather than downstream.
    if (mapping_selftest) {
        std::printf("case\tstatus\tn_alleles\tbypass_resolved\tbad_hap\tbad_value\tmapping"
                    "\tvector_size\n");
        const auto show = [&](const char* name, const std::vector<int>& av,
                              std::size_t na, int bypass) {
            const AlleleMapping m = build_allele_mapping(av, na, bypass);
            std::string mp;
            if (m.status == MappingStatus::Ok) {
                for (std::size_t k = 0; k < m.allele.size(); ++k) {
                    mp += (k ? "," : "") + std::to_string(m.allele[k]);
                }
            } else {
                mp = "-";
            }
            std::printf("%s\t%s\t%zu\t%zu\t%zu\t%d\t%s\t%zu\n", name,
                        mapping_status_name(m.status), m.n_alleles, m.n_bypass_resolved,
                        m.first_bad_haplotype, m.first_bad_value, mp.c_str(), m.allele.size());
        };
        // A NORMAL mapping, and a deliberately PERMUTED one: 3 haplotypes over 2 alleles, where the
        // haplotype index is not the allele index anywhere.
        show("normal",        {0, 1, 0},    2, -1);
        show("permuted",      {1, 0, 1},    2, -1);
        // A BYPASSING haplotype resolves to the block's bypass allele, not to -1 and not to 0.
        show("bypass",        {0, -1, 1},   3,  2);
        // NO bypass allele to resolve to: refused, never cast.
        show("missing",       {0, -1, 1},   2, -1);
        // A large negative must not wrap; it is refused at the same gate.
        show("missing_large", {0, -999, 1}, 2, -1);
        // An allele index past the block's allele count is refused before it can index anything.
        show("out_of_range",  {0, 5, 1},    2, -1);
        // A bypass index that is itself out of range must not be trusted either.
        show("bad_bypass",    {0, -1, 1},   2,  7);
        return 0;
    }

    // ---- HYBRID CHAIN vs BRUTE-FORCE PATH ORACLE ------------------------------------------------
    // The forward-backward recursion is checked against an independent enumeration of EVERY ordered
    // diploid state path. That is the only comparison that catches linkage applied at the wrong
    // edge, either factor applied twice, an accidental transition row-normalisation, an
    // ordered-state mapping error, or a correct best call resting on wrong posterior mass -- a
    // best-path comparison catches none of them.
    //
    // THE FIXTURE USES 3 HAPLOTYPES BUT 2 ALLELES PER BLOCK, with a deliberately non-identity map
    // (allele_a = [0,1,1], allele_b = [0,0,1]). If anything ever treats the haplotype index as an
    // allele index, the tables are the wrong shape and the oracle disagrees immediately.
    if (hybrid_oracle) {
        const std::size_t nh = 3, nb = 3, ns = nh * nh;
        HybridChain ch;
        ch.n_hap = nh; ch.n_blocks = nb; ch.recomb = 0.10;
        // Deterministic, asymmetric emissions: nothing here may be symmetric by accident, or a
        // mapping error could cancel itself.
        ch.log_emission.assign(nb, std::vector<double>(ns, 0.0));
        for (std::size_t b = 0; b < nb; ++b) {
            for (std::size_t i = 0; i < nh; ++i)
            for (std::size_t j = 0; j < nh; ++j) {
                ch.log_emission[b][i * nh + j] =
                    -0.37 * static_cast<double>((b + 1) * (i + 1)) -
                     0.11 * static_cast<double>((j + 2) * (b + 3)) -
                     0.05 * static_cast<double>(i * j);
            }
        }
        ch.edges.assign(nb, HybridEdge{});
        // edges[1]: LINKAGE-FREE. Its transition must remain exactly Li-Stephens.
        ch.edges[1].has_linkage = false;
        // edges[2]: an INFORMATIVE linkage edge, built through the real aggregation path so the
        // oracle sees the same mean-one table the chain does.
        {
            LinkageGeometry g;
            g.block_a = 1; g.block_b = 2;
            g.alleles_a = {std::string(600, 'A'), std::string(600, 'C')};
            g.alleles_b = {std::string(600, 'G'), std::string(600, 'T')};
            g.context = std::string(40, 'A');
            g.lflank = std::string(600, 'A');
            g.rflank = std::string(600, 'A');
            g.window_len.assign(4, 2440);
            g.exposure.assign(4, 1000.0);
            g.exposure_affine = true;
            g.ok = true;
            LinkageEmission m;
            m.n_a = 2; m.n_b = 2; m.log_p_bg = -400.0; m.ok = true; m.informative = true;
            m.mass = {-100.0, -137.0, -142.0, -103.0};
            const LinkageEdge E = aggregate_linkage_edge({m}, g, 0.05,
                                                         std::log1p(-0.05), std::log(0.05));
            if (!E.usable()) { std::fprintf(stderr, "oracle: linkage edge unusable\n"); return 2; }
            ch.edges[2].has_linkage = true;
            ch.edges[2].n_a = E.n_a; ch.edges[2].n_b = E.n_b;
            ch.edges[2].log_psi = E.log_psi;
            ch.edges[2].allele_a = {0, 1, 1};   // 3 haplotypes -> 2 alleles, NOT the identity
            ch.edges[2].allele_b = {0, 0, 1};
        }
        const HybridPosterior fb = hybrid_forward_backward(ch);
        const HybridPosterior bf = hybrid_bruteforce(ch);
        if (!fb.ok || !bf.ok) { std::fprintf(stderr, "oracle: a run failed\n"); return 2; }
        double worst_marg = 0.0;
        for (std::size_t b = 0; b < nb; ++b)
            for (std::size_t s = 0; s < ns; ++s)
                worst_marg = std::max(worst_marg,
                                      std::abs(std::exp(fb.log_marginal[b][s]) -
                                               std::exp(bf.log_marginal[b][s])));
        // A linkage-FREE control: with every edge plain Li-Stephens the two must still agree, which
        // separates "the recursion is right" from "the linkage table is right".
        HybridChain plain = ch;
        for (auto& e : plain.edges) { e.has_linkage = false; }
        const HybridPosterior fb0 = hybrid_forward_backward(plain);
        const HybridPosterior bf0 = hybrid_bruteforce(plain);
        double worst_marg0 = 0.0;
        for (std::size_t b = 0; b < nb; ++b)
            for (std::size_t s = 0; s < ns; ++s)
                worst_marg0 = std::max(worst_marg0,
                                       std::abs(std::exp(fb0.log_marginal[b][s]) -
                                                std::exp(bf0.log_marginal[b][s])));
        // Marginals must be distributions, and linkage must actually MOVE them, or the comparison
        // would pass on a chain where the factor does nothing.
        double worst_sum = 0.0, linkage_effect = 0.0;
        for (std::size_t b = 0; b < nb; ++b) {
            double sum = 0.0;
            for (std::size_t s = 0; s < ns; ++s) {
                sum += std::exp(fb.log_marginal[b][s]);
                linkage_effect = std::max(linkage_effect,
                                          std::abs(std::exp(fb.log_marginal[b][s]) -
                                                   std::exp(fb0.log_marginal[b][s])));
            }
            worst_sum = std::max(worst_sum, std::abs(sum - 1.0));
        }
        std::printf("metric\tvalue\n");
        std::printf("log_weight_sum_fb\t%.17g\n", fb.log_partition_unnormalised);
        std::printf("log_weight_sum_bruteforce\t%.17g\n", bf.log_partition_unnormalised);
        std::printf("log_weight_sum_abs_diff\t%.17g\n",
                    std::abs(fb.log_partition_unnormalised - bf.log_partition_unnormalised));
        std::printf("worst_marginal_abs_diff\t%.17g\n", worst_marg);
        std::printf("log_weight_sum_abs_diff_no_linkage\t%.17g\n",
                    std::abs(fb0.log_partition_unnormalised - bf0.log_partition_unnormalised));
        std::printf("worst_marginal_abs_diff_no_linkage\t%.17g\n", worst_marg0);
        std::printf("worst_marginal_sum_dev\t%.17g\n", worst_sum);
        std::printf("linkage_marginal_effect\t%.17g\n", linkage_effect);
        // BOTH KERNEL PATHS MUST HAVE RUN, in this one chain. A branch that exists but is never
        // taken is not covered, however many assertions surround it -- and the whole point of the
        // shared kernel is that the factorised and linked paths cannot drift apart unnoticed.
        std::printf("factorised_edges\t%zu\n", fb.factorised_edges);
        std::printf("linked_edges\t%zu\n", fb.linked_edges);
        std::printf("factorised_edges_no_linkage\t%zu\n", fb0.factorised_edges);
        std::printf("linked_edges_no_linkage\t%zu\n", fb0.linked_edges);

        // ---- GROUPED SPARSE vs DENSE KERNEL -----------------------------------------------
        // The same chain run two ways -- the dense n_h^4 table and the grouped sparse correction --
        // across the regimes that a single recombination rate or a single sparse class would hide:
        // r = 0 (identity only), r = 1 (uniform switch only), an intermediate r where all four
        // expanded terms are active, several non-neutral classes on one edge, several haplotypes
        // per allele, TWO linked edges in one chain, and a near-total cancellation.
        {
            std::printf("sparse_case\tr\tn_a\tn_b\tlinked_edges\tclasses\tcorrections"
                        "\tlogw_absdiff\tmarg_absdiff\tswap_absdiff\tclamped\tworst_negative"
                        "\tblock1_max_marginal\n");
            const double lam = 0.05, mx = std::log1p(-0.05), bw = std::log(0.05);
            const auto build_geom = [&](std::size_t na, std::size_t nb) {
                LinkageGeometry g;
                g.block_a = 1; g.block_b = 2;
                g.alleles_a.assign(na, std::string(600, 'A'));
                g.alleles_b.assign(nb, std::string(600, 'C'));
                g.context = std::string(40, 'A');
                g.lflank = std::string(600, 'A');
                g.rflank = std::string(600, 'A');
                g.window_len.assign(na * nb, 2440);
                g.exposure.assign(na * nb, 1000.0);
                g.exposure_affine = true;
                g.ok = true;
                return g;
            };
            const auto to_sparse = [&](const SparseLinkageEdge& SE,
                                       const std::vector<std::uint32_t>& aa,
                                       const std::vector<std::uint32_t>& ab) {
                SparseEdgeLinkage se;
                se.active = true; se.n_a = SE.n_a; se.n_b = SE.n_b;
                se.allele_a = aa; se.allele_b = ab;
                se.group_a.assign(SE.n_a, {});
                se.group_b.assign(SE.n_b, {});
                for (std::uint32_t h = 0; h < aa.size(); ++h) {
                    se.group_a[aa[h]].push_back(h);
                    se.group_b[ab[h]].push_back(h);
                }
                for (std::size_t a1 = 0; a1 < SE.n_a; ++a1)
                for (std::size_t a2 = a1 + 1; a2 < SE.n_a; ++a2)
                for (std::size_t b1 = 0; b1 < SE.n_b; ++b1)
                for (std::size_t b2 = b1 + 1; b2 < SE.n_b; ++b2) {
                    const double ls = SE.log_psi(a1, b1, a2, b2);
                    const double lc = SE.log_psi(a1, b2, a2, b1);
                    if (ls == 0.0 && lc == 0.0) continue;
                    SparsePhaseClass c;
                    c.amin = static_cast<std::uint32_t>(a1); c.amax = static_cast<std::uint32_t>(a2);
                    c.bmin = static_cast<std::uint32_t>(b1); c.bmax = static_cast<std::uint32_t>(b2);
                    c.straight_m1 = std::expm1(ls);   // expm1: a psi near one keeps its correction
                    c.crossed_m1 = std::expm1(lc);
                    se.classes.push_back(c);
                }
                return se;
            };
            const auto run_case =
                [&](const char* name, double r, std::size_t na, std::size_t nb,
                    const std::vector<std::uint32_t>& aa, const std::vector<std::uint32_t>& ab,
                    const std::vector<double>& mass, bool two_edges,
                    const std::vector<std::size_t>& favour, bool symmetric_emissions = false) {
                const std::size_t nh2 = aa.size(), nb2 = 3, ns2 = nh2 * nh2;
                const LinkageGeometry g = build_geom(na, nb);
                LinkageEmission m;
                m.n_a = na; m.n_b = nb; m.log_p_bg = -400.0; m.ok = true; m.informative = true;
                m.mass = mass;
                const LinkageEdge DE = aggregate_linkage_edge({m}, g, lam, mx, bw);
                const SparseLinkageEdge SE = build_sparse_linkage_edge({m}, g, lam, mx, bw);
                if (!DE.usable() || !SE.usable()) {
                    std::printf("%s\tunusable\n", name); return;
                }
                std::vector<std::vector<double>> em(nb2, std::vector<double>(ns2, 0.0));
                for (std::size_t b = 0; b < nb2; ++b)
                    for (std::size_t i = 0; i < nh2; ++i)
                        for (std::size_t j = 0; j < nh2; ++j)
                            em[b][i * nh2 + j] = std::exp(-0.37 * (b + 1) * (i + 1) -
                                                          0.11 * (j + 2) * (b + 3) - 0.05 * i * j);
                // SWAP INVARIANCE IS A PROPERTY OF THE KERNEL, and can only be tested against
                // emissions that are themselves symmetric. The default emissions above are
                // deliberately asymmetric in i and j -- good for catching index errors, useless for
                // this -- and a swap check against them fails even with NO linkage.
                if (symmetric_emissions) {
                    for (std::size_t b = 0; b < nb2; ++b)
                        for (std::size_t i = 0; i < nh2; ++i)
                            for (std::size_t j = 0; j < nh2; ++j) {
                                const double v = std::exp(-0.21 * (b + 1) * (i + 1) * (j + 1) -
                                                          0.05 * (i + j));
                                em[b][i * nh2 + j] = v;
                                em[b][j * nh2 + i] = v;
                            }
                }
                if (!favour.empty()) {
                    for (std::size_t b = 1; b < nb2; ++b)
                        for (std::size_t k = 0; k < ns2; ++k) em[b][k] = 1e-12;
                    for (std::size_t k : favour)
                        if (k < ns2) { em[1][k] = 1.0; em[2][k] = 1.0; }
                }
                std::vector<ChainEdgeLinkage> de(nb2);
                std::vector<SparseEdgeLinkage> se(nb2);
                const SparseEdgeLinkage sp1 = to_sparse(SE, aa, ab);
                for (std::size_t eb : two_edges ? std::vector<std::size_t>{1, 2}
                                                : std::vector<std::size_t>{2}) {
                    de[eb].active = true; de[eb].n_a = DE.n_a; de[eb].n_b = DE.n_b;
                    de[eb].allele_a = aa; de[eb].allele_b = ab; de[eb].log_psi = DE.log_psi;
                    se[eb] = sp1;
                }
                const auto emit_fn = [&](std::size_t b, std::vector<double>& ev) { ev = em[b]; };
                std::vector<std::vector<double>> f1, b1v, f2, b2v;
                ChainKernelStats s1, s2;
                chain_forward_backward(nh2, nb2, r, emit_fn, &de, f1, b1v, &s1, nullptr);
                chain_forward_backward(nh2, nb2, r, emit_fn, nullptr, f2, b2v, &s2, &se);
                double worst = 0.0, swap_diff = 0.0, mx_marg = 0.0;
                for (std::size_t b = 0; b < nb2; ++b) {
                    double z1 = 0.0, z2 = 0.0;
                    for (std::size_t k = 0; k < ns2; ++k) {
                        z1 += f1[b][k] * b1v[b][k]; z2 += f2[b][k] * b2v[b][k];
                    }
                    for (std::size_t i = 0; i < nh2; ++i)
                    for (std::size_t j = 0; j < nh2; ++j) {
                        const std::size_t k = i * nh2 + j, ks = j * nh2 + i;
                        const double p1 = z1 > 0 ? f1[b][k] * b1v[b][k] / z1 : 0.0;
                        const double p2 = z2 > 0 ? f2[b][k] * b2v[b][k] / z2 : 0.0;
                        const double p2s = z2 > 0 ? f2[b][ks] * b2v[b][ks] / z2 : 0.0;
                        worst = std::max(worst, std::abs(p1 - p2));
                        // HOMOLOGUE-SWAP INVARIANCE of the posterior itself: (i,j) and (j,i) are one
                        // diploid state written twice and must carry equal mass.
                        swap_diff = std::max(swap_diff, std::abs(p2 - p2s));
                        if (b == 1) mx_marg = std::max(mx_marg, p2);
                    }
                }
                std::printf("%s\t%.2f\t%zu\t%zu\t%zu\t%zu\t%zu\t%.17g\t%.17g\t%.17g"
                            "\t%zu\t%.3g\t%.10f\n",
                            name, r, na, nb, two_edges ? 2u : 1u, sp1.classes.size(),
                            s2.corrections_applied,
                            std::abs(s1.log_weight_sum - s2.log_weight_sum), worst, swap_diff,
                            s2.clamped_negatives, s2.worst_negative, mx_marg);
            };
            const std::vector<std::uint32_t> aa2 = {0, 1, 1}, ab2 = {0, 0, 1};
            const std::vector<double> m2 = {-100.0, -137.0, -142.0, -103.0};
            // r = 0 leaves only the identity term; r = 1 only the uniform switch; 0.10 activates all
            // four. A single rate cannot distinguish a wrong term from a right one.
            run_case("r0_identity_only",   0.0,  2, 2, aa2, ab2, m2, false, {});
            run_case("r1_uniform_switch",  1.0,  2, 2, aa2, ab2, m2, false, {});
            run_case("r_intermediate",     0.10, 2, 2, aa2, ab2, m2, false, {});
            run_case("two_linked_edges",   0.10, 2, 2, aa2, ab2, m2, true,  {});
            run_case("no_linkage_control", 0.10, 2, 2, aa2, ab2,
                     {-400.0, -400.0, -400.0, -400.0}, false, {});
            // Several haplotypes per allele, and 3x3 alleles so SEVERAL het x het classes exist on
            // one edge rather than the single class a 2x2 edge can hold.
            const std::vector<std::uint32_t> aa3 = {0, 0, 1, 2, 2}, ab3 = {0, 1, 1, 2, 2};
            std::vector<double> m3(9, -std::numeric_limits<double>::infinity());
            m3[0 * 3 + 0] = -100.0; m3[1 * 3 + 1] = -104.0; m3[2 * 3 + 2] = -101.0;
            m3[0 * 3 + 2] = -150.0; m3[2 * 3 + 0] = -155.0;
            run_case("multi_class_3x3",    0.10, 3, 3, aa3, ab3, m3, false, {});
            run_case("multi_class_r0",     0.0,  3, 3, aa3, ab3, m3, false, {});
            run_case("multi_class_r1",     1.0,  3, 3, aa3, ab3, m3, false, {});
            run_case("adversarial_mass_in_zero_phase", 0.10, 2, 2, aa2, ab2,
                     {-100.0, -1000.0, -1000.0, -100.0}, false, {1, 2, 6, 7});
            // SYMMETRIC-EMISSION arms, the only ones where swap invariance is a statement about the
            // kernel rather than about the fixture.
            run_case("swapsym_no_linkage",  0.10, 2, 2, aa2, ab2,
                     {-400.0, -400.0, -400.0, -400.0}, false, {}, true);
            run_case("swapsym_linked",      0.10, 2, 2, aa2, ab2, m2, false, {}, true);
            run_case("swapsym_multi_class", 0.10, 3, 3, aa3, ab3, m3, false, {}, true);
        }
        return 0;
    }

    // ---- LINKAGE FACTOR SELF-TEST ---------------------------------------------------------------
    // The neutrality property cannot be produced reliably from a read fixture: it needs an edge with
    // NO fragments, and one whose emissions are identical across every configuration. Both are
    // constructed here directly, over content classes of all three cardinalities (1, 2 and 4), which
    // is exactly where a sum-one normalisation would have leaked a -log|C| content penalty.
    if (linkage_selftest) {
        LinkageGeometry g;
        g.block_a = 0; g.block_b = 1;
        // EQUAL-LENGTH alleles, so the edge's exposure is identical across configurations and the
        // test isolates the phase factor from the exposure term.
        g.alleles_a = {std::string(600, 'A'), std::string(600, 'C')};
        g.alleles_b = {std::string(600, 'G'), std::string(600, 'T')};
        g.context = std::string(40, 'A');
        g.lflank = std::string(600, 'A');
        g.rflank = std::string(600, 'A');
        g.window_len.assign(4, 2440);
        g.exposure.assign(4, 1000.0);
        g.exposure_affine = true;
        g.ok = true;
        const double lam = 0.05, mix = std::log1p(-0.05), bgw = std::log(0.05);
        const auto report = [&](const char* name, const std::vector<LinkageEmission>& ems,
                                const LinkageGeometry* gover = nullptr, double lam_over = 0.0) {
            const LinkageEdge E = aggregate_linkage_edge(ems, gover ? *gover : g,
                                                         lam_over > 0.0 ? lam_over : lam,
                                                         mix, bgw);
            // GUARD BEFORE INDEXING. A refused edge has EMPTY score/log_psi with n_a and n_b still
            // set, so looping over n_a*n_b reads out of bounds.
            if (!E.usable()) {
                std::printf("%s\t%zu\t%s\t%s\t%s\t%s\t%d\t%zu\t%s\t%.17g"
                            "\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
                            name, ems.size(), "0", "0", "0", "-", 0, E.n_invalid,
                            linkage_status_name(E.status), E.exposure_asymmetry,
                            "0", "0", "0", "0", "0", "0", "0", "0", "0");
                return;
            }
            const std::size_t na = E.n_a, nb = E.n_b;
            double maxabs = 0.0, worst_mean = 0.0, swap = 0.0;
            std::map<std::size_t, std::size_t> class_sizes;
            std::vector<char> seen(E.log_psi.size(), 0);
            for (std::size_t a1 = 0; a1 < na; ++a1)
            for (std::size_t b1 = 0; b1 < nb; ++b1)
            for (std::size_t a2 = 0; a2 < na; ++a2)
            for (std::size_t b2 = 0; b2 < nb; ++b2) {
                const std::size_t c = ((a1*nb + b1)*na + a2)*nb + b2;
                maxabs = std::max(maxabs, std::abs(E.log_psi[c]));
                const std::size_t cs = ((a2*nb + b2)*na + a1)*nb + b1;
                swap = std::max(swap, std::abs(E.log_psi[c] - E.log_psi[cs]));
                if (seen[c]) continue;
                const std::size_t A2[2] = {a1, a2}, B2[2] = {b1, b2};
                std::vector<std::size_t> cls;
                for (int q1 = 0; q1 < 2; ++q1) for (int q2 = 0; q2 < 2; ++q2)
                    cls.push_back(((A2[q1]*nb + B2[q2])*na + A2[1-q1])*nb + B2[1-q2]);
                std::sort(cls.begin(), cls.end());
                cls.erase(std::unique(cls.begin(), cls.end()), cls.end());
                double m = 0.0;
                for (std::size_t z : cls) { m += std::exp(E.log_psi[z]); seen[z] = 1; }
                m /= static_cast<double>(cls.size());
                worst_mean = std::max(worst_mean, std::abs(m - 1.0));
                ++class_sizes[cls.size()];
            }
            std::string sizes;
            for (const auto& kv : class_sizes) {
                sizes += (sizes.empty() ? "" : ",") + std::to_string(kv.first) + "x" +
                         std::to_string(kv.second);
            }
            // THE MEAN-ONE INVARIANT, stated as a bound rather than assumed. psi averages to one
            // over a class of size |C|, so no configuration can exceed |C| and
            //     max log psi <= log|C| <= log 4 = 1.3863.
            // The large magnitudes a linkage edge produces are all NEGATIVE -- losing
            // configurations underflowing toward zero -- so exp(log psi) cannot overflow. Measured
            // here rather than reasoned about, with non-finite values counted separately.
            double max_log = -std::numeric_limits<double>::infinity();
            double min_log = std::numeric_limits<double>::infinity();
            std::size_t nonfinite = 0, over_bound = 0;
            {
                std::vector<char> seen2(E.log_psi.size(), 0);
                for (std::size_t a1 = 0; a1 < na; ++a1)
                for (std::size_t b1 = 0; b1 < nb; ++b1)
                for (std::size_t a2 = 0; a2 < na; ++a2)
                for (std::size_t b2 = 0; b2 < nb; ++b2) {
                    const std::size_t c = ((a1*nb + b1)*na + a2)*nb + b2;
                    const double v = E.log_psi[c];
                    if (!std::isfinite(v)) ++nonfinite;
                    else { max_log = std::max(max_log, v); min_log = std::min(min_log, v); }
                    if (seen2[c]) continue;
                    const std::size_t A3[2] = {a1, a2}, B3[2] = {b1, b2};
                    std::vector<std::size_t> cl;
                    for (int q1 = 0; q1 < 2; ++q1) for (int q2 = 0; q2 < 2; ++q2)
                        cl.push_back(((A3[q1]*nb + B3[q2])*na + A3[1-q1])*nb + B3[1-q2]);
                    std::sort(cl.begin(), cl.end());
                    cl.erase(std::unique(cl.begin(), cl.end()), cl.end());
                    const double bound = std::log(static_cast<double>(cl.size()));
                    for (std::size_t z : cl) {
                        seen2[z] = 1;
                        if (std::isfinite(E.log_psi[z]) && E.log_psi[z] > bound + 1e-9) ++over_bound;
                    }
                }
            }
            // SPARSE vs DENSE, over EVERY configuration. The sparse form is exact, not an
            // approximation, so any difference at all is a defect -- including at the extremely
            // negative losing phases, where the reconstruction must not lose the branch.
            const SparseLinkageEdge SP = build_sparse_linkage_edge(ems, gover ? *gover : g,
                                                                   lam_over > 0.0 ? lam_over : lam,
                                                                   mix, bgw);
            double sparse_diff = 0.0;
            std::size_t sparse_checked = 0;
            if (SP.usable()) {
                for (std::size_t a1 = 0; a1 < na; ++a1)
                for (std::size_t b1 = 0; b1 < nb; ++b1)
                for (std::size_t a2 = 0; a2 < na; ++a2)
                for (std::size_t b2 = 0; b2 < nb; ++b2) {
                    const std::size_t c = ((a1 * nb + b1) * na + a2) * nb + b2;
                    sparse_diff = std::max(sparse_diff,
                                           std::abs(E.log_psi[c] - SP.log_psi(a1, b1, a2, b2)));
                    ++sparse_checked;
                }
            }
            std::printf("%s\t%zu\t%.17g\t%.17g\t%.17g\t%s\t%d\t%zu\t%s\t%.17g"
                        "\t%.17g\t%.17g\t%zu\t%zu\t%.17g\t%zu\t%zu\t%zu\t%zu\n",
                        name, ems.size(), maxabs, worst_mean, swap, sizes.c_str(),
                        E.usable() ? 1 : 0, E.n_invalid, linkage_status_name(E.status),
                        E.exposure_asymmetry, max_log, min_log, nonfinite, over_bound,
                        sparse_diff, sparse_checked, SP.stored_classes, SP.theoretical_configs,
                        SP.support_cells);
        };
        const auto mk = [&](const std::vector<double>& mass) {
            LinkageEmission m;
            m.n_a = 2; m.n_b = 2; m.mass = mass; m.log_p_bg = -400.0; m.ok = true;
            return m;
        };
        std::printf("case\tfragments\tmax_abs_log_psi\tworst_mean_dev\tswap_asym\tclass_sizes"
                    "\tusable\tn_invalid\tstatus\texposure_asym"
                    "\tmax_log_psi\tmin_log_psi\tnonfinite\tover_class_bound"
                    "\tsparse_vs_dense\tconfigs_checked\tstored_classes\ttheoretical_configs"
                    "\tsupport_cells\n");
        report("zero_fragments", {});
        report("flat_emissions", {mk({-100.0, -100.0, -100.0, -100.0}),
                                  mk({-100.0, -100.0, -100.0, -100.0})});
        report("all_unplaced",   {mk({-std::numeric_limits<double>::infinity(),
                                      -std::numeric_limits<double>::infinity(),
                                      -std::numeric_limits<double>::infinity(),
                                      -std::numeric_limits<double>::infinity()})});
        report("informative",    {mk({-100.0, -140.0, -140.0, -100.0})});
        // EXPOSURE THAT DOES NOT CANCEL. The equal-length alleles above make exposure
        // phase-invariant, so they CANNOT detect an edge whose exposure differs between phases --
        // which is exactly what unequal short alleles produce once a window drops below the insert
        // support. Set directly here so the check is tested rather than the fixture's luck:
        // E(0,0)+E(1,1) = 0+9 = 9 against E(0,1)+E(1,0) = 1+4 = 5.
        {
            LinkageGeometry gu = g;
            gu.exposure = {0.0, 1.0, 4.0, 9.0};
            gu.exposure_affine = false;
            report("unequal_exposure", {}, &gu);
            report("unequal_exposure_with_frags", {mk({-100.0, -140.0, -140.0, -100.0})}, &gu);
        }
        // AN OWNED FRAGMENT WITH NO EMISSION must make the edge INCOMPLETE, never vanish.
        {
            LinkageEmission bad_em;   // ok == false
            report("invalid_emission", {mk({-100.0, -140.0, -140.0, -100.0}), bad_em});
        }
        // LAMBDA MUST REACH THE FACTORS. It cancels out wherever every configuration sits far
        // above or far below the background -- there it is a constant per fragment and the mean-one
        // centering removes it, which is why the saturated phase fixture is insensitive to it. Near
        // the crossover it must move psi, and if it does not, the parameter is decoration.
        {
            const auto near_floor = [&](const std::vector<double>& mass) {
                LinkageEmission m;
                m.n_a = 2; m.n_b = 2; m.mass = mass; m.log_p_bg = -400.0; m.ok = true;
                m.informative = true;
                return m;
            };
            // Masses straddling the -400 background, so the mixture is genuinely in transition.
            const std::vector<LinkageEmission> e = {near_floor({-398.0, -402.0, -402.0, -398.0})};
            report("lambda_crossover_lo", e, nullptr, 0.001);
            report("lambda_crossover_hi", e, nullptr, 1.000);
        }
        // A C4-SCALE EDGE. The dense table refuses these outright (118 x 119 is 197 million
        // ordered configurations, 1.6 GB); the sparse form stores one contrast per class that any
        // fragment's support can actually reach. Nothing is truncated and no score enters -- the
        // classes that are absent are provably neutral.
        {
            std::printf("large_case\tn_a\tn_b\ttheoretical_configs\tstored_classes"
                        "\tsupport_cells\tstatus\tfragments\tsparse_bytes\tdense_gb"
                        "\tdense_status_at_tiny_cap\n");
            for (const auto& dims : {std::pair<std::size_t, std::size_t>{118, 119},
                                     std::pair<std::size_t, std::size_t>{457, 410}}) {
                LinkageGeometry gl;
                gl.block_a = 0; gl.block_b = 1;
                gl.alleles_a.assign(dims.first, std::string(400, 'A'));
                gl.alleles_b.assign(dims.second, std::string(400, 'C'));
                gl.context = std::string(40, 'A');
                gl.lflank.clear(); gl.rflank.clear();
                gl.window_len.assign(dims.first * dims.second, 840);
                gl.exposure.assign(dims.first * dims.second, 1000.0);
                gl.exposure_affine = true;
                gl.ok = true;
                // A handful of fragments, each placing on a few allele pairs -- which is what a real
                // fragment does: it does not place on every allele of a 118-allele block.
                std::vector<LinkageEmission> ems;
                for (int f = 0; f < 20; ++f) {
                    LinkageEmission m;
                    m.n_a = dims.first; m.n_b = dims.second;
                    m.mass.assign(dims.first * dims.second,
                                  -std::numeric_limits<double>::infinity());
                    m.log_p_bg = -400.0; m.ok = true; m.informative = true;
                    for (int k = 0; k < 3; ++k) {
                        const std::size_t al = static_cast<std::size_t>((f * 7 + k * 13) %
                                                                        dims.first);
                        const std::size_t be = static_cast<std::size_t>((f * 5 + k * 11) %
                                                                        dims.second);
                        m.mass[al * dims.second + be] = -390.0 - k;
                    }
                    ems.push_back(std::move(m));
                }
                // THE DENSE CAP MUST NOT BE CONSULTED HERE. This edge's theoretical size far
                // exceeds the dense limit, and a production build that still asked it would refuse.
                // Passing a deliberately tiny dense cap alongside proves the sparse path ignores it.
                const SparseLinkageEdge SL = build_sparse_linkage_edge(ems, gl, 0.05,
                                                                       std::log1p(-0.05),
                                                                       std::log(0.05));
                const LinkageEdge DENSE_REFUSED = aggregate_linkage_edge(
                    ems, gl, 0.05, std::log1p(-0.05), std::log(0.05), 1000u);
                // A CONSISTENT DENSE BASELINE: LinkageEdge holds TWO double tables (score and
                // log_psi), so the dense cost is 16 bytes per ordered configuration, not 8.
                const double dense_gb = static_cast<double>(SL.theoretical_configs) * 16.0 / 1e9;
                std::printf("%zux%zu\t%zu\t%zu\t%zu\t%zu\t%zu\t%s\t%zu\t%zu\t%.3f\t%s\n",
                            dims.first, dims.second, dims.first, dims.second,
                            SL.theoretical_configs, SL.stored_classes, SL.support_cells,
                            linkage_status_name(SL.status), ems.size(), SL.bytes_total(), dense_gb,
                            linkage_status_name(DENSE_REFUSED.status));
            }
        }
        // THE DENSE TABLE MUST REFUSE, NOT TRUNCATE. A locus-scale block pair (LPA: 457 x 410) is
        // 35.1 billion configurations and 561.7 GB, so the dense form cannot claim to handle all
        // loci. Refusing keeps the marker shortlist from becoming an uncertified linkage cutoff.
        {
            LinkageGeometry gb2 = g;
            gb2.alleles_a.assign(64, std::string(600, 'A'));
            gb2.alleles_b.assign(64, std::string(600, 'C'));
            gb2.window_len.assign(64 * 64, 2440);
            gb2.exposure.assign(64 * 64, 1000.0);
            report("too_many_configs", {}, &gb2);
        }
        return 0;
    }

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

        // ---- PATH -> BLOCK PROJECTION, for EVERY path -------------------------------------------
        // Every path, not a requested subset. A name list is another chance to emit one pair's
        // alleles under another's label, which is what --force-haplotypes silently did: it ADDS to
        // the shortlist rather than restricting it, so a "projection" of a named pair came back as
        // the called pair and every block compared equal. Consumers filter on exact path names.
        const std::string pbt = dump_sequences + ".path_blocks.tsv";
        const std::string pbf = dump_sequences + ".path_blocks.fa";
        std::ofstream pt(pbt), pf(pbf);
        if (!pt) throw std::runtime_error("genotype-frag: cannot write " + pbt);
        if (!pf) throw std::runtime_error("genotype-frag: cannot write " + pbf);
        pt << "group\tpath\tframe\tprojection_status\tblock\tkind\tbubble_id\twalk_begin"
              "\twalk_end\twalk_strand\tblock_bp\tblock_md5\tcanonical_md5\tcatalogue_allele"
              "\tcatalogue_representable\n";
        std::size_t n_unproj = 0, n_partial = 0, n_slices = 0;
        const auto project = [&](const char* group, const std::vector<BlockAlleles>& src,
                                 const std::string& name) {
            const auto wi = by_name.find(name);
            bool wok = false;
            std::string walk;
            if (wi != by_name.end() && wi->second != nullptr) {
                walk = spell_path_steps_sequence(graph, wi->second->steps, &wok);
            }
            // src derives the geometry; `blocks` -- the reduced calling panel -- decides
            // representability, for retained and held-out paths alike.
            const PathProjection pr =
                wok ? project_path_blocks(src, blocks, name, walk) : PathProjection{};
            const char* fr = !pr.ok ? "NA" : (pr.reverse_frame ? "rc" : "fwd");
            if (!pr.ok) {
                ++n_unproj;
                // NA rather than a best-effort segmentation. A path whose map cannot be verified
                // has no block coordinates, and inventing them puts arbitrary blocks into whatever
                // consumes this table.
                pt << group << '\t' << name << '\t' << fr << "\tunprojectable"
                   << "\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n";
                return;
            }
            if (pr.partial) ++n_partial;
            const char* st = pr.partial ? "partial" : "complete";
            for (const PathBlockSlice& sl : pr.blocks) {
                const char* kind = "backbone";
                std::size_t bub = 0;
                if (sl.block < chain.size()) {
                    kind = chain[sl.block].kind == BlockKind::Bubble ? "bubble"
                         : chain[sl.block].kind == BlockKind::Flank  ? "flank" : "backbone";
                    bub = chain[sl.block].bubble_id;
                }
                pt << group << '\t' << name << '\t' << fr << '\t' << st << '\t'
                   << sl.block << '\t' << kind << '\t' << bub << '\t'
                   << sl.walk_begin << '\t' << sl.walk_end << '\t' << (sl.reverse ? '-' : '+')
                   << '\t' << sl.seq.size() << '\t' << sl.md5 << '\t' << sl.canonical_md5 << '\t';
                if (sl.catalogue_allele < 0) pt << "NA"; else pt << sl.catalogue_allele;
                pt << '\t' << (sl.catalogue_representable ? 1 : 0) << '\n';
                pf << '>' << name << ' ' << group << " block=" << sl.block
                   << " walk=" << sl.walk_begin << '-' << sl.walk_end
                   << " strand=" << (sl.reverse ? '-' : '+') << '\n';
                for (std::size_t off = 0; off < sl.seq.size(); off += 60) {
                    pf << sl.seq.substr(off, 60) << '\n';
                }
                if (sl.seq.empty()) pf << '\n';
                ++n_slices;
            }
            // Explicit unmapped intervals: nothing is silently assigned to the nearest block.
            for (const auto& [lo, hi] : pr.unmapped) {
                pt << group << '\t' << name << '\t' << fr << "\tunmapped\tNA\tNA\tNA\t"
                   << lo << '\t' << hi << "\tNA\t" << (hi - lo) << "\tNA\tNA\tNA\t0\n";
            }
        };

        for (const PathRecord& p : panel_graph.paths) { emit("panel", blocks, p.name);
                                                        project("panel", blocks, p.name); }
        // Held-out paths are projected against the SAME fixed chain. Their block sequence can be
        // perfectly valid while catalogue_representable is 0 -- that is the interesting case, not
        // an error.
        for (const PathRecord& p : held_out) { emit("held_out", held_blocks, p.name);
                                               project("held_out", held_blocks, p.name); }
        pt.flush(); pf.flush();
        if (!pt) throw std::runtime_error("genotype-frag: write failed for " + pbt);
        if (!pf) throw std::runtime_error("genotype-frag: write failed for " + pbf);

        tf.flush();
        ff.flush();
        if (!tf) throw std::runtime_error("genotype-frag: write failed for " + tsv);
        if (!ff) throw std::runtime_error("genotype-frag: write failed for " + fa);
        log.info("dumped " + std::to_string(n_rows) + " scored sequences (" +
                 std::to_string(panel_graph.paths.size()) + " panel, " +
                 std::to_string(held_out.size()) + " held out); " +
                 std::to_string(n_mismatch) + " do not round-trip against the GFA path, " +
                 std::to_string(n_incomplete) + " have no complete GFA spelling");
        log.info("path blocks: " + std::to_string(n_slices) + " slices, " +
                 std::to_string(n_partial) + " partial frames, " +
                 std::to_string(n_unproj) + " unprojectable");
        log.wrote({tsv, fa, pbt, pbf});
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
                                                    hopt.placement_topk, scope_tol).panel_domain_scope;
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

    // ---- EXPOSURE PRECONDITION PROBE ------------------------------------------------------------
    // Linkage is a conditional phase score whose exposure cancels only while every constructed
    // window exceeds the insert support; below that max(0, n - L + 1) clips and the affine form
    // stops describing the exposure at all. This prints both forms so the precondition is a
    // measurement, not an assumption -- a deletion or bypass allele can easily produce a short
    // window, and then the cancellation the linkage factor relies on is simply false.
    if (!exposure_probe.empty()) {
        long min_len_e = 1;
        const std::vector<Fragment> efr = load_fragments(read_paths);
        for (const Fragment& F : efr) {
            min_len_e = std::max<long>(min_len_e,
                                       fragment_insert_floor(F, opt.allow_overlapping_pairs));
        }
        const InsertPrior ip_e = make_insert_prior(opt.fragment_len, opt.fragment_sd,
                                                   opt.discordant_rate, opt.insert_sigmas, min_len_e);
        std::printf("window\texact\taffine\tdiff\tin_regime\tinsert_lo\tinsert_hi\n");
        for (const std::string& w : split_commas(exposure_probe)) {
            const std::size_t n2 = static_cast<std::size_t>(std::stoul(w));
            const ExposureCheck c = check_exposure(n2, ip_e);
            std::printf("%zu\t%.10g\t%.10g\t%.10g\t%d\t%ld\t%ld\n",
                        n2, c.exact, c.affine, c.exact - c.affine, c.in_regime ? 1 : 0,
                        ip_e.lo, ip_e.hi);
        }
        log.done();
        return 0;
    }

    // ---- EVIDENCE OWNERSHIP ---------------------------------------------------------------------
    // Which factor each fragment belongs to, decided ONCE over the whole candidate set before any
    // genotype is scored. This is the rule that keeps marker unaries and fragment linkage factors
    // from counting the same read twice; see FragmentOwner in genotype_fragments.hpp.
    if (!ownership_table.empty()) {
        const std::vector<Fragment> ofr = load_fragments(read_paths);
        if (ofr.empty()) throw std::runtime_error("genotype-frag: --ownership-table needs reads");
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
        // THE SAME REFUSAL --reconcile-scope makes, for the same reason: a scope computed over a
        // REDUCED candidate set is a scope for a smaller, easier problem than the one a caller
        // scores, and candidate independence is the property that makes ownership meaningful.
        if (frames.size() != panel_graph.paths.size()) {
            throw std::runtime_error(
                "genotype-frag: --ownership-table has a verified walk-to-block map for only " +
                std::to_string(frames.size()) + " of " +
                std::to_string(panel_graph.paths.size()) + " panel candidates. Assigning ownership "
                "over the reduced set would make the factor topology depend on which candidates "
                "happened to decompose. Exclude them explicitly with --exclude-haplotypes, or fix "
                "the decomposition. Run --origin-universe to list them.");
        }
        long min_len_o = 1;
        for (const Fragment& F : ofr) {
            min_len_o = std::max<long>(min_len_o,
                                       fragment_insert_floor(F, opt.allow_overlapping_pairs));
        }
        const InsertPrior ip_o = make_insert_prior(opt.fragment_len, opt.fragment_sd,
                                                   opt.discordant_rate, opt.insert_sigmas, min_len_o);
        const double lep_o = std::log(opt.error_rate / 3.0);
        const double l1m_o = std::log1p(-opt.error_rate);
        std::vector<PieceIndex> opidx;
        {
            std::size_t piece = 0;
            for (const Fragment& F : ofr) {
                if (F.r1.empty()) continue;
                piece = F.r1.size() / (mate_band_edits(opt.max_divergence, F.r1.size()) + 1);
                break;
            }
            if (piece >= 12) {
                opidx.resize(frames.size());
                for (std::size_t h = 0; h < frames.size(); ++h) {
                    opidx[h] = build_piece_index(frames[h].seq, piece);
                }
            }
        }
        // WHICH BLOCKS ARE VARIABLES. A block whose sequence is identical across every candidate
        // holds no genotype state: it is context, and spanning it must not cost a fragment any
        // factor arity. n_alleles counts DISTINCT SEQUENCES and a bypassing haplotype ALREADY has
        // an allele of its own in that count, so the test is exactly n_alleles > 1. Adding
        // `|| bypass_allele >= 0` marked every degenerate flank variable, because an empty block
        // reports allele 0 as its bypass -- one state, not two.
        std::vector<char> block_variable(blocks.size(), 0);
        std::size_t n_var = 0;
        for (std::size_t b = 0; b < blocks.size(); ++b) {
            const bool var = blocks[b].n_alleles > 1;
            block_variable[b] = var ? 1 : 0;
            if (var) ++n_var;
        }
        {
            const std::string bp = ownership_table + ".blocks.tsv";
            std::ofstream bf(bp);
            if (!bf) throw std::runtime_error("genotype-frag: cannot write " + bp);
            bf << "block\tkind\tn_alleles\tbypass\tvariable\n";
            for (std::size_t b = 0; b < blocks.size() && b < chain.size(); ++b) {
                bf << b << '\t' << (chain[b].kind == BlockKind::Bubble ? "bubble" :
                                    chain[b].kind == BlockKind::Backbone ? "backbone" : "flank")
                   << '\t' << blocks[b].n_alleles << '\t' << blocks[b].bypass_allele
                   << '\t' << static_cast<int>(block_variable[b]) << '\n';
            }
            bf.flush();
            log.wrote({bp});
        }
        std::vector<FragmentOwner> owners(ofr.size());
        run_parallel(ofr.size(), opt.threads, [&](std::size_t fi) {
            owners[fi] = assign_fragment_owner_panel_domain(ofr[fi], frames, block_variable, ip_o,
                                               opt.max_divergence, lep_o, l1m_o, scope_tol,
                                               opidx.empty() ? nullptr : &opidx);
        });
        // ---- LINKAGE POTENTIALS for the edge-owned fragments -------------------------------
        // Emitted before the chain exists so the two normalisation properties can be gated on real
        // data rather than trusted: the background must stay inside the mixture, and exposure must
        // be exact rather than assumed to cancel.
        if (!linkage_potential_out.empty()) {
            // GEOMETRY COMES FROM THE SHARED BUILDER. Window, context and flank construction lives
            // in build_linkage_geometry() and nowhere else -- duplicating it here would force the
            // genotype command to re-derive the same rule, which is the pattern that produced the
            // earlier block-coordinate defects.
            std::vector<std::vector<std::string>> ballele(blocks.size());
            for (std::size_t b = 0; b < blocks.size(); ++b) ballele[b] = blocks[b].allele_seq;
            const std::size_t FLANK = static_cast<std::size_t>(ip_o.hi);
            const double lam_o = hopt.haploid_depth > 0.0 ? hopt.haploid_depth : 0.05;
            const double mix = std::log1p(-opt.outlier_mix), bgw = std::log(opt.outlier_mix);

            // Group edge-owned fragments by the edge they own.
            std::map<std::pair<std::uint32_t, std::uint32_t>, std::vector<std::size_t>> by_edge;
            for (std::size_t fi = 0; fi < ofr.size(); ++fi) {
                if (owners[fi].kind != OwnerKind::Linkage) continue;
                by_edge[{owners[fi].block_lo, owners[fi].block_hi}].push_back(fi);
            }
            std::ofstream lp(linkage_potential_out);
            if (!lp) throw std::runtime_error("genotype-frag: cannot write " + linkage_potential_out);
            lp.precision(10);
            lp << "fragment\tblock_a\tblock_b\tn_a\tn_b\tmin_window\texposure_affine"
                  "\tcis_with_bg\tcis_no_bg\tbg_effect\tinformative\n";
            const std::string ep = linkage_potential_out + ".edges.tsv";
            std::ofstream ef(ep);
            if (!ef) throw std::runtime_error("genotype-frag: cannot write " + ep);
            ef.precision(10);
            ef << "block_a\tblock_b\tn_a\tn_b\tfragments\tinformative\tinvalid\tconfigs"
                  "\texposure_affine\tpsi_mean_dev_in_class\tbest_phase_config\tmax_phase_spread"
                  "\tswap_asymmetry\tusable\tstatus\texposure_asym\n";
            std::size_t emitted = 0, refused = 0;
            for (const auto& kv : by_edge) {
                std::vector<char> gvar(blocks.size(), 0);
                for (std::size_t q = 0; q < blocks.size(); ++q) {
                    gvar[q] = blocks[q].n_alleles > 1 ? 1 : 0;
                }
                const LinkageGeometry geom = build_linkage_geometry(
                    frames, ballele, gvar, kv.first.first, kv.first.second, FLANK, ip_o);
                if (!geom.ok) {
                    refused += kv.second.size();
                    log.info("edge " + std::to_string(kv.first.first) + "-" +
                             std::to_string(kv.first.second) + " refused: " + geom.refusal);
                    continue;
                }
                std::vector<LinkageEmission> ems;
                ems.reserve(kv.second.size());
                for (std::size_t fi : kv.second) {
                    const std::size_t len = ofr[fi].bases();
                    const std::size_t be = static_cast<std::size_t>(opt.bg_divergence *
                                                                    static_cast<double>(len));
                    const double bgf = static_cast<double>(be) * lep_o +
                                       static_cast<double>(len - be) * l1m_o;
                    ems.push_back(linkage_emission(ofr[fi], geom, ip_o, opt.max_divergence,
                                                   lep_o, l1m_o, bgf));
                }
                const LinkageEdge E = aggregate_linkage_edge(ems, geom, lam_o, mix, bgw);
                // Per-fragment contrast, with and without the background, so its role stays visible.
                const auto ladd = [](double x, double y) {
                    const double ninf = -std::numeric_limits<double>::infinity();
                    if (x == ninf) return y;
                    if (y == ninf) return x;
                    const double m = std::max(x, y);
                    return m + std::log1p(std::exp(-std::abs(x - y)));
                };
                for (std::size_t k = 0; k < ems.size(); ++k) {
                    const LinkageEmission& m = ems[k];
                    if (!m.ok || m.n_a < 2 || m.n_b < 2) { ++refused; continue; }
                    const auto conf = [&](std::size_t x1, std::size_t y1, std::size_t x2,
                                          std::size_t y2, bool with_bg) {
                        const double ninf = -std::numeric_limits<double>::infinity();
                        const double p = m.mass[x1 * m.n_b + y1], q = m.mass[x2 * m.n_b + y2];
                        double sig = ninf;
                        if (p != ninf) sig = p;
                        if (q != ninf) sig = (sig == ninf) ? q : ladd(sig, q);
                        if (sig != ninf) sig += mix + std::log(lam_o);
                        if (!with_bg) return sig;
                        return (sig == ninf) ? (bgw + m.log_p_bg) : ladd(sig, bgw + m.log_p_bg);
                    };
                    const double cb = conf(0,0,1,1,true)  - conf(0,1,1,0,true);
                    const double cn = conf(0,0,1,1,false) - conf(0,1,1,0,false);
                    std::size_t wmin = std::numeric_limits<std::size_t>::max();
                    for (std::size_t w : geom.window_len) wmin = std::min(wmin, w);
                    lp << ofr[kv.second[k]].name << '\t' << geom.block_a << '\t' << geom.block_b
                       << '\t' << m.n_a << '\t' << m.n_b << '\t' << wmin << '\t'
                       << (geom.exposure_affine ? 1 : 0) << '\t' << cb << '\t' << cn << '\t'
                       << (cn - cb) << '\t' << (m.informative ? 1 : 0) << '\n';
                    ++emitted;
                }
                // EDGE LEVEL: psi must be a conditional distribution over phase within each
                // unordered-content class, so each class's psi sums to 1. Reported, not assumed.
                double worst = 0.0, swap_asym = 0.0;
                std::size_t bestc = 0; double bestv = -1e300;   // SPREAD within a class, not level
                if (E.usable()) {
                    const std::size_t na = E.n_a, nb = E.n_b;
                    std::vector<char> seen(E.log_psi.size(), 0);
                    for (std::size_t a1 = 0; a1 < na; ++a1)
                    for (std::size_t b1 = 0; b1 < nb; ++b1)
                    for (std::size_t a2 = 0; a2 < na; ++a2)
                    for (std::size_t b2 = 0; b2 < nb; ++b2) {
                        const std::size_t c = ((a1*nb + b1)*na + a2)*nb + b2;
                        if (seen[c]) continue;
                        const std::size_t A2[2] = {a1, a2}, B2[2] = {b1, b2};
                        std::vector<std::size_t> cls;
                        for (int q1 = 0; q1 < 2; ++q1) for (int q2 = 0; q2 < 2; ++q2) {
                            cls.push_back(((A2[q1]*nb + B2[q2])*na + A2[1-q1])*nb + B2[1-q2]);
                        }
                        std::sort(cls.begin(), cls.end());
                        cls.erase(std::unique(cls.begin(), cls.end()), cls.end());
                        // MEAN, not sum: psi is a mean-one likelihood ratio within the class.
                        // Measuring the sum here reported a deviation of |C|-1 on every class and
                        // was simply the stale form of this check.
                        double acc2 = 0.0;
                        for (std::size_t z : cls) { acc2 += std::exp(E.log_psi[z]); seen[z] = 1; }
                        worst = std::max(worst, std::abs(acc2 / static_cast<double>(cls.size()) - 1.0));
                        // THE PHASE SPREAD WITHIN A CLASS is the only thing psi can say. A class
                        // with one configuration has none by construction -- a homozygous endpoint
                        // has no phase to choose -- so reporting the global best log_psi would
                        // report 0 from a degenerate class and look like a decided edge.
                        if (cls.size() > 1) {
                            double hi = -1e300, lo2 = 1e300;
                            for (std::size_t z : cls) {
                                hi = std::max(hi, E.log_psi[z]);
                                lo2 = std::min(lo2, E.log_psi[z]);
                            }
                            if (hi - lo2 > bestv) {
                                bestv = hi - lo2;
                                for (std::size_t z : cls) if (E.log_psi[z] == hi) bestc = z;
                            }
                        }
                    }
                    if (bestv < -1e299) bestv = 0.0;
                    // GLOBAL HOMOLOGUE SWAP must leave psi unchanged: (a1,b1,a2,b2) and
                    // (a2,b2,a1,b1) are the same diploid state written two ways. Any asymmetry here
                    // means the ordered state leaked into the potential, and the caller would then
                    // prefer one labelling of the same genotype.
                    for (std::size_t a1 = 0; a1 < na; ++a1)
                    for (std::size_t b1 = 0; b1 < nb; ++b1)
                    for (std::size_t a2 = 0; a2 < na; ++a2)
                    for (std::size_t b2 = 0; b2 < nb; ++b2) {
                        const std::size_t c  = ((a1*nb + b1)*na + a2)*nb + b2;
                        const std::size_t cs = ((a2*nb + b2)*na + a1)*nb + b1;
                        swap_asym = std::max(swap_asym, std::abs(E.log_psi[c] - E.log_psi[cs]));
                    }
                }
                ef << geom.block_a << '\t' << geom.block_b << '\t' << E.n_a << '\t' << E.n_b
                   << '\t' << E.n_fragments << '\t' << E.n_informative << '\t' << E.n_invalid
                   << '\t' << E.log_psi.size() << '\t' << (geom.exposure_affine ? 1 : 0) << '\t'
                   << worst << '\t' << bestc << '\t' << bestv << '\t'
                   << swap_asym << '\t' << (E.usable() ? 1 : 0) << '\t' << linkage_status_name(E.status) << '\t'
                   << E.exposure_asymmetry << '\n';
                // An unusable edge must PROPAGATE, not be quietly omitted from the summary.
                if (!E.usable()) {
                    log.info("edge " + std::to_string(geom.block_a) + "-" +
                             std::to_string(geom.block_b) + " is UNSUPPORTED/INCOMPLETE: " +
                             std::string(linkage_status_name(E.status)) + " (its log psi is zeroed and must not be consumed)");
                }
            }
            lp.flush(); ef.flush();
            log.info("linkage: " + std::to_string(by_edge.size()) + " edge(s), " +
                     std::to_string(emitted) + " fragment emissions, " +
                     std::to_string(refused) + " refused");
            log.wrote({linkage_potential_out, ep});
        }

        const OwnershipLedger led = ownership_ledger(owners);
        write_ownership_table(ownership_table, ofr, owners);
        {
            const std::string lp = ownership_table + ".ledger.tsv";
            std::ofstream lf(lp);
            if (!lf) throw std::runtime_error("genotype-frag: cannot write " + lp);
            lf.precision(10);
            lf << "metric\tvalue\n"
               << "candidates\t" << frames.size() << '\n'
               << "blocks\t" << blocks.size() << '\n'
               << "blocks_variable\t" << n_var << '\n'
               << "fragments\t" << led.total << '\n'
               << "unary\t" << led.unary << '\n'
               << "linkage\t" << led.linkage << '\n'
               << "wide\t" << led.wide << '\n'
               << "invariant\t" << led.invariant << '\n'
               << "unusable\t" << led.unusable << '\n'
               // SIZE STATISTICS, not information. Pooled in-band placement mass per class; a
               // low-mass class can still hold decisive likelihood ratios, so none of these is a
               // measure of what the linkage exclusion costs. Only the block-content regression
               // can establish that. Every class is listed so none goes unaccounted.
               << "linkage_fragment_share\t" << led.linkage_fragment_share << '\n'
               << "unary_in_band_mass_share\t" << led.unary_in_band_mass_share << '\n'
               << "linkage_in_band_mass_share\t" << led.linkage_in_band_mass_share << '\n'
               << "wide_in_band_mass_share\t" << led.wide_in_band_mass_share << '\n'
               << "invariant_in_band_mass_share\t" << led.invariant_in_band_mass_share << '\n'
               << "unusable_in_band_mass_share\t" << led.unusable_in_band_mass_share << '\n';
            lf.flush();
            log.wrote({lp});
        }
        log.info("ownership over " + std::to_string(frames.size()) + " candidates, " +
                 std::to_string(n_var) + " of " + std::to_string(blocks.size()) +
                 " blocks variable: " +
                 std::to_string(led.unary) + " unary, " + std::to_string(led.linkage) +
                 " linkage, " + std::to_string(led.wide) + " wide, " +
                 std::to_string(led.invariant) + " invariant, " +
                 std::to_string(led.unusable) + " unusable");
        log.info("linkage-owned (content evidence leaves the unaries): " +
                 std::to_string(led.linkage) + " of " + std::to_string(led.total) +
                 " fragments (" + std::to_string(100.0 * led.linkage_fragment_share) +
                 "%), holding " + std::to_string(100.0 * led.linkage_in_band_mass_share) +
                 "% of pooled in-band mass -- a SIZE statistic, not the information cost");
        log.wrote({ownership_table});
        log.done();
        return 0;
    }

    // ---- INTERVAL SCORING over named candidates -------------------------------------------------
    // The real experiment path. No exhaustive A/B here: that is O(|hap|) per cell and exists to
    // CERTIFY the bounded search on fixtures, not to run on 23953 fragments. Correctness of the
    // bounded search is established by tests/genotype_bounded_search.sh; this consumes it.
    if (!interval_score.empty()) {
        if (interval_cands.empty()) {
            throw std::runtime_error("genotype-frag: --interval-score needs --interval-candidates");
        }
        const std::vector<std::string> want = split_commas(interval_cands);
        const std::vector<Fragment> ifr = load_fragments(read_paths);
        if (ifr.empty()) throw std::runtime_error("genotype-frag: --interval-score needs reads");
        const auto by_name_i = path_records_by_name(graph);
        std::vector<std::string> inames, iseqs;
        for (const std::string& nm : want) {
            const auto it = by_name_i.find(nm);
            if (it == by_name_i.end() || it->second == nullptr) {
                throw std::runtime_error("genotype-frag: --interval-candidates names a path not in "
                                         "the graph: " + nm);
            }
            bool ok = false;
            std::string w = spell_path_steps_sequence(graph, it->second->steps, &ok);
            if (!ok) throw std::runtime_error("genotype-frag: cannot spell " + nm);
            inames.push_back(nm); iseqs.push_back(std::move(w));
        }
        long min_len_i = 1;
        for (const Fragment& F : ifr) {
            min_len_i = std::max<long>(min_len_i,
                                       fragment_insert_floor(F, opt.allow_overlapping_pairs));
        }
        const InsertPrior ip_i = make_insert_prior(opt.fragment_len, opt.fragment_sd,
                                                   opt.discordant_rate, opt.insert_sigmas, min_len_i);
        const double lep_i = std::log(opt.error_rate / 3.0);
        const double l1m_i = std::log1p(-opt.error_rate);
        const double lambda_i = hopt.haploid_depth > 0.0 ? hopt.haploid_depth : 0.05;
        const double log_mix_i = std::log1p(-opt.outlier_mix);
        const double log_bgw_i = std::log(opt.outlier_mix);
        const double log_lam_i = std::log(lambda_i);
        const std::size_t nc = iseqs.size();

        // Per (fragment, candidate) mass intervals, from the adaptive tail.
        std::vector<std::vector<MassInterval>> M(ifr.size(), std::vector<MassInterval>(nc));
        // THREE DISTINCT PROPERTIES, emitted side by side so one implementation reproduces the
        // original audit partition AND identifies tail-only placements, instead of two tools
        // classifying subtly different things and being treated as interchangeable:
        //   production_band_exists -- a state exists at the FIXED production band d;
        //   adaptive_tail_exists   -- a state appears only after deeper refinement (D > d);
        //   the contribution interval -- the quantity actually used for scoring.
        // Measured: two c4 fragments differ between the two, placing on the truth-equivalent
        // haplotype ~100 nats down, reachable only at D = 4d. They contribute 0.00 to the swing.
        std::vector<std::vector<char>> PB(ifr.size(), std::vector<char>(nc, 0));
        // LAZY REFINEMENT LEVEL per cell. 1 = the production band alone -- the CHEAP COARSE BOUND,
        // exact mass within d plus the analytic bound on everything beyond it. Deepening raises the
        // lower bound and lowers the upper bound MONOTONICALLY, so a class eliminated at one level
        // can never re-enter at a deeper one. 0 = not yet computed.
        std::vector<std::vector<std::uint8_t>> LV(ifr.size(), std::vector<std::uint8_t>(nc, 0));
        // The depth cap is the same one the eager scorer uses; lazy must not win by refining less
        // deeply than the oracle was allowed to.
        const std::size_t max_mult = 4;
        // ONE INDEX PER CANDIDATE, built once and reused across all fragments. The piece length is
        // read.size()/(d+1) and reads here are uniform, so a single length covers the run; a
        // fragment whose length differs falls back to the scanning path, which returns the same
        // placements.
        std::size_t idx_piece = 0;
        for (const Fragment& F : ifr) {
            if (F.r1.empty()) continue;
            idx_piece = F.r1.size() / (mate_band_edits(opt.max_divergence, F.r1.size()) + 1);
            break;
        }
        std::vector<PieceIndex> pidx(nc);
        if (idx_piece >= 12) {
            for (std::size_t h = 0; h < nc; ++h) pidx[h] = build_piece_index(iseqs[h], idx_piece);
            log.info("piece index: " + std::to_string(nc) + " candidates at " +
                     std::to_string(idx_piece) + " bp");
        }
        std::atomic<std::size_t> done_frags{0};

        // ---- CHECKPOINT / RESUME -----------------------------------------------------------
        // A multi-hour run that loses everything on interruption is not usable. Results are
        // checkpointed per fragment BATCH and written atomically (temp + rename), so a kill at any
        // moment leaves either the previous complete checkpoint or the new one -- never a torn file.
        //
        // RESUME IS REFUSED unless the binary, reads, graph, candidate manifest and parameters all
        // match. Resuming across any of those would silently splice two different computations
        // together, which is worse than starting over.
        std::string ck_sig;
        {
            std::ostringstream sig;
            sig << "v2|" << (interval_lazy ? "lazy|" : "eager|")
                << md5_hex(std::string(reinterpret_cast<const char*>(&opt.max_divergence),
                                       sizeof(double)))
                << '|' << opt.fragment_len << '|' << opt.fragment_sd << '|' << opt.error_rate
                << '|' << opt.bg_divergence << '|' << opt.outlier_mix << '|' << lambda_i
                << '|' << interval_tol << '|' << ifr.size() << '|' << nc << '|';
            std::string names_cat;
            for (const std::string& n : inames) names_cat += n + ",";
            sig << md5_hex(names_cat);
            std::string seq_cat;
            for (const std::string& q : iseqs) seq_cat += md5_hex(q);
            sig << '|' << md5_hex(seq_cat);
            ck_sig = sig.str();
        }
        std::vector<char> frag_done(ifr.size(), 0);
        std::size_t resumed = 0;
        if (!interval_ckpt.empty()) {
            std::ifstream ck(interval_ckpt);
            if (ck) {
                std::string line;
                bool sig_ok = false;
                if (std::getline(ck, line) && line.rfind("#sig	", 0) == 0) {
                    sig_ok = (line.substr(5) == ck_sig);
                }
                if (!sig_ok) {
                    log.info("checkpoint present but its signature does not match this run "
                             "(binary, reads, graph, candidates or parameters differ); starting over");
                } else {
                    while (std::getline(ck, line)) {
                        if (line.empty() || line[0] == '#') continue;
                        std::istringstream ls(line);
                        std::size_t fi = 0;
                        if (!(ls >> fi) || fi >= ifr.size()) continue;
                        bool ok = true;
                        for (std::size_t h = 0; h < nc; ++h) {
                            double lo = 0, up = 0; int pb = 0, lv = 0;
                            if (!(ls >> lo >> up >> pb >> lv)) { ok = false; break; }
                            M[fi][h] = MassInterval{lo, up};
                            PB[fi][h] = static_cast<char>(pb);
                            LV[fi][h] = static_cast<std::uint8_t>(lv);
                        }
                        if (ok) { frag_done[fi] = 1; ++resumed; }
                    }
                    log.info("resumed " + std::to_string(resumed) + " of " +
                             std::to_string(ifr.size()) + " fragments from the checkpoint");
                }
            }
        }
        std::mutex ck_mu;
        std::vector<std::size_t> pending;
        auto ck_flush = [&]() {
            if (interval_ckpt.empty()) return;
            const std::string tmp = interval_ckpt + ".tmp";
            std::ofstream out(tmp, std::ios::trunc);
            if (!out) return;
            out.precision(17);
            out << "#sig\t" << ck_sig << '\n';
            for (std::size_t fi = 0; fi < ifr.size(); ++fi) {
                if (!frag_done[fi]) continue;
                out << fi;
                for (std::size_t h = 0; h < nc; ++h) {
                    out << ' ' << M[fi][h].lower << ' ' << M[fi][h].upper
                        << ' ' << static_cast<int>(PB[fi][h])
                        << ' ' << static_cast<int>(LV[fi][h]);
                }
                out << '\n';
            }
            out.flush();
            if (out) { out.close(); std::rename(tmp.c_str(), interval_ckpt.c_str()); }
        };
        auto hb_last = std::chrono::steady_clock::now();
        std::vector<std::size_t> depth_hist(8, 0);
        std::size_t refined_cells = 0, tol_ok_cells = 0;
        const auto t0 = std::chrono::steady_clock::now();
        run_parallel(ifr.size(), opt.threads, [&](std::size_t fi) {
            const Fragment& F = ifr[fi];
            if (F.r1.empty() || F.r2.empty()) { frag_done[fi] = 1; return; }
            if (frag_done[fi]) return;                       // already in the checkpoint
            const std::size_t d1 = mate_band_edits(opt.max_divergence, F.r1.size());
            const std::size_t d2 = mate_band_edits(opt.max_divergence, F.r2.size());
            const std::string a1 = reverse_complement(F.r1), a2 = reverse_complement(F.r2);
            for (std::size_t h = 0; h < nc; ++h) {
                // At the PRODUCTION band only -- no deepening. This is what the original audit
                // classified on.
                // The adaptive tail's FIRST iteration uses D = d1*1 = d1, which is exactly the
                // production band. Computing PB separately did that work twice per cell -- a
                // straight 2x on the dominant loop. adaptive_tail_interval now reports whether its
                // depth-1 state set was non-empty, so the production-band answer comes free.
                const PieceIndex* ix = pidx.empty() ? nullptr : &pidx[h];
                if (interval_lazy) {
                    // COARSE PASS. Every cell gets the production-band level and NOTHING more. This
                    // is a screen, not an approximation: [lower, upper] is a certified enclosure of
                    // the same quantity the eager scorer converges to, just a wider one.
                    const TailLevel lv = tail_interval_level(F.r1, F.r2, a1, a2, iseqs[h], d1, d2,
                                                             1, ip_i, lep_i, l1m_i, ix);
                    PB[fi][h] = lv.nonempty ? 1 : 0;
                    M[fi][h] = MassInterval{lv.lower, lv.upper};
                    LV[fi][h] = 1;
                } else {
                    const TailInterval ti = adaptive_tail_interval(F.r1, F.r2, iseqs[h], d1, d2, ip_i,
                                                                   lep_i, l1m_i, interval_tol, 4, ix);
                    PB[fi][h] = ti.depth1_nonempty ? 1 : 0;
                    M[fi][h] = MassInterval{ti.lower, ti.upper};
                    LV[fi][h] = static_cast<std::uint8_t>(ti.depth);
                }
            }
            // PROGRESS, so a long run is observable and can be judged rather than waited out. The
            // full-panel run emitted nothing for 12 hours and was killed with no partial result.
            frag_done[fi] = 1;
            const std::size_t n = ++done_frags;
            // TIME-BASED HEARTBEAT plus a batched atomic checkpoint. A count-based line alone can go
            // quiet for a long time when cells are slow, which is exactly when progress matters.
            {
                std::lock_guard<std::mutex> lk(ck_mu);
                pending.push_back(fi);
                const auto now = std::chrono::steady_clock::now();
                const bool due_time =
                    std::chrono::duration<double>(now - hb_last).count() >= 60.0;
                const bool due_batch = pending.size() >= interval_batch;
                if (due_time || due_batch) {
                    ck_flush();
                    pending.clear();
                    hb_last = now;
                    log.info("interval-score: " + std::to_string(n + resumed) + " / " +
                             std::to_string(ifr.size()) + " fragments" +
                             (interval_ckpt.empty() ? "" : " (checkpointed)"));
                }
            }
        });
        {
            std::lock_guard<std::mutex> lk(ck_mu);
            ck_flush();
        }
        // INCOMPLETE if any fragment is unfinished: an interrupted run must never be interpreted.
        std::size_t unfinished = 0;
        for (std::size_t fi = 0; fi < ifr.size(); ++fi) if (!frag_done[fi]) ++unfinished;

        // Background per fragment. Needed BEFORE lazy refinement, because what makes a cell worth
        // deepening is its effect on a diploid CONTRIBUTION -- which is mixed with this floor --
        // and not the width of its raw mass interval in isolation.
        std::vector<double> bg(ifr.size(), 0.0);
        for (std::size_t fi = 0; fi < ifr.size(); ++fi) {
            const std::size_t len = ifr[fi].bases();
            const std::size_t be = static_cast<std::size_t>(opt.bg_divergence *
                                                            static_cast<double>(len));
            bg[fi] = static_cast<double>(be) * lep_i +
                     static_cast<double>(len - be) * l1m_i;
        }

        // ---- CANDIDATE-LEVEL LAZY REFINEMENT --------------------------------------------------
        //
        // THE SAFETY RULE, and it is the whole design: a candidate is dropped only when EVERY
        // genotype class containing it satisfies U(C) < B - tau. A candidate can be hopeless with
        // one partner and decisive with another, so a per-candidate score cutoff -- a marker
        // shortlist, a point likelihood, any single number per haplotype -- would reopen exactly
        // the recruitment defect this whole line of work exists to characterise. Elimination is a
        // consequence of class elimination, never a judgement about a candidate on its own.
        //
        // Correctness rests on monotonicity: deepening raises lower bounds and lowers upper bounds,
        // so B never falls and an eliminated class can never re-enter. The tolerance in the pruning
        // test is the SAME tau as in the plausible-set test -- pruning at U < B would discard
        // classes that legitimately belong to the tolerance-expanded set.
        struct RoundStat {
            std::size_t round = 0, plausible = 0, candidates = 0, marked = 0, deepened = 0;
            double best_lower = 0.0, worst_other_upper = 0.0, secs = 0.0;
            std::string phase;
        };
        std::vector<RoundStat> round_log;
        std::size_t coarse_nontol = 0, lazy_rounds = 0, lazy_gap_passes = 0;
        if (interval_lazy && unfinished == 0) {
            // THE EAGER DEEPENING BASELINE, measured from this run's own coarse pass rather than
            // quoted from the oracle. The eager scorer deepens exactly those cells whose
            // production-band interval misses the tolerance, so counting them here gives the two
            // arms a common denominator instead of comparing against cells_widened, which is a
            // different quantity (a cell can end non-degenerate without ever having been deepened).
            for (std::size_t fi = 0; fi < ifr.size(); ++fi) {
                if (ifr[fi].r1.empty() || ifr[fi].r2.empty()) continue;
                for (std::size_t h = 0; h < nc; ++h) {
                    const double w = M[fi][h].upper - M[fi][h].lower;
                    if (!(w <= interval_tol)) ++coarse_nontol;      // inf/NaN-safe
                }
            }
            std::vector<std::pair<std::size_t, std::size_t>> cpair;
            for (std::size_t a = 0; a < nc; ++a)
                for (std::size_t b = a; b < nc; ++b) cpair.emplace_back(a, b);
            const std::size_t ncls = cpair.size();
            std::vector<double> exa(nc, 0.0);
            for (std::size_t h = 0; h < nc; ++h) exa[h] = ip_i.exposure(iseqs[h].size());
            std::vector<MassInterval> cls(ncls);
            auto build_classes = [&]() {
                run_parallel(ncls, opt.threads, [&](std::size_t ci) {
                    const std::size_t a = cpair[ci].first, b = cpair[ci].second;
                    const bool hom = (a == b);
                    MassInterval sum{0.0, 0.0};
                    for (std::size_t fi = 0; fi < ifr.size(); ++fi) {
                        if (ifr[fi].r1.empty() || ifr[fi].r2.empty()) continue;
                        const MassInterval c = fragment_contribution(M[fi][a], M[fi][b], hom,
                                                                     log_mix_i, log_lam_i,
                                                                     log_bgw_i, bg[fi]);
                        sum.lower += c.lower;
                        sum.upper += c.upper;
                    }
                    cls[ci] = representative_total(sum, exa[a], exa[b], hom, lambda_i, 0.0);
                });
            };
            // Deepen exactly the cells that can still move one of the named classes. "Can still
            // move" is a STRUCTURAL test, not a tuned threshold: if a cell's diploid contribution
            // interval has zero width, refining it cannot change the class interval at any depth;
            // if the cell is already exact (no omitted mass left to bound) or at the depth cap,
            // there is nothing to compute. No coefficient, no ranking, no top-k.
            auto refine_for = [&](const std::vector<std::size_t>& which,
                                  std::size_t& marked_out) -> std::size_t {
                std::vector<std::vector<char>> mk(ifr.size(), std::vector<char>(nc, 0));
                run_parallel(ifr.size(), opt.threads, [&](std::size_t fi) {
                    if (ifr[fi].r1.empty() || ifr[fi].r2.empty()) return;
                    for (std::size_t ci : which) {
                        const std::size_t a = cpair[ci].first, b = cpair[ci].second;
                        const MassInterval c = fragment_contribution(M[fi][a], M[fi][b], a == b,
                                                                     log_mix_i, log_lam_i,
                                                                     log_bgw_i, bg[fi]);
                        if (c.upper > c.lower) { mk[fi][a] = 1; mk[fi][b] = 1; }
                    }
                });
                std::atomic<std::size_t> marked{0}, deep{0};
                run_parallel(ifr.size(), opt.threads, [&](std::size_t fi) {
                    const Fragment& F = ifr[fi];
                    if (F.r1.empty() || F.r2.empty()) return;
                    const std::size_t d1 = mate_band_edits(opt.max_divergence, F.r1.size());
                    const std::size_t d2 = mate_band_edits(opt.max_divergence, F.r2.size());
                    std::string a1, a2;
                    std::size_t lm = 0, ld = 0;
                    for (std::size_t h = 0; h < nc; ++h) {
                        if (!mk[fi][h]) continue;
                        ++lm;
                        // THE SAME PER-CELL STOPPING RULE THE EAGER SCORER USES, on the same
                        // declared interval_tol -- not a second, lazy-only threshold. Without it the
                        // lazy arm could take a cell DEEPER than the oracle ever did, which would
                        // make it tighter than the thing it is being validated against and could
                        // cost more work than it saves. With it, every lazy cell depth is bounded by
                        // the eager depth for the same cell, so lazy work is a strict subset of
                        // eager work and the ledger below compares like with like.
                        if (M[fi][h].upper - M[fi][h].lower <= interval_tol) continue;
                        if (LV[fi][h] >= max_mult) continue;
                        if (M[fi][h].lower == M[fi][h].upper) continue;  // nothing omitted to bound
                        if (a1.empty()) {
                            a1 = reverse_complement(F.r1);
                            a2 = reverse_complement(F.r2);
                        }
                        const PieceIndex* ix = pidx.empty() ? nullptr : &pidx[h];
                        const std::size_t nl = static_cast<std::size_t>(LV[fi][h]) + 1;
                        const TailLevel lv = tail_interval_level(F.r1, F.r2, a1, a2, iseqs[h],
                                                                 d1, d2, nl, ip_i, lep_i, l1m_i, ix);
                        M[fi][h] = MassInterval{lv.lower, lv.upper};
                        LV[fi][h] = static_cast<std::uint8_t>(nl);
                        ++ld;
                    }
                    marked += lm;
                    deep += ld;
                });
                marked_out += marked.load();
                return deep.load();
            };
            auto elapsed = [&]() {
                return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            };
            std::size_t lead = 0;
            // SCREEN. Refine the leader (to raise B) and every class not yet safely prunable (to
            // lower its upper), until nothing contests the leader or nothing can be refined.
            for (std::size_t round = 1; ; ++round) {
                build_classes();
                lead = 0;
                for (std::size_t i = 1; i < ncls; ++i) if (cls[i].lower > cls[lead].lower) lead = i;
                const double B = cls[lead].lower;
                std::vector<std::size_t> contest;
                double wou = -std::numeric_limits<double>::infinity();
                for (std::size_t i = 0; i < ncls; ++i) {
                    if (i == lead) continue;
                    wou = std::max(wou, cls[i].upper);
                    if (!safely_prunable(cls[i], B, interval_tau)) contest.push_back(i);
                }
                std::vector<char> alive(nc, 0);
                alive[cpair[lead].first] = 1; alive[cpair[lead].second] = 1;
                for (std::size_t i : contest) {
                    alive[cpair[i].first] = 1; alive[cpair[i].second] = 1;
                }
                RoundStat rs;
                rs.round = round; rs.phase = "screen";
                rs.plausible = contest.size() + 1;
                rs.candidates = static_cast<std::size_t>(std::count(alive.begin(), alive.end(), 1));
                rs.best_lower = B; rs.worst_other_upper = wou;
                if (contest.empty()) {
                    rs.secs = elapsed();
                    round_log.push_back(rs);
                    lazy_rounds = round;
                    log.info("lazy screen: certified after " + std::to_string(round) +
                             " round(s); 1 plausible of " + std::to_string(ncls));
                    break;
                }
                // REPORTED BEFORE THE ROUND IS SPENT, not after. These counts describe the state
                // the screen is starting from, and a round costs minutes: printing them only on
                // completion hides exactly the number being measured -- how much the CHEAP bounds
                // eliminated -- behind the expensive work they were supposed to avoid.
                log.info("lazy round " + std::to_string(round) + ": " +
                         std::to_string(rs.plausible) + " plausible of " + std::to_string(ncls) +
                         ", " + std::to_string(rs.candidates) + " candidates alive, B " +
                         std::to_string(B) + ", worst rival upper " + std::to_string(wou));
                if (interval_screen_only) {
                    rs.secs = elapsed();
                    round_log.push_back(rs);
                    lazy_rounds = round;
                    log.info("screen-only: stopping before refinement");
                    break;
                }
                std::vector<std::size_t> which = contest;
                which.push_back(lead);
                std::size_t marked = 0;
                const std::size_t deep = refine_for(which, marked);
                rs.marked = marked; rs.deepened = deep; rs.secs = elapsed();
                round_log.push_back(rs);
                lazy_rounds = round;
                log.info("  round " + std::to_string(round) + " deepened " +
                         std::to_string(deep) + " of " + std::to_string(marked) +
                         " marked cells in " + std::to_string(rs.secs) + " s");
                {
                    std::lock_guard<std::mutex> lk(ck_mu);
                    ck_flush();
                }
                if (deep == 0) break;   // nothing left to refine; the verdict stands as it is
            }
            // GAP STAGE. Certification only needs every rival upper below B - tau. The REPORTED
            // robust gap is B minus the LARGEST rival upper, so reproducing the oracle's gap needs
            // that particular rival -- and the leader -- refined to the declared tolerance; every
            // other eliminated class may stay coarse. The argmax rival is recomputed each pass,
            // because lowering one upper can promote another. The leader cannot change here: a
            // rival was eliminated with upper < B - tau, so its lower is below B, and B only rises.
            for (std::size_t pass = 1; !interval_screen_only; ++pass) {
                build_classes();
                std::size_t nlead = 0;
                for (std::size_t i = 1; i < ncls; ++i) if (cls[i].lower > cls[nlead].lower) nlead = i;
                lead = nlead;
                std::size_t rival = (lead == 0 && ncls > 1) ? 1 : 0;
                for (std::size_t i = 0; i < ncls; ++i) {
                    if (i == lead) continue;
                    if (cls[i].upper > cls[rival].upper) rival = i;
                }
                if (rival == lead) break;
                const double wl = cls[lead].upper - cls[lead].lower;
                const double wr = cls[rival].upper - cls[rival].lower;
                RoundStat rs;
                rs.round = pass; rs.phase = "gap";
                rs.plausible = 1; rs.best_lower = cls[lead].lower;
                rs.worst_other_upper = cls[rival].upper;
                rs.candidates = 0;
                if (wl <= interval_tol && wr <= interval_tol) {
                    rs.secs = elapsed();
                    round_log.push_back(rs);
                    break;
                }
                std::size_t marked = 0;
                const std::size_t deep = refine_for({lead, rival}, marked);
                rs.marked = marked; rs.deepened = deep; rs.secs = elapsed();
                round_log.push_back(rs);
                lazy_gap_passes = pass;
                log.info("lazy gap pass " + std::to_string(pass) + ": leader width " +
                         std::to_string(wl) + ", rival width " + std::to_string(wr) +
                         ", deepened " + std::to_string(deep) + " cells");
                {
                    std::lock_guard<std::mutex> lk(ck_mu);
                    ck_flush();
                }
                if (deep == 0) break;
            }
        }

        // Cell-level accounting, over the FINAL state of every cell.
        std::size_t deepened_cells = 0, level_steps = 0, coarse_cells = 0;
        for (std::size_t fi = 0; fi < ifr.size(); ++fi) {
            for (std::size_t h = 0; h < nc; ++h) {
                if (M[fi][h].lower != M[fi][h].upper) ++refined_cells; else ++tol_ok_cells;
                if (ifr[fi].r1.empty() || ifr[fi].r2.empty()) continue;
                ++coarse_cells;
                if (LV[fi][h] > 1) {
                    ++deepened_cells;
                    level_steps += static_cast<std::size_t>(LV[fi][h]) - 1;
                }
            }
        }
        // PER-FRAGMENT, PER-CANDIDATE placement table. This is what a classification consumer
        // needs, and it costs nothing extra here -- the intervals are already computed. The audit
        // previously called --bounded-search for this, which runs the exhaustive A/B on every cell:
        // that is the tool that CERTIFIES the bounded search on fixtures, and on three 226 kb
        // haplotypes it did not finish in 10 minutes. Verification and use are different jobs.
        {
            const std::string fp = interval_score + ".placements.tsv";
            std::ofstream pf(fp);
            if (!pf) throw std::runtime_error("genotype-frag: cannot write " + fp);
            pf << "fragment\tcandidate\tproduction_band_exists\tadaptive_tail_exists"
                  "\ttail_only\tlower\tupper\n";
            for (std::size_t fi = 0; fi < ifr.size(); ++fi) {
                if (ifr[fi].r1.empty() || ifr[fi].r2.empty()) continue;
                for (std::size_t h = 0; h < nc; ++h) {
                    const bool ad =
                        M[fi][h].lower != -std::numeric_limits<double>::infinity();
                    const bool pb = PB[fi][h] != 0;
                    pf << ifr[fi].name << '\t' << inames[h] << '\t' << (pb ? 1 : 0) << '\t'
                       << (ad ? 1 : 0) << '\t' << ((ad && !pb) ? 1 : 0) << '\t'
                       << M[fi][h].lower << '\t' << M[fi][h].upper << '\n';
                }
            }
            pf.flush();
            log.wrote({fp});
        }
        // Six unordered diplotypes for three candidates, INCLUDING homozygotes.
        std::ofstream isf(interval_score);
        if (!isf) throw std::runtime_error("genotype-frag: cannot write " + interval_score);
        isf.precision(10);
        isf << "class\thap_a\thap_b\thomozygous\tlower\tupper\twidth\tnominal\n";
        // PER-FRAGMENT CONTRIBUTIONS per class, so a score swing can be ATTRIBUTED rather than
        // assumed. 12697 fragment-candidate cells widened under complete recruitment, including
        // fragments outside the audited 76 -- "the 20 reclassified fragments caused the reversal"
        // does not follow from the reversal alone and has to be decomposed.
        // OFF BY DEFAULT. One row per (fragment, class): at 131 candidates that is 8646 unordered
        // diplotypes x 23953 fragments = 207 million rows. It exists for ATTRIBUTION on a handful of
        // named candidates, not for a full shortlist.
        const std::string cfp = interval_score + ".contrib.tsv";
        std::ofstream cf;
        if (interval_contrib) cf.open(cfp);
        if (interval_contrib && !cf) throw std::runtime_error("genotype-frag: cannot write " + cfp);
        if (interval_contrib) { cf.precision(12); cf << "fragment\tclass\tcontrib_lower\n"; }
        std::vector<MassInterval> classes;
        std::vector<std::string> clabel;
        std::vector<double> nominal;
        for (std::size_t a = 0; a < nc; ++a) {
            for (std::size_t b = a; b < nc; ++b) {
                const bool hom = (a == b);
                MassInterval sum{0.0, 0.0};
                for (std::size_t fi = 0; fi < ifr.size(); ++fi) {
                    if (ifr[fi].r1.empty() || ifr[fi].r2.empty()) continue;
                    const MassInterval c = fragment_contribution(M[fi][a], M[fi][b], hom,
                                                                 log_mix_i, log_lam_i,
                                                                 log_bgw_i, bg[fi]);
                    sum.lower += c.lower;
                    sum.upper += c.upper;
                    if (interval_contrib) {
                        cf << ifr[fi].name << '\t' << a << '_' << b << '\t' << c.lower << '\n';
                    }
                }
                const double ea = ip_i.exposure(iseqs[a].size());
                const double eb = ip_i.exposure(iseqs[b].size());
                const MassInterval tot = representative_total(sum, ea, eb, hom, lambda_i, 0.0);
                classes.push_back(tot);
                clabel.push_back(inames[a] + " | " + inames[b]);
                nominal.push_back(tot.lower);
                isf << clabel.back() << '\t' << inames[a] << '\t' << inames[b] << '\t'
                    << (hom ? 1 : 0) << '\t' << tot.lower << '\t' << tot.upper << '\t'
                    << (tot.upper - tot.lower) << '\t' << tot.lower << '\n';
            }
        }
        if (interval_contrib) { cf.flush(); log.wrote({cfp}); }
        const auto t1 = std::chrono::steady_clock::now();
        const double secs = std::chrono::duration<double>(t1 - t0).count();
        const Certification cert = certify(classes, interval_tau);
        // ROBUST LEADER is argmax L(C), not the nominal point winner. Reported separately so a
        // divergence between them is visible rather than assumed away.
        std::size_t leader = 0;
        for (std::size_t i = 1; i < classes.size(); ++i) {
            if (classes[i].lower > classes[leader].lower) leader = i;
        }
        std::size_t nom_win = 0;
        for (std::size_t i = 1; i < nominal.size(); ++i) {
            if (nominal[i] > nominal[nom_win]) nom_win = i;
        }
        double worst_other_upper = -std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i < classes.size(); ++i) {
            if (i == leader) continue;
            worst_other_upper = std::max(worst_other_upper, classes[i].upper);
        }
        const double robust_gap = classes[leader].lower - worst_other_upper;
        // THREE OUTCOMES, and the third must not be reported as biological ambiguity. A run stopped
        // by the depth cap with intervals still overlapping is UNFINISHED, not unresolved.
        Verdict vd = verdict_of(classes, cert, interval_tol);
        if (unfinished > 0) vd = Verdict::Incomplete;
        const char* verdict = verdict_name(vd);
        const bool all_tight = (vd != Verdict::Incomplete);
        isf << "# robust_leader\t" << clabel[leader] << '\n';
        isf << "# nominal_winner\t" << clabel[nom_win] << '\n';
        isf << "# leaders_agree\t" << (leader == nom_win ? 1 : 0) << '\n';
        isf << "# robust_gap\t" << robust_gap << '\n';
        isf << "# tau\t" << interval_tau << '\n';
        isf << "# n_plausible\t" << cert.plausible.size() << '\n';
        isf << "# verdict\t" << verdict << '\n';
        isf << "# all_plausible_tight\t" << (all_tight ? 1 : 0) << '\n';
        isf << "# interval_tol\t" << interval_tol << '\n';
        isf << "# fragments\t" << ifr.size() << '\n';
        isf << "# unfinished\t" << unfinished << '\n';
        isf << "# resumed\t" << resumed << '\n';
        isf << "# candidates\t" << nc << '\n';
        isf << "# diplotypes\t" << classes.size() << '\n';
        isf << "# cells_widened\t" << refined_cells << '\n';
        isf << "# cells_exact\t" << tol_ok_cells << '\n';
        isf << "# seconds\t" << secs << '\n';
        // CPU time alongside wall time. On four workers they differ by up to 4x, and an
        // acceleration that only moves wall time by taking more cores is not an acceleration.
        double cpu_secs = 0.0;
        {
            struct rusage ru {};
            if (getrusage(RUSAGE_SELF, &ru) == 0) {
                cpu_secs = static_cast<double>(ru.ru_utime.tv_sec + ru.ru_stime.tv_sec) +
                           1e-6 * static_cast<double>(ru.ru_utime.tv_usec + ru.ru_stime.tv_usec);
            }
        }
        isf << "# cpu_seconds\t" << cpu_secs << '\n';
        // THE ACCELERATION LEDGER. Four different counts that have been conflated before, so each
        // is named for exactly what it measures:
        //   coarse_cells        every fragment x candidate cell, all of which get a cheap bound;
        //   coarse_nontolerant  cells whose production-band interval misses the tolerance -- the
        //                       set the EAGER scorer deepens, and therefore the baseline. Only the
        //                       lazy arm can report it directly, because only it keeps a cell at
        //                       level 1; for the eager arm the SAME set is deepened_cells, since
        //                       eager deepens a cell exactly when level 1 missed the tolerance.
        //                       So lazy's coarse_nontolerant and eager's deepened_cells are two
        //                       measurements of one quantity, and their agreement is a check;
        //   deepened_cells      cells this run actually took beyond the production band;
        //   level_steps         total depth increments, since a cell can be deepened repeatedly.
        // Comparing deepened_cells against the eager run's cells_widened would compare two
        // different quantities: a cell can end non-degenerate without ever having been deepened.
        isf << "# lazy\t" << (interval_lazy ? 1 : 0) << '\n';
        isf << "# coarse_cells\t" << coarse_cells << '\n';
        isf << "# coarse_nontolerant\t" << coarse_nontol << '\n';
        isf << "# deepened_cells\t" << deepened_cells << '\n';
        isf << "# level_steps\t" << level_steps << '\n';
        isf << "# lazy_screen_rounds\t" << lazy_rounds << '\n';
        isf << "# lazy_gap_passes\t" << lazy_gap_passes << '\n';
        isf.flush();
        if (!isf) throw std::runtime_error("genotype-frag: write failed for " + interval_score);
        log.info(std::string(verdict) + ": robust leader " + clabel[leader] + ", gap " +
                 std::to_string(robust_gap) + ", " + std::to_string(cert.plausible.size()) +
                 " plausible of " + std::to_string(classes.size()));
        if (!interval_rounds.empty()) {
            std::ofstream rf(interval_rounds);
            if (!rf) throw std::runtime_error("genotype-frag: cannot write " + interval_rounds);
            rf.precision(10);
            rf << "phase\tround\tplausible\tcandidates_alive\tcells_marked\tcells_deepened"
                  "\tbest_lower\tworst_other_upper\tseconds\n";
            for (const RoundStat& r : round_log) {
                rf << r.phase << '\t' << r.round << '\t' << r.plausible << '\t' << r.candidates
                   << '\t' << r.marked << '\t' << r.deepened << '\t' << r.best_lower << '\t'
                   << r.worst_other_upper << '\t' << r.secs << '\n';
            }
            rf.flush();
            log.wrote({interval_rounds});
        }
        log.wrote({interval_score});
        log.done();
        return 0;
    }

    // ---- BOUNDED-COMPLETE SINGLE-MATE SEARCH, against the exhaustive reference ------------------
    // Stage 1 is HAMMING only: reference_emission counts mismatches at a fixed offset and has no gap
    // model, so this is the emission it can certify. An indel-aware search needs its own oracle.
    if (!bounded_search.empty()) {
        const std::vector<Fragment> bfr = load_fragments(read_paths);
        if (bfr.empty()) throw std::runtime_error("genotype-frag: --bounded-search needs reads");
        std::vector<std::string> names;
        std::vector<std::string> seqs;
        // HOISTED. path_records_by_name returns BY VALUE, so calling it inside the loop left `it`
        // pointing into a destroyed temporary and compared it against an iterator from a second,
        // different temporary -- undefined behaviour that segfaulted immediately.
        const auto by_name_bs = path_records_by_name(graph);
        for (const PathRecord& p : panel_graph.paths) {
            const auto it = by_name_bs.find(p.name);
            if (it == by_name_bs.end() || it->second == nullptr) continue;
            bool ok = false;
            std::string w = spell_path_steps_sequence(graph, it->second->steps, &ok);
            if (!ok) continue;
            names.push_back(p.name); seqs.push_back(std::move(w));
        }
        if (seqs.empty()) throw std::runtime_error("genotype-frag: --bounded-search found no panel walks");
        std::ofstream bf(bounded_search);
        if (!bf) throw std::runtime_error("genotype-frag: cannot write " + bounded_search);
        bf << "fragment\tmate\tstrand\thaplotype\tbounded_n\texhaustive_n\tagree"
              "\tbest_edits\tmax_edits\tpieces\tcandidate_starts\tdistinct_starts\tverified\tfallback\n";
        std::size_t rows = 0, disagree = 0;
        std::uint64_t tot_cand = 0, tot_ver = 0, tot_exh = 0;
        // FRAGMENT STATES, built independently from the bounded and exhaustive placement vectors.
        // Comparing only single-mate placements cannot catch a lost library orientation, a wrong
        // combination of the four vectors, collapsed states, or lost haplotype identity.
        const std::string fs_path = bounded_search + ".states.tsv";
        std::ofstream sf(fs_path);
        if (!sf) throw std::runtime_error("genotype-frag: cannot write " + fs_path);
        sf << "fragment\thaplotype\tbounded_states\texhaustive_states\tagree\torientA\torientB"
              "\tinsert_lo\tinsert_hi\tstate_hap\tmass_bounded\tmass_exhaustive\tmass_reference"
              "\tomitted_bound\tupper\tcontains_ref\tinterval_nats\texposure"
              "\tad_depth\tad_lower\tad_upper\tad_width\tad_tol_ok\tad_contains"
              "\tad_lower_mono\tad_upper_mono\tad_inband_stable\tad_inband_at_d\n";
        std::size_t st_dis = 0, tot_a = 0, tot_b = 0, tot_states = 0;
        std::size_t mass_dis = 0, bound_viol = 0, mass_cells = 0, interval_miss = 0;
        std::size_t adaptive_miss = 0, adaptive_tight = 0;
        long bs_prior_lo = 0, bs_prior_hi = 0;
        InsertPrior bs_prior_g;
        {
            // THE PRODUCTION PRIOR, not a re-derived one. make_insert_prior uses 4 sigmas and
            // floors lo at the mates' combined length: for 120 bp mates that is [240,550], where a
            // hard-coded mean +/- 3sd gives [200,550-50] = [200,500]. Comparing states under the
            // wrong support both admits states the model rejects (insert 200-239) and omits valid
            // ones (501-550), so the audit would classify against a different fragment model than
            // the production C4 result it is meant to explain.
            // MAXIMUM, from 1, over EVERY loaded fragment -- singletons included. One shared insert
            // prior has to be valid for the LONGEST fragment, so the longest sets the lower support
            // bound; a minimum would put the bound below the length some fragments already occupy.
            // This mirrors the production loop exactly (std::max<long>, seeded at 1, no skipping).
            // The first version took the minimum AND skipped singletons, and a fixture of uniform
            // 120+120 pairs cannot tell the two apart -- min and max are both 240 there.
            long min_frag_len_bs = 1;
            for (const Fragment& F : bfr) {
                min_frag_len_bs = std::max<long>(min_frag_len_bs,
                                    fragment_insert_floor(F, opt.allow_overlapping_pairs));
            }
            const InsertPrior bs_prior = make_insert_prior(opt.fragment_len, opt.fragment_sd,
                                                           opt.discordant_rate, opt.insert_sigmas, min_frag_len_bs);
            const long ilo = bs_prior.lo, ihi = bs_prior.hi;
            bs_prior_lo = ilo; bs_prior_hi = ihi; bs_prior_g = bs_prior;
            log.info("bounded-search insert support [" + std::to_string(ilo) + "," +
                     std::to_string(ihi) + "] from the production prior");
            for (const Fragment& F : bfr) {
                if (F.r1.empty() || F.r2.empty()) continue;
                const std::string r1rc = reverse_complement(F.r1);
                const std::string r2rc = reverse_complement(F.r2);
                const std::size_t d1 = mate_band_edits(opt.max_divergence, F.r1.size());
                const std::size_t d2 = mate_band_edits(opt.max_divergence, F.r2.size());
                for (std::size_t h = 0; h < seqs.size(); ++h) {
                    const auto bf1 = bounded_mate_placements(F.r1, seqs[h], d1, nullptr);
                    const auto br1 = bounded_mate_placements(r1rc, seqs[h], d1, nullptr);
                    const auto bf2 = bounded_mate_placements(F.r2, seqs[h], d2, nullptr);
                    const auto br2 = bounded_mate_placements(r2rc, seqs[h], d2, nullptr);
                    const auto ef1 = exhaustive_mate_placements(F.r1, seqs[h], d1);
                    const auto er1 = exhaustive_mate_placements(r1rc, seqs[h], d1);
                    const auto ef2 = exhaustive_mate_placements(F.r2, seqs[h], d2);
                    const auto er2 = exhaustive_mate_placements(r2rc, seqs[h], d2);
                    const auto bs = enumerate_fragment_states(static_cast<std::uint32_t>(h),
                        bf1, br1, bf2, br2, F.r1.size(), F.r2.size(), ilo, ihi);
                    const auto es = enumerate_fragment_states(static_cast<std::uint32_t>(h),
                        ef1, er1, ef2, er2, F.r1.size(), F.r2.size(), ilo, ihi);
                    const bool same = bs.size() == es.size() &&
                        std::equal(bs.begin(), bs.end(), es.begin(),
                                   [](const FragmentState& x, const FragmentState& y) {
                                       return x == y && x.m1_edits == y.m1_edits &&
                                              x.m2_edits == y.m2_edits; });
                    if (!same) ++st_dis;
                    std::size_t na = 0, nb = 0;
                    for (const FragmentState& z : bs) { if (z.m1_fwd) ++na; else ++nb; }
                    tot_a += na; tot_b += nb; tot_states += bs.size();
                    // IN-BAND MASS, in the reference's terms, plus the untruncated reference for
                    // the same fragment and haplotype. The first two must be EQUAL (same states,
                    // same formula); the third must be >= them, because the reference integrates
                    // out-of-band states this search deliberately does not reach.
                    const double lep = std::log(opt.error_rate / 3.0);
                    const double l1m = std::log1p(-opt.error_rate);
                    const double mb = fragment_states_mass(bs, F.r1.size(), F.r2.size(),
                                                           bs_prior, lep, l1m);
                    const double me = fragment_states_mass(es, F.r1.size(), F.r2.size(),
                                                           bs_prior, lep, l1m);
                    ReferenceParams rp_bs;
                    rp_bs.error_rate = opt.error_rate; rp_bs.fragment_len = opt.fragment_len;
                    rp_bs.fragment_sd = opt.fragment_sd; rp_bs.bg_divergence = opt.bg_divergence;
                    const double mr = reference_fragment_on_haplotype(
                        F, seqs[h], rp_bs, bs_prior, reverse_complement(F.r2), lep, l1m);
                    const double kNegInfBs = -std::numeric_limits<double>::infinity();
                    const bool mass_eq = (mb == kNegInfBs && me == kNegInfBs) ||
                                         std::abs(mb - me) < 1e-9;
                    // 1e-9 absolute on log mass: these are the same additions in the same order, so
                    // only floating-point associativity can separate them.
                    if (!mass_eq) ++mass_dis;
                    if (mb != kNegInfBs && mr != kNegInfBs && mb > mr + 1e-9) ++bound_viol;
                    // THE INTERVAL. lower = in-band mass, upper = logadd(lower, omitted bound).
                    // The exact reference must lie inside it, or the bound does not bound.
                    const double ob = omitted_mass_bound(seqs[h].size(), F.r1.size(), F.r2.size(),
                                                         d1, d2, bs_prior, lep, l1m, bs);
                    const double up = (mb == kNegInfBs) ? ob
                                    : (ob == kNegInfBs ? mb
                                       : std::max(mb, ob) + std::log1p(std::exp(-std::abs(mb - ob))));
                    const bool contains = (mr == kNegInfBs) ||
                                          (mr <= up + 1e-9 &&
                                           (mb == kNegInfBs || mr >= mb - 1e-9));
                    if (!contains) ++interval_miss;
                    const double width = (up == kNegInfBs || mb == kNegInfBs) ? -1.0 : up - mb;
                    const double expo = bs_prior.exposure(seqs[h].size());
                    // ADAPTIVE TAIL: deepen D until the width meets tolerance, or report it as
                    // uncertifiable. 1 nat is a deliberately demanding target for a single
                    // fragment-haplotype cell.
                    const TailInterval ti = adaptive_tail_interval(
                        F.r1, F.r2, seqs[h], d1, d2, bs_prior, lep, l1m, 1.0, 6);
                    const double ti_width = (ti.upper == kNegInfBs || ti.lower == kNegInfBs)
                                          ? -1.0 : ti.upper - ti.lower;
                    const bool ti_ok = (mr == kNegInfBs) ||
                                       (mr <= ti.upper + 1e-9 &&
                                        (ti.lower == kNegInfBs || mr >= ti.lower - 1e-9));
                    if (!ti_ok) ++adaptive_miss;
                    if (ti.within_tolerance) ++adaptive_tight;
                    // EMIT EVEN WITH NO STATES. The tail-only case -- every placement out of band,
                    // so the in-band set is empty and lower = -inf -- is precisely where the
                    // omitted bound has to do work, and skipping empty cells hid it entirely.
                    {
                        // state_hap is FragmentState::hap, not the loop variable: printing the
                        // enclosing name would show the right answer even if the key lost it.
                        sf << F.name << '\t' << names[h] << '\t' << bs.size() << '\t' << es.size()
                           << '\t' << (same ? "yes" : "NO") << '\t' << na << '\t' << nb
                           << '\t' << ilo << '\t' << ihi << '\t'
                           << (bs.empty() ? h : bs.front().hap) << '\t'
                           << mb << '\t' << me << '\t' << mr << '\t'
                           << ob << '\t' << up << '\t' << (contains ? 1 : 0) << '\t'
                           << width << '\t' << expo << '\t'
                           << ti.depth << '\t' << ti.lower << '\t' << ti.upper << '\t'
                           << ti_width << '\t' << (ti.within_tolerance ? 1 : 0) << '\t'
                           << (ti_ok ? 1 : 0) << '\t'
                           << (ti.lower_monotone ? 1 : 0) << '\t'
                           << (ti.upper_monotone ? 1 : 0) << '\t'
                           << (ti.inband_stable ? 1 : 0) << '\t' << ti.inband_at_d << '\n';
                        ++mass_cells;
                    }
                }
            }
        // THE THREE BLOCKS BELOW ARE FIXTURE-SCALE and run only under --bounded-verify. Each is
        // O(|haplotype|) or worse PER CELL: the exhaustive A/B, the edit-class extraction at D+1,
        // and the D-boundary case at D = read length, which is a full exhaustive scan by
        // construction. Measured: leaving them unconditional made the 76-fragment c4 audit -- three
        // haplotypes of ~226 kb -- run for 90 minutes without finishing. They CERTIFY the bounded
        // search on small fixtures; they are not part of using it.
        if (bounded_verify) {
        // ---- MULTIPLICITY, isolated to one edit class -----------------------------------------
        // K distinct states at the SAME likelihood must sum to log K above a single one, and all K
        // coordinates must survive deduplication. Reported per (fragment, haplotype) at the exact
        // (d1+1, 0) class, so a repeat's omitted copies are visible as a COUNT rather than inferred
        // from a log-space subtraction.
        {
            const std::string mp = bounded_search + ".editclass.tsv";
            std::ofstream mf(mp);
            if (!mf) throw std::runtime_error("genotype-frag: cannot write " + mp);
            mf << "fragment\thaplotype\te1\te2\tcount\tclass_mass\tper_state_mass"
                  "\tmin_start\tmax_start\tspan\tdistinct_starts\n";
            const double lep2 = std::log(opt.error_rate / 3.0);
            const double l1m2 = std::log1p(-opt.error_rate);
            for (const Fragment& F : bfr) {
                if (F.r1.empty() || F.r2.empty()) continue;
                const std::size_t D1 = mate_band_edits(opt.max_divergence, F.r1.size()) + 1;
                const std::size_t D2 = mate_band_edits(opt.max_divergence, F.r2.size()) + 1;
                const std::string a1 = reverse_complement(F.r1), a2 = reverse_complement(F.r2);
                for (std::size_t h = 0; h < seqs.size(); ++h) {
                    const auto st = enumerate_fragment_states(static_cast<std::uint32_t>(h),
                        bounded_mate_placements(F.r1, seqs[h], D1, nullptr),
                        bounded_mate_placements(a1, seqs[h], D1, nullptr),
                        bounded_mate_placements(F.r2, seqs[h], D2, nullptr),
                        bounded_mate_placements(a2, seqs[h], D2, nullptr),
                        F.r1.size(), F.r2.size(), bs_prior_lo, bs_prior_hi);
                    for (std::uint32_t e1 = 0; e1 <= 1; ++e1) {
                        const EditClass ec = edit_class_mass(st, e1, 0, F.r1.size(), F.r2.size(),
                                                             bs_prior_g, lep2, l1m2);
                        if (ec.count == 0) continue;
                        // The COUNT alone is not enough: extra starts, inserts or orientations
                        // inside the first eight copies could reach 9 or 16 while later copies stay
                        // truncated. The coordinate SPAN shows which copies actually contributed.
                        long lo_s = 0, hi_s = 0; bool firstv = true;
                        std::set<long> starts;
                        for (const FragmentState& z : st) {
                            if (z.m1_edits != e1 || z.m2_edits != 0) continue;
                            starts.insert(z.frag_start);
                            if (firstv) { lo_s = hi_s = z.frag_start; firstv = false; }
                            lo_s = std::min(lo_s, z.frag_start);
                            hi_s = std::max(hi_s, z.frag_start);
                        }
                        mf << F.name << '\t' << names[h] << '\t' << e1 << "\t0\t"
                           << ec.count << '\t' << ec.mass << '\t'
                           << (ec.mass - std::log(static_cast<double>(ec.count))) << '\t'
                           << lo_s << '\t' << hi_s << '\t' << (hi_s - lo_s) << '\t'
                           << starts.size() << '\n';
                    }
                }
            }
            mf.flush();
            log.wrote({mp});
        }

        // ---- D-BOUNDARY: the terminal case ----------------------------------------------------
        // At D equal to the maximum possible mate mismatch count -- the read length -- every
        // (start, insert, orientation) in the universe is in band, so the residual set is EMPTY,
        // the bound is -inf, and lower, upper and the exhaustive mass must agree numerically. This
        // closes the induction: every intermediate D is bracketed between the production band and
        // an endpoint where the bound is provably exact.
        {
            const std::string bp = bounded_search + ".dboundary.tsv";
            std::ofstream df(bp);
            if (!df) throw std::runtime_error("genotype-frag: cannot write " + bp);
            df << "fragment\thaplotype\tn_omitted_bound\tlower\tupper\treference\tagree\n";
            const double lep3 = std::log(opt.error_rate / 3.0);
            const double l1m3 = std::log1p(-opt.error_rate);
            std::size_t dbad = 0;
            for (const Fragment& F : bfr) {
                if (F.r1.empty() || F.r2.empty()) continue;
                const std::string a1 = reverse_complement(F.r1), a2 = reverse_complement(F.r2);
                for (std::size_t h = 0; h < seqs.size(); ++h) {
                    const std::size_t D1 = F.r1.size(), D2 = F.r2.size();
                    const auto st = enumerate_fragment_states(static_cast<std::uint32_t>(h),
                        bounded_mate_placements(F.r1, seqs[h], D1, nullptr),
                        bounded_mate_placements(a1, seqs[h], D1, nullptr),
                        bounded_mate_placements(F.r2, seqs[h], D2, nullptr),
                        bounded_mate_placements(a2, seqs[h], D2, nullptr),
                        F.r1.size(), F.r2.size(), bs_prior_lo, bs_prior_hi);
                    const double m = fragment_states_mass(st, F.r1.size(), F.r2.size(),
                                                          bs_prior_g, lep3, l1m3);
                    const double b = omitted_mass_bound(seqs[h].size(), F.r1.size(), F.r2.size(),
                                                        D1, D2, bs_prior_g, lep3, l1m3, st);
                    ReferenceParams rp3;
                    rp3.error_rate = opt.error_rate; rp3.fragment_len = opt.fragment_len;
                    rp3.fragment_sd = opt.fragment_sd; rp3.bg_divergence = opt.bg_divergence;
                    const double r = reference_fragment_on_haplotype(F, seqs[h], rp3, bs_prior_g,
                                                                     a2, lep3, l1m3);
                    const bool agree = (b == -std::numeric_limits<double>::infinity()) &&
                                       ((m == r) || (m != -std::numeric_limits<double>::infinity() &&
                                                     r != -std::numeric_limits<double>::infinity() &&
                                                     std::abs(m - r) < 1e-9));
                    if (!agree) ++dbad;
                    df << F.name << '\t' << names[h] << '\t' << b << '\t' << m << '\t' << m
                       << '\t' << r << '\t' << (agree ? 1 : 0) << '\n';
                }
            }
            df.flush();
            log.wrote({bp});
            log.info("D-boundary: " + std::to_string(dbad) + " cell(s) where the residual is "
                     "non-empty or lower/upper/reference disagree at D = read length");
        }

        // ---- AGGREGATION GATES ------------------------------------------------------------------
        // Analytic fixtures, not data-driven: each isolates one property of the order
        // per-representative sum -> exposure -> class aggregation.
        {
            const std::string ap = bounded_search + ".aggregate.tsv";
            std::ofstream af(ap);
            if (!af) throw std::runtime_error("genotype-frag: cannot write " + ap);
            af << "case\tvalue_a\tvalue_b\tequal\n";
            const double kNegInfBs2 = -std::numeric_limits<double>::infinity();
            const double lmix = std::log(0.9), llam = std::log(0.05);
            const double lbgw = std::log(0.1), lpbg = -12.0;
            const auto eq = [](double x, double y) { return std::abs(x - y) < 1e-9; };
            std::size_t abad = 0;
            const auto emit = [&](const char* nm, double x, double y) {
                const bool e = eq(x, y);
                if (!e) ++abad;
                af << nm << '\t' << x << '\t' << y << '\t' << (e ? 1 : 0) << '\n';
            };
            const MassInterval A{-3.0, -2.5}, B{-4.0, -3.2};

            // 1. SWAPPING a and b leaves the contribution unchanged: the pair is unordered.
            const auto ab = fragment_contribution(A, B, false, lmix, llam, lbgw, lpbg);
            const auto ba = fragment_contribution(B, A, false, lmix, llam, lbgw, lpbg);
            emit("swap_ab_lower", ab.lower, ba.lower);
            emit("swap_ab_upper", ab.upper, ba.upper);

            // 2. HOMOZYGOUS DOUBLING is present EXACTLY ONCE: 2*M_a, and the background is mixed
            //    once. Compared against combining A with itself as if heterozygous, which gives the
            //    same 2*M_a -- the two must agree, proving the doubling is not applied twice.
            const auto hom = fragment_contribution(A, A, true, lmix, llam, lbgw, lpbg);
            const auto het_aa = fragment_contribution(A, A, false, lmix, llam, lbgw, lpbg);
            emit("hom_doubling_once", hom.lower, het_aa.lower);

            // 3. ADDING AN IDENTICAL ALIAS must not change the class: weights sum to one, so
            //    catalogue multiplicity cannot manufacture confidence.
            const MassInterval R{-5.0, -4.0};
            const auto c1 = aggregate_class({R}, {1.0});
            const auto c2 = aggregate_class({R, R}, {0.5, 0.5});
            emit("alias_lower", c1.lower, c2.lower);
            emit("alias_upper", c1.upper, c2.upper);

            // 4. SPLITTING one representative into two half-weight aliases is a no-op.
            const auto c3 = aggregate_class({R, R, R}, {0.5, 0.25, 0.25});
            emit("split_lower", c1.lower, c3.lower);

            // 5. NAIVE SUMMING (weights all 1.0) must NOT equal the weighted class -- otherwise the
            //    weighting is inert and duplicates would inflate the score. This one asserts a
            //    DIFFERENCE, so a degenerate implementation cannot pass everything above.
            const auto naive = aggregate_class({R, R}, {1.0, 1.0});
            af << "naive_differs\t" << naive.lower << '\t' << c1.lower << '\t'
               << (eq(naive.lower, c1.lower) ? 0 : 1) << '\n';
            if (eq(naive.lower, c1.lower)) ++abad;

            // 6. PER-FRAGMENT aggregation differs from per-representative aggregation whenever
            //    different representatives explain different fragments. Two fragments, two
            //    representatives, each strong on one: aggregating per fragment lets BOTH be
            //    explained well -- a mosaic -- while the correct order forces one representative to
            //    explain both. The two must NOT agree.
            const MassInterval P1{-1.0, -1.0}, P2{-9.0, -9.0};
            const auto rep1 = fragment_contribution(P1, P1, true, lmix, llam, lbgw, lpbg);
            const auto rep1b = fragment_contribution(P2, P2, true, lmix, llam, lbgw, lpbg);
            const auto rep2 = fragment_contribution(P2, P2, true, lmix, llam, lbgw, lpbg);
            const auto rep2b = fragment_contribution(P1, P1, true, lmix, llam, lbgw, lpbg);
            // correct: sum per representative, then aggregate the two totals
            const MassInterval t1{rep1.lower + rep1b.lower, rep1.upper + rep1b.upper};
            const MassInterval t2{rep2.lower + rep2b.lower, rep2.upper + rep2b.upper};
            const auto correct = aggregate_class({t1, t2}, {0.5, 0.5});
            // mosaic: aggregate per fragment, then sum
            const auto f1 = aggregate_class({rep1, rep2}, {0.5, 0.5});
            const auto f2 = aggregate_class({rep1b, rep2b}, {0.5, 0.5});
            const double mosaic = f1.lower + f2.lower;
            af << "mosaic_differs\t" << mosaic << '\t' << correct.lower << '\t'
               << (eq(mosaic, correct.lower) ? 0 : 1) << '\n';
            if (eq(mosaic, correct.lower)) ++abad;

            // ---- EXPOSURE ------------------------------------------------------------------
            const MassInterval ZERO{kNegInfBs2, kNegInfBs2};
            const double EA = 1000.0, EB = 600.0, LAM = 0.05;
            // Zero fragments must leave exactly -lambda(E_a + E_b).
            const auto z = representative_total(MassInterval{0.0, 0.0}, EA, EB, false, LAM, 0.0);
            emit("exposure_zero_frags", z.lower, -LAM * (EA + EB));
            // Adding fragments must not repeat the charge: the exposure term is identical.
            const auto withf = representative_total(MassInterval{-42.0, -41.0}, EA, EB, false, LAM, 0.0);
            emit("exposure_not_repeated", withf.lower - (-42.0), z.lower);
            // Homozygous exposure is -2*lambda*E_a, matching the 2*M_a chromosome-copy factor.
            const auto hz = representative_total(MassInterval{0.0, 0.0}, EA, EA, true, LAM, 0.0);
            emit("exposure_homozygous", hz.lower, -2.0 * LAM * EA);
            // Representatives of DIFFERENT length carry different exposure, so it must be added
            // before class aggregation -- equivalent output representatives are not interchangeable
            // at this step.
            const auto rA = representative_total(MassInterval{-10.0, -10.0}, EA, EA, true, LAM, 0.0);
            const auto rB = representative_total(MassInterval{-10.0, -10.0}, EB, EB, true, LAM, 0.0);
            af << "exposure_length_matters\t" << rA.lower << '\t' << rB.lower << '\t'
               << (eq(rA.lower, rB.lower) ? 0 : 1) << '\n';
            if (eq(rA.lower, rB.lower)) ++abad;

            // ---- CONDITION G --------------------------------------------------------------
            const std::string gp = bounded_search + ".certify.tsv";
            std::ofstream gf(gp);
            gf << "case\tn_plausible\tcertified\texpect_certified\tok\n";
            std::size_t gbad = 0;
            const auto gcase = [&](const char* nm, const std::vector<MassInterval>& cs,
                                   double tau, bool want_cert, std::size_t want_n) {
                const Certification c = certify(cs, tau);
                const bool good = (c.certified == want_cert) && (c.plausible.size() == want_n);
                if (!good) ++gbad;
                gf << nm << '\t' << c.plausible.size() << '\t' << (c.certified ? 1 : 0) << '\t'
                   << (want_cert ? 1 : 0) << '\t' << (good ? 1 : 0) << '\n';
            };
            // clear winner: its lower is above every other upper
            gcase("clear_winner", {{-10.0, -9.0}, {-20.0, -19.0}}, 0.0, true, 1);
            // point-estimate winner but OVERLAPPING intervals -> unresolved, not certified
            gcase("overlap_unresolved", {{-10.0, -5.0}, {-11.0, -4.0}}, 0.0, false, 2);
            // exact ties -> equivalence set
            gcase("exact_tie", {{-10.0, -10.0}, {-10.0, -10.0}}, 0.0, false, 2);
            // exact POINT intervals recover the ordinary argmax
            gcase("point_argmax", {{-10.0, -10.0}, {-12.0, -12.0}}, 0.0, true, 1);
            // tau widens the plausible set consistently
            gcase("tau_widens", {{-10.0, -10.0}, {-10.5, -10.5}}, 1.0, false, 2);
            // MONOTONICITY of tau: raising it can only ENLARGE the plausible set. Asserted rather
            // than assumed, because a sign error in the threshold would shrink it instead and every
            // single-tau case above would still pass.
            {
                const std::vector<MassInterval> cs{{-10.0, -10.0}, {-10.5, -10.5}, {-14.0, -14.0}};
                std::size_t prev = 0; bool mono = true;
                for (double t : {0.0, 0.25, 1.0, 5.0}) {
                    const std::size_t n = certify(cs, t).plausible.size();
                    if (n < prev) mono = false;
                    prev = n;
                }
                gf << "tau_monotone\t" << prev << "\t0\t0\t" << (mono ? 1 : 0) << '\n';
                if (!mono) ++gbad;
            }
            // ---- LAZY PRUNING, against the unpruned reference ------------------------------
            // The pruned arm is an OPTIMIZATION measured against the complete one; the unpruned
            // path stays runnable and is the reference in every case below.
            {
                // The fourth class is the DISCRIMINATING one: at tau=2, B=-10.0 so the plausible
                // threshold is -12.0. Its upper of -11.0 sits between -12.0 and B, so the correct
                // rule KEEPS it while the tighter U(C) < B rule would prune it. Without a class in
                // that band the unsafety of the tighter rule cannot be demonstrated -- measured,
                // the first version of this fixture reported 0 wrongly-pruned classes.
                const std::vector<MassInterval> cs{
                    {-10.0, -9.5}, {-10.4, -9.8}, {-13.0, -11.0},
                    {-14.0, -13.5}, {-30.0, -29.0}};
                for (double t : {0.0, 0.5, 2.0}) {
                    const Certification full = certify(cs, t);
                    // lazy: drop everything safely prunable at this tau, then certify the rest
                    std::vector<MassInterval> kept;
                    std::size_t pruned = 0;
                    for (const MassInterval& c : cs) {
                        if (safely_prunable(c, full.best_lower, t)) { ++pruned; continue; }
                        kept.push_back(c);
                    }
                    const Certification lazy = certify(kept, t);
                    const bool same = (lazy.plausible.size() == full.plausible.size()) &&
                                      (lazy.certified == full.certified);
                    // NON-VACUITY: pruning must actually remove something, or "agrees" is trivial.
                    const bool nonvac = pruned > 0;
                    gf << "lazy_tau" << t << '\t' << lazy.plausible.size() << '\t'
                       << (lazy.certified ? 1 : 0) << '\t' << (full.certified ? 1 : 0) << '\t'
                       << ((same && nonvac) ? 1 : 0) << '\n';
                    if (!same || !nonvac) ++gbad;
                }
                // A class pruned at tau=0 must stay pruned as B rises: B only ever increases during
                // refinement, so the threshold only ever moves away from a pruned class.
                const Certification c0 = certify(cs, 0.0);
                bool stays = true;
                for (const MassInterval& c : cs) {
                    if (!safely_prunable(c, c0.best_lower, 0.0)) continue;
                    if (!safely_prunable(c, c0.best_lower + 5.0, 0.0)) stays = false;
                }
                gf << "pruned_stays_pruned\t0\t0\t0\t" << (stays ? 1 : 0) << '\n';
                if (!stays) ++gbad;
                // And the tighter rule U(C) < B is UNSAFE: at tau > 0 it removes a class the
                // plausible set keeps. Asserted as a DIFFERENCE so the correction cannot regress.
                const double tt = 2.0;
                const Certification ft = certify(cs, tt);
                std::size_t bad_pruned = 0;
                for (std::size_t i = 0; i < cs.size(); ++i) {
                    const bool wrong_rule = cs[i].upper < ft.best_lower;          // no tau
                    const bool in_plausible =
                        std::find(ft.plausible.begin(), ft.plausible.end(), i) != ft.plausible.end();
                    if (wrong_rule && in_plausible) ++bad_pruned;
                }
                gf << "tighter_rule_is_unsafe\t" << bad_pruned << "\t0\t0\t"
                   << (bad_pruned > 0 ? 1 : 0) << '\n';
                if (bad_pruned == 0) ++gbad;
            }
            // THREE OUTCOMES, non-vacuously. The same overlapping case must read INCOMPLETE under a
            // restricted budget and CERTIFIED once refined, so "unresolved" can never be produced by
            // simply stopping early.
            {
                const std::vector<MassInterval> wide{{-10.0, -5.0}, {-11.0, -4.0}};
                const std::vector<MassInterval> tight{{-10.0, -9.99}, {-14.0, -13.99}};
                const std::vector<MassInterval> tie{{-10.0, -10.0}, {-10.0, -10.0}};
                struct VC { const char* nm; const std::vector<MassInterval>* cs; double tol;
                            Verdict want; };
                const VC vcs[] = {
                    {"verdict_incomplete", &wide,  0.5, Verdict::Incomplete},
                    {"verdict_certified",  &tight, 0.5, Verdict::Certified},
                    {"verdict_unresolved", &tie,   0.5, Verdict::Unresolved},
                };
                for (const VC& v : vcs) {
                    const Certification cc = certify(*v.cs, 0.0);
                    const Verdict got = verdict_of(*v.cs, cc, v.tol);
                    const bool good = got == v.want;
                    if (!good) ++gbad;
                    gf << v.nm << '\t' << cc.plausible.size() << '\t' << (cc.certified ? 1 : 0)
                       << "\t0\t" << (good ? 1 : 0) << '\n';
                }
            }
            gf.flush();
            log.wrote({gp});
            log.info("certification gates: " + std::to_string(gbad) + " failed");
            abad += gbad;

            af.flush();
            log.wrote({ap});
            log.info("aggregation gates: " + std::to_string(abad) + " failed");
        }

        }  // end --bounded-verify

        // HAPLOTYPE IDENTITY, exercised directly. States are enumerated per haplotype into separate
        // vectors, so removing `hap` from the key cannot collapse anything there -- every vector
        // holds one haplotype value. (An earlier mutation run reported this as "caught"; it was not,
        // something else had failed.) Here two otherwise-identical states from two haplotypes are
        // unioned, sorted and uniqued: the key must keep them apart.
        {
            std::vector<MatePlacement> f1{{100, 0}}, r1{}, f2{}, r2{{100 + ihi - 120, 0}};
            const auto s0 = enumerate_fragment_states(0, f1, r1, f2, r2, 120, 120, ilo, ihi);
            const auto s1 = enumerate_fragment_states(1, f1, r1, f2, r2, 120, 120, ilo, ihi);
            std::vector<FragmentState> both = s0;
            both.insert(both.end(), s1.begin(), s1.end());
            std::sort(both.begin(), both.end());
            both.erase(std::unique(both.begin(), both.end()), both.end());
            const std::string hp = bounded_search + ".hapkey.tsv";
            std::ofstream hf(hp);
            hf << "states_hap0\tstates_hap1\tunioned_unique\n"
               << s0.size() << '\t' << s1.size() << '\t' << both.size() << '\n';
            hf.flush();
            log.wrote({hp});
        }
        }
        sf.flush();
        log.wrote({fs_path});
        log.info("fragment states: " + std::to_string(tot_states) + " built, " +
                 std::to_string(st_dis) + " cell(s) disagree with exhaustive; orientation A " +
                 std::to_string(tot_a) + ", orientation B " + std::to_string(tot_b));
        log.info("in-band mass: " + std::to_string(mass_cells) + " cells, " +
                 std::to_string(mass_dis) + " bounded/exhaustive disagreements, " +
                 std::to_string(bound_viol) + " exceeding the untruncated reference");
        log.info("adaptive tail: " + std::to_string(adaptive_tight) + " of " +
                 std::to_string(mass_cells) + " cells within 1 nat, " +
                 std::to_string(adaptive_miss) + " not containing the reference");
        if (adaptive_miss != 0) {
            throw std::runtime_error("genotype-frag: the adaptive-tail interval does not contain "
                                     "the exact reference on " + std::to_string(adaptive_miss) +
                                     " cell(s)");
        }
        log.info("intervals: " + std::to_string(interval_miss) +
                 " cell(s) where the exact reference falls outside [lower, upper]");
        if (interval_miss != 0) {
            throw std::runtime_error("genotype-frag: the exact reference mass falls OUTSIDE the "
                                     "reported interval on " + std::to_string(interval_miss) +
                                     " cell(s); the omitted-mass bound does not bound");
        }
        if (mass_dis != 0) {
            throw std::runtime_error("genotype-frag: bounded and exhaustive in-band MASS differ on " +
                                     std::to_string(mass_dis) + " cell(s)");
        }
        if (bound_viol != 0) {
            throw std::runtime_error("genotype-frag: in-band mass EXCEEDS the untruncated reference "
                                     "on " + std::to_string(bound_viol) + " cell(s); the in-band sum "
                                     "cannot be larger than the integral that contains it");
        }
        if (st_dis != 0) {
            throw std::runtime_error("genotype-frag: bounded and exhaustive fragment-state sets "
                                     "differ on " + std::to_string(st_dis) + " cell(s)");
        }

        for (const Fragment& F : bfr) {
            for (int m = 0; m < 2; ++m) {
                const std::string& raw = m == 0 ? F.r1 : F.r2;
                if (raw.empty()) continue;
                for (int strand = 0; strand < 2; ++strand) {
                    const std::string r = strand == 0 ? raw : reverse_complement(raw);
                    // THE PRODUCTION BAND, from the shared helper. Re-deriving it here gave
                    // floor(div*len) and certified a band one edit NARROWER than production.
                    const std::size_t d = mate_band_edits(opt.max_divergence, r.size());
                    for (std::size_t h = 0; h < seqs.size(); ++h) {
                        SearchWork w;
                        const auto b = bounded_mate_placements(r, seqs[h], d, &w);
                        const auto e = exhaustive_mate_placements(r, seqs[h], d);
                        const bool same = b.size() == e.size() &&
                            std::equal(b.begin(), b.end(), e.begin(),
                                       [](const MatePlacement& x, const MatePlacement& y) {
                                           return x.start == y.start && x.edits == y.edits; });
                        if (!same) ++disagree;
                        tot_cand += w.candidate_starts; tot_ver += w.verified;
                        tot_exh += seqs[h].size() >= r.size() ? seqs[h].size() - r.size() + 1 : 0;
                        std::uint32_t best = 0xFFFFFFFFu;
                        for (const MatePlacement& p2 : b) best = std::min(best, p2.edits);
                        bf << F.name << '\t' << (m + 1) << '\t' << (strand == 0 ? '+' : '-')
                           << '\t' << names[h] << '\t' << b.size() << '\t' << e.size() << '\t'
                           << (same ? "yes" : "NO") << '\t'
                           << (best == 0xFFFFFFFFu ? -1 : static_cast<long>(best)) << '\t'
                           << d << '\t' << w.pieces << '\t' << w.candidate_starts << '\t' << w.distinct_starts
                           << '\t' << w.verified << '\t' << (w.exhaustive_fallback ? 1 : 0) << '\n';
                        ++rows;
                    }
                }
            }
        }
        bf.flush();
        if (!bf) throw std::runtime_error("genotype-frag: write failed for " + bounded_search);

        // ---- VALID-FR RULE, exercised directly ---------------------------------------------
        // Asserting the rule on synthetic coordinates rather than only through fragment counts:
        // a state-set comparison can agree while the rule is wrong in a way both sides share.
        {
            const std::string fr = bounded_search + ".fr.tsv";
            std::ofstream ff(fr);
            if (!ff) throw std::runtime_error("genotype-frag: cannot write " + fr);
            ff << "case\tfwd_start\trev_end\tinsert_lo\tinsert_hi\tinsert\tvalid\texpected\n";
            struct C { const char* name; long fs, re, lo, hi; bool want; };
            const C cases[] = {
                {"in_support",        100, 449, 200, 500, true},
                {"at_lo_edge",        100, 299, 200, 500, true},
                {"at_hi_edge",        100, 599, 200, 500, true},
                {"below_support",     100, 250, 200, 500, false},
                {"above_support",     100, 700, 200, 500, false},
                {"reverse_upstream",  500,  99, 200, 500, false},
                {"same_position",     100,  99, 200, 500, false},
            };
            std::size_t wrong = 0;
            for (const C& c : cases) {
                const bool v = valid_fr_coordinates(c.fs, c.re, c.lo, c.hi);
                if (v != c.want) ++wrong;
                ff << c.name << '\t' << c.fs << '\t' << c.re << '\t' << c.lo << '\t' << c.hi
                   << '\t' << (c.re - c.fs + 1) << '\t' << (v ? 1 : 0) << '\t' << (c.want ? 1 : 0)
                   << '\n';
            }
            ff.flush();
            log.wrote({fr});
            if (wrong != 0) {
                throw std::runtime_error("genotype-frag: the valid-FR rule is wrong on " +
                                         std::to_string(wrong) + " synthetic case(s)");
            }
            // The shared window must agree with the predicate, or the join and the test can drift.
            for (const C& c : cases) {
                const auto w2 = fr_reverse_end_window(c.fs, c.lo, c.hi);
                const bool inwin = c.re >= w2.first && c.re <= w2.second;
                if (inwin != valid_fr_coordinates(c.fs, c.re, c.lo, c.hi)) {
                    throw std::runtime_error("genotype-frag: fr_reverse_end_window disagrees with "
                                             "valid_fr_coordinates");
                }
            }
        }
        log.info("bounded search: " + std::to_string(rows) + " (mate,strand,haplotype) cells, " +
                 std::to_string(disagree) + " disagree with exhaustive; verified " +
                 std::to_string(tot_ver) + " starts against " + std::to_string(tot_exh) +
                 " exhaustive (" +
                 std::to_string(tot_exh ? 100.0 * static_cast<double>(tot_ver) /
                                          static_cast<double>(tot_exh) : 0.0) + "%)");
        log.wrote({bounded_search});
        if (disagree != 0) {
            throw std::runtime_error("genotype-frag: --bounded-search disagreed with the exhaustive "
                                     "reference on " + std::to_string(disagree) + " cell(s); the "
                                     "search is not complete within the band");
        }
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
              "\tpanel_domain_scope_blocks\tpanel_domain_scope\tachieved_bound\tbound_holds"
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
               << u.panel_domain_scope.size() << '\t';
            for (std::size_t q = 0; q < u.panel_domain_scope.size(); ++q) of << (q ? "," : "") << u.panel_domain_scope[q];
            if (u.panel_domain_scope.empty()) of << '.';
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
