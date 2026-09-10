#include "panvar/genotype_command.hpp"

#include "panvar/align.hpp"
#include "panvar/cli_utils.hpp"
#include "panvar/genotype_blocks.hpp"
#include "panvar/genotype_markers.hpp"
#include "panvar/genotype.hpp"
#include "panvar/md5.hpp"

#include "panvar/candidate_frame.hpp"
#include "panvar/genotype_fragments.hpp"
#include "panvar/genotype_index.hpp"
#include "panvar/graph_utils.hpp"
#include "panvar/genotype_reads.hpp"
#include "panvar/node_coverage.hpp"
#include "panvar/gfa.hpp"
#include "panvar/output.hpp"
#include "panvar/pangenie_model.hpp"

#include "panvar/parallel.hpp"
#include "panvar/syncmer.hpp"

#include <unordered_set>
#include <map>
#include <set>
#include <cstdio>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <atomic>
#include <chrono>
#include <fstream>
#include <sys/resource.h>
#include <iomanip>

#include <iostream>
#include <stdexcept>

namespace panvar {
namespace {

// Depth provenance, said out loud. A block that inherited the region's depth reads identically to one
// that measured its own, and at block_class=array that value is the denominator converting marker
// multiplicity into copy number -- so an inherited one carries none of that block's evidence into the
// number the block exists to produce.
void log_depth_provenance(cli::RunLog& log,
                          const std::vector<Block>& chain,
                          const std::vector<BlockDepth>& depth) {
    std::map<std::string, std::size_t> by_source;
    std::vector<std::size_t> fallback_bubbles;
    for (std::size_t bi = 0; bi < depth.size(); ++bi) {
        ++by_source[depth_source_name(depth[bi].source)];
        if (depth[bi].source == DepthSource::RegionFallback && bi < chain.size() &&
            chain[bi].kind == BlockKind::Bubble) {
            fallback_bubbles.push_back(chain[bi].bubble_id);
        }
    }
    std::string summary;
    for (const auto& [name, n] : by_source) {
        if (!summary.empty()) summary += ", ";
        summary += std::to_string(n) + " " + name;
    }
    log.info("depth provenance: " + summary);
    if (!fallback_bubbles.empty()) {
        std::string ids;
        for (std::size_t i = 0; i < fallback_bubbles.size() && i < 12; ++i) {
            if (!ids.empty()) ids += ",";
            ids += std::to_string(fallback_bubbles[i]);
        }
        if (fallback_bubbles.size() > 12) ids += ",...";
        log.info("  " + std::to_string(fallback_bubbles.size()) + " bubble block(s) have NO anchors of "
                 "their own and take the region's depth entirely (bubble " + ids + "); at a tandem "
                 "array that is the copy-number denominator");
    }
}

void print_genotype_help() {
    std::cout
        << "Usage:\n"
        << "  panvar genotype -i <graph.gfa> -b <bubble-prefix> -r <ref-path> -o <out-prefix>\n"
        << "                  -R <reads.fq.gz> [-R <reads2.fq.gz>]\n"
        << "\n"
        << "EXPERIMENTAL, and not part of the reviewed release. Build with\n"
        << "-DPANVAR_ENABLE_EXPERIMENTAL_GENOTYPE=ON. Interfaces and output columns may change.\n"
        << "\n"
        << "Genotype a short-read sample against a regional graph, per bubble. The locus is cut into an\n"
        << "alternating Flank / Bubble / Backbone chain, and a diploid Li-Stephens HMM over panel\n"
        << "haplotype pairs is run along it, with a negative-binomial emission over closed-syncmer\n"
        << "marker counts. Each block reports the pair of panel alleles the sample carries, a GQ, and\n"
        << "an independent read-derived length estimate (mass_bp).\n"
        << "\n"
        << "KNOWN LIMIT, at a tandem array: short reads cannot resolve the ORDER of repeat units --\n"
        << "measured, a fragment must span a unit junction to say anything about arrangement. Content\n"
        << "and total copy number are recoverable; arrangement is not. Read `mass_bp` as the copy-number\n"
        << "answer at `block_class=array`, not `called_bp`, which is quantised to panel allele lengths.\n"
        << "\n"
        << "--audit skips genotyping and reports only per-bubble marker feasibility: for every allele,\n"
        << "how many discriminative markers exist in each candidate unit -- syncmer nodes (present in\n"
        << "this allele, absent from its siblings) and syncmer adjacencies (the same rule on\n"
        << "consecutive-syncmer pairs, which is what carries a deletion junction).\n"
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
        << "                              Counts them against the block marker panel; without this the\n"
        << "                              run can only build an index or an --audit\n"
        << "      --min-anchors <N>       Below this many invariant markers a block is flagged\n"
        << "                              low_anchor (default 20). It no longer gates the local\n"
        << "                              estimate: every block with any anchors contributes,\n"
        << "                              shrunk toward the region by its own weight\n"
        << "      --marker-rule <r>       panvar (default) = keep markers whose multiplicity varies\n"
        << "                              across alleles; unique = carried by exactly one allele (any\n"
        << "                              copy number); pangenie = unique AND occurring once,\n"
        << "                              presence/absence only (PanGenie's actual rule);\n"
        << "                              mixed = presence/absence outside bubbles,\n"
        << "                              multiplicity inside them\n"
        << "      --fragment-len <N>      Library fragment length, used to discount correlated\n"
        << "                              markers when computing GQ (default 350; 0 disables)\n"
        << "      --depth-model <m>       How per-haplotype depth is estimated: joint (default,\n"
        << "                              a second pass that divides each block's anchors by how many\n"
        << "                              of its called alleles traverse it -- the only estimator that\n"
        << "                              is right when a haplotype's deletion removes whole blocks;\n"
        << "                              skipped entirely when no block has a bypass allele),\n"
        << "                              median (per-block anchor median/2 shrunk toward the region),\n"
        << "                              quantile (one region-wide value from a high quantile of the\n"
        << "                              block medians, so a deletion covering much of the locus\n"
        << "                              cannot drag it), or bases (total read bases over reference\n"
        << "                              length, independent of block structure)\n"
        << "      --depth-quantile <q>    Quantile for --depth-model quantile (default 0.75)\n"
        << "      --hybrid-call           Marker unaries plus certified fragment linkage. Linkage-\n"
        << "                              owned fragments leave the marker counts and feed edge\n"
        << "                              factors instead; the whole exchange is transactional, so\n"
        << "                              any refusal leaves the legacy call untouched.\n"
        << "      --hybrid-edges <path>   Per-edge sparse construction measurements, serialised\n"
        << "                              from the objects the caller builds.\n"
        << "      --hybrid-grouped        Build linkage edges through EXACT signature grouping:\n"
        << "                              alleles indistinguishable to this edge's fragments share\n"
        << "                              one delta, and the content-class table is never visited.\n"
        << "                              Alleles remain distinct HMM states throughout.\n"
        << "      --hybrid-grouped-report <path>  Per-edge grouped build measurements.\n"
        << "      --hybrid-edge-signature <a,b,path>  JOINT structural emission signatures on one\n"
        << "                              edge: per-fragment classes, their sum, and the joint\n"
        << "                              classes after intersecting every fragment's partition.\n"
        << "      --hybrid-context-dependence <ctx,b1,...,path>  Does block ctx affect ANY\n"
        << "                              fragment signature of the factor over b1..bk? Exact:\n"
        << "                              a state can touch ctx only via an accepted placement\n"
        << "                              there, so it enumerates those directly.\n"
        << "      --hybrid-wide-inventory <path>  Every Wide fragment's exact variable scope,\n"
        << "                              grouped, with the minimal factor interval each would\n"
        << "                              need. Wide evidence keeps the model INCOMPLETE, so this\n"
        << "                              says how many higher-order factors a COMPLETE call needs.\n"
        << "      --hybrid-super-ledger <b1,...,bk,path>  EVIDENCE ACCOUNTING for a candidate\n"
        << "                              super-factor over those blocks: which fragments would\n"
        << "                              enter it, each exactly once, and what else claims them.\n"
        << "      --hybrid-arity-probe <b1,...,bk,path>  Geometry over k consecutive variable\n"
        << "                              blocks: EVERY haploid window length by complete\n"
        << "                              enumeration, the affine boundary, and whether exposure\n"
        << "                              cancellation follows by algebra.\n"
        << "      --hybrid-triple-probe <...>  DEPRECATED alias for --hybrid-arity-probe.\n"
        << "      --hybrid-exposure-probe <path>  Per-configuration window length and exposure for\n"
        << "                              every candidate edge, with the affine regime boundary.\n"
        << "      --hybrid-edge-oracle <a,b,path>  Run ONE edge through BOTH the supported search\n"
        << "                              and the dense oracle and compare them cell by cell --\n"
        << "                              finite set, multiplicity and mass -- writing the full\n"
        << "                              counter set and the per-fragment proposal distribution.\n"
        << "      --hybrid-dry-run        Plan the transaction and report, then stop before calling.\n"
        << "      --hybrid-max-classes <n>  Operational budget: stored phase classes per edge.\n"
        << "      --hybrid-max-bytes <n>    Operational budget: measured bytes per edge.\n"
        << "      --hybrid-max-proposed-cells <n>  Operational budget: distinct allele-pair\n"
        << "                              cells retained by the support search, over the locus.\n"
        << "      --hybrid-max-full-read-verifications <n>  Operational budget: positional\n"
        << "                              starts subjected to whole-read Hamming verification.\n"
        << "                              Both are charged BEFORE the work, and exhaustion is a\n"
        << "                              model-level refusal, never a partial emission.\n"
        << "      --hybrid-max-verified-windows <n>  DEPRECATED alias for\n"
        << "                              --hybrid-max-proposed-cells (same quantity, honest name).\n"
        << "      --hybrid-max-dense-emission-windows <n>  Operational budget on DENSE emission\n"
        << "                              construction: fragments x n_A x n_B summed over edges.\n"
        << "                              Exceeding it gives INCOMPLETE:\n"
        << "                              dense-emission-window-limit, nothing subtracted.\n"
        << "      --hybrid-status <path>  Where the hybrid status report is written, including\n"
        << "                              every linkage parameter used.\n"
        << "      --hybrid-lambda-estimate  Estimate lambda candidate-independently as\n"
        << "                              N_fragments / (2 * median panel haplotype length). The\n"
        << "                              only option available on real data.\n"
        << "      --hybrid-lambda <x>     Fragment-start intensity for the linkage factors. NOT the\n"
        << "                              marker model's lambda_hap (different units), and never\n"
        << "                              fitted to the winning candidate. Default 0.05.\n"
        << "      --hybrid-bg-divergence <x>  Background per-base disagreement. Default 0.10.\n"
        << "      --hybrid-outlier-mix <x>    Background weight eta. Default 0.05.\n"
        << "      --hybrid-preflight <p>  Report candidate-frame coverage over the HMM's declared\n"
        << "                              state universe, before any hybrid wiring.\n"
        << "      --exclude-fragments <f> Fragment names (one per line) whose reads must NOT be\n"
        << "                              counted into the marker panel. Their occurrences are\n"
        << "                              subtracted; the markers themselves stay in the panel.\n"
        << "      --dump-markers <path>   Write every marker count, one row per (block, marker):\n"
        << "                              anchor or informative, k-mer GC, offset and clump, per-allele\n"
        << "                              multiplicity, and with --truth-haplotypes this sample's own\n"
        << "                              copy number and the count divided by it. Separating\n"
        << "                              efficiency from dosage needs that last pair\n"
        << "      --depth-estimator <e>   How anchor counts are reduced to one depth: median\n"
        << "                              (default, historical), mean, or trimmed (central 80%).\n"
        << "                              Anchor counts are small integers, so a median of them\n"
        << "                              is an integer however many are pooled -- the region\n"
        << "                              depth can only land on 11.0, 11.5, 12.0 and so on. At a\n"
        << "                              tandem array that is the copy-number denominator\n"
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
        << "      --model-pangenie        Genotype with a faithful port of PanGenie's model instead\n"
        << "                              of panvar's: their unique-once presence/absence marker rule,\n"
        << "                              their geometric+Poisson emission capped at copy number 2, and\n"
        << "                              their distance-scaled Li-Stephens transition. Same panel and\n"
        << "                              same read counts, so any difference is the model alone\n"
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
        << "      --all-kmers             Use every k-mer instead of the ~17% closed-syncmer sample\n"
        << "                              (~6x more markers, ~6x the memory and time)\n"
        << "      --no-region-unique      Skip region filtering entirely (diagnostic). This disables\n"
        << "                              BOTH rules -- markers varying in more than one block, and\n"
        << "                              markers the panel shows more often than the blocks account\n"
        << "                              for -- so it cannot attribute a loss to either one. It also\n"
        << "                              turns off the audit that counts them, which then reports\n"
        << "                              nothing rather than zero. Measured to break leave-zero-out\n"
        << "      --ledger-block <N>      Write every candidate marker of block N -- syncmer nodes and\n"
        << "                              2-syncmer adjacencies alike -- with the reason it was kept\n"
        << "                              or dropped, taken before the filter erases them\n"
        << "                              (<prefix>.blockN.ledger.tsv). Needs no reads: the ledger is\n"
        << "                              a property of the panel. One block, since the pre-filter set\n"
        << "                              is several times the retained one\n"
        << "      --dump-haplotype-alleles <f>  Write which allele every panel haplotype carries at\n"
        << "                              every block. The panel's basic fact, and the only way to\n"
        << "                              recover TRANSITIONS between blocks: per-block dumps cannot\n"
        << "                              say which allele pairs co-occur on one haplotype\n"
        << "      --max-alleles-block <B:N>  Cap block B at N candidates, overriding --max-alleles\n"
        << "                              there only. Raising the cap locus-wide to study one block\n"
        << "                              expands every other block too, so the neighbouring\n"
        << "                              emissions change and a difference at the block under study\n"
        << "                              cannot be attributed to it. Repeatable\n"
        << "      --oracle-called-only    With --certified-oracle, align ONLY the called pair rather\n"
        << "                              than searching all 2A. Reports called_total_edits; best_*\n"
        << "                              and excess_total_edits are NA. Use when the certified best\n"
        << "                              is already known from an earlier run -- the search is 914\n"
        << "                              alignments of ~50 kb at an array, about half an hour a sample\n"
        << "      --probe-pair <BLK:A,B>  Report the emission rank, tie count and delta for allele\n"
        << "                              pair (A,B) at block BLK instead of the truth's. The\n"
        << "                              certified oracle names the best REACHABLE pair, which at an\n"
        << "                              unrepresentable block is not the truth; this asks which\n"
        << "                              stage lost it. truth_rank -2 means pruned before scoring.\n"
        << "                              DIAGNOSTIC: the truth_* columns then describe the probe.\n"
        << "                              REPEATABLE: every occurrence adds a row to\n"
        << "                              <prefix>.probe_pairs.tsv (rank, ties, delta, raw score);\n"
        << "                              the FIRST also redirects the truth_* columns as above, so\n"
        << "                              adding a probe cannot change what an earlier one reported.\n"
        << "                              Probing costs nothing -- the emission matrix already exists\n"
        << "                              -- which is how a FIXED reference pair can be scored in\n"
        << "                              every arm of an experiment. Needed because `delta` is\n"
        << "                              measured against the block optimum and the optimum MOVES\n"
        << "                              when the counts do: deltas from different runs are not on a\n"
        << "                              common scale, but two pairs probed in the same run are,\n"
        << "                              since the block baseline is common to both and cancels\n"
        << "      --noiseless-counts <BLK|BLK:A,B>\n"
        << "                              Replace block BLK's observed marker counts with the counts a\n"
        << "                              pair would produce with NO read noise: lambda*(m1+m2)+mu at\n"
        << "                              every marker, rounded. With BLK alone the source is this\n"
        << "                              sample's truth haplotypes (needs --truth-haplotypes), which\n"
        << "                              under leave-one-out is OFF-PANEL -- that is the point, it is\n"
        << "                              what tests the emission's geometry on a projection. With\n"
        << "                              BLK:A,B it is panel alleles A and B, an on-panel control the\n"
        << "                              emission must rank first if it is a proper likelihood.\n"
        << "                              Only the TARGET BLOCK's counts are made noiseless: depth,\n"
        << "                              lambda and dispersion still come from the reads, and the\n"
        << "                              injected counts are built from that same lambda. So this\n"
        << "                              removes acquisition and sampling noise AT THAT BLOCK,\n"
        << "                              conditional on the estimated nuisance parameters -- it does\n"
        << "                              not make the run noiseless. Rejected\n"
        << "                              with --edge-weight or --evidence other than syncmer, since\n"
        << "                              adjacency and coverage counts are NOT synthesized and mixing\n"
        << "                              them with noiseless nodes would attribute nothing.\n"
        << "                              DIAGNOSTIC: the run no longer describes these reads\n"
        << "      --noiseless-scope <all|present|absent>\n"
        << "                              Which half of the block's markers --noiseless-counts\n"
        << "                              rewrites: those the source pair carries, those it does not,\n"
        << "                              or both (default all). The other half keeps its observed\n"
        << "                              counts. Splitting it is what makes the diagnosis causal --\n"
        << "                              if `absent` alone reproduces the damage, the defect is in\n"
        << "                              how unsupported markers are scored. Note `absent` writes 0,\n"
        << "                              not mu: counts are integers and mu is well under 1, so this\n"
        << "                              is the SHARPEST possible absent reading, sharper than real\n"
        << "                              reads which scatter stray counts there\n"
        << "      --max-linkage-emission-loss <F|inf>\n"
        << "                              How much block-local emission a state may give up and\n"
        << "                              still be reachable by the chain. States losing more than F\n"
        << "                              to the block optimum are excluded BEFORE forward-backward,\n"
        << "                              so posterior and GQ are computed under the constraint.\n"
        << "                              Default inf = unrestricted. 0 is NOT off: it admits states\n"
        << "                              TIED with the optimum, so linkage can still resolve a tie.\n"
        << "                              Measured over 6 loci x 20 donors, linkage moved off a\n"
        << "                              unique local optimum 93 times: 20 rescues (median move\n"
        << "                              0.15) against 73 overrides (median 1.57)\n"
        << "      --edge-weight <F>       Weight on 2-syncmer adjacency evidence (default 0 = off).\n"
        << "                              Adjacencies come from the same reads as the nodes, so this\n"
        << "                              double-counts; measured worth 0.5-4% at 0.25-0.5, and 1.0\n"
        << "                              helps one locus while hurting another. Cannot be combined\n"
        << "                              with --compositional (different likelihood scales)\n"
        << "      --max-dense-alleles <N> Above N alleles in a block, score pairwise separation with\n"
        << "                              an O(n) sparse accumulator instead of the O(n^2) dense one\n"
        << "                              (default 2048; identical results, ~40% slower, but the\n"
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
    std::map<std::size_t, std::size_t> max_alleles_override;
    double fragment_len = 350.0;
    bool fragment_len_set = false;
    double recomb_rate = 1.0;
    double carrier_weight = 0.0;
    DepthModel depth_model = DepthModel::Joint;
    DepthEstimator depth_estimator = DepthEstimator::Median;
    std::string dump_markers;
    std::string exclude_fragments_path;
    std::string hybrid_preflight;
    bool hybrid_call = false;
    std::string hybrid_geometry_probe;
    std::string hybrid_orientation_probe;
    HybridLinkageParameters hyb_params;
    // OPERATIONAL limits: bytes and counts, never a likelihood or score. A resource refusal must be
    // attributable to the machine, not to which genotype looked good.
    std::size_t hybrid_max_classes = 50000000;
    std::size_t hybrid_max_bytes = 2000000000;
    // DENSE-EMISSION window alignments summed over edges: fragments x n_A x n_B. Named for the
    // construction it bounds, so that when the bounded-complete support search lands it reports its
    // OWN counters -- seed hits, proposed states, verified states -- and this budget cannot silently
    // acquire a new meaning. C4 needs 1,154,642 and does not finish, so the default refuses it.
    std::size_t hybrid_max_dense_emission_windows = 200000;
    // WINDOWS ACTUALLY VERIFIED, summed over edges -- the work production does, as opposed to the
    // dense cross-product it avoids. Its own name and its own budget, so neither figure can be
    // mistaken for the other.
    // SEPARATE OPERATIONAL COUNTERS. "verified windows" named a quantity that, after direct
    // positional verification, is neither verified nor a window.
    std::uint64_t hybrid_max_proposed_cells = 5000000;
    std::uint64_t hybrid_max_full_read_verifications = 200000000;
    bool hybrid_verified_windows_alias_used = false;
    std::size_t verified_windows_total = 0;
    std::uint64_t hyb_work_proposed = 0, hyb_work_verifications = 0, hyb_work_bases = 0;
    std::string hybrid_edges_path;
    bool hybrid_dry_run = false;
    std::string hybrid_edge_oracle;      // "<a>,<b>:<path>" -- one edge, both paths, compared
    bool hybrid_grouped = false;         // build edges through exact signature grouping
    std::string hybrid_grouped_report;   // per-edge grouped build measurements
    std::string hybrid_edge_signature;   // "<a>,<b>,<path>" -- joint emission signatures
    std::string hybrid_exposure_probe;   // per-configuration window length and exposure
    std::string hybrid_triple_probe;     // "<b1>,...,<bk>,<path>" -- k-variable geometry
    std::string hybrid_super_ledger;     // "<b1>,...,<bk>,<path>" -- evidence accounting
    std::string hybrid_wide_inventory;   // every Wide fragment's scope, grouped
    std::string hybrid_context_dependence;  // "<ctx>,<b1>,...,<path>" -- is ctx emission-relevant?
    std::string hybrid_interval_probe;      // "<b1>,...,<path>" -- the generic k-block geometry
    // REPEATABLE. Two overlapping factors have to exist AT THE SAME TIME for the higher-order
    // recurrence to be run on them, and building them in separate processes cannot show that.
    std::vector<std::string> hybrid_factor_runs;   // each "<b1>,...,<path>"
    std::string hybrid_higher_bench;               // where to write the recurrence benchmark
    std::uint64_t hybrid_max_message_entries = 0;  // 0 = unbounded; exceeding it is a REFUSAL
    std::size_t hybrid_plan_threads = 1;
    bool hybrid_plan_only = false;                 // report requirements, allocate nothing
    bool hybrid_higher = false;                    // route the call through the higher-order chain
    std::string hybrid_higher_report;              // realised counts, written AFTER the call
    // THE BLOCK CATALOGUE, for an evaluator that must not reconstruct blocks from raw graph
    // coordinates. Emits exactly what the caller itself indexes: each block's allele sequences, and
    // each panel path's allele index per block. An accuracy harness can then take truth and called
    // sequences from the same projection the call was made in, rather than a parallel one.
    std::string dump_block_catalogue;
    // COMPATIBILITY ONLY. The model's insert floor is max(|r1|, |r2|); this restores |r1| + |r2|,
    // which declares every overlapping pair impossible.
    bool hybrid_no_overlap_pairs = false;
    bool hybrid_factor_oracle = false;      // also run the exhaustive oracle and compare
    // REPEATABLE AND POSITIONAL, matched to --hybrid-factor-run in order. Two factors do not
    // supersede the same edges -- F1 over {2,3,4,5} replaces 3-4 and 4-5, F2 over {4,5,6}
    // replaces nothing -- so one global list cannot describe both, and applying F1's list to F2
    // would hand F2 evidence that belongs to F1.
    std::vector<std::string> hybrid_factor_supersedes;   // "a-b,c-d" per factor run
    bool hybrid_triple_alias_used = false;
    struct EdgeRow {
        std::uint32_t a = 0, b = 0;
        std::size_t na = 0, nb = 0, owned = 0, informative = 0, support = 0, stored = 0,
                    predicted = 0, theoretical = 0, bytes = 0;
        double build_s = 0.0;
        std::string status;
        // The SUPPORT SEARCH's own counters, distinct from the dense baseline.
        std::size_t seed_hits = 0, proposed = 0, verified = 0, dense_windows = 0, fallbacks = 0;
        // Direct verification's own work, which is what production now pays: whole reads checked
        // at seeded starts, and the emission cells that came out finite.
        std::size_t read_verifications = 0, finite_cells = 0;
    };
    std::vector<EdgeRow> edge_rows;
    std::string hybrid_status_path;
    double depth_quantile = 0.75;
    long dump_block = -1;
    long ledger_block = -1;
    std::string alleles_out;
    // Infinity = unrestricted, and it is the default so current output cannot move. NOT 0: tau of 0
    // admits every state tied with the block optimum, which is a real and useful setting.
    double max_linkage_loss = std::numeric_limits<double>::infinity();
    long probe_block = -1; long probe_a = -1, probe_b = -1;
    // --noiseless-counts. `noiseless_a/b` stay -1 when the source is the sample's own truth.
    long noiseless_block = -1; long noiseless_a = -1, noiseless_b = -1;
    // Which half of the marker set the injection rewrites. Splitting it is what turns "noiseless is
    // worse" into a causal statement: if rewriting only the markers the source LACKS reproduces the
    // damage, the defect is in how absent markers are scored.
    std::string noiseless_scope = "all";       // all | present | absent
    // Where to write the source pair's multiplicity vector for the injected block. It is the only
    // quantity in the decomposition that cannot be recovered from --dump-block: under leave-one-out
    // the truth is off-panel, so no allele column carries it.
    std::string injection_out;
    // Whole-locus marker multiplicity for the truth haplotypes. --dump-injection and
    // --noiseless-counts both measure the truth's marker content WITHIN one block; this counts every
    // panel marker across the truth's ENTIRE walk, which is what says whether a block decomposition
    // has missed places that generate a marker's counts.
    std::string truth_markers_out;
    // Every --probe-pair, in order. The first also redirects the truth_* columns, which is what the
    // single-probe form has always done; all of them appear in the probe table.
    std::vector<std::array<std::size_t, 3>> probe_pairs;
    bool oracle_called_only = false;
    // Extra allele pairs to price against the truth, beyond the certified best and the called pair.
    // The oracle's cost is the 2A alignment cache; once it exists, any pair is two lookups. Without
    // this, asking "how far from truth is the pair this new score picked" means re-running the whole
    // 914-alignment search.
    std::vector<std::array<std::size_t, 2>> oracle_pairs;
    bool depth_calibration = false;
    double mass_weight = 0.0;
    bool nearest_rank = false;
    bool oracle_rank = false;
    bool certified_oracle = false;
    std::size_t certified_oracle_max = 1024;
    long certified_oracle_block = -1;
    std::string explain_pair;
    double marker_outlier = 0.0;
    bool restore_stripped = false;
    bool compositional = false;
    double edge_weight = 0.0;
    double robust_c = 0.0;
    double mass_window = 0.0;
    double scale_weight = 1.0;
    long deconvolve = -1;
    long cosine_block = -1;
    bool node_coverage = false;
    long coverage_block = -1;
    std::string coverage_preset = "sr";
    std::string evidence = "syncmer";
    bool model_pangenie = false;
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
        else if (arg == "--max-alleles-block") {
            const std::string v = require_value(arg);
            const std::size_t c = v.find(':');
            if (c == std::string::npos || v.find(':', c + 1) != std::string::npos)
                throw std::runtime_error("--max-alleles-block wants exactly BLOCK:N (e.g. 13:512)");
            max_alleles_override[cli::parse_size_arg(arg, v.substr(0, c))] =
                cli::parse_size_arg(arg, v.substr(c + 1));
        }
        else if (arg == "--fragment-len") { fragment_len = std::stod(require_value(arg)); fragment_len_set = true; }
        else if (arg == "--recomb-rate") recomb_rate = std::stod(require_value(arg));
        else if (arg == "--depth-model") {
            const std::string v = require_value(arg);
            if (v == "median") depth_model = DepthModel::Median;
            else if (v == "quantile") depth_model = DepthModel::Quantile;
            else if (v == "bases") depth_model = DepthModel::Bases;
            else if (v == "joint") depth_model = DepthModel::Joint;
            else throw std::runtime_error("genotype: --depth-model must be median|quantile|bases|joint");
        }
        else if (arg == "--depth-quantile") depth_quantile = std::stod(require_value(arg));
        else if (arg == "--exclude-fragments") exclude_fragments_path = require_value(arg);
        else if (arg == "--hybrid-preflight") hybrid_preflight = require_value(arg);
        else if (arg == "--hybrid-call") hybrid_call = true;
        else if (arg == "--hybrid-geometry-probe") hybrid_geometry_probe = require_value(arg);
        else if (arg == "--hybrid-orientation-probe") hybrid_orientation_probe = require_value(arg);
        else if (arg == "--hybrid-status") hybrid_status_path = require_value(arg);
        else if (arg == "--hybrid-fragment-sd") hyb_params.fragment_sd = std::stod(require_value(arg));
        else if (arg == "--hybrid-divergence") hyb_params.max_divergence = std::stod(require_value(arg));
        else if (arg == "--hybrid-error-rate") hyb_params.error_rate = std::stod(require_value(arg));
        else if (arg == "--hybrid-bg-divergence") hyb_params.bg_divergence = std::stod(require_value(arg));
        else if (arg == "--hybrid-outlier-mix") hyb_params.outlier_mix = std::stod(require_value(arg));
        else if (arg == "--hybrid-max-classes") hybrid_max_classes = std::stoull(require_value(arg));
        else if (arg == "--hybrid-max-bytes") hybrid_max_bytes = std::stoull(require_value(arg));
        else if (arg == "--hybrid-max-proposed-cells")
            hybrid_max_proposed_cells = std::stoull(require_value(arg));
        else if (arg == "--hybrid-max-full-read-verifications")
            hybrid_max_full_read_verifications = std::stoull(require_value(arg));
        else if (arg == "--hybrid-max-verified-windows") {
            // DEPRECATED ALIAS, kept because its old meaning -- allele-pair cells retained -- is
            // exactly what --hybrid-max-proposed-cells now names. Accepted rather than silently
            // reinterpreted, and it says so.
            hybrid_max_proposed_cells = std::stoull(require_value(arg));
            hybrid_verified_windows_alias_used = true;
        }
        else if (arg == "--hybrid-max-dense-emission-windows")
            hybrid_max_dense_emission_windows = std::stoull(require_value(arg));
        else if (arg == "--hybrid-edges") hybrid_edges_path = require_value(arg);
        else if (arg == "--hybrid-dry-run") hybrid_dry_run = true;
        else if (arg == "--hybrid-edge-oracle") hybrid_edge_oracle = require_value(arg);
        else if (arg == "--hybrid-grouped") hybrid_grouped = true;
        else if (arg == "--hybrid-grouped-report") hybrid_grouped_report = require_value(arg);
        else if (arg == "--hybrid-edge-signature") hybrid_edge_signature = require_value(arg);
        else if (arg == "--hybrid-exposure-probe") hybrid_exposure_probe = require_value(arg);
        else if (arg == "--hybrid-arity-probe") hybrid_triple_probe = require_value(arg);
        else if (arg == "--hybrid-super-ledger") hybrid_super_ledger = require_value(arg);
        else if (arg == "--hybrid-wide-inventory") hybrid_wide_inventory = require_value(arg);
        else if (arg == "--hybrid-context-dependence")
            hybrid_context_dependence = require_value(arg);
        else if (arg == "--hybrid-interval-probe") hybrid_interval_probe = require_value(arg);
        else if (arg == "--hybrid-factor-run") hybrid_factor_runs.push_back(require_value(arg));
        else if (arg == "--hybrid-higher-bench") hybrid_higher_bench = require_value(arg);
        else if (arg == "--hybrid-max-message-entries")
            hybrid_max_message_entries = std::stoull(require_value(arg));
        else if (arg == "--hybrid-plan-threads")
            hybrid_plan_threads = static_cast<std::size_t>(std::stoul(require_value(arg)));
        else if (arg == "--hybrid-plan-only") hybrid_plan_only = true;
        else if (arg == "--hybrid-higher") {
            hybrid_higher = true;
            // GROUPED CONSTRUCTION IS NOT OPTIONAL HERE. Ungrouped, edge 6-7 predicts 1.97e8
            // classes against a 5e7 cap and refuses on resources; grouped, the same edge is a few
            // thousand. Leaving that to a second flag makes a forgotten switch look like a
            // modelling limit, which is exactly what it did.
            hybrid_grouped = true;
        }
        else if (arg == "--hybrid-higher-report") hybrid_higher_report = require_value(arg);
        else if (arg == "--dump-block-catalogue") dump_block_catalogue = require_value(arg);
        else if (arg == "--no-overlap-pairs") hybrid_no_overlap_pairs = true;
        else if (arg == "--hybrid-factor-oracle") hybrid_factor_oracle = true;
        else if (arg == "--hybrid-factor-supersede")
            hybrid_factor_supersedes.push_back(require_value(arg));
        else if (arg == "--hybrid-triple-probe") {
            // DEPRECATED: the probe takes any number of consecutive blocks now, so "triple" names
            // a case rather than the feature. Accepted, and said out loud.
            hybrid_triple_probe = require_value(arg);
            hybrid_triple_alias_used = true;
        }
        else if (arg == "--hybrid-lambda-estimate") {
            hyb_params.lambda_source = HybridLinkageParameters::LambdaSource::Estimated;
        }
        else if (arg == "--hybrid-lambda") {
            // FRAGMENT-START INTENSITY, supplied. Never derived from the marker model's lambda_hap
            // (different units) and never fitted to the winning candidate.
            hyb_params.lambda = std::stod(require_value(arg));
            hyb_params.lambda_source = HybridLinkageParameters::LambdaSource::Supplied;
        }
        else if (arg == "--dump-markers") dump_markers = require_value(arg);
        else if (arg == "--dump-anchors") dump_markers = require_value(arg);  // former name
        else if (arg == "--depth-estimator") {
            const std::string v = require_value(arg);
            if (v == "median") depth_estimator = DepthEstimator::Median;
            else if (v == "mean") depth_estimator = DepthEstimator::Mean;
            else if (v == "trimmed") depth_estimator = DepthEstimator::TrimmedMean;
            else throw std::runtime_error("genotype: --depth-estimator must be median|mean|trimmed");
        }
        else if (arg == "--dump-block") dump_block = std::stol(require_value(arg));
        else if (arg == "--ledger-block") ledger_block = std::stol(require_value(arg));
        else if (arg == "--dump-haplotype-alleles") alleles_out = require_value(arg);
        else if (arg == "--probe-pair") {
            // BLOCK:A,B -- point the emission-rank diagnostics at a chosen allele pair instead of
            // the truth's. The certified oracle names the best REACHABLE pair, which at an
            // unrepresentable block is not the truth and so has no rank reported; this is how to
            // ask which stage lost it. truth_rank == -2 already means "pruned before scoring".
            const std::string v = require_value(arg);
            const std::size_t c = v.find(':'), m = v.find(',');
            if (c == std::string::npos || m == std::string::npos || m < c ||
                v.find(':', c + 1) != std::string::npos || v.find(',', m + 1) != std::string::npos) {
                throw std::runtime_error("--probe-pair wants exactly BLOCK:A,B (e.g. 13:271,304)");
            }
            const auto whole = [&](const std::string& t) {
                if (t.empty()) throw std::runtime_error("--probe-pair has an empty field in '" + v + "'");
                std::size_t used = 0;
                long n = 0;
                try {
                    n = std::stol(t, &used);
                } catch (const std::exception&) {
                    throw std::runtime_error("--probe-pair: '" + t + "' is not a whole number");
                }
                if (used != t.size())
                    throw std::runtime_error("--probe-pair: '" + t + "' is not a whole number");
                if (n < 0) throw std::runtime_error("--probe-pair indices must be >= 0");
                return n;
            };
            const long pb = whole(v.substr(0, c));
            const long pa = whole(v.substr(c + 1, m - c - 1));
            const long pbb = whole(v.substr(m + 1));
            // Repeatable. The first occurrence keeps the historical behaviour of pointing the truth_*
            // columns at the pair; later ones only add rows to the probe table, so adding a second
            // probe cannot silently change what the first one reported.
            if (probe_block < 0) { probe_block = pb; probe_a = pa; probe_b = pbb; }
            probe_pairs.push_back({static_cast<std::size_t>(pb), static_cast<std::size_t>(pa),
                                   static_cast<std::size_t>(pbb)});
        }
        else if (arg == "--noiseless-counts") {
            // BLOCK, or BLOCK:A,B -- replace this block's observed marker counts with the counts a
            // stated allele pair would produce with no read noise, then run the ordinary pipeline.
            // The point is to score the SAME emission on data it cannot blame: if it still prefers
            // the wrong pair, no amount of read work can help.
            const std::string v = require_value(arg);
            const std::size_t c = v.find(':');
            const auto whole = [&](const std::string& t) {
                if (t.empty())
                    throw std::runtime_error("--noiseless-counts has an empty field in '" + v + "'");
                std::size_t used = 0;
                long n = 0;
                try {
                    n = std::stol(t, &used);
                } catch (const std::exception&) {
                    throw std::runtime_error("--noiseless-counts: '" + t + "' is not a whole number");
                }
                if (used != t.size())
                    throw std::runtime_error("--noiseless-counts: '" + t + "' is not a whole number");
                if (n < 0) throw std::runtime_error("--noiseless-counts indices must be >= 0");
                return n;
            };
            if (c == std::string::npos) {
                noiseless_block = whole(v);          // source = this sample's truth haplotypes
            } else {
                const std::size_t m = v.find(',');
                if (m == std::string::npos || m < c ||
                    v.find(':', c + 1) != std::string::npos ||
                    v.find(',', m + 1) != std::string::npos) {
                    throw std::runtime_error("--noiseless-counts wants BLOCK or BLOCK:A,B "
                                             "(e.g. 13, or 13:271,304)");
                }
                noiseless_block = whole(v.substr(0, c));
                noiseless_a = whole(v.substr(c + 1, m - c - 1));
                noiseless_b = whole(v.substr(m + 1));
            }
        }
        else if (arg == "--dump-injection") injection_out = require_value(arg);
        else if (arg == "--dump-truth-marker-counts") truth_markers_out = require_value(arg);
        else if (arg == "--noiseless-scope") {
            noiseless_scope = require_value(arg);
            if (noiseless_scope != "all" && noiseless_scope != "present" && noiseless_scope != "absent")
                throw std::runtime_error("genotype: --noiseless-scope must be all|present|absent");
        }
        else if (arg == "--max-linkage-emission-loss") {
            const std::string v = require_value(arg);
            if (v == "inf" || v == "INF" || v == "infinity") {
                max_linkage_loss = std::numeric_limits<double>::infinity();
            } else {
                max_linkage_loss = std::stod(v);
                if (!(max_linkage_loss >= 0.0)) {
                    throw std::runtime_error(arg + " must be >= 0 or 'inf'; a negative loss would "
                                             "exclude the block optimum itself");
                }
            }
        }
        else if (arg == "--depth-calibration") depth_calibration = true;
        else if (arg == "--mass-weight") mass_weight = std::stod(require_value(arg));
        else if (arg == "--nearest-emission-rank") nearest_rank = true;
        else if (arg == "--oracle-emission-rank") oracle_rank = true;
        else if (arg == "--certified-oracle") certified_oracle = true;
        // The 2A search is the whole cost: 914 alignments of ~50 kb sequences at LPA's array is
        // half an hour per sample. When the certified best is already known from an earlier run,
        // only the CALLED pair needs aligning, which is 2.
        else if (arg == "--oracle-called-only") oracle_called_only = true;
        else if (arg == "--oracle-pair") {
            const std::string v = require_value(arg);
            const std::size_t c = v.find(',');
            if (c == std::string::npos || v.find(',', c + 1) != std::string::npos)
                throw std::runtime_error("--oracle-pair wants exactly A,B (e.g. 271,304)");
            const auto whole = [&](const std::string& t) {
                if (t.empty()) throw std::runtime_error("--oracle-pair has an empty field in '" + v + "'");
                std::size_t used = 0; long n = 0;
                try { n = std::stol(t, &used); }
                catch (const std::exception&) {
                    throw std::runtime_error("--oracle-pair: '" + t + "' is not a whole number"); }
                if (used != t.size())
                    throw std::runtime_error("--oracle-pair: '" + t + "' is not a whole number");
                if (n < 0) throw std::runtime_error("--oracle-pair indices must be >= 0");
                return static_cast<std::size_t>(n);
            };
            oracle_pairs.push_back({whole(v.substr(0, c)), whole(v.substr(c + 1))});
        }
        else if (arg == "--certified-oracle-block")
            certified_oracle_block = std::stol(require_value(arg));
        else if (arg == "--certified-oracle-max-alleles")
            certified_oracle_max = static_cast<std::size_t>(std::stoul(require_value(arg)));
        else if (arg == "--explain-pair") explain_pair = require_value(arg);
        else if (arg == "--marker-outlier") marker_outlier = std::stod(require_value(arg));
        else if (arg == "--restore-stripped-alleles") restore_stripped = true;
        else if (arg == "--compositional") compositional = true;
        else if (arg == "--edge-weight") {
            edge_weight = std::stod(require_value(arg));
            if (!(edge_weight >= 0.0) || !std::isfinite(edge_weight))
                throw std::runtime_error("--edge-weight must be a finite value >= 0");
        }
        else if (arg == "--robust") robust_c = std::stod(require_value(arg));
        else if (arg == "--mass-window") mass_window = std::stod(require_value(arg));
        else if (arg == "--scale-weight") scale_weight = std::stod(require_value(arg));
        else if (arg == "--deconvolve") deconvolve = std::stol(require_value(arg));
        else if (arg == "--cosine-block") cosine_block = std::stol(require_value(arg));
        else if (arg == "--node-coverage") node_coverage = true;
        else if (arg == "--coverage-block") coverage_block = std::stol(require_value(arg));
        else if (arg == "--coverage-preset") coverage_preset = require_value(arg);
        else if (arg == "--evidence") {
            evidence = require_value(arg);
            if (evidence != "syncmer" && evidence != "coverage" && evidence != "auto" &&
                evidence != "combined") {
                throw std::runtime_error("genotype: --evidence must be syncmer|coverage|auto|combined");
            }
        }
        else if (arg == "--model-pangenie") model_pangenie = true;
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
    // AFTER the parse loop, not before it. Assigned at declaration this reads the default and the flag
    // never reaches the marker builder -- the exact failure this module has hit repeatedly.
    if (fragment_len < 0.0) {
        throw std::runtime_error("genotype: --fragment-len must be >= 0");
    }
    {
        options.fragment_len = fragment_len;
    }
    if (!oracle_pairs.empty() && !certified_oracle) {
        throw std::runtime_error("genotype: --oracle-pair prices a pair using the certified oracle's "
                                 "alignment cache; it does nothing without --certified-oracle");
    }
    if (oracle_called_only && !certified_oracle) {
        throw std::runtime_error("genotype: --oracle-called-only only means anything with "
                                 "--certified-oracle; on its own it reports nothing");
    }
    // A diagnostic that quietly does nothing is worse than one that refuses: the run looks like an
    // arm of the experiment and is actually the baseline, which is exactly how a null result gets
    // manufactured. Both combinations below were silently inert.
    if (!truth_markers_out.empty() && truth_haplotypes.empty()) {
        throw std::runtime_error("genotype: --dump-truth-marker-counts counts panel markers across "
                                 "the TRUTH haplotypes' walks; pass --truth-haplotypes");
    }
    if (!injection_out.empty() && noiseless_block < 0) {
        throw std::runtime_error("genotype: --dump-injection writes the vector --noiseless-counts "
                                 "injects; without it there is nothing to write");
    }
    if (noiseless_scope != "all" && noiseless_block < 0) {
        throw std::runtime_error("genotype: --noiseless-scope selects which half of the markers "
                                 "--noiseless-counts rewrites; on its own it rewrites nothing");
    }
    if (!index_in.empty() && (noiseless_block >= 0 || !probe_pairs.empty())) {
        throw std::runtime_error("genotype: --probe-pair, --noiseless-counts and --noiseless-scope "
                                 "are not implemented on the indexed route (--index) and would be "
                                 "ignored; run the direct route with -i/--gfa");
    }
    if (noiseless_block >= 0) {
        // It replaces counts, so there must be counts to replace: the panel, the depth model and phi
        // are all built from the reads and only the named block's node counts are overwritten.
        if (read_paths.empty()) {
            throw std::runtime_error("genotype: --noiseless-counts replaces the counts at one block "
                                     "and needs the rest of the run to be real; pass -R/--reads");
        }
        if (noiseless_a < 0 && truth_haplotypes.empty()) {
            throw std::runtime_error("genotype: --noiseless-counts <BLK> takes its counts from this "
                                     "sample's truth; pass --truth-haplotypes, or name a panel pair "
                                     "with --noiseless-counts <BLK:A,B>");
        }
        if (edge_weight != 0.0) {
            throw std::runtime_error("genotype: --noiseless-counts does not synthesize adjacency "
                                     "counts, so with --edge-weight the emission would mix noiseless "
                                     "nodes with observed edges and attribute nothing");
        }
        if (evidence != "syncmer") {
            throw std::runtime_error("genotype: --noiseless-counts does not synthesize alignment "
                                     "coverage, so with --evidence " + evidence + " the emission "
                                     "would mix noiseless counts with observed coverage");
        }
    }
    if (!audit && !audit_linkage && read_paths.empty() && index_out.empty() && ledger_block < 0 &&
        alleles_out.empty()) {
        throw std::runtime_error(
            "genotype: pass --audit, --audit-linkage, --build-index, --ledger-block, "
            "--dump-haplotype-alleles, or -R/--reads");
    }
    // Reads need the block chain and the marker panel, but not the linkage/novelty audits or the
    // separation statistics -- those are diagnostics and were roughly half the runtime.
    const bool need_blocks =
        audit_linkage || !read_paths.empty() || !index_out.empty() || ledger_block >= 0 ||
        !alleles_out.empty();
    // Every option the model takes, in one place, so the indexed and direct paths cannot diverge.
    auto make_genotype_options = [&]() {
        GenotypeOptions g;
        g.threads = options.threads;
        g.max_alleles_per_block = max_alleles;
        g.max_alleles_override = max_alleles_override;
        g.max_linkage_emission_loss = max_linkage_loss;
        g.fragment_len = fragment_len;
        g.provenance = provenance;
        g.recomb_rate = recomb_rate;
        g.carrier_weight = carrier_weight;
        g.mass_weight = mass_weight;
        g.marker_outlier = marker_outlier;
        g.compositional = compositional;
        // The compositional emission is a multinomial SHAPE score; the adjacency term is a negative
        // binomial over absolute counts. Adding them is the incommensurable-scale error that already
        // cost this module once, and the edge term currently reaches only the non-compositional
        // branch, so combining them would silently do nothing. Refuse rather than mislead.
        if (edge_weight > 0.0 && compositional)
            throw std::runtime_error("--edge-weight cannot be combined with --compositional: the "
                                     "compositional emission is a multinomial shape score and the "
                                     "adjacency term is a negative binomial over counts, so the two "
                                     "are not on the same scale");
        g.edge_weight = edge_weight;
        g.scale_weight = scale_weight;
        g.robust_c = robust_c;
        g.mass_window = mass_window;
        return g;
    };
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
        // The index's clumps were computed at ONE fragment length. Building the emission from a
        // different runtime value mixes the two: the clump count reflects the index's length while
        // every other term uses the caller's. Inherit when unspecified, refuse when contradicted --
        // silently preferring either would reintroduce exactly the indexed-versus-direct drift v5
        // exists to remove.
        {
            // No `> 0` guard: 0 is a MEANINGFUL value, it disables clumping. Skipping inheritance
            // there sent an index built at 0 back to the 350 default and produced GQ 33.03 against
            // the direct route's 99 on identical reads. Format v5 always writes the field, so it is
            // always present and always authoritative.
            if (!fragment_len_set) {
                fragment_len = idx.panel.fragment_len;
            } else if (std::abs(fragment_len - idx.panel.fragment_len) > 1e-9) {
                throw std::runtime_error(
                    "genotype: --fragment-len " + std::to_string(fragment_len) +
                    " but the index was built at " + std::to_string(idx.panel.fragment_len) +
                    "; rebuild the index or drop the flag");
            }
        }
        log.info("index " + index_in + ": " + std::to_string(idx.chain.size()) + " blocks, " +
                 std::to_string(idx.haplotype_names.size()) + " haplotypes, " +
                 std::to_string(idx.panel.node_codes.size()) + " markers");
        const ReadCounts rc = count_reads(read_paths, idx.panel, options.threads);
        log.info("reads: " + std::to_string(rc.reads) + " (" + std::to_string(rc.bases / 1000) +
                 " kb); " + std::to_string(rc.syncmers) + " syncmers, " +
                 std::to_string(100 * rc.matched_syncmers / std::max<std::uint64_t>(1, rc.syncmers)) +
                 "% matched a panel marker");
        DepthRegionStats idx_region_stats;
        const std::vector<BlockDepth> depth =
            estimate_depth(idx.panel, rc, min_anchors, uneven_tolerance, depth_model,
                           depth_quantile, 0, depth_estimator, &idx_region_stats);
        log.info("region anchors: " + std::to_string(idx_region_stats.n_anchor) + "; median " +
                 std::to_string(idx_region_stats.median) + ", mean " + std::to_string(idx_region_stats.mean) +
                 ", trimmed " + std::to_string(idx_region_stats.trimmed_mean) +
                 "; selected_anchor_center " + std::to_string(idx_region_stats.used));
        log_depth_provenance(log, idx.chain, depth);
        // Same options as the direct path. Assembling them twice let the two drift: the indexed path
        // silently ignored the recombination rate, carrier weight, provenance, compositional and robust
        // scoring, the mass window and the scale weight, so --index and --bubble-prefix-in did not mean
        // the same thing.
        GenotypeOptions gopt = make_genotype_options();
        GenotypeSummary gsum;
        const std::vector<BlockCall> calls = genotype_sample(idx.chain, idx.blocks, idx.panel, rc,
                                                             depth, idx.haplotype_names, gopt, &gsum);
        log.info("calls: " + std::to_string(gsum.called) + " PASS, " + std::to_string(gsum.no_calls) +
                 " no-call, " + std::to_string(gsum.off_panel) + " off-panel; mean GQ " +
                 std::to_string(gsum.mean_gq));
        {
            std::size_t n_arr = 0;
            for (const BlockCall& c : calls) if (c.is_array) ++n_arr;
            if (n_arr > 0) {
                log.info(std::to_string(n_arr) + " block(s) are tandem arrays (block_class=array): "
                         "there the called allele pair is the closest panel allele BY CONTENT, and "
                         "copy number is mass_bp +- mass_bp_sd, not called_bp");
            }
        }
        write_read_audit(out_prefix, idx.chain, idx.panel, rc, depth);
        if (!dump_markers.empty()) {
            write_marker_dump(dump_markers, idx.chain, idx.panel, rc, depth);
            log.wrote({dump_markers});
        }
        write_genotypes(out_prefix, idx.chain, idx.blocks, calls, idx.haplotype_names);
        log.wrote({out_prefix + ".reads.depth.tsv", out_prefix + ".genotypes.tsv"});
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

        // Haplotype x block allele assignment. Held internally by the linkage audit and never
        // written out, though it is the panel's basic fact: which allele each haplotype carries at
        // each block. Anything reasoning about TRANSITIONS between blocks -- boundary-spanning
        // markers, mosaic structure, linkage -- needs the pairs (allele at bi, allele at bi+1),
        // which cannot be recovered from per-block dumps.
        if (!alleles_out.empty()) {
            std::ofstream af(alleles_out);
            if (!af) throw std::runtime_error("genotype: cannot write " + alleles_out);
            af << "haplotype\tblock_index\tblock_kind\tallele\n";
            std::size_t cells = 0;
            for (const PathRecord& p : panel_graph.paths) {
                for (std::size_t bi = 0; bi < blocks.size() && bi < chain.size(); ++bi) {
                    const auto it = blocks[bi].allele_of.find(p.name);
                    if (it == blocks[bi].allele_of.end()) continue;   // does not traverse this block
                    af << p.name << '\t' << bi << '\t'
                       << (chain[bi].kind == BlockKind::Bubble     ? "bubble"
                           : chain[bi].kind == BlockKind::Flank    ? "flank"
                                                                   : "backbone")
                       << '\t' << it->second << '\n';
                    ++cells;
                }
            }
            log.info("haplotype-allele matrix: " + std::to_string(cells) + " cells over " +
                     std::to_string(panel_graph.paths.size()) + " haplotypes and " +
                     std::to_string(blocks.size()) + " blocks");
            log.wrote({alleles_out});
        }

        // A block index nobody has is a typo, not a request for nothing. Silently writing no ledger
        // and exiting 0 is how an experiment reports "the filters dropped nothing" when in fact it
        // never looked.
        // Every fate the ledger reports is a filter decision. With the filter off there are no
        // decisions to report, and a ledger of "retained, retained, retained" would read as evidence
        // the filter kept everything rather than that it never ran.
        if (ledger_block >= 0 && !options.require_region_unique) {
            throw std::runtime_error(
                "genotype: --ledger-block cannot be combined with --no-region-unique; the ledger "
                "reports why each marker was kept or dropped, and that flag removes the decision");
        }
        if (ledger_block >= 0 && static_cast<std::size_t>(ledger_block) >= chain.size()) {
            throw std::runtime_error("genotype: --ledger-block " + std::to_string(ledger_block) +
                                     " is out of range; the chain has " +
                                     std::to_string(chain.size()) + " blocks (0.." +
                                     std::to_string(chain.size() - 1) + ")");
        }
        ReadPanel read_panel;
        options.restore_stripped_alleles = restore_stripped;
        options.ledger_block = static_cast<int>(ledger_block);
        // The ledger is a property of the PANEL, so it must not require reads. Without this the only
        // way to get one was to ask for an index as well, which is a 40 MB side effect of a question
        // about markers.
        const bool need_panel =
            !read_paths.empty() || !index_out.empty() || ledger_block >= 0;
        const std::vector<BlockMarkerStats> bstats =
            build_block_marker_panel(chain, blocks, options,
                                     need_panel ? &read_panel : nullptr,
                                     need_panel ? &panel_graph : nullptr,
                                     want_audit_stats);
        // Written before anything reads: the ledger is a property of the panel, not of a sample.
        if (!read_panel.ledger.empty()) {
            const std::string lp = out_prefix + ".block" + std::to_string(ledger_block) + ".ledger.tsv";
            std::ofstream lf(lp);
            if (!lf) throw std::runtime_error("genotype: cannot write " + lp);
            lf << "unit\tallele\tslot\tcode\tmult\tvary_blocks\tvary_where\tactual\texpected\tfate\n";
            std::size_t kept = 0, mb = 0, oe = 0, both = 0, edges = 0;
            for (const MarkerLedgerRow& r : read_panel.ledger) {
                const char* fate = r.fate == MarkerFate::Retained     ? "retained"
                                   : r.fate == MarkerFate::MultiBlock ? "multi_block"
                                   : r.fate == MarkerFate::OverExpected ? "over_expected"
                                                                        : "both";
                switch (r.fate) {
                    case MarkerFate::Retained: ++kept; break;
                    case MarkerFate::MultiBlock: ++mb; break;
                    case MarkerFate::OverExpected: ++oe; break;
                    case MarkerFate::Both: ++both; break;
                }
                if (r.unit == MarkerUnit::Edge) ++edges;
                std::string where;
                for (const std::uint32_t b : r.where) {
                    if (!where.empty()) where += ',';
                    where += std::to_string(b);
                }
                lf << (r.unit == MarkerUnit::Edge ? "edge" : "node") << '\t' << r.allele << '\t'
                   << r.slot << '\t' << r.code << '\t' << r.mult << '\t' << r.vary_blocks << '\t'
                   << (where.empty() ? "." : where) << '\t' << r.actual << '\t' << r.expected << '\t'
                   << fate << '\n';
            }
            log.info("block " + std::to_string(ledger_block) + " marker ledger: " +
                     std::to_string(read_panel.ledger.size()) + " candidate occurrences (" +
                     std::to_string(read_panel.ledger.size() - edges) + " node, " +
                     std::to_string(edges) + " edge), " +
                     std::to_string(kept) + " retained, " + std::to_string(mb) + " multi-block, " +
                     std::to_string(oe) + " over-expected, " + std::to_string(both) + " both");
            log.wrote({lp});
        }
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
            // --no-region-unique skips the pass that fills these, so they stay 0. Printed as a
            // percentage that reads "0% confined" -- the opposite of the truth, and indistinguishable
            // from a panel where nothing IS confined. Say "not measured" instead.
            if (!options.require_region_unique) {
                log.info("confinement by context: not measured (--no-region-unique skips the pass)");
            } else {
                log.info("confinement by context: 1-syncmer " +
                         std::to_string(100 * read_panel.confined_vary_nodes /
                                        std::max<std::size_t>(1, read_panel.vary_nodes)) +
                         "% (" + std::to_string(read_panel.confined_vary_nodes) + "/" +
                         std::to_string(read_panel.vary_nodes) + "), 2-syncmer " +
                         std::to_string(100 * read_panel.confined_vary_edges /
                                        std::max<std::size_t>(1, read_panel.vary_edges)) +
                         "% (" + std::to_string(read_panel.confined_vary_edges) + "/" +
                         std::to_string(read_panel.vary_edges) + ")");
            }
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
            const ReadCounts rc_pre = dump_block >= 0 ? count_reads(read_paths, read_panel, options.threads)
                                                     : ReadCounts{};
            // Per-allele diagnostic for one block: does marker multiplicity scale with allele length?
            // At a tandem array it must, or the emission has no copy-number signal to work with and no
            // amount of reweighting will recover one.
            if (dump_block >= 0 && static_cast<std::size_t>(dump_block) < read_panel.by_block.size()) {
                const std::size_t bi = static_cast<std::size_t>(dump_block);
                const std::string dpath = out_prefix + ".block" + std::to_string(bi) + ".tsv";
                std::ofstream d(dpath);
                // obs_sum lets the reader solve for the depth each allele implies. For a marker at
                // multiplicity m in allele a, the reads carry lambda*(m_truth1 + m_truth2), so
                // obs_sum/sum_multiplicity is about 2*lambda for an allele whose copy number matches
                // the sample, lower for one that is too long and higher for one too short.
                // Per-marker confinement audit for this block: why is each surviving marker still
                // here, and does its panel-wide occurrence match what the block accounts for?
                {
                    const std::string apath = out_prefix + ".block" + std::to_string(bi) + ".conf.tsv";
                    std::ofstream cf(apath);
                    cf << "allele\tslot\tmult\tvary_blocks\tocc_blocks\tactual\texpected\tobs\n";
                    for (std::size_t ai = 0; ai < read_panel.by_block[bi].size(); ++ai) {
                        for (const auto& [slot, m] : read_panel.by_block[bi][ai].nodes) {
                            cf << ai << '\t' << slot << '\t' << m << '\t'
                               << (slot < read_panel.dbg_vary.size() ? read_panel.dbg_vary[slot] : 0) << '\t'
                               << (slot < read_panel.dbg_occ.size() ? read_panel.dbg_occ[slot] : 0) << '\t'
                               << (slot < read_panel.dbg_actual.size() ? read_panel.dbg_actual[slot] : 0) << '\t'
                               << (slot < read_panel.dbg_expected.size() ? read_panel.dbg_expected[slot] : 0) << '\t'
                               << rc_pre.node[slot] << '\n';
                        }
                    }
                    log.wrote({apath});
                }
                // Allele sequences for the block, so an external experiment uses the SAME sequence the
                // model does. Reconstructing a block interval by node-id range is not equivalent to the
                // interval walk and silently produces different sequences.
                {
                    const std::string fp = out_prefix + ".block" + std::to_string(bi) + ".fa";
                    std::ofstream ff(fp);
                    // An EMPTY allele is a real allele -- the bypass a deletion takes -- and skipping
                    // it made the dump disagree with every other per-block output about how many
                    // alleles the block has. A consumer counting FASTA records then silently drops
                    // the one allele a deletion is about. Emitted with an empty sequence line.
                    for (std::size_t ai = 0; ai < blocks[bi].allele_seq.size(); ++ai) {
                        ff << '>' << ai << ' ' << blocks[bi].allele_seq[ai].size() << '\n'
                           << blocks[bi].allele_seq[ai] << '\n';
                    }
                    log.wrote({fp});
                }
                d << "allele\tn_haplotypes\tbp\tn_markers\tsum_multiplicity\tobs_sum\n";
                for (std::size_t ai = 0; ai < read_panel.by_block[bi].size(); ++ai) {
                    std::uint64_t mult = 0;
                    std::uint64_t obs = 0;
                    for (const auto& [slot, m] : read_panel.by_block[bi][ai].nodes) {
                        mult += m;
                        obs += rc_pre.node[slot];
                    }
                    d << ai << '\t'
                      << (ai < blocks[bi].allele_haplotypes.size() ? blocks[bi].allele_haplotypes[ai] : 0)
                      << '\t' << (ai < blocks[bi].allele_bp.size() ? blocks[bi].allele_bp[ai] : 0)
                      << '\t' << read_panel.by_block[bi][ai].nodes.size() << '\t' << mult
                      << '\t' << obs << '\n';
                }
                log.wrote({dpath});
            }

            // Why does the emission prefer one candidate pair over another? The likelihood is a sum
            // over thousands of markers, so a 100-unit gap can be one marker screaming or ten thousand
            // whispering, and those call for opposite fixes. This splits the difference by multiplicity
            // band and by whether the two pairs actually disagree about that marker.
            if (!explain_pair.empty()) {
                std::size_t ebi = 0, a1 = 0, a2 = 0, b1 = 0, b2 = 0;
                if (std::sscanf(explain_pair.c_str(), "%zu:%zu,%zu:%zu,%zu", &ebi, &a1, &a2, &b1, &b2) != 5) {
                    throw std::runtime_error("genotype: --explain-pair wants block:a1,a2:b1,b2");
                }
                const ReadCounts rce = count_reads(read_paths, read_panel, options.threads);
                const std::vector<BlockDepth> de =
                    estimate_depth(read_panel, rce, min_anchors, uneven_tolerance, DepthModel::Median,
                                   depth_quantile, 0);
                const double lam = ebi < de.size() ? de[ebi].lambda_hap : 0.0;
                const double mub = 0.02 * lam;
                // Both marker units, because H4 asks exactly this: junction markers are the only
                // local evidence of unit ORDER inside a tandem array, and an allele with extra copies
                // creates junctions the sample may lack -- a penalty node counts cannot express. If
                // adjacencies carry the same presence bias as nodes they will vote the same way, and
                // that settles it without rebuilding the emission's edge path.
                const std::string ep = out_prefix + ".explain.tsv";
                std::ofstream e(ep);
                for (int unit = 0; unit < 2; ++unit) {
                    const bool use_edges = unit == 1;
                    auto profile = [&](std::size_t x, std::size_t y) {
                        std::unordered_map<std::uint64_t, std::uint32_t> t;
                        const auto& B0 = read_panel.by_block[ebi];
                        if (x < B0.size()) for (const auto& [s2, m] : (use_edges ? B0[x].edges : B0[x].nodes)) t[s2] += m;
                        if (y < B0.size()) for (const auto& [s2, m] : (use_edges ? B0[y].edges : B0[y].nodes)) t[s2] += m;
                        return t;
                    };
                    const auto A = profile(a1, a2);
                    const auto B = profile(b1, b2);
                    std::set<std::uint64_t> all;
                    for (const auto& [s2, m] : A) { (void)m; all.insert(s2); }
                    for (const auto& [s2, m] : B) { (void)m; all.insert(s2); }
                    auto lp = [&](double o, double mean) {
                        if (mean <= 0.0) mean = 1e-9;
                        return -mean + o * std::log(mean) - std::lgamma(o + 1.0);
                    };
                    auto obs_of = [&](std::uint64_t slot) {
                        return static_cast<double>(use_edges ? rce.edge[slot] : rce.node[slot]);
                    };
                    struct Cls { std::size_t n = 0; double obs = 0, ll = 0; };
                    Cls only_a, only_b, both_differ, both_same;
                    for (const std::uint64_t s2 : all) {
                        const auto ia = A.find(s2); const auto ib = B.find(s2);
                        const double ma = ia == A.end() ? 0.0 : ia->second;
                        const double mb = ib == B.end() ? 0.0 : ib->second;
                        const double o = obs_of(s2);
                        const double d = lp(o, lam * ma + mub) - lp(o, lam * mb + mub);
                        Cls& c2 = (ma > 0 && mb == 0) ? only_a : (mb > 0 && ma == 0) ? only_b
                                : (ma != mb) ? both_differ : both_same;
                        ++c2.n; c2.obs += o; c2.ll += d;
                    }
                    const char* u = use_edges ? "edges" : "nodes";
                    e << u << "\tcarried_by_A_only\t" << only_a.n << '\t' << only_a.obs << '\t' << only_a.ll << '\n';
                    e << u << "\tcarried_by_B_only\t" << only_b.n << '\t' << only_b.obs << '\t' << only_b.ll << '\n';
                    e << u << "\tboth_differing_copies\t" << both_differ.n << '\t' << both_differ.obs << '\t' << both_differ.ll << '\n';
                    e << u << "\tboth_same_copies\t" << both_same.n << '\t' << both_same.obs << '\t' << both_same.ll << '\n';
                    e << u << "\tTOTAL\t" << all.size() << '\t' << (only_a.obs+only_b.obs+both_differ.obs+both_same.obs)
                      << '\t' << (only_a.ll+only_b.ll+both_differ.ll+both_same.ll) << '\n';
                }
                log.wrote({ep});
            }

            // H5, as a measurement before it is a feature. The defect found by --explain-pair is that
            // scoring WHOLE panel alleles rewards whichever allele contains most of the sample's unit
            // variants, regardless of how many copies that costs. The proposed cure is to stop choosing
            // two alleles and instead infer how much of each is present -- a non-negative mixture whose
            // weights sum to 2 (a diploid) rather than a 0/1/2 indicator.
            //
            // Panel alleles are used as the basis rather than decomposed unit variants, because that
            // needs no repeat decomposition and answers the question that matters first: can ANY
            // non-negative combination explain these counts better than the best pair, and does its
            // implied length match the truth? If it cannot, unit-level decomposition will not either.
            if (deconvolve >= 0) {
                const std::size_t dbi = static_cast<std::size_t>(deconvolve);
                const ReadCounts rcd = count_reads(read_paths, read_panel, options.threads);
                const std::vector<BlockDepth> dd =
                    estimate_depth(read_panel, rcd, min_anchors, uneven_tolerance, DepthModel::Median,
                                   depth_quantile, 0);
                const double lam = dbi < dd.size() ? dd[dbi].lambda_hap : 0.0;
                const auto& BB = read_panel.by_block[dbi];
                std::vector<std::uint32_t> uni;
                for (const auto& ms : BB) for (const auto& [s2, m] : ms.nodes) { (void)m; uni.push_back(s2); }
                std::sort(uni.begin(), uni.end());
                uni.erase(std::unique(uni.begin(), uni.end()), uni.end());
                std::unordered_map<std::uint32_t, std::size_t> row;
                for (std::size_t i = 0; i < uni.size(); ++i) row[uni[i]] = i;
                std::vector<double> obs(uni.size(), 0.0);
                for (std::size_t i = 0; i < uni.size(); ++i) obs[i] = static_cast<double>(rcd.node[uni[i]]);
                // Multiplicative update for non-negative least squares, renormalised to sum 2 each
                // step. Cheap, monotone, and it needs no solver dependency.
                std::vector<double> w(BB.size(), 2.0 / static_cast<double>(std::max<std::size_t>(1, BB.size())));
                std::vector<double> pred(uni.size(), 0.0);
                for (int it = 0; it < 300; ++it) {
                    std::fill(pred.begin(), pred.end(), 0.0);
                    for (std::size_t a = 0; a < BB.size(); ++a) {
                        if (w[a] <= 0.0) continue;
                        for (const auto& [s2, m] : BB[a].nodes) pred[row[s2]] += w[a] * lam * m;
                    }
                    std::vector<double> num(BB.size(), 0.0), den(BB.size(), 0.0);
                    for (std::size_t a = 0; a < BB.size(); ++a) {
                        for (const auto& [s2, m] : BB[a].nodes) {
                            const std::size_t r = row[s2];
                            num[a] += lam * m * obs[r];
                            den[a] += lam * m * pred[r];
                        }
                    }
                    double sum = 0.0;
                    for (std::size_t a = 0; a < BB.size(); ++a) {
                        if (den[a] > 0.0) w[a] *= num[a] / den[a];
                        sum += w[a];
                    }
                    if (sum > 0.0) for (double& x : w) x *= 2.0 / sum;
                }
                double bp_mix = 0.0;
                std::vector<std::pair<double, std::size_t>> top;
                for (std::size_t a = 0; a < BB.size(); ++a) {
                    bp_mix += w[a] * static_cast<double>(a < blocks[dbi].allele_bp.size() ? blocks[dbi].allele_bp[a] : 0);
                    if (w[a] > 0.01) top.emplace_back(w[a], a);
                }
                std::sort(top.rbegin(), top.rend());
                const std::string dp = out_prefix + ".deconv.tsv";
                std::ofstream d2(dp);
                d2 << "#block\t" << dbi << "\tn_alleles\t" << BB.size() << "\tlambda\t" << lam
                   << "\timplied_bp\t" << bp_mix << "\tn_alleles_with_weight_over_0.01\t" << top.size() << '\n';
                d2 << "allele\tweight\tallele_bp\n";
                for (const auto& [ww, a] : top) {
                    d2 << a << '\t' << ww << '\t'
                       << (a < blocks[dbi].allele_bp.size() ? blocks[dbi].allele_bp[a] : 0) << '\n';
                }
                log.wrote({dp});
            }

            // cosigt's score, on our vectors. cosigt compares the sample's per-node coverage against
            // the summed coverage of every candidate haplotype pair by cosine similarity; here the same
            // score runs over one block's marker counts, which is the cheapest way to find out whether
            // the SCORE is what matters before deciding whether the alignment-derived COVERAGE is.
            //
            // The reason to expect a difference: our likelihood predicts a marker the candidate lacks
            // at the error background, so real counts there cost tens of log units and veto the
            // candidate outright. Cosine has no such term -- a marker the candidate lacks contributes
            // nothing to the dot product and only enters the sample's magnitude -- so it should not
            // carry the bias toward longer arrays that --explain-pair exposed.
            if (cosine_block >= 0) {
                const std::size_t cbi = static_cast<std::size_t>(cosine_block);
                const ReadCounts rcc = count_reads(read_paths, read_panel, options.threads);
                const auto& BB = read_panel.by_block[cbi];
                const std::size_t na = BB.size();
                std::vector<std::unordered_map<std::uint32_t, double>> vec(na);
                std::vector<double> self_dot(na, 0.0), obs_dot(na, 0.0);
                double obs_norm2 = 0.0;
                std::unordered_map<std::uint32_t, double> obs;
                for (const auto& ms : BB) {
                    for (const auto& [s2, m] : ms.nodes) {
                        if (obs.find(s2) == obs.end()) {
                            const double o = static_cast<double>(rcc.node[s2]);
                            obs[s2] = o;
                            obs_norm2 += o * o;
                        }
                    }
                }
                for (std::size_t a = 0; a < na; ++a) {
                    for (const auto& [s2, m] : BB[a].nodes) {
                        const double v = static_cast<double>(m);
                        vec[a][s2] = v;
                        self_dot[a] += v * v;
                        obs_dot[a] += v * obs[s2];
                    }
                }
                const double obs_norm = std::sqrt(std::max(1e-12, obs_norm2));
                auto cos_of = [&](std::size_t a, std::size_t b) {
                    double cross = 0.0;
                    const auto& S = vec[a].size() <= vec[b].size() ? vec[a] : vec[b];
                    const auto& L = vec[a].size() <= vec[b].size() ? vec[b] : vec[a];
                    for (const auto& [s2, v] : S) {
                        const auto it = L.find(s2);
                        if (it != L.end()) cross += v * it->second;
                    }
                    const double n2 = self_dot[a] + self_dot[b] + 2.0 * cross;
                    if (n2 <= 0.0) return 0.0;
                    return (obs_dot[a] + obs_dot[b]) / (std::sqrt(n2) * obs_norm);
                };
                std::vector<std::pair<double, std::pair<std::size_t, std::size_t>>> all;
                all.reserve(na * (na + 1) / 2);
                for (std::size_t a = 0; a < na; ++a) {
                    for (std::size_t b = a; b < na; ++b) all.emplace_back(cos_of(a, b), std::make_pair(a, b));
                }
                std::sort(all.rbegin(), all.rend());
                const std::string cp = out_prefix + ".cosine.tsv";
                std::ofstream cf(cp);
                cf << "rank\tallele1\tallele2\tcosine\tbp1\tbp2\ttotal_bp\n";
                auto bpv = [&](std::size_t a) {
                    return a < blocks[cbi].allele_bp.size() ? blocks[cbi].allele_bp[a] : 0;
                };
                for (std::size_t r = 0; r < std::min<std::size_t>(20, all.size()); ++r) {
                    cf << (r + 1) << '\t' << all[r].second.first << '\t' << all[r].second.second << '\t'
                       << all[r].first << '\t' << bpv(all[r].second.first) << '\t'
                       << bpv(all[r].second.second) << '\t'
                       << (bpv(all[r].second.first) + bpv(all[r].second.second)) << '\n';
                }
                log.wrote({cp});
            }

            // Gate A for the coverage evidence path: build the panel's per-node traversal vectors,
            // project the reads onto nodes by alignment, and run the four invariants. Deliberately
            // separate from the marker path so the two can be compared on the same sample rather than
            // replacing each other before there is evidence to justify it.
            if (node_coverage) {
                const NodeIndex nidx = build_node_index(panel_graph);
                const PanelCoverage pcov = build_panel_coverage(panel_graph, nidx);
                std::size_t short_nodes = 0;
                for (const std::uint32_t l : nidx.length) if (l < options.kmer_size) ++short_nodes;
                log.info("node coverage: " + std::to_string(nidx.size()) + " nodes (" +
                         std::to_string(100 * short_nodes / std::max<std::size_t>(1, nidx.size())) +
                         "% shorter than k=" + std::to_string(options.kmer_size) + ", so unreachable by "
                         "any k-mer), " + std::to_string(pcov.path_names.size()) + " panel paths");
                CoverageOptions copt;
                copt.threads = options.threads;
                copt.preset = coverage_preset;
                const SampleCoverage scov = inject_reads(panel_graph, nidx, pcov, read_paths, copt);
                log.info("node coverage: " + std::to_string(scov.reads) + " reads, " +
                         std::to_string(100 * scov.aligned / std::max<std::uint64_t>(1, scov.reads)) +
                         "% placed, " + std::to_string(scov.placements) + " placements (" +
                         std::to_string(scov.placements / std::max<std::uint64_t>(1, scov.aligned)) +
                         " per placed read)");
                std::vector<std::string> tnames;
                if (!truth_haplotypes.empty()) {
                    const auto comma = truth_haplotypes.find(',');
                    if (comma != std::string::npos) {
                        tnames.push_back(truth_haplotypes.substr(0, comma));
                        tnames.push_back(truth_haplotypes.substr(comma + 1));
                    }
                }
                const CoverageAudit ca = audit_coverage(panel_graph, nidx, pcov, scov, tnames, copt);
                log.info("GATE A  probe " + std::to_string(ca.probe_exact) + "/" +
                         std::to_string(ca.probe_total) + " intervals project onto exactly the right "
                         "nodes | mass ratio " + std::to_string(ca.mass_ratio) +
                         " (1.0 = no bases lost or double counted) | self correlation " +
                         std::to_string(ca.self_pearson) + ", slope " + std::to_string(ca.self_slope) +
                         " (slope = per-haplotype read depth)");
                // Phase B at one block: score every allele pair on the coverage vectors by the
                // three candidate scores at once, so likelihood, cosine and Pearson are compared on
                // identical evidence rather than across runs.
                if (coverage_block >= 0 && static_cast<std::size_t>(coverage_block) < chain.size()) {
                    const std::size_t cb = static_cast<std::size_t>(coverage_block);
                    const auto avec = block_allele_node_vectors(panel_graph, path_indexes, bubbles,
                                                                chain[cb], blocks[cb], nidx);
                    // Region-wide depth per haplotype copy: total coverage over the panel's typical
                    // total traversal, halved for the diploid. Sample-independent apart from the reads
                    // themselves, and it needs no anchors -- every node counts.
                    // Depth from ANCHOR nodes: those every panel path traverses exactly once, so the
                    // sample must carry two copies whatever its genotype and their coverage measures
                    // depth alone. Dividing total coverage by the panel's median total traversal was
                    // tried first and is biased -- it assumes the sample's total length equals the
                    // panel median, and at a CNV locus that is exactly the quantity in question. It
                    // read 16.14 where the true depth was 15.0, and a 7% depth error becomes a 7%
                    // copy-number error: 17 kb on a 252 kb array.
                    std::vector<double> anchor_cov;
                    for (std::size_t n = 0; n < nidx.size(); ++n) {
                        bool invariant = !pcov.by_path.empty();
                        for (const auto& v : pcov.by_path) if (v[n] != 1) { invariant = false; break; }
                        if (invariant) anchor_cov.push_back(scov.node[n]);
                    }
                    double lam = 1.0;
                    if (anchor_cov.size() >= 20) {
                        std::sort(anchor_cov.begin(), anchor_cov.end());
                        lam = anchor_cov[anchor_cov.size() / 2] / 2.0;
                    } else {
                        std::vector<double> tot;
                        for (const auto& v : pcov.by_path) {
                            double t = 0.0;
                            for (const std::uint32_t m : v) t += m;
                            tot.push_back(t);
                        }
                        std::sort(tot.begin(), tot.end());
                        double cov_sum = 0.0;
                        for (const double c : scov.node) cov_sum += c;
                        if (!tot.empty() && tot[tot.size() / 2] > 0.0) lam = cov_sum / (2.0 * tot[tot.size() / 2]);
                    }
                    log.info("coverage depth: lambda " + std::to_string(lam) + " from " +
                             std::to_string(anchor_cov.size()) + " invariant nodes");
                    auto sc = score_block_by_coverage(avec, blocks[cb].allele_bp, scov, nidx, lam, 0.0);
                    const std::string cp = out_prefix + ".covscore.tsv";
                    std::ofstream cf(cp);
                    cf << "#block\t" << cb << "\tlambda\t" << lam << "\tn_alleles\t" << avec.size() << "\n";
                    cf << "score\trank\tallele1\tallele2\ttotal_bp\tvalue\n";
                    auto dump = [&](const char* nm, auto key) {
                        auto v = sc;
                        std::sort(v.begin(), v.end(), [&](const CoverageScore& x, const CoverageScore& y) {
                            return key(x) > key(y);
                        });
                        for (std::size_t r = 0; r < std::min<std::size_t>(8, v.size()); ++r) {
                            cf << nm << '\t' << (r + 1) << '\t' << v[r].allele1 << '\t' << v[r].allele2
                               << '\t' << v[r].bp << '\t' << key(v[r]) << '\n';
                        }
                    };
                    dump("loglik", [](const CoverageScore& x) { return x.loglik; });
                    dump("cosine", [](const CoverageScore& x) { return x.cosine; });
                    dump("pearson", [](const CoverageScore& x) { return x.pearson; });
                    log.wrote({cp});
                }
                write_node_coverage(out_prefix, nidx, pcov, scov);
                log.wrote({out_prefix + ".nodecov.sample.tsv", out_prefix + ".nodecov.panel.tsv"});
            }

            // THE HMM'S DECLARED STATE UNIVERSE, built here because the hybrid transaction must be
            // planned BEFORE marker counting: its exclusion set is an input to counting.
            std::vector<std::string> hap_names;
            hap_names.reserve(panel_graph.paths.size());
            for (const PathRecord& p : panel_graph.paths) hap_names.push_back(p.name);

            // ---- ORIENTATION PROBE --------------------------------------------------------
            // Every candidate's CHAIN-ORIENTED sequence and block spans. An exact reverse-complement
            // duplicate must be indistinguishable here from its original: same spans, same bytes.
            // Comparing raw walk coordinates instead is what refused every C4 edge.
            if (!hybrid_orientation_probe.empty()) {
                const FrameCoverage ocov = assess_frame_coverage(graph, blocks, hap_names);
                std::ofstream op(hybrid_orientation_probe);
                if (!op) throw std::runtime_error("genotype: cannot write " +
                                                  hybrid_orientation_probe);
                op << "candidate\treverse_frame\tchain_seq_md5\tchain_len\tblock_spans\n";
                for (std::size_t h = 0; h < ocov.frames.size(); ++h) {
                    const CandidateFrame& F = ocov.frames[h];
                    const std::string cs = chain_oriented_sequence(F);
                    std::string spans;
                    for (std::size_t b = 0; b < blocks.size(); ++b) {
                        std::size_t lo = 0, hi = 0;
                        if (!chain_oriented_block_span(F, static_cast<std::uint32_t>(b), lo, hi)) {
                            continue;
                        }
                        spans += (spans.empty() ? "" : ",") + std::to_string(b) + ":" +
                                 std::to_string(lo) + "-" + std::to_string(hi);
                    }
                    op << ocov.framed_names[h] << '\t' << (F.reverse_frame ? 1 : 0) << '\t'
                       << md5_hex(cs) << '\t' << cs.size() << '\t' << spans << '\n';
                }
                op.flush();
                log.wrote({hybrid_orientation_probe});
            }

            // ---- LINKAGE GEOMETRY PROBE ---------------------------------------------------
            // Geometry depends only on the frames, the block alleles and the flank width -- NOT on
            // any fragment -- so every candidate edge can be tested in seconds instead of behind a
            // ten-minute ownership pass. Written because the first C4 run refused all ten edges and
            // the status recorded only "not-computed": the reason string was dropped, so the
            // failure could not be diagnosed without guessing.
            // ---- THREE-VARIABLE GEOMETRY FEASIBILITY ------------------------------------------
            // A single window below hi-1 kills the whole option, so it is checked first and
            // cheaply, over EVERY haploid combination rather than from the arithmetic
            // 721 - min|A4| -- that shortcut assumes the three minimising alleles co-occur, which
            // is exactly the kind of assumption that has had to be retracted here before.
            if (!hybrid_triple_probe.empty()) {
                if (hybrid_triple_alias_used) {
                    log.info("hybrid: --hybrid-triple-probe is deprecated; it is "
                             "--hybrid-arity-probe, which takes any number of consecutive blocks");
                }
                std::vector<std::string> parts;
                std::string cur;
                for (char c : hybrid_triple_probe) {
                    if (c == ',') { parts.push_back(cur); cur.clear(); } else cur.push_back(c);
                }
                parts.push_back(cur);
                if (parts.size() < 4) {
                    throw std::runtime_error(
                        "genotype: --hybrid-triple-probe expects <b1>,<b2>,...,<path>");
                }
                // GENERALISED to any number of consecutive variable blocks. The three-block form
                // was hardcoded, and its cancellation gate SAMPLED eight alleles per block -- which
                // inverted the conclusion, because the single deciding window used an allele
                // outside the sample. Complete enumeration is now structural rather than optional.
                std::vector<std::size_t> bl;
                for (std::size_t q = 0; q + 1 < parts.size(); ++q) bl.push_back(std::stoul(parts[q]));
                const std::size_t ba = bl.front(), bb = bl.size() > 1 ? bl[1] : bl.front(),
                                  bc = bl.back();
                // THE PATH IS THE LAST FIELD, not the fourth. With three blocks those coincide;
                // with four the fourth field is a BLOCK, and the probe silently wrote a file named
                // "5" and reported nothing.
                std::ofstream tp(parts.back());
                if (!tp) throw std::runtime_error("genotype: cannot write " + parts.back());
                if (ba >= blocks.size() || bb >= blocks.size() || bc >= blocks.size()) {
                    throw std::runtime_error("genotype: --hybrid-triple-probe block out of range");
                }
                const auto& A = blocks[ba].allele_seq;
                const auto& Bv = blocks[bb].allele_seq;
                const auto& Cv = blocks[bc].allele_seq;
                // The two intervening contexts, taken from the pairwise geometries so this probe
                // measures the SAME sequence the pairwise factors do rather than a second opinion.
                const FrameCoverage tcov = assess_frame_coverage(graph, blocks, hap_names);
                std::vector<std::vector<std::string>> tball(blocks.size());
                std::vector<char> tvar(blocks.size(), 0);
                for (std::size_t q = 0; q < blocks.size(); ++q) {
                    tball[q] = blocks[q].allele_seq;
                    tvar[q] = blocks[q].n_alleles > 1 ? 1 : 0;
                }
                const InsertPrior tip = make_insert_prior(fragment_len, hyb_params.fragment_sd,
                                                          hyb_params.discordant_rate,
                                                          hyb_params.insert_sigmas, 1);
                const std::size_t tflank = static_cast<std::size_t>(tip.hi);
                const std::size_t affine_from =
                    static_cast<std::size_t>(std::max<long>(0, tip.hi - 1));
                tp << "field\tvalue\n";
                {
                    std::string bs, as;
                    std::size_t prod = 1;
                    for (std::size_t q = 0; q < bl.size(); ++q) {
                        bs += (q ? "," : "") + std::to_string(bl[q]);
                        as += (q ? "x" : "") + std::to_string(blocks[bl[q]].n_alleles);
                        prod *= blocks[bl[q]].n_alleles;
                    }
                    tp << "blocks\t" << bs << '\n';
                    tp << "alleles\t" << as << '\n';
                    tp << "haploid_windows\t" << prod << '\n';
                }
                // CONSECUTIVE geometries only. The three-block form used g1 = (b0,b1) and
                // g2 = (b1,b2), which for four blocks asks for a geometry between blocks that are
                // not adjacent -- and it refused with "intervening context differs" while reporting
                // the wrong block list.
                bool chain_ok = true;
                std::string chain_refusal;
                for (std::size_t q = 0; q + 1 < bl.size(); ++q) {
                    const LinkageGeometry gq = build_linkage_geometry(
                        tcov.frames, tball, tvar, static_cast<std::uint32_t>(bl[q]),
                        static_cast<std::uint32_t>(bl[q + 1]), tflank, tip);
                    tp << "pair_" << bl[q] << "_" << bl[q + 1] << "_ok\t" << (gq.ok ? 1 : 0)
                       << "\t" << (gq.refusal.empty() ? "-" : gq.refusal) << '\n';
                    if (!gq.ok) { chain_ok = false; chain_refusal = gq.refusal; }
                }
                if (!chain_ok) {
                    tp << "REFUSED\ta consecutive pairwise geometry is unavailable: "
                       << chain_refusal << '\n';
                } else if (bl.size() > 3) {
                    // ---- GENERAL ARITY: complete enumeration over the whole haploid product ----
                    // Contexts come from the consecutive pairwise geometries, so this measures the
                    // same sequence the pairwise factors do.
                    std::vector<std::size_t> ctx;
                    std::size_t lfl = 0, rfl = 0;
                    for (std::size_t q = 0; q + 1 < bl.size(); ++q) {
                        const LinkageGeometry gq = build_linkage_geometry(
                            tcov.frames, tball, tvar, static_cast<std::uint32_t>(bl[q]),
                            static_cast<std::uint32_t>(bl[q + 1]), tflank, tip);
                        ctx.push_back(gq.context.size());
                        if (q == 0) lfl = gq.lflank.size();
                        if (q + 2 == bl.size()) rfl = gq.rflank.size();
                    }
                    {
                        std::size_t fixed = lfl + rfl;
                        for (std::size_t c : ctx) fixed += c;
                        std::vector<std::vector<std::size_t>> lens;
                        std::size_t total = 1;
                        for (std::size_t b : bl) {
                            std::vector<std::size_t> v;
                            for (const std::string& x : blocks[b].allele_seq) v.push_back(x.size());
                            total *= v.size();
                            lens.push_back(std::move(v));
                        }
                        tp << "haploid_windows_total\t" << total << '\n';
                        tp << "fixed_context_bp\t" << fixed << '\n';
                        // Odometer over the complete product: no sampling, no recursion depth cap.
                        std::vector<std::size_t> idx(bl.size(), 0);
                        std::size_t minw = SIZE_MAX, maxw = 0, below = 0, visited = 0;
                        double worst_affine_gap = 0.0;
                        bool done = false;
                        while (!done) {
                            std::size_t w = fixed;
                            for (std::size_t q = 0; q < bl.size(); ++q) w += lens[q][idx[q]];
                            ++visited;
                            minw = std::min(minw, w); maxw = std::max(maxw, w);
                            if (static_cast<long>(w) < tip.hi - 1) ++below;
                            const ExposureCheck ec = check_exposure(w, tip);
                            worst_affine_gap = std::max(worst_affine_gap,
                                                        std::abs(ec.exact - ec.affine));
                            for (std::size_t q = 0; ; ++q) {
                                if (q == bl.size()) { done = true; break; }
                                if (++idx[q] < lens[q].size()) break;
                                idx[q] = 0;
                            }
                        }
                        tp << "haploid_windows_visited\t" << visited << '\n';
                        tp << "complete_enumeration\t" << (visited == total ? 1 : 0) << '\n';
                        tp << "affine_from\t" << (tip.hi - 1) << '\n';
                        tp << "min_window\t" << minw << '\n';
                        tp << "max_window\t" << maxw << '\n';
                        tp << "windows_below_affine\t" << below << '\n';
                        tp << "all_windows_affine\t" << (below == 0 ? 1 : 0) << '\n';
                        // THE ALGEBRAIC INVARIANT, gated exhaustively rather than sampled. If exact
                        // equals affine at EVERY window, then E(w) = w + 1 - E[L] identically;
                        // a diploid's total is E(w1) + E(w2) = (w1 + w2) + 2(1 - E[L]); and
                        // w1 + w2 sums each block's two allele lengths, which no re-pairing of the
                        // homologues can change. Cancellation is then exact by algebra, over ALL
                        // relative phase arrangements, without enumerating any of them.
                        // THE INVARIANT IS PROVED, NOT MEASURED. In regime, max(0, n - L + 1) is
                        // n - L + 1 for every L in the support, so exact = sum p(L)(n - L + 1) =
                        // n + 1 - E[L] = affine IDENTICALLY in real arithmetic. Hence
                        // below == 0 alone establishes it; a diploid's E(w1) + E(w2) is then
                        // (w1 + w2) + 2(1 - E[L]), and w1 + w2 sums each block's two allele
                        // lengths, which no re-pairing of the homologues changes. Cancellation is
                        // exact over every relative phase arrangement without enumerating one.
                        //
                        // The residual below is therefore a NUMERICAL cross-check, not the gate:
                        // `exact` accumulates one term per insert length while `affine` is closed
                        // form, so they differ by double rounding. Reported with its relative size
                        // so a real discrepancy could not hide inside it -- and NOT used to relax
                        // the gate, because a threshold chosen to make a result pass is not a gate.
                        tp << "worst_exact_minus_affine_abs\t" << worst_affine_gap << '\n';
                        tp << "worst_exact_minus_affine_rel\t"
                           << (maxw ? worst_affine_gap / static_cast<double>(maxw) : 0.0) << '\n';
                        tp << "all_in_affine_regime\t" << (below == 0 ? 1 : 0) << '\n';
                        tp << "cancellation_follows_by_algebra\t" << (below == 0 ? 1 : 0) << '\n';
                        log.info("hybrid arity probe over " + std::to_string(bl.size()) +
                                 " blocks: " + std::to_string(visited) + " of " +
                                 std::to_string(total) + " windows, " + std::to_string(minw) + "-" +
                                 std::to_string(maxw) + " bp, " + std::to_string(below) +
                                 " below affine, worst exact-affine gap " +
                                 std::to_string(worst_affine_gap));
                    }
                } else {
                    // THE THREE-BLOCK FORM, and the only place these two geometries are built.
                    // Constructing g2 = (second, last) unconditionally meant that for four blocks
                    // it asked for a geometry between non-adjacent blocks -- dead, malformed, and
                    // a standing hazard even while unused.
                    const LinkageGeometry g1 = build_linkage_geometry(
                        tcov.frames, tball, tvar, static_cast<std::uint32_t>(ba),
                        static_cast<std::uint32_t>(bb), tflank, tip);
                    const LinkageGeometry g2 = build_linkage_geometry(
                        tcov.frames, tball, tvar, static_cast<std::uint32_t>(bb),
                        static_cast<std::uint32_t>(bc), tflank, tip);
                    const std::size_t c1 = g1.context.size(), c2 = g2.context.size();
                    tp << "context_ab\t" << c1 << '\n';
                    tp << "context_bc\t" << c2 << '\n';
                    tp << "lflank\t" << g1.lflank.size() << '\n';
                    tp << "rflank\t" << g2.rflank.size() << '\n';
                    std::size_t minw = SIZE_MAX, maxw = 0, below = 0;
                    for (const std::string& x : A)
                    for (const std::string& y : Bv)
                    for (const std::string& z : Cv) {
                        const std::size_t w = g1.lflank.size() + x.size() + c1 + y.size() + c2 +
                                              z.size() + g2.rflank.size();
                        minw = std::min(minw, w); maxw = std::max(maxw, w);
                        if (w < affine_from) ++below;
                    }
                    tp << "insert_hi\t" << tip.hi << '\n';
                    tp << "affine_from\t" << affine_from << '\n';
                    tp << "min_window\t" << minw << '\n';
                    tp << "max_window\t" << maxw << '\n';
                    tp << "windows_below_affine\t" << below << '\n';
                    tp << "all_windows_affine\t" << (below == 0 ? 1 : 0) << '\n';
                    // EXPOSURE ACROSS ALL FOUR ARRANGEMENTS. With three heterozygous blocks,
                    // fixing the A assignment leaves two choices at B and two at C: four relative
                    // phase arrangements, not straight versus crossed. Cancellation must hold over
                    // all of them, and is asserted rather than argued from affineness.
                    const auto expo = [&](const std::string& x, const std::string& y,
                                          const std::string& z) {
                        const std::size_t w = g1.lflank.size() + x.size() + c1 + y.size() + c2 +
                                              z.size() + g2.rflank.size();
                        return check_exposure(w, tip).exact;
                    };
                    // EVERY arrangement, not a sample. Capping the loops at eight alleles per
                    // block would have left the one sub-affine window possibly outside the checked
                    // set -- and that window is the entire question, so a sampled answer about it
                    // is no answer. 16^2 * 10^2 * 8^2 is 1.6M tuples and costs nothing.
                    double asym = 0.0;
                    const std::size_t la = A.size(), lb = Bv.size(), lc = Cv.size();
                    for (std::size_t i1 = 0; i1 < la; ++i1)
                    for (std::size_t i2 = 0; i2 < la; ++i2)
                    for (std::size_t j1 = 0; j1 < lb; ++j1)
                    for (std::size_t j2 = 0; j2 < lb; ++j2)
                    for (std::size_t k1 = 0; k1 < lc; ++k1)
                    for (std::size_t k2 = 0; k2 < lc; ++k2) {
                        const double base = expo(A[i1], Bv[j1], Cv[k1]) +
                                            expo(A[i2], Bv[j2], Cv[k2]);
                        const double alt1 = expo(A[i1], Bv[j2], Cv[k1]) +
                                            expo(A[i2], Bv[j1], Cv[k2]);
                        const double alt2 = expo(A[i1], Bv[j1], Cv[k2]) +
                                            expo(A[i2], Bv[j2], Cv[k1]);
                        const double alt3 = expo(A[i1], Bv[j2], Cv[k2]) +
                                            expo(A[i2], Bv[j1], Cv[k1]);
                        asym = std::max(asym, std::abs(base - alt1));
                        asym = std::max(asym, std::abs(base - alt2));
                        asym = std::max(asym, std::abs(base - alt3));
                    }
                    // THE SAME GATE AS THE GENERAL BRANCH. Reporting a numeric-equality verdict
                    // here while k > 3 uses the proved invariant would be two different answers to
                    // one question, and the numeric one reads as a failure on pure rounding.
                    tp << "max_exposure_asymmetry_all_arrangements\t" << asym << '\n';
                    tp << "all_in_affine_regime\t" << (below == 0 ? 1 : 0) << '\n';
                    tp << "cancellation_follows_by_algebra\t" << (below == 0 ? 1 : 0) << '\n';
                    // WHICH window is short, and by how much. One window below the boundary is a
                    // different situation from many, and the size of the shortfall says whether
                    // the residual asymmetry is a rounding artefact or a real term.
                    std::size_t sa = 0, sb = 0, sc = 0;
                    std::size_t worst_short = 0;
                    for (std::size_t i = 0; i < A.size(); ++i)
                    for (std::size_t j = 0; j < Bv.size(); ++j)
                    for (std::size_t k = 0; k < Cv.size(); ++k) {
                        const std::size_t w = g1.lflank.size() + A[i].size() + c1 + Bv[j].size() +
                                              c2 + Cv[k].size() + g2.rflank.size();
                        if (w < affine_from && (worst_short == 0 || w < worst_short)) {
                            worst_short = w; sa = i; sb = j; sc = k;
                        }
                    }
                    if (worst_short != 0) {
                        tp << "shortest_window_alleles\t" << sa << "," << sb << "," << sc << '\n';
                        tp << "shortest_window_bp\t" << worst_short << '\n';
                        tp << "shortfall_bp\t" << (affine_from - worst_short) << '\n';
                        tp << "allele_lengths\t" << A[sa].size() << "," << Bv[sb].size() << ","
                           << Cv[sc].size() << '\n';
                    }
                    log.info("hybrid triple probe " + std::to_string(ba) + "," +
                             std::to_string(bb) + "," + std::to_string(bc) + ": " +
                             std::to_string(A.size()) + "x" + std::to_string(Bv.size()) + "x" +
                             std::to_string(Cv.size()) + ", windows " + std::to_string(minw) +
                             "-" + std::to_string(maxw) + ", affine from " +
                             std::to_string(affine_from) + ", " + std::to_string(below) +
                             " below, asymmetry " + std::to_string(asym));
                }
            }
            if (!hybrid_geometry_probe.empty()) {
                const FrameCoverage gcov = assess_frame_coverage(graph, blocks, hap_names);
                std::vector<std::vector<std::string>> ball(blocks.size());
                for (std::size_t b = 0; b < blocks.size(); ++b) ball[b] = blocks[b].allele_seq;
                long minl = 1;
                for (const std::string& rp : read_paths) { (void)rp; }
                const InsertPrior gip = make_insert_prior(fragment_len, hyb_params.fragment_sd,
                                                          hyb_params.discordant_rate,
                                                          hyb_params.insert_sigmas,
                                                          std::max<long>(1, minl));
                std::ofstream gp(hybrid_geometry_probe);
                if (!gp) throw std::runtime_error("genotype: cannot write " + hybrid_geometry_probe);
                gp << "block_a\tblock_b\tn_alleles_a\tn_alleles_b\tvariable_a\tvariable_b"
                      "\tflank_bp\tlflank_derived\trflank_derived\tok\treason\n";
                // ADJACENT VARIABLE PAIRS, which is what ownership actually produces -- not
                // adjacent physical pairs. They coincide only where every block is variable (C4);
                // where a fixed backbone separates two bubbles the real edge skips over it, and
                // probing physical neighbours would miss exactly the edge that carries phase.
                std::vector<std::size_t> varb;
                for (std::size_t b = 0; b < blocks.size(); ++b) {
                    if (blocks[b].n_alleles > 1) varb.push_back(b);
                }
                for (std::size_t vi = 1; vi < varb.size(); ++vi) {
                    const std::size_t a = varb[vi - 1], b = varb[vi];
                    // Three widths: the production reach, an intermediate one where the REQUEST
                    // is smaller than the available invariant context (so the cap must bite
                    // exactly), and zero.
                    for (std::size_t flank : {static_cast<std::size_t>(gip.hi),
                                              static_cast<std::size_t>(100),
                                              static_cast<std::size_t>(0)}) {
                        std::vector<char> gvar(blocks.size(), 0);
                        for (std::size_t q = 0; q < blocks.size(); ++q) {
                            gvar[q] = blocks[q].n_alleles > 1 ? 1 : 0;
                        }
                        const LinkageGeometry g2 = build_linkage_geometry(
                            gcov.frames, ball, gvar, static_cast<std::uint32_t>(a),
                            static_cast<std::uint32_t>(b), flank, gip);
                        gp << a << '\t' << b << '\t' << blocks[a].n_alleles << '\t'
                           << blocks[b].n_alleles << '\t' << (blocks[a].n_alleles > 1 ? 1 : 0)
                           << '\t' << (blocks[b].n_alleles > 1 ? 1 : 0) << '\t' << flank << '\t'
                           << g2.lflank_bp << '\t' << g2.rflank_bp << '\t'
                           << (g2.ok ? 1 : 0) << '\t'
                           << (g2.refusal.empty() ? "-" : g2.refusal) << '\n';
                    }
                }
                gp.flush();
                log.wrote({hybrid_geometry_probe});
            }

            // ---- HYBRID TRANSACTION -------------------------------------------------------
            // Planned in full before a single occurrence is subtracted. Every failure path leaves
            // the marker counts untouched, so the legacy call remains exactly what it would be.
            FrameCoverage hyb_cov;
            HybridActivation hyb_act;
            std::unordered_set<std::string> hyb_exclusions;
            // THE HIGHER-ORDER MATERIAL, declared where the genotyper can still see it. The
            // tables must outlive every pointer into them, so they live here rather than in the
            // block that builds them.
            std::vector<IntervalFactorTable> built_tables;
            std::vector<std::vector<std::uint32_t>> built_blocks;
            std::vector<std::string> built_supersede;
            std::vector<HybridHigherFactor> higher_factors;
            std::vector<std::vector<std::uint32_t>> higher_hap_allele;
            HigherOrderPlan higher_plan;
            HigherOrderStats higher_stats;
            std::string higher_refusal;
            bool higher_active = false;
            // Rendered where `owners` is in scope; printed with the status file later.
            std::vector<std::string> unusable_lines;
            std::size_t neutral_pairwise_reported = 0, superseded_reported = 0;
            std::vector<HigherFactorScope> higher_scopes;      // scopes CREDITED to the ledger
            std::vector<HigherFactorScope> planned_only_stored; // what the planner proposed
            std::size_t planned_only_reported = 0;
            std::vector<std::size_t> plan_wide_consumed;       // Wide fragments each factor takes
            long min_len_recorded = 0, ip_lo_recorded = 0, ip_hi_recorded = 0;
            double ip_residual_lo = 0.0, ip_residual_hi = 0.0;
            std::size_t hyb_fragments_loaded = 0, hyb_edges_considered = 0;
            if (hybrid_call) {
                hyb_cov = assess_frame_coverage(graph, blocks, hap_names);
                if (!hyb_cov.all_states_usable) {
                    hyb_act.refusal = "candidate-frame-coverage";
                } else {
                    const std::vector<Fragment> hf = load_fragments(read_paths);
                    hyb_fragments_loaded = hf.size();
                    std::vector<char> block_variable(blocks.size(), 0);
                    for (std::size_t b = 0; b < blocks.size(); ++b) {
                        block_variable[b] = blocks[b].n_alleles > 1 ? 1 : 0;
                    }
                    // THE INSERT SUPPORT, and the assumption hidden in it.
                    //
                    // r1 + r2 as the minimum says "an insert shorter than the two mates is not a
                    // state at all", which is only true if the mates cannot OVERLAP. They can:
                    // a 350 +- 50 library with 150 bp mates puts 15.9% of fragments under 300 bp,
                    // and those pairs are physically real -- the two reads simply share sequence.
                    // Their emissions still multiply correctly conditional on the template,
                    // because the sequencing errors are separate observations.
                    //
                    // --hybrid-overlap-pairs uses max(r1, r2) instead, which is the true floor:
                    // an insert cannot be shorter than its longest mate. It is a FLAG rather than
                    // the default because widening the support changes the prior normalisation,
                    // the exposure, every fragment mass, ownership, and therefore which fragments
                    // F1 and F2 own -- so the ledger and the real factors have to be rebuilt after
                    // it, not merely re-read.
                    hyb_params.allow_overlapping_pairs = !hybrid_no_overlap_pairs;
                    const long min_len =
                        fragment_insert_floor(hf, hyb_params.allow_overlapping_pairs);
                    hyb_params.fragment_len = fragment_len;   // the library geometry in force
                    // ALWAYS COMPUTED, even when a supplied lambda is used, so the ratio between
                    // the two is reportable and the arms are comparable.
                    std::vector<std::size_t> panel_lengths;
                    panel_lengths.reserve(hyb_cov.frames.size());
                    for (const CandidateFrame& F : hyb_cov.frames) {
                        panel_lengths.push_back(F.seq.size());
                    }
                    hyb_params.n_fragments = hf.size();
                    hyb_params.lambda_estimated = estimate_fragment_lambda(
                        hf.size(), panel_lengths, &hyb_params.median_panel_length);
                    if (hyb_params.lambda_source ==
                        HybridLinkageParameters::LambdaSource::Estimated) {
                        hyb_params.lambda = hyb_params.lambda_estimated;
                    }
                    min_len_recorded = min_len;
                    const InsertPrior ip = make_insert_prior(hyb_params.fragment_len,
                                                             hyb_params.fragment_sd,
                                                             hyb_params.discordant_rate,
                                                             hyb_params.insert_sigmas, min_len);
                    ip_lo_recorded = ip.lo; ip_hi_recorded = ip.hi;
                    ip_residual_lo = ip.log_residual_below;
                    ip_residual_hi = ip.log_residual_above;
                    const double lep = std::log(hyb_params.error_rate / 3.0);
                    const double l1m = std::log1p(-hyb_params.error_rate);
                    // OWNERSHIP over the SAME frames the preflight reported: the object is passed,
                    // not recomputed, so report and decision cannot describe different runs.
                    // THE PIECE INDEX, built once over the panel and reused for every fragment.
                    // Without it bounded_mate_placements scans each candidate linearly: on C4 that
                    // is 131 candidates x 226 kb x 4 mates x 23953 fragments, and the run did not
                    // finish in 55 minutes. The index is the same pigeonhole structure the interval
                    // scorer uses; omitting it here was a wiring omission, not a design choice.
                    std::vector<PieceIndex> pidx;
                    {
                        std::size_t piece = 0;
                        for (const Fragment& F : hf) {
                            if (F.r1.empty()) continue;
                            piece = F.r1.size() /
                                    (mate_band_edits(hyb_params.max_divergence, F.r1.size()) + 1);
                            break;
                        }
                        if (piece >= 12) {
                            pidx.resize(hyb_cov.frames.size());
                            run_parallel(hyb_cov.frames.size(), options.threads,
                                         [&](std::size_t h) {
                                             pidx[h] = build_piece_index(hyb_cov.frames[h].seq,
                                                                         piece);
                                         });
                            log.info("hybrid: piece index over " +
                                     std::to_string(hyb_cov.frames.size()) + " candidates at " +
                                     std::to_string(piece) + " bp");
                        }
                    }
                    std::vector<FragmentOwner> owners(hf.size());
                    std::atomic<std::size_t> own_done{0};
                    const auto own_t0 = std::chrono::steady_clock::now();
                    run_parallel(hf.size(), options.threads, [&](std::size_t fi) {
                        owners[fi] = assign_fragment_owner(hf[fi], hyb_cov.frames, block_variable,
                                                           ip, hyb_params.max_divergence, lep, l1m,
                                                           1e-6, pidx.empty() ? nullptr : &pidx);
                        // PROGRESS. A silent multi-minute pass is exactly what made the first C4
                        // attempt uninterpretable until it was killed.
                        const std::size_t n = ++own_done;
                        if (n % 5000 == 0) {
                            log.info("hybrid: ownership " + std::to_string(n) + " / " +
                                     std::to_string(hf.size()) + " fragments");
                        }
                    });
                    log.info("hybrid: ownership over " + std::to_string(hf.size()) +
                             " fragments in " +
                             std::to_string(std::chrono::duration<double>(
                                 std::chrono::steady_clock::now() - own_t0).count()) + " s");
                    std::map<std::pair<std::uint32_t, std::uint32_t>,
                             std::vector<std::size_t>> by_edge;
                    for (std::size_t fi = 0; fi < hf.size(); ++fi) {
                        if (owners[fi].kind == OwnerKind::Linkage) {
                            by_edge[{owners[fi].block_lo, owners[fi].block_hi}].push_back(fi);
                        }
                    }
                    hyb_edges_considered = by_edge.size();
                    // ---- RUN A REAL INTERVAL FACTOR ------------------------------------------
                    // Its evidence is the pairwise owners inside the block span plus the Wide
                    // fragments whose whole variable scope lies inside it -- the same rule the
                    // ledger used, so the counts must reconcile with it exactly.
                    for (std::size_t frun = 0; frun < hybrid_factor_runs.size(); ++frun) {
                        const std::string& hybrid_factor_run = hybrid_factor_runs[frun];
                        const std::string hybrid_factor_supersede =
                            frun < hybrid_factor_supersedes.size()
                                ? hybrid_factor_supersedes[frun] : std::string();
                        std::vector<std::string> fp;
                        std::string fcu;
                        for (char c : hybrid_factor_run) {
                            if (c == ',') { fp.push_back(fcu); fcu.clear(); } else fcu.push_back(c);
                        }
                        fp.push_back(fcu);
                        std::vector<std::uint32_t> fbl;
                        for (std::size_t q = 0; q + 1 < fp.size(); ++q)
                            fbl.push_back(static_cast<std::uint32_t>(std::stoul(fp[q])));
                        std::set<std::uint32_t> fbs(fbl.begin(), fbl.end());
                        std::vector<std::vector<std::string>> fball(blocks.size());
                        for (std::size_t q = 0; q < blocks.size(); ++q)
                            fball[q] = blocks[q].allele_seq;
                        const IntervalGeometry FG = build_interval_geometry(
                            hyb_cov.frames, fball, block_variable, fbl,
                            static_cast<std::size_t>(ip.hi), ip);
                        std::ofstream fo(fp.back());
                        if (!fo) throw std::runtime_error("genotype: cannot write " + fp.back());
                        fo << "field\tvalue\n";
                        std::string fbstr;
                        for (std::size_t q = 0; q < fbl.size(); ++q)
                            fbstr += (q ? "," : "") + std::to_string(fbl[q]);
                        fo << "blocks\t" << fbstr << '\n';
                        fo << "geometry_ok\t" << (FG.ok ? 1 : 0) << '\n';
                        if (!FG.ok) {
                            fo << "refusal\t" << FG.refusal << '\n';
                        } else {
                            // THE EVIDENCE SET. SCOPE IS NOT OWNERSHIP: a factor depends on the
                            // blocks in its span, but it consumes only the fragments assigned to
                            // it. Taking every pairwise owner inside the span is wrong -- for
                            // {4,5,6} that sweeps in edge 4-5's owners, which belong to F1, and
                            // edge 5-6's, which stay with the pairwise factor that is RETAINED.
                            // A factor consumes a pairwise edge's owners only when it SUPERSEDES
                            // that edge, which is a stated fact about the model and not something
                            // to infer from block membership.
                            std::set<std::pair<std::uint32_t, std::uint32_t>> sup_edges;
                            {
                                std::string tok;
                                std::vector<std::string> toks;
                                for (char c : hybrid_factor_supersede) {
                                    if (c == ',') { toks.push_back(tok); tok.clear(); }
                                    else tok.push_back(c);
                                }
                                if (!tok.empty()) toks.push_back(tok);
                                for (const std::string& t : toks) {
                                    const std::size_t dash = t.find('-');
                                    if (dash == std::string::npos) {
                                        throw std::runtime_error(
                                            "genotype: --hybrid-factor-supersede expects a-b,c-d");
                                    }
                                    sup_edges.insert({
                                        static_cast<std::uint32_t>(std::stoul(t.substr(0, dash))),
                                        static_cast<std::uint32_t>(std::stoul(t.substr(dash + 1)))});
                                }
                            }
                            std::vector<std::size_t> ev;
                            std::size_t n_from_superseded = 0;
                            for (const auto& kv : by_edge) {
                                if (!sup_edges.count(kv.first)) continue;
                                if (!fbs.count(kv.first.first) || !fbs.count(kv.first.second)) {
                                    throw std::runtime_error(
                                        "genotype: superseded edge " +
                                        std::to_string(kv.first.first) + "-" +
                                        std::to_string(kv.first.second) +
                                        " is not inside the factor's own scope");
                                }
                                ev.insert(ev.end(), kv.second.begin(), kv.second.end());
                                n_from_superseded += kv.second.size();
                            }
                            std::size_t n_wide_in = 0;
                            for (std::size_t fi = 0; fi < hf.size(); ++fi) {
                                if (owners[fi].kind != OwnerKind::Wide) continue;
                                bool inside = !owners[fi].var_scope.empty();
                                for (std::uint32_t b : owners[fi].var_scope)
                                    if (!fbs.count(b)) { inside = false; break; }
                                if (inside) { ev.push_back(fi); ++n_wide_in; }
                            }
                            std::sort(ev.begin(), ev.end());
                            const std::size_t before_u = ev.size();
                            ev.erase(std::unique(ev.begin(), ev.end()), ev.end());
                            fo << "haploid_cells\t" << FG.cells() << '\n';
                            fo << "min_window\t" << FG.min_window << '\n';
                            fo << "exposure_affine\t" << (FG.exposure_affine ? 1 : 0) << '\n';
                            fo << "superseded_edges\t"
                               << (hybrid_factor_supersede.empty() ? "-" : hybrid_factor_supersede)
                               << '\n';
                            fo << "evidence_from_superseded_edges\t" << n_from_superseded << '\n';
                            fo << "evidence_fragments\t" << ev.size() << '\n';
                            fo << "evidence_unique\t" << (before_u == ev.size() ? 1 : 0) << '\n';
                            fo << "evidence_wide\t" << n_wide_in << '\n';
                            const std::size_t piece0 = ev.empty() ? 16
                                : hf[ev.front()].r1.size() /
                                  (mate_band_edits(hyb_params.max_divergence,
                                                   hf[ev.front()].r1.size()) + 1);
                            const auto t_ix = std::chrono::steady_clock::now();
                            const IntervalSeedIndex FX = build_interval_seed_index(FG, piece0);
                            const double ix_s = std::chrono::duration<double>(
                                std::chrono::steady_clock::now() - t_ix).count();
                            fo << "piece\t" << piece0 << '\n';
                            fo << "index_ok\t" << (FX.ok ? 1 : 0) << '\n';
                            fo << "index_complete\t" << (FX.complete ? 1 : 0) << '\n';
                            fo << "shortest_allele\t" << FX.shortest_allele << '\n';
                            fo << "longest_invariant\t" << FX.longest_context << '\n';
                            fo << "index_build_seconds\t" << ix_s << '\n';
                            fo << "index_inside_keys\t" << FX.inside.size() << '\n';
                            fo << "index_boundary_keys\t" << FX.boundary.size() << '\n';
                            std::size_t contributed = 0, refused = 0, fin_tot = 0;
                            std::uint64_t seedh = 0, symst = 0, joins = 0, tup = 0, ver = 0,
                                          acc = 0, frst = 0, writes = 0, pbd = 0, pad = 0;
                            std::size_t cells_differ = 0, org_differ = 0, sig_differ2 = 0,
                                        con_differ2 = 0, oracle_cells = 0, oracle_finite = 0,
                                        oracle_absent = 0;
                            double worst_mass = 0.0;
                            std::ofstream mf;
                            if (hybrid_factor_oracle) {
                                mf.open(fp.back() + ".manifest.tsv");
                                if (mf) mf << "fragment\tcell\tallele_tuple\n";
                            }
                            const auto t_run = std::chrono::steady_clock::now();
                            // PROGRESS, per fragment. The first attempt at this run went 40
                            // minutes with no output: unable to say whether it was a tenth done or
                            // nearly finished, and therefore unable to justify either waiting or
                            // killing it. A long silent loop is not a measurement.
                            std::size_t done_n = 0;
                            double score_s = 0.0, oracle_s = 0.0;
                            std::vector<std::vector<std::string>> frag_sigs;
                            std::vector<IntervalEmission> ems_kept;
                            const auto t_prog = std::chrono::steady_clock::now();
                            for (std::size_t fi : ev) {
                                const std::size_t len = hf[fi].bases();
                                const std::size_t be = static_cast<std::size_t>(
                                    hyb_params.bg_divergence * static_cast<double>(len));
                                const double bgf = static_cast<double>(be) * lep +
                                                   static_cast<double>(len - be) * l1m;
                                const auto t_sc0 = std::chrono::steady_clock::now();
                                const IntervalEmission E = interval_emission(
                                    hf[fi], FG, ip, hyb_params.max_divergence, lep, l1m, bgf,
                                    &FX, nullptr, true, hybrid_factor_oracle);
                                score_s += std::chrono::duration<double>(
                                    std::chrono::steady_clock::now() - t_sc0).count();
                                if (E.work_refused || !E.ok) { ++refused; continue; }
                                ++contributed;
                                frag_sigs.push_back(E.cell_signature);
                                ems_kept.push_back(E);
                                seedh += E.seed_hits; symst += E.symbolic_states;
                                pbd += E.placements_before_dedup; pad += E.placements_after_dedup;
                                joins += E.joined_pairs; tup += E.tuple_expansions;
                                ver += E.full_read_verifications; acc += E.accepted_placements;
                                frst += E.verified_fr_states; fin_tot += E.finite_cells;
                                for (std::uint32_t cs : E.cell_states) writes += cs;
                                const auto t_or0 = std::chrono::steady_clock::now();
                                if (hybrid_factor_oracle) {
                                    // THE FROZEN MANIFEST, selected WITHOUT looking at any result:
                                    // a deterministic hash of the fragment name strides through the
                                    // cell space, plus the length extremes, the all-non-reference
                                    // tuple, and -- for a four-block factor -- a cross grid over
                                    // the OUTER blocks, which is where the 122x fan-out lives. The
                                    // selection is written to disk before a single comparison runs
                                    // and is never adjusted afterwards.
                                    std::uint64_t hsh = 1469598103934665603ull;
                                    for (char ch : hf[fi].name) {
                                        hsh ^= static_cast<std::uint64_t>(
                                            static_cast<unsigned char>(ch));
                                        hsh *= 1099511628211ull;
                                    }
                                    std::set<std::size_t> pick;
                                    const std::size_t NCELL = FG.cells();
                                    for (std::size_t q = 0; q < 32; ++q)
                                        pick.insert((hsh + q * 2654435761ull) % NCELL);
                                    // Length extremes and the all-non-reference tuple.
                                    std::vector<std::uint32_t> tmin(FG.alleles.size(), 0),
                                                               tmax(FG.alleles.size(), 0),
                                                               tnr(FG.alleles.size(), 0);
                                    for (std::size_t j = 0; j < FG.alleles.size(); ++j) {
                                        std::size_t lo2 = SIZE_MAX, hi2 = 0;
                                        for (std::uint32_t a = 0; a < FG.alleles[j].size(); ++a) {
                                            const std::size_t L = FG.alleles[j][a].size();
                                            if (L < lo2) { lo2 = L; tmin[j] = a; }
                                            if (L > hi2) { hi2 = L; tmax[j] = a; }
                                        }
                                        tnr[j] = FG.alleles[j].size() > 1 ? 1 : 0;
                                    }
                                    pick.insert(FG.cell_index(tmin));
                                    pick.insert(FG.cell_index(tmax));
                                    pick.insert(FG.cell_index(tnr));
                                    // The outer-block cross grid, on ONE fragment, capped.
                                    if (FG.alleles.size() == 4 && done_n == 0) {
                                        std::vector<std::uint32_t> t(4, 0);
                                        for (std::uint32_t a0 = 0;
                                             a0 < FG.alleles[0].size() && a0 < 25; ++a0)
                                        for (std::uint32_t a3 = 0; a3 < FG.alleles[3].size(); ++a3) {
                                            t[0] = a0; t[1] = 0; t[2] = 0; t[3] = a3;
                                            pick.insert(FG.cell_index(t));
                                        }
                                    }
                                    for (std::size_t c : pick) {
                                        if (mf) {
                                            std::vector<std::uint32_t> tc;
                                            FG.cell_choice(c, tc);
                                            mf << hf[fi].name << '\t' << c << '\t';
                                            for (std::size_t j = 0; j < tc.size(); ++j)
                                                mf << (j ? "," : "") << tc[j];
                                            mf << '\n';
                                        }
                                        ++oracle_cells;
                                        const IntervalOracleCell Oc = interval_oracle_cell(
                                            hf[fi], FG, ip, hyb_params.max_divergence, lep, l1m, c);
                                        const bool a1 = E.mass[c] != -INFINITY;
                                        const bool b1 = Oc.mass != -INFINITY;
                                        if (a1) ++oracle_finite; else ++oracle_absent;
                                        if (a1 != b1) { ++cells_differ; continue; }
                                        if (E.cell_origin[c] != Oc.origin) ++org_differ;
                                        if (E.cell_signature[c] != Oc.signature) ++sig_differ2;
                                        if (E.cell_contrib[c].size() != Oc.contrib.size()) ++con_differ2;
                                        else for (std::size_t q = 0; q < Oc.contrib.size(); ++q)
                                            if (E.cell_contrib[c][q] != Oc.contrib[q]) {
                                                ++con_differ2; break;
                                            }
                                        if (a1) worst_mass = std::max(worst_mass,
                                                                      std::abs(E.mass[c] - Oc.mass));
                                    }
                                }
                                oracle_s += std::chrono::duration<double>(
                                    std::chrono::steady_clock::now() - t_or0).count();
                                if (++done_n % 2 == 0 || done_n == ev.size()) {
                                    const double el = std::chrono::duration<double>(
                                        std::chrono::steady_clock::now() - t_prog).count();
                                    log.info("hybrid factor " + fbstr + ": " +
                                             std::to_string(done_n) + "/" +
                                             std::to_string(ev.size()) + " fragments, " +
                                             std::to_string(el) + " s elapsed, projected " +
                                             std::to_string(el / static_cast<double>(done_n) *
                                                            static_cast<double>(ev.size())) + " s");
                                }
                            }
                            // ---- K-DIMENSIONAL GROUPING, on this factor's real signatures ----
                            // The content-class space over RAW alleles is the product of
                            // C(n_j + 1, 2) and is far too large to enumerate; over signature
                            // classes it is the product of C(R_j + 1, 2). Whether the factor is
                            // constructible at all turns on that ratio.
                            const auto t_grp = std::chrono::steady_clock::now();
                            const IntervalGrouping GR =
                                build_interval_grouping(FG, frag_sigs);
                            const double grp_s = std::chrono::duration<double>(
                                std::chrono::steady_clock::now() - t_grp).count();
                            fo << "grouping_ok\t" << (GR.ok ? 1 : 0) << '\n';
                            if (GR.ok) {
                                std::string cpb;
                                for (std::size_t j = 0; j < GR.classes_per_block.size(); ++j)
                                    cpb += (j ? "x" : "") +
                                           std::to_string(GR.classes_per_block[j]);
                                std::string apb;
                                for (std::size_t j = 0; j < FG.alleles.size(); ++j)
                                    apb += (j ? "x" : "") + std::to_string(FG.alleles[j].size());
                                fo << "alleles_per_block\t" << apb << '\n';
                                fo << "signature_classes_per_block\t" << cpb << '\n';
                                fo << "distinct_cell_signatures\t" << GR.distinct_cell_signatures
                                   << '\n';
                                fo << "content_classes_raw\t" << GR.content_classes_raw << '\n';
                                fo << "content_classes_grouped\t" << GR.content_classes_grouped
                                   << '\n';
                                fo << "factor_lookup_configurations\t" << GR.factor_lookup_configurations
                                   << '\n';
                                fo << "classes_m_le_1_neutral\t" << GR.classes_m_le_1 << '\n';
                                fo << "classes_carrying_a_value\t"
                                   << (GR.content_classes_grouped - GR.classes_m_le_1) << '\n';
                                fo << "ordered_in_m_le_1\t" << GR.ordered_in_m_le_1 << '\n';
                                fo << "stored_canonical_values\t" << GR.stored_canonical_values
                                   << '\n';
                                fo << "stored_ordered_values\t" << GR.stored_ordered_values
                                   << '\n';
                                fo << "predicted_bytes_canonical\t" << GR.predicted_bytes << '\n';
                                fo << "cells_checked\t" << GR.cells_checked << '\n';
                                fo << "cells_disagreeing_with_representative\t"
                                   << GR.cells_disagreeing_with_representative << '\n';
                                fo << "joint_equality_verified\t"
                                   << (GR.joint_equality_verified ? 1 : 0) << '\n';
                                // THE ACTUAL ALLELE -> CLASS VECTORS, per block. Equal class
                                // COUNTS on a shared block do not prove identical partitions, and
                                // a refinement claim between two factors must be checked on the
                                // mappings themselves.
                                for (std::size_t j = 0; j < GR.allele_class.size(); ++j) {
                                    fo << "allele_class_block_" << FG.blocks[j] << '\t';
                                    for (std::size_t a = 0; a < GR.allele_class[j].size(); ++a)
                                        fo << (a ? "," : "") << GR.allele_class[j][a];
                                    fo << '\n';
                                }
                                fo << "content_class_reduction\t"
                                   << (GR.content_classes_grouped
                                        ? static_cast<double>(GR.content_classes_raw) /
                                          static_cast<double>(GR.content_classes_grouped) : 0.0)
                                   << '\n';
                                fo << "grouping_seconds\t" << grp_s << '\n';
                                log.info("hybrid factor " + fbstr + " grouping: " + apb +
                                         " alleles -> " + cpb + " signature classes; content "
                                         "classes " + std::to_string(GR.content_classes_raw) +
                                         " -> " + std::to_string(GR.content_classes_grouped));
                            }
                            // ---- BUILD THE ACTUAL FACTOR -------------------------------------
                            if (GR.ok) {
                                const auto t_fac = std::chrono::steady_clock::now();
                                IntervalFactorTable FT = build_interval_factor(
                                    FG, GR, ems_kept, hyb_params.lambda,
                                    std::log1p(-hyb_params.outlier_mix),
                                    std::log(hyb_params.outlier_mix));
                                const double fac_s = std::chrono::duration<double>(
                                    std::chrono::steady_clock::now() - t_fac).count();
                                fo << "factor_ok\t" << (FT.ok ? 1 : 0) << '\n';
                                fo << "factor_refusal\t"
                                   << (FT.refusal.empty() ? "-" : FT.refusal) << '\n';
                                fo << "factor_table_empty_on_refusal\t"
                                   << ((FT.ok || (FT.phase_value.empty() &&
                                                  FT.class_offset.empty())) ? 1 : 0) << '\n';
                                fo << "factor_classes_stored\t" << FT.classes_stored << '\n';
                                fo << "factor_phase_values\t" << FT.phase_values_stored << '\n';
                                fo << "factor_classes_neutral\t" << FT.classes_neutral << '\n';
                                fo << "factor_swap_failures\t" << FT.swap_partner_failures << '\n';
                                fo << "factor_max_log_psi\t" << FT.max_log_psi << '\n';
                                fo << "factor_worst_bound_slack\t" << FT.worst_bound_slack << '\n';
                                fo << "factor_payload_bytes\t" << FT.canonical_payload_bytes
                                   << '\n';
                                fo << "factor_payload_predicted\t" << FT.predicted_payload_bytes
                                   << '\n';
                                fo << "factor_total_bytes\t" << FT.total_factor_bytes << '\n';
                                fo << "factor_total_predicted\t" << FT.predicted_total_bytes
                                   << '\n';
                                fo << "factor_bytes_match_prediction\t"
                                   << ((FT.canonical_payload_bytes == FT.predicted_payload_bytes &&
                                        FT.total_factor_bytes == FT.predicted_total_bytes) ? 1 : 0)
                                   << '\n';
                                fo << "factor_build_seconds\t" << fac_s << '\n';
                                log.info("hybrid factor " + fbstr + " table: " +
                                         (FT.ok ? "usable" : "REFUSED " + FT.refusal) + ", " +
                                         std::to_string(FT.classes_stored) + " classes, " +
                                         std::to_string(FT.phase_values_stored) + " phase values, " +
                                         std::to_string(FT.total_factor_bytes) + " bytes, " +
                                         std::to_string(fac_s) + " s");
                                // RETAINED so the recurrence can be run on the real tables rather
                                // than on a reconstruction of them.
                                if (FT.ok) {
                                    built_tables.push_back(std::move(FT));
                                    built_blocks.push_back(fbl);
                                    built_supersede.push_back(hybrid_factor_supersede);
                                }
                            }
                            const double run_s = std::chrono::duration<double>(
                                std::chrono::steady_clock::now() - t_run).count();
                            fo << "fragments_contributed\t" << contributed << '\n';
                            fo << "fragments_refused\t" << refused << '\n';
                            fo << "all_contributed\t"
                               << ((contributed == ev.size() && refused == 0) ? 1 : 0) << '\n';
                            fo << "seed_hits\t" << seedh << '\n';
                            fo << "symbolic_states\t" << symst << '\n';
                            fo << "joined_pairs\t" << joins << '\n';
                            fo << "tuple_expansions\t" << tup << '\n';
                            fo << "full_read_verifications\t" << ver << '\n';
                            fo << "accepted_placements\t" << acc << '\n';
                            fo << "placements_before_dedup\t" << pbd << '\n';
                            fo << "placements_after_dedup\t" << pad << '\n';
                            fo << "verified_fr_states\t" << frst << '\n';
                            fo << "cell_state_writes\t" << writes << '\n';
                            fo << "finite_fragment_cells\t" << fin_tot << '\n';
                            // THREE DIFFERENT QUANTITIES, and conflating them would make a
                            // production figure look thirty times worse than it is. Symbolic
                            // scoring is what production would pay; the manifest, oracle and
                            // comparison exist only to certify it and would never run in a call.
                            fo << "symbolic_scoring_seconds\t" << score_s << '\n';
                            fo << "manifest_oracle_compare_seconds\t" << oracle_s << '\n';
                            fo << "factor_total_seconds\t" << run_s << '\n';
                            if (hybrid_factor_oracle) {
                                fo << "oracle_cells_checked\t" << oracle_cells << '\n';
                                fo << "oracle_cells_finite\t" << oracle_finite << '\n';
                                fo << "oracle_cells_absent\t" << oracle_absent << '\n';
                                fo << "oracle_cells_differ\t" << cells_differ << '\n';
                                fo << "oracle_origin_differ\t" << org_differ << '\n';
                                fo << "oracle_signature_differ\t" << sig_differ2 << '\n';
                                fo << "oracle_contribution_differ\t" << con_differ2 << '\n';
                                fo << "oracle_worst_mass\t" << worst_mass << '\n';
                            }
                            struct rusage ru {};
                            double pk = 0.0;
                            if (getrusage(RUSAGE_SELF, &ru) == 0) {
#ifdef __APPLE__
                                pk = static_cast<double>(ru.ru_maxrss) / 1048576.0;
#else
                                pk = static_cast<double>(ru.ru_maxrss) / 1024.0;
#endif
                            }
                            // WHOLE-PROCESS RSS, not this factor's. What the difference from the
                            // pre-factor baseline is attributable to has NOT been measured, and
                            // naming a cause without phase-level measurement would be a guess.
                            fo << "process_peak_rss_mb\t" << pk << '\n';
                            fo << "process_peak_rss_attribution\tUNMEASURED\n";
                            log.info("hybrid factor " + fbstr + ": " + std::to_string(contributed) +
                                     "/" + std::to_string(ev.size()) + " fragments, " +
                                     std::to_string(FG.cells()) + " cells, index complete=" +
                                     (FX.complete ? "yes" : "NO") + ", " + std::to_string(run_s) +
                                     " s" + (hybrid_factor_oracle
                                        ? ", oracle differ " + std::to_string(cells_differ) + "/" +
                                          std::to_string(org_differ) : ""));
                        }
                    }
                    // ---- THE RECURRENCE, ON THE REAL TABLES -----------------------------------
                    // Not a synthetic loop: the factors just built, the real panel, the real block
                    // chain. EMISSIONS ARE UNIFORM, and deliberately so -- nothing here prunes on
                    // emission value, so every reachable state is enumerated whatever the emissions
                    // are, and uniform ones give the same work as real ones. A degenerate emission
                    // could only remove states, never add them, so this is the upper end.
                    if (!hybrid_higher_bench.empty()) {
                        std::ofstream bo(hybrid_higher_bench);
                        if (!bo) throw std::runtime_error("genotype: cannot write " +
                                                          hybrid_higher_bench);
                        bo << "field\tvalue\n";
                        bo << "factors_built\t" << built_tables.size() << '\n';
                        std::size_t lo_b = blocks.size(), hi_b = 0;
                        for (const auto& bl : built_blocks)
                            for (std::uint32_t q : bl) {
                                lo_b = std::min<std::size_t>(lo_b, q);
                                hi_b = std::max<std::size_t>(hi_b, q);
                            }
                        if (built_tables.empty()) {
                            bo << "status\tNO_FACTORS\n";
                        } else {
                            // THE SLICE is exactly the blocks the factors span. Extending it would
                            // add plain Li-Stephens steps that measure the existing kernel, not
                            // this recurrence.
                            bool continue_after_plan = false;
                            const std::size_t NB = hi_b - lo_b + 1;
                            std::vector<std::size_t> keep_h;
                            for (std::size_t h = 0; h < hap_names.size(); ++h) {
                                bool full = true;
                                for (std::size_t q = lo_b; q <= hi_b && full; ++q)
                                    if (!blocks[q].allele_of.count(hap_names[h])) full = false;
                                if (full) keep_h.push_back(h);
                            }
                            const std::size_t NH = keep_h.size();
                            bo << "slice_blocks\t" << lo_b << "-" << hi_b << '\n';
                            bo << "slice_n_blocks\t" << NB << '\n';
                            bo << "panel_haplotypes\t" << hap_names.size() << '\n';
                            bo << "haplotypes_spanning_slice\t" << NH << '\n';
                            HybridChain HC;
                            // ANY r IN (0,1) GIVES THE SAME WORK. The enumeration depends on r
                            // only through whether the stay and switch components are active, not
                            // through their weights, so the reachable state set -- and therefore
                            // every count below -- is identical for any interior r. Only r = 0 or
                            // r = 1 would change it, by deleting a component.
                            HC.n_hap = NH; HC.n_blocks = NB; HC.recomb = 0.5;
                            HC.edges.assign(NB, HybridEdge{});
                            HC.log_emission.assign(NB, std::vector<double>(NH * NH, 0.0));
                            HC.hap_allele.assign(NH, std::vector<std::uint32_t>(NB, 0));
                            for (std::size_t i = 0; i < NH; ++i)
                                for (std::size_t q = 0; q < NB; ++q)
                                    HC.hap_allele[i][q] = static_cast<std::uint32_t>(
                                        blocks[lo_b + q].allele_of.at(hap_names[keep_h[i]]));
                            for (std::size_t f = 0; f < built_tables.size(); ++f) {
                                HybridHigherFactor HF;
                                for (std::uint32_t q : built_blocks[f])
                                    HF.blocks.push_back(static_cast<std::uint32_t>(q - lo_b));
                                HF.table = &built_tables[f];
                                HC.higher.push_back(HF);
                                bo << "factor_" << f << "_blocks\t";
                                for (std::size_t j = 0; j < built_blocks[f].size(); ++j)
                                    bo << (j ? "," : "") << built_blocks[f][j];
                                bo << '\n';
                                bo << "factor_" << f << "_classes\t"
                                   << built_tables[f].classes_stored << '\n';
                            }
                            // ---- THE PLAN, ALWAYS, BEFORE ANYTHING IS ALLOCATED ----------
                            const auto t_pl = std::chrono::steady_clock::now();
                            const HigherOrderPlan PL =
                                plan_higher_order(HC, hybrid_plan_threads);
                            const double pl_s = std::chrono::duration<double>(
                                std::chrono::steady_clock::now() - t_pl).count();
                            bo << "plan_ok\t" << (PL.ok ? 1 : 0) << '\n';
                            bo << "plan_refusal\t" << (PL.refusal.empty() ? "-" : PL.refusal)
                               << '\n';
                            bo << "plan_seconds\t" << pl_s << '\n';
                            if (PL.ok) {
                                for (std::size_t q = 0; q < PL.haploid_states.size(); ++q) {
                                    bo << "plan_block_" << (lo_b + q) << "_refinement\t"
                                       << PL.refinement_classes[q] << '\n';
                                    bo << "plan_block_" << (lo_b + q) << "_haploid_states\t"
                                       << PL.haploid_states[q] << '\n';
                                    bo << "plan_block_" << (lo_b + q) << "_message_entries\t"
                                       << PL.message_entries[q] << '\n';
                                }
                                bo << "plan_peak_message_entries\t"
                                   << PL.peak_message_entries << '\n';
                                bo << "plan_total_message_entries\t"
                                   << PL.total_message_entries << '\n';
                                bo << "plan_forward_updates\t" << PL.forward_updates << '\n';
                                bo << "plan_adjoint_updates\t" << PL.adjoint_updates << '\n';
                                bo << "plan_total_updates\t"
                                   << (PL.forward_updates + PL.adjoint_updates) << '\n';
                                bo << "plan_dense_equivalent_updates\t"
                                   << PL.dense_equivalent_updates << '\n';
                                bo << "plan_multiplier_reconstructions\t"
                                   << PL.multiplier_reconstructions << '\n';
                                bo << "plan_factor_lookups\t" << PL.factor_lookups << '\n';
                                bo << "plan_history_ops\t" << PL.history_ops << '\n';
                                bo << "plan_grouping_ops\t" << PL.grouping_ops << '\n';
                                bo << "plan_payload_bytes\t" << PL.payload_bytes << '\n';
                                bo << "plan_container_bytes\t" << PL.container_bytes << '\n';
                                bo << "plan_temporary_bytes\t" << PL.temporary_bytes << '\n';
                                bo << "plan_threads\t" << PL.threads << '\n';
                                bo << "plan_per_thread_bytes\t" << PL.per_thread_bytes << '\n';
                                bo << "plan_total_bytes\t" << PL.total_bytes << '\n';
                                bo << "plan_total_gb\t"
                                   << (static_cast<double>(PL.total_bytes) / 1073741824.0) << '\n';
                            }
                            if (hybrid_plan_only) {
                                bo << "status\tPLAN_ONLY\n";
                                log.info("higher-order PLAN ONLY on blocks " +
                                         std::to_string(lo_b) + "-" + std::to_string(hi_b) + ": " +
                                         (PL.ok ? (std::to_string(PL.forward_updates) +
                                                   " forward updates, " +
                                                   std::to_string(PL.total_bytes / 1048576) +
                                                   " MB")
                                                : ("REFUSED " + PL.refusal)));
                                continue_after_plan = true;
                            }
                            if (!continue_after_plan) {
                            const auto t_ho = std::chrono::steady_clock::now();
                            HigherOrderStats HS;
                            const HybridPosterior HP = hybrid_higher_order(
                                HC, hybrid_max_message_entries, &HS);
                            const double ho_s = std::chrono::duration<double>(
                                std::chrono::steady_clock::now() - t_ho).count();
                            bo << "status\t" << (HP.ok ? "COMPLETE" : "REFUSED") << '\n';
                            bo << "refusal\t" << (HS.refusal.empty() ? "-" : HS.refusal) << '\n';
                            bo << "recomb\t" << HC.recomb << '\n';
                            bo << "recomb_note\tany interior r gives the same enumeration\n";
                            bo << "forward_updates\t" << HS.forward_updates << '\n';
                            bo << "adjoint_updates\t" << HS.adjoint_updates << '\n';
                            bo << "total_updates\t"
                               << (HS.forward_updates + HS.adjoint_updates) << '\n';
                            bo << "forward_equals_adjoint\t"
                               << (HS.forward_updates == HS.adjoint_updates ? 1 : 0) << '\n';
                            // SEPARATE CATEGORIES. The tape-free adjoint buys memory with
                            // recomputation, so these are not interchangeable with the updates.
                            bo << "multiplier_reconstructions\t"
                               << HS.multiplier_reconstructions << '\n';
                            bo << "factor_lookups\t" << HS.factor_lookups << '\n';
                            bo << "history_ops\t" << HS.history_ops << '\n';
                            bo << "grouping_ops\t" << HS.grouping_ops << '\n';
                            bo << "peak_message_entries\t" << HS.peak_message_entries << '\n';
                            bo << "total_message_entries\t" << HS.total_message_entries << '\n';
                            bo << "predicted_payload_bytes\t"
                               << HS.predicted_payload_bytes << '\n';
                            bo << "predicted_payload_note\tmessages and adjoints only; excludes "
                                  "hash index, allocator overhead, factor tables, temporaries and "
                                  "any reduction buffers\n";
                            bo << "forward_seconds\t" << HS.forward_seconds << '\n';
                            bo << "adjoint_seconds\t" << HS.adjoint_seconds << '\n';
                            bo << "total_seconds\t" << ho_s << '\n';
                            struct rusage bru {};
                            double bpk = 0.0;
                            if (getrusage(RUSAGE_SELF, &bru) == 0) {
#ifdef __APPLE__
                                bpk = static_cast<double>(bru.ru_maxrss) / 1048576.0;
#else
                                bpk = static_cast<double>(bru.ru_maxrss) / 1024.0;
#endif
                            }
                            // WHOLE-PROCESS peak, which includes everything built before this
                            // point -- the factor tables above all. Not the recurrence's own.
                            bo << "process_peak_rss_mb\t" << bpk << '\n';
                            bo << "process_peak_rss_attribution\tUNMEASURED\n";
                            log.info("higher-order recurrence on blocks " + std::to_string(lo_b) +
                                     "-" + std::to_string(hi_b) + ": " +
                                     (HP.ok ? "COMPLETE" : "REFUSED " + HS.refusal) + ", " +
                                     std::to_string(HS.forward_updates) + " forward + " +
                                     std::to_string(HS.adjoint_updates) + " adjoint updates, " +
                                     std::to_string(ho_s) + " s");
                            }
                        }
                    }

                    // ---- INTERVAL GEOMETRY CROSS-CHECK ----------------------------------------
                    // The generic k-block geometry, reported through the type the factor will
                    // actually use. Its window statistics must reproduce the arity probe's, which
                    // computed them by an independent path -- a new type agreeing with an existing
                    // measurement is worth more than the same code reporting itself twice.
                    if (!hybrid_interval_probe.empty()) {
                        std::vector<std::string> ipp;
                        std::string icu;
                        for (char c : hybrid_interval_probe) {
                            if (c == ',') { ipp.push_back(icu); icu.clear(); } else icu.push_back(c);
                        }
                        ipp.push_back(icu);
                        std::vector<std::uint32_t> ib;
                        for (std::size_t q = 0; q + 1 < ipp.size(); ++q)
                            ib.push_back(static_cast<std::uint32_t>(std::stoul(ipp[q])));
                        std::vector<std::vector<std::string>> iball(blocks.size());
                        for (std::size_t q = 0; q < blocks.size(); ++q)
                            iball[q] = blocks[q].allele_seq;
                        const IntervalGeometry IG = build_interval_geometry(
                            hyb_cov.frames, iball, block_variable, ib,
                            static_cast<std::size_t>(ip.hi), ip);
                        std::ofstream io(ipp.back());
                        if (!io) throw std::runtime_error("genotype: cannot write " + ipp.back());
                        io << "field\tvalue\n";
                        std::string ibs;
                        for (std::size_t q = 0; q < ib.size(); ++q)
                            ibs += (q ? "," : "") + std::to_string(ib[q]);
                        io << "blocks\t" << ibs << '\n';
                        io << "ok\t" << (IG.ok ? 1 : 0) << '\n';
                        if (!IG.ok) {
                            io << "refusal\t" << IG.refusal << '\n';
                        } else {
                            io << "arity\t" << IG.arity() << '\n';
                            io << "cells\t" << IG.cells() << '\n';
                            std::string as;
                            for (std::size_t q = 0; q < IG.alleles.size(); ++q)
                                as += (q ? "x" : "") + std::to_string(IG.alleles[q].size());
                            io << "alleles\t" << as << '\n';
                            io << "lflank\t" << IG.lflank.size() << '\n';
                            io << "rflank\t" << IG.rflank.size() << '\n';
                            std::string cs;
                            for (std::size_t q = 0; q < IG.contexts.size(); ++q)
                                cs += (q ? "," : "") + std::to_string(IG.contexts[q].size());
                            io << "contexts\t" << cs << '\n';
                            io << "min_window\t" << IG.min_window << '\n';
                            io << "max_window\t" << IG.max_window << '\n';
                            io << "windows_below_affine\t" << IG.windows_below_affine << '\n';
                            io << "exposure_affine\t" << (IG.exposure_affine ? 1 : 0) << '\n';
                            // ROUND TRIP: cell_index and cell_choice must invert each other over
                            // the whole product, or every tensor index derived from them is wrong.
                            std::vector<std::uint32_t> ch;
                            std::size_t bad = 0;
                            for (std::size_t k2 = 0; k2 < IG.cells(); ++k2) {
                                IG.cell_choice(k2, ch);
                                if (IG.cell_index(ch) != k2) ++bad;
                            }
                            io << "index_roundtrip_failures\t" << bad << '\n';
                        }
                        log.info("hybrid interval geometry " + ibs + ": cells " +
                                 std::to_string(IG.ok ? IG.cells() : 0) + ", windows " +
                                 std::to_string(IG.min_window) + "-" +
                                 std::to_string(IG.max_window) + ", affine=" +
                                 (IG.exposure_affine ? "yes" : "no"));
                    }
                    // ---- CONTEXT DEPENDENCE ---------------------------------------------------
                    // Does a neighbouring block change any fragment's STRUCTURAL SIGNATURE, or is
                    // it only sequence that makes the window long enough for exposure to be affine?
                    //
                    // The question is exact and does not need a four-dimensional search. A state
                    // lying wholly to the right of block ctx shifts uniformly when ctx's allele
                    // changes: its start moves, its mismatch counts, insert length and multiplicity
                    // do not, so its signature is invariant by construction. Only a state TOUCHING
                    // ctx can differ -- and a state touches ctx only if one of the four mate
                    // variants has an accepted placement overlapping ctx's span. So enumerate
                    // exactly those: placements wholly inside each ctx allele, and placements
                    // straddling its right boundary, which reach at most one read length past it.
                    //
                    // If there are none, ctx cannot affect any signature, and the factor's
                    // statistical scope excludes it however much sequence it contributes.
                    if (!hybrid_context_dependence.empty()) {
                        std::vector<std::string> cp;
                        std::string cu;
                        for (char c : hybrid_context_dependence) {
                            if (c == ',') { cp.push_back(cu); cu.clear(); } else cu.push_back(c);
                        }
                        cp.push_back(cu);
                        const std::uint32_t ctx_b = static_cast<std::uint32_t>(std::stoul(cp[0]));
                        std::vector<std::uint32_t> fb;
                        for (std::size_t q = 1; q + 1 < cp.size(); ++q) {
                            fb.push_back(static_cast<std::uint32_t>(std::stoul(cp[q])));
                        }
                        std::set<std::uint32_t> fbset(fb.begin(), fb.end());
                        std::ofstream cd(cp.back());
                        if (!cd) throw std::runtime_error("genotype: cannot write " + cp.back());
                        // THE FACTOR'S EVIDENCE: the pairwise owners inside the block set plus the
                        // Wide fragments whose whole variable scope lies inside it.
                        std::vector<std::size_t> fset;
                        for (const auto& kv : by_edge) {
                            if (fbset.count(kv.first.first) && fbset.count(kv.first.second))
                                fset.insert(fset.end(), kv.second.begin(), kv.second.end());
                        }
                        for (std::size_t fi = 0; fi < hf.size(); ++fi) {
                            if (owners[fi].kind != OwnerKind::Wide) continue;
                            bool inside = !owners[fi].var_scope.empty();
                            for (std::uint32_t b : owners[fi].var_scope)
                                if (!fbset.count(b)) { inside = false; break; }
                            if (inside) fset.push_back(fi);
                        }
                        std::sort(fset.begin(), fset.end());
                        fset.erase(std::unique(fset.begin(), fset.end()), fset.end());
                        // The downstream sequence immediately after ctx, from ctx's own pairwise
                        // geometry so it is the same sequence the factor would see.
                        std::vector<std::vector<std::string>> cball(blocks.size());
                        for (std::size_t q = 0; q < blocks.size(); ++q)
                            cball[q] = blocks[q].allele_seq;
                        // THE CONTEXT MAY SIT ON EITHER SIDE. Building (ctx, factor.front())
                        // unconditionally asks for a geometry in the wrong order when ctx is to the
                        // RIGHT, and the probe then refused with "blocks overlap or are out of
                        // order" -- reporting nothing rather than testing block 7 at all.
                        const bool ctx_left = ctx_b < fb.front();
                        const std::uint32_t ga = ctx_left ? ctx_b : fb.back();
                        const std::uint32_t gb = ctx_left ? fb.front() : ctx_b;
                        const LinkageGeometry gctx = build_linkage_geometry(
                            hyb_cov.frames, cball, block_variable, ga, gb,
                            static_cast<std::size_t>(ip.hi), ip);
                        cd << "field\tvalue\n";
                        cd << "context_block\t" << ctx_b << '\n';
                        cd << "factor_blocks\t" << cp[1];
                        for (std::size_t q = 2; q + 1 < cp.size(); ++q) cd << "," << cp[q];
                        cd << '\n';
                        cd << "evidence_fragments\t" << fset.size() << '\n';
                        cd << "context_alleles\t" << blocks[ctx_b].n_alleles << '\n';
                        cd << "downstream_geometry_ok\t" << (gctx.ok ? 1 : 0) << '\n';
                        if (!gctx.ok) {
                            cd << "REFUSED\t" << gctx.refusal << '\n';
                        } else {
                            std::size_t touching = 0, placements = 0;
                            std::vector<std::size_t> touch_frags;
                            const auto& ctx_alleles = blocks[ctx_b].allele_seq;
                            // Which neighbour sequence abuts the context on the side facing the
                            // factor: its own alleles when the context is on the right.
                            const std::vector<std::string>& abut =
                                ctx_left ? gctx.alleles_b : gctx.alleles_a;
                            // One index per context allele, reused across every fragment.
                            for (std::size_t ai = 0; ai < ctx_alleles.size(); ++ai) {
                                const std::string& A2 = ctx_alleles[ai];
                                if (A2.empty()) continue;
                                PieceIndex px;
                                bool have_px = false;
                                for (std::size_t k = 0; k < fset.size(); ++k) {
                                    const Fragment& f = hf[fset[k]];
                                    const std::size_t d1 = mate_band_edits(
                                        hyb_params.max_divergence, f.r1.size());
                                    const std::size_t d2 = mate_band_edits(
                                        hyb_params.max_divergence, f.r2.size());
                                    const std::size_t p1 = f.r1.size() / (d1 + 1);
                                    if (!have_px && p1 >= 12 && A2.size() > 8 * p1) {
                                        px = build_piece_index(A2, p1); have_px = true;
                                    }
                                    const std::string a1 = reverse_complement(f.r1);
                                    const std::string a2s = reverse_complement(f.r2);
                                    const std::string* mv[4] = {&f.r1, &a1, &f.r2, &a2s};
                                    const std::size_t bd[4] = {d1, d1, d2, d2};
                                    std::size_t hits = 0;
                                    for (int m = 0; m < 4; ++m) {
                                        // Wholly inside the context allele.
                                        const auto pl = bounded_mate_placements(
                                            *mv[m], A2, bd[m], nullptr,
                                            have_px ? &px : nullptr);
                                        hits += pl.size();
                                        // Straddling its right boundary: at most one read length
                                        // of the allele's tail plus the same of what follows.
                                        const std::size_t L = mv[m]->size();
                                        if (L > 1) {
                                            const std::string tailA =
                                                A2.size() <= L - 1 ? A2 : A2.substr(A2.size() - (L - 1));
                                            const std::string headA =
                                                A2.size() <= L - 1 ? A2 : A2.substr(0, L - 1);
                                            for (const std::string& nx : abut) {
                                                // Facing side: for a left context the factor is
                                                // downstream, for a right context it is upstream.
                                                const std::string near =
                                                    ctx_left
                                                        ? (nx.size() <= L - 1 ? nx : nx.substr(0, L - 1))
                                                        : (nx.size() <= L - 1 ? nx
                                                                              : nx.substr(nx.size() - (L - 1)));
                                                const std::string bstr =
                                                    ctx_left ? tailA + gctx.context + near
                                                             : near + gctx.context + headA;
                                                const std::size_t ctx_lo =
                                                    ctx_left ? 0 : near.size() + gctx.context.size();
                                                const std::size_t ctx_hi =
                                                    ctx_left ? tailA.size() : bstr.size();
                                                const auto pb = bounded_mate_placements(
                                                    *mv[m], bstr, bd[m], nullptr, nullptr);
                                                for (const MatePlacement& q : pb) {
                                                    const long e = q.start + static_cast<long>(L);
                                                    if (q.start < static_cast<long>(ctx_hi) &&
                                                        e > static_cast<long>(ctx_lo)) ++hits;
                                                }
                                            }
                                        }
                                    }
                                    if (hits > 0) {
                                        placements += hits;
                                        touch_frags.push_back(fset[k]);
                                    }
                                }
                            }
                            std::sort(touch_frags.begin(), touch_frags.end());
                            touch_frags.erase(std::unique(touch_frags.begin(), touch_frags.end()),
                                              touch_frags.end());
                            touching = touch_frags.size();
                            // WHOSE fragments they are. The evidence set here is every pairwise
                            // owner inside the block span plus the Wide fragments scoped to it --
                            // which for {4,5,6} includes edge 4-5's owners, and those belong to a
                            // DIFFERENT factor. Without this breakdown a touch by someone else's
                            // fragment reads as a dependency of this one.
                            std::size_t touch_wide = 0, touch_link = 0;
                            std::map<std::pair<std::uint32_t,std::uint32_t>, std::size_t> touch_edge;
                            for (std::size_t fi : touch_frags) {
                                if (owners[fi].kind == OwnerKind::Wide) ++touch_wide;
                                else if (owners[fi].kind == OwnerKind::Linkage) {
                                    ++touch_link;
                                    touch_edge[{owners[fi].block_lo, owners[fi].block_hi}] += 1;
                                }
                            }
                            std::size_t ev_wide = 0, ev_link = 0;
                            for (std::size_t fi : fset) {
                                if (owners[fi].kind == OwnerKind::Wide) ++ev_wide; else ++ev_link;
                            }
                            cd << "evidence_wide\t" << ev_wide << '\n';
                            cd << "evidence_linkage\t" << ev_link << '\n';
                            cd << "touching_wide\t" << touch_wide << '\n';
                            cd << "touching_linkage\t" << touch_link << '\n';
                            for (const auto& te : touch_edge) {
                                cd << "touching_from_edge_" << te.first.first << "_"
                                   << te.first.second << "\t" << te.second << '\n';
                            }
                            cd << "fragments_with_a_placement_touching_context\t" << touching
                               << '\n';
                            cd << "total_touching_placements\t" << placements << '\n';
                            cd << "context_is_emission_relevant\t" << (touching ? 1 : 0) << '\n';
                            cd << "effective_statistical_scope\t"
                               << (touching ? "includes block " + std::to_string(ctx_b)
                                            : "EXCLUDES block " + std::to_string(ctx_b)) << '\n';
                            log.info("hybrid context dependence: block " + std::to_string(ctx_b) +
                                     " over " + std::to_string(fset.size()) + " evidence "
                                     "fragments -- " + std::to_string(touching) +
                                     " have a placement touching it, " +
                                     std::to_string(placements) + " placements");
                        }
                    }
                    // ---- WIDE INVENTORY -------------------------------------------------------
                    // Wide fragments are unconsumed evidence, and under the locus-wide transaction
                    // ANY unconsumed Wide keeps the model INCOMPLETE. So repairing the exposure
                    // edges is necessary and NOT sufficient: without this inventory the project
                    // could finish that work and find C4 still unreachable for a reason already
                    // visible now. Grouped by exact variable scope, because that is what decides
                    // which factor could ever consume them.
                    if (!hybrid_wide_inventory.empty()) {
                        std::map<std::vector<std::uint32_t>, std::size_t> by_scope;
                        std::size_t n_wide = 0;
                        for (std::size_t fi = 0; fi < hf.size(); ++fi) {
                            if (owners[fi].kind != OwnerKind::Wide) continue;
                            ++n_wide;
                            by_scope[owners[fi].var_scope] += 1;
                        }
                        std::ofstream wi(hybrid_wide_inventory);
                        if (!wi) throw std::runtime_error("genotype: cannot write " +
                                                          hybrid_wide_inventory);
                        wi << "var_scope\tarity\tfragments\tmin_block\tmax_block"
                              "\tminimal_factor_interval\tinterval_blocks\tconsecutive_variable"
                              "\tcovered_by_existing_edge\n";
                        std::size_t need_higher = 0, covered = 0;
                        for (const auto& kv : by_scope) {
                            const auto& sc = kv.first;
                            std::string ss;
                            for (std::size_t q = 0; q < sc.size(); ++q)
                                ss += (q ? "," : "") + std::to_string(sc[q]);
                            const std::uint32_t lo = sc.empty() ? 0 : sc.front();
                            const std::uint32_t hi = sc.empty() ? 0 : sc.back();
                            // The minimal factor interval is the variable run from lo to hi: a
                            // factor cannot depend on a subset of the variables a fragment reads.
                            std::size_t interval_blocks = 0;
                            for (std::size_t b = lo; b <= hi && b < blocks.size(); ++b)
                                if (blocks[b].n_alleles > 1) ++interval_blocks;
                            const bool consec = interval_blocks == sc.size();
                            // An existing PAIRWISE edge covers it only when the scope is exactly
                            // two variable blocks with nothing variable between them.
                            const bool by_pair = sc.size() == 2 && consec;
                            if (by_pair) covered += kv.second; else need_higher += kv.second;
                            wi << ss << '\t' << sc.size() << '\t' << kv.second << '\t' << lo
                               << '\t' << hi << '\t' << lo << "-" << hi << '\t'
                               << interval_blocks << '\t' << (consec ? 1 : 0) << '\t'
                               << (by_pair ? 1 : 0) << '\n';
                        }
                        wi << "#total_wide_fragments\t" << n_wide << '\n';
                        wi << "#distinct_scopes\t" << by_scope.size() << '\n';
                        wi << "#fragments_needing_higher_order\t" << need_higher << '\n';
                        wi << "#fragments_a_pairwise_edge_could_take\t" << covered << '\n';
                        log.info("hybrid wide inventory: " + std::to_string(n_wide) +
                                 " Wide fragments over " + std::to_string(by_scope.size()) +
                                 " distinct scopes; " + std::to_string(need_higher) +
                                 " need a higher-order factor");
                    }
                    // ---- EVIDENCE ACCOUNTING for a candidate super-factor ----------------------
                    // Defined BEFORE any search is generalised: a factor whose input set is not
                    // pinned is a factor whose cost cannot be judged and whose double-counting
                    // cannot be excluded. Every fragment that would enter, exactly once, and every
                    // competing claim on it, enumerated rather than assumed.
                    if (!hybrid_super_ledger.empty()) {
                        std::vector<std::string> lp;
                        std::string cu;
                        for (char c : hybrid_super_ledger) {
                            if (c == ',') { lp.push_back(cu); cu.clear(); } else cu.push_back(c);
                        }
                        lp.push_back(cu);
                        std::vector<std::uint32_t> sb;
                        for (std::size_t q = 0; q + 1 < lp.size(); ++q) {
                            sb.push_back(static_cast<std::uint32_t>(std::stoul(lp[q])));
                        }
                        std::set<std::uint32_t> sbset(sb.begin(), sb.end());
                        std::ofstream lg(lp.back());
                        if (!lg) throw std::runtime_error("genotype: cannot write " + lp.back());
                        lg << "field\tvalue\n";
                        std::string bs;
                        for (std::size_t q = 0; q < sb.size(); ++q)
                            bs += (q ? "," : "") + std::to_string(sb[q]);
                        lg << "blocks\t" << bs << '\n';
                        // The pairwise edges the super-factor would subsume: consecutive pairs
                        // inside the block set.
                        std::vector<std::size_t> uni;
                        std::size_t pair_total = 0;
                        std::map<std::pair<std::uint32_t,std::uint32_t>, std::size_t> per_pair;
                        for (const auto& kv : by_edge) {
                            if (sbset.count(kv.first.first) && sbset.count(kv.first.second)) {
                                per_pair[kv.first] = kv.second.size();
                                pair_total += kv.second.size();
                                uni.insert(uni.end(), kv.second.begin(), kv.second.end());
                            }
                        }
                        for (const auto& pp : per_pair) {
                            lg << "owners_" << pp.first.first << "_" << pp.first.second << "\t"
                               << pp.second << '\n';
                        }
                        std::sort(uni.begin(), uni.end());
                        const std::size_t before = uni.size();
                        uni.erase(std::unique(uni.begin(), uni.end()), uni.end());
                        lg << "sum_of_pair_owner_counts\t" << pair_total << '\n';
                        lg << "union_size\t" << uni.size() << '\n';
                        // DISJOINTNESS is a property to check, not to hope for: a fragment owned by
                        // two pairwise edges would be counted twice by a naive union.
                        lg << "pairwise_sets_disjoint\t" << (before == uni.size() ? 1 : 0) << '\n';
                        // WIDE fragments whose variable scope lies inside the block set: they are
                        // currently unconsumed, and a super-factor spanning their scope is exactly
                        // what could consume them. Listed either way.
                        std::size_t wide_in = 0, wide_out = 0;
                        for (std::size_t fi = 0; fi < hf.size(); ++fi) {
                            if (owners[fi].kind != OwnerKind::Wide) continue;
                            bool inside = !owners[fi].var_scope.empty();
                            for (std::uint32_t b : owners[fi].var_scope)
                                if (!sbset.count(b)) { inside = false; break; }
                            if (inside) ++wide_in; else ++wide_out;
                        }
                        lg << "wide_scope_inside_blocks\t" << wide_in << '\n';
                        lg << "wide_scope_outside_blocks\t" << wide_out << '\n';
                        // NO OTHER CLAIM: none of these fragments may be Unary-owned or already
                        // consumed elsewhere. Ownership is a partition, so this must hold by
                        // construction -- which is why it is worth asserting rather than assuming.
                        std::size_t misclassified = 0;
                        for (std::size_t fi : uni)
                            if (owners[fi].kind != OwnerKind::Linkage) ++misclassified;
                        lg << "non_linkage_in_union\t" << misclassified << '\n';
                        lg << "every_fragment_exactly_once\t"
                           << ((before == uni.size() && misclassified == 0) ? 1 : 0) << '\n';
                        log.info("hybrid super ledger over blocks " + bs + ": union " +
                                 std::to_string(uni.size()) + " of " + std::to_string(pair_total) +
                                 " pairwise owners, disjoint=" +
                                 std::string(before == uni.size() ? "yes" : "NO") + ", wide inside " +
                                 std::to_string(wide_in));
                    }
                    // AN OPERATIONAL WORK GUARD on dense emission construction, until the
                    // support search exists. The cost is fragments x n_A x n_B window alignments
                    // per edge; C4 needs 1,154,642 and does not finish. This is a WORK count -- no
                    // likelihood, no score -- so exceeding it is a machine refusal, and the command
                    // says so quickly instead of appearing to hang.
                    std::size_t dense_emission_window_alignments = 0;
                    bool work_overflow = false;
                    for (const auto& kv : by_edge) {
                        const std::size_t na2 = blocks[kv.first.first].n_alleles;
                        const std::size_t nb2 = blocks[kv.first.second].n_alleles;
                        std::size_t w = 0;
                        if (__builtin_mul_overflow(na2, nb2, &w) ||
                            __builtin_mul_overflow(w, kv.second.size(), &w) ||
                            __builtin_add_overflow(dense_emission_window_alignments, w,
                                                   &dense_emission_window_alignments)) {
                            work_overflow = true;
                            break;
                        }
                    }
                    // THE DENSE FIGURE IS NOW A BASELINE, NOT A GATE. Production verifies only
                    // the windows the support search proposes, so bounding work by the dense
                    // cross-product would refuse work that is never done. It is still computed --
                    // the reduction is measured against it -- and an estimate that OVERFLOWS is
                    // still a refusal, because a wrapped baseline is not a baseline.
                    if (work_overflow) {
                        hyb_act = HybridActivation{};
                        hyb_act.refusal = "dense-emission-window-limit: the baseline estimate "
                                          "overflowed";
                        log.info("hybrid: " + hyb_act.refusal +
                                 "; no edge built, nothing subtracted");
                    } else {
                    // ONE BUDGET FOR THE LOCUS, charged incrementally inside every emission.
                    HybridWorkBudget work;
                    work.max_proposed_cells = hybrid_max_proposed_cells;
                    work.max_full_read_verifications = hybrid_max_full_read_verifications;
                    if (hybrid_verified_windows_alias_used) {
                        log.info("hybrid: --hybrid-max-verified-windows is deprecated; it sets "
                                 "--hybrid-max-proposed-cells (" +
                                 std::to_string(hybrid_max_proposed_cells) + "), the same "
                                 "quantity under an honest name");
                    }
                    std::vector<std::vector<std::string>> ballele(blocks.size());
                    for (std::size_t b = 0; b < blocks.size(); ++b) {
                        ballele[b] = blocks[b].allele_seq;
                    }
                    // SPARSE IS THE PRODUCTION BUILDER. The dense LinkageEdge is reachable only
                    // from fixtures and oracles now; production never consults its configuration
                    // cap, which is why a C4-sized edge no longer refuses here.
                    std::map<std::pair<std::uint32_t, std::uint32_t>, SparseLinkageEdge> edge_map;
                    SparseResourceLimits rlim;
                    rlim.max_classes = hybrid_max_classes;
                    rlim.max_bytes = hybrid_max_bytes;
                    std::vector<EdgeStatusEntry> edge_status;
                    const std::size_t FLANK = static_cast<std::size_t>(ip.hi);
                    for (const auto& kv : by_edge) {
                        EdgeStatusEntry es;
                        es.block_a = kv.first.first; es.block_b = kv.first.second;
                        es.n_fragments = kv.second.size();
                        const LinkageGeometry geom = build_linkage_geometry(
                            hyb_cov.frames, ballele, block_variable, es.block_a, es.block_b,
                            FLANK, ip);
                        if (!geom.ok) {
                            // THE REASON, kept. Recording only "not-computed" made the first C4
                            // refusal undiagnosable without re-reading the source.
                            es.status = LinkageStatus::NotComputed;
                            es.detail = geom.refusal;
                            edge_status.push_back(es);
                            continue;
                        }
                        // ---- EXPOSURE GEOMETRY, no emissions built -------------------------
                        // Cancellation is exact when exposure is AFFINE in the window length,
                        // because window_len = const + |A_alpha| + |B_beta| makes it additive and
                        // the crossed and straight sums coincide. Below hi-1 it is not affine, and
                        // no amount of support work repairs that -- so the question is purely how
                        // many configurations fall short, and by how much.
                        if (!hybrid_exposure_probe.empty()) {
                            static bool exp_hdr = false;
                            std::ofstream xp(hybrid_exposure_probe,
                                             exp_hdr ? std::ios::app : std::ios::trunc);
                            if (!exp_hdr) {
                                xp << "block_a\tblock_b\tn_a\tn_b\tconfigs\tinsert_hi"
                                      "\taffine_from\tmin_window\tmax_window\tconfigs_below"
                                      "\tshortfall_bp\tmin_exposure\tmax_exposure"
                                      "\tmax_asymmetry\tlflank\trflank\texposure_affine"
                                      "\tleft_block_common_suffix\tright_block_common_prefix"
                                      "\treachable_min_window\tcould_become_affine\n";
                                exp_hdr = true;
                            }
                            const std::size_t na2 = geom.alleles_a.size();
                            const std::size_t nb2 = geom.alleles_b.size();
                            std::size_t below = 0, minw = SIZE_MAX, maxw = 0;
                            double mine = 1e300, maxe = -1e300, asym = 0.0;
                            const std::size_t affine_from =
                                static_cast<std::size_t>(std::max<long>(0, ip.hi - 1));
                            for (std::size_t k = 0; k < geom.window_len.size(); ++k) {
                                minw = std::min(minw, geom.window_len[k]);
                                maxw = std::max(maxw, geom.window_len[k]);
                                mine = std::min(mine, geom.exposure[k]);
                                maxe = std::max(maxe, geom.exposure[k]);
                                if (geom.window_len[k] < affine_from) ++below;
                            }
                            for (std::size_t a1 = 0; a1 < na2; ++a1)
                            for (std::size_t b1 = 0; b1 < nb2; ++b1)
                            for (std::size_t a2 = 0; a2 < na2; ++a2)
                            for (std::size_t b2 = 0; b2 < nb2; ++b2) {
                                asym = std::max(asym, std::abs(
                                    (geom.exposure[a1 * nb2 + b1] + geom.exposure[a2 * nb2 + b2]) -
                                    (geom.exposure[a1 * nb2 + b2] + geom.exposure[a2 * nb2 + b1])));
                            }
                            // INVARIANT CONTEXT INSIDE A VARIABLE BLOCK. A zero derived flank
                            // means no complete fixed block sits outside the edge -- it does NOT
                            // mean the neighbouring block offers no candidate-independent
                            // sequence. Every allele of the block to the left may share a common
                            // SUFFIX, and every allele of the block to the right a common PREFIX;
                            // those are invariant even though the block is variable, and they are
                            // legitimate flank. Measured before concluding the refusal is genuine.
                            const auto common_suffix = [](const std::vector<std::string>& v) {
                                if (v.empty()) return static_cast<std::size_t>(0);
                                std::size_t n = v[0].size();
                                for (const std::string& x : v) n = std::min(n, x.size());
                                std::size_t k = 0;
                                while (k < n) {
                                    const char c = v[0][v[0].size() - 1 - k];
                                    bool same = true;
                                    for (const std::string& x : v)
                                        if (x[x.size() - 1 - k] != c) { same = false; break; }
                                    if (!same) break;
                                    ++k;
                                }
                                return k;
                            };
                            const auto common_prefix = [](const std::vector<std::string>& v) {
                                if (v.empty()) return static_cast<std::size_t>(0);
                                std::size_t n = v[0].size();
                                for (const std::string& x : v) n = std::min(n, x.size());
                                std::size_t k = 0;
                                while (k < n) {
                                    const char c = v[0][k];
                                    bool same = true;
                                    for (const std::string& x : v)
                                        if (x[k] != c) { same = false; break; }
                                    if (!same) break;
                                    ++k;
                                }
                                return k;
                            };
                            std::size_t left_inv = 0, right_inv = 0;
                            if (es.block_a > 0)
                                left_inv = common_suffix(blocks[es.block_a - 1].allele_seq);
                            if (es.block_b + 1 < blocks.size())
                                right_inv = common_prefix(blocks[es.block_b + 1].allele_seq);
                            const std::size_t reachable = minw + left_inv + right_inv;
                            xp << es.block_a << '\t' << es.block_b << '\t' << na2 << '\t' << nb2
                               << '\t' << geom.window_len.size() << '\t' << ip.hi << '\t'
                               << affine_from << '\t' << minw << '\t' << maxw << '\t' << below
                               << '\t' << (minw < affine_from ? affine_from - minw : 0) << '\t'
                               << mine << '\t' << maxe << '\t' << asym << '\t'
                               << geom.lflank.size() << '\t' << geom.rflank.size() << '\t'
                               << (geom.exposure_affine ? 1 : 0) << '\t'
                               << left_inv << '\t' << right_inv << '\t' << reachable << '\t'
                               << (reachable >= affine_from ? 1 : 0) << '\n';
                        }
                        // THE ALLELE INDEX, ONCE PER EDGE. Built from this edge's own piece
                        // length; every fragment then does hash lookups instead of re-scanning
                        // alleles up to 26 kb. Scanning per fragment is 1,274,940 full-allele scans
                        // on C4 and does not finish.
                        AlleleProductIndex aidx;
                        if (!kv.second.empty()) {
                            const Fragment& f0 = hf[kv.second.front()];
                            const std::size_t d0 =
                                mate_band_edits(hyb_params.max_divergence, f0.r1.size());
                            const std::size_t p0 = f0.r1.size() / (d0 + 1);
                            if (p0 >= 8 && p0 <= 32) aidx = build_allele_product_index(geom, p0);
                        }
                        // ---- ONE-EDGE ORACLE COMPARISON ------------------------------------
                        // Counters can show speed; only this can show the likelihood survived.
                        // The supported search and the dense oracle must agree on WHICH cells are
                        // finite, on HOW MANY states built each one, and on the mass exactly.
                        if (!hybrid_edge_oracle.empty()) {
                            std::size_t c1 = hybrid_edge_oracle.find(',');
                            std::size_t c2 = hybrid_edge_oracle.find(',', c1 + 1);
                            if (c1 != std::string::npos && c2 != std::string::npos &&
                                std::stoul(hybrid_edge_oracle.substr(0, c1)) == es.block_a &&
                                std::stoul(hybrid_edge_oracle.substr(c1 + 1, c2 - c1 - 1)) ==
                                    es.block_b) {
                                const std::string path = hybrid_edge_oracle.substr(c2 + 1);
                                std::ofstream eo(path);
                                const auto t0 = std::chrono::steady_clock::now();
                                HybridWorkBudget wb;
                                wb.max_proposed_cells = hybrid_max_proposed_cells;
                                wb.max_full_read_verifications = hybrid_max_full_read_verifications;
                                std::size_t cells_differ = 0, mult_differ = 0, fin_s = 0, fin_d = 0,
                                            fb = 0, refused = 0;
                                double worst = 0.0;
                                std::uint64_t states_s = 0, states_d = 0;
                                eo << "distinct_masses\tmass_spread\t"
                                      "fragment\tdense_pairs\tproposed\tunique_seed_starts"
                                      "\tfull_read_verifications\taccepted_placements"
                                      "\tverified_fr_states\tfinite_cells\tfallback\n";
                                for (std::size_t fi : kv.second) {
                                    const std::size_t len = hf[fi].bases();
                                    const std::size_t bee = static_cast<std::size_t>(
                                        hyb_params.bg_divergence * static_cast<double>(len));
                                    const double bgf = static_cast<double>(bee) * lep +
                                                       static_cast<double>(len - bee) * l1m;
                                    AlleleProductSupport sp;
                                    const LinkageEmission S = linkage_emission_supported(
                                        hf[fi], geom, ip, hyb_params.max_divergence, lep, l1m, bgf,
                                        &sp, aidx.ok ? &aidx : nullptr, &wb);
                                    const LinkageEmission D = linkage_emission(
                                        hf[fi], geom, ip, hyb_params.max_divergence, lep, l1m, bgf);
                                    if (S.work_refused) { ++refused; continue; }
                                    if (sp.exhaustive_fallback) ++fb;
                                    for (std::size_t k = 0; k < D.mass.size(); ++k) {
                                        const bool fx = D.mass[k] != -INFINITY;
                                        const bool fy = S.mass[k] != -INFINITY;
                                        if (fx) ++fin_d;
                                        if (fy) ++fin_s;
                                        states_d += D.cell_states[k];
                                        states_s += S.cell_states[k];
                                        if (fx != fy) { ++cells_differ; continue; }
                                        if (D.cell_states[k] != S.cell_states[k]) ++mult_differ;
                                        if (fx) worst = std::max(worst,
                                                                 std::abs(D.mass[k] - S.mass[k]));
                                    }
                                    // PER FRAGMENT, not only totals: a mean hides whether a few
                                    // fragments propose everything or all of them propose a lot.
                                    // HOW MANY DISTINCT VALUES the mass actually takes over the
                                    // finite cells. If 126 beta cells carry two values, the support
                                    // is not reducible but the ALLELE SET is -- and that is a
                                    // different optimisation from a better seed filter.
                                    std::vector<double> vals;
                                    vals.reserve(D.mass.size());
                                    for (double m : D.mass)
                                        if (m != -INFINITY) vals.push_back(m);
                                    std::sort(vals.begin(), vals.end());
                                    const std::size_t n_distinct =
                                        static_cast<std::size_t>(
                                            std::unique(vals.begin(), vals.end(),
                                                        [](double x, double y) {
                                                            return std::abs(x - y) < 1e-9;
                                                        }) - vals.begin());
                                    const double spread = vals.empty() ? 0.0
                                                        : vals.back() - vals.front();
                                    eo << n_distinct << '\t' << spread << '\t'
                                       << hf[fi].name << '\t' << sp.dense_pairs << '\t'
                                       << sp.proposals.size() << '\t' << sp.unique_seed_starts
                                       << '\t' << S.full_read_verifications << '\t'
                                       << S.accepted_mate_placements << '\t'
                                       << S.verified_fr_states << '\t' << S.finite_emission_cells
                                       << '\t' << (sp.exhaustive_fallback ? 1 : 0) << '\n';
                                }
                                // ONE FRAGMENT'S CELLS IN FULL, against the B allele LENGTH.
                                // If every beta is finite but the mass still varies, the question
                                // is what beta is varying THROUGH. Length moves the insert; content
                                // moves the Hamming term. Dumping both makes that separable
                                // instead of a matter of opinion.
                                if (!kv.second.empty()) {
                                    const std::size_t f0i = kv.second.front();
                                    const std::size_t len0 = hf[f0i].bases();
                                    const std::size_t be0 = static_cast<std::size_t>(
                                        hyb_params.bg_divergence * static_cast<double>(len0));
                                    const double bgf0 = static_cast<double>(be0) * lep +
                                                        static_cast<double>(len0 - be0) * l1m;
                                    const LinkageEmission C = linkage_emission(
                                        hf[f0i], geom, ip, hyb_params.max_divergence, lep, l1m,
                                        bgf0);
                                    std::ofstream cf(path + ".cells.tsv");
                                    cf << "alpha\tbeta\ta_len\tb_len\tmass\tstates\n";
                                    cf.setf(std::ios::fixed);
                                    cf.precision(12);
                                    for (std::size_t al = 0; al < C.n_a; ++al) {
                                        for (std::size_t be = 0; be < C.n_b; ++be) {
                                            cf << al << '\t' << be << '\t'
                                               << geom.alleles_a[al].size() << '\t'
                                               << geom.alleles_b[be].size() << '\t'
                                               << C.mass[al * C.n_b + be] << '\t'
                                               << C.cell_states[al * C.n_b + be] << '\n';
                                        }
                                    }
                                }
                                const double secs = std::chrono::duration<double>(
                                    std::chrono::steady_clock::now() - t0).count();
                                eo << "#edge\t" << es.block_a << '-' << es.block_b << '\n'
                                   << "#owned_fragments\t" << kv.second.size() << '\n'
                                   << "#alleles\t" << blocks[es.block_a].n_alleles << 'x'
                                   << blocks[es.block_b].n_alleles << '\n'
                                   << "#total_possible_cells\t"
                                   << static_cast<std::uint64_t>(blocks[es.block_a].n_alleles) *
                                      blocks[es.block_b].n_alleles * kv.second.size() << '\n'
                                   << "#proposed_cells\t" << wb.proposed_cells << '\n'
                                   << "#full_read_verifications\t" << wb.full_read_verifications
                                   << '\n'
                                   << "#bases_compared_upper_bound\t"
                                   << wb.bases_compared_upper_bound << '\n'
                                   << "#fallbacks\t" << fb << '\n'
                                   << "#work_refused_fragments\t" << refused << '\n'
                                   << "#index_ok\t" << (aidx.ok ? 1 : 0) << '\n'
                                   << "#index_complete\t" << (aidx.complete ? 1 : 0) << '\n'
                                   << "#finite_cells_supported\t" << fin_s << '\n'
                                   << "#finite_cells_dense\t" << fin_d << '\n'
                                   << "#states_supported\t" << states_s << '\n'
                                   << "#states_dense\t" << states_d << '\n'
                                   << "#cells_differ\t" << cells_differ << '\n'
                                   << "#multiplicity_differ\t" << mult_differ << '\n'
                                   << "#worst_mass_diff\t" << worst << '\n'
                                   << "#seconds\t" << secs << '\n'
                                   << "#peak_rss_mb\t" << [] {
                                          struct rusage r {};
                                          if (getrusage(RUSAGE_SELF, &r) != 0) return 0.0;
#ifdef __APPLE__
                                          return static_cast<double>(r.ru_maxrss) / 1048576.0;
#else
                                          return static_cast<double>(r.ru_maxrss) / 1024.0;
#endif
                                      }() << '\n';
                                log.info("hybrid edge oracle " + std::to_string(es.block_a) + "-" +
                                         std::to_string(es.block_b) + ": " +
                                         std::to_string(cells_differ) + " cells differ, " +
                                         std::to_string(mult_differ) + " multiplicities differ, "
                                         "worst mass " + std::to_string(worst) + ", " +
                                         std::to_string(fin_s) + " of " +
                                         std::to_string(fin_d) + " finite, " +
                                         std::to_string(secs) + " s");
                            }
                        }
                        // ---- JOINT EMISSION SIGNATURES -------------------------------------
                        // Per-fragment classes can CUT ACROSS one another, so their sum is not the
                        // achievable compression and pooled equal masses are not shareable at all.
                        // Two allele pairs may be collapsed only when the WHOLE vector of
                        // per-fragment signatures is identical, which is what this intersects.
                        if (!hybrid_edge_signature.empty()) {
                            const std::size_t s1 = hybrid_edge_signature.find(',');
                            const std::size_t s2 = s1 == std::string::npos ? s1
                                : hybrid_edge_signature.find(',', s1 + 1);
                            if (s1 != std::string::npos && s2 != std::string::npos &&
                                std::stoul(hybrid_edge_signature.substr(0, s1)) == es.block_a &&
                                std::stoul(hybrid_edge_signature.substr(s1 + 1, s2 - s1 - 1)) ==
                                    es.block_b) {
                                const std::string spath = hybrid_edge_signature.substr(s2 + 1);
                                const auto ts = std::chrono::steady_clock::now();
                                const std::size_t ncell = static_cast<std::size_t>(
                                    blocks[es.block_a].n_alleles) * blocks[es.block_b].n_alleles;
                                std::vector<std::string> joint(ncell);
                                std::ofstream sf(spath);
                                sf << "fragment\tfinite_cells\tdistinct_signatures\n";
                                HybridWorkBudget sb;
                                sb.max_proposed_cells = hybrid_max_proposed_cells;
                                sb.max_full_read_verifications =
                                    hybrid_max_full_read_verifications;
                                std::uint64_t sum_per_fragment = 0;
                                std::size_t n_fb = 0, n_refused = 0, n_badsize = 0,
                                            n_contributed = 0;
                                std::vector<std::size_t> per_frag;
                                for (std::size_t fi : kv.second) {
                                    const std::size_t len = hf[fi].bases();
                                    const std::size_t bee = static_cast<std::size_t>(
                                        hyb_params.bg_divergence * static_cast<double>(len));
                                    const double bgf = static_cast<double>(bee) * lep +
                                                       static_cast<double>(len - bee) * l1m;
                                    AlleleProductSupport sp;
                                    std::vector<std::string> cs;
                                    const LinkageEmission S = linkage_emission_supported(
                                        hf[fi], geom, ip, hyb_params.max_divergence, lep, l1m, bgf,
                                        &sp, aidx.ok ? &aidx : nullptr, &sb, &cs);
                                    // FAIL CLOSED. A fragment that falls back, runs out of budget
                                    // or returns a malformed vector is COUNTED, never skipped: a
                                    // joint class computed over a subset of the fragments is not a
                                    // smaller answer, it is a wrong one -- and silently smaller in
                                    // exactly the flattering direction.
                                    if (sp.exhaustive_fallback) { ++n_fb; continue; }
                                    if (S.work_refused) { ++n_refused; continue; }
                                    if (cs.size() != ncell) { ++n_badsize; continue; }
                                    ++n_contributed;
                                    std::unordered_set<std::string> own(cs.begin(), cs.end());
                                    sum_per_fragment += own.size();
                                    per_frag.push_back(own.size());
                                    std::size_t fin = 0;
                                    for (std::size_t k = 0; k < ncell; ++k) {
                                        if (!cs[k].empty()) ++fin;
                                        // The joint key is the CONCATENATION over fragments, with a
                                        // length prefix so two different splits cannot alias.
                                        const std::uint32_t n =
                                            static_cast<std::uint32_t>(cs[k].size());
                                        joint[k].append(reinterpret_cast<const char*>(&n), 4);
                                        joint[k].append(cs[k]);
                                    }
                                    sf << hf[fi].name << '\t' << fin << '\t' << own.size() << '\n';
                                }
                                const bool sound = n_contributed == kv.second.size() &&
                                                   n_fb == 0 && n_refused == 0 && n_badsize == 0;
                                // SIGNATURE IDs, so the phase pattern can be keyed by them.
                                std::unordered_map<std::string, std::uint32_t> sid;
                                std::vector<std::uint32_t> qid(ncell, 0);
                                for (std::size_t k = 0; k < ncell; ++k) {
                                    const auto it = sid.emplace(joint[k],
                                        static_cast<std::uint32_t>(sid.size())).first;
                                    qid[k] = it->second;
                                }
                                const std::uint64_t jc = sid.size();
                                // THE FOUR-SIGNATURE PHASE PATTERN, which is what psi actually
                                // stores. K joint haploid signatures bound the RAW diploid score
                                // pairs by K^2, but a class is keyed by an unordered A-pair and an
                                // unordered B-pair and its delta is built from FOUR corners --
                                // straight {q11,q22} against crossed {q12,q21}, each unordered
                                // because the mixture is symmetric in its two arguments. Two
                                // content classes can share straight corners and differ on crossed,
                                // so K^2 is not the count of production factor states.
                                const std::size_t NA = blocks[es.block_a].n_alleles;
                                const std::size_t NB = blocks[es.block_b].n_alleles;
                                // THE TRAVERSAL QUESTION, which is separate from the scoring one.
                                // A quadruple's pattern depends on q only through four cells, so
                                // two A alleles with IDENTICAL signature rows are interchangeable
                                // in every pattern, and likewise two B alleles with identical
                                // columns. If the matrix has RA distinct rows and RB distinct
                                // columns, the enumeration is C(RA,2)*C(RB,2)-shaped instead of
                                // C(na,2)*C(nb,2) -- and nothing ever visits the full table.
                                // Reducing distinct DELTA VALUES to a small number would not by
                                // itself achieve that; this is what does.
                                std::unordered_map<std::string, std::uint32_t> rowid, colid;
                                std::vector<std::uint32_t> arow(NA), bcol(NB);
                                for (std::size_t a = 0; a < NA; ++a) {
                                    std::string key(NB * 4, '\0');
                                    for (std::size_t b = 0; b < NB; ++b)
                                        std::memcpy(&key[b * 4], &qid[a * NB + b], 4);
                                    arow[a] = rowid.emplace(key,
                                        static_cast<std::uint32_t>(rowid.size())).first->second;
                                }
                                for (std::size_t b = 0; b < NB; ++b) {
                                    std::string key(NA * 4, '\0');
                                    for (std::size_t a = 0; a < NA; ++a)
                                        std::memcpy(&key[a * 4], &qid[a * NB + b], 4);
                                    bcol[b] = colid.emplace(key,
                                        static_cast<std::uint32_t>(colid.size())).first->second;
                                }
                                const std::size_t RA = rowid.size(), RB = colid.size();
                                // TWO representatives per class where the class has two members.
                                // A quadruple whose two A alleles come from the SAME row class is a
                                // real class -- and a FLAT one, since identical rows force
                                // q11==q12 and q21==q22, hence straight == crossed and delta 0.
                                // Visiting only distinct classes would silently drop every one of
                                // them, which is what the first version of this check did.
                                std::vector<std::size_t> arep(RA, SIZE_MAX), arep2(RA, SIZE_MAX);
                                std::vector<std::size_t> brep(RB, SIZE_MAX), brep2(RB, SIZE_MAX);
                                for (std::size_t a = 0; a < NA; ++a) {
                                    if (arep[arow[a]] == SIZE_MAX) arep[arow[a]] = a;
                                    else if (arep2[arow[a]] == SIZE_MAX) arep2[arow[a]] = a;
                                }
                                for (std::size_t b = 0; b < NB; ++b) {
                                    if (brep[bcol[b]] == SIZE_MAX) brep[bcol[b]] = b;
                                    else if (brep2[bcol[b]] == SIZE_MAX) brep2[bcol[b]] = b;
                                }
                                // ORIENTATION MATTERS, and representatives alone cannot express
                                // it. The full scan keys a quadruple with a1<a2 and b1<b2, so which
                                // corner counts as STRAIGHT depends on the index order between the
                                // two classes' members -- and both orders occur, because class
                                // members are scattered through index space. Visiting one
                                // representative pair per unordered class pair therefore realises
                                // only one orientation and silently drops delta's mirror image.
                                // So: enumerate ORDERED class pairs, and keep one only when some
                                // real member pair actually realises that order.
                                std::vector<std::size_t> amin(RA, SIZE_MAX), amax(RA, 0),
                                                         acount(RA, 0);
                                std::vector<std::size_t> bmin(RB, SIZE_MAX), bmax(RB, 0),
                                                         bcount(RB, 0);
                                for (std::size_t a = 0; a < NA; ++a) {
                                    amin[arow[a]] = std::min(amin[arow[a]], a);
                                    amax[arow[a]] = std::max(amax[arow[a]], a);
                                    ++acount[arow[a]];
                                }
                                for (std::size_t b = 0; b < NB; ++b) {
                                    bmin[bcol[b]] = std::min(bmin[bcol[b]], b);
                                    bmax[bcol[b]] = std::max(bmax[bcol[b]], b);
                                    ++bcount[bcol[b]];
                                }
                                // Ordered class pair (i,j) is realisable as a1<a2 iff some member
                                // of i precedes some member of j; within one class it needs two.
                                std::vector<std::pair<std::size_t, std::size_t>> apairs, bpairs;
                                for (std::size_t i = 0; i < RA; ++i)
                                for (std::size_t j = 0; j < RA; ++j) {
                                    const bool okp = (i == j) ? acount[i] >= 2 : amin[i] < amax[j];
                                    if (okp) apairs.emplace_back(i, j);
                                }
                                for (std::size_t u = 0; u < RB; ++u)
                                for (std::size_t v = 0; v < RB; ++v) {
                                    const bool okp = (u == v) ? bcount[u] >= 2 : bmin[u] < bmax[v];
                                    if (okp) bpairs.emplace_back(u, v);
                                }
                                std::unordered_set<std::uint64_t> patterns;
                                std::size_t affected = 0, flat_full = 0;
                                for (std::size_t a1 = 0; a1 + 1 < NA; ++a1)
                                for (std::size_t a2 = a1 + 1; a2 < NA; ++a2)
                                for (std::size_t b1 = 0; b1 + 1 < NB; ++b1)
                                for (std::size_t b2 = b1 + 1; b2 < NB; ++b2) {
                                    const std::uint32_t q11 = qid[a1 * NB + b1];
                                    const std::uint32_t q22 = qid[a2 * NB + b2];
                                    const std::uint32_t q12 = qid[a1 * NB + b2];
                                    const std::uint32_t q21 = qid[a2 * NB + b1];
                                    // A class matters only if some corner carries mass; an
                                    // all-empty quadruple contributes nothing to psi.
                                    if (joint[a1 * NB + b1].find_first_not_of('\0') ==
                                            std::string::npos &&
                                        joint[a2 * NB + b2].find_first_not_of('\0') ==
                                            std::string::npos &&
                                        joint[a1 * NB + b2].find_first_not_of('\0') ==
                                            std::string::npos &&
                                        joint[a2 * NB + b1].find_first_not_of('\0') ==
                                            std::string::npos) continue;
                                    ++affected;
                                    const std::uint64_t sA = std::min(q11, q22);
                                    const std::uint64_t sB = std::max(q11, q22);
                                    const std::uint64_t cA = std::min(q12, q21);
                                    const std::uint64_t cB = std::max(q12, q21);
                                    patterns.insert((sA << 48) | (sB << 32) | (cA << 16) | cB);
                                    if (sA == cA && sB == cB) ++flat_full;
                                }
                                // THE COLLAPSED ENUMERATION, over class REPRESENTATIVES only.
                                // Its pattern set must equal the one from full enumeration, or the
                                // row/column collapse is not equivalence-preserving. Distinct
                                // classes may still share a representative pair when RA<NA, so
                                // ordered representative pairs are used, both orders, to reach the
                                // patterns a strict a1<a2 scan over representatives would miss.
                                std::unordered_set<std::uint64_t> patterns_collapsed;
                                std::size_t collapsed_visits = 0;
                                for (const auto& ap : apairs)
                                for (const auto& bp : bpairs) {
                                    // CLASS indices now, in the realised order: the first element
                                    // plays the role the smaller-indexed allele plays in the full
                                    // scan. No re-sorting by allele index, which is what discarded
                                    // the orientation before.
                                    const std::size_t A1 = arep[ap.first], A2 = arep[ap.second];
                                    const std::size_t B1 = brep[bp.first], B2 = brep[bp.second];
                                    ++collapsed_visits;
                                    const std::uint32_t q11 = qid[A1 * NB + B1];
                                    const std::uint32_t q22 = qid[A2 * NB + B2];
                                    const std::uint32_t q12 = qid[A1 * NB + B2];
                                    const std::uint32_t q21 = qid[A2 * NB + B1];
                                    // THE SAME emptiness test as the full scan, byte for byte. A
                                    // joint key is never zero-length -- it always carries one
                                    // length prefix per fragment -- so testing .empty() here while
                                    // the full scan tests for all-zero bytes would compare two
                                    // different sets and call the difference a finding.
                                    if (joint[A1 * NB + B1].find_first_not_of('\0') ==
                                            std::string::npos &&
                                        joint[A2 * NB + B2].find_first_not_of('\0') ==
                                            std::string::npos &&
                                        joint[A1 * NB + B2].find_first_not_of('\0') ==
                                            std::string::npos &&
                                        joint[A2 * NB + B1].find_first_not_of('\0') ==
                                            std::string::npos) continue;
                                    const std::uint64_t sA = std::min(q11, q22);
                                    const std::uint64_t sB = std::max(q11, q22);
                                    const std::uint64_t cA = std::min(q12, q21);
                                    const std::uint64_t cB = std::max(q12, q21);
                                    patterns_collapsed.insert(
                                        (sA << 48) | (sB << 32) | (cA << 16) | cB);
                                }
                                std::sort(per_frag.begin(), per_frag.end());
                                const double ssec = std::chrono::duration<double>(
                                    std::chrono::steady_clock::now() - ts).count();
                                sf << "#edge\t" << es.block_a << '-' << es.block_b << '\n'
                                   << "#alleles\t" << blocks[es.block_a].n_alleles << 'x'
                                   << blocks[es.block_b].n_alleles << '\n'
                                   << "#owned_fragments\t" << kv.second.size() << '\n'
                                   << "#fragments_contributed\t" << n_contributed << '\n'
                                   << "#fallback_fragments\t" << n_fb << '\n'
                                   << "#work_refused_fragments\t" << n_refused << '\n'
                                   << "#bad_signature_size\t" << n_badsize << '\n'
                                   << "#sound\t" << (sound ? 1 : 0) << '\n'
                                   << "#allele_pair_cells\t" << ncell << '\n'
                                   << "#per_fragment_classes_min\t"
                                   << (per_frag.empty() ? 0 : per_frag.front()) << '\n'
                                   << "#per_fragment_classes_median\t"
                                   << (per_frag.empty() ? 0 : per_frag[per_frag.size() / 2]) << '\n'
                                   << "#per_fragment_classes_max\t"
                                   << (per_frag.empty() ? 0 : per_frag.back()) << '\n'
                                   << "#sum_per_fragment_classes\t" << sum_per_fragment << '\n'
                                   << "#joint_classes\t" << jc << '\n'
                                   << "#joint_compression\t"
                                   << (jc ? static_cast<double>(ncell) / jc : 0.0) << '\n'
                                   << "#raw_diploid_signature_pairs_upper_bound\t" << jc * jc
                                   << '\n'
                                   << "#affected_classes\t" << affected << '\n'
                                   << "#distinct_phase_patterns\t" << patterns.size() << '\n'
                                   << "#signature_rows_A\t" << RA << " of " << NA << '\n'
                                   << "#signature_cols_B\t" << RB << " of " << NB << '\n'
                                   << "#collapsed_visits\t" << collapsed_visits << '\n'
                                   << "#collapsed_patterns\t" << patterns_collapsed.size() << '\n'
                                   << "#flat_quadruples_full\t" << flat_full << '\n'
                                   << "#collapse_equivalent\t"
                                   << (patterns_collapsed == patterns ? 1 : 0) << '\n'
                                   << "#traversal_reduction\t"
                                   << (collapsed_visits ? static_cast<double>(affected) /
                                                          collapsed_visits : 0.0) << '\n'
                                   << "#full_read_verifications\t" << sb.full_read_verifications
                                   << '\n'
                                   << "#seconds\t" << ssec << '\n';
                                if (!sound) {
                                    sf << "#REFUSED\tjoint classes not emitted: "
                                          "contributed " << n_contributed << " of "
                                       << kv.second.size() << ", fallbacks " << n_fb
                                       << ", work-refused " << n_refused
                                       << ", malformed " << n_badsize << '\n';
                                }
                                if (!sound) {
                                    throw std::runtime_error(
                                        "genotype: --hybrid-edge-signature refuses on edge " +
                                        std::to_string(es.block_a) + "-" +
                                        std::to_string(es.block_b) + ": " +
                                        std::to_string(n_contributed) + " of " +
                                        std::to_string(kv.second.size()) + " fragments "
                                        "contributed (fallbacks " + std::to_string(n_fb) +
                                        ", work-refused " + std::to_string(n_refused) +
                                        ", malformed " + std::to_string(n_badsize) + "). A joint "
                                        "class over a subset of the fragments is not a result.");
                                }
                                log.info("hybrid edge signature " + std::to_string(es.block_a) +
                                         "-" + std::to_string(es.block_b) + ": " +
                                         std::to_string(ncell) + " cells, per-fragment classes " +
                                         std::to_string(per_frag.empty() ? 0 : per_frag.front()) +
                                         "/" + std::to_string(per_frag.empty() ? 0
                                                : per_frag[per_frag.size() / 2]) + "/" +
                                         std::to_string(per_frag.empty() ? 0 : per_frag.back()) +
                                         " (min/med/max), JOINT " + std::to_string(jc) +
                                         ", phase patterns " + std::to_string(patterns.size()) +
                                         " over " + std::to_string(affected) + " affected classes; "
                                         "signature rows/cols " + std::to_string(RA) + "/" +
                                         std::to_string(RB) + ", collapsed visits " +
                                         std::to_string(collapsed_visits) + ", equivalent=" +
                                         (patterns_collapsed == patterns ? "yes" : "NO") + ", " +
                                         std::to_string(ssec) + " s");
                            }
                        }
                        std::vector<LinkageEmission> ems;
                        std::vector<std::vector<std::string>> edge_sigs;
                        ems.reserve(kv.second.size());
                        std::size_t n_inform = 0;
                        std::size_t e_seed_hits = 0, e_proposed = 0, e_verified = 0,
                                    e_dense = 0, e_fallback = 0, e_reads = 0, e_finite = 0;
                        for (std::size_t fi : kv.second) {
                            const std::size_t len = hf[fi].bases();
                            const std::size_t be = static_cast<std::size_t>(
                                hyb_params.bg_divergence * static_cast<double>(len));
                            const double bgf = static_cast<double>(be) * lep +
                                               static_cast<double>(len - be) * l1m;
                            AlleleProductSupport sup;
                            std::vector<std::string> csig;
                            ems.push_back(linkage_emission_supported(
                                hf[fi], geom, ip, hyb_params.max_divergence, lep, l1m, bgf, &sup,
                                aidx.ok ? &aidx : nullptr, &work,
                                hybrid_grouped ? &csig : nullptr));
                            if (hybrid_grouped) edge_sigs.push_back(std::move(csig));
                            if (ems.back().work_refused) break;
                            if (ems.back().informative) ++n_inform;
                            // THE SEARCH'S OWN COUNTERS, not the dense budget's. Seed hits,
                            // proposed states and verified windows describe what production
                            // actually did; the dense figure is kept only as the baseline it is
                            // measured against.
                            e_seed_hits += sup.seed_hits;
                            e_proposed += sup.proposed_states;
                            e_verified += sup.exhaustive_fallback ? sup.dense_pairs
                                                                  : sup.proposals.size();
                            e_dense += sup.dense_pairs;
                            e_reads += ems.back().full_read_verifications;
                            e_finite += ems.back().finite_emission_cells;
                            if (sup.exhaustive_fallback) ++e_fallback;
                        }
                        const auto t_build = std::chrono::steady_clock::now();
                        GroupedBuildStats gstats;
                        SparseLinkageEdge E =
                            hybrid_grouped
                                ? build_sparse_linkage_edge_grouped(
                                      ems, edge_sigs, geom, hyb_params.lambda,
                                      std::log1p(-hyb_params.outlier_mix),
                                      std::log(hyb_params.outlier_mix), rlim, &gstats)
                                : build_sparse_linkage_edge(
                                      ems, geom, hyb_params.lambda,
                                      std::log1p(-hyb_params.outlier_mix),
                                      std::log(hyb_params.outlier_mix), rlim);
                        if (hybrid_grouped && !hybrid_grouped_report.empty()) {
                            static bool gh = false;
                            std::ofstream gr(hybrid_grouped_report,
                                             gh ? std::ios::app : std::ios::trunc);
                            if (!gh) {
                                gr << "block_a\tblock_b\tn_a\tn_b\towned_fragments"
                                      "\tfragments_contributed\toracle_visits"
                                      "\trepresentative_visits\tpattern_evaluations"
                                      "\trow_classes\tcol_classes\tstored_class_quadruples"
                                      "\tpredicted_classes\tsupport_cells\tbytes_cell_signatures"
                                      "\tbytes_matrix\tbytes_members\tbytes_delta_class"
                                      "\tbytes_total\testimates_checked\tbuild_seconds"
                                      "\tstatus\n";
                                gh = true;
                            }
                            gr << es.block_a << '\t' << es.block_b << '\t'
                               << blocks[es.block_a].n_alleles << '\t'
                               << blocks[es.block_b].n_alleles << '\t' << kv.second.size() << '\t'
                               << edge_sigs.size() << '\t' << gstats.oracle_visits << '\t'
                               << gstats.representative_visits << '\t'
                               << gstats.pattern_evaluations << '\t' << gstats.row_classes << '\t'
                               << gstats.col_classes << '\t' << E.delta_class.size() << '\t'
                               << E.predicted_classes << '\t' << E.support_cells << '\t'
                               << gstats.bytes_cell_signatures << '\t' << gstats.bytes_matrix
                               << '\t' << gstats.bytes_members << '\t'
                               << gstats.bytes_delta_class << '\t' << E.bytes_total() << '\t'
                               << (gstats.estimates_checked ? 1 : 0) << '\t'
                               << gstats.build_seconds << '\t'
                               << linkage_status_name(E.status) << '\n';
                        }
                        const double build_s = std::chrono::duration<double>(
                            std::chrono::steady_clock::now() - t_build).count();
                        // PER EDGE, AS IT COMPLETES. A silent multi-minute loop is what made every
                        // previous attempt uninterpretable until it was killed.
                        log.info("hybrid: edge " + std::to_string(es.block_a) + "-" +
                                 std::to_string(es.block_b) + " " +
                                 std::to_string(blocks[es.block_a].n_alleles) + "x" +
                                 std::to_string(blocks[es.block_b].n_alleles) + ", " +
                                 std::to_string(kv.second.size()) + " frags, " +
                                 std::to_string(e_verified) + " verified of " +
                                 std::to_string(e_dense) + " dense, " +
                                 std::to_string(e_fallback) + " fallbacks, " +
                                 std::to_string(e_reads) + " read checks -> " +
                                 std::to_string(e_finite) + " finite cells, idx=" +
                                 (aidx.ok ? (aidx.complete ? "ok" : "incomplete") : "none") + ", " +
                                 std::to_string(build_s) + " s");
                        es.status = E.status;
                        // PER-EDGE MEASUREMENT, from the production build itself rather than a
                        // second profiler that would rebuild the support independently and could
                        // measure something the caller never used.
                        edge_rows.push_back(EdgeRow{
                            es.block_a, es.block_b, blocks[es.block_a].n_alleles,
                            blocks[es.block_b].n_alleles, kv.second.size(), n_inform,
                            E.support_cells, E.stored_classes, E.predicted_classes,
                            E.theoretical_configs, E.bytes_total(), build_s,
                            linkage_status_name(E.status), e_seed_hits, e_proposed, e_verified,
                            e_dense, e_fallback, e_reads, e_finite});
                        edge_status.push_back(es);
                        if (E.usable()) edge_map.emplace(kv.first, std::move(E));
                        // A BUDGET ON THE WORK PRODUCTION ACTUALLY DOES: windows verified, summed
                        // over edges. Checked as it accumulates so an expensive locus refuses
                        // early rather than after paying for every edge.
                        verified_windows_total += e_verified;
                        // THE BUDGET IS ALREADY ENFORCED INSIDE THE SEARCH. What is left here is
                        // only to propagate its refusal transactionally: an exhausted budget means
                        // some emission carries no mass, and a model missing an emission is not a
                        // smaller model.
                        if (work.exhausted) { work_overflow = true; break; }
                    }
                    hyb_work_proposed = work.proposed_cells;
                    hyb_work_verifications = work.full_read_verifications;
                    hyb_work_bases = work.bases_compared_upper_bound;
                    if (work_overflow && !edge_status.empty()) {
                        hyb_act = HybridActivation{};
                        hyb_act.refusal = work.exhausted
                            ? work.reason + ": " + std::to_string(work.proposed_cells) +
                              " proposed cells, " +
                              std::to_string(work.full_read_verifications) +
                              " full-read verifications (limits " +
                              std::to_string(hybrid_max_proposed_cells) + " / " +
                              std::to_string(hybrid_max_full_read_verifications) + ")"
                            : "verified-window-limit: the baseline estimate overflowed";
                        log.info("hybrid: " + hyb_act.refusal +
                                 "; no edge activated, nothing subtracted");
                    } else {
                    // Haplotype -> allele per block, validated at the int -> unsigned boundary.
                    std::vector<AlleleMapping> maps(blocks.size());
                    for (std::size_t b = 0; b < blocks.size(); ++b) {
                        std::vector<int> av(hap_names.size(), -1);
                        for (std::size_t h = 0; h < hap_names.size(); ++h) {
                            const auto it = blocks[b].allele_of.find(hap_names[h]);
                            if (it != blocks[b].allele_of.end()) {
                                av[h] = static_cast<int>(it->second);
                            }
                        }
                        maps[b] = build_allele_mapping(av, blocks[b].n_alleles,
                                                       blocks[b].bypass_allele);
                    }
                    // ---- WHY EACH UNUSABLE FRAGMENT IS UNUSABLE ---------------------------
                    // One policy cannot serve all of these. A fragment with NO IN-BAND ORIGINS is
                    // candidate-INDEPENDENT and cannot rank genotypes; one whose scope merely
                    // failed to certify has attributable mass and needs a deeper search; one with
                    // dominant UNMAPPED mass is candidate-DEPENDENT evidence the block chain
                    // cannot express, and only that last category may legitimately block the call.
                    {
                        struct RCat {
                            std::size_t n = 0, bases = 0;
                            double in_band = -std::numeric_limits<double>::infinity();
                            double omitted = -std::numeric_limits<double>::infinity();
                            double unmapped = -std::numeric_limits<double>::infinity();
                            std::size_t origins = 0, scope_blocks = 0, max_span = 0;
                        };
                        std::map<std::string, RCat> cats;
                        const auto la = [](double x, double y) {
                            const double ninf = -std::numeric_limits<double>::infinity();
                            if (x == ninf) return y;
                            if (y == ninf) return x;
                            const double m = std::max(x, y);
                            return m + std::log(std::exp(x - m) + std::exp(y - m));
                        };
                        for (std::size_t i = 0; i < owners.size(); ++i) {
                            if (owners[i].kind != OwnerKind::Unusable) continue;
                            RCat& c = cats[unusable_reason_name(owners[i].why)];
                            ++c.n;
                            c.in_band = la(c.in_band, owners[i].in_band);
                            c.omitted = la(c.omitted, owners[i].omitted_bound);
                            c.unmapped = la(c.unmapped, owners[i].unmapped);
                            c.origins += owners[i].origins;
                            c.scope_blocks += owners[i].scope.size();
                            if (!owners[i].scope.empty())
                                c.max_span = std::max<std::size_t>(
                                    c.max_span,
                                    owners[i].scope.back() - owners[i].scope.front() + 1);
                            if (i < hf.size()) c.bases += hf[i].r1.size() + hf[i].r2.size();
                        }
                        const auto num = [](double v) {
                            char b[40];
                            if (v == -std::numeric_limits<double>::infinity())
                                return std::string("-inf");
                            std::snprintf(b, sizeof b, "%.6f", v); return std::string(b);
                        };
                        for (const auto& kv : cats) {
                            const std::string k = "unusable_" + kv.first + "_";
                            unusable_lines.push_back(k + "count\t" +
                                                     std::to_string(kv.second.n) + "\n");
                            unusable_lines.push_back(k + "read_bases\t" +
                                                     std::to_string(kv.second.bases) + "\n");
                            unusable_lines.push_back(k + "log_in_band\t" +
                                                     num(kv.second.in_band) + "\n");
                            unusable_lines.push_back(k + "log_omitted_bound\t" +
                                                     num(kv.second.omitted) + "\n");
                            unusable_lines.push_back(k + "log_unmapped\t" +
                                                     num(kv.second.unmapped) + "\n");
                            unusable_lines.push_back(k + "origins\t" +
                                                     std::to_string(kv.second.origins) + "\n");
                            unusable_lines.push_back(k + "scope_blocks\t" +
                                                     std::to_string(kv.second.scope_blocks) + "\n");
                            unusable_lines.push_back(k + "max_block_span\t" +
                                                     std::to_string(kv.second.max_span) + "\n");
                        }
                        // ---- THE ONLY QUESTION THAT MATTERS FOR THESE FRAGMENTS ----------
                        // Not "is the raw tail mass small" -- every read probability is tiny --
                        // but "how far can the tail move this fragment's CONTRIBUTION", which is
                        //
                        //     log( (1-eta) * lambda * (M_a + M_b) + eta * P_bg )
                        //
                        // evaluated at M = 0 and at M = the certified omitted bound. The width of
                        // that interval is what a tail can do to the log-likelihood, and it is
                        // governed by the tail RELATIVE to eta * P_bg, not by its absolute size.
                        //
                        // A per-fragment width is not the guarantee either: 3,740 individually
                        // negligible widths can sum to something that reorders genotypes. So the
                        // widths are SUMMED and the aggregate is what must fit the declared
                        // tolerance.
                        //
                        // THIS IS AN UPPER BOUND ON THE WIDTH, NOT AN ESTIMATED EFFECT.
                        // The omitted bound used here is aggregated over ALL candidates, and both
                        // homologues are given it, so every fragment is charged more tail than any
                        // single candidate pair could actually carry. A number produced this way
                        // says "no more than this"; it does not say what the biological effect is,
                        // and it must not be quoted as one. Per-candidate bounds would tighten it,
                        // and that is what adaptive deepening should use.
                        //
                        // The counts are also not a statement about the genotyper as a whole:
                        // these fragments REMAIN IN THE MARKER COUNTS and can still rank genotypes
                        // through the block unaries. "No in-band placement" is a statement about
                        // the certified linkage-placement model only.
                        {
                            const double lmix = std::log1p(-hyb_params.outlier_mix);
                            const double lbgw = std::log(hyb_params.outlier_mix);
                            const double llam = std::log(hyb_params.lambda);
                            double sum_w = 0.0, worst_w = 0.0, worst_bg = -1e308;
                            std::size_t n_audited = 0, n_infinite = 0;
                            for (std::size_t i = 0; i < owners.size() && i < hf.size(); ++i) {
                                if (owners[i].kind != OwnerKind::Unusable) continue;
                                if (owners[i].why != UnusableReason::NoInBandOrigins) continue;
                                const std::size_t len = hf[i].bases();
                                const std::size_t bee = static_cast<std::size_t>(
                                    hyb_params.bg_divergence * static_cast<double>(len));
                                const double bgf = static_cast<double>(bee) * lep +
                                                   static_cast<double>(len - bee) * l1m;
                                MassInterval m{-std::numeric_limits<double>::infinity(),
                                               owners[i].omitted_bound};
                                const MassInterval ci =
                                    fragment_contribution(m, m, false, lmix, llam, lbgw, bgf);
                                ++n_audited;
                                if (!std::isfinite(ci.upper) || !std::isfinite(ci.lower)) {
                                    ++n_infinite; continue;
                                }
                                const double w = ci.upper - ci.lower;
                                sum_w += w;
                                worst_w = std::max(worst_w, w);
                                worst_bg = std::max(worst_bg, bgf);
                            }
                            unusable_lines.push_back(
                                "audit_no_in_band_audited\t" + std::to_string(n_audited) + "\n");
                            unusable_lines.push_back(
                                "audit_no_in_band_nonfinite\t" + std::to_string(n_infinite)+"\n");
                            char nb[64];
                            std::snprintf(nb, sizeof nb, "%.6f", worst_w);
                            unusable_lines.push_back(
                                std::string("audit_worst_contribution_width_upper_bound_nats\t") + nb + "\n");
                            std::snprintf(nb, sizeof nb, "%.6f", sum_w);
                            unusable_lines.push_back(
                                std::string("audit_summed_contribution_width_upper_bound_nats\t") + nb+"\n");
                            std::snprintf(nb, sizeof nb, "%.6f", worst_bg);
                            unusable_lines.push_back(
                                std::string("audit_worst_log_p_bg\t") + nb + "\n");
                            // DECLARED, not inferred. A total of 1e-3 nats cannot reorder a call
                            // whose margins are measured in nats; anything above it is not
                            // dismissible and the fragments stay Unusable.
                            const double kGlobalToleranceNats = 1e-3;
                            std::snprintf(nb, sizeof nb, "%.6f", kGlobalToleranceNats);
                            unusable_lines.push_back(
                                std::string("audit_global_tolerance_nats\t") + nb + "\n");
                            unusable_lines.push_back(
                                "audit_bound_kind\tall-candidate aggregate, both homologues; "
                                "UPPER BOUND on the width, not an effect estimate\n");
                            // ---- WHICH DIMENSION ACTUALLY REJECTED THEM ----------------------
                            // Before deepening anything, ask whether the edit band was ever the
                            // binding constraint. The lower insert tail has already accounted for
                            // 3,746 of these; the UPPER tail predicts about 0.76 fragments beyond
                            // mean + 4 sd for this library, so a survivor is more likely a long
                            // insert than a hard read. Deepening edits would attack the wrong
                            // dimension a second time.
                            for (std::size_t i = 0; i < owners.size() && i < hf.size(); ++i) {
                                if (owners[i].kind != OwnerKind::Unusable) continue;
                                if (owners[i].why != UnusableReason::NoInBandOrigins) continue;
                                const Fragment& F = hf[i];
                                long best = -1; std::size_t e1 = 0, e2 = 0; bool both = false;
                                const std::size_t d1 =
                                    mate_band_edits(hyb_params.max_divergence, F.r1.size());
                                const std::size_t d2 =
                                    mate_band_edits(hyb_params.max_divergence, F.r2.size());
                                const std::string rc2 = F.r2.empty() ? std::string()
                                                                     : reverse_complement(F.r2);
                                for (const CandidateFrame& fr : hyb_cov.frames) {
                                    if (!fr.ok || fr.seq.empty()) continue;
                                    const auto p1 = bounded_mate_placements(F.r1, fr.seq, d1,
                                                                            nullptr, nullptr);
                                    const auto p2 = bounded_mate_placements(rc2, fr.seq, d2,
                                                                            nullptr, nullptr);
                                    if (p1.empty() || p2.empty()) continue;
                                    both = true;
                                    for (const auto& a : p1)
                                        for (const auto& b : p2) {
                                            const long ins = static_cast<long>(b.start) +
                                                             static_cast<long>(F.r2.size()) -
                                                             static_cast<long>(a.start);
                                            if (ins <= 0) continue;
                                            if (best < 0 || ins < best) {
                                                best = ins; e1 = a.edits; e2 = b.edits;
                                            }
                                        }
                                }
                                unusable_lines.push_back(
                                    "residual_fragment\t" + F.name + "\tmates_place_individually=" +
                                    (both ? "1" : "0") + "\tbest_unrestricted_insert=" +
                                    std::to_string(best) + "\tedits=" + std::to_string(e1) + "+" +
                                    std::to_string(e2) + "\tsupport=" +
                                    std::to_string(ip_lo_recorded) + "-" +
                                    std::to_string(ip_hi_recorded) + "\tabove_upper=" +
                                    ((best > ip_hi_recorded) ? "1" : "0") + "\n");
                            }
                            unusable_lines.push_back(
                                "audit_marker_occurrences\tNOT MEASURED -- these fragments remain "
                                "in the marker counts and may still rank genotypes\n");
                            unusable_lines.push_back(
                                std::string("audit_aggregate_within_tolerance\t") +
                                ((n_infinite == 0 && sum_w <= kGlobalToleranceNats) ? "1" : "0") +
                                "\n");
                        }
                        // THE ONE CATEGORY THAT MUST BLOCK.
                        unusable_lines.push_back(
                            "unusable_blocking_count\t" +
                            std::to_string(cats.count("unmapped_mass_dominant")
                                               ? cats["unmapped_mass_dominant"].n : 0) + "\n");
                    }
                    // THE FACTOR SCOPES, stated before the transaction so completeness is
                    // assessed against the model that will actually run. Without --hybrid-higher
                    // this is empty and the assessment is the pairwise one, unchanged.
                    // THE PLANNER, when no factor list was supplied. A hand-written list is fitted
                    // to the reads that produced it; this derives the factors from THIS sample's
                    // own ledger, so a donor whose fragments span blocks nobody anticipated still
                    // gets a factor that consumes them.
                    if (hybrid_higher && hybrid_factor_runs.empty()) {
                        higher_scopes = plan_higher_factors(owners, edge_status, blocks.size());
                        std::string desc;
                        for (const HigherFactorScope& f : higher_scopes) {
                            desc += " {";
                            for (std::size_t q = 0; q < f.blocks.size(); ++q)
                                desc += (q ? "," : "") + std::to_string(f.blocks[q]);
                            desc += "}";
                            if (!f.superseded.empty()) {
                                desc += "<-";
                                for (std::size_t q = 0; q < f.superseded.size(); ++q)
                                    desc += (q ? "," : "") +
                                            std::to_string(f.superseded[q].first) + "-" +
                                            std::to_string(f.superseded[q].second);
                            }
                        }
                        log.info("factor plan (derived from this sample's ledger):" +
                                 (desc.empty() ? std::string(" none") : desc));
                    }
                    // The ledger, counted where `owners` is live.
                    plan_wide_consumed.assign(higher_scopes.size(), 0);
                    for (std::size_t i = 0; i < higher_scopes.size(); ++i)
                        for (const FragmentOwner& o2 : owners)
                            if (o2.kind == OwnerKind::Wide &&
                                higher_scopes[i].covers(o2.var_scope)) ++plan_wide_consumed[i];
                    planned_only_stored = higher_scopes;
                    planned_only_reported = higher_scopes.size();
                    // THE LEDGER CREDITS ONLY FACTORS THAT EXIST.
                    //
                    // higher_scopes above is what the planner WOULD build. Handing it to the
                    // completeness assessment tells the ledger those Wide and superseded fragments
                    // have consumers, and under the planner they do not: no table is built, so
                    // nothing consumes them. The assessment must be given the scopes of factors
                    // that were actually BUILT -- which is `built_blocks` -- so an unbuilt plan
                    // shows up as unconsumed evidence and refuses, instead of as a complete model.
                    //
                    // The planner's output stays in the report as the plan, clearly separate from
                    // what was realised.
                    std::vector<HigherFactorScope> planned_only = planned_only_stored;
                    higher_scopes.clear();
                    if (!hybrid_factor_runs.empty())
                        for (std::size_t f = 0; f < built_blocks.size(); ++f) {
                            HigherFactorScope hsc;
                            hsc.blocks = built_blocks[f];
                            std::string tok;
                            const std::string spec =
                                f < built_supersede.size() ? built_supersede[f] : std::string();
                            for (char ch : spec + ",") {
                                if (ch != ',') { tok.push_back(ch); continue; }
                                const std::size_t dash = tok.find('-');
                                if (dash != std::string::npos)
                                    hsc.superseded.push_back({
                                        static_cast<std::uint32_t>(std::stoul(tok.substr(0, dash))),
                                        static_cast<std::uint32_t>(std::stoul(tok.substr(dash+1)))});
                                tok.clear();
                            }
                            higher_scopes.push_back(hsc);
                        }
                    hyb_act = plan_hybrid_activation_sparse(hf, owners, edge_status, edge_map,
                                                            maps, blocks.size(), higher_scopes);
                    // ---- THE HIGHER-ORDER TRANSACTION ------------------------------------
                    // ORDER MATTERS AND IS THE POINT. The plan is computed on the ACTUAL
                    // production chain -- every block, the real panel -- BEFORE a single read is
                    // excluded and before any message is allocated. If anything here refuses, the
                    // whole transaction is abandoned: no factors, no exclusions, and the legacy
                    // call proceeds untouched.
                    if (hybrid_higher && !built_tables.empty()) {
                        higher_hap_allele.assign(hap_names.size(),
                                                 std::vector<std::uint32_t>(blocks.size(), 0));
                        bool spanned = true;
                        for (std::size_t h = 0; h < hap_names.size() && spanned; ++h)
                            for (std::size_t b = 0; b < blocks.size(); ++b) {
                                const auto it = blocks[b].allele_of.find(hap_names[h]);
                                if (it == blocks[b].allele_of.end()) { spanned = false; break; }
                                higher_hap_allele[h][b] =
                                    static_cast<std::uint32_t>(it->second);
                            }
                        // EVERY superseded edge, from every factor's own list.
                        std::set<std::pair<std::uint32_t, std::uint32_t>> superseded;
                        for (const std::string& spec : built_supersede) {
                            std::string tok;
                            for (char ch : spec + ",") {
                                if (ch != ',') { tok.push_back(ch); continue; }
                                const std::size_t dash = tok.find('-');
                                if (dash != std::string::npos)
                                    superseded.insert({
                                        static_cast<std::uint32_t>(std::stoul(tok.substr(0, dash))),
                                        static_cast<std::uint32_t>(std::stoul(tok.substr(dash+1)))});
                                tok.clear();
                            }
                        }
                        // ONE LIST. Retained pairwise edges first, then the interval factors --
                        // a superseded edge is never added, so it cannot be applied at all.
                        std::size_t retained_pairwise = 0, dropped_superseded = 0,
                                    neutral_pairwise = 0;
                        for (std::size_t b = 1; b < hyb_act.sparse_kernel_edges.size(); ++b) {
                            const SparseEdgeLinkage& e = hyb_act.sparse_kernel_edges[b];
                            // ACTIVE IS NOT THE SAME AS CARRYING A CORRECTION. An edge that built
                            // but stored no non-neutral phase class contributes exactly 1 and is
                            // not added -- correct, but it means active_edges is NOT evidence that
                            // that many pairwise factors reached inference. Counted separately so
                            // the two numbers can never be read as the same claim.
                            if (e.active && e.classes.empty()) ++neutral_pairwise;
                            if (!e.has_corrections()) continue;
                            const std::pair<std::uint32_t, std::uint32_t> key{
                                static_cast<std::uint32_t>(b - 1), static_cast<std::uint32_t>(b)};
                            if (superseded.count(key)) { ++dropped_superseded; continue; }
                            HybridHigherFactor HF;
                            HF.blocks = {key.first, key.second};
                            HF.pairwise = &hyb_act.sparse_kernel_edges[b];
                            higher_factors.push_back(HF);
                            ++retained_pairwise;
                        }
                        for (std::size_t f = 0; f < built_tables.size(); ++f) {
                            HybridHigherFactor HF;
                            HF.blocks = built_blocks[f];
                            HF.table = &built_tables[f];
                            higher_factors.push_back(HF);
                        }
                        HybridChain PC;
                        PC.n_hap = hap_names.size(); PC.n_blocks = blocks.size();
                        PC.recomb = 0.5;   // the plan does not depend on r beyond r in (0,1)
                        PC.edges.assign(blocks.size(), HybridEdge{});
                        PC.log_emission.assign(blocks.size(),
                                               std::vector<double>(PC.n_hap * PC.n_hap, 0.0));
                        PC.hap_allele = higher_hap_allele;
                        PC.higher = higher_factors;
                        higher_plan = spanned ? plan_higher_order(PC, 1) : HigherOrderPlan{};
                        if (!spanned) higher_plan.refusal =
                            "not every panel haplotype spans every block";
                        const bool fits = higher_plan.ok &&
                            (hybrid_max_message_entries == 0 ||
                             higher_plan.peak_message_entries <= hybrid_max_message_entries);
                        if (!fits) {
                            higher_refusal = higher_plan.ok
                                ? ("planned peak message of " +
                                   std::to_string(higher_plan.peak_message_entries) +
                                   " entries exceeds --hybrid-max-message-entries")
                                : higher_plan.refusal;
                            higher_factors.clear(); higher_hap_allele.clear();
                            log.info("hybrid higher-order: REFUSED before allocation (" +
                                     higher_refusal + "); nothing excluded, legacy call stands");
                        } else if (!hyb_act.hybrid_activated) {
                            // THE TRANSACTION IS ONE TRANSACTION. If the pairwise side did not
                            // activate, the model is INCOMPLETE and the higher factors must not
                            // run either -- a partly-higher-order call is not a defined object.
                            higher_refusal = "the hybrid transaction did not activate (" +
                                             hyb_act.refusal + ")";
                            higher_factors.clear(); higher_hap_allele.clear();
                            log.info("hybrid higher-order: STOOD DOWN because the transaction is "
                                     "incomplete; nothing excluded, legacy call stands");
                        } else {
                            higher_active = true;
                            neutral_pairwise_reported = neutral_pairwise;
                            superseded_reported = dropped_superseded;
                            log.info("hybrid higher-order: " +
                                     std::to_string(retained_pairwise) + " retained pairwise (" +
                                     std::to_string(neutral_pairwise) + " active but exactly "
                                     "neutral, contributing nothing) + " +
                                     std::to_string(built_tables.size()) + " interval factors, " +
                                     std::to_string(dropped_superseded) + " superseded edges "
                                     "excluded; planned peak " +
                                     std::to_string(higher_plan.peak_message_entries) +
                                     " entries, " +
                                     std::to_string(higher_plan.total_bytes / 1048576) + " MB, " +
                                     std::to_string(higher_plan.forward_updates) +
                                     " forward updates");
                        }
                    }
                    // EXCLUSIONS ONLY NOW, AND ONLY FROM ACTIVE CONSUMERS.
                    //
                    // `built_tables.empty()` used to appear in this condition as a permissive
                    // case, and with the planner it was ALWAYS true: the planner supplies scopes
                    // but no CLI factor runs, so no table is ever built, higher_order_active stays
                    // 0 -- and the ledger still counted the Wide and superseded fragments as
                    // consumed. Exclusions were then applied on the strength of consumers that do
                    // not exist. On one pilot donor that removed 316 fragments from the marker
                    // counts while zero interval factors reached inference, and the run still
                    // reported COMPLETE. That is evidence loss, not an incomplete call.
                    //
                    // THE INVARIANT: a plan with factors in it must have produced factors that
                    // reached inference, or nothing is excluded and the call is INCOMPLETE.
                    const bool plan_is_honoured =
                        planned_only.empty() ? true : higher_active;
                    if (hybrid_higher && !plan_is_honoured && higher_refusal.empty()) {
                        higher_refusal = "the plan plans " + std::to_string(planned_only.size()) +
                                         " factor(s) but none reached inference";
                    }
                    if (hyb_act.hybrid_activated && (!hybrid_higher || plan_is_honoured)) {
                        hyb_exclusions.insert(hyb_act.excluded_fragments.begin(),
                                              hyb_act.excluded_fragments.end());
                    } else if (hyb_act.hybrid_activated) {
                        hyb_act.hybrid_activated = false;
                        hyb_act.call_status = HybridCallStatus::Incomplete;
                        hyb_act.refusal = "higher-order: " + higher_refusal;
                    }
                    }   // end of the verified-window budget else
                    }   // end of the overflow guard else
                }
                // PER-EDGE MEASUREMENT, serialised from the objects the caller actually built.
                if (!hybrid_edges_path.empty()) {
                    std::ofstream ef(hybrid_edges_path);
                    if (!ef) throw std::runtime_error("genotype: cannot write " + hybrid_edges_path);
                    ef.precision(6);
                    ef << "block_a\tblock_b\tn_alleles_a\tn_alleles_b\towned_fragments"
                          "\tinformative_fragments\tsupport_cells\tstored_classes"
                          "\tpredicted_classes\ttheoretical_configs\tsparse_bytes"
                          "\tbuild_seconds\tstatus\tseed_hits\tproposed_states"
                          "\tverified_windows\tdense_windows\tfallbacks"
                          "\tread_verifications\tfinite_cells\treduction\n";
                    std::size_t tot_stored = 0, tot_bytes = 0, tot_support = 0, usable = 0;
                    double tot_build = 0.0;
                    for (const EdgeRow& r : edge_rows) {
                        ef << r.a << '\t' << r.b << '\t' << r.na << '\t' << r.nb << '\t'
                           << r.owned << '\t' << r.informative << '\t' << r.support << '\t'
                           << r.stored << '\t' << r.predicted << '\t' << r.theoretical << '\t'
                           << r.bytes << '\t' << r.build_s << '\t' << r.status << '\t'
                           << r.seed_hits << '\t' << r.proposed << '\t' << r.verified << '\t'
                           << r.dense_windows << '\t' << r.fallbacks << '\t'
                           << r.read_verifications << '\t' << r.finite_cells << '\t'
                           << (r.verified ? static_cast<double>(r.dense_windows) / r.verified : 0.0)
                           << '\n';
                        tot_stored += r.stored; tot_bytes += r.bytes; tot_support += r.support;
                        tot_build += r.build_s;
                        if (r.status == std::string("ok")) ++usable;
                    }
                    // AGGREGATE, because production holds every edge at once. A sum of isolated
                    // per-edge estimates is not what the process actually needs.
                    struct rusage ru {};
                    double peak_mb = 0.0;
                    if (getrusage(RUSAGE_SELF, &ru) == 0) {
#ifdef __APPLE__
                        peak_mb = static_cast<double>(ru.ru_maxrss) / 1048576.0;   // bytes
#else
                        peak_mb = static_cast<double>(ru.ru_maxrss) / 1024.0;      // kilobytes
#endif
                    }
                    ef << "#total_edges\t" << edge_rows.size() << '\n';
                    ef << "#usable_edges\t" << usable << '\n';
                    ef << "#refused_edges\t" << (edge_rows.size() - usable) << '\n';
                    ef << "#total_support_cells\t" << tot_support << '\n';
                    ef << "#total_stored_classes\t" << tot_stored << '\n';
                    ef << "#total_sparse_bytes\t" << tot_bytes << '\n';
                    std::size_t tot_seed = 0, tot_prop = 0, tot_ver = 0, tot_dense = 0,
                                tot_fb = 0;
                    for (const EdgeRow& r : edge_rows) {
                        tot_seed += r.seed_hits; tot_prop += r.proposed; tot_ver += r.verified;
                        tot_dense += r.dense_windows; tot_fb += r.fallbacks;
                    }
                    ef << "#total_build_seconds\t" << tot_build << '\n';
                    ef << "#total_seed_hits\t" << tot_seed << '\n';
                    ef << "#total_proposed_states\t" << tot_prop << '\n';
                    ef << "#total_verified_windows\t" << tot_ver << '\n';
                    ef << "#total_dense_windows\t" << tot_dense << '\n';
                    ef << "#total_exhaustive_fallbacks\t" << tot_fb << '\n';
                    ef << "#overall_reduction\t"
                       << (tot_ver ? static_cast<double>(tot_dense) / tot_ver : 0.0) << '\n';
                    ef << "#peak_rss_mb\t" << peak_mb << '\n';
                    ef << "#hybrid_status\t" << hybrid_call_status_name(hyb_act.call_status) << '\n';
                    ef << "#active_edges\t" << hyb_act.active_edges << '\n';
                    ef << "#fragments_excluded\t" << hyb_exclusions.size() << '\n';
                    ef << "#max_classes_limit\t" << hybrid_max_classes << '\n';
                    ef << "#max_bytes_limit\t" << hybrid_max_bytes << '\n';
                    ef.flush();
                    if (!ef) throw std::runtime_error("genotype: write failed for " +
                                                      hybrid_edges_path);
                    log.wrote({hybrid_edges_path});
                }
                // PROVENANCE. Every parameter the linkage factors used, written with the result,
                // so a real-data number is interpretable later without reading the source. lambda is
                // marked supplied or default because "which lambda" is the first question anyone
                // will ask of a C4 result.
                if (!hybrid_status_path.empty()) {
                    std::ofstream hs(hybrid_status_path);
                    if (!hs) throw std::runtime_error("genotype: cannot write " + hybrid_status_path);
                    hs.precision(10);
                    hs << "field\tvalue\n";
                    hs << "hybrid_status\t" << hybrid_call_status_name(hyb_act.call_status) << '\n';
                    // PROVENANCE: which fragment universe this run scored. Two runs under
                    // different insert floors are not comparable, and the difference is invisible
                    // in every other field.
                    hs << "insert_floor_policy\t"
                       << (hyb_params.allow_overlapping_pairs
                               ? "max(|r1|,|r2|) -- overlapping pairs are valid"
                               : "|r1|+|r2| -- COMPATIBILITY, overlapping pairs excluded")
                       << '\n';
                    hs << "insert_floor_bp\t" << min_len_recorded << '\n';
                    hs << "insert_support\t" << ip_lo_recorded << "-" << ip_hi_recorded << '\n';
                    hs << "insert_sigmas\t" << hyb_params.insert_sigmas << '\n';
                    // NOT CALLED ZERO. Any finite support leaves Gaussian tail mass outside it;
                    // reporting it is what lets a later caller certify the residual through the
                    // same contribution-width machinery instead of assuming it away.
                    hs << "insert_residual_log_below\t" << ip_residual_lo << '\n';
                    hs << "insert_residual_log_above\t" << ip_residual_hi << '\n';
                    hs << "ownership_complete\t" << (hyb_act.ownership_complete ? 1 : 0) << '\n';
                    hs << "factors_buildable\t" << (hyb_act.factors_buildable ? 1 : 0) << '\n';
                    hs << "hybrid_activated\t" << (hyb_act.hybrid_activated ? 1 : 0) << '\n';
                    hs << "legacy_call_status\tAVAILABLE\n";
                    hs << "reason\t" << (hyb_act.refusal.empty() ? "-" : hyb_act.refusal) << '\n';
                    hs << "all_states_usable\t" << (hyb_cov.all_states_usable ? 1 : 0) << '\n';
                    hs << "fragments_loaded\t" << hyb_fragments_loaded << '\n';
                    hs << "candidate_edges\t" << hyb_edges_considered << '\n';
                    hs << "active_edges\t" << hyb_act.active_edges << '\n';
                    hs << "fragments_excluded\t" << hyb_exclusions.size() << '\n';
                    // ---- THE HIGHER-ORDER ACCOUNTING ---------------------------------------
                    // Every claim a C4 acceptance would rest on, emitted as a number rather than
                    // left to be inferred from the log.
                    hs << "higher_order_requested\t" << (hybrid_higher ? 1 : 0) << '\n';
                    // THE PLAN AND ITS EVIDENCE LEDGER. Which factors the rule derived, what each
                    // one consumes, and what it supersedes -- so "supersession follows evidence"
                    // is a checkable statement rather than a described intention.
                    hs << "factor_plan_source\t"
                       << (hybrid_factor_runs.empty() ? "planner" : "explicit --hybrid-factor-run")
                       << '\n';
                    hs << "factor_plan_size\t" << planned_only_reported << '\n';
                    hs << "factor_scopes_credited_to_ledger\t" << higher_scopes.size() << '\n';
                    // THE INVARIANT, REPORTED. planned == built == reached inference, or the call
                    // is not a higher-order call at all.
                    hs << "factor_plan_honoured\t"
                       << ((higher_scopes.empty() || higher_active) ? 1 : 0) << '\n';
                    for (const HigherFactorScope& f : planned_only_stored) {
                        hs << "factor_plan\t";
                        for (std::size_t q = 0; q < f.blocks.size(); ++q)
                            hs << (q ? "," : "") << f.blocks[q];
                        hs << '\t';
                        const std::size_t pi =
                            static_cast<std::size_t>(&f - planned_only_stored.data());
                        hs << (pi < plan_wide_consumed.size() ? plan_wide_consumed[pi] : 0)
                           << '\t';
                        if (f.superseded.empty()) hs << '-';
                        for (std::size_t q = 0; q < f.superseded.size(); ++q)
                            hs << (q ? "," : "") << f.superseded[q].first << '-'
                               << f.superseded[q].second;
                        hs << '\n';
                    }
                    hs << "higher_order_active\t" << (higher_active ? 1 : 0) << '\n';
                    hs << "higher_order_refusal\t"
                       << (higher_refusal.empty() ? "-" : higher_refusal) << '\n';
                    if (hybrid_higher) {
                        std::size_t n_pw = 0, n_iv = 0;
                        std::set<std::string> factor_keys;
                        bool dup = false;
                        for (const HybridHigherFactor& f : higher_factors) {
                            std::string k;
                            for (std::uint32_t q : f.blocks) k += std::to_string(q) + ".";
                            k += f.table != nullptr ? "I" : "P";
                            if (!factor_keys.insert(k).second) dup = true;
                            if (f.table != nullptr) ++n_iv; else ++n_pw;
                        }
                        hs << "higher_retained_pairwise_factors\t" << n_pw << '\n';
                        hs << "higher_neutral_pairwise_edges\t" << neutral_pairwise_reported
                           << '\n';
                        hs << "higher_pairwise_accounted\t"
                           << ((n_pw + neutral_pairwise_reported + superseded_reported ==
                                hyb_act.active_edges + superseded_reported) ? 1 : 0) << '\n';
                        hs << "higher_interval_factors\t" << n_iv << '\n';
                        hs << "higher_each_factor_once\t" << (dup ? 0 : 1) << '\n';
                        hs << "higher_ls_transitions\t"
                           << (blocks.empty() ? 0 : blocks.size() - 1) << '\n';
                        for (const HybridHigherFactor& f : higher_factors) {
                            hs << "higher_factor\t";
                            for (std::size_t j = 0; j < f.blocks.size(); ++j)
                                hs << (j ? "," : "") << f.blocks[j];
                            hs << (f.table != nullptr ? "\tinterval\t" : "\tpairwise\t")
                               << (f.table != nullptr ? f.table->classes_stored
                                                      : f.pairwise->classes.size()) << '\n';
                        }
                        hs << "higher_plan_ok\t" << (higher_plan.ok ? 1 : 0) << '\n';
                        hs << "higher_plan_peak_message_entries\t"
                           << higher_plan.peak_message_entries << '\n';
                        hs << "higher_plan_forward_updates\t"
                           << higher_plan.forward_updates << '\n';
                        hs << "higher_plan_total_bytes\t" << higher_plan.total_bytes << '\n';
                        hs << "higher_planned_before_exclusions\t1\n";
                        // The REALISED counts are not known here -- this file is written before
                        // the genotyper runs. They go to --hybrid-higher-report, afterwards.
                        hs << "higher_actuals\tsee --hybrid-higher-report\n";
                        if (false) {
                            hs << "higher_actual_forward_updates\t"
                               << higher_stats.forward_updates << '\n';
                            hs << "higher_actual_adjoint_updates\t"
                               << higher_stats.adjoint_updates << '\n';
                            hs << "higher_actual_peak_message_entries\t"
                               << higher_stats.peak_message_entries << '\n';
                            hs << "higher_all_initial_emissions_finite\t"
                               << (higher_stats.all_initial_emissions_finite ? 1 : 0) << '\n';
                            hs << "higher_initial_states_dropped\t"
                               << higher_stats.initial_states_dropped << '\n';
                            // EXACT when every initial emission is finite; a BOUND otherwise. Both
                            // are checked, and which one applied is reported.
                            const bool exact = higher_stats.all_initial_emissions_finite;
                            const bool holds = exact
                                ? (higher_stats.forward_updates == higher_plan.forward_updates &&
                                   higher_stats.peak_message_entries ==
                                       higher_plan.peak_message_entries)
                                : (higher_stats.forward_updates <= higher_plan.forward_updates &&
                                   higher_stats.peak_message_entries <=
                                       higher_plan.peak_message_entries);
                            hs << "higher_plan_contract\t" << (exact ? "EQUAL" : "BOUND") << '\n';
                            hs << "higher_plan_contract_holds\t" << (holds ? 1 : 0) << '\n';
                            hs << "higher_forward_seconds\t"
                               << higher_stats.forward_seconds << '\n';
                            hs << "higher_adjoint_seconds\t"
                               << higher_stats.adjoint_seconds << '\n';
                        }
                    }
                    for (const std::string& ln : unusable_lines) hs << ln;
                    hs << "consumed_higher_wide\t"
                       << hyb_act.report.consumed_higher_wide << '\n';
                    hs << "consumed_higher_superseded\t"
                       << hyb_act.report.consumed_higher_superseded << '\n';
                    hs << "consumed_fragments\t" << hyb_act.consumed_fragments << '\n';
                    hs << "excluded_equals_consumed\t"
                       << ((hyb_act.excluded_fragments.size() == hyb_act.consumed_fragments &&
                            hyb_exclusions.size() == hyb_act.consumed_fragments) ? 1 : 0) << '\n';
                    hs << "owned_unary\t" << hyb_act.report.consumed_unary << '\n';
                    hs << "owned_linkage\t" << hyb_act.report.consumed_linkage << '\n';
                    hs << "owned_invariant\t" << hyb_act.report.invariant << '\n';
                    hs << "unconsumed_wide\t" << hyb_act.report.unconsumed_wide << '\n';
                    // THE SCOPES THEMSELVES. "1 wide" says a fragment has no consumer; it does not
                    // say which blocks it depends on, which is the only thing that tells you
                    // whether a factor set covers a cohort or was fitted to one donor's reads.
                    for (const std::vector<std::uint32_t>& sc : hyb_act.report.wide_scopes) {
                        hs << "unconsumed_wide_scope\t";
                        for (std::size_t q = 0; q < sc.size(); ++q)
                            hs << (q ? "," : "") << sc[q];
                        hs << '\n';
                    }
                    hs << "unconsumed_refused_edge\t" << hyb_act.report.unconsumed_refused_edge << '\n';
                    hs << "unconsumed_unusable\t" << hyb_act.report.unconsumed_unusable << '\n';
                    for (const EdgeRefusal& r : hyb_act.report.refusals) {
                        hs << "edge_refusal\t" << r.block_a << '-' << r.block_b << ':'
                           << linkage_status_name(r.status) << '(' << r.n_fragments << ')'
                           << (r.detail.empty() ? "" : " " + r.detail) << '\n';
                    }
                    // THE PARAMETER CONTRACT. These do not all come from one place, so each is named
                    // for what it is rather than grouped as "the model".
                    // OPERATIONAL WORK ACTUALLY DONE, beside the limits it was held to.
                    hs << "work_proposed_cells\t" << hyb_work_proposed << '\n';
                    hs << "work_full_read_verifications\t" << hyb_work_verifications << '\n';
                    hs << "work_bases_compared_upper_bound\t" << hyb_work_bases << '\n';
                    hs << "limit_proposed_cells\t" << hybrid_max_proposed_cells << '\n';
                    hs << "limit_full_read_verifications\t"
                       << hybrid_max_full_read_verifications << '\n';
                    hs << "param_lambda\t" << hyb_params.lambda << '\n';
                    hs << "param_lambda_source\t"
                       << (hyb_params.lambda_source ==
                               HybridLinkageParameters::LambdaSource::Supplied  ? "supplied"
                         : hyb_params.lambda_source ==
                               HybridLinkageParameters::LambdaSource::Estimated ? "estimated"
                                                                                : "default")
                       << '\n';
                    hs << "param_lambda_estimated\t" << hyb_params.lambda_estimated << '\n';
                    hs << "param_median_panel_length\t" << hyb_params.median_panel_length << '\n';
                    hs << "param_n_fragments\t" << hyb_params.n_fragments << '\n';
                    hs << "param_lambda_ratio_used_over_estimated\t"
                       << (hyb_params.lambda_estimated > 0.0
                               ? hyb_params.lambda / hyb_params.lambda_estimated : 0.0) << '\n';
                    hs << "param_bg_divergence\t" << hyb_params.bg_divergence << '\n';
                    hs << "param_outlier_mix\t" << hyb_params.outlier_mix << '\n';
                    hs << "param_error_rate\t" << hyb_params.error_rate << '\n';
                    hs << "param_max_divergence\t" << hyb_params.max_divergence << '\n';
                    hs << "param_fragment_len\t" << hyb_params.fragment_len << '\n';
                    hs << "param_fragment_sd\t" << hyb_params.fragment_sd << '\n';
                    hs << "param_discordant_rate\t" << hyb_params.discordant_rate << '\n';
                    hs << "param_insert_sigmas\t" << hyb_params.insert_sigmas << '\n';
                    hs.flush();
                    if (!hs) throw std::runtime_error("genotype: write failed for " +
                                                      hybrid_status_path);
                    log.wrote({hybrid_status_path});
                }
                log.info(std::string("hybrid: ") +
                         hybrid_call_status_name(hyb_act.call_status) + ", " +
                         std::to_string(hyb_fragments_loaded) + " fragments, " +
                         std::to_string(hyb_edges_considered) + " candidate edge(s), " +
                         std::to_string(hyb_act.active_edges) + " active, " +
                         std::to_string(hyb_exclusions.size()) + " excluded from markers" +
                         (hyb_act.refusal.empty() ? "" : "; " + hyb_act.refusal));
            }

            // A DRY RUN stops after transactional planning, having serialised the exact edge
            // objects and counters normal inference would consume -- not a second code path that
            // rebuilds them.
            if (hybrid_call && hybrid_dry_run) {
                log.info("hybrid dry run: planning complete, stopping before inference");
                log.done();
                return 0;
            }

            // MARKER OCCURRENCE EXCLUSION. Fragments owned by a linkage edge contribute their
            // sequence there and must leave the marker counts, or the same read is counted twice.
            // Only their OCCURRENCES are subtracted: a marker they share with a unary-owned
            // fragment keeps that fragment's counts and stays in the panel.
            std::unordered_set<std::string> excl;
            if (!exclude_fragments_path.empty()) {
                std::ifstream ef(exclude_fragments_path);
                if (!ef) throw std::runtime_error("genotype: cannot read " + exclude_fragments_path);
                std::string line;
                while (std::getline(ef, line)) {
                    while (!line.empty() && (line.back() == '\r' || line.back() == '\n')) line.pop_back();
                    if (!line.empty()) excl.insert(line);
                }
            }
            // The hybrid's exclusions are the COMMITTED ones: empty unless the transaction
            // activated, so a refusal cannot subtract anything.
            for (const std::string& nm : hyb_exclusions) excl.insert(nm);
            std::size_t excluded_reads = 0;
            ReadCounts rc = count_reads(read_paths, read_panel, options.threads,
                                        excl.empty() ? nullptr : &excl, &excluded_reads);
            if (!excl.empty()) {
                log.info("marker exclusion: " + std::to_string(excl.size()) +
                         " fragment(s) named, " + std::to_string(excluded_reads) +
                         " reads left the marker counts (their markers remain in the panel)");
            }
            log.info("reads: " + std::to_string(rc.reads) + " (" + std::to_string(rc.bases / 1000) +
                     " kb); " + std::to_string(rc.syncmers) + " syncmers, " +
                     std::to_string(100 * rc.matched_syncmers / std::max<std::uint64_t>(1, rc.syncmers)) +
                     "% matched a panel marker; " + std::to_string(rc.novel_adjacencies) +
                     " novel adjacencies");
            std::size_t region_bp_hint = 0;
            for (const PathRecord& p : panel_graph.paths) {
                if (p.name == reference_path) {
                    region_bp_hint = spell_path_steps_sequence(panel_graph, p.steps).size();
                    break;
                }
            }
            // The joint refinement can only change anything when some haplotype bypasses a block --
            // otherwise every block is traversed twice, the traversal count is a constant 2, and the
            // second pass reproduces the first. Gate on that, so a panel without bypass alleles pays
            // nothing for it.
            bool panel_has_bypass = false;
            for (const BlockAlleles& b : blocks) {
                if (b.bypass_allele >= 0) { panel_has_bypass = true; break; }
            }
            const bool run_joint = depth_model == DepthModel::Joint && panel_has_bypass;
            DepthRegionStats region_stats;
            std::vector<BlockDepth> depth =
                estimate_depth(read_panel, rc, min_anchors, uneven_tolerance,
                               // Joint refines a first pass, so its starting point matters: seeded
                               // with Median it reproduces Median's own homozygous answer and never
                               // escapes. Bases does not depend on the genotype at all, which is
                               // exactly what a starting point needs.
                               run_joint ? (region_bp_hint > 0 ? DepthModel::Bases : DepthModel::Quantile)
                                         : (depth_model == DepthModel::Joint ? DepthModel::Median
                                                                             : depth_model),
                               depth_quantile, region_bp_hint, depth_estimator, &region_stats);
            std::size_t uneven = 0;
            std::vector<double> lam;
            for (const BlockDepth& d : depth) { if (d.uneven) ++uneven; if (d.usable) lam.push_back(d.lambda_hap); }
            std::sort(lam.begin(), lam.end());
            log.info("per-haplotype depth (lambda): median " +
                     std::to_string(lam.empty() ? 0.0 : lam[lam.size() / 2]) + " over " +
                     std::to_string(lam.size()) + " usable blocks; " + std::to_string(uneven) +
                     " flagged UNEVEN");
            // All three side by side, because the choice is not neutral and the median cannot express
            // a value between half-integers however many anchors are pooled.
            // Deliberately NOT reported as lambda. Under Bases, Quantile or the final Joint pass the
            // fitted depth comes from elsewhere entirely, so printing this as the lambda in force
            // would misstate what the model used -- the same defect the raw-versus-fitted split above
            // exists to prevent, one level up.
            log.info("region anchors: " + std::to_string(region_stats.n_anchor) + "; median " +
                     std::to_string(region_stats.median) + ", mean " + std::to_string(region_stats.mean) +
                     ", trimmed " + std::to_string(region_stats.trimmed_mean) +
                     "; selected_anchor_center " + std::to_string(region_stats.used));
            log_depth_provenance(log, chain, depth);
            // The audit is written AFTER the joint pass below, not here: joint replaces every block's
            // fitted depth, so an audit written at this point would describe a state the emission
            // never used. `raw_anchor_*` are untouched by any model, so nothing is lost by waiting.

            // hap_names is built earlier now: the hybrid transaction has to be planned BEFORE
            // marker counting, because its exclusion set is an input to counting.
            // ---- HYBRID PREFLIGHT ---------------------------------------------------------
            // THE REQUIREMENT IS OVER THE HMM'S DECLARED STATE UNIVERSE, not the raw panel paths.
            // A recorded state reduction is legitimate; silently dropping candidates because their
            // frame construction failed is not -- it would make the linkage topology and the
            // posterior depend on an undocumented change of state space.
            //
            // Reported before any wiring so a disagreement with the earlier frame-coverage results
            // (131/131 on C4, 127/127 on CYP2D6) surfaces as a reconciliation problem rather than
            // as an inference difference nobody can attribute.
            if (!hybrid_preflight.empty()) {
                // ONE structured result. The writer below only SERIALISES it; hybrid activation
                // will consume the same object rather than recomputing coverage, so the report and
                // the decision cannot describe different runs.
                const FrameCoverage cov = assess_frame_coverage(graph, blocks, hap_names);
                std::ofstream pf(hybrid_preflight);
                if (!pf) throw std::runtime_error("genotype: cannot write " + hybrid_preflight);
                pf << "metric\tvalue\n";
                pf << "raw_panel_paths\t" << panel_graph.paths.size() << '\n';
                pf << "hmm_states\t" << hap_names.size() << '\n';
                pf << "framed_states\t" << cov.framed_names.size() << '\n';
                pf << "complete_frames\t" << cov.complete_names.size() << '\n';
                // ACCEPTED, not missing: a path ending inside a block gives a correct prefix whose
                // bytes agree, with only the remainder unmapped. Reported so the case stays visible.
                pf << "accepted_partial_frames\t" << cov.partial_names.size() << '\n';
                pf << "missing_states\t" << cov.missing_names.size() << '\n';
                pf << "names_unique\t" << (cov.names_unique ? 1 : 0) << '\n';
                pf << "all_states_usable\t" << (cov.all_states_usable ? 1 : 0) << '\n';
                std::vector<std::string> states = hap_names;
                std::sort(states.begin(), states.end());
                for (const std::string& nm : states) pf << "hmm_state\t" << nm << '\n';
                for (const std::string& nm : cov.partial_names) pf << "partial\t" << nm << '\n';
                for (std::size_t i2 = 0; i2 < cov.missing_names.size(); ++i2) {
                    pf << "missing\t" << cov.missing_names[i2] << '\t'
                       << cov.missing_reasons[i2] << '\n';
                }
                pf.flush();
                if (!pf) throw std::runtime_error("genotype: write failed for " + hybrid_preflight);
                log.info("hybrid preflight: " + std::to_string(hap_names.size()) +
                         " HMM states, " + std::to_string(cov.complete_names.size()) +
                         " complete + " + std::to_string(cov.partial_names.size()) +
                         " accepted partial, " + std::to_string(cov.missing_names.size()) +
                         " missing; coverage " +
                         (cov.all_states_usable ? "all states usable" : "NOT all states usable"));
                if (!cov.all_states_usable) {
                    // A SUBSTRATE result, not a biological one, and named so it cannot be read as
                    // genotype ambiguity.
                    log.info("hybrid_status INCOMPLETE / reason candidate-frame-coverage / "
                             "legacy_call_status AVAILABLE / hybrid_call NA");
                }
                log.wrote({hybrid_preflight});
            }

            GenotypeOptions gopt = make_genotype_options();
            gopt.scale_weight = scale_weight;
            GenotypeSummary gsum;
            std::vector<int> ta1;
            std::vector<int> ta2;
            // An EMPTY truth sequence is a real genotype -- the deletion/bypass allele -- so the
            // empty string cannot double as "no truth here". Conflating them made the evaluation
            // silently skip blocks where both haplotypes carry the deletion, which is exactly the
            // case a deletion caller has to get right.
            std::vector<char> tav1, tav2;   // does this haplotype traverse the block at all?
            std::vector<std::string> ts1;   // true spelled sequence per block, even when unrepresentable
            std::vector<std::string> ts2;
            if (!truth_haplotypes.empty()) {
                const auto comma = truth_haplotypes.find(',');
                if (comma == std::string::npos) {
                    throw std::runtime_error("genotype: --truth-haplotypes needs two comma-separated names");
                }
                const std::string a = truth_haplotypes.substr(0, comma);
                const std::string b = truth_haplotypes.substr(comma + 1);

                // Whole-locus marker multiplicity for the truth. Every other truth measurement here
                // is scoped to one block, which is exactly what cannot answer "does the block
                // decomposition account for everywhere this marker occurs". Spelled from the truth's
                // own full walk and counted with the panel's own syncmer code, so the two cannot
                // drift apart.
                if (!truth_markers_out.empty()) {
                    std::unordered_map<std::uint64_t, std::uint32_t> tot;
                    // Where each occurrence sits, per haplotype. Counts alone cannot test whether
                    // lambda * multiplicity is the right observation model: an occurrence near a
                    // sequence end has fewer possible fragment starts than one in the middle, so
                    // equal multiplicity does not mean equal expected coverage.
                    std::unordered_map<std::uint64_t, std::vector<std::pair<int, std::size_t>>> pos;
                    std::vector<std::size_t> hap_len;
                    std::size_t nfound = 0;
                    for (const std::string& nm : {a, b}) {
                        const PathRecord* pr = nullptr;
                        for (const PathRecord& p : graph.paths) if (p.name == nm) pr = &p;
                        if (pr == nullptr) continue;
                        ++nfound;
                        const std::string seq = spell_path_steps_sequence(graph, pr->steps);
                        const std::vector<KmerOccurrence> sy =
                            read_panel.all_kmers
                                ? collect_canonical_kmer_occurrences(seq, read_panel.kmer_size)
                                : collect_syncmers(seq, read_panel.kmer_size, read_panel.syncmer_s);
                        for (const KmerOccurrence& o : sy) {
                            ++tot[o.code];
                            pos[o.code].emplace_back(static_cast<int>(nfound - 1), o.start);
                        }
                        hap_len.push_back(seq.size());
                    }
                    if (nfound != 2) {
                        throw std::runtime_error("genotype: --dump-truth-marker-counts could not find "
                                                 "both truth haplotypes in the graph");
                    }
                    std::ofstream tf(truth_markers_out);
                    if (!tf) throw std::runtime_error("genotype: cannot write " + truth_markers_out);
                    tf << "slot\tcode\ttruth_multiplicity_whole_locus\toccurrences\n";
                    for (std::size_t sl = 0; sl < read_panel.node_codes.size(); ++sl) {
                        const auto it = tot.find(read_panel.node_codes[sl]);
                        tf << sl << '\t' << read_panel.node_codes[sl] << '\t'
                           << (it == tot.end() ? 0u : it->second) << '\t';
                        const auto pit = pos.find(read_panel.node_codes[sl]);
                        if (pit == pos.end()) { tf << ".\n"; continue; }
                        bool first = true;
                        for (const auto& [h, st] : pit->second) {
                            if (!first) tf << ',';
                            tf << h << ':' << st;
                            first = false;
                        }
                        tf << '\n';
                    }
                    tf << "# haplotype_lengths";
                    for (const std::size_t L : hap_len) tf << '\t' << L;
                    tf << '\n';
                    tf.flush();
                    if (!tf) throw std::runtime_error("genotype: write failed for " + truth_markers_out);
                    log.wrote({truth_markers_out});
                }
                ta1.assign(chain.size(), -1);
                ta2.assign(chain.size(), -1);
                // Under leave-one-out the truth haplotypes are not in the panel, so their allele is
                // resolved by spelling their own walk and matching it against the reduced panel's
                // allele sequences. No match means the panel simply cannot represent that allele --
                // the mosaic ceiling -- and the block is scored as unrepresentable rather than wrong.
                // Alleles of the held-out haplotypes, enumerated exactly as the panel's were.
                std::vector<BlockAlleles> held_blocks(chain.size());
                if (!held_out.empty()) {
                    Graph held_graph = graph;
                    held_graph.paths = held_out;
                    std::vector<BubblePathIndex> held_idx(held_graph.paths.size());
                    for (std::size_t k = 0; k < held_graph.paths.size(); ++k) {
                        held_idx[k] = build_bubble_path_index(held_graph.paths[k]);
                    }
                    for (std::size_t bi = 0; bi < chain.size(); ++bi) {
                        held_blocks[bi] = enumerate_block_alleles(held_graph, held_idx, bubbles,
                                                                 chain[bi], options.threads);
                    }
                }
                auto resolve = [&](const std::string& name, std::vector<int>& out_alleles,
                                   std::vector<std::string>& out_seq, std::vector<char>& out_avail) {
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
                            if (out_alleles[bi] >= 0) out_avail[bi] = 1;
                        }
                        return;
                    }
                    // Resolve a held-out haplotype's truth by running the SAME allele enumeration
                    // that built the panel, over a graph whose paths are the held-out ones. The two
                    // used to be separate code paths -- the panel enumerated alleles while the truth
                    // re-spelled the walk by hand -- and they drifted, so leave-one-out scores were
                    // measured against a truth the panel itself disagreed with. One routine used for
                    // both is the only way to keep that from happening again.
                    for (std::size_t bi = 0; bi < chain.size(); ++bi) {
                        const auto it = held_blocks[bi].allele_of.find(name);
                        if (it == held_blocks[bi].allele_of.end()) continue;
                        out_avail[bi] = 1;   // it traverses the block; the sequence may still be empty
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
                            if (blocks[bi].allele_seq[ai] == seq) { out_alleles[bi] = static_cast<int>(ai); break; }
                        }
                        // Kept even when no allele matches: an unrepresentable block still has a true
                        // sequence, and that is exactly the case the graded score exists to measure.
                        out_seq[bi] = seq;
                    }
                };
                ts1.assign(chain.size(), std::string());
                ts2.assign(chain.size(), std::string());
                tav1.assign(chain.size(), 0);
                tav2.assign(chain.size(), 0);
                resolve(a, ta1, ts1, tav1);
                resolve(b, ta2, ts2, tav2);

                // Depth calibration: with the truth known, every marker's expected count is
                // lambda*(m_truth1 + m_truth2), so the reads measure lambda directly. Regressing the
                // observed counts on the true multiplicity gives the depth the DATA implies, against
                // the depth the model is using. The two must agree, and a copy-number call is exactly
                // as wrong as they disagree: at an array of N copies a 1/N relative error in lambda
                // moves the call by one whole repeat unit.
                //
                // Stratified by multiplicity because a constant offset and a wrong slope look the same
                // in the aggregate but mean different things -- the first is a background term, the
                // second is the depth itself.
                if (depth_calibration) {
                    const std::string cpath = out_prefix + ".depth.calibration.tsv";
                    std::ofstream c(cpath);
                    if (!c) throw std::runtime_error("genotype: cannot write " + cpath);
                    // fano = observed variance of the residual divided by the predicted mean. A
                    // Poisson count gives 1. The emission uses a negative binomial whose dispersion is
                    // fitted on the depth ANCHORS, which are single-copy by construction, so this is
                    // the check that the fit still holds where the copy-number signal actually lives:
                    // at markers carried twenty or thirty times over.
                    c << "block_index\tblock_kind\tmult_class\tn_markers\tsum_mult\tsum_obs"
                         "\tlambda_implied\tlambda_used\tratio\tmean_pred\tfano\n";
                    const char* names[] = {"1", "2", "3-5", "6-10", ">10", "all"};
                    for (std::size_t bi = 0; bi < chain.size() && bi < read_panel.by_block.size(); ++bi) {
                        if (bi >= ta1.size() || ta1[bi] < 0 || ta2[bi] < 0) continue;
                        const auto& ba = read_panel.by_block[bi];
                        const std::size_t i1 = static_cast<std::size_t>(ta1[bi]);
                        const std::size_t i2 = static_cast<std::size_t>(ta2[bi]);
                        if (i1 >= ba.size() || i2 >= ba.size()) continue;
                        std::unordered_map<std::uint32_t, std::uint64_t> m;
                        for (const auto& [slot, k] : ba[i1].nodes) m[slot] += k;
                        for (const auto& [slot, k] : ba[i2].nodes) m[slot] += k;
                        double sm[6] = {0}, so[6] = {0}, sres[6] = {0}, spred[6] = {0};
                        std::size_t nm[6] = {0};
                        const double lam_b = bi < depth.size() ? depth[bi].lambda_hap : 0.0;
                        for (const auto& [slot, mult] : m) {
                            const std::size_t cls = mult == 1 ? 0 : mult == 2 ? 1 : mult <= 5 ? 2
                                                  : mult <= 10 ? 3 : 4;
                            const double obs = static_cast<double>(rc.node[slot]);
                            const double pred = lam_b * static_cast<double>(mult);
                            const double res = (obs - pred) * (obs - pred);
                            sm[cls] += static_cast<double>(mult); so[cls] += obs; ++nm[cls];
                            sm[5] += static_cast<double>(mult); so[5] += obs;  ++nm[5];
                            sres[cls] += res;   spred[cls] += pred;
                            sres[5] += res;     spred[5] += pred;
                        }
                        for (std::size_t cl = 0; cl < 6; ++cl) {
                            if (nm[cl] == 0) continue;
                            // Poisson ML for lambda in obs ~ Poisson(lambda*m): sum(obs)/sum(m).
                            const double implied = so[cl] / std::max(1.0, sm[cl]);
                            const double used = bi < depth.size() ? depth[bi].lambda_hap : 0.0;
                            c << bi << '\t'
                              << (chain[bi].kind == BlockKind::Bubble ? "bubble"
                                  : chain[bi].kind == BlockKind::Flank ? "flank" : "backbone")
                              << '\t' << names[cl] << '\t' << nm[cl] << '\t' << sm[cl] << '\t' << so[cl]
                              << '\t' << implied << '\t' << used << '\t'
                              << (used > 0.0 ? implied / used : 0.0) << '\t'
                              << spred[cl] / static_cast<double>(nm[cl]) << '\t'
                              << (spred[cl] > 0.0 ? sres[cl] / spred[cl] : 0.0) << '\n';
                        }
                    }
                    log.wrote({cpath});
                }
            }
            // Diagnostic: under leave-one-out the truth allele is often absent, so the emission-rank
            // diagnostic (which asks where the target pair sits by emission alone, ignoring linkage) has
            // nothing to aim at. Substitute the NEAREST AVAILABLE allele by length. If the nearest ranks
            // first by emission but is not called, something after the emission overrode it; if it ranks
            // poorly, the likelihood itself prefers a worse allele. Scoring is untouched -- the truth
            // check reads ta1/ta2 directly, not what is handed to genotype_sample.
            std::vector<int> pa1 = ta1;
            std::vector<int> pa2 = ta2;
            if (nearest_rank && !ta1.empty()) {
                auto nearest = [&](std::size_t bi, const std::string& truth) {
                    long best = -1;
                    long bestd = -1;
                    for (std::size_t ai = 0; ai < blocks[bi].allele_seq.size(); ++ai) {
                        const long d = std::labs(static_cast<long>(blocks[bi].allele_seq[ai].size()) -
                                                 static_cast<long>(truth.size()));
                        if (best < 0 || d < bestd) { best = static_cast<long>(ai); bestd = d; }
                    }
                    return static_cast<int>(best);
                };
                for (std::size_t bi = 0; bi < chain.size(); ++bi) {
                    if (bi < ts1.size() && !ts1[bi].empty() && pa1[bi] < 0) pa1[bi] = nearest(bi, ts1[bi]);
                    if (bi < ts2.size() && !ts2[bi].empty() && pa2[bi] < 0) pa2[bi] = nearest(bi, ts2[bi]);
                }
            }
            // Target the most IDENTICAL available allele per haplotype, and ask where the EMISSION
            // ranks that pair. This is the gate on every "we picked the wrong allele" claim: the
            // identity oracle chooses each haplotype's best independently, so the pair it names has
            // never been scored as a pair, and a target that no likelihood could prefer is not a defect.
            //
            // Length was the previous target and it was the wrong one -- it ignores sequence entirely,
            // and forcing the call to match it made identity worse. Identity is the meaningful oracle,
            // but it needs the same check before it is trusted.
            if (oracle_rank && !ta1.empty()) {
                const std::size_t kk = read_panel.kmer_size;
                const std::size_t ss = read_panel.syncmer_s != 0 ? read_panel.syncmer_s
                                                                 : default_syncmer_s(kk);
                auto sset = [&](const std::string& s) {
                    std::unordered_set<std::uint64_t> out;
                    for (const KmerOccurrence& o : collect_syncmers(s, kk, ss)) out.insert(o.code);
                    return out;
                };
                // Jaccard over every allele to shortlist, exact alignment on the shortlist. Jaccard
                // alone will not do: it is set-valued, so it ignores multiplicity, and inside a tandem
                // array two alleles with the same unit repertoire and different copy numbers score
                // alike -- blind to the one quantity the call has to get right.
                auto most_identical = [&](std::size_t bi, const std::string& truth) {
                    if (blocks[bi].allele_seq.empty() || blocks[bi].allele_seq.size() > 1024) return -1;
                    const auto tset = sset(truth);
                    std::vector<std::pair<double, std::size_t>> jac;
                    for (std::size_t ai = 0; ai < blocks[bi].allele_seq.size(); ++ai) {
                        const auto A = sset(blocks[bi].allele_seq[ai]);
                        std::size_t inter = 0;
                        for (const std::uint64_t c : (A.size() <= tset.size() ? A : tset)) {
                            if ((A.size() <= tset.size() ? tset : A).count(c) != 0) ++inter;
                        }
                        const std::size_t uni = A.size() + tset.size() - inter;
                        jac.emplace_back(uni == 0 ? 1.0 : static_cast<double>(inter) / static_cast<double>(uni), ai);
                    }
                    std::sort(jac.begin(), jac.end(), [](const auto& x, const auto& y) {
                        return x.first != y.first ? x.first > y.first : x.second < y.second;
                    });
                    long best = -1;
                    double bid = -1.0;
                    for (std::size_t r = 0; r < std::min<std::size_t>(16, jac.size()); ++r) {
                        const NwAlign n2 = nw_edit_distance(blocks[bi].allele_seq[jac[r].second], truth);
                        const double id = 1.0 - static_cast<double>(n2.edits) /
                                                static_cast<double>(std::max<std::size_t>(1, n2.aln_len));
                        if (id > bid) { bid = id; best = static_cast<long>(jac[r].second); }
                    }
                    return static_cast<int>(best);
                };
                for (std::size_t bi = 0; bi < chain.size(); ++bi) {
                    if (bi < ts1.size() && !ts1[bi].empty()) pa1[bi] = most_identical(bi, ts1[bi]);
                    if (bi < ts2.size() && !ts2[bi].empty()) pa2[bi] = most_identical(bi, ts2[bi]);
                }
            }
            // Coverage evidence, built once and handed to the same chain the marker path uses.
            CoverageEvidence cev;
            if (evidence != "syncmer") {
                const NodeIndex nidx = build_node_index(panel_graph);
                const PanelCoverage pcov = build_panel_coverage(panel_graph, nidx);
                CoverageOptions copt;
                copt.threads = options.threads;
                const SampleCoverage scov = inject_reads(panel_graph, nidx, pcov, read_paths, copt);
                // Depth from nodes every panel path traverses exactly once: the sample carries two
                // copies of those whatever its genotype, so their coverage measures depth alone.
                std::vector<double> anchor_cov;
                for (std::size_t n = 0; n < nidx.size(); ++n) {
                    bool inv = !pcov.by_path.empty();
                    for (const auto& v : pcov.by_path) if (v[n] != 1) { inv = false; break; }
                    if (inv) anchor_cov.push_back(scov.node[n]);
                }
                std::sort(anchor_cov.begin(), anchor_cov.end());
                cev.lambda = anchor_cov.size() >= 20 ? anchor_cov[anchor_cov.size() / 2] / 2.0 : 0.0;
                cev.node = scov.node;
                cev.block_allele_nodes.resize(chain.size());
                cev.use_block.assign(chain.size(), 0);
                cev.target_bp.assign(chain.size(), 0.0);
                cev.target_sd.assign(chain.size(), 0.0);
                std::size_t n_used = 0;
                for (std::size_t bi = 0; bi < chain.size(); ++bi) {
                    auto av = block_allele_node_vectors(panel_graph, path_indexes, bubbles,
                                                        chain[bi], blocks[bi], nidx);
                    std::uint32_t mx = 0;
                    for (const auto& v : av) for (const std::uint32_t m : v) mx = std::max(mx, m);
                    // `auto` routes only the tandem arrays to coverage, on the argument that the marker
                    // path is validated everywhere else and a marker rule cannot express copy number.
                    // Whether routing beats using coverage everywhere is a measurement, not a belief.
                    const bool is_array = mx >= 3;
                    const bool want = evidence == "coverage" ||
                                      ((evidence == "auto" || evidence == "combined") && is_array);
                    if (want) { cev.use_block[bi] = 1; ++n_used; }
                    // How much sequence the coverage says is here, across both haplotypes: the observed
                    // coverage over the block's nodes, divided by the depth, times each node's length.
                    // Alignment-derived, so independent of the markers it will be used to constrain.
                    if (want && cev.lambda > 0.0) {
                        std::vector<char> seen(nidx.size(), 0);
                        double bp = 0.0;
                        for (const auto& v : av) {
                            for (std::size_t n = 0; n < v.size(); ++n) {
                                if (v[n] > 0 && !seen[n]) {
                                    seen[n] = 1;
                                    bp += (scov.node[n] / cev.lambda) * static_cast<double>(nidx.length[n]);
                                }
                            }
                        }
                        cev.target_bp[bi] = bp;
                        // The estimator is good to about 1%; 2% is deliberately loose,
                        // so the constraint rules out whole repeat units without arbitrating between
                        // pairs that are both plausible -- that arbitration is the markers' job.
                        cev.target_sd[bi] = 0.02 * bp;
                    }
                    cev.block_allele_nodes[bi] = std::move(av);
                }
                // In combined mode the markers keep the emission and coverage contributes only the
                // total-length constraint. The emission fires on the per-block allele vectors, so
                // dropping those switches it off while use_block stays set and the target still
                // applies.
                if (evidence == "combined") {
                    for (auto& v : cev.block_allele_nodes) v.clear();
                    for (std::size_t bi = 0; bi < cev.use_block.size(); ++bi) {
                        if (cev.target_bp[bi] <= 0.0) cev.use_block[bi] = 0;
                    }
                }
                log.info("evidence " + evidence + ": coverage at " + std::to_string(n_used) + " of " +
                         std::to_string(chain.size()) + " blocks, lambda " + std::to_string(cev.lambda) +
                         " from " + std::to_string(anchor_cov.size()) + " invariant nodes");
            }
            // ---- --noiseless-counts -----------------------------------------------------------
            // Overwrite one block's observed marker counts with what a stated pair would produce at
            // this sample's depth with no read noise. Everything else -- panel, depth model, pruning,
            // emission, chain -- is the production path untouched, so a wrong call afterwards cannot
            // be charged to read sampling, sequencing error or marker acquisition.
            //
            // The counts are written as lambda*(m1+m2)+mu, which is exactly the emission's own mean
            // for that pair. A per-marker negative binomial is maximised at mean == observation, so a
            // proper likelihood MUST rank the source pair first when it is in the panel. That is why
            // the BLK:A,B form exists: it is a control the model cannot fail for an honest reason.
            if (noiseless_block >= 0) {
                const std::size_t nbi = static_cast<std::size_t>(noiseless_block);
                if (nbi >= chain.size())
                    throw std::runtime_error("--noiseless-counts block " +
                                             std::to_string(noiseless_block) +
                                             " is out of range (chain has " +
                                             std::to_string(chain.size()) + " blocks)");
                // lambda exactly as the emission reads it (genotype.cpp), so the injected counts and
                // the means they are compared against come from one number, not two.
                const double lambda = depth[nbi].usable ? depth[nbi].lambda_hap : 1.0;
                double lam = 0.0;
                std::size_t nlam = 0;
                for (const BlockDepth& d : depth) { if (d.usable) { lam += d.lambda_hap; ++nlam; } }
                lam = nlam ? lam / static_cast<double>(nlam) : 1.0;
                const double mu = gopt.error_background > 0.0 ? gopt.error_background
                                                              : std::max(0.01, 0.02 * lam);

                // The two source sequences.
                std::string src1, src2;
                std::string src_desc;
                if (noiseless_a >= 0) {
                    const std::size_t na = blocks[nbi].allele_seq.size();
                    if (static_cast<std::size_t>(noiseless_a) >= na ||
                        static_cast<std::size_t>(noiseless_b) >= na)
                        throw std::runtime_error("--noiseless-counts allele out of range: block " +
                                                 std::to_string(noiseless_block) + " has " +
                                                 std::to_string(na) + " alleles, asked for (" +
                                                 std::to_string(noiseless_a) + "," +
                                                 std::to_string(noiseless_b) + ")");
                    src1 = blocks[nbi].allele_seq[static_cast<std::size_t>(noiseless_a)];
                    src2 = blocks[nbi].allele_seq[static_cast<std::size_t>(noiseless_b)];
                    src_desc = "panel pair (" + std::to_string(noiseless_a) + "," +
                               std::to_string(noiseless_b) + ")";
                } else {
                    // The truth's own spelled sequence, which under leave-one-out is off-panel.
                    // Traversal is what decides availability: an empty string is a legitimate
                    // deletion allele and must not read as missing data.
                    if (nbi >= tav1.size() || (!tav1[nbi] && !tav2[nbi]))
                        throw std::runtime_error("--noiseless-counts <BLK>: neither truth haplotype "
                                                 "traverses block " + std::to_string(noiseless_block));
                    src1 = ts1[nbi];
                    src2 = ts2[nbi];
                    src_desc = "the truth haplotypes";
                }

                // Marker multiplicity of a sequence, extracted exactly as the panel extracted its
                // alleles'. A code the panel never kept as a marker has no slot and is dropped, which
                // is the same thing that happens to it in a real read.
                std::unordered_map<std::uint64_t, std::vector<std::uint32_t>> slot_of;
                slot_of.reserve(read_panel.node_codes.size() * 2);
                for (std::size_t s = 0; s < read_panel.node_codes.size(); ++s) {
                    slot_of[read_panel.node_codes[s]].push_back(static_cast<std::uint32_t>(s));
                }
                std::unordered_map<std::uint32_t, std::uint32_t> mult;
                const auto add = [&](const std::string& seq) {
                    if (seq.empty()) return;
                    const std::vector<KmerOccurrence> sy =
                        read_panel.all_kmers
                            ? collect_canonical_kmer_occurrences(seq, read_panel.kmer_size)
                            : collect_syncmers(seq, read_panel.kmer_size, read_panel.syncmer_s);
                    for (const KmerOccurrence& o : sy) {
                        const auto it = slot_of.find(o.code);
                        if (it == slot_of.end()) continue;
                        for (const std::uint32_t s : it->second) ++mult[s];
                    }
                };
                add(src1);
                add(src2);

                // Every marker of the block, not only the ones the source carries: a marker the pair
                // lacks is predicted at mu, and writing mu there is what makes the source the
                // per-marker maximum. Leaving the observed count would keep the veto the probe exists
                // to remove.
                std::vector<std::uint32_t> universe;
                for (const auto& mset : read_panel.by_block[nbi]) {
                    for (const auto& [slot, m] : mset.nodes) { (void)m; universe.push_back(slot); }
                }
                std::sort(universe.begin(), universe.end());
                universe.erase(std::unique(universe.begin(), universe.end()), universe.end());
                // Split by whether the SOURCE carries the marker, so each half can be rewritten on
                // its own. A marker the source lacks is predicted at mu, and mu is well under 1, so
                // rounding sends it to exactly 0 -- the sharpest possible "absent" reading, and
                // sharper than any real read pile-up, which scatters stray counts there. That is the
                // asymmetry the `absent` scope exists to test, and it is worth being precise about:
                // the treatment is "set to zero", not "set to mu", because counts are integers.
                if (!injection_out.empty()) {
                    std::ofstream jf(injection_out);
                    if (!jf) throw std::runtime_error("genotype: cannot write " + injection_out);
                    // first_pos and clump make the per-class decomposition answerable at the
                    // right unit. Adjacent syncmers a few bp apart sit inside one fragment and rise
                    // and fall together, so N markers disagreeing is not N independent pieces of
                    // evidence -- it may be one sequence difference counted many times. The block's
                    // rho discount already scales the total by clumps/markers, but it applies AFTER
                    // correlated markers have each voted, so a class total in marker units cannot say
                    // how many real differences are behind it.
                    const double frag = read_panel.fragment_len > 0.0 ? read_panel.fragment_len : 350.0;
                    jf << "slot\tsource_mult\tobserved_count\tfirst_pos\tclump\n";
                    std::vector<std::uint32_t> uni;
                    for (const auto& mset : read_panel.by_block[nbi])
                        for (const auto& [slot, m] : mset.nodes) { (void)m; uni.push_back(slot); }
                    std::sort(uni.begin(), uni.end());
                    uni.erase(std::unique(uni.begin(), uni.end()), uni.end());
                    for (const std::uint32_t sl : uni) {
                        const auto it = mult.find(sl);
                        const std::uint32_t fp = sl < read_panel.node_first_pos.size()
                                                     ? read_panel.node_first_pos[sl] : UINT32_MAX;
                        jf << sl << '\t' << (it == mult.end() ? 0u : it->second) << '\t'
                           << rc.node[sl] << '\t';
                        // No recorded position means no clump can be assigned; NA, never 0, which
                        // would silently pool every such marker into one enormous clump.
                        if (fp == UINT32_MAX) jf << "NA\tNA\n";
                        else jf << fp << '\t' << static_cast<long>(fp / frag) << '\n';
                    }
                    log.wrote({injection_out});
                }
                std::size_t n_pres = 0, n_abs = 0;
                double pres_before = 0.0, pres_after = 0.0, abs_before = 0.0, abs_after = 0.0;
                for (const std::uint32_t s : universe) {
                    const auto it = mult.find(s);
                    const double m = it == mult.end() ? 0.0 : static_cast<double>(it->second);
                    const bool present = m > 0.0;
                    if (present) { ++n_pres; pres_before += rc.node[s]; }
                    else         { ++n_abs;  abs_before  += rc.node[s]; }
                    const bool rewrite = noiseless_scope == "all" ||
                                         (noiseless_scope == "present" && present) ||
                                         (noiseless_scope == "absent" && !present);
                    if (rewrite) {
                        const double want = lambda * m + mu;
                        // Counts are integers, so the injection is noiseless only up to rounding --
                        // at most half a count per marker, against a lambda of tens.
                        rc.node[s] =
                            static_cast<std::uint32_t>(std::max<long long>(0, std::llround(want)));
                    }
                    if (present) pres_after += rc.node[s];
                    else         abs_after  += rc.node[s];
                }
                const auto mass = [](double v) {
                    return std::to_string(static_cast<long long>(v));
                };
                log.info("NOISELESS: block " + std::to_string(noiseless_block) + ", scope " +
                         noiseless_scope + ", source " + src_desc + ", lambda " +
                         std::to_string(lambda) + ", mu " + std::to_string(mu));
                log.info("NOISELESS:   markers the source carries: " + std::to_string(n_pres) +
                         ", mass " + mass(pres_before) + " -> " + mass(pres_after));
                log.info("NOISELESS:   markers it does not:        " + std::to_string(n_abs) +
                         ", mass " + mass(abs_before) + " -> " + mass(abs_after) +
                         (noiseless_scope == "all" || noiseless_scope == "absent"
                              ? "  (rewritten to 0, not to mu: counts are integers and mu < 1)" : ""));
                log.info("NOISELESS: this run no longer describes the supplied reads at block " +
                         std::to_string(noiseless_block));
            }
            // Applied before the FIRST call, which always runs: the joint-depth refinement below
            // fires only under some depth conditions, so patching it there made the flag inert
            // exactly where the diagnosis was needed.
            // EVERY probe is range-checked, not only the first. rank -2 means "pruned before
            // scoring", so an out-of-range index silently becomes a pruning finding -- which is the
            // conclusion the flag exists to establish. Validating only the first left the later ones
            // of a repeated probe able to forge exactly that.
            for (const auto& pp : probe_pairs) {
                if (pp[0] >= chain.size())
                    throw std::runtime_error("--probe-pair block " + std::to_string(pp[0]) +
                                             " is out of range (chain has " +
                                             std::to_string(chain.size()) + " blocks)");
                const std::size_t na = blocks[pp[0]].allele_seq.size();
                if (pp[1] >= na || pp[2] >= na)
                    throw std::runtime_error("--probe-pair allele out of range: block " +
                                             std::to_string(pp[0]) + " has " + std::to_string(na) +
                                             " alleles, asked for (" + std::to_string(pp[1]) + "," +
                                             std::to_string(pp[2]) + ")");
            }
            if (probe_block >= 0) {
                if (pa1.size() < chain.size()) pa1.resize(chain.size(), -1);
                if (pa2.size() < chain.size()) pa2.resize(chain.size(), -1);
                pa1[static_cast<std::size_t>(probe_block)] = static_cast<int>(probe_a);
                pa2[static_cast<std::size_t>(probe_block)] = static_cast<int>(probe_b);
                log.info("PROBE: block " + std::to_string(probe_block) + " truth columns now describe "
                         "the pair (" + std::to_string(probe_a) + "," + std::to_string(probe_b) +
                         "), NOT this sample's truth");
            }
            gopt.probe_pairs = probe_pairs;
            // ACTIVE EDGES ONLY, and only from a committed transaction. Null otherwise, which takes
            // the factorised path at every edge and reproduces the legacy chain exactly.
            if (hyb_act.hybrid_activated && !higher_active)
                gopt.sparse_linkage_edges = &hyb_act.sparse_kernel_edges;
            if (higher_active) {
                // ONE LIST OWNS EVERY FACTOR. The sparse edges stay null precisely so nothing can
                // be applied twice.
                gopt.higher_factors = &higher_factors;
                gopt.hap_allele = &higher_hap_allele;
                gopt.higher_max_message_entries = hybrid_max_message_entries;
                gopt.higher_stats = &higher_stats;
            }
            std::vector<ProbePairResult> probe_rows;
            std::vector<BlockCall> calls =
                genotype_sample(chain, blocks, read_panel, rc, depth, hap_names, gopt, &gsum,
                                pa1.empty() ? nullptr : &pa1, pa2.empty() ? nullptr : &pa2,
                                evidence == "syncmer" ? nullptr : &cev, &probe_rows);

            // ---- THE BLOCK CATALOGUE ---------------------------------------------------------
            // Exactly the projection the call was made in: allele_seq is what the reported allele
            // indices index into, and allele_of is how a panel path maps onto them. An evaluator
            // reading this cannot drift from the caller the way one rebuilding blocks from graph
            // coordinates can. The catalogue fingerprint is emitted so a harness can assert that
            // two runs indexed the same catalogue before comparing any allele number between them.
            if (!dump_block_catalogue.empty()) {
                std::ofstream bc(dump_block_catalogue);
                if (!bc) throw std::runtime_error("genotype: cannot write " +
                                                  dump_block_catalogue);
                bc << "# catalogue_fingerprint\t" << allele_catalogue_fingerprint(blocks) << '\n';
                bc << "# blocks\t" << blocks.size() << '\n';
                bc << "# paths\t" << hap_names.size() << '\n';
                bc << "kind\tblock\tkey\tvalue\n";
                for (std::size_t b = 0; b < blocks.size(); ++b) {
                    bc << "block_kind\t" << b << "\t-\t"
                       << (chain[b].kind == BlockKind::Bubble ? "bubble"
                           : chain[b].kind == BlockKind::Backbone ? "backbone" : "flank")
                       << '\n';
                    for (std::size_t ai = 0; ai < blocks[b].allele_seq.size(); ++ai)
                        bc << "allele\t" << b << '\t' << ai << '\t'
                           << blocks[b].allele_seq[ai] << '\n';
                    for (const std::string& nm : hap_names) {
                        const auto it = blocks[b].allele_of.find(nm);
                        bc << "path\t" << b << '\t' << nm << '\t'
                           << (it == blocks[b].allele_of.end()
                                   ? std::string("-") : std::to_string(it->second)) << '\n';
                    }
                }
                log.info("block catalogue: " + std::to_string(blocks.size()) + " blocks, " +
                         std::to_string(hap_names.size()) + " paths -> " + dump_block_catalogue);
            }

            // ---- WHAT THE RECURRENCE ACTUALLY DID, after it has done it -----------------------
            if (!hybrid_higher_report.empty()) {
                std::ofstream hr(hybrid_higher_report);
                if (!hr) throw std::runtime_error("genotype: cannot write " +
                                                  hybrid_higher_report);
                hr << "field\tvalue\n";
                hr << "higher_order_active\t" << (higher_active ? 1 : 0) << '\n';
                hr << "refusal\t" << (higher_refusal.empty() ? "-" : higher_refusal) << '\n';
                if (higher_active) {
                    const bool exact = higher_stats.all_initial_emissions_finite;
                    hr << "plan_contract\t" << (exact ? "EQUAL" : "BOUND") << '\n';
                    hr << "all_initial_emissions_finite\t" << (exact ? 1 : 0) << '\n';
                    hr << "initial_states_dropped\t"
                       << higher_stats.initial_states_dropped << '\n';
                    hr << "planned_forward_updates\t" << higher_plan.forward_updates << '\n';
                    hr << "actual_forward_updates\t" << higher_stats.forward_updates << '\n';
                    hr << "actual_adjoint_updates\t" << higher_stats.adjoint_updates << '\n';
                    hr << "planned_peak_message_entries\t"
                       << higher_plan.peak_message_entries << '\n';
                    hr << "actual_peak_message_entries\t"
                       << higher_stats.peak_message_entries << '\n';
                    const bool holds = exact
                        ? (higher_stats.forward_updates == higher_plan.forward_updates &&
                           higher_stats.peak_message_entries == higher_plan.peak_message_entries)
                        : (higher_stats.forward_updates <= higher_plan.forward_updates &&
                           higher_stats.peak_message_entries <=
                               higher_plan.peak_message_entries);
                    hr << "plan_contract_holds\t" << (holds ? 1 : 0) << '\n';
                    hr << "factor_lookups\t" << higher_stats.factor_lookups << '\n';
                    hr << "history_ops\t" << higher_stats.history_ops << '\n';
                    hr << "grouping_ops\t" << higher_stats.grouping_ops << '\n';
                    hr << "forward_seconds\t" << higher_stats.forward_seconds << '\n';
                    hr << "adjoint_seconds\t" << higher_stats.adjoint_seconds << '\n';
                    hr << "blocks_called\t" << calls.size() << '\n';
                }
            }

            if (model_pangenie) {
                // Same panel, same counts, same depth -- only the model differs, which is the whole
                // point: any difference is attributable to the emission, transition and marker rule
                // rather than to inputs.
                PanGenieOptions po;
                po.threads = options.threads;
                const std::vector<PanGenieCall> pg =
                    genotype_pangenie(chain, blocks, read_panel, rc, depth, hap_names, po);
                std::size_t undef = 0;
                for (std::size_t bi = 0; bi < calls.size() && bi < pg.size(); ++bi) {
                    calls[bi].allele1 = pg[bi].allele1;
                    calls[bi].allele2 = pg[bi].allele2;
                    calls[bi].gq = pg[bi].gq;
                    calls[bi].n_markers = pg[bi].n_markers;
                    calls[bi].evidence = pg[bi].undefined ? "none" : "local";
                    calls[bi].filter = pg[bi].undefined ? "NOMARKERS"
                                                        : (pg[bi].gq < gopt.min_gq ? "LOWGQ" : "PASS");
                    if (pg[bi].undefined) ++undef;
                }
                log.info("pangenie model: " + std::to_string(calls.size()) + " blocks, " +
                         std::to_string(undef) + " with no usable marker under their unique-once rule");
            }

            if (run_joint) {
                // The anchors of a block record lambda times the number of the sample's haplotypes
                // that traverse it, and that count is part of the genotype. Take it from the first
                // pass: divide each block's anchor median by how many of its called alleles are not
                // the bypass allele, then use the median of those as one region-wide lambda and call
                // again. One refinement is enough -- the count only ever takes the values 0, 1 or 2.
                std::vector<double> per_hap;
                for (std::size_t bi = 0; bi < depth.size() && bi < calls.size(); ++bi) {
                    // `local_center`, not `anchor_median`: the latter is deliberately always the
                    // median so the audit column is comparable across runs, so reading it here would
                    // silently ignore --depth-estimator in the DEFAULT model. `local_available`
                    // rather than a >0 test, because 0 is a legitimate anchor count.
                    if (!depth[bi].local_available || depth[bi].local_center <= 0.0) continue;
                    const int bp = bi < blocks.size() ? blocks[bi].bypass_allele : -1;
                    int traversing = 0;
                    if (static_cast<int>(calls[bi].allele1) != bp) ++traversing;
                    if (static_cast<int>(calls[bi].allele2) != bp) ++traversing;
                    if (traversing > 0) per_hap.push_back(depth[bi].local_center / traversing);
                }
                if (!per_hap.empty()) {
                    std::sort(per_hap.begin(), per_hap.end());
                    const double lambda = per_hap[per_hap.size() / 2];
                    log.info("joint depth: lambda " + std::to_string(lambda) +
                             " from " + std::to_string(per_hap.size()) +
                             " blocks, using the first pass's traversal counts");
                    for (BlockDepth& d : depth) {
                        d.median = lambda * 2.0;
                        d.lambda_hap = lambda;
                        d.usable = true;
                        d.source = DepthSource::Joint;
                        d.region_shrink_weight = 1.0;   // region-wide; raw and local_* stay untouched
                    }
                    // Same arguments as the first pass. Dropping the coverage evidence here silently
                    // reverted --evidence to the marker model whenever a bypass allele existed, which
                    // is exactly the case the joint pass exists for -- so the flag appeared to have no
                    // effect at precisely those loci. Passing ta instead of pa did the same to the
                    // emission-rank diagnostics.
                    // The second pass supersedes the first, so the probe table must describe THIS
                    // run and not hold both. Appending would silently double every row.
                    probe_rows.clear();
                    calls = genotype_sample(chain, blocks, read_panel, rc, depth, hap_names, gopt,
                                            &gsum, pa1.empty() ? nullptr : &pa1,
                                            pa2.empty() ? nullptr : &pa2,
                                            evidence == "syncmer" ? nullptr : &cev, &probe_rows);
                    log_depth_provenance(log, chain, depth);
                }
            }
            // Deferred from the first-pass logging so it describes the depth the emission actually
            // used, joint refinement included.
            write_read_audit(out_prefix, chain, read_panel, rc, depth);
            if (!dump_markers.empty()) {
                write_marker_dump(dump_markers, chain, read_panel, rc, depth,
                                  ts1.empty() ? nullptr : &ts1, ts2.empty() ? nullptr : &ts2,
                                  ta1.empty() ? nullptr : &ta1, ta2.empty() ? nullptr : &ta2);
                log.wrote({dump_markers});
            }
            log.info("model: lambda " + std::to_string(gsum.lambda_hap) + ", overdispersion phi " +
                     std::to_string(gsum.overdispersion) + ", error background " +
                     std::to_string(gsum.error_background));
            log.info("calls: " + std::to_string(gsum.called) + " PASS, " + std::to_string(gsum.no_calls) +
                     " no-call, " + std::to_string(gsum.off_panel) + " off-panel; mean GQ " +
                     std::to_string(gsum.mean_gq));

            // Say it in the log, not only in a column. At a tandem array the allele pair is the closest
            // panel allele BY CONTENT and its length is not the copy-number answer -- mass_bp is -- and
            // a reader who takes called_bp for the copy number gets it wrong by up to a whole repeat
            // unit. That is easy to do silently, so it is stated wherever such a block is called.
            {
                std::size_t n_arr = 0;
                for (const BlockCall& c : calls) if (c.is_array) ++n_arr;
                if (n_arr > 0) {
                    log.info(std::to_string(n_arr) + " block(s) are tandem arrays (block_class=array): "
                             "there the called allele pair is the closest panel allele BY CONTENT, and "
                             "copy number is mass_bp +- mass_bp_sd, not called_bp");
                }
            }
            write_genotypes(out_prefix, chain, blocks, calls, hap_names, probe_block >= 0);

            // The probe table. One row per --probe-pair, in the order given. `score` is raw and only
            // comparable within this run; `delta` is against this run's block optimum, which moves
            // when the counts do -- so to compare arms of an experiment, probe a FIXED reference pair
            // in every arm and take the difference of the two deltas, where the baseline cancels.
            if (!probe_pairs.empty()) {
                const std::string pp_path = out_prefix + ".probe_pairs.tsv";
                std::ofstream pf(pp_path);
                if (!pf) throw std::runtime_error("genotype: cannot write " + pp_path);
                pf << "block_index\tallele1\tallele2\trank\tties\tn_scored\tdelta\tscore\n";
                pf.setf(std::ios::fixed);
                for (const ProbePairResult& r : probe_rows) {
                    pf << r.block_index << '\t' << r.allele1 << '\t' << r.allele2 << '\t'
                       << r.rank << '\t' << r.ties << '\t' << r.n_scored << '\t';
                    // rank -2 means pruned before scoring: there is no emission, so a number in the
                    // score columns would be read as one.
                    if (r.rank == -2) pf << "NA\tNA\n";
                    // Ten decimals, not six. These numbers exist to be SUBTRACTED across runs -- a
                    // factorial design takes four of them -- and at scores in the thousands, six
                    // decimals leaves rounding noise the size of the effect being measured.
                    else pf << std::setprecision(10) << r.delta << '\t'
                            << std::setprecision(10) << r.score << '\n';
                }
                log.wrote({pp_path});
            }

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

                // ---- certified nearest-pair oracle -------------------------------------------
                // The existing oracle shortlists by syncmer-set Jaccard and aligns only the top 16, so
                // it certifies nothing: Jaccard is order- and multiplicity-blind, which is exactly
                // wrong at a tandem array. This one is exact and costs 2A alignments rather than
                // A(A+1)/2 pairs -- cache each allele against each truth haplotype, then every pair is
                // read from the cache. At 457 alleles that is 914 alignments instead of 104,653.
                //
                // Three criteria, because they do not agree. Minimum summed edit distance is not the
                // same ordering as maximum mean identity when the two alignments differ in length, and
                // neither is the same as minimum length error. Whichever is reported must be the one
                // optimised, so all three are.
                //
                // Per-haplotype residuals, not their sum: a summed length error of ~0 is equally
                // consistent with both haplotypes right and with one too long by exactly what the
                // other is too short. That distinction is the whole question of representability.
                if (certified_oracle && !blocks.empty() && !blocks[0].allele_seq.empty()) {
                    const std::string op = out_prefix + ".oracle.tsv";
                    std::ofstream orc(op);
                    if (!orc) throw std::runtime_error("genotype: cannot write " + op);
                    // Opened once, like orc. Opening it inside the per-block loop meant every block
                    // truncated the previous one's rows and the file described only the last block.
                    const std::string pp = out_prefix + ".oracle_pairs.tsv";
                    std::ofstream opf;
                    if (!oracle_pairs.empty()) {
                        opf.open(pp);
                        if (!opf) throw std::runtime_error("genotype: cannot write " + pp);
                        opf << "block_index\tallele1\tallele2\th1_id\th2_id\th1_lenerr"
                               "\th2_lenerr\ttotal_edits\tmean_id\texcess_total_edits\n";
                    }
                    orc << "block_index\tblock_kind\tbubble_id\tn_alleles\tcriterion"
                           "\tbest_a\tbest_b\tbest_h1_id\tbest_h2_id\tbest_h1_lenerr\tbest_h2_lenerr"
                           "\tbest_total_edits\tbest_mean_id\tn_tied"
                           "\tcalled_a\tcalled_b\tcalled_rank\tcalled_h1_id\tcalled_h2_id"
                           "\tcalled_h1_lenerr\tcalled_h2_lenerr\tcalled_mean_id"
                           // The quantity to prioritise on: how many edits the call costs OVER the
                           // best pair the reduced panel could have offered. Identity is a mean of
                           // two per-haplotype RATES, each normalised by its own alignment length,
                           // so scaling it by a block's bp only approximates this -- and the earlier
                           // `best_identity` is a top-16 Jaccard shortlist, a lower bound. These are
                           // exact and certified over all 2A alignments.
                           "\tcalled_total_edits\texcess_total_edits\n";
                    for (std::size_t bi = 0; bi < chain.size() && bi < blocks.size(); ++bi) {
                        // Skip only when NEITHER haplotype traverses the block. Both carrying the
                        // deletion is a genotype, not a missing value.
                        if (bi >= ts1.size() || (!tav1[bi] && !tav2[bi])) continue;
                        if (certified_oracle_block >= 0 &&
                            static_cast<long>(chain[bi].index) != certified_oracle_block) continue;
                        const std::size_t A = blocks[bi].allele_seq.size();
                        if (A == 0 || A > certified_oracle_max) continue;
                        // Progress, because 2A exact alignments of whole block sequences is minutes to
                        // hours and the output stream buffers -- an empty file says nothing about how
                        // far it has got, which is how the first run became unobservable.
                        log.info("certified oracle: block " + std::to_string(chain[bi].index) + ", " +
                                 std::to_string(A) + " alleles, " +
                                 std::to_string(oracle_called_only ? 4 : 2 * A) + " alignments" +
                                 (oracle_called_only ? " (called pair only; best_* will be NA)" : ""));
                        const std::string* tv[2] = {&ts1[bi], &ts2[bi]};

                        // 2A alignments, the entire cost of the exercise -- unless only the called
                        // pair is wanted, in which case it is 2.
                        const std::size_t ca0 = calls[bi].allele1, cb0 = calls[bi].allele2;
                        std::vector<std::array<NwAlign, 2>> al(A);
                        std::vector<std::array<long, 2>> dl(A);
                        // Under --oracle-called-only the cache is trimmed to the called pair; any
                        // allele named by --oracle-pair must survive that trim or it would be priced
                        // from an alignment that was never computed.
                        std::vector<char> want_allele(A, 0);
                        for (const auto& pr : oracle_pairs) {
                            if (pr[0] < A) want_allele[pr[0]] = 1;
                            if (pr[1] < A) want_allele[pr[1]] = 1;
                        }
                        for (std::size_t a = 0; a < A; ++a) {
                            if (oracle_called_only && a != ca0 && a != cb0 && !want_allele[a]) continue;
                            for (int h = 0; h < 2; ++h) {
                                al[a][static_cast<std::size_t>(h)] =
                                    nw_edit_distance(blocks[bi].allele_seq[a], *tv[h]);
                                dl[a][static_cast<std::size_t>(h)] =
                                    static_cast<long>(blocks[bi].allele_seq[a].size()) -
                                    static_cast<long>(tv[h]->size());
                            }
                        }
                        auto ident = [&](std::size_t a, int h) {
                            const NwAlign& n = al[a][static_cast<std::size_t>(h)];
                            return 1.0 - static_cast<double>(n.edits) /
                                             static_cast<double>(std::max<std::size_t>(1, n.aln_len));
                        };
                        // For a pair, take the better of the two haplotype assignments UNDER THE
                        // CRITERION BEING SCORED, not under a single fixed one.
                        struct PairScore { double v; int swap; };
                        auto score = [&](std::size_t a, std::size_t b, int crit) -> PairScore {
                            const double d0 = static_cast<double>(al[a][0].edits + al[b][1].edits);
                            const double d1 = static_cast<double>(al[a][1].edits + al[b][0].edits);
                            const double i0 = (ident(a, 0) + ident(b, 1)) / 2.0;
                            const double i1 = (ident(a, 1) + ident(b, 0)) / 2.0;
                            const double l0 = static_cast<double>(std::labs(dl[a][0]) + std::labs(dl[b][1]));
                            const double l1 = static_cast<double>(std::labs(dl[a][1]) + std::labs(dl[b][0]));
                            if (crit == 0) return d0 <= d1 ? PairScore{d0, 0} : PairScore{d1, 1};
                            if (crit == 1) return i0 >= i1 ? PairScore{-i0, 0} : PairScore{-i1, 1};
                            return l0 <= l1 ? PairScore{l0, 0} : PairScore{l1, 1};
                        };
                        const char* names[3] = {"edit_distance", "mean_identity", "length_error"};
                        const std::size_t ca = calls[bi].allele1, cb = calls[bi].allele2;
                        if (oracle_called_only) {
                            // Only the called pair was aligned; every other al[] entry is a
                            // default-constructed zero, which the search below would happily report
                            // as a perfect best pair. Emit the called edits and nothing else.
                            const long ce = (ca < A && cb < A)
                                ? static_cast<long>(al[ca][0].edits + al[cb][1].edits)
                                : -1;
                            const long ce_sw = (ca < A && cb < A)
                                ? static_cast<long>(al[ca][1].edits + al[cb][0].edits)
                                : -1;
                            const long cbest = (ce < 0) ? -1 : std::min(ce, ce_sw);
                            orc << chain[bi].index << '\t'
                                << (chain[bi].kind == BlockKind::Bubble ? "bubble"
                                    : chain[bi].kind == BlockKind::Flank ? "flank" : "backbone")
                                << '\t' << chain[bi].bubble_id << '\t' << A << '\t'
                                << "edit_distance\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t"
                                << ca << '\t' << cb << "\tNA\tNA\tNA\tNA\tNA\tNA\t";
                            if (cbest >= 0) orc << cbest << "\tNA\n"; else orc << "NA\tNA\n";
                            orc.flush();
                            // The named pairs still have to be priced here. This branch returns
                            // early, and the emit below it never ran under --oracle-called-only --
                            // which is exactly the combination the pricing pass uses, so the table
                            // came out with a header and no rows.
                            if (!oracle_pairs.empty()) {
                                for (const auto& pr : oracle_pairs) {
                                    opf << chain[bi].index << '\t' << pr[0] << '\t' << pr[1] << '\t';
                                    if (pr[0] >= A || pr[1] >= A) {
                                        opf << "NA\tNA\tNA\tNA\tNA\tNA\tNA\n"; continue;
                                    }
                                    const long e0 = static_cast<long>(al[pr[0]][0].edits + al[pr[1]][1].edits);
                                    const long e1 = static_cast<long>(al[pr[0]][1].edits + al[pr[1]][0].edits);
                                    const int sw = (e1 < e0) ? 1 : 0;
                                    const std::size_t x = sw ? pr[1] : pr[0], y = sw ? pr[0] : pr[1];
                                    opf << ident(x, 0) << '\t' << ident(y, 1) << '\t'
                                        << dl[x][0] << '\t' << dl[y][1] << '\t' << std::min(e0, e1)
                                        << '\t' << ((ident(x, 0) + ident(y, 1)) / 2.0)
                                        // The cache is trimmed here by construction, so the certified
                                        // optimum is not visible and no excess can be honest.
                                        << "\tNA\n";
                                }
                                opf.flush();
                            }
                            continue;
                        }
                        for (int crit = 0; crit < 3; ++crit) {
                            double best = std::numeric_limits<double>::infinity();
                            std::size_t ba = 0, bb = 0;
                            int bswap = 0;
                            std::size_t tied = 0;
                            const double cval = (ca < A && cb < A) ? score(ca, cb, crit).v
                                                                   : std::numeric_limits<double>::infinity();
                            std::size_t rank = 1;
                            for (std::size_t a = 0; a < A; ++a) {
                                for (std::size_t b = a; b < A; ++b) {
                                    const PairScore s = score(a, b, crit);
                                    if (s.v < best - 1e-12) { best = s.v; ba = a; bb = b; bswap = s.swap; tied = 1; }
                                    else if (s.v <= best + 1e-12) ++tied;
                                    if (s.v < cval - 1e-12) ++rank;
                                }
                            }
                            const std::size_t b1 = bswap ? bb : ba, b2 = bswap ? ba : bb;
                            const PairScore cs = (ca < A && cb < A) ? score(ca, cb, crit) : PairScore{0, 0};
                            const std::size_t c1 = cs.swap ? cb : ca, c2 = cs.swap ? ca : cb;
                            orc << chain[bi].index << '\t'
                                << (chain[bi].kind == BlockKind::Bubble ? "bubble"
                                    : chain[bi].kind == BlockKind::Flank ? "flank" : "backbone") << '\t'
                                << chain[bi].bubble_id << '\t' << A << '\t' << names[crit] << '\t'
                                << ba << '\t' << bb << '\t'
                                << ident(b1, 0) << '\t' << ident(b2, 1) << '\t'
                                << dl[b1][0] << '\t' << dl[b2][1] << '\t'
                                << (al[b1][0].edits + al[b2][1].edits) << '\t'
                                << ((ident(b1, 0) + ident(b2, 1)) / 2.0) << '\t' << tied << '\t';
                            if (ca < A && cb < A) {
                                orc << ca << '\t' << cb << '\t' << rank << '\t'
                                    << ident(c1, 0) << '\t' << ident(c2, 1) << '\t'
                                    << dl[c1][0] << '\t' << dl[c2][1] << '\t'
                                    << ((ident(c1, 0) + ident(c2, 1)) / 2.0) << '\t';
                                // Only the edit_distance row's `best` minimises EDITS, so only there
                                // is (called - best) the certified edit excess. Writing a number on
                                // the other two rows would invite a later script to aggregate all
                                // three, or to read a length-optimal pair's edit count as the
                                // certified optimum.
                                if (crit == 0) {
                                    const long ce = static_cast<long>(al[c1][0].edits + al[c2][1].edits);
                                    const long be = static_cast<long>(al[b1][0].edits + al[b2][1].edits);
                                    orc << ce << '\t' << (ce - be) << '\n';
                                } else {
                                    orc << "NA\tNA\n";
                                }
                            } else {
                                orc << ca << '\t' << cb << "\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n";
                            }
                        }
                        orc.flush();
                        // Named pairs, priced from the same cache and against the same truth. Written
                        // to their own file rather than as extra rows of oracle.tsv, whose schema is
                        // one row per criterion and is consumed by existing scripts.
                        if (!oracle_pairs.empty()) {
                            std::ofstream& of = opf;
                            // best_* on the edit_distance criterion, recomputed here so the excess is
                            // against the certified optimum and not against whatever ran last.
                            //
                            // Under --oracle-called-only the cache holds only the called pair and the
                            // alleles named here, so the minimum over it is NOT the certified optimum
                            // over all A alleles -- it is a minimum over a handful. Reporting that
                            // difference as excess would understate it, silently, by however much the
                            // real optimum beats the named set. Excess is NA there; total_edits is
                            // still exact, and the caller can subtract a best_total_edits obtained
                            // from a full run.
                            const bool excess_computable = !oracle_called_only;
                            std::size_t eb_a = 0, eb_b = 0; long eb = -1;
                            for (std::size_t a = 0; a < A; ++a) {
                                if (oracle_called_only && !want_allele[a] && a != ca0 && a != cb0) continue;
                                for (std::size_t b = a; b < A; ++b) {
                                    if (oracle_called_only && !want_allele[b] && b != ca0 && b != cb0) continue;
                                    const long e0 = static_cast<long>(al[a][0].edits + al[b][1].edits);
                                    const long e1 = static_cast<long>(al[a][1].edits + al[b][0].edits);
                                    const long e = std::min(e0, e1);
                                    if (eb < 0 || e < eb) { eb = e; eb_a = a; eb_b = b; }
                                }
                            }
                            (void)eb_a; (void)eb_b;
                            for (const auto& pr : oracle_pairs) {
                                of << chain[bi].index << '\t' << pr[0] << '\t' << pr[1] << '\t';
                                if (pr[0] >= A || pr[1] >= A) { of << "NA\tNA\tNA\tNA\tNA\tNA\tNA\n"; continue; }
                                // The better of the two haplotype assignments, by edits, matching how
                                // the oracle scores the called pair.
                                const long e0 = static_cast<long>(al[pr[0]][0].edits + al[pr[1]][1].edits);
                                const long e1 = static_cast<long>(al[pr[0]][1].edits + al[pr[1]][0].edits);
                                const int sw = (e1 < e0) ? 1 : 0;
                                const std::size_t x = sw ? pr[1] : pr[0], y = sw ? pr[0] : pr[1];
                                const long tot = std::min(e0, e1);
                                of << ident(x, 0) << '\t' << ident(y, 1) << '\t'
                                   << dl[x][0] << '\t' << dl[y][1] << '\t' << tot << '\t'
                                   << ((ident(x, 0) + ident(y, 1)) / 2.0) << '\t'
                                   << ((excess_computable && eb >= 0) ? std::to_string(tot - eb)
                                                                      : std::string("NA")) << '\n';
                            }
                            of.flush();
                            if (!of) throw std::runtime_error("genotype: write failed for " + pp);
                        }
                    }
                    if (!orc) throw std::runtime_error("genotype: write failed for " + op);
                    log.wrote({op});
                    if (!oracle_pairs.empty()) log.wrote({pp});
                }

                if (!blocks.empty() && !blocks[0].allele_seq.empty()) {
                    const std::string acc_path = out_prefix + ".accuracy.tsv";
                    std::ofstream acc(acc_path);
                    acc << "block_index\tblock_kind\tbubble_id\tn_alleles\trepresentable\texact\t"
                           "dbp\tbest_dbp\tidentity\tbest_identity\toracle_rank\tid_h1\tid_h2\trank_h1\trank_h2"
                           "\tlenrank_h1\tlenrank_h2"
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
                        if (!tav1[bi] && !tav2[bi]) continue;
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
                        // Rank by LENGTH as well as by identity. At a tandem array the two disagree
                        // sharply: alleles differing by several copies are still ~99% identical,
                        // because identity is dominated by the shared repeat unit. An identity oracle
                        // therefore reports success on a call whose copy number is badly wrong, which
                        // is precisely the quantity a CNV caller exists to get right.
                        long lenrank_h0 = -1, lenrank_h1 = -1;
                        if (blocks[bi].allele_seq.size() <= 1024) {
                            std::vector<std::unordered_set<std::uint64_t>> asets;
                            asets.reserve(blocks[bi].allele_seq.size());
                            for (const std::string& s : blocks[bi].allele_seq) asets.push_back(sync_set(s));
                            long lrank_h[2] = {-1, -1};
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
                                // Where the called allele ranks by |length - truth length|.
                                {
                                    std::vector<std::pair<long, std::size_t>> byd;
                                    byd.reserve(blocks[bi].allele_seq.size());
                                    for (std::size_t ai = 0; ai < blocks[bi].allele_seq.size(); ++ai) {
                                        byd.emplace_back(std::labs(
                                            static_cast<long>(blocks[bi].allele_seq[ai].size()) -
                                            static_cast<long>(truth.size())), ai);
                                    }
                                    std::sort(byd.begin(), byd.end());
                                    for (std::size_t r = 0; r < byd.size(); ++r) {
                                        if (byd[r].second == cv[h ^ swap]) {
                                            lrank_h[h] = static_cast<long>(r) + 1;
                                            break;
                                        }
                                    }
                                }
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
                                lenrank_h0 = lrank_h[0]; lenrank_h1 = lrank_h[1];
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
                            << lenrank_h0 << '\t' << lenrank_h1 << '\t'
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
            log.wrote({out_prefix + ".reads.depth.tsv", out_prefix + ".genotypes.tsv"});
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
