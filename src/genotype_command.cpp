#include "panvar/genotype_command.hpp"

#include "panvar/align.hpp"
#include "panvar/cli_utils.hpp"
#include "panvar/genotype_blocks.hpp"
#include "panvar/genotype_markers.hpp"
#include "panvar/genotype.hpp"
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
#include <fstream>
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
    // Every --probe-pair, in order. The first also redirects the truth_* columns, which is what the
    // single-probe form has always done; all of them appear in the probe table.
    std::vector<std::array<std::size_t, 3>> probe_pairs;
    bool oracle_called_only = false;
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
    if (oracle_called_only && !certified_oracle) {
        throw std::runtime_error("genotype: --oracle-called-only only means anything with "
                                 "--certified-oracle; on its own it reports nothing");
    }
    // A diagnostic that quietly does nothing is worse than one that refuses: the run looks like an
    // arm of the experiment and is actually the baseline, which is exactly how a null result gets
    // manufactured. Both combinations below were silently inert.
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

            ReadCounts rc = count_reads(read_paths, read_panel, options.threads);
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

            std::vector<std::string> hap_names;
            hap_names.reserve(panel_graph.paths.size());
            for (const PathRecord& p : panel_graph.paths) hap_names.push_back(p.name);
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
            std::vector<ProbePairResult> probe_rows;
            std::vector<BlockCall> calls =
                genotype_sample(chain, blocks, read_panel, rc, depth, hap_names, gopt, &gsum,
                                pa1.empty() ? nullptr : &pa1, pa2.empty() ? nullptr : &pa2,
                                evidence == "syncmer" ? nullptr : &cev, &probe_rows);

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
                        for (std::size_t a = 0; a < A; ++a) {
                            if (oracle_called_only && a != ca0 && a != cb0) continue;
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
                    }
                    if (!orc) throw std::runtime_error("genotype: write failed for " + op);
                    log.wrote({op});
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
