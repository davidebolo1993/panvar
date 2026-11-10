#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace panvar {

// `rebuild`: re-induce a fragmented low-complexity locus graph before bubble decomposition. A degree
// gate decides whether the graph is pathological; if so the locus is rebuilt by progressive graph
// generation, driving minigraph as a library. Ours are the two ends: the haplotype order fed to it
// (k-mer richness, since the seed becomes the backbone) and the emitted GFA, which carries per-haplotype
// P lines and preserves link orientation. Healthy graphs pass through untouched.
struct RebuildOptions {
    std::string gfa_path;            // input GFA
    std::string out_path;            // output GFA
    std::size_t kmer = 21;           // Criteria B: k for the k-mer richness metric (distinct k-mers
                                     // first, total k-mers breaking ties)
    std::size_t hub_degree = 50;     // Criteria A: a node with >= this many distinct neighbours is a hub
    std::size_t min_hubs = 10;       // Criteria A: >= this many hubs => pathological
    std::size_t min_var = 50;        // minigraph -L: minimum variant length to augment into the graph
    std::size_t min_align_len = 0;   // minigraph -l: chains shorter than this augment nothing. 0 = auto,
                                     // i.e. scaled to the locus (minigraph's default assumes chromosomes)
    std::string tmp_dir;             // parent for the per-haplotype FASTA scratch (empty = beside --out);
                                     // a dedicated subfolder is created under it and removed on exit
    std::size_t threads = 0;         // worker threads (0 = auto)
    bool force = false;              // rebuild even if the gate says healthy (testing / small inputs)
    bool base_align = false;         // minigraph -c: base-level alignment during graph generation
    bool clip_to_reference = false;  // with -r: drop haplotype sequence running past the reference's ends
    // ---- acceptance contract ----
    // minigraph augments variation above --min-var, so differences below it are collapsed and a walk need not
    // be byte-identical to its haplotype (the seed's is). The contract is structural plus a threshold: every
    // path comes back, spelling what it spelled before within these bounds, or the rebuild is rejected and
    // the original graph passed through.
    double min_recovered_identity = 0.98;  // per-path recovered-walk identity
    double min_matched_cover = 0.95;       // per-path matching bases / haplotype length
    std::string reference_path;            // if set, it seeds the graph and must come back exactly
    bool allow_loss = false;               // accept a rebuild that fails the contract (records why)
    // Fold copy-number modules (module_fold.hpp): the generator keeps a duplication's copies apart, so
    // without it a copy gain is a long insertion and nothing counts copies.
    bool fold_modules = true;
    double fold_identity = 0.95;           // identity of the repeat that makes a module (--fold-identity)
    std::size_t fold_min_copy_bp = 3000;   // shortest copy that makes a module (--fold-min-copy-bp)
    bool quiet = false;
    // Set by rebuild for its second pass: the graph the first pass generated, read instead of generating it
    // again, and the haplotypes whose copy gains and losses are written into their walks, which the first
    // pass left to a fold that did not take them.
    std::string generated_gfa;
    std::vector<std::string> copies_for;
};

struct RebuildSummary {
    std::size_t raw_nodes = 0;
    std::size_t raw_hubs = 0;        // #nodes with degree >= hub_degree
    std::size_t raw_maxdeg = 0;      // largest single-handle degree (not the two ends pooled)
    std::size_t raw_selfloops = 0;   // nodes with a self-loop: a folded tandem array, not pathology
    double raw_density = 0.0;        // nodes per kb of the longest haplotype span
    bool pathological = false;
    bool ran = false;                // did we actually rebuild (pathological or forced)?
    bool degenerate = false;         // a rebuild ran but recovered no variation; input passed through
    std::size_t haplotypes = 0;
    std::string seed;                // the reference if -r was given, else the richest haplotype
    std::size_t out_nodes = 0;
    std::size_t out_edges = 0;
    std::size_t out_maxdeg = 0;
    std::size_t out_hubs = 0;        // #nodes with degree >= hub_degree in the rebuilt graph
    std::size_t out_selfloops = 0;
    std::size_t paths_recovered = 0; // haplotypes whose walk mg_map resolved
    // (qe-qs)/len, the chain's outer envelope: where the alignment starts and ends, so a chain with a large
    // internal gap still reads high.
    double mean_query_cover = 0.0;
    // Matching query bases over query length, from the base-level alignment: an internal gap shows here.
    double mean_matched_cover = 0.0;
    double min_matched_cover = 0.0;
    // Matches over alignment-block length, within the aligned region.
    double mean_chain_identity = 0.0;
    // Identity of the recovered walk re-spelled from the rebuilt graph against the original haplotype: what a
    // caller gets back.
    double mean_walk_identity = 0.0;
    double min_walk_identity = 0.0;
    std::size_t walk_identity_checked = 0;   // haplotypes the walk check ran on
    std::size_t walks_below_99 = 0;          // recovered at under 99% identity
    // ---- acceptance ----
    bool accepted = false;                   // did the rebuild satisfy the contract?
    bool overridden = false;                 // emitted anyway under --allow-loss (contract not met)
    bool reference_exact = true;             // -r given: its walk spells it byte for byte (not
                                             // overridable by --allow-loss)
    std::string reject_reason;               // why not, when it did not
    std::size_t paths_failing = 0;           // paths below the identity/cover bounds
    std::size_t dangling_steps = 0;          // consecutive P steps with no emitted L
    std::size_t end_cuts = 0;                // segment ends cut where a haplotype stops
    std::size_t inversions_joined = 0;       // haplotypes whose large inversion was joined into one walk
    std::size_t gaps_kept = 0;               // assembly gaps emitted as N pieces on their haplotype's walk
    std::size_t tips_added = 0;              // pieces holding a haplotype's own start/end past the graph
    std::size_t clipped_haplotypes = 0;      // --clip-to-reference: haplotypes shortened
    std::size_t clipped_bp = 0;              //   ... and the bases dropped
    std::size_t pruned_segments = 0;         // segments (or cut pieces) no walk uses, not emitted
    std::size_t pruned_links = 0;            // links no walk uses, not emitted
    std::size_t indels_added = 0;            // insertions or deletions the generator left out, written into walks
    bool report_written = false;
    std::size_t modules_found = 0;           // copy-number modules found in the rebuilt graph
    std::size_t modules_folded = 0;          // ... folded (the others are reported with why not)
};

RebuildSummary rebuild_graph(const RebuildOptions& options);

// The report beside an output graph: its name without .gfa or .gfa.gz, plus .rebuild.report.tsv.
std::string rebuild_report_path(const std::string& out_path);

} // namespace panvar
