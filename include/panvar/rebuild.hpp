#pragma once

#include <cstddef>
#include <string>

namespace panvar {

// `rebuild`: re-induce a fragmented low-complexity locus graph before bubble decomposition. A degree
// gate decides whether the graph is pathological; if so the locus is rebuilt by progressive graph
// generation, driving minigraph as a library. Ours are the two ends: the haplotype order fed to it
// (k-mer richness, since the seed becomes the backbone) and the emitted GFA, which carries per-haplotype
// P lines and preserves link orientation. Healthy graphs pass through untouched.
struct RebuildOptions {
    std::string gfa_path;            // input GFA (the pggb graph for the locus)
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
    // ---- acceptance contract ----
    // minigraph augments variation ABOVE --min-var, so sub-threshold differences (SNPs, short indels)
    // are collapsed by construction and a recovered walk is never byte-identical to its haplotype.
    // The contract is therefore structural plus a threshold, not losslessness: every path must come
    // back, spelling what it spelled before to within these bounds, or the rebuild is rejected and the
    // original graph is passed through unchanged.
    double min_recovered_identity = 0.98;  // per-path recovered-walk identity
    double min_matched_cover = 0.95;       // per-path matching bases / haplotype length
    std::string reference_path;            // if set, must be recovered and meet the bounds; else reject
    bool allow_loss = false;               // accept a rebuild that fails the contract (records why)
    std::string audit_path;                // per-path audit TSV (empty = <out>.rebuild_audit.tsv)
    bool quiet = false;
};

struct RebuildSummary {
    std::size_t raw_nodes = 0;
    std::size_t raw_hubs = 0;        // #nodes with degree >= hub_degree
    std::size_t raw_maxdeg = 0;      // largest single-handle degree (not the two ends pooled)
    std::size_t raw_selfloops = 0;   // nodes with a self-loop: a folded tandem array, not pathology
    double raw_density = 0.0;        // nodes per kb of the longest haplotype span
    bool pathological = false;
    bool ran = false;                // did we actually rebuild (pathological or forced)?
    std::size_t haplotypes = 0;
    std::string seed;                // richest haplotype (Criteria B)
    std::size_t out_nodes = 0;
    std::size_t out_edges = 0;
    std::size_t out_maxdeg = 0;
    std::size_t out_hubs = 0;        // #nodes with degree >= hub_degree in the rebuilt graph
    std::size_t out_selfloops = 0;
    std::size_t paths_recovered = 0; // haplotypes whose walk mg_map resolved
    // (qe-qs)/len -- the chain's OUTER envelope. It says where the alignment starts and ends, not how
    // much of the query is actually aligned, so a chain with a large internal gap still reads high.
    double mean_query_cover = 0.0;
    // Matching query bases over query length: the fraction genuinely aligned, so an internal gap
    // shows up here where the envelope hides it. This is the number an acceptance threshold wants.
    double mean_matched_cover = 0.0;
    double min_matched_cover = 0.0;
    // Matches over alignment-block length, within the aligned region.
    double mean_chain_identity = 0.0;
    // Identity of the RECOVERED WALK re-spelled from the rebuilt graph against the original haplotype.
    // The end-to-end check: it is what a caller actually gets back, and it catches an error anywhere in
    // mapping, chain choice or emission rather than only in the alignment.
    double mean_walk_identity = 0.0;
    double min_walk_identity = 0.0;
    std::size_t walk_identity_checked = 0;   // haplotypes the walk check ran on
    std::size_t walks_below_99 = 0;          // recovered at under 99% identity
    // ---- acceptance ----
    bool accepted = false;                   // did the rebuild satisfy the contract?
    std::string reject_reason;               // why not, when it did not
    std::size_t paths_failing = 0;           // paths below the identity/cover bounds
    std::size_t dangling_steps = 0;          // consecutive P steps with no emitted L
    bool audit_written = false;
};

RebuildSummary rebuild_graph(const RebuildOptions& options);

} // namespace panvar
