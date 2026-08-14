#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace panvar {

// POA-realigns the interiors of selected bubbles on a panphorte graph, so graph-builder artifacts
// (e.g. one indel split into a spurious INS+DEL) collapse into a clean call. Folded REP nodes are held
// fixed and only their flanks realigned; bubbles with an unfolded copy-number revisit are skipped.
// Emits panphorte's output family, so the result drops straight into `call`/`describe`/`inspect`.
struct RefineOptions {
    std::string gfa_path;                    // input panphorte normalized.sorted.gfa
    std::string bubbles_csv_in;              // input bubbles CSV (panphorte.bubbles.csv)
    std::string reference_path;              // reference path name or unique case-insensitive substring
    std::string out_prefix;                  // output prefix
    std::string gtf_path;                    // optional GTF -> <prefix>.bandage_genes.csv
    std::size_t max_poa_bp = 5000;           // primary cost guard: skip a residual segment longer than this
    std::size_t max_walks = 500;             // skip a residual segment with more distinct walks than this
    std::size_t min_bubbles = 1;             // only rebuild regions fusing >= this many bubbles
    std::vector<std::string> only_bubble_ids; // targeted mode: rebuild only these bubble ids (empty = auto)
    // Interior-span filter for the re-snarled call-ready CSV. Was BubbleCallOptions' own default,
    // applied silently to an input that may have been produced with a different one.
    std::size_t resnarl_min_variant_bp = 50;
    // A region that some paths traverse only partly: retaining their old nodes also retains the old
    // EDGES between them, so the pre-refinement topology survives beside the refined one and a walk can
    // still take it. Sequence losslessness cannot detect that, because every path still spells the same
    // bases. Skipping such a region is the conservative default.
    // Estimated abPOA work, longest x total bases over the DISTINCT sequences, in DP cells.
    // 0 disables. An independent budget rather than a function of the other two bounds, which would be
    // redundant given they already cap the longest sequence and the sequence count.
    std::size_t max_poa_work = 0;
    bool partial_path_policy_skip = true;
    bool no_flip = false;
    bool quiet = false;
};

struct RefineSummary {
    bool wrote_gene_annotation = false;  // <prefix>.bandage_genes.csv, when --gtf projected genes
    std::size_t regions_rebuilt = 0;
    std::size_t regions_skipped = 0;
    std::size_t nodes_added = 0;
    std::size_t nodes_removed = 0;
    // Interior nodes kept back from deletion because a path that does not span both anchors still
    // walks them. Non-zero means the naive removal would have corrupted the graph.
    std::size_t nodes_retained_referenced = 0;
    std::size_t bubbles_after = 0;
};

// Run refine end to end, writing <prefix>.normalized.sorted.gfa + <prefix>.bubbles.csv +
// <prefix>.bandage_nodes.csv (+ <prefix>.bandage_genes.csv when --gtf). Returns a summary.
RefineSummary refine_graph(const RefineOptions& options);

} // namespace panvar
