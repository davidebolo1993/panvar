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
    bool no_flip = false;
    bool quiet = false;
};

struct RefineSummary {
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
