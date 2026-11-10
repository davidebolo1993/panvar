#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace panvar {

// POA-realigns the interiors of selected bubbles on a panphorte graph, so graph-builder artifacts
// (e.g. one indel split into a spurious INS+DEL) collapse into a clean call. Folded repeats (the nodes
// panphorte's unit graph table names, and any node with a self-loop) are held fixed and only their flanks
// realigned; a bubble with an unfolded revisit of --min-dup-bp or more is skipped. Emits panphorte's
// output family, repeat tables included, so the result drops straight into `call`/`describe`/`inspect`.
struct RefineOptions {
    std::string gfa_path;                    // input panphorte normalized.sorted.gfa
    std::string bubbles_csv_in;              // input bubbles CSV (panphorte.bubbles.csv)
    std::string reference_path;              // reference path, as -r resolves it
    std::string out_prefix;                  // output prefix
    std::string gtf_path;                    // optional GTF -> <prefix>.bandage_genes.csv
    // Where panphorte's unit graph table is read from: <pfx>.panphorte.unit_graph.tsv, when present (empty =
    // none). Written back under out_prefix.
    std::string panphorte_prefix;
    std::size_t max_poa_bp = 5000;           // skip a region with a residual sequence longer than this
    std::size_t max_walks = 500;             // skip a residual segment with more distinct sequences than this
    std::size_t min_bubbles = 1;             // only rebuild regions fusing >= this many bubbles
    std::vector<std::string> only_bubble_ids; // targeted mode: rebuild only these bubble ids (empty = auto)
    // Interior-span filter for the re-snarled call-ready CSV (bubble's --min-variant-bp), set here so it
    // can match the threshold the input was built with.
    std::size_t resnarl_min_variant_bp = 50;
    // The re-snarl's fusion by placement and (off by default) proximity fusion, as in panphorte (see there).
    std::size_t resnarl_merge_nearby_bp = 0;
    bool resnarl_fuse_ambiguous = true;
    // Estimated abPOA work per residual segment, longest x total bases over the distinct sequences, in DP
    // cells; 0 disables. A budget of its own: a bound derived from the other two would be redundant.
    std::size_t max_poa_work = 0;
    // A crossing whose repeat visits of non-unit nodes (every visit after a node's first) carry fewer
    // bases than this is realigned rather than treated as an unfolded duplication. 0, the default, makes
    // any revisit skip the region: short loops are copy-number signal for `call`, and realigning them
    // turns it into plain indels.
    std::size_t min_dup_bp = 0;
    bool no_flip = false;
    bool quiet = false;
};

struct RefineSummary {
    bool wrote_gene_annotation = false;  // <prefix>.bandage_genes.csv, when --gtf projected genes
    std::size_t regions_rebuilt = 0;
    std::size_t regions_skipped = 0;
    std::size_t nodes_added = 0;
    std::size_t nodes_removed = 0;
    std::size_t bubbles_after = 0;
    std::size_t twins_merged = 0;   // nodes merged into an identical twin
    std::size_t twin_paths_walking_several = 0;   // paths that walked two nodes of one merged group
    std::size_t rep_tables_read = 0;              // panphorte repeat tables read
    std::vector<std::string> rep_tables_written;  // ... and written back under the output prefix
    std::vector<std::string> stale_removed;       // optional outputs of an earlier run at the prefix, removed
};

// Run refine end to end, writing <prefix>.normalized.sorted.gfa + <prefix>.bubbles.csv +
// <prefix>.bandage_nodes.csv + the report table (+ <prefix>.bandage_genes.csv when --gtf, + the
// panphorte repeat tables it read, renumbered). Returns a summary.
RefineSummary refine_graph(const RefineOptions& options);

} // namespace panvar
