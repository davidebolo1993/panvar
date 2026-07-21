#pragma once

// Reference-guided node ordering for a GFA, internal so the pipeline does not depend on odgi. Restores
// the invariant the rest of panvar relies on -- numeric node id == reference order -- in three steps:
// flip nodes the reference traverses in reverse, sort into reference-guided topological order (each
// bubble's interior contiguous after its source), renumber 1..N. Idempotent on an already-sorted graph.

#include "panvar/gfa_io.hpp"

#include <string>
#include <unordered_map>

namespace panvar {

struct GraphSortOptions {
    std::string reference_path; // path name or unique case-insensitive substring (required)
    bool flip = true;           // reorient nodes the reference traverses in reverse
    bool renumber = true;       // assign 1..N ids in sorted order
};

struct GraphSortResult {
    std::unordered_map<std::string, std::string> id_remap; // old id -> new id (when renumber)
    std::string resolved_reference; // the path name actually used
    std::size_t flipped_nodes = 0;
    std::size_t reference_steps = 0;
};

// Flip + reference-guided sort + (optional) renumber, in place on `model`.
// Throws if the reference path cannot be resolved uniquely.
GraphSortResult sort_graph_reference(GfaModel& model, const GraphSortOptions& options);

} // namespace panvar
