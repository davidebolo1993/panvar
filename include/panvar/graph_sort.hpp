#pragma once

// Reference-guided node ORDERING for a GFA, done internally so the pipeline does not
// depend on `odgi sort`/`odgi flip`. Operates on the editable GfaModel (so both `bubble`
// and `panphorte` can reuse it) and restores the invariant the rest of panvar relies on:
// numeric node id == reference order. Three steps:
//   1. flip   - reverse-complement nodes the reference path traverses in reverse so the
//               reference reads all-forward (optional; a no-op when already forward).
//   2. sort   - reference-guided topological order (reference backbone in path order, each
//               bubble's interior contiguous right after its source).
//   3. renumber - assign ids 1..N in sorted order (rewriting S/L/P/W).
// Idempotent on a graph that is already sorted+flipped along the reference.

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
