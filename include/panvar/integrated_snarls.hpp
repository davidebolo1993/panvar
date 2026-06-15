#pragma once

// Faithful port of vg's IntegratedSnarlFinder cactus decomposition, restricted to emitting
// the TOP-LEVEL (source,sink) boundary pairs that `bubble` consumes. The algorithm is
// vendored near-verbatim (integrated_snarl_finder.cpp); only the graph backing is swapped
// for a tiny in-process handle shim, and the begin/end traversal callbacks collect the
// outermost snarls. See snarls.cpp for the caller.

#include "panvar/gfa_io.hpp"

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace panvar {

// Minimal bidirected graph description: nodes are ranks 0..N-1; edges carry GFA-style
// orientation (true = the link uses the reverse strand of that endpoint).
struct SnarlGraphInput {
    std::vector<std::string> node_ids;  // index = node rank
    std::vector<std::size_t> node_len;  // bp length per rank (for rooting by sequence)
    struct Edge {
        std::size_t from_rank;
        bool from_rev;
        std::size_t to_rank;
        bool to_rev;
    };
    std::vector<Edge> edges;
};

// Top-level snarl boundary pairs (node ids), matching vg's `vg snarls` top-level set.
std::vector<std::pair<std::string, std::string>> find_top_level_snarls_cactus(
    const SnarlGraphInput& in);

// Build the cactus finder input from an editable GFA model (ranks = node_order index,
// edge orientation from L-line orients). Use after an internal sort so ids are stable.
SnarlGraphInput snarl_input_from_model(const GfaModel& model);

} // namespace panvar
