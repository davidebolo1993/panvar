#pragma once

// Port of vg's IntegratedSnarlFinder cactus decomposition (integrated_snarls.cpp), emitting the top-level
// (source, sink) boundary pairs `bubble` consumes. The algorithm is vendored near-verbatim; only the graph
// backing is a small in-process handle shim, and the traversal callbacks collect the outermost snarls.

#include "panvar/gfa_io.hpp"

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace panvar {

// Minimal bidirected graph description: nodes are ranks 0..N-1; edges carry GFA-style orientation (true:
// the link uses the reverse strand of that endpoint).
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
    // First and last node rank of every path, for recognising the snarl that closes a root cycle (see
    // find_top_level_snarls_cactus). Empty: no such recognition.
    std::vector<std::pair<std::size_t, std::size_t>> path_ends;
    // Per node rank, whether paths walk it around: one path twice, or some forward and others backward (each
    // path read on the strand its end nodes' ranks give). A cycle none of whose nodes is walked around is a
    // site's alleles, not a repeat or an inversion, and does not root the decomposition (see
    // traverse_decomposition). Empty: no node is walked around, so no cycle roots it; vg's rooting, where
    // every cycle may, is every entry 1.
    std::vector<char> walked_twice;
};

// Which side of each boundary the snarl's interior hangs off: true = the node's end side.
struct SnarlInward {
    bool source_end = true;
    bool sink_end = false;
};

// Top-level snarl boundary pairs (node ids), as vg's `vg snarls` top-level set, with two differences. A cycle
// no path walks around does not root the decomposition (see `walked_twice`), so where vg roots on such a
// cycle the snarls differ. And where a component's decomposition roots in a cycle (an inverted duplication, a locus-wide rearrangement), one
// top-level snarl closes the cycle and wraps everything outside its chain, the locus ends included, so every
// real site there is nested in it. A snarl in a cycle-rooted component holding an end (first or last step)
// of more than half of that component's paths may be such a wrapper: it is returned with the snarls inside
// it, recursively, and bubble's overlap pass keeps the smaller sites. `inward`, when given, receives each
// pair's inward sides in the same order, since the pair alone carries no orientation.
std::vector<std::pair<std::string, std::string>> find_top_level_snarls_cactus(
    const SnarlGraphInput& in, std::vector<SnarlInward>* inward = nullptr);

// Build the cactus finder input from an editable GFA model (ranks = node_order index, edge orientation from
// L-line orients). Use after an internal sort so ids are stable.
SnarlGraphInput snarl_input_from_model(const GfaModel& model);

} // namespace panvar
