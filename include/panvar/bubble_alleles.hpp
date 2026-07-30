#pragma once

#include "panvar/bubble_path.hpp"
#include "panvar/bubbles.hpp"
#include "panvar/gfa.hpp"

#include <cstddef>
#include <string>
#include <unordered_set>
#include <vector>

namespace panvar {

// One distinct allele of a bubble: the walk itself plus the haplotypes that carry it.
// Alleles are grouped by canonical-walk signature, so two haplotypes share an allele only when
// they traverse the same nodes in the same orientations.
struct BubbleAllele {
    std::size_t id = 0;                     // index within the set, in first-encounter order
    std::string signature;                  // grouping key ("12+,13-,14+")
    std::vector<PathStep> steps;            // the walk, canonical source->sink
    std::vector<std::string> members;       // panel path names carrying it
};

struct BubbleAlleleSet {
    std::size_t bubble_id = 0;
    std::vector<BubbleAllele> alleles;
    std::unordered_set<std::string> traversing;   // path names that cross the bubble at all
    std::vector<PathStep> reference_steps;        // empty when the reference does not cross
    std::string reference_signature;
    bool has_reference = false;
};

// Enumerate the distinct alleles of `bubble` across every path in `graph`. `path_indexes` is
// parallel to `graph.paths`. Uses `bubble_steps`, so pure-deletion alleles (no interior node) are
// represented rather than dropped. Allele order is the order haplotypes are first seen, which keeps
// it stable for a given graph.
BubbleAlleleSet enumerate_bubble_alleles(
    const Graph& graph,
    const std::vector<BubblePathIndex>& path_indexes,
    const Bubble& bubble,
    const std::string& reference_path_name);

} // namespace panvar
