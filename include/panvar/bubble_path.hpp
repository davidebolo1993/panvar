#pragma once

#include "panvar/bubbles.hpp"
#include "panvar/gfa.hpp"

#include <cstddef>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

namespace panvar {

struct BubblePathIndex {
    std::unordered_map<std::string, std::vector<std::size_t>> positions;
};

struct BubblePathInterval {
    std::size_t left = 0;
    std::size_t right = 0;
    std::size_t inside_count = 0;
    bool source_to_sink = true;
};

BubblePathIndex build_bubble_path_index(const PathRecord& path);
std::optional<BubblePathInterval> find_best_bubble_path_interval(
    const BubblePathIndex& index,
    const Bubble& bubble);
// `out_step_indices`, when given, receives the index in `path.steps` each returned step came from.
// This is the only safe way to know which OCCURRENCE of a repeated node a walk position refers to.
std::vector<PathStep> canonical_bubble_path_steps(
    const PathRecord& path,
    const Bubble& bubble,
    const BubblePathInterval& interval,
    std::vector<std::size_t>* out_step_indices = nullptr);

// Steps of `path` across `bubble` (canonical source->sink). Falls back to an empty interior
// ([source, sink]) for paths that cross with no inside node (a pure deletion / short side of an
// insertion), which the inside-node-only interval finder above would otherwise drop. Prefer this
// over `canonical_bubble_path_steps` whenever every allele of a bubble must be represented.
// `used_interval`, when non-null, receives the interval the steps came from -- including for the
// empty-interior fallback, where `inside_count` is 0. Callers that report the interval as metadata
// need it for every allele, not only the ones with an inside node.
std::optional<std::vector<PathStep>> bubble_steps(
    const PathRecord& path,
    const BubblePathIndex& index,
    const Bubble& bubble,
    BubblePathInterval* used_interval = nullptr,
    std::vector<std::size_t>* used_step_indices = nullptr);

// Steps of `path` from `from_node` to `to_node`, oriented from->to. Unlike `bubble_steps` this keys
// on the two boundary nodes alone, so it works for a stretch with no declared interior — the
// backbone between two consecutive bubbles, which carries the sub-threshold variation that
// `--min-variant-bp` kept out of the bubble list. Picks the shortest spanning interval; returns
// nullopt when the path does not traverse both boundaries in either orientation.
std::optional<std::vector<PathStep>> interval_steps(
    const PathRecord& path,
    const BubblePathIndex& index,
    const std::string& from_node,
    const std::string& to_node);

// The same interval with both boundary nodes removed. A bubble owns its source and sink -- they are
// its anchors, and `bubble_steps` returns them -- so a backbone stretch that also included them would
// put the same sequence in two blocks. That is not cosmetic: blocks are meant to TILE the haplotype,
// and marker confinement drops anything varying in more than one block, so a shared boundary node
// silently destroys the markers of both blocks that touch it. Returns nullopt when nothing lies
// strictly between the boundaries.
std::optional<std::vector<PathStep>> interval_interior_steps(
    const PathRecord& path,
    const BubblePathIndex& index,
    const std::string& from_node,
    const std::string& to_node);

// Steps of `path` on one side of a single boundary node: the stretch before the first bubble or
// after the last. Orientation is taken from how the path traverses the boundary itself, so a path
// running through the graph backwards still yields a reference-oriented flank. `leading` selects the
// reference-upstream side.
std::optional<std::vector<PathStep>> flank_steps(
    const PathRecord& path,
    const BubblePathIndex& index,
    const std::string& boundary_node,
    bool leading);

} // namespace panvar
