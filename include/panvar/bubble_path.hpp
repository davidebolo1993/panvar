#pragma once

#include "panvar/bubbles.hpp"
#include "panvar/gfa.hpp"

#include <cstddef>
#include <optional>
#include <string>
#include <unordered_map>
#include <unordered_set>
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
// Only the positions of `keep`'s nodes: the boundaries and interiors of the sites a run reads. A whole-path
// index over every path of a panel holds every step of every walk, far more than a run needs.
BubblePathIndex build_bubble_path_index(const PathRecord& path, const std::unordered_set<std::string>& keep);
std::optional<BubblePathInterval> find_best_bubble_path_interval(
    const BubblePathIndex& index,
    const Bubble& bubble);
// `out_step_indices`, when given, receives the index in `path.steps` each returned step came from: the only
// safe way to know which occurrence of a repeated node a walk position refers to.
std::vector<PathStep> canonical_bubble_path_steps(
    const PathRecord& path,
    const Bubble& bubble,
    const BubblePathInterval& interval,
    std::vector<std::size_t>* out_step_indices = nullptr);

// Steps of `path` across `bubble` (canonical source->sink). Falls back to an empty interior ([source, sink])
// for paths that cross with no inside node (a pure deletion, the short side of an insertion), which the
// inside-node interval finder above would drop; prefer this whenever every allele must be represented.
// `used_interval`, when non-null, receives the interval the steps came from, the fallback's included. A path
// crossing with an empty interior more than once gets one window over all those crossings, as with an
// interior. `inward_only` admits a window only when it leaves its first boundary and reaches its last through
// the inward sides as `source_reverse`/`sink_reverse` give them, so a walk past the outer sides (`3+ 1+`
// beside a `1+..3+` site) is not a deletion, and one that enters, turns back and reaches the far boundary from
// outside (`1+ 2+ 2- 1- 4+ 3-`) is no crossing; such a window is narrowed to the crossings inside it (those
// `bubble_crossings` finds), so a real crossing elsewhere on the path is kept. Use it only when the orientation is known to be right
// (bubble's own sites): a site table's columns can disagree with the walks, and there those walks and a
// crossing look the same.
std::optional<std::vector<PathStep>> bubble_steps(
    const PathRecord& path,
    const BubblePathIndex& index,
    const Bubble& bubble,
    BubblePathInterval* used_interval = nullptr,
    std::vector<std::size_t>* used_step_indices = nullptr,
    bool inward_only = false);

// Every crossing of `bubble` by `path`, in path order: consecutive boundary occurrences, source then sink or
// sink then source, that leave the first and reach the second through the sides facing the interior (as the
// reference walks them, or mirrored for a path walking the site backwards). The stretch between two
// crossings leaves through an outer side and is not one; `bubble_steps` instead takes one window spanning
// both. With `fallback`, a path none of whose boundary pairs is a crossing gets that one window: for a site
// table whose orientation columns do not describe the walks. A table without orientation columns gives every
// path that window, with no crossings tried. Where the columns are right, that window is a
// path touching the site from outside (`1- 2+ 3-`, or `3+ 1+`), not an allele, so callers decide per site
// with `bubble_orientation_holds`. Turn an interval into steps with `canonical_bubble_path_steps`.
std::vector<BubblePathInterval> bubble_crossings(
    const PathRecord& path,
    const BubblePathIndex& index,
    const Bubble& bubble,
    bool fallback = true);

// Whether the site table's boundary orientations describe this graph: the table has orientation columns and
// some path crosses through both inward sides as they give them. Then `bubble_crossings` is called without
// its fallback; otherwise the fallback is the reading. A table without the columns reads every boundary
// '+', and one path that happens to cross that way must not make every other path strict.
bool bubble_orientation_holds(const std::vector<PathRecord>& paths, const std::vector<BubblePathIndex>& indexes,
                              const Bubble& bubble);

} // namespace panvar
