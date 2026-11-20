#pragma once

#include "panvar/bubble_path.hpp"
#include "panvar/bubbles.hpp"
#include "panvar/gfa.hpp"

#include <cstddef>
#include <vector>

namespace panvar {

// The sites call reads in place of each site of `bubbles`: empty where the site is read as it is. A site some
// path, the reference included, crosses more than once lies in copies the graph collapsed. It is read one
// reference crossing at a time, each path from the crossing matched to it (Bubble::windows), and the stretch
// between two matched crossings where a path holds other copies than the reference is a site of its own
// (Bubble::between_copies), listed under the site before it. `ref_idx` indexes the reference in `graph.paths`.
std::vector<std::vector<Bubble>> copy_sites(const Graph& graph, const std::vector<BubblePathIndex>& path_indexes,
                                            const std::vector<Bubble>& bubbles, std::size_t ref_idx);

} // namespace panvar
