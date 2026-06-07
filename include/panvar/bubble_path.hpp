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
std::vector<PathStep> canonical_bubble_path_steps(
    const PathRecord& path,
    const Bubble& bubble,
    const BubblePathInterval& interval);

} // namespace panvar
