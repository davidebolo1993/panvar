#include "panvar/bubble_path.hpp"

#include <algorithm>

namespace panvar {
namespace {

std::vector<std::size_t> collect_inside_positions(
    const BubblePathIndex& index,
    const Bubble& bubble) {

    std::vector<std::size_t> inside;
    inside.reserve(bubble.inside.size());
    for (const auto& node_id : bubble.inside) {
        const auto it = index.positions.find(node_id);
        if (it == index.positions.end()) {
            continue;
        }
        inside.insert(inside.end(), it->second.begin(), it->second.end());
    }
    std::sort(inside.begin(), inside.end());
    inside.erase(std::unique(inside.begin(), inside.end()), inside.end());
    return inside;
}

std::size_t count_inside_between(
    const std::vector<std::size_t>& inside_positions,
    std::size_t left,
    std::size_t right) {

    if (inside_positions.empty() || left >= right) {
        return 0;
    }
    const auto begin_it = std::upper_bound(inside_positions.begin(), inside_positions.end(), left);
    const auto end_it = std::lower_bound(inside_positions.begin(), inside_positions.end(), right);
    return static_cast<std::size_t>(std::distance(begin_it, end_it));
}

bool better_interval(const BubblePathInterval& a, const BubblePathInterval& b) {
    if (a.inside_count != b.inside_count) {
        return a.inside_count > b.inside_count;
    }
    const std::size_t a_span = a.right - a.left;
    const std::size_t b_span = b.right - b.left;
    if (a_span != b_span) {
        return a_span < b_span;
    }
    return a.left < b.left;
}

void evaluate_direction(
    const std::vector<std::size_t>& start_positions,
    const std::vector<std::size_t>& end_positions,
    bool source_to_sink,
    const std::vector<std::size_t>& inside_positions,
    std::optional<BubblePathInterval>& best) {

    for (const std::size_t start_pos : start_positions) {
        const auto end_it = std::upper_bound(end_positions.begin(), end_positions.end(), start_pos);
        if (end_it == end_positions.end()) {
            continue;
        }
        BubblePathInterval candidate;
        candidate.left = start_pos;
        candidate.right = *end_it;
        candidate.inside_count = count_inside_between(inside_positions, candidate.left, candidate.right);
        candidate.source_to_sink = source_to_sink;
        if (candidate.inside_count == 0) {
            continue;
        }
        if (!best.has_value() || better_interval(candidate, *best)) {
            best = candidate;
        }
    }
}

} // namespace

BubblePathIndex build_bubble_path_index(const PathRecord& path) {
    BubblePathIndex out;
    out.positions.reserve(path.steps.size());
    for (std::size_t i = 0; i < path.steps.size(); ++i) {
        out.positions[path.steps[i].node_id].push_back(i);
    }
    return out;
}

std::optional<BubblePathInterval> find_best_bubble_path_interval(
    const BubblePathIndex& index,
    const Bubble& bubble) {

    const auto src_it = index.positions.find(bubble.source);
    const auto sink_it = index.positions.find(bubble.sink);
    if (src_it == index.positions.end() || sink_it == index.positions.end()) {
        return std::nullopt;
    }
    const auto inside_positions = collect_inside_positions(index, bubble);
    if (inside_positions.empty()) {
        return std::nullopt;
    }

    std::optional<BubblePathInterval> best;
    evaluate_direction(src_it->second, sink_it->second, true, inside_positions, best);
    evaluate_direction(sink_it->second, src_it->second, false, inside_positions, best);
    return best;
}

std::vector<PathStep> canonical_bubble_path_steps(
    const PathRecord& path,
    const Bubble& bubble,
    const BubblePathInterval& interval) {

    if (interval.left >= path.steps.size() || interval.right >= path.steps.size() || interval.left > interval.right) {
        return {};
    }

    std::vector<PathStep> out;
    out.reserve(interval.right - interval.left + 1);
    if (interval.source_to_sink) {
        for (std::size_t i = interval.left; i <= interval.right; ++i) {
            out.push_back(path.steps[i]);
        }
        return out;
    }

    for (std::size_t i = interval.right + 1; i > interval.left; --i) {
        const PathStep& step = path.steps[i - 1];
        out.push_back(PathStep{step.node_id, !step.reverse});
    }
    if (!out.empty() && (out.front().node_id != bubble.source || out.back().node_id != bubble.sink)) {
        std::reverse(out.begin(), out.end());
        for (auto& step : out) {
            step.reverse = !step.reverse;
        }
    }
    return out;
}

} // namespace panvar
