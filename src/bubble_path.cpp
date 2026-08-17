#include "panvar/bubble_path.hpp"

#include <algorithm>
#include <unordered_set>

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

    // Every later endpoint, not just the nearest. The rule is "the interval with the most interior
    // steps wins", but taking only the first endpoint after each start made that unreachable: where a
    // boundary recurs -- a tandem array, a duplication, any path revisiting a boundary -- the enclosing
    // traversal ends at a LATER occurrence, and the nearest one closes a short interval containing
    // almost none of the snarl. Bounded per start, since only the widest few can ever win.
    constexpr std::size_t kMaxEndsPerStart = 64;
    for (const std::size_t start_pos : start_positions) {
        const auto first = std::upper_bound(end_positions.begin(), end_positions.end(), start_pos);
        std::size_t examined = 0;
        for (auto it = first; it != end_positions.end() && examined < kMaxEndsPerStart; ++it, ++examined) {
            BubblePathInterval candidate;
            candidate.left = start_pos;
            candidate.right = *it;
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
    const BubblePathInterval& interval,
    std::vector<std::size_t>* out_step_indices) {

    if (interval.left >= path.steps.size() || interval.right >= path.steps.size() || interval.left > interval.right) {
        if (out_step_indices != nullptr) out_step_indices->clear();
        return {};
    }

    // The index each returned step came from in `path.steps`, carried alongside rather than derived by
    // the caller: the reverse branch below may flip the walk a second time, so `left + j` and
    // `right - j` are both wrong in one of the four combinations. A caller that needs to know WHICH
    // occurrence of a repeated node it is looking at cannot reconstruct this safely from outside.
    std::vector<PathStep> out;
    std::vector<std::size_t> idx;
    out.reserve(interval.right - interval.left + 1);
    idx.reserve(interval.right - interval.left + 1);
    if (interval.source_to_sink) {
        for (std::size_t i = interval.left; i <= interval.right; ++i) {
            out.push_back(path.steps[i]);
            idx.push_back(i);
        }
        if (out_step_indices != nullptr) *out_step_indices = std::move(idx);
        return out;
    }

    for (std::size_t i = interval.right + 1; i > interval.left; --i) {
        const PathStep& step = path.steps[i - 1];
        out.push_back(PathStep{step.node_id, !step.reverse});
        idx.push_back(i - 1);
    }
    if (!out.empty() && (out.front().node_id != bubble.source || out.back().node_id != bubble.sink)) {
        std::reverse(out.begin(), out.end());
        std::reverse(idx.begin(), idx.end());
        for (auto& step : out) {
            step.reverse = !step.reverse;
        }
    }
    if (out_step_indices != nullptr) *out_step_indices = std::move(idx);
    return out;
}

// Steps of `path` across `bubble` (canonical source->sink). Falls back to an empty interior
// ([source, sink]) for paths that cross with no inside node (a pure deletion / short side of an
// insertion), which the inside-node-only interval finder would otherwise drop.
std::optional<std::vector<PathStep>> bubble_steps(
    const PathRecord& path, const BubblePathIndex& index, const Bubble& bubble,
    BubblePathInterval* used_interval, std::vector<std::size_t>* used_step_indices) {
    const auto iv = find_best_bubble_path_interval(index, bubble);
    if (iv.has_value()) {
        std::vector<PathStep> s = canonical_bubble_path_steps(path, bubble, *iv, used_step_indices);
        if (!s.empty()) {
            if (used_interval != nullptr) *used_interval = *iv;
            return s;
        }
    }
    const auto si = index.positions.find(bubble.source);
    const auto ki = index.positions.find(bubble.sink);
    if (si == index.positions.end() || ki == index.positions.end()) return std::nullopt;
    const auto record = [&](std::size_t l, std::size_t r, bool s2s) {
        if (used_step_indices != nullptr) *used_step_indices = s2s ? std::vector<std::size_t>{l, r}
                                                                   : std::vector<std::size_t>{r, l};
        if (used_interval == nullptr) return;
        used_interval->left = l;
        used_interval->right = r;
        used_interval->inside_count = 0;
        used_interval->source_to_sink = s2s;
    };
    const std::unordered_set<std::size_t> sink_pos(ki->second.begin(), ki->second.end());
    for (const std::size_t p : si->second) {                 // forward: source then sink
        if (sink_pos.count(p + 1)) {
            record(p, p + 1, true);
            return std::vector<PathStep>{ path.steps[p], path.steps[p + 1] };
        }
    }
    const std::unordered_set<std::size_t> src_pos(si->second.begin(), si->second.end());
    for (const std::size_t p : ki->second) {                 // reverse: sink then source -> flip
        if (src_pos.count(p + 1)) {
            record(p, p + 1, false);
            return std::vector<PathStep>{
                PathStep{ path.steps[p + 1].node_id, !path.steps[p + 1].reverse },
                PathStep{ path.steps[p].node_id, !path.steps[p].reverse } };
        }
    }
    return std::nullopt;
}

std::optional<std::vector<PathStep>> interval_steps(
    const PathRecord& path,
    const BubblePathIndex& index,
    const std::string& from_node,
    const std::string& to_node) {

    const auto fi = index.positions.find(from_node);
    const auto ti = index.positions.find(to_node);
    if (fi == index.positions.end() || ti == index.positions.end()) return std::nullopt;

    // Shortest interval spanning the two boundaries, in either direction.
    bool found = false;
    bool forward = true;
    std::size_t best_left = 0;
    std::size_t best_right = 0;
    auto consider = [&](std::size_t l, std::size_t r, bool fwd) {
        if (r <= l) return;
        if (!found || (r - l) < (best_right - best_left)) {
            found = true;
            forward = fwd;
            best_left = l;
            best_right = r;
        }
    };
    for (const std::size_t a : fi->second) {
        for (const std::size_t b : ti->second) {
            if (b > a) consider(a, b, true);
        }
    }
    if (!found) {
        for (const std::size_t b : ti->second) {
            for (const std::size_t a : fi->second) {
                if (a > b) consider(b, a, false);
            }
        }
    }
    if (!found) return std::nullopt;

    std::vector<PathStep> out;
    out.reserve(best_right - best_left + 1);
    if (forward) {
        for (std::size_t i = best_left; i <= best_right; ++i) out.push_back(path.steps[i]);
    } else {
        for (std::size_t i = best_right + 1; i > best_left; --i) {
            const PathStep& s = path.steps[i - 1];
            out.push_back(PathStep{s.node_id, !s.reverse});
        }
    }
    return out;
}

std::optional<std::vector<PathStep>> interval_interior_steps(
    const PathRecord& path,
    const BubblePathIndex& index,
    const std::string& from_node,
    const std::string& to_node) {

    auto full = interval_steps(path, index, from_node, to_node);
    if (!full.has_value() || full->size() <= 2) return std::nullopt;
    return std::vector<PathStep>(full->begin() + 1, full->end() - 1);
}

std::optional<std::vector<PathStep>> flank_steps(
    const PathRecord& path,
    const BubblePathIndex& index,
    const std::string& boundary_node,
    bool leading) {

    const auto it = index.positions.find(boundary_node);
    if (it == index.positions.end() || it->second.empty()) return std::nullopt;

    // A path may run through the graph in either direction; the boundary node's own orientation says
    // which. When it is reversed, the reference-upstream flank sits AFTER the boundary in path order.
    const std::size_t pos = leading ? it->second.front() : it->second.back();
    const bool rev = path.steps[pos].reverse;
    const bool take_prefix = leading != rev;

    std::vector<PathStep> out;
    if (take_prefix) {
        if (pos == 0) return std::nullopt;
        out.assign(path.steps.begin(), path.steps.begin() + static_cast<long>(pos));
    } else {
        if (pos + 1 >= path.steps.size()) return std::nullopt;
        out.assign(path.steps.begin() + static_cast<long>(pos) + 1, path.steps.end());
    }
    if (out.empty()) return std::nullopt;
    if (rev) {
        std::reverse(out.begin(), out.end());
        for (PathStep& st : out) st.reverse = !st.reverse;
    }
    return out;
}

} // namespace panvar
