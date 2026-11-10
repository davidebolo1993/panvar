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

    // The interval with the most interior steps wins, shortest span then leftmost breaking a tie. No
    // enumeration of start/end pairs is needed (quadratic, and any bound clips the answer where a boundary
    // recurs): the interior count is monotone under containment, so the widest interval attains the
    // maximum, and the shortest interval attaining it is the tightest window around those interior steps,
    // the last start before the first interior step and the first end after the last. Two binary searches.
    if (start_positions.empty() || end_positions.empty()) {
        return;
    }
    const std::size_t widest_left = start_positions.front();
    const std::size_t widest_right = end_positions.back();
    if (widest_right <= widest_left) {
        return;
    }
    const std::size_t inside_count = count_inside_between(inside_positions, widest_left, widest_right);
    if (inside_count == 0) {
        return;
    }

    // inside_count > 0 makes every dereference below valid: there is at least one interior position
    // strictly inside the widest interval, so each bound is bracketed by an endpoint on its own side.
    const std::size_t first_inside =
        *std::upper_bound(inside_positions.begin(), inside_positions.end(), widest_left);
    const std::size_t last_inside =
        *(std::lower_bound(inside_positions.begin(), inside_positions.end(), widest_right) - 1);

    BubblePathInterval candidate;
    candidate.left =
        *(std::lower_bound(start_positions.begin(), start_positions.end(), first_inside) - 1);
    candidate.right =
        *std::upper_bound(end_positions.begin(), end_positions.end(), last_inside);
    candidate.inside_count = inside_count;
    candidate.source_to_sink = source_to_sink;
    if (!best.has_value() || better_interval(candidate, *best)) {
        best = candidate;
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

BubblePathIndex build_bubble_path_index(const PathRecord& path, const std::unordered_set<std::string>& keep) {
    BubblePathIndex out;
    for (std::size_t i = 0; i < path.steps.size(); ++i)
        if (keep.count(path.steps[i].node_id)) out.positions[path.steps[i].node_id].push_back(i);
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

    // The index each returned step came from in `path.steps`, carried alongside, so a caller knows which
    // occurrence of a repeated node it is looking at. A reversed interval has its source at `right` and
    // its sink at `left`, so reading it backwards gives source to sink.
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
    if (out_step_indices != nullptr) *out_step_indices = std::move(idx);
    return out;
}

// Steps of `path` across `bubble` (canonical source->sink), with an empty-interior fallback; see the
// header.
std::optional<std::vector<PathStep>> bubble_steps(
    const PathRecord& path, const BubblePathIndex& index, const Bubble& bubble,
    BubblePathInterval* used_interval, std::vector<std::size_t>* used_step_indices, bool inward_only) {
    // With `inward_only`, a window must leave its first boundary and reach its last through the inward sides:
    // a walk into the site that turns back and reaches the far boundary from outside does not cross it.
    const auto inward = [&](const BubblePathInterval& w) {
        if (!inward_only || bubble.source == bubble.sink) return true;
        if (w.left >= path.steps.size() || w.right >= path.steps.size()) return false;
        const bool rl = path.steps[w.left].reverse, rr = path.steps[w.right].reverse;
        return w.source_to_sink ? rl == bubble.source_reverse && rr == bubble.sink_reverse
                                : rl != bubble.sink_reverse && rr != bubble.source_reverse;
    };
    auto iv = find_best_bubble_path_interval(index, bubble);
    if (iv.has_value() && !inward(*iv)) {
        // The widest window runs past a turn at an end; keep the crossings inside it instead. Per direction,
        // one window from the first crossing with an interior to the last, and the better of the two.
        iv.reset();
        const auto inside_positions = collect_inside_positions(index, bubble);
        const auto crossings = bubble_crossings(path, index, bubble, false);
        for (const bool fwd : {true, false}) {
            std::optional<BubblePathInterval> w;
            for (const BubblePathInterval& c : crossings) {
                if (c.source_to_sink != fwd || c.right <= c.left + 1) continue;
                if (!w) w = c;
                w->right = c.right;
            }
            if (!w) continue;
            w->inside_count = count_inside_between(inside_positions, w->left, w->right);
            if (w->inside_count > 0 && (!iv || better_interval(*w, *iv))) iv = w;
        }
    }
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
    // A crossing with nothing between the boundaries is an adjacent pair: source then sink, or sink
    // then source (flipped). With `inward_only` it must be read through both inward sides as well (see
    // the header). A path crossing this way more than once gets one window over all its crossings, as
    // with an interior: `1 3 1 3` is two copies of the deletion, not the deletion.
    constexpr std::size_t kNone = static_cast<std::size_t>(-1);
    std::size_t f_first = kNone, f_last = kNone, r_first = kNone, r_last = kNone;
    const std::unordered_set<std::size_t> sink_pos(ki->second.begin(), ki->second.end());
    for (const std::size_t p : si->second) {                 // forward: source then sink
        if (!sink_pos.count(p + 1)) continue;
        if (inward_only && (path.steps[p].reverse != bubble.source_reverse ||
                            path.steps[p + 1].reverse != bubble.sink_reverse)) continue;
        if (f_first == kNone) f_first = p;
        f_last = p;
    }
    const std::unordered_set<std::size_t> src_pos(si->second.begin(), si->second.end());
    for (const std::size_t p : ki->second) {                 // reverse: sink then source -> flip
        if (!src_pos.count(p + 1)) continue;
        if (inward_only && (path.steps[p].reverse == bubble.sink_reverse ||
                            path.steps[p + 1].reverse == bubble.source_reverse)) continue;
        if (r_first == kNone) r_first = p;
        r_last = p;
    }
    BubblePathInterval ev;
    if (f_first != kNone) {
        ev.left = f_first;
        ev.right = f_last + 1;
        ev.source_to_sink = true;
    } else if (r_first != kNone) {
        ev.left = r_first;
        ev.right = r_last + 1;
        ev.source_to_sink = false;
    } else {
        return std::nullopt;
    }
    ev.inside_count = 0;
    std::vector<PathStep> s = canonical_bubble_path_steps(path, bubble, ev, used_step_indices);
    if (s.empty()) return std::nullopt;
    if (used_interval != nullptr) *used_interval = ev;
    return s;
}

// Every crossing, in path order; see the header for the rule. Callers that edit or spell each crossing need
// them apart: one window over two crossings would also take in the stretch between them.
std::vector<BubblePathInterval> bubble_crossings(const PathRecord& path, const BubblePathIndex& index,
                                                 const Bubble& bubble, bool fallback) {
    std::vector<BubblePathInterval> out;
    const auto si = index.positions.find(bubble.source);
    const auto ki = index.positions.find(bubble.sink);
    // A table without orientations gives none to test: its default would pick crossings at random.
    const bool oriented = bubble.orientation_known || !fallback;
    if (oriented && bubble.source != bubble.sink && si != index.positions.end() && ki != index.positions.end()) {
        std::vector<std::pair<std::size_t, bool>> hits;   // (step, is_source), in path order
        hits.reserve(si->second.size() + ki->second.size());
        for (const std::size_t p : si->second) hits.emplace_back(p, true);
        for (const std::size_t p : ki->second) hits.emplace_back(p, false);
        std::sort(hits.begin(), hits.end());
        for (std::size_t k = 0; k + 1 < hits.size(); ++k) {
            const auto [pa, a_src] = hits[k];
            const auto [pb, b_src] = hits[k + 1];
            if (a_src == b_src) continue;
            const bool ra = path.steps[pa].reverse, rb = path.steps[pb].reverse;
            const bool fwd = a_src && ra == bubble.source_reverse && rb == bubble.sink_reverse;
            const bool bwd = !a_src && ra != bubble.sink_reverse && rb != bubble.source_reverse;
            if (!fwd && !bwd) continue;
            BubblePathInterval iv;
            iv.left = pa;
            iv.right = pb;
            iv.inside_count = pb - pa - 1;
            iv.source_to_sink = fwd;
            out.push_back(iv);
        }
    }
    if (out.empty() && fallback) {
        BubblePathInterval used{};
        if (bubble_steps(path, index, bubble, &used).has_value()) out.push_back(used);
    }
    return out;
}

bool bubble_orientation_holds(const std::vector<PathRecord>& paths, const std::vector<BubblePathIndex>& indexes,
                              const Bubble& bubble) {
    if (!bubble.orientation_known) return false;
    for (std::size_t pi = 0; pi < paths.size(); ++pi)
        if (!bubble_crossings(paths[pi], indexes[pi], bubble, false).empty()) return true;
    return false;
}

} // namespace panvar
