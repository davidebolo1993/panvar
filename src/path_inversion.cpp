#include "panvar/path_inversion.hpp"

#include <algorithm>
#include <map>
#include <tuple>
#include <unordered_map>
#include <cstdint>

namespace panvar {

PathInversionScan find_path_inversions(const Graph& graph, const PathRecord& reference, std::size_t min_bp) {
    PathInversionScan out;
    out.span.assign(graph.paths.size(), {SIZE_MAX, 0});
    if (min_bp == 0) return out;
    const auto len_of = [&](const std::string& id) -> std::size_t {
        const auto it = graph.nodes.find(id);
        return it == graph.nodes.end() ? 0 : it->second.sequence.size();
    };
    struct URef { std::size_t start = 0, len = 0; bool reverse = false; };
    std::unordered_map<std::string, std::size_t> ref_count;
    for (const PathStep& st : reference.steps) ++ref_count[st.node_id];
    std::unordered_map<std::string, URef> uref;
    std::size_t at = 0;
    for (const PathStep& st : reference.steps) {
        const std::size_t len = len_of(st.node_id);
        if (ref_count[st.node_id] == 1) uref[st.node_id] = URef{at, len, st.reverse};
        at += len;
    }
    const std::size_t ref_total = at;
    struct Run {
        bool fwd = true;
        std::size_t bp = 0, lo = SIZE_MAX, hi = 0, falls = 0, pairs = 0;
        std::size_t first_start = 0, last_start = 0;   // reference starts of the first and last steps it walks
        // Same for its outer runs of its own (not absorbed), preferring those of min_bp or more: these do not
        // depend on which way the haplotype is stored.
        std::size_t outer_first = 0, outer_last = 0;
        bool outer_long = false;
        std::string lo_node, hi_node;
    };
    struct Cand { std::size_t pi, left, a, b, right; std::string left_node, right_node; };
    std::vector<Cand> cands;
    for (std::size_t pi = 0; pi < graph.paths.size(); ++pi) {
        if (&graph.paths[pi] == &reference) { out.span[pi] = {0, ref_total == 0 ? 0 : ref_total - 1}; continue; }
        struct Hit { std::size_t start, len; bool same; const std::string* node; };
        std::vector<Hit> hits;
        std::unordered_map<std::string, std::size_t> visits;
        for (const PathStep& st : graph.paths[pi].steps) if (uref.count(st.node_id)) ++visits[st.node_id];
        for (const PathStep& st : graph.paths[pi].steps) {
            const auto it = uref.find(st.node_id);
            if (it == uref.end() || visits[st.node_id] != 1) continue;
            hits.push_back({it->second.start, it->second.len, st.reverse == it->second.reverse, &it->first});
        }
        if (hits.empty()) continue;
        for (const Hit& h : hits) {
            out.span[pi].first = std::min(out.span[pi].first, h.start);
            out.span[pi].second = std::max(out.span[pi].second, h.start + h.len - 1);
        }
        // Runs of hits walked the same way, in walk order. A run under min_bp between two others is absorbed
        // into the one before it: its bases do not move that run's reference range, but its steps count in
        // its falling order. A short run at either end is kept, as an end's flank.
        const auto runs_of = [&](const std::vector<Hit>& hs) {
            std::vector<Run> runs;
            for (const Hit& h : hs) {
                if (runs.empty() || runs.back().fwd != h.same) {
                    runs.push_back(Run{h.same});
                    runs.back().first_start = h.start;
                } else {
                    Run& r = runs.back();
                    ++r.pairs; if (h.start < r.last_start) ++r.falls;
                }
                Run& r = runs.back();
                r.bp += h.len;
                r.last_start = h.start;
                if (h.start < r.lo) { r.lo = h.start; r.lo_node = *h.node; }
                if (h.start + h.len - 1 >= r.hi) { r.hi = h.start + h.len - 1; r.hi_node = *h.node; }
            }
            // The step pair joining a run to the one it extends is tested like any other.
            const auto extend = [](Run& m, const Run& r) {
                m.pairs += r.pairs + 1; m.falls += r.falls + (r.first_start < m.last_start ? 1 : 0);
                m.last_start = r.last_start;
            };
            std::vector<Run> merged;
            for (std::size_t i = 0; i < runs.size(); ++i) {
                const Run& r = runs[i];
                if (r.bp < min_bp && i > 0 && i + 1 < runs.size()) { extend(merged.back(), r); continue; }
                if (!merged.empty() && merged.back().fwd == r.fwd) {
                    Run& m = merged.back();
                    extend(m, r);
                    m.bp += r.bp;
                    if (r.lo < m.lo) { m.lo = r.lo; m.lo_node = r.lo_node; }
                    if (r.hi > m.hi) { m.hi = r.hi; m.hi_node = r.hi_node; }
                    if (r.bp >= min_bp && !m.outer_long) { m.outer_first = r.first_start; m.outer_long = true; }
                    if (r.bp >= min_bp || !m.outer_long) m.outer_last = r.last_start;
                } else {
                    merged.push_back(r);
                    Run& m = merged.back();
                    m.outer_first = r.first_start; m.outer_last = r.last_start; m.outer_long = r.bp >= min_bp;
                }
            }
            return merged;
        };
        std::vector<Run> merged = runs_of(hits);
        // A haplotype stored on the other strand is read the reference's way: reversed, its steps turned.
        // Its strand comes from where its ends lie on the reference, not from how they are walked, which an
        // inversion near an end, a short inverted end node or one longer than both flanks would mislead. Its
        // ends are the outer steps of its outermost runs of min_bp or more, short runs absorbed into them not
        // counted, so a short node out of place near an end does not decide it; with fewer than two such runs,
        // its first and last steps decide.
        std::size_t first_long = SIZE_MAX, last_long = SIZE_MAX;
        for (std::size_t k = 0; k < merged.size(); ++k)
            if (merged[k].bp >= min_bp) { if (first_long == SIZE_MAX) first_long = k; last_long = k; }
        const bool reversed = first_long != last_long
                                  ? merged[last_long].outer_last < merged[first_long].outer_first
                                  : hits.back().start < hits.front().start;
        if (reversed) {
            std::reverse(hits.begin(), hits.end());
            for (Hit& h : hits) h.same = !h.same;
            merged = runs_of(hits);
        }
        // Only an end run can be under min_bp: an inversion near the end of a haplotype keeps its short flank.
        for (std::size_t k = 1; k + 1 < merged.size(); ++k) {
            const Run& L = merged[k - 1]; const Run& M = merged[k]; const Run& R = merged[k + 1];
            if (M.fwd || !L.fwd || !R.fwd || M.bp < min_bp) continue;
            if (M.pairs > 0 && 10 * M.falls < 8 * M.pairs) continue;   // not walked in falling order
            if (!(L.hi < M.lo && M.hi < R.lo)) continue;              // flanks must lie outside it
            cands.push_back(Cand{pi, L.hi, M.lo, M.hi, R.lo, L.hi_node, R.lo_node});
        }
    }
    std::sort(cands.begin(), cands.end(), [](const Cand& x, const Cand& y) {
        return std::tie(x.a, x.b, x.pi) < std::tie(y.a, y.b, y.pi);
    });
    std::vector<std::vector<Cand>> groups;
    for (const Cand& c : cands) {
        bool placed = false;
        for (auto& g : groups) {
            const Cand& r = g.front();
            const std::size_t ov_lo = std::max(r.a, c.a), ov_hi = std::min(r.b, c.b);
            const std::size_t ov = ov_hi >= ov_lo ? ov_hi - ov_lo + 1 : 0;
            if (10 * ov >= 9 * std::max(r.b - r.a + 1, c.b - c.a + 1)) { g.push_back(c); placed = true; break; }
        }
        if (!placed) groups.push_back({c});
    }
    for (const auto& g : groups) {
        std::map<std::tuple<std::size_t, std::size_t, std::size_t, std::size_t>, std::size_t> votes;
        for (const Cand& c : g) ++votes[{c.left, c.a, c.b, c.right}];
        const auto best = std::max_element(votes.begin(), votes.end(),
            [](const auto& x, const auto& y) { return x.second < y.second; })->first;
        PathInversion inv;
        std::tie(inv.left_end, inv.inv_lo, inv.inv_hi, inv.right_start) = best;
        for (const Cand& c : g) {
            if (std::tie(c.left, c.a, c.b, c.right) == best) { inv.left_node = c.left_node; inv.right_node = c.right_node; }
            inv.carriers.push_back(c.pi);
        }
        std::sort(inv.carriers.begin(), inv.carriers.end());
        inv.carriers.erase(std::unique(inv.carriers.begin(), inv.carriers.end()), inv.carriers.end());
        out.inversions.push_back(std::move(inv));
    }
    return out;
}

} // namespace panvar
