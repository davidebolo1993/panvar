#include "panvar/bubbles.hpp"

#include "panvar/bubble_path.hpp"

#include "panvar/cli_utils.hpp"

#include <algorithm>
#include <atomic>
#include <cstdlib>
#include <iostream>
#include <cstdint>
#include <fstream>
#include <functional>
#include <limits>
#include <optional>
#include <queue>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace panvar {
namespace {

// Compact node id <-> index mapping (indices key the per-endpoint dedup map).
struct PackedGraph {
    std::vector<std::string> node_ids;
    std::unordered_map<std::string, std::uint32_t> node_idx_of;
};

struct BubbleCandidateIdx {
    std::uint32_t source = 0;
    std::uint32_t sink = 0;
    std::vector<std::uint32_t> inside;
};

struct EndpointKey {
    std::uint32_t a = 0;
    std::uint32_t b = 0;

    bool operator==(const EndpointKey& other) const {
        return a == other.a && b == other.b;
    }
};

struct EndpointKeyHash {
    std::size_t operator()(const EndpointKey& key) const {
        const std::uint64_t packed = (static_cast<std::uint64_t>(key.a) << 32U) | key.b;
        return std::hash<std::uint64_t>{}(packed);
    }
};

using CandidateMap = std::unordered_map<EndpointKey, BubbleCandidateIdx, EndpointKeyHash>;

// The bubble-path index and the walk extraction live in bubble_path.{hpp,cpp} and are what every other
// module uses. bubbles.cpp used to carry a second, structurally identical index with its own interval
// search, and the two had drifted: the local one required at least one declared interior node between
// the boundaries, so a direct source->sink allele -- a pure deletion -- was invisible to bubble scoring
// while `bubble_steps` handled it correctly everywhere else. Aliasing rather than duplicating removes
// the class of divergence, not just this instance.
using PathIndex = BubblePathIndex;

// How many endpoint pairs fell back from the graph-derived interior to the path-derived one because
// the handle traversal hit its cap. The fallback is a silent downgrade of a stated contract -- the
// interior stops being "what the graph carries between the boundaries" and becomes "what this panel
// happened to walk" -- so the run must say when it happened. Atomic because discovery is threaded.
std::atomic<std::size_t> g_interior_traversal_truncated{0};

// Defined below, next to the rest of the interior machinery; declared here because merging needs to
// recompute a fused bubble's interior with exactly the rule discovery used.
std::unordered_set<std::uint32_t> interior_indices_for_pair(
    const Graph& graph,
    const PackedGraph& packed,
    const std::vector<PathIndex>& path_indexes,
    const std::string& source_id,
    const std::string& sink_id,
    std::uint32_t src_idx,
    std::uint32_t sink_idx);
bool reachable_interior(
    const Graph& graph,
    const PackedGraph& packed,
    std::uint64_t start_handle,
    std::uint32_t stop_idx,
    std::uint32_t origin_idx,
    bool backward,
    std::unordered_set<std::uint32_t>& out);

PackedGraph pack_graph(const Graph& graph) {
    PackedGraph packed;
    packed.node_ids.reserve(graph.nodes.size());
    for (const auto& [node_id, _node] : graph.nodes) {
        packed.node_ids.push_back(node_id);
    }
    std::sort(packed.node_ids.begin(), packed.node_ids.end());

    packed.node_idx_of.reserve(packed.node_ids.size() * 2);
    for (std::size_t i = 0; i < packed.node_ids.size(); ++i) {
        packed.node_idx_of[packed.node_ids[i]] = static_cast<std::uint32_t>(i);
    }
    return packed;
}

EndpointKey normalize_endpoint_key(std::uint32_t source, std::uint32_t sink) {
    if (source < sink) {
        return EndpointKey{source, sink};
    }
    return EndpointKey{sink, source};
}

// Keep one candidate per unordered endpoint pair, preferring the larger inside set.
bool add_or_replace_candidate(CandidateMap& dedup, BubbleCandidateIdx candidate) {
    const EndpointKey key = normalize_endpoint_key(candidate.source, candidate.sink);
    const auto it = dedup.find(key);
    if (it == dedup.end() || candidate.inside.size() > it->second.inside.size()) {
        dedup[key] = std::move(candidate);
        return true;
    }
    return false;
}

std::string endpoint_key(const std::string& a, const std::string& b) {
    if (a < b) {
        return a + "|" + b;
    }
    return b + "|" + a;
}

bool is_unsigned_decimal(const std::string& s) {
    if (s.empty()) {
        return false;
    }
    for (const char c : s) {
        if (c < '0' || c > '9') {
            return false;
        }
    }
    return true;
}

std::string_view strip_leading_zeros(std::string_view s) {
    std::size_t i = 0;
    while (i + 1 < s.size() && s[i] == '0') {
        ++i;
    }
    return s.substr(i);
}

bool node_id_less(const std::string& a, const std::string& b) {
    const bool a_num = is_unsigned_decimal(a);
    const bool b_num = is_unsigned_decimal(b);
    if (a_num && b_num) {
        const std::string_view av = strip_leading_zeros(a);
        const std::string_view bv = strip_leading_zeros(b);
        if (av.size() != bv.size()) {
            return av.size() < bv.size();
        }
        if (av != bv) {
            return av < bv;
        }
        return a < b;
    }
    if (a_num != b_num) {
        return a_num;
    }
    return a < b;
}

// Orient every bubble so `source` is the reference-LEFT boundary and `sink` the reference-RIGHT one.
//
// The cactus decomposition names the two boundaries of a snarl, but nothing in it fixes which is which:
// on a three-node graph whose reference is 1,2,3 it reports source=3, sink=1. That is not wrong as a
// snarl -- the pair is unordered -- but every consumer downstream reads them as an interval in
// reference order: `call` derives POS from the source, merging joins one bubble's sink to the next
// bubble's source, and a coordinate written from a reversed pair lands at the wrong end.
//
// A bubble the reference does not traverse has no reference order to take, and is left alone.
void orient_bubbles_to_reference(const Graph& graph, const std::string& reference_path,
                                 std::vector<Bubble>& bubbles) {
    if (reference_path.empty() || bubbles.empty()) return;
    std::size_t ref_idx = graph.paths.size();
    for (std::size_t i = 0; i < graph.paths.size(); ++i)
        if (graph.paths[i].name == reference_path) { ref_idx = i; break; }
    if (ref_idx == graph.paths.size()) {
        std::size_t hits = 0, last = 0;
        for (std::size_t i = 0; i < graph.paths.size(); ++i)
            if (graph.paths[i].name.find(reference_path) != std::string::npos) { ++hits; last = i; }
        if (hits != 1) return;                      // ambiguous or absent: no order to impose
        ref_idx = last;
    }
    // First occurrence of each node on the reference walk. A boundary the reference visits more than
    // once has no single coordinate; the first is the one `call` already anchors on.
    std::unordered_map<std::string, std::size_t> first_pos;
    const PathRecord& ref = graph.paths[ref_idx];
    for (std::size_t i = 0; i < ref.steps.size(); ++i)
        first_pos.emplace(ref.steps[i].node_id, i);
    for (Bubble& b : bubbles) {
        const auto s_it = first_pos.find(b.source);
        const auto k_it = first_pos.find(b.sink);
        if (s_it == first_pos.end() || k_it == first_pos.end()) continue;
        if (k_it->second < s_it->second) std::swap(b.source, b.sink);
        // Each boundary is a HANDLE: record which way the reference reads it, so a consumer knows
        // which side of the node faces into the bubble.
        const auto so = first_pos.find(b.source);
        const auto ko = first_pos.find(b.sink);
        if (so != first_pos.end()) b.source_reverse = ref.steps[so->second].reverse;
        if (ko != first_pos.end()) b.sink_reverse = ref.steps[ko->second].reverse;
    }
}

bool bubble_endpoint_less(const Bubble& a, const Bubble& b) {
    if (a.source != b.source) {
        return node_id_less(a.source, b.source);
    }
    if (a.sink != b.sink) {
        return node_id_less(a.sink, b.sink);
    }
    return a.inside.size() < b.inside.size();
}




// Reference coordinates of a bubble's boundaries: the bp offset at which each boundary node starts on
// the reference walk, and its length. Bubbles the reference does not span cannot be ordered against
// each other and are never merged.
struct RefSpan {
    std::size_t src_start = 0;   // bp offset of the source node's first base
    std::size_t snk_start = 0;   // bp offset of the sink node's first base
    std::size_t snk_len = 0;
    std::size_t src_len = 0;
    bool ok = false;
};

std::unordered_map<std::string, std::pair<std::size_t, std::size_t>> reference_offsets(
    const Graph& graph, const std::string& reference_path) {
    std::unordered_map<std::string, std::pair<std::size_t, std::size_t>> out;   // node -> (bp start, len)
    if (reference_path.empty()) return out;
    std::size_t ref_idx = graph.paths.size();
    for (std::size_t i = 0; i < graph.paths.size(); ++i)
        if (graph.paths[i].name == reference_path) { ref_idx = i; break; }
    if (ref_idx == graph.paths.size()) {
        std::size_t hits = 0, last = 0;
        for (std::size_t i = 0; i < graph.paths.size(); ++i)
            if (graph.paths[i].name.find(reference_path) != std::string::npos) { ++hits; last = i; }
        if (hits != 1) return out;
        ref_idx = last;
    }
    std::size_t off = 0;
    for (const PathStep& st : graph.paths[ref_idx].steps) {
        const auto nit = graph.nodes.find(st.node_id);
        const std::size_t len = nit == graph.nodes.end() ? 0 : nit->second.sequence.size();
        out.emplace(st.node_id, std::make_pair(off, len));   // first occurrence wins
        off += len;
    }
    return out;
}

// Merge bubbles that sit close together on the REFERENCE.
//
// The old version sorted by source node id and joined `current.sink` to `next.source` through an
// undirected shortest path, which assumed cactus order follows reference order. It does not: on the
// bundled synthetic graph the reference traverses 2..4 while cactus reports source=4, sink=2, and with
// a large threshold the fused bubble ended up containing nodes before its own source and after its own
// sink. It also seeded the distance with the whole length of the starting boundary node, so two
// bubbles sharing a boundary had a nominal gap of zero yet failed to merge when that node was long.
//
// Ordering and distance now both come from reference coordinates. The gap between two bubbles is the
// sequence strictly BETWEEN the facing boundaries -- after the left bubble's sink ends, before the
// right bubble's source begins -- which is zero for bubbles that abut or share a boundary. A fused
// bubble is validated: every interior node must lie inside the new span, or the merge is refused.
std::vector<Bubble> merge_nearby_bubbles(
    const Graph& graph,
    const PackedGraph& packed,
    const std::vector<PathIndex>& path_indexes,
    const std::vector<Bubble>& input,
    std::size_t max_bp,
    const std::string& reference_path) {

    if (max_bp == 0 || input.size() < 2) {
        return input;
    }
    const auto ref_off = reference_offsets(graph, reference_path);
    if (ref_off.empty()) return input;   // no reference to order against: merging would be guesswork

    auto span_of = [&](const Bubble& b) {
        RefSpan r;
        const auto s = ref_off.find(b.source);
        const auto k = ref_off.find(b.sink);
        if (s == ref_off.end() || k == ref_off.end()) return r;
        r.src_start = s->second.first;  r.src_len = s->second.second;
        r.snk_start = k->second.first;  r.snk_len = k->second.second;
        r.ok = r.src_start <= r.snk_start;      // boundaries are reference-ordered by now
        return r;
    };

    // Reference order, and only bubbles the reference actually spans take part.
    std::vector<Bubble> spanned, unspanned;
    for (const Bubble& b : input) (span_of(b).ok ? spanned : unspanned).push_back(b);
    if (std::getenv("PANVAR_MERGE_DEBUG"))
        std::cerr << "[merge] spanned=" << spanned.size() << " unspanned=" << unspanned.size() << '\n';
    std::sort(spanned.begin(), spanned.end(), [&](const Bubble& a, const Bubble& b) {
        const RefSpan ra = span_of(a), rb = span_of(b);
        if (ra.src_start != rb.src_start) return ra.src_start < rb.src_start;
        return ra.snk_start < rb.snk_start;
    });

    std::vector<Bubble> merged;
    merged.reserve(spanned.size());
    for (std::size_t i = 0; i < spanned.size(); ++i) {
        if (merged.empty()) { merged.push_back(spanned[i]); continue; }
        Bubble& cur = merged.back();
        const RefSpan rc = span_of(cur), rn = span_of(spanned[i]);
        const std::size_t cur_end = rc.snk_start + rc.snk_len;   // one past the left bubble's last base
        // Strictly between the facing boundaries; overlapping or abutting bubbles have a gap of 0.
        const std::size_t gap = rn.src_start > cur_end ? rn.src_start - cur_end : 0;
        if (std::getenv("PANVAR_MERGE_DEBUG"))
            std::cerr << "[merge] " << cur.source << ".." << cur.sink << " -> " << spanned[i].source
                      << ".." << spanned[i].sink << " gap=" << gap << " max=" << max_bp << '\n';
        if (gap > max_bp) { merged.push_back(spanned[i]); continue; }

        Bubble fused = cur;
        fused.sink = spanned[i].sink;
        fused.sink_reverse = spanned[i].sink_reverse;
        // The fused interior is recomputed from the graph between the new outer boundaries, not
        // assembled from the parts. Taking the union of the two interiors plus whatever the REFERENCE
        // carries across the connector describes only the reference route between them: an alternate
        // connector branch -- and any cycle riding on it, which is what --superbubbles must see -- is
        // in neither part's interior and is not on the reference, so it went missing exactly when the
        // fused bubble started to span it.
        std::unordered_set<std::string> inside_set(cur.inside.begin(), cur.inside.end());
        inside_set.insert(spanned[i].inside.begin(), spanned[i].inside.end());
        inside_set.insert(cur.sink);                 // old boundaries become interior
        inside_set.insert(spanned[i].source);
        {
            const auto s_it = packed.node_idx_of.find(fused.source);
            const auto k_it = packed.node_idx_of.find(fused.sink);
            if (s_it != packed.node_idx_of.end() && k_it != packed.node_idx_of.end()) {
                for (const std::uint32_t idx : interior_indices_for_pair(
                         graph, packed, path_indexes, fused.source, fused.sink,
                         s_it->second, k_it->second)) {
                    inside_set.insert(packed.node_ids[idx]);
                }
            }
        }

        const RefSpan rf = span_of(fused);
        bool contained = rf.ok;
        if (contained) {
            const std::size_t lo = rf.src_start, hi = rf.snk_start + rf.snk_len;
            for (const std::string& n : inside_set) {
                if (n == fused.source || n == fused.sink) continue;
                const auto it = ref_off.find(n);
                if (it == ref_off.end()) continue;   // off-reference alleles legitimately sit inside
                if (it->second.first < lo || it->second.first >= hi) { contained = false; break; }
            }
        }
        if (!contained) {
            if (std::getenv("PANVAR_MERGE_DEBUG"))
                std::cerr << "[merge]   refused: fused span would not contain every interior node\n";
            merged.push_back(spanned[i]); continue;
        }

        fused.inside.clear();
        fused.inside.reserve(inside_set.size());
        for (const std::string& n : inside_set)
            if (n != fused.source && n != fused.sink) fused.inside.push_back(n);
        std::sort(fused.inside.begin(), fused.inside.end());
        fused.path_support = 0; fused.min_inside_bp = 0; fused.max_inside_bp = 0;
        fused.long_path_support = 0; fused.inversion_signal = false;
        fused.distinct_alleles = 0; fused.ref_allele_support = 0;
        fused.alt_allele_support_max = 0; fused.alt_allele_support_min = 0;
        cur = std::move(fused);
    }
    for (Bubble& b : unspanned) merged.push_back(std::move(b));
    std::sort(merged.begin(), merged.end(), bubble_endpoint_less);
    return merged;
}

// ---- Path indexing + snarl JSONL import ----------------------------------------------

PathIndex build_path_index(const PathRecord& path) { return build_bubble_path_index(path); }

std::optional<std::string> extract_node_id_for_key(
    const std::string& line,
    const std::string& key,
    bool use_last_key) {

    const std::string token = "\"" + key + "\"";
    const std::size_t key_pos = use_last_key ? line.rfind(token) : line.find(token);
    if (key_pos == std::string::npos) {
        return std::nullopt;
    }
    const std::size_t node_key_pos = line.find("\"node_id\"", key_pos);
    if (node_key_pos == std::string::npos) {
        return std::nullopt;
    }
    const std::size_t colon_pos = line.find(':', node_key_pos);
    if (colon_pos == std::string::npos) {
        return std::nullopt;
    }
    const std::size_t q1 = line.find('"', colon_pos + 1);
    if (q1 == std::string::npos) {
        return std::nullopt;
    }
    const std::size_t q2 = line.find('"', q1 + 1);
    if (q2 == std::string::npos || q2 <= q1 + 1) {
        return std::nullopt;
    }
    return line.substr(q1 + 1, q2 - q1 - 1);
}

// Pick, for one snarl endpoint pair, the path interval with the most inside steps
// (shortest span as a tiebreak), in either source->sink or sink->source direction.
struct EndpointOnlyInterval {
    std::size_t left = 0;
    std::size_t right = 0;
    std::size_t inside_steps = 0;
    bool source_to_sink = true;
};

bool better_endpoint_interval(const EndpointOnlyInterval& a, const EndpointOnlyInterval& b) {
    if (a.inside_steps != b.inside_steps) {
        return a.inside_steps > b.inside_steps;
    }
    const std::size_t a_span = a.right - a.left;
    const std::size_t b_span = b.right - b.left;
    if (a_span != b_span) {
        return a_span < b_span;
    }
    return a.left < b.left;
}

void evaluate_endpoint_interval_direction(
    const std::vector<std::size_t>& start_positions,
    const std::vector<std::size_t>& end_positions,
    bool source_to_sink,
    std::optional<EndpointOnlyInterval>& best) {

    // Where a boundary recurs, the enclosing traversal ends at a LATER occurrence and the nearest one
    // closes an interval holding almost none of the snarl. This used to enumerate the endpoints after
    // each start, bounded at 64 to keep it affordable, which silently clipped any path revisiting a
    // boundary more often than that.
    //
    // The enumeration was never needed. This interval is ranked on `inside_steps`, which here is just
    // `right - left - 1` -- the span, not a count of interior nodes -- so "most interior steps" and
    // "widest" are the same requirement, and better_endpoint_interval's shortest-span tie-break is
    // unreachable because equal inside_steps means equal span. The widest pair is the earliest start
    // with the latest end, and it is unique: any other pair has a strictly smaller span.
    //
    // An ADJACENT pair (end_pos == start_pos + 1) is a real crossing with an empty interior, and it was
    // once rejected here. That looks harmless -- it contributes no interior nodes -- but this interval
    // is also what tells the graph-derived search which SIDE of each boundary faces inward. With every
    // stored path taking the direct source->sink route, no interval survived, the graph search was
    // never seeded, and an alternate allele that exists in the graph produced no bubble at all. It is
    // admitted here (right > left), and it wins only when there is nothing wider.
    if (start_positions.empty() || end_positions.empty()) {
        return;
    }
    EndpointOnlyInterval candidate;
    candidate.left = start_positions.front();
    candidate.right = end_positions.back();
    if (candidate.right <= candidate.left) {
        return;
    }
    candidate.inside_steps = candidate.right - candidate.left - 1;
    candidate.source_to_sink = source_to_sink;
    if (!best.has_value() || better_endpoint_interval(candidate, *best)) {
        best = candidate;
    }
}

std::optional<EndpointOnlyInterval> find_best_endpoint_interval(
    const PathIndex& index,
    const std::string& source,
    const std::string& sink) {

    const auto src_it = index.positions.find(source);
    const auto sink_it = index.positions.find(sink);
    if (src_it == index.positions.end() || sink_it == index.positions.end()) {
        return std::nullopt;
    }

    std::optional<EndpointOnlyInterval> best;
    evaluate_endpoint_interval_direction(src_it->second, sink_it->second, true, best);
    evaluate_endpoint_interval_direction(sink_it->second, src_it->second, false, best);
    return best;
}

// Parse top-level (source,sink) boundary pairs from a vg snarls JSONL (skip nested snarls).
std::vector<std::pair<std::string, std::string>> read_snarl_pairs_jsonl(const std::string& path) {
    std::ifstream snarls_in(path);
    if (!snarls_in) {
        throw std::runtime_error("Failed to open snarl JSONL input: " + path);
    }
    std::vector<std::pair<std::string, std::string>> pairs;
    std::string line;
    while (std::getline(snarls_in, line)) {
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        if (line.find("\"start\"") == std::string::npos || line.find("\"end\"") == std::string::npos) {
            continue;
        }
        if (line.find("\"parent\"") != std::string::npos) {
            continue; // top-level snarls only
        }
        const auto source_opt = extract_node_id_for_key(line, "start", true);
        const auto sink_opt = extract_node_id_for_key(line, "end", false);
        if (!source_opt.has_value() || !sink_opt.has_value()) {
            continue;
        }
        pairs.emplace_back(*source_opt, *sink_opt);
    }
    return pairs;
}

// Build bubble candidates from top-level (source,sink) pairs: for each pair, union the
// inside nodes seen on every path that crosses source->sink. Shared by the internal snarl
// finder and the --snarls-in override.

// The interior of one endpoint pair: everything the panel walks between the boundaries, unioned with
// everything the GRAPH carries between them. One function, so discovery and merging cannot drift --
// a fused bubble used to take the union of its parts' interiors plus the nodes the REFERENCE carries
// across the connector, which omits any alternate connector branch and any cycle riding on one.
std::unordered_set<std::uint32_t> interior_indices_for_pair(
    const Graph& graph,
    const PackedGraph& packed,
    const std::vector<PathIndex>& path_indexes,
    const std::string& source_id,
    const std::string& sink_id,
    std::uint32_t src_idx,
    std::uint32_t sink_idx) {

    std::unordered_set<std::uint32_t> inside_idx_set;
    std::optional<EndpointOnlyInterval> best_iv;
    std::size_t best_iv_path = 0;
    for (std::size_t p_idx = 0; p_idx < path_indexes.size(); ++p_idx) {
        const auto interval = find_best_endpoint_interval(path_indexes[p_idx], source_id, sink_id);
        if (!interval.has_value()) {
            continue;
        }
        if (!best_iv.has_value() || better_endpoint_interval(*interval, *best_iv)) {
            best_iv = interval;
            best_iv_path = p_idx;
        }
        const auto& steps = graph.paths[p_idx].steps;
        for (std::size_t i = interval->left + 1; i < interval->right; ++i) {
            const auto idx_it = packed.node_idx_of.find(steps[i].node_id);
            if (idx_it == packed.node_idx_of.end()) {
                continue;
            }
            if (idx_it->second == src_idx || idx_it->second == sink_idx) {
                continue;
            }
            inside_idx_set.insert(idx_it->second);
        }
    }

    // B5: add the branches the graph carries between these boundaries that no stored path takes.
    // The walked interval supplies the orientation -- which side of each boundary faces inward --
    // which the unordered cactus pair does not.
    if (best_iv.has_value()) {
        const auto& steps = graph.paths[best_iv_path].steps;
        const PathStep& near = steps[best_iv->left];
        const PathStep& far = steps[best_iv->right];
        const auto near_it = packed.node_idx_of.find(near.node_id);
        const auto far_it = packed.node_idx_of.find(far.node_id);
        if (near_it != packed.node_idx_of.end() && far_it != packed.node_idx_of.end()) {
            const std::uint64_t near_handle =
                (static_cast<std::uint64_t>(near_it->second) << 1) | (near.reverse ? 1u : 0u);
            const std::uint64_t far_handle =
                (static_cast<std::uint64_t>(far_it->second) << 1) | (far.reverse ? 1u : 0u);
            std::unordered_set<std::uint32_t> forward;
            std::unordered_set<std::uint32_t> backward;
            if (reachable_interior(graph, packed, near_handle, far_it->second, near_it->second,
                                   false, forward) &&
                reachable_interior(graph, packed, far_handle, near_it->second, far_it->second,
                                   true, backward)) {
                for (const std::uint32_t idx : forward) {
                    if (backward.count(idx) == 0) continue;   // not on a boundary-to-boundary route
                    if (idx == src_idx || idx == sink_idx) continue;
                    inside_idx_set.insert(idx);
                }
            } else {
                g_interior_traversal_truncated.fetch_add(1, std::memory_order_relaxed);
            }
        }
    }
    return inside_idx_set;
}

// B5. The interior of a snarl is a property of the GRAPH, not of the panel that happens to walk it.
//
// Collecting only the nodes some stored path visits between the boundaries understates the site
// wherever the graph carries a branch this panel does not use -- an allele nobody in the cohort
// carries is still part of the variant. Everything computed from `inside` inherits that: the
// interior-span filters, the --superbubbles acyclicity search (a cycle on an unwalked branch is
// invisible), and the fused-span validation in merging.
//
// So the interior is derived here from the graph, in the standard snarl sense: the nodes reachable
// from the near boundary's inner handle that also reach the far boundary, without leaving through a
// boundary. The traversal is over oriented HANDLES for the reason interior_has_cycle gives -- leaving
// a node forward departs its end, reversed departs its start.
//
// Two deliberate conservatisms. The result is UNIONED with the path-derived set rather than replacing
// it, so a pair that is not really a snarl cannot lose nodes the panel proves are between its
// boundaries; and the search is capped, because a leaky pair floods the whole component, in which case
// the path-derived answer stands alone.
constexpr std::size_t kMaxInteriorHandles = 1u << 20;

// Successors of an encoded handle (node index * 2 + orientation; 0 forward, 1 reverse).
// `backward` walks predecessors instead, using the bidirected identity
// pred(n, o) == { reverse(h) : h in succ(n, 1 - o) }, which needs no reverse index.
void handle_successors(
    const Graph& graph,
    const PackedGraph& packed,
    std::uint64_t handle,
    bool backward,
    std::vector<std::uint64_t>& out) {

    out.clear();
    const std::uint32_t idx = static_cast<std::uint32_t>(handle >> 1);
    const int orient = static_cast<int>(handle & 1u);
    const int depart = backward ? (1 - orient) : orient;

    const auto nit = graph.nodes.find(packed.node_ids[idx]);
    if (nit == graph.nodes.end()) return;
    const std::vector<Neighbor>& arcs = depart ? nit->second.start : nit->second.end;
    for (const Neighbor& nb : arcs) {
        const auto jt = packed.node_idx_of.find(nb.node_id);
        if (jt == packed.node_idx_of.end()) continue;
        // side 0 = entering the neighbour's start = arriving on its forward handle
        std::uint64_t next = (static_cast<std::uint64_t>(jt->second) << 1) | (nb.side == 0 ? 0u : 1u);
        if (backward) next ^= 1u;   // the predecessor is that handle reversed
        out.push_back(next);
    }
}

// Nodes reachable from `start_handle`, stopping at `stop_idx` and never passing back through
// `origin_idx`. Returns false if the cap is hit, in which case `out` is meaningless.
bool reachable_interior(
    const Graph& graph,
    const PackedGraph& packed,
    std::uint64_t start_handle,
    std::uint32_t stop_idx,
    std::uint32_t origin_idx,
    bool backward,
    std::unordered_set<std::uint32_t>& out) {

    out.clear();
    std::unordered_set<std::uint64_t> seen;
    std::vector<std::uint64_t> stack{start_handle};
    std::vector<std::uint64_t> succ;
    seen.insert(start_handle);

    while (!stack.empty()) {
        const std::uint64_t h = stack.back();
        stack.pop_back();
        handle_successors(graph, packed, h, backward, succ);
        for (const std::uint64_t nh : succ) {
            const std::uint32_t nidx = static_cast<std::uint32_t>(nh >> 1);
            if (nidx == origin_idx) continue;   // would leave through the near boundary
            if (nidx == stop_idx) continue;     // reached the far boundary; do not expand past it
            out.insert(nidx);
            if (seen.size() >= kMaxInteriorHandles) return false;
            if (seen.insert(nh).second) stack.push_back(nh);
        }
    }
    return true;
}

void collect_candidates_for_pairs(
    const Graph& graph,
    const PackedGraph& packed,
    const std::vector<std::pair<std::string, std::string>>& pairs,
    const std::vector<PathIndex>& path_indexes,
    CandidateMap& dedup,
    std::vector<SnarlDebugEntry>* debug_entries) {

    std::size_t candidate_id = 0;
    for (const auto& [source_id, sink_id] : pairs) {
        ++candidate_id;
        SnarlDebugEntry debug;
        debug.candidate_id = candidate_id;

        if (source_id == sink_id) {
            if (debug_entries != nullptr) debug_entries->push_back(std::move(debug));
            continue;
        }

        const auto src_idx_it = packed.node_idx_of.find(source_id);
        const auto sink_idx_it = packed.node_idx_of.find(sink_id);
        if (src_idx_it == packed.node_idx_of.end() || sink_idx_it == packed.node_idx_of.end()) {
            if (debug_entries != nullptr) {
                debug.source = source_id;
                debug.sink = sink_id;
                debug_entries->push_back(std::move(debug)); // endpoint not in graph
            }
            continue;
        }

        std::unordered_set<std::uint32_t> inside_idx_set = interior_indices_for_pair(
            graph, packed, path_indexes, source_id, sink_id, src_idx_it->second, sink_idx_it->second);

        std::vector<std::uint32_t> inside;
        inside.reserve(inside_idx_set.size());
        for (const std::uint32_t idx : inside_idx_set) {
            inside.push_back(idx);
        }
        std::sort(inside.begin(), inside.end());

        debug.source = source_id;
        debug.sink = sink_id;
        debug.inside_node_count = inside.size();
        if (inside.empty()) {
            if (debug_entries != nullptr) {
                debug_entries->push_back(std::move(debug)); // empty inside
            }
            continue;
        }

        BubbleCandidateIdx candidate;
        candidate.source = src_idx_it->second;
        candidate.sink = sink_idx_it->second;
        candidate.inside = std::move(inside);
        add_or_replace_candidate(dedup, std::move(candidate));
        if (debug_entries != nullptr) {
            debug_entries->push_back(std::move(debug)); // acceptance decided after filtering
        }
    }
}

// ---- Per-bubble metrics on supporting path intervals ----------------------------------



struct CandidateInterval {
    std::size_t left = 0;
    std::size_t right = 0;
    std::size_t inside_count = 0;
    bool source_to_sink = true;
};



// Does the snarl's INTERIOR contain a directed cycle?
//
// Cyclicity was previously inferred from what the stored paths happen to do -- a self-loop, a path
// revisiting an interior node, a node used in both orientations. All three are real signals, but none
// of them looks at the graph: an interior cycle that no stored path traverses was reported as acyclic
// and survived --superbubbles, which exists precisely to exclude it. A superbubble's interior must be a
// DAG whether or not anybody walked it.
//
// The search is over oriented HANDLES, not node names. A bidirected node is two vertices: leaving it
// forward departs its end, leaving it reversed departs its start, so `A.end` gives the successors of
// (A,+) and `A.start` those of (A,-). Collapsing to node names would call `A+ -> B+ -> A-` a cycle when
// it is a perfectly ordinary hairpin-free traversal of two distinct handles.
bool interior_has_cycle(const Graph& graph, const std::vector<std::string>& inside) {
    if (inside.size() < 2) return false;
    std::unordered_map<std::string, std::size_t> index;
    index.reserve(inside.size());
    for (std::size_t i = 0; i < inside.size(); ++i) index.emplace(inside[i], i);

    // vertex id = 2*node_index + orientation (0 = forward, 1 = reverse)
    const std::size_t n = inside.size() * 2;
    std::vector<std::vector<std::size_t>> adj(n);
    for (std::size_t i = 0; i < inside.size(); ++i) {
        const auto nit = graph.nodes.find(inside[i]);
        if (nit == graph.nodes.end()) continue;
        for (int orient = 0; orient < 2; ++orient) {
            const std::vector<Neighbor>& out = orient ? nit->second.start : nit->second.end;
            for (const Neighbor& nb : out) {
                const auto jt = index.find(nb.node_id);
                if (jt == index.end()) continue;              // leaves the interior: not an interior cycle
                // side 0 = entering the neighbour's start = its forward handle
                adj[i * 2 + static_cast<std::size_t>(orient)]
                    .push_back(jt->second * 2 + (nb.side == 0 ? 0u : 1u));
            }
        }
    }

    // Iterative three-colour DFS; a grey successor is a back edge and therefore a cycle.
    std::vector<char> colour(n, 0);   // 0 white, 1 grey, 2 black
    std::vector<std::pair<std::size_t, std::size_t>> stack;   // (vertex, next successor to visit)
    for (std::size_t start = 0; start < n; ++start) {
        if (colour[start] != 0) continue;
        stack.clear();
        stack.emplace_back(start, 0);
        colour[start] = 1;
        while (!stack.empty()) {
            auto& [v, next] = stack.back();
            if (next < adj[v].size()) {
                const std::size_t w = adj[v][next++];
                if (colour[w] == 1) return true;              // back edge
                if (colour[w] == 0) { colour[w] = 1; stack.emplace_back(w, 0); }
            } else {
                colour[v] = 2;
                stack.pop_back();
            }
        }
    }
    return false;
}



void assign_ids(std::vector<Bubble>& bubbles) {
    for (std::size_t i = 0; i < bubbles.size(); ++i) {
        bubbles[i].id = i + 1;
    }
}

} // namespace

std::unordered_set<std::string> nodes_in_bubble(const Bubble& bubble) {
    std::unordered_set<std::string> out;
    out.insert(bubble.source);
    out.insert(bubble.sink);
    for (const auto& n : bubble.inside) {
        out.insert(n);
    }
    return out;
}

BubbleCallReport call_bubbles_report(const Graph& graph, const BubbleCallOptions& options) {
    const PackedGraph packed = pack_graph(graph);
    CandidateMap dedup;
    dedup.reserve(packed.node_ids.size());

    std::vector<SnarlDebugEntry> snarl_debug;
    auto* debug_ptr = options.collect_snarl_debug ? &snarl_debug : nullptr;

    if (graph.paths.empty()) {
        throw std::runtime_error(
            "Bubble refinement requires paths (P/W) loaded; rerun with path-aware parsing enabled");
    }

    std::vector<PathIndex> path_indexes;
    path_indexes.reserve(graph.paths.size());
    for (const auto& p : graph.paths) {
        path_indexes.push_back(build_path_index(p));
    }

    // Snarl boundary pairs: pre-computed cactus pairs (the internal default, supplied by the
    // command after sorting) win; otherwise an external vg snarls JSONL.
    std::vector<std::pair<std::string, std::string>> snarl_pairs;
    if (options.snarl_source_supplied || !options.snarl_pairs_override.empty()) {
        // An empty override from a supplied source is a valid answer -- a linear graph has no snarl --
        // and must produce an empty table, not an error about a missing source.
        snarl_pairs = options.snarl_pairs_override;
    } else if (!options.snarls_input_path.empty()) {
        snarl_pairs = read_snarl_pairs_jsonl(options.snarls_input_path);
    } else {
        throw std::runtime_error(
            "call_bubbles_report: no snarl source — provide snarl_pairs_override (internal cactus "
            "finder) or snarls_input_path (--snarls-in)");
    }
    collect_candidates_for_pairs(graph, packed, snarl_pairs, path_indexes, dedup, debug_ptr);

    std::vector<Bubble> bubbles;
    bubbles.reserve(dedup.size());
    for (auto& [key, cand] : dedup) {
        (void)key;
        Bubble bubble;
        bubble.source = packed.node_ids[cand.source];
        bubble.sink = packed.node_ids[cand.sink];
        bubble.inside.reserve(cand.inside.size());
        for (const std::uint32_t inside_idx : cand.inside) {
            bubble.inside.push_back(packed.node_ids[inside_idx]);
        }
        bubbles.push_back(std::move(bubble));
    }

    // Per-bubble support: number of crossing paths, internal bp span (min/max), the count of
    // paths whose internal span reaches --min-variant-bp, and an inversion signal (an internal
    // node seen in both orientations across paths).
    // Which path is the reference, so allele support can be split into reference and alternate. Exact
    // name first, then a unique substring; anything ambiguous is left unresolved rather than guessed,
    // in which case every allele is reported as an alternate.
    std::size_t ref_idx = graph.paths.size();
    if (!options.reference_path.empty()) {
        for (std::size_t i = 0; i < graph.paths.size(); ++i)
            if (graph.paths[i].name == options.reference_path) { ref_idx = i; break; }
        if (ref_idx == graph.paths.size()) {
            std::size_t hits = 0, last = 0;
            for (std::size_t i = 0; i < graph.paths.size(); ++i)
                if (graph.paths[i].name.find(options.reference_path) != std::string::npos) { ++hits; last = i; }
            if (hits == 1) ref_idx = last;
        }
    }

    const auto compute_bubble_metrics = [&](std::vector<Bubble>& target, const char* progress_label) {
        // An empty label suppresses the bar (honors --quiet).
        cli::ProgressBar progress(options.quiet ? "" : progress_label, target.size());
        for (auto& bubble : target) {
            std::size_t supported_paths = 0;
            std::size_t min_inside_bp = std::numeric_limits<std::size_t>::max();
            std::size_t max_inside_bp = 0;
            std::size_t long_path_support = 0;
            bool has_inside_bp = false;
            bool inversion_signal = false;
            bool cyclic = false;

            std::unordered_set<std::string> inside_nodes(bubble.inside.begin(), bubble.inside.end());
            std::unordered_map<std::string, unsigned char> orientation_mask;
            orientation_mask.reserve(bubble.inside.size() * 2);

            // A directed cycle anywhere in the induced oriented interior, whether or not a stored path
            // walks it. This is the structural test; the path-derived signals below are additional.
            if (interior_has_cycle(graph, bubble.inside)) cyclic = true;

            // An interior node carrying a self-loop edge makes the snarl cyclic (tandem unit).
            for (const std::string& id : bubble.inside) {
                const auto nit = graph.nodes.find(id);
                if (nit == graph.nodes.end()) continue;
                bool loop = false;
                for (const Neighbor& nb : nit->second.start) if (nb.node_id == id) { loop = true; break; }
                if (!loop) for (const Neighbor& nb : nit->second.end) if (nb.node_id == id) { loop = true; break; }
                if (loop) { cyclic = true; break; }
            }

            std::unordered_map<std::string, std::size_t> allele_counts;
            std::string ref_signature;
            for (std::size_t p_idx = 0; p_idx < path_indexes.size(); ++p_idx) {
                // `bubble_steps` returns the canonical source->sink walk and, when a path crosses with no
                // interior at all, falls back to the adjacent source/sink pair. That fallback is the
                // pure-deletion allele, and it is why this must not reimplement the search: requiring an
                // interior node made a deletion count as no support and reported min_inside_bp as the
                // shortest INSERTION rather than 0.
                const auto steps_opt = bubble_steps(graph.paths[p_idx], path_indexes[p_idx], bubble);
                if (!steps_opt.has_value() || steps_opt->size() < 2) {
                    continue;
                }
                ++supported_paths;
                const std::vector<PathStep>& steps = *steps_opt;
                // The walk itself, so support can be reported per ALLELE rather than per traversal.
                {
                    std::string sig;
                    sig.reserve(steps.size() * 8);
                    for (const PathStep& st : steps) {
                        sig += st.node_id;
                        sig += st.reverse ? '-' : '+';
                        sig += ',';
                    }
                    ++allele_counts[sig];
                    if (p_idx == ref_idx) ref_signature = sig;
                }

                std::size_t inside_bp = 0;
                std::unordered_set<std::string> seen_this_path;
                for (std::size_t i = 1; i + 1 < steps.size(); ++i) {
                    const auto& step = steps[i];
                    if (inside_nodes.find(step.node_id) == inside_nodes.end()) {
                        continue;
                    }
                    const auto node_it = graph.nodes.find(step.node_id);
                    if (node_it != graph.nodes.end()) {
                        inside_bp += node_it->second.sequence.size();
                    }
                    // A haplotype revisiting an interior node within the snarl = a cycle.
                    if (!seen_this_path.insert(step.node_id).second) {
                        cyclic = true;
                    }
                    auto& mask = orientation_mask[step.node_id];
                    mask |= static_cast<unsigned char>(step.reverse ? 0x2 : 0x1);
                    if (mask == 0x3) {
                        inversion_signal = true;
                    }
                }

                min_inside_bp = std::min(min_inside_bp, inside_bp);
                max_inside_bp = std::max(max_inside_bp, inside_bp);
                has_inside_bp = true;

                if (options.min_variant_bp > 0 && inside_bp >= options.min_variant_bp) {
                    ++long_path_support;
                }
            }

            // path_support counts TRAVERSALS. On a fully-typed panel nearly every haplotype crosses
            // nearly every bubble, so it says little about any particular allele; these say what the
            // traversals actually contain.
            bubble.distinct_alleles = allele_counts.size();
            bubble.ref_allele_support = 0;
            std::size_t alt_max = 0, alt_min = 0;
            bool any_alt = false;
            for (const auto& [sig, n] : allele_counts) {
                if (!ref_signature.empty() && sig == ref_signature) { bubble.ref_allele_support = n; continue; }
                alt_max = std::max(alt_max, n);
                alt_min = any_alt ? std::min(alt_min, n) : n;
                any_alt = true;
            }
            bubble.alt_allele_support_max = alt_max;
            bubble.alt_allele_support_min = any_alt ? alt_min : 0;
            bubble.path_support = supported_paths;
            bubble.min_inside_bp = has_inside_bp ? min_inside_bp : 0;
            bubble.max_inside_bp = has_inside_bp ? max_inside_bp : 0;
            bubble.long_path_support = long_path_support;
            bubble.inversion_signal = inversion_signal;
            bubble.cyclic = cyclic || inversion_signal;
            progress.tick();
        }
    };

    compute_bubble_metrics(bubbles, "Scoring bubbles");

    std::vector<Bubble> non_snp_bubbles;
    non_snp_bubbles.reserve(bubbles.size());
    for (const auto& bubble : bubbles) {
        if (bubble.max_inside_bp > 1 || bubble.inversion_signal) {
            non_snp_bubbles.push_back(bubble);
        }
    }

    // Candidate metrics by endpoint, captured before filtering so the debug TSV can report
    // them for rejected candidates too.
    struct DebugMetrics {
        std::size_t n_paths = 0;
        std::size_t min_inside_bp = 0;
        std::size_t inside_node_count = 0;
        std::size_t long_path_support = 0;
        bool inversion_signal = false;
    };
    std::unordered_map<std::string, DebugMetrics> metrics_by_endpoint;
    metrics_by_endpoint.reserve(bubbles.size() * 2);
    for (const auto& bubble : bubbles) {
        metrics_by_endpoint[endpoint_key(bubble.source, bubble.sink)] =
            DebugMetrics{bubble.path_support, bubble.min_inside_bp, bubble.inside.size(),
                         bubble.long_path_support, bubble.inversion_signal};
    }

    // Boundaries carry reference order BEFORE anything reads them as an interval. Merging in
    // particular is defined in reference coordinates, so it has to run on oriented pairs; when this
    // ran after the merge instead, every bubble looked reference-inverted (source after sink) and
    // merging refused them all as unorderable.
    orient_bubbles_to_reference(graph, options.reference_path, bubbles);

    // Every filter, applied to whatever the current bubble set is. Merging creates NEW bubbles whose
    // metrics are not those of their parts -- two 1 bp bubbles 10 bp apart fuse into one spanning 12 bp
    // -- so a set filtered before merging is not a filtered set afterwards. Measured: with
    // `--max-variant-bp 1 --merge-nearby-bp 50` the emitted bubble had max_inside_bp 12.
    auto apply_filters = [&](std::vector<Bubble>& bs) {
        if (options.superbubbles_only) {
            bs.erase(std::remove_if(bs.begin(), bs.end(), [&](const Bubble& b) { return b.cyclic; }),
                     bs.end());
        }
        if (options.min_path_support > 0) {
            bs.erase(std::remove_if(bs.begin(), bs.end(),
                                    [&](const Bubble& b) { return b.path_support < options.min_path_support; }),
                     bs.end());
        }
        if (options.min_alt_support > 0) {
            bs.erase(std::remove_if(bs.begin(), bs.end(),
                                    [&](const Bubble& b) {
                                        return b.alt_allele_support_max < options.min_alt_support;
                                    }),
                     bs.end());
        }
        if (options.min_variant_bp > 0) {
            bs.erase(std::remove_if(bs.begin(), bs.end(),
                                    [&](const Bubble& b) {
                                        return b.long_path_support == 0 && !b.inversion_signal;
                                    }),
                     bs.end());
        }
        if (options.max_variant_bp > 0) {
            // Bounds the per-bubble walk cost so panphorte/call stay tractable on hypervariable
            // regions. Trades away SVs larger than the cap.
            bs.erase(std::remove_if(bs.begin(), bs.end(),
                                    [&](const Bubble& b) { return b.max_inside_bp > options.max_variant_bp; }),
                     bs.end());
        }
    };

    apply_filters(bubbles);

    // Always reported, quiet or not: the interior of every affected site is panel-derived rather than
    // graph-derived, so an allele nobody in this panel carries is missing from it, and everything
    // computed from `inside` -- the span filters, the --superbubbles acyclicity search, the emitted
    // interior -- silently inherits that. A user cannot tell from the output that it happened.
    if (const std::size_t truncated = g_interior_traversal_truncated.exchange(0); truncated > 0) {
        std::cerr << "[bubble] WARNING: the graph-derived interior search hit its "
                  << kMaxInteriorHandles << "-handle cap on " << truncated
                  << " endpoint pair(s); their interiors are PANEL-derived, so a branch no stored path"
                     " walks is not counted as part of those sites\n";
    }

    if (options.merge_nearby_bp > 0 && bubbles.size() > 1) {
        const std::size_t before = bubbles.size();
        bubbles = merge_nearby_bubbles(graph, packed, path_indexes, bubbles, options.merge_nearby_bp,
                                       options.reference_path);
        compute_bubble_metrics(bubbles, "Rescoring merged");
        const std::size_t merged = bubbles.size();
        apply_filters(bubbles);   // a fused bubble has to clear the same bars its parts did
        if (bubbles.size() != merged && !options.quiet) {
            std::cerr << "[bubble] " << (merged - bubbles.size())
                      << " merged bubble(s) dropped by re-applying the filters (a fusion can exceed "
                         "a bound its parts satisfied)\n";
        }
        (void)before;
    }

    // Overlapping retained sites: two bubbles claiming the same interior node describe the same
    // sequence twice. Downstream that is not a nuisance but a hard stop -- panphorte refuses to fold
    // both, because rewriting one span inside the other would corrupt the enclosing one, so a locus
    // with an overlap cannot be processed at all (ANKRD36C fails there today).
    //
    // Project policy is to keep the ENCLOSING site and drop the smaller ones it contains, so the
    // conflict resolves toward more sequence rather than less. Reported, because dropping a site
    // loses finer resolution and that should be visible rather than silent.
    if (bubbles.size() > 1) {
        std::sort(bubbles.begin(), bubbles.end(), [](const Bubble& a, const Bubble& b) {
            if (a.inside.size() != b.inside.size()) return a.inside.size() < b.inside.size();  // EXPERIMENT: prefer smaller
            return bubble_endpoint_less(a, b);   // deterministic among equals
        });
        std::unordered_set<std::string> claimed;
        std::vector<Bubble> kept;
        std::size_t dropped = 0;
        for (Bubble& b : bubbles) {
            const std::string* clash = nullptr;
            for (const std::string& n : b.inside) {
                if (claimed.count(n)) { clash = &n; break; }
            }
            if (clash != nullptr) {
                ++dropped;
                if (!options.quiet) {
                    std::cerr << "[bubble] dropped site " << b.source << ".." << b.sink << " ("
                              << b.inside.size() << " interior node(s)): its interior node " << *clash
                              << " is already claimed by a smaller retained site\n";
                }
                continue;
            }
            for (const std::string& n : b.inside) claimed.insert(n);
            kept.push_back(std::move(b));
        }
        if (dropped > 0 && !options.quiet) {
            std::cerr << "[bubble] " << dropped << " overlapping site(s) dropped so the emitted set is "
                         "pairwise disjoint; the smaller sites are retained, which keeps the finer "
                         "resolution an enclosing site would have combined\n";
        }
        bubbles = std::move(kept);
    }

    std::sort(bubbles.begin(), bubbles.end(), bubble_endpoint_less);
    assign_ids(bubbles);

    std::unordered_set<std::string> final_endpoint_keys;
    final_endpoint_keys.reserve(bubbles.size() * 2);
    for (const auto& bubble : bubbles) {
        final_endpoint_keys.insert(endpoint_key(bubble.source, bubble.sink));
    }

    if (options.collect_snarl_debug) {
        std::size_t next_debug_id = 0;
        std::unordered_set<std::string> debug_endpoint_keys;
        debug_endpoint_keys.reserve(snarl_debug.size() * 2 + bubbles.size() * 2);
        for (const auto& debug : snarl_debug) {
            next_debug_id = std::max(next_debug_id, debug.candidate_id);
            if (!debug.source.empty() && !debug.sink.empty()) {
                debug_endpoint_keys.insert(endpoint_key(debug.source, debug.sink));
            }
        }

        for (auto& debug : snarl_debug) {
            if (debug.source.empty() || debug.sink.empty()) {
                continue;
            }
            const std::string key = endpoint_key(debug.source, debug.sink);
            const auto metrics_it = metrics_by_endpoint.find(key);
            if (metrics_it != metrics_by_endpoint.end()) {
                debug.n_paths = metrics_it->second.n_paths;
                debug.min_inside_bp = metrics_it->second.min_inside_bp;
                if (debug.inside_node_count == 0) {
                    debug.inside_node_count = metrics_it->second.inside_node_count;
                }
                debug.long_path_support = metrics_it->second.long_path_support;
                debug.inversion_signal = metrics_it->second.inversion_signal;
            }
            debug.accepted = final_endpoint_keys.find(key) != final_endpoint_keys.end();
        }

        // Final bubbles produced by merging have no original candidate row; add one each.
        for (const auto& bubble : bubbles) {
            const std::string key = endpoint_key(bubble.source, bubble.sink);
            if (debug_endpoint_keys.find(key) != debug_endpoint_keys.end()) {
                continue;
            }
            SnarlDebugEntry synthetic;
            synthetic.candidate_id = ++next_debug_id;
            synthetic.source = bubble.source;
            synthetic.sink = bubble.sink;
            synthetic.inside_node_count = bubble.inside.size();
            synthetic.n_paths = bubble.path_support;
            synthetic.min_inside_bp = bubble.min_inside_bp;
            synthetic.long_path_support = bubble.long_path_support;
            synthetic.inversion_signal = bubble.inversion_signal;
            synthetic.accepted = true;
            snarl_debug.push_back(std::move(synthetic));
        }
    }

    BubbleCallReport report;
    report.bubbles = std::move(bubbles);
    report.non_snp_bubbles = std::move(non_snp_bubbles);
    if (options.collect_snarl_debug) {
        report.snarl_debug = std::move(snarl_debug);
    }
    return report;
}

std::vector<Bubble> call_bubbles(const Graph& graph, const BubbleCallOptions& options) {
    return call_bubbles_report(graph, options).bubbles;
}

} // namespace panvar
