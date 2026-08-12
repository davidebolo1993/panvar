#include "panvar/bubbles.hpp"

#include "panvar/bubble_path.hpp"

#include "panvar/cli_utils.hpp"

#include <algorithm>
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

bool bubble_endpoint_less(const Bubble& a, const Bubble& b) {
    if (a.source != b.source) {
        return node_id_less(a.source, b.source);
    }
    if (a.sink != b.sink) {
        return node_id_less(a.sink, b.sink);
    }
    return a.inside.size() < b.inside.size();
}

// ---- Nearby-bubble merge (graph-distance based) --------------------------------------

struct NodeAdjacency {
    std::unordered_map<std::string, std::vector<std::string>> neighbors;
    std::unordered_map<std::string, std::size_t> node_len_bp;
};

NodeAdjacency build_node_adjacency(const Graph& graph) {
    NodeAdjacency out;
    out.neighbors.reserve(graph.nodes.size() * 2);
    out.node_len_bp.reserve(graph.nodes.size() * 2);

    for (const auto& [node_id, node] : graph.nodes) {
        out.node_len_bp[node_id] = std::max<std::size_t>(1, node.sequence.size());
    }

    auto add_neighbor = [&](const std::string& from, const std::string& to) {
        out.neighbors[from].push_back(to);
    };

    for (const auto& [node_id, node] : graph.nodes) {
        for (const auto& n : node.start) {
            if (graph.nodes.find(n.node_id) == graph.nodes.end()) {
                continue;
            }
            add_neighbor(node_id, n.node_id);
            add_neighbor(n.node_id, node_id);
        }
        for (const auto& n : node.end) {
            if (graph.nodes.find(n.node_id) == graph.nodes.end()) {
                continue;
            }
            add_neighbor(node_id, n.node_id);
            add_neighbor(n.node_id, node_id);
        }
    }

    for (auto& [node_id, vec] : out.neighbors) {
        (void)node_id;
        std::sort(vec.begin(), vec.end());
        vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
    }

    return out;
}

std::optional<std::vector<std::string>> shortest_node_path_within_bp(
    const NodeAdjacency& adj,
    const std::string& from_node,
    const std::string& to_node,
    std::size_t max_bp) {

    const auto it_from = adj.node_len_bp.find(from_node);
    const auto it_to = adj.node_len_bp.find(to_node);
    if (it_from == adj.node_len_bp.end() || it_to == adj.node_len_bp.end()) {
        return std::nullopt;
    }

    const std::size_t start_cost = it_from->second;
    if (start_cost > max_bp) {
        return std::nullopt;
    }

    using HeapItem = std::pair<std::size_t, std::string>;
    std::priority_queue<HeapItem, std::vector<HeapItem>, std::greater<HeapItem>> heap;
    std::unordered_map<std::string, std::size_t> dist;
    std::unordered_map<std::string, std::string> prev;
    dist.reserve(adj.node_len_bp.size() / 2 + 16);
    prev.reserve(adj.node_len_bp.size() / 2 + 16);

    dist[from_node] = start_cost;
    heap.push({start_cost, from_node});

    while (!heap.empty()) {
        auto [d, node] = heap.top();
        heap.pop();
        const auto dist_it = dist.find(node);
        if (dist_it == dist.end() || d != dist_it->second) {
            continue;
        }
        if (d > max_bp) {
            continue;
        }
        if (node == to_node) {
            break;
        }

        const auto neigh_it = adj.neighbors.find(node);
        if (neigh_it == adj.neighbors.end()) {
            continue;
        }
        for (const auto& nxt : neigh_it->second) {
            const auto len_it = adj.node_len_bp.find(nxt);
            if (len_it == adj.node_len_bp.end()) {
                continue;
            }
            const std::size_t nd = d + len_it->second;
            if (nd > max_bp) {
                continue;
            }
            const auto old_it = dist.find(nxt);
            if (old_it == dist.end() || nd < old_it->second) {
                dist[nxt] = nd;
                prev[nxt] = node;
                heap.push({nd, nxt});
            }
        }
    }

    const auto to_dist_it = dist.find(to_node);
    if (to_dist_it == dist.end() || to_dist_it->second > max_bp) {
        return std::nullopt;
    }

    std::vector<std::string> path_rev;
    path_rev.push_back(to_node);
    while (!path_rev.empty() && path_rev.back() != from_node) {
        const auto p_it = prev.find(path_rev.back());
        if (p_it == prev.end()) {
            return std::nullopt;
        }
        path_rev.push_back(p_it->second);
    }
    if (path_rev.empty() || path_rev.back() != from_node) {
        return std::nullopt;
    }
    std::reverse(path_rev.begin(), path_rev.end());
    return path_rev;
}

std::vector<Bubble> merge_nearby_bubbles(
    const Graph& graph,
    const std::vector<Bubble>& input,
    std::size_t max_bp) {

    if (max_bp == 0 || input.size() < 2) {
        return input;
    }

    std::vector<Bubble> bubbles = input;
    std::sort(bubbles.begin(), bubbles.end(), bubble_endpoint_less);

    const NodeAdjacency adj = build_node_adjacency(graph);

    std::vector<Bubble> merged;
    merged.reserve(bubbles.size());
    Bubble current = bubbles.front();
    for (std::size_t i = 1; i < bubbles.size(); ++i) {
        const Bubble& next = bubbles[i];
        const auto path_opt = shortest_node_path_within_bp(adj, current.sink, next.source, max_bp);
        if (!path_opt.has_value()) {
            merged.push_back(std::move(current));
            current = next;
            continue;
        }

        std::unordered_set<std::string> inside_set;
        inside_set.reserve((current.inside.size() + next.inside.size() + path_opt->size()) * 2);
        for (const auto& n : current.inside) {
            inside_set.insert(n);
        }
        for (const auto& n : next.inside) {
            inside_set.insert(n);
        }
        // Old boundary nodes become internal when two nearby bubbles are fused.
        inside_set.insert(current.sink);
        inside_set.insert(next.source);
        // Include any connector nodes on the shortest path between boundaries.
        for (std::size_t pi = 1; pi + 1 < path_opt->size(); ++pi) {
            inside_set.insert((*path_opt)[pi]);
        }

        Bubble fused = current;
        fused.sink = next.sink;
        fused.path_support = 0;
        fused.min_inside_bp = 0;
        fused.max_inside_bp = 0;
        fused.long_path_support = 0;
        fused.inversion_signal = false;
        fused.inside.clear();
        fused.inside.reserve(inside_set.size());
        for (const auto& node_id : inside_set) {
            if (node_id == fused.source || node_id == fused.sink) {
                continue;
            }
            fused.inside.push_back(node_id);
        }
        std::sort(fused.inside.begin(), fused.inside.end());
        current = std::move(fused);
    }
    merged.push_back(std::move(current));
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

    for (const std::size_t start_pos : start_positions) {
        const auto end_it = std::upper_bound(end_positions.begin(), end_positions.end(), start_pos);
        if (end_it == end_positions.end()) {
            continue;
        }
        const std::size_t end_pos = *end_it;
        if (end_pos <= start_pos + 1) {
            continue;
        }

        EndpointOnlyInterval candidate;
        candidate.left = start_pos;
        candidate.right = end_pos;
        candidate.inside_steps = end_pos - start_pos - 1;
        candidate.source_to_sink = source_to_sink;
        if (!best.has_value() || better_endpoint_interval(candidate, *best)) {
            best = candidate;
        }
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

        std::unordered_set<std::uint32_t> inside_idx_set;
        for (std::size_t p_idx = 0; p_idx < path_indexes.size(); ++p_idx) {
            const auto interval = find_best_endpoint_interval(path_indexes[p_idx], source_id, sink_id);
            if (!interval.has_value()) {
                continue;
            }
            const auto& steps = graph.paths[p_idx].steps;
            for (std::size_t i = interval->left + 1; i < interval->right; ++i) {
                const auto idx_it = packed.node_idx_of.find(steps[i].node_id);
                if (idx_it == packed.node_idx_of.end()) {
                    continue;
                }
                if (idx_it->second == src_idx_it->second || idx_it->second == sink_idx_it->second) {
                    continue;
                }
                inside_idx_set.insert(idx_it->second);
            }
        }

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

std::vector<std::size_t> collect_inside_positions(const PathIndex& path, const Bubble& bubble) {
    std::vector<std::size_t> inside_positions;
    inside_positions.reserve(bubble.inside.size());
    for (const auto& node : bubble.inside) {
        const auto it = path.positions.find(node);
        if (it == path.positions.end()) {
            continue;
        }
        inside_positions.insert(inside_positions.end(), it->second.begin(), it->second.end());
    }
    if (inside_positions.empty()) {
        return inside_positions;
    }
    std::sort(inside_positions.begin(), inside_positions.end());
    inside_positions.erase(std::unique(inside_positions.begin(), inside_positions.end()), inside_positions.end());
    return inside_positions;
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

struct CandidateInterval {
    std::size_t left = 0;
    std::size_t right = 0;
    std::size_t inside_count = 0;
    bool source_to_sink = true;
};

bool better_candidate(const CandidateInterval& a, const CandidateInterval& b) {
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

void evaluate_direction_candidates(
    const std::vector<std::size_t>& start_positions,
    const std::vector<std::size_t>& end_positions,
    bool source_to_sink,
    const std::vector<std::size_t>& inside_positions,
    std::optional<CandidateInterval>& best) {

    for (const std::size_t start_pos : start_positions) {
        const auto end_it = std::upper_bound(end_positions.begin(), end_positions.end(), start_pos);
        if (end_it == end_positions.end()) {
            continue;
        }
        const std::size_t end_pos = *end_it;
        const std::size_t inside_count = count_inside_between(inside_positions, start_pos, end_pos);
        if (inside_count == 0) {
            continue;
        }

        CandidateInterval candidate;
        candidate.left = start_pos;
        candidate.right = end_pos;
        candidate.inside_count = inside_count;
        candidate.source_to_sink = source_to_sink;

        if (!best.has_value() || better_candidate(candidate, *best)) {
            best = candidate;
        }
    }
}

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

std::optional<CandidateInterval> find_best_interval(const PathIndex& index, const Bubble& bubble) {
    const auto src_it = index.positions.find(bubble.source);
    const auto sink_it = index.positions.find(bubble.sink);
    if (src_it == index.positions.end() || sink_it == index.positions.end()) {
        return std::nullopt;
    }

    const auto inside_positions = collect_inside_positions(index, bubble);
    if (inside_positions.empty()) {
        return std::nullopt;
    }

    std::optional<CandidateInterval> best;
    evaluate_direction_candidates(src_it->second, sink_it->second, true, inside_positions, best);
    evaluate_direction_candidates(sink_it->second, src_it->second, false, inside_positions, best);
    return best;
}

std::vector<PathStep> canonical_steps_for_bubble(
    const PathRecord& path,
    const Bubble& bubble,
    const CandidateInterval& interval) {

    if (interval.left >= path.steps.size() || interval.right >= path.steps.size() || interval.left > interval.right) {
        return {};
    }

    std::vector<PathStep> out;
    out.reserve(interval.right - interval.left + 1);

    if (interval.source_to_sink) {
        for (std::size_t i = interval.left; i <= interval.right; ++i) {
            out.push_back(path.steps[i]);
        }
    } else {
        for (std::size_t i = interval.right + 1; i > interval.left; --i) {
            const PathStep& s = path.steps[i - 1];
            out.push_back(PathStep{s.node_id, !s.reverse});
        }
        if (!out.empty() && (out.front().node_id != bubble.source || out.back().node_id != bubble.sink)) {
            std::reverse(out.begin(), out.end());
            for (auto& step : out) {
                step.reverse = !step.reverse;
            }
        }
    }

    return out;
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

    if (options.merge_nearby_bp > 0 && bubbles.size() > 1) {
        const std::size_t before = bubbles.size();
        bubbles = merge_nearby_bubbles(graph, bubbles, options.merge_nearby_bp);
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
