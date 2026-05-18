#include "panvar/bubbles.hpp"

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <functional>
#include <limits>
#include <optional>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace panvar {
namespace {

struct PackedNeighbor {
    std::uint32_t node_idx = 0;
    int side = 0;
};

struct PackedGraph {
    std::vector<std::string> node_ids;
    std::unordered_map<std::string, std::uint32_t> node_idx_of;
    std::vector<std::vector<PackedNeighbor>> start;
    std::vector<std::vector<PackedNeighbor>> end;
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

struct PathIndex {
    std::unordered_map<std::string, std::vector<std::size_t>> positions;
};

PathIndex build_path_index(const PathRecord& path);

struct CandidateInterval {
    std::size_t left = 0;
    std::size_t right = 0;
    std::size_t inside_count = 0;
    bool source_to_sink = true;
};

struct EndpointOnlyInterval {
    std::size_t left = 0;
    std::size_t right = 0;
    std::size_t inside_steps = 0;
    bool source_to_sink = true;
};

struct SearchWorkspace {
    std::vector<std::uint64_t> visited_epoch;
    std::vector<std::uint64_t> seen_epoch;
    std::vector<std::uint64_t> open_epoch;
    std::vector<std::uint64_t> inside_epoch;
    std::uint64_t search_epoch = 0;
    std::uint64_t inside_epoch_id = 0;

    explicit SearchWorkspace(std::size_t node_count)
        : visited_epoch(node_count, 0),
          seen_epoch(node_count * 2, 0),
          open_epoch(node_count * 2, 0),
          inside_epoch(node_count, 0) {}

    std::uint64_t next_search_epoch() {
        ++search_epoch;
        return search_epoch;
    }

    std::uint64_t next_inside_epoch() {
        ++inside_epoch_id;
        return inside_epoch_id;
    }
};

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

    packed.start.resize(packed.node_ids.size());
    packed.end.resize(packed.node_ids.size());

    for (std::size_t i = 0; i < packed.node_ids.size(); ++i) {
        const std::string& node_id = packed.node_ids[i];
        const auto it = graph.nodes.find(node_id);
        if (it == graph.nodes.end()) {
            continue;
        }
        const Node& node = it->second;

        auto& start_vec = packed.start[i];
        start_vec.reserve(node.start.size());
        for (const auto& n : node.start) {
            const auto idx_it = packed.node_idx_of.find(n.node_id);
            if (idx_it == packed.node_idx_of.end()) {
                continue;
            }
            start_vec.push_back(PackedNeighbor{idx_it->second, n.side});
        }

        auto& end_vec = packed.end[i];
        end_vec.reserve(node.end.size());
        for (const auto& n : node.end) {
            const auto idx_it = packed.node_idx_of.find(n.node_id);
            if (idx_it == packed.node_idx_of.end()) {
                continue;
            }
            end_vec.push_back(PackedNeighbor{idx_it->second, n.side});
        }
    }

    return packed;
}

const std::vector<PackedNeighbor>& side_children(
    const PackedGraph& graph,
    std::uint32_t node_idx,
    int direction) {

    if (node_idx >= graph.node_ids.size()) {
        throw std::runtime_error("Invalid node index in bubble traversal");
    }
    return direction == 0 ? graph.start[node_idx] : graph.end[node_idx];
}

EndpointKey normalize_endpoint_key(std::uint32_t source, std::uint32_t sink) {
    if (source < sink) {
        return EndpointKey{source, sink};
    }
    return EndpointKey{sink, source};
}

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

std::string join_strings(const std::vector<std::string>& values, char sep) {
    std::ostringstream oss;
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (i > 0) {
            oss << sep;
        }
        oss << values[i];
    }
    return oss.str();
}

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
        auto& vec = out.neighbors[from];
        vec.push_back(to);
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
        fused.type = BubbleType::Super;
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

std::optional<std::uint32_t> find_single_open_handle(
    const std::vector<std::uint64_t>& open_epoch,
    std::uint64_t epoch,
    const std::vector<std::uint32_t>& open_stack) {

    for (auto it = open_stack.rbegin(); it != open_stack.rend(); ++it) {
        if (open_epoch[*it] == epoch) {
            return *it;
        }
    }
    return std::nullopt;
}

std::optional<BubbleCandidateIdx> find_superbubble_from_source(
    const PackedGraph& graph,
    std::uint32_t source,
    int direction,
    std::size_t max_nodes_per_search,
    SearchWorkspace& ws) {

    const std::uint64_t epoch = ws.next_search_epoch();
    std::vector<std::uint32_t> open_stack;
    std::vector<std::uint32_t> nodes_inside;
    open_stack.reserve(64);
    nodes_inside.reserve(64);

    std::size_t open_count = 0;
    std::size_t seen_count = 0;
    std::size_t visited_count = 0;

    auto mark_seen = [&](std::uint32_t handle) {
        if (ws.seen_epoch[handle] != epoch) {
            ws.seen_epoch[handle] = epoch;
            ++seen_count;
        }
    };
    auto unmark_seen = [&](std::uint32_t handle) {
        if (ws.seen_epoch[handle] == epoch) {
            ws.seen_epoch[handle] = 0;
            --seen_count;
        }
    };
    auto push_open = [&](std::uint32_t handle) {
        if (ws.open_epoch[handle] != epoch) {
            ws.open_epoch[handle] = epoch;
            open_stack.push_back(handle);
            ++open_count;
        }
    };

    const std::uint32_t start_handle = (source * 2U) + static_cast<std::uint32_t>(direction);
    mark_seen(start_handle);
    push_open(start_handle);

    while (open_count > 0) {
        std::optional<std::uint32_t> current_handle;
        while (!open_stack.empty()) {
            const std::uint32_t h = open_stack.back();
            open_stack.pop_back();
            if (ws.open_epoch[h] != epoch) {
                continue;
            }
            ws.open_epoch[h] = 0;
            --open_count;
            current_handle = h;
            break;
        }
        if (!current_handle.has_value()) {
            break;
        }

        const std::uint32_t current_node = *current_handle / 2U;
        const int current_direction = static_cast<int>(*current_handle & 1U);

        if (ws.visited_epoch[current_node] == epoch) {
            continue;
        }
        ws.visited_epoch[current_node] = epoch;
        ++visited_count;
        nodes_inside.push_back(current_node);
        unmark_seen(*current_handle);

        const auto& children = side_children(graph, current_node, current_direction);
        if (children.empty()) {
            break;
        }

        bool loop_to_source = false;
        for (const auto& child : children) {
            if (child.node_idx == source) {
                loop_to_source = true;
                break;
            }

            const int child_direction = (child.side == 0) ? 1 : 0;
            const std::uint32_t child_handle =
                (child.node_idx * 2U) + static_cast<std::uint32_t>(child_direction);
            mark_seen(child_handle);

            const auto& parent_side = side_children(graph, child.node_idx, child.side);
            bool all_parents_visited = true;
            for (const auto& p : parent_side) {
                if (ws.visited_epoch[p.node_idx] != epoch) {
                    all_parents_visited = false;
                    break;
                }
            }
            if (all_parents_visited) {
                push_open(child_handle);
            }
        }

        if (loop_to_source) {
            break;
        }

        if (open_count == 1 && seen_count == 1) {
            const auto sink_handle = find_single_open_handle(ws.open_epoch, epoch, open_stack);
            if (!sink_handle.has_value()) {
                break;
            }
            const std::uint32_t sink_node = *sink_handle / 2U;
            nodes_inside.push_back(sink_node);

            if (nodes_inside.size() == 2) {
                break;
            }

            const std::uint64_t inside_epoch = ws.next_inside_epoch();
            std::vector<std::uint32_t> inside;
            inside.reserve(nodes_inside.size());
            for (const std::uint32_t node_idx : nodes_inside) {
                if (node_idx == source || node_idx == sink_node) {
                    continue;
                }
                if (ws.inside_epoch[node_idx] == inside_epoch) {
                    continue;
                }
                ws.inside_epoch[node_idx] = inside_epoch;
                inside.push_back(node_idx);
            }

            if (inside.empty()) {
                break;
            }

            return BubbleCandidateIdx{source, sink_node, std::move(inside)};
        }

        if (visited_count > max_nodes_per_search) {
            break;
        }
    }

    return std::nullopt;
}

[[maybe_unused]] void collect_superbubble_candidates(
    const PackedGraph& packed,
    std::size_t max_nodes_per_search,
    CandidateMap& dedup) {

    SearchWorkspace ws(packed.node_ids.size());
    for (std::uint32_t node_idx = 0; node_idx < packed.node_ids.size(); ++node_idx) {
        for (int direction = 0; direction <= 1; ++direction) {
            auto candidate =
                find_superbubble_from_source(packed, node_idx, direction, max_nodes_per_search, ws);
            if (!candidate.has_value()) {
                continue;
            }
            add_or_replace_candidate(dedup, std::move(*candidate));
        }
    }
}

struct HandleGraph {
    std::vector<std::vector<std::uint32_t>> out;
    std::vector<std::vector<std::uint32_t>> in;
};

HandleGraph build_handle_graph(const PackedGraph& packed) {
    HandleGraph hg;
    const std::size_t handle_count = packed.node_ids.size() * 2;
    hg.out.resize(handle_count);
    hg.in.resize(handle_count);

    for (std::uint32_t node_idx = 0; node_idx < packed.node_ids.size(); ++node_idx) {
        for (int side = 0; side <= 1; ++side) {
            const std::uint32_t handle = (node_idx * 2U) + static_cast<std::uint32_t>(side);
            const auto& children = side_children(packed, node_idx, side);
            auto& out_vec = hg.out[handle];
            out_vec.reserve(children.size());
            for (const auto& child : children) {
                const int child_direction = (child.side == 0) ? 1 : 0;
                const std::uint32_t child_handle =
                    (child.node_idx * 2U) + static_cast<std::uint32_t>(child_direction);
                out_vec.push_back(child_handle);
                hg.in[child_handle].push_back(handle);
            }
        }
    }

    return hg;
}

struct SccResult {
    std::vector<int> component_of_handle;
    int component_count = 0;
};

SccResult compute_handle_sccs(const HandleGraph& hg) {
    const std::size_t n = hg.out.size();
    std::vector<int> index(n, -1);
    std::vector<int> lowlink(n, -1);
    std::vector<bool> on_stack(n, false);
    std::vector<std::uint32_t> stack;
    stack.reserve(n);
    std::vector<int> component_of_handle(n, -1);

    int next_index = 0;
    int component_count = 0;

    std::function<void(std::uint32_t)> strongconnect = [&](std::uint32_t v) {
        index[v] = next_index;
        lowlink[v] = next_index;
        ++next_index;
        stack.push_back(v);
        on_stack[v] = true;

        for (const std::uint32_t w : hg.out[v]) {
            if (index[w] < 0) {
                strongconnect(w);
                lowlink[v] = std::min(lowlink[v], lowlink[w]);
            } else if (on_stack[w]) {
                lowlink[v] = std::min(lowlink[v], index[w]);
            }
        }

        if (lowlink[v] != index[v]) {
            return;
        }

        while (!stack.empty()) {
            const std::uint32_t w = stack.back();
            stack.pop_back();
            on_stack[w] = false;
            component_of_handle[w] = component_count;
            if (w == v) {
                break;
            }
        }
        ++component_count;
    };

    for (std::uint32_t v = 0; v < n; ++v) {
        if (index[v] < 0) {
            strongconnect(v);
        }
    }

    return SccResult{std::move(component_of_handle), component_count};
}

bool has_degree_one_each_side(const Graph& graph, const std::string& node_id) {
    const auto it = graph.nodes.find(node_id);
    if (it == graph.nodes.end()) {
        return false;
    }
    return it->second.start.size() == 1 && it->second.end.size() == 1;
}

bool is_simple_bubble(const Graph& graph, const Bubble& bubble) {
    if (bubble.inside.size() != 2) {
        return false;
    }

    const std::string& a = bubble.inside[0];
    const std::string& b = bubble.inside[1];

    if (!has_degree_one_each_side(graph, a) || !has_degree_one_each_side(graph, b)) {
        return false;
    }

    const auto neigh_a = graph.neighbors_of(a);
    const auto neigh_b = graph.neighbors_of(b);
    if (neigh_a != neigh_b) {
        return false;
    }

    if (graph.has_edge_between(bubble.source, bubble.sink) || graph.has_edge_between(bubble.sink, bubble.source)) {
        return false;
    }

    return true;
}

bool is_insertion_bubble(const Graph& graph, const Bubble& bubble) {
    if (bubble.inside.size() != 1) {
        return false;
    }

    const std::string& mid = bubble.inside.front();
    if (!has_degree_one_each_side(graph, mid)) {
        return false;
    }

    auto neighbors = graph.neighbors_of(mid);
    std::vector<std::string> expected{bubble.source, bubble.sink};
    std::sort(neighbors.begin(), neighbors.end());
    std::sort(expected.begin(), expected.end());
    return neighbors == expected;
}

PathIndex build_path_index(const PathRecord& path) {
    PathIndex out;
    out.positions.reserve(path.steps.size());

    for (std::size_t i = 0; i < path.steps.size(); ++i) {
        const auto& step = path.steps[i];
        out.positions[step.node_id].push_back(i);
    }

    return out;
}

std::string lowercase_ascii(std::string value) {
    for (char& c : value) {
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }
    return value;
}

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

void collect_snarl_jsonl_candidates(
    const Graph& graph,
    const PackedGraph& packed,
    const BubbleCallOptions& options,
    CandidateMap& dedup,
    std::vector<SnarlDebugEntry>* debug_entries) {

    if (options.snarls_input_path.empty()) {
        throw std::runtime_error("Bubble refinement requires --snarls-in");
    }
    if (graph.paths.empty()) {
        throw std::runtime_error(
            "Bubble refinement requires paths (P/W) loaded; rerun with path-aware parsing enabled");
    }

    std::vector<PathIndex> path_indexes;
    path_indexes.reserve(graph.paths.size());
    for (const auto& path : graph.paths) {
        path_indexes.push_back(build_path_index(path));
    }

    std::ifstream snarls_in(options.snarls_input_path);
    if (!snarls_in) {
        throw std::runtime_error("Failed to open snarl JSONL input: " + options.snarls_input_path);
    }

    std::string line;
    std::size_t candidate_id = 0;
    while (std::getline(snarls_in, line)) {
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        if (line.find("\"start\"") == std::string::npos || line.find("\"end\"") == std::string::npos) {
            continue;
        }
        // Keep top-level snarls only.
        if (line.find("\"parent\"") != std::string::npos) {
            continue;
        }

        ++candidate_id;
        SnarlDebugEntry debug;
        debug.candidate_id = candidate_id;
        debug.mode = "snarl-jsonl";

        const auto source_opt = extract_node_id_for_key(line, "start", true);
        const auto sink_opt = extract_node_id_for_key(line, "end", false);
        if (!source_opt.has_value() || !sink_opt.has_value() || *source_opt == *sink_opt) {
            if (debug_entries != nullptr) {
                debug.accepted = false;
                debug.reason = "reject:parse-start-end";
                debug_entries->push_back(std::move(debug));
            }
            continue;
        }

        const auto src_idx_it = packed.node_idx_of.find(*source_opt);
        const auto sink_idx_it = packed.node_idx_of.find(*sink_opt);
        if (src_idx_it == packed.node_idx_of.end() || sink_idx_it == packed.node_idx_of.end()) {
            if (debug_entries != nullptr) {
                debug.accepted = false;
                debug.reason = "reject:endpoint-not-in-graph";
                debug.source = *source_opt;
                debug.sink = *sink_opt;
                debug_entries->push_back(std::move(debug));
            }
            continue;
        }

        std::unordered_set<std::uint32_t> inside_idx_set;
        for (std::size_t p_idx = 0; p_idx < path_indexes.size(); ++p_idx) {
            const auto interval = find_best_endpoint_interval(path_indexes[p_idx], *source_opt, *sink_opt);
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

        debug.source = *source_opt;
        debug.sink = *sink_opt;
        debug.inside_node_count = inside.size();
        if (inside.empty()) {
            if (debug_entries != nullptr) {
                debug.accepted = false;
                debug.reason = "reject:empty-inside";
                debug_entries->push_back(std::move(debug));
            }
            continue;
        }

        BubbleCandidateIdx candidate;
        candidate.source = src_idx_it->second;
        candidate.sink = sink_idx_it->second;
        candidate.inside = std::move(inside);
        const bool kept = add_or_replace_candidate(dedup, std::move(candidate));
        if (debug_entries != nullptr) {
            debug.accepted = kept;
            debug.reason = kept ? "accept:snarl-jsonl" : "reject:duplicate-endpoints-smaller-inside";
            debug_entries->push_back(std::move(debug));
        }
    }

}

std::optional<std::size_t> pick_reference_path_index_for_snarl_macro(
    const Graph& graph,
    [[maybe_unused]] const BubbleCallOptions& options) {

    if (graph.paths.empty()) {
        return std::nullopt;
    }

    auto score_name = [](const std::string& path_name, std::size_t step_count) -> long long {
        const std::string lowered = lowercase_ascii(path_name);
        long long score = static_cast<long long>(step_count);
        if (lowered.find("grch38_1") != std::string::npos || lowered.find("grch38#1") != std::string::npos) {
            score += 1000000000LL;
        } else if (lowered.find("grch38") != std::string::npos) {
            score += 900000000LL;
        }
        if (lowered.find("chm13_1") != std::string::npos || lowered.find("chm13#1") != std::string::npos) {
            score += 800000000LL;
        } else if (lowered.find("chm13") != std::string::npos) {
            score += 700000000LL;
        }
        if (lowered.find("reference") != std::string::npos) {
            score += 500000000LL;
        } else if (lowered.find("ref") != std::string::npos) {
            score += 100000000LL;
        }
        return score;
    };

    std::size_t best_idx = 0;
    long long best_score = std::numeric_limits<long long>::min();
    std::string best_name;
    for (std::size_t i = 0; i < graph.paths.size(); ++i) {
        const auto& path = graph.paths[i];
        const long long score = score_name(path.name, path.steps.size());
        if (score > best_score || (score == best_score && path.name < best_name)) {
            best_idx = i;
            best_score = score;
            best_name = path.name;
        }
    }
    return std::optional<std::size_t>(best_idx);
}

std::vector<std::vector<std::size_t>> build_reference_positions_by_node(
    const Graph& graph,
    const PackedGraph& packed,
    const std::optional<std::size_t>& ref_path_idx) {

    std::vector<std::vector<std::size_t>> ref_positions_by_node(packed.node_ids.size());
    if (!ref_path_idx.has_value()) {
        return ref_positions_by_node;
    }
    if (*ref_path_idx >= graph.paths.size()) {
        return ref_positions_by_node;
    }

    const auto& reference_path = graph.paths[*ref_path_idx];
    for (std::size_t ref_pos = 0; ref_pos < reference_path.steps.size(); ++ref_pos) {
        const std::string& node_id = reference_path.steps[ref_pos].node_id;
        const auto idx_it = packed.node_idx_of.find(node_id);
        if (idx_it == packed.node_idx_of.end()) {
            continue;
        }
        ref_positions_by_node[idx_it->second].push_back(ref_pos);
    }
    return ref_positions_by_node;
}

std::optional<std::pair<std::uint32_t, std::uint32_t>> choose_macro_boundary_pair(
    const std::vector<std::uint32_t>& boundaries,
    const std::vector<std::vector<std::size_t>>& ref_positions_by_node) {

    if (boundaries.size() < 2) {
        return std::nullopt;
    }

    struct AnchoredBoundary {
        std::size_t pos = 0;
        std::uint32_t node = 0;
    };
    std::vector<AnchoredBoundary> anchored;
    anchored.reserve(boundaries.size());
    for (const std::uint32_t node_idx : boundaries) {
        if (node_idx >= ref_positions_by_node.size()) {
            continue;
        }
        const auto& positions = ref_positions_by_node[node_idx];
        if (positions.empty()) {
            continue;
        }
        anchored.push_back(AnchoredBoundary{positions.front(), node_idx});
    }

    if (anchored.size() >= 2) {
        const auto min_it = std::min_element(
            anchored.begin(), anchored.end(), [](const AnchoredBoundary& a, const AnchoredBoundary& b) {
                return a.pos < b.pos;
            });
        const auto max_it = std::max_element(
            anchored.begin(), anchored.end(), [](const AnchoredBoundary& a, const AnchoredBoundary& b) {
                return a.pos < b.pos;
            });
        if (min_it != anchored.end() && max_it != anchored.end() && min_it->node != max_it->node) {
            return std::make_pair(min_it->node, max_it->node);
        }
    }

    std::vector<std::uint32_t> sorted = boundaries;
    std::sort(sorted.begin(), sorted.end());
    sorted.erase(std::unique(sorted.begin(), sorted.end()), sorted.end());
    if (sorted.size() < 2 || sorted.front() == sorted.back()) {
        return std::nullopt;
    }
    return std::make_pair(sorted.front(), sorted.back());
}

[[maybe_unused]] void collect_cycle_snarl_candidates(
    const Graph& graph,
    const PackedGraph& packed,
    const BubbleCallOptions& options,
    CandidateMap& dedup,
    std::vector<SnarlDebugEntry>* debug_entries) {

    if (packed.node_ids.empty()) {
        return;
    }

    const auto reference_path_idx = pick_reference_path_index_for_snarl_macro(graph, options);
    const auto ref_positions_by_node =
        build_reference_positions_by_node(graph, packed, reference_path_idx);

    const HandleGraph hg = build_handle_graph(packed);
    const SccResult scc = compute_handle_sccs(hg);
    if (scc.component_count <= 0) {
        return;
    }

    std::vector<std::vector<std::uint32_t>> handles_by_comp(static_cast<std::size_t>(scc.component_count));
    for (std::uint32_t h = 0; h < scc.component_of_handle.size(); ++h) {
        const int comp = scc.component_of_handle[h];
        if (comp < 0) {
            continue;
        }
        handles_by_comp[static_cast<std::size_t>(comp)].push_back(h);
    }

    std::size_t candidate_id = 0;
    for (const auto& comp_handles : handles_by_comp) {
        if (comp_handles.empty()) {
            continue;
        }
        ++candidate_id;

        SnarlDebugEntry debug;
        debug.candidate_id = candidate_id;
        debug.mode = "snarl";
        debug.component_handle_count = comp_handles.size();

        bool has_cycle = comp_handles.size() > 1;
        if (!has_cycle) {
            const std::uint32_t h = comp_handles.front();
            for (const std::uint32_t nxt : hg.out[h]) {
                if (nxt == h) {
                    has_cycle = true;
                    break;
                }
            }
        }
        if (!has_cycle) {
            if (debug_entries != nullptr) {
                debug.accepted = false;
                debug.reason = "reject:not-cyclic";
                debug_entries->push_back(std::move(debug));
            }
            continue;
        }

        std::unordered_set<std::uint32_t> comp_nodes;
        comp_nodes.reserve(comp_handles.size() * 2);
        std::unordered_set<std::uint32_t> boundary_nodes;
        boundary_nodes.reserve(8);

        for (const std::uint32_t h : comp_handles) {
            const int comp_id = scc.component_of_handle[h];
            const std::uint32_t node_idx = h / 2U;
            comp_nodes.insert(node_idx);

            bool touches_outside = false;
            for (const std::uint32_t v : hg.out[h]) {
                if (scc.component_of_handle[v] != comp_id) {
                    touches_outside = true;
                    break;
                }
            }
            if (!touches_outside) {
                for (const std::uint32_t u : hg.in[h]) {
                    if (scc.component_of_handle[u] != comp_id) {
                        touches_outside = true;
                        break;
                    }
                }
            }
            if (touches_outside) {
                boundary_nodes.insert(node_idx);
            }
        }
        debug.component_node_count = comp_nodes.size();
        debug.boundary_node_count = boundary_nodes.size();

        std::vector<std::uint32_t> boundaries;
        boundaries.reserve(boundary_nodes.size());
        for (const std::uint32_t n : boundary_nodes) {
            boundaries.push_back(n);
        }
        std::sort(boundaries.begin(), boundaries.end());
        std::vector<std::string> boundary_node_ids;
        boundary_node_ids.reserve(boundaries.size());
        for (const std::uint32_t idx : boundaries) {
            boundary_node_ids.push_back(packed.node_ids[idx]);
        }
        debug.boundary_nodes = join_strings(boundary_node_ids, ';');

        if (boundaries.size() < 2) {
            if (debug_entries != nullptr) {
                debug.accepted = false;
                debug.reason = "reject:boundary-count-lt2";
                debug_entries->push_back(std::move(debug));
            }
            continue;
        }

        std::uint32_t endpoint_a = boundaries[0];
        std::uint32_t endpoint_b = boundaries[1];
        bool used_macro_pair = false;
        if (boundaries.size() > 2) {
            const auto macro_pair = choose_macro_boundary_pair(boundaries, ref_positions_by_node);
            if (!macro_pair.has_value()) {
                if (debug_entries != nullptr) {
                    debug.accepted = false;
                    debug.reason = "reject:boundary-count-gt2-no-pair";
                    debug_entries->push_back(std::move(debug));
                }
                continue;
            }
            endpoint_a = macro_pair->first;
            endpoint_b = macro_pair->second;
            used_macro_pair = true;
        }

        if (endpoint_a == endpoint_b) {
            if (debug_entries != nullptr) {
                debug.accepted = false;
                debug.reason = "reject:degenerate-endpoints";
                debug_entries->push_back(std::move(debug));
            }
            continue;
        }

        std::vector<std::uint32_t> inside;
        inside.reserve(comp_nodes.size());
        for (const std::uint32_t node_idx : comp_nodes) {
            if (node_idx == endpoint_a || node_idx == endpoint_b) {
                continue;
            }
            inside.push_back(node_idx);
        }
        debug.inside_node_count = inside.size();
        debug.source = packed.node_ids[endpoint_a];
        debug.sink = packed.node_ids[endpoint_b];

        if (inside.empty()) {
            if (debug_entries != nullptr) {
                debug.accepted = false;
                debug.reason = "reject:empty-inside";
                debug_entries->push_back(std::move(debug));
            }
            continue;
        }
        std::sort(inside.begin(), inside.end());
        inside.erase(std::unique(inside.begin(), inside.end()), inside.end());
        debug.inside_node_count = inside.size();

        BubbleCandidateIdx candidate;
        candidate.source = endpoint_a;
        candidate.sink = endpoint_b;
        candidate.inside = std::move(inside);
        const bool kept = add_or_replace_candidate(dedup, std::move(candidate));
        if (debug_entries != nullptr) {
            debug.accepted = kept;
            if (kept) {
                debug.reason = used_macro_pair ? "accept:macro-boundary-pair" : "accept";
            } else {
                debug.reason = "reject:duplicate-endpoints-smaller-inside";
            }
            debug_entries->push_back(std::move(debug));
        }
    }
}

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

std::string bubble_type_to_string(BubbleType type) {
    switch (type) {
        case BubbleType::Simple:
            return "simple";
        case BubbleType::Insertion:
            return "insertion";
        case BubbleType::Super:
        default:
            return "super";
    }
}

BubbleCallReport call_bubbles_report(const Graph& graph, const BubbleCallOptions& options) {
    const PackedGraph packed = pack_graph(graph);
    CandidateMap dedup;
    dedup.reserve(packed.node_ids.size());

    std::vector<SnarlDebugEntry> snarl_debug;
    auto* debug_ptr = options.collect_snarl_debug ? &snarl_debug : nullptr;

    std::vector<PathIndex> path_indexes;
    path_indexes.reserve(graph.paths.size());
    for (const auto& p : graph.paths) {
        path_indexes.push_back(build_path_index(p));
    }

    collect_snarl_jsonl_candidates(graph, packed, options, dedup, debug_ptr);

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

    const auto compute_bubble_metrics = [&](std::vector<Bubble>& target) {
        for (auto& bubble : target) {
            std::size_t supported_paths = 0;
            std::size_t min_inside_bp = std::numeric_limits<std::size_t>::max();
            std::size_t max_inside_bp = 0;
            std::size_t long_path_support = 0;
            bool has_inside_bp = false;
            bool inversion_signal = false;

            std::unordered_set<std::string> inside_nodes;
            inside_nodes.reserve(bubble.inside.size() * 2);
            for (const auto& n : bubble.inside) {
                inside_nodes.insert(n);
            }
            std::unordered_map<std::string, unsigned char> orientation_mask;
            orientation_mask.reserve(bubble.inside.size() * 2);

            for (std::size_t p_idx = 0; p_idx < path_indexes.size(); ++p_idx) {
                const auto interval = find_best_interval(path_indexes[p_idx], bubble);
                if (!interval.has_value()) {
                    continue;
                }
                ++supported_paths;

                const auto steps = canonical_steps_for_bubble(graph.paths[p_idx], bubble, *interval);
                if (steps.size() < 2) {
                    continue;
                }

                std::size_t inside_bp = 0;
                for (std::size_t i = 1; i + 1 < steps.size(); ++i) {
                    const auto& step = steps[i];
                    if (inside_nodes.find(step.node_id) == inside_nodes.end()) {
                        continue;
                    }

                    const auto node_it = graph.nodes.find(step.node_id);
                    if (node_it != graph.nodes.end()) {
                        inside_bp += node_it->second.sequence.size();
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

            bubble.path_support = supported_paths;
            if (has_inside_bp) {
                bubble.min_inside_bp = min_inside_bp;
                bubble.max_inside_bp = max_inside_bp;
            } else {
                bubble.min_inside_bp = 0;
                bubble.max_inside_bp = 0;
            }
            bubble.long_path_support = long_path_support;
            bubble.inversion_signal = inversion_signal;
        }
    };

    compute_bubble_metrics(bubbles);

    std::vector<Bubble> non_snp_bubbles;
    non_snp_bubbles.reserve(bubbles.size());
    for (const auto& bubble : bubbles) {
        if (bubble.max_inside_bp > 1 || bubble.inversion_signal) {
            non_snp_bubbles.push_back(bubble);
        }
    }

    struct DebugMetrics {
        std::size_t n_paths = 0;
        std::size_t min_inside_bp = 0;
        std::size_t inside_node_count = 0;
    };
    std::unordered_map<std::string, DebugMetrics> metrics_by_endpoint;
    metrics_by_endpoint.reserve(bubbles.size() * 2);
    for (const auto& bubble : bubbles) {
        DebugMetrics metrics;
        metrics.n_paths = bubble.path_support;
        metrics.min_inside_bp = bubble.min_inside_bp;
        metrics.inside_node_count = bubble.inside.size();
        metrics_by_endpoint[endpoint_key(bubble.source, bubble.sink)] = metrics;
    }

    std::unordered_map<std::string, std::string> rejection_reason_by_endpoint;
    rejection_reason_by_endpoint.reserve(bubbles.size() * 2);

    if (options.min_path_support > 0) {
        for (const auto& bubble : bubbles) {
            if (bubble.path_support < options.min_path_support) {
                const std::string key = endpoint_key(bubble.source, bubble.sink);
                if (rejection_reason_by_endpoint.find(key) == rejection_reason_by_endpoint.end()) {
                    rejection_reason_by_endpoint[key] = "reject:min-path-support";
                }
            }
        }
        bubbles.erase(
            std::remove_if(
                bubbles.begin(),
                bubbles.end(),
                [&](const Bubble& b) { return b.path_support < options.min_path_support; }),
            bubbles.end());
    }

    if (options.min_variant_bp > 0) {
        for (const auto& bubble : bubbles) {
            if (bubble.long_path_support == 0 && !bubble.inversion_signal) {
                const std::string key = endpoint_key(bubble.source, bubble.sink);
                if (rejection_reason_by_endpoint.find(key) == rejection_reason_by_endpoint.end()) {
                    rejection_reason_by_endpoint[key] = "reject:min-variant-bp";
                }
            }
        }
        bubbles.erase(
            std::remove_if(
                bubbles.begin(),
                bubbles.end(),
                [&](const Bubble& b) {
                    return b.long_path_support == 0 && !b.inversion_signal;
                }),
            bubbles.end());
    }

    if (options.merge_nearby_bp > 0 && bubbles.size() > 1) {
        std::unordered_set<std::string> premerge_keys;
        premerge_keys.reserve(bubbles.size() * 2);
        for (const auto& bubble : bubbles) {
            premerge_keys.insert(endpoint_key(bubble.source, bubble.sink));
        }

        bubbles = merge_nearby_bubbles(graph, bubbles, options.merge_nearby_bp);
        compute_bubble_metrics(bubbles);

        std::unordered_set<std::string> postmerge_keys;
        postmerge_keys.reserve(bubbles.size() * 2);
        for (const auto& bubble : bubbles) {
            postmerge_keys.insert(endpoint_key(bubble.source, bubble.sink));
        }
        for (const auto& key : premerge_keys) {
            if (postmerge_keys.find(key) != postmerge_keys.end()) {
                continue;
            }
            if (rejection_reason_by_endpoint.find(key) == rejection_reason_by_endpoint.end()) {
                rejection_reason_by_endpoint[key] = "reject:merged-nearby";
            }
        }
    }

    for (auto& bubble : bubbles) {
        if (is_simple_bubble(graph, bubble)) {
            bubble.type = BubbleType::Simple;
        } else if (is_insertion_bubble(graph, bubble)) {
            bubble.type = BubbleType::Insertion;
        } else {
            bubble.type = BubbleType::Super;
        }
    }

    std::sort(bubbles.begin(), bubbles.end(), bubble_endpoint_less);

    assign_ids(bubbles);

    std::unordered_set<std::string> final_endpoint_keys;
    final_endpoint_keys.reserve(bubbles.size() * 2);
    for (const auto& bubble : bubbles) {
        const std::string key = endpoint_key(bubble.source, bubble.sink);
        final_endpoint_keys.insert(key);
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
            }
            const bool in_final = final_endpoint_keys.find(key) != final_endpoint_keys.end();
            if (in_final) {
                if (debug.reason.rfind("accept:", 0) == 0 || debug.reason.empty()) {
                    debug.accepted = true;
                    debug.reason = "accept:final";
                }
                continue;
            }

            if (debug.reason.rfind("accept:", 0) == 0 || debug.reason.empty()) {
                debug.accepted = false;
                const auto reject_it = rejection_reason_by_endpoint.find(key);
                if (reject_it != rejection_reason_by_endpoint.end()) {
                    debug.reason = reject_it->second;
                } else {
                    debug.reason = "reject:filtered-out";
                }
            }
        }

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
            synthetic.accepted = true;
            synthetic.reason = "accept:final-merged";
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
