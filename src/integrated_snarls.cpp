#include "panvar/integrated_snarls.hpp"

// Port of vg's integrated_snarl_finder (MIT) - just the traverse_decomposition path, on an in-process
// HandleShim. Cactus / 3-edge-connected decomposition unchanged, so snarls match `vg snarls`.
// Method: Paten et al. 2018, Superbubbles, Ultrabubbles and Cacti (https://doi.org/10.1089/cmb.2017.0251).
// Upstream: https://github.com/vgteam/vg/blob/master/src/integrated_snarl_finder.cpp

#include "panvar/three_edge_cc.hpp"
#include "panvar/union_find.hpp"

#include <algorithm>
#include <unordered_map>
#include <cstdint>
#include <functional>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace panvar {
namespace {

// A handle is a union-find index: node_rank * 2 + reverse_bit. This makes uf_rank(handle)
// == handle, exactly matching vg's RankedHandleGraph indexing convention.
using handle_t = std::uint64_t;
using edge_t = std::pair<handle_t, handle_t>;

// Thin stand-in for vg's RankedHandleGraph over SnarlGraphInput.
struct HandleShim {
    const SnarlGraphInput* in;

    explicit HandleShim(const SnarlGraphInput* in) : in(in) {}

    handle_t get_handle(std::size_t rank, bool rev) const {
        return (static_cast<handle_t>(rank) << 1) | (rev ? 1u : 0u);
    }
    bool get_is_reverse(handle_t h) const { return (h & 1u) != 0; }
    handle_t flip(handle_t h) const { return h ^ 1u; }
    handle_t forward(handle_t h) const { return h & ~static_cast<handle_t>(1u); }
    std::size_t get_node_count() const { return in->node_ids.size(); }
    std::size_t get_length(handle_t h) const { return in->node_len[h >> 1]; }
    const std::string& get_id(handle_t h) const { return in->node_ids[h >> 1]; }
    // 1-based rank in union-find space, matching vg (uf_rank subtracts 1).
    std::size_t handle_to_rank(handle_t h) const { return static_cast<std::size_t>(h) + 1; }
    handle_t rank_to_handle(std::size_t rank) const { return static_cast<handle_t>(rank - 1); }

    void for_each_edge(const std::function<void(const edge_t&)>& iteratee) const {
        for (const auto& e : in->edges) {
            iteratee(edge_t{get_handle(e.from_rank, e.from_rev), get_handle(e.to_rank, e.to_rev)});
        }
    }
};

// ---- MergedAdjacencyGraph (ported) -------------------------------------------------------
class MergedAdjacencyGraph {
protected:
    const HandleShim* graph;
    mutable UnionFind union_find;

    std::size_t uf_rank(handle_t into) const { return graph->handle_to_rank(into) - 1; }
    handle_t uf_handle(std::size_t rank) const { return graph->rank_to_handle(rank + 1); }
    std::size_t get_weighted_length(const handle_t& handle) const { return graph->get_length(handle); }

public:
    explicit MergedAdjacencyGraph(const HandleShim* graph)
        : graph(graph), union_find(graph->get_node_count() * 2, true) {
        graph->for_each_edge([&](const edge_t& e) {
            auto into_b = graph->flip(e.second);
            merge(e.first, into_b);
        });
    }

    MergedAdjacencyGraph(const MergedAdjacencyGraph& other) : MergedAdjacencyGraph(other.graph) {
        other.for_each_membership([&](handle_t head, handle_t member) { merge(head, member); });
    }

    void merge(handle_t into_a, handle_t into_b) {
        union_find.union_groups(uf_rank(into_a), uf_rank(into_b));
    }

    handle_t find(handle_t into) const { return uf_handle(union_find.find_group(uf_rank(into))); }

    void for_each_head(const std::function<void(handle_t)>& iteratee) const {
        std::vector<bool> seen_heads(union_find.size(), false);
        for (std::size_t i = 0; i < union_find.size(); i++) {
            if (!seen_heads[i]) {
                std::size_t head = union_find.find_group(i);
                if (!seen_heads[head]) {
                    seen_heads[head] = true;
                    iteratee(uf_handle(head));
                }
            }
        }
    }

    void for_each_other_member(handle_t head, const std::function<void(handle_t)>& iteratee) const {
        std::size_t head_rank = uf_rank(head);
        std::vector<std::size_t> group = union_find.group(head_rank);
        for (auto& member_rank : group) {
            if (member_rank != head_rank) iteratee(uf_handle(member_rank));
        }
    }

    void for_each_member(handle_t head, const std::function<void(handle_t)>& iteratee) const {
        std::size_t head_rank = uf_rank(head);
        std::vector<std::size_t> group = union_find.group(head_rank);
        for (auto& member_rank : group) iteratee(uf_handle(member_rank));
    }

    void for_each_membership(const std::function<void(handle_t, handle_t)>& iteratee) const {
        std::vector<std::vector<std::size_t>> uf_components = union_find.all_groups();
        for (auto& component : uf_components) {
            for (std::size_t i = 1; i < component.size(); i++) {
                iteratee(uf_handle(component[0]), uf_handle(component[i]));
            }
        }
    }

    std::pair<std::vector<std::pair<std::size_t, handle_t>>, std::unordered_map<handle_t, handle_t>>
    cycles_in_cactus() const {
        std::pair<std::vector<std::pair<std::size_t, handle_t>>, std::unordered_map<handle_t, handle_t>>
            to_return;
        auto& longest_cycles = to_return.first;
        auto& next_edge = to_return.second;

        std::unordered_map<handle_t, std::size_t> visited_frame;

        struct DFSFrame {
            handle_t here;
            std::vector<handle_t> todo;
        };
        std::vector<DFSFrame> stack;

        for_each_head([&](handle_t component_root) {
            if (!visited_frame.count(component_root)) {
                stack.emplace_back();
                stack.back().here = component_root;

                longest_cycles.emplace_back();
                auto& longest_cycle = longest_cycles.back();

                while (!stack.empty()) {
                    auto& frame = stack.back();
                    auto frame_head = find(frame.here);

                    auto frame_it = visited_frame.find(frame_head);
                    if (frame_it == visited_frame.end()) {
                        frame_it = visited_frame.emplace_hint(frame_it, frame_head, stack.size() - 1);
                        for_each_member(frame_head, [&](handle_t member) {
                            if (member != frame.here || stack.size() == 1) {
                                frame.todo.push_back(graph->flip(member));
                            }
                        });
                    }

                    if (!frame.todo.empty()) {
                        handle_t edge_into = frame.todo.back();
                        handle_t connected_head = find(edge_into);
                        frame.todo.pop_back();

                        auto connected_it = visited_frame.find(connected_head);
                        if (connected_it == visited_frame.end()) {
                            stack.emplace_back();
                            stack.back().here = edge_into;
                        } else {
                            if (frame_it->second > connected_it->second) {
                                std::size_t cycle_length_bp = get_weighted_length(edge_into);
                                handle_t prev_edge = edge_into;
                                for (std::size_t i = connected_it->second + 1; i < stack.size(); i++) {
                                    cycle_length_bp += get_weighted_length(stack[i].here);
                                    next_edge[prev_edge] = stack[i].here;
                                    prev_edge = stack[i].here;
                                }
                                next_edge[prev_edge] = edge_into;

                                if (cycle_length_bp > longest_cycle.first) {
                                    longest_cycle.first = cycle_length_bp;
                                    longest_cycle.second = edge_into;
                                }
                            }
                        }
                    } else {
                        stack.pop_back();
                    }
                }

                if (longest_cycle.first == 0) {
                    longest_cycles.pop_back();
                }
            }
        });

        return to_return;
    }

    std::vector<handle_t> find_cycle_path_in_cactus(
        const std::unordered_map<handle_t, handle_t>& next_along_cycle, handle_t start_head,
        handle_t end_head) const {
        std::vector<handle_t> cycle_path;
        std::vector<std::tuple<handle_t, std::vector<handle_t>, bool>> cycle_stack;

        std::vector<handle_t> roots;
        for_each_member(start_head, [&](handle_t inbound) {
            if (next_along_cycle.count(inbound)) roots.push_back(inbound);
        });

        for (auto& root : roots) {
            cycle_stack.emplace_back(root, std::vector<handle_t>(), false);
            while (!cycle_stack.empty()) {
                auto& cycle_frame = cycle_stack.back();
                if (!std::get<2>(cycle_frame)) {
                    std::get<2>(cycle_frame) = true;
                    for (auto it = next_along_cycle.find(std::get<0>(cycle_frame));
                         it->second != std::get<0>(cycle_frame); it = next_along_cycle.find(it->second)) {
                        handle_t node = find(it->second);
                        if (node == end_head) {
                            cycle_path.reserve(cycle_stack.size());
                            for (auto& f : cycle_stack) cycle_path.push_back(std::get<0>(f));
                            return cycle_path;
                        }
                        for_each_member(node, [&](handle_t inbound) {
                            if (inbound != it->second && next_along_cycle.count(inbound)) {
                                std::get<1>(cycle_frame).push_back(inbound);
                            }
                        });
                    }
                }
                if (!std::get<1>(cycle_frame).empty()) {
                    handle_t child = std::get<1>(cycle_frame).back();
                    std::get<1>(cycle_frame).pop_back();
                    cycle_stack.emplace_back(child, std::vector<handle_t>(), false);
                } else {
                    cycle_stack.pop_back();
                }
            }
        }
        throw std::runtime_error("Could not find cycle path!");
    }

    std::pair<std::vector<std::pair<std::size_t, std::vector<handle_t>>>, std::unordered_map<handle_t, handle_t>>
    longest_paths_in_forest(const std::vector<std::pair<std::size_t, handle_t>>& longest_simple_cycles) const {
        std::pair<std::vector<std::pair<std::size_t, std::vector<handle_t>>>, std::unordered_map<handle_t, handle_t>>
            to_return;
        auto& longest_tree_paths = to_return.first;
        auto& deepest_child_edge = to_return.second;

        struct DFSRecord {
            handle_t parent_edge;
            std::size_t leaf_path_length = 0;
            handle_t second_deepest_child_edge;
            bool has_second_deepest_child = false;
            handle_t longest_subtree_path_root;
            std::size_t longest_subtree_path_length = 0;
        };
        std::unordered_map<handle_t, DFSRecord> records;

        struct DFSFrame {
            handle_t here;
            std::vector<handle_t> todo;
        };
        std::vector<DFSFrame> stack;

        auto try_root = [&](handle_t traversal_root, std::size_t root_cycle_length) {
            if (!records.count(traversal_root)) {
                stack.emplace_back();
                stack.back().here = traversal_root;

                while (!stack.empty()) {
                    auto& frame = stack.back();
                    auto frame_head = find(frame.here);

                    auto frame_it = records.find(frame_head);
                    if (frame_it == records.end()) {
                        frame_it = records.emplace_hint(frame_it, frame_head, DFSRecord());
                        frame_it->second.parent_edge = graph->flip(frame.here);
                        frame_it->second.longest_subtree_path_root = frame_head;

                        for_each_member(frame_head, [&](handle_t member) {
                            auto flipped = graph->flip(member);
                            if (find(flipped) != frame_head) {
                                frame.todo.push_back(flipped);
                            }
                        });
                    }

                    auto& record = frame_it->second;

                    if (!frame.todo.empty()) {
                        handle_t edge_into = frame.todo.back();
                        handle_t connected_head = find(edge_into);
                        frame.todo.pop_back();

                        if (!records.count(connected_head)) {
                            stack.emplace_back();
                            stack.back().here = edge_into;
                        }
                    } else {
                        auto deepest_child_edge_it = deepest_child_edge.find(frame_head);

                        if (stack.size() > 1) {
                            auto& parent_frame = stack[stack.size() - 2];
                            auto parent_head = find(parent_frame.here);
                            auto& parent_record = records[parent_head];

                            record.leaf_path_length = get_weighted_length(frame.here);
                            if (deepest_child_edge_it != deepest_child_edge.end()) {
                                record.leaf_path_length +=
                                    records[find(deepest_child_edge_it->second)].leaf_path_length;
                            }

                            auto parent_deepest_child_it = deepest_child_edge.find(parent_head);
                            if (parent_deepest_child_it == deepest_child_edge.end()) {
                                deepest_child_edge.emplace_hint(parent_deepest_child_it, parent_head,
                                                               frame.here);
                            } else if (records[find(parent_deepest_child_it->second)].leaf_path_length <
                                       record.leaf_path_length) {
                                parent_record.second_deepest_child_edge = parent_deepest_child_it->second;
                                parent_record.has_second_deepest_child = true;
                                parent_deepest_child_it->second = frame.here;
                            } else if (!parent_record.has_second_deepest_child) {
                                parent_record.second_deepest_child_edge = frame.here;
                                parent_record.has_second_deepest_child = true;
                            } else if (records[find(parent_record.second_deepest_child_edge)]
                                           .leaf_path_length < record.leaf_path_length) {
                                parent_record.second_deepest_child_edge = frame.here;
                            }
                        }

                        if (record.has_second_deepest_child || stack.size() == 1) {
                            std::size_t longest_here_path_length = 0;
                            if (deepest_child_edge_it != deepest_child_edge.end()) {
                                longest_here_path_length +=
                                    records[find(deepest_child_edge_it->second)].leaf_path_length;
                            }
                            if (record.has_second_deepest_child) {
                                longest_here_path_length +=
                                    records[find(record.second_deepest_child_edge)].leaf_path_length;
                            }

                            if (record.longest_subtree_path_root == frame_head ||
                                longest_here_path_length > record.longest_subtree_path_length) {
                                record.longest_subtree_path_root = frame_head;
                                record.longest_subtree_path_length = longest_here_path_length;
                            }
                        }

                        if (stack.size() > 1 && record.longest_subtree_path_length > 0) {
                            auto& parent_frame = stack[stack.size() - 2];
                            auto parent_head = find(parent_frame.here);
                            auto& parent_record = records[parent_head];

                            if (parent_record.longest_subtree_path_root == parent_head ||
                                parent_record.longest_subtree_path_length <
                                    record.longest_subtree_path_length) {
                                parent_record.longest_subtree_path_root = record.longest_subtree_path_root;
                                parent_record.longest_subtree_path_length =
                                    record.longest_subtree_path_length;
                            }
                        }

                        if (stack.size() == 1) {
                            if (record.longest_subtree_path_length >= root_cycle_length) {
                                longest_tree_paths.emplace_back();
                                longest_tree_paths.back().first = record.longest_subtree_path_length;
                                auto& path = longest_tree_paths.back().second;

                                auto& path_root_frame = records[record.longest_subtree_path_root];

                                if (path_root_frame.has_second_deepest_child) {
                                    path.push_back(path_root_frame.second_deepest_child_edge);
                                    auto path_trace_it = deepest_child_edge.find(find(path.back()));
                                    while (path_trace_it != deepest_child_edge.end()) {
                                        path.push_back(path_trace_it->second);
                                        path_trace_it = deepest_child_edge.find(find(path.back()));
                                    }
                                    std::vector<handle_t> flipped;
                                    flipped.reserve(path.size());
                                    for (auto path_it = path.rbegin(); path_it != path.rend(); ++path_it) {
                                        flipped.push_back(graph->flip(*path_it));
                                    }
                                    path = std::move(flipped);
                                }

                                if (deepest_child_edge.count(record.longest_subtree_path_root)) {
                                    path.push_back(deepest_child_edge[record.longest_subtree_path_root]);
                                    auto path_trace_it = deepest_child_edge.find(find(path.back()));
                                    while (path_trace_it != deepest_child_edge.end()) {
                                        path.push_back(path_trace_it->second);
                                        path_trace_it = deepest_child_edge.find(find(path.back()));
                                    }
                                }

                                handle_t cursor = record.longest_subtree_path_root;
                                std::vector<handle_t> convergence_to_old_root;
                                while (cursor != frame_head) {
                                    auto& cursor_record = records[cursor];
                                    convergence_to_old_root.push_back(cursor_record.parent_edge);
                                    cursor = find(cursor_record.parent_edge);
                                }

                                while (!convergence_to_old_root.empty()) {
                                    handle_t parent_child_edge = convergence_to_old_root.back();
                                    handle_t child_head = find(parent_child_edge);
                                    handle_t parent_head = find(graph->flip(parent_child_edge));

                                    auto& child_record = records[child_head];
                                    auto& parent_record = records[parent_head];

                                    deepest_child_edge_it = deepest_child_edge.find(child_head);
                                    if (deepest_child_edge_it != deepest_child_edge.end() &&
                                        find(deepest_child_edge_it->second) == parent_head) {
                                        if (child_record.has_second_deepest_child) {
                                            deepest_child_edge_it->second =
                                                child_record.second_deepest_child_edge;
                                            child_record.has_second_deepest_child = false;
                                        } else {
                                            deepest_child_edge.erase(deepest_child_edge_it);
                                            deepest_child_edge_it = deepest_child_edge.end();
                                        }
                                    }

                                    child_record.leaf_path_length = get_weighted_length(parent_child_edge);
                                    if (deepest_child_edge_it != deepest_child_edge.end()) {
                                        child_record.leaf_path_length +=
                                            records[find(deepest_child_edge_it->second)].leaf_path_length;
                                    }

                                    auto parent_deepest_child_it = deepest_child_edge.find(parent_head);
                                    if (parent_deepest_child_it == deepest_child_edge.end()) {
                                        deepest_child_edge.emplace_hint(parent_deepest_child_it,
                                                                       parent_head, parent_child_edge);
                                    } else if (records[find(parent_deepest_child_it->second)]
                                                   .leaf_path_length < child_record.leaf_path_length) {
                                        parent_record.second_deepest_child_edge =
                                            parent_deepest_child_it->second;
                                        parent_record.has_second_deepest_child = true;
                                        parent_deepest_child_it->second = parent_child_edge;
                                    } else if (!parent_record.has_second_deepest_child) {
                                        parent_record.second_deepest_child_edge = parent_child_edge;
                                        parent_record.has_second_deepest_child = true;
                                    } else if (records[find(parent_record.second_deepest_child_edge)]
                                                   .leaf_path_length < child_record.leaf_path_length) {
                                        parent_record.second_deepest_child_edge = parent_child_edge;
                                    }

                                    convergence_to_old_root.pop_back();
                                }

                                if (path.empty()) {
                                    path.push_back(traversal_root);
                                }
                            }
                        }

                        stack.pop_back();
                    }
                }
            }
        };

        for (auto it = longest_simple_cycles.begin(); it != longest_simple_cycles.end(); ++it) {
            try_root(find(it->second), it->first);
        }
        for_each_head([&](handle_t head) { try_root(head, 0); });

        return to_return;
    }
};

// ---- node set ignoring orientation -------------------------------------------------------
class HandleGraphNodeSet {
private:
    std::unordered_set<handle_t> visited;
    const HandleShim* graph;

public:
    explicit HandleGraphNodeSet(const HandleShim* graph) : graph(graph) {}
    std::size_t size() const { return visited.size(); }
    void insert(const handle_t& here) { visited.insert(graph->forward(here)); }
    bool count(const handle_t& here) const { return visited.count(graph->forward(here)); }
};

// ---- core traversal (ported) -------------------------------------------------------------
void traverse_computed_decomposition(
    const HandleShim* graph, MergedAdjacencyGraph& cactus, const MergedAdjacencyGraph& forest,
    std::vector<std::pair<std::size_t, std::vector<handle_t>>>& longest_paths,
    std::unordered_map<handle_t, handle_t>& towards_deepest_leaf,
    std::vector<std::pair<std::size_t, handle_t>>& longest_cycles,
    std::unordered_map<handle_t, handle_t>& next_along_cycle,
    const std::function<void(handle_t)>& begin_chain, const std::function<void(handle_t)>& end_chain,
    const std::function<void(handle_t)>& begin_snarl, const std::function<void(handle_t)>& end_snarl) {

    HandleGraphNodeSet visited(graph);
    std::size_t to_decompose = graph->get_node_count();

    while (visited.size() < to_decompose) {
        struct SnarlChainFrame {
            bool is_snarl = true;
            bool saw_children = false;
            std::pair<handle_t, handle_t> bounds;
            std::vector<handle_t> todo;
        };
        std::vector<SnarlChainFrame> stack;

        if (longest_cycles.empty() ||
            (!longest_paths.empty() && longest_cycles.back().first <= longest_paths.back().first)) {
            if (!visited.count(longest_paths.back().second.front())) {
                handle_t first_edge = longest_paths.back().second.front();

                if (longest_paths.back().first == 0) {
                    cactus.for_each_member(cactus.find(first_edge), [&](handle_t inbound) {
                        if (!graph->get_is_reverse(inbound)) {
                            begin_chain(inbound);
                            end_chain(inbound);
                            visited.insert(inbound);
                        }
                    });
                } else {
                    for (std::size_t i = 1; i < longest_paths.back().second.size(); i++) {
                        handle_t prev_path_edge = longest_paths.back().second[i - 1];
                        handle_t prev_head = forest.find(prev_path_edge);
                        handle_t next_path_edge = longest_paths.back().second[i];
                        towards_deepest_leaf[prev_head] = next_path_edge;
                    }

                    stack.emplace_back();
                    stack.back().is_snarl = true;
                    stack.back().todo.push_back(graph->flip(first_edge));

                    cactus.for_each_member(cactus.find(graph->flip(first_edge)), [&](handle_t inbound) {
                        if (inbound == graph->flip(first_edge)) return;
                        if (next_along_cycle.count(inbound)) {
                            stack.back().todo.push_back(inbound);
                        } else if (cactus.find(inbound) == cactus.find(graph->flip(inbound)) &&
                                   !graph->get_is_reverse(inbound)) {
                            begin_chain(inbound);
                            end_chain(inbound);
                            visited.insert(inbound);
                        }
                    });
                }
            }
            longest_paths.pop_back();
        } else {
            if (!visited.count(longest_cycles.back().second)) {
                stack.emplace_back();
                stack.back().is_snarl = true;

                stack.emplace_back();
                stack.back().is_snarl = false;
                stack.back().bounds =
                    std::make_pair(longest_cycles.back().second, longest_cycles.back().second);
            }
            longest_cycles.pop_back();
        }

        while (!stack.empty()) {
            auto& frame = stack.back();

            if (stack.size() > 1 && !frame.saw_children) {
                frame.saw_children = true;
                (frame.is_snarl ? begin_snarl : begin_chain)(frame.bounds.first);

                if (frame.is_snarl) {
                    visited.insert(frame.bounds.first);
                    visited.insert(frame.bounds.second);

                    cactus.for_each_member(cactus.find(frame.bounds.first), [&](handle_t inbound) {
                        if (inbound == frame.bounds.first || graph->flip(inbound) == frame.bounds.second) {
                            // boundary; stay inside
                        } else if (forest.find(graph->flip(inbound)) != forest.find(inbound)) {
                            frame.todo.push_back(inbound);
                        } else if (next_along_cycle.count(inbound)) {
                            frame.todo.push_back(inbound);
                        } else if (cactus.find(graph->flip(inbound)) == cactus.find(inbound) &&
                                   !graph->get_is_reverse(inbound)) {
                            begin_chain(inbound);
                            end_chain(inbound);
                            visited.insert(inbound);
                        }
                    });
                } else {
                    handle_t here = frame.bounds.first;
                    std::unordered_set<handle_t> seen;
                    std::size_t region_start = frame.todo.size();
                    do {
                        seen.insert(here);
                        frame.todo.push_back(here);
                        here = next_along_cycle.at(here);
                    } while (here != frame.bounds.second);
                    std::reverse(frame.todo.begin() + region_start, frame.todo.end());
                }
            }

            if (!frame.todo.empty()) {
                handle_t task = frame.todo.back();
                frame.todo.pop_back();

                if (frame.is_snarl) {
                    auto next_along_cycle_it = next_along_cycle.find(task);
                    if (next_along_cycle_it != next_along_cycle.end()) {
                        handle_t outgoing = next_along_cycle_it->second;
                        stack.emplace_back();
                        stack.back().is_snarl = false;
                        stack.back().bounds = std::make_pair(outgoing, task);
                    } else {
                        handle_t edge = graph->flip(task);
                        handle_t cactus_head = cactus.find(edge);
                        auto deepest_it = towards_deepest_leaf.find(forest.find(cactus_head));
                        while (deepest_it != towards_deepest_leaf.end()) {
                            handle_t next_back_head = cactus.find(graph->flip(deepest_it->second));

                            if (cactus_head != next_back_head) {
                                std::vector<handle_t> cycle_path = cactus.find_cycle_path_in_cactus(
                                    next_along_cycle, cactus_head, next_back_head);

                                while (!cycle_path.empty()) {
                                    auto through_path_member = next_along_cycle.find(cycle_path.back());
                                    auto through_end = through_path_member;
                                    do {
                                        through_end = next_along_cycle.find(through_end->second);
                                    } while (cactus.find(through_end->first) != cactus.find(next_back_head));

                                    cactus.merge(cycle_path.back(), next_back_head);
                                    std::swap(through_path_member->second, through_end->second);

                                    if (through_path_member->first == through_path_member->second) {
                                        next_along_cycle.erase(through_path_member);
                                    }
                                    if (through_end->first == through_end->second) {
                                        next_along_cycle.erase(through_end);
                                    }
                                    cycle_path.pop_back();
                                }
                            }

                            next_along_cycle[edge] = deepest_it->second;
                            edge = deepest_it->second;
                            cactus_head = cactus.find(edge);
                            deepest_it = towards_deepest_leaf.find(forest.find(cactus_head));
                        }

                        if (edge == graph->flip(task)) {
                            visited.insert(edge);
                            begin_chain(graph->forward(edge));
                            end_chain(graph->forward(edge));
                        } else {
                            next_along_cycle[edge] = graph->flip(task);
                        }

                        cactus.for_each_member(cactus_head, [&](handle_t inbound) {
                            if (next_along_cycle.count(inbound)) {
                                frame.todo.push_back(inbound);
                            } else if (cactus.find(graph->flip(inbound)) == cactus.find(inbound) &&
                                       !graph->get_is_reverse(inbound)) {
                                begin_chain(inbound);
                                end_chain(inbound);
                                visited.insert(inbound);
                            }
                        });

                        cactus.merge(edge, task);
                    }
                } else {
                    handle_t out_edge = next_along_cycle.at(task);
                    stack.emplace_back();
                    stack.back().is_snarl = true;
                    stack.back().bounds = std::make_pair(task, out_edge);
                }
            } else {
                if (stack.size() > 1) {
                    (frame.is_snarl ? end_snarl : end_chain)(frame.bounds.second);
                }
                stack.pop_back();
            }
        }
    }
}

// ---- orchestration (ported from traverse_decomposition, bdsg overlay dropped) ------------
void traverse_decomposition(const HandleShim* graph, const std::function<void(handle_t)>& begin_chain,
                            const std::function<void(handle_t)>& end_chain,
                            const std::function<void(handle_t)>& begin_snarl,
                            const std::function<void(handle_t)>& end_snarl) {
    MergedAdjacencyGraph cactus(graph);

    std::vector<std::pair<handle_t, handle_t>> merge_list;
    three_edge_connected_component_merges<handle_t>(
        [&](const std::function<void(handle_t)>& emit_node) {
            cactus.for_each_head([&](handle_t head) { emit_node(head); });
        },
        [&](handle_t node, const std::function<void(handle_t)>& emit_edge) {
            cactus.for_each_member(node, [&](handle_t other_member) {
                handle_t member_connected_head = cactus.find(graph->flip(other_member));
                if (member_connected_head == node && graph->get_is_reverse(other_member)) {
                    return;
                }
                emit_edge(member_connected_head);
            });
        },
        [&](handle_t a, handle_t b) { merge_list.emplace_back(a, b); });

    for (auto& ab : merge_list) cactus.merge(ab.first, ab.second);
    merge_list.clear();

    MergedAdjacencyGraph forest(cactus);

    auto cycles = cactus.cycles_in_cactus();
    auto& longest_cycles = cycles.first;
    auto& next_along_cycle = cycles.second;

    for (auto& kv : next_along_cycle) forest.merge(kv.first, kv.second);

    auto forest_paths = forest.longest_paths_in_forest(longest_cycles);
    auto& longest_paths = forest_paths.first;
    auto& towards_deepest_leaf = forest_paths.second;

    std::sort(longest_cycles.begin(), longest_cycles.end());
    std::sort(longest_paths.begin(), longest_paths.end());

    traverse_computed_decomposition(graph, cactus, forest, longest_paths, towards_deepest_leaf,
                                    longest_cycles, next_along_cycle, begin_chain, end_chain,
                                    begin_snarl, end_snarl);
}

} // namespace

std::vector<std::pair<std::string, std::string>> find_top_level_snarls_cactus(
    const SnarlGraphInput& in) {
    std::vector<std::pair<std::string, std::string>> out;
    if (in.node_ids.empty()) return out;

    HandleShim shim(&in);

    // Collect the outermost snarls: a snarl with no enclosing snarl is top-level.
    int open_snarls = 0;
    std::string pending_start;
    auto begin_chain = [&](handle_t) {};
    auto end_chain = [&](handle_t) {};
    auto begin_snarl = [&](handle_t h) {
        if (open_snarls == 0) pending_start = shim.get_id(h);
        ++open_snarls;
    };
    auto end_snarl = [&](handle_t h) {
        --open_snarls;
        if (open_snarls == 0) out.emplace_back(pending_start, shim.get_id(h));
    };

    traverse_decomposition(&shim, begin_chain, end_chain, begin_snarl, end_snarl);
    return out;
}

SnarlGraphInput snarl_input_from_model(const GfaModel& model) {
    SnarlGraphInput sgi;
    sgi.node_ids = model.node_order;
    sgi.node_len.resize(model.node_order.size());
    std::unordered_map<std::string, std::size_t> rank_of;
    rank_of.reserve(model.node_order.size() * 2);
    for (std::size_t i = 0; i < model.node_order.size(); ++i) {
        rank_of[model.node_order[i]] = i;
        const auto it = model.seq.find(model.node_order[i]);
        sgi.node_len[i] = (it != model.seq.end()) ? it->second.size() : 0;
    }
    for (const GfaEdge& e : model.edges) {
        const auto fi = rank_of.find(e.from);
        const auto ti = rank_of.find(e.to);
        if (fi == rank_of.end() || ti == rank_of.end()) continue;
        sgi.edges.push_back({fi->second, e.from_orient == '-', ti->second, e.to_orient == '-'});
    }
    return sgi;
}

} // namespace panvar
