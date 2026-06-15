#include "panvar/three_edge_cc.hpp"

// Vendored, near-verbatim, from vg src/algorithms/three_edge_connected_components.cpp
// (three_edge_connected_component_merges_dense). Debug output removed; namespace changed.
// Independent implementation of Norouzi and Tsin (2014) "A simple 3-edge connected
// component algorithm revisited".

#include <cassert>
#include <limits>
#include <vector>

namespace panvar {

void three_edge_connected_component_merges_dense(
    std::size_t node_count,
    std::size_t first_root,
    const std::function<void(std::size_t, const std::function<void(std::size_t)>&)>& for_each_connected_node,
    const std::function<void(std::size_t, std::size_t)>& same_component) {

    using number_t = std::size_t;
    assert(node_count < std::numeric_limits<number_t>::max());

    struct TsinNode {
        number_t dfs_counter;
        number_t dfs_exit;
        number_t low_point;
        number_t effective_degree = 0;
        number_t path_tail;
        bool is_on_path;
        bool visited = false;
    };

    std::vector<TsinNode> nodes(node_count);

    // Absorb-eject along a (sub)path; see the original for the full derivation.
    auto absorb_all_along_path = [&](number_t into, number_t path_start, number_t path_past_end) {
        [[maybe_unused]] bool path_null = true;
        number_t here = path_start;
        while (here != path_past_end) {
            if (here == std::numeric_limits<number_t>::max()) {
                assert(path_null);
                break;
            }
            auto& here_node = nodes[here];
            if (here_node.is_on_path) {
                if (into == std::numeric_limits<number_t>::max()) {
                    into = here;
                } else {
                    path_null = false;
                    nodes[into].effective_degree =
                        (nodes[into].effective_degree + here_node.effective_degree - 2);
                    same_component(into, here);
                }
            }
            here = here_node.path_tail;
        }
    };

    struct DFSStackFrame {
        number_t current;
        std::vector<number_t> neighbors;
        bool saw_parent_tree_edge = false;
        bool recursing = false;
    };

    std::vector<DFSStackFrame> stack;
    number_t next_unvisited = 0;
    number_t dfs_counter = 1;

    while (next_unvisited != node_count) {
        if (!nodes[first_root].visited) {
            stack.emplace_back();
            stack.back().current = first_root;
        } else {
            stack.emplace_back();
            stack.back().current = next_unvisited;
        }

        while (!stack.empty()) {
            auto& frame = stack.back();
            auto& node = nodes[frame.current];

            if (!node.visited) {
                node.visited = true;

                if (frame.current == next_unvisited) {
                    do {
                        next_unvisited++;
                    } while (next_unvisited != node_count && nodes[next_unvisited].visited);
                }

                node.dfs_counter = dfs_counter;
                dfs_counter++;
                node.low_point = node.dfs_counter;
                node.path_tail = std::numeric_limits<number_t>::max();
                node.is_on_path = true;

                for_each_connected_node(frame.current, [&](std::size_t connected) {
                    frame.neighbors.push_back(connected);
                });
                continue;
            } else {
                if (!frame.neighbors.empty()) {
                    number_t neighbor_number = frame.neighbors.back();
                    auto& neighbor = nodes[neighbor_number];

                    if (!frame.recursing) {
                        node.effective_degree++;

                        if (!neighbor.visited) {
                            frame.recursing = true;
                            stack.emplace_back();
                            stack.back().current = neighbor_number;
                        } else {
                            if (stack.size() > 1 && neighbor_number == stack[stack.size() - 2].current &&
                                !frame.saw_parent_tree_edge) {
                                // tree edge in from parent
                                frame.saw_parent_tree_edge = true;
                            } else if (neighbor.dfs_counter < node.dfs_counter) {
                                // outgoing back-edge (step 1.2)
                                if (neighbor.dfs_counter < node.low_point) {
                                    absorb_all_along_path(std::numeric_limits<number_t>::max(),
                                                          frame.current,
                                                          std::numeric_limits<number_t>::max());
                                    node.low_point = neighbor.dfs_counter;
                                    node.is_on_path = true;
                                    node.path_tail = std::numeric_limits<number_t>::max();
                                }
                            } else if (node.dfs_counter < neighbor.dfs_counter) {
                                // incoming back-edge (step 1.3)
                                node.effective_degree -= 2;
                                number_t replacement_neighbor_number = frame.current;
                                number_t candidate = nodes[replacement_neighbor_number].path_tail;
                                while (candidate != std::numeric_limits<number_t>::max() &&
                                       nodes[candidate].dfs_counter <= neighbor.dfs_counter &&
                                       nodes[candidate].dfs_exit >= neighbor.dfs_exit) {
                                    replacement_neighbor_number = candidate;
                                    candidate = nodes[replacement_neighbor_number].path_tail;
                                }
                                auto& replacement_neighbor = nodes[replacement_neighbor_number];
                                absorb_all_along_path(std::numeric_limits<number_t>::max(),
                                                      frame.current,
                                                      replacement_neighbor.path_tail);
                                node.path_tail = replacement_neighbor.path_tail;
                            } else {
                                // self loop: censor the edge
                                node.effective_degree--;
                            }
                            frame.neighbors.pop_back();
                        }
                    } else {
                        // returned from recursion on neighbor
                        if (neighbor.low_point == neighbor.dfs_counter) {
                            // bridge edge: hide it
                            neighbor.effective_degree--;
                            node.effective_degree--;
                        } else {
                            if (neighbor.effective_degree == 2) {
                                neighbor.is_on_path = false;
                            }
                            assert(neighbor.effective_degree != 1);

                            if (node.low_point <= neighbor.low_point) {
                                absorb_all_along_path(frame.current,
                                                      neighbor_number,
                                                      std::numeric_limits<number_t>::max());
                            } else {
                                node.low_point = neighbor.low_point;
                                absorb_all_along_path(std::numeric_limits<number_t>::max(),
                                                      frame.current,
                                                      std::numeric_limits<number_t>::max());
                                node.is_on_path = true;
                                node.path_tail = neighbor_number;
                            }
                        }
                        frame.recursing = false;
                        frame.neighbors.pop_back();
                    }
                } else {
                    node.dfs_exit = dfs_counter;
                    stack.pop_back();
                }
            }
        }
    }
}

} // namespace panvar
