#pragma once

// Vendored from vg (src/algorithms/three_edge_connected_components.cpp), the dense
// "merges" entry point only. Independent implementation of Norouzi & Tsin (2014),
// "A simple 3-edge connected component algorithm revisited". Self-contained: it reports,
// via `same_component`, pairs of nodes that belong to the same 3-edge-connected component
// (a spanning set per component); feed those into a union-find to recover the components.
// We use this to build the cactus graph for internal snarl finding (see snarls.cpp).

#include <cstddef>
#include <functional>
#include <unordered_map>
#include <vector>

namespace panvar {

// node ids are dense integers 0..node_count-1. `for_each_connected_node(n, emit)` must emit
// every node adjacent to n (each undirected edge from both endpoints). `same_component(a,b)`
// is called for pairs to be unioned.
void three_edge_connected_component_merges_dense(
    std::size_t node_count,
    std::size_t first_root,
    const std::function<void(std::size_t, const std::function<void(std::size_t)>&)>& for_each_connected_node,
    const std::function<void(std::size_t, std::size_t)>& same_component);

// Templated wrapper (vendored from vg's header): ranks arbitrary hashable nodes into a
// dense space, runs the dense algorithm, and reports same-component pairs in node space.
template <typename TECCNode>
void three_edge_connected_component_merges(
    const std::function<void(const std::function<void(TECCNode)>&)>& for_each_node,
    const std::function<void(TECCNode, const std::function<void(TECCNode)>&)>& for_each_connected_node,
    const std::function<void(TECCNode, TECCNode)>& same_component) {

    std::vector<TECCNode> rank_to_node;
    std::unordered_map<TECCNode, std::size_t> node_to_rank;
    for_each_node([&](TECCNode node) {
        node_to_rank[node] = rank_to_node.size();
        rank_to_node.push_back(node);
    });

    three_edge_connected_component_merges_dense(
        rank_to_node.size(), 0,
        [&](std::size_t rank, const std::function<void(std::size_t)>& visit_connected) {
            for_each_connected_node(rank_to_node[rank], [&](TECCNode connected) {
                visit_connected(node_to_rank[connected]);
            });
        },
        [&](std::size_t a, std::size_t b) { same_component(rank_to_node[a], rank_to_node[b]); });
}

} // namespace panvar
