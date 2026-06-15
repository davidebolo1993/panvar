#include "panvar/union_find.hpp"

// Vendored near-verbatim from vg's `structures` dependency (Apache-2.0).

namespace panvar {

UnionFind::UnionFind(std::size_t size, bool include_children) : include_children(include_children) {
    uf_nodes.reserve(size);
    for (std::size_t i = 0; i < size; i++) {
        uf_nodes.emplace_back(i);
    }
}

std::size_t UnionFind::size() { return uf_nodes.size(); }

void UnionFind::resize(std::size_t size) {
    std::size_t start_size = uf_nodes.size();
    if (size <= start_size) {
        return;
    }
    for (std::size_t i = start_size; i < size; i++) {
        uf_nodes.emplace_back(i);
    }
}

std::size_t UnionFind::find_group(std::size_t i) {
    std::vector<std::size_t> path;
    while (uf_nodes[i].head != i) {
        path.push_back(i);
        i = uf_nodes[i].head;
    }
    std::unordered_set<std::size_t>& head_children = uf_nodes[i].children;
    for (std::size_t p = 1; p < path.size(); p++) {
        std::size_t j = path[p - 1];
        uf_nodes[j].head = i;
        if (include_children) {
            uf_nodes[path[p]].children.erase(j);
            head_children.insert(j);
        }
    }
    return i;
}

std::size_t UnionFind::union_groups(std::size_t i, std::size_t j) {
    std::size_t head_i = find_group(i);
    std::size_t head_j = find_group(j);
    if (head_i == head_j) {
        return head_i;
    }
    UFNode& node_i = uf_nodes[head_i];
    UFNode& node_j = uf_nodes[head_j];
    if (node_i.rank > node_j.rank) {
        node_j.head = head_i;
        node_i.size += node_j.size;
        if (include_children) {
            node_i.children.insert(head_j);
        }
        return head_i;
    } else {
        node_i.head = head_j;
        node_j.size += node_i.size;
        if (include_children) {
            node_j.children.insert(head_i);
        }
        if (node_j.rank == node_i.rank) {
            node_j.rank++;
        }
        return head_j;
    }
}

std::size_t UnionFind::group_size(std::size_t i) { return uf_nodes[find_group(i)].size; }

std::vector<std::size_t> UnionFind::group(std::size_t i) {
    assert(include_children);
    std::vector<std::size_t> to_return;
    std::vector<std::size_t> stack{find_group(i)};
    while (!stack.empty()) {
        std::size_t curr = stack.back();
        stack.pop_back();
        to_return.push_back(curr);
        for (std::size_t child : uf_nodes[curr].children) {
            stack.push_back(child);
        }
    }
    return to_return;
}

std::vector<std::vector<std::size_t>> UnionFind::all_groups() {
    std::vector<std::vector<std::size_t>> to_return(uf_nodes.size());
    for (std::size_t i = 0; i < uf_nodes.size(); i++) {
        to_return[find_group(i)].push_back(i);
    }
    auto new_end = std::remove_if(to_return.begin(), to_return.end(),
                                  [](const std::vector<std::size_t>& grp) { return grp.empty(); });
    to_return.resize(new_end - to_return.begin());
    return to_return;
}

} // namespace panvar
