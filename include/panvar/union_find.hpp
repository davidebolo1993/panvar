#pragma once

// Vendored from vg's dependency `structures` (deps/structures, Apache-2.0):
// a union-find that also supports enumerating the members of a group
// (include_children) and listing all groups. Used by the vendored cactus snarl
// finder (integrated_snarls.cpp). Kept in namespace panvar to avoid collisions.

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <unordered_set>
#include <vector>

namespace panvar {

class UnionFind {
public:
    // group()/all_groups() require include_children = true.
    explicit UnionFind(std::size_t size, bool include_children = true);

    std::size_t size();
    void resize(std::size_t size);
    std::size_t find_group(std::size_t i);
    std::size_t union_groups(std::size_t i, std::size_t j);
    std::size_t group_size(std::size_t i);
    std::vector<std::size_t> group(std::size_t i);
    std::vector<std::vector<std::size_t>> all_groups();

private:
    struct UFNode {
        explicit UFNode(std::size_t index) : rank(0), size(1), head(index) {}
        std::size_t rank;
        std::size_t size;
        std::size_t head;
        std::unordered_set<std::size_t> children;
    };
    std::vector<UFNode> uf_nodes;
    bool include_children;
};

} // namespace panvar
