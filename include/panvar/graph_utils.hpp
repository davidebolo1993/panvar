#pragma once

#include "panvar/gfa.hpp"

#include <cstddef>
#include <cstdint>
#include <numeric>
#include <utility>
#include <string>
#include <unordered_map>
#include <vector>

namespace panvar {

// Union-find with path compression + union by rank. Shared by the walk clustering in
// `inspect` and the connected-component event merging in `call`.
struct DisjointSet {
    explicit DisjointSet(std::size_t n) : parent(n), rank(n, 0) {
        std::iota(parent.begin(), parent.end(), std::size_t{0});
    }
    std::size_t find(std::size_t x) {
        while (parent[x] != x) {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        return x;
    }
    void unite(std::size_t a, std::size_t b) {
        a = find(a);
        b = find(b);
        if (a == b) {
            return;
        }
        if (rank[a] < rank[b]) {
            std::swap(a, b);
        }
        parent[b] = a;
        if (rank[a] == rank[b]) {
            ++rank[a];
        }
    }
    std::vector<std::size_t> parent;
    std::vector<unsigned char> rank;
};

std::unordered_map<std::string, const PathRecord*> path_records_by_name(const Graph& graph);
std::vector<std::size_t> path_prefix_bp(
    const PathRecord& path,
    const std::unordered_map<std::string, Node>& nodes);
std::string reverse_complement(const std::string& sequence);
std::string spell_path_steps_sequence(
    const Graph& graph,
    const std::vector<PathStep>& steps,
    bool* complete = nullptr);

// Walk identity helpers shared across modules: a stable per-step hash token, the
// token vector for a walk, and a human-readable "id+/id-" signature string.
std::uint64_t hash_step_token(const PathStep& step);
std::vector<std::uint64_t> build_walk_tokens(const std::vector<PathStep>& steps);
std::string build_walk_signature(const std::vector<PathStep>& steps);

} // namespace panvar
