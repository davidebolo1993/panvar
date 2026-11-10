#pragma once

#include "panvar/gfa.hpp"

#include <cstddef>
#include <cstdint>
#include <numeric>
#include <utility>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace panvar {

// Union-find with path compression + union by rank, for the connected-component event merging in `call`.
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

std::vector<std::size_t> path_prefix_bp(
    const PathRecord& path,
    const std::unordered_map<std::string, Node>& nodes);
std::string reverse_complement(const std::string& sequence);
std::string spell_path_steps_sequence(
    const Graph& graph,
    const std::vector<PathStep>& steps,
    bool* complete = nullptr);

// Walk identity helpers shared across modules: a stable per-step hash token and a human-readable
// "id+/id-" signature string.
std::uint64_t hash_step_token(const PathStep& step);
std::string build_walk_signature(const std::vector<PathStep>& steps);

// Nodes carrying a self-loop edge (a panphorte unit, or any tandem-unit node): these are
// copy-number loci even when a haplotype traverses them only once. Read from the graph's edges, since
// `bubbles.csv` does not record it.
std::unordered_set<std::string> self_loop_nodes(const Graph& graph);

// Resolve a reference path QUERY to an exact path name by match_reference_name (ref_path.hpp), or a W line
// named by its own coordinates. Ambiguity and absence are errors, never resolved by file order.
//
// Every consumer compares path names exactly, so a module that accepts a loose query and then passes
// the QUERY on has silently resolved nothing: boundaries come back unoriented, reference allele support
// reads 0, and nothing says so. Resolve once, at the edge, and pass the resolved name inward.
std::string resolve_reference_path_name(const Graph& graph, const std::string& query,
                                        const std::string& module);

// Refuses a graph with no P/W paths, as every module reads the haplotypes' walks; a file with no segments
// either is named as empty.
void require_graph_paths(std::size_t segments, std::size_t paths, const std::string& module);

// One authoritative check that a graph is safe to reason about, shared by every module that reads
// paths. Each condition, left unchecked, gives a silently WRONG answer rather than a failure: a step
// naming a missing node spells a shorter sequence than the file describes; a step pair with no link
// describes a traversal the graph does not permit; duplicate path names make "which path" depend on
// file order; a non-zero or unknown overlap double-counts bases in every span and identity figure,
// since spelling concatenates whole segments.
//
// `module` prefixes the message. `require_sequences` is for callers that spell or measure bp;
// `require_zero_overlaps` for callers that concatenate.
void validate_graph_paths(const Graph& graph, const std::string& module,
                          bool require_sequences, bool require_zero_overlaps);

} // namespace panvar
