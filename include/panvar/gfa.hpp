#pragma once

#include <cstddef>
#include <string>
#include <unordered_map>
#include <vector>

namespace panvar {

struct Neighbor {
    std::string node_id;
    // 0 means we enter the neighbor from its start side, 1 from its end side.
    int side = 0;
    int overlap = 0;
};

struct Node {
    std::string id;
    std::string sequence;
    std::vector<Neighbor> start;
    std::vector<Neighbor> end;
};

struct PathStep {
    std::string node_id;
    bool reverse = false;
};

struct PathRecord {
    std::string name;
    char type = 'P';
    std::vector<PathStep> steps;
};

struct Graph {
    std::unordered_map<std::string, Node> nodes;
    std::vector<PathRecord> paths;

    std::vector<std::string> neighbors_of(const std::string& node_id) const;
    bool has_edge_between(const std::string& a, const std::string& b) const;
};

struct ParseGfaOptions {
    bool include_paths = true;
    bool include_sequences = true;
};

Graph parse_gfa(const std::string& gfa_path, const ParseGfaOptions& options = {});

// The two halves of L-line handling, exposed so a Graph built from an in-memory GfaModel cannot drift
// from one parsed out of the file that model would have written. Both `parse_gfa` and
// `graph_from_model` (gfa_io.hpp) go through these, so there is one definition of what a link means.

// '*' on an L line means the overlap is UNKNOWN, not zero: it yields -1. "10M" and "10" both yield 10.
int parse_gfa_overlap(const std::string& field);

// Attach one oriented link, applying the side mapping and the duplicate check. Silently ignores a link
// whose endpoints are not both present, which is what parse_gfa does for an L line that names a
// segment with no S line.
void add_gfa_edge(Graph& graph, const std::string& from, char from_orient,
                  const std::string& to, char to_orient, int overlap);

} // namespace panvar
