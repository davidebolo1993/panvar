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

} // namespace panvar
