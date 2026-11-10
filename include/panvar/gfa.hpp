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

// A P line's overlap column: '*', or one CIGAR per link that is a zero-length match (0M). panvar spells a
// path by concatenating its segments, so an overlap would duplicate bases, and is refused here as on an L
// line. Throws with the line number.
void check_p_overlaps(const std::string& field, std::size_t line_no);

// The name a W line's walk goes by: sample#hap#seqid:start-end, with the interval written as every
// other path name in panvar is read (ref_path.cpp): 1-based, inclusive. A W line's own start is
// 0-based and its end exclusive (GFA 1.1), so the start is shifted by one; read as written, a walk
// starting at 0 would not be a PanSN interval and every other position would be one base early. A start
// or end that is not a number ('*') leaves the interval out.
std::string w_path_name(const std::string& sample, const std::string& hap, const std::string& seqid,
                        const std::string& start, const std::string& end);

// A query naming a W line by its own coordinates ("sample#hap#seq:1000-1221", the W line's 0-based
// start) as panvar names it ("...:1001-1221"); empty when the query ends in no :start-end.
std::string w_literal_as_path_name(const std::string& query);

} // namespace panvar
