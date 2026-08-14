#pragma once

// Editable in-memory GFA model with faithful read/write, used by the panphorte
// normalization module (which mutates nodes/edges/paths and writes a new GFA).
// `parse_gfa` (gfa.hpp) builds a read-only adjacency Graph for traversal; this
// model instead preserves the raw S/L/P/W records so the graph can be rewritten
// and serialized back without losing edge orientation, overlaps, or W metadata.

#include "panvar/gfa.hpp"

#include <cstddef>
#include <string>
#include <unordered_map>
#include <vector>

namespace panvar {

struct GfaEdge {
    std::string from;
    char from_orient = '+';
    std::string to;
    char to_orient = '+';
    std::string overlap = "0M"; // raw overlap field (e.g. "0M" or "*")
};

// The one canonical name for a path, matching what parse_gfa() records in Graph::paths: a W line is
// sample#hap#seqid:start-end, a P line is its name. Modules kept private copies of this rule and one of
// them dropped the :start-end suffix, so a W-line graph resolved a reference that then matched no path.
std::string gfa_path_name(const struct GfaPath& p);

struct GfaPath {
    char type = 'P';                 // 'P' or 'W'
    std::vector<PathStep> steps;     // oriented node walk
    // P-line fields:
    std::string name;                // path name (P only)
    std::string overlaps = "*";      // optional P overlaps column
    // W-line fields:
    std::string sample;
    std::string hap;
    std::string seqid;
    std::string start;
    std::string end;
};

struct GfaModel {
    std::vector<std::string> header_lines;             // raw H lines (verbatim)
    std::vector<std::string> node_order;               // node ids in input order
    std::unordered_map<std::string, std::string> seq;  // id -> sequence ("" means '*')
    std::vector<GfaEdge> edges;
    std::vector<GfaPath> paths;

    bool has_node(const std::string& id) const { return seq.find(id) != seq.end(); }
};

// Read a GFA into the editable model (S/L/P/W + H header lines).
GfaModel read_gfa_model(const std::string& gfa_path);

// Serialize the model: H lines, then S (node_order), then L (edges), then P/W.
void write_gfa_model(const std::string& out_path, const GfaModel& model);

// Allocate a fresh numeric node id not present in the model, registering it with
// the given sequence (appended to node_order). Returns the new id as a string.
std::string add_new_node(GfaModel& model, const std::string& sequence);

} // namespace panvar
