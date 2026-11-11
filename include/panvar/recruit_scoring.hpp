#pragma once

// The scoring of one fragment against one locus, shared by recruit's read scorer and the locus index's
// calibration: each read's k-mers labelled (P, S, R, D, N), its locus segment, independent sites, and the
// decision certain / ambiguous / unrecruited with its reason.

#include "panvar/flat_map.hpp"
#include "panvar/locus_index.hpp"
#include "panvar/read_input.hpp"
#include "panvar/recruit_reads.hpp"

#include <cstdint>
#include <string>
#include <vector>

namespace panvar::recruit::scoring {

constexpr std::size_t kK = 31, kS = 11;
constexpr std::uint8_t kLabelP = 0, kLabelS = 1, kLabelR = 2, kLabelD = 3, kLabelN = 4, kMultiBit = 8, kBackboneBit = 16;

// One locus, as reads see it: each k-mer's label (P/S/R, D for decoys), with bits for multi-copy P and
// backbone k-mers (single-copy, in every haplotype: their depth is the sample's depth at the locus); the
// graph nodes each private k-mer lies on.
struct LocusTable {
    std::string name;
    FlatMap<std::uint8_t> label;
    FlatMap<std::uint32_t> nodeset;   // private k-mer -> its node set
    std::vector<std::uint32_t> nodeset_first, nodeset_nodes;
    std::vector<std::string> nodes;
    std::uint64_t backbone = 0;       // backbone k-mers in the index
    std::uint32_t min_sites = 0;      // the private sites certain needs here (0: the run's)
};

struct Row {
    std::string fragment;
    std::uint32_t locus = 0;
    Evidence ev;
    std::string nodes;           // graph nodes of the private k-mers in the segments, "." if none
    std::uint32_t backbone = 0;  // backbone k-mers in the segments
    std::string seg[2];          // each mate's locus segment in read bases, "." when it has none
    const char* mates = "";      // each mate's status: in, part (a segment covering part of it), out
    const char* fate = "";
    const char* reason = "";
};

// What a read's segment search and edge test cost at an error rate.
struct Costs {
    double gap_cost = 0, flank_share = 0;
};
Costs costs_for(double error_rate);

// A locus index as reads see it.
LocusTable table_of(const LocusIndex& x);

// One fragment against one locus.
Row score_fragment(const Fragment& fr, const LocusTable& t, std::uint32_t li, const ReadOptions& o, double gap_cost, double flank_share);

}  // namespace panvar::recruit::scoring
