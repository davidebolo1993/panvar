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
constexpr std::size_t kLongRead = 1000;   // a read this long or longer is held to the long-read decoy share
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
    double private_share = 0;         // private k-mers among the locus's k-mers
    // The long-read rules (see long_read_rules), set by a run for its read class; off in calibration, which the
    // cuts were set by.
    bool dense_conflicts = false;     // a conflict read with dense private sites is certain
    std::uint32_t min_end_bases = 0;  // a read certain on a segment that leaves the locus needs this many bases of it
    double run_error = 0;             // the error rate the dense rule holds mismatch runs to
};

// The long-read rules. Dense: a conflict read (private and decoy sites) is certain when private k-mers are at least
// kDenseShare of its k-mers, with kDenseSites private sites or more, decoy sites in its segment at most twice its
// private sites plus two, and no more mismatch runs than the run's error rate explains; hifi and ont runs, at loci
// whose k-mers are at least kDenseLocusShare private (below that, paralog reads are as dense as the locus's own).
// End: in a hifi run, a read certain on a segment that leaves the locus needs kMinEndBases of it (a short end stretch
// carries chance sites).
constexpr double kDenseShare = 0.015;
constexpr std::uint32_t kDenseSites = 5;
constexpr double kDenseLocusShare = 0.5;
constexpr std::uint32_t kMinEndBases = 500;

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

// One fragment against one locus. With `certain_through`, also the largest private-site cut (up to
// kMaxCut) at which the fragment is certain, 0 if none: certain at a cut means certain at every lower one.
Row score_fragment(const Fragment& fr, const LocusTable& t, std::uint32_t li, const ReadOptions& o, double gap_cost, double flank_share,
                   std::uint32_t* certain_through = nullptr);

}  // namespace panvar::recruit::scoring
