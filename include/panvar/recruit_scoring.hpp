#pragma once

// The scoring of one fragment against one locus, shared by recruit's read scorer and the locus index's
// calibration: each read's k-mers labelled (P, S, R, D, N), its locus segment, independent sites, the score,
// and the decision certain / likely / ambiguous / unrecruited with its reason.

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

// The score of a judged stretch with n callable k-mers, p private and d decoy sites: the log-likelihood ratio of
// two Poisson processes, a locus read against an off-target one, with private sites at r_in and r_out and decoy
// sites at q_in and q_out per k k-mers. Length, density and decoys against private evidence in one number;
// sites cluster on paralog variants, so it ranks fragments and the thresholds come from calibration.
struct Rates {
    double r_in = 0, r_out = 0, q_in = 0, q_out = 0;
    bool known() const { return r_in > 0 && r_out > 0 && q_in > 0 && q_out > 0; }
};
struct Counts {
    std::uint32_t p = 0, d = 0, n = 0;
};
double score(const Counts& c, const Rates& r);
std::int32_t score_bin(double s);

// What decides a fragment under the score, before rates: the pooled stretch (a pair's mates pooled, sites counted
// once) and each mate alone, and which of them may call it. A stretch may call it with a private site, callable
// and not divergent; the pooled one with a mate wholly in the locus (or a single read), a mate alone unless the
// other mate leans to a homolog away from a locus edge (the rescue).
struct Judged {
    Counts pooled, mate[2];
    bool pooled_ok = false, mate_ok[2] = {false, false};
    bool pooled_mates_out = false;   // the pooled stretch could call it, but both mates leave the locus
};
// The deciding score: the best of the stretches that may call the fragment; false when none may.
bool deciding_score(const Judged& j, const Rates& r, double& s);

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
    // The score: the rates of the run's read class, and the tiers' bins. Fragments are decided by the score when
    // `by_score`, else by the private-site cut (the first rules, an uncalibrated index or --min-sites).
    Rates rates;
    bool by_score = false;
    std::int32_t certain_bin = 0, likely_bin = 0;
    bool tiers = false;                // the bins are set (the tier column)
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
    bool scored = false;         // the locus has rates: `score` holds the deciding score (or the pooled one)
    double score = 0;
    const char* tier = ".";      // the tier the deciding score reaches: certain, likely, or "."
};

// What a read's segment search and edge test cost at an error rate.
struct Costs {
    double gap_cost = 0, flank_share = 0;
};
Costs costs_for(double error_rate);

// A locus index as reads see it.
LocusTable table_of(const LocusIndex& x);

// One fragment against one locus. With `certain_through`, also the largest private-site cut (up to
// kMaxCut) at which the fragment is certain by the first rules, 0 if none: certain at a cut means certain at
// every lower one. With `judged`, what the score decides on.
Row score_fragment(const Fragment& fr, const LocusTable& t, std::uint32_t li, const ReadOptions& o, double gap_cost, double flank_share,
                   std::uint32_t* certain_through = nullptr, Judged* judged = nullptr);

}  // namespace panvar::recruit::scoring
