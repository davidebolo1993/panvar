#pragma once

#include "panvar/gfa.hpp"

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace panvar {

// Inversions read from the haplotypes' walks: a run of the nodes the reference holds once, walked against the
// reference, between two runs walked with it. Only nodes a haplotype walks once count (an inverted repeat at
// the breakpoints is walked both ways). A run under `min_bp` between two others is absorbed into its
// neighbours without moving their reference range; a short run at a haplotype's end is kept as a flank, since
// an inversion can lie near the end of the locus. Each haplotype is read in the reference's orientation: one
// whose outermost runs of at least `min_bp` (its first and last steps when it has fewer than two) lie in
// falling reference order is stored on the other strand and is read reversed, so both strands give the same
// result. A reversed run of at least `min_bp`, in falling reference order (eight step pairs in ten), between
// two forward runs is an inversion; haplotypes whose reversed runs overlap by nine tenths are one, placed at
// the flanks most of them share. Positions are 0-based offsets along the reference path.
struct PathInversion {
    std::size_t left_end = 0;      // last base of the forward run before the inversion
    std::size_t inv_lo = 0;        // first base of the reversed run
    std::size_t inv_hi = 0;        // last base of the reversed run
    std::size_t right_start = 0;   // first base of the forward run after it
    std::string left_node, right_node;   // the reference nodes holding left_end and right_start
    std::vector<std::size_t> carriers;   // path indices
};

struct PathInversionScan {
    std::vector<PathInversion> inversions;
    // Per path, the reference offsets its once-walked reference nodes reach ({SIZE_MAX, 0} when none);
    // the reference path itself spans everything.
    std::vector<std::pair<std::size_t, std::size_t>> span;
};

// The shortest inversion and flank read when the caller's size floor is off (0).
constexpr std::size_t kPathInversionMinBp = 50;

PathInversionScan find_path_inversions(const Graph& graph, const PathRecord& reference, std::size_t min_bp);

} // namespace panvar
