#pragma once

// Banded fitting alignment (unit-cost edit distance), backed by edlib's bit-parallel
// (Myers) kernel. Used by panphorte's approximate tandem-copy collapse and by call's
// cross-haplotype event-merge sequence identity. See src/align.cpp for the edlib adapter.

#include <cstddef>
#include <string>

namespace panvar {

struct FitAlignResult {
    bool ok = false;            // false if the optimal path left the band
    std::size_t edits = 0;      // edit distance (mismatches + indels)
    std::size_t query_len = 0;  // |query| (all of query is consumed)
    std::size_t target_consumed = 0; // target bases aligned (copy length in target)
    double identity = 0.0;      // 1 - edits / max(query_len, target_consumed)
};

// Fit the entire `query` into `target` starting at target[0] (no leading target
// gap), with a free end gap on the target, within +/- `band` of the diagonal.
// Returns the best end position on the target (the copy length) and identity.
// Linear space (two DP rows).
FitAlignResult fit_align(const std::string& query, const std::string& target, std::size_t band);

struct NwAlign {
    std::size_t edits = 0;    // global edit distance (delta)
    std::size_t aln_len = 0;  // total alignment length in columns (S)
};

// Global (Needleman-Wunsch) edit distance with the alignment length, via edlib
// EDLIB_MODE_NW + EDLIB_TASK_PATH. Used by benchmark's round-trip QV:
// QV = -10*log10(max(0.5, edits) / aln_len).
NwAlign nw_edit_distance(const std::string& a, const std::string& b);

} // namespace panvar
