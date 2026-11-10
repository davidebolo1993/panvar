#pragma once

// Banded fitting alignment (unit-cost edit distance), backed by edlib's bit-parallel
// (Myers) kernel. Used by panphorte's approximate tandem-copy collapse and by call's
// cross-haplotype event-merge sequence identity. See src/align.cpp for the edlib adapter.

#include <cstddef>
#include <cstdint>
#include <string>

namespace panvar {

struct FitAlignResult {
    bool ok = false;            // false if the best fit needs more than `band` edits
    std::size_t edits = 0;      // edit distance (mismatches + indels)
    std::size_t query_len = 0;  // |query| (all of query is consumed)
    std::size_t target_consumed = 0; // target bases aligned (copy length in target)
    double identity = 0.0;      // 1 - edits / max(query_len, target_consumed)
};

// Fit the entire `query` into `target` starting at target[0] (no leading target
// gap), with a free end gap on the target, using at most `band` edits (edlib SHW mode).
// Returns the best end position on the target (the copy length) and identity; ok is false
// when more than `band` edits are needed.
FitAlignResult fit_align(const std::string& query, const std::string& target, std::size_t band);

// Global (Needleman-Wunsch) alignment of `target` and `query` as edit operations, one per alignment
// column: 'M' match, 'X' mismatch, 'I' a base of `query` only, 'D' a base of `target` only. Empty
// sequences give the obvious run. Returns the edit distance. Used by benchmark to cut the residual
// between a reconstruction and a haplotype into events.
std::size_t nw_edit_ops(const std::string& target, const std::string& query, std::string& ops);

// The fewest edits placing all of `query` somewhere inside `target` (edlib HW mode): SIZE_MAX when more
// than `max_k`.
std::size_t infix_distance(const std::string& query, const std::string& target, std::size_t max_k);

// The best placement of all of `query` inside `target` (edlib HW mode) within `max_k` edits: where it lies in
// `target` ([start, end], 0-based) and its gap-compressed identity, each insertion or deletion counting once
// whatever its length (1 - (mismatches + gaps) / aligned columns, gaps counted as runs). -1 when none.
double infix_identity(const std::string& query, const std::string& target, std::size_t max_k, std::size_t& start,
                      std::size_t& end);

// As nw_edit_ops, but only when the edit distance is at most `max_k`: SIZE_MAX (and `ops` empty)
// otherwise. edlib bands the search to max_k, so a small bound keeps a long alignment cheap.
std::size_t nw_edit_ops_bounded(const std::string& target, const std::string& query, std::size_t max_k,
                                std::string& ops);

// Global edit distance only (EDLIB_TASK_DISTANCE): no traceback, so no alignment length. About 2x
// faster than nw_edit_ops on long sequences; for callers that need only the count. With max_k,
// edlib stops once the distance must exceed it, and any distance over max_k is returned as max_k + 1.
std::size_t nw_distance(const std::string& a, const std::string& b, std::size_t max_k = SIZE_MAX);

} // namespace panvar
