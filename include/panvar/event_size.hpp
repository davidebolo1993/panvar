#pragma once

// How call and benchmark size a change between two sequences: one rule, so an event call leaves under
// --min-sv-bp is one benchmark counts as small, and one call reports is one benchmark counts as large.

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace panvar {

// A banded alignment of two allele sequences above this many cells is not made here (one byte of traceback
// per cell, so 32 MB per aligning thread).
constexpr std::size_t kSeqAlignCellCap = 32000000;
// A replaced run is sized by aligning its two sides up to this larger side; above it, by the larger side.
constexpr std::size_t kSubEditMaxBp = 20000;

// One change between two alleles found by aligning them: [t0, t1) of the reference allele replaced by
// [q0, q1) of the haplotype's. An empty side makes it a deletion or an insertion.
struct AlignedChange {
    std::size_t t0 = 0, t1 = 0, q0 = 0, q1 = 0;
};

// The changes an alignment of hap to ref makes, from its CIGAR ((op, length) pairs, op 0 aligned columns
// compared base by base, 7 '=', 8 'X', 1 hap bases, 2 ref bases). The changes are the gaps. Two gaps are
// one change, a run replaced, unless the sequence between them aligns well enough to pay for opening a
// new long gap: some stretch of it scoring 24 (the long-gap open cost) at 2 per match and -4 per
// mismatch, as twelve matching bases do. False if the CIGAR does not cover both.
bool changes_from_cigar(const std::vector<std::pair<int, int>>& cigar, const std::string& ref,
                        const std::string& hap, std::vector<AlignedChange>& out);

// Global alignment of hap to ref, banded to w diagonals either side (widened where the band could hide a
// better placement), with minimap2's default scoring: match 2, mismatch 4, and two gap costs, 4 + 2 per base
// and 24 + 1 per base, the cheaper applying. Its changes as changes_from_cigar cuts them. False when the
// alignment failed.
bool aligned_changes(const std::string& ref, const std::string& hap, int w, std::vector<AlignedChange>& out);

// The size of a change from `ref` to `hap`: their shared ends trimmed, the rest aligned (aligned_changes),
// and the largest single change, each counting its larger side. A one-sided change is its length. Above
// kSubEditMaxBp a side, or kSeqAlignCellCap cells, the larger side whole. Sized by the larger side alone, one
// base changed on a long run would be a long event; by edit distance, an unrelated replacement would be about
// half its length; summed over the changes, scattered small edits in a repeat would add up to one large one.
std::size_t change_size_bp(const std::string& ref, const std::string& hap);

} // namespace panvar
