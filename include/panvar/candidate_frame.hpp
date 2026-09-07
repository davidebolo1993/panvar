#pragma once

#include "panvar/genotype_blocks.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace panvar {

// ---------------------------------------------------------------------------------------------
// THE CANDIDATE FRAME: a panel haplotype's walk bytes together with its VERIFIED block map.
//
// Extracted into its own header so the genotype command can use it without pulling in
// genotype_fragments.hpp, which declares a second, unrelated BlockCall and therefore cannot coexist
// with genotype.hpp. Splitting the shared type is the right fix; making one of the callers work
// around the collision would not be.
struct CandidateFrame {
    std::string seq;                    // walk bytes, authoritative
    std::vector<std::size_t> offsets;   // segment start offsets IN WALK COORDINATES, ascending
    // The BLOCK INDEX of each segment, in walk order. For an antiparallel candidate the blocks run
    // backwards along the walk, so the k-th segment is block (nblocks-1-k) -- and a coordinate
    // lookup that returned the segment index would report mirrored block numbers. Measured: an
    // antiparallel duplicate of an existing path changed a fragment's scope from {1} to {1,2},
    // which is impossible for a strand-symmetric scorer and was exactly this.
    std::vector<std::uint32_t> block_at;
    bool reverse_frame = false;
    bool ok = false;                    // false: no verified coordinate map for this candidate
    // A PATH MAY END INSIDE A BLOCK. Its assembly contig, or the interval that cut it, can stop
    // partway through a block, so the concatenated alleles are a correct PREFIX (or suffix) of the
    // walk rather than the whole of it. That is not a spelling error -- the bytes agree -- but the
    // "concatenation equals the walk" test cannot certify a map over the uncovered remainder.
    // Refusing the candidate outright loses it from the panel; inventing a block for the remainder
    // would be worse. So the map covers [mapped_lo, mapped_hi) in WALK coordinates and the rest is
    // reported as unmapped, where a position has no block rather than a guessed one.
    // Measured: cyp2d6 NA18989#1#haplotype1 ends 1978 bp past bubble 8's sink, short of bubble 9.
    bool partial = false;
    std::size_t mapped_lo = 0, mapped_hi = 0;
};

// Returned by ordered_block_span for a position outside [mapped_lo, mapped_hi). Such a position
// belongs to no block, so its mass is unattributable: it can never be dropped by restricting the
// scope, and it can never justify adding a block to one.
inline constexpr std::uint32_t kUnmappedBlock = 0xFFFFFFFFu;


// The block index range a walk interval touches, in CHAIN order. Returns kUnmappedBlock at either
// end for a position outside the frame's verified window.
std::pair<std::uint32_t, std::uint32_t> ordered_block_span(
    const CandidateFrame& frame, long start, long end);

// Build a candidate's frame by checking the concatenated block alleles against the walk bytes.
CandidateFrame build_candidate_frame(
    const std::vector<BlockAlleles>& blocks,
    const std::string& name,
    const std::string& walk);

}  // namespace panvar
