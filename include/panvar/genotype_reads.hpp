#pragma once

#include "panvar/genotype_blocks.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace panvar {

// The markers one allele carries, as (global slot, expected multiplicity). Sparse on purpose: a
// dense allele-by-marker matrix would be tens of MB per block on the larger loci, while the sparse
// form costs only the total number of marker occurrences.
struct AlleleMarkerSet {
    std::vector<std::pair<std::uint32_t, std::uint32_t>> nodes;
    std::vector<std::pair<std::uint32_t, std::uint32_t>> edges;
};

// Everything reads are counted against: the marker universe with dense slots, each allele's expected
// multiplicities, and the invariant markers that serve as depth anchors.
struct ReadPanel {
    std::size_t kmer_size = 31;
    std::size_t syncmer_s = 0;
    bool all_kmers = false;
    std::vector<std::uint64_t> node_codes;                  // slot -> canonical syncmer code
    std::vector<std::uint64_t> edge_keys;                   // slot -> adjacency key
    // Every adjacency any panel allele carries, BEFORE confinement and the marker rules thin it down.
    // `edge_keys` is only the retained informative subset, so judging a read's adjacencies against it
    // counts constant and filtered panel adjacencies as novel -- which is most of them. Only this set
    // can say whether an arrangement is genuinely off-panel.
    std::vector<std::uint64_t> all_edge_keys;
    std::vector<std::vector<AlleleMarkerSet>> by_block;     // [block][allele]
    // Markers carried by every allele of a block at the same multiplicity. They cannot separate
    // anything, which is exactly why they are the right depth reference: their observed count is a
    // clean read-out of coverage, uncontaminated by which allele the sample carries.
    std::vector<std::vector<std::uint32_t>> anchor_slots;   // [block] -> node slots
    // Per-slot confinement audit, filled only when requested. `vary` is the number of blocks whose
    // count depends on the genotype, `occ` the number of blocks the marker appears in at all, `actual`
    // its panel-wide occurrence count from re-spelling every path, and `expected` what the blocks that
    // retained it account for. A marker survives only when vary <= 1 and actual <= expected, so these
    // four numbers say exactly why any given marker is still here.
    std::vector<std::uint32_t> dbg_vary;
    std::vector<std::uint32_t> dbg_occ;
    std::vector<std::uint64_t> dbg_actual;
    std::vector<std::uint64_t> dbg_expected;
    // Per block: how many distinct fragment-length windows the surviving markers occupy. This, not the
    // marker count, is how many independent observations the block really supplies -- markers inside one
    // fragment are carried by the same reads and rise and fall together.
    std::vector<double> marker_clumps;
    std::size_t blocks_restored = 0;          // blocks whose markers were put back for allele balance
    std::size_t region_filtered_markers = 0;   // dropped for occurring elsewhere in the region
    std::size_t dropped_multi_block = 0;       // ...because they appear in more than one block
    std::size_t dropped_over_expected = 0;     // ...because the panel shows more copies than blocks own
    std::size_t informative_before_filter = 0;
    std::size_t dropped_adjacent_blocks = 0;   // multi-block markers confined to neighbouring blocks
    std::size_t dropped_distant_blocks = 0;    // ...spread across non-adjacent blocks
    // Does longer context buy confinement? A single syncmer's count mixes every block it occurs in,
    // but the observed PAIR (a -> b) may occur in only one of them. Comparing these two rates says
    // whether read-level context can reach blocks where no single marker is usable.
    std::size_t vary_nodes = 0;                // syncmers whose count varies with genotype somewhere
    std::size_t confined_vary_nodes = 0;       // ...and varies in exactly one block
    std::size_t vary_edges = 0;                // same, for 2-syncmer context
    std::size_t confined_vary_edges = 0;
    // Per block, the up-to-3 other blocks its informative markers most often also occur in.
    std::vector<std::vector<std::pair<std::uint32_t, std::uint32_t>>> block_overlap;
};

struct ReadCounts {
    std::vector<std::uint32_t> node;        // parallel to ReadPanel::node_codes
    std::vector<std::uint32_t> edge;        // parallel to ReadPanel::edge_keys
    std::uint64_t reads = 0;
    std::uint64_t bases = 0;
    std::uint64_t syncmers = 0;
    std::uint64_t matched_syncmers = 0;     // syncmers landing on a known marker
    std::uint64_t novel_adjacencies = 0;    // consecutive pair whose adjacency is in no allele:
                                            // positive evidence of sequence the panel lacks
};

// Project reads onto the panel exactly as syng's -noAddK does: take each read's closed syncmers in
// order, map them to known markers, and treat anything unmatched as absent rather than as new.
ReadCounts count_reads(
    const std::vector<std::string>& read_paths,
    const ReadPanel& panel,
    std::size_t threads);

struct BlockDepth {
    std::size_t block_index = 0;
    std::size_t n_anchor = 0;
    // The raw anchor-count median, never rewritten by a depth model. `median` below is the model's
    // fitted value; this is the observation, which is what a second pass has to reason from.
    double anchor_median = 0.0;
    double median = 0.0;        // anchor count median
    double mad = 0.0;
    double lambda_hap = 0.0;    // expected count for one haplotype copy (median / 2, diploid)
    bool usable = false;
    bool uneven = false;        // MAD/median above the tolerance: coverage too ragged to trust
};

// How the per-haplotype depth lambda is estimated. All of them must answer the same question -- what
// read count does ONE copy of a marker produce -- and they differ in what they assume about the sample.
//
//   Median   per-block anchor median / 2, shrunk toward the region. Assumes BOTH of the sample's
//            haplotypes traverse every block. Where one bypasses (a deletion spanning the site) only
//            one copy is present, the anchors read lambda rather than 2*lambda, and the estimate halves
//            -- after which the model needs two copies to explain the data and calls homozygous.
//   Quantile a single region-wide lambda from a high quantile of the per-block anchor medians, so a
//            deletion covering most of the locus cannot drag it. Robust as long as some reasonable
//            fraction of blocks is diploid; the quantile is where that assumption is stated.
//   Bases    lambda from total read bases over the reference length, independent of block structure
//            entirely. Biased when the sample's own length differs from the reference, which is exactly
//            what a large indel does.
//   Joint    two passes: call once with Median, then re-estimate lambda using how many of each
//            block's CALLED alleles actually traverse it, and call again. This is the only one that
//            can be right when a haplotype bypasses many blocks, because then "one haplotype at full
//            depth" and "two at half depth" produce identical anchor counts and lambda is simply not
//            identifiable from the anchors alone -- it needs the genotype, which needs lambda.
enum class DepthModel { Median, Quantile, Bases, Joint };

// Per-block depth from that block's own invariant markers, with a region-wide fallback for blocks
// that have too few of their own.
std::vector<BlockDepth> estimate_depth(
    const ReadPanel& panel,
    const ReadCounts& counts,
    std::size_t min_anchors,
    double uneven_tolerance,
    DepthModel model = DepthModel::Median,
    double depth_quantile = 0.75,
    std::size_t region_bp = 0);

void write_read_audit(
    const std::string& out_prefix,
    const std::vector<Block>& chain,
    const ReadPanel& panel,
    const ReadCounts& counts,
    const std::vector<BlockDepth>& depth);

} // namespace panvar
