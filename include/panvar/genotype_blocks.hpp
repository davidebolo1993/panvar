#pragma once

#include "panvar/bubble_path.hpp"
#include "panvar/bubbles.hpp"
#include "panvar/gfa.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace panvar {

// A locus is genotyped as an alternating chain along the reference: each called bubble is a block,
// and the stretch between two consecutive bubbles is also a block. The backbone blocks are not
// filler -- they carry the sub-threshold variation `--min-variant-bp` kept out of the bubble list,
// which is 60-83% of each haplotype and is what actually identifies which haplotypes a sample is on.
enum class BlockKind { Bubble, Backbone, Flank };

struct Block {
    std::size_t index = 0;        // position along the chain
    BlockKind kind = BlockKind::Bubble;
    std::size_t bubble_id = 0;    // 0 for a backbone block
    std::string source;
    std::string sink;
};

// Alleles of one block, grouped by the sequence they spell: two haplotypes taking different graph
// routes through identical bases are one allele, because no read-based method can separate them.
struct BlockAlleles {
    std::size_t block_index = 0;
    std::size_t n_alleles = 0;
    std::size_t n_walk_alleles = 0;                   // before sequence grouping
    std::vector<std::size_t> allele_haplotypes;       // per allele: carrier count
    std::vector<std::size_t> allele_bp;               // per allele: spelled length
    std::unordered_map<std::string, std::size_t> allele_of;  // haplotype name -> allele index
    std::vector<std::string> allele_seq;                     // representative sequence per allele
};

// For alleles carried by exactly one haplotype -- the ones a panel-based model cannot supply -- how
// much of their content is nonetheless present elsewhere in the same block? An allele whose syncmers
// and adjacencies all occur in some sibling is a novel *arrangement* of known material, recoverable
// by matching at finer granularity. One carrying syncmers seen nowhere is novel *sequence*, which
// only local assembly can reach. The split decides whether "guess the unrepresented structure" is
// feasible at all.
struct NoveltyReport {
    std::size_t private_alleles = 0;
    double mean_syncmer_reuse = 0.0;      // fraction of a private allele's syncmers seen in a sibling
    double mean_adjacency_reuse = 0.0;    // same for consecutive-syncmer adjacencies
    std::size_t fully_reusable = 0;       // private alleles that are pure rearrangements
    std::size_t with_novel_sequence = 0;  // private alleles carrying syncmers found nowhere else
};

NoveltyReport measure_novelty(
    const std::vector<Block>& chain,
    const std::vector<BlockAlleles>& blocks,
    std::size_t kmer_size,
    std::size_t syncmer_s);

// Bubbles must arrive in reference order (bubble id order, which `graph_sort` guarantees by
// renumbering nodes into reference order). Emits bubble blocks interleaved with the backbone
// stretches between them.
std::vector<Block> build_block_chain(const std::vector<Bubble>& bubbles);

BlockAlleles enumerate_block_alleles(
    const Graph& graph,
    const std::vector<BubblePathIndex>& path_indexes,
    const std::vector<Bubble>& bubbles,
    const Block& block,
    std::size_t threads);

struct LinkageReport {
    std::size_t n_haplotypes = 0;
    std::size_t n_blocks = 0;
    std::size_t uniquely_identified = 0;
    // Mean size of the compatible-haplotype set after consuming the first N blocks of the chain.
    std::vector<double> collapse_curve;
    std::vector<std::size_t> class_size;   // per haplotype, final compatible-set size

    // Leave-one-out ceiling. Identifying a haplotype that is *in* the panel is easy and not the real
    // task: a real sample is not in the panel. For each haplotype in turn, drop it, take the
    // remaining haplotype that agrees with it at the most blocks, and ask whether that stand-in
    // carries the same allele at each BUBBLE block. That is the best any panel-based method could do
    // with perfect read evidence, so it bounds the whole approach.
    std::vector<double> loo_bubble_agreement;   // per haplotype, fraction of bubble blocks recovered
    double loo_mean_agreement = 0.0;
    std::size_t loo_perfect = 0;                // haplotypes whose bubble alleles are fully recovered

    // The mosaic ceiling, which is the one that actually bounds an HMM: a recombination-aware model
    // may use a different panel haplotype at every block, so what limits it is not whether one
    // neighbour matches everywhere but whether SOME remaining haplotype carries the held-out
    // haplotype's allele at each bubble. Anything below this is a modelling loss, not a panel limit.
    double mosaic_ceiling = 0.0;                // mean over haplotypes of recoverable bubble alleles
    double singleton_mass = 0.0;                // fraction of (hap, bubble) cells with a private allele
};

// How far the block chain narrows "which panel haplotype is this?". A sample that is a panel
// haplotype can only be pinned down as tightly as this allows, so it bounds the whole model.
LinkageReport measure_linkage(
    const Graph& graph,
    const std::vector<BlockAlleles>& blocks);

// Fills the leave-one-out fields of `rep`; call after measure_linkage.
void measure_leave_one_out(
    const Graph& graph,
    const std::vector<Block>& chain,
    const std::vector<BlockAlleles>& blocks,
    LinkageReport& rep);

void write_linkage_audit(
    const std::string& out_prefix,
    const Graph& graph,
    const std::vector<Block>& chain,
    const std::vector<BlockAlleles>& blocks,
    const LinkageReport& report);

} // namespace panvar
