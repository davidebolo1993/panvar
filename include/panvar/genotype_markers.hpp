#pragma once

#include "panvar/bubble_alleles.hpp"
#include "panvar/bubbles.hpp"
#include "panvar/genotype_blocks.hpp"
#include "panvar/genotype_reads.hpp"
#include "panvar/gfa.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace panvar {

// Which rule decides what counts as a usable marker. `Panvar` keeps any marker whose multiplicity
// varies across the block's alleles. `PanGenie` reproduces the rule in PanGenie's
// src/uniquekmercomputer.cpp: a k-mer must occur EXACTLY ONCE in an allele (`entry.second == 1`) and
// on EXACTLY ONE allele (`local_count > 1` is skipped), and is stored as presence/absence
// (`bool kmer_on_allele`). Holding everything else fixed, this isolates the marker rule as the single
// variable -- a cleaner comparison than two tools that also differ in model, VCF and parameters.
// PanGenie's rule bundles two independent constraints; `Unique` isolates the first so their
// contributions can be told apart:
//   Panvar   - multiplicity varies across the block's alleles (permissive)
//   Unique   - carried by exactly ONE allele, any multiplicity
//   PanGenie - carried by exactly one allele AND occurring exactly once in it, presence/absence
enum class MarkerRule { Panvar, Unique, PanGenie };

struct MarkerOptions {
    MarkerRule rule = MarkerRule::Panvar;
    std::size_t kmer_size = 31;
    std::size_t syncmer_s = 0;              // 0 = auto (default_syncmer_s)
    std::size_t min_markers = 10;           // per-allele "green" threshold
    std::size_t max_multiplicity = 0;       // 0 = no cap; multiplicity IS the copy-number signal
    bool require_region_unique = true;
    // Pairwise allele separation uses a dense n^2 scratch matrix, which is the fastest form but costs
    // n_alleles^2 * 8 bytes per block (1.7 MB at today's largest block of 463 alleles, but 200 MB at
    // 5,000 and 800 MB at 10,000 -- times the worker count, so it OOMs rather than merely slowing).
    // Above this many alleles a sparse accumulator is used instead: identical results, O(n) memory,
    // roughly 40% slower. 0 disables the dense path entirely.
    std::size_t max_dense_alleles = 2048;
    // Last-resort path for panels where even the sparse accumulator is too slow: score each allele
    // only against its K most similar siblings, found by MinHash sketch rather than by touching every
    // marker carrier. Approximate -- the hardest sibling can be missed -- so 0 (exact) is the default
    // and any run using it should be compared against an exact run first. An exact bound,
    // differ(a,b) >= |carried[a]-carried[b]|, is applied alongside to recover cheap certain misses.
    std::size_t sep_top_k = 0;
    std::size_t sketch_size = 64;
    std::size_t threads = 0;
};

// Marker inventory for one allele of one bubble. Both candidate units are carried side by side so
// the audit can compare them: syncmer NODES (present in this allele, absent from every sibling) and
// syncmer ADJACENCIES (the same rule applied to consecutive-syncmer pairs, which is what carries a
// deletion junction or a shared-sequence rearrangement).
struct AlleleMarkers {
    std::size_t bubble_id = 0;
    std::size_t allele_id = 0;
    std::size_t n_haplotypes = 0;
    std::size_t allele_bp = 0;
    bool is_reference = false;

    std::size_t n_syncmers_total = 0;       // syncmer occurrences in the allele walk
    std::size_t n_syncmers_distinct = 0;

    // A marker separates two alleles when its MULTIPLICITY differs between them. Presence/absence is
    // the special case (0 vs >0); equal presence at different copy number is the tandem case, which a
    // strict private-set rule misses entirely.
    std::size_t n_informative_nodes = 0;    // non-constant multiplicity across this bubble's alleles
    std::size_t n_informative_edges = 0;
    std::size_t n_carried_nodes = 0;        // informative markers this allele actually carries
    std::size_t n_carried_edges = 0;

    // The binding constraint: separation from the single hardest sibling to tell apart. An allele is
    // genotypable only if this is large enough, however many markers it has in total.
    std::size_t min_separating_nodes = 0;
    std::size_t min_separating_edges = 0;
    std::size_t hardest_sibling = 0;
    double nearest_sibling_jaccard = 0.0;   // node-set Jaccard against the most similar sibling

    std::size_t n_nodes_lost_region = 0;    // dropped because they occur elsewhere in the region
    std::size_t n_edges_lost_region = 0;

    std::vector<std::uint64_t> node_codes;  // final marker sets
    std::vector<std::uint64_t> edge_keys;
};

struct BubbleMarkerReport {
    std::size_t bubble_id = 0;
    std::string source;
    std::string sink;
    std::size_t n_inside = 0;
    std::size_t n_alleles = 0;        // distinct SEQUENCES (the genotyping unit)
    std::size_t n_alleles_walk = 0;   // distinct walk signatures (what `call` groups by)
    std::size_t n_traversing = 0;
    std::size_t singleton_alleles = 0;      // alleles carried by exactly one haplotype
    double singleton_hap_frac = 0.0;        // fraction of haplotypes on a singleton allele
    bool folded = false;                    // carries a REP self-loop node
    bool has_reference = false;
};

struct MarkerPanel {
    std::vector<BubbleMarkerReport> bubbles;
    std::vector<AlleleMarkers> markers;     // flattened over (bubble, allele)
    std::size_t kmer_size = 31;
    std::size_t syncmer_s = 0;
};

// Per-block marker sufficiency: the panel can *represent* a block's alleles (Phase A1), but that is a
// different question from whether short reads can *tell them apart*. This measures the second.
struct BlockMarkerStats {
    std::size_t block_index = 0;
    bool is_bubble = true;
    std::size_t bubble_id = 0;
    std::size_t n_alleles = 0;
    std::size_t n_haplotypes = 0;
    std::size_t block_bp = 0;               // median allele length

    std::size_t n_informative_nodes = 0;
    std::size_t n_informative_edges = 0;
    std::size_t median_sep_nodes = 0;       // median over alleles of separation from the hardest sibling
    std::size_t median_sep_edges = 0;
    std::size_t unseparable_alleles = 0;    // alleles with zero separation from some sibling

    // Haplotype-weighted: the share of haplotypes sitting on an allele that reads could resolve.
    double separable_mass_nodes = 0.0;
    double separable_mass_edges = 0.0;

    // A3: the multiplicity cap exists to drop repeat-derived markers, but inside a tandem array
    // multiplicity IS the copy-number signal, so this reports what the cap would remove.
    std::size_t markers_over_cap = 0;
    std::size_t max_marker_multiplicity_seen = 0;

    double seconds = 0.0;                   // A3: where the time actually goes
};

// Score marker sufficiency over a whole block chain. Alleles come from `BlockAlleles::allele_seq`,
// so this handles bubble and backbone blocks identically.
// `out_panel`, when non-null, is filled with the marker universe reads are counted against: dense
// slots, each allele's expected multiplicities, and each block's invariant depth anchors.
std::vector<BlockMarkerStats> build_block_marker_panel(
    const std::vector<Block>& chain,
    const std::vector<BlockAlleles>& blocks,
    const MarkerOptions& options,
    ReadPanel* out_panel = nullptr,
    const Graph* graph_for_region_uniqueness = nullptr,
    bool want_separation_stats = true);

void write_block_marker_audit(
    const std::string& out_prefix,
    const std::vector<BlockMarkerStats>& stats);

// Build the per-bubble syncmer graphs and score both marker units on them.
MarkerPanel build_marker_panel(
    const Graph& graph,
    const std::vector<Bubble>& bubbles,
    const std::string& reference_path,
    const MarkerOptions& options);

// Audit tables: <out_prefix>.audit.bubbles.tsv and <out_prefix>.audit.alleles.tsv.
void write_marker_audit(const std::string& out_prefix, const MarkerPanel& panel);

} // namespace panvar
