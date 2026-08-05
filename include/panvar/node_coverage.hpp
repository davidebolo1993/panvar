#pragma once

// A second evidence path for `genotype`, alongside the syncmer marker counts: per-node coverage, as
// cosigt derives it from gfainject + gafpack.
//
// Why a second path rather than a replacement. The marker path is validated everywhere except tandem
// arrays, where it calls the array up to two repeat units too long. The cause is that a marker the
// candidate allele does not carry is predicted at the error background, so real counts there veto the
// candidate -- and under leave-one-out every candidate lacks some of the sample's sequence, while a
// longer array carries more distinct unit variants and so lacks less. Nodes do not have that problem:
// the graph holds an array as a cycle, so alleles differ in how many TIMES they cross a node rather
// than in which nodes they contain. At lpa's KIV-2 block the haplotypes use 2319-2347 of 3791 nodes
// (a 1% spread) while their step counts vary threefold.
//
// Why alignment rather than syncmers here. The median node in that graph is 1 bp and 86% of nodes are
// shorter than k=31, so no k-mer can be assigned to them at all -- and those short nodes are exactly
// the SNP and indel bubbles that separate haplotypes. An alignment carries base-level coverage across
// node boundaries and reaches them.

#include "panvar/gfa.hpp"
#include "panvar/genotype_blocks.hpp"
#include "panvar/bubble_path.hpp"
#include "panvar/bubbles.hpp"

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace panvar {

// Dense node numbering shared by the panel and the sample, so the two vectors are directly comparable
// and neither needs a hash lookup per access.
struct NodeIndex {
    std::vector<std::string> id;                            // dense index -> GFA node name
    std::vector<std::uint32_t> length;                      // dense index -> node sequence length
    std::unordered_map<std::string, std::uint32_t> of;      // GFA node name -> dense index
    std::size_t size() const { return id.size(); }
};

NodeIndex build_node_index(const Graph& graph);

// How many times each panel path traverses each node. This is `odgi paths -H`, and at a cyclic array
// it is the copy number of each repeat-unit variant.
struct PanelCoverage {
    std::vector<std::string> path_names;
    std::vector<std::vector<std::uint32_t>> by_path;        // [path][node]
    std::vector<std::string> path_seq;                      // spelled, for the aligner index
};

PanelCoverage build_panel_coverage(const Graph& graph, const NodeIndex& index);

struct CoverageOptions {
    std::string preset = "sr";
    // Alternative placements kept per read. A read from inside a tandem array matches every copy and
    // every haplotype carrying it, so this has to be generous or the array is silently under-covered.
    std::size_t best_n = 64;
    // Divide each node's accumulated bases by its length, so coverage is in read-depth units and a
    // long node does not outweigh a short one purely by being long.
    bool len_scale = true;
    // How a read that aligns in several places is charged to nodes.
    //
    //   UnionMax  -- take the largest overlap each NODE receives across the read's alignments, and add
    //                that once. The panel holds hundreds of near-identical haplotypes, so a read from a
    //                unique region aligns to all of them at the SAME graph location: those are one
    //                placement seen many times, not competing placements. Likewise inside a cyclic
    //                array, where every copy is the same node. This keeps coverage in read-depth units.
    //   UniformSplit -- divide by the number of alignments, which is cosigt's `--weight-queries`. It is
    //                scale-invariant nonsense for a likelihood but harmless for cosine, which is
    //                itself scale-invariant. Kept as the baseline to compare against.
    //   All       -- charge every alignment in full. Only useful for showing what the other two fix.
    //
    // Measured on the synthetic locus with 17 panel paths: UniformSplit gives self-correlation 0.767
    // and a depth slope of 9.4 where the truth is 15, because a node crossed by one path keeps 1/17 of
    // the read while a node crossed by all seventeen keeps all of it.
    enum class Placement { UnionMax, UniformSplit, All };
    Placement placement = Placement::UnionMax;
    std::size_t threads = 0;
};

struct SampleCoverage {
    std::vector<double> node;
    std::uint64_t reads = 0;          // reads seen
    std::uint64_t aligned = 0;        // reads with at least one placement
    std::uint64_t placements = 0;     // alignments kept across all reads
    double bases_placed = 0.0;        // sum of aligned target bases, before any weighting
};

// Map reads to the panel's path sequences and project each alignment onto the nodes that path
// traverses over the aligned interval.
SampleCoverage inject_reads(
    const Graph& graph,
    const NodeIndex& index,
    const PanelCoverage& panel,
    const std::vector<std::string>& read_paths,
    const CoverageOptions& options);

// Gate A. Each is exact and each catches a different way the projection can be wrong; a correlation
// that is merely high is not a pass, because a systematic offset or a flipped orientation still
// correlates well.
struct CoverageAudit {
    // A read taken from a known path at a known offset must land on exactly the nodes that path
    // traverses there -- catches off-by-one and orientation errors in the projection.
    std::size_t probe_total = 0;
    std::size_t probe_exact = 0;
    // sum(coverage * node length) against total placed bases -- catches double counting and loss.
    double mass_ratio = 0.0;
    // A haplotype's own reads against its own traversal vector -- the end-to-end check, since it
    // exercises alignment, projection, weighting and node indexing together.
    double self_pearson = 0.0;
    double self_slope = 0.0;          // regression slope; should be the per-haplotype read depth
};

CoverageAudit audit_coverage(
    const Graph& graph,
    const NodeIndex& index,
    const PanelCoverage& panel,
    const SampleCoverage& sample,
    const std::vector<std::string>& sample_paths,
    const CoverageOptions& options);

// Per-allele node multiplicity inside one block, using the same step extraction the marker panel
// uses, so an allele means the same thing on both evidence paths and the two are comparable.
std::vector<std::vector<std::uint32_t>> block_allele_node_vectors(
    const Graph& graph,
    const std::vector<BubblePathIndex>& path_indexes,
    const std::vector<Bubble>& bubbles,
    const Block& block,
    const BlockAlleles& alleles,
    const NodeIndex& index);

// Score every allele pair at one block by the coverage the reads put on its nodes. Negative binomial
// on node coverage, the direct analogue of the marker emission -- the difference being what the
// evidence is, not how it is weighed.
struct CoverageScore {
    std::size_t allele1 = 0;
    std::size_t allele2 = 0;
    double loglik = 0.0;
    double cosine = 0.0;
    double pearson = 0.0;
    std::size_t bp = 0;
};

std::vector<CoverageScore> score_block_by_coverage(
    const std::vector<std::vector<std::uint32_t>>& allele_vec,
    const std::vector<std::size_t>& allele_bp,
    const SampleCoverage& sample,
    const NodeIndex& index,
    double lambda,
    double overdispersion);

void write_node_coverage(
    const std::string& out_prefix,
    const NodeIndex& index,
    const PanelCoverage& panel,
    const SampleCoverage& sample);

} // namespace panvar
