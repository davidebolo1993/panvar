#pragma once

#include "panvar/bubbles.hpp"
#include "panvar/genotype_model.hpp"
#include "panvar/gfa.hpp"
#include "panvar/mosaic_spelling.hpp"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace panvar {

// The graph-facing adapter. It is the only place where graph structure becomes statistical input:
// every later stage reads `StructuralLocusInput` and never touches nodes, edges or walks. The
// adapter therefore owns exactly one decision — where the block boundaries are — and must leave the
// panel sequences byte-for-byte intact, which `verify_complete_panel_round_trip` then proves.
struct LocusInputOptions {
    std::string locus_id;
    // Reference path name, already resolved against the graph. It fixes the chain's coordinate frame
    // and therefore which templates count as antiparallel; it is not otherwise privileged.
    std::string reference_path_name;
};

struct LocusBlockReport {
    std::string block_id;
    GenotypeBlockKind kind = GenotypeBlockKind::Variant;
    std::size_t bubble_id = 0;             // 0 for a context block
    std::size_t source_alleles = 0;        // distinct walks, before sequence collapse
    std::size_t canonical_alleles = 0;     // distinct sequences, after collapse
    std::size_t min_bp = 0;
    std::size_t max_bp = 0;
    std::size_t bypass_templates = 0;      // templates that do not walk this block at all
    // Bytes the divergence trim moved out of this variant block into its left and right neighbours.
    std::size_t trimmed_prefix_bp = 0;
    std::size_t trimmed_suffix_bp = 0;
    // True when this block's boundary is defined by where the graph cut a node rather than by where
    // alleles diverge, so a node split would change the model. See require_representation_independent.
    bool representation_dependent = false;
};

// How one template covers one block. This replaces a boolean `traverses`, which conflated two states
// that are not alike: a template that walks PAST a block still has the block's neighbours on either
// side and has genuinely deleted the region, while a template whose path ENDS before the block never
// observed it at all. Both spell zero bases, so the second silently drops the bases the path does
// carry there, and the blocks stop tiling the haplotype.
enum class BlockCoverage {
    Traversing,      // walks the block and spells its allele
    EmptyInterior,   // walks both boundaries with nothing strictly between them
    BridgedBypass,   // does not walk it, but is observed on BOTH sides: an admissible deletion
    Uncovered,       // does not walk it and is not observed on one side: absent coverage, not a deletion
};

struct TemplateBlockCoverage {
    std::string template_name;
    std::string block_id;
    BlockCoverage coverage = BlockCoverage::Traversing;
    // Meaningful for Uncovered only: the side the template has no observation on.
    bool uncovered_on_left = false;
};

// "uncovered_left" / "uncovered_right" for Uncovered, the classification name otherwise.
const char* block_coverage_label(BlockCoverage coverage, bool uncovered_on_left);

struct LocusInputReport {
    std::size_t variant_blocks = 0;
    std::size_t context_blocks = 0;
    std::size_t total_source_alleles = 0;
    std::size_t total_canonical_alleles = 0;
    std::size_t forward_templates = 0;
    std::size_t reverse_templates = 0;
    // Templates that do not span the chain in a consistent direction. They are reported as partial
    // frames and refused by `prepare_genotype_model`, never silently completed.
    std::vector<std::string> refused_templates;
    std::size_t contained_bubbles_dropped = 0;
    std::size_t abutting_boundaries_trimmed = 0;
    // Bases moved out of variant blocks into their context neighbours by the divergence trim. This
    // is what makes the decomposition independent of where the graph cut its boundary nodes.
    std::size_t boundary_bases_moved_to_context = 0;
    // Blocks whose boundaries depend on graph representation, in chain order. A variant block that
    // carries a bypass allele cannot be divergence-trimmed -- its traversing alleles share boundary
    // bytes the bypassing template does not have -- so where the graph cut its boundary nodes decides
    // how many bytes sit in it and in both neighbours.
    std::vector<std::string> representation_dependent_blocks;
    std::vector<LocusBlockReport> blocks;
    // Every (template, block) pair the template does not simply traverse, in chain order. Only the
    // non-traversing pairs are kept: they are the interesting ones and the full cross product is
    // templates x blocks. Reporting only -- no allele, refusal or fingerprint depends on it.
    std::vector<TemplateBlockCoverage> non_traversing_coverage;
};

struct LocusInputResult {
    StructuralLocusInput input;
    // The bytes each panel path spells in its own frame, for the round-trip verifier only. These are
    // deliberately kept outside the prepared model so no authoritative sequence can enter inference.
    std::vector<AuthoritativePanelSequence> authoritative;
    LocusInputReport report;
};

// `bubbles` must come from a reference-sorted graph, so bubble id order is reference order. Bubbles
// contained in another bubble are dropped: a chain is a partition, and a contained bubble would put
// the same sequence in two blocks.
LocusInputResult build_locus_input(
    const Graph& graph,
    const std::vector<Bubble>& bubbles,
    const LocusInputOptions& options);

// Fail closed. Version one refuses a model whose block boundaries depend on graph representation,
// naming the blocks, rather than scoring a model that a node split would change.
void require_representation_independent(const LocusInputReport& report);

// EVALUATION SUPPORT, never inference. Spells an arbitrary path of `graph` onto the block chain the
// adapter builds from `bubbles`, in canonical chain orientation, applying the same boundary-node
// handling and the recorded divergence trims positionally, so a held-out truth path is expressed in
// exactly the block frame of a model that never saw it. A block is nullopt when the path's walk is
// shorter than the bytes the trim moves out of it, where no positional frame exists. Returns nullopt
// when the path does not span the chain in a consistent direction. `report` must be the report of the
// model built from the same bubbles.
struct PathProjection {
    ChainOrientation orientation = ChainOrientation::Forward;
    std::vector<std::optional<std::string>> blocks;
};

std::optional<PathProjection> project_path_onto_locus(
    const Graph& graph,
    const std::vector<Bubble>& bubbles,
    const LocusInputReport& report,
    const std::string& path_name);

// A digest of the statistical model alone: block order, kind, allele sequences, panel rows and
// template names. It deliberately excludes node identifiers, walk signatures and input record order,
// so splitting, merging or renaming graph nodes leaves it unchanged while changing the block
// boundaries or the panel composition does not.
std::string model_fingerprint(const PreparedGenotypeModel& model);

} // namespace panvar
