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
    // Panel templates removed by an explicit curation decision, by graph path name. A curated template
    // contributes no panel row and no authoritative sequence, so the round trip verifies exactly what
    // the model scores. This is the ONLY way a template may leave the panel: a partial frame still
    // fails model construction closed, because silently completing it is the defect this guards.
    // A name that matches no path in the graph is refused rather than ignored -- a manifest that has
    // drifted from its graph must say so, not quietly curate nothing.
    std::vector<std::string> excluded_templates;
};

struct LocusBlockReport {
    std::string block_id;
    GenotypeBlockKind kind = GenotypeBlockKind::Variant;
    std::size_t bubble_id = 0;             // 0 for a context block
    std::size_t source_alleles = 0;        // distinct walks, before sequence collapse
    std::size_t canonical_alleles = 0;     // distinct sequences, after collapse
    // Curation: source alleles (and distinct sequences) carried only by curated templates. They shaped the block's
    // structure (frame, divergence trim) but are not callable: absent from the block's alleles above.
    std::size_t curation_inactive_source_alleles = 0;
    std::size_t curation_inactive_canonical_alleles = 0;
    std::size_t min_bp = 0;
    std::size_t max_bp = 0;
    std::size_t bypass_templates = 0;      // templates that do not walk this block at all
    // The chain nodes bounding this block. A leading flank has no source and a trailing flank no
    // sink; these are what a template must reach for the block to be covered at all.
    std::string source_node;
    std::string sink_node;
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

// A template whose locus frame is partial, and why. Coverage and representation are reported as
// separate concepts and must not be conflated: `uncovered_blocks` are blocks this template never
// observed, which is a property of the template alone, while the blocks a bypass makes
// representation-dependent are a property of the chain and are reported in the locus report.
struct PartialTemplateReport {
    std::string name;
    std::string reason;
    std::vector<std::string> uncovered_blocks;
    // "left", "right" or "both": the terminal anchor of the chain this template fails to reach.
    // Empty when the frame is partial for a reason other than coverage.
    std::string missing_terminal_anchor;
};

struct LocusInputReport {
    std::size_t variant_blocks = 0;
    std::size_t context_blocks = 0;
    std::size_t total_source_alleles = 0;
    std::size_t total_canonical_alleles = 0;
    std::size_t total_curation_inactive_source_alleles = 0;
    std::size_t total_curation_inactive_canonical_alleles = 0;
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
    // Every template refused a complete frame, with the reason. Parallel in content to
    // `refused_templates`, which keeps only the names.
    std::vector<PartialTemplateReport> partial_templates;
    // Templates removed from the panel by an explicit curation decision, each with the reason it would
    // otherwise have carried. A curated panel has to be auditable: the record says what left and why,
    // never merely that the panel is smaller than the graph.
    struct CuratedExclusion {
        std::string template_name;
        std::string reason;
    };
    std::vector<CuratedExclusion> curated_exclusions;
    // Digest over the curated names, so a status can identify which curation produced a model without
    // quoting the whole list. Empty when nothing was curated. The model fingerprint changes
    // independently of this, because the panel rows it hashes no longer contain the excluded template.
    std::string curation_digest;
};

// M4-B (EXPERIMENTAL, consumed only by copying_transition distance_ref_core_v1_experimental): the reference-coordinate
// geometry of the block chain, computed by the adapter against the FRAME -- every allele a complete-frame template carries
// after the divergence trim, curated templates included, before curation removes any allele from the callable model.
// Positions are on the reference path's axis in canonical chain order (invariant segment 0, reference allele of block 0,
// segment 1, ...). A block's core is its reference interval with the maximal prefix and suffix shared by every frame
// allele removed; a block whose frame has one distinct sequence has an empty core (a NEUTRAL state). Never derived from
// reads, truth, candidates or the callable (curated) allele set. Unavailable geometry is explicit, never substituted.
struct LocusGeometry {
    std::string status = "unavailable: not computed";   // "available", or "unavailable: <reason>"
    std::string reference_path_name;
    std::vector<std::string> block_ids;
    std::vector<long long> ref_start, ref_end;        // the reference allele's interval, bp
    std::vector<long long> core_start, core_end;      // the polymorphic core (core_start == core_end: empty)
    std::vector<long long> position_x2;               // twice the copying state's position (core midpoint), exact
    std::vector<std::size_t> frame_distinct_alleles;  // distinct frame sequences per block (1 = neutral state)
    std::string frame_digest;                         // block ids + every frame allele sequence
    std::string digest;                               // status, reference, intervals, frame digest
    bool available() const noexcept { return status == "available"; }
    // Distance between consecutive core midpoints (bp; exact halves). Throws when unavailable.
    std::vector<double> distances() const;
};

// The geometry rule itself, separated from the graph walk so it can be gated on fixtures. `frame_alleles[b]` are block
// b's frame sequences; `reference_alleles` (one per block) or nullptr with `unavailable_reason`.
LocusGeometry compute_locus_geometry(const std::vector<std::string>& block_ids,
                                     const std::vector<std::string>& invariant_segments,
                                     const std::vector<std::vector<std::string>>& frame_alleles,
                                     const std::vector<std::string>* reference_alleles,
                                     const std::string& reference_path_name,
                                     const std::string& unavailable_reason = "");

struct LocusInputResult {
    StructuralLocusInput input;
    LocusGeometry geometry;   // M4-B
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
