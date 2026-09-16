#pragma once

#include <cstddef>
#include <map>
#include <string>
#include <vector>

namespace panvar {

// A read exactly as observed. Qualities may be empty when the input format has none; when present,
// they must have one byte per base. The replacement caller does not assume one read length per run.
struct SequencingRead {
    std::string name;
    std::string bases;
    std::string qualities;
};

struct ReadPair {
    SequencingRead first;
    SequencingRead second;
};

// The declared paired-end library model. Estimation and explicit overrides will both produce this
// object, so every downstream scorer receives the same validated policy and records its provenance.
// There is one of these, not one per scorer: a private copy of any of these parameters is how two
// stages start disagreeing about the same library.
struct LibraryModel {
    double insert_mean_bp = 0.0;
    double insert_sd_bp = 0.0;
    double support_sigmas = 0.0;
    // Per-base substitution probability. Version one is substitutions only, with a fixed-position
    // Hamming emission, so this is the whole error model.
    double substitution_error = 0.0;
    // Robust outlier mixture weight, eta in the fragment contract. A fragment that belongs to no
    // candidate -- a chimera, a mismapped pair, contamination -- is absorbed here instead of driving
    // the score of every candidate to minus infinity together.
    double outlier_weight = 0.0;
    // Fragment-start intensity per haplotype base: lambda in the fragment contract. This is a
    // fragment-start rate, not a marker depth, and the two must not be conflated.
    double fragment_start_intensity = 0.0;
    std::string provenance;
};

// The candidate-independent background probability of observing this fragment's bases at all, under
// a uniform base model over both mates. It is the P_background(f) of the fragment contract: it
// depends on the fragment's own read lengths and on nothing about any candidate, so it can never
// reorder candidates -- which is exactly why it is safe to include.
double log_background_probability(const ReadPair& fragment);

struct InsertSupport {
    std::size_t min_bp = 1;
    std::size_t max_bp = 0;

    bool empty() const noexcept { return min_bp > max_bp; }
};

struct DiscreteInsertDistribution {
    InsertSupport support;
    std::vector<double> log_probability_by_offset;

    double log_probability(std::size_t insert_bp) const noexcept;
};

// The one definition of log-space addition. It sits in this lowest layer because both the insert
// normalisation here and the fragment scorer above need it, and two copies of it is exactly how two
// stages start disagreeing about the same sum.
double log_add_exp(double lhs, double rhs);

void validate_read_pair(const ReadPair& fragment);
void validate_library_model(const LibraryModel& model);

// Insert length is the template span from the first outer read boundary to the second. Overlapping
// mates are physical states, so the lower limit is the longer observed mate, not their length sum.
// The Gaussian support is intersected with that fragment-specific physical limit.
InsertSupport fragment_insert_support(
    const ReadPair& fragment,
    const LibraryModel& model);

// The library's own statistical insert support, before any fragment's physical floor is applied.
InsertSupport library_insert_support(const LibraryModel& model);

// The discrete Gaussian normalised over the library's statistical support. This is the insert
// probability of the fragment likelihood: the physical floor of an individual fragment restricts
// which placements are geometrically possible, but it does not renormalise the probability, because
// a template shorter than the longer mate is a state of the library that cannot produce that
// fragment, not a state that was never generated. Exposure uses the same distribution.
DiscreteInsertDistribution build_library_insert_distribution(const LibraryModel& model);

// Build once per fragment (or per equal mate-length bucket), then reuse for every candidate and
// placement. This prevents normalization over the support from entering the scoring inner loop.
// Probabilities are a discrete Gaussian normalized over this fragment's declared insert support.
DiscreteInsertDistribution build_fragment_insert_distribution(
    const ReadPair& fragment,
    const LibraryModel& model);

// Panel template h carries template_alleles[h][b] at block b. Allele identifiers are block-local.
struct PanelAlleleMatrix {
    std::vector<std::string> template_names;
    std::vector<std::vector<std::size_t>> template_alleles;

    std::size_t template_count() const noexcept { return template_alleles.size(); }
    std::size_t block_count() const noexcept {
        return template_alleles.empty() ? 0 : template_alleles.front().size();
    }
};

// Structural inputs are intentionally separate from truth and experiment metadata. A graph-facing
// adapter may populate these records later, but the statistical model only accepts this narrow,
// auditable representation.
enum class ChainOrientation { Forward, Reverse };
enum class TemplateFrameStatus { Complete, Partial };

struct BlockAlleleInput {
    // Block-local identifier of the graph walk/record carrying this sequence. Several source records
    // may spell identical sequence; they are retained as aliases of one canonical allele.
    std::string source_id;
    std::string sequence;
};

// A locus is genotyped as an alternating chain: a variant block is a called bubble, and a context
// block is the stretch between two consecutive bubbles (or the flanks beyond the outermost ones).
// Context blocks are not filler. They carry the sub-threshold variation the SV size filter kept out
// of the bubble list, which is what identifies which panel templates a sample sits on. Only variant
// blocks are reporting units; both kinds are modelled identically.
enum class GenotypeBlockKind { Variant, Context };

struct GenotypeBlockInput {
    std::string block_id;
    std::vector<BlockAlleleInput> alleles;
    GenotypeBlockKind kind = GenotypeBlockKind::Variant;
};

struct PanelTemplateInput {
    std::string name;
    TemplateFrameStatus frame_status = TemplateFrameStatus::Complete;
    // Exactly one block-local source allele identifier for every ordered block.
    std::vector<std::string> allele_source_ids;
    // Orientation of this template's authoritative source path relative to the canonical block
    // chain. It is retained per template because one panel may contain both parallel and
    // antiparallel paths.
    ChainOrientation source_orientation = ChainOrientation::Forward;
    // Why the frame is partial, for the refusal message; empty when the frame is complete. Last,
    // because this is an aggregate that fixtures initialise positionally: a new field inserted
    // ahead of `allele_source_ids` silently swallows the allele list and compiles clean.
    std::string partial_reason;
};

struct StructuralLocusInput {
    std::string locus_id;
    // Blocks, alleles and invariant segments below are already expressed in canonical chain order
    // and orientation. This records how that canonical chain relates to the external reference
    // coordinate frame; it must not cause statistical code to reinterpret their bytes.
    ChainOrientation chain_orientation = ChainOrientation::Forward;
    std::vector<GenotypeBlockInput> blocks;
    // There are blocks+1 invariant segments: before block 0, between each adjacent pair, and after
    // the final block. Empty segments are valid and are not missing data.
    std::vector<std::string> invariant_segments;
    std::vector<PanelTemplateInput> templates;
};

struct CanonicalBlockAllele {
    std::size_t id = 0;
    std::string stable_id;
    std::string sequence;
    std::vector<std::string> source_ids;
    std::size_t panel_carrier_count = 0;
};

struct PreparedGenotypeBlock {
    std::string block_id;
    GenotypeBlockKind kind = GenotypeBlockKind::Variant;
    std::vector<CanonicalBlockAllele> alleles;
    std::map<std::string, std::size_t> allele_id_by_source;

    std::size_t allele_id_for_source(const std::string& source_id) const;
};

// Immutable after construction: only const accessors expose its components. Construction validates
// the complete panel frame first, canonicalizes sequence-identical alleles, and assigns IDs from
// sequence order rather than input record order.
class PreparedGenotypeModel {
public:
    const std::string& locus_id() const noexcept { return locus_id_; }
    ChainOrientation chain_orientation() const noexcept { return chain_orientation_; }
    TemplateFrameStatus frame_status() const noexcept { return TemplateFrameStatus::Complete; }
    const std::vector<PreparedGenotypeBlock>& blocks() const noexcept { return blocks_; }
    const std::vector<std::string>& invariant_segments() const noexcept {
        return invariant_segments_;
    }
    const PanelAlleleMatrix& panel() const noexcept { return panel_; }
    const std::vector<ChainOrientation>& template_source_orientations() const noexcept {
        return template_source_orientations_;
    }

    const CanonicalBlockAllele& allele(
        std::size_t block_id,
        std::size_t allele_id) const;

private:
    friend PreparedGenotypeModel prepare_genotype_model(const StructuralLocusInput& input);

    PreparedGenotypeModel(
        std::string locus_id,
        ChainOrientation chain_orientation,
        std::vector<PreparedGenotypeBlock> blocks,
        std::vector<std::string> invariant_segments,
        PanelAlleleMatrix panel,
        std::vector<ChainOrientation> template_source_orientations);

    std::string locus_id_;
    ChainOrientation chain_orientation_ = ChainOrientation::Forward;
    std::vector<PreparedGenotypeBlock> blocks_;
    std::vector<std::string> invariant_segments_;
    PanelAlleleMatrix panel_;
    std::vector<ChainOrientation> template_source_orientations_;
};

PreparedGenotypeModel prepare_genotype_model(const StructuralLocusInput& input);

struct MosaicPath {
    std::vector<std::size_t> alleles;

    bool operator==(const MosaicPath& other) const noexcept { return alleles == other.alleles; }
    bool operator!=(const MosaicPath& other) const noexcept { return !(*this == other); }
};

struct MosaicDiplotype {
    MosaicPath first;
    MosaicPath second;
};

void validate_panel(const PanelAlleleMatrix& panel);

} // namespace panvar
