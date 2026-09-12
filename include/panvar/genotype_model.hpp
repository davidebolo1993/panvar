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
struct LibraryModel {
    double insert_mean_bp = 0.0;
    double insert_sd_bp = 0.0;
    double support_sigmas = 0.0;
    std::string provenance;
};

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

void validate_read_pair(const ReadPair& fragment);
void validate_library_model(const LibraryModel& model);

// Insert length is the template span from the first outer read boundary to the second. Overlapping
// mates are physical states, so the lower limit is the longer observed mate, not their length sum.
// The Gaussian support is intersected with that fragment-specific physical limit.
InsertSupport fragment_insert_support(
    const ReadPair& fragment,
    const LibraryModel& model);

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

struct GenotypeBlockInput {
    std::string block_id;
    std::vector<BlockAlleleInput> alleles;
};

struct PanelTemplateInput {
    std::string name;
    TemplateFrameStatus frame_status = TemplateFrameStatus::Complete;
    // Exactly one block-local source allele identifier for every ordered block.
    std::vector<std::string> allele_source_ids;
};

struct StructuralLocusInput {
    std::string locus_id;
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
        PanelAlleleMatrix panel);

    std::string locus_id_;
    ChainOrientation chain_orientation_ = ChainOrientation::Forward;
    std::vector<PreparedGenotypeBlock> blocks_;
    std::vector<std::string> invariant_segments_;
    PanelAlleleMatrix panel_;
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
