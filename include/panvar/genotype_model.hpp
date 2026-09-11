#pragma once

#include <cstddef>
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
