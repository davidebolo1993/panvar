#pragma once

#include "panvar/genotype_model.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace panvar {

struct MarkerModelParameters {
    std::size_t kmer_size = 0;
    std::size_t syncmer_s = 0;
    double expected_observations_per_copy = 0.0;
    double background_observations = 0.0;
    std::string provenance;
};

struct BlockMarkerInventory {
    std::string block_id;
    std::size_t allele_count = 0;
    // Local marker index -> slot in PreparedMarkerPanel::codes and fixed copies elsewhere on one
    // haplotype. The varying block's copies are stored sparsely per allele below.
    std::vector<std::uint32_t> marker_slots;
    std::vector<std::uint32_t> fixed_haploid_multiplicity;
    std::vector<std::vector<std::pair<std::uint32_t, std::uint32_t>>> allele_markers;
};

struct PreparedMarkerPanel {
    std::size_t kmer_size = 0;
    std::size_t syncmer_s = 0;
    std::vector<std::uint64_t> codes;
    std::vector<BlockMarkerInventory> blocks;
    std::vector<std::uint32_t> anchor_slots;
    std::vector<std::uint32_t> anchor_haploid_multiplicity;
    std::size_t total_alleles = 0;
    std::size_t scored_alleles = 0;
    std::size_t dropped_multiblock_markers = 0;
};

struct MarkerObservations {
    std::vector<std::uint64_t> counts;
    std::uint64_t reads = 0;
    std::uint64_t bases = 0;
    std::uint64_t syncmer_occurrences = 0;
    std::uint64_t matched_occurrences = 0;
};

struct BlockProposalEmissions {
    std::size_t allele_count = 0;
    std::size_t scored_alleles = 0;
    std::vector<double> log_score;
    // Full Poisson score at zero variable copies. It is omitted from log_score because it is the
    // same for every allele pair in this block, but retained so that omission is auditable.
    double candidate_independent_log_constant = 0.0;

    double at(std::size_t allele1, std::size_t allele2) const;
};

struct PreparedProposalEmissions {
    std::vector<BlockProposalEmissions> blocks;
    std::size_t total_alleles = 0;
    std::size_t scored_alleles = 0;

    bool every_allele_scored() const noexcept { return total_alleles == scored_alleles; }
};

void validate_marker_parameters(const MarkerModelParameters& parameters);

PreparedMarkerPanel build_marker_panel(
    const PreparedGenotypeModel& model,
    std::size_t kmer_size,
    std::size_t syncmer_s);

MarkerObservations count_marker_observations(
    const PreparedMarkerPanel& panel,
    const std::vector<ReadPair>& fragments);

PreparedProposalEmissions build_proposal_emissions(
    const PreparedMarkerPanel& panel,
    const MarkerObservations& observations,
    const MarkerModelParameters& parameters);

} // namespace panvar
