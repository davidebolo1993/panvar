#pragma once

#include "panvar/genotype_model.hpp"

#include <cstddef>
#include <cstdint>
#include <optional>
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

// Per-copy marker depth, estimated from the invariant anchors. Every panel haplotype carries an
// anchor's fixed copies, so its diploid multiplicity is known whatever the sample's genotype, and
// count_k ~ Poisson(background + rate * 2 * m_k). A variable marker cannot be used for this: its
// expected count depends on the genotype being inferred.
//
// With a declared background the rate is the one-dimensional Poisson maximum-likelihood root. Without
// one, both are fitted jointly by damped Newton, identifiable only when the anchors carry at least two
// distinct multiplicities; with a single multiplicity background and rate trade off exactly, so the
// fit refuses and the background must be declared.
struct MarkerDepthEstimate {
    double rate_per_copy = 0.0;
    double background = 0.0;
    bool background_estimated = false;
    bool converged = false;
    // The fitted background reached its positive floor: the data prefer no background at all.
    bool background_at_floor = false;
    // Why the fit stopped: interior_stationary, boundary_kkt, closed_form, declared_background_root,
    // no_ascent, or iteration_cap. Only the first four are convergence.
    std::string termination;
    // Score components and log-likelihood at the returned estimate, so a status file can show a fit is
    // stationary rather than assert it. On the boundary the background gradient is negative by design.
    double gradient_background = 0.0;
    double gradient_rate = 0.0;
    double log_likelihood = 0.0;
    std::size_t anchors_used = 0;
    std::size_t distinct_multiplicities = 0;
    std::size_t iterations = 0;
    std::string provenance;
};

MarkerDepthEstimate estimate_marker_depth(
    const PreparedMarkerPanel& panel,
    const MarkerObservations& observations,
    std::optional<double> declared_background);

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
