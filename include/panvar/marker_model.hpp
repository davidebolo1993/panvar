#pragma once

#include "panvar/genotype_model.hpp"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace panvar {

// The marker model's one background floor, in expected observations per marker. The depth fit keeps
// its background at or above it, a declared background below it is refused, and proposal emissions
// refuse a background below it, so the background the fit reports is the background proposal scoring
// uses. It is not numerically negligible: at an observed marker that a candidate gives zero copies,
// each observation contributes log(background) to that candidate's score.
inline constexpr double kMarkerBackgroundFloor = 1e-6;

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
// With a declared background (at least kMarkerBackgroundFloor) the rate is the one-dimensional Poisson
// maximum-likelihood root. Without one, both are fitted jointly by active-set Newton, identifiable only
// when the anchors carry at least two distinct multiplicities; with a single multiplicity background and
// rate trade off exactly, so the fit refuses and the background must be declared.
struct MarkerDepthEstimate {
    double rate_per_copy = 0.0;
    double background = 0.0;
    bool background_estimated = false;
    bool converged = false;
    // The fitted background reached its positive floor: the data prefer no background at all.
    bool background_at_floor = false;
    // Why the fit stopped: interior_stationary, boundary_kkt, declared_background_root, no_ascent, or
    // iteration_cap. Only the first three are convergence.
    std::string termination;
    // Score components and log-likelihood at the returned estimate, so a status file can show a fit is
    // stationary rather than assert it. On the boundary the background gradient is negative by design.
    double gradient_background = 0.0;
    double gradient_rate = 0.0;
    double log_likelihood = 0.0;
    // Joint fit only: the floating-point resolution of the compensated log-likelihood, and the gain the
    // Newton model still predicted when the fit stopped (on the boundary, the rate-only face's). Convergence
    // means the second is no larger than the first.
    double likelihood_resolution = 0.0;
    double newton_decrement = 0.0;
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

// Proposal parameters from a converged depth fit, with the fitted or declared background unchanged.
MarkerModelParameters proposal_marker_parameters(
    const MarkerDepthEstimate& depth,
    std::size_t kmer_size,
    std::size_t syncmer_s);

// I3-D: completed-or-attempted build_marker_panel executions in this process (tests prove a route never built one).
std::uint64_t marker_panel_builds();

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

// ---- I2-E (EXPERIMENTAL): exact reuse of marker-proposal contributions -------------------------------------------------
// Within one block of one call, a marker's contribution to an allele pair -- poisson_log_kernel(count, mean) minus the
// marker's baseline kernel, mean = background + per_copy * (2 * fixed + copies) -- depends only on the marker and the
// pair's combined variable copies, so the allele-pair loop recomputes the same term many times. With reuse on, each
// block holds a lazily filled table: per marker, entries for copies 1 .. 2 * (the marker's largest copy count among the
// block's alleles), each evaluated once by the SAME helper the direct path uses and then served. The traversal, the
// sparse merge and the summation order are the direct path's own (one shared loop). The table is per call and per block
// (built from this call's observations and parameters, released after the block); its budget bounds the RETAINED
// storage only -- the allocated capacity of the table and its two per-marker index vectors -- and markers whose range
// does not fit (admitted in marker order) are evaluated directly. It does not bound the temporary planning arrays (each
// marker's largest copy count and its tentative offset) that every enabled block allocates to decide admission, even
// when nothing is admitted; those are reported apart (peak_planning_bytes). Off (the default and the three-argument
// overload): every term evaluated directly, and nothing is allocated for reuse.
struct ProposalReuseOptions {
    bool enabled = false;
    std::size_t table_budget_bytes = 0;   // per block
};
struct ProposalReuseReport {
    bool enabled = false;
    std::uint64_t table_budget_bytes = 0;
    std::uint64_t visits = 0;     // marker visits in the allele-pair merge: one contribution each
    std::uint64_t executed = 0;   // contributions evaluated (table fills plus direct evaluations)
    std::uint64_t reused = 0;     // contributions served from a table
    std::uint64_t fallback = 0;   // of executed: direct because the marker's range was not tabled
    std::uint64_t markers_tabled = 0, markers_not_tabled = 0;
    std::uint64_t peak_table_bytes = 0;     // the largest block's retained table and index, allocated capacity (budgeted)
    std::uint64_t peak_planning_bytes = 0;  // the largest block's temporary planning arrays (not budgeted)
    double seconds_table_build = 0.0;     // range scan and allocation, summed over blocks
    double seconds_total = 0.0;
};
PreparedProposalEmissions build_proposal_emissions(
    const PreparedMarkerPanel& panel,
    const MarkerObservations& observations,
    const MarkerModelParameters& parameters,
    const ProposalReuseOptions& reuse,
    ProposalReuseReport* report = nullptr);

} // namespace panvar
