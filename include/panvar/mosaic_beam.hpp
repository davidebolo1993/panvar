#pragma once

#include "panvar/genotype_model.hpp"
#include "panvar/marker_model.hpp"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

namespace panvar {

struct MosaicBeamParameters {
    std::size_t unique_candidate_target = 0;
    std::size_t raw_history_limit = 0;
    std::size_t expansion_chunk_size = 262144;
    std::uint64_t max_extension_evaluations = std::numeric_limits<std::uint64_t>::max();
};

struct MosaicCandidate {
    std::size_t candidate_id = 0;
    std::size_t unique_rank = 0;
    MosaicDiplotype diplotype;
    double best_proposal_score = 0.0;
    double retained_logsum_proposal_score = 0.0;
    double proposal_emission_component = 0.0;
    double proposal_transition_component = 0.0;
    std::size_t equivalent_retained_histories = 0;
    bool every_allele_scored = false;
};

struct MosaicBeamResult {
    std::vector<MosaicCandidate> candidates;
    std::size_t beam_width_used = 0;
    std::size_t final_raw_histories = 0;
    std::uint64_t extensions_evaluated = 0;
    std::uint64_t extension_upper_bound = 0;
    bool unique_target_reached = false;
    bool raw_history_limit_reached = false;
};

MosaicBeamResult propose_mosaic_diplotypes(
    const PanelAlleleMatrix& panel,
    const PreparedProposalEmissions& emissions,
    const std::vector<double>& switch_probability,
    const MosaicBeamParameters& parameters);

} // namespace panvar
