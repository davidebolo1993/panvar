#pragma once

#include "panvar/genotype_model.hpp"
#include "panvar/marker_model.hpp"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

namespace panvar {

// Deterministic top-K mosaic proposal with equivalence merging in every layer.
//
// Two template histories that have spelled the same allele prefix on a homologue and sit on the same
// current template have identical futures: every later emission and transition depends only on the
// prefix already spelled and the current template. The beam therefore keys a state by the UNORDERED pair
// of (allele-prefix id, current template), merges every history reaching that key before pruning, and
// carries the log-summed proposal mass forward. The pairing of each prefix with its own template is part
// of the key: canonicalising the template pair and the prefix pair separately would exchange homologues.
//
// Without this, hundreds of histories describing one mosaic occupied the width and crowded out distinct
// mosaics, so the truth was never proposed at widths that a single merged state would have fitted in.
struct MosaicBeamParameters {
    std::size_t unique_candidate_target = 0;
    // Merged states retained per block. A resource cap, not a model parameter: a search that prunes and
    // still misses the target is incomplete, however large this is.
    std::size_t state_limit = 0;
    std::uint64_t max_extension_evaluations = std::numeric_limits<std::uint64_t>::max();
};

struct MosaicCandidate {
    std::size_t candidate_id = 0;
    std::size_t unique_rank = 0;
    MosaicDiplotype diplotype;
    // The best single retained history spelling this diplotype, and its two components.
    double best_proposal_score = 0.0;
    double proposal_emission_component = 0.0;
    double proposal_transition_component = 0.0;
    // log sum over every retained history spelling this diplotype.
    double retained_logsum_proposal_score = 0.0;
    // Physical ordered template-pair histories merged into this candidate. Saturates at 2^64-1 rather
    // than wrapping; the flag records when it did.
    std::uint64_t equivalent_retained_histories = 0;
    bool history_count_saturated = false;
    bool every_allele_scored = false;
};

struct MosaicBeamResult {
    std::vector<MosaicCandidate> candidates;   // at most unique_candidate_target, best retained mass first
    std::size_t state_limit_used = 0;
    std::size_t final_states = 0;
    std::size_t unique_candidates_found = 0;   // distinct diplotypes among final states, before truncation
    std::uint64_t states_pruned = 0;           // merged states discarded, summed over every block
    std::uint64_t extensions_evaluated = 0;
    std::uint64_t extension_upper_bound = 0;
    // No merged state was ever pruned, so every diplotype with positive proposal mass was found.
    bool search_exhaustive = false;
    bool unique_target_reached = false;
    // The returned set omits candidates -- because states were pruned, or because more diplotypes were
    // found than the target -- so any ranking over it is conditional on the candidate set.
    bool candidate_set_conditional = false;
};

MosaicBeamResult propose_mosaic_diplotypes(
    const PanelAlleleMatrix& panel,
    const PreparedProposalEmissions& emissions,
    const std::vector<double>& switch_probability,
    const MosaicBeamParameters& parameters);

} // namespace panvar
