#pragma once

#include "panvar/genotype_model.hpp"
#include "panvar/marker_model.hpp"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <utility>
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
    // Diagnostic only: follow one declared diplotype through every layer. Its prefixes are looked up and
    // never interned, so prefix ids, ordering and pruning are identical to an untraced search.
    std::optional<MosaicDiplotype> trace_target = std::nullopt;
    // Optionally, the template pair (first path's template, second path's template) whose own state is
    // also followed -- in leave-zero-out, the donor's two panel rows.
    std::optional<std::pair<std::size_t, std::size_t>> trace_templates = std::nullopt;
};

// One layer of a trace: what was offered to pruning, what survived, and where the target stood.
struct MosaicBeamTraceLayer {
    std::size_t block = 0;
    std::uint64_t states_offered = 0;      // merged states offered to the layer's selection
    std::size_t states_retained = 0;
    // Distinct allele-prefix pairs among the retained states: how many diplotype prefixes the width holds.
    std::size_t prefix_pairs_retained = 0;
    double retained_min_log_mass = -std::numeric_limits<double>::infinity();
    // Both target prefixes exist as interned prefixes at this layer.
    bool target_prefixes_present = false;
    std::uint64_t target_states_offered = 0;   // states spelling the target prefix pair, any templates
    std::uint64_t target_states_retained = 0;
    double target_best_log_mass = -std::numeric_limits<double>::infinity();
    double target_logsum_log_mass = -std::numeric_limits<double>::infinity();
    std::uint32_t target_best_template1 = 0;   // templates of the best target state, as stored
    std::uint32_t target_best_template2 = 0;
    // Offered states with strictly more mass than the best target state: its rank minus one.
    std::uint64_t offered_above_target_best = 0;
    double declared_templates_log_mass = -std::numeric_limits<double>::infinity();
    bool declared_templates_offered = false;
    bool declared_templates_retained = false;
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
    std::vector<MosaicBeamTraceLayer> trace;   // one per block when a trace target was declared
};

MosaicBeamResult propose_mosaic_diplotypes(
    const PanelAlleleMatrix& panel,
    const PreparedProposalEmissions& emissions,
    const std::vector<double>& switch_probability,
    const MosaicBeamParameters& parameters);

} // namespace panvar
