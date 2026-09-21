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
    // Grouped proposer only: haplotype-prefix vector cells held at once. Checked before anything is
    // allocated, so a refused search allocates nothing.
    std::uint64_t max_prefix_vector_cells = std::numeric_limits<std::uint64_t>::max();
    // Diagnostic only: follow one declared diplotype through every layer. Its prefixes are looked up and
    // never interned, so prefix ids, ordering and pruning are identical to an untraced search.
    std::optional<MosaicDiplotype> trace_target = std::nullopt;
    // Optionally, the template pair (first path's template, second path's template) whose own state is
    // also followed -- in leave-zero-out, the donor's two panel rows.
    std::optional<std::pair<std::size_t, std::size_t>> trace_templates = std::nullopt;
};

// The beam's pre-allocation budget: templates^2 for the first block, plus state_limit * templates^2 for
// every block after it, in integers with every multiplication overflow-checked. Production compares
// MosaicBeamParameters::max_extension_evaluations against this, so anything that PREDICTS the same work
// must call it rather than restate the formula -- a second version in floating point is not the same
// function near the cap, which is the only place the answer decides anything. Throws std::length_error
// when the bound itself overflows uint64.
std::uint64_t extension_upper_bound(std::size_t templates, std::size_t blocks, std::size_t state_limit);


// One layer of a trace: what was offered to pruning, what survived, and where the target stood.
struct MosaicBeamTraceLayer {
    std::size_t block = 0;
    std::uint64_t states_offered = 0;      // merged states offered to the layer's selection
    std::size_t states_retained = 0;
    // Distinct allele-prefix pairs among the states OFFERED to this layer's selection: how many
    // diplotype prefixes the layer had to choose between. Read against prefix_pairs_retained, it
    // separates a layer that was offered little diversity from one that was offered plenty and kept
    // little.
    std::size_t prefix_pairs_offered = 0;
    // Distinct allele-prefix pairs among the retained states: how many diplotype prefixes the width holds.
    std::size_t prefix_pairs_retained = 0;
    // The most retained states sharing a single prefix pair. The width is spent on DIPLOTYPE PREFIXES
    // only to the extent this stays small: a layer whose width is largely one prefix's template
    // variants is exploring one history many times over, which is what prefix_pairs_retained alone
    // cannot show -- it counts the prefixes present, not how unevenly the width is divided among them.
    std::size_t max_states_per_prefix = 0;
    // States this layer's selection discarded: offered minus retained. The result's `states_pruned` is
    // the sum over blocks and says nothing about WHERE the search narrowed.
    std::uint64_t states_pruned = 0;
    double retained_min_log_mass = -std::numeric_limits<double>::infinity();
    // Both target prefixes exist as interned prefixes at this layer.
    bool target_prefixes_present = false;
    std::uint64_t target_states_offered = 0;   // states spelling the target prefix pair, any templates
    std::uint64_t target_states_retained = 0;
    double target_best_log_mass = -std::numeric_limits<double>::infinity();
    double target_logsum_log_mass = -std::numeric_limits<double>::infinity();
    // Diagnostic-only target accounting: the target prefix group's mass as the layer RETAINED it, and
    // what the layer dropped -- log(offered - retained). Equal offered and retained masses report -inf
    // lost, which is the honest reading of "nothing was lost".
    double target_retained_logsum_log_mass = -std::numeric_limits<double>::infinity();
    double target_lost_log_mass = -std::numeric_limits<double>::infinity();
    // Rank of the target's prefix GROUP among every group offered, by aggregate group mass -- the
    // quantity prefix-first selection actually ranks on. 0 when the target was not offered. A rank
    // among states (offered_above_target_best) cannot stand in for it.
    std::size_t target_group_rank_by_offered_mass = 0;
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

// The protected panel baseline: every unordered pair of panel templates, deduplicated by allele paths, with
// the first path not lexicographically greater than the second and ids in enumeration order. In
// leave-zero-out the truth is one of these, and it is a candidate whatever the beam pruned.
std::vector<MosaicCandidate> panel_diplotypes(const PanelAlleleMatrix& panel);

// The beam's diplotypes that are not panel pairs, in proposal order, ranked 1.. with ids following the panel
// baseline; unique_candidates_found counts them. The beam must have returned every diplotype it found.
MosaicBeamResult off_panel_mosaics(const MosaicBeamResult& beam, const std::vector<MosaicCandidate>& panel);

MosaicBeamResult propose_mosaic_diplotypes(
    const PanelAlleleMatrix& panel,
    const PreparedProposalEmissions& emissions,
    const std::vector<double>& switch_probability,
    const MosaicBeamParameters& parameters);
// The grouped-prefix proposer: an EXPERIMENTAL second implementation, retaining allele-prefix pairs and
// carrying each one's whole distribution over current template pairs instead of a single attribution of
// it. propose_mosaic_diplotypes above is untouched and is the oracle this is tested against; nothing
// here is wired into calling. Same inputs, same MosaicBeamResult, same candidate ordering.
MosaicBeamResult propose_mosaic_diplotypes_grouped(const PanelAlleleMatrix& panel,
                                                   const PreparedProposalEmissions& emissions,
                                                   const std::vector<double>& switch_probability,
                                                   const MosaicBeamParameters& parameters);

// The two pre-allocation counts, taken from the exact loops the grouped proposer runs: one template
// scan per retained haplotype prefix (at most two per retained pair) plus one child per ordered allele
// pair per retained pair, and the vector cells held across a layer boundary. Both are integer and
// overflow-checked, and both are checked before anything is allocated.
std::uint64_t grouped_extension_upper_bound(std::size_t templates, std::size_t state_limit,
                                            const std::vector<std::size_t>& allele_counts);
std::uint64_t grouped_prefix_vector_cells(std::size_t templates, std::size_t state_limit);

} // namespace panvar
