#pragma once

// M3-D (EXPERIMENTAL, default off): one bounded, likelihood-guided joint-block pass from a converged incumbent.
//
// This is a bounded heuristic neighbourhood, not an exhaustive search, not A*, and not a certificate of global
// optimality. The frozen rule (M3-C):
//   1. score every single-block alternative of every block on each homologue, the other homologue fixed, by exact block
//      replacement from the incumbent's census-holding slots;
//   2. shortlist `top` (4) alternatives per block and homologue by fragment + exposure (the prior excluded), ties by
//      allele index;
//   3. enumerate the shortlisted singles and every same-homologue combination of two different blocks (all block pairs),
//      whatever the single-block intermediates score;
//   4. rank the resulting unordered sequence pairs, with the incumbent, by the complete score (fragment + exposure + the
//      unchanged sequence-level pair prior) under the repair round's own certified move rule.
// Both bounds (single-block alternatives, catalogue candidates; 2,048 each) are checked from the model before anything is
// scored; a pass over either bound is skipped, never truncated. Shortlisting scores are not candidates: only catalogue
// members (and the incumbent) are reported and ranked.
//
// Everything runs in the caller's PreparedScoringContext (never a second live context), through its existing exact
// machinery: add_repaired_haplotypes (planned, batched), add_repaired_haplotype(..., retain_census = true) to promote a
// shortlisted single (census by delta), score_candidates. Caps are the context's own cumulative caps; a cap crossed
// throws ResourceExhausted out of the pass, which the caller handles under its incomplete contract.

#include "panvar/fragment_likelihood.hpp"
#include "panvar/genotype_model.hpp"

#include <cstddef>
#include <cstdint>
#include <functional>
#include <string>
#include <vector>

namespace panvar {

struct JointSearchOptions {
    std::size_t top = 4;
    std::size_t single_limit = 2048;
    std::size_t catalogue_limit = 2048;
    double switch_probability = 0.01;
    // GATE HOOK, never production: permute the single and joint request lists (results must not change).
    std::uint64_t request_order_seed = 0;
};

struct JointSearchBounds {
    std::uint64_t single_alternatives = 0;   // sum over blocks of 2 * (alleles - 1)
    std::uint64_t catalogue_raw = 0;         // shortlisted singles + same-homologue block pairs, before deduplication
    bool within = false;
};
// Model-only: computable before any read is touched.
JointSearchBounds joint_search_bounds(const PreparedGenotypeModel& model, const JointSearchOptions& options);

struct JointCandidate {
    std::vector<std::size_t> first, second;   // allele paths
    std::size_t slot_first = 0, slot_second = 0;
    DiplotypeScore score;
    double prior = 0.0;
    std::string move;                         // first generating move; "single:h<k>:<block>:<from>><to>:r<rank>" or "pair:..."
    std::size_t equivalent_moves = 1;
};

struct JointSearchResult {
    std::string status = "not_run";           // skipped_bounds | complete
    std::string decision = "none";            // moved | no_improvement | tie_at_best | rescore_incomplete | skipped
    JointSearchBounds bounds;
    std::size_t singles_scored = 0, promoted = 0, joint_requests = 0;
    std::vector<JointCandidate> catalogue;    // unique sequence pairs, incumbent excluded, generation order
    double incumbent_total = 0.0;
    bool moved = false;
    std::size_t winner = 0;                   // index into catalogue when moved
    std::size_t winner_set_size = 0;
    bool incumbent_in_top_set = false;
    double certified_gap = 0.0;
    // work and resources
    std::uint64_t origins_reused = 0, origins_recomputed = 0, dp_executions = 0, window_bases = 0;
    std::uint64_t full_placements = 0;        // during the pass (expected 0)
    std::uint64_t census_bytes_peak = 0;      // retained census bytes after promotion
    std::uint64_t transient_bytes_peak = 0;   // transient payload after the joint batch
    std::uint64_t transient_bytes_after = 0;  // after the pass released its slots
    std::size_t slots_created = 0, live_transient_after = 0;
    double seconds_singles = 0.0, seconds_promotion = 0.0, seconds_joint = 0.0, seconds_scoring = 0.0,
           seconds_priors = 0.0, seconds_total = 0.0;
};

// The pass. `slot_first` / `slot_second` must hold censuses (a converged repair search's incumbent slots). `on_scored`
// receives every catalogue member once the whole catalogue is scored (never a partial ranking), with its spelled bytes.
// `keep` lists slots the release must preserve (the caller's pins). When the pass moves, the winner's slots keep their
// censuses and are preserved too; every other slot the pass created is released before it returns.
JointSearchResult run_joint_block_search(
    PreparedScoringContext& context, const PreparedGenotypeModel& model, const std::vector<std::size_t>& first,
    const std::vector<std::size_t>& second, std::size_t slot_first, std::size_t slot_second, double incumbent_total,
    std::size_t fragments, std::size_t mismatch_band, const JointSearchOptions& options,
    const std::function<void(const JointCandidate&, const std::string&, const std::string&)>& on_scored,
    const std::vector<std::size_t>& keep);

}  // namespace panvar
