#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace panvar {

// The only public surface of the mosaic caller. Debug probes live in unit tests or a developer
// executable, not here: `genotype-frag`, `--hybrid-call` and `--mosaic-rescore` are deliberately
// gone and must not come back.
struct GenotypeCommandOptions {
    std::string gfa_path;
    // Optional. A bubbles CSV from `panvar bubble`; when absent the bubbles are called internally
    // with the same SV threshold, so the block alphabet is the same either way.
    std::string bubbles_path;
    std::string reads1_path;
    std::string reads2_path;
    std::string reference_path;
    std::string out_prefix = "genotype";
    // Declared library parameters. Version one does not estimate the insert distribution from the
    // data, so this is required and its provenance is recorded in the status file.
    std::string library_path;
    // Optional. A curation manifest: one panel template name per line, '#' comments ignored. Named
    // templates are removed from the panel before the model is built, and the run records which and
    // why. Curation is the only sanctioned way a template leaves the panel -- a partial frame still
    // fails model construction closed -- and a manifest naming an absent template is refused.
    std::string curation_manifest_path;
    std::size_t threads = 1;
    // Off-panel mosaic diplotypes added to the protected panel baseline; 0 scores the panel alone. Every
    // unordered panel diplotype is always a candidate and is not counted here.
    // FROZEN at 128 by the declared selection rule (width ladder + K ladder, five frozen C4 plans):
    // the smallest K that exactly scores every development truth, gives the same block calls or
    // tied-winner block consensus as K = 256, and is stable at the next larger proposal width. K = 8,
    // 32 and 64 all fail the first clause alone, because every_block's truth sits at proposal rank 87.
    std::size_t max_mosaic_candidates = 128;

    std::size_t kmer_size = 31;
    std::size_t syncmer_s = 15;
    std::size_t min_variant_bp = 50;
    double switch_probability = 0.01;
    // The declared mismatch band of the truncated fragment likelihood. It is part of the model, not a
    // search depth, so it is fixed for the whole run and recorded in the status file.
    std::size_t mismatch_band = 6;
    std::size_t minimum_piece_bp = 12;
    // Declared resource budgets. Reaching one ends the run incomplete with no ranking and no call; the
    // values are recorded in the status file. Not command-line options: the public surface is fixed.
    std::uint64_t max_beam_extensions = 2000000000ULL;
    // The fixed production ceiling on proposer-owned memory, over the WHOLE predicted peak: rolling
    // vectors, persistent ancestry, pair records and every material temporary and index. 256 MiB is
    // well above what the frozen widths predict (about 18 MB at 1,024 and 72 MB at 4,096 on C4) and
    // far below the process budget, so it refuses a configuration that would be absurd rather than
    // one that is merely large. Checked before the proposal is constructed, so a refusal writes no
    // ranking and no call.
    std::uint64_t max_proposal_payload_bytes = 268435456ULL;   // 256 MiB
    // Merged proposal states retained per block: a resource cap on the grouped-state search. A search
    // that prunes and finds fewer than max_mosaic_candidates off-panel mosaics still calls -- every panel
    // diplotype was evaluated -- and reports mosaic_target_reached false and candidate_scope CONDITIONAL.
    // FROZEN at 256, the smallest plausible width: on all five frozen C4 plans the truth is retained at
    // every layer and its proposal rank is unchanged across 256, 512, 1,024, 2,048 and 4,096 (ranks
    // 6, 6, 1, 8, 87), so width is not the binding constraint there. Counted in allele-prefix PAIRS.
    // This is measured on five synthetic plans at one locus; it is not evidence that 256 suffices
    // everywhere.
    std::size_t proposal_state_limit = 256;
    // Panel pairs times fragments: the predicted diploid-combination work of the protected panel baseline.
    // A baseline over this is refused whole, never scored partially.
    std::uint64_t max_panel_pair_fragment_evaluations = 20000000000ULL;
    std::uint64_t max_piece_hits = 2000000000ULL;
    std::uint64_t max_proposed_starts_per_haplotype = 50000000ULL;
    std::uint64_t max_physical_origins = 500000000ULL;
    std::uint64_t max_state_classes = 200000000ULL;
    std::uint64_t max_mass_cells = 200000000ULL;
    // Bytes the exact emission cache may spend retaining prefilter refusals; zero disables it and
    // reproduces the uncached path exactly. Accepted fits are retained regardless, being far fewer
    // and far dearer to recompute. Bound into the run signature, like every cap above.
    std::size_t emission_cache_bytes = 0;
    bool quiet = false;
    // Development harness only; not a command-line option. When non-empty (strictly increasing, the
    // last equal to max_mosaic_candidates), the model, markers and proposal are computed once, and each prefix
    // k is then refused, ranked and written under <out_prefix>.K<k> by the same code a standalone run
    // with max_mosaic_candidates = k uses: K only truncates the final proposal ordering. Every candidate a
    // prefix contains is scored; the production refusal rule still applies to each k.
    std::vector<std::size_t> development_candidate_prefixes;
    // Not a command-line option. When positive (and at least as large as every primary candidate set), the same proposal
    // and scores are also ranked with this many mosaics, and read-supported phase is reported READS_AUDITED only
    // if that larger set agrees on the winner, the block genotypes and read-supported phase. Without an audit a
    // read-supported phase is READS_WITHIN_CANDIDATES.
    std::size_t phase_audit_candidates = 0;

    // ---- panel-anchored adaptive completion (the production algorithm; set by a profile) --------
    // When repair_rounds > 0 the call is not made from the proposal's candidate set alone. The panel
    // and grouped-mosaic union is exactly scored, the best exactly scored PANEL diplotype anchors an
    // exact one-block coordinate repair, the union and every candidate the repair scored are pooled,
    // and the pool is exactly reranked. The decision is then made over the pool's numerical top set:
    // one candidate calls normally; several with one unordered genotype are
    // CONDITIONAL_GENOTYPE_EQUIVALENT (unphased, no FASTA); several that disagree anywhere refuse.
    // A search that does not converge -- round cap, an unresolved tie among improving neighbours, a
    // cycle, a resource refusal or an interruption -- emits no genotype. 0 keeps the unrepaired call,
    // which is a development view only.
    std::size_t repair_rounds = 0;
    std::size_t repair_local_pairs = 32;
    std::size_t repair_origin_cache_bytes = 402653184;   // 384 MiB, planned retention
    // One round's neighbourhood above this is refused before any allocation.
    std::uint64_t repair_max_neighbours = 670;
    // CUMULATIVE caps over the whole repair search, a different declared quantity from the per-pass
    // caps above, which bound the union's single scoring pass.
    std::uint64_t repair_max_physical_origins = 4900000000ULL;
    std::uint64_t repair_max_state_classes = 4200000000ULL;
    std::uint64_t repair_max_piece_hits = 14000000000ULL;
    std::uint64_t repair_max_latent_pairs = 14000000000ULL;
    // v2: exact branch-and-pool at an improving-move tie (off by default, so a v1 profile runs unchanged). Every
    // tied neighbour is continued by a warm branch over the production search's scoring context, seeded by the
    // exact block delta of a promoted repair; every candidate joins the pool and the pooled rule decides.
    bool repair_tie_branching = false;
    std::size_t repair_tie_max_branches = 16;
    bool test_reverse_tie_branch_order = false;   // gate hook only: the pooled result must not depend on it
    // Process limits, enforced by a watchdog over the whole run (0: none). Crossing one stops the
    // process at once: committed repair rounds stay as <prefix>.repair_rounds.partial.tsv, the staged
    // graph and every final output are removed, and the status says TIME_LIMIT or RESOURCE_REFUSAL.
    std::uint64_t footprint_limit_bytes = 0;
    double time_limit_seconds = 0.0;
    // Provenance of the profile the settings came from, and of every setting a command-line flag
    // overrode after it. Recorded; never interpreted.
    std::string profile_path;
    std::string profile_sha256;
    std::string profile_name;
    std::vector<std::string> overridden_settings;
    // Development harness only; not a command-line option. Stops the process through the time-limit
    // path after this many committed repair rounds that moved, so the interruption contract -- partial
    // rounds kept, no final output, TIME_LIMIT -- is testable deterministically. 0: never.
    std::size_t test_interrupt_after_repair_rounds = 0;

    // ---- DIAGNOSTIC ONLY (v2/post-hoc): improving-move-tie branching. Not a command-line option, off by
    // default, and never part of a call. When on, a run that stops on an improving-move tie first writes
    // every production output exactly as it would otherwise; only then does it reconstruct the terminal
    // tied round, repair from every tied best neighbour on its own branch, pool every scored candidate
    // and rank the pool with the production rule. Everything it finds goes to <prefix>.diag.* files, and
    // its outcome is DIAGNOSTIC_UNIQUE, DIAGNOSTIC_GENOTYPE_EQUIVALENT, DIAGNOSTIC_DIVERGENT or
    // DIAGNOSTIC_INCOMPLETE -- a conditional diagnostic over a pruned search, not a certified call.
    bool diagnostic_tie_branching = false;
    std::size_t diagnostic_max_branches = 16;          // branches started, nested ties included
    double diagnostic_time_limit_seconds = 3600.0;     // wall spent in the diagnostic, checked per round
    std::uint64_t diagnostic_footprint_limit_bytes = 3758096384ULL;   // checked per round
    std::uint64_t diagnostic_max_total_physical_origins = 19600000000ULL;   // over every branch
    bool diagnostic_reverse_branch_order = false;      // gate hook: the result must not depend on it
    // MEASUREMENT ONLY: per search (production, then each branch) its rounds and origin-demand ledger, to
    // price what a warm cache shared across branches would avoid. Changes no score and no decision.
    bool diagnostic_reuse_trace = false;
    // DIAGNOSTIC branch-speed pilot. Warm branches share the production search's scoring context: each
    // tied neighbour is created by the exact block delta of a promoted repair from its parent incumbent's
    // census (pinned, never mutated), instead of a fresh context placing the start in full. The digest
    // memo keys sequence digests by allele path -- immutable content under the fixed model, never a
    // slot or a buffer -- and the check recomputes every memo hit and counts disagreements.
    bool diagnostic_warm_branches = false;
    bool diagnostic_digest_memo = false;
    bool diagnostic_digest_memo_check = false;
    std::uint64_t diagnostic_branch_request_order_seed = 0;   // gate hook: permutes branch requests
    bool diagnostic_pool_dump = false;                        // gate hook: the whole pool, bit-exact
};

// A run that cannot complete emits its diagnostics and says so. It never presents a partial ranking
// as a result.
enum class GenotypeRunStatus {
    Complete,
    ModelIncomplete,
    CandidateIncomplete,
    SpellingIncomplete,
    RescoreIncomplete,
    // Scoring finished, but no single whole-locus diplotype could be separated from its rivals by more
    // than floating-point error, and every candidate that could not be separated agrees on at least one
    // block. Those blocks are called by consensus; blocks where the tied candidates disagree are
    // no-called individually. Deliberately NOT Complete: no unique whole-locus winner was certified, no
    // phased haplotype is emitted, and phase is UNRESOLVED. A caller must be able to tell this from a
    // certified winner without reading the status file.
    BlockConsensus,
    // The pooled top set holds several candidates, numerically tied on total score, that carry the
    // same unordered allele pair at every block: the block genotype is emitted, phase is UNRESOLVED,
    // the equivalence set is reported and no haplotype FASTA is written.
    ConditionalGenotypeEquivalent,
    // The repair search stopped without converging (round cap, a tie among improving neighbours, a
    // cycle, or an incumbent that did not reproduce). The trajectory is kept; no genotype.
    RepairIncomplete,
    // Written by the watchdog as the process stops; a run never returns these.
    TimeLimit,
    ResourceRefusal,
};

const char* genotype_run_status_name(GenotypeRunStatus status);

// The phase claim after an audit. A read-supported phase within the primary candidates becomes READS_AUDITED
// only when the audit set was ranked to completion, agrees on the winner and on every block genotype, and is
// itself read-supported; with the same call but a phase tie there it is PRIOR_ONLY; otherwise UNRESOLVED. Any
// other primary claim is unchanged. Exposed for its gate.
std::string audited_phase_support(const std::string& within_candidates, bool audit_complete, bool same_winner,
                                  bool same_block_genotypes, const std::string& audit_phase_support);

// Per-block genotype confidence within a retained candidate set. Candidate 0 is the winner; scores are total
// log scores and error_bounds their floating-point bounds; genotypes[c][b] is candidate c's unordered allele
// pair at block b. The posterior is the set mass of candidates sharing the winner's genotype there; the quality
// is -10 log10 of the rest, summed in log space so a posterior indistinguishable from one keeps a finite GQ;
// numerically_separated says every candidate with another genotype is below the winner beyond error.
// Exposed for its gate.
struct BlockGenotypeConfidence {
    double posterior = 1.0;
    double quality = 0.0;
    bool numerically_separated = true;
};
std::vector<BlockGenotypeConfidence> block_genotype_confidence(
    const std::vector<double>& scores,
    const std::vector<double>& error_bounds,
    const std::vector<std::vector<std::pair<std::size_t, std::size_t>>>& genotypes);

// Exposed for its gate. The block-consensus rule, over the candidates that could not be separated from
// one another: a block is called when every tied candidate carries the same unordered allele pair there,
// and no-called when they disagree. Deliberately a pure function of the tied genotypes alone -- it reads
// no score, no prior and no order, so reversing the candidate list cannot change its answer, and the
// rule can be tested without contriving two candidates whose scores coincide to within 1e-9.
struct BlockConsensusCall {
    bool agreed = false;
    std::size_t allele1 = 0;
    std::size_t allele2 = 0;
};
// `tied_genotypes[c][b]` is candidate c's unordered allele pair at block b, first <= second. An empty
// set agrees on nothing: with no candidates there is no consensus to report.
std::vector<BlockConsensusCall> block_consensus(
    const std::vector<std::vector<std::pair<std::size_t, std::size_t>>>& tied_genotypes);

GenotypeRunStatus run_genotype(const GenotypeCommandOptions& options);

// Reads a genotype profile: one `key<TAB>value` per line, '#' comment lines. Every key sets one
// option; an unknown key, a repeated key or an unparsable value refuses. Records the profile's path,
// name and SHA-256 in the options so every run states what it was configured by.
void apply_genotype_profile(const std::string& path, GenotypeCommandOptions& options);

int run_genotype_command(const std::vector<std::string>& args);

} // namespace panvar
