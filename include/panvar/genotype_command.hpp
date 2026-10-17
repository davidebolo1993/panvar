#pragma once

#include "panvar/genotype_model.hpp"

#include <cstddef>
#include <cstdint>
#include <iosfwd>
#include <limits>
#include <optional>
#include <map>
#include <string>
#include <unordered_map>
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
    // Bytes of pending tie-branch seeds kept as census-holding slots (built by the exact block delta when their tie is
    // recorded, admitted in canonical tie order while they fit, released when their branch runs); every other seed is
    // materialised just in time. 0 (the default): no seed is kept, the census-free frontier.
    std::uint64_t repair_tie_seed_cache_bytes = 0;
    bool test_reverse_tie_branch_order = false;   // gate hook only: the pooled result must not depend on it
    // Representation-conflict quality guard (an annotation; it never changes a call, a score or a phase).
    // A called block is LOWQUAL_REPRESENTATION_CONFLICT when the marker (depth) channel prefers another
    // available pair at that block, and substituting that pair into the call shows the likelihood's exposure
    // (depth) term preferring it too while its fragment term carries the call -- each by at least the
    // declared threshold. Truth-free and locus-free; see representation_conflict_rule_v1. Off by default: it is
    // bound explicitly by the v3 profiles, so the v1 and v2 profiles run as they were accepted.
    bool representation_guard = false;
    // Exact union placement reuse (performance only; every score is the reference's bit for bit). 0: off.
    // The cap bounds cached placements CONSUMED over the union pass (0: none); crossing it refuses.
    std::size_t union_placement_cache_bytes = 0;
    std::uint64_t max_union_cached_placements_consumed = 0;
    // Exact union origin reuse, on top of placement reuse (performance only; every score is the reference's bit for
    // bit). 0: off. When on, both consumption caps must be positive: served origins and served state aggregates
    // are not executed work and are bounded apart, and crossing either refuses with no ranking.
    std::size_t union_origin_cache_bytes = 0;
    std::uint64_t max_union_cached_origins_consumed = 0;
    std::uint64_t max_union_cached_states_consumed = 0;
    // Exact union prefilter reuse (performance only; every score and logical counter is the reference's). 0: off;
    // when on, the consumption cap must be positive.
    std::size_t union_prefilter_cache_bytes = 0;
    std::uint64_t max_union_cached_prefilters_consumed = 0;
    // Scratch for the union's out-of-place start sorter, over all workers (performance only; the list is the same).
    // 0: in place always.
    std::size_t start_sort_scratch_bytes = 0;
    // I2-A (EXPERIMENTAL, profile keys only, never public settings): the backend of the initial one-shot union pass.
    // "reference" (the default) is score_diplotypes, unchanged. "indexed_fine_experimental" scores the same pass
    // through the I1.4 content index (score_diplotypes_indexed_union) with the same scores and the same refusal
    // semantics; indexed_union_store_bytes (required positive then) bounds its placement store. Only the union pass
    // changes: the adaptive repair (PreparedScoringContext) and the dosage stage keep their own scorers.
    std::string union_scoring_backend = "reference";
    std::size_t indexed_union_store_bytes = 0;
    // I2-C: the indexed union's placement-store layout: "exact" (default; capacity-charged, exact-sized rows) or
    // "legacy_i2a" (the I2-A layout, for comparison only). Valid only with the indexed backend.
    std::string indexed_union_store_layout = "exact";
    // I4-B (EXPERIMENTAL, profile key only): the one-shot union store's evidence representation: "reconstructed"
    // (default; coarse interior rows copied from their fine pieces) or "fine_views_experimental" (coarse units served as
    // views over the fine units' shared rows). Valid only with the indexed backend and the exact store layout.
    std::string indexed_union_evidence_layout = "reconstructed";
    // I2-C (measurement key): report the detailed stage attribution block (seconds_detail_*). Timers always run.
    bool timing_detail = false;
    // I2-D (EXPERIMENTAL, profile key only): seed the adaptive search's anchor pair -- scored columns and censuses -- in
    // one placement pass. Off (default): score, then retain each census by a second full placement. Setting the key
    // explicitly (either value) reports adaptive_fused_seed and adaptive_full_placements in the status.
    bool adaptive_fused_seed = false;
    // I2-E (EXPERIMENTAL, profile key only): exact reuse of marker-proposal contributions in build_proposal_emissions,
    // with this per-block table budget in bytes (allocated capacity). 0 (default): off, every term evaluated directly.
    // When set, the status reports the marker_proposal_reuse_* block. The dosage stage's emissions are not affected.
    std::size_t marker_proposal_reuse_table_bytes = 0;
    // I2-F (EXPERIMENTAL, profile key only): candidate combination uses the run's library terms prepared once per
    // combination call (FragmentScoringOptions::prepared_library_terms) in the initial union pass and the adaptive search;
    // the dosage stage keeps the scalar path. Setting the key (either value) reports the combination_* block.
    bool combination_prepared_library_terms = false;
    // I2-G (MEASUREMENT key, never for timing runs): sample every N-th verification request into the verify_profile_*
    // status block (union pass and whole run). 0 (default): off.
    std::uint32_t verification_profile_period = 0;
    // I2-I (MEASUREMENT key, never for timing runs): the verification-request census (verification_census.hpp), sampling
    // 1 in N read contents into <prefix>.verification_census.tsv; 64 MiB of census storage at most. 0 (default): off.
    std::uint32_t verification_census_period = 0;
    // I2-H (EXPERIMENTAL, profile key only): `gap_recursion_backend scaled_reference|fixed_geometry_experimental`, the
    // emitting gapped recursion's implementation for the run's scoring (FragmentScoringOptions::gap_recursion_backend):
    // the union pass, the adaptive search and dosage. Exact, bit for bit; setting the key reports it in the status.
    bool gap_recursion_fixed_geometry = false;
    // I2-K (EXPERIMENTAL, profile key only): `adaptive_emission_handoff 0|1` -- the indexed union pass's emission results
    // handed to the adaptive search's scoring context (FragmentScoringOptions::emission_handoff_out/_in). Exact; setting
    // the key reports the handoff in the status.
    bool adaptive_emission_handoff = false;
    // I2-B (EXPERIMENTAL, profile keys only, never public settings): how the adaptive search's window enumeration builds
    // its placement rows. "reference" (the default) scans each enumerated window whole, unchanged. "fine_experimental"
    // builds the same rows from a round-local, per-worker fine store (PreparedScoringContext::
    // enable_fine_window_discovery): adaptive_fine_store_bytes bounds its allocated capacity and
    // adaptive_fine_max_served_placements caps what it serves over the search (both required positive then). Setting
    // the key explicitly, either value, adds the adaptive_window_* report block to the status.
    std::string adaptive_window_discovery = "reference";
    std::size_t adaptive_fine_store_bytes = 0;
    std::uint64_t adaptive_fine_max_served_placements = 0;
    // Per-lookup emission-cache timing (two clock reads a lookup). Off: no clock is read and lookup time is reported
    // as not measured.
    bool emission_lookup_timing = false;
    // Gate hooks only, never settings: score every competitor in the live context even when the pool holds
    // it (the outside-pool path), and replace the declared threshold (NaN: kRepresentationBoundNats).
    bool test_representation_guard_force_context = false;
    // Gate hook only, never a setting: keep the pre-fix branch lifetime (every parent incumbent pinned for the whole
    // branching search, no release at a branch's end), so a gate can compare the two lifetimes on one binary.
    bool test_retain_branch_parents = false;
    // Gate hook only, never a setting: build every pending branch seed eagerly at its tie, as a census-holding
    // slot pair derived by block delta from the parent incumbent (the pre-census-free frontier). Production keeps
    // pending seeds as compact descriptors and materialises a branch's start just before the branch runs.
    // test_retain_branch_parents implies it.
    bool test_eager_branch_seeds = false;
    // Gate hook only: the branch with this id declares a start score one ulp off, so its just-in-time
    // reproduction must fail closed (-1: off).
    long test_perturb_branch_seed = -1;
    double test_representation_guard_bound_nats = std::numeric_limits<double>::quiet_NaN();
    // Gate hook: where the marker's best pair is the called pair, take the best other pair as the competitor.
    bool test_representation_guard_runner_up = false;
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
    // Where every setting a preset, profile or flag set came from ("preset short-read", "profile <path>",
    // "flag --kmer-size"); a key not listed took its built-in default. Written to <prefix>.settings.tsv.
    std::map<std::string, std::string> setting_sources;
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
    // MEASUREMENT ONLY: <prefix>.union_haplotype_work.tsv, the union pass's placement work per unique haplotype
    // (panel or mosaic), written whether the pass completes or refuses. Changes no score, order or counter.
    bool diagnostic_union_work_trace = false;
    // DIAGNOSTIC ONLY: stop after the union pass is scored and ranked. Writes <prefix>.union_diagnostic_ranking.tsv
    // and a status UNION_DIAGNOSTIC_ONLY; no repair, no genotype artefact, never a call.
    bool diagnostic_union_only = false;
    // DIAGNOSTIC ONLY: after a completed union pass, replay its seeding path single-threaded over this many evenly
    // spaced panel and as many mosaic unique haplotypes, timing each phase and exact alternatives in isolation;
    // writes <prefix>.union_seeding_profile.tsv. 0: off. Changes no score, order or counter.
    std::size_t diagnostic_union_seeding_profile = 0;
    // DIAGNOSTIC ONLY: after the union is spelled, take the exact census of its direct verification stream
    // (census_union_verification) with this scratch directory, instead of scoring; writes
    // <prefix>.union_verification_census.tsv. Empty: off. A union-only diagnostic, never a call.
    std::string diagnostic_union_verification_census;
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
    // docs/modules/genotype.md, Dosage modes: "off" (default), "audit" or "repair". Audit runs the
    // dosage stage after every existing output is committed and writes only <prefix>.dosage_*; it never changes a call,
    // a score, a status value or an existing file. Repair also applies COMPLETE_DIPLOTYPE_REPAIR blocks to the reported
    // genotype, keeping the caller's outputs as <prefix>.original.*. `--dosage-audit shadow` is audit. Last, so
    // positional fixtures stay aligned.
    std::string dosage_mode = "off";
    // SHADOW ONLY: the allocation observability test inside the shadow audit (off until its cost is reviewed).
    bool dosage_allocation_test = false;
    // I2-J (EXPERIMENTAL, profile key only): `dosage_census_backend reference|interned_experimental`, the dosage block
    // census's implementation (dosage::CensusBackend). Exact; setting the key records it with the census timers.
    bool dosage_census_interned = false;
    // Diagnostic only (docs/algorithms/genotype.md, stage 3): after the run, estimate the library from the reads and the
    // model and write <prefix>.library_estimate.tsv beside the declared library. Never used by the call. Off by default.
    bool library_diagnostic = false;
    // library_sensitivity_audit_v2 (docs/algorithms/genotype.md, stage 9): under --library auto, the call is
    // followed at 0.9x and 1.1x the estimated fragment-start intensity and a genotype that changes is annotated. On by
    // default; the profile key library_sensitivity_audit 0 turns it off (measurement only). No effect on --library FILE.
    bool library_sensitivity_audit = true;
    // --library was omitted on the command line, which means --library auto (recorded as library_invocation).
    bool library_omitted = false;
    // TEST HOOK (gates only): the audit's endpoint 1 (0.9x) or 2 (1.1x) fails before any work, to gate that an
    // incomplete endpoint is reported NOT_EVALUATED and changes nothing. 0: off.
    std::size_t test_library_sensitivity_fail_view = 0;
    // The technology preset applied before the profile and the flags (docs/modules/genotype.md, Required inputs), and
    // the keys it set. Empty: none. A preset never names a locus.
    std::string preset_name;
    bool preset_explicit = false;   // false: the command line's default preset (an omitted --preset means short-read)
    std::vector<std::string> preset_keys;
};

// What the dosage stage and the locus table need from a run: the prepared model and the sorted graph's node sequences.
struct GenotypeDosageCapture {
    std::optional<PreparedGenotypeModel> model;
    std::unordered_map<std::string, std::string> node_sequences;
    // The library the run scored with: the declared one, or the one --library auto built.
    std::optional<LibraryModel> library;
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
    // DIAGNOSTIC ONLY (diagnostic_union_only): the union pass was scored and ranked and the run stopped there --
    // no repair, no genotype, no haplotype. Never a call.
    UnionDiagnosticOnly,
    // --library auto could not validate its estimate (too few decisive pairs, several insert modes, or a contaminated
    // pair distribution). Nothing is scored; the message names --library FILE.
    LibraryRefused,
};

const char* genotype_run_status_name(GenotypeRunStatus status);

// The combined per-block quality (block_quality_rule_v1), truth-free, from two facts the call already reports:
// the representation filter and the genotype support. PASS only for representation PASS with READS support;
// otherwise LOW_QUALITY, with every applying reason in canonical order -- representation_conflict (the filter is
// LOWQUAL_REPRESENTATION_CONFLICT), representation_not_evaluated (any other non-PASS filter: fail closed),
// no_read_support (support other than READS), library_sensitive (library_sensitivity_audit_v1 found this block's
// unordered genotype different at a declared alternative fragment-start intensity) -- joined by ','. It annotates;
// it never changes a call.
// conditional_GQ is not part of it: that quality is conditional on the scored candidate set and is not enough
// for high-quality status.
struct BlockQuality {
    std::string filter;   // PASS or LOW_QUALITY
    std::string reason;   // none, or the reasons in canonical order
};
BlockQuality block_quality(const std::string& representation_filter, const std::string& genotype_support,
                           bool library_sensitive = false);
// The locus filter with a library-sensitivity flag added (library_sensitivity_audit_v2): LOWQUAL_LIBRARY_SENSITIVE (some
// block's genotype differs at an endpoint) or NOT_EVALUATED_LIBRARY_SENSITIVITY (the audit did not complete; a
// locus-level warning, no block is accused). PASS, NA and empty become the flag; LOWQUAL_LIBRARY_SENSITIVE also
// replaces NOT_EVALUATED (a flag wins); anything else keeps its value and gains ";<flag>". Adding a present flag is a
// no-op.
std::string with_library_flag(const std::string& locus_filter, const std::string& flag);

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
// The same run; when capture is non-null it also receives the prepared model and the sorted graph's node sequences.
GenotypeRunStatus run_genotype(const GenotypeCommandOptions& options, GenotypeDosageCapture* capture);

// Reads a genotype profile: one `key<TAB>value` per line, '#' comment lines. Every key sets one
// option; an unknown key, a repeated key or an unparsable value refuses. Records the profile's path,
// name and SHA-256 in the options so every run states what it was configured by.
void apply_genotype_profile(const std::string& path, GenotypeCommandOptions& options);

// One named setting, validated exactly as a profile line is; context names the source in every error.
void apply_genotype_setting(const std::string& context, const std::string& key, const std::string& value,
                            GenotypeCommandOptions& options);
// The public settings (every profile key except test and diagnostic hooks), in a fixed order, and a setting's effective
// value as the run records it. Cross-field constraints (syncmer_s <= kmer_size) are checked once every source has been
// applied, before any work.
const std::vector<std::string>& genotype_public_settings();
std::string genotype_setting_value(const GenotypeCommandOptions& options, const std::string& key);
void validate_genotype_settings(const GenotypeCommandOptions& options);

// The technology presets (docs/modules/genotype.md, Required inputs). "short-read": the production algorithm and the
// short-read technology constants, as (key, value) settings, and the short-read library-model constants. The insert
// is not a preset value (NaN): it comes from --library auto or --library FILE. "long-read" is reserved and refused as
// not implemented; any other name is refused.
std::vector<std::pair<std::string, std::string>> genotype_preset_settings(const std::string& name);
LibraryModel genotype_preset_library(const std::string& name);

// A declared library file (--library FILE): one `key value` per line; blank lines and lines whose first non-space
// character is '#' are skipped, and a value may be followed only by a '#' comment. The six required keys are
// insert_mean_bp, insert_sd_bp, support_sigmas, substitution_error, outlier_weight and fragment_start_intensity;
// marker_background is optional. An unknown key, a repeated key, a missing key, trailing content or a value that is not
// a finite number is refused with the file, line and key named. The gap model is fixed and cannot be declared.
LibraryModel read_library(const std::string& path, std::optional<double>& marker_background);

// Resolves the command line into options: the preset, then every --profile, then the flags, independent of argument
// order. --adaptive-repair and --repair-rounds are resolved together after every flag is read; an explicit
// contradiction between them is refused. Returns the exit status when the arguments asked for help (printed to out),
// and nothing otherwise.
std::optional<int> parse_genotype_arguments(const std::vector<std::string>& args, GenotypeCommandOptions& options,
                                            std::ostream& out);

int run_genotype_command(const std::vector<std::string>& args);

} // namespace panvar
