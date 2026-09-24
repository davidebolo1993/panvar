#pragma once

// Panel-anchored adaptive completion: the exact coordinate repair search, and the decision made over
// the pool of every candidate it and the proposal exactly scored. One implementation, used by
// `panvar genotype` and by the developer probe; neither carries its own copy of the loop.

#include "panvar/fragment_likelihood.hpp"
#include "panvar/genotype_model.hpp"
#include "panvar/marker_model.hpp"

#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace panvar {

// ---- process resources ------------------------------------------------------------------------
// macOS physical footprint: what the kernel charges the process, compressed pages included. Reported
// BESIDE resident size, not instead of it -- resident size falls when the compressor takes pages, so
// under pressure it can read far below the live data it is supposed to bound. Zero where unavailable.
struct FootprintSample {
    std::uint64_t current = 0;
    std::uint64_t peak = 0;
};
FootprintSample physical_footprint();
// CURRENT resident size. ru_maxrss is a high-water mark that can only rise, so it cannot tell a
// plateau from a leak; this can.
std::uint64_t current_rss_bytes();
double peak_rss_gib();
// DIAGNOSTIC ONLY: asks the allocator to return free pages and reports what it released.
std::uint64_t release_free_pages();
// Whole-file replacement by rename, so a reader -- or a run killed mid-write -- only ever sees the
// previous complete version or the new complete version, never a torn one.
void atomic_write(const std::string& path, const std::string& content);
// SHA-256 of a file's bytes, lowercase hex. Provenance only.
std::string sha256_file(const std::string& path);

// ---- the search --------------------------------------------------------------------------------
// OneShot scores every round cold; Context keeps one PreparedScoringContext; Delta additionally
// scores each neighbour as a single-block replacement of the incumbent. Production uses Delta with
// the planned origin cache; the other modes are the references the gates compare it against.
enum class RepairScoringMode { OneShot, Context, Delta };

struct RepairSearchOptions {
    std::size_t local_pairs = 32;        // unordered allele pairs per block, by proposal emission
    std::size_t max_rounds = 11;
    std::size_t first_round = 1;         // a resumed leg continues the round numbering
    RepairScoringMode mode = RepairScoringMode::Delta;
    bool origin_cache = true;            // Delta only
    bool planned_cache = true;           // planned retention; false: first_fit
    std::size_t origin_cache_bytes = 0;
    std::uint64_t origin_cache_records = 0;
    double switch_probability = 0.01;
    // 0: none. A round's neighbourhood above this is refused before any allocation.
    std::uint64_t declared_neighbour_cap = 0;
    // Per-round memory accounts at four fixed points (context modes only).
    bool memory_log = true;
    // Gate hooks, never production: permute the neighbour and repair-request lists every round, and
    // re-score the final neighbourhood N more times without promoting (hold cycles).
    std::uint64_t request_order_seed = 0;
    std::size_t hold_cycles = 0;
    bool pressure_relief = false;
    // DIAGNOSTIC (tie branching): the search builds no context of its own; the caller attaches one that
    // another search owns and keeps alive (attach_context).
    bool external_context = false;
};

// One committed round, with every work and resource column the rounds table reports.
struct RepairRound {
    std::size_t index = 0;
    double start_total = 0.0;
    bool moved = false;
    std::size_t block = 0;
    std::size_t rank = 0;
    int assignment = 0;
    std::size_t allele_x = 0;
    std::size_t allele_y = 0;
    std::vector<std::size_t> first;
    std::vector<std::size_t> second;
    double fragment = 0.0, exposure = 0.0, prior = 0.0, total = 0.0;
    std::uint64_t generated = 0, duplicate = 0, rejected_prior = 0, suppressed = 0, self_skipped = 0;
    std::size_t neighbours = 0, unique_haplotypes = 0, scored_pairs = 0;
    bool complete = false, unique_best = true, order_invariant = true, swap_invariant = true;
    // Whether the incumbent is in the round's numerical top set. A round that ties at the top with
    // the incumbent inside it has converged: no neighbour improves beyond floating-point error.
    bool incumbent_in_top_set = true;
    double separation = 0.0;
    double seconds = 0.0, rss = 0.0;
    std::size_t haps_new = 0, haps_reused = 0;
    std::uint64_t bases_scanned = 0, piece_hits = 0;
    std::uint64_t cache_requests = 0, cache_hits = 0, cache_misses = 0, executed_dps = 0;
    std::uint64_t physical_origins = 0, state_classes = 0;
    double placement_seconds = 0.0, combination_seconds = 0.0;
    std::uint64_t origins_reused = 0, origins_recomputed = 0, window_bases = 0;
    std::uint64_t transient_before = 0, transient_after = 0;
    std::uint64_t full_bytes = 0, window_bytes = 0, key_bytes = 0;
    double seconds_origins = 0.0, seconds_merge = 0.0, seconds_hashing = 0.0;
    double seconds_lookup = 0.0;
    std::uint64_t entry_metadata_bytes = 0, content_arena_bytes = 0;
    std::uint64_t row_constructions = 0, accumulator_growths = 0, origin_growths = 0;
    std::uint64_t oc_hits = 0, oc_misses = 0, oc_stored = 0, oc_entries = 0, oc_bytes = 0;
    std::uint64_t origins_constructed = 0, origins_loaded = 0;
    std::uint64_t fragments_reduced = 0, fragments_copied = 0, merge_working_bytes = 0;
    double seconds_enumeration = 0.0, seconds_filter = 0.0, seconds_aggregate = 0.0;
    double seconds_reduce = 0.0;
    double round_wall_seconds = 0.0;
    std::uint64_t oc_hits_all = 0, oc_misses_all = 0, oc_admitted = 0, oc_admitted_bytes = 0;
    std::uint64_t oc_declined = 0, oc_declined_bytes = 0, origins_constructed_all = 0;
    double seconds_enumeration_all = 0.0, seconds_reenumeration_all = 0.0;
    std::uint64_t cum_piece_hits = 0, cum_physical_origins = 0, cum_state_classes = 0;
    std::uint64_t cum_latent_pairs = 0, mass_cells = 0, records_consumed = 0;
    std::uint64_t census_count = 0, census_bytes = 0, oc_bytes_end = 0;
    std::uint64_t phys_footprint = 0, phys_footprint_peak = 0, rss_current = 0;
    std::uint64_t plan_distinct = 0, plan_known = 0, plan_selected = 0, plan_selected_bytes = 0;
    std::uint64_t plan_retained = 0, plan_evicted = 0, plan_evicted_bytes = 0, oc_not_selected = 0;
    std::uint64_t oc_group_reuses = 0, oc_evicted = 0, oc_evicted_bytes = 0;
};

// A candidate's complete identity: its unordered pair of allele paths.
using DiplotypeIdentity = std::pair<std::vector<std::size_t>, std::vector<std::size_t>>;
DiplotypeIdentity diplotype_identity(const std::vector<std::size_t>& a, const std::vector<std::size_t>& b);

struct RepairHooks {
    // Every candidate a round exactly scored -- the incumbent first, then its neighbours -- with its
    // spelled haplotypes, its score and its sequence-pair prior. Called after the incumbent guard.
    std::function<void(const std::vector<std::size_t>& first, const std::vector<std::size_t>& second,
                       const std::string& first_bases, const std::string& second_bases,
                       const DiplotypeScore& score, double prior)>
        on_scored;
    // After each round is committed to rounds(). May throw (a caller's own resource bound); the
    // exception leaves the search where it stands.
    std::function<void(const RepairRound& round)> on_round;
};

// The declared bound on ONE round's neighbourhood, from the block shape: per block, both assignments
// of up to `local_pairs` unordered allele pairs. Computable before anything is allocated.
std::uint64_t repair_neighbour_bound(const PreparedGenotypeModel& model, std::size_t local_pairs);

class CoordinateRepairSearch {
public:
    CoordinateRepairSearch(const PreparedGenotypeModel& model, const PreparedProposalEmissions& emissions,
                           const std::vector<ReadPair>& fragments, const LibraryModel& library,
                           const FragmentScoringOptions& scoring, const RepairSearchOptions& options);
    ~CoordinateRepairSearch();
    CoordinateRepairSearch(const CoordinateRepairSearch&) = delete;
    CoordinateRepairSearch& operator=(const CoordinateRepairSearch&) = delete;

    // A fresh scoring context (context modes), configured by the options. Resets the seeding count.
    void new_context();
    PreparedScoringContext* context() { return shared_ != nullptr ? shared_ : context_.get(); }

    // ---- DIAGNOSTIC tie branching over one shared context -----------------------------------------
    // Borrows a context another search owns (options.external_context). The context is not owned: it
    // must outlive this search, and its transient round is already open.
    void attach_context(PreparedScoringContext* shared);
    // The census-holding slots of the current incumbent (valid once a delta search has seeded).
    std::size_t incumbent_slot_first() const { return inc_slot_a_; }
    std::size_t incumbent_slot_second() const { return inc_slot_b_; }
    // A neighbour of the current incumbent (one block replaced, side-aligned as the round generated it),
    // created by the same exact block delta a promoted repair uses, census retained. Returns its two
    // slots. A side equal to an incumbent haplotype resolves to that haplotype's slot, as a round does.
    std::pair<std::size_t, std::size_t> seed_neighbour_by_delta(const std::vector<std::size_t>& first,
                                                                const std::vector<std::size_t>& second,
                                                                BlockDeltaWork* work = nullptr);
    // start(), then adopt two census-holding slots of the attached context as the incumbent, so round
    // one deltas from them instead of placing the start in full.
    void start_from_slots(const std::vector<std::size_t>& first, const std::vector<std::size_t>& second,
                          double total, std::size_t slot_first, std::size_t slot_second);
    // Slots every release must keep (other branches' parents and pending seeds). Never mutated here.
    void set_pinned(std::vector<std::size_t> pinned) { pinned_ = std::move(pinned); }
    // Opens the context's transient round once; everything created after it is released between
    // rounds unless retained.
    void open_transient_round();
    void add_seeding_origins(std::uint64_t origins) { seeding_origins_ += origins; }
    std::uint64_t seeding_origins() const { return seeding_origins_; }
    // The memory-row kind of the rounds that follow ("round", or "seed<k>" in a multi-start).
    void set_round_kind(const std::string& kind) { round_kind_ = kind; }

    // Begins a search from a start diplotype whose total score is declared: round one refuses unless
    // it reproduces. Resets the trajectory, the visited set (to the start) and the stop reason.
    void start(const std::vector<std::size_t>& first, const std::vector<std::size_t>& second, double total);
    // Pre-occupies a candidate, as a resumed leg restores the earlier leg's incumbents.
    void add_visited(const std::vector<std::size_t>& first, const std::vector<std::size_t>& second);
    // Runs rounds until a stop. False only when the incumbent did not reproduce; failed_round(),
    // incumbent_expected() and incumbent_scored() then describe it. Throws ResourceExhausted at a cap.
    bool run(const RepairHooks& hooks);

    const std::vector<RepairRound>& rounds() const { return rounds_; }
    const std::string& stop_reason() const { return stop_reason_; }
    const std::vector<std::size_t>& incumbent_first() const { return cur_first_; }
    const std::vector<std::size_t>& incumbent_second() const { return cur_second_; }
    double incumbent_total() const { return cur_total_; }
    std::size_t failed_round() const { return failed_round_; }
    double incumbent_expected() const { return failed_expected_; }
    double incumbent_scored() const { return failed_scored_; }
    std::size_t visited_size() const { return visited_.size(); }
    // The stop is a converged local optimum: no neighbour improves beyond floating-point error. True
    // for no_improvement and no_unvisited_neighbour, and for a round tie that includes the incumbent.
    bool converged() const;
    std::size_t hold_index() const { return hold_index_; }
    double cumulative_seconds() const { return cumulative_seconds_; }
    std::uint64_t cumulative_pairs() const { return cumulative_pairs_; }
    std::uint64_t neighbour_bound() const { return neighbour_bound_; }
    std::size_t last_round() const { return options_.first_round + options_.max_rounds - 1; }
    bool memory_on() const { return memory_on_; }
    const std::string& memory_text() const { return memory_text_; }
    const RepairSearchOptions& options() const { return options_; }
    const FragmentScoringOptions& scoring() const { return scoring_; }
    const std::vector<std::vector<std::pair<std::size_t, std::size_t>>>& local_pairs() const { return local_; }

    // The rounds table (truth-free), one row per round, with every column the gates compare.
    std::string render_rounds(const std::vector<RepairRound>& rounds) const;

private:
    void snapshot(const char* kind, std::size_t index, const char* point, bool result_alive, bool relieve);

    const PreparedGenotypeModel& model_;
    const std::vector<ReadPair>& fragments_;
    const LibraryModel& library_;
    FragmentScoringOptions scoring_;
    RepairSearchOptions options_;
    std::vector<double> switches_;
    std::vector<std::vector<std::pair<std::size_t, std::size_t>>> local_;
    std::uint64_t neighbour_bound_ = 0;
    std::unique_ptr<PreparedScoringContext> context_;
    PreparedScoringContext* shared_ = nullptr;
    std::vector<std::size_t> pinned_;
    bool transient_opened_ = false;
    std::uint64_t seeding_origins_ = 0;
    bool memory_on_ = false;
    std::string memory_text_;
    std::string round_kind_ = "round";
    std::map<std::size_t, std::string> role_of_;
    std::size_t inc_slot_a_ = 0;
    std::size_t inc_slot_b_ = 0;
    bool have_inc_slots_ = false;
    std::vector<RepairRound> rounds_;
    std::string stop_reason_ = "round_cap";
    std::vector<std::size_t> cur_first_;
    std::vector<std::size_t> cur_second_;
    double cur_total_ = 0.0;
    std::set<DiplotypeIdentity> visited_;
    std::size_t hold_remaining_ = 0;
    std::size_t hold_index_ = 0;
    bool holding_ = false;
    double cumulative_seconds_ = 0.0;
    std::uint64_t cumulative_pairs_ = 0;
    std::size_t failed_round_ = 0;
    double failed_expected_ = 0.0;
    double failed_scored_ = 0.0;
};

// ---- the final decision over the pool -----------------------------------------------------------
// Spelled haplotypes are grouped by a 128-bit content digest plus length; allele paths are compared
// exactly. The pool can hold thousands of 200 kb haplotypes, which it must not keep as strings.
struct SequenceDigest {
    std::uint64_t a = 0;
    std::uint64_t b = 0;
    std::size_t length = 0;
    bool operator<(const SequenceDigest& o) const;
    bool operator==(const SequenceDigest& o) const;
};
SequenceDigest sequence_digest(const std::string& bases);

// One exactly scored candidate. `score` is the fragment likelihood; `sequence_prior` the unordered
// sequence-pair prior it is ranked with; `path_prior` its exact allele-path prior, which chooses the
// representative among candidates spelling the same two sequences.
struct PoolCandidate {
    std::vector<std::size_t> first;
    std::vector<std::size_t> second;
    SequenceDigest first_sequence;
    SequenceDigest second_sequence;
    DiplotypeScore score;
    double sequence_prior = 0.0;
    double path_prior = 0.0;
};

// Step 3-5 of the frozen contract, over the numerical top-score equivalence set:
//   Unique             -- one candidate; the normal genotype decision is made from it;
//   GenotypeEquivalent -- several, carrying the same unordered allele pair at every block: the
//                         unphased genotype is emitted, phase is UNRESOLVED, no haplotype FASTA;
//   DistinctTie        -- tied candidates disagree at some block: refused.
// Lower-scoring candidates of the winner's genotype never turn a unique winner into a tie.
enum class TopSetDecision { Unique, GenotypeEquivalent, DistinctTie };
const char* top_set_decision_name(TopSetDecision decision);

struct PooledDecision {
    // The pool in canonical order (by identity); every index below points into it.
    std::vector<PoolCandidate> candidates;
    // Sequence-level candidates: members spelling one unordered sequence pair, in canonical order.
    std::vector<std::vector<std::size_t>> groups;
    std::vector<std::size_t> representative;   // per group: max path prior, ties to canonical order
    std::vector<bool> projection_ambiguous;    // per group: members disagree at some block
    CertifiedRanking ranking;                  // over groups; candidate_id is a group index
    std::uint64_t score_mismatches = 0;        // members of one group that did not score identically
    TopSetDecision decision = TopSetDecision::DistinctTie;
    std::vector<std::size_t> top_groups;       // the numerical top set, in rank order
    std::vector<std::size_t> disagreeing_blocks;   // DistinctTie only
    // Unordered allele pair per block of the decided genotype (Unique or GenotypeEquivalent).
    std::vector<std::pair<std::size_t, std::size_t>> genotype;
    // Separation of the genotype: against the best group carrying any other genotype.
    bool has_other_genotype = false;
    std::size_t best_other_genotype_group = 0;
    double genotype_separation_gap = 0.0;
    bool genotype_separation_certified = true;
    // The best candidate below the top set with the SAME genotype: evidence for PRIOR_ONLY.
    bool has_same_genotype_alternative = false;
    std::size_t best_same_genotype_group = 0;
    double same_genotype_score_gap = 0.0;
    double same_genotype_fragment_gap = 0.0;
    double same_genotype_prior_gap = 0.0;
    bool same_genotype_fragment_identical = false;
    // Phase support, a separate axis: READS_WITHIN_CANDIDATES when no rival is a fragment contender,
    // PRIOR_ONLY when some rival's fragment term is not below the winner's beyond error and only the
    // prior chose it, UNRESOLVED when total scores stay tied within the top set.
    std::string phase_support;
    std::size_t fragment_contenders = 0;
    double fragment_term_gap = 0.0;
    std::vector<bool> block_genotype_read_supported;
    std::vector<double> block_genotype_posterior;
    std::vector<double> block_genotype_quality;
    std::vector<bool> block_numerically_separated;
};

// Ranks the pool exactly and makes the decision. Identities must be unique. Input order never
// changes the result: the pool is put in canonical order first.
PooledDecision decide_pooled_genotype(std::vector<PoolCandidate> pool, std::size_t block_count,
                                      std::size_t fragments, std::size_t mismatch_band);

// ---- phase-switch symmetry of two tied repair branches (diagnostic) ------------------------------
// Two tied improving moves that carry one unordered genotype differ by a LOCAL phase switch at one block:
// the second is the first with its two haplotypes' alleles at that block exchanged. (A global homologue
// swap maps an unordered diplotype onto itself, so two distinct tied candidates are never global swaps of
// each other.) The two branches continued from them are exactly equivalent when, at every depth, the
// candidates branch B scored are the local switches of those branch A scored, with bit-identical scores
// and priors -- and only then is continuing one of them provably enough. Anything else -- a candidate one
// branch scores and the other does not (visited sets differ), a score or prior that differs in a single bit,
// a different number of rounds -- breaks the equivalence at that depth. Convergence alone proves nothing.
// Bit identity is the exact criterion, but two haplotype orders sum the same terms in different orders, so
// an exactly symmetric switch can still differ in the last bits. The second verdict therefore takes the
// certify rule's own equality: counterparts whose totals differ by no more than the sum of their numerical
// error bounds (the bound certify_scored_ranking computes, from `fragment_count`) are indistinguishable, the
// same way the caller's tie is. Neighbourhoods and round counts must still match exactly for either verdict.
struct SymmetryTraceItem {
    std::size_t round = 0;
    std::vector<std::size_t> first;
    std::vector<std::size_t> second;
    DiplotypeScore score;
    double prior = 0.0;
};
struct SymmetryDepth {
    std::size_t round = 0;
    std::size_t a_candidates = 0;
    std::size_t b_candidates = 0;
    std::size_t matched = 0;             // B candidates whose local switch A also scored
    std::size_t bit_identical = 0;       // of those, identical score components and prior, bit for bit
    std::size_t indistinguishable = 0;   // of those, totals within the certify rule's summed error bounds
    double max_total_difference = 0.0;   // over matched candidates
    double max_difference_over_bound = 0.0;   // max |total difference| / summed error bound
};
struct SymmetryReport {
    bool equivalent = false;             // bit-exact
    std::size_t first_break_round = 0;   // 0: none
    std::string first_break;             // what broke there
    bool numerically_equivalent = false;          // under the certify rule's numerical equality
    std::size_t numeric_first_break_round = 0;
    std::string numeric_first_break;
    std::vector<SymmetryDepth> depths;
};
SymmetryReport check_phase_switch_symmetry(const std::vector<SymmetryTraceItem>& a,
                                           const std::vector<SymmetryTraceItem>& b, std::size_t block,
                                           std::size_t fragment_count);

} // namespace panvar
