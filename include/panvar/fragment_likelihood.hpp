#pragma once

#include "panvar/genotype_model.hpp"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <tuple>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace panvar {

// The single definition of read emission, insert probability, valid-FR geometry, exposure, fragment
// contribution and log-addition. Every caller uses these; a private copy of any of them is how the
// accelerated scorer and its oracle start disagreeing about the model they are supposed to share.

// Which mate is read off the haplotype's forward strand. The other mate matches the haplotype's
// reverse strand, which is the same as its own reverse complement matching forward, so every
// comparison below is a forward comparison once the reads are pre-complemented.
enum class MateOrientation { FirstForward, SecondForward };

// A single physical origin of a fragment on one haplotype.
//
// A physical origin is a PLACE, described by the two endpoints of the template and the library
// orientation: the forward mate is anchored at the left endpoint, the reverse mate at the right
// endpoint, and `insert_bp` is their separation. Two distinct repeat copies produce two of these
// even when every statistic agrees, because they are two places the fragment could physically have
// come from.
//
// An ALIGNMENT is not an origin. Once indels are admitted, one origin is explained by several latent
// alignments -- the same left and right endpoints reached with different numbers of haplotype bases
// consumed by each mate -- and those are competing explanations of one physical event, not extra
// physical events. Their masses are therefore SUMMED into `log_latent_mass`, never counted as
// additional multiplicity. Getting that backwards inflates a repeat-free locus with phantom copies.
struct FragmentPlacement {
    // The template's left endpoint: where the forward mate is anchored.
    std::size_t forward_start = 0;
    // Right endpoint minus left endpoint. The right endpoint is forward_start + insert_bp.
    std::size_t insert_bp = 0;
    MateOrientation orientation = MateOrientation::FirstForward;
    // log sum over every latent alignment pair explaining this origin, of the two mates' normalised
    // emissions. This is the whole of what the origin contributes, before the orientation prior and
    // the insert probability -- both of which apply once per origin, because there is one template.
    double log_latent_mass = -std::numeric_limits<double>::infinity();
    // How many latent alignment pairs were summed. One under a gapless model; up to (2d+1)^2 with a
    // displacement band of d. Reported so the pairing cost is measured rather than assumed.
    std::size_t latent_alignments = 0;

    // ---- diagnostics: the single highest-mass latent alignment of this origin -------------------
    // None of these may be used to recompute a score. Under the indel-aware emission a mismatch count
    // does not determine an emission, and one alignment's consumed span does not carry the origin's
    // mass -- attaching the summed mass to the best span's geometry is exactly the defect this
    // structure exists to prevent.
    std::size_t reverse_start = 0;
    std::size_t forward_consumed_bp = 0;
    std::size_t reverse_consumed_bp = 0;
    std::size_t forward_mismatches = 0;
    std::size_t reverse_mismatches = 0;
};

// Fixed-position Hamming emission of `read` against `haplotype` starting at `start`. Version one is
// substitutions only. Returns nullopt when the read runs off the end; a non-ACGT base in either
// sequence refuses, because silently scoring it as a mismatch would turn missing data into evidence.
std::optional<std::size_t> count_mismatches(
    const std::string& haplotype,
    const std::string& read,
    std::size_t start,
    std::size_t mismatch_limit = std::numeric_limits<std::size_t>::max());

double log_emission(
    std::size_t read_length,
    std::size_t mismatches,
    double substitution_error);

// ---------------------------------------------------------------------------------------------
// The declared fragment model: a truncated likelihood.
//
// The mismatch band d is part of the model and there is exactly one definition of it. A mate is
// generated with independent substitutions at rate e, CONDITIONED on carrying at most d of them, so
// its emission is (1-e)^(L-m) (e/3)^m / F(d,L) for m <= d and zero beyond, where F(d,L) is the
// binomial probability of at most d substitutions. Under that model the in-band mass is not an
// approximation of anything: it is the exact likelihood, a fragment with no in-band placement on a
// haplotype has zero mass there, and neither fact depends on which other candidates are scored.
//
// Each physical template is sequenced in one of two library orientations with probability 1/2, the
// insert length follows the library distribution normalised over the library's own support, and
// the fragment-start intensity lambda counts templates per base per chromosome copy summed over both
// orientations. With those three definitions the signal density integrates exactly to the exposure:
// summing M_f,h over every possible observation with fixed read lengths gives E_h. That identity is
// gated directly on a tiny locus.
// ---------------------------------------------------------------------------------------------

// log F(d, L): the probability that a mate of length L carries at most d substitutions.
double log_band_mass(std::size_t read_length, std::size_t mismatch_band, double substitution_error);

// Emission of a mate under the truncated model. -inf beyond the band.
double log_truncated_emission(
    std::size_t read_length,
    std::size_t mismatches,
    double substitution_error,
    std::size_t mismatch_band);

// ---------------------------------------------------------------------------------------------
// Indel-aware emission, version two.
//
// Version one aligned a read at a fixed start and counted substitutions. A single inserted or
// deleted base shifts every downstream base, so a read spanning a one-base indel accumulates
// roughly L/2 mismatches, leaves the declared band, and gets NO placement at all -- it falls to the
// outlier background on the haplotype that is actually correct, while a haplotype that happens to
// match by chance keeps it. That is not a tolerance problem, it is a misspecified emission, and it
// is what cost the two-edit allele 2,871 nats to a sixteen-edit one on NA19240.
//
// The replacement is a banded affine-gap pair-HMM, conditioned on the observed read length:
//
//   states     M consumes one read base and one haplotype base, I one read base, D one haplotype
//              base. No I<->D transition, so a gap has one kind.
//   emissions  M emits the read base with 1-e when it matches and e/3 otherwise; I emits a read
//              base uniformly at 1/4. Both sum to one over the emitted base, and D emits nothing.
//   band       the total edit budget `mismatch_band` counts substitutions AND gap bases together,
//              and net displacement is bounded by `max_gap_bp`.
//
// Because every emitting state's distribution sums to one over the base it emits, summing the
// unnormalised path mass over all 4^L reads of length L leaves exactly the transition mass of the
// admissible alignment SHAPES -- which depends on the read length and the band, never on the
// haplotype or the read. That sum is the normaliser below, so the emission is a proper distribution
// over reads of the observed length, the exposure identity carries over unchanged, and no score can
// depend on which other candidates are present.
//
// At max_gap_bp == 0 this reduces exactly to log_truncated_emission, which is gated as an equality.

// The shape mass log Z: the same recursion with every emission set to one. Independent of the
// sequences' CONTENT, but not of how much haplotype is available: `available_bp` bounds the bases an
// alignment may consume from the start in question, and the normaliser must honour the same bound the
// numerator does. Normalising against an unbounded haplotype instead leaks mass near the ends --
// shapes dropped from the numerator would still divide the denominator, and the emission would sum to
// less than one there, breaking the exposure identity by a small, load-bearing amount.
double log_gap_shape_mass(
    std::size_t read_length,
    std::size_t available_bp,
    std::size_t mismatch_band,
    const LibraryModel& library);

// One terminal consumed-span class of the alignments of a read at a fixed start.
//
// A gapped read does not occupy read_length haplotype bases, and it does not occupy ONE number of
// them either: with a displacement band of d it may end anywhere in read_length +/- d, and each of
// those spans places the mate's far endpoint somewhere different. Because the spans imply different
// fragment geometries and different insert probabilities, the mass of each has to be carried
// separately as far as the geometry -- summing them first and attaching the total to the best span's
// insert length weights real alignment mass with the wrong insert probability, which is what broke
// the exposure identity.
struct GappedComponent {
    // Haplotype bases this class of alignments consumes from `start`.
    std::size_t consumed_bp = 0;
    // Normalised log mass of every admissible alignment ending with that consumed span.
    double log_emission = -std::numeric_limits<double>::infinity();
    // Fewest total edits -- substitutions plus gap bases -- within this class. Diagnostic.
    std::size_t edits = 0;
};

// The components of one fit, held INLINE.
//
// A fit has at most 2*max_gap_bp+1 terminal spans -- seven at the production gap width -- so eight
// slots cover the declared model outright, and anything wider spills to the heap so correctness never
// depends on that bound. The reason it is not simply a std::vector: a cache hit copies components
// into a freshly constructed fit, and with a vector that was one heap allocation per hit, measured at
// 64 allocations over 64 hits and 640 over 640. The arena removed the key allocation; this removes
// the last one on the hit path.
class ComponentList {
public:
    using value_type = GappedComponent;
    static constexpr std::size_t kInline = 8;

    std::size_t size() const noexcept { return size_; }
    bool empty() const noexcept { return size_ == 0; }
    void clear() noexcept {
        size_ = 0;
        heap_.clear();
    }

    GappedComponent* begin() noexcept { return data(); }
    GappedComponent* end() noexcept { return data() + size_; }
    const GappedComponent* begin() const noexcept { return data(); }
    const GappedComponent* end() const noexcept { return data() + size_; }
    GappedComponent& operator[](std::size_t at) noexcept { return data()[at]; }
    const GappedComponent& operator[](std::size_t at) const noexcept { return data()[at]; }
    GappedComponent& back() noexcept { return data()[size_ - 1]; }
    const GappedComponent& back() const noexcept { return data()[size_ - 1]; }

    void push_back(const GappedComponent& value) {
        if (heap_.empty() && size_ < kInline) {
            inline_[size_++] = value;
            return;
        }
        if (heap_.empty()) heap_.assign(inline_, inline_ + size_);
        heap_.push_back(value);
        size_ = heap_.size();
    }

    template <typename Iterator>
    void assign(Iterator first, Iterator last) {
        clear();
        for (; first != last; ++first) push_back(*first);
    }

    void resize(std::size_t count) {
        while (size_ > count) {
            if (!heap_.empty()) heap_.pop_back();
            --size_;
        }
        while (size_ < count) push_back(GappedComponent{});
    }

    void reserve(std::size_t) noexcept {}   // inline storage needs none; the spill handles the rest

private:
    GappedComponent* data() noexcept { return heap_.empty() ? inline_ : heap_.data(); }
    const GappedComponent* data() const noexcept { return heap_.empty() ? inline_ : heap_.data(); }

    GappedComponent inline_[kInline]{};
    std::vector<GappedComponent> heap_;
    std::size_t size_ = 0;
};

// Every admissible alignment of one read at one start, grouped by what it consumes.
struct GappedFit {
    // One entry per terminal consumed span carrying mass, ascending. At most 2*max_gap_bp+1 of them.
    ComponentList components;
    // log sum of the components: the whole read's normalised emission at this start, which is what
    // must sum to one over every read of this length. The components share ONE normaliser -- the
    // shape mass of the read length -- so summing them reproduces this exactly. Normalising each
    // span against its own shape mass would make every span a distribution in its own right and
    // their sum several times one.
    double log_emission = -std::numeric_limits<double>::infinity();
    // The highest-mass component's span and edit count. DIAGNOSTIC ONLY: it is not where the
    // emission lives, and using it as the geometry of the whole emission is the defect above.
    std::size_t consumed_bp = 0;
    std::size_t edits = 0;
};

// Sums over every admissible alignment of the read at this start -- it is a likelihood, not a best
// path -- and returns that sum decomposed by consumed haplotype span, so each part can be given the
// fragment geometry it actually implies. Returns nullopt when the haplotype cannot hold the read at
// this start at all.
// One oriented read form's Myers/Hyyro bit-parallel pattern, built once during fragment preparation
// instead of rebuilt over the read on every verification. The masks are carried as 64-bit halves so
// that the compiler extension they are assembled into stays inside the implementation.
//
// `usable` is false for an empty read or one longer than the word, where the caller falls back to
// edlib. That fallback is the same exact decision by a different route, not an approximation.
struct MyersPattern {
    std::uint64_t peq_lo[4] = {0, 0, 0, 0};
    std::uint64_t peq_hi[4] = {0, 0, 0, 0};
    std::uint32_t length = 0;
    bool usable = false;
};

// Deterministic from the read alone; holds no reference to it.
MyersPattern prepare_myers_pattern(const std::string& read);

// `pattern`, when given, is the read's prepared Myers pattern and must correspond to exactly this
// read; passing nullptr rebuilds it, which is what the reference and test paths do. The decision it
// yields is identical either way -- preparation moves work, not meaning.
std::optional<GappedFit> gapped_emission(
    const std::string& haplotype,
    const std::string& read,
    std::size_t start,
    const LibraryModel& library,
    std::size_t mismatch_band,
    const MyersPattern* pattern = nullptr,
    // True only when the caller has established that the whole haplotype and the whole read are ACGT, so the
    // per-window ACGT refusal cannot fire. Never changes a result; a caller unsure passes false.
    bool acgt_known = false);

// The committed log-space recursion, retained verbatim as the REFERENCE the optimised path is
// certified against. It is never called in production: gapped_emission computes the same quantities
// in scaled probability space, and the oracle suite compares the two component by component on
// random and adversarial reads. Keeping the reference in the library rather than snapshotting its
// numbers means the comparison keeps working when either side is touched again.
std::optional<GappedFit> gapped_emission_reference(
    const std::string& haplotype,
    const std::string& read,
    std::size_t start,
    const LibraryModel& library,
    std::size_t mismatch_band);

// I2-H: which implementation computes the EMITTING gapped recursion. The scaled recursion is the reference and
// the control; the experimental backend must agree with it bit for bit, and falls back to it for any
// configuration it does not support. The normaliser (the non-emitting pass) always uses the reference.
enum class GapRecursionBackend { scaled_reference, fixed_geometry_experimental };

// I2-H (benchmark and gates only; not a production entry point).
namespace gapped_dp_test {
// The emitting recursion alone: the raw (unnormalised) total in log space, and its consumed-span components
// when `components` is given. `available` is how many haplotype bases may be consumed from `start`.
double recursion(const std::string& haplotype, const std::string& read, std::size_t start, std::size_t available,
                 std::size_t band, const LibraryModel& library, ComponentList* components, GapRecursionBackend backend);
// The complete emission operation for a request whose prefilter decision is already known to be feasible:
// recursion, components, normaliser lookup, normalisation, and the emission cache's key, lookup and store
// whenever the cache is armed. The haplotype window and read must be ACGT (not re-checked).
std::optional<GappedFit> emission_known_feasible(const std::string& haplotype, const std::string& read,
                                                 std::size_t start, const LibraryModel& library, std::size_t band,
                                                 GapRecursionBackend backend);
// The public emission (prechecks, refusal of non-ACGT input, prefilter, cache) with the backend chosen.
std::optional<GappedFit> emission(const std::string& haplotype, const std::string& read, std::size_t start,
                                  const LibraryModel& library, std::size_t band, GapRecursionBackend backend);
// Whether the experimental backend has a fixed-geometry recursion for this band and gap allowance (otherwise it
// computes with the reference), and how many fixed-geometry recursions this thread has executed.
bool fixed_geometry_supported(std::size_t band, const LibraryModel& library);
std::uint64_t fixed_geometry_executions();
// The exact feasibility prefilter, as gapped_emission applies it.
bool feasible(const std::string& haplotype, const std::string& read, std::size_t start, const LibraryModel& library,
              std::size_t band);
}  // namespace gapped_dp_test

// Work performed by the emission on THIS thread since the counters were last taken. Reading them
// resets them, so a caller accumulates its own share without coordinating with any other thread.
struct GappedEmissionCounters {
    // Placements where the recursion actually ran.
    std::uint64_t dp_evaluations = 0;
    // Placements the exact feasibility prefilter proved to have zero mass, so the recursion did not.
    std::uint64_t prefilter_rejections = 0;
    // Purely observational work counts, drained from the same thread-local rather than incremented
    // on shared atomics. At 1.19 billion prefilters a contended fetch_add is itself a measurable
    // share of the run, and none of these decides anything: they are merged deterministically in
    // thread-index order after placement, so their totals do not depend on the schedule.
    std::uint64_t executed_prefilters = 0;
    std::uint64_t dp_cache_requests = 0;
    std::uint64_t dp_cache_hits = 0;
    std::uint64_t dp_cache_misses = 0;
    std::uint64_t executed_dps = 0;
};


GappedEmissionCounters take_gapped_emission_counters();

// A census of the verification stream, for deciding whether exact reuse is worth building at all.
// Off by default; when off it costs one relaxed bool test per emission.
//
// An emission KEY is exactly (oriented read bytes, the target window the recursion can touch). Three
// things deliberately do NOT enter it:
//
//   * the haplotype identity or the start offset -- two windows with the same bytes give the same
//     answer, which is the whole point, and is what lets reuse cross haplotypes rather than only
//     repeat within one;
//   * the available length at a boundary -- the window's LENGTH already carries it. The recursion
//     consumes at most read_length + max_gap_bp bases, so every start with at least that much
//     remaining yields a full-length window and an identical computation, while a start nearer the
//     end yields a shorter window and is a genuinely different key;
//   * direction -- reverse-role forms are already reversed against a reversed haplotype before the
//     emission sees them, so the convention is normalised by the sequence itself.
//
// Band and library parameters are fixed within a scoring run and so are not part of the key.
struct EmissionCensus {
    std::uint64_t requests = 0;                    // total verification requests
    std::uint64_t distinct = 0;                    // distinct exact keys
    std::uint64_t prefilter_refusal_requests = 0;  // answered by the prefilter, no recursion
    std::uint64_t prefilter_refusal_distinct = 0;
    std::uint64_t dp_requests = 0;                 // requests that executed the recursion
    std::uint64_t dp_distinct = 0;
    std::uint64_t key_bytes = 0;                   // predicted cache payload: keys ...
    std::uint64_t value_bytes = 0;                 // ... and the values they would hold
};

// Census hashing is 128-bit. A collision would MERGE two distinct keys and therefore UNDER-count
// distinctness, which biases the census against building the cache rather than for it.
void begin_emission_census();
EmissionCensus end_emission_census();

// ---------------------------------------------------------------------------------------------
// Exact reuse of identical emission calculations.
//
// The census measured 81x repetition among DP-evaluated requests and 55x among prefilter refusals
// on NA19240 2x, because the panel's 129 haplotypes share most local windows at homologous
// positions. Reuse is therefore CROSS-HAPLOTYPE and long range: a key recurs about once per
// haplotype pass. That rules out per-thread caches (each worker owns a disjoint slice of the
// haplotypes, so it would see a fraction of the repetition) and rules out eviction under a cyclic
// access pattern, which is LRU's worst case.
//
// The two streams are budgeted differently because the census says they differ by two orders of
// magnitude in both size and value: at 30x the DP set projects to ~1.0M keys costing ~26us each to
// recompute, while the refusal set projects to ~8.7M keys costing ~1.6us each. So DP results are
// always retained, and refusals are retained while the declared budget allows. Exhaustion falls
// back to uncached exact computation -- never eviction -- which changes nothing observable because
// every cached value is an exact function of its key.
// The cache holds ONLY what the recursion produced: accepted fits and refusals discovered after the
// recursion ran. Prefilter refusals are deliberately not retained. Measured at 30x, retaining them
// cost 733,497,742 declined insertions, drove the hit rate from 97.3% down to 37.8%, and spent the
// budget on the cheap stream -- a prefilter refusal costs ~1.6us to recompute against ~26us for a
// recursion, so holding one in the same hash-table representation as the other is a poor trade at
// any budget. The prefilter therefore runs BEFORE a key is formed, and a refusal returns without
// ever touching the table.
struct EmissionCacheCounters {
    std::uint64_t executed_prefilters = 0;   // prefilters actually run (every request pays one)
    std::uint64_t dp_cache_requests = 0;     // requests that survived the prefilter and were looked up
    std::uint64_t dp_cache_hits = 0;
    std::uint64_t dp_cache_misses = 0;
    std::uint64_t executed_dps = 0;          // recursions actually run
    std::uint64_t payload_bytes = 0;         // key and value bytes retained
    std::uint64_t entries = 0;
    std::uint64_t entries_not_retained = 0;  // entries the declared budget declined to hold
    // Bytes charged against the declared budget, which governs REFUSAL entries only -- accepted fits
    // are retained regardless. payload_bytes mixes the two, so it cannot be held to the budget; this
    // can, and is what an exact budget boundary has to be measured against.
    std::uint64_t budgeted_bytes_used = 0;
    // Where the retained bytes actually are. The table holds fixed-size entries carrying arena
    // offsets, and the content arena holds the read and window bytes each entry compares against, so
    // these two are the whole of what a retained entry costs besides its components.
    std::uint64_t entry_metadata_bytes = 0;
    std::uint64_t content_arena_bytes = 0;
};

// Explicitly owned scorer memory, counted rather than inferred. RSS is not evidence about which
// buffer grew, and the number of haplotypes does not by itself imply more simultaneously live
// placement vectors: each worker handles one haplotype at a time, so haplotype count drives
// allocator CHURN, not residency. These accounts separate the two.
struct ScorerMemoryAccounts {
    std::uint64_t mass_matrix_bytes = 0;
    std::uint64_t haplotype_sequence_bytes = 0;      // forward and reversed copies
    std::uint64_t form_sequence_bytes = 0;
    std::uint64_t piece_index_bytes = 0;
    std::uint64_t emission_cache_table_bytes = 0;
    std::uint64_t emission_cache_arena_bytes = 0;
    // Summed over workers, at each worker's own peak. Rows are placements actually held; capacity is
    // what the reused buffers retain between haplotypes.
    std::uint64_t placement_rows_peak = 0;
    std::uint64_t placement_capacity_bytes_peak = 0;
    std::uint64_t proposed_starts_peak = 0;
    std::uint64_t proposed_starts_capacity_bytes_peak = 0;
};

// `budget_bytes` bounds what refusal entries may occupy; zero disables the cache entirely, which is
// the default and reproduces the uncached implementation exactly. The cache lives for one scoring
// pass: keys reference the caller's form and haplotype buffers, which must outlive it.
// TEST ONLY. Forces every emission key's hash to zero, so every distinct key collides into one
// shard and byte equality is the sole thing separating them. A real 128-bit collision cannot be
// manufactured to order; collapsing the hash exercises the same code path, harder, on every key
// rather than on one contrived pair. Never set on a production path.
void set_emission_cache_hash_collapse(bool collapse);
// TEST ONLY: off, prepared forms hash their whole read into every emission key as the public path does (the
// reference); on (the default), they resume from the hash state carried with their bytes. The keys are equal.
void set_emission_key_prefix(bool enabled);
// Gates only: the resumed key equals key_view_of's for this read and window, in every field.
bool emission_key_prefix_matches(const std::string& read, const std::string& window);

void begin_emission_cache(std::size_t budget_bytes);

// The live counters, without disarming or discarding anything. Separate from end_emission_cache so
// that a caller can record what the cache did while the cache is still in use -- the teardown
// belongs to whichever scope owns the buffers the keys point into, and doing both in one call meant
// whoever read first silently disarmed it for everyone else.
EmissionCacheCounters read_emission_cache_counters();

EmissionCacheCounters end_emission_cache();

// I2-K (EXPERIMENTAL, off by default): the union pass's emission results, handed to the adaptive search.
//
// Not a second cache: the process-global emission cache's own shards (entry tables, key-byte slabs, component arenas)
// are MOVED into this object when union discovery would otherwise discard them, together with the cache's budget
// accounting and the scoring context the results were computed under. The caller owns it; destroying it frees them.
// The adaptive PreparedScoringContext adopts it -- moves the shards back -- only under an exactly compatible scoring
// context (library error and gap model, mismatch band, recursion backend) and only when the inherited budgeted bytes
// fit its own unchanged cache budget; otherwise it discards it and starts empty. Adoption consumes it. Entries keep
// their exact identity (read and window bytes, compared byte for byte); placement rows, origins and decisions are never
// part of it. Inherited entries are marked so that hits on them are counted.
class EmissionCacheHandoff {
public:
    EmissionCacheHandoff();
    ~EmissionCacheHandoff();
    EmissionCacheHandoff(const EmissionCacheHandoff&) = delete;
    EmissionCacheHandoff& operator=(const EmissionCacheHandoff&) = delete;
    bool holds() const;
    std::uint64_t entries() const;
    std::uint64_t allocated_bytes() const;   // table capacity + slabs + component-arena capacity
    std::uint64_t budgeted_bytes() const;
    void clear();
    struct Contents;
    std::unique_ptr<Contents> contents;
};

// What the most recent handoff did (process-wide; reset with reset_emission_handoff_report, at the start of a run).
struct EmissionHandoffReport {
    bool detached = false, offered = false, adopted = false;
    std::string refusal = "none";   // none, not_offered, context, budget, cache_disabled
    std::uint64_t entries = 0, allocated_bytes = 0, budgeted_bytes = 0;
    std::uint64_t inherited_hits = 0;   // cache hits served by inherited entries (live while the adopting context lives)
};
EmissionHandoffReport emission_handoff_report();
void reset_emission_handoff_report();

// I3-A.1: process-global storage a finished sample can leave behind (capacity of process-global containers, plus the
// calling thread's own scratch), as allocated capacity in bytes. Worker threads' thread-locals end with the workers.
struct RetainedProcessStorage {
    std::uint64_t emission_cache_table_bytes = 0;      // entry-table capacity, kept by end_emission_cache's clear()
    std::uint64_t emission_cache_component_bytes = 0;  // component-arena capacity, likewise
    std::uint64_t emission_cache_slab_bytes = 0;       // key-byte slabs (freed by end_emission_cache; 0 between samples)
    std::uint64_t emission_census_bytes = 0;           // the measurement-only emission census table (buckets and nodes)
    std::uint64_t calling_thread_scratch_bytes = 0;    // this thread's gapped-recursion scratch
    bool emission_cache_armed = false;
    std::uint64_t total() const {
        return emission_cache_table_bytes + emission_cache_component_bytes + emission_cache_slab_bytes +
               emission_census_bytes + calling_thread_scratch_bytes;
    }
};
RetainedProcessStorage retained_process_storage();
// Releases that storage (returns what it held before). Refuses -- std::logic_error -- while the emission cache is armed
// or pinned by a live PreparedScoringContext: it is only for the gap between finished samples.
RetainedProcessStorage release_sample_retained_storage();

// What the most recently completed scoring pass's cache did. end_emission_cache() tears the table
// down before any status writer runs, so the figures are retained here rather than recovered by
// threading the scoring result through every caller that wants to report them.
EmissionCacheCounters last_emission_cache_counters();

// Runs the previous window enumeration alongside the current one and refuses if they disagree on any
// placement, origin, counter or ordering. Off by default and never on a production path: it doubles
// the work. The reference is kept for exactly this reason, as gap_recursion_logspace is.
void set_window_enumeration_self_check(bool enabled);

// Stage 6 certification. When on, every origin-path delta also runs a fresh stage-5 enumeration and
// the stage-5 std::map merge, and refuses on any divergence in retained origins, merged StateKeys,
// multiplicities, reduced masses, or a promoted census. Off by default: it more than doubles the work.
void set_origin_cache_self_check(bool enabled);
// TEST ONLY. Collapses every origin-cache key hash to one value, so byte comparison alone separates
// entries. Never set on a production path.
void set_origin_cache_hash_collapse(bool collapse);

// Seconds spent inside cache lookups, accumulated per thread and folded in as placement drains its
// counters. Reported separately from origin construction because the two answer different questions:
// whether the cache is cheap to consult, and whether the work it saves was worth consulting it for.
// Per-lookup timing is off unless enabled (diagnostics); emission_cache_lookup_seconds() is then not a measurement.
void set_emission_lookup_timing(bool enabled);
bool emission_lookup_timing_enabled();
double emission_cache_lookup_seconds();
void reset_emission_cache_lookup_seconds();

// The memory accounts of the most recently completed scoring pass, retained for the same reason the
// cache counters are: write_status is reached from eleven call sites and carries no scoring result,
// so threading a diagnostic through all of them would be a worse change than keeping the last
// pass's figures here.
ScorerMemoryAccounts last_scorer_memory_accounts();

// The same emission, anchored at the read's RIGHT endpoint instead of at its start.
//
// An FR pair is sequenced from both ends of ONE template: the forward mate is read rightward from
// the left endpoint, the reverse mate leftward from the right endpoint. Each mate's emission has to
// be a proper distribution given ITS OWN anchor, and the reverse mate's anchor is the right
// endpoint -- exactly as FragmentPlacement above describes a physical origin.
//
// This is not gapped_emission with the arguments moved about. The recursion is DIRECTIONAL: an
// alignment may end in a gap but may not begin in one, so a path and its reversal do not carry the
// same mass (measured at 3.8 to 11.2 nats apart), and the normaliser is conditioned on the sequence
// available in the direction the mate is actually read.
//
// Left-anchoring the reverse mate is what broke the exposure identity. Several starts share one
// right endpoint -- rs = R - consumed, one per consumed span -- and each was normalised against a
// different amount of remaining haplotype, so gathering one consumed-span class from each summed to
// MORE than one wherever the haplotype end bound the normaliser. The excess appeared only near the
// right boundary and walked inward as max_gap_bp widened, which is precisely what was measured.
// Anchored here, the identity holds exactly at every (mismatch_band, max_gap_bp) combination, with
// no precondition relating the two, and at max_gap_bp 0 it is identical to gapped_emission.
//
// `end` is one past the last haplotype base the mate may consume; a component's consumed span is
// measured leftward from it, so that alignment occupies [end - consumed_bp, end).
std::optional<GappedFit> gapped_emission_ending_at(
    const std::string& haplotype,
    const std::string& read,
    std::size_t end,
    const LibraryModel& library,
    std::size_t mismatch_band);

// log(1/2): each template is read in one of two library orientations.
inline constexpr double kLogMateOrientationPrior = -0.69314718055994530942;

// The template span from the first outer boundary to the second. Overlapping mates are physical
// states, so the floor is the longer mate, never the sum. Returns nullopt when the geometry is not a
// valid FR pair or the span falls outside this fragment's declared support.
std::optional<std::size_t> fr_insert_length(
    std::size_t forward_start,
    std::size_t forward_length,
    std::size_t reverse_start,
    std::size_t reverse_length,
    const InsertSupport& support);

// The same admission rule stated on the template's two ENDPOINTS, which is how an origin is named
// once the reverse mate is anchored at its right endpoint. Equivalent to fr_insert_length with
// reverse_start = reverse_end - reverse_consumed_bp, but it cannot underflow when a gapped mate
// consumes more than its endpoint offset, and it does not ask the caller to reconstruct a start it
// no longer has.
std::optional<std::size_t> fr_insert_from_endpoints(
    std::size_t forward_start,
    std::size_t forward_consumed_bp,
    std::size_t reverse_end,
    std::size_t reverse_consumed_bp,
    const InsertSupport& support);

// Expected number of templates on one chromosome copy of this haplotype that can produce a fragment
// whose longer mate is `physical_floor_bp`, under the library insert distribution:
// sum over i >= floor of P(i) * (L - i + 1)^+. Analytic, and summed over both orientations because
// lambda is.
double haplotype_exposure(
    std::size_t haplotype_length,
    const DiscreteInsertDistribution& insert,
    std::size_t physical_floor_bp = 0);

// log((1-eta) * lambda * mass + eta * P_background(f)). `log_mass` is the log of the summed
// placement mass over BOTH homologues, so a homozygote has already counted two chromosome copies.
double log_fragment_contribution(
    double log_mass,
    double log_background,
    const LibraryModel& library);

// The candidate-dependent part of the integrated intensity, which is what the -E term of the score
// must be. The per-fragment intensity is the superposition of a signal process,
// (1-eta) * lambda * (M_a + M_b), and a background process, eta * P_background(f). Integrated over
// every possible observation the signal gives (1-eta) * lambda * (E_a + E_b) and the background gives
// eta per read-length class, the same for every candidate and therefore dropped. The (1-eta) factor
// is NOT candidate-independent: scoring with lambda * (E_a + E_b) instead ranks haplotypes of different
// lengths against a density that does not integrate to its own exposure.
double expected_signal_fragments(double exposure_sum, const LibraryModel& library);

// ---------------------------------------------------------------------------------------------
// Reference implementation.
// ---------------------------------------------------------------------------------------------

// Every physical origin of `fragment` on `haplotype`, found by trying every start pair in both
// library orientations. This is the oracle: it is O(L^2) per fragment and exists to certify the
// accelerated scorer on small fixtures, never to run on a real locus.
//
// `mismatch_limit` is the declared band of the truncated model. Placements beyond it have zero
// emission under that model, so they are absent from the enumeration by definition.
std::vector<FragmentPlacement> enumerate_placements_exhaustive(
    const std::string& haplotype,
    const ReadPair& fragment,
    const InsertSupport& support,
    const LibraryModel& library,
    std::size_t mismatch_band);

// Summed log mass of a fragment on one chromosome copy under the truncated model: sum over physical
// origins of (1/2) * P_library(insert) * the origin's summed latent alignment mass. `insert` must be
// the library distribution. Returns -inf when no placement is accepted.
//
// The orientation prior and the insert probability multiply the origin ONCE, because one origin is
// one template with one insert length. The latent alignments inside it were already summed when it
// was built, at the geometry each of them implied.
// There is deliberately no band parameter. Under the substitution-only model the emission was
// recoverable from a mismatch count, so the same placement list could be re-normalised at scoring
// time by passing a different band. A gapped alignment is not summarised by a count: the band selects
// AND normalises the placements when they are enumerated, and this can only sum what they carry.
// Keeping a band argument here would let a caller ask for a renormalisation that silently does not
// happen, which is worse than not offering it.
double log_placement_mass(
    const std::vector<FragmentPlacement>& placements,
    const ReadPair& fragment,
    const DiscreteInsertDistribution& insert,
    const LibraryModel& library);

// ---------------------------------------------------------------------------------------------
// Accelerated scorer.
// ---------------------------------------------------------------------------------------------

// K diplotype candidates contain at most 2K haploid mosaics, and mosaics that differ in their allele
// path can still spell identical bytes. Scoring the unique sequences once and combining two columns
// per candidate is what keeps the cost linear in unique haplotypes rather than quadratic in
// candidates.
struct InternedHaplotypes {
    std::vector<std::string> sequences;
    // Parallel to the input: which unique sequence each spelled haplotype resolved to.
    std::vector<std::size_t> index_of;
    std::size_t distinct_allele_paths = 0;
};

InternedHaplotypes intern_haplotypes(
    const std::vector<MosaicPath>& allele_paths,
    const std::vector<std::string>& sequences);

struct HaplotypeCandidate {
    std::size_t first = 0;    // index into InternedHaplotypes::sequences
    std::size_t second = 0;
};

// Thrown when a declared resource cap is reached. Scoring stops and returns nothing: a ranking built
// from the placements found before the cap would silently favour whichever haplotypes were scanned
// first, which is a partial ranking presented as a result.
class ResourceExhausted : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

// MEASUREMENT ONLY: one unique haplotype's placement work in a one-shot scoring pass -- deltas of the worker's own
// counters around that haplotype, and the time of its two phases. Filled only when requested; never read by a score.
struct HaplotypeWork {
    bool recorded = false;
    std::uint64_t haplotype_bp = 0;
    std::uint64_t piece_hits = 0, proposed_starts = 0, unique_starts = 0, full_read_verifications = 0;
    std::uint64_t prefilter_rejections = 0, dp_evaluations = 0, executed_dps = 0, accepted_mate_placements = 0;
    std::uint64_t physical_origins = 0, statistical_state_classes = 0;
    // Phase 2 (union origin reuse): origins and state aggregates served from the cache, and the StateKey classes
    // of the merged census -- executed plus served reconcile with the reference's physical_origins and classes.
    std::uint64_t reused_physical_origins = 0, reused_state_aggregates = 0, merged_state_classes = 0;
    double seconds_seeding_and_emission = 0.0;   // exact seeds, prefilter, gapped emission, both strands
    double seconds_pairing_and_reduction = 0.0;  // mate pairing, physical origins, StateKey reduction
    double seconds_cached_merge = 0.0;           // of which merging served state aggregates (phase 2)
    // Of seconds_seeding_and_emission, both strands: piece scan and start generation, start ordering and
    // deduplication, per-start verification (window checks, prefilter, emission lookup, DP, row storage), phase-1
    // replay, and the final row sort. The rest is the reverse rows' coordinate turn.
    double seconds_seed_scan = 0.0, seconds_start_sort = 0.0, seconds_verification = 0.0;
    double seconds_replay = 0.0, seconds_row_sort = 0.0;
};

// Exact union placement reuse (one-shot passes). A start whose whole dependency window -- the bases its seeds
// can come from, [start - max_gap_bp, start + R + max_gap_bp) with R the longest read form, and the bases its
// verification and emission read -- lies inside one segment of a haplotype has placements that depend on that
// segment's bytes alone (the emission normaliser is saturated there). Such placements are computed once per
// distinct (strand, segment bytes) occurring at least twice, stored in local coordinates, and replayed at every
// occurrence; every other start is placed directly as before. Rows are sorted by (start, consumed) at the end
// either way, so the result is the reference's, bit for bit. The segmentation is a hint: any segmentation is
// exact, a poor one only reuses less. What was reused is reported here.
struct UnionPlacementCacheReport {
    bool active = false;
    std::size_t read_span_bp = 0;                 // R
    std::size_t max_gap_bp = 0;                   // the seed and verification margin
    std::uint64_t budget_bytes = 0;
    std::uint64_t max_placements_consumed = 0;    // 0: none
    std::uint64_t segment_occurrences = 0;        // (haplotype, strand, segment) with a non-empty interior zone
    std::uint64_t distinct_segments = 0;
    std::uint64_t eligible_segments = 0;          // distinct segments occurring at least twice
    std::uint64_t admitted_segments = 0;
    std::uint64_t declined_segments = 0;          // did not fit the byte budget: placed directly instead
    std::uint64_t evictions = 0;                  // nothing admitted is ever evicted
    std::uint64_t hits = 0;                       // occurrences served from the cache
    std::uint64_t misses = 0;                     // occurrences placed directly
    std::uint64_t zone_starts_requested = 0;      // start positions inside interior zones
    std::uint64_t zone_starts_served = 0;         // of those, served from the cache
    std::uint64_t placements_enumerated = 0;      // accepted placements computed while filling (whole segments)
    std::uint64_t placements_cached = 0;          // interior placements retained
    std::uint64_t placements_replayed = 0;        // consumed: replayed into haplotypes (held to the cap)
    std::uint64_t payload_bytes = 0;              // retained placements
    std::uint64_t index_bytes = 0;                // segment keys held while filling
    std::uint64_t peak_bytes = 0;                 // payload + index + the batch being filled, at its peak
    double seconds_fill = 0.0;
};
// TEST ONLY: collapse every union segment key's hash, so byte comparison alone separates segments.
void set_union_placement_cache_hash_collapse(bool collapse);

// Exact union ORIGIN reuse (phase 2; one-shot passes, on top of phase 1). A physical origin is (fragment, arm,
// forward start s, insert I), and its latent mass sums the forward placements starting at s against the reverse
// placements ending at e = s + I -- nothing else. When s lies in the forward interior zone of one segment
// occurrence, [o + g, o + L - R - g], and e in the same occurrence's reverse interior zone, [o + R + g, o + L - g],
// every one of those placements is a phase-1 replayed placement of that segment, so the origin -- insert,
// orientation, latent mass, summed in the same order -- is a function of the segment's bytes alone. Such origins
// are constructed once per distinct segment (both strands admitted by phase 1) as sorted per-fragment StateKey
// multiplicities, and merged at every occurrence; each occurrence contributes its own copies, so multiplicity is
// preserved. Every other pair -- across a boundary, near the locus ends, in a short or declined segment -- is
// paired directly. The partition is by (s, e) alone, so it is disjoint and exhaustive; the reduction runs over the
// merged sorted census in ascending StateKey order, exactly as the reference's map does.
//
// Resource semantics: physical_origins, statistical_state_classes, latent_alignment_pairs and valid_fr_states
// count EXECUTED work (direct pairing plus the one-time construction), and max_physical_origins and
// max_state_classes bind those. Served origins and aggregates are counted apart (reused_* in the counters) and
// bound by max_union_cached_origins_consumed and max_union_cached_states_consumed; crossing either refuses before
// any ranking. The reference's totals are recovered exactly as executed - constructed + reused.
struct UnionOriginCacheReport {
    bool active = false;
    std::uint64_t budget_bytes = 0;
    std::uint64_t max_origins_consumed = 0;
    std::uint64_t max_states_consumed = 0;
    std::uint64_t segment_occurrences = 0;        // forward segment occurrences with a non-empty interior zone
    std::uint64_t eligible_segments = 0;          // distinct segments with both strands admitted by phase 1
    std::uint64_t admitted_segments = 0;
    std::uint64_t declined_segments = 0;          // did not fit the byte budget: paired directly
    std::uint64_t evictions = 0;                  // nothing admitted is ever evicted
    std::uint64_t hits = 0;                       // occurrences served
    std::uint64_t misses = 0;                     // occurrences paired directly
    std::uint64_t origins_constructed = 0;        // physical origins built once, while filling (every eligible segment)
    std::uint64_t state_aggregates_constructed = 0;
    std::uint64_t latent_pairs_constructed = 0;
    std::uint64_t valid_states_constructed = 0;
    std::uint64_t origin_records_held = 0;        // aggregates retained (admitted segments)
    std::uint64_t origins_held = 0;               // physical origins those stand for, per occurrence
    std::uint64_t origins_consumed = 0;           // served over the pass (held to max_origins_consumed)
    std::uint64_t states_consumed = 0;            // aggregates served over the pass (held to max_states_consumed)
    std::uint64_t payload_bytes = 0;
    std::uint64_t index_bytes = 0;                // the per-haplotype zone tables
    std::uint64_t peak_bytes = 0;                 // payload + index + the batch being filled, at its peak
    double seconds_fill = 0.0;                    // cached-origin construction (wall)
};

// DIAGNOSTIC ONLY: a deterministic single-thread replay of the union pass's seeding path over sampled unique
// haplotypes, after the production phase-1/2 fills, splitting its cost into phases by repeating each phase in
// isolation over the same inputs; and timing exact alternatives against it -- each alternative's output is
// compared with the production path's and any difference throws. Nanoseconds, summed over the sampled
// haplotypes' two strands. Nothing here enters a score.
struct UnionSeedingPhases {
    std::uint64_t haplotypes = 0;
    std::uint64_t production_ns = 0;          // the production collect_form_placements calls, as one
    std::uint64_t scan_ns = 0;                // seed scan and index lookups alone (no start generated)
    std::uint64_t scan_and_generate_ns = 0;   // the same scan generating every proposed start
    std::uint64_t start_sort_ns = 0;          // std::sort + std::unique of the proposals
    std::uint64_t start_sort_bucketed_ns = 0; // the counting-by-form alternative (exact)
    std::uint64_t window_check_ns = 0;        // eligibility, window and ACGT validation per unique start
    std::uint64_t prefilter_ns = 0;           // the exact feasibility prefilter per unique start
    std::uint64_t prefilter_early_ns = 0;     // the early-refusing kernel (exact)
    std::uint64_t emission_call_ns = 0;       // gapped_emission per unique start, every DP already cached
    std::uint64_t dp_timed_ns = 0;            // gap recursions timed directly (a prefix of the survivors)
    std::uint64_t storage_ns = 0;             // FormPlacement construction and row storage
    std::uint64_t replay_ns = 0;              // phase-1 cached placements appended
    std::uint64_t row_sort_ns = 0;
    std::uint64_t piece_hits = 0, proposals = 0, unique_starts = 0, prefilters = 0, prefilter_rejections = 0;
    std::uint64_t dp_requests = 0, executed_dps = 0, dp_timed = 0, placements = 0, replayed = 0;
    // The emission key path over the prefilter survivors, every key already cached: key construction by the full
    // hash and by the carried read-side state (identical h0/h1, checked), the production lookup, and its parts
    // without the clock reads, without the shard lock, without the component copy, and the probe alone.
    std::uint64_t survivors = 0;
    std::uint64_t key_full_ns = 0, key_prefix_ns = 0, lookup_ns = 0, lookup_unclocked_ns = 0;
    std::uint64_t lookup_nolock_ns = 0, lookup_nocopy_ns = 0, probe_only_ns = 0, clock_pairs_ns = 0;
    // The early-refusing prefilter over the same eligible starts, streamed through 2, 4 and 8 lanes per read form
    // (decisions checked identical to the one-start kernel's).
    std::uint64_t batched_starts = 0, lanes2_ns = 0, lanes4_ns = 0, lanes8_ns = 0;
};
struct UnionSeedingProfile {
    UnionSeedingPhases panel, mosaic;
    std::vector<std::size_t> sampled;              // unique haplotype indices, in replay order
    std::vector<std::uint64_t> sampled_hashes;     // FNV-1a of each sampled haplotype's bytes
    std::uint64_t bucketed_sort_checked = 0;       // proposal lists compared, all identical
    std::uint64_t early_prefilter_checked = 0;     // prefilter decisions compared, all identical
    std::uint64_t lane_decisions_checked = 0;      // multi-lane decisions compared, all identical
    std::uint64_t prefix_keys_checked = 0;         // prefix-state keys compared, all identical
    // The placement-cache fill, per eligible segment in canonical order: wall seconds and length, with the batch
    // width it ran at -- enough to price another assignment of the same work.
    std::vector<std::pair<double, std::size_t>> fill_segments;
    std::size_t fill_workers = 0;
    double fill_seconds = 0.0;
};

// Exact union prefilter reuse (one-shot passes). The feasibility prefilter's decision depends on the read's bytes and
// the window's bytes alone (the window's length carries the end condition), and the union's direct stream repeats
// the same pairs many times. Read contents are interned at preparation and window contents per pass -- both by
// complete byte equality -- and each (read, window) decision is remembered in a fixed table keyed by the two exact
// IDs; a remembered decision stands in for the kernel. Scores, placements and every logical counter are unchanged;
// executed_prefilters counts only kernels that ran. Hits are held to their own consumption cap.
struct UnionPrefilterCacheReport {
    bool active = false;
    std::uint64_t budget_bytes = 0;
    std::uint64_t max_consumed = 0;
    std::size_t read_length = 0, window_length = 0;
    std::uint64_t window_contents = 0;       // distinct window contents interned
    std::uint64_t window_index_bytes = 0;    // their copies and tables
    std::uint64_t positions_indexed = 0;     // eligible starts given a window identity
    std::uint64_t decision_table_bytes = 0;
    std::uint64_t decisions_held = 0;
    std::uint64_t decisions_not_retained = 0;   // the table was full: the kernel ran and nothing was stored
    std::uint64_t hits = 0;                  // consumed: decisions served instead of a kernel
    std::uint64_t misses = 0;                // the kernel ran
    double worker_seconds_index = 0.0;       // window identities, summed over workers
};
void set_union_prefilter_cache_hash_collapse(bool collapse);   // TEST ONLY: one hash for every window

// The start sorter over one pass. With a scratch budget, each of `workers` holds at most budget / workers bytes
// (4 per proposal), kept across haplotypes; a list too large for that share is ordered in place instead -- the
// same list either way -- and counted as a fallback. Without a budget every list is ordered in place.
struct StartSortReport {
    std::uint64_t budget_bytes = 0;
    std::size_t workers = 0;
    std::uint64_t share_bytes = 0;         // the most any worker may hold
    std::uint64_t out_of_place = 0;
    std::uint64_t in_place = 0;
    std::uint64_t fallbacks = 0;           // in place because the share was too small
    std::uint64_t peak_worker_bytes = 0;   // the largest scratch a worker held (never above share_bytes)
};

struct FragmentScoringOptions {
    // The declared band of the truncated model. Fixed for a run: changing it changes the model.
    std::size_t mismatch_band = 6;
    // Pigeonhole pieces shorter than this are too unspecific to index, so the scorer refuses rather
    // than silently degrading into a quadratic scan.
    std::size_t minimum_piece_bp = 10;
    // Operational caps, checked as the work happens rather than estimated beforehand. The defaults
    // leave two orders of magnitude of headroom over C4 at 30x and K=64 (14.8M piece hits, 1.9M
    // proposed starts per haplotype, 0.88M physical origins) while stopping a repeat-rich locus such
    // as LPA before it exhausts memory. Reaching one throws ResourceExhausted.
    //
    // Under parallel placement these are PROCESS-WIDE budgets, not per worker. Two mechanisms, with
    // different jobs, and the distinction matters:
    //
    //   * In flight, max_piece_hits and max_physical_origins are shared across workers and bound the
    //     WORK actually performed. Crossing one cancels every worker promptly rather than letting each
    //     consume a full cap -- which is what a circuit breaker is for, and what a per-worker guard
    //     silently gave up. physical_origins is an exact shared counter; piece_hits is published in
    //     batches of piece_hit_publish_batch because it is the hottest counter in the program (320M
    //     increments on C4 at audit 256), and batching was chosen to avoid unmeasured contention on a
    //     shared cache line. Whether an exact atomic would in fact be affordable here has NOT been
    //     measured. The overshoot that batching admits is bounded and declared: at most
    //     threads * piece_hit_publish_batch beyond the cap, which is 8,192 against a 2,000,000,000
    //     default.
    //
    //   * ACCEPTANCE is decided once, after the counters are merged, on the exact totals. A REFUSAL may
    //     also be raised earlier, by whichever worker observes the shared budget cross the cap -- so it
    //     is not true that every verdict waits for the merge. Both are sound, and the accept-or-refuse
    //     answer still does not depend on the thread count: a published total never exceeds the exact
    //     total, so a run whose exact total is within the cap cannot trip an early refusal, and a run
    //     that trips one would have exceeded the cap on the merged totals as well. A total equal to the
    //     cap accepts, one above refuses, however the work was partitioned.
    //
    // max_state_classes is deliberately only a post-merge acceptance cap. It is counted after a
    // fragment's states already exist, so bounding it in flight would limit neither memory nor work.
    std::uint64_t max_piece_hits = 2000000000ULL;
    // Increments a worker accumulates locally before publishing to the shared piece-hit budget. Larger
    // is cheaper and admits proportionally more overshoot; the bound is threads * this value.
    std::uint64_t piece_hit_publish_batch = 4096;
    std::uint64_t max_proposed_starts_per_haplotype = 50000000ULL;
    std::uint64_t max_physical_origins = 500000000ULL;
    std::uint64_t max_state_classes = 200000000ULL;
    std::uint64_t max_mass_cells = 200000000ULL;
    // Admitting indels multiplies the pairing work by up to (2*max_gap_bp+1)^2 per physical origin,
    // so it gets a declared bound of its own rather than being absorbed into the origin cap. Like
    // max_state_classes this is an acceptance cap checked on the merged totals: the pairs are counted
    // as the geometry is tested, and a running per-worker crossing point would move with the
    // partition. The default holds the same headroom over C4 that the origin cap does.
    std::uint64_t max_latent_alignment_pairs = 2000000000ULL;
    // Bytes the exact emission cache may spend retaining prefilter refusals. Zero -- the default --
    // disables the cache and reproduces the uncached path exactly, so no existing caller changes
    // behaviour by upgrading. DP results are retained regardless of this budget: the census shows
    // them to be both far fewer and far dearer to recompute.
    std::size_t emission_cache_bytes = 0;
    // Diagnostics only; nothing here enters a score or a ranking. When set, the result carries for
    // every (fragment, haplotype) the accepted physical-origin count and the single highest-mass origin,
    // and for the listed fragments every accepted origin.
    bool record_origin_summary = false;
    std::vector<std::size_t> record_origins_for_fragments;
    // MEASUREMENT ONLY (one-shot passes): when set, sized to the unique haplotypes and filled per haplotype. Owned by
    // the caller, so a pass refused at a cap still leaves the haplotypes it finished.
    std::vector<HaplotypeWork>* record_haplotype_work = nullptr;
    // Exact union placement reuse (see UnionPlacementCacheReport). One-shot passes only; off unless both a
    // segmentation and a byte budget are given. The segmentation is per unique haplotype (index order of the
    // interned set): segment lengths summing to the haplotype's length.
    const std::vector<std::vector<std::uint32_t>>* union_segment_lengths = nullptr;
    std::size_t union_placement_cache_bytes = 0;
    // A cap on cached placements CONSUMED over the pass (0: none). Replayed placements are not executed work
    // and are not charged as such, so reuse needs a bound of its own: exceeding it refuses, with no ranking.
    std::uint64_t max_union_cached_placements_consumed = 0;
    // Exact union origin reuse (see UnionOriginCacheReport): off unless phase 1 is on and a byte budget is given.
    // Both consumption caps are then required to be positive -- a served origin is not free.
    std::size_t union_origin_cache_bytes = 0;
    std::uint64_t max_union_cached_origins_consumed = 0;
    std::uint64_t max_union_cached_states_consumed = 0;
    // Exact union prefilter reuse (see UnionPrefilterCacheReport): off unless given a byte budget, which bounds the
    // decision table; the consumption cap bounds decisions served and must then be positive. Needs one read length.
    std::size_t union_prefilter_cache_bytes = 0;
    std::uint64_t max_union_cached_prefilters_consumed = 0;
    // Scratch for the out-of-place start sorter, over all workers (see StartSortReport). 0: in place always.
    std::size_t start_sort_scratch_bytes = 0;
    // Workers for the diploid combination phase. Candidates are partitioned into fixed contiguous
    // ranges and each worker writes only its own complete DiplotypeScore entries, so the result does
    // not depend on the thread count: no cross-thread floating-point reduction exists to reorder.
    // Defaults to serial, so every existing caller keeps its behaviour without being touched.
    std::size_t threads = 1;
    // The library-intensity sensitivity audit (library_sensitivity_audit_v1): further libraries whose scores every
    // diploid combination also computes, from the same mass columns and exposures in the same pass. Nothing is placed
    // again: the library enters only the combination. The run's own scores are unchanged, and each audit score is,
    // bit for bit, what a run declaring that library computes. Empty: no audit views.
    std::vector<LibraryModel> audit_libraries;
    // I2-F (EXPERIMENTAL, off by default): the run's own library's combination constants -- log1p(-w) + log(lambda),
    // log(w), (1 - w) * lambda -- are prepared ONCE per combination call (validated there, through the audit views' own
    // preparation) instead of being re-validated and re-evaluated per (candidate, fragment) term by the scalar
    // log_fragment_contribution / expected_signal_fragments. Each prepared constant is the exact subexpression those
    // functions evaluate, in the same association, so every score component is unchanged bit for bit.
    bool prepared_library_terms = false;
    // I2-H (EXPERIMENTAL, off by default): the emitting gapped recursion's implementation for every placement this
    // run's scoring performs. fixed_geometry_experimental computes the same cells with compile-time geometry and
    // single assignment, bit for bit; geometries it does not support, and the normaliser, use the scaled reference.
    GapRecursionBackend gap_recursion_backend = GapRecursionBackend::scaled_reference;
    // I2-K (EXPERIMENTAL, off by default): the indexed union batch detaches its emission cache into this handoff
    // instead of discarding it; an adaptive PreparedScoringContext adopts it under an exactly compatible context.
    EmissionCacheHandoff* emission_handoff_out = nullptr;
    EmissionCacheHandoff* emission_handoff_in = nullptr;
};

// I2-F: what preparing the main library's terms removed, summed over combination calls since the last reset (counted
// once per call from its candidate and fragment counts, never per term). Process-wide; measurement only.
// I2-G correction: the scalar reference evaluates log1p(-w) and log(lambda) only for terms whose combined mass is finite
// (fragments_with_mass, already computed per candidate), and log(w) for every term when w > 0; it validates once per
// term and once per candidate. Preparation evaluates those logarithms and one validation once per call that combines any
// term. The counts are taken once per call after the combination -- never per term -- and the net is SIGNED: preparation
// can exceed an empty or nearly empty reference.
struct PreparedLibraryTermsTally {
    std::uint64_t combinations = 0;                  // combination calls that used prepared terms
    std::uint64_t terms = 0;                         // (candidate, fragment) terms combined with them
    std::uint64_t hypothetical_scalar_logs = 0;      // logarithms the scalar reference would have evaluated
    std::uint64_t prepared_logs = 0;                 // logarithms preparation evaluated
    std::int64_t net_logs_avoided = 0;               // hypothetical minus prepared (may be negative)
    std::uint64_t hypothetical_scalar_validations = 0;
    std::uint64_t prepared_validations = 0;
    std::int64_t net_validations_avoided = 0;
};
PreparedLibraryTermsTally prepared_library_terms_tally();

void reset_prepared_library_terms_tally();

struct OriginSummary {
    std::uint64_t origins = 0;
    FragmentPlacement best;
    double best_log_term = -std::numeric_limits<double>::infinity();
};

// Measured work, not estimated viability. `candidate_bases_scanned` counts haplotype bases actually
// walked, which is the honest denominator; counting hash lookups instead flatters the index.
struct FragmentScoringCounters {
    std::uint64_t candidate_bases_scanned = 0;
    std::uint64_t piece_hits = 0;
    std::uint64_t proposed_starts = 0;
    std::uint64_t unique_starts = 0;
    std::uint64_t full_read_verifications = 0;
    // Of those verifications, how many reached the pair-HMM and how many the exact prefilter proved
    // impossible first. Their sum is full_read_verifications, and the rejected fraction is the
    // headline number for whether the prefilter is earning its place.
    std::uint64_t dp_evaluations = 0;
    std::uint64_t prefilter_rejections = 0;
    std::uint64_t accepted_mate_placements = 0;
    std::uint64_t valid_fr_states = 0;
    // Latent alignment pairs TESTED: one per (forward consumed span, reverse consumed span) offered
    // to the FR geometry, whether or not it was accepted. This is the pairing work that admitting
    // indels multiplies, which is why it counts attempts rather than successes -- `valid_fr_states`
    // is the accepted subset. Measured rather than assumed negligible: a mate has at most
    // 2*max_gap_bp+1 terminal spans, so one mate pair can present up to (2*max_gap_bp+1)^2 of them --
    // 49 at the production max_gap_bp of 3 -- before the collapse onto endpoints reduces them.
    std::uint64_t latent_alignment_pairs = 0;
    // Fragments whose declared insert support lets the geometry's floor -- insert >= the larger
    // consumed span -- actually bind. Only those need their consumed-span components carried
    // separately through the pairing: for every other fragment the mass at an anchor pair
    // factorises exactly as (sum over forward spans) * (sum over reverse spans), so one
    // whole-emission row per anchor yields the IDENTICAL likelihood at a fraction of the cost.
    // Counted because a fast path nobody can see is one that can diverge unnoticed; the scorer's
    // agreement with the always-componentised oracle is what certifies it.
    std::uint64_t fragments_with_floor_binding = 0;
    std::uint64_t physical_origins = 0;
    std::uint64_t statistical_state_classes = 0;
    // Phase-2 union origin reuse only (zero otherwise): what was SERVED rather than executed. The reference's
    // physical_origins is physical_origins - constructed + reused_physical_origins (likewise the latent pairs and
    // valid states), and its statistical_state_classes is merged_state_classes.
    std::uint64_t reused_physical_origins = 0;
    std::uint64_t reused_state_aggregates = 0;
    std::uint64_t reused_latent_alignment_pairs = 0;
    std::uint64_t reused_valid_fr_states = 0;
    std::uint64_t merged_state_classes = 0;
};

// Exact relative to the declared truncated model. `absolute_term_sum` is the sum of the magnitudes
// of every term added into the score, which is what bounds its floating-point accumulation error.
struct DiplotypeScore {
    double log_likelihood = 0.0;
    double fragment_term = 0.0;
    double exposure_term = 0.0;
    double absolute_term_sum = 0.0;
    std::size_t fragments_with_mass = 0;
};

struct FragmentScoringResult {
    UnionPlacementCacheReport union_placement_cache;
    UnionOriginCacheReport union_origin_cache;
    UnionPrefilterCacheReport union_prefilter_cache;
    StartSortReport start_sort;
    std::vector<DiplotypeScore> scores;     // parallel to the requested candidates
    // [audit library][candidate]: the same candidates under FragmentScoringOptions::audit_libraries.
    std::vector<std::vector<DiplotypeScore>> audit_scores;
    // Row-major [fragment * unique_haplotypes + haplotype]: the summed log placement mass of one
    // fragment on one chromosome copy of one haplotype. Retained because it is what the exhaustive
    // oracle certifies directly; recovering it from a candidate score would compare the two paths
    // through the mixture rather than at the mass itself.
    std::vector<double> log_mass;
    // Per unique haplotype: expected templates on one chromosome copy, weighted by the observed
    // distribution of fragment physical floors.
    std::vector<double> exposure;
    std::size_t unique_haplotypes = 0;
    FragmentScoringCounters counters;
    // Wall clock of the three phases inside one call, so a caller can attribute rescoring instead of
    // reporting it whole. Preparation is per fragment; placement and mass construction is per
    // (fragment, haplotype); diploid combination is per (candidate, fragment) and is the term that
    // grows with the panel. Which of the last two dominates decides how this should be parallelised.
    double seconds_fragment_preparation = 0.0;
    double seconds_placement_and_mass = 0.0;
    double seconds_diploid_combination = 0.0;
    // What was actually used, which is not what was asked for when there are fewer candidates than
    // workers. Reported so a run can never claim parallelism it did not use.
    std::size_t threads_used_combination = 1;
    // Workers that actually ran the placement and mass phase. One when a single worker was asked for,
    // when there are fewer haplotypes than workers, and when full-origin recording is requested --
    // that diagnostic writes a shared map from the innermost loop and deliberately stays serial rather
    // than being made deterministic for a path nothing statistical depends on.
    std::size_t threads_used_placement = 1;
    // candidates x fragments: the work the combination phase performs, independent of thread count.
    std::uint64_t candidate_fragment_evaluations = 0;
    // Work the exact emission cache did, if it was armed. Deliberately NOT part of
    // FragmentScoringCounters: those count logical events and must stay identical between a cached
    // and an uncached run of the same input, which is what makes the two comparable at all.
    EmissionCacheCounters emission_cache;
    // What this pass actually owned, measured at its peak. Reported so a memory question is answered
    // by an account rather than by reasoning backwards from RSS.
    ScorerMemoryAccounts memory;
    // Filled only when requested in FragmentScoringOptions; same row-major layout as log_mass.
    std::vector<OriginSummary> origin_summary;
    // (fragment, haplotype) -> every accepted origin, for record_origins_for_fragments only.
    std::map<std::pair<std::size_t, std::size_t>, std::vector<FragmentPlacement>> recorded_origins;

    double mass_of(std::size_t fragment, std::size_t haplotype) const;
};

// Scores exactly the requested candidates and nothing else. The raw score of a candidate does not
// depend on which other candidates were requested, on their order, or on any proposal score.
// What a refused run actually did, as it stood when a cap was reached. A refusal that reports nothing
// is unfalsifiable -- a breaker that stopped at its cap cannot be told from one that stopped long past
// it -- and reporting only counters loses the rest: which workers ran, and where the time went. Both
// were previously left at their success-path defaults on a refusal, so a refused run claimed one
// placement thread however many it used, and attributed its whole duration to no phase at all.
struct FragmentScoringObservation {
    FragmentScoringCounters counters;
    // Workers that actually ran placement, not the number requested.
    std::size_t threads_used_placement = 1;
    // Zero when the phase never ran, which is the usual case for a placement refusal. Deliberately not
    // defaulted to one: claiming a worker for a phase that did not happen is the defect this replaces.
    std::size_t threads_used_combination = 0;
    // The phases up to the refusal. Preparation is complete by then, placement is however far it got,
    // and combination has usually not started.
    double seconds_fragment_preparation = 0.0;
    double seconds_placement_and_mass = 0.0;
    double seconds_diploid_combination = 0.0;
    UnionPlacementCacheReport union_placement_cache;   // what the reuse did before the refusal
    UnionOriginCacheReport union_origin_cache;
    UnionPrefilterCacheReport union_prefilter_cache;
    StartSortReport start_sort;
};

// A measurement, not an optimisation. Oriented read forms repeat -- four forms per fragment, and at
// depth many fragments carry identical sequence -- so their placements could in principle be computed
// once per distinct sequence and reused. Whether that is worth building depends on the saving in PIECE
// HITS, not in forms: duplicate forms are plausibly the repetitive ones that contribute
// disproportionately many hits, so the form-level collapse is not the answer. This reports both, and
// the exact weighted prediction, so the decision rests on the quantity that actually costs.
struct PieceHitCensus {
    std::size_t total_forms = 0;
    std::size_t unique_forms = 0;
    // Seed entries across every piece-length bucket, and how many survive a collapse to one form per
    // distinct sequence.
    std::uint64_t total_entries = 0;
    std::uint64_t unique_entries = 0;
    // The scorer's own scan, counted rather than placed: sum over scanned positions of the entries
    // matching that seed. `unique` is the same sum counting only the first form of each distinct
    // sequence, which is exactly what a collapse would leave. Their ratio is the predicted saving.
    std::uint64_t predicted_piece_hits_all = 0;
    std::uint64_t predicted_piece_hits_unique = 0;
    // Distribution, so a mean saving cannot hide a single pathological seed or haplotype.
    std::uint64_t max_piece_hits_at_one_seed = 0;
    std::uint64_t max_piece_hits_one_haplotype = 0;
};

// Builds the same forms and piece indexes the scorer builds and scans the same haplotypes, counting
// hits instead of placing them. Changes nothing and is never called from scoring.
PieceHitCensus census_piece_hits(
    const InternedHaplotypes& haplotypes,
    const std::vector<ReadPair>& fragments,
    const FragmentScoringOptions& options);

// `observed_on_failure`, when given, receives the work performed when a cap was reached, so a refused
// run can report what it actually did rather than nothing at all.
// DIAGNOSTIC ONLY: an exact census of the union's direct verification stream (phase-1 zones applied, as the pass
// runs it). Every request is identified by (read content, window content): read contents interned by complete
// byte equality, window contents -- every eligible window of every haplotype strand, its length min(available,
// read + max_gap_bp) carrying the end condition the emission's normaliser depends on -- by complete byte
// equality against a representative, hashes only choosing where to look. Every request's early-refusing
// prefilter decision and the columns it scanned are recorded; every repetition of an identity must reproduce
// both. Requests are aggregated through `scratch_dir` (removed afterwards). Requires one read length.
struct UnionCensusClass {
    std::uint64_t haplotypes = 0, requests = 0, ineligible = 0, accepted = 0, rejected = 0;
    std::uint64_t first_accepted = 0, first_rejected = 0, repeated_accepted = 0, repeated_rejected = 0;
    std::uint64_t columns_accepted = 0, columns_rejected = 0, columns_repeated_rejected = 0;
    std::uint64_t columns_repeated_accepted = 0;
    std::uint64_t prefilter_accepted_ns = 0, prefilter_rejected_ns = 0;   // production kernel, timed by decision
    std::uint64_t identities_first_here = 0;
};
struct UnionVerificationCensus {
    UnionCensusClass panel, mosaic;
    std::size_t read_length = 0, window_length = 0;
    std::uint64_t read_forms = 0, read_contents = 0;
    std::uint64_t window_positions = 0, window_contents = 0, window_table_bytes = 0, window_index_bytes = 0;
    std::uint64_t identities = 0;
    std::map<std::uint64_t, std::uint64_t> multiplicity;   // occurrences -> identities with that many
    // Cumulative distinct window contents after each haplotype of a class, in interned order.
    std::vector<std::uint64_t> accumulation_panel, accumulation_mosaic;
    std::uint64_t fill_executed_prefilters = 0, fill_verifications = 0;
    std::uint64_t decisions_checked = 0;
    std::uint64_t record_bytes = 0;
    double seconds_windows = 0.0, seconds_stream = 0.0, seconds_aggregate = 0.0;
};
UnionVerificationCensus census_union_verification(const InternedHaplotypes& haplotypes,
                                                  const std::vector<ReadPair>& fragments, const LibraryModel& library,
                                                  const FragmentScoringOptions& options,
                                                  const std::vector<char>& is_panel, const std::string& scratch_dir);

// See UnionSeedingProfile. `options` must carry the union segmentation and both reuse budgets as the pass ran.
UnionSeedingProfile profile_union_seeding(const InternedHaplotypes& haplotypes, const std::vector<ReadPair>& fragments,
                                          const LibraryModel& library, const FragmentScoringOptions& options,
                                          const std::vector<std::size_t>& sample, const std::vector<char>& is_panel);

FragmentScoringResult score_diplotypes(
    const InternedHaplotypes& haplotypes,
    const std::vector<HaplotypeCandidate>& candidates,
    const std::vector<ReadPair>& fragments,
    const LibraryModel& library,
    const FragmentScoringOptions& options,
    FragmentScoringObservation* observed_on_failure = nullptr);

// ---------------------------------------------------------------------------------------------
// Certification.
// ---------------------------------------------------------------------------------------------

enum class RescoreStatus { Complete, RescoreIncomplete };

// Every score is exact relative to the declared model, so the only uncertainty left in a comparison
// is floating-point accumulation. `numeric_error_bound` is a conservative first-order bound on it:
// (terms + 16) * DBL_EPSILON * sum |term| for a recursively summed score. It is not a model interval
// and it never widens with the band.
struct CertifiedCandidate {
    std::size_t candidate_id = 0;      // index into the requested candidate list
    double score = 0.0;                // fragment term + exposure term + exact prior
    double fragment_term = 0.0;
    double exposure_term = 0.0;
    double prior = 0.0;
    double numeric_error_bound = 0.0;
    std::size_t fragments_with_mass = 0;
    // exp(score - log sum_c exp(score_c)) over the retained candidates. A probability within the
    // candidate set only: omitted candidates carry mass this does not see, so it is not a posterior.
    double candidate_set_probability = 0.0;
};

struct CertificationOptions {
    // A floor under the numerical tie tolerance, for scores so small that the computed bound is below
    // any meaningful difference.
    double minimum_tie_tolerance_nats = 1e-9;
};

struct CertifiedRanking {
    RescoreStatus status = RescoreStatus::RescoreIncomplete;
    std::vector<CertifiedCandidate> candidates;   // sorted by score, descending
    // One id when a winner is certified; every candidate that cannot be separated from the best when
    // it is not. An undecided ranking reports the set rather than picking from it.
    std::vector<std::size_t> winner_set;
    // min over rivals of (best - its error bound) - (rival + its error bound). Positive exactly when
    // the status is Complete.
    double certified_gap = 0.0;
    // Raw score difference between the best and the runner-up.
    double score_gap_to_runner_up = 0.0;
    double candidate_set_log_weight = 0.0;
    std::size_t mismatch_band = 0;
    FragmentScoringCounters counters;
};

// Scores every requested candidate once under the declared truncated model, adds the exact prior,
// and certifies a winner when its score minus its numerical error bound exceeds every rival's score
// plus theirs. There is no refinement loop: nothing is approximated, so nothing can be deepened.
// A candidate's score depends only on its own two haplotypes, the fragments and the library.
//
// `exact_prior` is added to the matching candidate's score. It must be exact, which is why it is a
// separate argument. The proposal score is not, and may never be, passed here: it comes from the
// same reads as the fragment term, and adding both would count that evidence twice.
// Certification of scores already computed by score_diplotypes: scores[c] is candidate c. A diplotype's
// score does not depend on the other candidates it is scored with, so a subset of one scoring pass is
// ranked exactly as scoring that subset alone would rank it.
CertifiedRanking certify_scored_ranking(
    const std::vector<DiplotypeScore>& scores,
    std::size_t fragment_count,
    std::size_t mismatch_band,
    const FragmentScoringCounters& counters,
    const CertificationOptions& certification,
    const std::vector<double>& exact_prior = {});

CertifiedRanking certify_diplotype_ranking(
    const InternedHaplotypes& haplotypes,
    const std::vector<HaplotypeCandidate>& candidates,
    const std::vector<ReadPair>& fragments,
    const LibraryModel& library,
    const FragmentScoringOptions& options,
    const CertificationOptions& certification,
    const std::vector<double>& exact_prior = {});


// ---------------------------------------------------------------------------------------------
// Persistent scoring context: exact incremental haplotype scoring.
// ---------------------------------------------------------------------------------------------
//
// score_diplotypes prepares everything, scores everything, and throws it all away. A repair search
// calls it once per round over a candidate set that is almost entirely the previous round's, so it
// re-derives the same read forms, the same piece indexes and the same haplotype columns every time.
// This context owns that work for its whole lifetime instead, and scores only what is new.
//
// It is EXACT, not an approximation of the one-shot path: a haplotype's mass column and exposure
// depend on that haplotype, the fragments and the library alone -- never on which other haplotypes
// are present -- so a column computed in batch one is the column batch seven would have computed.
// The gates require the same haplotypes supplied in 1, 2, 3 and 7 batches to agree cell for cell.
//
// What it owns, for its whole lifetime, and why each has to be owned rather than borrowed:
//
//   read forms and piece indexes   rebuilt per call otherwise, and the emission cache's keys point
//                                  into the form bytes.
//   insert supports, backgrounds   per fragment, fixed for the search.
//   haplotype sequences, forward   in std::deque, whose element addresses never move. The cache key
//   and reversed                   holds raw pointers into these; a std::vector would relocate them
//                                  on growth, and only the accident that a 226kb string is not
//                                  short-string-optimised keeps that from being a dangling read
//                                  today. Requirement: do not retain pointers into relocatable
//                                  vectors.
//   emission cache storage         armed once for the context, never between batches, which is the
//                                  only reason a key may outlive the batch that created it.
//   mass columns and exposures     retained so combination can be rerun cheaply, and so a later
//                                  batch reuses them exactly rather than recomputing them.
//   cumulative counters, refusal   caps bound the WHOLE search, not each round separately.
//
// A refusal invalidates the context permanently: every later call throws, and no partial ranking is
// ever produced. That is the same contract score_diplotypes has within one call, extended to the
// search, because a search that continued past a cap would rank whatever it had managed to score.

// What one incremental batch actually did. Per batch, and cumulative, because the caps are
// cumulative and a per-batch figure alone cannot be checked against them.
struct IncrementalBatchReport {
    std::size_t haplotypes_offered = 0;
    std::size_t haplotypes_new = 0;        // sequences not already held, so actually placed
    std::size_t haplotypes_reused = 0;     // already had a column, reused exactly
    std::size_t unique_haplotypes_total = 0;
    // This batch alone: what placing the new haplotypes cost.
    FragmentScoringCounters counters;
    EmissionCacheCounters emission_cache;
    double seconds_placement_and_mass = 0.0;
    std::size_t threads_used_placement = 1;
    ScorerMemoryAccounts memory;
    // The whole context so far, which is what the caps are held against.
    FragmentScoringCounters cumulative_counters;
    EmissionCacheCounters cumulative_emission_cache;
    double cumulative_seconds_placement_and_mass = 0.0;
};

// ---------------------------------------------------------------------------------------------
// Stage 2: exact haploid block-delta scoring.
// ---------------------------------------------------------------------------------------------
//
// Placing a repaired haplotype from scratch rescans the whole locus to learn what one block changed.
// A delta reuses the rest. It is exact, not an approximation, and the reason is structural: a
// fragment's mass on one haplotype is a reduction over a StateKey -> multiplicity census, and an
// origin the replacement cannot touch keeps its key unchanged.
//
// "Cannot touch" is decided on the half-open template interval [forward_start, forward_start +
// insert_bp): an origin is affected exactly when that interval intersects the replaced block. An
// origin ending precisely at the block start, or starting precisely at its end, is NOT affected.
//
// One subtlety decides whether reuse is sound at the ends, and it is not about which bases an origin
// reads. The emission normaliser log_gap_shape_mass clamps available_bp at read_length + max_gap_bp
// and is bit-identical beyond it, so an unaffected origin keeps its normaliser only while it is
// saturated. A length-changing replacement shifts the forward mate's available length (L - start)
// for origins BEFORE the block, and the reverse mate's (its right endpoint) for origins AFTER it.
// Both are saturated unless the origin lies within read_length + max_gap_bp of the corresponding
// end. So the affected set widens to a whole end exactly when the replaced block sits within that
// margin of it -- which is why repairing a first or last block costs more than repairing an interior
// one, and is measured rather than assumed.

// One accepted physical origin, carrying only what the reduction and the affected test need. The
// diagnostic fields of FragmentPlacement are deliberately absent: they are not inputs to a score,
// and at hundreds of thousands of origins per haplotype they would dominate the census.
struct HaploidOrigin {
    std::uint32_t forward_start = 0;   // template left endpoint
    std::uint32_t insert_bp = 0;       // template span; the right endpoint is start + insert_bp
    std::uint32_t orientation = 0;     // which library arm produced it
    double log_latent_mass = 0.0;      // summed over this origin's latent alignments
};

// The complete origin census of ONE haplotype. Retained only for haplotypes a later round may delta
// from -- retaining it for every haplotype ever scored would cost more than the placement it saves.
struct HaploidOriginCensus {
    std::vector<std::vector<HaploidOrigin>> origins_by_fragment;
    std::size_t haplotype_length = 0;
    std::uint64_t payload_bytes() const;
};

// What one delta actually did. Reuse is reported, never assumed: origins_reused against
// origins_recomputed is the claim this whole stage rests on, and a delta that quietly fell back to
// full placement would otherwise look like a fast one.
struct BlockDeltaWork {
    std::uint64_t origins_reused = 0;
    std::uint64_t origins_recomputed = 0;
    std::uint64_t local_window_bases = 0;
    std::uint64_t piece_hits = 0;
    std::uint64_t dp_executions = 0;
    std::uint64_t cache_misses = 0;
    // Requests and hits as well as misses. Reporting misses alone hides the probe VOLUME, which is
    // what decides whether an owned key's copy costs anything -- and a delta path that reported
    // hits as zero because nobody filled the field is worse than one that reported nothing.
    std::uint64_t cache_requests = 0;
    std::uint64_t cache_hits = 0;
    // Predicted from the census shape before it is built, observed from the containers afterwards.
    // Two numbers rather than one, because a predicted payload nobody checks is a guess.
    std::uint64_t census_payload_bytes_predicted = 0;
    std::uint64_t census_payload_bytes_observed = 0;
    std::uint64_t local_window_payload_bytes = 0;
    // Stage 3 accounts. A repair is held as SEGMENTS, so the bytes it materialises are the window it
    // actually looked at, not a 226kb haplotype and its reverse. full_haplotype_bytes_materialised is
    // therefore expected to be zero for a neighbour and non-zero only when one is promoted to
    // incumbent -- if it is not, the segmented representation is not doing its job.
    std::uint64_t full_haplotype_bytes_materialised = 0;
    std::uint64_t local_window_bytes_materialised = 0;
    std::uint64_t sequence_key_bytes = 0;
    // Where the time goes, split three ways rather than reported whole.
    // Container GROWTH events, which is where the allocations are. Persistent scratch means these
    // should fall to zero once warmed: a non-zero figure after the first round says something is
    // still being built per candidate rather than reused.
    // Stage 6: the window-origin cache and the preaggregated merge. Work that was LOADED is reported
    // apart from work that was PERFORMED, because a cache hit performs none of it.
    std::uint64_t origin_cache_hits = 0;
    std::uint64_t origin_cache_misses = 0;
    std::uint64_t origin_cache_stored = 0;
    // Demand, not just outcome: a miss whose records were built and then NOT retained because the
    // budget was full is a decline, and every later request for that window pays the enumeration again.
    std::uint64_t origin_cache_declined = 0;
    // Planned retention only: a miss on a window the round's plan did not select. Not a decline -- the
    // budget had room; the window's priority did not earn it.
    std::uint64_t origin_cache_not_selected = 0;
    // Planned retention within a round. A group reuse is a request served from the records the
    // immediately preceding enumeration built for the SAME window -- the canonical batch order makes
    // every window's requests adjacent -- so it executes no work and holds nothing extra. Evictions here
    // are the in-round replacements that made room for a higher-priority window.
    std::uint64_t origin_group_reuses = 0;
    std::uint64_t origin_cache_evicted = 0;
    std::uint64_t origin_cache_evicted_bytes = 0;
    std::uint64_t origin_cache_admitted_bytes = 0;
    std::uint64_t origin_cache_declined_bytes = 0;
    double seconds_reenumeration = 0.0;   // misses on a window already missed before
    std::uint64_t cached_origin_records_consumed = 0;
    std::uint64_t physical_origins_constructed = 0;
    std::uint64_t fragments_reduced = 0;
    std::uint64_t fragments_copied = 0;
    double seconds_enumeration = 0.0;
    double seconds_affected_filter = 0.0;
    double seconds_aggregate_merge = 0.0;
    double seconds_reduce = 0.0;
    std::uint64_t placement_row_constructions = 0;
    std::uint64_t accumulator_growths = 0;
    std::uint64_t origin_output_growths = 0;
    double seconds_placement_and_origins = 0.0;
    double seconds_statekey_merge = 0.0;
    double seconds_hashing_and_windows = 0.0;
    // True when the replacement sat within the end margin, or no census was held, so the haplotype
    // was placed in full. Still exact; simply not a saving.
    bool fell_back_to_full = false;
    double seconds = 0.0;
};

// A repair as segments, exposed so the extraction can be compared against the materialised
// equivalent directly. The context uses exactly these operations, so a gate that drives them is
// testing the real path rather than a parallel copy of it.
//
// All three describe the haplotype formed by replacing [block_begin, block_end) of `base` with
// `replacement`, and none of them materialises it.
std::size_t segmented_length(const std::string& base, std::size_t block_begin, std::size_t block_end,
                             const std::string& replacement);
std::uint64_t segmented_hash(const std::string& base, std::size_t block_begin, std::size_t block_end,
                             const std::string& replacement);
// `forward_out` receives [begin, end); `reversed_out` receives that same range reversed, which is
// what the reverse-role scan consumes.
void segmented_window(const std::string& base, std::size_t block_begin, std::size_t block_end,
                      const std::string& replacement, std::size_t begin, std::size_t end,
                      std::string& forward_out, std::string& reversed_out);

// What a context holds at one instant, measured from its containers rather than inferred from RSS.
// LOGICAL is what is in use; CAPACITY is what the containers keep reserved. Node-based containers
// (maps, hashed maps) do not expose their internals, so those are charged a declared per-node
// allowance and say so. The purpose is to compare the sum with RSS and see what RSS holds that the
// context cannot account for -- not to assign that remainder to anything without evidence.
// One row per DISTINCT window the origin path was asked for, in first-request order. Identity is the
// origin cache's own: window bytes, both end flags and the relative replacement interval. It records
// demand whether or not the window was admitted, so what a budget declines -- and what declining it
// costs on every later request -- is measured rather than estimated.
struct OriginDemandRow {
    std::uint64_t order = 0;
    std::uint64_t key_hash = 0;
    std::uint64_t window_bytes = 0;
    bool at_start = false;
    bool at_end = false;
    std::uint32_t rel_block_begin = 0;
    std::uint32_t rel_new_end = 0;
    std::uint64_t requests = 0;
    std::uint64_t hits = 0;
    std::uint64_t misses = 0;
    bool admitted = false;
    std::uint64_t declined = 0;               // misses built and not retained: the budget was full
    std::uint64_t entry_bytes = 0;            // declared cost of one entry; the same at every miss
    std::uint64_t records = 0;                // retained origin records in one entry
    std::uint64_t origins_constructed = 0;    // executed, summed over every miss
    double seconds_enumeration = 0.0;         // enumeration + record building, summed over misses
    double seconds_reenumeration = 0.0;       // the same, for misses after the first
    // MEASUREMENT ONLY, at the first request: the exact sorted StateKey -> multiplicity aggregate of
    // the affected retained origins -- the array the merge consumes -- and what storing it instead of
    // the raw records would cost. For an interior window every retained record is affected
    // (retained_not_affected must be 0), so the aggregate depends on the key alone. For a window at an
    // end the affected subset reads the haplotype length, so its figure describes one call only.
    std::uint64_t group_reuses = 0;               // served by the adjacent enumeration, no work
    std::uint64_t evictions = 0;                  // times this window's entry was evicted
    bool aggregate_measured = false;
    std::uint64_t aggregate_states = 0;
    std::uint64_t aggregate_fragments = 0;
    std::uint64_t aggregate_multiplicity = 0;
    std::uint64_t retained_not_affected = 0;
    std::uint64_t aggregate_entry_bytes = 0;         // StateCount records (24 bytes each)
    std::uint64_t aggregate_entry_bytes_packed = 0;  // the same states at 16 bytes each
};

// The identity of the window a repair WOULD request, computed without enumerating anything: the same
// bounds, bytes and key add_repaired_haplotype uses. already_held means the repaired haplotype is held,
// so the repair would request no window at all.
struct PlannedRepairWindow {
    std::uint64_t haplotype_hash = 0;
    bool already_held = false;
    std::uint64_t key_hash = 0;
    std::string bytes;
    bool at_start = false;
    bool at_end = false;
    std::uint32_t rel_block_begin = 0;
    std::uint32_t rel_new_end = 0;
};

// One repair of a batch: replace [block_begin, block_end) of `incumbent` with `replacement`.
struct RepairRequest {
    std::size_t incumbent = 0;
    std::size_t block_begin = 0;
    std::size_t block_end = 0;
    std::string replacement;
};

struct OriginCachePlanReport {
    std::uint64_t planned_repairs = 0;
    std::uint64_t window_requests = 0;     // distinct new haplotypes, each requesting one window
    std::uint64_t distinct_windows = 0;
    std::uint64_t known_windows = 0;       // exact entry size known from an earlier enumeration
    std::uint64_t selected_entries = 0;
    std::uint64_t selected_bytes = 0;
    std::uint64_t retained_entries = 0;    // selected and already held
    std::uint64_t evicted_entries = 0;
    std::uint64_t evicted_bytes = 0;
};

// ---- MEASUREMENT ONLY: the necessary-condition piece screen, in shadow ---------------------------
//
// Lemma. Partition a read into d disjoint pieces (offsets o_j, lengths l_j). Let an alignment of the
// read against a prefix of haplotype[s, s+W) use at most k unit-cost edits (substitutions, inserted
// and deleted bases) -- the existing feasibility prefilter's acceptance condition. A piece containing
// no edit is EXACT: it matches haplotype[s + o_j + delta_j, + l_j) with delta_j the alignment's
// displacement over it (haplotype bases consumed minus read bases consumed). Then:
//   (a) at least d - k pieces are exact -- each broken piece holds an edit, and pieces are disjoint;
//   (b) the exact pieces j_1 < ... < j_m lie in increasing, non-overlapping haplotype positions inside
//       the window: delta_b - delta_a >= -(read bases strictly between a and b);
//   (c) |delta_{j_1}| + sum |delta_{j_{t+1}} - delta_{j_t}| <= k, because the alignment starts at s
//       (delta 0 before the first read base) and every unit of displacement change is an indel base.
// (b) and (c) are the offset consistency. (a) and (c) are separate bounds on the same k edits: an
// indel inside a broken piece is counted by both, so they are NOT added.
//
// The shadow evaluates the screen for d = k+1, k+2 and 2k+1 on every start that reaches the
// prefilter, rejects nothing, and counts what it WOULD reject -- including any start the prefilter or
// the exact scorer accepted, which must be zero.
struct PieceScreenShadowRow {
    std::uint64_t pieces = 0;
    std::uint64_t entering = 0;                    // starts that reached the existing prefilter
    std::uint64_t count_pass = 0;                  // passed (a) alone
    std::uint64_t admitted = 0;                    // passed (a), (b) and (c)
    std::uint64_t rejected = 0;
    std::uint64_t rejected_prefilter_passed = 0;   // MUST be 0
    std::uint64_t rejected_accepted = 0;           // MUST be 0
    std::uint64_t screen_ns = 0;                   // query time
    std::uint64_t rejected_emission_ns = 0;        // gapped_emission time on starts it would reject
    std::uint64_t rejected_prefilter_ns = 0;       // of which inside the existing prefilter
};
struct PieceScreenShadowReport {
    std::uint64_t band = 0;
    std::uint64_t starts = 0;                      // verification calls observed
    std::uint64_t entering = 0;
    std::uint64_t prefilter_passed = 0;
    std::uint64_t accepted = 0;
    std::uint64_t emission_ns = 0;                 // all gapped_emission calls observed
    std::uint64_t prefilter_ns = 0;                // all prefilter executions observed
    std::uint64_t code_build_ns = 0;               // substring codes of the scanned haplotypes/windows
    std::uint64_t layout_ns = 0;                   // read piece layouts
    PieceScreenShadowRow rows[3];
};
void set_piece_screen_shadow(bool enabled);   // enabling also resets the totals
PieceScreenShadowReport piece_screen_shadow_report();
// The screen on one placement, the independent brute-force reference used to gate it, the count
// condition alone, and the existing prefilter under the same window rule. For gates only.
bool piece_screen_admits(const std::string& haplotype, std::size_t start, const std::string& read,
                         std::size_t band, std::size_t max_gap_bp, std::size_t pieces);
bool piece_screen_admits_reference(const std::string& haplotype, std::size_t start,
                                   const std::string& read, std::size_t band, std::size_t max_gap_bp,
                                   std::size_t pieces);
bool piece_screen_count_passes(const std::string& haplotype, std::size_t start, const std::string& read,
                               std::size_t band, std::size_t max_gap_bp, std::size_t pieces);
bool feasibility_prefilter_admits(const std::string& haplotype, std::size_t start,
                                  const std::string& read, std::size_t band, std::size_t max_gap_bp);
// Gates only: the prefilter's decision by the reference (every column) and the production (early-refusing)
// kernels, each with the read's masks prepared ahead and built on the spot. Exactness means all four agree.
struct PrefilterDecisions {
    bool reference_unprepared = false, reference_prepared = false;
    bool production_unprepared = false, production_prepared = false;
};
PrefilterDecisions feasibility_prefilter_decisions(const std::string& haplotype, std::size_t start,
                                                   const std::string& read, std::size_t band, std::size_t max_gap_bp);
// Gates only: proposed (form, start) pairs ordered and deduplicated by std::sort + std::unique (bucketed false) or
// by the production counting-by-form path (true). Exactness means the two lists are equal.
std::vector<std::pair<std::uint32_t, std::uint32_t>> sort_unique_proposals_scatter_for_test(
    std::vector<std::pair<std::uint32_t, std::uint32_t>> proposals, std::size_t form_count);
std::vector<std::pair<std::uint32_t, std::uint32_t>> sort_unique_proposals_for_test(
    std::vector<std::pair<std::uint32_t, std::uint32_t>> proposals, std::size_t form_count, bool bucketed);

// ---- MEASUREMENT ONLY: one window's enumeration split across fragment ranges ---------------------
struct WindowEnumerationBench {
    std::size_t workers = 0;
    std::size_t fragments = 0;
    std::size_t window_bytes = 0;
    double reference_cold_seconds = 0.0;       // production path, one thread, emission cache cold
    double production_warm_seconds = 0.0;      // production path, one thread, emission cache warm
    std::vector<double> wall_seconds;          // partitioned runs, one per repeat
    std::vector<double> worker_busy_seconds;   // last repeat
    double prefilter_seconds = 0.0;            // summed over workers, last repeat
    double emission_seconds = 0.0;             // gapped_emission, summed over workers, last repeat
    bool origins_identical = false;            // against the cold production reference
    bool counters_identical = false;           // against the warm production pass
    std::string first_difference;
    std::uint64_t origins = 0;
    std::uint64_t piece_hits = 0;
    std::uint64_t latent_pairs = 0;
    std::uint64_t executed_dps = 0;
    std::uint64_t executed_prefilters = 0;
    std::uint64_t verifications = 0;
    std::uint64_t duplicate_scan_bases = 0;    // extra window scanning the split performs
};

struct DemandWindow {
    bool found = false;
    std::string bytes;
    bool at_start = false;
    bool at_end = false;
    std::uint32_t rel_block_begin = 0;
    std::uint32_t rel_new_end = 0;
};

// Every thread doing scoring work -- placement, combination, window enumeration -- holds one slot of
// this process-wide gauge while it works. The peak is how a gate proves the total never exceeds the
// requested thread count, whichever stages run.
std::size_t peak_active_workers();
void reset_peak_active_workers();

// Window enumeration split across the context's workers. The work counters and every result are the
// one-thread values exactly; these are what parallelism added on top.
struct ParallelEnumerationStats {
    std::size_t workers = 1;                    // per enumeration: min(threads, fragments)
    std::uint64_t parallel_enumerations = 0;
    std::uint64_t duplicate_dps = 0;            // recursions two workers both ran for one new key
    std::uint64_t duplicate_scan_bases = 0;     // window scanning repeated by workers 1..W-1
};

// I2-B (EXPERIMENTAL): round-local fine discovery inside the adaptive search's window enumeration
// (PreparedScoringContext::enable_fine_window_discovery). Tallies are cumulative over the search. seconds_rows /
// seconds_pairing time every production window enumeration, either path, summed over enumeration workers (worker
// elapsed); `enumerations` counts windows enumerated.
struct FineWindowReport {
    bool enabled = false;
    std::uint64_t budget_bytes = 0, served_cap = 0;
    std::uint64_t rounds = 0;
    std::uint64_t worker_windows = 0, fallback_windows = 0;   // per worker store: windows built, of which by the reference
    std::uint64_t units_discovered = 0, units_served = 0, units_not_retained = 0;
    std::uint64_t placements_discovered = 0, placements_served = 0, hash_collisions = 0;
    // A container ledger, not process memory. RETAINED (held, what the budget bounds): the held units' byte capacity, their
    // rows' allocated capacity, an ESTIMATED 48 bytes per hash node and the bucket array. TEMPORARY: the reusable discovery
    // scratch (rows, proposals, touched lists), units used once and not retained, and each window's candidate unit copies.
    // Each is the max over rounds of the workers' summed per-round peaks (the workers' peaks need not coincide in time).
    std::uint64_t round_peak_bytes = 0;            // retained + temporary
    std::uint64_t round_peak_held_bytes = 0;       // retained only
    std::uint64_t round_peak_temporary_bytes = 0;  // temporary only
    std::uint64_t enumerations = 0;
    double seconds_rows = 0.0, seconds_pairing = 0.0;
};

struct ContextMemoryAccount {
    std::uint64_t slots_total = 0;
    std::uint64_t slots_live = 0;
    std::uint64_t slots_released = 0;
    std::uint64_t slot_table_bytes = 0;          // per-slot bookkeeping, released slots included
    std::uint64_t sequence_bytes = 0;            // forward + reverse haplotype bytes in use
    std::uint64_t sequence_capacity = 0;
    std::uint64_t segment_bytes = 0;             // segmented repairs: node + replacement
    std::uint64_t mass_column_bytes = 0;
    std::uint64_t mass_column_capacity = 0;
    std::uint64_t census_count = 0;
    std::uint64_t census_bytes = 0;
    std::uint64_t census_capacity = 0;
    std::uint64_t origin_cache_entries = 0;
    std::uint64_t origin_cache_payload = 0;      // the declared, budgeted figure
    std::uint64_t origin_cache_capacity = 0;
    std::uint64_t merge_working_bytes = 0;       // group preaggregates + full aggregates
    std::uint64_t merge_working_capacity = 0;
    std::uint64_t origin_demand_bytes = 0;       // the demand ledger: one row + key bytes per window
    std::uint64_t window_scratch_capacity = 0;   // persistent placement rows and flat origin output
    std::uint64_t origin_scratch_capacity = 0;   // miss/check records, candidate states, merge buffer
    std::uint64_t emission_cache_table = 0;
    std::uint64_t emission_cache_arena = 0;      // components + content slabs
    std::uint64_t prepared_capacity = 0;         // read forms, piece indexes, fragments: constant
    std::uint64_t last_scoring_result_bytes = 0; // flattened matrix of the latest score_candidates
    std::uint64_t logical_total = 0;
    std::uint64_t capacity_total = 0;
};

class PreparedScoringContext {
public:
    // Fragments, library and options are fixed for the context's life: they define the model, and a
    // context whose model changed between batches would retain columns computed under the old one.
    PreparedScoringContext(std::vector<ReadPair> fragments,
                           LibraryModel library,
                           FragmentScoringOptions options);
    ~PreparedScoringContext();
    PreparedScoringContext(const PreparedScoringContext&) = delete;
    PreparedScoringContext& operator=(const PreparedScoringContext&) = delete;

    // Places and scores only the sequences not already held; returns one haplotype index per offered
    // sequence, in the order offered. Duplicate bytes within the batch resolve to one column, exactly
    // as intern_haplotypes collapses them in the one-shot path.
    // `reuse_released_slots` places new sequences into released repair slots before appending any, as a repair
    // does, so a caller that places and releases repeatedly (a tie-branch start materialised just in time) does
    // not grow the slot table. Off by default: every other caller's slot numbering is unchanged.
    // `retain_censuses` (I2-D): the batch's resulting slots also keep their origin censuses, captured during the SAME
    // placement pass for every slot this batch places (one placement per distinct new sequence; byte-identical
    // sequences share one slot and one census). A slot that was already held without a census is censused as
    // retain_census does; one already holding a census is left as it is. Off: unchanged.
    std::vector<std::size_t> add_haplotypes(const std::vector<std::string>& sequences,
                                            IncrementalBatchReport* report = nullptr,
                                            bool reuse_released_slots = false,
                                            bool retain_censuses = false);
    // I2-D: full haplotype placements this context has performed (one per slot add_haplotypes places, one per
    // retain_census that places). Measurement and gates only.
    std::uint64_t full_placement_count() const;

    // Combination from the retained columns. No placement happens here, so this is cheap and may be
    // rerun over any candidate set the context already holds haplotypes for.
    FragmentScoringResult score_candidates(const std::vector<HaplotypeCandidate>& candidates) const;

    std::size_t unique_haplotypes() const;
    // Haplotypes that have ever occupied a slot. Released repair slots are reused, so this is not the
    // same as unique_haplotypes(): that is how many slots exist, this is how many haplotypes have held
    // one. Reuse is what keeps the slot count -- and everything sized by it -- bounded across rounds.
    std::size_t haplotypes_created() const;
    double mass_of(std::size_t fragment, std::size_t haplotype) const;
    double exposure_of(std::size_t haplotype) const;
    const FragmentScoringCounters& cumulative_counters() const;
    const EmissionCacheCounters& cumulative_emission_cache() const;
    // True once a cap has been crossed. The context is then permanently unusable.
    bool refused() const;

    // ---- stage 2 ------------------------------------------------------------------------------
    //
    // There is no freeze any more, and removing it is the point of stage 3. A cache key used to hold
    // raw pointers into the read and the haplotype window it was built from, so anything transient
    // had to be kept out of the table -- which cost a recomputed recursion on every later miss, about
    // 171,000 of them per round. Keys now OWN the bytes they compare, so a transient repair haplotype
    // may be destroyed while every entry it produced stays valid and keeps serving.
    //
    // Marks the point after which added haplotypes are transient and may be released. It may be called
    // before anything is added, which makes every haplotype -- the seed included -- releasable once it
    // is no longer kept. Repair search does that: a seed that has stopped being an incumbent is never
    // delta'd from again, so its full census is pure retention.
    void begin_transient_round();

    // Scores the incumbent haplotype with [block_begin, block_end) replaced by `replacement`, reusing
    // every origin the replacement cannot touch. Returns the new haplotype's index. Requires a
    // retained census for `incumbent`; refuses otherwise rather than silently placing in full,
    // because a silent fallback is indistinguishable from a delta that does not work.
    //
    // `retain_census` decides whether the new haplotype keeps a census of its own. It defaults to
    // false because retaining one per neighbour is not affordable -- a census runs to hundreds of
    // thousands of origins, and a round offers hundreds of neighbours -- so the caller scores a
    // round's neighbours without censuses and asks for one only for the winner it moves to.
    std::size_t add_repaired_haplotype(std::size_t incumbent,
                                       std::size_t block_begin,
                                       std::size_t block_end,
                                       const std::string& replacement,
                                       BlockDeltaWork* work = nullptr,
                                       bool retain_census = false);

    // ---- stage 6 (diagnostic path only; stage 5 remains the default and the reference) -------
    //
    // Routes add_repaired_haplotype through the window-origin cache and the preaggregated sorted
    // merge. Off by default, in which case the stage-5 code runs unchanged.
    void enable_origin_cache_path(bool enabled);
    // I2-B (EXPERIMENTAL, off by default): within each add_repaired_haplotypes round, every window enumeration builds its
    // raw placement rows from a round-local, per-enumeration-worker fine store (FineWindowReport) instead of scanning the
    // whole window; pairing, filtering, merging and reduction are unchanged, and so is every row. `budget_bytes` bounds
    // the stores' allocated capacity (split evenly across workers); `max_served_placements` caps placements served from
    // them over the search (crossing it refuses). Single add_repaired_haplotype calls outside a round, full placements
    // and every self-check keep the reference path.
    void enable_fine_window_discovery(std::size_t budget_bytes, std::uint64_t max_served_placements);
    FineWindowReport fine_window_report() const;
    // Bytes the origin cache may retain (0: retain nothing, recompute every window exactly), and a cap
    // on cached origin records CONSUMED over the search (0: none). Crossing the budget stops
    // retention only; crossing the consumption cap refuses, like every other declared cap.
    void set_origin_cache_limits(std::size_t budget_bytes, std::uint64_t max_records_consumed);
    std::size_t origin_cache_entries() const;
    std::uint64_t origin_cache_bytes() const;
    std::uint64_t origin_cache_not_retained() const;
    std::uint64_t origin_records_consumed() const;
    // Working storage for the preaggregated merge: per-group unaffected StateKey arrays and lazily
    // built full aggregates. Reported apart from the cache budget, because it is per-round working
    // state released with the incumbent, not retained across rounds.
    std::uint64_t merge_working_bytes() const;
    ContextMemoryAccount memory_account() const;
    std::vector<OriginDemandRow> origin_demand() const;
    // Planned selective retention. plan_repaired_window describes the window a repair would request;
    // plan_origin_cache then chooses the round's cache set from those windows and evicts every entry
    // outside it, releasing its storage. Selection is deterministic and independent of request order:
    // a window's priority is the number of distinct new haplotypes requesting it divided by its EXACT
    // entry bytes (compared as integers), ties broken by key hash and then the full key. Only windows
    // whose exact size is already known -- enumerated in an earlier round -- are candidates, so no
    // entry is ever retained speculatively. Once planned, a miss is retained only if its window was
    // selected; the selected total never exceeds the budget, so neither does the cache.
    PlannedRepairWindow plan_repaired_window(std::size_t incumbent, std::size_t block_begin,
                                             std::size_t block_end, const std::string& replacement) const;
    OriginCachePlanReport plan_origin_cache(const std::vector<PlannedRepairWindow>& windows);
    // A whole round's repairs at once, under planned retention. The batch is planned, then executed in
    // a CANONICAL order -- by window identity, then haplotype, then the request itself -- so the order
    // the caller enumerated its candidates in cannot change which windows are enumerated, reused,
    // retained or evicted, nor any score. Requests for one window are adjacent: the first enumerates,
    // the rest are group reuses. A window whose exact size was unknown before this round, or that the
    // plan selected, is retained on its miss by priority (distinct requesting haplotypes per exact
    // entry byte, the plan's total order), evicting strictly lower-ranked entries first so the payload
    // never exceeds the budget. Returns one slot per request, in the caller's order.
    std::vector<std::size_t> add_repaired_haplotypes(const std::vector<RepairRequest>& requests,
                                                     std::vector<BlockDeltaWork>* works = nullptr,
                                                     OriginCachePlanReport* plan = nullptr);
    std::vector<std::uint64_t> origin_cache_key_hashes() const;   // live entries, ascending
    // MEASUREMENT ONLY. The recorded window with this key hash, if the demand ledger holds it; and one
    // window's enumeration timed with its fragments divided among `workers` contiguous ranges.
    DemandWindow demand_window(std::uint64_t key_hash) const;
    ParallelEnumerationStats parallel_enumeration_stats() const;
    // A digest of the latest window enumeration's origins and offsets, bit for bit. For gates only.
    std::uint64_t last_window_digest() const;
    WindowEnumerationBench benchmark_window_enumeration(const std::string& bases, bool at_start, bool at_end,
                                                        std::size_t workers, std::size_t repeats);
    std::uint64_t origin_cache_slots() const;                     // entry slots, live or reusable
    // Slots currently holding a full origin census, ascending. A caller that knows what each slot is
    // (an incumbent, a seed, a winner being promoted) can then say why every one of them is live.
    std::vector<std::size_t> census_slots() const;

    // Retains an origin census for this haplotype, so a later round may delta from it.
    void retain_census(std::size_t haplotype);
    bool holds_census(std::size_t haplotype) const;
    // The retained census of `haplotype`, or null when none is held. Read-only; never consulted by a score, a delta or
    // a cache (regression tests compare censuses through it).
    const HaploidOriginCensus* census_view(std::size_t haplotype) const;

    // Releases every haplotype added since the last checkpoint that is not named in `keep`: its
    // sequence, its reversed copy, its mass column and its census. Safe because every cache entry
    // owns the bytes it compares (the emission cache in its arena, the origin cache in its entries),
    // so nothing retained points into what is being freed.
    void release_transient(const std::vector<std::size_t>& keep);
    // Bytes currently held by transient sequences, mass columns and censuses, so a caller can check
    // that storage returns to its previous level after every round instead of trusting that it does.
    std::uint64_t transient_payload_bytes() const;
    // What one slot holds right now -- sequence, reversed copy, segment, mass column and census -- counted exactly as
    // transient_payload_bytes counts it; 0 for a released slot. What a caller keeping the slot alive retains.
    std::uint64_t slot_payload_bytes(std::size_t slot) const;
    std::uint64_t retained_census_bytes() const;
    // I4-C: the transient slots that currently hold anything (sequence, mass column, segment or census), ascending.
    std::vector<std::size_t> live_transient_slots() const;
    // I4-C: releases exactly the named transient slots (each must be a live slot added since the last checkpoint), with
    // the same per-slot release as release_transient -- and the group and full aggregates keyed by them as incumbent.
    // Nothing else is touched. Used to give back a rejected tie-branch seed's own slots at once.
    void release_slots(const std::vector<std::size_t>& slots);

private:
    void release_slot_storage(std::size_t slot);   // I4-C: the one per-slot release both public releases use
    struct State;
    std::unique_ptr<State> state_;
};

// ---- EXPERIMENTAL (I1): shared-content placement index -------------------------------------------------------------
//
// Outside every production route. A start s of a read form on a strand text T is proposed and verified from the bytes of
// D(s) = [s - g, s + L + g) alone (piece hits implying s lie in it; the verification window is [s, s + min(available,
// L + g)), with the normaliser clamped at L + g), plus clipping where D(s) meets an end of T. So each strand of each
// haplotype decomposes into content UNITS, each evaluated once per read batch with the scorer's own proposal and
// verification (collect_form_placements) and then expanded to every occurrence:
//   * a segment interior (strand, bytes): it keeps the starts whose whole D(s) lies inside the segment;
//   * a boundary window (strand, bytes of T[b - W, b + W) clipped to T, whether that clip reached either end of T, and b's
//     offset in the window), W = max read length + 2g, for every segment boundary b and both ends of T: it keeps the
//     starts whose D(s) crosses b -- all of D(s) lies in the window, so the result is the text's own.
// Every start falls in one of the two (D(s) crosses some boundary or end, or lies inside one segment); a start owned by
// several windows is the same placement and is kept once. Copies are separate occurrences at separate starts, so copy
// multiplicity is kept. Unit keys are compared byte for byte (the hash only selects). Pairing, physical origins,
// reduction, exposure and the diploid combination are the scorer's own.
struct ContentIndexReport {
    std::uint64_t haplotypes = 0;
    std::uint64_t segment_occurrences = 0, window_occurrences = 0;   // registered, both strands
    std::uint64_t interior_units = 0, window_units = 0;              // distinct
    std::uint64_t unit_bytes = 0;                                    // distinct unit bytes held
    std::uint64_t occurrence_bytes = 0;
    std::uint64_t hash_collisions = 0;                               // equal hash, different key: separated by bytes
    double seconds_registration = 0.0;                               // cumulative add_haplotype wall time
    // I1.4 fine discovery: distinct fine units, their text bytes and owned core bytes; fine pieces registered over all
    // distinct coarse interiors; distinct coarse interior and window bytes.
    std::uint64_t fine_units = 0, fine_unit_bytes = 0, fine_core_bytes = 0, fine_occurrences = 0;
    std::uint64_t coarse_interior_bytes = 0, window_bytes = 0;
};

class ContentPlacementIndex {
public:
    // Bound to a maximum read length and max_gap_bp: windows are sized from both, and a batch with longer reads or a
    // different gap is refused.
    // fine_discovery (I1.4): coarse interior segments are discovered through content-defined fine units (owned cores
    // with their dependency context) and their placement rows reconstructed; everything downstream is unchanged.
    ContentPlacementIndex(std::size_t max_read_bp, std::size_t max_gap_bp, bool fine_discovery = false);
    ~ContentPlacementIndex();
    ContentPlacementIndex(const ContentPlacementIndex&) = delete;
    ContentPlacementIndex& operator=(const ContentPlacementIndex&) = delete;
    // Registers one distinct haplotype sequence with a segmentation (lengths summing to its length; zeros allowed) and
    // returns its slot. Any segmentation is exact; a finer one shares more.
    std::size_t add_haplotype(const std::string& bases, const std::vector<std::uint32_t>& segment_lengths);
    std::size_t haplotype_count() const;
    std::size_t unit_count() const;
    ContentIndexReport report() const;
    struct State;
    std::unique_ptr<State> state_;
};

void set_content_index_hash_collapse(bool collapse);   // TEST ONLY: one hash for every content unit
// I2-C TEST ONLY: request spare capacity per stored vector (unused capacity the ledger must charge), and the negative
// control that charges occupied rows only.
void set_store_spare_capacity(std::size_t slots);
void set_store_charge_rows_only(bool on);
// I2-B TEST ONLY: compare every fine window's raw rows with the reference's two scans (throws std::logic_error on a
// difference); and deliberately break ownership or context (0: none; 1 overlap; 2 left context; 3 right context; 4 end
// window) so the comparison can be shown to catch it.
void set_fine_window_row_check(bool on);
void set_fine_window_mutation(int mutation);

struct ContentBatchReport {
    std::uint64_t units_evaluated = 0, unit_bases_scanned = 0, discovery_calls = 0;
    std::uint64_t proposed_starts = 0, unique_starts = 0, full_read_verifications = 0, dp_evaluations = 0;
    std::uint64_t unit_placements_accepted = 0, unit_placements_kept = 0, placement_store_bytes = 0;
    // I2-K (measurement): the emission cache as discovery left it, just before it was discarded -- entries held and the
    // storage actually allocated (table capacity, key-byte slabs, component-arena capacity) and the budgeted bytes.
    std::uint64_t emission_cache_entries_held = 0, emission_cache_table_bytes = 0, emission_cache_slab_bytes = 0,
                  emission_cache_component_bytes = 0, emission_cache_budgeted_bytes = 0;
    double seconds_discovery = 0.0;                       // wall, cumulative over discover() calls
    double worker_seconds_scan = 0.0, worker_seconds_sort = 0.0, worker_seconds_verification = 0.0;
    double seconds_expansion_and_pairing = 0.0;           // wall, last score()
    double worker_seconds_expansion = 0.0, worker_seconds_pairing = 0.0;
    std::uint64_t expanded_placements = 0, duplicate_placements_removed = 0, physical_origins = 0, state_classes = 0;
    double seconds_candidate_scoring = 0.0;
    // I1.1: the read batch's exact prefilter-decision cache (FragmentScoringOptions::union_prefilter_cache_bytes):
    // "on", "off", or an explicit fallback to uncached discovery; its report is cumulative over discover() calls.
    std::string prefilter_status = "off";
    UnionPrefilterCacheReport prefilter;
    std::uint64_t executed_prefilters = 0, executed_dps = 0, prefilter_rejections = 0;   // executed work    // I1.2: interior-origin reuse (FragmentScoringOptions::union_origin_cache_bytes as the batch's declared budget):
    // aggregates built from this batch's verified unit placements by the reference's own phase-2 builder, merged at
    // every occurrence through guarded zones. physical_origins is then direct + served (logical).
    std::string origin_reuse_status = "off";
    std::uint64_t origin_budget_bytes = 0, origin_held_bytes = 0;
    std::uint64_t origin_entities_eligible = 0, origin_entities_built = 0, origin_entities_admitted = 0, origin_entities_declined = 0;
    std::uint64_t origins_constructed = 0, direct_physical_origins = 0, served_physical_origins = 0, zones_applied = 0;
    double seconds_origin_construction = 0.0;
    // I1.3: retention (append-only, within this context). *_last: the last score() call; *_total: cumulative.
    // physical_origins, state_classes, expanded_placements and the origin split are totals over every completed slot.
    std::uint64_t columns_computed_last = 0, columns_computed_total = 0, columns_retained_last = 0;
    std::uint64_t pairs_requested_last = 0, pairs_computed_last = 0, pairs_computed_total = 0, retained_pairs = 0;
    std::uint64_t retained_column_bytes = 0, retained_pair_bytes = 0;
    double seconds_pair_lookup_and_combination = 0.0, seconds_materialisation = 0.0;
    // I1.4: discovery work under fine units
    std::uint64_t units_discovered_directly = 0, owned_starts_skipped_before_verification = 0;
    std::uint64_t coarse_units_reconstructed = 0, reconstructed_placements = 0;
    double seconds_reconstruction = 0.0;
    // I2-A: guards and executed work, in the reference's terms. Discovery: piece hits and haplotype bases walked
    // (cumulative), the placement store charged against its declared budget (stored rows, coarse plus fine, x row size;
    // 0: no budget) and the workers that ran. score(): `counters_last` is the last call's EXECUTED work -- discovery's
    // (cumulative) seeding counters, and the origins, state classes, latent pairs and valid states actually enumerated
    // by direct pairing and by aggregate construction -- with what was SERVED in its reused_* fields and the merged
    // (logical) census in merged_state_classes, as the reference reports its phase-2 pass. The acceptance caps are held
    // to these executed totals; the served and expanded quantities are held to the consumption caps.
    std::uint64_t piece_hits = 0, candidate_bases_scanned = 0;
    std::uint64_t store_budget_bytes = 0, store_charged_bytes = 0;
    std::size_t threads_used_discovery = 0, threads_used_pairing = 0;
    FragmentScoringCounters counters_last;
    std::uint64_t placements_consumed_last = 0;   // expanded placements served from the store (cap: cached placements)
    std::uint64_t origins_consumed_last = 0, states_consumed_last = 0;   // served at zones (caps: cached origins/states)
    std::uint64_t physical_origins_observed = 0;  // the shared in-flight origin count (what a refusal saw)
    EmissionCacheCounters emission_cache;         // cumulative over discover() calls
    // I2-C: the placement-store ledger (see ContentIndexBatch::State): per-unit rows' allocated capacity plus the per-unit
    // container's capacity -- a container ledger, NOT process memory. live: now; peak: the most held at once (fine and
    // coarse rows overlap during reconstruction); allocated / released: cumulative. worker_scratch_bytes: the workers'
    // temporary discovery buffers (rows, proposals, the per-unit buffer) held at once, plus reconstruction's use counts --
    // temporary, outside the budget. store_charged_bytes is what the budget was charged: allocated capacity (exact
    // layout) or occupied rows (legacy).
    std::string store_layout = "exact";
    std::uint64_t store_live_bytes = 0, store_peak_bytes = 0, store_allocated_bytes = 0, store_released_bytes = 0;
    std::uint64_t reconstruction_peak_bytes = 0, worker_scratch_bytes = 0;
    // I4-B (EXPERIMENTAL): "reconstructed" (coarse interior rows copied from their fine pieces) or "fine_views" (coarse
    // interior units served as views over the fine units' shared rows; nothing copied). Under views: the coarse units
    // served as views, the rows they represent (what reconstruction would have copied), the view table's charged bytes
    // (store metadata), and the origin-zone gathers -- the representation's transition-time copies, charged to the store
    // ledger while they exist -- their peak simultaneous bytes, cumulative bytes and gather time.
    std::string evidence_layout = "reconstructed";
    std::uint64_t view_coarse_units = 0, view_rows_represented = 0, view_table_bytes = 0;
    std::uint64_t view_zone_temporary_peak_bytes = 0, view_zone_temporary_bytes_total = 0;
    double seconds_view_zone_gather = 0.0;
};

class ContentIndexBatch {
public:
    // One read batch over a prepared index (the index is shared, not copied; it must outlive the batch).
    ContentIndexBatch(const ContentPlacementIndex& index, std::vector<ReadPair> fragments, LibraryModel library,
                      FragmentScoringOptions options);
    ~ContentIndexBatch();
    ContentIndexBatch(const ContentIndexBatch&) = delete;
    ContentIndexBatch& operator=(const ContentIndexBatch&) = delete;
    // Proposes and verifies every unit of the index this batch has not evaluated yet (all, first; only new units after
    // the index grew).
    void discover();
    // Expands, pairs and scores; candidate members are index slots. Requires every unit discovered.
    FragmentScoringResult score(const std::vector<HaplotypeCandidate>& candidates);
    // I1.3: completed slot columns and pair likelihoods are retained across score() calls (the default). false
    // reproduces I1.2: every call re-expands and re-pairs every slot and recombines every requested pair.
    void set_retention(bool retain);
    // I2-A: a declared byte budget for the placement store (stored rows x row size, fine and coarse); crossing it
    // refuses with ResourceExhausted and no ranking. 0 (the default): no budget, as before.
    void set_placement_store_budget(std::size_t bytes);
    // I2-A: one discovery pass only. After it the fine units' rows -- needed only to reconstruct coarse units added
    // later -- are released, and a later discover() that finds new units is refused (std::logic_error).
    void set_one_shot(bool one_shot);
    // I2-C: the I2-A layout (push_back growth, the budget charging occupied rows), for comparison only. Default: exact.
    void set_legacy_store_layout(bool legacy);
    // I4-B (EXPERIMENTAL): serve coarse interior units as views over their fine pieces' rows instead of reconstructing
    // copied coarse rows (fine discovery and the exact layout only). The fine rows are then kept until the batch ends.
    // Masses, exposures, scores and every cap's executed and served work are unchanged; the store holds the distinct
    // evidence (fine and window rows), the view table and in-flight origin-zone gathers.
    void set_fine_views(bool on);
    // TEST ONLY: an independent recount of the placement containers' allocated capacity (what the ledger must equal).
    std::uint64_t store_capacity_recount() const;
    const ContentBatchReport& report() const;
    // TEST ONLY: every discovered coarse interior unit's placement rows, keyed by strand and bytes.
    std::map<std::string, std::vector<std::tuple<std::uint32_t, std::uint32_t, std::uint32_t, std::uint32_t, double>>>
    interior_placement_rows() const;
    struct State;
    std::unique_ptr<State> state_;
};

// ---- I2-A (EXPERIMENTAL): the I1.4 indexed route at the caller's one-shot union boundary ----------------------------
// A thin adapter with score_diplotypes' contract: the same inputs, the same result structure (scores, audit scores,
// masses, exposures, counters, phase timings, workers, cache reports) and the same refusal semantics -- a declared cap
// throws ResourceExhausted with no ranking, after filling `observed_on_failure`. It builds a fine-discovery
// ContentPlacementIndex over `segments` (per unique haplotype, summing to its length) and one ContentIndexBatch,
// discovers once and scores once, then releases both before returning or rethrowing: nothing indexed outlives the call.
//
// Guards, each the reference's own where one applies:
//   * input validation and positive caps (score_diplotypes' rules); per-origin and per-haplotype diagnostics refused;
//   * max_mass_cells, before any column is allocated;
//   * max_piece_hits in flight (shared, batched publication) during discovery; max_proposed_starts_per_haplotype per
//     discovered TEXT (a unit is a substring of a haplotype, so it bounds the same per-call buffer, but a haplotype the
//     reference would refuse may pass: no refusal equivalence is claimed);
//   * max_physical_origins in flight over direct pairing AND aggregate construction (both executed work);
//   * the acceptance caps (piece hits, physical origins, state classes, latent pairs) on the merged executed totals;
//   * union_origin_cache_bytes bounds the held aggregates; served origins and state aggregates are charged to
//     max_union_cached_origins_consumed / max_union_cached_states_consumed (both required positive when it is on);
//   * every expanded placement is served from the store, so all of them are charged to
//     max_union_cached_placements_consumed (required positive);
//   * union_prefilter_cache_bytes / max_union_cached_prefilters_consumed: the shared decision cache and its in-flight
//     cap; mixed read lengths are refused exactly as the reference refuses them (no silent uncached fallback);
//   * `store_budget_bytes` (required positive) bounds the placement store, fine and reconstructed coarse rows.
// Counters are EXECUTED work (ContentBatchReport::counters_last); they need not equal the reference's, and a cap can
// refuse on one route and not on the other. Score identity is not resource-policy equivalence.
struct IndexedUnionReport {
    std::size_t max_read_bp = 0, store_budget_bytes = 0;
    double seconds_index_build = 0.0, seconds_preparation = 0.0, seconds_release = 0.0;
    std::string refusal;   // the ResourceExhausted message, if the call was refused
    ContentIndexReport index;
    ContentBatchReport batch;
};
FragmentScoringResult score_diplotypes_indexed_union(const InternedHaplotypes& haplotypes,
                                                     const std::vector<HaplotypeCandidate>& candidates,
                                                     const std::vector<ReadPair>& fragments, const LibraryModel& library,
                                                     const FragmentScoringOptions& options,
                                                     const std::vector<std::vector<std::uint32_t>>& segments,
                                                     std::size_t store_budget_bytes, IndexedUnionReport* report = nullptr,
                                                     FragmentScoringObservation* observed_on_failure = nullptr,
                                                     bool legacy_store_layout = false, bool fine_views = false);

} // namespace panvar
