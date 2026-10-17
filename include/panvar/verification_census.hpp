#pragma once

// I2-I (MEASUREMENT ONLY): a bounded census of verification requests -- how much verification work repeats for
// equivalent (scoring context, read content, target-window content), and where the repetition comes from.
//
// Exact counters cover every request. Repetition is measured on a prospectively fixed, content-based sample: a
// request is sampled exactly when its READ CONTENT's existing hash (the emission key prefix, computed once per read
// form at preparation) falls in 1 of `period` residues, so every occurrence of a sampled read is in the sample and
// its repeats are counted exactly. Hashes only select the sample and address the table; identity is established by
// the scoring context's exact parameters and byte comparison of the read and of the window. The window is the
// haplotype bases from the start, at most read length + max_gap_bp -- exactly what the feasibility decision and the
// emission read (the clamped available length is the window's length).
//
// The census never changes a verdict, a row or a counter of inference. It does not merge physical observations: it
// counts them. Storage (tables, window and read arenas, metadata) is charged against a hard byte budget BEFORE it is
// allocated; when the next allocation would exceed it, coverage is marked incomplete and new keys are dropped (and
// counted), while inference continues unchanged.

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace panvar {

// Which part of the run is verifying: set by the caller (genotype_command) at stage boundaries.
enum class CensusStage : std::uint8_t { other = 0, union_discovery = 1, adaptive = 2, dosage = 3 };
// Which placement routine is verifying: set by that routine for its own thread.
enum class CensusSite : std::uint8_t { other = 0, indexed_unit = 1, full_haplotype = 2, window = 3 };

// Report groups: the union pass; the adaptive search's full-haplotype placements (anchor seeding and any other full
// placement) and its window enumeration; dosage's full-haplotype placements and windows; anything else.
enum class CensusGroup : std::uint8_t {
    union_discovery = 0, adaptive_full = 1, adaptive_window = 2, dosage_full = 3, dosage_window = 4, other = 5
};
constexpr std::size_t kCensusGroups = 6;
const char* census_group_name(std::size_t group);

// What a request's result depends on besides its read and window content.
struct CensusContext {
    double substitution_error = 0.0, gap_open = 0.0, gap_extend = 0.0;
    std::size_t max_gap_bp = 0, band = 0;
};

// One verification request, as the placement loop saw it.
struct CensusOutcome {
    bool consulted = false;        // the union decision cache was consulted (union only)
    bool served = false;           // ... and held the decision
    int retained = -1;             // after a miss: 1 stored, 0 not retained (admission stopped), 2 already present
    bool kernel = false;           // the feasibility kernel executed
    bool feasible = false;         // the decision (remembered or computed) admitted the request
    bool emission_hit = false;     // the emission cache (or a deferred earlier result) answered it
    bool recursion = false;        // the recursion executed
    bool fixed_recursion = false;  // ... on the fixed-geometry backend (I2-H)
    bool accepted = false;         // a fit with finite emission
    std::uint32_t rows = 0;        // placement rows emitted for it
};

struct CensusExact {   // every request
    std::uint64_t requests = 0, units = 0, consulted = 0, served = 0, missed = 0, missed_not_retained = 0,
                  missed_already_present = 0, kernels = 0, feasibility_rejections = 0, emission_hits = 0, recursions = 0,
                  fixed_recursions = 0, refused_after_recursion = 0, accepted = 0, rows = 0;
};

struct CensusSampled {   // requests whose read content is in the sample
    std::uint64_t requests = 0, new_keys = 0, repeats = 0;
    // where a repeat's PREVIOUS occurrence of the same key came from
    std::uint64_t repeat_from_other_group = 0, repeat_within_unit = 0, repeat_across_units = 0;
    std::uint64_t repeat_other_fragment = 0, repeat_same_fragment_same_unit = 0, repeat_same_fragment_other_unit = 0;
    std::uint64_t served = 0, missed = 0, repeated_misses = 0, repeated_misses_after_not_retained = 0;
    std::uint64_t kernels = 0, repeated_kernels = 0, repeated_kernels_cross_group = 0;
    std::uint64_t recursions = 0, repeated_recursions = 0, repeated_recursions_cross_group = 0;
    std::uint64_t emission_hits = 0, accepted = 0, rows = 0, dropped = 0;
    // distinct keys requested in this group, and how many requests each received here: 1, 2, 3-4, 5-8, 9-16, 17+
    std::uint64_t keys_seen = 0;
    std::array<std::uint64_t, 6> requests_per_key{};
};

struct VerificationCensusReport {
    bool active = false;
    std::uint32_t period = 0;
    std::uint64_t budget_bytes = 0, peak_bytes = 0;
    bool complete = true;                  // false: the storage budget refused a key at least once
    std::uint64_t sampled_read_contents = 0, distinct_keys = 0, contexts = 0;
    std::array<CensusExact, kCensusGroups> exact{};
    std::array<CensusSampled, kCensusGroups> sampled{};
    // distinct keys by the set of groups that requested them (bit g = group g)
    std::vector<std::pair<std::uint32_t, std::uint64_t>> keys_by_group_mask;
};

// Lifecycle (process-wide; one census at a time). period 0 is off.
void begin_verification_census(std::uint32_t period, std::uint64_t budget_bytes);
VerificationCensusReport end_verification_census();
bool verification_census_active();
void set_verification_census_stage(CensusStage stage);
CensusStage verification_census_stage();
// Tests: every table address collides, so identity rests on the exact comparison alone.
void set_verification_census_hash_collapse(bool collapse);

// The placement routine's identity for this thread, restored on scope exit.
class CensusSiteScope {
public:
    explicit CensusSiteScope(CensusSite site);
    ~CensusSiteScope();
    CensusSiteScope(const CensusSiteScope&) = delete;
    CensusSiteScope& operator=(const CensusSiteScope&) = delete;
private:
    CensusSite previous_;
};

// One discovery unit (one placement call over one text) on this thread: assigns it a fresh unit id. Its exact
// counters are published to the census when the scope closes.
class CensusUnitScope {
public:
    CensusUnitScope();
    ~CensusUnitScope();
    CensusUnitScope(const CensusUnitScope&) = delete;
    CensusUnitScope& operator=(const CensusUnitScope&) = delete;
private:
    bool on_;
};

// Whether a read with this content hash is in the sample (the census must be active).
bool census_sampled(std::uint64_t read_hash_a, std::uint64_t read_hash_b);

// Records one request. `read_hash_*` is the read's existing content hash (its emission key prefix). The read and the
// window are compared byte for byte; `fragment` is the request's fragment index in the run's read list.
void census_record(const CensusContext& context, std::uint64_t read_hash_a, std::uint64_t read_hash_b,
                   const char* read, std::size_t read_length, const char* window, std::size_t window_length,
                   std::uint32_t fragment, const CensusOutcome& outcome);

// The census as a compact key/value table (group, metric, value).
std::string verification_census_table(const VerificationCensusReport& report);

}  // namespace panvar
