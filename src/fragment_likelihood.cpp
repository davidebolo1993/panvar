#include "panvar/fragment_likelihood.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <exception>
#include <functional>
#include <map>
#include <numeric>
#include <memory>
#include <deque>
#include <thread>
#include <set>
#include <tuple>
#include <type_traits>
#include <utility>
#include <array>
#include <cmath>
#include <cstring>
#include <mutex>
#include <condition_variable>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <string_view>
#include <cstdio>
#include <filesystem>

#include "edlib.h"

namespace panvar {
namespace {

constexpr double kNegInf = -std::numeric_limits<double>::infinity();

bool is_acgt(char base) noexcept {
    return base == 'A' || base == 'C' || base == 'G' || base == 'T';
}

char complement(char base) {
    switch (base) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        default:
            throw std::invalid_argument(
                "fragment likelihood is substitution-only over ACGT; a non-ACGT base must refuse "
                "rather than be scored as a mismatch");
    }
}

std::string reverse_complement_read(const std::string& bases);

std::string reverse_complement_read(const std::string& bases) {
    std::string out(bases.size(), 'N');
    for (std::size_t i = 0; i < bases.size(); ++i) {
        out[bases.size() - 1 - i] = complement(bases[i]);
    }
    return out;
}

// A reverse-role mate is compared against the haplotype's forward strand as rc(mate), and is then
// anchored at its RIGHT endpoint by running the recursion over reversed sequences. The reversed form
// of rc(mate) is just its complement, since rc(X) = reverse(complement(X)). Spelling that out here
// keeps the scorer from reversing a string per verification in its hottest loop.
std::string complement_read(const std::string& bases) {
    std::string out(bases.size(), 'N');
    for (std::size_t i = 0; i < bases.size(); ++i) out[i] = complement(bases[i]);
    return out;
}

// One physical origin under construction, while the latent alignments that explain it are summed.
// Shared by the oracle and the accelerated scorer so the two cannot drift in how they collapse.
struct OriginAccumulator {
    double log_latent_mass = kNegInf;
    std::size_t latent_alignments = 0;
    // The best single latent alignment, kept only to describe the origin in diagnostics.
    double best_pair = kNegInf;
    std::size_t reverse_start = 0;
    std::size_t forward_consumed_bp = 0;
    std::size_t reverse_consumed_bp = 0;
    std::size_t forward_mismatches = 0;
    std::size_t reverse_mismatches = 0;

    void add(double pair_mass, std::size_t reverse, std::size_t forward_consumed,
             std::size_t reverse_consumed, std::size_t forward_edits, std::size_t reverse_edits) {
        log_latent_mass = log_add_exp(log_latent_mass, pair_mass);
        ++latent_alignments;
        if (pair_mass <= best_pair) return;
        best_pair = pair_mass;
        reverse_start = reverse;
        forward_consumed_bp = forward_consumed;
        reverse_consumed_bp = reverse_consumed;
        forward_mismatches = forward_edits;
        reverse_mismatches = reverse_edits;
    }
};

} // namespace

std::optional<std::size_t> count_mismatches(
    const std::string& haplotype,
    const std::string& read,
    std::size_t start,
    std::size_t mismatch_limit) {
    if (read.empty()) {
        throw std::invalid_argument("fragment likelihood requires non-empty reads");
    }
    if (start > haplotype.size() || haplotype.size() - start < read.size()) {
        return std::nullopt;
    }
    std::size_t mismatches = 0;
    for (std::size_t i = 0; i < read.size(); ++i) {
        const char reference = haplotype[start + i];
        const char observed = read[i];
        if (!is_acgt(reference) || !is_acgt(observed)) {
            throw std::invalid_argument(
                "fragment likelihood is substitution-only over ACGT; a non-ACGT base must refuse "
                "rather than be scored as a mismatch");
        }
        if (reference != observed && ++mismatches > mismatch_limit) return std::nullopt;
    }
    return mismatches;
}

double log_emission(
    std::size_t read_length,
    std::size_t mismatches,
    double substitution_error) {
    if (mismatches > read_length) {
        throw std::invalid_argument("mismatch count exceeds the read length");
    }
    if (!(substitution_error > 0.0) || !(substitution_error < 0.75)) {
        throw std::invalid_argument("substitution error must be in (0, 0.75)");
    }
    // A substitution lands on one of three other bases, so a specific observed mismatch carries
    // error/3. Folding the /3 in here rather than at the call site keeps the oracle and the
    // accelerated path on the same constant.
    const double matches = static_cast<double>(read_length - mismatches);
    return matches * std::log(1.0 - substitution_error) +
           static_cast<double>(mismatches) * std::log(substitution_error / 3.0);
}

double log_band_mass(std::size_t read_length, std::size_t mismatch_band, double substitution_error) {
    if (!(substitution_error > 0.0) || !(substitution_error < 0.75)) {
        throw std::invalid_argument("substitution error must be in (0, 0.75)");
    }
    if (read_length == 0) throw std::invalid_argument("band mass needs a non-empty read");
    const std::size_t top = std::min(mismatch_band, read_length);
    const double n = static_cast<double>(read_length);
    double total = kNegInf;
    for (std::size_t m = 0; m <= top; ++m) {
        const double k = static_cast<double>(m);
        const double log_choose = std::lgamma(n + 1.0) - std::lgamma(k + 1.0) - std::lgamma(n - k + 1.0);
        total = log_add_exp(total, log_choose + k * std::log(substitution_error) +
                                       (n - k) * std::log1p(-substitution_error));
    }
    return total;
}

double log_truncated_emission(
    std::size_t read_length,
    std::size_t mismatches,
    double substitution_error,
    std::size_t mismatch_band) {
    if (mismatches > mismatch_band) return kNegInf;
    return log_emission(read_length, mismatches, substitution_error) -
           log_band_mass(read_length, mismatch_band, substitution_error);
}

// -------------------------------------------------------------------------------------------
// Indel-aware emission: a banded affine-gap pair-HMM, conditioned on the observed read length.
// -------------------------------------------------------------------------------------------

namespace {

enum GapState { kMatch = 0, kInsert = 1, kDelete = 2, kStates = 3 };

// The recursion, shared by the emission and its normaliser. `emitting` false replaces every
// emission with one, which is exactly the sum over all 4^L reads and therefore the shape mass.
// Cells are indexed [read position][displacement][state][edits used].
// `available` is how many haplotype bases may be consumed from `start`. It bounds BOTH modes: the
// numerator must not walk off the haplotype, and the normaliser must not count shapes that would.
// Bounding only the emitting pass is what leaked mass near the haplotype ends -- shapes dropped from
// the numerator were still dividing the denominator, so the emission summed to less than one there.
//
// THIS IS THE REFERENCE. It is kept exactly as committed and is no longer on the production path;
// gap_recursion_scaled below computes the same quantities far more cheaply, and the oracle suite
// compares the two component by component. Do not optimise this one -- its value is that it is the
// implementation the model was gated on.
double gap_recursion_logspace(
    const std::string* haplotype,
    const std::string* read,
    std::size_t start,
    std::size_t read_length,
    std::size_t available,
    std::size_t band,
    const LibraryModel& library,
    bool emitting,
    ComponentList* components) {
    const std::size_t delta = library.max_gap_bp;
    const std::size_t width = 2 * delta + 1;
    const std::size_t budget = band + 1;                 // edits 0..band
    const double log_match = std::log(1.0 - library.substitution_error);
    const double log_mismatch = std::log(library.substitution_error / 3.0);
    const double log_insert_emission = std::log(0.25);
    const double log_continue = library.gap_open > 0.0
        ? std::log(1.0 - 2.0 * library.gap_open) : 0.0;
    const double log_open = library.gap_open > 0.0 ? std::log(library.gap_open) : kNegInf;
    const double log_extend = library.gap_extend > 0.0 ? std::log(library.gap_extend) : kNegInf;
    const double log_close = library.gap_extend < 1.0 ? std::log(1.0 - library.gap_extend) : kNegInf;

    const std::size_t cells = width * kStates * budget;
    std::vector<double> previous(cells, kNegInf);
    std::vector<double> current(cells, kNegInf);

    const auto at = [&](std::size_t displacement, std::size_t state, std::size_t edits) {
        return (displacement * kStates + state) * budget + edits;
    };
    // Displacement k is stored offset by delta, so index 0 means k = -delta.
    previous[at(delta, kMatch, 0)] = 0.0;

    for (std::size_t i = 1; i <= read_length; ++i) {
        std::fill(current.begin(), current.end(), kNegInf);
        for (std::size_t d = 0; d < width; ++d) {
            const long k = static_cast<long>(d) - static_cast<long>(delta);
            // Haplotype bases consumed before this read base under M: i - 1 + k.
            const long consumed_before = static_cast<long>(i) - 1 + k;
            for (std::size_t e = 0; e < budget; ++e) {
                // --- M: consume one read base and one haplotype base -------------------------
                if (consumed_before >= 0) {
                    const std::size_t index = start + static_cast<std::size_t>(consumed_before);
                    // The match step, as (weight, edit cost) branches.
                    //
                    // Emitting, the read base is known, so there is exactly one branch. Summing over
                    // ALL 4^L reads there are two, and they must both be walked: a match position
                    // carries (1-e) over the one base that matches, and a mismatch carries 3*(e/3) = e
                    // over the three that do not, at a cost of one edit. Collapsing that to a single
                    // zero-cost branch -- which is what computing the normaliser with "no emission"
                    // does -- makes Z count only the all-match path, and the emission then reduces to
                    // the UNTRUNCATED substitution model instead of the declared truncated one.
                    double weights[2] = {0.0, 0.0};
                    std::size_t costs[2] = {0, 1};
                    std::size_t branches = 0;
                    // The consumption bound applies whether or not we are emitting.
                    bool usable = static_cast<std::size_t>(consumed_before) < available;
                    if (emitting && usable) {
                        {
                            const char reference = (*haplotype)[index];
                            const char observed = (*read)[i - 1];
                            if (!is_acgt(reference) || !is_acgt(observed)) {
                                throw std::invalid_argument(
                                    "gapped emission is defined over ACGT; a non-ACGT base must "
                                    "refuse rather than be scored as a mismatch");
                            }
                            const bool same = reference == observed;
                            weights[0] = same ? log_match : log_mismatch;
                            costs[0] = same ? 0 : 1;
                            branches = 1;
                        }
                    } else {
                        weights[0] = std::log(1.0 - library.substitution_error);
                        costs[0] = 0;
                        weights[1] = std::log(library.substitution_error);
                        costs[1] = 1;
                        branches = 2;
                    }
                    for (std::size_t b = 0; usable && b < branches; ++b) {
                        const double emission = weights[b];
                        const std::size_t cost = costs[b];
                        if (e < cost) continue;
                        const std::size_t from = e - cost;
                        const double arrivals[kStates] = {
                            previous[at(d, kMatch, from)] + log_continue,
                            previous[at(d, kInsert, from)] + log_close,
                            previous[at(d, kDelete, from)] + log_close,
                        };
                        double sum = kNegInf;
                        for (std::size_t s = 0; s < kStates; ++s) {
                            sum = log_add_exp(sum, arrivals[s]);
                        }
                        const std::size_t cell = at(d, kMatch, e);
                        current[cell] = log_add_exp(current[cell], sum + emission);
                    }
                }
                // --- I: consume one read base only, displacement falls by one ----------------
                if (d + 1 < width && e >= 1) {
                    const double emission = emitting ? log_insert_emission : 0.0;
                    const std::size_t from = e - 1;
                    const double opened = previous[at(d + 1, kMatch, from)] + log_open;
                    const double extended = previous[at(d + 1, kInsert, from)] + log_extend;
                    const std::size_t cell = at(d, kInsert, e);
                    current[cell] = log_add_exp(current[cell], log_add_exp(opened, extended) + emission);
                }
            }
        }
        // --- D: consume haplotype only. No read base advances, so it is resolved within the
        // current read position, sweeping displacement upward because D raises k by one.
        for (std::size_t d = 1; d < width; ++d) {
            const long k = static_cast<long>(d) - static_cast<long>(delta);
            const long consumed_before = static_cast<long>(i) + k - 1;
            if (consumed_before < 0) continue;
            if (static_cast<std::size_t>(consumed_before) >= available) continue;
            for (std::size_t e = 1; e < budget; ++e) {
                const std::size_t from = e - 1;
                const double opened = current[at(d - 1, kMatch, from)] + log_open;
                const double extended = current[at(d - 1, kDelete, from)] + log_extend;
                const std::size_t cell = at(d, kDelete, e);
                current[cell] = log_add_exp(current[cell], log_add_exp(opened, extended));
            }
        }
        previous.swap(current);
    }

    double total = kNegInf;
    // The displacement band bounds consumption from BELOW as well as above. Without the floor an
    // alignment may terminate having consumed almost no haplotype -- a read absorbed entirely by
    // insertions -- which is not a placement in any useful sense: near the end of a haplotype every
    // start would then admit the read, and those degenerate origins put a whole extra unit of mass
    // into the exposure sum. With it, a read that cannot fit yields no placement, exactly as the
    // substitution-only model behaved when count_mismatches ran off the end.
    const std::size_t floor_bp = read_length - std::min(read_length, delta);
    // Each terminal displacement is one CONSUMED-SPAN CLASS, and they are kept apart here rather than
    // summed into one number. They have to be: every class ends the mate at a different haplotype
    // coordinate, so each implies a different fragment geometry and a different insert probability.
    // Collapsing them and handing the total to the best class's span weights the other classes' mass
    // with an insert probability that is not theirs -- which is what made the signal density fail to
    // integrate to the exposure. The caller applies the geometry to each class separately.
    for (std::size_t d = 0; d < width; ++d) {
        const long k = static_cast<long>(d) - static_cast<long>(delta);
        const long consumed = static_cast<long>(read_length) + k;
        if (consumed < 0 || static_cast<std::size_t>(consumed) > available) continue;
        if (static_cast<std::size_t>(consumed) < floor_bp) continue;
        double span = kNegInf;
        // Fewest edits reaching this span. A cell is INDEXED by the edits spent to reach it, so the
        // lowest finite e in the class is that minimum exactly -- there is nothing to track alongside
        // the recursion, and the arrays that used to do it could only restate their own index.
        std::size_t span_edits = budget;
        for (std::size_t s = 0; s < kStates; ++s) {
            for (std::size_t e = 0; e < budget; ++e) {
                const double cell = previous[at(d, s, e)];
                if (cell == kNegInf) continue;
                span = log_add_exp(span, cell);
                if (e < span_edits) span_edits = e;
            }
        }
        if (span == kNegInf) continue;
        total = log_add_exp(total, span);
        if (components != nullptr) {
            components->push_back({static_cast<std::size_t>(consumed), span, span_edits});
        }
    }
    return total;
}

// ---------------------------------------------------------------------------------------------
// Optimisation only. Everything below computes exactly what the reference above computes: same
// probability model, same gap parameters, same band, same normalisation, same physical origins,
// same resource semantics. What changes is arithmetic and memory, nothing declared.
// ---------------------------------------------------------------------------------------------

// Per-thread DP scratch. The recursion used to allocate two vectors on every call, and there are two
// calls per placement, which at C4 scale is tens of millions of allocations for buffers of 147
// doubles. Sized once per worker and reused; growth is monotone, so a worker settles on its first
// placement and never allocates again.
struct GapScratch {
    std::vector<double> first;
    std::vector<double> second;

    void ensure(std::size_t cells) {
        if (first.size() < cells) {
            first.assign(cells, 0.0);
            second.assign(cells, 0.0);
        }
    }
};

GapScratch& gap_scratch() {
    static thread_local GapScratch scratch;
    return scratch;
}

struct EmissionWork {
    std::uint64_t dp_evaluations = 0;
    std::uint64_t prefilter_rejections = 0;
    // Observational only, and deliberately not atomic: at 1.19 billion prefilters a shared
    // fetch_add is contention with no decision attached to it.
    std::uint64_t executed_prefilters = 0;
    std::uint64_t dp_cache_requests = 0;
    std::uint64_t dp_cache_hits = 0;
    std::uint64_t dp_cache_misses = 0;
    std::uint64_t executed_dps = 0;
    std::uint64_t cache_lookup_nanos = 0;
};

EmissionWork& emission_work() {
    static thread_local EmissionWork work;
    return work;
}

// The process-wide worker gauge. A thread holds a slot for exactly as long as it does scoring work.
std::atomic<std::size_t>& active_worker_count() {
    static std::atomic<std::size_t> count{0};
    return count;
}
std::atomic<std::size_t>& active_worker_peak() {
    static std::atomic<std::size_t> peak{0};
    return peak;
}
struct ActiveWorker {
    ActiveWorker() {
        const std::size_t now = active_worker_count().fetch_add(1, std::memory_order_relaxed) + 1;
        std::size_t peak = active_worker_peak().load(std::memory_order_relaxed);
        while (now > peak && !active_worker_peak().compare_exchange_weak(peak, now, std::memory_order_relaxed)) {
        }
    }
    ~ActiveWorker() { active_worker_count().fetch_sub(1, std::memory_order_relaxed); }
    ActiveWorker(const ActiveWorker&) = delete;
    ActiveWorker& operator=(const ActiveWorker&) = delete;
};

// MEASUREMENT ONLY: per-thread time inside the feasibility prefilter and inside gapped_emission. Off
// unless enabled; when off the prefilter pays one relaxed load.
struct EmissionTimingDiag {
    std::uint64_t prefilter_ns = 0;
    std::uint64_t emission_ns = 0;
};
std::atomic<bool>& emission_timing_active() {
    static std::atomic<bool> active{false};
    return active;
}
EmissionTimingDiag& emission_timing() {
    static thread_local EmissionTimingDiag timing;
    return timing;
}
std::uint64_t nanos_since(std::chrono::steady_clock::time_point from) {
    return static_cast<std::uint64_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - from).count());
}

// ---- verification-stream census (opt-in; see EmissionCensus in the header) --------------------
struct CensusKey {
    std::uint64_t a = 0;
    std::uint64_t b = 0;
    bool operator==(const CensusKey& other) const noexcept {
        return a == other.a && b == other.b;
    }
};
struct CensusKeyHash {
    std::size_t operator()(const CensusKey& key) const noexcept {
        return static_cast<std::size_t>(key.a ^ (key.b * 0x9e3779b97f4a7c15ULL));
    }
};
struct CensusEntry {
    std::uint32_t requests = 0;
    std::uint16_t value_bytes = 0;
    bool executed_dp = false;
};

std::mutex& census_mutex() { static std::mutex guard; return guard; }
std::atomic<bool>& census_active() { static std::atomic<bool> active{false}; return active; }
std::unordered_map<CensusKey, CensusEntry, CensusKeyHash>& census_map() {
    static std::unordered_map<CensusKey, CensusEntry, CensusKeyHash> table;
    return table;
}
std::uint64_t& census_key_bytes() { static std::uint64_t total = 0; return total; }

// 128 bits over (read bytes, separator, window bytes). Two independent mixers, so a collision needs
// both to agree; at ~4e7 keys that probability is negligible, and it would under-count distinctness.
CensusKey census_key_of(const std::string& read, const char* window, std::size_t window_length) {
    std::uint64_t a = 1469598103934665603ULL;
    std::uint64_t b = 0x2545f4914f6cdd1dULL;
    const auto mix = [&](unsigned char value) {
        a = (a ^ value) * 1099511628211ULL;
        b = (b + value) * 0x9e3779b97f4a7c15ULL;
        b ^= b >> 29;
    };
    for (const char c : read) mix(static_cast<unsigned char>(c));
    mix(0xff);                                    // separator: (read, window) is not (read||window)
    for (std::size_t i = 0; i < window_length; ++i) mix(static_cast<unsigned char>(window[i]));
    return CensusKey{a, b};
}

// ---- exact emission cache (see EmissionCacheCounters in the header) --------------------------
//
// Keys hold POINTERS into the caller's form and haplotype buffers rather than copies: storing the
// ~203 key bytes per entry would cost 1.96GB at 30x, against ~48 by reference. The referenced
// buffers must outlive the scoring pass, which is why the reversed haplotypes are precomputed once
// below rather than rebuilt per haplotype iteration.
// The key OWNS its bytes. It used to hold pointers into the caller's read and haplotype buffers,
// which is why anything transient had to be kept out of the table. `content` is the read followed by
// the window, and `read_length` records where one ends and the other begins -- so two different
// (read, window) pairs cannot collide by concatenating to the same bytes, and identity stays exact.
//
// Hashes select a bucket and nothing else. Equality compares the COMPLETE content, so a collision
// costs a comparison and never a wrong answer. Nothing here interprets a base, so the key admits
// whatever the emission admits; today require_acgt_window refuses anything but ACGT before a key is
// ever formed, and the key is correct either way.
// A lookup's key, as a VIEW. Non-owning, stack-only, valid for exactly the duration of the lookup:
// the pointers are the caller's read and the haplotype window it is being compared against. Nothing
// here allocates, and nothing here is retained.
//
// The previous revision owned its bytes in a std::string, which fixed the dangling-pointer problem
// and introduced a worse one: the call site copy-assigned a fresh key per request, so 8.2M lookups a
// round became 8.2M heap allocations. Ownership belongs to the CACHE, not to the key.
struct KeyView {
    std::uint64_t h0 = 0;
    std::uint64_t h1 = 0;
    const char* read = nullptr;
    const char* window = nullptr;
    std::uint32_t read_length = 0;
    std::uint32_t window_length = 0;
};

enum class CacheOutcome : std::uint8_t { Accepted, RefusedAfterRecursion };

struct CacheValue {
    double log_emission = 0.0;
    std::uint32_t consumed_bp = 0;
    std::uint32_t edits = 0;
    std::uint32_t component_offset = 0;
    std::uint16_t component_count = 0;
    // Both store paths set this explicitly; the default is the conservative one, since a value that
    // somehow escaped without being set describes a recursion that ran and found nothing rather than
    // an acceptance nobody computed. The distinction is carried in the VALUE because the uncached
    // path counts a prefilter rejection and a post-recursion refusal differently, and a logical
    // counter that moved when the cache was armed would make cached and uncached runs incomparable.
    CacheOutcome outcome = CacheOutcome::RefusedAfterRecursion;
};

// A retained entry. Fixed size, stored inline in an open-addressed table, and carrying OFFSETS into
// the content arena rather than pointers -- so a transient haplotype may be destroyed while every
// entry it produced keeps serving, and so arena growth cannot invalidate anything.
//
// `state` is an explicit occupancy flag rather than a hash sentinel, because the forced-collision
// test drives every key to h0 == h1 == 0 and a zero hash must remain a legal key.
struct ArenaEntry {
    std::uint64_t h0 = 0;
    std::uint64_t h1 = 0;
    std::uint32_t read_length = 0;
    std::uint32_t window_length = 0;
    std::uint32_t slab = 0;
    std::uint32_t offset = 0;
    CacheValue value;
    std::uint8_t state = 0;   // 0 empty, 1 occupied
};

// Slabs are append-only and never reallocated, which is what makes an offset permanent. A request
// larger than a slab gets one of its own rather than being split.
constexpr std::size_t kSlabBytes = 1u << 20;
constexpr std::size_t kNoSlot = static_cast<std::size_t>(-1);

struct CacheShard {
    std::mutex guard;
    std::vector<ArenaEntry> table;                 // power-of-two capacity, linear probing
    std::size_t occupied = 0;
    std::vector<std::unique_ptr<char[]>> slabs;    // stable: never moved, never resized
    std::size_t slab_used = 0;
    std::vector<GappedComponent> arena;            // components, addressed by index
};

// Measured, not assumed: raising this to 1024 moved a repair round from 82.20 s to 82.73 s, so the
// 8.2M lookups a round performs are NOT queueing behind these mutexes and finer sharding buys
// nothing. Left at the value the 560M-lookup one-shot pass was sized for.
constexpr std::size_t kCacheShards = 64;

struct EmissionCache {
    std::array<CacheShard, kCacheShards> shards;
    std::atomic<bool> active{false};
    // Governs EVERY retained entry: its table slots, its arena bytes and its components. Declining
    // one costs a recomputation, never a wrong answer.
    std::atomic<std::size_t> budget_bytes{0};
    std::atomic<std::size_t> budgeted_bytes_used{0};
    std::atomic<std::uint64_t> entries{0};
    std::atomic<std::uint64_t> payload_bytes{0};
    std::atomic<std::uint64_t> entries_not_retained{0};
    // I2-K: hits served by entries inherited from the union pass (state 2), counted under the shard lock.
    std::atomic<std::uint64_t> inherited_hits{0};
    // Nanoseconds inside lookups, folded in per haplotype rather than per lookup: at 8.2M lookups a
    // round a shared counter touched per request would be contention with no decision attached.
    std::atomic<std::uint64_t> lookup_nanos{0};
};

EmissionCache& emission_cache() {
    static EmissionCache cache;
    return cache;
}

std::atomic<bool>& cache_hash_collapse() {
    static std::atomic<bool> collapse{false};
    return collapse;
}

// Builds the view on the stack. No allocation, no copy, nothing retained.
KeyView key_view_of(const std::string& read, const char* window, std::size_t window_length) {
    std::uint64_t a = 1469598103934665603ULL;
    std::uint64_t b = 0x2545f4914f6cdd1dULL;
    const auto mix = [&](unsigned char value) {
        a = (a ^ value) * 1099511628211ULL;
        b = (b + value) * 0x9e3779b97f4a7c15ULL;
        b ^= b >> 29;
    };
    for (const char c : read) mix(static_cast<unsigned char>(c));
    mix(0xff);
    for (std::size_t i = 0; i < window_length; ++i) mix(static_cast<unsigned char>(window[i]));
    KeyView key;
    // Under the test collapse every key hashes alike, so the table degenerates to one probe chain
    // and byte comparison alone separates them. The answers may not change.
    const bool collapse = cache_hash_collapse().load(std::memory_order_relaxed);
    key.h0 = collapse ? 0 : a;
    key.h1 = collapse ? 0 : b;
    key.read = read.data();
    key.read_length = static_cast<std::uint32_t>(read.size());
    key.window = window;
    key.window_length = static_cast<std::uint32_t>(window_length);
    return key;
}

// The key hash is byte-serial, so its state after "read bytes + 0xff" is a function of the read alone: carried per
// prepared form, the window's bytes resume from it and give exactly key_view_of's h0 and h1. Identity is still
// decided by complete byte comparison; only the hashing of the read side is not repeated.
struct KeyPrefix {
    std::uint64_t a = 0;
    std::uint64_t b = 0;
};
KeyPrefix key_prefix_of(const std::string& read) {
    KeyPrefix prefix{1469598103934665603ULL, 0x2545f4914f6cdd1dULL};
    const auto mix = [&](unsigned char value) {
        prefix.a = (prefix.a ^ value) * 1099511628211ULL;
        prefix.b = (prefix.b + value) * 0x9e3779b97f4a7c15ULL;
        prefix.b ^= prefix.b >> 29;
    };
    for (const char c : read) mix(static_cast<unsigned char>(c));
    mix(0xff);
    return prefix;
}
KeyView key_view_from_prefix(const KeyPrefix& prefix, const std::string& read, const char* window,
                             std::size_t window_length) {
    std::uint64_t a = prefix.a;
    std::uint64_t b = prefix.b;
    for (std::size_t i = 0; i < window_length; ++i) {
        const unsigned char value = static_cast<unsigned char>(window[i]);
        a = (a ^ value) * 1099511628211ULL;
        b = (b + value) * 0x9e3779b97f4a7c15ULL;
        b ^= b >> 29;
    }
    KeyView key;
    const bool collapse = cache_hash_collapse().load(std::memory_order_relaxed);
    key.h0 = collapse ? 0 : a;
    key.h1 = collapse ? 0 : b;
    key.read = read.data();
    key.read_length = static_cast<std::uint32_t>(read.size());
    key.window = window;
    key.window_length = static_cast<std::uint32_t>(window_length);
    return key;
}

std::size_t slot_hash(std::uint64_t h0, std::uint64_t h1) {
    return static_cast<std::size_t>(h0 ^ (h1 * 0x9e3779b97f4a7c15ULL));
}

// Content bytes are 8-aligned inside a slab, so the padding a retained entry occupies is charged
// rather than ignored.
std::size_t aligned_content_bytes(std::size_t read_length, std::size_t window_length) {
    const std::size_t raw = read_length + window_length;
    return (raw + 7u) & ~static_cast<std::size_t>(7u);
}

// What a retained entry declares against the budget: two table slots, because the table is held at
// or below half load so every entry reserves the slack that keeps probing short; its aligned content
// bytes; and the slab pointer's own share. Components are added by the fit path.
std::size_t arena_entry_bytes(std::size_t read_length, std::size_t window_length) {
    return 2 * sizeof(ArenaEntry) + aligned_content_bytes(read_length, window_length) +
           sizeof(std::unique_ptr<char[]>);
}

const char* entry_bytes(const CacheShard& shard, const ArenaEntry& entry) {
    return shard.slabs[entry.slab].get() + entry.offset;
}

// Exact identity: lengths first, then the complete retained content. The hash only chose where to
// start probing.
bool entry_matches(const CacheShard& shard, const ArenaEntry& entry, const KeyView& key) {
    if (entry.h0 != key.h0 || entry.h1 != key.h1) return false;
    if (entry.read_length != key.read_length || entry.window_length != key.window_length) return false;
    const char* bytes = entry_bytes(shard, entry);
    return std::memcmp(bytes, key.read, key.read_length) == 0 &&
           std::memcmp(bytes + key.read_length, key.window, key.window_length) == 0;
}

std::size_t find_occupied_slot(const CacheShard& shard, const KeyView& key) {
    if (shard.table.empty()) return kNoSlot;
    const std::size_t mask = shard.table.size() - 1;
    std::size_t at = slot_hash(key.h0, key.h1) & mask;
    for (std::size_t probes = 0; probes <= mask; ++probes) {
        const ArenaEntry& entry = shard.table[at];
        if (entry.state == 0) return kNoSlot;
        if (entry_matches(shard, entry, key)) return at;
        at = (at + 1) & mask;
    }
    return kNoSlot;
}

void grow_table(CacheShard& shard) {
    const std::size_t wanted = shard.table.empty() ? 1024 : shard.table.size() * 2;
    std::vector<ArenaEntry> next(wanted);
    const std::size_t mask = wanted - 1;
    for (const ArenaEntry& entry : shard.table) {
        if (entry.state == 0) continue;
        std::size_t at = slot_hash(entry.h0, entry.h1) & mask;
        while (next[at].state != 0) at = (at + 1) & mask;
        next[at] = entry;
    }
    shard.table.swap(next);
}

// Copies the key's bytes into append-only slab storage exactly once, on the insertion path only.
void store_content(CacheShard& shard, const KeyView& key, ArenaEntry& entry) {
    const std::size_t need = aligned_content_bytes(key.read_length, key.window_length);
    if (shard.slabs.empty() || shard.slab_used + need > kSlabBytes) {
        shard.slabs.push_back(std::make_unique<char[]>(std::max(kSlabBytes, need)));
        shard.slab_used = 0;
    }
    char* at = shard.slabs.back().get() + shard.slab_used;
    std::memcpy(at, key.read, key.read_length);
    std::memcpy(at + key.read_length, key.window, key.window_length);
    entry.slab = static_cast<std::uint32_t>(shard.slabs.size() - 1);
    entry.offset = static_cast<std::uint32_t>(shard.slab_used);
    shard.slab_used += need;
}

enum class CacheProbe { Disabled, Miss, HitAccepted, HitRefusedAfterRecursion };

bool same_key(const KeyView& a, const KeyView& b) {
    return a.h0 == b.h0 && a.h1 == b.h1 && a.read_length == b.read_length &&
           a.window_length == b.window_length && std::memcmp(a.read, b.read, a.read_length) == 0 &&
           std::memcmp(a.window, b.window, a.window_length) == 0;
}

// A parallel window enumeration freezes the shared emission cache: workers only read it. A key that
// misses there is resolved from the worker's own earlier result, or computed and recorded; after the
// join the caller replays every record in single-thread order, inserting into the shared cache exactly
// where one thread would have. That keeps hits, misses, executed recursions and retention decisions
// the one-thread values, whatever the schedule.
struct DeferredEmissions {
    struct Record {
        KeyView key;
        bool refused = false;
        GappedFit fit;
    };
    std::vector<Record> records;   // every shared-cache miss, in this worker's request order
    std::unordered_multimap<std::uint64_t, std::size_t> first_of;   // h0 -> first computing record
    const Record* find(const KeyView& key) const {
        const auto range = first_of.equal_range(key.h0);
        for (auto it = range.first; it != range.second; ++it) {
            if (same_key(records[it->second].key, key)) return &records[it->second];
        }
        return nullptr;
    }
    void add(const KeyView& key, bool refused, const GappedFit* fit) {
        Record record;
        record.key = key;
        record.refused = refused;
        if (fit != nullptr) record.fit = *fit;
        first_of.emplace(key.h0, records.size());
        records.push_back(std::move(record));
    }
    void clear() {
        records.clear();
        first_of.clear();
    }
};
DeferredEmissions*& deferred_emissions() {
    static thread_local DeferredEmissions* deferred = nullptr;
    return deferred;
}

// The hit path: one mutex, one probe chain, one memcmp pair, and a component copy into a buffer the
// caller already owns. No allocation.
// Lookup timing is diagnostic work -- two clock reads per lookup, measured at ~4.6% of union placement -- so it runs
// only when asked for; otherwise no clock is read and the time is reported as not measured, never as zero.
std::atomic<bool>& emission_lookup_timing() {
    static std::atomic<bool> enabled{false};
    return enabled;
}

CacheProbe cache_lookup(const KeyView& key, GappedFit* out, EmissionWork& work) {
    EmissionCache& cache = emission_cache();
    if (!cache.active.load(std::memory_order_relaxed)) return CacheProbe::Disabled;
    const bool timed = emission_lookup_timing().load(std::memory_order_relaxed);
    const auto mark = timed ? std::chrono::steady_clock::now() : std::chrono::steady_clock::time_point{};
    CacheShard& shard = cache.shards[key.h0 % kCacheShards];
    CacheProbe answer = CacheProbe::Miss;
    {
        const std::lock_guard<std::mutex> guard(shard.guard);
        const std::size_t at = find_occupied_slot(shard, key);
        if (at == kNoSlot) {
            ++work.dp_cache_misses;
        } else {
            ++work.dp_cache_hits;
            if (shard.table[at].state == 2) cache.inherited_hits.fetch_add(1, std::memory_order_relaxed);   // I2-K
            const CacheValue& value = shard.table[at].value;
            if (value.outcome == CacheOutcome::RefusedAfterRecursion) {
                answer = CacheProbe::HitRefusedAfterRecursion;
            } else {
                out->log_emission = value.log_emission;
                out->consumed_bp = value.consumed_bp;
                out->edits = value.edits;
                // assign into the caller's existing capacity; the arena may be appended to by a
                // concurrent insert, so it is read under the lock.
                out->components.assign(shard.arena.begin() + value.component_offset,
                                       shard.arena.begin() + value.component_offset +
                                           value.component_count);
                answer = CacheProbe::HitAccepted;
            }
        }
    }
    if (timed) {
        work.cache_lookup_nanos += static_cast<std::uint64_t>(
            std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - mark).count());
    }
    return answer;
}

// MEASUREMENT ONLY: the hit path's parts, without the clock reads -- with or without the shard lock, with or
// without the component copy. Single-threaded replay only (without the lock it races any writer).
CacheProbe cache_probe_measured(const KeyView& key, GappedFit* out, bool lock, bool copy) {
    EmissionCache& cache = emission_cache();
    CacheShard& shard = cache.shards[key.h0 % kCacheShards];
    std::unique_lock<std::mutex> guard(shard.guard, std::defer_lock);
    if (lock) guard.lock();
    const std::size_t at = find_occupied_slot(shard, key);
    if (at == kNoSlot) return CacheProbe::Miss;
    const CacheValue& value = shard.table[at].value;
    if (value.outcome == CacheOutcome::RefusedAfterRecursion) return CacheProbe::HitRefusedAfterRecursion;
    if (copy) {
        out->log_emission = value.log_emission;
        out->consumed_bp = value.consumed_bp;
        out->edits = value.edits;
        out->components.assign(shard.arena.begin() + value.component_offset,
                               shard.arena.begin() + value.component_offset + value.component_count);
    }
    return CacheProbe::HitAccepted;
}

// Shared by both store paths: charge the budget, take the slot, copy the bytes once.
bool insert_entry(const KeyView& key, const CacheValue& value, std::size_t extra_bytes) {
    EmissionCache& cache = emission_cache();
    const std::size_t bytes = arena_entry_bytes(key.read_length, key.window_length) + extra_bytes;
    const std::size_t used = cache.budgeted_bytes_used.fetch_add(bytes, std::memory_order_relaxed) + bytes;
    if (used > cache.budget_bytes.load(std::memory_order_relaxed)) {
        cache.budgeted_bytes_used.fetch_sub(bytes, std::memory_order_relaxed);
        cache.entries_not_retained.fetch_add(1, std::memory_order_relaxed);
        return false;
    }
    CacheShard& shard = cache.shards[key.h0 % kCacheShards];
    const std::lock_guard<std::mutex> guard(shard.guard);
    if (find_occupied_slot(shard, key) != kNoSlot) {
        cache.budgeted_bytes_used.fetch_sub(bytes, std::memory_order_relaxed);
        return false;   // another worker got there first; the budget is not spent twice
    }
    if (shard.table.empty() || (shard.occupied + 1) * 2 > shard.table.size()) grow_table(shard);
    const std::size_t mask = shard.table.size() - 1;
    std::size_t at = slot_hash(key.h0, key.h1) & mask;
    while (shard.table[at].state != 0) at = (at + 1) & mask;
    ArenaEntry& entry = shard.table[at];
    entry.h0 = key.h0;
    entry.h1 = key.h1;
    entry.read_length = key.read_length;
    entry.window_length = key.window_length;
    entry.value = value;
    store_content(shard, key, entry);
    entry.state = 1;
    ++shard.occupied;
    cache.entries.fetch_add(1, std::memory_order_relaxed);
    cache.payload_bytes.fetch_add(bytes, std::memory_order_relaxed);
    return true;
}

bool cache_contains(const KeyView& key) {
    EmissionCache& cache = emission_cache();
    CacheShard& shard = cache.shards[key.h0 % kCacheShards];
    const std::lock_guard<std::mutex> guard(shard.guard);
    return find_occupied_slot(shard, key) != kNoSlot;
}

void cache_store_refusal(const KeyView& key, CacheOutcome outcome) {
    if (!emission_cache().active.load(std::memory_order_relaxed)) return;
    CacheValue value;
    value.outcome = outcome;
    (void)insert_entry(key, value, 0);
}

void cache_store_fit(const KeyView& key, const GappedFit& fit) {
    EmissionCache& cache = emission_cache();
    if (!cache.active.load(std::memory_order_relaxed)) return;
    const std::size_t extra = fit.components.size() * sizeof(GappedComponent);
    // The components must land in the arena under the same lock that takes the slot, so the offset an
    // entry records is the offset its components actually occupy.
    const std::size_t bytes = arena_entry_bytes(key.read_length, key.window_length) + extra;
    const std::size_t used = cache.budgeted_bytes_used.fetch_add(bytes, std::memory_order_relaxed) + bytes;
    if (used > cache.budget_bytes.load(std::memory_order_relaxed)) {
        cache.budgeted_bytes_used.fetch_sub(bytes, std::memory_order_relaxed);
        cache.entries_not_retained.fetch_add(1, std::memory_order_relaxed);
        return;
    }
    CacheShard& shard = cache.shards[key.h0 % kCacheShards];
    const std::lock_guard<std::mutex> guard(shard.guard);
    if (find_occupied_slot(shard, key) != kNoSlot) {
        cache.budgeted_bytes_used.fetch_sub(bytes, std::memory_order_relaxed);
        return;
    }
    if (shard.table.empty() || (shard.occupied + 1) * 2 > shard.table.size()) grow_table(shard);
    const std::size_t mask = shard.table.size() - 1;
    std::size_t at = slot_hash(key.h0, key.h1) & mask;
    while (shard.table[at].state != 0) at = (at + 1) & mask;
    ArenaEntry& entry = shard.table[at];
    entry.h0 = key.h0;
    entry.h1 = key.h1;
    entry.read_length = key.read_length;
    entry.window_length = key.window_length;
    entry.value.outcome = CacheOutcome::Accepted;
    entry.value.log_emission = fit.log_emission;
    entry.value.consumed_bp = static_cast<std::uint32_t>(fit.consumed_bp);
    entry.value.edits = static_cast<std::uint32_t>(fit.edits);
    entry.value.component_offset = static_cast<std::uint32_t>(shard.arena.size());
    entry.value.component_count = static_cast<std::uint16_t>(fit.components.size());
    shard.arena.insert(shard.arena.end(), fit.components.begin(), fit.components.end());
    store_content(shard, key, entry);
    entry.state = 1;
    ++shard.occupied;
    cache.entries.fetch_add(1, std::memory_order_relaxed);
    cache.payload_bytes.fetch_add(bytes, std::memory_order_relaxed);
}

void census_record(const std::string& read, const char* window, std::size_t window_length,
                   bool executed_dp, std::size_t value_bytes) {
    if (!census_active().load(std::memory_order_relaxed)) return;
    const CensusKey key = census_key_of(read, window, window_length);
    const std::lock_guard<std::mutex> guard(census_mutex());
    CensusEntry& entry = census_map()[key];
    if (entry.requests == 0) {
        entry.executed_dp = executed_dp;
        entry.value_bytes = static_cast<std::uint16_t>(value_bytes);
        census_key_bytes() += read.size() + window_length;
    }
    ++entry.requests;
}

// The contract REFUSES a non-ACGT base rather than scoring it as a mismatch, and the prefilter must
// not quietly turn that refusal into an absence -- edlib would simply find an N-bearing read
// infeasible and the placement would vanish instead of throwing. So the window the recursion can
// touch is validated before anything else looks at it, which is also what lets the hot loop below
// stop testing every base it visits.
void require_acgt_window(const std::string& haplotype, std::size_t start, std::size_t window,
                         const std::string& read) {
    for (const char base : read) {
        if (!is_acgt(base)) {
            throw std::invalid_argument(
                "gapped emission is defined over ACGT; a non-ACGT base must refuse rather than be "
                "scored as a mismatch");
        }
    }
    for (std::size_t i = 0; i < window; ++i) {
        if (!is_acgt(haplotype[start + i])) {
            throw std::invalid_argument(
                "gapped emission is defined over ACGT; a non-ACGT base must refuse rather than be "
                "scored as a mismatch");
        }
    }
}

// The EXACT feasibility prefilter.
//
// An admissible alignment consumes read_length +/- max_gap_bp haplotype bases from the anchor and
// spends at most `band` edits, counting substitutions and gap bases alike -- which is precisely
// unit-cost Levenshtein. edlib's SHW (prefix) mode fits the WHOLE query from target[0] with a free
// gap at the target's end, so it minimises edits over EVERY target prefix. Those prefixes are a
// superset of the consumed spans this model admits, so the SHW minimum can only be smaller than the
// constrained one. If it already exceeds the band, no admissible alignment exists and the mass is
// exactly zero.
//
// Therefore: NO false negatives, by construction rather than by tuning. False positives are allowed
// and expected -- the minimising prefix may consume a span outside read_length +/- max_gap_bp -- and
// they simply proceed to the recursion, which decides exactly. This narrows nothing: not the band,
// not the gap width, not the model, and it never consults truth.
// Myers/Hyyro bit-parallel, inline and allocation-free. edlib computes exactly this quantity, but at
// 3.4us per call on a 100x103 problem it cannot reach the budget: that cost is per-call setup and
// internal allocation, not the kernel. This form is 709ns, and it was checked to make the IDENTICAL
// admit/reject decision as edlib on 4,000 cases spanning zero through band+3 planted edits,
// indel-bearing reads and unrelated sequence, with zero disagreements. A read longer than the word
// falls back to edlib rather than guessing.
constexpr std::size_t kMyersWidth = 128;
using MyersWord = unsigned __int128;

int two_bit_code(char base) noexcept {
    switch (base) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        default: return 3;   // 'T'; the window is ACGT-validated before this runs
    }
}

// The kernel, once the pattern is in hand. Identical arithmetic whether the masks were prepared
// ahead of time or built on the spot.
bool myers_shw_run(const MyersWord peq[4], std::size_t n, const char* target, std::size_t m,
                   std::size_t band) {
    MyersWord vp = ~MyersWord{0};
    MyersWord vn = 0;
    std::size_t score = n;                       // the empty target prefix costs n deletions
    if (score <= band) return true;
    const MyersWord top = MyersWord{1} << (n - 1);
    for (std::size_t j = 0; j < m; ++j) {
        const MyersWord x0 = peq[two_bit_code(target[j])] | vn;
        const MyersWord d0 = (((x0 & vp) + vp) ^ vp) | x0;
        const MyersWord hn = vp & d0;
        const MyersWord hp = vn | ~(vp | d0);
        if (hp & top) ++score; else if (hn & top) --score;
        const MyersWord x1 = (hp << 1) | 1;
        vn = x1 & d0;
        vp = (hn << 1) | ~(x1 | d0);
        if (score <= band) return true;          // some prefix already fits the budget
    }
    return false;
}

// The same kernel with an exact early refusal. The last-row score is the full query's edit distance to the
// target prefix consumed so far, and adjacent prefixes differ by at most one edit, so after column j no later
// prefix can score below score - (m - 1 - j). Once that bound exceeds the band no prefix can fit it and the
// answer is already false -- the same answer the full scan returns, reached earlier. Acceptances are untouched:
// they return at the first prefix within the band exactly as before.
bool myers_shw_run_early(const MyersWord peq[4], std::size_t n, const char* target, std::size_t m,
                         std::size_t band) {
    MyersWord vp = ~MyersWord{0};
    MyersWord vn = 0;
    std::size_t score = n;
    if (score <= band) return true;
    if (score > band + m) return false;
    const MyersWord top = MyersWord{1} << (n - 1);
    for (std::size_t j = 0; j < m; ++j) {
        const MyersWord x0 = peq[two_bit_code(target[j])] | vn;
        const MyersWord d0 = (((x0 & vp) + vp) ^ vp) | x0;
        const MyersWord hn = vp & d0;
        const MyersWord hp = vn | ~(vp | d0);
        if (hp & top) ++score; else if (hn & top) --score;
        const MyersWord x1 = (hp << 1) | 1;
        vn = x1 & d0;
        vp = (hn << 1) | ~(x1 | d0);
        if (score <= band) return true;
        if (score > band + (m - 1 - j)) return false;   // no remaining prefix can come back within the band
    }
    return false;
}

// MEASUREMENT ONLY: myers_shw_run_early, also reporting the columns it scanned before deciding (0 when it decided
// before the first column) -- a per-request cost proxy for the census.
bool myers_shw_run_early_columns(const MyersWord peq[4], std::size_t n, const char* target, std::size_t m,
                                 std::size_t band, std::uint32_t& columns) {
    MyersWord vp = ~MyersWord{0};
    MyersWord vn = 0;
    std::size_t score = n;
    columns = 0;
    if (score <= band) return true;
    if (score > band + m) return false;
    const MyersWord top = MyersWord{1} << (n - 1);
    for (std::size_t j = 0; j < m; ++j) {
        const MyersWord x0 = peq[two_bit_code(target[j])] | vn;
        const MyersWord d0 = (((x0 & vp) + vp) ^ vp) | x0;
        const MyersWord hn = vp & d0;
        const MyersWord hp = vn | ~(vp | d0);
        if (hp & top) ++score; else if (hn & top) --score;
        const MyersWord x1 = (hp << 1) | 1;
        vn = x1 & d0;
        vp = (hn << 1) | ~(x1 | d0);
        columns = static_cast<std::uint32_t>(j + 1);
        if (score <= band) return true;
        if (score > band + (m - 1 - j)) return false;
    }
    return false;
}

// Many starts of ONE read form (one pattern) through the early-refusing kernel, `Lanes` at a time: each lane
// runs exactly myers_shw_run_early's steps on its own start, lanes advance one column per round so their
// independent dependency chains overlap, and a lane that decides takes the next start at once. Decisions land at
// their own index; nothing a lane does reads another lane. `windows[i] == 0` is decided false, as the
// single-start path decides it.
template <std::size_t Lanes>
void myers_early_lanes(const MyersWord peq[4], std::size_t n, const char* text, const std::uint32_t* starts,
                       const std::uint32_t* windows, std::size_t count, std::size_t band, char* out) {
    MyersWord vp[Lanes], vn[Lanes];
    std::size_t score[Lanes], column[Lanes], width[Lanes], index[Lanes];
    const char* target[Lanes];
    bool live[Lanes];
    const MyersWord top = MyersWord{1} << (n - 1);
    std::size_t next = 0;
    const auto load = [&](std::size_t lane) {
        while (next < count) {
            const std::size_t i = next++;
            const std::size_t m = windows[i];
            if (m == 0) { out[i] = 0; continue; }
            if (n <= band) { out[i] = 1; continue; }
            if (n > band + m) { out[i] = 0; continue; }
            vp[lane] = ~MyersWord{0};
            vn[lane] = 0;
            score[lane] = n;
            column[lane] = 0;
            width[lane] = m;
            index[lane] = i;
            target[lane] = text + starts[i];
            return true;
        }
        return false;
    };
    std::size_t active = 0;
    for (std::size_t lane = 0; lane < Lanes; ++lane) {
        live[lane] = load(lane);
        active += live[lane] ? 1 : 0;
    }
    while (active != 0) {
        for (std::size_t lane = 0; lane < Lanes; ++lane) {
            if (!live[lane]) continue;
            const std::size_t j = column[lane];
            const MyersWord x0 = peq[two_bit_code(target[lane][j])] | vn[lane];
            const MyersWord d0 = (((x0 & vp[lane]) + vp[lane]) ^ vp[lane]) | x0;
            const MyersWord hn = vp[lane] & d0;
            const MyersWord hp = vn[lane] | ~(vp[lane] | d0);
            if (hp & top) ++score[lane]; else if (hn & top) --score[lane];
            const MyersWord x1 = (hp << 1) | 1;
            vn[lane] = x1 & d0;
            vp[lane] = (hn << 1) | ~(x1 | d0);
            column[lane] = j + 1;
            int decided = -1;
            if (score[lane] <= band) decided = 1;
            else if (score[lane] > band + (width[lane] - 1 - j)) decided = 0;
            if (decided < 0) continue;
            out[index[lane]] = static_cast<char>(decided);
            live[lane] = load(lane);
            if (!live[lane]) --active;
        }
    }
}

bool myers_shw_within(const char* query, std::size_t n, const char* target, std::size_t m,
                      std::size_t band) {
    MyersWord peq[4] = {0, 0, 0, 0};
    for (std::size_t i = 0; i < n; ++i) peq[two_bit_code(query[i])] |= MyersWord{1} << i;
    return myers_shw_run(peq, n, target, m, band);
}

MyersWord myers_word_of(std::uint64_t low, std::uint64_t high) {
    return (static_cast<MyersWord>(high) << 64) | static_cast<MyersWord>(low);
}

// The prefilter as it was: every column scanned. Kept as the reference the early-refusing kernel is gated
// against (and timed against in the seeding profile).
bool feasible_within_band_reference(const std::string& haplotype, std::size_t start, std::size_t window,
                                    const std::string& read, std::size_t band, const MyersPattern* pattern) {
    if (window == 0) return false;
    if (pattern != nullptr && pattern->usable) {
        MyersWord peq[4];
        for (int i = 0; i < 4; ++i) peq[i] = myers_word_of(pattern->peq_lo[i], pattern->peq_hi[i]);
        return myers_shw_run(peq, pattern->length, haplotype.data() + start, window, band);
    }
    if (read.size() <= kMyersWidth) {
        return myers_shw_within(read.data(), read.size(), haplotype.data() + start, window, band);
    }
    const EdlibAlignConfig config = edlibNewAlignConfig(
        static_cast<int>(band), EDLIB_MODE_SHW, EDLIB_TASK_DISTANCE, nullptr, 0);
    EdlibAlignResult result =
        edlibAlign(read.data(), static_cast<int>(read.size()), haplotype.data() + start,
                   static_cast<int>(window), config);
    const bool feasible = result.status == EDLIB_STATUS_OK && result.editDistance >= 0;
    edlibFreeAlignResult(result);
    return feasible;
}

bool feasible_within_band(const std::string& haplotype, std::size_t start, std::size_t window,
                          const std::string& read, std::size_t band,
                          const MyersPattern* pattern) {
    if (window == 0) return false;
    // Prepared: the masks were built once for this form during fragment preparation. Rebuilding
    // them here cost a pass over the read on every one of 1.19 billion verifications. The kernel refuses
    // early -- exactly (see myers_shw_run_early): most starts are refusals, and most refusals are decided
    // long before the last column.
    if (pattern != nullptr && pattern->usable) {
        MyersWord peq[4];
        for (int i = 0; i < 4; ++i) peq[i] = myers_word_of(pattern->peq_lo[i], pattern->peq_hi[i]);
        return myers_shw_run_early(peq, pattern->length, haplotype.data() + start, window, band);
    }
    if (read.size() <= kMyersWidth) {
        MyersWord peq[4] = {0, 0, 0, 0};
        for (std::size_t i = 0; i < read.size(); ++i) peq[two_bit_code(read[i])] |= MyersWord{1} << i;
        return myers_shw_run_early(peq, read.size(), haplotype.data() + start, window, band);
    }
    const EdlibAlignConfig config = edlibNewAlignConfig(
        static_cast<int>(band), EDLIB_MODE_SHW, EDLIB_TASK_DISTANCE, nullptr, 0);
    EdlibAlignResult result =
        edlibAlign(read.data(), static_cast<int>(read.size()), haplotype.data() + start,
                   static_cast<int>(window), config);
    const bool feasible = result.status == EDLIB_STATUS_OK && result.editDistance >= 0;
    edlibFreeAlignResult(result);
    return feasible;
}

// The same recursion in SCALED PROBABILITY SPACE.
//
// The log-space form evaluates log_add_exp -- an exp and a log1p behind a call that does not inline
// across translation units -- about 39,000 times for a 100 bp read at band 6, which measured 77us
// per placement and is what put a C4 donor into hours. Here a transition is a multiply and a merge
// is an add. Dynamic range is held by renormalising each completed read row by its maximum and
// carrying log(scale) alongside, which is the standard HMM scaling and is exact up to rounding.
//
// The row is scaled only AFTER the deletion closure: D resolves within the row and reads cells that
// scaling would otherwise have rewritten underneath it.
double gap_recursion_scaled(
    const std::string* haplotype,
    const std::string* read,
    std::size_t start,
    std::size_t read_length,
    std::size_t available,
    std::size_t band,
    const LibraryModel& library,
    bool emitting,
    ComponentList* components) {
    const std::size_t delta = library.max_gap_bp;
    const std::size_t width = 2 * delta + 1;
    const std::size_t budget = band + 1;
    const double p_match = 1.0 - library.substitution_error;
    const double p_mismatch = library.substitution_error / 3.0;
    // Summed over all 4^L reads a mismatch position carries 3 * (e/3) = e.
    const double p_any_mismatch = library.substitution_error;
    const double p_insert_emission = 0.25;
    const double t_continue = library.gap_open > 0.0 ? 1.0 - 2.0 * library.gap_open : 1.0;
    const double t_open = library.gap_open;
    const double t_extend = library.gap_extend;
    const double t_close = 1.0 - library.gap_extend;

    const std::size_t cells = width * kStates * budget;
    GapScratch& scratch = gap_scratch();
    scratch.ensure(cells);
    double* previous = scratch.first.data();
    double* current = scratch.second.data();
    std::fill(previous, previous + cells, 0.0);

    const auto at = [&](std::size_t displacement, std::size_t state, std::size_t edits) {
        return (displacement * kStates + state) * budget + edits;
    };
    previous[at(delta, kMatch, 0)] = 1.0;
    double log_scale = 0.0;

    // Row base pointers, hoisted. The index arithmetic (d*kStates+s)*budget+e was recomputed on
    // every one of the ~15,000 cell accesses per read, and the innermost loop strided across
    // non-contiguous cells. Walking `e` contiguously off a hoisted base is the same arithmetic in a
    // cache-friendly order -- the measured 2.3ns per cell update was memory traffic and loop
    // overhead, not the multiply-adds.
    const auto row = [&](double* plane, std::size_t displacement, std::size_t state) {
        return plane + (displacement * kStates + state) * budget;
    };

    for (std::size_t i = 1; i <= read_length; ++i) {
        std::fill(current, current + cells, 0.0);
        // The largest value written this row, tracked as we go instead of rescanning all 147 cells.
        double best = 0.0;
        for (std::size_t d = 0; d < width; ++d) {
            const long k = static_cast<long>(d) - static_cast<long>(delta);
            // Reaching displacement k costs at least |k| gap bases, so every lower edit level at
            // this displacement is unreachable. Exact bookkeeping, not a pruning heuristic.
            const std::size_t lowest = static_cast<std::size_t>(k < 0 ? -k : k);
            if (lowest >= budget) continue;
            const long consumed_before = static_cast<long>(i) - 1 + k;
            if (consumed_before >= 0 && static_cast<std::size_t>(consumed_before) < available) {
                const std::size_t index = start + static_cast<std::size_t>(consumed_before);
                double weight[2] = {0.0, 0.0};
                std::size_t cost[2] = {0, 1};
                std::size_t branches = 0;
                if (emitting) {
                    const bool same = (*haplotype)[index] == (*read)[i - 1];
                    weight[0] = same ? p_match : p_mismatch;
                    cost[0] = same ? 0 : 1;
                    branches = 1;
                } else {
                    weight[0] = p_match;         cost[0] = 0;
                    weight[1] = p_any_mismatch;  cost[1] = 1;
                    branches = 2;
                }
                const double* from_match = row(previous, d, kMatch);
                const double* from_insert = row(previous, d, kInsert);
                const double* from_delete = row(previous, d, kDelete);
                double* into = row(current, d, kMatch);
                for (std::size_t b = 0; b < branches; ++b) {
                    const double w = weight[b];
                    if (w == 0.0) continue;
                    const std::size_t c = cost[b];
                    for (std::size_t e = std::max(lowest, c); e < budget; ++e) {
                        const std::size_t f = e - c;
                        const double arrived = from_match[f] * t_continue +
                                               (from_insert[f] + from_delete[f]) * t_close;
                        if (arrived == 0.0) continue;
                        const double value = into[e] + arrived * w;
                        into[e] = value;
                        if (value > best) best = value;
                    }
                }
            }
            // I: consume one read base only; displacement falls by one.
            if (d + 1 < width) {
                const double emission = emitting ? p_insert_emission : 1.0;
                const double* from_match = row(previous, d + 1, kMatch);
                const double* from_insert = row(previous, d + 1, kInsert);
                double* into = row(current, d, kInsert);
                for (std::size_t e = std::max<std::size_t>(lowest, 1); e < budget; ++e) {
                    const std::size_t f = e - 1;
                    const double arrived = from_match[f] * t_open + from_insert[f] * t_extend;
                    if (arrived == 0.0) continue;
                    const double value = into[e] + arrived * emission;
                    into[e] = value;
                    if (value > best) best = value;
                }
            }
        }
        // D: consume haplotype only, resolved within this read position, sweeping displacement
        // upward because D raises k by one.
        for (std::size_t d = 1; d < width; ++d) {
            const long k = static_cast<long>(d) - static_cast<long>(delta);
            const std::size_t lowest = static_cast<std::size_t>(k < 0 ? -k : k);
            if (lowest >= budget) continue;
            const long consumed_before = static_cast<long>(i) + k - 1;
            if (consumed_before < 0) continue;
            if (static_cast<std::size_t>(consumed_before) >= available) continue;
            const double* from_match = row(current, d - 1, kMatch);
            const double* from_delete = row(current, d - 1, kDelete);
            double* into = row(current, d, kDelete);
            for (std::size_t e = std::max<std::size_t>(lowest, 1); e < budget; ++e) {
                const std::size_t f = e - 1;
                const double arrived = from_match[f] * t_open + from_delete[f] * t_extend;
                if (arrived == 0.0) continue;
                const double value = into[e] + arrived;
                into[e] = value;
                if (value > best) best = value;
            }
        }
        if (best == 0.0) return kNegInf;   // nothing admissible survives; the mass is exactly zero
        // Rescale only when the row drifts far from unity. Doubles carry ~300 decades of headroom,
        // so a scale every few dozen rows is ample, and this removes ~99 of the 100 log() calls and
        // the 147-cell multiply that went with each. Correctness does not depend on WHEN we rescale,
        // only that the carried log_scale stays consistent with the row.
        if (best < 1e-100 || best > 1e100) {
            const double inverse = 1.0 / best;
            for (std::size_t c = 0; c < cells; ++c) current[c] *= inverse;
            log_scale += std::log(best);
        }
        std::swap(previous, current);
    }

    double total = 0.0;
    const std::size_t floor_bp = read_length - std::min(read_length, delta);
    for (std::size_t d = 0; d < width; ++d) {
        const long k = static_cast<long>(d) - static_cast<long>(delta);
        const long consumed = static_cast<long>(read_length) + k;
        if (consumed < 0 || static_cast<std::size_t>(consumed) > available) continue;
        if (static_cast<std::size_t>(consumed) < floor_bp) continue;
        double span = 0.0;
        std::size_t span_edits = budget;
        for (std::size_t s = 0; s < kStates; ++s) {
            for (std::size_t e = 0; e < budget; ++e) {
                const double value = previous[at(d, s, e)];
                if (value == 0.0) continue;
                span += value;
                if (e < span_edits) span_edits = e;
            }
        }
        if (span == 0.0) continue;
        total += span;
        if (components != nullptr) {
            components->push_back(
                {static_cast<std::size_t>(consumed), std::log(span) + log_scale, span_edits});
        }
    }
    if (total == 0.0) return kNegInf;
    return std::log(total) + log_scale;
}

template <typename F, std::size_t... I>
inline void for_each_index_impl(F&& f, std::index_sequence<I...>) {
    (f(std::integral_constant<std::size_t, I>{}), ...);
}

// Calls f(std::integral_constant<std::size_t, i>) for i = 0 .. N-1, in order.
template <std::size_t N, typename F>
inline void for_each_index(F&& f) {
    for_each_index_impl(std::forward<F>(f), std::make_index_sequence<N>{});
}

// I2-H (EXPERIMENTAL, default off): the EMITTING recursion with its geometry fixed at compile time.
//
// The same arithmetic as gap_recursion_scaled(emitting = true), cell for cell; only the bookkeeping differs. In
// the emitting pass every cell of a row has exactly one writer -- M and I from the previous row, D from this row's
// M and D one displacement lower -- so each cell is ASSIGNED once (its value, or the zero the reference's fill
// would have left) instead of zero-filled and then accumulated into. The accumulations removed were all into exact
// zeros: 0 + x*w rounds once to x*w whether or not it is fused, and 0 + arrived is arrived; every value is
// non-negative, so no signed zero arises. The data-dependent skips (arrived == 0, w == 0) left exactly the zeros
// now written, and the row maximum is a maximum over the same non-negative values in another order. Every
// `arrived` expression is the reference's source text, so floating-point contraction is the same. Rescaling, the
// terminal reduction and the components are the reference's code on the same layout. Two planes of at most
// 7 x 3 x 9 doubles live on the stack for one call; no thread-local scratch is touched.
template <std::size_t Delta, std::size_t Budget>
double gap_recursion_fixed(const std::string& haplotype, const std::string& read, std::size_t start,
                           std::size_t read_length, std::size_t available, const LibraryModel& library,
                           ComponentList* components) {
    constexpr std::size_t width = 2 * Delta + 1;
    constexpr std::size_t cells = width * kStates * Budget;
    const std::size_t delta = Delta;
    const std::size_t budget = Budget;
    const double p_match = 1.0 - library.substitution_error;
    const double p_mismatch = library.substitution_error / 3.0;
    const double p_insert_emission = 0.25;
    const double t_continue = library.gap_open > 0.0 ? 1.0 - 2.0 * library.gap_open : 1.0;
    const double t_open = library.gap_open;
    const double t_extend = library.gap_extend;
    const double t_close = 1.0 - library.gap_extend;

    double planes[2][cells];
    double* previous = planes[0];
    double* current = planes[1];
    std::fill(previous, previous + cells, 0.0);
    const auto at = [](std::size_t displacement, std::size_t state, std::size_t edits) {
        return (displacement * kStates + state) * Budget + edits;
    };
    previous[at(Delta, kMatch, 0)] = 1.0;
    double log_scale = 0.0;

    for (std::size_t i = 1; i <= read_length; ++i) {
        double best = 0.0;
        // This row's M and D one displacement lower, carried for the deletion step.
        double match_below[Budget];
        double delete_below[Budget];
        for_each_index<width>([&](auto displacement) {
            constexpr std::size_t d = decltype(displacement)::value;
            constexpr long k = static_cast<long>(d) - static_cast<long>(Delta);
            constexpr std::size_t lowest = static_cast<std::size_t>(k < 0 ? -k : k);
            constexpr std::size_t lowest_gap = lowest > 1 ? lowest : 1;
            double m[Budget];
            double ins[Budget];
            double del[Budget];
            for (std::size_t e = 0; e < Budget; ++e) m[e] = ins[e] = del[e] = 0.0;
            if constexpr (lowest < Budget) {
                const long consumed_before = static_cast<long>(i) - 1 + k;
                const bool in_haplotype =
                    consumed_before >= 0 && static_cast<std::size_t>(consumed_before) < available;
                // M: consume one read base and one haplotype base.
                if (in_haplotype) {
                    const std::size_t index = start + static_cast<std::size_t>(consumed_before);
                    const double* from_match = previous + (d * kStates + kMatch) * Budget;
                    const double* from_insert = previous + (d * kStates + kInsert) * Budget;
                    const double* from_delete = previous + (d * kStates + kDelete) * Budget;
                    if (haplotype[index] == read[i - 1]) {
                        const double w = p_match;
                        for (std::size_t e = lowest; e < Budget; ++e) {
                            const std::size_t f = e;
                            const double arrived = from_match[f] * t_continue +
                                                   (from_insert[f] + from_delete[f]) * t_close;
                            m[e] = arrived * w;
                        }
                    } else {
                        const double w = p_mismatch;
                        for (std::size_t e = lowest_gap; e < Budget; ++e) {
                            const std::size_t f = e - 1;
                            const double arrived = from_match[f] * t_continue +
                                                   (from_insert[f] + from_delete[f]) * t_close;
                            m[e] = arrived * w;
                        }
                    }
                }
                // I: consume one read base only; displacement falls by one.
                if constexpr (d + 1 < width) {
                    const double* from_match = previous + ((d + 1) * kStates + kMatch) * Budget;
                    const double* from_insert = previous + ((d + 1) * kStates + kInsert) * Budget;
                    for (std::size_t e = lowest_gap; e < Budget; ++e) {
                        const std::size_t f = e - 1;
                        const double arrived = from_match[f] * t_open + from_insert[f] * t_extend;
                        ins[e] = arrived * p_insert_emission;
                    }
                }
                // D: consume haplotype only, from this row's M and D one displacement lower.
                if constexpr (d >= 1) {
                    if (in_haplotype) {
                        for (std::size_t e = lowest_gap; e < Budget; ++e) {
                            const std::size_t f = e - 1;
                            const double arrived = match_below[f] * t_open + delete_below[f] * t_extend;
                            del[e] = arrived;
                        }
                    }
                }
            }
            double* into_match = current + (d * kStates + kMatch) * Budget;
            double* into_insert = current + (d * kStates + kInsert) * Budget;
            double* into_delete = current + (d * kStates + kDelete) * Budget;
            for (std::size_t e = 0; e < Budget; ++e) {
                into_match[e] = m[e];
                into_insert[e] = ins[e];
                into_delete[e] = del[e];
                best = std::max(best, m[e]);
                best = std::max(best, ins[e]);
                best = std::max(best, del[e]);
                match_below[e] = m[e];
                delete_below[e] = del[e];
            }
        });
        if (best == 0.0) return kNegInf;   // nothing admissible survives; the mass is exactly zero
        if (best < 1e-100 || best > 1e100) {
            const double inverse = 1.0 / best;
            for (std::size_t c = 0; c < cells; ++c) current[c] *= inverse;
            log_scale += std::log(best);
        }
        std::swap(previous, current);
    }

    double total = 0.0;
    const std::size_t floor_bp = read_length - std::min(read_length, delta);
    for (std::size_t d = 0; d < width; ++d) {
        const long k = static_cast<long>(d) - static_cast<long>(delta);
        const long consumed = static_cast<long>(read_length) + k;
        if (consumed < 0 || static_cast<std::size_t>(consumed) > available) continue;
        if (static_cast<std::size_t>(consumed) < floor_bp) continue;
        double span = 0.0;
        std::size_t span_edits = budget;
        for (std::size_t s = 0; s < kStates; ++s) {
            for (std::size_t e = 0; e < budget; ++e) {
                const double value = previous[at(d, s, e)];
                if (value == 0.0) continue;
                span += value;
                if (e < span_edits) span_edits = e;
            }
        }
        if (span == 0.0) continue;
        total += span;
        if (components != nullptr) {
            components->push_back(
                {static_cast<std::size_t>(consumed), std::log(span) + log_scale, span_edits});
        }
    }
    if (total == 0.0) return kNegInf;
    return std::log(total) + log_scale;
}

using FixedRecursion = double (*)(const std::string&, const std::string&, std::size_t, std::size_t, std::size_t,
                                  const LibraryModel&, ComponentList*);

template <std::size_t Delta, std::size_t... Index>
constexpr std::array<FixedRecursion, sizeof...(Index)> fixed_recursions(std::index_sequence<Index...>) {
    return {&gap_recursion_fixed<Delta, Index + 2>...};   // mismatch band Index + 1, so Budget Index + 2
}

// The supported geometries: max_gap_bp 1..3, mismatch band 1..8. Anything else keeps the scaled recursion.
FixedRecursion fixed_recursion_for(std::size_t max_gap_bp, std::size_t band) {
    static constexpr std::array<std::array<FixedRecursion, 8>, 3> table = {
        fixed_recursions<1>(std::make_index_sequence<8>{}), fixed_recursions<2>(std::make_index_sequence<8>{}),
        fixed_recursions<3>(std::make_index_sequence<8>{})};
    if (max_gap_bp < 1 || max_gap_bp > 3 || band < 1 || band > 8) return nullptr;
    return table[max_gap_bp - 1][band - 1];
}

// Fixed-geometry recursions executed on this thread (gates only: proves the experimental path actually ran).
thread_local std::uint64_t t_fixed_recursions = 0;

// The emitting recursion, by backend. The experimental backend falls back to the reference for any geometry it
// does not support; the normaliser never comes here.
double gap_recursion_emitting(const std::string& haplotype, const std::string& read, std::size_t start,
                              std::size_t read_length, std::size_t available, std::size_t band,
                              const LibraryModel& library, ComponentList* components, GapRecursionBackend backend) {
    if (backend == GapRecursionBackend::fixed_geometry_experimental) {
        if (const FixedRecursion fixed = fixed_recursion_for(library.max_gap_bp, band)) {
            ++t_fixed_recursions;
            return fixed(haplotype, read, start, read_length, available, library, components);
        }
    }
    return gap_recursion_scaled(&haplotype, &read, start, read_length, available, band, library, true, components);
}

} // namespace

void set_emission_cache_hash_collapse(bool collapse) {
    cache_hash_collapse().store(collapse, std::memory_order_relaxed);
}

namespace scoring_detail {
// A persistent context arms the cache for its whole life. While it does, begin_emission_cache must
// refuse rather than clear the table out from under it: the context's keys point into ITS buffers,
// and a one-shot pass that wiped them would silently recompute everything while reporting hits.
std::atomic<bool>& emission_cache_pinned() {
    static std::atomic<bool> pinned{false};
    return pinned;
}
}  // namespace scoring_detail

void begin_emission_cache(std::size_t budget_bytes) {
    if (scoring_detail::emission_cache_pinned().load(std::memory_order_relaxed)) {
        throw std::runtime_error(
            "the emission cache is held by a live PreparedScoringContext; destroy it before starting "
            "a one-shot scoring pass");
    }
    EmissionCache& cache = emission_cache();
    cache.active.store(false, std::memory_order_relaxed);
    // Discard whatever this thread accumulated BEFORE the cache was armed, so that "what the cache
    // did" means work done during its lifetime rather than everything the thread ever did. Without
    // this the end-of-life fold reports the caller's earlier uncached work as the cache's own:
    // measured at executed_dps 54 against 18 genuine misses, the difference being exactly a 36-read
    // pre-pass. Same shape as the residue defect fixed in placement -- a shared accumulator whose
    // contents depend on what ran first.
    (void)take_gapped_emission_counters();
    for (CacheShard& shard : cache.shards) {
        const std::lock_guard<std::mutex> guard(shard.guard);
        shard.table.clear();
        shard.occupied = 0;
        shard.slabs.clear();
        shard.slab_used = 0;
        shard.arena.clear();
    }
    cache.lookup_nanos.store(0, std::memory_order_relaxed);
    cache.budget_bytes.store(budget_bytes, std::memory_order_relaxed);
    cache.budgeted_bytes_used.store(0, std::memory_order_relaxed);
    cache.entries.store(0, std::memory_order_relaxed);
    cache.payload_bytes.store(0, std::memory_order_relaxed);
    cache.entries_not_retained.store(0, std::memory_order_relaxed);
    // Zero is off, and off reproduces the uncached path exactly rather than approximately.
    cache.active.store(budget_bytes > 0, std::memory_order_relaxed);
}

EmissionCacheCounters read_emission_cache_counters() {
    EmissionCache& cache = emission_cache();
    // Only the atomic-backed half. The observational counters are worker-local and are merged by
    // the caller, which is why this must not zero them.
    EmissionCacheCounters out;
    out.payload_bytes = cache.payload_bytes.load(std::memory_order_relaxed);
    out.entries = cache.entries.load(std::memory_order_relaxed);
    out.entries_not_retained = cache.entries_not_retained.load(std::memory_order_relaxed);
    // What the declared budget is actually holding, separate from the fits it does not govern.
    out.budgeted_bytes_used = cache.budgeted_bytes_used.load(std::memory_order_relaxed);
    for (CacheShard& shard : cache.shards) {
        const std::lock_guard<std::mutex> guard(shard.guard);
        out.entry_metadata_bytes +=
            static_cast<std::uint64_t>(shard.table.capacity()) * sizeof(ArenaEntry);
        out.content_arena_bytes +=
            static_cast<std::uint64_t>(shard.slabs.size()) * kSlabBytes;
    }
    return out;
}

// Table and arena bytes actually held, summed across shards. Read before teardown, under each
// shard's own lock. Counted rather than estimated: kCacheEntryBytes is the budget's declared unit,
// which is not the same claim as what the container really holds.
void emission_cache_storage_bytes(std::uint64_t* table_bytes, std::uint64_t* arena_bytes) {
    EmissionCache& cache = emission_cache();
    std::uint64_t table = 0;
    std::uint64_t arena = 0;
    for (CacheShard& shard : cache.shards) {
        const std::lock_guard<std::mutex> guard(shard.guard);
        // Per entry now, not a count times a constant: the key owns a variable number of bytes, so
        // the only honest figure is the sum of what each one actually holds.
        table += static_cast<std::uint64_t>(shard.table.capacity()) * sizeof(ArenaEntry);
        arena += static_cast<std::uint64_t>(shard.arena.capacity()) * sizeof(GappedComponent) +
                 static_cast<std::uint64_t>(shard.slabs.size()) * kSlabBytes;
    }
    if (table_bytes != nullptr) *table_bytes = table;
    if (arena_bytes != nullptr) *arena_bytes = arena;
}

EmissionCacheCounters& last_emission_cache_slot() {
    static EmissionCacheCounters last;
    return last;
}

EmissionCacheCounters last_emission_cache_counters() { return last_emission_cache_slot(); }

ScorerMemoryAccounts& last_memory_slot() {
    static ScorerMemoryAccounts last;
    return last;
}

ScorerMemoryAccounts last_scorer_memory_accounts() { return last_memory_slot(); }

void set_emission_lookup_timing(bool enabled) { emission_lookup_timing().store(enabled, std::memory_order_relaxed); }
bool emission_lookup_timing_enabled() { return emission_lookup_timing().load(std::memory_order_relaxed); }

double emission_cache_lookup_seconds() {
    return static_cast<double>(emission_cache().lookup_nanos.load(std::memory_order_relaxed)) / 1e9;
}

void reset_emission_cache_lookup_seconds() {
    emission_cache().lookup_nanos.store(0, std::memory_order_relaxed);
}

struct EmissionCacheHandoff::Contents {
    struct Shard {
        std::vector<ArenaEntry> table;
        std::size_t occupied = 0;
        std::vector<std::unique_ptr<char[]>> slabs;
        std::size_t slab_used = 0;
        std::vector<GappedComponent> arena;
    };
    std::array<Shard, kCacheShards> shards;
    std::uint64_t entries = 0, budgeted_bytes = 0, allocated_bytes = 0;
    // every parameter a cached fit depends on besides its key's bytes (read and window, whose lengths carry the
    // normaliser's read length and clamped available length)
    double substitution_error = 0.0, gap_open = 0.0, gap_extend = 0.0;
    std::size_t max_gap_bp = 0, band = 0;
    GapRecursionBackend backend = GapRecursionBackend::scaled_reference;
};

EmissionCacheHandoff::EmissionCacheHandoff() = default;
EmissionCacheHandoff::~EmissionCacheHandoff() = default;
bool EmissionCacheHandoff::holds() const { return contents != nullptr; }
std::uint64_t EmissionCacheHandoff::entries() const { return contents ? contents->entries : 0; }
std::uint64_t EmissionCacheHandoff::allocated_bytes() const { return contents ? contents->allocated_bytes : 0; }
std::uint64_t EmissionCacheHandoff::budgeted_bytes() const { return contents ? contents->budgeted_bytes : 0; }
void EmissionCacheHandoff::clear() { contents.reset(); }

namespace {
std::mutex& handoff_report_mutex() {
    static std::mutex m;
    return m;
}
EmissionHandoffReport& handoff_report_slot() {
    static EmissionHandoffReport r;
    return r;
}

bool same_double_bits(double a, double b) { return std::memcmp(&a, &b, sizeof(double)) == 0; }

// The emission cache's teardown. With `keep`, the shards' storage is moved into it instead of being discarded.
EmissionCacheCounters finish_emission_cache(EmissionCacheHandoff::Contents* keep);
}  // namespace

EmissionHandoffReport emission_handoff_report() {
    const std::lock_guard<std::mutex> lock(handoff_report_mutex());
    EmissionHandoffReport r = handoff_report_slot();
    if (r.adopted) r.inherited_hits = emission_cache().inherited_hits.load(std::memory_order_relaxed);
    return r;
}

void reset_emission_handoff_report() {
    const std::lock_guard<std::mutex> lock(handoff_report_mutex());
    handoff_report_slot() = EmissionHandoffReport{};
    emission_cache().inherited_hits.store(0, std::memory_order_relaxed);
}

// I2-K: ends the cache as end_emission_cache does, moving its storage into `out` (replacing anything it held).
void detach_emission_cache(EmissionCacheHandoff& out, const LibraryModel& library, std::size_t band,
                           GapRecursionBackend backend, bool record = true) {
    auto contents = std::make_unique<EmissionCacheHandoff::Contents>();
    contents->budgeted_bytes = emission_cache().budgeted_bytes_used.load(std::memory_order_relaxed);
    (void)finish_emission_cache(contents.get());
    for (const auto& shard : contents->shards) {
        contents->entries += shard.occupied;
        contents->allocated_bytes += static_cast<std::uint64_t>(shard.table.capacity()) * sizeof(ArenaEntry) +
                                     static_cast<std::uint64_t>(shard.slabs.size()) * kSlabBytes +
                                     static_cast<std::uint64_t>(shard.arena.capacity()) * sizeof(GappedComponent);
    }
    contents->substitution_error = library.substitution_error;
    contents->gap_open = library.gap_open;
    contents->gap_extend = library.gap_extend;
    contents->max_gap_bp = library.max_gap_bp;
    contents->band = band;
    contents->backend = backend;
    if (record) {
        const std::lock_guard<std::mutex> lock(handoff_report_mutex());
        EmissionHandoffReport& r = handoff_report_slot();
        r.detached = true;
        r.entries = contents->entries;
        r.allocated_bytes = contents->allocated_bytes;
        r.budgeted_bytes = contents->budgeted_bytes;
    }
    out.contents = std::move(contents);
}

// I2-K: called right after begin_emission_cache by the adopting context. Consumes the handoff whatever happens.
std::string adopt_emission_cache(EmissionCacheHandoff* handoff, const LibraryModel& library, std::size_t band,
                                 GapRecursionBackend backend, std::size_t budget_bytes, bool record = true) {
    if (handoff == nullptr) return "not_offered";
    std::unique_ptr<EmissionCacheHandoff::Contents> c = std::move(handoff->contents);
    const std::lock_guard<std::mutex> lock(handoff_report_mutex());
    EmissionHandoffReport scratch;   // S2-D: a batch-owned adoption leaves the process-wide I2-K record untouched
    EmissionHandoffReport& r = record ? handoff_report_slot() : scratch;
    if (!c) {   // nothing to adopt; a context created after an adoption leaves that record as it is
        if (!r.offered) r.refusal = "not_offered";
        return "not_offered";
    }
    r.offered = true;
    if (!same_double_bits(c->substitution_error, library.substitution_error) || !same_double_bits(c->gap_open, library.gap_open) ||
        !same_double_bits(c->gap_extend, library.gap_extend) || c->max_gap_bp != library.max_gap_bp || c->band != band ||
        c->backend != backend) {
        r.refusal = "context";
        return r.refusal;   // the contents are freed here; the adaptive search computes as without a handoff
    }
    EmissionCache& cache = emission_cache();
    if (budget_bytes == 0 || !cache.active.load(std::memory_order_relaxed)) { r.refusal = "cache_disabled"; return r.refusal; }
    if (c->budgeted_bytes > budget_bytes) { r.refusal = "budget"; return r.refusal; }
    for (std::size_t i = 0; i < kCacheShards; ++i) {
        CacheShard& shard = cache.shards[i];
        auto& from = c->shards[i];
        const std::lock_guard<std::mutex> guard(shard.guard);
        shard.table.swap(from.table);
        shard.occupied = from.occupied;
        shard.slabs.swap(from.slabs);
        shard.slab_used = from.slab_used;
        shard.arena.swap(from.arena);
        for (ArenaEntry& entry : shard.table) {
            if (entry.state != 0) entry.state = 2;   // inherited: hits on it are counted
        }
    }
    // The budget governs every retained entry, inherited ones included; the stage's own counters start at zero.
    cache.budgeted_bytes_used.store(c->budgeted_bytes, std::memory_order_relaxed);
    cache.inherited_hits.store(0, std::memory_order_relaxed);
    r.adopted = true;
    r.refusal = "none";
    return r.refusal;
}

EmissionCacheCounters end_emission_cache() { return finish_emission_cache(nullptr); }

RetainedProcessStorage retained_process_storage() {
    RetainedProcessStorage r;
    EmissionCache& cache = emission_cache();
    r.emission_cache_armed = cache.active.load(std::memory_order_relaxed) ||
                             scoring_detail::emission_cache_pinned().load(std::memory_order_relaxed);
    for (CacheShard& shard : cache.shards) {
        const std::lock_guard<std::mutex> guard(shard.guard);
        r.emission_cache_table_bytes += static_cast<std::uint64_t>(shard.table.capacity()) * sizeof(ArenaEntry);
        r.emission_cache_component_bytes += static_cast<std::uint64_t>(shard.arena.capacity()) * sizeof(GappedComponent);
        r.emission_cache_slab_bytes += static_cast<std::uint64_t>(shard.slabs.size()) * kSlabBytes;
    }
    {
        const std::lock_guard<std::mutex> guard(census_mutex());
        auto& table = census_map();
        r.emission_census_bytes = static_cast<std::uint64_t>(table.bucket_count()) * sizeof(void*) +
                                  static_cast<std::uint64_t>(table.size()) * (sizeof(CensusKey) + sizeof(CensusEntry) + 2 * sizeof(void*));
    }
    const GapScratch& scratch = gap_scratch();
    r.calling_thread_scratch_bytes = (scratch.first.capacity() + scratch.second.capacity()) * sizeof(double);
    return r;
}

RetainedProcessStorage release_sample_retained_storage() {
    const RetainedProcessStorage before = retained_process_storage();
    if (before.emission_cache_armed) {
        throw std::logic_error("release_sample_retained_storage: the emission cache is armed or pinned; release only "
                               "between finished samples");
    }
    for (CacheShard& shard : emission_cache().shards) {
        const std::lock_guard<std::mutex> guard(shard.guard);
        std::vector<ArenaEntry>().swap(shard.table);
        shard.occupied = 0;
        std::vector<std::unique_ptr<char[]>>().swap(shard.slabs);
        shard.slab_used = 0;
        std::vector<GappedComponent>().swap(shard.arena);
    }
    {
        const std::lock_guard<std::mutex> guard(census_mutex());
        std::unordered_map<CensusKey, CensusEntry, CensusKeyHash>().swap(census_map());
    }
    GapScratch& scratch = gap_scratch();
    std::vector<double>().swap(scratch.first);
    std::vector<double>().swap(scratch.second);
    return before;
}

namespace {
EmissionCacheCounters finish_emission_cache(EmissionCacheHandoff::Contents* keep) {
    EmissionCache& cache = emission_cache();
    cache.active.store(false, std::memory_order_relaxed);
    // One aggregation path, shared with the non-destructive read, so the two cannot drift.
    EmissionCacheCounters out = read_emission_cache_counters();
    // Fold in whatever THIS thread accumulated and nobody drained. The observational counters moved
    // into a thread-local, and only collect_form_placements drains it -- so a caller that uses
    // gapped_emission directly (the probe's diagnostics, the oracle, any future one) would be told
    // the cache did nothing at all while it was plainly filling. Reporting zero work beside a
    // non-zero entry count is worse than not reporting: it reads as a broken cache rather than a
    // broken instrument.
    const GappedEmissionCounters pending = take_gapped_emission_counters();
    out.executed_prefilters += pending.executed_prefilters;
    out.dp_cache_requests += pending.dp_cache_requests;
    out.dp_cache_hits += pending.dp_cache_hits;
    out.dp_cache_misses += pending.dp_cache_misses;
    out.executed_dps += pending.executed_dps;
    // Deliberately does NOT write the retained slot. Teardown runs after score_diplotypes has
    // already stored the MERGED worker totals there, and this view has none of them: every worker's
    // thread-local was drained by collect_form_placements, so folding here yields zeros that would
    // stamp over the real figures. Only the atomic-backed fields survived that overwrite, which is
    // exactly the symptom seen -- payload_bytes and entries correct beside five zeroed counters.
    // Direct callers use this function's RETURN value, which is complete for their single thread.
    // Keys own their bytes now, so this is reclaiming memory rather than avoiding dangling reads.
    for (std::size_t i = 0; i < kCacheShards; ++i) {
        CacheShard& shard = cache.shards[i];
        const std::lock_guard<std::mutex> guard(shard.guard);
        if (keep != nullptr) {   // I2-K: the storage moves to the handoff; the shard is left empty
            auto& to = keep->shards[i];
            to.table.swap(shard.table);
            to.occupied = shard.occupied;
            to.slabs.swap(shard.slabs);
            to.slab_used = shard.slab_used;
            to.arena.swap(shard.arena);
        }
        shard.table.clear();
        shard.occupied = 0;
        shard.slabs.clear();
        shard.slab_used = 0;
        shard.arena.clear();
    }
    return out;
}
}  // namespace

void begin_emission_census() {
    const std::lock_guard<std::mutex> guard(census_mutex());
    census_map().clear();
    census_key_bytes() = 0;
    census_active().store(true, std::memory_order_relaxed);
}

EmissionCensus end_emission_census() {
    census_active().store(false, std::memory_order_relaxed);
    const std::lock_guard<std::mutex> guard(census_mutex());
    EmissionCensus out;
    out.distinct = census_map().size();
    out.key_bytes = census_key_bytes();
    for (const auto& item : census_map()) {
        const CensusEntry& entry = item.second;
        out.requests += entry.requests;
        out.value_bytes += entry.value_bytes;
        if (entry.executed_dp) {
            out.dp_requests += entry.requests;
            ++out.dp_distinct;
        } else {
            out.prefilter_refusal_requests += entry.requests;
            ++out.prefilter_refusal_distinct;
        }
    }
    census_map().clear();
    census_key_bytes() = 0;
    return out;
}

GappedEmissionCounters take_gapped_emission_counters() {
    EmissionWork& work = emission_work();
    GappedEmissionCounters out;
    out.dp_evaluations = work.dp_evaluations;
    out.prefilter_rejections = work.prefilter_rejections;
    out.executed_prefilters = work.executed_prefilters;
    out.dp_cache_requests = work.dp_cache_requests;
    out.dp_cache_hits = work.dp_cache_hits;
    out.dp_cache_misses = work.dp_cache_misses;
    out.executed_dps = work.executed_dps;
    // Folded in here rather than per lookup: once per haplotype instead of 8.2M times a round.
    emission_cache().lookup_nanos.fetch_add(work.cache_lookup_nanos, std::memory_order_relaxed);
    work = EmissionWork{};
    return out;
}

double log_gap_shape_mass(
    std::size_t read_length,
    std::size_t available_bp,
    std::size_t mismatch_band,
    const LibraryModel& library) {
    if (read_length == 0) throw std::invalid_argument("gapped emission requires a non-empty read");
    // ONE normaliser, summed over every terminal span. Deliberately not per span: the components are
    // parts of a single distribution over reads of this length, not several distributions, and
    // dividing each by its own shape mass would make each of them sum to one on its own and their
    // total several times that.
    //
    // MEMOISED, because this is a pure function of its arguments and gapped_emission calls it on
    // every verification -- which made each verification run the whole recursion TWICE, once for the
    // numerator and once for a denominator that depends on neither sequence. Measured at 175us of a
    // 255us call: 69% of the cost of every placement in the program.
    //
    // The key space is tiny once `available_bp` is clamped. The haplotype end cuts admissible shapes
    // only while it is within reach; once available_bp reaches read_length + max_gap_bp, no terminal
    // span and no intermediate path is excluded, and every larger value yields a BIT-IDENTICAL
    // result -- measured across available_bp from 103 to 200,000 at read_length 100, max_gap_bp 3.
    // So the clamp leaves exactly max_gap_bp + 1 distinct entries per read length, and it is an
    // identity rather than a tolerance.
    const std::size_t key_available = std::min(available_bp, read_length + library.max_gap_bp);

    struct ShapeKey {
        std::size_t read_length;
        std::size_t available;
        std::size_t band;
        std::size_t max_gap;
        double substitution_error;
        double gap_open;
        double gap_extend;

        bool operator<(const ShapeKey& other) const noexcept {
            return std::tie(read_length, available, band, max_gap, substitution_error, gap_open,
                            gap_extend) <
                   std::tie(other.read_length, other.available, other.band, other.max_gap,
                            other.substitution_error, other.gap_open, other.gap_extend);
        }
    };
    // Every parameter the value depends on is in the key, so a stale hit is impossible: change the
    // band or any gap parameter and it is a different entry, not a reused one. Thread-local because
    // the placement phase runs on several workers and a shared cache would need a lock in the
    // hottest path in the program; each worker paying for four entries is not worth contending for.
    static thread_local std::map<ShapeKey, double> cache;
    const ShapeKey key{read_length,   key_available,       mismatch_band,     library.max_gap_bp,
                       library.substitution_error, library.gap_open, library.gap_extend};
    const auto found = cache.find(key);
    if (found != cache.end()) return found->second;
    const double value = gap_recursion_scaled(nullptr, nullptr, 0, read_length, key_available,
                                              mismatch_band, library, false, nullptr);
    cache.emplace(key, value);
    return value;
}

MyersPattern prepare_myers_pattern(const std::string& read) {
    MyersPattern pattern;
    if (read.empty() || read.size() > kMyersWidth) return pattern;   // usable stays false
    MyersWord peq[4] = {0, 0, 0, 0};
    for (std::size_t i = 0; i < read.size(); ++i) {
        peq[two_bit_code(read[i])] |= MyersWord{1} << i;
    }
    for (int i = 0; i < 4; ++i) {
        pattern.peq_lo[i] = static_cast<std::uint64_t>(peq[i]);
        pattern.peq_hi[i] = static_cast<std::uint64_t>(peq[i] >> 64);
    }
    pattern.length = static_cast<std::uint32_t>(read.size());
    pattern.usable = true;
    return pattern;
}

namespace {
std::atomic<bool>& emission_key_prefix_enabled() {
    static std::atomic<bool> enabled{true};
    return enabled;
}
}  // namespace

void set_emission_key_prefix(bool enabled) { emission_key_prefix_enabled().store(enabled, std::memory_order_relaxed); }

// The one emission body. `key_prefix`, when given, is the hash state after `read` and its separator (taken from
// the read's own bytes by key_prefix_of, and carried by the prepared form that owns those bytes); the key it
// resumes is key_view_of's exactly. The public entry point passes none and hashes the read as before.
std::optional<GappedFit> gapped_emission_core(
    const std::string& haplotype,
    const std::string& read,
    std::size_t start,
    const LibraryModel& library,
    std::size_t mismatch_band,
    const MyersPattern* pattern,
    bool acgt_known,
    const KeyPrefix* key_prefix,
    // The prefilter's decision for exactly this (read, window) content, when the caller already holds it (the
    // union prefilter cache): -1 unknown. A known decision skips only the kernel; every logical counter moves
    // as if it had run, and executed_prefilters counts only kernels that ran.
    int known_feasible = -1,
    bool* feasible_out = nullptr,
    // I2-H: which implementation computes the emitting recursion (FragmentScoringOptions::gap_recursion_backend).
    GapRecursionBackend backend = GapRecursionBackend::scaled_reference) {
    if (read.empty()) throw std::invalid_argument("gapped emission requires a non-empty read");
    if (start > haplotype.size()) return std::nullopt;
    // ELIGIBILITY is by read length, exactly as the substitution-only model required: the haplotype
    // must be able to hold the read from here. The alignment inside that window is still free to
    // consume read_length +/- max_gap_bp, which is what admits indels.
    //
    // Admitting starts with only read_length - max_gap_bp bases left looks more permissive but is
    // wrong twice over. Such a start can be reached only by absorbing the overhang as insertions, and
    // because the normaliser is conditioned on what is available there, those gapped shapes are the
    // ONLY admissible ones and normalise to near certainty -- a gap costing 1e-4 in the interior
    // arrives at ~0.2 at the boundary. Summed over reads that manufactured exactly one extra unit of
    // density at the haplotype end, which is why the signal integrated to exposure + 1 regardless of
    // haplotype length or gap width.
    if (haplotype.size() - start < read.size()) return std::nullopt;
    // Both passes see the same bound on what the haplotype can supply from here.
    const std::size_t available = haplotype.size() - start;
    // The recursion can touch at most this much haplotype from the anchor, so this is the window the
    // refusal contract and the prefilter both apply to.
    const std::size_t window = std::min(available, read.size() + library.max_gap_bp);
    // A caller that has already established the whole haplotype and the whole read are ACGT may say so: the
    // check below then cannot fail, and repeating it cost ~11% of the union's seeding on 203 bytes per start.
    if (!acgt_known) require_acgt_window(haplotype, start, window, read);

    EmissionWork& work = emission_work();
    EmissionCache& cache = emission_cache();

    // ORDER MATTERS, and this is the correction. The exact feasibility prefilter runs FIRST, before
    // any key is constructed or looked up. It costs ~1.6us; forming a key hashes ~203 bytes and then
    // probes a shared table. Refusals are the overwhelming majority of requests (993M of 1,189M at
    // 30x), so making them pay for a key -- and then retaining them in the same representation as a
    // recursion result -- spent the whole budget on the cheap stream and collapsed the hit rate to
    // 37.8%. A refusal now returns here, having touched no table at all.
    bool feasible = known_feasible > 0;
    if (known_feasible < 0) {
        ++work.executed_prefilters;
        const bool timing = emission_timing_active().load(std::memory_order_relaxed);
        const auto prefilter_mark = timing ? std::chrono::steady_clock::now() : std::chrono::steady_clock::time_point{};
        feasible = feasible_within_band(haplotype, start, window, read, mismatch_band, pattern);
        if (timing) emission_timing().prefilter_ns += nanos_since(prefilter_mark);
    }
    if (feasible_out != nullptr) *feasible_out = feasible;
    if (!feasible) {
        ++work.prefilter_rejections;
        census_record(read, haplotype.data() + start, window, false, sizeof(double));
        return std::nullopt;
    }
    // LOGICAL: this request reached the recursion stage, whether or not the recursion runs for it.
    // Identical with and without the cache, which is what keeps the two comparable.
    ++work.dp_evaluations;

    // Only survivors get a key, and only recursion results are ever retained.
    const bool caching = cache.active.load(std::memory_order_relaxed);
    KeyView key;
    if (caching) {
        key = key_prefix != nullptr && emission_key_prefix_enabled().load(std::memory_order_relaxed)
                  ? key_view_from_prefix(*key_prefix, read, haplotype.data() + start, window)
                  : key_view_of(read, haplotype.data() + start, window);
        ++work.dp_cache_requests;
        GappedFit cached;
        switch (cache_lookup(key, &cached, work)) {
            case CacheProbe::HitRefusedAfterRecursion:
                census_record(read, haplotype.data() + start, window, true, sizeof(double));
                return std::nullopt;
            case CacheProbe::HitAccepted:
                census_record(read, haplotype.data() + start, window, true, sizeof(double));
                return cached;
            default:
                break;
        }
        if (DeferredEmissions* deferred = deferred_emissions(); deferred != nullptr) {
            if (const DeferredEmissions::Record* earlier = deferred->find(key)) {
                DeferredEmissions::Record repeat;
                repeat.key = key;
                repeat.refused = earlier->refused;
                repeat.fit = earlier->fit;
                deferred->records.push_back(repeat);
                if (repeat.refused) return std::nullopt;
                return repeat.fit;
            }
        }
    }

    ++work.executed_dps;

    GappedFit fit;
    // Filled in place rather than through a scratch vector, so an accepted placement costs one
    // allocation instead of two.
    const double raw_total = gap_recursion_emitting(haplotype, read, start, read.size(), available, mismatch_band,
                                                    library, &fit.components, backend);
    if (raw_total == kNegInf || fit.components.empty()) {
        // The recursion RAN and found nothing. That is still work a cache would avoid, so it counts
        // against the DP stream rather than the prefilter's -- and it is retained as a refusal,
        // because recomputing it costs a whole recursion.
        if (caching) {
            if (DeferredEmissions* deferred = deferred_emissions(); deferred != nullptr) {
                deferred->add(key, true, nullptr);
            } else {
                cache_store_refusal(key, CacheOutcome::RefusedAfterRecursion);
            }
        }
        census_record(read, haplotype.data() + start, window, true, sizeof(double));
        return std::nullopt;
    }
    const double shapes = log_gap_shape_mass(read.size(), available, mismatch_band, library);

    // Every component is divided by the SAME normaliser, so their sum is the whole read's emission
    // and the decomposition changes where the mass is attached without changing how much there is.
    fit.log_emission = raw_total - shapes;
    double best = kNegInf;
    for (GappedComponent& component : fit.components) {
        component.log_emission -= shapes;
        if (component.log_emission > best) {
            best = component.log_emission;
            fit.consumed_bp = component.consumed_bp;
            fit.edits = component.edits;
        }
    }
    // What a cache would hold for an accepted fit: the whole-read emission, the best consumed span
    // and edit count, and one entry per consumed-span component.
    census_record(read, haplotype.data() + start, window, true,
                  3 * sizeof(double) + fit.components.size() * (2 * sizeof(std::size_t) +
                                                                sizeof(double)));
    if (caching) {
        if (DeferredEmissions* deferred = deferred_emissions(); deferred != nullptr) {
            deferred->add(key, false, &fit);
        } else {
            cache_store_fit(key, fit);
        }
    }
    return fit;
}

std::optional<GappedFit> gapped_emission(
    const std::string& haplotype,
    const std::string& read,
    std::size_t start,
    const LibraryModel& library,
    std::size_t mismatch_band,
    const MyersPattern* pattern,
    bool acgt_known) {
    return gapped_emission_core(haplotype, read, start, library, mismatch_band, pattern, acgt_known, nullptr);
}

namespace gapped_dp_test {
double recursion(const std::string& haplotype, const std::string& read, std::size_t start, std::size_t available,
                 std::size_t band, const LibraryModel& library, ComponentList* components, GapRecursionBackend backend) {
    return gap_recursion_emitting(haplotype, read, start, read.size(), available, band, library, components, backend);
}

std::optional<GappedFit> emission_known_feasible(const std::string& haplotype, const std::string& read,
                                                 std::size_t start, const LibraryModel& library, std::size_t band,
                                                 GapRecursionBackend backend) {
    return gapped_emission_core(haplotype, read, start, library, band, nullptr, true, nullptr, 1, nullptr, backend);
}

std::optional<GappedFit> emission(const std::string& haplotype, const std::string& read, std::size_t start,
                                  const LibraryModel& library, std::size_t band, GapRecursionBackend backend) {
    return gapped_emission_core(haplotype, read, start, library, band, nullptr, false, nullptr, -1, nullptr, backend);
}

bool fixed_geometry_supported(std::size_t band, const LibraryModel& library) {
    return fixed_recursion_for(library.max_gap_bp, band) != nullptr;
}

std::uint64_t fixed_geometry_executions() { return t_fixed_recursions; }

bool feasible(const std::string& haplotype, const std::string& read, std::size_t start, const LibraryModel& library,
              std::size_t band) {
    if (read.empty() || start > haplotype.size() || haplotype.size() - start < read.size()) return false;
    const std::size_t window = std::min(haplotype.size() - start, read.size() + library.max_gap_bp);
    return feasible_within_band(haplotype, start, window, read, band, nullptr);
}
}  // namespace gapped_dp_test

bool emission_key_prefix_matches(const std::string& read, const std::string& window) {
    const KeyView full = key_view_of(read, window.data(), window.size());
    const KeyView resumed = key_view_from_prefix(key_prefix_of(read), read, window.data(), window.size());
    return full.h0 == resumed.h0 && full.h1 == resumed.h1 && full.read == resumed.read &&
           full.read_length == resumed.read_length && full.window == resumed.window &&
           full.window_length == resumed.window_length;
}

std::optional<GappedFit> gapped_emission_reference(
    const std::string& haplotype,
    const std::string& read,
    std::size_t start,
    const LibraryModel& library,
    std::size_t mismatch_band) {
    if (read.empty()) throw std::invalid_argument("gapped emission requires a non-empty read");
    if (start > haplotype.size()) return std::nullopt;
    if (haplotype.size() - start < read.size()) return std::nullopt;
    const std::size_t available = haplotype.size() - start;
    // No prefilter and no memoised normaliser: this path exists to be believed without depending on
    // any of the machinery it is used to check, so it recomputes the shape mass in log space too,
    // unclamped. It keeps the reference recursion's own non-ACGT refusal.
    ComponentList raw;
    const double raw_total = gap_recursion_logspace(&haplotype, &read, start, read.size(), available,
                                                    mismatch_band, library, true, &raw);
    if (raw_total == kNegInf || raw.empty()) return std::nullopt;
    const double shapes = gap_recursion_logspace(nullptr, nullptr, 0, read.size(), available,
                                                 mismatch_band, library, false, nullptr);

    GappedFit fit;
    fit.log_emission = raw_total - shapes;
    fit.components.reserve(raw.size());
    double best = kNegInf;
    for (GappedComponent component : raw) {
        component.log_emission -= shapes;
        if (component.log_emission > best) {
            best = component.log_emission;
            fit.consumed_bp = component.consumed_bp;
            fit.edits = component.edits;
        }
        fit.components.push_back(component);
    }
    return fit;
}

std::optional<GappedFit> gapped_emission_ending_at(
    const std::string& haplotype,
    const std::string& read,
    std::size_t end,
    const LibraryModel& library,
    std::size_t mismatch_band) {
    if (read.empty()) throw std::invalid_argument("gapped emission requires a non-empty read");
    if (end > haplotype.size()) return std::nullopt;
    // Eligibility by read length, measured leftward: the haplotype must be able to hold the read
    // ending here, exactly as the left-anchored form requires it to hold the read starting there.
    if (end < read.size()) return std::nullopt;
    // Reversing BOTH sequences turns "read leftward from `end`" into "read rightward from
    // haplotype.size() - end", which is the recursion already written. Reversing both preserves the
    // match/mismatch pattern exactly; what changes is the direction the gap model runs in, and that
    // is the whole point -- the mate really is sequenced from this end.
    const std::string reversed_haplotype(haplotype.rbegin(), haplotype.rend());
    const std::string reversed_read(read.rbegin(), read.rend());
    return gapped_emission(reversed_haplotype, reversed_read, haplotype.size() - end, library,
                           mismatch_band);
}

std::optional<std::size_t> fr_insert_from_endpoints(
    std::size_t forward_start,
    std::size_t forward_consumed_bp,
    std::size_t reverse_end,
    std::size_t reverse_consumed_bp,
    const InsertSupport& support) {
    if (reverse_end <= forward_start) return std::nullopt;
    const std::size_t insert = reverse_end - forward_start;
    // Overlapping mates are a physical state, so the floor is the longer mate's CONSUMED span,
    // never the sum. A gapped mate occupies a different stretch than its read length, so the floor
    // has to be stated on what was consumed.
    if (insert < std::max(forward_consumed_bp, reverse_consumed_bp)) return std::nullopt;
    if (support.empty() || insert < support.min_bp || insert > support.max_bp) return std::nullopt;
    return insert;
}

std::optional<std::size_t> fr_insert_length(
    std::size_t forward_start,
    std::size_t forward_length,
    std::size_t reverse_start,
    std::size_t reverse_length,
    const InsertSupport& support) {
    // The reverse mate's own right edge closes the template. A reverse mate that ends before the
    // forward mate begins is not an FR pair at all.
    const std::size_t reverse_end = reverse_start + reverse_length;
    if (reverse_end <= forward_start) return std::nullopt;
    const std::size_t insert = reverse_end - forward_start;
    // Overlapping mates are a physical state, so the floor is the longer mate, not their sum. A
    // run-wide maximum read length here would wrongly raise the floor for shorter fragments.
    if (insert < std::max(forward_length, reverse_length)) return std::nullopt;
    if (support.empty() || insert < support.min_bp || insert > support.max_bp) return std::nullopt;
    return insert;
}

double haplotype_exposure(
    std::size_t haplotype_length,
    const DiscreteInsertDistribution& insert,
    std::size_t physical_floor_bp) {
    if (insert.support.empty()) return 0.0;
    double exposure = 0.0;
    const std::size_t from = std::max(insert.support.min_bp, physical_floor_bp);
    for (std::size_t length = from; length <= insert.support.max_bp; ++length) {
        if (haplotype_length < length) continue;
        const double probability = std::exp(insert.log_probability(length));
        exposure += probability * static_cast<double>(haplotype_length - length + 1);
    }
    return exposure;
}

double log_fragment_contribution(
    double log_mass,
    double log_background,
    const LibraryModel& library) {
    validate_library_model(library);
    const double signal = log_mass == kNegInf
        ? kNegInf
        : std::log1p(-library.outlier_weight) +
              std::log(library.fragment_start_intensity) + log_mass;
    if (library.outlier_weight == 0.0) {
        if (signal == kNegInf) {
            throw std::invalid_argument(
                "a fragment with no accepted placement and no outlier mixture makes every candidate "
                "impossible; declare a positive outlier weight or refuse the fragment upstream");
        }
        return signal;
    }
    return log_add_exp(signal, std::log(library.outlier_weight) + log_background);
}

double expected_signal_fragments(double exposure_sum, const LibraryModel& library) {
    validate_library_model(library);
    return (1.0 - library.outlier_weight) * library.fragment_start_intensity * exposure_sum;
}

namespace {
// One audit library's combination constants: each is the exact subexpression log_fragment_contribution and
// expected_signal_fragments evaluate, in the same association, so an audit score is bit for bit the score a run
// declaring that library computes.
struct AuditCombination {
    double log_signal_scale = 0.0;   // log1p(-w) + log(lambda)
    double log_outlier = 0.0;        // log(w)
    bool outlier_free = false;
    double exposure_scale = 0.0;     // (1 - w) * lambda
};
constexpr std::size_t kMaxAuditLibraries = 4;
std::vector<AuditCombination> audit_combinations(const std::vector<LibraryModel>& libraries) {
    if (libraries.size() > kMaxAuditLibraries) throw std::invalid_argument("at most four audit libraries");
    std::vector<AuditCombination> out;
    for (const LibraryModel& library : libraries) {
        validate_library_model(library);
        AuditCombination c;
        c.log_signal_scale = std::log1p(-library.outlier_weight) + std::log(library.fragment_start_intensity);
        c.outlier_free = library.outlier_weight == 0.0;
        c.log_outlier = c.outlier_free ? 0.0 : std::log(library.outlier_weight);
        c.exposure_scale = (1.0 - library.outlier_weight) * library.fragment_start_intensity;
        out.push_back(c);
    }
    return out;
}
// One fragment's contribution from prepared constants: exactly log_fragment_contribution's expressions, in its
// association (signal = (log1p(-w) + log(lambda)) + mass; log(w) + background; the zero-outlier refusal).
inline double prepared_contribution(const AuditCombination& c, double combined, double background) {
    const double signal = combined == kNegInf ? kNegInf : c.log_signal_scale + combined;
    if (c.outlier_free) {
        if (signal == kNegInf) {
            throw std::invalid_argument(
                "a fragment with no accepted placement and no outlier mixture makes every candidate "
                "impossible; declare a positive outlier weight or refuse the fragment upstream");
        }
        return signal;
    }
    return log_add_exp(signal, c.log_outlier + background);
}
// I2-F: the tally of what prepared main-library terms removed (one update per combination call)
struct PreparedTally {
    std::atomic<std::uint64_t> combinations{0}, terms{0}, hypothetical_logs{0}, prepared_logs{0};
    std::atomic<std::uint64_t> hypothetical_validations{0}, prepared_validations{0};
};
PreparedTally& prepared_tally() {
    static PreparedTally t;
    return t;
}
// I2-G correction (see PreparedLibraryTermsTally). Called once per prepared call, after the combination, with the call's
// candidate and fragment counts and the summed fragments_with_mass of its candidates.
void note_prepared_combination(std::uint64_t candidates, std::uint64_t fragments, std::uint64_t with_mass,
                               const AuditCombination& main) {
    const std::uint64_t terms = candidates * fragments;
    PreparedTally& t = prepared_tally();
    t.combinations.fetch_add(1, std::memory_order_relaxed);
    t.terms.fetch_add(terms, std::memory_order_relaxed);
    // scalar: log1p(-w) + log(lambda) where the combined mass is finite; log(w) on every term when w > 0
    t.hypothetical_logs.fetch_add(2 * with_mass + (main.outlier_free ? 0 : terms), std::memory_order_relaxed);
    t.prepared_logs.fetch_add(main.outlier_free ? 2 : 3, std::memory_order_relaxed);
    // scalar: one validation per term (log_fragment_contribution) and one per candidate (expected_signal_fragments)
    t.hypothetical_validations.fetch_add(terms + candidates, std::memory_order_relaxed);
    t.prepared_validations.fetch_add(1, std::memory_order_relaxed);
}
// One candidate's diploid combination under the run's library and every audit library, in one pass over its two
// mass columns. The run's own score is computed by exactly the expressions of the audit-free loop. `main` (I2-F): the
// run's library prepared once; null: the scalar reference functions per term.
void combine_with_audits(const std::vector<double>& mass_a, const std::vector<double>& mass_b,
                         const std::vector<double>& background, std::size_t fragments, double exposure_sum,
                         const LibraryModel& library, const std::vector<AuditCombination>& audits, DiplotypeScore& score,
                         std::vector<std::vector<DiplotypeScore>>& audit_scores, std::size_t candidate,
                         const AuditCombination* main = nullptr) {
    const std::size_t views = audits.size();
    score.exposure_term = main != nullptr ? -(main->exposure_scale * exposure_sum) : -expected_signal_fragments(exposure_sum, library);
    score.absolute_term_sum = std::abs(score.exposure_term);
    double term[kMaxAuditLibraries] = {};
    double absolute[kMaxAuditLibraries] = {};
    for (std::size_t k = 0; k < views; ++k) {
        DiplotypeScore& audit = audit_scores[k][candidate];
        audit.exposure_term = -(audits[k].exposure_scale * exposure_sum);
        absolute[k] = std::abs(audit.exposure_term);
    }
    double fragment_term = 0.0;
    for (std::size_t f = 0; f < fragments; ++f) {
        const double combined = log_add_exp(mass_a[f], mass_b[f]);
        if (combined != kNegInf) ++score.fragments_with_mass;
        const double contribution = main != nullptr ? prepared_contribution(*main, combined, background[f])
                                                    : log_fragment_contribution(combined, background[f], library);
        fragment_term += contribution;
        score.absolute_term_sum += std::abs(contribution);
        for (std::size_t k = 0; k < views; ++k) {
            const double value = prepared_contribution(audits[k], combined, background[f]);
            term[k] += value;
            absolute[k] += std::abs(value);
        }
    }
    score.fragment_term = fragment_term;
    score.log_likelihood = score.exposure_term + fragment_term;
    for (std::size_t k = 0; k < views; ++k) {
        DiplotypeScore& audit = audit_scores[k][candidate];
        audit.fragment_term = term[k];
        audit.log_likelihood = audit.exposure_term + term[k];
        audit.absolute_term_sum = absolute[k];
        audit.fragments_with_mass = score.fragments_with_mass;
    }
}
} // namespace

PreparedLibraryTermsTally prepared_library_terms_tally() {
    const PreparedTally& t = prepared_tally();
    PreparedLibraryTermsTally out;
    out.combinations = t.combinations.load();
    out.terms = t.terms.load();
    out.hypothetical_scalar_logs = t.hypothetical_logs.load();
    out.prepared_logs = t.prepared_logs.load();
    out.net_logs_avoided = static_cast<std::int64_t>(out.hypothetical_scalar_logs) - static_cast<std::int64_t>(out.prepared_logs);
    out.hypothetical_scalar_validations = t.hypothetical_validations.load();
    out.prepared_validations = t.prepared_validations.load();
    out.net_validations_avoided =
        static_cast<std::int64_t>(out.hypothetical_scalar_validations) - static_cast<std::int64_t>(out.prepared_validations);
    return out;
}
void reset_prepared_library_terms_tally() {
    PreparedTally& t = prepared_tally();
    t.combinations = 0;
    t.terms = 0;
    t.hypothetical_logs = 0;
    t.prepared_logs = 0;
    t.hypothetical_validations = 0;
    t.prepared_validations = 0;
}

std::vector<FragmentPlacement> enumerate_placements_exhaustive(
    const std::string& haplotype,
    const ReadPair& fragment,
    const InsertSupport& support,
    const LibraryModel& library,
    std::size_t mismatch_band) {
    validate_read_pair(fragment);
    validate_library_model(library);
    std::vector<FragmentPlacement> placements;
    if (haplotype.empty() || support.empty()) return placements;

    const std::string first_rc = reverse_complement_read(fragment.first.bases);
    const std::string second_rc = reverse_complement_read(fragment.second.bases);

    for (int arm = 0; arm < 2; ++arm) {
        const MateOrientation orientation =
            arm == 0 ? MateOrientation::FirstForward : MateOrientation::SecondForward;
        // Whichever mate reads off the forward strand is compared as given; the other is compared
        // pre-complemented, which makes both comparisons forward.
        const std::string& forward_bases =
            arm == 0 ? fragment.first.bases : fragment.second.bases;
        const std::string& reverse_bases = arm == 0 ? second_rc : first_rc;

        // Each anchor's components once per form. The pairing below reads them O(L^2) times, and
        // recomputing the recursion inside that loop made the oracle quadratic in DP evaluations
        // rather than merely in pairings.
        //
        // The two mates are anchored at OPPOSITE ends, because that is how the template was read:
        // the forward mate rightward from the left endpoint, the reverse mate leftward from the
        // right endpoint. Anchoring both at a start is what broke the exposure identity.
        std::vector<std::optional<GappedFit>> forward_fits(haplotype.size() + 1);
        std::vector<std::optional<GappedFit>> reverse_fits(haplotype.size() + 1);
        for (std::size_t at = 0; at <= haplotype.size(); ++at) {
            if (at < haplotype.size()) {
                forward_fits[at] = gapped_emission(haplotype, forward_bases, at, library,
                                                   mismatch_band);
            }
            reverse_fits[at] =
                gapped_emission_ending_at(haplotype, reverse_bases, at, library, mismatch_band);
        }

        // Keyed by (left endpoint, insert). The right endpoint is their sum, so this names the
        // template's two ends -- which is what a physical origin is. Every latent alignment pair
        // arriving at the same pair of ends is one more explanation of ONE template, so it is summed
        // into that origin rather than producing a second one.
        std::map<std::pair<std::size_t, std::size_t>, OriginAccumulator> origins;
        for (std::size_t forward_start = 0; forward_start < haplotype.size(); ++forward_start) {
            const std::optional<GappedFit>& forward_fit = forward_fits[forward_start];
            if (!forward_fit.has_value()) continue;
            for (std::size_t reverse_end = 0; reverse_end <= haplotype.size(); ++reverse_end) {
                const std::optional<GappedFit>& reverse_fit = reverse_fits[reverse_end];
                if (!reverse_fit.has_value()) continue;
                // Geometry follows the bases each mate CONSUMES, not the bases it carries, and each
                // consumed-span class is given its own geometry: a gapped mate spans a different
                // stretch of haplotype than its read length, and a different one per class.
                for (const GappedComponent& forward_span : forward_fit->components) {
                    if (forward_span.log_emission == kNegInf) continue;
                    for (const GappedComponent& reverse_span : reverse_fit->components) {
                        if (reverse_span.log_emission == kNegInf) continue;
                        const std::optional<std::size_t> insert = fr_insert_from_endpoints(
                            forward_start, forward_span.consumed_bp,
                            reverse_end, reverse_span.consumed_bp, support);
                        if (!insert.has_value()) continue;
                        origins[{forward_start, *insert}].add(
                            forward_span.log_emission + reverse_span.log_emission,
                            reverse_end - reverse_span.consumed_bp,
                            forward_span.consumed_bp, reverse_span.consumed_bp,
                            forward_span.edits, reverse_span.edits);
                    }
                }
            }
        }

        for (const auto& item : origins) {
            FragmentPlacement placement;
            placement.orientation = orientation;
            placement.forward_start = item.first.first;
            placement.insert_bp = item.first.second;
            placement.log_latent_mass = item.second.log_latent_mass;
            placement.latent_alignments = item.second.latent_alignments;
            placement.reverse_start = item.second.reverse_start;
            placement.forward_consumed_bp = item.second.forward_consumed_bp;
            placement.reverse_consumed_bp = item.second.reverse_consumed_bp;
            placement.forward_mismatches = item.second.forward_mismatches;
            placement.reverse_mismatches = item.second.reverse_mismatches;
            placements.push_back(placement);
        }
    }
    return placements;
}

double log_placement_mass(
    const std::vector<FragmentPlacement>& placements,
    const ReadPair& fragment,
    const DiscreteInsertDistribution& insert,
    const LibraryModel& library) {
    validate_read_pair(fragment);
    validate_library_model(library);
    double mass = kNegInf;
    for (const FragmentPlacement& placement : placements) {
        // One template, so one orientation prior and one insert probability. The latent alignments
        // were summed into log_latent_mass at the geometry each of them implied, which is the whole
        // point of carrying them separately that far.
        const double term =
            kLogMateOrientationPrior +
            insert.log_probability(placement.insert_bp) +
            placement.log_latent_mass;
        mass = log_add_exp(mass, term);
    }
    return mass;
}


// -------------------------------------------------------------------------------------------
// Accelerated scorer.
// -------------------------------------------------------------------------------------------

// Named rather than anonymous: PreparedScoringContext::State holds these, and a class with external
// linkage may not have members whose types come from the anonymous namespace.
namespace scoring_detail {

constexpr std::uint64_t kNoCode = std::numeric_limits<std::uint64_t>::max();

std::uint64_t two_bit(char base) {
    switch (base) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return kNoCode;
    }
}

// One oriented copy of one mate, already turned so that it matches the haplotype's forward strand.
struct ReadForm {
    std::size_t fragment = 0;
    std::string bases;
    bool is_forward_mate_of_first_arm = false;   // bookkeeping only; arms are resolved by role below
    // Built once during fragment preparation. Rebuilding the masks inside the hot feasibility test
    // cost a pass over the read on every one of 1.19 billion verifications.
    MyersPattern pattern;
    bool all_acgt = false;   // every base is ACGT, so no window check can fail on the read's side
    // The emission key's hash state after these bytes: lives and dies with them, so it can never describe
    // other bytes. The cache keeps its own copy of every key's bytes, so its lifetime does not depend on this.
    KeyPrefix key_prefix{};
    // Equal for two forms exactly when their bytes are equal (interned at preparation by byte comparison).
    std::uint32_t content_id = 0;
};

// piece code -> (form, offset of the piece inside that form)
struct PieceEntry {
    std::uint64_t code = 0;
    std::uint32_t form = 0;
    std::uint32_t offset = 0;
};

struct PieceIndex {
    std::size_t piece_bp = 0;
    std::vector<PieceEntry> entries;   // sorted by code
};

// One verified consumed-span class of one read form at one start on one haplotype. There is a row
// per class, not per start: the classes end the mate at different coordinates and so pair with
// different reverse mates, and merging them before the geometry is applied is exactly the defect
// this file was corrected for. The emission is carried rather than a mismatch count, because under
// the indel-aware model the count no longer determines it.
struct FormPlacement {
    std::uint32_t start = 0;
    std::uint32_t consumed = 0;        // haplotype bases this alignment class spans
    std::uint32_t mismatches = 0;      // diagnostic only
    double emission = 0.0;             // normalised log mass of this class
};

// The exact statistical class of a physical origin. Origins sharing every field may be stored once
// with a multiplicity, but ONLY because every term of their mass is then identical. Keying on
// mismatch counts was sufficient when the count determined the emission; it no longer does, so the
// summed latent mass is itself part of the key. Two origins with equal counts but different gapped
// alignments are different classes and must not be collapsed.
//
// The multiplicity here counts PLACES -- distinct left endpoints carrying the same insert and the
// same latent mass, which is what repeat copies look like. It never counts alternative alignments of
// one place: those were summed into log_latent_mass before this key was formed.
struct StateKey {
    std::uint32_t insert_bp = 0;
    std::uint32_t orientation = 0;
    double log_latent_mass = 0.0;

    bool operator<(const StateKey& other) const noexcept {
        return std::tie(insert_bp, orientation, log_latent_mass) <
               std::tie(other.insert_bp, other.orientation, other.log_latent_mass);
    }
};

// Templated on the form container so the one-shot path, the census and the persistent context can
// each hold forms in whatever container their lifetime requires -- a vector where the set is built
// once and discarded, a deque where the emission cache keeps pointers into the bytes.
template <typename Forms>
std::vector<PieceIndex> build_piece_indexes(
    const Forms& forms,
    std::size_t band,
    std::size_t minimum_piece_bp) {
    // Reads of different lengths yield pieces of different lengths, and one rolling scan can only
    // carry one length. Bucket by piece length and scan once per bucket.
    std::map<std::size_t, std::vector<PieceEntry>> by_length;
    for (std::size_t f = 0; f < forms.size(); ++f) {
        const std::size_t length = forms[f].bases.size();
        // The band guarantees some contiguous 1/(band+1) chunk is exact. Keying on a PREFIX of that
        // chunk keeps the guarantee, because an exact chunk has an exact prefix, and caps the key at
        // what a 2-bit code can hold.
        const std::size_t chunk = length / (band + 1);
        const std::size_t piece = std::min<std::size_t>(31, chunk);
        if (piece < minimum_piece_bp) {
            throw std::invalid_argument(
                "the declared mismatch band splits a read into pieces too short to index; widen the "
                "reads, narrow the band, or lower minimum_piece_bp deliberately");
        }
        for (std::size_t p = 0; p <= band; ++p) {
            const std::size_t offset = p * chunk;
            std::uint64_t code = 0;
            bool usable = true;
            for (std::size_t i = 0; i < piece; ++i) {
                const std::uint64_t bits = two_bit(forms[f].bases[offset + i]);
                if (bits == kNoCode) { usable = false; break; }
                code = (code << 2) | bits;
            }
            if (!usable) continue;
            by_length[piece].push_back({code, static_cast<std::uint32_t>(f),
                                        static_cast<std::uint32_t>(offset)});
        }
    }
    std::vector<PieceIndex> indexes;
    for (auto& item : by_length) {
        std::sort(item.second.begin(), item.second.end(),
                  [](const PieceEntry& lhs, const PieceEntry& rhs) {
                      return std::tie(lhs.code, lhs.form, lhs.offset) <
                             std::tie(rhs.code, rhs.form, rhs.offset);
                  });
        indexes.push_back({item.first, std::move(item.second)});
    }
    return indexes;
}

// Process-wide budget for the two caps that bound WORK rather than merely decide acceptance. Shared by
// every placement worker, so T workers cannot each consume a full cap before anyone refuses -- which is
// what a per-worker guard silently allowed, and what makes a 64-thread request nearly breaker-free.
//
// Crossing a cap here DOES decide the answer: the worker that observes it refuses immediately, which is
// sound because a shared total can only undercount the exact one. Acceptance is the half decided later,
// on the exact merged totals. `cancelled` is what makes the refusal prompt: the worker that crosses a
// cap raises it, and every other worker abandons its range at the next check -- in the piece-hit
// publication, at each physical origin, and at the haplotype, fragment, sort and verification
// boundaries -- rather than running to completion.
struct PlacementBudget {
    std::atomic<std::uint64_t> piece_hits{0};
    std::atomic<std::uint64_t> physical_origins{0};
    std::atomic<bool> cancelled{false};
    std::uint64_t max_piece_hits = 0;
    std::uint64_t max_physical_origins = 0;
    std::uint64_t publish_batch = 4096;

    bool stop() const { return cancelled.load(std::memory_order_relaxed); }
    void cancel() { cancelled.store(true, std::memory_order_relaxed); }
};

// Raised by a worker abandoning its range because another worker already crossed a cap. It never
// escapes score_diplotypes: the cap that was actually crossed is reported instead, so the message a
// caller sees does not depend on which worker noticed first.
class PlacementCancelled : public std::exception {};

// ---- MEASUREMENT ONLY: the necessary-condition piece screen (see the lemma in the header) ----------
constexpr std::uint64_t kNoPieceCode = std::numeric_limits<std::uint64_t>::max();

struct PieceLayout {
    std::size_t pieces = 0;
    bool usable = false;   // false: the screen cannot decide, so it admits
    std::vector<std::uint32_t> offset;
    std::vector<std::uint32_t> length;
    std::vector<std::uint64_t> code;
};

PieceLayout piece_layout(const std::string& read, std::size_t pieces) {
    PieceLayout layout;
    layout.pieces = pieces;
    const std::size_t n = read.size();
    if (pieces == 0 || pieces > n) return layout;
    for (std::size_t j = 0; j < pieces; ++j) {
        const std::size_t begin = j * n / pieces;
        const std::size_t end = (j + 1) * n / pieces;
        if (end - begin == 0 || end - begin > 31) return layout;
        std::uint64_t code = 0;
        for (std::size_t i = begin; i < end; ++i) {
            const std::uint64_t bits = two_bit(read[i]);
            if (bits == kNoCode) return layout;
            code = (code << 2) | bits;
        }
        layout.offset.push_back(static_cast<std::uint32_t>(begin));
        layout.length.push_back(static_cast<std::uint32_t>(end - begin));
        layout.code.push_back(code);
    }
    layout.usable = true;
    return layout;
}

// Exact 2-bit codes of every `length`-substring of the text, built once per length and text.
struct HaplotypeCodes {
    const std::string* text = nullptr;
    std::map<std::size_t, std::vector<std::uint64_t>> by_length;
    std::uint64_t build_ns = 0;
    const std::vector<std::uint64_t>& of(std::size_t length) {
        auto found = by_length.find(length);
        if (found != by_length.end()) return found->second;
        const auto mark = std::chrono::steady_clock::now();
        std::vector<std::uint64_t> codes(text->size(), kNoPieceCode);
        const std::uint64_t mask = (std::uint64_t{1} << (2 * length)) - 1;
        std::uint64_t code = 0;
        std::size_t run = 0;
        for (std::size_t at = 0; at < text->size(); ++at) {
            const std::uint64_t bits = two_bit((*text)[at]);
            if (bits == kNoCode) { run = 0; code = 0; continue; }
            code = ((code << 2) | bits) & mask;
            if (++run >= length) codes[at + 1 - length] = code;
        }
        build_ns += nanos_since(mark);
        return by_length.emplace(length, std::move(codes)).first->second;
    }
};

enum class ScreenVerdict { CountReject, ChainReject, Admit };

// The screen for a start whose alignment may consume haplotype[start, start + window).
ScreenVerdict piece_screen(HaplotypeCodes& codes, std::size_t start, std::size_t window,
                           const PieceLayout& layout, std::size_t k) {
    if (!layout.usable) return ScreenVerdict::Admit;
    const std::size_t d = layout.pieces;
    if (d <= k) return ScreenVerdict::Admit;   // the lemma says nothing
    const std::size_t need = d - k;
    struct Node {
        std::uint32_t piece;
        int delta;
    };
    static thread_local std::vector<Node> nodes;
    static thread_local std::vector<int> best;
    nodes.clear();
    const long K = static_cast<long>(k);
    std::size_t with_match = 0;
    for (std::size_t j = 0; j < d; ++j) {
        const std::vector<std::uint64_t>& table = codes.of(layout.length[j]);
        bool any = false;
        for (long delta = -K; delta <= K; ++delta) {
            const long at = static_cast<long>(start) + static_cast<long>(layout.offset[j]) + delta;
            if (at < static_cast<long>(start)) continue;
            if (static_cast<std::size_t>(at) + layout.length[j] > start + window) continue;
            if (table[static_cast<std::size_t>(at)] != layout.code[j]) continue;
            nodes.push_back({static_cast<std::uint32_t>(j), static_cast<int>(delta)});
            any = true;
        }
        if (any) ++with_match;
        if (with_match + (d - j - 1) < need) return ScreenVerdict::CountReject;
    }
    if (with_match < need) return ScreenVerdict::CountReject;
    // Longest offset-consistent chain within the displacement budget: best[n * (k+1) + c] is the most
    // exact pieces in a chain ending at node n having spent c units of displacement.
    const std::size_t width = k + 1;
    best.assign(nodes.size() * width, -1);
    for (std::size_t b = 0; b < nodes.size(); ++b) {
        const std::size_t first = static_cast<std::size_t>(std::abs(nodes[b].delta));
        if (first <= k) best[b * width + first] = 1;
        for (std::size_t a = 0; a < b; ++a) {
            if (nodes[a].piece >= nodes[b].piece) continue;
            const long between = static_cast<long>(layout.offset[nodes[b].piece]) -
                                 static_cast<long>(layout.offset[nodes[a].piece] + layout.length[nodes[a].piece]);
            const long change = static_cast<long>(nodes[b].delta) - static_cast<long>(nodes[a].delta);
            if (change < -between) continue;   // positions would overlap or go backwards
            const std::size_t step = static_cast<std::size_t>(change < 0 ? -change : change);
            for (std::size_t c = 0; c + step <= k; ++c) {
                const int from = best[a * width + c];
                if (from < 0) continue;
                int& to = best[b * width + c + step];
                if (from + 1 > to) to = from + 1;
            }
        }
        for (std::size_t c = 0; c <= k; ++c) {
            if (best[b * width + c] >= static_cast<int>(need)) return ScreenVerdict::Admit;
        }
    }
    return ScreenVerdict::ChainReject;
}

std::atomic<bool>& piece_screen_shadow_active() {
    static std::atomic<bool> active{false};
    return active;
}
std::mutex& piece_screen_mutex() {
    static std::mutex mutex;
    return mutex;
}
PieceScreenShadowReport& piece_screen_totals() {
    static PieceScreenShadowReport totals;
    return totals;
}

// ---- exact union placement reuse -------------------------------------------------------------------
struct CachedPlacement {
    std::uint32_t form = 0;
    std::uint32_t start = 0;   // local to the segment
    std::uint32_t consumed = 0;
    std::uint32_t mismatches = 0;
    double emission = 0.0;
};
struct UnionSegmentEntry {
    std::vector<CachedPlacement> placements;   // the interior zone's, in (form, start) order
};
// Starts [lo, hi] of one haplotype strand are served by `entry`, translated by `shift`.
struct ReplayZone {
    std::uint32_t lo = 0;
    std::uint32_t hi = 0;
    std::uint32_t shift = 0;
    const UnionSegmentEntry* entry = nullptr;
};
struct UnionReplay {
    std::vector<std::vector<ReplayZone>> forward, reverse;   // per haplotype
    std::atomic<std::uint64_t> consumed{0};
    std::uint64_t max_consumed = 0;   // 0: none
};
std::atomic<bool>& union_hash_collapse() {
    static std::atomic<bool> collapse{false};
    return collapse;
}

// Proposed (form, start) pairs sorted ascending and deduplicated, without a comparison sort over the whole list:
// an in-place counting permutation by form (the forms are dense indices; American-flag cycle leading, so no
// second copy of the list is held), then each form's few starts sorted and made unique in place. The result is
// exactly std::sort + std::unique's -- the same pairs in the same order -- and depends on nothing but the
// multiset of proposals. `cursor` and `bound` are worker storage, one entry per form.
void sort_unique_proposals_bucketed(std::vector<std::pair<std::uint32_t, std::uint32_t>>& proposed,
                                    std::size_t form_count, std::vector<std::uint32_t>& cursor,
                                    std::vector<std::uint32_t>& bound) {
    cursor.assign(form_count, 0);
    for (const auto& item : proposed) ++cursor[item.first];
    bound.resize(form_count);
    std::uint32_t running = 0;
    for (std::size_t form = 0; form < form_count; ++form) {
        const std::uint32_t size = cursor[form];
        cursor[form] = running;   // the next unfilled slot of this form's bucket
        running += size;
        bound[form] = running;    // its end
    }
    for (std::size_t form = 0; form < form_count; ++form) {
        while (cursor[form] < bound[form]) {
            auto value = proposed[cursor[form]];
            if (value.first == form) {
                ++cursor[form];
                continue;
            }
            // Carry the displaced pair to its own bucket, and whatever it displaces onward, until one lands here.
            do {
                std::swap(value, proposed[cursor[value.first]++]);
            } while (value.first != form);
            proposed[cursor[form]++] = value;
        }
    }
    std::size_t out = 0;
    std::size_t begin = 0;
    for (std::size_t form = 0; form < form_count; ++form) {
        const std::size_t end = bound[form];
        if (begin == end) continue;
        auto* first = proposed.data() + begin;
        auto* last = proposed.data() + end;
        if (end - begin <= 24) {
            for (auto* it = first + 1; it < last; ++it) {   // insertion sort: a form holds a few dozen starts
                const auto value = *it;
                auto* hole = it;
                while (hole > first && value.second < (hole - 1)->second) {
                    *hole = *(hole - 1);
                    --hole;
                }
                *hole = value;
            }
        } else {
            std::sort(first, last, [](const auto& a, const auto& b) { return a.second < b.second; });
        }
        // Compaction: out <= begin always, so a form's pairs never overwrite ones not yet read.
        for (auto* it = first; it < last; ++it) {
            if (out != 0 && proposed[out - 1] == *it) continue;
            proposed[out++] = *it;
        }
        begin = end;
    }
    proposed.resize(out);
}

// ---- the union prefilter cache ---------------------------------------------------------------------------
//
// The prefilter's decision is a function of the read's bytes and the window's bytes alone (the window's length,
// min(available, read + max_gap_bp), carries the end condition). The census found the union's direct stream
// made of few such pairs, each requested ~50 times: remembered once, a pair's decision is exact every time.
// Identity is by ID: read contents interned at preparation, window contents interned here -- each distinct
// window's bytes copied into this index's own storage, found by a hash that only chooses where to look and
// confirmed by length and complete byte comparison. The decision table is keyed by the two IDs, exact integers.
struct UnionPrefilterCache {
    struct Stored {
        std::uint64_t hash = 0;
        const char* bytes = nullptr;
        std::uint32_t length = 0;
        std::uint32_t id = 0;
    };
    struct Shard {
        std::mutex guard;
        std::vector<Stored> table;   // open addressing, at most half full
        std::size_t occupied = 0;
        std::vector<std::unique_ptr<char[]>> slabs;
        std::size_t slab_used = 0;
        std::uint64_t bytes = 0;
    };
    static constexpr std::size_t kShards = 256;
    static constexpr std::size_t kWindowBits = 27;
    std::array<Shard, kShards> shards;
    std::atomic<std::uint32_t> next_window{0};
    std::size_t read_length = 0;
    std::size_t window_length = 0;
    // Decisions: ((key + 1) << 1) | feasible, key = read content << kWindowBits | window content; 0 is empty.
    std::unique_ptr<std::atomic<std::uint64_t>[]> slots;
    std::uint64_t slot_bits = 0;
    std::atomic<std::uint64_t> held{0};
    std::uint64_t max_held = 0;   // half the slots: admission stops there, and misses run the kernel
    std::atomic<std::uint64_t> not_retained{0};
    std::atomic<std::uint64_t> consumed{0};   // hits, published in batches of publish_batch
    std::uint64_t max_consumed = 0;
    std::uint64_t publish_batch = 4096;
    std::atomic<std::uint64_t> misses{0};
    std::atomic<std::uint64_t> positions_indexed{0};
    std::atomic<std::uint64_t> index_nanos{0};
};
std::atomic<bool>& union_prefilter_hash_collapse() {
    static std::atomic<bool> collapse{false};
    return collapse;
}

std::uint32_t window_content_id(UnionPrefilterCache& cache, const char* bytes, std::size_t length, std::uint64_t hash) {
    if (union_prefilter_hash_collapse().load(std::memory_order_relaxed)) hash = 0;
    UnionPrefilterCache::Shard& shard = cache.shards[hash >> 56];
    const std::lock_guard<std::mutex> guard(shard.guard);
    if (shard.table.empty() || (shard.occupied + 1) * 2 > shard.table.size()) {
        std::vector<UnionPrefilterCache::Stored> next(shard.table.empty() ? 1024 : shard.table.size() * 2);
        const std::size_t mask = next.size() - 1;
        for (const UnionPrefilterCache::Stored& stored : shard.table) {
            if (stored.bytes == nullptr) continue;
            std::size_t at = (stored.hash * 0x9e3779b97f4a7c15ULL) & mask;
            while (next[at].bytes != nullptr) at = (at + 1) & mask;
            next[at] = stored;
        }
        shard.table.swap(next);
    }
    const std::size_t mask = shard.table.size() - 1;
    std::size_t at = (hash * 0x9e3779b97f4a7c15ULL) & mask;
    for (;;) {
        UnionPrefilterCache::Stored& stored = shard.table[at];
        if (stored.bytes == nullptr) break;
        if (stored.hash == hash && stored.length == length && std::memcmp(stored.bytes, bytes, length) == 0) {
            return stored.id;
        }
        at = (at + 1) & mask;
    }
    constexpr std::size_t kSlab = 1u << 16;   // 256 shards: a large slab each would idle ~256 MB on a small locus
    if (shard.slabs.empty() || shard.slab_used + length > kSlab) {
        shard.slabs.push_back(std::make_unique<char[]>(kSlab));
        shard.slab_used = 0;
        shard.bytes += kSlab;
    }
    char* copy = shard.slabs.back().get() + shard.slab_used;
    std::memcpy(copy, bytes, length);
    shard.slab_used += length;
    const std::uint32_t id = cache.next_window.fetch_add(1, std::memory_order_relaxed);
    if (id >= (std::uint32_t{1} << UnionPrefilterCache::kWindowBits)) {
        throw ResourceExhausted("the union prefilter cache's window index outgrew its 27-bit identities");
    }
    shard.table[at] = {hash, copy, static_cast<std::uint32_t>(length), id};
    ++shard.occupied;
    return id;
}

// Every eligible start's window content ID on one text (eligible: at least read_length bases remain).
void window_content_ids(UnionPrefilterCache& cache, const std::string& text, std::vector<std::uint32_t>& ids) {
    const auto mark = std::chrono::steady_clock::now();
    const std::size_t r = cache.read_length;
    const std::size_t W = cache.window_length;
    const std::size_t N = text.size();
    ids.assign(N >= r ? N - r + 1 : 0, 0);
    constexpr std::uint64_t kBase = 0x100000001b3ULL;
    std::uint64_t base_power = 1;
    for (std::size_t i = 1; i < W; ++i) base_power *= kBase;
    const auto code = [](char c) -> std::uint64_t { return static_cast<unsigned char>(c); };
    std::uint64_t rolling = 0;
    if (N >= W) for (std::size_t i = 0; i < W; ++i) rolling = rolling * kBase + code(text[i]);
    for (std::size_t pos = 0; pos + r <= N; ++pos) {
        const std::size_t length = std::min(N - pos, W);
        std::uint64_t hash = 0;
        if (length == W) {
            if (pos > 0) rolling = (rolling - code(text[pos - 1]) * base_power) * kBase + code(text[pos + W - 1]);
            hash = rolling;
        } else {
            for (std::size_t i = 0; i < length; ++i) hash = hash * kBase + code(text[pos + i]);
        }
        hash ^= static_cast<std::uint64_t>(length) * 0xbf58476d1ce4e5b9ULL;
        ids[pos] = window_content_id(cache, text.data() + pos, length, hash);
    }
    cache.positions_indexed.fetch_add(ids.size(), std::memory_order_relaxed);
    cache.index_nanos.fetch_add(static_cast<std::uint64_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(
                                    std::chrono::steady_clock::now() - mark).count()),
                                std::memory_order_relaxed);
}

std::size_t decision_slot(const UnionPrefilterCache& cache, std::uint64_t key) {
    return static_cast<std::size_t>((key * 0x9e3779b97f4a7c15ULL) >> (64 - cache.slot_bits));
}

// -1 unknown, else the remembered decision.
int prefilter_decision(const UnionPrefilterCache& cache, std::uint64_t key) {
    const std::uint64_t want = key + 1;
    const std::size_t mask = (std::size_t{1} << cache.slot_bits) - 1;
    for (std::size_t at = decision_slot(cache, key);; at = (at + 1) & mask) {
        const std::uint64_t entry = cache.slots[at].load(std::memory_order_acquire);
        if (entry == 0) return -1;
        if ((entry >> 1) == want) return static_cast<int>(entry & 1);
    }
}

void remember_prefilter(UnionPrefilterCache& cache, std::uint64_t key, bool feasible) {
    if (cache.held.load(std::memory_order_relaxed) >= cache.max_held) {
        cache.not_retained.fetch_add(1, std::memory_order_relaxed);
        return;
    }
    const std::uint64_t packed = ((key + 1) << 1) | (feasible ? 1u : 0u);
    const std::size_t mask = (std::size_t{1} << cache.slot_bits) - 1;
    for (std::size_t at = decision_slot(cache, key);; at = (at + 1) & mask) {
        std::uint64_t entry = cache.slots[at].load(std::memory_order_acquire);
        if (entry == 0) {
            if (cache.slots[at].compare_exchange_strong(entry, packed, std::memory_order_acq_rel)) {
                cache.held.fetch_add(1, std::memory_order_relaxed);
                return;
            }
        }
        if ((entry >> 1) == key + 1) return;   // another worker remembered it first: the same decision
    }
}

// The same list by an out-of-place scatter: each proposal's start is written into its form's bucket of `scratch`
// (4 bytes per proposal; the form is the bucket), each bucket sorted and made unique, and the pairs written back in
// order. Exactly std::sort + std::unique's result. `scratch` must already hold proposed.size() entries of
// capacity -- the caller's budget decides that, so nothing here allocates it.
void sort_unique_proposals_scatter(std::vector<std::pair<std::uint32_t, std::uint32_t>>& proposed,
                                   std::size_t form_count, std::vector<std::uint32_t>& cursor,
                                   std::vector<std::uint32_t>& bound, std::vector<std::uint32_t>& scratch) {
    cursor.assign(form_count, 0);
    for (const auto& item : proposed) ++cursor[item.first];
    bound.resize(form_count);
    std::uint32_t running = 0;
    for (std::size_t form = 0; form < form_count; ++form) {
        const std::uint32_t size = cursor[form];
        cursor[form] = running;
        running += size;
        bound[form] = running;
    }
    scratch.resize(proposed.size());   // within the capacity the caller reserved
    for (const auto& item : proposed) scratch[cursor[item.first]++] = item.second;
    std::size_t out = 0;
    std::size_t begin = 0;
    for (std::size_t form = 0; form < form_count; ++form) {
        const std::size_t end = bound[form];
        if (begin == end) continue;
        std::uint32_t* first = scratch.data() + begin;
        std::uint32_t* last = scratch.data() + end;
        if (end - begin <= 24) {
            for (std::uint32_t* it = first + 1; it < last; ++it) {
                const std::uint32_t value = *it;
                std::uint32_t* hole = it;
                while (hole > first && value < *(hole - 1)) {
                    *hole = *(hole - 1);
                    --hole;
                }
                *hole = value;
            }
        } else {
            std::sort(first, last);
        }
        for (std::uint32_t* it = first; it < last; ++it) {
            if (it != first && *it == *(it - 1)) continue;
            proposed[out++] = {static_cast<std::uint32_t>(form), *it};
        }
        begin = end;
    }
    proposed.resize(out);
}

// What the start sorter did over one pass: reset when a pass begins, read when it ends.
struct StartSortTally {
    std::atomic<std::uint64_t> out_of_place{0};
    std::atomic<std::uint64_t> in_place{0};
    std::atomic<std::uint64_t> fallbacks{0};            // scratch was budgeted but too small for this list
    std::atomic<std::uint64_t> peak_worker_bytes{0};    // the largest scratch any worker held
};
StartSortTally& start_sort_tally() {
    static StartSortTally tally;
    return tally;
}

// MEASUREMENT ONLY: where one collect_form_placements call spent its time, in nanoseconds. A handful of clock
// reads per call; filled only when a caller passes it.
struct SeedingTimes {
    std::uint64_t scan_ns = 0;        // piece scan, index lookups and start generation
    std::uint64_t sort_ns = 0;        // start ordering and deduplication
    std::uint64_t verify_ns = 0;      // per start: window checks, prefilter, emission lookup, DP, row storage
    std::uint64_t replay_ns = 0;      // phase-1 cached placements appended
    std::uint64_t row_sort_ns = 0;    // the final (start, consumed) row sort
};

// Scan one haplotype once per piece-length bucket, turning every piece hit into an implied read
// start. This is the inversion the design calls for: the haplotype is walked, not queried once per
// read against a position hash the size of the locus.
template <typename Forms>
void collect_form_placements(
    const std::string& haplotype,
    const Forms& forms,
    // Parallel to `forms`: whether this form's consumed-span classes must be carried separately.
    // False only when the fragment's insert support makes the geometry's floor unreachable, where
    // emitting the whole emission once per anchor is exactly equal and far cheaper.
    const std::vector<char>& componentise,
    const std::vector<PieceIndex>& indexes,
    PlacementBudget& budget,
    std::size_t band,
    const LibraryModel& library,
    std::vector<std::vector<FormPlacement>>& out,
    // Worker-local scratch, reused across haplotypes. This was a function-local vector, so every
    // haplotype paid a fresh allocation and every one of its growth steps again -- churn, not
    // residency, and the kind an allocator holds onto rather than returning.
    std::vector<std::pair<std::uint32_t, std::uint32_t>>& proposed,
    FragmentScoringCounters& counters,
    // Worker-local observational totals. Summed in thread-index order after placement, so the
    // reported figures do not depend on the schedule even though no atomic guards them.
    EmissionCacheCounters& cache_work,
    const FragmentScoringOptions& options,
    // When given, receives the index of every row this call writes to for the first time, so a
    // caller holding persistent rows can clear and sort just those instead of walking all 22,754.
    std::vector<std::uint32_t>* touched = nullptr,
    // Exact union reuse: starts inside these zones are not proposed here; their cached placements are
    // replayed instead, and counted against `replay`'s consumption cap. `skip` is worker scratch.
    const std::vector<ReplayZone>* zones = nullptr,
    std::vector<std::uint8_t>* skip = nullptr,
    UnionReplay* replay = nullptr,
    SeedingTimes* times = nullptr,
    UnionPrefilterCache* prefilter_cache = nullptr,
    // EXPERIMENTAL (I1.4): only starts in [first, second) are verified -- the discovery unit's owned core; the rest
    // are dropped after deduplication and before any verification, and counted in `owned_skipped`.
    const std::pair<std::uint32_t, std::uint32_t>* owned = nullptr,
    std::uint64_t* owned_skipped = nullptr) {
    auto phase_mark = times != nullptr ? std::chrono::steady_clock::now() : std::chrono::steady_clock::time_point{};
    const auto lap = [&](std::uint64_t SeedingTimes::*field) {
        if (times == nullptr) return;
        const auto now = std::chrono::steady_clock::now();
        times->*field += static_cast<std::uint64_t>(
            std::chrono::duration_cast<std::chrono::nanoseconds>(now - phase_mark).count());
        phase_mark = now;
    };
    proposed.clear();   // keeps capacity: the point of hoisting it
    const bool replaying = zones != nullptr && !zones->empty();
    if (replaying) {
        skip->assign(haplotype.size() + library.max_gap_bp + 2, 0);
        for (const ReplayZone& zone : *zones) {
            std::fill(skip->begin() + zone.lo, skip->begin() + zone.hi + 1, static_cast<std::uint8_t>(1));
        }
    }
    // Piece hits accumulated since the last publication to the shared budget. Drained before this
    // function returns, so nothing a worker counted goes unpublished and the merged total stays exact.
    std::uint64_t unpublished_piece_hits = 0;
    const auto publish_piece_hits = [&]() {
        if (unpublished_piece_hits == 0) return;
        const std::uint64_t total =
            budget.piece_hits.fetch_add(unpublished_piece_hits, std::memory_order_relaxed) +
            unpublished_piece_hits;
        unpublished_piece_hits = 0;
        if (total > budget.max_piece_hits) {
            budget.cancel();
            throw ResourceExhausted(
                "fragment scoring reached its piece-hit cap; the locus is too repetitive for the "
                "declared caps and no partial ranking is returned");
        }
    };

    for (const PieceIndex& index : indexes) {
        const std::size_t piece = index.piece_bp;
        if (haplotype.size() < piece) continue;
        counters.candidate_bases_scanned += haplotype.size();

        const std::uint64_t mask = piece == 32
            ? std::numeric_limits<std::uint64_t>::max()
            : ((std::uint64_t{1} << (2 * piece)) - 1);
        std::uint64_t code = 0;
        std::size_t run = 0;
        for (std::size_t at = 0; at < haplotype.size(); ++at) {
            const std::uint64_t bits = two_bit(haplotype[at]);
            if (bits == kNoCode) { run = 0; code = 0; continue; }
            code = ((code << 2) | bits) & mask;
            if (++run < piece) continue;
            const std::size_t piece_start = at + 1 - piece;

            const auto range = std::equal_range(
                index.entries.begin(), index.entries.end(), code,
                [](const auto& lhs, const auto& rhs) {
                    if constexpr (std::is_same_v<std::decay_t<decltype(lhs)>, PieceEntry>) {
                        return lhs.code < rhs;
                    } else {
                        return lhs < rhs.code;
                    }
                });
            for (auto it = range.first; it != range.second; ++it) {
                // The hottest counter in the program -- 320M increments on C4 at audit 256 -- so it is
                // accumulated locally and published to the shared budget in batches, to avoid unmeasured
                // contention on a shared cache line; whether an exact atomic would be affordable here
                // has not been measured. The overshoot batching admits is bounded by
                // threads * publish_batch and declared in the header. Refusing here is sound the moment
                // the published global total exceeds the cap, because a published total can only
                // undercount the exact one; acceptance is still decided later, on the merged totals.
                ++counters.piece_hits;
                if (++unpublished_piece_hits >= budget.publish_batch) {
                    const std::uint64_t total =
                        budget.piece_hits.fetch_add(unpublished_piece_hits,
                                                    std::memory_order_relaxed) +
                        unpublished_piece_hits;
                    unpublished_piece_hits = 0;
                    if (total > budget.max_piece_hits) {
                        budget.cancel();
                        throw ResourceExhausted(
                            "fragment scoring reached its piece-hit cap; the locus is too repetitive "
                            "for the declared caps and no partial ranking is returned");
                    }
                    if (budget.stop()) throw PlacementCancelled();
                }
                if (it->offset > piece_start) continue;   // implied start would be negative
                // The seed implies a start only under a gapless alignment. An indel BEFORE the seed
                // shifts the true start by up to max_gap_bp in either direction, so the whole
                // neighbourhood is proposed. Pigeonhole is unaffected -- d edits still break at most
                // d pieces, so no read the substitution-only scan found is lost -- and the extra
                // starts are deduplicated below before any verification happens.
                const std::size_t implied = piece_start - it->offset;
                const std::size_t lo = implied > library.max_gap_bp ? implied - library.max_gap_bp : 0;
                const std::size_t hi = implied + library.max_gap_bp;
                for (std::size_t shifted = lo; shifted <= hi; ++shifted) {
                    if (replaying && shifted < skip->size() && (*skip)[shifted] != 0) continue;   // served by the cache
                    proposed.emplace_back(it->form, static_cast<std::uint32_t>(shifted));
                    ++counters.proposed_starts;
                }
                if (proposed.size() > options.max_proposed_starts_per_haplotype) {
                    throw ResourceExhausted(
                        "fragment scoring reached its per-haplotype proposed-start cap; no partial "
                        "ranking is returned");
                }
            }
        }
    }

    // Drain whatever this haplotype counted below the last batch boundary, so the shared total is the
    // true total rather than a multiple of publish_batch, and the merged totals stay exact.
    publish_piece_hits();
    lap(&SeedingTimes::scan_ns);

    // Sorting and full-read verification are the rest of this haplotype's cost and neither touches the
    // budget, so without these checks a worker cancelled by someone else's refusal would still pay for
    // both before noticing. Prompt cancellation has to mean prompt in the expensive loops too.
    if (budget.stop()) throw PlacementCancelled();

    // Several pieces of one read can rediscover one physical placement. Deduplicate before
    // verifying, so a placement found three ways is verified once and counted once. Counting by form and
    // sorting each form's few starts gives exactly std::sort + std::unique's list (gated), at a fraction of
    // the cost; the storage is the worker's own.
    {
        // Out of place when the worker's share of the declared scratch budget holds this list (4 bytes a
        // proposal), in place otherwise; either gives exactly the same list. The worker's scratch only grows, to
        // at most its share, and is kept across haplotypes.
        thread_local std::vector<std::uint32_t> form_cursor;
        thread_local std::vector<std::uint32_t> form_bound;
        thread_local std::vector<std::uint32_t> start_scratch;
        const std::size_t share = options.start_sort_scratch_bytes / std::max<std::size_t>(1, options.threads);
        StartSortTally& tally = start_sort_tally();
        const bool fits = !proposed.empty() && share > 0 && proposed.size() * sizeof(std::uint32_t) <= share;
        if (fits) {
            if (start_scratch.capacity() < proposed.size()) {
                std::vector<std::uint32_t>().swap(start_scratch);
                start_scratch.reserve(proposed.size());
                if (start_scratch.capacity() * sizeof(std::uint32_t) > share) {
                    // The allocator rounded past the share: give it back rather than exceed the budget.
                    std::vector<std::uint32_t>().swap(start_scratch);
                }
                const std::uint64_t held = start_scratch.capacity() * sizeof(std::uint32_t);
                std::uint64_t seen = tally.peak_worker_bytes.load(std::memory_order_relaxed);
                while (held > seen && !tally.peak_worker_bytes.compare_exchange_weak(seen, held)) {
                }
            }
        }
        if (fits && start_scratch.capacity() >= proposed.size()) {
            sort_unique_proposals_scatter(proposed, forms.size(), form_cursor, form_bound, start_scratch);
            tally.out_of_place.fetch_add(1, std::memory_order_relaxed);
        } else if (!proposed.empty()) {
            sort_unique_proposals_bucketed(proposed, forms.size(), form_cursor, form_bound);
            tally.in_place.fetch_add(1, std::memory_order_relaxed);
            if (options.start_sort_scratch_bytes > 0) tally.fallbacks.fetch_add(1, std::memory_order_relaxed);
        }
    }
    // Whether this haplotype strand is all ACGT, established once: with an all-ACGT form, no window check
    // can then fail, and gapped_emission is told so.
    const bool haplotype_acgt = !proposed.empty() && std::all_of(haplotype.begin(), haplotype.end(), is_acgt);
    // The prefilter cache serves ACGT texts only (its window identities are taken over whole texts).
    UnionPrefilterCache* const remembered = haplotype_acgt ? prefilter_cache : nullptr;
    thread_local std::vector<std::uint32_t> window_ids;
    if (remembered != nullptr) window_content_ids(*remembered, haplotype, window_ids);
    std::uint64_t unpublished_hits = 0;
    std::uint64_t local_misses = 0;
    const auto publish_hits = [&]() {
        if (remembered == nullptr || unpublished_hits == 0) return;
        const std::uint64_t total =
            remembered->consumed.fetch_add(unpublished_hits, std::memory_order_relaxed) + unpublished_hits;
        unpublished_hits = 0;
        if (total > remembered->max_consumed) {
            budget.cancel();
            throw ResourceExhausted(
                "fragment scoring reached its cached-prefilter consumption cap; no partial ranking is returned");
        }
    };
    if (owned != nullptr) {
        const auto keep_end = std::remove_if(proposed.begin(), proposed.end(), [&](const std::pair<std::uint32_t, std::uint32_t>& item) {
            return item.second < owned->first || item.second >= owned->second;
        });
        if (owned_skipped != nullptr) *owned_skipped += static_cast<std::uint64_t>(proposed.end() - keep_end);
        proposed.erase(keep_end, proposed.end());
    }
    counters.unique_starts += proposed.size();
    if (budget.stop()) throw PlacementCancelled();
    lap(&SeedingTimes::sort_ns);

    // MEASUREMENT ONLY: the shadow piece screen and emission timing. Neither changes a verdict.
    const bool shadow = piece_screen_shadow_active().load(std::memory_order_relaxed);
    const bool timing = shadow || emission_timing_active().load(std::memory_order_relaxed);
    PieceScreenShadowReport shadow_local;
    HaplotypeCodes shadow_codes;
    shadow_codes.text = &haplotype;
    const std::size_t shadow_pieces[3] = {band + 1, band + 2, 2 * band + 1};
    PieceLayout shadow_layouts[3];
    std::uint32_t shadow_form = std::numeric_limits<std::uint32_t>::max();
    for (const auto& item : proposed) {
        // One relaxed load against a full-read comparison: free at this granularity, and verification is
        // where a cancelled worker would otherwise spend the longest.
        if (budget.stop()) throw PlacementCancelled();
        ++counters.full_read_verifications;
        ScreenVerdict verdict[3] = {ScreenVerdict::Admit, ScreenVerdict::Admit, ScreenVerdict::Admit};
        if (shadow) {
            const std::string& read = forms[item.first].bases;
            if (item.first != shadow_form) {
                const auto mark = std::chrono::steady_clock::now();
                for (int i = 0; i < 3; ++i) shadow_layouts[i] = piece_layout(read, shadow_pieces[i]);
                shadow_local.layout_ns += nanos_since(mark);
                shadow_form = item.first;
                for (int i = 0; i < 3; ++i) {
                    for (const std::uint32_t length : shadow_layouts[i].length) (void)shadow_codes.of(length);
                }
            }
            if (item.second <= haplotype.size() && haplotype.size() - item.second >= read.size()) {
                const std::size_t window = std::min(haplotype.size() - item.second, read.size() + library.max_gap_bp);
                for (int i = 0; i < 3; ++i) {
                    const auto mark = std::chrono::steady_clock::now();
                    verdict[i] = piece_screen(shadow_codes, item.second, window, shadow_layouts[i], band);
                    shadow_local.rows[i].screen_ns += nanos_since(mark);
                }
            }
        }
        EmissionWork work_before;
        std::uint64_t prefilter_before = 0;
        if (shadow) {
            work_before = emission_work();
            prefilter_before = emission_timing().prefilter_ns;
        }
        const auto emission_mark = timing ? std::chrono::steady_clock::now() : std::chrono::steady_clock::time_point{};
        // A remembered decision stands in for the kernel on exactly the same (read, window) content.
        int known = -1;
        std::uint64_t decision_key = 0;
        bool consult = false;
        if (remembered != nullptr && forms[item.first].all_acgt && forms[item.first].bases.size() == remembered->read_length &&
            item.second <= haplotype.size() && haplotype.size() - item.second >= remembered->read_length) {
            consult = true;
            decision_key = (static_cast<std::uint64_t>(forms[item.first].content_id) << UnionPrefilterCache::kWindowBits) |
                           window_ids[item.second];
            known = prefilter_decision(*remembered, decision_key);
            if (known >= 0) {
                if (++unpublished_hits >= remembered->publish_batch) publish_hits();
            } else {
                ++local_misses;
            }
        }
        bool computed = false;
        const std::optional<GappedFit> fit =
            gapped_emission_core(haplotype, forms[item.first].bases, item.second, library, band,
                                 &forms[item.first].pattern, haplotype_acgt && forms[item.first].all_acgt,
                                 &forms[item.first].key_prefix, known, consult && known < 0 ? &computed : nullptr,
                                 options.gap_recursion_backend);
        if (consult && known < 0) remember_prefilter(*remembered, decision_key, computed);
        if (timing) {
            const std::uint64_t spent = nanos_since(emission_mark);
            emission_timing().emission_ns += spent;
            if (shadow) {
                const EmissionWork& work_after = emission_work();
                const bool entered = work_after.executed_prefilters != work_before.executed_prefilters;
                const bool passed = entered && work_after.prefilter_rejections == work_before.prefilter_rejections;
                const bool accepted = fit.has_value() && fit->log_emission != kNegInf;
                const std::uint64_t prefilter_spent = emission_timing().prefilter_ns - prefilter_before;
                ++shadow_local.starts;
                shadow_local.emission_ns += spent;
                if (entered) {
                    ++shadow_local.entering;
                    shadow_local.prefilter_ns += prefilter_spent;
                    if (passed) ++shadow_local.prefilter_passed;
                    if (accepted) ++shadow_local.accepted;
                    for (int i = 0; i < 3; ++i) {
                        PieceScreenShadowRow& row = shadow_local.rows[i];
                        ++row.entering;
                        if (verdict[i] != ScreenVerdict::CountReject) ++row.count_pass;
                        if (verdict[i] == ScreenVerdict::Admit) {
                            ++row.admitted;
                            continue;
                        }
                        ++row.rejected;
                        row.rejected_emission_ns += spent;
                        row.rejected_prefilter_ns += prefilter_spent;
                        if (passed) ++row.rejected_prefilter_passed;
                        if (accepted) ++row.rejected_accepted;
                    }
                }
            }
        }
        if (!fit.has_value() || fit->log_emission == kNegInf) continue;
        if (componentise[item.first] != 0) {
            // One row per consumed-span class. The scorer must see the same decomposition the oracle
            // sees, or the two agree only where no alignment is gapped.
            for (const GappedComponent& component : fit->components) {
                if (component.log_emission == kNegInf) continue;
                FormPlacement placement;
                placement.start = item.second;
                placement.consumed = static_cast<std::uint32_t>(component.consumed_bp);
                placement.mismatches = static_cast<std::uint32_t>(component.edits);
                placement.emission = component.log_emission;
                if (touched != nullptr && out[item.first].empty()) {
                    touched->push_back(item.first);
                }
                out[item.first].push_back(placement);
                ++counters.accepted_mate_placements;
            }
            continue;
        }
        // The exact fast path, not an approximation. The insert is the distance between the two
        // ANCHORS, so a consumed span changes the likelihood only through the floor
        // insert >= max(consumed spans). Where the declared support cannot reach that floor, every
        // span pair at an anchor pair is admitted, the sum factorises into the two mates' whole
        // emissions, and carrying the classes separately would compute the same number through
        // (2*max_gap_bp+1)^2 times as many terms. Measured: identical to the last bit above the
        // floor, and genuinely different below it -- which is why the predicate, not a tolerance,
        // decides. The consumed span recorded here is the highest-mass class, used only for the
        // floor check that cannot bind and for diagnostics.
        FormPlacement placement;
        placement.start = item.second;
        placement.consumed = static_cast<std::uint32_t>(fit->consumed_bp);
        placement.mismatches = static_cast<std::uint32_t>(fit->edits);
        placement.emission = fit->log_emission;
        if (touched != nullptr && out[item.first].empty()) touched->push_back(item.first);
        out[item.first].push_back(placement);
        ++counters.accepted_mate_placements;
    }
    if (shadow) {
        shadow_local.code_build_ns += shadow_codes.build_ns;
        const std::lock_guard<std::mutex> lock(piece_screen_mutex());
        PieceScreenShadowReport& totals = piece_screen_totals();
        totals.starts += shadow_local.starts;
        totals.entering += shadow_local.entering;
        totals.prefilter_passed += shadow_local.prefilter_passed;
        totals.accepted += shadow_local.accepted;
        totals.emission_ns += shadow_local.emission_ns;
        totals.prefilter_ns += shadow_local.prefilter_ns;
        totals.code_build_ns += shadow_local.code_build_ns;
        totals.layout_ns += shadow_local.layout_ns;
        for (int i = 0; i < 3; ++i) {
            PieceScreenShadowRow& to = totals.rows[i];
            const PieceScreenShadowRow& from = shadow_local.rows[i];
            to.pieces = shadow_pieces[i];
            to.entering += from.entering;
            to.count_pass += from.count_pass;
            to.admitted += from.admitted;
            to.rejected += from.rejected;
            to.rejected_prefilter_passed += from.rejected_prefilter_passed;
            to.rejected_accepted += from.rejected_accepted;
            to.screen_ns += from.screen_ns;
            to.rejected_emission_ns += from.rejected_emission_ns;
            to.rejected_prefilter_ns += from.rejected_prefilter_ns;
        }
        totals.band = band;
    }
    // What the emission actually did on this thread since it was last asked. Taken here rather than
    // threaded through gapped_emission's signature, so the model's entry point keeps its shape.
    const GappedEmissionCounters emission = take_gapped_emission_counters();
    counters.dp_evaluations += emission.dp_evaluations;
    counters.prefilter_rejections += emission.prefilter_rejections;
    cache_work.executed_prefilters += emission.executed_prefilters;
    cache_work.dp_cache_requests += emission.dp_cache_requests;
    cache_work.dp_cache_hits += emission.dp_cache_hits;
    cache_work.dp_cache_misses += emission.dp_cache_misses;
    cache_work.executed_dps += emission.executed_dps;
    publish_hits();
    if (remembered != nullptr) remembered->misses.fetch_add(local_misses, std::memory_order_relaxed);
    lap(&SeedingTimes::verify_ns);

    // Exact union reuse: the cached placements of every zone, translated to this haplotype. Appended before
    // the sort below, which puts every row in (start, consumed) order -- unique within a row -- so the rows
    // are the reference's exactly. Each occurrence is replayed, so repeat multiplicity is preserved.
    if (replaying) {
        std::uint64_t replayed = 0;
        for (const ReplayZone& zone : *zones) {
            for (const CachedPlacement& cached : zone.entry->placements) {
                FormPlacement placement;
                placement.start = cached.start + zone.shift;
                placement.consumed = cached.consumed;
                placement.mismatches = cached.mismatches;
                placement.emission = cached.emission;
                out[cached.form].push_back(placement);
                ++replayed;
            }
        }
        counters.accepted_mate_placements += replayed;
        const std::uint64_t total = replay->consumed.fetch_add(replayed, std::memory_order_relaxed) + replayed;
        if (replay->max_consumed != 0 && total > replay->max_consumed) {
            budget.cancel();
            throw ResourceExhausted(
                "fragment scoring reached its cached-placement consumption cap; no partial ranking is returned");
        }
    }
    lap(&SeedingTimes::replay_ns);

    // Sorted by start first, so the pairing below can walk one left endpoint's classes as a group,
    // and by consumed span within it so the row order is fixed rather than insertion-dependent.
    // Sorting only the touched rows is identical: an untouched row is empty, and sorting an empty
    // range is a no-op.
    const auto sort_row = [](std::vector<FormPlacement>& list) {
        std::sort(list.begin(), list.end(), [](const FormPlacement& lhs, const FormPlacement& rhs) {
            return std::tie(lhs.start, lhs.consumed) < std::tie(rhs.start, rhs.consumed);
        });
    };
    if (touched != nullptr) {
        for (const std::uint32_t row : *touched) sort_row(out[row]);
    } else {
        for (std::vector<FormPlacement>& list : out) sort_row(list);
    }
    lap(&SeedingTimes::row_sort_ns);
}

} // namespace scoring_detail

using namespace scoring_detail;

PieceHitCensus census_piece_hits(
    const InternedHaplotypes& haplotypes,
    const std::vector<ReadPair>& fragments,
    const FragmentScoringOptions& options) {
    PieceHitCensus census;
    const std::size_t band = options.mismatch_band;

    // The same four oriented forms per fragment the scorer builds: two forward-role mates scanned
    // against the haplotype, and two reverse-role mates stored reversed -- complement(mate) -- and
    // scanned against the reversed haplotype. Both are indexed together here because a piece hit
    // costs the same whichever direction its form is scanned in, and the census exists to count
    // hits; the haplotype is scanned in both orientations below so the total matches the scorer's.
    std::vector<ReadForm> forms;
    forms.reserve(fragments.size() * 4);
    for (std::size_t f = 0; f < fragments.size(); ++f) {
        validate_read_pair(fragments[f]);
        // The census counts piece hits and never runs a feasibility test, so the pattern is left
        // default-constructed -- stated explicitly rather than left to a missing-initialiser
        // warning, since an unset field that happens not to be read is exactly the kind of gap that
        // becomes a defect the moment someone reads it.
        forms.push_back({f, fragments[f].first.bases, true, MyersPattern{}});
        forms.push_back({f, complement_read(fragments[f].second.bases), false, MyersPattern{}});
        forms.push_back({f, fragments[f].second.bases, false, MyersPattern{}});
        forms.push_back({f, complement_read(fragments[f].first.bases), true, MyersPattern{}});
    }
    census.total_forms = forms.size();
    if (forms.empty()) return census;

    // A collapse keeps exactly the first form of each distinct sequence; every other form would be
    // answered from it. Marking them is what makes the prediction exact rather than proportional.
    std::set<std::string> seen;
    std::vector<char> first_of_sequence(forms.size(), 0);
    for (std::size_t i = 0; i < forms.size(); ++i) {
        if (seen.insert(forms[i].bases).second) first_of_sequence[i] = 1;
    }
    census.unique_forms = seen.size();

    const std::vector<PieceIndex> indexes = build_piece_indexes(forms, band, options.minimum_piece_bp);
    for (const PieceIndex& index : indexes) {
        census.total_entries += index.entries.size();
        for (const PieceEntry& entry : index.entries) {
            census.unique_entries += first_of_sequence[entry.form];
        }
    }

    for (const std::string& forward_haplotype : haplotypes.sequences) {
      const std::string reversed_haplotype(forward_haplotype.rbegin(), forward_haplotype.rend());
      std::uint64_t hits_here = 0;
      for (int direction = 0; direction < 2; ++direction) {
        const std::string& haplotype = direction == 0 ? forward_haplotype : reversed_haplotype;
        for (const PieceIndex& index : indexes) {
            const std::size_t piece = index.piece_bp;
            if (piece == 0 || haplotype.size() < piece) continue;
            const std::uint64_t mask = piece >= 32 ? ~0ULL : ((1ULL << (2 * piece)) - 1);
            std::uint64_t code = 0;
            std::size_t valid = 0;
            for (std::size_t at = 0; at < haplotype.size(); ++at) {
                const std::uint64_t bits = two_bit(haplotype[at]);
                if (bits == kNoCode) {
                    valid = 0;
                    code = 0;
                    continue;
                }
                code = ((code << 2) | bits) & mask;
                if (++valid < piece) continue;
                const auto range = std::equal_range(
                    index.entries.begin(), index.entries.end(), code,
                    [](const auto& lhs, const auto& rhs) {
                        if constexpr (std::is_same_v<std::decay_t<decltype(lhs)>, PieceEntry>) {
                            return lhs.code < rhs;
                        } else {
                            return lhs < rhs.code;
                        }
                    });
                std::uint64_t all = 0;
                std::uint64_t unique = 0;
                for (auto it = range.first; it != range.second; ++it) {
                    ++all;
                    unique += first_of_sequence[it->form];
                }
                census.predicted_piece_hits_all += all;
                census.predicted_piece_hits_unique += unique;
                hits_here += all;
                census.max_piece_hits_at_one_seed = std::max(census.max_piece_hits_at_one_seed, all);
            }
        }
      }
      census.max_piece_hits_one_haplotype =
          std::max(census.max_piece_hits_one_haplotype, hits_here);
    }
    return census;
}

double FragmentScoringResult::mass_of(std::size_t fragment, std::size_t haplotype) const {
    if (unique_haplotypes == 0 || haplotype >= unique_haplotypes ||
        fragment * unique_haplotypes + haplotype >= log_mass.size()) {
        throw std::out_of_range("fragment/haplotype mass index is out of range");
    }
    return log_mass[fragment * unique_haplotypes + haplotype];
}

InternedHaplotypes intern_haplotypes(
    const std::vector<MosaicPath>& allele_paths,
    const std::vector<std::string>& sequences) {
    if (allele_paths.size() != sequences.size()) {
        throw std::invalid_argument("allele paths and spelled sequences differ in size");
    }
    InternedHaplotypes result;
    result.index_of.reserve(sequences.size());
    std::map<std::string, std::size_t> slot_by_sequence;
    std::set<std::vector<std::size_t>> distinct_paths;
    for (std::size_t i = 0; i < sequences.size(); ++i) {
        distinct_paths.insert(allele_paths[i].alleles);
        const auto found = slot_by_sequence.find(sequences[i]);
        if (found != slot_by_sequence.end()) {
            result.index_of.push_back(found->second);
            continue;
        }
        const std::size_t slot = result.sequences.size();
        slot_by_sequence.emplace(sequences[i], slot);
        result.sequences.push_back(sequences[i]);
        result.index_of.push_back(slot);
    }
    // Two different allele paths can spell identical bytes. They are one haplotype to any read-based
    // method, so they are scored once -- but the count of distinct paths is retained, because
    // collapsing them silently would hide that the candidate set is narrower than it looks.
    result.distinct_allele_paths = distinct_paths.size();
    return result;
}

namespace scoring_detail {

// ---- fragment preparation, owned once ----------------------------------------------------------
//
// Phase one of the one-shot scorer, lifted out whole so the persistent context prepares fragments by
// the same definition rather than a parallel copy of it. Nothing here depends on a haplotype, which
// is exactly why it can be hoisted out of a per-round call.
struct PreparedFragmentSet {
    std::deque<ReadForm> forward_forms;
    std::deque<ReadForm> reverse_forms;
    std::vector<char> forward_componentise;
    std::vector<char> reverse_componentise;
    std::vector<PieceIndex> forward_indexes;
    std::vector<PieceIndex> reverse_indexes;
    std::vector<InsertSupport> support_by_fragment;
    std::vector<double> background_by_fragment;
    std::vector<char> floor_binds;
    std::map<std::size_t, std::size_t> fragments_by_floor;
    std::uint64_t fragments_with_floor_binding = 0;
};

PreparedFragmentSet prepare_fragment_set(
    const std::vector<ReadPair>& fragments,
    const LibraryModel& library,
    const FragmentScoringOptions& options) {
    PreparedFragmentSet out;
    const std::size_t band = options.mismatch_band;
    out.support_by_fragment.reserve(fragments.size());
    out.background_by_fragment.reserve(fragments.size());
    out.floor_binds.reserve(fragments.size());
    for (std::size_t f = 0; f < fragments.size(); ++f) {
        validate_read_pair(fragments[f]);
        const auto add_form = [](std::deque<ReadForm>& into, std::size_t fragment, std::string bases,
                                 bool first_arm) {
            ReadForm form;
            form.fragment = fragment;
            form.bases = std::move(bases);
            form.is_forward_mate_of_first_arm = first_arm;
            form.pattern = prepare_myers_pattern(form.bases);
            form.all_acgt = std::all_of(form.bases.begin(), form.bases.end(), is_acgt);
            form.key_prefix = key_prefix_of(form.bases);
            into.push_back(std::move(form));
        };
        add_form(out.forward_forms, f, fragments[f].first.bases, true);
        add_form(out.reverse_forms, f, complement_read(fragments[f].second.bases), false);
        add_form(out.forward_forms, f, fragments[f].second.bases, false);
        add_form(out.reverse_forms, f, complement_read(fragments[f].first.bases), true);
        out.support_by_fragment.push_back(fragment_insert_support(fragments[f], library));
        {
            const std::size_t widest_span =
                std::max(fragments[f].first.bases.size(), fragments[f].second.bases.size()) +
                library.max_gap_bp;
            const bool binds = out.support_by_fragment.back().empty() ||
                               out.support_by_fragment.back().min_bp < widest_span;
            out.floor_binds.push_back(binds ? 1 : 0);
            if (binds) ++out.fragments_with_floor_binding;
        }
        out.background_by_fragment.push_back(log_background_probability(fragments[f]));
        ++out.fragments_by_floor[std::max(fragments[f].first.bases.size(),
                                          fragments[f].second.bases.size())];
    }
    out.forward_componentise.assign(out.forward_forms.size(), 0);
    out.reverse_componentise.assign(out.reverse_forms.size(), 0);
    for (std::size_t i = 0; i < out.forward_forms.size(); ++i) {
        out.forward_componentise[i] = out.floor_binds[out.forward_forms[i].fragment];
    }
    for (std::size_t i = 0; i < out.reverse_forms.size(); ++i) {
        out.reverse_componentise[i] = out.floor_binds[out.reverse_forms[i].fragment];
    }
    {
        std::unordered_map<std::string_view, std::uint32_t> content_of;
        for (std::deque<ReadForm>* forms : {&out.forward_forms, &out.reverse_forms}) {
            for (ReadForm& form : *forms) {
                form.content_id = content_of.emplace(std::string_view(form.bases),
                                                     static_cast<std::uint32_t>(content_of.size())).first->second;
            }
        }
    }
    out.forward_indexes = out.forward_forms.empty()
        ? std::vector<PieceIndex>{}
        : build_piece_indexes(out.forward_forms, band, options.minimum_piece_bp);
    out.reverse_indexes = out.reverse_forms.empty()
        ? std::vector<PieceIndex>{}
        : build_piece_indexes(out.reverse_forms, band, options.minimum_piece_bp);
    return out;
}

// One haplotype's exposure, averaged over the observed distribution of fragment physical floors.
// Summed in the map's own order, so it does not matter when a haplotype is scored: the same length
// gives the same value in batch one and batch seven.
double exposure_for_length(
    std::size_t length,
    std::size_t fragment_count,
    const std::map<std::size_t, std::size_t>& fragments_by_floor,
    const DiscreteInsertDistribution& library_insert) {
    if (fragment_count == 0) return haplotype_exposure(length, library_insert, 0);
    double total = 0.0;
    for (const auto& item : fragments_by_floor) {
        total += static_cast<double>(item.second) / static_cast<double>(fragment_count) *
                 haplotype_exposure(length, library_insert, item.first);
    }
    return total;
}

struct UnionOriginReuse;

// Everything placement reads and never writes.
struct PlacementInputs {
    const PreparedFragmentSet* prepared = nullptr;
    const std::vector<ReadPair>* fragments = nullptr;
    const LibraryModel* library = nullptr;
    const DiscreteInsertDistribution* library_insert = nullptr;
    const FragmentScoringOptions* options = nullptr;
    const std::deque<std::string>* sequences = nullptr;
    const std::deque<std::string>* reversed = nullptr;
    std::size_t band = 0;
    UnionReplay* replay = nullptr;   // exact union reuse, one-shot passes only
    UnionOriginReuse* origin_reuse = nullptr;   // exact union origin reuse (phase 2), one-shot passes only
    UnionPrefilterCache* prefilter_cache = nullptr;   // exact union prefilter reuse, one-shot passes only
    // The slots to place, when they are not one contiguous range (released slots reused): the placement range
    // then indexes this list. Null: the range is the slots themselves.
    const std::vector<std::size_t>* slot_list = nullptr;
};

// Where a haplotype's complete origin census is written, when one is wanted. Distinct from
// OriginRecording below: that is a diagnostic carrying alignment detail no score reads, this is the
// exact input a later block-delta reduces from.
struct CensusCapture {
    // [haplotype slot] -> the census to fill, or null for haplotypes nobody will delta from.
    std::map<std::size_t, HaploidOriginCensus*> wanted;
};

// The ONE reduction from a StateKey census to a fragment's mass. Shared by the full scorer and the
// delta so the two cannot reduce in different orders: std::map orders by (insert, orientation,
// latent mass), and log_add_exp is not associative in floating point, so a second copy of this loop
// is how bit-identical quietly becomes merely close.
double reduce_states(const std::map<StateKey, std::uint64_t>& states,
                     const DiscreteInsertDistribution& library_insert) {
    double total = kNegInf;
    for (const auto& item : states) {
        const double term = std::log(static_cast<double>(item.second)) + kLogMateOrientationPrior +
                            library_insert.log_probability(item.first.insert_bp) +
                            item.first.log_latent_mass;
        total = log_add_exp(total, term);
    }
    return total;
}

// Everything one window enumeration needs, owned across calls. The reference built a vector per read
// form on every candidate -- about 57,000 container constructions each, 36 million a round -- which is
// what this exists to stop. Rows persist; only the ones the previous window touched are cleared.
// ---- I1.4 content units (shared with the adaptive window path, I2-B) -------------------------------------------
struct IndexUnit {
    std::uint8_t strand = 0;
    bool interior = false;
    bool at_start = false, at_end = false;   // windows: the clip reached the text's start / end
    std::uint32_t boundary = 0;             // windows: b's offset in the window
    // I1.4 fine discovery unit: owns the starts in [core_lo, core_hi) of its text whose dependency span lies in it
    bool fine = false;
    std::uint32_t core_lo = 0, core_hi = 0;
    std::string bytes;
};
struct IndexOccurrence {
    std::uint32_t unit = 0;
    std::uint32_t offset = 0;   // the unit's first byte, in its strand's own coordinates
};
using StoredPlacement = CachedPlacement;   // (form, start within the unit, consumed, mismatches, emission)
std::atomic<bool>& content_index_hash_collapse() {
    static std::atomic<bool> collapse{false};
    return collapse;
}
std::uint64_t unit_hash(const IndexUnit& u) {
    if (content_index_hash_collapse().load(std::memory_order_relaxed)) return 0;
    std::uint64_t h = 1469598103934665603ULL;
    for (const char c : u.bytes) h = (h ^ static_cast<unsigned char>(c)) * 1099511628211ULL;
    const std::uint64_t tag = (static_cast<std::uint64_t>(u.strand) << 1) ^ (static_cast<std::uint64_t>(u.interior) << 2) ^
                              (static_cast<std::uint64_t>(u.at_start) << 3) ^ (static_cast<std::uint64_t>(u.at_end) << 4) ^
                              (static_cast<std::uint64_t>(u.boundary) << 8) ^ (static_cast<std::uint64_t>(u.fine) << 5) ^
                              (static_cast<std::uint64_t>(u.core_lo) << 24) ^ (static_cast<std::uint64_t>(u.core_hi) << 44);
    return (h ^ tag) * 1099511628211ULL;
}
bool same_unit(const IndexUnit& a, const IndexUnit& b) {
    return a.strand == b.strand && a.interior == b.interior && a.at_start == b.at_start && a.at_end == b.at_end &&
           a.boundary == b.boundary && a.fine == b.fine && a.core_lo == b.core_lo && a.core_hi == b.core_hi && a.bytes == b.bytes;
}
// I1.4: content-defined cores of one interior segment's bytes (start positions [0, len) partitioned). A gear rolling hash
// (FastCDC-style: h = (h << 1) + G[byte], so h depends on the last 64 bytes) cuts after position p when the core holds
// at least kMinCore starts and bits 40..48 of h are zero, or at kMaxCore. Fixed before any read-based measurement. The
// hash's history is finite (64 bytes), but whether a position may cut also depends on the distance to the PRECEDING cut
// (the minimum and maximum core sizes), so an edit can move a cut and, through it, later cuts farther away than its hash
// window: locality is typical, not guaranteed. Exactness never depends on where the cuts fall.
constexpr std::size_t kMinCore = 256, kMaxCore = 2048;
constexpr std::uint64_t kCutMask = 511;
const std::array<std::uint64_t, 256>& gear_table() {
    static const std::array<std::uint64_t, 256> table = [] {
        std::array<std::uint64_t, 256> t{};
        std::uint64_t x = 0x9e3779b97f4a7c15ULL;
        for (auto& v : t) {   // splitmix64
            x += 0x9e3779b97f4a7c15ULL;
            std::uint64_t z = x;
            z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
            z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
            v = z ^ (z >> 31);
        }
        return t;
    }();
    return table;
}
std::vector<std::pair<std::size_t, std::size_t>> content_cores(const std::string& bytes) {
    const auto& G = gear_table();
    std::vector<std::pair<std::size_t, std::size_t>> cores;
    std::uint64_t h = 0;
    std::size_t last = 0;
    for (std::size_t p = 0; p < bytes.size(); ++p) {
        h = (h << 1) + G[static_cast<unsigned char>(bytes[p])];
        const std::size_t size = p + 1 - last;
        if ((size >= kMinCore && ((h >> 40) & kCutMask) == 0) || size >= kMaxCore) {
            cores.push_back({last, p + 1});
            last = p + 1;
        }
    }
    if (last < bytes.size()) cores.push_back({last, bytes.size()});
    return cores;
}


// ---- shared discovery helpers: the I1.4 batch and the I2-B adaptive window path use these same functions ----------
// The fine units of one interior text: each content-defined core with g bytes before it and max_read + g after it
// (clipped to the text), so every owned start's dependency span [s - g, s + L + g) that lies in the text lies in the
// fine text. `emit(unit, offset of its text)`.
template <typename Emit>
void for_each_fine_piece(const std::string& bytes, std::uint8_t strand, std::size_t gap, std::size_t max_read, Emit&& emit) {
    for (const auto& [lo, hi] : content_cores(bytes)) {
        const std::size_t t0 = lo >= gap ? lo - gap : 0;
        const std::size_t t1 = std::min(bytes.size(), hi + max_read + gap);
        IndexUnit f;
        f.strand = strand;
        f.fine = true;
        f.core_lo = static_cast<std::uint32_t>(lo - t0);
        f.core_hi = static_cast<std::uint32_t>(hi - t0);
        f.bytes = bytes.substr(t0, t1 - t0);
        emit(std::move(f), static_cast<std::uint32_t>(t0));
    }
}
// The boundary window of `text` at b: [b - window, b + window) clipped to the text; it keeps the starts whose dependency
// span crosses b. `lo` receives the window's first byte.
IndexUnit boundary_window_unit(const std::string& text, std::uint8_t strand, std::size_t b, std::size_t window, std::size_t& lo) {
    const std::size_t n = text.size();
    lo = b > window ? b - window : 0;
    const std::size_t hi = std::min(n, b + window);
    IndexUnit u;
    u.strand = strand;
    u.at_start = lo == 0;
    u.at_end = hi == n;
    u.boundary = static_cast<std::uint32_t>(b - lo);
    u.bytes = text.substr(lo, hi - lo);
    return u;
}
// Whether a unit keeps a verified start x of a read form of length L: interior, all of D(x) = [x - g, x + L + g) inside
// the text; window, D(x) crosses b; fine (I1.4), x in the owned core and all of D(x) inside the fine text.
bool keep_unit_start(const IndexUnit& unit, std::int64_t x, std::int64_t L, std::int64_t g) {
    const std::int64_t len = static_cast<std::int64_t>(unit.bytes.size());
    const std::int64_t b = unit.boundary;
    return unit.fine ? (x >= unit.core_lo && x < unit.core_hi && x - g >= 0 && x + L + g <= len)
           : unit.interior ? (x - g >= 0 && x + L + g <= len)
                           : (x - g < b && b < x + L + g);
}
// One unit's placements: proposed and verified on the unit's own bytes (a fine unit verifies only the starts of its
// owned core), kept by keep_unit_start, and appended to `store` in (form, start, consumed) order. `rows` ([2], one per
// strand), `proposed` and `touched` are the caller's reusable scratch. Returns the placements accepted before keeping.
std::uint64_t discover_unit_placements(const IndexUnit& unit, const PreparedFragmentSet& prepared,
                                       const std::vector<PieceIndex>& forward_indexes,
                                       const std::vector<PieceIndex>& reverse_indexes, PlacementBudget& budget,
                                       std::size_t band, const LibraryModel& library, const FragmentScoringOptions& options,
                                       std::vector<std::vector<FormPlacement>>* rows,
                                       std::vector<std::pair<std::uint32_t, std::uint32_t>>& proposed,
                                       std::vector<std::uint32_t>& touched, FragmentScoringCounters& counters,
                                       EmissionCacheCounters& cache, SeedingTimes* times, UnionPrefilterCache* prefilter,
                                       std::uint64_t* skipped, std::vector<StoredPlacement>& store) {
    const bool forward = unit.strand == 0;
    auto& out = rows[unit.strand];
    if (out.size() < (forward ? prepared.forward_forms.size() : prepared.reverse_forms.size())) {
        out.resize(forward ? prepared.forward_forms.size() : prepared.reverse_forms.size());
    }
    touched.clear();
    const std::pair<std::uint32_t, std::uint32_t> core{unit.core_lo, unit.core_hi};
    const std::pair<std::uint32_t, std::uint32_t>* owned = unit.fine ? &core : nullptr;
    if (forward) {
        collect_form_placements(unit.bytes, prepared.forward_forms, prepared.forward_componentise, forward_indexes, budget,
                                band, library, out, proposed, counters, cache, options, &touched, nullptr, nullptr, nullptr,
                                times, prefilter, owned, skipped);
    } else {
        collect_form_placements(unit.bytes, prepared.reverse_forms, prepared.reverse_componentise, reverse_indexes, budget,
                                band, library, out, proposed, counters, cache, options, &touched, nullptr, nullptr, nullptr,
                                times, prefilter, owned, skipped);
    }
    std::sort(touched.begin(), touched.end());
    const std::int64_t g = static_cast<std::int64_t>(library.max_gap_bp);
    std::uint64_t accepted = 0;
    for (const std::uint32_t form : touched) {
        const std::int64_t L = static_cast<std::int64_t>(
            (forward ? prepared.forward_forms[form] : prepared.reverse_forms[form]).bases.size());
        for (const FormPlacement& p : out[form]) {
            ++accepted;
            if (keep_unit_start(unit, static_cast<std::int64_t>(p.start), L, g)) {
                store.push_back({form, p.start, p.consumed, p.mismatches, p.emission});
            }
        }
        out[form].clear();
    }
    return accepted;
}

// ---- I2-B (EXPERIMENTAL): round-local fine discovery for the adaptive search's window enumeration ----------------
// A window's raw placement rows -- what collect_form_placements returns on the window's bytes and on its reversal --
// are the union, over the window treated as one text with one segment, of its two end windows' and its fine units'
// kept placements (the I1/I1.4 dependency argument: a start depends only on D(s) = [s - g, s + L + g), clipped where it
// meets an end, and the emission normaliser saturates within it). Units are content keys compared byte for byte, so a
// unit discovered for one window serves every later window of the round that contains it. The store belongs to one
// enumeration worker (its fragment range: placements for its own read forms only), lives for one round
// (add_repaired_haplotypes) and is released at the round's end, refusal or not. Its budget counts allocated capacity
// and metadata; a unit that does not fit is used once and dropped, and closes the store, after which a window not
// fully served is enumerated by the reference path (recorded). Served placements are charged to a cumulative cap.
struct WindowFineStore {
    bool active = false;
    std::size_t max_read = 0, gap = 0, window = 0;
    std::uint64_t budget = 0, served_cap = 0;
    std::atomic<std::uint64_t>* served_total = nullptr;   // shared by every worker's store, cumulative over the search
    std::deque<IndexUnit> units;
    std::deque<std::vector<StoredPlacement>> placements;
    std::unordered_multimap<std::uint64_t, std::uint32_t> by_hash;
    std::uint64_t bytes = 0;     // allocated capacity of the held units, their rows and their index entries
    bool closed = false;
    // round-local scratch (its capacity is in the memory account, not the budget)
    std::vector<std::vector<FormPlacement>> rows[2];
    std::vector<std::uint32_t> touched;
    std::vector<std::pair<std::uint32_t, std::uint32_t>> proposed;
    std::deque<std::vector<StoredPlacement>> transient;   // units used once (not retained)
    struct Piece {
        const std::vector<StoredPlacement>* rows = nullptr;
        std::uint32_t offset = 0;
        std::uint8_t strand = 0;
        bool skip_left_crossing = false;
    };
    std::vector<Piece> pieces;
    // tallies since the last fold (see PreparedScoringContext::State::end_fine_round)
    std::uint64_t windows = 0, fallback_windows = 0, units_discovered = 0, units_served = 0, units_not_retained = 0;
    std::uint64_t placements_served = 0, placements_discovered = 0, hash_collisions = 0, peak_bytes = 0, peak_held = 0;
    std::uint64_t peak_temporary = 0;   // I2-C: scratch + the window's candidate unit copies, at their largest

    static std::uint64_t entry_bytes(const IndexUnit& unit, const std::vector<StoredPlacement>& rows) {
        return sizeof(IndexUnit) + unit.bytes.capacity() + sizeof(std::vector<StoredPlacement>) +
               rows.capacity() * sizeof(StoredPlacement) + 48;   // 48: one hash node, generously
    }
    std::uint64_t scratch_bytes() const {
        std::uint64_t b = touched.capacity() * sizeof(std::uint32_t) + proposed.capacity() * 8 + pieces.capacity() * sizeof(Piece);
        for (const auto& r : rows) {
            b += r.capacity() * sizeof(std::vector<FormPlacement>);
            for (const auto& v : r) b += v.capacity() * sizeof(FormPlacement);
        }
        for (const auto& t : transient) b += t.capacity() * sizeof(StoredPlacement) + sizeof(t);
        return b;
    }
    std::uint64_t held_bytes() const { return bytes + by_hash.bucket_count() * sizeof(void*); }
    std::uint64_t account_bytes() const { return held_bytes() + scratch_bytes(); }
    void release() {
        std::deque<IndexUnit>().swap(units);
        std::deque<std::vector<StoredPlacement>>().swap(placements);
        std::unordered_multimap<std::uint64_t, std::uint32_t>().swap(by_hash);
        for (auto& r : rows) std::vector<std::vector<FormPlacement>>().swap(r);
        std::vector<std::uint32_t>().swap(touched);
        std::vector<std::pair<std::uint32_t, std::uint32_t>>().swap(proposed);
        std::deque<std::vector<StoredPlacement>>().swap(transient);
        std::vector<Piece>().swap(pieces);
        bytes = 0;
        closed = false;
        active = false;
    }
};

struct WindowScratch {
    std::vector<std::vector<FormPlacement>> forward_rows;
    std::vector<std::vector<FormPlacement>> reverse_rows;
    std::vector<std::uint32_t> touched_forward;
    std::vector<std::uint32_t> touched_reverse;
    std::vector<std::pair<std::uint32_t, std::uint32_t>> proposed;
    // Contiguous equal-endpoint groups of the reverse range, replacing the per-forward-start map.
    std::vector<std::size_t> group_begin;
    std::vector<std::size_t> group_end;
    std::vector<OriginAccumulator> accumulators;
    std::vector<char> accepted;
    // Flat output: fragment f owns [offsets[f], offsets[f + 1]).
    std::vector<HaploidOrigin> origins;
    std::vector<std::size_t> offsets;
    // Growth events, which is where the allocations are.
    std::uint64_t row_constructions = 0;
    std::uint64_t accumulator_growths = 0;
    std::uint64_t origin_growths = 0;
    // I2-B: the worker's round-local fine store (inactive unless the context enabled it for the round), and the time this
    // scratch's windows spent building raw rows (discovery: the reference's two scans, or fine reconstruction) and pairing
    // them, worker elapsed, both paths. Measurement only.
    WindowFineStore fine;
    double seconds_rows = 0.0, seconds_pairing = 0.0;
    std::uint64_t windows_enumerated = 0;

    void ensure_rows(std::size_t forward_count, std::size_t reverse_count) {
        if (forward_rows.size() < forward_count) {
            row_constructions += forward_count - forward_rows.size();
            forward_rows.resize(forward_count);
        }
        if (reverse_rows.size() < reverse_count) {
            row_constructions += reverse_count - reverse_rows.size();
            reverse_rows.resize(reverse_count);
        }
    }
    // Only what the previous window wrote.
    void clear_touched() {
        for (const std::uint32_t row : touched_forward) forward_rows[row].clear();
        for (const std::uint32_t row : touched_reverse) reverse_rows[row].clear();
        touched_forward.clear();
        touched_reverse.clear();
    }
};

std::atomic<bool>& window_self_check() {
    static std::atomic<bool> enabled{false};
    return enabled;
}

// ---- stage 6: the window-origin cache and the preaggregated merge ------------------------------

struct StateCount {
    StateKey key;
    std::uint64_t count = 0;
};

// The map's own notion of "same key": neither orders before the other.
inline bool same_state(const StateKey& lhs, const StateKey& rhs) {
    return !(lhs < rhs) && !(rhs < lhs);
}

// One retained origin, relative to its window, packed to 16 bytes. insert and arm share a word as
// (insert << 1 | arm), which orders exactly as StateKey orders (insert, orientation).
struct CachedOrigin {
    std::uint32_t rel_start = 0;
    std::uint32_t insert_and_arm = 0;
    double log_latent_mass = 0.0;
    std::uint32_t insert() const { return insert_and_arm >> 1; }
    std::uint32_t arm() const { return insert_and_arm & 1u; }
};
static_assert(sizeof(CachedOrigin) == 16, "a cached origin is declared at 16 bytes");

// Records per fragment, sparse: only fragments with at least one record appear. Within a fragment,
// records are sorted by StateKey and then by relative start -- so a filtered subsequence is still
// StateKey-sorted, and its multiplicities are run lengths. That is the "sorted StateKey ->
// multiplicity contribution" and the raw relative origins in one array rather than two.
struct OriginRecordSet {
    std::vector<std::uint32_t> fragments;
    std::vector<std::uint32_t> starts;   // fragments.size() + 1
    std::vector<CachedOrigin> origins;
    void clear() {
        fragments.clear();
        starts.clear();
        origins.clear();
    }
    bool operator==(const OriginRecordSet& other) const {
        if (fragments != other.fragments || starts != other.starts) return false;
        if (origins.size() != other.origins.size()) return false;
        for (std::size_t i = 0; i < origins.size(); ++i) {
            const CachedOrigin& a = origins[i];
            const CachedOrigin& b = other.origins[i];
            if (a.rel_start != b.rel_start || a.insert_and_arm != b.insert_and_arm ||
                !(a.log_latent_mass == b.log_latent_mass)) {
                return false;
            }
        }
        return true;
    }
};

// Sparse per-fragment sorted StateKey -> multiplicity arrays, for the merge.
struct StateAggregateSet {
    std::vector<std::uint32_t> fragments;
    std::vector<std::uint32_t> starts;   // fragments.size() + 1 once closed
    std::vector<StateCount> states;
    void clear() {
        fragments.clear();
        starts.clear();
        states.clear();
    }
    std::uint64_t bytes() const {
        return static_cast<std::uint64_t>(fragments.capacity() + starts.capacity()) *
                   sizeof(std::uint32_t) +
               static_cast<std::uint64_t>(states.capacity()) * sizeof(StateCount);
    }
};

// Identity: exact window bytes, the two end flags, and the relative replacement interval that
// decided which origins were retained. window.begin is deliberately NOT part of it -- records are
// relative and translated on use, so the same bytes at another locus offset reuse correctly.
struct OriginCacheEntry {
    std::string bytes;
    bool at_start = false;
    bool at_end = false;
    std::uint32_t rel_block_begin = 0;
    std::uint32_t rel_new_end = 0;
    OriginRecordSet records;
    std::uint64_t payload = 0;
};

// Declared per-entry cost: the struct, a hash-node allowance, the key bytes, the origin array and
// both offset arrays. Size-based, so what a budget is held to does not depend on how a vector grew.
constexpr std::uint64_t kOriginCacheNodeBytes = 64;
std::uint64_t origin_entry_bytes(const std::string& bytes, const OriginRecordSet& records) {
    return sizeof(OriginCacheEntry) + kOriginCacheNodeBytes + bytes.size() +
           static_cast<std::uint64_t>(records.origins.size()) * sizeof(CachedOrigin) +
           static_cast<std::uint64_t>(records.fragments.size() + records.starts.size()) *
               sizeof(std::uint32_t);
}

struct OriginCache {
    std::unordered_multimap<std::uint64_t, std::size_t> index;
    std::deque<OriginCacheEntry> entries;   // stable addresses
    std::uint64_t bytes = 0;
    std::uint64_t not_retained = 0;
    // Planned retention evicts. An evicted slot's storage is released at once and the slot is reused by
    // the next admission, so slots never outnumber the most entries ever live together.
    std::vector<char> live;
    std::vector<std::uint64_t> demand_of;   // the demand-ledger row each slot holds
    std::vector<std::uint64_t> multiplicity;   // its requesters in the current plan: its priority
    std::vector<std::size_t> free_slots;
};

std::uint64_t origin_key_hash(const std::string& bytes, bool at_start, bool at_end,
                              std::uint32_t rel_block_begin, std::uint32_t rel_new_end);

// A window's full identity, in the order the plan ranks ties by.
using OriginIdentity = std::tuple<std::uint64_t, std::string, bool, bool, std::uint32_t, std::uint32_t>;

OriginIdentity entry_identity(const OriginCacheEntry& entry) {
    return OriginIdentity{origin_key_hash(entry.bytes, entry.at_start, entry.at_end, entry.rel_block_begin,
                                          entry.rel_new_end),
                          entry.bytes, entry.at_start, entry.at_end, entry.rel_block_begin,
                          entry.rel_new_end};
}

// The retention order, shared by the plan and by in-round replacement: more distinct requesters per
// exact byte first (compared exactly, as integers), then the identity ascending. A total order, so
// every decision taken from it is a function of the windows alone.
bool ranks_before(std::uint64_t m_a, std::uint64_t bytes_a, const OriginIdentity& id_a, std::uint64_t m_b,
                  std::uint64_t bytes_b, const OriginIdentity& id_b) {
    const unsigned __int128 lhs = static_cast<unsigned __int128>(m_a) * bytes_b;
    const unsigned __int128 rhs = static_cast<unsigned __int128>(m_b) * bytes_a;
    if (lhs != rhs) return lhs > rhs;
    return id_a < id_b;
}

// Evicts one live slot: the index entry goes, the bytes are released now, and the slot is reused by
// the next admission. Returns the payload freed.
std::uint64_t evict_origin_slot(OriginCache& cache, std::size_t at) {
    OriginCacheEntry& entry = cache.entries[at];
    const std::uint64_t hash =
        origin_key_hash(entry.bytes, entry.at_start, entry.at_end, entry.rel_block_begin, entry.rel_new_end);
    const auto range = cache.index.equal_range(hash);
    for (auto it = range.first; it != range.second; ++it) {
        if (it->second == at) {
            cache.index.erase(it);
            break;
        }
    }
    const std::uint64_t freed = entry.payload;
    cache.bytes -= freed;
    entry = OriginCacheEntry{};   // releases the key bytes and every record array
    cache.live[at] = 0;
    cache.multiplicity[at] = 0;
    cache.free_slots.push_back(at);
    return freed;
}

std::atomic<bool>& origin_self_check() {
    static std::atomic<bool> enabled{false};
    return enabled;
}
std::atomic<bool>& origin_hash_collapse() {
    static std::atomic<bool> collapse{false};
    return collapse;
}

std::uint64_t origin_key_hash(const std::string& bytes, bool at_start, bool at_end,
                              std::uint32_t rel_block_begin, std::uint32_t rel_new_end) {
    if (origin_hash_collapse().load(std::memory_order_relaxed)) return 0;
    std::uint64_t h = 1469598103934665603ULL;
    const auto mix = [&h](std::uint64_t value) { h = (h ^ value) * 1099511628211ULL; };
    for (const char c : bytes) mix(static_cast<unsigned char>(c));
    mix(0x100u | (at_start ? 1u : 0u) | (at_end ? 2u : 0u));
    mix(rel_block_begin);
    mix(rel_new_end);
    return h;
}

const OriginCacheEntry* find_origin_entry(const OriginCache& cache, std::uint64_t hash,
                                          const std::string& bytes, bool at_start, bool at_end,
                                          std::uint32_t rel_block_begin, std::uint32_t rel_new_end) {
    const auto range = cache.index.equal_range(hash);
    for (auto it = range.first; it != range.second; ++it) {
        const OriginCacheEntry& entry = cache.entries[it->second];
        // Hash chose the candidates; this decides. Flags and interval first, then every byte.
        if (entry.at_start == at_start && entry.at_end == at_end &&
            entry.rel_block_begin == rel_block_begin && entry.rel_new_end == rel_new_end &&
            entry.bytes == bytes) {
            return &entry;
        }
    }
    return nullptr;
}

bool cached_order(const CachedOrigin& lhs, const CachedOrigin& rhs) {
    if (lhs.insert_and_arm != rhs.insert_and_arm) return lhs.insert_and_arm < rhs.insert_and_arm;
    if (lhs.log_latent_mass != rhs.log_latent_mass) return lhs.log_latent_mass < rhs.log_latent_mass;
    return lhs.rel_start < rhs.rel_start;
}

// The complete window enumeration, relative and pruned to what COULD be affected. An interior
// window (neither end flag) keeps exactly the origins intersecting the relative replacement; a
// window touching an end keeps everything, because there the end-margin conditions can fire. The
// exact absolute filter still runs at every use, so correctness rests only on "retained is a superset
// of affected".
void build_origin_records(const WindowScratch& scratch, std::size_t fragment_count, std::size_t base,
                          std::uint32_t rel_block_begin, std::uint32_t rel_new_end, bool keep_all,
                          OriginRecordSet& out) {
    out.clear();
    for (std::size_t f = 0; f < fragment_count; ++f) {
        const std::size_t first = out.origins.size();
        for (std::size_t at = scratch.offsets[f]; at < scratch.offsets[f + 1]; ++at) {
            const HaploidOrigin& origin = scratch.origins[at];
            const std::uint32_t rel = static_cast<std::uint32_t>(origin.forward_start - base);
            if (!keep_all && !(rel < rel_new_end && rel_block_begin < rel + origin.insert_bp)) continue;
            CachedOrigin record;
            record.rel_start = rel;
            record.insert_and_arm = (origin.insert_bp << 1) | (origin.orientation & 1u);
            record.log_latent_mass = origin.log_latent_mass;
            out.origins.push_back(record);
        }
        if (out.origins.size() == first) continue;
        std::sort(out.origins.begin() + static_cast<std::ptrdiff_t>(first), out.origins.end(),
                  cached_order);
        out.fragments.push_back(static_cast<std::uint32_t>(f));
        out.starts.push_back(static_cast<std::uint32_t>(first));
    }
    out.starts.push_back(static_cast<std::uint32_t>(out.origins.size()));
}

// Appends one fragment's aggregate to `out`. `first` is where this fragment's states begin, so a
// run never merges into the previous fragment's last state when two fragments share one array.
void sort_and_count(std::vector<StateKey>& keys, std::vector<StateCount>& out, std::size_t first) {
    std::sort(keys.begin(), keys.end());
    for (const StateKey& key : keys) {
        if (out.size() > first && same_state(out.back().key, key)) {
            ++out.back().count;
        } else {
            out.push_back({key, 1});
        }
    }
}

// The SAME arithmetic as reduce_states, in the same order, over a sorted array instead of a map:
// ascending StateKey, one log_add_exp per key. Bit-identical because every operation is.
double reduce_sorted(const StateCount* first, const StateCount* last,
                     const DiscreteInsertDistribution& library_insert) {
    double total = kNegInf;
    for (; first != last; ++first) {
        const double term = std::log(static_cast<double>(first->count)) + kLogMateOrientationPrior +
                            library_insert.log_probability(first->key.insert_bp) +
                            first->key.log_latent_mass;
        total = log_add_exp(total, term);
    }
    return total;
}

// Linear merge of two sorted multiplicity arrays; equal keys add, exactly as ++states[key] would.
void merge_sorted(const StateCount* a, const StateCount* a_end, const StateCount* b,
                  const StateCount* b_end, std::vector<StateCount>& out) {
    out.clear();
    while (a != a_end && b != b_end) {
        if (a->key < b->key) {
            out.push_back(*a++);
        } else if (b->key < a->key) {
            out.push_back(*b++);
        } else {
            out.push_back({a->key, a->count + b->count});
            ++a;
            ++b;
        }
    }
    out.insert(out.end(), a, a_end);
    out.insert(out.end(), b, b_end);
}

// Unaffected incumbent origins, preaggregated once per exact (incumbent, block interval, length
// shift): the unaffected predicate reads the shift, so two shifts may not share a preaggregation.
struct GroupKey {
    std::size_t incumbent = 0;
    std::size_t block_begin = 0;
    std::size_t block_end = 0;
    std::ptrdiff_t shift = 0;
    bool operator<(const GroupKey& other) const {
        return std::tie(incumbent, block_begin, block_end, shift) <
               std::tie(other.incumbent, other.block_begin, other.block_end, other.shift);
    }
};
struct GroupEntry {
    std::uint64_t unaffected_origins = 0;
    // Only fragments with at least one AFFECTED incumbent origin; every other fragment's unaffected
    // set is its whole census.
    StateAggregateSet unaffected;
};

// Diagnostics only; the context leaves every pointer null.
struct OriginRecording {
    bool summary = false;
    const std::vector<char>* record_origins = nullptr;
    std::size_t unique_stride = 0;
    std::vector<OriginSummary>* summaries = nullptr;
    std::map<std::pair<std::size_t, std::size_t>, std::vector<FragmentPlacement>>* recorded = nullptr;
};

// ---- exact union origin reuse (phase 2) ------------------------------------------------------------
//
// One distinct segment's interior origins -- forward start in its forward zone, reverse end in its reverse
// zone -- as sorted per-fragment StateKey multiplicities over both arms: the census of ONE occurrence.
struct UnionOriginEntry {
    StateAggregateSet states;
    std::uint64_t origins = 0;        // physical origins one occurrence stands for (the counts' sum)
    std::uint64_t latent_pairs = 0;   // latent alignment pairs tested inside one occurrence
    std::uint64_t valid_states = 0;   // of those, valid FR states
};
// One served occurrence, in forward haplotype coordinates: a pair is served exactly when its forward start is in
// [start_lo, start_hi] and its reverse end in [end_lo, end_hi] of the SAME zone.
struct OriginZone {
    std::uint32_t start_lo = 0, start_hi = 0;
    std::uint32_t end_lo = 0, end_hi = 0;
    const UnionOriginEntry* entry = nullptr;
};
struct UnionOriginReuse {
    std::vector<std::vector<OriginZone>> zones;   // per haplotype, ascending and disjoint
    std::atomic<std::uint64_t> origins_consumed{0};
    std::atomic<std::uint64_t> states_consumed{0};
    std::uint64_t max_origins = 0;
    std::uint64_t max_states = 0;
};

// A contiguous row of placements, sorted by (start, consumed).
struct RowView {
    const FormPlacement* data = nullptr;
    std::size_t size = 0;
    const FormPlacement* begin() const { return data; }
    const FormPlacement* end() const { return data + size; }
    bool empty() const { return size == 0; }
};

// The ONE pairing of a fragment's mate placements into physical origins, shared by the reference path, the
// phase-2 direct path, the phase-2 construction and its self-check, so none of them can pair differently. Per
// arm and per forward start (ascending), every reverse placement ending in [start + min_bp, start + max_bp] is
// tested, forward rows outer and reverse rows inner, and the valid pairs are summed per insert; `emit` receives
// each origin once, in ascending insert. With `zones`, the pairs a zone serves are skipped -- the start inside
// the zone's forward range and the end inside its reverse range -- and nothing else changes: an origin fixes
// its end (e = start + insert), so it is skipped whole or not at all, and a kept origin sums the same pairs in
// the same order.
template <typename Emit>
void pair_fragment(const RowView (&forward_rows)[2], const RowView (&reverse_rows)[2], const InsertSupport& support,
                   const OriginZone* zones, std::size_t zone_count, std::map<std::uint32_t, OriginAccumulator>& origins,
                   std::uint64_t& latent_pairs, std::uint64_t& valid_states, Emit&& emit) {
    const auto start_below = [](const FormPlacement& lhs, std::size_t value) { return lhs.start < value; };
    const auto start_above = [](std::size_t value, const FormPlacement& rhs) { return value < rhs.start; };
    for (std::size_t arm = 0; arm < 2; ++arm) {
        const RowView& forward = forward_rows[arm];
        const RowView& reverse = reverse_rows[arm];
        if (forward.empty() || reverse.empty()) continue;
        std::size_t zone = 0;
        for (std::size_t first_row = 0; first_row < forward.size;) {
            const std::uint32_t forward_start = forward.data[first_row].start;
            std::size_t last_row = first_row;
            while (last_row < forward.size && forward.data[last_row].start == forward_start) ++last_row;
            const std::size_t lo_end = static_cast<std::size_t>(forward_start) + support.min_bp;
            const std::size_t hi_end = static_cast<std::size_t>(forward_start) + support.max_bp;
            const FormPlacement* begin = std::lower_bound(reverse.begin(), reverse.end(), lo_end, start_below);
            const FormPlacement* end = std::upper_bound(begin, reverse.end(), hi_end, start_above);
            const FormPlacement* skip_begin = end;
            const FormPlacement* skip_end = end;
            while (zone < zone_count && zones[zone].start_hi < forward_start) ++zone;
            if (zone < zone_count && zones[zone].start_lo <= forward_start) {
                skip_begin = std::lower_bound(begin, end, static_cast<std::size_t>(zones[zone].end_lo), start_below);
                skip_end = std::upper_bound(skip_begin, end, static_cast<std::size_t>(zones[zone].end_hi), start_above);
            }

            origins.clear();
            for (std::size_t row = first_row; row < last_row; ++row) {
                const FormPlacement& mate = forward.data[row];
                const auto pair_range = [&](const FormPlacement* from, const FormPlacement* to) {
                    for (const FormPlacement* it = from; it != to; ++it) {
                        ++latent_pairs;
                        const std::optional<std::size_t> insert = fr_insert_from_endpoints(
                            forward_start, mate.consumed, it->start, it->consumed, support);
                        if (!insert.has_value()) continue;
                        ++valid_states;
                        origins[static_cast<std::uint32_t>(*insert)].add(
                            mate.emission + it->emission, it->start - it->consumed, mate.consumed, it->consumed,
                            mate.mismatches, it->mismatches);
                    }
                };
                pair_range(begin, skip_begin);
                pair_range(skip_end, end);
            }
            for (const auto& item : origins) emit(arm, forward_start, item.first, item.second);
            first_row = last_row;
        }
    }
}

// Phase 2: every serving zone's cached census of fragment f merged into `merged_counts` (the direct census on entry),
// in StateKey order, copies kept apart by adding their multiplicities. `zone_cursor` (one per zone, reset per
// haplotype) advances monotonically with f. Shared by place_haplotype_range and the EXPERIMENTAL content index (I1.2).
void merge_served_states(std::size_t f, const OriginZone* zones, std::size_t zone_count, std::vector<std::size_t>& zone_cursor,
                         std::vector<StateCount>& merged_counts, std::vector<StateCount>& merge_scratch) {
    for (std::size_t z = 0; z < zone_count; ++z) {
        const StateAggregateSet& cached = zones[z].entry->states;
        std::size_t& at = zone_cursor[z];
        while (at < cached.fragments.size() && cached.fragments[at] < f) ++at;
        if (at == cached.fragments.size() || cached.fragments[at] != f) continue;
        const StateCount* first = cached.states.data() + cached.starts[at];
        const StateCount* last = cached.states.data() + cached.starts[at + 1];
        merge_sorted(merged_counts.data(), merged_counts.data() + merged_counts.size(), first, last, merge_scratch);
        merged_counts.swap(merge_scratch);
    }
}

// Phase 2: one distinct segment's interior origins -- forward placements (strand-local starts, (form, start, consumed)
// order) paired with reverse placements (reversed-strand-local starts, turned here to forward ends, length - start) --
// as sorted per-fragment StateKey multiplicities over both arms: the census of ONE occurrence. Every placement given
// must lie in its zone. Returns the physical origins constructed. `budget`, when given, charges and caps them.
std::uint64_t build_origin_entry(const std::vector<CachedPlacement>& forward_cached,
                                 const std::vector<CachedPlacement>& reverse_cached, std::size_t length,
                                 std::size_t form_count, std::size_t fragment_count, const PreparedFragmentSet& prepared,
                                 PlacementBudget* budget, StateAggregateSet& built_states, UnionOriginEntry& entry) {
    // Rows per form: forward placements as cached (already (start, consumed) within a form); reverse ones
    // turned to forward coordinates -- the end, length - start -- and re-sorted, as the direct path does.
    std::vector<FormPlacement> forward_rows(forward_cached.size());
    std::vector<std::size_t> forward_at(form_count + 1, 0);
    for (std::size_t i = 0; i < forward_cached.size(); ++i) {
        const CachedPlacement& c = forward_cached[i];
        forward_rows[i] = {c.start, c.consumed, c.mismatches, c.emission};
        ++forward_at[c.form + 1];
    }
    std::vector<FormPlacement> reverse_rows(reverse_cached.size());
    std::vector<std::size_t> reverse_at(form_count + 1, 0);
    for (std::size_t i = 0; i < reverse_cached.size(); ++i) {
        const CachedPlacement& c = reverse_cached[i];
        reverse_rows[i] = {static_cast<std::uint32_t>(length - c.start), c.consumed, c.mismatches, c.emission};
        ++reverse_at[c.form + 1];
    }
    for (std::size_t form = 0; form < form_count; ++form) {
        forward_at[form + 1] += forward_at[form];
        reverse_at[form + 1] += reverse_at[form];
    }
    for (std::size_t form = 0; form < form_count; ++form) {
        std::sort(reverse_rows.begin() + static_cast<std::ptrdiff_t>(reverse_at[form]),
                  reverse_rows.begin() + static_cast<std::ptrdiff_t>(reverse_at[form + 1]),
                  [](const FormPlacement& lhs, const FormPlacement& rhs) {
                      return std::tie(lhs.start, lhs.consumed) < std::tie(rhs.start, rhs.consumed);
                  });
    }
    const auto row = [](const std::vector<FormPlacement>& rows, const std::vector<std::size_t>& at, std::size_t form) {
        return RowView{rows.data() + at[form], at[form + 1] - at[form]};
    };
    built_states.clear();
    std::uint64_t constructed = 0;
    std::map<std::uint32_t, OriginAccumulator> origins;
    std::vector<StateKey> keys;
    for (std::size_t f = 0; f < fragment_count; ++f) {
        if (budget != nullptr && budget->stop()) throw PlacementCancelled();
        const InsertSupport& support = prepared.support_by_fragment[f];
        if (support.empty()) continue;
        const RowView forward_view[2] = {row(forward_rows, forward_at, 2 * f), row(forward_rows, forward_at, 2 * f + 1)};
        const RowView reverse_view[2] = {row(reverse_rows, reverse_at, 2 * f), row(reverse_rows, reverse_at, 2 * f + 1)};
        if ((forward_view[0].empty() || reverse_view[0].empty()) && (forward_view[1].empty() || reverse_view[1].empty())) {
            continue;
        }
        keys.clear();
        pair_fragment(forward_view, reverse_view, support, nullptr, 0, origins, entry.latent_pairs, entry.valid_states,
                      [&](std::size_t arm, std::uint32_t, std::uint32_t insert, const OriginAccumulator& origin) {
                          if (budget != nullptr) {
                              const std::uint64_t seen = budget->physical_origins.fetch_add(1, std::memory_order_relaxed) + 1;
                              if (seen > budget->max_physical_origins) {
                                  budget->cancel();
                                  throw ResourceExhausted(
                                      "fragment scoring reached its physical-origin cap; no partial ranking is returned");
                              }
                          }
                          StateKey key;
                          key.insert_bp = insert;
                          key.orientation = static_cast<std::uint32_t>(arm);
                          key.log_latent_mass = origin.log_latent_mass;
                          keys.push_back(key);
                      });
        if (keys.empty()) continue;
        constructed += keys.size();
        const std::size_t first = built_states.states.size();
        sort_and_count(keys, built_states.states, first);
        built_states.fragments.push_back(static_cast<std::uint32_t>(f));
        built_states.starts.push_back(static_cast<std::uint32_t>(first));
    }
    built_states.starts.push_back(static_cast<std::uint32_t>(built_states.states.size()));
    entry.states.fragments.assign(built_states.fragments.begin(), built_states.fragments.end());
    entry.states.starts.assign(built_states.starts.begin(), built_states.starts.end());
    entry.states.states.assign(built_states.states.begin(), built_states.states.end());
    entry.origins = constructed;
    return constructed;
}

// Places haplotypes [first_h, last_h) and writes their mass COLUMNS. Column storage rather than the
// row-major matrix is what makes a later batch appendable: a new haplotype adds a column and moves
// nothing, whereas a row-major matrix restrides every existing cell when the haplotype count grows.
// The arithmetic per cell is untouched, so the flattened result is identical to the matrix the
// one-shot path used to build in place.
void place_haplotype_range(
    const PlacementInputs& in,
    std::size_t first_h,
    std::size_t last_h,
    std::size_t stride,
    PlacementBudget& budget,
    std::vector<std::vector<double>>& mass_column,
    FragmentScoringCounters& counters,
    ScorerMemoryAccounts& memory,
    EmissionCacheCounters& cache_work,
    const OriginRecording& origins_out,
    const CensusCapture* census_out) {
    (void)take_gapped_emission_counters();
    const PreparedFragmentSet& prepared = *in.prepared;
    const std::vector<ReadPair>& fragments = *in.fragments;
    const LibraryModel& library = *in.library;
    const DiscreteInsertDistribution& library_insert = *in.library_insert;
    const FragmentScoringOptions& options = *in.options;
    const std::size_t band = in.band;
    const auto form_slot = [](std::size_t fragment, std::size_t arm) { return fragment * 2 + arm; };

    std::map<StateKey, std::uint64_t> states;
    std::map<std::uint32_t, OriginAccumulator> origins;
    std::vector<std::vector<FormPlacement>> forward_placements(prepared.forward_forms.size());
    std::vector<std::vector<FormPlacement>> reverse_placements(prepared.reverse_forms.size());
    std::vector<std::pair<std::uint32_t, std::uint32_t>> proposed_scratch;
    std::vector<std::uint8_t> replay_skip;
    // Phase 2 worker scratch: the direct origins' keys, their census, the merged census, and per-zone cursors.
    std::vector<StateKey> direct_keys;
    std::vector<StateCount> direct_counts, merged_counts, merge_scratch;
    std::vector<std::size_t> zone_cursor;
    const bool reuse_self_check = in.origin_reuse != nullptr && origin_self_check().load(std::memory_order_relaxed);

    for (std::size_t at = first_h; at < last_h; at += stride) {
        const std::size_t h = in.slot_list != nullptr ? (*in.slot_list)[at] : at;
        if (budget.stop()) throw PlacementCancelled();
        const std::string& haplotype = (*in.sequences)[h];
        if (prepared.forward_forms.empty()) continue;
        // MEASUREMENT ONLY: this worker's counters before the haplotype (read, never changed).
        const bool trace_work = options.record_haplotype_work != nullptr && census_out == nullptr;
        const FragmentScoringCounters work_before = trace_work ? counters : FragmentScoringCounters{};
        const std::uint64_t executed_before = trace_work ? cache_work.executed_dps : 0;
        const auto work_mark = std::chrono::steady_clock::now();

        for (std::vector<FormPlacement>& row : forward_placements) row.clear();
        const std::uint64_t rows_before = counters.accepted_mate_placements;
        SeedingTimes seeding_times;
        SeedingTimes* const timed = trace_work ? &seeding_times : nullptr;
        collect_form_placements(haplotype, prepared.forward_forms, prepared.forward_componentise,
                                prepared.forward_indexes, budget, band, library, forward_placements,
                                proposed_scratch, counters, cache_work, options, nullptr,
                                in.replay != nullptr ? &in.replay->forward[h] : nullptr, &replay_skip, in.replay,
                                timed, in.prefilter_cache);
        memory.proposed_starts_peak =
            std::max<std::uint64_t>(memory.proposed_starts_peak, proposed_scratch.size());
        memory.proposed_starts_capacity_bytes_peak = std::max<std::uint64_t>(
            memory.proposed_starts_capacity_bytes_peak,
            proposed_scratch.capacity() * sizeof(std::pair<std::uint32_t, std::uint32_t>));

        const std::string& reversed_haplotype = (*in.reversed)[h];
        for (std::vector<FormPlacement>& row : reverse_placements) row.clear();
        collect_form_placements(reversed_haplotype, prepared.reverse_forms,
                                prepared.reverse_componentise, prepared.reverse_indexes, budget, band,
                                library, reverse_placements, proposed_scratch, counters, cache_work,
                                options, nullptr, in.replay != nullptr ? &in.replay->reverse[h] : nullptr,
                                &replay_skip, in.replay, timed, in.prefilter_cache);
        memory.proposed_starts_peak =
            std::max<std::uint64_t>(memory.proposed_starts_peak, proposed_scratch.size());
        memory.proposed_starts_capacity_bytes_peak = std::max<std::uint64_t>(
            memory.proposed_starts_capacity_bytes_peak,
            proposed_scratch.capacity() * sizeof(std::pair<std::uint32_t, std::uint32_t>));
        memory.placement_rows_peak = std::max<std::uint64_t>(
            memory.placement_rows_peak, counters.accepted_mate_placements - rows_before);
        const auto work_seeded = std::chrono::steady_clock::now();

        for (std::vector<FormPlacement>& list : reverse_placements) {
            for (FormPlacement& placement : list) {
                placement.start = static_cast<std::uint32_t>(haplotype.size() - placement.start);
            }
            std::sort(list.begin(), list.end(),
                      [](const FormPlacement& lhs, const FormPlacement& rhs) {
                          return std::tie(lhs.start, lhs.consumed) <
                                 std::tie(rhs.start, rhs.consumed);
                      });
        }

        HaploidOriginCensus* census = nullptr;
        if (census_out != nullptr) {
            const auto found = census_out->wanted.find(h);
            if (found != census_out->wanted.end()) {
                census = found->second;
                census->haplotype_length = haplotype.size();
                census->origins_by_fragment.assign(fragments.size(), {});
            }
        }

        // Phase 2: the occurrences this haplotype is served, charged whole before any of them is merged, so a
        // consumption cap refuses before a partial census exists.
        const std::vector<OriginZone>* origin_zones =
            in.origin_reuse != nullptr ? &in.origin_reuse->zones[h] : nullptr;
        const bool reusing = origin_zones != nullptr && !origin_zones->empty();
        double merge_seconds = 0.0;
        if (reusing) {
            std::uint64_t served_origins = 0, served_states = 0;
            for (const OriginZone& zone : *origin_zones) {
                served_origins += zone.entry->origins;
                served_states += zone.entry->states.states.size();
                counters.reused_latent_alignment_pairs += zone.entry->latent_pairs;
                counters.reused_valid_fr_states += zone.entry->valid_states;
            }
            counters.reused_physical_origins += served_origins;
            counters.reused_state_aggregates += served_states;
            const std::uint64_t origins_total =
                in.origin_reuse->origins_consumed.fetch_add(served_origins, std::memory_order_relaxed) + served_origins;
            const std::uint64_t states_total =
                in.origin_reuse->states_consumed.fetch_add(served_states, std::memory_order_relaxed) + served_states;
            if (origins_total > in.origin_reuse->max_origins) {
                budget.cancel();
                throw ResourceExhausted("fragment scoring reached its cached-origin consumption cap; no partial ranking "
                                        "is returned");
            }
            if (states_total > in.origin_reuse->max_states) {
                budget.cancel();
                throw ResourceExhausted("fragment scoring reached its cached-state consumption cap; no partial ranking "
                                        "is returned");
            }
            zone_cursor.assign(origin_zones->size(), 0);
        }

        for (std::size_t f = 0; f < fragments.size(); ++f) {
            if (budget.stop()) throw PlacementCancelled();
            const InsertSupport& support = prepared.support_by_fragment[f];
            if (support.empty()) continue;
            const RowView forward_rows[2] = {
                {forward_placements[form_slot(f, 0)].data(), forward_placements[form_slot(f, 0)].size()},
                {forward_placements[form_slot(f, 1)].data(), forward_placements[form_slot(f, 1)].size()}};
            const RowView reverse_rows[2] = {
                {reverse_placements[form_slot(f, 0)].data(), reverse_placements[form_slot(f, 0)].size()},
                {reverse_placements[form_slot(f, 1)].data(), reverse_placements[form_slot(f, 1)].size()}};

            if (reusing) {
                // Direct: every pair no zone serves, charged as executed work.
                direct_keys.clear();
                pair_fragment(forward_rows, reverse_rows, support, origin_zones->data(), origin_zones->size(), origins,
                              counters.latent_alignment_pairs, counters.valid_fr_states,
                              [&](std::size_t arm, std::uint32_t, std::uint32_t insert, const OriginAccumulator& origin) {
                                  ++counters.physical_origins;
                                  const std::uint64_t seen =
                                      budget.physical_origins.fetch_add(1, std::memory_order_relaxed) + 1;
                                  if (seen > budget.max_physical_origins) {
                                      budget.cancel();
                                      throw ResourceExhausted(
                                          "fragment scoring reached its physical-origin cap; no partial "
                                          "ranking is returned");
                                  }
                                  if (budget.stop()) throw PlacementCancelled();
                                  StateKey key;
                                  key.insert_bp = insert;
                                  key.orientation = static_cast<std::uint32_t>(arm);
                                  key.log_latent_mass = origin.log_latent_mass;
                                  direct_keys.push_back(key);
                              });
                direct_counts.clear();
                sort_and_count(direct_keys, direct_counts, 0);
                counters.statistical_state_classes += direct_counts.size();
                // Served: each occurrence's census added in, copies kept apart by adding their multiplicities.
                const auto merge_mark = std::chrono::steady_clock::now();
                merged_counts.assign(direct_counts.begin(), direct_counts.end());
                merge_served_states(f, origin_zones->data(), origin_zones->size(), zone_cursor, merged_counts, merge_scratch);
                merge_seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - merge_mark).count();
                if (reuse_self_check) {
                    // Fresh enumeration of the whole fragment, no zone skipped, reduced as the reference reduces.
                    std::uint64_t unused_pairs = 0, unused_states = 0;
                    states.clear();
                    pair_fragment(forward_rows, reverse_rows, support, nullptr, 0, origins, unused_pairs, unused_states,
                                  [&](std::size_t arm, std::uint32_t, std::uint32_t insert, const OriginAccumulator& origin) {
                                      StateKey key;
                                      key.insert_bp = insert;
                                      key.orientation = static_cast<std::uint32_t>(arm);
                                      key.log_latent_mass = origin.log_latent_mass;
                                      ++states[key];
                                  });
                    bool same = states.size() == merged_counts.size();
                    std::size_t k = 0;
                    for (auto it = states.begin(); same && it != states.end(); ++it, ++k) {
                        const StateCount& merged = merged_counts[k];
                        same = merged.key.insert_bp == it->first.insert_bp &&
                               merged.key.orientation == it->first.orientation &&
                               std::memcmp(&merged.key.log_latent_mass, &it->first.log_latent_mass, sizeof(double)) == 0 &&
                               merged.count == it->second;
                    }
                    if (!same) {
                        throw std::logic_error("union origin reuse self-check: the merged census of fragment " +
                                               std::to_string(f) + " on haplotype " + std::to_string(h) +
                                               " differs from a fresh enumeration");
                    }
                }
                if (merged_counts.empty()) continue;
                counters.merged_state_classes += merged_counts.size();
                mass_column[h][f] = reduce_sorted(merged_counts.data(), merged_counts.data() + merged_counts.size(),
                                                  library_insert);
                continue;
            }

            states.clear();
            pair_fragment(forward_rows, reverse_rows, support, nullptr, 0, origins, counters.latent_alignment_pairs,
                          counters.valid_fr_states,
                          [&](std::size_t arm, std::uint32_t forward_start, std::uint32_t insert,
                              const OriginAccumulator& origin) {
                        ++counters.physical_origins;
                        const std::uint64_t seen =
                            budget.physical_origins.fetch_add(1, std::memory_order_relaxed) + 1;
                        if (seen > budget.max_physical_origins) {
                            budget.cancel();
                            throw ResourceExhausted(
                                "fragment scoring reached its physical-origin cap; no partial "
                                "ranking is returned");
                        }
                        if (budget.stop()) throw PlacementCancelled();
                        StateKey key;
                        key.insert_bp = insert;
                        key.orientation = static_cast<std::uint32_t>(arm);
                        key.log_latent_mass = origin.log_latent_mass;
                        ++states[key];
                        if (census != nullptr) {
                            census->origins_by_fragment[f].push_back(
                                {forward_start, insert, static_cast<std::uint32_t>(arm), origin.log_latent_mass});
                        }
                        const bool want_named =
                            origins_out.record_origins != nullptr &&
                            (*origins_out.record_origins)[f] != 0;
                        if (origins_out.summary || want_named) {
                            FragmentPlacement placement;
                            placement.forward_start = forward_start;
                            placement.insert_bp = insert;
                            placement.log_latent_mass = origin.log_latent_mass;
                            placement.latent_alignments = origin.latent_alignments;
                            placement.reverse_start = origin.reverse_start;
                            placement.forward_mismatches = origin.forward_mismatches;
                            placement.reverse_mismatches = origin.reverse_mismatches;
                            placement.forward_consumed_bp = origin.forward_consumed_bp;
                            placement.reverse_consumed_bp = origin.reverse_consumed_bp;
                            placement.orientation = arm == 0 ? MateOrientation::FirstForward
                                                             : MateOrientation::SecondForward;
                            if (origins_out.summary && origins_out.summaries != nullptr) {
                                OriginSummary& summary =
                                    (*origins_out.summaries)[f * origins_out.unique_stride + h];
                                ++summary.origins;
                                const double term = kLogMateOrientationPrior +
                                                    library_insert.log_probability(insert) +
                                                    origin.log_latent_mass;
                                if (term > summary.best_log_term) {
                                    summary.best_log_term = term;
                                    summary.best = placement;
                                }
                            }
                            if (want_named && origins_out.recorded != nullptr) {
                                (*origins_out.recorded)[{f, h}].push_back(placement);
                            }
                        }
                    });
            if (states.empty()) continue;
            counters.statistical_state_classes += states.size();
            counters.merged_state_classes += states.size();
            mass_column[h][f] = reduce_states(states, library_insert);
        }
        if (trace_work && h < options.record_haplotype_work->size()) {
            HaplotypeWork& work = (*options.record_haplotype_work)[h];   // one worker per haplotype: no race
            const auto work_done = std::chrono::steady_clock::now();
            work.recorded = true;
            work.haplotype_bp = haplotype.size();
            work.piece_hits = counters.piece_hits - work_before.piece_hits;
            work.proposed_starts = counters.proposed_starts - work_before.proposed_starts;
            work.unique_starts = counters.unique_starts - work_before.unique_starts;
            work.full_read_verifications = counters.full_read_verifications - work_before.full_read_verifications;
            work.prefilter_rejections = counters.prefilter_rejections - work_before.prefilter_rejections;
            work.dp_evaluations = counters.dp_evaluations - work_before.dp_evaluations;
            work.executed_dps = cache_work.executed_dps - executed_before;
            work.accepted_mate_placements = counters.accepted_mate_placements - work_before.accepted_mate_placements;
            work.physical_origins = counters.physical_origins - work_before.physical_origins;
            work.statistical_state_classes = counters.statistical_state_classes - work_before.statistical_state_classes;
            work.reused_physical_origins = counters.reused_physical_origins - work_before.reused_physical_origins;
            work.reused_state_aggregates = counters.reused_state_aggregates - work_before.reused_state_aggregates;
            work.merged_state_classes = counters.merged_state_classes - work_before.merged_state_classes;
            work.seconds_cached_merge = merge_seconds;
            work.seconds_seed_scan = static_cast<double>(seeding_times.scan_ns) * 1e-9;
            work.seconds_start_sort = static_cast<double>(seeding_times.sort_ns) * 1e-9;
            work.seconds_verification = static_cast<double>(seeding_times.verify_ns) * 1e-9;
            work.seconds_replay = static_cast<double>(seeding_times.replay_ns) * 1e-9;
            work.seconds_row_sort = static_cast<double>(seeding_times.row_sort_ns) * 1e-9;
            work.seconds_seeding_and_emission = std::chrono::duration<double>(work_seeded - work_mark).count();
            work.seconds_pairing_and_reduction = std::chrono::duration<double>(work_done - work_seeded).count();
        }
    }
    {
        std::uint64_t capacity = 0;
        for (const std::vector<FormPlacement>& row : forward_placements) capacity += row.capacity();
        for (const std::vector<FormPlacement>& row : reverse_placements) capacity += row.capacity();
        memory.placement_capacity_bytes_peak = std::max<std::uint64_t>(
            memory.placement_capacity_bytes_peak, capacity * sizeof(FormPlacement));
    }
}

// Exact union reuse, stage one: every distinct (strand, segment bytes) with a non-empty interior zone that
// occurs at least twice, filled in a canonical order -- most occurrences first, then shorter, then strand and
// first occurrence -- and admitted in that order while it fits the byte budget, so what is retained never
// depends on which worker finished first. A declined segment is simply placed directly. Filling counts its
// executed seeding, verification and emission work; its placements are not accepted placements until replayed.
void fill_union_cache(const PlacementInputs& in, const std::vector<std::vector<std::uint32_t>>& lengths,
                      std::size_t budget_bytes, std::uint64_t max_consumed, std::size_t threads,
                      PlacementBudget& budget, FragmentScoringCounters& counters, EmissionCacheCounters& cache_work,
                      std::vector<UnionSegmentEntry>& entries, UnionReplay& replay,
                      UnionPlacementCacheReport& report,
                      // Phase 2, when `origin_reuse` is given: origins constructed from the admitted placements.
                      std::size_t origin_budget_bytes = 0, UnionOriginReuse* origin_reuse = nullptr,
                      std::deque<UnionOriginEntry>* origin_entries = nullptr,
                      UnionOriginCacheReport* origin_report = nullptr,
                      // MEASUREMENT ONLY: per eligible segment (canonical order), its fill's wall seconds and length.
                      std::vector<std::pair<double, std::size_t>>* segment_seconds = nullptr) {
    const auto fill_mark = std::chrono::steady_clock::now();
    const PreparedFragmentSet& prepared = *in.prepared;
    const LibraryModel& library = *in.library;
    std::size_t R = 0;
    for (const auto& form : prepared.forward_forms) R = std::max(R, form.bases.size());
    for (const auto& form : prepared.reverse_forms) R = std::max(R, form.bases.size());
    const std::size_t g = library.max_gap_bp;
    report.active = true;
    report.read_span_bp = R;
    report.max_gap_bp = g;
    report.budget_bytes = budget_bytes;
    report.max_placements_consumed = max_consumed;
    const std::size_t unique = in.sequences->size();
    if (lengths.size() != unique) throw std::invalid_argument("the union segmentation needs one entry per unique haplotype");
    replay.forward.assign(unique, {});
    replay.reverse.assign(unique, {});
    replay.max_consumed = max_consumed;

    struct Occurrence {
        std::size_t haplotype = 0;
        std::uint32_t offset = 0;   // in the strand's own coordinates
    };
    struct Distinct {
        int strand = 0;
        std::string bytes;
        std::vector<Occurrence> occurrences;
        std::size_t partner = static_cast<std::size_t>(-1);   // a forward segment's reverse strand
    };
    std::vector<Distinct> distinct;
    std::unordered_multimap<std::uint64_t, std::size_t> by_hash;
    const bool collapse = union_hash_collapse().load(std::memory_order_relaxed);
    const std::size_t minimum = R + 2 * g + 1;   // a segment this long has at least one interior start
    for (std::size_t h = 0; h < unique; ++h) {
        const std::string& forward = (*in.sequences)[h];
        std::uint64_t sum = 0;
        for (const std::uint32_t length : lengths[h]) sum += length;
        if (sum != forward.size()) throw std::invalid_argument("union segment lengths must sum to the haplotype length");
        std::uint32_t offset = 0;
        for (const std::uint32_t length : lengths[h]) {
            if (length >= minimum) {
                std::size_t forward_distinct = 0;
                for (int strand = 0; strand < 2; ++strand) {
                    const std::uint32_t at = strand == 0 ? offset : static_cast<std::uint32_t>(forward.size() - offset - length);
                    const std::string& text = strand == 0 ? forward : (*in.reversed)[h];
                    std::string bytes = text.substr(at, length);
                    ++report.segment_occurrences;
                    report.zone_starts_requested += length - minimum + 1;
                    std::uint64_t key = 0;
                    if (!collapse) {
                        key = 1469598103934665603ULL;
                        for (const char base : bytes) key = (key ^ static_cast<unsigned char>(base)) * 1099511628211ULL;
                        key = (key ^ static_cast<std::uint64_t>(strand + 1)) * 1099511628211ULL;
                    }
                    std::size_t found = distinct.size();
                    const auto range = by_hash.equal_range(key);
                    for (auto it = range.first; it != range.second; ++it) {
                        if (distinct[it->second].strand == strand && distinct[it->second].bytes == bytes) found = it->second;
                    }
                    if (found == distinct.size()) {
                        distinct.push_back({strand, std::move(bytes), {}});
                        by_hash.emplace(key, found);
                    }
                    distinct[found].occurrences.push_back({h, at});
                    // The reverse strand's bytes are a function of the forward bytes, so a forward segment has
                    // exactly one partner, and the two share their occurrences mirrored.
                    if (strand == 0) forward_distinct = found;
                    else distinct[forward_distinct].partner = found;
                }
            }
            offset += length;
        }
    }
    report.distinct_segments = distinct.size();
    for (const Distinct& d : distinct) report.index_bytes += d.bytes.size() + sizeof(Distinct) + d.occurrences.size() * sizeof(Occurrence);

    std::vector<std::size_t> eligible;
    for (std::size_t i = 0; i < distinct.size(); ++i) {
        if (distinct[i].occurrences.size() >= 2) eligible.push_back(i);
    }
    std::sort(eligible.begin(), eligible.end(), [&](std::size_t x, std::size_t y) {
        const Distinct& a = distinct[x];
        const Distinct& b = distinct[y];
        if (a.occurrences.size() != b.occurrences.size()) return a.occurrences.size() > b.occurrences.size();
        if (a.bytes.size() != b.bytes.size()) return a.bytes.size() < b.bytes.size();
        if (a.strand != b.strand) return a.strand < b.strand;
        return std::tie(a.occurrences.front().haplotype, a.occurrences.front().offset) <
               std::tie(b.occurrences.front().haplotype, b.occurrences.front().offset);
    });
    report.eligible_segments = eligible.size();
    entries.assign(eligible.size(), UnionSegmentEntry{});
    if (segment_seconds != nullptr) segment_seconds->assign(eligible.size(), {0.0, 0});
    std::vector<char> admitted(eligible.size(), 0);

    // Bounded batches of canonical segments. Within one, the workers claim segments longest first from a shared
    // cursor, so a long segment no longer holds seven idle workers; each result lands in its canonical slot, the
    // counters are merged in slot order and admission runs in canonical order after the batch, so what is
    // computed, counted and retained does not depend on which worker took what. The batch bound caps what is
    // held before admission decides.
    constexpr std::size_t kFillBatch = 64;
    const std::size_t workers = std::max<std::size_t>(1, std::min(threads, eligible.size()));
    std::uint64_t held = 0;
    report.peak_bytes = report.index_bytes;
    for (std::size_t begin = 0; begin < eligible.size(); begin += kFillBatch) {
        const std::size_t end = std::min(begin + kFillBatch, eligible.size());
        std::vector<FragmentScoringCounters> local(end - begin);
        std::vector<EmissionCacheCounters> local_cache(end - begin);
        std::vector<std::exception_ptr> failures(end - begin);
        const auto fill_one = [&](std::size_t slot) {
            const ActiveWorker active;
            const auto segment_mark = std::chrono::steady_clock::now();
            struct Record {
                std::vector<std::pair<double, std::size_t>>* into;
                std::size_t at;
                std::size_t length;
                std::chrono::steady_clock::time_point mark;
                ~Record() {
                    if (into != nullptr) {
                        (*into)[at] = {std::chrono::duration<double>(std::chrono::steady_clock::now() - mark).count(), length};
                    }
                }
            } record{segment_seconds, begin + slot, distinct[eligible[begin + slot]].bytes.size(), segment_mark};
            const Distinct& d = distinct[eligible[begin + slot]];
            const bool forward = d.strand == 0;
            const auto& forms = forward ? prepared.forward_forms : prepared.reverse_forms;
            std::vector<std::vector<FormPlacement>> rows(forms.size());
            std::vector<std::pair<std::uint32_t, std::uint32_t>> proposed;
            (void)take_gapped_emission_counters();
            collect_form_placements(d.bytes, forms, forward ? prepared.forward_componentise : prepared.reverse_componentise,
                                    forward ? prepared.forward_indexes : prepared.reverse_indexes, budget, in.band,
                                    library, rows, proposed, local[slot], local_cache[slot], *in.options, nullptr,
                                    nullptr, nullptr, nullptr, nullptr, in.prefilter_cache);
            const std::size_t length = d.bytes.size();
            std::vector<CachedPlacement>& out = entries[begin + slot].placements;
            for (std::size_t f = 0; f < rows.size(); ++f) {
                for (const FormPlacement& p : rows[f]) {
                    if (p.start < g || p.start + R + g > length) continue;   // outside the interior zone
                    out.push_back({static_cast<std::uint32_t>(f), p.start, p.consumed, p.mismatches, p.emission});
                }
            }
        };
        std::vector<std::size_t> claim_order(end - begin);
        for (std::size_t slot = 0; slot < claim_order.size(); ++slot) claim_order[slot] = slot;
        std::stable_sort(claim_order.begin(), claim_order.end(), [&](std::size_t x, std::size_t y) {
            return distinct[eligible[begin + x]].bytes.size() > distinct[eligible[begin + y]].bytes.size();
        });
        const std::size_t batch_workers = std::min(workers, end - begin);
        if (batch_workers == 1) {
            for (const std::size_t slot : claim_order) fill_one(slot);
        } else {
            std::atomic<std::size_t> cursor{0};
            std::vector<std::thread> pool;
            for (std::size_t w = 0; w < batch_workers; ++w) {
                pool.emplace_back([&]() {
                    for (std::size_t at = cursor.fetch_add(1); at < claim_order.size(); at = cursor.fetch_add(1)) {
                        const std::size_t slot = claim_order[at];
                        try {
                            fill_one(slot);
                        } catch (const PlacementCancelled&) {
                            return;
                        } catch (...) {
                            failures[slot] = std::current_exception();
                            budget.cancel();
                            return;
                        }
                    }
                });
            }
            for (std::thread& worker : pool) worker.join();
            for (const std::exception_ptr& failure : failures) {
                if (failure) std::rethrow_exception(failure);
            }
            if (budget.stop()) throw PlacementCancelled();
        }
        std::uint64_t transient = 0;
        for (std::size_t slot = 0; slot < end - begin; ++slot) {
            const FragmentScoringCounters& part = local[slot];
            counters.candidate_bases_scanned += part.candidate_bases_scanned;
            counters.piece_hits += part.piece_hits;
            counters.proposed_starts += part.proposed_starts;
            counters.unique_starts += part.unique_starts;
            counters.full_read_verifications += part.full_read_verifications;
            counters.dp_evaluations += part.dp_evaluations;
            counters.prefilter_rejections += part.prefilter_rejections;
            report.placements_enumerated += part.accepted_mate_placements;
            const EmissionCacheCounters& e = local_cache[slot];
            cache_work.executed_prefilters += e.executed_prefilters;
            cache_work.dp_cache_requests += e.dp_cache_requests;
            cache_work.dp_cache_hits += e.dp_cache_hits;
            cache_work.dp_cache_misses += e.dp_cache_misses;
            cache_work.executed_dps += e.executed_dps;
            transient += entries[begin + slot].placements.size() * sizeof(CachedPlacement);
        }
        report.peak_bytes = std::max<std::uint64_t>(report.peak_bytes, report.index_bytes + held + transient);
        for (std::size_t slot = 0; slot < end - begin; ++slot) {
            std::vector<CachedPlacement>& placements = entries[begin + slot].placements;
            const std::uint64_t bytes = placements.size() * sizeof(CachedPlacement);
            if (held + bytes <= budget_bytes) {
                admitted[begin + slot] = 1;
                held += bytes;
                ++report.admitted_segments;
                report.placements_cached += placements.size();
            } else {
                ++report.declined_segments;
                std::vector<CachedPlacement>().swap(placements);
            }
        }
    }
    report.payload_bytes = held;

    for (std::size_t e = 0; e < eligible.size(); ++e) {
        if (admitted[e] == 0) continue;
        const Distinct& d = distinct[eligible[e]];
        const std::uint32_t length = static_cast<std::uint32_t>(d.bytes.size());
        for (const Occurrence& occurrence : d.occurrences) {
            ReplayZone zone;
            zone.lo = occurrence.offset + static_cast<std::uint32_t>(g);
            zone.hi = occurrence.offset + length - static_cast<std::uint32_t>(R + g);
            zone.shift = occurrence.offset;
            zone.entry = &entries[e];
            (d.strand == 0 ? replay.forward : replay.reverse)[occurrence.haplotype].push_back(zone);
            ++report.hits;
            report.zone_starts_served += length - minimum + 1;
        }
    }
    for (auto* plan : {&replay.forward, &replay.reverse}) {
        for (std::vector<ReplayZone>& zones : *plan) {
            std::sort(zones.begin(), zones.end(), [](const ReplayZone& x, const ReplayZone& y) { return x.lo < y.lo; });
        }
    }
    report.misses = report.segment_occurrences - report.hits;
    if (origin_reuse == nullptr) {
        report.seconds_fill = std::chrono::duration<double>(std::chrono::steady_clock::now() - fill_mark).count();
        return;
    }
    // ---- phase 2: interior origins, constructed from the admitted placements of both strands ---------------
    const auto origin_mark = std::chrono::steady_clock::now();
    UnionOriginCacheReport& oreport = *origin_report;
    oreport.active = true;
    oreport.budget_bytes = origin_budget_bytes;
    oreport.max_origins_consumed = origin_reuse->max_origins;
    oreport.max_states_consumed = origin_reuse->max_states;
    origin_reuse->zones.assign(unique, {});
    std::vector<std::size_t> position_of(distinct.size(), static_cast<std::size_t>(-1));
    for (std::size_t e = 0; e < eligible.size(); ++e) position_of[eligible[e]] = e;
    // Candidates in the canonical order phase 1 used, forward strands only; both strands must be admitted.
    std::vector<std::size_t> candidates;   // positions in `eligible`
    for (std::size_t e = 0; e < eligible.size(); ++e) {
        const Distinct& d = distinct[eligible[e]];
        if (d.strand != 0) continue;
        oreport.segment_occurrences += d.occurrences.size();
        if (admitted[e] == 0 || d.partner == static_cast<std::size_t>(-1)) continue;
        const std::size_t partner = position_of[d.partner];
        if (partner == static_cast<std::size_t>(-1) || admitted[partner] == 0) continue;
        candidates.push_back(e);
    }
    for (const Distinct& d : distinct) {
        if (d.strand == 0 && d.occurrences.size() < 2) oreport.segment_occurrences += d.occurrences.size();
    }
    oreport.eligible_segments = candidates.size();
    origin_entries->clear();
    std::vector<UnionOriginEntry*> built(candidates.size(), nullptr);
    std::vector<char> kept(candidates.size(), 0);
    const std::size_t form_count = prepared.forward_forms.size();
    const std::size_t fragment_count = in.fragments->size();
    const std::size_t oworkers = std::max<std::size_t>(1, std::min(threads, candidates.size()));
    std::uint64_t oheld = 0;
    oreport.peak_bytes = 0;
    const auto entry_bytes = [](const UnionOriginEntry& entry) -> std::uint64_t {
        return sizeof(UnionOriginEntry) +
               static_cast<std::uint64_t>(entry.states.fragments.size() + entry.states.starts.size()) *
                   sizeof(std::uint32_t) +
               static_cast<std::uint64_t>(entry.states.states.size()) * sizeof(StateCount);
    };
    // Per-slot scratch, reused across batches: an entry is built here and copied out once at its exact size, so
    // the fill leaves no trail of doubled and abandoned arrays behind the retained ones.
    std::vector<StateAggregateSet> slot_scratch(oworkers);
    for (std::size_t begin = 0; begin < candidates.size(); begin += oworkers) {
        const std::size_t end = std::min(begin + oworkers, candidates.size());
        std::vector<UnionOriginEntry> batch(end - begin);
        std::vector<std::uint64_t> batch_origins(end - begin, 0);
        std::vector<std::exception_ptr> failures(end - begin);
        const auto build_one = [&](std::size_t slot) {
            const ActiveWorker active;
            const std::size_t e = candidates[begin + slot];
            const Distinct& d = distinct[eligible[e]];
            const std::size_t length = d.bytes.size();
            const std::vector<CachedPlacement>& forward_cached = entries[e].placements;
            const std::vector<CachedPlacement>& reverse_cached = entries[position_of[d.partner]].placements;
            UnionOriginEntry& entry = batch[slot];
            batch_origins[slot] = build_origin_entry(forward_cached, reverse_cached, length, form_count, fragment_count,
                                                     prepared, &budget, slot_scratch[slot], entry);
        };
        if (end - begin == 1) {
            build_one(0);
        } else {
            std::vector<std::thread> pool;
            for (std::size_t slot = 0; slot < end - begin; ++slot) {
                pool.emplace_back([&, slot]() {
                    try {
                        build_one(slot);
                    } catch (const PlacementCancelled&) {
                    } catch (...) {
                        failures[slot] = std::current_exception();
                    }
                });
            }
            for (std::thread& worker : pool) worker.join();
            for (const std::exception_ptr& failure : failures) {
                if (failure) std::rethrow_exception(failure);
            }
            if (budget.stop()) throw PlacementCancelled();
        }
        // Executed work, merged in slot order; then admission in canonical order.
        std::uint64_t transient = 0;
        for (std::size_t slot = 0; slot < end - begin; ++slot) {
            const UnionOriginEntry& entry = batch[slot];
            counters.physical_origins += entry.origins;
            counters.statistical_state_classes += entry.states.states.size();
            counters.latent_alignment_pairs += entry.latent_pairs;
            counters.valid_fr_states += entry.valid_states;
            oreport.origins_constructed += entry.origins;
            oreport.state_aggregates_constructed += entry.states.states.size();
            oreport.latent_pairs_constructed += entry.latent_pairs;
            oreport.valid_states_constructed += entry.valid_states;
            transient += entry_bytes(entry);
        }
        oreport.peak_bytes = std::max<std::uint64_t>(oreport.peak_bytes, oheld + transient);
        for (std::size_t slot = 0; slot < end - begin; ++slot) {
            UnionOriginEntry& entry = batch[slot];
            const std::uint64_t bytes = entry_bytes(entry);
            if (oheld + bytes <= origin_budget_bytes) {
                oheld += bytes;
                kept[begin + slot] = 1;
                ++oreport.admitted_segments;
                oreport.origin_records_held += entry.states.states.size();
                oreport.origins_held += entry.origins;
                origin_entries->push_back(std::move(entry));
                built[begin + slot] = &origin_entries->back();
            } else {
                ++oreport.declined_segments;
            }
        }
    }
    oreport.payload_bytes = oheld;
    for (std::size_t c = 0; c < candidates.size(); ++c) {
        if (kept[c] == 0) continue;
        const Distinct& d = distinct[eligible[candidates[c]]];
        const std::uint32_t length = static_cast<std::uint32_t>(d.bytes.size());
        for (const Occurrence& occurrence : d.occurrences) {
            OriginZone zone;
            zone.start_lo = occurrence.offset + static_cast<std::uint32_t>(g);
            zone.start_hi = occurrence.offset + length - static_cast<std::uint32_t>(R + g);
            zone.end_lo = occurrence.offset + static_cast<std::uint32_t>(R + g);
            zone.end_hi = occurrence.offset + length - static_cast<std::uint32_t>(g);
            zone.entry = built[c];
            origin_reuse->zones[occurrence.haplotype].push_back(zone);
            ++oreport.hits;
        }
    }
    for (std::vector<OriginZone>& zones : origin_reuse->zones) {
        std::sort(zones.begin(), zones.end(),
                  [](const OriginZone& x, const OriginZone& y) { return x.start_lo < y.start_lo; });
        oreport.index_bytes += zones.size() * sizeof(OriginZone);
    }
    oreport.peak_bytes = std::max<std::uint64_t>(oreport.peak_bytes, oheld + oreport.index_bytes);
    oreport.misses = oreport.segment_occurrences - oreport.hits;
    oreport.seconds_fill = std::chrono::duration<double>(std::chrono::steady_clock::now() - origin_mark).count();
    report.seconds_fill = std::chrono::duration<double>(origin_mark - fill_mark).count();
}

// Runs placement over [0, unique) with the declared worker count, merging in thread-index order so
// the totals are the serial totals exactly. Shared by the one-shot path and every incremental batch,
// which is what makes their counters comparable rather than merely similar.
std::size_t run_placement(
    const PlacementInputs& in,
    std::size_t first_h,
    std::size_t last_h,
    std::size_t requested_threads,
    bool force_serial,
    PlacementBudget& budget,
    std::vector<std::vector<double>>& mass_column,
    FragmentScoringCounters& counters,
    ScorerMemoryAccounts& memory,
    EmissionCacheCounters& cache_work,
    const OriginRecording& origins_out,
    const CensusCapture* census_out = nullptr) {
    const std::size_t span_count = last_h - first_h;
    const std::size_t threads =
        (span_count == 0 || force_serial)
            ? 1
            : std::max<std::size_t>(1, std::min(requested_threads, span_count));
    if (threads <= 1) {
        const ActiveWorker active;
        place_haplotype_range(in, first_h, last_h, 1, budget, mass_column, counters, memory, cache_work,
                              origins_out, census_out);
        return threads;
    }
    const std::size_t span = (span_count + threads - 1) / threads;
    // Under union origin reuse the haplotypes are dealt round-robin rather than in contiguous spans: served
    // haplotypes cost a fraction of unserved ones, and the interned order groups them. Each column is still
    // computed whole by one worker, and every counter is an integer sum or a maximum, so the assignment changes
    // the schedule and nothing else.
    const bool interleave = in.origin_reuse != nullptr;
    std::vector<FragmentScoringCounters> local(threads);
    std::vector<ScorerMemoryAccounts> local_memory(threads);
    std::vector<EmissionCacheCounters> local_cache(threads);
    std::vector<std::exception_ptr> failures(threads);
    std::vector<std::thread> workers;
    workers.reserve(threads);
    for (std::size_t t = 0; t < threads; ++t) {
        const std::size_t lo = interleave ? std::min(first_h + t, last_h) : std::min(first_h + t * span, last_h);
        const std::size_t hi = interleave ? last_h : std::min(lo + span, last_h);
        const std::size_t stride = interleave ? threads : 1;
        workers.emplace_back([&, t, lo, hi, stride]() {
            const ActiveWorker active;
            try {
                place_haplotype_range(in, lo, hi, stride, budget, mass_column, local[t], local_memory[t],
                                      local_cache[t], origins_out, census_out);
            } catch (const PlacementCancelled&) {
            } catch (...) {
                failures[t] = std::current_exception();
            }
        });
    }
    for (std::thread& worker : workers) worker.join();
    for (const std::exception_ptr& failure : failures) {
        if (!failure) continue;
        std::rethrow_exception(failure);
    }
    for (const EmissionCacheCounters& part : local_cache) {
        cache_work.executed_prefilters += part.executed_prefilters;
        cache_work.dp_cache_requests += part.dp_cache_requests;
        cache_work.dp_cache_hits += part.dp_cache_hits;
        cache_work.dp_cache_misses += part.dp_cache_misses;
        cache_work.executed_dps += part.executed_dps;
    }
    for (const ScorerMemoryAccounts& part : local_memory) {
        memory.placement_rows_peak = std::max(memory.placement_rows_peak, part.placement_rows_peak);
        memory.placement_capacity_bytes_peak =
            std::max(memory.placement_capacity_bytes_peak, part.placement_capacity_bytes_peak);
        memory.proposed_starts_peak =
            std::max(memory.proposed_starts_peak, part.proposed_starts_peak);
        memory.proposed_starts_capacity_bytes_peak = std::max(
            memory.proposed_starts_capacity_bytes_peak, part.proposed_starts_capacity_bytes_peak);
    }
    for (const FragmentScoringCounters& part : local) {
        counters.candidate_bases_scanned += part.candidate_bases_scanned;
        counters.piece_hits += part.piece_hits;
        counters.proposed_starts += part.proposed_starts;
        counters.unique_starts += part.unique_starts;
        counters.full_read_verifications += part.full_read_verifications;
        counters.dp_evaluations += part.dp_evaluations;
        counters.prefilter_rejections += part.prefilter_rejections;
        counters.accepted_mate_placements += part.accepted_mate_placements;
        counters.valid_fr_states += part.valid_fr_states;
        counters.latent_alignment_pairs += part.latent_alignment_pairs;
        counters.physical_origins += part.physical_origins;
        counters.statistical_state_classes += part.statistical_state_classes;
        counters.reused_physical_origins += part.reused_physical_origins;
        counters.reused_state_aggregates += part.reused_state_aggregates;
        counters.reused_latent_alignment_pairs += part.reused_latent_alignment_pairs;
        counters.reused_valid_fr_states += part.reused_valid_fr_states;
        counters.merged_state_classes += part.merged_state_classes;
    }
    return threads;
}

// ---- a repair, held as segments rather than as 226kb of bytes ----------------------------------
//
// A neighbour differs from its incumbent in one block, so materialising it means copying a whole
// haplotype, reversing it, and -- worst of all -- using it as a std::map key, where every comparison
// memcmps up to 226kb. Held as prefix + replacement + suffix, none of that happens: length, hash and
// equality all stream over the three segments, and only the window a delta actually reads is ever
// turned into bytes.
//
// `base` must outlive the segment. That is guaranteed by construction: a segment is only ever taken
// against a materialised INCUMBENT, and an incumbent is retained for the whole round in which its
// neighbours exist.
struct SegmentedSequence {
    const std::string* base = nullptr;
    std::size_t block_begin = 0;
    std::size_t block_end = 0;
    std::string replacement;

    std::size_t replaced_at_end() const { return block_begin + replacement.size(); }

    std::size_t length() const {
        return base->size() - (block_end - block_begin) + replacement.size();
    }

    char at(std::size_t i) const {
        if (i < block_begin) return (*base)[i];
        const std::size_t after = replaced_at_end();
        if (i < after) return replacement[i - block_begin];
        return (*base)[i - after + block_end];
    }

    // Appends [from, to) by SEGMENT, so each piece is one memcpy rather than a per-byte walk.
    void append_range(std::size_t from, std::size_t to, std::string& out) const {
        if (to <= from) return;
        const std::size_t after = replaced_at_end();
        // prefix
        if (from < block_begin) {
            const std::size_t stop = std::min(to, block_begin);
            out.append(base->data() + from, stop - from);
            from = stop;
        }
        // replacement
        if (from < after && from >= block_begin) {
            const std::size_t stop = std::min(to, after);
            out.append(replacement.data() + (from - block_begin), stop - from);
            from = stop;
        }
        // suffix, mapped back into the base's coordinates
        if (from < to) {
            const std::size_t base_from = from - after + block_end;
            out.append(base->data() + base_from, to - from);
        }
    }

    void append_range_reversed(std::size_t from, std::size_t to, std::string& out) const {
        if (to <= from) return;
        const std::size_t start = out.size();
        append_range(from, to, out);
        std::reverse(out.begin() + static_cast<std::ptrdiff_t>(start), out.end());
    }

    // FNV-1a streamed over the three segments. Selects a dedup bucket only: equality below decides.
    std::uint64_t hash() const {
        std::uint64_t h = 1469598103934665603ULL;
        const auto mix = [&h](const char* data, std::size_t count) {
            for (std::size_t i = 0; i < count; ++i) {
                h = (h ^ static_cast<unsigned char>(data[i])) * 1099511628211ULL;
            }
        };
        mix(base->data(), block_begin);
        mix(replacement.data(), replacement.size());
        mix(base->data() + block_end, base->size() - block_end);
        return h;
    }

    // Exact, after a hash match. The common case -- the same block of the same base replaced by the
    // same allele -- is settled without looking at a single base; anything else is compared for real,
    // because two different allele paths may spell the identical haplotype and must share one column.
    bool equals(const SegmentedSequence& other) const {
        if (base == other.base && block_begin == other.block_begin &&
            block_end == other.block_end && replacement == other.replacement) {
            return true;
        }
        const std::size_t n = length();
        if (n != other.length()) return false;
        constexpr std::size_t kChunk = 4096;
        std::string mine;
        std::string theirs;
        for (std::size_t at = 0; at < n; at += kChunk) {
            const std::size_t stop = std::min(n, at + kChunk);
            mine.clear();
            theirs.clear();
            append_range(at, stop, mine);
            other.append_range(at, stop, theirs);
            if (mine != theirs) return false;
        }
        return true;
    }

    bool equals_materialised(const std::string& other) const {
        const std::size_t n = length();
        if (n != other.size()) return false;
        constexpr std::size_t kChunk = 4096;
        std::string mine;
        for (std::size_t at = 0; at < n; at += kChunk) {
            const std::size_t stop = std::min(n, at + kChunk);
            mine.clear();
            append_range(at, stop, mine);
            if (std::memcmp(mine.data(), other.data() + at, mine.size()) != 0) return false;
        }
        return true;
    }

    void materialise(std::string& out) const {
        out.clear();
        out.reserve(length());
        append_range(0, length(), out);
    }
};

// Streamed hash of an already-materialised sequence, so the two representations share one bucket
// space. Must agree with SegmentedSequence::hash byte for byte.
std::uint64_t hash_materialised(const std::string& bases) {
    std::uint64_t h = 1469598103934665603ULL;
    for (const char c : bases) h = (h ^ static_cast<unsigned char>(c)) * 1099511628211ULL;
    return h;
}

// ---- the local window a replacement can disturb ------------------------------------------------
//
// Transient by construction: these bytes live only for one delta. That is safe because the emission
// cache copies whatever it retains into its own arena -- nothing it keeps points into `bases` below.
struct LocalWindow {
    std::string bases;
    std::string reversed;
    std::size_t begin = 0;     // offset of `bases` within the repaired haplotype
    bool at_start = false;     // the window's left edge IS the haplotype's start
    bool at_end = false;       // and likewise its right edge
};

// Enumerates every origin the window can support, in the repaired haplotype's own coordinates.
//
// The one trap here is the normaliser. gapped_emission divides by log_gap_shape_mass(read_length,
// haplotype.size() - start, ...), so a placement scored against a WINDOW is normalised against the
// window's remaining length, not the haplotype's. Where the window edge is synthetic those two
// differ, and the emission would silently disagree with the full scorer. The normaliser saturates at
// read_length + max_gap_bp, so a placement that keeps that much window on the relevant side is
// identical either way, and one that does not is discarded. Nothing is lost by discarding it: the
// halo is a whole maximum insert wider than the replacement, so no origin that touches the
// replacement can reach within a read length of a synthetic edge.
// THE REFERENCE. Kept verbatim so the optimised enumeration below can be certified against it, the
// way gap_recursion_logspace certifies the scaled recursion. It allocates a vector per read form on
// every call, which is exactly what the replacement exists to stop doing, so it is never on a
// production path -- only the self-check drives it.
void enumerate_window_origins_reference(
    const PreparedFragmentSet& prepared,
    const std::vector<ReadPair>& fragments,
    const LibraryModel& library,
    const FragmentScoringOptions& options,
    std::size_t band,
    const LocalWindow& window,
    PlacementBudget& budget,
    FragmentScoringCounters& counters,
    EmissionCacheCounters& cache_work,
    std::vector<std::vector<HaploidOrigin>>& out) {
    const auto form_slot = [](std::size_t fragment, std::size_t arm) { return fragment * 2 + arm; };
    std::vector<std::vector<FormPlacement>> forward_placements(prepared.forward_forms.size());
    std::vector<std::vector<FormPlacement>> reverse_placements(prepared.reverse_forms.size());
    std::vector<std::pair<std::uint32_t, std::uint32_t>> proposed_scratch;

    collect_form_placements(window.bases, prepared.forward_forms, prepared.forward_componentise,
                            prepared.forward_indexes, budget, band, library, forward_placements,
                            proposed_scratch, counters, cache_work, options);
    collect_form_placements(window.reversed, prepared.reverse_forms, prepared.reverse_componentise,
                            prepared.reverse_indexes, budget, band, library, reverse_placements,
                            proposed_scratch, counters, cache_work, options);

    const std::size_t width = window.bases.size();
    // Forward rows: drop anything that would have been normalised against a synthetic right edge,
    // then shift into the repaired haplotype's coordinates.
    for (std::size_t i = 0; i < forward_placements.size(); ++i) {
        const std::size_t saturating = prepared.forward_forms[i].bases.size() + library.max_gap_bp;
        std::vector<FormPlacement> kept;
        kept.reserve(forward_placements[i].size());
        for (const FormPlacement& placement : forward_placements[i]) {
            if (!window.at_end && width - placement.start < saturating) continue;
            FormPlacement moved = placement;
            moved.start = static_cast<std::uint32_t>(window.begin + placement.start);
            kept.push_back(moved);
        }
        forward_placements[i] = std::move(kept);
    }
    // Reverse rows carry a reversed-coordinate start; map to a right endpoint exactly as the full
    // scorer does, drop anything normalised against a synthetic LEFT edge, then shift and re-sort.
    for (std::size_t i = 0; i < reverse_placements.size(); ++i) {
        const std::size_t saturating = prepared.reverse_forms[i].bases.size() + library.max_gap_bp;
        std::vector<FormPlacement> kept;
        kept.reserve(reverse_placements[i].size());
        for (const FormPlacement& placement : reverse_placements[i]) {
            const std::size_t endpoint = width - placement.start;
            if (!window.at_start && endpoint < saturating) continue;
            FormPlacement moved = placement;
            moved.start = static_cast<std::uint32_t>(window.begin + endpoint);
            kept.push_back(moved);
        }
        std::sort(kept.begin(), kept.end(), [](const FormPlacement& lhs, const FormPlacement& rhs) {
            return std::tie(lhs.start, lhs.consumed) < std::tie(rhs.start, rhs.consumed);
        });
        reverse_placements[i] = std::move(kept);
    }

    // The pairing loop, identical to the full scorer's: one group per left endpoint, one origin per
    // reachable right endpoint, latent alignments summed into the origin rather than counted as
    // extra copies.
    out.assign(fragments.size(), {});
    std::map<std::uint32_t, OriginAccumulator> origins;
    for (std::size_t f = 0; f < fragments.size(); ++f) {
        const InsertSupport& support = prepared.support_by_fragment[f];
        if (support.empty()) continue;
        for (std::size_t arm = 0; arm < 2; ++arm) {
            const std::vector<FormPlacement>& forward = forward_placements[form_slot(f, arm)];
            const std::vector<FormPlacement>& reverse = reverse_placements[form_slot(f, arm)];
            if (forward.empty() || reverse.empty()) continue;
            for (std::size_t first_row = 0; first_row < forward.size();) {
                const std::uint32_t forward_start = forward[first_row].start;
                std::size_t last_row = first_row;
                while (last_row < forward.size() && forward[last_row].start == forward_start) ++last_row;
                const std::size_t lo_end = static_cast<std::size_t>(forward_start) + support.min_bp;
                const std::size_t hi_end = static_cast<std::size_t>(forward_start) + support.max_bp;
                const auto begin = std::lower_bound(
                    reverse.begin(), reverse.end(), lo_end,
                    [](const FormPlacement& lhs, std::size_t value) { return lhs.start < value; });
                const auto end = std::upper_bound(
                    begin, reverse.end(), hi_end,
                    [](std::size_t value, const FormPlacement& rhs) { return value < rhs.start; });
                origins.clear();
                for (std::size_t row = first_row; row < last_row; ++row) {
                    for (auto it = begin; it != end; ++it) {
                        ++counters.latent_alignment_pairs;
                        const std::optional<std::size_t> insert = fr_insert_from_endpoints(
                            forward_start, forward[row].consumed, it->start, it->consumed, support);
                        if (!insert.has_value()) continue;
                        ++counters.valid_fr_states;
                        origins[static_cast<std::uint32_t>(*insert)].add(
                            forward[row].emission + it->emission, it->start - it->consumed,
                            forward[row].consumed, it->consumed, forward[row].mismatches,
                            it->mismatches);
                    }
                }
                for (const auto& item : origins) {
                    ++counters.physical_origins;
                    const std::uint64_t seen =
                        budget.physical_origins.fetch_add(1, std::memory_order_relaxed) + 1;
                    if (seen > budget.max_physical_origins) {
                        budget.cancel();
                        throw ResourceExhausted(
                            "fragment scoring reached its physical-origin cap; no partial ranking is "
                            "returned");
                    }
                    out[f].push_back({forward_start, item.first, static_cast<std::uint32_t>(arm),
                                      item.second.log_latent_mass});
                }
                first_row = last_row;
            }
        }
    }
}

// I2-B TEST ONLY: an exact per-window comparison of the fine rows with the reference's own two scans, and deliberate
// ownership/context mutations (1: overlapping ownership, core_hi + 1; 2: left context one byte short; 3: right context
// two bytes short; 4: the left end window two bytes short) that the comparison must catch.
std::atomic<bool>& fine_window_row_check() {
    static std::atomic<bool> on{false};
    return on;
}
std::atomic<int>& fine_window_mutation() {
    static std::atomic<int> m{0};
    return m;
}

// I2-B: a window's raw rows -- exactly what the reference's two collect_form_placements calls write into scratch's
// forward rows (window coordinates) and reverse rows (reversed-window coordinates), each in (start, consumed) order, with
// the same touched lists -- built from the worker's round-local fine store. Executed work (the discovery of units not
// held) is counted in `counters`/`cache_work` and charged to the shared budget by collect_form_placements itself;
// served placements are tallied apart and charged to the cumulative served cap. `proposed_by_pass` receives the
// executed proposals of the forward and the reverse pieces. A closed store that does not hold every piece of the window
// hands it to the reference path (recorded); the result is the same either way.
void fine_window_rows(const PreparedFragmentSet& prepared, const LibraryModel& library, const FragmentScoringOptions& options,
                      std::size_t band, const LocalWindow& window, const std::vector<PieceIndex>& forward_indexes,
                      const std::vector<PieceIndex>& reverse_indexes, PlacementBudget& budget,
                      FragmentScoringCounters& counters, EmissionCacheCounters& cache_work, WindowScratch& scratch,
                      std::uint64_t* proposed_by_pass) {
    WindowFineStore& st = scratch.fine;
    ++st.windows;
    // the window as one text with one segment, per strand: its end windows (cuts 0 and n) and, when it can keep an
    // interior start at all (n > 2g, as add_haplotype requires), its fine units
    // Each start has exactly one owner: a start whose span crosses BOTH ends (a window shorter than a read's span) is
    // the left end window's, so the right one, when it reaches the text's start, skips starts below g.
    struct Candidate {
        IndexUnit unit;
        std::uint32_t offset;
        bool skip_left_crossing;
    };
    std::vector<Candidate> candidates;
    const int mutation = fine_window_mutation().load(std::memory_order_relaxed);
    for (std::uint8_t strand = 0; strand < 2; ++strand) {
        const std::string& text = strand == 0 ? window.bases : window.reversed;
        const std::size_t n = text.size();
        for (const std::size_t b : {std::size_t{0}, n}) {
            if (b == n && n == 0) break;
            std::size_t lo = 0;
            IndexUnit u = boundary_window_unit(text, strand, b, st.window, lo);
            if (mutation == 4 && b == 0 && n > 0 && u.bytes.size() > 2) u.bytes.resize(u.bytes.size() - 2);
            candidates.push_back({std::move(u), static_cast<std::uint32_t>(lo), b == n && n > 0 && lo == 0});
        }
        if (n > 2 * st.gap) {
            for_each_fine_piece(text, strand, st.gap, st.max_read, [&](IndexUnit&& u, std::uint32_t offset) {
                if (mutation == 1 && u.core_hi < u.bytes.size()) ++u.core_hi;
                if (mutation == 2 && offset > 0 && u.core_lo > 0) {
                    u.bytes.erase(0, 1);
                    --u.core_lo;
                    --u.core_hi;
                    ++offset;
                }
                if (mutation == 3 && offset + u.bytes.size() < n && u.bytes.size() >= u.core_hi + 2) {
                    u.bytes.resize(u.bytes.size() - 2);
                }
                candidates.push_back({std::move(u), offset, false});
            });
        }
    }
    const auto find = [&](const IndexUnit& unit, std::uint64_t key) -> long {
        const auto range = st.by_hash.equal_range(key);
        for (auto it = range.first; it != range.second; ++it) {
            if (same_unit(st.units[it->second], unit)) return static_cast<long>(it->second);
            ++st.hash_collisions;
        }
        return -1;
    };
    std::vector<std::uint64_t> keys(candidates.size());
    for (std::size_t i = 0; i < candidates.size(); ++i) keys[i] = unit_hash(candidates[i].unit);
    if (st.closed) {
        bool all_held = true;
        for (std::size_t i = 0; i < candidates.size() && all_held; ++i) all_held = find(candidates[i].unit, keys[i]) >= 0;
        if (!all_held) {
            // the reference path, unchanged: the two scans over the whole window
            ++st.fallback_windows;
            const std::uint64_t at_start = counters.proposed_starts;
            collect_form_placements(window.bases, prepared.forward_forms, prepared.forward_componentise, forward_indexes,
                                    budget, band, library, scratch.forward_rows, scratch.proposed, counters, cache_work,
                                    options, &scratch.touched_forward);
            const std::uint64_t after_forward = counters.proposed_starts;
            collect_form_placements(window.reversed, prepared.reverse_forms, prepared.reverse_componentise, reverse_indexes,
                                    budget, band, library, scratch.reverse_rows, scratch.proposed, counters, cache_work,
                                    options, &scratch.touched_reverse);
            if (proposed_by_pass != nullptr) {
                proposed_by_pass[0] = after_forward - at_start;
                proposed_by_pass[1] = counters.proposed_starts - after_forward;
            }
            return;
        }
    }
    st.pieces.clear();
    st.transient.clear();
    std::uint64_t proposed[2] = {0, 0}, served = 0;
    for (std::size_t i = 0; i < candidates.size(); ++i) {
        Candidate& c = candidates[i];
        const std::uint8_t strand = c.unit.strand;
        const long held = find(c.unit, keys[i]);
        if (held >= 0) {
            const std::vector<StoredPlacement>& rows = st.placements[static_cast<std::size_t>(held)];
            st.pieces.push_back({&rows, c.offset, strand, c.skip_left_crossing});
            ++st.units_served;
            served += rows.size();
            continue;
        }
        std::vector<StoredPlacement> found;
        const std::uint64_t before = counters.proposed_starts;
        discover_unit_placements(c.unit, prepared, forward_indexes, reverse_indexes, budget, band, library, options, st.rows,
                                 st.proposed, st.touched, counters, cache_work, nullptr, nullptr, nullptr, found);
        proposed[strand] += counters.proposed_starts - before;
        st.placements_discovered += found.size();
        std::vector<StoredPlacement> exact(found.begin(), found.end());   // exact capacity: no growth slack is held
        const std::uint64_t cost = WindowFineStore::entry_bytes(c.unit, exact);
        const std::uint64_t buckets_after =
            std::max<std::uint64_t>(st.by_hash.bucket_count(), (st.by_hash.size() + 1) * 2) * sizeof(void*);
        if (!st.closed && st.bytes + cost + buckets_after <= st.budget) {
            const std::uint32_t id = static_cast<std::uint32_t>(st.units.size());
            st.units.push_back(std::move(c.unit));
            st.placements.push_back(std::move(exact));
            st.by_hash.emplace(keys[i], id);
            st.bytes += cost;
            ++st.units_discovered;
            st.pieces.push_back({&st.placements.back(), c.offset, strand, c.skip_left_crossing});
        } else {
            st.closed = true;   // later windows not fully held take the reference path
            ++st.units_not_retained;
            st.transient.push_back(std::move(exact));
            st.pieces.push_back({&st.transient.back(), c.offset, strand, c.skip_left_crossing});
        }
    }
    if (proposed_by_pass != nullptr) {
        proposed_by_pass[0] = proposed[0];
        proposed_by_pass[1] = proposed[1];
    }
    if (served != 0) {
        st.placements_served += served;
        if (st.served_total != nullptr) {
            const std::uint64_t total = st.served_total->fetch_add(served, std::memory_order_relaxed) + served;
            if (st.served_cap != 0 && total > st.served_cap) {
                budget.cancel();
                throw ResourceExhausted("the adaptive fine store reached its served-placement cap; no partial ranking is "
                                        "returned");
            }
        }
    }
    // reconstruction: every piece's owned placements, shifted to the window, into the raw rows; then the reference's row
    // order. No deduplication: ownership is exclusive, so a duplicate would be a defect for the row check to catch.
    const std::uint32_t g32 = static_cast<std::uint32_t>(st.gap);
    for (const WindowFineStore::Piece& piece : st.pieces) {
        auto& rows = piece.strand == 0 ? scratch.forward_rows : scratch.reverse_rows;
        auto& touched = piece.strand == 0 ? scratch.touched_forward : scratch.touched_reverse;
        for (const StoredPlacement& p : *piece.rows) {
            if (piece.skip_left_crossing && p.start < g32) continue;   // the left end window owns it
            std::vector<FormPlacement>& row = rows[p.form];
            if (row.empty()) touched.push_back(p.form);
            FormPlacement q;
            q.start = piece.offset + p.start;
            q.consumed = p.consumed;
            q.mismatches = p.mismatches;
            q.emission = p.emission;
            row.push_back(q);
        }
    }
    for (int strand = 0; strand < 2; ++strand) {
        auto& rows = strand == 0 ? scratch.forward_rows : scratch.reverse_rows;
        for (const std::uint32_t form : strand == 0 ? scratch.touched_forward : scratch.touched_reverse) {
            std::vector<FormPlacement>& row = rows[form];
            std::sort(row.begin(), row.end(), [](const FormPlacement& a, const FormPlacement& b) {
                return std::tie(a.start, a.consumed) < std::tie(b.start, b.consumed);
            });
        }
    }
    {
        // I2-C: temporary storage of this window: the reusable scratch plus the candidate units' byte copies and keys
        std::uint64_t candidate_bytes = candidates.capacity() * sizeof(Candidate) + keys.capacity() * sizeof(std::uint64_t);
        for (const Candidate& c : candidates) candidate_bytes += c.unit.bytes.capacity();
        st.peak_temporary = std::max(st.peak_temporary, st.scratch_bytes() + candidate_bytes);
        st.peak_bytes = std::max(st.peak_bytes, st.held_bytes() + st.scratch_bytes() + candidate_bytes);
    }
    st.peak_held = std::max(st.peak_held, st.held_bytes());
    st.transient.clear();
    if (fine_window_row_check().load(std::memory_order_relaxed)) {
        // TEST ONLY: the reference's two scans over the whole window, against a throwaway budget, row for row
        PlacementBudget spare;
        spare.max_piece_hits = std::numeric_limits<std::uint64_t>::max();
        spare.max_physical_origins = std::numeric_limits<std::uint64_t>::max();
        FragmentScoringCounters spare_counters;
        EmissionCacheCounters spare_cache;
        std::vector<std::vector<FormPlacement>> ref[2] = {
            std::vector<std::vector<FormPlacement>>(prepared.forward_forms.size()),
            std::vector<std::vector<FormPlacement>>(prepared.reverse_forms.size())};
        std::vector<std::pair<std::uint32_t, std::uint32_t>> spare_proposed;
        collect_form_placements(window.bases, prepared.forward_forms, prepared.forward_componentise, forward_indexes, spare,
                                band, library, ref[0], spare_proposed, spare_counters, spare_cache, options);
        collect_form_placements(window.reversed, prepared.reverse_forms, prepared.reverse_componentise, reverse_indexes,
                                spare, band, library, ref[1], spare_proposed, spare_counters, spare_cache, options);
        for (int strand = 0; strand < 2; ++strand) {
            const auto& got = strand == 0 ? scratch.forward_rows : scratch.reverse_rows;
            for (std::size_t form = 0; form < ref[strand].size(); ++form) {
                const auto& a = ref[strand][form];
                const auto& b = got[form];
                bool same = a.size() == b.size();
                for (std::size_t i = 0; same && i < a.size(); ++i) {
                    same = a[i].start == b[i].start && a[i].consumed == b[i].consumed && a[i].mismatches == b[i].mismatches &&
                           std::memcmp(&a[i].emission, &b[i].emission, sizeof(double)) == 0;
                }
                if (!same) {
                    throw std::logic_error("fine window rows differ from the reference at strand " + std::to_string(strand) +
                                           " form " + std::to_string(form) + " (" + std::to_string(b.size()) + " vs " +
                                           std::to_string(a.size()) + " rows)");
                }
            }
        }
    }
}

// The optimised enumeration. Same placements, same origins, same counters, same order -- the only
// things that change are who owns the memory and how the per-endpoint sums are addressed.
//
// Two substitutions, both exact:
//
//   persistent rows   the reference built a vector per read form per candidate. These persist in the
//                     caller's scratch and only the rows the previous window wrote are cleared.
//
//   endpoint groups   the reference kept a std::map<insert, OriginAccumulator> per forward start.
//                     fr_insert_from_endpoints returns exactly reverse_end - forward_start when it
//                     accepts, so for a fixed forward start the insert is strictly increasing in the
//                     endpoint: contiguous equal-endpoint runs of the reverse range ARE the map's
//                     keys, already in ascending order. The loop nesting is unchanged -- forward row
//                     outer, reverse row inner, groups ascending -- so every accumulator receives its
//                     additions in the sequence the map's accumulator received them, and log_add_exp
//                     is applied to the same values in the same order. That is what keeps it
//                     bit-identical rather than merely equal in exact arithmetic.
//
// The pair COUNTER is order-sensitive too: latent_alignment_pairs counts every pair tested, rejected
// ones included, and the physical-origin budget is charged per emitted origin in ascending insert
// order. Reaching the same origins by a different route would keep the masses and move both.
void enumerate_window_origins(
    const PreparedFragmentSet& prepared,
    const std::vector<ReadPair>& fragments,
    const LibraryModel& library,
    const FragmentScoringOptions& options,
    std::size_t band,
    const LocalWindow& window,
    PlacementBudget& budget,
    FragmentScoringCounters& counters,
    EmissionCacheCounters& cache_work,
    WindowScratch& scratch,
    // MEASUREMENT ONLY (the fragment-range diagnostic): piece indexes restricted to one range's forms,
    // and that range. The defaults are the production path, unchanged.
    const std::vector<PieceIndex>* forward_indexes = nullptr,
    const std::vector<PieceIndex>* reverse_indexes = nullptr,
    std::size_t fragment_begin = 0,
    std::size_t fragment_end = std::numeric_limits<std::size_t>::max(),
    // Proposed starts of the forward and the reverse pass, for the exact post-merge cap check.
    std::uint64_t* proposed_by_pass = nullptr) {
    const auto form_slot = [](std::size_t fragment, std::size_t arm) { return fragment * 2 + arm; };
    scratch.ensure_rows(prepared.forward_forms.size(), prepared.reverse_forms.size());
    scratch.clear_touched();

    const auto rows_mark = std::chrono::steady_clock::now();
    ++scratch.windows_enumerated;
    if (scratch.fine.active) {
        // I2-B (experimental): the same raw rows, from the worker's round-local fine store
        fine_window_rows(prepared, library, options, band, window,
                         forward_indexes != nullptr ? *forward_indexes : prepared.forward_indexes,
                         reverse_indexes != nullptr ? *reverse_indexes : prepared.reverse_indexes, budget, counters,
                         cache_work, scratch, proposed_by_pass);
    } else {
    const std::uint64_t proposed_at_start = counters.proposed_starts;
    collect_form_placements(window.bases, prepared.forward_forms, prepared.forward_componentise,
                            forward_indexes != nullptr ? *forward_indexes : prepared.forward_indexes,
                            budget, band, library, scratch.forward_rows,
                            scratch.proposed, counters, cache_work, options, &scratch.touched_forward);
    const std::uint64_t proposed_after_forward = counters.proposed_starts;
    collect_form_placements(window.reversed, prepared.reverse_forms, prepared.reverse_componentise,
                            reverse_indexes != nullptr ? *reverse_indexes : prepared.reverse_indexes,
                            budget, band, library, scratch.reverse_rows,
                            scratch.proposed, counters, cache_work, options, &scratch.touched_reverse);
    if (proposed_by_pass != nullptr) {
        proposed_by_pass[0] = proposed_after_forward - proposed_at_start;
        proposed_by_pass[1] = counters.proposed_starts - proposed_after_forward;
    }
    }
    const auto pairing_mark = std::chrono::steady_clock::now();
    scratch.seconds_rows += std::chrono::duration<double>(pairing_mark - rows_mark).count();
    const std::size_t last_fragment = std::min(fragment_end, fragments.size());

    const std::size_t width = window.bases.size();
    // Filtered and shifted IN PLACE. resize() never grows here, so no row reallocates.
    for (const std::uint32_t i : scratch.touched_forward) {
        std::vector<FormPlacement>& row = scratch.forward_rows[i];
        const std::size_t saturating = prepared.forward_forms[i].bases.size() + library.max_gap_bp;
        std::size_t write = 0;
        for (std::size_t read = 0; read < row.size(); ++read) {
            if (!window.at_end && width - row[read].start < saturating) continue;
            FormPlacement moved = row[read];
            moved.start = static_cast<std::uint32_t>(window.begin + row[read].start);
            row[write++] = moved;
        }
        row.resize(write);
    }
    for (const std::uint32_t i : scratch.touched_reverse) {
        std::vector<FormPlacement>& row = scratch.reverse_rows[i];
        const std::size_t saturating = prepared.reverse_forms[i].bases.size() + library.max_gap_bp;
        std::size_t write = 0;
        for (std::size_t read = 0; read < row.size(); ++read) {
            const std::size_t endpoint = width - row[read].start;
            if (!window.at_start && endpoint < saturating) continue;
            FormPlacement moved = row[read];
            moved.start = static_cast<std::uint32_t>(window.begin + endpoint);
            row[write++] = moved;
        }
        row.resize(write);
        std::sort(row.begin(), row.end(), [](const FormPlacement& lhs, const FormPlacement& rhs) {
            return std::tie(lhs.start, lhs.consumed) < std::tie(rhs.start, rhs.consumed);
        });
    }

    scratch.origins.clear();
    scratch.offsets.assign(fragments.size() + 1, 0);
    for (std::size_t f = fragment_begin; f < last_fragment; ++f) {
        scratch.offsets[f] = scratch.origins.size();
        const InsertSupport& support = prepared.support_by_fragment[f];
        if (support.empty()) continue;
        for (std::size_t arm = 0; arm < 2; ++arm) {
            const std::vector<FormPlacement>& forward = scratch.forward_rows[form_slot(f, arm)];
            const std::vector<FormPlacement>& reverse = scratch.reverse_rows[form_slot(f, arm)];
            if (forward.empty() || reverse.empty()) continue;
            for (std::size_t first_row = 0; first_row < forward.size();) {
                const std::uint32_t forward_start = forward[first_row].start;
                std::size_t last_row = first_row;
                while (last_row < forward.size() && forward[last_row].start == forward_start) ++last_row;
                const std::size_t lo_end = static_cast<std::size_t>(forward_start) + support.min_bp;
                const std::size_t hi_end = static_cast<std::size_t>(forward_start) + support.max_bp;
                const auto begin = std::lower_bound(
                    reverse.begin(), reverse.end(), lo_end,
                    [](const FormPlacement& lhs, std::size_t value) { return lhs.start < value; });
                const auto end = std::upper_bound(
                    begin, reverse.end(), hi_end,
                    [](std::size_t value, const FormPlacement& rhs) { return value < rhs.start; });

                // Contiguous equal-endpoint runs, ascending: the map's keys without the map.
                scratch.group_begin.clear();
                scratch.group_end.clear();
                for (auto it = begin; it != end;) {
                    auto stop = it;
                    while (stop != end && stop->start == it->start) ++stop;
                    scratch.group_begin.push_back(static_cast<std::size_t>(it - reverse.begin()));
                    scratch.group_end.push_back(static_cast<std::size_t>(stop - reverse.begin()));
                    it = stop;
                }
                const std::size_t groups = scratch.group_begin.size();
                if (groups == 0) {
                    first_row = last_row;
                    continue;
                }
                if (scratch.accumulators.size() < groups) {
                    scratch.accumulator_growths += groups - scratch.accumulators.size();
                    scratch.accumulators.resize(groups);
                    scratch.accepted.resize(groups);
                }
                for (std::size_t g = 0; g < groups; ++g) {
                    scratch.accumulators[g] = OriginAccumulator{};
                    scratch.accepted[g] = 0;
                }

                // Forward row outer, reverse row inner, groups ascending -- the reference's order.
                for (std::size_t row = first_row; row < last_row; ++row) {
                    for (std::size_t g = 0; g < groups; ++g) {
                        for (std::size_t at = scratch.group_begin[g]; at < scratch.group_end[g]; ++at) {
                            const FormPlacement& rev = reverse[at];
                            ++counters.latent_alignment_pairs;
                            const std::optional<std::size_t> insert = fr_insert_from_endpoints(
                                forward_start, forward[row].consumed, rev.start, rev.consumed, support);
                            if (!insert.has_value()) continue;
                            ++counters.valid_fr_states;
                            scratch.accumulators[g].add(forward[row].emission + rev.emission,
                                                        rev.start - rev.consumed, forward[row].consumed,
                                                        rev.consumed, forward[row].mismatches,
                                                        rev.mismatches);
                            scratch.accepted[g] = 1;
                        }
                    }
                }

                for (std::size_t g = 0; g < groups; ++g) {
                    if (scratch.accepted[g] == 0) continue;   // the map held no key for this endpoint
                    ++counters.physical_origins;
                    const std::uint64_t seen =
                        budget.physical_origins.fetch_add(1, std::memory_order_relaxed) + 1;
                    if (seen > budget.max_physical_origins) {
                        budget.cancel();
                        throw ResourceExhausted(
                            "fragment scoring reached its physical-origin cap; no partial ranking is "
                            "returned");
                    }
                    const std::uint32_t insert =
                        reverse[scratch.group_begin[g]].start - forward_start;
                    if (scratch.origins.size() == scratch.origins.capacity()) ++scratch.origin_growths;
                    scratch.origins.push_back({forward_start, insert, static_cast<std::uint32_t>(arm),
                                               scratch.accumulators[g].log_latent_mass});
                }
                first_row = last_row;
            }
        }
    }
    scratch.offsets[last_fragment] = scratch.origins.size();
    scratch.seconds_pairing += std::chrono::duration<double>(std::chrono::steady_clock::now() - pairing_mark).count();
}

// The four acceptance caps, decided on whatever totals the caller holds them against: per call for
// the one-shot path, cumulative for a context. Throwing here is the same refusal either way.
void enforce_acceptance_caps(
    const FragmentScoringCounters& counters,
    const FragmentScoringOptions& options,
    const std::function<void()>& note) {
    if (counters.piece_hits > options.max_piece_hits) {
        note();
        throw ResourceExhausted(
            "fragment scoring reached its piece-hit cap; the locus is too repetitive for the declared "
            "caps and no partial ranking is returned");
    }
    if (counters.physical_origins > options.max_physical_origins) {
        note();
        throw ResourceExhausted(
            "fragment scoring reached its physical-origin cap; no partial ranking is returned");
    }
    if (counters.statistical_state_classes > options.max_state_classes) {
        note();
        throw ResourceExhausted(
            "fragment scoring reached its statistical-class cap; no partial ranking is returned");
    }
    if (counters.latent_alignment_pairs > options.max_latent_alignment_pairs) {
        note();
        throw ResourceExhausted(
            "fragment scoring reached its latent-alignment-pair cap; the declared gap width expands "
            "the pairing work beyond the cap and no partial ranking is returned");
    }
}

}  // namespace scoring_detail

// The exact union prefilter-decision cache (UnionPrefilterCacheReport), built once per read set. Shared by
// score_diplotypes and the EXPERIMENTAL content-index batch (I1.1) so the two cannot configure it differently. Needs one
// read length (throws otherwise; single_read_length() lets a caller fall back instead) and a positive consumption cap.
bool single_read_length(const PreparedFragmentSet& prepared, std::size_t& length) {
    length = 0;
    for (const auto* forms : {&prepared.forward_forms, &prepared.reverse_forms}) {
        for (const ReadForm& form : *forms) {
            if (length == 0) length = form.bases.size();
            if (form.bases.size() != length) return false;
        }
    }
    return true;
}
std::unique_ptr<UnionPrefilterCache> make_union_prefilter_cache(const PreparedFragmentSet& prepared, const LibraryModel& library,
                                                                const FragmentScoringOptions& options,
                                                                UnionPrefilterCacheReport& report) {
    if (options.max_union_cached_prefilters_consumed == 0) {
        throw std::invalid_argument("union prefilter reuse needs a positive consumption cap");
    }
    auto prefilter_cache = std::make_unique<UnionPrefilterCache>();
    std::size_t r = 0;
    if (!single_read_length(prepared, r)) throw std::invalid_argument("union prefilter reuse needs one read length");
    prefilter_cache->read_length = r;
    prefilter_cache->window_length = r + library.max_gap_bp;
    std::uint64_t slots = 1024;
    while (slots * 2 * sizeof(std::uint64_t) <= options.union_prefilter_cache_bytes) slots *= 2;
    prefilter_cache->slot_bits = 0;
    while ((std::uint64_t{1} << prefilter_cache->slot_bits) < slots) ++prefilter_cache->slot_bits;
    prefilter_cache->slots.reset(new std::atomic<std::uint64_t>[slots]());
    prefilter_cache->max_held = slots / 2;
    prefilter_cache->max_consumed = options.max_union_cached_prefilters_consumed;
    prefilter_cache->publish_batch = std::max<std::uint64_t>(1, options.piece_hit_publish_batch);
    report.active = true;
    report.budget_bytes = options.union_prefilter_cache_bytes;
    report.max_consumed = options.max_union_cached_prefilters_consumed;
    report.read_length = r;
    report.window_length = r + library.max_gap_bp;
    report.decision_table_bytes = slots * sizeof(std::uint64_t);
    return prefilter_cache;
}
void fill_union_prefilter_report(UnionPrefilterCache& cache, UnionPrefilterCacheReport& report) {
    report.window_contents = cache.next_window.load(std::memory_order_relaxed);
    report.window_index_bytes = 0;
    for (UnionPrefilterCache::Shard& shard : cache.shards) {
        report.window_index_bytes += shard.bytes + shard.table.capacity() * sizeof(UnionPrefilterCache::Stored);
    }
    report.positions_indexed = cache.positions_indexed.load(std::memory_order_relaxed);
    report.decisions_held = cache.held.load(std::memory_order_relaxed);
    report.decisions_not_retained = cache.not_retained.load(std::memory_order_relaxed);
    report.hits = cache.consumed.load(std::memory_order_relaxed);
    report.misses = cache.misses.load(std::memory_order_relaxed);
    report.worker_seconds_index = static_cast<double>(cache.index_nanos.load(std::memory_order_relaxed)) * 1e-9;
}

// The diploid combination of placed mass columns into candidate scores, shared by score_diplotypes and the
// EXPERIMENTAL content index (I1) so the two cannot combine differently. Fills result.scores (already sized),
// result.audit_scores, threads_used_combination, candidate_fragment_evaluations and seconds_diploid_combination.
void combine_candidate_scores(const std::vector<std::vector<double>>& mass_column, const std::vector<double>& exposure,
                              const std::vector<double>& background, std::size_t fragment_count,
                              const LibraryModel& library, const std::vector<HaplotypeCandidate>& candidates,
                              std::size_t threads, const std::vector<LibraryModel>& audit_libraries,
                              FragmentScoringResult& result, bool prepared_library_terms = false) {
    const auto combination_mark = std::chrono::steady_clock::now();
    const std::vector<AuditCombination> audits = audit_combinations(audit_libraries);
    // I2-F: the run's library prepared once (validated by audit_combinations), or the scalar reference per term
    // (prepared only when some term will be combined, so validation happens exactly where the scalar path would first
    // validate: never for an empty candidate set)
    const std::vector<AuditCombination> main_terms = prepared_library_terms && !candidates.empty()
                                                         ? audit_combinations({library})
                                                         : std::vector<AuditCombination>{};
    const AuditCombination* main = main_terms.empty() ? nullptr : &main_terms.front();
    result.audit_scores.assign(audits.size(), std::vector<DiplotypeScore>(candidates.size()));
    const auto combine_range = [&](std::size_t begin, std::size_t end) {
        for (std::size_t c = begin; c < end && !audits.empty(); ++c) {
            const std::size_t a = candidates[c].first;
            const std::size_t b = candidates[c].second;
            combine_with_audits(mass_column[a], mass_column[b], background, fragment_count,
                                exposure[a] + exposure[b], library, audits, result.scores[c], result.audit_scores, c, main);
        }
        for (std::size_t c = begin; c < end && audits.empty(); ++c) {
            const std::size_t a = candidates[c].first;
            const std::size_t b = candidates[c].second;
            DiplotypeScore& score = result.scores[c];
            score.exposure_term = main != nullptr ? -(main->exposure_scale * (exposure[a] + exposure[b]))
                                                  : -expected_signal_fragments(exposure[a] + exposure[b], library);
            score.absolute_term_sum = std::abs(score.exposure_term);
            double fragment_term = 0.0;
            for (std::size_t f = 0; f < fragment_count; ++f) {
                const double combined = log_add_exp(mass_column[a][f], mass_column[b][f]);
                if (combined != kNegInf) ++score.fragments_with_mass;
                const double contribution = main != nullptr ? prepared_contribution(*main, combined, background[f])
                                                            : log_fragment_contribution(combined, background[f], library);
                fragment_term += contribution;
                score.absolute_term_sum += std::abs(contribution);
            }
            score.fragment_term = fragment_term;
            score.log_likelihood = score.exposure_term + fragment_term;
        }
    };

    const std::size_t combination_threads =
        candidates.empty() ? 1
                           : std::max<std::size_t>(1, std::min(threads, candidates.size()));
    result.threads_used_combination = combination_threads;
    result.candidate_fragment_evaluations =
        static_cast<std::uint64_t>(candidates.size()) * static_cast<std::uint64_t>(fragment_count);

    if (combination_threads <= 1) {
        const ActiveWorker active;
        combine_range(0, candidates.size());
    } else {
        const std::size_t span = (candidates.size() + combination_threads - 1) / combination_threads;
        std::vector<std::exception_ptr> failures(combination_threads);
        std::vector<std::thread> workers;
        workers.reserve(combination_threads);
        for (std::size_t t = 0; t < combination_threads; ++t) {
            const std::size_t begin = std::min(t * span, candidates.size());
            const std::size_t end = std::min(begin + span, candidates.size());
            workers.emplace_back([&, t, begin, end]() {
                const ActiveWorker active;
                try {
                    combine_range(begin, end);
                } catch (...) {
                    failures[t] = std::current_exception();
                }
            });
        }
        for (std::thread& worker : workers) worker.join();
        for (const std::exception_ptr& failure : failures) {
            if (failure) std::rethrow_exception(failure);
        }
    }
    if (main != nullptr) {
        std::uint64_t with_mass = 0;
        for (const DiplotypeScore& score : result.scores) with_mass += score.fragments_with_mass;
        note_prepared_combination(candidates.size(), fragment_count, with_mass, *main);
    }
    result.seconds_diploid_combination =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - combination_mark).count();
}

FragmentScoringResult score_diplotypes(
    const InternedHaplotypes& haplotypes,
    const std::vector<HaplotypeCandidate>& candidates,
    const std::vector<ReadPair>& fragments,
    const LibraryModel& library,
    const FragmentScoringOptions& options,
    FragmentScoringObservation* observed_on_failure) {
    validate_library_model(library);
    if (options.minimum_piece_bp == 0) {
        throw std::invalid_argument("minimum indexable piece length must be positive");
    }
    if (options.max_piece_hits == 0 || options.max_proposed_starts_per_haplotype == 0 ||
        options.max_physical_origins == 0 || options.max_state_classes == 0 ||
        options.max_mass_cells == 0 || options.max_latent_alignment_pairs == 0) {
        throw std::invalid_argument("every fragment-scoring resource cap must be positive");
    }
    for (const HaplotypeCandidate& candidate : candidates) {
        if (candidate.first >= haplotypes.sequences.size() ||
            candidate.second >= haplotypes.sequences.size()) {
            throw std::invalid_argument("diplotype candidate names an unknown haplotype");
        }
    }

    FragmentScoringResult result;
    result.scores.resize(candidates.size());
    result.unique_haplotypes = haplotypes.sequences.size();
    FragmentScoringCounters& counters = result.counters;
    const std::size_t band = options.mismatch_band;
    const DiscreteInsertDistribution library_insert = build_library_insert_distribution(library);

    const auto preparation_mark = std::chrono::steady_clock::now();
    const PreparedFragmentSet prepared = prepare_fragment_set(fragments, library, options);
    counters.fragments_with_floor_binding = prepared.fragments_with_floor_binding;
    result.seconds_fragment_preparation =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - preparation_mark).count();

    const auto placement_mark = std::chrono::steady_clock::now();
    const std::size_t unique = haplotypes.sequences.size();
    if (!fragments.empty() && unique > options.max_mass_cells / fragments.size()) {
        throw ResourceExhausted(
            "the fragment-by-haplotype mass matrix would exceed its declared cell cap; no partial "
            "ranking is returned");
    }
    std::vector<std::vector<double>> mass_column(unique,
                                                 std::vector<double>(fragments.size(), kNegInf));
    std::vector<double> exposure(unique, 0.0);
    if (options.record_origin_summary) {
        result.origin_summary.assign(fragments.size() * unique, OriginSummary{});
    }
    std::vector<char> record_origins(fragments.size(), 0);
    for (const std::size_t f : options.record_origins_for_fragments) {
        if (f >= fragments.size()) {
            throw std::invalid_argument("origin recording names an unknown fragment");
        }
        record_origins[f] = 1;
    }
    for (std::size_t h = 0; h < unique; ++h) {
        exposure[h] = exposure_for_length(haplotypes.sequences[h].size(), fragments.size(),
                                          prepared.fragments_by_floor, library_insert);
    }

    // Stable storage, appended once here. A deque so the emission cache's raw window pointers stay
    // valid however the container grows -- the one-shot path never grows it, but the definition is
    // shared with the context, which does.
    std::deque<std::string> sequences;
    std::deque<std::string> reversed_sequences;
    for (const std::string& sequence : haplotypes.sequences) {
        sequences.push_back(sequence);
        reversed_sequences.emplace_back(sequence.rbegin(), sequence.rend());
    }

    struct EmissionCacheScope {
        explicit EmissionCacheScope(std::size_t bytes) { begin_emission_cache(bytes); }
        ~EmissionCacheScope() { counters = end_emission_cache(); }
        EmissionCacheScope(const EmissionCacheScope&) = delete;
        EmissionCacheScope& operator=(const EmissionCacheScope&) = delete;
        EmissionCacheCounters counters;
    } emission_cache_scope(options.emission_cache_bytes);

    PlacementBudget budget;
    budget.max_piece_hits = options.max_piece_hits;
    budget.max_physical_origins = options.max_physical_origins;
    budget.publish_batch = std::max<std::uint64_t>(1, options.piece_hit_publish_batch);

    PlacementInputs inputs;
    inputs.prepared = &prepared;
    inputs.fragments = &fragments;
    inputs.library = &library;
    inputs.library_insert = &library_insert;
    inputs.options = &options;
    inputs.sequences = &sequences;
    inputs.reversed = &reversed_sequences;
    inputs.band = band;

    OriginRecording origins_out;
    origins_out.summary = options.record_origin_summary;
    origins_out.record_origins = &record_origins;
    origins_out.unique_stride = unique;
    origins_out.summaries = &result.origin_summary;
    origins_out.recorded = &result.recorded_origins;

    const bool origins_force_serial = !options.record_origins_for_fragments.empty();
    UnionReplay union_replay;
    std::vector<UnionSegmentEntry> union_entries;
    UnionOriginReuse union_origins;
    std::deque<UnionOriginEntry> union_origin_entries;
    const bool union_phase_one = options.union_segment_lengths != nullptr && options.union_placement_cache_bytes > 0;
    const bool union_phase_two = union_phase_one && options.union_origin_cache_bytes > 0;
    if (options.union_origin_cache_bytes > 0) {
        if (!union_phase_one) {
            throw std::invalid_argument("union origin reuse needs union placement reuse (a segmentation and its budget)");
        }
        if (options.max_union_cached_origins_consumed == 0 || options.max_union_cached_states_consumed == 0) {
            throw std::invalid_argument("union origin reuse needs both of its consumption caps, positive");
        }
        if (options.record_origin_summary || !options.record_origins_for_fragments.empty()) {
            throw std::invalid_argument("union origin reuse does not record per-origin diagnostics");
        }
    }
    union_origins.max_origins = options.max_union_cached_origins_consumed;
    union_origins.max_states = options.max_union_cached_states_consumed;
    {
        StartSortTally& tally = start_sort_tally();
        tally.out_of_place.store(0);
        tally.in_place.store(0);
        tally.fallbacks.store(0);
        tally.peak_worker_bytes.store(0);
        result.start_sort.budget_bytes = options.start_sort_scratch_bytes;
        result.start_sort.workers = std::max<std::size_t>(1, options.threads);
        result.start_sort.share_bytes = options.start_sort_scratch_bytes / result.start_sort.workers;
    }
    const auto start_sort_report = [&](StartSortReport& report) {
        const StartSortTally& tally = start_sort_tally();
        report.out_of_place = tally.out_of_place.load();
        report.in_place = tally.in_place.load();
        report.fallbacks = tally.fallbacks.load();
        report.peak_worker_bytes = tally.peak_worker_bytes.load();
    };
    std::unique_ptr<UnionPrefilterCache> prefilter_cache;
    if (options.union_prefilter_cache_bytes > 0) {
        prefilter_cache = make_union_prefilter_cache(prepared, library, options, result.union_prefilter_cache);
        inputs.prefilter_cache = prefilter_cache.get();
    }
    const auto prefilter_report = [&](UnionPrefilterCacheReport& report) {
        if (prefilter_cache) fill_union_prefilter_report(*prefilter_cache, report);
    };

    const auto note_observed_work = [&]() {
        if (observed_on_failure == nullptr) return;
        observed_on_failure->counters = counters;
        observed_on_failure->counters.piece_hits =
            std::max(counters.piece_hits, budget.piece_hits.load(std::memory_order_relaxed));
        observed_on_failure->counters.physical_origins = std::max(
            counters.physical_origins, budget.physical_origins.load(std::memory_order_relaxed));
        observed_on_failure->threads_used_placement = result.threads_used_placement;
        observed_on_failure->threads_used_combination = 0;
        observed_on_failure->seconds_fragment_preparation = result.seconds_fragment_preparation;
        observed_on_failure->seconds_placement_and_mass =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - placement_mark).count();
        observed_on_failure->seconds_diploid_combination = 0.0;
        observed_on_failure->union_placement_cache = result.union_placement_cache;
        observed_on_failure->union_placement_cache.placements_replayed = union_replay.consumed.load(std::memory_order_relaxed);
        observed_on_failure->union_origin_cache = result.union_origin_cache;
        observed_on_failure->union_origin_cache.origins_consumed =
            union_origins.origins_consumed.load(std::memory_order_relaxed);
        observed_on_failure->union_origin_cache.states_consumed =
            union_origins.states_consumed.load(std::memory_order_relaxed);
        observed_on_failure->union_prefilter_cache = result.union_prefilter_cache;
        prefilter_report(observed_on_failure->union_prefilter_cache);
        observed_on_failure->start_sort = result.start_sort;
        start_sort_report(observed_on_failure->start_sort);
    };

    result.threads_used_placement =
        (unique == 0 || origins_force_serial)
            ? 1
            : std::max<std::size_t>(1, std::min(options.threads, unique));
    try {
        if (options.record_haplotype_work != nullptr) options.record_haplotype_work->assign(unique, HaplotypeWork{});
        if (union_phase_one && unique > 0 && !prepared.forward_forms.empty()) {
            fill_union_cache(inputs, *options.union_segment_lengths, options.union_placement_cache_bytes,
                             options.max_union_cached_placements_consumed, options.threads, budget, counters,
                             result.emission_cache, union_entries, union_replay, result.union_placement_cache,
                             options.union_origin_cache_bytes, union_phase_two ? &union_origins : nullptr,
                             &union_origin_entries, &result.union_origin_cache);
            inputs.replay = &union_replay;
            if (union_phase_two) inputs.origin_reuse = &union_origins;
        }
        result.threads_used_placement =
            run_placement(inputs, 0, unique, options.threads, origins_force_serial, budget,
                          mass_column, counters, result.memory, result.emission_cache, origins_out);
        result.union_placement_cache.placements_replayed = union_replay.consumed.load(std::memory_order_relaxed);
        result.union_origin_cache.origins_consumed = union_origins.origins_consumed.load(std::memory_order_relaxed);
        result.union_origin_cache.states_consumed = union_origins.states_consumed.load(std::memory_order_relaxed);
        prefilter_report(result.union_prefilter_cache);
        start_sort_report(result.start_sort);
    } catch (const ResourceExhausted&) {
        note_observed_work();
        throw;
    }

    enforce_acceptance_caps(counters, options, note_observed_work);

    result.seconds_placement_and_mass =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - placement_mark).count();

    combine_candidate_scores(mass_column, exposure, prepared.background_by_fragment, fragments.size(), library,
                             candidates, options.threads, options.audit_libraries, result, options.prepared_library_terms);

    // Flattened to the documented row-major layout. Column storage is an internal choice about what
    // can be appended; the published matrix keeps its shape.
    result.log_mass.assign(fragments.size() * unique, kNegInf);
    for (std::size_t h = 0; h < unique; ++h) {
        for (std::size_t f = 0; f < fragments.size(); ++f) {
            result.log_mass[f * unique + h] = mass_column[h][f];
        }
    }
    result.exposure = std::move(exposure);
    {
        const EmissionCacheCounters stored = read_emission_cache_counters();
        result.emission_cache.payload_bytes = stored.payload_bytes;
        result.emission_cache.entries = stored.entries;
        result.emission_cache.entries_not_retained = stored.entries_not_retained;
        result.emission_cache.budgeted_bytes_used = stored.budgeted_bytes_used;
    }
    emission_cache_storage_bytes(&result.memory.emission_cache_table_bytes,
                                 &result.memory.emission_cache_arena_bytes);
    result.memory.mass_matrix_bytes =
        static_cast<std::uint64_t>(result.log_mass.size()) * sizeof(double);
    for (const std::string& sequence : sequences) {
        result.memory.haplotype_sequence_bytes += sequence.capacity();
    }
    for (const std::string& sequence : reversed_sequences) {
        result.memory.haplotype_sequence_bytes += sequence.capacity();
    }
    for (const ReadForm& form : prepared.forward_forms) {
        result.memory.form_sequence_bytes += form.bases.capacity();
    }
    for (const ReadForm& form : prepared.reverse_forms) {
        result.memory.form_sequence_bytes += form.bases.capacity();
    }
    for (const PieceIndex& index : prepared.forward_indexes) {
        result.memory.piece_index_bytes += index.entries.capacity() * sizeof(PieceEntry);
    }
    for (const PieceIndex& index : prepared.reverse_indexes) {
        result.memory.piece_index_bytes += index.entries.capacity() * sizeof(PieceEntry);
    }
    last_memory_slot() = result.memory;
    last_emission_cache_slot() = result.emission_cache;
    return result;
}

UnionSeedingProfile profile_union_seeding(const InternedHaplotypes& haplotypes, const std::vector<ReadPair>& fragments,
                                          const LibraryModel& library, const FragmentScoringOptions& options,
                                          const std::vector<std::size_t>& sample, const std::vector<char>& is_panel) {
    using Clock = std::chrono::steady_clock;
    const auto ns_since = [](Clock::time_point mark) {
        return static_cast<std::uint64_t>(
            std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - mark).count());
    };
    validate_library_model(library);
    UnionSeedingProfile profile;
    const std::size_t band = options.mismatch_band;
    const DiscreteInsertDistribution library_insert = build_library_insert_distribution(library);
    const PreparedFragmentSet prepared = prepare_fragment_set(fragments, library, options);
    std::deque<std::string> sequences;
    std::deque<std::string> reversed_sequences;
    for (const std::string& sequence : haplotypes.sequences) {
        sequences.push_back(sequence);
        reversed_sequences.emplace_back(sequence.rbegin(), sequence.rend());
    }
    struct EmissionCacheScope {
        explicit EmissionCacheScope(std::size_t bytes) { begin_emission_cache(bytes); }
        ~EmissionCacheScope() { (void)end_emission_cache(); }
    } emission_cache_scope(options.emission_cache_bytes);
    PlacementBudget budget;
    budget.max_piece_hits = std::numeric_limits<std::uint64_t>::max();
    budget.max_physical_origins = std::numeric_limits<std::uint64_t>::max();
    PlacementInputs inputs;
    inputs.prepared = &prepared;
    inputs.fragments = &fragments;
    inputs.library = &library;
    inputs.library_insert = &library_insert;
    inputs.options = &options;
    inputs.sequences = &sequences;
    inputs.reversed = &reversed_sequences;
    inputs.band = band;
    UnionReplay replay;
    std::vector<UnionSegmentEntry> entries;
    UnionOriginReuse origin_reuse;
    origin_reuse.max_origins = std::numeric_limits<std::uint64_t>::max();
    origin_reuse.max_states = std::numeric_limits<std::uint64_t>::max();
    std::deque<UnionOriginEntry> origin_entries;
    UnionPlacementCacheReport placement_report;
    UnionOriginCacheReport origin_report;
    FragmentScoringCounters fill_counters;
    EmissionCacheCounters fill_cache;
    const bool reuse = options.union_segment_lengths != nullptr && options.union_placement_cache_bytes > 0;
    if (reuse) {
        fill_union_cache(inputs, *options.union_segment_lengths, options.union_placement_cache_bytes, 0, options.threads,
                         budget, fill_counters, fill_cache, entries, replay, placement_report,
                         0, nullptr, &origin_entries, &origin_report, &profile.fill_segments);   // phase 1 only
        profile.fill_workers = std::max<std::size_t>(1, std::min<std::size_t>(options.threads, placement_report.eligible_segments));
        profile.fill_seconds = placement_report.seconds_fill;
    }
    const auto code_of = [](char base) { return two_bit(base); };
    std::vector<std::uint8_t> skip;
    std::vector<std::pair<std::uint32_t, std::uint32_t>> proposed, generated, sorted;
    std::vector<std::uint32_t> counts, bounds;
    for (const std::size_t h : sample) {
        if (h >= sequences.size()) throw std::invalid_argument("seeding profile samples an unknown haplotype");
        profile.sampled.push_back(h);
        std::uint64_t hash = 1469598103934665603ULL;
        for (const char base : sequences[h]) hash = (hash ^ static_cast<unsigned char>(base)) * 1099511628211ULL;
        profile.sampled_hashes.push_back(hash);
        UnionSeedingPhases& phases = (h < is_panel.size() && is_panel[h] != 0) ? profile.panel : profile.mosaic;
        ++phases.haplotypes;
        for (int strand = 0; strand < 2; ++strand) {
            const std::string& text = strand == 0 ? sequences[h] : reversed_sequences[h];
            const auto& forms = strand == 0 ? prepared.forward_forms : prepared.reverse_forms;
            const auto& componentise = strand == 0 ? prepared.forward_componentise : prepared.reverse_componentise;
            const auto& indexes = strand == 0 ? prepared.forward_indexes : prepared.reverse_indexes;
            const std::vector<ReplayZone>* zones =
                reuse ? &(strand == 0 ? replay.forward : replay.reverse)[h] : nullptr;
            const bool replaying = zones != nullptr && !zones->empty();

            // The production call, as the pass makes it (it also warms the emission cache for these keys).
            std::vector<std::vector<FormPlacement>> rows(forms.size());
            FragmentScoringCounters counters;
            EmissionCacheCounters cache_work;
            (void)take_gapped_emission_counters();
            auto mark = Clock::now();
            collect_form_placements(text, forms, componentise, indexes, budget, band, library, rows, proposed, counters,
                                    cache_work, options, nullptr, zones, &skip, replaying ? &replay : nullptr);
            phases.production_ns += ns_since(mark);
            phases.piece_hits += counters.piece_hits;
            phases.proposals += counters.proposed_starts;
            phases.unique_starts += counters.unique_starts;
            phases.prefilters += cache_work.executed_prefilters;
            phases.prefilter_rejections += counters.prefilter_rejections;
            phases.dp_requests += counters.dp_evaluations;
            phases.executed_dps += cache_work.executed_dps;

            // Seed scan alone, then the same scan generating starts: the difference is generation.
            if (replaying) {
                skip.assign(text.size() + library.max_gap_bp + 2, 0);
                for (const ReplayZone& zone : *zones) {
                    std::fill(skip.begin() + zone.lo, skip.begin() + zone.hi + 1, static_cast<std::uint8_t>(1));
                }
            }
            for (int generate = 0; generate < 2; ++generate) {
                generated.clear();
                std::uint64_t sink = 0;
                mark = Clock::now();
                for (const PieceIndex& index : indexes) {
                    const std::size_t piece = index.piece_bp;
                    if (text.size() < piece) continue;
                    const std::uint64_t mask =
                        piece == 32 ? std::numeric_limits<std::uint64_t>::max() : ((std::uint64_t{1} << (2 * piece)) - 1);
                    std::uint64_t code = 0;
                    std::size_t run = 0;
                    for (std::size_t at = 0; at < text.size(); ++at) {
                        const std::uint64_t bits = code_of(text[at]);
                        if (bits == kNoCode) { run = 0; code = 0; continue; }
                        code = ((code << 2) | bits) & mask;
                        if (++run < piece) continue;
                        const std::size_t piece_start = at + 1 - piece;
                        const auto range = std::equal_range(
                            index.entries.begin(), index.entries.end(), code, [](const auto& lhs, const auto& rhs) {
                                if constexpr (std::is_same_v<std::decay_t<decltype(lhs)>, PieceEntry>) {
                                    return lhs.code < rhs;
                                } else {
                                    return lhs < rhs.code;
                                }
                            });
                        for (auto it = range.first; it != range.second; ++it) {
                            if (it->offset > piece_start) continue;
                            const std::size_t implied = piece_start - it->offset;
                            const std::size_t lo = implied > library.max_gap_bp ? implied - library.max_gap_bp : 0;
                            const std::size_t hi = implied + library.max_gap_bp;
                            if (generate == 0) {
                                sink += lo + hi;
                                continue;
                            }
                            for (std::size_t shifted = lo; shifted <= hi; ++shifted) {
                                if (replaying && shifted < skip.size() && skip[shifted] != 0) continue;
                                generated.emplace_back(it->form, static_cast<std::uint32_t>(shifted));
                            }
                        }
                    }
                }
                (generate == 0 ? phases.scan_ns : phases.scan_and_generate_ns) += ns_since(mark);
                if (generate == 0) {
                    volatile std::uint64_t keep = sink;   // keeps the scan-only loop from being elided
                    (void)keep;
                }
            }
            if (generated.size() != counters.proposed_starts) {
                throw std::logic_error("seeding profile: the replayed scan proposed a different number of starts");
            }
            // Start ordering: production's sort + unique, and the bucketed alternative, on the same proposals.
            sorted = generated;
            mark = Clock::now();
            std::sort(sorted.begin(), sorted.end());
            sorted.erase(std::unique(sorted.begin(), sorted.end()), sorted.end());
            phases.start_sort_ns += ns_since(mark);
            if (sorted != proposed) throw std::logic_error("seeding profile: the replayed sort differs from production");
            std::vector<std::pair<std::uint32_t, std::uint32_t>> bucketed = generated;
            mark = Clock::now();
            sort_unique_proposals_bucketed(bucketed, forms.size(), counts, bounds);
            phases.start_sort_bucketed_ns += ns_since(mark);
            if (bucketed != proposed) throw std::logic_error("seeding profile: the bucketed start set differs");
            ++profile.bucketed_sort_checked;

            // Per unique start: the window checks, the prefilter (both kernels), the emission call and the storage.
            std::vector<std::size_t> window_of(proposed.size(), 0);
            mark = Clock::now();
            for (std::size_t i = 0; i < proposed.size(); ++i) {
                const std::string& read = forms[proposed[i].first].bases;
                const std::size_t start = proposed[i].second;
                if (start > text.size() || text.size() - start < read.size()) continue;
                const std::size_t window = std::min(text.size() - start, read.size() + library.max_gap_bp);
                require_acgt_window(text, start, window, read);
                window_of[i] = window;
            }
            phases.window_check_ns += ns_since(mark);
            std::vector<char> feasible(proposed.size(), 0);
            mark = Clock::now();
            for (std::size_t i = 0; i < proposed.size(); ++i) {
                if (window_of[i] == 0) continue;
                const ReadForm& form = forms[proposed[i].first];
                feasible[i] = feasible_within_band_reference(text, proposed[i].second, window_of[i], form.bases, band,
                                                             &form.pattern) ? 1 : 0;
            }
            phases.prefilter_ns += ns_since(mark);
            std::vector<char> early(proposed.size(), 0);
            mark = Clock::now();
            for (std::size_t i = 0; i < proposed.size(); ++i) {
                if (window_of[i] == 0) continue;
                const ReadForm& form = forms[proposed[i].first];
                early[i] = feasible_within_band(text, proposed[i].second, window_of[i], form.bases, band, &form.pattern) ? 1 : 0;
            }
            phases.prefilter_early_ns += ns_since(mark);
            if (early != feasible) throw std::logic_error("seeding profile: the early-refusing prefilter differs");
            profile.early_prefilter_checked += proposed.size();

            // The same eligible starts streamed through 2, 4 and 8 lanes, one read form at a time.
            {
                std::vector<std::uint32_t> lane_starts, lane_windows;
                std::vector<char> lane_out;
                std::vector<std::pair<std::size_t, std::size_t>> groups;   // [begin, end) in the compacted arrays
                std::vector<std::size_t> origin_of;                          // compacted -> proposed index
                std::vector<const ReadForm*> group_form;
                for (std::size_t i = 0; i < proposed.size();) {
                    std::size_t k = i;
                    const std::size_t first = lane_starts.size();
                    while (k < proposed.size() && proposed[k].first == proposed[i].first) {
                        if (window_of[k] != 0) {
                            lane_starts.push_back(proposed[k].second);
                            lane_windows.push_back(static_cast<std::uint32_t>(window_of[k]));
                            origin_of.push_back(k);
                        }
                        ++k;
                    }
                    const ReadForm& form = forms[proposed[i].first];
                    if (lane_starts.size() > first) {
                        if (!form.pattern.usable) throw std::logic_error("seeding profile: a form without a pattern");
                        groups.push_back({first, lane_starts.size()});
                        group_form.push_back(&form);
                    }
                    i = k;
                }
                phases.batched_starts += lane_starts.size();
                const auto run_lanes = [&](auto lanes_tag, std::uint64_t& into) {
                    constexpr std::size_t kLanes = decltype(lanes_tag)::value;
                    lane_out.assign(lane_starts.size(), 2);
                    const auto lane_mark = Clock::now();
                    for (std::size_t g = 0; g < groups.size(); ++g) {
                        const MyersPattern& pattern = group_form[g]->pattern;
                        MyersWord peq[4];
                        for (int c = 0; c < 4; ++c) peq[c] = myers_word_of(pattern.peq_lo[c], pattern.peq_hi[c]);
                        myers_early_lanes<kLanes>(peq, pattern.length, text.data(), lane_starts.data() + groups[g].first,
                                                  lane_windows.data() + groups[g].first,
                                                  groups[g].second - groups[g].first, band,
                                                  lane_out.data() + groups[g].first);
                    }
                    into += ns_since(lane_mark);
                    for (std::size_t c = 0; c < lane_out.size(); ++c) {
                        if (lane_out[c] != early[origin_of[c]]) {
                            throw std::logic_error("seeding profile: a multi-lane prefilter decision differs");
                        }
                    }
                    profile.lane_decisions_checked += lane_out.size();
                };
                run_lanes(std::integral_constant<std::size_t, 2>{}, phases.lanes2_ns);
                run_lanes(std::integral_constant<std::size_t, 4>{}, phases.lanes4_ns);
                run_lanes(std::integral_constant<std::size_t, 8>{}, phases.lanes8_ns);
            }

            // The emission key path over the survivors: every key is cached by the production call above.
            {
                std::vector<std::size_t> survivors;
                for (std::size_t i = 0; i < proposed.size(); ++i) {
                    if (feasible[i] != 0) survivors.push_back(i);
                }
                phases.survivors += survivors.size();
                std::vector<KeyView> keys(survivors.size());
                mark = Clock::now();
                for (std::size_t k = 0; k < survivors.size(); ++k) {
                    const std::size_t i = survivors[k];
                    keys[k] = key_view_of(forms[proposed[i].first].bases, text.data() + proposed[i].second, window_of[i]);
                }
                phases.key_full_ns += ns_since(mark);
                std::vector<KeyPrefix> prefix_of(forms.size());
                std::vector<char> prefixed(forms.size(), 0);
                mark = Clock::now();
                for (std::size_t k = 0; k < survivors.size(); ++k) {
                    const std::size_t i = survivors[k];
                    const std::uint32_t f = proposed[i].first;
                    if (prefixed[f] == 0) {
                        prefix_of[f] = key_prefix_of(forms[f].bases);
                        prefixed[f] = 1;
                    }
                    const KeyView key = key_view_from_prefix(prefix_of[f], forms[f].bases, text.data() + proposed[i].second,
                                                             window_of[i]);
                    if (key.h0 != keys[k].h0 || key.h1 != keys[k].h1) {
                        throw std::logic_error("seeding profile: a prefix-state key differs from the full hash");
                    }
                }
                phases.key_prefix_ns += ns_since(mark);
                profile.prefix_keys_checked += survivors.size();
                GappedFit out;
                EmissionWork scratch_work;
                std::uint64_t sink = 0;
                mark = Clock::now();
                for (const KeyView& key : keys) sink += static_cast<std::uint64_t>(cache_lookup(key, &out, scratch_work));
                phases.lookup_ns += ns_since(mark);
                const auto parts = [&](bool lock, bool copy, std::uint64_t& into) {
                    const auto part_mark = Clock::now();
                    for (const KeyView& key : keys) sink += static_cast<std::uint64_t>(cache_probe_measured(key, &out, lock, copy));
                    into += ns_since(part_mark);
                };
                parts(true, true, phases.lookup_unclocked_ns);
                parts(false, true, phases.lookup_nolock_ns);
                parts(true, false, phases.lookup_nocopy_ns);
                parts(false, false, phases.probe_only_ns);
                mark = Clock::now();
                for (std::size_t k = 0; k < keys.size(); ++k) {
                    const auto a = Clock::now();
                    const auto b = Clock::now();
                    sink += static_cast<std::uint64_t>((b - a).count());
                }
                phases.clock_pairs_ns += ns_since(mark);
                volatile std::uint64_t keep = sink;
                (void)keep;
            }
            // Only accepted fits are kept (a fit carries its component array inline; millions of refusals would not).
            std::vector<std::pair<std::size_t, GappedFit>> fits;
            mark = Clock::now();
            for (std::size_t i = 0; i < proposed.size(); ++i) {
                const ReadForm& form = forms[proposed[i].first];
                std::optional<GappedFit> fit =
                    gapped_emission(text, form.bases, proposed[i].second, library, band, &form.pattern);
                if (fit.has_value() && fit->log_emission != kNegInf) fits.emplace_back(i, std::move(*fit));
            }
            phases.emission_call_ns += ns_since(mark);
            (void)take_gapped_emission_counters();
            // Direct recursions over a fixed prefix of the survivors, for a per-DP cost.
            std::size_t timed = 0;
            mark = Clock::now();
            for (std::size_t i = 0; i < proposed.size() && timed < 2000; ++i) {
                if (feasible[i] == 0) continue;
                const ReadForm& form = forms[proposed[i].first];
                GappedFit fit;
                (void)gap_recursion_scaled(&text, &form.bases, proposed[i].second, form.bases.size(),
                                           text.size() - proposed[i].second, band, library, true, &fit.components);
                ++timed;
            }
            phases.dp_timed_ns += ns_since(mark);
            phases.dp_timed += timed;
            // Storage, replay and the row sort, on fresh rows.
            for (std::vector<FormPlacement>& row : rows) row.clear();
            mark = Clock::now();
            for (const auto& [i, accepted] : fits) {
                const GappedFit* fit = &accepted;
                const std::uint32_t form = proposed[i].first;
                if (componentise[form] != 0) {
                    for (const GappedComponent& component : fit->components) {
                        if (component.log_emission == kNegInf) continue;
                        rows[form].push_back({proposed[i].second, static_cast<std::uint32_t>(component.consumed_bp),
                                              static_cast<std::uint32_t>(component.edits), component.log_emission});
                        ++phases.placements;
                    }
                    continue;
                }
                rows[form].push_back({proposed[i].second, static_cast<std::uint32_t>(fit->consumed_bp),
                                      static_cast<std::uint32_t>(fit->edits), fit->log_emission});
                ++phases.placements;
            }
            phases.storage_ns += ns_since(mark);
            mark = Clock::now();
            if (replaying) {
                for (const ReplayZone& zone : *zones) {
                    for (const CachedPlacement& cached : zone.entry->placements) {
                        rows[cached.form].push_back(
                            {cached.start + zone.shift, cached.consumed, cached.mismatches, cached.emission});
                        ++phases.replayed;
                    }
                }
            }
            phases.replay_ns += ns_since(mark);
            mark = Clock::now();
            for (std::vector<FormPlacement>& list : rows) {
                std::sort(list.begin(), list.end(), [](const FormPlacement& lhs, const FormPlacement& rhs) {
                    return std::tie(lhs.start, lhs.consumed) < std::tie(rhs.start, rhs.consumed);
                });
            }
            phases.row_sort_ns += ns_since(mark);
        }
    }
    return profile;
}

UnionVerificationCensus census_union_verification(const InternedHaplotypes& haplotypes,
                                                  const std::vector<ReadPair>& fragments, const LibraryModel& library,
                                                  const FragmentScoringOptions& options,
                                                  const std::vector<char>& is_panel, const std::string& scratch_dir) {
    using Clock = std::chrono::steady_clock;
    validate_library_model(library);
    UnionVerificationCensus census;
    const std::size_t band = options.mismatch_band;
    const std::size_t g = library.max_gap_bp;
    const DiscreteInsertDistribution library_insert = build_library_insert_distribution(library);
    const PreparedFragmentSet prepared = prepare_fragment_set(fragments, library, options);
    std::deque<std::string> sequences;
    std::deque<std::string> reversed_sequences;
    for (const std::string& sequence : haplotypes.sequences) {
        sequences.push_back(sequence);
        reversed_sequences.emplace_back(sequence.rbegin(), sequence.rend());
    }
    const std::size_t unique = sequences.size();
    // One read length: then every eligible start has exactly one window.
    std::size_t r = 0;
    for (const auto* forms : {&prepared.forward_forms, &prepared.reverse_forms}) {
        for (const ReadForm& form : *forms) {
            if (r == 0) r = form.bases.size();
            if (form.bases.size() != r) throw std::invalid_argument("the verification census needs one read length");
            if (!form.pattern.usable || !form.all_acgt) throw std::invalid_argument("the census needs prepared ACGT reads");
        }
    }
    for (const std::string& sequence : sequences) {
        if (sequence.size() >= (std::size_t{1} << 19)) throw std::invalid_argument("the census packs haplotypes below 2^19 bp");
        if (!std::all_of(sequence.begin(), sequence.end(), is_acgt)) {
            throw std::invalid_argument("the census needs ACGT haplotypes");
        }
    }
    if (unique >= (std::size_t{1} << 9)) throw std::invalid_argument("the census packs at most 511 haplotypes");
    if (unique >= (std::size_t{1} << 12)) throw std::invalid_argument("the census orders at most 4095 haplotypes");
    census.read_length = r;
    census.window_length = r + g;

    // Phase-1 zones, exactly as the pass builds them (the fill's own verification work is recorded apart).
    PlacementBudget budget;
    budget.max_piece_hits = std::numeric_limits<std::uint64_t>::max();
    budget.max_physical_origins = std::numeric_limits<std::uint64_t>::max();
    PlacementInputs inputs;
    inputs.prepared = &prepared;
    inputs.fragments = &fragments;
    inputs.library = &library;
    inputs.library_insert = &library_insert;
    inputs.options = &options;
    inputs.sequences = &sequences;
    inputs.reversed = &reversed_sequences;
    inputs.band = band;
    UnionReplay replay;
    std::vector<UnionSegmentEntry> entries;
    std::deque<UnionOriginEntry> origin_entries;
    UnionPlacementCacheReport placement_report;
    UnionOriginCacheReport origin_report;
    const bool reuse = options.union_segment_lengths != nullptr && options.union_placement_cache_bytes > 0;
    if (reuse) {
        FragmentScoringCounters fill_counters;
        EmissionCacheCounters fill_cache;
        // The emission cache changes how fast the fill runs, never which placements it keeps, so it stays off here:
        // its tables would otherwise stay resident through the census.
        begin_emission_cache(0);
        fill_union_cache(inputs, *options.union_segment_lengths, options.union_placement_cache_bytes, 0, options.threads,
                         budget, fill_counters, fill_cache, entries, replay, placement_report, 0, nullptr,
                         &origin_entries, &origin_report);
        (void)end_emission_cache();
        census.fill_executed_prefilters = fill_cache.executed_prefilters;
        census.fill_verifications = fill_counters.full_read_verifications;
    }

    // Read contents, by complete byte equality.
    std::unordered_map<std::string_view, std::uint32_t> read_id_of_bytes;
    std::vector<std::uint32_t> read_ids[2];
    for (int strand = 0; strand < 2; ++strand) {
        const auto& forms = strand == 0 ? prepared.forward_forms : prepared.reverse_forms;
        for (const ReadForm& form : forms) {
            const auto inserted = read_id_of_bytes.emplace(std::string_view(form.bases),
                                                           static_cast<std::uint32_t>(read_id_of_bytes.size()));
            read_ids[strand].push_back(inserted.first->second);
        }
        census.read_forms += forms.size();
    }
    census.read_contents = read_id_of_bytes.size();

    // Window contents: every eligible window (N - p >= r) of every strand, length min(N - p, r + g).
    auto mark = Clock::now();
    const std::size_t W = r + g;
    std::size_t positions = 0;
    for (const std::string& sequence : sequences) positions += 2 * (sequence.size() >= r ? sequence.size() - r + 1 : 0);
    census.window_positions = positions;
    std::size_t bits = 20;
    while ((std::size_t{1} << bits) < positions) ++bits;   // load <= 1/2 even if every window were distinct
    std::vector<std::uint64_t> table(std::size_t{1} << bits, 0);
    census.window_table_bytes = table.size() * sizeof(std::uint64_t);
    const std::uint64_t mask = (std::uint64_t{1} << bits) - 1;
    const auto text_of = [&](std::size_t h, int strand) -> const std::string& {
        return strand == 0 ? sequences[h] : reversed_sequences[h];
    };
    std::vector<std::vector<std::uint32_t>> window_ids(2 * unique);
    std::uint64_t next_id = 0;
    constexpr std::uint64_t kBase = 0x100000001b3ULL;
    std::uint64_t base_power = 1;   // kBase^(W-1)
    for (std::size_t i = 1; i < W; ++i) base_power *= kBase;
    const auto code = [](char c) -> std::uint64_t { return static_cast<unsigned char>(c); };
    for (std::size_t h = 0; h < unique; ++h) {
        for (int strand = 0; strand < 2; ++strand) {
            const std::string& text = text_of(h, strand);
            const std::size_t N = text.size();
            std::vector<std::uint32_t>& ids = window_ids[2 * h + strand];
            if (N < r) continue;
            ids.assign(N - r + 1, 0);
            std::uint64_t rolling = 0;
            if (N >= W) for (std::size_t i = 0; i < W; ++i) rolling = rolling * kBase + code(text[i]);
            for (std::size_t pos = 0; pos + r <= N; ++pos) {
                const std::size_t length = std::min(N - pos, W);
                std::uint64_t hash = 0;
                if (length == W) {
                    if (pos > 0) rolling = (rolling - code(text[pos - 1]) * base_power) * kBase + code(text[pos + W - 1]);
                    hash = rolling;
                } else {
                    for (std::size_t i = 0; i < length; ++i) hash = hash * kBase + code(text[pos + i]);
                }
                hash ^= static_cast<std::uint64_t>(length) * 0x9e3779b97f4a7c15ULL;
                std::uint64_t slot = (hash * 0xbf58476d1ce4e5b9ULL) >> (64 - bits);
                for (;;) {
                    const std::uint64_t entry = table[slot];
                    if (entry == 0) {
                        const std::uint64_t rep = (static_cast<std::uint64_t>(h) << 20) |
                                                  (static_cast<std::uint64_t>(strand) << 19) | pos;
                        table[slot] = ((next_id + 1) << 37) | (static_cast<std::uint64_t>(length) << 29) | rep;
                        ids[pos] = static_cast<std::uint32_t>(next_id++);
                        break;
                    }
                    const std::size_t entry_length = (entry >> 29) & 0xff;
                    if (entry_length == length) {
                        const std::uint64_t rep = entry & ((std::uint64_t{1} << 29) - 1);
                        const std::string& other = text_of(rep >> 20, static_cast<int>((rep >> 19) & 1));
                        if (std::memcmp(other.data() + (rep & ((1u << 19) - 1)), text.data() + pos, length) == 0) {
                            ids[pos] = static_cast<std::uint32_t>((entry >> 37) - 1);
                            break;
                        }
                    }
                    slot = (slot + 1) & mask;
                }
            }
        }
    }
    if (next_id >= (std::uint64_t{1} << 27) - 1) throw std::logic_error("the census window ids overflow 27 bits");
    census.window_contents = next_id;
    census.window_index_bytes = positions * sizeof(std::uint32_t);
    std::vector<std::uint64_t>().swap(table);
    {
        std::vector<char> seen(next_id, 0);
        std::uint64_t distinct[2] = {0, 0};
        for (std::size_t h = 0; h < unique; ++h) {
            const int cls = h < is_panel.size() && is_panel[h] != 0 ? 0 : 1;
            for (int strand = 0; strand < 2; ++strand) {
                for (const std::uint32_t id : window_ids[2 * h + strand]) {
                    if (seen[id] & (1 << cls)) continue;
                    seen[id] |= static_cast<char>(1 << cls);
                    ++distinct[cls];
                }
            }
            (cls == 0 ? census.accumulation_panel : census.accumulation_mosaic).push_back(distinct[cls]);
        }
    }
    census.seconds_windows = std::chrono::duration<double>(Clock::now() - mark).count();

    // The stream: every haplotype's unique starts, as the pass proposes and orders them, each request's decision
    // and columns written to one of kParts partitions by its identity.
    struct Record {
        std::uint64_t key;      // read id << 27 | window id
        std::uint32_t order;    // haplotype << 20 | strand << 19 | start: the canonical first occurrence is the least
        std::uint16_t columns;
        std::uint16_t flags;    // 1 accepted, 2 panel
    };
    static_assert(sizeof(Record) == 16, "a census record is 16 bytes");
    constexpr std::size_t kParts = 64;
    std::filesystem::create_directories(scratch_dir);
    std::vector<std::FILE*> files(kParts, nullptr);
    std::vector<std::mutex> file_guard(kParts);
    for (std::size_t part = 0; part < kParts; ++part) {
        files[part] = std::fopen((scratch_dir + "/part" + std::to_string(part) + ".bin").c_str(), "wb");
        if (files[part] == nullptr) throw std::runtime_error("the census cannot write its scratch partitions");
    }
    mark = Clock::now();
    std::atomic<std::size_t> next_haplotype{0};
    std::vector<UnionCensusClass> local_classes(2 * std::max<std::size_t>(1, options.threads));
    std::vector<std::uint64_t> local_checked(std::max<std::size_t>(1, options.threads), 0);
    std::vector<std::exception_ptr> failures(std::max<std::size_t>(1, options.threads));
    const auto work = [&](std::size_t worker) {
        try {
            std::vector<std::vector<Record>> buffers(kParts);
            const auto flush = [&](std::size_t part) {
                if (buffers[part].empty()) return;
                const std::lock_guard<std::mutex> guard(file_guard[part]);
                std::fwrite(buffers[part].data(), sizeof(Record), buffers[part].size(), files[part]);
                buffers[part].clear();
            };
            std::vector<std::pair<std::uint32_t, std::uint32_t>> proposed;
            std::vector<std::uint32_t> cursor, bound;
            std::vector<std::uint8_t> skip;
            std::vector<std::uint32_t> accepted_list, rejected_list;
            for (std::size_t h = next_haplotype.fetch_add(1); h < unique; h = next_haplotype.fetch_add(1)) {
                const int cls = h < is_panel.size() && is_panel[h] != 0 ? 0 : 1;
                UnionCensusClass& into = local_classes[2 * worker + cls];
                ++into.haplotypes;
                for (int strand = 0; strand < 2; ++strand) {
                    const std::string& text = text_of(h, strand);
                    const auto& forms = strand == 0 ? prepared.forward_forms : prepared.reverse_forms;
                    const auto& indexes = strand == 0 ? prepared.forward_indexes : prepared.reverse_indexes;
                    const std::vector<ReplayZone>* zones = reuse ? &(strand == 0 ? replay.forward : replay.reverse)[h] : nullptr;
                    const bool replaying = zones != nullptr && !zones->empty();
                    if (replaying) {
                        skip.assign(text.size() + g + 2, 0);
                        for (const ReplayZone& zone : *zones) {
                            std::fill(skip.begin() + zone.lo, skip.begin() + zone.hi + 1, static_cast<std::uint8_t>(1));
                        }
                    }
                    proposed.clear();
                    for (const PieceIndex& index : indexes) {   // collect_form_placements' scan, verbatim in effect
                        const std::size_t piece = index.piece_bp;
                        if (text.size() < piece) continue;
                        const std::uint64_t piece_mask =
                            piece == 32 ? std::numeric_limits<std::uint64_t>::max() : ((std::uint64_t{1} << (2 * piece)) - 1);
                        std::uint64_t rolling_code = 0;
                        std::size_t run = 0;
                        for (std::size_t at = 0; at < text.size(); ++at) {
                            const std::uint64_t bits2 = two_bit(text[at]);
                            if (bits2 == kNoCode) { run = 0; rolling_code = 0; continue; }
                            rolling_code = ((rolling_code << 2) | bits2) & piece_mask;
                            if (++run < piece) continue;
                            const std::size_t piece_start = at + 1 - piece;
                            const auto range = std::equal_range(
                                index.entries.begin(), index.entries.end(), rolling_code, [](const auto& lhs, const auto& rhs) {
                                    if constexpr (std::is_same_v<std::decay_t<decltype(lhs)>, PieceEntry>) {
                                        return lhs.code < rhs;
                                    } else {
                                        return lhs < rhs.code;
                                    }
                                });
                            for (auto it = range.first; it != range.second; ++it) {
                                if (it->offset > piece_start) continue;
                                const std::size_t implied = piece_start - it->offset;
                                const std::size_t lo = implied > g ? implied - g : 0;
                                for (std::size_t shifted = lo; shifted <= implied + g; ++shifted) {
                                    if (replaying && shifted < skip.size() && skip[shifted] != 0) continue;
                                    proposed.emplace_back(it->form, static_cast<std::uint32_t>(shifted));
                                }
                            }
                        }
                    }
                    sort_unique_proposals_bucketed(proposed, forms.size(), cursor, bound);
                    const std::vector<std::uint32_t>& ids = window_ids[2 * h + strand];
                    accepted_list.clear();
                    rejected_list.clear();
                    for (std::size_t i = 0; i < proposed.size(); ++i) {
                        const std::uint32_t form_index = proposed[i].first;
                        const std::size_t start = proposed[i].second;
                        ++into.requests;
                        if (start > text.size() || text.size() - start < r) {
                            ++into.ineligible;
                            continue;
                        }
                        const std::size_t window = std::min(text.size() - start, W);
                        const MyersPattern& pattern = forms[form_index].pattern;
                        MyersWord peq[4];
                        for (int c = 0; c < 4; ++c) peq[c] = myers_word_of(pattern.peq_lo[c], pattern.peq_hi[c]);
                        std::uint32_t columns = 0;
                        const bool accepted = myers_shw_run_early_columns(peq, pattern.length, text.data() + start, window,
                                                                          band, columns);
                        (accepted ? accepted_list : rejected_list).push_back(static_cast<std::uint32_t>(i));
                        (accepted ? into.accepted : into.rejected) += 1;
                        (accepted ? into.columns_accepted : into.columns_rejected) += columns;
                        Record record;
                        record.key = (static_cast<std::uint64_t>(read_ids[strand][form_index]) << 27) | ids[start];
                        record.order = static_cast<std::uint32_t>((h << 20) | (static_cast<std::size_t>(strand) << 19) | start);
                        record.columns = static_cast<std::uint16_t>(columns);
                        record.flags = static_cast<std::uint16_t>((accepted ? 1 : 0) | (cls == 0 ? 2 : 0));
                        const std::size_t part = (record.key * 0x9e3779b97f4a7c15ULL) >> 58;
                        buffers[part].push_back(record);
                        if (buffers[part].size() >= 8192) flush(part);
                    }
                    // The production kernel, timed by decision, and checked against the recorded decisions.
                    for (int accepted = 0; accepted < 2; ++accepted) {
                        const std::vector<std::uint32_t>& list = accepted != 0 ? accepted_list : rejected_list;
                        std::size_t agree = 0;
                        const auto timed = Clock::now();
                        for (const std::uint32_t i : list) {
                            const std::size_t start = proposed[i].second;
                            const ReadForm& form = forms[proposed[i].first];
                            agree += feasible_within_band(text, start, std::min(text.size() - start, W), form.bases, band,
                                                          &form.pattern) == (accepted != 0) ? 1 : 0;
                        }
                        (accepted != 0 ? into.prefilter_accepted_ns : into.prefilter_rejected_ns) += static_cast<std::uint64_t>(
                            std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - timed).count());
                        if (agree != list.size()) throw std::logic_error("the census kernel disagrees with the prefilter");
                        local_checked[worker] += list.size();
                    }
                }
            }
            for (std::size_t part = 0; part < kParts; ++part) flush(part);
        } catch (...) {
            failures[worker] = std::current_exception();
        }
    };
    {
        const std::size_t workers = std::max<std::size_t>(1, options.threads);
        std::vector<std::thread> pool;
        for (std::size_t w = 0; w < workers; ++w) pool.emplace_back(work, w);
        for (std::thread& worker : pool) worker.join();
        for (std::FILE* file : files) std::fclose(file);
        for (const std::exception_ptr& failure : failures) {
            if (failure) std::rethrow_exception(failure);
        }
    }
    for (std::size_t w = 0; w < local_classes.size(); ++w) {
        UnionCensusClass& into = (w % 2 == 0) ? census.panel : census.mosaic;
        const UnionCensusClass& from = local_classes[w];
        into.haplotypes += from.haplotypes;
        into.requests += from.requests;
        into.ineligible += from.ineligible;
        into.accepted += from.accepted;
        into.rejected += from.rejected;
        into.columns_accepted += from.columns_accepted;
        into.columns_rejected += from.columns_rejected;
        into.prefilter_accepted_ns += from.prefilter_accepted_ns;
        into.prefilter_rejected_ns += from.prefilter_rejected_ns;
    }
    for (const std::uint64_t checked : local_checked) census.decisions_checked += checked;
    census.seconds_stream = std::chrono::duration<double>(Clock::now() - mark).count();

    // Aggregation, one partition at a time: sort by (identity, order); the least order is the first occurrence.
    mark = Clock::now();
    for (std::size_t part = 0; part < kParts; ++part) {
        const std::string path = scratch_dir + "/part" + std::to_string(part) + ".bin";
        std::vector<Record> records;
        {
            std::FILE* file = std::fopen(path.c_str(), "rb");
            if (file == nullptr) throw std::runtime_error("the census cannot read its scratch partition");
            std::fseek(file, 0, SEEK_END);
            const long size = std::ftell(file);
            std::fseek(file, 0, SEEK_SET);
            records.resize(static_cast<std::size_t>(size) / sizeof(Record));
            if (std::fread(records.data(), sizeof(Record), records.size(), file) != records.size()) {
                throw std::runtime_error("the census read a short partition");
            }
            std::fclose(file);
            std::filesystem::remove(path);
        }
        census.record_bytes += records.size() * sizeof(Record);
        std::sort(records.begin(), records.end(), [](const Record& a, const Record& b) {
            return a.key != b.key ? a.key < b.key : a.order < b.order;
        });
        for (std::size_t i = 0; i < records.size();) {
            std::size_t j = i;
            while (j < records.size() && records[j].key == records[i].key) ++j;
            ++census.identities;
            ++census.multiplicity[j - i];
            ((records[i].flags & 2) != 0 ? census.panel : census.mosaic).identities_first_here += 1;
            for (std::size_t k = i; k < j; ++k) {
                const Record& record = records[k];
                if (record.columns != records[i].columns || (record.flags & 1) != (records[i].flags & 1)) {
                    throw std::logic_error("the census found one identity with two decisions or two costs");
                }
                UnionCensusClass& into = (record.flags & 2) != 0 ? census.panel : census.mosaic;
                const bool accepted = (record.flags & 1) != 0;
                if (k == i) {
                    (accepted ? into.first_accepted : into.first_rejected) += 1;
                } else if (accepted) {
                    ++into.repeated_accepted;
                    into.columns_repeated_accepted += record.columns;
                } else {
                    ++into.repeated_rejected;
                    into.columns_repeated_rejected += record.columns;
                }
            }
            i = j;
        }
    }
    std::filesystem::remove(scratch_dir);
    census.seconds_aggregate = std::chrono::duration<double>(Clock::now() - mark).count();
    return census;
}

// ---------------------------------------------------------------------------------------------
// PreparedScoringContext
// ---------------------------------------------------------------------------------------------

// A persistent pool for window enumeration: size-1 background threads, and the caller as worker 0.
// Created once per context, so no window pays for thread creation. Only one stage uses it at a time --
// the caller is blocked in run() while it works -- so it never runs beside placement or combination.
class EnumerationPool {
public:
    explicit EnumerationPool(std::size_t size) : size_(size) {
        for (std::size_t w = 1; w < size_; ++w) threads_.emplace_back([this, w]() { serve(w); });
    }
    ~EnumerationPool() {
        {
            const std::lock_guard<std::mutex> lock(mutex_);
            stop_ = true;
        }
        start_.notify_all();
        for (std::thread& thread : threads_) thread.join();
    }
    EnumerationPool(const EnumerationPool&) = delete;
    EnumerationPool& operator=(const EnumerationPool&) = delete;
    std::size_t size() const { return size_; }
    // task(w) for w in [0, count); count <= size(). Tasks must not throw: callers capture.
    void run(std::size_t count, const std::function<void(std::size_t)>& task) {
        {
            const std::lock_guard<std::mutex> lock(mutex_);
            task_ = &task;
            count_ = count;
            pending_ = count > 0 ? count - 1 : 0;
            ++generation_;
        }
        start_.notify_all();
        if (count > 0) {
            const ActiveWorker active;
            task(0);
        }
        std::unique_lock<std::mutex> lock(mutex_);
        done_.wait(lock, [this]() { return pending_ == 0; });
        task_ = nullptr;
    }

private:
    void serve(std::size_t w) {
        std::uint64_t seen = 0;
        for (;;) {
            const std::function<void(std::size_t)>* task = nullptr;
            {
                std::unique_lock<std::mutex> lock(mutex_);
                start_.wait(lock, [&]() { return stop_ || generation_ != seen; });
                if (stop_) return;
                seen = generation_;
                if (w >= count_) continue;
                task = task_;
            }
            {
                const ActiveWorker active;
                (*task)(w);
            }
            {
                const std::lock_guard<std::mutex> lock(mutex_);
                --pending_;
            }
            done_.notify_one();
        }
    }
    std::size_t size_ = 1;
    std::vector<std::thread> threads_;
    std::mutex mutex_;
    std::condition_variable start_;
    std::condition_variable done_;
    const std::function<void(std::size_t)>* task_ = nullptr;
    std::size_t count_ = 0;
    std::size_t pending_ = 0;
    std::uint64_t generation_ = 0;
    bool stop_ = false;
};

struct PreparedScoringContext::State {
    std::vector<ReadPair> fragments;
    LibraryModel library;
    FragmentScoringOptions options;
    DiscreteInsertDistribution library_insert;
    PreparedFragmentSet prepared;

    // Materialised bytes exist only for haplotypes that are actually PLACED in full: the initial
    // panel, and whichever repair has been promoted to incumbent. A transient neighbour leaves these
    // entries empty and lives in `segments` instead, so a round costs no 226kb copies, no reversed
    // copies, and no map keyed on whole haplotypes.
    std::deque<std::string> sequences;
    std::deque<std::string> reversed;
    std::map<std::size_t, SegmentedSequence> segments;
    // Dedup by streamed hash, then EXACT comparison on a match. A hash alone cannot decide identity,
    // and byte-identical haplotypes spelled by different allele paths must share one column.
    std::unordered_multimap<std::uint64_t, std::size_t> by_hash;
    std::vector<std::uint64_t> hash_by_slot;

    std::size_t slot_length(std::size_t slot) const {
        const auto found = segments.find(slot);
        return found != segments.end() ? found->second.length() : sequences[slot].size();
    }
    // The window a repair requests. One function for the request AND the plan, so a planned identity
    // can never differ from the one requested.
    struct RepairWindowBounds {
        std::size_t margin = 0;
        std::size_t shortest = 0;
        std::size_t begin = 0;
        std::size_t end = 0;
    };
    RepairWindowBounds repair_window_bounds(std::size_t block_begin, std::size_t block_end,
                                            std::size_t old_length, std::size_t new_length,
                                            std::size_t new_end, std::ptrdiff_t shift) const {
        // The margin beyond which the emission normaliser is saturated, so an unaffected origin keeps
        // its key under a length change. Taken at the LONGEST read, which only ever widens the set.
        std::size_t longest_read = 0;
        for (const ReadPair& fragment : fragments) {
            longest_read = std::max({longest_read, fragment.first.bases.size(), fragment.second.bases.size()});
        }
        RepairWindowBounds bounds;
        bounds.margin = longest_read + library.max_gap_bp;
        const std::size_t max_insert = library_insert_support(library).max_bp;
        const std::size_t halo = max_insert + bounds.margin;
        bounds.shortest = std::min(old_length, new_length);
        bounds.begin = block_begin > halo ? block_begin - halo : 0;
        bounds.end = std::min(new_length, new_end + halo);
        // A length change reaches the ends through the normaliser, so the window must cover whichever
        // end the replacement sits within `margin` of, or the origins there would be dropped by both
        // halves.
        if (shift != 0 && block_begin < bounds.margin) bounds.begin = 0;
        if (shift != 0 && bounds.shortest < block_end + bounds.margin) bounds.end = new_length;
        return bounds;
    }

    bool slot_equals(std::size_t slot, const SegmentedSequence& other) const {
        const auto found = segments.find(slot);
        if (found != segments.end()) return found->second.equals(other);
        return other.equals_materialised(sequences[slot]);
    }
    bool slot_equals_materialised(std::size_t slot, const std::string& other) const {
        const auto found = segments.find(slot);
        if (found != segments.end()) return found->second.equals_materialised(other);
        return sequences[slot] == other;
    }
    // Gives a transient repair real bytes. Done once, when it is promoted to incumbent, because the
    // next round deltas FROM it and full placement needs a string.
    void materialise(std::size_t slot) {
        const auto found = segments.find(slot);
        if (found == segments.end()) return;
        found->second.materialise(sequences[slot]);
        reversed[slot].assign(sequences[slot].rbegin(), sequences[slot].rend());
        segments.erase(found);
    }

    std::vector<std::vector<double>> mass_column;   // [haplotype][fragment]
    std::vector<double> exposure;

    FragmentScoringCounters cumulative;
    EmissionCacheCounters cumulative_cache;
    double cumulative_placement_seconds = 0.0;
    PlacementBudget budget;                          // cumulative across every batch
    bool refused = false;
    bool cache_pinned = false;
    // `transient_watermark` is the slot from which release_transient is allowed to free, so the
    // initial panel columns are structurally out of reach.
    std::size_t transient_watermark = 0;
    std::map<std::size_t, HaploidOriginCensus> censuses;
    // Reused by every delta: rows, groups, accumulators and the flat origin output.
    WindowScratch window_scratch;
    // Size of the matrix the latest score_candidates flattened for its caller, for memory accounts.
    std::uint64_t last_result_bytes = 0;
    // Released repair slots, reused before any new slot is appended. Without this the slot count grew
    // by a whole neighbourhood every round, and score_candidates -- which flattens one column per slot
    // -- grew with it: measured at 26 MiB per round on NA20509 10x, the only accounted item that grew.
    std::vector<std::size_t> free_slots;
    std::size_t slots_created = 0;
    std::uint64_t full_placements = 0;   // I2-D: full haplotype placements performed (measurement and gates)

    // ---- window enumeration across workers. Declared after the scratch and caches it writes into.
    struct EnumerationPartition {
        std::size_t begin = 0;
        std::size_t end = 0;
        std::vector<PieceIndex> forward;
        std::vector<PieceIndex> reverse;
        WindowScratch scratch;
    };
    std::vector<EnumerationPartition> partitions;
    std::unique_ptr<EnumerationPool> pool;
    ParallelEnumerationStats parallel_stats;

    // ---- I2-B (experimental): round-local fine window discovery. Stores live in the enumeration workers' scratch.
    bool fine_enabled = false, fine_round = false;
    std::uint64_t fine_budget = 0, fine_served_cap = 0;
    std::atomic<std::uint64_t> fine_served{0};
    FineWindowReport fine_report;
    std::uint64_t window_enumerations = 0;
    std::size_t longest_read() const {
        std::size_t longest = 0;
        for (const ReadPair& f : fragments) longest = std::max({longest, f.first.bases.size(), f.second.bases.size()});
        return longest;
    }
    void configure_fine(WindowFineStore& store, std::uint64_t budget_share) {
        store.max_read = longest_read();
        store.gap = library.max_gap_bp;
        store.window = store.max_read + 2 * store.gap;
        store.budget = budget_share;
        store.served_cap = fine_served_cap;
        store.served_total = &fine_served;
    }
    std::vector<WindowFineStore*> fine_stores() {
        std::vector<WindowFineStore*> out{&window_scratch.fine};
        for (EnumerationPartition& p : partitions) out.push_back(&p.scratch.fine);
        return out;
    }
    // A round's end (completed or not): tallies folded into the report, every store released.
    void end_fine_round() {
        std::uint64_t round_peak = 0, round_held = 0, round_temporary = 0;
        for (WindowFineStore* st : fine_stores()) {
            round_held += st->peak_held;
            round_temporary += st->peak_temporary;
            st->peak_held = 0;
            st->peak_temporary = 0;
            fine_report.worker_windows += st->windows;
            fine_report.fallback_windows += st->fallback_windows;
            fine_report.units_discovered += st->units_discovered;
            fine_report.units_served += st->units_served;
            fine_report.units_not_retained += st->units_not_retained;
            fine_report.placements_discovered += st->placements_discovered;
            fine_report.placements_served += st->placements_served;
            fine_report.hash_collisions += st->hash_collisions;
            round_peak += st->peak_bytes;
            st->windows = st->fallback_windows = st->units_discovered = st->units_served = st->units_not_retained = 0;
            st->placements_discovered = st->placements_served = st->hash_collisions = st->peak_bytes = 0;
            st->release();
        }
        fine_report.round_peak_bytes = std::max(fine_report.round_peak_bytes, round_peak);
        fine_report.round_peak_held_bytes = std::max(fine_report.round_peak_held_bytes, round_held);
        fine_report.round_peak_temporary_bytes = std::max(fine_report.round_peak_temporary_bytes, round_temporary);
        fine_round = false;
    }

    std::size_t enumeration_workers() const {
        return std::max<std::size_t>(1, std::min(options.threads, fragments.size()));
    }

    // Every production window enumeration goes through here. One worker: the reference function
    // itself, untouched. Several: contiguous fragment ranges, each with its own piece indexes, scratch
    // and counters, merged in ascending fragment order into `out` -- the same origins, offsets and work
    // counters one thread produces.
    void enumerate_window(const LocalWindow& window, PlacementBudget& shared_budget,
                          FragmentScoringCounters& counters, EmissionCacheCounters& cache_work,
                          WindowScratch& out);

    // A slot for a new repair: a released one if any, else a new one. A reused slot's column is reset
    // to -inf, which is what a fresh slot starts at and what stage 5 relies on for a fragment with no
    // states. Its sequence strings are already empty (released), its census, segment and hash entry
    // already gone, and anything keyed by its old identity as an incumbent was dropped when it was
    // released -- so no retained structure can mistake the new haplotype for the old one.
    std::size_t allocate_repair_slot(std::uint64_t hash) {
        std::size_t slot = 0;
        if (!free_slots.empty()) {
            slot = free_slots.back();
            free_slots.pop_back();
            hash_by_slot[slot] = hash;
            mass_column[slot].assign(fragments.size(), kNegInf);
        } else {
            slot = sequences.size();
            hash_by_slot.push_back(hash);
            sequences.emplace_back();
            reversed.emplace_back();
            mass_column.resize(slot + 1, std::vector<double>(fragments.size(), kNegInf));
            exposure.resize(slot + 1, 0.0);
        }
        by_hash.emplace(hash, slot);
        ++slots_created;
        return slot;
    }

    // Stage 6, off unless enabled. The cache is content-keyed and survives release_transient; the
    // group preaggregates and full aggregates are working state keyed by incumbent and are dropped
    // with it.
    bool origin_path = false;
    std::size_t origin_cache_budget = 0;
    std::uint64_t origin_record_cap = 0;
    std::uint64_t origin_records_consumed = 0;
    OriginCache origin_cache;
    // Planned retention: set by plan_origin_cache, indexed by demand-ledger row. Once planned, a miss is
    // retained only if its row was selected by the latest plan.
    bool origin_planned = false;
    std::vector<char> plan_selected;
    std::map<OriginIdentity, std::uint64_t> plan_multiplicity;
    // The group holder: the identity whose records miss_records holds from the latest enumeration.
    bool holder_valid = false;
    std::string holder_bytes;
    bool holder_at_start = false;
    bool holder_at_end = false;
    std::uint32_t holder_rel_block_begin = 0;
    std::uint32_t holder_rel_new_end = 0;
    // The demand ledger: every distinct window requested, admitted or not. Same identity as the cache.
    std::unordered_multimap<std::uint64_t, std::size_t> demand_index;
    std::deque<std::pair<std::string, OriginDemandRow>> demand;
    std::map<GroupKey, GroupEntry> groups;
    std::map<std::pair<std::size_t, std::uint32_t>, std::vector<StateCount>> full_aggregates;
    OriginRecordSet miss_records;
    OriginRecordSet check_records;
    StateAggregateSet candidate_states;
    std::vector<StateKey> key_scratch;
    std::vector<StateCount> merged;

    void require_usable() const {
        if (refused) {
            throw ResourceExhausted(
                "this scoring context refused a declared cap; the whole search is invalid and no "
                "partial ranking is available from it");
        }
    }
};

void PreparedScoringContext::State::enumerate_window(const LocalWindow& window, PlacementBudget& shared_budget,
                                                     FragmentScoringCounters& counters,
                                                     EmissionCacheCounters& cache_work, WindowScratch& out) {
    const std::size_t workers = enumeration_workers();
    const std::size_t band = options.mismatch_band;
    parallel_stats.workers = workers;
    ++window_enumerations;
    // I2-B: a store is active only for the windows of an open round, and only during this call
    struct Deactivate {
        std::vector<WindowFineStore*> stores;
        ~Deactivate() {
            for (WindowFineStore* s : stores) s->active = false;
        }
    } deactivate;
    if (workers <= 1) {
        const ActiveWorker active;
        if (fine_round) {
            if (out.fine.budget == 0) configure_fine(out.fine, fine_budget);
            out.fine.active = true;
            deactivate.stores.push_back(&out.fine);
        }
        enumerate_window_origins(prepared, fragments, library, options, band, window, shared_budget, counters,
                                 cache_work, out);
        return;
    }
    const std::size_t n = fragments.size();
    if (partitions.size() != workers) {
        partitions.clear();
        partitions.resize(workers);
        const auto restrict = [&](const std::vector<PieceIndex>& all, const std::deque<ReadForm>& forms,
                                  std::size_t begin, std::size_t end) {
            std::vector<PieceIndex> part;
            for (const PieceIndex& index : all) {
                PieceIndex sub;
                sub.piece_bp = index.piece_bp;
                for (const PieceEntry& entry : index.entries) {
                    const std::size_t fragment = forms[entry.form].fragment;
                    if (fragment >= begin && fragment < end) sub.entries.push_back(entry);
                }
                part.push_back(std::move(sub));
            }
            return part;
        };
        for (std::size_t w = 0; w < workers; ++w) {
            EnumerationPartition& p = partitions[w];
            p.begin = n * w / workers;
            p.end = n * (w + 1) / workers;
            p.forward = restrict(prepared.forward_indexes, prepared.forward_forms, p.begin, p.end);
            p.reverse = restrict(prepared.reverse_indexes, prepared.reverse_forms, p.begin, p.end);
        }
    }
    if (!pool || pool->size() != workers) pool = std::make_unique<EnumerationPool>(workers);
    if (fine_round) {
        for (EnumerationPartition& p : partitions) {
            if (p.scratch.fine.budget == 0) configure_fine(p.scratch.fine, fine_budget / workers);
            p.scratch.fine.active = true;
            deactivate.stores.push_back(&p.scratch.fine);
        }
    }

    struct Outcome {
        FragmentScoringCounters counters;
        EmissionCacheCounters cache;
        std::uint64_t proposed[2] = {0, 0};
        DeferredEmissions deferred;
        std::exception_ptr failure;
        bool refusal = false;
    };
    std::vector<Outcome> outcomes(workers);
    const bool caching = emission_cache().active.load(std::memory_order_relaxed);
    pool->run(workers, [&](std::size_t w) {
        EnumerationPartition& p = partitions[w];
        Outcome& o = outcomes[w];
        deferred_emissions() = caching ? &o.deferred : nullptr;
        try {
            enumerate_window_origins(prepared, fragments, library, options, band, window, shared_budget,
                                     o.counters, o.cache, p.scratch, &p.forward, &p.reverse, p.begin, p.end,
                                     o.proposed);
        } catch (const PlacementCancelled&) {
            // Another worker crossed a cap; its refusal is the one reported.
        } catch (const ResourceExhausted&) {
            o.refusal = true;
            o.failure = std::current_exception();
        } catch (...) {
            o.failure = std::current_exception();
        }
        deferred_emissions() = nullptr;
    });
    // The worker that crossed a cap is reported, never a cancellation it caused.
    for (const Outcome& o : outcomes) {
        if (o.refusal) std::rethrow_exception(o.failure);
    }
    for (const Outcome& o : outcomes) {
        if (o.failure) std::rethrow_exception(o.failure);
    }
    if (shared_budget.stop()) {
        throw std::logic_error("a window enumeration was cancelled without a refusal");
    }
    // The per-haplotype proposed-start cap bounds one PASS over the window; each worker saw only its
    // share, so it is decided exactly here, on the merged totals.
    for (int pass = 0; pass < 2; ++pass) {
        std::uint64_t proposed = 0;
        for (const Outcome& o : outcomes) proposed += o.proposed[pass];
        if (proposed > options.max_proposed_starts_per_haplotype) {
            shared_budget.cancel();
            throw ResourceExhausted(
                "fragment scoring reached its per-haplotype proposed-start cap; no partial ranking is "
                "returned");
        }
    }
    // Origins and offsets, in ascending fragment order.
    out.origins.clear();
    out.offsets.assign(n + 1, 0);
    for (std::size_t w = 0; w < workers; ++w) {
        const EnumerationPartition& p = partitions[w];
        for (std::size_t f = p.begin; f < p.end; ++f) {
            out.offsets[f] = out.origins.size();
            out.origins.insert(out.origins.end(),
                               p.scratch.origins.begin() + static_cast<std::ptrdiff_t>(p.scratch.offsets[f]),
                               p.scratch.origins.begin() + static_cast<std::ptrdiff_t>(p.scratch.offsets[f + 1]));
        }
    }
    out.offsets[n] = out.origins.size();
    // Work counters: sums, except the window scan every worker repeats -- counted once.
    for (std::size_t w = 0; w < workers; ++w) {
        const FragmentScoringCounters& c = outcomes[w].counters;
        if (w == 0) {
            counters.candidate_bases_scanned += c.candidate_bases_scanned;
        } else {
            parallel_stats.duplicate_scan_bases += c.candidate_bases_scanned;
        }
        counters.piece_hits += c.piece_hits;
        counters.proposed_starts += c.proposed_starts;
        counters.unique_starts += c.unique_starts;
        counters.full_read_verifications += c.full_read_verifications;
        counters.dp_evaluations += c.dp_evaluations;
        counters.prefilter_rejections += c.prefilter_rejections;
        counters.accepted_mate_placements += c.accepted_mate_placements;
        counters.valid_fr_states += c.valid_fr_states;
        counters.latent_alignment_pairs += c.latent_alignment_pairs;
        counters.physical_origins += c.physical_origins;
        counters.statistical_state_classes += c.statistical_state_classes;
    }
    // The emission cache, replayed in one thread's order: the forward pass of every range, then the
    // reverse pass of every range. A request the frozen cache already answered is a hit either way.
    std::uint64_t shared_hits = 0;
    std::uint64_t executed = 0;
    for (const Outcome& o : outcomes) {
        cache_work.executed_prefilters += o.cache.executed_prefilters;
        cache_work.dp_cache_requests += o.cache.dp_cache_requests;
        shared_hits += o.cache.dp_cache_hits;
        executed += o.cache.executed_dps;
    }
    std::uint64_t replay_hits = 0;
    std::uint64_t replay_misses = 0;
    if (caching) {
        const char* forward_begin = window.bases.data();
        const char* forward_end = forward_begin + window.bases.size();
        for (int pass = 0; pass < 2; ++pass) {
            for (const Outcome& o : outcomes) {
                for (const DeferredEmissions::Record& record : o.deferred.records) {
                    const bool forward = record.key.window >= forward_begin && record.key.window < forward_end;
                    if (forward != (pass == 0)) continue;
                    if (cache_contains(record.key)) {
                        ++replay_hits;
                        continue;
                    }
                    ++replay_misses;
                    if (record.refused) {
                        cache_store_refusal(record.key, CacheOutcome::RefusedAfterRecursion);
                    } else {
                        cache_store_fit(record.key, record.fit);
                    }
                }
            }
        }
        cache_work.dp_cache_hits += shared_hits + replay_hits;
        cache_work.dp_cache_misses += replay_misses;
        cache_work.executed_dps += replay_misses;
        parallel_stats.duplicate_dps += executed - replay_misses;
    } else {
        cache_work.dp_cache_hits += shared_hits;
        cache_work.executed_dps += executed;
    }
    ++parallel_stats.parallel_enumerations;
}

PreparedScoringContext::PreparedScoringContext(std::vector<ReadPair> fragments,
                                               LibraryModel library,
                                               FragmentScoringOptions options)
    : state_(std::make_unique<State>()) {
    validate_library_model(library);
    if (options.minimum_piece_bp == 0) {
        throw std::invalid_argument("minimum indexable piece length must be positive");
    }
    if (options.max_piece_hits == 0 || options.max_proposed_starts_per_haplotype == 0 ||
        options.max_physical_origins == 0 || options.max_state_classes == 0 ||
        options.max_mass_cells == 0 || options.max_latent_alignment_pairs == 0) {
        throw std::invalid_argument("every fragment-scoring resource cap must be positive");
    }
    if (options.record_origin_summary || !options.record_origins_for_fragments.empty()) {
        // Origin recording indexes by a stride that changes every time a haplotype is added, and it
        // forces the serial path. Refusing is honest; half-supporting it would give a diagnostic
        // whose rows silently referred to a different haplotype count than the one that produced them.
        throw std::invalid_argument(
            "origin recording is not available from an incremental scoring context");
    }
    state_->fragments = std::move(fragments);
    state_->library = std::move(library);
    state_->options = std::move(options);
    state_->library_insert = build_library_insert_distribution(state_->library);
    state_->prepared =
        prepare_fragment_set(state_->fragments, state_->library, state_->options);
    state_->cumulative.fragments_with_floor_binding = state_->prepared.fragments_with_floor_binding;
    state_->budget.max_piece_hits = state_->options.max_piece_hits;
    state_->budget.max_physical_origins = state_->options.max_physical_origins;
    state_->budget.publish_batch =
        std::max<std::uint64_t>(1, state_->options.piece_hit_publish_batch);

    // Armed once, for the whole life of the context. This is the only reason a cache key may outlive
    // the batch that created it: every buffer it points into is owned here and never moves.
    if (scoring_detail::emission_cache_pinned().exchange(true, std::memory_order_relaxed)) {
        throw std::runtime_error("another PreparedScoringContext already holds the emission cache");
    }
    state_->cache_pinned = true;
    scoring_detail::emission_cache_pinned().store(false, std::memory_order_relaxed);
    begin_emission_cache(state_->options.emission_cache_bytes);
    // I2-K (experimental): the union pass's emission results, under an exactly compatible context only
    (void)adopt_emission_cache(state_->options.emission_handoff_in, state_->library, state_->options.mismatch_band,
                         state_->options.gap_recursion_backend, state_->options.emission_cache_bytes);
    state_->options.emission_handoff_in = nullptr;
    scoring_detail::emission_cache_pinned().store(true, std::memory_order_relaxed);
}

PreparedScoringContext::~PreparedScoringContext() {
    if (state_ && state_->cache_pinned) {
        scoring_detail::emission_cache_pinned().store(false, std::memory_order_relaxed);
        state_->cumulative_cache = end_emission_cache();
    }
}

std::uint64_t PreparedScoringContext::full_placement_count() const { return state_->full_placements; }

std::vector<std::size_t> PreparedScoringContext::add_haplotypes(
    const std::vector<std::string>& sequences,
    IncrementalBatchReport* report,
    bool reuse_released_slots,
    bool retain_censuses) {
    State& state = *state_;
    state.require_usable();
    const std::size_t before = state.sequences.size();

    std::vector<std::size_t> out;
    out.reserve(sequences.size());
    std::vector<std::size_t> placed;   // the slots this batch places, appended or reused
    std::size_t reused = 0;
    for (const std::string& sequence : sequences) {
        const std::uint64_t hash = hash_materialised(sequence);
        std::size_t existing = state.sequences.size();
        const auto range = state.by_hash.equal_range(hash);
        for (auto it = range.first; it != range.second; ++it) {
            if (state.slot_equals_materialised(it->second, sequence)) {
                existing = it->second;
                break;
            }
        }
        if (existing != state.sequences.size()) {
            ++reused;
            out.push_back(existing);
            continue;
        }
        std::size_t slot = state.sequences.size();
        if (reuse_released_slots && !state.free_slots.empty()) {
            // A released slot owns nothing (see allocate_repair_slot): its column restarts at -inf.
            slot = state.free_slots.back();
            state.free_slots.pop_back();
            state.hash_by_slot[slot] = hash;
            state.sequences[slot] = sequence;
            state.reversed[slot].assign(sequence.rbegin(), sequence.rend());
            state.mass_column[slot].assign(state.fragments.size(), kNegInf);
        } else {
            state.hash_by_slot.push_back(hash);
            state.sequences.push_back(sequence);
            state.reversed.emplace_back(sequence.rbegin(), sequence.rend());
        }
        state.by_hash.emplace(hash, slot);
        ++state.slots_created;
        placed.push_back(slot);
        out.push_back(slot);
    }
    const std::size_t after = state.sequences.size();
    const std::size_t added = placed.size();
    const bool contiguous = placed.size() == after - before;   // nothing reused: exactly [before, after)

    if (!state.fragments.empty() &&
        after > state.options.max_mass_cells / state.fragments.size()) {
        state.refused = true;
        throw ResourceExhausted(
            "the fragment-by-haplotype mass matrix would exceed its declared cell cap; no partial "
            "ranking is returned");
    }

    state.mass_column.resize(after, std::vector<double>(state.fragments.size(), kNegInf));
    state.exposure.resize(after, 0.0);
    for (const std::size_t h : placed) {
        state.exposure[h] =
            exposure_for_length(state.sequences[h].size(), state.fragments.size(),
                                state.prepared.fragments_by_floor, state.library_insert);
    }

    FragmentScoringCounters batch;
    EmissionCacheCounters batch_cache;
    ScorerMemoryAccounts batch_memory;
    std::size_t threads_used = 1;
    // I2-D: censuses of the slots this batch places, filled by the same placement pass; owned by the context only once
    // the batch's work has been accepted
    std::map<std::size_t, HaploidOriginCensus> fresh_censuses;
    CensusCapture capture;
    if (retain_censuses) {
        for (const std::size_t h : placed) capture.wanted.emplace(h, &fresh_censuses[h]);
    }
    const auto mark = std::chrono::steady_clock::now();
    if (added != 0) {
        PlacementInputs inputs;
        inputs.prepared = &state.prepared;
        inputs.fragments = &state.fragments;
        inputs.library = &state.library;
        inputs.library_insert = &state.library_insert;
        inputs.options = &state.options;
        inputs.sequences = &state.sequences;
        inputs.reversed = &state.reversed;
        inputs.band = state.options.mismatch_band;
        if (!contiguous) inputs.slot_list = &placed;
        OriginRecording none;
        try {
            threads_used = run_placement(inputs, contiguous ? before : 0, contiguous ? after : placed.size(),
                                         state.options.threads, false, state.budget, state.mass_column, batch,
                                         batch_memory, batch_cache, none, retain_censuses ? &capture : nullptr);
        } catch (const ResourceExhausted&) {
            // A refusal invalidates the WHOLE search, not this batch: columns already retained were
            // computed under a budget that has now been exceeded, and ranking from them would present
            // a partial search as a complete one.
            state.refused = true;
            throw;
        }
    }
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - mark).count();
    state.full_placements += added;

    // Cumulative first, because the caps bound the whole search rather than any one batch.
    state.cumulative.candidate_bases_scanned += batch.candidate_bases_scanned;
    state.cumulative.piece_hits += batch.piece_hits;
    state.cumulative.proposed_starts += batch.proposed_starts;
    state.cumulative.unique_starts += batch.unique_starts;
    state.cumulative.full_read_verifications += batch.full_read_verifications;
    state.cumulative.dp_evaluations += batch.dp_evaluations;
    state.cumulative.prefilter_rejections += batch.prefilter_rejections;
    state.cumulative.accepted_mate_placements += batch.accepted_mate_placements;
    state.cumulative.valid_fr_states += batch.valid_fr_states;
    state.cumulative.latent_alignment_pairs += batch.latent_alignment_pairs;
    state.cumulative.physical_origins += batch.physical_origins;
    state.cumulative.statistical_state_classes += batch.statistical_state_classes;
    state.cumulative_cache.executed_prefilters += batch_cache.executed_prefilters;
    state.cumulative_cache.dp_cache_requests += batch_cache.dp_cache_requests;
    state.cumulative_cache.dp_cache_hits += batch_cache.dp_cache_hits;
    state.cumulative_cache.dp_cache_misses += batch_cache.dp_cache_misses;
    state.cumulative_cache.executed_dps += batch_cache.executed_dps;
    state.cumulative_placement_seconds += seconds;

    try {
        enforce_acceptance_caps(state.cumulative, state.options, []() {});
    } catch (const ResourceExhausted&) {
        state.refused = true;
        throw;
    }
    if (retain_censuses) {
        // I2-D: a slot this batch did not place (already held) is censused the existing way unless it already holds one.
        // I2-E fix: those fallback placements are part of the batch, so they run BEFORE the captured censuses are
        // published; if one is refused, every census this batch produced -- captured or fallback -- is withdrawn, and the
        // context stays refused (retain_census marked it). Nothing from a failed batch is ever owned by the context.
        std::vector<std::size_t> fallback_done;
        try {
            for (const std::size_t h : out) {
                if (state.censuses.count(h) != 0 || fresh_censuses.count(h) != 0) continue;
                if (std::find(fallback_done.begin(), fallback_done.end(), h) != fallback_done.end()) continue;
                retain_census(h);
                fallback_done.push_back(h);
            }
        } catch (...) {
            for (const std::size_t h : fallback_done) state.censuses.erase(h);
            state.refused = true;
            throw;
        }
        for (auto& [h, census] : fresh_censuses) state.censuses.emplace(h, std::move(census));
    }

    if (report != nullptr) {
        report->haplotypes_offered = sequences.size();
        report->haplotypes_new = added;
        report->haplotypes_reused = reused;
        report->unique_haplotypes_total = after;
        report->counters = batch;
        report->emission_cache = batch_cache;
        report->seconds_placement_and_mass = seconds;
        report->threads_used_placement = threads_used;
        report->memory = batch_memory;
        const EmissionCacheCounters stored = read_emission_cache_counters();
        report->emission_cache.payload_bytes = stored.payload_bytes;
        report->emission_cache.entries = stored.entries;
        report->emission_cache.entries_not_retained = stored.entries_not_retained;
        report->cumulative_counters = state.cumulative;
        report->cumulative_emission_cache = state.cumulative_cache;
        report->cumulative_emission_cache.payload_bytes = stored.payload_bytes;
        report->cumulative_emission_cache.entries = stored.entries;
        report->cumulative_emission_cache.entries_not_retained = stored.entries_not_retained;
        report->cumulative_seconds_placement_and_mass = state.cumulative_placement_seconds;
    }
    return out;
}

FragmentScoringResult PreparedScoringContext::score_candidates(
    const std::vector<HaplotypeCandidate>& candidates) const {
    State& state = *state_;
    state.require_usable();
    const std::size_t unique = state.sequences.size();
    for (const HaplotypeCandidate& candidate : candidates) {
        if (candidate.first >= unique || candidate.second >= unique) {
            throw std::invalid_argument("diplotype candidate names an unknown haplotype");
        }
    }
    // A released haplotype has no column. Naming one is a caller error, and it must be refused here
    // rather than discovered as a read through freed storage further down.
    for (const HaplotypeCandidate& candidate : candidates) {
        if (state.mass_column[candidate.first].empty() || state.mass_column[candidate.second].empty()) {
            throw std::invalid_argument(
                "diplotype candidate names a haplotype whose storage has been released");
        }
    }
    FragmentScoringResult result;
    result.scores.resize(candidates.size());
    result.unique_haplotypes = unique;
    result.counters = state.cumulative;
    result.emission_cache = state.cumulative_cache;

    const std::vector<ReadPair>& fragments = state.fragments;
    const auto combination_mark = std::chrono::steady_clock::now();
    const std::vector<AuditCombination> audits = audit_combinations(state.options.audit_libraries);
    // I2-F: the run's library prepared once (validated by audit_combinations), or the scalar reference per term
    const std::vector<AuditCombination> main_terms = state.options.prepared_library_terms && !candidates.empty()
                                                         ? audit_combinations({state.library})
                                                         : std::vector<AuditCombination>{};
    const AuditCombination* main = main_terms.empty() ? nullptr : &main_terms.front();
    result.audit_scores.assign(audits.size(), std::vector<DiplotypeScore>(candidates.size()));
    const auto combine_range = [&](std::size_t begin, std::size_t end) {
        for (std::size_t c = begin; c < end && !audits.empty(); ++c) {
            const std::size_t a = candidates[c].first;
            const std::size_t b = candidates[c].second;
            combine_with_audits(state.mass_column[a], state.mass_column[b], state.prepared.background_by_fragment,
                                fragments.size(), state.exposure[a] + state.exposure[b], state.library, audits,
                                result.scores[c], result.audit_scores, c, main);
        }
        for (std::size_t c = begin; c < end && audits.empty(); ++c) {
            const std::size_t a = candidates[c].first;
            const std::size_t b = candidates[c].second;
            DiplotypeScore& score = result.scores[c];
            score.exposure_term =
                main != nullptr ? -(main->exposure_scale * (state.exposure[a] + state.exposure[b]))
                                : -expected_signal_fragments(state.exposure[a] + state.exposure[b], state.library);
            score.absolute_term_sum = std::abs(score.exposure_term);
            double fragment_term = 0.0;
            for (std::size_t f = 0; f < fragments.size(); ++f) {
                const double combined =
                    log_add_exp(state.mass_column[a][f], state.mass_column[b][f]);
                if (combined != kNegInf) ++score.fragments_with_mass;
                const double contribution =
                    main != nullptr ? prepared_contribution(*main, combined, state.prepared.background_by_fragment[f])
                                    : log_fragment_contribution(combined, state.prepared.background_by_fragment[f],
                                                                state.library);
                fragment_term += contribution;
                score.absolute_term_sum += std::abs(contribution);
            }
            score.fragment_term = fragment_term;
            score.log_likelihood = score.exposure_term + fragment_term;
        }
    };
    const std::size_t threads =
        candidates.empty()
            ? 1
            : std::max<std::size_t>(1, std::min(state.options.threads, candidates.size()));
    result.threads_used_combination = threads;
    result.candidate_fragment_evaluations =
        static_cast<std::uint64_t>(candidates.size()) * static_cast<std::uint64_t>(fragments.size());
    if (threads <= 1) {
        const ActiveWorker active;
        combine_range(0, candidates.size());
    } else {
        const std::size_t span = (candidates.size() + threads - 1) / threads;
        std::vector<std::exception_ptr> failures(threads);
        std::vector<std::thread> workers;
        workers.reserve(threads);
        for (std::size_t t = 0; t < threads; ++t) {
            const std::size_t begin = std::min(t * span, candidates.size());
            const std::size_t end = std::min(begin + span, candidates.size());
            workers.emplace_back([&, t, begin, end]() {
                const ActiveWorker active;
                try {
                    combine_range(begin, end);
                } catch (...) {
                    failures[t] = std::current_exception();
                }
            });
        }
        for (std::thread& worker : workers) worker.join();
        for (const std::exception_ptr& failure : failures) {
            if (failure) std::rethrow_exception(failure);
        }
    }
    if (main != nullptr) {
        std::uint64_t with_mass = 0;
        for (const DiplotypeScore& score : result.scores) with_mass += score.fragments_with_mass;
        note_prepared_combination(candidates.size(), fragments.size(), with_mass, *main);
    }
    result.seconds_diploid_combination =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - combination_mark).count();
    result.seconds_placement_and_mass = state.cumulative_placement_seconds;

    result.log_mass.assign(fragments.size() * unique, kNegInf);
    for (std::size_t h = 0; h < unique; ++h) {
        // Released slots keep their index so nothing is renumbered, but hold no column. Their cells
        // stay at -inf: the haplotype no longer exists, and walking into the empty vector is what
        // turned the first release into a segfault.
        if (state.mass_column[h].empty()) continue;
        for (std::size_t f = 0; f < fragments.size(); ++f) {
            result.log_mass[f * unique + h] = state.mass_column[h][f];
        }
    }
    result.exposure = state.exposure;
    result.memory.mass_matrix_bytes =
        static_cast<std::uint64_t>(result.log_mass.size()) * sizeof(double);
    state.last_result_bytes =
        static_cast<std::uint64_t>(result.log_mass.capacity() + result.exposure.capacity()) *
            sizeof(double) +
        static_cast<std::uint64_t>(result.scores.capacity()) * sizeof(DiplotypeScore);
    return result;
}

std::size_t PreparedScoringContext::unique_haplotypes() const { return state_->sequences.size(); }

std::size_t PreparedScoringContext::haplotypes_created() const { return state_->slots_created; }

double PreparedScoringContext::mass_of(std::size_t fragment, std::size_t haplotype) const {
    state_->require_usable();
    return state_->mass_column.at(haplotype).at(fragment);
}

double PreparedScoringContext::exposure_of(std::size_t haplotype) const {
    state_->require_usable();
    return state_->exposure.at(haplotype);
}

const FragmentScoringCounters& PreparedScoringContext::cumulative_counters() const {
    return state_->cumulative;
}

const EmissionCacheCounters& PreparedScoringContext::cumulative_emission_cache() const {
    return state_->cumulative_cache;
}

bool PreparedScoringContext::refused() const { return state_->refused; }

void PreparedScoringContext::enable_fine_window_discovery(std::size_t budget_bytes, std::uint64_t max_served_placements) {
    if (budget_bytes == 0 || max_served_placements == 0) {
        throw std::invalid_argument("fine window discovery needs a positive store budget and a positive served-placement cap");
    }
    State& state = *state_;
    state.fine_enabled = true;
    state.fine_budget = budget_bytes;
    state.fine_served_cap = max_served_placements;
    state.fine_report.enabled = true;
    state.fine_report.budget_bytes = budget_bytes;
    state.fine_report.served_cap = max_served_placements;
}

FineWindowReport PreparedScoringContext::fine_window_report() const {
    const State& state = *state_;
    FineWindowReport r = state.fine_report;
    r.enumerations = state.window_enumerations;
    r.seconds_rows = state.window_scratch.seconds_rows;
    r.seconds_pairing = state.window_scratch.seconds_pairing;
    for (const auto& p : state.partitions) {
        r.seconds_rows += p.scratch.seconds_rows;
        r.seconds_pairing += p.scratch.seconds_pairing;
    }
    return r;
}

void PreparedScoringContext::enable_origin_cache_path(bool enabled) {
    state_->require_usable();
    state_->origin_path = enabled;
}

void PreparedScoringContext::set_origin_cache_limits(std::size_t budget_bytes,
                                                     std::uint64_t max_records_consumed) {
    state_->require_usable();
    state_->origin_cache_budget = budget_bytes;
    state_->origin_record_cap = max_records_consumed;
}

std::vector<OriginDemandRow> PreparedScoringContext::origin_demand() const {
    std::vector<OriginDemandRow> out;
    out.reserve(state_->demand.size());
    for (const auto& row : state_->demand) out.push_back(row.second);
    return out;
}

PlannedRepairWindow PreparedScoringContext::plan_repaired_window(std::size_t incumbent,
                                                                std::size_t block_begin,
                                                                std::size_t block_end,
                                                                const std::string& replacement) const {
    const State& state = *state_;
    if (incumbent >= state.sequences.size() || state.sequences[incumbent].empty()) {
        throw std::invalid_argument("the incumbent haplotype has not been materialised");
    }
    const std::string& old_sequence = state.sequences[incumbent];
    if (block_end < block_begin || block_end > old_sequence.size()) {
        throw std::invalid_argument("the replaced block is not an interval of the incumbent");
    }
    SegmentedSequence repaired;
    repaired.base = &old_sequence;
    repaired.block_begin = block_begin;
    repaired.block_end = block_end;
    repaired.replacement = replacement;
    PlannedRepairWindow out;
    out.haplotype_hash = repaired.hash();
    const auto range = state.by_hash.equal_range(out.haplotype_hash);
    for (auto it = range.first; it != range.second && !out.already_held; ++it) {
        out.already_held = state.slot_equals(it->second, repaired);
    }
    const std::size_t new_length = repaired.length();
    const std::size_t new_end = block_begin + replacement.size();
    const std::ptrdiff_t shift = static_cast<std::ptrdiff_t>(replacement.size()) -
                                 static_cast<std::ptrdiff_t>(block_end - block_begin);
    const State::RepairWindowBounds bounds = state.repair_window_bounds(
        block_begin, block_end, old_sequence.size(), new_length, new_end, shift);
    repaired.append_range(bounds.begin, bounds.end, out.bytes);
    out.at_start = bounds.begin == 0;
    out.at_end = bounds.end == new_length;
    out.rel_block_begin = static_cast<std::uint32_t>(block_begin - bounds.begin);
    out.rel_new_end = static_cast<std::uint32_t>(new_end - bounds.begin);
    out.key_hash = origin_key_hash(out.bytes, out.at_start, out.at_end, out.rel_block_begin, out.rel_new_end);
    return out;
}

OriginCachePlanReport PreparedScoringContext::plan_origin_cache(
    const std::vector<PlannedRepairWindow>& windows) {
    State& state = *state_;
    state.require_usable();
    OriginCachePlanReport report;
    report.planned_repairs = windows.size();
    // A window's multiplicity is the number of DISTINCT new haplotypes requesting it: a haplotype
    // offered twice is placed once, and a held one requests nothing.
    using Identity = OriginIdentity;
    std::map<Identity, std::set<std::uint64_t>> requesters;
    for (const PlannedRepairWindow& w : windows) {
        if (w.already_held) continue;
        requesters[Identity{w.key_hash, w.bytes, w.at_start, w.at_end, w.rel_block_begin, w.rel_new_end}]
            .insert(w.haplotype_hash);
    }
    report.distinct_windows = requesters.size();
    state.plan_multiplicity.clear();
    for (const auto& item : requesters) state.plan_multiplicity.emplace(item.first, item.second.size());
    state.holder_valid = false;
    struct Candidate {
        const Identity* identity = nullptr;
        std::uint64_t multiplicity = 0;
        std::uint64_t bytes = 0;
        std::uint64_t row = 0;
    };
    std::vector<Candidate> candidates;
    for (const auto& item : requesters) {
        report.window_requests += item.second.size();
        const Identity& id = item.first;
        const auto range = state.demand_index.equal_range(std::get<0>(id));
        for (auto it = range.first; it != range.second; ++it) {
            const auto& row = state.demand[it->second];
            if (row.second.at_start == std::get<2>(id) && row.second.at_end == std::get<3>(id) &&
                row.second.rel_block_begin == std::get<4>(id) && row.second.rel_new_end == std::get<5>(id) &&
                row.first == std::get<1>(id)) {
                // Known only once enumerated: entry_bytes is written at a miss, identical at every miss.
                if (row.second.misses != 0) {
                    candidates.push_back({&id, item.second.size(), row.second.entry_bytes, row.second.order});
                }
                break;
            }
        }
    }
    report.known_windows = candidates.size();
    // Priority m/bytes, descending, compared exactly as m_a * bytes_b > m_b * bytes_a; then key hash,
    // then the full key. A total order over identities, so the selection is a function of the SET of
    // planned windows and never of the order they were offered in.
    std::sort(candidates.begin(), candidates.end(), [](const Candidate& a, const Candidate& b) {
        return ranks_before(a.multiplicity, a.bytes, *a.identity, b.multiplicity, b.bytes, *b.identity);
    });
    state.plan_selected.assign(state.demand.size(), 0);
    std::uint64_t used = 0;
    for (const Candidate& c : candidates) {
        if (used + c.bytes > state.origin_cache_budget) continue;
        used += c.bytes;
        state.plan_selected[c.row] = 1;
        ++report.selected_entries;
    }
    report.selected_bytes = used;
    state.origin_planned = true;
    // Evict everything outside the selection, releasing its storage now: the cache may only ever hold
    // selected windows, whose total fits the budget.
    OriginCache& cache = state.origin_cache;
    for (std::size_t at = 0; at < cache.entries.size(); ++at) {
        if (cache.live[at] == 0) continue;
        if (state.plan_selected[cache.demand_of[at]] != 0) {
            ++report.retained_entries;
            const auto planned = state.plan_multiplicity.find(entry_identity(cache.entries[at]));
            cache.multiplicity[at] = planned != state.plan_multiplicity.end() ? planned->second : 0;
            continue;
        }
        ++state.demand[cache.demand_of[at]].second.evictions;
        report.evicted_bytes += evict_origin_slot(cache, at);
        ++report.evicted_entries;
    }
    return report;
}

std::vector<std::size_t> PreparedScoringContext::add_repaired_haplotypes(
    const std::vector<RepairRequest>& requests, std::vector<BlockDeltaWork>* works,
    OriginCachePlanReport* plan) {
    State& state = *state_;
    state.require_usable();
    std::vector<PlannedRepairWindow> planned;
    planned.reserve(requests.size());
    for (const RepairRequest& r : requests) {
        planned.push_back(plan_repaired_window(r.incumbent, r.block_begin, r.block_end, r.replacement));
    }
    if (state.origin_path) {
        const OriginCachePlanReport report = plan_origin_cache(planned);
        if (plan != nullptr) *plan = report;
    }
    // The canonical order: held haplotypes (no window) first, then by window identity so each window's
    // requests are adjacent, then by haplotype and by the request itself. Total, so ties cannot fall
    // back on the caller's order.
    std::vector<std::size_t> order(requests.size());
    std::iota(order.begin(), order.end(), std::size_t{0});
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
        const PlannedRepairWindow& x = planned[a];
        const PlannedRepairWindow& y = planned[b];
        const RepairRequest& p = requests[a];
        const RepairRequest& q = requests[b];
        return std::tie(x.already_held, x.key_hash, x.bytes, x.at_start, x.at_end, x.rel_block_begin,
                        x.rel_new_end, x.haplotype_hash, p.incumbent, p.block_begin, p.block_end,
                        p.replacement) <
               std::tie(y.already_held, y.key_hash, y.bytes, y.at_start, y.at_end, y.rel_block_begin,
                        y.rel_new_end, y.haplotype_hash, q.incumbent, q.block_begin, q.block_end,
                        q.replacement);
    });
    std::vector<std::size_t> slots(requests.size(), 0);
    if (works != nullptr) works->assign(requests.size(), BlockDeltaWork{});
    // I2-B: one round, one set of fine stores, released when the round ends however it ends
    struct FineRound {
        State* state = nullptr;
        ~FineRound() {
            if (state != nullptr) state->end_fine_round();
        }
    } fine_round;
    if (state.fine_enabled) {
        state.fine_round = true;
        ++state.fine_report.rounds;
        fine_round.state = &state;
    }
    for (const std::size_t i : order) {
        const RepairRequest& r = requests[i];
        slots[i] = add_repaired_haplotype(r.incumbent, r.block_begin, r.block_end, r.replacement,
                                          works != nullptr ? &(*works)[i] : nullptr);
    }
    return slots;
}

std::vector<std::uint64_t> PreparedScoringContext::origin_cache_key_hashes() const {
    const OriginCache& cache = state_->origin_cache;
    std::vector<std::uint64_t> out;
    for (std::size_t at = 0; at < cache.entries.size(); ++at) {
        if (cache.live[at] == 0) continue;
        const OriginCacheEntry& e = cache.entries[at];
        out.push_back(origin_key_hash(e.bytes, e.at_start, e.at_end, e.rel_block_begin, e.rel_new_end));
    }
    std::sort(out.begin(), out.end());
    return out;
}

std::uint64_t PreparedScoringContext::origin_cache_slots() const {
    return state_->origin_cache.entries.size();
}

void set_piece_screen_shadow(bool enabled) {
    if (enabled) {
        const std::lock_guard<std::mutex> lock(piece_screen_mutex());
        piece_screen_totals() = PieceScreenShadowReport{};
    }
    piece_screen_shadow_active().store(enabled, std::memory_order_relaxed);
}

PieceScreenShadowReport piece_screen_shadow_report() {
    const std::lock_guard<std::mutex> lock(piece_screen_mutex());
    return piece_screen_totals();
}

namespace {
// The window the prefilter would test, or 0 when gapped_emission refuses before it.
std::size_t prefilter_window(const std::string& haplotype, std::size_t start, const std::string& read,
                             std::size_t max_gap_bp) {
    if (start > haplotype.size() || haplotype.size() - start < read.size()) return 0;
    return std::min(haplotype.size() - start, read.size() + max_gap_bp);
}
}  // namespace

bool piece_screen_admits(const std::string& haplotype, std::size_t start, const std::string& read,
                         std::size_t band, std::size_t max_gap_bp, std::size_t pieces) {
    const std::size_t window = prefilter_window(haplotype, start, read, max_gap_bp);
    if (window == 0) return false;
    HaplotypeCodes codes;
    codes.text = &haplotype;
    return piece_screen(codes, start, window, piece_layout(read, pieces), band) == ScreenVerdict::Admit;
}

bool piece_screen_count_passes(const std::string& haplotype, std::size_t start, const std::string& read,
                               std::size_t band, std::size_t max_gap_bp, std::size_t pieces) {
    const std::size_t window = prefilter_window(haplotype, start, read, max_gap_bp);
    if (window == 0) return false;
    HaplotypeCodes codes;
    codes.text = &haplotype;
    return piece_screen(codes, start, window, piece_layout(read, pieces), band) != ScreenVerdict::CountReject;
}

// Written independently of piece_screen: plain substring comparison and a depth-first search over
// every choice of exact-or-skip per piece, pruned only by the same three conditions.
bool piece_screen_admits_reference(const std::string& haplotype, std::size_t start, const std::string& read,
                                   std::size_t band, std::size_t max_gap_bp, std::size_t pieces) {
    const std::size_t window = prefilter_window(haplotype, start, read, max_gap_bp);
    if (window == 0) return false;
    const std::size_t n = read.size();
    if (pieces == 0 || pieces > n || pieces <= band) return true;
    std::vector<std::size_t> offset(pieces), length(pieces);
    for (std::size_t j = 0; j < pieces; ++j) {
        offset[j] = j * n / pieces;
        length[j] = (j + 1) * n / pieces - offset[j];
        if (length[j] == 0 || length[j] > 31) return true;
    }
    const long k = static_cast<long>(band);
    std::vector<std::vector<long>> matches(pieces);
    for (std::size_t j = 0; j < pieces; ++j) {
        for (long delta = -k; delta <= k; ++delta) {
            const long at = static_cast<long>(start + offset[j]) + delta;
            if (at < static_cast<long>(start)) continue;
            if (static_cast<std::size_t>(at) + length[j] > start + window) continue;
            if (haplotype.compare(static_cast<std::size_t>(at), length[j], read, offset[j], length[j]) == 0) {
                matches[j].push_back(delta);
            }
        }
    }
    const std::size_t need = pieces - band;
    std::function<bool(std::size_t, bool, std::size_t, long, long, std::size_t)> search =
        [&](std::size_t j, bool have, std::size_t last, long last_delta, long spent, std::size_t count) {
            if (count >= need) return true;
            if (j == pieces || count + (pieces - j) < need) return false;
            for (const long delta : matches[j]) {
                long cost = spent;
                if (have) {
                    const long between = static_cast<long>(offset[j]) - static_cast<long>(offset[last] + length[last]);
                    if (delta - last_delta < -between) continue;
                    cost += std::labs(delta - last_delta);
                } else {
                    cost += std::labs(delta);
                }
                if (cost > k) continue;
                if (search(j + 1, true, j, delta, cost, count + 1)) return true;
            }
            return search(j + 1, have, last, last_delta, spent, count);
        };
    return search(0, false, 0, 0, 0, 0);
}

bool feasibility_prefilter_admits(const std::string& haplotype, std::size_t start, const std::string& read,
                                  std::size_t band, std::size_t max_gap_bp) {
    const std::size_t window = prefilter_window(haplotype, start, read, max_gap_bp);
    if (window == 0) return false;
    return feasible_within_band(haplotype, start, window, read, band, nullptr);
}

PrefilterDecisions feasibility_prefilter_decisions(const std::string& haplotype, std::size_t start,
                                                   const std::string& read, std::size_t band, std::size_t max_gap_bp) {
    PrefilterDecisions out;
    const std::size_t window = prefilter_window(haplotype, start, read, max_gap_bp);
    const MyersPattern pattern = prepare_myers_pattern(read);
    out.reference_unprepared = window != 0 && feasible_within_band_reference(haplotype, start, window, read, band, nullptr);
    out.reference_prepared = window != 0 && feasible_within_band_reference(haplotype, start, window, read, band, &pattern);
    out.production_unprepared = window != 0 && feasible_within_band(haplotype, start, window, read, band, nullptr);
    out.production_prepared = window != 0 && feasible_within_band(haplotype, start, window, read, band, &pattern);
    return out;
}

std::vector<std::pair<std::uint32_t, std::uint32_t>> sort_unique_proposals_scatter_for_test(
    std::vector<std::pair<std::uint32_t, std::uint32_t>> proposals, std::size_t form_count) {
    std::vector<std::uint32_t> cursor, bound, scratch;
    scratch.reserve(proposals.size());
    sort_unique_proposals_scatter(proposals, form_count, cursor, bound, scratch);
    return proposals;
}

std::vector<std::pair<std::uint32_t, std::uint32_t>> sort_unique_proposals_for_test(
    std::vector<std::pair<std::uint32_t, std::uint32_t>> proposals, std::size_t form_count, bool bucketed) {
    if (bucketed) {
        std::vector<std::uint32_t> cursor, bound;
        sort_unique_proposals_bucketed(proposals, form_count, cursor, bound);
    } else {
        std::sort(proposals.begin(), proposals.end());
        proposals.erase(std::unique(proposals.begin(), proposals.end()), proposals.end());
    }
    return proposals;
}

DemandWindow PreparedScoringContext::demand_window(std::uint64_t key_hash) const {
    DemandWindow out;
    const auto range = state_->demand_index.equal_range(key_hash);
    for (auto it = range.first; it != range.second; ++it) {
        const auto& row = state_->demand[it->second];
        out.found = true;
        out.bytes = row.first;
        out.at_start = row.second.at_start;
        out.at_end = row.second.at_end;
        out.rel_block_begin = row.second.rel_block_begin;
        out.rel_new_end = row.second.rel_new_end;
        break;
    }
    return out;
}

WindowEnumerationBench PreparedScoringContext::benchmark_window_enumeration(const std::string& bases,
                                                                            bool at_start, bool at_end,
                                                                            std::size_t workers,
                                                                            std::size_t repeats) {
    State& state = *state_;
    state.require_usable();
    if (workers == 0) workers = 1;
    WindowEnumerationBench out;
    out.workers = workers;
    out.fragments = state.fragments.size();
    out.window_bytes = bases.size();
    LocalWindow window;
    window.bases = bases;
    window.reversed.assign(bases.rbegin(), bases.rend());
    window.begin = 0;
    window.at_start = at_start;
    window.at_end = at_end;
    const std::size_t band = state.options.mismatch_band;
    const auto open_budget = [](PlacementBudget& budget) {
        budget.max_piece_hits = std::numeric_limits<std::uint64_t>::max();
        budget.max_physical_origins = std::numeric_limits<std::uint64_t>::max();
    };
    using Clock = std::chrono::steady_clock;
    const auto seconds_since = [](Clock::time_point from) {
        return std::chrono::duration<double>(Clock::now() - from).count();
    };
    // Reference: the production path, one thread, emission cache cold for this window.
    WindowScratch reference;
    FragmentScoringCounters reference_counters;
    EmissionCacheCounters reference_cache;
    {
        PlacementBudget budget;
        open_budget(budget);
        const auto mark = Clock::now();
        enumerate_window_origins(state.prepared, state.fragments, state.library, state.options, band, window,
                                 budget, reference_counters, reference_cache, reference);
        out.reference_cold_seconds = seconds_since(mark);
    }
    // The counter reference: the same path again, emission cache now warm.
    WindowScratch warm;
    FragmentScoringCounters warm_counters;
    EmissionCacheCounters warm_cache;
    {
        PlacementBudget budget;
        open_budget(budget);
        const auto mark = Clock::now();
        enumerate_window_origins(state.prepared, state.fragments, state.library, state.options, band, window,
                                 budget, warm_counters, warm_cache, warm);
        out.production_warm_seconds = seconds_since(mark);
    }
    // Contiguous fragment ranges, and each range's own piece indexes (built before any timing).
    const std::size_t n = state.fragments.size();
    struct Worker {
        std::size_t begin = 0, end = 0;
        std::vector<PieceIndex> forward, reverse;
        WindowScratch scratch;
        FragmentScoringCounters counters;
        EmissionCacheCounters cache;
        double busy = 0.0;
        EmissionTimingDiag timing;
    };
    std::vector<Worker> pool(workers);
    const auto restrict = [&](const std::vector<PieceIndex>& all, const std::deque<ReadForm>& forms,
                              std::size_t begin, std::size_t end) {
        std::vector<PieceIndex> part;
        for (const PieceIndex& index : all) {
            PieceIndex sub;
            sub.piece_bp = index.piece_bp;
            for (const PieceEntry& entry : index.entries) {
                const std::size_t fragment = forms[entry.form].fragment;
                if (fragment >= begin && fragment < end) sub.entries.push_back(entry);
            }
            part.push_back(std::move(sub));
        }
        return part;
    };
    for (std::size_t w = 0; w < workers; ++w) {
        pool[w].begin = n * w / workers;
        pool[w].end = n * (w + 1) / workers;
        pool[w].forward = restrict(state.prepared.forward_indexes, state.prepared.forward_forms, pool[w].begin, pool[w].end);
        pool[w].reverse = restrict(state.prepared.reverse_indexes, state.prepared.reverse_forms, pool[w].begin, pool[w].end);
    }
    emission_timing_active().store(true, std::memory_order_relaxed);
    for (std::size_t r = 0; r < std::max<std::size_t>(repeats, 1); ++r) {
        PlacementBudget budget;
        open_budget(budget);
        for (Worker& worker : pool) {
            worker.counters = FragmentScoringCounters{};
            worker.cache = EmissionCacheCounters{};
        }
        const auto mark = Clock::now();
        std::vector<std::thread> threads;
        for (Worker& worker : pool) {
            threads.emplace_back([&, target = &worker]() {
                const EmissionTimingDiag before = emission_timing();
                const auto started = Clock::now();
                enumerate_window_origins(state.prepared, state.fragments, state.library, state.options, band, window,
                                         budget, target->counters, target->cache, target->scratch,
                                         &target->forward, &target->reverse, target->begin, target->end);
                target->busy = seconds_since(started);
                target->timing.prefilter_ns = emission_timing().prefilter_ns - before.prefilter_ns;
                target->timing.emission_ns = emission_timing().emission_ns - before.emission_ns;
            });
        }
        for (std::thread& thread : threads) thread.join();
        out.wall_seconds.push_back(seconds_since(mark));
    }
    emission_timing_active().store(false, std::memory_order_relaxed);
    // Merge in ascending fragment order and compare.
    std::vector<HaploidOrigin> merged;
    std::vector<std::size_t> offsets(n + 1, 0);
    FragmentScoringCounters total;
    EmissionCacheCounters cache;
    for (std::size_t w = 0; w < workers; ++w) {
        const Worker& worker = pool[w];
        for (std::size_t f = worker.begin; f < worker.end; ++f) {
            offsets[f] = merged.size();
            merged.insert(merged.end(), worker.scratch.origins.begin() + static_cast<std::ptrdiff_t>(worker.scratch.offsets[f]),
                          worker.scratch.origins.begin() + static_cast<std::ptrdiff_t>(worker.scratch.offsets[f + 1]));
        }
        out.worker_busy_seconds.push_back(worker.busy);
        out.prefilter_seconds += static_cast<double>(worker.timing.prefilter_ns) * 1e-9;
        out.emission_seconds += static_cast<double>(worker.timing.emission_ns) * 1e-9;
        if (w == 0) {
            total.candidate_bases_scanned += worker.counters.candidate_bases_scanned;
        } else {
            out.duplicate_scan_bases += worker.counters.candidate_bases_scanned;
        }
        total.piece_hits += worker.counters.piece_hits;
        total.proposed_starts += worker.counters.proposed_starts;
        total.unique_starts += worker.counters.unique_starts;
        total.full_read_verifications += worker.counters.full_read_verifications;
        total.dp_evaluations += worker.counters.dp_evaluations;
        total.prefilter_rejections += worker.counters.prefilter_rejections;
        total.accepted_mate_placements += worker.counters.accepted_mate_placements;
        total.valid_fr_states += worker.counters.valid_fr_states;
        total.latent_alignment_pairs += worker.counters.latent_alignment_pairs;
        total.physical_origins += worker.counters.physical_origins;
        cache.executed_prefilters += worker.cache.executed_prefilters;
        cache.dp_cache_requests += worker.cache.dp_cache_requests;
        cache.dp_cache_hits += worker.cache.dp_cache_hits;
        cache.dp_cache_misses += worker.cache.dp_cache_misses;
        cache.executed_dps += worker.cache.executed_dps;
    }
    offsets[n] = merged.size();
    out.origins_identical = merged.size() == reference.origins.size();
    for (std::size_t i = 0; out.origins_identical && i < merged.size(); ++i) {
        const HaploidOrigin& a = merged[i];
        const HaploidOrigin& b = reference.origins[i];
        out.origins_identical = a.forward_start == b.forward_start && a.insert_bp == b.insert_bp &&
                                a.orientation == b.orientation &&
                                std::memcmp(&a.log_latent_mass, &b.log_latent_mass, sizeof(double)) == 0;
        if (!out.origins_identical) out.first_difference = "origin " + std::to_string(i);
    }
    if (out.origins_identical) {
        for (std::size_t f = 0; f <= n; ++f) {
            if (offsets[f] != reference.offsets[f]) {
                out.origins_identical = false;
                out.first_difference = "offset of fragment " + std::to_string(f);
                break;
            }
        }
    }
    const std::pair<const char*, std::pair<std::uint64_t, std::uint64_t>> checks[] = {
        {"candidate_bases_scanned", {total.candidate_bases_scanned, warm_counters.candidate_bases_scanned}},
        {"piece_hits", {total.piece_hits, warm_counters.piece_hits}},
        {"proposed_starts", {total.proposed_starts, warm_counters.proposed_starts}},
        {"unique_starts", {total.unique_starts, warm_counters.unique_starts}},
        {"full_read_verifications", {total.full_read_verifications, warm_counters.full_read_verifications}},
        {"dp_evaluations", {total.dp_evaluations, warm_counters.dp_evaluations}},
        {"prefilter_rejections", {total.prefilter_rejections, warm_counters.prefilter_rejections}},
        {"accepted_mate_placements", {total.accepted_mate_placements, warm_counters.accepted_mate_placements}},
        {"valid_fr_states", {total.valid_fr_states, warm_counters.valid_fr_states}},
        {"latent_alignment_pairs", {total.latent_alignment_pairs, warm_counters.latent_alignment_pairs}},
        {"physical_origins", {total.physical_origins, warm_counters.physical_origins}},
        {"executed_prefilters", {cache.executed_prefilters, warm_cache.executed_prefilters}},
        {"dp_cache_requests", {cache.dp_cache_requests, warm_cache.dp_cache_requests}},
        {"dp_cache_hits", {cache.dp_cache_hits, warm_cache.dp_cache_hits}},
        {"dp_cache_misses", {cache.dp_cache_misses, warm_cache.dp_cache_misses}},
        {"executed_dps", {cache.executed_dps, warm_cache.executed_dps}},
    };
    out.counters_identical = true;
    for (const auto& check : checks) {
        if (check.second.first != check.second.second) {
            out.counters_identical = false;
            if (out.first_difference.empty()) {
                out.first_difference = std::string(check.first) + " " + std::to_string(check.second.first) +
                                       " vs " + std::to_string(check.second.second);
            }
            break;
        }
    }
    out.origins = merged.size();
    out.piece_hits = total.piece_hits;
    out.latent_pairs = total.latent_alignment_pairs;
    out.executed_dps = cache.executed_dps;
    out.executed_prefilters = cache.executed_prefilters;
    out.verifications = total.full_read_verifications;
    return out;
}

std::size_t peak_active_workers() { return active_worker_peak().load(std::memory_order_relaxed); }

void reset_peak_active_workers() {
    active_worker_peak().store(active_worker_count().load(std::memory_order_relaxed), std::memory_order_relaxed);
}

ParallelEnumerationStats PreparedScoringContext::parallel_enumeration_stats() const {
    ParallelEnumerationStats out = state_->parallel_stats;
    out.workers = state_->enumeration_workers();
    return out;
}

std::uint64_t PreparedScoringContext::last_window_digest() const {
    const WindowScratch& w = state_->window_scratch;
    std::uint64_t h = 1469598103934665603ULL;
    const auto mix = [&h](const void* data, std::size_t bytes) {
        const unsigned char* p = static_cast<const unsigned char*>(data);
        for (std::size_t i = 0; i < bytes; ++i) h = (h ^ p[i]) * 1099511628211ULL;
    };
    for (const HaploidOrigin& o : w.origins) {
        mix(&o.forward_start, sizeof(o.forward_start));
        mix(&o.insert_bp, sizeof(o.insert_bp));
        mix(&o.orientation, sizeof(o.orientation));
        mix(&o.log_latent_mass, sizeof(o.log_latent_mass));
    }
    for (const std::size_t offset : w.offsets) mix(&offset, sizeof(offset));
    return h;
}

std::vector<std::size_t> PreparedScoringContext::census_slots() const {
    std::vector<std::size_t> out;
    for (const auto& item : state_->censuses) out.push_back(item.first);
    std::sort(out.begin(), out.end());
    return out;
}

std::size_t PreparedScoringContext::origin_cache_entries() const {
    return state_->origin_cache.index.size();
}
std::uint64_t PreparedScoringContext::origin_cache_bytes() const { return state_->origin_cache.bytes; }
std::uint64_t PreparedScoringContext::origin_cache_not_retained() const {
    return state_->origin_cache.not_retained;
}
std::uint64_t PreparedScoringContext::origin_records_consumed() const {
    return state_->origin_records_consumed;
}
std::uint64_t PreparedScoringContext::merge_working_bytes() const {
    std::uint64_t total = 0;
    for (const auto& item : state_->groups) total += sizeof(item) + item.second.unaffected.bytes();
    for (const auto& item : state_->full_aggregates) {
        total += sizeof(item) + static_cast<std::uint64_t>(item.second.capacity()) * sizeof(StateCount);
    }
    return total;
}

ContextMemoryAccount PreparedScoringContext::memory_account() const {
    const State& state = *state_;
    ContextMemoryAccount a;
    // A default string's capacity is its short-string buffer, inside the object itself: only the
    // excess is heap.
    const std::size_t inline_capacity = std::string().capacity();
    const auto heap = [inline_capacity](const std::string& text) -> std::uint64_t {
        return text.capacity() > inline_capacity ? text.capacity() : 0;
    };
    // Declared allowances for node-based containers, whose internals are not observable.
    constexpr std::uint64_t kTreeNode = 32;
    constexpr std::uint64_t kHashNode = 16;

    a.slots_total = state.sequences.size();
    for (std::size_t slot = 0; slot < state.sequences.size(); ++slot) {
        const bool live = !state.sequences[slot].empty() || !state.mass_column[slot].empty() ||
                          state.segments.count(slot) != 0;
        if (live) {
            ++a.slots_live;
        } else {
            ++a.slots_released;
        }
        a.sequence_bytes += state.sequences[slot].size() + state.reversed[slot].size();
        a.sequence_capacity += heap(state.sequences[slot]) + heap(state.reversed[slot]);
        a.mass_column_bytes += state.mass_column[slot].size() * sizeof(double);
        a.mass_column_capacity += state.mass_column[slot].capacity() * sizeof(double);
    }
    a.slot_table_bytes =
        static_cast<std::uint64_t>(state.sequences.size()) * 2 * sizeof(std::string) +
        static_cast<std::uint64_t>(state.mass_column.capacity()) * sizeof(std::vector<double>) +
        static_cast<std::uint64_t>(state.exposure.capacity()) * sizeof(double) +
        static_cast<std::uint64_t>(state.hash_by_slot.capacity()) * sizeof(std::uint64_t) +
        static_cast<std::uint64_t>(state.by_hash.size()) *
            (sizeof(std::pair<const std::uint64_t, std::size_t>) + kHashNode) +
        static_cast<std::uint64_t>(state.by_hash.bucket_count()) * sizeof(void*);
    for (const auto& segment : state.segments) {
        a.segment_bytes += sizeof(segment) + kTreeNode + heap(segment.second.replacement);
    }
    for (const auto& item : state.censuses) {
        ++a.census_count;
        for (const std::vector<HaploidOrigin>& list : item.second.origins_by_fragment) {
            a.census_bytes += list.size() * sizeof(HaploidOrigin);
        }
        a.census_capacity += item.second.payload_bytes() + sizeof(item) + kTreeNode;
    }
    a.origin_cache_entries = state.origin_cache.index.size();
    a.origin_cache_payload = state.origin_cache.bytes;
    for (const OriginCacheEntry& entry : state.origin_cache.entries) {
        a.origin_cache_capacity +=
            sizeof(OriginCacheEntry) + heap(entry.bytes) +
            static_cast<std::uint64_t>(entry.records.origins.capacity()) * sizeof(CachedOrigin) +
            static_cast<std::uint64_t>(entry.records.fragments.capacity() +
                                       entry.records.starts.capacity()) *
                sizeof(std::uint32_t);
    }
    a.origin_cache_capacity +=
        static_cast<std::uint64_t>(state.origin_cache.index.size()) *
            (sizeof(std::pair<const std::uint64_t, std::size_t>) + kHashNode) +
        static_cast<std::uint64_t>(state.origin_cache.index.bucket_count()) * sizeof(void*);
    for (const auto& item : state.groups) {
        a.merge_working_bytes +=
            item.second.unaffected.states.size() * sizeof(StateCount) +
            (item.second.unaffected.fragments.size() + item.second.unaffected.starts.size()) *
                sizeof(std::uint32_t);
        a.merge_working_capacity += item.second.unaffected.bytes() + sizeof(item) + kTreeNode;
    }
    for (const auto& item : state.full_aggregates) {
        a.merge_working_bytes += item.second.size() * sizeof(StateCount);
        a.merge_working_capacity +=
            item.second.capacity() * sizeof(StateCount) + sizeof(item) + kTreeNode;
    }
    {
        const WindowScratch& w = state.window_scratch;
        std::uint64_t cap = static_cast<std::uint64_t>(w.forward_rows.capacity() + w.reverse_rows.capacity()) *
                            sizeof(std::vector<FormPlacement>);
        for (const auto& row : w.forward_rows) cap += row.capacity() * sizeof(FormPlacement);
        for (const auto& row : w.reverse_rows) cap += row.capacity() * sizeof(FormPlacement);
        cap += static_cast<std::uint64_t>(w.touched_forward.capacity() + w.touched_reverse.capacity()) *
                   sizeof(std::uint32_t) +
               static_cast<std::uint64_t>(w.proposed.capacity()) *
                   sizeof(std::pair<std::uint32_t, std::uint32_t>) +
               static_cast<std::uint64_t>(w.group_begin.capacity() + w.group_end.capacity() +
                                          w.offsets.capacity()) *
                   sizeof(std::size_t) +
               static_cast<std::uint64_t>(w.accumulators.capacity()) * sizeof(OriginAccumulator) +
               w.accepted.capacity() +
               static_cast<std::uint64_t>(w.origins.capacity()) * sizeof(HaploidOrigin);
        // I2-B: live fine stores (held + scratch, capacity); zero outside a round, since a round releases them
        cap += w.fine.account_bytes();
        for (const auto& p : state.partitions) cap += p.scratch.fine.account_bytes();
        a.window_scratch_capacity = cap;
    }
    {
        const auto records = [](const OriginRecordSet& r) {
            return static_cast<std::uint64_t>(r.origins.capacity()) * sizeof(CachedOrigin) +
                   static_cast<std::uint64_t>(r.fragments.capacity() + r.starts.capacity()) *
                       sizeof(std::uint32_t);
        };
        a.origin_scratch_capacity =
            records(state.miss_records) + records(state.check_records) +
            state.candidate_states.bytes() +
            static_cast<std::uint64_t>(state.key_scratch.capacity()) * sizeof(StateKey) +
            static_cast<std::uint64_t>(state.merged.capacity()) * sizeof(StateCount);
    }
    emission_cache_storage_bytes(&a.emission_cache_table, &a.emission_cache_arena);
    {
        std::uint64_t prepared = 0;
        for (const ReadForm& form : state.prepared.forward_forms) prepared += sizeof(ReadForm) + heap(form.bases);
        for (const ReadForm& form : state.prepared.reverse_forms) prepared += sizeof(ReadForm) + heap(form.bases);
        for (const PieceIndex& index : state.prepared.forward_indexes) {
            prepared += index.entries.capacity() * sizeof(PieceEntry);
        }
        for (const PieceIndex& index : state.prepared.reverse_indexes) {
            prepared += index.entries.capacity() * sizeof(PieceEntry);
        }
        for (const ReadPair& fragment : state.fragments) {
            prepared += sizeof(ReadPair) + heap(fragment.first.bases) + heap(fragment.second.bases);
        }
        a.prepared_capacity = prepared;
    }
    a.last_scoring_result_bytes = state.last_result_bytes;
    for (const auto& row : state.demand) {
        a.origin_demand_bytes += sizeof(row) + heap(row.first);
    }
    a.origin_demand_bytes += static_cast<std::uint64_t>(state.demand_index.size()) *
                                 (sizeof(std::pair<const std::uint64_t, std::size_t>) + kHashNode) +
                             static_cast<std::uint64_t>(state.demand_index.bucket_count()) * sizeof(void*);
    a.logical_total = a.slot_table_bytes + a.sequence_bytes + a.segment_bytes + a.mass_column_bytes +
                      a.origin_demand_bytes +
                      a.census_bytes + a.origin_cache_payload + a.merge_working_bytes +
                      a.emission_cache_table + a.emission_cache_arena + a.prepared_capacity +
                      a.last_scoring_result_bytes;
    a.capacity_total = a.slot_table_bytes + a.sequence_capacity + a.segment_bytes + a.origin_demand_bytes +
                       a.mass_column_capacity + a.census_capacity + a.origin_cache_capacity +
                       a.merge_working_capacity + a.window_scratch_capacity +
                       a.origin_scratch_capacity + a.emission_cache_table + a.emission_cache_arena +
                       a.prepared_capacity + a.last_scoring_result_bytes;
    return a;
}

void set_origin_cache_self_check(bool enabled) {
    scoring_detail::origin_self_check().store(enabled, std::memory_order_relaxed);
}

void set_union_prefilter_cache_hash_collapse(bool collapse) {
    union_prefilter_hash_collapse().store(collapse, std::memory_order_relaxed);
}

void set_union_placement_cache_hash_collapse(bool collapse) {
    scoring_detail::union_hash_collapse().store(collapse, std::memory_order_relaxed);
}

void set_origin_cache_hash_collapse(bool collapse) {
    scoring_detail::origin_hash_collapse().store(collapse, std::memory_order_relaxed);
}

void set_window_enumeration_self_check(bool enabled) {
    scoring_detail::window_self_check().store(enabled, std::memory_order_relaxed);
}

namespace {
SegmentedSequence segment_of(const std::string& base, std::size_t block_begin, std::size_t block_end,
                             const std::string& replacement) {
    if (block_end < block_begin || block_end > base.size()) {
        throw std::invalid_argument("the replaced block is not an interval of the base");
    }
    SegmentedSequence segment;
    segment.base = &base;
    segment.block_begin = block_begin;
    segment.block_end = block_end;
    segment.replacement = replacement;
    return segment;
}
}  // namespace

std::size_t segmented_length(const std::string& base, std::size_t block_begin, std::size_t block_end,
                             const std::string& replacement) {
    return segment_of(base, block_begin, block_end, replacement).length();
}

std::uint64_t segmented_hash(const std::string& base, std::size_t block_begin, std::size_t block_end,
                             const std::string& replacement) {
    return segment_of(base, block_begin, block_end, replacement).hash();
}

void segmented_window(const std::string& base, std::size_t block_begin, std::size_t block_end,
                      const std::string& replacement, std::size_t begin, std::size_t end,
                      std::string& forward_out, std::string& reversed_out) {
    const SegmentedSequence segment = segment_of(base, block_begin, block_end, replacement);
    if (end > segment.length() || begin > end) {
        throw std::invalid_argument("the requested window is not a range of the repaired haplotype");
    }
    forward_out.clear();
    segment.append_range(begin, end, forward_out);
    reversed_out.clear();
    segment.append_range_reversed(begin, end, reversed_out);
}

std::uint64_t HaploidOriginCensus::payload_bytes() const {
    std::uint64_t total = sizeof(HaploidOriginCensus);
    for (const std::vector<HaploidOrigin>& list : origins_by_fragment) {
        total += sizeof(list) + static_cast<std::uint64_t>(list.capacity()) * sizeof(HaploidOrigin);
    }
    return total;
}

void PreparedScoringContext::begin_transient_round() {
    State& state = *state_;
    state.require_usable();
    // Nothing is frozen. Owned cache keys mean a transient haplotype may be destroyed while the
    // entries it produced keep serving, so all this records is the slot from which release is allowed.
    state.transient_watermark = state.sequences.size();
}

const HaploidOriginCensus* PreparedScoringContext::census_view(std::size_t haplotype) const {
    const auto found = state_->censuses.find(haplotype);
    return found == state_->censuses.end() ? nullptr : &found->second;
}

bool PreparedScoringContext::holds_census(std::size_t haplotype) const {
    return state_->censuses.count(haplotype) != 0;
}

void PreparedScoringContext::retain_census(std::size_t haplotype) {
    State& state = *state_;
    state.require_usable();
    if (haplotype >= state.sequences.size()) {
        throw std::invalid_argument("census requested for an unknown haplotype");
    }
    if (state.censuses.count(haplotype) != 0) return;
    // Placed in full, because nothing is held to delta from. One of these per round is the price of
    // being able to delta the next round's several hundred neighbours.
    HaploidOriginCensus census;
    PlacementInputs inputs;
    inputs.prepared = &state.prepared;
    inputs.fragments = &state.fragments;
    inputs.library = &state.library;
    inputs.library_insert = &state.library_insert;
    inputs.options = &state.options;
    inputs.sequences = &state.sequences;
    inputs.reversed = &state.reversed;
    inputs.band = state.options.mismatch_band;
    CensusCapture capture;
    capture.wanted.emplace(haplotype, &census);
    OriginRecording none;
    FragmentScoringCounters counters;
    ScorerMemoryAccounts memory;
    EmissionCacheCounters cache_work;
    std::vector<std::vector<double>> scratch = state.mass_column;
    ++state.full_placements;
    try {
        place_haplotype_range(inputs, haplotype, haplotype + 1, 1, state.budget, scratch, counters,
                              memory, cache_work, none, &capture);
    } catch (const ResourceExhausted&) {
        state.refused = true;
        throw;
    }
    state.mass_column[haplotype] = scratch[haplotype];
    state.censuses.emplace(haplotype, std::move(census));
    // This placement is search work like any other. Leaving it out of the cumulative counters while
    // the shared budget atomic already counted it is what made a cap measured from those counters
    // refuse a run that did exactly the same work.
    state.cumulative.candidate_bases_scanned += counters.candidate_bases_scanned;
    state.cumulative.piece_hits += counters.piece_hits;
    state.cumulative.proposed_starts += counters.proposed_starts;
    state.cumulative.unique_starts += counters.unique_starts;
    state.cumulative.full_read_verifications += counters.full_read_verifications;
    state.cumulative.dp_evaluations += counters.dp_evaluations;
    state.cumulative.prefilter_rejections += counters.prefilter_rejections;
    state.cumulative.accepted_mate_placements += counters.accepted_mate_placements;
    state.cumulative.valid_fr_states += counters.valid_fr_states;
    state.cumulative.latent_alignment_pairs += counters.latent_alignment_pairs;
    state.cumulative.physical_origins += counters.physical_origins;
    state.cumulative.statistical_state_classes += counters.statistical_state_classes;
    state.cumulative_cache.executed_prefilters += cache_work.executed_prefilters;
    state.cumulative_cache.dp_cache_requests += cache_work.dp_cache_requests;
    state.cumulative_cache.dp_cache_hits += cache_work.dp_cache_hits;
    state.cumulative_cache.dp_cache_misses += cache_work.dp_cache_misses;
    state.cumulative_cache.executed_dps += cache_work.executed_dps;
    try {
        enforce_acceptance_caps(state.cumulative, state.options, []() {});
    } catch (const ResourceExhausted&) {
        state.refused = true;
        throw;
    }
}

std::size_t PreparedScoringContext::add_repaired_haplotype(
    std::size_t incumbent,
    std::size_t block_begin,
    std::size_t block_end,
    const std::string& replacement,
    BlockDeltaWork* work,
    bool retain_census_for_it) {
    State& state = *state_;
    state.require_usable();
    // Refusing is deliberate. Placing in full when no census is held would report a delta that never
    // happened, and the whole point of this path is that its reuse is measurable.
    const auto held = state.censuses.find(incumbent);
    if (held == state.censuses.end()) {
        throw std::invalid_argument("no origin census is held for the incumbent haplotype");
    }
    const HaploidOriginCensus& census = held->second;
    // The incumbent must have real bytes: a delta reads its prefix and suffix directly.
    if (state.sequences[incumbent].empty()) {
        throw std::invalid_argument("the incumbent haplotype has not been materialised");
    }
    const std::string& old_sequence = state.sequences[incumbent];
    if (block_end < block_begin || block_end > old_sequence.size()) {
        throw std::invalid_argument("the replaced block is not an interval of the incumbent");
    }

    const auto mark = std::chrono::steady_clock::now();
    BlockDeltaWork local;
    const auto hash_mark = std::chrono::steady_clock::now();
    // The repair, as three immutable segments. No haplotype is copied and none is reversed.
    SegmentedSequence repaired;
    repaired.base = &old_sequence;
    repaired.block_begin = block_begin;
    repaired.block_end = block_end;
    repaired.replacement = replacement;
    const std::uint64_t repaired_hash = repaired.hash();
    local.sequence_key_bytes = replacement.size() + sizeof(SegmentedSequence);
    const std::size_t old_length = old_sequence.size();
    const std::size_t new_length = repaired.length();
    const std::size_t new_end = block_begin + replacement.size();
    const std::ptrdiff_t shift =
        static_cast<std::ptrdiff_t>(replacement.size()) -
        static_cast<std::ptrdiff_t>(block_end - block_begin);

    // Already held? Hand back the same slot. Found by streamed hash, then confirmed byte for byte --
    // two different allele paths may spell the identical haplotype and must share one column. The one
    // exception is a caller that now wants a census for it (the round's winner, scored earlier
    // without one), which must fill the EXISTING slot rather than mint a second for identical bytes.
    std::size_t held_slot = state.sequences.size();
    {
        const auto range = state.by_hash.equal_range(repaired_hash);
        for (auto it = range.first; it != range.second; ++it) {
            if (state.slot_equals(it->second, repaired)) {
                held_slot = it->second;
                break;
            }
        }
    }
    local.seconds_hashing_and_windows =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - hash_mark).count();
    const bool already_held = held_slot != state.sequences.size();
    if (already_held && (!retain_census_for_it || state.censuses.count(held_slot) != 0)) {
        if (work != nullptr) {
            local.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - mark).count();
            *work = local;
        }
        return held_slot;
    }

    const State::RepairWindowBounds bounds =
        state.repair_window_bounds(block_begin, block_end, old_length, new_length, new_end, shift);
    const std::size_t margin = bounds.margin;
    const std::size_t shortest = bounds.shortest;

    // Affected on the OLD haplotype. Half-open: an origin ending exactly at block_begin, or starting
    // exactly at block_end, does not intersect and is reused.
    const auto affected = [&](const HaploidOrigin& origin) {
        const std::size_t begin = origin.forward_start;
        const std::size_t end = begin + origin.insert_bp;
        if (begin < block_end && block_begin < end) return true;
        if (shift == 0) return false;
        if (end <= block_begin) {
            // Before the block: its forward mate's available length is L - begin, which the length
            // change moves. Unsaturated under EITHER length means the key is not safe to reuse.
            return shortest < begin + margin;
        }
        // After the block: its reverse mate's available length is the right endpoint, which shifts.
        const std::size_t moved_end = static_cast<std::size_t>(static_cast<std::ptrdiff_t>(end) + shift);
        return std::min(end, moved_end) < margin;
    };

    LocalWindow window;
    window.begin = bounds.begin;
    const std::size_t window_end = bounds.end;
    window.at_start = window.begin == 0;
    window.at_end = window_end == new_length;
    {
        const auto window_mark = std::chrono::steady_clock::now();
        window.bases.clear();
        repaired.append_range(window.begin, window_end, window.bases);
        window.reversed.clear();
        repaired.append_range_reversed(window.begin, window_end, window.reversed);
        local.seconds_hashing_and_windows +=
            std::chrono::duration<double>(std::chrono::steady_clock::now() - window_mark).count();
    }
    local.local_window_bases = window.bases.size();
    local.local_window_bytes_materialised =
        static_cast<std::uint64_t>(window.bases.size() + window.reversed.size());
    local.local_window_payload_bytes =
        static_cast<std::uint64_t>(window.bases.capacity() + window.reversed.capacity());

    // Affected on the NEW haplotype, mirroring the test above so the two halves partition exactly.
    const auto affected_new = [&](const HaploidOrigin& origin) {
        const std::size_t begin = origin.forward_start;
        const std::size_t end = begin + origin.insert_bp;
        if (begin < new_end && block_begin < end) return true;
        if (shift == 0) return false;
        if (end <= block_begin) return shortest < begin + margin;
        const std::size_t before = static_cast<std::size_t>(static_cast<std::ptrdiff_t>(end) - shift);
        return std::min(before, end) < margin;
    };

    // ==== stage 6: window-origin cache + preaggregated sorted merge. Off unless enabled; when off,
    // ==== everything below this block is the stage-5 path, unchanged.
    if (state.origin_path) {
        using Clock = std::chrono::steady_clock;
        const auto seconds_since = [](Clock::time_point from) {
            return std::chrono::duration<double>(Clock::now() - from).count();
        };
        const std::size_t fragment_count = state.fragments.size();
        const std::size_t base_offset = window.begin;
        const std::uint32_t rel_block_begin = static_cast<std::uint32_t>(block_begin - base_offset);
        const std::uint32_t rel_new_end = static_cast<std::uint32_t>(new_end - base_offset);
        // At an end the margin conditions can fire, so nothing may be pruned there.
        const bool keep_all = window.at_start || window.at_end;
        const bool checking = origin_self_check().load(std::memory_order_relaxed);

        FragmentScoringCounters counters;
        EmissionCacheCounters cache_work;
        const std::uint64_t rows_before = state.window_scratch.row_constructions;
        const std::uint64_t accumulators_before = state.window_scratch.accumulator_growths;
        const std::uint64_t origins_before = state.window_scratch.origin_growths;

        // ---- 1. the window's origins: loaded from the cache, or enumerated and retained if it fits
        const auto enumeration_mark = Clock::now();
        const std::uint64_t key_hash = origin_key_hash(window.bases, window.at_start, window.at_end,
                                                       rel_block_begin, rel_new_end);
        const OriginCacheEntry* entry =
            find_origin_entry(state.origin_cache, key_hash, window.bases, window.at_start,
                              window.at_end, rel_block_begin, rel_new_end);
        OriginDemandRow* demand = nullptr;
        {
            const auto range = state.demand_index.equal_range(key_hash);
            for (auto it = range.first; it != range.second && demand == nullptr; ++it) {
                auto& row = state.demand[it->second];
                if (row.second.at_start == window.at_start && row.second.at_end == window.at_end &&
                    row.second.rel_block_begin == rel_block_begin &&
                    row.second.rel_new_end == rel_new_end && row.first == window.bases) {
                    demand = &row.second;
                }
            }
            if (demand == nullptr) {
                OriginDemandRow fresh;
                fresh.order = state.demand.size();
                fresh.key_hash = key_hash;
                fresh.window_bytes = window.bases.size();
                fresh.at_start = window.at_start;
                fresh.at_end = window.at_end;
                fresh.rel_block_begin = rel_block_begin;
                fresh.rel_new_end = rel_new_end;
                state.demand.emplace_back(window.bases, fresh);
                state.demand_index.emplace(key_hash, state.demand.size() - 1);
                demand = &state.demand.back().second;   // deque: stable address
            }
            ++demand->requests;
        }
        const OriginRecordSet* source = nullptr;
        bool scratch_is_fresh = false;
        if (entry != nullptr) {
            ++local.origin_cache_hits;
            ++demand->hits;
            const std::uint64_t records = entry->records.origins.size();
            // A separate cap on what REUSE consumes. Executed-work caps are not charged for a hit --
            // that work was not performed -- so reuse needs a bound of its own.
            if (state.origin_record_cap != 0 &&
                state.origin_records_consumed + records > state.origin_record_cap) {
                state.refused = true;
                throw ResourceExhausted(
                    "the origin cache reached its declared consumption cap; no partial ranking is "
                    "returned");
            }
            state.origin_records_consumed += records;
            local.cached_origin_records_consumed = records;
            source = &entry->records;
        } else if (state.origin_planned && state.holder_valid && state.holder_at_start == window.at_start &&
                   state.holder_at_end == window.at_end &&
                   state.holder_rel_block_begin == rel_block_begin &&
                   state.holder_rel_new_end == rel_new_end && state.holder_bytes == window.bases) {
            // The same window the previous enumeration built: its records are still in miss_records.
            // Like a hit, it executes nothing and is bounded by the consumption cap instead.
            ++local.origin_group_reuses;
            ++demand->group_reuses;
            const std::uint64_t records = state.miss_records.origins.size();
            if (state.origin_record_cap != 0 &&
                state.origin_records_consumed + records > state.origin_record_cap) {
                state.refused = true;
                throw ResourceExhausted(
                    "the origin cache reached its declared consumption cap; no partial ranking is "
                    "returned");
            }
            state.origin_records_consumed += records;
            local.cached_origin_records_consumed = records;
            source = &state.miss_records;
        } else {
            ++local.origin_cache_misses;
            const bool size_known_before = demand->misses != 0;
            const auto miss_mark = Clock::now();
            try {
                state.enumerate_window(window, state.budget, counters, cache_work, state.window_scratch);
            } catch (const ResourceExhausted&) {
                state.refused = true;
                throw;
            }
            scratch_is_fresh = true;
            local.physical_origins_constructed = counters.physical_origins;
            build_origin_records(state.window_scratch, fragment_count, base_offset, rel_block_begin,
                                 rel_new_end, keep_all, state.miss_records);
            source = &state.miss_records;
            state.holder_valid = true;
            state.holder_bytes = window.bases;
            state.holder_at_start = window.at_start;
            state.holder_at_end = window.at_end;
            state.holder_rel_block_begin = rel_block_begin;
            state.holder_rel_new_end = rel_new_end;
            const std::uint64_t cost = origin_entry_bytes(window.bases, state.miss_records);
            {
                const double spent = seconds_since(miss_mark);
                if (demand->misses != 0) {
                    demand->seconds_reenumeration += spent;
                    local.seconds_reenumeration += spent;
                }
                ++demand->misses;
                demand->seconds_enumeration += spent;
                demand->origins_constructed += counters.physical_origins;
                demand->entry_bytes = cost;
                demand->records = state.miss_records.origins.size();
            }
            // First fit keeps any miss that fits. Planned retention considers a window the plan selected,
            // or one whose exact size was unknown until this miss; a known window the plan left out
            // stays out. Both planned cases may displace strictly lower-ranked entries.
            bool selected = true;
            bool replace = false;
            std::uint64_t multiplicity = 1;
            OriginIdentity identity;
            if (state.origin_planned) {
                selected = (demand->order < state.plan_selected.size() &&
                            state.plan_selected[demand->order] != 0) ||
                           !size_known_before;
                replace = selected;
                identity = OriginIdentity{key_hash, window.bases, window.at_start, window.at_end,
                                          rel_block_begin, rel_new_end};
                const auto planned = state.plan_multiplicity.find(identity);
                if (planned != state.plan_multiplicity.end()) multiplicity = planned->second;
            }
            bool fits = state.origin_cache.bytes + cost <= state.origin_cache_budget;
            if (state.origin_cache_budget != 0 && selected && !fits && replace) {
                OriginCache& cache = state.origin_cache;
                std::vector<std::pair<std::size_t, OriginIdentity>> below;
                for (std::size_t at = 0; at < cache.entries.size(); ++at) {
                    if (cache.live[at] == 0) continue;
                    OriginIdentity other = entry_identity(cache.entries[at]);
                    if (ranks_before(multiplicity, cost, identity, cache.multiplicity[at],
                                     cache.entries[at].payload, other)) {
                        below.emplace_back(at, std::move(other));
                    }
                }
                // Worst first.
                std::sort(below.begin(), below.end(), [&](const auto& a, const auto& b) {
                    return ranks_before(cache.multiplicity[b.first], cache.entries[b.first].payload, b.second,
                                        cache.multiplicity[a.first], cache.entries[a.first].payload, a.second);
                });
                std::uint64_t freeable = 0;
                std::size_t take = 0;
                while (take < below.size() && cache.bytes - freeable + cost > state.origin_cache_budget) {
                    freeable += cache.entries[below[take].first].payload;
                    ++take;
                }
                if (cache.bytes - freeable + cost <= state.origin_cache_budget) {
                    for (std::size_t k = 0; k < take; ++k) {
                        ++state.demand[cache.demand_of[below[k].first]].second.evictions;
                        local.origin_cache_evicted_bytes += evict_origin_slot(cache, below[k].first);
                        ++local.origin_cache_evicted;
                    }
                    fits = true;
                }
            }
            if (state.origin_cache_budget != 0 && selected && fits) {
                OriginCacheEntry stored;
                stored.bytes = window.bases;
                stored.at_start = window.at_start;
                stored.at_end = window.at_end;
                stored.rel_block_begin = rel_block_begin;
                stored.rel_new_end = rel_new_end;
                stored.records = state.miss_records;
                stored.payload = cost;
                OriginCache& cache = state.origin_cache;
                std::size_t at = 0;
                if (!cache.free_slots.empty()) {
                    at = cache.free_slots.back();
                    cache.free_slots.pop_back();
                    cache.entries[at] = std::move(stored);
                    cache.live[at] = 1;
                    cache.demand_of[at] = demand->order;
                    cache.multiplicity[at] = multiplicity;
                } else {
                    cache.entries.push_back(std::move(stored));
                    cache.live.push_back(1);
                    cache.demand_of.push_back(demand->order);
                    cache.multiplicity.push_back(multiplicity);
                    at = cache.entries.size() - 1;
                }
                cache.index.emplace(key_hash, at);
                cache.bytes += cost;
                ++local.origin_cache_stored;
                local.origin_cache_admitted_bytes += cost;
                demand->admitted = true;
                source = &cache.entries[at].records;   // deque: stable address
                if (cache.bytes > state.origin_cache_budget) {
                    throw std::logic_error("origin cache exceeded its declared budget");
                }
            } else if (state.origin_cache_budget != 0 && !selected) {
                ++local.origin_cache_not_selected;
            } else if (state.origin_cache_budget != 0) {
                // Exhaustion stops RETENTION. The records just built are still used, exactly.
                ++state.origin_cache.not_retained;
                ++local.origin_cache_declined;
                local.origin_cache_declined_bytes += cost;
                ++demand->declined;
            }
        }
        // Certification of a hit: enumerate fresh against a throwaway budget and require the
        // retained records to be exactly what a fresh enumeration would retain.
        if (checking && !scratch_is_fresh) {
            PlacementBudget spare;
            spare.max_piece_hits = std::numeric_limits<std::uint64_t>::max();
            spare.max_physical_origins = std::numeric_limits<std::uint64_t>::max();
            spare.publish_batch = std::numeric_limits<std::uint64_t>::max();
            FragmentScoringCounters spare_counters;
            EmissionCacheCounters spare_cache;
            enumerate_window_origins(state.prepared, state.fragments, state.library, state.options,
                                     state.options.mismatch_band, window, spare, spare_counters,
                                     spare_cache, state.window_scratch);
            scratch_is_fresh = true;
            build_origin_records(state.window_scratch, fragment_count, base_offset, rel_block_begin,
                                 rel_new_end, keep_all, state.check_records);
            if (!(state.check_records == *source)) {
                throw std::runtime_error(
                    "origin-cache self-check: retained origins differ from a fresh enumeration");
            }
        }
        local.seconds_enumeration = seconds_since(enumeration_mark);
        local.seconds_placement_and_origins = local.seconds_enumeration;
        local.placement_row_constructions = state.window_scratch.row_constructions - rows_before;
        local.accumulator_growths = state.window_scratch.accumulator_growths - accumulators_before;
        local.origin_output_growths = state.window_scratch.origin_growths - origins_before;

        // ---- 2. the EXPLICIT affected filter, on absolute coordinates, aggregated per fragment.
        // Records are StateKey-sorted within a fragment, so the kept subsequence is too, and each
        // multiplicity is a run length.
        const auto filter_mark = Clock::now();
        StateAggregateSet& candidate = state.candidate_states;
        candidate.clear();
        for (std::size_t i = 0; i < source->fragments.size(); ++i) {
            const std::size_t first = candidate.states.size();
            for (std::uint32_t at = source->starts[i]; at < source->starts[i + 1]; ++at) {
                const CachedOrigin& record = source->origins[at];
                HaploidOrigin origin;
                origin.forward_start = static_cast<std::uint32_t>(base_offset + record.rel_start);
                origin.insert_bp = record.insert();
                origin.orientation = record.arm();
                origin.log_latent_mass = record.log_latent_mass;
                if (!affected_new(origin)) continue;
                ++local.origins_recomputed;
                StateKey key;
                key.insert_bp = origin.insert_bp;
                key.orientation = origin.orientation;
                key.log_latent_mass = origin.log_latent_mass;
                if (candidate.states.size() > first && same_state(candidate.states.back().key, key)) {
                    ++candidate.states.back().count;
                } else {
                    candidate.states.push_back({key, 1});
                }
            }
            if (candidate.states.size() > first) {
                candidate.fragments.push_back(source->fragments[i]);
                candidate.starts.push_back(static_cast<std::uint32_t>(first));
            }
        }
        candidate.starts.push_back(static_cast<std::uint32_t>(candidate.states.size()));
        local.seconds_affected_filter = seconds_since(filter_mark);
        if (!demand->aggregate_measured) {
            // Measurement only: nothing below reads these, and nothing is retained differently.
            demand->aggregate_measured = true;
            demand->aggregate_states = candidate.states.size();
            demand->aggregate_fragments = candidate.fragments.size();
            std::uint64_t multiplicity = 0;
            for (const StateCount& item : candidate.states) multiplicity += item.count;
            demand->aggregate_multiplicity = multiplicity;
            demand->retained_not_affected = source->origins.size() - multiplicity;
            const std::uint64_t shared = sizeof(OriginCacheEntry) + kOriginCacheNodeBytes +
                                         window.bases.size() +
                                         static_cast<std::uint64_t>(2 * candidate.fragments.size() + 1) *
                                             sizeof(std::uint32_t);
            demand->aggregate_entry_bytes =
                shared + static_cast<std::uint64_t>(candidate.states.size()) * sizeof(StateCount);
            demand->aggregate_entry_bytes_packed =
                shared + static_cast<std::uint64_t>(candidate.states.size()) * 16u;
        }

        // ---- 3. the unaffected incumbent origins, preaggregated once per exact (block, shift)
        const auto merge_mark = Clock::now();
        const GroupKey group_key{incumbent, block_begin, block_end, shift};
        auto group_it = state.groups.find(group_key);
        if (group_it == state.groups.end()) {
            GroupEntry built;
            for (std::size_t f = 0; f < fragment_count; ++f) {
                const std::vector<HaploidOrigin>& list = census.origins_by_fragment[f];
                std::size_t hit = 0;
                for (const HaploidOrigin& origin : list) {
                    if (affected(origin)) ++hit;
                }
                built.unaffected_origins += list.size() - hit;
                if (hit == 0) continue;   // this fragment's unaffected set is its whole census
                state.key_scratch.clear();
                for (const HaploidOrigin& origin : list) {
                    if (affected(origin)) continue;
                    StateKey key;
                    key.insert_bp = origin.insert_bp;
                    key.orientation = origin.orientation;
                    key.log_latent_mass = origin.log_latent_mass;
                    state.key_scratch.push_back(key);
                }
                const std::size_t first = built.unaffected.states.size();
                built.unaffected.fragments.push_back(static_cast<std::uint32_t>(f));
                built.unaffected.starts.push_back(static_cast<std::uint32_t>(first));
                sort_and_count(state.key_scratch, built.unaffected.states, first);
            }
            built.unaffected.starts.push_back(
                static_cast<std::uint32_t>(built.unaffected.states.size()));
            group_it = state.groups.emplace(group_key, std::move(built)).first;
        }
        const GroupEntry& group = group_it->second;
        local.origins_reused = group.unaffected_origins;

        // ---- 4. the slot, exactly as stage 5 makes it
        std::size_t slot = 0;
        if (already_held) {
            slot = held_slot;
        } else {
            slot = state.allocate_repair_slot(repaired_hash);
            state.segments.emplace(slot, repaired);
            state.exposure[slot] = exposure_for_length(new_length, fragment_count,
                                                       state.prepared.fragments_by_floor,
                                                       state.library_insert);
        }
        if (retain_census_for_it) {
            const bool was_segmented = state.segments.count(slot) != 0;
            state.materialise(slot);
            local.full_haplotype_bytes_materialised =
                was_segmented ? static_cast<std::uint64_t>(state.sequences[slot].size() +
                                                           state.reversed[slot].size())
                              : 0;
        }

        // ---- 5. the column. A fragment with no affected origin on EITHER side has the incumbent's
        // StateKey multiset, so its reduction is the incumbent's value bit for bit: copied. Every
        // other fragment is merged linearly and reduced in ascending StateKey order.
        std::vector<double>& column = state.mass_column[slot];
        if (slot != incumbent) column = state.mass_column[incumbent];
        std::map<std::uint32_t, std::vector<StateCount>> merged_for_check;
        std::uint64_t reduced = 0;
        double reduce_seconds = 0.0;
        const std::vector<std::uint32_t>& group_fragments = group.unaffected.fragments;
        const std::vector<std::uint32_t>& candidate_fragments = candidate.fragments;
        std::size_t gi = 0;
        std::size_t ci = 0;
        while (gi < group_fragments.size() || ci < candidate_fragments.size()) {
            bool in_group = false;
            bool in_candidate = false;
            std::uint32_t f = 0;
            if (ci == candidate_fragments.size() ||
                (gi < group_fragments.size() && group_fragments[gi] < candidate_fragments[ci])) {
                f = group_fragments[gi];
                in_group = true;
            } else if (gi == group_fragments.size() || candidate_fragments[ci] < group_fragments[gi]) {
                f = candidate_fragments[ci];
                in_candidate = true;
            } else {
                f = group_fragments[gi];
                in_group = true;
                in_candidate = true;
            }
            const StateCount* base_first = nullptr;
            const StateCount* base_last = nullptr;
            if (in_group) {
                base_first = group.unaffected.states.data() + group.unaffected.starts[gi];
                base_last = group.unaffected.states.data() + group.unaffected.starts[gi + 1];
            } else {
                const auto aggregate_key = std::make_pair(incumbent, f);
                auto aggregate = state.full_aggregates.find(aggregate_key);
                if (aggregate == state.full_aggregates.end()) {
                    state.key_scratch.clear();
                    for (const HaploidOrigin& origin : census.origins_by_fragment[f]) {
                        StateKey key;
                        key.insert_bp = origin.insert_bp;
                        key.orientation = origin.orientation;
                        key.log_latent_mass = origin.log_latent_mass;
                        state.key_scratch.push_back(key);
                    }
                    std::vector<StateCount> whole;
                    sort_and_count(state.key_scratch, whole, 0);
                    aggregate = state.full_aggregates.emplace(aggregate_key, std::move(whole)).first;
                }
                base_first = aggregate->second.data();
                base_last = aggregate->second.data() + aggregate->second.size();
            }
            const StateCount* add_first = nullptr;
            const StateCount* add_last = nullptr;
            if (in_candidate) {
                add_first = candidate.states.data() + candidate.starts[ci];
                add_last = candidate.states.data() + candidate.starts[ci + 1];
            }
            merge_sorted(base_first, base_last, add_first, add_last, state.merged);
            if (state.merged.empty()) {
                column[f] = kNegInf;
            } else {
                const auto reduce_mark = Clock::now();
                column[f] = reduce_sorted(state.merged.data(),
                                          state.merged.data() + state.merged.size(),
                                          state.library_insert);
                reduce_seconds += seconds_since(reduce_mark);
                // Executed work only: copied fragments formed no states.
                state.cumulative.statistical_state_classes += state.merged.size();
            }
            ++reduced;
            if (checking) merged_for_check.emplace(f, state.merged);
            if (in_group) ++gi;
            if (in_candidate) ++ci;
        }
        local.fragments_reduced = reduced;
        local.fragments_copied = fragment_count - reduced;
        local.seconds_reduce = reduce_seconds;
        local.seconds_aggregate_merge = seconds_since(merge_mark) - reduce_seconds;
        local.seconds_statekey_merge = local.seconds_aggregate_merge + reduce_seconds;

        // ---- 6. certification against the stage-5 std::map merge over a FRESH, unpruned enumeration
        if (checking) {
            std::map<StateKey, std::uint64_t> states;
            for (std::size_t f = 0; f < fragment_count; ++f) {
                states.clear();
                for (const HaploidOrigin& origin : census.origins_by_fragment[f]) {
                    if (affected(origin)) continue;
                    StateKey key;
                    key.insert_bp = origin.insert_bp;
                    key.orientation = origin.orientation;
                    key.log_latent_mass = origin.log_latent_mass;
                    ++states[key];
                }
                for (std::size_t at = state.window_scratch.offsets[f];
                     at < state.window_scratch.offsets[f + 1]; ++at) {
                    const HaploidOrigin& origin = state.window_scratch.origins[at];
                    if (!affected_new(origin)) continue;
                    StateKey key;
                    key.insert_bp = origin.insert_bp;
                    key.orientation = origin.orientation;
                    key.log_latent_mass = origin.log_latent_mass;
                    ++states[key];
                }
                const double expected = states.empty() ? kNegInf : reduce_states(states, state.library_insert);
                if (!(column[f] == expected)) {
                    throw std::runtime_error("origin-path self-check: mass diverged at fragment " +
                                             std::to_string(f));
                }
                const auto merged = merged_for_check.find(static_cast<std::uint32_t>(f));
                if (merged == merged_for_check.end()) continue;   // copied: its mass was checked
                if (merged->second.size() != states.size()) {
                    throw std::runtime_error("origin-path self-check: state count diverged at fragment " +
                                             std::to_string(f));
                }
                auto expected_state = states.begin();
                for (const StateCount& mine : merged->second) {
                    if (!same_state(mine.key, expected_state->first) ||
                        mine.count != expected_state->second) {
                        throw std::runtime_error(
                            "origin-path self-check: key or multiplicity diverged at fragment " +
                            std::to_string(f));
                    }
                    ++expected_state;
                }
            }
        }

        // ---- 7. the promoted winner's census: unaffected (moved) + affected (translated records)
        local.census_payload_bytes_predicted = census.payload_bytes();
        if (retain_census_for_it) {
            HaploidOriginCensus fresh;
            fresh.haplotype_length = new_length;
            fresh.origins_by_fragment.assign(fragment_count, {});
            for (std::size_t f = 0; f < fragment_count; ++f) {
                for (const HaploidOrigin& origin : census.origins_by_fragment[f]) {
                    if (affected(origin)) continue;
                    HaploidOrigin moved = origin;
                    if (origin.forward_start >= block_end) {
                        moved.forward_start = static_cast<std::uint32_t>(
                            static_cast<std::ptrdiff_t>(origin.forward_start) + shift);
                    }
                    fresh.origins_by_fragment[f].push_back(moved);
                }
            }
            for (std::size_t i = 0; i < source->fragments.size(); ++i) {
                const std::uint32_t f = source->fragments[i];
                for (std::uint32_t at = source->starts[i]; at < source->starts[i + 1]; ++at) {
                    const CachedOrigin& record = source->origins[at];
                    HaploidOrigin origin;
                    origin.forward_start = static_cast<std::uint32_t>(base_offset + record.rel_start);
                    origin.insert_bp = record.insert();
                    origin.orientation = record.arm();
                    origin.log_latent_mass = record.log_latent_mass;
                    if (affected_new(origin)) fresh.origins_by_fragment[f].push_back(origin);
                }
            }
            // Census ORDER carries no meaning -- every consumer aggregates it -- so the stage-5 census
            // is compared as a per-fragment multiset.
            if (checking) {
                const auto order = [](const HaploidOrigin& lhs, const HaploidOrigin& rhs) {
                    return std::tie(lhs.forward_start, lhs.insert_bp, lhs.orientation,
                                    lhs.log_latent_mass) <
                           std::tie(rhs.forward_start, rhs.insert_bp, rhs.orientation,
                                    rhs.log_latent_mass);
                };
                for (std::size_t f = 0; f < fragment_count; ++f) {
                    std::vector<HaploidOrigin> expect;
                    for (const HaploidOrigin& origin : census.origins_by_fragment[f]) {
                        if (affected(origin)) continue;
                        HaploidOrigin moved = origin;
                        if (origin.forward_start >= block_end) {
                            moved.forward_start = static_cast<std::uint32_t>(
                                static_cast<std::ptrdiff_t>(origin.forward_start) + shift);
                        }
                        expect.push_back(moved);
                    }
                    for (std::size_t at = state.window_scratch.offsets[f];
                         at < state.window_scratch.offsets[f + 1]; ++at) {
                        const HaploidOrigin& origin = state.window_scratch.origins[at];
                        if (affected_new(origin)) expect.push_back(origin);
                    }
                    std::vector<HaploidOrigin> mine = fresh.origins_by_fragment[f];
                    std::sort(expect.begin(), expect.end(), order);
                    std::sort(mine.begin(), mine.end(), order);
                    bool same = expect.size() == mine.size();
                    for (std::size_t k = 0; same && k < mine.size(); ++k) {
                        same = mine[k].forward_start == expect[k].forward_start &&
                               mine[k].insert_bp == expect[k].insert_bp &&
                               mine[k].orientation == expect[k].orientation &&
                               mine[k].log_latent_mass == expect[k].log_latent_mass;
                    }
                    if (!same) {
                        throw std::runtime_error(
                            "origin-path self-check: promoted census diverged at fragment " +
                            std::to_string(f));
                    }
                }
            }
            local.census_payload_bytes_observed = fresh.payload_bytes();
            state.censuses.emplace(slot, std::move(fresh));
        }

        // ---- 8. work counters: ONLY what was executed. A hit leaves every one of these at zero.
        local.piece_hits = counters.piece_hits;
        local.dp_executions = cache_work.executed_dps;
        local.cache_misses = cache_work.dp_cache_misses;
        local.cache_requests = cache_work.dp_cache_requests;
        local.cache_hits = cache_work.dp_cache_hits;
        state.cumulative.piece_hits += counters.piece_hits;
        state.cumulative.candidate_bases_scanned += counters.candidate_bases_scanned;
        state.cumulative.physical_origins += counters.physical_origins;
        state.cumulative.latent_alignment_pairs += counters.latent_alignment_pairs;
        state.cumulative.valid_fr_states += counters.valid_fr_states;
        state.cumulative.full_read_verifications += counters.full_read_verifications;
        state.cumulative.proposed_starts += counters.proposed_starts;
        state.cumulative.unique_starts += counters.unique_starts;
        state.cumulative.accepted_mate_placements += counters.accepted_mate_placements;
        state.cumulative.dp_evaluations += counters.dp_evaluations;
        state.cumulative.prefilter_rejections += counters.prefilter_rejections;
        state.cumulative_cache.executed_prefilters += cache_work.executed_prefilters;
        state.cumulative_cache.dp_cache_requests += cache_work.dp_cache_requests;
        state.cumulative_cache.dp_cache_hits += cache_work.dp_cache_hits;
        state.cumulative_cache.dp_cache_misses += cache_work.dp_cache_misses;
        state.cumulative_cache.executed_dps += cache_work.executed_dps;
        try {
            enforce_acceptance_caps(state.cumulative, state.options, []() {});
        } catch (const ResourceExhausted&) {
            state.refused = true;
            throw;
        }
        local.seconds = seconds_since(mark);
        state.cumulative_placement_seconds += local.seconds;
        if (work != nullptr) *work = local;
        return slot;
    }

    FragmentScoringCounters counters;
    EmissionCacheCounters cache_work;
    const std::uint64_t rows_before = state.window_scratch.row_constructions;
    const std::uint64_t accumulators_before = state.window_scratch.accumulator_growths;
    const std::uint64_t origins_before = state.window_scratch.origin_growths;
    const auto origins_mark = std::chrono::steady_clock::now();
    try {
        state.enumerate_window(window, state.budget, counters, cache_work, state.window_scratch);
    } catch (const ResourceExhausted&) {
        state.refused = true;
        throw;
    }
    // Certification against the implementation this replaced. Off unless asked for: it doubles the
    // work. The reference gets a throwaway budget so its origins are not charged twice, and its
    // counters start from zero so they can be compared as deltas.
    if (window_self_check().load(std::memory_order_relaxed)) {
        PlacementBudget spare;
        spare.max_piece_hits = std::numeric_limits<std::uint64_t>::max();
        spare.max_physical_origins = std::numeric_limits<std::uint64_t>::max();
        spare.publish_batch = std::numeric_limits<std::uint64_t>::max();
        FragmentScoringCounters reference_counters;
        EmissionCacheCounters reference_cache;
        std::vector<std::vector<HaploidOrigin>> reference;
        enumerate_window_origins_reference(state.prepared, state.fragments, state.library,
                                           state.options, state.options.mismatch_band, window, spare,
                                           reference_counters, reference_cache, reference);
        if (reference_counters.latent_alignment_pairs != counters.latent_alignment_pairs ||
            reference_counters.valid_fr_states != counters.valid_fr_states ||
            reference_counters.physical_origins != counters.physical_origins ||
            reference_counters.accepted_mate_placements != counters.accepted_mate_placements) {
            throw std::runtime_error("window enumeration self-check: counters diverged");
        }
        for (std::size_t f = 0; f < state.fragments.size(); ++f) {
            const std::size_t lo = state.window_scratch.offsets[f];
            const std::size_t hi = state.window_scratch.offsets[f + 1];
            if (hi - lo != reference[f].size()) {
                throw std::runtime_error("window enumeration self-check: origin count diverged");
            }
            for (std::size_t i = 0; i < reference[f].size(); ++i) {
                const HaploidOrigin& mine = state.window_scratch.origins[lo + i];
                const HaploidOrigin& theirs = reference[f][i];
                if (mine.forward_start != theirs.forward_start || mine.insert_bp != theirs.insert_bp ||
                    mine.orientation != theirs.orientation ||
                    !(mine.log_latent_mass == theirs.log_latent_mass)) {
                    throw std::runtime_error(
                        "window enumeration self-check: origin diverged at fragment " +
                        std::to_string(f) + " index " + std::to_string(i));
                }
            }
        }
    }
    local.placement_row_constructions = state.window_scratch.row_constructions - rows_before;
    local.accumulator_growths = state.window_scratch.accumulator_growths - accumulators_before;
    local.origin_output_growths = state.window_scratch.origin_growths - origins_before;
    local.seconds_placement_and_origins =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - origins_mark).count();

    std::size_t slot = 0;
    if (already_held) {
        slot = held_slot;   // refilling an existing column, and its census
    } else {
        slot = state.allocate_repair_slot(repaired_hash);
        state.segments.emplace(slot, repaired);
        state.exposure[slot] = exposure_for_length(new_length, state.fragments.size(),
                                                   state.prepared.fragments_by_floor,
                                                   state.library_insert);
    }
    // Promotion to incumbent is the only thing that earns real bytes, and a census request is the
    // signal for it: the next round reads this haplotype's prefix and suffix directly.
    if (retain_census_for_it) {
        const bool was_segmented = state.segments.count(slot) != 0;
        state.materialise(slot);
        local.full_haplotype_bytes_materialised =
            was_segmented ? static_cast<std::uint64_t>(state.sequences[slot].size() +
                                                       state.reversed[slot].size())
                          : 0;
    }

    local.census_payload_bytes_predicted = census.payload_bytes();
    HaploidOriginCensus fresh;
    if (retain_census_for_it) {
        fresh.haplotype_length = new_length;
        fresh.origins_by_fragment.assign(state.fragments.size(), {});
    }

    const auto merge_mark = std::chrono::steady_clock::now();
    std::map<StateKey, std::uint64_t> states;
    for (std::size_t f = 0; f < state.fragments.size(); ++f) {
        states.clear();
        for (const HaploidOrigin& origin : census.origins_by_fragment[f]) {
            if (affected(origin)) continue;
            ++local.origins_reused;
            StateKey key;
            key.insert_bp = origin.insert_bp;
            key.orientation = origin.orientation;
            key.log_latent_mass = origin.log_latent_mass;
            ++states[key];
            if (retain_census_for_it) {
                HaploidOrigin moved = origin;
                if (origin.forward_start >= block_end) {
                    moved.forward_start = static_cast<std::uint32_t>(
                        static_cast<std::ptrdiff_t>(origin.forward_start) + shift);
                }
                fresh.origins_by_fragment[f].push_back(moved);
            }
        }
        for (std::size_t at = state.window_scratch.offsets[f];
             at < state.window_scratch.offsets[f + 1]; ++at) {
            const HaploidOrigin& origin = state.window_scratch.origins[at];
            if (!affected_new(origin)) continue;
            ++local.origins_recomputed;
            StateKey key;
            key.insert_bp = origin.insert_bp;
            key.orientation = origin.orientation;
            key.log_latent_mass = origin.log_latent_mass;
            ++states[key];
            if (retain_census_for_it) fresh.origins_by_fragment[f].push_back(origin);
        }
        if (states.empty()) continue;
        state.cumulative.statistical_state_classes += states.size();
        state.mass_column[slot][f] = reduce_states(states, state.library_insert);
    }
    local.seconds_statekey_merge =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - merge_mark).count();

    if (retain_census_for_it) {
        local.census_payload_bytes_observed = fresh.payload_bytes();
        state.censuses.emplace(slot, std::move(fresh));
    }
    local.piece_hits = counters.piece_hits;
    local.dp_executions = cache_work.executed_dps;
    local.cache_misses = cache_work.dp_cache_misses;
    local.cache_requests = cache_work.dp_cache_requests;
    local.cache_hits = cache_work.dp_cache_hits;
    state.cumulative.piece_hits += counters.piece_hits;
    state.cumulative.candidate_bases_scanned += counters.candidate_bases_scanned;
    state.cumulative.physical_origins += counters.physical_origins;
    state.cumulative.latent_alignment_pairs += counters.latent_alignment_pairs;
    state.cumulative.valid_fr_states += counters.valid_fr_states;
    state.cumulative.full_read_verifications += counters.full_read_verifications;
    state.cumulative.proposed_starts += counters.proposed_starts;
    state.cumulative.unique_starts += counters.unique_starts;
    state.cumulative.accepted_mate_placements += counters.accepted_mate_placements;
    state.cumulative.dp_evaluations += counters.dp_evaluations;
    state.cumulative.prefilter_rejections += counters.prefilter_rejections;
    state.cumulative_cache.executed_prefilters += cache_work.executed_prefilters;
    state.cumulative_cache.dp_cache_requests += cache_work.dp_cache_requests;
    state.cumulative_cache.dp_cache_hits += cache_work.dp_cache_hits;
    state.cumulative_cache.dp_cache_misses += cache_work.dp_cache_misses;
    state.cumulative_cache.executed_dps += cache_work.executed_dps;
    try {
        enforce_acceptance_caps(state.cumulative, state.options, []() {});
    } catch (const ResourceExhausted&) {
        state.refused = true;
        throw;
    }
    local.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - mark).count();
    state.cumulative_placement_seconds += local.seconds;
    if (work != nullptr) *work = local;
    return slot;
}

void PreparedScoringContext::release_transient(const std::vector<std::size_t>& keep) {
    State& state = *state_;
    const std::set<std::size_t> retained(keep.begin(), keep.end());
    // Group preaggregates and full aggregates are keyed by incumbent and are pure working state:
    // an incumbent that is not kept will never be delta'd from again. The origin cache is content-
    // keyed and deliberately survives -- that is what lets the next round reuse this one's windows.
    for (auto it = state.groups.begin(); it != state.groups.end();) {
        if (retained.count(it->first.incumbent) == 0) {
            it = state.groups.erase(it);
        } else {
            ++it;
        }
    }
    for (auto it = state.full_aggregates.begin(); it != state.full_aggregates.end();) {
        if (retained.count(it->first.first) == 0) {
            it = state.full_aggregates.erase(it);
        } else {
            ++it;
        }
    }
    for (std::size_t slot = state.transient_watermark; slot < state.sequences.size(); ++slot) {
        if (retained.count(slot) != 0) continue;
        if (state.sequences[slot].empty() && state.mass_column[slot].empty() &&
            state.segments.count(slot) == 0) {
            continue;
        }
        release_slot_storage(slot);
    }
}

void PreparedScoringContext::release_slot_storage(std::size_t slot) {
    State& state = *state_;
    // Safe because cache keys OWN their bytes now: entries this haplotype produced keep serving
    // after it is gone. That is what removing the freeze bought.
    const auto range = state.by_hash.equal_range(state.hash_by_slot[slot]);
    for (auto it = range.first; it != range.second; ++it) {
        if (it->second == slot) { state.by_hash.erase(it); break; }
    }
    state.segments.erase(slot);
    std::string().swap(state.sequences[slot]);
    std::string().swap(state.reversed[slot]);
    std::vector<double>().swap(state.mass_column[slot]);
    state.censuses.erase(slot);
    state.free_slots.push_back(slot);
}

std::vector<std::size_t> PreparedScoringContext::live_transient_slots() const {
    const State& state = *state_;
    std::vector<std::size_t> out;
    for (std::size_t slot = state.transient_watermark; slot < state.sequences.size(); ++slot) {
        if (state.sequences[slot].empty() && state.mass_column[slot].empty() && state.segments.count(slot) == 0) continue;
        out.push_back(slot);
    }
    return out;
}

void PreparedScoringContext::release_slots(const std::vector<std::size_t>& slots) {
    State& state = *state_;
    const std::set<std::size_t> named(slots.begin(), slots.end());
    for (const std::size_t slot : named) {
        if (slot < state.transient_watermark || slot >= state.sequences.size() ||
            (state.sequences[slot].empty() && state.mass_column[slot].empty() && state.segments.count(slot) == 0)) {
            throw std::invalid_argument("release_slots: a named slot is not a live transient haplotype");
        }
    }
    for (auto it = state.groups.begin(); it != state.groups.end();) {
        it = named.count(it->first.incumbent) != 0 ? state.groups.erase(it) : std::next(it);
    }
    for (auto it = state.full_aggregates.begin(); it != state.full_aggregates.end();) {
        it = named.count(it->first.first) != 0 ? state.full_aggregates.erase(it) : std::next(it);
    }
    for (const std::size_t slot : named) release_slot_storage(slot);
}

std::uint64_t PreparedScoringContext::transient_payload_bytes() const {
    const State& state = *state_;
    std::uint64_t total = 0;
    for (std::size_t slot = state.transient_watermark; slot < state.sequences.size(); ++slot) {
        // A released slot keeps its shell but owns no heap. std::string::capacity() still reports the
        // short-string buffer for an empty string, so counting it would report storage that was in
        // fact returned -- and this figure exists precisely to show that it was.
        if (state.sequences[slot].empty() && state.mass_column[slot].empty() &&
            state.segments.count(slot) == 0) {
            continue;
        }
        total += state.sequences[slot].capacity() + state.reversed[slot].capacity();
        const auto segment = state.segments.find(slot);
        if (segment != state.segments.end()) {
            total += segment->second.replacement.capacity() + sizeof(SegmentedSequence);
        }
        total += static_cast<std::uint64_t>(state.mass_column[slot].capacity()) * sizeof(double);
        const auto census = state.censuses.find(slot);
        if (census != state.censuses.end()) total += census->second.payload_bytes();
    }
    return total;
}

std::uint64_t PreparedScoringContext::slot_payload_bytes(std::size_t slot) const {
    const State& state = *state_;
    if (slot >= state.sequences.size()) throw std::invalid_argument("payload requested for an unknown haplotype");
    if (state.sequences[slot].empty() && state.mass_column[slot].empty() && state.segments.count(slot) == 0) return 0;
    std::uint64_t total = state.sequences[slot].capacity() + state.reversed[slot].capacity();
    const auto segment = state.segments.find(slot);
    if (segment != state.segments.end()) total += segment->second.replacement.capacity() + sizeof(SegmentedSequence);
    total += static_cast<std::uint64_t>(state.mass_column[slot].capacity()) * sizeof(double);
    const auto census = state.censuses.find(slot);
    if (census != state.censuses.end()) total += census->second.payload_bytes();
    return total;
}

std::uint64_t PreparedScoringContext::retained_census_bytes() const {
    std::uint64_t total = 0;
    for (const auto& item : state_->censuses) total += item.second.payload_bytes();
    return total;
}


CertifiedRanking certify_diplotype_ranking(
    const InternedHaplotypes& haplotypes,
    const std::vector<HaplotypeCandidate>& candidates,
    const std::vector<ReadPair>& fragments,
    const LibraryModel& library,
    const FragmentScoringOptions& options,
    const CertificationOptions& certification,
    const std::vector<double>& exact_prior) {
    if (!exact_prior.empty() && exact_prior.size() != candidates.size()) {
        throw std::invalid_argument("exact prior must be parallel to the candidate list");
    }
    if (candidates.empty()) {
        throw std::invalid_argument("certification needs at least one candidate to rank");
    }
    const FragmentScoringResult scored =
        score_diplotypes(haplotypes, candidates, fragments, library, options);
    return certify_scored_ranking(scored.scores, fragments.size(), options.mismatch_band, scored.counters,
                                  certification, exact_prior);
}

CertifiedRanking certify_scored_ranking(
    const std::vector<DiplotypeScore>& scores,
    std::size_t fragment_count,
    std::size_t mismatch_band,
    const FragmentScoringCounters& counters,
    const CertificationOptions& certification,
    const std::vector<double>& exact_prior) {
    if (!exact_prior.empty() && exact_prior.size() != scores.size()) {
        throw std::invalid_argument("exact prior must be parallel to the candidate list");
    }
    if (scores.empty()) {
        throw std::invalid_argument("certification needs at least one candidate to rank");
    }
    if (!std::isfinite(certification.minimum_tie_tolerance_nats) ||
        certification.minimum_tie_tolerance_nats < 0.0) {
        throw std::invalid_argument("minimum tie tolerance must be finite and non-negative");
    }
    for (double prior : exact_prior) {
        if (!std::isfinite(prior)) {
            throw std::invalid_argument("an exact prior must be finite to rank its candidate");
        }
    }

    CertifiedRanking result;
    result.counters = counters;
    result.mismatch_band = mismatch_band;
    result.candidates.resize(scores.size());
    const double terms = static_cast<double>(fragment_count + 16);
    double log_weight = kNegInf;
    for (std::size_t c = 0; c < scores.size(); ++c) {
        CertifiedCandidate& entry = result.candidates[c];
        entry.candidate_id = c;
        entry.prior = exact_prior.empty() ? 0.0 : exact_prior[c];
        entry.fragment_term = scores[c].fragment_term;
        entry.exposure_term = scores[c].exposure_term;
        entry.score = scores[c].log_likelihood + entry.prior;
        entry.fragments_with_mass = scores[c].fragments_with_mass;
        entry.numeric_error_bound = terms * std::numeric_limits<double>::epsilon() *
                                    (scores[c].absolute_term_sum + std::abs(entry.prior));
        entry.numeric_error_bound =
            std::max(entry.numeric_error_bound, 0.5 * certification.minimum_tie_tolerance_nats);
        log_weight = log_add_exp(log_weight, entry.score);
    }
    result.candidate_set_log_weight = log_weight;
    for (CertifiedCandidate& entry : result.candidates) {
        entry.candidate_set_probability = std::exp(entry.score - log_weight);
    }

    std::sort(result.candidates.begin(), result.candidates.end(),
              [](const CertifiedCandidate& lhs, const CertifiedCandidate& rhs) {
                  if (lhs.score != rhs.score) return lhs.score > rhs.score;
                  return lhs.candidate_id < rhs.candidate_id;
              });

    const CertifiedCandidate& best = result.candidates.front();
    if (result.candidates.size() == 1) {
        result.status = RescoreStatus::Complete;
        result.winner_set = {best.candidate_id};
        result.certified_gap = std::numeric_limits<double>::infinity();
        result.score_gap_to_runner_up = std::numeric_limits<double>::infinity();
        return result;
    }

    result.score_gap_to_runner_up = best.score - result.candidates[1].score;
    const double best_lower = best.score - best.numeric_error_bound;
    double certified_gap = std::numeric_limits<double>::infinity();
    std::vector<std::size_t> inseparable = {best.candidate_id};
    for (std::size_t i = 1; i < result.candidates.size(); ++i) {
        const CertifiedCandidate& rival = result.candidates[i];
        const double margin = best_lower - (rival.score + rival.numeric_error_bound);
        certified_gap = std::min(certified_gap, margin);
        if (margin <= 0.0) inseparable.push_back(rival.candidate_id);
    }
    if (inseparable.size() == 1) {
        result.status = RescoreStatus::Complete;
        result.winner_set = inseparable;
        result.certified_gap = certified_gap;
    } else {
        // The best candidate cannot be separated from these by more than floating-point error. Taking
        // the highest point estimate here would present an uncertified pick as a call.
        result.status = RescoreStatus::RescoreIncomplete;
        result.winner_set = inseparable;
        result.certified_gap = 0.0;
    }
    return result;
}

// ---- EXPERIMENTAL (I1): shared-content placement index (see the header) ----------------------------------------
namespace {
}  // namespace

struct ContentPlacementIndex::State {
    std::size_t max_read = 0, gap = 0, window = 0;
    std::deque<IndexUnit> units;
    std::unordered_multimap<std::uint64_t, std::uint32_t> by_hash;
    std::vector<std::uint32_t> lengths;
    std::vector<std::vector<IndexOccurrence>> occurrences[2];   // [strand][slot]
    ContentIndexReport report;
    // I1.2: segment entities -- a distinct interior segment as its forward and reverse interior units -- and their
    // occurrences per slot (forward offset, ascending), for interior-origin reuse.
    struct Entity {
        std::uint32_t forward_unit = 0, reverse_unit = 0, length = 0, occurrences = 0;
    };
    std::vector<Entity> entities;
    std::unordered_map<std::uint32_t, std::uint32_t> entity_of_forward_unit;
    std::vector<std::vector<std::pair<std::uint32_t, std::uint32_t>>> entity_occurrences;   // [slot] -> (entity, offset)
    // I1.4: fine discovery. Each coarse interior unit is discovered through its content-defined fine units (fine unit,
    // offset of its text in the segment), and its own placement rows are reconstructed from theirs.
    bool fine_discovery = false;
    std::unordered_map<std::uint32_t, std::vector<std::pair<std::uint32_t, std::uint32_t>>> fine_of_coarse;
};

ContentPlacementIndex::ContentPlacementIndex(std::size_t max_read_bp, std::size_t max_gap_bp, bool fine_discovery)
    : state_(std::make_unique<State>()) {
    if (max_read_bp == 0) throw std::invalid_argument("the content index needs a positive maximum read length");
    state_->fine_discovery = fine_discovery;
    state_->max_read = max_read_bp;
    state_->gap = max_gap_bp;
    state_->window = max_read_bp + 2 * max_gap_bp;
}
ContentPlacementIndex::~ContentPlacementIndex() = default;
void set_content_index_hash_collapse(bool collapse) { content_index_hash_collapse().store(collapse); }
void set_fine_window_row_check(bool on) { fine_window_row_check().store(on); }
void set_fine_window_mutation(int mutation) { fine_window_mutation().store(mutation); }
std::size_t ContentPlacementIndex::haplotype_count() const { return state_->lengths.size(); }
std::size_t ContentPlacementIndex::unit_count() const { return state_->units.size(); }
ContentIndexReport ContentPlacementIndex::report() const { return state_->report; }

std::size_t ContentPlacementIndex::add_haplotype(const std::string& bases, const std::vector<std::uint32_t>& segment_lengths) {
    State& st = *state_;
    const auto mark = std::chrono::steady_clock::now();
    std::uint64_t sum = 0;
    for (const std::uint32_t l : segment_lengths) sum += l;
    if (sum != bases.size()) throw std::invalid_argument("content index: segment lengths must sum to the haplotype length");
    const std::size_t slot = st.lengths.size();
    st.lengths.push_back(static_cast<std::uint32_t>(bases.size()));
    const std::size_t n = bases.size();
    const auto intern_unit = [&](IndexUnit&& unit, bool& fresh) {
        fresh = false;
        const std::uint64_t key = unit_hash(unit);
        const auto range = st.by_hash.equal_range(key);
        for (auto it = range.first; it != range.second; ++it) {
            if (same_unit(st.units[it->second], unit)) return it->second;
            ++st.report.hash_collisions;
        }
        fresh = true;
        const std::uint32_t id = static_cast<std::uint32_t>(st.units.size());
        if (unit.fine) {
            ++st.report.fine_units;
            st.report.fine_unit_bytes += unit.bytes.size();
            st.report.fine_core_bytes += unit.core_hi - unit.core_lo;
        } else {
            (unit.interior ? st.report.interior_units : st.report.window_units) += 1;
            (unit.interior ? st.report.coarse_interior_bytes : st.report.window_bytes) += unit.bytes.size();
        }
        st.report.unit_bytes += unit.bytes.size();
        st.units.push_back(std::move(unit));
        st.by_hash.emplace(key, id);
        return id;
    };
    // I1.4: a fresh coarse interior unit's fine units: each content-defined core with g bytes before it and
    // max_read + g after it (clipped to the segment), so every owned start's dependency span [s - g, s + L + g) that
    // lies in the segment lies in the fine text.
    const auto partition = [&](std::uint32_t coarse) {
        const std::string bytes = st.units[coarse].bytes;
        const std::uint8_t strand = st.units[coarse].strand;
        std::vector<std::pair<std::uint32_t, std::uint32_t>> pieces;
        for_each_fine_piece(bytes, strand, st.gap, st.max_read, [&](IndexUnit&& f, std::uint32_t t0) {
            bool fresh = false;
            pieces.push_back({intern_unit(std::move(f), fresh), t0});
        });
        st.report.fine_occurrences += pieces.size();
        st.fine_of_coarse[coarse] = std::move(pieces);
    };
    const auto intern = [&](IndexUnit&& unit) {
        const bool interior = unit.interior;
        bool fresh = false;
        const std::uint32_t id = intern_unit(std::move(unit), fresh);
        if (fresh && interior && st.fine_discovery) partition(id);
        return id;
    };
    std::vector<std::pair<std::uint32_t, std::uint32_t>> entity_occ;
    for (std::uint8_t strand = 0; strand < 2; ++strand) {
        const std::string text = strand == 0 ? bases : std::string(bases.rbegin(), bases.rend());
        std::vector<std::uint32_t> lengths = segment_lengths;
        if (strand == 1) std::reverse(lengths.begin(), lengths.end());
        std::vector<IndexOccurrence> occ;
        std::vector<std::size_t> cuts{0};
        std::size_t at = 0;
        for (const std::uint32_t l : lengths) {
            // an interior can keep a start only when the segment holds g + 1 + g bytes; shorter ones are covered by windows
            if (l > 2 * st.gap) {
                IndexUnit u;
                u.strand = strand;
                u.interior = true;
                u.bytes = text.substr(at, l);
                const std::string forward_bytes = strand == 0 ? u.bytes : std::string();
                const std::uint32_t id = intern(std::move(u));
                occ.push_back({id, static_cast<std::uint32_t>(at)});
                ++st.report.segment_occurrences;
                if (strand == 0) {
                    // the segment's reverse interior unit: the same key the strand-1 pass interns, so the same id
                    IndexUnit r;
                    r.strand = 1;
                    r.interior = true;
                    r.bytes.assign(forward_bytes.rbegin(), forward_bytes.rend());
                    const std::uint32_t rid = intern(std::move(r));
                    auto [it, fresh] = st.entity_of_forward_unit.emplace(id, static_cast<std::uint32_t>(st.entities.size()));
                    if (fresh) st.entities.push_back({id, rid, l, 0});
                    ++st.entities[it->second].occurrences;
                    entity_occ.emplace_back(it->second, static_cast<std::uint32_t>(at));
                }
            }
            at += l;
            cuts.push_back(at);
        }
        cuts.push_back(n);
        std::sort(cuts.begin(), cuts.end());
        cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());
        for (const std::size_t b : cuts) {
            std::size_t lo = 0;
            IndexUnit u = boundary_window_unit(text, strand, b, st.window, lo);
            occ.push_back({intern(std::move(u)), static_cast<std::uint32_t>(lo)});
            ++st.report.window_occurrences;
        }
        st.report.occurrence_bytes += occ.size() * sizeof(IndexOccurrence);
        st.occurrences[strand].push_back(std::move(occ));
    }
    st.entity_occurrences.push_back(std::move(entity_occ));
    ++st.report.haplotypes;
    st.report.seconds_registration += std::chrono::duration<double>(std::chrono::steady_clock::now() - mark).count();
    return slot;
}

struct ContentIndexBatch::State {
    const ContentPlacementIndex* index = nullptr;
    std::vector<ReadPair> fragments;
    LibraryModel library;
    FragmentScoringOptions options;
    DiscreteInsertDistribution library_insert;
    PreparedFragmentSet prepared;
    std::vector<std::vector<StoredPlacement>> placements;   // per unit, sorted by (form, start, consumed)
    std::size_t evaluated = 0;
    ContentBatchReport report;
    // I1.1: the exact prefilter-decision cache, owned by this read batch and kept across discover() calls
    std::unique_ptr<UnionPrefilterCache> prefilter;
    // I1.2: interior-origin aggregates, one per admitted segment entity, owned by this batch (built from its own
    // verified unit placements); entity_entry: -1 not built, -2 declined (budget), -3 ineligible, else the entry index.
    std::deque<UnionOriginEntry> origin_entries;
    std::vector<long> entity_entry;
    std::uint64_t origin_held = 0;
    std::size_t zone_read = 0;   // R: the batch's longest read; zones need it
    // I1.3: completed calculations, retained within this immutable read/library/scoring context (append-only).
    // A slot's mass column and exposure depend only on its own sequence, the reads, the library and the options --
    // never on which other slots exist -- so a completed slot is never recomputed. Pair likelihood components (no
    // priors, normalisers or selections) are kept per unordered slot pair: log_add_exp and the exposure sum are
    // symmetric, so (a, b) and (b, a) are the same bits.
    struct SlotAccounting {
        std::uint64_t expanded = 0, removed = 0, plain_origins = 0, direct = 0, served = 0, classes = 0, zones = 0;
    };
    struct PairResult {
        DiplotypeScore score;
        std::vector<DiplotypeScore> audits;
    };
    std::vector<std::vector<double>> mass_column;
    std::vector<double> exposure;
    std::vector<char> completed;
    std::vector<SlotAccounting> accounting;
    std::unordered_map<std::uint64_t, PairResult> pairs;
    bool retain = true;
    // I2-A: the placement store's declared byte budget (0: none), the one-shot lifecycle, and discovery's executed
    // seeding work (cumulative), which every score() reports with its own pairing work.
    std::uint64_t store_budget = 0;
    bool one_shot = false, sealed = false;
    FragmentScoringCounters discovery_counters;
    // I2-C: the placement-store ledger. It covers exactly the per-unit placement rows' ALLOCATED capacity and the
    // per-unit container's own capacity (the `placements` vector of vectors) -- nothing else: not the index's unit bytes
    // and occurrence tables, the prepared reads, the prefilter or emission caches, origin aggregates, mass columns, or the
    // workers' temporary scratch (reported apart). It is a container ledger, not process memory. Exact layout (default):
    // every growth's REQUESTED capacity is reserved against the budget before the storage is allocated (lock-free,
    // overflow-checked; a reservation that would cross the budget refuses). Capacity the vector returns beyond its request
    // is charged right AFTER the allocation exists (and refuses then), so the ledger may briefly trail the allocator by
    // that excess; none was observed with exact reserves. The tested guarantee: at every observation point the ledger
    // equals the recounted capacity, and the budget accepts at the exact peak and refuses one byte below it. Rows are
    // copied into exact-capacity vectors, coarse rows are
    // built from their known counts, and a one-shot batch releases each fine unit's rows once the last coarse unit using
    // it is built. Legacy layout (the I2-A behaviour, kept for comparison): rows grow by push_back and the budget charges
    // occupied rows; the ledger then only measures capacity.
    bool exact_store = true;
    // I4-B (EXPERIMENTAL): coarse interior units served as views over their fine pieces (index-owned fine_of_coarse
    // entries); view_pieces[u] is null for every other unit. Its capacity is charged to the store ledger (metadata).
    bool fine_views = false;
    std::vector<const std::vector<std::pair<std::uint32_t, std::uint32_t>>*> view_pieces;
    std::uint64_t view_table_charged = 0;
    std::atomic<std::uint64_t> zone_live{0}, zone_peak{0}, zone_total{0};
    std::atomic<std::uint64_t> store_live{0}, store_peak{0}, store_allocated{0}, store_released{0};
    std::uint64_t outer_charged = 0;
    // S2-C (EXPERIMENTAL): the demand-driven route. A batch is either bulk (discover/score) or selective
    // (discover_for/score_selected), never both. unit_done: per unit, discovered (directly, or a coarse interior built
    // or viewed from its fine units). Cumulative charges across every selective call, so no cap resets between requests.
    bool selective = false;
    std::vector<char> unit_done;
    std::uint64_t selective_origins = 0, selective_placements_consumed = 0, selective_origins_consumed = 0,
                  selective_states_consumed = 0;
    FragmentScoringCounters selective_score_counters;   // executed pairing and aggregate work, summed over calls
    // S2-D (EXPERIMENTAL): bounded emission reuse across the selective discovery calls of this batch. Between calls the
    // cache is detached into `carried` (owned here, destroyed with the batch -- never shared across samples); each call
    // adopts it, under the same byte budget, only if the library, band and recursion backend match (the existing
    // handoff's identity check; keys carry the read and window bytes). poisoned: a refused or failed selective call
    // ends the batch -- no later call may continue on incomplete evidence.
    bool persistent_emissions = false;
    EmissionCacheHandoff carried;
    bool poisoned = false;
    // S2-E: the running totals of the blocks before this batch, and build_columns' no-export mode
    bool baseline_set = false;
    UnionCumulativeCharges baseline;
    bool materialise = true;
    bool reserve_store(std::uint64_t bytes, bool enforce) {
        std::uint64_t cur = store_live.load(std::memory_order_relaxed);
        for (;;) {
            if (cur > std::numeric_limits<std::uint64_t>::max() - bytes) return false;   // overflow
            if (enforce && store_budget != 0 && (bytes > store_budget || cur > store_budget - bytes)) return false;
            if (store_live.compare_exchange_weak(cur, cur + bytes, std::memory_order_relaxed)) break;
        }
        store_allocated.fetch_add(bytes, std::memory_order_relaxed);
        const std::uint64_t now = cur + bytes;
        std::uint64_t peak = store_peak.load(std::memory_order_relaxed);
        while (now > peak && !store_peak.compare_exchange_weak(peak, now, std::memory_order_relaxed)) {
        }
        return true;
    }
    void release_store(std::uint64_t bytes) {
        store_live.fetch_sub(bytes, std::memory_order_relaxed);
        store_released.fetch_add(bytes, std::memory_order_relaxed);
    }
};
// TEST ONLY: the exact layout charges occupied rows instead of allocated capacity (a negative control the capacity gate
// must catch).
std::atomic<bool>& store_charge_rows_only() {
    static std::atomic<bool> on{false};
    return on;
}
// TEST ONLY: request this many spare slots per stored vector (unused capacity the ledger must charge).
std::atomic<std::size_t>& store_spare_capacity() {
    static std::atomic<std::size_t> spare{0};
    return spare;
}
// The capacity a stored vector of n rows is requested with, and what the ledger charges for it before allocation.
std::uint64_t requested_slots(std::uint64_t n) { return n + store_spare_capacity().load(std::memory_order_relaxed); }
std::uint64_t charged_bytes(std::uint64_t n) {
    return (store_charge_rows_only().load(std::memory_order_relaxed) ? n : requested_slots(n)) * sizeof(StoredPlacement);
}
void ContentIndexBatch::set_retention(bool retain) { state_->retain = retain; }
void ContentIndexBatch::set_placement_store_budget(std::size_t bytes) {
    state_->store_budget = bytes;
    state_->report.store_budget_bytes = bytes;
}
void ContentIndexBatch::set_one_shot(bool one_shot) { state_->one_shot = one_shot; }
void ContentIndexBatch::set_legacy_store_layout(bool legacy) { state_->exact_store = !legacy; }
void ContentIndexBatch::set_fine_views(bool on) {
    state_->fine_views = on;
    state_->report.evidence_layout = on ? "fine_views" : "reconstructed";
}
void set_store_charge_rows_only(bool on) { store_charge_rows_only().store(on); }
void set_store_spare_capacity(std::size_t slots) { store_spare_capacity().store(slots); }
std::uint64_t ContentIndexBatch::store_capacity_recount() const {
    std::uint64_t b = state_->placements.capacity() * sizeof(std::vector<StoredPlacement>);
    for (const auto& v : state_->placements) b += v.capacity() * sizeof(StoredPlacement);
    b += state_->view_pieces.capacity() * sizeof(state_->view_pieces[0]);   // I4-B: the view table
    return b;
}

namespace {
// I2-A: a worker that stopped because another crossed a cap throws PlacementCancelled; the cap actually crossed is what
// propagates, so the message does not depend on which worker noticed first.
void rethrow_primary_failure(const std::vector<std::exception_ptr>& failures) {
    std::exception_ptr cancelled;
    for (const std::exception_ptr& f : failures) {
        if (!f) continue;
        try {
            std::rethrow_exception(f);
        } catch (const PlacementCancelled&) {
            cancelled = f;
        } catch (...) {
            throw;
        }
    }
    if (cancelled) std::rethrow_exception(cancelled);
}
void add_counters(FragmentScoringCounters& into, const FragmentScoringCounters& c) {
    into.candidate_bases_scanned += c.candidate_bases_scanned;
    into.piece_hits += c.piece_hits;
    into.proposed_starts += c.proposed_starts;
    into.unique_starts += c.unique_starts;
    into.full_read_verifications += c.full_read_verifications;
    into.dp_evaluations += c.dp_evaluations;
    into.prefilter_rejections += c.prefilter_rejections;
    into.accepted_mate_placements += c.accepted_mate_placements;
    into.valid_fr_states += c.valid_fr_states;
    into.latent_alignment_pairs += c.latent_alignment_pairs;
    into.physical_origins += c.physical_origins;
    into.statistical_state_classes += c.statistical_state_classes;
    into.reused_physical_origins += c.reused_physical_origins;
    into.reused_state_aggregates += c.reused_state_aggregates;
    into.reused_latent_alignment_pairs += c.reused_latent_alignment_pairs;
    into.reused_valid_fr_states += c.reused_valid_fr_states;
    into.merged_state_classes += c.merged_state_classes;
}
}  // namespace

ContentIndexBatch::ContentIndexBatch(const ContentPlacementIndex& index, std::vector<ReadPair> fragments, LibraryModel library,
                                     FragmentScoringOptions options)
    : state_(std::make_unique<State>()) {
    validate_library_model(library);
    if (library.max_gap_bp != index.state_->gap) throw std::invalid_argument("content index: max_gap_bp differs from the index's");
    for (const ReadPair& f : fragments) {
        if (f.first.bases.size() > index.state_->max_read || f.second.bases.size() > index.state_->max_read) {
            throw std::invalid_argument("content index: a read is longer than the index's maximum read length");
        }
    }
    if (!options.record_origins_for_fragments.empty() || options.record_origin_summary || options.union_segment_lengths != nullptr) {
        throw std::invalid_argument("content index: origin recording and union caches are not supported");
    }
    State& st = *state_;
    st.index = &index;
    st.fragments = std::move(fragments);
    st.library = std::move(library);
    st.options = std::move(options);
    st.library_insert = build_library_insert_distribution(st.library);
    st.prepared = prepare_fragment_set(st.fragments, st.library, st.options);
    for (const auto* forms : {&st.prepared.forward_forms, &st.prepared.reverse_forms})
        for (const ReadForm& form : *forms) st.zone_read = std::max(st.zone_read, form.bases.size());
    st.report.origin_reuse_status = st.options.union_origin_cache_bytes > 0 ? "on" : "off";
    st.report.origin_budget_bytes = st.options.union_origin_cache_bytes;
    if (st.options.union_prefilter_cache_bytes > 0) {
        std::size_t length = 0;
        if (single_read_length(st.prepared, length)) {
            st.prefilter = make_union_prefilter_cache(st.prepared, st.library, st.options, st.report.prefilter);
            st.report.prefilter_status = "on";
        } else {
            st.report.prefilter_status = "fallback_uncached_mixed_read_lengths";   // the cache needs one read length
        }
    } else {
        st.report.prefilter_status = "off";
    }
}
ContentIndexBatch::~ContentIndexBatch() = default;
const ContentBatchReport& ContentIndexBatch::report() const { return state_->report; }

void ContentIndexBatch::discover() {
    if (state_->selective) throw std::logic_error("content index: bulk discovery on a selective batch");
    discover_impl(nullptr);
}
void ContentIndexBatch::set_persistent_emission_cache(bool on) {
    if (on && state_->options.emission_handoff_out != nullptr) {
        throw std::invalid_argument("content index: persistent emissions and the adaptive emission handoff are exclusive");
    }
    state_->persistent_emissions = on;
    state_->report.persistent_emission_status = on ? "on" : "off";
}
void ContentIndexBatch::set_cumulative_baseline(const UnionCumulativeCharges& baseline) {
    if (state_->selective || state_->evaluated != 0) throw std::logic_error("content index: a baseline is set before bulk discovery");
    state_->baseline_set = true;
    state_->baseline = baseline;
    if (state_->prefilter) state_->prefilter->consumed.store(baseline.prefilters_consumed, std::memory_order_relaxed);
}
void ContentIndexBatch::build_columns() {
    if (state_->selective) throw std::logic_error("content index: build_columns() on a selective batch");
    state_->materialise = false;
    (void)score_impl({}, false);
    state_->materialise = true;
}
const std::vector<std::vector<double>>& ContentIndexBatch::mass_columns() const { return state_->mass_column; }
UnionCumulativeCharges ContentIndexBatch::charges() const {
    const State& st = *state_;
    UnionCumulativeCharges c;
    c.executed = st.report.counters_last;   // discovery's and the last score call's executed work (this batch only)
    c.piece_hits = st.discovery_counters.piece_hits;
    c.physical_origins = st.report.counters_last.physical_origins;
    c.placements_consumed = st.report.placements_consumed_last - (st.baseline_set ? st.baseline.placements_consumed : 0);
    c.origins_consumed = st.report.origins_consumed_last - (st.baseline_set ? st.baseline.origins_consumed : 0);
    c.states_consumed = st.report.states_consumed_last - (st.baseline_set ? st.baseline.states_consumed : 0);
    if (st.prefilter) {
        c.prefilters_consumed = st.prefilter->consumed.load(std::memory_order_relaxed) - (st.baseline_set ? st.baseline.prefilters_consumed : 0);
    }
    return c;
}
void ContentIndexBatch::release_carried_emissions() {
    state_->report.persistent_carried_bytes = 0;
    state_->report.persistent_carried_entries = 0;
    state_->carried.clear();
}
void ContentIndexBatch::discover_for(const std::vector<std::size_t>& slots) {
    State& st = *state_;
    if (st.poisoned) throw std::logic_error("content index: a refused or failed selective batch cannot continue");
    if (!st.selective && st.evaluated != 0) throw std::logic_error("content index: selective discovery on a bulk batch");
    if (st.one_shot) throw std::invalid_argument("content index: selective discovery needs a batch that is not one-shot");
    if (!st.exact_store) throw std::invalid_argument("content index: selective discovery needs the exact store layout");
    st.selective = true;
    try {
        discover_impl(&slots);
    } catch (...) {
        st.poisoned = true;   // S2-D: incomplete evidence never feeds a later call
        st.carried.clear();
        throw;
    }
}
std::size_t ContentIndexBatch::units_discovered() const {
    std::size_t n = 0;
    for (const char d : state_->unit_done) n += d != 0;
    return state_->selective ? n : state_->evaluated;
}
std::size_t ContentIndexBatch::columns_completed() const {
    std::size_t n = 0;
    for (const char d : state_->completed) n += d != 0;
    return n;
}
void ContentIndexBatch::discover_impl(const std::vector<std::size_t>* slots) {
    State& st = *state_;
    const ContentPlacementIndex::State& ix = *st.index->state_;
    const std::size_t first = st.evaluated, last = ix.units.size();
    // The units evaluated now: `direct` are discovered by proposal and verification; `interior` (fine discovery only)
    // are coarse interiors built or viewed from their fine units. Bulk: every unit added since the last call, in index
    // order (exactly as before). Selective (S2-C): the units the requested slots' occurrences name and, under fine
    // discovery, each coarse interior's fine units -- each at most once per batch.
    std::vector<std::size_t> direct, interior;
    if (slots == nullptr) {
        if (first == last) return;
        for (std::size_t u = first; u < last; ++u) {
            if (ix.fine_discovery && ix.units[u].interior) interior.push_back(u);
            else direct.push_back(u);
        }
    } else {
        st.unit_done.resize(last, 0);
        std::vector<char> want(last, 0);
        for (const std::size_t h : *slots) {
            if (h >= ix.lengths.size()) throw std::invalid_argument("content index: selective discovery names an unknown slot");
            for (int strand = 0; strand < 2; ++strand) {
                for (const IndexOccurrence& occ : ix.occurrences[strand][h]) want[occ.unit] = 1;
            }
        }
        std::vector<char> queued(last, 0);
        for (std::size_t u = 0; u < last; ++u) {
            if (!want[u] || st.unit_done[u]) continue;
            if (ix.fine_discovery && ix.units[u].interior) {
                interior.push_back(u);
                for (const auto& [fine, offset] : ix.fine_of_coarse.at(static_cast<std::uint32_t>(u))) {
                    if (!st.unit_done[fine] && !queued[fine]) {
                        queued[fine] = 1;
                        direct.push_back(fine);
                    }
                }
            } else if (!queued[u]) {
                queued[u] = 1;
                direct.push_back(u);
            }
        }
        std::sort(direct.begin(), direct.end());
        if (direct.empty() && interior.empty()) return;
    }
    if (st.sealed) throw std::logic_error("content index: a one-shot batch was asked to discover new units");
    if (st.fine_views && !st.exact_store) throw std::invalid_argument("content index: fine views need the exact store layout");
    const auto mark = std::chrono::steady_clock::now();
    const char* const kStoreRefusal =
        "the indexed union's placement store reached its declared byte budget; no partial ranking is returned";
    // I2-C: the per-unit container grows to exactly `last` entries, its capacity reserved first (and, in the legacy
    // layout, only measured)
    if (st.placements.capacity() < last) {
        const std::uint64_t grow = (static_cast<std::uint64_t>(last) - st.placements.capacity()) * sizeof(std::vector<StoredPlacement>);
        if (!st.reserve_store(grow, st.exact_store)) throw ResourceExhausted(kStoreRefusal);
        st.placements.reserve(last);
        st.outer_charged += grow;
        if (st.placements.capacity() * sizeof(std::vector<StoredPlacement>) > st.outer_charged) {
            const std::uint64_t extra = st.placements.capacity() * sizeof(std::vector<StoredPlacement>) - st.outer_charged;
            if (!st.reserve_store(extra, st.exact_store)) throw ResourceExhausted(kStoreRefusal);
            st.outer_charged += extra;
        }
    }
    st.placements.resize(last);
    // I2-A: stored rows charged against the declared store budget as each unit completes (a monotone shared total, so
    // the verdict does not depend on the schedule: it refuses exactly when the merged total exceeds the budget).
    // I2-C: the legacy layout only; the exact layout charges allocated capacity through the ledger instead.
    std::atomic<std::uint64_t> store_charged{st.report.store_charged_bytes};
    const auto charge_store = [&](std::uint64_t bytes) {
        const std::uint64_t total = store_charged.fetch_add(bytes, std::memory_order_relaxed) + bytes;
        if (st.store_budget != 0 && total > st.store_budget) {
            throw ResourceExhausted("the indexed union's placement store reached its declared byte budget; no partial "
                                    "ranking is returned");
        }
    };
    begin_emission_cache(st.options.emission_cache_bytes);
    // S2-D: whatever happens below, the process-global cache is never left armed after this call (an exception between
    // here and the explicit end/detach ends it here, discarding its contents)
    struct CacheGuard {
        bool open = true;
        ~CacheGuard() {
            if (open) (void)end_emission_cache();
        }
    } cache_guard;
    const bool persist = slots != nullptr && st.persistent_emissions;
    if (persist && st.carried.holds()) {
        const std::string verdict = adopt_emission_cache(&st.carried, st.library, st.options.mismatch_band,
                                                         st.options.gap_recursion_backend, st.options.emission_cache_bytes, false);
        if (verdict == "none") ++st.report.persistent_adoptions;
        else ++st.report.persistent_adoption_refusals;
    }
    PlacementBudget budget;
    budget.max_piece_hits = st.options.max_piece_hits;
    budget.max_physical_origins = st.options.max_physical_origins;
    budget.publish_batch = std::max<std::uint64_t>(1, st.options.piece_hit_publish_batch);
    if (slots != nullptr) budget.piece_hits.store(st.report.piece_hits, std::memory_order_relaxed);   // S2-C: cumulative
    if (slots == nullptr && st.baseline_set) budget.piece_hits.store(st.baseline.piece_hits, std::memory_order_relaxed);   // S2-E
    const std::size_t work_items = slots == nullptr ? last - first : direct.size();
    const std::size_t threads = std::max<std::size_t>(1, std::min(st.options.threads, std::max<std::size_t>(1, work_items)));
    std::atomic<std::size_t> next{0};
    std::vector<FragmentScoringCounters> counters(threads);
    std::vector<EmissionCacheCounters> cache(threads);
    std::vector<SeedingTimes> times(threads);
    std::vector<std::uint64_t> accepted(threads, 0), kept(threads, 0), bases(threads, 0);
    std::vector<std::uint64_t> skipped(threads, 0), evaluated_units(threads, 0), scratch(threads, 0);
    std::vector<std::exception_ptr> failures(threads);
    const auto work = [&](std::size_t t) {
        const ActiveWorker active;
        std::vector<std::vector<FormPlacement>> rows[2] = {
            std::vector<std::vector<FormPlacement>>(st.prepared.forward_forms.size()),
            std::vector<std::vector<FormPlacement>>(st.prepared.reverse_forms.size())};
        std::vector<std::pair<std::uint32_t, std::uint32_t>> proposed;
        std::vector<std::uint32_t> touched;
        std::vector<StoredPlacement> found;   // I2-C: the worker's reused discovery buffer (temporary, reported apart)
        const auto scratch_bytes = [&]() {
            std::uint64_t b = found.capacity() * sizeof(StoredPlacement) + proposed.capacity() * 8 + touched.capacity() * 4;
            for (const auto& r : rows) {
                b += r.capacity() * sizeof(std::vector<FormPlacement>);
                for (const auto& v : r) b += v.capacity() * sizeof(FormPlacement);
            }
            return b;
        };
        try {
            for (std::size_t i = next++; i < direct.size(); i = next++) {
                if (budget.stop()) throw PlacementCancelled();
                const std::size_t u = direct[i];
                const IndexUnit& unit = ix.units[u];   // I1.4: coarse interiors are not here; reconstructed below
                std::vector<StoredPlacement>& store = st.placements[u];
                if (st.exact_store) {
                    found.clear();
                    accepted[t] += discover_unit_placements(unit, st.prepared, st.prepared.forward_indexes,
                                                            st.prepared.reverse_indexes, budget, st.options.mismatch_band,
                                                            st.library, st.options, rows, proposed, touched, counters[t],
                                                            cache[t], &times[t], st.prefilter.get(), &skipped[t], found);
                    const std::uint64_t n = found.size();
                    if (n > (std::numeric_limits<std::uint64_t>::max() / sizeof(StoredPlacement)) - store_spare_capacity().load()) {
                        throw ResourceExhausted(kStoreRefusal);
                    }
                    if (n != 0 || store_spare_capacity().load() != 0) {
                        // reserved before the persistent vector is allocated, at the capacity it is requested with
                        if (!st.reserve_store(charged_bytes(n), true)) {
                            budget.cancel();
                            throw ResourceExhausted(kStoreRefusal);
                        }
                        std::vector<StoredPlacement> exact;
                        exact.reserve(static_cast<std::size_t>(requested_slots(n)));
                        exact.insert(exact.end(), found.begin(), found.end());
                        store.swap(exact);
                        // anything allocated beyond the request is charged too (unless the negative control drops it)
                        if (!store_charge_rows_only().load(std::memory_order_relaxed) && store.capacity() > requested_slots(n) &&
                            !st.reserve_store((store.capacity() - requested_slots(n)) * sizeof(StoredPlacement), true)) {
                            budget.cancel();
                            throw ResourceExhausted(kStoreRefusal);
                        }
                    }
                } else {
                    accepted[t] += discover_unit_placements(unit, st.prepared, st.prepared.forward_indexes,
                                                            st.prepared.reverse_indexes, budget, st.options.mismatch_band,
                                                            st.library, st.options, rows, proposed, touched, counters[t],
                                                            cache[t], &times[t], st.prefilter.get(), &skipped[t], store);
                    charge_store(static_cast<std::uint64_t>(store.size()) * sizeof(StoredPlacement));
                    (void)st.reserve_store(store.capacity() * sizeof(StoredPlacement), false);   // measured only
                }
                kept[t] += store.size();
                ++evaluated_units[t];
                bases[t] += unit.bytes.size();
            }
        } catch (...) {
            failures[t] = std::current_exception();
            budget.cancel();
        }
        scratch[t] = scratch_bytes();
    };
    std::vector<std::thread> workers;
    for (std::size_t t = 0; t < threads; ++t) workers.emplace_back(work, t);
    for (std::thread& w : workers) w.join();
    ContentBatchReport& r = st.report;
    {
        const EmissionCacheCounters stored = read_emission_cache_counters();
        {   // I2-K (measurement): what discovery's cache actually holds before it is discarded
            EmissionCache& held = emission_cache();
            r.emission_cache_entries_held = r.emission_cache_table_bytes = r.emission_cache_slab_bytes = 0;
            r.emission_cache_component_bytes = 0;
            for (CacheShard& shard : held.shards) {
                const std::lock_guard<std::mutex> guard(shard.guard);
                r.emission_cache_entries_held += shard.occupied;
                r.emission_cache_table_bytes += static_cast<std::uint64_t>(shard.table.capacity()) * sizeof(ArenaEntry);
                r.emission_cache_slab_bytes += static_cast<std::uint64_t>(shard.slabs.size()) * kSlabBytes;
                r.emission_cache_component_bytes += static_cast<std::uint64_t>(shard.arena.capacity()) * sizeof(GappedComponent);
            }
            r.emission_cache_budgeted_bytes = held.budgeted_bytes_used.load(std::memory_order_relaxed);
        }
        st.report.persistent_inherited_hits += emission_cache().inherited_hits.load(std::memory_order_relaxed);   // S2-D
        cache_guard.open = false;
        const bool failed = std::any_of(failures.begin(), failures.end(), [](const std::exception_ptr& f) { return static_cast<bool>(f); });
        if (persist && !failed) {   // S2-D: carried to the batch's next selective discovery (same budget, same identity)
            detach_emission_cache(st.carried, st.library, st.options.mismatch_band, st.options.gap_recursion_backend, false);
            st.report.persistent_carried_bytes = st.carried.allocated_bytes();
            st.report.persistent_carried_entries = st.carried.entries();
        } else if (st.options.emission_handoff_out != nullptr) {   // I2-K (experimental): kept for the adaptive search
            detach_emission_cache(*st.options.emission_handoff_out, st.library, st.options.mismatch_band,
                                  st.options.gap_recursion_backend);
        } else {
            (void)end_emission_cache();
        }
        for (std::size_t t = 0; t < threads; ++t) {
            r.emission_cache.executed_prefilters += cache[t].executed_prefilters;
            r.emission_cache.dp_cache_requests += cache[t].dp_cache_requests;
            r.emission_cache.dp_cache_hits += cache[t].dp_cache_hits;
            r.emission_cache.dp_cache_misses += cache[t].dp_cache_misses;
            r.emission_cache.executed_dps += cache[t].executed_dps;
        }
        r.emission_cache.payload_bytes = stored.payload_bytes;
        r.emission_cache.entries = stored.entries;
        r.emission_cache.entries_not_retained = stored.entries_not_retained;
        r.emission_cache.budgeted_bytes_used = stored.budgeted_bytes_used;
    }
    // I2-A: executed work is merged BEFORE any failure propagates, so a refused pass reports what it did
    for (std::size_t t = 0; t < threads; ++t) add_counters(st.discovery_counters, counters[t]);
    r.piece_hits = std::max(st.discovery_counters.piece_hits, budget.piece_hits.load(std::memory_order_relaxed));
    r.candidate_bases_scanned = st.discovery_counters.candidate_bases_scanned;
    r.store_charged_bytes = st.exact_store ? st.store_allocated.load(std::memory_order_relaxed)
                                           : store_charged.load(std::memory_order_relaxed);
    r.threads_used_discovery = std::max(r.threads_used_discovery, threads);
    {
        std::uint64_t sum = 0;
        for (const std::uint64_t b : scratch) sum += b;
        r.worker_scratch_bytes = std::max(r.worker_scratch_bytes, sum);   // held at once while the workers ran
    }
    const auto note_ledger = [&]() {
        r.store_layout = st.exact_store ? "exact" : "legacy_i2a";
        r.store_live_bytes = st.store_live.load(std::memory_order_relaxed);
        r.store_peak_bytes = st.store_peak.load(std::memory_order_relaxed);
        r.store_allocated_bytes = st.store_allocated.load(std::memory_order_relaxed);
        r.store_released_bytes = st.store_released.load(std::memory_order_relaxed);
    };
    note_ledger();
    r.counters_last = st.discovery_counters;   // what a refusal before any score() reports
    r.counters_last.piece_hits = r.piece_hits;
    r.counters_last.fragments_with_floor_binding = st.prepared.fragments_with_floor_binding;
    rethrow_primary_failure(failures);
    for (std::size_t t = 0; t < threads; ++t) {
        r.proposed_starts += counters[t].proposed_starts;
        r.unique_starts += counters[t].unique_starts;
        r.full_read_verifications += counters[t].full_read_verifications;
        r.dp_evaluations += counters[t].dp_evaluations;
        r.executed_prefilters += cache[t].executed_prefilters;
        r.executed_dps += cache[t].executed_dps;
        r.prefilter_rejections += counters[t].prefilter_rejections;
        r.unit_placements_accepted += accepted[t];
        r.unit_placements_kept += kept[t];
        r.unit_bases_scanned += bases[t];
        r.worker_seconds_scan += static_cast<double>(times[t].scan_ns) * 1e-9;
        r.worker_seconds_sort += static_cast<double>(times[t].sort_ns) * 1e-9;
        r.worker_seconds_verification += static_cast<double>(times[t].verify_ns + times[t].row_sort_ns) * 1e-9;
    }
    r.units_evaluated += slots == nullptr ? last - first : direct.size() + interior.size();
    ++r.discovery_calls;
    for (std::size_t t = 0; t < threads; ++t) {
        r.owned_starts_skipped_before_verification += skipped[t];
        r.units_discovered_directly += evaluated_units[t];
    }
    if (ix.fine_discovery) {
        // I1.4: each new coarse interior unit's rows, exactly as coarse discovery would keep them: the union of its fine
        // units' owned placements, translated to the segment and put in (form, start, consumed) order
        const auto reconstruct_mark = std::chrono::steady_clock::now();
        const auto row_order = [](const StoredPlacement& a, const StoredPlacement& b) {
            return std::tie(a.form, a.start, a.consumed) < std::tie(b.form, b.start, b.consumed);
        };
        if (st.fine_views) {
            // I4-B: no coarse row is copied. Each new coarse interior unit is served as a view over its fine pieces (the
            // index's fine_of_coarse entry); the fine rows -- the distinct, separately verified evidence -- stay held until
            // the batch ends. The view table is store metadata: its capacity is charged before it is allocated.
            if (st.view_pieces.capacity() < last) {
                const std::uint64_t grow = (static_cast<std::uint64_t>(last) - st.view_pieces.capacity()) * sizeof(st.view_pieces[0]);
                if (!st.reserve_store(grow, true)) {
                    note_ledger();
                    throw ResourceExhausted(kStoreRefusal);
                }
                st.view_table_charged += grow;
                st.view_pieces.reserve(last);
                const std::uint64_t held = st.view_pieces.capacity() * sizeof(st.view_pieces[0]);
                if (held > st.view_table_charged) {
                    if (!st.reserve_store(held - st.view_table_charged, true)) {
                        note_ledger();
                        throw ResourceExhausted(kStoreRefusal);
                    }
                    st.view_table_charged = held;
                }
            }
            st.view_pieces.resize(last, nullptr);
            for (const std::size_t u : interior) {
                const auto& pieces = ix.fine_of_coarse.at(static_cast<std::uint32_t>(u));
                st.view_pieces[u] = &pieces;
                std::uint64_t n = 0;
                for (const auto& [fine, offset] : pieces) n += st.placements[fine].size();
                ++r.view_coarse_units;
                r.view_rows_represented += n;
            }
            r.view_table_bytes = st.view_table_charged;
        } else if (st.exact_store) {
            // I2-C: each coarse unit sized from its known count and its capacity reserved before it is allocated; in a
            // one-shot batch a fine unit's rows are released once the last new coarse unit using it is built. Peak = the
            // fine rows still needed + the coarse rows built so far, never the whole store twice.
            std::unordered_map<std::uint32_t, std::uint32_t> uses;
            if (st.one_shot) {
                for (const std::size_t u : interior) {
                    for (const auto& [fine, offset] : ix.fine_of_coarse.at(static_cast<std::uint32_t>(u))) ++uses[fine];
                }
                r.worker_scratch_bytes = std::max<std::uint64_t>(r.worker_scratch_bytes, uses.bucket_count() * sizeof(void*) + uses.size() * 32);
            }
            for (const std::size_t u : interior) {
                const auto& pieces = ix.fine_of_coarse.at(static_cast<std::uint32_t>(u));
                std::uint64_t n = 0;
                for (const auto& [fine, offset] : pieces) n += st.placements[fine].size();
                if (n > (std::numeric_limits<std::uint64_t>::max() / sizeof(StoredPlacement)) - store_spare_capacity().load() ||
                    !st.reserve_store(charged_bytes(n), true)) {
                    note_ledger();
                    throw ResourceExhausted(kStoreRefusal);
                }
                std::vector<StoredPlacement> built;
                built.reserve(static_cast<std::size_t>(requested_slots(n)));
                for (const auto& [fine, offset] : pieces) {
                    for (const StoredPlacement& p : st.placements[fine]) {
                        StoredPlacement q = p;
                        q.start = p.start + offset;
                        built.push_back(q);
                    }
                }
                std::sort(built.begin(), built.end(), row_order);
                if (!store_charge_rows_only().load(std::memory_order_relaxed) && built.capacity() > requested_slots(n) &&
                    !st.reserve_store((built.capacity() - requested_slots(n)) * sizeof(StoredPlacement), true)) {
                    note_ledger();
                    throw ResourceExhausted(kStoreRefusal);
                }
                const std::uint64_t old_capacity = st.placements[u].capacity() * sizeof(StoredPlacement);
                st.placements[u].swap(built);
                if (old_capacity != 0) st.release_store(old_capacity);
                ++r.coarse_units_reconstructed;
                r.reconstructed_placements += n;
                if (st.one_shot) {
                    for (const auto& [fine, offset] : pieces) {
                        if (--uses[fine] != 0) continue;
                        const std::uint64_t bytes = st.placements[fine].capacity() * sizeof(StoredPlacement);
                        std::vector<StoredPlacement>().swap(st.placements[fine]);
                        st.release_store(bytes);
                    }
                }
            }
        } else {
            {
                // I2-A: the reconstructed rows are charged whole before any is built
                std::uint64_t rows = 0;
                for (std::size_t u = first; u < last; ++u) {
                    if (!ix.units[u].interior) continue;
                    for (const auto& [fine, offset] : ix.fine_of_coarse.at(static_cast<std::uint32_t>(u))) rows += st.placements[fine].size();
                }
                charge_store(rows * sizeof(StoredPlacement));
                r.store_charged_bytes = store_charged.load(std::memory_order_relaxed);
            }
            for (std::size_t u = first; u < last; ++u) {
                if (!ix.units[u].interior) continue;
                std::vector<StoredPlacement>& store = st.placements[u];
                store.clear();
                for (const auto& [fine, offset] : ix.fine_of_coarse.at(static_cast<std::uint32_t>(u))) {
                    for (const StoredPlacement& p : st.placements[fine]) {
                        StoredPlacement q = p;
                        q.start = p.start + offset;
                        store.push_back(q);
                    }
                }
                std::sort(store.begin(), store.end(), row_order);
                (void)st.reserve_store(store.capacity() * sizeof(StoredPlacement), false);   // measured only
                ++r.coarse_units_reconstructed;
                r.reconstructed_placements += store.size();
            }
            if (st.one_shot) {
                // I2-A: fine rows exist only to reconstruct coarse units; a one-shot batch never adds any, so release them
                for (std::size_t u = 0; u < last; ++u) {
                    if (!ix.units[u].fine) continue;
                    const std::uint64_t bytes = st.placements[u].capacity() * sizeof(StoredPlacement);
                    std::vector<StoredPlacement>().swap(st.placements[u]);
                    st.release_store(bytes);
                }
            }
        }
        r.reconstruction_peak_bytes = std::max(r.reconstruction_peak_bytes, st.store_peak.load(std::memory_order_relaxed));
        r.seconds_reconstruction += std::chrono::duration<double>(std::chrono::steady_clock::now() - reconstruct_mark).count();
    }
    if (st.one_shot) st.sealed = true;
    if (st.prefilter) fill_union_prefilter_report(*st.prefilter, r.prefilter);
    if (st.exact_store) r.store_charged_bytes = st.store_allocated.load(std::memory_order_relaxed);
    note_ledger();
    r.placement_store_bytes = 0;
    for (const auto& v : st.placements) r.placement_store_bytes += v.capacity() * sizeof(StoredPlacement);
    if (slots == nullptr) {
        st.evaluated = last;
    } else {   // S2-C
        for (const std::size_t u : direct) st.unit_done[u] = 1;
        for (const std::size_t u : interior) st.unit_done[u] = 1;
    }
    r.seconds_discovery += std::chrono::duration<double>(std::chrono::steady_clock::now() - mark).count();
}

FragmentScoringResult ContentIndexBatch::score(const std::vector<HaplotypeCandidate>& candidates) {
    if (state_->selective) throw std::logic_error("content index: bulk score() on a selective batch");
    return score_impl(candidates, false);
}
FragmentScoringResult ContentIndexBatch::score_selected(const std::vector<HaplotypeCandidate>& candidates) {
    if (!state_->selective) throw std::logic_error("content index: score_selected() before any selective discovery");
    if (state_->poisoned) throw std::logic_error("content index: a refused or failed selective batch cannot continue");
    {   // preconditions, checked before any work: a violated one changes nothing and does not poison the batch
        const ContentPlacementIndex::State& ix = *state_->index->state_;
        for (const HaplotypeCandidate& c : candidates) {
            for (const std::size_t h : {c.first, c.second}) {
                if (h >= ix.lengths.size()) throw std::invalid_argument("content index: candidate names an unknown slot");
                if (h < state_->completed.size() && state_->completed[h]) continue;
                for (int strand = 0; strand < 2; ++strand) {
                    for (const IndexOccurrence& occ : ix.occurrences[strand][h]) {
                        if (occ.unit >= state_->unit_done.size() || !state_->unit_done[occ.unit]) {
                            throw std::logic_error("content index: a requested slot needs evidence not yet discovered; an "
                                                   "incomplete column is never scored");
                        }
                    }
                }
            }
        }
    }
    try {
        return score_impl(candidates, true);
    } catch (...) {
        state_->poisoned = true;
        state_->carried.clear();
        throw;
    }
}
FragmentScoringResult ContentIndexBatch::score_impl(const std::vector<HaplotypeCandidate>& candidates, bool selective) {
    State& st = *state_;
    const ContentPlacementIndex::State& ix = *st.index->state_;
    if (!selective && st.evaluated != ix.units.size()) throw std::logic_error("content index: score() before every unit was discovered");
    const std::size_t unique = ix.lengths.size(), nf = st.fragments.size();
    for (const HaplotypeCandidate& c : candidates) {
        if (c.first >= unique || c.second >= unique) throw std::invalid_argument("content index: candidate names an unknown slot");
    }
    // S2-C: the requested slots; each needs every unit its occurrences name discovered (a coarse interior: built or
    // viewed). Unknown evidence is never zero mass: an incomplete column is refused before anything is computed.
    std::vector<char> requested(selective ? unique : 0, 0);
    std::size_t requested_new = 0;
    if (selective) {
        st.completed.resize(unique, 0);
        for (const HaplotypeCandidate& c : candidates) {
            for (const std::size_t h : {c.first, c.second}) {
                if (requested[h]) continue;
                requested[h] = 1;
                if (st.completed[h]) continue;
                ++requested_new;
                for (int strand = 0; strand < 2; ++strand) {
                    for (const IndexOccurrence& occ : ix.occurrences[strand][h]) {
                        if (occ.unit >= st.unit_done.size() || !st.unit_done[occ.unit]) {
                            throw std::logic_error("content index: a requested slot needs evidence not yet discovered; an "
                                                   "incomplete column is never scored");
                        }
                    }
                }
            }
        }
    }
    // I2-A: the reference's mass-cell cap, before any column is allocated (S2-C: the columns held after this call)
    const std::size_t columns_after = selective ? std::count(st.completed.begin(), st.completed.end(), 1) + requested_new : unique;
    if (nf != 0 && columns_after > st.options.max_mass_cells / nf) {
        throw ResourceExhausted(
            "the fragment-by-haplotype mass matrix would exceed its declared cell cap; no partial ranking is returned");
    }
    FragmentScoringResult result;
    result.unique_haplotypes = unique;
    result.scores.resize(candidates.size());
    // I2-A: executed pairing work is charged to the shared physical-origin budget as it happens (direct pairing and
    // aggregate construction alike), and served quantities to their consumption caps (0: no cap, the experimental
    // benchmark's setting; the caller's adapter requires them positive).
    PlacementBudget budget;
    budget.max_piece_hits = st.options.max_piece_hits;
    budget.max_physical_origins = st.options.max_physical_origins;
    std::atomic<std::uint64_t> placements_consumed{0}, origins_consumed{0}, states_consumed{0};
    if (!selective && st.baseline_set) {   // S2-E: a block's batch starts from the blocks before it
        budget.physical_origins.store(st.baseline.physical_origins, std::memory_order_relaxed);
        placements_consumed.store(st.baseline.placements_consumed, std::memory_order_relaxed);
        origins_consumed.store(st.baseline.origins_consumed, std::memory_order_relaxed);
        states_consumed.store(st.baseline.states_consumed, std::memory_order_relaxed);
    }
    if (selective) {   // S2-C: every cap charged cumulatively over the batch's calls, never reset between requests
        budget.physical_origins.store(st.selective_origins, std::memory_order_relaxed);
        placements_consumed.store(st.selective_placements_consumed, std::memory_order_relaxed);
        origins_consumed.store(st.selective_origins_consumed, std::memory_order_relaxed);
        states_consumed.store(st.selective_states_consumed, std::memory_order_relaxed);
    }
    FragmentScoringCounters built;   // aggregate construction: executed
    const auto note_consumption = [&]() {
        st.report.placements_consumed_last = placements_consumed.load(std::memory_order_relaxed);
        st.report.origins_consumed_last = origins_consumed.load(std::memory_order_relaxed);
        st.report.states_consumed_last = states_consumed.load(std::memory_order_relaxed);
        st.report.physical_origins_observed = budget.physical_origins.load(std::memory_order_relaxed);
    };
    // I1.3: retained columns; only slots never completed are computed (all of them with retention off: I1.2's behaviour)
    if (!st.retain) {
        std::fill(st.completed.begin(), st.completed.end(), 0);
        st.pairs.clear();
    }
    st.mass_column.resize(unique);
    st.exposure.resize(unique, 0.0);
    st.completed.resize(unique, 0);
    st.accounting.resize(unique);
    std::vector<std::size_t> todo_slots;
    for (std::size_t h = 0; h < unique; ++h) {
        if (st.completed[h]) continue;
        if (selective && !requested[h]) continue;   // S2-C: only the requested slots' columns
        todo_slots.push_back(h);
        st.mass_column[h].assign(nf, kNegInf);
        st.exposure[h] = exposure_for_length(ix.lengths[h], nf, st.prepared.fragments_by_floor, st.library_insert);
    }
    std::vector<std::vector<double>>& mass_column = st.mass_column;
    std::vector<double>& exposure = st.exposure;
    // ---- I1.2: interior-origin aggregates for every eligible entity not yet built (eligible: two or more occurrences
    // and a nonempty zone), in canonical order (most occurrences, then shorter, then entity id), admitted while they
    // fit the declared budget; a declined entity is paired directly.
    const std::uint32_t R = static_cast<std::uint32_t>(st.zone_read), g = static_cast<std::uint32_t>(ix.gap);
    ContentBatchReport& rr = st.report;
    if (st.options.union_origin_cache_bytes > 0) {
        const auto build_mark = std::chrono::steady_clock::now();
        st.entity_entry.resize(ix.entities.size(), -1);
        std::vector<std::uint32_t> todo;
        // S2-C: only the entities the slots computed now contain (each built at most once per batch)
        std::vector<char> relevant;
        if (selective) {
            relevant.assign(ix.entities.size(), 0);
            for (const std::size_t h : todo_slots) for (const auto& [entity, offset] : ix.entity_occurrences[h]) relevant[entity] = 1;
        }
        for (std::uint32_t e = 0; e < ix.entities.size(); ++e) {
            if (selective && !relevant[e]) continue;
            if (st.entity_entry[e] >= 0 || st.entity_entry[e] == -2) continue;
            const auto& en = ix.entities[e];
            const bool eligible = en.occurrences >= 2 && en.length >= R + 2 * g + 1;
            st.entity_entry[e] = eligible ? -1 : -3;
            if (eligible) todo.push_back(e);
        }
        std::sort(todo.begin(), todo.end(), [&](std::uint32_t a, std::uint32_t b) {
            const auto& x = ix.entities[a];
            const auto& y = ix.entities[b];
            return std::make_tuple(-static_cast<long>(x.occurrences), x.length, a) < std::make_tuple(-static_cast<long>(y.occurrences), y.length, b);
        });
        rr.origin_entities_eligible += todo.size();
        const std::size_t workers = std::max<std::size_t>(1, std::min(st.options.threads, todo.size()));
        std::vector<StateAggregateSet> scratch(workers);
        std::vector<double> gather_ns(workers, 0.0);   // I4-B: seconds gathering view zones, per worker slot
        for (std::size_t begin = 0; begin < todo.size(); begin += workers) {
            const std::size_t end = std::min(begin + workers, todo.size());
            std::vector<UnionOriginEntry> made(end - begin);
            std::vector<std::exception_ptr> failures(end - begin);
            std::vector<std::pair<std::uint64_t, std::uint64_t>> zone_rows(end - begin, {0, 0});   // I4-B: per-entity counts
            const auto build_one = [&](std::size_t slot) {
                const ActiveWorker active;
                const auto& en = ix.entities[todo[begin + slot]];
                // the zone's placements only: starts in [g, length - R - g] on each strand
                const auto in_zone = [&](const std::vector<StoredPlacement>& all) {
                    std::vector<CachedPlacement> kept;
                    for (const StoredPlacement& p : all) if (p.start >= g && p.start + R + g <= en.length) kept.push_back(p);
                    return kept;
                };
                if (!st.fine_views) {
                    build_origin_entry(in_zone(st.placements[en.forward_unit]), in_zone(st.placements[en.reverse_unit]),
                                       en.length, st.prepared.forward_forms.size(), nf, st.prepared, &budget, scratch[slot],
                                       made[slot]);
                    return;
                }
                // I4-B: the zone rows gathered from the unit's view (each fine piece's rows at piece offset + start), in
                // reconstruction's (form, start, consumed) order -- keys are unique within a coarse unit (disjoint owned
                // cores). These gathers are the views' transition-time copies; the whole batch's were counted and charged
                // to the store ledger before any worker started (zone_rows), and are released after the batch joins.
                const auto gather_mark = std::chrono::steady_clock::now();
                const auto view_zone = [&](std::uint32_t unit, std::uint64_t n) {
                    const auto* pieces = unit < st.view_pieces.size() ? st.view_pieces[unit] : nullptr;
                    if (pieces == nullptr) return in_zone(st.placements[unit]);
                    std::vector<CachedPlacement> kept;
                    kept.reserve(static_cast<std::size_t>(n));
                    for (const auto& [fine, offset] : *pieces) {
                        for (const StoredPlacement& p : st.placements[fine]) {
                            const std::uint32_t at = p.start + offset;
                            if (at < g || at + R + g > en.length) continue;
                            CachedPlacement q = p;
                            q.start = at;
                            kept.push_back(q);
                        }
                    }
                    if (kept.size() != n || kept.capacity() != n) {
                        throw std::logic_error("content index: a view zone gather differs from its charged count");
                    }
                    std::sort(kept.begin(), kept.end(), [](const CachedPlacement& x, const CachedPlacement& y) {
                        return std::tie(x.form, x.start, x.consumed) < std::tie(y.form, y.start, y.consumed);
                    });
                    return kept;
                };
                const std::vector<CachedPlacement> forward = view_zone(en.forward_unit, zone_rows[slot].first);
                const std::vector<CachedPlacement> reverse = view_zone(en.reverse_unit, zone_rows[slot].second);
                gather_ns[slot] += std::chrono::duration<double>(std::chrono::steady_clock::now() - gather_mark).count();
                build_origin_entry(forward, reverse, en.length, st.prepared.forward_forms.size(), nf, st.prepared, &budget,
                                   scratch[slot], made[slot]);
            };
            // I4-B: under views, each entity's zone rows are counted first (in parallel, no allocation) and the batch's
            // total is charged to the store ledger before any gather is allocated: the ledger's peak, and so the store
            // verdict, depend on the thread count only, never on the schedule. Released after the batch joins.
            std::uint64_t batch_zone_bytes = 0;
            if (st.fine_views) {
                std::vector<std::thread> counters_pool;
                for (std::size_t slot = 0; slot < end - begin; ++slot) {
                    counters_pool.emplace_back([&, slot]() {
                        const auto& en = ix.entities[todo[begin + slot]];
                        const auto count = [&](std::uint32_t unit) -> std::uint64_t {
                            const auto* pieces = unit < st.view_pieces.size() ? st.view_pieces[unit] : nullptr;
                            if (pieces == nullptr) return 0;   // not a view: the old uncharged path applies
                            std::uint64_t n = 0;
                            for (const auto& [fine, offset] : *pieces) {
                                for (const StoredPlacement& p : st.placements[fine]) {
                                    const std::uint32_t at = p.start + offset;
                                    n += (at >= g && at + R + g <= en.length) ? 1 : 0;
                                }
                            }
                            return n;
                        };
                        zone_rows[slot] = {count(en.forward_unit), count(en.reverse_unit)};
                    });
                }
                for (std::thread& w : counters_pool) w.join();
                for (std::size_t slot = 0; slot < end - begin; ++slot) {
                    batch_zone_bytes += (zone_rows[slot].first + zone_rows[slot].second) * sizeof(CachedPlacement);
                }
                if (batch_zone_bytes != 0 && !st.reserve_store(batch_zone_bytes, true)) {
                    note_consumption();
                    throw ResourceExhausted("the indexed union's placement store reached its declared byte budget; no partial "
                                            "ranking is returned");
                }
                const std::uint64_t now = st.zone_live.fetch_add(batch_zone_bytes, std::memory_order_relaxed) + batch_zone_bytes;
                st.zone_peak.store(std::max(st.zone_peak.load(std::memory_order_relaxed), now), std::memory_order_relaxed);
                st.zone_total.fetch_add(batch_zone_bytes, std::memory_order_relaxed);
            }
            struct ZoneRelease {
                State& s;
                std::uint64_t bytes;
                ~ZoneRelease() {
                    if (bytes == 0) return;
                    s.release_store(bytes);
                    s.zone_live.fetch_sub(bytes, std::memory_order_relaxed);
                }
            } zone_release{st, batch_zone_bytes};
            std::vector<std::thread> pool;
            for (std::size_t slot = 0; slot < end - begin; ++slot) {
                pool.emplace_back([&, slot]() {
                    try { build_one(slot); } catch (...) { failures[slot] = std::current_exception(); budget.cancel(); }
                });
            }
            for (std::thread& w : pool) w.join();
            if (std::any_of(failures.begin(), failures.end(), [](const std::exception_ptr& f) { return static_cast<bool>(f); })) {
                note_consumption();
                rethrow_primary_failure(failures);
            }
            for (std::size_t slot = 0; slot < end - begin; ++slot) {
                UnionOriginEntry& entry = made[slot];
                built.physical_origins += entry.origins;   // I2-A: construction is executed work, as in the reference
                built.statistical_state_classes += entry.states.states.size();
                built.latent_alignment_pairs += entry.latent_pairs;
                built.valid_fr_states += entry.valid_states;
                const std::uint64_t bytes = sizeof(UnionOriginEntry) +
                    static_cast<std::uint64_t>(entry.states.fragments.size() + entry.states.starts.size()) * sizeof(std::uint32_t) +
                    static_cast<std::uint64_t>(entry.states.states.size()) * sizeof(StateCount);
                rr.origins_constructed += entry.origins;
                ++rr.origin_entities_built;
                if (st.origin_held + bytes <= st.options.union_origin_cache_bytes) {
                    st.origin_held += bytes;
                    st.entity_entry[todo[begin + slot]] = static_cast<long>(st.origin_entries.size());
                    st.origin_entries.push_back(std::move(entry));
                    ++rr.origin_entities_admitted;
                } else {
                    st.entity_entry[todo[begin + slot]] = -2;
                    ++rr.origin_entities_declined;
                }
            }
        }
        rr.origin_held_bytes = st.origin_held;
        for (const double x : gather_ns) rr.seconds_view_zone_gather += x;
        rr.view_zone_temporary_peak_bytes = st.zone_peak.load(std::memory_order_relaxed);
        rr.view_zone_temporary_bytes_total = st.zone_total.load(std::memory_order_relaxed);
        rr.store_peak_bytes = st.store_peak.load(std::memory_order_relaxed);
        rr.store_live_bytes = st.store_live.load(std::memory_order_relaxed);
        rr.store_allocated_bytes = st.store_allocated.load(std::memory_order_relaxed);
        rr.store_released_bytes = st.store_released.load(std::memory_order_relaxed);
        rr.seconds_origin_construction += std::chrono::duration<double>(std::chrono::steady_clock::now() - build_mark).count();
    }
    const bool zoned = st.options.union_origin_cache_bytes > 0;

    const auto mark = std::chrono::steady_clock::now();
    const std::size_t threads = std::max<std::size_t>(1, std::min(st.options.threads, std::max<std::size_t>(1, todo_slots.size())));
    std::atomic<std::size_t> next{0};
    std::vector<double> expansion_ns(threads, 0), pairing_ns(threads, 0);
    std::vector<std::uint64_t> served(threads, 0), direct(threads, 0), zones_used(threads, 0);
    std::vector<std::uint64_t> expanded(threads, 0), removed(threads, 0), origins_count(threads, 0), classes(threads, 0);
    std::vector<FragmentScoringCounters> executed(threads);   // I2-A: the reference's counters, executed work only
    std::vector<std::exception_ptr> failures(threads);
    const auto charge_origin = [&]() {
        const std::uint64_t seen = budget.physical_origins.fetch_add(1, std::memory_order_relaxed) + 1;
        if (seen > budget.max_physical_origins) {
            budget.cancel();
            throw ResourceExhausted("fragment scoring reached its physical-origin cap; no partial ranking is returned");
        }
        if (budget.stop()) throw PlacementCancelled();
    };
    const auto consume = [&](std::atomic<std::uint64_t>& total, std::uint64_t amount, std::uint64_t cap, const char* message) {
        const std::uint64_t now = total.fetch_add(amount, std::memory_order_relaxed) + amount;
        if (cap != 0 && now > cap) {
            budget.cancel();
            throw ResourceExhausted(message);
        }
    };
    const auto work = [&](std::size_t t) {
        const ActiveWorker active;
        std::vector<std::vector<FormPlacement>> rows[2] = {
            std::vector<std::vector<FormPlacement>>(st.prepared.forward_forms.size()),
            std::vector<std::vector<FormPlacement>>(st.prepared.reverse_forms.size())};
        std::vector<std::uint32_t> touched[2];
        std::map<StateKey, std::uint64_t> states;
        std::map<std::uint32_t, OriginAccumulator> origins;
        std::vector<OriginZone> zones;
        std::vector<std::size_t> zone_cursor;
        std::vector<StateKey> direct_keys;
        std::vector<StateCount> direct_counts, merged_counts, merge_scratch;
        try {
            for (std::size_t i = next++; i < todo_slots.size(); i = next++) {
                if (budget.stop()) throw PlacementCancelled();
                const std::size_t h = todo_slots[i];
                const std::uint64_t b_expanded = expanded[t], b_removed = removed[t], b_plain = origins_count[t];
                const std::uint64_t b_direct = direct[t], b_served = served[t], b_classes = classes[t], b_zones = zones_used[t];
                const auto e0 = std::chrono::steady_clock::now();
                const std::uint32_t n = ix.lengths[h];
                for (int strand = 0; strand < 2; ++strand) {
                    const auto expand_rows = [&](const std::vector<StoredPlacement>& unit_rows, std::uint32_t offset) {
                        for (const StoredPlacement& p : unit_rows) {
                            auto& row = rows[strand][p.form];
                            if (row.empty()) touched[strand].push_back(p.form);
                            FormPlacement q;
                            q.start = offset + p.start;
                            q.consumed = p.consumed;
                            q.mismatches = p.mismatches;
                            q.emission = p.emission;
                            row.push_back(q);
                            ++expanded[t];
                        }
                    };
                    for (const IndexOccurrence& occ : ix.occurrences[strand][h]) {
                        // I4-B: a coarse unit served as a view expands each fine piece's rows at its piece offset: the
                        // same multiset of placements as its reconstructed rows (the per-form sort below orders them)
                        const auto* pieces = occ.unit < st.view_pieces.size() ? st.view_pieces[occ.unit] : nullptr;
                        if (pieces != nullptr) {
                            for (const auto& [fine, offset] : *pieces) expand_rows(st.placements[fine], occ.offset + offset);
                        } else {
                            expand_rows(st.placements[occ.unit], occ.offset);
                        }
                    }
                    for (const std::uint32_t form : touched[strand]) {
                        auto& row = rows[strand][form];
                        if (strand == 1) for (FormPlacement& q : row) q.start = n - q.start;   // reverse rows: forward end
                        std::sort(row.begin(), row.end(), [](const FormPlacement& a, const FormPlacement& b) {
                            return std::tie(a.start, a.consumed) < std::tie(b.start, b.consumed);
                        });
                        const auto end = std::unique(row.begin(), row.end(), [](const FormPlacement& a, const FormPlacement& b) {
                            return a.start == b.start && a.consumed == b.consumed;
                        });
                        removed[t] += static_cast<std::uint64_t>(row.end() - end);
                        row.erase(end, row.end());
                    }
                }
                const auto e1 = std::chrono::steady_clock::now();
                // I2-A: every expanded placement is served from the store, not executed: charged whole for the slot
                consume(placements_consumed, expanded[t] - b_expanded, st.options.max_union_cached_placements_consumed,
                        "fragment scoring reached its cached-placement consumption cap; no partial ranking is returned");
                // I1.2: the slot's served zones, one per admitted entity occurrence (ascending, disjoint): a pair is
                // served exactly when its forward start and reverse end lie in the SAME occurrence's interior zone
                zones.clear();
                if (zoned) {
                    for (const auto& [entity, offset] : ix.entity_occurrences[h]) {
                        const long entry = st.entity_entry[entity];
                        if (entry < 0) continue;
                        const std::uint32_t len = ix.entities[entity].length;
                        OriginZone z;
                        z.start_lo = offset + g;
                        z.start_hi = offset + len - (R + g);
                        z.end_lo = offset + R + g;
                        z.end_hi = offset + len - g;
                        z.entry = &st.origin_entries[static_cast<std::size_t>(entry)];
                        zones.push_back(z);
                        served[t] += z.entry->origins;
                    }
                    zone_cursor.assign(zones.size(), 0);
                    zones_used[t] += zones.size();
                    // I2-A: the occurrences served, charged whole before any is merged (the reference's phase-2 rule)
                    std::uint64_t served_origins = 0, served_states = 0;
                    for (const OriginZone& z : zones) {
                        served_origins += z.entry->origins;
                        served_states += z.entry->states.states.size();
                        executed[t].reused_latent_alignment_pairs += z.entry->latent_pairs;
                        executed[t].reused_valid_fr_states += z.entry->valid_states;
                    }
                    executed[t].reused_physical_origins += served_origins;
                    executed[t].reused_state_aggregates += served_states;
                    consume(origins_consumed, served_origins, st.options.max_union_cached_origins_consumed,
                            "fragment scoring reached its cached-origin consumption cap; no partial ranking is returned");
                    consume(states_consumed, served_states, st.options.max_union_cached_states_consumed,
                            "fragment scoring reached its cached-state consumption cap; no partial ranking is returned");
                }
                std::uint64_t latent_pairs = 0, valid_states = 0;
                for (std::size_t f = 0; f < nf; ++f) {
                    const InsertSupport& support = st.prepared.support_by_fragment[f];
                    if (support.empty()) continue;
                    if (!zones.empty()) {
                        const RowView forward_rows[2] = {{rows[0][2 * f].data(), rows[0][2 * f].size()},
                                                         {rows[0][2 * f + 1].data(), rows[0][2 * f + 1].size()}};
                        const RowView reverse_rows[2] = {{rows[1][2 * f].data(), rows[1][2 * f].size()},
                                                         {rows[1][2 * f + 1].data(), rows[1][2 * f + 1].size()}};
                        direct_keys.clear();
                        pair_fragment(forward_rows, reverse_rows, support, zones.data(), zones.size(), origins, latent_pairs,
                                      valid_states,
                                      [&](std::size_t arm, std::uint32_t, std::uint32_t insert, const OriginAccumulator& origin) {
                                          ++direct[t];
                                          charge_origin();
                                          StateKey key;
                                          key.insert_bp = insert;
                                          key.orientation = static_cast<std::uint32_t>(arm);
                                          key.log_latent_mass = origin.log_latent_mass;
                                          direct_keys.push_back(key);
                                      });
                        direct_counts.clear();
                        sort_and_count(direct_keys, direct_counts, 0);
                        executed[t].statistical_state_classes += direct_counts.size();
                        merged_counts.assign(direct_counts.begin(), direct_counts.end());
                        merge_served_states(f, zones.data(), zones.size(), zone_cursor, merged_counts, merge_scratch);
                        if (merged_counts.empty()) continue;
                        executed[t].merged_state_classes += merged_counts.size();
                        classes[t] += merged_counts.size();
                        mass_column[h][f] = reduce_sorted(merged_counts.data(), merged_counts.data() + merged_counts.size(),
                                                          st.library_insert);
                        continue;
                    }
                    const RowView forward_rows[2] = {{rows[0][2 * f].data(), rows[0][2 * f].size()},
                                                     {rows[0][2 * f + 1].data(), rows[0][2 * f + 1].size()}};
                    const RowView reverse_rows[2] = {{rows[1][2 * f].data(), rows[1][2 * f].size()},
                                                     {rows[1][2 * f + 1].data(), rows[1][2 * f + 1].size()}};
                    states.clear();
                    pair_fragment(forward_rows, reverse_rows, support, nullptr, 0, origins, latent_pairs, valid_states,
                                  [&](std::size_t arm, std::uint32_t, std::uint32_t insert, const OriginAccumulator& origin) {
                                      ++origins_count[t];
                                      charge_origin();
                                      StateKey key;
                                      key.insert_bp = insert;
                                      key.orientation = static_cast<std::uint32_t>(arm);
                                      key.log_latent_mass = origin.log_latent_mass;
                                      ++states[key];
                                  });
                    if (states.empty()) continue;
                    classes[t] += states.size();
                    executed[t].statistical_state_classes += states.size();
                    mass_column[h][f] = reduce_states(states, st.library_insert);
                }
                executed[t].latent_alignment_pairs += latent_pairs;
                executed[t].valid_fr_states += valid_states;
                for (int strand = 0; strand < 2; ++strand) {
                    for (const std::uint32_t form : touched[strand]) rows[strand][form].clear();
                    touched[strand].clear();
                }
                const auto e2 = std::chrono::steady_clock::now();
                expansion_ns[t] += std::chrono::duration<double>(e1 - e0).count();
                pairing_ns[t] += std::chrono::duration<double>(e2 - e1).count();
                st.accounting[h] = {expanded[t] - b_expanded, removed[t] - b_removed, origins_count[t] - b_plain,
                                    direct[t] - b_direct, served[t] - b_served, classes[t] - b_classes, zones_used[t] - b_zones};
            }
        } catch (...) {
            failures[t] = std::current_exception();
            budget.cancel();
        }
    };
    std::vector<std::thread> workers;
    for (std::size_t t = 0; t < threads; ++t) workers.emplace_back(work, t);
    for (std::thread& w : workers) w.join();
    note_consumption();
    rr.threads_used_pairing = threads;
    {
        // I2-A: this call's executed work (discovery's seeding counters are cumulative), merged before any refusal
        // propagates; then the reference's acceptance caps on the merged totals
        FragmentScoringCounters c = st.discovery_counters;
        c.fragments_with_floor_binding = st.prepared.fragments_with_floor_binding;
        FragmentScoringCounters call;   // this call's executed pairing and aggregate work
        add_counters(call, built);
        for (const FragmentScoringCounters& e : executed) add_counters(call, e);
        for (std::size_t t = 0; t < threads; ++t) call.physical_origins += origins_count[t] + direct[t];
        if (selective) {   // S2-C: the acceptance caps hold to the batch's cumulative executed totals
            add_counters(st.selective_score_counters, call);
            add_counters(c, st.selective_score_counters);
            st.selective_origins = budget.physical_origins.load(std::memory_order_relaxed);
            st.selective_placements_consumed = placements_consumed.load(std::memory_order_relaxed);
            st.selective_origins_consumed = origins_consumed.load(std::memory_order_relaxed);
            st.selective_states_consumed = states_consumed.load(std::memory_order_relaxed);
        } else {
            add_counters(c, call);
        }
        rr.counters_last = c;
        rethrow_primary_failure(failures);
        if (!selective && st.baseline_set) {   // S2-E: the acceptance caps hold to the whole union's executed totals
            FragmentScoringCounters whole = st.baseline.executed;
            add_counters(whole, c);
            whole.piece_hits = st.baseline.executed.piece_hits + c.piece_hits;
            enforce_acceptance_caps(whole, st.options, [] {});
        } else {
            enforce_acceptance_caps(c, st.options, [] {});
        }
    }
    ContentBatchReport& r = st.report;
    r.seconds_expansion_and_pairing = std::chrono::duration<double>(std::chrono::steady_clock::now() - mark).count();
    r.worker_seconds_expansion = r.worker_seconds_pairing = 0;
    for (std::size_t t = 0; t < threads; ++t) {
        r.worker_seconds_expansion += expansion_ns[t];
        r.worker_seconds_pairing += pairing_ns[t];
    }
    for (const std::size_t h : todo_slots) st.completed[h] = 1;
    r.columns_computed_last = todo_slots.size();
    r.columns_computed_total += todo_slots.size();
    r.columns_retained_last = selective ? static_cast<std::uint64_t>(std::count(requested.begin(), requested.end(), 1)) - todo_slots.size()
                                        : unique - todo_slots.size();
    // totals over every completed slot (each counted once, whenever it was computed)
    r.expanded_placements = r.duplicate_placements_removed = r.physical_origins = r.state_classes = 0;
    r.direct_physical_origins = r.served_physical_origins = r.zones_applied = 0;
    for (std::size_t h = 0; h < unique; ++h) {
        const auto& a = st.accounting[h];
        r.expanded_placements += a.expanded;
        r.duplicate_placements_removed += a.removed;
        r.direct_physical_origins += a.direct;
        r.served_physical_origins += a.served;
        r.zones_applied += a.zones;
        r.physical_origins += a.plain_origins + a.direct + a.served;   // logical: executed plus served
        r.state_classes += a.classes;
    }
    result.seconds_placement_and_mass = r.seconds_discovery + r.seconds_expansion_and_pairing;
    // pair likelihoods: each unordered slot pair combined once (shared combination), then served from the cache
    const auto combination_mark = std::chrono::steady_clock::now();
    const auto key_of = [](std::size_t a, std::size_t b) {
        return (static_cast<std::uint64_t>(std::min(a, b)) << 32) | static_cast<std::uint64_t>(std::max(a, b));
    };
    std::vector<HaplotypeCandidate> fresh;
    std::vector<std::uint64_t> fresh_keys;
    {
        std::unordered_set<std::uint64_t> seen;
        for (const HaplotypeCandidate& c : candidates) {
            const std::uint64_t key = key_of(c.first, c.second);
            if (st.pairs.count(key) != 0 || !seen.insert(key).second) continue;
            fresh.push_back({std::min(c.first, c.second), std::max(c.first, c.second)});
            fresh_keys.push_back(key);
        }
    }
    r.pairs_requested_last = candidates.size();
    r.pairs_computed_last = fresh.size();
    r.pairs_computed_total += fresh.size();
    r.seconds_candidate_scoring = 0.0;
    if (!fresh.empty()) {
        FragmentScoringResult computed;
        computed.scores.resize(fresh.size());
        combine_candidate_scores(mass_column, exposure, st.prepared.background_by_fragment, nf, st.library, fresh,
                                 st.options.threads, st.options.audit_libraries, computed, st.options.prepared_library_terms);
        r.seconds_candidate_scoring = computed.seconds_diploid_combination;
        result.threads_used_combination = computed.threads_used_combination;
        result.candidate_fragment_evaluations = computed.candidate_fragment_evaluations;
        for (std::size_t i = 0; i < fresh.size(); ++i) {
            State::PairResult pr;
            pr.score = computed.scores[i];
            for (const auto& audit : computed.audit_scores) pr.audits.push_back(audit[i]);
            st.pairs.emplace(fresh_keys[i], std::move(pr));
        }
    }
    const std::size_t audit_count = audit_combinations(st.options.audit_libraries).size();
    result.audit_scores.assign(audit_count, std::vector<DiplotypeScore>(candidates.size()));
    for (std::size_t c = 0; c < candidates.size(); ++c) {
        const State::PairResult& pr = st.pairs.at(key_of(candidates[c].first, candidates[c].second));
        result.scores[c] = pr.score;
        for (std::size_t k = 0; k < audit_count; ++k) result.audit_scores[k][c] = pr.audits[k];
    }
    result.seconds_diploid_combination = std::chrono::duration<double>(std::chrono::steady_clock::now() - combination_mark).count();
    r.seconds_pair_lookup_and_combination = result.seconds_diploid_combination;
    result.counters = r.counters_last;
    result.threads_used_placement = std::max<std::size_t>(1, std::max(r.threads_used_discovery, threads));
    // output materialisation: the flattened mass matrix and exposures of every slot (S2-C: none on the selective route --
    // the scores are returned, the columns stay retained, and nothing is exported)
    const auto materialise_mark = std::chrono::steady_clock::now();
    if (!selective && st.materialise) {
        result.log_mass.assign(nf * unique, kNegInf);
        for (std::size_t h = 0; h < unique; ++h)
            for (std::size_t f = 0; f < nf; ++f) result.log_mass[f * unique + h] = mass_column[h][f];
        result.exposure = exposure;
    }
    r.seconds_materialisation = std::chrono::duration<double>(std::chrono::steady_clock::now() - materialise_mark).count();
    r.retained_column_bytes = 0;
    for (const auto& column : mass_column) r.retained_column_bytes += column.capacity() * sizeof(double);
    r.retained_column_bytes += exposure.capacity() * sizeof(double) + st.accounting.capacity() * sizeof(State::SlotAccounting);
    r.retained_pairs = st.pairs.size();
    r.retained_pair_bytes = st.pairs.size() * (sizeof(std::uint64_t) + sizeof(State::PairResult) + 32 + audit_count * sizeof(DiplotypeScore)) +
                            st.pairs.bucket_count() * sizeof(void*);
    return result;
}

std::map<std::string, std::vector<std::tuple<std::uint32_t, std::uint32_t, std::uint32_t, std::uint32_t, double>>>
ContentIndexBatch::interior_placement_rows() const {
    const State& st = *state_;
    const ContentPlacementIndex::State& ix = *st.index->state_;
    std::map<std::string, std::vector<std::tuple<std::uint32_t, std::uint32_t, std::uint32_t, std::uint32_t, double>>> out;
    for (std::size_t u = 0; u < st.evaluated; ++u) {
        if (!ix.units[u].interior) continue;
        auto& rows = out[std::string(1, static_cast<char>('0' + ix.units[u].strand)) + ix.units[u].bytes];
        const auto* pieces = u < st.view_pieces.size() ? st.view_pieces[u] : nullptr;
        if (pieces != nullptr) {   // I4-B: a view, materialised on demand in reconstruction's order
            std::vector<StoredPlacement> built;
            for (const auto& [fine, offset] : *pieces) {
                for (StoredPlacement p : st.placements[fine]) {
                    p.start += offset;
                    built.push_back(p);
                }
            }
            std::sort(built.begin(), built.end(), [](const StoredPlacement& a, const StoredPlacement& b) {
                return std::tie(a.form, a.start, a.consumed) < std::tie(b.form, b.start, b.consumed);
            });
            for (const StoredPlacement& p : built) rows.emplace_back(p.form, p.start, p.consumed, p.mismatches, p.emission);
            continue;
        }
        for (const StoredPlacement& p : st.placements[u]) rows.emplace_back(p.form, p.start, p.consumed, p.mismatches, p.emission);
    }
    return out;
}

// ---- I2-A (EXPERIMENTAL): the indexed route at the caller's one-shot union boundary (see the header) -------------------
namespace {
UnionOriginCacheReport indexed_origin_report(const ContentBatchReport& b, const FragmentScoringOptions& o) {
    UnionOriginCacheReport u;
    u.active = o.union_origin_cache_bytes > 0;
    u.budget_bytes = o.union_origin_cache_bytes;
    u.max_origins_consumed = o.max_union_cached_origins_consumed;
    u.max_states_consumed = o.max_union_cached_states_consumed;
    u.eligible_segments = b.origin_entities_eligible;   // segment entities with two or more occurrences and a zone
    u.admitted_segments = b.origin_entities_admitted;
    u.declined_segments = b.origin_entities_declined;
    u.hits = b.zones_applied;                            // occurrences served
    u.origins_constructed = b.origins_constructed;
    u.origins_consumed = b.origins_consumed_last;
    u.states_consumed = b.states_consumed_last;
    u.payload_bytes = b.origin_held_bytes;
    u.peak_bytes = b.origin_held_bytes;
    u.seconds_fill = b.seconds_origin_construction;
    return u;
}
}  // namespace

FragmentScoringResult score_diplotypes_indexed_union(const InternedHaplotypes& haplotypes,
                                                     const std::vector<HaplotypeCandidate>& candidates,
                                                     const std::vector<ReadPair>& fragments, const LibraryModel& library,
                                                     const FragmentScoringOptions& options,
                                                     const std::vector<std::vector<std::uint32_t>>& segments,
                                                     std::size_t store_budget_bytes, IndexedUnionReport* report,
                                                     FragmentScoringObservation* observed_on_failure,
                                                     bool legacy_store_layout, bool fine_views) {
    using Clock = std::chrono::steady_clock;
    const auto since = [](Clock::time_point mark) { return std::chrono::duration<double>(Clock::now() - mark).count(); };
    // score_diplotypes' own validation, then the adapter's declared requirements
    validate_library_model(library);
    if (options.minimum_piece_bp == 0) throw std::invalid_argument("minimum indexable piece length must be positive");
    if (options.max_piece_hits == 0 || options.max_proposed_starts_per_haplotype == 0 || options.max_physical_origins == 0 ||
        options.max_state_classes == 0 || options.max_mass_cells == 0 || options.max_latent_alignment_pairs == 0) {
        throw std::invalid_argument("every fragment-scoring resource cap must be positive");
    }
    const std::size_t unique = haplotypes.sequences.size();
    for (const HaplotypeCandidate& c : candidates) {
        if (c.first >= unique || c.second >= unique) throw std::invalid_argument("diplotype candidate names an unknown haplotype");
    }
    if (options.record_origin_summary || !options.record_origins_for_fragments.empty() ||
        options.record_haplotype_work != nullptr) {
        throw std::invalid_argument("the indexed union backend records no per-origin or per-haplotype diagnostics");
    }
    if (segments.size() != unique) throw std::invalid_argument("the indexed union backend needs one segmentation per unique haplotype");
    if (store_budget_bytes == 0) throw std::invalid_argument("the indexed union backend needs a positive placement-store budget");
    if (options.max_union_cached_placements_consumed == 0) {
        throw std::invalid_argument("the indexed union backend needs a positive cached-placement consumption cap: every "
                                    "expanded placement is served from its store");
    }
    if (options.union_origin_cache_bytes > 0 &&
        (options.max_union_cached_origins_consumed == 0 || options.max_union_cached_states_consumed == 0)) {
        throw std::invalid_argument("union origin reuse needs both of its consumption caps, positive");
    }
    IndexedUnionReport local;
    IndexedUnionReport& rep = report != nullptr ? *report : local;
    rep = IndexedUnionReport{};
    rep.store_budget_bytes = store_budget_bytes;
    if (!fragments.empty() && unique > options.max_mass_cells / fragments.size()) {
        // before any work, as in the reference (whose observation is likewise left at its defaults here)
        rep.refusal = "the fragment-by-haplotype mass matrix would exceed its declared cell cap; no partial ranking is returned";
        throw ResourceExhausted(rep.refusal);
    }
    rep.max_read_bp = 1;
    for (const ReadPair& f : fragments) rep.max_read_bp = std::max({rep.max_read_bp, f.first.bases.size(), f.second.bases.size()});
    FragmentScoringOptions batch_options = options;
    batch_options.union_segment_lengths = nullptr;     // the index carries the segmentation; no phase-1 replay
    batch_options.union_placement_cache_bytes = 0;
    {
        StartSortTally& tally = start_sort_tally();
        tally.out_of_place.store(0);
        tally.in_place.store(0);
        tally.fallbacks.store(0);
        tally.peak_worker_bytes.store(0);
    }
    const auto start_sort_report = [&]() {
        StartSortReport q;
        const StartSortTally& tally = start_sort_tally();
        q.budget_bytes = options.start_sort_scratch_bytes;
        q.workers = std::max<std::size_t>(1, options.threads);
        q.share_bytes = options.start_sort_scratch_bytes / q.workers;
        q.out_of_place = tally.out_of_place.load();
        q.in_place = tally.in_place.load();
        q.fallbacks = tally.fallbacks.load();
        q.peak_worker_bytes = tally.peak_worker_bytes.load();
        return q;
    };
    std::unique_ptr<ContentPlacementIndex> index;
    std::unique_ptr<ContentIndexBatch> batch;
    const auto release = [&]() {
        const auto mark = Clock::now();
        batch.reset();
        index.reset();
        rep.seconds_release = since(mark);
    };
    const auto started = Clock::now();
    FragmentScoringResult result;
    try {
        auto mark = Clock::now();
        index = std::make_unique<ContentPlacementIndex>(rep.max_read_bp, library.max_gap_bp, true);
        for (std::size_t h = 0; h < unique; ++h) {
            if (index->add_haplotype(haplotypes.sequences[h], segments[h]) != h) {
                throw std::logic_error("the indexed union backend's slots must follow the interned order");
            }
        }
        rep.seconds_index_build = since(mark);
        mark = Clock::now();
        batch = std::make_unique<ContentIndexBatch>(*index, fragments, library, batch_options);
        rep.seconds_preparation = since(mark);
        // the reference refuses mixed read lengths under the prefilter cache; so does this route (no silent fallback)
        if (options.union_prefilter_cache_bytes > 0 && batch->report().prefilter_status != "on") {
            throw std::invalid_argument("union prefilter reuse needs one read length");
        }
        batch->set_placement_store_budget(store_budget_bytes);
        batch->set_legacy_store_layout(legacy_store_layout);
        batch->set_fine_views(fine_views);   // I4-B (EXPERIMENTAL)
        batch->set_one_shot(true);
        batch->discover();
        result = batch->score(candidates);
    } catch (const ResourceExhausted& error) {
        rep.refusal = error.what();
        if (index) rep.index = index->report();
        if (batch) rep.batch = batch->report();
        if (observed_on_failure != nullptr) {
            FragmentScoringObservation& o = *observed_on_failure;
            o = FragmentScoringObservation{};
            o.counters = rep.batch.counters_last;
            o.counters.physical_origins = std::max(o.counters.physical_origins, rep.batch.physical_origins_observed);
            o.threads_used_placement =
                std::max<std::size_t>(1, std::max(rep.batch.threads_used_discovery, rep.batch.threads_used_pairing));
            o.threads_used_combination = 0;
            o.seconds_fragment_preparation = rep.seconds_preparation;
            o.seconds_placement_and_mass = std::max(0.0, since(started) - rep.seconds_preparation);
            o.seconds_diploid_combination = 0.0;
            o.union_origin_cache = indexed_origin_report(rep.batch, options);
            o.union_prefilter_cache = rep.batch.prefilter;
            o.start_sort = start_sort_report();
        }
        release();
        throw;
    } catch (...) {
        release();
        throw;
    }
    rep.index = index->report();
    rep.batch = batch->report();
    const ContentBatchReport& b = rep.batch;
    result.union_origin_cache = indexed_origin_report(b, options);
    result.union_prefilter_cache = b.prefilter;
    result.start_sort = start_sort_report();
    result.emission_cache = b.emission_cache;
    result.seconds_fragment_preparation = rep.seconds_preparation;
    result.seconds_placement_and_mass =
        rep.seconds_index_build + b.seconds_discovery + b.seconds_origin_construction + b.seconds_expansion_and_pairing;
    // memory accounts, counted from the containers (placement peaks are not tracked on this route and stay 0)
    result.memory.mass_matrix_bytes = static_cast<std::uint64_t>(result.log_mass.size()) * sizeof(double);
    result.memory.haplotype_sequence_bytes = rep.index.unit_bytes;
    for (const auto* forms : {&batch->state_->prepared.forward_forms, &batch->state_->prepared.reverse_forms}) {
        for (const ReadForm& form : *forms) result.memory.form_sequence_bytes += form.bases.capacity();
    }
    for (const auto* indexes : {&batch->state_->prepared.forward_indexes, &batch->state_->prepared.reverse_indexes}) {
        for (const PieceIndex& pi : *indexes) result.memory.piece_index_bytes += pi.entries.capacity() * sizeof(PieceEntry);
    }
    last_memory_slot() = result.memory;
    last_emission_cache_slot() = result.emission_cache;
    release();   // nothing indexed outlives the call
    return result;
}

BlockedUnionResult score_diplotypes_indexed_blocked(const InternedHaplotypes& haplotypes,
                                                    const std::vector<HaplotypeCandidate>& candidates,
                                                    const std::vector<ReadPair>& fragments, const LibraryModel& library,
                                                    const FragmentScoringOptions& options,
                                                    const std::vector<std::vector<std::uint32_t>>& segments,
                                                    std::size_t store_budget_bytes, bool fine_views, std::size_t block_fragments,
                                                    const std::function<void(BlockedBlockReport&, bool before)>& on_block) {
    using Clock = std::chrono::steady_clock;
    const auto since = [](Clock::time_point mark) { return std::chrono::duration<double>(Clock::now() - mark).count(); };
    // ---- the bulk adapter's validation, over the WHOLE run (splitting must never admit what bulk refuses)
    validate_library_model(library);
    if (block_fragments == 0) throw std::invalid_argument("a fragment block needs at least one read pair");
    if (options.minimum_piece_bp == 0) throw std::invalid_argument("minimum indexable piece length must be positive");
    if (options.max_piece_hits == 0 || options.max_proposed_starts_per_haplotype == 0 || options.max_physical_origins == 0 ||
        options.max_state_classes == 0 || options.max_mass_cells == 0 || options.max_latent_alignment_pairs == 0) {
        throw std::invalid_argument("every fragment-scoring resource cap must be positive");
    }
    const std::size_t unique = haplotypes.sequences.size(), nf = fragments.size();
    for (const HaplotypeCandidate& c : candidates) {
        if (c.first >= unique || c.second >= unique) throw std::invalid_argument("diplotype candidate names an unknown haplotype");
    }
    if (options.record_origin_summary || !options.record_origins_for_fragments.empty() || options.record_haplotype_work != nullptr) {
        throw std::invalid_argument("the indexed union backend records no per-origin or per-haplotype diagnostics");
    }
    if (segments.size() != unique) throw std::invalid_argument("the indexed union backend needs one segmentation per unique haplotype");
    if (store_budget_bytes == 0) throw std::invalid_argument("the indexed union backend needs a positive placement-store budget");
    if (options.max_union_cached_placements_consumed == 0) {
        throw std::invalid_argument("the indexed union backend needs a positive cached-placement consumption cap");
    }
    if (options.union_origin_cache_bytes > 0 &&
        (options.max_union_cached_origins_consumed == 0 || options.max_union_cached_states_consumed == 0)) {
        throw std::invalid_argument("union origin reuse needs both of its consumption caps, positive");
    }
    if (options.emission_handoff_out != nullptr || options.emission_handoff_in != nullptr) {
        throw std::invalid_argument("the blocked union pass takes no emission handoff");
    }
    BlockedUnionResult out;
    if (nf != 0 && unique > options.max_mass_cells / nf) {   // the whole run's matrix, exactly as bulk refuses it
        out.refusal = "the fragment-by-haplotype mass matrix would exceed its declared cell cap; no partial ranking is returned";
        return out;
    }
    if (options.union_prefilter_cache_bytes > 0) {   // the cache needs one read length over the WHOLE run, as bulk requires
        std::size_t length = 0;
        bool single = true;
        for (const ReadPair& f : fragments) {
            for (const std::size_t l : {f.first.bases.size(), f.second.bases.size()}) {
                if (length == 0) length = l;
                else if (l != length) single = false;
            }
        }
        if (!single) throw std::invalid_argument("union prefilter reuse needs one read length");
    }
    std::size_t max_read_bp = 1;
    for (const ReadPair& f : fragments) max_read_bp = std::max({max_read_bp, f.first.bases.size(), f.second.bases.size()});
    FragmentScoringOptions batch_options = options;
    batch_options.union_segment_lengths = nullptr;
    batch_options.union_placement_cache_bytes = 0;
    // ---- one index over every haplotype
    auto mark = Clock::now();
    ContentPlacementIndex index(max_read_bp, library.max_gap_bp, true);
    for (std::size_t h = 0; h < unique; ++h) {
        if (index.add_haplotype(haplotypes.sequences[h], segments[h]) != h) throw std::logic_error("slots must follow the interned order");
    }
    out.seconds_index_build = since(mark);
    out.index = index.report();
    out.index_unit_bytes = out.index.unit_bytes;
    // ---- whole-run quantities: exposure from every fragment's floor; background per fragment
    const DiscreteInsertDistribution library_insert = build_library_insert_distribution(library);
    std::map<std::size_t, std::size_t> by_floor;
    std::vector<double> background(nf);
    for (std::size_t f = 0; f < nf; ++f) {
        ++by_floor[std::max(fragments[f].first.bases.size(), fragments[f].second.bases.size())];
        background[f] = log_background_probability(fragments[f]);
    }
    std::vector<double> exposure(unique);
    for (std::size_t h = 0; h < unique; ++h) exposure[h] = exposure_for_length(haplotypes.sequences[h].size(), nf, by_floor, library_insert);
    // ---- per-candidate accumulators, seeded exactly as the bulk combination seeds them
    const std::vector<AuditCombination> audits = audit_combinations(options.audit_libraries);
    const std::vector<AuditCombination> main_terms = options.prepared_library_terms && !candidates.empty()
                                                         ? audit_combinations({library})
                                                         : std::vector<AuditCombination>{};
    const AuditCombination* main = main_terms.empty() ? nullptr : &main_terms.front();
    const std::size_t nc = candidates.size(), views = audits.size();
    std::vector<DiplotypeScore> scores(nc);
    std::vector<double> fragment_term(nc, 0.0);
    std::vector<std::vector<DiplotypeScore>> audit_scores(views, std::vector<DiplotypeScore>(nc));
    std::vector<std::vector<double>> audit_term(views, std::vector<double>(nc, 0.0));
    for (std::size_t c = 0; c < nc; ++c) {
        const std::size_t a = candidates[c].first, b = candidates[c].second;
        DiplotypeScore& score = scores[c];
        score.exposure_term = main != nullptr ? -(main->exposure_scale * (exposure[a] + exposure[b]))
                                              : -expected_signal_fragments(exposure[a] + exposure[b], library);
        score.absolute_term_sum = std::abs(score.exposure_term);
        for (std::size_t k = 0; k < views; ++k) {
            audit_scores[k][c].exposure_term = -(audits[k].exposure_scale * (exposure[a] + exposure[b]));
            audit_scores[k][c].absolute_term_sum = std::abs(audit_scores[k][c].exposure_term);
        }
    }
    out.accumulator_bytes = nc * (sizeof(DiplotypeScore) + sizeof(double)) * (1 + views) + exposure.capacity() * sizeof(double) +
                            background.capacity() * sizeof(double);
    const std::size_t threads = std::max<std::size_t>(1, std::min(options.threads, std::max<std::size_t>(1, nc)));
    UnionCumulativeCharges cumulative;
    FragmentScoringResult& result = out.result;
    std::uint64_t floor_binding = 0;
    try {
        for (std::size_t begin = 0; begin < nf; begin += block_fragments) {
            const std::size_t end = std::min(nf, begin + block_fragments);
            BlockedBlockReport rep;
            rep.first_fragment = begin;
            rep.fragments = end - begin;
            if (on_block) on_block(rep, true);
            const auto block_mark = Clock::now();
            {
                ContentIndexBatch batch(index, std::vector<ReadPair>(fragments.begin() + begin, fragments.begin() + end), library,
                                        batch_options);
                batch.set_placement_store_budget(store_budget_bytes);
                batch.set_fine_views(fine_views);
                batch.set_one_shot(true);
                batch.set_cumulative_baseline(cumulative);
                batch.discover();
                batch.build_columns();
                const std::vector<std::vector<double>>& column = batch.mass_columns();
                // fold this block into every candidate, fragments in their original order (the bulk loop's operations)
                const auto combine_mark = Clock::now();
                const auto fold = [&](std::size_t lo, std::size_t hi) {
                    for (std::size_t c = lo; c < hi; ++c) {
                        const std::vector<double>& ma = column[candidates[c].first];
                        const std::vector<double>& mb = column[candidates[c].second];
                        DiplotypeScore& score = scores[c];
                        double term = fragment_term[c];
                        for (std::size_t i = 0; i < end - begin; ++i) {
                            const double combined = log_add_exp(ma[i], mb[i]);
                            if (combined != kNegInf) ++score.fragments_with_mass;
                            const double contribution = main != nullptr ? prepared_contribution(*main, combined, background[begin + i])
                                                                        : log_fragment_contribution(combined, background[begin + i], library);
                            term += contribution;
                            score.absolute_term_sum += std::abs(contribution);
                            for (std::size_t k = 0; k < views; ++k) {
                                const double value = prepared_contribution(audits[k], combined, background[begin + i]);
                                audit_term[k][c] += value;
                                audit_scores[k][c].absolute_term_sum += std::abs(value);
                            }
                        }
                        fragment_term[c] = term;
                    }
                };
                const std::size_t span = (nc + threads - 1) / threads;
                std::vector<std::thread> pool;
                for (std::size_t t = 0; t < threads; ++t) {
                    const std::size_t lo = std::min(t * span, nc), hi = std::min(lo + span, nc);
                    if (lo < hi) pool.emplace_back(fold, lo, hi);
                }
                for (std::thread& w : pool) w.join();
                out.seconds_combination += since(combine_mark);
                const ContentBatchReport& br = batch.report();
                const UnionCumulativeCharges own = batch.charges();
                add_counters(cumulative.executed, own.executed);
                cumulative.piece_hits += own.piece_hits;
                cumulative.physical_origins += own.physical_origins;
                cumulative.placements_consumed += own.placements_consumed;
                cumulative.origins_consumed += own.origins_consumed;
                cumulative.states_consumed += own.states_consumed;
                cumulative.prefilters_consumed += own.prefilters_consumed;
                floor_binding += batch.state_->prepared.fragments_with_floor_binding;
                rep.units_evaluated = br.units_evaluated;
                rep.store_peak_bytes = br.store_peak_bytes;
                rep.origin_held_bytes = br.origin_held_bytes;
                rep.retained_column_bytes = br.retained_column_bytes;
                rep.prefilter_held_bytes = br.prefilter.window_index_bytes + br.prefilter.decision_table_bytes;
                rep.emission_cache_bytes = br.emission_cache_table_bytes + br.emission_cache_slab_bytes + br.emission_cache_component_bytes;
                rep.emission_cache_entries = br.emission_cache_entries_held;
                rep.worker_scratch_bytes = br.worker_scratch_bytes;
                rep.executed_dps = br.executed_dps;
                rep.executed_prefilters = br.executed_prefilters;
                rep.direct_physical_origins = br.direct_physical_origins;
                rep.served_physical_origins = br.served_physical_origins;
                rep.origins_constructed = br.origins_constructed;
                rep.origin_entities_admitted = br.origin_entities_admitted;
                rep.origin_entities_declined = br.origin_entities_declined;
                rep.cumulative_physical_origins = cumulative.physical_origins;
                rep.cumulative_piece_hits = cumulative.piece_hits;
                result.union_origin_cache.origins_constructed += br.origins_constructed;
                result.union_origin_cache.origins_consumed += br.origins_consumed_last;
                result.union_origin_cache.hits += br.zones_applied;
                result.union_origin_cache.admitted_segments += br.origin_entities_admitted;
                result.union_origin_cache.declined_segments += br.origin_entities_declined;
                result.union_origin_cache.peak_bytes = std::max(result.union_origin_cache.peak_bytes, br.origin_held_bytes);
                result.emission_cache.executed_dps += br.emission_cache.executed_dps;
                result.emission_cache.executed_prefilters += br.emission_cache.executed_prefilters;
                result.emission_cache.dp_cache_hits += br.emission_cache.dp_cache_hits;
                result.emission_cache.dp_cache_misses += br.emission_cache.dp_cache_misses;
                result.emission_cache.dp_cache_requests += br.emission_cache.dp_cache_requests;
                result.threads_used_placement = std::max(result.threads_used_placement,
                                                         std::max<std::size_t>(br.threads_used_discovery, br.threads_used_pairing));
            }   // the block's batch -- its placements, aggregates, caches and columns -- is released here
            rep.seconds = since(block_mark);
            if (on_block) on_block(rep, false);
            out.blocks.push_back(rep);
        }
    } catch (const ResourceExhausted& error) {
        out.refusal = error.what();   // no union ranking: the accumulated partial scores are discarded
        result.scores.clear();
        result.audit_scores.clear();
        result.counters = cumulative.executed;
        return out;
    }
    for (std::size_t c = 0; c < nc; ++c) {
        scores[c].fragment_term = fragment_term[c];
        scores[c].log_likelihood = scores[c].exposure_term + fragment_term[c];
        for (std::size_t k = 0; k < views; ++k) {
            audit_scores[k][c].fragment_term = audit_term[k][c];
            audit_scores[k][c].log_likelihood = audit_scores[k][c].exposure_term + audit_term[k][c];
            audit_scores[k][c].fragments_with_mass = scores[c].fragments_with_mass;
        }
    }
    if (main != nullptr) {   // the I2-F tally, once for the whole combination, as bulk records it
        std::uint64_t with_mass = 0;
        for (const DiplotypeScore& score : scores) with_mass += score.fragments_with_mass;
        note_prepared_combination(nc, nf, with_mass, *main);
    }
    result.scores = std::move(scores);
    result.audit_scores = std::move(audit_scores);
    result.counters = cumulative.executed;
    result.counters.fragments_with_floor_binding = floor_binding;
    result.unique_haplotypes = unique;
    result.candidate_fragment_evaluations = static_cast<std::uint64_t>(nc) * static_cast<std::uint64_t>(nf);
    result.threads_used_combination = threads;
    result.seconds_diploid_combination = out.seconds_combination;
    out.completed = true;
    return out;
}

ProgressiveResult score_diplotypes_indexed_progressive(const InternedHaplotypes& haplotypes,
                                                       const std::vector<std::vector<HaplotypeCandidate>>& stages,
                                                       const std::vector<ReadPair>& fragments, const LibraryModel& library,
                                                       const FragmentScoringOptions& options,
                                                       const std::vector<std::vector<std::uint32_t>>& segments,
                                                       std::size_t store_budget_bytes, bool fine_views,
                                                       const std::function<void(ProgressiveStageReport&, const FragmentScoringResult&)>& on_stage,
                                                       bool persistent_emissions) {
    using Clock = std::chrono::steady_clock;
    const auto since = [](Clock::time_point mark) { return std::chrono::duration<double>(Clock::now() - mark).count(); };
    // the bulk adapter's validation, unchanged
    validate_library_model(library);
    if (options.minimum_piece_bp == 0) throw std::invalid_argument("minimum indexable piece length must be positive");
    if (options.max_piece_hits == 0 || options.max_proposed_starts_per_haplotype == 0 || options.max_physical_origins == 0 ||
        options.max_state_classes == 0 || options.max_mass_cells == 0 || options.max_latent_alignment_pairs == 0) {
        throw std::invalid_argument("every fragment-scoring resource cap must be positive");
    }
    const std::size_t unique = haplotypes.sequences.size();
    for (const auto& stage : stages) {
        for (const HaplotypeCandidate& c : stage) {
            if (c.first >= unique || c.second >= unique) throw std::invalid_argument("diplotype candidate names an unknown haplotype");
        }
    }
    if (options.record_origin_summary || !options.record_origins_for_fragments.empty() || options.record_haplotype_work != nullptr) {
        throw std::invalid_argument("the indexed union backend records no per-origin or per-haplotype diagnostics");
    }
    if (segments.size() != unique) throw std::invalid_argument("the indexed union backend needs one segmentation per unique haplotype");
    if (store_budget_bytes == 0) throw std::invalid_argument("the indexed union backend needs a positive placement-store budget");
    if (options.max_union_cached_placements_consumed == 0) {
        throw std::invalid_argument("the indexed union backend needs a positive cached-placement consumption cap");
    }
    if (options.union_origin_cache_bytes > 0 &&
        (options.max_union_cached_origins_consumed == 0 || options.max_union_cached_states_consumed == 0)) {
        throw std::invalid_argument("union origin reuse needs both of its consumption caps, positive");
    }
    std::size_t max_read_bp = 1;
    for (const ReadPair& f : fragments) max_read_bp = std::max({max_read_bp, f.first.bases.size(), f.second.bases.size()});
    FragmentScoringOptions batch_options = options;
    batch_options.union_segment_lengths = nullptr;
    batch_options.union_placement_cache_bytes = 0;
    ProgressiveResult out;
    const auto started = Clock::now();
    auto mark = Clock::now();
    ContentPlacementIndex index(max_read_bp, library.max_gap_bp, true);
    for (std::size_t h = 0; h < unique; ++h) {
        if (index.add_haplotype(haplotypes.sequences[h], segments[h]) != h) throw std::logic_error("slots must follow the interned order");
    }
    out.seconds_index_build = since(mark);
    mark = Clock::now();
    ContentIndexBatch batch(index, fragments, library, batch_options);
    out.seconds_preparation = since(mark);
    if (options.union_prefilter_cache_bytes > 0 && batch.report().prefilter_status != "on") {
        throw std::invalid_argument("union prefilter reuse needs one read length");
    }
    batch.set_placement_store_budget(store_budget_bytes);
    batch.set_fine_views(fine_views);
    batch.set_persistent_emission_cache(persistent_emissions);   // S2-D (off: S2-C's per-call cache)
    {
        StartSortTally& tally = start_sort_tally();
        tally.out_of_place.store(0);
        tally.in_place.store(0);
        tally.fallbacks.store(0);
        tally.peak_worker_bytes.store(0);
    }
    double seconds_combination = 0.0, seconds_placement = 0.0;
    try {
        for (const auto& stage : stages) {
            ProgressiveStageReport rep;
            rep.candidates = stage.size();
            std::vector<std::size_t> slots;
            for (const HaplotypeCandidate& c : stage) { slots.push_back(c.first); slots.push_back(c.second); }
            std::sort(slots.begin(), slots.end());
            slots.erase(std::unique(slots.begin(), slots.end()), slots.end());
            rep.slots = slots.size();
            const ContentBatchReport before = batch.report();
            const std::size_t units_before = batch.units_discovered(), columns_before = batch.columns_completed();
            std::set<std::uint32_t> needed;
            for (const std::size_t h : slots) {
                for (int strand = 0; strand < 2; ++strand) for (const IndexOccurrence& occ : index.state_->occurrences[strand][h]) needed.insert(occ.unit);
            }
            for (const std::uint32_t u : needed) {   // the stage's occurrence units already held before it
                rep.units_reused += u < batch.state_->unit_done.size() && batch.state_->unit_done[u] ? 1 : 0;
            }
            const auto d0 = Clock::now();
            batch.discover_for(slots);
            rep.seconds_discovery = since(d0);
            rep.units_discovered = batch.units_discovered() - units_before;   // coarse units and fine dependencies
            const auto s0 = Clock::now();
            FragmentScoringResult r = batch.score_selected(stage);
            rep.seconds_scoring = since(s0);
            const ContentBatchReport& after = batch.report();
            rep.columns_computed = batch.columns_completed() - columns_before;
            rep.columns_reused = slots.size() - rep.columns_computed;
            rep.pairs_computed = after.pairs_computed_last;
            rep.pairs_served = stage.size() - std::min<std::uint64_t>(stage.size(), after.pairs_computed_last);
            rep.executed_prefilters = after.executed_prefilters - before.executed_prefilters;
            rep.executed_dps = after.executed_dps - before.executed_dps;
            rep.full_read_verifications = after.full_read_verifications - before.full_read_verifications;
            rep.piece_hits = after.piece_hits - before.piece_hits;
            rep.direct_physical_origins = after.direct_physical_origins - before.direct_physical_origins;
            rep.served_physical_origins = after.served_physical_origins - before.served_physical_origins;
            rep.origins_constructed = after.origins_constructed - before.origins_constructed;
            rep.expanded_placements = after.expanded_placements - before.expanded_placements;
            rep.store_live_bytes = after.store_live_bytes;
            rep.store_peak_bytes = after.store_peak_bytes;
            rep.retained_column_bytes = after.retained_column_bytes;
            rep.origin_held_bytes = after.origin_held_bytes;
            rep.worker_scratch_bytes = after.worker_scratch_bytes;
            rep.dp_cache_hits = after.emission_cache.dp_cache_hits - before.emission_cache.dp_cache_hits;
            rep.dp_cache_misses = after.emission_cache.dp_cache_misses - before.emission_cache.dp_cache_misses;
            rep.carried_hits = after.persistent_inherited_hits - before.persistent_inherited_hits;
            rep.carried_bytes = after.persistent_carried_bytes;
            rep.carried_entries = after.persistent_carried_entries;
            rep.candidate_fragment_evaluations = r.candidate_fragment_evaluations;
            // footprint_* are left to the caller's on_stage hook (this layer does not sample process memory)
            rep.seconds_since_start = since(started);
            seconds_combination += r.seconds_diploid_combination;
            seconds_placement += rep.seconds_discovery + rep.seconds_scoring - r.seconds_diploid_combination;
            out.summary.candidate_fragment_evaluations += r.candidate_fragment_evaluations;
            out.summary.threads_used_combination = std::max(out.summary.threads_used_combination, r.threads_used_combination);
            if (on_stage) on_stage(rep, r);
            out.scores.push_back(std::move(r.scores));
            out.audit_scores.push_back(std::move(r.audit_scores));
            out.stages.push_back(rep);
        }
        out.completed = true;
    } catch (const ResourceExhausted& error) {
        out.refusal = error.what();   // the completed stages keep their scores; nothing partial is scored
    }
    batch.release_carried_emissions();   // S2-D: the carried cache ends with the schedule
    out.index = index.report();
    out.batch = batch.report();
    const ContentBatchReport& b = out.batch;
    FragmentScoringResult& sum = out.summary;
    sum.counters = b.counters_last;
    sum.unique_haplotypes = unique;
    sum.union_origin_cache = indexed_origin_report(b, options);
    sum.union_prefilter_cache = b.prefilter;
    {
        const StartSortTally& tally = start_sort_tally();
        sum.start_sort.budget_bytes = options.start_sort_scratch_bytes;
        sum.start_sort.workers = std::max<std::size_t>(1, options.threads);
        sum.start_sort.share_bytes = options.start_sort_scratch_bytes / sum.start_sort.workers;
        sum.start_sort.out_of_place = tally.out_of_place.load();
        sum.start_sort.in_place = tally.in_place.load();
        sum.start_sort.fallbacks = tally.fallbacks.load();
        sum.start_sort.peak_worker_bytes = tally.peak_worker_bytes.load();
    }
    sum.emission_cache = b.emission_cache;
    sum.seconds_fragment_preparation = out.seconds_preparation;
    sum.seconds_placement_and_mass = out.seconds_index_build + seconds_placement;
    sum.seconds_diploid_combination = seconds_combination;
    sum.threads_used_placement = std::max<std::size_t>(1, std::max(b.threads_used_discovery, b.threads_used_pairing));
    return out;
}

} // namespace panvar
