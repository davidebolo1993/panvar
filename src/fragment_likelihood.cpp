#include "panvar/fragment_likelihood.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <exception>
#include <functional>
#include <map>
#include <memory>
#include <deque>
#include <thread>
#include <set>
#include <tuple>
#include <type_traits>
#include <array>
#include <cmath>
#include <cstring>
#include <mutex>
#include <stdexcept>
#include <unordered_map>

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

// The hit path: one mutex, one probe chain, one memcmp pair, and a component copy into a buffer the
// caller already owns. No allocation.
CacheProbe cache_lookup(const KeyView& key, GappedFit* out, EmissionWork& work) {
    EmissionCache& cache = emission_cache();
    if (!cache.active.load(std::memory_order_relaxed)) return CacheProbe::Disabled;
    const auto mark = std::chrono::steady_clock::now();
    CacheShard& shard = cache.shards[key.h0 % kCacheShards];
    CacheProbe answer = CacheProbe::Miss;
    {
        const std::lock_guard<std::mutex> guard(shard.guard);
        const std::size_t at = find_occupied_slot(shard, key);
        if (at == kNoSlot) {
            ++work.dp_cache_misses;
        } else {
            ++work.dp_cache_hits;
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
    work.cache_lookup_nanos += static_cast<std::uint64_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - mark)
            .count());
    return answer;
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

bool myers_shw_within(const char* query, std::size_t n, const char* target, std::size_t m,
                      std::size_t band) {
    MyersWord peq[4] = {0, 0, 0, 0};
    for (std::size_t i = 0; i < n; ++i) peq[two_bit_code(query[i])] |= MyersWord{1} << i;
    return myers_shw_run(peq, n, target, m, band);
}

MyersWord myers_word_of(std::uint64_t low, std::uint64_t high) {
    return (static_cast<MyersWord>(high) << 64) | static_cast<MyersWord>(low);
}

bool feasible_within_band(const std::string& haplotype, std::size_t start, std::size_t window,
                          const std::string& read, std::size_t band,
                          const MyersPattern* pattern) {
    if (window == 0) return false;
    // Prepared: the masks were built once for this form during fragment preparation. Rebuilding
    // them here cost a pass over the read on every one of 1.19 billion verifications.
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

double emission_cache_lookup_seconds() {
    return static_cast<double>(emission_cache().lookup_nanos.load(std::memory_order_relaxed)) / 1e9;
}

void reset_emission_cache_lookup_seconds() {
    emission_cache().lookup_nanos.store(0, std::memory_order_relaxed);
}

EmissionCacheCounters end_emission_cache() {
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
    for (CacheShard& shard : cache.shards) {
        const std::lock_guard<std::mutex> guard(shard.guard);
        shard.table.clear();
        shard.occupied = 0;
        shard.slabs.clear();
        shard.slab_used = 0;
        shard.arena.clear();
    }
    return out;
}

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

std::optional<GappedFit> gapped_emission(
    const std::string& haplotype,
    const std::string& read,
    std::size_t start,
    const LibraryModel& library,
    std::size_t mismatch_band,
    const MyersPattern* pattern) {
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
    require_acgt_window(haplotype, start, window, read);

    EmissionWork& work = emission_work();
    EmissionCache& cache = emission_cache();

    // ORDER MATTERS, and this is the correction. The exact feasibility prefilter runs FIRST, before
    // any key is constructed or looked up. It costs ~1.6us; forming a key hashes ~203 bytes and then
    // probes a shared table. Refusals are the overwhelming majority of requests (993M of 1,189M at
    // 30x), so making them pay for a key -- and then retaining them in the same representation as a
    // recursion result -- spent the whole budget on the cheap stream and collapsed the hit rate to
    // 37.8%. A refusal now returns here, having touched no table at all.
    ++work.executed_prefilters;
    if (!feasible_within_band(haplotype, start, window, read, mismatch_band, pattern)) {
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
        key = key_view_of(read, haplotype.data() + start, window);
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
    }

    ++work.executed_dps;

    GappedFit fit;
    // Filled in place rather than through a scratch vector, so an accepted placement costs one
    // allocation instead of two.
    const double raw_total = gap_recursion_scaled(&haplotype, &read, start, read.size(), available,
                                                  mismatch_band, library, true, &fit.components);
    if (raw_total == kNegInf || fit.components.empty()) {
        // The recursion RAN and found nothing. That is still work a cache would avoid, so it counts
        // against the DP stream rather than the prefilter's -- and it is retained as a refusal,
        // because recomputing it costs a whole recursion.
        if (caching) cache_store_refusal(key, CacheOutcome::RefusedAfterRecursion);
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
    if (caching) cache_store_fit(key, fit);
    return fit;
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
    std::vector<std::uint32_t>* touched = nullptr) {
    proposed.clear();   // keeps capacity: the point of hoisting it
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

    // Sorting and full-read verification are the rest of this haplotype's cost and neither touches the
    // budget, so without these checks a worker cancelled by someone else's refusal would still pay for
    // both before noticing. Prompt cancellation has to mean prompt in the expensive loops too.
    if (budget.stop()) throw PlacementCancelled();

    // Several pieces of one read can rediscover one physical placement. Deduplicate before
    // verifying, so a placement found three ways is verified once and counted once.
    std::sort(proposed.begin(), proposed.end());
    proposed.erase(std::unique(proposed.begin(), proposed.end()), proposed.end());
    counters.unique_starts += proposed.size();
    if (budget.stop()) throw PlacementCancelled();

    for (const auto& item : proposed) {
        // One relaxed load against a full-read comparison: free at this granularity, and verification is
        // where a cancelled worker would otherwise spend the longest.
        if (budget.stop()) throw PlacementCancelled();
        ++counters.full_read_verifications;
        const std::optional<GappedFit> fit =
            gapped_emission(haplotype, forms[item.first].bases, item.second, library, band,
                            &forms[item.first].pattern);
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

// Diagnostics only; the context leaves every pointer null.
struct OriginRecording {
    bool summary = false;
    const std::vector<char>* record_origins = nullptr;
    std::size_t unique_stride = 0;
    std::vector<OriginSummary>* summaries = nullptr;
    std::map<std::pair<std::size_t, std::size_t>, std::vector<FragmentPlacement>>* recorded = nullptr;
};

// Places haplotypes [first_h, last_h) and writes their mass COLUMNS. Column storage rather than the
// row-major matrix is what makes a later batch appendable: a new haplotype adds a column and moves
// nothing, whereas a row-major matrix restrides every existing cell when the haplotype count grows.
// The arithmetic per cell is untouched, so the flattened result is identical to the matrix the
// one-shot path used to build in place.
void place_haplotype_range(
    const PlacementInputs& in,
    std::size_t first_h,
    std::size_t last_h,
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

    for (std::size_t h = first_h; h < last_h; ++h) {
        if (budget.stop()) throw PlacementCancelled();
        const std::string& haplotype = (*in.sequences)[h];
        if (prepared.forward_forms.empty()) continue;

        for (std::vector<FormPlacement>& row : forward_placements) row.clear();
        const std::uint64_t rows_before = counters.accepted_mate_placements;
        collect_form_placements(haplotype, prepared.forward_forms, prepared.forward_componentise,
                                prepared.forward_indexes, budget, band, library, forward_placements,
                                proposed_scratch, counters, cache_work, options);
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
                                options);
        memory.proposed_starts_peak =
            std::max<std::uint64_t>(memory.proposed_starts_peak, proposed_scratch.size());
        memory.proposed_starts_capacity_bytes_peak = std::max<std::uint64_t>(
            memory.proposed_starts_capacity_bytes_peak,
            proposed_scratch.capacity() * sizeof(std::pair<std::uint32_t, std::uint32_t>));
        memory.placement_rows_peak = std::max<std::uint64_t>(
            memory.placement_rows_peak, counters.accepted_mate_placements - rows_before);

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

        for (std::size_t f = 0; f < fragments.size(); ++f) {
            if (budget.stop()) throw PlacementCancelled();
            const InsertSupport& support = prepared.support_by_fragment[f];
            if (support.empty()) continue;
            states.clear();

            for (std::size_t arm = 0; arm < 2; ++arm) {
                const std::vector<FormPlacement>& forward = forward_placements[form_slot(f, arm)];
                const std::vector<FormPlacement>& reverse = reverse_placements[form_slot(f, arm)];
                if (forward.empty() || reverse.empty()) continue;

                for (std::size_t first_row = 0; first_row < forward.size();) {
                    const std::uint32_t forward_start = forward[first_row].start;
                    std::size_t last_row = first_row;
                    while (last_row < forward.size() && forward[last_row].start == forward_start) {
                        ++last_row;
                    }
                    const std::size_t lo_end =
                        static_cast<std::size_t>(forward_start) + support.min_bp;
                    const std::size_t hi_end =
                        static_cast<std::size_t>(forward_start) + support.max_bp;
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
                                forward_start, forward[row].consumed, it->start, it->consumed,
                                support);
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
                                "fragment scoring reached its physical-origin cap; no partial "
                                "ranking is returned");
                        }
                        if (budget.stop()) throw PlacementCancelled();
                        StateKey key;
                        key.insert_bp = item.first;
                        key.orientation = static_cast<std::uint32_t>(arm);
                        key.log_latent_mass = item.second.log_latent_mass;
                        ++states[key];
                        if (census != nullptr) {
                            census->origins_by_fragment[f].push_back(
                                {forward_start, item.first, static_cast<std::uint32_t>(arm),
                                 item.second.log_latent_mass});
                        }
                        const bool want_named =
                            origins_out.record_origins != nullptr &&
                            (*origins_out.record_origins)[f] != 0;
                        if (origins_out.summary || want_named) {
                            FragmentPlacement placement;
                            placement.forward_start = forward_start;
                            placement.insert_bp = item.first;
                            placement.log_latent_mass = item.second.log_latent_mass;
                            placement.latent_alignments = item.second.latent_alignments;
                            placement.reverse_start = item.second.reverse_start;
                            placement.forward_mismatches = item.second.forward_mismatches;
                            placement.reverse_mismatches = item.second.reverse_mismatches;
                            placement.forward_consumed_bp = item.second.forward_consumed_bp;
                            placement.reverse_consumed_bp = item.second.reverse_consumed_bp;
                            placement.orientation = arm == 0 ? MateOrientation::FirstForward
                                                             : MateOrientation::SecondForward;
                            if (origins_out.summary && origins_out.summaries != nullptr) {
                                OriginSummary& summary =
                                    (*origins_out.summaries)[f * origins_out.unique_stride + h];
                                ++summary.origins;
                                const double term = kLogMateOrientationPrior +
                                                    library_insert.log_probability(item.first) +
                                                    item.second.log_latent_mass;
                                if (term > summary.best_log_term) {
                                    summary.best_log_term = term;
                                    summary.best = placement;
                                }
                            }
                            if (want_named && origins_out.recorded != nullptr) {
                                (*origins_out.recorded)[{f, h}].push_back(placement);
                            }
                        }
                    }
                    first_row = last_row;
                }
            }
            if (states.empty()) continue;
            counters.statistical_state_classes += states.size();
            mass_column[h][f] = reduce_states(states, library_insert);
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
        place_haplotype_range(in, first_h, last_h, budget, mass_column, counters, memory, cache_work,
                              origins_out, census_out);
        return threads;
    }
    const std::size_t span = (span_count + threads - 1) / threads;
    std::vector<FragmentScoringCounters> local(threads);
    std::vector<ScorerMemoryAccounts> local_memory(threads);
    std::vector<EmissionCacheCounters> local_cache(threads);
    std::vector<std::exception_ptr> failures(threads);
    std::vector<std::thread> workers;
    workers.reserve(threads);
    for (std::size_t t = 0; t < threads; ++t) {
        const std::size_t lo = std::min(first_h + t * span, last_h);
        const std::size_t hi = std::min(lo + span, last_h);
        workers.emplace_back([&, t, lo, hi]() {
            try {
                place_haplotype_range(in, lo, hi, budget, mass_column, local[t], local_memory[t],
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
// Transient by construction: these bytes live only for one delta. That is precisely why the cache
// must be frozen before any of this runs -- a key built here would point into `bases` below, which
// is about to go out of scope.
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
    WindowScratch& scratch) {
    const auto form_slot = [](std::size_t fragment, std::size_t arm) { return fragment * 2 + arm; };
    scratch.ensure_rows(prepared.forward_forms.size(), prepared.reverse_forms.size());
    scratch.clear_touched();

    collect_form_placements(window.bases, prepared.forward_forms, prepared.forward_componentise,
                            prepared.forward_indexes, budget, band, library, scratch.forward_rows,
                            scratch.proposed, counters, cache_work, options, &scratch.touched_forward);
    collect_form_placements(window.reversed, prepared.reverse_forms, prepared.reverse_componentise,
                            prepared.reverse_indexes, budget, band, library, scratch.reverse_rows,
                            scratch.proposed, counters, cache_work, options, &scratch.touched_reverse);

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
    for (std::size_t f = 0; f < fragments.size(); ++f) {
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
    scratch.offsets[fragments.size()] = scratch.origins.size();
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
    };

    result.threads_used_placement =
        (unique == 0 || origins_force_serial)
            ? 1
            : std::max<std::size_t>(1, std::min(options.threads, unique));
    try {
        result.threads_used_placement =
            run_placement(inputs, 0, unique, options.threads, origins_force_serial, budget,
                          mass_column, counters, result.memory, result.emission_cache, origins_out);
    } catch (const ResourceExhausted&) {
        note_observed_work();
        throw;
    }

    enforce_acceptance_caps(counters, options, note_observed_work);

    result.seconds_placement_and_mass =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - placement_mark).count();

    const auto combination_mark = std::chrono::steady_clock::now();
    const auto combine_range = [&](std::size_t begin, std::size_t end) {
        for (std::size_t c = begin; c < end; ++c) {
            const std::size_t a = candidates[c].first;
            const std::size_t b = candidates[c].second;
            DiplotypeScore& score = result.scores[c];
            score.exposure_term = -expected_signal_fragments(exposure[a] + exposure[b], library);
            score.absolute_term_sum = std::abs(score.exposure_term);
            double fragment_term = 0.0;
            for (std::size_t f = 0; f < fragments.size(); ++f) {
                const double combined = log_add_exp(mass_column[a][f], mass_column[b][f]);
                if (combined != kNegInf) ++score.fragments_with_mass;
                const double contribution = log_fragment_contribution(
                    combined, prepared.background_by_fragment[f], library);
                fragment_term += contribution;
                score.absolute_term_sum += std::abs(contribution);
            }
            score.fragment_term = fragment_term;
            score.log_likelihood = score.exposure_term + fragment_term;
        }
    };

    const std::size_t combination_threads =
        candidates.empty() ? 1
                           : std::max<std::size_t>(1, std::min(options.threads, candidates.size()));
    result.threads_used_combination = combination_threads;
    result.candidate_fragment_evaluations =
        static_cast<std::uint64_t>(candidates.size()) * static_cast<std::uint64_t>(fragments.size());

    if (combination_threads <= 1) {
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
    result.seconds_diploid_combination =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - combination_mark).count();

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

// ---------------------------------------------------------------------------------------------
// PreparedScoringContext
// ---------------------------------------------------------------------------------------------

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

    void require_usable() const {
        if (refused) {
            throw ResourceExhausted(
                "this scoring context refused a declared cap; the whole search is invalid and no "
                "partial ranking is available from it");
        }
    }
};

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
    scoring_detail::emission_cache_pinned().store(true, std::memory_order_relaxed);
}

PreparedScoringContext::~PreparedScoringContext() {
    if (state_ && state_->cache_pinned) {
        scoring_detail::emission_cache_pinned().store(false, std::memory_order_relaxed);
        state_->cumulative_cache = end_emission_cache();
    }
}

std::vector<std::size_t> PreparedScoringContext::add_haplotypes(
    const std::vector<std::string>& sequences,
    IncrementalBatchReport* report) {
    State& state = *state_;
    state.require_usable();
    const std::size_t before = state.sequences.size();

    std::vector<std::size_t> out;
    out.reserve(sequences.size());
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
        const std::size_t slot = state.sequences.size();
        state.by_hash.emplace(hash, slot);
        state.hash_by_slot.push_back(hash);
        state.sequences.push_back(sequence);
        state.reversed.emplace_back(sequence.rbegin(), sequence.rend());
        out.push_back(slot);
    }
    const std::size_t after = state.sequences.size();
    const std::size_t added = after - before;

    if (!state.fragments.empty() &&
        after > state.options.max_mass_cells / state.fragments.size()) {
        state.refused = true;
        throw ResourceExhausted(
            "the fragment-by-haplotype mass matrix would exceed its declared cell cap; no partial "
            "ranking is returned");
    }

    state.mass_column.resize(after, std::vector<double>(state.fragments.size(), kNegInf));
    state.exposure.resize(after, 0.0);
    for (std::size_t h = before; h < after; ++h) {
        state.exposure[h] =
            exposure_for_length(state.sequences[h].size(), state.fragments.size(),
                                state.prepared.fragments_by_floor, state.library_insert);
    }

    FragmentScoringCounters batch;
    EmissionCacheCounters batch_cache;
    ScorerMemoryAccounts batch_memory;
    std::size_t threads_used = 1;
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
        OriginRecording none;
        try {
            threads_used = run_placement(inputs, before, after, state.options.threads, false,
                                         state.budget, state.mass_column, batch, batch_memory,
                                         batch_cache, none);
        } catch (const ResourceExhausted&) {
            // A refusal invalidates the WHOLE search, not this batch: columns already retained were
            // computed under a budget that has now been exceeded, and ranking from them would present
            // a partial search as a complete one.
            state.refused = true;
            throw;
        }
    }
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - mark).count();

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
    const auto combine_range = [&](std::size_t begin, std::size_t end) {
        for (std::size_t c = begin; c < end; ++c) {
            const std::size_t a = candidates[c].first;
            const std::size_t b = candidates[c].second;
            DiplotypeScore& score = result.scores[c];
            score.exposure_term =
                -expected_signal_fragments(state.exposure[a] + state.exposure[b], state.library);
            score.absolute_term_sum = std::abs(score.exposure_term);
            double fragment_term = 0.0;
            for (std::size_t f = 0; f < fragments.size(); ++f) {
                const double combined =
                    log_add_exp(state.mass_column[a][f], state.mass_column[b][f]);
                if (combined != kNegInf) ++score.fragments_with_mass;
                const double contribution = log_fragment_contribution(
                    combined, state.prepared.background_by_fragment[f], state.library);
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
    return result;
}

std::size_t PreparedScoringContext::unique_haplotypes() const { return state_->sequences.size(); }

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
    try {
        place_haplotype_range(inputs, haplotype, haplotype + 1, state.budget, scratch, counters,
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

    // The margin beyond which the emission normaliser is saturated, so an unaffected origin keeps its
    // key under a length change. Taken at the LONGEST read, which only ever widens the affected set.
    std::size_t longest_read = 0;
    for (const ReadPair& fragment : state.fragments) {
        longest_read = std::max({longest_read, fragment.first.bases.size(), fragment.second.bases.size()});
    }
    const std::size_t margin = longest_read + state.library.max_gap_bp;
    const std::size_t max_insert = library_insert_support(state.library).max_bp;
    const std::size_t halo = max_insert + margin;
    const std::size_t shortest = std::min(old_length, new_length);

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
    window.begin = block_begin > halo ? block_begin - halo : 0;
    std::size_t window_end = std::min(new_length, new_end + halo);
    // A length change reaches the ends through the normaliser, so the window must cover whichever end
    // the replacement sits within `margin` of, or the origins there would be dropped by both halves.
    if (shift != 0 && block_begin < margin) window.begin = 0;
    if (shift != 0 && shortest < block_end + margin) window_end = new_length;
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

    FragmentScoringCounters counters;
    EmissionCacheCounters cache_work;
    const std::uint64_t rows_before = state.window_scratch.row_constructions;
    const std::uint64_t accumulators_before = state.window_scratch.accumulator_growths;
    const std::uint64_t origins_before = state.window_scratch.origin_growths;
    const auto origins_mark = std::chrono::steady_clock::now();
    try {
        enumerate_window_origins(state.prepared, state.fragments, state.library, state.options,
                                 state.options.mismatch_band, window, state.budget, counters,
                                 cache_work, state.window_scratch);
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
        slot = state.sequences.size();
        state.by_hash.emplace(repaired_hash, slot);
        state.hash_by_slot.push_back(repaired_hash);
        // Empty placeholders keep every slot-indexed container aligned; the bytes live in `segments`.
        state.sequences.emplace_back();
        state.reversed.emplace_back();
        state.segments.emplace(slot, repaired);
        state.mass_column.resize(slot + 1, std::vector<double>(state.fragments.size(), kNegInf));
        state.exposure.resize(slot + 1, 0.0);
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
    for (std::size_t slot = state.transient_watermark; slot < state.sequences.size(); ++slot) {
        if (retained.count(slot) != 0) continue;
        if (state.sequences[slot].empty() && state.mass_column[slot].empty() &&
            state.segments.count(slot) == 0) {
            continue;
        }
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
    }
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

} // namespace panvar
