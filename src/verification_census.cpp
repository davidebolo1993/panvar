#include "panvar/verification_census.hpp"

#include <algorithm>
#include <atomic>
#include <cstring>
#include <map>
#include <memory>
#include <mutex>
#include <sstream>

namespace panvar {

namespace {

constexpr std::size_t kShards = 64;
constexpr std::size_t kSlabBytes = std::size_t{16} << 10;

std::uint64_t splitmix(std::uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

std::uint64_t read_hash_of(std::uint64_t a, std::uint64_t b) { return splitmix(a ^ ((b << 31) | (b >> 33))); }

// One distinct (context, read, window) key and what its requests did.
struct Entry {
    std::uint64_t hash = 0;
    const char* window = nullptr;   // in the shard's arena; null marks an empty slot
    std::uint32_t read = 0;         // read-store id
    std::uint16_t window_length = 0;
    std::uint8_t context = 0;
    std::uint8_t groups = 0;        // bit g: requested by group g
    std::uint32_t requests[kCensusGroups] = {};
    std::uint32_t last_fragment = 0, last_unit = 0;
    std::uint8_t last_group = 0, kernel_group = 0, recursion_group = 0;
    bool kernel_seen = false, recursion_seen = false, missed_before = false, not_retained_before = false;
};

struct Arena {
    std::vector<std::unique_ptr<char[]>> slabs;
    std::size_t used = kSlabBytes;
};

struct Shard {
    std::mutex guard;
    std::vector<Entry> table;   // open addressing, at most half full
    std::size_t occupied = 0;
    Arena arena;
};

struct ReadSlot {
    std::uint64_t hash = 0;
    const char* bytes = nullptr;
    std::uint32_t length = 0;
    std::uint32_t id = 0;
};

struct Census {
    std::atomic<bool> active{false};
    std::atomic<std::uint64_t> generation{0};
    std::uint32_t period = 0;
    std::uint64_t budget = 0;
    std::atomic<std::uint64_t> charged{0};
    std::atomic<std::uint64_t> peak{0};
    std::atomic<bool> complete{true};
    std::atomic<bool> collapse{false};
    std::atomic<std::uint8_t> stage{0};
    std::atomic<std::uint32_t> next_unit{0};
    std::array<Shard, kShards> shards;
    std::mutex reads_guard;
    std::vector<ReadSlot> reads;
    std::size_t reads_occupied = 0;
    Arena read_arena;
    std::mutex contexts_guard;
    std::vector<CensusContext> contexts;
    std::mutex totals_guard;
    std::array<CensusExact, kCensusGroups> exact{};
    std::array<CensusSampled, kCensusGroups> sampled{};   // under each shard's lock, merged here: see record
};

Census& census() {
    static Census c;
    return c;
}

// Storage is charged BEFORE it is allocated; a refusal marks coverage incomplete.
bool charge(std::uint64_t bytes) {
    Census& c = census();
    const std::uint64_t now = c.charged.fetch_add(bytes, std::memory_order_relaxed) + bytes;
    if (now > c.budget) {
        c.charged.fetch_sub(bytes, std::memory_order_relaxed);
        c.complete.store(false, std::memory_order_relaxed);
        return false;
    }
    std::uint64_t seen = c.peak.load(std::memory_order_relaxed);
    while (now > seen && !c.peak.compare_exchange_weak(seen, now)) {
    }
    return true;
}

void uncharge(std::uint64_t bytes) { census().charged.fetch_sub(bytes, std::memory_order_relaxed); }

const char* arena_copy(Arena& arena, const char* bytes, std::size_t length) {
    if (length > kSlabBytes) return nullptr;
    if (arena.used + length > kSlabBytes) {
        if (!charge(kSlabBytes + sizeof(std::unique_ptr<char[]>))) return nullptr;
        arena.slabs.push_back(std::make_unique<char[]>(kSlabBytes));
        arena.used = 0;
    }
    char* at = arena.slabs.back().get() + arena.used;
    std::memcpy(at, bytes, length);
    arena.used += length;
    return at;
}

template <typename Slot>
bool grow_if_needed(std::vector<Slot>& table, std::size_t occupied, std::uint64_t (*hash_of)(const Slot&),
                    bool (*empty)(const Slot&)) {
    if (!table.empty() && (occupied + 1) * 2 <= table.size()) return true;
    const std::size_t next_size = table.empty() ? 256 : table.size() * 2;
    if (!charge(next_size * sizeof(Slot))) return false;
    std::vector<Slot> next(next_size);
    const std::size_t mask = next_size - 1;
    for (const Slot& slot : table) {
        if (empty(slot)) continue;
        std::size_t at = (hash_of(slot) * 0x9e3779b97f4a7c15ULL >> 7) & mask;
        while (!empty(next[at])) at = (at + 1) & mask;
        next[at] = slot;
    }
    uncharge(table.size() * sizeof(Slot));
    table.swap(next);
    return true;
}

std::uint64_t entry_hash(const Entry& e) { return e.hash; }
bool entry_empty(const Entry& e) { return e.window == nullptr; }
std::uint64_t read_slot_hash(const ReadSlot& s) { return s.hash; }
bool read_slot_empty(const ReadSlot& s) { return s.bytes == nullptr; }

// The read's store id, by exact bytes (hash only addresses). -1: refused for budget.
long read_id(Census& c, std::uint64_t hash, const char* bytes, std::size_t length) {
    const std::lock_guard<std::mutex> lock(c.reads_guard);
    if (!c.reads.empty()) {
        const std::size_t mask = c.reads.size() - 1;
        for (std::size_t at = (hash * 0x9e3779b97f4a7c15ULL >> 7) & mask; c.reads[at].bytes != nullptr; at = (at + 1) & mask) {
            const ReadSlot& s = c.reads[at];
            if (s.hash == hash && s.length == length && std::memcmp(s.bytes, bytes, length) == 0) return s.id;
        }
    }
    if (!grow_if_needed(c.reads, c.reads_occupied, &read_slot_hash, &read_slot_empty)) return -1;
    const char* copy = arena_copy(c.read_arena, bytes, length);
    if (copy == nullptr) return -1;
    const std::size_t mask = c.reads.size() - 1;
    std::size_t at = (hash * 0x9e3779b97f4a7c15ULL >> 7) & mask;
    while (c.reads[at].bytes != nullptr) at = (at + 1) & mask;
    c.reads[at] = {hash, copy, static_cast<std::uint32_t>(length), static_cast<std::uint32_t>(c.reads_occupied)};
    return static_cast<long>(c.reads_occupied++);
}

int context_id(Census& c, const CensusContext& x) {
    const std::lock_guard<std::mutex> lock(c.contexts_guard);
    for (std::size_t i = 0; i < c.contexts.size(); ++i) {
        const CensusContext& y = c.contexts[i];
        if (std::memcmp(&x.substitution_error, &y.substitution_error, sizeof(double)) == 0 &&
            std::memcmp(&x.gap_open, &y.gap_open, sizeof(double)) == 0 &&
            std::memcmp(&x.gap_extend, &y.gap_extend, sizeof(double)) == 0 && x.max_gap_bp == y.max_gap_bp &&
            x.band == y.band) {
            return static_cast<int>(i);
        }
    }
    if (c.contexts.size() >= 255) return -1;
    c.contexts.push_back(x);
    return static_cast<int>(c.contexts.size() - 1);
}

thread_local CensusSite t_site = CensusSite::other;
thread_local std::uint32_t t_unit = 0;
struct LocalExact {
    std::uint64_t generation = 0;
    std::array<CensusExact, kCensusGroups> cells{};
    bool dirty = false;
};
thread_local LocalExact t_exact;

std::size_t current_group() {
    const CensusStage stage = static_cast<CensusStage>(census().stage.load(std::memory_order_relaxed));
    switch (stage) {
        case CensusStage::union_discovery: return static_cast<std::size_t>(CensusGroup::union_discovery);
        case CensusStage::adaptive:
            return static_cast<std::size_t>(t_site == CensusSite::window || t_site == CensusSite::indexed_unit
                                                ? CensusGroup::adaptive_window : CensusGroup::adaptive_full);
        case CensusStage::dosage:
            return static_cast<std::size_t>(t_site == CensusSite::window ? CensusGroup::dosage_window : CensusGroup::dosage_full);
        default: return static_cast<std::size_t>(CensusGroup::other);
    }
}

LocalExact& local_exact() {
    const std::uint64_t generation = census().generation.load(std::memory_order_relaxed);
    if (t_exact.generation != generation) {
        t_exact = LocalExact{};
        t_exact.generation = generation;
    }
    return t_exact;
}

void flush_local() {
    LocalExact& local = local_exact();
    if (!local.dirty) return;
    Census& c = census();
    const std::lock_guard<std::mutex> lock(c.totals_guard);
    for (std::size_t g = 0; g < kCensusGroups; ++g) {
        CensusExact& to = c.exact[g];
        const CensusExact& from = local.cells[g];
        to.requests += from.requests; to.units += from.units; to.consulted += from.consulted; to.served += from.served;
        to.missed += from.missed; to.missed_not_retained += from.missed_not_retained;
        to.missed_already_present += from.missed_already_present; to.kernels += from.kernels;
        to.feasibility_rejections += from.feasibility_rejections; to.emission_hits += from.emission_hits;
        to.recursions += from.recursions; to.fixed_recursions += from.fixed_recursions;
        to.refused_after_recursion += from.refused_after_recursion; to.accepted += from.accepted; to.rows += from.rows;
    }
    local.cells = {};
    local.dirty = false;
}

void bucket(std::array<std::uint64_t, 6>& histogram, std::uint32_t n) {
    if (n == 0) return;
    const std::size_t i = n == 1 ? 0 : n == 2 ? 1 : n <= 4 ? 2 : n <= 8 ? 3 : n <= 16 ? 4 : 5;
    ++histogram[i];
}

}  // namespace

const char* census_group_name(std::size_t group) {
    static const char* names[kCensusGroups] = {"union_discovery", "adaptive_full_haplotype", "adaptive_window",
                                               "dosage_full_haplotype", "dosage_window", "other"};
    return group < kCensusGroups ? names[group] : "unknown";
}

void begin_verification_census(std::uint32_t period, std::uint64_t budget_bytes) {
    Census& c = census();
    c.active.store(false, std::memory_order_relaxed);
    for (Shard& shard : c.shards) {
        const std::lock_guard<std::mutex> lock(shard.guard);
        std::vector<Entry>().swap(shard.table);
        shard.occupied = 0;
        shard.arena = Arena{};
    }
    {
        const std::lock_guard<std::mutex> lock(c.reads_guard);
        std::vector<ReadSlot>().swap(c.reads);
        c.reads_occupied = 0;
        c.read_arena = Arena{};
    }
    {
        const std::lock_guard<std::mutex> lock(c.contexts_guard);
        c.contexts.clear();
    }
    {
        const std::lock_guard<std::mutex> lock(c.totals_guard);
        c.exact = {};
        c.sampled = {};
    }
    c.period = period;
    c.budget = budget_bytes;
    c.charged.store(0);
    c.peak.store(0);
    c.complete.store(true);
    c.stage.store(0);
    c.next_unit.store(0);
    c.generation.fetch_add(1);
    c.active.store(period > 0, std::memory_order_relaxed);
}

bool verification_census_active() { return census().active.load(std::memory_order_relaxed); }

void set_verification_census_stage(CensusStage stage) { census().stage.store(static_cast<std::uint8_t>(stage)); }
CensusStage verification_census_stage() { return static_cast<CensusStage>(census().stage.load()); }
void set_verification_census_hash_collapse(bool collapse) { census().collapse.store(collapse); }

CensusSiteScope::CensusSiteScope(CensusSite site) : previous_(t_site) { t_site = site; }
CensusSiteScope::~CensusSiteScope() { t_site = previous_; }

CensusUnitScope::CensusUnitScope() : on_(verification_census_active()) {
    if (!on_) return;
    t_unit = census().next_unit.fetch_add(1, std::memory_order_relaxed);
    LocalExact& local = local_exact();
    ++local.cells[current_group()].units;
    local.dirty = true;
}

CensusUnitScope::~CensusUnitScope() {
    if (on_) flush_local();
}

bool census_sampled(std::uint64_t a, std::uint64_t b) {
    const std::uint32_t period = census().period;
    return period != 0 && read_hash_of(a, b) % period == 0;
}

void census_record(const CensusContext& context, std::uint64_t read_hash_a, std::uint64_t read_hash_b, const char* read,
                   std::size_t read_length, const char* window, std::size_t window_length, std::uint32_t fragment,
                   const CensusOutcome& o) {
    Census& c = census();
    if (!c.active.load(std::memory_order_relaxed)) return;
    const std::size_t g = current_group();
    LocalExact& local = local_exact();
    CensusExact& x = local.cells[g];
    local.dirty = true;
    ++x.requests;
    if (o.consulted) {
        ++x.consulted;
        if (o.served) ++x.served;
        else {
            ++x.missed;
            if (o.retained == 0) ++x.missed_not_retained;
            if (o.retained == 2) ++x.missed_already_present;
        }
    }
    if (o.kernel) ++x.kernels;
    if (!o.feasible) ++x.feasibility_rejections;
    if (o.emission_hit) ++x.emission_hits;
    if (o.recursion) ++x.recursions;
    if (o.fixed_recursion) ++x.fixed_recursions;
    if (o.recursion && !o.accepted) ++x.refused_after_recursion;
    if (o.accepted) ++x.accepted;
    x.rows += o.rows;

    const std::uint64_t rh = read_hash_of(read_hash_a, read_hash_b);
    if (rh % c.period != 0) return;
    // Sampled: exact identity by context parameters and bytes.
    std::uint64_t wh = 1469598103934665603ULL;
    for (std::size_t i = 0; i < window_length; ++i) wh = (wh ^ static_cast<unsigned char>(window[i])) * 1099511628211ULL;
    const int ctx = context_id(c, context);
    const long rid = read_id(c, c.collapse.load() ? 0 : rh, read, read_length);
    std::uint64_t hash = splitmix(rh ^ splitmix(wh ^ (static_cast<std::uint64_t>(window_length) << 48) ^
                                                static_cast<std::uint64_t>(ctx < 0 ? 255 : ctx)));
    if (c.collapse.load()) hash = 0;
    Shard& shard = c.shards[hash >> 58];
    const std::lock_guard<std::mutex> lock(shard.guard);
    const auto drop = [&]() {
        const std::lock_guard<std::mutex> totals(c.totals_guard);
        ++c.sampled[g].dropped;
        c.complete.store(false);
    };
    if (ctx < 0 || rid < 0 || window_length > 0xffff) { drop(); return; }
    Entry* e = nullptr;
    if (!shard.table.empty()) {
        const std::size_t mask = shard.table.size() - 1;
        for (std::size_t at = (hash * 0x9e3779b97f4a7c15ULL >> 7) & mask; shard.table[at].window != nullptr; at = (at + 1) & mask) {
            Entry& candidate = shard.table[at];
            if (candidate.hash == hash && candidate.read == static_cast<std::uint32_t>(rid) &&
                candidate.context == static_cast<std::uint8_t>(ctx) && candidate.window_length == window_length &&
                std::memcmp(candidate.window, window, window_length) == 0) {
                e = &candidate;
                break;
            }
        }
    }
    const bool fresh = e == nullptr;
    if (fresh) {
        if (!grow_if_needed(shard.table, shard.occupied, &entry_hash, &entry_empty)) { drop(); return; }
        const char* copy = window_length == 0 ? "" : arena_copy(shard.arena, window, window_length);
        if (copy == nullptr) { drop(); return; }
        const std::size_t mask = shard.table.size() - 1;
        std::size_t at = (hash * 0x9e3779b97f4a7c15ULL >> 7) & mask;
        while (shard.table[at].window != nullptr) at = (at + 1) & mask;
        e = &shard.table[at];
        *e = Entry{};
        e->hash = hash;
        e->window = copy;
        e->read = static_cast<std::uint32_t>(rid);
        e->window_length = static_cast<std::uint16_t>(window_length);
        e->context = static_cast<std::uint8_t>(ctx);
        ++shard.occupied;
    }
    // Per-group sampled totals: accumulated into a local copy and merged under the totals lock, so concurrent shards
    // never write the same counters unguarded.
    CensusSampled d;
    ++d.requests;
    if (fresh) {
        ++d.new_keys;
    } else {
        ++d.repeats;
        if (e->last_group != g) ++d.repeat_from_other_group;
        else if (e->last_unit == t_unit) ++d.repeat_within_unit;
        else ++d.repeat_across_units;
        if (e->last_fragment != fragment) ++d.repeat_other_fragment;
        else if (e->last_unit == t_unit) ++d.repeat_same_fragment_same_unit;
        else ++d.repeat_same_fragment_other_unit;
    }
    if (o.consulted) {
        if (o.served) ++d.served;
        else {
            ++d.missed;
            if (e->missed_before) ++d.repeated_misses;
            if (e->not_retained_before) ++d.repeated_misses_after_not_retained;
            e->missed_before = true;
            if (o.retained == 0) e->not_retained_before = true;
        }
    }
    if (o.kernel) {
        ++d.kernels;
        if (e->kernel_seen) {
            ++d.repeated_kernels;
            if (e->kernel_group != g) ++d.repeated_kernels_cross_group;
        } else {
            e->kernel_seen = true;
            e->kernel_group = static_cast<std::uint8_t>(g);
        }
    }
    if (o.recursion) {
        ++d.recursions;
        if (e->recursion_seen) {
            ++d.repeated_recursions;
            if (e->recursion_group != g) ++d.repeated_recursions_cross_group;
        } else {
            e->recursion_seen = true;
            e->recursion_group = static_cast<std::uint8_t>(g);
        }
    }
    if (o.emission_hit) ++d.emission_hits;
    if (o.accepted) ++d.accepted;
    d.rows += o.rows;
    if (e->requests[g] != 0xffffffffu) ++e->requests[g];
    e->groups = static_cast<std::uint8_t>(e->groups | (1u << g));
    e->last_group = static_cast<std::uint8_t>(g);
    e->last_unit = t_unit;
    e->last_fragment = fragment;
    const std::lock_guard<std::mutex> totals(c.totals_guard);
    CensusSampled& to = c.sampled[g];
    to.requests += d.requests; to.new_keys += d.new_keys; to.repeats += d.repeats;
    to.repeat_from_other_group += d.repeat_from_other_group; to.repeat_within_unit += d.repeat_within_unit;
    to.repeat_across_units += d.repeat_across_units; to.repeat_other_fragment += d.repeat_other_fragment;
    to.repeat_same_fragment_same_unit += d.repeat_same_fragment_same_unit;
    to.repeat_same_fragment_other_unit += d.repeat_same_fragment_other_unit;
    to.served += d.served; to.missed += d.missed; to.repeated_misses += d.repeated_misses;
    to.repeated_misses_after_not_retained += d.repeated_misses_after_not_retained;
    to.kernels += d.kernels; to.repeated_kernels += d.repeated_kernels;
    to.repeated_kernels_cross_group += d.repeated_kernels_cross_group;
    to.recursions += d.recursions; to.repeated_recursions += d.repeated_recursions;
    to.repeated_recursions_cross_group += d.repeated_recursions_cross_group;
    to.emission_hits += d.emission_hits; to.accepted += d.accepted; to.rows += d.rows;
}

VerificationCensusReport end_verification_census() {
    Census& c = census();
    VerificationCensusReport r;
    r.active = c.active.load();
    flush_local();
    c.active.store(false, std::memory_order_relaxed);
    r.period = c.period;
    r.budget_bytes = c.budget;
    r.peak_bytes = c.peak.load();
    r.complete = c.complete.load();
    std::map<std::uint32_t, std::uint64_t> masks;
    {
        const std::lock_guard<std::mutex> lock(c.totals_guard);
        r.exact = c.exact;
        r.sampled = c.sampled;
    }
    for (Shard& shard : c.shards) {
        const std::lock_guard<std::mutex> lock(shard.guard);
        for (const Entry& e : shard.table) {
            if (e.window == nullptr) continue;
            ++r.distinct_keys;
            ++masks[e.groups];
            for (std::size_t g = 0; g < kCensusGroups; ++g) {
                if (e.requests[g] == 0) continue;
                ++r.sampled[g].keys_seen;
                bucket(r.sampled[g].requests_per_key, e.requests[g]);
            }
        }
        std::vector<Entry>().swap(shard.table);
        shard.occupied = 0;
        shard.arena = Arena{};
    }
    {
        const std::lock_guard<std::mutex> lock(c.reads_guard);
        r.sampled_read_contents = c.reads_occupied;
        std::vector<ReadSlot>().swap(c.reads);
        c.reads_occupied = 0;
        c.read_arena = Arena{};
    }
    {
        const std::lock_guard<std::mutex> lock(c.contexts_guard);
        r.contexts = c.contexts.size();
        c.contexts.clear();
    }
    r.keys_by_group_mask.assign(masks.begin(), masks.end());
    c.charged.store(0);
    c.generation.fetch_add(1);
    return r;
}

std::string verification_census_table(const VerificationCensusReport& r) {
    std::ostringstream out;
    out << "group\tmetric\tvalue\n";
    out << "all\tcensus_period\t" << r.period << "\nall\tcensus_budget_bytes\t" << r.budget_bytes
        << "\nall\tcensus_peak_bytes\t" << r.peak_bytes << "\nall\tcensus_coverage\t" << (r.complete ? "complete" : "incomplete")
        << "\nall\tsampled_read_contents\t" << r.sampled_read_contents << "\nall\tsampled_distinct_keys\t" << r.distinct_keys
        << "\nall\tscoring_contexts\t" << r.contexts << '\n';
    for (const auto& item : r.keys_by_group_mask) {
        std::string groups;
        for (std::size_t g = 0; g < kCensusGroups; ++g) {
            if ((item.first >> g) & 1u) groups += (groups.empty() ? "" : "+") + std::string(census_group_name(g));
        }
        out << "all\tsampled_keys_requested_by:" << groups << '\t' << item.second << '\n';
    }
    for (std::size_t g = 0; g < kCensusGroups; ++g) {
        const CensusExact& x = r.exact[g];
        const CensusSampled& s = r.sampled[g];
        if (x.requests == 0 && x.units == 0 && s.requests == 0) continue;
        const std::string n = census_group_name(g);
        const auto row = [&](const char* metric, std::uint64_t value) { out << n << '\t' << metric << '\t' << value << '\n'; };
        row("exact_units", x.units);
        row("exact_requests", x.requests);
        row("exact_decision_consulted", x.consulted);
        row("exact_decision_served", x.served);
        row("exact_decision_missed", x.missed);
        row("exact_decision_missed_not_retained", x.missed_not_retained);
        row("exact_decision_missed_already_present", x.missed_already_present);
        row("exact_kernels", x.kernels);
        row("exact_feasibility_rejections", x.feasibility_rejections);
        row("exact_emission_hits", x.emission_hits);
        row("exact_recursions", x.recursions);
        row("exact_fixed_geometry_recursions", x.fixed_recursions);
        row("exact_refused_after_recursion", x.refused_after_recursion);
        row("exact_accepted", x.accepted);
        row("exact_rows", x.rows);
        row("sampled_requests", s.requests);
        row("sampled_keys_seen", s.keys_seen);
        row("sampled_new_keys", s.new_keys);
        row("sampled_repeats", s.repeats);
        row("sampled_repeat_from_other_group", s.repeat_from_other_group);
        row("sampled_repeat_within_unit", s.repeat_within_unit);
        row("sampled_repeat_across_units", s.repeat_across_units);
        row("sampled_repeat_other_fragment", s.repeat_other_fragment);
        row("sampled_repeat_same_fragment_same_unit", s.repeat_same_fragment_same_unit);
        row("sampled_repeat_same_fragment_other_unit", s.repeat_same_fragment_other_unit);
        row("sampled_decision_served", s.served);
        row("sampled_decision_missed", s.missed);
        row("sampled_repeated_misses", s.repeated_misses);
        row("sampled_repeated_misses_after_not_retained", s.repeated_misses_after_not_retained);
        row("sampled_kernels", s.kernels);
        row("sampled_repeated_kernels", s.repeated_kernels);
        row("sampled_repeated_kernels_cross_group", s.repeated_kernels_cross_group);
        row("sampled_recursions", s.recursions);
        row("sampled_repeated_recursions", s.repeated_recursions);
        row("sampled_repeated_recursions_cross_group", s.repeated_recursions_cross_group);
        row("sampled_emission_hits", s.emission_hits);
        row("sampled_accepted", s.accepted);
        row("sampled_rows", s.rows);
        row("sampled_dropped_for_budget", s.dropped);
        static const char* bins[6] = {"1", "2", "3-4", "5-8", "9-16", "17+"};
        for (std::size_t i = 0; i < 6; ++i) {
            out << n << "\tsampled_keys_with_requests_" << bins[i] << '\t' << s.requests_per_key[i] << '\n';
        }
    }
    return out.str();
}

}  // namespace panvar
