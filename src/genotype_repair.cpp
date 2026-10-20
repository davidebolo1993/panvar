#include "panvar/genotype_repair.hpp"

#include "panvar/genotype_command.hpp"
#include "panvar/ls_hmm.hpp"
#include "panvar/mosaic_spelling.hpp"

#include <sys/resource.h>
#if defined(__APPLE__)
#include <mach/mach.h>
#include <malloc/malloc.h>
#endif
#include <unistd.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <optional>
#include <random>
#include <sstream>
#include <stdexcept>

namespace panvar {

// ---- process resources ------------------------------------------------------------------------

std::uint64_t current_rss_bytes() {
#if defined(__APPLE__)
    mach_task_basic_info_data_t info{};
    mach_msg_type_number_t count = MACH_TASK_BASIC_INFO_COUNT;
    if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO, reinterpret_cast<task_info_t>(&info),
                  &count) != KERN_SUCCESS) {
        return 0;
    }
    return info.resident_size;
#else
    std::ifstream statm("/proc/self/statm");
    std::uint64_t size = 0;
    std::uint64_t resident = 0;
    statm >> size >> resident;
    return resident * static_cast<std::uint64_t>(sysconf(_SC_PAGESIZE));
#endif
}

std::optional<std::uint64_t> allocator_bytes_in_use() {
#if defined(__APPLE__)
    malloc_statistics_t stats{};
    malloc_zone_statistics(nullptr, &stats);
    return static_cast<std::uint64_t>(stats.size_in_use);
#else
    return std::nullopt;
#endif
}

std::uint64_t release_free_pages() {
#if defined(__APPLE__)
    return static_cast<std::uint64_t>(malloc_zone_pressure_relief(nullptr, 0));
#else
    return 0;
#endif
}

FootprintSample physical_footprint() {
    FootprintSample sample;
#if defined(__APPLE__)
    task_vm_info_data_t info{};
    mach_msg_type_number_t count = TASK_VM_INFO_COUNT;
    if (task_info(mach_task_self(), TASK_VM_INFO, reinterpret_cast<task_info_t>(&info), &count) !=
        KERN_SUCCESS) {
        return sample;
    }
    if (count >= TASK_VM_INFO_REV1_COUNT) sample.current = info.phys_footprint;
    if (count >= TASK_VM_INFO_REV3_COUNT && info.ledger_phys_footprint_peak > 0) {
        sample.peak = static_cast<std::uint64_t>(info.ledger_phys_footprint_peak);
    }
#endif
    return sample;
}

void atomic_write(const std::string& path, const std::string& content) {
    const std::string staging = path + ".writing";
    {
        std::ofstream file(staging, std::ios::binary | std::ios::trunc);
        if (!file) throw std::runtime_error("cannot write " + staging);
        file << content;
        file.flush();
        if (!file) throw std::runtime_error("short write to " + staging);
    }
    std::filesystem::rename(staging, path);
}

double peak_rss_gib() {
    rusage usage{};
    getrusage(RUSAGE_SELF, &usage);
#if defined(__APPLE__)
    return static_cast<double>(usage.ru_maxrss) / (1024.0 * 1024.0 * 1024.0);   // Darwin reports bytes
#else
    return static_cast<double>(usage.ru_maxrss) / (1024.0 * 1024.0);            // Linux reports KiB
#endif
}

namespace {

// FIPS 180-4 SHA-256, streaming. Provenance of inputs only; nothing statistical depends on it.
class Sha256 {
public:
    void update(const unsigned char* data, std::size_t size) {
        for (std::size_t i = 0; i < size; ++i) {
            block_[used_++] = data[i];
            if (used_ == 64) {
                compress();
                used_ = 0;
            }
        }
        bits_ += static_cast<std::uint64_t>(size) * 8;
    }
    std::string hex() {
        const std::uint64_t bits = bits_;
        const unsigned char pad = 0x80;
        update(&pad, 1);
        const unsigned char zero = 0;
        while (used_ != 56) update(&zero, 1);
        unsigned char length[8];
        for (int i = 0; i < 8; ++i) length[i] = static_cast<unsigned char>(bits >> (56 - 8 * i));
        update(length, 8);
        std::ostringstream out;
        for (const std::uint32_t word : state_) out << std::hex << std::setw(8) << std::setfill('0') << word;
        return out.str();
    }

private:
    static std::uint32_t rotr(std::uint32_t x, int n) { return (x >> n) | (x << (32 - n)); }
    void compress() {
        static const std::uint32_t k[64] = {
            0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
            0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
            0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
            0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
            0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
            0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
            0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
            0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2};
        std::uint32_t w[64];
        for (int i = 0; i < 16; ++i) {
            w[i] = (static_cast<std::uint32_t>(block_[4 * i]) << 24) |
                   (static_cast<std::uint32_t>(block_[4 * i + 1]) << 16) |
                   (static_cast<std::uint32_t>(block_[4 * i + 2]) << 8) | static_cast<std::uint32_t>(block_[4 * i + 3]);
        }
        for (int i = 16; i < 64; ++i) {
            const std::uint32_t s0 = rotr(w[i - 15], 7) ^ rotr(w[i - 15], 18) ^ (w[i - 15] >> 3);
            const std::uint32_t s1 = rotr(w[i - 2], 17) ^ rotr(w[i - 2], 19) ^ (w[i - 2] >> 10);
            w[i] = w[i - 16] + s0 + w[i - 7] + s1;
        }
        std::uint32_t a = state_[0], b = state_[1], c = state_[2], d = state_[3];
        std::uint32_t e = state_[4], f = state_[5], g = state_[6], h = state_[7];
        for (int i = 0; i < 64; ++i) {
            const std::uint32_t s1 = rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25);
            const std::uint32_t ch = (e & f) ^ (~e & g);
            const std::uint32_t t1 = h + s1 + ch + k[i] + w[i];
            const std::uint32_t s0 = rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22);
            const std::uint32_t maj = (a & b) ^ (a & c) ^ (b & c);
            const std::uint32_t t2 = s0 + maj;
            h = g; g = f; f = e; e = d + t1; d = c; c = b; b = a; a = t1 + t2;
        }
        state_[0] += a; state_[1] += b; state_[2] += c; state_[3] += d;
        state_[4] += e; state_[5] += f; state_[6] += g; state_[7] += h;
    }
    std::array<std::uint32_t, 8> state_ = {0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
                                           0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19};
    unsigned char block_[64] = {};
    std::size_t used_ = 0;
    std::uint64_t bits_ = 0;
};

bool same_score(const DiplotypeScore& x, const DiplotypeScore& y) {
    return std::memcmp(&x.log_likelihood, &y.log_likelihood, sizeof(double)) == 0 &&
           std::memcmp(&x.fragment_term, &y.fragment_term, sizeof(double)) == 0 &&
           std::memcmp(&x.exposure_term, &y.exposure_term, sizeof(double)) == 0 &&
           std::memcmp(&x.absolute_term_sum, &y.absolute_term_sum, sizeof(double)) == 0 &&
           x.fragments_with_mass == y.fragments_with_mass;
}

std::string join(const std::vector<std::size_t>& v) {
    std::string text;
    for (std::size_t i = 0; i < v.size(); ++i) {
        if (i != 0) text += ',';
        text += std::to_string(v[i]);
    }
    return text;
}

std::string margin(std::uint64_t cap_value, std::uint64_t used) {
    if (cap_value == 0) return "none";
    return std::to_string(static_cast<std::int64_t>(cap_value) - static_cast<std::int64_t>(used));
}

} // namespace

std::string sha256_file(const std::string& path) {
    std::ifstream file(path, std::ios::binary);
    if (!file) throw std::runtime_error("cannot open for hashing: " + path);
    Sha256 hash;
    std::vector<char> buffer(1 << 20);
    while (file) {
        file.read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
        const std::streamsize got = file.gcount();
        if (got > 0) hash.update(reinterpret_cast<const unsigned char*>(buffer.data()), static_cast<std::size_t>(got));
    }
    return hash.hex();
}

std::string sha256_bytes(const std::string& bytes) {
    Sha256 hash;
    hash.update(reinterpret_cast<const unsigned char*>(bytes.data()), bytes.size());
    return hash.hex();
}

struct Sha256Stream::State {
    Sha256 hash;
};
Sha256Stream::Sha256Stream() : state_(std::make_unique<State>()) {}
Sha256Stream::~Sha256Stream() = default;
void Sha256Stream::update(const void* data, std::size_t size) {
    state_->hash.update(static_cast<const unsigned char*>(data), size);
}
std::string Sha256Stream::hex() { return state_->hash.hex(); }

// ---- the search --------------------------------------------------------------------------------

DiplotypeIdentity diplotype_identity(const std::vector<std::size_t>& a, const std::vector<std::size_t>& b) {
    return a <= b ? std::make_pair(a, b) : std::make_pair(b, a);
}

std::uint64_t repair_neighbour_bound(const PreparedGenotypeModel& model, std::size_t local_pairs) {
    std::uint64_t bound = 0;
    for (std::size_t b = 0; b < model.blocks().size(); ++b) {
        const std::uint64_t a = model.blocks()[b].alleles.size();
        bound += std::min<std::uint64_t>(local_pairs, a * (a + 1) / 2) * 2;
    }
    return bound;
}

CoordinateRepairSearch::CoordinateRepairSearch(const PreparedGenotypeModel& model,
                                               const PreparedProposalEmissions& emissions,
                                               const std::vector<ReadPair>& fragments,
                                               const LibraryModel& library,
                                               const FragmentScoringOptions& scoring,
                                               const RepairSearchOptions& options)
    : model_(model), fragments_(fragments), library_(library), scoring_(scoring), options_(options) {
    const std::size_t blocks = model.blocks().size();
    if (blocks < 2) throw std::invalid_argument("the repair search needs at least two blocks");
    if (emissions.blocks.size() != blocks) throw std::invalid_argument("one proposal emission per block");
    if (!options.switch_vector.empty()) {   // M4-B: the run's shared transition vector
        if (options.switch_vector.size() + 1 != blocks) throw std::invalid_argument("the repair's transition vector does not match the block count");
        switches_ = options.switch_vector;
    } else {
        switches_.assign(blocks - 1, options.switch_probability);
    }
    // The per-block local pairs, ranked once. Deterministic: emission descending, then indices.
    local_.resize(blocks);
    for (std::size_t b = 0; b < blocks; ++b) {
        const BlockProposalEmissions& emission = emissions.blocks[b];
        std::vector<std::pair<double, std::pair<std::size_t, std::size_t>>> ranked;
        for (std::size_t x = 0; x < emission.allele_count; ++x) {
            for (std::size_t y = x; y < emission.allele_count; ++y) ranked.push_back({emission.at(x, y), {x, y}});
        }
        std::sort(ranked.begin(), ranked.end(), [](const auto& lhs, const auto& rhs) {
            if (lhs.first != rhs.first) return lhs.first > rhs.first;
            return lhs.second < rhs.second;
        });
        const std::size_t take = std::min(options.local_pairs, ranked.size());
        for (std::size_t r = 0; r < take; ++r) local_[b].push_back(ranked[r].second);
    }
    // M9-G: junction-eligible blocks (every pair emission tied, more pairs than the limit), only when the experiment is on
    junction_block_.assign(blocks, 0);
    if (options.junction_reads) {
        for (std::size_t b = 0; b < blocks; ++b) {
            const BlockProposalEmissions& emission = emissions.blocks[b];
            std::vector<double> all;
            for (std::size_t x = 0; x < emission.allele_count; ++x) {
                for (std::size_t y = x; y < emission.allele_count; ++y) all.push_back(emission.at(x, y));
            }
            if (junction_block_eligible(all, options.local_pairs)) {
                junction_block_[b] = 1;
                ++junction_report_.eligible_blocks;
            }
        }
        junction_report_.block_eligible.assign(junction_block_.begin(), junction_block_.end());   // M11-C: reported, never read by the search
    }
    neighbour_bound_ = repair_neighbour_bound(model, options.local_pairs);
    if (options_.mode != RepairScoringMode::OneShot && !options_.external_context) new_context();
    memory_on_ = context_ != nullptr && options_.memory_log;
    if (memory_on_) {
        memory_text_ += "kind\tindex\tpoint\trss_current\trss_peak\trelief_freed\trss_after_relief\t"
                        "slots_total\tslots_live\tslots_released\tslot_table\tsequence_bytes\t"
                        "sequence_capacity\tsegment_bytes\tmass_bytes\tmass_capacity\tcensus_count\t"
                        "census_bytes\tcensus_capacity\torigin_cache_entries\torigin_cache_payload\t"
                        "origin_cache_capacity\tmerge_bytes\tmerge_capacity\twindow_scratch\t"
                        "origin_scratch\temission_table\temission_arena\tprepared\tscoring_result\t"
                        "logical_total\tcapacity_total\trss_minus_capacity\tphys_footprint\t"
                        "phys_footprint_peak\torigin_demand\tcensus_slots\tcensus_roles\n";
    }
    hold_remaining_ = options_.hold_cycles;
}

CoordinateRepairSearch::~CoordinateRepairSearch() = default;

void CoordinateRepairSearch::new_context() {
    context_.reset();
    context_ = std::make_unique<PreparedScoringContext>(fragments_, library_, scoring_);
    if (options_.origin_cache) {
        context_->enable_origin_cache_path(true);
        context_->set_origin_cache_limits(options_.origin_cache_bytes, options_.origin_cache_records);
    }
    if (options_.fine_window_bytes != 0) {
        context_->enable_fine_window_discovery(options_.fine_window_bytes, options_.fine_window_served_cap);
    }
    transient_opened_ = false;
    seeding_origins_ = 0;   // seeding work is counted within one context
}

void CoordinateRepairSearch::open_transient_round() {
    if (transient_opened_ || !context()) return;
    context()->begin_transient_round();
    transient_opened_ = true;
}

void CoordinateRepairSearch::attach_context(PreparedScoringContext* shared) {
    if (!options_.external_context || shared == nullptr) {
        throw std::logic_error("attach_context needs options.external_context and a live context");
    }
    shared_ = shared;
    transient_opened_ = true;   // the owner opened it
}

std::pair<std::size_t, std::size_t> CoordinateRepairSearch::seed_neighbour_by_delta(
    const std::vector<std::size_t>& first, const std::vector<std::size_t>& second, BlockDeltaWork* work) {
    PreparedScoringContext* context = this->context();
    if (context == nullptr || !have_inc_slots_) throw std::logic_error("seeding needs a seeded delta search");
    const std::size_t blocks = model_.blocks().size();
    const HaplotypeSequence incumbent_a = spell_mosaic(model_, MosaicPath{cur_first_});
    const HaplotypeSequence incumbent_b = spell_mosaic(model_, MosaicPath{cur_second_});
    // As a round resolves a local path: an incumbent haplotype is its slot; otherwise a delta from the
    // same side's incumbent at the one block where they differ, census retained (a promotion).
    const auto resolve = [&](const std::vector<std::size_t>& path, int side) {
        if (path == cur_first_) return inc_slot_a_;
        if (path == cur_second_) return inc_slot_b_;
        const std::vector<std::size_t>& base_path = side == 0 ? cur_first_ : cur_second_;
        const HaplotypeSequence& base = side == 0 ? incumbent_a : incumbent_b;
        std::size_t block = blocks;
        for (std::size_t b = 0; b < blocks; ++b) {
            if (path[b] == base_path[b]) continue;
            if (block != blocks) throw std::invalid_argument("a seeded neighbour differs from its side at two blocks");
            block = b;
        }
        const SequenceInterval span = base.block_intervals[block];
        BlockDeltaWork local;
        const std::size_t slot = context->add_repaired_haplotype(
            side == 0 ? inc_slot_a_ : inc_slot_b_, span.begin, span.end,
            model_.blocks()[block].alleles[path[block]].sequence, &local, true);
        if (work != nullptr) {
            work->origins_reused += local.origins_reused;
            work->origins_recomputed += local.origins_recomputed;
            work->physical_origins_constructed += local.physical_origins_constructed;
            work->seconds += local.seconds;
        }
        return slot;
    };
    const std::size_t a = resolve(first, 0);
    const std::size_t b = second == first ? a : resolve(second, 1);
    return {a, b};
}

void CoordinateRepairSearch::start_from_slots(const std::vector<std::size_t>& first,
                                              const std::vector<std::size_t>& second, double total,
                                              std::size_t slot_first, std::size_t slot_second) {
    PreparedScoringContext* context = this->context();
    if (context == nullptr || options_.mode != RepairScoringMode::Delta) {
        throw std::logic_error("start_from_slots needs a delta search with a context");
    }
    if (!context->holds_census(slot_first) || !context->holds_census(slot_second)) {
        throw std::invalid_argument("a seeded start must hold a census on both haplotypes");
    }
    start(first, second, total);
    inc_slot_a_ = slot_first;
    inc_slot_b_ = slot_second;
    role_of_[inc_slot_a_] = inc_slot_a_ == inc_slot_b_ ? "seed_hap1+2" : "seed_hap1";
    if (inc_slot_b_ != inc_slot_a_) role_of_[inc_slot_b_] = "seed_hap2";
    have_inc_slots_ = true;
}

void CoordinateRepairSearch::start(const std::vector<std::size_t>& first, const std::vector<std::size_t>& second,
                                   double total) {
    const std::size_t blocks = model_.blocks().size();
    if (first.size() != blocks || second.size() != blocks) {
        throw std::invalid_argument("a repair start needs one allele per block on both haplotypes");
    }
    cur_first_ = first;
    cur_second_ = second;
    cur_total_ = total;
    visited_.clear();
    visited_.insert(diplotype_identity(first, second));
    rounds_.clear();
    stop_reason_ = "round_cap";
    have_inc_slots_ = false;
    role_of_.clear();
    failed_round_ = 0;
}

void CoordinateRepairSearch::add_visited(const std::vector<std::size_t>& first, const std::vector<std::size_t>& second) {
    visited_.insert(diplotype_identity(first, second));
}

bool same_decision(const RoundVerdict& a, const RoundVerdict& b) {
    return a.complete == b.complete && a.unique_best == b.unique_best &&
           a.incumbent_in_top_set == b.incumbent_in_top_set && a.moved == b.moved && a.stop == b.stop &&
           a.winner == b.winner;
}

const std::vector<std::pair<std::size_t, std::size_t>>& CoordinateRepairSearch::local_for(
    std::size_t b, const std::vector<std::size_t>& first, const std::vector<std::size_t>& second) {
    if (b >= junction_block_.size() || junction_block_[b] == 0) return local_[b];
    const JunctionProposalCache::Key key = JunctionProposalCache::key(b, first, second);
    if (const JunctionProposal* hit = junction_cache_.find(key)) {
        ++junction_report_.cache_hits;
        return hit->pairs;
    }
    const auto mark = std::chrono::steady_clock::now();
    const std::size_t blocks = model_.blocks().size();
    std::vector<std::string> alleles;
    alleles.reserve(model_.blocks()[b].alleles.size());
    for (const auto& allele : model_.blocks()[b].alleles) alleles.push_back(allele.sequence);
    JunctionContext context;
    const std::vector<std::size_t>* haps[2] = {&first, &second};
    for (int h = 0; h < 2; ++h) {
        for (std::size_t i = 0; i < b; ++i) context.prefix[h] += model_.blocks()[i].alleles[(*haps[h])[i]].sequence;
        for (std::size_t i = b + 1; i < blocks; ++i) context.suffix[h] += model_.blocks()[i].alleles[(*haps[h])[i]].sequence;
        context.incumbent[h] = (*haps[h])[b];
    }
    JunctionProposal proposal = junction_local_pairs(alleles, context, *options_.junction_reads, options_.junction_rate,
                                                     options_.junction_background, options_.local_pairs);
    ++junction_report_.invocations;
    junction_report_.cutoff_tie_crossings += proposal.cutoff_tie_crosses ? 1 : 0;
    junction_report_.seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - mark).count();
    return junction_cache_.insert(key, std::move(proposal)).pairs;
}

bool CoordinateRepairSearch::converged() const {
    if (stop_reason_ == "no_improvement" || stop_reason_ == "no_unvisited_neighbour") return true;
    // A round tie with the incumbent inside the top set: nothing improves on it beyond error.
    if (stop_reason_ == "rescore_incomplete" && !rounds_.empty()) return rounds_.back().incumbent_in_top_set;
    return false;
}

void CoordinateRepairSearch::snapshot(const char* kind, std::size_t index, const char* point, bool result_alive,
                                      bool relieve) {
    if (!memory_on_) return;
    const ContextMemoryAccount a = context_->memory_account();
    const FootprintSample footprint = physical_footprint();
    std::string slots_text;
    std::string roles_text;
    for (const std::size_t slot : context_->census_slots()) {
        if (!slots_text.empty()) {
            slots_text += ',';
            roles_text += ',';
        }
        slots_text += std::to_string(slot);
        const auto role = role_of_.find(slot);
        roles_text += role == role_of_.end() ? std::string("unknown") : role->second;
        if (slot == inc_slot_a_ || slot == inc_slot_b_) roles_text += "*";
    }
    const std::uint64_t rss = current_rss_bytes();
    const std::uint64_t peak = static_cast<std::uint64_t>(peak_rss_gib() * 1073741824.0);
    std::uint64_t freed = 0;
    std::uint64_t after = rss;
    if (relieve) {
        freed = release_free_pages();
        after = current_rss_bytes();
    }
    // The flattened result belongs to the caller: count it only while the caller still holds it.
    const std::uint64_t result = result_alive ? a.last_scoring_result_bytes : 0;
    const std::uint64_t logical = a.logical_total - a.last_scoring_result_bytes + result;
    const std::uint64_t capacity = a.capacity_total - a.last_scoring_result_bytes + result;
    std::ostringstream memory_log;
    memory_log << kind << '\t' << index << '\t' << point << '\t' << rss << '\t' << peak << '\t' << freed << '\t'
               << after << '\t' << a.slots_total << '\t' << a.slots_live << '\t' << a.slots_released << '\t'
               << a.slot_table_bytes << '\t' << a.sequence_bytes << '\t' << a.sequence_capacity << '\t'
               << a.segment_bytes << '\t' << a.mass_column_bytes << '\t' << a.mass_column_capacity << '\t'
               << a.census_count << '\t' << a.census_bytes << '\t' << a.census_capacity << '\t'
               << a.origin_cache_entries << '\t' << a.origin_cache_payload << '\t' << a.origin_cache_capacity
               << '\t' << a.merge_working_bytes << '\t' << a.merge_working_capacity << '\t'
               << a.window_scratch_capacity << '\t' << a.origin_scratch_capacity << '\t' << a.emission_cache_table
               << '\t' << a.emission_cache_arena << '\t' << a.prepared_capacity << '\t' << result << '\t' << logical
               << '\t' << capacity << '\t'
               << (static_cast<std::int64_t>(rss) - static_cast<std::int64_t>(capacity)) << '\t'
               << footprint.current << '\t' << footprint.peak << '\t' << a.origin_demand_bytes << '\t'
               << (slots_text.empty() ? "-" : slots_text) << '\t' << (roles_text.empty() ? "-" : roles_text)
               << '\n';
    memory_text_ += memory_log.str();
}

void CoordinateRepairSearch::live_probe(std::size_t round, const char* point, std::size_t done, std::size_t total) {
    if (options_.memory_live_path.empty()) return;
    PreparedScoringContext* context = this->context();
    if (context == nullptr) return;
    std::ofstream out(options_.memory_live_path, std::ios::app);
    if (!live_header_written_) {
        out << "round\tpoint\trequests_done\trequests_total\tphys_footprint\tphys_footprint_peak\trss_current\trss_peak\t"
               "slots_total\tslots_live\tslots_released\tslot_table\tsequence_bytes\tsequence_capacity\tsegment_bytes\t"
               "mass_bytes\tmass_capacity\tcensus_count\tcensus_bytes\tcensus_capacity\torigin_cache_entries\t"
               "origin_cache_payload\torigin_cache_capacity\tmerge_bytes\tmerge_capacity\torigin_demand\twindow_scratch\t"
               "origin_scratch\temission_table\temission_arena\temission_entries\temission_entry_metadata\t"
               "emission_content_arena\tprepared\tlast_scoring_result\tlogical_total\tcapacity_total\t"
               "footprint_minus_capacity\thaplotypes_created\tunique_haplotypes\tcum_physical_origins\tcum_state_classes\t"
               "cum_piece_hits\tallocator_in_use\tpartitions\tpartition_scratch_bytes\tpartition_scratch_capacity\t"
               "partition_piece_index_bytes\tpartition_piece_index_capacity\tcapacity_total_extended\t"
               "footprint_minus_capacity_extended\treuse_only_admission\tadmission_restricted_plans\t"
               "admission_singleton_windows\tadmission_eligible_windows\tadmission_misses_while_restricted\t"
               "admission_singleton_misses_not_retained\tadmission_retained_while_restricted\t"
               "admission_hits_while_restricted\tadmission_group_reuses_while_restricted\n";
        live_header_written_ = true;
    }
    const ContextMemoryAccount a = context->memory_account();
    const FootprintSample fp = physical_footprint();
    const EmissionCacheCounters cache = read_emission_cache_counters();
    const FragmentScoringCounters& cum = context->cumulative_counters();
    out << round << '\t' << point << '\t' << done << '\t' << total << '\t' << fp.current << '\t' << fp.peak << '\t'
        << current_rss_bytes() << '\t' << static_cast<std::uint64_t>(peak_rss_gib() * 1073741824.0) << '\t' << a.slots_total
        << '\t' << a.slots_live << '\t' << a.slots_released << '\t' << a.slot_table_bytes << '\t' << a.sequence_bytes << '\t'
        << a.sequence_capacity << '\t' << a.segment_bytes << '\t' << a.mass_column_bytes << '\t' << a.mass_column_capacity
        << '\t' << a.census_count << '\t' << a.census_bytes << '\t' << a.census_capacity << '\t' << a.origin_cache_entries
        << '\t' << a.origin_cache_payload << '\t' << a.origin_cache_capacity << '\t' << a.merge_working_bytes << '\t'
        << a.merge_working_capacity << '\t' << a.origin_demand_bytes << '\t' << a.window_scratch_capacity << '\t'
        << a.origin_scratch_capacity << '\t' << a.emission_cache_table << '\t' << a.emission_cache_arena << '\t' << cache.entries
        << '\t' << cache.entry_metadata_bytes << '\t' << cache.content_arena_bytes << '\t' << a.prepared_capacity << '\t'
        << a.last_scoring_result_bytes << '\t' << a.logical_total << '\t' << a.capacity_total << '\t'
        << (static_cast<std::int64_t>(fp.current) - static_cast<std::int64_t>(a.capacity_total)) << '\t'
        << context->haplotypes_created() << '\t' << context->unique_haplotypes() << '\t' << cum.physical_origins << '\t'
        << cum.statistical_state_classes << '\t' << cum.piece_hits << '\t';
    // S2-I: the allocator's own in-use total (a different quantity from the footprint), and the workers' storage
    const std::optional<std::uint64_t> in_use = allocator_bytes_in_use();
    out << (in_use ? std::to_string(*in_use) : std::string("unavailable")) << '\t' << a.partitions << '\t'
        << a.partition_scratch_bytes << '\t' << a.partition_scratch_capacity << '\t' << a.partition_piece_index_bytes << '\t'
        << a.partition_piece_index_capacity << '\t' << a.capacity_total_extended << '\t'
        << (static_cast<std::int64_t>(fp.current) - static_cast<std::int64_t>(a.capacity_total_extended)) << '\t';
    const OriginAdmissionReport ad = context->origin_admission_report();   // S2-J
    out << (context->reuse_only_admission() ? 1 : 0) << '\t' << ad.restricted_plans << '\t' << ad.singleton_windows << '\t'
        << ad.eligible_windows << '\t' << ad.misses_while_restricted << '\t' << ad.singleton_misses_not_retained << '\t'
        << ad.retained_while_restricted << '\t' << ad.hits_while_restricted << '\t' << ad.group_reuses_while_restricted << '\n';
}

void CoordinateRepairSearch::write_seed_report(const IncrementalBatchReport& r) {
    std::string file = options_.memory_live_path;
    const std::string live = "refinement_memory_live.tsv";
    if (file.size() >= live.size() && file.compare(file.size() - live.size(), live.size(), live) == 0) {
        file.replace(file.size() - live.size(), live.size(), "refinement_seed_report.tsv");
    } else {
        file += ".seed_report.tsv";
    }
    const ScorerMemoryAccounts& m = r.memory;
    std::ofstream out(file, std::ios::trunc);
    out << "key\tvalue\nschema\trefinement_seed_report_v1\nnote\tDIAGNOSTIC ONLY: the seed placement's own batch report "
           "(add_haplotypes). The *_peak values are per-worker transient peaks summed over workers at each worker's own peak: "
           "not simultaneously live, never part of a live total\n"
        << "haplotypes_offered\t" << r.haplotypes_offered << "\nhaplotypes_new\t" << r.haplotypes_new << "\nthreads_used_placement\t"
        << r.threads_used_placement << "\nseconds_placement_and_mass\t" << r.seconds_placement_and_mass
        << "\nmass_matrix_bytes\t" << m.mass_matrix_bytes << "\nhaplotype_sequence_bytes\t" << m.haplotype_sequence_bytes
        << "\nform_sequence_bytes\t" << m.form_sequence_bytes << "\npiece_index_bytes\t" << m.piece_index_bytes
        << "\nemission_cache_table_bytes\t" << m.emission_cache_table_bytes << "\nemission_cache_arena_bytes\t"
        << m.emission_cache_arena_bytes << "\nplacement_rows_peak\t" << m.placement_rows_peak
        << "\nplacement_capacity_bytes_peak\t" << m.placement_capacity_bytes_peak << "\nproposed_starts_peak\t"
        << m.proposed_starts_peak << "\nproposed_starts_capacity_bytes_peak\t" << m.proposed_starts_capacity_bytes_peak
        << "\nphysical_origins\t" << r.counters.physical_origins << "\nstatistical_state_classes\t"
        << r.counters.statistical_state_classes << "\nemission_entries\t" << r.emission_cache.entries << '\n';
}

bool CoordinateRepairSearch::run(const RepairHooks& hooks) {
    const std::size_t blocks = model_.blocks().size();
    const bool use_delta = options_.mode == RepairScoringMode::Delta;
    const bool planned_cache = options_.planned_cache;
    const std::uint64_t request_order_seed = options_.request_order_seed;
    PreparedScoringContext* context = this->context();
    const std::size_t first_round = options_.first_round;
    const std::size_t final_round = last_round();
    const auto key_of = [](const std::vector<std::size_t>& a, const std::vector<std::size_t>& b) {
        return diplotype_identity(a, b);
    };

    // I2-C: measurement only -- disjoint wall laps, accumulated into stage_seconds_
    using StageClock = std::chrono::steady_clock;
    const auto since = [](StageClock::time_point t) { return std::chrono::duration<double>(StageClock::now() - t).count(); };
    RepairStageSeconds& lap = stage_seconds_;
    // S2-H (diagnostic): the context's batch-repair safe points feed the live file; removed when run() returns
    std::size_t probe_round = first_round;
    struct ProbeGuard {
        PreparedScoringContext* context = nullptr;
        ~ProbeGuard() {
            if (context != nullptr) context->set_progress_probe({});
        }
    } probe_guard;
    if (!options_.memory_live_path.empty() && context != nullptr) {
        context->set_progress_probe([this, &probe_round](const char* point, std::size_t done, std::size_t total) {
            live_probe(probe_round, point, done, total);
        });
        probe_guard.context = context;
    }
    for (std::size_t round = first_round; round <= final_round; ++round) {
        probe_round = round;
        live_probe(round, "round_entry");
        struct RoundClock {
            double& total;
            StageClock::time_point start = StageClock::now();
            ~RoundClock() { total += std::chrono::duration<double>(StageClock::now() - start).count(); }
        } round_clock{lap.rounds_wall};
        const auto timed_snapshot = [&](const char* kind, std::size_t index, const char* point, bool alive, bool relief) {
            const auto t = StageClock::now();
            snapshot(kind, index, point, alive, relief);
            lap.memory_snapshots += since(t);
        };
        RepairRound rec;
        rec.index = round;
        const auto round_mark = std::chrono::steady_clock::now();
        const char* kind = holding_ ? "hold" : round_kind_.c_str();
        const std::size_t snapshot_index = holding_ ? hold_index_ : round;

        // The incumbent, spelled once per round: its block_intervals are the coordinates every
        // replacement below is expressed in.
        auto lap_mark = StageClock::now();
        const HaplotypeSequence incumbent_a = spell_mosaic(model_, MosaicPath{cur_first_});
        const HaplotypeSequence incumbent_b = spell_mosaic(model_, MosaicPath{cur_second_});
        lap.incumbent_spelling += since(lap_mark);
        lap_mark = StageClock::now();
        if (use_delta && !have_inc_slots_) {
            // The initial pass: place the incumbent in full, then mark everything after it
            // transient. Owned cache keys mean a released haplotype's entries stay valid and keep
            // serving. The transient round opens BEFORE the seed is added, so a seed haplotype is
            // released like any other former incumbent: once a move replaces it nothing deltas from
            // it again, and cycle prevention needs only its allele path, which the visited set holds.
            open_transient_round();
            const std::uint64_t seeding_before = context->cumulative_counters().physical_origins;
            live_probe(round, "before_seed_placement");
            // I2-D: fused, the censuses come from the same pass and the retain_census calls below place nothing
            IncrementalBatchReport seed_report;   // S2-I: an output-only report, captured only for the diagnostic
            const std::vector<std::size_t> seeded =
                context->add_haplotypes({incumbent_a.bases, incumbent_b.bases},
                                        options_.memory_live_path.empty() ? nullptr : &seed_report, false, options_.fused_seed);
            live_probe(round, "after_seed_placement");
            if (!options_.memory_live_path.empty()) write_seed_report(seed_report);
            if (hooks.on_checkpoint) hooks.on_checkpoint("initial_incumbent_placed");
            inc_slot_a_ = seeded.front();
            inc_slot_b_ = seeded.back();
            role_of_[inc_slot_a_] = inc_slot_a_ == inc_slot_b_ ? "seed_hap1+2" : "seed_hap1";
            if (inc_slot_b_ != inc_slot_a_) role_of_[inc_slot_b_] = "seed_hap2";
            live_probe(round, "before_census_retention");
            context->retain_census(inc_slot_a_);
            context->retain_census(inc_slot_b_);
            live_probe(round, "after_census_retention");
            if (hooks.on_checkpoint) hooks.on_checkpoint("initial_incumbent_census_retained");
            seeding_origins_ += context->cumulative_counters().physical_origins - seeding_before;
            have_inc_slots_ = true;
        }
        lap.seed_placement_census += since(lap_mark);
        if (use_delta) rec.transient_before = context->transient_payload_bytes();
        timed_snapshot(kind, snapshot_index, "before_candidates", false, false);
        const double lookup_before = emission_cache_lookup_seconds();
        const auto generation_mark = StageClock::now();
        double spelling_inside = 0.0, priors_inside = 0.0;

        // S2-J: reuse-only origin-cache admission in round 1 only, restored when this round ends however it ends
        struct AdmissionScope {
            PreparedScoringContext* context = nullptr;
            bool previous = false;
            ~AdmissionScope() {
                if (context != nullptr) context->set_reuse_only_admission(previous);
            }
        } admission_scope;
        if (round == 1 && options_.round1_reuse_only_admission && use_delta && context != nullptr) {
            admission_scope.context = context;
            admission_scope.previous = context->reuse_only_admission();
            context->set_reuse_only_admission(true);
        }
        live_probe(round, "before_neighbourhood_construction");
        // ---- neighbourhood construction. NO TRUTH IN SCOPE.
        std::vector<MosaicPath> paths;
        std::vector<std::string> sequences;
        std::map<std::vector<std::size_t>, std::size_t> slot_of_path;
        const auto slot_for = [&](const std::vector<std::size_t>& alleles) {
            const auto found = slot_of_path.find(alleles);
            if (found != slot_of_path.end()) return found->second;
            const std::size_t slot = paths.size();
            slot_of_path.emplace(alleles, slot);
            paths.push_back(MosaicPath{alleles});
            const auto spell_mark = StageClock::now();
            sequences.push_back(spell_mosaic(model_, MosaicPath{alleles}).bases);
            spelling_inside += since(spell_mark);
            return slot;
        };
        struct Neighbour {
            std::size_t s1 = 0, s2 = 0, block = 0, rank = 0, x = 0, y = 0;
            int assignment = 0;
            std::vector<std::size_t> first, second;
        };
        // How each local haplotype relates to the incumbent. A neighbour differs from it at exactly
        // one block on one side, which is precisely what the delta needs to be told.
        struct Provenance {
            bool is_incumbent = false;
            int side = 0;            // 0: derived from the incumbent's first haplotype, 1: second
            std::size_t block = 0;
            std::size_t allele = 0;
        };
        std::vector<Provenance> provenance;
        const auto slot_with = [&](const std::vector<std::size_t>& alleles, Provenance note) {
            const std::size_t before = paths.size();
            const std::size_t slot = slot_for(alleles);
            if (paths.size() > before) {
                provenance.resize(paths.size());
                provenance[slot] = note;
            }
            return slot;
        };
        std::vector<Neighbour> neighbours;
        std::set<DiplotypeIdentity> seen_round;
        const std::size_t cur_s1 = slot_with(cur_first_, Provenance{true, 0, 0, 0});
        const std::size_t cur_s2 = slot_with(cur_second_, Provenance{true, 1, 0, 0});
        seen_round.insert(key_of(cur_first_, cur_second_));
        for (std::size_t b = 0; b < blocks; ++b) {
            const std::vector<std::pair<std::size_t, std::size_t>>& local = local_for(b, cur_first_, cur_second_);   // M9-G
            for (std::size_t r = 0; r < local.size(); ++r) {
                const std::size_t x = local[r].first;
                const std::size_t y = local[r].second;
                for (int assignment = 0; assignment < 2; ++assignment) {
                    if (assignment == 1 && x == y) continue;   // the swap is the same diplotype
                    std::vector<std::size_t> first = cur_first_;
                    std::vector<std::size_t> second = cur_second_;
                    first[b] = assignment == 0 ? x : y;
                    second[b] = assignment == 0 ? y : x;
                    ++rec.generated;
                    const auto key = key_of(first, second);
                    if (key == key_of(cur_first_, cur_second_)) {
                        ++rec.self_skipped;
                        continue;
                    }
                    if (visited_.count(key) != 0) {
                        ++rec.suppressed;
                        continue;
                    }
                    if (!seen_round.insert(key).second) {
                        ++rec.duplicate;
                        continue;
                    }
                    const std::size_t s1 = slot_with(first, Provenance{false, 0, b, first[b]});
                    const std::size_t s2 = slot_with(second, Provenance{false, 1, b, second[b]});
                    const auto prior_mark = StageClock::now();
                    const double prior = log_unordered_sequence_pair_prior(model_, sequences[s1], sequences[s2], switches_);
                    priors_inside += since(prior_mark);
                    if (!std::isfinite(prior)) {
                        ++rec.rejected_prior;
                        continue;
                    }
                    neighbours.push_back({s1, s2, b, r + 1, x, y, assignment, first, second});
                }
            }
        }
        lap.neighbour_spelling += spelling_inside;
        lap.neighbour_priors += priors_inside;
        lap.neighbour_generation += since(generation_mark) - spelling_inside - priors_inside;
        rec.neighbours = neighbours.size();
        live_probe(round, "after_neighbourhood_construction", 0, neighbours.size());
        if (request_order_seed != 0) {
            std::mt19937_64 rng(request_order_seed * 1000003ULL + round);
            std::shuffle(neighbours.begin(), neighbours.end(), rng);
        }
        if (neighbours.empty()) {
            stop_reason_ = "no_unvisited_neighbour";
            rounds_.push_back(rec);
            break;
        }

        // ---- ONE scoring pass: the current candidate, then its neighbours
        const auto scoring_mark = StageClock::now();
        std::map<std::pair<std::size_t, std::size_t>, std::size_t> scored_index;
        std::vector<HaplotypeCandidate> scored_pairs;
        std::vector<std::size_t> pair_of_item;
        std::vector<std::size_t> slot_of_local;
        std::map<std::size_t, const std::string*> sequence_of_slot;
        IncrementalBatchReport batch;
        std::optional<InternedHaplotypes> interned;
        BlockDeltaWork round_delta;
        const auto accumulate = [&](const BlockDeltaWork& work) {
            round_delta.origins_reused += work.origins_reused;
            round_delta.origins_recomputed += work.origins_recomputed;
            round_delta.local_window_bases += work.local_window_bases;
            round_delta.piece_hits += work.piece_hits;
            round_delta.dp_executions += work.dp_executions;
            round_delta.cache_misses += work.cache_misses;
            round_delta.cache_requests += work.cache_requests;
            round_delta.cache_hits += work.cache_hits;
            round_delta.full_haplotype_bytes_materialised += work.full_haplotype_bytes_materialised;
            round_delta.local_window_bytes_materialised += work.local_window_bytes_materialised;
            round_delta.sequence_key_bytes += work.sequence_key_bytes;
            round_delta.seconds_placement_and_origins += work.seconds_placement_and_origins;
            round_delta.seconds_statekey_merge += work.seconds_statekey_merge;
            round_delta.seconds_hashing_and_windows += work.seconds_hashing_and_windows;
            round_delta.placement_row_constructions += work.placement_row_constructions;
            round_delta.accumulator_growths += work.accumulator_growths;
            round_delta.origin_output_growths += work.origin_output_growths;
            round_delta.origin_cache_hits += work.origin_cache_hits;
            round_delta.origin_cache_misses += work.origin_cache_misses;
            round_delta.origin_cache_stored += work.origin_cache_stored;
            round_delta.origin_cache_declined += work.origin_cache_declined;
            round_delta.origin_cache_not_selected += work.origin_cache_not_selected;
            round_delta.origin_cache_admitted_bytes += work.origin_cache_admitted_bytes;
            round_delta.origin_cache_declined_bytes += work.origin_cache_declined_bytes;
            round_delta.seconds_reenumeration += work.seconds_reenumeration;
            round_delta.cached_origin_records_consumed += work.cached_origin_records_consumed;
            round_delta.physical_origins_constructed += work.physical_origins_constructed;
            round_delta.fragments_reduced += work.fragments_reduced;
            round_delta.fragments_copied += work.fragments_copied;
            round_delta.seconds_enumeration += work.seconds_enumeration;
            round_delta.seconds_affected_filter += work.seconds_affected_filter;
            round_delta.seconds_aggregate_merge += work.seconds_aggregate_merge;
            round_delta.seconds_reduce += work.seconds_reduce;
            round_delta.seconds += work.seconds;
            round_delta.origin_group_reuses += work.origin_group_reuses;
            round_delta.origin_cache_evicted += work.origin_cache_evicted;
            round_delta.origin_cache_evicted_bytes += work.origin_cache_evicted_bytes;
        };
        const auto resolve_by_delta = [&](std::size_t local, bool keep_census) {
            const Provenance& note = provenance[local];
            if (note.is_incumbent) return note.side == 0 ? inc_slot_a_ : inc_slot_b_;
            const HaplotypeSequence& base = note.side == 0 ? incumbent_a : incumbent_b;
            const std::size_t from = note.side == 0 ? inc_slot_a_ : inc_slot_b_;
            const SequenceInterval span = base.block_intervals[note.block];
            BlockDeltaWork work;
            const std::size_t slot = context->add_repaired_haplotype(
                from, span.begin, span.end, model_.blocks()[note.block].alleles[note.allele].sequence, &work,
                keep_census);
            accumulate(work);
            return slot;
        };

        if (use_delta) {
            provenance.resize(paths.size());
            slot_of_local.resize(paths.size());
            const std::size_t held_before = context->haplotypes_created();
            if (planned_cache) {
                // The round's repairs as ONE batch: the context plans them, then executes them in its
                // canonical order, so how this loop enumerated them cannot matter.
                std::vector<std::size_t> requested;
                for (std::size_t i = 0; i < paths.size(); ++i) {
                    const Provenance& note = provenance[i];
                    if (note.is_incumbent) {
                        slot_of_local[i] = note.side == 0 ? inc_slot_a_ : inc_slot_b_;
                    } else {
                        requested.push_back(i);
                    }
                }
                if (request_order_seed != 0) {
                    std::mt19937_64 rng(request_order_seed * 7919ULL + round);
                    std::shuffle(requested.begin(), requested.end(), rng);
                }
                std::vector<RepairRequest> requests;
                for (const std::size_t i : requested) {
                    const Provenance& note = provenance[i];
                    const HaplotypeSequence& base = note.side == 0 ? incumbent_a : incumbent_b;
                    const SequenceInterval span = base.block_intervals[note.block];
                    requests.push_back({note.side == 0 ? inc_slot_a_ : inc_slot_b_, span.begin, span.end,
                                        model_.blocks()[note.block].alleles[note.allele].sequence});
                }
                std::vector<BlockDeltaWork> works;
                OriginCachePlanReport plan;
                const std::vector<std::size_t> slots = context->add_repaired_haplotypes(requests, &works, &plan);
                for (std::size_t k = 0; k < requested.size(); ++k) {
                    slot_of_local[requested[k]] = slots[k];
                    accumulate(works[k]);
                }
                rec.plan_distinct = plan.distinct_windows;
                rec.plan_known = plan.known_windows;
                rec.plan_selected = plan.selected_entries;
                rec.plan_selected_bytes = plan.selected_bytes;
                rec.plan_retained = plan.retained_entries;
                rec.plan_evicted = plan.evicted_entries;
                rec.plan_evicted_bytes = plan.evicted_bytes;
            } else {
                for (std::size_t i = 0; i < paths.size(); ++i) slot_of_local[i] = resolve_by_delta(i, false);
            }
            for (std::size_t i = 0; i < sequences.size(); ++i) sequence_of_slot.emplace(slot_of_local[i], &sequences[i]);
            rec.haps_new = context->haplotypes_created() - held_before;
            rec.haps_reused = paths.size() - rec.haps_new;
            rec.unique_haplotypes = context->haplotypes_created();
        } else if (context) {
            slot_of_local = context->add_haplotypes(sequences, &batch);
            for (std::size_t i = 0; i < sequences.size(); ++i) sequence_of_slot.emplace(slot_of_local[i], &sequences[i]);
            rec.haps_new = batch.haplotypes_new;
            rec.haps_reused = batch.haplotypes_reused;
            rec.unique_haplotypes = batch.unique_haplotypes_total;
        } else {
            interned.emplace(intern_haplotypes(paths, sequences));
            slot_of_local = interned->index_of;
            for (std::size_t i = 0; i < sequences.size(); ++i) {
                sequence_of_slot.emplace(slot_of_local[i], &interned->sequences[slot_of_local[i]]);
            }
            rec.haps_new = interned->sequences.size();
            rec.unique_haplotypes = interned->sequences.size();
        }
        const auto add_item = [&](std::size_t a, std::size_t b2) {
            std::size_t x = slot_of_local[a];
            std::size_t y = slot_of_local[b2];
            if (y < x) std::swap(x, y);
            const auto inserted = scored_index.emplace(std::make_pair(x, y), scored_pairs.size());
            if (inserted.second) scored_pairs.push_back({x, y});
            pair_of_item.push_back(inserted.first->second);
        };
        add_item(cur_s1, cur_s2);                                  // item 0 is always the incumbent
        for (const Neighbour& n : neighbours) add_item(n.s1, n.s2);
        rec.scored_pairs = scored_pairs.size();

        live_probe(round, "before_candidate_combination", 0, scored_pairs.size());
        const auto mark = std::chrono::steady_clock::now();
        const FragmentScoringResult scored = context ? context->score_candidates(scored_pairs)
                                                     : score_diplotypes(*interned, scored_pairs, fragments_, library_, scoring_);
        live_probe(round, "after_candidate_combination", 0, scored_pairs.size());
        const std::size_t view_count = 1 + scored.audit_scores.size();
        if (options_.score_view >= view_count) throw std::logic_error("the repair search's score view was not scored");
        const auto view_of = [&](std::size_t v) -> const std::vector<DiplotypeScore>& {
            return v == 0 ? scored.scores : scored.audit_scores[v - 1];
        };
        const std::vector<DiplotypeScore>& own = view_of(options_.score_view);
        rec.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - mark).count();
        if (context) rec.seconds += batch.seconds_placement_and_mass;
        // The deltas ran while slots were being resolved, which is before `mark`.
        if (use_delta) rec.seconds += round_delta.seconds;
        lap.scoring += since(scoring_mark);
        timed_snapshot(kind, snapshot_index, "after_scoring", true, false);
        rec.rss = peak_rss_gib();
        {
            const FragmentScoringCounters& work = context ? batch.counters : scored.counters;
            const EmissionCacheCounters& cache = context ? batch.emission_cache : scored.emission_cache;
            rec.bases_scanned = work.candidate_bases_scanned;
            rec.piece_hits = work.piece_hits;
            rec.cache_requests = cache.dp_cache_requests;
            rec.cache_hits = cache.dp_cache_hits;
            rec.cache_misses = cache.dp_cache_misses;
            rec.executed_dps = cache.executed_dps;
            rec.physical_origins = work.physical_origins;
            rec.state_classes = work.statistical_state_classes;
            if (use_delta) {
                rec.bases_scanned = round_delta.local_window_bases;
                rec.piece_hits = round_delta.piece_hits;
                rec.cache_misses = round_delta.cache_misses;
                rec.cache_requests = round_delta.cache_requests;
                rec.cache_hits = round_delta.cache_hits;
                rec.executed_dps = round_delta.dp_executions;
                rec.origins_reused = round_delta.origins_reused;
                rec.origins_recomputed = round_delta.origins_recomputed;
                rec.window_bases = round_delta.local_window_bases;
                rec.full_bytes = round_delta.full_haplotype_bytes_materialised;
                rec.window_bytes = round_delta.local_window_bytes_materialised;
                rec.key_bytes = round_delta.sequence_key_bytes;
                rec.seconds_origins = round_delta.seconds_placement_and_origins;
                rec.seconds_merge = round_delta.seconds_statekey_merge;
                rec.seconds_hashing = round_delta.seconds_hashing_and_windows;
                rec.row_constructions = round_delta.placement_row_constructions;
                rec.accumulator_growths = round_delta.accumulator_growths;
                rec.origin_growths = round_delta.origin_output_growths;
                rec.oc_hits = round_delta.origin_cache_hits;
                rec.oc_misses = round_delta.origin_cache_misses;
                rec.oc_stored = round_delta.origin_cache_stored;
                rec.origins_constructed = round_delta.physical_origins_constructed;
                rec.origins_loaded = round_delta.cached_origin_records_consumed;
                rec.fragments_reduced = round_delta.fragments_reduced;
                rec.fragments_copied = round_delta.fragments_copied;
                rec.seconds_enumeration = round_delta.seconds_enumeration;
                rec.seconds_filter = round_delta.seconds_affected_filter;
                rec.seconds_aggregate = round_delta.seconds_aggregate_merge;
                rec.seconds_reduce = round_delta.seconds_reduce;
                rec.oc_entries = context->origin_cache_entries();
                rec.oc_bytes = context->origin_cache_bytes();
                rec.merge_working_bytes = context->merge_working_bytes();
            }
            // Not measured unless lookup timing is on: NaN, written NA, never a zero that reads as a measurement.
            rec.seconds_lookup = emission_lookup_timing_enabled() ? emission_cache_lookup_seconds() - lookup_before
                                                                  : std::numeric_limits<double>::quiet_NaN();
            {
                const EmissionCacheCounters live = read_emission_cache_counters();
                rec.entry_metadata_bytes = live.entry_metadata_bytes;
                rec.content_arena_bytes = live.content_arena_bytes;
            }
            rec.placement_seconds = use_delta ? round_delta.seconds
                                              : (context ? batch.seconds_placement_and_mass
                                                         : scored.seconds_placement_and_mass);
            rec.combination_seconds = scored.seconds_diploid_combination;
        }
        if (!holding_) {
            cumulative_seconds_ += rec.seconds;
            cumulative_pairs_ += scored_pairs.size();
        }

        lap_mark = StageClock::now();
        std::vector<double> prior_of_pair(scored_pairs.size(), 0.0);
        for (std::size_t slot = 0; slot < scored_pairs.size(); ++slot) {
            prior_of_pair[slot] = log_unordered_sequence_pair_prior(
                model_, *sequence_of_slot.at(scored_pairs[slot].first), *sequence_of_slot.at(scored_pairs[slot].second),
                switches_);
        }
        lap.scored_pair_priors += since(lap_mark);

        // ---- the guard: the incumbent must still score what it scored before any neighbour is read.
        // Round 1 checks against the start's declared score, which is what ties the search to the
        // ranking it continues rather than to one it produced itself.
        const double incumbent = own[pair_of_item[0]].log_likelihood + prior_of_pair[pair_of_item[0]];
        if (std::fabs(incumbent - cur_total_) > 1e-6) {
            failed_round_ = round;
            failed_expected_ = cur_total_;
            failed_scored_ = incumbent;
            stop_reason_ = "incumbent_not_reproduced";
            return false;
        }
        // The scored value is the authoritative one; a declared start score is only an assertion about
        // it, possibly through text at finite precision. Adopting it makes a resumed leg byte-identical.
        cur_total_ = incumbent;
        rec.start_total = incumbent;
        lap_mark = StageClock::now();
        if (hooks.on_scored && !holding_) {
            for (std::size_t item = 0; item < pair_of_item.size(); ++item) {
                const std::size_t at = pair_of_item[item];
                const std::size_t s1 = item == 0 ? cur_s1 : neighbours[item - 1].s1;
                const std::size_t s2 = item == 0 ? cur_s2 : neighbours[item - 1].s2;
                hooks.on_scored(item == 0 ? cur_first_ : neighbours[item - 1].first,
                                item == 0 ? cur_second_ : neighbours[item - 1].second, sequences[s1], sequences[s2],
                                own[at], prior_of_pair[at]);
                if (hooks.on_scored_views) {
                    std::vector<DiplotypeScore> views;
                    for (std::size_t v = 0; v < view_count; ++v) views.push_back(view_of(v)[at]);
                    hooks.on_scored_views(views);
                }
            }
        }

        lap.pool_hook += since(lap_mark);
        const auto certification_mark = StageClock::now();
        // ---- certify the round with production's own certifier, in presentation order
        const auto certify = [&](bool reversed, const std::vector<DiplotypeScore>& view_scores) {
            std::vector<std::size_t> order(pair_of_item.size());
            for (std::size_t i = 0; i < order.size(); ++i) order[i] = i;
            if (reversed) std::reverse(order.begin(), order.end());
            std::vector<DiplotypeScore> scores;
            std::vector<double> priors;
            std::vector<std::size_t> item_of_local;
            std::set<std::size_t> taken;
            for (const std::size_t i : order) {
                const std::size_t slot = pair_of_item[i];
                if (!taken.insert(slot).second) continue;
                scores.push_back(view_scores[slot]);
                priors.push_back(prior_of_pair[slot]);
                item_of_local.push_back(i);
            }
            const CertifiedRanking ranking = certify_scored_ranking(
                scores, fragments_.size(), scoring_.mismatch_band, scored.counters, CertificationOptions{}, priors);
            return std::make_pair(ranking, item_of_local);
        };
        const auto forward = certify(false, own);
        const CertifiedRanking& ranking = forward.first;
        const std::size_t winner_item = forward.second[ranking.candidates.front().candidate_id];
        rec.complete = ranking.status == RescoreStatus::Complete;
        rec.unique_best = ranking.winner_set.size() == 1;
        rec.separation = ranking.certified_gap;
        rec.incumbent_in_top_set = false;
        for (const std::size_t id : ranking.winner_set) {
            if (forward.second[id] == 0) rec.incumbent_in_top_set = true;
        }

        // Order invariance, from the scores already computed: no new likelihood call.
        {
            const auto back = certify(true, own);
            const std::size_t back_item = back.second[back.first.candidates.front().candidate_id];
            const auto alleles_of = [&](std::size_t item) {
                if (item == 0) return std::minmax(cur_first_, cur_second_);
                return std::minmax(neighbours[item - 1].first, neighbours[item - 1].second);
            };
            if (alleles_of(back_item) != alleles_of(winner_item) ||
                back.first.candidates.front().score != ranking.candidates.front().score) {
                rec.order_invariant = false;
            }
        }
        // Homologue swap: rebuilding the neighbourhood with the incumbent's haplotypes exchanged must
        // give the identical key set, because both assignments are generated for every local pair.
        {
            std::set<DiplotypeIdentity> swapped;
            for (std::size_t b = 0; b < blocks; ++b) {
                // M9-G: the same incumbent (unordered), so the same memoised list as the neighbourhood above
                const std::vector<std::pair<std::size_t, std::size_t>>& local = local_for(b, cur_second_, cur_first_);
                for (std::size_t r = 0; r < local.size(); ++r) {
                    const std::size_t x = local[r].first;
                    const std::size_t y = local[r].second;
                    for (int assignment = 0; assignment < 2; ++assignment) {
                        if (assignment == 1 && x == y) continue;
                        std::vector<std::size_t> first = cur_second_;   // homologues exchanged
                        std::vector<std::size_t> second = cur_first_;
                        first[b] = assignment == 0 ? x : y;
                        second[b] = assignment == 0 ? y : x;
                        swapped.insert(key_of(first, second));
                    }
                }
            }
            for (const Neighbour& n : neighbours) {
                if (swapped.count(key_of(n.first, n.second)) == 0) rec.swap_invariant = false;
            }
        }

        // ---- the move rule
        if (winner_item != 0) {
            const Neighbour& best = neighbours[winner_item - 1];
            const std::size_t slot = pair_of_item[winner_item];
            rec.block = best.block;
            rec.rank = best.rank;
            rec.assignment = best.assignment;
            rec.allele_x = best.x;
            rec.allele_y = best.y;
            rec.first = best.first;
            rec.second = best.second;
            rec.fragment = own[slot].fragment_term;
            rec.exposure = own[slot].exposure_term;
            rec.prior = prior_of_pair[slot];
            rec.total = own[slot].log_likelihood + prior_of_pair[slot];
        } else {
            rec.first = cur_first_;
            rec.second = cur_second_;
            const std::size_t slot = pair_of_item[0];
            rec.fragment = own[slot].fragment_term;
            rec.exposure = own[slot].exposure_term;
            rec.prior = prior_of_pair[slot];
            rec.total = incumbent;
        }

        if (!rec.complete) {
            stop_reason_ = "rescore_incomplete";
        } else if (winner_item == 0) {
            stop_reason_ = "no_improvement";
        } else if (!rec.unique_best) {
            stop_reason_ = "tie_at_best";
        } else if (rec.total <= cur_total_) {
            stop_reason_ = "no_improvement";
        } else if (visited_.count(key_of(rec.first, rec.second)) != 0) {
            stop_reason_ = "cycle";                                  // suppression should preclude this
        } else {
            rec.moved = true;
        }
        // ---- the same round at every scored view (library_sensitivity_audit_v1). The neighbourhood is already scored
        // under each audit library, so each view's move rule is certification and comparison only; the search itself
        // follows its own view. A view whose rule disagrees with the search's own decision is where its trajectory forks.
        if (view_count > 1 && !holding_) {
            rec.verdicts.resize(view_count);
            for (std::size_t v = 0; v < view_count; ++v) {
                const std::vector<DiplotypeScore>& at = view_of(v);
                const auto certified = certify(false, at);
                const CertifiedRanking& r = certified.first;
                const std::size_t w = certified.second[r.candidates.front().candidate_id];
                RoundVerdict& verdict = rec.verdicts[v];
                verdict.complete = r.status == RescoreStatus::Complete;
                verdict.unique_best = r.winner_set.size() == 1;
                for (const std::size_t id : r.winner_set) {
                    if (certified.second[id] == 0) verdict.incumbent_in_top_set = true;
                }
                verdict.incumbent_total = at[pair_of_item[0]].log_likelihood + prior_of_pair[pair_of_item[0]];
                verdict.winner_total = at[pair_of_item[w]].log_likelihood + prior_of_pair[pair_of_item[w]];
                verdict.winner = w == 0 ? key_of(cur_first_, cur_second_)
                                        : key_of(neighbours[w - 1].first, neighbours[w - 1].second);
                if (!verdict.complete) {
                    verdict.stop = "rescore_incomplete";
                } else if (w == 0) {
                    verdict.stop = "no_improvement";
                } else if (!verdict.unique_best) {
                    verdict.stop = "tie_at_best";
                } else if (verdict.winner_total <= verdict.incumbent_total) {
                    verdict.stop = "no_improvement";
                } else if (visited_.count(verdict.winner) != 0) {
                    verdict.stop = "cycle";
                } else {
                    verdict.moved = true;
                    verdict.stop = "moved";
                }
            }
            // The search's own view must reproduce the decision the search just made, or the replay is not the rule.
            const RoundVerdict& self = rec.verdicts[options_.score_view];
            const bool agrees = self.moved == rec.moved && self.complete == rec.complete &&
                                self.unique_best == rec.unique_best &&
                                self.incumbent_in_top_set == rec.incumbent_in_top_set &&
                                (rec.moved ? self.winner == key_of(rec.first, rec.second) : self.stop == stop_reason_);
            if (!agrees) throw std::logic_error("a round's verdict at the search's own view does not reproduce its decision");
        }
        lap.certification_and_move += since(certification_mark);
        const auto promotion_mark = StageClock::now();
        const double snapshots_before_promotion = lap.memory_snapshots;
        if (use_delta) {
            std::vector<std::size_t> keep;
            if (rec.moved && winner_item != 0) {
                const Neighbour& best = neighbours[winner_item - 1];
                // A delta, not a full placement: the winner's census is built from the incumbent's.
                const std::size_t won_a = resolve_by_delta(best.s1, true);
                const std::size_t won_b = resolve_by_delta(best.s2, true);
                keep = {won_a, won_b};
                const std::string promoted = "promoted_r" + std::to_string(round);
                if (won_a != inc_slot_a_ && won_a != inc_slot_b_) role_of_[won_a] = promoted + "_hap1";
                if (won_b != inc_slot_a_ && won_b != inc_slot_b_ && won_b != won_a) role_of_[won_b] = promoted + "_hap2";
                inc_slot_a_ = won_a;
                inc_slot_b_ = won_b;
            } else {
                keep = {inc_slot_a_, inc_slot_b_};
            }
            timed_snapshot(kind, snapshot_index, "after_promotion", true, false);
            live_probe(round, "after_promotion");
            std::vector<std::size_t> kept = keep;
            kept.insert(kept.end(), pinned_.begin(), pinned_.end());
            context->release_transient(kept);
            for (auto it = role_of_.begin(); it != role_of_.end();) {
                if (std::find(keep.begin(), keep.end(), it->first) == keep.end()) {
                    it = role_of_.erase(it);
                } else {
                    ++it;
                }
            }
            rec.transient_after = context->transient_payload_bytes();
            timed_snapshot(kind, snapshot_index, "after_release", true,
                           options_.pressure_relief && holding_ && hold_remaining_ == 0);
            live_probe(round, "after_release");
            // Re-read AFTER the promotion: the one materialisation a round performs happens there.
            rec.full_bytes = round_delta.full_haplotype_bytes_materialised;
            rec.key_bytes = round_delta.sequence_key_bytes;
            rec.seconds_origins = round_delta.seconds_placement_and_origins;
            rec.seconds_merge = round_delta.seconds_statekey_merge;
            rec.seconds_hashing = round_delta.seconds_hashing_and_windows;
            rec.oc_hits_all = round_delta.origin_cache_hits;
            rec.oc_misses_all = round_delta.origin_cache_misses;
            rec.oc_admitted = round_delta.origin_cache_stored;
            rec.oc_admitted_bytes = round_delta.origin_cache_admitted_bytes;
            rec.oc_declined = round_delta.origin_cache_declined;
            rec.oc_declined_bytes = round_delta.origin_cache_declined_bytes;
            rec.oc_not_selected = round_delta.origin_cache_not_selected;
            rec.oc_group_reuses = round_delta.origin_group_reuses;
            rec.oc_evicted = round_delta.origin_cache_evicted;
            rec.oc_evicted_bytes = round_delta.origin_cache_evicted_bytes;
            rec.origins_constructed_all = round_delta.physical_origins_constructed;
            rec.seconds_enumeration_all = round_delta.seconds_enumeration;
            rec.seconds_reenumeration_all = round_delta.seconds_reenumeration;
        }
        lap.promotion_release += since(promotion_mark) - (lap.memory_snapshots - snapshots_before_promotion);
        if (context) {
            const FragmentScoringCounters& cumulative = context->cumulative_counters();
            rec.cum_piece_hits = cumulative.piece_hits;
            rec.cum_physical_origins = cumulative.physical_origins;
            rec.cum_state_classes = cumulative.statistical_state_classes;
            rec.cum_latent_pairs = cumulative.latent_alignment_pairs;
            rec.mass_cells = static_cast<std::uint64_t>(fragments_.size()) * context->unique_haplotypes();
            rec.records_consumed = context->origin_records_consumed();
            rec.census_count = context->census_slots().size();
            rec.census_bytes = context->retained_census_bytes();
            rec.oc_bytes_end = context->origin_cache_bytes();
        }
        {
            const FootprintSample footprint = physical_footprint();
            rec.phys_footprint = footprint.current;
            rec.phys_footprint_peak = footprint.peak;
            rec.rss_current = current_rss_bytes();
        }
        rec.round_wall_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - round_mark).count();
        if (!holding_) {
            rounds_.push_back(rec);
            live_probe(round, "round_committed");
            const auto hook_mark = StageClock::now();
            if (hooks.on_round) hooks.on_round(rounds_.back());
            lap.trajectory_hook += since(hook_mark);
        }
        if (!rec.moved) {
            if (hold_remaining_ > 0) {
                // Not a round: the counter does not advance, and nothing is recorded as one.
                holding_ = true;
                --hold_remaining_;
                ++hold_index_;
                --round;
                continue;
            }
            break;
        }
        visited_.insert(key_of(rec.first, rec.second));
        cur_first_ = rec.first;
        cur_second_ = rec.second;
        cur_total_ = rec.total;
        if (round == final_round) stop_reason_ = "round_cap";
    }
    return true;
}

std::string CoordinateRepairSearch::render_rounds(const std::vector<RepairRound>& rounds) const {
    const std::size_t final_round = last_round();
    std::ostringstream table;
    table << std::setprecision(12)
          << "round\tstart_total\tmoved\trepaired_block\tlocal_pair_rank\tassignment\tallele_x\t"
             "allele_y\thap1_alleles\thap2_alleles\tfragment_term\texposure_term\tprior\ttotal\t"
             "delta_total\tstatus\tunique_best\tseparation\tneighbours_generated\tneighbours_scored\t"
             "duplicate\tsuppressed_visited\tself_skipped\trejected_infinite_prior\t"
             "unique_haplotypes\tscored_pairs\torder_invariant\tswap_invariant\tseconds\tpeak_rss_gib\t"
             "haplotypes_new\thaplotypes_reused\tbases_scanned\tpiece_hits\tcache_requests\t"
             "cache_hits\tcache_misses\texecuted_dps\tphysical_origins\tstate_classes\t"
             "placement_seconds\tcombination_seconds\torigins_reused\torigins_recomputed\t"
             "window_bases\ttransient_before\ttransient_after\tfull_haplotype_bytes\t"
             "window_bytes\tsequence_key_bytes\tseconds_origins\tseconds_merge\t"
             "seconds_hashing\tseconds_lookup\tentry_metadata_bytes\tcontent_arena_bytes\t"
             "row_constructions\taccumulator_growths\torigin_growths\t"
             "oc_hits\toc_misses\toc_stored\toc_entries\toc_bytes\torigins_constructed\t"
             "origins_loaded\tfragments_reduced\tfragments_copied\tmerge_working_bytes\t"
             "seconds_enumeration\tseconds_filter\tseconds_aggregate\tseconds_reduce\t"
             "round_wall_seconds\toc_hits_all\toc_misses_all\toc_admitted\toc_admitted_bytes\t"
             "oc_declined\toc_declined_bytes\torigins_constructed_all\tseconds_enumeration_all\t"
             "seconds_reenumeration_all\tcum_piece_hits\tmargin_piece_hits\tcum_physical_origins\t"
             "margin_physical_origins\tcum_state_classes\tmargin_state_classes\tcum_latent_pairs\t"
             "margin_latent_pairs\tmass_cells\tmargin_mass_cells\tmargin_origin_cache_bytes\t"
             "records_consumed\tmargin_records_consumed\tmargin_neighbours\tmargin_rounds\t"
             "census_count\tcensus_bytes\tphys_footprint\tphys_footprint_peak\trss_current\t"
             "plan_distinct_windows\tplan_known_windows\tplan_selected_entries\tplan_selected_bytes\t"
             "plan_retained_entries\tplan_evicted_entries\tplan_evicted_bytes\toc_not_selected\t"
             "oc_group_reuses\toc_evicted_in_round\toc_evicted_in_round_bytes\n";
    for (const RepairRound& r : rounds) {
        table << r.index << '\t' << r.start_total << '\t' << (r.moved ? "true" : "false") << '\t'
              << (r.moved ? model_.blocks()[r.block].block_id : std::string("NA")) << '\t'
              << (r.moved ? std::to_string(r.rank) : std::string("NA")) << '\t'
              << (r.moved ? std::to_string(r.assignment) : std::string("NA")) << '\t'
              << (r.moved ? std::to_string(r.allele_x) : std::string("NA")) << '\t'
              << (r.moved ? std::to_string(r.allele_y) : std::string("NA")) << '\t' << join(r.first) << '\t'
              << join(r.second) << '\t' << r.fragment << '\t' << r.exposure << '\t' << r.prior << '\t' << r.total
              << '\t' << (r.total - r.start_total) << '\t' << (r.complete ? "COMPLETE" : "RESCORE_INCOMPLETE")
              << '\t' << (r.unique_best ? "true" : "false") << '\t' << r.separation << '\t' << r.generated << '\t'
              << r.neighbours << '\t' << r.duplicate << '\t' << r.suppressed << '\t' << r.self_skipped << '\t'
              << r.rejected_prior << '\t' << r.unique_haplotypes << '\t' << r.scored_pairs << '\t'
              << (r.order_invariant ? "true" : "false") << '\t' << (r.swap_invariant ? "true" : "false") << '\t'
              << r.seconds << '\t' << r.rss << '\t' << r.haps_new << '\t' << r.haps_reused << '\t'
              << r.bases_scanned << '\t' << r.piece_hits << '\t' << r.cache_requests << '\t' << r.cache_hits
              << '\t' << r.cache_misses << '\t' << r.executed_dps << '\t' << r.physical_origins << '\t'
              << r.state_classes << '\t' << r.placement_seconds << '\t' << r.combination_seconds << '\t'
              << r.origins_reused << '\t' << r.origins_recomputed << '\t' << r.window_bases << '\t'
              << r.transient_before << '\t' << r.transient_after << '\t' << r.full_bytes << '\t' << r.window_bytes
              << '\t' << r.key_bytes << '\t' << r.seconds_origins << '\t' << r.seconds_merge << '\t'
              << r.seconds_hashing << '\t';
        if (std::isnan(r.seconds_lookup)) {
            table << "NA";
        } else {
            table << r.seconds_lookup;
        }
        table << '\t' << r.entry_metadata_bytes << '\t'
              << r.content_arena_bytes << '\t' << r.row_constructions << '\t' << r.accumulator_growths << '\t'
              << r.origin_growths << '\t' << r.oc_hits << '\t' << r.oc_misses << '\t' << r.oc_stored << '\t'
              << r.oc_entries << '\t' << r.oc_bytes << '\t' << r.origins_constructed << '\t' << r.origins_loaded
              << '\t' << r.fragments_reduced << '\t' << r.fragments_copied << '\t' << r.merge_working_bytes << '\t'
              << r.seconds_enumeration << '\t' << r.seconds_filter << '\t' << r.seconds_aggregate << '\t'
              << r.seconds_reduce << '\t' << r.round_wall_seconds << '\t' << r.oc_hits_all << '\t'
              << r.oc_misses_all << '\t' << r.oc_admitted << '\t' << r.oc_admitted_bytes << '\t' << r.oc_declined
              << '\t' << r.oc_declined_bytes << '\t' << r.origins_constructed_all << '\t'
              << r.seconds_enumeration_all << '\t' << r.seconds_reenumeration_all << '\t' << r.cum_piece_hits
              << '\t' << margin(scoring_.max_piece_hits, r.cum_piece_hits) << '\t' << r.cum_physical_origins
              << '\t' << margin(scoring_.max_physical_origins, r.cum_physical_origins) << '\t'
              << r.cum_state_classes << '\t' << margin(scoring_.max_state_classes, r.cum_state_classes) << '\t'
              << r.cum_latent_pairs << '\t' << margin(scoring_.max_latent_alignment_pairs, r.cum_latent_pairs)
              << '\t' << r.mass_cells << '\t' << margin(scoring_.max_mass_cells, r.mass_cells) << '\t'
              << (options_.origin_cache ? margin(options_.origin_cache_bytes, r.oc_bytes_end) : std::string("none"))
              << '\t' << r.records_consumed << '\t' << margin(options_.origin_cache_records, r.records_consumed)
              << '\t' << margin(options_.declared_neighbour_cap, r.neighbours) << '\t' << (final_round - r.index)
              << '\t' << r.census_count << '\t' << r.census_bytes << '\t' << r.phys_footprint << '\t'
              << r.phys_footprint_peak << '\t' << r.rss_current << '\t' << r.plan_distinct << '\t' << r.plan_known
              << '\t' << r.plan_selected << '\t' << r.plan_selected_bytes << '\t' << r.plan_retained << '\t'
              << r.plan_evicted << '\t' << r.plan_evicted_bytes << '\t' << r.oc_not_selected << '\t'
              << r.oc_group_reuses << '\t' << r.oc_evicted << '\t' << r.oc_evicted_bytes << '\n';
    }
    return table.str();
}

// ---- the final decision -------------------------------------------------------------------------

bool SequenceDigest::operator<(const SequenceDigest& o) const {
    return std::tie(a, b, length) < std::tie(o.a, o.b, o.length);
}
bool SequenceDigest::operator==(const SequenceDigest& o) const {
    return a == o.a && b == o.b && length == o.length;
}

SequenceDigest sequence_digest(const std::string& bases) {
    SequenceDigest id;
    id.length = bases.size();
    std::uint64_t h1 = 1469598103934665603ULL;
    std::uint64_t h2 = 0x9e3779b97f4a7c15ULL;
    for (const char c : bases) {
        const std::uint64_t v = static_cast<unsigned char>(c);
        h1 = (h1 ^ v) * 1099511628211ULL;
        h2 = (h2 + v + 1) * 0xbf58476d1ce4e5b9ULL;
        h2 ^= h2 >> 31;
    }
    id.a = h1;
    id.b = h2;
    return id;
}

const char* top_set_decision_name(TopSetDecision decision) {
    switch (decision) {
        case TopSetDecision::Unique: return "UNIQUE";
        case TopSetDecision::GenotypeEquivalent: return "GENOTYPE_EQUIVALENT";
        case TopSetDecision::DistinctTie: return "DISTINCT_TIE";
    }
    return "DISTINCT_TIE";
}

PooledDecision decide_pooled_genotype(std::vector<PoolCandidate> pool, std::size_t block_count,
                                      std::size_t fragments, std::size_t mismatch_band) {
    if (pool.empty()) throw std::invalid_argument("the pooled decision needs at least one candidate");
    for (PoolCandidate& c : pool) {
        if (c.first.size() != block_count || c.second.size() != block_count) {
            throw std::invalid_argument("every pooled candidate needs one allele per block on both haplotypes");
        }
        // Canonical orientation of the identity; the sequences follow their haplotypes.
        if (c.second < c.first) {
            std::swap(c.first, c.second);
            std::swap(c.first_sequence, c.second_sequence);
        }
    }
    std::sort(pool.begin(), pool.end(), [](const PoolCandidate& x, const PoolCandidate& y) {
        return std::tie(x.first, x.second) < std::tie(y.first, y.second);
    });
    for (std::size_t i = 1; i < pool.size(); ++i) {
        if (pool[i].first == pool[i - 1].first && pool[i].second == pool[i - 1].second) {
            throw std::invalid_argument("pooled candidate identities must be unique");
        }
    }
    PooledDecision out;
    out.candidates = std::move(pool);
    const std::vector<PoolCandidate>& candidates = out.candidates;
    const auto genotype_of = [&](std::size_t i) {
        std::vector<std::pair<std::size_t, std::size_t>> row(block_count);
        for (std::size_t b = 0; b < block_count; ++b) {
            row[b] = std::minmax(candidates[i].first[b], candidates[i].second[b]);
        }
        return row;
    };

    // Sequence-level candidates, as the caller ranks: members spelling one unordered sequence pair are
    // one state. Groups are ordered by their sequence pair, members by identity: both canonical.
    std::map<std::pair<SequenceDigest, SequenceDigest>, std::size_t> group_of;
    for (std::size_t i = 0; i < candidates.size(); ++i) {
        SequenceDigest x = candidates[i].first_sequence;
        SequenceDigest y = candidates[i].second_sequence;
        if (y < x) std::swap(x, y);
        group_of.emplace(std::make_pair(x, y), 0);
    }
    {
        std::size_t next = 0;
        for (auto& item : group_of) item.second = next++;
    }
    out.groups.assign(group_of.size(), {});
    for (std::size_t i = 0; i < candidates.size(); ++i) {
        SequenceDigest x = candidates[i].first_sequence;
        SequenceDigest y = candidates[i].second_sequence;
        if (y < x) std::swap(x, y);
        out.groups[group_of.at({x, y})].push_back(i);
    }
    std::vector<DiplotypeScore> scores;
    std::vector<double> priors;
    for (const std::vector<std::size_t>& members : out.groups) {
        std::size_t best = members.front();
        bool ambiguous = false;
        for (const std::size_t m : members) {
            if (candidates[m].path_prior > candidates[best].path_prior) best = m;
            if (!same_score(candidates[m].score, candidates[members.front()].score) ||
                std::memcmp(&candidates[m].sequence_prior, &candidates[members.front()].sequence_prior,
                            sizeof(double)) != 0) {
                ++out.score_mismatches;
            }
        }
        const auto chosen = genotype_of(best);
        for (const std::size_t m : members) {
            if (genotype_of(m) != chosen) ambiguous = true;
        }
        out.representative.push_back(best);
        out.projection_ambiguous.push_back(ambiguous);
        scores.push_back(candidates[members.front()].score);
        priors.push_back(candidates[members.front()].sequence_prior);
    }
    out.ranking = certify_scored_ranking(scores, fragments, mismatch_band, FragmentScoringCounters{},
                                         CertificationOptions{}, priors);

    // ---- the numerical top set, in rank order
    const std::set<std::size_t> winners(out.ranking.winner_set.begin(), out.ranking.winner_set.end());
    for (const CertifiedCandidate& c : out.ranking.candidates) {
        if (winners.count(c.candidate_id) != 0) out.top_groups.push_back(c.candidate_id);
    }
    const CertifiedCandidate& best = out.ranking.candidates.front();
    out.genotype = genotype_of(out.representative[best.candidate_id]);
    {
        std::set<std::size_t> disagreeing;
        for (const std::size_t g : out.top_groups) {
            for (const std::size_t m : out.groups[g]) {
                const auto row = genotype_of(m);
                for (std::size_t b = 0; b < block_count; ++b) {
                    if (row[b] != out.genotype[b]) disagreeing.insert(b);
                }
            }
        }
        out.disagreeing_blocks.assign(disagreeing.begin(), disagreeing.end());
    }
    if (!out.disagreeing_blocks.empty()) {
        out.decision = TopSetDecision::DistinctTie;
    } else {
        out.decision = out.top_groups.size() == 1 ? TopSetDecision::Unique : TopSetDecision::GenotypeEquivalent;
    }

    // ---- separation of the genotype, and the best same-genotype alternative below the top set
    const auto carries_only = [&](std::size_t g) {
        for (const std::size_t m : out.groups[g]) {
            if (genotype_of(m) != out.genotype) return false;
        }
        return true;
    };
    for (const CertifiedCandidate& c : out.ranking.candidates) {
        if (winners.count(c.candidate_id) != 0) continue;
        if (!carries_only(c.candidate_id)) {
            if (!out.has_other_genotype) {
                out.has_other_genotype = true;
                out.best_other_genotype_group = c.candidate_id;
                out.genotype_separation_gap = best.score - c.score;
                out.genotype_separation_certified =
                    best.score - best.numeric_error_bound > c.score + c.numeric_error_bound;
            }
        } else if (!out.has_same_genotype_alternative) {
            out.has_same_genotype_alternative = true;
            out.best_same_genotype_group = c.candidate_id;
            out.same_genotype_score_gap = best.score - c.score;
            out.same_genotype_fragment_gap = best.fragment_term - c.fragment_term;
            out.same_genotype_prior_gap = best.prior - c.prior;
            out.same_genotype_fragment_identical =
                std::memcmp(&best.fragment_term, &c.fragment_term, sizeof(double)) == 0;
        }
        if (out.has_other_genotype && out.has_same_genotype_alternative) break;
    }
    if (!out.has_other_genotype) out.genotype_separation_gap = std::numeric_limits<double>::infinity();

    // ---- phase support and per-block read support, exactly as the unpooled caller reports them
    out.block_genotype_read_supported.assign(block_count, true);
    out.fragment_term_gap = out.ranking.candidates.size() > 1
                                ? best.fragment_term - out.ranking.candidates[1].fragment_term
                                : std::numeric_limits<double>::infinity();
    for (std::size_t i = 1; i < out.ranking.candidates.size(); ++i) {
        const CertifiedCandidate& rival = out.ranking.candidates[i];
        // A fragment contender: the reads alone do not put it below the winner.
        if (rival.fragment_term + rival.numeric_error_bound < best.fragment_term - best.numeric_error_bound) continue;
        ++out.fragment_contenders;
        const auto row = genotype_of(out.representative[rival.candidate_id]);
        for (std::size_t b = 0; b < block_count; ++b) {
            if (row[b] != out.genotype[b]) out.block_genotype_read_supported[b] = false;
        }
    }
    if (out.decision != TopSetDecision::Unique) {
        out.phase_support = "UNRESOLVED";
    } else {
        out.phase_support = out.fragment_contenders == 0 ? "READS_WITHIN_CANDIDATES" : "PRIOR_ONLY";
    }
    {
        std::vector<double> ranked_scores;
        std::vector<double> ranked_bounds;
        std::vector<std::vector<std::pair<std::size_t, std::size_t>>> ranked_genotypes;
        for (const CertifiedCandidate& entry : out.ranking.candidates) {
            ranked_scores.push_back(entry.score);
            ranked_bounds.push_back(entry.numeric_error_bound);
            ranked_genotypes.push_back(genotype_of(out.representative[entry.candidate_id]));
        }
        for (const BlockGenotypeConfidence& block :
             block_genotype_confidence(ranked_scores, ranked_bounds, ranked_genotypes)) {
            out.block_genotype_posterior.push_back(block.posterior);
            out.block_genotype_quality.push_back(block.quality);
            out.block_numerically_separated.push_back(block.numerically_separated);
        }
    }
    return out;
}

SymmetryReport check_phase_switch_symmetry(const std::vector<SymmetryTraceItem>& a,
                                           const std::vector<SymmetryTraceItem>& b, std::size_t block,
                                           std::size_t fragment_count) {
    SymmetryReport report;
    // certify_scored_ranking's per-candidate bound, with its default floor
    const auto error_bound = [&](const SymmetryTraceItem& item) {
        const double terms = static_cast<double>(fragment_count + 16);
        return std::max(terms * std::numeric_limits<double>::epsilon() *
                            (item.score.absolute_term_sum + std::fabs(item.prior)),
                        0.5 * CertificationOptions{}.minimum_tie_tolerance_nats);
    };
    const auto key_of = [](const SymmetryTraceItem& item, bool switch_block, std::size_t at) {
        std::vector<std::size_t> x = item.first;
        std::vector<std::size_t> y = item.second;
        if (switch_block) std::swap(x[at], y[at]);
        return diplotype_identity(x, y);
    };
    std::size_t rounds_a = 0, rounds_b = 0;
    for (const auto& item : a) rounds_a = std::max(rounds_a, item.round);
    for (const auto& item : b) rounds_b = std::max(rounds_b, item.round);
    const std::size_t rounds = std::max(rounds_a, rounds_b);
    for (std::size_t r = 1; r <= rounds; ++r) {
        std::map<DiplotypeIdentity, const SymmetryTraceItem*> in_a;
        SymmetryDepth depth;
        depth.round = r;
        for (const auto& item : a) {
            if (item.round == r && in_a.emplace(key_of(item, false, block), &item).second) ++depth.a_candidates;
        }
        std::set<DiplotypeIdentity> seen_b;
        for (const auto& item : b) {
            if (item.round != r || !seen_b.insert(key_of(item, true, block)).second) continue;
            ++depth.b_candidates;
            const auto found = in_a.find(key_of(item, true, block));
            if (found == in_a.end()) continue;
            ++depth.matched;
            const SymmetryTraceItem& other = *found->second;
            if (same_score(item.score, other.score) &&
                std::memcmp(&item.prior, &other.prior, sizeof(double)) == 0) {
                ++depth.bit_identical;
            }
            const double difference = std::fabs((item.score.log_likelihood + item.prior) -
                                                (other.score.log_likelihood + other.prior));
            const double bound = error_bound(item) + error_bound(other);
            if (difference <= bound) ++depth.indistinguishable;
            depth.max_total_difference = std::max(depth.max_total_difference, difference);
            depth.max_difference_over_bound = std::max(depth.max_difference_over_bound, difference / bound);
        }
        report.depths.push_back(depth);
        if (report.first_break_round == 0) {
            std::string why;
            if (r > rounds_a || r > rounds_b) why = "one branch stopped at an earlier round";
            else if (depth.a_candidates != depth.b_candidates || depth.matched != depth.b_candidates)
                why = "the switched neighbourhoods differ (" + std::to_string(depth.b_candidates - depth.matched) +
                      " candidates without a switched counterpart)";
            else if (depth.bit_identical != depth.matched)
                why = std::to_string(depth.matched - depth.bit_identical) + " switched candidates score differently";
            if (!why.empty()) {
                report.first_break_round = r;
                report.first_break = why;
            }
        }
        if (report.numeric_first_break_round == 0) {
            std::string why;
            if (r > rounds_a || r > rounds_b) why = "one branch stopped at an earlier round";
            else if (depth.a_candidates != depth.b_candidates || depth.matched != depth.b_candidates)
                why = "the switched neighbourhoods differ (" + std::to_string(depth.b_candidates - depth.matched) +
                      " candidates without a switched counterpart)";
            else if (depth.indistinguishable != depth.matched)
                why = std::to_string(depth.matched - depth.indistinguishable) +
                      " switched candidates differ beyond the certified numerical bound";
            if (!why.empty()) {
                report.numeric_first_break_round = r;
                report.numeric_first_break = why;
            }
        }
    }
    report.equivalent = rounds > 0 && report.first_break_round == 0;
    report.numerically_equivalent = rounds > 0 && report.numeric_first_break_round == 0;
    return report;
}


// ---- improving-move tie branching ----------------------------------------------------------------------

TieSet terminal_tie(const std::vector<TieItem>& items, std::size_t round, std::size_t fragments,
                    std::size_t mismatch_band) {
    TieSet tie;
    std::vector<const TieItem*> here;
    for (const TieItem& item : items) {
        if (item.round == round) here.push_back(&item);
    }
    if (here.empty()) throw std::invalid_argument("no scored candidate in the tied round");
    tie.incumbent = *here.front();
    std::map<std::pair<SequenceDigest, SequenceDigest>, std::size_t> slot;
    std::vector<std::vector<const TieItem*>> groups;
    std::vector<DiplotypeScore> scores;
    std::vector<double> priors;
    for (const TieItem* item : here) {
        SequenceDigest x = item->first_sequence;
        SequenceDigest y = item->second_sequence;
        if (y < x) std::swap(x, y);
        const auto found = slot.emplace(std::make_pair(x, y), groups.size());
        if (found.second) {
            groups.emplace_back();
            scores.push_back(item->score);
            priors.push_back(item->prior);
        }
        groups[found.first->second].push_back(item);
    }
    const CertifiedRanking ranking =
        certify_scored_ranking(scores, fragments, mismatch_band, FragmentScoringCounters{}, CertificationOptions{}, priors);
    const DiplotypeIdentity incumbent_key = diplotype_identity(tie.incumbent.first, tie.incumbent.second);
    for (const std::size_t g : ranking.winner_set) {
        for (const TieItem* item : groups[g]) {
            if (diplotype_identity(item->first, item->second) != incumbent_key) tie.tied.push_back(*item);
        }
    }
    std::sort(tie.tied.begin(), tie.tied.end(), [](const TieItem& x, const TieItem& y) {
        return diplotype_identity(x.first, x.second) < diplotype_identity(y.first, y.second);
    });
    return tie;
}

TieBranchingReport run_tie_branching(const TieSet& root, const std::vector<DiplotypeIdentity>& root_history,
                                     std::size_t max_branches, std::size_t fragments, std::size_t mismatch_band,
                                     const std::function<BranchResult(const BranchStart&)>& run_branch,
                                     const std::function<void(const TieSet&)>& on_tie, bool reverse,
                                     const std::function<void(long)>& on_settled) {
    TieBranchingReport report;
    const auto settled = [&](long branch) {
        if (on_settled) on_settled(branch);
    };
    std::vector<BranchStart> worklist;
    std::set<DiplotypeIdentity> started;
    const auto record = [&](TieSet tie, const std::vector<DiplotypeIdentity>& history) {
        tie.id = report.ties.size();
        if (tie.tied.empty()) throw std::runtime_error("a split without an improving tie");
        report.ties.push_back(tie);
        on_tie(report.ties.back());
        std::vector<TieItem> order = report.ties.back().tied;
        if (reverse) std::reverse(order.begin(), order.end());
        for (const TieItem& n : order) {
            if (!started.insert(diplotype_identity(n.first, n.second)).second) continue;   // one branch per start
            BranchStart branch;
            branch.parent_tie = static_cast<long>(tie.id);
            branch.start = n;
            branch.history = history;
            worklist.push_back(branch);
        }
    };
    try {
        record(root, root_history);
    } catch (const std::exception& error) {
        report.incomplete = std::string("the production tie could not be branched: ") + error.what();
        return report;
    }
    settled(-1);
    for (std::size_t next = 0; next < worklist.size(); ++next) {
        // Every unique start in the worklist -- run or waiting -- needs a branch of its own, so once there are more of
        // them than the cap allows, completing within the cap is impossible: fail closed now, before any further
        // seed is materialised or scored. With the worklist at or below the cap the search proceeds unchanged.
        if (worklist.size() > max_branches) {
            report.incomplete = "branch cap " + std::to_string(max_branches) + " binds (at least " +
                                std::to_string(worklist.size()) + " branches required; proven after " +
                                std::to_string(report.branches.size()) + " branches)";
            break;
        }
        if (report.branches.size() >= max_branches) {
            report.incomplete = "branch cap " + std::to_string(max_branches) + " binds (" +
                                std::to_string(worklist.size()) + " branches required)";
            break;
        }
        BranchRecord entry;
        entry.start = worklist[next];
        entry.start.id = report.branches.size();
        try {
            entry.result = run_branch(entry.start);
        } catch (const std::exception& error) {
            entry.result = BranchResult{};
            entry.result.outcome = BranchOutcome::Incomplete;
            entry.result.why = error.what();
        }
        report.branches.push_back(entry);
        const BranchRecord& done = report.branches.back();
        const long done_id = static_cast<long>(done.start.id);
        if (done.result.outcome == BranchOutcome::Incomplete) {
            report.incomplete = "branch " + std::to_string(done.start.id) + ": " + done.result.why;
            settled(done_id);
            break;
        }
        if (done.result.outcome == BranchOutcome::Converged) {
            std::string text;
            for (std::size_t b = 0; b < done.result.terminal_first.size(); ++b) {
                const auto pair = std::minmax(done.result.terminal_first[b], done.result.terminal_second[b]);
                text += (b == 0 ? "" : ",") + std::to_string(pair.first) + "/" + std::to_string(pair.second);
            }
            report.terminal_genotypes.insert(text);
            settled(done_id);
            continue;
        }
        std::vector<DiplotypeIdentity> history = done.start.history;
        history.push_back(diplotype_identity(done.start.start.first, done.start.start.second));
        history.insert(history.end(), done.result.moved_to.begin(), done.result.moved_to.end());
        try {
            TieSet nested = terminal_tie(done.result.items, done.result.last_round, fragments, mismatch_band);
            nested.branch = static_cast<long>(done.start.id);
            record(nested, history);
        } catch (const std::exception& error) {
            report.incomplete = "branch " + std::to_string(done.start.id) + ": " + error.what();
            settled(done_id);
            break;
        }
        settled(done_id);
    }
    return report;
}

LibrarySensitivity classify_library_sensitivity(const std::vector<std::pair<std::size_t, std::size_t>>& primary,
                                                const std::vector<LibraryAuditView>& views) {
    using PairSet = std::set<std::pair<std::size_t, std::size_t>>;
    LibrarySensitivity out;
    const std::size_t blocks = primary.size();
    out.projections.assign(views.size(), {});
    out.independent_projections.assign(views.size(), {});
    out.call_class.assign(views.size(), "incomplete");
    out.independent_call_class.assign(views.size(), "NA");
    out.top_set_size.assign(views.size(), 0);
    out.top_set_digest.assign(views.size(), "NA");
    out.representative_digest.assign(views.size(), "NA");
    out.independent_digest.assign(views.size(), "NA");
    std::set<std::string> anchors;
    std::string incomplete;
    const auto identity_line = [](const PoolCandidate& c) {
        const DiplotypeIdentity key = diplotype_identity(c.first, c.second);
        std::string line;
        for (const std::size_t a : key.first) line += std::to_string(a) + ",";
        line += "|";
        for (const std::size_t a : key.second) line += std::to_string(a) + ",";
        return line;
    };
    const auto digest_of = [](std::vector<std::string> lines) {
        std::sort(lines.begin(), lines.end());
        std::string text;
        for (const std::string& line : lines) text += line + "\n";
        return sha256_bytes(text).substr(0, 16);
    };
    const auto class_of = [](const PooledDecision& d) {
        return d.decision == TopSetDecision::Unique               ? "strict_call"
               : d.decision == TopSetDecision::GenotypeEquivalent ? "genotype_equivalent_call"
                                                                  : "no_call";
    };
    // The unordered pairs the top set's representatives carry at each block, and the digests of its identities.
    const auto project = [&](const PooledDecision& d, std::vector<PairSet>& at, std::string* members_digest,
                             std::string* representatives_digest, std::size_t* size) {
        at.assign(blocks, {});
        std::vector<std::string> members, representatives;
        for (const std::size_t g : d.top_groups) {
            const PoolCandidate& rep = d.candidates[d.representative[g]];
            representatives.push_back(identity_line(rep));
            for (std::size_t b = 0; b < blocks; ++b) at[b].insert(std::minmax(rep.first[b], rep.second[b]));
            for (const std::size_t m : d.groups[g]) members.push_back(identity_line(d.candidates[m]));
        }
        if (members_digest != nullptr) *members_digest = digest_of(members);
        if (representatives_digest != nullptr) *representatives_digest = digest_of(representatives);
        if (size != nullptr) *size = members.size();
    };
    for (std::size_t v = 0; v < views.size(); ++v) {
        const LibraryAuditView& view = views[v];
        anchors.insert(view.anchor);
        const bool missing = view.decision == nullptr || (v > 0 && view.independent == nullptr);
        if (!view.incomplete.empty() || missing) {
            if (incomplete.empty()) {
                incomplete = "view_" + std::to_string(v) + ":" + (view.incomplete.empty() ? "no_decision" : view.incomplete);
            }
            continue;
        }
        out.call_class[v] = class_of(*view.decision);
        project(*view.decision, out.projections[v], &out.top_set_digest[v], &out.representative_digest[v],
                &out.top_set_size[v]);
        if (v > 0) {
            out.independent_call_class[v] = class_of(*view.independent);
            project(*view.independent, out.independent_projections[v], &out.independent_digest[v], nullptr, nullptr);
        }
    }
    out.distinct_anchors = anchors.size();
    if (!views.empty() && !out.projections.front().empty()) {
        out.central_matches_primary = true;
        for (std::size_t b = 0; b < blocks; ++b) {
            if (out.projections.front()[b] != PairSet{primary[b]}) out.central_matches_primary = false;
        }
    }
    if (!incomplete.empty()) {
        out.incomplete = incomplete;
        return out;   // NOT_EVALUATED, LIBRARY_SENSITIVITY_NOT_EVALUATED: nothing flagged
    }
    out.block_reason.assign(blocks, "STABLE");
    for (std::size_t b = 0; b < blocks; ++b) {
        const PairSet reference{primary[b]};
        bool rank = false, pipeline = false;
        for (std::size_t v = 1; v < views.size(); ++v) {
            if (out.projections[v][b] != reference) rank = true;
            if (out.independent_projections[v][b] != reference) pipeline = true;
        }
        if (rank) out.rank_blocks.push_back(b);
        if (pipeline) out.pipeline_blocks.push_back(b);
        if (rank || pipeline) out.sensitive_blocks.push_back(b);
        out.block_reason[b] = rank && pipeline ? "BOTH" : rank ? "LIKELIHOOD_RANK" : pipeline ? "SEARCH_PATH" : "STABLE";
    }
    const bool any_rank = !out.rank_blocks.empty(), any_pipeline = !out.pipeline_blocks.empty();
    const bool any_search_path_only =
        std::find(out.block_reason.begin(), out.block_reason.end(), "SEARCH_PATH") != out.block_reason.end();
    const bool any_rank_only =
        std::find(out.block_reason.begin(), out.block_reason.end(), "LIKELIHOOD_RANK") != out.block_reason.end();
    const bool any_both = std::find(out.block_reason.begin(), out.block_reason.end(), "BOTH") != out.block_reason.end();
    out.status = out.sensitive_blocks.empty() ? "STABLE" : "SENSITIVE";
    out.reason = any_both || (any_rank_only && any_search_path_only) ? "BOTH"
                 : any_rank                                          ? "LIKELIHOOD_RANK"
                 : any_pipeline                                      ? "SEARCH_PATH"
                                                                     : "STABLE";
    if (out.sensitive_blocks.empty()) {
        for (std::size_t v = 1; v < views.size(); ++v) {
            if (out.top_set_digest[v] != out.top_set_digest[0] ||
                out.representative_digest[v] != out.representative_digest[0] ||
                out.independent_digest[v] != out.top_set_digest[0]) {
                out.phase_or_representative_only = true;
            }
        }
    }
    return out;
}

TieBranchVerdict tie_branch_verdict(bool incomplete, TopSetDecision decision) {
    if (incomplete) return {"REPAIR_INCOMPLETE", "no_call", "tie_branching_incomplete", "tie_branched"};
    switch (decision) {
    case TopSetDecision::Unique: return {"COMPLETE", "strict_call", "none", "tie_branched"};
    case TopSetDecision::GenotypeEquivalent:
        return {"CONDITIONAL_GENOTYPE_EQUIVALENT", "genotype_equivalent_call", "none", "tie_branched"};
    case TopSetDecision::DistinctTie: break;
    }
    return {"RESCORE_INCOMPLETE", "no_call", "divergent_tie", "tie_branched"};
}


// ---- representation-conflict quality guard -------------------------------------------------------------

std::pair<std::size_t, std::size_t> marker_best_pair(const BlockProposalEmissions& emission) {
    if (emission.allele_count == 0) throw std::invalid_argument("a block with no allele has no marker pair");
    std::pair<std::size_t, std::size_t> top{0, 0};
    double top_value = -std::numeric_limits<double>::infinity();
    for (std::size_t x = 0; x < emission.allele_count; ++x) {
        for (std::size_t y = x; y < emission.allele_count; ++y) {
            if (emission.at(x, y) > top_value) {   // strict: the first pair in index order wins a tie
                top_value = emission.at(x, y);
                top = {x, y};
            }
        }
    }
    return top;
}

std::pair<std::size_t, std::size_t> marker_best_pair(const BlockProposalEmissions& emission,
                                                     std::pair<std::size_t, std::size_t> excluded) {
    std::pair<std::size_t, std::size_t> top{0, 0};
    double top_value = -std::numeric_limits<double>::infinity();
    bool found = false;
    for (std::size_t x = 0; x < emission.allele_count; ++x) {
        for (std::size_t y = x; y < emission.allele_count; ++y) {
            if (std::make_pair(x, y) == excluded) continue;
            if (!found || emission.at(x, y) > top_value) {
                top_value = emission.at(x, y);
                top = {x, y};
                found = true;
            }
        }
    }
    if (!found) throw std::invalid_argument("a block with one allele pair has no runner-up");
    return top;
}

RepresentationEvidence representation_conflict(double marker_margin, const DiplotypeScore& called, double called_prior,
                                               std::vector<RepresentationAssignment> competitor_assignments,
                                               double bound) {
    if (competitor_assignments.empty()) throw std::invalid_argument("a representation competitor needs a scored assignment");
    std::sort(competitor_assignments.begin(), competitor_assignments.end(),
              [](const RepresentationAssignment& x, const RepresentationAssignment& y) {
                  const double tx = x.score.log_likelihood + x.prior;
                  const double ty = y.score.log_likelihood + y.prior;
                  if (tx != ty) return tx > ty;
                  if (x.score.fragment_term != y.score.fragment_term) return x.score.fragment_term > y.score.fragment_term;
                  return x.score.exposure_term > y.score.exposure_term;
              });
    const RepresentationAssignment& best = competitor_assignments.front();
    RepresentationEvidence out;
    out.marker_margin = marker_margin;
    out.fragment_diff = called.fragment_term - best.score.fragment_term;
    out.exposure_diff = called.exposure_term - best.score.exposure_term;
    out.likelihood_margin = (called.log_likelihood + called_prior) - (best.score.log_likelihood + best.prior);
    out.statistic = std::min({marker_margin, out.fragment_diff, -out.exposure_diff});
    out.conflict = out.statistic >= bound;
    return out;
}

} // namespace panvar
