// S2-C (MEASUREMENT ONLY): the progressive indexed kernel against the bulk adapter on real union inputs
// (<prefix>.union_kernel_inputs.bin from the caller, plus the reads). Each mode runs in its own process so peak memory is
// its own.
//   s2c_progressive_harness bulk        <inputs.bin> <reads_1.fq> <reads_2.fq> <out prefix>
//   s2c_progressive_harness progressive <inputs.bin> <reads_1.fq> <reads_2.fq> <out prefix> [prefix_size=31 [footprint_limit_bytes
//     [per_call|persist]]]   (S2-D: persist carries the emission cache across the batch's discovery calls)
//     a footprint guard (4 GiB by default, the local envelope) records the completed stages and exits 125 when crossed
//   s2c_progressive_harness blocked     <inputs.bin> <reads_1.fq> <reads_2.fq> <out prefix> [block=4096 [footprint_limit_bytes]]
//     S2-E: the fragment-blocked union pass; the footprint guard (4 GiB default) also covers bulk (argument 6 there)
//   s2c_progressive_harness compare     <inputs.bin> <bulk prefix> <progressive prefix> <out.tsv>
// The schedule (declared here, before any timing): stage 1 the S2-A early pair; stage 2 the next `prefix_size` pairs by
// the truth-free proposal score (descending; ties to pool order); stage 3 the rest of the pool in pool order. No pair is
// excluded; an interrupted schedule is provisional. Outputs: <prefix>.scores.bin (pool order) and <prefix>.tsv.
#include "panvar/fragment_likelihood.hpp"
#include "panvar/genotype_model.hpp"
#include "panvar/genotype_repair.hpp"

#include <sys/resource.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace {
using Clock = std::chrono::steady_clock;
double since(Clock::time_point t) { return std::chrono::duration<double>(Clock::now() - t).count(); }
double peak_rss_bytes() {
    rusage u{};
    getrusage(RUSAGE_SELF, &u);
    return static_cast<double>(u.ru_maxrss);
}

struct Inputs {
    panvar::InternedHaplotypes haplotypes;
    std::vector<std::vector<std::uint32_t>> segments;
    std::vector<panvar::HaplotypeCandidate> pool;
    std::vector<double> proposal, prior;
    std::uint64_t early_first = 0, early_second = 0;
    panvar::LibraryModel library;
    panvar::FragmentScoringOptions options;
    std::size_t store_bytes = 0;
    bool legacy = false, views = false;
};
Inputs read_inputs(const std::string& file) {
    std::ifstream in(file, std::ios::binary);
    if (!in) throw std::runtime_error("cannot read " + file);
    const std::string all((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    if (all.compare(0, 5, "PVKI1") != 0) throw std::runtime_error("not a PVKI1 file");
    std::size_t at = 8;
    const auto get = [&](auto& v) {
        std::memcpy(&v, all.data() + at, sizeof v);
        at += sizeof v;
    };
    const auto get_library = [&](panvar::LibraryModel& l) {
        get(l.insert_mean_bp); get(l.insert_sd_bp); get(l.support_sigmas); get(l.substitution_error);
        get(l.outlier_weight); get(l.fragment_start_intensity); get(l.gap_open); get(l.gap_extend);
        std::uint64_t g = 0;
        get(g);
        l.max_gap_bp = g;
        l.provenance = "union kernel inputs";
    };
    Inputs x;
    std::uint64_t H = 0;
    get(H);
    for (std::uint64_t h = 0; h < H; ++h) {
        std::uint64_t n = 0, k = 0;
        get(n);
        x.haplotypes.sequences.push_back(all.substr(at, n));
        at += n;
        get(k);
        std::vector<std::uint32_t> seg(k);
        for (auto& v : seg) get(v);
        x.segments.push_back(seg);
        x.haplotypes.index_of.push_back(h);
    }
    std::uint64_t P = 0;
    get(P);
    for (std::uint64_t i = 0; i < P; ++i) {
        std::uint64_t a = 0, b = 0;
        double p = 0, q = 0;
        get(a); get(b); get(p); get(q);
        x.pool.push_back({a, b});
        x.proposal.push_back(p);
        x.prior.push_back(q);
    }
    get(x.early_first);
    get(x.early_second);
    get_library(x.library);
    std::uint64_t v[24];
    for (auto& y : v) get(y);
    panvar::FragmentScoringOptions& o = x.options;
    o.threads = v[0]; o.mismatch_band = v[1]; o.minimum_piece_bp = v[2]; o.max_piece_hits = v[3]; o.piece_hit_publish_batch = v[4];
    o.max_proposed_starts_per_haplotype = v[5]; o.max_physical_origins = v[6]; o.max_state_classes = v[7]; o.max_mass_cells = v[8];
    o.max_latent_alignment_pairs = v[9]; o.emission_cache_bytes = v[10]; o.max_union_cached_placements_consumed = v[11];
    o.union_origin_cache_bytes = v[12]; o.max_union_cached_origins_consumed = v[13]; o.max_union_cached_states_consumed = v[14];
    o.union_prefilter_cache_bytes = v[15]; o.max_union_cached_prefilters_consumed = v[16]; o.start_sort_scratch_bytes = v[17];
    o.prepared_library_terms = v[18] != 0;
    o.gap_recursion_backend = v[19] == 0 ? panvar::GapRecursionBackend::scaled_reference : panvar::GapRecursionBackend::fixed_geometry_experimental;
    x.store_bytes = v[20]; x.legacy = v[21] != 0; x.views = v[22] != 0;
    for (std::uint64_t i = 0; i < v[23]; ++i) {
        panvar::LibraryModel l;
        get_library(l);
        o.audit_libraries.push_back(l);
    }
    return x;
}
std::vector<panvar::ReadPair> load_pairs(const std::string& r1, const std::string& r2) {
    const auto load = [](const std::string& file) {
        std::vector<panvar::SequencingRead> out;
        std::ifstream in(file);
        if (!in) throw std::runtime_error("cannot read " + file);
        for (std::string head, bases, plus, quals; std::getline(in, head) && std::getline(in, bases) && std::getline(in, plus) &&
                                                   std::getline(in, quals);) {
            panvar::SequencingRead r;
            const std::size_t sp = head.find_first_of(" \t");
            r.name = head.substr(1, sp == std::string::npos ? std::string::npos : sp - 1);
            for (char& c : bases) c = static_cast<char>(std::toupper(c));
            r.bases = bases;
            r.qualities = quals;
            out.push_back(std::move(r));
        }
        return out;
    };
    const auto a = load(r1), b = load(r2);
    std::vector<panvar::ReadPair> out;
    for (std::size_t i = 0; i < a.size() && i < b.size(); ++i) out.push_back({a[i], b[i]});
    return out;
}
void write_scores(const std::string& file, const std::vector<panvar::DiplotypeScore>& s, const std::vector<char>& have) {
    std::ofstream out(file, std::ios::binary);
    for (std::size_t i = 0; i < s.size(); ++i) {
        out.write(&have[i], 1);
        out.write(reinterpret_cast<const char*>(&s[i]), sizeof(panvar::DiplotypeScore));
    }
}
void read_scores(const std::string& file, std::vector<panvar::DiplotypeScore>& s, std::vector<char>& have, std::size_t n) {
    std::ifstream in(file, std::ios::binary);
    if (!in) throw std::runtime_error("cannot read " + file);
    s.resize(n);
    have.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        in.read(&have[i], 1);
        in.read(reinterpret_cast<char*>(&s[i]), sizeof(panvar::DiplotypeScore));
    }
}
struct Row {
    std::vector<std::pair<std::string, std::string>> kv;
    template <class T> void put(const std::string& k, const T& v) {
        std::ostringstream o;
        o << std::setprecision(10) << v;
        kv.push_back({k, o.str()});
    }
    void write(const std::string& file) const {
        std::ofstream out(file);
        out << "key\tvalue\n";
        for (const auto& [k, v] : kv) out << k << '\t' << v << '\n';
        for (const auto& [k, v] : kv) std::cout << k << '\t' << v << '\n';
    }
};
}  // namespace

int main(int argc, char** argv) {
    try {
        const std::string mode = argc > 1 ? argv[1] : "";
        if ((mode == "bulk" || mode == "progressive" || mode == "blocked") && argc >= 6) {
            const auto t0 = Clock::now();
            Inputs x = read_inputs(argv[2]);
            const auto fragments = load_pairs(argv[3], argv[4]);
            const std::string out = argv[5];
            Row row;
            row.put("mode", mode);
            row.put("haplotypes", x.haplotypes.sequences.size());
            row.put("pool", x.pool.size());
            row.put("fragments", fragments.size());
            row.put("threads", x.options.threads);
            row.put("fine_views", x.views ? 1 : 0);
            row.put("seconds_load", since(t0));
            row.put("footprint_after_load", panvar::physical_footprint().current);
            std::vector<panvar::DiplotypeScore> scores(x.pool.size());
            std::vector<char> have(x.pool.size(), 0);
            // S2-E: the footprint guard for bulk and blocked runs (the progressive mode keeps its own below)
            static std::mutex sample_mutex;
            static std::vector<std::string> sample_rows;
            static bool sample_done = false;
            std::thread guard_thread;
            const auto t_guard = Clock::now();
            if (mode == "bulk" || mode == "blocked") {
                const std::uint64_t limit = mode == "bulk" ? (argc > 6 ? std::stoull(argv[6]) : (std::uint64_t{4} << 30))
                                                           : (argc > 7 ? std::stoull(argv[7]) : (std::uint64_t{4} << 30));
                row.put("footprint_limit_bytes", limit);
                guard_thread = std::thread([&, limit]() {
                    for (;;) {
                        {
                            std::lock_guard<std::mutex> lock(sample_mutex);
                            if (sample_done) return;
                            const panvar::FootprintSample fp = panvar::physical_footprint();
                            if (fp.current > limit) {
                                std::ofstream f(out + ".tsv");
                                f << "key\tvalue\n";
                                for (const auto& [k, v] : row.kv) f << k << '\t' << v << '\n';
                                for (const std::string& line : sample_rows) f << line;
                                f << "completed\tno\nrefusal\tharness footprint guard: " << fp.current << " bytes > " << limit
                                  << "\nseconds_at_refusal\t" << since(t_guard) << "\nfootprint_peak\t" << fp.peak << '\n';
                                f.close();
                                std::cerr << "s2c_progressive_harness: footprint guard " << fp.current << " > " << limit << '\n';
                                std::_Exit(125);
                            }
                        }
                        std::this_thread::sleep_for(std::chrono::milliseconds(100));
                    }
                });
            }
            const auto stop_guard = [&]() {
                {
                    std::lock_guard<std::mutex> lock(sample_mutex);
                    sample_done = true;
                }
                if (guard_thread.joinable()) guard_thread.join();
            };
            if (mode == "blocked") {
                const std::size_t block = argc > 6 ? std::stoul(argv[6]) : 4096;
                row.put("block_fragments", block);
                const auto t = Clock::now();
                std::size_t number = 0;
                const panvar::BlockedUnionResult r = panvar::score_diplotypes_indexed_blocked(
                    x.haplotypes, x.pool, fragments, x.library, x.options, x.segments, x.store_bytes, x.views, block,
                    [&](panvar::BlockedBlockReport& rep, bool before) {
                        const panvar::FootprintSample fp = panvar::physical_footprint();
                        if (before) { rep.footprint_before = fp.current; return; }
                        rep.footprint_after = fp.current;
                        std::lock_guard<std::mutex> lock(sample_mutex);
                        const std::string k = "block" + std::to_string(++number) + "_";
                        std::ostringstream o;
                        o << k << "fragments\t" << rep.fragments << '\n' << k << "units_evaluated\t" << rep.units_evaluated << '\n'
                          << k << "store_peak_bytes\t" << rep.store_peak_bytes << '\n' << k << "origin_held_bytes\t" << rep.origin_held_bytes
                          << '\n' << k << "retained_column_bytes\t" << rep.retained_column_bytes << '\n' << k << "prefilter_held_bytes\t"
                          << rep.prefilter_held_bytes << '\n' << k << "emission_cache_bytes\t" << rep.emission_cache_bytes << '\n' << k
                          << "emission_cache_entries\t" << rep.emission_cache_entries << '\n' << k << "executed_dps\t" << rep.executed_dps
                          << '\n' << k << "direct_physical_origins\t" << rep.direct_physical_origins << '\n' << k
                          << "served_physical_origins\t" << rep.served_physical_origins << '\n' << k << "origins_constructed\t"
                          << rep.origins_constructed << '\n' << k << "origin_entities_admitted\t" << rep.origin_entities_admitted << '\n'
                          << k << "origin_entities_declined\t" << rep.origin_entities_declined << '\n' << k
                          << "cumulative_physical_origins\t" << rep.cumulative_physical_origins << '\n' << k << "footprint_before\t"
                          << rep.footprint_before << '\n' << k << "footprint_after\t" << rep.footprint_after << '\n' << k
                          << "footprint_peak_so_far\t" << fp.peak << '\n' << k << "seconds\t" << rep.seconds << '\n' << k
                          << "seconds_since_start\t" << since(t) << '\n';
                        sample_rows.push_back(o.str());
                    });
                stop_guard();
                row.put("seconds_total", since(t));
                row.put("completed", r.completed ? "yes" : "no");
                row.put("refusal", r.refusal.empty() ? "none" : r.refusal);
                row.put("seconds_index_build", r.seconds_index_build);
                row.put("seconds_combination", r.seconds_combination);
                row.put("blocks", r.blocks.size());
                row.put("accumulator_bytes", r.accumulator_bytes);
                row.put("index_unit_bytes", r.index_unit_bytes);
                row.put("executed_physical_origins", r.result.counters.physical_origins);
                row.put("executed_dps", r.result.emission_cache.executed_dps);
                row.put("executed_prefilters", r.result.emission_cache.executed_prefilters);
                row.put("origins_constructed", r.result.union_origin_cache.origins_constructed);
                row.put("origin_segments_admitted", r.result.union_origin_cache.admitted_segments);
                row.put("origin_segments_declined", r.result.union_origin_cache.declined_segments);
                row.put("full_read_verifications", r.result.counters.full_read_verifications);
                row.put("state_classes", r.result.counters.statistical_state_classes);
                for (const std::string& line : sample_rows) {
                    const std::size_t tab = line.find('\t');
                    (void)tab;
                }
                if (r.completed) {
                    scores = r.result.scores;
                    std::fill(have.begin(), have.end(), 1);
                }
                row.put("seconds_first_evaluated_candidate", r.completed ? since(t) : -1.0);
            } else if (mode == "bulk") {
                const auto t = Clock::now();
                panvar::IndexedUnionReport rep;
                std::string refusal;
                try {
                    const panvar::FragmentScoringResult r = panvar::score_diplotypes_indexed_union(
                        x.haplotypes, x.pool, fragments, x.library, x.options, x.segments, x.store_bytes, &rep, nullptr, x.legacy, x.views);
                    scores = r.scores;
                    std::fill(have.begin(), have.end(), 1);
                } catch (const panvar::ResourceExhausted& e) {
                    refusal = e.what();
                }
                row.put("seconds_total", since(t));
                row.put("seconds_first_evaluated_candidate", refusal.empty() ? since(t) : -1.0);
                row.put("refusal", refusal.empty() ? "none" : refusal);
                row.put("units_discovered", rep.batch.units_evaluated);
                row.put("index_units", rep.index.interior_units + rep.index.window_units + rep.index.fine_units);
                row.put("columns_computed", rep.batch.columns_computed_total);
                row.put("executed_prefilters", rep.batch.executed_prefilters);
                row.put("executed_dps", rep.batch.executed_dps);
                row.put("full_read_verifications", rep.batch.full_read_verifications);
                row.put("direct_physical_origins", rep.batch.direct_physical_origins);
                row.put("served_physical_origins", rep.batch.served_physical_origins);
                row.put("origins_constructed", rep.batch.origins_constructed);
                row.put("expanded_placements", rep.batch.expanded_placements);
                row.put("store_peak_bytes", rep.batch.store_peak_bytes);
                row.put("origin_held_bytes", rep.batch.origin_held_bytes);
                row.put("retained_column_bytes", rep.batch.retained_column_bytes);
                row.put("mass_matrix_bytes", static_cast<std::uint64_t>(x.pool.empty() ? 0 : fragments.size() * x.haplotypes.sequences.size() * 8));
                row.put("executed_physical_origins", rep.batch.counters_last.physical_origins);
                row.put("origin_segments_admitted", rep.batch.origin_entities_admitted);
                row.put("origin_segments_declined", rep.batch.origin_entities_declined);
                row.put("state_classes", rep.batch.counters_last.statistical_state_classes);
                stop_guard();
            } else {
                const std::size_t prefix = argc > 6 ? std::stoul(argv[6]) : 31;
                // the declared schedule
                std::vector<std::vector<panvar::HaplotypeCandidate>> stages(3);
                std::vector<std::vector<std::size_t>> stage_index(3);
                std::size_t early = x.pool.size();
                for (std::size_t c = 0; c < x.pool.size(); ++c) {
                    if (x.pool[c].first == x.early_first && x.pool[c].second == x.early_second) early = c;
                }
                if (early == x.pool.size()) throw std::runtime_error("the early pair is not in the pool");
                std::vector<std::size_t> order(x.pool.size());
                std::iota(order.begin(), order.end(), 0);
                std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return x.proposal[a] > x.proposal[b]; });
                stage_index[0].push_back(early);
                for (const std::size_t c : order) {
                    if (c == early) continue;
                    if (stage_index[1].size() < prefix) stage_index[1].push_back(c);
                }
                std::vector<char> placed(x.pool.size(), 0);
                for (std::size_t s = 0; s < 2; ++s) for (const std::size_t c : stage_index[s]) placed[c] = 1;
                for (std::size_t c = 0; c < x.pool.size(); ++c) if (!placed[c]) stage_index[2].push_back(c);
                for (std::size_t s = 0; s < 3; ++s) for (const std::size_t c : stage_index[s]) stages[s].push_back(x.pool[c]);
                row.put("prefix_size", prefix);
                const std::uint64_t limit = argc > 7 ? std::stoull(argv[7]) : (std::uint64_t{4} << 30);
                row.put("footprint_limit_bytes", limit);
                const bool persist = argc > 8 && std::string(argv[8]) == "persist";
                row.put("emission_cache_lifetime", persist ? "persist_across_discovery_calls" : "per_discovery_call");
                const auto t = Clock::now();
                // the footprint guard: completed stages are the checkpoint it keeps (no work after the limit)
                static std::mutex guard_mutex;
                static std::vector<std::string> done_rows;
                static bool finished = false;
                std::thread guard([&, limit]() {
                    for (;;) {
                        {
                            std::lock_guard<std::mutex> lock(guard_mutex);
                            if (finished) return;
                            const panvar::FootprintSample fp = panvar::physical_footprint();
                            if (fp.current > limit) {
                                std::ofstream f(out + ".tsv");
                                f << "key\tvalue\n";
                                for (const auto& [k, v] : row.kv) f << k << '\t' << v << '\n';
                                for (const std::string& line : done_rows) f << line;
                                f << "completed\tno (provisional: completed stages only)\nrefusal\tharness footprint guard: "
                                  << fp.current << " bytes > " << limit << "\nseconds_at_refusal\t" << since(t) << '\n';
                                f.close();
                                std::cerr << "s2c_progressive_harness: footprint guard " << fp.current << " > " << limit << '\n';
                                std::_Exit(125);
                            }
                        }
                        std::this_thread::sleep_for(std::chrono::milliseconds(100));
                    }
                });
                std::vector<panvar::ProgressiveStageReport> reports;
                const panvar::ProgressiveResult r = panvar::score_diplotypes_indexed_progressive(
                    x.haplotypes, stages, fragments, x.library, x.options, x.segments, x.store_bytes, x.views,
                    [&](panvar::ProgressiveStageReport& rep, const panvar::FragmentScoringResult&) {
                        const panvar::FootprintSample fp = panvar::physical_footprint();
                        rep.footprint_bytes = fp.current;
                        rep.footprint_peak_bytes = fp.peak;
                        std::lock_guard<std::mutex> lock(guard_mutex);
                        const std::string k = "stage" + std::to_string(done_rows.size() + 1) + "_";
                        std::ostringstream o;
                        o << std::setprecision(10) << k << "candidates\t" << rep.candidates << '\n' << k << "units_discovered\t"
                          << rep.units_discovered << '\n' << k << "columns_computed\t" << rep.columns_computed << '\n' << k
                          << "executed_dps\t" << rep.executed_dps << '\n' << k << "store_live_bytes\t" << rep.store_live_bytes << '\n'
                          << k << "origin_held_bytes\t" << rep.origin_held_bytes << '\n' << k << "footprint_bytes\t" << fp.current
                          << '\n' << k << "seconds_since_start\t" << rep.seconds_since_start << '\n';
                        done_rows.push_back(o.str());
                    }, persist);
                {
                    std::lock_guard<std::mutex> lock(guard_mutex);
                    finished = true;
                }
                guard.join();
                row.put("seconds_total", since(t));
                row.put("completed", r.completed ? "yes" : "no (provisional: completed stages only)");
                row.put("refusal", r.refusal.empty() ? "none" : r.refusal);
                row.put("seconds_index_build", r.seconds_index_build);
                row.put("seconds_preparation", r.seconds_preparation);
                row.put("index_units", r.index.interior_units + r.index.window_units + r.index.fine_units);
                for (std::size_t s = 0; s < r.stages.size(); ++s) {
                    const auto& st = r.stages[s];
                    const std::string k = "stage" + std::to_string(s + 1) + "_";
                    row.put(k + "candidates", st.candidates);
                    row.put(k + "slots", st.slots);
                    row.put(k + "units_discovered", st.units_discovered);
                    row.put(k + "units_reused", st.units_reused);
                    row.put(k + "columns_computed", st.columns_computed);
                    row.put(k + "columns_reused", st.columns_reused);
                    row.put(k + "pairs_computed", st.pairs_computed);
                    row.put(k + "executed_prefilters", st.executed_prefilters);
                    row.put(k + "executed_dps", st.executed_dps);
                    row.put(k + "full_read_verifications", st.full_read_verifications);
                    row.put(k + "direct_physical_origins", st.direct_physical_origins);
                    row.put(k + "served_physical_origins", st.served_physical_origins);
                    row.put(k + "origins_constructed", st.origins_constructed);
                    row.put(k + "expanded_placements", st.expanded_placements);
                    row.put(k + "store_live_bytes", st.store_live_bytes);
                    row.put(k + "store_peak_bytes", st.store_peak_bytes);
                    row.put(k + "origin_held_bytes", st.origin_held_bytes);
                    row.put(k + "retained_column_bytes", st.retained_column_bytes);
                    row.put(k + "worker_scratch_bytes", st.worker_scratch_bytes);
                    row.put(k + "footprint_bytes", st.footprint_bytes);
                    row.put(k + "footprint_peak_bytes", st.footprint_peak_bytes);
                    row.put(k + "seconds_discovery", st.seconds_discovery);
                    row.put(k + "seconds_scoring", st.seconds_scoring);
                    row.put(k + "seconds_since_start", st.seconds_since_start);
                    row.put(k + "dp_cache_hits", st.dp_cache_hits);
                    row.put(k + "dp_cache_misses", st.dp_cache_misses);
                    row.put(k + "carried_hits", st.carried_hits);
                    row.put(k + "carried_bytes", st.carried_bytes);
                    row.put(k + "carried_entries", st.carried_entries);
                    for (std::size_t i = 0; i < stage_index[s].size(); ++i) {
                        scores[stage_index[s][i]] = r.scores[s][i];
                        have[stage_index[s][i]] = 1;
                    }
                }
                row.put("seconds_first_evaluated_candidate", r.stages.empty() ? -1.0 : r.stages[0].seconds_since_start);
                row.put("dp_cache_hits", r.batch.emission_cache.dp_cache_hits);
                row.put("dp_cache_misses", r.batch.emission_cache.dp_cache_misses);
                row.put("persistent_adoptions", r.batch.persistent_adoptions);
                row.put("persistent_adoption_refusals", r.batch.persistent_adoption_refusals);
                row.put("persistent_carried_hits", r.batch.persistent_inherited_hits);
                row.put("units_discovered", r.batch.units_evaluated);
                row.put("columns_computed", r.batch.columns_computed_total);
                row.put("executed_prefilters", r.batch.executed_prefilters);
                row.put("executed_dps", r.batch.executed_dps);
                row.put("full_read_verifications", r.batch.full_read_verifications);
                row.put("direct_physical_origins", r.batch.direct_physical_origins);
                row.put("served_physical_origins", r.batch.served_physical_origins);
                row.put("origins_constructed", r.batch.origins_constructed);
                row.put("expanded_placements", r.batch.expanded_placements);
                row.put("store_peak_bytes", r.batch.store_peak_bytes);
                row.put("origin_held_bytes", r.batch.origin_held_bytes);
                row.put("retained_column_bytes", r.batch.retained_column_bytes);
            }
            row.put("pairs_scored", static_cast<std::uint64_t>(std::count(have.begin(), have.end(), 1)));
            row.put("footprint_end", panvar::physical_footprint().current);
            row.put("footprint_peak", panvar::physical_footprint().peak);
            row.put("peak_rss_bytes", peak_rss_bytes());
            write_scores(out + ".scores.bin", scores, have);
            row.write(out + ".tsv");
            if (mode == "blocked") {   // the per-block rows after the summary
                std::ofstream f(out + ".tsv", std::ios::app);
                for (const std::string& line : sample_rows) f << line;
            }
            return 0;
        }
        if (mode == "compare" && argc >= 6) {
            const Inputs x = read_inputs(argv[2]);
            std::vector<panvar::DiplotypeScore> a, b;
            std::vector<char> ha, hb;
            read_scores(std::string(argv[3]) + ".scores.bin", a, ha, x.pool.size());
            read_scores(std::string(argv[4]) + ".scores.bin", b, hb, x.pool.size());
            std::uint64_t both = 0, differ = 0;
            for (std::size_t i = 0; i < x.pool.size(); ++i) {
                if (!ha[i] || !hb[i]) continue;
                ++both;
                if (std::memcmp(&a[i].log_likelihood, &b[i].log_likelihood, 8) || std::memcmp(&a[i].fragment_term, &b[i].fragment_term, 8) ||
                    std::memcmp(&a[i].exposure_term, &b[i].exposure_term, 8) || std::memcmp(&a[i].absolute_term_sum, &b[i].absolute_term_sum, 8) ||
                    a[i].fragments_with_mass != b[i].fragments_with_mass) ++differ;
            }
            // the numerical top set by total score (likelihood + sequence prior): the exact maximum and every pair tied with it
            const auto top = [&](const std::vector<panvar::DiplotypeScore>& s, const std::vector<char>& h) {
                double best = -INFINITY;
                for (std::size_t i = 0; i < s.size(); ++i) if (h[i]) best = std::max(best, s[i].log_likelihood + x.prior[i]);
                std::set<std::size_t> out;
                for (std::size_t i = 0; i < s.size(); ++i) if (h[i] && s[i].log_likelihood + x.prior[i] == best) out.insert(i);
                return out;
            };
            Row row;
            row.put("pairs_scored_bulk", static_cast<std::uint64_t>(std::count(ha.begin(), ha.end(), 1)));
            row.put("pairs_scored_progressive", static_cast<std::uint64_t>(std::count(hb.begin(), hb.end(), 1)));
            row.put("pairs_compared", both);
            row.put("pairs_bitwise_different", differ);
            const auto ta = top(a, ha), tb = top(b, hb);
            row.put("top_set_bulk_size", ta.size());
            row.put("top_set_identical", ta == tb ? "yes" : "NO");
            row.write(argv[5]);
            return differ == 0 && ta == tb ? 0 : 1;
        }
        std::cerr << "usage: see the header comment\n";
        return 2;
    } catch (const std::exception& e) {
        std::cerr << "s2c_progressive_harness: " << e.what() << '\n';
        return 1;
    }
}
