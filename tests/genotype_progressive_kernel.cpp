// S2-C (EXPERIMENTAL): the demand-driven route of the indexed scorer against its bulk control, on deterministic
// fixtures with repeated copies, a distant copy, indels, empty and short segments, a mosaic junction, a sequence-identical
// labelling, a short haplotype and reads in both orientations. Every evaluated score must equal the bulk adapter's bit
// for bit; reordered and repeated requests must reuse work and change nothing; an incomplete column is refused; a cap
// crossed mid-schedule keeps the completed stages; 1 and 4 threads agree; the selective route exports no mass matrix.
#include "panvar/fragment_likelihood.hpp"
#include "panvar/genotype_model.hpp"

#include <algorithm>
#include <cstring>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using panvar::DiplotypeScore;
using panvar::HaplotypeCandidate;
using panvar::ReadPair;

void require(bool ok, const std::string& what) {
    if (!ok) throw std::runtime_error(what);
}
bool bits(double a, double b) { return std::memcmp(&a, &b, sizeof(double)) == 0; }
bool same(const DiplotypeScore& a, const DiplotypeScore& b) {
    return bits(a.log_likelihood, b.log_likelihood) && bits(a.fragment_term, b.fragment_term) &&
           bits(a.exposure_term, b.exposure_term) && bits(a.absolute_term_sum, b.absolute_term_sum) &&
           a.fragments_with_mass == b.fragments_with_mass;
}
std::mt19937_64 rng(20261018);
std::string random_bases(std::size_t n) {
    std::string s(n, 'A');
    for (char& c : s) c = "ACGT"[rng() & 3];
    return s;
}
std::string revcomp(const std::string& s) {
    std::string out(s.rbegin(), s.rend());
    for (char& c : out) c = c == 'A' ? 'T' : c == 'C' ? 'G' : c == 'G' ? 'C' : 'A';
    return out;
}

struct Fixture {
    panvar::InternedHaplotypes haplotypes;
    std::vector<std::vector<std::uint32_t>> segments;
    std::vector<ReadPair> fragments;
    panvar::LibraryModel library;
    std::vector<HaplotypeCandidate> pool;
};

Fixture build() {
    Fixture fx;
    const std::string a = random_bases(1800);
    std::string b = a;
    for (const std::size_t at : {200u, 650u, 1100u, 1500u}) b[at] = b[at] == 'A' ? 'C' : 'A';           // SNPs
    const std::string c = a.substr(0, 900) + a.substr(500, 300) + a.substr(900);                       // tandem copy
    const std::string d = a + random_bases(400) + a.substr(100, 350);                                 // distant copy
    std::string e = a;
    e.insert(700, "GATTA");                                                                            // insertion
    e.erase(1300, 3);                                                                                  // deletion
    const std::string f = a.substr(0, 1000) + b.substr(1000);                                          // mosaic junction
    const std::string g = a.substr(0, 600) + a.substr(680);                                            // an empty block joined across
    const std::string h = a.substr(400, 420);                                                          // short
    const std::vector<std::string> paths_bytes = {a, b, c, d, e, f, g, h, b};   // the last: a sequence-identical labelling of b
    std::vector<panvar::MosaicPath> paths;
    for (std::size_t i = 0; i < paths_bytes.size(); ++i) paths.push_back(panvar::MosaicPath{{i}});
    fx.haplotypes = panvar::intern_haplotypes(paths, paths_bytes);
    // segmentations: block-like cuts with zero-length and short segments, a cut at the junction
    for (const std::string& s : fx.haplotypes.sequences) {
        std::vector<std::uint32_t> seg;
        std::size_t left = s.size();
        const std::size_t cuts[] = {0, 300, 7, 0, 400, 13, 600};
        for (const std::size_t cut : cuts) {
            const std::size_t take = std::min(cut, left);
            seg.push_back(static_cast<std::uint32_t>(take));
            left -= take;
        }
        seg.push_back(static_cast<std::uint32_t>(left));
        fx.segments.push_back(seg);
    }
    // reads: 100 bp mates from every haplotype, both orientations, a few substitutions
    for (std::size_t i = 0; i < 90; ++i) {
        const std::string& src = fx.haplotypes.sequences[i % fx.haplotypes.sequences.size()];
        const std::size_t insert = 240 + rng() % 60;
        if (src.size() < insert + 1) continue;
        const std::size_t start = rng() % (src.size() - insert);
        std::string m1 = src.substr(start, 100), m2 = revcomp(src.substr(start + insert - 100, 100));
        if (i % 3 == 0) m1[37] = m1[37] == 'A' ? 'G' : 'A';
        ReadPair p;
        p.first.bases = i % 2 ? m2 : m1;
        p.second.bases = i % 2 ? m1 : m2;
        p.first.name = p.second.name = "r" + std::to_string(i);
        p.first.qualities = p.second.qualities = std::string(100, 'I');
        fx.fragments.push_back(p);
    }
    fx.library.insert_mean_bp = 270;
    fx.library.insert_sd_bp = 25;
    fx.library.support_sigmas = 4;
    fx.library.substitution_error = 0.01;
    fx.library.outlier_weight = 0.002;
    fx.library.fragment_start_intensity = 0.05;
    fx.library.provenance = "s2c fixture";
    for (std::size_t x = 0; x < fx.haplotypes.sequences.size(); ++x) {
        for (std::size_t y = x; y < fx.haplotypes.sequences.size(); ++y) fx.pool.push_back({x, y});
    }
    return fx;
}

panvar::FragmentScoringOptions options(std::size_t threads) {
    panvar::FragmentScoringOptions o;
    o.threads = threads;
    o.emission_cache_bytes = 1 << 24;
    o.union_origin_cache_bytes = 1 << 24;
    o.max_union_cached_origins_consumed = 2000000000ULL;
    o.max_union_cached_states_consumed = 2000000000ULL;
    o.max_union_cached_placements_consumed = 2000000000ULL;
    o.union_prefilter_cache_bytes = 1 << 22;
    o.max_union_cached_prefilters_consumed = 2000000000ULL;
    o.start_sort_scratch_bytes = 1 << 22;
    o.prepared_library_terms = true;
    return o;
}

std::vector<DiplotypeScore> flatten(const panvar::ProgressiveResult& r, const std::vector<std::vector<HaplotypeCandidate>>& stages,
                                    const std::vector<HaplotypeCandidate>& pool) {
    std::vector<DiplotypeScore> out(pool.size());
    std::vector<char> filled(pool.size(), 0);
    for (std::size_t s = 0; s < r.scores.size(); ++s) {
        for (std::size_t i = 0; i < stages[s].size(); ++i) {
            for (std::size_t c = 0; c < pool.size(); ++c) {
                if (pool[c].first == stages[s][i].first && pool[c].second == stages[s][i].second) {
                    out[c] = r.scores[s][i];
                    filled[c] = 1;
                }
            }
        }
    }
    for (const char f : filled) require(f, "progressive: every pool candidate scored");
    return out;
}
}  // namespace

int main() {
    try {
        const Fixture fx = build();
        const std::size_t store = std::size_t{1} << 30;
        for (const bool views : {true, false}) {
            const std::string tag = views ? "fine_views" : "reconstructed";
            const panvar::FragmentScoringResult bulk = panvar::score_diplotypes_indexed_union(
                fx.haplotypes, fx.pool, fx.fragments, fx.library, options(2), fx.segments, store, nullptr, nullptr, false, views);
            // the declared schedule shape: one candidate, a small prefix, the remainder
            std::vector<std::vector<HaplotypeCandidate>> stages(3);
            stages[0].push_back(fx.pool[5]);
            for (std::size_t c = 0; c < fx.pool.size(); ++c) {
                if (c == 5) continue;
                (stages[1].size() < 4 ? stages[1] : stages[2]).push_back(fx.pool[c]);
            }
            const panvar::ProgressiveResult prog = panvar::score_diplotypes_indexed_progressive(
                fx.haplotypes, stages, fx.fragments, fx.library, options(2), fx.segments, store, views);
            require(prog.completed, tag + ": the schedule completes");
            for (std::size_t i = 0; i < prog.stages.size(); ++i) {
                const auto& st = prog.stages[i];
                std::cout << "  " << tag << " stage " << i + 1 << ": candidates " << st.candidates << ", slots " << st.slots
                          << ", units discovered " << st.units_discovered << " (reused " << st.units_reused << "), columns "
                          << st.columns_computed << " (reused " << st.columns_reused << "), executed dps " << st.executed_dps
                          << ", direct origins " << st.direct_physical_origins << ", served " << st.served_physical_origins
                          << ", constructed " << st.origins_constructed << '\n';
            }
            std::cout << "  " << tag << " index units: interior " << prog.index.interior_units << ", window " << prog.index.window_units
                      << ", fine " << prog.index.fine_units << "; bulk served origins " << bulk.union_origin_cache.origins_consumed << '\n';
            const auto flat = flatten(prog, stages, fx.pool);
            for (std::size_t c = 0; c < fx.pool.size(); ++c) require(same(flat[c], bulk.scores[c]), tag + ": progressive == bulk, bit for bit");
            require(prog.stages[0].slots <= 2 && prog.stages[0].columns_computed <= 2, tag + ": stage 1 builds only its own columns");
            require(prog.stages[0].units_discovered < prog.index.interior_units + prog.index.window_units + prog.index.fine_units,
                    tag + ": stage 1 does not discover the whole index");
            // reordered and repeated: the pool reversed in one stage, then the first stage again
            std::vector<std::vector<HaplotypeCandidate>> other(2);
            other[0].assign(fx.pool.rbegin(), fx.pool.rend());
            other[1] = stages[0];
            const panvar::ProgressiveResult again = panvar::score_diplotypes_indexed_progressive(
                fx.haplotypes, other, fx.fragments, fx.library, options(2), fx.segments, store, views);
            require(again.completed, tag + ": reordered schedule completes");
            for (std::size_t i = 0; i < other[0].size(); ++i) require(same(again.scores[0][i], bulk.scores[fx.pool.size() - 1 - i]), tag + ": reordered == bulk");
            require(same(again.scores[1][0], bulk.scores[5]), tag + ": a repeated request returns the same bits");
            require(again.stages[1].units_discovered == 0 && again.stages[1].columns_computed == 0 && again.stages[1].pairs_computed == 0,
                    tag + ": a repeated request reuses every unit, column and pair");
            // duplicates: a homozygous pair and the same pair twice in one request
            std::vector<std::vector<HaplotypeCandidate>> dup = {{{3, 3}, {1, 4}, {1, 4}, {4, 1}}};
            const panvar::ProgressiveResult d = panvar::score_diplotypes_indexed_progressive(
                fx.haplotypes, dup, fx.fragments, fx.library, options(2), fx.segments, store, views);
            const auto index_of = [&](std::size_t x, std::size_t y) {
                for (std::size_t c = 0; c < fx.pool.size(); ++c) if (fx.pool[c].first == std::min(x, y) && fx.pool[c].second == std::max(x, y)) return c;
                throw std::runtime_error("no pool entry");
            };
            require(same(d.scores[0][0], bulk.scores[index_of(3, 3)]) && same(d.scores[0][1], bulk.scores[index_of(1, 4)]) &&
                        same(d.scores[0][2], bulk.scores[index_of(1, 4)]) && same(d.scores[0][3], bulk.scores[index_of(1, 4)]),
                    tag + ": homozygous and duplicate requests == bulk");
            // 1 vs 4 threads
            const panvar::ProgressiveResult t1 = panvar::score_diplotypes_indexed_progressive(
                fx.haplotypes, stages, fx.fragments, fx.library, options(1), fx.segments, store, views);
            const panvar::ProgressiveResult t4 = panvar::score_diplotypes_indexed_progressive(
                fx.haplotypes, stages, fx.fragments, fx.library, options(4), fx.segments, store, views);
            const auto f1 = flatten(t1, stages, fx.pool), f4 = flatten(t4, stages, fx.pool);
            for (std::size_t c = 0; c < fx.pool.size(); ++c) require(same(f1[c], f4[c]) && same(f1[c], flat[c]), tag + ": 1 == 4 threads");
            // a cap crossed mid-schedule: the completed stages keep their scores, nothing partial is scored
            panvar::FragmentScoringOptions capped = options(2);
            capped.max_physical_origins = prog.stages[0].direct_physical_origins + prog.stages[0].origins_constructed + 1;
            const panvar::ProgressiveResult cut = panvar::score_diplotypes_indexed_progressive(
                fx.haplotypes, stages, fx.fragments, fx.library, capped, fx.segments, store, views);
            std::cout << "  " << tag << " capped schedule: " << cut.scores.size() << " stage(s) kept; refusal: " << cut.refusal.substr(0, 60) << '\n';
            require(!cut.completed && !cut.refusal.empty() && !cut.scores.empty() && cut.scores.size() < stages.size(),
                    tag + ": a cap crossed mid-schedule stops it with the earlier stages kept (" + std::to_string(cut.scores.size()) + ")");
            for (std::size_t s = 0; s < cut.scores.size(); ++s) {
                for (std::size_t i = 0; i < stages[s].size(); ++i) {
                    require(same(cut.scores[s][i], flat[index_of(stages[s][i].first, stages[s][i].second)]), tag + ": kept stage scores exact");
                }
            }
        }
        // S2-D: bounded emission reuse across the batch's discovery calls -- identical scores (default budget, a tiny
        // budget that evicts/declines, forced hash collisions), real reuse, poisoning after a refusal, no armed cache left
        {
            const panvar::FragmentScoringResult bulk = panvar::score_diplotypes_indexed_union(
                fx.haplotypes, fx.pool, fx.fragments, fx.library, options(2), fx.segments, store, nullptr, nullptr, false, true);
            std::vector<std::vector<HaplotypeCandidate>> stages(3);
            stages[0].push_back(fx.pool[5]);
            for (std::size_t c = 0; c < fx.pool.size(); ++c) {
                if (c == 5) continue;
                (stages[1].size() < 4 ? stages[1] : stages[2]).push_back(fx.pool[c]);
            }
            for (const int variant : {0, 1, 2}) {
                panvar::FragmentScoringOptions o = options(2);
                if (variant == 1) o.emission_cache_bytes = 4096;   // tiny: entries declined or evicted
                if (variant == 2) panvar::set_emission_cache_hash_collapse(true);
                const panvar::ProgressiveResult off = panvar::score_diplotypes_indexed_progressive(
                    fx.haplotypes, stages, fx.fragments, fx.library, o, fx.segments, store, true, {}, false);
                const panvar::ProgressiveResult on = panvar::score_diplotypes_indexed_progressive(
                    fx.haplotypes, stages, fx.fragments, fx.library, o, fx.segments, store, true, {}, true);
                panvar::set_emission_cache_hash_collapse(false);
                const auto a = flatten(off, stages, fx.pool), b = flatten(on, stages, fx.pool);
                for (std::size_t c = 0; c < fx.pool.size(); ++c) {
                    require(same(a[c], b[c]) && same(b[c], bulk.scores[c]), "S2-D: persistent emissions keep every score (variant " + std::to_string(variant) + ")");
                }
                std::cout << "  S2-D variant " << variant << ": executed dps per-call " << off.batch.executed_dps << ", persistent "
                          << on.batch.executed_dps << "; adoptions " << on.batch.persistent_adoptions << ", carried hits "
                          << on.batch.persistent_inherited_hits << '\n';
                if (variant == 0) {
                    require(on.batch.persistent_adoptions == 2 && on.batch.persistent_inherited_hits > 0 &&
                                on.batch.executed_dps <= off.batch.executed_dps,
                            "S2-D: the carried cache is adopted by each later discovery and serves hits");
                }
                require(!panvar::retained_process_storage().emission_cache_armed, "S2-D: no emission cache is left armed after a schedule");
            }
            // a refused selective call poisons the batch, and the global cache is not left armed
            panvar::ContentPlacementIndex index(100, fx.library.max_gap_bp, true);
            for (std::size_t h = 0; h < fx.haplotypes.sequences.size(); ++h) index.add_haplotype(fx.haplotypes.sequences[h], fx.segments[h]);
            panvar::FragmentScoringOptions o = options(2);
            o.max_piece_hits = 1;   // refuses inside discovery
            panvar::ContentIndexBatch batch(index, fx.fragments, fx.library, o);
            batch.set_fine_views(true);
            batch.set_persistent_emission_cache(true);
            bool refused = false, poisoned = false;
            try {
                batch.discover_for({0, 1});
                (void)batch.score_selected({{0, 1}});   // the cap may bind at publication or at the acceptance check
            } catch (const panvar::ResourceExhausted&) {
                refused = true;
            }
            require(refused && !panvar::retained_process_storage().emission_cache_armed, "S2-D: a refusal leaves no armed cache");
            try { batch.discover_for({2}); } catch (const std::logic_error&) { poisoned = true; }
            require(poisoned, "S2-D: a refused batch cannot continue on incomplete evidence");
        }
        // an incomplete column is refused; the selective route exports no matrix; bulk and selective never mix
        {
            panvar::ContentPlacementIndex index(100, fx.library.max_gap_bp, true);
            for (std::size_t h = 0; h < fx.haplotypes.sequences.size(); ++h) index.add_haplotype(fx.haplotypes.sequences[h], fx.segments[h]);
            panvar::ContentIndexBatch batch(index, fx.fragments, fx.library, options(2));
            batch.set_fine_views(true);
            batch.discover_for({0});
            bool refused = false;
            try {
                (void)batch.score_selected({{0, 3}});
            } catch (const std::logic_error&) {
                refused = true;
            }
            require(refused, "a slot needing undiscovered evidence is refused, never scored as zero mass");
            const panvar::FragmentScoringResult r = batch.score_selected({{0, 0}});
            require(r.log_mass.empty() && r.exposure.empty() && r.scores.size() == 1, "the selective route exports no mass matrix");
            require(batch.columns_completed() == 1, "only the requested column was built");
            bool mixed = false;
            try {
                batch.discover();
            } catch (const std::logic_error&) {
                mixed = true;
            }
            require(mixed, "a selective batch refuses bulk discovery");
        }
        std::cout << "genotype progressive kernel: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "genotype progressive kernel: FAIL: " << error.what() << '\n';
        return 1;
    }
}
