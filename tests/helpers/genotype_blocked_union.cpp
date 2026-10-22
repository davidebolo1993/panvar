// S2-E (EXPERIMENTAL): fragment-blocked union scoring against the bulk indexed adapter. Every candidate's score
// components (and every audit view) must be bit-identical for every block size -- one fragment, awkward boundaries, a final
// short block, more than the whole run -- per-fragment masses must equal bulk's, mixed read lengths must be supported
// without the prefilter cache and refused with it exactly as bulk refuses, 1 and 4 threads must agree, and a cumulative
// cap reached in a later block must return no scores.
#include "panvar/fragment_likelihood.hpp"
#include "panvar/genotype_model.hpp"

#include <algorithm>
#include <cmath>
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

// mixed: 100 bp mates, except the second half of the fragments uses 80 bp (the first block of 7 stays uniform)
Fixture build(bool mixed) {
    Fixture fx;
    const std::string a = random_bases(1800);
    std::string b = a;
    for (const std::size_t at : {200u, 650u, 1100u, 1500u}) b[at] = b[at] == 'A' ? 'C' : 'A';
    const std::string c = a.substr(0, 900) + a.substr(500, 300) + a.substr(900);            // repeated copy
    const std::string d = a + random_bases(400) + a.substr(100, 350);                      // distant copy
    std::string e = a;
    e.insert(700, "GATTA");                                                                 // indels
    e.erase(1300, 3);
    const std::string f = a.substr(0, 1000) + b.substr(1000);                               // junction
    const std::string g = a.substr(0, 600) + a.substr(680);                                 // empty block joined across
    const std::string h = a.substr(400, 420);                                               // short
    const std::vector<std::string> bytes = {a, b, c, d, e, f, g, h, b};   // last: an identical labelling of b
    std::vector<panvar::MosaicPath> paths;
    for (std::size_t i = 0; i < bytes.size(); ++i) paths.push_back(panvar::MosaicPath{{i}});
    fx.haplotypes = panvar::intern_haplotypes(paths, bytes);
    for (const std::string& s : fx.haplotypes.sequences) {
        std::vector<std::uint32_t> seg;
        std::size_t left = s.size();
        for (const std::size_t cut : {0u, 300u, 7u, 0u, 400u, 13u, 600u}) {
            const std::size_t take = std::min<std::size_t>(cut, left);
            seg.push_back(static_cast<std::uint32_t>(take));
            left -= take;
        }
        seg.push_back(static_cast<std::uint32_t>(left));
        fx.segments.push_back(seg);
    }
    for (std::size_t i = 0; i < 90; ++i) {
        const std::string& src = fx.haplotypes.sequences[i % fx.haplotypes.sequences.size()];
        const std::size_t len = mixed && i >= 45 ? 80 : 100;
        const std::size_t insert = 240 + rng() % 60;
        if (src.size() < insert + 1) continue;
        const std::size_t start = rng() % (src.size() - insert);
        std::string m1 = src.substr(start, len), m2 = revcomp(src.substr(start + insert - len, len));
        if (i % 3 == 0) m1[37] = m1[37] == 'A' ? 'G' : 'A';
        ReadPair p;
        p.first.bases = i % 2 ? m2 : m1;
        p.second.bases = i % 2 ? m1 : m2;
        p.first.name = p.second.name = "r" + std::to_string(i);
        p.first.qualities = std::string(p.first.bases.size(), 'I');
        p.second.qualities = std::string(p.second.bases.size(), 'I');
        fx.fragments.push_back(p);
    }
    fx.library.insert_mean_bp = 270;
    fx.library.insert_sd_bp = 25;
    fx.library.support_sigmas = 4;
    fx.library.substitution_error = 0.01;
    fx.library.outlier_weight = 0.002;
    fx.library.fragment_start_intensity = 0.05;
    fx.library.provenance = "s2e fixture";
    for (std::size_t x = 0; x < fx.haplotypes.sequences.size(); ++x) {
        for (std::size_t y = x; y < fx.haplotypes.sequences.size(); ++y) fx.pool.push_back({x, y});
    }
    return fx;
}

panvar::FragmentScoringOptions options(std::size_t threads, bool prefilter) {
    panvar::FragmentScoringOptions o;
    o.threads = threads;
    o.emission_cache_bytes = 1 << 24;
    o.union_origin_cache_bytes = 1 << 24;
    o.max_union_cached_origins_consumed = 2000000000ULL;
    o.max_union_cached_states_consumed = 2000000000ULL;
    o.max_union_cached_placements_consumed = 2000000000ULL;
    if (prefilter) {
        o.union_prefilter_cache_bytes = 1 << 22;
        o.max_union_cached_prefilters_consumed = 2000000000ULL;
    }
    o.start_sort_scratch_bytes = 1 << 22;
    o.prepared_library_terms = true;
    return o;
}
}  // namespace

int main() {
    try {
        const std::size_t store = std::size_t{1} << 30;
        for (const bool mixed : {false, true}) {
            const Fixture fx = build(mixed);
            const std::size_t nf = fx.fragments.size();
            const std::string tag = mixed ? "mixed read lengths" : "uniform read lengths";
            for (const bool audits : {false, true}) {
                panvar::FragmentScoringOptions o = options(2, !mixed);
                if (audits) {
                    for (const double factor : {0.9, 1.1}) {
                        panvar::LibraryModel l = fx.library;
                        l.fragment_start_intensity *= factor;
                        o.audit_libraries.push_back(l);
                    }
                }
                const panvar::FragmentScoringResult bulk = panvar::score_diplotypes_indexed_union(
                    fx.haplotypes, fx.pool, fx.fragments, fx.library, o, fx.segments, store, nullptr, nullptr, false, true);
                for (const std::size_t block : {std::size_t{1}, std::size_t{7}, std::size_t{13}, nf - 1, nf, nf + 5}) {
                    const panvar::BlockedUnionResult r = panvar::score_diplotypes_indexed_blocked(
                        fx.haplotypes, fx.pool, fx.fragments, fx.library, o, fx.segments, store, true, block);
                    require(r.completed, tag + ": blocked completes (block " + std::to_string(block) + ")");
                    require(r.blocks.size() == (nf + block - 1) / block, tag + ": every block processed");
                    for (std::size_t c = 0; c < fx.pool.size(); ++c) {
                        require(same(r.result.scores[c], bulk.scores[c]),
                                tag + ": blocked == bulk bit for bit (block " + std::to_string(block) + ", candidate " + std::to_string(c) + ")");
                        for (std::size_t k = 0; k < bulk.audit_scores.size(); ++k) {
                            require(same(r.result.audit_scores[k][c], bulk.audit_scores[k][c]), tag + ": audit view identical");
                        }
                    }
                    require(r.result.log_mass.empty(), tag + ": no matrix exported");
                }
            }
            // per-fragment masses of each block's batch == bulk's matrix (block size 7; seeded cumulative charges)
            {
                panvar::FragmentScoringOptions o = options(2, !mixed);
                const panvar::FragmentScoringResult bulk = panvar::score_diplotypes_indexed_union(
                    fx.haplotypes, fx.pool, fx.fragments, fx.library, o, fx.segments, store, nullptr, nullptr, false, true);
                panvar::ContentPlacementIndex index(100, fx.library.max_gap_bp, true);
                for (std::size_t h = 0; h < fx.haplotypes.sequences.size(); ++h) index.add_haplotype(fx.haplotypes.sequences[h], fx.segments[h]);
                const std::size_t H = fx.haplotypes.sequences.size();
                for (std::size_t begin = 0; begin < nf; begin += 7) {
                    const std::size_t end = std::min(nf, begin + 7);
                    panvar::ContentIndexBatch batch(index, std::vector<ReadPair>(fx.fragments.begin() + begin, fx.fragments.begin() + end),
                                                    fx.library, o);
                    batch.set_fine_views(true);
                    batch.set_one_shot(true);
                    batch.set_placement_store_budget(store);
                    batch.discover();
                    batch.build_columns();
                    for (std::size_t h = 0; h < H; ++h) {
                        for (std::size_t i = 0; i < end - begin; ++i) {
                            require(bits(batch.mass_columns()[h][i], bulk.log_mass[(begin + i) * H + h]), tag + ": per-fragment mass identical");
                        }
                    }
                }
            }
            // 1 vs 4 threads
            {
                const auto r1 = panvar::score_diplotypes_indexed_blocked(fx.haplotypes, fx.pool, fx.fragments, fx.library,
                                                                         options(1, !mixed), fx.segments, store, true, 13);
                const auto r4 = panvar::score_diplotypes_indexed_blocked(fx.haplotypes, fx.pool, fx.fragments, fx.library,
                                                                         options(4, !mixed), fx.segments, store, true, 13);
                for (std::size_t c = 0; c < fx.pool.size(); ++c) require(same(r1.result.scores[c], r4.result.scores[c]), tag + ": 1 == 4 threads");
            }
            // a cumulative cap reached in a later block: no scores; bulk refuses the same cap
            {
                const auto ref = panvar::score_diplotypes_indexed_blocked(fx.haplotypes, fx.pool, fx.fragments, fx.library,
                                                                          options(2, !mixed), fx.segments, store, true, 13);
                panvar::FragmentScoringOptions capped = options(2, !mixed);
                capped.max_physical_origins = ref.blocks[1].cumulative_physical_origins;   // blocks 1-2 fit, a later one does not
                const auto cut = panvar::score_diplotypes_indexed_blocked(fx.haplotypes, fx.pool, fx.fragments, fx.library, capped,
                                                                          fx.segments, store, true, 13);
                require(!cut.completed && !cut.refusal.empty() && cut.result.scores.empty() && cut.blocks.size() >= 2 &&
                            cut.blocks.size() < ref.blocks.size(),
                        tag + ": a cumulative cap in a later block returns no scores (" + std::to_string(cut.blocks.size()) + " blocks done)");
                bool bulk_refused = false;
                try {
                    (void)panvar::score_diplotypes_indexed_union(fx.haplotypes, fx.pool, fx.fragments, fx.library, capped, fx.segments,
                                                                 store, nullptr, nullptr, false, true);
                } catch (const panvar::ResourceExhausted&) {
                    bulk_refused = true;
                }
                require(bulk_refused, tag + ": bulk refuses the same cumulative cap");
                std::cout << "  " << tag << ": cap " << capped.max_physical_origins << " refused in block " << cut.blocks.size() + 1
                          << " of " << ref.blocks.size() << "; " << cut.refusal.substr(0, 50) << '\n';
            }
        }
        // the prefilter cache over mixed read lengths: refused by bulk and by blocked, even though the first block of 7
        // is uniform
        {
            const Fixture fx = build(true);
            bool bulk = false, blocked = false;
            try {
                (void)panvar::score_diplotypes_indexed_union(fx.haplotypes, fx.pool, fx.fragments, fx.library, options(2, true),
                                                             fx.segments, std::size_t{1} << 30, nullptr, nullptr, false, true);
            } catch (const std::invalid_argument&) {
                bulk = true;
            }
            try {
                (void)panvar::score_diplotypes_indexed_blocked(fx.haplotypes, fx.pool, fx.fragments, fx.library, options(2, true),
                                                               fx.segments, std::size_t{1} << 30, true, 7);
            } catch (const std::invalid_argument&) {
                blocked = true;
            }
            require(bulk && blocked, "mixed read lengths under the prefilter cache: refused by both routes");
        }
        std::cout << "genotype blocked union: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "genotype blocked union: FAIL: " << error.what() << '\n';
        return 1;
    }
}
