// Gates for PreparedScoringContext: exact incremental haplotype scoring.
//
// The one-shot scorer is the ORACLE here. A context that scores the same haplotypes in batches must
// reproduce it cell for cell, not closely -- every arithmetic operation is the same operation in the
// same order, so these compare with ==, and a tolerance would only hide the day that stops being
// true. What batching changes is WHEN a column is computed, and a column depends on its own
// haplotype, the fragments and the library alone.
//
// The gates, and what each would catch:
//
//   1  one-shot against 2, 3 and 7 batches -- a column that depended on the haplotype count, the
//      restriding defect that row-major storage would have introduced.
//   2  batch-order invariance -- a column that depended on what was scored before it.
//   3  1/2/4-thread identity -- a reduction whose order moved with the partition.
//   4  cache keys after container growth -- the real point of the deque. Adding one sequence at a
//      time forces many growths; if a key pointed into relocatable storage, a later hit would read
//      moved bytes and the columns would silently diverge.
//   5  cumulative caps -- a cap applied per batch instead of to the whole search would accept a run
//      that exceeds it several times over.
//   6  teardown -- a cache left armed, or holding keys into a destroyed context's buffers, would be
//      read by the next pass.

#include "panvar/fragment_likelihood.hpp"
#include "panvar/genotype_model.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr double kNegInf = -std::numeric_limits<double>::infinity();

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

void require_identical(double observed, double expected, const std::string& message) {
    if (observed == expected) return;
    if (observed == kNegInf && expected == kNegInf) return;
    throw std::runtime_error(message + ": observed=" + std::to_string(observed) +
                             " expected=" + std::to_string(expected));
}

template <typename Exception>
void require_throws(const std::function<void()>& operation, const std::string& message) {
    try {
        operation();
    } catch (const Exception&) {
        return;
    } catch (...) {
        throw std::runtime_error(message + " (threw, but not the expected type)");
    }
    throw std::runtime_error(message);
}

std::string revcomp(const std::string& bases) {
    std::string out(bases.size(), 'N');
    for (std::size_t i = 0; i < bases.size(); ++i) {
        const char base = bases[bases.size() - 1 - i];
        out[i] = base == 'A' ? 'T' : base == 'C' ? 'G' : base == 'G' ? 'C' : 'A';
    }
    return out;
}

panvar::SequencingRead read_of(const std::string& name, const std::string& bases) {
    return {name, bases, std::string(bases.size(), 'I')};
}

panvar::LibraryModel library() {
    panvar::LibraryModel model;
    model.insert_mean_bp = 260.0;
    model.insert_sd_bp = 20.0;
    model.support_sigmas = 3.0;
    model.substitution_error = 0.02;
    model.outlier_weight = 0.005;
    model.fragment_start_intensity = 0.01;
    model.provenance = "fixture: declared";
    return model;
}

panvar::FragmentScoringOptions options(std::size_t threads = 1) {
    panvar::FragmentScoringOptions out;
    out.mismatch_band = 6;
    out.minimum_piece_bp = 10;
    out.threads = threads;
    // Exact publication, so a cap crosses where the totals say it does rather than up to a batch late.
    out.piece_hit_publish_batch = 1;
    return out;
}

std::string random_bases(std::mt19937_64& rng, std::size_t length) {
    static const char* alphabet = "ACGT";
    std::string out(length, 'A');
    for (char& base : out) base = alphabet[rng() % 4];
    return out;
}

// A panel that shares most of its sequence, which is what makes the emission cache do real work:
// distinct haplotypes with common windows are exactly the case cross-haplotype reuse exists for.
std::vector<std::string> panel(std::mt19937_64& rng, std::size_t count, std::size_t length) {
    const std::string backbone = random_bases(rng, length);
    std::vector<std::string> out;
    for (std::size_t i = 0; i < count; ++i) {
        std::string sequence = backbone;
        // Two small edited windows per haplotype, so they differ locally and agree elsewhere.
        for (int edit = 0; edit < 2; ++edit) {
            const std::size_t at = (rng() % (length - 40)) + 10;
            const std::string patch = random_bases(rng, 12);
            sequence.replace(at, patch.size(), patch);
        }
        out.push_back(sequence);
    }
    return out;
}

std::vector<panvar::ReadPair> fragments_from(
    std::mt19937_64& rng, const std::vector<std::string>& sources, std::size_t count,
    std::size_t read_length, std::size_t insert_bp) {
    std::vector<panvar::ReadPair> out;
    for (std::size_t i = 0; i < count; ++i) {
        const std::string& source = sources[rng() % sources.size()];
        const std::size_t start = rng() % (source.size() - insert_bp);
        const std::string forward = source.substr(start, read_length);
        const std::string reverse = revcomp(source.substr(start + insert_bp - read_length, read_length));
        out.push_back({read_of("f" + std::to_string(i) + "/1", forward),
                       read_of("f" + std::to_string(i) + "/2", reverse)});
    }
    return out;
}

panvar::InternedHaplotypes intern(const std::vector<std::string>& sequences) {
    std::vector<panvar::MosaicPath> paths;
    for (std::size_t i = 0; i < sequences.size(); ++i) paths.push_back(panvar::MosaicPath{{i}});
    return panvar::intern_haplotypes(paths, sequences);
}

std::vector<panvar::HaplotypeCandidate> all_pairs(std::size_t unique) {
    std::vector<panvar::HaplotypeCandidate> out;
    for (std::size_t a = 0; a < unique; ++a) {
        for (std::size_t b = a; b < unique; ++b) out.push_back({a, b});
    }
    return out;
}

// Split [0, n) into `batches` contiguous chunks, the way a repair search hands over a round at a time.
std::vector<std::vector<std::string>> split(const std::vector<std::string>& all, std::size_t batches) {
    std::vector<std::vector<std::string>> out(batches);
    for (std::size_t i = 0; i < all.size(); ++i) out[i * batches / all.size()].push_back(all[i]);
    return out;
}

struct Fixture {
    std::vector<std::string> haplotypes;
    std::vector<panvar::ReadPair> fragments;
};

Fixture fixture(std::uint64_t seed, std::size_t haplotype_count = 12, std::size_t fragment_count = 60) {
    std::mt19937_64 rng(seed);
    Fixture out;
    out.haplotypes = panel(rng, haplotype_count, 900);
    out.fragments = fragments_from(rng, out.haplotypes, fragment_count, 100, 260);
    return out;
}

void compare_to_one_shot(const Fixture& f, std::size_t batches, std::size_t threads) {
    const panvar::InternedHaplotypes interned = intern(f.haplotypes);
    const std::vector<panvar::HaplotypeCandidate> candidates = all_pairs(interned.sequences.size());
    // One-shot FIRST: while a context lives it owns the emission cache and this would refuse.
    const panvar::FragmentScoringResult oracle = panvar::score_diplotypes(
        interned, candidates, f.fragments, library(), options(threads));

    panvar::FragmentScoringResult batched;
    std::vector<std::size_t> mapped(interned.sequences.size(), 0);
    {
        panvar::PreparedScoringContext context(f.fragments, library(), options(threads));
        std::size_t next = 0;
        for (const std::vector<std::string>& group : split(interned.sequences, batches)) {
            panvar::IncrementalBatchReport report;
            const std::vector<std::size_t> slots = context.add_haplotypes(group, &report);
            for (const std::size_t slot : slots) mapped[next++] = slot;
            require(report.haplotypes_new + report.haplotypes_reused == report.haplotypes_offered,
                    "every offered haplotype is either new or reused");
        }
        require(next == interned.sequences.size(), "every haplotype was offered exactly once");
        require(context.unique_haplotypes() == interned.sequences.size(),
                "the context holds one column per unique haplotype");
        // Contiguous batches of an already-interned list preserve order, so the two index spaces agree.
        for (std::size_t h = 0; h < mapped.size(); ++h) {
            require(mapped[h] == h, "batching must not renumber the haplotypes");
        }
        for (std::size_t h = 0; h < interned.sequences.size(); ++h) {
            require_identical(context.exposure_of(h), oracle.exposure[h],
                              "exposure of haplotype " + std::to_string(h));
            for (std::size_t frag = 0; frag < f.fragments.size(); ++frag) {
                require_identical(context.mass_of(frag, h),
                                  oracle.mass_of(frag, h),
                                  "mass cell f" + std::to_string(frag) + " h" + std::to_string(h));
            }
        }
        batched = context.score_candidates(candidates);
    }

    for (std::size_t c = 0; c < candidates.size(); ++c) {
        const panvar::DiplotypeScore& a = batched.scores[c];
        const panvar::DiplotypeScore& b = oracle.scores[c];
        const std::string at = " at candidate " + std::to_string(c);
        require_identical(a.fragment_term, b.fragment_term, "fragment term" + at);
        require_identical(a.exposure_term, b.exposure_term, "exposure term" + at);
        require_identical(a.log_likelihood, b.log_likelihood, "total" + at);
        require_identical(a.absolute_term_sum, b.absolute_term_sum, "absolute term sum" + at);
        require(a.fragments_with_mass == b.fragments_with_mass, "fragments with mass" + at);
    }

    const panvar::CertifiedRanking ranked_batched = panvar::certify_scored_ranking(
        batched.scores, f.fragments.size(), options().mismatch_band, batched.counters,
        panvar::CertificationOptions{});
    const panvar::CertifiedRanking ranked_oracle = panvar::certify_scored_ranking(
        oracle.scores, f.fragments.size(), options().mismatch_band, oracle.counters,
        panvar::CertificationOptions{});
    require(ranked_batched.candidates.size() == ranked_oracle.candidates.size(),
            "the two rankings retain the same candidates");
    for (std::size_t i = 0; i < ranked_oracle.candidates.size(); ++i) {
        require(ranked_batched.candidates[i].candidate_id == ranked_oracle.candidates[i].candidate_id,
                "ranking order differs at position " + std::to_string(i));
        require_identical(ranked_batched.candidates[i].score, ranked_oracle.candidates[i].score,
                          "ranked score at position " + std::to_string(i));
    }
    require(ranked_batched.winner_set == ranked_oracle.winner_set, "the winner set differs");
    require(ranked_batched.status == ranked_oracle.status, "the certification status differs");
}

// -------------------------------------------------------------------------------------------
void test_batched_matches_one_shot() {
    const Fixture f = fixture(20260923);
    for (const std::size_t batches : {1u, 2u, 3u, 7u}) {
        compare_to_one_shot(f, batches, 1);
        std::cout << "  ok: " << batches << " batch(es) reproduce the one-shot scorer exactly\n";
    }
}

void test_batch_order_invariance() {
    const Fixture f = fixture(772211);
    const panvar::InternedHaplotypes interned = intern(f.haplotypes);

    // Score the same haplotypes in two different batch orders and compare COLUMNS BY SEQUENCE, which
    // is the claim: a haplotype's column may not depend on what was scored before it.
    const auto columns_for = [&](const std::vector<std::size_t>& order) {
        std::vector<std::vector<double>> out(interned.sequences.size());
        panvar::PreparedScoringContext context(f.fragments, library(), options(1));
        for (const std::size_t h : order) {
            const std::vector<std::size_t> slots = context.add_haplotypes({interned.sequences[h]});
            std::vector<double> column;
            for (std::size_t frag = 0; frag < f.fragments.size(); ++frag) {
                column.push_back(context.mass_of(frag, slots.front()));
            }
            out[h] = std::move(column);
        }
        return out;
    };

    std::vector<std::size_t> forward(interned.sequences.size());
    for (std::size_t i = 0; i < forward.size(); ++i) forward[i] = i;
    std::vector<std::size_t> shuffled = forward;
    std::mt19937_64 rng(99);
    std::shuffle(shuffled.begin(), shuffled.end(), rng);
    require(shuffled != forward, "the fixture must actually reorder, or nothing is tested");

    const std::vector<std::vector<double>> a = columns_for(forward);
    const std::vector<std::vector<double>> b = columns_for(shuffled);
    for (std::size_t h = 0; h < a.size(); ++h) {
        for (std::size_t frag = 0; frag < a[h].size(); ++frag) {
            require_identical(b[h][frag], a[h][frag],
                              "column of haplotype " + std::to_string(h) + " moved with batch order");
        }
    }
    std::cout << "  ok: a haplotype's column does not depend on the order batches arrived in\n";
}

void test_thread_identity() {
    const Fixture f = fixture(31337);
    for (const std::size_t threads : {2u, 4u}) {
        compare_to_one_shot(f, 3, threads);
        std::cout << "  ok: " << threads << " threads reproduce the serial result exactly\n";
    }
}

void test_cache_keys_survive_container_growth() {
    // One haplotype per batch, with the cache armed. A deque never relocates an element, so keys
    // stay valid; a vector would have moved the strings the keys point into many times over.
    const Fixture f = fixture(5150, 24);
    panvar::FragmentScoringOptions cached = options(1);
    cached.emission_cache_bytes = 1u << 20;

    const panvar::InternedHaplotypes interned = intern(f.haplotypes);
    const std::vector<panvar::HaplotypeCandidate> candidates = all_pairs(interned.sequences.size());
    const panvar::FragmentScoringResult oracle =
        panvar::score_diplotypes(interned, candidates, f.fragments, library(), cached);

    panvar::FragmentScoringResult batched;
    std::uint64_t hits = 0;
    {
        panvar::PreparedScoringContext context(f.fragments, library(), cached);
        for (const std::string& sequence : interned.sequences) {
            panvar::IncrementalBatchReport report;
            context.add_haplotypes({sequence}, &report);
            hits = report.cumulative_emission_cache.dp_cache_hits;
        }
        batched = context.score_candidates(candidates);
    }
    require(hits > 0, "the cache must actually serve hits across batches, or growth is untested");
    for (std::size_t c = 0; c < candidates.size(); ++c) {
        require_identical(batched.scores[c].log_likelihood, oracle.scores[c].log_likelihood,
                          "score after container growth at candidate " + std::to_string(c));
    }
    std::cout << "  ok: " << interned.sequences.size() << " separate batches, " << hits
              << " cross-batch cache hits, results unchanged\n";
}

void test_cumulative_caps() {
    const Fixture f = fixture(6060842);
    const panvar::InternedHaplotypes interned = intern(f.haplotypes);

    // What the whole search actually costs, measured rather than assumed.
    std::uint64_t total_piece_hits = 0;
    {
        panvar::PreparedScoringContext context(f.fragments, library(), options(1));
        for (const std::vector<std::string>& group : split(interned.sequences, 4)) {
            context.add_haplotypes(group);
        }
        total_piece_hits = context.cumulative_counters().piece_hits;
    }
    require(total_piece_hits > 0, "the fixture must do real work");

    // Exactly at the cumulative cap: accepted.
    {
        panvar::FragmentScoringOptions at_cap = options(1);
        at_cap.max_piece_hits = total_piece_hits;
        panvar::PreparedScoringContext context(f.fragments, library(), at_cap);
        for (const std::vector<std::string>& group : split(interned.sequences, 4)) {
            context.add_haplotypes(group);
        }
        require(!context.refused(), "a search exactly at its cumulative cap must be accepted");
        require(context.cumulative_counters().piece_hits == total_piece_hits,
                "the cumulative total must reproduce");
    }
    // One below: refused. A per-batch cap would have let this through, because no single batch
    // reaches the total.
    {
        panvar::FragmentScoringOptions under = options(1);
        under.max_piece_hits = total_piece_hits - 1;
        panvar::PreparedScoringContext context(f.fragments, library(), under);
        bool refused = false;
        try {
            for (const std::vector<std::string>& group : split(interned.sequences, 4)) {
                context.add_haplotypes(group);
            }
        } catch (const panvar::ResourceExhausted&) {
            refused = true;
        }
        require(refused, "one below the cumulative cap must refuse");
        require(context.refused(), "the refusal must invalidate the context");
        // And the whole search is invalid, not merely this batch: nothing may be read back.
        require_throws<panvar::ResourceExhausted>(
            [&]() { context.score_candidates({{0, 0}}); },
            "a refused context must not produce a partial ranking");
        require_throws<panvar::ResourceExhausted>(
            [&]() { (void)context.mass_of(0, 0); },
            "a refused context must not hand back retained columns");
        require_throws<panvar::ResourceExhausted>(
            [&]() { (void)context.add_haplotypes({interned.sequences.front()}); },
            "a refused context must not accept more work");
    }
    std::cout << "  ok: cumulative cap " << total_piece_hits
              << " accepts, one below refuses and invalidates the search\n";
}

void test_teardown_leaves_no_cache() {
    const Fixture f = fixture(424242);
    const panvar::InternedHaplotypes interned = intern(f.haplotypes);
    const std::vector<panvar::HaplotypeCandidate> candidates = all_pairs(interned.sequences.size());
    panvar::FragmentScoringOptions cached = options(1);
    cached.emission_cache_bytes = 1u << 20;

    const panvar::FragmentScoringResult before =
        panvar::score_diplotypes(interned, candidates, f.fragments, library(), cached);

    {
        panvar::PreparedScoringContext context(f.fragments, library(), cached);
        context.add_haplotypes(interned.sequences);
        // While the context is alive it OWNS the cache: a one-shot pass would clear the table its
        // keys live in, so it must refuse rather than silently recompute.
        require_throws<std::runtime_error>(
            [&]() {
                (void)panvar::score_diplotypes(interned, candidates, f.fragments, library(), cached);
            },
            "a one-shot pass must refuse while a context holds the cache");
    }

    // After teardown the cache is unpinned AND empty. If a key into the destroyed context's buffers
    // had survived, this pass would read moved bytes; if the table had merely been left armed, its
    // entries would answer for windows that no longer exist.
    const panvar::FragmentScoringResult after =
        panvar::score_diplotypes(interned, candidates, f.fragments, library(), cached);
    for (std::size_t c = 0; c < candidates.size(); ++c) {
        require_identical(after.scores[c].log_likelihood, before.scores[c].log_likelihood,
                          "score after context teardown at candidate " + std::to_string(c));
    }

    // An exception thrown mid-search must tear down just as cleanly as a return.
    {
        panvar::FragmentScoringOptions under = options(1);
        under.emission_cache_bytes = 1u << 20;
        under.max_piece_hits = 1;
        try {
            panvar::PreparedScoringContext context(f.fragments, library(), under);
            context.add_haplotypes(interned.sequences);
            require(false, "the capped context should have refused");
        } catch (const panvar::ResourceExhausted&) {
        }
    }
    const panvar::FragmentScoringResult after_throw =
        panvar::score_diplotypes(interned, candidates, f.fragments, library(), cached);
    for (std::size_t c = 0; c < candidates.size(); ++c) {
        require_identical(after_throw.scores[c].log_likelihood, before.scores[c].log_likelihood,
                          "score after an exception teardown at candidate " + std::to_string(c));
    }
    std::cout << "  ok: the cache is owned, released and left empty on both return and throw\n";
}

void test_reuse_is_actually_reuse() {
    // Requirement 2 stated as a measurement: offering the same haplotypes again must place nothing.
    const Fixture f = fixture(8899);
    const panvar::InternedHaplotypes interned = intern(f.haplotypes);
    panvar::PreparedScoringContext context(f.fragments, library(), options(1));
    panvar::IncrementalBatchReport first;
    context.add_haplotypes(interned.sequences, &first);
    panvar::IncrementalBatchReport again;
    context.add_haplotypes(interned.sequences, &again);
    require(first.haplotypes_new == interned.sequences.size(), "the first batch scores everything");
    require(again.haplotypes_new == 0, "a repeated batch must score nothing");
    require(again.haplotypes_reused == interned.sequences.size(), "and must reuse everything");
    require(again.counters.piece_hits == 0, "a reused batch must perform no placement work");
    require(again.cumulative_counters.piece_hits == first.cumulative_counters.piece_hits,
            "and must not add to the cumulative totals");
    std::cout << "  ok: re-offering " << interned.sequences.size()
              << " haplotypes places nothing and reuses every column\n";
}

}  // namespace

int main() {
    try {
        std::cout << "== 1. batched scoring reproduces the one-shot scorer ==\n";
        test_batched_matches_one_shot();
        std::cout << "== 2. batch-order invariance ==\n";
        test_batch_order_invariance();
        std::cout << "== 3. thread identity ==\n";
        test_thread_identity();
        std::cout << "== 4. cache keys survive container growth ==\n";
        test_cache_keys_survive_container_growth();
        std::cout << "== 5. cumulative caps ==\n";
        test_cumulative_caps();
        std::cout << "== 6. teardown on return and on throw ==\n";
        test_teardown_leaves_no_cache();
        std::cout << "== 7. reuse is reuse ==\n";
        test_reuse_is_actually_reuse();
        std::cout << "\nscoring-context gates passed\n";
    } catch (const std::exception& error) {
        std::cerr << "\nSCORING-CONTEXT GATE FAILED: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
