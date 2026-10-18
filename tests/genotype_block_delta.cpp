// Stage 2 gates: the exact haploid block-delta scorer against the full scorer.
//
// The full scorer is the ORACLE. A delta reuses the StateKey census of every origin the replacement
// cannot touch and re-enumerates only a local window, so the claim is equality, not closeness: the
// reused keys carry the identical doubles, the recomputed ones come from the identical emission, and
// both halves reduce through the same std::map order. These compare with ==.
//
// What each gate would catch:
//   substitution / insertion / deletion / empty allele -- a suffix translation that is wrong for
//     length-changing replacements, or an affected test that assumes equal lengths.
//   first and last block -- the end-margin defect. log_gap_shape_mass saturates at read_length +
//     max_gap_bp, so a length change moves the normaliser of origins near an end even though they
//     read none of the replaced bases. Reuse there is only sound if those origins are called
//     affected and the window is widened to cover them.
//   boundary-adjacent fragments -- the half-open test. An origin ending exactly at block_begin, or
//     starting exactly at block_end, must be REUSED. Get the interval closed and it is dropped from
//     the reused half while the window half also declines it, and the mass silently falls.
//   tandem repeats -- identical StateKeys must collapse with multiplicity, not be counted as copies.
//   gapped reads -- a placement carries several consumed-span components; merging them before the
//     geometry is applied is the defect this scorer was corrected for once already.
//   random sweep -- every haploid mass cell and exposure, over intervals nobody chose by hand.

#include "panvar/fragment_likelihood.hpp"
#include "panvar/genotype_model.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <new>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

// Every heap allocation in this binary passes through here, so "no allocation on the hit path" is a
// counted fact rather than a reading of the code.
std::atomic<std::uint64_t>& allocation_count() {
    static std::atomic<std::uint64_t> count{0};
    return count;
}
void* operator new(std::size_t size) {
    allocation_count().fetch_add(1, std::memory_order_relaxed);
    void* memory = std::malloc(size == 0 ? 1 : size);
    if (memory == nullptr) throw std::bad_alloc();
    return memory;
}
void operator delete(void* memory) noexcept { std::free(memory); }
void operator delete(void* memory, std::size_t) noexcept { std::free(memory); }
void* operator new[](std::size_t size) {
    allocation_count().fetch_add(1, std::memory_order_relaxed);
    void* memory = std::malloc(size == 0 ? 1 : size);
    if (memory == nullptr) throw std::bad_alloc();
    return memory;
}
void operator delete[](void* memory) noexcept { std::free(memory); }
void operator delete[](void* memory, std::size_t) noexcept { std::free(memory); }

namespace {

constexpr double kNegInf = -std::numeric_limits<double>::infinity();

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

void require_identical(double observed, double expected, const std::string& message) {
    if (observed == expected) return;
    if (observed == kNegInf && expected == kNegInf) return;
    throw std::runtime_error(message + ": delta=" + std::to_string(observed) +
                             " full=" + std::to_string(expected) +
                             " difference=" + std::to_string(observed - expected));
}

template <typename Exception>
void require_throws(const std::function<void()>& operation, const std::string& message) {
    try {
        operation();
    } catch (const Exception&) {
        return;
    } catch (...) {
        throw std::runtime_error(message + " (threw the wrong type)");
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

panvar::FragmentScoringOptions options() {
    panvar::FragmentScoringOptions out;
    out.mismatch_band = 6;
    out.minimum_piece_bp = 10;
    out.threads = 1;
    out.piece_hit_publish_batch = 1;
    out.emission_cache_bytes = 1u << 20;
    return out;
}

std::string random_bases(std::mt19937_64& rng, std::size_t length) {
    static const char* alphabet = "ACGT";
    std::string out(length, 'A');
    for (char& base : out) base = alphabet[rng() % 4];
    return out;
}

panvar::InternedHaplotypes intern(const std::vector<std::string>& sequences) {
    std::vector<panvar::MosaicPath> paths;
    for (std::size_t i = 0; i < sequences.size(); ++i) paths.push_back(panvar::MosaicPath{{i}});
    return panvar::intern_haplotypes(paths, sequences);
}

// A fragment whose TEMPLATE is exactly [start, start + insert). Used both for ordinary coverage and,
// deliberately, to place templates flush against a block boundary.
panvar::ReadPair fragment_at(const std::string& source, std::size_t start, std::size_t insert,
                             std::size_t read_length, const std::string& name) {
    const std::string forward = source.substr(start, read_length);
    const std::string reverse = revcomp(source.substr(start + insert - read_length, read_length));
    return {read_of(name + "/1", forward), read_of(name + "/2", reverse)};
}

struct Case {
    std::string name;
    std::size_t begin = 0;
    std::size_t end = 0;
    std::string replacement;
};

// Scores every case both ways and compares every haploid mass cell and the exposure.
//
// The oracle runs FIRST and alone: a live context owns the emission cache, so a one-shot pass while
// one exists would refuse.
void compare(const std::string& base, const std::vector<panvar::ReadPair>& fragments,
             const std::vector<Case>& cases, const std::string& label) {
    std::vector<std::string> variants{base};
    for (const Case& item : cases) {
        std::string repaired = base;
        repaired.replace(item.begin, item.end - item.begin, item.replacement);
        require(repaired != base, label + "/" + item.name + ": the case must actually change the haplotype");
        variants.push_back(repaired);
    }
    const panvar::InternedHaplotypes interned = intern(variants);
    std::vector<panvar::HaplotypeCandidate> candidates;
    for (std::size_t i = 0; i < interned.sequences.size(); ++i) candidates.push_back({i, i});
    const panvar::FragmentScoringResult oracle =
        panvar::score_diplotypes(interned, candidates, fragments, library(), options());

    std::uint64_t reused = 0;
    std::uint64_t recomputed = 0;
    {
        panvar::PreparedScoringContext context(fragments, library(), options());
        context.add_haplotypes({base});
        context.begin_transient_round();
        context.retain_census(0);
        for (std::size_t c = 0; c < cases.size(); ++c) {
            panvar::BlockDeltaWork work;
            const std::size_t slot = context.add_repaired_haplotype(
                0, cases[c].begin, cases[c].end, cases[c].replacement, &work);
            reused += work.origins_reused;
            recomputed += work.origins_recomputed;
            const std::string where = label + "/" + cases[c].name;
            require(!work.fell_back_to_full, where + ": the delta path must actually be taken");
            require_identical(context.exposure_of(slot), oracle.exposure[interned.index_of[c + 1]],
                              where + " exposure");
            for (std::size_t f = 0; f < fragments.size(); ++f) {
                require_identical(context.mass_of(f, slot),
                                  oracle.mass_of(f, interned.index_of[c + 1]),
                                  where + " mass cell f" + std::to_string(f));
            }
        }
    }
    require(reused > 0, label + ": no origin was reused, so the delta proved nothing");
    require(recomputed > 0, label + ": no origin was recomputed, so the window proved nothing");
    std::cout << "  ok: " << label << " -- " << cases.size() << " cases, " << reused
              << " origins reused, " << recomputed << " recomputed, every cell identical\n";
}

// -------------------------------------------------------------------------------------------
struct Locus {
    std::string base;
    std::vector<panvar::ReadPair> fragments;
};

// Dense coverage plus templates deliberately flush against the block boundaries used below.
Locus locus(std::uint64_t seed, const std::vector<std::size_t>& flush_boundaries) {
    std::mt19937_64 rng(seed);
    Locus out;
    out.base = random_bases(rng, 900);
    const std::size_t read_length = 100;
    const std::size_t insert = 260;
    std::size_t n = 0;
    for (std::size_t start = 0; start + insert <= out.base.size(); start += 17) {
        out.fragments.push_back(fragment_at(out.base, start, insert, read_length,
                                            "f" + std::to_string(n++)));
    }
    // Half-open boundary cases: a template ENDING exactly at a boundary, and one STARTING exactly at
    // it. Both must be reused when that boundary is a block edge.
    for (const std::size_t at : flush_boundaries) {
        if (at >= insert) {
            out.fragments.push_back(fragment_at(out.base, at - insert, insert, read_length,
                                                "endsAt" + std::to_string(at)));
        }
        if (at + insert <= out.base.size()) {
            out.fragments.push_back(fragment_at(out.base, at, insert, read_length,
                                                "startsAt" + std::to_string(at)));
        }
    }
    return out;
}

// A branch start seeded warm -- the tied neighbour promoted by delta from its parent incumbent's census,
// census retained -- against the same start seeded cold (placed in full, census from full placement). Both
// the start and a neighbour the branch then deltas FROM that census must equal the one-shot oracle in every
// mass cell and in exposure, bit for bit, with and without the origin-cache path.
void test_delta_promoted_seed_matches_cold_seed() {
    const std::size_t b1 = 300;
    const std::size_t b2 = 600;
    const Locus l = locus(71, {b1, b1 + 30, b2, b2 + 30});
    std::mt19937_64 rng(9);
    const std::string r1 = random_bases(rng, 30);
    const std::string r2 = random_bases(rng, 45);
    std::string seed = l.base;
    seed.replace(b1, 30, r1);
    std::string next = seed;
    next.replace(b2, 30, r2);
    const panvar::InternedHaplotypes interned = intern({seed, next});
    const panvar::FragmentScoringResult oracle =
        panvar::score_diplotypes(interned, {{0, 0}, {1, 1}}, l.fragments, library(), options());
    const auto check = [&](panvar::PreparedScoringContext& context, std::size_t slot, std::size_t oracle_index,
                           const std::string& where) {
        require_identical(context.exposure_of(slot), oracle.exposure[oracle_index], where + " exposure");
        for (std::size_t f = 0; f < l.fragments.size(); ++f) {
            require_identical(context.mass_of(f, slot), oracle.mass_of(f, oracle_index),
                              where + " mass cell f" + std::to_string(f));
        }
    };
    for (const bool origin_path : {false, true}) {
        const std::string mode = origin_path ? "origin-cache path" : "stage-5 path";
        {
            panvar::PreparedScoringContext context(l.fragments, library(), options());
            if (origin_path) {
                context.enable_origin_cache_path(true);
                context.set_origin_cache_limits(1u << 24, 0);
            }
            context.add_haplotypes({l.base});
            context.begin_transient_round();
            context.retain_census(0);
            panvar::BlockDeltaWork promote;
            const std::size_t start = context.add_repaired_haplotype(0, b1, b1 + 30, r1, &promote, true);
            require(!promote.fell_back_to_full, mode + ": the promotion must be a delta");
            require(promote.origins_reused > 0 && promote.origins_recomputed > 0,
                    mode + ": the promotion must reuse and recompute origins");
            require(context.holds_census(start), mode + ": a promoted seed must hold a census");
            check(context, start, interned.index_of[0], mode + " warm seed");
            panvar::BlockDeltaWork work;
            const std::size_t onward = context.add_repaired_haplotype(start, b2, b2 + 30, r2, &work);
            require(!work.fell_back_to_full, mode + ": the onward repair must delta from the promoted census");
            check(context, onward, interned.index_of[1], mode + " warm onward");
        }
        {
            panvar::PreparedScoringContext context(l.fragments, library(), options());
            if (origin_path) {
                context.enable_origin_cache_path(true);
                context.set_origin_cache_limits(1u << 24, 0);
            }
            context.add_haplotypes({seed});
            context.begin_transient_round();
            context.retain_census(0);
            check(context, 0, interned.index_of[0], mode + " cold seed");
            panvar::BlockDeltaWork work;
            const std::size_t onward = context.add_repaired_haplotype(0, b2, b2 + 30, r2, &work);
            check(context, onward, interned.index_of[1], mode + " cold onward");
        }
    }
    std::cout << "  ok: a delta-promoted seed and its onward repair equal cold full seeding, every cell, both paths\n";
}

void test_replacement_kinds() {
    const std::size_t b = 400;
    const Locus l = locus(11, {b, b + 30});
    std::mt19937_64 rng(5);
    compare(l.base, l.fragments,
            {{"substitution_same_length", b, b + 30, random_bases(rng, 30)},
             {"insertion_longer", b, b + 30, random_bases(rng, 55)},
             {"deletion_shorter", b, b + 30, random_bases(rng, 8)},
             {"empty_allele", b, b + 30, ""}},
            "replacement kinds");
}

void test_block_positions() {
    const Locus l = locus(23, {0, 40, 860, 900});
    std::mt19937_64 rng(6);
    compare(l.base, l.fragments,
            {{"first_block", 0, 40, random_bases(rng, 40)},
             {"first_block_length_changing", 0, 40, random_bases(rng, 61)},
             {"last_block", 860, 900, random_bases(rng, 40)},
             {"last_block_length_changing", 860, 900, random_bases(rng, 19)},
             {"interior", 430, 470, random_bases(rng, 40)}},
            "block positions");
}

void test_adjacent_short_blocks() {
    const Locus l = locus(37, {300, 306, 312});
    std::mt19937_64 rng(7);
    compare(l.base, l.fragments,
            {{"short_left", 300, 306, random_bases(rng, 6)},
             {"short_right", 306, 312, random_bases(rng, 6)},
             {"short_left_grows", 300, 306, random_bases(rng, 14)},
             {"short_right_shrinks", 306, 312, random_bases(rng, 2)}},
            "adjacent short blocks");
}

void test_tandem_repeat_multiplicity() {
    // A tandem array makes many origins share an insert AND a latent mass, so they collapse into one
    // StateKey with multiplicity. If the delta counted places instead of merging counts, or merged
    // two genuinely different keys, the mass would move.
    std::mt19937_64 rng(41);
    const std::string unit = random_bases(rng, 60);
    std::string base = random_bases(rng, 240);
    for (int i = 0; i < 6; ++i) base += unit;
    base += random_bases(rng, 240);
    Locus l;
    l.base = base;
    std::size_t n = 0;
    for (std::size_t start = 0; start + 260 <= base.size(); start += 13) {
        l.fragments.push_back(fragment_at(base, start, 260, 100, "t" + std::to_string(n++)));
    }
    compare(l.base, l.fragments,
            {{"inside_the_array", 300, 360, random_bases(rng, 60)},
             {"array_shortened", 300, 360, random_bases(rng, 20)},
             {"array_extended", 300, 360, random_bases(rng, 95)}},
            "tandem repeat multiplicity");
}

void test_gapped_components() {
    // Reads carrying an indel relative to the haplotype place through the pair-HMM with several
    // terminal consumed spans, each implying a different geometry.
    std::mt19937_64 rng(53);
    Locus l;
    l.base = random_bases(rng, 900);
    std::size_t n = 0;
    for (std::size_t start = 0; start + 260 <= l.base.size(); start += 19) {
        panvar::ReadPair pair = fragment_at(l.base, start, 260, 100, "g" + std::to_string(n++));
        // A one- or two-base indel inside one mate, alternating which mate carries it so both
        // orientations exercise the gapped path.
        std::string bases = pair.first.bases;
        if (n % 2 == 0) {
            bases.erase(50, 1);
            bases += random_bases(rng, 1);
            pair.first = read_of(pair.first.name, bases);
        } else {
            std::string other = pair.second.bases;
            other.insert(40, random_bases(rng, 2));
            other.resize(100);
            pair.second = read_of(pair.second.name, other);
        }
        l.fragments.push_back(pair);
    }
    compare(l.base, l.fragments,
            {{"substitution", 420, 460, random_bases(rng, 40)},
             {"length_changing", 420, 460, random_bases(rng, 73)}},
            "gapped multi-component spans");
}

void test_random_sweep() {
    for (std::uint64_t seed = 100; seed < 106; ++seed) {
        std::mt19937_64 rng(seed);
        const std::size_t begin = 120 + (rng() % 600);
        const std::size_t span = 5 + (rng() % 80);
        const std::size_t end = std::min<std::size_t>(begin + span, 890);
        const std::size_t new_len = rng() % 90;
        const Locus l = locus(seed, {begin, end});
        compare(l.base, l.fragments,
                {{"sweep", begin, end, random_bases(rng, new_len)}},
                "random sweep seed " + std::to_string(seed));
    }
}

void test_preconditions_refuse() {
    const Locus l = locus(77, {400});
    std::mt19937_64 rng(9);
    const std::string replacement = random_bases(rng, 40);
    panvar::PreparedScoringContext context(l.fragments, library(), options());
    context.add_haplotypes({l.base});
    context.begin_transient_round();
    // No census: refusing is the point. Placing in full here would report a delta that never
    // happened, which is indistinguishable from a delta that does not work.
    require_throws<std::invalid_argument>(
        [&]() { (void)context.add_repaired_haplotype(0, 400, 440, replacement); },
        "a missing census must refuse rather than silently place in full");
    context.retain_census(0);
    panvar::BlockDeltaWork work;
    (void)context.add_repaired_haplotype(0, 400, 440, replacement, &work);
    require(work.origins_reused > 0, "the delta should run once its preconditions hold");
    std::cout << "  ok: a missing census refuses rather than placing in full\n";
}

void test_transient_destruction_then_cache_hits() {
    // The gate the freeze used to make impossible. Score some repairs, destroy them, then score a
    // haplotype that shares their windows: the entries they inserted must still be there and still
    // serve. With pointer-keyed entries this would read freed bytes; with owned keys it is simply a
    // hit. A rising hit count after destruction is the whole claim of stage 3.
    const Locus l = locus(64, {400});
    std::mt19937_64 rng(21);
    panvar::PreparedScoringContext context(l.fragments, library(), options());
    context.add_haplotypes({l.base});
    context.begin_transient_round();
    context.retain_census(0);

    std::vector<std::size_t> doomed;
    std::vector<std::string> replacements;
    for (int i = 0; i < 4; ++i) replacements.push_back(random_bases(rng, 40));
    for (const std::string& replacement : replacements) {
        doomed.push_back(context.add_repaired_haplotype(0, 400, 440, replacement));
    }
    const std::uint64_t hits_before = context.cumulative_emission_cache().dp_cache_hits;
    // Every one of them goes, along with its sequence bytes.
    context.release_transient({});
    require(context.transient_payload_bytes() == 0, "the transient haplotypes must really be gone");

    // Re-score the identical repairs. Their windows are byte-identical to the destroyed ones, so
    // these must be served from entries whose source buffers no longer exist.
    panvar::BlockDeltaWork work;
    for (const std::string& replacement : replacements) {
        (void)context.add_repaired_haplotype(0, 400, 440, replacement, &work);
    }
    const std::uint64_t hits_after = context.cumulative_emission_cache().dp_cache_hits;
    require(hits_after > hits_before,
            "entries produced by destroyed haplotypes must still serve hits");
    std::cout << "  ok: after destroying every transient haplotype, the cache served "
              << (hits_after - hits_before) << " further hits\n";
}

void test_transient_release() {
    const Locus l = locus(88, {400});
    std::mt19937_64 rng(10);
    panvar::PreparedScoringContext context(l.fragments, library(), options());
    context.add_haplotypes({l.base});
    context.begin_transient_round();
    context.retain_census(0);
    const std::uint64_t before = context.transient_payload_bytes();
    require(before == 0, "nothing is transient before the first repair");

    std::vector<std::size_t> slots;
    for (int i = 0; i < 6; ++i) {
        slots.push_back(context.add_repaired_haplotype(0, 400, 440, random_bases(rng, 40)));
    }
    const std::uint64_t peak = context.transient_payload_bytes();
    require(peak > 0, "six repaired haplotypes must occupy something");
    // Keep one winner; everything else goes back.
    context.release_transient({slots.front()});
    const std::uint64_t after = context.transient_payload_bytes();
    require(after < peak, "releasing the losers must return storage");
    context.release_transient({});
    require(context.transient_payload_bytes() == 0,
            "releasing every transient haplotype must return to the previous level exactly");
    std::cout << "  ok: transient storage peaks at " << peak << " bytes and returns to 0\n";
}

void test_scoring_after_release() {
    // The gate that was missing. Releasing and scoring were each covered, but never scoring AFTER a
    // release -- and a released slot keeps its index while losing its column, so the flattening loop
    // walked into an empty vector. A repair search hits this on its second round, every time.
    const Locus l = locus(91, {400});
    std::mt19937_64 rng(12);
    panvar::PreparedScoringContext context(l.fragments, library(), options());
    context.add_haplotypes({l.base});
    context.begin_transient_round();
    context.retain_census(0);

    std::vector<std::size_t> slots;
    for (int i = 0; i < 4; ++i) {
        slots.push_back(context.add_repaired_haplotype(0, 400, 440, random_bases(rng, 40)));
    }
    const std::size_t keep = slots.front();
    const panvar::FragmentScoringResult before = context.score_candidates({{0, keep}});
    context.release_transient({keep});

    // The survivor still scores, and identically: releasing the losers may not disturb it.
    const panvar::FragmentScoringResult after = context.score_candidates({{0, keep}});
    require_identical(after.scores[0].log_likelihood, before.scores[0].log_likelihood,
                      "the surviving candidate's score after releasing its rivals");
    for (std::size_t f = 0; f < l.fragments.size(); ++f) {
        require_identical(context.mass_of(f, keep), before.mass_of(f, keep),
                          "surviving mass cell f" + std::to_string(f));
    }
    // And a released haplotype must be refused, not read through.
    const std::size_t gone = slots.back();
    require_throws<std::invalid_argument>(
        [&]() { (void)context.score_candidates({{0, gone}}); },
        "scoring a released haplotype must be refused");
    std::cout << "  ok: scoring after a release survives, and released slots are refused\n";
}

// Every mass cell of a full delta run, so two configurations can be compared exactly.
std::vector<double> delta_cells(const Locus& l, const std::vector<std::string>& replacements) {
    std::vector<double> out;
    panvar::PreparedScoringContext context(l.fragments, library(), options());
    context.add_haplotypes({l.base});
    context.begin_transient_round();
    context.retain_census(0);
    for (const std::string& replacement : replacements) {
        const std::size_t slot = context.add_repaired_haplotype(0, 400, 440, replacement);
        for (std::size_t f = 0; f < l.fragments.size(); ++f) out.push_back(context.mass_of(f, slot));
    }
    return out;
}

void test_forced_hash_collisions() {
    // Collapsing every key's hash to zero puts every distinct key in one bucket of one shard, so
    // byte equality alone separates them. A key that decided identity by hash -- or that compared
    // only a prefix -- would start returning one entry's answer for another's question, and the
    // masses would move. They may not.
    const Locus l = locus(55, {400});
    std::mt19937_64 rng(31);
    std::vector<std::string> replacements;
    for (int i = 0; i < 5; ++i) replacements.push_back(random_bases(rng, 40));

    const std::vector<double> ordinary = delta_cells(l, replacements);
    panvar::set_emission_cache_hash_collapse(true);
    std::vector<double> collided;
    try {
        collided = delta_cells(l, replacements);
    } catch (...) {
        panvar::set_emission_cache_hash_collapse(false);
        throw;
    }
    panvar::set_emission_cache_hash_collapse(false);

    require(ordinary.size() == collided.size(), "both runs must produce the same number of cells");
    for (std::size_t i = 0; i < ordinary.size(); ++i) {
        require_identical(collided[i], ordinary[i],
                          "mass cell " + std::to_string(i) + " under forced hash collision");
    }
    std::cout << "  ok: every key colliding into one bucket changes nothing across "
              << ordinary.size() << " cells\n";
}

void test_ambiguous_bases_refuse() {
    // The emission is defined over ACGT and require_acgt_window refuses anything else BEFORE a key is
    // formed, so an ambiguous base can never reach the cache. The identity question is therefore
    // settled by refusal, and that is what this asserts -- in both the read and the haplotype. The
    // owned key copies raw bytes and would admit a wider alphabet unchanged if the model ever did.
    std::mt19937_64 rng(77);
    const std::string clean = random_bases(rng, 600);
    std::string ambiguous = clean;
    ambiguous[300] = 'N';

    const panvar::ReadPair good = fragment_at(clean, 100, 260, 100, "clean");
    panvar::ReadPair ambiguous_read = good;
    std::string bases = ambiguous_read.first.bases;
    bases[40] = 'N';
    ambiguous_read.first = read_of(ambiguous_read.first.name, bases);

    require_throws<std::invalid_argument>(
        [&]() {
            (void)panvar::gapped_emission(clean, ambiguous_read.first.bases, 100, library(), 6);
        },
        "an ambiguous base in the READ must refuse");
    require_throws<std::invalid_argument>(
        [&]() { (void)panvar::gapped_emission(ambiguous, good.first.bases, 250, library(), 6); },
        "an ambiguous base in the HAPLOTYPE window must refuse");
    std::cout << "  ok: ambiguous bases refuse in both the read and the window, before any key\n";
}

void test_identical_after_transient_destruction() {
    // Stronger than "a hit still happens": the ANSWERS served from entries whose source buffers were
    // destroyed must be the same answers. A key that still referenced freed bytes could compare equal
    // by luck and return a neighbouring window's fit.
    const Locus l = locus(66, {400});
    std::mt19937_64 rng(41);
    std::vector<std::string> replacements;
    for (int i = 0; i < 4; ++i) replacements.push_back(random_bases(rng, 40));

    panvar::PreparedScoringContext context(l.fragments, library(), options());
    context.add_haplotypes({l.base});
    context.begin_transient_round();
    context.retain_census(0);

    std::vector<std::vector<double>> before;
    for (const std::string& replacement : replacements) {
        const std::size_t slot = context.add_repaired_haplotype(0, 400, 440, replacement);
        std::vector<double> column;
        for (std::size_t f = 0; f < l.fragments.size(); ++f) column.push_back(context.mass_of(f, slot));
        before.push_back(std::move(column));
    }
    context.release_transient({});
    require(context.transient_payload_bytes() == 0, "every transient haplotype must be gone");

    for (std::size_t i = 0; i < replacements.size(); ++i) {
        const std::size_t slot = context.add_repaired_haplotype(0, 400, 440, replacements[i]);
        for (std::size_t f = 0; f < l.fragments.size(); ++f) {
            require_identical(context.mass_of(f, slot), before[i][f],
                              "cell f" + std::to_string(f) + " of repair " + std::to_string(i) +
                                  " after its source buffers were destroyed");
        }
    }
    std::cout << "  ok: answers served after destruction are identical, not merely present\n";
}

void test_segmented_windows_match_materialised() {
    // The segmented extraction against the obvious materialised answer, for every kind of
    // replacement and at both ends of the locus. A segment walk that mishandled the suffix mapping
    // after a length change, or reversed the wrong range, would show here and nowhere else -- mass
    // equality can hide it whenever the affected set happens not to reach the discrepancy.
    std::mt19937_64 rng(4242);
    const std::string base = random_bases(rng, 900);
    struct Case { std::string name; std::size_t begin, end; std::string replacement; };
    const std::vector<Case> cases = {
        {"substitution", 400, 440, random_bases(rng, 40)},
        {"insertion", 400, 440, random_bases(rng, 77)},
        {"deletion", 400, 440, random_bases(rng, 9)},
        {"empty_allele", 400, 440, ""},
        {"first_block", 0, 40, random_bases(rng, 55)},
        {"first_block_empty", 0, 40, ""},
        {"last_block", 860, 900, random_bases(rng, 12)},
        {"last_block_empty", 860, 900, ""},
        {"whole_locus", 0, 900, random_bases(rng, 300)},
    };
    std::size_t checked = 0;
    for (const Case& item : cases) {
        std::string materialised = base;
        materialised.replace(item.begin, item.end - item.begin, item.replacement);
        require(panvar::segmented_length(base, item.begin, item.end, item.replacement) ==
                    materialised.size(),
                item.name + ": segmented length must match");
        require(panvar::segmented_hash(base, item.begin, item.end, item.replacement) ==
                    panvar::segmented_hash(materialised, 0, 0, ""),
                item.name + ": segmented hash must match the materialised hash");

        const std::size_t n = materialised.size();
        std::vector<std::pair<std::size_t, std::size_t>> ranges = {
            {0, n}, {0, 0}, {0, std::min<std::size_t>(n, 50)},
            {n > 50 ? n - 50 : 0, n},
        };
        if (n > item.begin) ranges.push_back({item.begin, n});
        if (item.begin > 10) ranges.push_back({item.begin - 10, std::min(n, item.begin + 10)});
        if (n > 200) ranges.push_back({100, 200});
        for (const auto& range : ranges) {
            if (range.second > n || range.first > range.second) continue;
            std::string forward;
            std::string backward;
            panvar::segmented_window(base, item.begin, item.end, item.replacement, range.first,
                                     range.second, forward, backward);
            const std::string expect = materialised.substr(range.first, range.second - range.first);
            require(forward == expect,
                    item.name + ": forward window [" + std::to_string(range.first) + "," +
                        std::to_string(range.second) + ")");
            require(backward == std::string(expect.rbegin(), expect.rend()),
                    item.name + ": reversed window [" + std::to_string(range.first) + "," +
                        std::to_string(range.second) + ")");
            ++checked;
        }
    }
    std::cout << "  ok: " << checked << " segmented windows match their materialised equivalents\n";
}

void test_path_distinct_haplotypes_share_one_column() {
    // Two DIFFERENT replacements that spell the same bytes. Held as segments they are distinct
    // objects with distinct block bounds, so only an exact comparison after the hash match can see
    // that they are one haplotype -- and they must share one column, or the same haplotype would be
    // placed twice and counted twice.
    const Locus l = locus(31337, {400});
    std::mt19937_64 rng(13);
    const std::string replacement = random_bases(rng, 40);

    panvar::PreparedScoringContext context(l.fragments, library(), options());
    context.add_haplotypes({l.base});
    context.begin_transient_round();
    context.retain_census(0);

    // [400,440) <- R, and [400,441) <- R + the base byte at 440: identical haplotypes, different paths.
    const std::size_t first = context.add_repaired_haplotype(0, 400, 440, replacement);
    const std::size_t second =
        context.add_repaired_haplotype(0, 400, 441, replacement + l.base[440]);
    require(first == second,
            "byte-identical haplotypes reached by different replacements must share one slot");

    // And a genuinely different one must NOT collapse onto it.
    const std::size_t other = context.add_repaired_haplotype(0, 400, 440, random_bases(rng, 40));
    require(other != first, "a different haplotype must not be deduplicated onto the same slot");
    std::cout << "  ok: path-distinct but byte-identical repairs share one column\n";
}

void test_cache_budget_exact_boundary() {
    // The budget governs retained REFUSAL entries, and now counts the bytes their keys own. At
    // exactly the measured figure everything is retained; one byte below, something is declined.
    const Locus l = locus(808, {400});
    const panvar::InternedHaplotypes interned = intern({l.base});
    const std::vector<panvar::HaplotypeCandidate> candidates = {{0, 0}};

    panvar::FragmentScoringOptions generous = options();
    generous.emission_cache_bytes = 1u << 30;
    (void)panvar::score_diplotypes(interned, candidates, l.fragments, library(), generous);
    const std::uint64_t used = panvar::last_emission_cache_counters().budgeted_bytes_used;
    // Non-zero because the budget now governs every retained entry, keys included. While it governed
    // only post-recursion refusals it was very nearly vacuous: the prefilter is exact, so almost
    // nothing reaches that stream, and this figure was 0 on any clean fixture.
    require(used > 0, "the fixture must retain some budgeted entries");

    panvar::FragmentScoringOptions at_bound = options();
    at_bound.emission_cache_bytes = static_cast<std::size_t>(used);
    (void)panvar::score_diplotypes(interned, candidates, l.fragments, library(), at_bound);
    const panvar::EmissionCacheCounters exact = panvar::last_emission_cache_counters();
    require(exact.entries_not_retained == 0,
            "exactly at the budget, every budgeted entry must be retained");

    panvar::FragmentScoringOptions below = options();
    below.emission_cache_bytes = static_cast<std::size_t>(used - 1);
    (void)panvar::score_diplotypes(interned, candidates, l.fragments, library(), below);
    const panvar::EmissionCacheCounters short_by_one = panvar::last_emission_cache_counters();
    require(short_by_one.entries_not_retained > 0,
            "one byte below the budget, at least one entry must be declined");
    std::cout << "  ok: budget " << used << " retains everything, " << (used - 1) << " declines "
              << short_by_one.entries_not_retained << "\n";
}

void test_no_dangling_reads_after_destruction() {
    // "No retained pointer addresses" stated as something observable: destroy every transient, then
    // churn the heap hard enough that freed pages are very likely reused, and require the cache to
    // return the SAME answers. A key still pointing at released bytes would now be comparing against
    // whatever overwrote them.
    const Locus l = locus(99991, {400});
    std::mt19937_64 rng(17);
    std::vector<std::string> replacements;
    for (int i = 0; i < 4; ++i) replacements.push_back(random_bases(rng, 40));

    panvar::PreparedScoringContext context(l.fragments, library(), options());
    context.add_haplotypes({l.base});
    context.begin_transient_round();
    context.retain_census(0);

    std::vector<std::vector<double>> before;
    for (const std::string& replacement : replacements) {
        const std::size_t slot = context.add_repaired_haplotype(0, 400, 440, replacement);
        std::vector<double> column;
        for (std::size_t f = 0; f < l.fragments.size(); ++f) column.push_back(context.mass_of(f, slot));
        before.push_back(std::move(column));
    }
    context.release_transient({});

    // Churn: a few hundred megabytes of short-lived allocations over the freed range.
    {
        std::vector<std::string> rubbish;
        for (int i = 0; i < 2000; ++i) rubbish.emplace_back(100000, static_cast<char>('A' + (i % 4)));
        rubbish.clear();
        rubbish.shrink_to_fit();
    }

    for (std::size_t i = 0; i < replacements.size(); ++i) {
        const std::size_t slot = context.add_repaired_haplotype(0, 400, 440, replacements[i]);
        for (std::size_t f = 0; f < l.fragments.size(); ++f) {
            require_identical(context.mass_of(f, slot), before[i][f],
                              "cell f" + std::to_string(f) + " of repair " + std::to_string(i) +
                                  " after the heap was churned over its released bytes");
        }
    }
    std::cout << "  ok: answers survive destruction followed by heavy heap reuse\n";
}

void test_arena_slabs_grow_and_compare_across_slabs() {
    // Slabs are 1 MB and append-only. Forcing every key into one shard (the hash collapse sends them
    // all to shard 0) makes a few thousand entries span several slabs, which is the only way to
    // exercise an entry whose bytes live in a slab other than the newest one. Every key must still
    // hit, and hit with the value it was stored with: a lookup that read from the wrong slab, or that
    // assumed the latest slab, would surface here and essentially nowhere else.
    // Reads are drawn FROM the haplotype at their own start, so each is an exact match that survives
    // the prefilter and is retained. Random reads would all be rejected before a key is even formed
    // -- prefilter refusals never reach the table -- and nothing would be stored to span a slab.
    std::mt19937_64 rng(20260923);
    const std::string haplotype = random_bases(rng, 20000);
    const std::size_t count = 6000;
    std::vector<std::size_t> starts;
    for (std::size_t i = 0; i < count; ++i) starts.push_back(i);

    panvar::set_emission_cache_hash_collapse(true);
    panvar::begin_emission_cache(1u << 30);
    std::vector<double> stored;
    try {
        for (const std::size_t start : starts) {
            const std::string read = haplotype.substr(start, 100);
            const auto fit = panvar::gapped_emission(haplotype, read, start, library(), 6);
            stored.push_back(fit.has_value() ? fit->log_emission : kNegInf);
        }
        const std::uint64_t arena = panvar::read_emission_cache_counters().content_arena_bytes;
        require(arena > (1u << 20),
                "the fixture must span more than one slab, or cross-slab lookup is untested");
        // Second pass: every one of these is a hit, and entries are spread across the slabs.
        for (std::size_t i = 0; i < starts.size(); ++i) {
            const std::string read = haplotype.substr(starts[i], 100);
            const auto fit = panvar::gapped_emission(haplotype, read, starts[i], library(), 6);
            const double again = fit.has_value() ? fit->log_emission : kNegInf;
            require_identical(again, stored[i],
                              "cross-slab hit " + std::to_string(i) + " returned a different value");
        }
        std::cout << "  ok: " << starts.size() << " entries across " << (arena >> 20)
                  << " slabs, every cross-slab lookup identical\n";
    } catch (...) {
        panvar::end_emission_cache();
        panvar::set_emission_cache_hash_collapse(false);
        throw;
    }
    panvar::end_emission_cache();
    panvar::set_emission_cache_hash_collapse(false);
}

void test_hits_allocate_nothing() {
    // The microbenchmark the requirement asks for, counted through the replaced operator new above.
    //
    // The property is that allocations do not GROW WITH HITS. A raw zero is not measurable here: any
    // loop that drives the cache carries some fixed cost of its own, and an earlier version of this
    // gate divided exactly that by the hit count and reported "0.0156 per hit" as if the hit path
    // were leaking. Two passes of very different lengths settle it -- if the count is the same for 64
    // hits and for 640, nothing in the hit path is allocating.
    // Each read is scored at ITS OWN start, so every call is a genuine accepted hit. Scoring them
    // all at start 0 -- as this did -- made 63 of 64 calls prefilter rejections, which never form a
    // key and never touch the table, so the "per hit" rate was measuring almost no hits at all.
    std::mt19937_64 rng(5150);
    const std::string haplotype = random_bases(rng, 2000);
    std::vector<std::size_t> starts;
    for (std::size_t i = 0; i < 64; ++i) starts.push_back(i * 7);
    std::vector<std::string> reads;
    for (const std::size_t start : starts) reads.push_back(haplotype.substr(start, 100));

    panvar::begin_emission_cache(1u << 28);
    std::uint64_t on_insert = 0;
    std::uint64_t on_hit = 0;
    std::uint64_t single_pass = 0;
    std::uint64_t wide = 0;
    std::uint64_t narrow = 0;
    try {
        const std::uint64_t before_insert = allocation_count().load();
        for (std::size_t i = 0; i < starts.size(); ++i) {
            (void)panvar::gapped_emission(haplotype, reads[i], starts[i], library(), 6);
        }
        on_insert = allocation_count().load() - before_insert;

        // Hoisted: constructing a LibraryModel per call is the test's own cost, not the cache's.
        const panvar::LibraryModel model = library();

        // ONE loop entry each, 64 hits against 640. Holding the number of passes fixed is the whole
        // point: an earlier version varied the pass count instead, so a fixed per-pass cost looked
        // like growth. Built up front so its own allocations are outside the counted region.
        std::vector<std::size_t> long_run;
        long_run.reserve(starts.size() * 10);
        for (int repeat = 0; repeat < 10; ++repeat) {
            for (const std::size_t start : starts) long_run.push_back(start);
        }

        const std::uint64_t before_single = allocation_count().load();
        for (std::size_t i = 0; i < starts.size(); ++i) {
            (void)panvar::gapped_emission(haplotype, reads[i], starts[i], model, 6);
        }
        single_pass = allocation_count().load() - before_single;

        const std::uint64_t before_wide = allocation_count().load();
        for (const std::size_t start : long_run) {
            (void)panvar::gapped_emission(haplotype, reads[start / 7], start, model, 6);
        }
        wide = allocation_count().load() - before_wide;
        narrow = 0;
        on_hit = wide > single_pass ? wide - single_pass : 0;
    } catch (...) {
        panvar::end_emission_cache();
        throw;
    }
    panvar::end_emission_cache();

    std::cout << "  measured: " << on_insert << " allocations over " << reads.size()
              << " insertions\n";
    std::cout << "  measured: one pass, " << reads.size() << " hits -> " << single_pass
              << " allocations\n";
    std::cout << "  measured: one pass, " << (reads.size() * 10) << " hits -> " << wide
              << " allocations\n";
    require(on_hit == 0,
            "allocations grew with hits: " + std::to_string(single_pass) + " at " +
                std::to_string(reads.size()) + " hits, " + std::to_string(wide) + " at " +
                std::to_string(reads.size() * 10) + " -- the hit path is allocating");
    std::cout << "  ok: ten times the hits, the same " << wide
              << " allocations -- nothing on the hit path scales with hits\n";
}

void test_window_enumeration_matches_reference() {
    // The old-versus-new gate. With the self-check armed, every delta runs BOTH enumerations and
    // refuses on any divergence in origins, their order, or the counters -- including
    // latent_alignment_pairs, which counts rejected pairs too, so an implementation reaching the same
    // origins by a different test sequence would be caught even though every mass matched.
    //
    // Driven through the real call path over the cases the suite already covers: substitutions,
    // insertions, deletions, empty alleles, first and last block with normaliser widening, adjacent
    // short blocks, tandem repeats (duplicate endpoints with several consumed-span classes), gapped
    // multi-component spans, and both orientations.
    panvar::set_window_enumeration_self_check(true);
    try {
        const std::size_t b = 400;
        const Locus l = locus(11, {b, b + 30});
        std::mt19937_64 rng(5);
        compare(l.base, l.fragments,
                {{"substitution", b, b + 30, random_bases(rng, 30)},
                 {"insertion", b, b + 30, random_bases(rng, 55)},
                 {"deletion", b, b + 30, random_bases(rng, 8)},
                 {"empty_allele", b, b + 30, ""}},
                "self-checked replacement kinds");

        const Locus ends = locus(23, {0, 40, 860, 900});
        compare(ends.base, ends.fragments,
                {{"first_block", 0, 40, random_bases(rng, 61)},
                 {"last_block", 860, 900, random_bases(rng, 19)}},
                "self-checked locus ends");

        // A tandem array: many endpoints repeat, so groups carry several rows apiece.
        std::mt19937_64 tandem(41);
        const std::string unit = random_bases(tandem, 60);
        std::string base = random_bases(tandem, 240);
        for (int i = 0; i < 6; ++i) base += unit;
        base += random_bases(tandem, 240);
        Locus repeats;
        repeats.base = base;
        std::size_t n = 0;
        for (std::size_t start = 0; start + 260 <= base.size(); start += 13) {
            repeats.fragments.push_back(fragment_at(base, start, 260, 100, "t" + std::to_string(n++)));
        }
        compare(repeats.base, repeats.fragments,
                {{"inside_the_array", 300, 360, random_bases(tandem, 60)},
                 {"array_shortened", 300, 360, random_bases(tandem, 20)}},
                "self-checked tandem repeats");
    } catch (...) {
        panvar::set_window_enumeration_self_check(false);
        throw;
    }
    panvar::set_window_enumeration_self_check(false);
    std::cout << "  ok: grouped enumeration matches the reference on every case, origins and counters\n";
}

void test_scratch_does_not_allocate_per_candidate() {
    // Persistent scratch means the per-form row vectors, the endpoint accumulators and the flat
    // origin output are built once and reused. After a warm-up candidate, scoring more candidates
    // must not allocate in proportion to them -- that is the whole point of the change, and the
    // reference allocated roughly 57,000 containers per candidate.
    const Locus l = locus(4242, {400});
    std::mt19937_64 rng(19);
    std::vector<std::string> replacements;
    for (int i = 0; i < 40; ++i) replacements.push_back(random_bases(rng, 40));

    panvar::PreparedScoringContext context(l.fragments, library(), options());
    context.add_haplotypes({l.base});
    context.begin_transient_round();
    context.retain_census(0);

    // Two properties, and they are different. CONSTRUCTION of per-form row containers must stop
    // outright after the first candidate -- that is the reference's 57,000-containers-per-candidate
    // defect, and it is measured directly. What remains is row GROWTH: each window sits at a
    // different locus position, so a different subset of read forms matches, and a row touched for
    // the first time at candidate five still has no capacity. That cost is a warm-up, not a per-form
    // cost, so the test is that it DECLINES -- not an absolute threshold, which would just be a
    // number chosen to make the result pass.
    panvar::BlockDeltaWork warm;
    (void)context.add_repaired_haplotype(0, 400, 440, replacements[0], &warm);

    panvar::BlockDeltaWork second;
    (void)context.add_repaired_haplotype(0, 400, 440, replacements[1], &second);

    const std::size_t batch = 15;
    const std::uint64_t before_early = allocation_count().load();
    for (std::size_t i = 2; i < 2 + batch; ++i) {
        panvar::BlockDeltaWork work;
        (void)context.add_repaired_haplotype(0, 400, 440, replacements[i], &work);
    }
    const double early = static_cast<double>(allocation_count().load() - before_early) / batch;

    std::uint64_t later_rows = 0;
    std::uint64_t later_accumulators = 0;
    std::uint64_t later_origins = 0;
    const std::uint64_t before_late = allocation_count().load();
    for (std::size_t i = 2 + batch; i < 2 + 2 * batch; ++i) {
        panvar::BlockDeltaWork work;
        (void)context.add_repaired_haplotype(0, 400, 440, replacements[i], &work);
        later_rows += work.placement_row_constructions;
        later_accumulators += work.accumulator_growths;
        later_origins += work.origin_output_growths;
    }
    const double late = static_cast<double>(allocation_count().load() - before_late) / batch;

    std::cout << "  measured: warm-up candidate built " << warm.placement_row_constructions
              << " placement rows; the next built " << second.placement_row_constructions << "\n";
    std::cout << "  measured: allocations per candidate, first " << batch << " -> " << early
              << ", next " << batch << " -> " << late << "\n";
    std::cout << "  measured: over " << batch << " warmed candidates -- row constructions "
              << later_rows << ", accumulator growths " << later_accumulators
              << ", origin-output growths " << later_origins << "\n";
    // What this batch is responsible for: the enumeration's own containers. Zero on all three means
    // no per-form row is built and neither the endpoint accumulators nor the flat origin output grow.
    require(second.placement_row_constructions == 0,
            "a warmed candidate must construct no placement rows, built " +
                std::to_string(second.placement_row_constructions));
    require(later_rows == 0,
            "later candidates must construct no placement rows, built " + std::to_string(later_rows));
    require(later_accumulators == 0,
            "endpoint accumulators must be reused, grew " + std::to_string(later_accumulators));
    require(later_origins == 0,
            "the flat origin output must be reused, grew " + std::to_string(later_origins));
    // The residual is REPORTED, not asserted. It does not fall with warm-up -- measured flat at about
    // 1170 a candidate -- so it is not row growth, and with the three counters above at zero it is
    // not the enumeration at all. The StateKey merge clears and refills a std::map per fragment,
    // allocating a node per state across every fragment, and that is a different call site than the
    // one this batch was asked to change.
    std::cout << "  reported: " << early << " then " << late
              << " total allocations per candidate, flat -- outside the enumeration (see the merge)\n";
    std::cout << "  ok: the enumeration builds no per-form or per-origin container after warm-up\n";
}

// ================================================================================================
// Stage 6: the window-origin cache and the preaggregated sorted merge.
// ================================================================================================

struct Repair {
    std::size_t begin = 0;
    std::size_t end = 0;
    std::string allele;
};

struct SearchResult {
    std::vector<double> cells;
    std::vector<double> exposures;
    std::vector<double> scores;
    std::uint64_t hits = 0, misses = 0, stored = 0, copied = 0, reduced = 0;
    std::uint64_t executed_origins = 0, consumed = 0, not_retained = 0, cache_bytes = 0;
    std::uint64_t hit_calls_with_work = 0;   // a hit that nevertheless charged executed work
};

struct SearchOptions {
    bool origin_path = false;
    std::size_t budget = 0;
    std::uint64_t record_cap = 0;
    std::size_t threads = 1;
    std::uint64_t max_physical_origins = 0;
};

// Round A scores every repair against the base; one same-length repair is promoted to incumbent and
// the rest released; round B scores the same repairs against the winner. Windows away from the
// promoted block are byte-identical in round B, which is what a cross-round cache is for.
SearchResult mini_search(const Locus& l, const std::vector<Repair>& repairs, std::size_t promote,
                         const SearchOptions& how) {
    panvar::FragmentScoringOptions opts = options();
    opts.threads = how.threads;
    if (how.max_physical_origins != 0) opts.max_physical_origins = how.max_physical_origins;
    panvar::PreparedScoringContext context(l.fragments, library(), opts);
    context.add_haplotypes({l.base});
    context.begin_transient_round();
    context.retain_census(0);
    if (how.origin_path) {
        context.enable_origin_cache_path(true);
        context.set_origin_cache_limits(how.budget, how.record_cap);
    }
    SearchResult out;
    const auto record = [&](std::size_t incumbent, std::size_t slot, const panvar::BlockDeltaWork& w) {
        for (std::size_t f = 0; f < l.fragments.size(); ++f) out.cells.push_back(context.mass_of(f, slot));
        out.exposures.push_back(context.exposure_of(slot));
        out.scores.push_back(context.score_candidates({{incumbent, slot}}).scores[0].log_likelihood);
        out.hits += w.origin_cache_hits;
        out.misses += w.origin_cache_misses;
        out.stored += w.origin_cache_stored;
        out.copied += w.fragments_copied;
        out.reduced += w.fragments_reduced;
        if (w.origin_cache_hits != 0 &&
            (w.piece_hits != 0 || w.dp_executions != 0 || w.physical_origins_constructed != 0)) {
            ++out.hit_calls_with_work;
        }
    };
    for (const Repair& r : repairs) {
        panvar::BlockDeltaWork w;
        const std::size_t slot = context.add_repaired_haplotype(0, r.begin, r.end, r.allele, &w);
        record(0, slot, w);
    }
    panvar::BlockDeltaWork promoted;
    const Repair& chosen = repairs[promote];
    const std::size_t winner =
        context.add_repaired_haplotype(0, chosen.begin, chosen.end, chosen.allele, &promoted, true);
    record(0, winner, promoted);
    context.release_transient({winner});
    for (std::size_t i = 0; i < repairs.size(); ++i) {
        if (i == promote) continue;
        panvar::BlockDeltaWork w;
        const std::size_t slot =
            context.add_repaired_haplotype(winner, repairs[i].begin, repairs[i].end, repairs[i].allele, &w);
        record(winner, slot, w);
    }
    out.executed_origins = context.cumulative_counters().physical_origins;
    out.consumed = context.origin_records_consumed();
    out.not_retained = context.origin_cache_not_retained();
    out.cache_bytes = context.origin_cache_bytes();
    return out;
}

void require_same_search(const SearchResult& got, const SearchResult& want, const std::string& label) {
    require(got.cells.size() == want.cells.size(), label + ": cell counts differ");
    for (std::size_t i = 0; i < want.cells.size(); ++i) {
        require_identical(got.cells[i], want.cells[i], label + ": mass cell " + std::to_string(i));
    }
    for (std::size_t i = 0; i < want.exposures.size(); ++i) {
        require_identical(got.exposures[i], want.exposures[i], label + ": exposure " + std::to_string(i));
    }
    for (std::size_t i = 0; i < want.scores.size(); ++i) {
        require_identical(got.scores[i], want.scores[i], label + ": score " + std::to_string(i));
    }
}

// Five blocks: one near the start (an at_start window), three interior, one near the end (at_end).
struct MultiBlock {
    Locus locus;
    std::vector<Repair> repairs;
    std::size_t promote = 0;   // a same-length substitution, so coordinates hold in round B
};

MultiBlock multi_block(std::uint64_t seed) {
    std::mt19937_64 rng(seed);
    MultiBlock out;
    out.locus.base = random_bases(rng, 6000);
    std::size_t n = 0;
    for (std::size_t start = 0; start + 260 <= out.locus.base.size(); start += 11) {
        out.locus.fragments.push_back(fragment_at(out.locus.base, start, 260, 100, "m" + std::to_string(n++)));
    }
    const std::vector<std::size_t> blocks = {50, 1500, 3000, 4500, 5900};
    for (const std::size_t b : blocks) {
        const std::size_t e = std::min<std::size_t>(b + 40, out.locus.base.size());
        out.repairs.push_back({b, e, random_bases(rng, e - b)});           // same length
        out.repairs.push_back({b, e, random_bases(rng, (e - b) + 23)});    // longer
        out.repairs.push_back({b, e, random_bases(rng, (e - b) / 3)});     // shorter
    }
    out.promote = 6;   // the same-length substitution at block 3000
    return out;
}

void test_origin_path_matches_stage5() {
    const MultiBlock m = multi_block(606);
    const SearchResult stage5 = mini_search(m.locus, m.repairs, m.promote, {});
    SearchOptions cached;
    cached.origin_path = true;
    cached.budget = std::size_t{1} << 30;
    panvar::set_origin_cache_self_check(true);
    SearchResult got;
    try {
        got = mini_search(m.locus, m.repairs, m.promote, cached);
    } catch (...) {
        panvar::set_origin_cache_self_check(false);
        throw;
    }
    panvar::set_origin_cache_self_check(false);
    require_same_search(got, stage5, "origin path vs stage 5");
    require(got.hits > 0, "round B must reuse windows, or the cache path was never exercised");
    require(got.copied > 0, "some fragments must be copied, or the untouched-fragment path is untested");
    require(got.hit_calls_with_work == 0, "a cache hit must charge no executed work");
    std::cout << "  ok: bit-identical to stage 5 over " << got.cells.size() << " cells; " << got.hits
              << " hits, " << got.misses << " misses; " << got.copied << " fragments copied, "
              << got.reduced << " reduced; self-check certified keys, multiplicities and census\n";
}

void test_origin_budgets_identical() {
    const MultiBlock m = multi_block(707);
    SearchOptions disabled;
    disabled.origin_path = true;
    disabled.budget = 0;
    SearchOptions full = disabled;
    full.budget = std::size_t{1} << 30;
    const SearchResult none = mini_search(m.locus, m.repairs, m.promote, disabled);
    const SearchResult all = mini_search(m.locus, m.repairs, m.promote, full);
    SearchOptions small = disabled;
    small.budget = static_cast<std::size_t>(all.cache_bytes / 3);
    const SearchResult some = mini_search(m.locus, m.repairs, m.promote, small);
    require_same_search(all, none, "fully retained vs disabled");
    require_same_search(some, none, "undersized vs disabled");
    require(none.hits == 0, "a disabled cache must not hit");
    require(all.hits > some.hits, "an undersized cache must hit less than a full one");
    require(some.not_retained > 0, "the undersized budget must actually decline entries");
    std::cout << "  ok: disabled (" << none.hits << " hits), undersized (" << some.hits << " hits, "
              << some.not_retained << " declined) and full (" << all.hits << " hits, " << all.cache_bytes
              << " bytes) are identical\n";
}

void test_origin_forced_collisions() {
    const MultiBlock m = multi_block(808);
    const SearchResult stage5 = mini_search(m.locus, m.repairs, m.promote, {});
    SearchOptions cached;
    cached.origin_path = true;
    cached.budget = std::size_t{1} << 30;
    panvar::set_origin_cache_hash_collapse(true);
    panvar::set_origin_cache_self_check(true);
    SearchResult got;
    try {
        got = mini_search(m.locus, m.repairs, m.promote, cached);
    } catch (...) {
        panvar::set_origin_cache_hash_collapse(false);
        panvar::set_origin_cache_self_check(false);
        throw;
    }
    panvar::set_origin_cache_hash_collapse(false);
    panvar::set_origin_cache_self_check(false);
    require_same_search(got, stage5, "forced collisions vs stage 5");
    require(got.hits > 0, "collisions must still let genuine hits through");
    std::cout << "  ok: every key in one bucket, " << got.hits << " hits, results unchanged\n";
}

// halo = max supported insert + longest read + max_gap_bp, exactly as the delta computes it.
std::size_t fixture_halo() {
    return panvar::library_insert_support(library()).max_bp + 100 + library().max_gap_bp;
}

void test_origin_translation_across_offsets() {
    // Two copies of one segment at different offsets. The window around the same relative block in
    // each is byte-identical, so the second must be a HIT, served from records stored relative to the
    // first window and translated to the second's offset.
    std::mt19937_64 rng(909);
    const std::size_t halo = fixture_halo();
    const std::string zone = random_bases(rng, halo + 40 + halo);
    const std::string base = random_bases(rng, 700) + zone + random_bases(rng, 700) + zone + random_bases(rng, 700);
    Locus l;
    l.base = base;
    std::size_t n = 0;
    for (std::size_t start = 0; start + 260 <= base.size(); start += 9) {
        l.fragments.push_back(fragment_at(base, start, 260, 100, "z" + std::to_string(n++)));
    }
    const std::string allele = random_bases(rng, 40);
    const std::size_t first_block = 700 + halo;
    const std::size_t second_block = 700 + zone.size() + 700 + halo;
    const std::vector<Repair> repairs = {{first_block, first_block + 40, allele},
                                         {second_block, second_block + 40, allele}};

    const auto run = [&](bool origin) {
        panvar::PreparedScoringContext context(l.fragments, library(), options());
        context.add_haplotypes({l.base});
        context.begin_transient_round();
        context.retain_census(0);
        if (origin) {
            context.enable_origin_cache_path(true);
            context.set_origin_cache_limits(std::size_t{1} << 30, 0);
        }
        std::vector<double> cells;
        std::vector<panvar::BlockDeltaWork> works;
        for (const Repair& r : repairs) {
            panvar::BlockDeltaWork w;
            const std::size_t slot = context.add_repaired_haplotype(0, r.begin, r.end, r.allele, &w);
            for (std::size_t f = 0; f < l.fragments.size(); ++f) cells.push_back(context.mass_of(f, slot));
            works.push_back(w);
        }
        return std::make_pair(cells, works);
    };
    const auto stage5 = run(false);
    panvar::set_origin_cache_self_check(true);
    std::pair<std::vector<double>, std::vector<panvar::BlockDeltaWork>> cached;
    try {
        cached = run(true);
    } catch (...) {
        panvar::set_origin_cache_self_check(false);
        throw;
    }
    panvar::set_origin_cache_self_check(false);
    require(cached.second[0].origin_cache_misses == 1, "the first copy must be enumerated");
    require(cached.second[1].origin_cache_hits == 1,
            "the second copy's identical window must be a hit at a different offset");
    for (std::size_t i = 0; i < stage5.first.size(); ++i) {
        require_identical(cached.first[i], stage5.first[i], "translated cell " + std::to_string(i));
    }
    std::cout << "  ok: identical window " << (second_block - first_block)
              << " bp further along is a hit, translated, and every cell matches stage 5\n";
}

void test_origin_end_flags_never_alias() {
    // The same bytes twice, the second copy flush with the locus END, so its window carries at_end and
    // the first does not. at_end changes which placements the enumeration keeps, so a hit here would
    // be a wrong answer: it must be a miss.
    std::mt19937_64 rng(1010);
    const std::size_t halo = fixture_halo();
    const std::string zone = random_bases(rng, halo + 40 + halo);
    const std::string base = random_bases(rng, 700) + zone + random_bases(rng, 700) + zone;
    Locus l;
    l.base = base;
    std::size_t n = 0;
    for (std::size_t start = 0; start + 260 <= base.size(); start += 9) {
        l.fragments.push_back(fragment_at(base, start, 260, 100, "e" + std::to_string(n++)));
    }
    const std::string allele = random_bases(rng, 40);
    const std::size_t first_block = 700 + halo;
    const std::size_t second_block = 700 + zone.size() + 700 + halo;

    const auto run = [&](bool origin) {
        panvar::PreparedScoringContext context(l.fragments, library(), options());
        context.add_haplotypes({l.base});
        context.begin_transient_round();
        context.retain_census(0);
        if (origin) {
            context.enable_origin_cache_path(true);
            context.set_origin_cache_limits(std::size_t{1} << 30, 0);
        }
        std::vector<double> cells;
        std::vector<panvar::BlockDeltaWork> works;
        for (const std::size_t b : {first_block, second_block}) {
            panvar::BlockDeltaWork w;
            const std::size_t slot = context.add_repaired_haplotype(0, b, b + 40, allele, &w);
            for (std::size_t f = 0; f < l.fragments.size(); ++f) cells.push_back(context.mass_of(f, slot));
            works.push_back(w);
        }
        return std::make_pair(cells, works);
    };
    const auto stage5 = run(false);
    const auto cached = run(true);
    require(cached.second[1].origin_cache_misses == 1,
            "an at_end window must never be served by an interior window's entry");
    for (std::size_t i = 0; i < stage5.first.size(); ++i) {
        require_identical(cached.first[i], stage5.first[i], "end-flag cell " + std::to_string(i));
    }
    std::cout << "  ok: byte-identical windows differing only in at_end do not alias\n";
}

void test_origin_length_shifts_keep_groups_apart() {
    // Three lengths at one near-start block: three shifts, three unaffected sets. Near an end the
    // unaffected predicate reads the shift, so a preaggregation shared across shifts would be wrong
    // for two of them; the self-check compares every merged key against the stage-5 map.
    std::mt19937_64 rng(1111);
    Locus l;
    l.base = random_bases(rng, 2500);
    std::size_t n = 0;
    for (std::size_t start = 0; start + 260 <= l.base.size(); start += 7) {
        l.fragments.push_back(fragment_at(l.base, start, 260, 100, "s" + std::to_string(n++)));
    }
    std::vector<Repair> repairs;
    for (const std::size_t length : {40u, 5u, 71u, 40u, 5u, 71u}) {   // each shift twice: groups reused
        repairs.push_back({30, 70, random_bases(rng, length)});
    }
    const SearchResult stage5 = mini_search(l, repairs, 0, {});
    SearchOptions cached;
    cached.origin_path = true;
    cached.budget = std::size_t{1} << 30;
    panvar::set_origin_cache_self_check(true);
    SearchResult got;
    try {
        got = mini_search(l, repairs, 0, cached);
    } catch (...) {
        panvar::set_origin_cache_self_check(false);
        throw;
    }
    panvar::set_origin_cache_self_check(false);
    require_same_search(got, stage5, "length shifts vs stage 5");
    std::cout << "  ok: shifts -35, 0 and +31 at a near-start block keep separate preaggregations\n";
}

void test_origin_caps_count_executed_work() {
    const MultiBlock m = multi_block(1212);
    SearchOptions none;
    none.origin_path = true;
    SearchOptions full = none;
    full.budget = std::size_t{1} << 30;
    const SearchResult enumerated = mini_search(m.locus, m.repairs, m.promote, none);
    const SearchResult reused = mini_search(m.locus, m.repairs, m.promote, full);
    require(reused.executed_origins < enumerated.executed_origins,
            "hits must not be charged as executed physical origins");

    // Executed-work cap: exactly what the cached search executed is enough for it.
    SearchOptions exact = full;
    exact.max_physical_origins = reused.executed_origins;
    (void)mini_search(m.locus, m.repairs, m.promote, exact);

    // Consumption cap: at the measured figure accepted, one below refused with no ranking.
    SearchOptions at_cap = full;
    at_cap.record_cap = reused.consumed;
    const SearchResult accepted = mini_search(m.locus, m.repairs, m.promote, at_cap);
    require_same_search(accepted, reused, "consumption cap at the bound");
    SearchOptions under = full;
    under.record_cap = reused.consumed - 1;
    bool refused = false;
    try {
        (void)mini_search(m.locus, m.repairs, m.promote, under);
    } catch (const panvar::ResourceExhausted&) {
        refused = true;
    }
    require(refused, "one record below the consumption cap must refuse");
    std::cout << "  ok: executed origins " << reused.executed_origins << " (cached) vs "
              << enumerated.executed_origins << " (enumerated); consumption cap " << reused.consumed
              << " accepts, one below refuses\n";
}

void test_origin_threads_identical() {
    const MultiBlock m = multi_block(1313);
    SearchOptions base;
    base.origin_path = true;
    base.budget = std::size_t{1} << 30;
    const SearchResult one = mini_search(m.locus, m.repairs, m.promote, base);
    for (const std::size_t threads : {2u, 4u, 8u}) {
        SearchOptions how = base;
        how.threads = threads;
        require_same_search(mini_search(m.locus, m.repairs, m.promote, how), one,
                            std::to_string(threads) + " threads");
    }
    std::cout << "  ok: 1/2/4/8 threads bit-identical on the origin path\n";
}

void test_released_slots_are_reused() {
    // The memory cause measured on NA20509: released slots were never reused, so every round appended
    // a neighbourhood of slots, and score_candidates -- one column per slot -- grew by 26 MiB a round.
    // Ten release cycles must now hold the slot count and the scoring result constant, on both paths,
    // while every mass cell stays identical to the first cycle's.
    const Locus l = locus(2718, {400});
    std::mt19937_64 rng(27);
    std::vector<std::string> replacements;
    for (int i = 0; i < 8; ++i) replacements.push_back(random_bases(rng, 40));
    for (const bool origin : {false, true}) {
        panvar::PreparedScoringContext context(l.fragments, library(), options());
        context.add_haplotypes({l.base});
        context.begin_transient_round();
        context.retain_census(0);
        if (origin) {
            context.enable_origin_cache_path(true);
            context.set_origin_cache_limits(std::size_t{1} << 30, 0);
        }
        std::vector<double> first;
        std::size_t slots = 0;
        std::uint64_t result_bytes = 0;
        for (int cycle = 0; cycle < 10; ++cycle) {
            std::vector<double> cells;
            std::size_t probe_slot = 0;
            for (const std::string& replacement : replacements) {
                probe_slot = context.add_repaired_haplotype(0, 400, 440, replacement);
                for (std::size_t f = 0; f < l.fragments.size(); ++f) cells.push_back(context.mass_of(f, probe_slot));
            }
            (void)context.score_candidates({{0, probe_slot}});
            const std::uint64_t result = context.memory_account().last_scoring_result_bytes;
            if (cycle == 0) {
                first = cells;
                slots = context.unique_haplotypes();
                result_bytes = result;
            } else {
                for (std::size_t i = 0; i < first.size(); ++i) {
                    require_identical(cells[i], first[i],
                                      "cycle " + std::to_string(cycle) + " cell " + std::to_string(i));
                }
                require(context.unique_haplotypes() == slots,
                        "released slots must be reused, not appended: " +
                            std::to_string(context.unique_haplotypes()) + " after " + std::to_string(slots));
                require(result == result_bytes, "the scoring result must not grow across cycles");
            }
            context.release_transient({});
        }
        require(context.haplotypes_created() == 1 + replacements.size() * 10,
                "haplotypes_created must stay cumulative");
        std::cout << "  ok: " << (origin ? "origin path" : "stage-5 path") << ": 10 cycles, " << slots
                  << " slots held throughout, scoring result constant at " << result_bytes
                  << " bytes, masses identical\n";
    }
}

void test_released_seed_is_respelled_identically() {
    // The census-lifetime correction: repair search now opens the transient round BEFORE the seed, so a
    // seed that stops being an incumbent is released with its full census. Two things must then hold.
    // Only the promoted incumbent keeps a census. And a later neighbour that spells the seed's haplotype
    // again -- which the search can generate, since it may restore one block on one homologue -- is
    // rebuilt by delta from the new incumbent and must come back bit-identical to the seed's own full
    // placement, cell for cell, in exposure, and in the diplotype score. If it did not, releasing the
    // seed would change a later round's ranking.
    const MultiBlock m = multi_block(4242);
    const Repair& chosen = m.repairs[m.promote];   // same length, so the winner's coordinates hold
    const std::string original = m.locus.base.substr(chosen.begin, chosen.end - chosen.begin);
    for (const bool origin : {false, true}) {
        const std::string label = origin ? "origin path" : "stage-5 path";
        panvar::PreparedScoringContext context(m.locus.fragments, library(), options());
        if (origin) {
            context.enable_origin_cache_path(true);
            context.set_origin_cache_limits(std::size_t{1} << 30, 0);
        }
        context.begin_transient_round();
        const std::size_t seed = context.add_haplotypes({m.locus.base}).front();
        context.retain_census(seed);
        std::vector<double> seed_cells;
        for (std::size_t f = 0; f < m.locus.fragments.size(); ++f) seed_cells.push_back(context.mass_of(f, seed));
        const double seed_exposure = context.exposure_of(seed);
        panvar::BlockDeltaWork promoted;
        const std::size_t winner =
            context.add_repaired_haplotype(seed, chosen.begin, chosen.end, chosen.allele, &promoted, true);
        const double seed_score = context.score_candidates({{seed, winner}}).scores[0].log_likelihood;
        require(context.census_slots().size() == 2, label + ": seed and winner both hold a census");
        context.release_transient({winner});
        require(context.census_slots() == std::vector<std::size_t>{winner},
                label + ": after promotion only the incumbent may keep a full census");
        panvar::BlockDeltaWork work;
        const std::size_t again = context.add_repaired_haplotype(
            winner, chosen.begin, chosen.begin + chosen.allele.size(), original, &work);
        require(!work.fell_back_to_full, label + ": the respelled seed must be built by delta");
        for (std::size_t f = 0; f < seed_cells.size(); ++f) {
            require_identical(context.mass_of(f, again), seed_cells[f],
                              label + ": respelled seed, mass cell f" + std::to_string(f));
        }
        require_identical(context.exposure_of(again), seed_exposure, label + ": respelled seed exposure");
        require_identical(context.score_candidates({{again, winner}}).scores[0].log_likelihood, seed_score,
                          label + ": respelled seed, diplotype score");
        std::cout << "  ok: " << label << ": seed released at promotion (census slots {" << winner
                  << "}), respelled by delta bit-identically over " << seed_cells.size() << " cells\n";
    }
}

// ---- planned selective retention -----------------------------------------------------------------
struct PlannedResult {
    std::vector<double> cells;                              // canonical: by repair, then by round
    std::vector<std::vector<std::uint64_t>> cache_after;    // live key hashes after each round
    std::vector<panvar::OriginCachePlanReport> plans;
    std::uint64_t hits = 0, misses = 0, not_selected = 0, max_bytes = 0, last_round_misses = 0;
};

// Round A scores every repair against the seed (no size known yet, so nothing can be selected); the
// promoted winner is kept; rounds B-D score the rest against it. B learns nothing new about sizes but
// admits what its plan selected; C and D should then hit. `order_seed` permutes the request order.
PlannedResult planned_search(const MultiBlock& m, std::size_t budget, bool origin, std::size_t threads,
                             unsigned order_seed) {
    panvar::FragmentScoringOptions opts = options();
    opts.threads = threads;
    panvar::PreparedScoringContext context(m.locus.fragments, library(), opts);
    context.add_haplotypes({m.locus.base});
    context.begin_transient_round();
    context.retain_census(0);
    if (origin) {
        context.enable_origin_cache_path(true);
        context.set_origin_cache_limits(budget, 0);
    }
    std::vector<std::size_t> order(m.repairs.size());
    std::iota(order.begin(), order.end(), std::size_t{0});
    if (order_seed != 0) {
        std::mt19937_64 rng(order_seed);
        std::shuffle(order.begin(), order.end(), rng);
    }
    PlannedResult out;
    const auto round = [&](std::size_t incumbent, bool skip_promoted) {
        if (origin) {
            std::vector<panvar::PlannedRepairWindow> planned;
            for (const std::size_t i : order) {
                if (skip_promoted && i == m.promote) continue;
                const Repair& r = m.repairs[i];
                planned.push_back(context.plan_repaired_window(incumbent, r.begin, r.end, r.allele));
            }
            out.plans.push_back(context.plan_origin_cache(planned));
            require(out.plans.back().selected_bytes <= budget, "a plan selected more than its budget");
        }
        std::vector<std::vector<double>> by_repair(m.repairs.size());
        std::uint64_t misses = 0;
        for (const std::size_t i : order) {
            if (skip_promoted && i == m.promote) continue;
            const Repair& r = m.repairs[i];
            panvar::BlockDeltaWork w;
            const std::size_t slot = context.add_repaired_haplotype(incumbent, r.begin, r.end, r.allele, &w);
            for (std::size_t f = 0; f < m.locus.fragments.size(); ++f) by_repair[i].push_back(context.mass_of(f, slot));
            by_repair[i].push_back(context.exposure_of(slot));
            out.hits += w.origin_cache_hits;
            misses += w.origin_cache_misses;
            out.not_selected += w.origin_cache_not_selected;
            if (origin) {
                out.max_bytes = std::max<std::uint64_t>(out.max_bytes, context.origin_cache_bytes());
                require(context.origin_cache_bytes() <= budget, "the cache exceeded its declared budget");
            }
        }
        out.misses += misses;
        out.last_round_misses = misses;
        for (const auto& cells : by_repair) out.cells.insert(out.cells.end(), cells.begin(), cells.end());
        if (origin) out.cache_after.push_back(context.origin_cache_key_hashes());
    };
    round(0, false);
    const Repair& chosen = m.repairs[m.promote];
    const std::size_t winner = context.add_repaired_haplotype(0, chosen.begin, chosen.end, chosen.allele, nullptr, true);
    context.release_transient({winner});
    for (int r = 0; r < 3; ++r) {
        round(winner, true);
        context.release_transient({winner});
    }
    return out;
}

void require_same_cells(const PlannedResult& got, const PlannedResult& want, const std::string& label) {
    require(got.cells.size() == want.cells.size(), label + ": cell counts differ");
    for (std::size_t i = 0; i < want.cells.size(); ++i) {
        require_identical(got.cells[i], want.cells[i], label + ": cell " + std::to_string(i));
    }
}

void test_planned_cache_scores_identically() {
    const MultiBlock m = multi_block(2027);
    const PlannedResult stage5 = planned_search(m, 0, false, 1, 0);
    const PlannedResult full = planned_search(m, std::size_t{1} << 30, true, 1, 0);
    require_same_cells(full, stage5, "planned, full budget vs stage 5");
    require(full.plans[0].known_windows == 0 && full.plans[0].selected_entries == 0,
            "round A knows no exact size, so its plan must select nothing");
    require(full.plans[1].selected_entries > 0, "round B must select the windows round A measured");
    require(full.last_round_misses == 0,
            "with room for everything, the last round must be all hits: a planned identity that differed "
            "from the requested one would miss here");
    const std::size_t small = static_cast<std::size_t>(full.plans[1].selected_bytes / 2);
    const PlannedResult some = planned_search(m, small, true, 1, 0);
    require_same_cells(some, stage5, "planned, half budget vs stage 5");
    require(some.max_bytes <= small, "the undersized cache exceeded its budget");
    require(some.hits < full.hits && some.hits > 0, "half the budget must hit less, but still hit");
    const PlannedResult none = planned_search(m, 0, true, 1, 0);
    require_same_cells(none, stage5, "planned, disabled vs stage 5");
    require(none.hits == 0, "a disabled cache must not hit");
    std::cout << "  ok: disabled, half (" << small << " bytes, " << some.hits << " hits) and full ("
              << full.hits << " hits, last round " << full.last_round_misses
              << " misses) are bit-identical to stage 5\n";
}

void test_planned_cache_order_independent() {
    const MultiBlock m = multi_block(2028);
    const PlannedResult base = planned_search(m, std::size_t{1} << 30, true, 1, 0);
    for (const std::size_t budget : {std::size_t{1} << 30, static_cast<std::size_t>(base.plans[1].selected_bytes / 3)}) {
        const PlannedResult reference = planned_search(m, budget, true, 1, 0);
        for (const unsigned seed : {7u, 11u, 29u}) {
            const PlannedResult shuffled = planned_search(m, budget, true, 1, seed);
            const std::string label = "budget " + std::to_string(budget) + ", order " + std::to_string(seed);
            require_same_cells(shuffled, reference, label);
            require(shuffled.cache_after == reference.cache_after, label + ": cache contents depend on request order");
            require(shuffled.hits == reference.hits && shuffled.misses == reference.misses &&
                        shuffled.not_selected == reference.not_selected,
                    label + ": hit/miss accounting depends on request order");
            for (std::size_t r = 0; r < reference.plans.size(); ++r) {
                const auto& a = shuffled.plans[r];
                const auto& b = reference.plans[r];
                require(a.selected_entries == b.selected_entries && a.selected_bytes == b.selected_bytes &&
                            a.evicted_entries == b.evicted_entries && a.known_windows == b.known_windows,
                        label + ": the plan depends on request order");
            }
        }
    }
    for (const std::size_t threads : {2u, 4u, 8u}) {
        const PlannedResult t = planned_search(m, static_cast<std::size_t>(base.plans[1].selected_bytes / 3), true, threads, 0);
        const PlannedResult one = planned_search(m, static_cast<std::size_t>(base.plans[1].selected_bytes / 3), true, 1, 0);
        require_same_cells(t, one, std::to_string(threads) + " threads");
        require(t.cache_after == one.cache_after, std::to_string(threads) + " threads: cache contents differ");
    }
    std::cout << "  ok: three request orders x two budgets give identical cells, plans, cache contents and "
                 "accounting; 1/2/4/8 threads identical\n";
}

void test_planned_eviction_releases_storage() {
    // Alternate the planned set between two disjoint halves of the repairs. Every switch must evict the
    // other half, release its storage and reuse its slots: slots and retained capacity may not grow.
    const MultiBlock m = multi_block(2029);
    panvar::PreparedScoringContext context(m.locus.fragments, library(), options());
    context.add_haplotypes({m.locus.base});
    context.begin_transient_round();
    context.retain_census(0);
    context.enable_origin_cache_path(true);
    context.set_origin_cache_limits(std::size_t{1} << 30, 0);
    const auto request = [&](int parity, bool plan) {
        if (plan) {
            std::vector<panvar::PlannedRepairWindow> planned;
            for (std::size_t i = 0; i < m.repairs.size(); ++i) {
                if (parity >= 0 && static_cast<int>(i % 2) != parity) continue;
                const Repair& r = m.repairs[i];
                planned.push_back(context.plan_repaired_window(0, r.begin, r.end, r.allele));
            }
            (void)context.plan_origin_cache(planned);
        }
        for (std::size_t i = 0; i < m.repairs.size(); ++i) {
            if (parity >= 0 && static_cast<int>(i % 2) != parity) continue;
            const Repair& r = m.repairs[i];
            (void)context.add_repaired_haplotype(0, r.begin, r.end, r.allele);
        }
        context.release_transient({});
    };
    request(-1, false);   // every size becomes known
    std::uint64_t slots_after_warm = 0;
    std::uint64_t capacity[2] = {0, 0};
    for (int cycle = 0; cycle < 10; ++cycle) {
        request(cycle % 2, true);
        const std::uint64_t cap = context.memory_account().origin_cache_capacity;
        if (cycle == 1) slots_after_warm = context.origin_cache_slots();
        if (cycle >= 2) {
            require(context.origin_cache_slots() == slots_after_warm, "evicted slots must be reused, not appended");
            require(cap == capacity[cycle % 2], "retained capacity must return to the same level each cycle");
        }
        capacity[cycle % 2] = cap;
    }
    std::cout << "  ok: 10 alternating plans, " << slots_after_warm << " slots throughout, capacity "
              << capacity[0] << " / " << capacity[1] << " bytes per half, no growth\n";
}

// ---- planned retention within a round: the canonical batch ----------------------------------------
struct BatchResult {
    std::vector<double> cells;                               // canonical: by request index, then round
    std::vector<double> scores;                              // every (A-request, B-request) diplotype
    std::vector<std::vector<std::uint64_t>> cache_after;     // live key hashes after each round
    std::vector<std::uint64_t> accounting;                   // hits, misses, reuses, stored, evicted, ...
    std::uint64_t max_bytes = 0, slots_after_second = 0, slots_after_last = 0, group_reuses = 0, evicted = 0;
};

// Two incumbents that differ only at one point far from every repair, so the same allele repaired on
// either spells the same window: multiplicity 2, and a group reuse for the second. Three rounds of the
// same batch; `order_seed` permutes the requests the caller offers.
BatchResult batch_search(const MultiBlock& m, std::size_t budget, bool origin, std::size_t threads,
                         unsigned order_seed) {
    panvar::FragmentScoringOptions opts = options();
    opts.threads = threads;
    panvar::PreparedScoringContext context(m.locus.fragments, library(), opts);
    std::string other = m.locus.base;
    const std::size_t point = 2250;
    other[point] = other[point] == 'A' ? 'C' : 'A';
    context.begin_transient_round();
    const std::vector<std::size_t> seeds = context.add_haplotypes({m.locus.base, other});
    context.retain_census(seeds[0]);
    context.retain_census(seeds[1]);
    if (origin) {
        context.enable_origin_cache_path(true);
        context.set_origin_cache_limits(budget, 0);
    }
    std::vector<panvar::RepairRequest> requests;
    for (const std::size_t seed_slot : seeds) {
        for (const Repair& r : m.repairs) requests.push_back({seed_slot, r.begin, r.end, r.allele});
    }
    std::vector<std::size_t> order(requests.size());
    std::iota(order.begin(), order.end(), std::size_t{0});
    if (order_seed != 0) {
        std::mt19937_64 rng(order_seed);
        std::shuffle(order.begin(), order.end(), rng);
    }
    BatchResult out;
    std::uint64_t hits = 0, misses = 0, stored = 0, evicted = 0, not_selected = 0, declined = 0;
    for (int round = 0; round < 3; ++round) {
        std::vector<panvar::RepairRequest> offered;
        for (const std::size_t i : order) offered.push_back(requests[i]);
        std::vector<std::size_t> slots(requests.size());
        if (origin) {
            std::vector<panvar::BlockDeltaWork> works;
            const std::vector<std::size_t> got = context.add_repaired_haplotypes(offered, &works);
            for (std::size_t k = 0; k < order.size(); ++k) {
                slots[order[k]] = got[k];
                const panvar::BlockDeltaWork& w = works[k];
                hits += w.origin_cache_hits;
                misses += w.origin_cache_misses;
                out.group_reuses += w.origin_group_reuses;
                stored += w.origin_cache_stored;
                evicted += w.origin_cache_evicted;
                not_selected += w.origin_cache_not_selected;
                declined += w.origin_cache_declined;
            }
            out.max_bytes = std::max<std::uint64_t>(out.max_bytes, context.origin_cache_bytes());
            require(context.origin_cache_bytes() <= budget, "the cache exceeded its declared budget");
        } else {
            for (std::size_t i = 0; i < requests.size(); ++i) {
                const panvar::RepairRequest& r = requests[i];
                slots[i] = context.add_repaired_haplotype(r.incumbent, r.block_begin, r.block_end, r.replacement);
            }
        }
        for (std::size_t i = 0; i < requests.size(); ++i) {
            for (std::size_t f = 0; f < m.locus.fragments.size(); ++f) out.cells.push_back(context.mass_of(f, slots[i]));
            out.cells.push_back(context.exposure_of(slots[i]));
        }
        std::vector<panvar::HaplotypeCandidate> pairs;
        const std::size_t half = m.repairs.size();
        for (std::size_t i = 0; i < half; ++i) pairs.push_back({slots[i], slots[half + i]});
        const panvar::FragmentScoringResult scored = context.score_candidates(pairs);
        for (const auto& score : scored.scores) out.scores.push_back(score.log_likelihood);
        if (origin) out.cache_after.push_back(context.origin_cache_key_hashes());
        if (round == 1) out.slots_after_second = context.origin_cache_slots();
        context.release_transient(seeds);
    }
    out.slots_after_last = context.origin_cache_slots();
    out.evicted = evicted;
    out.accounting = {hits, misses, out.group_reuses, stored, evicted, not_selected, declined,
                      context.origin_cache_bytes(), context.origin_records_consumed(),
                      context.cumulative_counters().physical_origins};
    return out;
}

void test_batch_retention_is_order_independent() {
    const MultiBlock m = multi_block(3031);
    const BatchResult stage5 = batch_search(m, 0, false, 1, 0);
    const BatchResult full = batch_search(m, std::size_t{1} << 30, true, 1, 0);
    // With nothing retained, a window's second requester can only be served by the group holder.
    const BatchResult unretained = batch_search(m, 0, true, 1, 0);
    require(unretained.group_reuses > 0, "the fixture must give some window two requesters in one round");
    const std::size_t tight = static_cast<std::size_t>(full.accounting[7] / 3);
    std::uint64_t evictions_seen = 0;
    for (const std::size_t budget : {std::size_t{1} << 30, tight, std::size_t{0}}) {
        const BatchResult reference = batch_search(m, budget, true, 1, 0);
        const std::string at = "budget " + std::to_string(budget);
        for (std::size_t i = 0; i < stage5.cells.size(); ++i) {
            require_identical(reference.cells[i], stage5.cells[i], at + ": cell " + std::to_string(i) + " vs stage 5");
        }
        for (std::size_t i = 0; i < stage5.scores.size(); ++i) {
            require_identical(reference.scores[i], stage5.scores[i], at + ": score " + std::to_string(i) + " vs stage 5");
        }
        require(reference.max_bytes <= budget, at + ": payload exceeded the budget");
        require(reference.slots_after_last == reference.slots_after_second,
                at + ": cache slots grew after the selection settled");
        evictions_seen += reference.evicted;
        for (const unsigned seed : {5u, 17u, 1234567u}) {
            const BatchResult shuffled = batch_search(m, budget, true, 1, seed);
            const std::string label = at + ", order " + std::to_string(seed);
            for (std::size_t i = 0; i < reference.cells.size(); ++i) {
                require_identical(shuffled.cells[i], reference.cells[i], label + ": cell " + std::to_string(i));
            }
            for (std::size_t i = 0; i < reference.scores.size(); ++i) {
                require_identical(shuffled.scores[i], reference.scores[i], label + ": score " + std::to_string(i));
            }
            require(shuffled.cache_after == reference.cache_after, label + ": cache identities depend on order");
            require(shuffled.accounting == reference.accounting, label + ": cache accounting depends on order");
        }
        for (const std::size_t threads : {2u, 4u, 8u}) {
            const BatchResult t = batch_search(m, budget, true, threads, 0);
            require(t.cells == reference.cells && t.scores == reference.scores, at + ": threads changed a value");
            require(t.cache_after == reference.cache_after && t.accounting == reference.accounting,
                    at + ": threads changed the cache");
        }
    }
    require(evictions_seen > 0, "the tight budget must force in-round replacement, or replacement is untested");
    std::cout << "  ok: 4 request orders x 3 budgets x 1/2/4/8 threads: identical cells, scores, cache identities "
                 "and accounting; " << unretained.group_reuses << " group reuses unretained, " << evictions_seen
              << " in-round evictions at the tight budget; matches stage 5\n";
}

// ---- the necessary-condition piece screen (measurement only) ------------------------------------
// Builds a read from haplotype[start ...] under an explicit edit script, then keeps its first `length`
// bases. Kinds: 's' substitution, 'i' insertion (a read base the haplotype lacks), 'd' deletion (a
// haplotype base the read skips). Positions are read coordinates.
std::string edited_read(const std::string& hap, std::size_t start, std::size_t length,
                        std::vector<std::pair<char, std::size_t>> script, std::mt19937_64& rng) {
    std::sort(script.begin(), script.end(), [](const auto& a, const auto& b) { return a.second < b.second; });
    std::string read;
    std::size_t h = start;
    std::size_t e = 0;
    const std::string bases = "ACGT";
    while (read.size() < length && h < hap.size()) {
        if (e < script.size() && script[e].second <= read.size()) {
            const char kind = script[e++].first;
            if (kind == 'i') {
                read.push_back(bases[rng() % 4]);
            } else if (kind == 'd') {
                ++h;
            } else {
                char b = bases[rng() % 4];
                while (b == hap[h]) b = bases[rng() % 4];
                read.push_back(b);
                ++h;
            }
            continue;
        }
        read.push_back(hap[h++]);
    }
    return read;
}

void test_piece_screen_is_a_necessary_condition() {
    const std::size_t band = options().mismatch_band;
    const std::size_t gap = library().max_gap_bp;
    const std::vector<std::size_t> pieces = {band + 1, band + 2, 2 * band + 1};
    const std::size_t length = 100;
    std::mt19937_64 rng(4141);
    std::uint64_t checked = 0, feasible = 0, screened = 0, consistency_only = 0;
    const auto check = [&](const std::string& hap, std::size_t start, const std::string& read,
                           const std::string& label) {
        if (read.size() != length) return;
        const bool prefilter = panvar::feasibility_prefilter_admits(hap, start, read, band, gap);
        if (prefilter) ++feasible;
        for (const std::size_t d : pieces) {
            const bool admit = panvar::piece_screen_admits(hap, start, read, band, gap, d);
            const bool reference = panvar::piece_screen_admits_reference(hap, start, read, band, gap, d);
            const std::string at = label + " (d=" + std::to_string(d) + ", start " + std::to_string(start) + ")";
            require(admit == reference, at + ": the screen disagrees with its brute-force reference");
            if (prefilter) require(admit, at + ": FALSE NEGATIVE -- the prefilter admits a start the screen rejects");
            if (!admit) {
                ++screened;
                if (panvar::piece_screen_count_passes(hap, start, read, band, gap, d)) ++consistency_only;
            }
            ++checked;
        }
    };
    const std::string unique = random_bases(rng, 1200);
    std::string repeat;
    {
        const std::string unit = random_bases(rng, 7);
        while (repeat.size() < 1200) repeat += unit;
        for (int v = 0; v < 12; ++v) repeat[rng() % repeat.size()] = "ACGT"[rng() % 4];
    }
    std::string period14;
    {
        const std::string unit = random_bases(rng, 14);
        while (period14.size() < 1200) period14 += unit;
    }
    // Directed fixtures: each is feasible, and each breaks a naive offset rule.
    const std::size_t s0 = 300;
    const std::size_t chunk = length / (band + 1);
    struct Directed {
        std::string name;
        std::vector<std::pair<char, std::size_t>> script;
        std::size_t start_shift;   // read taken from hap[s0 + shift] but tested at s0 (leading deletions)
    };
    const std::vector<Directed> directed = {
        {"insertion inside piece 1: displacement falls across a broken piece", {{'i', chunk + 3}}, 0},
        {"deletion exactly at a piece boundary: displacement rises between adjacent exact pieces", {{'d', 2 * chunk}}, 0},
        {"two leading deletions: the first exact piece sits at displacement +2", {}, 2},
        {"insertion at the first read base", {{'i', 0}}, 0},
        {"deletion just before the last read base", {{'d', length - 1}}, 0},
        {"k substitutions, one per piece boundary", {}, 0},
        {"mixed k edits: insertion, deletion, substitutions", {{'i', 10}, {'d', 40}, {'s', 55}, {'i', 70}, {'d', 85}, {'s', 95}}, 0},
        {"k insertions clustered, then k deletions would exceed", {{'i', 20}, {'i', 21}, {'i', 22}}, 0},
    };
    for (const auto& hap : {unique, repeat, period14}) {
        for (Directed d : directed) {
            if (d.name.rfind("k substitutions", 0) == 0) {
                d.script.clear();
                for (std::size_t j = 0; j < band; ++j) d.script.push_back({'s', (j + 1) * chunk});
            }
            const std::string read = edited_read(hap, s0 + d.start_shift, length, d.script, rng);
            require(read.size() == length, d.name + ": fixture read has the wrong length");
            // Leading deletions are the alignment consuming haplotype bases before the first read base.
            std::string at_s0 = read;
            check(hap, s0, at_s0, d.name);
            require(panvar::feasibility_prefilter_admits(hap, s0, at_s0, band, gap),
                    d.name + ": the directed fixture must be feasible, or it tests nothing");
        }
    }
    // Random edit scripts: kinds mixed or pure, counts 0..k+3, positions random or on piece boundaries,
    // tested at the true start and at neighbouring starts.
    for (int trial = 0; trial < 3000; ++trial) {
        const std::string& hap = trial % 3 == 0 ? unique : (trial % 3 == 1 ? repeat : period14);
        const std::size_t start = 200 + rng() % 700;
        const std::size_t edits = rng() % (band + 4);
        const int mode = static_cast<int>(rng() % 4);   // 0 mixed, 1 subs, 2 ins, 3 dels
        const std::size_t d = pieces[rng() % pieces.size()];
        std::vector<std::pair<char, std::size_t>> script;
        for (std::size_t e = 0; e < edits; ++e) {
            const char kind = mode == 1 ? 's' : mode == 2 ? 'i' : mode == 3 ? 'd' : "sid"[rng() % 3];
            std::size_t at = rng() % length;
            if (rng() % 2 == 0) at = std::min<std::size_t>(length - 1, (rng() % d) * length / d);   // a boundary
            script.push_back({kind, at});
        }
        const std::string read = edited_read(hap, start, length, script, rng);
        for (int shift = -3; shift <= 3; ++shift) {
            check(hap, static_cast<std::size_t>(static_cast<long>(start) + shift), read, "random trial " + std::to_string(trial));
        }
    }
    // Scrambled reads: each piece copied from its own displacement in [-k, k]. The count condition
    // passes for every one, so any rejection here is the offset consistency at work.
    for (int trial = 0; trial < 1500; ++trial) {
        const std::string& hap = trial % 2 == 0 ? unique : period14;
        const std::size_t start = 200 + rng() % 700;
        const std::size_t d = pieces[rng() % pieces.size()];
        std::string read;
        for (std::size_t j = 0; j < d; ++j) {
            const std::size_t o = j * length / d;
            const std::size_t l = (j + 1) * length / d - o;
            const long delta = static_cast<long>(rng() % (2 * band + 1)) - static_cast<long>(band);
            const long at = std::max<long>(static_cast<long>(start + o) + delta, static_cast<long>(start));
            read += hap.substr(static_cast<std::size_t>(at), l);
        }
        check(hap, start, read, "scrambled trial " + std::to_string(trial));
    }
    require(feasible > 1000, "the fixtures must include many feasible placements");
    require(screened > 1000, "the fixtures must include placements the screen rejects");
    require(consistency_only > 0, "some rejection must come from offset consistency alone, or it is untested");
    std::cout << "  ok: " << checked << " screen evaluations (" << feasible
              << " prefilter-feasible placements, every one admitted); " << screened << " rejections, "
              << consistency_only << " by offset consistency alone; screen == brute-force reference throughout\n";
}

// ---- window enumeration across workers ------------------------------------------------------------
std::vector<std::uint64_t> work_fields(const panvar::BlockDeltaWork& w) {
    // Every work and accounting field. Excluded: timings, and the three allocation-growth diagnostics,
    // which count each worker's own scratch growth and so necessarily differ with the worker count.
    return {w.origins_reused, w.origins_recomputed, w.local_window_bases, w.piece_hits, w.dp_executions,
            w.cache_misses, w.cache_requests, w.cache_hits, w.census_payload_bytes_predicted,
            w.census_payload_bytes_observed, w.local_window_payload_bytes, w.full_haplotype_bytes_materialised,
            w.local_window_bytes_materialised, w.sequence_key_bytes, w.origin_cache_hits, w.origin_cache_misses,
            w.origin_cache_stored, w.origin_cache_declined, w.origin_cache_not_selected, w.origin_group_reuses,
            w.origin_cache_evicted, w.origin_cache_evicted_bytes, w.origin_cache_admitted_bytes,
            w.origin_cache_declined_bytes, w.cached_origin_records_consumed, w.physical_origins_constructed,
            w.fragments_reduced, w.fragments_copied, static_cast<std::uint64_t>(w.fell_back_to_full)};
}

struct ParallelTrace {
    std::vector<double> values;           // masses, exposures and diplotype scores
    std::vector<std::uint64_t> work;      // every field of every call, then the cumulative counters
    std::vector<std::uint64_t> digests;   // the enumerated origins and offsets, bit for bit, per call
};

// Substitutions, insertions, deletions and empty alleles at interior, start and end windows; a round
// against the seed, a promotion, a round against the winner, then one planned batch. `mode`: 0 the
// stage-5 path, 1 the origin cache with room for everything, 2 a tight planned budget that evicts.
ParallelTrace parallel_trace(const Locus& l, const std::vector<Repair>& repairs, std::size_t promote,
                             std::size_t threads, int mode, std::size_t tight_budget,
                             std::uint64_t max_origins = 0, std::uint64_t max_piece_hits = 0) {
    panvar::FragmentScoringOptions opts = options();
    opts.threads = threads;
    if (max_origins != 0) opts.max_physical_origins = max_origins;
    if (max_piece_hits != 0) opts.max_piece_hits = max_piece_hits;
    panvar::PreparedScoringContext context(l.fragments, library(), opts);
    context.begin_transient_round();
    context.add_haplotypes({l.base});
    context.retain_census(0);
    if (mode != 0) {
        context.enable_origin_cache_path(true);
        context.set_origin_cache_limits(mode == 1 ? std::size_t{1} << 30 : tight_budget, 0);
    }
    ParallelTrace out;
    const auto record = [&](std::size_t slot, const panvar::BlockDeltaWork& w) {
        for (std::size_t f = 0; f < l.fragments.size(); ++f) out.values.push_back(context.mass_of(f, slot));
        out.values.push_back(context.exposure_of(slot));
        const std::vector<std::uint64_t> fields = work_fields(w);
        out.work.insert(out.work.end(), fields.begin(), fields.end());
        out.digests.push_back(context.last_window_digest());
    };
    std::vector<std::size_t> round_a;
    for (const Repair& r : repairs) {
        panvar::BlockDeltaWork w;
        const std::size_t slot = context.add_repaired_haplotype(0, r.begin, r.end, r.allele, &w);
        record(slot, w);
        round_a.push_back(slot);
    }
    for (std::size_t i = 0; i + 1 < round_a.size(); i += 2) {
        out.values.push_back(context.score_candidates({{round_a[i], round_a[i + 1]}}).scores[0].log_likelihood);
    }
    panvar::BlockDeltaWork promoted;
    const Repair& chosen = repairs[promote];
    const std::size_t winner = context.add_repaired_haplotype(0, chosen.begin, chosen.end, chosen.allele, &promoted, true);
    record(winner, promoted);
    context.release_transient({winner});
    for (std::size_t i = 0; i < repairs.size(); ++i) {
        if (i == promote) continue;
        panvar::BlockDeltaWork w;
        const std::size_t slot = context.add_repaired_haplotype(winner, repairs[i].begin, repairs[i].end, repairs[i].allele, &w);
        record(slot, w);
    }
    context.release_transient({winner});
    if (mode != 0) {
        std::vector<panvar::RepairRequest> batch;
        for (std::size_t i = 0; i < repairs.size(); ++i) {
            if (i != promote) batch.push_back({winner, repairs[i].begin, repairs[i].end, repairs[i].allele});
        }
        std::vector<panvar::BlockDeltaWork> works;
        const std::vector<std::size_t> slots = context.add_repaired_haplotypes(batch, &works);
        for (std::size_t k = 0; k < slots.size(); ++k) record(slots[k], works[k]);
        context.release_transient({winner});
    }
    const panvar::FragmentScoringCounters& c = context.cumulative_counters();
    out.work.insert(out.work.end(), {c.candidate_bases_scanned, c.piece_hits, c.proposed_starts, c.unique_starts,
                                     c.full_read_verifications, c.dp_evaluations, c.prefilter_rejections,
                                     c.accepted_mate_placements, c.valid_fr_states, c.latent_alignment_pairs,
                                     c.physical_origins, c.statistical_state_classes,
                                     context.origin_records_consumed(), context.origin_cache_bytes(),
                                     static_cast<std::uint64_t>(context.origin_cache_entries())});
    return out;
}

void require_same_trace(const ParallelTrace& got, const ParallelTrace& want, const std::string& label) {
    require(got.values.size() == want.values.size(), label + ": value counts differ");
    for (std::size_t i = 0; i < want.values.size(); ++i) {
        require_identical(got.values[i], want.values[i], label + ": value " + std::to_string(i));
    }
    require(got.work == want.work, label + ": a work or accounting counter differs");
    require(got.digests == want.digests, label + ": enumerated origins or their order differ");
}

Locus tandem_locus(std::uint64_t seed) {
    std::mt19937_64 rng(seed);
    const std::string unit = random_bases(rng, 10);
    Locus l;
    l.base = random_bases(rng, 900);
    for (int i = 0; i < 60; ++i) l.base += unit;
    l.base += random_bases(rng, 900);
    std::size_t n = 0;
    for (std::size_t start = 0; start + 260 <= l.base.size(); start += 9) {
        l.fragments.push_back(fragment_at(l.base, start, 260, 100, "t" + std::to_string(n++)));
    }
    return l;
}

void test_parallel_enumeration_is_exact() {
    MultiBlock m = multi_block(5151);
    std::mt19937_64 rng(5152);
    for (const std::size_t b : {50u, 3000u, 5900u}) {   // empty alleles at start, interior and end blocks
        m.repairs.push_back({b, std::min<std::size_t>(b + 40, m.locus.base.size()), std::string()});
    }
    const Locus t = tandem_locus(5153);
    const std::vector<Repair> tandem_repairs = {
        {1000, 1040, random_bases(rng, 40)}, {1100, 1140, random_bases(rng, 55)}, {1200, 1240, random_bases(rng, 12)},
        {1300, 1340, std::string()}, {1000, 1040, t.base.substr(1010, 40)}};   // the last is a shifted copy
    std::size_t traces = 0;
    for (int mode = 0; mode <= 2; ++mode) {
        using Fixture = std::pair<const Locus*, const std::vector<Repair>*>;
        for (const Fixture& fixture : {Fixture{&m.locus, &m.repairs}, Fixture{&t, &tandem_repairs}}) {
            const std::size_t promote = fixture.second == &m.repairs ? m.promote : 0;
            const ParallelTrace one = parallel_trace(*fixture.first, *fixture.second, promote, 1, mode, 60000);
            for (const std::size_t threads : {2u, 4u, 8u}) {
                const ParallelTrace many = parallel_trace(*fixture.first, *fixture.second, promote, threads, mode, 60000);
                require_same_trace(many, one, "mode " + std::to_string(mode) + ", " +
                                                  (fixture.second == &m.repairs ? "multi-block" : "tandem") + ", " +
                                                  std::to_string(threads) + " threads");
                ++traces;
            }
        }
    }
    std::cout << "  ok: " << traces << " traces (stage 5, full cache, evicting cache x multi-block, tandem x "
                 "2/4/8 threads) bit-identical to one thread: masses, exposures, scores, origins and order, "
                 "every work and cache counter\n";
}

void test_parallel_caps_exact() {
    const MultiBlock m = multi_block(6161);
    const ParallelTrace reference = parallel_trace(m.locus, m.repairs, m.promote, 1, 0, 0);
    // The cumulative physical origins and piece hits are the 11th and 2nd of the trailing counters.
    const std::size_t tail = reference.work.size() - 15;
    const std::uint64_t origins = reference.work[tail + 10];
    const std::uint64_t piece_hits = reference.work[tail + 1];
    for (const std::size_t threads : {1u, 2u, 4u, 8u}) {
        const std::string at = std::to_string(threads) + " threads";
        require_same_trace(parallel_trace(m.locus, m.repairs, m.promote, threads, 0, 0, origins, 0), reference,
                           at + ", physical-origin cap exactly at the total");
        require_same_trace(parallel_trace(m.locus, m.repairs, m.promote, threads, 0, 0, 0, piece_hits), reference,
                           at + ", piece-hit cap exactly at the total");
        for (const int which : {0, 1}) {
            bool refused = false;
            std::string message;
            try {
                (void)parallel_trace(m.locus, m.repairs, m.promote, threads, 0, 0, which == 0 ? origins - 1 : 0,
                                     which == 1 ? piece_hits - 1 : 0);
            } catch (const panvar::ResourceExhausted& e) {
                refused = true;
                message = e.what();
            }
            const std::string cap = which == 0 ? "physical-origin cap" : "piece-hit cap";
            require(refused, at + ": one below the " + cap + " must refuse");
            require(message.find(cap) != std::string::npos,
                    at + ": the refusal must name the " + cap + ", not a cancellation: " + message);
        }
    }
    std::cout << "  ok: 1/2/4/8 threads accept exactly at the physical-origin (" << origins << ") and piece-hit ("
              << piece_hits << ") totals, identically, and refuse one below naming the crossed cap\n";
}

void test_active_workers_never_exceed_threads() {
    const MultiBlock m = multi_block(7171);
    for (const std::size_t threads : {1u, 2u, 4u, 8u}) {
        panvar::reset_peak_active_workers();
        panvar::FragmentScoringOptions opts = options();
        opts.threads = threads;
        panvar::PreparedScoringContext context(m.locus.fragments, library(), opts);
        context.begin_transient_round();
        std::vector<std::string> seeds{m.locus.base};
        for (int v = 0; v < 5; ++v) {   // several haplotypes, so placement runs in parallel too
            std::string other = m.locus.base;
            other[700 + 900 * v] = other[700 + 900 * v] == 'A' ? 'C' : 'A';
            seeds.push_back(other);
        }
        const std::vector<std::size_t> slots = context.add_haplotypes(seeds);   // placement stage
        context.retain_census(slots[0]);
        std::vector<panvar::HaplotypeCandidate> pairs;
        for (const Repair& r : m.repairs) {                                     // enumeration stage
            pairs.push_back({slots[0], context.add_repaired_haplotype(slots[0], r.begin, r.end, r.allele)});
        }
        (void)context.score_candidates(pairs);                                  // combination stage
        const std::size_t peak = panvar::peak_active_workers();
        require(peak <= threads, std::to_string(threads) + " threads requested but " + std::to_string(peak) +
                                     " workers were active at once");
        if (threads > 1) require(peak >= 2, std::to_string(threads) + " threads: no stage ran in parallel");
        std::cout << "  ok: --threads " << threads << ": peak " << peak
                  << " active workers across placement, enumeration and combination\n";
    }
}

void test_cumulative_caps_through_the_delta() {
    const Locus l = locus(99, {400});
    std::mt19937_64 rng(11);
    std::vector<std::string> replacements;
    for (int i = 0; i < 4; ++i) replacements.push_back(random_bases(rng, 40));

    std::uint64_t total = 0;
    {
        panvar::PreparedScoringContext context(l.fragments, library(), options());
        context.add_haplotypes({l.base});
        context.begin_transient_round();
        context.retain_census(0);
        for (const std::string& replacement : replacements) {
            (void)context.add_repaired_haplotype(0, 400, 440, replacement);
        }
        total = context.cumulative_counters().physical_origins;
    }
    require(total > 0, "the fixture must produce origins");

    {
        panvar::FragmentScoringOptions at_cap = options();
        at_cap.max_physical_origins = total;
        panvar::PreparedScoringContext context(l.fragments, library(), at_cap);
        context.add_haplotypes({l.base});
        context.begin_transient_round();
        context.retain_census(0);
        for (const std::string& replacement : replacements) {
            (void)context.add_repaired_haplotype(0, 400, 440, replacement);
        }
        require(!context.refused(), "exactly at the cumulative cap must be accepted");
    }
    {
        panvar::FragmentScoringOptions under = options();
        under.max_physical_origins = total - 1;
        panvar::PreparedScoringContext context(l.fragments, library(), under);
        bool refused = false;
        try {
            context.add_haplotypes({l.base});
            context.begin_transient_round();
            context.retain_census(0);
            for (const std::string& replacement : replacements) {
                (void)context.add_repaired_haplotype(0, 400, 440, replacement);
            }
        } catch (const panvar::ResourceExhausted&) {
            refused = true;
        }
        require(refused, "one below the cumulative cap must refuse");
        require(context.refused(), "and must invalidate the whole search");
        require_throws<panvar::ResourceExhausted>(
            [&]() { (void)context.score_candidates({{0, 0}}); },
            "a refused context must not emit a partial ranking");
    }
    std::cout << "  ok: cumulative cap " << total << " accepts, one below refuses with no ranking\n";
}

}  // namespace


// ---- I2-B (experimental): round-local fine discovery in the adaptive window enumeration ------------------------------
// Rounds go through add_repaired_haplotypes (the production batch); a winner is promoted by delta with its census kept,
// the rest released, and later rounds deltas from it into reused slots. With the per-window row check on, every window
// the fine path builds is compared row for row (start, consumed, mismatches, emission bits) with the reference's own two
// scans over the whole window; the cells, exposures and score components (audit views included) must equal the
// reference context's bit for bit.
struct FineHow {
    bool fine = false;
    std::size_t budget = std::size_t{1} << 28;
    std::uint64_t served_cap = std::uint64_t{1} << 40;
    std::size_t threads = 1;
    unsigned order_seed = 0;
    bool origin = true;
    std::size_t origin_budget = std::size_t{1} << 20;
    bool audit = false;
    bool repeat = false;
    std::uint64_t max_physical_origins = 0;
    panvar::LibraryModel library_model = library();
};
struct FineResult {
    std::vector<double> cells;
    std::vector<panvar::DiplotypeScore> scores;
    std::vector<std::vector<panvar::DiplotypeScore>> audits;
    std::vector<std::size_t> slots;
    panvar::FineWindowReport report;
    std::size_t unique = 0, created = 0;
};

FineResult fine_search(const Locus& l, const std::vector<Repair>& repairs, std::size_t promote, const FineHow& how) {
    panvar::FragmentScoringOptions opts = options();
    opts.threads = how.threads;
    if (how.max_physical_origins != 0) opts.max_physical_origins = how.max_physical_origins;
    if (how.audit) {
        panvar::LibraryModel a = how.library_model, b = how.library_model;
        a.fragment_start_intensity *= 2.0;
        b.outlier_weight = 0.02;
        opts.audit_libraries = {a, b};
    }
    panvar::PreparedScoringContext context(l.fragments, how.library_model, opts);
    context.add_haplotypes({l.base});
    context.begin_transient_round();
    context.retain_census(0);
    if (how.origin) {
        context.enable_origin_cache_path(true);
        context.set_origin_cache_limits(how.origin_budget, 0);
    }
    if (how.fine) context.enable_fine_window_discovery(how.budget, how.served_cap);
    std::vector<std::size_t> order(repairs.size());
    std::iota(order.begin(), order.end(), std::size_t{0});
    if (how.order_seed != 0) {
        std::mt19937_64 rng(how.order_seed);
        std::shuffle(order.begin(), order.end(), rng);
    }
    FineResult out;
    const auto round = [&](std::size_t incumbent, bool skip_promoted) {
        std::vector<panvar::RepairRequest> requests;
        std::vector<std::size_t> which;
        for (const std::size_t i : order) {
            if (skip_promoted && i == promote) continue;
            requests.push_back({incumbent, repairs[i].begin, repairs[i].end, repairs[i].allele});
            which.push_back(i);
        }
        if (how.repeat && !requests.empty()) {   // the same request twice in one round
            requests.push_back(requests.front());
            which.push_back(which.front());
        }
        const std::vector<std::size_t> slots = context.add_repaired_haplotypes(requests);
        // recorded in repair order, so permuted requests compare position by position
        std::vector<std::size_t> slot_of(repairs.size(), SIZE_MAX);
        for (std::size_t k = 0; k < which.size(); ++k) slot_of[which[k]] = slots[k];
        for (std::size_t i = 0; i < repairs.size(); ++i) {
            if (slot_of[i] == SIZE_MAX) continue;
            const std::size_t slot = slot_of[i];
            for (std::size_t f = 0; f < l.fragments.size(); ++f) out.cells.push_back(context.mass_of(f, slot));
            out.cells.push_back(context.exposure_of(slot));
            const panvar::FragmentScoringResult r = context.score_candidates({{incumbent, slot}});
            out.scores.push_back(r.scores[0]);
            for (std::size_t k = 0; k < r.audit_scores.size(); ++k) {
                if (out.audits.size() <= k) out.audits.resize(k + 1);
                out.audits[k].push_back(r.audit_scores[k][0]);
            }
        }
    };
    round(0, false);
    const Repair& chosen = repairs[promote];
    const std::size_t winner = context.add_repaired_haplotype(0, chosen.begin, chosen.end, chosen.allele, nullptr, true);
    context.release_transient({winner});
    for (int r = 0; r < 2; ++r) {
        round(winner, true);
        context.release_transient({winner});
    }
    out.report = context.fine_window_report();
    out.unique = context.unique_haplotypes();
    out.created = context.haplotypes_created();
    return out;
}

void require_same_fine(const FineResult& got, const FineResult& want, const std::string& label) {
    require(got.cells.size() == want.cells.size() && got.scores.size() == want.scores.size() &&
                got.audits.size() == want.audits.size(),
            label + ": result shapes differ");
    for (std::size_t i = 0; i < want.cells.size(); ++i) require_identical(got.cells[i], want.cells[i], label + ": cell " + std::to_string(i));
    const auto same = [&](const panvar::DiplotypeScore& a, const panvar::DiplotypeScore& b, const std::string& what) {
        require_identical(a.log_likelihood, b.log_likelihood, what + " log likelihood");
        require_identical(a.fragment_term, b.fragment_term, what + " fragment term");
        require_identical(a.exposure_term, b.exposure_term, what + " exposure term");
        require_identical(a.absolute_term_sum, b.absolute_term_sum, what + " term sum");
        require(a.fragments_with_mass == b.fragments_with_mass, what + " fragments with mass");
    };
    for (std::size_t i = 0; i < want.scores.size(); ++i) same(got.scores[i], want.scores[i], label + ": score " + std::to_string(i));
    for (std::size_t k = 0; k < want.audits.size(); ++k)
        for (std::size_t i = 0; i < want.audits[k].size(); ++i)
            same(got.audits[k][i], want.audits[k][i], label + ": audit " + std::to_string(k) + "/" + std::to_string(i));
    require(got.unique == want.unique && got.created == want.created, label + ": slot reuse differs");
}

// Repeated copies, both mate orientations, mixed read lengths, blocks at both locus ends, tiny and empty blocks.
MultiBlock fine_fixture(std::uint64_t seed) {
    std::mt19937_64 rng(seed);
    MultiBlock out;
    const std::string unit = random_bases(rng, 60);
    std::string base = random_bases(rng, 1400);
    for (int i = 0; i < 8; ++i) base += unit;
    base += random_bases(rng, 1400);
    out.locus.base = base;
    std::size_t n = 0;
    for (std::size_t start = 0; start + 260 <= base.size(); start += 9) {
        panvar::ReadPair pair = fragment_at(base, start, 260, 100, "x" + std::to_string(n));
        if (n % 2 == 1) std::swap(pair.first, pair.second);              // both orientations
        if (n % 7 == 0) {                                                 // mixed read lengths
            pair.second.bases.resize(93);
            pair.second.qualities.resize(93);
        }
        out.locus.fragments.push_back(pair);
        ++n;
    }
    const std::size_t array = 1400;
    out.repairs = {
        {0, 6, random_bases(rng, 6)},                       // the locus start
        {2, 3, ""},                                          // a tiny block emptied
        {array + 60, array + 120, random_bases(rng, 60)},   // a copy of the repeat replaced
        {array + 60, array + 120, ""},                       // one copy deleted
        {array + 60, array + 120, unit + unit},              // one copy duplicated
        {700, 700, random_bases(rng, 37)},                   // an empty block filled (a new junction)
        {900, 950, random_bases(rng, 50)},
        {900, 950, random_bases(rng, 3)},
        {base.size() - 5, base.size(), random_bases(rng, 9)},   // the locus end
    };
    out.promote = 6;   // same length
    return out;
}

void test_fine_window_discovery() {
    panvar::set_fine_window_row_check(true);
    std::size_t checked = 0;
    for (const std::uint64_t seed : {2101u, 2102u}) {
        const MultiBlock m = seed == 2101u ? multi_block(seed) : fine_fixture(seed);
        const std::string fixture = seed == 2101u ? "multi-block" : "repeat/orientation/ends";
        for (const bool origin : {true, false}) {
            FineHow ref;
            ref.origin = origin;
            ref.audit = true;
            const FineResult want = fine_search(m.locus, m.repairs, m.promote, ref);
            require(!want.report.enabled && want.report.enumerations > 0, "fixture: the reference enumerates windows");
            for (const std::size_t threads : {std::size_t{1}, std::size_t{4}}) {
                FineHow how = ref;
                how.fine = true;
                how.threads = threads;
                FineHow ref_t = ref;
                ref_t.threads = threads;
                const FineResult base_t = fine_search(m.locus, m.repairs, m.promote, ref_t);
                require_same_fine(base_t, want, fixture + ": reference, " + std::to_string(threads) + " threads");
                const FineResult got = fine_search(m.locus, m.repairs, m.promote, how);
                const std::string label = fixture + (origin ? ", origin path" : ", stage 5") + ", " + std::to_string(threads) + " threads";
                require_same_fine(got, want, label);
                const panvar::FineWindowReport& r = got.report;
                require(r.enabled && r.rounds == 3 && r.units_discovered > 0 && r.units_served > 0 && r.placements_served > 0 &&
                            r.fallback_windows == 0 && r.units_not_retained == 0 && r.round_peak_held_bytes <= how.budget,
                        label + ": the fine path must discover, then serve shared units, within its budget");
                checked += r.worker_windows;
                if (threads == 1 && origin) {
                    std::cout << "  ok: " << fixture << " -- " << r.worker_windows << " windows, " << r.units_discovered
                              << " units discovered, " << r.units_served << " served (" << r.placements_served
                              << " placements served of " << r.placements_served + r.placements_discovered
                              << "), round peak " << r.round_peak_bytes << " bytes\n";
                }
            }
            // permuted and repeated requests
            for (const unsigned order_seed : {5u, 17u}) {
                FineHow how = ref;
                how.fine = true;
                how.order_seed = order_seed;
                how.repeat = true;
                FineHow ref_p = ref;
                ref_p.order_seed = order_seed;
                ref_p.repeat = true;
                require_same_fine(fine_search(m.locus, m.repairs, m.promote, how), fine_search(m.locus, m.repairs, m.promote, ref_p),
                                  fixture + ": permuted and repeated requests");
            }
        }
    }
    std::cout << "  ok: fine = reference (cells, exposures, components, audit views, slot reuse) at 1 and 4 threads, origin "
                 "path and stage 5, permuted and repeated requests; " << checked << " windows checked row for row\n";

    const MultiBlock m = fine_fixture(2103);
    FineHow ref;
    const FineResult want = fine_search(m.locus, m.repairs, m.promote, ref);
    // storage admission: a store too small to hold anything useful closes, and later windows take the reference path
    {
        FineHow how;
        how.fine = true;
        how.budget = 2048;
        const FineResult got = fine_search(m.locus, m.repairs, m.promote, how);
        require_same_fine(got, want, "undersized store");
        require(got.report.units_not_retained > 0 && got.report.fallback_windows > 0 &&
                    got.report.round_peak_held_bytes <= how.budget,
                "an undersized store must decline, fall back explicitly and stay within its budget");
        std::cout << "  ok: undersized store -- " << got.report.units_not_retained << " units not retained, "
                  << got.report.fallback_windows << " windows by the reference path, identical\n";
    }
    // the served-placement cap refuses, typed, and leaves the context refused
    {
        FineHow how;
        how.fine = true;
        how.served_cap = 1;
        bool typed = false;
        try {
            (void)fine_search(m.locus, m.repairs, m.promote, how);
        } catch (const panvar::ResourceExhausted&) {
            typed = true;
        } catch (...) {
        }
        require(typed, "the served-placement refusal must be ResourceExhausted");
    }
    // the existing cumulative caps still refuse, typed, with the fine path on
    {
        FineHow how;
        how.fine = true;
        how.max_physical_origins = 2000;
        bool typed = false;
        try {
            (void)fine_search(m.locus, m.repairs, m.promote, how);
        } catch (const panvar::ResourceExhausted&) {
            typed = true;
        } catch (...) {
        }
        require(typed, "the physical-origin cap must still refuse with the fine path on");
    }
    // independent read/library contexts
    {
        const MultiBlock other = multi_block(2104);
        FineHow a, b;
        b.library_model.insert_sd_bp = 15.0;
        const FineResult ref_a = fine_search(m.locus, m.repairs, m.promote, a);
        const FineResult ref_b = fine_search(other.locus, other.repairs, other.promote, b);
        a.fine = b.fine = true;
        require_same_fine(fine_search(m.locus, m.repairs, m.promote, a), ref_a, "context A");
        require_same_fine(fine_search(other.locus, other.repairs, other.promote, b), ref_b, "context B (other reads and library)");
    }
    // negative controls: broken ownership or context must be caught by the row check
    for (const int mutation : {1, 2, 3, 4}) {
        panvar::set_fine_window_mutation(mutation);
        FineHow how;
        how.fine = true;
        bool caught = false;
        try {
            (void)fine_search(m.locus, m.repairs, m.promote, how);
        } catch (const std::logic_error&) {
            caught = true;
        } catch (...) {
        }
        panvar::set_fine_window_mutation(0);
        require(caught, "mutation " + std::to_string(mutation) + " must be caught by the row check");
    }
    panvar::set_fine_window_row_check(false);
    std::cout << "  ok: served cap and physical-origin cap refuse typed; independent contexts; mutations 1-4 (overlap, "
                 "left context, right context, end window) caught\n";
}


// ---- I2-D: anchor scores and censuses in one placement pass -----------------------------------------------------------
// The fused route (add_haplotypes with retain_censuses) against score-then-retain: identical slots, cells, exposures,
// score components (audit views included) and censuses (every origin's coordinates, orientation and latent mass, in
// order); one full placement per distinct new sequence, none for the later retain_census; executed work halved
// honestly; and everything downstream (delta rounds with and without fine windows, promotion, release, slot reuse)
// identical, at 1 and 4 threads.
struct SeedResult {
    std::vector<std::size_t> slots;
    std::vector<double> cells;
    std::vector<panvar::DiplotypeScore> scores;
    std::vector<std::vector<panvar::HaploidOrigin>> census;   // [slot-ordered][origins]
    std::vector<std::size_t> census_lengths;
    std::uint64_t placements_after_seed = 0, placements_after_retain = 0;
    panvar::FragmentScoringCounters counters_after_retain;
    std::vector<double> downstream;
};

SeedResult seed_route(const Locus& l, const std::string& a, const std::string& b, bool fused, std::size_t threads,
                      const std::vector<Repair>& repairs, bool fine_windows) {
    panvar::FragmentScoringOptions opts = options();
    opts.threads = threads;
    panvar::LibraryModel x = library(), y = library();
    x.fragment_start_intensity *= 2.0;
    y.outlier_weight = 0.02;
    opts.audit_libraries = {x, y};
    panvar::PreparedScoringContext context(l.fragments, library(), opts);
    context.begin_transient_round();
    context.enable_origin_cache_path(true);
    context.set_origin_cache_limits(std::size_t{1} << 20, 0);
    if (fine_windows) context.enable_fine_window_discovery(std::size_t{1} << 26, std::uint64_t{1} << 40);
    SeedResult out;
    out.slots = context.add_haplotypes({a, b}, nullptr, false, fused);
    out.placements_after_seed = context.full_placement_count();
    for (const std::size_t s : out.slots) context.retain_census(s);
    out.placements_after_retain = context.full_placement_count();
    out.counters_after_retain = context.cumulative_counters();
    std::set<std::size_t> distinct(out.slots.begin(), out.slots.end());
    for (const std::size_t s : distinct) {
        for (std::size_t f = 0; f < l.fragments.size(); ++f) out.cells.push_back(context.mass_of(f, s));
        out.cells.push_back(context.exposure_of(s));
        const panvar::HaploidOriginCensus* c = context.census_view(s);
        require(c != nullptr, "every seed slot must hold a census");
        out.census_lengths.push_back(c->haplotype_length);
        for (const auto& per_fragment : c->origins_by_fragment) out.census.push_back(per_fragment);
    }
    const panvar::FragmentScoringResult r = context.score_candidates({{out.slots[0], out.slots[1]}});
    out.scores.push_back(r.scores[0]);
    for (const auto& audit : r.audit_scores) out.scores.push_back(audit[0]);
    // downstream: a round of deltas from the first anchor, a promotion with census, release, a second round into reused slots
    const auto round = [&](std::size_t incumbent) {
        std::vector<panvar::RepairRequest> requests;
        for (const Repair& rep : repairs) requests.push_back({incumbent, rep.begin, rep.end, rep.allele});
        const std::vector<std::size_t> got = context.add_repaired_haplotypes(requests);
        for (const std::size_t s : got) {
            for (std::size_t f = 0; f < l.fragments.size(); ++f) out.downstream.push_back(context.mass_of(f, s));
            out.downstream.push_back(context.score_candidates({{incumbent, s}}).scores[0].log_likelihood);
        }
    };
    round(out.slots[0]);
    const std::size_t winner = context.add_repaired_haplotype(out.slots[0], repairs[0].begin, repairs[0].end, repairs[0].allele, nullptr, true);
    context.release_transient({winner, out.slots[1]});
    std::vector<Repair> later(repairs.begin() + 1, repairs.end());
    std::vector<panvar::RepairRequest> requests;
    for (const Repair& rep : later) requests.push_back({winner, rep.begin, rep.end, rep.allele});
    const std::vector<std::size_t> got = context.add_repaired_haplotypes(requests);
    for (const std::size_t s : got) {
        for (std::size_t f = 0; f < l.fragments.size(); ++f) out.downstream.push_back(context.mass_of(f, s));
    }
    return out;
}

void require_same_seed(const SeedResult& got, const SeedResult& want, const std::string& label) {
    require(got.slots == want.slots, label + ": slots differ");
    require(got.cells.size() == want.cells.size() && got.downstream.size() == want.downstream.size(), label + ": shapes differ");
    for (std::size_t i = 0; i < want.cells.size(); ++i) require_identical(got.cells[i], want.cells[i], label + ": cell " + std::to_string(i));
    for (std::size_t i = 0; i < want.scores.size(); ++i) {
        require_identical(got.scores[i].log_likelihood, want.scores[i].log_likelihood, label + ": score");
        require_identical(got.scores[i].fragment_term, want.scores[i].fragment_term, label + ": fragment term");
        require_identical(got.scores[i].exposure_term, want.scores[i].exposure_term, label + ": exposure term");
        require_identical(got.scores[i].absolute_term_sum, want.scores[i].absolute_term_sum, label + ": term sum");
        require(got.scores[i].fragments_with_mass == want.scores[i].fragments_with_mass, label + ": fragments with mass");
    }
    require(got.census_lengths == want.census_lengths && got.census.size() == want.census.size(), label + ": census shapes differ");
    for (std::size_t i = 0; i < want.census.size(); ++i) {
        require(got.census[i].size() == want.census[i].size(), label + ": census origin counts differ");
        for (std::size_t k = 0; k < want.census[i].size(); ++k) {
            const panvar::HaploidOrigin& p = got.census[i][k];
            const panvar::HaploidOrigin& q = want.census[i][k];
            require(p.forward_start == q.forward_start && p.insert_bp == q.insert_bp && p.orientation == q.orientation &&
                        std::memcmp(&p.log_latent_mass, &q.log_latent_mass, sizeof(double)) == 0,
                    label + ": a census origin differs");
        }
    }
    for (std::size_t i = 0; i < want.downstream.size(); ++i) {
        require_identical(got.downstream[i], want.downstream[i], label + ": downstream value " + std::to_string(i));
    }
}

void test_fused_seed() {
    std::size_t compared = 0;
    for (const std::uint64_t seed : {2201u, 2202u}) {
        const MultiBlock m = seed == 2201u ? multi_block(seed) : fine_fixture(seed);
        const std::string fixture = seed == 2201u ? "multi-block" : "repeat/orientation/mixed lengths/ends/tiny";
        std::string other = m.locus.base;
        other.replace(m.repairs[1].begin, m.repairs[1].end - m.repairs[1].begin, m.repairs[1].allele);
        for (const bool identical : {false, true}) {
            const std::string& b = identical ? m.locus.base : other;
            for (const bool fine_windows : {false, true}) {
                const SeedResult want = seed_route(m.locus, m.locus.base, b, false, 1, m.repairs, fine_windows);
                for (const std::size_t threads : {std::size_t{1}, std::size_t{4}}) {
                    const std::string label = fixture + (identical ? ", homozygous" : ", distinct") +
                                              (fine_windows ? ", fine windows" : "") + ", " + std::to_string(threads) + " threads";
                    const SeedResult ref = seed_route(m.locus, m.locus.base, b, false, threads, m.repairs, fine_windows);
                    const SeedResult got = seed_route(m.locus, m.locus.base, b, true, threads, m.repairs, fine_windows);
                    require_same_seed(ref, want, label + " (reference, threads)");
                    require_same_seed(got, want, label);
                    const std::uint64_t distinct = identical ? 1 : 2;
                    require(got.slots[0] == got.slots[1] ? identical : !identical, label + ": slot sharing must follow byte identity");
                    require(ref.placements_after_seed == distinct && ref.placements_after_retain == 2 * distinct,
                            label + ": the reference route places each distinct anchor twice");
                    require(got.placements_after_seed == distinct && got.placements_after_retain == distinct,
                            label + ": the fused route places each distinct anchor once; retain_census then places nothing");
                    require(2 * got.counters_after_retain.physical_origins == ref.counters_after_retain.physical_origins &&
                                2 * got.counters_after_retain.piece_hits == ref.counters_after_retain.piece_hits,
                            label + ": executed work must halve honestly");
                    ++compared;
                }
            }
        }
    }
    // an already-held slot without a census, beside a new one, with released-slot reuse (the JIT materialisation shape)
    {
        const MultiBlock m = fine_fixture(2203);
        std::string c = m.locus.base;
        c.replace(m.repairs[6].begin, m.repairs[6].end - m.repairs[6].begin, m.repairs[6].allele);
        const auto run = [&](bool fused) {
            panvar::PreparedScoringContext context(m.locus.fragments, library(), options());
            context.begin_transient_round();
            const std::size_t held = context.add_haplotypes({m.locus.base}).front();
            const std::size_t scratch = context.add_haplotypes({c}).front();
            context.release_transient({held});   // frees `scratch` for reuse
            std::vector<std::size_t> slots;
            if (fused) {
                slots = context.add_haplotypes({m.locus.base, c}, nullptr, true, true);
            } else {
                slots = context.add_haplotypes({m.locus.base, c}, nullptr, true);
                for (const std::size_t s : slots) context.retain_census(s);
            }
            require(slots[0] == held && slots[1] == scratch, "an existing slot is reused by identity, a released one by reuse");
            std::vector<double> cells;
            for (const std::size_t s : slots) {
                for (std::size_t f = 0; f < m.locus.fragments.size(); ++f) cells.push_back(context.mass_of(f, s));
                require(context.holds_census(s), "every requested slot must hold a census");
                for (const auto& v : context.census_view(s)->origins_by_fragment) cells.push_back(static_cast<double>(v.size()));
            }
            cells.push_back(static_cast<double>(context.full_placement_count()));
            return cells;
        };
        const std::vector<double> a = run(false), b = run(true);
        require(a.size() == b.size(), "existing/reused slots: shapes differ");
        for (std::size_t i = 0; i + 1 < a.size(); ++i) require_identical(b[i], a[i], "existing/reused slots: value " + std::to_string(i));
        require(a.back() == 5 && b.back() == 4, "the held slot is censused once by placement; the reused slot's census comes from its own pass");
    }
    // I2-E: a batch with a NEW slot and an already-HELD slot without a census, where the held slot's fallback census
    // placement crosses the cap: typed, the context permanently refused, and no census from the failed batch owned --
    // neither the new slot's captured census nor any fallback census.
    {
        const MultiBlock m = fine_fixture(2205);
        std::string fresh = m.locus.base;
        fresh.replace(m.repairs[6].begin, m.repairs[6].end - m.repairs[6].begin, m.repairs[6].allele);
        // measure the work: the held slot's placement (X) and the new slot's (Y), uncapped
        std::uint64_t x = 0, y = 0;
        {
            panvar::PreparedScoringContext probe(m.locus.fragments, library(), options());
            (void)probe.add_haplotypes({m.locus.base});
            x = probe.cumulative_counters().physical_origins;
            (void)probe.add_haplotypes({fresh});
            y = probe.cumulative_counters().physical_origins - x;
        }
        require(x > 0 && y > 0, "fixture: both placements must enumerate origins");
        panvar::FragmentScoringOptions opts = options();
        opts.max_physical_origins = x + y + x - 1;   // the first pass fits (x + y), the fallback census (+x) does not
        panvar::PreparedScoringContext context(m.locus.fragments, library(), opts);
        context.begin_transient_round();
        const std::size_t held = context.add_haplotypes({m.locus.base}).front();
        bool typed = false;
        try {
            (void)context.add_haplotypes({m.locus.base, fresh}, nullptr, false, true);
        } catch (const panvar::ResourceExhausted&) {
            typed = true;
        }
        require(typed && context.refused(), "a fallback census placement crossing the cap must refuse, typed, permanently");
        require(!context.holds_census(held), "the failed fallback must leave no census for the held slot");
        for (std::size_t slot = 0; slot < context.unique_haplotypes(); ++slot) {
            require(!context.holds_census(slot), "no census from the failed batch may be owned (slot " + std::to_string(slot) + ")");
        }
        bool refused_again = false;
        try {
            (void)context.add_haplotypes({fresh});
        } catch (const panvar::ResourceExhausted&) {
            refused_again = true;
        }
        require(refused_again, "the refused context must refuse every later call");
    }
    // a cap refusal during the fused pass: typed, the context refused, no census left behind
    {
        const MultiBlock m = fine_fixture(2204);
        panvar::FragmentScoringOptions opts = options();
        opts.max_physical_origins = 1000;
        panvar::PreparedScoringContext context(m.locus.fragments, library(), opts);
        bool typed = false;
        try {
            (void)context.add_haplotypes({m.locus.base}, nullptr, false, true);
        } catch (const panvar::ResourceExhausted&) {
            typed = true;
        }
        require(typed && context.refused(), "a fused refusal must be typed and refuse the context");
        require(!context.holds_census(0), "a refused fused pass must leave no census");
        bool refused_again = false;
        try {
            context.retain_census(0);
        } catch (const panvar::ResourceExhausted&) {
            refused_again = true;
        }
        require(refused_again, "a refused context must refuse every later call");
    }
    std::cout << "  ok: fused seeding = score-then-retain in " << compared << " configurations (cells, exposures, components, "
                 "audit views, every census origin, downstream rounds with and without fine windows, promotion, release, slot "
                 "reuse); one placement per distinct anchor, none at retain; executed work halved; held/reused slots; refusal "
                 "leaves no census\n";
}

// S2-I: the partition storage memory_account adds (each worker's ordinary WindowScratch and its forward/reverse
// PieceIndex vectors; the fine stores stay in the existing window-scratch ledger) equals an independent recount, with
// and without fine window discovery, and the extended total adds exactly the two capacities.
void test_partition_storage_accounts() {
    const Locus t = tandem_locus(5153);
    std::mt19937_64 rng(5154);
    for (const bool fine : {false, true}) {
        panvar::FragmentScoringOptions opts = options();
        opts.threads = 4;
        panvar::PreparedScoringContext context(t.fragments, library(), opts);
        context.begin_transient_round();
        context.add_haplotypes({t.base});
        context.retain_census(0);
        context.enable_origin_cache_path(true);
        context.set_origin_cache_limits(std::size_t{1} << 30, 0);
        if (fine) context.enable_fine_window_discovery(std::size_t{1} << 26, 2000000000ULL);
        std::vector<panvar::RepairRequest> requests = {{0, 1000, 1040, random_bases(rng, 40)}, {0, 1100, 1140, random_bases(rng, 55)},
                                                      {0, 1300, 1340, std::string()}};
        std::vector<panvar::BlockDeltaWork> works;
        (void)context.add_repaired_haplotypes(requests, &works, nullptr);
        const panvar::ContextMemoryAccount a = context.memory_account();
        const auto recount = context.recount_partition_storage_for_testing();
        const std::string at = fine ? "fine discovery" : "reference discovery";
        require(a.partitions > 1, at + ": the repairs enumerated across partitions (" + std::to_string(a.partitions) + ")");
        require(a.partition_scratch_capacity == recount.first && a.partition_piece_index_capacity == recount.second,
                at + ": partition capacities equal the independent recount (" + std::to_string(a.partition_scratch_capacity) + "/" +
                    std::to_string(recount.first) + ", " + std::to_string(a.partition_piece_index_capacity) + "/" +
                    std::to_string(recount.second) + ")");
        require(a.partition_piece_index_bytes > 0 && a.partition_piece_index_bytes <= a.partition_piece_index_capacity &&
                    a.partition_scratch_bytes <= a.partition_scratch_capacity && a.partition_scratch_capacity > 0,
                at + ": sizes never exceed capacities");
        require(a.capacity_total_extended == a.capacity_total + a.partition_scratch_capacity + a.partition_piece_index_capacity,
                at + ": the extended total adds exactly the two partition capacities (fine stores not counted twice)");
        std::cout << "  " << at << ": partitions " << a.partitions << ", scratch " << a.partition_scratch_bytes << " / "
                  << a.partition_scratch_capacity << " B, piece indexes " << a.partition_piece_index_bytes << " / "
                  << a.partition_piece_index_capacity << " B\n";
    }
}

int main() {
    if (const char* only = std::getenv("PANVAR_BLOCK_DELTA_ONLY"); only != nullptr && std::string(only) == "partition_accounts") {
        try {
            test_partition_storage_accounts();
            std::cout << "block delta (partition_accounts only): PASS\n";
            return 0;
        } catch (const std::exception& e) {
            std::cerr << "block delta (partition_accounts only): FAIL: " << e.what() << '\n';
            return 1;
        }
    }
    try {
        std::cout << "== I2-D: fused anchor seeding ==\n";
        test_fused_seed();
        std::cout << "== I2-B: fine window discovery in the adaptive search ==\n";
        test_fine_window_discovery();
        std::cout << "== substitutions, insertions, deletions, empty alleles ==\n";
        test_replacement_kinds();
        std::cout << "== first block, last block, length-changing ==\n";
        test_block_positions();
        std::cout << "== adjacent short blocks ==\n";
        test_adjacent_short_blocks();
        std::cout << "== tandem-repeat multiplicity and StateKey collapse ==\n";
        test_tandem_repeat_multiplicity();
        std::cout << "== gapped alignments with multiple consumed-span components ==\n";
        test_gapped_components();
        std::cout << "== random small-locus sweep ==\n";
        test_random_sweep();
        std::cout << "== preconditions ==\n";
        test_preconditions_refuse();
        std::cout << "== transient destruction, then cache hits ==\n";
        test_transient_destruction_then_cache_hits();
        std::cout << "== transient release ==\n";
        test_transient_release();
        std::cout << "== scoring after a release ==\n";
        test_scoring_after_release();
        std::cout << "== released slots are reused ==\n";
        test_released_slots_are_reused();
        std::cout << "== a released seed census, respelled by delta, is bit-identical ==\n";
        test_released_seed_is_respelled_identically();
        std::cout << "== planned retention: scoring identical at every budget ==\n";
        test_planned_cache_scores_identically();
        std::cout << "== planned retention: independent of request order and threads ==\n";
        test_planned_cache_order_independent();
        std::cout << "== planned retention: eviction releases and reuses storage ==\n";
        test_planned_eviction_releases_storage();
        std::cout << "== planned retention within a round: adversarial request order ==\n";
        test_batch_retention_is_order_independent();
        std::cout << "== piece screen: a necessary condition, with offset consistency ==\n";
        test_piece_screen_is_a_necessary_condition();
        std::cout << "== window enumeration across workers: bit-identical to one thread ==\n";
        test_parallel_enumeration_is_exact();
        std::cout << "== window enumeration across workers: caps exact at 1/2/4/8 threads ==\n";
        test_parallel_caps_exact();
        std::cout << "== active workers never exceed --threads ==\n";
        test_active_workers_never_exceed_threads();
        std::cout << "== stage 6: origin path versus stage 5 (self-checked) ==\n";
        test_origin_path_matches_stage5();
        std::cout << "== stage 6: disabled, full and undersized origin caches ==\n";
        test_origin_budgets_identical();
        std::cout << "== stage 6: forced origin-cache collisions ==\n";
        test_origin_forced_collisions();
        std::cout << "== stage 6: identical windows at different offsets ==\n";
        test_origin_translation_across_offsets();
        std::cout << "== stage 6: end flags never alias ==\n";
        test_origin_end_flags_never_alias();
        std::cout << "== stage 6: length shifts keep preaggregations apart ==\n";
        test_origin_length_shifts_keep_groups_apart();
        std::cout << "== stage 6: caps count executed work; consumption is capped separately ==\n";
        test_origin_caps_count_executed_work();
        std::cout << "== stage 6: thread identity ==\n";
        test_origin_threads_identical();
        std::cout << "== warm branch seeding: delta-promoted census equals cold full seeding ==\n";
        test_delta_promoted_seed_matches_cold_seed();
        std::cout << "== grouped enumeration versus the reference ==\n";
        test_window_enumeration_matches_reference();
        std::cout << "== scratch reuse and per-candidate allocation ==\n";
        test_scratch_does_not_allocate_per_candidate();
        std::cout << "== arena slab growth and cross-slab equality ==\n";
        test_arena_slabs_grow_and_compare_across_slabs();
        std::cout << "== allocations on the hit path ==\n";
        test_hits_allocate_nothing();
        std::cout << "== segmented versus materialised windows ==\n";
        test_segmented_windows_match_materialised();
        std::cout << "== path-distinct, byte-identical haplotypes ==\n";
        test_path_distinct_haplotypes_share_one_column();
        std::cout << "== exact cache-budget boundary ==\n";
        test_cache_budget_exact_boundary();
        std::cout << "== no dangling reads after destruction ==\n";
        test_no_dangling_reads_after_destruction();
        std::cout << "== forced cache-hash collisions ==\n";
        test_forced_hash_collisions();
        std::cout << "== ambiguous bases ==\n";
        test_ambiguous_bases_refuse();
        std::cout << "== identity after transient destruction ==\n";
        test_identical_after_transient_destruction();
        std::cout << "== cumulative caps through the delta ==\n";
        test_cumulative_caps_through_the_delta();
        std::cout << "\nblock-delta gates passed\n";
    } catch (const std::exception& error) {
        std::cerr << "\nBLOCK-DELTA GATE FAILED: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
