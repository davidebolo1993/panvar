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
#include <cstdlib>
#include <functional>
#include <iostream>
#include <new>
#include <random>
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

int main() {
    try {
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
