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
        context.freeze_emission_cache();
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
    // Not frozen yet: scoring a transient haplotype would let its bytes into the pointer cache.
    require_throws<std::invalid_argument>(
        [&]() { (void)context.add_repaired_haplotype(0, 400, 440, replacement); },
        "an unfrozen cache must refuse a transient repair haplotype");
    context.freeze_emission_cache();
    // Frozen, but no census: refusing is the point. Placing in full here would report a delta that
    // never happened.
    require_throws<std::invalid_argument>(
        [&]() { (void)context.add_repaired_haplotype(0, 400, 440, replacement); },
        "a missing census must refuse rather than silently place in full");
    context.retain_census(0);
    panvar::BlockDeltaWork work;
    (void)context.add_repaired_haplotype(0, 400, 440, replacement, &work);
    require(work.origins_reused > 0, "the delta should run once its preconditions hold");
    std::cout << "  ok: an unfrozen cache and a missing census both refuse\n";
}

void test_transient_release() {
    const Locus l = locus(88, {400});
    std::mt19937_64 rng(10);
    panvar::PreparedScoringContext context(l.fragments, library(), options());
    context.add_haplotypes({l.base});
    context.freeze_emission_cache();
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
    context.freeze_emission_cache();
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

void test_cumulative_caps_through_the_delta() {
    const Locus l = locus(99, {400});
    std::mt19937_64 rng(11);
    std::vector<std::string> replacements;
    for (int i = 0; i < 4; ++i) replacements.push_back(random_bases(rng, 40));

    std::uint64_t total = 0;
    {
        panvar::PreparedScoringContext context(l.fragments, library(), options());
        context.add_haplotypes({l.base});
        context.freeze_emission_cache();
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
        context.freeze_emission_cache();
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
            context.freeze_emission_cache();
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
        std::cout << "== transient release ==\n";
        test_transient_release();
        std::cout << "== scoring after a release ==\n";
        test_scoring_after_release();
        std::cout << "== cumulative caps through the delta ==\n";
        test_cumulative_caps_through_the_delta();
        std::cout << "\nblock-delta gates passed\n";
    } catch (const std::exception& error) {
        std::cerr << "\nBLOCK-DELTA GATE FAILED: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
