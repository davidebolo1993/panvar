// Phase 6c. The accelerated scorer must reproduce the exhaustive oracle's mass exactly, not
// approximately: every approximation the archived scorer made -- occurrence caps, placement top-K,
// start binning, midpoint deduplication, best-placement instead of summed-placement -- is a way of
// being slightly wrong here, and each one is unsafe exactly where repeat-origin multiplicity is
// evidence.

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

void require_close(double observed, double expected, double tolerance, const std::string& message) {
    if (observed == kNegInf && expected == kNegInf) return;
    if (!std::isfinite(observed) || !std::isfinite(expected) ||
        std::abs(observed - expected) > tolerance) {
        throw std::runtime_error(message + ": observed=" + std::to_string(observed) +
                                 " expected=" + std::to_string(expected));
    }
}

void require_throws(const std::function<void()>& operation, const std::string& message) {
    try {
        operation();
    } catch (const std::invalid_argument&) {
        return;
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

panvar::SequencingRead read(const std::string& name, const std::string& bases) {
    return {name, bases, std::string(bases.size(), 'I')};
}

panvar::LibraryModel library() {
    panvar::LibraryModel model;
    model.insert_mean_bp = 120.0;
    model.insert_sd_bp = 15.0;
    model.support_sigmas = 3.0;
    model.substitution_error = 0.02;
    model.outlier_weight = 0.005;
    model.fragment_start_intensity = 0.01;
    model.provenance = "fixture: declared";
    return model;
}

panvar::InternedHaplotypes intern(const std::vector<std::string>& sequences) {
    std::vector<panvar::MosaicPath> paths;
    for (std::size_t i = 0; i < sequences.size(); ++i) {
        paths.push_back(panvar::MosaicPath{{i}});
    }
    return panvar::intern_haplotypes(paths, sequences);
}

std::string random_bases(std::mt19937_64& rng, std::size_t length) {
    static const char* alphabet = "ACGT";
    std::string out(length, 'A');
    for (char& base : out) base = alphabet[rng() % 4];
    return out;
}

// -------------------------------------------------------------------------------------------
// Gate 1: exact agreement with the exhaustive oracle.
// -------------------------------------------------------------------------------------------

void certify_against_oracle(
    const std::vector<std::string>& haplotypes,
    const std::vector<panvar::ReadPair>& fragments,
    std::size_t band,
    const std::string& label) {
    const panvar::LibraryModel model = library();
    const panvar::InternedHaplotypes interned = intern(haplotypes);

    std::vector<panvar::HaplotypeCandidate> candidates;
    for (std::size_t a = 0; a < interned.sequences.size(); ++a) {
        for (std::size_t b = a; b < interned.sequences.size(); ++b) {
            candidates.push_back({a, b});
        }
    }
    panvar::FragmentScoringOptions options;
    options.mismatch_band = band;
    options.minimum_piece_bp = 8;
    const panvar::FragmentScoringResult actual =
        panvar::score_diplotypes(interned, candidates, fragments, model, options);

    std::size_t compared = 0;
    std::size_t finite = 0;
    for (std::size_t f = 0; f < fragments.size(); ++f) {
        const panvar::InsertSupport support =
            panvar::fragment_insert_support(fragments[f], model);
        const panvar::DiscreteInsertDistribution insert =
            panvar::build_fragment_insert_distribution(fragments[f], model);
        for (std::size_t h = 0; h < interned.sequences.size(); ++h) {
            const auto placements = panvar::enumerate_placements_exhaustive(
                interned.sequences[h], fragments[f], support, band);
            const double expected =
                panvar::log_placement_mass(placements, fragments[f], insert, model);
            require_close(actual.mass_of(f, h), expected, 1e-9,
                          label + ": accelerated mass must equal the exhaustive mass");
            ++compared;
            if (std::isfinite(expected)) ++finite;
        }
    }
    // A green comparison over rows that are all -inf proves nothing: the fixture must actually
    // place fragments.
    require(compared > 0 && finite >= compared / 4,
            label + ": fixture must place fragments on most haplotypes, or it certifies nothing");
    std::cout << "  " << label << ": " << compared << " (fragment, haplotype) masses agree, "
              << finite << " of them finite\n";
}

// One exact planted fragment, one mismatching, on two related haplotypes.
void test_planted() {
    std::mt19937_64 rng(11);
    const std::string backbone = random_bases(rng, 400);
    std::string variant = backbone;
    variant[180] = variant[180] == 'A' ? 'C' : 'A';
    variant[240] = variant[240] == 'G' ? 'T' : 'G';

    std::vector<panvar::ReadPair> fragments;
    for (std::size_t start : {40, 150, 210, 300}) {
        const std::size_t insert = 120;
        fragments.push_back({read("f" + std::to_string(start) + "/1", backbone.substr(start, 60)),
                             read("f" + std::to_string(start) + "/2",
                                  revcomp(backbone.substr(start + insert - 60, 60)))});
    }
    certify_against_oracle({backbone, variant}, fragments, 2, "planted substitutions");
}

// Repeated units: a fragment has several physical origins, whose multiplicity is evidence.
void test_repeats() {
    std::mt19937_64 rng(29);
    const std::string unit = random_bases(rng, 130);
    const std::string spacer = random_bases(rng, 40);
    const std::string two_copies = spacer + unit + spacer + unit + spacer;
    const std::string three_copies = spacer + unit + spacer + unit + spacer + unit + spacer;

    std::vector<panvar::ReadPair> fragments;
    fragments.push_back({read("r/1", unit.substr(0, 55)),
                         read("r/2", revcomp(unit.substr(65, 55)))});
    fragments.push_back({read("s/1", spacer.substr(0, 30) + unit.substr(0, 25)),
                         read("s/2", revcomp(unit.substr(70, 50)))});
    certify_against_oracle({two_copies, three_copies}, fragments, 1, "tandem repeat origins");

    // The mass must actually grow with copy number, or the fixture is not exercising multiplicity.
    const panvar::InternedHaplotypes interned = intern({two_copies, three_copies});
    panvar::FragmentScoringOptions options;
    options.mismatch_band = 1;
    options.minimum_piece_bp = 8;
    const panvar::FragmentScoringResult result = panvar::score_diplotypes(
        interned, {{0, 0}}, fragments, library(), options);
    require(result.mass_of(0, 1) > result.mass_of(0, 0) + 0.3,
            "three repeat copies must carry more placement mass than two: collapsing repeated "
            "origins to one would erase exactly this signal");
    require(result.counters.physical_origins > result.counters.statistical_state_classes,
            "the repeat fixture must exercise state compression, or the counters prove nothing");
}

// Mixed read lengths in one run: piece-length bucketing must not lose placements.
void test_mixed_read_lengths() {
    std::mt19937_64 rng(97);
    const std::string haplotype = random_bases(rng, 500);
    std::vector<panvar::ReadPair> fragments;
    fragments.push_back({read("a/1", haplotype.substr(50, 70)),
                         read("a/2", revcomp(haplotype.substr(100, 70)))});
    fragments.push_back({read("b/1", haplotype.substr(200, 40)),
                         read("b/2", revcomp(haplotype.substr(280, 40)))});
    // Overlapping mates: a physical state whose floor is the longer mate, not the sum.
    fragments.push_back({read("c/1", haplotype.substr(350, 90)),
                         read("c/2", revcomp(haplotype.substr(350, 90)))});
    certify_against_oracle({haplotype}, fragments, 2, "mixed read lengths");
}

// Randomised sweep: the strongest form of the gate, since it does not depend on which cases the
// author happened to think of.
void test_random_sweep() {
    std::size_t total = 0;
    for (std::uint64_t seed = 1; seed <= 12; ++seed) {
        std::mt19937_64 rng(seed * 7919);
        const std::size_t length = 260 + (rng() % 240);
        std::string haplotype = random_bases(rng, length);
        // A duplicated block gives some fragments genuinely several origins.
        const std::size_t copy_from = rng() % (length / 3);
        haplotype += haplotype.substr(copy_from, 110);
        std::string other = haplotype;
        for (int edit = 0; edit < 6; ++edit) {
            const std::size_t at = rng() % other.size();
            other[at] = "ACGT"[rng() % 4];
        }

        std::vector<panvar::ReadPair> fragments;
        for (int f = 0; f < 5; ++f) {
            const std::size_t insert = 95 + (rng() % 50);
            if (haplotype.size() < insert + 1) continue;
            const std::size_t start = rng() % (haplotype.size() - insert);
            const std::size_t mate = 45 + (rng() % 20);
            if (insert < mate) continue;
            std::string first = haplotype.substr(start, mate);
            std::string second = haplotype.substr(start + insert - mate, mate);
            // Sprinkle substitutions so the band is genuinely exercised.
            for (int e = 0; e < 2; ++e) {
                if (rng() % 2) first[rng() % first.size()] = "ACGT"[rng() % 4];
                if (rng() % 2) second[rng() % second.size()] = "ACGT"[rng() % 4];
            }
            fragments.push_back({read("f/1", first), read("f/2", revcomp(second))});
        }
        if (fragments.empty()) continue;
        certify_against_oracle({haplotype, other}, fragments, 2,
                               "random seed " + std::to_string(seed));
        ++total;
    }
    require(total >= 10, "the random sweep must actually run its seeds");
}

// -------------------------------------------------------------------------------------------
// Gates 2-10 and 13: the scorer's contract, independent of the oracle.
// -------------------------------------------------------------------------------------------

void test_candidate_independence() {
    std::mt19937_64 rng(4242);
    const std::string a = random_bases(rng, 400);
    // The decoys must be divergent enough that reads genuinely stop placing on them, or the
    // candidate-independence gates below would hold trivially because every candidate scores alike.
    std::string b = a;
    for (int edit = 0; edit < 40; ++edit) b[rng() % b.size()] = "ACGT"[rng() % 4];
    std::string c = a;
    for (int edit = 0; edit < 80; ++edit) c[rng() % c.size()] = "ACGT"[rng() % 4];

    std::vector<panvar::ReadPair> fragments;
    for (std::size_t start : {30, 120, 260}) {
        fragments.push_back({read("f/1", a.substr(start, 60)),
                             read("f/2", revcomp(a.substr(start + 60, 60)))});
    }
    const panvar::InternedHaplotypes interned = intern({a, b, c});
    panvar::FragmentScoringOptions options;
    options.mismatch_band = 2;
    options.minimum_piece_bp = 8;
    const panvar::LibraryModel model = library();

    const std::vector<panvar::HaplotypeCandidate> full = {{0, 1}, {1, 2}, {0, 0}, {0, 2}};
    const panvar::FragmentScoringResult all =
        panvar::score_diplotypes(interned, full, fragments, model, options);

    // Gate 8: reordering the request must not change any raw score.
    std::vector<panvar::HaplotypeCandidate> shuffled = full;
    std::reverse(shuffled.begin(), shuffled.end());
    const panvar::FragmentScoringResult reordered =
        panvar::score_diplotypes(interned, shuffled, fragments, model, options);
    for (std::size_t i = 0; i < full.size(); ++i) {
        require_close(all.scores[i].log_likelihood,
                      reordered.scores[full.size() - 1 - i].log_likelihood, 1e-12,
                      "candidate order must not change a raw score");
    }

    // Gate 9: one candidate alone must score identically to the same candidate in a larger list.
    for (std::size_t i = 0; i < full.size(); ++i) {
        const panvar::FragmentScoringResult alone =
            panvar::score_diplotypes(interned, {full[i]}, fragments, model, options);
        require_close(alone.scores[0].log_likelihood, all.scores[i].log_likelihood, 1e-12,
                      "a candidate scored alone must match its score from a larger list");
    }

    // Non-vacuity: the candidates must not all score the same, or the gates above are empty.
    require(std::abs(all.scores[0].log_likelihood - all.scores[2].log_likelihood) > 1e-6,
            "the fixture must separate candidates, or candidate-independence proves nothing");
}

void test_homozygote_counts_two_copies() {
    std::mt19937_64 rng(555);
    const std::string a = random_bases(rng, 380);
    std::vector<panvar::ReadPair> fragments;
    for (std::size_t start : {40, 200}) {
        fragments.push_back({read("f/1", a.substr(start, 55)),
                             read("f/2", revcomp(a.substr(start + 60, 55)))});
    }
    const panvar::InternedHaplotypes interned = intern({a});
    panvar::FragmentScoringOptions options;
    options.mismatch_band = 2;
    options.minimum_piece_bp = 8;
    const panvar::LibraryModel model = library();
    const panvar::FragmentScoringResult result =
        panvar::score_diplotypes(interned, {{0, 0}}, fragments, model, options);

    const double one_copy = panvar::haplotype_exposure(
        a.size(), [&] {
            panvar::DiscreteInsertDistribution insert;
            const double lo = std::ceil(model.insert_mean_bp -
                                        model.support_sigmas * model.insert_sd_bp);
            const double hi = std::floor(model.insert_mean_bp +
                                         model.support_sigmas * model.insert_sd_bp);
            insert.support = {static_cast<std::size_t>(lo), static_cast<std::size_t>(hi)};
            double log_z = kNegInf;
            for (std::size_t l = insert.support.min_bp; l <= insert.support.max_bp; ++l) {
                const double z = (static_cast<double>(l) - model.insert_mean_bp) / model.insert_sd_bp;
                log_z = panvar::log_add_exp(log_z, -0.5 * z * z);
            }
            for (std::size_t l = insert.support.min_bp; l <= insert.support.max_bp; ++l) {
                const double z = (static_cast<double>(l) - model.insert_mean_bp) / model.insert_sd_bp;
                insert.log_probability_by_offset.push_back(-0.5 * z * z - log_z);
            }
            return insert;
        }());
    require_close(result.scores[0].exposure_term,
                  -model.fragment_start_intensity * 2.0 * one_copy, 1e-9,
                  "a homozygote is exposed as two chromosome copies");
    require(one_copy > 0.0, "exposure must be positive, or the gate above is vacuous");
}

void test_refusals() {
    std::mt19937_64 rng(17);
    const std::string haplotype = random_bases(rng, 300);
    const panvar::InternedHaplotypes interned = intern({haplotype});
    const panvar::LibraryModel model = library();
    panvar::FragmentScoringOptions options;
    options.mismatch_band = 2;
    options.minimum_piece_bp = 8;

    // Gate 13: version one is substitutions over ACGT. An ambiguous base is refused by name, not
    // absorbed as a mismatch, which would turn missing data into evidence.
    const std::vector<panvar::ReadPair> ambiguous = {
        {read("n/1", haplotype.substr(10, 50)),
         read("n/2", std::string(20, 'N') + revcomp(haplotype.substr(90, 30)))}};
    require_throws([&] { panvar::score_diplotypes(interned, {{0, 0}}, ambiguous, model, options); },
                   "a non-ACGT read base must refuse with a named reason");

    // Gate 14: a band too wide to index refuses outright rather than returning a partial ranking.
    const std::vector<panvar::ReadPair> ordinary = {
        {read("f/1", haplotype.substr(10, 50)), read("f/2", revcomp(haplotype.substr(90, 50)))}};
    panvar::FragmentScoringOptions unindexable = options;
    unindexable.mismatch_band = 40;
    require_throws(
        [&] { panvar::score_diplotypes(interned, {{0, 0}}, ordinary, model, unindexable); },
        "a band that shreds reads into unindexable pieces must refuse, not degrade silently");

    require_throws(
        [&] { panvar::score_diplotypes(interned, {{0, 5}}, ordinary, model, options); },
        "a candidate naming an unknown haplotype must refuse");
}

void test_interning() {
    const std::string shared = "ACGTACGTGGCCAATTACGTACGTGGCCAATT";
    // Two different allele paths that spell identical bytes are one haplotype to any read-based
    // method, but the caller must still be able to see that the candidate set is narrower than the
    // number of paths suggests.
    const panvar::InternedHaplotypes interned = panvar::intern_haplotypes(
        {panvar::MosaicPath{{0, 1}}, panvar::MosaicPath{{1, 0}}, panvar::MosaicPath{{2, 2}}},
        {shared, shared, shared + "TTTT"});
    require(interned.sequences.size() == 2, "identical bytes must be scored once");
    require(interned.index_of[0] == interned.index_of[1] &&
                interned.index_of[2] != interned.index_of[0],
            "interning must map every input to its unique sequence");
    require(interned.distinct_allele_paths == 3,
            "distinct allele paths must stay visible after byte interning");
    require_throws([&] { panvar::intern_haplotypes({panvar::MosaicPath{{0}}}, {shared, shared}); },
                   "mismatched path and sequence counts must refuse");
}

void test_counters_are_measured() {
    std::mt19937_64 rng(8080);
    const std::string haplotype = random_bases(rng, 600);
    std::vector<panvar::ReadPair> fragments;
    for (std::size_t start : {50, 150, 300, 420}) {
        fragments.push_back({read("f/1", haplotype.substr(start, 55)),
                             read("f/2", revcomp(haplotype.substr(start + 60, 55)))});
    }
    const panvar::InternedHaplotypes interned = intern({haplotype});
    panvar::FragmentScoringOptions options;
    options.mismatch_band = 2;
    options.minimum_piece_bp = 8;
    const panvar::FragmentScoringResult result =
        panvar::score_diplotypes(interned, {{0, 0}}, fragments, library(), options);

    const panvar::FragmentScoringCounters& counters = result.counters;
    require(counters.candidate_bases_scanned >= haplotype.size(),
            "the scan must walk the haplotype it scored");
    require(counters.piece_hits > 0 && counters.proposed_starts > 0,
            "the piece index must actually propose starts");
    // Gate 5: one placement rediscovered through several pieces is verified and counted once.
    require(counters.unique_starts <= counters.proposed_starts,
            "deduplication must not increase the number of starts");
    require(counters.proposed_starts > counters.unique_starts,
            "the fixture must contain a placement found through several pieces, or gate 5 is empty");
    require(counters.full_read_verifications == counters.unique_starts,
            "every unique proposed start is verified exactly once");
    require(counters.accepted_mate_placements >= fragments.size() * 2,
            "both mates of every planted fragment must place");
    require(counters.physical_origins >= counters.statistical_state_classes,
            "compression can only ever reduce the class count below the origin count");
    std::cout << "  counters: scanned=" << counters.candidate_bases_scanned
              << " piece_hits=" << counters.piece_hits
              << " proposed=" << counters.proposed_starts
              << " unique=" << counters.unique_starts
              << " verified=" << counters.full_read_verifications
              << " origins=" << counters.physical_origins
              << " classes=" << counters.statistical_state_classes << '\n';
}

} // namespace

int main() {
    try {
        test_interning();
        test_planted();
        test_repeats();
        test_mixed_read_lengths();
        test_random_sweep();
        test_candidate_independence();
        test_homozygote_counts_two_copies();
        test_refusals();
        test_counters_are_measured();
        std::cout << "genotype fragment scorer: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "genotype fragment scorer: FAIL: " << error.what() << '\n';
        return 1;
    }
}
