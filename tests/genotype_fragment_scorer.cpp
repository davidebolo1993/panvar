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
            panvar::build_library_insert_distribution(model);
        for (std::size_t h = 0; h < interned.sequences.size(); ++h) {
            const auto placements = panvar::enumerate_placements_exhaustive(
                interned.sequences[h], fragments[f], support, model, band);
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

// The consumed-span decomposition changes a likelihood ONLY where the geometry's floor
// -- insert >= the larger consumed span -- can be reached. Above it the mass at an anchor pair
// factorises into the two mates' whole emissions, so the scorer emits one row per anchor instead of
// one per span, and the answer is identical. That shortcut is worth roughly (2*max_gap_bp+1)^2 in
// the pairing loop, which is the difference between a 60-second donor and a 20-minute one.
//
// Both branches must therefore be certified against the oracle, which always componentises: if only
// the componentised branch were exercised the shortcut would be certified by nothing, and if only
// the fast branch were, the decomposition itself would be.
void test_both_span_paths_agree_with_the_oracle() {
    std::mt19937_64 rng(13);
    const std::string haplotype = random_bases(rng, 420);
    const panvar::LibraryModel model = library();   // library support [75, 165]
    require(model.max_gap_bp > 0, "fixture: this gate is empty with gaps closed");
    const panvar::InternedHaplotypes interned = intern({haplotype});
    panvar::FragmentScoringOptions options;
    options.mismatch_band = 2;
    options.minimum_piece_bp = 8;
    const panvar::DiscreteInsertDistribution insert =
        panvar::build_library_insert_distribution(model);

    // 60 bp mates: widest consumed span 63, below the support floor of 75, so the floor is
    // unreachable and the scorer must take the exact fast path.
    const std::vector<panvar::ReadPair> unreachable = {
        {read("n/1", haplotype.substr(40, 60)), read("n/2", revcomp(haplotype.substr(100, 60)))}};
    // 90 bp fully overlapping mates: the fragment's own floor lifts the support to [90, 165] while
    // the widest consumed span is 93, so the floor CAN bind and the classes must be carried.
    const std::vector<panvar::ReadPair> reachable = {
        {read("w/1", haplotype.substr(200, 90)), read("w/2", revcomp(haplotype.substr(200, 90)))}};

    struct Case { const std::vector<panvar::ReadPair>* fragments; std::uint64_t binding; const char* name; };
    const Case cases[] = {{&unreachable, 0, "floor unreachable (fast path)"},
                          {&reachable, 1, "floor reachable (componentised)"}};
    for (const Case& item : cases) {
        const panvar::FragmentScoringResult scored =
            panvar::score_diplotypes(interned, {{0, 0}}, *item.fragments, model, options);
        require(scored.counters.fragments_with_floor_binding == item.binding,
                std::string("fixture: ") + item.name +
                    " must select the branch it names, or this gate certifies the wrong path");
        const panvar::InsertSupport support =
            panvar::fragment_insert_support((*item.fragments)[0], model);
        const auto placements = panvar::enumerate_placements_exhaustive(
            haplotype, (*item.fragments)[0], support, model, options.mismatch_band);
        const double expected =
            panvar::log_placement_mass(placements, (*item.fragments)[0], insert, model);
        require(std::isfinite(expected), "fixture: the fragment must place, or nothing is compared");
        require_close(scored.mass_of(0, 0), expected, 1e-9,
                      std::string("the scorer must equal the always-componentised oracle on the ") +
                          item.name + " branch");
        std::cout << "  " << item.name << ": scorer mass " << scored.mass_of(0, 0)
                  << " = oracle mass " << expected << '\n';
    }
}

// The exact emission cache must be exactly that: a run with it armed and a run without it must
// agree on every score, every mass and every LOGICAL counter, at any thread count. The cache is
// shared across haplotypes and filled concurrently, so this is also where a dangling window
// pointer or a shallow component copy would surface -- both would be silent otherwise.
void test_emission_cache_is_exact() {
    std::mt19937_64 rng(606);
    // Repeated units, so one read form meets identical windows on several haplotypes and the cache
    // is genuinely exercised rather than merely present.
    const std::string unit = random_bases(rng, 180);
    const std::string spacer = random_bases(rng, 60);
    const std::string a = spacer + unit + spacer + unit + spacer;
    std::string b = a;
    for (int edit = 0; edit < 15; ++edit) b[rng() % b.size()] = "ACGT"[rng() % 4];
    std::string c = a;
    c.insert(c.begin() + 250, 'A');          // a one-base indel between panel members
    c.resize(a.size());

    std::vector<panvar::ReadPair> fragments;
    for (std::size_t start : {30, 120, 210, 300, 390}) {
        if (start + 130 > a.size()) continue;
        fragments.push_back({read("f/1", a.substr(start, 60)),
                             read("f/2", revcomp(a.substr(start + 70, 60)))});
    }
    require(!fragments.empty(), "fixture: the sweep must produce fragments");

    const panvar::InternedHaplotypes interned = intern({a, b, c});
    const panvar::LibraryModel model = library();
    std::vector<panvar::HaplotypeCandidate> candidates;
    for (std::size_t x = 0; x < 3; ++x) {
        for (std::size_t y = x; y < 3; ++y) candidates.push_back({x, y});
    }

    panvar::FragmentScoringOptions uncached;
    uncached.mismatch_band = 3;
    uncached.minimum_piece_bp = 8;
    uncached.threads = 1;
    uncached.emission_cache_bytes = 0;
    const panvar::FragmentScoringResult reference =
        panvar::score_diplotypes(interned, candidates, fragments, model, uncached);
    require(reference.emission_cache.requests == 0,
            "a zero budget must leave the cache disarmed, not merely empty");

    std::uint64_t total_hits = 0;
    double worst_mass = 0.0;
    double worst_score = 0.0;
    for (std::size_t threads : {1, 2, 4, 8}) {
        panvar::FragmentScoringOptions cached = uncached;
        cached.threads = threads;
        cached.emission_cache_bytes = 1u << 22;   // 4 MB: ample for this fixture
        const panvar::FragmentScoringResult actual =
            panvar::score_diplotypes(interned, candidates, fragments, model, cached);

        for (std::size_t f = 0; f < fragments.size(); ++f) {
            for (std::size_t h = 0; h < interned.sequences.size(); ++h) {
                const double want = reference.mass_of(f, h);
                const double got = actual.mass_of(f, h);
                if (want == -std::numeric_limits<double>::infinity() ||
                    got == -std::numeric_limits<double>::infinity()) {
                    require(want == got,
                            "a cached run must place exactly the fragments an uncached run places");
                    continue;
                }
                worst_mass = std::max(worst_mass, std::abs(want - got));
            }
        }
        for (std::size_t k = 0; k < candidates.size(); ++k) {
            worst_score = std::max(worst_score, std::abs(reference.scores[k].log_likelihood -
                                                         actual.scores[k].log_likelihood));
            require(reference.scores[k].fragments_with_mass == actual.scores[k].fragments_with_mass,
                    "fragment disposition must not depend on the cache");
        }
        // The LOGICAL counters describe what the model did, not what the machine avoided, so they
        // must be untouched by caching or by the thread count. This is the check that would catch a
        // cache that quietly changed which placements were found.
        const panvar::FragmentScoringCounters& want = reference.counters;
        const panvar::FragmentScoringCounters& got = actual.counters;
        // One require per counter, not a ten-way conjunction. A bundled assertion names only the
        // bundle when it fails, which sends the reader hunting for which member broke -- the same
        // "evidence you cannot read" fault this suite rejects elsewhere.
        const auto same = [&](const char* name, std::uint64_t expected, std::uint64_t actual_value) {
            require(expected == actual_value,
                    std::string("logical counter ") + name + " must be identical with and without "
                    "the cache at threads=" + std::to_string(threads) + ": uncached=" +
                    std::to_string(expected) + " cached=" + std::to_string(actual_value));
        };
        same("full_read_verifications", want.full_read_verifications, got.full_read_verifications);
        same("dp_evaluations", want.dp_evaluations, got.dp_evaluations);
        same("prefilter_rejections", want.prefilter_rejections, got.prefilter_rejections);
        same("accepted_mate_placements", want.accepted_mate_placements, got.accepted_mate_placements);
        same("valid_fr_states", want.valid_fr_states, got.valid_fr_states);
        same("latent_alignment_pairs", want.latent_alignment_pairs, got.latent_alignment_pairs);
        same("physical_origins", want.physical_origins, got.physical_origins);
        same("statistical_state_classes", want.statistical_state_classes,
             got.statistical_state_classes);
        same("proposed_starts", want.proposed_starts, got.proposed_starts);
        same("unique_starts", want.unique_starts, got.unique_starts);
        // ... while the WORK counters must show the cache actually did something.
        require(actual.emission_cache.hits > 0,
                "fixture: the cache must actually hit, or this gate proves only that an unused "
                "cache is harmless");
        require(actual.emission_cache.executed_dps < got.dp_evaluations,
                "executed recursions must fall below the logical count once reuse is captured");
        total_hits += actual.emission_cache.hits;
    }
    require(worst_mass <= 1e-9 && worst_score <= 1e-9,
            "cached and uncached scores must agree within the declared tolerance");
    std::cout << "  emission cache exact at 1/2/4/8 threads: worst mass delta " << worst_mass
              << ", worst score delta " << worst_score << ", " << total_hits << " hits\n";
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
        a.size(), panvar::build_library_insert_distribution(model), 55);
    require_close(result.scores[0].exposure_term,
                  -(1.0 - model.outlier_weight) * model.fragment_start_intensity * 2.0 * one_copy, 1e-9,
                  "a homozygote is exposed as two chromosome copies");
    require(one_copy > 0.0, "exposure must be positive, or the gate above is vacuous");
}

// Exposure averages over the observed physical floors, because a template shorter than a fragment's
// longer mate cannot produce that fragment. With two floor classes in equal number the exposure must
// be the mean of the two analytic counts, not the count at either one.
void test_exposure_uses_read_length_distribution() {
    std::mt19937_64 rng(777);
    const std::string a = random_bases(rng, 400);
    panvar::LibraryModel model = library();
    model.insert_mean_bp = 110.0;
    model.insert_sd_bp = 25.0;
    model.support_sigmas = 3.0;          // support [35, 185]: floors of 40 and 90 both truncate it
    const std::vector<panvar::ReadPair> fragments = {
        {read("s/1", a.substr(10, 40)), read("s/2", revcomp(a.substr(80, 40)))},
        {read("l/1", a.substr(150, 90)), read("l/2", revcomp(a.substr(170, 90)))},
    };
    panvar::FragmentScoringOptions options;
    options.mismatch_band = 2;
    options.minimum_piece_bp = 8;
    const panvar::FragmentScoringResult result =
        panvar::score_diplotypes(intern({a}), {{0, 0}}, fragments, model, options);
    const panvar::DiscreteInsertDistribution insert = panvar::build_library_insert_distribution(model);
    const double short_floor = panvar::haplotype_exposure(a.size(), insert, 40);
    const double long_floor = panvar::haplotype_exposure(a.size(), insert, 90);
    require(short_floor > long_floor + 1.0,
            "fixture: the two floors must give genuinely different exposures");
    require_close(result.exposure[0], 0.5 * (short_floor + long_floor), 1e-9,
                  "exposure must average the analytic count over the observed floor distribution");
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
    // The pairing work admitting indels creates, measured rather than assumed negligible.
    require(counters.latent_alignment_pairs >= counters.valid_fr_states,
            "tested latent alignment pairs cannot be fewer than the accepted ones");
    require(counters.valid_fr_states >= counters.physical_origins,
            "accepted pairs collapse onto origins, so they cannot be the fewer of the two");
    std::cout << "  counters: scanned=" << counters.candidate_bases_scanned
              << " piece_hits=" << counters.piece_hits
              << " proposed=" << counters.proposed_starts
              << " unique=" << counters.unique_starts
              << " verified=" << counters.full_read_verifications
              << " tested_pairs=" << counters.latent_alignment_pairs
              << " accepted_pairs=" << counters.valid_fr_states
              << " origins=" << counters.physical_origins
              << " classes=" << counters.statistical_state_classes << '\n';
}


// -------------------------------------------------------------------------------------------
// Certification under the truncated model.
// -------------------------------------------------------------------------------------------

// The reviewer's decisive gate. A fragment with no in-band placement on A but an in-band placement
// on B must leave A's score and fragment disposition identical whether A is scored alone, alongside
// B, or ranked with B. This failed on the previous implementation, where A's upper bound was -42.2
// scored alone and -186.9 alongside B.
void test_candidate_neutral_disposition() {
    std::mt19937_64 rng(5150);
    const std::string a = random_bases(rng, 600);
    std::string b = a;
    b.replace(200, 200, random_bases(rng, 200));

    std::vector<panvar::ReadPair> fragments;
    for (std::size_t start : {20, 430}) {
        fragments.push_back({read("s/1", a.substr(start, 55)),
                             read("s/2", revcomp(a.substr(start + 65, 55)))});
    }
    fragments.push_back({read("d/1", b.substr(230, 55)), read("d/2", revcomp(b.substr(295, 55)))});

    const panvar::LibraryModel model = library();
    panvar::FragmentScoringOptions options;
    options.mismatch_band = 2;
    options.minimum_piece_bp = 8;

    const panvar::InternedHaplotypes only_a = intern({a});
    const panvar::InternedHaplotypes both = intern({a, b});
    const panvar::FragmentScoringResult alone =
        panvar::score_diplotypes(only_a, {{0, 0}}, fragments, model, options);
    const panvar::FragmentScoringResult together =
        panvar::score_diplotypes(both, {{0, 0}, {1, 1}}, fragments, model, options);
    require(alone.mass_of(2, 0) == kNegInf && together.mass_of(2, 1) != kNegInf,
            "fixture: the discriminating fragment must place on B and not on A");

    const auto same_score = [&](const panvar::DiplotypeScore& lhs, const panvar::DiplotypeScore& rhs,
                                const std::string& where) {
        require(lhs.log_likelihood == rhs.log_likelihood && lhs.fragment_term == rhs.fragment_term &&
                    lhs.exposure_term == rhs.exposure_term,
                "A's score must be bit-identical " + where);
        require(lhs.fragments_with_mass == rhs.fragments_with_mass,
                "A's fragment disposition must be identical " + where);
    };
    same_score(alone.scores[0], together.scores[0], "alone and alongside B");

    const panvar::CertifiedRanking ranked = panvar::certify_diplotype_ranking(
        both, {{0, 0}, {1, 1}}, fragments, model, options, {});
    const panvar::CertifiedRanking solo = panvar::certify_diplotype_ranking(
        only_a, {{0, 0}}, fragments, model, options, {});
    const auto entry_for = [](const panvar::CertifiedRanking& ranking, std::size_t id) {
        for (const panvar::CertifiedCandidate& entry : ranking.candidates) {
            if (entry.candidate_id == id) return entry;
        }
        throw std::runtime_error("candidate missing from ranking");
    };
    const panvar::CertifiedCandidate a_ranked = entry_for(ranked, 0);
    const panvar::CertifiedCandidate a_solo = entry_for(solo, 0);
    require(a_ranked.score == a_solo.score && a_ranked.fragment_term == a_solo.fragment_term &&
                a_ranked.numeric_error_bound == a_solo.numeric_error_bound,
            "A's certified score and error bound must not depend on the rest of the candidate set");
    require(a_ranked.fragments_with_mass == a_solo.fragments_with_mass,
            "A's certified disposition must not depend on the rest of the candidate set");
    // Non-vacuity: B must genuinely be different, or the equalities above prove nothing.
    require(entry_for(ranked, 1).fragments_with_mass != a_ranked.fragments_with_mass,
            "fixture: B must place the discriminating fragment that A cannot");
}

// A winner is emitted only when its score minus its error bound exceeds every rival's score plus
// theirs, and a tie is reported as a set rather than resolved.
void test_winner_requires_proof() {
    std::mt19937_64 rng(646);
    const std::string truth = random_bases(rng, 450);
    std::string decoy = truth;
    for (int edit = 0; edit < 60; ++edit) decoy[rng() % decoy.size()] = "ACGT"[rng() % 4];

    std::vector<panvar::ReadPair> fragments;
    for (std::size_t start = 20; start + 130 < truth.size(); start += 25) {
        fragments.push_back({read("f/1", truth.substr(start, 55)),
                             read("f/2", revcomp(truth.substr(start + 70, 55)))});
    }
    const panvar::InternedHaplotypes interned = intern({truth, decoy});
    const panvar::LibraryModel model = library();
    panvar::FragmentScoringOptions options;
    options.mismatch_band = 3;
    options.minimum_piece_bp = 8;

    const panvar::CertifiedRanking ranking = panvar::certify_diplotype_ranking(
        interned, {{0, 0}, {0, 1}, {1, 1}}, fragments, model, options, {});
    require(ranking.status == panvar::RescoreStatus::Complete,
            "a clearly separated fixture must certify a winner");
    require(ranking.winner_set.size() == 1 && ranking.winner_set.front() == 0,
            "the homozygous truth must win");
    require(ranking.certified_gap > 0.0, "a certified winner must have a positive certified gap");
    for (std::size_t i = 1; i < ranking.candidates.size(); ++i) {
        require(ranking.candidates[0].score - ranking.candidates[0].numeric_error_bound >
                    ranking.candidates[i].score + ranking.candidates[i].numeric_error_bound,
                "the winner's lower edge must clear every rival's upper edge");
    }
    require(ranking.candidates[0].numeric_error_bound > 0.0 &&
                ranking.candidates[0].numeric_error_bound < 1e-6,
            "the numeric error bound must be positive and far below any real separation");

    // Two identical candidates cannot be separated, and must not be.
    const panvar::CertifiedRanking tied = panvar::certify_diplotype_ranking(
        interned, {{0, 0}, {0, 0}}, fragments, model, options, {});
    require(tied.status == panvar::RescoreStatus::RescoreIncomplete,
            "two identical candidates must not resolve to a winner");
    require(tied.winner_set.size() == 2, "an undecided ranking must report its whole set");
    require(tied.certified_gap == 0.0, "an undecided ranking has no certified gap");

    // The exact prior enters the certified score, and only through its own candidate.
    const panvar::CertifiedRanking with_prior = panvar::certify_diplotype_ranking(
        interned, {{0, 0}, {1, 1}}, fragments, model, options, {}, {-5.0, -1.0});
    const panvar::CertifiedRanking without_prior = panvar::certify_diplotype_ranking(
        interned, {{0, 0}, {1, 1}}, fragments, model, options, {});
    for (const auto& entry : with_prior.candidates) {
        for (const auto& other : without_prior.candidates) {
            if (entry.candidate_id != other.candidate_id) continue;
            require_close(entry.score - other.score, entry.candidate_id == 0 ? -5.0 : -1.0, 1e-9,
                          "the exact prior must shift only its own candidate's score");
        }
    }
}

// The candidate-set probability is a normalised weight over the retained set, never a posterior.
void test_candidate_set_probability() {
    std::mt19937_64 rng(2024);
    const std::string truth = random_bases(rng, 400);
    std::vector<std::string> pool = {truth};
    for (int d = 0; d < 4; ++d) {
        std::string decoy = truth;
        for (int edit = 0; edit < 30 + 10 * d; ++edit) decoy[rng() % decoy.size()] = "ACGT"[rng() % 4];
        pool.push_back(decoy);
    }
    std::vector<panvar::ReadPair> fragments;
    for (std::size_t start = 20; start + 130 < truth.size(); start += 30) {
        fragments.push_back({read("f/1", truth.substr(start, 55)),
                             read("f/2", revcomp(truth.substr(start + 70, 55)))});
    }
    std::vector<panvar::HaplotypeCandidate> candidates;
    for (std::size_t h = 0; h < pool.size(); ++h) candidates.push_back({h, h});
    panvar::FragmentScoringOptions options;
    options.mismatch_band = 3;
    options.minimum_piece_bp = 8;
    const panvar::CertifiedRanking ranking = panvar::certify_diplotype_ranking(
        intern(pool), candidates, fragments, library(), options, {});

    double total = 0.0;
    for (const auto& entry : ranking.candidates) total += entry.candidate_set_probability;
    require_close(total, 1.0, 1e-9, "candidate-set probabilities must sum to one over the set");
    require(ranking.winner_set.front() == 0 &&
                ranking.candidates.front().candidate_set_probability > 0.5,
            "the planted truth must win and dominate the retained set");
    require(ranking.candidates.size() == candidates.size(),
            "every requested candidate must be reported");
    std::cout << "  certification: winner=" << ranking.winner_set.front()
              << " certified_gap=" << ranking.certified_gap
              << " p_set=" << ranking.candidates.front().candidate_set_probability << '\n';
}

// A proposal score cannot reach the fragment score, because the scorer never sees one.
void test_proposal_cannot_reach_fragment_score() {
    std::mt19937_64 rng(9001);
    const std::string a = random_bases(rng, 350);
    std::string b = a;
    for (int edit = 0; edit < 30; ++edit) b[rng() % b.size()] = "ACGT"[rng() % 4];
    std::vector<panvar::ReadPair> fragments;
    for (std::size_t start : {25, 160}) {
        fragments.push_back({read("f/1", a.substr(start, 55)),
                             read("f/2", revcomp(a.substr(start + 70, 55)))});
    }
    const panvar::InternedHaplotypes interned = intern({a, b});
    panvar::FragmentScoringOptions options;
    options.mismatch_band = 2;
    options.minimum_piece_bp = 8;
    const panvar::FragmentScoringResult forward = panvar::score_diplotypes(
        interned, {{0, 0}, {0, 1}}, fragments, library(), options);
    const panvar::FragmentScoringResult backward = panvar::score_diplotypes(
        interned, {{0, 1}, {0, 0}}, fragments, library(), options);
    require_close(forward.scores[0].log_likelihood, backward.scores[1].log_likelihood, 1e-12,
                  "a candidate's fragment score must not depend on its proposal rank");
    require_close(forward.scores[1].log_likelihood, backward.scores[0].log_likelihood, 1e-12,
                  "a candidate's fragment score must not depend on its proposal rank");
}


// Every declared resource cap stops scoring with ResourceExhausted and returns nothing -- neither a
// score from score_diplotypes nor a ranking from certification. A ranking assembled from the
// placements found before a cap would favour whichever haplotypes happened to be scanned first.
void test_resource_caps_refuse_without_partial_output() {
    std::mt19937_64 rng(31);
    const std::string a = random_bases(rng, 500);
    std::vector<panvar::ReadPair> fragments;
    for (std::size_t start : {20, 150, 300}) {
        fragments.push_back({read("f/1", a.substr(start, 55)),
                             read("f/2", revcomp(a.substr(start + 65, 55)))});
    }
    const panvar::InternedHaplotypes interned = intern({a});
    const panvar::LibraryModel model = library();
    panvar::FragmentScoringOptions base;
    base.mismatch_band = 2;
    base.minimum_piece_bp = 8;

    const panvar::FragmentScoringResult baseline =
        panvar::score_diplotypes(interned, {{0, 0}}, fragments, model, base);
    require(baseline.counters.piece_hits > 1 && baseline.counters.proposed_starts > 1 &&
                baseline.counters.physical_origins > 1 &&
                baseline.counters.statistical_state_classes > 1,
            "fixture: the uncapped run must exceed a cap of one on every counter");

    const auto expect_exhausted = [&](const panvar::FragmentScoringOptions& options,
                                      const std::string& which) {
        bool thrown = false;
        try {
            panvar::score_diplotypes(interned, {{0, 0}}, fragments, model, options);
        } catch (const panvar::ResourceExhausted&) {
            thrown = true;
        }
        require(thrown, "the " + which + " cap must throw ResourceExhausted from scoring");
        thrown = false;
        try {
            panvar::certify_diplotype_ranking(interned, {{0, 0}}, fragments, model, options, {});
        } catch (const panvar::ResourceExhausted&) {
            thrown = true;
        }
        require(thrown, "the " + which + " cap must stop certification without a ranking");
    };
    panvar::FragmentScoringOptions capped = base;
    capped.max_piece_hits = 1;
    expect_exhausted(capped, "piece-hit");
    capped = base;
    capped.max_proposed_starts_per_haplotype = 1;
    expect_exhausted(capped, "proposed-start");
    capped = base;
    capped.max_physical_origins = 1;
    expect_exhausted(capped, "physical-origin");
    capped = base;
    capped.max_state_classes = 1;
    expect_exhausted(capped, "state-class");
    capped = base;
    capped.max_latent_alignment_pairs = 1;
    expect_exhausted(capped, "latent-alignment-pair");
    capped = base;
    capped.max_mass_cells = 1;
    expect_exhausted(capped, "mass-cell");

    capped = base;
    capped.max_physical_origins = 0;
    require_throws([&] { panvar::score_diplotypes(interned, {{0, 0}}, fragments, model, capped); },
                   "a zero resource cap is a malformed declaration and must refuse");
    capped = base;
    capped.max_latent_alignment_pairs = 0;
    require_throws([&] { panvar::score_diplotypes(interned, {{0, 0}}, fragments, model, capped); },
                   "a zero latent-alignment-pair cap is a malformed declaration and must refuse");
}

// GATE 4. Reads simulated from either candidate must prefer the candidate they came from -- with the
// two candidates separated by a one-base homopolymer insertion, which is precisely the difference a
// substitution-only emission could not see. Under that emission a read spanning the indel lost every
// placement on its OWN source and fell to the background, which is how NA19240's two-edit allele lost
// to a sixteen-edit one.
void test_reads_prefer_their_source_candidate() {
    std::mt19937_64 rng(31337);
    std::string a = random_bases(rng, 520);
    a.replace(200, 8, "TTTTTTTT");          // a homopolymer for the indel to live in
    std::string b = a;
    b.insert(204, "T");                     // one extra base: every downstream base shifts
    b[330] = b[330] == 'A' ? 'C' : 'A';     // and one plain substitution

    const auto reads_from = [&](const std::string& source, const std::string& tag) {
        std::vector<panvar::ReadPair> out;
        for (std::size_t start = 120; start + 125 < source.size() && start < 300; start += 8) {
            out.push_back({read(tag + "/1", source.substr(start, 60)),
                           read(tag + "/2", revcomp(source.substr(start + 60, 60)))});
        }
        return out;
    };

    const panvar::InternedHaplotypes interned = intern({a, b});
    panvar::FragmentScoringOptions options;
    options.mismatch_band = 6;
    options.minimum_piece_bp = 8;
    const panvar::LibraryModel model = library();
    require(model.max_gap_bp > 0, "fixture: this gate is meaningless with gaps closed");
    const std::vector<panvar::HaplotypeCandidate> candidates = {{0, 0}, {1, 1}};

    const std::vector<panvar::ReadPair> a_reads = reads_from(a, "a");
    const std::vector<panvar::ReadPair> b_reads = reads_from(b, "b");
    require(!a_reads.empty() && a_reads.size() == b_reads.size(),
            "fixture: both read sets must be non-empty and comparable");

    const panvar::FragmentScoringResult from_a =
        panvar::score_diplotypes(interned, candidates, a_reads, model, options);
    require(from_a.scores[0].log_likelihood > from_a.scores[1].log_likelihood,
            "reads simulated from A must prefer A");
    const panvar::FragmentScoringResult from_b =
        panvar::score_diplotypes(interned, candidates, b_reads, model, options);
    require(from_b.scores[1].log_likelihood > from_b.scores[0].log_likelihood,
            "reads simulated from B must prefer B");

    // Non-vacuity, and the heart of the matter: the reads that span the indel must actually PLACE on
    // their own source. Under the substitution-only emission they placed on neither.
    std::size_t spanning = 0;
    std::size_t placed = 0;
    for (std::size_t f = 0; f < b_reads.size(); ++f) {
        const std::size_t start = 120 + f * 8;
        if (start > 204 || start + 120 < 205) continue;   // does this fragment straddle the indel?
        ++spanning;
        if (from_b.mass_of(f, 1) != kNegInf) ++placed;
    }
    require(spanning > 0, "fixture: some fragment must span the planted indel");
    require(placed == spanning,
            "every fragment spanning the indel must place on the haplotype it came from; losing "
            "them to the background is the NA19240 defect");
    std::cout << "  generator preference: A-reads favour A by "
              << from_a.scores[0].log_likelihood - from_a.scores[1].log_likelihood
              << " nats, B-reads favour B by "
              << from_b.scores[1].log_likelihood - from_b.scores[0].log_likelihood << " nats; "
              << placed << "/" << spanning << " indel-spanning fragments placed\n";
}

// GATE 6. One physical template read in either library orientation is the same template. Swapping
// which mate is named first must therefore leave the mass exactly unchanged -- the two arms each
// carry half, so exchanging them is a relabelling and nothing else.
void test_orientation_roles_are_symmetric() {
    std::mt19937_64 rng(2718);
    const std::string haplotype = random_bases(rng, 420);
    const panvar::InternedHaplotypes interned = intern({haplotype});
    panvar::FragmentScoringOptions options;
    options.mismatch_band = 3;
    options.minimum_piece_bp = 8;
    const panvar::LibraryModel model = library();

    std::vector<panvar::ReadPair> direct;
    std::vector<panvar::ReadPair> swapped;
    for (std::size_t start : {30, 140, 250}) {
        const std::string left = haplotype.substr(start, 55);
        const std::string right = revcomp(haplotype.substr(start + 65, 55));
        direct.push_back({read("f/1", left), read("f/2", right)});
        // The same physical fragment, with the mate roles exchanged.
        swapped.push_back({read("f/1", right), read("f/2", left)});
    }
    const panvar::FragmentScoringResult a =
        panvar::score_diplotypes(interned, {{0, 0}}, direct, model, options);
    const panvar::FragmentScoringResult b =
        panvar::score_diplotypes(interned, {{0, 0}}, swapped, model, options);
    for (std::size_t f = 0; f < direct.size(); ++f) {
        require(a.mass_of(f, 0) != kNegInf, "fixture: the planted fragments must place");
        require_close(a.mass_of(f, 0), b.mass_of(f, 0), 1e-12,
                      "exchanging the mate roles must not change a fragment's mass: each "
                      "orientation carries exactly half, so the swap is a relabelling");
    }
    require_close(a.scores[0].log_likelihood, b.scores[0].log_likelihood, 1e-12,
                  "nor may it change the candidate score");
}

// GATES 7 and 10. Multiplicity counts PLACES. Alternative gapped alignments of one place are latent
// explanations and must be summed into it, never counted as extra physical copies -- otherwise a
// repeat-free locus acquires phantom copy number the moment gaps are opened.
//
// The pairing work this creates is bounded and measured here rather than assumed negligible: for one
// left endpoint and one insert, the forward mate has at most 2d+1 consumed spans and the reverse mate
// at most 2d+1 that reach the same right endpoint, so accepted geometries per physical origin cannot
// exceed (2d+1)^2 -- 49 at the production max_gap_bp of 3.
void test_multiplicity_counts_places_not_alignments() {
    std::mt19937_64 rng(29);
    const std::string unit = random_bases(rng, 130);
    const std::string spacer = random_bases(rng, 40);
    const std::string two = spacer + unit + spacer + unit + spacer;
    const std::string three = spacer + unit + spacer + unit + spacer + unit + spacer;

    const panvar::InternedHaplotypes interned = intern({two, three});
    panvar::FragmentScoringOptions options;
    options.mismatch_band = 3;
    options.minimum_piece_bp = 8;

    // Production-like geometry: the insert support starts at 75 while the widest consumed span is
    // 58, so the geometry's floor is unreachable and the scorer takes the exact fast path. A real
    // library lives here -- 350 bp inserts against 100 bp reads -- and there is NO pairing expansion
    // in this regime, which is the whole point of the shortcut.
    const std::vector<panvar::ReadPair> wide_insert = {
        {read("r/1", unit.substr(0, 55)), read("r/2", revcomp(unit.substr(65, 55)))}};
    // Tight geometry: a short insert support brings the floor within reach of a gapped span, so the
    // consumed-span classes must be carried separately and the pairing genuinely expands.
    panvar::LibraryModel tight_library = library();
    tight_library.insert_mean_bp = 70.0;
    tight_library.insert_sd_bp = 5.0;
    tight_library.support_sigmas = 3.0;            // support [55, 85] against a widest span of 58
    const std::vector<panvar::ReadPair> tight_insert = {
        {read("t/1", unit.substr(0, 55)), read("t/2", revcomp(unit.substr(15, 55)))}};

    std::size_t expanded = 0;
    std::size_t collapsed = 0;
    for (std::size_t max_gap : {0, 3}) {
        for (int which = 0; which < 2; ++which) {
            panvar::LibraryModel model = which == 0 ? library() : tight_library;
            model.max_gap_bp = max_gap;
            const std::vector<panvar::ReadPair>& fragments =
                which == 0 ? wide_insert : tight_insert;
            const panvar::FragmentScoringResult result =
                panvar::score_diplotypes(interned, {{0, 0}}, fragments, model, options);

            // Multiplicity counts PLACES. Three copies against two must carry log(3/2) more mass,
            // whether or not each place is explained by several latent alignments.
            require_close(result.mass_of(0, 1) - result.mass_of(0, 0), std::log(3.0 / 2.0), 0.02,
                          "placement mass must scale with copy number, not with the number of "
                          "latent alignments per copy");

            const panvar::FragmentScoringCounters& c = result.counters;
            require(c.physical_origins > 0 && c.valid_fr_states > 0,
                    "fixture: the repeat must produce origins");
            const std::uint64_t width = 2 * static_cast<std::uint64_t>(max_gap) + 1;
            require(c.valid_fr_states <= c.physical_origins * width * width,
                    "accepted geometries per physical origin cannot exceed (2*max_gap_bp+1)^2");
            require(c.latent_alignment_pairs >= c.valid_fr_states,
                    "tested pairs cannot be fewer than accepted ones");
            if (max_gap == 0 || c.fragments_with_floor_binding == 0) {
                // Either no gapped span exists at all, or the floor is unreachable and the scorer
                // collapses the spans exactly. One latent alignment per origin, by construction.
                require(c.valid_fr_states == c.physical_origins,
                        "with no reachable insert floor every origin must carry exactly one latent "
                        "alignment: that is what the exact fast path means");
                ++collapsed;
            } else {
                require(c.valid_fr_states > c.physical_origins,
                        "with the floor reachable and gaps open, some origin must be explained by "
                        "several latent alignments, or the componentised branch is untested");
                ++expanded;
            }
            std::cout << "  multiplicity max_gap=" << max_gap
                      << (which == 0 ? " wide-insert " : " tight-insert")
                      << ": floor_binding=" << c.fragments_with_floor_binding
                      << " origins=" << c.physical_origins
                      << " accepted_pairs=" << c.valid_fr_states
                      << " tested_pairs=" << c.latent_alignment_pairs << '\n';
        }
    }
    // Both branches must actually run here. If only one did, the other is measured by nothing --
    // which is how the previous version of this gate came to assert a premise the fast path had
    // already made false.
    require(expanded > 0 && collapsed > 0,
            "both the expanding and the collapsed branch must be exercised by this gate");
}

} // namespace

int main() {
    try {
        test_candidate_neutral_disposition();
        test_interning();
        test_planted();
        test_repeats();
        test_mixed_read_lengths();
        test_random_sweep();
        test_candidate_independence();
        test_reads_prefer_their_source_candidate();
        test_orientation_roles_are_symmetric();
        test_multiplicity_counts_places_not_alignments();
        test_both_span_paths_agree_with_the_oracle();
        test_emission_cache_is_exact();
        test_homozygote_counts_two_copies();
        test_exposure_uses_read_length_distribution();
        test_refusals();
        test_resource_caps_refuse_without_partial_output();
        test_counters_are_measured();
        test_winner_requires_proof();
        test_candidate_set_probability();
        test_proposal_cannot_reach_fragment_score();
        std::cout << "genotype fragment scorer: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "genotype fragment scorer: FAIL: " << error.what() << '\n';
        return 1;
    }
}
