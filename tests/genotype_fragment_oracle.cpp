// Phase 6b. The exhaustive placement enumerator is the reference the accelerated scorer will be
// certified against, so it must first be certified itself — against values computed by hand here,
// not against another implementation of the same idea.

#include "panvar/fragment_likelihood.hpp"
#include "panvar/genotype_model.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <map>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace {

constexpr double kNegInf = -std::numeric_limits<double>::infinity();

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

void require_close(double observed, double expected, double tolerance, const std::string& message) {
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

// The substitution-only fixture library: gaps are CLOSED here deliberately. Every gate below that
// asserts an exact placement count, an exact mismatch value or an exact normalisation identity was
// written against the substitution-only model, and those statements stay true only while the
// displacement band is zero. Gates that mean to open gaps use gap_library() and say so.
panvar::LibraryModel library(double error = 0.01, double outlier = 0.001) {
    panvar::LibraryModel model;
    model.insert_mean_bp = 20.0;
    model.insert_sd_bp = 4.0;
    model.support_sigmas = 3.0;
    model.substitution_error = error;
    model.outlier_weight = outlier;
    model.fragment_start_intensity = 0.05;
    model.max_gap_bp = 0;
    model.provenance = "fixture: declared";
    return model;
}

// -------------------------------------------------------------------------------------------
// Primitives, against hand arithmetic.
// -------------------------------------------------------------------------------------------

void test_emission() {
    // Eight matching bases at 1% error, plus two mismatches each worth error/3.
    require_close(panvar::log_emission(10, 2, 0.01),
                  8.0 * std::log(0.99) + 2.0 * std::log(0.01 / 3.0), 1e-12,
                  "emission must be (1-e)^matches * (e/3)^mismatches");
    require_close(panvar::log_emission(10, 0, 0.01), 10.0 * std::log(0.99), 1e-12,
                  "a perfect read emits (1-e)^L");
    require(panvar::log_emission(10, 0, 0.01) > panvar::log_emission(10, 1, 0.01),
            "an extra mismatch must cost likelihood");
    require_throws([] { return panvar::log_emission(5, 6, 0.01); },
                   "more mismatches than bases must refuse");
    require_throws([] { return panvar::log_emission(5, 1, 0.75); },
                   "an uninformative error rate must refuse");

    const std::string haplotype = "ACGTACGTAC";
    require(panvar::count_mismatches(haplotype, "CGTA", 1).value() == 0,
            "an exact interior match has no mismatch");
    require(panvar::count_mismatches(haplotype, "CATA", 1).value() == 1,
            "one substituted base is one mismatch");
    require(!panvar::count_mismatches(haplotype, "ACGTA", 7).has_value(),
            "a read running off the end has no placement");
    // The band is part of the model: a placement beyond it is omitted, not scored.
    require(!panvar::count_mismatches(haplotype, "CATA", 1, 0).has_value(),
            "a placement outside the declared band must be omitted");
    require(panvar::count_mismatches(haplotype, "CATA", 1, 1).value() == 1,
            "a placement exactly at the band must be kept");
    require_throws([&] { return panvar::count_mismatches(haplotype, "CNTA", 1); },
                   "a non-ACGT read base must refuse, not count as a mismatch");
    require_throws([&] { return panvar::count_mismatches("ACNTACGT", "ACNT", 0); },
                   "a non-ACGT haplotype base must refuse");
}

void test_fr_geometry() {
    panvar::InsertSupport support{10, 30};
    // Forward mate at 5 of length 8; reverse mate ends at 25 -> insert 20.
    require(panvar::fr_insert_length(5, 8, 17, 8, support).value() == 20,
            "insert spans the forward start to the reverse end");
    // Fully overlapping mates: the floor is the longer mate, never the sum.
    require(panvar::fr_insert_length(5, 12, 5, 12, support).value() == 12,
            "overlapping mates have an insert equal to the longer mate, not their sum");
    require(!panvar::fr_insert_length(5, 12, 5, 6, support).has_value(),
            "a span shorter than the longer mate is not a physical state");
    // A shorter fragment must not inherit a longer one's floor.
    require(panvar::fr_insert_length(0, 10, 0, 10, support).value() == 10,
            "a short fragment keeps its own floor");
    require(!panvar::fr_insert_length(20, 8, 5, 8, support).has_value(),
            "a reverse mate ending before the forward start is not an FR pair");
    require(!panvar::fr_insert_length(5, 8, 40, 8, support).has_value(),
            "a span outside the declared support is rejected");
}

void test_exposure() {
    panvar::DiscreteInsertDistribution insert;
    insert.support = {4, 6};
    // Hand distribution: 0.25, 0.5, 0.25 over inserts 4, 5, 6.
    insert.log_probability_by_offset = {std::log(0.25), std::log(0.5), std::log(0.25)};
    // L = 10 gives 7, 6 and 5 start positions for inserts 4, 5 and 6.
    require_close(panvar::haplotype_exposure(10, insert),
                  0.25 * 7.0 + 0.5 * 6.0 + 0.25 * 5.0, 1e-12,
                  "exposure is the insert-weighted count of valid start positions");
    // A haplotype shorter than an insert offers no start for it.
    require_close(panvar::haplotype_exposure(5, insert), 0.25 * 2.0 + 0.5 * 1.0, 1e-12,
                  "an insert longer than the haplotype contributes no exposure");
    require_close(panvar::haplotype_exposure(3, insert), 0.0, 1e-12,
                  "a haplotype shorter than every insert has zero exposure");
    require(panvar::haplotype_exposure(20, insert) > panvar::haplotype_exposure(10, insert),
            "a longer haplotype is exposed to more fragment starts");
}

void test_fragment_contribution() {
    const panvar::LibraryModel model = library();
    const double background = -50.0;
    const double mass = -8.0;
    require_close(
        panvar::log_fragment_contribution(mass, background, model),
        std::log((1.0 - model.outlier_weight) * model.fragment_start_intensity * std::exp(mass) +
                 model.outlier_weight * std::exp(background)),
        1e-9,
        "fragment contribution must be the declared (1-eta)*lambda*mass + eta*background mixture");

    // A fragment that fits no candidate is absorbed by the outlier arm rather than annihilating
    // every candidate together.
    const double absent = panvar::log_fragment_contribution(kNegInf, background, model);
    require_close(absent, std::log(model.outlier_weight) + background, 1e-9,
                  "an unplaceable fragment falls back to the outlier arm exactly");
    // Without an outlier arm that situation is a modelling error, not a score of minus infinity.
    panvar::LibraryModel strict = model;
    strict.outlier_weight = 0.0;
    require_throws([&] { panvar::log_fragment_contribution(kNegInf, background, strict); },
                   "an unplaceable fragment with no outlier mixture must refuse");
}

// -------------------------------------------------------------------------------------------
// The oracle, against hand enumeration.
// -------------------------------------------------------------------------------------------

// A deliberately tiny locus with one exact planted origin.
void test_single_origin() {
    const std::string haplotype = "TTTTACGTACGTGGCCAAAGTCTGACCATTTT";
    const std::string mate1 = "ACGTACGT";               // matches at 4
    const std::string mate2 = revcomp("TGACCATT");      // its revcomp matches at 20, ends at 28
    const panvar::ReadPair fragment{read("f/1", mate1), read("f/2", mate2)};
    const panvar::InsertSupport support{16, 30};

    const auto placements =
        panvar::enumerate_placements_exhaustive(haplotype, fragment, support, library(), 0);
    require(placements.size() == 1, "one planted origin must yield exactly one placement");
    const panvar::FragmentPlacement& only = placements.front();
    require(only.orientation == panvar::MateOrientation::FirstForward,
            "mate one reads off the forward strand here");
    require(only.forward_start == 4 && only.reverse_start == 22,
            "hand-placed mates must be found at their planted coordinates");
    require(only.insert_bp == 26, "insert is the reverse end minus the forward start");
    require(only.forward_mismatches == 0 && only.reverse_mismatches == 0,
            "an exact planted origin has no mismatch");

    panvar::DiscreteInsertDistribution insert;
    insert.support = support;
    insert.log_probability_by_offset.assign(support.max_bp - support.min_bp + 1,
                                            std::log(1.0 / 15.0));
    const panvar::LibraryModel model = library();
    // At band zero a perfect mate's truncated emission is (1-e)^8 / (1-e)^8 = 1, so the mass is the
    // insert probability times the 1/2 orientation prior, and nothing else.
    require_close(
        panvar::log_placement_mass(placements, fragment, insert, model),
        std::log(1.0 / 15.0) + std::log(0.5), 1e-12,
        "single-origin mass is the insert probability times the orientation prior");
    // A wider band still renormalises each mate by the binomial mass it admits -- but the band now
    // acts where the placement is made, not where it is scored, so the claim is asserted on the
    // emission itself. Re-enumerating at band one would also admit other one-mismatch origins and
    // change the mass for a second, unrelated reason, which would test two things at once.
    const double e = model.substitution_error;
    const double admitted =
        std::log(std::pow(1.0 - e, 8) + 8.0 * e * std::pow(1.0 - e, 7));   // F(1, 8)
    const auto exact_at_band_zero =
        panvar::gapped_emission(haplotype, fragment.first.bases, 4, model, 0);
    const auto exact_at_band_one =
        panvar::gapped_emission(haplotype, fragment.first.bases, 4, model, 1);
    require(exact_at_band_zero.has_value() && exact_at_band_one.has_value(),
            "the planted mate must place under both bands");
    require_close(exact_at_band_zero->log_emission, 0.0, 1e-12,
                  "at band zero a perfect mate is the only admitted sequence, so its emission is one");
    require_close(exact_at_band_one->log_emission, 8.0 * std::log(1.0 - e) - admitted, 1e-12,
                  "a wider band renormalises each mate by the binomial mass it admits");
}

// The other library orientation must be populated, not quietly skipped. The physical template spans
// [4,30): its left end reads ACGTACGT off the forward strand, its right end reads AATGGTCA off the
// reverse strand. Which mate got which end is a property of the fragment, not of the locus, so both
// mate orderings must find the same physical origin with the same insert.
void test_both_orientations() {
    const std::string haplotype = "TTTTACGTACGTGGCCAAAGTCTGACCATTTT";
    const std::string left_forward = "ACGTACGT";               // haplotype[4..12)
    const std::string right_reverse = revcomp("TGACCATT");     // revcomp of haplotype[22..30)
    const panvar::InsertSupport support{16, 30};

    const panvar::ReadPair first_reads_left{
        read("f/1", left_forward), read("f/2", right_reverse)};
    const panvar::ReadPair second_reads_left{
        read("f/1", right_reverse), read("f/2", left_forward)};

    const auto arm_a =
        panvar::enumerate_placements_exhaustive(haplotype, first_reads_left, support,
                                                library(), 0);
    const auto arm_b =
        panvar::enumerate_placements_exhaustive(haplotype, second_reads_left, support,
                                                library(), 0);
    require(arm_a.size() == 1 && arm_b.size() == 1,
            "the same physical fragment must be found whichever mate read the left end");
    require(arm_a.front().orientation == panvar::MateOrientation::FirstForward,
            "mate one reading the forward strand is the FirstForward arm");
    require(arm_b.front().orientation == panvar::MateOrientation::SecondForward,
            "mate two reading the forward strand is the SecondForward arm");
    require(arm_a.front().insert_bp == 26 && arm_b.front().insert_bp == 26,
            "one physical template has one insert length whichever mate named its left end");
    require(arm_a.front().forward_start == arm_b.front().forward_start &&
                arm_a.front().reverse_start == arm_b.front().reverse_start,
            "both arms must resolve to the same coordinates on the haplotype");
}

// Repeat-copy multiplicity is evidence. Two identical copies must both survive, with identical
// statistics, because they are two places the fragment could physically have come from.
void test_repeat_multiplicity() {
    const std::string unit = "ACGTACGTAAGGTTCCACCATGGT";
    const std::string haplotype = "TT" + unit + "TT" + unit + "TT";
    const std::string mate1 = "ACGTACGT";
    const std::string mate2 = revcomp("ACCATGGT");
    const panvar::ReadPair fragment{read("f/1", mate1), read("f/2", mate2)};
    const panvar::InsertSupport support{20, 28};

    const auto placements =
        panvar::enumerate_placements_exhaustive(haplotype, fragment, support, library(), 0);
    require(placements.size() == 2, "two identical repeat copies must yield two physical origins");
    require(placements[0].forward_start != placements[1].forward_start,
            "the two origins must be at distinct coordinates");
    require(placements[0].insert_bp == placements[1].insert_bp &&
                placements[0].forward_mismatches == placements[1].forward_mismatches &&
                placements[0].reverse_mismatches == placements[1].reverse_mismatches,
            "identical repeat copies have identical statistics, which is why they may be compressed "
            "only with an exact log-multiplicity");

    panvar::DiscreteInsertDistribution insert;
    insert.support = support;
    insert.log_probability_by_offset.assign(support.max_bp - support.min_bp + 1,
                                            std::log(1.0 / 9.0));
    const panvar::LibraryModel model = library();
    const double both = panvar::log_placement_mass(placements, fragment, insert, model);
    const double one = panvar::log_placement_mass({placements[0]}, fragment, insert, model);
    require_close(both, one + std::log(2.0), 1e-12,
                  "two origins with equal statistics must contribute exactly twice the mass: "
                  "collapsing them to one would discard real evidence");
}

// The band is part of the declared model. A mismatching origin enters the mass only when the band
// admits it, and both the oracle and any accelerated path must use the same band.
void test_band_is_model() {
    const std::string haplotype = "TTTTACGTACGTGGCCAAAGTCTGACCATTTT";
    const std::string mate1 = "ACGTTCGT";            // one substitution against position 4
    const std::string mate2 = revcomp("TGACCATT");
    const panvar::ReadPair fragment{read("f/1", mate1), read("f/2", mate2)};
    const panvar::InsertSupport support{16, 30};

    require(panvar::enumerate_placements_exhaustive(haplotype, fragment, support,
                                                    library(), 0).empty(),
            "a band of zero omits the mismatching origin entirely");
    const auto banded =
        panvar::enumerate_placements_exhaustive(haplotype, fragment, support, library(), 1);
    require(banded.size() == 1 && banded.front().forward_mismatches == 1,
            "widening the band to one admits exactly that origin");

    panvar::DiscreteInsertDistribution insert;
    insert.support = support;
    insert.log_probability_by_offset.assign(15, std::log(1.0 / 15.0));
    const panvar::LibraryModel model = library();
    const double narrow = panvar::log_placement_mass({}, fragment, insert, model);
    const double wide = panvar::log_placement_mass(banded, fragment, insert, model);
    require(narrow == kNegInf && std::isfinite(wide),
            "the band changes the mass, so it is part of the model and not an optimisation");
}

// A homozygote counts two chromosome copies in both mass and exposure.
void test_homozygote_counts_two_copies() {
    const std::string haplotype = "TTTTACGTACGTGGCCAAAGTCTGACCATTTT";
    const panvar::ReadPair fragment{read("f/1", "ACGTACGT"), read("f/2", revcomp("TGACCATT"))};
    const panvar::InsertSupport support{16, 30};
    panvar::DiscreteInsertDistribution insert;
    insert.support = support;
    insert.log_probability_by_offset.assign(15, std::log(1.0 / 15.0));
    const panvar::LibraryModel model = library();

    const auto placements =
        panvar::enumerate_placements_exhaustive(haplotype, fragment, support, model, 0);
    const double one_copy = panvar::log_placement_mass(placements, fragment, insert, model);
    const double two_copies = panvar::log_add_exp(one_copy, one_copy);
    require_close(two_copies, one_copy + std::log(2.0), 1e-12,
                  "a homozygote's mass is twice one chromosome copy's");
    require_close(2.0 * panvar::haplotype_exposure(haplotype.size(), insert),
                  panvar::haplotype_exposure(haplotype.size(), insert) +
                      panvar::haplotype_exposure(haplotype.size(), insert),
                  1e-12,
                  "a homozygote's exposure is twice one chromosome copy's");
}

// Enumeration order must not change the mass, and a placement rediscovered by several routes must
// not be counted twice. The oracle enumerates coordinates directly, so it is the definition of
// "counted once"; this pins that definition down for the accelerated path to match.
void test_order_and_uniqueness() {
    const std::string unit = "ACGTACGTAAGGTTCCACCATGGT";
    const std::string haplotype = "TT" + unit + "TT" + unit + "TT";
    const panvar::ReadPair fragment{read("f/1", "ACGTACGT"), read("f/2", revcomp("ACCATGGT"))};
    const panvar::InsertSupport support{20, 28};
    panvar::DiscreteInsertDistribution insert;
    insert.support = support;
    insert.log_probability_by_offset.assign(9, std::log(1.0 / 9.0));
    const panvar::LibraryModel model = library();

    auto placements = panvar::enumerate_placements_exhaustive(haplotype, fragment, support, model, 0);
    const double forward = panvar::log_placement_mass(placements, fragment, insert, model);
    std::reverse(placements.begin(), placements.end());
    require_close(panvar::log_placement_mass(placements, fragment, insert, model), forward, 1e-12,
                  "placement order must not change the summed mass");

    // An origin is its two template ends plus the orientation. Keying on the reverse START instead
    // would be wrong now: several reverse starts reach the same right endpoint by consuming
    // different numbers of bases, and those are one place, not several.
    std::set<std::tuple<std::size_t, std::size_t, int>> distinct;
    for (const panvar::FragmentPlacement& placement : placements) {
        distinct.emplace(placement.forward_start, placement.insert_bp,
                         static_cast<int>(placement.orientation));
    }
    require(distinct.size() == placements.size(),
            "each physical origin must appear exactly once in the enumeration");
}


// The truncated emission is a distribution: summed over every mate sequence of length L it is one,
// at any band. C(L,m) 3^m sequences carry exactly m mismatches against a fixed reference.
void test_truncated_emission_normalises() {
    for (double e : {0.001, 0.02, 0.2}) {
        for (std::size_t length : {1, 5, 12}) {
            for (std::size_t band = 0; band <= length + 1; ++band) {
                double total = 0.0;
                for (std::size_t m = 0; m <= std::min(band, length); ++m) {
                    const double sequences = std::exp(std::lgamma(length + 1.0) -
                                                      std::lgamma(m + 1.0) -
                                                      std::lgamma(length - m + 1.0)) *
                                             std::pow(3.0, static_cast<double>(m));
                    total += sequences *
                             std::exp(panvar::log_truncated_emission(length, m, e, band));
                }
                require_close(total, 1.0, 1e-9,
                              "the truncated emission must sum to one over every mate sequence");
            }
            require(panvar::log_truncated_emission(length, 2, e, 1) == kNegInf,
                    "a mate beyond the band has zero emission under the truncated model");
        }
    }
}

// The normalisation gate now lives after gap_library(), because it must be run with gaps genuinely
// open across a grid of displacement bands. See test_signal_density_integrates_to_exposure below.


// The complete intensity -- signal plus background -- summed over every observation must equal the
// score's integrated intensity: (1-eta) * lambda * E for the signal, plus eta for the background of this
// single read-length class. The outlier weight is large here so that the pre-review exposure term,
// lambda * E without the (1-eta) factor, is visibly wrong and restoring it fails this gate.
void test_complete_intensity_integrates() {
    const std::string haplotype = "ACGTTGCAAGTCCA";
    panvar::LibraryModel model;
    model.insert_mean_bp = 8.0;
    model.insert_sd_bp = 2.0;
    model.support_sigmas = 2.0;
    model.substitution_error = 0.1;
    model.outlier_weight = 0.25;
    model.fragment_start_intensity = 0.4;
    model.provenance = "fixture: declared";
    const panvar::DiscreteInsertDistribution insert = panvar::build_library_insert_distribution(model);
    const std::size_t mate = 3;
    const std::size_t band = 1;
    const double exposure = panvar::haplotype_exposure(haplotype.size(), insert, mate);

    // One chromosome copy: the signal density of a single haplotype.
    double total = 0.0;
    for (std::size_t code1 = 0; code1 < 64; ++code1) {
        for (std::size_t code2 = 0; code2 < 64; ++code2) {
            std::string r1(mate, 'A');
            std::string r2(mate, 'A');
            for (std::size_t i = 0; i < mate; ++i) {
                r1[i] = "ACGT"[(code1 >> (2 * i)) & 3];
                r2[i] = "ACGT"[(code2 >> (2 * i)) & 3];
            }
            const panvar::ReadPair fragment{read("o/1", r1), read("o/2", r2)};
            const auto placements = panvar::enumerate_placements_exhaustive(
                haplotype, fragment, panvar::fragment_insert_support(fragment, model), model, band);
            const double mass = panvar::log_placement_mass(placements, fragment, insert, model);
            total += std::exp(panvar::log_fragment_contribution(
                mass, panvar::log_background_probability(fragment), model));
        }
    }
    const double signal = panvar::expected_signal_fragments(exposure, model);
    require_close(total, signal + model.outlier_weight, 1e-9,
                  "signal plus background intensity must integrate to (1-eta)*lambda*E + eta");
    const double unnormalised = model.fragment_start_intensity * exposure;
    require(std::abs(unnormalised + model.outlier_weight - total) > 1e-3,
            "fixture: lambda*E without (1-eta) must disagree with the integrated intensity");

    // The scorer's exposure term is exactly the candidate-dependent part of that integral.
    const std::vector<panvar::ReadPair> fragments = {
        {read("f/1", haplotype.substr(0, mate)), read("f/2", revcomp(haplotype.substr(6, mate)))}};
    panvar::FragmentScoringOptions options;
    options.mismatch_band = band;
    options.minimum_piece_bp = 1;
    const panvar::FragmentScoringResult scored = panvar::score_diplotypes(
        panvar::intern_haplotypes({panvar::MosaicPath{{0}}}, {haplotype}), {{0, 0}}, fragments, model, options);
    require_close(scored.scores[0].exposure_term, -2.0 * signal, 1e-12,
                  "a homozygote's exposure term must be minus the integrated signal of two copies");
    std::cout << "  complete intensity: summed " << total << " = (1-eta)*lambda*E + eta "
              << signal + model.outlier_weight << '\n';
}

} // namespace

// A library whose indel parameters are stated here rather than inherited, so a gate cannot silently
// change meaning when a default moves.
panvar::LibraryModel gap_library(double substitution_error, std::size_t max_gap_bp) {
    panvar::LibraryModel model;
    model.insert_mean_bp = 350.0;
    model.insert_sd_bp = 40.0;
    model.support_sigmas = 4.0;
    model.substitution_error = substitution_error;
    model.outlier_weight = 0.001;
    model.fragment_start_intensity = 0.05;
    model.gap_open = 1e-4;
    model.gap_extend = 0.5;
    model.max_gap_bp = max_gap_bp;
    model.provenance = "fixture: declared";
    return model;
}

// GATE 3, as an equality rather than a control: with the displacement band closed the indel-aware
// emission must reproduce the substitution-only truncated emission exactly, including its refusal
// beyond the band. If this drifts, every number the old model produced has silently changed.
void test_gapped_reduces_to_substitution_model() {
    std::string haplotype(24, 'A');
    for (std::size_t i = 0; i < haplotype.size(); ++i) haplotype[i] = "ACGT"[(i * 7 + 3) % 4];
    std::size_t compared = 0;
    for (double e : {0.002, 0.02, 0.2}) {
        const panvar::LibraryModel model = gap_library(e, 0);
        for (std::size_t band : {0, 1, 2, 3}) {
            for (std::size_t length : {5, 8, 12}) {
                for (std::size_t m = 0; m <= std::min(band + 1, length); ++m) {
                    std::string bases = haplotype.substr(1, length);
                    for (std::size_t j = 0; j < m; ++j) {
                        // Distinct positions, and a base that really differs.
                        const char was = bases[j];
                        bases[j] = was == 'A' ? 'C' : 'A';
                    }
                    const auto fit = panvar::gapped_emission(haplotype, bases, 1, model, band);
                    const double expected = panvar::log_truncated_emission(length, m, e, band);
                    if (m > band) {
                        require(!fit.has_value() || fit->log_emission == kNegInf,
                                "beyond the band the gapped emission must refuse, as the "
                                "substitution-only model does");
                        continue;
                    }
                    require(fit.has_value(), "within the band the gapped emission must exist");
                    require_close(fit->log_emission, expected, 1e-9,
                                  "with max_gap_bp 0 the gapped emission must equal the "
                                  "substitution-only truncated emission exactly");
                    require(fit->consumed_bp == length,
                            "an ungapped alignment consumes exactly the read length");
                    ++compared;
                }
            }
        }
    }
    std::cout << "  gapped emission reduces to the substitution model on " << compared
              << " (error, band, length, mismatch) combinations\n";
}

// The property the normaliser exists for, and the precondition for the exposure identity: at a fixed
// start the emission is a distribution over the reads of that length, so summing it over every one of
// them is exactly one. Checked with gaps closed AND open, because an unnormalised gap model would
// still pass the closed case.
void test_gapped_emission_normalises() {
    const std::string haplotype = "ACGTTGCAAGTCCATGA";
    const std::size_t length = 4;
    for (std::size_t max_gap : {0, 1, 2}) {
        for (std::size_t band : {1, 2, 4}) {
            const panvar::LibraryModel model = gap_library(0.1, max_gap);
            double total = 0.0;
            std::size_t reads = 0;
            for (std::size_t code = 0; code < 256; ++code) {
                std::string bases(length, 'A');
                for (std::size_t i = 0; i < length; ++i) bases[i] = "ACGT"[(code >> (2 * i)) & 3];
                const auto fit = panvar::gapped_emission(haplotype, bases, 4, model, band);
                if (fit.has_value() && fit->log_emission != kNegInf) total += std::exp(fit->log_emission);
                ++reads;
            }
            require(reads == 256, "the gate must sum over every read of this length");
            require_close(total, 1.0, 1e-9,
                          "the gapped emission must sum to one over every read of the observed "
                          "length, or it is not a probability and the exposure identity fails");
        }
    }
    std::cout << "  gapped emission normalises over all 256 reads at every band and gap width\n";
}

// GATE 1. The NA19240 mechanism in miniature: a read carrying one extra base in a homopolymer. Under
// a fixed-position emission every base after the run is shifted, the mismatches leave the band, and
// the read gets NO placement on the haplotype it actually came from -- which is how a two-edit allele
// lost to a sixteen-edit one. The fixture asserts that failure rather than assuming it, so the gate
// cannot quietly stop testing anything.
void test_one_base_homopolymer_insertion() {
    const std::string prefix = "CGTACGT";
    const std::string tail = "TCTCTGCAAAAC";
    const std::string haplotype = prefix + std::string(8, 'A') + tail + "GGTTCA";
    const std::string bases = prefix + std::string(9, 'A') + tail;   // one extra A
    const std::size_t band = 6;

    const auto ungapped = panvar::count_mismatches(haplotype, bases, 0, band);
    require(!ungapped.has_value(),
            "fixture: the substitution-only model must FAIL on this read, or the gate proves nothing");

    const panvar::LibraryModel model = gap_library(0.002, 2);
    const auto fit = panvar::gapped_emission(haplotype, bases, 0, model, band);
    require(fit.has_value() && std::isfinite(fit->log_emission),
            "a one-base homopolymer insertion must receive a valid gapped placement");
    require(fit->edits == 1, "it must cost exactly one edit, not a shifted run of substitutions");
    require(fit->consumed_bp == bases.size() - 1,
            "an inserted read base consumes one fewer haplotype base");
    std::cout << "  one-base homopolymer insertion: substitution model refuses, gapped emission "
              << fit->log_emission << " at 1 edit\n";
}

// GATE 2. The same, with the read missing a base the haplotype carries.
void test_one_base_deletion() {
    const std::string prefix = "CGTACGT";
    const std::string tail = "TCTCTGCAAAAC";
    const std::string haplotype = prefix + std::string(9, 'A') + tail + "GGTTCA";
    const std::string bases = prefix + std::string(8, 'A') + tail;   // one A short
    const std::size_t band = 6;

    const auto ungapped = panvar::count_mismatches(haplotype, bases, 0, band);
    require(!ungapped.has_value(),
            "fixture: the substitution-only model must FAIL on this read, or the gate proves nothing");

    const panvar::LibraryModel model = gap_library(0.002, 2);
    const auto fit = panvar::gapped_emission(haplotype, bases, 0, model, band);
    require(fit.has_value() && std::isfinite(fit->log_emission),
            "a one-base deletion must receive a valid gapped placement");
    require(fit->edits == 1, "it must cost exactly one edit");
    require(fit->consumed_bp == bases.size() + 1,
            "a deleted read base consumes one extra haplotype base");
    std::cout << "  one-base deletion: substitution model refuses, gapped emission "
              << fit->log_emission << " at 1 edit\n";
}

// -------------------------------------------------------------------------------------------
// Consumed-span components: the correction, and the defect it replaces.
//
// gapped_emission sums alignment mass over every terminal consumed span. Each span ends the mate at
// a different haplotype coordinate, so each implies a different fragment geometry and a different
// insert probability. The previous implementation summed the mass first and attached the total to
// the single best span, weighting every other span's mass with an insert probability that was not
// its own. The gates below pin down the corrected behaviour AND reject the old one.
// -------------------------------------------------------------------------------------------

// A second, genuinely different insert distribution over the same support: strongly skewed toward
// long inserts, where the library's own distribution is near-symmetric. Misattributing mass between
// neighbouring consumed spans costs something quite different under the two, so an identity that
// holds under both is not holding by an accident of shape.
panvar::DiscreteInsertDistribution skewed_insert(const panvar::InsertSupport& support) {
    panvar::DiscreteInsertDistribution out;
    out.support = support;
    const std::size_t span = support.max_bp - support.min_bp + 1;
    std::vector<double> weight(span);
    double total = 0.0;
    for (std::size_t i = 0; i < span; ++i) {
        const double x = static_cast<double>(i + 1);
        weight[i] = x * x * x;
        total += weight[i];
    }
    for (std::size_t i = 0; i < span; ++i) {
        out.log_probability_by_offset.push_back(std::log(weight[i] / total));
    }
    return out;
}

// A near-symmetric discrete distribution over an EXPLICIT support, normalised over it.
//
// The gate builds its own rather than taking the library's, because the exposure identity constrains
// where the support's floor may sit. A template of length i can host a mate whose alignment consumes
// c haplotype bases only if c <= i, and consumed spans reach read_length + max_gap_bp -- while
// haplotype_exposure counts EVERY template at or above its floor. If the floor were below
// max_mate_length + max_gap_bp, templates between the two would be counted by the exposure but could
// not host the widest admissible alignments, and the identity would fail for a reason that has
// nothing to do with normalisation. Every support below is therefore floored at that value. The
// production library satisfies it with room to spare: inserts near 350 bp against 100 + 3.
panvar::DiscreteInsertDistribution gaussian_insert(const panvar::InsertSupport& support,
                                                   double mean, double sd) {
    panvar::DiscreteInsertDistribution out;
    out.support = support;
    std::vector<double> weight;
    double total = 0.0;
    for (std::size_t i = support.min_bp; i <= support.max_bp; ++i) {
        const double z = (static_cast<double>(i) - mean) / sd;
        weight.push_back(std::exp(-0.5 * z * z));
        total += weight.back();
    }
    for (double w : weight) out.log_probability_by_offset.push_back(std::log(w / total));
    return out;
}

// The support this fixture must use for a given fragment shape, by the rule above.
panvar::InsertSupport safe_support(std::size_t first_len, std::size_t second_len,
                                   std::size_t max_gap_bp, std::size_t haplotype_length) {
    panvar::InsertSupport support;
    support.min_bp = std::max(first_len, second_len) + max_gap_bp;
    support.max_bp = haplotype_length;
    return support;
}

// The fixture the whole grid runs on: short enough to enumerate every observation exhaustively, long
// enough that interior starts exist alongside the two boundaries.
panvar::LibraryModel exposure_model(std::size_t max_gap_bp) {
    panvar::LibraryModel model;
    model.insert_mean_bp = 6.0;
    model.insert_sd_bp = 1.5;
    model.support_sigmas = 2.0;
    model.substitution_error = 0.1;
    model.outlier_weight = 0.001;
    model.fragment_start_intensity = 0.05;
    // Deliberately loose, so gapped alignments carry mass worth measuring rather than 1e-12 of it:
    // a collapse defect hidden below the tolerance would be a gate that cannot fail.
    model.gap_open = 5e-3;
    model.gap_extend = 0.4;
    model.max_gap_bp = max_gap_bp;
    model.provenance = "fixture: declared";
    return model;
}

// THE MUTATION, kept here in the suite rather than in the library: exactly what the implementation
// did before the consumed-span correction. gapped_emission sums mass over every terminal span, and
// this attaches that whole sum to the single BEST span's geometry -- so every other span's mass is
// weighted with an insert probability that is not its own. The normalisation gate must reject it.
std::vector<panvar::FragmentPlacement> collapse_to_best_span(
    const std::string& haplotype,
    const panvar::ReadPair& fragment,
    const panvar::InsertSupport& support,
    const panvar::LibraryModel& model,
    std::size_t band) {
    std::vector<panvar::FragmentPlacement> out;
    const std::string rc_first = revcomp(fragment.first.bases);
    const std::string rc_second = revcomp(fragment.second.bases);
    for (int arm = 0; arm < 2; ++arm) {
        const std::string& forward = arm == 0 ? fragment.first.bases : fragment.second.bases;
        const std::string& reverse = arm == 0 ? rc_second : rc_first;
        for (std::size_t fs = 0; fs < haplotype.size(); ++fs) {
            const auto forward_fit = panvar::gapped_emission(haplotype, forward, fs, model, band);
            if (!forward_fit.has_value() || forward_fit->log_emission == kNegInf) continue;
            for (std::size_t re = 0; re <= haplotype.size(); ++re) {
                // Right-anchored, exactly as the corrected implementation anchors it. The ONLY thing
                // this function mutates is the best-span collapse; if it also left-anchored the
                // reverse mate it would be rejected for two reasons at once and would prove nothing
                // about the collapse.
                const auto reverse_fit =
                    panvar::gapped_emission_ending_at(haplotype, reverse, re, model, band);
                if (!reverse_fit.has_value() || reverse_fit->log_emission == kNegInf) continue;
                const auto length = panvar::fr_insert_from_endpoints(
                    fs, forward_fit->consumed_bp, re, reverse_fit->consumed_bp, support);
                if (!length.has_value()) continue;
                panvar::FragmentPlacement placement;
                placement.orientation = arm == 0 ? panvar::MateOrientation::FirstForward
                                                 : panvar::MateOrientation::SecondForward;
                placement.forward_start = fs;
                placement.insert_bp = *length;
                placement.reverse_start = re - reverse_fit->consumed_bp;
                placement.forward_consumed_bp = forward_fit->consumed_bp;
                placement.reverse_consumed_bp = reverse_fit->consumed_bp;
                placement.forward_mismatches = forward_fit->edits;
                placement.reverse_mismatches = reverse_fit->edits;
                placement.latent_alignments = 1;
                // The defect in one line: each mate's WHOLE emission, at the best span's geometry.
                placement.log_latent_mass = forward_fit->log_emission + reverse_fit->log_emission;
                out.push_back(placement);
            }
        }
    }
    return out;
}

struct IntegratedSignal {
    double total = 0.0;
    double first_forward = 0.0;
    std::size_t observations = 0;
    std::size_t origins = 0;
    std::size_t multi_alignment_origins = 0;
    std::size_t boundary_left = 0;    // origins anchored at the very first base
    std::size_t boundary_right = 0;   // origins whose template ends at the very last base
    std::size_t interior = 0;
};

// Sums the signal mass over EVERY observation with the given mate lengths. `collapse` selects the
// MUTATION: when set, each mate's whole emission is attached to its single best consumed span, which
// is exactly what the implementation did before this correction.
IntegratedSignal integrate_signal(
    const std::string& haplotype,
    const panvar::LibraryModel& model,
    const panvar::InsertSupport& support,
    const panvar::DiscreteInsertDistribution& insert,
    std::size_t band,
    std::size_t first_len,
    std::size_t second_len,
    bool collapse) {
    IntegratedSignal out;
    const std::size_t first_codes = std::size_t{1} << (2 * first_len);
    const std::size_t second_codes = std::size_t{1} << (2 * second_len);
    for (std::size_t code1 = 0; code1 < first_codes; ++code1) {
        for (std::size_t code2 = 0; code2 < second_codes; ++code2) {
            std::string r1(first_len, 'A');
            std::string r2(second_len, 'A');
            for (std::size_t i = 0; i < first_len; ++i) r1[i] = "ACGT"[(code1 >> (2 * i)) & 3];
            for (std::size_t i = 0; i < second_len; ++i) r2[i] = "ACGT"[(code2 >> (2 * i)) & 3];
            const panvar::ReadPair fragment{read("o/1", r1), read("o/2", r2)};

            std::vector<panvar::FragmentPlacement> placements;
            if (collapse) {
                placements = collapse_to_best_span(haplotype, fragment, support, model, band);
            } else {
                placements = panvar::enumerate_placements_exhaustive(haplotype, fragment, support,
                                                                     model, band);
            }
            const double mass = panvar::log_placement_mass(placements, fragment, insert, model);
            if (mass != kNegInf) out.total += std::exp(mass);

            std::vector<panvar::FragmentPlacement> arm;
            for (const panvar::FragmentPlacement& p : placements) {
                if (p.orientation == panvar::MateOrientation::FirstForward) arm.push_back(p);
                ++out.origins;
                if (p.latent_alignments > 1) ++out.multi_alignment_origins;
                const bool left = p.forward_start == 0;
                const bool right = p.forward_start + p.insert_bp == haplotype.size();
                if (left) ++out.boundary_left;
                if (right) ++out.boundary_right;
                if (!left && !right) ++out.interior;
            }
            const double arm_mass = panvar::log_placement_mass(arm, fragment, insert, model);
            if (arm_mass != kNegInf) out.first_forward += std::exp(arm_mass);
            ++out.observations;
        }
    }
    return out;
}

// GATES 6 and 8. Summing the signal mass over EVERY possible observation must reproduce the exposure
// exactly, and each library orientation must carry exactly half of it. This is the identity that ties
// the -lambda*E term to the fragment density: if it fails, the two describe different generated
// processes and every likelihood ratio between haplotypes of different lengths is off.
//
// Run over band < max_gap_bp, band == max_gap_bp and band > max_gap_bp, with NO precondition
// relating the two. An earlier reading of this file's failures proposed requiring band > max_gap_bp
// as a declared model precondition; that was a symptom being mistaken for a contract. A correctly
// normalised truncated model owes the identity at every combination -- gap shapes the edit budget
// cannot pay for simply carry zero mass, which costs the identity nothing.
struct GridCell {
    std::size_t max_gap;
    std::size_t band;
    std::size_t first_len;
    std::size_t second_len;
    bool skewed;
};

void test_signal_density_integrates_to_exposure() {
    const std::string haplotype = "ACGTTGCAAG";   // L = 10
    // Every cell keeps max_gap_bp < the shorter mate. At max_gap_bp >= read length the consumption
    // floor read_length - min(read_length, max_gap_bp) collapses to zero and a read may be absorbed
    // entirely into insertions, consuming NO haplotype at all -- which is not a placement in any
    // useful sense and is a regime production never enters (100 bp reads against max_gap_bp 3).
    // Testing the identity there would be testing a degenerate model, not this one.
    const GridCell grid[] = {
        // gaps closed: the substitution-only control, which must stay exact.
        {0, 1, 2, 2, false}, {0, 3, 3, 3, false},
        // band == max_gap, band > max_gap, a skewed insert distribution, and unequal mate lengths.
        {1, 1, 2, 2, false}, {1, 2, 2, 2, false}, {1, 3, 2, 2, true}, {1, 1, 2, 3, false},
        {1, 2, 2, 3, true},
        // band < max_gap, band == max_gap, band > max_gap at a wider displacement band.
        {2, 1, 3, 3, false}, {2, 2, 3, 3, false}, {2, 3, 3, 3, true},
    };
    std::size_t cells = 0;
    std::size_t below = 0, equal = 0, above = 0;
    std::size_t multi_alignment_origins = 0;
    std::size_t boundary_left = 0;
    std::size_t boundary_right = 0;
    std::size_t interior = 0;

    for (const GridCell& cell : grid) {
        const panvar::LibraryModel model = exposure_model(cell.max_gap);
        const panvar::InsertSupport support =
            safe_support(cell.first_len, cell.second_len, cell.max_gap, haplotype.size());
        require(support.min_bp < support.max_bp, "fixture: the support must be non-degenerate");
        const panvar::DiscreteInsertDistribution insert =
            cell.skewed ? skewed_insert(support)
                        : gaussian_insert(support, 0.5 * (support.min_bp + support.max_bp), 2.0);
        const double exposure =
            panvar::haplotype_exposure(haplotype.size(), insert, support.min_bp);
        require(exposure > 0.0, "fixture: exposure must be positive, or the identity is vacuous");

        const IntegratedSignal signal = integrate_signal(
            haplotype, model, support, insert, cell.band, cell.first_len, cell.second_len, false);
        const std::string where =
            "max_gap=" + std::to_string(cell.max_gap) + " band=" + std::to_string(cell.band) +
            " mates=" + std::to_string(cell.first_len) + "/" + std::to_string(cell.second_len) +
            (cell.skewed ? " skewed-insert" : " gaussian-insert");
        require_close(signal.total, exposure, 1e-9,
                      "the signal density summed over every observation must equal the exposure ["
                      + where + "]");
        require_close(signal.first_forward, 0.5 * exposure, 1e-9,
                      "each library orientation must carry exactly half of the exposure [" + where
                      + "]");
        if (cell.max_gap > 0) {
            if (cell.band < cell.max_gap) ++below;
            else if (cell.band == cell.max_gap) ++equal;
            else ++above;
        }
        multi_alignment_origins += signal.multi_alignment_origins;
        boundary_left += signal.boundary_left;
        boundary_right += signal.boundary_right;
        interior += signal.interior;
        ++cells;
    }
    require(cells == sizeof(grid) / sizeof(grid[0]), "the grid must run every declared cell");
    // The identity is required at every relation between the band and the displacement width, with
    // no precondition connecting them. An earlier reading of this file proposed band > max_gap_bp as
    // a declared precondition; it was a symptom mistaken for a contract, and these three counters
    // exist so that reading cannot quietly come back.
    require(below > 0 && equal > 0 && above > 0,
            "the grid must cover band < max_gap, band == max_gap and band > max_gap");
    // Non-vacuity, all four of them. Without multi-alignment origins the componentisation is never
    // exercised; without both boundaries and the interior the identity is only being checked where
    // the haplotype does not constrain the alignment.
    require(multi_alignment_origins > 0,
            "the grid must produce origins explained by several latent alignments, or the "
            "componentisation is untested");
    require(boundary_left > 0 && boundary_right > 0 && interior > 0,
            "the grid must cover boundary and interior starts, or the identity is only checked "
            "where the haplotype end does not bind");
    std::cout << "  exposure identity holds on " << cells
              << " cells (" << below << " band<max_gap, " << equal << " band==max_gap, " << above
              << " band>max_gap); " << multi_alignment_origins
              << " multi-alignment origins, boundaries " << boundary_left << "/" << boundary_right
              << ", interior " << interior << '\n';
}

// GATE 7 (the mutation). Collapsing every consumed span's mass onto the best span is the defect this
// correction removes, so the normalisation gate must REJECT it. With gaps closed the mutation is a
// no-op and the identity must still hold -- which is what proves the gate is detecting the collapse
// rather than merely detecting that something changed.
// How many (forward span, reverse span) pairs the insert floor excludes at admissible anchor pairs.
// Zero exactly when support.min_bp >= max mate length + max_gap_bp, which is what makes the collapse
// inert there; positive is what makes the mutation detectable.
std::size_t floor_exclusions(const std::string& haplotype, const panvar::ReadPair& fragment,
                             const panvar::InsertSupport& support,
                             const panvar::LibraryModel& model, std::size_t band) {
    std::size_t excluded = 0;
    const std::string reverse = revcomp(fragment.second.bases);
    for (std::size_t fs = 0; fs < haplotype.size(); ++fs) {
        const auto forward_fit =
            panvar::gapped_emission(haplotype, fragment.first.bases, fs, model, band);
        if (!forward_fit.has_value()) continue;
        for (std::size_t re = fs + 1; re <= haplotype.size(); ++re) {
            const auto reverse_fit =
                panvar::gapped_emission_ending_at(haplotype, reverse, re, model, band);
            if (!reverse_fit.has_value()) continue;
            const std::size_t insert = re - fs;
            if (insert < support.min_bp || insert > support.max_bp) continue;
            for (const auto& a : forward_fit->components) {
                for (const auto& b : reverse_fit->components) {
                    if (insert < std::max(a.consumed_bp, b.consumed_bp)) ++excluded;
                }
            }
        }
    }
    return excluded;
}

void test_best_span_collapse_is_rejected() {
    const std::string haplotype = "ACGTTGCAAGTC";   // L = 12
    const panvar::ReadPair probe_fragment{read("p/1", haplotype.substr(1, 3)),
                                          read("p/2", revcomp(haplotype.substr(6, 3)))};

    // (a) Gaps closed. Every mate then has exactly one consumed span, so the collapse is a no-op
    // and the identity must still hold. This is what shows the gate is detecting the COLLAPSE and
    // not merely detecting that something changed.
    {
        const panvar::LibraryModel closed = exposure_model(0);
        const panvar::InsertSupport support = safe_support(2, 2, 0, haplotype.size());
        const panvar::DiscreteInsertDistribution insert =
            gaussian_insert(support, 0.5 * (support.min_bp + support.max_bp), 2.0);
        const double exposure =
            panvar::haplotype_exposure(haplotype.size(), insert, support.min_bp);
        const IntegratedSignal collapsed =
            integrate_signal(haplotype, closed, support, insert, 2, 2, 2, true);
        require_close(collapsed.total, exposure, 1e-9,
                      "with gaps closed every mate has one consumed span, so the collapse is a "
                      "no-op and the identity must still hold");
    }

    // (b) Gaps open, insert floor NOT binding. Here the collapse is provably inert, and saying so
    // is part of the contract rather than a hole in it. With the reverse mate anchored at its right
    // endpoint the insert is re - fs: a function of the two ANCHORS alone. The consumed spans enter
    // only through the floor insert >= max(consumed spans), so where that floor cannot bind, the
    // correct mass at an anchor pair factorises as (sum over forward spans) * (sum over reverse
    // spans) -- which is exactly what the collapse computes. The exposure identity therefore CANNOT
    // reject the collapse in this regime, and a gate claiming it does would be asserting something
    // false. Measured: identical to the last bit, with zero span pairs excluded.
    const panvar::LibraryModel model = exposure_model(2);
    {
        const panvar::InsertSupport support = safe_support(3, 3, 2, haplotype.size());
        const panvar::DiscreteInsertDistribution insert =
            gaussian_insert(support, 0.5 * (support.min_bp + support.max_bp), 2.0);
        const double exposure =
            panvar::haplotype_exposure(haplotype.size(), insert, support.min_bp);
        require(floor_exclusions(haplotype, probe_fragment, support, model, 3) == 0,
                "fixture: at or above the safe floor nothing may be excluded by the insert floor");
        for (std::size_t band : {1, 2, 3}) {
            const IntegratedSignal correct =
                integrate_signal(haplotype, model, support, insert, band, 3, 3, false);
            const IntegratedSignal collapsed =
                integrate_signal(haplotype, model, support, insert, band, 3, 3, true);
            require_close(correct.total, exposure, 1e-9,
                          "the componentised path must satisfy the identity here");
            require_close(collapsed.total, correct.total, 1e-12,
                          "where the insert floor cannot bind, the best-span collapse is inert by "
                          "construction and must agree exactly");
        }
    }

    // (c) Gaps open, insert floor BINDING. This is where the componentisation is observable and
    // where the mutation must be rejected: admission is decided per consumed SPAN, and the collapse
    // decides it once from the best span, so it keeps mass the model excludes and drops mass the
    // model keeps.
    std::size_t rejected = 0;
    {
        const panvar::InsertSupport support{3, haplotype.size()};   // below 3 + max_gap_bp = 5
        const panvar::DiscreteInsertDistribution insert =
            gaussian_insert(support, 0.5 * (support.min_bp + support.max_bp), 2.0);
        const std::size_t excluded = floor_exclusions(haplotype, probe_fragment, support, model, 3);
        require(excluded > 0,
                "fixture: below the safe floor the insert floor must exclude span pairs, or there "
                "is nothing for the collapse to get wrong");
        for (std::size_t band : {2, 3}) {
            const IntegratedSignal correct =
                integrate_signal(haplotype, model, support, insert, band, 3, 3, false);
            const IntegratedSignal collapsed =
                integrate_signal(haplotype, model, support, insert, band, 3, 3, true);
            require(std::abs(collapsed.total - correct.total) > 1e-9,
                    "with the insert floor binding, collapsing consumed-span mass onto the best "
                    "span must give a different answer at band=" + std::to_string(band) +
                    "; a gate that both versions pass is not a gate");
            ++rejected;
        }
        std::cout << "  best-span collapse: inert where the insert floor cannot bind (exactly, by "
                     "factorisation), rejected on " << rejected << " floor-binding cells with "
                  << excluded << " span pairs excluded\n";
    }
    require(rejected == 2, "the mutation must be rejected on every floor-binding cell");
}

// GATE 9. Two consumed spans of one mate, carrying non-zero mass, landing on inserts with DIFFERENT
// probabilities. The joint likelihood must equal the component-wise reference -- which applies each
// span's own insert probability -- and must differ from the best-span collapse.
void test_component_factorisation_and_collapse() {
    const std::string haplotype = "ACGTTGCAAGTCAG";
    const panvar::LibraryModel model = exposure_model(2);
    // Deliberately BELOW the safe floor of 3 + max_gap_bp = 5, because that is the only regime in
    // which the componentisation is observable at all.
    //
    // Requirement 9 as originally framed asked for two consumed spans of one mate carrying
    // DIFFERENT insert probabilities. That cannot arise once the reverse mate is anchored at its
    // right endpoint: the insert is re - fs, so every consumed span at one anchor pair shares one
    // insert, and spans differing in insert are different origins. The premise belonged to the
    // left-anchored model this batch removed. What survives of the requirement, and is checked
    // here, is the part that still has content: each span is ADMITTED on its own geometry, and
    // grouping latent alignments onto physical origins does not change the likelihood.
    const panvar::InsertSupport support{3, haplotype.size()};
    // Skewed hard, so mass attributed to the wrong geometry cannot hide inside the tolerance.
    const panvar::DiscreteInsertDistribution insert = skewed_insert(support);
    const std::size_t band = 3;

    const panvar::ReadPair fragment{read("t/1", haplotype.substr(2, 3)),
                                    read("t/2", revcomp(haplotype.substr(6, 3)))};
    require(floor_exclusions(haplotype, fragment, support, model, band) > 0,
            "fixture: the insert floor must exclude some span pairs, or admission per span is "
            "indistinguishable from admission per best span");

    // An independent reference: every (start, span) pair summed directly, with no grouping at all.
    // Summation is associative, so the production path -- which groups latent alignments onto
    // physical origins first -- must agree with it exactly. That equality is the whole of the claim
    // that collapsing by endpoint does not change the likelihood.
    double reference = 0.0;
    std::set<std::size_t> inserts_used;
    std::size_t span_pairs = 0;
    const std::string rc_first = revcomp(fragment.first.bases);
    const std::string rc_second = revcomp(fragment.second.bases);
    for (int arm = 0; arm < 2; ++arm) {
        const std::string& forward = arm == 0 ? fragment.first.bases : fragment.second.bases;
        const std::string& reverse = arm == 0 ? rc_second : rc_first;
        for (std::size_t fs = 0; fs < haplotype.size(); ++fs) {
            const auto forward_fit = panvar::gapped_emission(haplotype, forward, fs, model, band);
            if (!forward_fit.has_value()) continue;
            for (std::size_t re = 0; re <= haplotype.size(); ++re) {
                const auto reverse_fit =
                    panvar::gapped_emission_ending_at(haplotype, reverse, re, model, band);
                if (!reverse_fit.has_value()) continue;
                for (const auto& fspan : forward_fit->components) {
                    for (const auto& rspan : reverse_fit->components) {
                        const auto length = panvar::fr_insert_from_endpoints(
                            fs, fspan.consumed_bp, re, rspan.consumed_bp, support);
                        if (!length.has_value()) continue;
                        inserts_used.insert(*length);
                        ++span_pairs;
                        reference += std::exp(std::log(0.5) + insert.log_probability(*length) +
                                              fspan.log_emission + rspan.log_emission);
                    }
                }
            }
        }
    }
    require(inserts_used.size() >= 2,
            "fixture: at least two distinct inserts must be reachable, or the insert probability "
            "plays no part in the comparison");
    require(span_pairs > inserts_used.size(),
            "fixture: several latent span pairs must share endpoints, or the grouping claim is "
            "untested");

    const auto placements =
        panvar::enumerate_placements_exhaustive(haplotype, fragment, support, model, band);
    const double produced =
        std::exp(panvar::log_placement_mass(placements, fragment, insert, model));
    require_close(produced, reference, 1e-12,
                  "the joint likelihood must equal the component-wise reference: grouping latent "
                  "alignments onto physical origins may not change the total");

    // At least one origin must genuinely be explained by more than one latent alignment, and the two
    // inserts involved must carry different probabilities, or the gate would pass on a fixture where
    // the distinction does not arise.
    std::size_t multi = 0;
    for (const panvar::FragmentPlacement& p : placements) {
        if (p.latent_alignments > 1) ++multi;
    }
    require(multi > 0, "fixture: some origin must have several latent alignments");
    const std::size_t lo = *inserts_used.begin();
    const std::size_t hi = *inserts_used.rbegin();
    require(std::abs(insert.log_probability(lo) - insert.log_probability(hi)) > 0.1,
            "fixture: the reachable inserts must carry visibly different probabilities");

    const auto collapsed = collapse_to_best_span(haplotype, fragment, support, model, band);
    const double collapsed_mass =
        std::exp(panvar::log_placement_mass(collapsed, fragment, insert, model));
    require(std::abs(collapsed_mass - produced) > 1e-12 * std::max(1.0, produced),
            "with the insert floor binding, the best-span collapse must give a different answer, "
            "or this fixture does not exercise the defect");
    // Printed as DIFFERENCES, not as three values. At six significant figures these masses are
    // indistinguishable, so printing them side by side made a passing gate look exactly like a gate
    // that detects nothing -- and a gate whose evidence cannot be read is one nobody can check.
    // Measured margins on this fixture: the collapse differs by ~3e-9 absolute (4e-7 relative)
    // against a 1e-12 threshold, while at or above the safe floor it differs by ~4e-18, i.e. zero.
    std::cout << "  component factorisation over " << inserts_used.size() << " distinct inserts: "
              << "componentised - reference = " << std::abs(produced - reference)
              << " (must be 0), collapse - componentised = " << std::abs(collapsed_mass - produced)
              << " (must be > 0)\n";
}

// -------------------------------------------------------------------------------------------
// The optimised emission against the committed log-space reference.
//
// This gate carries weight that none of the others can. enumerate_placements_exhaustive calls
// gapped_emission, which is now the SCALED path, so the exposure identity, the factorisation gate
// and the scorer-versus-oracle comparison all exercise the optimised arithmetic on both sides and
// cannot tell it apart from the reference. gapped_emission_reference is kept in the library so that
// this comparison is against the implementation the model was actually gated on.
//
// Identical REFUSAL is the prefilter's contract. The prefilter minimises edits over every target
// prefix, a superset of the consumed spans the model admits, so its minimum can only undershoot the
// constrained one -- it may over-admit, never falsely reject. The fixtures below therefore walk the
// edit count from zero through band + 2, so both sides of that boundary are exercised.
// -------------------------------------------------------------------------------------------
void test_scaled_emission_matches_the_reference() {
    std::mt19937_64 rng(90210);
    std::string haplotype(400, 'A');
    for (char& base : haplotype) base = "ACGT"[rng() % 4];
    // A homopolymer, so one-base indels are genuinely ambiguous rather than unique.
    haplotype.replace(150, 10, std::string(10, 'T'));

    std::size_t compared = 0;
    std::size_t accepted = 0;
    std::size_t refused = 0;
    std::size_t multi_component = 0;
    double worst_component = 0.0;
    double worst_total = 0.0;

    for (double error : {0.002, 0.02}) {
        for (std::size_t max_gap : {0, 1, 3}) {
            const panvar::LibraryModel model = gap_library(error, max_gap);
            for (std::size_t band : {2, 6}) {
                // Interior starts, a boundary start, and one inside the homopolymer.
                for (std::size_t start : {0, 37, 148, 290, 400 - 60}) {
                    const std::size_t length = 60;
                    if (start + length > haplotype.size()) continue;
                    for (std::size_t edits = 0; edits <= band + 2; ++edits) {
                        std::string read = haplotype.substr(start, length);
                        // Substitutions at spread positions, then a one-base indel for the widest
                        // cases so gapped shapes carry real mass rather than only substitutions.
                        for (std::size_t e = 0; e < edits && e < length; ++e) {
                            const std::size_t at = (e * 7 + 3) % length;
                            read[at] = read[at] == 'A' ? 'C' : 'A';
                        }
                        if (max_gap > 0 && edits > 0 && start == 148) {
                            read.insert(read.begin() + 5, 'T');   // homopolymer insertion
                        }

                        const auto fast =
                            panvar::gapped_emission(haplotype, read, start, model, band);
                        const auto reference =
                            panvar::gapped_emission_reference(haplotype, read, start, model, band);
                        ++compared;
                        require(fast.has_value() == reference.has_value(),
                                "the optimised emission must accept and refuse exactly where the "
                                "reference does: a prefilter that rejects a placement carrying mass "
                                "is a false negative, which this model cannot tolerate");
                        if (!fast.has_value()) { ++refused; continue; }
                        ++accepted;
                        if (fast->components.size() > 1) ++multi_component;
                        require(fast->components.size() == reference->components.size(),
                                "the two must report the same consumed-span classes");
                        for (std::size_t c = 0; c < fast->components.size(); ++c) {
                            require(fast->components[c].consumed_bp ==
                                        reference->components[c].consumed_bp,
                                    "component consumed spans must match exactly");
                            require(fast->components[c].edits == reference->components[c].edits,
                                    "component edit counts must match exactly");
                            const double gap = std::abs(fast->components[c].log_emission -
                                                        reference->components[c].log_emission);
                            worst_component = std::max(worst_component, gap);
                        }
                        worst_total =
                            std::max(worst_total, std::abs(fast->log_emission -
                                                           reference->log_emission));
                        require(fast->consumed_bp == reference->consumed_bp,
                                "the highest-mass consumed span must agree");
                    }
                }
            }
        }
    }
    // A declared tolerance, not a fitted one: scaled probability space differs from log space only
    // by floating-point rounding, and 1e-9 in a log mass is far below any separation that decides a
    // call while still being orders of magnitude above the rounding actually observed.
    require(worst_component < 1e-9,
            "every consumed-span component must match the reference within the declared tolerance");
    require(worst_total < 1e-9, "the whole-read emission must match the reference likewise");
    // Non-vacuity, all three: the sweep must accept, must refuse, and must produce gapped classes.
    require(accepted > 0 && refused > 0 && multi_component > 0,
            "the sweep must exercise acceptance, refusal and multi-component placements, or it "
            "certifies far less than it appears to");
    std::cout << "  scaled emission matches the log-space reference on " << compared
              << " placements (" << accepted << " accepted, " << refused << " refused, "
              << multi_component << " multi-component); worst component delta " << worst_component
              << ", worst total delta " << worst_total << '\n';
}

// The cache must answer exactly what the reference computes, including when every key collides.
// Byte equality is the only thing that separates keys under the collapse, so if the comparison were
// wrong -- or if a component copy were a view into an arena another insert reallocated -- this is
// where it shows. Each key is requested twice so the second request is served from the cache.
void test_emission_cache_matches_the_reference_under_hash_collapse() {
    std::string haplotype(300, 'A');
    for (std::size_t i = 0; i < haplotype.size(); ++i) haplotype[i] = "ACGT"[(i * 5 + 1) % 4];
    haplotype.replace(120, 9, std::string(9, 'T'));       // a homopolymer for real indel ambiguity

    // Held alive for the whole cached section: keys reference these buffers.
    std::vector<std::string> reads;
    std::vector<std::size_t> starts;
    for (std::size_t start : {0, 44, 118, 121, 200, 240}) {
        if (start + 50 > haplotype.size()) continue;
        for (std::size_t edits = 0; edits <= 3; ++edits) {
            std::string bases = haplotype.substr(start, 50);
            for (std::size_t e = 0; e < edits; ++e) {
                const std::size_t at = (e * 11 + 2) % bases.size();
                bases[at] = bases[at] == 'A' ? 'C' : 'A';
            }
            reads.push_back(bases);
            starts.push_back(start);
        }
        // The NA19240 mechanism itself, through the cache: one extra base, and one missing.
        std::string inserted = haplotype.substr(start, 50);
        inserted.insert(inserted.begin() + 5, inserted[5]);
        inserted.resize(50);
        reads.push_back(inserted);
        starts.push_back(start);
        std::string deleted = haplotype.substr(start, 51);
        deleted.erase(deleted.begin() + 5);
        deleted.resize(50);
        reads.push_back(deleted);
        starts.push_back(start);
        // Unrelated sequence, far outside any band. Without these the sweep refuses nothing once
        // the displacement band opens -- the indel-bearing reads that a gapless model turned away
        // now place -- and the CACHED REFUSAL path would go untested at exactly the settings that
        // matter. Derived deterministically from the start so the fixture stays reproducible.
        std::string unrelated(50, 'A');
        for (std::size_t i = 0; i < unrelated.size(); ++i) {
            unrelated[i] = "ACGT"[(start * 7 + i * 3 + 1) % 4];
        }
        reads.push_back(unrelated);
        starts.push_back(start);
    }

    for (std::size_t max_gap : {0, 2, 3}) {
        const panvar::LibraryModel model = gap_library(0.002, max_gap);
        const std::size_t band = 6;
        std::vector<std::optional<panvar::GappedFit>> expected;
        expected.reserve(reads.size());
        for (std::size_t i = 0; i < reads.size(); ++i) {
            expected.push_back(
                panvar::gapped_emission_reference(haplotype, reads[i], starts[i], model, band));
        }

        panvar::set_emission_cache_hash_collapse(true);
        panvar::begin_emission_cache(1u << 20);
        std::size_t accepted = 0;
        std::size_t refused = 0;
        double worst = 0.0;
        for (int pass = 0; pass < 2; ++pass) {     // second pass must come from the cache
            for (std::size_t i = 0; i < reads.size(); ++i) {
                const auto got =
                    panvar::gapped_emission(haplotype, reads[i], starts[i], model, band);
                require(got.has_value() == expected[i].has_value(),
                        "a totally collided cache must still accept and refuse exactly as the "
                        "reference does");
                if (!got.has_value()) { if (pass == 0) ++refused; continue; }
                if (pass == 0) ++accepted;
                require(got->components.size() == expected[i]->components.size(),
                        "consumed-span components must survive the cache intact");
                for (std::size_t c = 0; c < got->components.size(); ++c) {
                    require(got->components[c].consumed_bp == expected[i]->components[c].consumed_bp &&
                                got->components[c].edits == expected[i]->components[c].edits,
                            "component spans and edit counts must be identical");
                    worst = std::max(worst, std::abs(got->components[c].log_emission -
                                                     expected[i]->components[c].log_emission));
                }
                require(got->consumed_bp == expected[i]->consumed_bp,
                        "the highest-mass span must be identical");
                worst = std::max(worst, std::abs(got->log_emission - expected[i]->log_emission));
            }
        }
        const panvar::EmissionCacheCounters counters = panvar::end_emission_cache();
        panvar::set_emission_cache_hash_collapse(false);

        require(worst < 1e-9, "a cached emission must match the reference within the tolerance");
        require(accepted > 0 && refused > 0,
                "fixture: the sweep must both accept and refuse, or the collision path is only "
                "half exercised");
        // Two corrections over the obvious threshold, both measured rather than assumed. Only
        // recursion results are cached, so the unrelated reads here -- turned away before a key
        // exists -- can never hit. And reads DUPLICATE: at max_gap 0 the 24 accepted reads collapse
        // to 12 distinct keys, at max_gap 2 and 3 the 36 collapse to 18. So the bound is one hit
        // per distinct accepted KEY, counted here, not one per accepted read.
        std::set<std::pair<std::string, std::string>> distinct_accepted;
        for (std::size_t i = 0; i < reads.size(); ++i) {
            if (!expected[i].has_value()) continue;
            const std::size_t window =
                std::min(haplotype.size() - starts[i], reads[i].size() + max_gap);
            distinct_accepted.insert({reads[i], haplotype.substr(starts[i], window)});
        }
        require(!distinct_accepted.empty(), "fixture: some read must be accepted");
        require(counters.dp_cache_hits >= static_cast<std::uint64_t>(distinct_accepted.size()),
                "every distinct accepted key must be served from the cache on the second pass, or "
                "nothing was reused");
        // Every request pays a prefilter, so the old "executed work < requests" form is false by
        // construction. What reuse actually buys is recursions avoided among the requests that got
        // as far as a lookup.
        require(counters.executed_dps < counters.dp_cache_requests,
                "reuse must have avoided recursions, not merely recorded them");
        std::cout << "  hash-collapsed cache matches the reference at max_gap=" << max_gap << ": "
                  << accepted << " accepted, " << refused << " refused, " << counters.dp_cache_hits
                  << " hits, worst delta " << worst << '\n';
    }
}

// The prepared Myers pattern must make the IDENTICAL admit/reject decision as building the masks on
// the spot, and both must agree with edlib. Preparation moves work out of the hot path; if it moved
// meaning as well, every placement downstream would be wrong in a way no score comparison would
// localise. The sweep walks the edit count through band + 3 so both sides of the boundary are
// exercised, and includes indel-bearing and unrelated reads.
void test_prepared_myers_pattern_matches_the_unprepared_decision() {
    std::mt19937_64 rng(31337);
    std::string haplotype(4000, 'A');
    for (char& base : haplotype) base = "ACGT"[rng() % 4];
    haplotype.replace(1500, 12, std::string(12, 'G'));

    const std::size_t length = 100;
    const std::size_t band = 6;
    const panvar::LibraryModel model = gap_library(0.002, 3);

    std::size_t compared = 0;
    std::size_t admitted = 0;
    std::size_t refused = 0;
    for (int trial = 0; trial < 4000; ++trial) {
        const std::size_t start = rng() % (haplotype.size() - length - 8);
        std::string read;
        if (trial % 9 == 8) {
            read.assign(length, 'A');                       // unrelated: far outside any band
            for (char& base : read) base = "ACGT"[rng() % 4];
        } else {
            read = haplotype.substr(start, length);
            const std::size_t edits = trial % (band + 4);
            for (std::size_t e = 0; e < edits; ++e) {
                const std::size_t at = (e * 17 + 5) % length;
                read[at] = read[at] == 'A' ? 'C' : 'A';
            }
            if (trial % 5 == 0 && edits > 0) {
                read.insert(read.begin() + 20, read[20]);    // one-base insertion
                read.resize(length);
            }
        }

        const panvar::MyersPattern pattern = panvar::prepare_myers_pattern(read);
        require(pattern.usable && pattern.length == length,
                "fixture: a 100 bp read must yield a usable prepared pattern");
        // nullptr rebuilds the masks; the pattern uses the prepared ones. Same decision required.
        const auto unprepared =
            panvar::gapped_emission(haplotype, read, start, model, band, nullptr);
        const auto prepared =
            panvar::gapped_emission(haplotype, read, start, model, band, &pattern);
        require(unprepared.has_value() == prepared.has_value(),
                "the prepared pattern must accept and refuse exactly where rebuilding does");
        ++compared;
        if (!prepared.has_value()) { ++refused; continue; }
        ++admitted;
        require(prepared->components.size() == unprepared->components.size(),
                "preparation must not change the consumed-span decomposition");
        require(prepared->log_emission == unprepared->log_emission,
                "preparation must not change the emission by even a bit");
    }
    require(compared == 4000, "the mutation set must run every case");
    require(admitted > 0 && refused > 0,
            "the sweep must both admit and refuse, or one side of the decision is untested");
    std::cout << "  prepared Myers pattern matches the unprepared decision on " << compared
              << " cases (" << admitted << " admitted, " << refused << " refused)\n";
}

int main() {
    try {
        test_prepared_myers_pattern_matches_the_unprepared_decision();
        test_scaled_emission_matches_the_reference();
        test_emission_cache_matches_the_reference_under_hash_collapse();
        test_gapped_reduces_to_substitution_model();
        test_gapped_emission_normalises();
        test_one_base_homopolymer_insertion();
        test_one_base_deletion();
        test_signal_density_integrates_to_exposure();
        test_best_span_collapse_is_rejected();
        test_component_factorisation_and_collapse();
        test_emission();
        test_fr_geometry();
        test_exposure();
        test_fragment_contribution();
        test_single_origin();
        test_both_orientations();
        test_repeat_multiplicity();
        test_band_is_model();
        test_homozygote_counts_two_copies();
        test_order_and_uniqueness();
        test_truncated_emission_normalises();
        test_complete_intensity_integrates();
        std::cout << "genotype fragment oracle: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "genotype fragment oracle: FAIL: " << error.what() << '\n';
        return 1;
    }
}
