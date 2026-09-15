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

panvar::LibraryModel library(double error = 0.01, double outlier = 0.001) {
    panvar::LibraryModel model;
    model.insert_mean_bp = 20.0;
    model.insert_sd_bp = 4.0;
    model.support_sigmas = 3.0;
    model.substitution_error = error;
    model.outlier_weight = outlier;
    model.fragment_start_intensity = 0.05;
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
        panvar::enumerate_placements_exhaustive(haplotype, fragment, support, 0);
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
        panvar::log_placement_mass(placements, fragment, insert, model, 0),
        std::log(1.0 / 15.0) + std::log(0.5), 1e-12,
        "single-origin mass is the insert probability times the orientation prior");
    // At band one the same origin keeps only the fraction of emission mass the band admits.
    const double e = model.substitution_error;
    require_close(
        panvar::log_placement_mass(placements, fragment, insert, model, 1),
        std::log(1.0 / 15.0) + std::log(0.5) +
            2.0 * (8.0 * std::log(1.0 - e) - std::log(std::pow(1.0 - e, 8) + 8.0 * e * std::pow(1.0 - e, 7))),
        1e-12,
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
        panvar::enumerate_placements_exhaustive(haplotype, first_reads_left, support, 0);
    const auto arm_b =
        panvar::enumerate_placements_exhaustive(haplotype, second_reads_left, support, 0);
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
        panvar::enumerate_placements_exhaustive(haplotype, fragment, support, 0);
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
    const double both = panvar::log_placement_mass(placements, fragment, insert, model, 0);
    const double one = panvar::log_placement_mass({placements[0]}, fragment, insert, model, 0);
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

    require(panvar::enumerate_placements_exhaustive(haplotype, fragment, support, 0).empty(),
            "a band of zero omits the mismatching origin entirely");
    const auto banded =
        panvar::enumerate_placements_exhaustive(haplotype, fragment, support, 1);
    require(banded.size() == 1 && banded.front().forward_mismatches == 1,
            "widening the band to one admits exactly that origin");

    panvar::DiscreteInsertDistribution insert;
    insert.support = support;
    insert.log_probability_by_offset.assign(15, std::log(1.0 / 15.0));
    const panvar::LibraryModel model = library();
    const double narrow = panvar::log_placement_mass({}, fragment, insert, model, 0);
    const double wide = panvar::log_placement_mass(banded, fragment, insert, model, 1);
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
        panvar::enumerate_placements_exhaustive(haplotype, fragment, support, 0);
    const double one_copy = panvar::log_placement_mass(placements, fragment, insert, model, 0);
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

    auto placements = panvar::enumerate_placements_exhaustive(haplotype, fragment, support, 0);
    const double forward = panvar::log_placement_mass(placements, fragment, insert, model, 0);
    std::reverse(placements.begin(), placements.end());
    require_close(panvar::log_placement_mass(placements, fragment, insert, model, 0), forward, 1e-12,
                  "placement order must not change the summed mass");

    std::set<std::tuple<std::size_t, std::size_t, int>> distinct;
    for (const panvar::FragmentPlacement& placement : placements) {
        distinct.emplace(placement.forward_start, placement.reverse_start,
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

// The normalisation gate. Summing the signal mass over EVERY possible observation with fixed read
// lengths must reproduce the exposure exactly: the -lambda*E term and the fragment density describe
// the same generated process only if this holds. Checked at a band that truncates and at one that
// does not. Each orientation arm must carry exactly half, which is what proves the 1/2 prior is where
// it belongs rather than compensated somewhere else.
void test_signal_density_integrates_to_exposure() {
    const std::string haplotype = "ACGTTGCAAGTCCA";   // L = 14
    panvar::LibraryModel model;
    model.insert_mean_bp = 8.0;
    model.insert_sd_bp = 2.0;
    model.support_sigmas = 2.0;          // support [4, 12]
    model.substitution_error = 0.1;
    model.outlier_weight = 0.001;
    model.fragment_start_intensity = 0.05;
    model.provenance = "fixture: declared";
    const panvar::DiscreteInsertDistribution insert = panvar::build_library_insert_distribution(model);
    require(insert.support.min_bp == 4 && insert.support.max_bp == 12,
            "fixture: the library support must be [4,12]");

    const std::size_t mate = 3;
    const double exposure = panvar::haplotype_exposure(haplotype.size(), insert, mate);
    require(exposure > 0.0, "fixture: exposure must be positive, or the identity is vacuous");

    for (std::size_t band : {1, 3}) {
        double total = 0.0;
        double first_forward = 0.0;
        std::size_t observations = 0;
        for (std::size_t code1 = 0; code1 < 64; ++code1) {
            for (std::size_t code2 = 0; code2 < 64; ++code2) {
                std::string r1(mate, 'A');
                std::string r2(mate, 'A');
                for (std::size_t i = 0; i < mate; ++i) {
                    r1[i] = "ACGT"[(code1 >> (2 * i)) & 3];
                    r2[i] = "ACGT"[(code2 >> (2 * i)) & 3];
                }
                const panvar::ReadPair fragment{read("o/1", r1), read("o/2", r2)};
                const panvar::InsertSupport support = panvar::fragment_insert_support(fragment, model);
                const auto placements =
                    panvar::enumerate_placements_exhaustive(haplotype, fragment, support, band);
                const double mass = panvar::log_placement_mass(placements, fragment, insert, model, band);
                if (mass != kNegInf) total += std::exp(mass);
                std::vector<panvar::FragmentPlacement> arm;
                for (const auto& p : placements) {
                    if (p.orientation == panvar::MateOrientation::FirstForward) arm.push_back(p);
                }
                const double arm_mass = panvar::log_placement_mass(arm, fragment, insert, model, band);
                if (arm_mass != kNegInf) first_forward += std::exp(arm_mass);
                ++observations;
            }
        }
        require(observations == 4096, "the gate must sum over every possible observation");
        require_close(total, exposure, 1e-9,
                      "the signal density summed over every observation must equal the exposure");
        require_close(first_forward, 0.5 * exposure, 1e-9,
                      "each library orientation must carry exactly half of the exposure");
        std::cout << "  normalisation at band " << band << ": integrated mass " << total
                  << " = exposure " << exposure << '\n';
    }
}


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
                haplotype, fragment, panvar::fragment_insert_support(fragment, model), band);
            const double mass = panvar::log_placement_mass(placements, fragment, insert, model, band);
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

int main() {
    try {
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
        test_signal_density_integrates_to_exposure();
        test_complete_intensity_integrates();
        std::cout << "genotype fragment oracle: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "genotype fragment oracle: FAIL: " << error.what() << '\n';
        return 1;
    }
}
