#pragma once

#include "panvar/genotype_model.hpp"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <string>
#include <vector>

namespace panvar {

// The single definition of read emission, insert probability, valid-FR geometry, exposure, fragment
// contribution and log-addition. Every caller uses these; a private copy of any of them is how the
// accelerated scorer and its oracle start disagreeing about the model they are supposed to share.

double log_add_exp(double lhs, double rhs);

// Which mate is read off the haplotype's forward strand. The other mate matches the haplotype's
// reverse strand, which is the same as its own reverse complement matching forward, so every
// comparison below is a forward comparison once the reads are pre-complemented.
enum class MateOrientation { FirstForward, SecondForward };

// A single physical origin of a fragment on one haplotype: where each mate sits and how far apart.
// Two distinct repeat copies produce two of these even when every statistic agrees, because they are
// two places the fragment could physically have come from.
struct FragmentPlacement {
    std::size_t forward_start = 0;
    std::size_t reverse_start = 0;
    std::size_t insert_bp = 0;
    std::size_t forward_mismatches = 0;
    std::size_t reverse_mismatches = 0;
    MateOrientation orientation = MateOrientation::FirstForward;
};

// Fixed-position Hamming emission of `read` against `haplotype` starting at `start`. Version one is
// substitutions only. Returns nullopt when the read runs off the end; a non-ACGT base in either
// sequence refuses, because silently scoring it as a mismatch would turn missing data into evidence.
std::optional<std::size_t> count_mismatches(
    const std::string& haplotype,
    const std::string& read,
    std::size_t start,
    std::size_t mismatch_limit = std::numeric_limits<std::size_t>::max());

double log_emission(
    std::size_t read_length,
    std::size_t mismatches,
    double substitution_error);

// The template span from the first outer boundary to the second. Overlapping mates are physical
// states, so the floor is the longer mate, never the sum. Returns nullopt when the geometry is not a
// valid FR pair or the span falls outside this fragment's declared support.
std::optional<std::size_t> fr_insert_length(
    std::size_t forward_start,
    std::size_t forward_length,
    std::size_t reverse_start,
    std::size_t reverse_length,
    const InsertSupport& support);

// Expected number of valid fragment placements on one chromosome copy of this haplotype, under the
// same normalized insert distribution the mass term uses. Analytic: sum_i P(i) * (L - i + 1)^+.
double haplotype_exposure(
    std::size_t haplotype_length,
    const DiscreteInsertDistribution& insert);

// log((1-eta) * lambda * mass + eta * P_background(f)). `log_mass` is the log of the summed
// placement mass over BOTH homologues, so a homozygote has already counted two chromosome copies.
double log_fragment_contribution(
    double log_mass,
    double log_background,
    const LibraryModel& library);

// ---------------------------------------------------------------------------------------------
// Reference implementation.
// ---------------------------------------------------------------------------------------------

// Every physical origin of `fragment` on `haplotype`, found by trying every start pair in both
// library orientations. This is the oracle: it is O(L^2) per fragment and exists to certify the
// accelerated scorer on small fixtures, never to run on a real locus.
//
// `mismatch_limit` is the declared band. It is part of the model, not an optimisation: placements
// beyond it are omitted from the mass by both this and the accelerated path, so both must apply the
// same band or their masses are not comparable.
std::vector<FragmentPlacement> enumerate_placements_exhaustive(
    const std::string& haplotype,
    const ReadPair& fragment,
    const InsertSupport& support,
    std::size_t mismatch_limit);

// Summed log mass of a fragment on one chromosome copy: sum over placements of
// P(insert) * P(mate bases | haplotype). Returns -inf when no placement is accepted.
double log_placement_mass(
    const std::vector<FragmentPlacement>& placements,
    const ReadPair& fragment,
    const DiscreteInsertDistribution& insert,
    const LibraryModel& library);

} // namespace panvar
