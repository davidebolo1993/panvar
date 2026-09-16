#pragma once

#include "panvar/genotype_model.hpp"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace panvar {

// The single definition of read emission, insert probability, valid-FR geometry, exposure, fragment
// contribution and log-addition. Every caller uses these; a private copy of any of them is how the
// accelerated scorer and its oracle start disagreeing about the model they are supposed to share.

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

// ---------------------------------------------------------------------------------------------
// The declared fragment model: a truncated likelihood.
//
// The mismatch band d is part of the model and there is exactly one definition of it. A mate is
// generated with independent substitutions at rate e, CONDITIONED on carrying at most d of them, so
// its emission is (1-e)^(L-m) (e/3)^m / F(d,L) for m <= d and zero beyond, where F(d,L) is the
// binomial probability of at most d substitutions. Under that model the in-band mass is not an
// approximation of anything: it is the exact likelihood, a fragment with no in-band placement on a
// haplotype has zero mass there, and neither fact depends on which other candidates are scored.
//
// Each physical template is sequenced in one of two library orientations with probability 1/2, the
// insert length follows the library distribution normalised over the library's own support, and
// the fragment-start intensity lambda counts templates per base per chromosome copy summed over both
// orientations. With those three definitions the signal density integrates exactly to the exposure:
// summing M_f,h over every possible observation with fixed read lengths gives E_h. That identity is
// gated directly on a tiny locus.
// ---------------------------------------------------------------------------------------------

// log F(d, L): the probability that a mate of length L carries at most d substitutions.
double log_band_mass(std::size_t read_length, std::size_t mismatch_band, double substitution_error);

// Emission of a mate under the truncated model. -inf beyond the band.
double log_truncated_emission(
    std::size_t read_length,
    std::size_t mismatches,
    double substitution_error,
    std::size_t mismatch_band);

// log(1/2): each template is read in one of two library orientations.
inline constexpr double kLogMateOrientationPrior = -0.69314718055994530942;

// The template span from the first outer boundary to the second. Overlapping mates are physical
// states, so the floor is the longer mate, never the sum. Returns nullopt when the geometry is not a
// valid FR pair or the span falls outside this fragment's declared support.
std::optional<std::size_t> fr_insert_length(
    std::size_t forward_start,
    std::size_t forward_length,
    std::size_t reverse_start,
    std::size_t reverse_length,
    const InsertSupport& support);

// Expected number of templates on one chromosome copy of this haplotype that can produce a fragment
// whose longer mate is `physical_floor_bp`, under the library insert distribution:
// sum over i >= floor of P(i) * (L - i + 1)^+. Analytic, and summed over both orientations because
// lambda is.
double haplotype_exposure(
    std::size_t haplotype_length,
    const DiscreteInsertDistribution& insert,
    std::size_t physical_floor_bp = 0);

// log((1-eta) * lambda * mass + eta * P_background(f)). `log_mass` is the log of the summed
// placement mass over BOTH homologues, so a homozygote has already counted two chromosome copies.
double log_fragment_contribution(
    double log_mass,
    double log_background,
    const LibraryModel& library);

// The candidate-dependent part of the integrated intensity, which is what the -E term of the score
// must be. The per-fragment intensity is the superposition of a signal process,
// (1-eta) * lambda * (M_a + M_b), and a background process, eta * P_background(f). Integrated over
// every possible observation the signal gives (1-eta) * lambda * (E_a + E_b) and the background gives
// eta per read-length class, the same for every candidate and therefore dropped. The (1-eta) factor
// is NOT candidate-independent: scoring with lambda * (E_a + E_b) instead ranks haplotypes of different
// lengths against a density that does not integrate to its own exposure.
double expected_signal_fragments(double exposure_sum, const LibraryModel& library);

// ---------------------------------------------------------------------------------------------
// Reference implementation.
// ---------------------------------------------------------------------------------------------

// Every physical origin of `fragment` on `haplotype`, found by trying every start pair in both
// library orientations. This is the oracle: it is O(L^2) per fragment and exists to certify the
// accelerated scorer on small fixtures, never to run on a real locus.
//
// `mismatch_limit` is the declared band of the truncated model. Placements beyond it have zero
// emission under that model, so they are absent from the enumeration by definition.
std::vector<FragmentPlacement> enumerate_placements_exhaustive(
    const std::string& haplotype,
    const ReadPair& fragment,
    const InsertSupport& support,
    std::size_t mismatch_limit);

// Summed log mass of a fragment on one chromosome copy under the truncated model: sum over
// placements of (1/2) * P_library(insert) * truncated emission of both mates. `insert` must be the
// library distribution. Returns -inf when no placement is accepted.
double log_placement_mass(
    const std::vector<FragmentPlacement>& placements,
    const ReadPair& fragment,
    const DiscreteInsertDistribution& insert,
    const LibraryModel& library,
    std::size_t mismatch_band);

// ---------------------------------------------------------------------------------------------
// Accelerated scorer.
// ---------------------------------------------------------------------------------------------

// K diplotype candidates contain at most 2K haploid mosaics, and mosaics that differ in their allele
// path can still spell identical bytes. Scoring the unique sequences once and combining two columns
// per candidate is what keeps the cost linear in unique haplotypes rather than quadratic in
// candidates.
struct InternedHaplotypes {
    std::vector<std::string> sequences;
    // Parallel to the input: which unique sequence each spelled haplotype resolved to.
    std::vector<std::size_t> index_of;
    std::size_t distinct_allele_paths = 0;
};

InternedHaplotypes intern_haplotypes(
    const std::vector<MosaicPath>& allele_paths,
    const std::vector<std::string>& sequences);

struct HaplotypeCandidate {
    std::size_t first = 0;    // index into InternedHaplotypes::sequences
    std::size_t second = 0;
};

// Thrown when a declared resource cap is reached. Scoring stops and returns nothing: a ranking built
// from the placements found before the cap would silently favour whichever haplotypes were scanned
// first, which is a partial ranking presented as a result.
class ResourceExhausted : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

struct FragmentScoringOptions {
    // The declared band of the truncated model. Fixed for a run: changing it changes the model.
    std::size_t mismatch_band = 6;
    // Pigeonhole pieces shorter than this are too unspecific to index, so the scorer refuses rather
    // than silently degrading into a quadratic scan.
    std::size_t minimum_piece_bp = 10;
    // Operational caps, checked as the work happens rather than estimated beforehand. The defaults
    // leave two orders of magnitude of headroom over C4 at 30x and K=64 (14.8M piece hits, 1.9M
    // proposed starts per haplotype, 0.88M physical origins) while stopping a repeat-rich locus such
    // as LPA before it exhausts memory. Reaching one throws ResourceExhausted.
    std::uint64_t max_piece_hits = 2000000000ULL;
    std::uint64_t max_proposed_starts_per_haplotype = 50000000ULL;
    std::uint64_t max_physical_origins = 500000000ULL;
    std::uint64_t max_state_classes = 200000000ULL;
    std::uint64_t max_mass_cells = 200000000ULL;
    // Diagnostics only; nothing here enters a score or a ranking. When set, the result carries for
    // every (fragment, haplotype) the accepted physical-origin count and the single highest-mass origin,
    // and for the listed fragments every accepted origin.
    bool record_origin_summary = false;
    std::vector<std::size_t> record_origins_for_fragments;
    // Workers for the diploid combination phase. Candidates are partitioned into fixed contiguous
    // ranges and each worker writes only its own complete DiplotypeScore entries, so the result does
    // not depend on the thread count: no cross-thread floating-point reduction exists to reorder.
    // Defaults to serial, so every existing caller keeps its behaviour without being touched.
    std::size_t threads = 1;
};

struct OriginSummary {
    std::uint64_t origins = 0;
    FragmentPlacement best;
    double best_log_term = -std::numeric_limits<double>::infinity();
};

// Measured work, not estimated viability. `candidate_bases_scanned` counts haplotype bases actually
// walked, which is the honest denominator; counting hash lookups instead flatters the index.
struct FragmentScoringCounters {
    std::uint64_t candidate_bases_scanned = 0;
    std::uint64_t piece_hits = 0;
    std::uint64_t proposed_starts = 0;
    std::uint64_t unique_starts = 0;
    std::uint64_t full_read_verifications = 0;
    std::uint64_t accepted_mate_placements = 0;
    std::uint64_t valid_fr_states = 0;
    std::uint64_t physical_origins = 0;
    std::uint64_t statistical_state_classes = 0;
};

// Exact relative to the declared truncated model. `absolute_term_sum` is the sum of the magnitudes
// of every term added into the score, which is what bounds its floating-point accumulation error.
struct DiplotypeScore {
    double log_likelihood = 0.0;
    double fragment_term = 0.0;
    double exposure_term = 0.0;
    double absolute_term_sum = 0.0;
    std::size_t fragments_with_mass = 0;
};

struct FragmentScoringResult {
    std::vector<DiplotypeScore> scores;     // parallel to the requested candidates
    // Row-major [fragment * unique_haplotypes + haplotype]: the summed log placement mass of one
    // fragment on one chromosome copy of one haplotype. Retained because it is what the exhaustive
    // oracle certifies directly; recovering it from a candidate score would compare the two paths
    // through the mixture rather than at the mass itself.
    std::vector<double> log_mass;
    // Per unique haplotype: expected templates on one chromosome copy, weighted by the observed
    // distribution of fragment physical floors.
    std::vector<double> exposure;
    std::size_t unique_haplotypes = 0;
    FragmentScoringCounters counters;
    // Wall clock of the three phases inside one call, so a caller can attribute rescoring instead of
    // reporting it whole. Preparation is per fragment; placement and mass construction is per
    // (fragment, haplotype); diploid combination is per (candidate, fragment) and is the term that
    // grows with the panel. Which of the last two dominates decides how this should be parallelised.
    double seconds_fragment_preparation = 0.0;
    double seconds_placement_and_mass = 0.0;
    double seconds_diploid_combination = 0.0;
    // What was actually used, which is not what was asked for when there are fewer candidates than
    // workers. Reported so a run can never claim parallelism it did not use.
    std::size_t threads_used_combination = 1;
    // candidates x fragments: the work the combination phase performs, independent of thread count.
    std::uint64_t candidate_fragment_evaluations = 0;
    // Filled only when requested in FragmentScoringOptions; same row-major layout as log_mass.
    std::vector<OriginSummary> origin_summary;
    // (fragment, haplotype) -> every accepted origin, for record_origins_for_fragments only.
    std::map<std::pair<std::size_t, std::size_t>, std::vector<FragmentPlacement>> recorded_origins;

    double mass_of(std::size_t fragment, std::size_t haplotype) const;
};

// Scores exactly the requested candidates and nothing else. The raw score of a candidate does not
// depend on which other candidates were requested, on their order, or on any proposal score.
FragmentScoringResult score_diplotypes(
    const InternedHaplotypes& haplotypes,
    const std::vector<HaplotypeCandidate>& candidates,
    const std::vector<ReadPair>& fragments,
    const LibraryModel& library,
    const FragmentScoringOptions& options);

// ---------------------------------------------------------------------------------------------
// Certification.
// ---------------------------------------------------------------------------------------------

enum class RescoreStatus { Complete, RescoreIncomplete };

// Every score is exact relative to the declared model, so the only uncertainty left in a comparison
// is floating-point accumulation. `numeric_error_bound` is a conservative first-order bound on it:
// (terms + 16) * DBL_EPSILON * sum |term| for a recursively summed score. It is not a model interval
// and it never widens with the band.
struct CertifiedCandidate {
    std::size_t candidate_id = 0;      // index into the requested candidate list
    double score = 0.0;                // fragment term + exposure term + exact prior
    double fragment_term = 0.0;
    double exposure_term = 0.0;
    double prior = 0.0;
    double numeric_error_bound = 0.0;
    std::size_t fragments_with_mass = 0;
    // exp(score - log sum_c exp(score_c)) over the retained candidates. A probability within the
    // candidate set only: omitted candidates carry mass this does not see, so it is not a posterior.
    double candidate_set_probability = 0.0;
};

struct CertificationOptions {
    // A floor under the numerical tie tolerance, for scores so small that the computed bound is below
    // any meaningful difference.
    double minimum_tie_tolerance_nats = 1e-9;
};

struct CertifiedRanking {
    RescoreStatus status = RescoreStatus::RescoreIncomplete;
    std::vector<CertifiedCandidate> candidates;   // sorted by score, descending
    // One id when a winner is certified; every candidate that cannot be separated from the best when
    // it is not. An undecided ranking reports the set rather than picking from it.
    std::vector<std::size_t> winner_set;
    // min over rivals of (best - its error bound) - (rival + its error bound). Positive exactly when
    // the status is Complete.
    double certified_gap = 0.0;
    // Raw score difference between the best and the runner-up.
    double score_gap_to_runner_up = 0.0;
    double candidate_set_log_weight = 0.0;
    std::size_t mismatch_band = 0;
    FragmentScoringCounters counters;
};

// Scores every requested candidate once under the declared truncated model, adds the exact prior,
// and certifies a winner when its score minus its numerical error bound exceeds every rival's score
// plus theirs. There is no refinement loop: nothing is approximated, so nothing can be deepened.
// A candidate's score depends only on its own two haplotypes, the fragments and the library.
//
// `exact_prior` is added to the matching candidate's score. It must be exact, which is why it is a
// separate argument. The proposal score is not, and may never be, passed here: it comes from the
// same reads as the fragment term, and adding both would count that evidence twice.
// Certification of scores already computed by score_diplotypes: scores[c] is candidate c. A diplotype's
// score does not depend on the other candidates it is scored with, so a subset of one scoring pass is
// ranked exactly as scoring that subset alone would rank it.
CertifiedRanking certify_scored_ranking(
    const std::vector<DiplotypeScore>& scores,
    std::size_t fragment_count,
    std::size_t mismatch_band,
    const FragmentScoringCounters& counters,
    const CertificationOptions& certification,
    const std::vector<double>& exact_prior = {});

CertifiedRanking certify_diplotype_ranking(
    const InternedHaplotypes& haplotypes,
    const std::vector<HaplotypeCandidate>& candidates,
    const std::vector<ReadPair>& fragments,
    const LibraryModel& library,
    const FragmentScoringOptions& options,
    const CertificationOptions& certification,
    const std::vector<double>& exact_prior = {});

} // namespace panvar
