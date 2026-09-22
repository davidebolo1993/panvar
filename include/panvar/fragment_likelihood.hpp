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

// A single physical origin of a fragment on one haplotype.
//
// A physical origin is a PLACE, described by the two endpoints of the template and the library
// orientation: the forward mate is anchored at the left endpoint, the reverse mate at the right
// endpoint, and `insert_bp` is their separation. Two distinct repeat copies produce two of these
// even when every statistic agrees, because they are two places the fragment could physically have
// come from.
//
// An ALIGNMENT is not an origin. Once indels are admitted, one origin is explained by several latent
// alignments -- the same left and right endpoints reached with different numbers of haplotype bases
// consumed by each mate -- and those are competing explanations of one physical event, not extra
// physical events. Their masses are therefore SUMMED into `log_latent_mass`, never counted as
// additional multiplicity. Getting that backwards inflates a repeat-free locus with phantom copies.
struct FragmentPlacement {
    // The template's left endpoint: where the forward mate is anchored.
    std::size_t forward_start = 0;
    // Right endpoint minus left endpoint. The right endpoint is forward_start + insert_bp.
    std::size_t insert_bp = 0;
    MateOrientation orientation = MateOrientation::FirstForward;
    // log sum over every latent alignment pair explaining this origin, of the two mates' normalised
    // emissions. This is the whole of what the origin contributes, before the orientation prior and
    // the insert probability -- both of which apply once per origin, because there is one template.
    double log_latent_mass = -std::numeric_limits<double>::infinity();
    // How many latent alignment pairs were summed. One under a gapless model; up to (2d+1)^2 with a
    // displacement band of d. Reported so the pairing cost is measured rather than assumed.
    std::size_t latent_alignments = 0;

    // ---- diagnostics: the single highest-mass latent alignment of this origin -------------------
    // None of these may be used to recompute a score. Under the indel-aware emission a mismatch count
    // does not determine an emission, and one alignment's consumed span does not carry the origin's
    // mass -- attaching the summed mass to the best span's geometry is exactly the defect this
    // structure exists to prevent.
    std::size_t reverse_start = 0;
    std::size_t forward_consumed_bp = 0;
    std::size_t reverse_consumed_bp = 0;
    std::size_t forward_mismatches = 0;
    std::size_t reverse_mismatches = 0;
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

// ---------------------------------------------------------------------------------------------
// Indel-aware emission, version two.
//
// Version one aligned a read at a fixed start and counted substitutions. A single inserted or
// deleted base shifts every downstream base, so a read spanning a one-base indel accumulates
// roughly L/2 mismatches, leaves the declared band, and gets NO placement at all -- it falls to the
// outlier background on the haplotype that is actually correct, while a haplotype that happens to
// match by chance keeps it. That is not a tolerance problem, it is a misspecified emission, and it
// is what cost the two-edit allele 2,871 nats to a sixteen-edit one on NA19240.
//
// The replacement is a banded affine-gap pair-HMM, conditioned on the observed read length:
//
//   states     M consumes one read base and one haplotype base, I one read base, D one haplotype
//              base. No I<->D transition, so a gap has one kind.
//   emissions  M emits the read base with 1-e when it matches and e/3 otherwise; I emits a read
//              base uniformly at 1/4. Both sum to one over the emitted base, and D emits nothing.
//   band       the total edit budget `mismatch_band` counts substitutions AND gap bases together,
//              and net displacement is bounded by `max_gap_bp`.
//
// Because every emitting state's distribution sums to one over the base it emits, summing the
// unnormalised path mass over all 4^L reads of length L leaves exactly the transition mass of the
// admissible alignment SHAPES -- which depends on the read length and the band, never on the
// haplotype or the read. That sum is the normaliser below, so the emission is a proper distribution
// over reads of the observed length, the exposure identity carries over unchanged, and no score can
// depend on which other candidates are present.
//
// At max_gap_bp == 0 this reduces exactly to log_truncated_emission, which is gated as an equality.

// The shape mass log Z: the same recursion with every emission set to one. Independent of the
// sequences' CONTENT, but not of how much haplotype is available: `available_bp` bounds the bases an
// alignment may consume from the start in question, and the normaliser must honour the same bound the
// numerator does. Normalising against an unbounded haplotype instead leaks mass near the ends --
// shapes dropped from the numerator would still divide the denominator, and the emission would sum to
// less than one there, breaking the exposure identity by a small, load-bearing amount.
double log_gap_shape_mass(
    std::size_t read_length,
    std::size_t available_bp,
    std::size_t mismatch_band,
    const LibraryModel& library);

// One terminal consumed-span class of the alignments of a read at a fixed start.
//
// A gapped read does not occupy read_length haplotype bases, and it does not occupy ONE number of
// them either: with a displacement band of d it may end anywhere in read_length +/- d, and each of
// those spans places the mate's far endpoint somewhere different. Because the spans imply different
// fragment geometries and different insert probabilities, the mass of each has to be carried
// separately as far as the geometry -- summing them first and attaching the total to the best span's
// insert length weights real alignment mass with the wrong insert probability, which is what broke
// the exposure identity.
struct GappedComponent {
    // Haplotype bases this class of alignments consumes from `start`.
    std::size_t consumed_bp = 0;
    // Normalised log mass of every admissible alignment ending with that consumed span.
    double log_emission = -std::numeric_limits<double>::infinity();
    // Fewest total edits -- substitutions plus gap bases -- within this class. Diagnostic.
    std::size_t edits = 0;
};

// Every admissible alignment of one read at one start, grouped by what it consumes.
struct GappedFit {
    // One entry per terminal consumed span carrying mass, ascending. At most 2*max_gap_bp+1 of them.
    std::vector<GappedComponent> components;
    // log sum of the components: the whole read's normalised emission at this start, which is what
    // must sum to one over every read of this length. The components share ONE normaliser -- the
    // shape mass of the read length -- so summing them reproduces this exactly. Normalising each
    // span against its own shape mass would make every span a distribution in its own right and
    // their sum several times one.
    double log_emission = -std::numeric_limits<double>::infinity();
    // The highest-mass component's span and edit count. DIAGNOSTIC ONLY: it is not where the
    // emission lives, and using it as the geometry of the whole emission is the defect above.
    std::size_t consumed_bp = 0;
    std::size_t edits = 0;
};

// Sums over every admissible alignment of the read at this start -- it is a likelihood, not a best
// path -- and returns that sum decomposed by consumed haplotype span, so each part can be given the
// fragment geometry it actually implies. Returns nullopt when the haplotype cannot hold the read at
// this start at all.
std::optional<GappedFit> gapped_emission(
    const std::string& haplotype,
    const std::string& read,
    std::size_t start,
    const LibraryModel& library,
    std::size_t mismatch_band);

// The same emission, anchored at the read's RIGHT endpoint instead of at its start.
//
// An FR pair is sequenced from both ends of ONE template: the forward mate is read rightward from
// the left endpoint, the reverse mate leftward from the right endpoint. Each mate's emission has to
// be a proper distribution given ITS OWN anchor, and the reverse mate's anchor is the right
// endpoint -- exactly as FragmentPlacement above describes a physical origin.
//
// This is not gapped_emission with the arguments moved about. The recursion is DIRECTIONAL: an
// alignment may end in a gap but may not begin in one, so a path and its reversal do not carry the
// same mass (measured at 3.8 to 11.2 nats apart), and the normaliser is conditioned on the sequence
// available in the direction the mate is actually read.
//
// Left-anchoring the reverse mate is what broke the exposure identity. Several starts share one
// right endpoint -- rs = R - consumed, one per consumed span -- and each was normalised against a
// different amount of remaining haplotype, so gathering one consumed-span class from each summed to
// MORE than one wherever the haplotype end bound the normaliser. The excess appeared only near the
// right boundary and walked inward as max_gap_bp widened, which is precisely what was measured.
// Anchored here, the identity holds exactly at every (mismatch_band, max_gap_bp) combination, with
// no precondition relating the two, and at max_gap_bp 0 it is identical to gapped_emission.
//
// `end` is one past the last haplotype base the mate may consume; a component's consumed span is
// measured leftward from it, so that alignment occupies [end - consumed_bp, end).
std::optional<GappedFit> gapped_emission_ending_at(
    const std::string& haplotype,
    const std::string& read,
    std::size_t end,
    const LibraryModel& library,
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

// The same admission rule stated on the template's two ENDPOINTS, which is how an origin is named
// once the reverse mate is anchored at its right endpoint. Equivalent to fr_insert_length with
// reverse_start = reverse_end - reverse_consumed_bp, but it cannot underflow when a gapped mate
// consumes more than its endpoint offset, and it does not ask the caller to reconstruct a start it
// no longer has.
std::optional<std::size_t> fr_insert_from_endpoints(
    std::size_t forward_start,
    std::size_t forward_consumed_bp,
    std::size_t reverse_end,
    std::size_t reverse_consumed_bp,
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
    const LibraryModel& library,
    std::size_t mismatch_band);

// Summed log mass of a fragment on one chromosome copy under the truncated model: sum over physical
// origins of (1/2) * P_library(insert) * the origin's summed latent alignment mass. `insert` must be
// the library distribution. Returns -inf when no placement is accepted.
//
// The orientation prior and the insert probability multiply the origin ONCE, because one origin is
// one template with one insert length. The latent alignments inside it were already summed when it
// was built, at the geometry each of them implied.
// There is deliberately no band parameter. Under the substitution-only model the emission was
// recoverable from a mismatch count, so the same placement list could be re-normalised at scoring
// time by passing a different band. A gapped alignment is not summarised by a count: the band selects
// AND normalises the placements when they are enumerated, and this can only sum what they carry.
// Keeping a band argument here would let a caller ask for a renormalisation that silently does not
// happen, which is worse than not offering it.
double log_placement_mass(
    const std::vector<FragmentPlacement>& placements,
    const ReadPair& fragment,
    const DiscreteInsertDistribution& insert,
    const LibraryModel& library);

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
    //
    // Under parallel placement these are PROCESS-WIDE budgets, not per worker. Two mechanisms, with
    // different jobs, and the distinction matters:
    //
    //   * In flight, max_piece_hits and max_physical_origins are shared across workers and bound the
    //     WORK actually performed. Crossing one cancels every worker promptly rather than letting each
    //     consume a full cap -- which is what a circuit breaker is for, and what a per-worker guard
    //     silently gave up. physical_origins is an exact shared counter; piece_hits is published in
    //     batches of piece_hit_publish_batch because it is the hottest counter in the program (320M
    //     increments on C4 at audit 256), and batching was chosen to avoid unmeasured contention on a
    //     shared cache line. Whether an exact atomic would in fact be affordable here has NOT been
    //     measured. The overshoot that batching admits is bounded and declared: at most
    //     threads * piece_hit_publish_batch beyond the cap, which is 8,192 against a 2,000,000,000
    //     default.
    //
    //   * ACCEPTANCE is decided once, after the counters are merged, on the exact totals. A REFUSAL may
    //     also be raised earlier, by whichever worker observes the shared budget cross the cap -- so it
    //     is not true that every verdict waits for the merge. Both are sound, and the accept-or-refuse
    //     answer still does not depend on the thread count: a published total never exceeds the exact
    //     total, so a run whose exact total is within the cap cannot trip an early refusal, and a run
    //     that trips one would have exceeded the cap on the merged totals as well. A total equal to the
    //     cap accepts, one above refuses, however the work was partitioned.
    //
    // max_state_classes is deliberately only a post-merge acceptance cap. It is counted after a
    // fragment's states already exist, so bounding it in flight would limit neither memory nor work.
    std::uint64_t max_piece_hits = 2000000000ULL;
    // Increments a worker accumulates locally before publishing to the shared piece-hit budget. Larger
    // is cheaper and admits proportionally more overshoot; the bound is threads * this value.
    std::uint64_t piece_hit_publish_batch = 4096;
    std::uint64_t max_proposed_starts_per_haplotype = 50000000ULL;
    std::uint64_t max_physical_origins = 500000000ULL;
    std::uint64_t max_state_classes = 200000000ULL;
    std::uint64_t max_mass_cells = 200000000ULL;
    // Admitting indels multiplies the pairing work by up to (2*max_gap_bp+1)^2 per physical origin,
    // so it gets a declared bound of its own rather than being absorbed into the origin cap. Like
    // max_state_classes this is an acceptance cap checked on the merged totals: the pairs are counted
    // as the geometry is tested, and a running per-worker crossing point would move with the
    // partition. The default holds the same headroom over C4 that the origin cap does.
    std::uint64_t max_latent_alignment_pairs = 2000000000ULL;
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
    // Latent alignment pairs TESTED: one per (forward consumed span, reverse consumed span) offered
    // to the FR geometry, whether or not it was accepted. This is the pairing work that admitting
    // indels multiplies, which is why it counts attempts rather than successes -- `valid_fr_states`
    // is the accepted subset. Measured rather than assumed negligible: a mate has at most
    // 2*max_gap_bp+1 terminal spans, so one mate pair can present up to (2*max_gap_bp+1)^2 of them --
    // 49 at the production max_gap_bp of 3 -- before the collapse onto endpoints reduces them.
    std::uint64_t latent_alignment_pairs = 0;
    // Fragments whose declared insert support lets the geometry's floor -- insert >= the larger
    // consumed span -- actually bind. Only those need their consumed-span components carried
    // separately through the pairing: for every other fragment the mass at an anchor pair
    // factorises exactly as (sum over forward spans) * (sum over reverse spans), so one
    // whole-emission row per anchor yields the IDENTICAL likelihood at a fraction of the cost.
    // Counted because a fast path nobody can see is one that can diverge unnoticed; the scorer's
    // agreement with the always-componentised oracle is what certifies it.
    std::uint64_t fragments_with_floor_binding = 0;
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
    // Workers that actually ran the placement and mass phase. One when a single worker was asked for,
    // when there are fewer haplotypes than workers, and when full-origin recording is requested --
    // that diagnostic writes a shared map from the innermost loop and deliberately stays serial rather
    // than being made deterministic for a path nothing statistical depends on.
    std::size_t threads_used_placement = 1;
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
// What a refused run actually did, as it stood when a cap was reached. A refusal that reports nothing
// is unfalsifiable -- a breaker that stopped at its cap cannot be told from one that stopped long past
// it -- and reporting only counters loses the rest: which workers ran, and where the time went. Both
// were previously left at their success-path defaults on a refusal, so a refused run claimed one
// placement thread however many it used, and attributed its whole duration to no phase at all.
struct FragmentScoringObservation {
    FragmentScoringCounters counters;
    // Workers that actually ran placement, not the number requested.
    std::size_t threads_used_placement = 1;
    // Zero when the phase never ran, which is the usual case for a placement refusal. Deliberately not
    // defaulted to one: claiming a worker for a phase that did not happen is the defect this replaces.
    std::size_t threads_used_combination = 0;
    // The phases up to the refusal. Preparation is complete by then, placement is however far it got,
    // and combination has usually not started.
    double seconds_fragment_preparation = 0.0;
    double seconds_placement_and_mass = 0.0;
    double seconds_diploid_combination = 0.0;
};

// A measurement, not an optimisation. Oriented read forms repeat -- four forms per fragment, and at
// depth many fragments carry identical sequence -- so their placements could in principle be computed
// once per distinct sequence and reused. Whether that is worth building depends on the saving in PIECE
// HITS, not in forms: duplicate forms are plausibly the repetitive ones that contribute
// disproportionately many hits, so the form-level collapse is not the answer. This reports both, and
// the exact weighted prediction, so the decision rests on the quantity that actually costs.
struct PieceHitCensus {
    std::size_t total_forms = 0;
    std::size_t unique_forms = 0;
    // Seed entries across every piece-length bucket, and how many survive a collapse to one form per
    // distinct sequence.
    std::uint64_t total_entries = 0;
    std::uint64_t unique_entries = 0;
    // The scorer's own scan, counted rather than placed: sum over scanned positions of the entries
    // matching that seed. `unique` is the same sum counting only the first form of each distinct
    // sequence, which is exactly what a collapse would leave. Their ratio is the predicted saving.
    std::uint64_t predicted_piece_hits_all = 0;
    std::uint64_t predicted_piece_hits_unique = 0;
    // Distribution, so a mean saving cannot hide a single pathological seed or haplotype.
    std::uint64_t max_piece_hits_at_one_seed = 0;
    std::uint64_t max_piece_hits_one_haplotype = 0;
};

// Builds the same forms and piece indexes the scorer builds and scans the same haplotypes, counting
// hits instead of placing them. Changes nothing and is never called from scoring.
PieceHitCensus census_piece_hits(
    const InternedHaplotypes& haplotypes,
    const std::vector<ReadPair>& fragments,
    const FragmentScoringOptions& options);

// `observed_on_failure`, when given, receives the work performed when a cap was reached, so a refused
// run can report what it actually did rather than nothing at all.
FragmentScoringResult score_diplotypes(
    const InternedHaplotypes& haplotypes,
    const std::vector<HaplotypeCandidate>& candidates,
    const std::vector<ReadPair>& fragments,
    const LibraryModel& library,
    const FragmentScoringOptions& options,
    FragmentScoringObservation* observed_on_failure = nullptr);

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
