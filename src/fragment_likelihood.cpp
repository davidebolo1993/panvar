#include "panvar/fragment_likelihood.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <exception>
#include <map>
#include <thread>
#include <set>
#include <tuple>
#include <type_traits>
#include <cmath>
#include <stdexcept>

namespace panvar {
namespace {

constexpr double kNegInf = -std::numeric_limits<double>::infinity();

bool is_acgt(char base) noexcept {
    return base == 'A' || base == 'C' || base == 'G' || base == 'T';
}

char complement(char base) {
    switch (base) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        default:
            throw std::invalid_argument(
                "fragment likelihood is substitution-only over ACGT; a non-ACGT base must refuse "
                "rather than be scored as a mismatch");
    }
}

std::string reverse_complement_read(const std::string& bases);

std::string reverse_complement_read(const std::string& bases) {
    std::string out(bases.size(), 'N');
    for (std::size_t i = 0; i < bases.size(); ++i) {
        out[bases.size() - 1 - i] = complement(bases[i]);
    }
    return out;
}

} // namespace

std::optional<std::size_t> count_mismatches(
    const std::string& haplotype,
    const std::string& read,
    std::size_t start,
    std::size_t mismatch_limit) {
    if (read.empty()) {
        throw std::invalid_argument("fragment likelihood requires non-empty reads");
    }
    if (start > haplotype.size() || haplotype.size() - start < read.size()) {
        return std::nullopt;
    }
    std::size_t mismatches = 0;
    for (std::size_t i = 0; i < read.size(); ++i) {
        const char reference = haplotype[start + i];
        const char observed = read[i];
        if (!is_acgt(reference) || !is_acgt(observed)) {
            throw std::invalid_argument(
                "fragment likelihood is substitution-only over ACGT; a non-ACGT base must refuse "
                "rather than be scored as a mismatch");
        }
        if (reference != observed && ++mismatches > mismatch_limit) return std::nullopt;
    }
    return mismatches;
}

double log_emission(
    std::size_t read_length,
    std::size_t mismatches,
    double substitution_error) {
    if (mismatches > read_length) {
        throw std::invalid_argument("mismatch count exceeds the read length");
    }
    if (!(substitution_error > 0.0) || !(substitution_error < 0.75)) {
        throw std::invalid_argument("substitution error must be in (0, 0.75)");
    }
    // A substitution lands on one of three other bases, so a specific observed mismatch carries
    // error/3. Folding the /3 in here rather than at the call site keeps the oracle and the
    // accelerated path on the same constant.
    const double matches = static_cast<double>(read_length - mismatches);
    return matches * std::log(1.0 - substitution_error) +
           static_cast<double>(mismatches) * std::log(substitution_error / 3.0);
}

double log_band_mass(std::size_t read_length, std::size_t mismatch_band, double substitution_error) {
    if (!(substitution_error > 0.0) || !(substitution_error < 0.75)) {
        throw std::invalid_argument("substitution error must be in (0, 0.75)");
    }
    if (read_length == 0) throw std::invalid_argument("band mass needs a non-empty read");
    const std::size_t top = std::min(mismatch_band, read_length);
    const double n = static_cast<double>(read_length);
    double total = kNegInf;
    for (std::size_t m = 0; m <= top; ++m) {
        const double k = static_cast<double>(m);
        const double log_choose = std::lgamma(n + 1.0) - std::lgamma(k + 1.0) - std::lgamma(n - k + 1.0);
        total = log_add_exp(total, log_choose + k * std::log(substitution_error) +
                                       (n - k) * std::log1p(-substitution_error));
    }
    return total;
}

double log_truncated_emission(
    std::size_t read_length,
    std::size_t mismatches,
    double substitution_error,
    std::size_t mismatch_band) {
    if (mismatches > mismatch_band) return kNegInf;
    return log_emission(read_length, mismatches, substitution_error) -
           log_band_mass(read_length, mismatch_band, substitution_error);
}

std::optional<std::size_t> fr_insert_length(
    std::size_t forward_start,
    std::size_t forward_length,
    std::size_t reverse_start,
    std::size_t reverse_length,
    const InsertSupport& support) {
    // The reverse mate's own right edge closes the template. A reverse mate that ends before the
    // forward mate begins is not an FR pair at all.
    const std::size_t reverse_end = reverse_start + reverse_length;
    if (reverse_end <= forward_start) return std::nullopt;
    const std::size_t insert = reverse_end - forward_start;
    // Overlapping mates are a physical state, so the floor is the longer mate, not their sum. A
    // run-wide maximum read length here would wrongly raise the floor for shorter fragments.
    if (insert < std::max(forward_length, reverse_length)) return std::nullopt;
    if (support.empty() || insert < support.min_bp || insert > support.max_bp) return std::nullopt;
    return insert;
}

double haplotype_exposure(
    std::size_t haplotype_length,
    const DiscreteInsertDistribution& insert,
    std::size_t physical_floor_bp) {
    if (insert.support.empty()) return 0.0;
    double exposure = 0.0;
    const std::size_t from = std::max(insert.support.min_bp, physical_floor_bp);
    for (std::size_t length = from; length <= insert.support.max_bp; ++length) {
        if (haplotype_length < length) continue;
        const double probability = std::exp(insert.log_probability(length));
        exposure += probability * static_cast<double>(haplotype_length - length + 1);
    }
    return exposure;
}

double log_fragment_contribution(
    double log_mass,
    double log_background,
    const LibraryModel& library) {
    validate_library_model(library);
    const double signal = log_mass == kNegInf
        ? kNegInf
        : std::log1p(-library.outlier_weight) +
              std::log(library.fragment_start_intensity) + log_mass;
    if (library.outlier_weight == 0.0) {
        if (signal == kNegInf) {
            throw std::invalid_argument(
                "a fragment with no accepted placement and no outlier mixture makes every candidate "
                "impossible; declare a positive outlier weight or refuse the fragment upstream");
        }
        return signal;
    }
    return log_add_exp(signal, std::log(library.outlier_weight) + log_background);
}

double expected_signal_fragments(double exposure_sum, const LibraryModel& library) {
    validate_library_model(library);
    return (1.0 - library.outlier_weight) * library.fragment_start_intensity * exposure_sum;
}

std::vector<FragmentPlacement> enumerate_placements_exhaustive(
    const std::string& haplotype,
    const ReadPair& fragment,
    const InsertSupport& support,
    std::size_t mismatch_limit) {
    validate_read_pair(fragment);
    std::vector<FragmentPlacement> placements;
    if (haplotype.empty() || support.empty()) return placements;

    const std::string first_rc = reverse_complement_read(fragment.first.bases);
    const std::string second_rc = reverse_complement_read(fragment.second.bases);

    for (int arm = 0; arm < 2; ++arm) {
        const MateOrientation orientation =
            arm == 0 ? MateOrientation::FirstForward : MateOrientation::SecondForward;
        // Whichever mate reads off the forward strand is compared as given; the other is compared
        // pre-complemented, which makes both comparisons forward.
        const std::string& forward_bases =
            arm == 0 ? fragment.first.bases : fragment.second.bases;
        const std::string& reverse_bases = arm == 0 ? second_rc : first_rc;

        for (std::size_t forward_start = 0; forward_start < haplotype.size(); ++forward_start) {
            const std::optional<std::size_t> forward_mismatches =
                count_mismatches(haplotype, forward_bases, forward_start, mismatch_limit);
            if (!forward_mismatches.has_value()) continue;
            for (std::size_t reverse_start = 0; reverse_start < haplotype.size(); ++reverse_start) {
                const std::optional<std::size_t> insert = fr_insert_length(
                    forward_start, forward_bases.size(),
                    reverse_start, reverse_bases.size(), support);
                if (!insert.has_value()) continue;
                const std::optional<std::size_t> reverse_mismatches =
                    count_mismatches(haplotype, reverse_bases, reverse_start, mismatch_limit);
                if (!reverse_mismatches.has_value()) continue;

                FragmentPlacement placement;
                placement.orientation = orientation;
                placement.forward_start = forward_start;
                placement.reverse_start = reverse_start;
                placement.insert_bp = *insert;
                placement.forward_mismatches = *forward_mismatches;
                placement.reverse_mismatches = *reverse_mismatches;
                placements.push_back(placement);
            }
        }
    }
    return placements;
}

double log_placement_mass(
    const std::vector<FragmentPlacement>& placements,
    const ReadPair& fragment,
    const DiscreteInsertDistribution& insert,
    const LibraryModel& library,
    std::size_t mismatch_band) {
    validate_read_pair(fragment);
    validate_library_model(library);
    double mass = kNegInf;
    for (const FragmentPlacement& placement : placements) {
        const std::size_t forward_length = placement.orientation == MateOrientation::FirstForward
            ? fragment.first.bases.size() : fragment.second.bases.size();
        const std::size_t reverse_length = placement.orientation == MateOrientation::FirstForward
            ? fragment.second.bases.size() : fragment.first.bases.size();
        const double term =
            kLogMateOrientationPrior +
            insert.log_probability(placement.insert_bp) +
            log_truncated_emission(forward_length, placement.forward_mismatches,
                                   library.substitution_error, mismatch_band) +
            log_truncated_emission(reverse_length, placement.reverse_mismatches,
                                   library.substitution_error, mismatch_band);
        mass = log_add_exp(mass, term);
    }
    return mass;
}


// -------------------------------------------------------------------------------------------
// Accelerated scorer.
// -------------------------------------------------------------------------------------------

namespace {

constexpr std::uint64_t kNoCode = std::numeric_limits<std::uint64_t>::max();

std::uint64_t two_bit(char base) {
    switch (base) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return kNoCode;
    }
}

// One oriented copy of one mate, already turned so that it matches the haplotype's forward strand.
struct ReadForm {
    std::size_t fragment = 0;
    std::string bases;
    bool is_forward_mate_of_first_arm = false;   // bookkeeping only; arms are resolved by role below
};

// piece code -> (form, offset of the piece inside that form)
struct PieceEntry {
    std::uint64_t code = 0;
    std::uint32_t form = 0;
    std::uint32_t offset = 0;
};

struct PieceIndex {
    std::size_t piece_bp = 0;
    std::vector<PieceEntry> entries;   // sorted by code
};

// A verified placement of one read form on one haplotype.
struct FormPlacement {
    std::uint32_t start = 0;
    std::uint32_t mismatches = 0;
};

// The exact statistical class of a physical origin. Origins sharing all four fields may be stored
// once with a multiplicity, but only because every term of their mass is identical.
struct StateKey {
    std::uint32_t insert_bp = 0;
    std::uint32_t forward_mismatches = 0;
    std::uint32_t reverse_mismatches = 0;
    std::uint32_t orientation = 0;

    bool operator<(const StateKey& other) const noexcept {
        return std::tie(insert_bp, forward_mismatches, reverse_mismatches, orientation) <
               std::tie(other.insert_bp, other.forward_mismatches, other.reverse_mismatches,
                        other.orientation);
    }
};

std::vector<PieceIndex> build_piece_indexes(
    const std::vector<ReadForm>& forms,
    std::size_t band,
    std::size_t minimum_piece_bp) {
    // Reads of different lengths yield pieces of different lengths, and one rolling scan can only
    // carry one length. Bucket by piece length and scan once per bucket.
    std::map<std::size_t, std::vector<PieceEntry>> by_length;
    for (std::size_t f = 0; f < forms.size(); ++f) {
        const std::size_t length = forms[f].bases.size();
        // The band guarantees some contiguous 1/(band+1) chunk is exact. Keying on a PREFIX of that
        // chunk keeps the guarantee, because an exact chunk has an exact prefix, and caps the key at
        // what a 2-bit code can hold.
        const std::size_t chunk = length / (band + 1);
        const std::size_t piece = std::min<std::size_t>(31, chunk);
        if (piece < minimum_piece_bp) {
            throw std::invalid_argument(
                "the declared mismatch band splits a read into pieces too short to index; widen the "
                "reads, narrow the band, or lower minimum_piece_bp deliberately");
        }
        for (std::size_t p = 0; p <= band; ++p) {
            const std::size_t offset = p * chunk;
            std::uint64_t code = 0;
            bool usable = true;
            for (std::size_t i = 0; i < piece; ++i) {
                const std::uint64_t bits = two_bit(forms[f].bases[offset + i]);
                if (bits == kNoCode) { usable = false; break; }
                code = (code << 2) | bits;
            }
            if (!usable) continue;
            by_length[piece].push_back({code, static_cast<std::uint32_t>(f),
                                        static_cast<std::uint32_t>(offset)});
        }
    }
    std::vector<PieceIndex> indexes;
    for (auto& item : by_length) {
        std::sort(item.second.begin(), item.second.end(),
                  [](const PieceEntry& lhs, const PieceEntry& rhs) {
                      return std::tie(lhs.code, lhs.form, lhs.offset) <
                             std::tie(rhs.code, rhs.form, rhs.offset);
                  });
        indexes.push_back({item.first, std::move(item.second)});
    }
    return indexes;
}

// Process-wide budget for the two caps that bound WORK rather than merely decide acceptance. Shared by
// every placement worker, so T workers cannot each consume a full cap before anyone refuses -- which is
// what a per-worker guard silently allowed, and what makes a 64-thread request nearly breaker-free.
//
// Crossing a cap here DOES decide the answer: the worker that observes it refuses immediately, which is
// sound because a shared total can only undercount the exact one. Acceptance is the half decided later,
// on the exact merged totals. `cancelled` is what makes the refusal prompt: the worker that crosses a
// cap raises it, and every other worker abandons its range at the next check -- in the piece-hit
// publication, at each physical origin, and at the haplotype, fragment, sort and verification
// boundaries -- rather than running to completion.
struct PlacementBudget {
    std::atomic<std::uint64_t> piece_hits{0};
    std::atomic<std::uint64_t> physical_origins{0};
    std::atomic<bool> cancelled{false};
    std::uint64_t max_piece_hits = 0;
    std::uint64_t max_physical_origins = 0;
    std::uint64_t publish_batch = 4096;

    bool stop() const { return cancelled.load(std::memory_order_relaxed); }
    void cancel() { cancelled.store(true, std::memory_order_relaxed); }
};

// Raised by a worker abandoning its range because another worker already crossed a cap. It never
// escapes score_diplotypes: the cap that was actually crossed is reported instead, so the message a
// caller sees does not depend on which worker noticed first.
class PlacementCancelled : public std::exception {};

// Scan one haplotype once per piece-length bucket, turning every piece hit into an implied read
// start. This is the inversion the design calls for: the haplotype is walked, not queried once per
// read against a position hash the size of the locus.
void collect_form_placements(
    const std::string& haplotype,
    const std::vector<ReadForm>& forms,
    const std::vector<PieceIndex>& indexes,
    PlacementBudget& budget,
    std::size_t band,
    std::vector<std::vector<FormPlacement>>& out,
    FragmentScoringCounters& counters,
    const FragmentScoringOptions& options) {
    std::vector<std::pair<std::uint32_t, std::uint32_t>> proposed;   // (form, start)
    // Piece hits accumulated since the last publication to the shared budget. Drained before this
    // function returns, so nothing a worker counted goes unpublished and the merged total stays exact.
    std::uint64_t unpublished_piece_hits = 0;
    const auto publish_piece_hits = [&]() {
        if (unpublished_piece_hits == 0) return;
        const std::uint64_t total =
            budget.piece_hits.fetch_add(unpublished_piece_hits, std::memory_order_relaxed) +
            unpublished_piece_hits;
        unpublished_piece_hits = 0;
        if (total > budget.max_piece_hits) {
            budget.cancel();
            throw ResourceExhausted(
                "fragment scoring reached its piece-hit cap; the locus is too repetitive for the "
                "declared caps and no partial ranking is returned");
        }
    };

    for (const PieceIndex& index : indexes) {
        const std::size_t piece = index.piece_bp;
        if (haplotype.size() < piece) continue;
        counters.candidate_bases_scanned += haplotype.size();

        const std::uint64_t mask = piece == 32
            ? std::numeric_limits<std::uint64_t>::max()
            : ((std::uint64_t{1} << (2 * piece)) - 1);
        std::uint64_t code = 0;
        std::size_t run = 0;
        for (std::size_t at = 0; at < haplotype.size(); ++at) {
            const std::uint64_t bits = two_bit(haplotype[at]);
            if (bits == kNoCode) { run = 0; code = 0; continue; }
            code = ((code << 2) | bits) & mask;
            if (++run < piece) continue;
            const std::size_t piece_start = at + 1 - piece;

            const auto range = std::equal_range(
                index.entries.begin(), index.entries.end(), code,
                [](const auto& lhs, const auto& rhs) {
                    if constexpr (std::is_same_v<std::decay_t<decltype(lhs)>, PieceEntry>) {
                        return lhs.code < rhs;
                    } else {
                        return lhs < rhs.code;
                    }
                });
            for (auto it = range.first; it != range.second; ++it) {
                // The hottest counter in the program -- 320M increments on C4 at audit 256 -- so it is
                // accumulated locally and published to the shared budget in batches, to avoid unmeasured
                // contention on a shared cache line; whether an exact atomic would be affordable here
                // has not been measured. The overshoot batching admits is bounded by
                // threads * publish_batch and declared in the header. Refusing here is sound the moment
                // the published global total exceeds the cap, because a published total can only
                // undercount the exact one; acceptance is still decided later, on the merged totals.
                ++counters.piece_hits;
                if (++unpublished_piece_hits >= budget.publish_batch) {
                    const std::uint64_t total =
                        budget.piece_hits.fetch_add(unpublished_piece_hits,
                                                    std::memory_order_relaxed) +
                        unpublished_piece_hits;
                    unpublished_piece_hits = 0;
                    if (total > budget.max_piece_hits) {
                        budget.cancel();
                        throw ResourceExhausted(
                            "fragment scoring reached its piece-hit cap; the locus is too repetitive "
                            "for the declared caps and no partial ranking is returned");
                    }
                    if (budget.stop()) throw PlacementCancelled();
                }
                if (it->offset > piece_start) continue;   // implied start would be negative
                proposed.emplace_back(it->form,
                                      static_cast<std::uint32_t>(piece_start - it->offset));
                ++counters.proposed_starts;
                if (proposed.size() > options.max_proposed_starts_per_haplotype) {
                    throw ResourceExhausted(
                        "fragment scoring reached its per-haplotype proposed-start cap; no partial "
                        "ranking is returned");
                }
            }
        }
    }

    // Drain whatever this haplotype counted below the last batch boundary, so the shared total is the
    // true total rather than a multiple of publish_batch, and the merged totals stay exact.
    publish_piece_hits();

    // Sorting and full-read verification are the rest of this haplotype's cost and neither touches the
    // budget, so without these checks a worker cancelled by someone else's refusal would still pay for
    // both before noticing. Prompt cancellation has to mean prompt in the expensive loops too.
    if (budget.stop()) throw PlacementCancelled();

    // Several pieces of one read can rediscover one physical placement. Deduplicate before
    // verifying, so a placement found three ways is verified once and counted once.
    std::sort(proposed.begin(), proposed.end());
    proposed.erase(std::unique(proposed.begin(), proposed.end()), proposed.end());
    counters.unique_starts += proposed.size();
    if (budget.stop()) throw PlacementCancelled();

    for (const auto& item : proposed) {
        // One relaxed load against a full-read comparison: free at this granularity, and verification is
        // where a cancelled worker would otherwise spend the longest.
        if (budget.stop()) throw PlacementCancelled();
        ++counters.full_read_verifications;
        const std::optional<std::size_t> mismatches =
            count_mismatches(haplotype, forms[item.first].bases, item.second, band);
        if (!mismatches.has_value()) continue;
        out[item.first].push_back({item.second, static_cast<std::uint32_t>(*mismatches)});
        ++counters.accepted_mate_placements;
    }
    for (std::vector<FormPlacement>& list : out) {
        std::sort(list.begin(), list.end(), [](const FormPlacement& lhs, const FormPlacement& rhs) {
            return lhs.start < rhs.start;
        });
    }
}

} // namespace

double FragmentScoringResult::mass_of(std::size_t fragment, std::size_t haplotype) const {
    if (unique_haplotypes == 0 || haplotype >= unique_haplotypes ||
        fragment * unique_haplotypes + haplotype >= log_mass.size()) {
        throw std::out_of_range("fragment/haplotype mass index is out of range");
    }
    return log_mass[fragment * unique_haplotypes + haplotype];
}

InternedHaplotypes intern_haplotypes(
    const std::vector<MosaicPath>& allele_paths,
    const std::vector<std::string>& sequences) {
    if (allele_paths.size() != sequences.size()) {
        throw std::invalid_argument("allele paths and spelled sequences differ in size");
    }
    InternedHaplotypes result;
    result.index_of.reserve(sequences.size());
    std::map<std::string, std::size_t> slot_by_sequence;
    std::set<std::vector<std::size_t>> distinct_paths;
    for (std::size_t i = 0; i < sequences.size(); ++i) {
        distinct_paths.insert(allele_paths[i].alleles);
        const auto found = slot_by_sequence.find(sequences[i]);
        if (found != slot_by_sequence.end()) {
            result.index_of.push_back(found->second);
            continue;
        }
        const std::size_t slot = result.sequences.size();
        slot_by_sequence.emplace(sequences[i], slot);
        result.sequences.push_back(sequences[i]);
        result.index_of.push_back(slot);
    }
    // Two different allele paths can spell identical bytes. They are one haplotype to any read-based
    // method, so they are scored once -- but the count of distinct paths is retained, because
    // collapsing them silently would hide that the candidate set is narrower than it looks.
    result.distinct_allele_paths = distinct_paths.size();
    return result;
}

FragmentScoringResult score_diplotypes(
    const InternedHaplotypes& haplotypes,
    const std::vector<HaplotypeCandidate>& candidates,
    const std::vector<ReadPair>& fragments,
    const LibraryModel& library,
    const FragmentScoringOptions& options,
    FragmentScoringObservation* observed_on_failure) {
    validate_library_model(library);
    if (options.minimum_piece_bp == 0) {
        throw std::invalid_argument("minimum indexable piece length must be positive");
    }
    if (options.max_piece_hits == 0 || options.max_proposed_starts_per_haplotype == 0 ||
        options.max_physical_origins == 0 || options.max_state_classes == 0 ||
        options.max_mass_cells == 0) {
        throw std::invalid_argument("every fragment-scoring resource cap must be positive");
    }
    for (const HaplotypeCandidate& candidate : candidates) {
        if (candidate.first >= haplotypes.sequences.size() ||
            candidate.second >= haplotypes.sequences.size()) {
            throw std::invalid_argument("diplotype candidate names an unknown haplotype");
        }
    }

    FragmentScoringResult result;
    result.scores.resize(candidates.size());
    result.unique_haplotypes = haplotypes.sequences.size();
    FragmentScoringCounters& counters = result.counters;
    const std::size_t band = options.mismatch_band;
    const DiscreteInsertDistribution library_insert = build_library_insert_distribution(library);

    // Filled completely before any worker starts, and read-only thereafter. A lazily-filled shared
    // cache would be a data race under placement parallelism, and giving each worker its own would
    // recompute the same handful of values per thread; the set of lengths it can ever be asked for is
    // exactly the distinct mate lengths of the fragments, which are known here.
    std::map<std::size_t, double> band_mass_by_length;
    for (const ReadPair& fragment : fragments) {
        for (const std::size_t length : {fragment.first.bases.size(), fragment.second.bases.size()}) {
            if (band_mass_by_length.find(length) == band_mass_by_length.end()) {
                band_mass_by_length.emplace(length, log_band_mass(length, band, library.substitution_error));
            }
        }
    }
    const auto band_mass = [&](std::size_t length) -> double {
        const auto found = band_mass_by_length.find(length);
        if (found != band_mass_by_length.end()) return found->second;
        // Only reachable if a length appears that no fragment mate carries, which cannot happen for
        // the calls below; computing it rather than asserting keeps the function total.
        return log_band_mass(length, band, library.substitution_error);
    };

    // Phase 1 of three: per-fragment preparation. Timed so a caller can attribute rescoring rather
    // than report it whole -- it is 57-93% of a run, and which phase dominates decides how it should
    // be parallelised.
    const auto preparation_mark = std::chrono::steady_clock::now();

    // Four oriented forms per fragment: each mate as given, and each pre-complemented. Whichever
    // mate reads off the forward strand is compared as given, so every comparison is forward.
    std::vector<ReadForm> forms;
    forms.reserve(fragments.size() * 4);
    std::vector<InsertSupport> support_by_fragment;
    std::vector<double> background_by_fragment;
    std::map<std::size_t, std::size_t> fragments_by_floor;
    support_by_fragment.reserve(fragments.size());
    background_by_fragment.reserve(fragments.size());
    for (std::size_t f = 0; f < fragments.size(); ++f) {
        validate_read_pair(fragments[f]);
        forms.push_back({f, fragments[f].first.bases, true});
        forms.push_back({f, reverse_complement_read(fragments[f].first.bases), false});
        forms.push_back({f, fragments[f].second.bases, false});
        forms.push_back({f, reverse_complement_read(fragments[f].second.bases), true});
        support_by_fragment.push_back(fragment_insert_support(fragments[f], library));
        background_by_fragment.push_back(log_background_probability(fragments[f]));
        ++fragments_by_floor[std::max(fragments[f].first.bases.size(),
                                      fragments[f].second.bases.size())];
    }
    const auto form_slot = [](std::size_t fragment, std::size_t which) {
        return fragment * 4 + which;   // 0: mate1, 1: rc(mate1), 2: mate2, 3: rc(mate2)
    };

    const std::vector<PieceIndex> indexes =
        forms.empty() ? std::vector<PieceIndex>{}
                      : build_piece_indexes(forms, band, options.minimum_piece_bp);

    result.seconds_fragment_preparation =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - preparation_mark).count();
    // Phase 2: placement and mass construction, per (fragment, haplotype).
    const auto placement_mark = std::chrono::steady_clock::now();

    const std::size_t unique = haplotypes.sequences.size();
    // Checked before the allocation it guards, so a refusal costs nothing.
    if (!fragments.empty() && unique > options.max_mass_cells / fragments.size()) {
        throw ResourceExhausted(
            "the fragment-by-haplotype mass matrix would exceed its declared cell cap; no partial "
            "ranking is returned");
    }
    std::vector<double> mass(fragments.size() * unique, kNegInf);
    std::vector<double> exposure(unique, 0.0);
    if (options.record_origin_summary) result.origin_summary.assign(fragments.size() * unique, OriginSummary{});
    std::vector<char> record_origins(fragments.size(), 0);
    for (const std::size_t f : options.record_origins_for_fragments) {
        if (f >= fragments.size()) throw std::invalid_argument("origin recording names an unknown fragment");
        record_origins[f] = 1;
    }

    // The fragment likelihood conditions on the observed read lengths, and a template shorter than a
    // fragment's longer mate cannot produce that fragment. Exposure therefore averages the analytic
    // count over the observed distribution of physical floors, so the -lambda*E term and the signal
    // density describe the same generated process.
    for (std::size_t h = 0; h < unique; ++h) {
        const std::size_t length = haplotypes.sequences[h].size();
        if (fragments.empty()) {
            exposure[h] = haplotype_exposure(length, library_insert, 0);
            continue;
        }
        for (const auto& item : fragments_by_floor) {
            exposure[h] += static_cast<double>(item.second) /
                           static_cast<double>(fragments.size()) *
                           haplotype_exposure(length, library_insert, item.first);
        }
    }

    // Placement and mass construction, partitioned over haplotypes. Every write below is indexed by
    // this worker's own h -- mass[f * unique + h], exposure[h], origin_summary[f * unique + h] -- so no
    // two workers touch the same cell, and `forms`, `indexes`, the per-fragment supports and the
    // sequences are read-only by the time we get here.
    //
    // Counters are per worker and merged in thread-index order, so the totals a run reports do not
    // depend on the schedule. The two running-total caps inside this phase -- physical_origins and
    // statistical_state_classes -- are therefore checked AFTER the merge rather than mid-loop: their
    // crossing point moves with the partition, and the declared contract is that a given input either
    // completes or refuses regardless of thread count. `max_proposed_starts_per_haplotype` is unaffected
    // and still throws locally, because it tests a per-haplotype vector rather than a running total.
    // One budget for every worker: these two caps bound the work performed, not merely the answer, so
    // they cannot be per worker without letting T workers each consume a full cap.
    PlacementBudget budget;
    budget.max_piece_hits = options.max_piece_hits;
    budget.max_physical_origins = options.max_physical_origins;
    budget.publish_batch = std::max<std::uint64_t>(1, options.piece_hit_publish_batch);

    const auto place_range = [&](std::size_t first_h, std::size_t last_h,
                                 FragmentScoringCounters& counters) {
    std::map<StateKey, std::uint64_t> states;
    for (std::size_t h = first_h; h < last_h; ++h) {
        // Haplotype and fragment boundaries are the natural places to notice someone else's refusal:
        // both are cheap relative to the work they precede, and a worker that has just finished one
        // haplotype would otherwise start another that can never be used.
        if (budget.stop()) throw PlacementCancelled();
        const std::string& haplotype = haplotypes.sequences[h];
        if (forms.empty()) continue;

        std::vector<std::vector<FormPlacement>> placements(forms.size());
        collect_form_placements(haplotype, forms, indexes, budget, band, placements, counters, options);

        for (std::size_t f = 0; f < fragments.size(); ++f) {
            if (budget.stop()) throw PlacementCancelled();
            const std::size_t first_length = fragments[f].first.bases.size();
            const std::size_t second_length = fragments[f].second.bases.size();
            const InsertSupport& support = support_by_fragment[f];
            if (support.empty()) continue;
            states.clear();

            for (int arm = 0; arm < 2; ++arm) {
                // FirstForward: mate one as given is the forward mate, rc(mate two) the reverse one.
                const std::vector<FormPlacement>& forward =
                    placements[form_slot(f, arm == 0 ? 0 : 2)];
                const std::vector<FormPlacement>& reverse =
                    placements[form_slot(f, arm == 0 ? 3 : 1)];
                const std::size_t forward_length = arm == 0 ? first_length : second_length;
                const std::size_t reverse_length = arm == 0 ? second_length : first_length;
                if (forward.empty() || reverse.empty()) continue;

                for (const FormPlacement& start : forward) {
                    // reverse_start + reverse_length - forward_start must lie in the support, so the
                    // reverse starts form one contiguous window. Two-pointer bounds, no map.
                    const std::size_t lo_end = static_cast<std::size_t>(start.start) +
                                               std::max(support.min_bp, std::max(forward_length,
                                                                                 reverse_length));
                    const std::size_t hi_end =
                        static_cast<std::size_t>(start.start) + support.max_bp;
                    if (hi_end < reverse_length) continue;
                    const std::size_t lo_start =
                        lo_end < reverse_length ? 0 : lo_end - reverse_length;
                    const std::size_t hi_start = hi_end - reverse_length;

                    const auto begin = std::lower_bound(
                        reverse.begin(), reverse.end(), lo_start,
                        [](const FormPlacement& lhs, std::size_t value) {
                            return lhs.start < value;
                        });
                    const auto end = std::upper_bound(
                        begin, reverse.end(), hi_start,
                        [](std::size_t value, const FormPlacement& rhs) {
                            return value < rhs.start;
                        });
                    for (auto it = begin; it != end; ++it) {
                        const std::optional<std::size_t> insert = fr_insert_length(
                            start.start, forward_length, it->start, reverse_length, support);
                        if (!insert.has_value()) continue;
                        ++counters.valid_fr_states;
                        // Thirty times cooler than piece_hits (9.4M against 320M on C4), so this one
                        // carries an exact shared counter: the work bound is tight, with no batching
                        // overshoot at all. Crossing it cancels every worker rather than letting each
                        // consume a full cap. Refusing here is exact and immediate; acceptance is still
                        // decided later, on the merged totals.
                        ++counters.physical_origins;
                        const std::uint64_t origins =
                            budget.physical_origins.fetch_add(1, std::memory_order_relaxed) + 1;
                        if (origins > budget.max_physical_origins) {
                            budget.cancel();
                            throw ResourceExhausted(
                                "fragment scoring reached its physical-origin cap; no partial "
                                "ranking is returned");
                        }
                        if (budget.stop()) throw PlacementCancelled();
                        StateKey key;
                        key.insert_bp = static_cast<std::uint32_t>(*insert);
                        key.forward_mismatches = start.mismatches;
                        key.reverse_mismatches = it->mismatches;
                        key.orientation = static_cast<std::uint32_t>(arm);
                        ++states[key];
                        if (options.record_origin_summary || record_origins[f] != 0) {
                            FragmentPlacement placement;
                            placement.forward_start = start.start;
                            placement.reverse_start = it->start;
                            placement.insert_bp = *insert;
                            placement.forward_mismatches = start.mismatches;
                            placement.reverse_mismatches = it->mismatches;
                            placement.orientation =
                                arm == 0 ? MateOrientation::FirstForward : MateOrientation::SecondForward;
                            if (options.record_origin_summary) {
                                OriginSummary& summary = result.origin_summary[f * unique + h];
                                ++summary.origins;
                                const double term =
                                    kLogMateOrientationPrior + library_insert.log_probability(*insert) +
                                    log_emission(forward_length, start.mismatches, library.substitution_error) -
                                    band_mass(forward_length) +
                                    log_emission(reverse_length, it->mismatches, library.substitution_error) -
                                    band_mass(reverse_length);
                                if (term > summary.best_log_term) {
                                    summary.best_log_term = term;
                                    summary.best = placement;
                                }
                            }
                            if (record_origins[f] != 0) result.recorded_origins[{f, h}].push_back(placement);
                        }
                    }
                }
            }
            if (states.empty()) continue;
            // A pure running total with no memory role, so it is accumulated here and decided once,
            // after the ordered merge, where the sum is the same at any thread count.
            counters.statistical_state_classes += states.size();

            double total = kNegInf;
            for (const auto& item : states) {
                const std::size_t forward_length =
                    item.first.orientation == 0 ? first_length : second_length;
                const std::size_t reverse_length =
                    item.first.orientation == 0 ? second_length : first_length;
                // Exact log-multiplicity: every origin in this class has an identical mass term, so
                // multiplying by the count is not an approximation of summing them, it IS the sum.
                const double term =
                    std::log(static_cast<double>(item.second)) + kLogMateOrientationPrior +
                    library_insert.log_probability(item.first.insert_bp) +
                    log_emission(forward_length, item.first.forward_mismatches,
                                 library.substitution_error) - band_mass(forward_length) +
                    log_emission(reverse_length, item.first.reverse_mismatches,
                                 library.substitution_error) - band_mass(reverse_length);
                total = log_add_exp(total, term);
            }
            mass[f * unique + h] = total;
        }
    }
    };

    // Full-origin recording writes result.recorded_origins, a shared map, from the innermost loop.
    // Rather than make that deterministic for a diagnostic no score depends on, it forces the serial
    // path -- and says so in threads_used_placement instead of claiming workers it did not use.
    const bool origins_force_serial = !options.record_origins_for_fragments.empty();
    const std::size_t placement_threads =
        (unique == 0 || origins_force_serial)
            ? 1
            : std::max<std::size_t>(1, std::min(options.threads, unique));
    result.threads_used_placement = placement_threads;

    // A refusal reports the work it performed, so the breaker can be held to a bound and not merely to
    // a verdict. After a worker threw, `counters` is still empty -- the merge below never ran -- and the
    // budget's published totals are the only global record; after a post-merge check `counters` is exact
    // and includes the remainders no batch published. The larger of the two is right in both cases. On
    // the worker path the other fields stay zero: they are per-range sums that were never merged.
    const auto note_observed_work = [&]() {
        if (observed_on_failure == nullptr) return;
        observed_on_failure->counters = counters;
        observed_on_failure->counters.piece_hits =
            std::max(counters.piece_hits, budget.piece_hits.load(std::memory_order_relaxed));
        observed_on_failure->counters.physical_origins = std::max(
            counters.physical_origins, budget.physical_origins.load(std::memory_order_relaxed));
        // The workers that actually ran, and the phases as they stand at the refusal: preparation is
        // complete, placement is however far it got, combination has not started. Reporting the
        // placement elapsed rather than zero is the point -- a refused run's time is not
        // unattributable, and it is exactly when knowing where the time went matters most.
        observed_on_failure->threads_used_placement = placement_threads;
        observed_on_failure->threads_used_combination = 0;
        observed_on_failure->seconds_fragment_preparation = result.seconds_fragment_preparation;
        observed_on_failure->seconds_placement_and_mass =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - placement_mark).count();
        observed_on_failure->seconds_diploid_combination = 0.0;
    };

    if (placement_threads <= 1) {
        try {
            place_range(0, unique, counters);
        } catch (const ResourceExhausted&) {
            note_observed_work();
            throw;
        }
    } else {
        const std::size_t span = (unique + placement_threads - 1) / placement_threads;
        std::vector<FragmentScoringCounters> local(placement_threads);
        std::vector<std::exception_ptr> failures(placement_threads);
        std::vector<std::thread> workers;
        workers.reserve(placement_threads);
        for (std::size_t t = 0; t < placement_threads; ++t) {
            const std::size_t first_h = std::min(t * span, unique);
            const std::size_t last_h = std::min(first_h + span, unique);
            workers.emplace_back([&, t, first_h, last_h]() {
                try {
                    place_range(first_h, last_h, local[t]);
                } catch (const PlacementCancelled&) {
                    // Another worker crossed a cap and this one abandoned its range. That is the
                    // breaker working, not a failure of its own: the worker that actually crossed
                    // reports the cap, so the message does not depend on who noticed first.
                } catch (...) {
                    failures[t] = std::current_exception();
                }
            });
        }
        for (std::thread& worker : workers) worker.join();
        // Rethrown in thread-index order, so which worker's message surfaces does not depend on which
        // one crossed the cap first. A refusal records the work performed before it stopped.
        for (const std::exception_ptr& failure : failures) {
            if (!failure) continue;
            try {
                std::rethrow_exception(failure);
            } catch (const ResourceExhausted&) {
                note_observed_work();
                throw;
            }
        }
        // Summed in thread-index order, so the totals a run reports are the serial totals exactly.
        for (const FragmentScoringCounters& part : local) {
            counters.candidate_bases_scanned += part.candidate_bases_scanned;
            counters.piece_hits += part.piece_hits;
            counters.proposed_starts += part.proposed_starts;
            counters.unique_starts += part.unique_starts;
            counters.full_read_verifications += part.full_read_verifications;
            counters.accepted_mate_placements += part.accepted_mate_placements;
            counters.valid_fr_states += part.valid_fr_states;
            counters.physical_origins += part.physical_origins;
            counters.statistical_state_classes += part.statistical_state_classes;
        }
    }

    // The ACCEPTANCE decision, taken once on the merged totals. A refusal may already have been raised
    // in flight, by a worker that watched the shared budget cross a cap -- so this is not the only place
    // a verdict is reached, and claiming otherwise would be wrong. That early refusal is sound because a
    // shared total can only undercount the exact one. What this block settles is the other half: a run
    // no worker stopped is accepted or refused on exact totals, so the answer is the same at any thread
    // count -- the crossing point moves with the partition, the total does not.
    if (counters.piece_hits > options.max_piece_hits) {
        note_observed_work();
        throw ResourceExhausted(
            "fragment scoring reached its piece-hit cap; the locus is too repetitive for the declared "
            "caps and no partial ranking is returned");
    }
    if (counters.physical_origins > options.max_physical_origins) {
        note_observed_work();
        throw ResourceExhausted(
            "fragment scoring reached its physical-origin cap; no partial ranking is returned");
    }
    // Unlike the two above, this one bounds nothing in flight: state classes are counted after the
    // per-fragment merge, so the work and the memory are already spent by the time it is consulted. It
    // is an acceptance cap on the result, not a breaker, and so it needs no shared budget.
    if (counters.statistical_state_classes > options.max_state_classes) {
        note_observed_work();
        throw ResourceExhausted(
            "fragment scoring reached its statistical-class cap; no partial ranking is returned");
    }

    result.seconds_placement_and_mass =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - placement_mark).count();
    // Phase 3: diploid combination, per (candidate, fragment). This is the term that grows with the
    // panel -- roughly 108,000 candidates on the large loci against 8,515 on c4.
    const auto combination_mark = std::chrono::steady_clock::now();

    // One candidate's score depends on nothing but its own two haplotypes and the shared, now
    // immutable, mass matrix -- so candidates partition cleanly. Each worker takes a fixed contiguous
    // range and writes only its own complete DiplotypeScore entries; within a candidate the fragment
    // order is exactly the serial order, and no value is accumulated across threads. There is
    // therefore no floating-point reduction whose order could change, which is what makes this
    // bit-for-bit identical at any thread count rather than merely close.
    const auto combine_range = [&](std::size_t begin, std::size_t end) {
        for (std::size_t c = begin; c < end; ++c) {
            const std::size_t a = candidates[c].first;
            const std::size_t b = candidates[c].second;
            DiplotypeScore& score = result.scores[c];
            // A homozygote names the same haplotype twice and therefore counts two chromosome copies
            // in both terms, which is what the fragment contract requires.
            score.exposure_term = -expected_signal_fragments(exposure[a] + exposure[b], library);
            score.absolute_term_sum = std::abs(score.exposure_term);
            double fragment_term = 0.0;
            for (std::size_t f = 0; f < fragments.size(); ++f) {
                const double combined =
                    log_add_exp(mass[f * unique + a], mass[f * unique + b]);
                if (combined != kNegInf) ++score.fragments_with_mass;
                const double contribution =
                    log_fragment_contribution(combined, background_by_fragment[f], library);
                fragment_term += contribution;
                score.absolute_term_sum += std::abs(contribution);
            }
            score.fragment_term = fragment_term;
            score.log_likelihood = score.exposure_term + fragment_term;
        }
    };

    const std::size_t combination_threads =
        candidates.empty() ? 1
                           : std::max<std::size_t>(1, std::min(options.threads, candidates.size()));
    result.threads_used_combination = combination_threads;
    result.candidate_fragment_evaluations =
        static_cast<std::uint64_t>(candidates.size()) * static_cast<std::uint64_t>(fragments.size());

    if (combination_threads <= 1) {
        combine_range(0, candidates.size());
    } else {
        const std::size_t span = (candidates.size() + combination_threads - 1) / combination_threads;
        std::vector<std::exception_ptr> failures(combination_threads);
        std::vector<std::thread> workers;
        workers.reserve(combination_threads);
        for (std::size_t t = 0; t < combination_threads; ++t) {
            const std::size_t begin = std::min(t * span, candidates.size());
            const std::size_t end = std::min(begin + span, candidates.size());
            workers.emplace_back([&, t, begin, end]() {
                try {
                    combine_range(begin, end);
                } catch (...) {
                    failures[t] = std::current_exception();
                }
            });
        }
        for (std::thread& worker : workers) worker.join();
        // Every worker is joined before anything is rethrown, and the rethrow walks thread order, so
        // which worker happened to fail first cannot change which refusal the caller sees.
        for (const std::exception_ptr& failure : failures) {
            if (failure) std::rethrow_exception(failure);
        }
    }
    result.seconds_diploid_combination =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - combination_mark).count();
    result.log_mass = std::move(mass);
    result.exposure = std::move(exposure);
    return result;
}

CertifiedRanking certify_diplotype_ranking(
    const InternedHaplotypes& haplotypes,
    const std::vector<HaplotypeCandidate>& candidates,
    const std::vector<ReadPair>& fragments,
    const LibraryModel& library,
    const FragmentScoringOptions& options,
    const CertificationOptions& certification,
    const std::vector<double>& exact_prior) {
    if (!exact_prior.empty() && exact_prior.size() != candidates.size()) {
        throw std::invalid_argument("exact prior must be parallel to the candidate list");
    }
    if (candidates.empty()) {
        throw std::invalid_argument("certification needs at least one candidate to rank");
    }
    const FragmentScoringResult scored =
        score_diplotypes(haplotypes, candidates, fragments, library, options);
    return certify_scored_ranking(scored.scores, fragments.size(), options.mismatch_band, scored.counters,
                                  certification, exact_prior);
}

CertifiedRanking certify_scored_ranking(
    const std::vector<DiplotypeScore>& scores,
    std::size_t fragment_count,
    std::size_t mismatch_band,
    const FragmentScoringCounters& counters,
    const CertificationOptions& certification,
    const std::vector<double>& exact_prior) {
    if (!exact_prior.empty() && exact_prior.size() != scores.size()) {
        throw std::invalid_argument("exact prior must be parallel to the candidate list");
    }
    if (scores.empty()) {
        throw std::invalid_argument("certification needs at least one candidate to rank");
    }
    if (!std::isfinite(certification.minimum_tie_tolerance_nats) ||
        certification.minimum_tie_tolerance_nats < 0.0) {
        throw std::invalid_argument("minimum tie tolerance must be finite and non-negative");
    }
    for (double prior : exact_prior) {
        if (!std::isfinite(prior)) {
            throw std::invalid_argument("an exact prior must be finite to rank its candidate");
        }
    }

    CertifiedRanking result;
    result.counters = counters;
    result.mismatch_band = mismatch_band;
    result.candidates.resize(scores.size());
    const double terms = static_cast<double>(fragment_count + 16);
    double log_weight = kNegInf;
    for (std::size_t c = 0; c < scores.size(); ++c) {
        CertifiedCandidate& entry = result.candidates[c];
        entry.candidate_id = c;
        entry.prior = exact_prior.empty() ? 0.0 : exact_prior[c];
        entry.fragment_term = scores[c].fragment_term;
        entry.exposure_term = scores[c].exposure_term;
        entry.score = scores[c].log_likelihood + entry.prior;
        entry.fragments_with_mass = scores[c].fragments_with_mass;
        entry.numeric_error_bound = terms * std::numeric_limits<double>::epsilon() *
                                    (scores[c].absolute_term_sum + std::abs(entry.prior));
        entry.numeric_error_bound =
            std::max(entry.numeric_error_bound, 0.5 * certification.minimum_tie_tolerance_nats);
        log_weight = log_add_exp(log_weight, entry.score);
    }
    result.candidate_set_log_weight = log_weight;
    for (CertifiedCandidate& entry : result.candidates) {
        entry.candidate_set_probability = std::exp(entry.score - log_weight);
    }

    std::sort(result.candidates.begin(), result.candidates.end(),
              [](const CertifiedCandidate& lhs, const CertifiedCandidate& rhs) {
                  if (lhs.score != rhs.score) return lhs.score > rhs.score;
                  return lhs.candidate_id < rhs.candidate_id;
              });

    const CertifiedCandidate& best = result.candidates.front();
    if (result.candidates.size() == 1) {
        result.status = RescoreStatus::Complete;
        result.winner_set = {best.candidate_id};
        result.certified_gap = std::numeric_limits<double>::infinity();
        result.score_gap_to_runner_up = std::numeric_limits<double>::infinity();
        return result;
    }

    result.score_gap_to_runner_up = best.score - result.candidates[1].score;
    const double best_lower = best.score - best.numeric_error_bound;
    double certified_gap = std::numeric_limits<double>::infinity();
    std::vector<std::size_t> inseparable = {best.candidate_id};
    for (std::size_t i = 1; i < result.candidates.size(); ++i) {
        const CertifiedCandidate& rival = result.candidates[i];
        const double margin = best_lower - (rival.score + rival.numeric_error_bound);
        certified_gap = std::min(certified_gap, margin);
        if (margin <= 0.0) inseparable.push_back(rival.candidate_id);
    }
    if (inseparable.size() == 1) {
        result.status = RescoreStatus::Complete;
        result.winner_set = inseparable;
        result.certified_gap = certified_gap;
    } else {
        // The best candidate cannot be separated from these by more than floating-point error. Taking
        // the highest point estimate here would present an uncertified pick as a call.
        result.status = RescoreStatus::RescoreIncomplete;
        result.winner_set = inseparable;
        result.certified_gap = 0.0;
    }
    return result;
}

} // namespace panvar
