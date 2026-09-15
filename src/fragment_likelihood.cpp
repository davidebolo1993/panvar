#include "panvar/fragment_likelihood.hpp"

#include <algorithm>
#include <map>
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

// Scan one haplotype once per piece-length bucket, turning every piece hit into an implied read
// start. This is the inversion the design calls for: the haplotype is walked, not queried once per
// read against a position hash the size of the locus.
void collect_form_placements(
    const std::string& haplotype,
    const std::vector<ReadForm>& forms,
    const std::vector<PieceIndex>& indexes,
    std::size_t band,
    std::vector<std::vector<FormPlacement>>& out,
    FragmentScoringCounters& counters) {
    std::vector<std::pair<std::uint32_t, std::uint32_t>> proposed;   // (form, start)

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
                ++counters.piece_hits;
                if (it->offset > piece_start) continue;   // implied start would be negative
                proposed.emplace_back(it->form,
                                      static_cast<std::uint32_t>(piece_start - it->offset));
                ++counters.proposed_starts;
            }
        }
    }

    // Several pieces of one read can rediscover one physical placement. Deduplicate before
    // verifying, so a placement found three ways is verified once and counted once.
    std::sort(proposed.begin(), proposed.end());
    proposed.erase(std::unique(proposed.begin(), proposed.end()), proposed.end());
    counters.unique_starts += proposed.size();

    for (const auto& item : proposed) {
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
    const FragmentScoringOptions& options) {
    validate_library_model(library);
    if (options.minimum_piece_bp == 0) {
        throw std::invalid_argument("minimum indexable piece length must be positive");
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

    std::map<std::size_t, double> band_mass_by_length;
    const auto band_mass = [&](std::size_t length) {
        const auto found = band_mass_by_length.find(length);
        if (found != band_mass_by_length.end()) return found->second;
        const double value = log_band_mass(length, band, library.substitution_error);
        band_mass_by_length.emplace(length, value);
        return value;
    };

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

    const std::size_t unique = haplotypes.sequences.size();
    std::vector<double> mass(fragments.size() * unique, kNegInf);
    std::vector<double> exposure(unique, 0.0);

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

    std::map<StateKey, std::uint64_t> states;
    for (std::size_t h = 0; h < unique; ++h) {
        const std::string& haplotype = haplotypes.sequences[h];
        if (forms.empty()) continue;

        std::vector<std::vector<FormPlacement>> placements(forms.size());
        collect_form_placements(haplotype, forms, indexes, band, placements, counters);

        for (std::size_t f = 0; f < fragments.size(); ++f) {
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
                        ++counters.physical_origins;
                        StateKey key;
                        key.insert_bp = static_cast<std::uint32_t>(*insert);
                        key.forward_mismatches = start.mismatches;
                        key.reverse_mismatches = it->mismatches;
                        key.orientation = static_cast<std::uint32_t>(arm);
                        ++states[key];
                    }
                }
            }
            if (states.empty()) continue;
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

    for (std::size_t c = 0; c < candidates.size(); ++c) {
        const std::size_t a = candidates[c].first;
        const std::size_t b = candidates[c].second;
        DiplotypeScore& score = result.scores[c];
        // A homozygote names the same haplotype twice and therefore counts two chromosome copies in
        // both terms, which is what the fragment contract requires.
        score.exposure_term = -library.fragment_start_intensity * (exposure[a] + exposure[b]);
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
    if (!std::isfinite(certification.minimum_tie_tolerance_nats) ||
        certification.minimum_tie_tolerance_nats < 0.0) {
        throw std::invalid_argument("minimum tie tolerance must be finite and non-negative");
    }
    for (double prior : exact_prior) {
        if (!std::isfinite(prior)) {
            throw std::invalid_argument("an exact prior must be finite to rank its candidate");
        }
    }

    const FragmentScoringResult scored =
        score_diplotypes(haplotypes, candidates, fragments, library, options);

    CertifiedRanking result;
    result.counters = scored.counters;
    result.mismatch_band = options.mismatch_band;
    result.candidates.resize(candidates.size());
    const double terms = static_cast<double>(fragments.size() + 16);
    double log_weight = kNegInf;
    for (std::size_t c = 0; c < candidates.size(); ++c) {
        CertifiedCandidate& entry = result.candidates[c];
        entry.candidate_id = c;
        entry.prior = exact_prior.empty() ? 0.0 : exact_prior[c];
        entry.fragment_term = scored.scores[c].fragment_term;
        entry.exposure_term = scored.scores[c].exposure_term;
        entry.score = scored.scores[c].log_likelihood + entry.prior;
        entry.fragments_with_mass = scored.scores[c].fragments_with_mass;
        entry.numeric_error_bound = terms * std::numeric_limits<double>::epsilon() *
                                    (scored.scores[c].absolute_term_sum + std::abs(entry.prior));
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
