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

// A reverse-role mate is compared against the haplotype's forward strand as rc(mate), and is then
// anchored at its RIGHT endpoint by running the recursion over reversed sequences. The reversed form
// of rc(mate) is just its complement, since rc(X) = reverse(complement(X)). Spelling that out here
// keeps the scorer from reversing a string per verification in its hottest loop.
std::string complement_read(const std::string& bases) {
    std::string out(bases.size(), 'N');
    for (std::size_t i = 0; i < bases.size(); ++i) out[i] = complement(bases[i]);
    return out;
}

// One physical origin under construction, while the latent alignments that explain it are summed.
// Shared by the oracle and the accelerated scorer so the two cannot drift in how they collapse.
struct OriginAccumulator {
    double log_latent_mass = kNegInf;
    std::size_t latent_alignments = 0;
    // The best single latent alignment, kept only to describe the origin in diagnostics.
    double best_pair = kNegInf;
    std::size_t reverse_start = 0;
    std::size_t forward_consumed_bp = 0;
    std::size_t reverse_consumed_bp = 0;
    std::size_t forward_mismatches = 0;
    std::size_t reverse_mismatches = 0;

    void add(double pair_mass, std::size_t reverse, std::size_t forward_consumed,
             std::size_t reverse_consumed, std::size_t forward_edits, std::size_t reverse_edits) {
        log_latent_mass = log_add_exp(log_latent_mass, pair_mass);
        ++latent_alignments;
        if (pair_mass <= best_pair) return;
        best_pair = pair_mass;
        reverse_start = reverse;
        forward_consumed_bp = forward_consumed;
        reverse_consumed_bp = reverse_consumed;
        forward_mismatches = forward_edits;
        reverse_mismatches = reverse_edits;
    }
};

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

// -------------------------------------------------------------------------------------------
// Indel-aware emission: a banded affine-gap pair-HMM, conditioned on the observed read length.
// -------------------------------------------------------------------------------------------

namespace {

enum GapState { kMatch = 0, kInsert = 1, kDelete = 2, kStates = 3 };

// The recursion, shared by the emission and its normaliser. `emitting` false replaces every
// emission with one, which is exactly the sum over all 4^L reads and therefore the shape mass.
// Cells are indexed [read position][displacement][state][edits used].
// `available` is how many haplotype bases may be consumed from `start`. It bounds BOTH modes: the
// numerator must not walk off the haplotype, and the normaliser must not count shapes that would.
// Bounding only the emitting pass is what leaked mass near the haplotype ends -- shapes dropped from
// the numerator were still dividing the denominator, so the emission summed to less than one there.
double gap_recursion(
    const std::string* haplotype,
    const std::string* read,
    std::size_t start,
    std::size_t read_length,
    std::size_t available,
    std::size_t band,
    const LibraryModel& library,
    bool emitting,
    std::vector<GappedComponent>* components) {
    const std::size_t delta = library.max_gap_bp;
    const std::size_t width = 2 * delta + 1;
    const std::size_t budget = band + 1;                 // edits 0..band
    const double log_match = std::log(1.0 - library.substitution_error);
    const double log_mismatch = std::log(library.substitution_error / 3.0);
    const double log_insert_emission = std::log(0.25);
    const double log_continue = library.gap_open > 0.0
        ? std::log(1.0 - 2.0 * library.gap_open) : 0.0;
    const double log_open = library.gap_open > 0.0 ? std::log(library.gap_open) : kNegInf;
    const double log_extend = library.gap_extend > 0.0 ? std::log(library.gap_extend) : kNegInf;
    const double log_close = library.gap_extend < 1.0 ? std::log(1.0 - library.gap_extend) : kNegInf;

    const std::size_t cells = width * kStates * budget;
    std::vector<double> previous(cells, kNegInf);
    std::vector<double> current(cells, kNegInf);

    const auto at = [&](std::size_t displacement, std::size_t state, std::size_t edits) {
        return (displacement * kStates + state) * budget + edits;
    };
    // Displacement k is stored offset by delta, so index 0 means k = -delta.
    previous[at(delta, kMatch, 0)] = 0.0;

    for (std::size_t i = 1; i <= read_length; ++i) {
        std::fill(current.begin(), current.end(), kNegInf);
        for (std::size_t d = 0; d < width; ++d) {
            const long k = static_cast<long>(d) - static_cast<long>(delta);
            // Haplotype bases consumed before this read base under M: i - 1 + k.
            const long consumed_before = static_cast<long>(i) - 1 + k;
            for (std::size_t e = 0; e < budget; ++e) {
                // --- M: consume one read base and one haplotype base -------------------------
                if (consumed_before >= 0) {
                    const std::size_t index = start + static_cast<std::size_t>(consumed_before);
                    // The match step, as (weight, edit cost) branches.
                    //
                    // Emitting, the read base is known, so there is exactly one branch. Summing over
                    // ALL 4^L reads there are two, and they must both be walked: a match position
                    // carries (1-e) over the one base that matches, and a mismatch carries 3*(e/3) = e
                    // over the three that do not, at a cost of one edit. Collapsing that to a single
                    // zero-cost branch -- which is what computing the normaliser with "no emission"
                    // does -- makes Z count only the all-match path, and the emission then reduces to
                    // the UNTRUNCATED substitution model instead of the declared truncated one.
                    double weights[2] = {0.0, 0.0};
                    std::size_t costs[2] = {0, 1};
                    std::size_t branches = 0;
                    // The consumption bound applies whether or not we are emitting.
                    bool usable = static_cast<std::size_t>(consumed_before) < available;
                    if (emitting && usable) {
                        {
                            const char reference = (*haplotype)[index];
                            const char observed = (*read)[i - 1];
                            if (!is_acgt(reference) || !is_acgt(observed)) {
                                throw std::invalid_argument(
                                    "gapped emission is defined over ACGT; a non-ACGT base must "
                                    "refuse rather than be scored as a mismatch");
                            }
                            const bool same = reference == observed;
                            weights[0] = same ? log_match : log_mismatch;
                            costs[0] = same ? 0 : 1;
                            branches = 1;
                        }
                    } else {
                        weights[0] = std::log(1.0 - library.substitution_error);
                        costs[0] = 0;
                        weights[1] = std::log(library.substitution_error);
                        costs[1] = 1;
                        branches = 2;
                    }
                    for (std::size_t b = 0; usable && b < branches; ++b) {
                        const double emission = weights[b];
                        const std::size_t cost = costs[b];
                        if (e < cost) continue;
                        const std::size_t from = e - cost;
                        const double arrivals[kStates] = {
                            previous[at(d, kMatch, from)] + log_continue,
                            previous[at(d, kInsert, from)] + log_close,
                            previous[at(d, kDelete, from)] + log_close,
                        };
                        double sum = kNegInf;
                        for (std::size_t s = 0; s < kStates; ++s) {
                            sum = log_add_exp(sum, arrivals[s]);
                        }
                        const std::size_t cell = at(d, kMatch, e);
                        current[cell] = log_add_exp(current[cell], sum + emission);
                    }
                }
                // --- I: consume one read base only, displacement falls by one ----------------
                if (d + 1 < width && e >= 1) {
                    const double emission = emitting ? log_insert_emission : 0.0;
                    const std::size_t from = e - 1;
                    const double opened = previous[at(d + 1, kMatch, from)] + log_open;
                    const double extended = previous[at(d + 1, kInsert, from)] + log_extend;
                    const std::size_t cell = at(d, kInsert, e);
                    current[cell] = log_add_exp(current[cell], log_add_exp(opened, extended) + emission);
                }
            }
        }
        // --- D: consume haplotype only. No read base advances, so it is resolved within the
        // current read position, sweeping displacement upward because D raises k by one.
        for (std::size_t d = 1; d < width; ++d) {
            const long k = static_cast<long>(d) - static_cast<long>(delta);
            const long consumed_before = static_cast<long>(i) + k - 1;
            if (consumed_before < 0) continue;
            if (static_cast<std::size_t>(consumed_before) >= available) continue;
            for (std::size_t e = 1; e < budget; ++e) {
                const std::size_t from = e - 1;
                const double opened = current[at(d - 1, kMatch, from)] + log_open;
                const double extended = current[at(d - 1, kDelete, from)] + log_extend;
                const std::size_t cell = at(d, kDelete, e);
                current[cell] = log_add_exp(current[cell], log_add_exp(opened, extended));
            }
        }
        previous.swap(current);
    }

    double total = kNegInf;
    // The displacement band bounds consumption from BELOW as well as above. Without the floor an
    // alignment may terminate having consumed almost no haplotype -- a read absorbed entirely by
    // insertions -- which is not a placement in any useful sense: near the end of a haplotype every
    // start would then admit the read, and those degenerate origins put a whole extra unit of mass
    // into the exposure sum. With it, a read that cannot fit yields no placement, exactly as the
    // substitution-only model behaved when count_mismatches ran off the end.
    const std::size_t floor_bp = read_length - std::min(read_length, delta);
    // Each terminal displacement is one CONSUMED-SPAN CLASS, and they are kept apart here rather than
    // summed into one number. They have to be: every class ends the mate at a different haplotype
    // coordinate, so each implies a different fragment geometry and a different insert probability.
    // Collapsing them and handing the total to the best class's span weights the other classes' mass
    // with an insert probability that is not theirs -- which is what made the signal density fail to
    // integrate to the exposure. The caller applies the geometry to each class separately.
    for (std::size_t d = 0; d < width; ++d) {
        const long k = static_cast<long>(d) - static_cast<long>(delta);
        const long consumed = static_cast<long>(read_length) + k;
        if (consumed < 0 || static_cast<std::size_t>(consumed) > available) continue;
        if (static_cast<std::size_t>(consumed) < floor_bp) continue;
        double span = kNegInf;
        // Fewest edits reaching this span. A cell is INDEXED by the edits spent to reach it, so the
        // lowest finite e in the class is that minimum exactly -- there is nothing to track alongside
        // the recursion, and the arrays that used to do it could only restate their own index.
        std::size_t span_edits = budget;
        for (std::size_t s = 0; s < kStates; ++s) {
            for (std::size_t e = 0; e < budget; ++e) {
                const double cell = previous[at(d, s, e)];
                if (cell == kNegInf) continue;
                span = log_add_exp(span, cell);
                if (e < span_edits) span_edits = e;
            }
        }
        if (span == kNegInf) continue;
        total = log_add_exp(total, span);
        if (components != nullptr) {
            components->push_back({static_cast<std::size_t>(consumed), span, span_edits});
        }
    }
    return total;
}

} // namespace

double log_gap_shape_mass(
    std::size_t read_length,
    std::size_t available_bp,
    std::size_t mismatch_band,
    const LibraryModel& library) {
    if (read_length == 0) throw std::invalid_argument("gapped emission requires a non-empty read");
    // ONE normaliser, summed over every terminal span. Deliberately not per span: the components are
    // parts of a single distribution over reads of this length, not several distributions, and
    // dividing each by its own shape mass would make each of them sum to one on its own and their
    // total several times that.
    //
    // MEMOISED, because this is a pure function of its arguments and gapped_emission calls it on
    // every verification -- which made each verification run the whole recursion TWICE, once for the
    // numerator and once for a denominator that depends on neither sequence. Measured at 175us of a
    // 255us call: 69% of the cost of every placement in the program.
    //
    // The key space is tiny once `available_bp` is clamped. The haplotype end cuts admissible shapes
    // only while it is within reach; once available_bp reaches read_length + max_gap_bp, no terminal
    // span and no intermediate path is excluded, and every larger value yields a BIT-IDENTICAL
    // result -- measured across available_bp from 103 to 200,000 at read_length 100, max_gap_bp 3.
    // So the clamp leaves exactly max_gap_bp + 1 distinct entries per read length, and it is an
    // identity rather than a tolerance.
    const std::size_t key_available = std::min(available_bp, read_length + library.max_gap_bp);

    struct ShapeKey {
        std::size_t read_length;
        std::size_t available;
        std::size_t band;
        std::size_t max_gap;
        double substitution_error;
        double gap_open;
        double gap_extend;

        bool operator<(const ShapeKey& other) const noexcept {
            return std::tie(read_length, available, band, max_gap, substitution_error, gap_open,
                            gap_extend) <
                   std::tie(other.read_length, other.available, other.band, other.max_gap,
                            other.substitution_error, other.gap_open, other.gap_extend);
        }
    };
    // Every parameter the value depends on is in the key, so a stale hit is impossible: change the
    // band or any gap parameter and it is a different entry, not a reused one. Thread-local because
    // the placement phase runs on several workers and a shared cache would need a lock in the
    // hottest path in the program; each worker paying for four entries is not worth contending for.
    static thread_local std::map<ShapeKey, double> cache;
    const ShapeKey key{read_length,   key_available,       mismatch_band,     library.max_gap_bp,
                       library.substitution_error, library.gap_open, library.gap_extend};
    const auto found = cache.find(key);
    if (found != cache.end()) return found->second;
    const double value = gap_recursion(nullptr, nullptr, 0, read_length, key_available,
                                       mismatch_band, library, false, nullptr);
    cache.emplace(key, value);
    return value;
}

std::optional<GappedFit> gapped_emission(
    const std::string& haplotype,
    const std::string& read,
    std::size_t start,
    const LibraryModel& library,
    std::size_t mismatch_band) {
    if (read.empty()) throw std::invalid_argument("gapped emission requires a non-empty read");
    if (start > haplotype.size()) return std::nullopt;
    // ELIGIBILITY is by read length, exactly as the substitution-only model required: the haplotype
    // must be able to hold the read from here. The alignment inside that window is still free to
    // consume read_length +/- max_gap_bp, which is what admits indels.
    //
    // Admitting starts with only read_length - max_gap_bp bases left looks more permissive but is
    // wrong twice over. Such a start can be reached only by absorbing the overhang as insertions, and
    // because the normaliser is conditioned on what is available there, those gapped shapes are the
    // ONLY admissible ones and normalise to near certainty -- a gap costing 1e-4 in the interior
    // arrives at ~0.2 at the boundary. Summed over reads that manufactured exactly one extra unit of
    // density at the haplotype end, which is why the signal integrated to exposure + 1 regardless of
    // haplotype length or gap width.
    if (haplotype.size() - start < read.size()) return std::nullopt;
    // Both passes see the same bound on what the haplotype can supply from here.
    const std::size_t available = haplotype.size() - start;
    std::vector<GappedComponent> raw;
    const double raw_total = gap_recursion(&haplotype, &read, start, read.size(), available,
                                           mismatch_band, library, true, &raw);
    if (raw_total == kNegInf || raw.empty()) return std::nullopt;
    const double shapes = log_gap_shape_mass(read.size(), available, mismatch_band, library);

    GappedFit fit;
    // Every component is divided by the SAME normaliser, so their sum is the whole read's emission
    // and the decomposition changes where the mass is attached without changing how much there is.
    fit.log_emission = raw_total - shapes;
    fit.components.reserve(raw.size());
    double best = kNegInf;
    for (GappedComponent component : raw) {
        component.log_emission -= shapes;
        if (component.log_emission > best) {
            best = component.log_emission;
            fit.consumed_bp = component.consumed_bp;
            fit.edits = component.edits;
        }
        fit.components.push_back(component);
    }
    return fit;
}

std::optional<GappedFit> gapped_emission_ending_at(
    const std::string& haplotype,
    const std::string& read,
    std::size_t end,
    const LibraryModel& library,
    std::size_t mismatch_band) {
    if (read.empty()) throw std::invalid_argument("gapped emission requires a non-empty read");
    if (end > haplotype.size()) return std::nullopt;
    // Eligibility by read length, measured leftward: the haplotype must be able to hold the read
    // ending here, exactly as the left-anchored form requires it to hold the read starting there.
    if (end < read.size()) return std::nullopt;
    // Reversing BOTH sequences turns "read leftward from `end`" into "read rightward from
    // haplotype.size() - end", which is the recursion already written. Reversing both preserves the
    // match/mismatch pattern exactly; what changes is the direction the gap model runs in, and that
    // is the whole point -- the mate really is sequenced from this end.
    const std::string reversed_haplotype(haplotype.rbegin(), haplotype.rend());
    const std::string reversed_read(read.rbegin(), read.rend());
    return gapped_emission(reversed_haplotype, reversed_read, haplotype.size() - end, library,
                           mismatch_band);
}

std::optional<std::size_t> fr_insert_from_endpoints(
    std::size_t forward_start,
    std::size_t forward_consumed_bp,
    std::size_t reverse_end,
    std::size_t reverse_consumed_bp,
    const InsertSupport& support) {
    if (reverse_end <= forward_start) return std::nullopt;
    const std::size_t insert = reverse_end - forward_start;
    // Overlapping mates are a physical state, so the floor is the longer mate's CONSUMED span,
    // never the sum. A gapped mate occupies a different stretch than its read length, so the floor
    // has to be stated on what was consumed.
    if (insert < std::max(forward_consumed_bp, reverse_consumed_bp)) return std::nullopt;
    if (support.empty() || insert < support.min_bp || insert > support.max_bp) return std::nullopt;
    return insert;
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
    const LibraryModel& library,
    std::size_t mismatch_band) {
    validate_read_pair(fragment);
    validate_library_model(library);
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

        // Each anchor's components once per form. The pairing below reads them O(L^2) times, and
        // recomputing the recursion inside that loop made the oracle quadratic in DP evaluations
        // rather than merely in pairings.
        //
        // The two mates are anchored at OPPOSITE ends, because that is how the template was read:
        // the forward mate rightward from the left endpoint, the reverse mate leftward from the
        // right endpoint. Anchoring both at a start is what broke the exposure identity.
        std::vector<std::optional<GappedFit>> forward_fits(haplotype.size() + 1);
        std::vector<std::optional<GappedFit>> reverse_fits(haplotype.size() + 1);
        for (std::size_t at = 0; at <= haplotype.size(); ++at) {
            if (at < haplotype.size()) {
                forward_fits[at] = gapped_emission(haplotype, forward_bases, at, library,
                                                   mismatch_band);
            }
            reverse_fits[at] =
                gapped_emission_ending_at(haplotype, reverse_bases, at, library, mismatch_band);
        }

        // Keyed by (left endpoint, insert). The right endpoint is their sum, so this names the
        // template's two ends -- which is what a physical origin is. Every latent alignment pair
        // arriving at the same pair of ends is one more explanation of ONE template, so it is summed
        // into that origin rather than producing a second one.
        std::map<std::pair<std::size_t, std::size_t>, OriginAccumulator> origins;
        for (std::size_t forward_start = 0; forward_start < haplotype.size(); ++forward_start) {
            const std::optional<GappedFit>& forward_fit = forward_fits[forward_start];
            if (!forward_fit.has_value()) continue;
            for (std::size_t reverse_end = 0; reverse_end <= haplotype.size(); ++reverse_end) {
                const std::optional<GappedFit>& reverse_fit = reverse_fits[reverse_end];
                if (!reverse_fit.has_value()) continue;
                // Geometry follows the bases each mate CONSUMES, not the bases it carries, and each
                // consumed-span class is given its own geometry: a gapped mate spans a different
                // stretch of haplotype than its read length, and a different one per class.
                for (const GappedComponent& forward_span : forward_fit->components) {
                    if (forward_span.log_emission == kNegInf) continue;
                    for (const GappedComponent& reverse_span : reverse_fit->components) {
                        if (reverse_span.log_emission == kNegInf) continue;
                        const std::optional<std::size_t> insert = fr_insert_from_endpoints(
                            forward_start, forward_span.consumed_bp,
                            reverse_end, reverse_span.consumed_bp, support);
                        if (!insert.has_value()) continue;
                        origins[{forward_start, *insert}].add(
                            forward_span.log_emission + reverse_span.log_emission,
                            reverse_end - reverse_span.consumed_bp,
                            forward_span.consumed_bp, reverse_span.consumed_bp,
                            forward_span.edits, reverse_span.edits);
                    }
                }
            }
        }

        for (const auto& item : origins) {
            FragmentPlacement placement;
            placement.orientation = orientation;
            placement.forward_start = item.first.first;
            placement.insert_bp = item.first.second;
            placement.log_latent_mass = item.second.log_latent_mass;
            placement.latent_alignments = item.second.latent_alignments;
            placement.reverse_start = item.second.reverse_start;
            placement.forward_consumed_bp = item.second.forward_consumed_bp;
            placement.reverse_consumed_bp = item.second.reverse_consumed_bp;
            placement.forward_mismatches = item.second.forward_mismatches;
            placement.reverse_mismatches = item.second.reverse_mismatches;
            placements.push_back(placement);
        }
    }
    return placements;
}

double log_placement_mass(
    const std::vector<FragmentPlacement>& placements,
    const ReadPair& fragment,
    const DiscreteInsertDistribution& insert,
    const LibraryModel& library) {
    validate_read_pair(fragment);
    validate_library_model(library);
    double mass = kNegInf;
    for (const FragmentPlacement& placement : placements) {
        // One template, so one orientation prior and one insert probability. The latent alignments
        // were summed into log_latent_mass at the geometry each of them implied, which is the whole
        // point of carrying them separately that far.
        const double term =
            kLogMateOrientationPrior +
            insert.log_probability(placement.insert_bp) +
            placement.log_latent_mass;
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

// One verified consumed-span class of one read form at one start on one haplotype. There is a row
// per class, not per start: the classes end the mate at different coordinates and so pair with
// different reverse mates, and merging them before the geometry is applied is exactly the defect
// this file was corrected for. The emission is carried rather than a mismatch count, because under
// the indel-aware model the count no longer determines it.
struct FormPlacement {
    std::uint32_t start = 0;
    std::uint32_t consumed = 0;        // haplotype bases this alignment class spans
    std::uint32_t mismatches = 0;      // diagnostic only
    double emission = 0.0;             // normalised log mass of this class
};

// The exact statistical class of a physical origin. Origins sharing every field may be stored once
// with a multiplicity, but ONLY because every term of their mass is then identical. Keying on
// mismatch counts was sufficient when the count determined the emission; it no longer does, so the
// summed latent mass is itself part of the key. Two origins with equal counts but different gapped
// alignments are different classes and must not be collapsed.
//
// The multiplicity here counts PLACES -- distinct left endpoints carrying the same insert and the
// same latent mass, which is what repeat copies look like. It never counts alternative alignments of
// one place: those were summed into log_latent_mass before this key was formed.
struct StateKey {
    std::uint32_t insert_bp = 0;
    std::uint32_t orientation = 0;
    double log_latent_mass = 0.0;

    bool operator<(const StateKey& other) const noexcept {
        return std::tie(insert_bp, orientation, log_latent_mass) <
               std::tie(other.insert_bp, other.orientation, other.log_latent_mass);
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
    // Parallel to `forms`: whether this form's consumed-span classes must be carried separately.
    // False only when the fragment's insert support makes the geometry's floor unreachable, where
    // emitting the whole emission once per anchor is exactly equal and far cheaper.
    const std::vector<char>& componentise,
    const std::vector<PieceIndex>& indexes,
    PlacementBudget& budget,
    std::size_t band,
    const LibraryModel& library,
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
                // The seed implies a start only under a gapless alignment. An indel BEFORE the seed
                // shifts the true start by up to max_gap_bp in either direction, so the whole
                // neighbourhood is proposed. Pigeonhole is unaffected -- d edits still break at most
                // d pieces, so no read the substitution-only scan found is lost -- and the extra
                // starts are deduplicated below before any verification happens.
                const std::size_t implied = piece_start - it->offset;
                const std::size_t lo = implied > library.max_gap_bp ? implied - library.max_gap_bp : 0;
                const std::size_t hi = implied + library.max_gap_bp;
                for (std::size_t shifted = lo; shifted <= hi; ++shifted) {
                    proposed.emplace_back(it->form, static_cast<std::uint32_t>(shifted));
                    ++counters.proposed_starts;
                }
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
        const std::optional<GappedFit> fit =
            gapped_emission(haplotype, forms[item.first].bases, item.second, library, band);
        if (!fit.has_value() || fit->log_emission == kNegInf) continue;
        if (componentise[item.first] != 0) {
            // One row per consumed-span class. The scorer must see the same decomposition the oracle
            // sees, or the two agree only where no alignment is gapped.
            for (const GappedComponent& component : fit->components) {
                if (component.log_emission == kNegInf) continue;
                FormPlacement placement;
                placement.start = item.second;
                placement.consumed = static_cast<std::uint32_t>(component.consumed_bp);
                placement.mismatches = static_cast<std::uint32_t>(component.edits);
                placement.emission = component.log_emission;
                out[item.first].push_back(placement);
                ++counters.accepted_mate_placements;
            }
            continue;
        }
        // The exact fast path, not an approximation. The insert is the distance between the two
        // ANCHORS, so a consumed span changes the likelihood only through the floor
        // insert >= max(consumed spans). Where the declared support cannot reach that floor, every
        // span pair at an anchor pair is admitted, the sum factorises into the two mates' whole
        // emissions, and carrying the classes separately would compute the same number through
        // (2*max_gap_bp+1)^2 times as many terms. Measured: identical to the last bit above the
        // floor, and genuinely different below it -- which is why the predicate, not a tolerance,
        // decides. The consumed span recorded here is the highest-mass class, used only for the
        // floor check that cannot bind and for diagnostics.
        FormPlacement placement;
        placement.start = item.second;
        placement.consumed = static_cast<std::uint32_t>(fit->consumed_bp);
        placement.mismatches = static_cast<std::uint32_t>(fit->edits);
        placement.emission = fit->log_emission;
        out[item.first].push_back(placement);
        ++counters.accepted_mate_placements;
    }
    // Sorted by start first, so the pairing below can walk one left endpoint's classes as a group,
    // and by consumed span within it so the row order is fixed rather than insertion-dependent.
    for (std::vector<FormPlacement>& list : out) {
        std::sort(list.begin(), list.end(), [](const FormPlacement& lhs, const FormPlacement& rhs) {
            return std::tie(lhs.start, lhs.consumed) < std::tie(rhs.start, rhs.consumed);
        });
    }
}

} // namespace

PieceHitCensus census_piece_hits(
    const InternedHaplotypes& haplotypes,
    const std::vector<ReadPair>& fragments,
    const FragmentScoringOptions& options) {
    PieceHitCensus census;
    const std::size_t band = options.mismatch_band;

    // The same four oriented forms per fragment the scorer builds: two forward-role mates scanned
    // against the haplotype, and two reverse-role mates stored reversed -- complement(mate) -- and
    // scanned against the reversed haplotype. Both are indexed together here because a piece hit
    // costs the same whichever direction its form is scanned in, and the census exists to count
    // hits; the haplotype is scanned in both orientations below so the total matches the scorer's.
    std::vector<ReadForm> forms;
    forms.reserve(fragments.size() * 4);
    for (std::size_t f = 0; f < fragments.size(); ++f) {
        validate_read_pair(fragments[f]);
        forms.push_back({f, fragments[f].first.bases, true});
        forms.push_back({f, complement_read(fragments[f].second.bases), false});
        forms.push_back({f, fragments[f].second.bases, false});
        forms.push_back({f, complement_read(fragments[f].first.bases), true});
    }
    census.total_forms = forms.size();
    if (forms.empty()) return census;

    // A collapse keeps exactly the first form of each distinct sequence; every other form would be
    // answered from it. Marking them is what makes the prediction exact rather than proportional.
    std::set<std::string> seen;
    std::vector<char> first_of_sequence(forms.size(), 0);
    for (std::size_t i = 0; i < forms.size(); ++i) {
        if (seen.insert(forms[i].bases).second) first_of_sequence[i] = 1;
    }
    census.unique_forms = seen.size();

    const std::vector<PieceIndex> indexes = build_piece_indexes(forms, band, options.minimum_piece_bp);
    for (const PieceIndex& index : indexes) {
        census.total_entries += index.entries.size();
        for (const PieceEntry& entry : index.entries) {
            census.unique_entries += first_of_sequence[entry.form];
        }
    }

    for (const std::string& forward_haplotype : haplotypes.sequences) {
      const std::string reversed_haplotype(forward_haplotype.rbegin(), forward_haplotype.rend());
      std::uint64_t hits_here = 0;
      for (int direction = 0; direction < 2; ++direction) {
        const std::string& haplotype = direction == 0 ? forward_haplotype : reversed_haplotype;
        for (const PieceIndex& index : indexes) {
            const std::size_t piece = index.piece_bp;
            if (piece == 0 || haplotype.size() < piece) continue;
            const std::uint64_t mask = piece >= 32 ? ~0ULL : ((1ULL << (2 * piece)) - 1);
            std::uint64_t code = 0;
            std::size_t valid = 0;
            for (std::size_t at = 0; at < haplotype.size(); ++at) {
                const std::uint64_t bits = two_bit(haplotype[at]);
                if (bits == kNoCode) {
                    valid = 0;
                    code = 0;
                    continue;
                }
                code = ((code << 2) | bits) & mask;
                if (++valid < piece) continue;
                const auto range = std::equal_range(
                    index.entries.begin(), index.entries.end(), code,
                    [](const auto& lhs, const auto& rhs) {
                        if constexpr (std::is_same_v<std::decay_t<decltype(lhs)>, PieceEntry>) {
                            return lhs.code < rhs;
                        } else {
                            return lhs < rhs.code;
                        }
                    });
                std::uint64_t all = 0;
                std::uint64_t unique = 0;
                for (auto it = range.first; it != range.second; ++it) {
                    ++all;
                    unique += first_of_sequence[it->form];
                }
                census.predicted_piece_hits_all += all;
                census.predicted_piece_hits_unique += unique;
                hits_here += all;
                census.max_piece_hits_at_one_seed = std::max(census.max_piece_hits_at_one_seed, all);
            }
        }
      }
      census.max_piece_hits_one_haplotype =
          std::max(census.max_piece_hits_one_haplotype, hits_here);
    }
    return census;
}

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
        options.max_mass_cells == 0 || options.max_latent_alignment_pairs == 0) {
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

    // The per-length band-mass cache is gone with the substitution-only emission that needed it. The
    // indel-aware emission normalises itself -- gapped_emission divides by the shape mass of the
    // read's own length before returning -- so there is nothing left for the scorer to divide by, and
    // a cache here would only be a second place for the normaliser to disagree with the oracle.

    // Phase 1 of three: per-fragment preparation. Timed so a caller can attribute rescoring rather
    // than report it whole -- it is 57-93% of a run, and which phase dominates decides how it should
    // be parallelised.
    const auto preparation_mark = std::chrono::steady_clock::now();

    // Four oriented forms per fragment: each mate as given, and each pre-complemented. Whichever
    // mate reads off the forward strand is compared as given, so every comparison is forward.
    // Two form sets, because the two mates of an FR pair are anchored at OPPOSITE ends of the
    // template and their emissions are therefore computed by scanning in opposite directions.
    //
    //   forward-role, slot 2f+arm: the mate that reads rightward off the forward strand. Scanned
    //                              against the haplotype, anchored at its start.
    //   reverse-role, slot 2f+arm: the mate that reads leftward from the right endpoint, compared
    //                              as rc(mate) and stored already reversed -- which is complement
    //                              (mate). Scanned against the REVERSED haplotype, so the same
    //                              left-anchored machinery anchors it at its right endpoint once
    //                              its coordinates are mapped back.
    //
    // Anchoring both mates at a start is what broke the exposure identity, and the two paths must
    // make the same choice or the scorer and the oracle agree only where no alignment is gapped.
    std::vector<ReadForm> forward_forms;
    std::vector<ReadForm> reverse_forms;
    forward_forms.reserve(fragments.size() * 2);
    reverse_forms.reserve(fragments.size() * 2);
    std::vector<InsertSupport> support_by_fragment;
    std::vector<double> background_by_fragment;
    std::vector<char> floor_binds;
    std::map<std::size_t, std::size_t> fragments_by_floor;
    support_by_fragment.reserve(fragments.size());
    background_by_fragment.reserve(fragments.size());
    floor_binds.reserve(fragments.size());
    for (std::size_t f = 0; f < fragments.size(); ++f) {
        validate_read_pair(fragments[f]);
        // arm 0: mate one reads forward, mate two is the reverse mate.
        forward_forms.push_back({f, fragments[f].first.bases, true});
        reverse_forms.push_back({f, complement_read(fragments[f].second.bases), false});
        // arm 1: the other way round.
        forward_forms.push_back({f, fragments[f].second.bases, false});
        reverse_forms.push_back({f, complement_read(fragments[f].first.bases), true});
        support_by_fragment.push_back(fragment_insert_support(fragments[f], library));
        // Can the geometry's floor reach this fragment at all? The widest span either mate can
        // consume is its read length plus the displacement band; if the support's shortest insert
        // already exceeds that, no admissible insert can fall below a consumed span, and the
        // consumed-span decomposition cannot change this fragment's likelihood.
        {
            const std::size_t widest_span =
                std::max(fragments[f].first.bases.size(), fragments[f].second.bases.size()) +
                library.max_gap_bp;
            const bool binds = support_by_fragment.back().empty() ||
                               support_by_fragment.back().min_bp < widest_span;
            floor_binds.push_back(binds ? 1 : 0);
            if (binds) ++counters.fragments_with_floor_binding;
        }
        background_by_fragment.push_back(log_background_probability(fragments[f]));
        ++fragments_by_floor[std::max(fragments[f].first.bases.size(),
                                      fragments[f].second.bases.size())];
    }
    const auto form_slot = [](std::size_t fragment, std::size_t arm) { return fragment * 2 + arm; };

    // Parallel to each form set: a form takes the exact fast path exactly when its fragment's insert
    // floor is unreachable. Decided per fragment, once, from declared quantities only -- never from
    // the data, and never from a tolerance.
    std::vector<char> forward_componentise(forward_forms.size(), 0);
    std::vector<char> reverse_componentise(reverse_forms.size(), 0);
    for (std::size_t i = 0; i < forward_forms.size(); ++i) {
        forward_componentise[i] = floor_binds[forward_forms[i].fragment];
    }
    for (std::size_t i = 0; i < reverse_forms.size(); ++i) {
        reverse_componentise[i] = floor_binds[reverse_forms[i].fragment];
    }

    const std::vector<PieceIndex> forward_indexes =
        forward_forms.empty() ? std::vector<PieceIndex>{}
                              : build_piece_indexes(forward_forms, band, options.minimum_piece_bp);
    const std::vector<PieceIndex> reverse_indexes =
        reverse_forms.empty() ? std::vector<PieceIndex>{}
                              : build_piece_indexes(reverse_forms, band, options.minimum_piece_bp);

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
    // Right endpoints reachable from one left endpoint, reused across groups to keep the allocation
    // out of the inner loop.
    std::map<std::uint32_t, OriginAccumulator> origins;
    for (std::size_t h = first_h; h < last_h; ++h) {
        // Haplotype and fragment boundaries are the natural places to notice someone else's refusal:
        // both are cheap relative to the work they precede, and a worker that has just finished one
        // haplotype would otherwise start another that can never be used.
        if (budget.stop()) throw PlacementCancelled();
        const std::string& haplotype = haplotypes.sequences[h];
        if (forward_forms.empty()) continue;

        std::vector<std::vector<FormPlacement>> forward_placements(forward_forms.size());
        collect_form_placements(haplotype, forward_forms, forward_componentise, forward_indexes,
                                budget, band, library, forward_placements, counters, options);

        // The reverse-role scan runs against the reversed haplotype, which anchors those mates at
        // their right endpoints. One reversal per haplotype, not per verification.
        const std::string reversed_haplotype(haplotype.rbegin(), haplotype.rend());
        std::vector<std::vector<FormPlacement>> reverse_placements(reverse_forms.size());
        collect_form_placements(reversed_haplotype, reverse_forms, reverse_componentise,
                                reverse_indexes, budget, band, library, reverse_placements, counters,
                                options);
        // Map reversed coordinates back: a class occupying reversed [s, s+consumed) occupies forward
        // [L-s-consumed, L-s), so its forward-coordinate RIGHT ENDPOINT is L - s. `start` carries
        // that endpoint from here on, and the list is re-sorted because reversing the axis reverses
        // the order the two-pointer window below depends on.
        for (std::vector<FormPlacement>& list : reverse_placements) {
            for (FormPlacement& placement : list) {
                placement.start =
                    static_cast<std::uint32_t>(haplotype.size() - placement.start);
            }
            std::sort(list.begin(), list.end(),
                      [](const FormPlacement& lhs, const FormPlacement& rhs) {
                          return std::tie(lhs.start, lhs.consumed) <
                                 std::tie(rhs.start, rhs.consumed);
                      });
        }

        for (std::size_t f = 0; f < fragments.size(); ++f) {
            if (budget.stop()) throw PlacementCancelled();
            const InsertSupport& support = support_by_fragment[f];
            if (support.empty()) continue;
            states.clear();

            for (std::size_t arm = 0; arm < 2; ++arm) {
                // FirstForward: mate one as given is the forward mate, rc(mate two) the reverse one.
                // `forward` rows carry a START, `reverse` rows carry a RIGHT ENDPOINT.
                const std::vector<FormPlacement>& forward = forward_placements[form_slot(f, arm)];
                const std::vector<FormPlacement>& reverse = reverse_placements[form_slot(f, arm)];
                if (forward.empty() || reverse.empty()) continue;

                // Forward rows are grouped by start, because one start is one LEFT ENDPOINT and all
                // of its consumed-span classes are latent explanations of templates anchored there.
                for (std::size_t first_row = 0; first_row < forward.size();) {
                    const std::uint32_t forward_start = forward[first_row].start;
                    std::size_t last_row = first_row;
                    while (last_row < forward.size() && forward[last_row].start == forward_start) {
                        ++last_row;
                    }

                    // insert = right endpoint - left endpoint must lie in the support, so the
                    // admissible right endpoints form one contiguous window. Two-pointer bounds, no
                    // map. The window depends only on the left endpoint, so it is derived once per
                    // group.
                    //
                    // No max_gap_bp slack is needed any more, and adding it back would be a mistake.
                    // The old widening existed only because a reverse START had to be turned into an
                    // endpoint through a consumed span that was not yet known; the reverse rows now
                    // carry the endpoint itself, so these bounds are exact. fr_insert_from_endpoints
                    // still decides admission on the spans actually consumed.
                    const std::size_t lo_end =
                        static_cast<std::size_t>(forward_start) + support.min_bp;
                    const std::size_t hi_end =
                        static_cast<std::size_t>(forward_start) + support.max_bp;

                    const auto begin = std::lower_bound(
                        reverse.begin(), reverse.end(), lo_end,
                        [](const FormPlacement& lhs, std::size_t value) {
                            return lhs.start < value;
                        });
                    const auto end = std::upper_bound(
                        begin, reverse.end(), hi_end,
                        [](std::size_t value, const FormPlacement& rhs) {
                            return value < rhs.start;
                        });

                    // One entry per RIGHT ENDPOINT reachable from this left endpoint. Every latent
                    // alignment pair that reaches the same right endpoint -- whatever spans the two
                    // mates consumed on the way -- is another explanation of ONE template, so the
                    // masses are summed and the origin is counted once. Letting each pair be its own
                    // origin would manufacture physical copies out of alignment ambiguity.
                    origins.clear();
                    for (std::size_t row = first_row; row < last_row; ++row) {
                        for (auto it = begin; it != end; ++it) {
                            // Counted whether or not the geometry accepts: this is the pairing WORK
                            // that admitting indels multiplies, and it has a declared cap of its own.
                            ++counters.latent_alignment_pairs;
                            const std::optional<std::size_t> insert = fr_insert_from_endpoints(
                                forward_start, forward[row].consumed, it->start, it->consumed,
                                support);
                            if (!insert.has_value()) continue;
                            ++counters.valid_fr_states;
                            origins[static_cast<std::uint32_t>(*insert)].add(
                                forward[row].emission + it->emission, it->start - it->consumed,
                                forward[row].consumed, it->consumed, forward[row].mismatches,
                                it->mismatches);
                        }
                    }

                    for (const auto& item : origins) {
                        // Thirty times cooler than piece_hits (9.4M against 320M on C4), so this one
                        // carries an exact shared counter: the work bound is tight, with no batching
                        // overshoot at all. Crossing it cancels every worker rather than letting each
                        // consume a full cap. Refusing here is exact and immediate; acceptance is still
                        // decided later, on the merged totals.
                        ++counters.physical_origins;
                        const std::uint64_t seen =
                            budget.physical_origins.fetch_add(1, std::memory_order_relaxed) + 1;
                        if (seen > budget.max_physical_origins) {
                            budget.cancel();
                            throw ResourceExhausted(
                                "fragment scoring reached its physical-origin cap; no partial "
                                "ranking is returned");
                        }
                        if (budget.stop()) throw PlacementCancelled();
                        StateKey key;
                        key.insert_bp = item.first;
                        key.orientation = static_cast<std::uint32_t>(arm);
                        key.log_latent_mass = item.second.log_latent_mass;
                        ++states[key];
                        if (options.record_origin_summary || record_origins[f] != 0) {
                            FragmentPlacement placement;
                            placement.forward_start = forward_start;
                            placement.insert_bp = item.first;
                            placement.log_latent_mass = item.second.log_latent_mass;
                            placement.latent_alignments = item.second.latent_alignments;
                            placement.reverse_start = item.second.reverse_start;
                            placement.forward_mismatches = item.second.forward_mismatches;
                            placement.reverse_mismatches = item.second.reverse_mismatches;
                            placement.forward_consumed_bp = item.second.forward_consumed_bp;
                            placement.reverse_consumed_bp = item.second.reverse_consumed_bp;
                            placement.orientation =
                                arm == 0 ? MateOrientation::FirstForward : MateOrientation::SecondForward;
                            if (options.record_origin_summary) {
                                OriginSummary& summary = result.origin_summary[f * unique + h];
                                ++summary.origins;
                                const double term =
                                    kLogMateOrientationPrior +
                                    library_insert.log_probability(item.first) +
                                    item.second.log_latent_mass;
                                if (term > summary.best_log_term) {
                                    summary.best_log_term = term;
                                    summary.best = placement;
                                }
                            }
                            if (record_origins[f] != 0) result.recorded_origins[{f, h}].push_back(placement);
                        }
                    }
                    first_row = last_row;
                }
            }
            if (states.empty()) continue;
            // A pure running total with no memory role, so it is accumulated here and decided once,
            // after the ordered merge, where the sum is the same at any thread count.
            counters.statistical_state_classes += states.size();

            double total = kNegInf;
            for (const auto& item : states) {
                // Exact log-multiplicity: every origin in this class has an identical mass term, so
                // multiplying by the count is not an approximation of summing them, it IS the sum.
                // The emissions are part of the key, so two origins collapse only when their masses
                // agree exactly -- which is what makes that still true under gapped alignments.
                const double term =
                    std::log(static_cast<double>(item.second)) + kLogMateOrientationPrior +
                    library_insert.log_probability(item.first.insert_bp) +
                    item.first.log_latent_mass;
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
            counters.latent_alignment_pairs += part.latent_alignment_pairs;
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
    // Like the class cap, an acceptance cap rather than a breaker: the pairs are counted as the
    // geometry is tested, so a running crossing point would move with the partition. It is declared
    // separately from the origin cap because indels multiply the pairing work by up to
    // (2*max_gap_bp+1)^2 per origin, and a bound that cannot be seen is a bound nobody is holding.
    if (counters.latent_alignment_pairs > options.max_latent_alignment_pairs) {
        note_observed_work();
        throw ResourceExhausted(
            "fragment scoring reached its latent-alignment-pair cap; the declared gap width expands "
            "the pairing work beyond the cap and no partial ranking is returned");
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
