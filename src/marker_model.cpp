#include "panvar/marker_model.hpp"

#include "panvar/syncmer.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
#include <set>
#include <stdexcept>
#include <tuple>

namespace panvar {
namespace {

using CountedMarkers = std::vector<std::pair<std::uint64_t, std::uint32_t>>;

struct VariationEvent {
    std::uint64_t code = 0;
    std::uint32_t block = 0;
};

struct FixedContribution {
    std::uint64_t code = 0;
    std::uint64_t count = 0;
};

std::uint64_t mix64(std::uint64_t value) {
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31);
}

class CompactCodeIndex {
public:
    explicit CompactCodeIndex(const std::vector<std::uint64_t>& codes) {
        std::size_t capacity = 1;
        while (capacity < codes.size() * 2 + 1) {
            if (capacity > std::numeric_limits<std::size_t>::max() / 2) {
                throw std::length_error("marker lookup capacity overflows size_t");
            }
            capacity *= 2;
        }
        keys_.assign(capacity, kEmpty);
        slots_.assign(capacity, 0);
        mask_ = capacity - 1;
        for (std::size_t slot = 0; slot < codes.size(); ++slot) {
            if (slot > std::numeric_limits<std::uint32_t>::max()) {
                throw std::length_error("marker universe exceeds uint32 slot range");
            }
            std::size_t at = static_cast<std::size_t>(mix64(codes[slot])) & mask_;
            while (keys_[at] != kEmpty) at = (at + 1) & mask_;
            keys_[at] = codes[slot];
            slots_[at] = static_cast<std::uint32_t>(slot);
        }
    }

    bool find(std::uint64_t code, std::uint32_t& slot) const noexcept {
        std::size_t at = static_cast<std::size_t>(mix64(code)) & mask_;
        while (keys_[at] != kEmpty) {
            if (keys_[at] == code) {
                slot = slots_[at];
                return true;
            }
            at = (at + 1) & mask_;
        }
        return false;
    }

private:
    static constexpr std::uint64_t kEmpty = std::numeric_limits<std::uint64_t>::max();
    std::vector<std::uint64_t> keys_;
    std::vector<std::uint32_t> slots_;
    std::size_t mask_ = 0;
};

CountedMarkers count_sequence_markers(
    const std::string& sequence,
    std::size_t kmer_size,
    std::size_t syncmer_s) {
    const std::vector<KmerOccurrence> occurrences =
        collect_syncmers(sequence, kmer_size, syncmer_s);
    std::vector<std::uint64_t> codes;
    codes.reserve(occurrences.size());
    for (const KmerOccurrence& occurrence : occurrences) codes.push_back(occurrence.code);
    std::sort(codes.begin(), codes.end());

    CountedMarkers result;
    for (std::size_t begin = 0; begin < codes.size();) {
        std::size_t end = begin + 1;
        while (end < codes.size() && codes[end] == codes[begin]) ++end;
        const std::size_t count = end - begin;
        if (count > std::numeric_limits<std::uint32_t>::max()) {
            throw std::length_error("one sequence has more than uint32 marker copies");
        }
        result.emplace_back(codes[begin], static_cast<std::uint32_t>(count));
        begin = end;
    }
    return result;
}

std::vector<FixedContribution> aggregate_fixed(std::vector<FixedContribution> values) {
    std::sort(values.begin(), values.end(), [](const auto& lhs, const auto& rhs) {
        return lhs.code < rhs.code;
    });
    std::vector<FixedContribution> result;
    for (const FixedContribution& value : values) {
        if (result.empty() || result.back().code != value.code) {
            result.push_back(value);
        } else {
            if (value.count > std::numeric_limits<std::uint64_t>::max() - result.back().count) {
                throw std::length_error("fixed marker multiplicity overflows uint64");
            }
            result.back().count += value.count;
        }
    }
    return result;
}

std::uint64_t fixed_count(
    const std::vector<FixedContribution>& fixed,
    std::uint64_t code) {
    const auto found = std::lower_bound(
        fixed.begin(), fixed.end(), code,
        [](const FixedContribution& value, std::uint64_t key) { return value.code < key; });
    return found != fixed.end() && found->code == code ? found->count : 0;
}

double poisson_log_kernel(std::uint64_t observed, double mean) {
    if (!(mean > 0.0) || !std::isfinite(mean)) {
        throw std::invalid_argument("marker Poisson mean must be finite and positive");
    }
    return static_cast<double>(observed) * std::log(mean) - mean;
}

double poisson_log_probability(std::uint64_t observed, double mean) {
    return poisson_log_kernel(observed, mean) -
           std::lgamma(static_cast<double>(observed) + 1.0);
}

} // namespace

double BlockProposalEmissions::at(std::size_t allele1, std::size_t allele2) const {
    if (allele1 >= allele_count || allele2 >= allele_count ||
        log_score.size() != allele_count * allele_count) {
        throw std::out_of_range("proposal emission allele pair is out of range");
    }
    return log_score[allele1 * allele_count + allele2];
}

void validate_marker_parameters(const MarkerModelParameters& parameters) {
    if (parameters.kmer_size == 0 || parameters.kmer_size > 31) {
        throw std::invalid_argument("marker k-mer size must be in [1,31]");
    }
    if (parameters.syncmer_s == 0 || parameters.syncmer_s > parameters.kmer_size) {
        throw std::invalid_argument("marker syncmer s must satisfy 1 <= s <= k");
    }
    if (!std::isfinite(parameters.expected_observations_per_copy) ||
        parameters.expected_observations_per_copy <= 0.0) {
        throw std::invalid_argument("marker expected observations per copy must be positive");
    }
    if (!std::isfinite(parameters.background_observations) ||
        parameters.background_observations <= 0.0) {
        throw std::invalid_argument("marker background observations must be positive");
    }
    if (parameters.provenance.empty()) {
        throw std::invalid_argument("marker model provenance must not be empty");
    }
}

PreparedMarkerPanel build_marker_panel(
    const PreparedGenotypeModel& model,
    std::size_t kmer_size,
    std::size_t syncmer_s) {
    if (kmer_size == 0 || kmer_size > 31 || syncmer_s == 0 || syncmer_s > kmer_size) {
        throw std::invalid_argument("marker parameters must satisfy 1 <= s <= k <= 31");
    }
    if (model.blocks().size() > std::numeric_limits<std::uint32_t>::max()) {
        throw std::length_error("genotype model has more than uint32 blocks");
    }

    std::vector<std::vector<CountedMarkers>> allele_counts(model.blocks().size());
    std::vector<VariationEvent> variations;
    std::vector<FixedContribution> fixed_values;

    for (std::size_t block_id = 0; block_id < model.blocks().size(); ++block_id) {
        const PreparedGenotypeBlock& block = model.blocks()[block_id];
        if (block.alleles.size() > std::numeric_limits<std::uint32_t>::max()) {
            throw std::length_error("genotype block has more than uint32 alleles");
        }
        allele_counts[block_id].reserve(block.alleles.size());
        for (std::size_t allele_id = 0; allele_id < block.alleles.size(); ++allele_id) {
            CountedMarkers counts = count_sequence_markers(
                block.alleles[allele_id].sequence, kmer_size, syncmer_s);
            allele_counts[block_id].push_back(std::move(counts));
        }

        // Each allele list is already sorted by code. Merge those lists directly rather than
        // copying every entry into one second N-element array and sorting it again. At C4 this keeps
        // construction proportional to total sparse markers * log(number of alleles), with only one
        // cursor per allele as additional storage.
        struct Cursor {
            std::uint64_t code = 0;
            std::uint32_t allele = 0;
            std::size_t offset = 0;
        };
        const auto later = [](const Cursor& lhs, const Cursor& rhs) {
            return std::tie(lhs.code, lhs.allele) > std::tie(rhs.code, rhs.allele);
        };
        std::priority_queue<Cursor, std::vector<Cursor>, decltype(later)> queue(later);
        for (std::size_t allele_id = 0; allele_id < allele_counts[block_id].size(); ++allele_id) {
            if (!allele_counts[block_id][allele_id].empty()) {
                queue.push({allele_counts[block_id][allele_id][0].first,
                            static_cast<std::uint32_t>(allele_id), 0});
            }
        }
        while (!queue.empty()) {
            const std::uint64_t code = queue.top().code;
            std::uint32_t min_count = std::numeric_limits<std::uint32_t>::max();
            std::uint32_t max_count = 0;
            std::size_t carriers = 0;
            while (!queue.empty() && queue.top().code == code) {
                Cursor cursor = queue.top();
                queue.pop();
                const CountedMarkers& source = allele_counts[block_id][cursor.allele];
                const std::uint32_t count = source[cursor.offset].second;
                min_count = std::min(min_count, count);
                max_count = std::max(max_count, count);
                ++carriers;
                ++cursor.offset;
                if (cursor.offset < source.size()) {
                    cursor.code = source[cursor.offset].first;
                    queue.push(cursor);
                }
            }
            if (carriers < block.alleles.size()) min_count = 0;
            if (min_count == max_count) {
                fixed_values.push_back({code, min_count});
            } else {
                variations.push_back({code, static_cast<std::uint32_t>(block_id)});
            }
        }
    }

    for (const std::string& context : model.invariant_segments()) {
        for (const auto& marker : count_sequence_markers(context, kmer_size, syncmer_s)) {
            fixed_values.push_back({marker.first, marker.second});
        }
    }
    const std::vector<FixedContribution> fixed = aggregate_fixed(std::move(fixed_values));

    std::sort(variations.begin(), variations.end(), [](const auto& lhs, const auto& rhs) {
        return std::tie(lhs.code, lhs.block) < std::tie(rhs.code, rhs.block);
    });
    std::vector<VariationEvent> retained;
    std::vector<std::uint64_t> variable_anywhere;
    std::size_t dropped_multiblock = 0;
    for (std::size_t begin = 0; begin < variations.size();) {
        std::size_t end = begin + 1;
        while (end < variations.size() && variations[end].code == variations[begin].code) ++end;
        variable_anywhere.push_back(variations[begin].code);
        if (end - begin == 1) {
            retained.push_back(variations[begin]);
        } else {
            ++dropped_multiblock;
        }
        begin = end;
    }

    std::vector<std::uint64_t> anchor_codes;
    anchor_codes.reserve(fixed.size());
    for (const FixedContribution& value : fixed) {
        if (!std::binary_search(variable_anywhere.begin(), variable_anywhere.end(), value.code)) {
            anchor_codes.push_back(value.code);
        }
    }
    std::vector<std::uint64_t> codes = anchor_codes;
    codes.reserve(anchor_codes.size() + retained.size());
    for (const VariationEvent& value : retained) codes.push_back(value.code);
    std::sort(codes.begin(), codes.end());
    codes.erase(std::unique(codes.begin(), codes.end()), codes.end());
    if (codes.size() > std::numeric_limits<std::uint32_t>::max()) {
        throw std::length_error("marker universe exceeds uint32 slot range");
    }

    PreparedMarkerPanel panel;
    panel.kmer_size = kmer_size;
    panel.syncmer_s = syncmer_s;
    panel.codes = std::move(codes);
    panel.dropped_multiblock_markers = dropped_multiblock;
    panel.blocks.resize(model.blocks().size());
    for (std::size_t block_id = 0; block_id < model.blocks().size(); ++block_id) {
        BlockMarkerInventory& out = panel.blocks[block_id];
        out.block_id = model.blocks()[block_id].block_id;
        out.allele_count = model.blocks()[block_id].alleles.size();
        out.allele_markers.resize(out.allele_count);
        panel.total_alleles += out.allele_count;
        panel.scored_alleles += out.allele_count;

        std::vector<std::uint64_t> local_codes;
        for (const VariationEvent& value : retained) {
            if (value.block == block_id) local_codes.push_back(value.code);
        }
        for (const std::uint64_t code : local_codes) {
            const auto global = std::lower_bound(panel.codes.begin(), panel.codes.end(), code);
            out.marker_slots.push_back(
                static_cast<std::uint32_t>(global - panel.codes.begin()));
            const std::uint64_t copies = fixed_count(fixed, code);
            if (copies > std::numeric_limits<std::uint32_t>::max()) {
                throw std::length_error("fixed haploid marker multiplicity exceeds uint32");
            }
            out.fixed_haploid_multiplicity.push_back(static_cast<std::uint32_t>(copies));
        }
        for (std::size_t allele_id = 0; allele_id < out.allele_count; ++allele_id) {
            const CountedMarkers& raw = allele_counts[block_id][allele_id];
            std::size_t local = 0;
            std::size_t source = 0;
            while (local < local_codes.size() && source < raw.size()) {
                if (local_codes[local] < raw[source].first) {
                    ++local;
                } else if (raw[source].first < local_codes[local]) {
                    ++source;
                } else {
                    out.allele_markers[allele_id].emplace_back(
                        static_cast<std::uint32_t>(local), raw[source].second);
                    ++local;
                    ++source;
                }
            }
        }
    }

    for (const std::uint64_t code : anchor_codes) {
        const auto global = std::lower_bound(panel.codes.begin(), panel.codes.end(), code);
        panel.anchor_slots.push_back(static_cast<std::uint32_t>(global - panel.codes.begin()));
        const std::uint64_t copies = fixed_count(fixed, code);
        if (copies > std::numeric_limits<std::uint32_t>::max()) {
            throw std::length_error("anchor haploid marker multiplicity exceeds uint32");
        }
        panel.anchor_haploid_multiplicity.push_back(static_cast<std::uint32_t>(copies));
    }
    return panel;
}

MarkerObservations count_marker_observations(
    const PreparedMarkerPanel& panel,
    const std::vector<ReadPair>& fragments) {
    if (panel.kmer_size == 0 || panel.kmer_size > 31 || panel.syncmer_s == 0 ||
        panel.syncmer_s > panel.kmer_size) {
        throw std::invalid_argument("prepared marker panel has invalid syncmer parameters");
    }
    const CompactCodeIndex index(panel.codes);
    MarkerObservations result;
    result.counts.assign(panel.codes.size(), 0);
    for (const ReadPair& fragment : fragments) {
        validate_read_pair(fragment);
        for (const SequencingRead* read : {&fragment.first, &fragment.second}) {
            ++result.reads;
            result.bases += read->bases.size();
            const std::vector<KmerOccurrence> occurrences =
                collect_syncmers(read->bases, panel.kmer_size, panel.syncmer_s);
            result.syncmer_occurrences += occurrences.size();
            for (const KmerOccurrence& occurrence : occurrences) {
                std::uint32_t slot = 0;
                if (!index.find(occurrence.code, slot)) continue;
                if (result.counts[slot] == std::numeric_limits<std::uint64_t>::max()) {
                    throw std::length_error("marker observation count overflows uint64");
                }
                ++result.counts[slot];
                ++result.matched_occurrences;
            }
        }
    }
    return result;
}

PreparedProposalEmissions build_proposal_emissions(
    const PreparedMarkerPanel& panel,
    const MarkerObservations& observations,
    const MarkerModelParameters& parameters) {
    validate_marker_parameters(parameters);
    if (parameters.kmer_size != panel.kmer_size || parameters.syncmer_s != panel.syncmer_s) {
        throw std::invalid_argument("marker parameters do not match the prepared panel");
    }
    if (observations.counts.size() != panel.codes.size()) {
        throw std::invalid_argument("marker observations and panel universe differ in size");
    }

    PreparedProposalEmissions result;
    result.total_alleles = panel.total_alleles;
    result.scored_alleles = panel.scored_alleles;
    result.blocks.reserve(panel.blocks.size());
    for (const BlockMarkerInventory& block : panel.blocks) {
        if (block.marker_slots.size() != block.fixed_haploid_multiplicity.size() ||
            block.allele_markers.size() != block.allele_count) {
            throw std::invalid_argument("prepared block marker inventory is inconsistent");
        }
        BlockProposalEmissions emission;
        emission.allele_count = block.allele_count;
        emission.scored_alleles = block.allele_count;
        emission.log_score.assign(block.allele_count * block.allele_count, 0.0);

        std::vector<double> baseline_kernel(block.marker_slots.size(), 0.0);
        for (std::size_t marker = 0; marker < block.marker_slots.size(); ++marker) {
            const std::uint32_t slot = block.marker_slots[marker];
            if (slot >= observations.counts.size()) {
                throw std::invalid_argument("block marker slot is outside the observation universe");
            }
            const double baseline_mean = parameters.background_observations +
                parameters.expected_observations_per_copy * 2.0 *
                    static_cast<double>(block.fixed_haploid_multiplicity[marker]);
            baseline_kernel[marker] = poisson_log_kernel(observations.counts[slot], baseline_mean);
            emission.candidate_independent_log_constant +=
                poisson_log_probability(observations.counts[slot], baseline_mean);
        }

        for (std::size_t allele1 = 0; allele1 < block.allele_count; ++allele1) {
            for (std::size_t allele2 = allele1; allele2 < block.allele_count; ++allele2) {
                const auto& first = block.allele_markers[allele1];
                const auto& second = block.allele_markers[allele2];
                std::size_t i = 0;
                std::size_t j = 0;
                double score = 0.0;
                while (i < first.size() || j < second.size()) {
                    std::uint32_t marker = 0;
                    std::uint64_t copies = 0;
                    if (j == second.size() ||
                        (i < first.size() && first[i].first < second[j].first)) {
                        marker = first[i].first;
                        copies = first[i].second;
                        ++i;
                    } else if (i == first.size() || second[j].first < first[i].first) {
                        marker = second[j].first;
                        copies = second[j].second;
                        ++j;
                    } else {
                        marker = first[i].first;
                        copies = static_cast<std::uint64_t>(first[i].second) + second[j].second;
                        ++i;
                        ++j;
                    }
                    if (marker >= block.marker_slots.size()) {
                        throw std::invalid_argument("allele marker index is outside its block");
                    }
                    const std::uint32_t slot = block.marker_slots[marker];
                    const double mean = parameters.background_observations +
                        parameters.expected_observations_per_copy *
                            (2.0 * static_cast<double>(block.fixed_haploid_multiplicity[marker]) +
                             static_cast<double>(copies));
                    score += poisson_log_kernel(observations.counts[slot], mean) -
                             baseline_kernel[marker];
                }
                emission.log_score[allele1 * block.allele_count + allele2] = score;
                emission.log_score[allele2 * block.allele_count + allele1] = score;
            }
        }
        result.blocks.push_back(std::move(emission));
    }
    return result;
}

MarkerDepthEstimate estimate_marker_depth(
    const PreparedMarkerPanel& panel,
    const MarkerObservations& observations,
    std::optional<double> declared_background) {
    if (panel.anchor_slots.size() != panel.anchor_haploid_multiplicity.size()) {
        throw std::invalid_argument("anchor slots and multiplicities differ in size");
    }
    if (observations.counts.size() != panel.codes.size()) {
        throw std::invalid_argument("marker observations and panel universe differ in size");
    }
    std::vector<double> counts;
    std::vector<double> copies;   // diploid: 2 * haploid multiplicity
    std::set<std::uint32_t> multiplicities;
    for (std::size_t k = 0; k < panel.anchor_slots.size(); ++k) {
        const std::uint32_t haploid = panel.anchor_haploid_multiplicity[k];
        if (haploid == 0) continue;
        if (panel.anchor_slots[k] >= observations.counts.size()) {
            throw std::invalid_argument("anchor slot is outside the observation universe");
        }
        counts.push_back(static_cast<double>(observations.counts[panel.anchor_slots[k]]));
        copies.push_back(2.0 * static_cast<double>(haploid));
        multiplicities.insert(haploid);
    }
    MarkerDepthEstimate estimate;
    estimate.anchors_used = counts.size();
    estimate.distinct_multiplicities = multiplicities.size();
    if (counts.empty()) {
        throw std::invalid_argument("no invariant anchors are available to estimate marker depth");
    }
    double count_sum = 0.0;
    double copy_sum = 0.0;
    for (std::size_t k = 0; k < counts.size(); ++k) {
        count_sum += counts[k];
        copy_sum += copies[k];
    }
    if (!(count_sum > 0.0)) {
        throw std::invalid_argument(
            "no read observation matched any invariant anchor, so there is no depth to estimate");
    }

    if (declared_background.has_value()) {
        const double background = *declared_background;
        if (!std::isfinite(background) || background < 0.0) {
            throw std::invalid_argument("declared marker background must be finite and non-negative");
        }
        estimate.background = background;
        estimate.background_estimated = false;
        if (background == 0.0) {
            estimate.rate_per_copy = count_sum / copy_sum;
            estimate.converged = true;
            estimate.termination = "closed_form";
            estimate.provenance = "anchor Poisson MLE, closed form, background declared 0";
            return estimate;
        }
        // d/dr log L = sum c_k d_k / (bg + r d_k) - sum d_k, strictly decreasing in r.
        const auto score = [&](double rate) {
            double value = -copy_sum;
            for (std::size_t k = 0; k < counts.size(); ++k) {
                value += counts[k] * copies[k] / (background + rate * copies[k]);
            }
            return value;
        };
        if (score(0.0) <= 0.0) {
            throw std::invalid_argument(
                "invariant anchors carry no signal above the declared background");
        }
        double lo = 0.0;
        double hi = std::max(1.0, count_sum / copy_sum);
        while (score(hi) > 0.0) hi *= 2.0;
        for (int i = 0; i < 200 && hi - lo > 1e-14 * hi; ++i) {
            const double mid = 0.5 * (lo + hi);
            (score(mid) > 0.0 ? lo : hi) = mid;
            estimate.iterations = static_cast<std::size_t>(i + 1);
        }
        estimate.rate_per_copy = 0.5 * (lo + hi);
        estimate.converged = true;
        estimate.termination = "declared_background_root";
        estimate.provenance = "anchor Poisson MLE, one-dimensional, background declared";
        return estimate;
    }

    if (multiplicities.size() < 2) {
        throw std::invalid_argument(
            "every invariant anchor has the same multiplicity, so background and per-copy rate are "
            "not separately identifiable; declare the marker background");
    }
    // Joint maximum likelihood by damped Newton. The log-likelihood
    //   sum_k c_k log(bg + r d_k) - N bg - r D
    // is concave in (bg, r), so Newton with a backtracking line search reaches the unique maximum
    // quadratically. EM on the same model is only linearly convergent and did not reach tolerance on
    // C4's 25,167 anchors in 20,000 iterations.
    const double n = static_cast<double>(counts.size());
    const double background_floor = 1e-9;
    const auto log_likelihood = [&](double bg, double r) {
        double value = -n * bg - r * copy_sum;
        for (std::size_t k = 0; k < counts.size(); ++k) {
            if (counts[k] > 0.0) value += counts[k] * std::log(bg + r * copies[k]);
        }
        return value;
    };
    // Start from least squares on c = bg + r d, pulled into the interior.
    const double mean_copies = copy_sum / n;
    const double mean_count = count_sum / n;
    double sxx = 0.0;
    double sxy = 0.0;
    for (std::size_t k = 0; k < counts.size(); ++k) {
        sxx += (copies[k] - mean_copies) * (copies[k] - mean_copies);
        sxy += (copies[k] - mean_copies) * (counts[k] - mean_count);
    }
    double rate = sxx > 0.0 ? sxy / sxx : mean_count / mean_copies;
    if (!(rate > 0.0)) rate = mean_count / mean_copies;
    double background = mean_count - rate * mean_copies;
    if (!(background > background_floor)) background = std::max(10.0 * background_floor, 0.01 * mean_count);
    double current = log_likelihood(background, rate);

    const auto derivatives = [&](double bg, double r, double& g_bg, double& g_rate, double& h_bb,
                                 double& h_br, double& h_rr) {
        g_bg = -n;
        g_rate = -copy_sum;
        h_bb = 0.0;
        h_br = 0.0;
        h_rr = 0.0;
        for (std::size_t k = 0; k < counts.size(); ++k) {
            const double mean = bg + r * copies[k];
            const double first = counts[k] / mean;
            const double second = counts[k] / (mean * mean);
            g_bg += first;
            g_rate += first * copies[k];
            h_bb -= second;
            h_br -= second * copies[k];
            h_rr -= second * copies[k] * copies[k];
        }
    };
    // The rate that maximises the likelihood with the background held fixed: the rate score is
    // strictly decreasing, so its root is found by bisection.
    const auto rate_given_background = [&](double bg) {
        const auto score = [&](double r) {
            double value = -copy_sum;
            for (std::size_t k = 0; k < counts.size(); ++k) value += counts[k] * copies[k] / (bg + r * copies[k]);
            return value;
        };
        double lo = 0.0;
        double hi = std::max(1.0, count_sum / copy_sum);
        while (score(hi) > 0.0) hi *= 2.0;
        for (int i = 0; i < 300 && hi - lo > 1e-15 * hi; ++i) {
            const double mid = 0.5 * (lo + hi);
            (score(mid) > 0.0 ? lo : hi) = mid;
        }
        return 0.5 * (lo + hi);
    };

    estimate.termination = "iteration_cap";
    double g_bg = 0.0;
    double g_rate = 0.0;
    double h_bb = 0.0;
    double h_br = 0.0;
    double h_rr = 0.0;
    for (std::size_t iteration = 1; iteration <= 200; ++iteration) {
        derivatives(background, rate, g_bg, g_rate, h_bb, h_br, h_rr);
        estimate.iterations = iteration;
        if (std::abs(g_bg) <= 1e-9 * n && std::abs(g_rate) <= 1e-9 * copy_sum) {
            estimate.converged = true;
            estimate.termination = "interior_stationary";
            break;
        }
        // Active constraint. When the background sits on its floor and still wants to fall, the
        // maximum is on that face: maximise over the rate alone there, then check the Karush-Kuhn-
        // Tucker condition. Projecting a joint Newton step instead never makes the rate stationary --
        // on NA19240 that stopped at the iteration cap with a rate 3.4% high.
        if (background <= background_floor && g_bg <= 0.0) {
            background = background_floor;
            rate = rate_given_background(background_floor);
            current = log_likelihood(background, rate);
            derivatives(background, rate, g_bg, g_rate, h_bb, h_br, h_rr);
            if (g_bg <= 0.0 && std::abs(g_rate) <= 1e-9 * copy_sum) {
                estimate.converged = true;
                estimate.background_at_floor = true;
                estimate.termination = "boundary_kkt";
                break;
            }
        }
        const double determinant = h_bb * h_rr - h_br * h_br;
        double step_bg = 0.0;
        double step_rate = 0.0;
        if (determinant > 0.0 && h_bb < 0.0) {
            step_bg = -(h_rr * g_bg - h_br * g_rate) / determinant;
            step_rate = -(-h_br * g_bg + h_bb * g_rate) / determinant;
        } else {
            step_bg = g_bg / std::max(1.0, -h_bb);
            step_rate = g_rate / std::max(1.0, -h_rr);
        }
        bool moved = false;
        double scale = 1.0;
        for (int halving = 0; halving < 60; ++halving, scale *= 0.5) {
            const double next_bg = std::max(background_floor, background + scale * step_bg);
            const double next_rate = rate + scale * step_rate;
            if (!(next_rate > 0.0)) continue;
            const double candidate = log_likelihood(next_bg, next_rate);
            if (candidate > current) {
                background = next_bg;
                rate = next_rate;
                current = candidate;
                moved = true;
                break;
            }
        }
        if (!moved) {
            estimate.termination = "no_ascent";
            break;
        }
    }
    derivatives(background, rate, g_bg, g_rate, h_bb, h_br, h_rr);
    estimate.gradient_background = g_bg;
    estimate.gradient_rate = g_rate;
    estimate.log_likelihood = log_likelihood(background, rate);
    estimate.background = background;
    estimate.rate_per_copy = rate;
    estimate.background_estimated = true;
    estimate.provenance = "anchor Poisson MLE, background and rate jointly by active-set Newton";
    return estimate;
}

} // namespace panvar
