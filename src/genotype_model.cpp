#include "panvar/genotype_model.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace panvar {
namespace {

double log_add_exp(double a, double b) {
    if (a == -std::numeric_limits<double>::infinity()) return b;
    if (b == -std::numeric_limits<double>::infinity()) return a;
    const double hi = std::max(a, b);
    const double lo = std::min(a, b);
    return hi + std::log1p(std::exp(lo - hi));
}

double gaussian_log_weight(std::size_t insert_bp, const LibraryModel& model) {
    const double z = (static_cast<double>(insert_bp) - model.insert_mean_bp) /
                     model.insert_sd_bp;
    return -0.5 * z * z;
}

void validate_read(const SequencingRead& read, const char* which) {
    if (read.bases.empty()) {
        throw std::invalid_argument(std::string(which) + " read is empty");
    }
    if (!read.qualities.empty() && read.qualities.size() != read.bases.size()) {
        throw std::invalid_argument(std::string(which) +
                                    " read has a different number of bases and qualities");
    }
}

} // namespace

void validate_read_pair(const ReadPair& fragment) {
    validate_read(fragment.first, "first");
    validate_read(fragment.second, "second");
}

void validate_library_model(const LibraryModel& model) {
    if (!std::isfinite(model.insert_mean_bp) || model.insert_mean_bp <= 0.0) {
        throw std::invalid_argument("insert mean must be finite and positive");
    }
    if (!std::isfinite(model.insert_sd_bp) || model.insert_sd_bp <= 0.0) {
        throw std::invalid_argument("insert standard deviation must be finite and positive");
    }
    if (!std::isfinite(model.support_sigmas) || model.support_sigmas <= 0.0) {
        throw std::invalid_argument("insert support sigmas must be finite and positive");
    }
    if (model.provenance.empty()) {
        throw std::invalid_argument("library model provenance must not be empty");
    }
    const double hi = model.insert_mean_bp + model.support_sigmas * model.insert_sd_bp;
    if (!std::isfinite(hi) || hi > static_cast<double>(std::numeric_limits<std::size_t>::max())) {
        throw std::invalid_argument("insert support exceeds the addressable size range");
    }
}

InsertSupport fragment_insert_support(
    const ReadPair& fragment,
    const LibraryModel& model) {
    validate_read_pair(fragment);
    validate_library_model(model);

    const std::size_t physical_floor =
        std::max(fragment.first.bases.size(), fragment.second.bases.size());
    const double raw_lo = std::ceil(
        model.insert_mean_bp - model.support_sigmas * model.insert_sd_bp);
    const double raw_hi = std::floor(
        model.insert_mean_bp + model.support_sigmas * model.insert_sd_bp);

    const std::size_t statistical_floor =
        raw_lo <= 1.0 ? 1 : static_cast<std::size_t>(raw_lo);
    const std::size_t hi = raw_hi < 1.0 ? 0 : static_cast<std::size_t>(raw_hi);
    return {std::max(physical_floor, statistical_floor), hi};
}

double DiscreteInsertDistribution::log_probability(std::size_t insert_bp) const noexcept {
    if (support.empty() || insert_bp < support.min_bp || insert_bp > support.max_bp) {
        return -std::numeric_limits<double>::infinity();
    }
    return log_probability_by_offset[insert_bp - support.min_bp];
}

DiscreteInsertDistribution build_fragment_insert_distribution(
    const ReadPair& fragment,
    const LibraryModel& model) {
    const InsertSupport support = fragment_insert_support(fragment, model);
    if (support.empty()) return {support, {}};

    double log_z = -std::numeric_limits<double>::infinity();
    for (std::size_t length = support.min_bp;; ++length) {
        log_z = log_add_exp(log_z, gaussian_log_weight(length, model));
        if (length == support.max_bp) break;
    }
    DiscreteInsertDistribution result;
    result.support = support;
    result.log_probability_by_offset.reserve(support.max_bp - support.min_bp + 1);
    for (std::size_t length = support.min_bp;; ++length) {
        result.log_probability_by_offset.push_back(gaussian_log_weight(length, model) - log_z);
        if (length == support.max_bp) break;
    }
    return result;
}

void validate_panel(const PanelAlleleMatrix& panel) {
    if (panel.template_alleles.empty()) {
        throw std::invalid_argument("panel has no templates");
    }
    if (panel.template_names.size() != panel.template_alleles.size()) {
        throw std::invalid_argument("panel template names and allele rows differ in size");
    }
    const std::size_t blocks = panel.template_alleles.front().size();
    if (blocks == 0) {
        throw std::invalid_argument("panel has no blocks");
    }
    for (const auto& row : panel.template_alleles) {
        if (row.size() != blocks) {
            throw std::invalid_argument("panel allele matrix is not rectangular");
        }
    }
}

} // namespace panvar
