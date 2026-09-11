#include "panvar/ls_hmm.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace panvar {
namespace {

constexpr double kNegInf = -std::numeric_limits<double>::infinity();

double log_add_exp(double a, double b) {
    if (a == kNegInf) return b;
    if (b == kNegInf) return a;
    const double hi = std::max(a, b);
    const double lo = std::min(a, b);
    return hi + std::log1p(std::exp(lo - hi));
}

double log_sum_exp(const std::vector<double>& values) {
    double result = kNegInf;
    for (double value : values) result = log_add_exp(result, value);
    return result;
}

void validate_path_inputs(
    const PanelAlleleMatrix& panel,
    const MosaicPath& mosaic,
    const std::vector<double>& switch_probability) {
    validate_panel(panel);
    if (mosaic.alleles.size() != panel.block_count()) {
        throw std::invalid_argument("mosaic and panel have different block counts");
    }
    if (switch_probability.size() + 1 != panel.block_count()) {
        throw std::invalid_argument("one switch probability is required per block boundary");
    }
    for (double r : switch_probability) {
        if (!std::isfinite(r) || r < 0.0 || r > 1.0) {
            throw std::invalid_argument("switch probabilities must be finite and in [0,1]");
        }
    }
}

} // namespace

double log_haploid_mosaic_prior(
    const PanelAlleleMatrix& panel,
    const MosaicPath& mosaic,
    const std::vector<double>& switch_probability) {
    validate_path_inputs(panel, mosaic, switch_probability);

    const std::size_t templates = panel.template_count();
    const double log_templates = std::log(static_cast<double>(templates));
    std::vector<double> previous(templates, kNegInf);
    std::vector<double> next(templates, kNegInf);

    for (std::size_t h = 0; h < templates; ++h) {
        if (panel.template_alleles[h][0] == mosaic.alleles[0]) {
            previous[h] = -log_templates;
        }
    }

    for (std::size_t block = 1; block < panel.block_count(); ++block) {
        const double r = switch_probability[block - 1];
        const double total = log_sum_exp(previous);
        const double log_stay = r < 1.0 ? std::log1p(-r) : kNegInf;
        const double log_switch = r > 0.0 ? std::log(r) - log_templates : kNegInf;

        std::fill(next.begin(), next.end(), kNegInf);
        for (std::size_t h = 0; h < templates; ++h) {
            if (panel.template_alleles[h][block] != mosaic.alleles[block]) continue;
            const double stay = previous[h] == kNegInf || log_stay == kNegInf
                                    ? kNegInf
                                    : previous[h] + log_stay;
            const double switched = total == kNegInf || log_switch == kNegInf
                                        ? kNegInf
                                        : total + log_switch;
            next[h] = log_add_exp(stay, switched);
        }
        previous.swap(next);
    }
    return log_sum_exp(previous);
}

double log_unordered_diplotype_prior(
    const PanelAlleleMatrix& panel,
    const MosaicDiplotype& diplotype,
    const std::vector<double>& switch_probability) {
    const double first =
        log_haploid_mosaic_prior(panel, diplotype.first, switch_probability);
    const double second =
        log_haploid_mosaic_prior(panel, diplotype.second, switch_probability);
    if (first == kNegInf || second == kNegInf) return kNegInf;
    return first + second + (diplotype.first != diplotype.second ? std::log(2.0) : 0.0);
}

} // namespace panvar
