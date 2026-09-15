#include "panvar/ls_hmm.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>

namespace panvar {
namespace {

constexpr double kNegInf = -std::numeric_limits<double>::infinity();

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

namespace {

bool matches_at(const std::string& sequence, std::size_t offset, const std::string& piece) {
    return offset <= sequence.size() && sequence.size() - offset >= piece.size() &&
           sequence.compare(offset, piece.size(), piece) == 0;
}

} // namespace

double log_haploid_sequence_prior(
    const PreparedGenotypeModel& model,
    const std::string& sequence,
    const std::vector<double>& switch_probability) {
    const PanelAlleleMatrix& panel = model.panel();
    validate_panel(panel);
    const std::size_t blocks = model.blocks().size();
    if (blocks != panel.block_count()) {
        throw std::logic_error("prepared model blocks and panel columns differ");
    }
    if (switch_probability.size() + 1 != blocks) {
        throw std::invalid_argument("one switch probability is required per block boundary");
    }
    for (double r : switch_probability) {
        if (!std::isfinite(r) || r < 0.0 || r > 1.0) {
            throw std::invalid_argument("switch probabilities must be finite and in [0,1]");
        }
    }
    const std::vector<std::string>& contexts = model.invariant_segments();
    if (contexts.size() != blocks + 1) {
        throw std::logic_error("prepared model has an invalid invariant-segment frame");
    }

    const std::size_t templates = panel.template_count();
    const double log_templates = std::log(static_cast<double>(templates));
    if (!matches_at(sequence, 0, contexts[0])) return kNegInf;

    // frontier[offset][h]: log of the summed prior of every path prefix that spells sequence[0,offset)
    // through the current block and sits on template h there. Different prefixes reaching the same
    // offset are different allele paths spelling the same bytes; summing them is the aggregation.
    std::map<std::size_t, std::vector<double>> frontier;
    for (std::size_t b = 0; b < blocks; ++b) {
        std::map<std::size_t, std::vector<double>> next;
        const PreparedGenotypeBlock& block = model.blocks()[b];

        std::map<std::size_t, std::vector<double>> sources;
        if (b == 0) {
            sources.emplace(contexts[0].size(), std::vector<double>());
        } else {
            sources.swap(frontier);
        }
        const double r = b == 0 ? 0.0 : switch_probability[b - 1];
        const double log_stay = r < 1.0 ? std::log1p(-r) : kNegInf;
        const double log_switch = r > 0.0 ? std::log(r) - log_templates : kNegInf;

        for (const auto& source : sources) {
            const std::size_t offset = source.first;
            const double total = b == 0 ? kNegInf : log_sum_exp(source.second);
            for (std::size_t a = 0; a < block.alleles.size(); ++a) {
                const std::string& allele = block.alleles[a].sequence;
                if (!matches_at(sequence, offset, allele)) continue;
                const std::size_t after_allele = offset + allele.size();
                if (!matches_at(sequence, after_allele, contexts[b + 1])) continue;
                const std::size_t after = after_allele + contexts[b + 1].size();

                std::vector<double>& slot = next[after];
                if (slot.empty()) slot.assign(templates, kNegInf);
                for (std::size_t h = 0; h < templates; ++h) {
                    if (panel.template_alleles[h][b] != a) continue;
                    double value = kNegInf;
                    if (b == 0) {
                        value = -log_templates;
                    } else {
                        const double stay = source.second[h] == kNegInf || log_stay == kNegInf
                                                ? kNegInf
                                                : source.second[h] + log_stay;
                        const double switched = total == kNegInf || log_switch == kNegInf
                                                    ? kNegInf
                                                    : total + log_switch;
                        value = log_add_exp(stay, switched);
                    }
                    slot[h] = log_add_exp(slot[h], value);
                }
            }
        }
        frontier.swap(next);
        if (frontier.empty()) return kNegInf;
    }
    const auto complete = frontier.find(sequence.size());
    return complete == frontier.end() ? kNegInf : log_sum_exp(complete->second);
}

double log_unordered_sequence_pair_prior(
    const PreparedGenotypeModel& model,
    const std::string& first,
    const std::string& second,
    const std::vector<double>& switch_probability) {
    const double a = log_haploid_sequence_prior(model, first, switch_probability);
    if (a == kNegInf) return kNegInf;
    const double b = first == second ? a : log_haploid_sequence_prior(model, second, switch_probability);
    if (b == kNegInf) return kNegInf;
    return a + b + (first != second ? std::log(2.0) : 0.0);
}

} // namespace panvar
