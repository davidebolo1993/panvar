#include "panvar/genotype_model.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <utility>

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
    std::set<std::string> template_names;
    for (const std::string& name : panel.template_names) {
        if (name.empty()) {
            throw std::invalid_argument("panel template name is empty");
        }
        if (!template_names.insert(name).second) {
            throw std::invalid_argument("duplicate panel template name: " + name);
        }
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

std::size_t PreparedGenotypeBlock::allele_id_for_source(const std::string& source_id) const {
    const auto found = allele_id_by_source.find(source_id);
    if (found == allele_id_by_source.end()) {
        throw std::out_of_range(
            "unknown source allele '" + source_id + "' in block '" + block_id + "'");
    }
    return found->second;
}

PreparedGenotypeModel::PreparedGenotypeModel(
    std::string locus_id,
    ChainOrientation chain_orientation,
    std::vector<PreparedGenotypeBlock> blocks,
    std::vector<std::string> invariant_segments,
    PanelAlleleMatrix panel,
    std::vector<ChainOrientation> template_source_orientations)
    : locus_id_(std::move(locus_id)),
      chain_orientation_(chain_orientation),
      blocks_(std::move(blocks)),
      invariant_segments_(std::move(invariant_segments)),
      panel_(std::move(panel)),
      template_source_orientations_(std::move(template_source_orientations)) {}

const CanonicalBlockAllele& PreparedGenotypeModel::allele(
    std::size_t block_id,
    std::size_t allele_id) const {
    if (block_id >= blocks_.size()) {
        throw std::out_of_range("prepared genotype block index is out of range");
    }
    if (allele_id >= blocks_[block_id].alleles.size()) {
        throw std::out_of_range(
            "prepared genotype allele index is out of range for block '" +
            blocks_[block_id].block_id + "'");
    }
    return blocks_[block_id].alleles[allele_id];
}

PreparedGenotypeModel prepare_genotype_model(const StructuralLocusInput& input) {
    if (input.locus_id.empty()) {
        throw std::invalid_argument("genotype locus ID is empty");
    }
    if (input.blocks.empty()) {
        throw std::invalid_argument("genotype locus has no blocks");
    }
    if (input.invariant_segments.size() != input.blocks.size() + 1) {
        throw std::invalid_argument(
            "genotype locus must have exactly one more invariant segment than blocks");
    }
    if (input.templates.empty()) {
        throw std::invalid_argument("panel has no templates");
    }

    std::set<std::string> block_ids;
    std::vector<PreparedGenotypeBlock> blocks;
    blocks.reserve(input.blocks.size());
    for (const GenotypeBlockInput& input_block : input.blocks) {
        if (input_block.block_id.empty()) {
            throw std::invalid_argument("genotype block ID is empty");
        }
        if (!block_ids.insert(input_block.block_id).second) {
            throw std::invalid_argument("duplicate genotype block ID: " + input_block.block_id);
        }
        if (input_block.alleles.empty()) {
            throw std::invalid_argument(
                "genotype block '" + input_block.block_id + "' has no alleles");
        }

        std::set<std::string> source_ids;
        std::map<std::string, std::vector<std::string>> sources_by_sequence;
        for (const BlockAlleleInput& input_allele : input_block.alleles) {
            if (input_allele.source_id.empty()) {
                throw std::invalid_argument(
                    "source allele ID is empty in block '" + input_block.block_id + "'");
            }
            if (!source_ids.insert(input_allele.source_id).second) {
                throw std::invalid_argument(
                    "duplicate source allele ID '" + input_allele.source_id +
                    "' in block '" + input_block.block_id + "'");
            }
            sources_by_sequence[input_allele.sequence].push_back(input_allele.source_id);
        }

        PreparedGenotypeBlock block;
        block.block_id = input_block.block_id;
        block.alleles.reserve(sources_by_sequence.size());
        for (auto& sequence_sources : sources_by_sequence) {
            std::sort(sequence_sources.second.begin(), sequence_sources.second.end());
            CanonicalBlockAllele allele;
            allele.id = block.alleles.size();
            allele.stable_id =
                input_block.block_id + ":allele:" + std::to_string(allele.id);
            allele.sequence = sequence_sources.first;
            allele.source_ids = std::move(sequence_sources.second);
            for (const std::string& source_id : allele.source_ids) {
                block.allele_id_by_source.emplace(source_id, allele.id);
            }
            block.alleles.push_back(std::move(allele));
        }
        blocks.push_back(std::move(block));
    }

    std::vector<const PanelTemplateInput*> ordered_templates;
    ordered_templates.reserve(input.templates.size());
    for (const PanelTemplateInput& input_template : input.templates) {
        ordered_templates.push_back(&input_template);
    }
    std::sort(
        ordered_templates.begin(),
        ordered_templates.end(),
        [](const PanelTemplateInput* lhs, const PanelTemplateInput* rhs) {
            return lhs->name < rhs->name;
        });

    PanelAlleleMatrix panel;
    panel.template_names.reserve(ordered_templates.size());
    panel.template_alleles.reserve(ordered_templates.size());
    std::vector<ChainOrientation> template_source_orientations;
    template_source_orientations.reserve(ordered_templates.size());
    std::string previous_name;
    for (const PanelTemplateInput* input_template : ordered_templates) {
        if (input_template->name.empty()) {
            throw std::invalid_argument("panel template name is empty");
        }
        if (!previous_name.empty() && input_template->name == previous_name) {
            throw std::invalid_argument("duplicate panel template name: " + input_template->name);
        }
        previous_name = input_template->name;
        if (input_template->frame_status != TemplateFrameStatus::Complete) {
            throw std::invalid_argument(
                "panel template '" + input_template->name + "' has a partial chain frame");
        }
        if (input_template->allele_source_ids.size() != blocks.size()) {
            throw std::invalid_argument(
                "panel template '" + input_template->name +
                "' does not map exactly one allele at every block");
        }

        std::vector<std::size_t> row;
        row.reserve(blocks.size());
        for (std::size_t block_id = 0; block_id < blocks.size(); ++block_id) {
            const std::string& source_id = input_template->allele_source_ids[block_id];
            const auto found = blocks[block_id].allele_id_by_source.find(source_id);
            if (found == blocks[block_id].allele_id_by_source.end()) {
                throw std::invalid_argument(
                    "panel template '" + input_template->name + "' maps unknown allele '" +
                    source_id + "' at block '" + blocks[block_id].block_id + "'");
            }
            row.push_back(found->second);
            ++blocks[block_id].alleles[found->second].panel_carrier_count;
        }
        panel.template_names.push_back(input_template->name);
        panel.template_alleles.push_back(std::move(row));
        template_source_orientations.push_back(input_template->source_orientation);
    }
    validate_panel(panel);

    return PreparedGenotypeModel(
        input.locus_id,
        input.chain_orientation,
        std::move(blocks),
        input.invariant_segments,
        std::move(panel),
        std::move(template_source_orientations));
}

} // namespace panvar
