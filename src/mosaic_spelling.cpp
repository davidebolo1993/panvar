#include "panvar/mosaic_spelling.hpp"

#include "panvar/graph_utils.hpp"

#include <algorithm>
#include <limits>
#include <map>
#include <stdexcept>

namespace panvar {
namespace {

void checked_add(std::size_t& total, std::size_t amount) {
    if (amount > std::numeric_limits<std::size_t>::max() - total) {
        throw std::length_error("mosaic sequence length overflows size_t");
    }
    total += amount;
}

std::size_t first_difference(const std::string& lhs, const std::string& rhs) {
    const std::size_t common = std::min(lhs.size(), rhs.size());
    std::size_t at = 0;
    while (at < common && lhs[at] == rhs[at]) ++at;
    return at;
}

} // namespace

HaplotypeSequence spell_mosaic(
    const PreparedGenotypeModel& model,
    const MosaicPath& mosaic) {
    const std::size_t block_count = model.blocks().size();
    if (mosaic.alleles.size() != block_count) {
        throw std::invalid_argument(
            "mosaic path does not contain exactly one allele for every prepared block");
    }
    if (model.invariant_segments().size() != block_count + 1) {
        throw std::logic_error("prepared model has an invalid invariant-segment frame");
    }

    std::size_t sequence_size = 0;
    for (const std::string& context : model.invariant_segments()) {
        checked_add(sequence_size, context.size());
    }
    for (std::size_t block_id = 0; block_id < block_count; ++block_id) {
        if (mosaic.alleles[block_id] >= model.blocks()[block_id].alleles.size()) {
            throw std::invalid_argument(
                "mosaic allele index is out of range at block '" +
                model.blocks()[block_id].block_id + "'");
        }
        checked_add(sequence_size, model.allele(block_id, mosaic.alleles[block_id]).sequence.size());
    }

    HaplotypeSequence result;
    result.mosaic = mosaic;
    result.bases.reserve(sequence_size);
    result.block_intervals.reserve(block_count);
    result.invariant_intervals.reserve(block_count + 1);

    for (std::size_t block_id = 0; block_id < block_count; ++block_id) {
        const std::size_t context_begin = result.bases.size();
        result.bases += model.invariant_segments()[block_id];
        result.invariant_intervals.push_back({context_begin, result.bases.size()});

        const std::size_t block_begin = result.bases.size();
        result.bases += model.allele(block_id, mosaic.alleles[block_id]).sequence;
        result.block_intervals.push_back({block_begin, result.bases.size()});
    }
    const std::size_t context_begin = result.bases.size();
    result.bases += model.invariant_segments().back();
    result.invariant_intervals.push_back({context_begin, result.bases.size()});

    if (result.bases.size() != sequence_size) {
        throw std::logic_error("mosaic spelling length disagrees with its checked size");
    }
    return result;
}

HaplotypeSequence spell_panel_template(
    const PreparedGenotypeModel& model,
    std::size_t template_id) {
    if (template_id >= model.panel().template_count()) {
        throw std::out_of_range("panel template index is out of range");
    }
    return spell_mosaic(model, MosaicPath{model.panel().template_alleles[template_id]});
}

std::string sequence_in_source_frame(
    const std::string& chain_bases,
    ChainOrientation source_orientation) {
    return source_orientation == ChainOrientation::Forward
        ? chain_bases
        : reverse_complement(chain_bases);
}

SequenceInterval interval_in_source_frame(
    SequenceInterval chain_interval,
    std::size_t sequence_size,
    ChainOrientation source_orientation) {
    if (chain_interval.begin > chain_interval.end || chain_interval.end > sequence_size) {
        throw std::out_of_range("chain interval lies outside the spelled sequence");
    }
    if (source_orientation == ChainOrientation::Forward) return chain_interval;
    return {sequence_size - chain_interval.end, sequence_size - chain_interval.begin};
}

void verify_complete_panel_round_trip(
    const PreparedGenotypeModel& model,
    const std::vector<AuthoritativePanelSequence>& authoritative_sequences) {
    if (model.frame_status() != TemplateFrameStatus::Complete) {
        throw std::invalid_argument("cannot verify a partial prepared panel frame");
    }
    if (model.template_source_orientations().size() != model.panel().template_count()) {
        throw std::logic_error("prepared template orientations and panel rows differ in size");
    }

    std::map<std::string, const std::string*> authoritative_by_name;
    for (const AuthoritativePanelSequence& sequence : authoritative_sequences) {
        if (sequence.template_name.empty()) {
            throw std::invalid_argument("authoritative panel sequence has an empty template name");
        }
        if (!authoritative_by_name.emplace(sequence.template_name, &sequence.source_frame_bases)
                 .second) {
            throw std::invalid_argument(
                "duplicate authoritative panel sequence: " + sequence.template_name);
        }
    }
    if (authoritative_by_name.size() != model.panel().template_count()) {
        throw std::invalid_argument(
            "authoritative sequence set and prepared panel differ in template count");
    }

    for (std::size_t template_id = 0; template_id < model.panel().template_count(); ++template_id) {
        const std::string& name = model.panel().template_names[template_id];
        const auto found = authoritative_by_name.find(name);
        if (found == authoritative_by_name.end()) {
            throw std::invalid_argument(
                "authoritative sequence is missing panel template '" + name + "'");
        }
        const HaplotypeSequence spelled = spell_panel_template(model, template_id);
        const std::string source_frame = sequence_in_source_frame(
            spelled.bases,
            model.template_source_orientations()[template_id]);
        if (source_frame != *found->second) {
            const std::size_t at = first_difference(source_frame, *found->second);
            throw std::invalid_argument(
                "panel template '" + name + "' fails sequence round-trip at byte " +
                std::to_string(at) + " (spelled " + std::to_string(source_frame.size()) +
                " bytes, authoritative " + std::to_string(found->second->size()) + " bytes)");
        }
    }
}

} // namespace panvar
