#pragma once

#include "panvar/genotype_model.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace panvar {

struct SequenceInterval {
    std::size_t begin = 0;
    std::size_t end = 0;

    std::size_t size() const noexcept { return end - begin; }
    bool operator==(const SequenceInterval& other) const noexcept {
        return begin == other.begin && end == other.end;
    }
    bool operator!=(const SequenceInterval& other) const noexcept { return !(*this == other); }
};

// Always stored in canonical block-chain orientation. A bypass allele has a zero-width block
// interval; it is still an explicit allele state in `mosaic`.
struct HaplotypeSequence {
    MosaicPath mosaic;
    std::string bases;
    std::vector<SequenceInterval> block_intervals;
    std::vector<SequenceInterval> invariant_intervals;
};

struct AuthoritativePanelSequence {
    std::string template_name;
    // Exact bytes spelled by the source graph path in that path's own orientation.
    std::string source_frame_bases;
};

HaplotypeSequence spell_mosaic(
    const PreparedGenotypeModel& model,
    const MosaicPath& mosaic);

HaplotypeSequence spell_panel_template(
    const PreparedGenotypeModel& model,
    std::size_t template_id);

std::string sequence_in_source_frame(
    const std::string& chain_bases,
    ChainOrientation source_orientation);

SequenceInterval interval_in_source_frame(
    SequenceInterval chain_interval,
    std::size_t sequence_size,
    ChainOrientation source_orientation);

// Refuses missing, duplicate, extra, partial or byte-disagreeing panel paths. Authoritative
// sequences stay outside PreparedGenotypeModel so truth/evaluation data cannot enter inference.
void verify_complete_panel_round_trip(
    const PreparedGenotypeModel& model,
    const std::vector<AuthoritativePanelSequence>& authoritative_sequences);

} // namespace panvar
