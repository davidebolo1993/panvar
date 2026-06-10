#pragma once

#include <string>
#include <vector>

#include "panvar/bubbles.hpp"

namespace panvar {

// Split a single CSV line into fields, honoring double-quoted fields and
// escaped (doubled) quotes. Shared by the CSV readers across modules.
std::vector<std::string> split_csv_line(const std::string& line);

void write_bubbles_csv(
    const std::string& output_path,
    const std::vector<Bubble>& bubbles);

std::vector<Bubble> read_bubbles_csv(const std::string& input_path);

void write_bandage_node_colors_csv(
    const std::string& output_path,
    const std::vector<Bubble>& retained_bubbles,
    const std::vector<Bubble>& non_snp_bubbles = {});

void write_snarl_debug_tsv(
    const std::string& output_path,
    const std::vector<SnarlDebugEntry>& entries);

} // namespace panvar
