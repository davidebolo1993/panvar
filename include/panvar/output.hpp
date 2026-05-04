#pragma once

#include <string>
#include <unordered_map>
#include <vector>

#include "panvar/bubbles.hpp"

namespace panvar {

void write_bubbles_csv(
    const std::string& output_path,
    const std::vector<Bubble>& bubbles);

std::vector<Bubble> read_bubbles_csv(const std::string& input_path);

void write_bandage_node_colors_csv(
    const std::string& output_path,
    const std::vector<Bubble>& bubbles);

void write_snarl_debug_tsv(
    const std::string& output_path,
    const std::vector<SnarlDebugEntry>& entries);

std::unordered_map<std::size_t, std::string> level_palette();

} // namespace panvar
