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

// Inversion sites (Bubble::inversion_parent) enclose other sites; they are read only when asked for, so
// a consumer that needs disjoint sites -- folding, refining -- sees the ordinary ones alone.
std::vector<Bubble> read_bubbles_csv(const std::string& input_path, bool with_inversion_sites = false);

void write_bandage_node_colors_csv(
    const std::string& output_path,
    const std::vector<Bubble>& retained_bubbles,
    const std::vector<Bubble>& non_snp_bubbles = {});

void write_snarl_debug_tsv(
    const std::string& output_path,
    const std::vector<SnarlDebugEntry>& entries);

} // namespace panvar
