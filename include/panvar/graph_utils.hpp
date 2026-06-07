#pragma once

#include "panvar/gfa.hpp"

#include <cstddef>
#include <string>
#include <unordered_map>
#include <vector>

namespace panvar {

std::unordered_map<std::string, const PathRecord*> path_records_by_name(const Graph& graph);
std::vector<std::size_t> path_prefix_bp(
    const PathRecord& path,
    const std::unordered_map<std::string, Node>& nodes);
std::string reverse_complement(const std::string& sequence);
std::string spell_path_steps_sequence(
    const Graph& graph,
    const std::vector<PathStep>& steps,
    bool* complete = nullptr);

} // namespace panvar
