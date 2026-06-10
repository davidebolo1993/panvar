#pragma once

#include "panvar/gfa.hpp"

#include <cstddef>
#include <cstdint>
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

// Walk identity helpers shared across modules: a stable per-step hash token, the
// token vector for a walk, and a human-readable "id+/id-" signature string.
std::uint64_t hash_step_token(const PathStep& step);
std::vector<std::uint64_t> build_walk_tokens(const std::vector<PathStep>& steps);
std::string build_walk_signature(const std::vector<PathStep>& steps);

} // namespace panvar
