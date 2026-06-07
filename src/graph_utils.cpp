#include "panvar/graph_utils.hpp"

#include <algorithm>
#include <stdexcept>

namespace panvar {

std::unordered_map<std::string, const PathRecord*> path_records_by_name(const Graph& graph) {
    std::unordered_map<std::string, const PathRecord*> out;
    out.reserve(graph.paths.size() * 2);
    for (const auto& path : graph.paths) {
        out[path.name] = &path;
    }
    return out;
}

std::vector<std::size_t> path_prefix_bp(
    const PathRecord& path,
    const std::unordered_map<std::string, Node>& nodes) {

    std::vector<std::size_t> pref(path.steps.size() + 1, 0);
    for (std::size_t i = 0; i < path.steps.size(); ++i) {
        const auto it = nodes.find(path.steps[i].node_id);
        const std::size_t len = (it == nodes.end()) ? 1 : std::max<std::size_t>(1, it->second.sequence.size());
        pref[i + 1] = pref[i] + len;
    }
    return pref;
}

std::string reverse_complement(const std::string& sequence) {
    std::string rc;
    rc.reserve(sequence.size());
    for (auto it = sequence.rbegin(); it != sequence.rend(); ++it) {
        switch (*it) {
            case 'A':
            case 'a':
                rc.push_back('T');
                break;
            case 'C':
            case 'c':
                rc.push_back('G');
                break;
            case 'G':
            case 'g':
                rc.push_back('C');
                break;
            case 'T':
            case 't':
                rc.push_back('A');
                break;
            default:
                rc.push_back('N');
                break;
        }
    }
    return rc;
}

std::string spell_path_steps_sequence(
    const Graph& graph,
    const std::vector<PathStep>& steps,
    bool* complete) {

    if (complete != nullptr) {
        *complete = true;
    }

    std::size_t total_len = 0;
    for (const auto& step : steps) {
        const auto node_it = graph.nodes.find(step.node_id);
        if (node_it == graph.nodes.end() || node_it->second.sequence.empty()) {
            if (complete != nullptr) {
                *complete = false;
                return "";
            }
            throw std::runtime_error("Cannot spell sequence: missing sequence for node " + step.node_id);
        }
        total_len += node_it->second.sequence.size();
    }

    std::string seq;
    seq.reserve(total_len);
    for (const auto& step : steps) {
        const auto& node_seq = graph.nodes.at(step.node_id).sequence;
        seq += step.reverse ? reverse_complement(node_seq) : node_seq;
    }
    return seq;
}

} // namespace panvar
