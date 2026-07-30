#include "panvar/graph_utils.hpp"

#include <algorithm>
#include <cstdint>
#include <stdexcept>

namespace panvar {

std::uint64_t hash_step_token(const PathStep& step) {
    // 64-bit FNV-1a over node id + orientation.
    std::uint64_t h = 1469598103934665603ULL;
    for (const unsigned char c : step.node_id) {
        h ^= static_cast<std::uint64_t>(c);
        h *= 1099511628211ULL;
    }
    h ^= (step.reverse ? 0xF0ULL : 0x0FULL);
    h *= 1099511628211ULL;
    return h;
}

std::vector<std::uint64_t> build_walk_tokens(const std::vector<PathStep>& steps) {
    std::vector<std::uint64_t> tokens;
    tokens.reserve(steps.size());
    for (const auto& step : steps) {
        tokens.push_back(hash_step_token(step));
    }
    return tokens;
}

std::string build_walk_signature(const std::vector<PathStep>& steps) {
    std::string sig;
    sig.reserve(steps.size() * 8);
    for (std::size_t i = 0; i < steps.size(); ++i) {
        if (i > 0) {
            sig.push_back(',');
        }
        sig += steps[i].node_id;
        sig.push_back(steps[i].reverse ? '-' : '+');
    }
    return sig;
}

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

std::unordered_set<std::string> self_loop_nodes(const Graph& graph) {
    std::unordered_set<std::string> out;
    for (const auto& [id, node] : graph.nodes) {
        bool loop = false;
        for (const Neighbor& nb : node.start) if (nb.node_id == id) { loop = true; break; }
        if (!loop) for (const Neighbor& nb : node.end) if (nb.node_id == id) { loop = true; break; }
        if (loop) out.insert(id);
    }
    return out;
}

} // namespace panvar
