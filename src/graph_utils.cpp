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

void validate_graph_paths(const Graph& graph, const std::string& module,
                          bool require_sequences, bool require_zero_overlaps) {
    if (graph.paths.empty()) throw std::runtime_error(module + ": input GFA has no P/W paths");
    std::unordered_set<std::string> seen;
    for (const PathRecord& p : graph.paths) {
        if (p.name.empty()) throw std::runtime_error(module + ": input GFA has a path with no name");
        if (!seen.insert(p.name).second)
            throw std::runtime_error(module + ": input GFA has a duplicate path name: " + p.name +
                                     " (which of the two a result refers to would be undefined)");
        if (p.steps.empty()) throw std::runtime_error(module + ": path has no steps: " + p.name);
        for (const PathStep& step : p.steps) {
            const auto it = graph.nodes.find(step.node_id);
            if (it == graph.nodes.end())
                throw std::runtime_error(module + ": path " + p.name +
                                         " references a node that is not in the graph: " + step.node_id);
            if (require_sequences && (it->second.sequence.empty() || it->second.sequence == "*"))
                throw std::runtime_error(module + ": node " + step.node_id + " on path " + p.name +
                                         " has no sequence (S line is '*')");
        }
        // Every consecutive pair of steps must be joined by an ORIENTED link that exists. A path
        // describing a traversal the graph does not permit spells a sequence no walk could produce,
        // and every later comparison would be made against that.
        for (std::size_t i = 1; i < p.steps.size(); ++i) {
            const PathStep& a = p.steps[i - 1];
            const PathStep& b = p.steps[i];
            const auto ita = graph.nodes.find(a.node_id);
            if (ita == graph.nodes.end()) continue;               // already reported above
            // Leave `a` by its end when forward, by its start when reverse; enter `b` at its start
            // when forward, at its end when reverse.
            const std::vector<Neighbor>& side = a.reverse ? ita->second.start : ita->second.end;
            const int want = b.reverse ? 1 : 0;
            bool linked = false;
            for (const Neighbor& n : side)
                if (n.node_id == b.node_id && n.side == want) { linked = true; break; }
            if (!linked)
                throw std::runtime_error(module + ": path " + p.name + " steps from " + a.node_id +
                                         (a.reverse ? "-" : "+") + " to " + b.node_id +
                                         (b.reverse ? "-" : "+") + " but the graph has no such link");
        }
    }
    if (!require_zero_overlaps) return;
    // Spelling and span measurement concatenate whole segments, so a non-zero overlap double-counts
    // the overlapping bases in every length and identity figure. '*' means UNKNOWN, not zero, so it
    // cannot be assumed away either.
    for (const auto& kv : graph.nodes) {
        for (const std::vector<Neighbor>* side : {&kv.second.start, &kv.second.end})
            for (const Neighbor& n : *side)
                if (n.overlap != 0)
                    throw std::runtime_error(
                        module + ": link " + kv.first + " -- " + n.node_id +
                        (n.overlap < 0 ? std::string(" has an UNKNOWN overlap ('*')")
                                       : " has a non-zero overlap (" + std::to_string(n.overlap) + "M)") +
                        "; this module measures by concatenation, which is only correct when every "
                        "overlap is a verified 0M");
    }
}

} // namespace panvar
