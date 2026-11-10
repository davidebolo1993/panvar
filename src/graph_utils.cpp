#include "panvar/graph_utils.hpp"
#include "panvar/ref_path.hpp"

#include <cctype>

#include <algorithm>
#include <array>
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

// Every IUPAC code has a complement (R<->Y, K<->M, B<->V, D<->H; S, W and N are their own), and case is
// kept, so a reverse-oriented segment spells what the path holds, ambiguity codes and soft-masking
// included. A character that is not a base (e.g. '-') becomes N.
std::string reverse_complement(const std::string& sequence) {
    static const auto table = [] {
        std::array<char, 256> t{};
        t.fill('N');
        const char* from = "ACGTRYKMBVDHSWN";
        const char* to   = "TGCAYRMKVBHDSWN";
        for (int i = 0; from[i] != '\0'; ++i) {
            t[static_cast<unsigned char>(from[i])] = to[i];
            t[static_cast<unsigned char>(std::tolower(static_cast<unsigned char>(from[i])))] =
                static_cast<char>(std::tolower(static_cast<unsigned char>(to[i])));
        }
        t[static_cast<unsigned char>('U')] = 'A';
        t[static_cast<unsigned char>('u')] = 'a';
        return t;
    }();
    std::string rc;
    rc.reserve(sequence.size());
    for (auto it = sequence.rbegin(); it != sequence.rend(); ++it) {
        rc.push_back(table[static_cast<unsigned char>(*it)]);
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

std::string resolve_reference_path_name(const Graph& graph, const std::string& query,
                                        const std::string& module) {
    if (query.empty()) return query;
    std::vector<std::string> names;
    for (const PathRecord& p : graph.paths) names.push_back(p.name);
    std::vector<std::string> hits;
    for (const std::size_t i : match_reference_name(names, query)) hits.push_back(names[i]);
    if (hits.size() == 1) return hits[0];
    if (hits.empty()) {
        // A W line named by its own (0-based) coordinates, which panvar names from start + 1.
        const std::string w = w_literal_as_path_name(query);
        for (const PathRecord& p : graph.paths)
            if (!w.empty() && p.type == 'W' && p.name == w) return p.name;
        throw std::runtime_error(module + ": reference path not found: " + query);
    }
    std::string msg = module + ": reference path ambiguous: " + query + " matches:";
    for (const std::string& h : hits) msg += "\n  " + h;
    throw std::runtime_error(msg);
}

void require_graph_paths(std::size_t segments, std::size_t paths, const std::string& module) {
    if (paths > 0) return;
    if (segments == 0)
        throw std::runtime_error(module + ": input GFA has no segments and no paths (an empty file, or not GFA)");
    throw std::runtime_error(module + ": input GFA has no P/W paths; panvar reads the haplotypes' walks");
}

void validate_graph_paths(const Graph& graph, const std::string& module,
                          bool require_sequences, bool require_zero_overlaps) {
    require_graph_paths(graph.nodes.size(), graph.paths.size(), module);
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
