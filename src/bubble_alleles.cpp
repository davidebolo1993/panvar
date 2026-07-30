#include "panvar/bubble_alleles.hpp"

#include "panvar/graph_utils.hpp"

#include <unordered_map>
#include <utility>

namespace panvar {

BubbleAlleleSet enumerate_bubble_alleles(
    const Graph& graph,
    const std::vector<BubblePathIndex>& path_indexes,
    const Bubble& bubble,
    const std::string& reference_path_name) {

    BubbleAlleleSet out;
    out.bubble_id = bubble.id;

    std::size_t ref_idx = graph.paths.size();
    for (std::size_t pi = 0; pi < graph.paths.size(); ++pi) {
        if (graph.paths[pi].name == reference_path_name) { ref_idx = pi; break; }
    }
    if (ref_idx < graph.paths.size()) {
        auto ref_opt = bubble_steps(graph.paths[ref_idx], path_indexes[ref_idx], bubble);
        if (ref_opt.has_value() && !ref_opt->empty()) {
            out.reference_steps = std::move(*ref_opt);
            out.reference_signature = build_walk_signature(out.reference_steps);
            out.has_reference = true;
        }
    }

    std::unordered_map<std::string, std::size_t> sig_to_allele;
    for (std::size_t pi = 0; pi < graph.paths.size(); ++pi) {
        const auto steps_opt = bubble_steps(graph.paths[pi], path_indexes[pi], bubble);
        if (!steps_opt.has_value() || steps_opt->empty()) continue;
        const std::vector<PathStep>& steps = *steps_opt;
        out.traversing.insert(graph.paths[pi].name);
        std::string sig = build_walk_signature(steps);
        auto it = sig_to_allele.find(sig);
        if (it == sig_to_allele.end()) {
            sig_to_allele.emplace(sig, out.alleles.size());
            BubbleAllele a;
            a.id = out.alleles.size();
            a.signature = std::move(sig);
            a.steps = steps;
            a.members.push_back(graph.paths[pi].name);
            out.alleles.push_back(std::move(a));
        } else {
            out.alleles[it->second].members.push_back(graph.paths[pi].name);
        }
    }

    return out;
}

} // namespace panvar
