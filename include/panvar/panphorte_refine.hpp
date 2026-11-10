#pragma once

// panphorte --refine: at each folded site, realign the flanks around every carrier's array block -- one
// alignment for the flanks before it, one for the flanks after -- and keep the result only when it is
// simpler. The array block (first to last unit step) is kept verbatim; haplotypes that walk no unit node
// are left as folded. Internal to panphorte; separate from the standalone `refine` module.

#include <cstddef>
#include <string>
#include <unordered_set>
#include <vector>

#include "panvar/gfa_io.hpp"

namespace panvar {

struct PanphorteRefineSite {
    std::size_t bubble_id = 0;
    std::string source;
    std::string sink;
    // The site's interior in the folded graph: the input bubble's surviving interior plus every node
    // folding spliced in for it (unit nodes and fragments).
    std::unordered_set<std::string> interior;
};

struct PanphorteRefineLimits {
    std::size_t max_poa_bp = 5000;
    std::size_t max_walks = 500;
    std::size_t max_poa_work = 0;   // 0: no limit
};

struct PanphorteRefineRow {
    std::size_t bubble_id = 0;
    std::string source;
    std::string sink;
    std::string status;   // refined | kept_original:not_simpler | skipped:<reason>
    std::string detail;
    std::size_t paths = 0;           // traversing paths
    std::size_t carriers = 0;        // traversing paths that walk a unit node
    std::size_t distinct_pre = 0;    // distinct non-empty flanks before the array block
    std::size_t distinct_post = 0;   // and after it
    // nodes / steps / single-path nodes: non-unit, outside the array blocks, over every traversing path
    std::size_t nodes_before = 0, nodes_after = 0;
    std::size_t steps_before = 0, steps_after = 0;
    std::size_t single_path_nodes_before = 0, single_path_nodes_after = 0;
    bool measured = false;           // before/after figures exist
};

struct PanphorteRefineResult {
    std::vector<PanphorteRefineRow> rows;
    std::size_t refined = 0;
    std::size_t kept = 0;
    std::size_t skipped = 0;
    std::size_t nodes_added = 0;
    std::size_t nodes_removed = 0;
    std::size_t edges_added = 0;
    std::size_t edges_removed = 0;
};

// Rewrites `model` in place. Every acceptance check (lossless against the folded graph, unit
// invariance, links, identical carriers on identical routes, no surviving old route) runs before
// returning; a failure throws and the caller writes nothing.
// `rep_nodes` are the unit nodes folding created. A unit is not recognised by a self-loop: one no haplotype
// walks twice in a row has none, and would be realigned away as flank.
PanphorteRefineResult panphorte_refine_sites(GfaModel& model, const std::vector<PanphorteRefineSite>& sites,
                                             const std::string& resolved_reference,
                                             const PanphorteRefineLimits& limits,
                                             const std::unordered_set<std::string>& rep_nodes);

} // namespace panvar
