#pragma once

#include <cstddef>
#include <string>

namespace panvar {

// `rebuild`: re-induce a fragmented low-complexity locus graph before bubble decomposition. A degree
// gate decides whether the graph is pathological; if so the locus is rebuilt by progressive graph
// generation, driving minigraph as a library. Ours are the two ends: the haplotype order fed to it
// (k-mer richness, since the seed becomes the backbone) and the emitted GFA, which carries per-haplotype
// P lines and preserves link orientation. Healthy graphs pass through untouched.
struct RebuildOptions {
    std::string gfa_path;            // input GFA (the pggb graph for the locus)
    std::string out_path;            // output GFA
    std::size_t kmer = 21;           // Criteria B: k for the k-mer richness metric (distinct k-mers
                                     // first, total k-mers breaking ties)
    std::size_t hub_degree = 50;     // Criteria A: a node with >= this many distinct neighbours is a hub
    std::size_t min_hubs = 10;       // Criteria A: >= this many hubs => pathological
    std::size_t min_var = 50;        // minigraph -L: minimum variant length to augment into the graph
    std::string tmp_dir;             // parent for the per-haplotype FASTA scratch (empty = beside --out);
                                     // a dedicated subfolder is created under it and removed on exit
    std::size_t threads = 0;         // worker threads (0 = auto)
    bool force = false;              // rebuild even if the gate says healthy (testing / small inputs)
    bool quiet = false;
};

struct RebuildSummary {
    std::size_t raw_nodes = 0;
    std::size_t raw_hubs = 0;        // #nodes with degree >= hub_degree
    std::size_t raw_maxdeg = 0;
    double raw_density = 0.0;        // nodes per kb of the longest haplotype span
    bool pathological = false;
    bool ran = false;                // did we actually rebuild (pathological or forced)?
    std::size_t haplotypes = 0;
    std::string seed;                // richest haplotype (Criteria B)
    std::size_t out_nodes = 0;
    std::size_t out_edges = 0;
    std::size_t out_maxdeg = 0;
    std::size_t out_hubs = 0;        // #nodes with degree >= hub_degree in the rebuilt graph
    std::size_t paths_recovered = 0; // haplotypes whose walk mg_map resolved
    double mean_query_cover = 0.0;   // mean fraction of each haplotype covered by its recovered walk
};

RebuildSummary rebuild_graph(const RebuildOptions& options);

} // namespace panvar
