#pragma once

#include <cstddef>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include "panvar/gfa.hpp"

namespace panvar {

struct Bubble {
    std::size_t id = 0;
    std::string source;
    std::string sink;
    std::vector<std::string> inside;
    std::size_t path_support = 0;
    std::size_t min_inside_bp = 0;
    std::size_t max_inside_bp = 0;
    // Supporting paths with inside span >= --min-variant-bp (drives the size filter).
    std::size_t long_path_support = 0;
    // True when an internal node is seen in both orientations across supporting paths.
    bool inversion_signal = false;
    // True when the snarl interior is cyclic (a haplotype revisits an interior node within
    // the snarl, an interior node has a self-loop, or there is an inversion signal). An
    // acyclic, single-entry/single-exit snarl is a superbubble (ultrabubble); --superbubbles
    // keeps only `!cyclic` sites.
    bool cyclic = false;
};

struct BubbleCallOptions {
    // Optional override: JSON-lines from `vg view -R -j <snarls.pb>`. When empty, snarls are
    // found internally (cactus finder, integrated_snarls.hpp); when set, these top-level
    // pairs are used instead.
    std::string snarls_input_path;
    // Reference path (name or unique substring) used to order the internal snarl finder.
    std::string reference_path;
    // Emit only acyclic superbubbles from the internal finder (default: all snarls).
    bool superbubbles_only = false;
    // Pre-computed top-level (source,sink) pairs. When non-empty these are used directly
    // (e.g. the cactus finder run on the sorted GfaModel), overriding both --snarls-in and
    // the order-based fallback finder.
    std::vector<std::pair<std::string, std::string>> snarl_pairs_override;
    bool collect_snarl_debug = false;
    // Keep bubbles only if at least one supporting path has inside length >= this threshold
    // (unless inversion_signal is detected). Set to 0 to disable size filtering.
    std::size_t min_variant_bp = 50;
    std::size_t min_path_support = 0;
    // Optional merge of nearby bubbles on graph distance (bp) between sink/source boundaries.
    // 0 disables merging.
    std::size_t merge_nearby_bp = 0;
    // Suppress the stderr progress bar (output files are unaffected).
    bool quiet = false;
};

struct SnarlDebugEntry {
    std::size_t candidate_id = 0;
    std::string source;
    std::string sink;
    std::size_t inside_node_count = 0;
    std::size_t n_paths = 0;
    std::size_t min_inside_bp = 0;
    std::size_t long_path_support = 0;
    bool inversion_signal = false;
    bool accepted = false;
};

struct BubbleCallReport {
    std::vector<Bubble> bubbles;
    // Candidate non-SNP bubbles before final filters (used for overview coloring).
    std::vector<Bubble> non_snp_bubbles;
    std::vector<SnarlDebugEntry> snarl_debug;
};

BubbleCallReport call_bubbles_report(const Graph& graph, const BubbleCallOptions& options);
std::vector<Bubble> call_bubbles(const Graph& graph, const BubbleCallOptions& options);

std::unordered_set<std::string> nodes_in_bubble(const Bubble& bubble);

} // namespace panvar
