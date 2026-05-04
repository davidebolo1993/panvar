#pragma once

#include <cstddef>
#include <string>
#include <unordered_set>
#include <vector>

#include "panvar/gfa.hpp"

namespace panvar {

enum class BubbleType {
    Simple,
    Super,
    Insertion
};

enum class SiteMode {
    Superbubble,
    Snarl
};

struct Bubble {
    std::size_t id = 0;
    std::string source;
    std::string sink;
    std::vector<std::string> inside;
    SiteMode site_mode = SiteMode::Snarl;
    BubbleType type = BubbleType::Super;
    std::size_t parent_id = 0;
    std::size_t nesting_level = 1;
    std::size_t path_support = 0;
    std::size_t min_inside_bp = 0;
    std::size_t max_inside_bp = 0;
    std::size_t long_path_support = 0;
    bool inversion_signal = false;
};

struct BubbleCallOptions {
    // Required: JSON-lines from `vg view -R -j <snarls.pb>`.
    std::string snarls_input_path;
    bool collect_snarl_debug = false;
    // Keep bubbles only if at least one supporting path has inside length >= this threshold
    // (unless inversion_signal is detected). Set to 0 to disable size filtering.
    std::size_t min_variant_bp = 50;
    std::size_t min_path_support = 0;
    std::size_t max_nesting_level = 0;
};

struct SnarlDebugEntry {
    std::size_t candidate_id = 0;
    std::string mode;
    std::string source;
    std::string sink;
    std::size_t component_handle_count = 0;
    std::size_t component_node_count = 0;
    std::size_t boundary_node_count = 0;
    std::size_t inside_node_count = 0;
    std::size_t n_paths = 0;
    std::size_t min_inside_bp = 0;
    bool nested = false;
    bool directed_acyclic_net_graph = false;
    bool accepted = false;
    std::string reason;
    std::string boundary_nodes;
};

struct BubbleCallReport {
    std::vector<Bubble> bubbles;
    std::vector<SnarlDebugEntry> snarl_debug;
};

BubbleCallReport call_bubbles_report(const Graph& graph, const BubbleCallOptions& options);
std::vector<Bubble> call_bubbles(const Graph& graph, const BubbleCallOptions& options);

std::string bubble_type_to_string(BubbleType type);
std::string site_mode_to_string(SiteMode mode);

std::unordered_set<std::string> nodes_in_bubble(const Bubble& bubble);

} // namespace panvar
