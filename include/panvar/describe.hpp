#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace panvar {

enum class DescribeFeatureMode {
    AllKmers,
    Minimizer,
    Syncmer
};

struct DescribeOptions {
    std::string gfa_path;
    std::string bubbles_csv_in;
    std::string out_dir = "describe_out";
    std::vector<std::size_t> bubble_ids;
    std::size_t kmer_size = 31;
    // Syncmer sampling is the default: a compact, evenly distributed,
    // substitution-robust marker set. Use AllKmers for the exhaustive set.
    DescribeFeatureMode feature_mode = DescribeFeatureMode::Syncmer;
    std::size_t minimizer_window = 15;
    // 0 selects an automatic closed-syncmer s-mer size from k.
    std::size_t syncmer_s = 0;
    // Symmetric minor-presence (MAF-style) filter: drop a feature when
    // min(present_paths, path_count - present_paths) <= min_feature_paths,
    // unless it varies in copy number across paths. 0 keeps all discriminative
    // features (legacy behavior); 1 (default) drops singletons and all-but-one.
    std::size_t min_feature_paths = 1;
    // 0 disables the safety cap. The sparse JSONL/features outputs are always written.
    std::size_t max_wide_features = 250000;
    bool force_wide_matrix = false;
    bool write_wide_matrix = true;
    bool quiet = false;
};

struct DescribeSummary {
    std::size_t bubbles_processed = 0;
    std::size_t bubbles_with_paths = 0;
    std::size_t paths_written = 0;
    std::size_t features_written = 0;       // discriminative k-mers kept
    std::size_t features_candidates = 0;    // k-mer candidates before filter
    std::size_t matrix_files_written = 0;
    std::size_t jsonl_files_written = 0;
    std::size_t node_edge_features_written = 0;    // node+edge kept
    std::size_t node_edge_candidates = 0;          // node+edge candidates before filter
    std::size_t graph_matrix_files_written = 0;
    std::size_t files_written = 0;
};

void describe_kmers_from_graph(
    const DescribeOptions& options,
    DescribeSummary* summary_out = nullptr);

} // namespace panvar
