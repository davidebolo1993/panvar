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
    DescribeFeatureMode feature_mode = DescribeFeatureMode::AllKmers;
    std::size_t minimizer_window = 15;
    // 0 selects an automatic closed-syncmer s-mer size from k.
    std::size_t syncmer_s = 0;
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
    std::size_t features_written = 0;
    std::size_t matrix_files_written = 0;
    std::size_t jsonl_files_written = 0;
    std::size_t files_written = 0;
};

void describe_kmers_from_graph(
    const DescribeOptions& options,
    DescribeSummary* summary_out = nullptr);

} // namespace panvar
