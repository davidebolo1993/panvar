#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace panvar {

enum class DescribeFeatureMode {
    AllKmers,
    Syncmer
};

struct DescribeOptions {
    std::string gfa_path;
    std::string bubbles_csv_in;
    std::string out_dir = "describe_out";
    std::vector<std::size_t> bubble_ids;
    std::size_t kmer_size = 31;
    DescribeFeatureMode feature_mode = DescribeFeatureMode::Syncmer;
    std::size_t syncmer_s = 0;   // 0 = from k
    // Drop a feature with min(present, absent) <= min_feature_paths, unless its count varies.
    std::size_t min_feature_paths = 1;
    std::size_t max_wide_features = 250000;   // 0 = no cap on a site's dense matrix
    bool force_wide_matrix = false;
    bool write_wide_matrix = true;
    // The pooled BIMBAM matrices, <out_dir>/<level>/<substrate>/, with their sidecars and column order.
    bool bimbam = true;
    // Substrates; --only-* select. The variant one also needs variant_vcf_path.
    bool emit_kmers = true;
    bool emit_graph = true;
    bool emit_variant = true;
    // call's <prefix>.variant_nodes.tsv: only its sites, and of them only the variant nodes and the
    // boundaries, give k-mer and graph features.
    std::string variant_nodes_path;
    // With variant_nodes_path: also the bases this far from a variant node (k-mers), and the nodes they
    // reach (graph features).
    std::size_t variant_flank_bp = 0;
    // cosigt sample table: also the sample level, each sample the sum of its haplotypes.
    std::string samples_path;
    // call's region VCF: the variant substrate, one row per record.
    std::string variant_vcf_path;
    // The reference path: k-mer and graph rows get its coordinates (chrom, pos, end) in their sidecars.
    std::string reference_path;
    // Rescale each BIMBAM feature to 0..2 (per-feature min-max), for tools expecting a diploid dosage.
    bool scale_dosage = false;
    std::size_t threads = 0;   // 0 = hardware concurrency; output does not depend on it
    bool quiet = false;
};

struct DescribeSummary {
    std::size_t bubbles_processed = 0;
    std::size_t bubbles_with_paths = 0;
    std::size_t paths_written = 0;
    std::size_t features_written = 0;       // k-mers kept
    std::size_t features_candidates = 0;    // k-mer candidates before filter
    std::size_t matrix_files_written = 0;
    std::size_t jsonl_files_written = 0;
    std::size_t node_edge_features_written = 0;    // node+edge kept
    std::size_t node_edge_candidates = 0;          // node+edge candidates before filter
    std::size_t graph_matrix_files_written = 0;
    std::size_t files_written = 0;
    std::size_t haplotypes = 0, nodes = 0;         // the graph read
    std::size_t variant_rows = 0;                  // records in the variant matrix
    std::size_t samples = 0;                       // with --samples
};

void describe_kmers_from_graph(
    const DescribeOptions& options,
    DescribeSummary* summary_out = nullptr);

// The variant substrate, from call's region VCF alone (no GFA needed).
void describe_variant_from_vcf(
    const DescribeOptions& options,
    DescribeSummary* summary_out = nullptr);

} // namespace panvar
