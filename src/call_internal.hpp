#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

#include "panvar/allele.hpp"
#include "panvar/gfa.hpp"

namespace panvar {

struct PathAssignment {
    std::size_t unique_idx = 0;
    std::string path_name;
    std::size_t interval_start = 0;
    std::size_t interval_end = 0;
    bool source_to_sink = true;
};

struct UniqueAllele {
    std::size_t allele_id = 0;
    std::vector<PathStep> steps;
    std::string signature;
    std::string sequence;
    std::string compare_token;
    std::vector<std::uint64_t> compare_steps;
    std::vector<int> compare_step_weights;
    std::size_t compare_steps_weight_sum = 0;
    bool uses_sequence_similarity = false;
    std::size_t path_support = 0;
    std::size_t sequence_length = 0;
};

struct Cluster {
    std::size_t cluster_id = 0;
    std::size_t representative_idx = 0;
    std::vector<std::size_t> member_unique_idxs;
    std::size_t total_path_support = 0;
};

struct PrecomputedClusterRow {
    std::size_t bubble_id = 0;
    std::string source;
    std::string sink;
    std::size_t cluster_id = 0;
    std::size_t representative_allele_id = 0;
    std::size_t total_path_support = 0;
    std::vector<std::size_t> member_alleles;
};

struct PrecomputedAssignmentRow {
    std::size_t bubble_id = 0;
    std::string source;
    std::string sink;
    std::string path_name;
    std::size_t cluster_id = 0;
    std::size_t allele_id = 0;
    std::size_t allele_length = 0;
    std::size_t interval_start = 0;
    std::size_t interval_end = 0;
    bool source_to_sink = true;
};

using PrecomputedClustersByBubble = std::unordered_map<std::size_t, std::vector<PrecomputedClusterRow>>;
using PrecomputedAssignmentsByBubble = std::unordered_map<std::size_t, std::vector<PrecomputedAssignmentRow>>;

struct VariantBubbleReport {
    struct ClusterDebugStatus {
        std::size_t cluster_id = 0;
        bool is_reference_cluster = false;
        std::size_t representative_allele_id = 0;
        std::string representative_haplotype = ".";
        std::size_t reference_length_bp = 0;
        std::size_t cluster_length_bp = 0;
        std::size_t paf_records = 0;
        bool minimap_best_ok = false;
        std::string orientation = ".";
        double best_edit_distance_norm = 0.0;
        std::size_t event_count = 0;
        bool dotplot_written = false;
        bool pairwise_vcf_written = false;
        std::string status = "pending";
    };

    std::size_t reference_cluster_id = 0;
    std::size_t dotplots_written = 0;
    std::size_t debug_reports_written = 0;
    std::string status = "ok";
    bool has_reference_assignment = false;
    std::size_t reference_interval_start = 0;
    std::size_t reference_interval_end = 0;
    std::size_t clusters_total = 0;
    std::size_t clusters_with_minimap_hit = 0;
    std::size_t clusters_with_events = 0;
    std::unordered_map<std::string, std::size_t> cluster_by_path;
    std::unordered_map<std::string, std::size_t> allele_length_by_path;
    std::vector<ClusterDebugStatus> cluster_statuses;
    struct VariantRecordForVcf {
        std::size_t cluster_id = 0;
        std::size_t event_index = 0;
        std::size_t representative_allele_id = 0;
        std::string representative_haplotype = ".";
        std::size_t cluster_path_support = 0;
        std::string event_type;
        std::string event_subtype;
        std::string orientation;
        std::size_t reference_length = 0;
        std::size_t cluster_length = 0;
        long long length_delta = 0;
        double best_edit_distance_norm = 0.0;
        std::size_t inserted_bp = 0;
        std::string member_alleles;
        std::size_t ref_offset_start_bp = 0;
        std::size_t ref_offset_end_bp = 0;
        std::size_t alt_offset_start_bp = 0;
        std::size_t alt_offset_end_bp = 0;
        std::size_t event_bp = 0;
        std::string event_sequence;
        int cn_delta = 0;
        bool preserve_ins_svlen = false;
        bool has_dup_evidence = false;
        double dup_best_similarity = 0.0;
        std::size_t dup_ref_start_bp = 0;
        std::size_t dup_ref_end_bp = 0;
        char dup_orientation = '+';
        std::size_t dup_unit_bp = 0;
        std::size_t dup_ref_copy_number = 0;
        std::size_t dup_added_copies = 0;
        std::size_t dup_alt_copy_number = 0;
        double dup_copy_ratio = 0.0;
        std::string pangene_cn_delta;
        std::string pangene_gain_genes;
        std::string pangene_loss_genes;
        std::size_t pangene_gain_copies = 0;
        std::size_t pangene_loss_copies = 0;
    };
    std::vector<VariantRecordForVcf> vcf_rows;
};

struct ParsedReferencePath {
    std::string chrom;
    std::size_t region_start_1based = 1;
    bool has_interval = false;
};

struct RegionVcfBubble {
    std::size_t bubble_id = 0;
    std::string source;
    std::string sink;
    std::size_t reference_cluster_id = 0;
    std::size_t reference_interval_start = 0;
    std::size_t reference_interval_end = 0;
    std::unordered_map<std::string, std::size_t> cluster_by_path;
    std::unordered_map<std::string, std::size_t> allele_length_by_path;
    std::vector<VariantBubbleReport::VariantRecordForVcf> variant_rows;
};

struct GeneAnnotation {
    std::size_t start0 = 0;
    std::size_t end0 = 0;
    std::string name;
};

using GeneAnnotationIndex = std::unordered_map<std::string, std::vector<GeneAnnotation>>;

ParsedReferencePath parse_reference_path_label(const std::string& path_name);
GeneAnnotationIndex load_gene_annotations_from_gtf(const std::string& gtf_path);
VariantBubbleReport write_variant_reports_for_bubble(
    const Bubble& bubble,
    const std::vector<UniqueAllele>& unique_alleles,
    const std::vector<Cluster>& clusters,
    const std::vector<std::size_t>& cluster_of_unique,
    const std::vector<PathAssignment>& assignments,
    const std::string& reference_path,
    const std::string& reference_sequence_global,
    std::size_t min_sv_bp,
    bool write_debug_files,
    const std::string& debug_output_dir,
    const GeneAnnotationIndex* dotplot_gene_index,
    const std::vector<std::string>& dotplot_gene_matches,
    const std::string& minimap_preset,
    std::size_t minimap_best_n,
    bool minimap_emit_secondary,
    bool split_ins_use_geometric_svlen,
    bool classify_ins);
void write_region_level_vcf(
    const std::string& output_path,
    const std::string& reference_path,
    const ParsedReferencePath& reference_meta,
    const std::string& reference_sequence,
    const std::vector<std::size_t>& reference_prefix_bp,
    const std::vector<RegionVcfBubble>& bubbles,
    const std::vector<std::string>& sample_names,
    std::size_t region_start_1based,
    std::size_t merge_window_bp,
    const std::string& merge_mode,
    std::size_t lenient_window_bp,
    double lenient_min_ref_jaccard,
    double min_seq_similarity,
    double max_seq_edit_fraction,
    std::size_t* records_written_out);

void validate_call_merge_options(const AlleleCallOptions& options);

void call_variants_from_precomputed_grouped_impl(
    const Graph& graph,
    const AlleleCallOptions& options,
    const std::vector<Bubble>& bubbles,
    const PrecomputedClustersByBubble& clusters_by_bubble,
    const PrecomputedAssignmentsByBubble& assignments_by_bubble,
    AlleleCallSummary* summary_out);

} // namespace panvar
