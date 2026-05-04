#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "panvar/bubbles.hpp"
#include "panvar/gfa.hpp"

namespace panvar {

enum class ClusterMode {
    Sequence,
    SequenceFast,
    Walk
};

struct AlleleCallOptions {
    // Required module-1 handoff: precomputed bubbles CSV from `panvar bubble`.
    std::string bubbles_csv_in;
    // Similarity threshold in [0, 1].
    double min_similarity = 0.90;
    ClusterMode cluster_mode = ClusterMode::Sequence;
    // Worker threads for pairwise distance computations (0 = auto).
    std::size_t threads = 0;
    // Emit progress logs on stderr while processing bubbles.
    bool show_progress = true;
    // Use fast sketch-guided + bounded distance in auto mode.
    bool fast_distance = true;
    // Legacy knob kept for compatibility (not used in current auto path).
    std::size_t exact_distance_max_bp = 4000;
    // Auto-switch from UPGMA to threshold-graph clustering when unique alleles
    // per bubble exceed this value (0 disables auto-switch).
    std::size_t max_upgma_alleles = 256;
    // Optional representative sequence export for each bubble-cluster.
    bool write_cluster_sequences = false;
    std::string cluster_sequences_csv_path;

    // Optional ODGI viz export (per-bubble ordered path list + cluster colors).
    bool write_odgi_viz_inputs = false;
    std::string odgi_viz_out_dir;
    std::string reference_path;
    std::size_t flank_nodes = 1;
    std::size_t max_paths_per_bubble = 0;
    std::string odgi_bin = "odgi";
    std::string odgi_input_path;
    std::size_t odgi_width = 2400;
    std::size_t odgi_path_height = 8;
    bool run_odgi = false;

    // Optional per-bubble similarity diagnostics (distance summaries + dendrogram/tree).
    bool write_similarity_reports = false;
    std::string similarity_out_dir;

    std::size_t min_sv_bp = 50;
    // Optional merged region-level multi-sample VCF (one record per non-reference cluster event).
    bool write_region_vcf = false;
    std::string region_vcf_path;
    // Optional per-cluster debug artifacts (FASTA, PAF, dotplot).
    bool write_debug_reports = false;
    std::string debug_out_dir;
    // Optional GTF annotations to overlay genes on debug dotplots.
    std::string dotplot_gtf_path;
    // Optional regex/term filters for genes highlighted on dotplots.
    // Repeat --dotplot-gene-match to pass multiple patterns.
    std::vector<std::string> dotplot_gene_matches;
    // Minimap2 preset used for per-cluster reference/allele alignment.
    // Supported: asm5, asm10, asm20 (default).
    std::string minimap_preset = "asm20";
    // Number of top chains to extend/alignment in minimap2 (best_n).
    std::size_t minimap_best_n = 8;
    // Keep secondary/supplementary alignments in PAF/debug artifacts.
    bool minimap_emit_secondary = true;
    // INS length for split-chain events:
    //  - true: keep geometric deviation (SVIM-like inter-segment span delta)
    //  - false: derive from realized ALT span in sequence coordinates.
    bool split_ins_use_geometric_svlen = true;
    // Optional INS subtype/CN annotation pass (NOVEL vs DUP-like classes).
    bool classify_ins = false;
    // Merge radius (bp) for collapsing equivalent events across clusters in
    // final region-level VCF.
    std::size_t vcf_merge_window_bp = 20;
    // Merge mode for cross-cluster event compaction in the final VCF:
    //  - strict: require start/end proximity within --vcf-merge-window-bp
    //  - lenient: allow overlapping non-INS events and a wider INS anchor window
    //             (while still requiring sequence similarity).
    std::string vcf_merge_mode = "strict";
    // In lenient mode, event starts can merge within this wider window.
    std::size_t vcf_merge_lenient_window_bp = 100;
    // In lenient mode (non-INS), minimum reference-interval Jaccard overlap
    // required when strict start/end window checks fail.
    double vcf_merge_lenient_min_ref_jaccard = 0.30;
    // Minimum sequence similarity for cross-cluster event merge in final VCF.
    double vcf_merge_min_seq_similarity = 0.80;
    // Bounded edit-distance cap used during sequence-similarity calculation.
    // This is an upper bound on allowed edit fraction while comparing
    // candidate events for merge.
    double vcf_merge_max_seq_edit_fraction = 0.35;
};

struct AlleleCallSummary {
    std::size_t bubbles_processed = 0;
    std::size_t bubbles_with_assignments = 0;
    std::size_t unique_alleles = 0;
    std::size_t clusters = 0;
    std::size_t assignments = 0;
    std::size_t cluster_sequences_written = 0;
    std::size_t bubbles_with_odgi_exports = 0;
    std::size_t odgi_runs_ok = 0;
    std::size_t odgi_runs_failed = 0;
    std::size_t similarity_reports_written = 0;
    std::size_t dotplots_written = 0;
    std::size_t region_vcf_records = 0;
    std::size_t debug_reports_written = 0;
};

void call_alleles_to_csv(
    const Graph& graph,
    const AlleleCallOptions& options,
    const std::string& clusters_csv_path,
    const std::string& assignments_csv_path,
    AlleleCallSummary* summary_out = nullptr);

} // namespace panvar
