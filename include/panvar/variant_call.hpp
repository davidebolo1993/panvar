#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "panvar/gfa.hpp"

namespace panvar {

struct VariantCallOptions {
    std::string gfa_path;                // the graph's own path, so an output cannot alias an input
    std::string bubbles_csv_in;          // module-1 bubbles CSV
    std::string reference_path;          // required: graph path used as the diff baseline
    std::string out_prefix;              // writes <prefix>.bubble_<id>.vcf + <prefix>.region.vcf
    std::size_t min_sv_bp = 50;          // minimum event size for a reported (merged) record
    std::size_t merge_distance_bp = 100; // within-bubble nearby same-type coalescing distance
    double merge_jaccard = 0.80;         // cross-haplotype length-weighted node-set Jaccard threshold
    double merge_seq_identity = 0.80;    // cross-haplotype event-sequence identity threshold (alt merge key)
    double merge_size_ratio = 0.0;       // length-ratio floor for the sequence merge (0 -> use merge_seq_identity)
    std::size_t min_haplotypes = 1;      // drop a record carried by fewer than this many haplotypes
    double min_maf = 0.0;                // drop a record whose written AF (AC/AN) is below this (0 = off)
    // One column per PanSN sample (SAMPLE#HAP#...), its haplotypes joined into a phased GT, the
    // reference path left out -- the layout vg deconstruct writes. Default: one column per path.
    bool diploid = false;
    // Compare each haplotype's allele with the reference's by aligning their sequences, where the bubble
    // carries no copy-number record, instead of node by node (--no-align-alleles).
    bool align_alleles = true;
    std::size_t rescue_min_bp = 0;       // floor for events kept for merge/rescue (0 -> min_sv_bp/2)
    bool classify_ins = false;           // minimap2 INS subtype refinement (NOVEL vs DUP)
    bool cn = false;                     // copy-number calling: count a collapsed module's copies (MODULE).
                                         // The UNIT route (a folded repeat's copies) is always on.
    std::string minimap_preset = "asm20";
    std::size_t minimap_best_n = 8;
    double ins_dup_min_identity = 0.90;  // identity for an INS to be subtyped DUP
    std::string gtf_path;                // optional reference-coordinate GTF: annotate variants with the
                                         // genes they touch (INFO GENES), write <prefix>.node_genes.tsv,
                                         // and a per-gene DUP copy-number table; needs a PanSN reference
    // panphorte's <prefix>.panphorte.unit_graph.tsv, when present: each folded repeat's unit graph is one
    // copy-number unit, counted at its anchor nodes (one pass per copy), as a unit node is counted.
    std::string unit_graph_tsv;
    std::size_t max_seq_bp = 20000;      // INSSEQ, DELSEQ and INVSEQ left out above this many bases (0: none)
    std::vector<std::size_t> bubble_ids; // if non-empty, restrict to these bubbles
    bool write_per_bubble_vcf = true;    // also keep each <prefix>.bubble_<id>.vcf
    bool write_variant_nodes = true;     // write <prefix>.variant_nodes.tsv (per-variant node set)
    std::size_t threads = 0;             // worker threads inside each bubble (0 = auto); output is
                                         // identical regardless of thread count
    bool quiet = false;
};

struct VariantCallSummary {
    std::size_t bubbles_seen = 0;
    std::size_t bubbles_with_reference = 0;
    std::size_t bubbles_with_calls = 0;
    std::size_t records_written = 0;
    // Divergent blocks the node-token DP abandoned (over its cell cap) and no base-level alignment of the
    // allele then covered: changes there are not called, so the count is reported.
    std::size_t skipped_large_segments = 0;
    // Module copy changes (--cn), per haplotype: carriers of a module record, haplotypes whose copy change is a
    // larger event of their own, carriers none of whose changes holds a whole anchor pass (no form), and
    // carriers whose forms gain or lose other anchor passes than their copy number says.
    std::size_t copy_carriers = 0, copy_own_event = 0, copy_no_form = 0, copy_pass_off = 0;
    std::size_t del = 0;
    std::size_t ins = 0;
    std::size_t inv = 0;
    std::size_t dup = 0;
    std::size_t sub = 0;   // SVTYPE=SUB: replacements whose net change is under --min-sv-bp
    // Module copy numbers not called because the graph is too tangled there to read copies from it:
    // a site holding many hub nodes, or a module spanning most of the reference.
    std::vector<std::size_t> tangle_ids;      // their bubble ids, a site with many hubs
    std::vector<std::size_t> oversized_ids;   // a module spanning most of the reference
    // Bubbles with at least one traversing path that traverses them more than once (two or more runs
    // from one boundary to the other through interior nodes only, anywhere on the path). Its walk
    // either spans the sequence between them or types only one; every record there carries
    // MULTI_CROSSING.
    std::size_t multi_crossing_bubbles = 0;
    std::size_t multi_crossing_paths = 0;      // those paths, summed over the bubbles
    // --diploid: records not written because no sample column carries an alternate (every carrier's slot is
    // a haplotype whose contigs disagree at the site).
    std::size_t diploid_records_dropped = 0;
    // Of the multi-crossing paths, those with a traversal outside the scored walk that differs from every
    // reference traversal by --min-sv-bp or more in length: a change no record describes
    // (MULTI_CROSSING_UNSCORED).
    std::size_t multi_crossing_unscored_bubbles = 0;
    std::size_t multi_crossing_unscored_paths = 0;
    std::vector<std::size_t> multi_crossing_unscored_ids;   // those bubbles
    // --gtf found genes in the region and <prefix>.node_genes.tsv / .dup_gene_cn.tsv were written.
    bool gene_tables_written = false;
    std::size_t per_bubble_vcfs = 0;   // <prefix>.bubble_<id>.vcf files written, one per site with records
    // Haplotype-site pairs set missing because the walk crosses an assembly gap (a run of >= 10 N).
    std::size_t gap_missing = 0;
    std::size_t ref_gap_bubbles = 0;   // sites skipped: the reference's own allele holds an assembly gap
    std::vector<std::size_t> ref_gap_ids;   // their bubble ids
    // Bubbles where a haplotype walk held a stretch over the alignment cap, and the haplotype-records
    // there written GT '.' instead of 0 because that haplotype could not be compared in full.
    std::vector<std::size_t> unaligned_bubbles;
    std::size_t unaligned_missing = 0;
    // Haplotype-records written GT '.' because the haplotype's own larger deletion takes in the record's bases.
    std::size_t covered_missing = 0;
};

void call_variants(
    const Graph& graph,
    const VariantCallOptions& options,
    VariantCallSummary* summary_out = nullptr);

} // namespace panvar
