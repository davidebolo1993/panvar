#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "panvar/gfa.hpp"

namespace panvar {

struct VariantCallOptions {
    std::string bubbles_csv_in;          // module-1 bubbles CSV
    std::string reference_path;          // required: graph path used as the diff baseline
    std::string out_prefix;              // writes <prefix>.bubble_<id>.vcf + <prefix>.region.vcf
    std::size_t min_sv_bp = 50;          // minimum event size for a reported (merged) record
    std::size_t merge_distance_bp = 100; // within-bubble nearby same-type coalescing distance
    double merge_jaccard = 0.80;         // cross-haplotype length-weighted node-set Jaccard threshold
    double merge_seq_identity = 0.80;    // cross-haplotype event-sequence identity threshold (alt merge key)
    double merge_size_ratio = 0.0;       // length-ratio floor for the sequence merge (0 -> use merge_seq_identity)
    std::size_t min_haplotypes = 1;      // drop a record carried by fewer than this many haplotypes
    double min_maf = 0.0;                 // drop a record with AF (carriers/traversing-haps) below this (0=off)
    bool multiallelic_loci = false;       // opt-in: collapse a bounded locus into ONE multiallelic record
                                         // (REF + ALT1,ALT2,... explicit seqs; GT indexes the allele)
    std::size_t multiallelic_max_bp = 5000; // skip multiallelic collapse if any allele seq exceeds this
    // Lossless companion VCF: one record per bubble carrying every distinct allele as explicit
    // sequence, each haplotype's GT indexing its own allele. Written alongside the region VCF.
    bool allele_vcf = false;
    std::size_t allele_vcf_max_bp = 0;    // skip a bubble whose alleles exceed this (0 = no limit)
    std::size_t rescue_min_bp = 0;       // floor for events kept for merge/rescue (0 -> min_sv_bp/2)
    bool classify_ins = false;           // minimap2 INS subtype refinement (NOVEL vs DUP)
    bool cn = false;                     // copy-number calling: enable all CN routes, resolved by topology
                                         // (self-loop REP > coverage bp/unit > peak multiplicity). The
                                         // always-on self-loop REP DUP fires regardless; the two
                                         // folded-cluster routes fire only under --cn.
    std::string minimap_preset = "asm20";
    std::size_t minimap_best_n = 8;
    double ins_dup_min_identity = 0.90;  // identity for an INS to be subtyped DUP
    // Low-complexity tangle guard: a bubble whose interior is a densely-interconnected tangle (many
    // high-degree hub nodes, reached from all over the graph) is a low-complexity region, not a
    // copy-number module. Its node multiplicity is meaningless as copy number, so the coverage/peak
    // DUP routes are suppressed there. A real paralog cluster -- even a large one -- is chain-like
    // (low-degree nodes) and has ~0 hubs, so it is unaffected.
    std::size_t tangle_hub_degree = 20;  // an interior node with >= this many distinct neighbours is a "hub"
    std::size_t tangle_min_hubs = 10;    // a bubble with >= this many hubs is a low-complexity tangle (0=off)
    // Physical implausibility guard for the peak/coverage DUP routes: a copy-number event whose size
    // spans more than this fraction of the whole reference isn't a coherent duplication -- it is a
    // tangle summing diffuse revisits (real peak/coverage DUPs span a few % of the locus). Suppressed.
    // The self-loop REP route (genuine folded tandems) is not gated by this. 0 = off.
    double max_dup_region_frac = 0.80;
    // Refuse a MODULE_BP copy-number call whose integer model fits badly: CN_ROUND_RESIDUAL above this
    // means the calibrated unit does not divide real walks cleanly and the reported integer came from
    // rounding a poor fit. 0 disables. Default OFF, on evidence: the fit is worst exactly where CN is
    // exact against pangene truth (GSTM1 fold residual 0.83, CN 466/466), because the unit is a
    // calibration constant for a heterogeneous paralog module and never claimed to be one real copy.
    double max_cn_model_residual = 0.0;
    std::string gtf_path;                // optional reference-coordinate GTF: annotate variants with the
                                         // genes they touch (INFO GENES), write <prefix>.node_genes.tsv,
                                         // and a per-gene DUP copy-number table; needs a PanSN reference
    std::vector<std::size_t> bubble_ids; // if non-empty, restrict to these bubbles
    bool write_per_bubble_vcf = true;    // also keep each <prefix>.bubble_<id>.vcf
    bool write_variant_nodes = true;     // write <prefix>.variant_nodes.tsv (per-variant node set)
    std::size_t threads = 0;             // worker threads for the per-bubble loop (0 = auto); output is
                                         // identical regardless of thread count
    bool quiet = false;
};

struct VariantCallSummary {
    std::size_t bubbles_seen = 0;
    std::size_t bubbles_with_reference = 0;
    std::size_t bubbles_with_calls = 0;
    std::size_t records_written = 0;
    // Divergent blocks abandoned because the node-token DP would exceed its cell cap. Each is a
    // guaranteed false negative, so it belongs in the summary rather than behind a debug env var.
    std::size_t skipped_large_segments = 0;
    std::size_t del = 0;
    std::size_t ins = 0;
    std::size_t inv = 0;
    std::size_t dup = 0;
    std::size_t multi = 0;  // multiallelic-locus records (--multiallelic-loci)
    std::size_t allele_records = 0;   // bubbles written to the allele VCF
    std::size_t allele_skipped = 0;   // bubbles skipped there (over --allele-vcf-max-bp, or no anchor)
    std::size_t tangle_bubbles = 0;   // bubbles flagged low-complexity tangles (CN routes suppressed)
    std::size_t oversized_dups = 0;   // peak/coverage DUPs suppressed for spanning too much reference
    // MODULE_BP calls refused by --max-cn-model-residual (0 unless that gate is enabled).
    std::size_t declined_cn_model = 0;
};

void call_variants(
    const Graph& graph,
    const VariantCallOptions& options,
    VariantCallSummary* summary_out = nullptr);

} // namespace panvar
