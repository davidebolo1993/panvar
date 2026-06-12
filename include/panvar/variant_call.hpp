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
    std::size_t rescue_min_bp = 0;       // floor for events kept for merge/rescue (0 -> min_sv_bp/2)
    bool classify_ins = false;           // minimap2 INS subtype refinement (NOVEL vs DUP)
    bool cn_from_multiplicity = false;    // emit DUP from peak node multiplicity for folded bubbles
                                         // with no self-loop (e.g. GSTM1) that panphorte could not collapse
    std::string minimap_preset = "asm20";
    std::size_t minimap_best_n = 8;
    double ins_dup_min_identity = 0.90;  // identity for an INS to be subtyped DUP
    std::vector<std::size_t> bubble_ids; // if non-empty, restrict to these bubbles
    bool write_per_bubble_vcf = true;    // also keep each <prefix>.bubble_<id>.vcf
    bool write_variant_paths = true;     // write <prefix>.variant_paths.tsv (per-variant carrier sub-walks)
    bool quiet = false;
};

struct VariantCallSummary {
    std::size_t bubbles_seen = 0;
    std::size_t bubbles_with_reference = 0;
    std::size_t bubbles_with_calls = 0;
    std::size_t records_written = 0;
    std::size_t del = 0;
    std::size_t ins = 0;
    std::size_t inv = 0;
    std::size_t dup = 0;
};

void call_variants(
    const Graph& graph,
    const VariantCallOptions& options,
    VariantCallSummary* summary_out = nullptr);

} // namespace panvar
