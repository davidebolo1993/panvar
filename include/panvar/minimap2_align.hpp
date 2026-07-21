#pragma once

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace panvar {

// Minimal best-alignment summary from minimap2, used for INS subtype refinement
// (does an inserted sequence map back to a local reference window?).
struct Minimap2Hit {
    bool ok = false;
    bool reverse = false;
    std::size_t query_start_bp = 0;
    std::size_t query_end_bp = 0;
    std::size_t target_start_bp = 0;
    std::size_t target_end_bp = 0;
    std::size_t n_matches = 0;     // mlen: matching bases in the alignment
    std::size_t aln_block_len = 0; // blen: alignment block length (match+mismatch+indel)
    // Gap-compressed identity (each indel run counts as ONE event, not by length), like minimap2's
    // `de`. Robust when a long query maps a copy missing a large insertion (e.g. C4 short form lacks
    // a large internal indel): the structural indel does not deflate identity. 0 when no cigar available.
    double gc_identity = 0.0;

    // Extended (EQX) CIGAR of the alignment as (op, length) pairs, op in {7='=', 8='X', 1='I', 2='D'},
    // walking query_start_bp/target_start_bp forward. Empty when no cigar was produced. Used to locate the
    // divergent columns between two near-identical paralog CDS (per-site copy-number resolution).
    std::vector<std::pair<int, int>> cigar;

    // Fraction of the alignment block that matched; 0 when there is no alignment.
    double identity() const {
        return aln_block_len == 0 ? 0.0
                                  : static_cast<double>(n_matches) / static_cast<double>(aln_block_len);
    }
};

// Align query_seq against a single target_seq and return the best primary hit.
// preset is a minimap2 asm preset (asm5|asm10|asm20). Returns ok=false on any
// failure or when either sequence is empty / too large.
Minimap2Hit minimap2_best_hit(
    const std::string& query_name,
    const std::string& query_seq,
    const std::string& target_name,
    const std::string& target_seq,
    const std::string& preset = "asm20",
    std::size_t best_n = 8);

// Align query_seq against target_seq and return ALL hit regions (primary + secondary),
// each as a Minimap2Hit. Used to count paralog copies (one hit per copy in the target).
// preset/best_n as above. Returns empty on failure or empty/oversized input.
std::vector<Minimap2Hit> minimap2_hits(
    const std::string& query_name,
    const std::string& query_seq,
    const std::string& target_name,
    const std::string& target_seq,
    const std::string& preset = "asm20",
    std::size_t best_n = 64);

} // namespace panvar
