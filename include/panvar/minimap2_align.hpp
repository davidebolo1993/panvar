#pragma once

#include <cstddef>
#include <string>

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

} // namespace panvar
