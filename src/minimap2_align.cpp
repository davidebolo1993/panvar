#include "panvar/minimap2_align.hpp"

// Thin wrapper over the minimap2 C API (statically linked submodule) for INS-subtype realignment and
// per-gene DUP resolution. Li 2018, https://doi.org/10.1093/bioinformatics/bty191

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <limits>

#include <minimap.h>

namespace panvar {
namespace {

std::size_t clamp_i32_to_size(std::int32_t value, std::size_t cap) {
    if (value <= 0) {
        return 0;
    }
    return std::min<std::size_t>(static_cast<std::size_t>(value), cap);
}

Minimap2Hit reg_to_hit(const mm_reg1_t& reg, std::size_t qlen, std::size_t tlen) {
    Minimap2Hit out;
    out.ok = true;
    out.reverse = reg.rev != 0;
    out.query_start_bp = clamp_i32_to_size(reg.qs, qlen);
    out.query_end_bp = clamp_i32_to_size(reg.qe, qlen);
    if (out.query_end_bp < out.query_start_bp) std::swap(out.query_start_bp, out.query_end_bp);
    out.target_start_bp = clamp_i32_to_size(reg.rs, tlen);
    out.target_end_bp = clamp_i32_to_size(reg.re, tlen);
    if (out.target_end_bp < out.target_start_bp) std::swap(out.target_start_bp, out.target_end_bp);
    out.n_matches = reg.mlen > 0 ? static_cast<std::size_t>(reg.mlen) : 0;
    out.aln_block_len = reg.blen > 0 ? static_cast<std::size_t>(reg.blen) : 0;
    // Gap-compressed identity from the (EQX) cigar: matches / (matches + mismatches + indel events).
    // Each I/D run is one event regardless of length, so a large structural indel barely affects it.
    if (reg.p != nullptr && reg.p->n_cigar > 0) {
        std::size_t m = 0, x = 0, gaps = 0;
        out.cigar.reserve(reg.p->n_cigar);
        for (uint32_t i = 0; i < reg.p->n_cigar; ++i) {
            const uint32_t op = reg.p->cigar[i] & 0xf;
            const uint32_t len = reg.p->cigar[i] >> 4;
            out.cigar.emplace_back(static_cast<int>(op), static_cast<int>(len));
            if (op == 7) m += len;          // '=' match
            else if (op == 8) x += len;     // 'X' mismatch
            else if (op == 1 || op == 2) ++gaps;  // 'I'/'D' indel run (one event)
        }
        const std::size_t denom = m + x + gaps;
        out.gc_identity = denom == 0 ? 0.0 : static_cast<double>(m) / static_cast<double>(denom);
    } else {
        out.gc_identity = out.identity();  // no cigar: fall back to block identity
    }
    return out;
}

} // namespace

Minimap2Hit minimap2_best_hit(
    const std::string& query_name,
    const std::string& query_seq,
    const std::string& target_name,
    const std::string& target_seq,
    const std::string& preset,
    std::size_t best_n) {

    Minimap2Hit out;
    if (query_seq.empty() || target_seq.empty()) {
        return out;
    }
    if (query_seq.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
        target_seq.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        return out;
    }

    mm_idxopt_t idx_opt;
    mm_mapopt_t map_opt;
    if (mm_set_opt(nullptr, &idx_opt, &map_opt) < 0) {
        return out;
    }
    const std::string chosen = preset.empty() ? std::string("asm20") : preset;
    if (mm_set_opt(chosen.c_str(), &idx_opt, &map_opt) < 0) {
        return out;
    }
    map_opt.flag |= (MM_F_CIGAR | MM_F_EQX);
    map_opt.flag &= ~MM_F_OUT_SAM;
    const int requested_best_n = static_cast<int>(
        std::min<std::size_t>(
            static_cast<std::size_t>(std::numeric_limits<int>::max()),
            std::max<std::size_t>(static_cast<std::size_t>(1), best_n)));
    map_opt.best_n = std::max(map_opt.best_n, requested_best_n);

    const char* seq_ptrs[1] = {target_seq.c_str()};
    const char* name_ptrs[1] = {target_name.c_str()};
    mm_idx_t* idx = mm_idx_str(
        idx_opt.w,
        idx_opt.k,
        (idx_opt.flag & MM_I_HPC) != 0 ? 1 : 0,
        idx_opt.bucket_bits,
        1,
        seq_ptrs,
        name_ptrs);
    if (idx == nullptr) {
        return out;
    }
    mm_mapopt_update(&map_opt, idx);

    mm_tbuf_t* tbuf = mm_tbuf_init();
    if (tbuf == nullptr) {
        mm_idx_destroy(idx);
        return out;
    }

    int n_regs = 0;
    mm_reg1_t* regs = mm_map(
        idx,
        static_cast<int>(query_seq.size()),
        query_seq.c_str(),
        &n_regs,
        tbuf,
        &map_opt,
        query_name.c_str());
    if (regs == nullptr || n_regs <= 0) {
        mm_tbuf_destroy(tbuf);
        mm_idx_destroy(idx);
        return out;
    }

    int best_idx = -1;
    for (int i = 0; i < n_regs; ++i) {
        if (best_idx < 0) {
            best_idx = i;
            continue;
        }
        const bool curr_primary = regs[i].parent == regs[i].id;
        const bool best_primary = regs[best_idx].parent == regs[best_idx].id;
        if (curr_primary != best_primary) {
            if (curr_primary) {
                best_idx = i;
            }
            continue;
        }
        if (regs[i].score > regs[best_idx].score ||
            (regs[i].score == regs[best_idx].score && regs[i].blen > regs[best_idx].blen) ||
            (regs[i].score == regs[best_idx].score &&
             regs[i].blen == regs[best_idx].blen &&
             regs[i].mlen > regs[best_idx].mlen)) {
            best_idx = i;
        }
    }

    if (best_idx >= 0) {
        out = reg_to_hit(regs[best_idx], query_seq.size(), target_seq.size());
    }

    for (int i = 0; i < n_regs; ++i) {
        std::free(regs[i].p);
    }
    std::free(regs);
    mm_tbuf_destroy(tbuf);
    mm_idx_destroy(idx);
    return out;
}

std::vector<Minimap2Hit> minimap2_hits(
    const std::string& query_name,
    const std::string& query_seq,
    const std::string& target_name,
    const std::string& target_seq,
    const std::string& preset,
    std::size_t best_n) {

    std::vector<Minimap2Hit> hits;
    if (query_seq.empty() || target_seq.empty()) return hits;
    if (query_seq.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
        target_seq.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        return hits;
    }

    mm_idxopt_t idx_opt;
    mm_mapopt_t map_opt;
    if (mm_set_opt(nullptr, &idx_opt, &map_opt) < 0) return hits;
    const std::string chosen = preset.empty() ? std::string("asm20") : preset;
    if (mm_set_opt(chosen.c_str(), &idx_opt, &map_opt) < 0) return hits;
    map_opt.flag |= (MM_F_CIGAR | MM_F_EQX);
    map_opt.flag &= ~MM_F_OUT_SAM;
    // Keep secondary alignments (one per paralog copy in the target) and raise the cap.
    const int want = static_cast<int>(std::min<std::size_t>(
        static_cast<std::size_t>(std::numeric_limits<int>::max()),
        std::max<std::size_t>(static_cast<std::size_t>(1), best_n)));
    map_opt.best_n = std::max(map_opt.best_n, want);
    // Keep moderately low-scoring secondaries, so a partial or truncated paralog copy still surfaces as
    // its own hit -- without admitting so many weak hits that they steal loci in the competition and
    // blur the separation between neighbouring paralogs.
    map_opt.pri_ratio = 0.5f;

    const char* seq_ptrs[1] = {target_seq.c_str()};
    const char* name_ptrs[1] = {target_name.c_str()};
    mm_idx_t* idx = mm_idx_str(idx_opt.w, idx_opt.k,
                               (idx_opt.flag & MM_I_HPC) != 0 ? 1 : 0,
                               idx_opt.bucket_bits, 1, seq_ptrs, name_ptrs);
    if (idx == nullptr) return hits;
    mm_mapopt_update(&map_opt, idx);

    mm_tbuf_t* tbuf = mm_tbuf_init();
    if (tbuf == nullptr) { mm_idx_destroy(idx); return hits; }

    int n_regs = 0;
    mm_reg1_t* regs = mm_map(idx, static_cast<int>(query_seq.size()), query_seq.c_str(),
                             &n_regs, tbuf, &map_opt, query_name.c_str());
    for (int i = 0; i < n_regs; ++i) {
        hits.push_back(reg_to_hit(regs[i], query_seq.size(), target_seq.size()));
        std::free(regs[i].p);
    }
    std::free(regs);
    mm_tbuf_destroy(tbuf);
    mm_idx_destroy(idx);
    return hits;
}

} // namespace panvar
