#include "panvar/minimap2_align.hpp"

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
        const mm_reg1_t& reg = regs[best_idx];
        out.ok = true;
        out.reverse = reg.rev != 0;
        out.query_start_bp = clamp_i32_to_size(reg.qs, query_seq.size());
        out.query_end_bp = clamp_i32_to_size(reg.qe, query_seq.size());
        if (out.query_end_bp < out.query_start_bp) {
            std::swap(out.query_start_bp, out.query_end_bp);
        }
        out.target_start_bp = clamp_i32_to_size(reg.rs, target_seq.size());
        out.target_end_bp = clamp_i32_to_size(reg.re, target_seq.size());
        if (out.target_end_bp < out.target_start_bp) {
            std::swap(out.target_start_bp, out.target_end_bp);
        }
        out.n_matches = reg.mlen > 0 ? static_cast<std::size_t>(reg.mlen) : 0;
        out.aln_block_len = reg.blen > 0 ? static_cast<std::size_t>(reg.blen) : 0;
    }

    for (int i = 0; i < n_regs; ++i) {
        std::free(regs[i].p);
    }
    std::free(regs);
    mm_tbuf_destroy(tbuf);
    mm_idx_destroy(idx);
    return out;
}

} // namespace panvar
