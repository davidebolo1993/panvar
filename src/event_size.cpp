#include "panvar/event_size.hpp"

#include <ksw2.h>

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cstring>

namespace panvar {

// The changes an alignment of hap to ref makes, from its CIGAR ((op, length) pairs, op 0 aligned columns
// compared base by base, 7 '=', 8 'X', 1 hap bases, 2 ref bases). The changes are the gaps. Two gaps are
// one change, a run replaced, unless the sequence between them aligns well enough to pay for opening a
// new long gap: some stretch of it scoring 24 (the long-gap open cost) at 2 per match and -4 per
// mismatch, as twelve matching bases do. Unrelated sequence matches here and there by chance, and those
// bases do not split a replaced run; a mismatch outside a change is a small variant and not reported.
// False if the CIGAR does not cover both.
bool changes_from_cigar(const std::vector<std::pair<int, int>>& cigar, const std::string& ref,
                        const std::string& hap, std::vector<AlignedChange>& out) {
    out.clear();
    std::size_t ti = 0, qi = 0;
    bool open = false;
    AlignedChange cur;
    auto same = [&](std::size_t a, std::size_t b) {
        const char x = static_cast<char>(ref[a] & 0xDF), y = static_cast<char>(hap[b] & 0xDF);
        return x == y && x != 'N';
    };
    long long run = 0;   // best-scoring stretch ending here since the last gap (kept at >= 0)
    for (const auto& [op, len_i] : cigar) {
        const std::size_t len = static_cast<std::size_t>(len_i);
        if (op == 0 || op == 7 || op == 8) {             // aligned columns
            if (ti + len > ref.size() || qi + len > hap.size()) return false;
            for (std::size_t c = 0; c < len; ++c, ++ti, ++qi) {
                if (!open) continue;
                run = std::max(0LL, run + ((op != 8 && same(ti, qi)) ? 2 : -4));
                if (run >= 24) { out.push_back(cur); open = false; }
            }
        } else if (op == 1 || op == 2) {                 // 1: hap bases (insertion), 2: ref bases
            if (!open) { cur = {ti, ti, qi, qi}; open = true; }
            if (op == 1) qi += len; else ti += len;
            cur.t1 = ti; cur.q1 = qi;                     // through the sequence since the last gap
            run = 0;
        } else {
            return false;
        }
    }
    if (ti != ref.size() || qi != hap.size()) { out.clear(); return false; }
    if (open) out.push_back(cur);
    return true;
}

// Global alignment of hap to ref, banded to w diagonals either side, with minimap2's default scoring:
// match 2, mismatch 4, and two gap costs, 4 + 2 per base and 24 + 1 per base, the cheaper applying, so
// one long indel costs less than the same bases scattered. False when the alignment failed.
bool aligned_changes(const std::string& ref, const std::string& hap, int w, std::vector<AlignedChange>& out) {
    out.clear();
    if (ref.empty() || hap.empty()) {
        if (!ref.empty() || !hap.empty()) out.push_back({0, ref.size(), 0, hap.size()});
        return true;
    }
    auto code = [](char c) -> std::uint8_t {
        switch (c & 0xDF) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return 4;
        }
    };
    std::vector<std::uint8_t> t(ref.size()), q(hap.size());
    for (std::size_t i = 0; i < ref.size(); ++i) t[i] = code(ref[i]);
    for (std::size_t i = 0; i < hap.size(); ++i) q[i] = code(hap[i]);
    std::int8_t mat[25];
    for (int i = 0; i < 5; ++i)
        for (int j = 0; j < 5; ++j) mat[i * 5 + j] = static_cast<std::int8_t>((i == 4 || j == 4) ? -1 : i == j ? 2 : -4);
    std::vector<std::pair<int, int>> cigar;
    auto align = [&](int band, int& score) {
        ksw_extz_t ez;
        std::memset(&ez, 0, sizeof(ez));
        ksw_extd2_sse(nullptr, static_cast<int>(q.size()), q.data(), static_cast<int>(t.size()), t.data(), 5, mat,
                      4, 2, 24, 1, band, -1, 0, 0, &ez);
        if (ez.n_cigar <= 0 || ez.cigar == nullptr) { std::free(ez.cigar); return false; }
        cigar.clear();
        for (int k = 0; k < ez.n_cigar; ++k)
            cigar.emplace_back(static_cast<int>(ez.cigar[k] & 0xf), static_cast<int>(ez.cigar[k] >> 4));
        std::free(ez.cigar);
        score = ez.score;
        return true;
    };
    int score = 0;
    if (!align(w, score)) return false;
    // The band can hold a worse placement than one outside it. An alignment reaching D diagonals from the main
    // one has at least 2D - |n - m| gap bases, each costing 1 or more and taking a column of 2 from the matches,
    // so it scores at most n + m - 2(2D - |n - m|). Where that bound for D = w + 1 is above the banded score,
    // the alignment is redone with the band beyond which the bound falls under it, and kept if it scores more.
    if (w >= 0) {
        const long long n = static_cast<long long>(ref.size()), m = static_cast<long long>(hap.size());
        const long long d = n > m ? n - m : m - n;
        const long long need = std::min(std::max(n, m), (n + m + 2 * d - static_cast<long long>(score) + 3) / 4 - 1);
        const std::size_t cells = static_cast<std::size_t>(std::max(n, m)) * static_cast<std::size_t>(2 * need + 1);
        if (need > w && cells <= kSeqAlignCellCap) {
            const std::vector<std::pair<int, int>> banded = cigar;
            int wide = 0;
            if (!align(static_cast<int>(need), wide) || wide <= score) cigar = banded;
        }
    }
    return changes_from_cigar(cigar, ref, hap, out);
}

// See the header: the largest single change of the aligned core, by its larger side.
std::size_t change_size_bp(const std::string& r, const std::string& h) {
    const std::size_t whole = std::max(r.size(), h.size());
    if (whole > kSubEditMaxBp) return whole;
    std::size_t p = 0, sfx = 0;
    while (p < r.size() && p < h.size() && r[p] == h[p]) ++p;
    while (sfx < r.size() - p && sfx < h.size() - p && r[r.size() - 1 - sfx] == h[h.size() - 1 - sfx]) ++sfx;
    const std::size_t n = r.size() - p - sfx, m = h.size() - p - sfx;
    if (n == 0 || m == 0) return std::max(n, m);
    const std::size_t band = (n > m ? n - m : m - n) + 64;
    if (std::max(n, m) * (2 * band + 1) > kSeqAlignCellCap) return whole;
    const std::string rc = r.substr(p, n), hc = h.substr(p, m);
    std::vector<AlignedChange> changes;
    if (!aligned_changes(rc, hc, static_cast<int>(band), changes)) return whole;
    std::size_t best = 0;
    for (AlignedChange c : changes) {
        while (c.t0 < c.t1 && c.q0 < c.q1 && rc[c.t0] == hc[c.q0]) { ++c.t0; ++c.q0; }
        while (c.t1 > c.t0 && c.q1 > c.q0 && rc[c.t1 - 1] == hc[c.q1 - 1]) { --c.t1; --c.q1; }
        best = std::max(best, std::max(c.t1 - c.t0, c.q1 - c.q0));
    }
    return best;
}

} // namespace panvar
