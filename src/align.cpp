#include "panvar/align.hpp"

#include <algorithm>
#include <limits>
#include <vector>

namespace panvar {

FitAlignResult fit_align(const std::string& query, const std::string& target, std::size_t band) {
    FitAlignResult out;
    out.query_len = query.size();
    const std::size_t n = query.size();
    const std::size_t m = target.size();
    if (n == 0) {
        out.ok = true;
        out.edits = 0;
        out.target_consumed = 0;
        out.identity = 1.0;
        return out;
    }

    constexpr std::size_t kInf = std::numeric_limits<std::size_t>::max() / 4;

    // Banded DP: dp[i][j] = edit distance aligning query[0:i] to target[0:j],
    // with start fixed at target[0] (no leading target gap) and |i-j| <= band.
    std::vector<std::size_t> prev(m + 1, kInf);
    std::vector<std::size_t> cur(m + 1, kInf);

    // i = 0: only j = 0 is reachable (no leading target gap).
    prev[0] = 0;

    for (std::size_t i = 1; i <= n; ++i) {
        std::fill(cur.begin(), cur.end(), kInf);
        const std::size_t jlo = (i > band) ? i - band : 0;
        const std::size_t jhi = std::min(m, i + band);
        for (std::size_t j = jlo; j <= jhi; ++j) {
            std::size_t best = kInf;
            if (j > 0 && prev[j - 1] != kInf) {
                const std::size_t sub = prev[j - 1] + (query[i - 1] == target[j - 1] ? 0 : 1);
                best = std::min(best, sub);
            }
            if (prev[j] != kInf) {
                best = std::min(best, prev[j] + 1); // query base, target gap (deletion in target)
            }
            if (j > 0 && cur[j - 1] != kInf) {
                best = std::min(best, cur[j - 1] + 1); // target base, query gap (insertion in target)
            }
            cur[j] = best;
        }
        std::swap(prev, cur);
    }

    // Final row is in `prev` after the swap; free end gap -> min over reachable j.
    std::size_t best_edits = kInf;
    std::size_t best_j = 0;
    const std::size_t jlo = (n > band) ? n - band : 0;
    const std::size_t jhi = std::min(m, n + band);
    for (std::size_t j = jlo; j <= jhi; ++j) {
        if (prev[j] < best_edits) {
            best_edits = prev[j];
            best_j = j;
        }
    }
    if (best_edits == kInf) {
        out.ok = false;
        return out;
    }
    out.ok = true;
    out.edits = best_edits;
    out.target_consumed = best_j;
    const std::size_t denom = std::max(out.query_len, out.target_consumed);
    out.identity = denom == 0 ? 1.0 : 1.0 - static_cast<double>(best_edits) / static_cast<double>(denom);
    return out;
}

} // namespace panvar
