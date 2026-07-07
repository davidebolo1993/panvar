#include "panvar/align.hpp"

#include <algorithm>
#include <string>

#include "edlib.h"

namespace panvar {

// Bit-parallel (Myers) banded fitting alignment via edlib, same contract as the old scalar DP: fit the
// whole query into target from target[0] (free trailing target gap) within an edit budget `band`.
// edlib's SHW/prefix mode is exactly this; k = band caps the edit distance to the same budget, so
// dissimilar pairs bail early (the win on the merge).
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
    if (m == 0) {
        // Whole query is unmatched; identical to the old DP's j=0 column result.
        out.ok = true;
        out.edits = n;
        out.target_consumed = 0;
        out.identity = 0.0;
        return out;
    }

    const EdlibAlignConfig cfg = edlibNewAlignConfig(
        static_cast<int>(band), EDLIB_MODE_SHW, EDLIB_TASK_DISTANCE, nullptr, 0);
    EdlibAlignResult res = edlibAlign(query.data(), static_cast<int>(n),
                                      target.data(), static_cast<int>(m), cfg);

    if (res.status != EDLIB_STATUS_OK || res.editDistance < 0 || res.numLocations == 0) {
        // Optimal path exceeds the edit budget (the old DP's "left the band").
        edlibFreeAlignResult(res);
        out.ok = false;
        return out;
    }

    // edlib returns end locations in ascending order; the smallest matches the old DP,
    // which kept the first (smallest) target position achieving the minimum edit count.
    const std::size_t end_pos = static_cast<std::size_t>(res.endLocations[0]);
    out.ok = true;
    out.edits = static_cast<std::size_t>(res.editDistance);
    out.target_consumed = end_pos + 1;
    const std::size_t denom = std::max(out.query_len, out.target_consumed);
    out.identity = denom == 0 ? 1.0 : 1.0 - static_cast<double>(out.edits) / static_cast<double>(denom);
    edlibFreeAlignResult(res);
    return out;
}

} // namespace panvar
