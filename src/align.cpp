#include "panvar/align.hpp"

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <string>
#include <cstdint>

#include "edlib.h"

namespace panvar {
namespace {
// Bases compare case-insensitively: a soft-masked (lower-case) graph has the same sequence as an upper-case
// one, and graph_sort writes the nodes it flips in upper case, so a case-sensitive compare would fail a
// lower-case copy of a repeat unit against an upper-case one.
const EdlibEqualityPair kCaseless[] = {{'A', 'a'}, {'C', 'c'}, {'G', 'g'}, {'T', 't'}, {'N', 'n'}};
constexpr int kCaselessN = sizeof(kCaseless) / sizeof(kCaseless[0]);
} // namespace

// Bit-parallel (Myers) fitting alignment via edlib: fit the whole query into target from target[0]
// (free trailing target gap) within an edit budget `band`. edlib's SHW/prefix mode is exactly this;
// k = band caps the edit distance, so dissimilar pairs stop early.
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
        // The whole query is unmatched.
        out.ok = true;
        out.edits = n;
        out.target_consumed = 0;
        out.identity = 0.0;
        return out;
    }

    const EdlibAlignConfig cfg = edlibNewAlignConfig(
        static_cast<int>(band), EDLIB_MODE_SHW, EDLIB_TASK_DISTANCE, kCaseless, kCaselessN);
    EdlibAlignResult res = edlibAlign(query.data(), static_cast<int>(n),
                                      target.data(), static_cast<int>(m), cfg);

    if (res.status != EDLIB_STATUS_OK || res.editDistance < 0 || res.numLocations == 0) {
        // The best fit needs more edits than the budget.
        edlibFreeAlignResult(res);
        out.ok = false;
        return out;
    }

    // edlib returns end locations in ascending order; keep the smallest target position
    // achieving the minimum edit count.
    const std::size_t end_pos = static_cast<std::size_t>(res.endLocations[0]);
    out.ok = true;
    out.edits = static_cast<std::size_t>(res.editDistance);
    out.target_consumed = end_pos + 1;
    const std::size_t denom = std::max(out.query_len, out.target_consumed);
    out.identity = denom == 0 ? 1.0 : 1.0 - static_cast<double>(out.edits) / static_cast<double>(denom);
    edlibFreeAlignResult(res);
    return out;
}

std::size_t nw_edit_ops(const std::string& target, const std::string& query, std::string& ops) {
    ops.clear();
    if (target.empty() || query.empty()) {
        ops.assign(target.size(), 'D');
        ops.append(query.size(), 'I');
        return target.size() + query.size();
    }
    const EdlibAlignConfig cfg = edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, kCaseless, kCaselessN);
    // edlib aligns `query` against `target`: its INSERT is a base of the query only, DELETE of the target.
    EdlibAlignResult res = edlibAlign(query.data(), static_cast<int>(query.size()),
                                      target.data(), static_cast<int>(target.size()), cfg);
    if (res.status != EDLIB_STATUS_OK || res.editDistance < 0 || res.alignment == nullptr) {
        edlibFreeAlignResult(res);
        throw std::runtime_error("edlib failed to align two sequences of " + std::to_string(target.size()) +
                                 " and " + std::to_string(query.size()) + " bp");
    }
    ops.reserve(static_cast<std::size_t>(res.alignmentLength));
    for (int k = 0; k < res.alignmentLength; ++k) {
        switch (res.alignment[k]) {
            case EDLIB_EDOP_MATCH: ops.push_back('M'); break;
            case EDLIB_EDOP_MISMATCH: ops.push_back('X'); break;
            case EDLIB_EDOP_INSERT: ops.push_back('I'); break;
            default: ops.push_back('D'); break;
        }
    }
    const std::size_t d = static_cast<std::size_t>(res.editDistance);
    edlibFreeAlignResult(res);
    return d;
}

std::size_t infix_distance(const std::string& query, const std::string& target, std::size_t max_k) {
    if (query.empty()) return 0;
    if (target.empty()) return query.size() <= max_k ? query.size() : SIZE_MAX;
    const int k = max_k < static_cast<std::size_t>(std::numeric_limits<int>::max()) ? static_cast<int>(max_k) : -1;
    const EdlibAlignConfig cfg = edlibNewAlignConfig(k, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, kCaseless, kCaselessN);
    EdlibAlignResult res = edlibAlign(query.data(), static_cast<int>(query.size()),
                                      target.data(), static_cast<int>(target.size()), cfg);
    const std::size_t d = (res.status == EDLIB_STATUS_OK && res.editDistance >= 0)
                              ? static_cast<std::size_t>(res.editDistance) : SIZE_MAX;
    edlibFreeAlignResult(res);
    return d;
}

double infix_identity(const std::string& query, const std::string& target, std::size_t max_k, std::size_t& start,
                      std::size_t& end) {
    start = end = 0;
    if (query.empty() || target.empty()) return -1.0;
    const int k = max_k < static_cast<std::size_t>(std::numeric_limits<int>::max()) ? static_cast<int>(max_k) : -1;
    const EdlibAlignConfig cfg = edlibNewAlignConfig(k, EDLIB_MODE_HW, EDLIB_TASK_PATH, kCaseless, kCaselessN);
    EdlibAlignResult res = edlibAlign(query.data(), static_cast<int>(query.size()),
                                      target.data(), static_cast<int>(target.size()), cfg);
    double identity = -1.0;
    // A placement ending before it starts matched no target base: no placement.
    if (res.status == EDLIB_STATUS_OK && res.editDistance >= 0 && res.numLocations > 0 && res.alignment != nullptr &&
        res.startLocations[0] >= 0 && res.endLocations[0] >= res.startLocations[0]) {
        start = static_cast<std::size_t>(res.startLocations[0]);
        end = static_cast<std::size_t>(res.endLocations[0]);
        std::size_t columns = 0, differences = 0;
        unsigned char prev = EDLIB_EDOP_MATCH;
        for (int i = 0; i < res.alignmentLength; ++i) {
            const unsigned char op = res.alignment[i];
            if (op == EDLIB_EDOP_MATCH) ++columns;
            else if (op == EDLIB_EDOP_MISMATCH) { ++columns; ++differences; }
            else if (op != prev) { ++columns; ++differences; }   // a gap counts once, at its first base
            prev = op;
        }
        identity = columns ? 1.0 - static_cast<double>(differences) / static_cast<double>(columns) : -1.0;
    }
    edlibFreeAlignResult(res);
    return identity;
}

std::size_t nw_edit_ops_bounded(const std::string& target, const std::string& query, std::size_t max_k,
                                std::string& ops) {
    ops.clear();
    const std::size_t len_gap = target.size() > query.size() ? target.size() - query.size() : query.size() - target.size();
    if (len_gap > max_k) return SIZE_MAX;
    if (target.empty() || query.empty()) {
        ops.assign(target.size(), 'D');
        ops.append(query.size(), 'I');
        return target.size() + query.size();
    }
    const int k = max_k < static_cast<std::size_t>(std::numeric_limits<int>::max()) ? static_cast<int>(max_k) : -1;
    const EdlibAlignConfig cfg = edlibNewAlignConfig(k, EDLIB_MODE_NW, EDLIB_TASK_PATH, kCaseless, kCaselessN);
    EdlibAlignResult res = edlibAlign(query.data(), static_cast<int>(query.size()),
                                      target.data(), static_cast<int>(target.size()), cfg);
    if (res.status != EDLIB_STATUS_OK || res.editDistance < 0 || res.alignment == nullptr) {
        edlibFreeAlignResult(res);
        return SIZE_MAX;
    }
    ops.reserve(static_cast<std::size_t>(res.alignmentLength));
    for (int i = 0; i < res.alignmentLength; ++i) {
        switch (res.alignment[i]) {
            case EDLIB_EDOP_MATCH: ops.push_back('M'); break;
            case EDLIB_EDOP_MISMATCH: ops.push_back('X'); break;
            case EDLIB_EDOP_INSERT: ops.push_back('I'); break;
            default: ops.push_back('D'); break;
        }
    }
    const std::size_t d = static_cast<std::size_t>(res.editDistance);
    edlibFreeAlignResult(res);
    return d;
}

std::size_t nw_distance(const std::string& a, const std::string& b, std::size_t max_k) {
    const bool bounded = max_k < static_cast<std::size_t>(std::numeric_limits<int>::max());
    if (a.empty() || b.empty()) {
        const std::size_t d = std::max(a.size(), b.size());
        return bounded && d > max_k ? max_k + 1 : d;
    }
    const std::size_t len_gap = a.size() > b.size() ? a.size() - b.size() : b.size() - a.size();
    if (bounded && len_gap > max_k) return max_k + 1;
    const EdlibAlignConfig cfg = edlibNewAlignConfig(bounded ? static_cast<int>(max_k) : -1, EDLIB_MODE_NW,
                                                     EDLIB_TASK_DISTANCE, kCaseless, kCaselessN);
    EdlibAlignResult res = edlibAlign(a.data(), static_cast<int>(a.size()),
                                      b.data(), static_cast<int>(b.size()), cfg);
    const std::size_t d = (res.status == EDLIB_STATUS_OK && res.editDistance >= 0)
                              ? static_cast<std::size_t>(res.editDistance)
                              : bounded ? max_k + 1 : std::max(a.size(), b.size());
    edlibFreeAlignResult(res);
    return d;
}

} // namespace panvar
