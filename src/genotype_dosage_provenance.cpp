#include "panvar/genotype_dosage_repair.hpp"

#include "edlib.h"

#include <algorithm>
#include <cstdlib>
#include <map>
#include <sstream>
#include <stdexcept>

namespace panvar::dosage {

namespace {

long edit_distance_bounded(const std::string& q, const std::string& t, long bound) {
    const EdlibAlignConfig cfg = edlibNewAlignConfig(static_cast<int>(bound), EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, nullptr, 0);
    EdlibAlignResult r = edlibAlign(q.data(), static_cast<int>(q.size()), t.data(), static_cast<int>(t.size()), cfg);
    if (r.status != EDLIB_STATUS_OK) { edlibFreeAlignResult(r); throw std::runtime_error("edlib failed"); }
    const long d = r.editDistance;
    edlibFreeAlignResult(r);
    return d;   // -1 when above the bound
}

} // namespace

std::string panel_digest(const std::vector<PanelSequence>& panel) {
    std::vector<std::pair<std::string, std::string>> rows;
    for (const PanelSequence& p : panel) rows.push_back({p.id, sha256_hex(p.sequence)});
    std::sort(rows.begin(), rows.end());
    std::string text;
    for (const auto& [id, h] : rows) text += id + "\t" + h + "\n";
    return sha256_hex(text);
}

std::vector<NearestPanel> nearest_panel_sequences(const std::string& query, const std::vector<PanelSequence>& panel) {
    // canonical order first, so the search itself does not depend on the panel's order
    std::vector<std::pair<std::string, const PanelSequence*>> order;
    for (const PanelSequence& p : panel) order.push_back({sha256_hex(p.sequence), &p});
    std::sort(order.begin(), order.end(), [](const auto& a, const auto& b) {
        return a.first != b.first ? a.first < b.first : a.second->id < b.second->id;
    });
    long best = -1;
    std::vector<std::pair<long, std::size_t>> dist(order.size(), {-1, 0});
    std::map<std::string, long> by_sequence;   // identical sequences share one distance
    for (std::size_t i = 0; i < order.size(); ++i) {
        const auto memo = by_sequence.find(order[i].first);
        long d;
        if (memo != by_sequence.end()) d = memo->second;
        else {
            d = edit_distance_bounded(query, order[i].second->sequence, best);   // bound -1: unbounded
            by_sequence[order[i].first] = d;
        }
        if (d < 0) continue;
        dist[i] = {d, 0};
        if (best < 0 || d < best) best = d;
    }
    std::vector<NearestPanel> out;
    std::map<std::string, NearestPanel> aligned;   // per sequence digest
    for (std::size_t i = 0; i < order.size(); ++i) {
        if (dist[i].first != best) continue;
        auto it = aligned.find(order[i].first);
        if (it == aligned.end()) {
            const PanelSequence& p = *order[i].second;
            const EdlibAlignConfig cfg = edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, nullptr, 0);
            EdlibAlignResult r = edlibAlign(query.data(), static_cast<int>(query.size()), p.sequence.data(), static_cast<int>(p.sequence.size()), cfg);
            if (r.status != EDLIB_STATUS_OK) { edlibFreeAlignResult(r); throw std::runtime_error("edlib path failed"); }
            NearestPanel n;
            n.sha256 = order[i].first;
            n.distance = r.editDistance;
            n.alignment_columns = static_cast<std::size_t>(r.alignmentLength);
            n.normalized = r.alignmentLength > 0 ? static_cast<double>(r.editDistance) / r.alignmentLength : 0.0;
            char* cig = edlibAlignmentToCigar(r.alignment, r.alignmentLength, EDLIB_CIGAR_EXTENDED);
            n.cigar = cig ? cig : "";
            std::free(cig);
            edlibFreeAlignResult(r);
            if (n.distance != best) throw std::runtime_error("edlib path distance disagrees with the distance pass");
            it = aligned.emplace(order[i].first, n).first;
        }
        NearestPanel n = it->second;
        n.id = order[i].second->id;
        out.push_back(std::move(n));
    }
    return out;
}

} // namespace panvar::dosage
