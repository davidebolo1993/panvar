#pragma once

#include <cstddef>
#include <string>

namespace panvar {

// How many bases a change may cross that differ from what it displaces and still count as placed equally
// well. Under the scoring call aligns alleles with (minimap2 defaults: match 2, mismatch 4, long gaps 24 to
// open), each such base costs 6; three cost 18 and a fourth 24, the price of opening another gap, so a
// placement within this many bases of the aligner's is one it could equally have chosen.
constexpr std::size_t kPlacementMismatches = 3;

// The stretch of reference points a change can occupy. Points are the positions between bases, 0-based
// on `ref`. The change is a deletion of ref[s, e), or, with s == e, the insertion of `ins` at point s.
// It is moved left and right one base at a time, as VCF normalisation does, while the haplotype it
// makes differs from the original in at most kPlacementMismatches bases -- none, for an exact repeat.
// A replacement (a deletion and an insertion at once) does not move: its stretch is its own span.
struct PlacementRange {
    std::size_t first = 0, last = 0;
};

inline PlacementRange placement_range(const std::string& ref, std::size_t s, std::size_t e, const std::string& ins,
                                      std::size_t max_mismatch = kPlacementMismatches) {
    if (e > s && !ins.empty()) return {s, e};
    const std::size_t n = e > s ? e - s : ins.size();
    if (n == 0) return {s, e};
    // The haplotype the change makes, base by base: what aligns after the change at a new placement.
    auto hap_at = [&](std::size_t p) -> char {
        if (e > s) return p < s ? ref[p] : ref[p + n];     // deletion
        if (p < s) return ref[p];                          // insertion
        if (p < s + n) return ins[p - s];
        return ref[p - n];
    };
    PlacementRange r{s, e};
    std::size_t mm = 0;
    for (std::size_t y = s; y > 0; --y) {                  // placing it at y - 1
        const char h = e > s ? ref[y - 1 + n] : hap_at(y - 1 + n);
        if (h != ref[y - 1]) { if (++mm > max_mismatch) break; continue; }
        r.first = y - 1;   // only a placement whose last step matched ends the stretch
    }
    mm = 0;
    const std::size_t end_cap = ref.size() - (e > s ? n : 0);
    for (std::size_t y = s; y < end_cap; ++y) {            // placing it at y + 1
        const char h = e > s ? ref[y + n] : hap_at(y);
        if (h != ref[y]) { if (++mm > max_mismatch) break; continue; }
        r.last = (y + 1) + (e > s ? n : 0);
    }
    return r;
}

} // namespace panvar
