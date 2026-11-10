#include "panvar/bubbles.hpp"

#include "panvar/bubble_path.hpp"

#include "panvar/cli_utils.hpp"
#include "panvar/placement.hpp"
#include "panvar/path_inversion.hpp"

#include <algorithm>
#include <cctype>
#include <atomic>
#include <cstdlib>
#include <iostream>
#include <cstdint>
#include <fstream>
#include <functional>
#include <limits>
#include <map>
#include <numeric>
#include <optional>
#include <set>
#include <queue>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace panvar {
namespace {

// Compact node id <-> index mapping (indices key the per-endpoint dedup map).
struct PackedGraph {
    std::vector<std::string> node_ids;
    std::unordered_map<std::string, std::uint32_t> node_idx_of;
};

struct BubbleCandidateIdx {
    std::uint32_t source = 0;
    std::uint32_t sink = 0;
    std::vector<std::uint32_t> inside;
    // The sides the interior hangs off (true: the node's end side), as the decomposition gave them or the
    // walks voted them.
    bool source_end = true;
    bool sink_end = false;
};

struct EndpointKey {
    std::uint32_t a = 0;
    std::uint32_t b = 0;

    bool operator==(const EndpointKey& other) const {
        return a == other.a && b == other.b;
    }
};

struct EndpointKeyHash {
    std::size_t operator()(const EndpointKey& key) const {
        const std::uint64_t packed = (static_cast<std::uint64_t>(key.a) << 32U) | key.b;
        return std::hash<std::uint64_t>{}(packed);
    }
};

using CandidateMap = std::unordered_map<EndpointKey, BubbleCandidateIdx, EndpointKeyHash>;

// The bubble-path index and walk extraction are the ones every other module uses (bubble_path.hpp), so
// scoring sees the same alleles they do, a direct source->sink deletion included.
using PathIndex = BubblePathIndex;

// Endpoint pairs whose graph-derived interior hit the traversal cap and fell back to the path-derived one.
// Reported, since such an interior holds only what this panel walks.
std::atomic<std::size_t> g_interior_traversal_truncated{0};

// Defined below with the rest of the interior code; declared here because merging recomputes a fused
// bubble's interior with the rule discovery uses.
std::unordered_set<std::uint32_t> interior_indices_for_pair(
    const Graph& graph,
    const PackedGraph& packed,
    const std::vector<PathIndex>& path_indexes,
    const std::string& source_id,
    const std::string& sink_id,
    std::uint32_t src_idx,
    std::uint32_t sink_idx,
    const SnarlInward* fixed = nullptr,
    SnarlInward* used = nullptr);
bool reachable_interior(
    const Graph& graph,
    const PackedGraph& packed,
    std::uint64_t start_handle,
    std::uint32_t stop_idx,
    std::uint32_t origin_idx,
    bool backward,
    std::unordered_set<std::uint32_t>& out);

PackedGraph pack_graph(const Graph& graph) {
    PackedGraph packed;
    packed.node_ids.reserve(graph.nodes.size());
    for (const auto& [node_id, _node] : graph.nodes) {
        packed.node_ids.push_back(node_id);
    }
    std::sort(packed.node_ids.begin(), packed.node_ids.end());

    packed.node_idx_of.reserve(packed.node_ids.size() * 2);
    for (std::size_t i = 0; i < packed.node_ids.size(); ++i) {
        packed.node_idx_of[packed.node_ids[i]] = static_cast<std::uint32_t>(i);
    }
    return packed;
}

EndpointKey normalize_endpoint_key(std::uint32_t source, std::uint32_t sink) {
    if (source < sink) {
        return EndpointKey{source, sink};
    }
    return EndpointKey{sink, source};
}

// Keep one candidate per unordered endpoint pair, preferring the larger inside set. Two snarls can share both
// boundary nodes on different sides (a cycle through both); downstream a site is identified by its boundary
// nodes, so one is kept and the other's interior counted in `collapsed` for the run to report. Of two inside
// sets of one size, the one the reference crosses (`ref_crosses`, when given): its alleles are the ones `call`
// can place.
bool add_or_replace_candidate(CandidateMap& dedup, BubbleCandidateIdx candidate, std::size_t* collapsed = nullptr,
                              const std::function<bool(const BubbleCandidateIdx&)>& ref_crosses = nullptr) {
    const EndpointKey key = normalize_endpoint_key(candidate.source, candidate.sink);
    const auto it = dedup.find(key);
    if (it != dedup.end() && collapsed != nullptr && it->second.inside != candidate.inside) ++*collapsed;
    if (it == dedup.end() || candidate.inside.size() > it->second.inside.size() ||
        (candidate.inside.size() == it->second.inside.size() && candidate.inside != it->second.inside &&
         ref_crosses && ref_crosses(candidate) && !ref_crosses(it->second))) {
        dedup[key] = std::move(candidate);
        return true;
    }
    return false;
}

std::string endpoint_key(const std::string& a, const std::string& b) {
    if (a < b) {
        return a + "|" + b;
    }
    return b + "|" + a;
}

bool is_unsigned_decimal(const std::string& s) {
    if (s.empty()) {
        return false;
    }
    for (const char c : s) {
        if (c < '0' || c > '9') {
            return false;
        }
    }
    return true;
}

std::string_view strip_leading_zeros(std::string_view s) {
    std::size_t i = 0;
    while (i + 1 < s.size() && s[i] == '0') {
        ++i;
    }
    return s.substr(i);
}

bool node_id_less(const std::string& a, const std::string& b) {
    const bool a_num = is_unsigned_decimal(a);
    const bool b_num = is_unsigned_decimal(b);
    if (a_num && b_num) {
        const std::string_view av = strip_leading_zeros(a);
        const std::string_view bv = strip_leading_zeros(b);
        if (av.size() != bv.size()) {
            return av.size() < bv.size();
        }
        if (av != bv) {
            return av < bv;
        }
        return a < b;
    }
    if (a_num != b_num) {
        return a_num;
    }
    return a < b;
}

// Orient every bubble so `source` is the reference-left boundary and `sink` the reference-right one. The cactus
// pair is unordered, while every consumer reads it as an interval in reference order: `call` derives POS from
// the source, and merging joins one bubble's sink to the next one's source. A bubble the reference does not
// traverse has no order to take and is left alone.
void orient_bubbles_to_reference(const Graph& graph, const std::string& reference_path,
                                 std::vector<Bubble>& bubbles) {
    if (reference_path.empty() || bubbles.empty()) return;
    std::size_t ref_idx = graph.paths.size();
    for (std::size_t i = 0; i < graph.paths.size(); ++i)
        if (graph.paths[i].name == reference_path) { ref_idx = i; break; }
    if (ref_idx == graph.paths.size()) {
        std::size_t hits = 0, last = 0;
        for (std::size_t i = 0; i < graph.paths.size(); ++i)
            if (graph.paths[i].name.find(reference_path) != std::string::npos) { ++hits; last = i; }
        if (hits != 1) return;                      // ambiguous or absent: no order to impose
        ref_idx = last;
    }
    // First occurrence of each node on the reference walk. A boundary the reference visits more than
    // once has no single coordinate; the first is the one `call` already anchors on.
    std::unordered_map<std::string, std::size_t> first_pos;
    const PathRecord& ref = graph.paths[ref_idx];
    for (std::size_t i = 0; i < ref.steps.size(); ++i)
        first_pos.emplace(ref.steps[i].node_id, i);
    for (Bubble& b : bubbles) {
        const auto s_it = first_pos.find(b.source);
        const auto k_it = first_pos.find(b.sink);
        if (s_it == first_pos.end() || k_it == first_pos.end()) continue;
        // Each boundary is a handle whose orientation says which side of the node faces into the bubble.
        // Those sides are the site's own; the reference only sets the order. Read from the other end, the
        // pair swaps and each handle flips, which names the same two sides. Taking the reference's own
        // reading would mix two snarls on one node pair (on other sides) into one row.
        if (k_it->second < s_it->second) {
            std::swap(b.source, b.sink);
            const bool s_rev = b.source_reverse;
            b.source_reverse = !b.sink_reverse;
            b.sink_reverse = !s_rev;
        }
    }
}

bool bubble_endpoint_less(const Bubble& a, const Bubble& b) {
    if (a.source != b.source) {
        return node_id_less(a.source, b.source);
    }
    if (a.sink != b.sink) {
        return node_id_less(a.sink, b.sink);
    }
    return a.inside.size() < b.inside.size();
}




// Reference coordinates of a bubble's boundaries: the bp offset at which each boundary node starts on
// the reference walk, and its length. Bubbles the reference does not span cannot be ordered against
// each other and are never merged.
struct RefSpan {
    std::size_t src_start = 0;   // bp offset of the source node's first base
    std::size_t snk_start = 0;   // bp offset of the sink node's first base
    std::size_t snk_len = 0;
    std::size_t src_len = 0;
    bool ok = false;
};

std::unordered_map<std::string, std::pair<std::size_t, std::size_t>> reference_offsets(
    const Graph& graph, const std::string& reference_path) {
    std::unordered_map<std::string, std::pair<std::size_t, std::size_t>> out;   // node -> (bp start, len)
    if (reference_path.empty()) return out;
    std::size_t ref_idx = graph.paths.size();
    for (std::size_t i = 0; i < graph.paths.size(); ++i)
        if (graph.paths[i].name == reference_path) { ref_idx = i; break; }
    if (ref_idx == graph.paths.size()) {
        std::size_t hits = 0, last = 0;
        for (std::size_t i = 0; i < graph.paths.size(); ++i)
            if (graph.paths[i].name.find(reference_path) != std::string::npos) { ++hits; last = i; }
        if (hits != 1) return out;
        ref_idx = last;
    }
    std::size_t off = 0;
    for (const PathStep& st : graph.paths[ref_idx].steps) {
        const auto nit = graph.nodes.find(st.node_id);
        const std::size_t len = nit == graph.nodes.end() ? 0 : nit->second.sequence.size();
        out.emplace(st.node_id, std::make_pair(off, len));   // first occurrence wins
        off += len;
    }
    return out;
}

// The reference path spelled forward, for repeat detection. Empty when the path is not found.
std::string reference_sequence(const Graph& graph, const std::string& reference_path) {
    std::string out;
    if (reference_path.empty()) return out;
    const PathRecord* ref = nullptr;
    for (const PathRecord& p : graph.paths)
        if (p.name == reference_path) { ref = &p; break; }
    if (ref == nullptr) {
        std::size_t hits = 0;
        for (const PathRecord& p : graph.paths)
            if (p.name.find(reference_path) != std::string::npos) { ++hits; ref = &p; }
        if (hits != 1) return std::string();
    }
    for (const PathStep& st : ref->steps) {
        const auto nit = graph.nodes.find(st.node_id);
        if (nit == graph.nodes.end()) continue;
        const std::string& s = nit->second.sequence;
        if (!st.reverse) { out += s; continue; }
        for (auto it = s.rbegin(); it != s.rend(); ++it) {
            switch (*it) {
                case 'A': case 'a': out += 'T'; break;
                case 'C': case 'c': out += 'G'; break;
                case 'G': case 'g': out += 'C'; break;
                case 'T': case 't': out += 'A'; break;
                default: out += 'N';
            }
        }
    }
    return out;
}

// Merge bubbles that sit close together on the reference. Order and distance come from reference
// coordinates, not node ids: cactus order does not follow the reference. The gap is measured between the
// two interiors (see below). A fused bubble is validated: every interior node must lie inside the new span,
// or the merge is refused.
std::vector<Bubble> merge_nearby_bubbles(
    const Graph& graph,
    const PackedGraph& packed,
    const std::vector<PathIndex>& path_indexes,
    const std::vector<Bubble>& input,
    std::size_t max_bp,
    const std::string& reference_path,
    std::vector<std::vector<Bubble>>* parts_out,
    const std::unordered_map<std::string, std::pair<std::size_t, std::size_t>>* offsets = nullptr) {

    // `parts_out`, when given, receives for every returned bubble the input bubbles it was fused from (one
    // if not fused), so a fusion that later fails a filter can be undone without losing its valid parts.
    auto unfused = [&](const std::vector<Bubble>& bs) {
        if (parts_out != nullptr) {
            parts_out->clear();
            for (const Bubble& b : bs) parts_out->push_back({b});
        }
        return bs;
    };
    if (max_bp == 0 || input.size() < 2) {
        return unfused(input);
    }
    // `offsets`, when given, are the reference's (reference_offsets): callers merging many groups compute
    // them once.
    const auto computed = offsets == nullptr ? reference_offsets(graph, reference_path)
                                             : std::unordered_map<std::string, std::pair<std::size_t, std::size_t>>{};
    const auto& ref_off = offsets == nullptr ? computed : *offsets;
    if (ref_off.empty()) return unfused(input);   // no reference to order against: merging would be guesswork

    auto span_of = [&](const Bubble& b) {
        RefSpan r;
        const auto s = ref_off.find(b.source);
        const auto k = ref_off.find(b.sink);
        if (s == ref_off.end() || k == ref_off.end()) return r;
        r.src_start = s->second.first;  r.src_len = s->second.second;
        r.snk_start = k->second.first;  r.snk_len = k->second.second;
        r.ok = r.src_start <= r.snk_start;      // boundaries are reference-ordered by now
        return r;
    };

    // Reference order, and only bubbles the reference actually spans take part.
    std::vector<Bubble> spanned, unspanned;
    for (const Bubble& b : input) (span_of(b).ok ? spanned : unspanned).push_back(b);
    if (std::getenv("PANVAR_MERGE_DEBUG"))
        std::cerr << "[merge] spanned=" << spanned.size() << " unspanned=" << unspanned.size() << '\n';
    std::sort(spanned.begin(), spanned.end(), [&](const Bubble& a, const Bubble& b) {
        const RefSpan ra = span_of(a), rb = span_of(b);
        if (ra.src_start != rb.src_start) return ra.src_start < rb.src_start;
        return ra.snk_start < rb.snk_start;
    });

    std::vector<Bubble> merged;
    std::vector<std::vector<Bubble>> parts;   // parallel to `merged`
    merged.reserve(spanned.size());
    for (std::size_t i = 0; i < spanned.size(); ++i) {
        if (merged.empty()) { merged.push_back(spanned[i]); parts.push_back({spanned[i]}); continue; }
        Bubble& cur = merged.back();
        const RefSpan rc = span_of(cur), rn = span_of(spanned[i]);
        // The reference sequence between the two interiors: from the left bubble's sink to the end of
        // the right bubble's source, both boundaries included, so a shared boundary counts its own
        // length. Measured strictly between the boundaries, a long boundary node next to a short gap
        // would read as nearby.
        const std::size_t right_end = rn.src_start + rn.src_len;   // one past the right source's last base
        const std::size_t gap = rn.src_start >= rc.snk_start ? right_end - rc.snk_start : 0;   // overlap: 0
        if (std::getenv("PANVAR_MERGE_DEBUG"))
            std::cerr << "[merge] " << cur.source << ".." << cur.sink << " -> " << spanned[i].source
                      << ".." << spanned[i].sink << " gap=" << gap << " max=" << max_bp << '\n';
        if (gap > max_bp) { merged.push_back(spanned[i]); parts.push_back({spanned[i]}); continue; }

        Bubble fused = cur;
        fused.sink = spanned[i].sink;
        fused.sink_reverse = spanned[i].sink_reverse;
        // The fused interior is recomputed from the graph between the new outer boundaries, not assembled
        // from the parts: the parts plus the reference's route between them miss an alternate connector
        // branch.
        std::unordered_set<std::string> inside_set(cur.inside.begin(), cur.inside.end());
        inside_set.insert(spanned[i].inside.begin(), spanned[i].inside.end());
        inside_set.insert(cur.sink);                 // old boundaries become interior
        inside_set.insert(spanned[i].source);
        {
            const auto s_it = packed.node_idx_of.find(fused.source);
            const auto k_it = packed.node_idx_of.find(fused.sink);
            if (s_it != packed.node_idx_of.end() && k_it != packed.node_idx_of.end()) {
                for (const std::uint32_t idx : interior_indices_for_pair(
                         graph, packed, path_indexes, fused.source, fused.sink,
                         s_it->second, k_it->second)) {
                    inside_set.insert(packed.node_ids[idx]);
                }
            }
        }

        const RefSpan rf = span_of(fused);
        bool contained = rf.ok;
        if (contained) {
            const std::size_t lo = rf.src_start, hi = rf.snk_start + rf.snk_len;
            for (const std::string& n : inside_set) {
                if (n == fused.source || n == fused.sink) continue;
                const auto it = ref_off.find(n);
                if (it == ref_off.end()) continue;   // off-reference alleles legitimately sit inside
                if (it->second.first < lo || it->second.first >= hi) { contained = false; break; }
            }
        }
        if (!contained) {
            if (std::getenv("PANVAR_MERGE_DEBUG"))
                std::cerr << "[merge]   refused: fused span would not contain every interior node\n";
            merged.push_back(spanned[i]); parts.push_back({spanned[i]}); continue;
        }

        fused.inside.clear();
        fused.inside.reserve(inside_set.size());
        for (const std::string& n : inside_set)
            if (n != fused.source && n != fused.sink) fused.inside.push_back(n);
        std::sort(fused.inside.begin(), fused.inside.end());
        fused.path_support = 0; fused.min_inside_bp = 0; fused.max_inside_bp = 0;
        fused.long_path_support = 0; fused.inversion_signal = false;
        fused.distinct_alleles = 0; fused.ref_allele_support = 0;
        fused.alt_allele_support_max = 0; fused.alt_allele_support_min = 0;
        cur = std::move(fused);
        parts.back().push_back(spanned[i]);
    }
    for (Bubble& b : unspanned) { parts.push_back({b}); merged.push_back(std::move(b)); }
    std::vector<std::size_t> order(merged.size());
    for (std::size_t i = 0; i < order.size(); ++i) order[i] = i;
    std::sort(order.begin(), order.end(),
              [&](std::size_t a, std::size_t b) { return bubble_endpoint_less(merged[a], merged[b]); });
    std::vector<Bubble> sorted;
    sorted.reserve(merged.size());
    if (parts_out != nullptr) parts_out->clear();
    for (const std::size_t i : order) {
        sorted.push_back(std::move(merged[i]));
        if (parts_out != nullptr) parts_out->push_back(std::move(parts[i]));
    }
    return sorted;
}

// Would fusing `fs` out of the retained sites in `member_keys` lose a site in the overlap pass? That pass
// keeps, of two retained sites sharing an interior node, the one with fewer interior nodes, so a retained site
// outside the fusion that shares nodes with it would either drop it, with every site fused into it, or be
// dropped by it. One the pass would drop anyway, by a smaller retained site it overlaps, does not count.
// Compared by nodes, not reference coordinates: a nested site off the reference has no reference span.
bool fusion_loses_site(const Bubble& fs, const std::unordered_set<std::string>& member_keys,
                       const std::vector<Bubble>& kept) {
    auto key = [](const Bubble& b) { return b.source + "\t" + b.sink; };
    const std::unordered_set<std::string> fused_nodes(fs.inside.begin(), fs.inside.end());
    auto ranks_before = [](const Bubble& a, const Bubble& b) {
        return a.inside.size() != b.inside.size() ? a.inside.size() < b.inside.size() : bubble_endpoint_less(a, b);
    };
    auto shares = [](const Bubble& a, const std::unordered_set<std::string>& nodes) {
        for (const std::string& n : a.inside) if (nodes.count(n)) return true;
        return false;
    };
    for (const Bubble& b : kept) {
        if (member_keys.count(key(b)) || !shares(b, fused_nodes)) continue;
        if (ranks_before(b, fs)) return true;
        const std::unordered_set<std::string> b_nodes(b.inside.begin(), b.inside.end());
        bool dropped_anyway = false;
        for (const Bubble& c : kept)
            if (key(c) != key(b) && ranks_before(c, b) && shares(c, b_nodes)) { dropped_anyway = true; break; }
        if (!dropped_anyway) return true;
    }
    return false;
}

// A site's distinct walks, kept by the first scoring pass so fusing by placement need not walk every
// path through every candidate site again.
struct SiteWalks {
    std::optional<std::vector<PathStep>> ref;
    std::vector<std::vector<PathStep>> alts;
};

// The crossings in a walk across a site (bubble_steps): [from, to] step ranges, each from a boundary to a
// boundary, split wherever the walk leaves the site's nodes. Passes that follow each other directly stay one
// range: a tandem copy.
std::vector<std::pair<std::size_t, std::size_t>> site_passes(const std::vector<PathStep>& steps, const Bubble& b,
                                                             const std::unordered_set<std::string>& inside) {
    std::vector<std::pair<std::size_t, std::size_t>> out;
    const auto own = [&](std::size_t i) {
        const std::string& n = steps[i].node_id;
        return n == b.source || n == b.sink || inside.count(n) > 0;
    };
    std::size_t from = 0;
    for (std::size_t i = 1; i < steps.size(); ++i) {
        if (own(i)) continue;
        if (i - 1 > from) out.emplace_back(from, i - 1);
        while (i < steps.size() && !own(i)) ++i;
        from = i;
    }
    if (from + 1 < steps.size()) out.emplace_back(from, steps.size() - 1);
    if (out.empty() && steps.size() >= 2) out.emplace_back(0, steps.size() - 1);
    return out;
}

// Fuse sites whose variants could be placed in each other. In a repeat the graph can split one indel over
// several sites, the largest kept and the rest under the size floor; which sites belong together is read from
// the variants, not from a definition of a repeat. Each haplotype's allele is trimmed of what it shares with
// the reference allele at either end. An insertion or deletion left over is slid left and right as VCF
// normalisation does, across up to three differing bases (placement_range), possibly out of the site; a
// replacement stays where it is. A site's own stretch is where its changes sit, its reach where they can be
// placed; two sites are one locus when the reach of one meets the other's own stretch. Each group holding a
// retained site is fused as --merge-nearby-bp fuses; other sites are left as they were.
std::vector<Bubble> fuse_ambiguous_placements(
    const Graph& graph,
    const PackedGraph& packed,
    const std::vector<PathIndex>& path_indexes,
    const std::vector<Bubble>& kept,
    const std::vector<Bubble>& candidates,
    const std::string& reference_path,
    const std::unordered_map<std::string, SiteWalks>& walk_cache,
    std::vector<std::vector<Bubble>>* parts_out,
    std::size_t* fused_out) {

    if (fused_out != nullptr) *fused_out = 0;
    if (parts_out != nullptr) parts_out->clear();
    const auto ref_off = reference_offsets(graph, reference_path);
    std::string ref = reference_sequence(graph, reference_path);
    auto unchanged = [&]() {
        if (parts_out != nullptr) parts_out->assign(kept.size(), {});
        return kept;
    };
    if (ref_off.empty() || ref.empty()) return unchanged();
    for (char& c : ref) c = static_cast<char>(c & 0xDF);
    std::size_t ref_idx = graph.paths.size();
    for (std::size_t i = 0; i < graph.paths.size(); ++i)
        if (graph.paths[i].name == reference_path) { ref_idx = i; break; }
    if (ref_idx == graph.paths.size()) return unchanged();

    auto comp = [](char c) {
        switch (c) { case 'A': return 'T'; case 'C': return 'G'; case 'G': return 'C'; case 'T': return 'A'; default: return 'N'; }
    };
    // The allele between a walk's first and last step; with `own`, only its steps on those nodes, so what a path
    // walks outside the site between two passes is not part of it.
    auto spell = [&](const std::vector<PathStep>& steps, const std::unordered_set<std::string>* own = nullptr) {
        std::string out;
        for (std::size_t i = 1; i + 1 < steps.size(); ++i) {
            if (own != nullptr && !own->count(steps[i].node_id)) continue;
            const auto it = graph.nodes.find(steps[i].node_id);
            if (it == graph.nodes.end()) continue;
            const std::string& sq = it->second.sequence;
            if (!steps[i].reverse) for (const char c : sq) out += static_cast<char>(c & 0xDF);
            else for (std::size_t k = sq.size(); k-- > 0;) out += comp(static_cast<char>(sq[k] & 0xDF));
        }
        return out;
    };
    auto revcomp = [&](const std::string& x) {
        std::string r(x.rbegin(), x.rend());
        for (char& c : r) c = comp(c);
        return r;
    };
    // [lo, hi) of a site's interior on the reference, from its boundaries' first occurrences.
    auto interior = [&](const Bubble& b, std::size_t& lo, std::size_t& hi) {
        const auto s = ref_off.find(b.source), k = ref_off.find(b.sink);
        if (s == ref_off.end() || k == ref_off.end() || s->second.first > k->second.first) return false;
        lo = s->second.first + s->second.second;
        hi = k->second.first;
        return lo <= hi;
    };
    constexpr std::size_t kMaxInterior = 100000;   // larger sites keep their own span
    // The stretches of reference points (a point is the position between two bases) the site's variants
    // can occupy: `reach` within the placement tolerance, `own` exactly (the normalisation range, or a
    // replacement's span). False when the site carries no change or cannot be placed. A site whose changes
    // all keep the length (substitutions) links nothing: two indels that could each be placed onto one such
    // site would otherwise be joined, each having moved.
    struct Stretch { std::size_t reach_first = SIZE_MAX, reach_last = 0, own_first = SIZE_MAX, own_last = 0; };
    auto stretch = [&](const Bubble& b, Stretch& st_out) {
        std::size_t lo = 0, hi = 0;
        if (!interior(b, lo, hi)) return false;
        // Too large to place base by base: its own span, linking only if its alleles differ in length.
        if (hi - lo > kMaxInterior) { st_out = {lo, hi, lo, hi}; return b.max_inside_bp != b.min_inside_bp; }
        const std::string r = ref.substr(lo, hi - lo);
        SiteWalks fresh;
        const SiteWalks* sw = nullptr;
        if (const auto it = walk_cache.find(b.source + "\t" + b.sink); it != walk_cache.end()) {
            sw = &it->second;
        } else {   // a site made after scoring (a fusion): walk the paths
            fresh.ref = bubble_steps(graph.paths[ref_idx], path_indexes[ref_idx], b, nullptr, nullptr, true);
            std::unordered_set<std::string> walked;
            for (std::size_t p = 0; p < graph.paths.size(); ++p) {
                if (p == ref_idx) continue;
                auto st = bubble_steps(graph.paths[p], path_indexes[p], b, nullptr, nullptr, true);
                if (!st.has_value() || st->size() < 2) continue;
                std::string sig;
                for (const PathStep& x : *st) { sig += x.node_id; sig += x.reverse ? '-' : '+'; }
                if (walked.insert(std::move(sig)).second) fresh.alts.push_back(std::move(*st));
            }
            sw = &fresh;
        }
        if (!sw->ref.has_value()) return false;
        std::unordered_set<std::string> own(b.inside.begin(), b.inside.end());
        own.insert(b.source);
        own.insert(b.sink);
        const std::string rw = spell(*sw->ref, &own);
        bool rev = false;
        if (rw != r) { if (revcomp(rw) == r) rev = true; else return false; }
        std::unordered_set<std::string> seen;
        bool any = false;
        st_out = Stretch{};
        std::vector<std::vector<PathStep>> crossings;   // each crossing is an allele (see site_passes)
        const std::unordered_set<std::string> inside_set(b.inside.begin(), b.inside.end());
        for (const std::vector<PathStep>& st : sw->alts)
            for (const auto& [from, to] : site_passes(st, b, inside_set))
                crossings.emplace_back(st.begin() + static_cast<long>(from), st.begin() + static_cast<long>(to) + 1);
        for (const std::vector<PathStep>& st : crossings) {
            if (st.size() < 2) continue;
            std::string a = spell(st, &own);
            if (rev) a = revcomp(a);
            if (a == r || !seen.insert(a).second) continue;
            std::size_t pre = 0, suf = 0;
            while (pre < r.size() && pre < a.size() && r[pre] == a[pre]) ++pre;
            while (suf < r.size() - pre && suf < a.size() - pre && r[r.size() - 1 - suf] == a[a.size() - 1 - suf]) ++suf;
            const std::size_t dl = r.size() - pre - suf, il = a.size() - pre - suf;
            const std::size_t s0 = lo + pre;                  // deleted [s0, s0 + dl), or insertion point s0
            const std::string ins = a.substr(pre, il);
            const PlacementRange wide = placement_range(ref, s0, s0 + dl, ins);
            const PlacementRange exact = placement_range(ref, s0, s0 + dl, ins, 0);
            st_out.reach_first = std::min(st_out.reach_first, wide.first);
            st_out.reach_last = std::max(st_out.reach_last, wide.last);
            st_out.own_first = std::min(st_out.own_first, exact.first);
            st_out.own_last = std::max(st_out.own_last, exact.last);
            any = any || dl != il;
        }
        return any;
    };

    auto key = [](const Bubble& b) { return b.source + "\t" + b.sink; };
    struct Item { const Bubble* b; bool kept; Stretch st; };
    std::vector<Item> items;
    std::unordered_set<std::string> kept_keys;
    std::vector<std::pair<std::size_t, std::size_t>> kept_spans;
    for (const Bubble& b : kept) {
        kept_keys.insert(key(b));
        std::size_t lo = 0, hi = 0;
        if (interior(b, lo, hi)) kept_spans.emplace_back(lo, hi);
    }
    // A site enclosing a smaller retained site is dropped by the overlap pass whatever is fused, so it
    // takes no part here: as a member its stretch would be its whole span, and a group holding a site
    // and the sites inside it cannot be fused into one.
    // Spans sorted by start, so both tests below look only at the spans that can answer them.
    std::sort(kept_spans.begin(), kept_spans.end());
    auto encloses_kept = [&](std::size_t lo, std::size_t hi) {
        for (auto it = std::lower_bound(kept_spans.begin(), kept_spans.end(), std::make_pair(lo, std::size_t{0}));
             it != kept_spans.end() && it->first <= hi; ++it)
            if (it->second <= hi && (it->first != lo || it->second != hi)) return true;
        return false;
    };
    std::unordered_set<std::string> containers;
    std::vector<std::pair<std::size_t, std::size_t>> plain_spans;   // kept sites that enclose none
    for (const Bubble& b : kept) {
        std::size_t lo = 0, hi = 0;
        const bool placed = interior(b, lo, hi);
        if (placed && encloses_kept(lo, hi)) { containers.insert(key(b)); continue; }
        if (placed) plain_spans.emplace_back(lo, hi);
        Stretch st;
        if (stretch(b, st)) items.push_back({&b, true, st});
    }
    // Whether [lo, hi) lies inside one of them: the furthest end among the spans starting at or before lo.
    std::sort(plain_spans.begin(), plain_spans.end());
    std::vector<std::size_t> reach_to(plain_spans.size());
    for (std::size_t i = 0; i < plain_spans.size(); ++i)
        reach_to[i] = std::max(i ? reach_to[i - 1] : 0, plain_spans[i].second);
    auto inside_plain = [&](std::size_t lo, std::size_t hi) {
        const auto it = std::upper_bound(plain_spans.begin(), plain_spans.end(), std::make_pair(lo, SIZE_MAX));
        return it != plain_spans.begin() && reach_to[static_cast<std::size_t>(it - plain_spans.begin()) - 1] >= hi;
    };
    for (const Bubble& c : candidates) {
        if (kept_keys.count(key(c))) continue;
        std::size_t lo = 0, hi = 0;
        if (!interior(c, lo, hi) || encloses_kept(lo, hi)) continue;
        if (inside_plain(lo, hi)) continue;   // a part of a site --merge-nearby-bp already fused
        Stretch st;
        if (stretch(c, st)) items.push_back({&c, false, st});
    }
    std::sort(items.begin(), items.end(), [](const Item& x, const Item& y) {
        return x.st.reach_first != y.st.reach_first ? x.st.reach_first < y.st.reach_first
                                                    : x.st.reach_last < y.st.reach_last;
    });
    // Two sites are one locus when a change of one can be placed in the other's own stretch: the
    // tolerance is spent once, by the change that moves. Letting both move and meet in between would
    // allow twice the mismatches a single placement is allowed.
    std::vector<std::size_t> parent(items.size());
    std::iota(parent.begin(), parent.end(), std::size_t{0});
    std::function<std::size_t(std::size_t)> find = [&](std::size_t x) {
        while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; }
        return x;
    };
    auto meets = [](std::size_t f1, std::size_t l1, std::size_t f2, std::size_t l2) { return f1 <= l2 && f2 <= l1; };
    for (std::size_t i = 0; i < items.size(); ++i) {
        for (std::size_t j = i + 1; j < items.size() && items[j].st.reach_first <= items[i].st.reach_last; ++j) {
            const Stretch& a = items[i].st;
            const Stretch& b = items[j].st;
            if (meets(a.reach_first, a.reach_last, b.own_first, b.own_last) ||
                meets(b.reach_first, b.reach_last, a.own_first, a.own_last))
                parent[find(i)] = find(j);
        }
    }
    std::map<std::size_t, std::vector<std::size_t>> groups;   // root -> items, in stretch order
    for (std::size_t i = 0; i < items.size(); ++i) groups[find(i)].push_back(i);

    std::vector<Bubble> out;
    std::vector<std::vector<Bubble>> parts;
    std::unordered_set<std::string> replaced;   // kept sites fused into a group
    for (const auto& [root, members] : groups) {
        std::vector<Bubble> group;
        std::vector<Bubble> group_kept;
        std::unordered_set<std::string> member_keys;
        std::size_t span_lo = SIZE_MAX, span_hi = 0;
        for (const std::size_t k : members) {
            group.push_back(*items[k].b);
            member_keys.insert(key(*items[k].b));
            if (items[k].kept) group_kept.push_back(*items[k].b);
            std::size_t lo = 0, hi = 0;
            if (interior(*items[k].b, lo, hi)) { span_lo = std::min(span_lo, lo); span_hi = std::max(span_hi, hi); }
        }
        if (group.size() < 2 || group_kept.empty()) continue;
        // Groups chain (a reaches b, b reaches c), and nothing else bounds how far: an imperfect repeat
        // long enough could chain across a locus. Past the span a single site is scanned for, the
        // group is left as it was.
        if (span_lo != SIZE_MAX && span_hi > span_lo && span_hi - span_lo > kMaxInterior) continue;
        const std::vector<Bubble> fused = merge_nearby_bubbles(
            graph, packed, path_indexes, group, std::numeric_limits<std::size_t>::max() / 4, reference_path, nullptr,
            &ref_off);
        if (fused.size() != 1) continue;   // a fusion was refused: leave the sites
        // A fusion the overlap pass would undo at the cost of a site is not made.
        if (fusion_loses_site(fused.front(), member_keys, kept)) continue;
        out.push_back(fused.front());
        parts.push_back(group_kept);
        for (const Bubble& b : group_kept) replaced.insert(key(b));
        if (fused_out != nullptr) ++*fused_out;
    }
    for (const Bubble& b : kept)
        if (!replaced.count(key(b))) { out.push_back(b); parts.emplace_back(); }
    if (parts_out != nullptr) *parts_out = std::move(parts);
    return out;
}

// ---- Path indexing + snarl JSONL import ----------------------------------------------

PathIndex build_path_index(const PathRecord& path) { return build_bubble_path_index(path); }

std::optional<std::string> extract_node_id_for_key(
    const std::string& line,
    const std::string& key,
    bool use_last_key) {

    const std::string token = "\"" + key + "\"";
    const std::size_t key_pos = use_last_key ? line.rfind(token) : line.find(token);
    if (key_pos == std::string::npos) {
        return std::nullopt;
    }
    const std::size_t node_key_pos = line.find("\"node_id\"", key_pos);
    if (node_key_pos == std::string::npos) {
        return std::nullopt;
    }
    const std::size_t colon_pos = line.find(':', node_key_pos);
    if (colon_pos == std::string::npos) {
        return std::nullopt;
    }
    // vg writes ids as strings; an unquoted number is read as the same id.
    std::size_t v = colon_pos + 1;
    while (v < line.size() && (line[v] == ' ' || line[v] == '\t')) ++v;
    if (v < line.size() && std::isdigit(static_cast<unsigned char>(line[v]))) {
        std::size_t e = v;
        while (e < line.size() && std::isdigit(static_cast<unsigned char>(line[e]))) ++e;
        return line.substr(v, e - v);
    }
    const std::size_t q1 = line.find('"', colon_pos + 1);
    if (q1 == std::string::npos) {
        return std::nullopt;
    }
    const std::size_t q2 = line.find('"', q1 + 1);
    if (q2 == std::string::npos || q2 <= q1 + 1) {
        return std::nullopt;
    }
    return line.substr(q1 + 1, q2 - q1 - 1);
}

// The `backward` flag of the handle under `key` ("start" or "end"): absent means false, as vg omits a
// false field.
bool extract_backward_for_key(const std::string& line, const std::string& key, bool use_last_key) {
    const std::string token = "\"" + key + "\"";
    const std::size_t key_pos = use_last_key ? line.rfind(token) : line.find(token);
    if (key_pos == std::string::npos) return false;
    const std::size_t close = line.find('}', key_pos);
    const std::size_t b = line.find("\"backward\"", key_pos);
    if (b == std::string::npos || (close != std::string::npos && b > close)) return false;
    const std::size_t colon = line.find(':', b);
    if (colon == std::string::npos) return false;
    std::size_t v = colon + 1;
    while (v < line.size() && (line[v] == ' ' || line[v] == '\t')) ++v;
    return line.compare(v, 4, "true") == 0;
}

// One crossing of a snarl endpoint pair along one path: the window between two consecutive
// boundary occurrences of opposite kind (source then sink, or sink then source).
struct EndpointOnlyInterval {
    std::size_t left = 0;
    std::size_t right = 0;
    std::size_t inside_steps = 0;
};

// One occurrence of a boundary node on a path, with the orientation it is walked in.
struct BoundaryHit {
    std::size_t pos = 0;
    bool is_source = false;
    bool reverse = false;
};

// Every occurrence of either boundary along `path`, in step order.
std::vector<BoundaryHit> boundary_hits(
    const PathIndex& index,
    const PathRecord& path,
    const std::string& source,
    const std::string& sink) {

    const auto src_it = index.positions.find(source);
    const auto sink_it = index.positions.find(sink);
    if (src_it == index.positions.end() || sink_it == index.positions.end()) {
        return {};
    }
    const std::vector<std::size_t>& src_pos = src_it->second;    // sorted: built in step order
    const std::vector<std::size_t>& sink_pos = sink_it->second;
    std::vector<BoundaryHit> out;
    out.reserve(src_pos.size() + sink_pos.size());
    std::size_t i = 0, j = 0;
    while (i < src_pos.size() || j < sink_pos.size()) {
        const bool take_source =
            j == sink_pos.size() || (i < src_pos.size() && src_pos[i] < sink_pos[j]);
        const std::size_t pos = take_source ? src_pos[i++] : sink_pos[j++];
        out.push_back(BoundaryHit{pos, take_source, path.steps[pos].reverse});
    }
    return out;
}

// Which side of each boundary the site's interior hangs off (true: the node's end side). A step departs
// through the end side iff it is forward, and arrives through the end side iff it is reversed.
struct InwardSides {
    bool source_end = true;
    bool sink_end = false;
};

bool departs_inward(const BoundaryHit& h, const InwardSides& in) {
    return !h.reverse == (h.is_source ? in.source_end : in.sink_end);
}

bool arrives_inward(const BoundaryHit& h, const InwardSides& in) {
    return h.reverse == (h.is_source ? in.source_end : in.sink_end);
}

// The panel's inward sides for one endpoint pair, voted from the walks. Each consecutive source/sink
// occurrence pair implies one assignment: a crossing the true one, the stretch between two crossings (out
// through one boundary's outer side, back in through the other's) its complement. Excursions exist only on
// paths that cross the site more than once, so across the panel they are outvoted, where one path's own vote
// could tie. A tie keeps the reference-sorted orientation: the source left forward, the sink entered forward.
// Nullopt when no path crosses the pair.
std::optional<InwardSides> panel_inward_sides(const std::vector<std::vector<BoundaryHit>>& hits) {
    std::size_t votes[4] = {0, 0, 0, 0};   // index = source_end * 2 + sink_end
    bool any = false;
    for (const std::vector<BoundaryHit>& h : hits) {
        for (std::size_t k = 0; k + 1 < h.size(); ++k) {
            const BoundaryHit& a = h[k];
            const BoundaryHit& b = h[k + 1];
            if (a.is_source == b.is_source) continue;
            const bool a_end = !a.reverse;   // the side `a` is left through
            const bool b_end = b.reverse;    // the side `b` is entered through
            const bool s_end = a.is_source ? a_end : b_end;
            const bool k_end = a.is_source ? b_end : a_end;
            ++votes[(s_end ? 2 : 0) + (k_end ? 1 : 0)];
            any = true;
        }
    }
    if (!any) return std::nullopt;
    // Preference order among equal votes: the reference-sorted orientation first.
    const int order[4] = {2, 1, 3, 0};
    int best = order[0];
    for (const int c : order)
        if (votes[c] > votes[best]) best = c;
    return InwardSides{(best & 2) != 0, (best & 1) != 0};
}

bool better_endpoint_interval(const EndpointOnlyInterval& a, const EndpointOnlyInterval& b) {
    if (a.inside_steps != b.inside_steps) {
        return a.inside_steps > b.inside_steps;
    }
    const std::size_t a_span = a.right - a.left;
    const std::size_t b_span = b.right - b.left;
    if (a_span != b_span) {
        return a_span < b_span;
    }
    return a.left < b.left;
}

// The stretches of one path that lie inside the site, one window per stretch. A stretch, between two
// consecutive boundary occurrences, is inside when the path leaves the first through its inward side and
// arrives at the second through its inward side: a crossing in either direction (`S+ .. K+`, `K- .. S-`) or
// a hairpin walked inside (`S+ .. S-`), but not the excursion between two crossings (`K+ .. S+`). Judged per
// stretch, so a recurring boundary does not sweep excursions in, and a path starting or ending between
// crossings cannot tip the decision. An adjacent pair (sink at source + 1) is a crossing with an empty
// interior and is kept: when every path takes the direct route, it is what tells the graph search which side
// of each boundary faces inward.
std::vector<EndpointOnlyInterval> inside_stretches(
    const std::vector<BoundaryHit>& hits,
    const InwardSides& inward) {

    std::vector<EndpointOnlyInterval> out;
    for (std::size_t k = 0; k + 1 < hits.size(); ++k) {
        const BoundaryHit& a = hits[k];
        const BoundaryHit& b = hits[k + 1];
        if (!departs_inward(a, inward) || !arrives_inward(b, inward)) continue;
        EndpointOnlyInterval c;
        c.left = a.pos;
        c.right = b.pos;
        c.inside_steps = b.pos - a.pos - 1;
        out.push_back(c);
    }
    return out;
}

// Parse top-level (source, sink) boundary pairs from a vg snarls JSONL, skipping nested snarls. `inward`
// receives each pair's inward sides from the handles' `backward` flags, since the same two nodes can bound
// different sites on different sides. Only a file carrying the flag somewhere is read this way (vg writes it
// only where true); otherwise `inward` stays empty and the sides are voted from the walks.
std::vector<std::pair<std::string, std::string>> read_snarl_pairs_jsonl(const std::string& path,
                                                                        std::vector<SnarlInward>* inward) {
    std::ifstream snarls_in(path);
    if (!snarls_in) {
        throw std::runtime_error("Failed to open snarl JSONL input: " + path);
    }
    std::vector<std::pair<std::string, std::string>> pairs;
    bool any_backward = false;
    std::string line;
    while (std::getline(snarls_in, line)) {
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        if (line.find("\"start\"") == std::string::npos || line.find("\"end\"") == std::string::npos) {
            continue;
        }
        if (line.find("\"parent\"") != std::string::npos) {
            continue; // top-level snarls only
        }
        const auto source_opt = extract_node_id_for_key(line, "start", true);
        const auto sink_opt = extract_node_id_for_key(line, "end", false);
        if (!source_opt.has_value() || !sink_opt.has_value()) {
            continue;
        }
        pairs.emplace_back(*source_opt, *sink_opt);
        if (inward != nullptr)
            inward->push_back(SnarlInward{!extract_backward_for_key(line, "start", true),
                                          extract_backward_for_key(line, "end", false)});
        any_backward = any_backward || line.find("\"backward\"") != std::string::npos;
    }
    if (inward != nullptr && !any_backward) inward->clear();
    return pairs;
}

// The interior of one endpoint pair: everything the panel walks between the boundaries, unioned with
// everything the graph carries between them. Candidates from either snarl source and fused bubbles all take
// their interior from here, so they follow one rule.
std::unordered_set<std::uint32_t> interior_indices_for_pair(
    const Graph& graph,
    const PackedGraph& packed,
    const std::vector<PathIndex>& path_indexes,
    const std::string& source_id,
    const std::string& sink_id,
    std::uint32_t src_idx,
    std::uint32_t sink_idx,
    const SnarlInward* fixed,
    SnarlInward* used) {

    std::unordered_set<std::uint32_t> inside_idx_set;
    std::optional<EndpointOnlyInterval> best_iv;
    std::size_t best_iv_path = 0;
    std::vector<std::vector<BoundaryHit>> hits(path_indexes.size());
    for (std::size_t p_idx = 0; p_idx < path_indexes.size(); ++p_idx) {
        hits[p_idx] = boundary_hits(path_indexes[p_idx], graph.paths[p_idx], source_id, sink_id);
    }
    // The decomposition's own orientation when it gave one: voted from the walks, a snarl whose every path
    // walks its outside (one closing a root cycle around the locus ends) would take the outside as interior.
    const std::optional<InwardSides> inward =
        fixed != nullptr ? std::optional<InwardSides>(InwardSides{fixed->source_end, fixed->sink_end})
                         : panel_inward_sides(hits);
    if (!inward.has_value()) {
        return inside_idx_set;
    }
    if (used != nullptr) *used = SnarlInward{inward->source_end, inward->sink_end};
    for (std::size_t p_idx = 0; p_idx < path_indexes.size(); ++p_idx) {
        // Every inside stretch contributes its interior (see inside_stretches). Only a stretch that runs
        // boundary to other boundary seeds the graph search below: a hairpin returns to the boundary it
        // left, so its two ends do not orient the pair.
        const std::vector<EndpointOnlyInterval> crossings = inside_stretches(hits[p_idx], *inward);
        const auto& steps = graph.paths[p_idx].steps;
        for (const EndpointOnlyInterval& crossing : crossings) {
            const bool spans_pair = steps[crossing.left].node_id != steps[crossing.right].node_id;
            if (spans_pair && (!best_iv.has_value() || better_endpoint_interval(crossing, *best_iv))) {
                best_iv = crossing;
                best_iv_path = p_idx;
            }
            for (std::size_t i = crossing.left + 1; i < crossing.right; ++i) {
                const auto idx_it = packed.node_idx_of.find(steps[i].node_id);
                if (idx_it == packed.node_idx_of.end()) {
                    continue;
                }
                if (idx_it->second == src_idx || idx_it->second == sink_idx) {
                    continue;
                }
                inside_idx_set.insert(idx_it->second);
            }
        }
    }

    // Add the branches the graph carries between these boundaries that no stored path takes. The walked
    // interval supplies which side of each boundary faces inward, which the unordered cactus pair does not.
    if (best_iv.has_value()) {
        const auto& steps = graph.paths[best_iv_path].steps;
        const PathStep& near = steps[best_iv->left];
        const PathStep& far = steps[best_iv->right];
        const auto near_it = packed.node_idx_of.find(near.node_id);
        const auto far_it = packed.node_idx_of.find(far.node_id);
        if (near_it != packed.node_idx_of.end() && far_it != packed.node_idx_of.end()) {
            const std::uint64_t near_handle =
                (static_cast<std::uint64_t>(near_it->second) << 1) | (near.reverse ? 1u : 0u);
            const std::uint64_t far_handle =
                (static_cast<std::uint64_t>(far_it->second) << 1) | (far.reverse ? 1u : 0u);
            std::unordered_set<std::uint32_t> forward;
            std::unordered_set<std::uint32_t> backward;
            if (reachable_interior(graph, packed, near_handle, far_it->second, near_it->second,
                                   false, forward) &&
                reachable_interior(graph, packed, far_handle, near_it->second, far_it->second,
                                   true, backward)) {
                for (const std::uint32_t idx : forward) {
                    if (backward.count(idx) == 0) continue;   // not on a boundary-to-boundary route
                    if (idx == src_idx || idx == sink_idx) continue;
                    inside_idx_set.insert(idx);
                }
            } else {
                g_interior_traversal_truncated.fetch_add(1, std::memory_order_relaxed);
            }
        }
    }
    return inside_idx_set;
}

// The interior of a snarl is a property of the graph, not of the panel walking it: a branch no stored path
// takes is still part of the variant, and the span filters and merging's span check read `inside`. So the
// interior also takes, in the standard snarl sense, the nodes reachable from the near boundary's inner handle
// that reach the far boundary without passing a boundary, over oriented handles. It is unioned with the path-derived set, so a pair that is not really a
// snarl keeps the nodes the panel walks between its boundaries; and the search is capped, since a leaky pair
// floods the component, in which case the path-derived set stands alone.
constexpr std::size_t kMaxInteriorHandles = 1u << 20;

// Successors of an encoded handle (node index * 2 + orientation; 0 forward, 1 reverse).
// `backward` walks predecessors instead, using the bidirected identity
// pred(n, o) == { reverse(h) : h in succ(n, 1 - o) }, which needs no reverse index.
void handle_successors(
    const Graph& graph,
    const PackedGraph& packed,
    std::uint64_t handle,
    bool backward,
    std::vector<std::uint64_t>& out) {

    out.clear();
    const std::uint32_t idx = static_cast<std::uint32_t>(handle >> 1);
    const int orient = static_cast<int>(handle & 1u);
    const int depart = backward ? (1 - orient) : orient;

    const auto nit = graph.nodes.find(packed.node_ids[idx]);
    if (nit == graph.nodes.end()) return;
    const std::vector<Neighbor>& arcs = depart ? nit->second.start : nit->second.end;
    for (const Neighbor& nb : arcs) {
        const auto jt = packed.node_idx_of.find(nb.node_id);
        if (jt == packed.node_idx_of.end()) continue;
        // side 0 = entering the neighbour's start = arriving on its forward handle
        std::uint64_t next = (static_cast<std::uint64_t>(jt->second) << 1) | (nb.side == 0 ? 0u : 1u);
        if (backward) next ^= 1u;   // the predecessor is that handle reversed
        out.push_back(next);
    }
}

// Nodes reachable from `start_handle`, stopping at `stop_idx` and never passing back through
// `origin_idx`. Returns false if the cap is hit, in which case `out` is meaningless.
bool reachable_interior(
    const Graph& graph,
    const PackedGraph& packed,
    std::uint64_t start_handle,
    std::uint32_t stop_idx,
    std::uint32_t origin_idx,
    bool backward,
    std::unordered_set<std::uint32_t>& out) {

    out.clear();
    std::unordered_set<std::uint64_t> seen;
    std::vector<std::uint64_t> stack{start_handle};
    std::vector<std::uint64_t> succ;
    seen.insert(start_handle);

    while (!stack.empty()) {
        const std::uint64_t h = stack.back();
        stack.pop_back();
        handle_successors(graph, packed, h, backward, succ);
        for (const std::uint64_t nh : succ) {
            const std::uint32_t nidx = static_cast<std::uint32_t>(nh >> 1);
            if (nidx == origin_idx) continue;   // would leave through the near boundary
            if (nidx == stop_idx) continue;     // reached the far boundary; do not expand past it
            out.insert(nidx);
            if (seen.size() >= kMaxInteriorHandles) return false;
            if (seen.insert(nh).second) stack.push_back(nh);
        }
    }
    return true;
}

// Whether a path walks a node, at these steps, twice in the same direction.
bool repeats_one_way(const PathRecord& path, const std::vector<std::size_t>& steps) {
    if (steps.size() < 2) return false;
    std::size_t fwd = 0;
    for (const std::size_t k : steps) fwd += !path.steps[k].reverse;
    return fwd >= 2 || steps.size() - fwd >= 2;
}

std::size_t collect_candidates_for_pairs(
    const Graph& graph,
    const PackedGraph& packed,
    const std::vector<std::pair<std::string, std::string>>& pairs,
    const std::vector<SnarlInward>& inward,
    const std::vector<PathIndex>& path_indexes,
    CandidateMap& dedup,
    std::vector<SnarlDebugEntry>* debug_entries,
    const PathRecord* ref,
    std::vector<std::pair<std::string, std::string>>* loop_pairs) {

    std::size_t candidate_id = 0;
    // Whether some path walks a node twice in one direction: a snarl between two such nodes can be the inside
    // of a loop, a tandem copy the graph stores once. An inverted repeat is walked once each way.
    const auto walked_twice = [&](const std::string& id) {
        for (std::size_t p = 0; p < path_indexes.size(); ++p) {
            const auto it = path_indexes[p].positions.find(id);
            if (it != path_indexes[p].positions.end() && repeats_one_way(graph.paths[p], it->second)) return true;
        }
        return false;
    };
    std::size_t collapsed = 0;
    // Whether the reference walks the candidate from one boundary to the other through its inward sides,
    // in either direction, with nothing but its interior between them.
    const auto ref_crosses = [&](const BubbleCandidateIdx& c) {
        if (ref == nullptr) return false;
        const std::string& src = packed.node_ids[c.source];
        const std::string& snk = packed.node_ids[c.sink];
        std::unordered_set<std::string> in;
        for (const std::uint32_t idx : c.inside) in.insert(packed.node_ids[idx]);
        const std::vector<PathStep>& st = ref->steps;
        const auto from = [&](std::size_t i, const std::string& to, bool to_rev) {
            std::size_t j = i + 1;
            while (j < st.size() && in.count(st[j].node_id)) ++j;
            return j < st.size() && st[j].node_id == to && st[j].reverse == to_rev;
        };
        for (std::size_t i = 0; i < st.size(); ++i) {
            if (st[i].node_id == src && st[i].reverse == !c.source_end && from(i, snk, c.sink_end)) return true;
            if (st[i].node_id == snk && st[i].reverse == !c.sink_end && from(i, src, c.source_end)) return true;
        }
        return false;
    };
    for (const auto& [source_id, sink_id] : pairs) {
        ++candidate_id;
        SnarlDebugEntry debug;
        debug.candidate_id = candidate_id;

        if (source_id == sink_id) {
            if (debug_entries != nullptr) debug_entries->push_back(std::move(debug));
            continue;
        }

        const auto src_idx_it = packed.node_idx_of.find(source_id);
        const auto sink_idx_it = packed.node_idx_of.find(sink_id);
        if (src_idx_it == packed.node_idx_of.end() || sink_idx_it == packed.node_idx_of.end()) {
            if (debug_entries != nullptr) {
                debug.source = source_id;
                debug.sink = sink_id;
                debug_entries->push_back(std::move(debug)); // endpoint not in graph
            }
            continue;
        }

        const SnarlInward* fixed = inward.size() == pairs.size() ? &inward[candidate_id - 1] : nullptr;
        SnarlInward used;
        std::unordered_set<std::uint32_t> inside_idx_set = interior_indices_for_pair(
            graph, packed, path_indexes, source_id, sink_id, src_idx_it->second, sink_idx_it->second, fixed, &used);

        std::vector<std::uint32_t> inside;
        inside.reserve(inside_idx_set.size());
        for (const std::uint32_t idx : inside_idx_set) {
            inside.push_back(idx);
        }
        std::sort(inside.begin(), inside.end());

        debug.source = source_id;
        debug.sink = sink_id;
        debug.inside_node_count = inside.size();
        if (inside.empty()) {
            // A snarl inside a loop has nothing between its boundaries: the copies pass through them. It is
            // kept as the seed of a site around the loop (see widen_loop_sites).
            if (loop_pairs != nullptr && (walked_twice(source_id) || walked_twice(sink_id)))
                loop_pairs->emplace_back(source_id, sink_id);
            if (debug_entries != nullptr) {
                debug_entries->push_back(std::move(debug)); // empty inside
            }
            continue;
        }

        BubbleCandidateIdx candidate;
        candidate.source = src_idx_it->second;
        candidate.sink = sink_idx_it->second;
        candidate.inside = std::move(inside);
        candidate.source_end = used.source_end;
        candidate.sink_end = used.sink_end;
        add_or_replace_candidate(dedup, std::move(candidate), &collapsed, ref_crosses);
        if (debug_entries != nullptr) {
            debug_entries->push_back(std::move(debug)); // acceptance decided after filtering
        }
    }
    return collapsed;
}

void assign_ids(std::vector<Bubble>& bubbles) {
    for (std::size_t i = 0; i < bubbles.size(); ++i) {
        bubbles[i].id = i + 1;
    }
}

} // namespace

BubbleCallReport call_bubbles_report(const Graph& graph, const BubbleCallOptions& options) {
    const PackedGraph packed = pack_graph(graph);
    CandidateMap dedup;
    dedup.reserve(packed.node_ids.size());

    std::vector<SnarlDebugEntry> snarl_debug;
    auto* debug_ptr = options.collect_snarl_debug ? &snarl_debug : nullptr;
    BubbleRunCounts counts;

    if (graph.paths.empty()) {
        throw std::runtime_error(
            "Bubble refinement requires paths (P/W) loaded; rerun with path-aware parsing enabled");
    }

    std::vector<PathIndex> path_indexes;
    path_indexes.reserve(graph.paths.size());
    for (const auto& p : graph.paths) {
        path_indexes.push_back(build_path_index(p));
    }

    // Snarl boundary pairs: precomputed cactus pairs (the internal default, supplied by the command) win;
    // otherwise an external vg snarls JSONL.
    std::vector<std::pair<std::string, std::string>> snarl_pairs;
    std::vector<SnarlInward> snarl_inward = options.snarl_inward_override;
    if (options.snarl_source_supplied || !options.snarl_pairs_override.empty()) {
        // An empty override from a supplied source is a valid answer (a linear graph has no snarl): an
        // empty table, not an error.
        snarl_pairs = options.snarl_pairs_override;
    } else if (!options.snarls_input_path.empty()) {
        snarl_inward.clear();
        snarl_pairs = read_snarl_pairs_jsonl(options.snarls_input_path, &snarl_inward);
    } else {
        throw std::runtime_error(
            "call_bubbles_report: no snarl source — provide snarl_pairs_override (internal cactus "
            "finder) or snarls_input_path (--snarls)");
    }
    if (snarl_inward.size() != snarl_pairs.size()) snarl_inward.clear();   // no orientation: voted from the walks
    // The reference path, so allele support can be split into reference and alternate: exact name, then a
    // unique substring; an ambiguous name is left unresolved and every allele reported as alternate. It also
    // settles which of two snarls on one node pair is kept, when their interiors are the same size.
    std::size_t ref_idx = graph.paths.size();
    if (!options.reference_path.empty()) {
        for (std::size_t i = 0; i < graph.paths.size(); ++i)
            if (graph.paths[i].name == options.reference_path) { ref_idx = i; break; }
        if (ref_idx == graph.paths.size()) {
            std::size_t hits = 0, last = 0;
            for (std::size_t i = 0; i < graph.paths.size(); ++i)
                if (graph.paths[i].name.find(options.reference_path) != std::string::npos) { ++hits; last = i; }
            if (hits == 1) ref_idx = last;
        }
    }
    counts.snarls = snarl_pairs.size();
    std::vector<std::pair<std::string, std::string>> loop_pairs;
    const std::size_t collapsed = collect_candidates_for_pairs(
        graph, packed, snarl_pairs, snarl_inward, path_indexes, dedup, debug_ptr,
        ref_idx < graph.paths.size() ? &graph.paths[ref_idx] : nullptr, &loop_pairs);
    if (collapsed > 0) {
        std::cerr << cli::log_prefix(options.log_module) << "warning: snarls sharing both boundary nodes with another: "
                  << collapsed << "; a site is named by its boundary nodes, so the one with more interior nodes (on a "
                     "tie, the one the reference crosses) is kept and the other's variation is not reported\n";
    }

    std::vector<Bubble> bubbles;
    bubbles.reserve(dedup.size());
    for (auto& [key, cand] : dedup) {
        (void)key;
        Bubble bubble;
        bubble.source = packed.node_ids[cand.source];
        bubble.sink = packed.node_ids[cand.sink];
        // Read into the site: the source left through its inward side, the sink entered through its.
        // These sides are kept: orienting to the reference below only swaps the pair into its order.
        bubble.source_reverse = !cand.source_end;
        bubble.sink_reverse = cand.sink_end;
        bubble.inside.reserve(cand.inside.size());
        for (const std::uint32_t inside_idx : cand.inside) {
            bubble.inside.push_back(packed.node_ids[inside_idx]);
        }
        bubbles.push_back(std::move(bubble));
    }
    // Reference order and orientation before scoring: the crossings scored are read through each
    // boundary's inward side (see bubble_steps), and the walk cache is keyed by the oriented pair.
    orient_bubbles_to_reference(graph, options.reference_path, bubbles);

    // Sites around loops. The graph can store a tandem copy once, as a loop the haplotypes pass through once
    // per copy. A snarl found inside the loop has boundaries the copies pass through: a haplotype's extra
    // copies then fall outside what the site measures, or the snarl has no interior at all. Such a site is
    // widened along the reference, on each side whose boundary some path walks twice, to the first node that
    // no path walks twice, that every path reaching past it on the reference walks, and that lies in no other
    // site's interior. A path ending before the node (a partial assembly) does not count against it. The sites the widened site covers are taken into it; one it would only partly
    // overlap keeps the site as it was, as does a side with no such node within 100 kb.
    if (ref_idx < graph.paths.size()) {
        std::unordered_set<std::string> multi;                 // nodes some path walks twice in one direction
        std::unordered_map<std::string, std::size_t> walkers;  // paths walking each node
        for (std::size_t p = 0; p < path_indexes.size(); ++p)
            for (const auto& [n, pos] : path_indexes[p].positions) {
                ++walkers[n];
                if (repeats_one_way(graph.paths[p], pos)) multi.insert(n);
            }
        std::vector<std::pair<std::string, std::string>> seeds = loop_pairs;
        for (const Bubble& b : bubbles)
            if (multi.count(b.source) || multi.count(b.sink)) seeds.emplace_back(b.source, b.sink);
        const PathRecord& ref = graph.paths[ref_idx];
        std::unordered_map<std::string, std::size_t> ref_pos;   // first step of each node on the reference
        for (std::size_t k = 0; k < ref.steps.size(); ++k) ref_pos.emplace(ref.steps[k].node_id, k);
        // Paths reaching past each reference step: a path's reach runs from the first to the last reference
        // step among the nodes it walks.
        std::vector<std::size_t> reaching(ref.steps.size(), 0);
        {
            std::vector<long long> diff(ref.steps.size() + 1, 0);
            for (const PathIndex& ix : path_indexes) {
                std::size_t lo = SIZE_MAX, hi = 0;
                for (const auto& [n, pos] : ix.positions) {
                    (void)pos;
                    const auto it = ref_pos.find(n);
                    if (it == ref_pos.end()) continue;
                    lo = std::min(lo, it->second);
                    hi = std::max(hi, it->second);
                }
                if (lo > hi) continue;
                ++diff[lo];
                --diff[hi + 1];
            }
            long long run = 0;
            for (std::size_t k = 0; k < ref.steps.size(); ++k) reaching[k] = static_cast<std::size_t>(run += diff[k]);
        }
        const auto len_of = [&](const std::string& id) -> std::size_t {
            const auto it = graph.nodes.find(id);
            return it == graph.nodes.end() ? 0 : it->second.sequence.size();
        };
        constexpr std::size_t kMaxWiden = 100000;
        constexpr std::size_t kNone = SIZE_MAX;
        // Built once, so a locus where every site is a seed (a haplotype walking it twice) stays linear: the
        // reference's offsets, the nearest step on each side whose node could be a boundary (walked twice by no
        // path, walked by every path reaching it), and which sites hold each node.
        const std::size_t nref = ref.steps.size();
        std::vector<std::size_t> off(nref + 1, 0), prev_ok(nref, kNone), next_ok(nref, kNone);
        for (std::size_t k = 0; k < nref; ++k) off[k + 1] = off[k] + len_of(ref.steps[k].node_id);
        const auto could_bound = [&](std::size_t k) {
            const std::string& n = ref.steps[k].node_id;
            return !multi.count(n) && walkers[n] >= reaching[k];
        };
        for (std::size_t k = 0; k < nref; ++k) prev_ok[k] = could_bound(k) ? k : (k > 0 ? prev_ok[k - 1] : kNone);
        for (std::size_t k = nref; k-- > 0;) next_ok[k] = could_bound(k) ? k : (k + 1 < nref ? next_ok[k + 1] : kNone);
        std::vector<char> dead(bubbles.size(), 0);
        std::unordered_map<std::string, std::size_t> interior_count;          // node -> sites holding it inside
        std::unordered_map<std::string, std::vector<std::size_t>> sites_of;  // node -> sites holding it at all
        std::unordered_map<std::string, std::size_t> site_at;                 // "a\tb" (either order) -> site
        const auto index_site = [&](std::size_t i) {
            const Bubble& b = bubbles[i];
            for (const std::string& n : b.inside) { ++interior_count[n]; sites_of[n].push_back(i); }
            sites_of[b.source].push_back(i);
            sites_of[b.sink].push_back(i);
            site_at[b.source + "\t" + b.sink] = i;
            site_at[b.sink + "\t" + b.source] = i;
        };
        const auto drop_site = [&](std::size_t i) {
            dead[i] = 1;
            for (const std::string& n : bubbles[i].inside) --interior_count[n];
            site_at.erase(bubbles[i].source + "\t" + bubbles[i].sink);
            site_at.erase(bubbles[i].sink + "\t" + bubbles[i].source);
        };
        for (std::size_t i = 0; i < bubbles.size(); ++i) index_site(i);
        std::size_t widened = 0;
        for (const auto& [src, snk] : seeds) {
            const auto a = ref_pos.find(src), z = ref_pos.find(snk);
            if (a == ref_pos.end() || z == ref_pos.end()) continue;
            const std::size_t lo = std::min(a->second, z->second), hi = std::max(a->second, z->second);
            const std::string lo_node = ref.steps[lo].node_id;
            const std::string hi_node = ref.steps[hi].node_id;
            std::unordered_set<std::string> own;   // the seed's own interior, when it is a site
            if (const auto it = site_at.find(src + "\t" + snk); it != site_at.end())
                own.insert(bubbles[it->second].inside.begin(), bubbles[it->second].inside.end());
            const auto in_other_site = [&](const std::string& n) {
                const auto c = interior_count.find(n);
                return c != interior_count.end() && c->second > (own.count(n) ? 1u : 0u);
            };
            std::string left = lo_node, right = hi_node;
            if (multi.count(lo_node)) {
                left.clear();
                for (std::size_t j = lo > 0 ? prev_ok[lo - 1] : kNone; j != kNone;
                     j = j > 0 ? prev_ok[j - 1] : kNone) {
                    if (off[lo] - off[j + 1] > kMaxWiden) break;
                    if (!in_other_site(ref.steps[j].node_id)) { left = ref.steps[j].node_id; break; }
                }
            }
            if (multi.count(hi_node)) {
                right.clear();
                for (std::size_t j = hi + 1 < nref ? next_ok[hi + 1] : kNone; j != kNone;
                     j = j + 1 < nref ? next_ok[j + 1] : kNone) {
                    if (off[j] - off[hi + 1] > kMaxWiden) break;
                    if (!in_other_site(ref.steps[j].node_id)) { right = ref.steps[j].node_id; break; }
                }
            }
            if (left.empty() || right.empty() || left == right) continue;
            if (site_at.count(left + "\t" + right)) continue;   // widened already, from another snarl of the same loop
            const auto li = packed.node_idx_of.find(left), ri = packed.node_idx_of.find(right);
            if (li == packed.node_idx_of.end() || ri == packed.node_idx_of.end()) continue;
            SnarlInward used;
            const std::unordered_set<std::uint32_t> in_idx = interior_indices_for_pair(
                graph, packed, path_indexes, left, right, li->second, ri->second, nullptr, &used);
            if (in_idx.empty()) continue;
            Bubble w;
            w.source = left;
            w.sink = right;
            w.source_reverse = !used.source_end;
            w.sink_reverse = used.sink_end;
            for (const std::uint32_t i : in_idx) w.inside.push_back(packed.node_ids[i]);
            std::sort(w.inside.begin(), w.inside.end());
            std::unordered_set<std::string> covered(w.inside.begin(), w.inside.end());
            covered.insert(left);
            covered.insert(right);
            const std::unordered_set<std::string> w_inside(w.inside.begin(), w.inside.end());
            // Only the sites holding one of its nodes can be covered by it or overlap it.
            std::unordered_set<std::size_t> near;
            for (const std::string& n : covered)
                if (const auto it = sites_of.find(n); it != sites_of.end())
                    for (const std::size_t i : it->second) if (!dead[i]) near.insert(i);
            std::vector<std::size_t> taken;
            bool partial = false;
            for (const std::size_t i : near) {
                const Bubble& b = bubbles[i];
                const bool contained = covered.count(b.source) && covered.count(b.sink) &&
                    std::all_of(b.inside.begin(), b.inside.end(), [&](const std::string& n) { return covered.count(n) > 0; });
                if (contained) taken.push_back(i);
                else if (std::any_of(b.inside.begin(), b.inside.end(), [&](const std::string& n) { return w_inside.count(n) > 0; }))
                    partial = true;
            }
            if (partial) continue;
            for (const std::size_t i : taken) drop_site(i);
            bubbles.push_back(std::move(w));
            dead.push_back(0);
            index_site(bubbles.size() - 1);
            ++widened;
        }
        counts.widened = widened;
        if (widened > 0) {
            std::vector<Bubble> alive;
            alive.reserve(bubbles.size());
            for (std::size_t i = 0; i < bubbles.size(); ++i) if (!dead[i]) alive.push_back(std::move(bubbles[i]));
            bubbles = std::move(alive);
        }
        if (widened > 0) {
            orient_bubbles_to_reference(graph, options.reference_path, bubbles);
        }
    }

    // Per-bubble support: number of crossing paths, internal bp span (min/max), the count of
    // paths whose internal span reaches --min-variant-bp, and an inversion signal (an internal
    // node seen in both orientations across paths).

    // Kept by the first scoring pass, for fusing by placement (see fuse_ambiguous_placements).
    std::unordered_map<std::string, SiteWalks> walk_cache;
    bool collect_walks = options.fuse_ambiguous;
    const auto compute_bubble_metrics = [&](std::vector<Bubble>& target, const char* progress_label) {
        // An empty label suppresses the bar (honors --quiet).
        cli::ProgressBar progress(options.quiet ? "" : progress_label, target.size());
        for (auto& bubble : target) {
            std::size_t supported_paths = 0;
            std::size_t min_inside_bp = std::numeric_limits<std::size_t>::max();
            std::size_t max_inside_bp = 0;
            std::size_t long_path_support = 0;
            bool has_inside_bp = false;
            bool inversion_signal = false;

            std::unordered_set<std::string> inside_nodes(bubble.inside.begin(), bubble.inside.end());
            std::unordered_map<std::string, unsigned char> orientation_mask;
            orientation_mask.reserve(bubble.inside.size() * 2);

            std::unordered_map<std::string, std::size_t> allele_counts;
            std::string ref_signature;
            for (std::size_t p_idx = 0; p_idx < path_indexes.size(); ++p_idx) {
                // `bubble_steps` falls back to the adjacent source/sink pair when a path crosses with no
                // interior: the pure-deletion allele, which counts as support with an inside span of 0.
                //
                // A path crossing the site more than once (a tandem array whose boundaries recur, a
                // duplicated site) yields one allele spanning every crossing through the interior: the array
                // is the variant, and its length is what --min-variant-bp and --max-variant-bp have to see.
                // An empty crossing before the first or after the last interior crossing is left out; one
                // between two of them is inside the allele, and one alone is the deletion allele.
                const auto steps_opt =
                    bubble_steps(graph.paths[p_idx], path_indexes[p_idx], bubble, nullptr, nullptr, true);
                if (!steps_opt.has_value() || steps_opt->size() < 2) {
                    continue;
                }
                ++supported_paths;
                const std::vector<PathStep>& steps = *steps_opt;
                // The walk itself, so support can be reported per allele rather than per traversal. Between two
                // passes a path may leave the site; what it walks out there is not this site's allele, so each
                // such stretch reads as one mark.
                {
                    std::string sig;
                    sig.reserve(steps.size() * 8);
                    bool outside = false;
                    for (std::size_t i = 0; i < steps.size(); ++i) {
                        const PathStep& st = steps[i];
                        const bool own = i == 0 || i + 1 == steps.size() || st.node_id == bubble.source ||
                                         st.node_id == bubble.sink || inside_nodes.count(st.node_id) > 0;
                        if (!own) {
                            if (!outside) sig += "~,";
                            outside = true;
                            continue;
                        }
                        outside = false;
                        sig += st.node_id;
                        sig += st.reverse ? '-' : '+';
                        sig += ',';
                    }
                    const bool first_walk = ++allele_counts[sig] == 1;
                    if (p_idx == ref_idx) ref_signature = sig;
                    if (collect_walks) {
                        SiteWalks& w = walk_cache[bubble.source + "\t" + bubble.sink];
                        if (p_idx == ref_idx) w.ref = steps;
                        else if (first_walk) w.alts.push_back(steps);
                    }
                }

                // Each crossing is measured on its own: the passes that follow each other directly are one (a tandem
                // copy, which includes the boundaries the copies share), and a path that leaves the site between two
                // passes crosses it twice, a copy elsewhere. What it walks outside the site is not counted.
                bool long_crossing = false;
                for (const auto& [from, to] : site_passes(steps, bubble, inside_nodes)) {
                    std::size_t inside_bp = 0;
                    for (std::size_t i = from + 1; i < to; ++i) {
                        const auto& step = steps[i];
                        const auto node_it = graph.nodes.find(step.node_id);
                        if (node_it != graph.nodes.end()) {
                            inside_bp += node_it->second.sequence.size();
                        }
                        if (step.node_id == bubble.source || step.node_id == bubble.sink) continue;
                        auto& mask = orientation_mask[step.node_id];
                        mask |= static_cast<unsigned char>(step.reverse ? 0x2 : 0x1);
                        if (mask == 0x3) {
                            inversion_signal = true;
                        }
                    }
                    min_inside_bp = std::min(min_inside_bp, inside_bp);
                    max_inside_bp = std::max(max_inside_bp, inside_bp);
                    has_inside_bp = true;
                    long_crossing = long_crossing || (options.min_variant_bp > 0 && inside_bp >= options.min_variant_bp);
                }
                if (long_crossing) {
                    ++long_path_support;
                }
            }

            // path_support counts traversals, which on a full panel say little about any allele; these say
            // what the traversals contain.
            bubble.distinct_alleles = allele_counts.size();
            bubble.ref_allele_support = 0;
            std::size_t alt_max = 0, alt_min = 0;
            bool any_alt = false;
            for (const auto& [sig, n] : allele_counts) {
                if (!ref_signature.empty() && sig == ref_signature) { bubble.ref_allele_support = n; continue; }
                alt_max = std::max(alt_max, n);
                alt_min = any_alt ? std::min(alt_min, n) : n;
                any_alt = true;
            }
            bubble.alt_allele_support_max = alt_max;
            bubble.alt_allele_support_min = any_alt ? alt_min : 0;
            bubble.path_support = supported_paths;
            bubble.min_inside_bp = has_inside_bp ? min_inside_bp : 0;
            bubble.max_inside_bp = has_inside_bp ? max_inside_bp : 0;
            bubble.long_path_support = long_path_support;
            bubble.inversion_signal = inversion_signal;
            progress.tick();
        }
    };

    counts.candidates = bubbles.size();
    compute_bubble_metrics(bubbles, (cli::log_prefix(options.log_module) + "scoring sites").c_str());
    collect_walks = false;

    std::vector<Bubble> non_snp_bubbles;
    non_snp_bubbles.reserve(bubbles.size());
    for (const auto& bubble : bubbles) {
        if (bubble.max_inside_bp > 1 || bubble.inversion_signal) {
            non_snp_bubbles.push_back(bubble);
        }
    }

    // Candidate metrics by endpoint, captured before filtering so the debug TSV can report
    // them for rejected candidates too.
    struct DebugMetrics {
        std::size_t n_paths = 0;
        std::size_t min_inside_bp = 0;
        std::size_t inside_node_count = 0;
        std::size_t long_path_support = 0;
        bool inversion_signal = false;
    };
    std::unordered_map<std::string, DebugMetrics> metrics_by_endpoint;
    metrics_by_endpoint.reserve(bubbles.size() * 2);
    for (const auto& bubble : bubbles) {
        metrics_by_endpoint[endpoint_key(bubble.source, bubble.sink)] =
            DebugMetrics{bubble.path_support, bubble.min_inside_bp, bubble.inside.size(),
                         bubble.long_path_support, bubble.inversion_signal};
    }

    // Boundaries take reference order before anything reads them as an interval: merging is defined in
    // reference coordinates and needs oriented pairs.
    orient_bubbles_to_reference(graph, options.reference_path, bubbles);

    // Every filter, applied to whatever the current bubble set is: merging creates bubbles whose metrics are
    // not those of their parts (two 1 bp bubbles 10 bp apart fuse into one spanning 12 bp), so a set filtered
    // before merging is not a filtered set afterwards.
    auto apply_filters = [&](std::vector<Bubble>& bs, bool size_floor = true) {
        if (options.min_path_support > 0) {
            bs.erase(std::remove_if(bs.begin(), bs.end(),
                                    [&](const Bubble& b) { return b.path_support < options.min_path_support; }),
                     bs.end());
        }
        if (options.min_alt_support > 0) {
            bs.erase(std::remove_if(bs.begin(), bs.end(),
                                    [&](const Bubble& b) {
                                        return b.alt_allele_support_max < options.min_alt_support;
                                    }),
                     bs.end());
        }
        if (size_floor && options.min_variant_bp > 0) {
            bs.erase(std::remove_if(bs.begin(), bs.end(),
                                    [&](const Bubble& b) {
                                        return b.long_path_support == 0 && !b.inversion_signal;
                                    }),
                     bs.end());
        }
        if (options.max_variant_bp > 0) {
            // Bounds the per-bubble walk cost on hypervariable regions, at the price of SVs larger than the
            // cap.
            bs.erase(std::remove_if(bs.begin(), bs.end(),
                                    [&](const Bubble& b) { return b.max_inside_bp > options.max_variant_bp; }),
                     bs.end());
        }
    };

    // The reference's node offsets, for every merge below.
    const auto ref_offsets_all = reference_offsets(graph, options.reference_path);

    // An event the graph splits across adjacent snarls, none reaching the size floor on its own. The pieces
    // that change length are kept aside and fused at the end (below); equal-length pieces only fall inside.
    std::vector<Bubble> sub_floor;
    if (options.min_variant_bp > 0)
        for (const Bubble& b : bubbles)
            if (b.long_path_support == 0 && !b.inversion_signal && b.max_inside_bp > b.min_inside_bp)
                sub_floor.push_back(b);

    // Every candidate as scored that only the size floor removes: fusing by placement takes in sites under
    // the size floor whose variants could be placed in a retained site's. A site another filter removes
    // (the support floors, --max-variant-bp) stays out.
    std::vector<Bubble> all_candidates = options.fuse_ambiguous ? bubbles : std::vector<Bubble>{};
    apply_filters(all_candidates, false);
    apply_filters(bubbles);
    counts.passed = bubbles.size();

    if (options.merge_nearby_bp > 0 && bubbles.size() > 1) {
        std::vector<std::vector<Bubble>> parts;
        const std::vector<Bubble> unmerged = bubbles;
        bubbles = merge_nearby_bubbles(graph, packed, path_indexes, bubbles, options.merge_nearby_bp,
                                       options.reference_path, &parts, &ref_offsets_all);
        // Only the merged sites are scored again: the others are the sites already scored.
        {
            std::vector<Bubble> merged_only;
            for (std::size_t i = 0; i < bubbles.size(); ++i)
                if (parts[i].size() >= 2) merged_only.push_back(bubbles[i]);
            compute_bubble_metrics(merged_only, (cli::log_prefix(options.log_module) + "scoring merged sites").c_str());
            for (std::size_t i = 0, m = 0; i < bubbles.size(); ++i)
                if (parts[i].size() >= 2) bubbles[i] = std::move(merged_only[m++]);
        }
        // A merged site has to clear the same bars its parts did, and must not be one the overlap pass would
        // drop together with its parts, or that would drop a retained site it overlaps (see
        // fusion_loses_site). A merged site that fails is not undone whole: the longest run of its parts
        // from the left that passes is merged, and so on along the rest, each part that merges with none
        // kept as it was.
        const auto passes = [&](const Bubble& b, const std::vector<Bubble>& members, bool& lost) {
            std::vector<Bubble> one{b};
            apply_filters(one);
            if (one.empty()) return false;
            std::unordered_set<std::string> member_keys;
            for (const Bubble& p : members) member_keys.insert(p.source + "\t" + p.sink);
            lost = fusion_loses_site(b, member_keys, unmerged);
            return !lost;
        };
        constexpr std::size_t kMaxRetried = 32;   // longer chains are kept as their parts
        std::vector<Bubble> kept;
        std::size_t undone = 0, overlapping = 0;
        for (std::size_t i = 0; i < bubbles.size(); ++i) {
            if (parts[i].size() < 2) {   // not a fusion: kept if it passes
                std::vector<Bubble> one{bubbles[i]};
                apply_filters(one);
                if (!one.empty()) kept.push_back(std::move(one.front()));
                continue;
            }
            bool lost = false;
            if (passes(bubbles[i], parts[i], lost)) { kept.push_back(std::move(bubbles[i])); ++counts.merged; continue; }
            if (lost) ++overlapping; else ++undone;
            const std::vector<Bubble>& ps = parts[i];
            if (ps.size() > kMaxRetried) {
                for (const Bubble& p : ps) kept.push_back(p);
                continue;
            }
            for (std::size_t a = 0; a < ps.size();) {
                std::size_t end = a;   // last part merged with ps[a]
                for (std::size_t z = ps.size() - 1; z > a; --z) {
                    if (a == 0 && z == ps.size() - 1) continue;   // the whole chain failed above
                    const std::vector<Bubble> run(ps.begin() + static_cast<long>(a), ps.begin() + static_cast<long>(z) + 1);
                    std::vector<Bubble> f = merge_nearby_bubbles(graph, packed, path_indexes, run,
                                                                 std::numeric_limits<std::size_t>::max() / 4,
                                                                 options.reference_path, nullptr, &ref_offsets_all);
                    if (f.size() != 1) continue;
                    compute_bubble_metrics(f, "");
                    bool l = false;
                    if (!passes(f.front(), run, l)) continue;
                    kept.push_back(std::move(f.front()));
                    ++counts.merged;
                    end = z;
                    break;
                }
                if (end == a) kept.push_back(ps[a]);
                a = end + 1;
            }
        }
        if (overlapping > 0) {
            std::cerr << cli::log_prefix(options.log_module) << "warning: merged sites that would lose a site to the overlap "
                         "step: " << overlapping << "; merged only the runs of their parts that do not\n";
        }
        std::sort(kept.begin(), kept.end(), bubble_endpoint_less);
        bubbles = std::move(kept);
        if (undone > 0) {
            std::cerr << cli::log_prefix(options.log_module) << "warning: merged sites failing the filters their parts "
                         "passed: " << undone << "; merged only the runs of their parts that pass\n";
        }
    }

    if (options.fuse_ambiguous && !bubbles.empty()) {
        std::size_t fused = 0;
        std::vector<std::vector<Bubble>> parts;
        std::vector<Bubble> next = fuse_ambiguous_placements(graph, packed, path_indexes, bubbles, all_candidates,
                                                             options.reference_path, walk_cache, &parts, &fused);
        if (fused > 0) {
            // Only the fused sites changed; the others keep their scores. A fused site that fails a filter is
            // undone on its own, its retained parts put back: each passed the filters itself.
            bubbles.clear();
            std::size_t undone = 0;
            for (std::size_t i = 0; i < next.size(); ++i) {
                if (parts[i].empty()) { bubbles.push_back(std::move(next[i])); continue; }
                std::vector<Bubble> one{std::move(next[i])};
                compute_bubble_metrics(one, "");
                std::vector<Bubble> check = one;
                apply_filters(check);
                if (check.size() == 1) { bubbles.push_back(std::move(one.front())); continue; }
                ++undone;
                --fused;
                for (Bubble& p : parts[i]) bubbles.push_back(std::move(p));
            }
            if (undone > 0)
                std::cerr << cli::log_prefix(options.log_module) << "warning: sites fused over placements failing the "
                             "filters: " << undone << "; kept the parts unfused\n";
        }
        std::sort(bubbles.begin(), bubbles.end(), bubble_endpoint_less);
        counts.fused = fused;
    }

    // Always reported, quiet or not: the interior of every affected site is panel-derived rather than
    // graph-derived, so an allele no one in the panel carries is missing from it and from everything
    // computed from `inside`, with nothing in the output to show it. Counted after merging, so the fused
    // bubbles' searches are included.
    if (const std::size_t truncated = g_interior_traversal_truncated.exchange(0); truncated > 0) {
        std::cerr << cli::log_prefix(options.log_module) << "warning: sites whose interior search hit its "
                  << kMaxInteriorHandles << "-handle limit: " << truncated
                  << "; their interiors are what the haplotypes walk, so a branch no haplotype walks is not part"
                     " of them\n";
    }

    // Overlapping retained sites: two bubbles claiming the same interior node describe the same sequence
    // twice, and panphorte refuses to fold both. The smaller sites are kept and the one enclosing them
    // dropped, toward finer resolution; this also drops an oversized site overlapping real ones, such as the
    // closing snarl of a component whose decomposition roots in a cycle. Reported, since it loses sequence.
    if (bubbles.size() > 1) {
        std::sort(bubbles.begin(), bubbles.end(), [](const Bubble& a, const Bubble& b) {
            if (a.inside.size() != b.inside.size()) return a.inside.size() < b.inside.size();  // smaller first
            return bubble_endpoint_less(a, b);   // deterministic among equals
        });
        std::unordered_set<std::string> claimed;
        std::vector<Bubble> kept;
        std::size_t dropped = 0;
        for (Bubble& b : bubbles) {
            const std::string* clash = nullptr;
            for (const std::string& n : b.inside) {
                if (claimed.count(n)) { clash = &n; break; }
            }
            if (clash != nullptr) {
                ++dropped;
                // Shown even with -q: the site's variation is lost.
                std::cerr << cli::log_prefix(options.log_module) << "warning: dropped overlapping site: " << b.source
                          << ".." << b.sink << " (" << b.inside.size() << " interior nodes; node " << *clash
                          << " is in a smaller site)\n";
                continue;
            }
            for (const std::string& n : b.inside) claimed.insert(n);
            kept.push_back(std::move(b));
        }
        counts.overlaps_dropped = dropped;
        bubbles = std::move(kept);
    }

    // Rescue of split events, after every other step so it changes none of them. The pieces no retained
    // site holds are fused within --merge-nearby-bp (else --min-variant-bp, as refine re-snarls with merging
    // off), never across a retained site. A fusion is kept when an allele differs in length from the
    // reference's by the floor -- what call reports -- it passes the filters, and it holds no node of a
    // retained or earlier rescued site.
    const std::size_t rescue_bp = options.merge_nearby_bp > 0 ? options.merge_nearby_bp : options.min_variant_bp;
    if (sub_floor.size() > 1 && rescue_bp > 0) {
        std::unordered_set<std::string> held;
        for (const Bubble& b : bubbles) {
            held.insert(b.inside.begin(), b.inside.end());
            held.insert(b.source);
            held.insert(b.sink);
        }
        const auto& ref_off = ref_offsets_all;
        const auto start_of = [&](const std::string& n) -> std::optional<std::size_t> {
            const auto it = ref_off.find(n);
            return it == ref_off.end() ? std::nullopt : std::optional<std::size_t>(it->second.first);
        };
        std::vector<std::size_t> retained_starts;
        for (const Bubble& b : bubbles)
            if (const auto s = start_of(b.source)) retained_starts.push_back(*s);
        std::sort(retained_starts.begin(), retained_starts.end());
        // Pieces between the same two retained sites form one group: the retained sites starting before it.
        std::map<std::size_t, std::vector<Bubble>> groups;
        for (const Bubble& b : sub_floor) {
            if (std::any_of(b.inside.begin(), b.inside.end(), [&](const std::string& n) { return held.count(n) > 0; }))
                continue;
            const auto s = start_of(b.source);
            if (!s) continue;
            groups[std::lower_bound(retained_starts.begin(), retained_starts.end(), *s) - retained_starts.begin()]
                .push_back(b);
        }
        // The reference's interior length at a site, when the reference crosses it.
        const auto ref_inside_bp = [&](const Bubble& b) -> std::optional<std::size_t> {
            if (ref_idx >= graph.paths.size()) return std::nullopt;
            const auto steps = bubble_steps(graph.paths[ref_idx], path_indexes[ref_idx], b, nullptr, nullptr, true);
            if (!steps.has_value() || steps->size() < 2) return std::nullopt;
            const std::unordered_set<std::string> in(b.inside.begin(), b.inside.end());
            std::size_t bp = 0;
            for (std::size_t i = 1; i + 1 < steps->size(); ++i) {
                const auto it = graph.nodes.find((*steps)[i].node_id);
                if (in.count((*steps)[i].node_id) && it != graph.nodes.end()) bp += it->second.sequence.size();
            }
            return bp;
        };
        std::vector<Bubble> rescued;
        for (auto& [gap, pieces] : groups) {
            (void)gap;
            if (pieces.size() < 2) continue;
            std::sort(pieces.begin(), pieces.end(), bubble_endpoint_less);
            std::vector<std::vector<Bubble>> parts;
            std::vector<Bubble> fused = merge_nearby_bubbles(graph, packed, path_indexes, pieces, rescue_bp,
                                                             options.reference_path, &parts, &ref_offsets_all);
            for (std::size_t i = 0; i < fused.size(); ++i)
                if (parts[i].size() >= 2) rescued.push_back(std::move(fused[i]));
        }
        if (!rescued.empty()) {
            compute_bubble_metrics(rescued, "");
            apply_filters(rescued);   // the bars every other site cleared
            std::size_t kept = 0;
            for (Bubble& b : rescued) {
                const auto ref_bp = ref_inside_bp(b);
                if (!ref_bp) continue;
                if (std::max(b.max_inside_bp - std::min(b.max_inside_bp, *ref_bp),
                             *ref_bp - std::min(*ref_bp, b.min_inside_bp)) < options.min_variant_bp)
                    continue;
                if (std::any_of(b.inside.begin(), b.inside.end(), [&](const std::string& n) { return held.count(n) > 0; }))
                    continue;
                held.insert(b.inside.begin(), b.inside.end());
                bubbles.push_back(std::move(b));
                ++kept;
            }
            counts.rescued = kept;
            // Rescued sites then merge with their neighbours as any site does. Only a fusion holding a rescued
            // site is made; it must pass the filters and take no node of a site outside it.
            if (kept > 0 && options.merge_nearby_bp > 0 && bubbles.size() > 1) {
                std::unordered_set<std::string> rescued_keys;
                for (std::size_t i = bubbles.size() - kept; i < bubbles.size(); ++i)
                    rescued_keys.insert(endpoint_key(bubbles[i].source, bubbles[i].sink));
                std::sort(bubbles.begin(), bubbles.end(), bubble_endpoint_less);
                std::vector<std::vector<Bubble>> parts;
                std::vector<Bubble> merged = merge_nearby_bubbles(graph, packed, path_indexes, bubbles,
                                                                  options.merge_nearby_bp, options.reference_path, &parts,
                                                                  &ref_offsets_all);
                std::vector<Bubble> next;
                for (std::size_t i = 0; i < merged.size(); ++i) {
                    const bool has_rescued = parts[i].size() >= 2 &&
                        std::any_of(parts[i].begin(), parts[i].end(), [&](const Bubble& p) {
                            return rescued_keys.count(endpoint_key(p.source, p.sink)) > 0;
                        });
                    if (has_rescued) {
                        std::vector<Bubble> one{merged[i]};
                        compute_bubble_metrics(one, "");
                        apply_filters(one);
                        if (!one.empty()) {
                            std::unordered_set<std::string> member;
                            for (const Bubble& p : parts[i]) {
                                member.insert(p.inside.begin(), p.inside.end());
                                member.insert(p.source);
                                member.insert(p.sink);
                            }
                            const bool takes = std::any_of(one.front().inside.begin(), one.front().inside.end(),
                                [&](const std::string& n) { return held.count(n) > 0 && member.count(n) == 0; });
                            if (!takes) { next.push_back(std::move(one.front())); continue; }
                        }
                    }
                    for (Bubble& p : parts[i]) next.push_back(std::move(p));
                }
                bubbles = std::move(next);

            }
        }
    }

    std::sort(bubbles.begin(), bubbles.end(), bubble_endpoint_less);
    // Ids follow the reference. On the graph bubble sorts, node ids already do; with --snarls the graph is
    // taken as it is, so the sites the reference reaches take their places in node order sorted by where
    // the reference reaches them, and the others keep theirs.
    if (!options.snarls_input_path.empty() && !ref_offsets_all.empty()) {
        const auto at = [&](const Bubble& b) -> std::optional<std::pair<std::size_t, std::size_t>> {
            const auto s = ref_offsets_all.find(b.source), k = ref_offsets_all.find(b.sink);
            if (s == ref_offsets_all.end() || k == ref_offsets_all.end()) return std::nullopt;
            return std::make_pair(s->second.first, k->second.first);
        };
        std::vector<std::size_t> slots;
        std::vector<Bubble> reached;
        for (std::size_t i = 0; i < bubbles.size(); ++i)
            if (at(bubbles[i])) { slots.push_back(i); reached.push_back(bubbles[i]); }
        std::stable_sort(reached.begin(), reached.end(), [&](const Bubble& x, const Bubble& y) { return *at(x) < *at(y); });
        for (std::size_t j = 0; j < slots.size(); ++j) bubbles[slots[j]] = std::move(reached[j]);
    }
    assign_ids(bubbles);

    // Inversion sites. A site needs boundary nodes every path crosses once in one direction. A large
    // inversion whose breakpoints lie in inverted repeats the graph collapses has none there: every path
    // walks the repeat nodes both ways, so the snarls at the breakpoints have no interior a path crosses.
    // It is read from the walks instead (find_path_inversions), and its site runs between backbone flanks
    // outside it. The sites inside keep their own variants and name it as their parent; ids follow the
    // ordinary sites', which keep theirs.
    if (ref_idx < graph.paths.size()) {
        const PathRecord& ref = graph.paths[ref_idx];
        const PathInversionScan scan = find_path_inversions(
            graph, ref, options.min_variant_bp > 0 ? options.min_variant_bp : kPathInversionMinBp);
        const auto ref_off = reference_offsets(graph, ref.name);
        const auto len_of = [&](const std::string& id) -> std::size_t {
            const auto it = graph.nodes.find(id);
            return it == graph.nodes.end() ? 0 : it->second.sequence.size();
        };
        std::unordered_map<std::string, std::size_t> interior_of;   // an ordinary site's interior node -> that site
        for (std::size_t i = 0; i < bubbles.size(); ++i)
            for (const std::string& n : bubbles[i].inside) interior_of.emplace(n, i);
        std::unordered_map<std::string, std::size_t> ref_index;   // reference-unique nodes -> their step
        {
            std::unordered_map<std::string, std::size_t> seen;
            for (const PathStep& st : ref.steps) ++seen[st.node_id];
            for (std::size_t k = 0; k < ref.steps.size(); ++k)
                if (seen[ref.steps[k].node_id] == 1) ref_index.emplace(ref.steps[k].node_id, k);
        }
        // A flank must be backbone: no path walks it twice, and every path whose reference span covers it walks
        // it. A node some path skips is a variant's, and that path would not cross a site between such flanks.
        std::unordered_set<std::string> not_flank;
        {
            // Reference-unique nodes by offset: a path's span covers a run of them, counted with a difference
            // array, and the ones it walks are read from its own steps.
            std::vector<std::pair<std::size_t, std::string>> by_off;
            for (const auto& [n, k] : ref_index) {
                (void)k;
                const auto o = ref_off.find(n);
                if (o != ref_off.end()) by_off.emplace_back(o->second.first, n);
            }
            std::sort(by_off.begin(), by_off.end());
            std::vector<long long> cover_diff(by_off.size() + 1, 0);
            std::unordered_map<std::string, std::size_t> walkers;
            for (std::size_t pi = 0; pi < graph.paths.size(); ++pi) {
                std::unordered_map<std::string, std::size_t> c;
                for (const PathStep& st : graph.paths[pi].steps) if (++c[st.node_id] == 2) not_flank.insert(st.node_id);
                const auto [lo, hi] = scan.span[pi];
                if (lo == SIZE_MAX) continue;
                const auto a = std::lower_bound(by_off.begin(), by_off.end(), std::make_pair(lo, std::string()));
                const auto z = std::upper_bound(by_off.begin(), by_off.end(), hi,
                                                [](std::size_t v, const auto& e) { return v < e.first; });
                if (a < z) { ++cover_diff[a - by_off.begin()]; --cover_diff[z - by_off.begin()]; }
                for (const auto& [n, cnt] : c) {
                    (void)cnt;
                    if (!ref_index.count(n)) continue;
                    const auto o = ref_off.find(n);
                    if (o != ref_off.end() && o->second.first >= lo && o->second.first <= hi) ++walkers[n];
                }
            }
            std::unordered_map<std::string, std::size_t> covering;
            long long run = 0;
            for (std::size_t i = 0; i < by_off.size(); ++i) {
                run += cover_diff[i];
                covering[by_off[i].second] = static_cast<std::size_t>(run);
            }
            for (const PathStep& st : ref.steps)
                if (walkers[st.node_id] < covering[st.node_id]) not_flank.insert(st.node_id);
        }
        std::size_t next_id = bubbles.size();
        std::vector<Bubble> parents;
        for (const PathInversion& inv : scan.inversions) {
            // A reversed run inside one site's interior is that site's to type. One reaching into a boundary
            // is not: a carrier crossing the site end to end reversed reads as its reference allele there.
            bool held = false;
            for (const Bubble& b : bubbles) {
                const auto a = ref_off.find(b.source), z = ref_off.find(b.sink);
                if (a == ref_off.end() || z == ref_off.end()) continue;
                const auto& first = a->second.first <= z->second.first ? a->second : z->second;
                const auto& last = a->second.first <= z->second.first ? z->second : a->second;
                if (first.first + first.second <= inv.inv_lo && inv.inv_hi < last.first) { held = true; break; }
            }
            if (held || inv.left_node.empty() || inv.right_node.empty()) continue;
            // A flank inside a site -- the breakpoint's repeat often is one, and where the graph collapses the
            // two inverted repeats into one site its interior holds both breakpoints -- is moved outward along the
            // reference, left from the left flank and right from the right one, to the first node in no site's
            // interior that is backbone (below). The sites passed over are then inside the inversion site.
            std::string left, right;
            {
                const auto li = ref_index.find(inv.left_node), ri = ref_index.find(inv.right_node);
                if (li != ref_index.end() && ri != ref_index.end() && li->second < ri->second) {
                    for (std::size_t k = li->second + 1; k-- > 0;) {
                        const std::string& n = ref.steps[k].node_id;
                        if (!interior_of.count(n) && !not_flank.count(n)) { left = n; break; }
                    }
                    for (std::size_t k = ri->second; k < ref.steps.size(); ++k) {
                        const std::string& n = ref.steps[k].node_id;
                        if (!interior_of.count(n) && !not_flank.count(n)) { right = n; break; }
                    }
                }
            }
            if (left.empty() || right.empty() || left == right) {
                // Shown even with -q: the inversion is lost.
                std::cerr << cli::log_prefix(options.log_module) << "warning: inversion not made a site: "
                          << inv.left_node << ".." << inv.right_node << " (no flank outside the sites around it)\n";
                continue;
            }
            // Every path crossing the flanks crosses each once (no path walks a flank twice); the interior is
            // what lies between. A partial haplotype ending inside crosses nothing.
            Bubble p;
            p.source = left; p.sink = right;
            std::unordered_set<std::string> inside;
            std::size_t crossing = 0, min_bp = SIZE_MAX, max_bp = 0;
            std::unordered_set<std::size_t> crossed;   // paths crossing the site
            for (std::size_t pi = 0; pi < graph.paths.size(); ++pi) {
                const PathRecord& path = graph.paths[pi];
                std::size_t at_s = SIZE_MAX, at_k = SIZE_MAX;
                for (std::size_t k = 0; k < path.steps.size(); ++k) {
                    if (path.steps[k].node_id == p.source) at_s = k;
                    if (path.steps[k].node_id == p.sink) at_k = k;
                }
                if (at_s == SIZE_MAX || at_k == SIZE_MAX) continue;
                const std::size_t lo = std::min(at_s, at_k), hi = std::max(at_s, at_k);
                std::size_t bp = 0;
                for (std::size_t k = lo + 1; k < hi; ++k) { inside.insert(path.steps[k].node_id); bp += len_of(path.steps[k].node_id); }
                ++crossing; min_bp = std::min(min_bp, bp); max_bp = std::max(max_bp, bp);
                crossed.insert(pi);
            }
            // Support counts the carriers crossing the site, as path_support does: a partial haplotype carrying
            // the inversion but ending inside it crosses nothing.
            std::size_t carrying = 0;
            for (const std::size_t c : inv.carriers) carrying += crossed.count(c);
            inside.erase(p.source); inside.erase(p.sink);
            p.inside.assign(inside.begin(), inside.end());
            std::sort(p.inside.begin(), p.inside.end());
            for (const PathStep& st : ref.steps) {   // boundaries as the reference walks them
                if (st.node_id == p.source) p.source_reverse = st.reverse;
                if (st.node_id == p.sink) p.sink_reverse = st.reverse;
            }
            p.path_support = crossing;
            p.distinct_alleles = (carrying > 0 ? 1 : 0) + (crossing > carrying ? 1 : 0);
            p.alt_allele_support_max = p.alt_allele_support_min = carrying;
            p.ref_allele_support = crossing - carrying;
            p.min_inside_bp = min_bp; p.max_inside_bp = max_bp;
            p.long_path_support = crossing;
            p.inversion_signal = true;
            p.inversion_parent = true;
            parents.push_back(std::move(p));
        }
        // The support floors and --max-variant-bp apply as to any site; the size floor does not, an
        // inversion being a change whatever its length.
        {
            const std::size_t made = parents.size();
            apply_filters(parents);
            if (!options.quiet && parents.size() < made)
                std::cerr << cli::log_prefix(options.log_module) << "inversion sites filtered out: "
                          << made - parents.size() << '\n';
        }
        // Inversion sites are disjoint: of two sharing a node, the one more paths carry is kept. Two may share
        // a boundary, as ordinary sites do, but a boundary in the other's interior is an overlap.
        std::stable_sort(parents.begin(), parents.end(), [](const Bubble& x, const Bubble& y) {
            return x.alt_allele_support_max > y.alt_allele_support_max;
        });
        {
            std::unordered_set<std::string> taken_inside, taken_ends;
            std::vector<Bubble> disjoint;
            for (Bubble& p : parents) {
                const bool overlaps =
                    taken_inside.count(p.source) || taken_inside.count(p.sink) ||
                    std::any_of(p.inside.begin(), p.inside.end(), [&](const std::string& n) {
                        return taken_inside.count(n) > 0 || taken_ends.count(n) > 0;
                    });
                if (overlaps) {
                    // Shown even with -q: the inversion is lost.
                    std::cerr << cli::log_prefix(options.log_module) << "warning: inversion site not made: " << p.source
                              << ".." << p.sink << " (it overlaps one more haplotypes carry)\n";
                    continue;
                }
                taken_inside.insert(p.inside.begin(), p.inside.end());
                taken_ends.insert(p.source); taken_ends.insert(p.sink);
                disjoint.push_back(std::move(p));
            }
            parents = std::move(disjoint);
        }
        std::sort(parents.begin(), parents.end(), bubble_endpoint_less);
        for (Bubble& p : parents) p.id = ++next_id;
        for (Bubble& p : parents) {
            // A site with both boundaries in the inversion site is inside it, with every node it holds,
            // including branches no path crossing both flanks walks.
            const std::unordered_set<std::string> in(p.inside.begin(), p.inside.end());
            const auto within = [&](const std::string& n) { return in.count(n) || n == p.source || n == p.sink; };
            std::set<std::string> grown(p.inside.begin(), p.inside.end());
            std::size_t children = 0;
            for (Bubble& b : bubbles) {
                if (b.parent_id != 0 || !within(b.source) || !within(b.sink)) continue;
                b.parent_id = p.id; ++children;
                grown.insert(b.inside.begin(), b.inside.end());
            }
            grown.erase(p.source); grown.erase(p.sink);
            p.inside.assign(grown.begin(), grown.end());
            if (!options.quiet)
                std::cerr << cli::log_prefix(options.log_module) << "inversion site: " << p.source << ".." << p.sink
                          << "; carriers: " << p.alt_allele_support_max << " of " << p.path_support
                          << " crossing; sites inside: " << children << '\n';
        }
        counts.inversion_sites = parents.size();
        for (Bubble& p : parents) bubbles.push_back(std::move(p));
    }

    std::unordered_set<std::string> final_endpoint_keys;
    final_endpoint_keys.reserve(bubbles.size() * 2);
    for (const auto& bubble : bubbles) {
        final_endpoint_keys.insert(endpoint_key(bubble.source, bubble.sink));
    }

    if (options.collect_snarl_debug) {
        std::size_t next_debug_id = 0;
        std::unordered_set<std::string> debug_endpoint_keys;
        debug_endpoint_keys.reserve(snarl_debug.size() * 2 + bubbles.size() * 2);
        for (const auto& debug : snarl_debug) {
            next_debug_id = std::max(next_debug_id, debug.candidate_id);
            if (!debug.source.empty() && !debug.sink.empty()) {
                debug_endpoint_keys.insert(endpoint_key(debug.source, debug.sink));
            }
        }

        for (auto& debug : snarl_debug) {
            if (debug.source.empty() || debug.sink.empty()) {
                continue;
            }
            const std::string key = endpoint_key(debug.source, debug.sink);
            const auto metrics_it = metrics_by_endpoint.find(key);
            if (metrics_it != metrics_by_endpoint.end()) {
                debug.n_paths = metrics_it->second.n_paths;
                debug.min_inside_bp = metrics_it->second.min_inside_bp;
                if (debug.inside_node_count == 0) {
                    debug.inside_node_count = metrics_it->second.inside_node_count;
                }
                debug.long_path_support = metrics_it->second.long_path_support;
                debug.inversion_signal = metrics_it->second.inversion_signal;
            }
            debug.accepted = final_endpoint_keys.find(key) != final_endpoint_keys.end();
        }

        // Final bubbles produced by merging have no original candidate row; add one each.
        for (const auto& bubble : bubbles) {
            const std::string key = endpoint_key(bubble.source, bubble.sink);
            if (debug_endpoint_keys.find(key) != debug_endpoint_keys.end()) {
                continue;
            }
            SnarlDebugEntry synthetic;
            synthetic.candidate_id = ++next_debug_id;
            synthetic.source = bubble.source;
            synthetic.sink = bubble.sink;
            synthetic.inside_node_count = bubble.inside.size();
            synthetic.n_paths = bubble.path_support;
            synthetic.min_inside_bp = bubble.min_inside_bp;
            synthetic.long_path_support = bubble.long_path_support;
            synthetic.inversion_signal = bubble.inversion_signal;
            synthetic.accepted = true;
            snarl_debug.push_back(std::move(synthetic));
        }
    }

    BubbleCallReport report;
    report.counts = counts;
    report.bubbles = std::move(bubbles);
    report.non_snp_bubbles = std::move(non_snp_bubbles);
    if (options.collect_snarl_debug) {
        report.snarl_debug = std::move(snarl_debug);
    }
    return report;
}

void view_inversion_sites_on_own_nodes(std::vector<Bubble>& bubbles) {
    for (Bubble& parent : bubbles) {
        if (!parent.inversion_parent) continue;
        std::unordered_set<std::string> inner;
        parent.inner_bounds.clear();
        for (const Bubble& b : bubbles) {
            if (b.parent_id != parent.id) continue;
            inner.insert(b.inside.begin(), b.inside.end());
            parent.inner_bounds.emplace_back(b.source, b.sink);
        }
        std::vector<std::string> own;
        for (const std::string& n : parent.inside) if (!inner.count(n)) own.push_back(n);
        parent.inside = std::move(own);
    }
}

bool inner_site_step(const Bubble& bubble, const std::string& from, const std::string& to) {
    for (const auto& [a, b] : bubble.inner_bounds)
        if ((from == a && to == b) || (from == b && to == a)) return true;
    return false;
}

} // namespace panvar
