#include "panvar/panphorte.hpp"
#include "panvar/panphorte_refine.hpp"

#include "panvar/align.hpp"
#include "panvar/bubble_path.hpp"
#include "panvar/bubbles.hpp"
#include "panvar/cli_utils.hpp"
#include "panvar/gfa.hpp"
#include "panvar/gfa_io.hpp"
#include "panvar/graph_sort.hpp"
#include "panvar/graph_utils.hpp"
#include "panvar/gtf.hpp"
#include "panvar/integrated_snarls.hpp"
#include "panvar/output.hpp"
#include "panvar/parallel.hpp"

#include <algorithm>
#include <atomic>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <mutex>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <iterator>

namespace panvar {
namespace {

constexpr std::size_t kMaxGapSteps = 32; // bound interruption length (steps) between copies

std::string step_sequence(const Graph& graph, const PathStep& step) {
    const auto it = graph.nodes.find(step.node_id);
    if (it == graph.nodes.end()) {
        return {};
    }
    return step.reverse ? reverse_complement(it->second.sequence) : it->second.sequence;
}

std::uint64_t hash_sequence(const std::string& s) {
    std::uint64_t h = 1469598103934665603ULL;
    for (const unsigned char c : s) {
        h ^= static_cast<std::uint64_t>(c);
        h *= 1099511628211ULL;
    }
    return h;
}

char orient_char(bool reverse) { return reverse ? '-' : '+'; }

// Per-node sequence hashes (forward and reverse complement) and length, so a step's token is a lookup
// rather than a re-hash.
struct NodeTok {
    std::uint64_t fwd = 0;
    std::uint64_t rc = 0;
    std::size_t len = 0;
};

std::unordered_map<std::string, NodeTok> build_node_tokens(const Graph& graph) {
    std::unordered_map<std::string, NodeTok> out;
    out.reserve(graph.nodes.size() * 2);
    for (const auto& [id, node] : graph.nodes) {
        NodeTok nt;
        nt.fwd = hash_sequence(node.sequence);
        nt.rc = hash_sequence(reverse_complement(node.sequence));
        nt.len = node.sequence.size();
        out.emplace(id, nt);
    }
    return out;
}

// One element of a (possibly interrupted) tandem array, in emission order: a unit
// copy (replaced by the unit node) or an interruption kept as literal path steps.
struct ArrayElement {
    bool is_copy = false;
    std::vector<PathStep> steps; // literal steps for interruptions
    std::string unit;            // a copy spelling other bases than the array's unit (letter case): its own
};

struct TandemArray {
    std::size_t start = 0;       // absolute index into path.steps (array start)
    std::size_t length = 0;      // number of steps replaced
    std::size_t copies = 0;      // unit copies (unit traversals)
    std::size_t interruption_bp = 0;
    std::string unit_seq;        // one unit, native orientation
    bool backward = false;       // found on the path read in reverse: native orientation is sink -> source
    std::vector<ArrayElement> elements;                        // ordered emission plan
    std::vector<std::pair<std::size_t, std::size_t>> copy_ranges; // absolute [lo,hi) of copies (for removal)
};

// An array whose unit is itself a shorter unit repeated (a node holding two whole copies, UU, matched as one
// copy) is folded as that shorter unit walked once per copy, so both copies count. The shorter unit is the
// smallest whole repeat of the unit's period that is at least min_unit_bp and divides it. The period ignores
// letter case, as motifs do, so a soft-masked copy still counts; a piece whose case differs from the first
// keeps its own bases, so the pieces walked in order spell the same bases and nothing is lost.
TandemArray split_repeated_unit(TandemArray a, std::size_t min_unit_bp) {
    const std::string u = a.unit_seq;
    const std::size_t n = u.size();
    const auto same = [](char x, char y) {
        return std::toupper(static_cast<unsigned char>(x)) == std::toupper(static_cast<unsigned char>(y));
    };
    std::size_t period = n;
    for (std::size_t p = 1; p <= n / 2; ++p) {
        if (n % p != 0) continue;
        bool ok = true;
        for (std::size_t i = p; i < n && ok; ++i) ok = same(u[i], u[i - p]);
        if (ok) { period = p; break; }
    }
    std::size_t L = n;
    for (std::size_t m = period; m < n; m += period)
        if (m >= min_unit_bp && n % m == 0) { L = m; break; }
    if (L == n) return a;
    const std::size_t k = n / L;
    a.unit_seq = u.substr(0, L);
    a.copies *= k;
    std::vector<ArrayElement> el;
    for (const ArrayElement& e : a.elements) {
        if (!e.is_copy) { el.push_back(e); continue; }
        for (std::size_t c = 0; c < k; ++c) {
            ArrayElement piece{true, {}, u.substr(c * L, L)};
            if (piece.unit == a.unit_seq) piece.unit.clear();
            el.push_back(std::move(piece));
        }
    }
    a.elements = std::move(el);
    return a;
}

// A path read in the other direction: steps reversed, each flipped. It spells the reverse complement.
PathRecord reversed_path(const PathRecord& path) {
    PathRecord r = path;
    std::reverse(r.steps.begin(), r.steps.end());
    for (PathStep& st : r.steps) st.reverse = !st.reverse;
    return r;
}

// A crossing's interior steps read source -> sink: one walked the other way is reversed and each step flipped, so a
// haplotype stored reversed spells the same interior as one stored forward.
std::vector<PathStep> canonical_interior_steps(const PathRecord& path, const BubblePathInterval& used) {
    std::vector<PathStep> nat(path.steps.begin() + static_cast<std::ptrdiff_t>(used.left + 1),
                              path.steps.begin() + static_cast<std::ptrdiff_t>(used.right));
    if (!used.source_to_sink) {
        std::reverse(nat.begin(), nat.end());
        for (PathStep& st : nat) st.reverse = !st.reverse;
    }
    return nat;
}

// An array found on reversed_path(path) (N steps), placed back on path: positions mirrored, elements in
// the other order with their steps reversed and flipped, the unit as path reads it.
TandemArray mirror_array(const TandemArray& a, std::size_t N) {
    TandemArray m = a;
    m.start = N - (a.start + a.length);
    m.unit_seq = reverse_complement(a.unit_seq);
    m.backward = !a.backward;
    std::reverse(m.elements.begin(), m.elements.end());
    for (ArrayElement& el : m.elements) {
        std::reverse(el.steps.begin(), el.steps.end());
        for (PathStep& st : el.steps) st.reverse = !st.reverse;
        if (!el.unit.empty()) el.unit = reverse_complement(el.unit);
    }
    m.copy_ranges.clear();
    for (auto it = a.copy_ranges.rbegin(); it != a.copy_ranges.rend(); ++it)
        m.copy_ranges.emplace_back(N - it->second, N - it->first);
    return m;
}

// The identity of a repeat unit, independent of where the array is cut, its strand and how many copies it
// holds. Prevalence is counted per unit: only a repeat the haplotypes share, not any tandem in each,
// separates a population VNTR from unrelated arrays at one site.
std::string canonical_motif_key(const std::string& unit_in) {
    if (unit_in.empty()) return unit_in;
    // Upper case, so soft-masked and upper-case carriers of one repeat are one motif and count together.
    std::string unit = unit_in;
    for (char& c : unit) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    // Primitive root: if the unit is itself w repeated, the motif is w.
    std::string prim = unit;
    for (std::size_t p = 1; p <= unit.size() / 2; ++p) {
        if (unit.size() % p != 0) continue;
        bool ok = true;
        for (std::size_t i = p; i < unit.size() && ok; ++i) ok = (unit[i] == unit[i - p]);
        if (ok) { prim = unit.substr(0, p); break; }
    }
    // Least rotation on both strands, by Booth's algorithm: linear in the unit, since the exact-mode
    // single-copy search calls this on every candidate span.
    const auto least_rotation = [](const std::string& t) {
        const std::size_t n = t.size();
        if (n < 2) return t;
        const std::string d = t + t;
        std::vector<long> f(2 * n, -1);
        std::size_t k = 0;
        for (std::size_t j = 1; j < 2 * n; ++j) {
            long i = f[j - k - 1];
            while (i != -1 && d[j] != d[k + static_cast<std::size_t>(i) + 1]) {
                if (d[j] < d[k + static_cast<std::size_t>(i) + 1]) k = j - static_cast<std::size_t>(i) - 1;
                i = f[static_cast<std::size_t>(i)];
            }
            if (i == -1 && d[j] != d[k]) {
                if (d[j] < d[k]) k = j;
                f[j - k] = -1;
            } else {
                f[j - k] = i + 1;
            }
        }
        return d.substr(k, n);
    };
    const std::string f = least_rotation(prim);
    const std::string r = least_rotation(reverse_complement(prim));
    return f <= r ? f : r;
}


// Copies of a known motif lying on whole node boundaries, including a lone copy and one split across
// several nodes. detect_tandems needs two copies to find a period, so a single-copy haplotype has no array
// of its own; left literal while its neighbours fold, `call` would read it as CN 0. The motif is known, so
// spans of whole consecutive steps are matched against it. Steps in `skip` ([start, end) ranges, copies
// already found) are passed over, so a copy beside them is still found on its own.
std::vector<TandemArray> find_motif_copies_by_node(
    const Graph& graph,
    const PathRecord& path,
    std::size_t left,
    std::size_t right,
    const std::string& motif_key,
    std::size_t min_unit_bp,
    const std::vector<std::pair<std::size_t, std::size_t>>& skip) {

    std::vector<TandemArray> out;
    if (right <= left || right >= path.steps.size()) return out;

    const auto spell_span = [&](std::size_t i, std::size_t j) {
        std::string acc;
        for (std::size_t k = i; k < j; ++k) {
            const auto nit = graph.nodes.find(path.steps[k].node_id);
            if (nit == graph.nodes.end()) return std::string();
            acc += path.steps[k].reverse ? reverse_complement(nit->second.sequence)
                                         : nit->second.sequence;
        }
        return acc;
    };

    // Accepted spans, left to right. Both ends are step boundaries by construction, so no copy can be
    // accepted whose edge falls inside a node.
    constexpr std::size_t kMaxSpanSteps = 64;
    std::vector<std::pair<std::size_t, std::size_t>> spans;
    std::vector<std::string> span_seq;
    // A span matching the motif is whole copies of its primitive, so spans of other lengths are skipped
    // without computing a key. Spans grow one step at a time rather than being re-spelled from scratch.
    const std::size_t period = motif_key.size();
    for (std::size_t i = left + 1; i < right;) {
        std::size_t limit = right;   // a span ends before the next skipped range
        bool skipped = false;
        for (const auto& [lo, hi] : skip) {
            if (i >= lo && i < hi) { i = hi; skipped = true; break; }
            if (lo > i) limit = std::min(limit, lo);
        }
        if (skipped) continue;
        std::size_t end = 0;
        std::string seq_at_end;
        std::string seq;
        for (std::size_t j = i + 1; j <= limit && j - i <= kMaxSpanSteps; ++j) {
            seq += spell_span(j - 1, j);   // a zero-length step adds nothing
            if (seq.empty()) continue;
            if (seq.size() >= min_unit_bp && period > 0 && seq.size() % period == 0 &&
                canonical_motif_key(seq) == motif_key) {
                end = j;
                seq_at_end = seq;
                break;
            }
        }
        if (end == 0) { ++i; continue; }
        spans.emplace_back(i, end);
        span_seq.push_back(std::move(seq_at_end));
        i = end;
    }

    // Contiguous spans of the same sequence form one array; a gap, or a span spelling the motif
    // differently (its reverse complement, a rotation, two copies in one node), starts a new one: every
    // copy of an array is rewritten as the array's one unit, which must spell each copy exactly.
    for (std::size_t k = 0; k < spans.size();) {
        std::size_t j = k;
        while (j + 1 < spans.size() && spans[j + 1].first == spans[j].second && span_seq[j + 1] == span_seq[k]) ++j;
        TandemArray a;
        a.start = spans[k].first;
        a.length = spans[j].second - spans[k].first;
        a.copies = j - k + 1;
        a.unit_seq = span_seq[k];
        a.elements.assign(a.copies, ArrayElement{true, {}, {}});
        a.copy_ranges.emplace_back(spans[k].first, spans[j].second);
        out.push_back(std::move(a));
        k = j + 1;
    }
    return out;
}

// Copies of a site's repeat that every crossing walks just outside the site, as shared backbone, belong to
// its array: read from the interior alone, every haplotype's copy number would come out that many short. The
// site is widened over them: each boundary moves out past the copies whose steps every crossing shares, and
// the copies fold with the rest. `copy_steps(outward, left)` is given the shared steps read away from the
// boundary (the boundary first, each step as the site reads it) and returns how many of them hold the copy
// nearest the site, or 0. A site is widened only where nothing else can change: the steps taken belong to no
// other site (a side stops before the first copy holding another site's node), the new boundaries lie inside
// no other site, and every walk through the taken steps is a crossing of the widened site at the expected
// place. `sites_of_node` includes what earlier folded sites took by widening, so two neighbours never take
// the same nodes; `inside_of_node` holds the input interiors.
std::optional<Bubble> widen_over_shared_copies(
    const Graph& graph,
    const std::vector<BubblePathIndex>& path_indexes,
    const Bubble& bubble,
    const std::unordered_map<std::string, std::vector<std::size_t>>& sites_of_node,
    const std::unordered_map<std::string, std::vector<std::size_t>>& inside_of_node,
    const std::function<std::size_t(const std::vector<PathStep>&, bool)>& copy_steps) {

    constexpr std::size_t kMaxOutwardSteps = 512;
    constexpr std::size_t kMaxSharedCopies = 64;
    if (!bubble_orientation_holds(graph.paths, path_indexes, bubble)) return std::nullopt;
    std::vector<std::vector<BubblePathInterval>> crossings(graph.paths.size());
    std::size_t n_crossings = 0;
    for (std::size_t pi = 0; pi < graph.paths.size(); ++pi) {
        crossings[pi] = bubble_crossings(graph.paths[pi], path_indexes[pi], bubble, false);
        n_crossings += crossings[pi].size();
    }
    if (n_crossings == 0) return std::nullopt;

    // The steps every crossing shares outward from one boundary, as the site reads them.
    const auto shared_outward = [&](bool left) {
        std::vector<PathStep> shared;
        bool first = true;
        for (std::size_t pi = 0; pi < graph.paths.size(); ++pi) {
            const std::vector<PathStep>& steps = graph.paths[pi].steps;
            for (const BubblePathInterval& c : crossings[pi]) {
                // Outward along the path when the boundary is the crossing's right end, else backward.
                const bool at_right = left != c.source_to_sink;
                const std::size_t from = at_right ? c.right : c.left;
                std::vector<PathStep> out;
                for (std::size_t k = 0; k < kMaxOutwardSteps; ++k) {
                    if (at_right ? from + k >= steps.size() : k > from) break;
                    PathStep st = steps[at_right ? from + k : from - k];
                    if (!c.source_to_sink) st.reverse = !st.reverse;
                    out.push_back(st);
                }
                if (first) { shared = std::move(out); first = false; continue; }
                std::size_t m = 0;
                while (m < shared.size() && m < out.size() && shared[m].node_id == out[m].node_id &&
                       shared[m].reverse == out[m].reverse)
                    ++m;
                shared.resize(m);
                if (shared.size() < 2) return std::vector<PathStep>{};   // no copy and new boundary to share
            }
        }
        return shared;
    };

    std::unordered_set<std::string> taken;   // nodes the widened site takes in
    std::size_t taken_steps[2] = {0, 0};
    PathStep new_boundary[2];
    for (int side = 0; side < 2; ++side) {
        const bool left = side == 0;
        const std::vector<PathStep> shared = shared_outward(left);
        // A node of another site's ends the widening before the copy holding it: the copies past it are
        // that site's to take.
        const auto other_sites = [&](std::size_t i) {
            const auto it = sites_of_node.find(shared[i].node_id);
            if (it != sites_of_node.end())
                for (const std::size_t id : it->second)
                    if (id != bubble.id) return true;
            return false;
        };
        std::size_t k = 0;
        for (std::size_t copy = 0; copy < kMaxSharedCopies && k + 1 < shared.size(); ++copy) {
            const std::vector<PathStep> rest(shared.begin() + static_cast<std::ptrdiff_t>(k), shared.end());
            const std::size_t n = copy_steps(rest, left);
            if (n == 0 || k + n >= shared.size()) break;   // the new boundary must be shared too
            bool owned = false;
            for (std::size_t i = k; i < k + n && !owned; ++i) owned = other_sites(i);
            if (owned) break;
            k += n;
        }
        if (k == 0) continue;
        for (std::size_t i = 0; i < k; ++i) taken.insert(shared[i].node_id);
        const auto it = sites_of_node.find(shared[k].node_id);
        if (it != sites_of_node.end())
            for (const std::size_t id : it->second)
                if (id == bubble.id) return std::nullopt;   // the flank leads back into the site
        const auto in = inside_of_node.find(shared[k].node_id);
        if (in != inside_of_node.end())
            for (const std::size_t id : in->second)
                if (id != bubble.id) return std::nullopt;   // the new boundary is inside another site
        taken_steps[side] = k;
        new_boundary[side] = shared[k];
    }
    if (taken_steps[0] == 0 && taken_steps[1] == 0) return std::nullopt;

    Bubble w = bubble;
    if (taken_steps[0] > 0) { w.source = new_boundary[0].node_id; w.source_reverse = new_boundary[0].reverse; }
    if (taken_steps[1] > 0) { w.sink = new_boundary[1].node_id; w.sink_reverse = new_boundary[1].reverse; }
    if (w.source == w.sink || taken.count(w.source) || taken.count(w.sink)) return std::nullopt;
    for (const std::string& n : taken)
        if (std::find(w.inside.begin(), w.inside.end(), n) == w.inside.end()) w.inside.push_back(n);
    // Every crossing must come out as the same crossing, widened, and no walk may pass the nodes taken in
    // anywhere else.
    for (std::size_t pi = 0; pi < graph.paths.size(); ++pi) {
        const std::vector<BubblePathInterval> now = bubble_crossings(graph.paths[pi], path_indexes[pi], w, false);
        if (now.size() != crossings[pi].size()) return std::nullopt;
        for (std::size_t i = 0; i < now.size(); ++i) {
            const BubblePathInterval& c = crossings[pi][i];
            const std::size_t out_left = taken_steps[c.source_to_sink ? 0 : 1];
            const std::size_t out_right = taken_steps[c.source_to_sink ? 1 : 0];
            if (now[i].source_to_sink != c.source_to_sink || c.left < out_left ||
                now[i].left != c.left - out_left || now[i].right != c.right + out_right)
                return std::nullopt;
        }
        const std::vector<PathStep>& steps = graph.paths[pi].steps;
        for (std::size_t k = 0; k < steps.size(); ++k) {
            if (!taken.count(steps[k].node_id)) continue;
            bool within = false;
            for (const BubblePathInterval& c : now) within = within || (k > c.left && k < c.right);
            if (!within) return std::nullopt;
        }
    }
    return w;
}

// Non-overlapping tandem arrays within the bubble interval [left, right], comparing units by sequence
// token so copies with the same sequence match across distinct nodes. Short interruptions between copies
// (CGG (A) CGG CGG) are kept as literal steps.
std::vector<TandemArray> detect_tandems(
    const Graph& graph,
    const std::unordered_map<std::string, NodeTok>& node_tok,
    const PathRecord& path,
    std::size_t left,
    std::size_t right,
    std::size_t min_unit_bp,
    std::size_t min_copies,
    double max_interruption_frac) {

    std::vector<TandemArray> arrays;
    if (right < left || right >= path.steps.size()) {
        return arrays;
    }
    const std::size_t m = right - left + 1;

    std::vector<std::uint64_t> tok(m);
    std::vector<std::size_t> bp(m);
    for (std::size_t k = 0; k < m; ++k) {
        const PathStep& st = path.steps[left + k];
        const auto it = node_tok.find(st.node_id);
        const NodeTok nt = (it != node_tok.end()) ? it->second : NodeTok{};
        tok[k] = st.reverse ? nt.rc : nt.fwd;
        bp[k] = nt.len;
    }

    // Materialize an oriented step sequence on demand (only for accepted arrays).
    auto seq_at = [&](std::size_t rel) { return step_sequence(graph, path.steps[left + rel]); };

    std::unordered_map<std::uint64_t, std::vector<std::size_t>> token_positions;
    for (std::size_t k = 0; k < m; ++k) {
        token_positions[tok[k]].push_back(k);
    }

    // Fast equality by 64-bit sequence token during detection; accepted arrays are verified by sequence
    // below, and the sequence-preserving check backstops a token collision.
    auto block_equal = [&](std::size_t a, std::size_t b, std::size_t p) {
        for (std::size_t t = 0; t < p; ++t) {
            if (tok[a + t] != tok[b + t]) {
                return false;
            }
        }
        return true;
    };
    auto block_equal_seq = [&](std::size_t a, std::size_t b, std::size_t p) {
        for (std::size_t t = 0; t < p; ++t) {
            if (seq_at(a + t) != seq_at(b + t)) {
                return false;
            }
        }
        return true;
    };

    std::size_t lower_bound = 0; // do not extend an array left of the previous array

    // The array with unit [i, i+p) and its copies around i, or nothing if it fails a gate.
    auto try_period = [&](std::size_t i, std::size_t p) -> std::optional<TandemArray> {
        std::size_t unit_bp = 0;
        for (std::size_t t = 0; t < p; ++t) {
            unit_bp += bp[i + t];
        }

        // Collect copy start positions around the anchor, tolerating gaps.
        std::vector<std::size_t> copy_starts = {i};

        // Extend right.
        std::size_t last = i;
        while (true) {
            const std::size_t adj = last + p;
            if (adj + p <= m && block_equal(i, adj, p)) {
                copy_starts.push_back(adj);
                last = adj;
                continue;
            }
            bool found = false;
            std::size_t gap_bp = 0;
            for (std::size_t g = 1; g <= kMaxGapSteps; ++g) {
                const std::size_t cand = last + p + g;
                if (cand + p > m) {
                    break;
                }
                gap_bp += bp[last + p + (g - 1)];
                if (gap_bp > unit_bp) {
                    break;
                }
                if (block_equal(i, cand, p)) {
                    copy_starts.push_back(cand);
                    last = cand;
                    found = true;
                    break;
                }
            }
            if (!found) {
                break;
            }
        }

        // Extend left (not below lower_bound). Collected in reverse, prepended once.
        std::vector<std::size_t> left_starts;
        std::size_t first = i;
        while (first >= p && first - p >= lower_bound) {
            const std::size_t adj = first - p;
            if (block_equal(i, adj, p)) {
                left_starts.push_back(adj);
                first = adj;
                continue;
            }
            bool found = false;
            std::size_t gap_bp = 0;
            for (std::size_t g = 1; g <= kMaxGapSteps; ++g) {
                if (first < p + g || first - p - g < lower_bound) {
                    break;
                }
                const std::size_t cand = first - p - g;
                gap_bp += bp[first - g];
                if (gap_bp > unit_bp) {
                    break;
                }
                if (block_equal(i, cand, p)) {
                    left_starts.push_back(cand);
                    first = cand;
                    found = true;
                    break;
                }
            }
            if (!found) {
                break;
            }
        }
        if (!left_starts.empty()) {
            std::reverse(left_starts.begin(), left_starts.end());
            copy_starts.insert(copy_starts.begin(), left_starts.begin(), left_starts.end());
        }

        const std::size_t array_start = copy_starts.front();
        const std::size_t array_end = copy_starts.back() + p; // exclusive
        const std::size_t copies = copy_starts.size();

        // Build emission plan + measure interruptions.
        std::vector<ArrayElement> elements;
        std::vector<std::pair<std::size_t, std::size_t>> copy_ranges;
        std::size_t interruption_bp = 0;
        for (std::size_t k = 0; k < copy_starts.size(); ++k) {
            ArrayElement copy_el;
            copy_el.is_copy = true;
            elements.push_back(std::move(copy_el));
            copy_ranges.emplace_back(left + copy_starts[k], left + copy_starts[k] + p);
            if (k + 1 < copy_starts.size()) {
                const std::size_t gap_lo = copy_starts[k] + p;
                const std::size_t gap_hi = copy_starts[k + 1];
                if (gap_hi > gap_lo) {
                    ArrayElement gap_el;
                    gap_el.is_copy = false;
                    for (std::size_t s = gap_lo; s < gap_hi; ++s) {
                        gap_el.steps.push_back(path.steps[left + s]);
                        interruption_bp += bp[s];
                    }
                    elements.push_back(std::move(gap_el));
                }
            }
        }

        std::size_t array_bp = 0;
        for (std::size_t s = array_start; s < array_end; ++s) {
            array_bp += bp[s];
        }

        bool seq_ok = true;
        for (std::size_t k = 0; seq_ok && k < copy_starts.size(); ++k) {
            if (!block_equal_seq(i, copy_starts[k], p)) {
                seq_ok = false; // token collision: not a real identical copy
            }
        }

        const bool ok = seq_ok && copies >= min_copies && unit_bp >= min_unit_bp &&
                        (array_bp == 0 ||
                         static_cast<double>(interruption_bp) <= max_interruption_frac * static_cast<double>(array_bp));
        if (!ok) return std::nullopt;
        TandemArray arr;
        arr.start = left + array_start;
        arr.length = array_end - array_start;
        arr.copies = copies;
        arr.interruption_bp = interruption_bp;
        arr.elements = std::move(elements);
        arr.copy_ranges = std::move(copy_ranges);
        for (std::size_t t = 0; t < p; ++t) {
            arr.unit_seq += seq_at(i + t);
        }
        return arr;
    };

    std::size_t i = 0;
    while (i < m) {
        // Candidate periods at i, shortest first: every p with a clean adjacent pair. The shortest is not
        // always the unit (in a unit a,a,b the pair a,a gives p = 1), and moving on to i + 1 would find a
        // rotation and lose a copy, so each is tried until one yields an array that clears every gate. One
        // under --min-unit-bp is skipped before the costly extension.
        std::optional<TandemArray> found;
        const auto pit = token_positions.find(tok[i]);
        if (pit != token_positions.end()) {
            std::size_t tried = 0;
            for (const std::size_t j : pit->second) {
                if (j <= i) {
                    continue;
                }
                const std::size_t p = j - i;
                if (i + 2 * p > m) {
                    break;
                }
                if (++tried > 8) {
                    break;
                }
                if (!block_equal(i, i + p, p)) {
                    continue;
                }
                std::size_t unit_bp = 0;
                for (std::size_t t = 0; t < p; ++t) unit_bp += bp[i + t];
                if (unit_bp < min_unit_bp) {
                    continue;
                }
                found = try_period(i, p);
                if (found) break;
            }
        }
        if (!found) {
            ++i;
            continue;
        }
        const std::size_t array_end = (found->start - left) + found->length;
        arrays.push_back(std::move(*found));
        lower_bound = array_end;
        i = array_end;
    }
    return arrays;
}


// ---- unit graph (approximate mode, --min-similarity below 1) ----------------------------------------------
//
// Copies fold onto a unit graph that keeps their own differences: the unit is cut wherever a copy differs
// from it, stretches no copy changes are shared segment nodes, and each difference is a node of its own (a
// replacement, an insertion) or a skip over segments (a deletion). Each copy walks the segments and its own
// differences, so it spells exactly what it did, and copies sharing a difference share its node. Copies are
// counted at the anchor, a unit position no copy deletes, where every copy walks exactly one node: the
// longest segment no copy changes when there is one, else the position covered by the fewest nodes.

// unit [lo, hi) replaced by `alt` in one copy: a mismatch run, a deletion (alt empty) or an insertion
// (lo == hi), merged when they abut.
struct UnitEdit {
    std::size_t lo = 0, hi = 0;
    std::string alt;
};

// The edits that turn `unit` into `copy` (both in unit orientation), from a global alignment. A column the
// aligner calls a match but whose bytes differ (it ignores letter case) is an edit too: the segments spell
// the unit's bytes, and the copy must spell its own.
std::vector<UnitEdit> unit_edits(const std::string& unit, const std::string& copy) {
    std::string ops;
    nw_edit_ops(unit, copy, ops);
    std::vector<UnitEdit> out;
    std::size_t ti = 0, qi = 0;
    bool open = false;
    for (const char op : ops) {
        const bool same = op == 'M' && unit[ti] == copy[qi];
        if (same) {
            open = false;
            ++ti;
            ++qi;
            continue;
        }
        if (!open) {
            out.push_back(UnitEdit{ti, ti, std::string()});
            open = true;
        }
        UnitEdit& e = out.back();
        if (op == 'M' || op == 'X') { e.alt += copy[qi++]; ++ti; }
        else if (op == 'I') e.alt += copy[qi++];
        else ++ti;   // 'D'
        e.hi = ti;
    }
    // Normalise, so one difference is one node whichever copy carries it (the aligner places an indel
    // anywhere it fits in repeated sequence): each edit is trimmed of bases it shares with the unit at
    // either end, and a pure insertion or deletion is shifted as far left as the unit allows, as VCF
    // normalisation does, never into the copy's previous edit. The edits still spell the copy exactly.
    std::size_t floor_lo = 0;   // the previous edit's end
    std::vector<UnitEdit> norm;
    for (UnitEdit e : out) {
        while (e.lo < e.hi && !e.alt.empty() && unit[e.hi - 1] == e.alt.back()) { --e.hi; e.alt.pop_back(); }
        while (e.lo < e.hi && !e.alt.empty() && unit[e.lo] == e.alt.front()) { ++e.lo; e.alt.erase(0, 1); }
        if (e.lo == e.hi && e.alt.empty()) continue;   // nothing left: the columns spelled the unit after all
        if (e.alt.empty()) {                            // deletion of unit [lo, hi)
            while (e.lo > floor_lo && unit[e.lo - 1] == unit[e.hi - 1]) { --e.lo; --e.hi; }
        } else if (e.lo == e.hi) {                      // insertion of alt at lo
            while (e.lo > floor_lo && unit[e.lo - 1] == e.alt.back()) {
                e.alt.insert(e.alt.begin(), unit[e.lo - 1]);
                e.alt.pop_back();
                --e.lo;
                e.hi = e.lo;
            }
        }
        floor_lo = e.hi;
        norm.push_back(std::move(e));
    }
    return norm;
}

struct UnitGraph {
    std::string unit;
    std::vector<std::size_t> cuts;                                    // sorted, with 0 and |unit|
    std::unordered_map<std::string, std::vector<UnitEdit>> edits;     // copy (unit orientation) -> edits
    std::map<std::size_t, std::string> segment_node;                  // segment start -> node id
    std::map<std::string, std::string> edit_node;                     // "lo:hi:alt" -> node id
    std::size_t anchor_pos = SIZE_MAX;                                // unit position no copy deletes
    std::vector<std::string> nodes;                                   // every node created, in order
};

// Build the unit graph's shape from every distinct copy at a site; nodes are created as walks need them.
UnitGraph build_unit_graph(const std::string& unit, const std::vector<std::string>& copies, std::size_t threads) {
    UnitGraph g;
    g.unit = unit;
    std::vector<std::vector<UnitEdit>> per(copies.size());
    run_parallel(copies.size(), threads, [&](std::size_t i) { per[i] = unit_edits(unit, copies[i]); });
    std::set<std::size_t> cuts{0, unit.size()};
    std::vector<char> changed(unit.size(), 0);   // a unit base some copy replaces or deletes
    for (std::size_t i = 0; i < copies.size(); ++i) {
        for (const UnitEdit& e : per[i]) {
            cuts.insert(e.lo);
            cuts.insert(e.hi);
            for (std::size_t k = e.lo; k < e.hi; ++k) changed[k] = 1;
        }
        g.edits.emplace(copies[i], std::move(per[i]));
    }
    g.cuts.assign(cuts.begin(), cuts.end());
    std::size_t best = 0;
    for (std::size_t c = 0; c + 1 < g.cuts.size(); ++c) {
        const std::size_t lo = g.cuts[c], hi = g.cuts[c + 1];
        if (hi <= lo || changed[lo]) continue;   // a segment is wholly changed or wholly kept
        if (hi - lo > best) { best = hi - lo; g.anchor_pos = lo; }
    }
    if (g.anchor_pos == SIZE_MAX && !unit.empty()) {
        // No segment is kept by every copy. Per position: is it deleted by some copy, and how many
        // distinct replacements cover it (each one more node in the anchor set)?
        std::vector<char> deleted(unit.size(), 0);
        std::vector<std::set<std::string>> cover(unit.size());
        for (const auto& [copy, edits] : g.edits) {
            for (const UnitEdit& e : edits) {
                for (std::size_t k = e.lo; k < e.hi; ++k) {
                    if (e.alt.empty()) deleted[k] = 1;
                    else cover[k].insert(std::to_string(e.lo) + ":" + std::to_string(e.hi) + ":" + e.alt);
                }
            }
        }
        std::size_t fewest = SIZE_MAX;
        for (std::size_t k = 0; k < unit.size(); ++k) {
            if (!deleted[k] && cover[k].size() < fewest) { fewest = cover[k].size(); g.anchor_pos = k; }
        }
    }
    return g;
}

// The nodes covering the anchor position: the segment there and every replacement of a stretch holding
// it, among the nodes the walks created. A copy passes exactly one of them, so passes count copies.
std::vector<std::string> unit_graph_anchor_nodes(const UnitGraph& g) {
    std::vector<std::string> out;
    if (g.anchor_pos == SIZE_MAX) return out;
    const auto seg = std::prev(std::upper_bound(g.cuts.begin(), g.cuts.end(), g.anchor_pos));
    const auto s = g.segment_node.find(*seg);
    if (s != g.segment_node.end()) out.push_back(s->second);
    std::set<std::string> seen;
    for (const auto& [copy, edits] : g.edits) {
        for (const UnitEdit& e : edits) {
            if (e.alt.empty() || e.lo > g.anchor_pos || e.hi <= g.anchor_pos) continue;
            const auto n = g.edit_node.find(std::to_string(e.lo) + ":" + std::to_string(e.hi) + ":" + e.alt);
            if (n != g.edit_node.end() && seen.insert(n->second).second) out.push_back(n->second);
        }
    }
    std::sort(out.begin() + (s != g.segment_node.end() ? 1 : 0), out.end(), [](const std::string& a, const std::string& b) {
        return a.size() != b.size() ? a.size() < b.size() : a < b;
    });
    return out;
}

// The walk of one copy (unit orientation) through the unit graph, creating nodes as first needed.
std::vector<PathStep> unit_graph_walk(UnitGraph& g, GfaModel& model, const std::string& copy) {
    const std::vector<UnitEdit>& edits = g.edits.at(copy);
    std::vector<PathStep> walk;
    auto node = [&](std::string seq, auto& table, const auto& key) {
        const auto it = table.find(key);
        if (it != table.end()) return it->second;
        const std::string id = add_new_node(model, seq);
        table.emplace(key, id);
        g.nodes.push_back(id);
        return id;
    };
    std::size_t p = 0, ei = 0;
    const std::size_t L = g.unit.size();
    while (p < L || ei < edits.size()) {
        if (ei < edits.size() && edits[ei].lo == p) {
            const UnitEdit& e = edits[ei++];
            if (!e.alt.empty()) {
                const std::string key = std::to_string(e.lo) + ":" + std::to_string(e.hi) + ":" + e.alt;
                walk.push_back(PathStep{node(e.alt, g.edit_node, key), false});
            }
            p = e.hi;
            continue;
        }
        if (p >= L) break;
        const std::size_t next = *std::upper_bound(g.cuts.begin(), g.cuts.end(), p);
        walk.push_back(PathStep{node(g.unit.substr(p, next - p), g.segment_node, p), false});
        p = next;
    }
    return walk;
}

std::string spell_model_path(const GfaModel& model, const GfaPath& path) {
    std::string out;
    for (const PathStep& s : path.steps) {
        const auto it = model.seq.find(s.node_id);
        if (it == model.seq.end()) {
            throw std::runtime_error("Path references missing node: " + s.node_id);
        }
        out += s.reverse ? reverse_complement(it->second) : it->second;
    }
    return out;
}

// Ensure every consecutive step pair in every path has a backing L edge.
struct EdgeSet {
    std::unordered_set<std::string> keys;
    GfaModel& model;
    std::size_t added = 0;

    explicit EdgeSet(GfaModel& m) : model(m) {
        for (const GfaEdge& e : m.edges) {
            keys.insert(key(e.from, e.from_orient, e.to, e.to_orient));
        }
    }
    // One arc, read from either end: a+ -> b+ and b- -> a- are the same link, and a GFA may store either.
    // Keyed on the stored direction alone, a link recorded from one side would not be found from the other
    // (left behind when obsolete, or added twice by ensure()).
    static std::string key(const std::string& f, char fo, const std::string& t, char to) {
        const std::string fwd = f + fo + ">" + t + to;
        const char rf = fo == '+' ? '-' : '+';
        const char rt = to == '+' ? '-' : '+';
        const std::string rev = t + rt + ">" + f + rf;
        return fwd <= rev ? fwd : rev;
    }
    void ensure(const PathStep& a, const PathStep& b) {
        const std::string k = key(a.node_id, orient_char(a.reverse), b.node_id, orient_char(b.reverse));
        if (keys.insert(k).second) {
            model.edges.push_back(GfaEdge{a.node_id, orient_char(a.reverse), b.node_id, orient_char(b.reverse), "0M"});
            ++added;
        }
    }
};

struct PathArray {
    std::size_t start = 0;
    std::size_t length = 0;
    std::size_t bubble_id = 0;   // the site that produced this edit, for the acceptance check
    std::vector<PathStep> replacement;
};

// Approximate (similarity-based) folding: seed one repeat unit per bubble, find near-identical copies by
// banded alignment of the unit against each path, and fold them onto the site's unit graph, which keeps
// each copy's own differences. Detection is per path and multithreaded; the graph mutation is serial.

constexpr std::size_t kAnchorSeeds = 8;     // seeds spread across the unit for anchoring

std::string sample_of(const std::string& path_name) {
    const std::size_t h = path_name.find('#');
    return h == std::string::npos ? path_name : path_name.substr(0, h);
}

// Bases compared as the aligner compares them, without regard to letter case: a soft-masked copy is the
// same copy.
bool same_base(char a, char b) {
    return std::toupper(static_cast<unsigned char>(a)) == std::toupper(static_cast<unsigned char>(b));
}

int base2(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return -1;
    }
}

// Start positions where `kmer` (<=31 bp, ACGT) occurs in s, via a 2-bit roll.
std::vector<std::size_t> find_kmer_starts(const std::string& s, const std::string& kmer) {
    std::vector<std::size_t> out;
    const std::size_t k = kmer.size();
    if (k == 0 || k > 31 || s.size() < k) return out;
    std::uint64_t target = 0;
    for (const char c : kmer) {
        const int b = base2(c);
        if (b < 0) return out;
        target = (target << 2) | static_cast<std::uint64_t>(b);
    }
    const std::uint64_t mask = (k < 32) ? ((1ULL << (2 * k)) - 1) : ~0ULL;
    std::uint64_t code = 0;
    std::size_t valid = 0;
    for (std::size_t i = 0; i < s.size(); ++i) {
        const int b = base2(s[i]);
        if (b < 0) { valid = 0; code = 0; continue; }
        code = ((code << 2) | static_cast<std::uint64_t>(b)) & mask;
        if (++valid >= k && code == target) out.push_back(i - k + 1);
    }
    return out;
}

struct ApproxCopy {
    std::size_t bp_lo = 0;
    std::size_t bp_hi = 0;
    bool rev = false;
    double identity = 0.0;
};

// Near-identical copies of unit r (and its RC) in spelled sequence s: multi-seed anchors propose copy
// starts, then a banded global alignment of r into the window decides each copy, its extent and
// orientation. The band is (1-min_sim)*|r| + 8, uncapped, so divergent copies with a large indel still align.
std::vector<ApproxCopy> detect_copies_native(const std::string& s, const std::string& r, double min_sim) {
    std::vector<ApproxCopy> copies;
    const std::size_t k = 16;
    const std::size_t band =
        std::max<std::size_t>(8, static_cast<std::size_t>((1.0 - min_sim) * static_cast<double>(r.size())) + 8);
    // An interior a little shorter than the unit can still be one copy (one with a small deletion and
    // nothing else around it), so within the band the alignment and the identity threshold decide.
    if (r.size() < k || s.size() + band < r.size() || s.size() < k) return copies;
    const std::string rrc = reverse_complement(r);

    std::vector<std::size_t> anchors;
    auto add_seeds = [&](const std::string& ref) {
        const std::size_t span = ref.size() - k;
        for (std::size_t i = 0; i < kAnchorSeeds; ++i) {
            const std::size_t off = (kAnchorSeeds <= 1) ? 0 : span * i / (kAnchorSeeds - 1);
            for (const std::size_t pos : find_kmer_starts(s, ref.substr(off, k))) {
                if (pos >= off) anchors.push_back(pos - off);
                // A copy missing its first bases (a deletion at the unit start) implies a start before the
                // interior: within the band, it is tried at the interior's start instead of never at all.
                else if (off - pos <= band) anchors.push_back(0);
            }
        }
    };
    add_seeds(r);
    add_seeds(rrc);
    std::sort(anchors.begin(), anchors.end());
    anchors.erase(std::unique(anchors.begin(), anchors.end()), anchors.end());

    std::size_t last_end = 0;
    bool have_last = false;
    bool tried_last_end = false;   // a copy starting right where the previous one ended was tried
    auto try_at = [&](std::size_t st) {
        const std::size_t winlen = std::min(s.size() - st, r.size() + band);
        const std::string win = s.substr(st, winlen);
        const FitAlignResult ff = fit_align(r, win, band);
        const FitAlignResult fr = fit_align(rrc, win, band);
        const bool use_rev = fr.identity > ff.identity;
        const FitAlignResult& best = use_rev ? fr : ff;
        if (st == last_end) tried_last_end = true;
        if (!(best.ok && best.identity >= min_sim && best.target_consumed >= r.size() / 2)) return false;
        copies.push_back({st, st + best.target_consumed, use_rev, best.identity});
        last_end = st + best.target_consumed;
        have_last = true;
        tried_last_end = false;
        return true;
    };
    // An anchor is where the copy would start if it matched the unit up to its seed, so a copy with an indel
    // before every seed is anchored a few bases off the previous copy's end; within the band of that end, it
    // is tried at the end. A copy missing its last bases lets the free end gap run its alignment into the
    // next copy, which then starts short and can fail the identity test: when the copy at the previous end
    // fails, it is tried where its own seeds place it, and the previous copy is trimmed back to that start
    // if it still passes there.
    auto try_overlapping = [&](std::size_t st) {
        ApproxCopy& prev = copies.back();
        if (st <= prev.bp_lo + r.size() / 2) return false;
        const std::size_t winlen = std::min(s.size() - st, r.size() + band);
        const std::string win = s.substr(st, winlen);
        const FitAlignResult ff = fit_align(r, win, band);
        const FitAlignResult fr = fit_align(rrc, win, band);
        const bool use_rev = fr.identity > ff.identity;
        const FitAlignResult& best = use_rev ? fr : ff;
        if (!(best.ok && best.identity >= min_sim && best.target_consumed >= r.size() / 2)) return false;
        const FitAlignResult pt = fit_align(prev.rev ? rrc : r, s.substr(prev.bp_lo, st - prev.bp_lo), band);
        if (!(pt.ok && pt.identity >= min_sim && pt.target_consumed >= r.size() / 2)) return false;
        prev.bp_hi = prev.bp_lo + pt.target_consumed;
        prev.identity = pt.identity;
        copies.push_back({st, st + best.target_consumed, use_rev, best.identity});
        last_end = st + best.target_consumed;
        tried_last_end = false;
        return true;
    };
    for (const std::size_t st : anchors) {
        if (have_last && st < last_end) {
            if (st + band >= last_end && !tried_last_end && last_end < s.size() && !try_at(last_end)) try_overlapping(st);
            continue;
        }
        if (try_at(st)) continue;
        if (have_last && st <= last_end + band && !tried_last_end && last_end < s.size()) try_at(last_end);
    }
    return copies;
}


// Base-level tandem seeding: the smallest lag (>= min_unit_bp) at which the spelled interval is
// self-similar. Fallback for graphs whose node boundaries do not follow repeat-unit boundaries, where
// `detect_tandems` -- which measures its period in node steps -- sees no period at all. Cost is bounded
// and independent of unit size: 2-bit k-mer scans propose candidate lags, each settled by a sampled
// identity check rather than a full O(len^2) scan.
std::string seed_unit_by_period(const std::string& s, std::size_t min_unit_bp, double min_identity) {
    constexpr std::size_t k = 16;
    constexpr std::size_t kProbes = 16;        // seed offsets spread across the interval: the array is
                                               // usually interior (the interval carries the source/sink
                                               // flanks, often kb of backbone), and probing only the
                                               // start would seed inside the flank and find nothing
    constexpr std::size_t kMaxCandidates = 8;  // candidate lags examined per offset
    constexpr std::size_t kSamples = 2048;     // sampled identity keeps each candidate ~O(1)
    if (min_unit_bp == 0 || s.size() < 2 * min_unit_bp || s.size() < 2 * k) {
        return {};
    }
    const std::size_t stride = std::max<std::size_t>(k, s.size() / kProbes);
    for (std::size_t off = 0, tried_off = 0; tried_off < kProbes && off + k <= s.size();
         off += stride, ++tried_off) {
        const std::vector<std::size_t> occ = find_kmer_starts(s, s.substr(off, k));
        std::size_t tried = 0;
        for (const std::size_t j : occ) {
            if (j <= off) {
                continue;
            }
            const std::size_t p = j - off;
            if (p < min_unit_bp) {
                continue;
            }
            if (off + 2 * p > s.size()) {
                break; // need two full copies to call it a tandem
            }
            if (++tried > kMaxCandidates) {
                break;
            }
            // Verify locally, on adjacent copy pairs walking right from the seed: measured across the whole
            // interval, the non-periodic source/sink flanks would make a real array's lag score near-random.
            const std::size_t sample_step = std::max<std::size_t>(1, p / kSamples);
            std::size_t copies = 1;
            for (std::size_t a = off; a + 2 * p <= s.size(); a += p) {
                std::size_t seen = 0;
                std::size_t hit = 0;
                for (std::size_t i = 0; i < p; i += sample_step) {
                    ++seen;
                    if (same_base(s[a + i], s[a + p + i])) {
                        ++hit;
                    }
                }
                if (seen == 0 || static_cast<double>(hit) / static_cast<double>(seen) < min_identity) {
                    break;
                }
                ++copies;
            }
            if (copies >= 2) {
                return s.substr(off, p); // a rotation of the unit; any rotation folds equivalently
            }
        }
    }
    return {};
}

// Whether `longer` is only copies of `shorter`, as folding on `shorter` would read it: a whole multiple
// (k >= 2 copies, its length within --min-similarity of k units), and, taken as a cycle so that any
// rotation counts, covered by k adjacent copies the copy detector finds, each in either orientation,
// with their edits and the bases left between them within --min-similarity of its length. Copies
// alternating two variants of one unit (U V U V) show an exact period of two copies, and folding on the
// dimer UV would count half of them. A module holding a sub-repeat has unique sequence besides it and does
// not fit.
bool near_multiple(std::string longer, std::string shorter, double min_similarity) {
    const std::size_t l = shorter.size(), n = longer.size();
    if (l == 0 || n < l + l / 2) return false;
    const std::size_t k = (n + l / 2) / l;
    const double budget = (1.0 - min_similarity) * static_cast<double>(n);
    const std::size_t whole = k * l;
    if (k < 2 || static_cast<double>(whole > n ? whole - n : n - whole) > budget) return false;
    for (char& c : longer) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    for (char& c : shorter) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    // Copies over the cycle: `longer` twice holds every rotation of it.
    const std::vector<ApproxCopy> copies = detect_copies_native(longer + longer, shorter, min_similarity);
    for (std::size_t i = 0; i < copies.size() && copies[i].bp_lo < n; ++i) {
        const std::size_t lo = copies[i].bp_lo, hi = lo + n;
        double cost = 0.0;
        std::size_t covered = 0, count = 0;
        for (std::size_t j = i; j < copies.size() && copies[j].bp_lo < hi; ++j) {
            const std::size_t end = std::min(copies[j].bp_hi, hi);
            const std::size_t span = end - copies[j].bp_lo;
            covered += span;
            cost += (1.0 - copies[j].identity) * static_cast<double>(std::max(span, l));
            ++count;
        }
        cost += static_cast<double>(n - std::min(covered, n));
        if (count >= 2 && cost <= budget) return true;
    }
    return false;
}

std::string shortest_period(const std::string& unit, const PanphorteOptions& options);

// Candidate repeat units for a bubble, best first, seeded by the exact tandem detector: a unit comes from a
// clean adjacent identical pair, so it is the true period (a whole module, not a recurring sub-segment).
// The per-path scan is parallel; the tally merges serially in path order and ties go to the longer unit,
// so the result is deterministic across threads.
std::vector<std::string> pick_reference_units(
    const Graph& graph,
    const std::unordered_map<std::string, NodeTok>& node_tok,
    const std::vector<BubblePathIndex>& path_indexes,
    const Bubble& bubble,
    const PanphorteOptions& options) {

    // A path touching the site from outside is not a crossing where the table's orientations hold
    // (bubble_orientation_holds); the positional fallback is only for a table whose orientations do not.
    const bool crossing_fallback = !bubble_orientation_holds(graph.paths, path_indexes, bubble);
    using Cand = std::pair<std::string, std::size_t>; // (unit_seq, copies)
    std::vector<std::vector<Cand>> per_path(graph.paths.size());
    run_parallel(graph.paths.size(), options.threads, [&](std::size_t pi) {
        // The interior only, as detection reads it: a periodic backbone anchor must not win the seed, since
        // anchors are never rewritten. Read source -> sink (a crossing walked the other way on the reversed
        // path), so the seed has the orientation detection matches copies in, however haplotypes are stored.
        std::vector<Cand> cands;
        const std::size_t N = graph.paths[pi].steps.size();
        std::optional<PathRecord> rev_path;
        for (const BubblePathInterval& used : bubble_crossings(graph.paths[pi], path_indexes[pi], bubble, crossing_fallback)) {
            if (used.inside_count == 0 || used.right < used.left + 2) continue;
            if (!used.source_to_sink && !rev_path) rev_path = reversed_path(graph.paths[pi]);
            // Exact tandem arrays only need an adjacent identical pair (min_copies 2),
            // independent of the approximate --min-copies the caller will enforce.
            const auto arrays = used.source_to_sink
                ? detect_tandems(graph, node_tok, graph.paths[pi], used.left + 1, used.right - 1,
                                 options.min_unit_bp, 2, options.max_interruption_frac)
                : detect_tandems(graph, node_tok, *rev_path, N - used.right, N - 2 - used.left,
                                 options.min_unit_bp, 2, options.max_interruption_frac);
            for (const TandemArray& arr : arrays) {
                if (arr.unit_seq.size() >= options.min_unit_bp) {
                    cands.emplace_back(arr.unit_seq, arr.copies);
                }
            }
        }
        per_path[pi] = std::move(cands);
    });
    // Votes: one per haplotype per canonical motif, so a sub-repeat cannot outvote the module containing
    // it by its copy count. Within a motif, the unit sequence (its rotation) is the copy-weighted favourite.
    std::map<std::string, std::set<std::size_t>> motif_haps;               // motif -> haplotypes
    std::map<std::string, std::map<std::string, std::size_t>> seq_weight;  // motif -> unit -> copies
    auto tally_paths = [&](const std::vector<std::vector<Cand>>& per) {
        for (std::size_t pi = 0; pi < per.size(); ++pi) {
            for (const auto& [seq, w] : per[pi]) {
                const std::string key = canonical_motif_key(seq);
                motif_haps[key].insert(pi);
                seq_weight[key][seq] += w;
            }
        }
    };
    tally_paths(per_path);
    // Fallback when no node-step period is found anywhere in the bubble: retry at base level, where a unit
    // straddling node boundaries is still visible. Graphs split at repeat boundaries never reach it.
    if (motif_haps.empty()) {
        std::vector<std::vector<Cand>> per_path_bp(graph.paths.size());
        run_parallel(graph.paths.size(), options.threads, [&](std::size_t pi) {
            for (const BubblePathInterval& used : bubble_crossings(graph.paths[pi], path_indexes[pi], bubble, crossing_fallback)) {
                if (used.inside_count == 0) continue;
                if (used.right <= used.left + 1) continue;       // interior only, as detection uses
                const std::string snat = spell_path_steps_sequence(graph, canonical_interior_steps(graph.paths[pi], used));
                const std::string unit =
                    seed_unit_by_period(snat, options.min_unit_bp, options.min_similarity);
                if (!unit.empty()) {
                    per_path_bp[pi].emplace_back(unit, snat.size() / unit.size());
                }
            }
        });
        tally_paths(per_path_bp);
    }
    // The unit standing for a motif: its copy-weighted favourite sequence.
    auto rep_of = [&](const std::string& motif) {
        std::string best;
        std::size_t best_w = 0;
        for (const auto& [seq, w] : seq_weight[motif]) {
            const bool better = w > best_w || (w == best_w && seq.size() > best.size()) ||
                                (w == best_w && seq.size() == best.size() && seq < best);
            if (best.empty() || better) { best_w = w; best = seq; }
        }
        return best;
    };
    // A motif that is only copies of a shorter one (near_multiple) is that motif: its haplotypes vote for
    // the shorter one, the most primitive it is a multiple of. Longest first, so a k-mer of a k-mer folds
    // all the way down. The shorter one must have support of its own, at least --min-array-prevalence of
    // the longer one's, as the container rule below asks of a module: otherwise one haplotype with two
    // adjacent copies of a module's sub-unit would move the whole population's module votes onto it.
    {
        std::vector<std::string> by_len;
        for (const auto& [motif, haps] : motif_haps) by_len.push_back(motif);
        std::sort(by_len.begin(), by_len.end(), [&](const std::string& a, const std::string& b) {
            const std::size_t la = rep_of(a).size(), lb = rep_of(b).size();
            return la != lb ? la > lb : a < b;
        });
        for (const std::string& a : by_len) {
            const std::string ua = rep_of(a);
            std::string into;
            for (const std::string& b : by_len) {
                if (b == a || !motif_haps.count(b)) continue;
                const std::string ub = rep_of(b);
                if (ub.size() >= ua.size()) continue;
                if (static_cast<double>(motif_haps[b].size()) <
                    options.min_array_prevalence * static_cast<double>(motif_haps[a].size()))
                    continue;
                if (near_multiple(ua, ub, options.min_similarity) && (into.empty() || ub.size() < rep_of(into).size()))
                    into = b;
            }
            if (into.empty()) continue;
            motif_haps[into].insert(motif_haps[a].begin(), motif_haps[a].end());
            motif_haps.erase(a);
        }
    }
    // Near-identical motifs are variants of one repeat; keyed by exact sequence alone they would split its
    // votes. Motifs whose units align end to end (as cycles, either strand) within --min-similarity are
    // pooled: the group votes with all its haplotypes and is represented by its medoid, the member closest
    // to the others weighted by their haplotypes. Pairwise, so only the kMaxPooled best-supported motifs.
    std::map<std::string, std::vector<std::string>> pooled;   // group's motif -> member units, closest first
    // The seed the vote gives without pooling: always a candidate, and the first, so pooling changes a site's
    // seed only where another unit detects strictly more copies.
    std::string unpooled;
    {
        std::size_t v = 0;
        for (const auto& [motif, haps] : motif_haps) {
            const std::string r = rep_of(motif);
            if (unpooled.empty() || haps.size() > v || (haps.size() == v && r.size() > unpooled.size())) {
                unpooled = r;
                v = haps.size();
            }
        }
    }
    {
        constexpr std::size_t kMaxPooled = 64;
        std::vector<std::string> ms;
        for (const auto& [motif, haps] : motif_haps) ms.push_back(motif);
        std::stable_sort(ms.begin(), ms.end(), [&](const std::string& a, const std::string& b) {
            return motif_haps[a].size() > motif_haps[b].size();
        });
        if (ms.size() > kMaxPooled) ms.resize(kMaxPooled);
        const std::size_t m = ms.size();
        std::vector<std::string> unit(m);
        for (std::size_t i = 0; i < m; ++i) unit[i] = rep_of(ms[i]);
        // Identity of a against b read as a cycle, or 0: one copy of a covering all of b within the budget.
        const auto near_identical = [&](const std::string& a, const std::string& b) {
            const double budget = (1.0 - options.min_similarity) * static_cast<double>(std::max(a.size(), b.size()));
            const double dl = a.size() > b.size() ? double(a.size() - b.size()) : double(b.size() - a.size());
            if (dl > budget) return 0.0;
            std::string ua = a, ub = b;
            for (char& c : ua) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
            for (char& c : ub) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
            for (const ApproxCopy& c : detect_copies_native(ub + ub, ua, options.min_similarity)) {
                if (c.bp_lo >= ub.size()) break;
                const double span = static_cast<double>(c.bp_hi - c.bp_lo);
                if (std::abs(span - static_cast<double>(ub.size())) <= budget) return c.identity;
            }
            return 0.0;
        };
        std::vector<std::vector<double>> sim(m, std::vector<double>(m, 0.0));
        std::vector<std::size_t> parent(m);
        for (std::size_t i = 0; i < m; ++i) parent[i] = i;
        const auto find = [&](std::size_t x) { while (parent[x] != x) x = parent[x] = parent[parent[x]]; return x; };
        for (std::size_t i = 0; i < m; ++i)
            for (std::size_t j = i + 1; j < m; ++j) {
                sim[i][j] = sim[j][i] = near_identical(unit[i], unit[j]);
                if (sim[i][j] > 0.0) parent[find(i)] = find(j);
            }
        std::map<std::size_t, std::vector<std::size_t>> groups;
        for (std::size_t i = 0; i < m; ++i) groups[find(i)].push_back(i);
        for (const auto& [root, members] : groups) {
            if (members.size() < 2) continue;
            std::size_t best = members.front();
            double best_score = -1.0;
            for (const std::size_t i : members) {    // members are in support order: ties keep the better supported
                double score = static_cast<double>(motif_haps[ms[i]].size());
                for (const std::size_t j : members)
                    if (j != i) score += sim[i][j] * static_cast<double>(motif_haps[ms[j]].size());
                if (score > best_score) { best_score = score; best = i; }
            }
            std::vector<std::pair<double, std::size_t>> order;
            for (const std::size_t i : members) {
                double score = static_cast<double>(motif_haps[ms[i]].size());
                for (const std::size_t j : members)
                    if (j != i) score += sim[i][j] * static_cast<double>(motif_haps[ms[j]].size());
                order.emplace_back(-score, i);
            }
            std::stable_sort(order.begin(), order.end(),
                             [](const auto& a, const auto& b) { return a.first < b.first; });
            for (const auto& [score, i] : order) pooled[ms[best]].push_back(unit[i]);
            for (const std::size_t i : members) {
                if (i == best) continue;
                motif_haps[ms[best]].insert(motif_haps[ms[i]].begin(), motif_haps[ms[i]].end());
                motif_haps.erase(ms[i]);
            }
        }
    }
    std::map<std::string, std::size_t> motif_votes;                        // motif -> haplotypes
    for (const auto& [motif, haps] : motif_haps) motif_votes[motif] = haps.size();
    std::string win_motif;
    std::size_t win_votes = 0;
    for (const auto& [motif, votes] : motif_votes) {
        const bool better = votes > win_votes ||
                            (votes == win_votes && rep_of(motif).size() > rep_of(win_motif).size());
        if (win_motif.empty() || better) { win_motif = motif; win_votes = votes; }
    }
    if (win_motif.empty()) return {};
    std::string best = rep_of(win_motif);
    // A pooled group's members are all candidates: each is cut at its own rotation, and one off the copies'
    // phase loses a copy at every array edge. The caller tests them.
    std::vector<std::string> candidates{best};
    // The unpooled winner leads only when it belongs to the winning group, so an unrelated repeat cannot
    // outvote the group by copy count in the caller.
    if (const auto it = pooled.find(win_motif); it != pooled.end()) {
        const bool member = std::find(it->second.begin(), it->second.end(), unpooled) != it->second.end();
        candidates = member ? std::vector<std::string>{unpooled} : std::vector<std::string>{};
        for (const std::string& u : it->second)
            if (u != unpooled || !member) candidates.push_back(u);
    }
    constexpr std::size_t kMaxCandidates = 5;   // each costs one phase test
    if (candidates.size() > kMaxCandidates) candidates.resize(kMaxCandidates);
    // A longer unit that contains the winner is the module the winner is a sub-repeat of: prefer it when it
    // is itself well supported (at least --min-array-prevalence of the winner's haplotypes), so a lone
    // outlier duplicating a larger region cannot take the seed from a population VNTR.
    {
        const std::string w_rc = reverse_complement(best);
        std::string container;
        for (const auto& [motif, votes] : motif_votes) {
            if (motif == win_motif) continue;
            if (static_cast<double>(votes) < options.min_array_prevalence * static_cast<double>(win_votes)) continue;
            const std::string rep = rep_of(motif);
            if (rep.size() <= best.size()) continue;
            const std::string twice = rep + rep;
            if (twice.find(best) == std::string::npos && twice.find(w_rc) == std::string::npos) continue;
            if (rep.size() > container.size()) container = rep;
        }
        if (!container.empty()) candidates = {container};
    }
    std::vector<std::string> out;
    for (std::string unit : candidates) {
        unit = shortest_period(unit, options);
        if (std::find(out.begin(), out.end(), unit) == out.end()) out.push_back(std::move(unit));
    }
    return out;
}

std::string shortest_period(const std::string& best, const PanphorteOptions& options) {
    // The unit is the shortest period it is made of, as exact mode's motif is its primitive root: a node
    // holding two copies (UU), or copies alternating two variants (U V U V), shows a period of two copies,
    // and folding on it would count half of them. Exact roots first, then near ones within --min-similarity
    // (near_multiple), up to kMaxFold copies per seeded unit and never below --min-unit-bp.
    {
        constexpr std::size_t kMaxFold = 8;
        const std::size_t L = best.size();
        for (std::size_t p = 1; p <= L / 2; ++p) {
            if (L % p != 0) continue;
            bool ok = true;
            for (std::size_t i = p; i < L && ok; ++i) ok = same_base(best[i], best[i - p]);
            // A root below --min-unit-bp (a microsatellite) is not a unit this site folds on: the seed stands.
            if (ok) return p >= options.min_unit_bp ? best.substr(0, p) : best;
        }
        for (std::size_t k = std::min(kMaxFold, L / std::max<std::size_t>(1, options.min_unit_bp)); k >= 2; --k) {
            const std::string sub = best.substr(0, (L + k / 2) / k);
            if (sub.size() >= options.min_unit_bp && near_multiple(best, sub, options.min_similarity)) return sub;
        }
    }
    return best;
}

// The phase of an approximate unit. A seed is cut wherever its probe landed, so it is a rotation of the
// copies: mid-array any rotation folds the same copies, but at the edges it does not. Rotated d bases into a
// copy, k copies become a d-base head fragment, k-1 unit steps and a tail fragment: one copy lost per
// haplotype, and a single-copy haplotype gets no unit at all.
//
// The rotation is proposed from the data rather than searched: on the commonest distinct interiors, the
// bases just before the first copy are matched against the unit's tail and those just after the last copy
// against its head (a cumulative Hamming scan, linear in the unit); the matched length is the phase offset.
// It is kept only if it detects strictly more copies on the same interiors, so a site that already folds
// well is unchanged. Cost: at most kPhaseScan + kPhaseSample single-interior detections per bubble.
struct UnitPhase {
    std::string unit;
    std::size_t shift = 0;          // bases the seeded unit was rotated left by (0: unchanged)
    std::size_t copies_seeded = 0;  // copies detected on the sample (weighted by haplotypes) with the seeded unit
    std::size_t copies_phased = 0;  // ... with the proposed rotation (= copies_seeded when nothing was proposed)
};
UnitPhase phase_unit_at_array_edge(const Graph& graph, const std::vector<BubblePathIndex>& path_indexes,
                                   const Bubble& bubble, const std::string& unit, const PanphorteOptions& options) {
    constexpr std::size_t kPhaseSample = 8;     // distinct interiors that propose and test the rotation
    constexpr std::size_t kMinPartial = 8;     // a partial copy shorter than this is not evidence of phase
    UnitPhase out;
    out.unit = unit;
    const std::size_t L = unit.size();
    // A path touching the site from outside is not a crossing where the table's orientations hold
    // (bubble_orientation_holds); the positional fallback is only for a table whose orientations do not.
    const bool crossing_fallback = !bubble_orientation_holds(graph.paths, path_indexes, bubble);
    std::map<std::string, std::size_t> interiors;   // spelled interior -> haplotypes
    for (std::size_t pi = 0; pi < graph.paths.size(); ++pi) {
        for (const BubblePathInterval& used : bubble_crossings(graph.paths[pi], path_indexes[pi], bubble, crossing_fallback)) {
            if (used.inside_count == 0) continue;
            if (used.right <= used.left + 1) continue;
            // Source -> sink, as detection reads it: the rotation is proposed and tested on what will be folded.
            ++interiors[spell_path_steps_sequence(graph, canonical_interior_steps(graph.paths[pi], used))];
        }
    }
    // The sample: the commonest interiors in which the seeded unit finds a copy (one without the repeat says
    // nothing about its phase). Ranked by haplotypes, then length, then sequence (deterministic), and scanned
    // in batches until kPhaseSample carry a copy or kPhaseScan were tried.
    constexpr std::size_t kPhaseScan = 64;
    std::vector<std::pair<const std::string*, std::size_t>> ranked;
    for (const auto& [s, n] : interiors) ranked.emplace_back(&s, n);
    std::stable_sort(ranked.begin(), ranked.end(), [](const auto& a, const auto& b) {
        if (a.second != b.second) return a.second > b.second;
        return a.first->size() > b.first->size();
    });
    if (ranked.size() > kPhaseScan) ranked.resize(kPhaseScan);
    std::vector<std::pair<const std::string*, std::size_t>> sample;
    std::vector<std::vector<ApproxCopy>> seeded;
    for (std::size_t from = 0; from < ranked.size() && sample.size() < kPhaseSample; from += kPhaseSample) {
        const std::size_t to = std::min(ranked.size(), from + kPhaseSample);
        std::vector<std::vector<ApproxCopy>> batch(to - from);
        run_parallel(batch.size(), options.threads, [&](std::size_t i) {
            batch[i] = detect_copies_native(*ranked[from + i].first, unit, options.min_similarity);
        });
        for (std::size_t i = 0; i < batch.size() && sample.size() < kPhaseSample; ++i) {
            if (batch[i].empty()) continue;
            sample.push_back(ranked[from + i]);
            seeded.push_back(std::move(batch[i]));
        }
    }
    // Counted for every unit, so candidate seeds can be compared; a unit too short for a partial copy of
    // kMinPartial bases at each edge is not re-phased.
    if (L < 2 * kMinPartial) {
        for (std::size_t i = 0; i < sample.size(); ++i) out.copies_seeded += seeded[i].size() * sample[i].second;
        out.copies_phased = out.copies_seeded;
        return out;
    }
    std::map<std::size_t, std::size_t> votes;   // proposed left rotation -> haplotypes
    for (std::size_t i = 0; i < sample.size(); ++i) {
        out.copies_seeded += seeded[i].size() * sample[i].second;
        if (seeded[i].empty() || seeded[i].front().rev) continue;
        const std::string& s = *sample[i].first;
        const std::size_t c0 = seeded[i].front().bp_lo;
        std::size_t mism = 0, matched = 0;
        for (std::size_t j = 1; j <= std::min(c0, L - 1); ++j) {
            const bool same = same_base(s[c0 - j], unit[L - j]);
            if (!same) ++mism;
            // The partial copy must begin on a matching base: extending it over a mismatch would pull a flanking
            // base into the first copy (a spurious difference, or lost under --lossy-units).
            if (same && j >= kMinPartial &&
                static_cast<double>(j - mism) / static_cast<double>(j) >= options.min_similarity) {
                matched = j;
            }
        }
        if (matched > 0) votes[L - matched] += sample[i].second;
        // The other edge: the bases after the last copy against the unit's head. A copy split by the rotation is
        // a head partial plus a tail partial, so whichever is long enough proposes the same rotation; the partial
        // must end on a matching base, for the same reason as above.
        const ApproxCopy& last = seeded[i].back();
        if (last.rev) continue;
        const std::size_t e = last.bp_hi;
        mism = 0;
        matched = 0;
        for (std::size_t j = 1; e + j <= s.size() && j <= L - 1; ++j) {
            const bool same = same_base(s[e + j - 1], unit[j - 1]);
            if (!same) ++mism;
            if (same && j >= kMinPartial &&
                static_cast<double>(j - mism) / static_cast<double>(j) >= options.min_similarity) {
                matched = j;
            }
        }
        if (matched > 0) votes[matched] += sample[i].second;
    }
    out.copies_phased = out.copies_seeded;
    std::size_t shift = 0, best_votes = 0;
    for (const auto& [d, n] : votes) {
        if (n > best_votes) { best_votes = n; shift = d; }
    }
    if (shift == 0) return out;
    const std::string rotated = unit.substr(shift) + unit.substr(0, shift);
    std::vector<std::size_t> phased(sample.size(), 0);
    run_parallel(sample.size(), options.threads, [&](std::size_t i) {
        phased[i] = detect_copies_native(*sample[i].first, rotated, options.min_similarity).size();
    });
    std::size_t total = 0;
    for (std::size_t i = 0; i < sample.size(); ++i) total += phased[i] * sample[i].second;
    if (total > out.copies_seeded) {
        out.unit = rotated;
        out.shift = shift;
        out.copies_phased = total;
    }
    return out;
}

} // namespace

namespace {

// Thrown when a folded site's original route is still walkable: other walks -- a contig that starts or
// ends inside the site, a haplotype whose partial copies sit on the same nodes -- keep every node and
// link of the replaced span. The site would be in the graph twice, folded and unfolded. Caught by
// panphorte_normalize, which declines those sites and runs again.
struct SharedRouteSites : std::runtime_error {
    std::set<std::size_t> bubble_ids;
    explicit SharedRouteSites(std::set<std::size_t> ids, const std::string& what)
        : std::runtime_error(what), bubble_ids(std::move(ids)) {}
};

void normalize_pass(const PanphorteOptions& options, PanphorteSummary* summary_out,
                    const std::set<std::size_t>& declined) {
    if (options.gfa_path.empty()) {
        throw std::runtime_error("panphorte requires -i/--gfa <graph.gfa>");
    }
    if (options.out_prefix.empty()) {
        throw std::runtime_error("panphorte requires -o/--out-prefix <prefix>");
    }
    if (options.refine && options.reference_path.empty())
        throw std::runtime_error("panphorte: --refine needs --reference-path (sites are aligned in reference "
                                 "orientation and the reference must cross each one)");

    cli::require_regular_input("panphorte", "-i", options.gfa_path);   // read here and again as a model
    ParseGfaOptions parse_opts;
    parse_opts.include_paths = true;
    parse_opts.include_sequences = true;
    const Graph graph = parse_gfa(options.gfa_path, parse_opts);
    // The same contract every other module applies: panphorte rewrites the graph, so a malformed input
    // (a path step with no link behind it) would come out looking repaired without ever being validated.
    validate_graph_paths(graph, "panphorte", true, true);

    GfaModel model = read_gfa_model(options.gfa_path);
    if (model.paths.size() != graph.paths.size()) {
        throw std::runtime_error("Internal: path count mismatch between Graph and GfaModel");
    }

    // No output may name an input: the inputs are read throughout the run, and the commit would overwrite
    // them.
    {
        const std::string outs[] = {
            options.out_prefix + ".panphorte.report.tsv",
            options.out_prefix + ".panphorte.copies.tsv",
            options.out_prefix + ".normalized.gfa",
            options.out_prefix + ".normalized.sorted.gfa",
            options.out_prefix + ".bubbles.csv",
            options.out_prefix + ".bandage_nodes.csv",
            options.out_prefix + ".bandage_genes.csv",
            options.out_prefix + ".panphorte.refine.tsv",
            options.out_prefix + ".panphorte.unit_graph.tsv",
        };
        // The GTF is an input like any other, so it is checked against every output.
        for (const std::string& in : {options.gfa_path, options.bubbles_csv_in, options.gtf_path}) {
            if (in.empty()) continue;
            std::error_code ec;
            const auto ip = std::filesystem::weakly_canonical(in, ec);
            if (ec || ip.empty()) continue;
            for (const std::string& o : outs) {
                std::error_code e2;
                const auto op = std::filesystem::weakly_canonical(o, e2);
                if (!e2 && ip == op)
                    throw std::runtime_error("panphorte: output '" + o + "' is the same file as input '" +
                                             in + "'");
            }
        }
    }

    const std::vector<Bubble> bubbles = read_bubbles_csv(options.bubbles_csv_in);
    // The CSV and the graph are separate inputs with nothing tying them together: a CSV built from a
    // different graph, or from this one before a renumbering, would otherwise run against nodes that do not
    // exist.
    {
        std::vector<std::string> missing;
        std::unordered_set<std::size_t> seen_ids;
        for (const Bubble& b : bubbles) {
            if (!seen_ids.insert(b.id).second)
                throw std::runtime_error("panphorte: duplicate bubble id in the CSV: " +
                                         std::to_string(b.id));
            const auto check = [&](const std::string& n) {
                if (graph.nodes.find(n) == graph.nodes.end() && missing.size() < 8)
                    missing.push_back("bubble " + std::to_string(b.id) + " node " + n);
            };
            check(b.source);
            check(b.sink);
            for (const std::string& n : b.inside) check(n);
        }
        if (!missing.empty())
            throw std::runtime_error(
                "panphorte: the bubbles CSV describes nodes the GFA does not contain (" +
                cli::join_with_comma(missing) + "); the CSV and the graph are not the same graph");
    }
    std::unordered_set<std::size_t> bubble_filter(options.bubble_ids.begin(), options.bubble_ids.end());
    // A --bubble-id that names nothing is an error, not an empty selection.
    if (!bubble_filter.empty()) {
        std::unordered_set<std::size_t> present;
        for (const Bubble& b : bubbles) present.insert(b.id);
        std::vector<std::string> absent;
        for (const std::size_t want : bubble_filter)
            if (present.find(want) == present.end()) absent.push_back(std::to_string(want));
        if (!absent.empty()) {
            // Inversion sites are not read (the sites inside them are folded), so name them as such.
            std::unordered_set<std::size_t> inversion_ids;
            for (const Bubble& b : read_bubbles_csv(options.bubbles_csv_in, true))
                if (b.inversion_parent) inversion_ids.insert(b.id);
            std::vector<std::string> inversions, missing_ids;
            for (const std::string& id : absent)
                (inversion_ids.count(std::stoul(id)) ? inversions : missing_ids).push_back(id);
            if (!inversions.empty())
                throw std::runtime_error("panphorte: --bubble-id names an inversion site, which is not folded (the "
                                         "sites inside it are): " + cli::join_with_comma(inversions));
            throw std::runtime_error("panphorte: --bubble-id not present in the bubbles CSV: " +
                                     cli::join_with_comma(missing_ids));
        }
    }
    // Sites must be disjoint: two sharing interior nodes describe the same sequence twice, and folding both
    // rewrites one span inside another. Caught here, naming the pair, rather than as overlapping edits after
    // every haplotype has been aligned. Only among the selected bubbles, so selecting one member of an
    // overlapping pair still works.
    {
        std::unordered_map<std::string, std::size_t> owner;
        for (const Bubble& b : bubbles) {
            if (!bubble_filter.empty() && bubble_filter.find(b.id) == bubble_filter.end()) continue;
            for (const std::string& n : b.inside) {
                const auto it = owner.find(n);
                if (it != owner.end()) {
                    throw std::runtime_error(
                        "panphorte: bubbles " + std::to_string(it->second) + " and " +
                        std::to_string(b.id) + " both claim interior node " + n +
                        "; overlapping sites describe the same sequence twice and folding both would "
                        "rewrite one span inside the other");
                }
                owner.emplace(n, b.id);
            }
        }
    }

    const std::unordered_map<std::string, NodeTok> node_tok = build_node_tokens(graph);

    std::vector<BubblePathIndex> path_indexes;
    path_indexes.reserve(graph.paths.size());
    for (const auto& p : graph.paths) {
        path_indexes.push_back(build_bubble_path_index(p));
    }
    // Node ids alone do not tie the table to this graph (numeric ids nearly always exist, and in order), so
    // each site must be crossed by some haplotype walking only its own interior between the boundaries and
    // spelling an interior length within the table's range. A table from another graph, or from this one
    // before a renumbering, has few such sites: refused when most lack one.
    if (declined.empty()) {
        std::vector<std::size_t> uncrossed;
        std::size_t checked = 0;
        for (const Bubble& b : bubbles) {
            if (!bubble_filter.empty() && bubble_filter.find(b.id) == bubble_filter.end()) continue;
            ++checked;
            const std::unordered_set<std::string> inside(b.inside.begin(), b.inside.end());
            bool clean = false;
            for (std::size_t pi = 0; pi < graph.paths.size() && !clean; ++pi) {
                const std::vector<PathStep>& steps = graph.paths[pi].steps;
                const auto len = [&](std::size_t k) { return graph.nodes.at(steps[k].node_id).sequence.size(); };
                // Passes back to back are measured as one span, boundaries walked again included, as bubble
                // measures them (whether the crossings come apart or as one window); a run is clean while every
                // pass in it walks only the interior.
                std::size_t run_bp = 0, prev_right = 0;
                bool run_clean = false;
                for (const BubblePathInterval& c : bubble_crossings(graph.paths[pi], path_indexes[pi], b)) {
                    bool all_inside = true;
                    std::size_t bp = 0;
                    for (std::size_t k = c.left + 1; k < c.right && all_inside; ++k) {
                        const std::string& id = steps[k].node_id;
                        all_inside = inside.count(id) != 0 || id == b.source || id == b.sink;
                        if (all_inside) bp += len(k);
                    }
                    const bool joins = run_clean && c.left == prev_right + 1;
                    run_bp = joins ? run_bp + len(prev_right) + len(c.left) + bp : bp;
                    run_clean = all_inside;
                    prev_right = c.right;
                    // A table without interior lengths (both 0) is judged on the nodes alone.
                    const auto in_range = [&](std::size_t x) {
                        return (b.min_inside_bp == 0 && b.max_inside_bp == 0) ||
                               (x >= b.min_inside_bp && x <= b.max_inside_bp);
                    };
                    if (all_inside && (in_range(bp) || in_range(run_bp))) { clean = true; break; }
                }
            }
            if (!clean) uncrossed.push_back(b.id);
        }
        if (!uncrossed.empty()) {
            std::string ids;
            for (std::size_t i = 0; i < uncrossed.size() && i < 8; ++i) ids += (i ? ", " : "") + std::to_string(uncrossed[i]);
            if (uncrossed.size() > 8) ids += ", ...";
            if (2 * uncrossed.size() > checked)
                throw std::runtime_error("panphorte: " + std::to_string(uncrossed.size()) + " of " +
                                         std::to_string(checked) + " sites are crossed by no haplotype through "
                                         "their own interior (" + ids + "); the bubbles CSV does not describe "
                                         "this graph (run bubble on it and pass that table)");
            std::cerr << cli::log_prefix("panphorte") << "warning: sites crossed by no haplotype through their own "
                      << "interior: " << uncrossed.size() << " (" << ids << "); check that the table was made from this graph\n";
        }
    }
    // The reference path, whose sequence --lossy-units keeps (SIZE_MAX: none given).
    std::size_t ref_pi = SIZE_MAX;
    if (!options.reference_path.empty()) {
        const std::string name = resolve_reference_path_name(graph, options.reference_path, "panphorte");
        for (std::size_t i = 0; i < graph.paths.size(); ++i)
            if (graph.paths[i].name == name) { ref_pi = i; break; }
    }

    // Approximate mode (--min-similarity < 1): near-identical copies fold onto a unit graph that keeps each
    // copy's differences (only the reference's with --lossy-units). Exact mode folds identical copies onto
    // one looped node.
    const bool approximate = options.min_similarity < 1.0;
    const auto say = [&](const std::string& message) {
        if (!options.quiet) std::cerr << cli::log_prefix("panphorte") << message << '\n';
    };
    // A pass repeated after declining sites does not log its input again.
    if (declined.empty()) {
        say("input: haplotypes: " + std::to_string(model.paths.size()) + "; nodes: " +
            std::to_string(model.node_order.size()) + "; edges: " + std::to_string(model.edges.size()));
        say("reference: " + (ref_pi == SIZE_MAX ? std::string("none") : graph.paths[ref_pi].name));
        say("sites: " + std::to_string(bubbles.size()) +
            (bubble_filter.empty() ? std::string() : "; selected: " + std::to_string(bubble_filter.size())));
        std::ostringstream sim;
        sim << options.min_similarity;
        say(approximate ? "copies: near-identical (similarity at least " + sim.str() + ")" : "copies: identical");
    }
    // Nothing lands until the run succeeds: every output is staged and committed at the end. The output
    // directory is made here, after the input has been read and checked, so a refused run leaves none.
    cli::ensure_parent_dir_for_file(options.out_prefix);
    cli::StagedOutputs staged("panphorte");
    std::set<std::string> written;   // final paths this run stages, for the stale-sibling sweep at the end
    auto stage_out = [&](const std::string& final_path) {
        written.insert(final_path);
        return staged.stage(final_path);
    };

    // Every fold is a unit graph: the nodes a copy of the unit walks, and its anchor (the nodes covering a
    // unit position no copy deletes, where copies are counted). Exact mode's fold is one node looped once per
    // copy, its own anchor. Buffered, because --reference-path renumbers the graph and the table is mapped
    // through the sort when written. Two phase-rotated units at one site are two unit graphs; their shared
    // input_bubble_id says they are one site's copy number.
    struct UnitGraphRecord {
        std::size_t bubble_id = 0;
        std::string unit;
        std::string canonical_motif;   // the unit's primitive motif, least rotation over both strands
        std::size_t anchor_pos = SIZE_MAX;
        std::vector<std::string> anchor;
        std::vector<std::string> nodes;
    };
    std::vector<UnitGraphRecord> unit_graphs;

    std::ofstream copies_out;
    if (approximate) {
        copies_out.open(stage_out(options.out_prefix + ".panphorte.copies.tsv"));
        if (!copies_out) {
            throw std::runtime_error("Failed to write panphorte copies TSV");
        }
        // These ids describe the input graph: under --reference-path the sort renumbers nodes and the
        // re-snarl reassigns bubble ids, so the same numbers then refer to different sequence.
        copies_out << "path_name\tsample\tinput_bubble_id\tcopies\tunit_bp\torientations\t"
                   << "mean_identity\tregion_bp\tinput_from_node\tinput_to_node\n";
    }

    // Unit nodes are site-local, keyed by (bubble id, canonical unit): keyed by sequence alone, two loci
    // carrying the same repeat would share a node, joining them in the graph and corrupting sorting, snarl
    // decomposition and any copy number read from node multiplicity.
    std::unordered_map<std::string, std::string> rep_by_unit; // "<bubble>\t<canonical unit>" -> node
    // Fragment nodes (the literal bases around and between approximate copies) are interned, not minted
    // per occurrence, so identical haplotypes walk identical routes. The key is site- and context-local:
    // bubble, role (lead / between copies / tail), the oriented step it follows, and the sequence read
    // source -> sink.
    std::unordered_map<std::string, std::string> fragment_by_key;
    std::unordered_set<std::string> removal_candidates;
    // Oriented links a path traversed across a span this run replaced, and the site that replaced it: a
    // replaced node another walk still uses would otherwise keep its old arcs, leaving the replaced branch
    // reachable. The bubble id names the site when a link cannot be resolved.
    std::unordered_map<std::string, std::size_t> obsolete_edge_keys;
    std::vector<std::vector<PathArray>> per_path_arrays(graph.paths.size());
    // Every bubble the loop visited, in order, and whether it was folded: the candidate sites of --refine.
    struct SeenSite {
        std::size_t bubble_id;
        std::string source;
        std::string sink;
        bool folded;
    };
    std::vector<SeenSite> seen_sites;

    // Every consecutive pair the original path walks across [lo, hi), plus the arcs into and out of the
    // span: those are the ones the replacement is meant to supersede.
    const auto record_replaced_span = [&](const PathRecord& p, std::size_t bubble_id, std::size_t lo,
                                          std::size_t hi) {
        if (p.steps.empty()) return;
        const std::size_t first = lo == 0 ? 1 : lo;
        const std::size_t last = std::min(hi + 1, p.steps.size());
        for (std::size_t k = first; k < last; ++k) {
            obsolete_edge_keys.emplace(EdgeSet::key(p.steps[k - 1].node_id,
                                                    orient_char(p.steps[k - 1].reverse),
                                                    p.steps[k].node_id,
                                                    orient_char(p.steps[k].reverse)),
                                       bubble_id);
        }
    };

    PanphorteSummary summary;
    std::ofstream report(stage_out(options.out_prefix + ".panphorte.report.tsv"));
    if (!report) {
        throw std::runtime_error("Failed to write panphorte report: " + options.out_prefix + ".panphorte.report.tsv");
    }
    // Enough to explain a decision without re-running: how many haplotypes were at the site, how many
    // carried the winning motif, that fraction, and why the bubble was or was not folded.
    report << "bubble_id\tnormalized\tunit_bp\tpaths_normalized\tmin_copies\tmax_copies\t"
           << "interruptions_bp\tnodes_collapsed\tn_traversing\tn_motif_carriers\tprevalence\t"
           << "n_motifs\tstatus\n";
    // Rows are held until the graph is pruned: nodes_collapsed counts the nodes the fold removed, known only
    // then, since a replaced node another site or haplotype still uses survives the fold.
    struct ReportRow {
        std::string head;                   // bubble_id .. interruptions_bp, tab-terminated
        std::vector<std::string> collapsed; // nodes in the replaced ranges
        std::string tail;                   // n_traversing .. status, newline-terminated
    };
    std::vector<ReportRow> report_rows;

    std::size_t total_bubbles = 0;
    for (const Bubble& b : bubbles) {
        if (bubble_filter.empty() || bubble_filter.find(b.id) != bubble_filter.end()) ++total_bubbles;
    }
    // Exact mode: a bar over the sites. Approximate mode logs each seeded site and shows a bar over its
    // haplotypes instead (below), since one site's alignment can take most of the run.
    const bool bubble_bar = !options.quiet && !approximate;
    cli::ProgressBar progress(bubble_bar ? cli::log_prefix("panphorte") + "folding sites" : std::string(),
                              bubble_bar ? total_bubbles : 0);

    // The sites naming each node (as a boundary or inside), and the sites holding it inside, so a site is
    // widened only over nodes of its own. The nodes a widened site took are added once it folds: a site left
    // unfolded edits nothing and leaves them to its neighbours. Its new boundaries are not added, since a
    // boundary is not edited: a neighbour may take that node and keep its edits beside the site's.
    std::unordered_map<std::string, std::vector<std::size_t>> sites_of_node;
    std::unordered_map<std::string, std::vector<std::size_t>> inside_of_node;
    for (const Bubble& b : bubbles) {
        sites_of_node[b.source].push_back(b.id);
        sites_of_node[b.sink].push_back(b.id);
        for (const std::string& n : b.inside) {
            sites_of_node[n].push_back(b.id);
            inside_of_node[n].push_back(b.id);
        }
    }
    std::vector<std::string> widened_nodes;   // nodes the current site took by widening
    const auto widen = [&](const Bubble& site,
                           const std::function<std::size_t(const std::vector<PathStep>&, bool)>& copy_steps) {
        std::optional<Bubble> w =
            widen_over_shared_copies(graph, path_indexes, site, sites_of_node, inside_of_node, copy_steps);
        widened_nodes.clear();
        if (!w) return w;
        const std::unordered_set<std::string> had(site.inside.begin(), site.inside.end());
        for (const std::string& n : w->inside)
            if (!had.count(n)) widened_nodes.push_back(n);
        return w;
    };
    const auto claim_widened = [&](std::size_t id, bool folded) {
        if (folded)
            for (const std::string& n : widened_nodes) sites_of_node[n].push_back(id);
        widened_nodes.clear();
    };

    for (const Bubble& input_bubble : bubbles) {
        if (!bubble_filter.empty() && bubble_filter.find(input_bubble.id) == bubble_filter.end()) {
            continue;
        }
        Bubble bubble = input_bubble;   // widened over copies shared just outside it, when there are any
        ++summary.bubbles_seen;
        if (declined.count(bubble.id) != 0) {
            // Folding it left its original route walkable (see SharedRouteSites): left as it is.
            seen_sites.push_back({bubble.id, bubble.source, bubble.sink, false});
            ReportRow row;
            row.head = std::to_string(bubble.id) + "\tno\t0\t0\t0\t0\t0\t";
            row.tail = ".\t.\t.\t.\tdeclined_shared_route\n";
            ++summary.status_counts["declined_shared_route"];
            report_rows.push_back(std::move(row));
            progress.tick();
            continue;
        }

        std::size_t paths_norm = 0;
        std::size_t min_copies_seen = 0;
        std::size_t max_copies_seen = 0;
        std::size_t unit_bp_report = 0;
        std::size_t interruptions_bp = 0;
        std::unordered_set<std::string> collapsed_nodes;

        if (approximate) {
            // Seed one unit, find its near-identical copies per path (either orientation, banded alignment)
            // and fold them onto the site's unit graph; sequence between copies stays as literal fragments.
            // Every candidate seed is phased and tested on the same sample interiors; the one detecting the
            // most copies is kept (ties: the earlier, the group's closest member).
            UnitPhase phase;
            for (const std::string& cand : pick_reference_units(graph, node_tok, path_indexes, bubble, options)) {
                UnitPhase p = phase_unit_at_array_edge(graph, path_indexes, bubble, cand, options);
                if (phase.unit.empty() || p.copies_phased > phase.copies_phased) phase = std::move(p);
            }
            const std::string& ref_unit = phase.unit;
            if (ref_unit.size() >= options.min_unit_bp) {
                // One line per seeded site, printed before the costly per-haplotype alignment so a slow site
                // shows in non-interactive logs too. Sites with no foldable unit stay silent.
                say("site " + std::to_string(bubble.id) + ": unit: " + std::to_string(ref_unit.size()) + " bp" +
                    (phase.shift == 0 ? std::string()
                                      : "; rotated to start at the array edge: " + std::to_string(phase.shift) + " bp"));
            }
            if (ref_unit.size() >= options.min_unit_bp && ref_unit.size() >= 16) {
                // Copies of the unit every crossing shares just outside the site: widen it over them. A copy
                // counts when it ends (or, past the sink, starts) within an alignment band of the boundary.
                const std::size_t band =
                    static_cast<std::size_t>((1.0 - options.min_similarity) * static_cast<double>(ref_unit.size())) + 8;
                const auto copy_steps = [&](const std::vector<PathStep>& outward, bool left) -> std::size_t {
                    std::string seq;
                    std::vector<std::size_t> reach;   // bases from the boundary out to the end of each step
                    for (const PathStep& st : outward) {
                        const std::string b = step_sequence(graph, st);
                        seq = left ? b + seq : seq + b;
                        reach.push_back(seq.size());
                        if (seq.size() >= ref_unit.size() + 2 * band) break;
                    }
                    std::optional<ApproxCopy> nearest;
                    for (const ApproxCopy& c : detect_copies_native(seq, ref_unit, options.min_similarity)) {
                        if (c.rev) continue;
                        if (!nearest || (left ? c.bp_hi > nearest->bp_hi : c.bp_lo < nearest->bp_lo)) nearest = c;
                    }
                    if (!nearest) return 0;
                    const std::size_t gap = left ? seq.size() - nearest->bp_hi : nearest->bp_lo;
                    if (gap > band) return 0;
                    const std::size_t need = left ? seq.size() - nearest->bp_lo : nearest->bp_hi;
                    for (std::size_t n = 0; n < reach.size(); ++n)
                        if (reach[n] >= need) return n + 1;
                    return 0;
                };
                if (std::optional<Bubble> w = widen(bubble, copy_steps)) {
                    bubble = std::move(*w);
                    say("site " + std::to_string(bubble.id) + ": widened to " + bubble.source + " .. " + bubble.sink +
                        " over copies every haplotype shares beside it");
                }
            }
            // Hoisted so the report can state them whether or not the inner block runs.
            std::size_t approx_traversing = 0, approx_carriers = 0;
            const bool seeded = ref_unit.size() >= options.min_unit_bp;
            double approx_prevalence = 0.0;
            // A unit made of copies of a motif shorter than --min-unit-bp is a short tandem repeat the size
            // floor made look like a copy unit: folding it would count copies of the wrong thing, so the site
            // stays sequence for call. It is one when, for some period p, the most common base of each phase
            // covers 90% of it (bubble's tandem test), or, where an indel shifts the phase part-way, it fits a
            // run of one of its own p-mers at 90% identity with indels over at least three copies (two copies
            // of a skewed alphabet fit by chance). The cheap screen for the fit, each base against the one p
            // further on, also passes homopolymer runs, which no one p-mer fits.
            bool unit_not_primitive = false;
            if (approximate) {
                const std::size_t n = ref_unit.size();
                std::string unit_uc = ref_unit;
                for (char& c : unit_uc) c = static_cast<char>(c & 0xDF);
                for (std::size_t p = 1; 2 * p <= n && p < options.min_unit_bp && !unit_not_primitive; ++p) {
                    std::size_t same = 0;
                    for (std::size_t i = 0; i + p < n; ++i) same += unit_uc[i] == unit_uc[i + p];
                    if (3 * p <= n && 10 * same >= 6 * (n - p)) {   // two copies: the consensus test
                        const std::size_t band = n / 10 + 8;
                        for (std::size_t q = 0; q < 4 && !unit_not_primitive; ++q) {
                            const std::size_t o = (n - p) * q / 4;    // a p-mer from each quarter
                            const std::size_t ph = o % p;             // in the phase the unit starts in
                            std::string rot = unit_uc.substr(o, p);
                            rot = rot.substr((p - ph) % p) + rot.substr(0, (p - ph) % p);
                            std::string reps;
                            while (reps.size() < n + band + p) reps += rot;
                            const FitAlignResult fa = fit_align(unit_uc, reps, band);
                            unit_not_primitive = fa.ok && fa.identity >= 0.9;
                        }
                    }
                    if (unit_not_primitive) break;
                    std::size_t m = 0;
                    for (std::size_t k = 0; k < p; ++k) {
                        std::size_t c[4] = {0, 0, 0, 0};
                        for (std::size_t i = k; i < n; i += p) {
                            switch (ref_unit[i] & 0xDF) {
                            case 'A': ++c[0]; break;
                            case 'C': ++c[1]; break;
                            case 'G': ++c[2]; break;
                            case 'T': ++c[3]; break;
                            default: break;
                            }
                        }
                        m += std::max({c[0], c[1], c[2], c[3]});
                    }
                    unit_not_primitive = static_cast<double>(m) >= 0.9 * static_cast<double>(n);
                }
            }
            if (ref_unit.size() >= options.min_unit_bp) {
                // One replaced step range holding one or more copies. Several copies can share a range
                // -- two copies inside a single node have no step boundary between them -- so the range
                // is emitted as fragment, unit, fragment, unit, ... fragment: `frag` always holds one
                // more entry than there are copies, and the outer two are the bases before the first
                // copy and after the last.
                struct Mapped {
                    std::size_t off_lo = 0, off_hi = 0;
                    std::vector<bool> rev;
                    std::vector<double> id;
                    std::vector<std::string> frag;
                    std::vector<std::string> seq;   // each copy's own bases, as the path spells them
                    std::size_t last_bp_hi = 0;   // end of the most recent copy, in interval bp
                    std::size_t n_copies() const { return rev.size(); }
                };
                // One per crossing, in canonical (source -> sink) orientation: a haplotype stored the other
                // way round is read as it would be forward, so it gets the same copy boundaries, fragments and
                // unit-graph walk as one stored forward. `left`/`n` place the interior on the path; a backward
                // crossing's canonical offset k is the path's step left + n - 1 - k, flipped. `before` is the
                // canonical step before the interior (the source as the crossing enters it).
                struct Res {
                    std::size_t left = 0, n = 0;
                    bool backward = false;
                    PathStep before;
                    std::vector<Mapped> copies;
                };
                std::vector<std::vector<Res>> results(graph.paths.size());
                const bool crossing_fallback = !bubble_orientation_holds(graph.paths, path_indexes, bubble);

                std::atomic<std::size_t> n_traverse{0};   // haplotypes that traverse this bubble
                std::size_t approx_interruptions_bp = 0;    // real interrupting bases within arrays
                std::size_t group_carriers = 0;             // paths with an array reaching min_copies
                std::mutex detail_mtx;
                cli::ProgressBar align_bar(options.quiet ? std::string()
                                                         : cli::log_prefix("panphorte") + "site " +
                                                               std::to_string(bubble.id) + ": aligning haplotypes",
                                           graph.paths.size());

                auto detect_one = [&](std::size_t pi) {
                    // Every crossing of the site on its own (bubble_crossings): a haplotype crossing it twice
                    // has two sets of copies and two edits. A haplotype crossing with no interior carries
                    // zero copies and is still a haplotype at this site, so it counts as traversing.
                    const std::vector<BubblePathInterval> crossings =
                        bubble_crossings(graph.paths[pi], path_indexes[pi], bubble, crossing_fallback);
                    if (crossings.empty()) return;
                    n_traverse.fetch_add(1);
                    std::size_t path_interruptions = 0;
                    std::size_t best_group = 0;
                    std::vector<Res> found;
                    for (const BubblePathInterval& used : crossings) {
                        if (used.inside_count == 0) continue;      // crosses with no copies
                        // Search the interior only: the source and sink are backbone, not bubble content, and
                        // must not contribute copies.
                        const std::size_t left = used.left + 1;
                        if (used.right <= left) continue;
                        const std::size_t n = used.right - left;
                        std::vector<PathStep> nat(
                            graph.paths[pi].steps.begin() + static_cast<std::ptrdiff_t>(left),
                            graph.paths[pi].steps.begin() + static_cast<std::ptrdiff_t>(left + n));
                        const bool backward = !used.source_to_sink;
                        PathStep before = graph.paths[pi].steps[used.left];
                        if (backward) {
                            std::reverse(nat.begin(), nat.end());
                            for (PathStep& st : nat) st.reverse = !st.reverse;
                            before = graph.paths[pi].steps[used.right];
                            before.reverse = !before.reverse;
                        }
                        const std::string snat = spell_path_steps_sequence(graph, nat);
                        std::vector<std::size_t> prefix(n + 1, 0); // cumulative bp at step boundaries
                        for (std::size_t off = 0; off < n; ++off) {
                            const auto it = node_tok.find(nat[off].node_id);
                            prefix[off + 1] = prefix[off] + (it != node_tok.end() ? it->second.len : 0);
                        }
                        const auto bp_copies = detect_copies_native(snat, ref_unit, options.min_similarity);
                        // A copy boundary rarely lands on a node boundary. The containing step range is
                        // taken -- floor for the start, ceiling for the end -- and the bases inside it
                        // that fall outside the copy are kept as literal fragments, emitted around the
                        // unit step at rewrite time: nothing outside a copy is lost, and no copy is declined.
                        // Binary searches over the step prefix.
                        const auto floor_off = [&](std::size_t bp) {
                            // last index lo with prefix[lo] <= bp, capped so lo + 1 stays in range
                            const auto it = std::upper_bound(prefix.begin() + 1, prefix.end(), bp);
                            const std::size_t lo = static_cast<std::size_t>(it - prefix.begin()) - 1;
                            return std::min(lo, prefix.size() >= 2 ? prefix.size() - 2 : 0);
                        };
                        const auto ceil_off = [&](std::size_t bp) {
                            const auto it = std::lower_bound(prefix.begin(), prefix.end(), bp);
                            const std::size_t hi = static_cast<std::size_t>(it - prefix.begin());
                            return hi < prefix.size() ? hi : prefix.size() - 1;
                        };
                        // Copies are grouped into blocks by the step range that contains them: two copies
                        // inside one node share a range, and mapped independently the second one's floor would
                        // collide with the first one's ceiling and the copy be lost.
                        std::vector<Mapped> copies;
                        for (const ApproxCopy& c : bp_copies) {
                            const std::size_t lo = floor_off(c.bp_lo);
                            const std::size_t hi = ceil_off(c.bp_hi);
                            if (hi <= lo) continue;   // cannot happen for a copy of positive length
                            if (!copies.empty() && lo < copies.back().off_hi) {
                                // Shares the open block: extend it and record the gap since the last copy.
                                Mapped& blk = copies.back();
                                if (hi > blk.off_hi) blk.off_hi = hi;
                                const std::size_t gap_lo = blk.last_bp_hi;
                                const std::size_t gap_hi = std::max(c.bp_lo, gap_lo);
                                blk.frag.push_back(snat.substr(gap_lo, gap_hi - gap_lo));
                                blk.rev.push_back(c.rev);
                                blk.id.push_back(c.identity);
                                blk.last_bp_hi = std::max(c.bp_hi, gap_hi);
                                blk.seq.push_back(snat.substr(gap_hi, blk.last_bp_hi - gap_hi));
                                continue;
                            }
                            Mapped blk;
                            blk.off_lo = lo;
                            blk.off_hi = hi;
                            const std::size_t range_lo = prefix[lo];
                            const std::size_t copy_lo = std::max(c.bp_lo, range_lo);
                            blk.frag.push_back(snat.substr(range_lo, copy_lo - range_lo));
                            blk.rev.push_back(c.rev);
                            blk.id.push_back(c.identity);
                            blk.last_bp_hi = std::max(copy_lo, std::min(c.bp_hi, prefix[hi]));
                            blk.seq.push_back(snat.substr(copy_lo, blk.last_bp_hi - copy_lo));
                            copies.push_back(std::move(blk));
                        }
                        // Close each block with the bases between its last copy and the end of its range.
                        for (Mapped& blk : copies) {
                            const std::size_t range_hi = prefix[blk.off_hi];
                            const std::size_t tail_lo = std::min(blk.last_bp_hi, range_hi);
                            blk.frag.push_back(snat.substr(tail_lo, range_hi - tail_lo));
                        }
                        // Copies are grouped into tandem arrays by the bases between consecutive copies,
                        // whatever the node boundaries, and --max-interruption-frac decides where a group ends.
                        {
                            bool open = false;
                            std::size_t gstart_lo = 0, ginterrupt = 0, gcount = 0, prev_hi = 0;
                            for (const ApproxCopy& c : bp_copies) {
                                if (!open) {
                                    open = true;
                                    gstart_lo = c.bp_lo;
                                    ginterrupt = 0;
                                    gcount = 1;
                                    prev_hi = c.bp_hi;
                                    continue;
                                }
                                const std::size_t gap = c.bp_lo > prev_hi ? c.bp_lo - prev_hi : 0;
                                const std::size_t span = c.bp_hi - gstart_lo;
                                const double frac = span > 0
                                    ? static_cast<double>(ginterrupt + gap) / static_cast<double>(span) : 0.0;
                                if (frac > options.max_interruption_frac) {
                                    best_group = std::max(best_group, gcount);
                                    path_interruptions += ginterrupt;
                                    gstart_lo = c.bp_lo;
                                    ginterrupt = 0;
                                    gcount = 1;
                                } else {
                                    ginterrupt += gap;
                                    ++gcount;
                                }
                                prev_hi = std::max(prev_hi, c.bp_hi);
                            }
                            if (open) {
                                best_group = std::max(best_group, gcount);
                                path_interruptions += ginterrupt;
                            }
                        }
                        // Record any crossing with >=1 detected copy. The min_copies gate is per bubble, below:
                        // once an array is confirmed, single-copy haplotypes fold through the unit too, so
                        // their copy number is 1, not 0.
                        if (!copies.empty()) found.push_back(Res{left, n, backward, before, std::move(copies)});
                    }
                    {
                        std::lock_guard<std::mutex> lk(detail_mtx);
                        approx_interruptions_bp += path_interruptions;
                        if (best_group >= options.min_copies) ++group_carriers;
                    }
                    results[pi] = std::move(found);
                    {
                        std::lock_guard<std::mutex> lk(detail_mtx);
                        align_bar.tick();
                    }
                };

                run_parallel(graph.paths.size(), options.threads, detect_one);
                align_bar.done();

                // Confirm a real tandem array before folding: some path must carry >= min_copies copies, and
                // >= --min-array-prevalence of traversing haplotypes an array of them, so a rare private
                // duplication is not folded over genes `call` resolves. Carriers are counted by group, not
                // total copies: two copies further apart than --max-interruption-frac are two arrays of one.
                const auto path_copies = [&](std::size_t pi) {
                    std::size_t nc = 0;
                    for (const Res& r : results[pi])
                        for (const auto& blk : r.copies) nc += blk.n_copies();
                    return nc;
                };
                std::size_t max_copies_path = 0;
                for (std::size_t pi = 0; pi < results.size(); ++pi)
                    max_copies_path = std::max(max_copies_path, path_copies(pi));
                const std::size_t traversing = n_traverse.load();
                const double prevalence =
                    traversing > 0 ? static_cast<double>(group_carriers) / static_cast<double>(traversing) : 0.0;
                const bool array_confirmed =
                    max_copies_path >= options.min_copies && prevalence >= options.min_array_prevalence &&
                    !unit_not_primitive;
                approx_traversing = traversing;
                approx_carriers = group_carriers;
                // Interrupting bases within arrays, for the report.
                interruptions_bp = approx_interruptions_bp;
                approx_prevalence = prevalence;

                // One unit graph per bubble (see build_unit_graph), from the copies whose sequence is kept: every
                // copy by default; with --lossy-units the reference's only, every other copy walking the plain
                // unit (its copy number kept, its own differences dropped). The unit itself is always included,
                // for that plain walk; where no kept copy differs from it the graph is one node, the unit
                // looped once per copy, as a single-node fold.
                std::optional<UnitGraph> ugraph;
                if (array_confirmed) {
                    std::set<std::string> distinct{ref_unit};
                    for (std::size_t pi = 0; pi < results.size(); ++pi) {
                        if (!options.lossless_units && pi != ref_pi) continue;
                        for (const Res& r : results[pi])
                            for (const Mapped& m : r.copies)
                                for (std::size_t ci = 0; ci < m.n_copies(); ++ci)
                                    distinct.insert(m.rev[ci] ? reverse_complement(m.seq[ci]) : m.seq[ci]);
                    }
                    ugraph = build_unit_graph(ref_unit, std::vector<std::string>(distinct.begin(), distinct.end()),
                                              options.threads);
                }
                // Serial: one PathArray and one copies.tsv row per crossing.
                for (std::size_t pi = 0; array_confirmed && pi < results.size(); ++pi) {
                    if (results[pi].empty()) continue;
                    const bool keep_sequence = options.lossless_units || pi == ref_pi;
                    std::size_t path_occ = 0;   // the haplotype's copies, over all its crossings
                    for (const Res& res : results[pi]) {
                    const std::size_t left = res.left;
                    const std::vector<Mapped>& copies = res.copies;
                    // Canonical offset -> the step as the crossing reads it, and its path index.
                    const auto path_index = [&](std::size_t off) { return res.backward ? left + res.n - 1 - off : left + off; };
                    const auto step_at = [&](std::size_t off) {
                        PathStep st = graph.paths[pi].steps[path_index(off)];
                        if (res.backward) st.reverse = !st.reverse;
                        return st;
                    };
                    // Canonical [lo, hi) -> the path's [lo, hi).
                    const auto path_span = [&](std::size_t lo, std::size_t hi) {
                        return res.backward ? std::make_pair(left + res.n - hi, left + res.n - lo)
                                            : std::make_pair(left + lo, left + hi);
                    };

                    PathArray pa;
                    pa.bubble_id = bubble.id;
                    const std::size_t first_lo = copies.front().off_lo;
                    std::size_t prev_hi = copies.front().off_lo;
                    std::string orients;
                    double id_sum = 0.0;
                    for (const Mapped& m : copies) {
                        for (std::size_t off = prev_hi; off < m.off_lo; ++off) {
                            pa.replacement.push_back(step_at(off));
                        }
                        // fragment, unit, fragment, unit, ... fragment. The fragments carry every base of the
                        // replaced range outside a copy, verbatim, as interned nodes: nothing outside an
                        // accepted copy is lost, and copies that share a node all fold.
                        const auto fragment_node = [&](char role, const std::string& bases) {
                            const PathStep& prev = pa.replacement.empty() ? res.before : pa.replacement.back();
                            const std::string key = std::to_string(bubble.id) + '\t' + role + '\t' + prev.node_id +
                                                    (prev.reverse ? "-" : "+") + '\t' + bases;
                            const auto fit = fragment_by_key.find(key);
                            if (fit != fragment_by_key.end()) return fit->second;
                            const std::string id = add_new_node(model, bases);
                            fragment_by_key.emplace(key, id);
                            ++summary.nodes_added;
                            ++summary.fragment_nodes_added;
                            return id;
                        };
                        for (std::size_t ci = 0; ci < m.n_copies(); ++ci) {
                            if (!m.frag[ci].empty()) {
                                pa.replacement.push_back(PathStep{fragment_node(ci == 0 ? 'L' : 'G', m.frag[ci]), false});
                            }
                            // The copy's own walk through the unit graph, or the plain unit's for a copy whose
                            // sequence is not kept; a reverse copy walks it backwards.
                            const std::string& as_unit = !keep_sequence ? ref_unit
                                                         : m.rev[ci] ? reverse_complement(m.seq[ci]) : m.seq[ci];
                            const std::vector<PathStep> w = unit_graph_walk(*ugraph, model, as_unit);
                            if (!m.rev[ci]) {
                                pa.replacement.insert(pa.replacement.end(), w.begin(), w.end());
                            } else {
                                for (auto it = w.rbegin(); it != w.rend(); ++it)
                                    pa.replacement.push_back(PathStep{it->node_id, !it->reverse});
                            }
                        }
                        if (!m.frag.back().empty()) {
                            pa.replacement.push_back(PathStep{fragment_node('T', m.frag.back()), false});
                        }
                        // The replaced steps, and the arcs across them, become candidates for removal;
                        // whether they go is decided globally, since another site may still use them.
                        for (std::size_t off = m.off_lo; off < m.off_hi; ++off) {
                            const std::string& nid = graph.paths[pi].steps[path_index(off)].node_id;
                            removal_candidates.insert(nid);
                            collapsed_nodes.insert(nid);
                        }
                        const auto [span_lo, span_hi] = path_span(m.off_lo, m.off_hi);
                        record_replaced_span(graph.paths[pi], bubble.id, span_lo, span_hi);
                        for (std::size_t ci = 0; ci < m.n_copies(); ++ci) {
                            orients += (m.rev[ci] ? '-' : '+');
                            id_sum += m.id[ci];
                        }
                        prev_hi = m.off_hi;
                    }
                    // Into the path: a backward crossing takes the canonical replacement reversed, each step flipped.
                    pa.start = path_span(first_lo, prev_hi).first;
                    pa.length = prev_hi - first_lo;
                    if (res.backward) {
                        std::reverse(pa.replacement.begin(), pa.replacement.end());
                        for (PathStep& st : pa.replacement) st.reverse = !st.reverse;
                    }
                    per_path_arrays[pi].push_back(std::move(pa));

                    std::size_t occ = 0;
                    for (const Mapped& m : copies) occ += m.n_copies();
                    path_occ += occ;
                    if (unit_bp_report == 0) unit_bp_report = ref_unit.size();
                    std::size_t region_bp = 0;
                    for (std::size_t off = copies.front().off_lo; off < copies.back().off_hi; ++off) {
                        const auto it = node_tok.find(graph.paths[pi].steps[path_index(off)].node_id);
                        if (it != node_tok.end()) region_bp += it->second.len;
                    }
                    copies_out << graph.paths[pi].name << '\t' << sample_of(graph.paths[pi].name) << '\t'
                               << bubble.id << '\t' << occ << '\t' << ref_unit.size() << '\t'
                               << orients << '\t' << (id_sum / static_cast<double>(occ)) << '\t'
                               << region_bp << '\t'
                               << step_at(copies.front().off_lo).node_id << '\t'
                               << step_at(copies.back().off_hi - 1).node_id << '\n';
                    }
                    ++paths_norm;
                    if (min_copies_seen == 0 || path_occ < min_copies_seen) min_copies_seen = path_occ;
                    if (path_occ > max_copies_seen) max_copies_seen = path_occ;
                }
                if (ugraph && !ugraph->nodes.empty()) {
                    summary.nodes_added += ugraph->nodes.size();
                    UnitGraphRecord ug;
                    ug.bubble_id = bubble.id;
                    ug.unit = ugraph->unit;
                    ug.canonical_motif = canonical_motif_key(ugraph->unit);
                    ug.anchor_pos = ugraph->anchor_pos;
                    ug.anchor = unit_graph_anchor_nodes(*ugraph);
                    ug.nodes = ugraph->nodes;
                    unit_graphs.push_back(std::move(ug));
                }
            }

            const bool norm_ok = paths_norm > 0;
            if (norm_ok) ++summary.bubbles_normalized;
            // A site left unfolded reports the unit it tried and no interruptions, as exact mode does.
            if (!norm_ok) {
                if (seeded) unit_bp_report = ref_unit.size();
                interruptions_bp = 0;
            }
            // A site widened but left unfolded keeps its own boundaries.
            const Bubble& seen = norm_ok ? bubble : input_bubble;
            seen_sites.push_back({bubble.id, seen.source, seen.sink, norm_ok});
            claim_widened(bubble.id, norm_ok);
            {
                const char* status = norm_ok ? "normalized"
                                     : ref_unit.empty() || ref_unit.size() < options.min_unit_bp ? "no_seed"
                                     // The copy detector seeds on 16-mers of the unit, so a shorter one finds
                                     // no copy: said as such, not as a prevalence result.
                                     : ref_unit.size() < 16 ? "unit_under_16bp"
                                     : unit_not_primitive ? "short_motif" : "below_prevalence";
                ++summary.status_counts[status];
                ReportRow row;
                std::ostringstream h, t;
                h << bubble.id << '\t' << (norm_ok ? "yes" : "no") << '\t'
                  << unit_bp_report << '\t' << paths_norm << '\t'
                  << min_copies_seen << '\t' << max_copies_seen << '\t'
                  << interruptions_bp << '\t';
                // Without a seed nothing was counted, so these are unknown ('.') rather than zero.
                t << (seeded ? std::to_string(approx_traversing) : ".") << '\t'
                  << (seeded ? std::to_string(approx_carriers) : ".") << '\t'
                  << (seeded ? std::to_string(approx_prevalence) : ".") << '\t'
                  << (ref_unit.empty() ? 0 : 1) << '\t' << status << '\n';
                row.head = h.str();
                row.tail = t.str();
                row.collapsed.assign(collapsed_nodes.begin(), collapsed_nodes.end());
                report_rows.push_back(std::move(row));
            }
            progress.tick();
            continue; // bubble fully handled in approximate mode
        }

        // Detect per-path tandem arrays, then apply the same prevalence gate as approximate mode: fold only
        // if a >= min_copies array is carried by >= --min-array-prevalence of the traversing haplotypes, so a
        // rare private duplication of a gene or segmental module is not folded.
        std::vector<std::vector<TandemArray>> path_arrays(graph.paths.size());
        std::size_t n_traverse = 0;
        bool crossing_fallback = false;
        // Carriers of each candidate motif, so prevalence is asked per motif rather than of any tandem at all.
        std::unordered_map<std::string, std::size_t> carriers_by_motif;
        std::unordered_map<std::string, std::size_t> unit_bp_by_motif;   // the longest array unit of each motif
        std::vector<std::vector<BubblePathInterval>> crossings_of(graph.paths.size());
        const auto detect_site = [&]() {   // fills the above for the site as it stands
            path_arrays.assign(graph.paths.size(), {});
            n_traverse = 0;
            crossing_fallback = !bubble_orientation_holds(graph.paths, path_indexes, bubble);
            carriers_by_motif.clear();
            unit_bp_by_motif.clear();
            crossings_of.assign(graph.paths.size(), {});
            for (std::size_t pi = 0; pi < graph.paths.size(); ++pi) {
                // Every crossing (bubble_crossings): a haplotype crossing with no interior node carries zero
                // copies and still counts in the denominator, or a lone carrier among deletion alleles would read
                // as full prevalence.
                crossings_of[pi] = bubble_crossings(graph.paths[pi], path_indexes[pi], bubble, crossing_fallback);
                if (crossings_of[pi].empty()) continue;
                ++n_traverse;
                std::vector<TandemArray> arrays;
                // A crossing is read source -> sink: one walked the other way is searched on the path read in
                // reverse and its arrays mirrored back, so a haplotype stored reversed gets the same copy
                // boundaries (and leftover partial copy, at the same end) as one stored forward.
                const std::size_t N = graph.paths[pi].steps.size();
                std::optional<PathRecord> rev_path;
                for (const BubblePathInterval& used : crossings_of[pi]) {
                    if (used.inside_count == 0) continue;             // crosses with 0 copies
                    // The interior only, as every other detection path: the boundaries are backbone and must
                    // not be folded into the unit node.
                    if (used.right < used.left + 2) continue;
                    if (used.source_to_sink) {
                        auto found = detect_tandems(
                            graph, node_tok, graph.paths[pi], used.left + 1, used.right - 1,
                            options.min_unit_bp, options.min_copies, options.max_interruption_frac);
                        for (TandemArray& a : found) arrays.push_back(std::move(a));
                        continue;
                    }
                    if (!rev_path) rev_path = reversed_path(graph.paths[pi]);
                    auto found = detect_tandems(
                        graph, node_tok, *rev_path, N - used.right, N - 2 - used.left,
                        options.min_unit_bp, options.min_copies, options.max_interruption_frac);
                    for (const TandemArray& a : found) arrays.push_back(mirror_array(a, N));
                }
                if (arrays.empty()) continue;
                std::unordered_set<std::string> motifs_here;
                for (const TandemArray& a : arrays) {
                    const std::string m = canonical_motif_key(a.unit_seq);
                    motifs_here.insert(m);
                    std::size_t& bp = unit_bp_by_motif[m];
                    bp = std::max(bp, a.unit_seq.size());
                }
                for (const std::string& m : motifs_here) ++carriers_by_motif[m];
                path_arrays[pi] = std::move(arrays);
            }
        };
        detect_site();
        {   // Copies of an array's unit every crossing shares just outside the site: widen it over them.
            std::set<std::pair<std::string, std::size_t>> units;   // (motif, unit length) of the arrays
            std::size_t longest = 0;
            for (const auto& arrays : path_arrays)
                for (const TandemArray& a : arrays) {
                    const std::string m = canonical_motif_key(a.unit_seq);
                    if (m.size() < options.min_unit_bp) continue;
                    units.emplace(m, a.unit_seq.size());
                    longest = std::max(longest, a.unit_seq.size());
                }
            const auto copy_steps = [&](const std::vector<PathStep>& outward, bool left) -> std::size_t {
                std::string seq;
                for (std::size_t n = 1; n <= outward.size(); ++n) {
                    const std::string st = step_sequence(graph, outward[n - 1]);
                    seq = left ? st + seq : seq + st;
                    if (seq.size() > longest) return 0;
                    if (seq.size() < options.min_unit_bp) continue;
                    if (units.count({canonical_motif_key(seq), seq.size()})) return n;
                }
                return 0;
            };
            if (!units.empty())
                if (std::optional<Bubble> w = widen(bubble, copy_steps)) {
                    bubble = std::move(*w);
                    detect_site();
                }
        }
        // Prevalence is asked of each motif on its own: two unrelated repeats in one bubble do not lend each
        // other support, and a motif that clears the gate is folded into its own site-local unit node.
        std::set<std::string> confirmed_motifs;   // ordered: the search order fixes unit node ids
        double prevalence = 0.0;
        bool short_motif_only = !carriers_by_motif.empty();
        for (const auto& [m, c] : carriers_by_motif) {
            const double p = n_traverse > 0 ? static_cast<double>(c) / static_cast<double>(n_traverse) : 0.0;
            // A unit that is copies of a motif shorter than --min-unit-bp is a short tandem repeat the
            // length floor made look like a copy unit, and is left as sequence for `call`, as in approximate
            // mode: counted in whole units, a haplotype one motif longer would read a whole copy up. The key is
            // the primitive motif, so its length tells.
            if (m.size() < options.min_unit_bp) continue;
            short_motif_only = false;
            if (p >= options.min_array_prevalence) confirmed_motifs.insert(m);
            if (p > prevalence) prevalence = p;
        }
        const bool array_confirmed = !confirmed_motifs.empty();
        // Once the site is confirmed, every haplotype carrying at least one copy of a confirmed motif folds:
        // --min-copies confirms the site, not each haplotype, and a one-copy haplotype left off the unit node
        // would read as CN 0.
        if (array_confirmed) {
            for (std::size_t pi = 0; pi < graph.paths.size(); ++pi) {
                // Every haplotype is searched for single copies of each confirmed motif outside its arrays of
                // that motif, so a lone copy beyond an array's end folds too rather than staying literal (which
                // would also keep the old route walkable where another haplotype's copy there folds).
                std::vector<std::pair<std::size_t, std::size_t>> taken;   // [start, end) of folded arrays
                for (const TandemArray& a : path_arrays[pi])
                    if (confirmed_motifs.count(canonical_motif_key(a.unit_seq)))
                        taken.emplace_back(a.start, a.start + a.length);
                const std::size_t N = graph.paths[pi].steps.size();
                std::optional<PathRecord> rev_path;
                for (const BubblePathInterval& used : crossings_of[pi]) {
                    if (used.inside_count == 0) continue;
                    if (!used.source_to_sink && !rev_path) rev_path = reversed_path(graph.paths[pi]);
                    for (const std::string& m : confirmed_motifs) {
                        // Read source -> sink, as the arrays above, passing over the copies already taken.
                        std::vector<TandemArray> singles;
                        if (used.source_to_sink) {
                            singles = find_motif_copies_by_node(
                                graph, graph.paths[pi], used.left, used.right, m, options.min_unit_bp, taken);
                        } else {
                            std::vector<std::pair<std::size_t, std::size_t>> taken_rev;
                            for (const auto& [lo, hi] : taken) taken_rev.emplace_back(N - hi, N - lo);
                            for (const TandemArray& a : find_motif_copies_by_node(
                                     graph, *rev_path, N - 1 - used.right, N - 1 - used.left, m, options.min_unit_bp,
                                     taken_rev))
                                singles.push_back(mirror_array(a, N));
                        }
                        for (TandemArray& a : singles) {
                            taken.emplace_back(a.start, a.start + a.length);
                            path_arrays[pi].push_back(std::move(a));
                        }
                    }
                }
            }
        }
        // The arrays that fold, each split into its smallest whole units.
        std::vector<std::vector<TandemArray>> fold_arrays(graph.paths.size());
        for (std::size_t pi = 0; array_confirmed && pi < graph.paths.size(); ++pi)
            for (const TandemArray& found_arr : path_arrays[pi])
                // Only confirmed motifs fold: an unrelated array in the same bubble does not ride on their
                // prevalence.
                if (confirmed_motifs.count(canonical_motif_key(found_arr.unit_seq)))
                    fold_arrays[pi].push_back(split_repeated_unit(found_arr, options.min_unit_bp));
        // A unit node is keyed by the lesser strand of its unit, so copies walked either way share it, and is
        // stored the way its copies read along the reference's direction through the site: the reference's
        // own copies decide, then the majority of all copies. The sort turns around only nodes the reference
        // walks, so a unit stored on the other strand would have every copy walk it reversed and read as
        // inverted wherever the reference carries no copy.
        bool ref_backward = false;   // the reference crosses the site sink -> source
        if (ref_pi != SIZE_MAX && !crossings_of[ref_pi].empty()) {
            ref_backward = true;
            for (const BubblePathInterval& used : crossings_of[ref_pi]) ref_backward &= !used.source_to_sink;
        }
        const auto copy_unit = [](const TandemArray& arr, const ArrayElement& el) -> const std::string& {
            return el.unit.empty() ? arr.unit_seq : el.unit;
        };
        std::unordered_map<std::string, std::pair<long, long>> strand_votes;   // canonical -> (reference, all)
        for (std::size_t pi = 0; pi < graph.paths.size(); ++pi)
            for (const TandemArray& arr : fold_arrays[pi])
                for (const ArrayElement& el : arr.elements) {
                    if (!el.is_copy) continue;
                    const std::string& unit = copy_unit(arr, el);
                    const std::string rc = reverse_complement(unit);
                    const std::string& along = arr.backward == ref_backward ? unit : rc;
                    const long v = along == std::min(unit, rc) ? 1 : -1;
                    auto& [ref_vote, all_vote] = strand_votes[std::min(unit, rc)];
                    if (pi == ref_pi) ref_vote += v;
                    all_vote += v;
                }
        for (std::size_t pi = 0; pi < graph.paths.size(); ++pi) {
            if (fold_arrays[pi].empty()) continue;
            bool rewritten = false;   // counted only if one of its arrays is actually folded
            std::size_t path_occ = 0; // the haplotype's copies, summed over its arrays
            for (const TandemArray& arr : fold_arrays[pi]) {
                const std::string& unit = arr.unit_seq;
                // The unit node of each copy, and the strand the copy walks it on.
                std::vector<PathStep> copy_steps;
                for (const ArrayElement& el : arr.elements) {
                    if (!el.is_copy) continue;
                    const std::string& own = copy_unit(arr, el);
                    const std::string rc = reverse_complement(own);
                    const std::string canonical = std::min(own, rc);
                    const auto& [ref_vote, all_vote] = strand_votes.at(canonical);
                    const bool canonical_along = ref_vote != 0 ? ref_vote > 0 : all_vote >= 0;
                    const std::string stored = canonical_along ? canonical : reverse_complement(canonical);
                    const std::string rep_key = std::to_string(bubble.id) + "\t" + canonical;
                    auto rit = rep_by_unit.find(rep_key);
                    std::string rep_id;
                    if (rit == rep_by_unit.end()) {
                        rep_id = add_new_node(model, stored);
                        rep_by_unit.emplace(rep_key, rep_id);
                        ++summary.nodes_added;
                        UnitGraphRecord ug;
                        ug.bubble_id = bubble.id;
                        ug.unit = stored;
                        ug.canonical_motif = canonical_motif_key(stored);
                        ug.anchor_pos = 0;
                        ug.anchor = {rep_id};
                        ug.nodes = {rep_id};
                        unit_graphs.push_back(std::move(ug));
                    } else {
                        rep_id = rit->second;
                    }
                    copy_steps.push_back(PathStep{rep_id, own != stored});
                }

                PathArray pa;
                pa.bubble_id = bubble.id;
                pa.start = arr.start;
                pa.length = arr.length;
                std::size_t next_copy = 0;
                for (const ArrayElement& el : arr.elements) {
                    if (el.is_copy) {
                        pa.replacement.push_back(copy_steps[next_copy++]);
                    } else {
                        pa.replacement.insert(pa.replacement.end(), el.steps.begin(), el.steps.end());
                    }
                }
                per_path_arrays[pi].push_back(std::move(pa));
                rewritten = true;

                for (const auto& r : arr.copy_ranges) {
                    for (std::size_t s = r.first; s < r.second; ++s) {
                        removal_candidates.insert(graph.paths[pi].steps[s].node_id);
                        collapsed_nodes.insert(graph.paths[pi].steps[s].node_id);
                    }
                    record_replaced_span(graph.paths[pi], bubble.id, r.first, r.second);
                }
                path_occ += arr.copies;
                if (unit_bp_report == 0) unit_bp_report = unit.size();
                interruptions_bp += arr.interruption_bp;
            }
            if (rewritten) {
                ++paths_norm;
                // Per haplotype, as in approximate mode: copies summed over all its arrays.
                if (min_copies_seen == 0 || path_occ < min_copies_seen) min_copies_seen = path_occ;
                if (path_occ > max_copies_seen) max_copies_seen = path_occ;
            }
        }

        const bool normalized = paths_norm > 0;
        if (normalized) {
            ++summary.bubbles_normalized;
        }
        // The prevalence's numerator: short motifs are not candidates, unless no motif is long enough, when the
        // most carried one is reported as the unit tried.
        std::size_t top_carriers = 0;
        std::size_t top_motif_bp = 0;
        std::string top_motif;
        for (const auto& [m, c] : carriers_by_motif)
            if ((short_motif_only || m.size() >= options.min_unit_bp) &&
                (c > top_carriers || (c == top_carriers && (m.size() > top_motif_bp ||
                                                            (m.size() == top_motif_bp && m < top_motif))))) {
                top_carriers = c;
                top_motif_bp = m.size();
                top_motif = m;
            }
        // A site left unfolded reports the unit it tried: its most carried motif, or for a short motif the
        // array unit made of it.
        if (!normalized) unit_bp_report = short_motif_only ? unit_bp_by_motif[top_motif] : top_motif_bp;
        if (short_motif_only && n_traverse > 0)
            prevalence = static_cast<double>(top_carriers) / static_cast<double>(n_traverse);
        const char* status = normalized ? "normalized"
                             : carriers_by_motif.empty() ? "no_tandem_detected"
                             : short_motif_only ? "short_motif"
                             : confirmed_motifs.empty() ? "below_prevalence"
                                                        : "no_copies_rewritten";
        ++summary.status_counts[status];
        const Bubble& seen = normalized ? bubble : input_bubble;   // a site widened but left unfolded keeps its own
        seen_sites.push_back({bubble.id, seen.source, seen.sink, normalized});
        claim_widened(bubble.id, normalized);
        {
            ReportRow row;
            std::ostringstream h, t;
            h << bubble.id << '\t' << (normalized ? "yes" : "no") << '\t'
              << unit_bp_report << '\t' << paths_norm << '\t'
              << min_copies_seen << '\t' << max_copies_seen << '\t'
              << interruptions_bp << '\t';
            t << n_traverse << '\t' << top_carriers << '\t' << std::to_string(prevalence) << '\t'   // as approximate mode writes it
              << carriers_by_motif.size() << '\t' << status << '\n';
            row.head = h.str();
            row.tail = t.str();
            row.collapsed.assign(collapsed_nodes.begin(), collapsed_nodes.end());
            report_rows.push_back(std::move(row));
        }
        progress.tick();
    }
    progress.done();

    // Splice replacements into each path's steps.
    for (std::size_t pi = 0; pi < model.paths.size(); ++pi) {
        auto& arrays = per_path_arrays[pi];
        if (arrays.empty()) {
            continue;
        }
        std::sort(arrays.begin(), arrays.end(),
                  [](const PathArray& a, const PathArray& b) { return a.start < b.start; });

        // Acceptance, structural and in every mode: the edits are ordered, non-overlapping and in bounds, and
        // the spliced result is exactly the original with those spans, and only those, substituted. It holds
        // under --lossy-units too, where copies may change sequence; whole-sequence equality is checked below.
        std::string expected_after;
        {
            const std::vector<PathStep>& old_steps = model.paths[pi].steps;
            const auto spell_steps = [&](const std::vector<PathStep>& st, std::size_t lo, std::size_t hi) {
                std::string acc;
                for (std::size_t k = lo; k < hi; ++k) {
                    const auto it = model.seq.find(st[k].node_id);
                    if (it == model.seq.end()) continue;
                    acc += st[k].reverse ? reverse_complement(it->second) : it->second;
                }
                return acc;
            };
            std::size_t prev_end = 0;
            for (const PathArray& a : arrays) {
                if (a.start < prev_end)
                    throw std::runtime_error("panphorte: overlapping edits on path " + graph.paths[pi].name);
                if (a.start + a.length > old_steps.size())
                    throw std::runtime_error("panphorte: edit runs past the end of path " +
                                             graph.paths[pi].name);
                if (a.length == 0 || a.replacement.empty())
                    throw std::runtime_error("panphorte: empty edit on path " + graph.paths[pi].name);
                expected_after += spell_steps(old_steps, prev_end, a.start);          // outside, verbatim
                expected_after += spell_steps(a.replacement, 0, a.replacement.size()); // the substitution
                prev_end = a.start + a.length;
            }
            expected_after += spell_steps(old_steps, prev_end, old_steps.size());
        }

        const std::vector<PathStep>& old_steps = model.paths[pi].steps;
        std::vector<PathStep> new_steps;
        new_steps.reserve(old_steps.size());
        std::size_t cursor = 0;
        for (const PathArray& a : arrays) {
            for (; cursor < a.start; ++cursor) {
                new_steps.push_back(old_steps[cursor]);
            }
            new_steps.insert(new_steps.end(), a.replacement.begin(), a.replacement.end());
            cursor = a.start + a.length;
        }
        for (; cursor < old_steps.size(); ++cursor) {
            new_steps.push_back(old_steps[cursor]);
        }
        model.paths[pi].steps = std::move(new_steps);
        {
            std::string got;
            for (const PathStep& st : model.paths[pi].steps) {
                const auto it = model.seq.find(st.node_id);
                if (it == model.seq.end()) continue;
                got += st.reverse ? reverse_complement(it->second) : it->second;
            }
            if (got != expected_after)
                throw std::runtime_error(
                    "panphorte: the rewrite of path " + graph.paths[pi].name +
                    " is not the original with the accepted copies substituted (" +
                    std::to_string(expected_after.size()) + " bp expected, " +
                    std::to_string(got.size()) + " bp produced)");
        }
        summary.paths_rewritten += 1;
    }

    // Ensure every consecutive step pair in every (rewritten) path has an L edge.
    EdgeSet edges(model);
    for (const GfaPath& p : model.paths) {
        for (std::size_t k = 1; k < p.steps.size(); ++k) {
            edges.ensure(p.steps[k - 1], p.steps[k]);
        }
    }
    summary.edges_added = edges.added;

    // Every rewritten path must spell its original sequence, in exact and approximate mode alike; with
    // --lossy-units only the reference is held to it.
    if (!approximate || options.lossless_units || ref_pi != SIZE_MAX) {
        for (std::size_t i = 0; i < model.paths.size(); ++i) {
            if (approximate && !options.lossless_units && i != ref_pi) continue;
            if (per_path_arrays[i].empty()) {
                continue;
            }
            const std::string original = spell_path_steps_sequence(graph, graph.paths[i].steps);
            if (spell_model_path(model, model.paths[i]) != original) {
                throw std::runtime_error(
                    "panphorte sequence-preservation violated for path index " + std::to_string(i));
            }
        }
    }

    // Remove collapsed nodes no longer referenced by any path; drop their edges.
    std::unordered_set<std::string> referenced;
    // The arcs the final paths walk, under the same canonical key as everything else, so a step pair
    // b-,a- is recognised as the stored a+>b+ rather than treated as a different link.
    std::unordered_set<std::string> traversed_edges;
    for (const GfaPath& p : model.paths) {
        for (const PathStep& s : p.steps) {
            referenced.insert(s.node_id);
        }
        for (std::size_t k = 1; k < p.steps.size(); ++k) {
            const PathStep& a = p.steps[k - 1];
            const PathStep& b = p.steps[k];
            traversed_edges.insert(
                EdgeSet::key(a.node_id, orient_char(a.reverse), b.node_id, orient_char(b.reverse)));
        }
    }
    std::unordered_set<std::string> to_remove;
    for (const std::string& id : removal_candidates) {
        if (referenced.find(id) == referenced.end()) {
            to_remove.insert(id);
        }
    }
    for (const ReportRow& row : report_rows) {
        std::size_t removed = 0;
        for (const std::string& n : row.collapsed) removed += to_remove.count(n);
        report << row.head << removed << '\t' << row.tail;
    }
    {
        std::vector<std::string> kept_order;
        kept_order.reserve(model.node_order.size());
        for (const std::string& id : model.node_order) {
            if (to_remove.find(id) == to_remove.end()) {
                kept_order.push_back(id);
            } else {
                model.seq.erase(id);
                ++summary.nodes_removed;
            }
        }
        model.node_order = std::move(kept_order);
        std::vector<GfaEdge> kept_edges;
        kept_edges.reserve(model.edges.size());
        for (const GfaEdge& e : model.edges) {
            if (to_remove.find(e.from) != to_remove.end() || to_remove.find(e.to) != to_remove.end()) {
                continue;
            }
            // Removing the nodes is not enough: a link keeps an old route walkable, and the re-snarl reads
            // links, not paths. A link is dropped only when no final path walks it, since links are global
            // (inside a repeat an arc is crossed once per copy); whether a leftover route matters is asked
            // per rewritten span below.
            const std::string k = EdgeSet::key(e.from, e.from_orient, e.to, e.to_orient);
            if (obsolete_edge_keys.find(k) != obsolete_edge_keys.end() &&
                traversed_edges.find(k) == traversed_edges.end()) {
                ++summary.edges_removed;
                continue;
            }
            kept_edges.push_back(e);
        }
        model.edges = std::move(kept_edges);
    }
    // Acceptance: is a route a path used before the rewrite still walkable after it? A surviving link alone
    // says nothing (at an array most survive because other paths need them); what matters is whether the
    // surviving nodes and arcs still spell the old allele end to end. Asked per edit, not per path: a
    // haplotype folded at two sites has two independent old routes, and one being gone would hide the other.
    // Each replaced span is checked with its arc in and its arc out.
    {
        std::unordered_set<std::string> live_nodes(model.node_order.begin(), model.node_order.end());
        std::unordered_set<std::string> live_edges;
        live_edges.reserve(model.edges.size() * 2);
        for (const GfaEdge& e : model.edges)
            live_edges.insert(EdgeSet::key(e.from, e.from_orient, e.to, e.to_orient));
        std::vector<std::string> survivors;
        std::set<std::size_t> survivor_sites;
        for (std::size_t i = 0; i < model.paths.size(); ++i) {
            const std::vector<PathStep>& was = graph.paths[i].steps;
            const auto step_live = [&](std::size_t k) {
                return live_nodes.find(was[k].node_id) != live_nodes.end();
            };
            const auto arc_live = [&](std::size_t k) {   // the arc from step k-1 to step k
                return live_edges.find(EdgeSet::key(was[k - 1].node_id, orient_char(was[k - 1].reverse),
                                                    was[k].node_id, orient_char(was[k].reverse))) !=
                       live_edges.end();
            };
            for (const PathArray& a : per_path_arrays[i]) {
                const std::size_t lo = a.start, hi = a.start + a.length;   // [lo, hi) was replaced
                bool walkable = true;
                for (std::size_t k = lo; walkable && k < hi; ++k) {
                    if (!step_live(k)) walkable = false;
                    else if (k > lo && !arc_live(k)) walkable = false;
                }
                // The anchors: without an arc in and an arc out the old interior is stranded, not a
                // route anything can walk.
                if (walkable && lo > 0 && !(step_live(lo - 1) && arc_live(lo))) walkable = false;
                if (walkable && hi < was.size() && !(step_live(hi) && arc_live(hi))) walkable = false;
                // Part of the old route is as bad as all of it when it joins the new one: at a step the
                // replacement keeps literal (a spacer, a partial copy), an old prefix can hand over to the folded
                // rest of the span, or the folded start to an old suffix (a contig ending inside the array keeps
                // the first copies' nodes and links).
                if (!walkable && lo > 0 && hi < was.size()) {
                    std::vector<PathStep> now;   // the new route between the same anchors
                    now.reserve(a.replacement.size() + 2);
                    now.push_back(was[lo - 1]);
                    now.insert(now.end(), a.replacement.begin(), a.replacement.end());
                    now.push_back(was[hi]);
                    const auto same_step = [](const PathStep& x, const PathStep& y) {
                        return x.node_id == y.node_id && x.reverse == y.reverse;
                    };
                    std::unordered_map<std::string, std::vector<std::size_t>> at;   // step -> positions in `now`
                    for (std::size_t q = 1; q + 1 < now.size(); ++q)
                        at[now[q].node_id + orient_char(now[q].reverse)].push_back(q);
                    // Old prefix was[lo-1 .. p] walkable, joining `now` at q: a hybrid unless it is `now`'s own prefix.
                    std::size_t lcp = 0;
                    while (lcp < now.size() && lo - 1 + lcp <= hi && same_step(was[lo - 1 + lcp], now[lcp])) ++lcp;
                    for (std::size_t p = lo; !walkable && p < hi && step_live(p) && arc_live(p); ++p) {
                        const auto it = at.find(was[p].node_id + orient_char(was[p].reverse));
                        if (it == at.end()) continue;
                        for (const std::size_t q : it->second)
                            if (!(q == p - (lo - 1) && q < lcp)) { walkable = true; break; }
                    }
                    // Old suffix was[p .. hi] walkable, entered from `now` at q.
                    std::size_t lcs = 0;
                    while (lcs < now.size() && hi >= lcs + lo - 1 &&
                           same_step(was[hi - lcs], now[now.size() - 1 - lcs])) ++lcs;
                    for (std::size_t p = hi - 1; !walkable && p >= lo && step_live(p) && arc_live(p + 1); --p) {
                        const auto it = at.find(was[p].node_id + orient_char(was[p].reverse));
                        if (it != at.end())
                            for (const std::size_t q : it->second)
                                if (!(hi - p == now.size() - 1 - q && hi - p < lcs)) { walkable = true; break; }
                        if (p == lo) break;
                    }
                }
                if (walkable) {
                    ++summary.routes_surviving;
                    survivor_sites.insert(a.bubble_id);
                    if (survivors.size() < 4)
                        survivors.push_back(graph.paths[i].name + " (bubble " +
                                            std::to_string(a.bubble_id) + ")");
                }
            }
        }
        // A graph carrying the same site both folded and unfolded is wrong for every consumer. It arises where
        // other walks keep the replaced span (contigs ending inside the site, haplotypes whose partial copies
        // sit on the same nodes), and separating the two would need the shared context cloned. Those sites are
        // declined (panphorte_normalize runs again without them) rather than the whole run refused.
        if (summary.routes_surviving > 0) {
            throw SharedRouteSites(survivor_sites,
                "panphorte: " + std::to_string(summary.routes_surviving) +
                " rewritten path(s) can still be walked along their ORIGINAL route in the normalized "
                "graph (" + cli::join_with_comma(survivors) +
                "), so the site is represented twice -- once folded and once as the branch it replaced, "
                "which the re-snarl will pick up as a second allele. This happens when every node and "
                "link of the replaced span is also used somewhere else, and separating the two would "
                "need the reused context cloned before the rewrite, which is not implemented");
        }
    }

    say("folded: " + std::to_string(summary.bubbles_normalized) + " of " + std::to_string(summary.bubbles_seen) +
        " sites; paths rewritten: " + std::to_string(summary.paths_rewritten));
    {
        std::string not_folded;
        for (const auto& [status, n] : summary.status_counts)
            if (status != "normalized") not_folded += (not_folded.empty() ? "" : "; ") + status + ": " + std::to_string(n);
        if (!not_folded.empty()) say("not folded: " + not_folded);
    }
    say("nodes: " + std::to_string(model.node_order.size()) + "; added: " + std::to_string(summary.nodes_added) +
        " (fragments: " + std::to_string(summary.fragment_nodes_added) + "); removed: " +
        std::to_string(summary.nodes_removed));
    say("edges: " + std::to_string(model.edges.size()) + "; added: " + std::to_string(summary.edges_added) +
        "; removed: " + std::to_string(summary.edges_removed));

    // --refine: realign the residue around each folded site (see panphorte_refine.hpp). Runs on the
    // graph exactly as folded and checked above, before the sort and the re-snarl.
    if (options.refine) {
        const std::string resolved = resolve_reference_path_name(graph, options.reference_path, "panphorte");
        std::unordered_map<std::size_t, std::unordered_set<std::string>> spliced;
        for (const auto& arrays : per_path_arrays)
            for (const PathArray& a : arrays)
                for (const PathStep& st : a.replacement) spliced[a.bubble_id].insert(st.node_id);
        std::unordered_map<std::size_t, const Bubble*> by_id;
        for (const Bubble& b : bubbles) by_id[b.id] = &b;
        std::vector<PanphorteRefineSite> sites;
        for (const SeenSite& s : seen_sites) {
            if (!s.folded) continue;
            PanphorteRefineSite site;
            site.bubble_id = s.bubble_id;
            site.source = s.source;
            site.sink = s.sink;
            for (const std::string& n : by_id.at(s.bubble_id)->inside)
                if (model.has_node(n)) site.interior.insert(n);
            for (const std::string& n : spliced[s.bubble_id])
                if (model.has_node(n)) site.interior.insert(n);
            site.interior.erase(s.source);
            site.interior.erase(s.sink);
            sites.push_back(std::move(site));
        }
        PanphorteRefineLimits lim;
        lim.max_poa_bp = options.refine_max_poa_bp;
        lim.max_walks = options.refine_max_walks;
        lim.max_poa_work = options.refine_max_poa_work;
        std::unordered_set<std::string> rep_nodes;
        for (const UnitGraphRecord& u : unit_graphs) rep_nodes.insert(u.nodes.begin(), u.nodes.end());
        const PanphorteRefineResult rr = panphorte_refine_sites(model, sites, resolved, lim, rep_nodes);
        summary.refine_refined = rr.refined;
        summary.refine_kept = rr.kept;
        summary.refine_skipped = rr.skipped;
        summary.refine_nodes_added = rr.nodes_added;
        summary.refine_nodes_removed = rr.nodes_removed;
        summary.edges_added += rr.edges_added;
        summary.edges_removed += rr.edges_removed;
        say("refine: realigned: " + std::to_string(rr.refined) + "; kept as folded: " + std::to_string(rr.kept) +
            "; skipped: " + std::to_string(rr.skipped) + "; nodes added: " + std::to_string(rr.nodes_added) +
            "; nodes removed: " + std::to_string(rr.nodes_removed));

        std::ofstream rt(stage_out(options.out_prefix + ".panphorte.refine.tsv"));
        if (!rt) throw std::runtime_error("panphorte: cannot write the refine table");
        // Ids are the folded graph's (input numbering), as in the copies table: the sort renumbers nodes.
        rt << "input_bubble_id\tinput_source\tinput_sink\trefine_status\tdetail\tpaths\tcarriers\t"
              "distinct_flanks_before\tdistinct_flanks_after\tnodes_before\tnodes_after\tsteps_before\t"
              "steps_after\tsingle_path_nodes_before\tsingle_path_nodes_after\n";
        std::unordered_map<std::size_t, const PanphorteRefineRow*> row_of;
        for (const PanphorteRefineRow& r : rr.rows) row_of[r.bubble_id] = &r;
        for (const SeenSite& s : seen_sites) {
            rt << s.bubble_id << '\t' << s.source << '\t' << s.sink << '\t';
            const auto it = row_of.find(s.bubble_id);
            if (it == row_of.end()) {
                rt << "skipped:not_folded\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\n";
                continue;
            }
            const PanphorteRefineRow& r = *it->second;
            const auto num = [&](std::size_t v) { return r.measured ? std::to_string(v) : std::string("-"); };
            rt << r.status << '\t' << (r.detail.empty() ? "-" : r.detail) << '\t' << r.paths << '\t'
               << r.carriers << '\t' << r.distinct_pre << '\t' << r.distinct_post << '\t'
               << num(r.nodes_before) << '\t' << num(r.nodes_after) << '\t' << num(r.steps_before) << '\t'
               << num(r.steps_after) << '\t' << num(r.single_path_nodes_before) << '\t'
               << num(r.single_path_nodes_after) << '\n';
        }
        rt.close();
        if (!rt) throw std::runtime_error("panphorte: failed to write the refine table");
    }

    // old unit node id -> the id it is delivered under. Empty means the ids did not move.
    std::unordered_map<std::string, std::string> rep_id_remap;

    // With a reference the output is made call-ready with no external tools: sorted and flipped along the
    // reference (placing the new unit nodes in reference order) and re-snarled with the cactus finder,
    // writing the sorted GFA, a bubbles CSV and Bandage colours. Without one, the unsorted GFA is written.
    if (options.reference_path.empty()) {
        write_gfa_model(stage_out(options.out_prefix + ".normalized.gfa"), model);
    } else {
        GraphSortOptions sort_opts;
        sort_opts.reference_path = options.reference_path;
        sort_opts.flip = !options.no_flip;
        // The sorter resolves an exact name, a unique case-insensitive match or a unique substring;
        // everything after this compares path names exactly, so it uses the resolved name.
        const GraphSortResult sort_result = sort_graph_reference(model, sort_opts);
        const std::string resolved_reference = sort_result.resolved_reference;
        rep_id_remap = sort_result.id_remap;
        say("sorted along the reference: nodes turned around: " + std::to_string(sort_result.flipped_nodes) +
            "; loops placed as one block: " + std::to_string(sort_result.cycle_blocks));

        const std::string sorted_gfa = stage_out(options.out_prefix + ".normalized.sorted.gfa");
        write_gfa_model(sorted_gfa, model);

        BubbleCallOptions bopts;
        bopts.reference_path = resolved_reference;
        // The re-snarl's size filter is panphorte's own option, since the input CSV may have been made
        // under a different threshold; a change in the bubble count is reported below.
        bopts.min_variant_bp = options.resnarl_min_variant_bp;
        bopts.merge_nearby_bp = options.resnarl_merge_nearby_bp;
        bopts.fuse_ambiguous = options.resnarl_fuse_ambiguous;
        bopts.snarl_pairs_override = find_top_level_snarls_cactus(snarl_input_from_model(model), &bopts.snarl_inward_override);
    bopts.snarl_source_supplied = true;
        bopts.quiet = options.quiet;
        bopts.log_module = "panphorte";

        ParseGfaOptions parse_options;
        parse_options.include_paths = true;
        parse_options.include_sequences = true;
        const Graph sorted_graph = parse_gfa(sorted_gfa, parse_options);
        const BubbleCallReport report = call_bubbles_report(sorted_graph, bopts);
        write_bubbles_csv(stage_out(options.out_prefix + ".bubbles.csv"), report.bubbles);
        // Inversion sites are left out of the folded set (read_bubbles_csv), so out of this count too. The
        // re-snarl reassigns ids, so they do not match the input's.
        say("sites after folding: " +
            std::to_string(std::count_if(report.bubbles.begin(), report.bubbles.end(),
                                         [](const Bubble& b) { return !b.inversion_parent; })) +
            " (ids reassigned)");
        write_bandage_node_colors_csv(stage_out(options.out_prefix + ".bandage_nodes.csv"),
                                      report.bubbles, report.non_snp_bubbles);

        // Gene annotation runs inside the transaction, so a malformed GTF fails before any output lands.
        if (!options.gtf_path.empty()) {
            summary.genes_written = emit_gene_annotation(
                sorted_graph, resolved_reference, options.gtf_path,
                stage_out(options.out_prefix + ".bandage_genes.csv"), "", "panphorte");
        }

        summary.sorted = true;
        summary.resnarled_bubbles = report.bubbles.size();
    }

    // The unit graph of every fold, ids as delivered. `unit` is the unit as folded, `canonical_motif` its
    // primitive motif (the min over both strands, so the same whichever way the sort reads the nodes) and
    // `copy_quantum` the motif copies one copy of the unit stands for. `anchor_pos` is a unit position
    // (1-based) no copy deletes and `anchor_nodes` the nodes covering it: a haplotype's copies are its passes
    // through them (one per copy). `nodes` are every node of the unit graph; a one-node unit graph is the unit
    // looped once per copy, its own anchor.
    {
        std::ofstream ug(stage_out(options.out_prefix + ".panphorte.unit_graph.tsv"));
        if (!ug) throw std::runtime_error("panphorte: cannot write the unit graph table");
        ug << "input_bubble_id\tunit_bp\tcanonical_motif\tcopy_quantum\tunit\tanchor_pos\tanchor_nodes\tn_nodes\tnodes\n";
        auto out_id = [&](const std::string& id) {
            const auto it = rep_id_remap.find(id);
            return it == rep_id_remap.end() ? id : it->second;
        };
        auto upper = [](std::string x) { for (char& c : x) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c))); return x; };
        for (const UnitGraphRecord& u : unit_graphs) {
            std::string ids;
            for (const std::string& n : u.nodes) {
                if (model.seq.find(out_id(n)) == model.seq.end())
                    throw std::runtime_error("panphorte: unit graph node " + n + " (bubble " +
                                             std::to_string(u.bubble_id) + ") is missing from the normalized graph");
                ids += (ids.empty() ? "" : ";") + out_id(n);
            }
            // A one-node unit graph must still spell its unit (on either strand: the sort may flip it), or
            // this table would name a node that is no longer the unit.
            if (u.nodes.size() == 1) {
                const std::string got = upper(model.seq.at(out_id(u.nodes.front()))), want = upper(u.unit);
                if (got != want && got != upper(reverse_complement(u.unit)))
                    throw std::runtime_error("panphorte: unit node " + out_id(u.nodes.front()) + " (bubble " +
                                             std::to_string(u.bubble_id) + ") no longer spells its unit");
            }
            std::string anchor;
            for (const std::string& n : u.anchor) anchor += (anchor.empty() ? "" : ";") + out_id(n);
            const std::size_t quantum =
                u.canonical_motif.empty() ? 1 : std::max<std::size_t>(1, u.unit.size() / u.canonical_motif.size());
            ug << u.bubble_id << '\t' << u.unit.size() << '\t' << u.canonical_motif << '\t' << quantum << '\t'
               << u.unit << '\t' << (anchor.empty() ? "." : std::to_string(u.anchor_pos + 1)) << '\t'
               << (anchor.empty() ? "." : anchor) << '\t' << u.nodes.size() << '\t' << ids << '\n';
        }
        ug.close();
        if (!ug) throw std::runtime_error("panphorte: failed to write the unit graph table");
    }

    // Everything succeeded: move the family into place. Closed and checked, not just flushed -- a
    // failure that surfaces on close would otherwise be committed as a complete file.
    report.close();
    if (!report) throw std::runtime_error("panphorte: failed to write the report table");
    if (copies_out.is_open()) {
        copies_out.close();
        if (!copies_out) throw std::runtime_error("panphorte: failed to write the copies table");
    }
    staged.commit();

    // The family's default-named files this run did not write were left by an earlier run with other
    // options, and `call -b <prefix>` would read a stale site table beside this run's graph. Removed,
    // unless one is an input.
    for (const char* suffix : {".normalized.gfa", ".normalized.sorted.gfa", ".bubbles.csv", ".bandage_nodes.csv",
                               ".bandage_genes.csv", ".panphorte.copies.tsv", ".panphorte.report.tsv",
                               ".panphorte.refine.tsv", ".panphorte.unit_graph.tsv",
                               // not written: a leftover, whose rows the unit graph table holds
                               ".panphorte.rep_provenance.tsv"}) {
        const std::string path = options.out_prefix + suffix;
        if (written.count(path) != 0) continue;
        std::error_code ec;
        bool is_input = false;
        for (const std::string& in : {options.gfa_path, options.bubbles_csv_in, options.gtf_path})
            is_input = is_input || (!in.empty() && std::filesystem::equivalent(path, in, ec));
        if (!is_input && std::filesystem::remove(path, ec) && !options.quiet)
            std::cerr << cli::log_prefix("panphorte") << "removed " << path << ", left by an earlier run at this prefix\n";
    }

    if (summary_out != nullptr) {
        *summary_out = summary;
    }
}

} // namespace

void panphorte_normalize(const PanphorteOptions& options, PanphorteSummary* summary_out) {
    std::set<std::size_t> declined;
    while (true) {
        try {
            normalize_pass(options, summary_out, declined);
            return;
        } catch (const SharedRouteSites& e) {
            std::size_t added = 0;
            for (const std::size_t id : e.bubble_ids) added += declined.insert(id).second;
            if (added == 0) throw std::runtime_error(e.what());   // nothing left to decline: refuse
            if (!options.quiet) {
                std::cerr << cli::log_prefix("panphorte") << "declined_shared_route: " << added
                          << " (sites whose original route other walks still use; folding again without them):";
                for (const std::size_t id : e.bubble_ids) std::cerr << ' ' << id;
                std::cerr << '\n';
            }
        }
    }
}

} // namespace panvar
