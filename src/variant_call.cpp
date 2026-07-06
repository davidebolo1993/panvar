#include "panvar/variant_call.hpp"

#include "panvar/align.hpp"
#include "panvar/bubble_path.hpp"
#include "panvar/bubbles.hpp"
#include "panvar/cli_utils.hpp"
#include "panvar/gene_cn_kmer.hpp"
#include "panvar/graph_utils.hpp"
#include "panvar/gtf.hpp"
#include "panvar/minimap2_align.hpp"
#include "panvar/output.hpp"
#include "panvar/parallel.hpp"
#include "panvar/ref_path.hpp"

#include <algorithm>
#include <cctype>
#include <cstddef>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <mutex>
#include <numeric>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace panvar {
namespace {

constexpr std::size_t kAlignCellCap = 4000000; // skip walk-diff if |ref|*|hap| exceeds this

enum class EvType { Del, Ins, Inv, Dup };

const char* ev_svtype(EvType t) {
    switch (t) {
        case EvType::Del: return "DEL";
        case EvType::Ins: return "INS";
        case EvType::Inv: return "INV";
        case EvType::Dup: return "DUP";
    }
    return "DEL";
}

std::size_t node_len(const Graph& graph, const std::string& id) {
    const auto it = graph.nodes.find(id);
    return it == graph.nodes.end() ? 0 : it->second.sequence.size();
}

// Nodes carrying a self-loop edge (a panphorte REP, or any tandem-unit node):
// these are copy-number loci even when a haplotype traverses them only once.
std::unordered_set<std::string> self_loop_nodes(const Graph& graph) {
    std::unordered_set<std::string> out;
    for (const auto& [id, node] : graph.nodes) {
        bool loop = false;
        for (const Neighbor& nb : node.start) if (nb.node_id == id) { loop = true; break; }
        if (!loop) for (const Neighbor& nb : node.end) if (nb.node_id == id) { loop = true; break; }
        if (loop) out.insert(id);
    }
    return out;
}

// Steps of `path` across `bubble` (canonical source->sink). Falls back to an empty interior
// ([source, sink]) for paths that cross with no inside node (a pure deletion / short side of an
// insertion), which the inside-node-only interval finder would otherwise drop.
std::optional<std::vector<PathStep>> bubble_steps(
    const PathRecord& path, const BubblePathIndex& index, const Bubble& bubble) {
    const auto iv = find_best_bubble_path_interval(index, bubble);
    if (iv.has_value()) {
        std::vector<PathStep> s = canonical_bubble_path_steps(path, bubble, *iv);
        if (!s.empty()) return s;
    }
    const auto si = index.positions.find(bubble.source);
    const auto ki = index.positions.find(bubble.sink);
    if (si == index.positions.end() || ki == index.positions.end()) return std::nullopt;
    const std::unordered_set<std::size_t> sink_pos(ki->second.begin(), ki->second.end());
    for (const std::size_t p : si->second) {                 // forward: source then sink
        if (sink_pos.count(p + 1)) return std::vector<PathStep>{ path.steps[p], path.steps[p + 1] };
    }
    const std::unordered_set<std::size_t> src_pos(si->second.begin(), si->second.end());
    for (const std::size_t p : ki->second) {                 // reverse: sink then source -> flip
        if (src_pos.count(p + 1)) {
            return std::vector<PathStep>{
                PathStep{ path.steps[p + 1].node_id, !path.steps[p + 1].reverse },
                PathStep{ path.steps[p].node_id, !path.steps[p].reverse } };
        }
    }
    return std::nullopt;
}

// One node-token position in a walk (oriented), for the DEL/INS/INV alignment.
struct Tok {
    std::uint64_t token = 0;
    std::string node_id;
    bool reverse = false;
};

// Token walk for the DEL/INS/INV alignment. Copy-number nodes (REP self-loops)
// are dropped entirely — they are handled separately as count-based DUP/CN events,
// so a haplotype's extra REP copies are never mistyped as INS.
std::vector<Tok> collapse_walk(
    const std::vector<PathStep>& steps,
    const std::unordered_set<std::string>& cn_nodes) {

    std::vector<Tok> out;
    for (const PathStep& s : steps) {
        if (cn_nodes.count(s.node_id) != 0) {
            continue;
        }
        Tok t;
        t.token = hash_step_token(s);
        t.node_id = s.node_id;
        t.reverse = s.reverse;
        out.push_back(std::move(t));
    }
    return out;
}

// A typed event derived from one haplotype walk vs the reference walk.
struct Event {
    EvType type = EvType::Del;
    std::vector<std::string> nodes;   // variant node set (length-weighted Jaccard key)
    std::string start_node;           // graph anchors (first/last involved node)
    std::string end_node;
    std::string seq;                  // event sequence (DELSEQ/INSSEQ/INVSEQ)
    std::size_t ref_cn = 0;           // DUP only
    std::size_t alt_cn = 0;           // DUP only
    std::string ins_subtype;          // INS only: "", "NOVEL", "DUP"
    std::string link_id;              // shared id for a co-located DEL+INS substitution (EVENTID)
    // Reference anchoring for VCF coordinates.
    std::string anchor_node;          // ref node POS is taken from
    bool anchor_after = false;        // true: POS = last base of anchor (INS); false: first base
    std::size_t size_bp = 0;          // |event| for min_sv filtering / SVLEN magnitude
    std::size_t ref_pos = 0;          // reference genomic position of the anchor (merge window)
    bool cn_peak = false;             // DUP from peak multiplicity (size_bp is the duplicated bp)
    std::size_t ru_len = 0;           // DUP only: repeat-unit length in bp (RU_LEN; one copy)
};

// Spell a token run into sequence. Builds in place with a single reservation -- no throwaway
// vector<PathStep> and no per-node reverse_complement temporary -- because this is the caller's
// hot path (one call per gap block per allele) and the allocator becomes the bottleneck under the
// parallel allele loop. Reverse bytes use the same complement mapping as reverse_complement().
std::string spell_toks(const Graph& graph, const std::vector<const Tok*>& toks) {
    std::size_t total = 0;
    for (const Tok* t : toks) {
        const auto it = graph.nodes.find(t->node_id);
        if (it == graph.nodes.end() || it->second.sequence.empty()) {
            // Defer to the canonical speller so the missing-node error is identical.
            std::vector<PathStep> steps;
            for (const Tok* u : toks) steps.push_back(PathStep{u->node_id, u->reverse});
            return spell_path_steps_sequence(graph, steps);
        }
        total += it->second.sequence.size();
    }
    std::string seq;
    seq.reserve(total);
    for (const Tok* t : toks) {
        const std::string& ns = graph.nodes.at(t->node_id).sequence;
        if (!t->reverse) { seq += ns; continue; }
        for (auto it = ns.rbegin(); it != ns.rend(); ++it) {
            switch (*it) {
                case 'A': case 'a': seq.push_back('T'); break;
                case 'C': case 'c': seq.push_back('G'); break;
                case 'G': case 'g': seq.push_back('C'); break;
                case 'T': case 't': seq.push_back('A'); break;
                default: seq.push_back('N'); break;
            }
        }
    }
    return seq;
}

std::size_t toks_bp(const Graph& graph, const std::vector<const Tok*>& toks) {
    std::size_t bp = 0;
    for (const Tok* t : toks) bp += node_len(graph, t->node_id);
    return bp;
}

// Is rc(hap block) == ref block (an inversion)? Compares node-id + orientation
// of the reference block against the reverse-complemented haplotype block.
bool is_inversion(const std::vector<const Tok*>& ref_blk, const std::vector<const Tok*>& hap_blk) {
    if (ref_blk.size() != hap_blk.size() || ref_blk.empty()) return false;
    for (std::size_t i = 0; i < ref_blk.size(); ++i) {
        const Tok* r = ref_blk[i];
        const Tok* h = hap_blk[hap_blk.size() - 1 - i]; // reverse order
        if (r->node_id != h->node_id || r->reverse == h->reverse) return false; // orientation must flip
    }
    return true;
}

// Align one sub-range R[r0,r1) vs H[h0,h1) with the node-token DP (diagonal only on
// equal tokens; else gaps) and append the DEL/INS/INV events to `events`.
// `preceding_ref_node` anchors an INS that opens the segment (the last matched ref
// node before it). Bounded by kAlignCellCap; segments between shared anchors are small.
void diff_segment(
    const Graph& graph,
    const std::vector<Tok>& R, std::size_t r0, std::size_t r1,
    const std::vector<Tok>& H, std::size_t h0, std::size_t h1,
    const std::string& preceding_ref_node,
    std::vector<Event>& events) {

    const std::size_t m = r1 - r0;
    const std::size_t n = h1 - h0;
    if (m == 0 && n == 0) return;
    if (m * n > kAlignCellCap) {
        if (std::getenv("PANVAR_CALL_DEBUG")) {
            std::cerr << "[diff] SKIP segment m=" << m << " n=" << n << " (cap " << kAlignCellCap << ")\n";
        }
        return; // unanchored divergent block too large to align; skip
    }

    // DP: match (+2) only when tokens equal, gaps (-1). Mismatch diagonals forbidden
    // so substitutions surface as separate ref-only / hap-only blocks (for INV check).
    const int kMatch = 2;
    const int kGap = -1;
    const int kNeg = -1000000;
    // Flat (m+1)x(n+1) DP buffer: one allocation instead of (m+1) row vectors. The full matrix is
    // retained because the traceback below revisits arbitrary cells.
    const std::size_t W = n + 1;
    std::vector<int> dp((m + 1) * W, 0);
    auto D = [&](std::size_t i, std::size_t j) -> int& { return dp[i * W + j]; };
    for (std::size_t i = 1; i <= m; ++i) D(i, 0) = static_cast<int>(i) * kGap;
    for (std::size_t j = 1; j <= n; ++j) D(0, j) = static_cast<int>(j) * kGap;
    for (std::size_t i = 1; i <= m; ++i) {
        for (std::size_t j = 1; j <= n; ++j) {
            const int diag = (R[r0 + i - 1].token == H[h0 + j - 1].token) ? D(i - 1, j - 1) + kMatch : kNeg;
            const int up = D(i - 1, j) + kGap;
            const int left = D(i, j - 1) + kGap;
            D(i, j) = std::max(diag, std::max(up, left));
        }
    }

    // Traceback into columns (ri, hi) into R/H absolute indices; -1 means a gap on that side.
    struct Col { long long ri; long long hi; };
    std::vector<Col> cols;
    std::size_t i = m, j = n;
    while (i > 0 || j > 0) {
        if (i > 0 && j > 0 && R[r0 + i - 1].token == H[h0 + j - 1].token &&
            D(i, j) == D(i - 1, j - 1) + kMatch) {
            cols.push_back({static_cast<long long>(r0 + i - 1), static_cast<long long>(h0 + j - 1)});
            --i; --j;
        } else if (i > 0 && D(i, j) == D(i - 1, j) + kGap) {
            cols.push_back({static_cast<long long>(r0 + i - 1), -1});
            --i;
        } else {
            cols.push_back({-1, static_cast<long long>(h0 + j - 1)});
            --j;
        }
    }
    std::reverse(cols.begin(), cols.end());

    // Emit events from maximal gap blocks -> DEL / INS / INV. Track the last matched
    // ref node as the anchor for an INS that follows it.
    std::string last_ref_node = preceding_ref_node;
    auto flush_block = [&](std::vector<const Tok*>& ref_blk, std::vector<const Tok*>& hap_blk) {
        if (ref_blk.empty() && hap_blk.empty()) return;
        if (!ref_blk.empty() && !hap_blk.empty() && is_inversion(ref_blk, hap_blk)) {
            Event e;
            e.type = EvType::Inv;
            for (const Tok* t : ref_blk) e.nodes.push_back(t->node_id);
            e.start_node = ref_blk.front()->node_id;
            e.end_node = ref_blk.back()->node_id;
            e.seq = spell_toks(graph, ref_blk);
            e.size_bp = toks_bp(graph, ref_blk);
            e.anchor_node = ref_blk.front()->node_id;
            events.push_back(std::move(e));
        } else {
            // A gap-block with both ref and hap content is a substitution: emit DEL+INS
            // but link them with a shared EVENTID so downstream can pair them.
            const bool substitution = !ref_blk.empty() && !hap_blk.empty();
            const std::string link =
                substitution ? ("sub_" + ref_blk.front()->node_id + "_" + hap_blk.front()->node_id) : "";
            if (!ref_blk.empty()) {
                Event e;
                e.type = EvType::Del;
                for (const Tok* t : ref_blk) e.nodes.push_back(t->node_id);
                e.start_node = ref_blk.front()->node_id;
                e.end_node = ref_blk.back()->node_id;
                e.seq = spell_toks(graph, ref_blk);
                e.size_bp = toks_bp(graph, ref_blk);
                e.anchor_node = ref_blk.front()->node_id;
                e.link_id = link;
                events.push_back(std::move(e));
            }
            if (!hap_blk.empty()) {
                Event e;
                e.type = EvType::Ins;
                for (const Tok* t : hap_blk) e.nodes.push_back(t->node_id);
                e.start_node = hap_blk.front()->node_id;
                e.end_node = hap_blk.back()->node_id;
                e.seq = spell_toks(graph, hap_blk);
                e.size_bp = toks_bp(graph, hap_blk);
                e.anchor_node = last_ref_node;
                e.anchor_after = true;
                e.link_id = link;
                events.push_back(std::move(e));
            }
        }
        ref_blk.clear();
        hap_blk.clear();
    };

    std::vector<const Tok*> ref_blk;
    std::vector<const Tok*> hap_blk;
    for (const Col& c : cols) {
        if (c.ri >= 0 && c.hi >= 0) {
            flush_block(ref_blk, hap_blk);
            last_ref_node = R[static_cast<std::size_t>(c.ri)].node_id;
        } else if (c.ri >= 0) {
            ref_blk.push_back(&R[static_cast<std::size_t>(c.ri)]);
        } else {
            hap_blk.push_back(&H[static_cast<std::size_t>(c.hi)]);
        }
    }
    flush_block(ref_blk, hap_blk);
}

// DEL/INS/INV of haplotype vs reference. Split both walks at shared anchors (tokens unique in both,
// chained in monotonic order) so the quadratic DP runs only between anchors - bounds cost on huge
// bubbles and gives identical breakpoints across haplotypes so the same event merges. CN is separate.
std::vector<Event> diff_walks(
    const Graph& graph,
    const std::vector<Tok>& R,
    const std::vector<Tok>& H,
    const std::string& bubble_source) {

    std::vector<Event> events;
    const std::size_t m = R.size();
    const std::size_t n = H.size();

    // Tokens unique within R and within H are anchor candidates.
    std::unordered_map<std::uint64_t, std::size_t> r_count, h_count;
    std::unordered_map<std::uint64_t, std::size_t> r_pos, h_pos;
    for (std::size_t i = 0; i < m; ++i) { ++r_count[R[i].token]; r_pos[R[i].token] = i; }
    for (std::size_t j = 0; j < n; ++j) { ++h_count[H[j].token]; h_pos[H[j].token] = j; }

    // Candidate anchors ordered by reference index; chain a monotonic subsequence on
    // the haplotype index (LIS) so anchors are consistent in both walks.
    std::vector<std::pair<std::size_t, std::size_t>> cand; // (ref_idx, hap_idx)
    for (std::size_t i = 0; i < m; ++i) {
        const std::uint64_t tok = R[i].token;
        if (r_count[tok] == 1) {
            const auto it = h_count.find(tok);
            if (it != h_count.end() && it->second == 1) {
                cand.emplace_back(i, h_pos[tok]);
            }
        }
    }
    // LIS on hap index (strictly increasing) over cand (already sorted by ref index).
    std::vector<std::size_t> chain_ref, chain_hap;
    {
        std::vector<std::size_t> tails;       // tails[k] = index into cand of smallest tail of an LIS of length k+1
        std::vector<long long> prev(cand.size(), -1);
        std::vector<std::size_t> tail_at;     // parallel hap value for binary search
        for (std::size_t c = 0; c < cand.size(); ++c) {
            const std::size_t hv = cand[c].second;
            auto pos = std::lower_bound(tail_at.begin(), tail_at.end(), hv);
            const std::size_t k = static_cast<std::size_t>(pos - tail_at.begin());
            if (k > 0) prev[c] = static_cast<long long>(tails[k - 1]);
            if (pos == tail_at.end()) { tails.push_back(c); tail_at.push_back(hv); }
            else { tails[k] = c; tail_at[k] = hv; }
        }
        if (!tails.empty()) {
            long long cur = static_cast<long long>(tails.back());
            while (cur >= 0) {
                chain_ref.push_back(cand[static_cast<std::size_t>(cur)].first);
                chain_hap.push_back(cand[static_cast<std::size_t>(cur)].second);
                cur = prev[static_cast<std::size_t>(cur)];
            }
            std::reverse(chain_ref.begin(), chain_ref.end());
            std::reverse(chain_hap.begin(), chain_hap.end());
        }
    }

    if (std::getenv("PANVAR_CALL_DEBUG")) {
        std::cerr << "[diff] m=" << m << " n=" << n << " anchors=" << chain_ref.size() << "\n";
    }

    // Walk the anchor chain, aligning each inter-anchor segment independently. The
    // anchor itself is a match (resets the INS anchor to that ref node).
    std::size_t r_prev = 0, h_prev = 0;
    std::string preceding_ref = bubble_source;
    for (std::size_t a = 0; a <= chain_ref.size(); ++a) {
        const std::size_t r_end = (a < chain_ref.size()) ? chain_ref[a] : m;
        const std::size_t h_end = (a < chain_hap.size()) ? chain_hap[a] : n;
        diff_segment(graph, R, r_prev, r_end, H, h_prev, h_end, preceding_ref, events);
        if (a < chain_ref.size()) {
            preceding_ref = R[r_end].node_id; // the anchor node
            r_prev = r_end + 1;
            h_prev = h_end + 1;
        }
    }
    return events;
}

// Within a haplotype, coalesce consecutive same-type DEL/INS/INV events whose
// reference gap is <= merge_distance_bp into one event (union node sets). The gap
// is measured in reference bp via the bubble's node->genomic-start map.
void coalesce_events(
    const Graph& graph,
    std::vector<Event>& events,
    const std::unordered_map<std::string, std::size_t>& ref_node_pos,
    const std::unordered_map<std::string, std::size_t>& hap_node_pos,
    std::size_t merge_distance_bp) {

    if (events.size() < 2) return;
    auto pos_in = [](const std::unordered_map<std::string, std::size_t>& m,
                     const std::string& node) -> long long {
        const auto it = m.find(node);
        return it == m.end() ? -1 : static_cast<long long>(it->second);
    };
    // Reference bp interval [lo, hi) that an event occupies on the reference.
    auto ref_span = [&](const Event& e, long long& lo, long long& hi) {
        if (e.type == EvType::Ins) {
            const long long p = pos_in(ref_node_pos, e.anchor_node);
            const long long end = (p < 0) ? p : p + static_cast<long long>(node_len(graph, e.anchor_node));
            lo = end; hi = end;
            return;
        }
        lo = -1; hi = -1;
        for (const std::string& n : e.nodes) {
            const long long p = pos_in(ref_node_pos, n);
            if (p < 0) continue;
            const long long q = p + static_cast<long long>(node_len(graph, n));
            if (lo < 0 || p < lo) lo = p;
            if (hi < 0 || q > hi) hi = q;
        }
    };
    // Haplotype bp interval the event's own (inserted) nodes occupy in this haplotype's
    // walk. Well-defined for INS (its nodes are on the haplotype); for DEL/INV the nodes
    // are reference-only and absent from hap_node_pos, so this yields no span and the
    // merge falls back to the reference metric. This catches same-type events that are
    // far apart on the reference (a deletion sits between them) yet contiguous in the
    // sample's own sequence.
    auto hap_span = [&](const Event& e, long long& lo, long long& hi) {
        lo = -1; hi = -1;
        for (const std::string& n : e.nodes) {
            const long long p = pos_in(hap_node_pos, n);
            if (p < 0) continue;
            const long long q = p + static_cast<long long>(node_len(graph, n));
            if (lo < 0 || p < lo) lo = p;
            if (hi < 0 || q > hi) hi = q;
        }
    };
    const long long d = static_cast<long long>(merge_distance_bp);
    auto gap_ok = [&](long long prev_hi, long long cur_lo) {
        return prev_hi >= 0 && cur_lo >= 0 && cur_lo - prev_hi <= d && cur_lo - prev_hi >= -d;
    };
    std::vector<Event> out;
    long long prev_hi = -1, prev_hhi = -1;
    for (Event& e : events) {
        long long lo = -1, hi = -1, hlo = -1, hhi = -1;
        ref_span(e, lo, hi);
        hap_span(e, hlo, hhi);
        // Merge when the same-type predecessor is near EITHER in reference space OR in this
        // haplotype's own sequence space.
        if (e.type != EvType::Dup && !out.empty() && out.back().type == e.type &&
            (gap_ok(prev_hi, lo) || gap_ok(prev_hhi, hlo))) {
            Event& prev = out.back();
            for (const std::string& nd : e.nodes) prev.nodes.push_back(nd);
            prev.end_node = e.end_node;
            if (!e.seq.empty()) prev.seq += e.seq;
            prev.size_bp += e.size_bp;
            if (hi >= 0) prev_hi = std::max(prev_hi, hi);
            if (hhi >= 0) prev_hhi = std::max(prev_hhi, hhi);
        } else {
            out.push_back(std::move(e));
            prev_hi = hi;
            prev_hhi = hhi;
        }
    }
    events = std::move(out);
}

double weighted_jaccard(
    const Graph& graph,
    const std::vector<std::string>& a,
    const std::vector<std::string>& b) {

    std::unordered_set<std::string> sa(a.begin(), a.end());
    std::unordered_set<std::string> sb(b.begin(), b.end());
    std::size_t inter = 0;
    std::size_t uni = 0;
    std::unordered_set<std::string> seen;
    for (const std::string& n : sa) {
        seen.insert(n);
        const std::size_t l = node_len(graph, n);
        uni += l;
        if (sb.count(n)) inter += l;
    }
    for (const std::string& n : sb) {
        if (seen.count(n)) continue;
        uni += node_len(graph, n);
    }
    return uni == 0 ? 0.0 : static_cast<double>(inter) / static_cast<double>(uni);
}

// Banded-alignment identity between two event sequences, gated by len_ratio (min shorter/longer to
// even compare - decoupled from min_id so same-motif STR alleles of different size can still merge).
double seq_identity(const std::string& a, const std::string& b, double min_id, double len_ratio) {
    if (a.empty() || b.empty()) return 0.0;
    const std::size_t lo = std::min(a.size(), b.size());
    const std::size_t hi = std::max(a.size(), b.size());
    if (static_cast<double>(lo) < len_ratio * static_cast<double>(hi)) return 0.0;
    const std::string& shorter = a.size() <= b.size() ? a : b;
    const std::string& longer = a.size() <= b.size() ? b : a;
    const std::size_t band = std::max<std::size_t>(
        8, static_cast<std::size_t>((1.0 - min_id) * static_cast<double>(shorter.size())) + 8);
    const FitAlignResult fa = fit_align(shorter, longer, band);
    return fa.ok ? fa.identity : 0.0;
}

// A merged event across haplotypes: a representative event + its carriers.
struct MergedRecord {
    Event seed;                                             // representative (largest member)
    std::vector<std::string> carriers;                      // haplotype path names with GT=1
    std::unordered_set<std::size_t> member_alleles;         // allele indices already merged in
    std::set<std::string> member_nodes;                     // union of every merged member event's
                                                            // nodes (describe handoff; ordered set)
    std::unordered_map<std::string, std::size_t> sample_cn; // DUP per-sample copy number
    std::size_t min_size_bp = 0;                            // smallest merged member size (SVLEN_RANGE)
    std::size_t max_size_bp = 0;                            // largest merged member size
    double merge_max_jaccard = -1.0;                        // strongest node-set Jaccard that merged a member (-1 = none)
    double merge_max_seqid = -1.0;                          // strongest sequence identity that merged a member (-1 = none)
};

std::string upper_base(char c) {
    const char u = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    if (u == 'A' || u == 'C' || u == 'G' || u == 'T' || u == 'N') return std::string(1, u);
    return "N";
}

} // namespace

void call_variants(
    const Graph& graph,
    const VariantCallOptions& options,
    VariantCallSummary* summary_out) {

    if (options.reference_path.empty()) {
        throw std::runtime_error("call requires --reference-path");
    }
    if (options.out_prefix.empty()) {
        throw std::runtime_error("call requires -o/--out-prefix");
    }

    // Locate the reference path record. Accept either the full path name or a
    // case-insensitive substring (e.g. "grch38" -> "grch38#1#chr6:..."). An exact name
    // always wins; otherwise the substring must match exactly one path, else we error
    // (not found, or ambiguous with the candidate list).
    const PathRecord* ref_path = nullptr;
    for (const PathRecord& p : graph.paths) {
        if (p.name == options.reference_path) { ref_path = &p; break; }
    }
    if (ref_path == nullptr) {
        auto lower = [](std::string s) {
            std::transform(s.begin(), s.end(), s.begin(),
                           [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
            return s;
        };
        const std::string needle = lower(options.reference_path);
        std::vector<const PathRecord*> hits;
        for (const PathRecord& p : graph.paths) {
            if (lower(p.name).find(needle) != std::string::npos) hits.push_back(&p);
        }
        if (hits.size() == 1) {
            ref_path = hits.front();
        } else if (hits.empty()) {
            throw std::runtime_error("Reference path not found in GFA: " + options.reference_path);
        } else {
            std::string msg = "Reference path '" + options.reference_path + "' is ambiguous; matches " +
                              std::to_string(hits.size()) + " paths:";
            for (const PathRecord* p : hits) msg += "\n  " + p->name;
            throw std::runtime_error(msg);
        }
    }
    const std::string& ref_name = ref_path->name;

    const std::vector<Bubble> bubbles = read_bubbles_csv(options.bubbles_csv_in);
    std::unordered_set<std::size_t> bubble_filter(options.bubble_ids.begin(), options.bubble_ids.end());

    const ParsedReferencePath ref_meta = parse_reference_path_label(ref_name);
    const std::vector<std::size_t> ref_prefix = path_prefix_bp(*ref_path, graph.nodes);

    // ---- Optional GTF gene annotation. The GTF is in reference coordinates; we read the
    // reference path's chrom + absolute start (PanSN) and project genes onto reference nodes.
    // node_genes maps a reference node id -> indices into `genes`. Built once; const in the
    // parallel loop. Annotation is skipped (genes stays empty) when --gtf is unset or the
    // reference name is not PanSN.
    std::vector<GeneFeature> genes;
    std::unordered_map<std::string, std::vector<int>> node_genes;
    if (!options.gtf_path.empty()) {
        if (!is_pansn(ref_name)) {
            std::cerr << "warning: --gtf given but reference path '" << ref_name
                      << "' is not PanSN (sample#hap#contig:start-end); skipping gene annotation\n";
        } else {
            const std::size_t lo = ref_meta.region_start_1based;
            const std::size_t hi = lo + (ref_prefix.empty() ? 0 : ref_prefix.back()) - 1;
            genes = parse_gtf(options.gtf_path, ref_meta.chrom, lo, hi);
            node_genes = project_genes_to_nodes(graph, *ref_path, ref_meta, genes);
        }
    }
    // Sorted, unique gene names overlapping a set of node ids (for the INFO GENES field).
    auto genes_for_nodes = [&](const std::vector<std::string>& ids) {
        std::vector<std::string> names;
        std::unordered_set<int> seen;
        for (const std::string& id : ids) {
            const auto it = node_genes.find(id);
            if (it == node_genes.end()) continue;
            for (int gi : it->second) if (seen.insert(gi).second) names.push_back(genes[gi].gene_name);
        }
        std::sort(names.begin(), names.end());
        return names;
    };

    // Fixed sample list: every haplotype path (the reference included).
    std::vector<std::string> sample_names;
    sample_names.reserve(graph.paths.size());
    for (const PathRecord& p : graph.paths) sample_names.push_back(p.name);

    // Per-path bubble index, reused across bubbles.
    std::vector<BubblePathIndex> path_indexes;
    path_indexes.reserve(graph.paths.size());
    for (const PathRecord& p : graph.paths) path_indexes.push_back(build_bubble_path_index(p));

    const std::unordered_set<std::string> selfloops = self_loop_nodes(graph);

    cli::ensure_parent_dir_for_file(options.out_prefix + ".region.vcf");

    auto write_vcf_header = [&](std::ostream& out) {
        out << "##fileformat=VCFv4.2\n";
        out << "##source=panvar call\n";
        out << "##reference=" << ref_name << "\n";
        out << "##contig=<ID=" << ref_meta.chrom << ">\n";
        out << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant\">\n";
        out << "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Structural variant type\">\n";
        out << "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length difference ALT-REF\">\n";
        out << "##INFO=<ID=SVLEN_RANGE,Number=2,Type=Integer,Description=\"Min,max event size among merged members (when they differ)\">\n";
        out << "##INFO=<ID=BUBBLE_ID,Number=1,Type=Integer,Description=\"panvar bubble identifier\">\n";
        out << "##INFO=<ID=START_NODE,Number=1,Type=String,Description=\"First graph node of the event\">\n";
        out << "##INFO=<ID=END_NODE,Number=1,Type=String,Description=\"Last graph node of the event\">\n";
        out << "##INFO=<ID=EVENT_NODES,Number=.,Type=String,Description=\"Variant node set\">\n";
        out << "##INFO=<ID=INS_SUBTYPE,Number=1,Type=String,Description=\"INS subtype: NOVEL or DUP (minimap2 refined)\">\n";
        out << "##INFO=<ID=REF_CN,Number=1,Type=Integer,Description=\"Reference copy number of the repeat unit (DUP)\">\n";
        out << "##INFO=<ID=RU_LEN,Number=1,Type=Integer,Description=\"Repeat-unit length in bp, one copy (DUP)\">\n";
        if (!genes.empty())
            out << "##INFO=<ID=GENES,Number=.,Type=String,Description=\"Gene(s) overlapping the variant (from --gtf)\">\n";
        out << "##INFO=<ID=NMERGED,Number=1,Type=Integer,Description=\"Haplotype carriers merged into this record\">\n";
        out << "##INFO=<ID=MERGE_JACCARD,Number=1,Type=Float,Description=\"Strongest node-set Jaccard that merged a member into this record (cross-haplotype merge evidence)\">\n";
        out << "##INFO=<ID=MERGE_SEQID,Number=1,Type=Float,Description=\"Strongest sequence identity that merged a member into this record, when the Jaccard gate did not decide it\">\n";
        out << "##INFO=<ID=MERGE_SIZE_RATIO,Number=1,Type=Float,Description=\"Smallest/largest member size among merged members (min,max also in SVLEN_RANGE)\">\n";
        out << "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Allele number = haplotypes traversing the bubble\">\n";
        out << "##INFO=<ID=AC,Number=1,Type=Integer,Description=\"Allele count = carrier haplotypes\">\n";
        out << "##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele frequency = AC/AN (over traversing haplotypes)\">\n";
        out << "##INFO=<ID=NALLELES,Number=1,Type=Integer,Description=\"Number of alleles (REF+ALTs) at a multiallelic locus\">\n";
        out << "##INFO=<ID=EVENTID,Number=1,Type=String,Description=\"Shared id linking a co-located DEL+INS substitution\">\n";
        out << "##INFO=<ID=INSSEQ,Number=1,Type=String,Description=\"Inserted sequence\">\n";
        out << "##INFO=<ID=DELSEQ,Number=1,Type=String,Description=\"Deleted reference sequence\">\n";
        out << "##INFO=<ID=INVSEQ,Number=1,Type=String,Description=\"Inverted reference sequence\">\n";
        out << "##ALT=<ID=DEL,Description=\"Deletion\">\n";
        out << "##ALT=<ID=INS,Description=\"Insertion\">\n";
        out << "##ALT=<ID=INV,Description=\"Inversion\">\n";
        out << "##ALT=<ID=DUP,Description=\"Copy-number / tandem duplication locus\">\n";
        out << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"1=carrier, 0=reference-like, .=bubble not traversed\">\n";
        out << "##FORMAT=<ID=CN,Number=1,Type=Integer,Description=\"Copy number of the repeat unit (DUP records)\">\n";
        out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
        for (const std::string& s : sample_names) out << '\t' << s;
        out << '\n';
    };

    // Records are buffered and coordinate-sorted before writing, so the region VCF (and
    // per-bubble VCFs) are sorted and bgzip/tabix/bcftools-indexable.
    struct OutRecord {
        std::size_t pos = 0;
        std::size_t end = 0;
        std::size_t bubble_id = 0;
        std::string id;
        std::string line;                       // full VCF row (with trailing newline)
        std::vector<std::string> prov_lines;    // variant_paths.tsv rows for this record
    };
    std::vector<OutRecord> out_records;
    std::vector<std::string> variant_nodes_rows;     // <prefix>.variant_nodes.tsv (describe handoff)
    std::vector<std::string> dup_gene_cn_rows;       // <prefix>.dup_gene_cn.tsv body (per-gene DUP CN)

    VariantCallSummary summary;

    // Per-bubble outputs, computed in parallel and merged in bubble order so the result is identical
    // regardless of thread count: the region VCF is coordinate-sorted afterwards (total order on
    // unique ids), and variant_nodes inherits deterministic bubble order from the merge.
    // A DUP whose per-gene copy number will be resolved by realignment (post-pass).
    struct DupGeneTarget {
        std::size_t bubble_id = 0;
        std::string variant_id;
        std::vector<int> gene_idx;   // indices into `genes` of the genes overlapping the bubble
        std::unordered_map<std::string, std::size_t> sample_cn;  // module CN per haplotype (FORMAT:CN)
    };
    std::vector<DupGeneTarget> dup_targets;

    struct BubbleOut {
        VariantCallSummary sum;
        std::vector<OutRecord> records;
        std::vector<std::string> variant_nodes;
        std::vector<DupGeneTarget> dup_targets;   // DUPs needing per-gene CN (when --gtf is active)
    };

    // Genomic 1-based start of a bubble node from the reference path (first occurrence).
    auto build_ref_node_pos = [&](const Bubble& bubble) {
        std::unordered_map<std::string, std::size_t> pos; // node -> 1-based genomic start
        std::unordered_set<std::string> bnodes(bubble.inside.begin(), bubble.inside.end());
        bnodes.insert(bubble.source);
        bnodes.insert(bubble.sink);
        for (std::size_t k = 0; k < ref_path->steps.size(); ++k) {
            const std::string& id = ref_path->steps[k].node_id;
            if (bnodes.count(id) && !pos.count(id)) {
                pos[id] = ref_meta.region_start_1based + ref_prefix[k];
            }
        }
        return pos;
    };

    // Progress over bubbles (stderr, TTY-gated; empty label = suppressed under --quiet).
    cli::ProgressBar call_progress(options.quiet ? "" : "Calling variants", bubbles.size());
    std::vector<BubbleOut> bouts(bubbles.size());
    // Bubbles are processed SERIALLY; parallelism lives INSIDE each bubble, over its alleles (the
    // per-allele walk-diff DP is the hot path). Copy-number/SV loci are typically one dominant folded
    // bubble that bubble-level parallelism could not split, so this is where the cores are needed.
    // Each bubble still writes only to its own bouts[bubble_idx], so results are order-deterministic.
    auto process_bubble = [&](std::size_t bubble_idx) {
        const Bubble& bubble = bubbles[bubble_idx];
        // Shadow the shared sinks with this bubble's local buffers: the (unchanged) loop body below
        // writes only through these names, so each bubble accumulates independently. id_counts is
        // bubble-scoped (variant ids embed the bubble id), so a local map reproduces the ids exactly.
        BubbleOut& bout = bouts[bubble_idx];
        VariantCallSummary& summary = bout.sum;
        std::vector<OutRecord>& out_records = bout.records;
        std::vector<std::string>& variant_nodes_rows = bout.variant_nodes;
        std::vector<DupGeneTarget>& dup_targets = bout.dup_targets;
        std::unordered_map<std::string, int> id_counts;
        if (!bubble_filter.empty() && bubble_filter.find(bubble.id) == bubble_filter.end()) return;
        ++summary.bubbles_seen;

        // Reference walk through this bubble. path_indexes is parallel to graph.paths.
        std::size_t ref_idx = graph.paths.size();
        for (std::size_t pi = 0; pi < graph.paths.size(); ++pi) {
            if (graph.paths[pi].name == ref_name) { ref_idx = pi; break; }
        }
        std::optional<std::vector<PathStep>> ref_opt =
            (ref_idx < graph.paths.size())
                ? bubble_steps(graph.paths[ref_idx], path_indexes[ref_idx], bubble)
                : std::optional<std::vector<PathStep>>{};
        if (!ref_opt.has_value() || ref_opt->empty()) {
            return; // reference does not traverse this bubble; cannot type events
        }
        ++summary.bubbles_with_reference;

        const std::vector<PathStep> ref_steps = std::move(*ref_opt);

        // Distinct alleles: group paths by canonical-walk signature.
        struct Allele { std::vector<PathStep> steps; std::vector<std::string> members; };
        std::unordered_map<std::string, std::size_t> sig_to_allele;
        std::vector<Allele> alleles;
        std::unordered_set<std::string> traverses; // path names that cross the bubble
        const std::string ref_sig = build_walk_signature(ref_steps);
        for (std::size_t pi = 0; pi < graph.paths.size(); ++pi) {
            const auto steps_opt = bubble_steps(graph.paths[pi], path_indexes[pi], bubble);
            if (!steps_opt.has_value() || steps_opt->empty()) continue;
            const std::vector<PathStep>& steps = *steps_opt;
            traverses.insert(graph.paths[pi].name);
            const std::string sig = build_walk_signature(steps);
            auto it = sig_to_allele.find(sig);
            if (it == sig_to_allele.end()) {
                sig_to_allele.emplace(sig, alleles.size());
                Allele a; a.steps = steps; a.members.push_back(graph.paths[pi].name);
                alleles.push_back(std::move(a));
            } else {
                alleles[it->second].members.push_back(graph.paths[pi].name);
            }
        }

        // CN nodes: self-loop nodes present in this bubble (a panphorte REP / genuine
        // tandem unit). They are handled as count-based DUP/CN events and excluded from
        // the DEL/INS/INV alignment. Ordinary nodes that merely recur (e.g. shared module
        // content) are NOT copy-number loci and stay in the alignment.
        std::unordered_map<std::string, std::size_t> ref_count;
        for (const PathStep& s : ref_steps) ++ref_count[s.node_id];
        std::unordered_set<std::string> cn_nodes;
        for (const auto& [id, c] : ref_count) { (void)c; if (selfloops.count(id)) cn_nodes.insert(id); }
        for (const Allele& a : alleles) {
            for (const PathStep& s : a.steps) if (selfloops.count(s.node_id)) cn_nodes.insert(s.node_id);
        }

        const std::vector<Tok> Rtok = collapse_walk(ref_steps, cn_nodes);

        std::vector<MergedRecord> merged;
        const std::unordered_map<std::string, std::size_t> ref_node_pos = build_ref_node_pos(bubble);

        const std::size_t rescue_floor =
            options.rescue_min_bp != 0 ? options.rescue_min_bp : std::max<std::size_t>(1, options.min_sv_bp / 2);

        auto ev_ref_pos = [&](const Event& e) -> long long {
            const auto it = ref_node_pos.find(e.anchor_node);
            if (it == ref_node_pos.end()) return -1;
            const std::size_t glen = node_len(graph, e.anchor_node);
            return static_cast<long long>(it->second + (e.anchor_after && glen > 0 ? glen - 1 : 0));
        };
        // Two non-DUP events are the same site: same type, anchors within the window,
        // and either node sets overlap (Jaccard) OR sequences are similar.
        // Returns whether a and b are the same event. When out_jac/out_seq are given, they
        // receive the node-set Jaccard and (only if the Jaccard gate did not already decide it)
        // the sequence identity actually computed, so a merged record can report the evidence.
        auto events_match = [&](const Event& a, const Event& b,
                                double* out_jac = nullptr, double* out_seq = nullptr) {
            if (out_jac) *out_jac = -1.0;
            if (out_seq) *out_seq = -1.0;
            if (a.type != b.type) return false;
            if (a.ref_pos == 0 || b.ref_pos == 0) {
                // fall back to node/seq only when no coordinate
            } else {
                // Position window scales with event size: a large insertion/deletion can
                // attach at a fuzzy breakpoint that varies by kilobases across haplotypes,
                // yet still be the same event. The sequence/Jaccard gate below keeps this
                // from over-merging genuinely distinct events at nearby positions.
                const long long window = static_cast<long long>(options.merge_distance_bp) +
                    static_cast<long long>(std::min(a.size_bp, b.size_bp));
                const long long d = static_cast<long long>(a.ref_pos) - static_cast<long long>(b.ref_pos);
                if (d > window || d < -window) return false;
            }
            const double jac = weighted_jaccard(graph, a.nodes, b.nodes);
            if (out_jac) *out_jac = jac;
            if (jac >= options.merge_jaccard) return true;
            const double len_ratio =
                options.merge_size_ratio > 0.0 ? options.merge_size_ratio : options.merge_seq_identity;
            const double sid = seq_identity(a.seq, b.seq, options.merge_seq_identity, len_ratio);
            if (out_seq) *out_seq = sid;
            return sid >= options.merge_seq_identity;
        };

        // ---- bp-coverage CN (--cn-from-coverage): total module/cluster copy number on a
        // pggb-collapsed paralog cluster, where the reference itself traverses the module >=2x.
        // PGGB collapses IDENTICAL copies (e.g. the C4 long-long / short-short modules, or the
        // CYP2D6/2D7/2D8P paralogs) onto shared nodes, so a 2-copy haplotype re-traverses those
        // nodes twice -- the copy number is carried as NODE MULTIPLICITY, not as a tandem block in
        // the spelled sequence. Copy number = (full-walk bp through the bubble) / (unit bp), with
        // unit = ref_full_bp / ref_fold (ref's peak node multiplicity over its own full walk).
        //
        // The full walk is the WIDEST source..sink span (all repeats included). This is the crux:
        // canonical_bubble_path_steps (used everywhere else, incl. the alleles above and the
        // self-loop counts below) takes the MINIMAL span covering the distinct inside nodes and so
        // collapses the repeated copies onto one traversal -- which would flatten every haplotype to
        // the same bp. Counting over the full span preserves the multiplicity == copy number. When
        // it fires it is the authority for the bubble, so the self-loop / walk-diff paths are skipped.
        // Validated: C4 total CN 131/131 exact; CYP2D6 ~92% after the paralog baseline (2D7+2D8P)
        // offset the comparator subtracts.
        //
        // Coverage is the route for PGGB-folded paralogs. It is suppressed only when the bubble carries a
        // genuine panphorte REP self-loop (a repeat unit >= min_sv_bp, e.g. LPA's 5.5 kb KIV-2 node):
        // there the self-loop's integer traversal count is exact while coverage's bp/unit ratio only
        // estimates it, so the self-loop must win. An *incidental* tiny PGGB self-loop (e.g. C4's 22 bp
        // node 1962, well below min_sv_bp) is NOT a repeat unit, so coverage still fires there — keeping
        // C4/CYP2D6 on coverage while letting LPA use the exact self-loop. This makes the routes disjoint
        // by topology (self-loop > coverage > peak) so passing both CN flags is safe on either graph.
        bool has_rep_selfloop = false;
        for (const std::string& cn : cn_nodes)
            if (node_len(graph, cn) >= options.min_sv_bp) { has_rep_selfloop = true; break; }
        bool coverage_fired = false;
        double cn_unit_bp = 0.0;   // one-copy bp of the coverage module (set when coverage fires)
        if (options.cn && !has_rep_selfloop) {
            const std::unordered_set<std::string> inside_set(bubble.inside.begin(), bubble.inside.end());
            // Widest oriented source..sink span of a path over this bubble (all repeats included):
            // forward = first source .. last sink; reversed crossing = first sink .. last source.
            auto span_of = [&](std::size_t pi) -> std::pair<std::size_t, std::size_t> {
                const auto& idx = path_indexes[pi].positions;
                const auto sit = idx.find(bubble.source);
                const auto kit = idx.find(bubble.sink);
                if (sit == idx.end() || kit == idx.end()) return {1, 0};  // empty (lo>hi)
                const std::size_t s0 = sit->second.front(), s1 = sit->second.back();
                const std::size_t k0 = kit->second.front(), k1 = kit->second.back();
                const std::size_t lo = (s0 <= k1) ? s0 : k0;
                const std::size_t hi = (s0 <= k1) ? k1 : s1;
                return {lo, hi};
            };
            // The FOLDED set = inside nodes the REFERENCE visits >=2x, i.e. the actual repeated unit.
            // Copy number is measured over these nodes ONLY, so a haplotype's deletion/insertion of
            // UNIQUE (single-visit) content -- interstitial sequence, or paralog-specific nodes that do
            // not recur -- is not misread as a copy loss/gain (it stays a DEL/INS). This is
            // ratio-preserving on genuine folded modules (every folded node recurs with the module) yet
            // immune to unique-content edits inside the bubble.
            std::unordered_set<std::string> folded_set;
            std::size_t ref_fold = 0;
            std::size_t ref_bp = 0;
            if (ref_idx < graph.paths.size()) {
                const auto ref_span = span_of(ref_idx);
                const std::vector<PathStep>& steps = graph.paths[ref_idx].steps;
                std::unordered_map<std::string, std::size_t> cnt;
                for (std::size_t i = ref_span.first; i <= ref_span.second && i < steps.size(); ++i)
                    if (inside_set.count(steps[i].node_id)) ++cnt[steps[i].node_id];
                for (const auto& kv : cnt) {
                    if (kv.second >= 2) folded_set.insert(kv.first);
                    ref_fold = std::max(ref_fold, kv.second);
                }
                for (std::size_t i = ref_span.first; i <= ref_span.second && i < steps.size(); ++i)
                    if (folded_set.count(steps[i].node_id)) ref_bp += node_len(graph, steps[i].node_id);
            }
            auto full_walk_bp = [&](std::size_t pi) -> std::size_t {
                const auto sp = span_of(pi);
                const std::vector<PathStep>& steps = graph.paths[pi].steps;
                std::size_t bp = 0;
                for (std::size_t i = sp.first; i <= sp.second && i < steps.size(); ++i)
                    if (folded_set.count(steps[i].node_id)) bp += node_len(graph, steps[i].node_id);
                return bp;
            };
            // Require the one-copy unit to reach min_sv_bp: a sub-threshold fold (e.g. an 8 bp repeat the
            // reference happens to traverse twice) is not an SV-scale copy-number module and would only add
            // tiny noise DUPs, so coverage does not fire there.
            if (ref_fold >= 2 && ref_bp > 0 &&
                static_cast<double>(ref_bp) / static_cast<double>(ref_fold) >= static_cast<double>(options.min_sv_bp)) {
                const double unit = static_cast<double>(ref_bp) / static_cast<double>(ref_fold);
                const std::size_t ref_copies = ref_fold;
                MergedRecord mr;
                mr.seed.type = EvType::Dup;
                mr.seed.nodes.push_back(bubble.source);
                mr.seed.start_node = bubble.source; mr.seed.end_node = bubble.sink;
                mr.seed.ref_cn = ref_copies; mr.seed.alt_cn = ref_copies;
                mr.seed.anchor_node = bubble.source;
                mr.seed.size_bp = static_cast<std::size_t>(unit);
                mr.seed.ru_len = static_cast<std::size_t>(unit);
                mr.seed.cn_peak = true;
                for (std::size_t pi = 0; pi < graph.paths.size(); ++pi) {
                    if (!traverses.count(graph.paths[pi].name)) continue;
                    const std::size_t hbp = full_walk_bp(pi);
                    if (hbp == 0) continue;
                    const std::size_t copies =
                        static_cast<std::size_t>(std::llround(static_cast<double>(hbp) / unit));
                    // CN is reported for every traversing haplotype (absolute module count); GT marks a
                    // CARRIER only when the count differs from the reference's (a gain or a loss), so AC/AF
                    // stay meaningful instead of flagging every haplotype.
                    mr.sample_cn[graph.paths[pi].name] = copies;
                    if (copies != ref_copies) mr.carriers.push_back(graph.paths[pi].name);
                    if (copies > mr.seed.alt_cn) mr.seed.alt_cn = copies;
                }
                if (!mr.carriers.empty()) {  // a CN-invariant module is not a variant record
                    // describe handoff: a coverage DUP's copy-number signal lives in the folded module's
                    // INSIDE nodes (their per-walk multiplicity), not in the bubble source. EVENT_NODES /
                    // POS stay the compact source anchor, but variant_nodes.tsv must carry the inside
                    // nodes so `describe --variant-nodes` retains the dosage features (otherwise it masks
                    // to the source node, which every haplotype traverses once → no variance → dropped).
                    mr.member_nodes.insert(bubble.inside.begin(), bubble.inside.end());
                    coverage_fired = true;
                    cn_unit_bp = unit;   // one-copy size, used to drop CN-loss DELs below
                    merged.push_back(std::move(mr));
                }
            }
        }

        // ---- DUP/CN events: count self-loop traversals per allele vs reference. Merge
        // on shared REP node; per-sample CN. Independent of the walk-diff alignment.
        // Skipped when bp-coverage fired (it is the authority for that bubble's copy number).
        for (std::size_t ai = 0; !coverage_fired && ai < alleles.size(); ++ai) {
            std::unordered_map<std::string, std::size_t> alt_count;
            for (const PathStep& s : alleles[ai].steps) ++alt_count[s.node_id];
            for (const std::string& cn : cn_nodes) {
                // Only a genuine REP repeat unit (self-loop node >= min_sv_bp) anchors a self-loop DUP. An
                // incidental tiny self-loop (e.g. C4's 22 bp node the reference may not even traverse) would
                // otherwise emit a spurious DUP with REF_CN=0; with the guard it is simply skipped, so a
                // mis-chosen CN flag misses rather than corrupts.
                if (node_len(graph, cn) < options.min_sv_bp) continue;
                const std::size_t rc = ref_count.count(cn) ? ref_count.at(cn) : 0;
                const std::size_t ac = alt_count.count(cn) ? alt_count.at(cn) : 0;
                if (rc == ac) continue;
                const std::size_t unit = node_len(graph, cn);
                const std::size_t delta = rc > ac ? rc - ac : ac - rc;
                if (unit * delta < options.min_sv_bp) continue;
                MergedRecord* grp = nullptr;
                for (MergedRecord& mr : merged) {
                    if (mr.seed.type == EvType::Dup && mr.seed.nodes.front() == cn) { grp = &mr; break; }
                }
                if (grp == nullptr) {
                    MergedRecord mr;
                    mr.seed.type = EvType::Dup;
                    mr.seed.nodes.push_back(cn);
                    mr.seed.start_node = cn; mr.seed.end_node = cn;
                    mr.seed.ref_cn = rc; mr.seed.alt_cn = ac;
                    mr.seed.anchor_node = bubble.source;
                    mr.seed.size_bp = unit * delta;
                    mr.seed.ru_len = unit;
                    merged.push_back(std::move(mr));
                    grp = &merged.back();
                }
                for (const std::string& m : alleles[ai].members) {
                    grp->carriers.push_back(m);
                    grp->sample_cn[m] = ac;
                }
            }
        }

        // ---- Peak-multiplicity DUP (--cn-from-multiplicity): for a folded bubble that
        // panphorte could not collapse (no self-loop CN node, e.g. the GSTM1 segdup
        // cluster), copy number is the per-haplotype PEAK node traversal multiplicity. A
        // haplotype that folds an extra paralog copy onto a shared node traverses that node
        // once more than the reference does at its own peak. Firing on the PEAK (not on any
        // node exceeding the reference) is what rejects cluster background: scattered
        // per-node excesses reflect paralog presence/absence, the global peak reflects gene
        // dosage. Validated on GSTM1: only the 2 true dup haplotypes (peak 4 > ref peak 3)
        // fire, 0 false positives. Skipped entirely when self-loop CN nodes exist, so
        // panphorte-normalized C4/LPA keep their exact fit_align counts.
        // Gated on the absence of a REP self-loop >= min_sv_bp (not on cn_nodes being empty): an incidental
        // tiny self-loop (C4's 22 bp node) must NOT block the per-node/peak estimate. The self-loop DUP route
        // above already guards each loop on size, so the two routes stay disjoint by topology.
        if (options.cn && !has_rep_selfloop && !coverage_fired) {
            std::size_t ref_peak = 0;
            for (const std::string& id : bubble.inside) {
                const auto it = ref_count.find(id);
                if (it != ref_count.end() && it->second > ref_peak) ref_peak = it->second;
            }
            MergedRecord* peak_grp = nullptr;
            for (std::size_t ai = 0; ai < alleles.size(); ++ai) {
                std::unordered_map<std::string, std::size_t> alt_count;
                for (const PathStep& s : alleles[ai].steps) ++alt_count[s.node_id];
                // ac_peak / peak_node decide IF this allele duplicates (gene dosage). The
                // peak node is often tiny, so its length is not the event size; size is the
                // extra sequence the allele traverses vs the reference, summed over all
                // folded nodes (excess_bp) — the actual duplicated content (~18.5 kb on the
                // GSTM1 carriers vs <7 kb for non-firing haplotypes).
                std::string peak_node;
                std::size_t ac_peak = 0;
                std::size_t excess_bp = 0;
                for (const std::string& id : bubble.inside) {
                    const auto it = alt_count.find(id);
                    const std::size_t c = it == alt_count.end() ? 0 : it->second;
                    if (c > ac_peak) { ac_peak = c; peak_node = id; }
                    const std::size_t rc_id = ref_count.count(id) ? ref_count.at(id) : 0;
                    if (c > rc_id) excess_bp += node_len(graph, id) * (c - rc_id);
                }
                if (ac_peak < 2 || ac_peak <= ref_peak) continue;
                if (excess_bp < options.min_sv_bp) continue;
                if (peak_grp == nullptr) {
                    MergedRecord mr;
                    mr.seed.type = EvType::Dup;
                    mr.seed.nodes.push_back(peak_node);
                    mr.seed.start_node = peak_node; mr.seed.end_node = peak_node;
                    mr.seed.ref_cn = ref_peak; mr.seed.alt_cn = ac_peak;
                    mr.seed.anchor_node = bubble.source;
                    mr.seed.size_bp = excess_bp;
                    mr.seed.ru_len = excess_bp / (ac_peak - ref_peak);  // per-copy duplicated content
                    mr.seed.cn_peak = true;
                    // describe handoff: like the coverage route, a peak DUP's copy-number signal lives in
                    // the folded module's INSIDE nodes (their per-walk multiplicity), not in the single peak
                    // node. EVENT_NODES / POS stay the compact peak-node anchor, but variant_nodes.tsv must
                    // carry the inside nodes so `describe --variant-nodes` keeps the dosage features (else it
                    // masks to one node every haplotype traverses the same way → no variance → dropped).
                    mr.member_nodes.insert(bubble.inside.begin(), bubble.inside.end());
                    merged.push_back(std::move(mr));
                    peak_grp = &merged.back();
                } else if (ac_peak > peak_grp->seed.alt_cn) {
                    peak_grp->seed.alt_cn = ac_peak;
                    peak_grp->seed.size_bp = excess_bp;
                    peak_grp->seed.ru_len = excess_bp / (ac_peak - ref_peak);
                }
                for (const std::string& m : alleles[ai].members) {
                    peak_grp->carriers.push_back(m);
                    peak_grp->sample_cn[m] = ac_peak;
                }
            }
        }

        // ---- DEL/INS/INV events: derive per allele, keep ALL (down to the rescue floor)
        // for the re-scan, then merge events >= floor by position + sequence/node match.
        // These are emitted even when bp-coverage fired: the coverage DUP is a copy-number DOSAGE
        // for the folded module, but it does not represent sequence-resolved insertions/inversions/
        // deletions inside the bubble (e.g. a rare gene-unit INS), so suppressing them would silently
        // drop real calls. The CN-DUP routes (self-loop / peak above) stay gated on !coverage_fired so
        // copy number is reported once; the walk-diff events coexist with the coverage DUP.
        // The per-allele walk-diff (diff_walks' O(m*n) segment DP) dominates the whole caller on
        // large folded bubbles (e.g. the GSTM CNV: ~hundreds of distinct alleles, large inter-anchor
        // segments). It is pure and writes only to its own allele_events[ai], so run it across cores;
        // the cross-haplotype merge below consumes the results in the same (allele) order, so output
        // is independent of thread count.
        std::vector<std::vector<Event>> allele_events(alleles.size());
        run_parallel(alleles.size(), options.threads, [&](std::size_t ai) {
            if (build_walk_signature(alleles[ai].steps) == ref_sig) return;
            const std::vector<Tok> Htok = collapse_walk(alleles[ai].steps, cn_nodes);
            std::vector<Event> events = diff_walks(graph, Rtok, Htok, bubble.source);
            // Per-haplotype node -> first bp offset along this allele's walk, so coalescing
            // can also measure same-type-event gaps in the sample's own sequence space.
            std::unordered_map<std::string, std::size_t> hap_node_pos;
            {
                std::size_t off = 0;
                for (const PathStep& s : alleles[ai].steps) {
                    hap_node_pos.emplace(s.node_id, off);
                    off += node_len(graph, s.node_id);
                }
            }
            coalesce_events(graph, events, ref_node_pos, hap_node_pos, options.merge_distance_bp);
            for (Event& e : events) {
                const long long p = ev_ref_pos(e);
                e.ref_pos = p < 0 ? 0 : static_cast<std::size_t>(p);
            }
            allele_events[ai] = std::move(events);
        });

        // ---- Cross-haplotype merge by TRANSITIVE single-linkage clustering. Greedy
        // first-fit (each event joins the first matching seed) is not transitive: if A~B
        // and B~C but A!~C, then C spawns its own record and the same biological variant
        // fragments across haplotypes. Instead build a graph whose edges are events_match
        // and take connected components, so A-B-C collapse to one record. On clean sites
        // where greedy already merged everything into one component the result is identical;
        // it can only ever merge MORE, never split or drop a carrier.
        {
            struct Cand { std::size_t ai; const Event* e; };
            std::vector<Cand> cands;
            for (std::size_t ai = 0; ai < alleles.size(); ++ai)
                for (const Event& e : allele_events[ai])
                    if (e.size_bp >= rescue_floor) cands.push_back({ai, &e});

            DisjointSet dsu(cands.size());
            // Position-sorted windowed sweep keeps edge-building near-linear: events_match
            // can only fire within |dpos| <= merge_distance_bp + min(size) <= merge_distance_bp
            // + size_i, so once an earlier (smaller-pos) event is past that bound we stop.
            // Zero-coordinate events (rare; anchor missing from the reference) match on
            // node/sequence only, so they are never position-pruned.
            std::vector<std::size_t> ord(cands.size());
            std::iota(ord.begin(), ord.end(), std::size_t{0});
            std::sort(ord.begin(), ord.end(), [&](std::size_t a, std::size_t b) {
                return cands[a].e->ref_pos < cands[b].e->ref_pos;
            });
            // Per-candidate strongest merge evidence (over the edges that joined it), carried
            // into the record so it can report why members were merged.
            std::vector<double> best_jac(cands.size(), -1.0), best_seq(cands.size(), -1.0);
            for (std::size_t a = 0; a < ord.size(); ++a) {
                const std::size_t i = ord[a];
                const Event& ei = *cands[i].e;
                for (std::size_t b = a; b-- > 0;) {
                    const std::size_t j = ord[b];
                    const Event& ej = *cands[j].e;
                    if (ei.ref_pos != 0 && ej.ref_pos != 0 &&
                        static_cast<long long>(ei.ref_pos) - static_cast<long long>(ej.ref_pos) >
                            static_cast<long long>(options.merge_distance_bp) +
                                static_cast<long long>(ei.size_bp))
                        break;
                    if (dsu.find(i) == dsu.find(j)) continue;
                    double jac = -1.0, sid = -1.0;
                    if (events_match(ei, ej, &jac, &sid)) {
                        dsu.unite(i, j);
                        best_jac[i] = std::max(best_jac[i], jac); best_jac[j] = std::max(best_jac[j], jac);
                        best_seq[i] = std::max(best_seq[i], sid); best_seq[j] = std::max(best_seq[j], sid);
                    }
                }
            }

            // Collapse components into records (deterministic order: first cand per root).
            std::vector<std::vector<std::size_t>> comps;
            std::unordered_map<std::size_t, std::size_t> root_to_comp;
            for (std::size_t i = 0; i < cands.size(); ++i) {
                const std::size_t r = dsu.find(i);
                auto it = root_to_comp.find(r);
                if (it == root_to_comp.end()) { root_to_comp[r] = comps.size(); comps.push_back({i}); }
                else comps[it->second].push_back(i);
            }
            for (const std::vector<std::size_t>& comp : comps) {
                std::size_t best = comp.front();
                for (std::size_t k : comp)
                    if (cands[k].e->size_bp > cands[best].e->size_bp) best = k;
                MergedRecord mr;
                mr.seed = *cands[best].e;                 // largest member represents the record
                mr.min_size_bp = mr.seed.size_bp;
                mr.max_size_bp = mr.seed.size_bp;
                std::string link = mr.seed.link_id;
                for (std::size_t k : comp) {
                    const Event& e = *cands[k].e;
                    mr.min_size_bp = std::min(mr.min_size_bp, e.size_bp);
                    mr.max_size_bp = std::max(mr.max_size_bp, e.size_bp);
                    mr.merge_max_jaccard = std::max(mr.merge_max_jaccard, best_jac[k]);
                    mr.merge_max_seqid = std::max(mr.merge_max_seqid, best_seq[k]);
                    if (link.empty() && !e.link_id.empty()) link = e.link_id;
                    mr.member_nodes.insert(e.nodes.begin(), e.nodes.end());
                    if (mr.member_alleles.insert(cands[k].ai).second)
                        for (const std::string& m : alleles[cands[k].ai].members)
                            mr.carriers.push_back(m);
                }
                if (mr.seed.link_id.empty()) mr.seed.link_id = link;
                merged.push_back(std::move(mr));
            }
        }

        // ---- Graph-level force-call: interrogate every non-member haplotype at each called
        // locus and add it as a carrier if its OWN walk-diff supports the variant, even when
        // its event fell below the size threshold. allele_events keeps every event the walk
        // produced (it is filtered only when SEEDING records), so this re-reads each
        // haplotype's actual graph traversal against the record's representative — the fixed
        // largest member from the CC pass. Carrier-adding is monotone (the representative
        // never changes), so one pass reaches the fixpoint; no re-alignment, no re-iteration.
        // (A pure node-set containment test was tried and rejected: on folded paralog loci
        // the inserted nodes are shared across nearly all haplotype walks, so it force-called
        // ~every haplotype as a carrier. Requiring the haplotype's own diff to register the
        // event keeps the reference-relative meaning intact.)
        // Parallel over records: each MergedRecord is interrogated against every non-member allele
        // independently and mutates only itself, so this is safe and order-independent. This is the
        // other hot path on big folded bubbles -- events_match -> seq_identity -> fit_align (banded
        // O(L^2)) on large CNV event sequences, run records x alleles times.
        run_parallel(merged.size(), options.threads, [&](std::size_t ri) {
            MergedRecord& mr = merged[ri];
            if (mr.seed.type == EvType::Dup) return;
            for (std::size_t ai = 0; ai < alleles.size(); ++ai) {
                if (mr.member_alleles.count(ai)) continue;
                for (const Event& e : allele_events[ai]) {
                    double jac = -1.0, sid = -1.0;
                    if (events_match(mr.seed, e, &jac, &sid)) {
                        mr.member_alleles.insert(ai);
                        mr.min_size_bp = std::min(mr.min_size_bp, e.size_bp);
                        mr.max_size_bp = std::max(mr.max_size_bp, e.size_bp);
                        mr.merge_max_jaccard = std::max(mr.merge_max_jaccard, jac);
                        mr.merge_max_seqid = std::max(mr.merge_max_seqid, sid);
                        mr.member_nodes.insert(e.nodes.begin(), e.nodes.end());
                        for (const std::string& m : alleles[ai].members) mr.carriers.push_back(m);
                        break;
                    }
                }
            }
        });

        // ---- De-duplicate a folded duplication against its peak-multiplicity DUP. The
        // walk-diff and the peak DUP both see the extra copy: the DUP as a count, the
        // walk-diff as an INS of the duplicated content. (The extra copy may traverse
        // paralogous node ids that are not on this reference path, so it is not recognizable
        // as "reused reference nodes" — instead it is keyed by being carried by exactly the
        // DUP's carriers and spanning a comparable number of bp.) When a CN route emitted the DUP,
        // drop the matching (duplicated-content) INS so the extra copy is reported once, as the DUP.
        // Genuine novel insertions (carriers exceed the DUP's, or a very different size) are kept.
        if (options.cn) {
            std::vector<std::pair<std::unordered_set<std::string>, std::size_t>> peak_dups;
            for (const MergedRecord& mr : merged) {
                if (mr.seed.type == EvType::Dup && mr.seed.cn_peak) {
                    peak_dups.emplace_back(
                        std::unordered_set<std::string>(mr.carriers.begin(), mr.carriers.end()),
                        mr.seed.size_bp);
                }
            }
            if (!peak_dups.empty()) {
                std::vector<MergedRecord> kept;
                for (MergedRecord& mr : merged) {
                    bool drop = false;
                    if (mr.seed.type == EvType::Ins && !mr.carriers.empty() && mr.seed.size_bp > 0) {
                        for (const auto& [carriers, dup_bp] : peak_dups) {
                            bool subset = true;
                            for (const std::string& c : mr.carriers)
                                if (!carriers.count(c)) { subset = false; break; }
                            const double ratio = static_cast<double>(mr.seed.size_bp) /
                                                 static_cast<double>(dup_bp == 0 ? 1 : dup_bp);
                            if (subset && ratio >= 0.5 && ratio <= 2.0) { drop = true; break; }
                        }
                    }
                    if (!drop) kept.push_back(std::move(mr));
                }
                merged = std::move(kept);
            }
        }

        // ---- INS subtype refinement on the representative only (bounded minimap2 calls).
        // Serialized across bubbles: minimap2 is invoked per record and this path is off by
        // default, so guarding it keeps the parallel bubble loop safe at no cost to the common case.
        if (options.classify_ins) {
            static std::mutex minimap_mtx;
            std::lock_guard<std::mutex> minimap_lock(minimap_mtx);
            std::string ref_window;
            for (MergedRecord& mr : merged) {
                if (mr.seed.type != EvType::Ins || mr.seed.seq.empty()) continue;
                if (ref_window.empty()) ref_window = spell_path_steps_sequence(graph, ref_steps);
                const Minimap2Hit hit = minimap2_best_hit(
                    "ins", mr.seed.seq, "refwin", ref_window, options.minimap_preset, options.minimap_best_n);
                mr.seed.ins_subtype = (hit.ok && hit.identity() >= options.ins_dup_min_identity &&
                                       hit.query_end_bp - hit.query_start_bp >= mr.seed.seq.size() / 2)
                                          ? "DUP" : "NOVEL";
            }
        }

        // ---- Keep records reaching min_sv_bp (by representative) and --min-haplotypes.
        {
            // Substitution arms share an EVENTID (link_id); size-filter them as a UNIT so a co-located
            // DEL+INS is kept or dropped together (the event size is the larger arm) and EVENTID is never
            // left orphaned -- e.g. a 26 bp DEL + 50 bp INS replacement is kept whole, not split.
            std::unordered_map<std::string, std::size_t> link_max_bp;
            for (const MergedRecord& mr : merged)
                if (!mr.seed.link_id.empty())
                    link_max_bp[mr.seed.link_id] = std::max(link_max_bp[mr.seed.link_id], mr.seed.size_bp);
            std::vector<MergedRecord> kept;
            for (MergedRecord& mr : merged) {
                std::unordered_set<std::string> uniq(mr.carriers.begin(), mr.carriers.end());
                std::size_t support = uniq.size();
                if (mr.seed.type == EvType::Dup) {
                    support = 0;
                    for (const auto& [s, cn] : mr.sample_cn) { (void)s; if (cn != mr.seed.ref_cn) ++support; }
                } else {
                    const std::size_t eff_bp = mr.seed.link_id.empty()
                        ? mr.seed.size_bp : link_max_bp[mr.seed.link_id];
                    if (eff_bp < options.min_sv_bp) continue;
                }
                if (support < options.min_haplotypes) continue;
                // MAF filter over traversing haplotypes (AN). Drops common/rare by frequency,
                // complementing the count-based --min-haplotypes.
                if (options.min_maf > 0.0) {
                    const std::size_t an = traverses.size();
                    if (an == 0 || static_cast<double>(support) / static_cast<double>(an) < options.min_maf) continue;
                }
                kept.push_back(std::move(mr));
            }
            merged = std::move(kept);
        }

        // ---- enforce the EVENTID contract: a link_id must denote a genuine co-located DEL+INS pair.
        // A substitution arm can vanish (a sub-threshold micro-arm that never became a record) or a link_id
        // can have smeared onto a large merged cluster during merging. Strip any link_id that does not have
        // BOTH a DEL and an INS among the survivors, so EVENTID is never orphaned/spurious. A now-orphaned
        // arm that ONLY passed the size filter via the substitution unit-rule (its size_bp < min_sv_bp, kept
        // because its partner was >= min_sv_bp) no longer qualifies once the pair is broken, so it is dropped
        // -- otherwise a lone sub-threshold arm would surface as a bare INS/DEL.
        {
            std::unordered_map<std::string, std::pair<bool, bool>> link_arms;  // link -> (has_del, has_ins)
            for (const MergedRecord& mr : merged) {
                if (mr.seed.link_id.empty()) continue;
                auto& a = link_arms[mr.seed.link_id];
                if (mr.seed.type == EvType::Del) a.first = true;
                else if (mr.seed.type == EvType::Ins) a.second = true;
            }
            std::vector<MergedRecord> kept;
            kept.reserve(merged.size());
            for (MergedRecord& mr : merged) {
                if (!mr.seed.link_id.empty()) {
                    const auto& a = link_arms[mr.seed.link_id];
                    if (!(a.first && a.second)) {           // orphaned pair -> strip the link
                        mr.seed.link_id.clear();
                        if (mr.seed.type != EvType::Dup && mr.seed.size_bp < options.min_sv_bp)
                            continue;                        // and drop it if it only survived via the pair
                    }
                }
                kept.push_back(std::move(mr));
            }
            merged = std::move(kept);
        }

        // ---- coverage-CN bubble: drop copy-number-LOSS DELs (redundant with the coverage DUP).
        // A DEL spanning >= half a copy unit is a haplotype carrying fewer folded copies -- the same fact
        // the coverage DUP already reports as a reduced per-sample CN -- so emitting it too would
        // double-count the event as both a low CN and a big DEL. Sequence-novel INS/INV, small/local DELs
        // (< half a unit), and genuine substitution arms (link_id still set after the contract check above)
        // are kept, so a rare gene-unit insertion in the same module is never lost. Runs after the EVENTID
        // contract pass so an orphaned (link-cleared) DEL is correctly treated as a lone CN-loss DEL.
        if (coverage_fired && cn_unit_bp > 0.0) {
            std::vector<MergedRecord> kept;
            kept.reserve(merged.size());
            for (MergedRecord& mr : merged) {
                if (mr.seed.type == EvType::Del && mr.seed.link_id.empty() &&
                    static_cast<double>(mr.seed.size_bp) >= 0.5 * cn_unit_bp)
                    continue;
                kept.push_back(std::move(mr));
            }
            merged = std::move(kept);
        }
        if (merged.empty()) return;
        ++summary.bubbles_with_calls;

        // Map each sample to its allele (for per-carrier sub-walk provenance).
        std::unordered_map<std::string, std::size_t> sample_to_allele;
        for (std::size_t ai = 0; ai < alleles.size(); ++ai) {
            for (const std::string& m : alleles[ai].members) sample_to_allele[m] = ai;
        }

        // Genes overlapping this bubble's reference span (for the per-gene DUP table, resolved by
        // realignment in a post-pass). The GENES INFO field is graph-based (collapsed nodes tag
        // multiple genes); the per-gene COPY NUMBER is resolved later by competitive realignment,
        // which separates collapsed paralogs (CYP2D6 vs 2D7) the graph cannot.
        std::vector<int> bubble_gene_idx;
        if (!node_genes.empty()) {
            std::unordered_set<int> seen;
            std::vector<std::string> bnodes(bubble.inside.begin(), bubble.inside.end());
            bnodes.push_back(bubble.source); bnodes.push_back(bubble.sink);
            for (const std::string& id : bnodes) {
                const auto it = node_genes.find(id);
                if (it == node_genes.end()) continue;
                for (int gi : it->second) if (seen.insert(gi).second) bubble_gene_idx.push_back(gi);
            }
            std::sort(bubble_gene_idx.begin(), bubble_gene_idx.end());
        }

        // ---- Optional multiallelic-locus record (--multiallelic-loci): collapse a bounded
        // locus (e.g. an STR/VNTR) into ONE record with explicit-sequence alleles
        // (REF + ALT1,ALT2,...), GT indexing the allele a sample carries. Skipped (falls back
        // to per-event records) when any allele exceeds --multiallelic-max-bp, so large SVs
        // keep their typed representation. Opt-in; default behavior unchanged.
        // A bubble that yielded a copy-number record (coverage CN, self-loop, or peak-multiplicity DUP)
        // is left typed and NOT collapsed -- otherwise the multiallelic record would silently discard
        // REF_CN/FORMAT:CN. So multiallelic applies only to pure DEL/INS/INV bubbles.
        bool bubble_has_cn = coverage_fired;
        for (const MergedRecord& mr : merged)
            if (mr.seed.type == EvType::Dup) { bubble_has_cn = true; break; }
        if (options.multiallelic_loci && !bubble_has_cn) {
            const auto sit = ref_node_pos.find(bubble.source);
            const auto snode = graph.nodes.find(bubble.source);
            auto interior_seq = [&](const std::vector<PathStep>& steps) -> std::string {
                if (steps.size() < 2) return std::string();
                std::vector<PathStep> inner(steps.begin() + 1, steps.end() - 1);
                return spell_path_steps_sequence(graph, inner);
            };
            bool ok = sit != ref_node_pos.end() && snode != graph.nodes.end() && !snode->second.sequence.empty();
            std::string anchor_base, ref_seq;
            std::vector<std::string> alt_seqs;               // VCF ALT alleles (index 1..)
            std::unordered_map<std::string, int> seq_to_idx; // allele seq -> VCF index (0=ref)
            std::vector<int> allele_vcf_idx(alleles.size(), -1);
            std::size_t apos = 0;
            if (ok) {
                const std::size_t slen = snode->second.sequence.size();
                anchor_base = upper_base(snode->second.sequence[slen - 1]);
                apos = sit->second + (slen > 0 ? slen - 1 : 0);
                ref_seq = anchor_base + interior_seq(ref_steps);
                seq_to_idx[ref_seq] = 0;
                if (ref_seq.size() > options.multiallelic_max_bp) ok = false;
            }
            std::size_t max_alt_delta = 0;
            for (std::size_t ai = 0; ok && ai < alleles.size(); ++ai) {
                const std::string seq = anchor_base + interior_seq(alleles[ai].steps);
                if (seq.size() > options.multiallelic_max_bp) { ok = false; break; }
                auto it = seq_to_idx.find(seq);
                int idx;
                if (it == seq_to_idx.end()) {
                    idx = 1 + static_cast<int>(alt_seqs.size());
                    alt_seqs.push_back(seq);
                    seq_to_idx.emplace(seq, idx);
                    const std::size_t d = seq.size() > ref_seq.size() ? seq.size() - ref_seq.size()
                                                                      : ref_seq.size() - seq.size();
                    max_alt_delta = std::max(max_alt_delta, d);
                } else {
                    idx = it->second;
                }
                allele_vcf_idx[ai] = idx;
            }
            // Only collapse if there is real, large-enough variation; else fall through.
            if (ok && !alt_seqs.empty() && max_alt_delta >= options.min_sv_bp) {
                const std::size_t pos = apos;
                const std::size_t end = pos + ref_seq.size() - 1;
                const std::string id = "bubble" + std::to_string(bubble.id) + "_MULTI_" + bubble.source;
                // per-alt AC over traversing haplotypes
                std::vector<std::size_t> ac(alt_seqs.size(), 0);
                std::size_t an = 0;
                for (const std::string& s : sample_names) {
                    if (!traverses.count(s)) continue;
                    ++an;
                    const auto ait = sample_to_allele.find(s);
                    if (ait == sample_to_allele.end()) continue;
                    const int vi = allele_vcf_idx[ait->second];
                    if (vi >= 1) ++ac[static_cast<std::size_t>(vi - 1)];
                }
                std::ostringstream info;
                info << "BUBBLE_ID=" << bubble.id << ";END=" << end
                     << ";NALLELES=" << (alt_seqs.size() + 1) << ";AN=" << an << ";AC=";
                for (std::size_t k = 0; k < ac.size(); ++k) { if (k) info << ','; info << ac[k]; }
                info << ";SVLEN=";
                for (std::size_t k = 0; k < alt_seqs.size(); ++k) {
                    if (k) info << ',';
                    info << (static_cast<long long>(alt_seqs[k].size()) - static_cast<long long>(ref_seq.size()));
                }
                std::ostringstream row;
                row << ref_meta.chrom << '\t' << pos << '\t' << id << '\t' << ref_seq << '\t';
                for (std::size_t k = 0; k < alt_seqs.size(); ++k) { if (k) row << ','; row << alt_seqs[k]; }
                row << "\t.\t.\t" << info.str() << "\tGT";
                for (const std::string& s : sample_names) {
                    row << '\t';
                    if (!traverses.count(s)) { row << '.'; continue; }
                    const auto ait = sample_to_allele.find(s);
                    row << (ait != sample_to_allele.end() && allele_vcf_idx[ait->second] >= 0
                                ? std::to_string(allele_vcf_idx[ait->second]) : ".");
                }
                row << '\n';
                OutRecord rec; rec.pos = pos; rec.end = end; rec.bubble_id = bubble.id;
                rec.id = id; rec.line = row.str();
                out_records.push_back(std::move(rec));
                ++summary.records_written;
                ++summary.multi;
                return; // locus represented by the multiallelic record; skip per-event emission
            }
        }
        // A carrier's realized sub-walk through the event: its canonical bubble steps
        // between the reference node flanking the event upstream and the first reference
        // node past the event end, as a GFA-style >node / <node string.
        auto carrier_subwalk = [&](const std::vector<PathStep>& steps,
                                   std::size_t epos, std::size_t eend) -> std::string {
            long long up = -1;
            for (std::size_t i = 0; i < steps.size(); ++i) {
                const auto it = ref_node_pos.find(steps[i].node_id);
                if (it != ref_node_pos.end() && static_cast<long long>(it->second) <= static_cast<long long>(epos)) {
                    up = static_cast<long long>(i);
                }
            }
            long long down = -1;
            for (std::size_t i = static_cast<std::size_t>(up < 0 ? 0 : up + 1); i < steps.size(); ++i) {
                const auto it = ref_node_pos.find(steps[i].node_id);
                if (it != ref_node_pos.end() && static_cast<long long>(it->second) >= static_cast<long long>(eend)) {
                    down = static_cast<long long>(i);
                    break;
                }
            }
            const std::size_t lo = up >= 0 ? static_cast<std::size_t>(up) : 0;
            std::size_t hi = down >= 0 ? static_cast<std::size_t>(down)
                                       : (steps.empty() ? 0 : steps.size() - 1);
            if (hi < lo) hi = lo;
            std::string out;
            std::size_t emitted = 0;
            for (std::size_t i = lo; i <= hi && i < steps.size(); ++i) {
                if (emitted >= 2000) { out += "..."; break; }
                out += steps[i].reverse ? '<' : '>';
                out += steps[i].node_id;
                ++emitted;
            }
            return out;
        };

        for (const MergedRecord& mr : merged) {
            const Event& e = mr.seed;
            // POS / REF base from the reference node map.
            std::size_t pos = 1;
            std::string ref_base = "N";
            const auto ait = ref_node_pos.find(e.anchor_node);
            if (ait != ref_node_pos.end()) {
                const std::size_t glen = node_len(graph, e.anchor_node);
                if (e.anchor_after && glen > 0) {
                    pos = ait->second + glen - 1; // last base of preceding ref node
                } else {
                    pos = ait->second;            // first base of first event node
                }
                const auto nit = graph.nodes.find(e.anchor_node);
                if (nit != graph.nodes.end() && !nit->second.sequence.empty()) {
                    const std::size_t bi = (e.anchor_after && glen > 0) ? glen - 1 : 0;
                    ref_base = upper_base(nit->second.sequence[bi]);
                }
            }
            std::size_t end = pos;
            long long svlen = 0;
            const char* svt = ev_svtype(e.type);
            if (e.type == EvType::Del) { svlen = -static_cast<long long>(e.size_bp); end = pos + e.size_bp; }
            else if (e.type == EvType::Ins) { svlen = static_cast<long long>(e.size_bp); end = pos; }
            else if (e.type == EvType::Inv) { svlen = static_cast<long long>(e.size_bp); end = pos + e.size_bp; }
            else if (e.cn_peak) { // peak-multiplicity DUP: size_bp is the duplicated content bp
                svlen = static_cast<long long>(e.size_bp);
                end = pos + e.size_bp;
            }
            else { // self-loop DUP: copies x unit length (signed: CN gain or loss vs reference)
                svlen = (static_cast<long long>(e.alt_cn) - static_cast<long long>(e.ref_cn)) *
                        static_cast<long long>(node_len(graph, e.nodes.empty() ? std::string() : e.nodes.front()));
                end = pos + node_len(graph, e.nodes.empty() ? std::string() : e.nodes.front());
            }

            // Ordered, deduplicated event node set: reference nodes by genomic position,
            // haplotype-only nodes kept in representative-walk order after them. This gives
            // a readable START->END progression instead of the raw merge order.
            std::vector<std::string> ev_nodes;
            {
                std::unordered_set<std::string> seen;
                std::vector<std::pair<long long, std::string>> keyed;
                for (std::size_t k = 0; k < e.nodes.size(); ++k) {
                    const std::string& nd = e.nodes[k];
                    if (!seen.insert(nd).second) continue;
                    const auto it = ref_node_pos.find(nd);
                    const long long key = it != ref_node_pos.end()
                        ? static_cast<long long>(it->second)
                        : (1LL << 60) + static_cast<long long>(k);
                    keyed.emplace_back(key, nd);
                }
                std::stable_sort(keyed.begin(), keyed.end(),
                                 [](const auto& a, const auto& b) { return a.first < b.first; });
                for (auto& p : keyed) ev_nodes.push_back(p.second);
            }
            const std::string start_node = ev_nodes.empty() ? e.start_node : ev_nodes.front();
            const std::string end_node = ev_nodes.empty() ? e.end_node : ev_nodes.back();

            // The describe handoff (variant_nodes.tsv) needs every node any carrier traverses for
            // this event, so each carrier's own allele sequence is sketched (not just the
            // representative's). Use the union of all merged member events' nodes, ordered like
            // ev_nodes; DUP records carry only the single REP node, so fall back to ev_nodes.
            std::vector<std::string> var_nodes;
            if (mr.member_nodes.empty()) {
                var_nodes = ev_nodes;
            } else {
                std::vector<std::pair<long long, std::string>> keyed;
                long long tail = 1LL << 60;
                for (const std::string& nd : mr.member_nodes) {
                    const auto it = ref_node_pos.find(nd);
                    const long long key = it != ref_node_pos.end()
                        ? static_cast<long long>(it->second)
                        : tail++;
                    keyed.emplace_back(key, nd);
                }
                std::stable_sort(keyed.begin(), keyed.end(),
                                 [](const auto& a, const auto& b) { return a.first < b.first; });
                for (auto& p : keyed) var_nodes.push_back(p.second);
            }

            std::unordered_set<std::string> carrier_set(mr.carriers.begin(), mr.carriers.end());

            // Unique variant id.
            const std::string base_id = "bubble" + std::to_string(bubble.id) + "_" + svt + "_" + start_node;
            const int seen_n = ++id_counts[base_id];
            const std::string id = seen_n == 1 ? base_id : base_id + "_" + std::to_string(seen_n);

            std::ostringstream info;
            info << "END=" << end << ";SVTYPE=" << svt << ";SVLEN=" << svlen;
            if (e.type != EvType::Dup && mr.min_size_bp != mr.max_size_bp) {
                const long long sign = svlen < 0 ? -1 : 1;
                const long long a = sign * static_cast<long long>(mr.min_size_bp);
                const long long b = sign * static_cast<long long>(mr.max_size_bp);
                info << ";SVLEN_RANGE=" << std::min(a, b) << "," << std::max(a, b);
                if (mr.max_size_bp > 0) {
                    char buf[32];
                    std::snprintf(buf, sizeof(buf), "%.4g",
                                  static_cast<double>(mr.min_size_bp) / static_cast<double>(mr.max_size_bp));
                    info << ";MERGE_SIZE_RATIO=" << buf;
                }
            }
            // Cross-haplotype merge evidence (only set when ≥2 events were actually merged).
            if (e.type != EvType::Dup && mr.merge_max_jaccard >= 0.0) {
                char buf[32]; std::snprintf(buf, sizeof(buf), "%.4g", mr.merge_max_jaccard);
                info << ";MERGE_JACCARD=" << buf;
            }
            if (e.type != EvType::Dup && mr.merge_max_seqid >= 0.0) {
                char buf[32]; std::snprintf(buf, sizeof(buf), "%.4g", mr.merge_max_seqid);
                info << ";MERGE_SEQID=" << buf;
            }
            const std::size_t an = traverses.size();
            const std::size_t ac = carrier_set.size();
            const double af = an > 0 ? static_cast<double>(ac) / static_cast<double>(an) : 0.0;
            info << ";BUBBLE_ID=" << bubble.id
                 << ";START_NODE=" << start_node << ";END_NODE=" << end_node
                 << ";NMERGED=" << carrier_set.size()
                 << ";AN=" << an << ";AC=" << ac;
            { char buf[32]; std::snprintf(buf, sizeof(buf), "%.4g", af); info << ";AF=" << buf; }
            info << ";EVENT_NODES=";
            for (std::size_t k = 0; k < ev_nodes.size(); ++k) { if (k) info << ','; info << ev_nodes[k]; }
            if (!e.link_id.empty()) info << ";EVENTID=bubble" << bubble.id << "_" << e.link_id;
            if (e.type == EvType::Dup) info << ";REF_CN=" << e.ref_cn;
            if (e.type == EvType::Dup && e.ru_len > 0) info << ";RU_LEN=" << e.ru_len;
            if (!genes.empty()) {
                // Genes the variant touches: its reference event nodes; for a DUP (and as a
                // fallback) the whole bubble's reference span (the folded module).
                std::vector<std::string> gnames = (e.type == EvType::Dup)
                    ? std::vector<std::string>{} : genes_for_nodes(ev_nodes);
                if (gnames.empty()) {
                    std::vector<std::string> bn(bubble.inside.begin(), bubble.inside.end());
                    bn.push_back(bubble.source); bn.push_back(bubble.sink);
                    gnames = genes_for_nodes(bn);
                }
                if (!gnames.empty()) {
                    info << ";GENES=";
                    for (std::size_t k = 0; k < gnames.size(); ++k) { if (k) info << ','; info << gnames[k]; }
                }
            }
            if (e.type == EvType::Ins && !e.ins_subtype.empty()) info << ";INS_SUBTYPE=" << e.ins_subtype;
            if (!e.seq.empty() && e.seq.size() <= 20000) {
                if (e.type == EvType::Ins) info << ";INSSEQ=" << e.seq;
                else if (e.type == EvType::Del) info << ";DELSEQ=" << e.seq;
                else if (e.type == EvType::Inv) info << ";INVSEQ=" << e.seq;
            }

            std::ostringstream row;
            row << ref_meta.chrom << '\t' << pos << '\t' << id << '\t' << ref_base
                << "\t<" << svt << ">\t.\t.\t" << info.str() << "\tGT:CN";
            for (const std::string& s : sample_names) {
                row << '\t';
                if (!traverses.count(s)) { row << ".:."; continue; }
                const bool carrier = carrier_set.count(s) != 0;
                row << (carrier ? "1" : "0") << ':';
                if (e.type == EvType::Dup) {
                    const auto cit = mr.sample_cn.find(s);
                    row << (cit != mr.sample_cn.end() ? std::to_string(cit->second) : std::to_string(e.ref_cn));
                } else {
                    row << '.';
                }
            }
            row << '\n';

            OutRecord rec;
            rec.pos = pos;
            rec.end = end;
            rec.bubble_id = bubble.id;
            rec.id = id;
            rec.line = row.str();
            if (options.write_variant_paths) {
                for (const std::string& s : mr.carriers) {
                    const auto sit = sample_to_allele.find(s);
                    std::string walk;
                    if (sit != sample_to_allele.end()) {
                        walk = carrier_subwalk(alleles[sit->second].steps, pos, end);
                    }
                    rec.prov_lines.push_back(
                        id + '\t' + std::to_string(bubble.id) + '\t' + svt + '\t' + s + "\t1\t" + walk);
                }
            }
            out_records.push_back(std::move(rec));

            if (options.write_variant_paths) {
                std::string nodes;
                for (std::size_t k = 0; k < var_nodes.size(); ++k) { if (k) nodes += ','; nodes += var_nodes[k]; }
                variant_nodes_rows.push_back(
                    id + '\t' + std::to_string(bubble.id) + '\t' + svt + '\t' + nodes);
            }

            // Record a per-gene-CN target for this DUP (resolved by realignment post-pass). Carry the
            // per-haplotype module CN so unreliable (collapsed) gene groups can report the total.
            if (e.type == EvType::Dup && !bubble_gene_idx.empty()) {
                dup_targets.push_back(DupGeneTarget{bubble.id, id, bubble_gene_idx, mr.sample_cn});
            }

            ++summary.records_written;
            if (e.type == EvType::Del) ++summary.del;
            else if (e.type == EvType::Ins) ++summary.ins;
            else if (e.type == EvType::Inv) ++summary.inv;
            else ++summary.dup;
        }
    };
    for (std::size_t bubble_idx = 0; bubble_idx < bubbles.size(); ++bubble_idx) process_bubble(bubble_idx);

    // Merge per-bubble outputs in bubble order: deterministic regardless of thread count.
    for (std::size_t bi = 0; bi < bouts.size(); ++bi) {
        call_progress.tick();
        BubbleOut& bo = bouts[bi];
        summary.bubbles_seen += bo.sum.bubbles_seen;
        summary.bubbles_with_reference += bo.sum.bubbles_with_reference;
        summary.bubbles_with_calls += bo.sum.bubbles_with_calls;
        summary.records_written += bo.sum.records_written;
        summary.del += bo.sum.del;
        summary.ins += bo.sum.ins;
        summary.inv += bo.sum.inv;
        summary.dup += bo.sum.dup;
        summary.multi += bo.sum.multi;
        for (OutRecord& r : bo.records) out_records.push_back(std::move(r));
        for (std::string& s : bo.variant_nodes) variant_nodes_rows.push_back(std::move(s));
        for (DupGeneTarget& t : bo.dup_targets) dup_targets.push_back(std::move(t));
    }

    // ---- Coordinate-sort all records and write the (indexable) region + per-bubble VCFs.
    std::stable_sort(out_records.begin(), out_records.end(),
                     [](const OutRecord& a, const OutRecord& b) {
                         if (a.pos != b.pos) return a.pos < b.pos;
                         if (a.end != b.end) return a.end < b.end;
                         return a.id < b.id;
                     });

    std::ofstream region_out(options.out_prefix + ".region.vcf");
    if (!region_out) {
        throw std::runtime_error("Failed to write region VCF: " + options.out_prefix + ".region.vcf");
    }
    write_vcf_header(region_out);

    std::map<std::size_t, std::ofstream> bubble_files;
    for (const OutRecord& rec : out_records) {
        region_out << rec.line;
        if (options.write_per_bubble_vcf) {
            auto& f = bubble_files[rec.bubble_id];
            if (!f.is_open()) {
                const std::string path = options.out_prefix + ".bubble_" + std::to_string(rec.bubble_id) + ".vcf";
                f.open(path);
                if (!f) throw std::runtime_error("Failed to write bubble VCF: " + path);
                write_vcf_header(f);
            }
            f << rec.line;
        }
    }

    if (options.write_variant_paths) {
        std::ofstream prov_out(options.out_prefix + ".variant_paths.tsv");
        if (!prov_out) {
            throw std::runtime_error("Failed to write variant paths TSV: " + options.out_prefix + ".variant_paths.tsv");
        }
        prov_out << "variant_id\tbubble_id\tsvtype\tsample\tgt\tsub_walk\n";
        for (const OutRecord& rec : out_records) {
            for (const std::string& pl : rec.prov_lines) prov_out << pl << '\n';
        }

        // variant_nodes.tsv: per-variant node set, the bridge for `describe --variant-nodes`
        // (restrict k-mer features to nodes participating in called variation).
        std::ofstream vn_out(options.out_prefix + ".variant_nodes.tsv");
        if (!vn_out) {
            throw std::runtime_error("Failed to write variant nodes TSV: " + options.out_prefix + ".variant_nodes.tsv");
        }
        vn_out << "variant_id\tbubble_id\tsvtype\tnode_ids\n";
        for (const std::string& r : variant_nodes_rows) vn_out << r << '\n';
    }

    // GTF annotation sidecars (independent of --no-variant-paths): node->genes map (consumed by
    // describe/gwas) and the per-gene DUP copy-number table. Per-gene copy number is resolved by
    // PRIVATE-K-MER DOSAGE: each gene's discriminative reference sequence (its merged CDS, where paralogs
    // differ; the gene span when a gene has no CDS) yields a set of canonical k-mers unique to it vs its
    // paralogs, and a haplotype's per-copy count of those k-mers is the gene's copy number -- no
    // per-haplotype alignment. This is what separates collapsed paralogs (CYP2D6 vs 2D7) that share graph
    // nodes; graph multiplicity alone reports only the module total. Evidence (hits / private-set size /
    // dosage) is written per row so every call is auditable.
    if (!genes.empty()) {
        write_node_genes_tsv(options.out_prefix + ".node_genes.tsv", node_genes, genes);

        if (!dup_targets.empty()) {
            // Per-gene resolution (spelling every haplotype + building k-mer sets) is only meaningful when
            // a target folds >=2 genes; a single-gene DUP's copy number is just the module total.
            const bool any_multi = std::any_of(dup_targets.begin(), dup_targets.end(),
                [](const DupGeneTarget& t) { return t.gene_idx.size() >= 2; });

            // Module-total row (single gene, or a gene with no private k-mers). Evidence columns are "."
            // because no per-gene k-mer split was made.
            auto emit_total = [&](const DupGeneTarget& t, const std::string& names, const char* reliable) {
                for (const PathRecord& p : graph.paths) {
                    const auto cit = t.sample_cn.find(p.name);
                    const std::string total = cit != t.sample_cn.end() ? std::to_string(cit->second) : ".";
                    dup_gene_cn_rows.push_back(std::to_string(t.bubble_id) + '\t' + t.variant_id + '\t' +
                                               p.name + '\t' + names + '\t' + total + '\t' + reliable +
                                               "\t.\t.\t.");
                }
            };

            std::vector<std::string> hap_seq;                     // per graph.paths (spelled once)
            std::unordered_map<int, std::string> marker_cache;    // gene index -> discriminative marker
            if (any_multi) {
                const std::string ref_full_seq = spell_path_steps_sequence(graph, ref_path->steps);
                const std::size_t region_start = ref_meta.region_start_1based;
                auto gene_span_of = [&](int gi) -> std::string {
                    const GeneFeature& g = genes[gi];
                    if (g.start_1based < region_start) return std::string();
                    const std::size_t off = g.start_1based - region_start;
                    if (off >= ref_full_seq.size()) return std::string();
                    const std::size_t len = std::min(g.end_1based - g.start_1based + 1, ref_full_seq.size() - off);
                    return ref_full_seq.substr(off, len);
                };
                auto marker_of = [&](int gi) -> std::string {    // merged CDS, else span (CDS-less genes)
                    const GeneFeature& g = genes[gi];
                    if (g.cds.empty()) return gene_span_of(gi);
                    std::string out;
                    for (const auto& iv : g.cds) {
                        if (iv.first < region_start) continue;
                        const std::size_t off = iv.first - region_start;
                        if (off >= ref_full_seq.size()) continue;
                        const std::size_t len = std::min(iv.second - iv.first + 1, ref_full_seq.size() - off);
                        out += ref_full_seq.substr(off, len);
                    }
                    return out;
                };
                for (const DupGeneTarget& t : dup_targets)
                    if (t.gene_idx.size() >= 2)
                        for (int gi : t.gene_idx)
                            if (!marker_cache.count(gi)) marker_cache[gi] = marker_of(gi);
                hap_seq.assign(graph.paths.size(), std::string());
                run_parallel(graph.paths.size(), options.threads, [&](std::size_t pi) {
                    hap_seq[pi] = spell_path_steps_sequence(graph, graph.paths[pi].steps);
                });
            }

            for (const DupGeneTarget& t : dup_targets) {
                if (t.gene_idx.size() < 2) {   // single gene: the module total IS its copy number
                    emit_total(t, genes[t.gene_idx[0]].gene_name, "1");
                    continue;
                }
                std::vector<std::string> markers;
                markers.reserve(t.gene_idx.size());
                for (int gi : t.gene_idx) markers.push_back(marker_cache[gi]);
                // Module copy number per haplotype (the reliable coverage total), used to split a
                // near-identical paralog pair by its per-site allele fraction.
                std::vector<long> total(graph.paths.size(), 0);
                for (std::size_t pi = 0; pi < graph.paths.size(); ++pi) {
                    const auto cit = t.sample_cn.find(graph.paths[pi].name);
                    if (cit != t.sample_cn.end()) total[pi] = static_cast<long>(cit->second);
                }
                const std::vector<std::vector<GeneCnEvidence>> ev = resolve_gene_cn(markers, hap_seq, total);
                for (std::size_t pi = 0; pi < graph.paths.size(); ++pi) {
                    const std::string& sname = graph.paths[pi].name;
                    const std::string prefix =
                        std::to_string(t.bubble_id) + '\t' + t.variant_id + '\t' + sname + '\t';
                    for (std::size_t k = 0; k < t.gene_idx.size(); ++k) {
                        const GeneCnEvidence& e = ev[pi][k];
                        const std::string gname = genes[t.gene_idx[k]].gene_name;
                        if (e.separable) {
                            char dbuf[32];
                            std::snprintf(dbuf, sizeof(dbuf), "%.2f", e.dosage);
                            dup_gene_cn_rows.push_back(
                                prefix + gname + '\t' + std::to_string(e.cn) + "\t1\t" + dbuf + '\t' +
                                std::to_string(e.hits) + '\t' + std::to_string(e.priv_kmers));
                        } else {   // no private k-mers -> indistinguishable from a paralog; report total
                            const auto cit = t.sample_cn.find(sname);
                            const std::string total = cit != t.sample_cn.end() ? std::to_string(cit->second) : ".";
                            dup_gene_cn_rows.push_back(prefix + gname + '\t' + total + "\t0\t.\t.\t0");
                        }
                    }
                }
            }
        }

        std::ofstream dg_out(options.out_prefix + ".dup_gene_cn.tsv");
        if (!dg_out) {
            throw std::runtime_error("Failed to write per-gene DUP CN TSV: " + options.out_prefix + ".dup_gene_cn.tsv");
        }
        dg_out << "bubble_id\tvariant_id\tsample\tgenes\tcn\treliable\tdosage\thits\tpriv_kmers\n";
        for (const std::string& r : dup_gene_cn_rows) dg_out << r << '\n';
    }

    if (summary_out) *summary_out = summary;
    if (!options.quiet) {
        std::cerr << "[call] bubbles=" << summary.bubbles_seen
                  << " with_ref=" << summary.bubbles_with_reference
                  << " with_calls=" << summary.bubbles_with_calls
                  << " records=" << summary.records_written
                  << " (DEL=" << summary.del << " INS=" << summary.ins
                  << " INV=" << summary.inv << " DUP=" << summary.dup
                  << " MULTI=" << summary.multi << ")\n";
    }
}

} // namespace panvar
