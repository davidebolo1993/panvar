#include "panvar/variant_call.hpp"

#include "panvar/align.hpp"
#include "panvar/bubble_alleles.hpp"
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
#include <atomic>
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

// One node-token position in a walk (oriented), for the DEL/INS/INV alignment.
struct Tok {
    std::uint64_t token = 0;
    std::string node_id;
    bool reverse = false;
    // Index of this token in the walk it was collapsed from. For the reference walk that resolves,
    // via the allele set's step indices, to a position in the reference path -- and therefore to the
    // OCCURRENCE the event actually sits at, which a node id alone cannot identify.
    std::size_t src_idx = 0;
};

// Token walk for the DEL/INS/INV alignment. Copy-number nodes (REP self-loops)
// are dropped entirely — they are handled separately as count-based DUP/CN events,
// so a haplotype's extra REP copies are never mistyped as INS.
std::vector<Tok> collapse_walk(
    const std::vector<PathStep>& steps,
    const std::unordered_set<std::string>& cn_nodes) {

    std::vector<Tok> out;
    for (std::size_t i = 0; i < steps.size(); ++i) {
        const PathStep& s = steps[i];
        if (cn_nodes.count(s.node_id) != 0) {
            continue;   // dropped, so the token index is NOT the step index: carry it explicitly
        }
        Tok t;
        t.token = hash_step_token(s);
        t.node_id = s.node_id;
        t.reverse = s.reverse;
        t.src_idx = i;
        out.push_back(std::move(t));
    }
    return out;
}

// A typed event derived from one haplotype walk vs the reference walk.
// How a DUP record's CN was obtained. Three genuinely different measurements, and a consumer cannot
// interpret CN, RU_LEN or SVLEN without knowing which one it is looking at:
//   Rep       a panphorte REP self-loop -- an exact traversal count of a literal repeat unit
//   ModuleBp  folded-node bp divided by a reference-calibrated unit; the module may hold several
//             paralogs, so its "unit" is the SHARED per-copy content, not a whole copy
//   Peak      the highest interior-node traversal multiplicity, used where nothing folded
enum class CnMethod { None, Rep, ModuleBp, Peak };

const char* cn_method_name(CnMethod m) {
    switch (m) {
        case CnMethod::Rep: return "REP";
        case CnMethod::ModuleBp: return "MODULE_BP";
        case CnMethod::Peak: return "PEAK";
        default: return "";
    }
}

// Both of these count copies of a collapsed module rather than of a literal repeat unit, and every
// rule that is about "a module DUP rather than a REP DUP" must cover both. The single boolean these
// two shared before the split hid that distinction, and every use of it meant this predicate.
bool is_module_cn(CnMethod m) { return m == CnMethod::ModuleBp || m == CnMethod::Peak; }

// What the copy is a copy OF. REPEAT_UNIT means CN counts a literal repeat unit; COLLAPSED_MODULE
// means it counts copies of a module that may contain several distinct paralogs, so per-copy content
// is not one uniform sequence and (CN - REF_CN) x RU_LEN is not the haplotype's event size.
const char* cn_scope_name(CnMethod m) {
    switch (m) {
        case CnMethod::Rep: return "REPEAT_UNIT";
        case CnMethod::ModuleBp:
        case CnMethod::Peak: return "COLLAPSED_MODULE";
        default: return "";
    }
}

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
    // Positions in the bubble's REFERENCE walk that this event occupies (SIZE_MAX = none). These
    // identify the occurrence, which anchor_node alone cannot: a reference that revisits a node gives
    // the same name to two different places, and every node->position map records only the first.
    std::size_t ref_tok_first = SIZE_MAX;   // first affected reference token (DEL/INV)
    std::size_t ref_tok_last = SIZE_MAX;    // last affected reference token (DEL/INV)
    std::size_t ref_tok_anchor = SIZE_MAX;  // last MATCHED reference token before the event (INS)
    std::string anchor_node;          // ref node POS is taken from (the walk-order flank)
    bool anchor_after = false;        // true: POS = last base of anchor (INS); false: first base
    std::size_t size_bp = 0;          // |event| for min_sv filtering / SVLEN magnitude
    std::size_t ref_pos = 0;          // reference genomic position of the anchor (merge window)
    // Which of the three CN routes produced this record. A boolean could only say "peak or not", so
    // the coverage and peak routes were indistinguishable in the output and in the code that reads it
    // -- and they answer different questions: one counts copies of a literal repeat unit, the other
    // counts copies of a collapsed paralog module whose per-copy content is not one uniform sequence.
    CnMethod cn_method = CnMethod::None;
    std::size_t ru_len = 0;           // DUP only: repeat-unit length in bp (RU_LEN; one copy)
    // Instrumentation for the module routes, so the unit they calibrate against is inspectable rather
    // than inferred. Zero when the route did not set them.
    std::size_t shared_fold_bp = 0;   // reference bp in the folded (revisited) node set
    std::size_t ref_fold = 0;         // how many times the reference revisits that set
    std::size_t module_ref_bp = 0;    // reference bp across the whole bubble interior
    double fold_residual = -1.0;      // MODULE_BP: spread of folded bp around ref_fold (-1 = n/a)
    double max_support = -1.0;        // MODULE_BP: share of folded bp at ref_fold, the anchor's support
    double step_bp = 0.0;             // MODULE_BP: bp one copy adds, from the panel's cluster spacing
    bool dosage_spacing = false;      // MODULE_BP: CN came from the spacing model, not hbp/unit
    double round_residual = -1.0;     // MODULE_BP: mean distance from a whole number of units
    double round_ambiguous_frac = -1.0; // MODULE_BP: share of traversers rounding near a coin flip
    std::size_t cn_clamped_zero = 0;  // MODULE_BP: traversers whose modelled dosage was below zero
    std::size_t module_span_ambiguous = 0;  // traversers whose module span had a repeated boundary
    // What the spacing estimate actually rests on. A step is "one copy" only if the clusters it was
    // measured between are ADJACENT copy states, and that is an assumption the estimator cannot verify
    // from within -- these let a reader see how thin the evidence is instead of taking the number.
    std::size_t step_clusters = 0;    // clusters the panel's walk lengths fell into
    std::size_t step_gaps = 0;        // gaps the median was taken over (clusters kept, minus one)
    std::size_t step_dropped = 0;     // haplotypes in singleton clusters, which are discarded
    double step_max_offint = -1.0;    // worst |gap/step - nearest integer|; large means the clusters
                                      // are not evenly spaced, so they are not one copy apart
    double step_max_multiple = 0.0;   // largest gap/step; >= 2 means a copy state is missing
};

// Spell a token run in place (one reservation, no per-node temporaries): hot path under the
// parallel allele loop, so the allocator would otherwise dominate.
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
// Segments abandoned because the node-token DP would exceed its cell cap. Each one is a divergent
// block that produced NO variant records, so it is a false negative -- and one that used to vanish
// unless PANVAR_CALL_DEBUG happened to be set. Counted here and reported at the end of the run.
//
// Process-global because diff_segment runs under the per-bubble parallelism and threading a counter
// through it would touch every frame between. call_variants resets it on entry, so the figure is
// per-run rather than cumulative -- without that, a second call in one process reports the first
// call's skips as its own. That reset assumes call_variants is not run concurrently with itself in
// one process, which is true of every command path here.
std::atomic<std::size_t> g_skipped_segments{0};

void diff_segment(
    const Graph& graph,
    const std::vector<Tok>& R, std::size_t r0, std::size_t r1,
    const std::vector<Tok>& H, std::size_t h0, std::size_t h1,
    const std::string& preceding_ref_node,
    std::size_t preceding_ref_tok,
    std::vector<Event>& events) {

    const std::size_t m = r1 - r0;
    const std::size_t n = h1 - h0;
    if (m == 0 && n == 0) return;
    if (m * n > kAlignCellCap) {
        if (std::getenv("PANVAR_CALL_DEBUG")) {
            std::cerr << "[diff] SKIP segment m=" << m << " n=" << n << " (cap " << kAlignCellCap << ")\n";
        }
        ++g_skipped_segments;   // counted and reported, not silently dropped
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
    // The anchor for an INS that OPENS this segment is the shared anchor the previous segment ended
    // on, so its token index has to arrive with it. Starting at SIZE_MAX left every such insertion
    // without an occurrence-aware anchor -- which is most of them, since segments are cut at anchors.
    std::size_t last_ref_tok = preceding_ref_tok;
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
            e.ref_tok_first = ref_blk.front()->src_idx;
            e.ref_tok_last = ref_blk.back()->src_idx;
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
                e.ref_tok_first = ref_blk.front()->src_idx;
                e.ref_tok_last = ref_blk.back()->src_idx;
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
                e.ref_tok_anchor = last_ref_tok;
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
            last_ref_tok = R[static_cast<std::size_t>(c.ri)].src_idx;
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
    std::size_t preceding_ref_tok = SIZE_MAX;   // no token precedes the first segment
    for (std::size_t a = 0; a <= chain_ref.size(); ++a) {
        const std::size_t r_end = (a < chain_ref.size()) ? chain_ref[a] : m;
        const std::size_t h_end = (a < chain_hap.size()) ? chain_hap[a] : n;
        diff_segment(graph, R, r_prev, r_end, H, h_prev, h_end, preceding_ref, preceding_ref_tok, events);
        if (a < chain_ref.size()) {
            preceding_ref = R[r_end].node_id; // the anchor node
            preceding_ref_tok = R[r_end].src_idx;
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
    // Haplotype-space bp span of the event's own nodes: defined for INS; DEL/INV nodes are
    // reference-only (absent from hap_node_pos) so this is empty and the merge falls back to
    // reference coords. Catches same-type events far apart on the reference but contiguous in the sample.
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
            // The reference walk span has to grow with the content. Extending size_bp while leaving
            // the span at the first block made the record claim more deleted bases than the interval
            // it pointed at contained -- 85 bp of sequence over a 58 bp span at LPA bubble 8.
            if (e.ref_tok_first != SIZE_MAX)
                prev.ref_tok_first = prev.ref_tok_first == SIZE_MAX
                    ? e.ref_tok_first : std::min(prev.ref_tok_first, e.ref_tok_first);
            if (e.ref_tok_last != SIZE_MAX)
                prev.ref_tok_last = prev.ref_tok_last == SIZE_MAX
                    ? e.ref_tok_last : std::max(prev.ref_tok_last, e.ref_tok_last);
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
    // MODULE_BP: the raw hbp/unit before rounding. A record-level mean residual hides which samples are
    // ambiguous; this keeps the question per sample, where the answer is actionable.
    std::unordered_map<std::string, double> sample_dosage;
    std::size_t min_size_bp = 0;                            // smallest merged member size (SVLEN_RANGE)
    std::size_t max_size_bp = 0;                            // largest merged member size
    double merge_max_jaccard = -1.0;                        // strongest node-set Jaccard that merged a member (-1 = none)
    double merge_max_seqid = -1.0;                          // strongest sequence identity that merged a member (-1 = none)
    // Weakest pairwise node Jaccard inside the merged component (-1 = single member). Single-linkage
    // can chain A-B-C where A and C would never merge directly, and the strongest-edge fields above
    // cannot show that; this can. `exact` is false when the component was too large to compare all
    // pairs and the figure is an upper bound taken against the representative.
    double merge_diameter = -1.0;
    bool merge_diameter_exact = true;
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

    // Per-run, not cumulative: see the declaration.
    g_skipped_segments.store(0);

    // Refuse a graph this caller cannot describe truthfully, before any of it is read. `spell()` skips
    // a step whose node is absent, so an incomplete graph does not fail -- it produces shorter
    // sequences, and every coordinate, SVLEN and allele derived from them is confidently wrong. Call
    // spells by concatenating whole nodes, so a non-zero overlap would be double-counted the same way.
    validate_graph_paths(graph, "call", /*require_sequences=*/true, /*require_zero_overlaps=*/true);

    // Exact name wins, else a unique case-insensitive substring. Shared with bubble and inspect rather
    // than reimplemented: this rule previously lived in two places and only one of them was fixed.
    const std::string ref_name = resolve_reference_path_name(graph, options.reference_path, "call");
    const PathRecord* ref_path = nullptr;
    for (const PathRecord& p : graph.paths) {
        if (p.name == ref_name) { ref_path = &p; break; }
    }
    if (ref_path == nullptr) throw std::runtime_error("call: reference path not found: " + ref_name);

    const std::vector<Bubble> bubbles = read_bubbles_csv(options.bubbles_csv_in);
    if (bubbles.empty()) {
        throw std::runtime_error("call: bubbles CSV has no rows: " + options.bubbles_csv_in);
    }
    // The CSV and the graph must be the same pair. A CSV from another graph names nodes that do not
    // exist here; every such bubble is then skipped and the run still exits 0 with a header-only VCF,
    // which is indistinguishable from a locus with no variation.
    {
        std::unordered_set<std::size_t> seen_ids;
        for (const Bubble& b : bubbles) {
            if (!seen_ids.insert(b.id).second) {
                throw std::runtime_error(
                    "call: bubbles CSV has a duplicate bubble id: " + std::to_string(b.id) +
                    " (each id must name one site; the same bubble would otherwise be called twice)");
            }
            for (const std::string& n : {b.source, b.sink}) {
                if (graph.nodes.find(n) == graph.nodes.end()) {
                    throw std::runtime_error(
                        "call: bubble " + std::to_string(b.id) + " boundary node '" + n +
                        "' is not in the graph -- the bubbles CSV does not belong to this GFA");
                }
            }
            for (const std::string& n : b.inside) {
                if (graph.nodes.find(n) == graph.nodes.end()) {
                    throw std::runtime_error(
                        "call: bubble " + std::to_string(b.id) + " interior node '" + n +
                        "' is not in the graph -- the bubbles CSV does not belong to this GFA");
                }
            }
        }
    }
    std::unordered_set<std::size_t> bubble_filter(options.bubble_ids.begin(), options.bubble_ids.end());
    // A requested id that does not exist is a typo, not an empty result: without this the run reports
    // success over zero bubbles and writes a header-only VCF.
    if (!bubble_filter.empty()) {
        std::unordered_set<std::size_t> have;
        for (const Bubble& b : bubbles) have.insert(b.id);
        std::vector<std::size_t> missing;
        for (const std::size_t want : bubble_filter) if (!have.count(want)) missing.push_back(want);
        if (!missing.empty()) {
            std::sort(missing.begin(), missing.end());
            std::string msg = "call: --bubble-id names " + std::to_string(missing.size()) +
                              " id(s) not in " + options.bubbles_csv_in + ":";
            for (const std::size_t m : missing) msg += " " + std::to_string(m);
            throw std::runtime_error(msg);
        }
    }

    const ParsedReferencePath ref_meta = parse_reference_path_label(ref_name);
    const std::vector<std::size_t> ref_prefix = path_prefix_bp(*ref_path, graph.nodes);

    // Last reference coordinate: no record's END may exceed it (the reference is finite). Guards
    // against a copy-number event whose bp estimate runs past the region -- a low-complexity tangle
    // can make the peak/coverage DUP's excess-bp sum balloon past the graph.
    std::size_t ref_total_bp = 0;
    for (const auto& st : ref_path->steps) ref_total_bp += node_len(graph, st.node_id);
    const std::size_t ref_end_1based =
        ref_meta.region_start_1based + (ref_total_bp > 0 ? ref_total_bp - 1 : 0);
    // A peak/coverage DUP spanning more than this many bp is a tangle artifact, not a real duplication.
    const std::size_t max_dup_bp = options.max_dup_region_frac > 0.0
        ? static_cast<std::size_t>(options.max_dup_region_frac * static_cast<double>(ref_total_bp))
        : 0;

    // Reference base at a 1-based genomic coordinate, by binary search over the reference walk's
    // cumulative bp. Deliberately coordinate-driven rather than node-driven: anchoring a DEL/INV on
    // the base BEFORE the event needs that base, and the node holding it is only reachable by walk
    // position -- but every node->position map here records a node's FIRST occurrence, which is the
    // wrong occurrence whenever the reference revisits it. Arithmetic on the coordinate has no such
    // ambiguity. Returns an empty string when the coordinate is outside the region.
    auto ref_base_at = [&](std::size_t genomic_1based) -> std::string {
        if (genomic_1based < ref_meta.region_start_1based) return std::string();
        const std::size_t off = genomic_1based - ref_meta.region_start_1based;
        if (ref_prefix.empty() || off >= ref_prefix.back()) return std::string();
        const std::size_t k =
            static_cast<std::size_t>(std::upper_bound(ref_prefix.begin(), ref_prefix.end(), off) -
                                     ref_prefix.begin()) - 1;
        if (k >= ref_path->steps.size()) return std::string();
        const PathStep& st = ref_path->steps[k];
        const auto nit = graph.nodes.find(st.node_id);
        if (nit == graph.nodes.end() || nit->second.sequence.empty()) return std::string();
        const std::string& s = nit->second.sequence;
        const std::size_t idx = off - ref_prefix[k];
        if (idx >= s.size()) return std::string();
        if (!st.reverse) return upper_base(s[idx]);
        return upper_base(reverse_complement(std::string(1, s[s.size() - 1 - idx]))[0]);
    };

    // Optional GTF annotation: project reference-coordinate genes onto reference nodes via the PanSN
    // chrom+start. node_genes maps a ref node id -> gene indices. Built once, const in the parallel
    // loop; skipped when --gtf is unset or the reference name isn't PanSN.
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
    // name -> path index, for the per-sample FORMAT:CNBP walk-length lookup on DUP records.
    std::unordered_map<std::string, std::size_t> name_to_pi;
    for (std::size_t pi = 0; pi < graph.paths.size(); ++pi) name_to_pi.emplace(graph.paths[pi].name, pi);

    const std::unordered_set<std::string> selfloops = self_loop_nodes(graph);

    // Distinct-neighbour degree per node, for the low-complexity tangle guard. A hub reached from all
    // over the graph has high degree; a repeat/paralog unit sits in a chain (degree ~2-6). One pass.
    std::unordered_map<std::string, std::size_t> node_degree;
    if (options.tangle_min_hubs > 0) {
        for (const auto& kv : graph.nodes) {
            std::unordered_set<std::string> nb;
            for (const Neighbor& n : kv.second.start) nb.insert(n.node_id);
            for (const Neighbor& n : kv.second.end) nb.insert(n.node_id);
            node_degree[kv.first] = nb.size();
        }
    }

    cli::ensure_parent_dir_for_file(options.out_prefix + ".region.vcf");

    auto write_vcf_header = [&](std::ostream& out) {
        out << "##fileformat=VCFv4.2\n";
        out << "##source=panvar call\n";
        out << "##reference=" << ref_name << "\n";
        out << "##contig=<ID=" << ref_meta.chrom << ">\n";
        out << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Last reference base the variant spans, inclusive. POS is the base BEFORE the event (symbolic convention), so the event occupies POS+1..END and END-POS is its reference span: |SVLEN| for DEL/INV, 0 for INS (which spans no reference), and the module's own reference span for a CN_SCOPE=COLLAPSED_MODULE DUP -- there it equals INFO/CN_MODULE_REF_BP, the interval FORMAT:CNBP is measured over\">\n";
        out << "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Structural variant type\">\n";
        out << "##INFO=<ID=SVLEN,Number=A,Type=Integer,Description=\"Length difference ALT-REF. Absent on CN_SCOPE=COLLAPSED_MODULE records: their carriers both gain and lose copies, so no single record-level size exists -- read FORMAT:CNBP for the per-sample size\">\n";
        out << "##INFO=<ID=SVLEN_RANGE,Number=2,Type=Integer,Description=\"Min,max event size among merged members (when they differ)\">\n";
        out << "##INFO=<ID=BUBBLE_ID,Number=1,Type=Integer,Description=\"panvar bubble identifier\">\n";
        out << "##INFO=<ID=START_NODE,Number=1,Type=String,Description=\"First graph node of the event\">\n";
        out << "##INFO=<ID=END_NODE,Number=1,Type=String,Description=\"Last graph node of the event\">\n";
        out << "##INFO=<ID=EVENT_NODES,Number=.,Type=String,Description=\"Variant node set\">\n";
        out << "##INFO=<ID=INS_SUBTYPE,Number=1,Type=String,Description=\"INS subtype: NOVEL or DUP (minimap2 refined)\">\n";
        out << "##INFO=<ID=REF_CN,Number=1,Type=Integer,Description=\"Reference copy number of the repeat unit (DUP)\">\n";
        out << "##INFO=<ID=RU_LEN,Number=1,Type=Integer,Description=\"Repeat-unit length in bp, one copy. Emitted only for CN_METHOD=REP, where the unit is a literal repeat and (CN-REF_CN)*RU_LEN is the haplotype size\">\n";
        out << "##INFO=<ID=CN_UNIT_BP,Number=1,Type=Integer,Description=\"CN_METHOD=MODULE_BP only: the reference-calibrated unit CN was divided by (CN_SHARED_BP/CN_REF_FOLD). This is the SHARED per-copy content, not a whole copy, so it is not a per-haplotype event size -- read FORMAT:CNBP for that\">\n";
        out << "##INFO=<ID=CN_METHOD,Number=1,Type=String,Description=\"How CN was measured: REP (traversal count of a panphorte REP self-loop, exact), MODULE_BP (folded-node bp over a reference-calibrated unit), PEAK (highest interior-node traversal multiplicity)\">\n";
        out << "##INFO=<ID=CN_SCOPE,Number=1,Type=String,Description=\"What one copy is a copy of: REPEAT_UNIT (a literal repeat unit, so (CN-REF_CN)*RU_LEN is the haplotype's event size) or COLLAPSED_MODULE (a module that may hold several distinct paralogs, so RU_LEN is NOT emitted -- the shared per-copy content is reported as CN_UNIT_BP, and per-haplotype size must be read from FORMAT:CNBP)\">\n";
        out << "##INFO=<ID=CN_SHARED_BP,Number=1,Type=Integer,Description=\"Reference bp in the folded (revisited) node set the MODULE_BP unit was calibrated from\">\n";
        out << "##INFO=<ID=CN_REF_FOLD,Number=1,Type=Integer,Description=\"How many times the reference revisits that folded set; the MODULE_BP unit is CN_SHARED_BP/CN_REF_FOLD\">\n";
        out << "##INFO=<ID=CN_MODULE_REF_BP,Number=1,Type=Integer,Description=\"Reference bp across the whole bubble interior. Against CN_SHARED_BP this is the shared-versus-total question: CN_UNIT_BP describes the shared part only, and their ratio is how far (CN-REF_CN)*CN_UNIT_BP understates a carrier's real gain or loss\">\n";
        out << "##INFO=<ID=CN_REF_MULTIPLICITY_HETEROGENEITY,Number=1,Type=Float,Description=\"MODULE_BP only: length-weighted spread of the folded set's reference multiplicities around CN_REF_FOLD. 0 means one coherent unit repeated CN_REF_FOLD times. A diagnostic of graph structure, not a correctness test -- see docs/algorithms/call.md\">\n";
        out << "##INFO=<ID=CN_REF_MAX_SUPPORT,Number=1,Type=Float,Description=\"MODULE_BP only: share of multiplicity-weighted folded bp that actually sits at CN_REF_FOLD, the multiplicity REF_CN was taken from. Low values mean the anchor is supported by a small part of the module\">\n";
        out << "##INFO=<ID=CN_DOSAGE_MODEL,Number=1,Type=String,Description=\"MODULE_BP only: REFERENCE_RATIO (CN = hbp/CN_UNIT_BP, the default) or PANEL_SPACING (CN = REF_CN + (hbp-CN_SHARED_BP)/CN_STEP_BP, --cn-unit-spacing). Determines what FORMAT:CNR_RAW means\">\n";
        out << "##INFO=<ID=CN_STEP_BP,Number=1,Type=Integer,Description=\"MODULE_BP only: bp one copy adds, estimated from the spacing between the panel's own copy-state clusters. Independent of CN_UNIT_BP, which comes from ref_bp/ref_fold\">\n";
        out << "##INFO=<ID=CN_STEP_RATIO,Number=1,Type=Float,Description=\"CN_STEP_BP / CN_UNIT_BP. 1.0 means hbp is proportional to copy number, which is what the MODULE_BP integer assumes. Measured at 1.45 on both reference paralog modules, i.e. the dosage is AFFINE and its error grows about 0.45 per copy either side of the reference copy number -- see --cn-unit-spacing and docs/algorithms/call.md\">\n";
        out << "##INFO=<ID=CN_STEP_SUPPORT,Number=3,Type=Integer,Description=\"What CN_STEP_BP rests on: clusters,gaps,dropped. clusters = copy-state clusters the panel's walk lengths fell into; gaps = differences the median was taken over (one fewer than the populated clusters); dropped = haplotypes in singleton clusters, which are discarded. gaps=1 means the step is a single difference between two clusters, with nothing to cross-check it\">\n";
        out << "##INFO=<ID=CN_STEP_OFFINT,Number=1,Type=Float,Description=\"Worst distance from a whole number of copies among the gaps CN_STEP_BP was taken over (0 = every gap is an exact multiple of the step). Large values mean the clusters are not evenly spaced, so treating a gap as one copy is unsupported\">\n";
        out << "##INFO=<ID=CN_STEP_MAX_MULTIPLE,Number=1,Type=Float,Description=\"Largest gap as a multiple of CN_STEP_BP. At or above 2 a copy state is missing from the panel between two observed clusters; the estimator assumes adjacency and cannot detect this on its own\">\n";
        out << "##INFO=<ID=REF_CN_SOURCE,Number=1,Type=String,Description=\"How REF_CN was anchored: REP_TRAVERSAL (exact self-loop count) or MAX_NODE_MULTIPLICITY (a heuristic -- one short node visited N times can set it, so absolute CN on that route is heuristic even where relative dosage is sound)\">\n";
        out << "##INFO=<ID=CN_CONFIDENCE,Number=1,Type=String,Description=\"HEURISTIC on CN_METHOD=PEAK, which infers dosage from the highest interior-node traversal multiplicity and is not validated against external copy-number truth\">\n";
        out << "##INFO=<ID=CN_ROUND_RESIDUAL,Number=1,Type=Float,Description=\"MODULE_BP only: MEAN distance from a whole number of units across traversers, 0..0.5. Reported for continuity, but the per-sample residuals are BIMODAL at a paralog module -- most near 0 or near 0.5 -- so this mean describes no sample and moves with their ratio. Read CN_ROUND_AMBIGUOUS_FRAC for the record-level summary and FORMAT:CNR_MARGIN per sample\">\n";
        out << "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"The record describes an event CLUSTER rather than one exact event: coalescing joined several pieces that have reference sequence retained between them, so the affected interval POS+1..END is wider than the |SVLEN| bases actually removed. On a PRECISE record END-POS equals |SVLEN| exactly; this flag marks the records where it cannot, instead of silently reporting one of the two numbers\">\n";
        out << "##INFO=<ID=CN_SPAN_AMBIGUOUS,Number=1,Type=Integer,Description=\"Path measurements taken over a module span whose source or sink boundary occurs MORE THAN ONCE in that path. The span is then first-source..last-sink, which is a choice: it can enclose two separate visits and the content between them rather than one module. Every quantity measured over the span -- CN, CNBP, CN_MODULE_REF_BP -- inherits that choice on those paths. Absent when every boundary is visited once\">\n";
        out << "##INFO=<ID=CN_CLAMPED_ZERO,Number=1,Type=Integer,Description=\"MODULE_BP only: traversing haplotypes whose modelled dosage came out BELOW zero copies and whose reported CN was therefore floored at 0. Only the spacing model can produce this, when a walk is shorter than the reference's by more than REF_CN steps. Their CN is a boundary value rather than a measurement; FORMAT:CNR_RAW still carries the unclamped dosage. Absent when none occurred\">\n";
        out << "##INFO=<ID=CN_ROUND_AMBIGUOUS_FRAC,Number=1,Type=Float,Description=\"MODULE_BP only: share of traversing haplotypes whose dosage sat more than 0.4 units from a whole number, i.e. whose integer CN came from a near coin-flip rounding. This is what --max-cn-model-residual gates on. A high value does NOT by itself mean the CN is wrong: the reference locus with the worst value is exact against pangene truth on every haplotype, because the unit is a calibration constant for a heterogeneous module rather than one real copy\">\n";
        if (!genes.empty())
            out << "##INFO=<ID=GENES,Number=.,Type=String,Description=\"Gene(s) overlapping the variant (from --gtf)\">\n";
        out << "##INFO=<ID=NMERGED,Number=1,Type=Integer,Description=\"Haplotype carriers merged into this record\">\n";
        out << "##INFO=<ID=MERGE_DIAMETER,Number=1,Type=Float,Description=\"Weakest pairwise node Jaccard between any two events merged into this record. Merging is transitive single-linkage, so a record can span members that would never have merged directly (A-B-C with A and C dissimilar); MERGE_JACCARD reports the strongest edge and cannot show that, this reports the worst. Near 0 means the record's members share almost no nodes and the chain reached a long way. Absent on unmerged records\">\n";
        out << "##INFO=<ID=MERGE_DIAMETER_EXACT,Number=0,Type=Flag,Description=\"MERGE_DIAMETER compared every pair. Absent means the component exceeded the all-pairs bound and the value is an upper bound measured against the representative only\">\n";
        out << "##INFO=<ID=MERGE_JACCARD,Number=1,Type=Float,Description=\"Strongest node-set Jaccard that merged a member into this record (cross-haplotype merge evidence)\">\n";
        out << "##INFO=<ID=MERGE_SEQID,Number=1,Type=Float,Description=\"Strongest sequence identity that merged a member into this record, when the Jaccard gate did not decide it\">\n";
        out << "##INFO=<ID=MERGE_SIZE_RATIO,Number=1,Type=Float,Description=\"Smallest/largest member size among merged members (min,max also in SVLEN_RANGE)\">\n";
        out << "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Allele number = haplotypes traversing the bubble\">\n";
        out << "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count = carrier haplotypes\">\n";
        out << "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency = AC/AN (over traversing haplotypes)\">\n";
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
        out << "##FORMAT=<ID=CNBP,Number=1,Type=Integer,Description=\"Actual linear bp gained(+)/lost(-) by this haplotype across the copy-number module vs the reference, from the spelled walk length (sum of node length x traversal multiplicity over the bubble source->sink, minus the reference's). Recovers the linear SV size that the folded one-copy RU_LEN does not convey; DUP records only.\">\n";
        out << "##FORMAT=<ID=CNR_RAW,Number=1,Type=Float,Description=\"CN_METHOD=MODULE_BP only: this haplotype's raw dosage before rounding to CN. Which quantity that is depends on INFO/CN_DOSAGE_MODEL: hbp/CN_UNIT_BP under REFERENCE_RATIO, or REF_CN+(hbp-CN_SHARED_BP)/CN_STEP_BP under PANEL_SPACING\">\n";
        out << "##FORMAT=<ID=CNRESID,Number=1,Type=Integer,Description=\"CN_METHOD=REP only: CNBP - (CN-REF_CN)*RU_LEN, the bp this haplotype gained or lost that its repeat-copy change does NOT explain. 0 means the whole size difference is copy number; a large value means other sequence changed in the same bubble. Emitted only for a literal repeat unit, where RU_LEN is a real per-copy length -- on a collapsed module RU_LEN is a calibration constant and the difference would measure that, not biology\">\n";
        out << "##FORMAT=<ID=CNR_MARGIN,Number=1,Type=Float,Description=\"CN_METHOD=MODULE_BP only: 0.5 minus this haplotype's distance from a whole number of units. Near 0.5 the integer CN is unambiguous; near 0 the rounding was a coin flip and CN should be read with that in mind\">\n";
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
    };
    std::vector<OutRecord> out_records;
    std::vector<OutRecord> allele_records;
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
        std::vector<OutRecord> allele_records;
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
    // Bubbles run serially; parallelism is inside each bubble, over its alleles (the walk-diff DP is
    // the hot path, and an SV locus is usually one dominant folded bubble). Each bubble writes only to
    // its own slot, so results are order-deterministic.
    auto process_bubble = [&](std::size_t bubble_idx) {
        const Bubble& bubble = bubbles[bubble_idx];
        // Shadow the shared sinks with this bubble's local buffers: the (unchanged) loop body below
        // writes only through these names, so each bubble accumulates independently. id_counts is
        // bubble-scoped (variant ids embed the bubble id), so a local map reproduces the ids exactly.
        BubbleOut& bout = bouts[bubble_idx];
        VariantCallSummary& summary = bout.sum;
        std::vector<OutRecord>& out_records = bout.records;
        std::vector<OutRecord>& allele_records = bout.allele_records;
        std::vector<std::string>& variant_nodes_rows = bout.variant_nodes;
        std::vector<DupGeneTarget>& dup_targets = bout.dup_targets;
        std::unordered_map<std::string, int> id_counts;
        if (!bubble_filter.empty() && bubble_filter.find(bubble.id) == bubble_filter.end()) return;
        ++summary.bubbles_seen;

        // Reference walk + distinct alleles (grouped by canonical-walk signature).
        std::size_t ref_idx = graph.paths.size();
        for (std::size_t pi = 0; pi < graph.paths.size(); ++pi) {
            if (graph.paths[pi].name == ref_name) { ref_idx = pi; break; }
        }
        BubbleAlleleSet allele_set = enumerate_bubble_alleles(graph, path_indexes, bubble, ref_name);
        if (!allele_set.has_reference) {
            return; // reference does not traverse this bubble; cannot type events
        }
        ++summary.bubbles_with_reference;

        const std::vector<PathStep> ref_steps = std::move(allele_set.reference_steps);
        // Which step of the reference PATH each entry of ref_steps came from. This is what makes a
        // walk position identify an occurrence rather than just a node name.
        const std::vector<std::size_t> ref_step_idx = std::move(allele_set.reference_step_indices);
        // Genomic extent [first base, last base] of a run of reference-walk tokens, by full-path step
        // index rather than by node name. Returns {0,0} when unavailable. Takes min/max over the run
        // because a reverse-oriented bubble walks the reference in DECREASING coordinate, so the
        // walk-order first token is the genomically LAST one.
        auto tok_span = [&](std::size_t lo_tok, std::size_t hi_tok)
            -> std::pair<std::size_t, std::size_t> {
            if (ref_step_idx.empty() || lo_tok == SIZE_MAX || hi_tok == SIZE_MAX) return {0, 0};
            if (lo_tok >= ref_step_idx.size() || hi_tok >= ref_step_idx.size()) return {0, 0};
            if (lo_tok > hi_tok) std::swap(lo_tok, hi_tok);
            std::size_t lo_full = SIZE_MAX, hi_full = 0;
            for (std::size_t j = lo_tok; j <= hi_tok; ++j) {
                lo_full = std::min(lo_full, ref_step_idx[j]);
                hi_full = std::max(hi_full, ref_step_idx[j]);
            }
            if (lo_full == SIZE_MAX || hi_full + 1 >= ref_prefix.size()) return {0, 0};
            return {ref_meta.region_start_1based + ref_prefix[lo_full],
                    ref_meta.region_start_1based + ref_prefix[hi_full + 1] - 1};
        };
        using Allele = BubbleAllele;
        const std::vector<Allele>& alleles = allele_set.alleles;
        const std::unordered_set<std::string>& traverses = allele_set.traversing;
        const std::string& ref_sig = allele_set.reference_signature;

        // CN nodes: self-loop nodes in this bubble (a REP / tandem unit), handled as count-based DUP
        // events and excluded from the DEL/INS/INV alignment. Ordinary recurring nodes stay in it.
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

        // The module's oriented step span in one path, shared by every consumer that measures the
        // module: the folded-set construction, reference and haplotype module bp, CNBP and
        // CN_MODULE_REF_BP. Two copies of this arithmetic existed and could drift; that is the same
        // duplication that let bubble's reference-name rule be fixed in one place and not the other.
        //
        // Deliberately the WIDEST span (first source .. last sink), not the tight allele interval:
        // a module's copies are exactly what lies between the outermost boundaries. `ambiguous`
        // reports when a boundary occurs more than once, because the span is then a CHOICE and can
        // sweep in content between two unrelated visits rather than one module.
        auto module_span = [&](std::size_t pi, bool* ambiguous = nullptr)
            -> std::pair<std::size_t, std::size_t> {
            if (ambiguous != nullptr) *ambiguous = false;
            const auto& idx = path_indexes[pi].positions;
            const auto sit = idx.find(bubble.source);
            const auto kit = idx.find(bubble.sink);
            if (sit == idx.end() || kit == idx.end()) return {1, 0};   // empty (lo > hi)
            if (ambiguous != nullptr && (sit->second.size() > 1 || kit->second.size() > 1))
                *ambiguous = true;
            const std::size_t s0 = sit->second.front(), s1 = sit->second.back();
            const std::size_t k0 = kit->second.front(), k1 = kit->second.back();
            const std::size_t lo = (s0 <= k1) ? s0 : k0;
            const std::size_t hi = (s0 <= k1) ? k1 : s1;
            return {lo, hi};
        };

        const std::size_t rescue_floor =
            options.rescue_min_bp != 0 ? options.rescue_min_bp : std::max<std::size_t>(1, options.min_sv_bp / 2);

        auto ev_ref_pos = [&](const Event& e) -> long long {
            // Occurrence-aware first: this position seeds the merge window and the sort order, so
            // taking a node's FIRST occurrence lets two events at two visits to the same node look
            // co-located and merge into one record.
            const std::size_t tok = e.ref_tok_first != SIZE_MAX ? e.ref_tok_first : e.ref_tok_anchor;
            if (tok != SIZE_MAX) {
                const auto span = tok_span(tok, e.ref_tok_last != SIZE_MAX ? e.ref_tok_last : tok);
                if (span.first > 0)
                    return static_cast<long long>(e.anchor_after ? span.second : span.first);
            }
            const auto it = ref_node_pos.find(e.anchor_node);
            if (it == ref_node_pos.end()) return -1;
            const std::size_t glen = node_len(graph, e.anchor_node);
            return static_cast<long long>(it->second + (e.anchor_after && glen > 0 ? glen - 1 : 0));
        };
        // Two non-DUP events are the same site: same type, anchors within the window, and node-set
        // Jaccard OR sequence identity. out_jac/out_seq receive the evidence (seq computed only when
        // the Jaccard gate didn't already decide it).
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

        // bp-coverage CN for a collapsed paralog cluster: the reference revisits the module >=2x, so CN
        // is node multiplicity, not a tandem block. CN = full-walk bp / unit (unit = ref_full_bp /
        // ref_fold), over the WIDEST source..sink span since the minimal one collapses the repeats.
        // Yields to a genuine REP self-loop; routes stay disjoint (self-loop > coverage > peak).
        bool has_rep_selfloop = false;
        for (const std::string& cn : cn_nodes)
            if (node_len(graph, cn) >= options.min_sv_bp) { has_rep_selfloop = true; break; }
        // Low-complexity tangle: a bubble whose interior holds many high-degree hub nodes (reached from
        // all over the graph) is a low-complexity region, not a copy-number module -- its node
        // multiplicity is meaningless as copy number. Suppress the coverage/peak DUP routes there. A
        // real paralog cluster is chain-like (low degree) and has ~0 hubs, so it is unaffected.
        bool is_tangle = false;
        if (options.tangle_min_hubs > 0) {
            std::size_t hubs = 0;
            for (const std::string& id : bubble.inside) {
                const auto it = node_degree.find(id);
                if (it != node_degree.end() && it->second >= options.tangle_hub_degree &&
                    ++hubs >= options.tangle_min_hubs) { is_tangle = true; break; }
            }
            if (is_tangle) ++summary.tangle_bubbles;
        }
        // Three distinct states. Overloading one flag is what let a DECLINED module CN suppress
        // --multiallelic-loci: "do not try a weaker route" and "a CN record exists" are not the same
        // claim, and the writer needs the second.
        bool cn_route_consumed = false;   // a module route ran; weaker routes must not answer
        bool module_cn_declined = false;  // ...and it refused, so no CN record was emitted
        double cn_unit_bp = 0.0;   // one-copy bp of the coverage module (set when coverage fires)
        if (options.cn && !has_rep_selfloop) {
            const std::unordered_set<std::string> inside_set(bubble.inside.begin(), bubble.inside.end());
            // The shared module resolver, so the CN routes and CNBP measure the SAME interval.
            std::size_t span_ambiguous = 0;
            auto span_of = [&](std::size_t pi) {
                bool amb = false;
                const auto sp = module_span(pi, &amb);
                if (amb) ++span_ambiguous;
                return sp;
            };
            // Folded set = inside nodes the reference revisits (>=2x): the repeat unit. Measuring CN
            // only over these keeps unique-content edits (interstitial / single-visit nodes) out of it.
            std::unordered_set<std::string> folded_set;
            std::size_t ref_fold = 0;
            std::size_t ref_bp = 0;
            double fold_residual = 0.0;   // spread of folded bp around ref_fold (the MAXIMUM, not the mode)
            double max_support = 0.0;     // share of folded bp actually at ref_fold
            double round_residual = 0.0;  // mean |hbp/unit - nearest integer| over traversers, 0..0.5
            std::size_t round_ambiguous = 0;  // traversers whose rounding was near a coin flip (>0.4)
            std::size_t cn_clamped_zero = 0;  // traversers whose modelled dosage came out below zero
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
                // How coherent is "the folded set is ref_fold copies of one unit"? Taking the MAX
                // multiplicity as the fold count assumes every folded node is revisited the same number
                // of times. If they sit at 2 and 3, dividing total folded bp by 3 yields a unit that is
                // not one copy of anything -- and CN is that unit's divisor, so the dosage inherits the
                // incoherence. Measured length-weighted, since one long node at the wrong multiplicity
                // matters more than many short ones.
                std::size_t fit_bp = 0, off_bp = 0, all_bp = 0;
                for (const std::string& id : folded_set) {
                    const std::size_t m = cnt.count(id) ? cnt.at(id) : 0;
                    const std::size_t len = node_len(graph, id);
                    all_bp += len * m;
                    if (m == ref_fold) fit_bp += len * m;
                    else off_bp += len * (m > ref_fold ? m - ref_fold : ref_fold - m);
                }
                // Two numbers, because one alone misleads. The heterogeneity ratio mixes observed bp for
                // matching nodes with distance-to-ref_fold for the rest, so it is a spread measure and
                // not "the fraction of bp at the wrong multiplicity". max_support is the plain question:
                // how much of the folded signal actually sits at the multiplicity REF_CN was taken from.
                fold_residual = (fit_bp + off_bp) > 0
                    ? static_cast<double>(off_bp) / static_cast<double>(fit_bp + off_bp) : 0.0;
                max_support = all_bp > 0 ? static_cast<double>(fit_bp) / static_cast<double>(all_bp) : 0.0;
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
                // The unit from ref_bp/ref_fold assumes hbp is PROPORTIONAL to copy number. Measured
                // against pangene truth it is not: each real copy adds ~1.45 of those units, so the
                // dosage is affine and the error grows ~0.45 per copy either side of the reference.
                // The panel itself carries the missing constant -- haplotypes cluster by copy state,
                // and the gap between adjacent clusters IS one copy. Estimated here and reported; the
                // integer CN still comes from the division route unless --cn-unit spacing is given,
                // because switching the default would move a result validated at 466/466.
                std::vector<std::size_t> all_hbp;
                for (std::size_t pi = 0; pi < graph.paths.size(); ++pi)
                    if (traverses.count(graph.paths[pi].name)) all_hbp.push_back(full_walk_bp(pi));
                double step_bp = 0.0;
                std::size_t step_clusters = 0, step_gaps = 0, step_dropped = 0;
                double step_max_offint = -1.0, step_max_multiple = 0.0;
                if (all_hbp.size() >= 4) {
                    std::sort(all_hbp.begin(), all_hbp.end());
                    std::vector<std::vector<std::size_t>> cl{{all_hbp.front()}};
                    for (std::size_t k = 1; k < all_hbp.size(); ++k) {
                        if (static_cast<double>(all_hbp[k] - cl.back().back()) > 0.15 * unit)
                            cl.emplace_back();
                        cl.back().push_back(all_hbp[k]);
                    }
                    step_clusters = cl.size();
                    std::vector<double> centres;
                    for (const auto& c : cl) {
                        // A cluster of one is dropped: a lone walk length is as likely to be a
                        // mis-folded outlier as a real copy state, and it would drag a centre. The
                        // cost is that a genuinely rare copy state contributes nothing, which is
                        // reported rather than hidden.
                        if (c.size() >= 2) centres.push_back(static_cast<double>(c[c.size() / 2]));
                        else step_dropped += c.size();
                    }
                    std::vector<double> gaps;
                    for (std::size_t k = 1; k < centres.size(); ++k) gaps.push_back(centres[k] - centres[k - 1]);
                    if (!gaps.empty()) {
                        std::vector<double> sorted_gaps = gaps;
                        std::sort(sorted_gaps.begin(), sorted_gaps.end());
                        step_bp = sorted_gaps[sorted_gaps.size() / 2];
                        step_gaps = gaps.size();
                        // Each gap should be a whole number of copies of the chosen step. A gap at 2x
                        // says a copy state is absent from the panel and the estimator cannot see it;
                        // a gap that is not near ANY integer says the clusters are not copy states.
                        for (const double g : gaps) {
                            const double r = step_bp > 0.0 ? g / step_bp : 0.0;
                            step_max_multiple = std::max(step_max_multiple, r);
                            step_max_offint = std::max(step_max_offint, std::fabs(r - std::round(r)));
                        }
                    }
                }
                std::size_t round_n = 0;
                MergedRecord mr;
                mr.seed.type = EvType::Dup;
                mr.seed.nodes.push_back(bubble.source);
                mr.seed.start_node = bubble.source; mr.seed.end_node = bubble.sink;
                mr.seed.ref_cn = ref_copies; mr.seed.alt_cn = ref_copies;
                mr.seed.anchor_node = bubble.source;
                mr.seed.size_bp = static_cast<std::size_t>(unit);
                mr.seed.ru_len = static_cast<std::size_t>(unit);
                mr.seed.cn_method = CnMethod::ModuleBp;
                mr.seed.step_bp = step_bp;
                mr.seed.dosage_spacing = (options.cn_unit_spacing && step_bp > 0.0);
                mr.seed.step_clusters = step_clusters;
                mr.seed.step_gaps = step_gaps;
                mr.seed.step_dropped = step_dropped;
                mr.seed.step_max_offint = step_max_offint;
                mr.seed.step_max_multiple = step_max_multiple;
                // --cn-unit-spacing asked for a specific model. If the panel cannot supply it, the
                // honest outcome is to say so, not to quietly compute the biased ratio model under a
                // flag the user set to avoid exactly that. CN_DOSAGE_MODEL would have recorded the
                // substitution, but only for a reader who thought to check it.
                if (options.cn_unit_spacing) {
                    // Reporting that a model is unsupported while still using it is not enough. The
                    // step is only a copy when the clusters it was measured between are adjacent copy
                    // states and evenly spaced, so each of those is a precondition, not a caveat.
                    std::string why;
                    if (step_bp <= 0.0) {
                        why = std::to_string(all_hbp.size()) + " traversing haplotype(s) fell into " +
                              std::to_string(step_clusters) +
                              " cluster(s), leaving no gap between two populated ones";
                    } else if (step_gaps < 2) {
                        why = "the step rests on a single gap between two clusters, with nothing to "
                              "cross-check that they are one copy apart";
                    } else if (step_max_multiple >= 1.5) {
                        std::ostringstream m; m.setf(std::ios::fixed); m.precision(2);
                        m << step_max_multiple;
                        why = "the widest gap is " + m.str() +
                              "x the chosen step, so at least one copy state is missing from the "
                              "panel between two observed clusters";
                    } else if (step_max_offint > 0.15) {
                        std::ostringstream o; o.setf(std::ios::fixed); o.precision(3);
                        o << step_max_offint;
                        why = "a gap sits " + o.str() +
                              " of a step from any whole number of copies, so the clusters are not "
                              "evenly spaced and are not copy states";
                    }
                    if (!why.empty()) {
                        throw std::runtime_error(
                            "call: --cn-unit-spacing was given but bubble " + std::to_string(bubble.id) +
                            " cannot support it: " + why +
                            ". Rerun without the flag to use the reference-ratio model, which is "
                            "biased but always defined");
                    }
                }
                mr.seed.shared_fold_bp = ref_bp;
                mr.seed.ref_fold = ref_fold;
                mr.seed.fold_residual = fold_residual;
                mr.seed.max_support = max_support;
                for (std::size_t pi = 0; pi < graph.paths.size(); ++pi) {
                    if (!traverses.count(graph.paths[pi].name)) continue;
                    // A traverser with no folded bp carries ZERO copies -- a complete loss of the
                    // module, which is the largest event this record can describe. Skipping it left the
                    // sample out of sample_cn, and the writer then filled the gap with REF_CN, so the
                    // haplotype that lost the whole module read exactly like one that lost nothing.
                    const std::size_t hbp = full_walk_bp(pi);
                    const double exact = (options.cn_unit_spacing && step_bp > 0.0)
                        ? static_cast<double>(ref_copies) +
                          (static_cast<double>(hbp) - static_cast<double>(ref_bp)) / step_bp
                        : static_cast<double>(hbp) / unit;
                    // Spacing mode can put a haplotype BELOW zero copies when its walk is shorter
                    // than the reference's by more than REF_CN steps. llround then yields a negative
                    // long, and the cast to size_t wraps it to something astronomical.
                    const long long rounded = std::llround(exact);
                    const std::size_t copies = rounded > 0 ? static_cast<std::size_t>(rounded) : 0;
                    if (rounded < 0) ++cn_clamped_zero;
                    // How far each haplotype sits from a whole number of units. A calibrated unit should
                    // divide real walks nearly exactly; persistent halves mean the unit is wrong, and
                    // rounding then manufactures a confident integer out of a bad fit.
                    //
                    // Measured against the model's own rounded value, NOT the zero-clamped CN. Using
                    // the clamped one made a negative dosage report a residual above 0.5 and a
                    // CNR_MARGIN below 0, both of which their headers rule out. The clamp is a
                    // reporting floor and belongs in CN_CLAMPED_ZERO, not in the fit statistic.
                    const double resid = std::fabs(exact - static_cast<double>(rounded));
                    round_residual += resid;
                    ++round_n;
                    // Count the samples whose rounding was close to a coin flip. The mean above cannot
                    // stand in for this: at a paralog module the per-sample residuals are BIMODAL --
                    // most sit near 0 or near 0.5 and almost none in between -- so the mean describes
                    // neither population and moves with their ratio rather than with either one.
                    if (resid > 0.4) ++round_ambiguous;
                    // CN is reported for every traversing haplotype (absolute module count); GT marks a
                    // CARRIER only when the count differs from the reference's (a gain or a loss), so AC/AF
                    // stay meaningful instead of flagging every haplotype.
                    mr.sample_cn[graph.paths[pi].name] = copies;
                    mr.sample_dosage[graph.paths[pi].name] = exact;
                    if (copies != ref_copies) mr.carriers.push_back(graph.paths[pi].name);
                    if (copies > mr.seed.alt_cn) mr.seed.alt_cn = copies;
                }
                mr.seed.round_residual = round_n > 0 ? round_residual / static_cast<double>(round_n) : 0.0;
                mr.seed.round_ambiguous_frac =
                    round_n > 0 ? static_cast<double>(round_ambiguous) / static_cast<double>(round_n) : 0.0;
                mr.seed.cn_clamped_zero = cn_clamped_zero;
                mr.seed.module_span_ambiguous = span_ambiguous;
                // Gate on the SHARE of samples whose rounding was a coin flip, not on the mean residual.
                // The mean is a bad summary of a bimodal distribution -- at GSTM1 the median sample sits
                // at 0.045 and the 90th percentile at 0.496, so the mean of 0.197 describes no sample --
                // and it was also the statistic this record's own header warned against acting on.
                //
                // Opt-in strictness, and still OFF by default, because a bad fit does not mean a wrong
                // answer here: GSTM1 has the worst fit of any reference locus (over half its samples
                // ambiguous) and its CN is exact against pangene truth on all 466. The unit is a
                // calibration constant for a heterogeneous paralog module, not one real copy, so
                // landing near a half-integer is what it does when it is working.
                if (options.max_cn_model_residual > 0.0 &&
                    mr.seed.round_ambiguous_frac > options.max_cn_model_residual) {
                    ++summary.declined_cn_model;
                    // Declining means NO copy-number call here, not a different one. Leaving
                    // leaving the route unconsumed let the peak route answer instead, and it answers something
                    // else entirely -- at GSTM1 the MODULE_BP record covers 309 carriers and the peak
                    // record that replaced it covered 2. A refused measurement must not be silently
                    // substituted by a weaker one; the sequence-resolved events and the allele VCF are
                    // what remain.
                    cn_route_consumed = true;
                    module_cn_declined = true;
                    if (!options.quiet) {
                        std::cerr << "[call] bubble " << bubble.id << ": module CN declined, "
                                  << mr.seed.round_ambiguous_frac * 100.0
                                  << "% of traversers round ambiguously, over "
                                  << options.max_cn_model_residual * 100.0 << "%"
                                  << " (sequence-resolved events and the allele VCF are unaffected)\n";
                    }
                } else
                if (!mr.carriers.empty()) {  // a CN-invariant module is not a variant record
                    // describe handoff: a coverage DUP's CN signal lives in the folded module's inside
                    // nodes (their per-walk multiplicity), not the bubble source. Carry them in
                    // variant_nodes.tsv or `describe --variant-nodes` masks to the source and drops it.
                    mr.member_nodes.insert(bubble.inside.begin(), bubble.inside.end());
                    cn_route_consumed = true;
                    cn_unit_bp = unit;   // one-copy size, used to drop CN-loss DELs below
                    merged.push_back(std::move(mr));
                }
            }
        }

        // ---- DUP/CN events: count self-loop traversals per allele vs reference. Merge
        // on shared REP node; per-sample CN. Independent of the walk-diff alignment.
        // Skipped when bp-coverage fired (it is the authority for that bubble's copy number).
        // Node-outer, allele-inner: every allele's traversal count of a REP node is a copy number for
        // that node, including the alleles that match the reference. Allele-outer skipped those before
        // the record existed, so a reference-like traverser ended up with no CN at all -- which the
        // writer turned into REF_CN by accident rather than by measurement, and which left the per-gene
        // sidecar with nothing to split for hundreds of LPA and ANKRD36C haplotypes.
        for (const std::string& cn : cn_nodes) {
            if (cn_route_consumed) break;
            // Only a genuine REP unit (self-loop >= min_sv_bp) anchors a self-loop DUP; an incidental
            // tiny self-loop would emit a spurious REF_CN=0 DUP, so skip it.
            if (node_len(graph, cn) < options.min_sv_bp) continue;
            const std::size_t rc = ref_count.count(cn) ? ref_count.at(cn) : 0;
            const std::size_t unit = node_len(graph, cn);
            std::vector<std::size_t> ac_of(alleles.size(), 0);
            for (std::size_t ai = 0; ai < alleles.size(); ++ai) {
                std::unordered_map<std::string, std::size_t> alt_count;
                for (const PathStep& s : alleles[ai].steps) ++alt_count[s.node_id];
                ac_of[ai] = alt_count.count(cn) ? alt_count.at(cn) : 0;
            }
            // The record exists if ANY allele differs from the reference by a callable amount; the
            // qualifying test is unchanged, only where it sits.
            MergedRecord* grp = nullptr;
            for (std::size_t ai = 0; ai < alleles.size(); ++ai) {
                const std::size_t ac = ac_of[ai];
                if (rc == ac) continue;
                const std::size_t delta = rc > ac ? rc - ac : ac - rc;
                if (unit * delta < options.min_sv_bp) continue;
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
                    mr.seed.cn_method = CnMethod::Rep;
                    merged.push_back(std::move(mr));
                    grp = &merged.back();
                }
            }
            if (grp == nullptr) continue;
            // Now fill EVERY allele, reference-like ones included. Carriers stay the differing ones, so
            // AC/AF are unchanged; what changes is that a traverser always has a measured CN.
            for (std::size_t ai = 0; ai < alleles.size(); ++ai) {
                const std::size_t ac = ac_of[ai];
                const std::size_t delta = rc > ac ? rc - ac : ac - rc;
                // Carrier keeps exactly the rule that produced this record -- a callable difference --
                // so AC/AF are untouched. Completeness is about sample_cn, not about who counts as a
                // carrier: a sub-threshold difference is a real CN and not a called SV.
                const bool carrier = (ac != rc) && (unit * delta >= options.min_sv_bp);
                for (const std::string& m : alleles[ai].members) {
                    grp->sample_cn[m] = ac;
                    if (carrier) grp->carriers.push_back(m);
                }
            }
        }

        // Peak-multiplicity DUP: where panphorte could not fold the cluster (no self-loop node), CN is
        // the per-haplotype peak node traversal count. Using the peak rather than any node above the
        // reference rejects cluster background: scattered per-node excesses are paralog presence/absence,
        // the peak is dosage. Gated on the absence of a REP self-loop so routes stay disjoint.
        if (options.cn && !has_rep_selfloop && !cn_route_consumed) {
            std::size_t ref_peak = 0;
            for (const std::string& id : bubble.inside) {
                const auto it = ref_count.find(id);
                if (it != ref_count.end() && it->second > ref_peak) ref_peak = it->second;
            }
            MergedRecord* peak_grp = nullptr;
            std::vector<std::size_t> peak_of(alleles.size(), 0);   // every allele's peak, qualifying or not
            std::vector<char> qualified(alleles.size(), 0);         // ...and whether it cleared the rule
            for (std::size_t ai = 0; ai < alleles.size(); ++ai) {
                std::unordered_map<std::string, std::size_t> alt_count;
                for (const PathStep& s : alleles[ai].steps) ++alt_count[s.node_id];
                // ac_peak/peak_node decide IF the allele duplicates; the peak node is often tiny, so
                // event size is the excess bp the allele traverses over the reference across all folded
                // nodes (the actual duplicated content), not the node length.
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
                peak_of[ai] = ac_peak;
                if (ac_peak < 2 || ac_peak <= ref_peak) continue;
                if (excess_bp < options.min_sv_bp) continue;
                qualified[ai] = 1;
                if (peak_grp == nullptr) {
                    MergedRecord mr;
                    mr.seed.type = EvType::Dup;
                    mr.seed.nodes.push_back(peak_node);
                    mr.seed.start_node = peak_node; mr.seed.end_node = peak_node;
                    mr.seed.ref_cn = ref_peak; mr.seed.alt_cn = ac_peak;
                    mr.seed.anchor_node = bubble.source;
                    mr.seed.size_bp = excess_bp;
                    mr.seed.ru_len = excess_bp / (ac_peak - ref_peak);  // per-copy duplicated content
                    mr.seed.cn_method = CnMethod::Peak;
                    // describe handoff: a peak DUP's CN signal lives in the folded module's inside nodes
                    // (their per-walk multiplicity), not the peak node. Carry them in variant_nodes.tsv
                    // or `describe --variant-nodes` masks to one invariant node and drops the feature.
                    mr.member_nodes.insert(bubble.inside.begin(), bubble.inside.end());
                    merged.push_back(std::move(mr));
                    peak_grp = &merged.back();
                } else if (ac_peak > peak_grp->seed.alt_cn) {
                    peak_grp->seed.alt_cn = ac_peak;
                    peak_grp->seed.size_bp = excess_bp;
                    peak_grp->seed.ru_len = excess_bp / (ac_peak - ref_peak);
                }
            }
            // Fill every allele, not only the duplicating ones: a traverser at reference multiplicity
            // has a measured CN of ref_peak, and leaving it absent made the writer invent one.
            if (peak_grp != nullptr) {
                for (std::size_t ai = 0; ai < alleles.size(); ++ai) {
                    for (const std::string& m : alleles[ai].members) {
                        peak_grp->sample_cn[m] = peak_of[ai];
                        // Carrier keeps the original rule, excess_bp test included, so AC/AF do not move.
                        if (qualified[ai]) peak_grp->carriers.push_back(m);
                    }
                }
            }
        }

        // DEL/INS/INV per allele, kept down to the rescue floor, then merged. Emitted even when coverage
        // fired: a coverage DUP is dosage, not the sequence-resolved indels inside the bubble. The
        // per-allele walk-diff dominates on big folded bubbles and writes only its own slot, so it runs
        // across cores; the merge below consumes in allele order, independent of thread count.
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

        // Cross-haplotype merge by transitive single-linkage: connected components over the
        // events_match graph, so A~B~C collapse to one record even when A!~C (greedy first-fit would
        // fragment them). Can only ever merge more, never split or drop a carrier.
        {
            struct Cand { std::size_t ai; const Event* e; };
            std::vector<Cand> cands;
            for (std::size_t ai = 0; ai < alleles.size(); ++ai)
                for (const Event& e : allele_events[ai])
                    if (e.size_bp >= rescue_floor) cands.push_back({ai, &e});

            DisjointSet dsu(cands.size());
            // Position-sorted windowed sweep keeps edge-building near-linear: events_match can only
            // fire within |dpos| <= merge_distance_bp + size_i, so we stop once past that bound.
            // Zero-coordinate events (anchor missing from the reference) match on node/sequence only
            // and are never position-pruned.
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
                // Single-linkage merges A with C whenever some B links to both, so the record can span
                // members that would never have merged directly. MERGE_JACCARD reports the STRONGEST
                // edge, which says nothing about that. The diameter -- the WEAKEST pairwise similarity
                // in the component -- is what shows how far the chain reached. Bounded, because this is
                // quadratic and a component can hold hundreds of members.
                double diameter = -1.0;
                bool diameter_exact = true;
                if (comp.size() >= 2) {
                    constexpr std::size_t kDiameterMaxMembers = 128;
                    if (comp.size() <= kDiameterMaxMembers) {
                        diameter = 1.0;
                        for (std::size_t x = 0; x < comp.size(); ++x)
                            for (std::size_t y = x + 1; y < comp.size(); ++y)
                                diameter = std::min(diameter, weighted_jaccard(
                                    graph, cands[comp[x]].e->nodes, cands[comp[y]].e->nodes));
                    } else {
                        // Against the representative only: an upper bound on the true diameter, and
                        // labelled as such rather than passed off as the real thing.
                        diameter = 1.0;
                        diameter_exact = false;
                        for (std::size_t k : comp)
                            if (k != best)
                                diameter = std::min(diameter, weighted_jaccard(
                                    graph, cands[best].e->nodes, cands[k].e->nodes));
                    }
                }
                MergedRecord mr;
                mr.seed = *cands[best].e;                 // largest member represents the record
                mr.merge_diameter = diameter;
                mr.merge_diameter_exact = diameter_exact;
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

        // Force-call: re-read every non-member haplotype's walk-diff and add it as a carrier if it
        // supports the record's fixed representative, even below the size threshold. Monotone, so one
        // pass suffices. A node-set containment test was rejected: on folded loci the inserted nodes are
        // shared across nearly all walks, so it force-called everyone.
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

        // De-dup a folded duplication against its module DUP: both see the extra copy (the DUP as a
        // count, the walk-diff as a duplicated-content INS). Keyed by the same carriers + comparable
        // bp (the copy may use off-reference paralog nodes). Drop the matching INS; keep genuine novel
        // insertions (more carriers, or a very different size).
        if (options.cn) {
            std::vector<std::pair<std::unordered_set<std::string>, std::size_t>> peak_dups;
            for (const MergedRecord& mr : merged) {
                if (mr.seed.type == EvType::Dup && is_module_cn(mr.seed.cn_method)) {
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

        // Drop copy-number DUPs that are low-complexity-tangle or physically implausible (span > a
        // large fraction of the reference). Runs AFTER the folded-INS de-dup above, so the matching INS
        // stays dropped -- we remove only the bogus DUP record, without resurfacing its content. Genuine
        // self-loop REP DUPs are never touched.
        if (is_tangle || max_dup_bp > 0) {
            std::vector<MergedRecord> kept;
            for (MergedRecord& mr : merged) {
                const bool cn_dup = (mr.seed.type == EvType::Dup && is_module_cn(mr.seed.cn_method));
                if (cn_dup && (is_tangle || (max_dup_bp > 0 && mr.seed.size_bp > max_dup_bp))) {
                    ++summary.oversized_dups;
                    continue;
                }
                kept.push_back(std::move(mr));
            }
            merged = std::move(kept);
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
        // An arm can vanish, or a link_id can smear onto a merged cluster, so strip any link_id lacking
        // both a DEL and an INS among survivors. An arm orphaned that way is dropped rather than
        // surfaced as a bare INS/DEL.
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

        // coverage-CN bubble: drop copy-number-loss DELs (a DEL >= half a copy unit is just fewer folded
        // copies, already reported by the coverage DUP -- emitting both double-counts). Keep novel
        // INS/INV, small local DELs, and substitution arms. Runs after the EVENTID contract pass so an
        // orphaned DEL reads as a lone CN-loss DEL.
        if (module_cn_declined) ++summary.declined_cn_model_bubbles;
        if (cn_route_consumed && cn_unit_bp > 0.0) {
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
        // The allele VCF is the LOSSLESS record of what the graph holds, so it must not be gated on
        // the interpreted caller having found something. It was: a bubble whose events were all
        // filtered -- by --min-sv-bp, by support, or by the tangle/oversize suppression that can remove
        // a MODULE_BP record after its duplicated-content INS was already dropped as redundant with it
        // -- returned here and produced no allele record either, so the one output that could still
        // describe the site described nothing. The interpreted work below stays skipped; only the
        // allele record now survives an empty call set.
        const bool no_interpreted_calls = merged.empty();
        if (!no_interpreted_calls) ++summary.bubbles_with_calls;

        // Map each sample to its allele (for per-carrier sub-walk provenance).
        std::unordered_map<std::string, std::size_t> sample_to_allele;
        for (std::size_t ai = 0; ai < alleles.size(); ++ai) {
            for (const std::string& m : alleles[ai].members) sample_to_allele[m] = ai;
        }

        // Genes overlapping this bubble's reference span (for the per-gene DUP table, resolved by
        // realignment in a post-pass). The GENES INFO field is graph-based (collapsed nodes tag
        // multiple genes); the per-gene COPY NUMBER is resolved later by competitive realignment,
        // which separates collapsed paralogs the graph cannot.
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

        // Lossless companion output (--allele-vcf): one record per bubble carrying EVERY distinct
        // allele as explicit sequence, each haplotype's GT indexing its own allele. The merged records
        // in the region VCF are an interpretation -- several of them can describe one walk, and a
        // carrier is given the merged representative's length rather than its own -- so a consumer
        // reconstructing a specific sample needs this instead. Written alongside, never in place of,
        // the region VCF, and deliberately not gated on CN: the region VCF keeps REF_CN/CN semantics
        // and this file keeps the sequence.
        if (options.allele_vcf) {
            // Anchor orientation-independently, as the merged path does: when the reference crosses
            // this bubble sink->source the genomically-upstream flank is the SINK, and the allele
            // sequences -- spelled canonically source->sink -- have to be reverse-complemented to sit
            // in reference-forward coordinates alongside POS.
            const auto a_ps = ref_node_pos.find(bubble.source);
            const auto a_pk = ref_node_pos.find(bubble.sink);
            const bool a_rev = a_ps != ref_node_pos.end() && a_pk != ref_node_pos.end() &&
                               a_pk->second < a_ps->second;
            const std::string& a_anchor_id = a_rev ? bubble.sink : bubble.source;
            const auto asit = ref_node_pos.find(a_anchor_id);
            const auto asnode = graph.nodes.find(a_anchor_id);
            auto a_interior = [&](const std::vector<PathStep>& steps) -> std::string {
                if (steps.size() < 2) return std::string();
                std::vector<PathStep> inner(steps.begin() + 1, steps.end() - 1);
                std::string seq = spell_path_steps_sequence(graph, inner);
                return a_rev ? reverse_complement(seq) : seq;
            };
            if (asit == ref_node_pos.end() || asnode == graph.nodes.end() || asnode->second.sequence.empty()) {
                ++summary.allele_skipped;
            } else {
                const std::size_t aslen = asnode->second.sequence.size();
                const std::string a_anchor = upper_base(asnode->second.sequence[aslen - 1]);
                const std::size_t a_pos = asit->second + (aslen > 0 ? aslen - 1 : 0);
                const std::string a_ref = a_anchor + a_interior(ref_steps);
                const std::size_t cap = options.allele_vcf_max_bp;
                std::vector<std::string> a_alts;
                std::unordered_map<std::string, int> a_idx_of;
                a_idx_of.emplace(a_ref, 0);
                std::vector<int> a_vcf_idx(alleles.size(), -1);
                bool a_ok = !(cap && a_ref.size() > cap);
                for (std::size_t ai = 0; a_ok && ai < alleles.size(); ++ai) {
                    const std::string seq = a_anchor + a_interior(alleles[ai].steps);
                    if (cap && seq.size() > cap) { a_ok = false; break; }
                    auto it = a_idx_of.find(seq);
                    if (it == a_idx_of.end()) {
                        const int idx = 1 + static_cast<int>(a_alts.size());
                        a_alts.push_back(seq);
                        a_idx_of.emplace(seq, idx);
                        a_vcf_idx[ai] = idx;
                    } else {
                        a_vcf_idx[ai] = it->second;
                    }
                }
                if (!a_ok || a_alts.empty()) {
                    ++summary.allele_skipped;
                } else {
                    std::vector<std::size_t> a_ac(a_alts.size(), 0);
                    std::size_t a_an = 0;
                    for (const std::string& s : sample_names) {
                        if (!traverses.count(s)) continue;
                        ++a_an;
                        const auto ait = sample_to_allele.find(s);
                        if (ait == sample_to_allele.end()) continue;
                        const int vi = a_vcf_idx[ait->second];
                        if (vi >= 1) ++a_ac[static_cast<std::size_t>(vi - 1)];
                    }
                    const std::size_t a_end = a_pos + a_ref.size() - 1;
                    std::ostringstream info;
                    info << "BUBBLE_ID=" << bubble.id << ";END=" << a_end
                         << ";NALLELES=" << (a_alts.size() + 1) << ";AN=" << a_an << ";AC=";
                    for (std::size_t k = 0; k < a_ac.size(); ++k) { if (k) info << ','; info << a_ac[k]; }
                    info << ";SVLEN=";
                    for (std::size_t k = 0; k < a_alts.size(); ++k) {
                        if (k) info << ',';
                        info << (static_cast<long long>(a_alts[k].size()) - static_cast<long long>(a_ref.size()));
                    }
                    std::ostringstream row;
                    row << ref_meta.chrom << '\t' << a_pos << '\t'
                        << ("bubble" + std::to_string(bubble.id) + "_ALLELES") << '\t' << a_ref << '\t';
                    for (std::size_t k = 0; k < a_alts.size(); ++k) { if (k) row << ','; row << a_alts[k]; }
                    row << "\t.\t.\t" << info.str() << "\tGT";
                    for (const std::string& s : sample_names) {
                        row << '\t';
                        if (!traverses.count(s)) { row << '.'; continue; }
                        const auto ait = sample_to_allele.find(s);
                        row << (ait != sample_to_allele.end() && a_vcf_idx[ait->second] >= 0
                                    ? std::to_string(a_vcf_idx[ait->second]) : ".");
                    }
                    row << '\n';
                    OutRecord arec;
                    arec.pos = a_pos; arec.end = a_end; arec.bubble_id = bubble.id;
                    arec.id = "bubble" + std::to_string(bubble.id) + "_ALLELES";
                    arec.line = row.str();
                    allele_records.push_back(std::move(arec));
                    ++summary.allele_records;
                }
            }
        }

        if (no_interpreted_calls) return;   // allele record written above; nothing to interpret

        // Optional multiallelic record (--multiallelic-loci): collapse a bounded locus (STR/VNTR) into
        // one record with explicit-sequence alleles (REF + ALTs), GT indexing each sample's allele.
        // Skipped when an allele exceeds --multiallelic-max-bp, or the bubble carries a CN record (would
        // discard REF_CN/CN) -- so it applies to pure DEL/INS/INV bubbles only. Opt-in.
        // A DECLINED module CN emitted no record, so it must not suppress the multiallelic form --
        // that was the flag meaning two things at once.
        bool bubble_has_cn = false;
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
        for (const MergedRecord& mr : merged) {
            const Event& e = mr.seed;
            // POS / REF base from the reference node map. Anchor orientation-independently: a
            // reverse-oriented bubble (source genomically downstream of sink) is walked in decreasing
            // coordinate, so the walk-order anchor would point the wrong way. Take the genomically
            // upstream reference node of the event; END then extends in increasing coordinate.
            auto rpos = [&](const std::string& n) -> long long {
                const auto it = ref_node_pos.find(n);
                return it == ref_node_pos.end() ? -1 : static_cast<long long>(it->second);
            };
            std::string anchor = e.anchor_node;
            bool anchor_after = e.anchor_after;
            // A symbolic DEL/INV is anchored on the base PRECEDING the event, so REF carries one real
            // base and ALT is the symbol. Computed as a coordinate (first affected base minus one)
            // rather than by stepping back a node: the preceding node is only identifiable by walk
            // position, and every node->position map here holds first occurrences.
            long long del_inv_pos = -1;
            long long del_inv_end = -1;
            long long ins_pos = -1;
            if (e.type == EvType::Del || e.type == EvType::Inv) {
                // Preferred: the walk positions this event actually occupies, which name the
                // OCCURRENCE. The node-name fallback below takes each node's FIRST occurrence, so a
                // reference that revisits an event node anchors the record at the wrong copy -- on a
                // fixture where the reference visits the deleted node twice, 150 bp upstream of the
                // deletion, inside sequence the haplotype still carries.
                const auto span = tok_span(e.ref_tok_first, e.ref_tok_last);
                if (span.first > 0) {
                    anchor_after = false;
                    del_inv_end = static_cast<long long>(span.second);
                    if (span.first > ref_meta.region_start_1based)
                        del_inv_pos = static_cast<long long>(span.first) - 1;
                    // keep `anchor` pointing at a node of the event for START_NODE/REF fallback
                    anchor = e.nodes.empty() ? anchor : e.nodes.front();
                } else {
                    // First (genomically-lowest) reference node of the event.
                    long long best = -1;
                    for (const std::string& n : e.nodes) {
                        const long long p = rpos(n);
                        if (p >= 0 && (best < 0 || p < best)) { best = p; anchor = n; }
                    }
                    anchor_after = false;
                    if (best > static_cast<long long>(ref_meta.region_start_1based)) del_inv_pos = best - 1;
                }
            } else if (e.type == EvType::Ins) {
                // anchor_node is the walk-order flank of the insertion. In a forward bubble it is the
                // genomically-upstream flank (POS = its last base); in a reverse-oriented bubble it is
                // the downstream flank, so POS is its first base (the insertion sits just before it).
                const long long ps = rpos(bubble.source), pk = rpos(bubble.sink);
                const bool reverse_bubble = (ps >= 0 && pk >= 0 && pk < ps);
                anchor_after = !reverse_bubble;
                // Prefer the flank's OWN occurrence over its first one. Without this an insertion
                // after the second visit to a repeated node is placed at the first visit, exactly as
                // DEL was: the node name is the same in both places.
                const auto flank = tok_span(e.ref_tok_anchor, e.ref_tok_anchor);
                if (flank.first > 0) ins_pos = static_cast<long long>(
                    reverse_bubble ? flank.first : flank.second);
            } else {
                // DUP: anchor on the genomically upstream flank, POS = its LAST base (as INS does), since
                // the duplicated content starts where that flank ends. Using the first base only looks
                // right when flanks are short; against a multi-kb flank node POS lands kb upstream.
                const long long ps = rpos(bubble.source), pk = rpos(bubble.sink);
                if (pk >= 0 && (ps < 0 || pk < ps)) anchor = bubble.sink;
                else anchor = bubble.source;
                anchor_after = true;
            }
            std::size_t pos = 1;
            std::string ref_base = "N";
            const auto ait = ref_node_pos.find(anchor);
            if (ait != ref_node_pos.end()) {
                const std::size_t glen = node_len(graph, anchor);
                if (anchor_after && glen > 0) {
                    pos = ait->second + glen - 1; // last base of the upstream flank (INS)
                } else {
                    pos = ait->second;            // first base of the first event node
                }
                const auto nit = graph.nodes.find(anchor);
                if (nit != graph.nodes.end() && !nit->second.sequence.empty()) {
                    const std::size_t bi = (anchor_after && glen > 0) ? glen - 1 : 0;
                    ref_base = upper_base(nit->second.sequence[bi]);
                }
            }
            // Step POS back onto the preceding base. Done after the node lookup above so the REF base
            // comes from the coordinate rather than from the event's own first node. An event starting
            // at the very first reference base has no preceding base and keeps its original anchor.
            if (del_inv_pos > 0) {
                const std::string b = ref_base_at(static_cast<std::size_t>(del_inv_pos));
                if (!b.empty()) {
                    pos = static_cast<std::size_t>(del_inv_pos);
                    ref_base = b;
                }
            }
            if (ins_pos > 0) {
                const std::string b = ref_base_at(static_cast<std::size_t>(ins_pos));
                if (!b.empty()) {
                    pos = static_cast<std::size_t>(ins_pos);
                    ref_base = b;
                }
            }
            std::size_t end = pos;
            long long svlen = 0;
            const char* svt = ev_svtype(e.type);
            if (e.type == EvType::Del) { svlen = -static_cast<long long>(e.size_bp); end = pos + e.size_bp; }
            else if (e.type == EvType::Ins) { svlen = static_cast<long long>(e.size_bp); end = pos; }
            else if (e.type == EvType::Inv) { svlen = static_cast<long long>(e.size_bp); end = pos + e.size_bp; }
            // Both module routes carry their size in size_bp: MODULE_BP the calibrated one-copy unit,
            // PEAK the duplicated content. Only the REP route can derive it from a node length, because
            // only there is the node the literal repeat unit -- at a module the anchor is the bubble
            // source, which is whatever short node happens to bound the site. The flag these two shared
            // before the split was load-bearing exactly here.
            else if (is_module_cn(e.cn_method)) {
                svlen = static_cast<long long>(e.size_bp);
                end = pos + e.size_bp;
            }
            else { // self-loop DUP: copies x unit length (signed: CN gain or loss vs reference)
                svlen = (static_cast<long long>(e.alt_cn) - static_cast<long long>(e.ref_cn)) *
                        static_cast<long long>(node_len(graph, e.nodes.empty() ? std::string() : e.nodes.front()));
                end = pos + node_len(graph, e.nodes.empty() ? std::string() : e.nodes.front());
            }
            // A coalesced record can describe a deletion in SEVERAL pieces with reference retained
            // between them. That is not one exact deletion: the affected interval is wider than the
            // bases removed. Such a record is marked IMPRECISE and END describes the whole affected
            // span, while SVLEN stays the net deleted bases. END-POS == |SVLEN| therefore holds for
            // every PRECISE record and the flag announces the records where it cannot.
            bool imprecise = false;
            if (del_inv_end > 0 && (e.type == EvType::Del || e.type == EvType::Inv)) {
                const long long span_bp = del_inv_end - static_cast<long long>(pos);
                const long long mag = svlen < 0 ? -svlen : svlen;
                if (span_bp > mag) {
                    end = static_cast<std::size_t>(del_inv_end);
                    imprecise = true;
                }
            }
            // A COLLAPSED_MODULE record describes the whole site, and its "unit" is only the shared part
            // of one copy, so pos + unit named a reference interval that corresponds to nothing. END is
            // the module's own reference span -- the interval FORMAT:CNBP is measured over.
            if (is_module_cn(e.cn_method)) {
                const auto sit = ref_node_pos.find(bubble.sink);
                const auto bit = ref_node_pos.find(bubble.source);
                if (sit != ref_node_pos.end() && bit != ref_node_pos.end()) {
                    // build_ref_node_pos already returns 1-based genomic starts, so the span end is the
                    // far boundary's own start plus its length -- adding the region offset again put END
                    // past the locus and it silently clamped to the last base.
                    const bool sink_is_far = sit->second >= bit->second;
                    const std::size_t far_start = sink_is_far ? sit->second : bit->second;
                    const std::string& far = sink_is_far ? bubble.sink : bubble.source;
                    // END closes the interval FORMAT:CNBP is measured over, and that sum runs over the
                    // bubble's INTERIOR only -- both boundary nodes are excluded. POS is already the
                    // last base of the near boundary, so the interior is POS+1 .. far_start-1 and END
                    // is the base before the far boundary starts. Including the far boundary node (or,
                    // as before, the far boundary plus one) described an interval no other field means.
                    (void)far;
                    end = far_start > 0 ? far_start - 1 : 0;
                }
            }
            if (end > ref_end_1based) end = ref_end_1based;  // END is a reference coordinate; never past the graph

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

            // variant_nodes.tsv must cover every node any carrier traverses, so take the union of the
            // merged members' nodes and the representative's. This keeps VCF EVENT_NODES a subset of
            // variant_nodes; otherwise a coverage DUP anchored on a bubble boundary loses its POS node.
            std::vector<std::string> var_nodes;
            if (mr.member_nodes.empty()) {
                var_nodes = ev_nodes;
            } else {
                std::set<std::string> all(mr.member_nodes.begin(), mr.member_nodes.end());
                all.insert(ev_nodes.begin(), ev_nodes.end());
                std::vector<std::pair<long long, std::string>> keyed;
                long long tail = 1LL << 60;
                for (const std::string& nd : all) {
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
            info << "END=" << end << ";SVTYPE=" << svt;
            // At a collapsed module the carriers both gain and lose -- GSTM1 spans CN 1..4 around
            // REF_CN=3 -- so a single record-level size is not merely imprecise, it does not exist.
            // Reporting the shared unit there understated every carrier by the shared-to-total ratio.
            if (!is_module_cn(e.cn_method)) info << ";SVLEN=" << svlen;
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
            if (mr.merge_diameter >= 0.0) {
                std::ostringstream d; d.setf(std::ios::fixed); d.precision(4);
                d << mr.merge_diameter;
                info << ";MERGE_DIAMETER=" << d.str();
                if (mr.merge_diameter_exact) info << ";MERGE_DIAMETER_EXACT";
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
            if (e.type == EvType::Dup && e.ru_len > 0) {
                if (e.cn_method == CnMethod::Rep) info << ";RU_LEN=" << e.ru_len;
                else if (e.cn_method == CnMethod::ModuleBp) info << ";CN_UNIT_BP=" << e.ru_len;
                // PEAK claims no unit at all: its "unit" is excess content divided by a copy delta,
                // which is an average over whatever the duplication happened to include.
            }
            // Which measurement produced CN, and what a copy is a copy OF. Without these a consumer
            // cannot tell an exact repeat-unit traversal count from a paralog module's dosage, and
            // reads RU_LEN as a per-copy size in both cases -- which it only is for the first.
            if (e.type == EvType::Dup && e.cn_method != CnMethod::None) {
                info << ";CN_METHOD=" << cn_method_name(e.cn_method)
                     << ";CN_SCOPE=" << cn_scope_name(e.cn_method);
            }
            if (e.shared_fold_bp > 0) info << ";CN_SHARED_BP=" << e.shared_fold_bp;
            if (e.ref_fold > 0) info << ";CN_REF_FOLD=" << e.ref_fold;
            if (e.fold_residual >= 0.0) {
                std::ostringstream r; r.setf(std::ios::fixed); r.precision(4);
                r << e.fold_residual;
                info << ";CN_REF_MULTIPLICITY_HETEROGENEITY=" << r.str();
            }
            if (e.max_support >= 0.0) {
                std::ostringstream r; r.setf(std::ios::fixed); r.precision(4);
                r << e.max_support; info << ";CN_REF_MAX_SUPPORT=" << r.str();
            }
            if (e.cn_method == CnMethod::ModuleBp) {
                info << ";CN_DOSAGE_MODEL="
                     << (e.dosage_spacing ? "PANEL_SPACING" : "REFERENCE_RATIO");
            }
            if (e.step_bp > 0.0) {
                info << ";CN_STEP_BP=" << static_cast<long long>(e.step_bp + 0.5);
                std::ostringstream r; r.setf(std::ios::fixed); r.precision(3);
                r << (e.ru_len > 0 ? e.step_bp / static_cast<double>(e.ru_len) : 0.0);
                info << ";CN_STEP_RATIO=" << r.str();
                info << ";CN_STEP_SUPPORT=" << e.step_clusters << ',' << e.step_gaps << ','
                     << e.step_dropped;
                if (e.step_max_offint >= 0.0) {
                    std::ostringstream o; o.setf(std::ios::fixed); o.precision(3);
                    o << e.step_max_offint; info << ";CN_STEP_OFFINT=" << o.str();
                    std::ostringstream m; m.setf(std::ios::fixed); m.precision(2);
                    m << e.step_max_multiple; info << ";CN_STEP_MAX_MULTIPLE=" << m.str();
                }
            }
            // REF_CN is not measured the same way on every route, and absolute CN is only as good as
            // its anchor. Say which one produced it rather than leaving a bare integer.
            if (e.type == EvType::Dup && e.cn_method != CnMethod::None) {
                info << ";REF_CN_SOURCE="
                     << (e.cn_method == CnMethod::Rep ? "REP_TRAVERSAL" : "MAX_NODE_MULTIPLICITY");
            }
            if (e.cn_method == CnMethod::Peak) info << ";CN_CONFIDENCE=HEURISTIC";
            if (imprecise) info << ";IMPRECISE";
            if (e.module_span_ambiguous > 0)
                info << ";CN_SPAN_AMBIGUOUS=" << e.module_span_ambiguous;
            if (e.cn_clamped_zero > 0) info << ";CN_CLAMPED_ZERO=" << e.cn_clamped_zero;
            if (e.round_ambiguous_frac >= 0.0) {
                std::ostringstream r; r.setf(std::ios::fixed); r.precision(4);
                r << e.round_ambiguous_frac; info << ";CN_ROUND_AMBIGUOUS_FRAC=" << r.str();
            }
            if (e.round_residual >= 0.0) {
                std::ostringstream r; r.setf(std::ios::fixed); r.precision(4);
                r << e.round_residual; info << ";CN_ROUND_RESIDUAL=" << r.str();
            }
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

            // FORMAT:CNBP (DUP only): the linear bp a haplotype gains/loses through the module -- spelled
            // walk length (node length x traversal multiplicity) minus the reference's. Unlike RU_LEN
            // (one folded copy) this is the real SV size, since the walk still traverses shared nodes its
            // own number of times. Uses the widest source..sink span, not the repeat-collapsed one.
            const bool is_dup = (e.type == EvType::Dup);
            const std::unordered_set<std::string> cnbp_inside(bubble.inside.begin(), bubble.inside.end());
            auto cnbp_walk_bp = [&](std::size_t pi) -> long long {
                const auto sp = module_span(pi);
                if (sp.first > sp.second) return 0;
                const std::vector<PathStep>& steps = graph.paths[pi].steps;
                long long bp = 0;
                for (std::size_t i = sp.first; i <= sp.second && i < steps.size(); ++i)
                    if (cnbp_inside.count(steps[i].node_id)) bp += static_cast<long long>(node_len(graph, steps[i].node_id));
                return bp;
            };
            const long long cnbp_ref_bp = (is_dup && ref_idx < graph.paths.size()) ? cnbp_walk_bp(ref_idx) : 0;
            // The reference's TOTAL bp across this bubble, beside CN_SHARED_BP (its folded subset). At a
            // collapsed paralog module the two differ by the paralog-private content that travels with
            // each copy but that the reference does not revisit -- which is exactly why RU_LEN, derived
            // from the shared part, understates what a carrier gains or loses.
            if (is_dup && cnbp_ref_bp > 0) info << ";CN_MODULE_REF_BP=" << cnbp_ref_bp;
            auto sample_cnbp = [&](const std::string& s) -> long long {
                const auto pit = name_to_pi.find(s);
                return pit == name_to_pi.end() ? 0 : cnbp_walk_bp(pit->second) - cnbp_ref_bp;
            };

            std::ostringstream row;
            row << ref_meta.chrom << '\t' << pos << '\t' << id << '\t' << ref_base
                << "\t<" << svt << ">\t.\t.\t" << info.str()
                << (is_dup ? (e.cn_method == CnMethod::Rep && e.ru_len > 0
                                  ? "\tGT:CN:CNBP:CNRESID"
                                  : mr.sample_dosage.empty() ? "\tGT:CN:CNBP" : "\tGT:CN:CNBP:CNR_RAW:CNR_MARGIN")
                           : "\tGT:CN");
            for (const std::string& s : sample_names) {
                row << '\t';
                if (!traverses.count(s)) {
                    row << (is_dup ? (mr.sample_dosage.empty() ? ".:.:." : ".:.:.:.:.") : ".:.");
                    continue;
                }
                const bool carrier = carrier_set.count(s) != 0;
                row << (carrier ? "1" : "0") << ':';
                if (is_dup) {
                    // A traverser without a measured CN is a bug in whichever route produced this
                    // record, not a reference-like sample. Filling the gap with REF_CN made a complete
                    // module loss indistinguishable from carrying the reference count, and hid any
                    // failure to compute a CN at all.
                    const auto cit = mr.sample_cn.find(s);
                    if (cit == mr.sample_cn.end()) {
                        throw std::runtime_error(
                            "panvar call: no copy number for " + s + " at " + id +
                            ", which traverses the site (" + cn_method_name(e.cn_method) +
                            " route); a traversing haplotype must always have a measured CN");
                    }
                    row << cit->second;
                    const long long this_cnbp = sample_cnbp(s);
                    row << ':' << this_cnbp;
                    // CNRESID (REP only): what the repeat-copy change does NOT account for. A literal
                    // REP unit has a real per-copy length, so (CN-REF_CN)*RU_LEN is a genuine
                    // prediction of this haplotype's bp change and the remainder is other sequence
                    // change in the same bubble. Deliberately absent on a collapsed module, where
                    // RU_LEN is a calibration constant rather than one copy and the difference would
                    // measure the calibration rather than the biology.
                    if (e.cn_method == CnMethod::Rep && e.ru_len > 0) {
                        const long long predicted =
                            (static_cast<long long>(cit->second) - static_cast<long long>(e.ref_cn)) *
                            static_cast<long long>(e.ru_len);
                        row << ':' << (this_cnbp - predicted);
                    }
                    // Per sample, because that is where the ambiguity lives: how far this haplotype sat
                    // from a whole number of units, and how much room it had before the rounding would
                    // have gone the other way.
                    if (!mr.sample_dosage.empty()) {
                        const auto dit = mr.sample_dosage.find(s);
                        if (dit == mr.sample_dosage.end()) { row << ":.:."; }
                        else {
                            const double resid = std::fabs(dit->second - std::round(dit->second));
                            std::ostringstream d; d.setf(std::ios::fixed); d.precision(3);
                            d << dit->second; row << ':' << d.str();
                            std::ostringstream m; m.setf(std::ios::fixed); m.precision(3);
                            m << (0.5 - resid); row << ':' << m.str();
                        }
                    }
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
            out_records.push_back(std::move(rec));

            if (options.write_variant_nodes) {
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
        summary.tangle_bubbles += bo.sum.tangle_bubbles;
        summary.oversized_dups += bo.sum.oversized_dups;
        summary.skipped_large_segments = g_skipped_segments.load();
        summary.allele_records += bo.sum.allele_records;
        summary.allele_skipped += bo.sum.allele_skipped;
        for (OutRecord& r : bo.records) out_records.push_back(std::move(r));
        for (OutRecord& r : bo.allele_records) allele_records.push_back(std::move(r));
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

    // Every destination is resolved and checked BEFORE anything is opened, then written to a staged
    // sibling and renamed in only once the whole run has succeeded. Previously a failure part-way
    // through left a complete-looking region VCF beside a non-zero exit, and a run could overwrite its
    // own input because nothing compared the two.
    cli::StagedOutputs staged("call");
    std::vector<std::string> finals;
    finals.push_back(options.out_prefix + ".region.vcf");
    if (options.allele_vcf) finals.push_back(options.out_prefix + ".alleles.vcf");
    if (options.write_variant_nodes) finals.push_back(options.out_prefix + ".variant_nodes.tsv");
    std::set<std::size_t> bubble_ids_out;
    if (options.write_per_bubble_vcf) {
        for (const OutRecord& rec : out_records) bubble_ids_out.insert(rec.bubble_id);
        for (const std::size_t bid : bubble_ids_out)
            finals.push_back(options.out_prefix + ".bubble_" + std::to_string(bid) + ".vcf");
    }
    // The GTF sidecars belong to the same transaction as the VCFs. Committing the VCF family first and
    // then writing them directly meant a sidecar failure left the VCFs already installed, and it kept
    // them out of the collision check below -- so an input GTF named <prefix>.node_genes.tsv could be
    // overwritten by the run reading it.
    const bool gtf_active = !options.gtf_path.empty() && !genes.empty();
    if (gtf_active) {
        finals.push_back(options.out_prefix + ".node_genes.tsv");
        finals.push_back(options.out_prefix + ".dup_gene_cn.tsv");
    }
    {
        const std::vector<std::string> inputs = {options.gfa_path, options.bubbles_csv_in,
                                                 options.gtf_path};
        std::unordered_set<std::string> seen_out;
        for (const std::string& f : finals) {
            if (!seen_out.insert(f).second)
                throw std::runtime_error("call: two outputs would be written to the same file: " + f);
            for (const std::string& in : inputs) {
                if (in.empty()) continue;
                std::error_code ec;
                if (f == in || std::filesystem::equivalent(f, in, ec))
                    throw std::runtime_error("call: output '" + f + "' is also an input; refusing to "
                                             "overwrite the data being read");
            }
        }
    }
    // Per-bubble VCFs an earlier run left at this prefix and this run will not rewrite. Collected now
    // but REMOVED only after the transaction commits: deleting first meant a later failure destroyed
    // part of a previous successful result and put nothing in its place. Not gated on
    // --write-per-bubble-vcf either -- a rerun with --no-per-bubble-vcf leaves the old files looking
    // like current output, which is the same defect in a different shape.
    std::vector<std::filesystem::path> obsolete_bubble_vcfs;
    {
        std::error_code ec;
        const std::filesystem::path pfx(options.out_prefix);
        const std::filesystem::path dir = pfx.has_parent_path() ? pfx.parent_path() : std::filesystem::path(".");
        const std::string base = pfx.filename().string() + ".bubble_";
        for (const auto& de : std::filesystem::directory_iterator(dir, ec)) {
            if (ec) break;
            const std::string name = de.path().filename().string();
            if (name.rfind(base, 0) != 0 || de.path().extension() != ".vcf") continue;
            const std::string mid = name.substr(base.size(), name.size() - base.size() - 4);
            if (mid.empty() || mid.find_first_not_of("0123456789") != std::string::npos) continue;
            if (options.write_per_bubble_vcf &&
                bubble_ids_out.count(static_cast<std::size_t>(std::stoull(mid)))) continue;
            obsolete_bubble_vcfs.push_back(de.path());
        }
    }

    const std::string region_path = staged.stage(options.out_prefix + ".region.vcf");
    std::ofstream region_out(region_path);
    if (!region_out) {
        throw std::runtime_error("Failed to write region VCF: " + options.out_prefix + ".region.vcf");
    }
    write_vcf_header(region_out);

    if (options.allele_vcf) {
        std::stable_sort(allele_records.begin(), allele_records.end(),
                         [](const OutRecord& a, const OutRecord& b) {
                             if (a.pos != b.pos) return a.pos < b.pos;
                             if (a.end != b.end) return a.end < b.end;
                             return a.id < b.id;
                         });
        const std::string path = options.out_prefix + ".alleles.vcf";
        std::ofstream allele_out(staged.stage(path));
        if (!allele_out) throw std::runtime_error("Failed to write allele VCF: " + path);
        write_vcf_header(allele_out);
        for (const OutRecord& rec : allele_records) allele_out << rec.line;
        allele_out.close();
        if (!allele_out) throw std::runtime_error("Failed to finalize allele VCF: " + path);
    }

    std::map<std::size_t, std::ofstream> bubble_files;
    for (const OutRecord& rec : out_records) {
        region_out << rec.line;
        if (options.write_per_bubble_vcf) {
            auto& f = bubble_files[rec.bubble_id];
            if (!f.is_open()) {
                const std::string path = options.out_prefix + ".bubble_" + std::to_string(rec.bubble_id) + ".vcf";
                f.open(staged.stage(path));
                if (!f) throw std::runtime_error("Failed to write bubble VCF: " + path);
                write_vcf_header(f);
            }
            f << rec.line;
        }
    }
    region_out.close();
    if (!region_out) {
        throw std::runtime_error("Failed to finalize region VCF: " + options.out_prefix + ".region.vcf");
    }
    for (auto& [bid, f] : bubble_files) {
        f.close();
        if (!f) throw std::runtime_error("Failed to finalize bubble VCF for bubble " + std::to_string(bid));
    }

    if (options.write_variant_nodes) {
        // variant_nodes.tsv: per-variant node set, the bridge for `describe --variant-nodes` (restrict
        // k-mer features to nodes participating in called variation) and the `benchmark` round-trip
        // (which nodes a call explains). Ordered START->END; the caller does not emit per-haplotype walks
        // -- `benchmark` reconstructs those from the graph + this node set.
        const std::string vn_path = options.out_prefix + ".variant_nodes.tsv";
        std::ofstream vn_out(staged.stage(vn_path));
        if (!vn_out) {
            throw std::runtime_error("Failed to write variant nodes TSV: " + vn_path);
        }
        vn_out << "variant_id\tbubble_id\tsvtype\tnode_ids\n";
        for (const std::string& r : variant_nodes_rows) vn_out << r << '\n';
        vn_out.close();
        if (!vn_out) throw std::runtime_error("Failed to finalize variant nodes TSV: " + vn_path);
    }

    // GTF sidecars: node->genes map and per-gene DUP CN. CN comes from private-k-mer dosage -- each
    // gene's discriminative reference (merged CDS, else gene span) yields k-mers unique against its
    // paralogs, and a haplotype's per-copy count of those is its CN. Separates collapsed paralogs that
    // graph multiplicity cannot, without per-haplotype alignment.
    if (!genes.empty()) {
        write_node_genes_tsv(staged.stage(options.out_prefix + ".node_genes.tsv"), node_genes, genes);

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

        const std::string dg_path = options.out_prefix + ".dup_gene_cn.tsv";
        std::ofstream dg_out(staged.stage(dg_path));
        if (!dg_out) {
            throw std::runtime_error("Failed to write per-gene DUP CN TSV: " + dg_path);
        }
        dg_out << "bubble_id\tvariant_id\tsample\tgenes\tcn\treliable\tdosage\thits\tpriv_kmers\n";
        for (const std::string& r : dup_gene_cn_rows) dg_out << r << '\n';
        dg_out.close();
        if (!dg_out) throw std::runtime_error("Failed to finalize per-gene DUP CN TSV: " + dg_path);
    }

    // One commit for the whole family -- VCFs, variant nodes and the GTF sidecars together.
    staged.commit();
    // Only now are the superseded files safe to drop: until the new outputs are installed, removing
    // them would leave a failed run having destroyed the previous good result.
    if (!obsolete_bubble_vcfs.empty()) {
        std::size_t removed = 0;
        for (const std::filesystem::path& p : obsolete_bubble_vcfs) {
            std::error_code ec;
            if (std::filesystem::remove(p, ec) && !ec) ++removed;
            else if (ec) {
                std::cerr << "warning: could not remove stale per-bubble VCF " << p.string()
                          << ": " << ec.message() << " (it does not belong to this run)\n";
            }
        }
        if (removed > 0 && !options.quiet) {
            std::cerr << "[call] removed " << removed
                      << " stale per-bubble VCF(s) from an earlier run at this prefix\n";
        }
    }

    if (summary_out) *summary_out = summary;
    // The command layer reports this summary via RunLog; no duplicate print here.
}

} // namespace panvar
