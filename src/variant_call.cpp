#include "panvar/variant_call.hpp"
#include "panvar/path_inversion.hpp"

#include <ksw2.h>

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
#include <cerrno>
#include <charconv>
#include <cstddef>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <mutex>
#include <numeric>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <cstdio>

namespace panvar {
namespace {

// Skip a walk-diff segment if |ref|*|hap| exceeds this. One byte per cell (the traceback direction), so
// the cap is 64 MB per aligning thread.
constexpr std::size_t kAlignCellCap = 64000000;
// A banded alignment of two allele sequences above this many cells goes to minimap2 instead (one byte of
// traceback per cell, so 32 MB per aligning thread).
constexpr std::size_t kSeqAlignCellCap = 32000000;
// A node-compared substitution is sized by aligning its two runs up to this larger side; above it, by the
// larger side.
constexpr std::size_t kSubEditMaxBp = 20000;

// One change between two alleles found by aligning them: [t0, t1) of the reference allele replaced by
// [q0, q1) of the haplotype's. An empty side makes it a deletion or an insertion.
struct AlignedChange {
    std::size_t t0 = 0, t1 = 0, q0 = 0, q1 = 0;
};

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

enum class EvType { Del, Ins, Inv, Dup, Sub };

const char* ev_svtype(EvType t) {
    switch (t) {
        case EvType::Del: return "DEL";
        case EvType::Ins: return "INS";
        case EvType::Inv: return "INV";
        case EvType::Dup: return "CNV";
        case EvType::Sub: return "SUB";
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
    // occurrence the event sits at, which a node id alone cannot identify.
    std::size_t src_idx = 0;
};

// Token walk with the given copy-number nodes (unit self-loops) dropped entirely -- their copies are
// counted as DUP/CN events instead, so a haplotype's extra copies are never mistyped as INS. Used where
// every loop node must be out of the walk (placing a UNIT record with no reference copy).
std::vector<Tok> collapse_walk(
    const std::vector<PathStep>& steps,
    const std::unordered_set<std::string>& cn_nodes) {

    std::vector<Tok> out;
    for (std::size_t i = 0; i < steps.size(); ++i) {
        const PathStep& s = steps[i];
        if (cn_nodes.count(s.node_id) != 0) {
            continue;   // dropped, so the token index is not the step index: carry it explicitly
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

// Token walk for the DEL/INS/INV alignment. A step on a twin is compared as the loop node it folds into
// (its token is that node's), while the token keeps the step's own node id, so events still name the
// nodes walked; a path routed through the twin is then not a deletion of one node plus an insertion of
// the other. A step is dropped when the node it counts as is in `drop`.
std::vector<Tok> collapse_walk(
    const std::vector<PathStep>& steps,
    const std::unordered_set<std::string>& drop,
    const std::unordered_map<std::string, std::string>& unit_of) {

    std::vector<Tok> out;
    for (std::size_t i = 0; i < steps.size(); ++i) {
        const PathStep& s = steps[i];
        const auto u = unit_of.find(s.node_id);
        const std::string& unit = u == unit_of.end() ? s.node_id : u->second;
        if (drop.count(unit) != 0) continue;
        Tok t;
        PathStep as_unit = s;
        as_unit.node_id = unit;
        t.token = hash_step_token(as_unit);
        t.node_id = s.node_id;
        t.reverse = s.reverse;
        t.src_idx = i;
        out.push_back(std::move(t));
    }
    return out;
}

// How a DUP record's CN was obtained. Two measurements, and CN, RU_LEN and SVLEN cannot be interpreted
// without knowing which:
//   Unit      a panphorte unit graph -- an exact count of a repeat unit's copies, where they pass its anchor
//   Module    (--cn) a collapsed module's copies, where they pass its anchor; the module may hold
//             several paralogs, so its copies need not be one uniform sequence
enum class CnMethod { None, Unit, Module };

const char* cn_method_name(CnMethod m) {
    switch (m) {
        case CnMethod::Unit: return "UNIT";
        case CnMethod::Module: return "MODULE";
        default: return "";
    }
}

// A module DUP rather than a UNIT DUP: it counts copies of a collapsed module, not of a literal repeat unit.
bool is_module_cn(CnMethod m) { return m == CnMethod::Module; }

// A typed event derived from one haplotype walk vs the reference walk.
struct Event {
    EvType type = EvType::Del;
    std::vector<std::string> nodes;   // variant node set (length-weighted Jaccard key)
    std::string start_node;           // graph anchors (first/last involved node)
    std::string end_node;
    std::string seq;                  // event sequence (DELSEQ/INSSEQ/INVSEQ)
    std::size_t ref_cn = 0;           // DUP only
    std::size_t alt_cn = 0;           // DUP only
    bool inverted_dup = false;        // UNIT DUP: some extra copy is reverse-complemented
    // Sub only: reference sequence replaced by other sequence in one gap between anchors. `seq` is the
    // replaced reference, `alt_seq` what replaces it; size_bp is the larger of the two, the size the
    // --min-sv-bp filter judges a substitution by.
    std::string alt_seq;
    std::size_t del_bp = 0;
    std::size_t ins_bp = 0;
    // An event found by aligning the haplotype's allele to the reference's at base level (the default
    // where the bubble carries no copy-number record). It has no reference tokens: seq_before is the
    // genomic position of the base before it, and seq/alt_seq are in walk orientation as elsewhere.
    bool seq_event = false;
    std::size_t seq_before = 0;
    std::string seq_ins;              // seq_event: the inserted sequence in genomic orientation
    bool node_found = false;          // seq_event: a node replacement that is a plain insertion or deletion
    // A node event that is one exact change -- one stretch of reference deleted or replaced, or one
    // insertion point -- written the way an aligned event is (base before it, bases deleted, sequence
    // inserted in genomic orientation), so the two can be compared as the haplotypes they make.
    bool splice_ok = false;
    std::size_t splice_before = 0, splice_del = 0;
    std::string splice_ins;
    // Reference anchoring for VCF coordinates.
    // Positions in the bubble's reference walk that this event occupies (SIZE_MAX = none). These
    // identify the occurrence, which anchor_node alone cannot: a reference that revisits a node gives
    // the same name to two different places, and every node->position map records only the first.
    std::size_t ref_tok_first = SIZE_MAX;   // first affected reference token (DEL/INV)
    std::size_t ref_tok_last = SIZE_MAX;    // last affected reference token (DEL/INV)
    std::size_t ref_tok_anchor = SIZE_MAX;  // last matched reference token before the event (INS);
                                            // UNIT DUP with no reference copy: the step its copies follow
    std::string anchor_node;          // ref node POS is taken from (the walk-order flank)
    bool anchor_after = false;        // true: POS = last base of anchor (INS); false: first base
    // One of the pieces a gap block was split into where counted copies had left the walks between
    // them (see diff_segment): pieces of different places, so never coalesced back into one event.
    bool across_copies = false;
    // A deletion or insertion split out of a replacement at a collapsed module (split_copy_replacements): its
    // place on the reference is its own tokens, not the stretch between the matches around the replacement.
    bool split_piece = false;
    // At a collapsed module: how many of the haplotype's copies (split_copies) lie before this event on its walks.
    // Events with copies between them are not coalesced: the copies are reported apart, by the module record.
    std::size_t copies_before = 0;
    // A node event's place in the two compared walks, half-open token ranges (SIZE_MAX = unknown). An empty
    // side sits at the token after the last match before the event.
    std::size_t r_lo = SIZE_MAX, r_hi = SIZE_MAX, h_lo = SIZE_MAX, h_hi = SIZE_MAX;
    std::size_t size_bp = 0;          // |event| for min_sv filtering / SVLEN magnitude
    std::size_t ref_pos = 0;          // reference genomic position of the anchor (merge window)
    // Which CN route produced this record. They answer different questions: one counts copies of a
    // literal repeat unit, the other copies of a collapsed paralog module whose per-copy content is not
    // one uniform sequence.
    CnMethod cn_method = CnMethod::None;
    std::size_t ru_len = 0;           // DUP only: repeat-unit length in bp (RU_LEN; one copy)
    // Instrumentation for the module route, so what it counted is inspectable. Zero when not set.
    std::size_t shared_fold_bp = 0;   // reference bp over the module's nodes (each visit)
    std::size_t ref_fold = 0;         // the reference's copies of the module
    std::size_t module_ref_bp = 0;    // reference bp across the whole bubble interior
    std::size_t anchor_bp = 0;        // MODULE: bp of the anchor nodes
    std::size_t module_nodes_bp = 0;  // MODULE: bp of the module's nodes, one visit each
    std::size_t module_span_ambiguous = 0;  // traversers whose module span had a repeated boundary
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

// The size of a node-compared substitution as the alignment step would size it: the two runs aligned,
// the largest single change, a run replaced counting its larger side. Sized by its larger side alone,
// one base changed on a long node would be a long event; by edit distance, an unrelated replacement would
// be about half its length, as unrelated sequence matches here and there by chance.
std::size_t substitution_change_bp(const std::string& r, const std::string& h) {
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

// Align one sub-range R[r0,r1) vs H[h0,h1) with the node-token DP (diagonal only on equal tokens; else
// gaps) and append the DEL/INS/INV events to `events`. `preceding_ref_node` anchors an INS that opens
// the segment (the last matched ref node before it). Bounded by kAlignCellCap; segments between shared
// anchors are small. A segment over the cap yields no events and adds one to *skipped, so the caller
// can tell that this haplotype's walk was not fully compared.
void diff_segment(
    const Graph& graph,
    const std::vector<Tok>& R, std::size_t r0, std::size_t r1,
    const std::vector<Tok>& H, std::size_t h0, std::size_t h1,
    const std::string& preceding_ref_node,
    std::size_t preceding_ref_tok,
    std::vector<Event>& events,
    std::vector<std::size_t>* hap_match,
    std::size_t* skipped,
    bool split_at_removed) {

    const std::size_t m = r1 - r0;
    const std::size_t n = h1 - h0;
    if (m == 0 && n == 0) return;
    if (m * n > kAlignCellCap) {
        if (std::getenv("PANVAR_CALL_DEBUG")) {
            std::cerr << "[diff] SKIP segment m=" << m << " n=" << n << " (cap " << kAlignCellCap << ")\n";
        }
        if (skipped != nullptr) ++*skipped;
        return; // unanchored divergent block too large to align; skip
    }

    // DP: match (+2) only when tokens equal, gaps (-1). Mismatch diagonals forbidden
    // so substitutions surface as separate ref-only / hap-only blocks (for INV check).
    const int kMatch = 2;
    const int kGap = -1;
    // Forbidden-mismatch sentinel. It stays below every reachable gap score, or the traceback could take
    // it: INT_MIN/4 is below any score a capped segment can reach, with no overflow risk.
    const int kNeg = std::numeric_limits<int>::min() / 4;
    // Scores in two rolling rows; only the traceback step is kept per cell, one byte, which is what
    // allows a cell cap this large. The step prefers diagonal, else up, else left.
    enum : std::uint8_t { kLeft = 0, kUp = 1, kDiag = 2 };
    const std::size_t W = n + 1;
    std::vector<std::uint8_t> step((m + 1) * W, kLeft);
    std::vector<int> prev(W), cur(W);
    for (std::size_t j = 0; j <= n; ++j) prev[j] = static_cast<int>(j) * kGap;
    for (std::size_t i = 1; i <= m; ++i) {
        cur[0] = static_cast<int>(i) * kGap;
        step[i * W] = kUp;
        for (std::size_t j = 1; j <= n; ++j) {
            const bool eq = R[r0 + i - 1].token == H[h0 + j - 1].token;
            const int diag = eq ? prev[j - 1] + kMatch : kNeg;
            const int up = prev[j] + kGap;
            const int left = cur[j - 1] + kGap;
            const int best = std::max(diag, std::max(up, left));
            cur[j] = best;
            step[i * W + j] = (eq && best == prev[j - 1] + kMatch) ? kDiag : best == up ? kUp : kLeft;
        }
        std::swap(prev, cur);
    }

    // Traceback into columns (ri, hi) into R/H absolute indices; -1 means a gap on that side.
    struct Col { long long ri; long long hi; };
    std::vector<Col> cols;
    std::size_t i = m, j = n;
    while (i > 0 || j > 0) {
        const std::uint8_t st = step[i * W + j];
        if (st == kDiag) {
            cols.push_back({static_cast<long long>(r0 + i - 1), static_cast<long long>(h0 + j - 1)});
            --i; --j;
        } else if (st == kUp) {
            cols.push_back({static_cast<long long>(r0 + i - 1), -1});
            --i;
        } else {
            cols.push_back({-1, static_cast<long long>(h0 + j - 1)});
            --j;
        }
    }
    std::reverse(cols.begin(), cols.end());
    if (hap_match != nullptr) {
        for (const Col& c : cols)
            if (c.ri >= 0 && c.hi >= 0) (*hap_match)[static_cast<std::size_t>(c.hi)] = R[static_cast<std::size_t>(c.ri)].src_idx;
    }

    // Emit events from maximal gap blocks -> DEL / INS / INV. Track the last matched
    // ref node as the anchor for an INS that follows it.
    std::string last_ref_node = preceding_ref_node;
    // The anchor for an INS that opens this segment is the shared anchor the previous segment ended
    // on, so its token index arrives with it; segments are cut at anchors, so this is most insertions.
    std::size_t last_ref_tok = preceding_ref_tok;
    // The tokens after the last matched column, where a block's empty side sits.
    std::size_t next_r = r0, next_h = h0;
    // One gap block, typed as INV / SUB / DEL / INS. An insertion is anchored on ins_node / ins_tok.
    auto emit_block = [&](std::vector<const Tok*> ref_blk, std::vector<const Tok*> hap_blk,
                          const std::string& ins_node, std::size_t ins_tok) {
        if (ref_blk.empty() && hap_blk.empty()) return;
        const std::size_t n_before = events.size();
        const std::size_t rl = ref_blk.empty() ? next_r : static_cast<std::size_t>(ref_blk.front() - R.data());
        const std::size_t rh = ref_blk.empty() ? rl : static_cast<std::size_t>(ref_blk.back() - R.data()) + 1;
        const std::size_t hl = hap_blk.empty() ? next_h : static_cast<std::size_t>(hap_blk.front() - H.data());
        const std::size_t hh = hap_blk.empty() ? hl : static_cast<std::size_t>(hap_blk.back() - H.data()) + 1;
        struct Place { std::vector<Event>& ev; std::size_t from, rl, rh, hl, hh;
            ~Place() { for (std::size_t k = from; k < ev.size(); ++k) { ev[k].r_lo = rl; ev[k].r_hi = rh; ev[k].h_lo = hl; ev[k].h_hi = hh; } } }
            place{events, n_before, rl, rh, hl, hh};
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
            // A gap block with both ref and hap content is a substitution: one event replacing the
            // reference run with the haplotype's, since a separate DEL and INS would each misstate the
            // haplotype's length change.
            if (!ref_blk.empty() && !hap_blk.empty()) {
                // Different nodes spelling the same sequence (a graph that kept two copies of one
                // sequence apart) are no change at all.
                const std::string rs = spell_toks(graph, ref_blk);
                const std::string hs = spell_toks(graph, hap_blk);
                if (rs == hs) { ref_blk.clear(); hap_blk.clear(); return; }
                Event e;
                e.type = EvType::Sub;
                for (const Tok* t : ref_blk) e.nodes.push_back(t->node_id);
                for (const Tok* t : hap_blk) e.nodes.push_back(t->node_id);
                e.start_node = ref_blk.front()->node_id;
                e.end_node = ref_blk.back()->node_id;
                e.seq = rs;
                e.alt_seq = hs;
                e.del_bp = toks_bp(graph, ref_blk);
                e.ins_bp = toks_bp(graph, hap_blk);
                // Sized by what changed (see substitution_change_bp), as the alignment step sizes it.
                e.size_bp = substitution_change_bp(rs, hs);
                e.anchor_node = ref_blk.front()->node_id;
                e.ref_tok_first = ref_blk.front()->src_idx;
                e.ref_tok_last = ref_blk.back()->src_idx;
                events.push_back(std::move(e));
            } else if (!ref_blk.empty()) {
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
                events.push_back(std::move(e));
            } else {
                Event e;
                e.type = EvType::Ins;
                for (const Tok* t : hap_blk) e.nodes.push_back(t->node_id);
                e.start_node = hap_blk.front()->node_id;
                e.end_node = hap_blk.back()->node_id;
                e.seq = spell_toks(graph, hap_blk);
                e.size_bp = toks_bp(graph, hap_blk);
                e.anchor_node = ins_node;
                e.ref_tok_anchor = ins_tok;
                e.anchor_after = true;
                events.push_back(std::move(e));
            }
        }
    };
    // A gap block is typed as one event only when each side's tokens are one stretch of its walk. Copies a
    // record counts leave both walks before they are compared (collapse_walk), so tokens on either side of
    // removed copies sit next to each other in the block although they lie a copy or more apart, and as
    // one event they would be placed at the first junction and spell the bases of every place. Such a
    // block is split where copies were removed, and the pieces are paired
    // by their place relative to the copies: the pieces before the first copy with each other, those after
    // the last copy with each other; a piece between two copies, a junction, stands alone.
    auto flush_block = [&](std::vector<const Tok*>& ref_blk, std::vector<const Tok*>& hap_blk) {
        if (ref_blk.empty() && hap_blk.empty()) return;
        struct Piece { std::vector<const Tok*> toks; bool lead = false, tail = false; };
        // Split one side at removed steps; lead = no copy removed before it within the block's place,
        // tail = none after it. `any_gap` is set when a copy was removed before, inside or after the block.
        auto pieces_of = [](const std::vector<const Tok*>& blk, const std::vector<Tok>& W, bool* any_gap) {
            std::vector<Piece> out;
            if (blk.empty()) return out;
            const std::size_t i0 = static_cast<std::size_t>(blk.front() - W.data());
            const std::size_t i1 = static_cast<std::size_t>(blk.back() - W.data());
            const bool gap_before = i0 > 0 && W[i0].src_idx - W[i0 - 1].src_idx > 1;
            const bool gap_after = i1 + 1 < W.size() && W[i1 + 1].src_idx - W[i1].src_idx > 1;
            out.push_back({});
            for (std::size_t k = 0; k < blk.size(); ++k) {
                if (k > 0 && blk[k]->src_idx - blk[k - 1]->src_idx > 1) out.push_back({});
                out.back().toks.push_back(blk[k]);
            }
            for (std::size_t k = 0; k < out.size(); ++k) {
                out[k].lead = k == 0 && !gap_before;
                out[k].tail = k + 1 == out.size() && !gap_after;
            }
            *any_gap = gap_before || gap_after || out.size() > 1;
            return out;
        };
        bool r_gap = false, h_gap = false;
        std::vector<Piece> rp = pieces_of(ref_blk, R, &r_gap);
        std::vector<Piece> hp = pieces_of(hap_blk, H, &h_gap);
        const auto single_place = [](const std::vector<Piece>& p) { return p.size() <= 1; };
        const auto compatible = [](const Piece& a, const Piece& b) {
            return (a.lead && a.tail) || (b.lead && b.tail) || (a.lead == b.lead && a.tail == b.tail);
        };
        if (!split_at_removed || (single_place(rp) && single_place(hp) &&
                                  (rp.empty() || hp.empty() || compatible(rp.front(), hp.front())))) {
            emit_block(ref_blk, hap_blk, last_ref_node, last_ref_tok);
            ref_blk.clear();
            hap_blk.clear();
            return;
        }
        const std::size_t n_before = events.size();
        // Which piece of `side` a piece with no copies on its own side (lead and tail at once) goes with:
        // the only end piece there is, else the one closest in length.
        auto pick_end = [&](const std::vector<Piece>& side, const Piece& free_piece) -> long long {
            long long lead = -1, tail = -1;
            for (std::size_t k = 0; k < side.size(); ++k) {
                if (side[k].lead) lead = static_cast<long long>(k);
                if (side[k].tail) tail = static_cast<long long>(k);
            }
            if (lead < 0 || tail < 0) return lead >= 0 ? lead : tail;
            const long long n = static_cast<long long>(toks_bp(graph, free_piece.toks));
            const long long dl = std::llabs(static_cast<long long>(toks_bp(graph, side[static_cast<std::size_t>(lead)].toks)) - n);
            const long long dt = std::llabs(static_cast<long long>(toks_bp(graph, side[static_cast<std::size_t>(tail)].toks)) - n);
            return dl <= dt ? lead : tail;
        };
        std::vector<long long> r_mate(rp.size(), -1), h_mate(hp.size(), -1);
        for (std::size_t a = 0; a < rp.size(); ++a)
            for (std::size_t b = 0; b < hp.size(); ++b) {
                if (r_mate[a] >= 0 || h_mate[b] >= 0) continue;
                const bool r_free = rp[a].lead && rp[a].tail, h_free = hp[b].lead && hp[b].tail;
                if (r_free || h_free) continue;
                if ((rp[a].lead && hp[b].lead) || (rp[a].tail && hp[b].tail)) { r_mate[a] = b; h_mate[b] = a; }
            }
        for (std::size_t a = 0; a < rp.size(); ++a)
            if (rp[a].lead && rp[a].tail && !hp.empty()) {
                const long long b = pick_end(hp, rp[a]);
                if (b >= 0 && h_mate[static_cast<std::size_t>(b)] < 0) { r_mate[a] = b; h_mate[static_cast<std::size_t>(b)] = static_cast<long long>(a); }
            }
        for (std::size_t b = 0; b < hp.size(); ++b)
            if (hp[b].lead && hp[b].tail && h_mate[b] < 0 && !rp.empty()) {
                const long long a = pick_end(rp, hp[b]);
                if (a >= 0 && r_mate[static_cast<std::size_t>(a)] < 0) { h_mate[b] = a; r_mate[static_cast<std::size_t>(a)] = static_cast<long long>(b); }
            }
        // Walk order: unpaired haplotype pieces before the copies, the reference pieces with their mates,
        // then the remaining unpaired haplotype pieces. An unpaired insertion anchors where any insertion
        // beside a counted array does, on the step before the block (the array's start).
        for (std::size_t b = 0; b < hp.size(); ++b)
            if (h_mate[b] < 0 && hp[b].lead) emit_block({}, hp[b].toks, last_ref_node, last_ref_tok);
        for (std::size_t a = 0; a < rp.size(); ++a)
            emit_block(rp[a].toks, r_mate[a] >= 0 ? hp[static_cast<std::size_t>(r_mate[a])].toks : std::vector<const Tok*>(),
                       last_ref_node, last_ref_tok);
        for (std::size_t b = 0; b < hp.size(); ++b)
            if (h_mate[b] < 0 && !hp[b].lead) emit_block({}, hp[b].toks, last_ref_node, last_ref_tok);
        for (std::size_t k = n_before; k < events.size(); ++k) events[k].across_copies = true;
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
            next_r = static_cast<std::size_t>(c.ri) + 1;
            next_h = static_cast<std::size_t>(c.hi) + 1;
        } else if (c.ri >= 0) {
            ref_blk.push_back(&R[static_cast<std::size_t>(c.ri)]);
        } else {
            hap_blk.push_back(&H[static_cast<std::size_t>(c.hi)]);
        }
    }
    flush_block(ref_blk, hap_blk);
}

// One replacement the comparison cut in pieces where the two walks share a few short nodes by chance,
// such as a lost copy read as a deletion and a nearby insertion. A piece whose length change reaches
// `floor` is joined with the next such piece when that one changes length the other way, both lie in one
// place on each walk (no copy removed between them), and only matches shorter than `join_bp` and pieces
// under `floor` separate them; the result is one substitution. A small piece beside a large one does not
// join it, and pieces that add up (two deletions) stay apart.
void join_split_replacements(const Graph& graph, const std::vector<Tok>& R, const std::vector<Tok>& H,
                             std::vector<Event>& events, std::size_t join_bp, std::size_t floor) {
    if (join_bp == 0 || events.size() < 2) return;
    const auto joinable = [](const Event& e) {
        return !e.seq_event && !e.across_copies && e.r_lo != SIZE_MAX && e.h_lo != SIZE_MAX &&
               (e.type == EvType::Del || e.type == EvType::Ins || e.type == EvType::Sub);
    };
    const auto net = [](const Event& e) -> long long {
        if (e.type == EvType::Ins) return static_cast<long long>(e.size_bp);
        if (e.type == EvType::Del) return -static_cast<long long>(e.size_bp);
        return static_cast<long long>(e.ins_bp) - static_cast<long long>(e.del_bp);
    };
    const auto big = [&](const Event& e) { return std::llabs(net(e)) >= static_cast<long long>(std::max<std::size_t>(1, floor)); };
    // the tokens between two pieces are a match: the same tokens on both walks, shorter than join_bp
    const auto short_match = [&](const Event& a, const Event& b) {
        if (b.r_lo < a.r_hi || b.h_lo < a.h_hi || b.r_lo - a.r_hi != b.h_lo - a.h_hi) return false;
        std::size_t bp = 0;
        for (std::size_t k = 0; k < b.r_lo - a.r_hi; ++k) {
            if (R[a.r_hi + k].token != H[a.h_hi + k].token) return false;
            bp += node_len(graph, R[a.r_hi + k].node_id);
        }
        return bp < join_bp;
    };
    // from a's start to b's end, consecutive steps of each walk: no counted copy left the walks there
    const auto one_place = [](const std::vector<Tok>& W, std::size_t lo, std::size_t hi) {
        for (std::size_t k = lo; k + 1 < hi; ++k)
            if (W[k + 1].src_idx != W[k].src_idx + 1) return false;
        return true;
    };
    std::vector<Event> out;
    std::size_t i = 0;
    while (i < events.size()) {
        std::size_t m = SIZE_MAX;   // the partner piece
        if (joinable(events[i]) && big(events[i]))
            for (std::size_t j = i + 1; j < events.size(); ++j) {
                if (!joinable(events[j]) || !short_match(events[j - 1], events[j]) ||
                    !one_place(R, events[i].r_lo, events[j].r_hi) || !one_place(H, events[i].h_lo, events[j].h_hi))
                    break;
                if (!big(events[j])) continue;
                if ((net(events[j]) > 0) != (net(events[i]) > 0)) m = j;
                break;
            }
        if (m == SIZE_MAX) { out.push_back(std::move(events[i])); ++i; continue; }
        std::vector<const Tok*> rt, ht;
        for (std::size_t k = events[i].r_lo; k < events[m].r_hi; ++k) rt.push_back(&R[k]);
        for (std::size_t k = events[i].h_lo; k < events[m].h_hi; ++k) ht.push_back(&H[k]);
        const std::string rs = spell_toks(graph, rt), hs = spell_toks(graph, ht);
        if (rs != hs) {   // the same bases either way: no change at all
            Event e;
            e.type = EvType::Sub;
            for (const Tok* t : rt) e.nodes.push_back(t->node_id);
            for (const Tok* t : ht) e.nodes.push_back(t->node_id);
            e.start_node = rt.front()->node_id;
            e.end_node = rt.back()->node_id;
            e.seq = rs;
            e.alt_seq = hs;
            e.del_bp = toks_bp(graph, rt);
            e.ins_bp = toks_bp(graph, ht);
            e.size_bp = substitution_change_bp(rs, hs);
            e.anchor_node = rt.front()->node_id;
            e.ref_tok_first = rt.front()->src_idx;
            e.ref_tok_last = rt.back()->src_idx;
            e.r_lo = events[i].r_lo; e.r_hi = events[m].r_hi; e.h_lo = events[i].h_lo; e.h_hi = events[m].h_hi;
            out.push_back(std::move(e));
        }
        i = m + 1;
    }
    events = std::move(out);
}

// Two pieces of one change at a collapsed module at most this far apart are joined when each alone is under
// --min-sv-bp (join_small_pieces): the distance at which benchmark reads two differences as one.
constexpr std::size_t kSmallPieceGapBp = 10;

// A walk's copy span at a collapsed module: its first and last token on one of the module's nodes, or false.
bool copy_span(const std::vector<Tok>& W, const std::unordered_set<std::string>& module_nodes, std::size_t& a,
               std::size_t& b) {
    a = SIZE_MAX; b = 0;
    for (std::size_t k = 0; k < W.size(); ++k)
        if (module_nodes.count(W[k].node_id)) { if (a == SIZE_MAX) a = k; b = k; }
    return a != SIZE_MAX;
}

// Whether tokens [lo, hi) of a walk lie inside its copy span [a, b]; an empty range sits before token lo, so
// it is inside when it falls between two tokens of the span.
bool in_copy_span(std::size_t lo, std::size_t hi, std::size_t a, std::size_t b) {
    if (lo == SIZE_MAX || a == SIZE_MAX) return false;
    if (lo == hi) return a < lo && lo <= b;
    return a <= lo && hi - 1 <= b;
}

// One deletion or one insertion the node comparison cut in pieces inside a collapsed module's copies, where the
// walks meet at short nodes by chance: two or more consecutive deletions (or insertions), separated only by
// stretches both walks share that are shorter than `join_bp` and lie in one place on each walk, all inside the
// copies on both walks, are one replacement of the reference from the first piece to the last by the
// haplotype's bases there (the shared stretches included, so it spells the haplotype exactly), sized by its
// net change, when that reaches `min_bp`. Replacements are not joined: their two sides say nothing about which
// way the change goes.
void join_copy_pieces(const Graph& graph, const std::vector<Tok>& R, const std::vector<Tok>& H,
                      const std::unordered_set<std::string>& module_nodes, std::vector<Event>& events,
                      std::size_t join_bp, std::size_t min_bp) {
    if (join_bp == 0 || events.size() < 2) return;
    std::size_t ra = 0, rb = 0, ha = 0, hb = 0;
    if (!copy_span(R, module_nodes, ra, rb) || !copy_span(H, module_nodes, ha, hb)) return;
    const auto joinable = [&](const Event& e) {
        return !e.seq_event && !e.across_copies && !e.split_piece && e.r_lo != SIZE_MAX && e.h_lo != SIZE_MAX &&
               (e.type == EvType::Del || e.type == EvType::Ins) && in_copy_span(e.r_lo, e.r_hi, ra, rb) &&
               in_copy_span(e.h_lo, e.h_hi, ha, hb);
    };
    // the stretch between two pieces: the same tokens on both walks, shorter than join_bp
    const auto short_match = [&](const Event& a, const Event& b) {
        if (b.r_lo < a.r_hi || b.h_lo < a.h_hi || b.r_lo - a.r_hi != b.h_lo - a.h_hi) return false;
        std::size_t bp = 0;
        for (std::size_t k = 0; k < b.r_lo - a.r_hi; ++k) {
            if (R[a.r_hi + k].token != H[a.h_hi + k].token) return false;
            bp += node_len(graph, R[a.r_hi + k].node_id);
        }
        return bp < join_bp;
    };
    // consecutive steps of the walk the tokens were collapsed from: one place, not across dropped copies
    const auto one_place = [](const std::vector<Tok>& W, std::size_t lo, std::size_t hi) {
        for (std::size_t k = lo; k + 1 < hi; ++k)
            if (W[k + 1].src_idx != W[k].src_idx + 1) return false;
        return true;
    };
    std::vector<Event> out;
    std::size_t i = 0;
    while (i < events.size()) {
        std::size_t m = i;
        if (joinable(events[i]))
            for (std::size_t j = i + 1; j < events.size(); ++j) {
                if (!joinable(events[j]) || events[j].type != events[i].type || !short_match(events[j - 1], events[j]) ||
                    !one_place(R, events[i].r_lo, events[j].r_hi) || !one_place(H, events[i].h_lo, events[j].h_hi))
                    break;
                m = j;
            }
        std::size_t net = 0;
        for (std::size_t k = i; k <= m; ++k) net += events[k].size_bp;
        if (m == i || net < min_bp) { out.push_back(std::move(events[i])); ++i; continue; }
        std::vector<const Tok*> rt, ht;
        for (std::size_t k = events[i].r_lo; k < events[m].r_hi; ++k) rt.push_back(&R[k]);
        for (std::size_t k = events[i].h_lo; k < events[m].h_hi; ++k) ht.push_back(&H[k]);
        Event e;
        e.type = EvType::Sub;
        for (const Tok* t : rt) e.nodes.push_back(t->node_id);
        for (const Tok* t : ht) e.nodes.push_back(t->node_id);
        e.start_node = rt.front()->node_id;
        e.end_node = rt.back()->node_id;
        e.seq = spell_toks(graph, rt);
        e.alt_seq = spell_toks(graph, ht);
        e.del_bp = toks_bp(graph, rt);
        e.ins_bp = toks_bp(graph, ht);
        e.size_bp = e.del_bp > e.ins_bp ? e.del_bp - e.ins_bp : e.ins_bp - e.del_bp;
        e.anchor_node = rt.front()->node_id;
        e.ref_tok_first = rt.front()->src_idx;
        e.ref_tok_last = rt.back()->src_idx;
        e.r_lo = events[i].r_lo; e.r_hi = events[m].r_hi; e.h_lo = events[i].h_lo; e.h_hi = events[m].h_hi;
        e.split_piece = true;   // placed by its own tokens: a copy taken out beside it leaves no match between
        out.push_back(std::move(e));
        i = m + 1;
    }
    events = std::move(out);
}

// Whether `piece` is a version of part of `copy`: all of it aligns inside `copy` at 80% identity or more. Two
// versions of one segment can share few k-mers yet align well.
bool version_of(const std::string& piece, const std::string& copy) {
    if (piece.empty() || copy.size() < piece.size() / 2) return false;
    return infix_distance(piece, copy, piece.size() / 5) != SIZE_MAX;
}

// At a collapsed module, a replacement whose one side holds a copy of the module's anchor and whose other
// side is SV-sized new sequence is a copy change beside an insertion or deletion of its own, as when a lost
// copy and a novel insertion fall in one gap, or are joined across a few nodes the walks share by chance. It is
// split into a deletion and an insertion, cut at those shared nodes, so the copy-number record can account for
// the copy while the other side stays an event. Only for a haplotype whose copy number changes that way
// (`copy_delta`, its copies less the reference's), when the side holding the anchor is about the copies changed
// (half to one and a quarter times `unit_bp` each, as the copy-change step judges an event), and when the other
// side is new: under half of its 15-mers in the copy side, and not a version of part of it (version_of) about
// the size of what the copy side holds beyond the copies changed, so a copy side whose extra part the haplotype
// holds as another version stays one replacement.
void split_copy_replacements(const Graph& graph, const std::vector<Tok>& R, const std::vector<Tok>& H,
                             const std::unordered_set<std::string>& anchor, long long copy_delta, std::size_t unit_bp,
                             std::size_t min_sv_bp, std::vector<Event>& events) {
    if (copy_delta == 0) return;
    const auto kmers15 = [](const std::string& s) {
        std::unordered_set<std::uint64_t> out;
        std::uint64_t f = 0, r = 0;
        std::size_t len = 0;
        const std::uint64_t mask = (1ULL << 30) - 1;
        for (const char c : s) {
            const int b = c == 'A' ? 0 : c == 'C' ? 1 : c == 'G' ? 2 : c == 'T' ? 3 : -1;
            if (b < 0) { len = 0; continue; }
            f = ((f << 2) | static_cast<std::uint64_t>(b)) & mask;
            r = (r >> 2) | (static_cast<std::uint64_t>(3 - b) << 28);
            if (++len >= 15) out.insert(std::min(f, r));
        }
        return out;
    };
    std::vector<Event> out;
    out.reserve(events.size());
    for (Event& e : events) {
        const auto keep = [&] { out.push_back(std::move(e)); };
        if (e.type != EvType::Sub || e.seq_event || e.r_lo == SIZE_MAX || e.h_lo == SIZE_MAX || e.r_lo == 0 ||
            e.r_hi > R.size() || e.h_hi > H.size() || e.r_lo >= e.r_hi || e.h_lo >= e.h_hi) { keep(); continue; }
        const std::size_t m = e.r_hi - e.r_lo, n = e.h_hi - e.h_lo;
        if (m * n > 4000000) { keep(); continue; }
        // the tokens both sides share, in order (a longest common subsequence)
        std::vector<std::uint32_t> lcs((m + 1) * (n + 1), 0);
        for (std::size_t a = m; a-- > 0;)
            for (std::size_t b = n; b-- > 0;)
                lcs[a * (n + 1) + b] = R[e.r_lo + a].token == H[e.h_lo + b].token
                    ? lcs[(a + 1) * (n + 1) + b + 1] + 1
                    : std::max(lcs[(a + 1) * (n + 1) + b], lcs[a * (n + 1) + b + 1]);
        std::vector<char> r_shared(m, 0), h_shared(n, 0);
        std::vector<std::pair<std::size_t, std::size_t>> pairs;   // absolute (ref, hap) token indices
        for (std::size_t a = 0, b = 0; a < m && b < n;) {
            if (R[e.r_lo + a].token == H[e.h_lo + b].token) {
                r_shared[a] = h_shared[b] = 1; pairs.emplace_back(e.r_lo + a, e.h_lo + b); ++a; ++b;
            } else if (lcs[(a + 1) * (n + 1) + b] >= lcs[a * (n + 1) + b + 1]) ++a;
            else ++b;
        }
        // what is left on each side must be one run, the copy on one side and the new sequence on the other
        const auto run = [](const std::vector<char>& shared, std::size_t& lo, std::size_t& hi) {
            lo = SIZE_MAX; hi = 0;
            for (std::size_t k = 0; k < shared.size(); ++k)
                if (!shared[k]) { if (lo == SIZE_MAX) lo = k; hi = k + 1; }
            if (lo == SIZE_MAX) return false;
            for (std::size_t k = lo; k < hi; ++k) if (shared[k]) return false;
            return true;
        };
        std::size_t rl = 0, rh = 0, hl = 0, hh = 0;
        if (!run(r_shared, rl, rh) || !run(h_shared, hl, hh)) { keep(); continue; }
        rl += e.r_lo; rh += e.r_lo; hl += e.h_lo; hh += e.h_lo;
        std::vector<const Tok*> rb, hb;
        for (std::size_t k = rl; k < rh; ++k) rb.push_back(&R[k]);
        for (std::size_t k = hl; k < hh; ++k) hb.push_back(&H[k]);
        const auto holds_anchor = [&](const std::vector<const Tok*>& b) {
            for (const Tok* t : b) if (anchor.count(t->node_id)) return true;
            return false;
        };
        const bool ra = holds_anchor(rb), ha = holds_anchor(hb);
        const std::size_t rbp = toks_bp(graph, rb), hbp = toks_bp(graph, hb);
        const double copies_bp = static_cast<double>(std::llabs(copy_delta)) * static_cast<double>(unit_bp);
        const double side_bp = static_cast<double>(ra ? rbp : hbp);
        if (ra == ha || (ra ? copy_delta > 0 : copy_delta < 0) || (ra ? hbp : rbp) < std::max<std::size_t>(1, min_sv_bp) ||
            side_bp < 0.5 * copies_bp || side_bp > 1.25 * copies_bp) { keep(); continue; }
        const std::string rs = spell_toks(graph, rb), hs = spell_toks(graph, hb);
        {
            const auto copy_k = kmers15(ra ? rs : hs), other_k = kmers15(ra ? hs : rs);
            std::size_t in = 0;
            for (const std::uint64_t k : other_k) in += copy_k.count(k);
            const double excess = side_bp - copies_bp, obp = static_cast<double>(ra ? hbp : rbp);
            const bool fits = excess > 0 && obp >= 0.5 * excess && obp <= 1.25 * excess;
            if (other_k.empty() || 2 * in >= other_k.size() || (fits && version_of(ra ? hs : rs, ra ? rs : hs))) { keep(); continue; }
        }
        // the insertion sits after the last shared token before it, or where the replacement began
        std::size_t anchor_r = e.r_lo - 1, before_h = e.h_lo;
        for (const auto& [pr, ph] : pairs) if (ph < hl) anchor_r = pr;
        for (const auto& [pr, ph] : pairs) if (pr < rl) before_h = ph + 1;
        Event d;
        d.type = EvType::Del;
        for (const Tok* t : rb) d.nodes.push_back(t->node_id);
        d.start_node = rb.front()->node_id;
        d.end_node = rb.back()->node_id;
        d.seq = rs;
        d.size_bp = rbp;
        d.anchor_node = rb.front()->node_id;
        d.ref_tok_first = rb.front()->src_idx;
        d.ref_tok_last = rb.back()->src_idx;
        d.r_lo = rl; d.r_hi = rh; d.h_lo = d.h_hi = before_h;
        Event a;
        a.type = EvType::Ins;
        for (const Tok* t : hb) a.nodes.push_back(t->node_id);
        a.start_node = hb.front()->node_id;
        a.end_node = hb.back()->node_id;
        a.seq = hs;
        a.size_bp = hbp;
        a.anchor_node = R[anchor_r].node_id;
        a.ref_tok_anchor = R[anchor_r].src_idx;
        a.anchor_after = true;
        a.r_lo = a.r_hi = anchor_r + 1; a.h_lo = hl; a.h_hi = hh;
        d.across_copies = a.across_copies = e.across_copies;
        d.split_piece = a.split_piece = true;
        if (anchor_r + 1 <= rl) { out.push_back(std::move(a)); out.push_back(std::move(d)); }
        else { out.push_back(std::move(d)); out.push_back(std::move(a)); }
    }
    events = std::move(out);
}

// A change of a walk, tokens [rl, rh) of R replaced by [hl, hh) of H, as an event: a deletion, an insertion or
// a replacement, with the fields the node comparison gives its own.
Event token_event(const Graph& graph, const std::vector<Tok>& R, const std::vector<Tok>& H, std::size_t rl,
                  std::size_t rh, std::size_t hl, std::size_t hh) {
    std::vector<const Tok*> rt, ht;
    for (std::size_t k = rl; k < rh; ++k) rt.push_back(&R[k]);
    for (std::size_t k = hl; k < hh; ++k) ht.push_back(&H[k]);
    Event e;
    if (!rt.empty() && !ht.empty()) {
        e.type = EvType::Sub;
        e.seq = spell_toks(graph, rt);
        e.alt_seq = spell_toks(graph, ht);
        e.del_bp = toks_bp(graph, rt);
        e.ins_bp = toks_bp(graph, ht);
        e.size_bp = substitution_change_bp(e.seq, e.alt_seq);
    } else if (!rt.empty()) {
        e.type = EvType::Del;
        e.seq = spell_toks(graph, rt);
        e.size_bp = toks_bp(graph, rt);
    } else {
        e.type = EvType::Ins;
        e.seq = spell_toks(graph, ht);
        e.size_bp = toks_bp(graph, ht);
    }
    for (const Tok* t : rt) e.nodes.push_back(t->node_id);
    for (const Tok* t : ht) e.nodes.push_back(t->node_id);
    if (!rt.empty()) {
        e.start_node = rt.front()->node_id;
        e.end_node = rt.back()->node_id;
        e.anchor_node = rt.front()->node_id;
        e.ref_tok_first = rt.front()->src_idx;
        e.ref_tok_last = rt.back()->src_idx;
    } else {
        e.start_node = ht.front()->node_id;
        e.end_node = ht.back()->node_id;
        if (rl > 0) { e.anchor_node = R[rl - 1].node_id; e.ref_tok_anchor = R[rl - 1].src_idx; }
        e.anchor_after = true;
    }
    e.r_lo = rl; e.r_hi = rh; e.h_lo = hl; e.h_hi = hh;
    e.split_piece = true;   // placed by its own tokens: a piece cut from a change, no match beside it
    return e;
}

// Pieces of one change the node comparison cut apart at a collapsed module, each under `min_bp` and so dropped
// on its own: consecutive pieces of any kind, separated only by stretches both walks share that are `gap_bp`
// or shorter and lie in one place on each walk, with one of them touching the copies, are one replacement from
// the first piece to the last (the shared stretches included, so it spells the haplotype exactly), when it
// reaches `min_bp`.
void join_small_pieces(const Graph& graph, const std::vector<Tok>& R, const std::vector<Tok>& H,
                       const std::unordered_set<std::string>& module_nodes, std::vector<Event>& events,
                       std::size_t gap_bp, std::size_t min_bp) {
    if (events.size() < 2) return;
    std::size_t ra = 0, rb = 0, ha = 0, hb = 0;
    const bool rspan = copy_span(R, module_nodes, ra, rb), hspan = copy_span(H, module_nodes, ha, hb);
    const auto small = [&](const Event& e) {
        return !e.seq_event && e.r_lo != SIZE_MAX && e.h_lo != SIZE_MAX && e.size_bp < min_bp &&
               (e.type == EvType::Del || e.type == EvType::Ins || e.type == EvType::Sub);
    };
    const auto touches = [&](const Event& e) {
        const bool r = rspan && (e.r_lo < e.r_hi ? e.r_lo <= rb && e.r_hi > ra : in_copy_span(e.r_lo, e.r_hi, ra, rb));
        const bool h = hspan && (e.h_lo < e.h_hi ? e.h_lo <= hb && e.h_hi > ha : in_copy_span(e.h_lo, e.h_hi, ha, hb));
        return r || h;
    };
    const auto near = [&](const Event& a, const Event& b) {
        if (b.r_lo < a.r_hi || b.h_lo < a.h_hi || b.r_lo - a.r_hi != b.h_lo - a.h_hi) return false;
        std::size_t bp = 0;
        for (std::size_t k = 0; k < b.r_lo - a.r_hi; ++k) {
            if (R[a.r_hi + k].token != H[a.h_hi + k].token) return false;
            bp += node_len(graph, R[a.r_hi + k].node_id);
        }
        return bp <= gap_bp;
    };
    const auto one_place = [](const std::vector<Tok>& W, std::size_t lo, std::size_t hi) {
        for (std::size_t k = lo; k + 1 < hi; ++k)
            if (W[k + 1].src_idx != W[k].src_idx + 1) return false;
        return true;
    };
    std::vector<Event> out;
    std::size_t i = 0;
    while (i < events.size()) {
        std::size_t m = i;
        bool touch = small(events[i]) && touches(events[i]);
        if (small(events[i]))
            for (std::size_t j = i + 1; j < events.size(); ++j) {
                if (!small(events[j]) || !near(events[j - 1], events[j]) || !one_place(R, events[i].r_lo, events[j].r_hi) ||
                    !one_place(H, events[i].h_lo, events[j].h_hi))
                    break;
                m = j;
                touch = touch || touches(events[j]);
            }
        if (m > i && touch && events[i].r_lo < events[m].r_hi && events[i].h_lo < events[m].h_hi) {
            Event e = token_event(graph, R, H, events[i].r_lo, events[m].r_hi, events[i].h_lo, events[m].h_hi);
            if (e.size_bp >= min_bp) { out.push_back(std::move(e)); i = m + 1; continue; }
        }
        out.push_back(std::move(events[i]));
        ++i;
    }
    events = std::move(out);
}

// One haplotype's copies at a collapsed module, told apart from its other events (see call.md, step 7).
// A copy as an edit of the reference: [start, end] replaced by seq (reference orientation). `code` spells seq
// by the reference stretches it repeats (=A-B forward, ~A-B reverse-complemented) and literal bases between,
// joined by '+'; empty where that is no shorter. A gained copy repeats the reference's copies, so the code stays
// short however long the copy.
struct CopyForm { std::size_t start = 0, end = 0; std::string seq, code; std::size_t literal_bp = 0; };
struct CopySplit {
    enum Mode { None, Copies, OwnEvent } mode = None;
    std::vector<Event> rest;        // Copies: its events with the copies taken out, in walk order
    std::vector<CopyForm> forms;    // Copies: the copies, as edits of the reference
    long long passes = 0;           // Copies: anchor passes the forms gain (lose: negative)
};

// The copies a haplotype gained or lost, from its own events (`events`, the node comparison of its walk H with
// the reference's R). `delta` is its copies less the reference's, `module_bp` its change over the module's
// nodes, `unit` one copy (CN_UNIT_BP).
//   - The copy span of a walk runs from its first to its last token on a module node. An event touching it is
//     part of the copy change; S is their net change.
//   - The haplotype is no carrier, the copy change being an event of its own, when the change beside the copies
//     goes their way and is more than the copies: (S - module_bp) has delta's sign and exceeds both
//     |module_bp| and |delta| copies. A deletion that takes a copy with it is a deletion.
//   - Otherwise each touching event is cut at the copy span's ends (node boundaries); the parts outside are
//     events of their own. A replacement whose two sides are versions of each other is not cut: a copy
//     replaced by a diverged version of it is one change. A part inside is a copy when it holds anchor passes
//     (what CN counts), either way. A part with none is a copy only by its bases: a gain's insertion that is a
//     version of the reference's copies, its first to last anchor pass (a copy's part on nodes the reference
//     does not walk), or, beside a loss, an insertion that is a version of what the loss took, about the excess
//     it took beyond the copies. A deletion with no anchor pass is a change inside a remaining copy.
//   - The copies are bounded by twice the larger of |module_bp| and |delta| copies: past it, versions go back
//     to the events, largest first.
// `place` writes a part as an edit of the reference, false where it cannot (the part then stays an event).
// `h_first`, `h_last`: the walk's copy span widened over sequence in place of the reference's end module nodes.
CopySplit split_copies(const Graph& graph, const std::vector<Tok>& R, const std::vector<Tok>& H,
                       const std::vector<Event>& events, const std::unordered_set<std::string>& module_nodes,
                       const std::unordered_set<std::string>& anchor_nodes, long long delta, long long module_bp,
                       std::size_t unit, const std::function<bool(const Event&, CopyForm&)>& place,
                       std::size_t h_first = SIZE_MAX, std::size_t h_last = SIZE_MAX) {
    CopySplit out;
    std::size_t ra = 0, rb = 0, ha = 0, hb = 0;
    const bool rspan = copy_span(R, module_nodes, ra, rb), hspan = copy_span(H, module_nodes, ha, hb);
    // sequence in place of the copies' end module nodes is in the copies, as CNBP counts it
    if (hspan && h_first != SIZE_MAX && h_first < ha) ha = h_first;
    if (hspan && h_last != SIZE_MAX && h_last > hb && h_last < H.size()) hb = h_last;
    const auto clampi = [](std::size_t x, std::size_t lo, std::size_t hi) { return std::min(std::max(x, lo), hi); };
    const auto net_of = [](const Event& e) -> long long {
        if (e.type == EvType::Ins) return static_cast<long long>(e.size_bp);
        if (e.type == EvType::Del) return -static_cast<long long>(e.size_bp);
        if (e.type == EvType::Sub) return static_cast<long long>(e.ins_bp) - static_cast<long long>(e.del_bp);
        return 0;
    };
    // Each event's cut points on both walks: before [lo, c1), inside [c1, c2), after [c2, hi). A walk with no
    // module node is not cut.
    struct Cut { std::size_t r1, r2, h1, h2; bool touches; };
    std::vector<Cut> cuts(events.size());
    long long S = 0;
    for (std::size_t i = 0; i < events.size(); ++i) {
        const Event& e = events[i];
        Cut& c = cuts[i];
        c.touches = false;
        if (e.seq_event || e.type == EvType::Inv || e.type == EvType::Dup || e.r_lo == SIZE_MAX || e.h_lo == SIZE_MAX) continue;
        c.r1 = rspan ? clampi(ra, e.r_lo, e.r_hi) : e.r_lo;
        c.r2 = rspan ? clampi(rb + 1, e.r_lo, e.r_hi) : e.r_hi;
        c.h1 = hspan ? clampi(ha, e.h_lo, e.h_hi) : e.h_lo;
        c.h2 = hspan ? clampi(hb + 1, e.h_lo, e.h_hi) : e.h_hi;
        // It touches the copies by its tokens inside a walk's copy span. A walk with no module node has no copy
        // span and is not cut: the event's tokens there go with whatever it touches on the other walk, and do
        // not make it touch, so a change elsewhere in the site stays out of the copies.
        c.touches = (rspan && c.r1 < c.r2) || (hspan && c.h1 < c.h2);
        if (c.touches) S += net_of(e);
    }
    const long long off = S - module_bp;
    const long long ad = delta < 0 ? -delta : delta;
    const long long scale = std::max(module_bp < 0 ? -module_bp : module_bp, ad * static_cast<long long>(unit));
    if (off * delta > 0 && (off < 0 ? -off : off) > scale) { out.mode = CopySplit::OwnEvent; return out; }
    out.mode = CopySplit::Copies;
    // the anchor passes a part gains: the length-weighted median over the anchor's nodes of its visits there
    const auto passes_of = [&](const Event& p) -> long long {
        std::unordered_map<std::string, long long> d;
        for (std::size_t k = p.r_lo; k < p.r_hi; ++k) if (anchor_nodes.count(R[k].node_id)) --d[R[k].node_id];
        for (std::size_t k = p.h_lo; k < p.h_hi; ++k) if (anchor_nodes.count(H[k].node_id)) ++d[H[k].node_id];
        std::vector<std::pair<long long, std::size_t>> cw;
        std::size_t total = 0;
        for (const std::string& n : anchor_nodes) {
            const auto it = d.find(n);
            const std::size_t len = node_len(graph, n);
            cw.emplace_back(it == d.end() ? 0 : it->second, len);
            total += len;
        }
        if (total == 0) return 0;
        std::sort(cw.begin(), cw.end());
        std::size_t acc = 0;
        for (const auto& [v, w] : cw) { acc += w; if (2 * acc >= total) return v; }
        return cw.back().first;
    };
    struct Part { Event e; bool inside = false; long long passes = 0; int role = 0; };   // role: 1 copy, 2 version
    std::vector<Part> parts;
    for (std::size_t i = 0; i < events.size(); ++i) {
        const Event& e = events[i];
        const Cut& c = cuts[i];
        if (!c.touches) { parts.push_back({e, false, 0, 0}); continue; }
        const std::size_t rs[4] = {e.r_lo, c.r1, c.r2, e.r_hi}, hs[4] = {e.h_lo, c.h1, c.h2, e.h_hi};
        // A replacement of a copy by a version of it is one change, not a copy and an insertion: not cut.
        const bool versions = e.type == EvType::Sub &&
                              (e.seq.size() <= e.alt_seq.size() ? version_of(e.seq, e.alt_seq) : version_of(e.alt_seq, e.seq));
        if (versions || (c.r1 == e.r_lo && c.r2 == e.r_hi && c.h1 == e.h_lo && c.h2 == e.h_hi)) {
            Part pt{e, true, 0, 0};
            pt.passes = passes_of(e);
            parts.push_back(std::move(pt));
            continue;
        }
        for (int k = 0; k < 3; ++k) {
            if (rs[k] == rs[k + 1] && hs[k] == hs[k + 1]) continue;
            Part pt{token_event(graph, R, H, rs[k], rs[k + 1], hs[k], hs[k + 1]), k == 1, 0, 0};
            if (pt.inside) pt.passes = passes_of(pt.e);
            parts.push_back(std::move(pt));
        }
    }
    long long copy_net = 0;
    for (Part& pt : parts)
        if (pt.inside && pt.passes != 0) { pt.role = 1; copy_net += net_of(pt.e); }
    if (delta > 0) {
        // a gain's part on nodes the reference does not walk: a version of the reference's copies (from its
        // first to its last pass through the anchor), or with no reference copy of the gained ones
        std::string copies;
        std::size_t aa = SIZE_MAX, ab = 0;
        for (std::size_t k = 0; k < R.size(); ++k) if (anchor_nodes.count(R[k].node_id)) { if (aa == SIZE_MAX) aa = k; ab = k; }
        if (aa != SIZE_MAX) { std::vector<const Tok*> t; for (std::size_t k = aa; k <= ab; ++k) t.push_back(&R[k]); copies = spell_toks(graph, t); }
        else for (const Part& pt : parts) if (pt.role == 1) copies += pt.e.type == EvType::Sub ? pt.e.alt_seq : pt.e.seq;
        for (Part& pt : parts) {
            if (!pt.inside || pt.role != 0 || net_of(pt.e) <= 0 || copies.empty()) continue;
            const std::string& ins = pt.e.type == EvType::Sub ? pt.e.alt_seq : pt.e.seq;
            if (version_of(ins, copies)) { pt.role = 2; copy_net += net_of(pt.e); }
        }
    } else if (copy_net < 0) {
        // beside a loss, an insertion putting back another version of what the loss took beyond its copies
        std::string lost;
        for (const Part& pt : parts) if (pt.role == 1 && pt.e.type != EvType::Ins) lost += pt.e.seq;
        double excess = static_cast<double>(-copy_net) - static_cast<double>(ad) * static_cast<double>(unit);
        std::vector<Part*> cand;
        for (Part& pt : parts)
            if (pt.inside && pt.role == 0 && net_of(pt.e) > 0 && pt.e.type != EvType::Del) cand.push_back(&pt);
        std::stable_sort(cand.begin(), cand.end(), [&](const Part* a, const Part* b) { return net_of(a->e) > net_of(b->e); });
        for (Part* pt : cand) {
            const double b = static_cast<double>(net_of(pt->e));
            if (excess <= 0 || b < 0.5 * excess || b > 1.25 * excess) continue;
            if (!version_of(pt->e.type == EvType::Sub ? pt->e.alt_seq : pt->e.seq, lost)) continue;
            pt->role = 2; copy_net += net_of(pt->e); excess -= b;
        }
    }
    // the copies bounded by the walk's own module change: versions go back first, the largest first
    const long long bound = 2 * scale;
    while ((copy_net < 0 ? -copy_net : copy_net) > bound) {
        Part* big = nullptr;
        for (Part& pt : parts)
            if (pt.role == 2 && (big == nullptr || std::llabs(net_of(pt.e)) > std::llabs(net_of(big->e)))) big = &pt;
        if (big == nullptr) break;
        big->role = 0; copy_net -= net_of(big->e);
    }
    std::size_t copies = 0;
    for (Part& pt : parts) {
        CopyForm f;
        if (pt.role != 0 && place(pt.e, f)) { out.forms.push_back(std::move(f)); out.passes += pt.passes; ++copies; }
        else { pt.e.copies_before = copies; out.rest.push_back(std::move(pt.e)); }
    }
    return out;
}


// DEL/INS/INV of haplotype vs reference. Split both walks at shared anchors (tokens unique in both,
// chained in monotonic order) so the quadratic DP runs only between anchors - bounds cost on huge
// bubbles and gives identical breakpoints across haplotypes so the same event merges. CN is separate.
std::vector<Event> diff_walks(
    const Graph& graph,
    const std::vector<Tok>& R,
    const std::vector<Tok>& H,
    const std::string& bubble_source,
    std::vector<std::size_t>* hap_match = nullptr,
    std::size_t* skipped = nullptr,
    bool split_at_removed = false,
    std::size_t join_bp = 0,
    std::size_t join_floor = 1) {

    // hap_match, when given, receives for each H token the reference src_idx it was aligned to
    // (SIZE_MAX where it is a gap): the occurrence-aware position of a haplotype step on the reference.
    std::vector<Event> events;
    if (hap_match != nullptr) hap_match->assign(H.size(), SIZE_MAX);
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
        diff_segment(graph, R, r_prev, r_end, H, h_prev, h_end, preceding_ref, preceding_ref_tok, events, hap_match, skipped,
                     split_at_removed);
        if (a < chain_ref.size()) {
            if (hap_match != nullptr) (*hap_match)[h_end] = R[r_end].src_idx;
            preceding_ref = R[r_end].node_id; // the anchor node
            preceding_ref_tok = R[r_end].src_idx;
            r_prev = r_end + 1;
            h_prev = h_end + 1;
        }
    }
    join_split_replacements(graph, R, H, events, join_bp, join_floor);
    return events;
}

// Within a haplotype, coalesce consecutive same-type DEL/INS/INV events whose gap is <= merge_distance_bp
// into one event (union node sets). A DEL or INV gap is measured in reference bp (by walk occurrence where
// known, else via the node->genomic-start map); an insertion's also in the haplotype's own sequence.
void coalesce_events(
    const Graph& graph,
    std::vector<Event>& events,
    const std::unordered_map<std::string, std::size_t>& ref_node_pos,
    const std::unordered_map<std::string, std::size_t>& hap_node_pos,
    std::size_t merge_distance_bp,
    const std::function<bool(std::size_t, std::size_t, long long&, long long&)>& tok_bp) {

    if (events.size() < 2) return;
    auto pos_in = [](const std::unordered_map<std::string, std::size_t>& m,
                     const std::string& node) -> long long {
        const auto it = m.find(node);
        return it == m.end() ? -1 : static_cast<long long>(it->second);
    };
    // Reference bp interval [lo, hi) that an event occupies on the reference.
    auto ref_span = [&](const Event& e, long long& lo, long long& hi) {
        // Occurrence-aware where the walk positions are known, so two events at distant visits to a
        // repeated node do not look adjacent through the node's first occurrence.
        const std::size_t t0 = e.type == EvType::Ins ? e.ref_tok_anchor : e.ref_tok_first;
        const std::size_t t1 = e.type == EvType::Ins ? e.ref_tok_anchor : e.ref_tok_last;
        if (t0 != SIZE_MAX && t1 != SIZE_MAX && tok_bp != nullptr) {
            long long a = -1, b = -1;
            if (tok_bp(t0, t1, a, b)) { lo = a; hi = b; return; }
        }
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
    // Haplotype-space bp span of the event's own nodes, used for insertions only: two INS join through it
    // where ins_anchor_ok below allows, the same reference anchor or an insertion with no anchor token. Not
    // for a DEL/INV: its nodes found elsewhere in the walk are other copies, which say nothing about where
    // on the reference the two events lie.
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
    // Reference distance between the two intervals, whichever comes first: where the reference crosses
    // the bubble sink->source, walk order is decreasing genomic order and an event starts before its
    // predecessor, and cur_lo - prev_hi would then count both events' lengths into the gap.
    auto ref_gap_ok = [&](long long prev_lo, long long prev_hi, long long cur_lo, long long cur_hi) {
        if (prev_lo >= 0 && cur_hi >= 0 && cur_lo < prev_lo) return gap_ok(cur_hi, prev_lo);
        return gap_ok(prev_hi, cur_lo);
    };
    std::vector<Event> out;
    long long prev_lo = -1, prev_hi = -1, prev_hhi = -1;
    for (Event& e : events) {
        long long lo = -1, hi = -1, hlo = -1, hhi = -1;
        ref_span(e, lo, hi);
        hap_span(e, hlo, hhi);
        // Merge when the same-type predecessor is near in reference space, or for insertions in this
        // haplotype's own sequence space. Two insertions at different reference anchors stay two, since
        // one record keeps one anchor: they join only when they share an anchor token, or when either
        // has none (nothing then says the anchors differ).
        const bool ins_anchor_ok =
            e.type != EvType::Ins || out.empty() ||
            out.back().ref_tok_anchor == SIZE_MAX || e.ref_tok_anchor == SIZE_MAX ||
            out.back().ref_tok_anchor == e.ref_tok_anchor;
        if (e.type != EvType::Dup && e.type != EvType::Sub && !out.empty() && out.back().type == e.type && ins_anchor_ok &&
            !e.across_copies && !out.back().across_copies && e.copies_before == out.back().copies_before &&
            (ref_gap_ok(prev_lo, prev_hi, lo, hi) || (e.type == EvType::Ins && gap_ok(prev_hhi, hlo)))) {
            Event& prev = out.back();
            for (const std::string& nd : e.nodes) prev.nodes.push_back(nd);
            prev.end_node = e.end_node;
            if (!e.seq.empty()) prev.seq += e.seq;
            prev.size_bp += e.size_bp;
            // The reference walk span has to grow with the content, or the record claims more deleted
            // bases than the interval it points at contains.
            if (e.ref_tok_first != SIZE_MAX)
                prev.ref_tok_first = prev.ref_tok_first == SIZE_MAX
                    ? e.ref_tok_first : std::min(prev.ref_tok_first, e.ref_tok_first);
            if (e.ref_tok_last != SIZE_MAX)
                prev.ref_tok_last = prev.ref_tok_last == SIZE_MAX
                    ? e.ref_tok_last : std::max(prev.ref_tok_last, e.ref_tok_last);
            if (lo >= 0) prev_lo = prev_lo < 0 ? lo : std::min(prev_lo, lo);
            if (hi >= 0) prev_hi = std::max(prev_hi, hi);
            if (hhi >= 0) prev_hhi = std::max(prev_hhi, hhi);
        } else {
            out.push_back(std::move(e));
            prev_lo = lo;
            prev_hi = hi;
            prev_hhi = hhi;
        }
    }
    events = std::move(out);
}

// Length-weighted Jaccard of two events' node multisets: a node counts as many times as the event
// walks it, so different copy counts of one repeat unit are not the same event. Events without
// repeated nodes score exactly as over sets.
double weighted_jaccard(
    const Graph& graph,
    const std::vector<std::string>& a,
    const std::vector<std::string>& b) {

    std::unordered_map<std::string, std::size_t> ca, cb;
    for (const std::string& n : a) ++ca[n];
    for (const std::string& n : b) ++cb[n];
    std::size_t inter = 0;
    std::size_t uni = 0;
    for (const auto& [n, x] : ca) {
        const std::size_t l = node_len(graph, n);
        const auto it = cb.find(n);
        const std::size_t y = it == cb.end() ? 0 : it->second;
        inter += std::min(x, y) * l;
        uni += std::max(x, y) * l;
    }
    for (const auto& [n, y] : cb)
        if (!ca.count(n)) uni += y * node_len(graph, n);
    return uni == 0 ? 0.0 : static_cast<double>(inter) / static_cast<double>(uni);
}

// Banded-alignment identity between two event sequences, gated by len_ratio (min shorter/longer to
// even compare - decoupled from min_id so same-motif STR alleles of different size can still merge).
double seq_identity(const std::string& a, const std::string& b, double min_id, double len_ratio) {
    if (a.empty() || b.empty()) return 0.0;
    const std::size_t lo = std::min(a.size(), b.size());
    const std::size_t hi = std::max(a.size(), b.size());
    if (static_cast<double>(lo) < len_ratio * static_cast<double>(hi)) return 0.0;
    // Symmetric: of two of one length, the lesser is fitted into the other, so events_match(a, b) and
    // (b, a) agree and a merge cannot depend on the order events were met in.
    const bool a_short = a.size() < b.size() || (a.size() == b.size() && a <= b);
    const std::string& shorter = a_short ? a : b;
    const std::string& longer = a_short ? b : a;
    const std::size_t band = std::max<std::size_t>(
        8, static_cast<std::size_t>((1.0 - min_id) * static_cast<double>(shorter.size())) + 8);
    const FitAlignResult fa = fit_align(shorter, longer, band);
    return fa.ok ? fa.identity : 0.0;
}

// A merged event across haplotypes: a representative event + its carriers.
struct MergedRecord {
    // MODULE: each carrier allele's copies, as edits of the reference (split_copies), and the record's
    // distinct forms in the order written (CN_FORMS).
    std::unordered_map<std::size_t, std::vector<std::size_t>> allele_forms;
    std::vector<CopyForm> forms;
    Event seed;                                             // representative (most-carried member)
    // The largest member event's size, which the --min-sv-bp filter judges a record by, so the choice of
    // representative does not change which records exist.
    std::size_t filter_bp = 0;
    std::vector<std::string> carriers;                      // haplotype path names with GT=1
    std::unordered_set<std::size_t> member_alleles;         // allele indices already merged in
    std::set<std::string> member_nodes;                     // union of every merged member event's
                                                            // nodes (describe handoff; ordered set)
    std::unordered_map<std::string, std::size_t> sample_cn; // DUP per-sample copy number
    // UNIT only: how many of those copies run opposite to the reference's copy (FORMAT/CNINV).
    std::unordered_map<std::string, std::size_t> sample_cninv;
    // MODULE: each path's bp change over the module's own nodes against the reference's, what the copy-change
    // step weighs. The bubble's whole interior, which the other routes measure over, can be far more than
    // the module.
    std::unordered_map<std::string, long long> sample_bp;
    // MODULE: each path's bp change over its copies, first to last module node with everything walked between
    // (FORMAT:CNBP). Sequence between the copies' shared parts is walked once and so is no module node, yet it
    // goes with a copy gained or lost.
    std::unordered_map<std::string, long long> sample_copies_bp;
    // MODULE: the nodes CN counts passes through, and every node of the module's copies (repeats inside them too).
    std::unordered_set<std::string> anchor_nodes, module_nodes;
    std::size_t min_size_bp = 0;                            // smallest merged member size (SVLEN_RANGE)
    std::size_t max_size_bp = 0;                            // largest merged member size
    std::unordered_map<std::size_t, std::size_t> member_size;   // member allele -> its own event's size
    std::unordered_map<std::size_t, long long> member_net;      // Sub: member allele -> its own net bp change
    // Member allele -> the events of its own that joined this record: what each carrier was called for,
    // before merging chose one representative. benchmark applies these (carrier_events.tsv). Usually one;
    // an allele can bring two events into one merged component.
    std::unordered_map<std::size_t, std::vector<Event>> member_event;
};

// A string INFO value as VCF requires: the characters that separate fields or values percent-encoded.
std::string vcf_escape(const std::string& s) {
    static const char* hex = "0123456789ABCDEF";
    std::string out;
    for (const char c : s) {
        if (c == ';' || c == ',' || c == '=' || c == ':' || c == '%' || c == ' ' || c == '\t' || c == '\n' ||
            c == '\r') {
            out += '%';
            out += hex[(static_cast<unsigned char>(c) >> 4) & 0xF];
            out += hex[static_cast<unsigned char>(c) & 0xF];
        } else {
            out += c;
        }
    }
    return out;
}

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

    // Refuse a graph this caller cannot describe truthfully, before any of it is read. `spell()` skips
    // a step whose node is absent, so an incomplete graph does not fail -- it produces shorter
    // sequences, and every coordinate, SVLEN and allele derived from them is confidently wrong. Call
    // spells by concatenating whole nodes, so a non-zero overlap would be double-counted the same way.
    validate_graph_paths(graph, "call", /*require_sequences=*/true, /*require_zero_overlaps=*/true);

    // Exact name wins, else a unique substring as written, else a unique case-insensitive one; the rule is
    // shared with bubble and inspect.
    const std::string ref_name = resolve_reference_path_name(graph, options.reference_path, "call");
    const PathRecord* ref_path = nullptr;
    for (const PathRecord& p : graph.paths) {
        if (p.name == ref_name) { ref_path = &p; break; }
    }
    if (ref_path == nullptr) throw std::runtime_error("call: reference path not found: " + ref_name);

    const std::vector<Bubble> bubbles = read_bubbles_csv(options.bubbles_csv_in, /*with_inversion_sites=*/true);
    // A locus with no sites (bubble and panphorte write a header-only CSV for one) has nothing to call:
    // its VCF is the header alone, not an error, and the run's 'sites: 0' line says so.
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

    const std::vector<std::size_t> ref_prefix = path_prefix_bp(*ref_path, graph.nodes);
    // The name's interval convention follows from the path's length (see parse_reference_path_label).
    const ParsedReferencePath ref_meta =
        parse_reference_path_label(ref_name, ref_prefix.empty() ? 0 : ref_prefix.back(), ref_path->type == 'W');
    if (ref_meta.length_mismatch)
        std::cerr << cli::log_prefix("call") << "warning: the interval in reference path name '" << ref_name << "' does not match the path's "
                  << (ref_prefix.empty() ? 0 : ref_prefix.back()) << " bp; reading its start as "
                  << (ref_meta.half_open ? "0-based (half-open)" : "the W line's start, plus one") << "\n";
    if (!ref_meta.has_interval)
        std::cerr << cli::log_prefix("call") << "warning: reference path name '" << ref_name << "' has no :start-end interval; "
                  << "positions are counted from the path's first base (1) on CHROM " << ref_meta.chrom << "\n";
    // The reference spelled in genomic orientation, which alleles are aligned against and events compared on.
    std::string ref_genome = spell_path_steps_sequence(graph, ref_path->steps);
    for (char& c : ref_genome) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));

    // Last reference coordinate: no record's END may exceed it (the reference is finite). Guards
    // against a copy-number event whose bp estimate runs past the region -- a low-complexity tangle
    // can make a module DUP's span balloon past the graph.
    std::size_t ref_total_bp = 0;
    for (const auto& st : ref_path->steps) ref_total_bp += node_len(graph, st.node_id);
    const std::size_t ref_end_1based =
        ref_meta.region_start_1based + (ref_total_bp > 0 ? ref_total_bp - 1 : 0);
    // A module spanning more than 80% of the reference is a tangle summing scattered revisits, not a
    // duplication: real modules span a few percent of a locus.
    const std::size_t max_dup_bp = static_cast<std::size_t>(0.8 * static_cast<double>(ref_total_bp));

    // Reference base at a 1-based genomic coordinate, by binary search over the reference walk's
    // cumulative bp. Coordinate-driven rather than node-driven: every node->position map here records a
    // node's first occurrence, which is the wrong one where the reference revisits the node. Returns an
    // empty string when the coordinate is outside the region.
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
            std::cerr << cli::log_prefix("call") << "warning: --gtf given but reference path '" << ref_name
                      << "' is not PanSN (sample#hap#contig:start-end); skipping gene annotation\n";
        } else {
            const std::size_t lo = ref_meta.region_start_1based;
            const std::size_t hi = lo + (ref_prefix.empty() ? 0 : ref_prefix.back()) - 1;
            genes = parse_gtf(options.gtf_path, ref_meta.chrom, lo, hi, "call");
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
        names.erase(std::unique(names.begin(), names.end()), names.end());   // two genes may share a name
        return names;
    };
    // Sorted, unique gene names whose span overlaps reference bases lo..hi (1-based, inclusive).
    auto genes_in_interval = [&](std::size_t lo, std::size_t hi) {
        std::vector<std::string> names;
        for (const GeneFeature& g : genes)
            if (g.start_1based <= hi && g.end_1based >= lo) names.push_back(g.gene_name);
        std::sort(names.begin(), names.end());
        names.erase(std::unique(names.begin(), names.end()), names.end());
        return names;
    };

    // Fixed sample list: every haplotype path (the reference included).
    std::vector<std::string> sample_names;
    sample_names.reserve(graph.paths.size());
    for (const PathRecord& p : graph.paths) sample_names.push_back(p.name);

    // --diploid output columns: one per PanSN sample (text before the first '#'), with a slot per
    // haplotype number (the second field) in increasing order; a path without '#' is a sample of its
    // own. A haplotype assembled in several contigs is one slot, its paths merged per record. The
    // reference path is not a column, as in vg deconstruct.
    // A sample's column: its haplotype slots, each the paths of one haplotype. A numeric PanSN haplotype
    // field names a slot; a path whose field is not a number (such as 'U' for unassigned contigs) is a
    // haplotype of its own, unphased with the others, so the column's GT is then joined by '/'.
    struct OutColumn { std::string name; std::vector<std::vector<std::size_t>> slots; bool phased = true; };
    std::vector<OutColumn> out_columns;
    std::size_t ref_column = SIZE_MAX;
    for (std::size_t i = 0; i < sample_names.size(); ++i) if (sample_names[i] == ref_name) ref_column = i;
    if (options.diploid) {
        std::map<std::string, std::size_t> col_of;
        std::vector<std::map<long long, std::vector<std::size_t>>> haps;
        std::vector<std::vector<std::size_t>> unnumbered;
        for (std::size_t i = 0; i < sample_names.size(); ++i) {
            if (i == ref_column) continue;
            const std::string& n = sample_names[i];
            const std::size_t a = n.find('#');
            const std::string sample = a == std::string::npos ? n : n.substr(0, a);
            long long hap = -1;
            if (a != std::string::npos) {
                const std::size_t b = n.find('#', a + 1);
                const std::string h = n.substr(a + 1, b == std::string::npos ? std::string::npos : b - a - 1);
                if (!h.empty() && h.size() < 18 && h.find_first_not_of("0123456789") == std::string::npos) hap = std::stoll(h);
            } else {
                hap = 0;   // a plain path name: the sample's one haplotype
            }
            auto [it, fresh] = col_of.emplace(sample, out_columns.size());
            if (fresh) { out_columns.push_back({sample, {}}); haps.emplace_back(); unnumbered.emplace_back(); }
            if (hap >= 0) haps[it->second][hap].push_back(i);
            else unnumbered[it->second].push_back(i);
        }
        for (std::size_t c = 0; c < out_columns.size(); ++c) {
            for (auto& [h, idx] : haps[c]) out_columns[c].slots.push_back(idx);
            for (const std::size_t i : unnumbered[c]) out_columns[c].slots.push_back({i});
            if (!unnumbered[c].empty() && out_columns[c].slots.size() > 1) out_columns[c].phased = false;
        }
    }
    // --diploid: each path's haplotype slot, numbered over all columns; AN and AC count slots.
    std::unordered_map<std::string, std::size_t> slot_of_path;
    {
        std::size_t k = 0;
        for (const OutColumn& col : out_columns)
            for (const std::vector<std::size_t>& slot : col.slots) {
                for (const std::size_t i : slot) slot_of_path.emplace(sample_names[i], k);
                ++k;
            }
    }
    // A record line with per-path columns -> per-sample columns. GT: the slots joined by '|' ('/' when the
    // column is unphased). CN and the other per-haplotype amounts: the sample's sum, '.' if a slot is
    // unknown; HCN keeps the slots' CN. AN loses the reference path's allele, and AF follows.
    auto to_diploid = [&](const std::string& line) -> std::string {
        std::vector<std::string> f;
        {
            std::size_t a = 0;
            const std::size_t n = line.size() - (line.empty() || line.back() != '\n' ? 0 : 1);
            while (a <= n) {
                const std::size_t b = line.find('\t', a);
                const std::size_t e = (b == std::string::npos || b > n) ? n : b;
                f.push_back(line.substr(a, e - a));
                a = e + 1;
            }
        }
        if (f.size() < 10) return line;
        auto split = [](const std::string& x, char d) {
            std::vector<std::string> out; std::size_t a = 0;
            for (;;) { const std::size_t b = x.find(d, a); out.push_back(x.substr(a, b - a)); if (b == std::string::npos) break; a = b + 1; }
            return out;
        };
        const std::vector<std::string> keys = split(f[8], ':');
        std::vector<std::vector<std::string>> v(f.size() - 9);
        for (std::size_t i = 0; i < v.size(); ++i) { v[i] = split(f[9 + i], ':'); v[i].resize(keys.size(), "."); }
        std::string info = f[7];
        std::size_t cn_k = keys.size();
        for (std::size_t k = 0; k < keys.size(); ++k) if (keys[k] == "CN") cn_k = k;
        // HCN, each haplotype's copies, on a CNV record only: elsewhere CN is '.'
        const bool hcn = cn_k < keys.size() && std::find(keys.begin(), keys.end(), "CNBP") != keys.end();
        std::size_t gt_k = keys.size();
        for (std::size_t k = 0; k < keys.size(); ++k) if (keys[k] == "GT") gt_k = k;
        std::string cols;
        long long an = 0, ac = 0;   // over the alleles written, so AN/AC/AF match the GT columns
        std::vector<long long> per_allele;   // count per ALT index (1-based; index 0 unused)
        for (const OutColumn& col : out_columns) {
            // Each slot's values: those of the path of the haplotype that crosses the site. Where two contigs of
            // one haplotype both cross it (overlapping contigs), they are the haplotype there when they agree.
            // When they differ, a contig carrying the record makes the haplotype a carrier, with that contig's
            // values, as the per-path AC counts it; otherwise which one it is is unknown: '.'.
            std::vector<std::vector<std::string>> sv(col.slots.size(), std::vector<std::string>(keys.size(), "."));
            for (std::size_t s = 0; s < col.slots.size(); ++s) {
                std::size_t known = 0, which = 0, carrying = SIZE_MAX;
                bool agree = true;
                for (const std::size_t i : col.slots[s]) {
                    const bool k_any = gt_k < keys.size() ? v[i][gt_k] != "." : v[i][0] != ".";
                    if (!k_any) continue;
                    if (known > 0 && v[i] != v[which]) agree = false;
                    if (carrying == SIZE_MAX && gt_k < keys.size() && v[i][gt_k] != "0") carrying = i;
                    ++known; which = i;
                }
                if (known >= 1 && agree) sv[s] = v[which];
                else if (carrying != SIZE_MAX) sv[s] = v[carrying];
                if (gt_k < keys.size() && sv[s][gt_k] != ".") {
                    ++an;
                    if (sv[s][gt_k] != "0") ++ac;
                    const std::string& g = sv[s][gt_k];
                    if (!g.empty() && g.find_first_not_of("0123456789") == std::string::npos && g.size() < 9) {
                        const std::size_t idx = static_cast<std::size_t>(std::stoul(g));
                        if (idx > 0) { if (per_allele.size() <= idx) per_allele.resize(idx + 1, 0); ++per_allele[idx]; }
                    }
                }
            }
            std::string& out = cols;
            out += '\t';
            for (std::size_t k = 0; k < keys.size(); ++k) {
                if (k) out += ':';
                const std::string& key = keys[k];
                if (col.slots.size() == 1) { out += sv[0][k]; continue; }
                if (key == "GT" || key == "CNFORM") {   // per haplotype, joined as GT is
                    for (std::size_t s = 0; s < sv.size(); ++s) { if (s) out += col.phased ? '|' : '/'; out += sv[s][k]; }
                    continue;
                }
                bool known = true;
                for (const auto& x : sv) if (x[k] == ".") known = false;
                if (!known) { out += '.'; continue; }
                long long t = 0;
                for (const auto& x : sv) t += std::stoll(x[k]);
                out += std::to_string(t);
            }
            if (hcn) {
                out += ':';
                for (std::size_t s = 0; s < sv.size(); ++s) { if (s) out += col.phased ? '|' : '/'; out += sv[s][cn_k]; }
            }
        }
        // AN/AC/AF from the alleles written (the reference path is not a column). A single AC counts the
        // carriers; a list (a multiallelic record) counts each ALT index, as many as it lists.
        if (gt_k < keys.size()) {
            std::vector<std::string> kv = split(info, ';');
            std::size_t n_alt = 0;
            for (const std::string& t : kv)
                if (t.rfind("AC=", 0) == 0) n_alt = static_cast<std::size_t>(std::count(t.begin(), t.end(), ',')) + 1;
            auto fmt_af = [&](long long c) {
                char buf[32]; std::snprintf(buf, sizeof(buf), "%.4g", an > 0 ? static_cast<double>(c) / an : 0.0);
                return std::string(buf);
            };
            for (std::string& t : kv) {
                if (t.rfind("AN=", 0) == 0) { t = "AN=" + std::to_string(an); continue; }
                const bool is_ac = t.rfind("AC=", 0) == 0, is_af = t.rfind("AF=", 0) == 0;
                if (!is_ac && !is_af) continue;
                std::string v;
                if (n_alt <= 1) v = is_ac ? std::to_string(ac) : fmt_af(ac);
                else
                    for (std::size_t i = 1; i <= n_alt; ++i) {
                        const long long c = i < per_allele.size() ? per_allele[i] : 0;
                        if (i > 1) v += ',';
                        v += is_ac ? std::to_string(c) : fmt_af(c);
                    }
                t = (is_ac ? "AC=" : "AF=") + v;
            }
            info.clear();
            for (std::size_t k = 0; k < kv.size(); ++k) { if (k) info += ';'; info += kv[k]; }
            // No sample column carries any of its alternates: every carrier's slot is unknown (contigs that
            // disagree). The record says nothing about the samples written, so it is not written.
            if (ac == 0) return std::string();
        }
        std::string out;
        for (std::size_t k = 0; k < 7; ++k) { out += f[k]; out += '\t'; }
        out += info; out += '\t'; out += f[8];
        if (hcn) out += ":HCN";
        out += cols;
        out += '\n';
        return out;
    };
    // Per-path columns: AN/AC/AF over haplotypes, as --diploid counts slots. The reference path leaves AN, and a
    // haplotype assembled in several contigs (PanSN sample#haplotype, numbered) counts once: a carrier where one
    // of its contigs carries the record, so AC is the haplotypes with a carrying column. A path with no
    // numbered haplotype is one of its own.
    std::vector<std::size_t> path_slot(sample_names.size(), SIZE_MAX);
    {
        std::map<std::string, std::size_t> slot_of;
        for (std::size_t i = 0; i < sample_names.size(); ++i) {
            if (i == ref_column) continue;
            const std::string& n = sample_names[i];
            const std::size_t a = n.find('#');
            const std::size_t b = a == std::string::npos ? std::string::npos : n.find('#', a + 1);
            const std::string h = a == std::string::npos ? std::string() : n.substr(a + 1, b == std::string::npos ? std::string::npos : b - a - 1);
            const bool numbered = !h.empty() && h.size() < 18 && h.find_first_not_of("0123456789") == std::string::npos;
            const std::string key = numbered ? n.substr(0, b) : "\x01" + n;
            path_slot[i] = slot_of.emplace(key, slot_of.size()).first->second;
        }
    }
    auto per_path_counts = [&](const std::string& line) -> std::string {
        std::vector<std::string> f;
        {
            std::size_t a = 0;
            const std::size_t n = line.size() - (line.empty() || line.back() != '\n' ? 0 : 1);
            while (a <= n) {
                const std::size_t b = line.find('\t', a);
                const std::size_t e = (b == std::string::npos || b > n) ? n : b;
                f.push_back(line.substr(a, e - a));
                a = e + 1;
            }
        }
        if (f.size() < 10 || f.size() - 9 != sample_names.size()) return line;
        std::unordered_map<std::size_t, std::string> slot_gt;   // slot -> its GT: a carrying contig's, else 0
        for (std::size_t i = 0; i < sample_names.size(); ++i) {
            if (path_slot[i] == SIZE_MAX) continue;
            const std::string& c = f[9 + i];
            const std::string gt = c.substr(0, c.find(':'));
            if (gt == ".") continue;
            auto [it, fresh] = slot_gt.emplace(path_slot[i], gt);
            if (!fresh && it->second == "0") it->second = gt;   // a carrying contig makes the haplotype a carrier
        }
        long long an = 0, ac = 0;
        for (const auto& [k, g] : slot_gt) {
            ++an;
            if (g != "0") ++ac;
        }
        std::string info, sep;
        std::size_t a = 0;
        const std::string& in = f[7];
        while (a <= in.size()) {
            const std::size_t b = std::min(in.find(';', a), in.size());
            std::string kv = in.substr(a, b - a);
            if (kv.rfind("AN=", 0) == 0) kv = "AN=" + std::to_string(an);
            else if (kv.rfind("AC=", 0) == 0) kv = "AC=" + std::to_string(ac);
            else if (kv.rfind("AF=", 0) == 0) {
                char buf[32]; std::snprintf(buf, sizeof(buf), "%.4g", an > 0 ? static_cast<double>(ac) / an : 0.0);
                kv = std::string("AF=") + buf;
            }
            info += sep + kv; sep = ";";
            a = b + 1;
        }
        f[7] = info;
        std::string out;
        for (std::size_t k = 0; k < f.size(); ++k) { if (k) out += '\t'; out += f[k]; }
        return out + '\n';
    };
    // A record line as it is written: per sample with --diploid, else per path with AN over haplotypes.
    auto written = [&](const std::string& line) { return options.diploid ? to_diploid(line) : per_path_counts(line); };
    auto emit = [&](std::ostream& o, const std::string& line) { o << written(line); };
    // --min-alt-af against the AF a record line is written with, AC/AN (after to_diploid with --diploid).
    const auto af_below = [&](const std::string& line) {
        std::size_t from = 0;   // INFO is the eighth column
        for (int k = 0; k < 7; ++k) {
            const std::size_t t = line.find('\t', from);
            if (t == std::string::npos) return false;
            from = t + 1;
        }
        const std::size_t to = line.find('\t', from);
        const std::string info = line.substr(from, to == std::string::npos ? std::string::npos : to - from);
        long long an = -1, ac = -1;
        std::size_t a = 0;
        while (a <= info.size()) {
            const std::size_t b = std::min(info.find(';', a), info.size());
            const std::string kv = info.substr(a, b - a);
            if (kv.rfind("AN=", 0) == 0) an = std::atoll(kv.c_str() + 3);
            else if (kv.rfind("AC=", 0) == 0) ac = std::atoll(kv.c_str() + 3);
            a = b + 1;
        }
        if (an < 0 || ac < 0) return false;
        return an == 0 || static_cast<double>(ac) / static_cast<double>(an) < options.min_maf;
    };

    // Per-path bubble index, reused across bubbles: the sites' boundaries and interiors only, which is all a walk
    // reads, so the index does not hold every step of the whole panel.
    std::unordered_set<std::string> site_nodes;
    for (const Bubble& b : bubbles) {
        if (b.inversion_parent) continue;   // typed from the whole walks below, not from a site walk
        site_nodes.insert(b.source);
        site_nodes.insert(b.sink);
        site_nodes.insert(b.inside.begin(), b.inside.end());
    }
    std::vector<BubblePathIndex> path_indexes;
    path_indexes.reserve(graph.paths.size());
    for (const PathRecord& p : graph.paths) path_indexes.push_back(build_bubble_path_index(p, site_nodes));
    // name -> path index, for the per-sample FORMAT:CNBP walk-length lookup on DUP records.
    std::unordered_map<std::string, std::size_t> name_to_pi;
    for (std::size_t pi = 0; pi < graph.paths.size(); ++pi) name_to_pi.emplace(graph.paths[pi].name, pi);

    const std::unordered_set<std::string> selfloops = self_loop_nodes(graph);
    // Assembly-gap nodes: a run of 10 or more N, the threshold rebuild uses for its gap pieces. A
    // haplotype walking one at a site has an unknown allele there, not an allele of N bases.
    std::unordered_set<std::string> gap_nodes;
    for (const auto& [id, node] : graph.nodes) {
        std::size_t run = 0;
        for (const char c : node.sequence) {
            run = (c == 'N' || c == 'n') ? run + 1 : 0;
            if (run >= 10) { gap_nodes.insert(id); break; }
        }
    }
    // A gap the graph cut into nodes of fewer than 10 N each is found along each walk: the steps whose N, the N
    // ends of partly-N nodes included, make up a run of 10 or more on that walk cross a gap there. Marked per
    // walk step, not per node: another walk through one of those nodes may hold a shorter run, which is no gap.
    std::vector<std::vector<std::size_t>> split_gap_steps(graph.paths.size());   // per path, ascending
    {
        const auto is_n = [](char c) { return c == 'N' || c == 'n'; };
        struct NEnds { bool all = false; std::size_t lead = 0, trail = 0; };   // in the node's stored orientation
        std::unordered_map<std::string, NEnds> n_ends;
        for (const auto& [id, node] : graph.nodes) {
            const std::string& sq = node.sequence;
            NEnds e;
            while (e.lead < sq.size() && is_n(sq[e.lead])) ++e.lead;
            e.all = !sq.empty() && e.lead == sq.size();
            while (!e.all && e.trail < sq.size() && is_n(sq[sq.size() - 1 - e.trail])) ++e.trail;
            if (e.lead > 0 || e.trail > 0) n_ends.emplace(id, e);
        }
        if (!n_ends.empty())
            for (std::size_t pi = 0; pi < graph.paths.size(); ++pi) {
                const PathRecord& p = graph.paths[pi];
                std::vector<std::size_t>& marks = split_gap_steps[pi];
                std::size_t run = 0;
                std::vector<std::size_t> run_steps;   // the steps holding the current run's N
                const auto close_run = [&]() {
                    if (run >= 10) marks.insert(marks.end(), run_steps.begin(), run_steps.end());
                    run_steps.clear();
                    run = 0;
                };
                for (std::size_t k = 0; k < p.steps.size(); ++k) {
                    const PathStep& s = p.steps[k];
                    const std::size_t len = graph.nodes.at(s.node_id).sequence.size();
                    if (len == 0) continue;
                    const auto it = n_ends.find(s.node_id);
                    if (it == n_ends.end()) { close_run(); continue; }
                    if (it->second.all) { run += len; run_steps.push_back(k); continue; }
                    // the run ends in this node: its N on the walk's way in count, those on the way out start the next
                    const std::size_t lead = s.reverse ? it->second.trail : it->second.lead;
                    const std::size_t trail = s.reverse ? it->second.lead : it->second.trail;
                    if (lead > 0) { run += lead; run_steps.push_back(k); }
                    close_run();
                    if (trail > 0) {
                        run = trail;
                        run_steps.push_back(k);
                    }
                }
                close_run();
                marks.erase(std::unique(marks.begin(), marks.end()), marks.end());
            }
    }
    const auto in_split_gap = [&](std::size_t pi, std::size_t step) {
        const std::vector<std::size_t>& m = split_gap_steps[pi];
        return std::binary_search(m.begin(), m.end(), step);
    };
    // For the UNIT twin test, as refine's twin merge reads them: which kinds of self-link each looped node
    // carries ('t' tandem n+ n+, 'r' fold-back n+ n-, 'l' fold-back n- n+), and which loops some path
    // in the whole graph walks.
    std::unordered_map<std::string, std::set<char>> loop_kinds;
    for (const std::string& id : selfloops) {
        const Node& n = graph.nodes.at(id);
        for (const Neighbor& nb : n.end) if (nb.node_id == id) loop_kinds[id].insert(nb.side == 0 ? 't' : 'r');
        for (const Neighbor& nb : n.start) if (nb.node_id == id) loop_kinds[id].insert(nb.side == 1 ? 't' : 'l');
    }
    // Unit graphs (panphorte's default approximate fold): a folded repeat's copies walk several nodes, segments
    // every copy shares and each copy's own differences, so no node loops for the UNIT route to count. A
    // copy passes the unit graph's anchor exactly once, on one of the anchor nodes; counted there, a unit graph
    // is one unit, keyed by its first anchor node, of the table's unit_bp. Its other nodes leave the
    // DEL/INS walks wherever that unit's copies are counted, as a one-node unit graph's visits do.
    std::unordered_map<std::string, std::string> ug_anchor_of;              // anchor node -> unit key
    std::unordered_map<std::string, std::string> ug_member_of;              // any unit graph node -> unit key
    std::unordered_map<std::string, std::vector<std::string>> ug_members;   // unit key -> its nodes
    std::unordered_map<std::string, std::size_t> ug_unit_bp;                // unit key -> unit length
    // An anchor node the sort flipped is walked reversed by a forward copy; which ones were flipped against the
    // key is read below. A copy's orientation is then its anchor step's, corrected.
    std::unordered_map<std::string, bool> ug_flip;
    if (!options.unit_graph_tsv.empty()) {
        std::ifstream in(options.unit_graph_tsv);
        if (!in) throw std::runtime_error("call: cannot read " + options.unit_graph_tsv);
        std::string line;
        std::getline(in, line);
        std::vector<std::string> head;
        {
            std::stringstream hs(line);
            for (std::string f; std::getline(hs, f, '\t');) head.push_back(f);
        }
        const auto col = [&](const char* name) {
            const auto it = std::find(head.begin(), head.end(), name);
            if (it == head.end())
                throw std::runtime_error("call: " + options.unit_graph_tsv + " has no " + name + " column");
            return static_cast<std::size_t>(it - head.begin());
        };
        const std::size_t c_bp = col("unit_bp"), c_anchor = col("anchor_nodes"), c_nodes = col("nodes");
        // Rows of one site and one motif are one repeat where a path runs from one into the other, with less
        // than a copy of other sequence between them: panphorte folds phase-rotated copies of a motif onto
        // separate unit graphs (an unsplit node cannot hold two rotations), and counted apart a haplotype with
        // one copy on each would read as a loss on one and a gain on the other. Two arrays of one motif with
        // other sequence between them stay two repeats, so a gain on one and a loss on the other do not
        // cancel. Optional columns, so a table without them keeps one unit per row.
        const auto opt_col = [&](const char* name) {
            const auto it = std::find(head.begin(), head.end(), name);
            return it == head.end() ? SIZE_MAX : static_cast<std::size_t>(it - head.begin());
        };
        const std::size_t c_site = opt_col("input_bubble_id"), c_motif = opt_col("canonical_motif");
        const std::size_t c_unit = opt_col("unit"), c_n = opt_col("n_nodes");
        // Node ids alone do not tie a table to a graph (sorted graphs reuse small ids), so each row must also
        // agree with it as refine checks: its anchors among its nodes, its node count its list's, its unit_bp
        // its unit's length, and a one-node unit graph spelling its unit on either strand.
        const auto refuse = [&](const std::string& why) {
            throw std::runtime_error("call: " + options.unit_graph_tsv + " does not describe this graph: " + why +
                                     "; give the table panphorte wrote with this graph (--panphorte-prefix)");
        };
        const auto bases = [](std::string s) {
            for (char& c : s) {
                const char u = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
                c = (u == 'A' || u == 'C' || u == 'G' || u == 'T') ? u : 'N';
            }
            return s;
        };
        const auto split = [](const std::string& s) {
            std::vector<std::string> out;
            if (s.empty() || s == ".") return out;
            std::stringstream ss(s);
            for (std::string f; std::getline(ss, f, ';');) out.push_back(f);
            return out;
        };
        struct Row { std::vector<std::string> anchors, nodes; std::size_t unit_bp; std::string repeat, motif; };
        std::vector<Row> rows;
        while (std::getline(in, line)) {
            if (!line.empty() && line.back() == '\r') line.pop_back();
            if (line.empty()) continue;
            std::vector<std::string> f;
            std::stringstream ls(line);
            for (std::string x; std::getline(ls, x, '\t');) f.push_back(x);
            if (f.size() != head.size())
                throw std::runtime_error("call: malformed row in " + options.unit_graph_tsv + ": " + line);
            Row r{split(f[c_anchor]), split(f[c_nodes]), static_cast<std::size_t>(std::stoull(f[c_bp])), std::string(), std::string()};
            for (const std::string& n : r.nodes)
                if (!graph.nodes.count(n))
                    throw std::runtime_error("call: " + options.unit_graph_tsv + " names node " + n +
                                             ", which the GFA does not contain; the table and the graph are not "
                                             "the same run's");
            const std::string row_name = "the row of " + (r.nodes.empty() ? std::string("no node") : "node " + r.nodes.front());
            for (const std::string& a : r.anchors)
                if (std::find(r.nodes.begin(), r.nodes.end(), a) == r.nodes.end())
                    refuse("anchor node " + a + " of " + row_name + " is not among its nodes");
            if (c_n != SIZE_MAX && f[c_n] != std::to_string(r.nodes.size()))
                refuse(row_name + " lists " + std::to_string(r.nodes.size()) + " nodes, its n_nodes says " + f[c_n]);
            if (c_unit != SIZE_MAX) {
                const std::string unit = bases(f[c_unit]);
                if (r.unit_bp != unit.size())
                    refuse(row_name + " has a " + std::to_string(unit.size()) + " bp unit, its unit_bp says " + f[c_bp]);
                if (r.nodes.size() == 1) {
                    const std::string node = bases(graph.nodes.at(r.nodes.front()).sequence);
                    if (node != unit && node != reverse_complement(unit))
                        refuse("unit node " + r.nodes.front() + " (" + std::to_string(node.size()) +
                               " bp) does not spell its unit (" + std::to_string(unit.size()) + " bp), on either strand");
                }
            }
            if (r.anchors.empty()) continue;   // no position every copy passes: left to the other routes
            if (c_site != SIZE_MAX && c_motif != SIZE_MAX && !f[c_motif].empty()) {
                r.repeat = f[c_site] + "\t" + f[c_motif];
                r.motif = bases(f[c_motif]);
            }
            rows.push_back(std::move(r));
        }
        std::unordered_map<std::string, std::size_t> row_of;   // member node -> row
        for (std::size_t i = 0; i < rows.size(); ++i)
            for (const std::string& n : rows[i].nodes) row_of.emplace(n, i);
        std::vector<std::size_t> parent(rows.size());
        for (std::size_t i = 0; i < rows.size(); ++i) parent[i] = i;
        const std::function<std::size_t(std::size_t)> root = [&](std::size_t i) {
            return parent[i] == i ? i : parent[i] = root(parent[i]);
        };
        // One into the other: less than a copy of other sequence between them (the phase offset of a rotation).
        for (const PathRecord& p : graph.paths) {
            std::size_t last = SIZE_MAX, gap = 0;
            for (const PathStep& st : p.steps) {
                const auto b = row_of.find(st.node_id);
                if (b == row_of.end()) { gap += node_len(graph, st.node_id); continue; }
                const std::size_t a = last;
                last = b->second;
                const std::size_t was = gap;
                gap = 0;
                if (a == SIZE_MAX || a == b->second || rows[a].repeat.empty() || rows[a].repeat != rows[b->second].repeat ||
                    was >= std::max(rows[a].unit_bp, rows[b->second].unit_bp))
                    continue;
                const std::size_t ra = root(a), rb = root(b->second);
                if (ra != rb) parent[std::max(ra, rb)] = std::min(ra, rb);
            }
        }
        for (std::size_t i = 0; i < rows.size(); ++i) {
            const Row& r = rows[i];
            const std::string key = rows[root(i)].anchors.front();
            for (const std::string& a : r.anchors) ug_anchor_of[a] = key;
            for (const std::string& n : r.nodes) ug_member_of[n] = key;
            std::vector<std::string>& members = ug_members[key];
            members.insert(members.end(), r.nodes.begin(), r.nodes.end());
            ug_unit_bp.emplace(key, r.unit_bp);
        }
        // Read on the site walks, source to sink, as copies are counted: a haplotype stored on the other strand
        // walks every node reversed, which says nothing about how the node's sequence was stored.
        std::unordered_map<std::string, std::pair<std::size_t, std::size_t>> fwd_rev;   // anchor -> (+, -) steps
        for (const Bubble& b : bubbles) {
            if (b.inversion_parent ||
                std::none_of(b.inside.begin(), b.inside.end(), [&](const std::string& id) { return ug_anchor_of.count(id) != 0; }))
                continue;
            for (std::size_t pi = 0; pi < graph.paths.size(); ++pi) {
                const auto st = bubble_steps(graph.paths[pi], path_indexes[pi], b, nullptr, nullptr, false);
                if (!st.has_value()) continue;
                for (const PathStep& s : *st)
                    if (ug_anchor_of.count(s.node_id)) ++(s.reverse ? fwd_rev[s.node_id].second : fwd_rev[s.node_id].first);
            }
        }
        const auto usual_rev = [&](const std::string& id) {
            const auto it = fwd_rev.find(id);
            return it != fwd_rev.end() && it->second.second > it->second.first;
        };
        // A row's other anchors are alternatives for the same stretch of the unit as its first, so their
        // sequence says how they were stored: clearly closer to the first's reverse complement is flipped. The
        // walks cannot say it, as a haplotype whose only copy on one of them is inverted walks it reversed. The
        // sequence decides only when both anchors are 16 bp or more and one strand is within a quarter of the
        // longer's length while the other is a quarter further: a short anchor (a 1 bp C against a G) matches
        // the other strand by chance. Rows joined into one repeat hold other stretches of the motif, so a row's
        // first anchor is read against the motif repeated, by the same rule. Anchors the sequence leaves
        // unclear take the walks' usual orientation.
        // Which strand of the motif a node's sequence lies on: +1, -1, or 0 when it does not say clearly.
        const auto motif_strand = [&](const std::string& motif, const std::string& id) {
            const std::string& sq = graph.nodes.at(id).sequence;
            if (motif.empty() || sq.size() < 16) return 0;
            std::string text;
            while (text.size() < sq.size() + motif.size()) text += motif;
            const std::size_t along = infix_distance(sq, text, sq.size());
            const std::size_t against = infix_distance(sq, reverse_complement(text), sq.size());
            const std::size_t near = std::min(along, against), far = std::max(along, against);
            if (near > sq.size() || 4 * near > sq.size() || (far <= sq.size() && 4 * (far - near) < sq.size())) return 0;
            return along < against ? 1 : -1;
        };
        for (std::size_t i = 0; i < rows.size(); ++i) {
            const std::string& key = rows[root(i)].anchors.front();
            const std::string& first = rows[i].anchors.front();
            bool first_flip = usual_rev(first) != usual_rev(key);
            if (first != key) {
                const int o_first = motif_strand(rows[i].motif, first), o_key = motif_strand(rows[i].motif, key);
                if (o_first != 0 && o_key != 0) first_flip = o_first != o_key;
            }
            ug_flip[first] = first_flip;
            const std::string& fs = graph.nodes.at(first).sequence;
            const std::string fs_rc = reverse_complement(fs);
            for (std::size_t k = 1; k < rows[i].anchors.size(); ++k) {
                const std::string& a = rows[i].anchors[k];
                const std::string& as = graph.nodes.at(a).sequence;
                bool rel = usual_rev(a) != usual_rev(first);
                if (std::min(as.size(), fs.size()) >= 16) {
                    const std::size_t len = std::max(as.size(), fs.size());
                    const std::size_t same = nw_distance(as, fs), opposite = nw_distance(as, fs_rc);
                    const std::size_t near = std::min(same, opposite), far = std::max(same, opposite);
                    if (4 * near <= len && 4 * (far - near) >= len) rel = opposite < same;
                }
                ug_flip[a] = rel != first_flip;
            }
        }
    }
    // A step's orientation as its copy reads: reversed for an anchor node walked the other way from its key.
    const auto copy_rev = [&](const PathStep& s) {
        const auto it = ug_flip.find(s.node_id);
        return it != ug_flip.end() && it->second ? !s.reverse : s.reverse;
    };
    std::unordered_set<std::string> walked_loops;
    for (const PathRecord& p : graph.paths)
        for (std::size_t k = 1; k < p.steps.size(); ++k)
            if (p.steps[k].node_id == p.steps[k - 1].node_id) walked_loops.insert(p.steps[k].node_id);

    // Node degree for the tangle guard, as rebuild's gate counts it: distinct neighbours on each end, the
    // larger of the two, self-loops apart. A hub reached from all over the graph has a high degree on one
    // end; a repeat or paralog unit sits in a chain, and a clean branch spreads over both ends.
    std::unordered_map<std::string, std::size_t> node_degree;
    for (const auto& kv : graph.nodes) {
        std::unordered_set<std::string> left, right;
        for (const Neighbor& n : kv.second.start) if (n.node_id != kv.first) left.insert(n.node_id);
        for (const Neighbor& n : kv.second.end) if (n.node_id != kv.first) right.insert(n.node_id);
        node_degree[kv.first] = std::max(left.size(), right.size());
    }

    cli::ensure_parent_dir_for_file(options.out_prefix + ".region.vcf");

    auto write_vcf_header = [&](std::ostream& out) {
        // the --max-seq-bp this run wrote the sequences with
        const std::string seq_limit = options.max_seq_bp == 0 ? std::string()
                                    : "; left out above " + std::to_string(options.max_seq_bp) + " bp (--max-seq-bp)";
        out << "##fileformat=VCFv4.2\n";
        out << "##source=panvar call\n";
        out << "##reference=" << ref_name << "\n";
        out << "##contig=<ID=" << ref_meta.chrom << ">\n";
        out << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Last reference base the variant spans, inclusive. POS is the base BEFORE the event (symbolic convention), so the event occupies POS+1..END and END-POS is its reference span: |SVLEN| for DEL/INV, 0 for INS (which spans no reference), the replaced bases for a SUBSTITUTION (whose SVLEN is the net change), and for a CNV the reference's copies: on CN_METHOD=MODULE the module's own reference span, the bubble interior (FORMAT:CNBP is measured over each haplotype's copies within it), on UNIT the reference's array\">\n";
        out << "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Structural variant type\">\n";
        out << "##INFO=<ID=SVLEN,Number=A,Type=Integer,Description=\"ALT-REF length difference for DEL/INS (negative for a deletion). A DEL/INS/INV record merged across haplotypes reports POS, END, SVLEN and its sequence from one representative member: the exact event carried by the most haplotypes (ties: the larger). Members' sizes are in SVLEN_RANGE; the record is kept by --min-sv-bp on its largest member, so a representative can be smaller than that."
            << " For INV, the length of the inverted segment, which changes no length. On a CN_METHOD=MODULE CNV, one copy of the module, whose copies can hold several paralogs (each haplotype's own size is FORMAT/CNBP); absent on a CN_METHOD=UNIT record whose carriers both gain and lose copies: no single record-level size exists there -- read SVLEN_RANGE, or FORMAT:CN and CNBP for the per-sample size. A UNIT record whose carriers all change in one direction gives its largest change\">\n";
        out << "##INFO=<ID=SVLEN_RANGE,Number=2,Type=Integer,Description=\"Min,max length change among the events merged into this record (every event each member brought), signed as SVLEN, when they differ; a substitution counts by its net change. On an INV, the inverted lengths. On a CN_METHOD=UNIT CNV: the smallest and largest carrier change, (CN-REF_CN)*RU_LEN in bp\">\n";
        out << "##INFO=<ID=BUBBLE_ID,Number=1,Type=Integer,Description=\"panvar bubble identifier\">\n";
        out << "##INFO=<ID=EVENT_NODES,Number=.,Type=String,Description=\"Variant node set; absent on PATH_INVERSION records\">\n";
        out << "##INFO=<ID=INS_SUBTYPE,Number=1,Type=String,Description=\"--classify-ins, on every INS record (a replacement flagged SUBSTITUTION too), from the bases it inserts (INSSEQ, also where too long to write): DUP where minimap2 aligns half of them or more to the reference across the site at --ins-dup-min-identity, NOVEL otherwise\">\n";
        out << "##INFO=<ID=INVDUP,Number=0,Type=Flag,Description=\"Inverted duplication: the extra copy is the reverse complement of the sequence it copies. On a CN_METHOD=UNIT CNV, some carrier gains copies walked opposite to the reference's (see FORMAT/CNINV), or, with no reference copy, holds a copy and its reverse complement; on an INS with INS_SUBTYPE=DUP, the inserted sequence aligns to the reference across the site on the reverse strand\">\n";
        out << "##INFO=<ID=REF_CN,Number=1,Type=Integer,Description=\"Reference copy number (CNV): copies of the repeat unit (CN_METHOD=UNIT) or of the module (MODULE)\">\n";
        out << "##INFO=<ID=UNIT_TWIN_NODES,Number=.,Type=String,Description=\"CN_METHOD=UNIT only: nodes refine's twin merge would fold into the loop node: the same sequence and the same oriented neighbours, a looped one only if its loop is walked. Their traversals are counted as copies of the unit in CN, REF_CN and CNINV\">\n";
        out << "##INFO=<ID=RU_LEN,Number=1,Type=Integer,Description=\"Repeat-unit length in bp, one represented copy. Emitted only for CN_METHOD=UNIT; (CN-REF_CN)*RU_LEN is the repeat-copy component of the haplotype size, while FORMAT/CNBP is the full module-walk change and CNRESID is their difference\">\n";
        out << "##INFO=<ID=CN_UNIT_BP,Number=1,Type=Integer,Description=\"CN_METHOD=MODULE only: one copy of the module, the reference's bp over its nodes divided by REF_CN; with REF_CN=0, the median over the haplotypes holding copies of their bp per copy. Copies of different types differ, so a carrier's own change is FORMAT:CNBP\">\n";
        out << "##INFO=<ID=CN_METHOD,Number=1,Type=String,Description=\"How CN was measured: UNIT (copies of a panphorte unit graph, counted exactly where they pass its anchor; a one-node unit graph's are its traversals) or MODULE (--cn: passes through the anchor of a collapsed module, the module nodes the most paths revisit)\">\n";
        out << "##INFO=<ID=CN_FORMS,Number=.,Type=String,Description=\"CN_METHOD=MODULE: the copies carriers gained or lost, each an edit of the reference, START-END:SEQ (reference bases START..END, 1-based, replaced by SEQ in reference orientation; END=START-1 for an insertion, SEQ empty for a deletion, * where its literal bases are over --max-seq-bp). SEQ may name the reference stretches it repeats, =A-B (bases A..B) or ~A-B (their reverse complement), joined with literal bases by +. FORMAT/CNFORM names each carrier's; with its other records they rebuild its walk\">\n";
        out << "##INFO=<ID=CN_MODULE_REF_BP,Number=1,Type=Integer,Description=\"Reference bp across the bubble interior, each visit counted, over which CN_METHOD=UNIT measures CNBP. On a MODULE record END-POS equals it, unless CN_SPAN_AMBIGUOUS (the span then also encloses the boundary visits between the outermost ones); on a UNIT record END covers the reference's copies, which can lie inside the interior\">\n";
        out << "##INFO=<ID=CN_ANCHOR_BP,Number=1,Type=Integer,Description=\"CN_METHOD=MODULE: bp of the module's anchor, the nodes CN counts passes through\">\n";
        out << "##INFO=<ID=CN_MODULE_NODES_BP,Number=1,Type=Integer,Description=\"CN_METHOD=MODULE: bp of the module's nodes, those some path visits more than once, one visit each\">\n";
        out << "##INFO=<ID=SUBSTITUTION,Number=0,Type=Flag,Description=\"The DEL or INS replaces the reference bases POS+1..END with other bases, typed by its net change, which SVLEN gives. Where both sides fit --max-seq-bp and the breakpoints are exact, REF and ALT spell it after the base at POS, as on a SUB; otherwise ALT is symbolic, with DELSEQ and INSSEQ where they fit. A replacement whose net change is under --min-sv-bp is a SVTYPE=SUB record\">\n";
        out << "##INFO=<ID=HOMLEN,Number=1,Type=Integer,Description=\"Exact INS/DEL (not a substitution, not IMPRECISE): how many bases the event can slide right along the reference and stay the same haplotype. POS is left-aligned, so the equivalent placements are POS..POS+HOMLEN\">\n";
        out << "##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Interval around POS where the event can start. A left-aligned exact INS/DEL: 0,HOMLEN, its equivalent placements. An IMPRECISE record: from the start of an N run the base at POS lies in, to how far the removed bases can lie after POS in a cluster (END-POS less |SVLEN|). A PATH_INVERSION: 0 to the length of the sequence between the left flank and the inverted run (typically an inverted repeat)\">\n";
        out << "##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"Interval around END where the event can end. A left-aligned exact DEL: 0,HOMLEN. An IMPRECISE record: from how far the removed bases can lie before END in a cluster, to the end of an N run after END. A PATH_INVERSION: minus the length of the sequence between the inverted run and the right flank, to 0\">\n";
        out << "##INFO=<ID=PATH_INVERSION,Number=0,Type=Flag,Description=\"An inversion at an inversion site (bubble kind 'inversion'), typed from the haplotypes' walks: a run of reference nodes a haplotype walks reversed, between two runs it walks forward. Such a site is made where no snarl holds the inversion -- its breakpoints lie in inverted repeats the graph collapses, so every path crosses the repeat nodes in both orientations -- and encloses the sites inside it. POS is the last forward base before the inverted run, END the last base before the forward run after it; CIPOS/CIEND give the breakpoint ranges and IMPRECISE is set when they are not zero. Typed from the walks rather than from one site's nodes, it has no EVENT_NODES, START_NODE or END_NODE\">\n";
        out << "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"The breakpoints are not exact, for one of three reasons. The record describes an event CLUSTER rather than one exact event: coalescing joined several pieces that keep reference sequence between them, or the deletion spans loop copies a CNV record counts (their bases are left out of SVLEN), so POS+1..END is wider than the |SVLEN| bases removed. A breakpoint abuts an N in the reference (the base before POS+1 or after END), which places the event only to within the N run. Or a PATH_INVERSION whose breakpoints lie in repeats, with CIPOS/CIEND not zero. On a PRECISE DEL or INV, END-POS equals |SVLEN|\">\n";
        out << "##INFO=<ID=CN_SPAN_AMBIGUOUS,Number=1,Type=Integer,Description=\"Path measurements taken over a module span whose source or sink boundary occurs MORE THAN ONCE in that path. The span is then first-source..last-sink, which is a choice: it can enclose two separate visits and the content between them rather than one module. Every quantity measured over the span -- CN, CNBP, CN_MODULE_REF_BP -- inherits that choice on those paths. Absent when every boundary is visited once\">\n";
        out << "##INFO=<ID=MULTI_CROSSING,Number=1,Type=Integer,Description=\"Traversing paths that traverse this bubble more than once, counted over the whole path. A traversal is a run of consecutive steps from one boundary to the other with only bubble-interior nodes between (possibly none); a boundary touched from or towards outside the bubble is not one. The path's walk is one source..sink window, so it either spans the sequence between traversals, and the events and copy numbers taken from it count the interior once per traversal, or holds one traversal and leaves the others untyped. Reported, not corrected. Absent when every path traverses once\">\n";
        out << "##INFO=<ID=MULTI_CROSSING_UNSCORED,Number=1,Type=Integer,Description=\"Of the MULTI_CROSSING paths, those with a traversal the scored walk does not hold whose allele length differs from every reference traversal of this bubble by --min-sv-bp or more: an SV-sized difference no record here describes. Absent when none\">\n";
        if (!genes.empty())
            out << "##INFO=<ID=GENES,Number=.,Type=String,Description=\"Gene(s) overlapping the variant (from --gtf); for an event inside a gene of a copy-number site's copies, that gene's paralogs its carriers hold\">\n";
        out << "##INFO=<ID=MERGE_MIN_SIM,Number=1,Type=Float,Description=\"How alike the events merged into this record are: the weakest similarity of a member event to the representative, in the test that admits it (MERGE_TEST). Under that test's threshold (--merge-seq-identity, --merge-jaccard), the member joined only through other members. Absent where every member holds the representative's event\">\n";
        out << "##INFO=<ID=MERGE_TEST,Number=1,Type=String,Description=\"The test MERGE_MIN_SIM is measured in: spliced (sequence identity of the haplotypes the two events make, for events found by aligning alleles or each one exact change), jaccard (the events' node sets, length-weighted) or seq (their sequences)\">\n";
        out << "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Allele number: haplotypes with an allele at the site (GT not .). The reference path is left out, and a haplotype assembled in several contigs (PanSN sample#haplotype) counts once: a carrier where one of its contigs carries the allele\">\n";
        out << "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count: haplotypes of AN that carry the allele\">\n";
        out << "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency = AC/AN; --min-alt-af is tested on it\">\n";
        out << "##INFO=<ID=INSSEQ,Number=1,Type=String,Description=\"Inserted sequence of the representative member (see SVLEN)" << seq_limit << "\">\n";
        out << "##INFO=<ID=DELSEQ,Number=1,Type=String,Description=\"Deleted reference sequence of the representative member (see SVLEN)" << seq_limit << "\">\n";
        out << "##INFO=<ID=INVSEQ,Number=1,Type=String,Description=\"Inverted reference sequence of the representative member (see SVLEN)" << seq_limit << "\">\n";
        out << "##FILTER=<ID=PASS,Description=\"All filters passed: every record written passed --min-sv-bp, --min-haplotypes and --min-alt-af\">\n";
        out << "##ALT=<ID=DEL,Description=\"Deletion\">\n";
        out << "##ALT=<ID=INS,Description=\"Insertion\">\n";
        out << "##ALT=<ID=INV,Description=\"Inversion\">\n";
        out << "##ALT=<ID=CNV,Description=\"Copy-number variable region (CN_METHOD=MODULE): FORMAT/CN gives each haplotype's copy number. GT is 1 where CN differs from REF_CN and the module's own nodes change by --min-sv-bp or more (FORMAT/MODBP), a gain or a loss, unless the haplotype's change beside its copies goes their way and is more than the copies (a deletion that takes a copy with it), whose record then holds it. A carrier's copies are in CN_FORMS and FORMAT/CNFORM; its other changes are in records of their own\">\n";
        out << "##ALT=<ID=CNV:TR,Description=\"Tandem-repeat copy-number variation: copies of a repeat unit (RU_LEN), counted exactly per haplotype in FORMAT/CN (CN_METHOD=UNIT)\">\n";
        out << (options.diploid
                    ? "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Per-haplotype genotype of the sample's paths, numbered PanSN haplotypes in order joined by '|'; if the sample also has haplotypes with no number (e.g. HPRC 'U'), each of those is a slot too and every slot is joined by '/' (unphased): 1=carrier, 0=reference-like, .=no allele there (as per path), or two contigs of the haplotype cross the site with different values\">\n"
                    : "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"1=carrier, 0=reference-like, .=no allele there: the haplotype does not cross the site, crosses an assembly gap there, could not be aligned there, or its own larger deletion removes the record's bases\">\n");
        out << "##FORMAT=<ID=CN,Number=1,Type=Integer,Description=\"Copy number on a CNV record: copies of the repeat unit (CN_METHOD=UNIT) or of the module (MODULE). A MODULE non-carrier (GT 0) whose change a written record of its own holds reads REF_CN, so a dosage read off CN does not count it twice; any other haplotype reads its measured copies. . on other records\">\n";
        out << "##FORMAT=<ID=CNBP,Number=1,Type=Integer,Description=\"Actual linear bp gained(+)/lost(-) by this haplotype vs the reference, from the spelled walk (node length x traversal multiplicity minus the reference's): over the bubble interior on CN_METHOD=UNIT; on MODULE, from the haplotype's first to its last module node, everything walked between. Recovers the linear SV size that a one-copy unit does not convey; CNV records only.\">\n";
        out << "##FORMAT=<ID=CNRESID,Number=1,Type=Integer,Description=\"CNBP - (CN-REF_CN) copies of RU_LEN (UNIT) or CN_UNIT_BP (MODULE), with the measured CN (a MODULE non-carrier written as REF_CN keeps its own copies here): the bp this haplotype gained or lost that whole copies do NOT explain. 0 means the whole size difference is copy number; otherwise a partial copy or other sequence changed in the same site, and on MODULE also copies of another size than the reference's\">\n";
        out << "##FORMAT=<ID=CNFORM,Number=1,Type=String,Description=\"CN_METHOD=MODULE: the copies this haplotype gained or lost, as 1-based indices into INFO/CN_FORMS, comma-joined; . for a non-carrier, or a carrier whose copies are in no change of its walk\">\n";
        out << "##FORMAT=<ID=MODBP,Number=1,Type=Integer,Description=\"CN_METHOD=MODULE: this haplotype's bp change over the module's own nodes against the reference's, each visit counted; GT 1 needs it to reach --min-sv-bp (see ALT CNV). CNBP measures over the copies with everything walked between them\">\n";
        out << "##FORMAT=<ID=CNINV,Number=1,Type=Integer,Description=\"CN_METHOD=UNIT only: how many of this haplotype's CN copies run opposite to the reference's copy of the unit (with no reference copy, to the orientation most haplotypes' copies take; when both are equally common, the haplotype's copies in its own less common orientation). 0 means every copy is in the reference orientation; a gain with CNINV above the reference's is an inverted duplication (INFO/INVDUP), and with no reference copy a haplotype holding copies in both orientations\">\n";
        if (options.diploid) {
            out << "##FORMAT=<ID=HCN,Number=1,Type=String,Description=\"--diploid, on a CNV record: the copy number of each of the sample's haplotypes, in GT order, joined like GT ('|', or '/' when the sample's GT is unphased); CN is their sum, and CNBP, CNRESID, MODBP and CNINV are sums too\">\n";
            out << "##panvar_columns=diploid: one column per PanSN sample, its haplotype paths joined into a phased GT ('/' when a sample has unnumbered haplotypes); the reference path is not a column and leaves AN\n";
        }
        out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
        if (options.diploid) for (const OutColumn& c : out_columns) out << '\t' << c.name;
        else for (const std::string& s : sample_names) out << '\t' << s;
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
        std::size_t unaligned_missing = 0;      // its haplotypes GT '.' as unaligned, and as within a larger
        std::size_t covered_missing = 0;        // deletion of their own (the run's counts lose a dropped record's)
        // A module record's sample written as REF_CN because records of its own (holders) hold its change: the
        // field it takes instead, with its measured CN, when none of those records is written.
        struct CnHeld { std::size_t col = 0; std::string measured; std::vector<std::string> holders; };
        std::vector<CnHeld> cn_held;
    };
    std::vector<OutRecord> out_records;
    // A record's module samples written as REF_CN whose holders are all left out (`gone`) take their measured CN.
    const auto apply_measured_cn = [](OutRecord& rec, const std::function<bool(const std::string&)>& gone) {
        if (rec.cn_held.empty()) return;
        std::vector<std::string> cols;
        std::size_t a = 0;
        std::string line = rec.line;
        const bool nl = !line.empty() && line.back() == '\n';
        if (nl) line.pop_back();
        while (true) {
            const std::size_t b = line.find('\t', a);
            cols.push_back(line.substr(a, b == std::string::npos ? std::string::npos : b - a));
            if (b == std::string::npos) break;
            a = b + 1;
        }
        bool changed = false;
        for (OutRecord::CnHeld& h : rec.cn_held) {
            bool kept = false;
            for (const std::string& id : h.holders) kept = kept || !gone(id);
            if (kept || h.col >= cols.size() || h.measured.empty()) continue;
            cols[h.col] = h.measured;
            h.measured.clear();
            changed = true;
        }
        if (!changed) return;
        std::string out;
        for (std::size_t k = 0; k < cols.size(); ++k) { if (k) out += '\t'; out += cols[k]; }
        rec.line = out + (nl ? "\n" : "");
    };
    std::vector<std::string> variant_nodes_rows;     // <prefix>.variant_nodes.tsv (describe handoff)
    std::vector<std::string> carrier_events_rows;    // <prefix>.carrier_events.tsv (benchmark handoff)
    std::vector<std::string> dup_gene_cn_rows;       // <prefix>.dup_gene_cn.tsv body (per-gene DUP CN)

    VariantCallSummary summary;

    // A DUP whose per-gene copy number will be resolved from private k-mers after the bubbles are called.
    struct DupGeneTarget {
        std::size_t bubble_id = 0;
        std::string variant_id;
        std::vector<int> gene_idx;   // indices into `genes` of the genes overlapping the bubble
        std::unordered_map<std::string, std::size_t> sample_cn;  // module CN per haplotype (FORMAT:CN)
        std::unordered_set<std::string> module_nodes;            // the nodes of the module's copies
        std::vector<char> in_copies;                             // per gene_idx: lies in the copies
        bool module = false;                                     // a MODULE record (else UNIT)
    };
    std::vector<DupGeneTarget> dup_targets;

    // Per-bubble outputs, each bubble writing its own slot (parallelism is inside a bubble), merged in
    // bubble order so the result is identical regardless of thread count: the region VCF is
    // coordinate-sorted afterwards (total order on unique ids), and variant_nodes keeps bubble order.
    struct BubbleOut {
        VariantCallSummary sum;
        std::vector<OutRecord> records;
        std::vector<std::string> variant_nodes;
        std::vector<std::string> carrier_events;
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
    cli::ProgressBar call_progress(options.quiet ? "" : cli::log_prefix("call") + "calling sites", bubbles.size());
    std::vector<BubbleOut> bouts(bubbles.size());
    // Bubbles run serially; parallelism is inside each bubble, over its alleles (the walk-diff DP is
    // the hot path, and an SV locus is usually one dominant folded bubble). Each bubble writes only to
    // its own slot, so results are order-deterministic.
    auto process_bubble = [&](std::size_t bubble_idx) {
        const Bubble& bubble = bubbles[bubble_idx];
        // Shadow the shared sinks with this bubble's local buffers: the loop body below
        // writes only through these names, so each bubble accumulates independently. id_counts is
        // bubble-scoped (variant ids embed the bubble id), so a local map reproduces the ids exactly.
        BubbleOut& bout = bouts[bubble_idx];
        VariantCallSummary& summary = bout.sum;
        std::vector<OutRecord>& out_records = bout.records;
        std::vector<std::string>& variant_nodes_rows = bout.variant_nodes;
        std::vector<std::string>& carrier_events_rows = bout.carrier_events;
        std::vector<DupGeneTarget>& dup_targets = bout.dup_targets;
        std::unordered_map<std::string, int> id_counts;
        if (!bubble_filter.empty() && bubble_filter.find(bubble.id) == bubble_filter.end()) return;
        if (bubble.inversion_parent) return;   // typed below, from the haplotypes' walks
        ++summary.bubbles_seen;

        // Reference walk + distinct alleles (grouped by canonical-walk signature).
        std::size_t ref_idx = graph.paths.size();
        for (std::size_t pi = 0; pi < graph.paths.size(); ++pi) {
            if (graph.paths[pi].name == ref_name) { ref_idx = pi; break; }
        }
        BubbleAlleleSet allele_set = enumerate_bubble_alleles(graph, path_indexes, bubble, ref_name, &gap_nodes);
        if (!allele_set.has_reference) {
            return; // reference does not traverse this bubble; cannot type events
        }
        // An assembly gap on the reference's own allele: the reference sequence there is unknown, so no
        // allele can be typed against it and the whole site is skipped.
        for (std::size_t k = 1; k + 1 < allele_set.reference_steps.size(); ++k)
            if (gap_nodes.count(allele_set.reference_steps[k].node_id) ||
                in_split_gap(ref_idx, allele_set.reference_step_indices[k])) {
                ++summary.ref_gap_bubbles;
                summary.ref_gap_ids.push_back(bubble.id);
                return;
            }
        // A path whose own walk crosses a split gap inside the site leaves the alleles as a gap-node walker
        // does. The rest keep the order they are first seen in, now counted without the paths left out.
        {
            std::optional<bool> inward;
            bool dropped = false;
            for (BubbleAllele& a : allele_set.alleles) {
                std::vector<std::string> kept;
                for (const std::string& name : a.members) {
                    const std::size_t pi = name_to_pi.at(name);
                    bool gap = false;
                    if (pi != ref_idx && !split_gap_steps[pi].empty()) {
                        if (!inward) inward = bubble_orientation_holds(graph.paths, path_indexes, bubble);
                        std::vector<std::size_t> idx;
                        bubble_steps(graph.paths[pi], path_indexes[pi], bubble, nullptr, &idx, *inward);
                        for (std::size_t k = 1; k + 1 < idx.size() && !gap; ++k) gap = in_split_gap(pi, idx[k]);
                    }
                    if (gap) {
                        ++allele_set.gap_paths;
                        allele_set.traversing.erase(name);
                        dropped = true;
                    } else {
                        kept.push_back(name);
                    }
                }
                a.members = std::move(kept);
            }
            if (dropped) {
                std::vector<BubbleAllele>& al = allele_set.alleles;
                al.erase(std::remove_if(al.begin(), al.end(), [](const BubbleAllele& a) { return a.members.empty(); }),
                         al.end());
                std::stable_sort(al.begin(), al.end(), [&](const BubbleAllele& x, const BubbleAllele& y) {
                    return name_to_pi.at(x.members.front()) < name_to_pi.at(y.members.front());
                });
                for (std::size_t i = 0; i < al.size(); ++i) al[i].id = i;
            }
        }
        summary.gap_missing += allele_set.gap_paths;
        ++summary.bubbles_with_reference;

        const std::vector<PathStep> ref_steps = std::move(allele_set.reference_steps);
        // Which step of the reference path each entry of ref_steps came from. This is what makes a
        // walk position identify an occurrence rather than just a node name.
        const std::vector<std::size_t> ref_step_idx = std::move(allele_set.reference_step_indices);
        // Genomic extent [first base, last base] of a run of reference-walk tokens, by full-path step
        // index rather than by node name. Returns {0,0} when unavailable. Takes min/max over the run
        // because a reverse-oriented bubble walks the reference in decreasing coordinate, so the
        // walk-order first token is the genomically last one.
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
        // The bubble's orientation on the reference, read from the reference's own walk across it. The
        // walk is spelled source->sink and ref_step_idx says where on the reference path each step
        // sits, so the reference crosses sink->source exactly when its source step comes after its
        // sink step; the boundaries' first occurrences would misread a boundary that also occurs
        // elsewhere on the reference. The record anchors, POS and reverse-complemented sequences use
        // this, and module_span reads each path's direction from its own crossing the same way.
        const bool ref_rev = ref_step_idx.size() >= 2 && ref_step_idx.front() > ref_step_idx.back();
        // The bubble's own nodes. A haplotype crossing the bubble twice is read from its first source to
        // its last sink, so its events can take in nodes lying between the two crossings; those belong
        // to no record of this bubble, and the readers of variant_nodes.tsv refuse them.
        std::unordered_set<std::string> bubble_own(bubble.inside.begin(), bubble.inside.end());
        bubble_own.insert(bubble.source);
        bubble_own.insert(bubble.sink);
        // Genomic first base of the boundary occurrence that walk uses (the source is its first step,
        // the sink its last), not the boundary node's first occurrence anywhere. 0 when unavailable.
        auto walk_boundary_start = [&](bool sink) -> std::size_t {
            if (ref_step_idx.empty()) return 0;
            const std::size_t k = sink ? ref_step_idx.back() : ref_step_idx.front();
            if (k + 1 >= ref_prefix.size()) return 0;
            return ref_meta.region_start_1based + ref_prefix[k];
        };
        using Allele = BubbleAllele;
        const std::vector<Allele>& alleles = allele_set.alleles;
        const std::unordered_set<std::string>& traverses = allele_set.traversing;
        const std::string& ref_sig = allele_set.reference_signature;

        // CN nodes: self-loop nodes in this bubble (a one-node unit graph / tandem unit), handled as count-based DUP
        // events. They leave the DEL/INS/INV alignment when a surviving record counts their copies (see
        // walk_drop below); ordinary recurring nodes always stay in it.
        std::unordered_map<std::string, std::size_t> ref_count;
        for (const PathStep& s : ref_steps) ++ref_count[s.node_id];
        std::unordered_set<std::string> cn_nodes;
        for (const auto& [id, c] : ref_count) { (void)c; if (selfloops.count(id)) cn_nodes.insert(id); }
        for (const Allele& a : alleles) {
            for (const PathStep& s : a.steps) if (selfloops.count(s.node_id)) cn_nodes.insert(s.node_id);
        }
        // A unit graph's copies, through its anchor (see ug_anchor_of).
        if (!ug_anchor_of.empty()) {
            for (const PathStep& s : ref_steps)
                if (const auto it = ug_anchor_of.find(s.node_id); it != ug_anchor_of.end()) cn_nodes.insert(it->second);
            for (const Allele& a : alleles)
                for (const PathStep& s : a.steps)
                    if (const auto it = ug_anchor_of.find(s.node_id); it != ug_anchor_of.end()) cn_nodes.insert(it->second);
        }
        // A unit graph's nodes, mapped to its unit where that unit's copies leave the walks.
        const auto ug_units = [&](const std::unordered_set<std::string>& dropped) {
            std::unordered_map<std::string, std::string> m;
            for (const std::string& key : dropped)
                if (const auto it = ug_members.find(key); it != ug_members.end())
                    for (const std::string& n : it->second) m[n] = key;
            return m;
        };
        // The unit's length: a unit graph's from its table (its key is one segment of it), else the loop node's.
        const auto unit_len = [&](const std::string& cn) {
            const auto it = ug_unit_bp.find(cn);
            return it != ug_unit_bp.end() ? it->second : node_len(graph, cn);
        };

        const std::vector<Tok> Rtok = collapse_walk(ref_steps, cn_nodes, ug_units(cn_nodes));

        std::vector<MergedRecord> merged;
        const std::unordered_map<std::string, std::size_t> ref_node_pos = build_ref_node_pos(bubble);

        // The module's oriented step span in one path, shared by every consumer that measures the
        // module (folded-set construction, reference and haplotype module bp, CNBP, CN_MODULE_REF_BP)
        // so the arithmetic cannot drift between them.
        //
        // Deliberately the widest span (first source .. last sink), not the tight allele interval: a
        // module's copies are what lies between the outermost boundaries. `ambiguous` reports a
        // boundary occurring more than once, since the span is then a choice and can sweep in content
        // between two unrelated visits rather than one module.
        // `opens_on_source` receives which boundary sits at the span's low end.
        auto module_span = [&](std::size_t pi, bool* ambiguous = nullptr, bool* opens_on_source = nullptr)
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
            // The direction the path crosses the site in, from its own crossing: a path walking it sink to source
            // opens the span on the sink. Where the boundaries first occur would mislead when a reversed
            // crossing's source also occurs further up the path.
            bool forward = s0 <= k1;
            const std::vector<BubblePathInterval> cr =
                bubble_crossings(graph.paths[pi], path_indexes[pi], bubble, true);
            if (!cr.empty()) forward = cr.front().source_to_sink;
            if (forward && s0 > k1) forward = false;
            if (!forward && k0 > s1) forward = true;
            const std::size_t lo = forward ? s0 : k0;
            const std::size_t hi = forward ? k1 : s1;
            if (opens_on_source != nullptr) *opens_on_source = forward;
            return {lo, hi};
        };

        const std::size_t rescue_floor =
            options.rescue_min_bp != 0 ? options.rescue_min_bp : std::max<std::size_t>(1, options.min_sv_bp / 2);

        // Traversing paths whose module span rests on a repeated boundary. Counted once per path here,
        // since the span resolver runs more than once per path, and at bubble scope so UNIT records
        // carry it too.
        std::size_t module_ambiguous_paths = 0;
        for (std::size_t pi = 0; pi < graph.paths.size(); ++pi) {
            if (!traverses.count(graph.paths[pi].name)) continue;
            bool amb = false;
            (void)module_span(pi, &amb);
            if (amb) ++module_ambiguous_paths;
        }
        // Traversing paths that traverse this bubble more than once. A traversal is a run of
        // consecutive path steps from one boundary to the other with only interior nodes between
        // (none at all counts: S directly followed by K). Counted over the whole path, not only the
        // walk: a boundary merely touched -- reached from, or left towards, outside the bubble, or
        // returned to without reaching the other one -- is no traversal, and a second traversal
        // outside the walk is still one. The walk is a single source..sink window, so it either spans
        // the sequence between traversals (events and copy numbers then count the interior once per
        // traversal) or holds one of them and leaves the others untyped. Reported, not corrected:
        // MULTI_CROSSING on every record of the bubble and a run-summary counter.
        std::size_t multi_crossing_paths = 0;
        // ...and of those, the paths with a traversal the walk does not hold whose allele differs from every
        // reference traversal by --min-sv-bp or more in length: an SV-sized change no record of this bubble
        // can describe. Reported, not called (MULTI_CROSSING_UNSCORED and the run summary).
        std::size_t multi_crossing_unscored = 0;
        if (bubble.source != bubble.sink) {
            const std::unordered_set<std::string> interior(bubble.inside.begin(), bubble.inside.end());
            // Each traversal of path `pi`: its boundary step indices (low, high) and its interior bp.
            struct Traversal { std::size_t lo, hi, bp; };
            auto traversals_of = [&](std::size_t pi) {
                std::vector<Traversal> out;
                const auto& idx = path_indexes[pi].positions;
                const auto sit = idx.find(bubble.source);
                const auto kit = idx.find(bubble.sink);
                if (sit == idx.end() || kit == idx.end()) return out;
                const std::vector<PathStep>& steps = graph.paths[pi].steps;
                // From each boundary visit, skip interior steps; a traversal ends on the other
                // boundary. Each run is scanned once, so the cost is the path's steps in the bubble.
                auto from = [&](const std::vector<std::size_t>& at, const std::string& other) {
                    for (const std::size_t b : at) {
                        std::size_t q = b + 1, bp = 0;
                        while (q < steps.size() && interior.count(steps[q].node_id)) bp += node_len(graph, steps[q++].node_id);
                        if (q < steps.size() && steps[q].node_id == other) out.push_back({b, q, bp});
                    }
                };
                from(sit->second, bubble.sink);
                from(kit->second, bubble.source);
                return out;
            };
            const std::vector<Traversal> ref_trav =
                ref_idx < graph.paths.size() ? traversals_of(ref_idx) : std::vector<Traversal>{};
            for (std::size_t pi = 0; pi < graph.paths.size(); ++pi) {
                if (!traverses.count(graph.paths[pi].name)) continue;
                const std::vector<Traversal> trav = traversals_of(pi);
                if (trav.size() <= 1) continue;
                ++multi_crossing_paths;
                if (ref_trav.empty() || pi == ref_idx) continue;
                BubblePathInterval window;
                if (!bubble_steps(graph.paths[pi], path_indexes[pi], bubble, &window).has_value()) continue;
                for (const Traversal& t : trav) {
                    if (t.lo >= window.left && t.hi <= window.right) continue;   // held by the walk
                    std::size_t nearest = SIZE_MAX;
                    for (const Traversal& r : ref_trav)
                        nearest = std::min(nearest, t.bp > r.bp ? t.bp - r.bp : r.bp - t.bp);
                    if (nearest >= options.min_sv_bp) { ++multi_crossing_unscored; break; }
                }
            }
        }
        if (multi_crossing_unscored > 0) {
            ++summary.multi_crossing_unscored_bubbles;
            summary.multi_crossing_unscored_paths += multi_crossing_unscored;
            summary.multi_crossing_unscored_ids.push_back(bubble.id);
        }
        if (multi_crossing_paths > 0) {
            ++summary.multi_crossing_bubbles;
            summary.multi_crossing_paths += multi_crossing_paths;
        }

        // The reference's module span as genomic coordinates, from the same resolver the folded bp
        // and CNBP are measured over, so with a repeated boundary a record's coordinates and its CN
        // still describe one interval.
        auto module_ref_bounds = [&](bool* ambiguous, bool* opens_on_source = nullptr)
            -> std::pair<std::size_t, std::size_t> {
            if (ref_idx >= graph.paths.size()) return {0, 0};
            const auto sp = module_span(ref_idx, ambiguous, opens_on_source);
            if (sp.first > sp.second || sp.second + 1 >= ref_prefix.size()) return {0, 0};
            return {ref_meta.region_start_1based + ref_prefix[sp.first],
                    ref_meta.region_start_1based + ref_prefix[sp.second + 1] - 1};
        };

        auto ev_ref_pos = [&](const Event& e) -> long long {
            // Occurrence-aware first: this position seeds the merge window and the sort order, and a
            // node's first occurrence would make two events at two visits to the same node look
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
        // An event as the change it makes to the reference: the base before it, the bases it deletes, the
        // sequence it inserts (genomic orientation). Aligned events always have one; node events when exact.
        struct Splice { std::size_t before = 0, del = 0; const std::string* ins = nullptr; };
        auto splice_of = [](const Event& e, Splice& sp) {
            if (e.seq_event) { sp = {e.seq_before, e.del_bp, &e.seq_ins}; return true; }
            if (e.splice_ok) { sp = {e.splice_before, e.splice_del, &e.splice_ins}; return true; }
            return false;
        };
        // The identity of the two haplotypes events a and b make, as edit distance against the larger
        // event; 0 when they cannot be compared or the stretch where they differ is not alike.
        auto spliced_identity = [&](const Event& a, const Event& b) -> double {
            Splice A, B;
            if (!splice_of(a, A) || !splice_of(b, B)) return 0.0;
            if (A.before == B.before && A.del == B.del && *A.ins == *B.ins) return 1.0;
            const std::size_t base = ref_meta.region_start_1based;
            const std::size_t big = std::max(a.size_bp, b.size_bp), small = std::min(a.size_bp, b.size_bp);
            const double len_ratio = options.merge_size_ratio > 0.0 ? options.merge_size_ratio : options.merge_seq_identity;
            if (big == 0 || static_cast<double>(small) < len_ratio * static_cast<double>(big)) return 0.0;
            const std::size_t sa = A.before + 1, sb = B.before + 1;
            const std::size_t ea = sa + A.del, eb = sb + B.del;
            const std::size_t lo = std::min(sa, sb), hi = std::max(ea, eb);
            if (lo < base || hi - base > ref_genome.size()) return 0.0;
            // Each haplotype over [lo, hi): the reference with the event spliced in. Read base by base,
            // so the shared ends are found without building either; only the core between them is
            // built and aligned, and the alignment stops once the identity cannot be reached.
            const std::size_t nx = hi - lo - A.del + A.ins->size(), ny = hi - lo - B.del + B.ins->size();
            auto at = [&](const Splice& e, std::size_t st, std::size_t en, std::size_t i) -> char {
                if (i < st - lo) return ref_genome[lo - base + i];
                i -= st - lo;
                if (i < e.ins->size()) return (*e.ins)[i];
                return ref_genome[en - base + (i - e.ins->size())];
            };
            std::size_t p = 0, q = 0;
            while (p < nx && p < ny && at(A, sa, ea, p) == at(B, sb, eb, p)) ++p;
            while (q < nx - p && q < ny - p && at(A, sa, ea, nx - 1 - q) == at(B, sb, eb, ny - 1 - q)) ++q;
            std::string x, y;
            x.reserve(nx - p - q);
            y.reserve(ny - p - q);
            for (std::size_t i = p; i < nx - q; ++i) x += at(A, sa, ea, i);
            for (std::size_t i = p; i < ny - q; ++i) y += at(B, sb, eb, i);
            // The largest distance that can still reach the identity, with room for rounding: over it is
            // no match, whatever the identity the bound itself would give.
            const auto max_k = static_cast<std::size_t>((1.0 - options.merge_seq_identity) * static_cast<double>(big) + 1e-6);
            const std::size_t ed = nw_distance(x, y, max_k);
            if (ed > max_k) return 0.0;
            // Against the event, two placements of one size in unique sequence look alike: unrelated
            // sequence aligns at about half an edit per base, so two equal large deletions a few kb apart
            // can differ by less than the identity allows. Where the placements are a structural variant's
            // length apart and the haplotypes differ by as many edits, the sequence between them must be the
            // same, as it is where one event sits in either copy of a duplication: the edits are measured
            // against that distance. Base changes clustered in one placement are not a distance apart.
            const std::size_t apart = std::max(sa > sb ? sa - sb : sb - sa, ea > eb ? ea - eb : eb - ea);
            if (apart >= options.min_sv_bp && ed >= options.min_sv_bp &&
                1.0 - static_cast<double>(ed) / static_cast<double>(apart) < options.merge_seq_identity)
                return 0.0;
            return ed >= big ? 0.0 : 1.0 - static_cast<double>(ed) / static_cast<double>(big);
        };
        // Two non-DUP events are the same site: compatible types, anchors within the window, net length
        // changes less than --min-sv-bp apart, and then the haplotypes they make where both are exact
        // changes, else node-set Jaccard or sequence identity. out_jac/out_seq receive the evidence
        // (seq computed only when the Jaccard gate did not already decide it).
        // out_test: the test that decided, "spliced" (the haplotypes the two make), "jaccard" (their nodes) or
        // "seq" (their sequences); empty where none was reached.
        auto events_match = [&](const Event& a, const Event& b,
                                double* out_jac = nullptr, double* out_seq = nullptr, const char** out_test = nullptr) {
            if (out_jac) *out_jac = -1.0;
            if (out_seq) *out_seq = -1.0;
            if (out_test) *out_test = "";
            // Net length change: a deletion's is minus its size, an insertion's its size.
            auto net_of = [](const Event& e) -> long long {
                if (e.type == EvType::Ins) return static_cast<long long>(e.size_bp);
                if (e.type == EvType::Del) return -static_cast<long long>(e.size_bp);
                return static_cast<long long>(e.ins_bp) - static_cast<long long>(e.del_bp);
            };
            auto sgn = [&](const Event& e) { const long long v = net_of(e); return v > 0 ? 1 : v < 0 ? -1 : 0; };
            // Two aligned events are the same when the haplotypes they make are, so a substitution can
            // match a deletion or insertion it nearly is (a run replaced by a few bases against the same
            // run deleted), provided their net changes point the same way; spliced identity and the net gate
            // below decide. Node events keep to one type.
            if (a.type != b.type) {
                const bool sub_pair = a.seq_event && b.seq_event &&
                    (a.type == EvType::Sub || b.type == EvType::Sub) &&
                    a.type != EvType::Inv && b.type != EvType::Inv && a.type != EvType::Dup && b.type != EvType::Dup;
                if (!sub_pair || sgn(a) != sgn(b) || sgn(a) == 0) return false;
            }
            // Substitutions share the nodes of the run they replace whichever way their length moves,
            // so a loss and a gain of one size would pass the node test: their net changes have to
            // point the same way (or both be none).
            if (a.type == EvType::Sub && b.type == EvType::Sub && sgn(a) != sgn(b)) return false;
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
            // Two events found by aligning alleles are the same when the haplotypes they make are: each is
            // spliced into the reference over the stretch both span and the results compared, as edit
            // distance against the larger event. In an exact repeat that is left-normalisation; where the
            // same deletion is placed differently inside a segmental duplication, it measures the
            // duplicate's homology rather than how far the two placements are apart. Their nodes are
            // not used: an aligned event's nodes are those it overlaps, so in a repeat deletions of very
            // different sizes share them all.
            // And their length changes must differ by less than --min-sv-bp, so that a record's carriers
            // never differ in length by a structural variant; the identity alone allows that on large
            // events. Base changes between them are left to the identity, as one large insertion carries
            // ordinary divergence.
            if (a.seq_event && b.seq_event) {
                const long long na = static_cast<long long>(a.ins_bp) - static_cast<long long>(a.del_bp);
                const long long nb = static_cast<long long>(b.ins_bp) - static_cast<long long>(b.del_bp);
                const long long dn = na > nb ? na - nb : nb - na;
                const bool len_ok = dn < static_cast<long long>(options.min_sv_bp);
                if (out_jac) *out_jac = -1.0;
                if (out_test) *out_test = "spliced";
                // The length gate is the cheap test: where it fails and the identity is not asked for, it is
                // not measured.
                if (!len_ok && out_seq == nullptr) return false;
                const double sid = spliced_identity(a, b);
                if (out_seq) *out_seq = sid;
                return sid >= options.merge_seq_identity && len_ok;
            }
            // The node rules keep the aligned rules' length gate: a record's carriers never differ in
            // length by a structural variant, however many nodes they share (a deletion shares most of
            // its bases with a smaller one inside it).
            {
                const long long na = net_of(a), nb = net_of(b);
                if (a.type != EvType::Inv && (na > nb ? na - nb : nb - na) >= static_cast<long long>(options.min_sv_bp))
                    return false;
            }
            // Two events that are each one exact change compare as aligned ones do, by the haplotypes they
            // make, so the same inserted node at two places is two insertions, not one.
            {
                Splice A, B;
                if (a.type != EvType::Inv && splice_of(a, A) && splice_of(b, B)) {
                    const double sid = spliced_identity(a, b);
                    if (out_seq) *out_seq = sid;
                    if (out_test) *out_test = "spliced";
                    return sid >= options.merge_seq_identity;
                }
            }
            const double jac = weighted_jaccard(graph, a.nodes, b.nodes);
            if (out_jac) *out_jac = jac;
            if (out_test) *out_test = "jaccard";
            if (jac >= options.merge_jaccard) return true;
            const double len_ratio =
                options.merge_size_ratio > 0.0 ? options.merge_size_ratio : options.merge_seq_identity;
            // A substitution is compared on both sequences, replaced and replacing: two that replace
            // the same reference run with different sequence are different alleles.
            const double sid = a.type == EvType::Sub
                ? seq_identity(a.seq + "|" + a.alt_seq, b.seq + "|" + b.alt_seq, options.merge_seq_identity, len_ratio)
                : seq_identity(a.seq, b.seq, options.merge_seq_identity, len_ratio);
            if (out_seq) *out_seq = sid;
            if (out_test) *out_test = "seq";
            return sid >= options.merge_seq_identity;
        };

        // A self-loop the UNIT route calls is one whose unit is SV-sized, or whose copy change is: the
        // event is unit x |copy change|, so a short unit gaining many copies is an SV-sized event. Gating
        // on the unit alone would lose such events outright, as the loop is removed from the DEL/INS/INV
        // walks and no other route sees them.
        //
        // A loop inside a counted module's copies yields to the module (see the module route below), so a
        // short element repeated in each module copy does not replace the module's record.
        std::unordered_map<std::string, std::vector<std::size_t>> cn_allele_count;
        // Copies walked opposite to the reference's own copy -- a fold-back (`A+ A-`) is an inverted
        // duplication, which CN alone reports exactly like a tandem one. Where the reference carries no
        // copy, the orientation most haplotype copies take stands in for the reference's.
        std::unordered_map<std::string, std::vector<std::size_t>> cn_allele_inv;
        std::unordered_map<std::string, bool> cn_ref_reverse;
        // Twins of a loop node: the nodes refine's twin merge would fold into it. A graph builder can put
        // one element on two nodes, a looped node for some contigs and a plain twin for the others, and
        // counting the loop node alone then reads a path routed
        // through the twin as a copy loss, and the reverse as a gain. The test is refine's: the same
        // sequence and the same oriented neighbours apart from self-links; a looped member only if some
        // path walks its loop; the kept node is the one with the most kinds of self-link (then the lowest
        // id), and a member joins only if the kept node has every kind it has and it is linked to no
        // member already in. Sharing every neighbour is what makes them one element in one place, so the
        // same sequence elsewhere (a moved copy) is not a twin. Twins count only toward copy number; the
        // DEL/INS/INV walks keep their nodes as they are.
        std::unordered_map<std::string, std::vector<std::string>> rep_twins;   // kept loop node -> its twins
        std::unordered_map<std::string, std::string> twin_of;                  // twin -> its kept loop node
        {
            // Oriented neighbourhood without self-links, as (side of this node, neighbour, neighbour side).
            auto nbr_key = [&](const std::string& id, const Node& n) {
                std::vector<std::string> k;
                for (const Neighbor& nb : n.start) if (nb.node_id != id) k.push_back("s" + nb.node_id + ":" + std::to_string(nb.side));
                for (const Neighbor& nb : n.end) if (nb.node_id != id) k.push_back("e" + nb.node_id + ":" + std::to_string(nb.side));
                std::sort(k.begin(), k.end());
                k.erase(std::unique(k.begin(), k.end()), k.end());
                return k;
            };
            auto id_less = [](const std::string& a, const std::string& b) {
                return a.size() != b.size() ? a.size() < b.size() : a < b;
            };
            std::map<std::pair<std::string, std::vector<std::string>>, std::vector<std::string>> groups;
            for (const std::string& id : bubble.inside) {
                const auto nit = graph.nodes.find(id);
                if (nit == graph.nodes.end() || nit->second.sequence.empty()) continue;
                if (selfloops.count(id) && !walked_loops.count(id)) continue;   // an unwalked loop never joins
                std::vector<std::string> key = nbr_key(id, nit->second);
                if (key.empty()) continue;   // no place to share
                groups[{nit->second.sequence, std::move(key)}].push_back(id);
            }
            static const std::set<char> kNoKinds;
            auto kinds_of = [&](const std::string& id) -> const std::set<char>& {
                const auto it = loop_kinds.find(id);
                return it == loop_kinds.end() ? kNoKinds : it->second;
            };
            for (auto& kv : groups) {
                std::vector<std::string>& ids = kv.second;
                if (ids.size() < 2) continue;
                std::sort(ids.begin(), ids.end(), id_less);
                std::string keep = ids.front();
                for (const std::string& id : ids)
                    if (kinds_of(id).size() > kinds_of(keep).size()) keep = id;
                if (!cn_nodes.count(keep)) continue;   // no walked loop among them: nothing for UNIT to count
                const std::set<char>& kk = kinds_of(keep);
                std::vector<std::string> members{keep};
                for (const std::string& id : ids) {
                    if (id == keep) continue;
                    const std::set<char>& k = kinds_of(id);
                    if (!std::includes(kk.begin(), kk.end(), k.begin(), k.end())) continue;
                    bool independent = true;
                    for (const std::string& m : members)
                        if (graph.has_edge_between(id, m)) { independent = false; break; }
                    if (!independent) continue;
                    members.push_back(id);
                    rep_twins[keep].push_back(id);
                    twin_of[id] = keep;
                }
            }
        }
        // The loop node a step counts toward: itself, or the loop node it is a twin of.
        auto rep_unit_of = [&](const std::string& id) -> const std::string* {
            const auto t = twin_of.find(id);
            if (t != twin_of.end()) return &t->second;
            if (const auto a = ug_anchor_of.find(id); a != ug_anchor_of.end() && cn_nodes.count(a->second))
                return &a->second;
            return cn_nodes.count(id) ? &id : nullptr;
        };
        std::unordered_set<std::string> rep_consumed;   // loops whose copies a surviving UNIT record counts
        bool module_record = false;                     // a MODULE DUP was written for this bubble
        // MODULE records the copy-change step dropped (every carrier's copy change lies in a larger event of
        // its own): no VCF record, but their copy numbers still split the genes in the copies.
        std::vector<MergedRecord> dropped_modules;
        auto same_unit = [&](const std::string& a, const std::string& b) {
            const std::string* ua = rep_unit_of(a);
            const std::string* ub = rep_unit_of(b);
            return (ua != nullptr ? *ua : a) == (ub != nullptr ? *ub : b);
        };
        std::unordered_map<std::string, std::size_t> rep_ref_count;   // reference copies, twins included
        for (const PathStep& s : ref_steps)
            if (const std::string* u = rep_unit_of(s.node_id)) ++rep_ref_count[*u];
        for (const std::string& cn : cn_nodes) {
            cn_allele_count[cn].assign(alleles.size(), 0);
            cn_allele_inv[cn].assign(alleles.size(), 0);
            cn_ref_reverse[cn] = false;
        }
        for (auto it = ref_steps.rbegin(); it != ref_steps.rend(); ++it) {
            const std::string* u = rep_unit_of(it->node_id);
            if (u != nullptr) cn_ref_reverse[*u] = copy_rev(*it);   // ends on the first occurrence
        }
        // A unit the reference does not carry has no reference orientation, and the way its node is stored
        // says nothing about the sequence. Its direct orientation is the one most haplotype copies take.
        // When both are equally common nothing says which is inverted, except within one haplotype: one
        // holding copies in both orientations has its fewer ones inverted, whichever is called direct.
        std::unordered_set<std::string> cn_no_orientation;
        // per unit the reference does not carry, per allele: its copies in its own less common orientation,
        // the copies that reverse-complement another copy of the same haplotype
        std::unordered_map<std::string, std::vector<std::size_t>> cn_allele_minor;
        {
            std::unordered_map<std::string, std::pair<std::size_t, std::size_t>> hap_fwd_rev;
            for (std::size_t ai = 0; ai < alleles.size(); ++ai) {
                std::unordered_map<std::string, std::pair<std::size_t, std::size_t>> own;
                for (const PathStep& s : alleles[ai].steps) {
                    const std::string* u = rep_unit_of(s.node_id);
                    if (u == nullptr || rep_ref_count.count(*u)) continue;
                    auto& fr = own[*u];
                    ++(copy_rev(s) ? fr.second : fr.first);
                }
                for (const auto& [u, fr] : own) {
                    hap_fwd_rev[u].first += fr.first * alleles[ai].members.size();
                    hap_fwd_rev[u].second += fr.second * alleles[ai].members.size();
                    auto& minor = cn_allele_minor[u];
                    minor.resize(alleles.size(), 0);
                    minor[ai] = std::min(fr.first, fr.second);
                }
            }
            for (const auto& [u, fr] : hap_fwd_rev) {
                if (fr.first == fr.second) cn_no_orientation.insert(u);
                else cn_ref_reverse[u] = fr.second > fr.first;
            }
        }
        for (std::size_t ai = 0; ai < alleles.size(); ++ai)
            for (const PathStep& s : alleles[ai].steps) {
                const std::string* u = rep_unit_of(s.node_id);
                if (u == nullptr) continue;
                ++cn_allele_count.at(*u)[ai];
                if (!cn_no_orientation.count(*u) && copy_rev(s) != cn_ref_reverse.at(*u)) ++cn_allele_inv.at(*u)[ai];
            }
        for (const std::string& u : cn_no_orientation) cn_allele_inv.at(u) = cn_allele_minor.at(u);
        auto rep_callable = [&](const std::string& cn) {
            const std::size_t unit = unit_len(cn);
            if (unit >= options.min_sv_bp) return true;
            const std::size_t rc = rep_ref_count.count(cn) ? rep_ref_count.at(cn) : 0;
            for (const std::size_t ac : cn_allele_count.at(cn)) {
                const std::size_t d = ac > rc ? ac - rc : rc - ac;
                if (unit * d >= options.min_sv_bp) return true;
            }
            return false;
        };
        // A tangle: a site whose interior holds ten or more hubs (fifty neighbours or more on one end, the
        // hub rebuild's gate counts) is a low-complexity region, not a copy-number module, and its node
        // multiplicity is no copy number. The module route is not run there. A paralog cluster is a
        // chain of low-degree nodes and has no hubs.
        bool is_tangle = false;
        {
            std::size_t hubs = 0;
            for (const std::string& id : bubble.inside) {
                const auto it = node_degree.find(id);
                if (it != node_degree.end() && it->second >= 50 && ++hubs >= 10) { is_tangle = true; break; }
            }
        }
        // --min-haplotypes / --min-alt-af over a record's carriers, the haplotypes its AC counts, not every
        // haplotype whose CN differs from REF_CN: on UNIT that includes changes too small to make a carrier.
        // --min-alt-af is tested on the AF a record is written with. A copy-number record has no '.'
        // haplotype, so per path its AN is known here, counted as the writer counts it (per_path_counts): the
        // haplotypes crossing the site, the reference path left out and a haplotype's contigs once, a carrier
        // where one of them carries it. Other records lose '.' haplotypes from AN when written, and are tested then.
        // With --diploid a slot whose contigs both cross the site is '.' where their values differ, and some
        // of those values (CNBP) are only measured when the record is written. Contigs that differ in GT, CN
        // or CNINV leave the slot out here; a slot of agreeing non-carrier contigs is left out too, so the AF
        // tested is never below the written one, which the final pass then tests exactly.
        // lowest: AN over every slot a contig of which crosses the site, so the AF tested is never above the
        // written one (a carrier's copies are taken from its events only for a record that cannot then drop).
        auto support_ok = [&](const MergedRecord& mr, bool lowest = false) {
            const std::unordered_set<std::string> uniq(mr.carriers.begin(), mr.carriers.end());
            if (uniq.size() < options.min_haplotypes) return false;
            if (options.min_maf > 0.0 && mr.seed.type == EvType::Dup) {
                std::size_t an = 0, ac = 0;
                if (!options.diploid) {
                    std::unordered_map<std::size_t, int> slot_gt;   // 1 where a contig carries it, else 0
                    for (const std::string& t : traverses) {
                        const auto pit = name_to_pi.find(t);
                        if (pit == name_to_pi.end() || path_slot[pit->second] == SIZE_MAX) continue;
                        int& g = slot_gt[path_slot[pit->second]];
                        if (uniq.count(t) != 0) g = 1;
                    }
                    for (const auto& [slot, g] : slot_gt) {
                        ++an;
                        if (g == 1) ++ac;
                    }
                } else {
                    std::map<std::size_t, std::vector<std::string>> slot_paths;
                    for (const std::string& t : traverses)
                        if (const auto it = slot_of_path.find(t); it != slot_of_path.end()) slot_paths[it->second].push_back(t);
                    const auto value = [](const std::unordered_map<std::string, std::size_t>& m, const std::string& k) {
                        const auto it = m.find(k);
                        return it == m.end() ? SIZE_MAX : it->second;
                    };
                    an = ac = 0;
                    for (const auto& [slot, paths] : slot_paths) {
                        const std::string& p0 = paths.front();
                        bool carrier = false;
                        for (const std::string& p : paths) carrier = carrier || uniq.count(p) != 0;
                        bool agree = true;
                        for (const std::string& p : paths)
                            agree = agree && (uniq.count(p) != 0) == carrier && value(mr.sample_cn, p) == value(mr.sample_cn, p0) &&
                                    value(mr.sample_cninv, p) == value(mr.sample_cninv, p0);
                        // a carrying contig makes the slot a carrier; non-carrier contigs count only where alone
                        if (!carrier && !lowest && (!agree || paths.size() > 1)) continue;
                        ++an;
                        if (carrier) ++ac;
                    }
                }
                if (an == 0 || static_cast<double>(ac) / static_cast<double>(an) < options.min_maf)
                    return false;
            }
            return true;
        };
        // A module record claims the repeats inside its copies. Checked when the record is made, so only a
        // record that is written claims them; a record dropped afterwards (tangle, oversize,
        // --min-haplotypes) would leave the copy change reported nowhere. A tangle or an oversized module
        // has no copy number at all; one under the carrier filters still gives its genes' copy numbers.
        // The span bound is on the record's reference span, the module's, not on one copy's size.
        auto module_counts = [&]() {
            if (is_tangle) { summary.tangle_ids.push_back(bubble.id); return false; }
            {
                bool amb = false;
                const auto mb = module_ref_bounds(&amb);
                // the bases between the module's two boundaries, as the record's POS..END
                const std::size_t ends = node_len(graph, bubble.source) + node_len(graph, bubble.sink);
                const std::size_t whole = mb.second >= mb.first && mb.first > 0 ? mb.second - mb.first + 1 : 0;
                const std::size_t span = whole > ends ? whole - ends : 0;
                if (span > max_dup_bp) { summary.oversized_ids.push_back(bubble.id); return false; }
            }
            return true;
        };
        // ---- A collapsed module's copies (--cn), counted as UNIT counts passes through its anchor. Not read
        // off the module's bp over a reference unit, which needs the reference to hold two copies and reads
        // a copy of another type as part of one, nor off the busiest single node, which a short repeat or a
        // novel insertion can be. Counted:
        //   - the module is the bubble's nodes some path visits more than once (the reference or any
        //     haplotype), leaving out the repeat units UNIT counts;
        //   - its anchor is the module nodes the most paths revisit, weighed by length, which are the ones
        //     present once in every copy whatever the copy's type: a node only some copies carry, one between
        //     copies, or content a few haplotypes duplicate apart from the module is revisited by fewer;
        //   - a path's copies are the length-weighted median of its visits over the anchor;
        //   - the reference's copy number is counted by the same rule, one copy or several.
        // Repeat units outside the module keep their UNIT records; those inside its copies (between the
        // reference's first and last module visit) stay in the walks, so a change within one copy is an
        // insertion or deletion and a whole extra copy goes with the copy change.
        std::unordered_set<std::string> unified_internal;   // repeat units inside the module's copies
        if (options.cn && ref_idx < graph.paths.size()) {
            const std::unordered_set<std::string> inside_set(bubble.inside.begin(), bubble.inside.end());
            std::vector<std::size_t> pis;   // the traversing paths
            std::size_t ref_k = SIZE_MAX;
            for (std::size_t pi = 0; pi < graph.paths.size(); ++pi) {
                if (!traverses.count(graph.paths[pi].name)) continue;
                if (pi == ref_idx) ref_k = pis.size();
                pis.push_back(pi);
            }
            std::vector<std::unordered_map<std::string, std::size_t>> visits(pis.size());
            for (std::size_t k = 0; k < pis.size(); ++k) {
                const auto sp = module_span(pis[k]);
                const std::vector<PathStep>& st = graph.paths[pis[k]].steps;
                for (std::size_t i = sp.first; i <= sp.second && i < st.size(); ++i) {
                    const std::string& id = st[i].node_id;
                    if (inside_set.count(id) && rep_unit_of(id) == nullptr) ++visits[k][id];
                }
            }
            std::vector<std::string> mod;
            {
                std::set<std::string> m;
                for (const auto& v : visits) for (const auto& [id, c] : v) if (c >= 2) m.insert(id);
                mod.assign(m.begin(), m.end());
            }
            auto visits_of = [&](std::size_t k, const std::string& id) -> std::size_t {
                const auto it = visits[k].find(id);
                return it == visits[k].end() ? 0 : it->second;
            };
            // length-weighted median of path k's visits over nodes (ties to the lower count)
            auto wmed = [&](std::size_t k, const std::vector<std::string>& nodes) -> std::size_t {
                std::vector<std::pair<std::size_t, std::size_t>> cw;
                std::size_t total = 0;
                for (const std::string& id : nodes) {
                    const std::size_t len = node_len(graph, id);
                    cw.emplace_back(visits_of(k, id), len);
                    total += len;
                }
                std::sort(cw.begin(), cw.end());
                std::size_t acc = 0;
                for (const auto& [c, w] : cw) { acc += w; if (2 * acc >= total) return c; }
                return cw.empty() ? 0 : cw.back().first;
            };
            if (!mod.empty() && ref_k != SIZE_MAX) {
                // The anchor: the module nodes the most paths revisit. A node in every copy is revisited by every
                // path holding two copies or more; one only some copies carry, one between copies, or content
                // a few haplotypes duplicate apart from the module, by fewer.
                // Weighed by length, so that a short element in more copies than the module does not stand
                // for it: the anchor is the nodes within a tenth of the revisit count that gathers the most
                // revisited bp.
                // A node every path visits equally often is a duplication fixed in the panel: it says nothing
                // about copy number, and counted at it a module beside it would read no change.
                std::vector<std::size_t> revisits(mod.size(), 0);
                for (std::size_t j = 0; j < mod.size(); ++j) {
                    bool fixed = true;
                    for (std::size_t k = 0; k < pis.size(); ++k) {
                        if (visits_of(k, mod[j]) >= 2) ++revisits[j];
                        if (visits_of(k, mod[j]) != visits_of(0, mod[j])) fixed = false;
                    }
                    if (fixed) revisits[j] = 0;
                }
                auto near = [](std::size_t r, std::size_t r0) { return 10 * (r > r0 ? r - r0 : r0 - r) <= r0; };
                std::size_t r_best = 0, mass_best = 0;
                for (const std::size_t r0 : std::set<std::size_t>(revisits.begin(), revisits.end())) {
                    std::size_t mass = 0;
                    for (std::size_t j = 0; j < mod.size(); ++j)
                        if (near(revisits[j], r0)) mass += revisits[j] * node_len(graph, mod[j]);
                    if (r0 > 0 && mass > mass_best) { mass_best = mass; r_best = r0; }
                }
                std::vector<std::string> anchor;
                std::size_t anchor_bp = 0, mod_bp = 0;
                for (std::size_t j = 0; j < mod.size(); ++j) {
                    const std::size_t len = node_len(graph, mod[j]);
                    mod_bp += len;
                    if (r_best > 0 && near(revisits[j], r_best)) { anchor.push_back(mod[j]); anchor_bp += len; }
                }
                const std::vector<std::string>& by = anchor;
                std::vector<std::size_t> copies(pis.size());
                for (std::size_t k = 0; k < pis.size(); ++k) copies[k] = wmed(k, by);
                const std::size_t ref_copies = copies[ref_k];
                // The repeat units inside the module's copies, which a path walks between its first and last
                // module visit: part of each copy, so part of its bp, though their copies are not counted apart
                // (they stay in the walks). Read from every path, not the reference alone, which may hold no copy.
                std::unordered_set<std::string> internal;
                const std::unordered_set<std::string> mod_set(mod.begin(), mod.end());
                for (std::size_t k = 0; k < pis.size(); ++k) {
                    const auto sp = module_span(pis[k]);
                    const std::vector<PathStep>& st = graph.paths[pis[k]].steps;
                    std::size_t lo = SIZE_MAX, hi = 0;
                    for (std::size_t i = sp.first; i <= sp.second && i < st.size(); ++i)
                        if (mod_set.count(st[i].node_id)) { lo = std::min(lo, i); hi = std::max(hi, i); }
                    for (std::size_t i = lo; lo != SIZE_MAX && i <= hi; ++i)
                        if (const std::string* u = rep_unit_of(st[i].node_id)) internal.insert(*u);
                }
                std::vector<std::size_t> internal_bp(pis.size(), 0);
                if (!internal.empty())
                    for (std::size_t k = 0; k < pis.size(); ++k) {
                        const auto sp = module_span(pis[k]);
                        const std::vector<PathStep>& st = graph.paths[pis[k]].steps;
                        for (std::size_t i = sp.first; i <= sp.second && i < st.size(); ++i) {
                            const std::string* u = inside_set.count(st[i].node_id) ? rep_unit_of(st[i].node_id) : nullptr;
                            if (u != nullptr && internal.count(*u)) internal_bp[k] += node_len(graph, st[i].node_id);
                        }
                    }
                auto module_bp = [&](std::size_t k) {
                    std::size_t bp = internal_bp[k];
                    for (const std::string& id : mod) bp += visits_of(k, id) * node_len(graph, id);
                    return bp;
                };
                // One copy's bp: the reference's module over its copies; with no reference copy, the median
                // over the haplotypes that carry one.
                double unit = 0.0;
                if (ref_copies > 0) {
                    unit = static_cast<double>(module_bp(ref_k)) / static_cast<double>(ref_copies);
                } else {
                    std::vector<double> per;
                    for (std::size_t k = 0; k < pis.size(); ++k)
                        if (copies[k] > 0) per.push_back(static_cast<double>(module_bp(k)) / static_cast<double>(copies[k]));
                    if (!per.empty()) {
                        std::sort(per.begin(), per.end());
                        unit = per[per.size() / 2];
                    }
                }
                MergedRecord mr;
                mr.seed.type = EvType::Dup;
                mr.seed.nodes.push_back(bubble.source);
                mr.seed.start_node = bubble.source; mr.seed.end_node = bubble.sink;
                mr.seed.ref_cn = ref_copies; mr.seed.alt_cn = ref_copies;
                mr.seed.anchor_node = bubble.source;
                mr.seed.size_bp = static_cast<std::size_t>(unit);
                mr.seed.ru_len = static_cast<std::size_t>(unit);
                mr.seed.cn_method = CnMethod::Module;
                mr.seed.shared_fold_bp = module_bp(ref_k);
                mr.seed.ref_fold = ref_copies;
                mr.seed.anchor_bp = anchor_bp;
                mr.seed.module_nodes_bp = mod_bp;
                mr.seed.module_span_ambiguous = module_ambiguous_paths;
                // A carrier's copy number differs from the reference's and its module changes by --min-sv-bp
                // or more. Copies x unit overstates a short module's change wherever the reference visits a
                // node of it more often than once per copy.
                const long long ref_mbp = static_cast<long long>(module_bp(ref_k));
                auto copies_bp = [&](std::size_t k) {
                    const auto sp = module_span(pis[k]);
                    const std::vector<PathStep>& st = graph.paths[pis[k]].steps;
                    std::size_t lo = SIZE_MAX, hi = 0;
                    for (std::size_t i = sp.first; i <= sp.second && i < st.size(); ++i)
                        if (mod_set.count(st[i].node_id)) { lo = std::min(lo, i); hi = std::max(hi, i); }
                    long long bp = 0;
                    for (std::size_t i = lo; lo != SIZE_MAX && i <= hi; ++i)
                        if (inside_set.count(st[i].node_id)) bp += static_cast<long long>(node_len(graph, st[i].node_id));
                    return bp;
                };
                const long long ref_cbp = copies_bp(ref_k);
                for (std::size_t k = 0; k < pis.size(); ++k) {
                    const std::string& name = graph.paths[pis[k]].name;
                    mr.sample_cn[name] = copies[k];
                    const long long dbp = static_cast<long long>(module_bp(k)) - ref_mbp;
                    mr.sample_bp[name] = dbp;
                    mr.sample_copies_bp[name] = copies_bp(k) - ref_cbp;
                    if (copies[k] != ref_copies && std::llabs(dbp) >= static_cast<long long>(options.min_sv_bp)) {
                        mr.carriers.push_back(name);
                        if (copies[k] > mr.seed.alt_cn) mr.seed.alt_cn = copies[k];
                    }
                }
                // One copy has to reach --min-sv-bp: a module of a few bp is pieces of a short
                // repeat, which UNIT counts where it is folded, and its copy changes are insertions where not.
                // And a module that is a repeat UNIT counts already -- its copies track the repeat's in nine
                // paths of ten and the repeat is half a copy or more -- is that repeat: UNIT counts it
                // exactly. A short element inside each copy is not.
                bool is_unit = false;
                {
                    std::unordered_map<std::string, std::size_t> allele_of;
                    for (std::size_t ai = 0; ai < alleles.size(); ++ai)
                        for (const std::string& m : alleles[ai].members) allele_of.emplace(m, ai);
                    for (const std::string& cn : cn_nodes) {
                        if (twin_of.count(cn) || !rep_callable(cn)) continue;
                        if (2.0 * static_cast<double>(unit_len(cn)) < unit) continue;
                        const std::vector<std::size_t>& ac = cn_allele_count.at(cn);
                        std::size_t same = 0;
                        for (std::size_t k = 0; k < pis.size(); ++k) {
                            const auto it = allele_of.find(graph.paths[pis[k]].name);
                            if (it != allele_of.end() && ac[it->second] == copies[k]) ++same;
                        }
                        if (same * 10 >= pis.size() * 9) { is_unit = true; break; }
                    }
                }
                // And so does what it counts, the anchor: an anchor under --min-sv-bp is no SV-scale copy,
                // and its carriers can be haplotypes deleting the whole site.
                if (!is_unit && unit >= static_cast<double>(options.min_sv_bp) &&
                    anchor_bp >= options.min_sv_bp && !mr.carriers.empty() && module_counts()) {
                    mr.anchor_nodes.insert(anchor.begin(), anchor.end());
                    mr.module_nodes.insert(mod.begin(), mod.end());
                    for (const std::string& id : bubble.inside)
                        if (const std::string* u = rep_unit_of(id); u != nullptr && internal.count(*u)) mr.module_nodes.insert(id);
                    if (support_ok(mr)) {
                        unified_internal = std::move(internal);
                        mr.member_nodes.insert(bubble.inside.begin(), bubble.inside.end());
                        module_record = true;
                        merged.push_back(std::move(mr));
                    } else {
                        dropped_modules.push_back(std::move(mr));
                    }
                }
            }
        }
        // ---- DUP/CN events: count self-loop traversals per allele vs reference. Merge on shared unit
        // node; per-sample CN. Independent of the walk-diff alignment; a loop inside a counted module's copies
        // is the module's (unified_internal).
        //
        // Node-outer, allele-inner, so that every allele's traversal count of a unit node becomes a copy
        // number, including alleles matching the reference: a reference-like traverser gets REF_CN by
        // measurement rather than by default, and the per-gene sidecar has its count to split.
        for (const std::string& cn : cn_nodes) {
            if (twin_of.count(cn) || unified_internal.count(cn)) continue;
            // No unit-length gate: the per-allele size test below is the one that matters (see
            // rep_callable). A short loop with no SV-sized copy change still produces no record.
            const std::size_t rc = rep_ref_count.count(cn) ? rep_ref_count.at(cn) : 0;
            const std::size_t unit = unit_len(cn);
            const std::vector<std::size_t>& ac_of = cn_allele_count.at(cn);
            // The record exists if any allele differs from the reference by a callable amount. Its
            // ALT_CN / SVLEN describe the largest copy change (ties: the larger copy number), so they
            // do not depend on the order haplotypes were enumerated in; per-sample CN is exact anyway.
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
                    mr.seed.module_span_ambiguous = module_ambiguous_paths;
                    mr.seed.nodes.push_back(cn);
                    mr.seed.start_node = cn; mr.seed.end_node = cn;
                    mr.seed.ref_cn = rc; mr.seed.alt_cn = ac;
                    mr.seed.anchor_node = bubble.source;
                    if (rc == 0) {
                        // No reference copy to place the record on: anchor on the reference step the
                        // copies follow in a carrying walk, which is where they are inserted. That step
                        // comes from aligning the walk to the reference, as an INS anchor does: the last
                        // walk step before the copies that aligned to a reference step. Neither a node
                        // name nor a visit count identifies it -- the walk may skip an earlier visit.
                        const std::vector<PathStep>& w = alleles[ai].steps;
                        std::size_t first_copy = w.size();
                        for (std::size_t i = 0; i < w.size(); ++i) {
                            const std::string* u = rep_unit_of(w[i].node_id);
                            if (u != nullptr && *u == cn) { first_copy = i; break; }
                        }
                        std::unordered_set<std::string> w_drop = cn_nodes;
                        if (const auto tw = rep_twins.find(cn); tw != rep_twins.end())
                            w_drop.insert(tw->second.begin(), tw->second.end());
                        for (const auto& [n, key] : ug_units(cn_nodes)) w_drop.insert(n);
                        const std::vector<Tok> Hw = collapse_walk(w, w_drop);
                        std::vector<std::size_t> hm;
                        diff_walks(graph, Rtok, Hw, bubble.source, &hm);
                        for (std::size_t t = Hw.size(); t-- > 0;) {
                            if (Hw[t].src_idx >= first_copy || hm[t] == SIZE_MAX) continue;
                            // Copy-number steps are not aligned, so walk on through any the walk and
                            // the reference share (an extra copy of another looped node, on either of its
                            // twins) up to the copies.
                            std::size_t rs = hm[t];
                            for (std::size_t i = Hw[t].src_idx + 1;
                                 i < first_copy && rs + 1 < ref_steps.size() &&
                                 same_unit(ref_steps[rs + 1].node_id, w[i].node_id) && ref_steps[rs + 1].reverse == w[i].reverse;
                                 ++i) ++rs;
                            mr.seed.ref_tok_anchor = rs;
                            mr.seed.anchor_node = ref_steps[rs].node_id;
                            break;
                        }
                    }
                    mr.seed.size_bp = unit * delta;
                    mr.seed.ru_len = unit;
                    mr.seed.cn_method = CnMethod::Unit;
                    // variant_nodes.tsv lists every node a copy can sit on, so node-level dosage reads
                    // a path routed through a twin as the copies it is.
                    if (const auto tw = rep_twins.find(cn); tw != rep_twins.end()) {
                        mr.member_nodes.insert(cn);
                        mr.member_nodes.insert(tw->second.begin(), tw->second.end());
                    }
                    if (const auto um = ug_members.find(cn); um != ug_members.end())
                        mr.member_nodes.insert(um->second.begin(), um->second.end());
                    merged.push_back(std::move(mr));
                    grp = &merged.back();
                }
                const std::size_t cur = grp->seed.alt_cn;
                const std::size_t cur_delta = cur > rc ? cur - rc : rc - cur;
                if (delta > cur_delta || (delta == cur_delta && ac > cur)) {
                    grp->seed.alt_cn = ac;
                    grp->seed.size_bp = unit * delta;
                }
            }
            if (grp == nullptr) continue;
            const std::vector<std::size_t>& inv_of = cn_allele_inv.at(cn);
            std::size_t ref_inv = 0;
            for (const PathStep& s : ref_steps) {
                const std::string* u = rep_unit_of(s.node_id);
                if (u != nullptr && *u == cn && copy_rev(s) != cn_ref_reverse.at(cn)) ++ref_inv;
            }
            // Now fill every allele, reference-like ones included. Carriers stay the differing ones, so
            // AC/AF are unchanged; what changes is that a traverser always has a measured CN.
            for (std::size_t ai = 0; ai < alleles.size(); ++ai) {
                const std::size_t ac = ac_of[ai];
                const std::size_t delta = rc > ac ? rc - ac : ac - rc;
                // Carrier keeps exactly the rule that produced this record -- a callable difference --
                // so AC/AF are untouched. Completeness is about sample_cn, not about who counts as a
                // carrier: a sub-threshold difference is a real CN and not a called SV.
                const bool carrier = (ac != rc) && (unit * delta >= options.min_sv_bp);
                // An inverted duplication: a gain whose haplotype runs more copies against the
                // reference's orientation than the reference itself does. With no reference copy, a
                // haplotype's copies are all new: one is an inverted duplication only where it holds a copy
                // and its reverse complement, not where its copies merely run against the other haplotypes'.
                const bool inv_gain = rc > 0 ? inv_of[ai] > ref_inv
                                             : cn_allele_minor.count(cn) && cn_allele_minor.at(cn)[ai] > 0;
                if (carrier && ac > rc && inv_gain) grp->seed.inverted_dup = true;
                for (const std::string& m : alleles[ai].members) {
                    grp->sample_cn[m] = ac;
                    grp->sample_cninv[m] = inv_of[ai];
                    if (carrier) grp->carriers.push_back(m);
                }
            }
            // The record counts this loop's copies only if it survives, so apply the support filter the
            // final pass will (support_ok, over the record's carriers).
            if (support_ok(*grp)) {
                rep_consumed.insert(cn);
                if (const auto tw = rep_twins.find(cn); tw != rep_twins.end())
                    for (const std::string& t : tw->second) if (cn_nodes.count(t)) rep_consumed.insert(t);
            }
        }

        // DEL/INS/INV per allele, kept down to the rescue floor, then merged. Emitted beside a module record:
        // its DUP is dosage, not the sequence-resolved indels inside the bubble. The
        // per-allele walk-diff dominates on big folded bubbles and writes only its own slot, so it runs
        // across cores; the merge below consumes in allele order, independent of thread count.
        std::vector<std::vector<Event>> allele_events(alleles.size());
        // Per allele, the node comparison's matched pairs (allele step, reference step), in order: where
        // a node event's two walks meet again on either side. carrier_events.tsv writes each node event as
        // the stretch between those anchors on both walks, so a copy-number unit the comparison left out
        // goes with each side's own copies rather than being spanned on one side only.
        std::vector<std::vector<std::pair<std::size_t, std::size_t>>> allele_matched(alleles.size());
        // A loop node leaves the walks only where a surviving UNIT record counts its copies; dropped
        // elsewhere, a sub-threshold loop's copies would be missing from the events around them. Kept, it
        // is an ordinary node, and a sub-threshold change of it a sub-threshold INS/DEL that the size
        // filter removes.
        // A counted module leaves the repeat units inside its copies in the walks: only those with a UNIT
        // record of their own go.
        const std::unordered_set<std::string>& walk_drop = rep_consumed;
        // Twins compare as their loop node; a unit graph whose copies are counted leaves with all its nodes.
        std::unordered_map<std::string, std::string> walk_unit_of = twin_of;
        for (const auto& [n, key] : ug_units(walk_drop)) walk_unit_of.emplace(n, key);
        const std::vector<Tok> Rdiff = collapse_walk(ref_steps, walk_drop, walk_unit_of);
        std::unordered_map<std::size_t, std::size_t> rdiff_pos;   // reference step index -> token
        for (std::size_t k = 0; k < Rdiff.size(); ++k) rdiff_pos.emplace(Rdiff[k].src_idx, k);
        // An allele whose walk held a stretch too large to align is only partly compared: a variant
        // there produced no event, so its haplotypes cannot be called reference on this bubble's
        // INS/DEL/INV records (GT . where they carry nothing). One slot per allele, so the parallel
        // writers never share one.
        std::vector<char> allele_unaligned(alleles.size(), 0);
        std::vector<char> allele_seq(alleles.size(), 0);   // its events come from aligning its bases
        // Such an allele's alignment: every change it found, small ones too, in offsets of the two compared
        // sequences (the reference's and the haplotype's interiors with the counted copies left out, genomic
        // orientation), so a copy stretch can be cut where the alignment cuts the bases beside the copies.
        std::vector<std::vector<AlignedChange>> allele_align(alleles.size());
        // Blocks too large to align that no base-level comparison of the allele then covered.
        std::vector<std::size_t> allele_skipped_blocks(alleles.size(), 0);
        // Alleles compared by sequence. Each haplotype's interior is aligned to the reference's at base
        // level, in genomic orientation, and its indels are the events: the graph's node boundaries place
        // an indel arbitrarily inside a repeat, and can split one over several nodes, which the node
        // comparison then reports as pieces. Copies a record counts have left both walks (walk_drop), so
        // what is aligned is the rest, and the reference's rest is the pieces of the reference between
        // the removed copies, each one stretch of the genome. A collapsed module (MODULE) is compared
        // node by node: its copies stay in the
        // walks, and aligned, an extra copy and a novel insertion beside it are one gap, which the copy
        // change then removes whole; node by node they are different nodes.
        // offset in seq_ref, genomic first base there, and how far left an event in it may move: to the
        // start of the piece where removed copies lie before it, else as far as the reference allows
        struct SeqPiece { std::size_t off, g0, floor; };
        bool seq_ok = options.align_alleles && !module_record && ref_steps.size() >= 2 && !ref_genome.empty();
        std::string seq_ref;             // the reference interior left in the comparison, genomic orientation
        std::vector<SeqPiece> seq_pieces;
        std::vector<std::pair<std::size_t, std::string>> seq_ref_nodes;   // offset in seq_ref -> interior node
        if (seq_ok) {
            const std::size_t n = ref_steps.size();
            const std::size_t base = ref_meta.region_start_1based;
            std::vector<std::size_t> kept;   // interior steps left in the comparison, in genomic order
            for (std::size_t k = 1; k + 1 < n; ++k) if (rdiff_pos.count(k)) kept.push_back(k);
            if (ref_rev) std::reverse(kept.begin(), kept.end());
            std::size_t prev_last = 0;
            for (const std::size_t k : kept) {
                const auto sp = tok_span(k, k);
                const std::size_t len = node_len(graph, ref_steps[k].node_id);
                if (sp.first <= base || sp.second + 1 - sp.first != len) { seq_ok = false; break; }
                std::string nw = spell_path_steps_sequence(graph, std::vector<PathStep>{ref_steps[k]});
                for (char& c : nw) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
                if (ref_genome.compare(sp.first - base, len, ref_rev ? reverse_complement(nw) : nw) != 0) {
                    seq_ok = false;
                    break;
                }
                if (seq_pieces.empty() || sp.first != prev_last + 1) {
                    const bool first_base = seq_pieces.empty() && sp.first == tok_span(1, n - 2).first;
                    seq_pieces.push_back({seq_ref.size(), sp.first, first_base ? 0 : sp.first});
                }
                seq_ref_nodes.emplace_back(seq_ref.size(), ref_steps[k].node_id);
                seq_ref.append(ref_genome, sp.first - base, len);
                prev_last = sp.second;
            }
            // walk order, as the node comparison lists an event's nodes
            if (ref_rev) std::reverse(seq_ref_nodes.begin(), seq_ref_nodes.end());
            if (seq_ok && seq_pieces.empty()) {
                // Nothing of the interior left: the junction between the two boundaries, or the place the
                // removed copies start, which is where an insertion beside them goes.
                std::size_t g = 0;
                if (n > 2) {
                    g = tok_span(1, n - 2).first;
                } else {
                    const auto a = tok_span(0, 0), b = tok_span(n - 1, n - 1);
                    const auto& first = ref_rev ? b : a;
                    if (first.second != 0 && first.second >= base) g = first.second + 1;
                }
                if (g <= base) seq_ok = false;
                else seq_pieces.push_back({0, g, 0});
            }
        }
        // The events of allele ai by aligning it, or false to keep the node comparison's. Only the part
        // between the alleles' shared start and end is aligned. That core goes to a banded global
        // alignment when it fits: the band reaches as far from the diagonal as the node comparison's own
        // events take the walk, plus 64, so a placement no worse than the graph's is always inside it.
        // A core too large for that (a long indel widens the band) goes to minimap2, anchored by up to
        // 200 bp of the shared sequence either side, and is used only if aligned end to end.
        auto seq_events = [&](std::size_t ai, const std::vector<Tok>& htok, const std::vector<Event>& node_events,
                              bool node_skipped, std::vector<Event>& out) {
            out.clear();
            const std::vector<PathStep>& w = alleles[ai].steps;
            if (w.size() < 2 || htok.size() < 2 || htok.front().src_idx != 0 || htok.back().src_idx + 1 != w.size())
                return false;
            // The haplotype's interior left in the comparison, and where copies left it (a junction).
            std::vector<PathStep> hin;
            std::vector<std::size_t> hap_junctions;   // walk-orientation offsets in hw
            {
                std::size_t off = 0;
                for (std::size_t j = 1; j + 1 < htok.size(); ++j) {
                    if (j > 1 && htok[j].src_idx != htok[j - 1].src_idx + 1) hap_junctions.push_back(off);
                    hin.push_back(PathStep{htok[j].node_id, htok[j].reverse});
                    off += node_len(graph, htok[j].node_id);
                }
            }
            std::string hw = hin.empty() ? std::string() : spell_path_steps_sequence(graph, hin);
            for (char& c : hw) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
            const std::string hs = ref_rev ? reverse_complement(hw) : hw;
            if (ref_rev) for (std::size_t& j : hap_junctions) j = hw.size() - j;
            const std::string& rs = seq_ref;
            std::size_t p = 0, sfx = 0;
            while (p < rs.size() && p < hs.size() && rs[p] == hs[p]) ++p;
            while (sfx < rs.size() - p && sfx < hs.size() - p && rs[rs.size() - 1 - sfx] == hs[hs.size() - 1 - sfx]) ++sfx;
            const std::size_t n = rs.size() - p - sfx, m = hs.size() - p - sfx;
            if (n == 0 && m == 0) { allele_align[ai].clear(); return true; }
            // A substitution says nothing about where the walk drifts inside the run it replaces: one
            // node can hold an insertion and a deletion far apart. So it may reach as far as its larger side.
            long long drift = 0, reach = 0;
            for (const Event& e : node_events) {
                if (e.type == EvType::Sub) {
                    const long long side = static_cast<long long>(std::max(e.del_bp, e.ins_bp));
                    reach = std::max(reach, (drift < 0 ? -drift : drift) + side);
                }
                drift += e.type == EvType::Ins ? static_cast<long long>(e.size_bp)
                       : e.type == EvType::Del ? -static_cast<long long>(e.size_bp)
                       : e.type == EvType::Sub ? static_cast<long long>(e.ins_bp) - static_cast<long long>(e.del_bp) : 0;
                reach = std::max(reach, drift < 0 ? -drift : drift);
            }
            reach = std::min<long long>(reach, static_cast<long long>(std::max(n, m)));
            const std::size_t band = std::max<std::size_t>(n > m ? n - m : m - n, static_cast<std::size_t>(reach)) + 64;
            std::vector<AlignedChange> changes;
            std::size_t at = p;   // where the aligned stretch starts in both alleles
            // A core with one side empty is one insertion or deletion: nothing to align. Where the node
            // comparison skipped the stretch, its events say nothing about how far the walk drifts, so the
            // band is unknown: minimap2 aligns it instead.
            if (n == 0 || m == 0) {
                changes.push_back({0, n, 0, m});
            } else if (!node_skipped && std::max(n, m) * (2 * band + 1) <= kSeqAlignCellCap) {
                if (!aligned_changes(rs.substr(p, n), hs.substr(p, m), static_cast<int>(band), changes)) return false;
            } else {
                const std::size_t fl = std::min<std::size_t>(p, 200), fr = std::min<std::size_t>(sfx, 200);
                const std::string rw = rs.substr(p - fl, n + fl + fr), qw = hs.substr(p - fl, m + fl + fr);
                const Minimap2Hit hit = minimap2_best_hit("hap", qw, "ref", rw, "asm20", 1);
                if (!hit.ok || hit.reverse || hit.cigar.empty()) return false;
                // minimap2 aligns locally: sequence either allele leaves unaligned at an end is a gap
                // there (a core opening on a long deletion aligns the rest and leaves the deletion out).
                std::vector<std::pair<int, int>> cig;
                if (hit.query_start_bp > 0) cig.emplace_back(1, static_cast<int>(hit.query_start_bp));
                if (hit.target_start_bp > 0) cig.emplace_back(2, static_cast<int>(hit.target_start_bp));
                cig.insert(cig.end(), hit.cigar.begin(), hit.cigar.end());
                if (hit.query_end_bp < qw.size()) cig.emplace_back(1, static_cast<int>(qw.size() - hit.query_end_bp));
                if (hit.target_end_bp < rw.size()) cig.emplace_back(2, static_cast<int>(rw.size() - hit.target_end_bp));
                if (!changes_from_cigar(cig, rw, qw, changes)) return false;
                at = p - fl;   // window coordinates (a gap may start in a flank)
            }
            // Haplotype interior nodes by genomic-orientation offset, for each event's node set.
            std::vector<std::pair<std::size_t, std::string>> hap_nodes;
            {
                std::size_t off = 0;
                for (const PathStep& st : hin) {
                    const std::size_t len = node_len(graph, st.node_id);
                    hap_nodes.emplace_back(ref_rev ? hw.size() - off - len : off, st.node_id);
                    off += len;
                }
            }
            auto nodes_over = [&](const std::vector<std::pair<std::size_t, std::string>>& ns, std::size_t lo,
                                  std::size_t hi, std::vector<std::string>& into) {
                for (const auto& [o, id] : ns) {
                    const std::size_t len = node_len(graph, id);
                    if (o < hi && o + len > lo) into.push_back(id);
                }
            };
            for (AlignedChange c : changes) {
                // The bases the two sides share at either end are not part of the change (minimap2's
                // terminal gaps take in the shared flank), so an insertion or deletion is typed as one.
                while (c.t0 < c.t1 && c.q0 < c.q1 && rs[at + c.t0] == hs[at + c.q0]) { ++c.t0; ++c.q0; }
                while (c.t1 > c.t0 && c.q1 > c.q0 && rs[at + c.t1 - 1] == hs[at + c.q1 - 1]) { --c.t1; --c.q1; }
                const std::size_t dl = c.t1 - c.t0, il = c.q1 - c.q0;
                const std::size_t size = std::max(dl, il);
                if (size == 0) continue;
                const std::size_t g0 = at + c.t0;                // ref interior offset of the change
                // A change reaching across copies that left either walk is not one place: the allele
                // keeps the node comparison, which splits it there. Checked before the size floor, as
                // a small change is as much a sign that the pieces are not lined up.
                for (std::size_t k = 1; k < seq_pieces.size(); ++k)
                    if (g0 < seq_pieces[k].off && seq_pieces[k].off < g0 + dl) return false;
                for (const std::size_t j : hap_junctions)
                    if (at + c.q0 < j && j < at + c.q1) return false;
                if (size < rescue_floor) continue;
                Event e;
                e.seq_event = true;
                // The piece of the reference it lies in; an insertion where copies left the reference
                // goes on the piece before them, at the copies' start.
                std::size_t pc = 0;
                while (pc + 1 < seq_pieces.size() &&
                       (seq_pieces[pc + 1].off < g0 || (seq_pieces[pc + 1].off == g0 && dl > 0))) ++pc;
                const std::size_t gx = seq_pieces[pc].g0 + (g0 - seq_pieces[pc].off);   // genomic, at g0
                if (gx <= ref_meta.region_start_1based) return false;   // no base before it
                std::string d = rs.substr(g0, dl), a = hs.substr(at + c.q0, il);
                e.seq_before = gx - 1;
                // Where it sits once moved as far left as the reference allows, for merging: the same
                // indel in two haplotypes, placed differently by the graph, then has one position. Not
                // past the start of its piece, into copies the walks no longer hold.
                std::size_t x = gx;                               // first deleted base / insertion point
                if (dl == 0 || il == 0) {
                    std::string u = dl == 0 ? a : d;
                    const std::size_t base = ref_meta.region_start_1based;
                    while (x > base + 1 && x > seq_pieces[pc].floor && ref_genome[x - 1 - base] == u.back() &&
                           (dl == 0 || ref_genome[x - 1 - base] == ref_genome[x - 1 - base + dl])) {
                        u = u.back() + u.substr(0, u.size() - 1);
                        --x;
                    }
                }
                e.ref_pos = x;
                if (dl == 0) { e.type = EvType::Ins; e.seq = ref_rev ? reverse_complement(a) : a; e.size_bp = il; }
                else if (il == 0) { e.type = EvType::Del; e.seq = ref_rev ? reverse_complement(d) : d; e.size_bp = dl; }
                else {
                    e.type = EvType::Sub;
                    e.seq = ref_rev ? reverse_complement(d) : d;
                    e.alt_seq = ref_rev ? reverse_complement(a) : a;
                    e.size_bp = size;
                }
                e.del_bp = dl; e.ins_bp = il;
                e.seq_ins = a;
                nodes_over(seq_ref_nodes, g0, g0 + dl, e.nodes);
                nodes_over(hap_nodes, at + c.q0, at + c.q1, e.nodes);
                if (e.nodes.empty()) e.nodes.push_back(bubble.source);
                e.start_node = e.nodes.front();
                e.end_node = e.nodes.back();
                e.anchor_node = bubble.source;
                out.push_back(std::move(e));
            }
            allele_align[ai].clear();
            for (const AlignedChange& c : changes) allele_align[ai].push_back({at + c.t0, at + c.t1, at + c.q0, at + c.q1});
            return true;
        };
        // A node replacement whose two sides share bases at either end (different nodes spelling the same
        // bases) and that is then a plain insertion or deletion is made one, placed as an aligned change is
        // (left-aligned along the reference), so it merges with the same change typed as one and is
        // classified as one. Its size and nodes are unchanged. Not at a collapsed module, whose copy-change
        // accounting reads the node comparison's events as they are.
        auto plain_indel = [&](Event& e) {
            const std::size_t base = ref_meta.region_start_1based;
            if (e.splice_before < base || e.splice_before + 1 - base + e.splice_del > ref_genome.size()) return;
            std::string d = ref_genome.substr(e.splice_before + 1 - base, e.splice_del);
            std::string a = e.splice_ins;
            std::size_t pre = 0;
            while (pre < d.size() && pre < a.size() && d[pre] == a[pre]) ++pre;
            std::size_t suf = 0;
            while (suf < d.size() - pre && suf < a.size() - pre && d[d.size() - 1 - suf] == a[a.size() - 1 - suf]) ++suf;
            const std::size_t dl = d.size() - pre - suf, il = a.size() - pre - suf;
            if ((dl != 0) == (il != 0)) return;   // still a replacement, or no change at all
            d = d.substr(pre, dl);
            a = a.substr(pre, il);
            const std::size_t gx = e.splice_before + 1 + pre;
            // Not past the start of the stretch of reference the walks still hold, into copies they no longer
            // do, as an aligned change is bounded.
            std::size_t floor = 0;
            if (e.r_lo != SIZE_MAX && e.r_lo < Rdiff.size()) {
                std::size_t k = e.r_lo, j = e.r_lo;
                while (k > 0 && Rdiff[k - 1].src_idx + 1 == Rdiff[k].src_idx) --k;
                while (j + 1 < Rdiff.size() && Rdiff[j + 1].src_idx == Rdiff[j].src_idx + 1) ++j;
                // only where the walks' reference leaves off on the genomic left: before the stretch in walk order
                // on a forward reference, after it on one walked backwards
                const bool gap_left = ref_rev ? j + 1 < Rdiff.size() : k > 0;
                if (gap_left) floor = tok_span(Rdiff[k].src_idx, Rdiff[j].src_idx).first;
            }
            std::size_t x = gx;
            std::string u = dl == 0 ? a : d;
            while (x > base + 1 && x > floor && ref_genome[x - 1 - base] == u.back() &&
                   (dl == 0 || ref_genome[x - 1 - base] == ref_genome[x - 1 - base + dl])) {
                u = u.back() + u.substr(0, u.size() - 1);
                --x;
            }
            e.seq_event = true;
            e.node_found = true;
            e.ref_tok_first = e.ref_tok_last = SIZE_MAX;   // placed by its bases, as an aligned change is
            e.seq_before = gx - 1;
            e.ref_pos = x;
            e.type = dl == 0 ? EvType::Ins : EvType::Del;
            e.seq = dl == 0 ? (ref_rev ? reverse_complement(a) : a) : (ref_rev ? reverse_complement(d) : d);
            e.alt_seq.clear();
            e.del_bp = dl; e.ins_bp = il;
            e.size_bp = std::max(dl, il);
            e.seq_ins = a;
        };
        // The module record, when this site has one: a replacement of a whole copy is split below, and each
        // carrier's copies are told apart from its other events (split_copies).
        const MergedRecord* module_rec = nullptr;
        std::size_t module_idx = SIZE_MAX;
        for (std::size_t ri = 0; ri < merged.size(); ++ri)
            if (!merged[ri].anchor_nodes.empty()) { module_rec = &merged[ri]; module_idx = ri; }
        const std::unordered_set<std::string> module_gated =
            module_rec == nullptr ? std::unordered_set<std::string>()
                                  : std::unordered_set<std::string>(module_rec->carriers.begin(), module_rec->carriers.end());
        std::vector<CopySplit> copy_split(alleles.size());
        std::vector<long long> cnbp_fix(alleles.size(), 0);   // CNBP's change where a walk replaced the copies' end nodes
        const std::unordered_set<std::string> interior(bubble.inside.begin(), bubble.inside.end());
        std::unordered_map<std::string, std::size_t> ref_tok_of_node;   // a node the reference walks: one place it does
        for (std::size_t k = 0; k < Rdiff.size(); ++k) ref_tok_of_node.emplace(Rdiff[k].node_id, k);
        run_parallel(alleles.size(), options.threads, [&](std::size_t ai) {
            if (build_walk_signature(alleles[ai].steps) == ref_sig) return;
            std::vector<Event> events;
            {
                const std::vector<Tok> Htok = collapse_walk(alleles[ai].steps, walk_drop, walk_unit_of);
                std::size_t skipped = 0;
                std::vector<std::size_t> hm;
                // Only whole copies leave the walks, for a UNIT record (every node of the unit graph or loop),
                // so a gap is split where they left it, beside a module record too: the module takes nothing
                // out of the walks, but a UNIT array beside it still leaves junction pieces.
                events = diff_walks(graph, Rdiff, Htok, bubble.source, options.write_variant_nodes || module_rec != nullptr ? &hm : nullptr, &skipped,
                                    true, options.min_sv_bp, rescue_floor);
                for (std::size_t j = 0; j < hm.size(); ++j)
                    if (hm[j] != SIZE_MAX) allele_matched[ai].emplace_back(Htok[j].src_idx, hm[j]);
                // CNBP runs from each walk's first module node to its last. A walk whose copies start (end) at another
                // module node than the reference's, after (before) the step both walks share right before (after) the
                // reference's copies, holds other sequence in place of the reference's leading (trailing) module nodes:
                // that sequence is in its copies. An insertion before copies that start as the reference's is beside
                // them, and a deletion of module nodes with nothing in their place changes nothing here.
                std::size_t h_first = SIZE_MAX, h_last = SIZE_MAX;   // the copy span widened so (split_copies)
                if (module_rec != nullptr && !hm.empty()) {
                    std::size_t ra = 0, rb = 0, ha = 0, hb = 0;
                    if (copy_span(Rdiff, module_rec->module_nodes, ra, rb) && copy_span(Htok, module_rec->module_nodes, ha, hb)) {
                        std::vector<std::size_t> hr(Htok.size(), SIZE_MAX);   // the reference token each is matched to
                        for (std::size_t j = 0; j < hm.size() && j < Htok.size(); ++j) {
                            if (hm[j] == SIZE_MAX) continue;
                            const auto it = rdiff_pos.find(hm[j]);
                            if (it != rdiff_pos.end()) hr[j] = it->second;
                        }
                        const auto bp_of = [&](const std::vector<Tok>& W, std::size_t lo, std::size_t hi) {   // [lo, hi)
                            long long b = 0;
                            for (std::size_t k = lo; k < hi; ++k)
                                if (interior.count(W[k].node_id)) b += static_cast<long long>(node_len(graph, W[k].node_id));
                            return b;
                        };
                        long long fix = 0;
                        if (Htok[ha].token != Rdiff[ra].token) {
                            std::size_t j0 = SIZE_MAX;
                            for (std::size_t j = 0; j < ha; ++j) if (hr[j] != SIZE_MAX && hr[j] < ra) j0 = j;
                            if (j0 != SIZE_MAX && hr[j0] + 1 == ra && j0 + 1 < ha) { fix += bp_of(Htok, j0 + 1, ha); h_first = j0 + 1; }
                        }
                        if (Htok[hb].token != Rdiff[rb].token) {
                            std::size_t j1 = SIZE_MAX;
                            for (std::size_t j = hb + 1; j < Htok.size() && j1 == SIZE_MAX; ++j) if (hr[j] != SIZE_MAX && hr[j] > rb) j1 = j;
                            if (j1 != SIZE_MAX && hr[j1] == rb + 1 && hb + 1 < j1) { fix += bp_of(Htok, hb + 1, j1); h_last = j1 - 1; }
                        }
                        cnbp_fix[ai] = fix;
                    }
                }
                if (skipped > 0) allele_unaligned[ai] = 1;
                allele_skipped_blocks[ai] = skipped;
                if (module_rec != nullptr && !alleles[ai].members.empty()) {
                    const auto cit = module_rec->sample_cn.find(alleles[ai].members.front());
                    if (cit != module_rec->sample_cn.end())
                        split_copy_replacements(graph, Rdiff, Htok, module_rec->anchor_nodes,
                                                static_cast<long long>(cit->second) - static_cast<long long>(module_rec->seed.ref_cn),
                                                module_rec->seed.ru_len, options.min_sv_bp, events);
                    // A carrier by the gate (CN differs, its module changes by --min-sv-bp): its copies, as edits
                    // of the reference, and its other events.
                    const std::string& m0 = alleles[ai].members.front();
                    const auto bit = module_rec->sample_bp.find(m0);
                    if (cit != module_rec->sample_cn.end() && module_gated.count(m0) && bit != module_rec->sample_bp.end()) {
                        const auto place = [&](const Event& pc, CopyForm& f) -> bool {
                            if (pc.r_lo == SIZE_MAX || pc.h_lo == SIZE_MAX) return false;
                            std::vector<const Tok*> ht;
                            for (std::size_t k = pc.h_lo; k < pc.h_hi; ++k) ht.push_back(&Htok[k]);
                            std::string seq = ht.empty() ? std::string() : spell_toks(graph, ht);
                            if (ref_rev) seq = reverse_complement(seq);
                            for (char& ch : seq) ch = upper_base(ch)[0];
                            if (pc.r_lo < pc.r_hi) {
                                const auto sp = tok_span(Rdiff[pc.r_lo].src_idx, Rdiff[pc.r_hi - 1].src_idx);
                                std::size_t bp = 0;
                                for (std::size_t k = pc.r_lo; k < pc.r_hi; ++k) bp += node_len(graph, Rdiff[k].node_id);
                                if (sp.first == 0 || sp.second + 1 - sp.first != bp) return false;   // not one stretch
                                f.start = sp.first; f.end = sp.second;
                            } else {
                                if (pc.r_lo == 0) return false;
                                const auto sp = tok_span(Rdiff[pc.r_lo - 1].src_idx, Rdiff[pc.r_lo - 1].src_idx);
                                if (sp.first == 0) return false;
                                const std::size_t before = ref_rev ? sp.first - 1 : sp.second;   // the base it follows
                                f.start = before + 1; f.end = before;
                            }
                            // The inserted bases as the reference stretches they repeat: each token of a node the
                            // reference walks, at the place it walks it, in walk order; kept where it spells the same
                            // bases and is shorter.
                            struct Seg { std::size_t a = 0, b = 0; bool rc = false; std::string lit; };
                            std::vector<Seg> segs;
                            for (const Tok* t : ht) {
                                const auto it = ref_tok_of_node.find(t->node_id);
                                if (it == ref_tok_of_node.end()) { segs.push_back({0, 0, false, spell_toks(graph, {t})}); continue; }
                                const Tok& r = Rdiff[it->second];
                                const auto sp = tok_span(r.src_idx, r.src_idx);
                                if (sp.first == 0) { segs.push_back({0, 0, false, spell_toks(graph, {t})}); continue; }
                                segs.push_back({sp.first, sp.second, (t->reverse != r.reverse) != ref_rev, std::string()});
                            }
                            if (ref_rev) {   // the copy in reference orientation: the walk's pieces backwards, each turned
                                std::reverse(segs.begin(), segs.end());
                                for (Seg& g : segs) { if (g.a == 0) g.lit = reverse_complement(g.lit); else g.rc = !g.rc; }
                            }
                            std::vector<Seg> merged_segs;
                            for (Seg& g : segs) {
                                if (!merged_segs.empty() && g.a != 0 && merged_segs.back().a != 0 && merged_segs.back().rc == g.rc &&
                                    (g.rc ? g.b + 1 == merged_segs.back().a : merged_segs.back().b + 1 == g.a)) {
                                    if (g.rc) merged_segs.back().a = g.a; else merged_segs.back().b = g.b;
                                } else if (!merged_segs.empty() && g.a == 0 && merged_segs.back().a == 0) merged_segs.back().lit += g.lit;
                                else merged_segs.push_back(std::move(g));
                            }
                            std::string code, made;
                            std::size_t lit_bp = 0;
                            const std::size_t base = ref_meta.region_start_1based;
                            for (const Seg& g : merged_segs) {
                                if (!code.empty()) code += '+';
                                if (g.a == 0) {
                                    std::string l = g.lit;
                                    for (char& ch : l) ch = upper_base(ch)[0];
                                    code += l; made += l; lit_bp += l.size();
                                    continue;
                                }
                                code += (g.rc ? '~' : '=') + std::to_string(g.a) + '-' + std::to_string(g.b);
                                if (g.a < base || g.b + 1 - base > ref_genome.size()) { made.clear(); break; }
                                std::string x = ref_genome.substr(g.a - base, g.b + 1 - g.a);
                                for (char& ch : x) ch = upper_base(ch)[0];
                                made += g.rc ? reverse_complement(x) : x;
                            }
                            f.seq = std::move(seq);
                            if (made == f.seq && code.size() < f.seq.size()) { f.code = std::move(code); f.literal_bp = lit_bp; }
                            else f.literal_bp = f.seq.size();
                            return true;
                        };
                        copy_split[ai] = split_copies(graph, Rdiff, Htok, events, module_rec->module_nodes, module_rec->anchor_nodes,
                                                      static_cast<long long>(cit->second) - static_cast<long long>(module_rec->seed.ref_cn),
                                                      bit->second, module_rec->seed.ru_len, place, h_first, h_last);
                    }
                    // Pieces of one change are joined once the copies are decided, on the parts left, so joining
                    // never changes which parts are copies.
                    for (std::vector<Event>* evs : {&events, &copy_split[ai].rest}) {
                        join_copy_pieces(graph, Rdiff, Htok, module_rec->module_nodes, *evs, options.merge_distance_bp, options.min_sv_bp);
                        join_small_pieces(graph, Rdiff, Htok, module_rec->module_nodes, *evs, kSmallPieceGapBp, options.min_sv_bp);
                    }
                }
                // An inversion is a node-level finding (the walk turns round); such an allele keeps it.
                bool inv = false;
                for (const Event& e : events) inv = inv || e.type == EvType::Inv;
                std::vector<Event> sev;
                if (seq_ok && !inv && seq_events(ai, Htok, events, skipped > 0, sev)) {
                    events = std::move(sev);
                    allele_unaligned[ai] = 0;
                    allele_seq[ai] = 1;
                    allele_skipped_blocks[ai] = 0;
                }
            }
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
            // The events as written: placed on the reference and coalesced. A carrier's events with its copies
            // taken out go the same way.
            const auto finish = [&](std::vector<Event>& events) {
                // An insertion next to a counted loop array anchors on the reference step before the
                // array in walk order. In a forward bubble that is the array's genomic start, the
                // left-aligned place for content inside a tandem array and where the UNIT record sits; in
                // a reverse bubble it is the array's end, so the same insertion would move by the array's
                // length with the bubble's orientation. Anchor it on the far side of the array instead,
                // which is its genomic start there too.
                if (ref_rev) {
                    for (Event& e : events) {
                        if (e.type != EvType::Ins || e.ref_tok_anchor == SIZE_MAX) continue;
                        const auto it = rdiff_pos.find(e.ref_tok_anchor);
                        if (it == rdiff_pos.end() || it->second + 1 >= Rdiff.size()) continue;
                        const std::size_t next = Rdiff[it->second + 1].src_idx;
                        if (next <= e.ref_tok_anchor + 1) continue;   // nothing dropped in between
                        e.ref_tok_anchor = next - 1;
                        e.anchor_node = ref_steps[next - 1].node_id;
                    }
                }
                // The token->genomic resolver, so nearness is judged at the occurrence the events sit at.
                // Aligned events need none: the alignment already made each indel one gap.
                if (events.empty() || !events.front().seq_event) coalesce_events(graph, events, ref_node_pos, hap_node_pos, options.merge_distance_bp,
                                [&](std::size_t t0, std::size_t t1, long long& lo, long long& hi) {
                                    const auto sp = tok_span(t0, t1);
                                    if (sp.first == 0) return false;
                                    lo = static_cast<long long>(sp.first);
                                    hi = static_cast<long long>(sp.second) + 1;   // half-open
                                    return true;
                                });
                for (Event& e : events) {
                    if (e.seq_event) continue;   // already placed, left-aligned
                    const long long p = ev_ref_pos(e);
                    e.ref_pos = p < 0 ? 0 : static_cast<std::size_t>(p);
                    // One exact change: a deletion or substitution whose reference tokens are exactly the
                    // bases it removes (not pieces with reference kept between), or an insertion at one point.
                    if (e.type == EvType::Del || e.type == EvType::Sub) {
                        if (e.ref_tok_first == SIZE_MAX || e.ref_tok_last == SIZE_MAX) continue;
                        const auto sp = tok_span(e.ref_tok_first, e.ref_tok_last);
                        const std::size_t del = e.type == EvType::Sub ? e.del_bp : e.size_bp;
                        if (sp.first <= ref_meta.region_start_1based || sp.second + 1 - sp.first != del) continue;
                        e.splice_before = sp.first - 1;
                        e.splice_del = del;
                        if (e.type == EvType::Sub) e.splice_ins = ref_rev ? reverse_complement(e.alt_seq) : e.alt_seq;
                        for (char& c : e.splice_ins) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
                        e.splice_ok = !ref_genome.empty();
                        if (e.type == EvType::Sub && e.splice_ok && !e.across_copies && !e.split_piece && !module_record) plain_indel(e);
                    } else if (e.type == EvType::Ins && e.ref_tok_anchor != SIZE_MAX) {
                        const auto sp = tok_span(e.ref_tok_anchor, e.ref_tok_anchor);
                        if (sp.first <= ref_meta.region_start_1based) continue;
                        e.splice_before = ref_rev ? sp.first - 1 : sp.second;
                        e.splice_ins = ref_rev ? reverse_complement(e.seq) : e.seq;
                        for (char& c : e.splice_ins) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
                        e.splice_ok = !ref_genome.empty();
                    }
                }
            };
            finish(events);
            if (copy_split[ai].mode == CopySplit::Copies) finish(copy_split[ai].rest);
            allele_events[ai] = std::move(events);
        });
        std::unordered_set<std::string> unaligned_haps;
        for (std::size_t ai = 0; ai < alleles.size(); ++ai) {
            if (allele_unaligned[ai]) unaligned_haps.insert(alleles[ai].members.begin(), alleles[ai].members.end());
            summary.skipped_large_segments += allele_skipped_blocks[ai];
        }
        if (!unaligned_haps.empty()) summary.unaligned_bubbles.push_back(bubble.id);
        // A module record's carriers are the haplotypes whose copy change is copies; one whose change is a larger
        // event of its own leaves it (GT 0, its CN set where the record is written). The carriers' copies go to
        // the record, as its forms, and their other events on to the merge -- only when the record holds with
        // its carriers left: tested on the lowest AF it can be written with, so the filters run at writing
        // cannot drop it and the copies with it. Otherwise no carrier's copies are taken from its events.
        if (module_rec != nullptr) {
            MergedRecord& d = merged[module_idx];
            for (std::size_t ai = 0; ai < alleles.size(); ++ai)
                if (cnbp_fix[ai] != 0)
                    for (const std::string& m : alleles[ai].members)
                        if (const auto it = d.sample_copies_bp.find(m); it != d.sample_copies_bp.end()) it->second += cnbp_fix[ai];
            std::unordered_map<std::string, std::size_t> allele_of;
            for (std::size_t ai = 0; ai < alleles.size(); ++ai)
                for (const std::string& m : alleles[ai].members) allele_of.emplace(m, ai);
            MergedRecord probe;
            probe.seed.type = EvType::Dup;
            probe.sample_cn = d.sample_cn;
            for (const std::string& c : d.carriers) {
                const auto it = allele_of.find(c);
                if (it != allele_of.end() && copy_split[it->second].mode == CopySplit::OwnEvent) { ++summary.copy_own_event; continue; }
                probe.carriers.push_back(c);
            }
            if (!probe.carriers.empty() && support_ok(probe, true)) {
                d.carriers = probe.carriers;
                std::map<std::tuple<std::size_t, std::size_t, std::string>, std::size_t> form_idx;
                for (std::size_t ai = 0; ai < alleles.size(); ++ai) {
                    CopySplit& cs = copy_split[ai];
                    if (cs.mode != CopySplit::Copies) continue;
                    summary.copy_carriers += alleles[ai].members.size();
                    if (cs.forms.empty()) summary.copy_no_form += alleles[ai].members.size();
                    const auto cit = d.sample_cn.find(alleles[ai].members.front());
                    const long long delta = cit == d.sample_cn.end() ? 0 :
                        static_cast<long long>(cit->second) - static_cast<long long>(d.seed.ref_cn);
                    if (!cs.forms.empty() && cs.passes != delta) summary.copy_pass_off += alleles[ai].members.size();
                    std::vector<std::size_t>& fi = d.allele_forms[ai];
                    for (CopyForm& f : cs.forms) {
                        const auto [it, fresh] = form_idx.try_emplace({f.start, f.end, f.seq}, d.forms.size());
                        if (fresh) d.forms.push_back(std::move(f));
                        fi.push_back(it->second);
                    }
                    allele_events[ai] = std::move(cs.rest);
                }
                // forms in reference order, so CN_FORMS does not depend on the order alleles came in
                std::vector<std::size_t> order(d.forms.size()), rank(d.forms.size());
                for (std::size_t k = 0; k < order.size(); ++k) order[k] = k;
                std::sort(order.begin(), order.end(), [&](std::size_t x, std::size_t y) {
                    return std::tie(d.forms[x].start, d.forms[x].end, d.forms[x].seq) <
                           std::tie(d.forms[y].start, d.forms[y].end, d.forms[y].seq); });
                std::vector<CopyForm> sorted;
                for (std::size_t k = 0; k < order.size(); ++k) { rank[order[k]] = k; sorted.push_back(std::move(d.forms[order[k]])); }
                d.forms = std::move(sorted);
                for (auto& [ai, fi] : d.allele_forms) {
                    for (std::size_t& k : fi) k = rank[k];
                    std::sort(fi.begin(), fi.end());
                }
            } else {
                // no carrier's copies leave its events: the record is not written
                dropped_modules.push_back(std::move(d));
                merged.erase(merged.begin() + static_cast<std::ptrdiff_t>(module_idx));
                module_rec = nullptr;
            }
        }
        // The reference bases each allele deletes, [first, last] genomic, from its exact deletions and
        // substitutions: a haplotype whose own larger deletion takes in a record's bases has no allele at that
        // record, and is written GT '.' there rather than 0 (see covered_here). Not at a collapsed module,
        // whose copies the walks still hold: losing one copy is not losing the bases the others carry.
        auto event_ref_iv = [&](const Event& e) -> std::pair<std::size_t, std::size_t> {
            if (e.seq_event) {
                if (e.del_bp == 0) return {0, 0};
                return {e.seq_before + 1, e.seq_before + e.del_bp};
            }
            if (e.ref_tok_first == SIZE_MAX || e.ref_tok_last == SIZE_MAX) return {0, 0};
            const auto sp = tok_span(e.ref_tok_first, e.ref_tok_last);
            if (sp.first == 0) return {0, 0};
            const std::size_t del = e.type == EvType::Sub ? e.del_bp : e.size_bp;
            if (sp.second + 1 - sp.first != del) return {0, 0};   // pieces with reference kept between them
            return sp;
        };
        std::vector<std::vector<std::pair<std::size_t, std::size_t>>> allele_deleted(alleles.size());
        if (!module_record)
            for (std::size_t ai = 0; ai < alleles.size(); ++ai)
                for (const Event& e : allele_events[ai])
                    if (e.type == EvType::Del || e.type == EvType::Sub) {
                        const auto iv = event_ref_iv(e);
                        if (iv.first > 0) allele_deleted[ai].push_back(iv);
                    }
        std::unordered_map<std::string, std::size_t> member_allele;
        for (std::size_t ai = 0; ai < alleles.size(); ++ai)
            for (const std::string& m : alleles[ai].members) member_allele.emplace(m, ai);

        // Cross-haplotype merge by transitive single-linkage: connected components over the
        // events_match graph, so A~B~C collapse to one record even when A!~C (greedy first-fit would
        // fragment them). Components holding aligned events are then regrouped around representatives.
        // Content classes over every allele event: events equal in type, ref_pos, size_bp, nodes, seq and
        // alt_seq are taken to match the same events, so each content is compared once; most haplotype
        // events are such copies, and this saves the bulk of the alignment work.
        std::vector<std::vector<std::size_t>> ev_class(alleles.size());   // [ai][k] -> class
        std::vector<const Event*> class_rep;
        {
            auto same_content = [](const Event& a, const Event& b) {
                return a.type == b.type && a.ref_pos == b.ref_pos && a.size_bp == b.size_bp &&
                       a.seq == b.seq && a.alt_seq == b.alt_seq && a.nodes == b.nodes;
            };
            std::unordered_map<std::uint64_t, std::vector<std::size_t>> buckets;
            for (std::size_t ai = 0; ai < alleles.size(); ++ai) {
                ev_class[ai].reserve(allele_events[ai].size());
                for (const Event& e : allele_events[ai]) {
                    std::uint64_t h = std::hash<std::string>{}(e.seq) ^ (std::hash<std::string>{}(e.alt_seq) << 1);
                    auto mix = [&h](std::uint64_t v) { h ^= v + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2); };
                    mix(static_cast<std::uint64_t>(e.type));
                    mix(e.ref_pos);
                    mix(e.size_bp);
                    for (const std::string& n : e.nodes) mix(std::hash<std::string>{}(n));
                    std::vector<std::size_t>& b = buckets[h];
                    std::size_t c = SIZE_MAX;
                    for (const std::size_t k : b) if (same_content(*class_rep[k], e)) { c = k; break; }
                    if (c == SIZE_MAX) { c = class_rep.size(); class_rep.push_back(&e); b.push_back(c); }
                    ev_class[ai].push_back(c);
                }
            }
        }
        {
            struct Cand { std::size_t ai; const Event* e; std::size_t cls; };
            std::vector<Cand> cands;
            for (std::size_t ai = 0; ai < alleles.size(); ++ai)
                for (std::size_t k = 0; k < allele_events[ai].size(); ++k) {
                    const Event& e = allele_events[ai][k];
                    if (e.size_bp >= rescue_floor) cands.push_back({ai, &e, ev_class[ai][k]});
                }
            DisjointSet dsu(cands.size());
            // Position-sorted windowed sweep keeps edge-building near-linear: events_match can only
            // fire within |dpos| <= merge_distance_bp + size_i, so we stop once past that bound.
            // Zero-coordinate events (anchor missing from the reference) sort first and are not
            // position-pruned against each other, but a nonzero event's scan stops at its first
            // out-of-window predecessor and so rarely reaches them: in effect they are compared only
            // among themselves. Rare in practice: ev_ref_pos falls back to the anchor node, a reference
            // node of the bubble.
            // Swept over one candidate per content class: two members of a class always match (same
            // nodes, Jaccard 1), and whether two classes match depends only on their contents.
            std::vector<std::size_t> ord;
            {
                std::unordered_map<std::size_t, std::size_t> first_of_class;
                for (std::size_t i = 0; i < cands.size(); ++i) {
                    const auto [it, fresh] = first_of_class.try_emplace(cands[i].cls, i);
                    if (fresh) ord.push_back(i);
                    else dsu.unite(it->second, i);
                }
            }
            std::sort(ord.begin(), ord.end(), [&](std::size_t a, std::size_t b) {
                return cands[a].e->ref_pos < cands[b].e->ref_pos;
            });
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
                    if (events_match(ei, ej, nullptr, nullptr)) dsu.unite(i, j);
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
            // Events found by aligning alleles merge around a representative instead, as truvari's
            // collapse does: contents in order of how many haplotypes carry them, each joining the first
            // representative it matches, or becoming one. Single linkage would chain the length alleles of
            // a repeat, each close to the next, into one record whose members share nothing with its
            // representative. Node events keep single linkage, except at a collapsed module: there every
            // copy shares the same nodes, so single linkage chains different copies' events, a haplotype's
            // own among them, into one record. There a record also holds one event per haplotype, as the
            // VCF holds one allele per haplotype per record.
            {
                std::vector<std::vector<std::size_t>> split;
                for (std::vector<std::size_t>& comp : comps) {
                    bool aligned = false;
                    for (const std::size_t k : comp) aligned = aligned || cands[k].e->seq_event;
                    if (comp.size() < 2 || (!aligned && !module_record)) { split.push_back(std::move(comp)); continue; }
                    std::map<std::size_t, std::vector<std::size_t>> by_cls;
                    for (const std::size_t k : comp) by_cls[cands[k].cls].push_back(k);
                    struct ClsInfo { std::size_t cls, haps; };
                    std::vector<ClsInfo> order;
                    for (const auto& [cls, ks] : by_cls) {
                        std::unordered_set<std::size_t> ais;
                        std::size_t h = 0;
                        for (const std::size_t k : ks) if (ais.insert(cands[k].ai).second) h += alleles[cands[k].ai].members.size();
                        order.push_back({cls, h});
                    }
                    std::sort(order.begin(), order.end(), [&](const ClsInfo& x, const ClsInfo& y) {
                        if (x.haps != y.haps) return x.haps > y.haps;
                        const Event& a = *class_rep[x.cls];
                        const Event& b = *class_rep[y.cls];
                        if (a.size_bp != b.size_bp) return a.size_bp > b.size_bp;
                        if (module_record)   // the same content on other copies' nodes: their order fixed too
                            return std::tie(a.ref_pos, a.seq, a.alt_seq, a.nodes) < std::tie(b.ref_pos, b.seq, b.alt_seq, b.nodes);
                        return std::tie(a.ref_pos, a.seq, a.alt_seq) < std::tie(b.ref_pos, b.seq, b.alt_seq);
                    });
                    std::vector<std::size_t> seeds;          // class of each group's representative
                    std::vector<std::vector<std::size_t>> groups;
                    std::vector<std::unordered_set<std::size_t>> group_alleles;   // module site: the alleles in each
                    for (const ClsInfo& ci : order) {
                        const Event& ev = *class_rep[ci.cls];
                        if (module_record) {
                            // Each candidate on its own: one whose allele a group holds already goes on to the next
                            // representative it matches, or becomes one.
                            std::vector<signed char> ok(seeds.size(), -1);   // the class against each seed, once
                            for (const std::size_t k : by_cls.at(ci.cls)) {
                                std::size_t g = groups.size();
                                for (std::size_t x = 0; x < seeds.size() && g == groups.size(); ++x) {
                                    if (group_alleles[x].count(cands[k].ai)) continue;
                                    if (ok[x] < 0) {
                                        const Event& sd = *class_rep[seeds[x]];
                                        const bool sw = ev.seq < sd.seq || (ev.seq == sd.seq && ev.size_bp < sd.size_bp);
                                        ok[x] = (sw ? events_match(ev, sd, nullptr, nullptr) : events_match(sd, ev, nullptr, nullptr)) ? 1 : 0;
                                    }
                                    if (ok[x] == 1) g = x;
                                }
                                if (g == groups.size()) {
                                    seeds.push_back(ci.cls); groups.emplace_back(); group_alleles.emplace_back(); ok.push_back(1);
                                }
                                groups[g].push_back(k);
                                group_alleles[g].insert(cands[k].ai);
                            }
                            continue;
                        }
                        std::size_t g = groups.size();
                        for (std::size_t x = 0; x < seeds.size() && g == groups.size(); ++x) {
                            const Event& sd = *class_rep[seeds[x]];
                            const bool sw = ev.seq < sd.seq || (ev.seq == sd.seq && ev.size_bp < sd.size_bp);
                            if (sw ? events_match(ev, sd, nullptr, nullptr) : events_match(sd, ev, nullptr, nullptr)) g = x;
                        }
                        if (g == groups.size()) { seeds.push_back(ci.cls); groups.emplace_back(); }
                        const std::vector<std::size_t>& ks = by_cls.at(ci.cls);
                        groups[g].insert(groups[g].end(), ks.begin(), ks.end());
                    }
                    for (std::vector<std::size_t>& g : groups) {
                        std::sort(g.begin(), g.end());
                        split.push_back(std::move(g));
                    }
                }
                comps = std::move(split);
            }
            // The representative is the member event the most haplotypes carry as written, then the larger,
            // then by content: a record then reports what most of its carriers carry rather than one
            // outlier's size or sequence, and does not move when the GFA's paths are reordered.
            const std::size_t first_new_record = merged.size();
            // An event as the record would write it: its minimal change, the bases its two sides share at either
            // end trimmed and an insertion or deletion left-aligned. Two members that differ only in how the walks
            // placed them then write the same, and count as one. An event not placed on the reference by its bases
            // goes by its type, anchor and sequences.
            using WrittenKey = std::tuple<int, std::size_t, std::string, std::string>;
            const auto written_key = [&](const Event& e) -> WrittenKey {
                const std::size_t base = ref_meta.region_start_1based;
                std::size_t pos = 0;
                std::string d, a;
                if (e.seq_event && !ref_base_at(e.seq_before).empty()) {
                    pos = e.seq_before;
                    d = e.type == EvType::Ins ? std::string() : e.seq;
                    a = e.type == EvType::Ins ? e.seq : e.type == EvType::Sub ? e.alt_seq : std::string();
                    if (ref_rev) { d = reverse_complement(d); a = reverse_complement(a); }
                } else if (!e.seq_event && e.splice_ok && (e.type == EvType::Del || e.type == EvType::Ins || e.type == EvType::Sub) &&
                           e.splice_before >= base && e.splice_before + 1 - base + e.splice_del <= ref_genome.size()) {
                    pos = e.splice_before;
                    d = ref_genome.substr(e.splice_before + 1 - base, e.splice_del);
                    a = e.splice_ins;
                } else {
                    return {static_cast<int>(e.type), e.ref_pos, e.seq, e.alt_seq};
                }
                for (char& c : d) c = upper_base(c)[0];
                for (char& c : a) c = upper_base(c)[0];
                std::size_t pre = 0;
                while (pre < d.size() && pre < a.size() && d[pre] == a[pre]) ++pre;
                std::size_t suf = 0;
                while (suf < d.size() - pre && suf < a.size() - pre && d[d.size() - 1 - suf] == a[a.size() - 1 - suf]) ++suf;
                d = d.substr(pre, d.size() - pre - suf);
                a = a.substr(pre, a.size() - pre - suf);
                pos += pre;
                if (d.empty() != a.empty()) {
                    std::string& x = d.empty() ? a : d;
                    for (;;) {
                        const std::string b = ref_base_at(pos);
                        if (b.empty() || x.empty() || b[0] != x.back() || ref_base_at(pos - 1).empty()) break;
                        x = x.back() + x.substr(0, x.size() - 1);
                        --pos;
                    }
                }
                return {-1, pos, d, a};
            };
            // true: a ranks before b. The same content on different nodes (a graph that kept two copies of one
            // sequence apart) goes to the member whose first carrier comes first by name, which renumbering the
            // graph does not change; node ids decide only between two events of one allele.
            auto first_carrier = [&](std::size_t ai) {
                return *std::min_element(alleles[ai].members.begin(), alleles[ai].members.end());
            };
            auto content_better = [&](const Event& a, std::size_t a_ai, const Event& b, std::size_t b_ai) {
                if (a.ref_pos != b.ref_pos) return a.ref_pos < b.ref_pos;
                if (a.seq != b.seq) return a.seq < b.seq;
                if (a.alt_seq != b.alt_seq) return a.alt_seq < b.alt_seq;
                if (a_ai != b_ai) return first_carrier(a_ai) < first_carrier(b_ai);
                return a.nodes < b.nodes;
            };
            for (const std::vector<std::size_t>& comp : comps) {
                // haplotypes carrying each member's exact event, over the component's members
                std::vector<std::size_t> exact_haps(comp.size(), 0);
                {
                    // An "exact" member event: the same as written (written_key; nodes may differ). Grouped by
                    // that key rather than compared pairwise: each group's haplotypes are summed once over its
                    // distinct alleles.
                    std::map<WrittenKey, std::unordered_set<std::size_t>> by_event;
                    const auto& ekey = written_key;
                    for (const std::size_t k : comp) by_event[ekey(*cands[k].e)].insert(cands[k].ai);
                    auto haps = [&](const std::unordered_set<std::size_t>& ais) {
                        std::size_t n = 0;
                        for (const std::size_t ai : ais) n += alleles[ai].members.size();
                        return n;
                    };
                    for (std::size_t x = 0; x < comp.size(); ++x)
                        exact_haps[x] = haps(by_event.at(ekey(*cands[comp[x]].e)));
                }
                // A member under --min-sv-bp represents the record only when none reaches it: the record is
                // kept on its largest member, and a smaller representative would report a size under the
                // filter it passed.
                // A substitution merged with deletions or insertions is sized by its net change, as they are,
                // since the haplotypes they make are the same.
                bool comp_mixed = false;
                for (const std::size_t k : comp) comp_mixed = comp_mixed || cands[k].e->type != EvType::Sub;
                bool comp_has_sub = false;
                for (const std::size_t k : comp) comp_has_sub = comp_has_sub || cands[k].e->type == EvType::Sub;
                comp_mixed = comp_mixed && comp_has_sub;
                const auto eff_size = [&](const Event& e) -> std::size_t {
                    if (!comp_mixed || e.type != EvType::Sub) return e.size_bp;
                    const long long n = static_cast<long long>(e.ins_bp) - static_cast<long long>(e.del_bp);
                    return static_cast<std::size_t>(n < 0 ? -n : n);
                };
                const auto sv_sized = [&](const Event& e) { return eff_size(e) >= options.min_sv_bp; };
                std::size_t bx = 0;
                for (std::size_t x = 1; x < comp.size(); ++x) {
                    const Event& a = *cands[comp[x]].e;
                    const Event& b = *cands[comp[bx]].e;
                    bool better;
                    if (sv_sized(a) != sv_sized(b)) better = sv_sized(a);
                    else if (exact_haps[x] != exact_haps[bx]) better = exact_haps[x] > exact_haps[bx];
                    // an aligned member places the change at the base; a node one, at the nodes it replaced
                    else if (a.seq_event != b.seq_event) better = a.seq_event;
                    else if (a.size_bp != b.size_bp) better = a.size_bp > b.size_bp;
                    else better = content_better(a, cands[comp[x]].ai, b, cands[comp[bx]].ai);
                    if (better) bx = x;
                }
                const std::size_t best = comp[bx];
                MergedRecord mr;
                mr.seed = *cands[best].e;                 // most-carried member represents the record
                mr.min_size_bp = mr.seed.size_bp;
                mr.max_size_bp = mr.seed.size_bp;
                for (std::size_t k : comp) {
                    const Event& e = *cands[k].e;
                    mr.min_size_bp = std::min(mr.min_size_bp, e.size_bp);
                    mr.max_size_bp = std::max(mr.max_size_bp, e.size_bp);
                    { std::size_t& ms = mr.member_size[cands[k].ai]; ms = std::max(ms, e.size_bp); }
                    if (e.type == EvType::Sub)
                        mr.member_net[cands[k].ai] = static_cast<long long>(e.ins_bp) - static_cast<long long>(e.del_bp);
                    mr.filter_bp = std::max(mr.filter_bp, eff_size(e));
                    mr.member_nodes.insert(e.nodes.begin(), e.nodes.end());
                    mr.member_event[cands[k].ai].push_back(e);
                    if (mr.member_alleles.insert(cands[k].ai).second)
                        for (const std::string& m : alleles[cands[k].ai].members)
                            mr.carriers.push_back(m);
                }
                merged.push_back(std::move(mr));
            }
            // ID suffixes (_2, _3) follow emission order, so put the new records in content order first.
            std::stable_sort(merged.begin() + static_cast<std::ptrdiff_t>(first_new_record), merged.end(),
                                 [](const MergedRecord& a, const MergedRecord& b) {
                                     const Event& x = a.seed;
                                     const Event& y = b.seed;
                                     if (x.ref_pos != y.ref_pos) return x.ref_pos < y.ref_pos;
                                     if (x.size_bp != y.size_bp) return x.size_bp < y.size_bp;
                                     if (x.type != y.type) return x.type < y.type;
                                     if (x.seq != y.seq) return x.seq < y.seq;
                                     return x.nodes < y.nodes;
                                 });
        }

        // Force-call: re-read every non-member haplotype's walk-diff and add it as a carrier if it
        // supports the record's fixed representative, even below the size threshold. Monotone, so one
        // pass suffices. Not a node-set containment test: on folded loci the inserted nodes are shared
        // across nearly all walks, so it would force-call everyone.
        // At a collapsed module an event goes to one record only: one of --rescue-min-bp or more is a
        // candidate and found its record in the merge, and a smaller one goes to the first record (in
        // record order) it matches. Representatives that do not match each other can both match it, and
        // filed twice it would be counted twice.
        const auto add_member = [&](MergedRecord& mr, std::size_t ai, const Event& e) {
            mr.member_alleles.insert(ai);
            mr.member_event[ai].push_back(e);
            mr.member_size[ai] = e.size_bp;
            if (e.type == EvType::Sub)
                mr.member_net[ai] = static_cast<long long>(e.ins_bp) - static_cast<long long>(e.del_bp);
            mr.min_size_bp = std::min(mr.min_size_bp, e.size_bp);
            mr.max_size_bp = std::max(mr.max_size_bp, e.size_bp);
            mr.member_nodes.insert(e.nodes.begin(), e.nodes.end());
            for (const std::string& m : alleles[ai].members) mr.carriers.push_back(m);
        };
        if (module_record) {
            std::vector<std::vector<std::pair<std::size_t, std::vector<std::size_t>>>> hits(merged.size());
            run_parallel(merged.size(), options.threads, [&](std::size_t ri) {
                const MergedRecord& mr = merged[ri];
                if (mr.seed.type == EvType::Dup) return;
                std::vector<signed char> seen(class_rep.size(), -1);
                for (std::size_t ai = 0; ai < alleles.size(); ++ai) {
                    if (mr.member_alleles.count(ai)) continue;
                    std::vector<std::size_t> ks;
                    for (std::size_t k = 0; k < allele_events[ai].size(); ++k) {
                        const Event& e = allele_events[ai][k];
                        if (e.size_bp >= rescue_floor) continue;
                        signed char& ok = seen[ev_class[ai][k]];
                        if (ok < 0) ok = events_match(mr.seed, e) ? 1 : 0;
                        if (ok == 1) ks.push_back(k);
                    }
                    if (!ks.empty()) hits[ri].emplace_back(ai, std::move(ks));
                }
            });
            std::set<std::pair<std::size_t, std::size_t>> claimed;
            for (std::size_t ri = 0; ri < merged.size(); ++ri)
                for (const auto& [ai, ks] : hits[ri])
                    for (const std::size_t k : ks)
                        if (claimed.insert({ai, k}).second) { add_member(merged[ri], ai, allele_events[ai][k]); break; }
        } else
        run_parallel(merged.size(), options.threads, [&](std::size_t ri) {
            MergedRecord& mr = merged[ri];
            if (mr.seed.type == EvType::Dup) return;
            // The seed against each content class once (see ev_class).
            struct SeedMatch { signed char ok = -1; };
            std::vector<SeedMatch> seen(class_rep.size());
            for (std::size_t ai = 0; ai < alleles.size(); ++ai) {
                if (mr.member_alleles.count(ai)) continue;
                for (std::size_t k = 0; k < allele_events[ai].size(); ++k) {
                    const Event& e = allele_events[ai][k];
                    // An aligned event in a record already is not claimed by another: representatives
                    // that do not match each other can both match it.
                    if (mr.seed.seq_event && e.seq_event && e.size_bp >= rescue_floor) continue;
                    SeedMatch& sm = seen[ev_class[ai][k]];
                    if (sm.ok < 0) sm.ok = events_match(mr.seed, e) ? 1 : 0;
                    if (sm.ok == 1) { add_member(mr, ai, e); break; }
                }
            }
        });

        // ---- Keep records reaching min_sv_bp (by their largest member) and --min-haplotypes.
        {
            std::vector<MergedRecord> kept;
            for (MergedRecord& mr : merged) {
                if (mr.seed.type != EvType::Dup) {
                    if (mr.filter_bp < options.min_sv_bp) continue;
                }
                if (!support_ok(mr)) {
                    // a module record still gives its genes' copy numbers
                    if (!mr.anchor_nodes.empty()) dropped_modules.push_back(std::move(mr));
                    continue;
                }
                kept.push_back(std::move(mr));
            }
            merged = std::move(kept);
        }

        if (!merged.empty()) ++summary.bubbles_with_calls;

        // Map each sample to its allele (for per-carrier sub-walk provenance).
        std::unordered_map<std::string, std::size_t> sample_to_allele;
        for (std::size_t ai = 0; ai < alleles.size(); ++ai) {
            for (const std::string& m : alleles[ai].members) sample_to_allele[m] = ai;
        }

        // Genes on this bubble's interior nodes (for the per-gene DUP table, resolved by private
        // k-mer dosage after the bubbles are called). The GENES INFO field is graph-based (collapsed
        // nodes tag multiple genes); the per-gene copy number is resolved later from k-mers, which
        // separates collapsed paralogs the graph cannot. Interior only: the source and sink are the
        // flanks the site is anchored on, and a gene lying only there would get the module's copy
        // number though no copy of it varies.
        std::vector<int> bubble_gene_idx;
        if (!node_genes.empty()) {
            std::unordered_set<int> seen;
            for (const std::string& id : bubble.inside) {
                const auto it = node_genes.find(id);
                if (it == node_genes.end()) continue;
                for (int gi : it->second) if (seen.insert(gi).second) bubble_gene_idx.push_back(gi);
            }
            std::sort(bubble_gene_idx.begin(), bubble_gene_idx.end());
        }
        if (!bubble_gene_idx.empty())
            for (const MergedRecord& d : dropped_modules) {
                DupGeneTarget t{bubble.id, ".", bubble_gene_idx, d.sample_cn, d.module_nodes, {}, true};
                dup_targets.push_back(std::move(t));
            }

        if (merged.empty()) return;   // the dropped modules' genes are queued above

        // A module record's non-carriers whose change records of their own hold: their event in such a record
        // (their own, else the record's) is on the module's anchor and changes its length the way their copy
        // number does -- the test the copy accounting above applies to a carrier's events. They read REF_CN, so
        // a dosage read off CN does not count the change twice; a non-carrier whose change no record holds keeps
        // the copies its walk holds. Which of those records is written is known only once the write-time
        // filters ran, so the holders are kept with the line (OutRecord::cn_held).
        std::unordered_map<const MergedRecord*, std::unordered_map<std::string, std::vector<const MergedRecord*>>> cn_held;
        for (const MergedRecord& d : merged) {
            if (!is_module_cn(d.seed.cn_method) || d.seed.type != EvType::Dup || d.anchor_nodes.empty()) continue;
            const std::unordered_set<std::string> car(d.carriers.begin(), d.carriers.end());
            const auto net_of = [](const Event& ev) -> long long {
                if (ev.type == EvType::Ins) return static_cast<long long>(ev.size_bp);
                if (ev.type == EvType::Del) return -static_cast<long long>(ev.size_bp);
                return static_cast<long long>(ev.ins_bp) - static_cast<long long>(ev.del_bp);
            };
            for (const MergedRecord& r : merged) {
                if (&r == &d || r.seed.type == EvType::Dup) continue;
                for (const std::string& m : r.carriers) {
                    const auto cit = d.sample_cn.find(m);
                    if (car.count(m) || cit == d.sample_cn.end()) continue;
                    const long long delta = static_cast<long long>(cit->second) - static_cast<long long>(d.seed.ref_cn);
                    if (delta == 0) continue;
                    const auto ait = sample_to_allele.find(m);
                    const auto me = ait == sample_to_allele.end() ? r.member_event.end() : r.member_event.find(ait->second);
                    std::vector<const Event*> evs;
                    if (me != r.member_event.end()) for (const Event& ev : me->second) evs.push_back(&ev);
                    else evs.push_back(&r.seed);
                    long long net = 0;
                    bool anchor = false;
                    for (const Event* ev : evs) {
                        net += net_of(*ev);
                        for (const std::string& n : ev->nodes) anchor = anchor || d.anchor_nodes.count(n);
                    }
                    if (anchor && net * delta > 0) cn_held[&d][m].push_back(&r);
                }
            }
        }
        // The copy stretch of an aligned allele (see cn_body): its first and last reference base and the bases
        // that replace them; false where the alignment is not at hand.
        struct Stretch { std::size_t g_first = 0, g_last = 0; std::string alt; };
        std::unordered_map<const void*, std::vector<std::size_t>> ref_junctions;   // per record's copy-node set
        const auto aligned_stretch = [&](std::size_t ai, const std::unordered_set<std::string>& copy_nodes,
                                         Stretch& out) -> bool {
            if (!seq_ok || seq_pieces.empty() || ref_steps.size() < 3) return false;
            const std::size_t base = ref_meta.region_start_1based;
            const std::vector<AlignedChange>& ch = allele_align[ai];
            // Reference: genomic position of each compared base, and the junctions where this record's copies sit.
            const std::string& rs = seq_ref;
            const auto rgen = [&](std::size_t i) {   // genomic position of compared reference base i
                std::size_t pc = 0;
                while (pc + 1 < seq_pieces.size() && seq_pieces[pc + 1].off <= i) ++pc;
                return seq_pieces[pc].g0 + (i - seq_pieces[pc].off);
            };
            std::vector<std::size_t>& rj = ref_junctions[&copy_nodes];   // offsets in rs where a copy of this record lies
            if (rj.empty()) {
                std::vector<std::size_t> order;
                for (std::size_t k = 1; k + 1 < ref_steps.size(); ++k) order.push_back(k);
                if (ref_rev) std::reverse(order.begin(), order.end());
                std::size_t off = 0;
                for (const std::size_t k : order) {
                    if (rdiff_pos.count(k)) off += node_len(graph, ref_steps[k].node_id);
                    else if (copy_nodes.count(ref_steps[k].node_id) && (rj.empty() || rj.back() != off)) rj.push_back(off);
                }
                if (off != rs.size()) { rj.clear(); return false; }
            }
            if (rj.empty()) return false;
            // Haplotype: its interior with every step (hf) and with the counted copies left out (hs), genomic
            // orientation; for each compared base its place in hf; and the junctions where its copies sit.
            const std::vector<PathStep>& w = alleles[ai].steps;
            // a step is compared unless its unit's copies left the walks (as collapse_walk decides)
            const auto kept = [&](std::size_t k) {
                const auto u = walk_unit_of.find(w[k].node_id);
                return !walk_drop.count(u == walk_unit_of.end() ? w[k].node_id : u->second);
            };
            std::vector<std::size_t> order;
            for (std::size_t k = 1; k + 1 < w.size(); ++k) order.push_back(k);
            if (ref_rev) std::reverse(order.begin(), order.end());
            std::string hf;
            std::vector<std::size_t> hpos;   // compared hap base -> offset in hf
            std::vector<std::size_t> hj;     // offsets in hs where a copy of this record lies
            hf.reserve(rs.size() * 2);
            for (const std::size_t k : order) {
                const auto nit = graph.nodes.find(w[k].node_id);
                if (nit == graph.nodes.end()) return false;
                const std::string& nseq = nit->second.sequence;
                const std::size_t at = hf.size();
                // the step in genomic orientation: walked backwards on one strand of the two
                if (w[k].reverse != ref_rev) hf += reverse_complement(nseq);
                else hf += nseq;
                for (std::size_t q = at; q < hf.size(); ++q) hf[q] = static_cast<char>(std::toupper(static_cast<unsigned char>(hf[q])));
                if (kept(k)) for (std::size_t q = at; q < hf.size(); ++q) hpos.push_back(q);
                else if (copy_nodes.count(w[k].node_id) && (hj.empty() || hj.back() != hpos.size())) hj.push_back(hpos.size());
            }
            long long net = 0;
            for (const AlignedChange& c : ch)
                net += static_cast<long long>(c.q1 - c.q0) - static_cast<long long>(c.t1 - c.t0);
            if (static_cast<long long>(rs.size()) + net != static_cast<long long>(hpos.size())) return false;
            // The alignment as boundary pairs: between changes the two sides run base for base.
            // f(r, after): the hap boundary at reference boundary r, past an insertion there when `after`.
            const auto f = [&](std::size_t r, bool after) -> std::size_t {
                long long d = 0;
                for (const AlignedChange& c : ch) {
                    if (c.t1 < r || (c.t1 == r && (c.t0 < r || after)))
                        d += static_cast<long long>(c.q1 - c.q0) - static_cast<long long>(c.t1 - c.t0);
                    else if (c.t0 < r && r < c.t1) return SIZE_MAX;   // inside a change: no single place
                }
                return static_cast<std::size_t>(static_cast<long long>(r) + d);
            };
            // finv(h, lo): the reference boundary at hap boundary h; inside a change, its near edge on that side
            const auto finv = [&](std::size_t h, bool lo) -> std::size_t {
                long long d = 0;
                for (const AlignedChange& c : ch) {
                    const long long dc = static_cast<long long>(c.q1 - c.q0) - static_cast<long long>(c.t1 - c.t0);
                    if (h >= c.q1 && !(h == c.q1 && c.q0 == c.q1)) { d += dc; continue; }
                    if (h > c.q0) return lo ? c.t0 : c.t1;
                    break;
                }
                return static_cast<std::size_t>(static_cast<long long>(h) - d);
            };
            std::size_t rlo = rj.front(), rhi = rj.back();
            if (!hj.empty()) { rlo = std::min(rlo, finv(hj.front(), true)); rhi = std::max(rhi, finv(hj.back(), false)); }
            // never inside a change: widen to its edges
            for (const AlignedChange& c : ch) {
                if (c.t0 < rlo && rlo < c.t1) rlo = c.t0;
                if (c.t0 < rhi && rhi < c.t1) rhi = c.t1;
            }
            std::size_t hlo = f(rlo, true), hhi = f(rhi, false);
            if (hlo == SIZE_MAX || hhi == SIZE_MAX) return false;
            if (!hj.empty()) { hlo = std::min(hlo, hj.front()); hhi = std::max(hhi, hj.back()); }
            // An insertion at the stretch's first place goes in before it, and one at its last place after it
            // (or, where the reference's copies sit there, among them): outside the stretch only where that is
            // also its side of the haplotype's copies; otherwise its bases are in the stretch, in their place.
            for (const AlignedChange& c : ch) {
                if (c.t0 != c.t1 || (c.t0 != rlo && c.t0 != rhi)) continue;
                bool inside = false;
                if (c.t0 == rlo && !hj.empty() && c.q1 > hj.front()) inside = true;
                if (c.t0 == rhi && c.t0 != rlo) {
                    const bool among = std::find(rj.begin(), rj.end(), rhi) != rj.end();
                    if (among || (!hj.empty() && c.q0 < hj.back())) inside = true;
                }
                if (inside) { hlo = std::min(hlo, c.q0); hhi = std::max(hhi, c.q1); }
            }
            if (hhi < hlo) return false;
            // reference interval: the bases between compared base rlo-1 and compared base rhi (copies included)
            const std::size_t g_first = rlo == 0 ? tok_span(1, ref_steps.size() - 2).first : rgen(rlo - 1) + 1;
            const std::size_t g_last = rhi >= rs.size() ? tok_span(1, ref_steps.size() - 2).second : rgen(rhi) - 1;
            const std::size_t f_first = hlo == 0 ? 0 : hpos[hlo - 1] + 1;
            const std::size_t f_last = hhi >= hpos.size() ? hf.size() : hpos[hhi];
            if (g_first <= base || g_last + 1 < g_first || f_last < f_first) return false;
            out = {g_first, g_last, hf.substr(f_first, f_last - f_first)};
            return true;
        };
        std::unordered_map<const MergedRecord*, std::string> id_of;   // this site's record ids, as written
        std::vector<std::pair<std::size_t, std::vector<std::vector<const MergedRecord*>>>> held_rows;

        std::string ins_ref_window;   // the reference across the site, spelled once for --classify-ins
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
            // A symbolic DEL/INV is anchored on the base preceding the event, so REF carries one real
            // base and ALT is the symbol. Computed as a coordinate (first affected base minus one)
            // rather than by stepping back a node: the preceding node is only identifiable by walk
            // position, and every node->position map here holds first occurrences.
            long long del_inv_pos = -1;
            long long del_inv_end = -1;
            long long ins_pos = -1;
            long long rep_pos = -1;          // UNIT DUP: the base before the array, when resolved
            long long rep_last = -1;         // UNIT DUP: the array's last base on the reference
            std::size_t anchor_start = 0;    // genomic start of `anchor`'s occurrence, when known
            bool rep_no_ref_copy = false;    // UNIT DUP whose reference holds no copy: spans no reference
            if (e.type == EvType::Del || e.type == EvType::Inv || e.type == EvType::Sub) {
                // Preferred: the walk positions this event occupies, which name the occurrence. The
                // node-name fallback below takes each node's first occurrence, so a reference that
                // revisits an event node can anchor the record at the wrong copy.
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
                // the downstream flank, so the insertion sits just before its first base and POS, the
                // base before the insertion, is one less than that.
                const bool reverse_bubble = ref_rev;
                anchor_after = !reverse_bubble;
                // Prefer the flank's own occurrence over its first one, so an insertion after the second
                // visit to a repeated node is not placed at the first visit.
                const auto flank = tok_span(e.ref_tok_anchor, e.ref_tok_anchor);
                if (flank.first > 0) ins_pos = reverse_bubble ? static_cast<long long>(flank.first) - 1
                                                              : static_cast<long long>(flank.second);
                else if (reverse_bubble && rpos(e.anchor_node) > 0) ins_pos = rpos(e.anchor_node) - 1;
            } else {
                // DUP: anchor on the genomically upstream flank, POS = its last base (as INS does).
                // The module routes replace POS/END with the module's own span further down. The anchor
                // is the boundary occurrence the reference walk uses, not the node's first occurrence.
                const bool reverse_bubble = ref_rev;
                anchor = reverse_bubble ? bubble.sink : bubble.source;
                anchor_after = true;
                anchor_start = walk_boundary_start(reverse_bubble);
                // A UNIT record describes one node, so it is placed on that node, not on the bubble's
                // boundary: other sequence can sit between the two, and the boundary would put the record
                // upstream of its array. With a reference copy, POS is the base before the genomically first copy
                // in this bubble's reference walk. Without one, the copies are inserted next to the
                // reference occurrence they follow in a carrying walk: after it in a forward bubble,
                // before it in a reverse one. Both come from the walk's own steps; a node's first
                // occurrence can be a copy upstream of the bubble.
                if (e.cn_method == CnMethod::Unit && !e.nodes.empty()) {
                    long long up = -1;
                    for (std::size_t j = 0; j < ref_steps.size(); ++j) {
                        // A unit graph's copy starts on its first node, upstream of the anchor it is counted at.
                        const auto m = ug_member_of.find(ref_steps[j].node_id);
                        const std::string* u = m != ug_member_of.end() ? &m->second : rep_unit_of(ref_steps[j].node_id);
                        if (u == nullptr || *u != e.nodes.front()) continue;
                        const auto sp = tok_span(j, j);
                        if (sp.first > 0 && (up < 0 || static_cast<long long>(sp.first) < up))
                            up = static_cast<long long>(sp.first);
                        if (sp.first > 0) rep_last = std::max(rep_last, static_cast<long long>(sp.second));
                    }
                    if (up < 0 && ref_count.count(e.nodes.front())) up = rpos(e.nodes.front());
                    if (up >= 0) {
                        if (up > static_cast<long long>(ref_meta.region_start_1based)) rep_pos = up - 1;
                    } else {
                        rep_no_ref_copy = true;
                        const auto flank = tok_span(e.ref_tok_anchor, e.ref_tok_anchor);
                        if (flank.first > 0) {
                            rep_pos = reverse_bubble ? static_cast<long long>(flank.first) - 1
                                                     : static_cast<long long>(flank.second);
                        } else {
                            const long long a = rpos(e.anchor_node);
                            if (a >= 0) {
                                rep_pos = reverse_bubble
                                    ? a - 1 : a + static_cast<long long>(node_len(graph, e.anchor_node)) - 1;
                            }
                        }
                    }
                }
            }
            std::size_t pos = 1;
            std::string ref_base = "N";
            const auto ait = ref_node_pos.find(anchor);
            if (ait != ref_node_pos.end()) {
                const std::size_t glen = node_len(graph, anchor);
                const std::size_t astart = anchor_start > 0 ? anchor_start : ait->second;
                if (anchor_after && glen > 0) {
                    pos = astart + glen - 1;      // last base of the upstream flank (INS)
                } else {
                    pos = astart;                 // first base of the first event node
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
            bool del_inv_stepped_back = false;
            if (del_inv_pos > 0) {
                const std::string b = ref_base_at(static_cast<std::size_t>(del_inv_pos));
                if (!b.empty()) {
                    pos = static_cast<std::size_t>(del_inv_pos);
                    ref_base = b;
                    del_inv_stepped_back = true;
                }
            }
            if (ins_pos > 0) {
                const std::string b = ref_base_at(static_cast<std::size_t>(ins_pos));
                if (!b.empty()) {
                    pos = static_cast<std::size_t>(ins_pos);
                    ref_base = b;
                }
            }
            if (rep_pos > 0) {
                const std::string b = ref_base_at(static_cast<std::size_t>(rep_pos));
                if (!b.empty()) {
                    pos = static_cast<std::size_t>(rep_pos);
                    ref_base = b;
                }
            }
            std::size_t end = pos;
            long long svlen = 0;
            // A substitution is written by its net length change: one of --min-sv-bp or more makes it a
            // DEL or INS whose SVLEN is that net change, with DELSEQ (the replaced reference, POS+1..END)
            // and INSSEQ and the SUBSTITUTION flag; a smaller net change, a SUB record with the two
            // sequences as REF and ALT, as an aligner or vg would write it.
            const long long sub_net = e.type == EvType::Sub
                ? static_cast<long long>(e.ins_bp) - static_cast<long long>(e.del_bp) : 0;
            const bool sub_sized = e.type == EvType::Sub &&
                static_cast<std::size_t>(sub_net < 0 ? -sub_net : sub_net) >= options.min_sv_bp;
            const char* svt = !sub_sized ? ev_svtype(e.type) : sub_net < 0 ? "DEL" : "INS";
            // POS is the base before a DEL/INV unless none exists (the event starts at the region's
            // first base); POS is then the first affected base, and END must not add one more.
            const std::size_t no_pre_base =
                ((e.type == EvType::Del || e.type == EvType::Inv || e.type == EvType::Sub) && !del_inv_stepped_back && e.size_bp > 0) ? 1 : 0;
            if (e.type == EvType::Del) { svlen = -static_cast<long long>(e.size_bp); end = pos + e.size_bp - no_pre_base; }
            else if (e.type == EvType::Ins) { svlen = static_cast<long long>(e.size_bp); end = pos; }
            else if (e.type == EvType::Inv) { svlen = static_cast<long long>(e.size_bp); end = pos + e.size_bp - no_pre_base; }
            else if (e.type == EvType::Sub) { svlen = sub_net; end = pos + e.del_bp - no_pre_base; }
            // A module record carries its size in size_bp, one copy. Only the UNIT route can derive it from
            // a node length, because only there is the node the literal repeat unit -- at a module the
            // anchor is the bubble source, which is whatever short node happens to bound the site.
            else if (is_module_cn(e.cn_method)) {
                svlen = static_cast<long long>(e.size_bp);
                end = pos + e.size_bp;
            }
            else { // unit DUP: copies x unit length (signed: CN gain or loss vs reference)
                // The unit's length (RU_LEN): a one-node unit graph's node, or a unit graph's unit_bp. The node a
                // unit graph is keyed by is only one segment of it, so its length is not the unit's.
                const std::size_t unit =
                    e.ru_len > 0 ? e.ru_len : node_len(graph, e.nodes.empty() ? std::string() : e.nodes.front());
                svlen = (static_cast<long long>(e.alt_cn) - static_cast<long long>(e.ref_cn)) *
                        static_cast<long long>(unit);
                // END is the last base of the reference's last copy, so the record covers the whole array (one
                // unit long where the copies were not found on the walk). With no reference copy the record spans
                // no reference sequence, so END = POS, as for an INS.
                end = rep_no_ref_copy ? pos : rep_last > static_cast<long long>(pos) ? static_cast<std::size_t>(rep_last) : pos + unit;
            }
            // A coalesced record can describe a deletion in several pieces with reference retained
            // between them. That is not one exact deletion: the affected interval is wider than the
            // bases removed. Such a record is marked IMPRECISE and END describes the whole affected
            // span, while SVLEN stays the net deleted bases. END-POS == |SVLEN| therefore holds for
            // every precise record and the flag announces the records where it cannot.
            bool imprecise = false;
            long long cluster_slack = 0;   // a cluster's span beyond the bases it removes
            if (del_inv_end > 0 && (e.type == EvType::Del || e.type == EvType::Inv)) {
                // Bases from POS+1 to END; with no preceding base, POS is itself the first affected one.
                const long long span_bp =
                    del_inv_end - static_cast<long long>(pos) + static_cast<long long>(no_pre_base);
                const long long mag = svlen < 0 ? -svlen : svlen;
                if (span_bp > mag) {
                    end = static_cast<std::size_t>(del_inv_end);
                    imprecise = true;
                    cluster_slack = span_bp - mag;
                }
            }
            // A COLLAPSED_MODULE record describes the whole site, and its "unit" is only the shared part
            // of one copy, so pos + unit would name a reference interval that corresponds to nothing. END is
            // the module's own reference span -- the interval FORMAT:CNBP is measured over.
            if (is_module_cn(e.cn_method)) {
                // POS/END from the resolver that measures the module, so the record's interval and the
                // interval its CN/CNBP were taken over are the same one. The boundaries span
                // [lo, hi]; the module content is the interior, so POS is the last base of the near
                // boundary and END the base before the far one starts.
                bool amb = false;
                bool opens_on_source = true;
                const auto mb = module_ref_bounds(&amb, &opens_on_source);
                if (mb.first > 0 && mb.second > mb.first) {
                    const std::size_t near_len = node_len(graph, bubble.source);
                    const std::size_t far_len = node_len(graph, bubble.sink);
                    // Which boundary the span opens on is the resolver's own answer, not a comparison
                    // of first occurrences, which disagrees wherever a boundary recurs on both sides.
                    const std::size_t open_len = opens_on_source ? near_len : far_len;
                    const std::size_t close_len = opens_on_source ? far_len : near_len;
                    if (mb.first + open_len <= mb.second && mb.second >= close_len) {
                        pos = mb.first + (open_len > 0 ? open_len - 1 : 0);
                        end = mb.second - (close_len > 0 ? close_len - 1 : 0) - 1;
                        const std::string b = ref_base_at(pos);
                        if (!b.empty()) ref_base = b;
                    }
                }
            }
            if (end > ref_end_1based) end = ref_end_1based;  // END is a reference coordinate; never past the graph

            // Left-align an exact INS/DEL, as VCF normalisation does: inside a repeat the graph places
            // the event anywhere along it, and the same event then sits at different positions across
            // records and callers. While the anchor base equals the allele's last base, the event can
            // move one base left with the allele rotated. HOMLEN is how far it could slide in all,
            // measured from the left-aligned position (CIPOS 0,HOMLEN). A substitution is not moved, and
            // an IMPRECISE record has no single breakpoint to move.
            std::string aligned_seq;          // the allele in reference orientation, rotated if moved
            std::size_t homlen = 0;
            if ((e.type == EvType::Del || e.type == EvType::Ins) && e.cn_method == CnMethod::None &&
                !imprecise && !e.seq.empty() && !e.seq_event &&
                (e.type == EvType::Ins || del_inv_stepped_back)) {
                std::string a = ref_rev ? reverse_complement(e.seq) : e.seq;
                for (char& c : a) c = upper_base(c)[0];
                std::size_t right = 0;   // slide right from the graph's position
                {
                    std::string r = a;
                    std::size_t next = (e.type == EvType::Del ? end : pos) + 1;
                    for (;;) {
                        const std::string b = ref_base_at(next);
                        if (b.empty() || b[0] != r.front()) break;
                        r = r.substr(1) + r.front();
                        ++next; ++right;
                    }
                }
                std::size_t left = 0;
                for (;;) {
                    const std::string b = ref_base_at(pos);
                    if (b.empty() || b[0] != a.back() || ref_base_at(pos - 1).empty()) break;
                    a = a.back() + a.substr(0, a.size() - 1);
                    --pos; --end; ++left;
                }
                if (left > 0) { ref_base = ref_base_at(pos); aligned_seq = a; }
                homlen = left + right;
            }

            // An event found by aligning the alleles (seq_event) is written from its own position, the base
            // before it, as its minimal change: the sequence the two sides share at either end is trimmed,
            // as VCF normalisation does. What is left is an INS, a DEL, or -- both sides still holding
            // sequence -- a DEL/INS of the net change flagged SUBSTITUTION, or a SUB spelling both when the
            // net change is under --min-sv-bp. A replacement from the node comparison that is one exact change
            // is written the same way when its two sides share bases at either end, as different nodes spelling
            // the same bases do: a deletion with one base of the haplotype's own is then a plain deletion.
            bool node_minimal = false;
            if (!e.seq_event && e.type == EvType::Sub && e.splice_ok && !imprecise && !e.seq.empty() &&
                !e.alt_seq.empty()) {
                std::string d = ref_rev ? reverse_complement(e.seq) : e.seq;
                for (char& c : d) c = upper_base(c)[0];
                const std::string& a = e.splice_ins;
                // identical sides are no change; the node comparison drops them, and this keeps the guard here
                node_minimal = d != a && (upper_base(d.front()) == upper_base(a.front()) ||
                                          upper_base(d.back()) == upper_base(a.back()));
            }
            const std::size_t before = e.seq_event ? e.seq_before : e.splice_before;
            bool rp = false, rp_sub = false;
            std::string rp_del, rp_ins;   // reference orientation
            if ((e.seq_event || node_minimal) && !ref_base_at(before).empty()) {
                pos = before;
                const std::string& del_w = e.type == EvType::Ins ? std::string() : e.seq;
                const std::string& ins_w = e.type == EvType::Ins ? e.seq : e.type == EvType::Sub ? e.alt_seq : std::string();
                std::string d = ref_rev ? reverse_complement(del_w) : del_w;
                std::string a = ref_rev ? reverse_complement(ins_w) : ins_w;
                for (char& c : d) c = upper_base(c)[0];
                for (char& c : a) c = upper_base(c)[0];
                std::size_t pre = 0;
                while (pre < d.size() && pre < a.size() && d[pre] == a[pre]) ++pre;
                std::size_t suf = 0;
                while (suf < d.size() - pre && suf < a.size() - pre && d[d.size() - 1 - suf] == a[a.size() - 1 - suf]) ++suf;
                d = d.substr(pre, d.size() - pre - suf);
                a = a.substr(pre, a.size() - pre - suf);
                pos += pre;
                rp = true;
                svlen = static_cast<long long>(a.size()) - static_cast<long long>(d.size());
                rp_sub = !d.empty() && !a.empty() &&
                         static_cast<std::size_t>(svlen < 0 ? -svlen : svlen) < options.min_sv_bp;
                svt = rp_sub ? "SUB" : d.empty() ? "INS" : a.empty() ? "DEL" : svlen < 0 ? "DEL" : "INS";
                if (d.empty() || a.empty()) {         // left-align, and how far it could slide right
                    std::string x = d.empty() ? a : d;
                    std::size_t right = 0;
                    {
                        std::string r = x;
                        std::size_t next = pos + (d.empty() ? 0 : d.size()) + 1;
                        for (;;) {
                            const std::string b = ref_base_at(next);
                            if (b.empty() || x.empty() || b[0] != r.front()) break;
                            r = r.substr(1) + r.front();
                            ++next; ++right;
                        }
                    }
                    std::size_t left = 0;
                    for (;;) {
                        const std::string b = ref_base_at(pos);
                        if (b.empty() || x.empty() || b[0] != x.back() || ref_base_at(pos - 1).empty()) break;
                        x = x.back() + x.substr(0, x.size() - 1);
                        --pos; ++left;
                    }
                    homlen = left + right;
                    (d.empty() ? a : d) = x;
                }
                end = pos + d.size();
                ref_base = ref_base_at(pos);
                rp_del = d; rp_ins = a;
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
                    if (!bubble_own.count(nd) || !seen.insert(nd).second) continue;
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
            // variant_nodes; otherwise a module DUP anchored on a bubble boundary loses its POS node.
            std::vector<std::string> var_nodes;
            if (mr.member_nodes.empty()) {
                var_nodes = ev_nodes;
            } else {
                std::set<std::string> all;
                for (const std::string& nd : mr.member_nodes) if (bubble_own.count(nd)) all.insert(nd);
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

            // A replacement of --min-sv-bp or more net change (SUBSTITUTION) is spelled in REF and ALT as a SUB is,
            // when both sides fit --max-seq-bp: a VCF tool then reads the bases it replaces and what replaces
            // them. Above that, or IMPRECISE, it stays symbolic with DELSEQ and INSSEQ where they fit.
            const auto fits = [&options](std::size_t n) { return options.max_seq_bp == 0 || n <= options.max_seq_bp; };
            const bool explicit_sub = !imprecise &&
                (rp ? (!rp_sub && !rp_del.empty() && !rp_ins.empty() && fits(rp_del.size()) && fits(rp_ins.size()))
                    : (e.type == EvType::Sub && sub_sized && !e.seq.empty() && !e.alt_seq.empty() &&
                       fits(e.seq.size()) && fits(e.alt_seq.size())));
            std::ostringstream info;
            info << "END=" << end << ";SVTYPE=" << svt;
            if (rp ? (!rp_del.empty() && !rp_ins.empty() && !rp_sub) : sub_sized) info << ";SUBSTITUTION";
            // At a collapsed module carriers both gain and lose copies around REF_CN, so a single
            // record-level size is not merely imprecise -- it does not exist. Per-carrier CNBP carries
            // it instead. A UNIT record is the same when its carriers change in both directions: some
            // lose copies of the unit, others gain them, and the largest change describes only the
            // most extreme carrier. SVLEN is then left out there too, and SVLEN_RANGE gives the
            // smallest and largest carrier change, (CN - REF_CN) x RU_LEN. A record whose carriers all
            // gain, or all lose, keeps its SVLEN: the largest change, signed.
            long long rep_min = 0, rep_max = 0;
            bool rep_sized = false;
            if (e.type == EvType::Dup && e.cn_method == CnMethod::Unit) {
                for (const std::string& m : mr.carriers) {
                    const auto it = mr.sample_cn.find(m);
                    if (it == mr.sample_cn.end()) continue;
                    const long long d = (static_cast<long long>(it->second) - static_cast<long long>(e.ref_cn)) *
                                        static_cast<long long>(e.ru_len);
                    if (!rep_sized) { rep_min = rep_max = d; rep_sized = true; }
                    else { rep_min = std::min(rep_min, d); rep_max = std::max(rep_max, d); }
                }
            }
            const bool rep_mixed = rep_sized && rep_min < 0 && rep_max > 0;
            if (!is_module_cn(e.cn_method) && !rep_mixed) info << ";SVLEN=" << svlen;
            if (homlen > 0) {
                info << ";HOMLEN=" << homlen << ";CIPOS=0," << homlen;
                if (e.type == EvType::Del || (rp && rp_ins.empty() && !rp_del.empty())) info << ";CIEND=0," << homlen;
            }
            if (rep_sized && rep_min != rep_max) info << ";SVLEN_RANGE=" << rep_min << "," << rep_max;
            // Members' sizes by their net change, the quantity SVLEN gives: every event each remaining member
            // brought into the record (an allele can bring two), so the representative's is among them. A
            // substitution's size is its larger side, which would not compare with SVLEN.
            if (e.type == EvType::Del || e.type == EvType::Ins || e.type == EvType::Sub) {
                bool any = false;
                long long lo = 0, hi = 0;
                for (const std::size_t ai : mr.member_alleles) {
                    const auto me = mr.member_event.find(ai);
                    if (me == mr.member_event.end()) continue;
                    for (const Event& m : me->second) {
                        const long long n = m.type == EvType::Ins ? static_cast<long long>(m.size_bp)
                                          : m.type == EvType::Del ? -static_cast<long long>(m.size_bp)
                                          : static_cast<long long>(m.ins_bp) - static_cast<long long>(m.del_bp);
                        if (!any) { lo = hi = n; any = true; }
                        else { lo = std::min(lo, n); hi = std::max(hi, n); }
                    }
                }
                if (any && lo != hi) info << ";SVLEN_RANGE=" << lo << "," << hi;
            } else if (e.type != EvType::Dup && mr.min_size_bp != mr.max_size_bp) {
                info << ";SVLEN_RANGE=" << mr.min_size_bp << "," << mr.max_size_bp;
            }
            // How alike the merged events are: each remaining member event against the representative, in the test
            // that admits it (events_match), and the weakest of them. A member that joined only through others
            // (node events merge by single linkage) can fall under the test's threshold, which shows how far the
            // chain reached. Each content once.
            if (e.type != EvType::Dup && e.type != EvType::Inv && mr.member_alleles.size() > 0) {
                double min_sim = 2.0;
                const char* min_test = "";
                std::set<std::tuple<int, std::size_t, std::size_t, std::string, std::string, std::vector<std::string>>> seen;
                for (const std::size_t ai : mr.member_alleles) {
                    const auto me = mr.member_event.find(ai);
                    if (me == mr.member_event.end()) continue;
                    for (const Event& m : me->second) {
                        if (m.type == e.type && m.ref_pos == e.ref_pos && m.size_bp == e.size_bp && m.seq == e.seq &&
                            m.alt_seq == e.alt_seq && m.nodes == e.nodes) continue;
                        if (!seen.emplace(static_cast<int>(m.type), m.ref_pos, m.size_bp, m.seq, m.alt_seq, m.nodes).second) continue;
                        double jac = -1.0, sid = -1.0;
                        const char* test = "";
                        events_match(e, m, &jac, &sid, &test);
                        const std::string t = test;
                        const double v = t == "jaccard" ? jac : sid;
                        if (t.empty() || v < 0.0) continue;
                        if (v < min_sim) { min_sim = v; min_test = test; }
                    }
                }
                if (min_sim <= 1.0) {
                    char buf[32]; std::snprintf(buf, sizeof(buf), "%.4g", min_sim);
                    info << ";MERGE_MIN_SIM=" << buf << ";MERGE_TEST=" << min_test;
                }
            }
            // An unaligned haplotype that carries nothing here is unknown, not reference (see
            // allele_unaligned): it leaves AN and is written GT '.'. A carrier keeps its 1.
            // A haplotype that deleted more than this record's bases, all of them, has no allele here: GT '.'
            // and out of AN, as a spanning deletion is; written 0 it would read as carrying the reference.
            std::pair<std::size_t, std::size_t> rec_iv{0, 0};   // the record's bases; an insertion: first == last
            if (e.type == EvType::Del || e.type == EvType::Sub || e.type == EvType::Inv) {
                rec_iv = event_ref_iv(e);
                if (e.type == EvType::Inv && rec_iv.first == 0 && e.ref_tok_first != SIZE_MAX)
                    rec_iv = tok_span(e.ref_tok_first, e.ref_tok_last);
            } else if (e.type == EvType::Ins) {
                // the base before the insertion; a node insertion's from its splice, which knows the bubble's
                // orientation (its merge position is the anchor's last base, the far side in a reverse bubble)
                const std::size_t p = e.seq_event ? e.seq_before : e.splice_ok ? e.splice_before : e.ref_pos;
                rec_iv = {p, p};
            }
            // INSSEQ, DELSEQ and INVSEQ up to --max-seq-bp, so a record line stays a size VCF tools read
            const auto seq_fits = [&options](std::size_t n) { return options.max_seq_bp == 0 || n <= options.max_seq_bp; };
            auto covered_here = [&](const std::string& s) {
                if (rec_iv.first == 0 || carrier_set.count(s) != 0) return false;
                const auto it = member_allele.find(s);
                if (it == member_allele.end()) return false;
                for (const auto& [a, b] : allele_deleted[it->second]) {
                    if (e.type == EvType::Ins ? (a < rec_iv.first && b > rec_iv.second)
                                              : (a <= rec_iv.first && b >= rec_iv.second &&
                                                 b - a > rec_iv.second - rec_iv.first))
                        return true;
                }
                return false;
            };
            auto unknown_here = [&](const std::string& s) {
                if (e.type == EvType::Dup) return false;
                return (unaligned_haps.count(s) != 0 && carrier_set.count(s) == 0) || covered_here(s);
            };
            std::size_t n_unknown = 0, n_unaligned = 0;
            for (const std::string& s : sample_names)
                if (traverses.count(s) && unknown_here(s)) {
                    ++n_unknown;
                    if (unaligned_haps.count(s)) ++n_unaligned;
                }
            summary.unaligned_missing += n_unaligned;
            summary.covered_missing += n_unknown - n_unaligned;
            const std::size_t an = traverses.size() - n_unknown;
            const std::size_t ac = carrier_set.size();
            const double af = an > 0 ? static_cast<double>(ac) / static_cast<double>(an) : 0.0;
            info << ";BUBBLE_ID=" << bubble.id << ";AN=" << an << ";AC=" << ac;
            { char buf[32]; std::snprintf(buf, sizeof(buf), "%.4g", af); info << ";AF=" << buf; }
            info << ";EVENT_NODES=";
            for (std::size_t k = 0; k < ev_nodes.size(); ++k) { if (k) info << ','; info << ev_nodes[k]; }
            if (e.type == EvType::Dup) info << ";REF_CN=" << e.ref_cn;
            if (e.type == EvType::Dup && e.ru_len > 0) {
                if (e.cn_method == CnMethod::Unit) {
                    info << ";RU_LEN=" << e.ru_len;
                    const auto tw = rep_twins.find(e.nodes.front());
                    if (tw != rep_twins.end()) {
                        info << ";UNIT_TWIN_NODES=";
                        for (std::size_t k = 0; k < tw->second.size(); ++k) info << (k ? "," : "") << tw->second[k];
                    }
                }
                else if (e.cn_method == CnMethod::Module) info << ";CN_UNIT_BP=" << e.ru_len;
            }
            // Which measurement produced CN: an exact repeat-unit traversal count (UNIT, counted at the unit
            // graph's anchor) or a module's dosage (MODULE, passes through its anchor), whose copies can hold
            // several paralogs, so one copy's length is no carrier's size.
            if (e.type == EvType::Dup && e.cn_method != CnMethod::None) {
                info << ";CN_METHOD=" << cn_method_name(e.cn_method);
                if (e.cn_method == CnMethod::Module)
                    info << ";CN_ANCHOR_BP=" << e.anchor_bp << ";CN_MODULE_NODES_BP=" << e.module_nodes_bp;
            }
            // A breakpoint against an N run in the reference (unresolved bases, an assembly gap) is not known to
            // the base: the event may begin anywhere in the run, so the record is IMPRECISE.
            if (!imprecise && e.type != EvType::Dup && (ref_base == "N" || ref_base_at(end + 1) == "N"))
                imprecise = true;
            if (imprecise) {
                info << ";IMPRECISE";
                // How far each breakpoint can lie from where it is written. A cluster's |SVLEN| removed bases lie
                // somewhere in POS+1..END: its start up to the slack later, its end up to the slack earlier. A
                // breakpoint at an N run can lie anywhere in the run.
                if (homlen == 0) {
                    long long n_before = 0, n_after = 0;
                    for (std::size_t q = pos; q > 0 && ref_base_at(q) == "N"; --q) ++n_before;
                    for (std::size_t q = end + 1; ref_base_at(q) == "N"; ++q) ++n_after;
                    const bool has_end = e.type == EvType::Del || e.type == EvType::Inv || e.type == EvType::Sub;
                    if (!has_end) info << ";CIPOS=" << -n_before << "," << n_after;
                    else info << ";CIPOS=" << -n_before << "," << cluster_slack << ";CIEND=" << -cluster_slack << "," << n_after;
                }
            }
            if (e.module_span_ambiguous > 0)
                info << ";CN_SPAN_AMBIGUOUS=" << e.module_span_ambiguous;
            if (multi_crossing_paths > 0) info << ";MULTI_CROSSING=" << multi_crossing_paths;
                if (multi_crossing_unscored > 0) info << ";MULTI_CROSSING_UNSCORED=" << multi_crossing_unscored;
            if (!genes.empty()) {
                // Genes the variant touches: for a copy-number record the genes on the nodes it counts --
                // a MODULE's copies, a UNIT's repeat (not every gene in the site, which would name a repeat
                // between two genes as in both); for any other record the genes over its reference interval,
                // the bases it replaces or the two an insertion falls between. Not the nodes, which a
                // collapsed graph shares between places and an insertion does not have.
                std::vector<std::string> gnames;
                if (e.type == EvType::Dup && !mr.module_nodes.empty())
                    gnames = genes_for_nodes(std::vector<std::string>(mr.module_nodes.begin(), mr.module_nodes.end()));
                else if (e.type == EvType::Dup) {
                    std::unordered_set<std::string> keys;
                    for (const std::string& nd : ev_nodes) if (nd != bubble.source && nd != bubble.sink) keys.insert(nd);
                    std::vector<std::string> unit;
                    for (const std::string& nd : bubble.inside) {
                        const auto m = ug_member_of.find(nd);
                        const std::string* u = m != ug_member_of.end() ? &m->second : rep_unit_of(nd);
                        if (keys.count(nd) || (u != nullptr && keys.count(*u))) unit.push_back(nd);
                    }
                    gnames = genes_for_nodes(unit);
                } else if (e.type == EvType::Ins || end <= pos) {
                    gnames = genes_in_interval(pos, pos + 1);
                } else {
                    gnames = genes_in_interval(no_pre_base ? pos : pos + 1, end);
                }
                if (!gnames.empty()) {
                    info << ";GENES=";
                    for (std::size_t k = 0; k < gnames.size(); ++k) { if (k) info << ','; info << vcf_escape(gnames[k]); }
                }
            }
            // --classify-ins: every INS record by the bases it writes (INSSEQ, also above --max-seq-bp), a
            // replacement's inserted side included, against the reference across the site (bounded minimap2
            // calls, one at a time). Both are spelled along the site's walk, so a reverse-strand hit means the
            // inserted copy is the reverse complement of the sequence it copies.
            bool invdup = e.inverted_dup;
            if (options.classify_ins && std::string(svt) == "INS") {
                std::string w = rp ? rp_ins
                              : e.type == EvType::Ins ? (!aligned_seq.empty() ? aligned_seq : ref_rev ? reverse_complement(e.seq) : e.seq)
                              : e.type == EvType::Sub ? (ref_rev ? reverse_complement(e.alt_seq) : e.alt_seq) : std::string();
                if (ref_rev) w = reverse_complement(w);
                if (!w.empty()) {
                    static std::mutex minimap_mtx;
                    std::lock_guard<std::mutex> minimap_lock(minimap_mtx);
                    if (ins_ref_window.empty()) ins_ref_window = spell_path_steps_sequence(graph, ref_steps);
                    const Minimap2Hit hit = minimap2_best_hit("ins", w, "refwin", ins_ref_window, options.minimap_preset,
                                                              options.minimap_best_n);
                    const bool dup = hit.ok && hit.identity() >= options.ins_dup_min_identity &&
                                     hit.query_end_bp - hit.query_start_bp >= w.size() / 2;
                    info << ";INS_SUBTYPE=" << (dup ? "DUP" : "NOVEL");
                    invdup = invdup || (dup && hit.reverse);
                }
            }
            if (invdup) info << ";INVDUP";
            if (rp_sub || explicit_sub) {
                // REF/ALT spell it
            } else if (rp) {
                if (!rp_del.empty() && seq_fits(rp_del.size())) info << ";DELSEQ=" << rp_del;
                if (!rp_ins.empty() && seq_fits(rp_ins.size())) info << ";INSSEQ=" << rp_ins;
            } else if (!e.seq.empty() && seq_fits(e.seq.size())) {
                // Event sequences are spelled along the canonical source->sink walk, and stay that way
                // internally (merging compares them within one bubble). Where the reference crosses
                // this bubble sink->source they are written reverse-complemented, in reference-forward
                // orientation.
                const std::string oseq = !aligned_seq.empty() ? aligned_seq
                                         : ref_rev ? reverse_complement(e.seq) : e.seq;
                if (e.type == EvType::Ins) info << ";INSSEQ=" << oseq;
                else if (e.type == EvType::Del) info << ";DELSEQ=" << oseq;
                else if (e.type == EvType::Sub && sub_sized) info << ";DELSEQ=" << oseq;
                else if (e.type == EvType::Inv) info << ";INVSEQ=" << oseq;
            }
            // A substitution's INSSEQ is what reconstructs the haplotype, so it has its own size check
            // rather than the replaced sequence's.
            if (!rp && !explicit_sub && e.type == EvType::Sub && sub_sized && !e.alt_seq.empty() && seq_fits(e.alt_seq.size()))
                info << ";INSSEQ=" << (ref_rev ? reverse_complement(e.alt_seq) : e.alt_seq);

            // FORMAT:CNBP (DUP only): the linear bp a haplotype gains/loses through the module -- spelled
            // walk length (node length x traversal multiplicity) minus the reference's. Unlike RU_LEN
            // (one folded copy) this is the real SV size, since the walk still traverses shared nodes its
            // own number of times. Uses the widest source..sink span, not the repeat-collapsed one. A
            // MODULE record instead carries each path's change over its copies (sample_copies_bp).
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
            // The reference's total bp across this bubble. At a collapsed paralog module it holds the
            // paralog-private content that travels with each copy but that the reference does not revisit,
            // which one copy's shared part (CN_UNIT_BP) leaves out.
            if (is_dup && cnbp_ref_bp > 0) info << ";CN_MODULE_REF_BP=" << cnbp_ref_bp;
            const bool with_forms = is_dup && is_module_cn(e.cn_method) && !mr.forms.empty();
            if (with_forms) {
                info << ";CN_FORMS=";
                for (std::size_t k = 0; k < mr.forms.size(); ++k) {
                    const CopyForm& f = mr.forms[k];
                    info << (k ? "," : "") << f.start << '-' << f.end << ':'
                         << (!seq_fits(f.literal_bp) ? std::string("*") : f.code.empty() ? f.seq : f.code);
                }
            }
            auto sample_cnbp = [&](const std::string& s) -> long long {
                if (!mr.sample_copies_bp.empty()) {
                    const auto it = mr.sample_copies_bp.find(s);
                    return it == mr.sample_copies_bp.end() ? 0 : it->second;
                }
                const auto pit = name_to_pi.find(s);
                return pit == name_to_pi.end() ? 0 : cnbp_walk_bp(pit->second) - cnbp_ref_bp;
            };

            // REF/ALT: symbolic, except for a SUB, which spells the replaced and the replacing sequence
            // after the shared base at POS (none when the event starts the region).
            std::string ref_col = ref_base, alt_col = std::string("<") + svt + ">";
            if (is_dup && e.cn_method == CnMethod::Unit) alt_col = "<CNV:TR>";   // VCF 4.4 tandem-repeat CNV
            if (rp_sub || (explicit_sub && rp)) {
                ref_col = ref_base + rp_del;
                alt_col = ref_base + rp_ins;
            } else if (e.type == EvType::Sub && (!sub_sized || explicit_sub) && !rp) {
                const std::string pad = no_pre_base ? std::string() : ref_base;
                ref_col = pad + (ref_rev ? reverse_complement(e.seq) : e.seq);
                alt_col = pad + (ref_rev ? reverse_complement(e.alt_seq) : e.alt_seq);
                for (char& c : ref_col) c = upper_base(c)[0];
                for (char& c : alt_col) c = upper_base(c)[0];
            }
            std::ostringstream row;
            row << ref_meta.chrom << '\t' << pos << '\t' << id << '\t' << ref_col
                << '\t' << alt_col << "\t.\tPASS\t" << info.str()
                << (is_dup ? (e.cn_method == CnMethod::Unit && e.ru_len > 0
                                  ? "\tGT:CN:CNBP:CNRESID:CNINV"
                                  : e.cn_method == CnMethod::Module && e.ru_len > 0
                                      ? "\tGT:CN:CNBP:CNRESID:MODBP"
                                      : is_module_cn(e.cn_method) ? "\tGT:CN:CNBP:MODBP" : "\tGT:CN:CNBP")
                           : "\tGT:CN")
                << (with_forms ? ":CNFORM" : "");
            std::vector<OutRecord::CnHeld> held_fix;
            std::vector<std::vector<const MergedRecord*>> held_ptrs;
            for (const std::string& s : sample_names) {
                row << '\t';
                if (!traverses.count(s)) {
                    row << (is_dup ? ".:.:." : ".:.");
                    continue;
                }
                if (unknown_here(s)) { row << ".:."; continue; }
                const bool carrier = carrier_set.count(s) != 0;
                row << (carrier ? "1" : "0") << ':';
                if (is_dup) {
                    // A traverser without a measured CN is a bug in whichever route produced this
                    // record, not a reference-like sample: filling in REF_CN would make a complete
                    // module loss indistinguishable from carrying the reference count.
                    const auto cit = mr.sample_cn.find(s);
                    if (cit == mr.sample_cn.end()) {
                        throw std::runtime_error(
                            "panvar call: no copy number for " + s + " at " + id +
                            ", which traverses the site (" + cn_method_name(e.cn_method) +
                            " route); a traversing haplotype must always have a measured CN");
                    }
                    // A module record's non-carrier whose change a kept record of its own holds (cn_held) reads
                    // the reference's count, so a dosage read off CN does not count it twice. CNBP stays as measured.
                    const auto held = cn_held.find(&mr);
                    const auto hs_it = held == cn_held.end() ? decltype(held->second.end())() : held->second.find(s);
                    const bool held_here = is_module_cn(e.cn_method) && !carrier && held != cn_held.end() &&
                                           hs_it != held->second.end();
                    const std::size_t cn_out = held_here ? e.ref_cn : cit->second;
                    // MODBP: the module's own change, what decides GT at a MODULE record
                    const auto mbit = mr.sample_bp.find(s);
                    const std::string modbp = mbit == mr.sample_bp.end() ? std::string(".") : std::to_string(mbit->second);
                    if (held_here) {
                        std::ostringstream m;
                        m << "0:" << cit->second << ':' << sample_cnbp(s);
                        if (e.ru_len > 0)
                            m << ':' << (sample_cnbp(s) - (static_cast<long long>(cit->second) - static_cast<long long>(e.ref_cn)) *
                                                          static_cast<long long>(e.ru_len));
                        m << ':' << modbp;
                        if (with_forms) m << ":.";
                        OutRecord::CnHeld h;
                        h.col = 9 + static_cast<std::size_t>(&s - sample_names.data());
                        h.measured = m.str();
                        held_fix.push_back(std::move(h));
                        held_ptrs.push_back(hs_it->second);
                    }
                    row << cn_out;
                    const long long this_cnbp = sample_cnbp(s);
                    row << ':' << this_cnbp;
                    // CNRESID: the change in bases that whole copies do not explain, (CN-REF_CN) copies
                    // of RU_LEN (UNIT) or of CN_UNIT_BP, the reference's copy (MODULE), from the measured CN
                    // (a held non-carrier's CN is written as REF_CN, but its copies still explain its bases).
                    // On a module it also holds the size difference of copies of another type than the reference's.
                    if ((e.cn_method == CnMethod::Unit || e.cn_method == CnMethod::Module) && e.ru_len > 0) {
                        const long long predicted =
                            (static_cast<long long>(cit->second) - static_cast<long long>(e.ref_cn)) *
                            static_cast<long long>(e.ru_len);
                        row << ':' << (this_cnbp - predicted);
                    }
                    if (is_module_cn(e.cn_method)) row << ':' << modbp;
                    if (with_forms) {
                        std::string fs;
                        const auto ait = sample_to_allele.find(s);
                        const auto fit = !carrier || ait == sample_to_allele.end() ? mr.allele_forms.end() : mr.allele_forms.find(ait->second);
                        if (fit != mr.allele_forms.end())
                            for (std::size_t k = 0; k < fit->second.size(); ++k) fs += (k ? "," : "") + std::to_string(fit->second[k] + 1);
                        row << ':' << (fs.empty() ? std::string(".") : fs);
                    }
                    if (e.cn_method == CnMethod::Unit && e.ru_len > 0) {
                        const auto iit = mr.sample_cninv.find(s);
                        row << ':' << (iit == mr.sample_cninv.end() ? std::string(".")
                                                                     : std::to_string(iit->second));
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
            rec.unaligned_missing = n_unaligned;
            rec.covered_missing = n_unknown - n_unaligned;
            rec.cn_held = std::move(held_fix);
            if (!held_ptrs.empty()) held_rows.emplace_back(out_records.size(), std::move(held_ptrs));
            id_of[&mr] = id;
            out_records.push_back(std::move(rec));

            if (options.write_variant_nodes) {
                std::string nodes;
                for (std::size_t k = 0; k < var_nodes.size(); ++k) { if (k) nodes += ','; nodes += var_nodes[k]; }
                variant_nodes_rows.push_back(
                    id + '\t' + std::to_string(bubble.id) + '\t' + svt + '\t' + nodes);

                // carrier_events.tsv: for each carrier, the event it was called for -- its own, not the
                // record's representative -- as the reference interval it replaces and the bases that
                // replace it, forward strand. Carriers with the same event share a row. benchmark
                // applies these to the reference and compares the result with each haplotype.
                std::map<std::string, std::vector<std::string>> ce_rows;   // row body -> carriers
                const auto fwd = [&](const std::string& w) { return ref_rev ? reverse_complement(w) : w; };
                // event_bp: the event as call measured it -- its net change, an inversion's length -- which
                // a node event's stretch can exceed: the stretch also holds the copies between its anchors.
                const auto body = [](const char* kind, std::size_t rs, std::size_t re, const std::string& cn,
                                     const std::string& alt, const std::string& walk_cn = ".",
                                     const std::string& event_bp = ".") {
                    const std::string eb = event_bp != "." ? event_bp
                        : std::to_string(static_cast<long long>(alt.size()) -
                                         (re + 1 >= rs ? static_cast<long long>(re + 1 - rs) : 0LL));
                    return std::string(kind) + '\t' + std::to_string(rs) + '\t' + std::to_string(re) + '\t' +
                           cn + '\t' + walk_cn + '\t' + eb + '\t' + std::to_string(alt.size()) + '\t' +
                           (alt.empty() ? "." : alt);
                };
                const std::string unresolved = "unresolved\t0\t0\t.\t.\t.\t0\t.";
                // An insertion after reference walk token k, in genomic terms: after its last base on a
                // forward bubble, before its first on a reverse one (ref_start = the base it precedes).
                const auto insertion_at = [&](std::size_t k, std::size_t& rs, std::size_t& re) {
                    const auto sp = tok_span(k, k);
                    if (sp.first == 0) return false;
                    if (!ref_rev) { rs = sp.second + 1; re = sp.second; }
                    else { rs = sp.first; re = sp.first - 1; }
                    return true;
                };
                // A node event is the stretch between the pairs its two walks were matched at on either side:
                // the reference steps strictly between, replaced by the allele's. A copy-number unit the
                // comparison left out lies between them on both walks, so each side keeps its own copies --
                // spanning it on the reference alone would write a substitution beside an array as the array's loss.
                const auto event_body = [&](const Event& ev, std::size_t ai) -> std::string {
                    const char* kind = ev.type == EvType::Ins ? "ins" : ev.type == EvType::Del ? "del"
                                     : ev.type == EvType::Sub ? "sub" : ev.type == EvType::Inv ? "inv" : nullptr;
                    if (ev.seq_event)   // a node replacement made a plain insertion or deletion keeps its node kind
                        return body(ev.node_found && kind != nullptr ? kind : "seq", ev.seq_before + 1,
                                    ev.seq_before + ev.del_bp, ".", ev.seq_ins);
                    if (kind == nullptr || ai >= alleles.size()) return unresolved;
                    if (ev.split_piece) {
                        const std::string ev_bp = std::to_string(
                            ev.type == EvType::Ins ? static_cast<long long>(ev.size_bp)
                            : ev.type == EvType::Sub ? static_cast<long long>(ev.ins_bp) - static_cast<long long>(ev.del_bp)
                                                     : -static_cast<long long>(ev.size_bp));
                        std::size_t rs = 0, re = 0;
                        if (ev.type == EvType::Ins) {
                            if (ev.ref_tok_anchor == SIZE_MAX || !insertion_at(ev.ref_tok_anchor, rs, re)) return unresolved;
                            return body(kind, rs, re, ".", fwd(ev.seq), ".", ev_bp);
                        }
                        const auto sp = tok_span(std::min(ev.ref_tok_first, ev.ref_tok_last), std::max(ev.ref_tok_first, ev.ref_tok_last));
                        return sp.first == 0 ? unresolved
                                             : body(kind, sp.first, sp.second, ".", ev.type == EvType::Sub ? fwd(ev.alt_seq) : "", ".", ev_bp);
                    }
                    std::size_t lo_ref, hi_ref;   // the event's reference steps; an insertion has none
                    if (ev.type == EvType::Ins) {
                        if (ev.ref_tok_anchor == SIZE_MAX) return unresolved;
                        lo_ref = ev.ref_tok_anchor + 1; hi_ref = ev.ref_tok_anchor;
                    } else {
                        if (ev.ref_tok_first == SIZE_MAX || ev.ref_tok_last == SIZE_MAX) return unresolved;
                        lo_ref = std::min(ev.ref_tok_first, ev.ref_tok_last);
                        hi_ref = std::max(ev.ref_tok_first, ev.ref_tok_last);
                    }
                    const auto& mp = allele_matched[ai];   // (allele step, reference step), increasing
                    const std::pair<std::size_t, std::size_t>* before = nullptr;
                    const std::pair<std::size_t, std::size_t>* after = nullptr;
                    for (const auto& pr : mp) {
                        if (pr.second < lo_ref) before = &pr;
                        else if (pr.second > hi_ref && after == nullptr) after = &pr;
                    }
                    if (before == nullptr || after == nullptr || after->first <= before->first ||
                        after->second <= before->second)
                        return unresolved;
                    const std::string ev_bp = std::to_string(
                        ev.type == EvType::Inv ? static_cast<long long>(ev.size_bp)
                        : ev.type == EvType::Ins ? static_cast<long long>(ev.size_bp)
                        : ev.type == EvType::Del ? -static_cast<long long>(ev.size_bp)
                        : static_cast<long long>(ev.ins_bp) - static_cast<long long>(ev.del_bp));
                    const std::vector<PathStep>& hw = alleles[ai].steps;
                    const std::string alt = fwd(spell_path_steps_sequence(
                        graph, std::vector<PathStep>(hw.begin() + static_cast<std::ptrdiff_t>(before->first) + 1,
                                                     hw.begin() + static_cast<std::ptrdiff_t>(after->first))));
                    if (after->second == before->second + 1) {   // nothing of the reference between: insertion
                        std::size_t rs = 0, re = 0;
                        if (!insertion_at(before->second, rs, re)) return unresolved;
                        return body(kind, rs, re, ".", alt, ".", ev_bp);
                    }
                    const auto sp = tok_span(before->second + 1, after->second - 1);
                    return sp.first == 0 ? unresolved : body(kind, sp.first, sp.second, ".", alt, ".", ev_bp);
                };
                // A copy-number carrier was counted over the record's own nodes: the reference's run of
                // them is replaced by the carrier's run, which holds the copies it was counted with.
                // The copies are counted over the record's event (unit) nodes.
                std::unordered_set<std::string> cn_nodes, cn_count_nodes;
                if (e.type == EvType::Dup) {
                    for (const std::string& nd : var_nodes)
                        if (nd != bubble.source && nd != bubble.sink) cn_nodes.insert(nd);
                    for (const std::string& nd : ev_nodes)
                        if (nd != bubble.source && nd != bubble.sink) cn_count_nodes.insert(nd);
                }
                const auto run_of = [&](const std::vector<PathStep>& w, std::size_t& a, std::size_t& b) {
                    a = SIZE_MAX; b = SIZE_MAX;
                    for (std::size_t k = 1; k + 1 < w.size(); ++k)
                        if (cn_nodes.count(w[k].node_id)) { if (a == SIZE_MAX) a = k; b = k; }
                    return a != SIZE_MAX;
                };
                std::unordered_map<std::size_t, std::pair<bool, Stretch>> stretch_of;   // per allele, found once
                ref_junctions.clear();   // keyed by this record's copy-node set, which lives in this iteration
                const auto cn_body = [&](std::size_t ai, const std::string& cn) -> std::string {
                    if (cn_nodes.empty() || ai >= alleles.size()) return unresolved;
                    const std::vector<PathStep>& hw = alleles[ai].steps;
                    std::size_t ra = 0, rb = 0, ha = 0, hb = 0;
                    const bool in_ref = run_of(ref_steps, ra, rb), in_hap = run_of(hw, ha, hb);
                    // The copies the carrier's walk holds, by the record's own rule. MODULE: the length-weighted
                    // median of its visits over the anchor nodes. UNIT: its steps through the unit, its twins
                    // and the unit graph's anchors included -- the key node alone would read a twin as no copy.
                    std::string wcn = ".";
                    if (!mr.anchor_nodes.empty()) {
                        std::unordered_map<std::string, std::size_t> visits;
                        for (std::size_t k = 1; k + 1 < hw.size(); ++k)
                            if (mr.anchor_nodes.count(hw[k].node_id)) ++visits[hw[k].node_id];
                        std::vector<std::pair<std::size_t, std::size_t>> cw;   // (visits, length)
                        std::size_t total = 0;
                        for (const std::string& id : mr.anchor_nodes) {
                            const auto it = visits.find(id);
                            const std::size_t len = node_len(graph, id);
                            cw.emplace_back(it == visits.end() ? 0 : it->second, len);
                            total += len;
                        }
                        std::sort(cw.begin(), cw.end());
                        std::size_t acc = 0, med = cw.empty() ? 0 : cw.back().first;
                        for (const auto& [c, w] : cw) { acc += w; if (2 * acc >= total) { med = c; break; } }
                        wcn = std::to_string(med);
                    } else if (!cn_count_nodes.empty()) {
                        std::unordered_map<std::string, std::size_t> copies;
                        std::size_t walk_cn = 0;
                        for (std::size_t k = 1; k + 1 < hw.size(); ++k) {
                            const std::string* u = rep_unit_of(hw[k].node_id);
                            const std::string& key = u != nullptr ? *u : hw[k].node_id;
                            if (cn_count_nodes.count(key)) walk_cn = std::max(walk_cn, ++copies[key]);
                        }
                        wcn = std::to_string(walk_cn);
                    }
                    // An allele whose events come from aligning its bases (its copies left out of both walks) has
                    // every change beside its copies as an event of its own. Its stretch is cut where that
                    // alignment cuts the two walks: from the first copy on either walk to the last, each place on
                    // one walk taken to the other by the alignment, and a change of the alignment never cut in two.
                    // So no base of the haplotype is both in the stretch and in an event of its own outside it,
                    // and a base between copies on one walk is in the stretch on the other too.
                    if (allele_seq[ai] && in_ref) {
                        auto it = stretch_of.find(ai);
                        if (it == stretch_of.end()) {
                            Stretch st;
                            const bool ok = aligned_stretch(ai, cn_nodes, st);
                            it = stretch_of.emplace(ai, std::make_pair(ok, std::move(st))).first;
                        }
                        if (it->second.first)
                            return body("cn", it->second.second.g_first, it->second.second.g_last, cn, it->second.second.alt, wcn);
                    }
                    // Otherwise the stretch between the steps the carrier's walk and the reference were matched at
                    // on either side of the copies, as a node event is written: the carrier's own sequence beside
                    // or between its copies (an insertion, a deletion of flank) is in it, so whatever order the
                    // walk has them in is what benchmark applies. Its own events inside the stretch are part of
                    // it; benchmark checks the stretch's length against the CN plus those events.
                    const auto& mp = allele_matched[ai];
                    const std::pair<std::size_t, std::size_t>* before = nullptr;
                    const std::pair<std::size_t, std::size_t>* after = nullptr;
                    for (const auto& pr : mp) {
                        const bool lo = (!in_hap || pr.first < ha) && (!in_ref || pr.second < ra);
                        const bool hi = (!in_hap || pr.first > hb) && (!in_ref || pr.second > rb);
                        if (lo) before = &pr;
                        else if (hi && after == nullptr) after = &pr;
                    }
                    if (before == nullptr || after == nullptr || after->first <= before->first ||
                        after->second <= before->second)
                        return unresolved;
                    const std::string alt = fwd(spell_path_steps_sequence(
                        graph, std::vector<PathStep>(hw.begin() + static_cast<std::ptrdiff_t>(before->first) + 1,
                                                     hw.begin() + static_cast<std::ptrdiff_t>(after->first))));
                    if (after->second == before->second + 1) {   // no reference step between: an insertion
                        std::size_t rs = 0, re = 0;
                        if (!insertion_at(before->second, rs, re)) return unresolved;
                        return body("cn", rs, re, cn, alt, wcn);
                    }
                    const auto sp = tok_span(before->second + 1, after->second - 1);
                    return sp.first == 0 ? unresolved : body("cn", sp.first, sp.second, cn, alt, wcn);
                };
                for (const std::string& m : mr.carriers) {
                    const auto ait = sample_to_allele.find(m);
                    std::string b = unresolved;
                    if (ait != sample_to_allele.end()) {
                        if (e.type == EvType::Dup) {
                            const auto cit = mr.sample_cn.find(m);
                            b = cn_body(ait->second, cit == mr.sample_cn.end() ? std::string(".")
                                                                                  : std::to_string(cit->second));
                        } else if (const auto mit = mr.member_event.find(ait->second); mit != mr.member_event.end()) {
                            for (const Event& me : mit->second) ce_rows[event_body(me, ait->second)].push_back(m);
                            continue;
                        }
                    }
                    ce_rows[b].push_back(m);
                }
                for (auto& [b, haps] : ce_rows) {
                    std::sort(haps.begin(), haps.end());
                    std::string row_ce = id + '\t' + std::to_string(bubble.id) + '\t' + b + '\t';
                    for (std::size_t k = 0; k < haps.size(); ++k) { if (k) row_ce += ','; row_ce += haps[k]; }
                    carrier_events_rows.push_back(std::move(row_ce));
                }
            }

            // Record a per-gene-CN target for this copy-number record (resolved from k-mers after the
            // bubbles are called), with the per-haplotype CN so unseparable gene groups can report the
            // total, and the nodes of what it counts: a MODULE's copies, a UNIT's repeat unit (its key,
            // twins and unit-graph nodes). Which genes lie in those copies is decided at the end: a tandem
            // gene copy folded into a unit counts the gene, a unit that is only a piece of a gene does not.
            if (e.type == EvType::Dup && !bubble_gene_idx.empty()) {
                std::unordered_set<std::string> counted = mr.module_nodes;
                if (counted.empty()) {
                    std::unordered_set<std::string> keys;
                    for (const std::string& nd : ev_nodes) if (nd != bubble.source && nd != bubble.sink) keys.insert(nd);
                    for (const std::string& nd : bubble.inside) {
                        const auto m = ug_member_of.find(nd);
                        const std::string* u = m != ug_member_of.end() ? &m->second : rep_unit_of(nd);
                        if (keys.count(nd) || (u != nullptr && keys.count(*u))) counted.insert(nd);
                    }
                }
                DupGeneTarget t{bubble.id, id, bubble_gene_idx, mr.sample_cn, std::move(counted), {},
                                !mr.anchor_nodes.empty()};
                dup_targets.push_back(std::move(t));
            }

            ++summary.records_written;
            if (rp_sub) ++summary.sub;
            else if (rp) { if (svlen < 0) ++summary.del; else ++summary.ins; }
            else if (e.type == EvType::Sub && !sub_sized) ++summary.sub;
            else if (e.type == EvType::Del || (sub_sized && sub_net < 0)) ++summary.del;
            else if (e.type == EvType::Ins || sub_sized) ++summary.ins;
            else if (e.type == EvType::Inv) ++summary.inv;
            else ++summary.dup;
        }
        // The holders' ids, now that every record of the site is written; a sample none of whose holders was
        // written takes its measured CN at once.
        for (auto& [k, ptrs] : held_rows) {
            OutRecord& rec = out_records[k];
            for (std::size_t f = 0; f < ptrs.size(); ++f)
                for (const MergedRecord* r : ptrs[f])
                    if (const auto it = id_of.find(r); it != id_of.end()) rec.cn_held[f].holders.push_back(it->second);
            apply_measured_cn(rec, [](const std::string&) { return false; });
        }
    };
    // Ticked as bubbles are processed, not in the merge below, so the bar follows the run.
    for (std::size_t bubble_idx = 0; bubble_idx < bubbles.size(); ++bubble_idx) {
        process_bubble(bubble_idx);
        call_progress.tick();
    }

    // Merge per-bubble outputs in bubble order: deterministic regardless of thread count.
    for (std::size_t bi = 0; bi < bouts.size(); ++bi) {
        BubbleOut& bo = bouts[bi];
        summary.bubbles_seen += bo.sum.bubbles_seen;
        summary.bubbles_with_reference += bo.sum.bubbles_with_reference;
        summary.bubbles_with_calls += bo.sum.bubbles_with_calls;
        summary.records_written += bo.sum.records_written;
        summary.del += bo.sum.del;
        summary.ins += bo.sum.ins;
        summary.inv += bo.sum.inv;
        summary.dup += bo.sum.dup;
        summary.sub += bo.sum.sub;
        summary.tangle_ids.insert(summary.tangle_ids.end(), bo.sum.tangle_ids.begin(), bo.sum.tangle_ids.end());
        summary.oversized_ids.insert(summary.oversized_ids.end(), bo.sum.oversized_ids.begin(), bo.sum.oversized_ids.end());
        summary.skipped_large_segments += bo.sum.skipped_large_segments;
        summary.copy_carriers += bo.sum.copy_carriers;
        summary.copy_own_event += bo.sum.copy_own_event;
        summary.copy_no_form += bo.sum.copy_no_form;
        summary.copy_pass_off += bo.sum.copy_pass_off;
        // Every per-bubble counter has to be carried over here, or the run summary never reports it.
        summary.multi_crossing_bubbles += bo.sum.multi_crossing_bubbles;
        summary.multi_crossing_paths += bo.sum.multi_crossing_paths;
        summary.multi_crossing_unscored_bubbles += bo.sum.multi_crossing_unscored_bubbles;
        summary.multi_crossing_unscored_paths += bo.sum.multi_crossing_unscored_paths;
        summary.multi_crossing_unscored_ids.insert(summary.multi_crossing_unscored_ids.end(),
                                                   bo.sum.multi_crossing_unscored_ids.begin(),
                                                   bo.sum.multi_crossing_unscored_ids.end());
        summary.gap_missing += bo.sum.gap_missing;
        summary.ref_gap_bubbles += bo.sum.ref_gap_bubbles;
        summary.ref_gap_ids.insert(summary.ref_gap_ids.end(), bo.sum.ref_gap_ids.begin(), bo.sum.ref_gap_ids.end());
        summary.unaligned_missing += bo.sum.unaligned_missing;
        summary.covered_missing += bo.sum.covered_missing;
        summary.unaligned_bubbles.insert(summary.unaligned_bubbles.end(), bo.sum.unaligned_bubbles.begin(),
                                         bo.sum.unaligned_bubbles.end());
        for (OutRecord& r : bo.records) out_records.push_back(std::move(r));
        for (std::string& s : bo.variant_nodes) variant_nodes_rows.push_back(std::move(s));
        for (std::string& s : bo.carrier_events) carrier_events_rows.push_back(std::move(s));
        for (DupGeneTarget& t : bo.dup_targets) dup_targets.push_back(std::move(t));
    }

    // ---- Inversion sites. A large inversion whose breakpoints lie in inverted repeats the graph collapses
    // has no snarl of its own: every path, the reference too, walks the repeat nodes in both orientations,
    // so the snarls at the breakpoints have no interior a path crosses. bubble makes it a site from the
    // walks (Bubble::inversion_parent), enclosing the sites
    // inside it, and it is typed the same way: each haplotype, turned to the reference's orientation, as
    // runs of the nodes the reference holds once, walked with or against it. A reversed run between two
    // forward runs, its reference positions falling, is an inversion; it is written at the inversion site
    // whose flanks hold it.
    bool any_inversion_site = false;
    for (const Bubble& bb : bubbles) any_inversion_site = any_inversion_site || bb.inversion_parent;
    if (any_inversion_site) {
        const PathInversionScan scan = find_path_inversions(
            graph, *ref_path, options.min_sv_bp > 0 ? options.min_sv_bp : kPathInversionMinBp);
        const std::size_t g0 = ref_meta.region_start_1based;   // genomic, 1-based
        std::vector<std::pair<std::size_t, std::size_t>> span(graph.paths.size(), {SIZE_MAX, 0});
        for (std::size_t pi = 0; pi < graph.paths.size(); ++pi)
            if (scan.span[pi].first != SIZE_MAX) span[pi] = {g0 + scan.span[pi].first, g0 + scan.span[pi].second};
        // A reversed run inside one site's interior is that site's to type -- an inversion, or something
        // folded back on itself -- and only one no ordinary site's interior holds is written at an inversion
        // site (the rule bubble makes inversion sites by). site_span is an ordinary site's interior.
        std::vector<std::pair<std::size_t, std::size_t>> site_span;
        std::vector<std::pair<std::size_t, std::pair<std::size_t, std::size_t>>> parent_span;   // bubble index, span
        {
            std::unordered_map<std::string, std::size_t> first_at;
            for (std::size_t k = 0; k < ref_path->steps.size(); ++k)
                first_at.emplace(ref_path->steps[k].node_id, ref_meta.region_start_1based + ref_prefix[k]);
            for (std::size_t bi = 0; bi < bubbles.size(); ++bi) {
                const Bubble& bb = bubbles[bi];
                const auto a = first_at.find(bb.source), b = first_at.find(bb.sink);
                if (a == first_at.end() || b == first_at.end()) continue;
                const std::pair<std::size_t, std::size_t> sp{
                    std::min(a->second, b->second),
                    std::max(a->second + node_len(graph, a->first), b->second + node_len(graph, b->first))};
                if (bb.inversion_parent) { parent_span.emplace_back(bi, sp); continue; }
                const bool src_first = a->second <= b->second;
                site_span.emplace_back(src_first ? a->second + node_len(graph, a->first) : b->second + node_len(graph, b->first),
                                       src_first ? b->second : a->second);
            }
        }
        // A parent's own nodes: its interior less its children's -- the stretch the inversion turns.
        const auto own_nodes = [&](const Bubble& parent) {
            std::unordered_set<std::string> child;
            for (const Bubble& bb : bubbles)
                if (bb.parent_id == parent.id) child.insert(bb.inside.begin(), bb.inside.end());
            std::string out;
            for (const std::string& n : parent.inside)
                if (!child.count(n)) { if (!out.empty()) out += ','; out += n; }
            return out;
        };
        // An inversion a site already typed is not written again.
        std::vector<std::pair<std::size_t, std::size_t>> site_inv;
        for (const OutRecord& r : out_records) {
            if (r.line.find("\tSVTYPE=INV;") == std::string::npos && r.line.find(";SVTYPE=INV;") == std::string::npos)
                continue;
            if (options.min_maf > 0.0 && af_below(written(r.line))) continue;   // not written
            site_inv.emplace_back(r.pos, r.end);
        }
        std::size_t written = 0, no_site = 0;
        std::unordered_map<std::size_t, int> per_parent;
        for (const PathInversion& inv_site : scan.inversions) {
            // The representative breakpoints: the flanks most carriers share.
            const std::size_t left = g0 + inv_site.left_end, a = g0 + inv_site.inv_lo, b = g0 + inv_site.inv_hi,
                              right = g0 + inv_site.right_start;
            const std::size_t pos = left, end = right - 1;   // POS: last forward base; END: last before the right flank
            bool dup = false;
            for (const auto& [sp, se] : site_span) if (sp <= a && b < se) dup = true;
            for (const auto& [sp, se] : site_inv) {
                const std::size_t lo = std::max(sp, pos), hi = std::min(se, end);
                if (hi > lo && 2 * (hi - lo) >= end - pos) dup = true;
            }
            if (dup) continue;
            // The inversion site whose flanks hold the reversed run.
            const Bubble* parent = nullptr;
            for (const auto& [bi, sp] : parent_span)
                if (sp.first <= a && b < sp.second) { parent = &bubbles[bi]; break; }
            if (parent == nullptr) { ++no_site; continue; }
            if (!bubble_filter.empty() && !bubble_filter.count(parent->id)) continue;
            const std::unordered_set<std::size_t> carriers(inv_site.carriers.begin(), inv_site.carriers.end());
            std::string gts;
            std::size_t an = 0, ac = 0;
            std::vector<std::string> carrier_names;
            for (std::size_t pi = 0; pi < graph.paths.size(); ++pi) {
                std::string gt;
                if (carriers.count(pi)) { gt = "1"; ++an; ++ac; carrier_names.push_back(graph.paths[pi].name); }
                else if (span[pi].first <= left && span[pi].second >= right) { gt = "0"; ++an; }
                else gt = ".";
                gts += '\t' + gt + ":.";
            }
            // the same carrier filters as every other record (--min-alt-af on the written AF, below)
            if (ac == 0 || ac < options.min_haplotypes) continue;
            const long long cipos = static_cast<long long>(a - left - 1), ciend = static_cast<long long>(right - b - 1);
            const std::string base_id = "bubble" + std::to_string(parent->id) + "_INV_" + parent->source;
            const int seen_n = ++per_parent[parent->id];
            const std::string id = seen_n == 1 ? base_id : base_id + "_" + std::to_string(seen_n);
            char afb[32];
            std::snprintf(afb, sizeof(afb), "%.4g", static_cast<double>(ac) / static_cast<double>(an));
            std::ostringstream info;
            info << "END=" << end << ";SVTYPE=INV;SVLEN=" << (end - pos) << ";PATH_INVERSION;BUBBLE_ID=" << parent->id;
            if (cipos > 0 || ciend > 0) info << ";IMPRECISE;CIPOS=0," << cipos << ";CIEND=-" << ciend << ",0";
            info << ";AN=" << an << ";AC=" << ac << ";AF=" << afb;
            if (!genes.empty()) {
                const std::vector<std::string> gnames = genes_in_interval(pos + 1, end);
                if (!gnames.empty()) {
                    info << ";GENES=";
                    for (std::size_t k = 0; k < gnames.size(); ++k) { if (k) info << ','; info << vcf_escape(gnames[k]); }
                }
            }
            const std::size_t off = pos - ref_meta.region_start_1based;
            const std::string ref_base = off < ref_genome.size() ? std::string(1, ref_genome[off]) : std::string("N");
            OutRecord rec;
            rec.pos = pos; rec.end = end; rec.bubble_id = parent->id; rec.id = id;
            rec.line = ref_meta.chrom + '\t' + std::to_string(pos) + '\t' + id + '\t' + ref_base + "\t<INV>\t.\tPASS\t" +
                       info.str() + "\tGT:CN" + gts + '\n';
            out_records.push_back(std::move(rec));
            if (options.write_variant_nodes) {
                // benchmark applies the inverted reference stretch to each carrier; the repeats at the
                // breakpoints come out as small differences
                std::sort(carrier_names.begin(), carrier_names.end());
                const std::size_t g0 = pos + 1, g1 = end;
                const std::string inv = reverse_complement(ref_genome.substr(g0 - ref_meta.region_start_1based, g1 - g0 + 1));
                variant_nodes_rows.push_back(id + '\t' + std::to_string(parent->id) + "\tINV\t" + own_nodes(*parent));
                std::string row = id + '\t' + std::to_string(parent->id) + "\tinv\t" + std::to_string(g0) + '\t' + std::to_string(g1) + "\t.\t.\t" +
                                  std::to_string(g1 - g0 + 1) + '\t' + std::to_string(inv.size()) + '\t' + inv + '\t';
                for (std::size_t k = 0; k < carrier_names.size(); ++k) { if (k) row += ','; row += carrier_names[k]; }
                carrier_events_rows.push_back(std::move(row));
            }
            ++written;
            ++summary.records_written;
            ++summary.inv;
        }
        if (written > 0 && !options.quiet)
            std::cerr << cli::log_prefix("call") << "inversions typed at inversion sites: " << written << "\n";
        if (no_site > 0)
            std::cerr << cli::log_prefix("call") << "warning: " << no_site << " inversion(s) read from the walks lie in no inversion site of "
                      << options.bubbles_csv_in << " and are not written (made at a different --min-variant-bp?)\n";
    }

    // ---- Coordinate-sort all records and write the (indexable) region + per-bubble VCFs.
    std::stable_sort(out_records.begin(), out_records.end(),
                     [](const OutRecord& a, const OutRecord& b) {
                         if (a.pos != b.pos) return a.pos < b.pos;
                         if (a.end != b.end) return a.end < b.end;
                         return a.id < b.id;
                     });
    // A record is not written when no sample column carries it (--diploid), or when its written AF is under
    // --min-alt-af. Found once, before the outputs are listed, so the per-site files, the
    // node and carrier tables and the run's counts all leave it out.
    if (options.diploid || options.min_maf > 0.0) {
        std::unordered_set<std::string> dropped;
        std::size_t no_column = 0;
        // A module sample read as REF_CN because a dropped record held its change takes its measured CN. That
        // changes the line it is written with, and with --diploid can make a slot's contigs agree, which adds
        // to AN: the records are tested again until none is dropped (each restoration happens once).
        for (bool again = true; again;) {
            out_records.erase(std::remove_if(out_records.begin(), out_records.end(), [&](const OutRecord& r) {
                const std::string w = written(r.line);
                if (w.empty()) ++no_column;
                else if (options.min_maf <= 0.0 || !af_below(w)) return false;
                dropped.insert(r.id);
                summary.unaligned_missing -= r.unaligned_missing;
                summary.covered_missing -= r.covered_missing;
                return true;
            }), out_records.end());
            again = false;
            for (OutRecord& r : out_records) {
                const std::string before = r.line;
                apply_measured_cn(r, [&](const std::string& id) { return dropped.count(id) != 0; });
                if (r.line != before) again = true;
            }
        }
        summary.diploid_records_dropped = no_column;
        // a dropped copy-number record's genes keep their per-gene rows, under no record
        for (DupGeneTarget& t : dup_targets) if (dropped.count(t.variant_id)) t.variant_id = ".";
        if (!dropped.empty()) {
            const auto keep = [&](std::vector<std::string>& rows) {
                rows.erase(std::remove_if(rows.begin(), rows.end(), [&](const std::string& row) {
                    return dropped.count(row.substr(0, row.find('\t'))) != 0;
                }), rows.end());
            };
            keep(variant_nodes_rows);
            keep(carrier_events_rows);
            summary.records_written = out_records.size();
            summary.del = summary.ins = summary.inv = summary.dup = summary.sub = 0;
            // 'with records' counts ordinary sites, as 'sites' does: an inversion record carries its
            // inversion site's id
            std::unordered_set<std::size_t> inversion_sites;
            for (const Bubble& b : bubbles) if (b.inversion_parent) inversion_sites.insert(b.id);
            std::set<std::size_t> sites;
            for (const OutRecord& r : out_records) {
                const std::size_t a = r.line.find("SVTYPE=");
                const std::string type = a == std::string::npos ? "" : r.line.substr(a + 7, 3);
                if (type == "DEL") ++summary.del;
                else if (type == "INS") ++summary.ins;
                else if (type == "INV") ++summary.inv;
                else if (type == "CNV") ++summary.dup;
                else if (type == "SUB") ++summary.sub;
                if (r.bubble_id != 0 && !inversion_sites.count(r.bubble_id)) sites.insert(r.bubble_id);
            }
            summary.bubbles_with_calls = sites.size();
        }
    }

    // Every destination is resolved and checked before anything is opened, then written to a staged
    // sibling and renamed in only once the whole run has succeeded, so a failure part-way leaves no
    // complete-looking output and a run cannot overwrite its own input.
    cli::StagedOutputs staged("call");
    std::vector<std::string> finals;
    finals.push_back(options.out_prefix + ".region.vcf");
    if (options.write_variant_nodes) {
        finals.push_back(options.out_prefix + ".variant_nodes.tsv");
        finals.push_back(options.out_prefix + ".carrier_events.tsv");
    }
    std::set<std::size_t> bubble_ids_out;
    if (options.write_per_bubble_vcf) {
        for (const OutRecord& rec : out_records) if (rec.bubble_id != 0) bubble_ids_out.insert(rec.bubble_id);
        for (const std::size_t bid : bubble_ids_out)
            finals.push_back(options.out_prefix + ".bubble_" + std::to_string(bid) + ".vcf");
    }
    // The GTF sidecars belong to the same transaction as the VCFs, so a sidecar failure installs
    // nothing, and the collision check below covers them too (an input GTF named
    // <prefix>.node_genes.tsv is not overwritten by the run reading it).
    const bool gtf_active = !options.gtf_path.empty() && !genes.empty();
    if (gtf_active) {
        finals.push_back(options.out_prefix + ".node_genes.tsv");
        finals.push_back(options.out_prefix + ".dup_gene_cn.tsv");
    }
    // The files this run reads: no output may overwrite one, and the stale-file sweeps below never remove one.
    const std::vector<std::string> inputs = {options.gfa_path, options.bubbles_csv_in, options.gtf_path,
                                             options.unit_graph_tsv};
    const auto is_input = [&](const std::string& f) {
        for (const std::string& in : inputs) {
            if (in.empty()) continue;
            std::error_code ec;
            if (f == in || std::filesystem::equivalent(f, in, ec)) return true;
        }
        return false;
    };
    {
        std::unordered_set<std::string> seen_out;
        for (const std::string& f : finals) {
            if (!seen_out.insert(f).second)
                throw std::runtime_error("call: two outputs would be written to the same file: " + f);
            if (is_input(f))
                throw std::runtime_error("call: output '" + f + "' is also an input; refusing to "
                                         "overwrite the data being read");
        }
    }
    // Per-bubble VCFs an earlier run left at this prefix and this run will not rewrite. Collected now
    // but removed only after the transaction commits, so a later failure does not destroy part of a
    // previous result. Not gated on --write-per-bubble-vcf either: a rerun with --no-per-bubble-vcf
    // would otherwise leave the old files looking like current output.
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
            std::size_t id = 0;
            const auto parsed = std::from_chars(mid.data(), mid.data() + mid.size(), id);
            if (parsed.ec != std::errc()) continue;   // too large for a site id: no file of call's
            if (options.write_per_bubble_vcf && bubble_ids_out.count(id)) continue;
            obsolete_bubble_vcfs.push_back(de.path());
        }
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

        // carrier_events.tsv: per record, each carrier's own event, applied by `benchmark`. Coordinates
        // are 1-based on the reference; an insertion has ref_end = ref_start - 1 and goes in before
        // ref_start. `alt` is forward strand; kind says which route the event came from, and
        // `unresolved` rows are carriers whose event could not be placed on the reference. On a
        // copy-number record, cn is the CN reported and walk_cn the copies the carrier's walk holds,
        // counted by the record's own rule (see cn_body). event_bp is the event as call measured it (net
        // change, an inversion's length), which a node event's stretch can exceed: the stretch runs between
        // the anchors where its two walks meet again and holds whatever copies lie between them.
        const std::string ce_path = options.out_prefix + ".carrier_events.tsv";
        std::ofstream ce_out(staged.stage(ce_path));
        if (!ce_out) throw std::runtime_error("Failed to write carrier events TSV: " + ce_path);
        ce_out << "variant_id\tbubble_id\tkind\tref_start\tref_end\tcn\twalk_cn\tevent_bp\talt_bp\talt\thaplotypes\n";
        for (const std::string& r : carrier_events_rows) ce_out << r << '\n';
        ce_out.close();
        if (!ce_out) throw std::runtime_error("Failed to finalize carrier events TSV: " + ce_path);
    }

    // GTF sidecars: node->genes map and per-gene DUP CN. CN comes from private-k-mer dosage -- each
    // gene's discriminative reference (merged CDS, else gene span) yields k-mers unique against its
    // paralogs, and a haplotype's per-copy count of those is its CN. Separates collapsed paralogs that
    // graph multiplicity cannot, without per-haplotype alignment.
    if (!genes.empty()) {
        write_node_genes_tsv(staged.stage(options.out_prefix + ".node_genes.tsv"), node_genes, genes);

        // Which of a target's genes lie in the module's copies: most of the gene's coding sequence (its
        // span if it has none) lies in the reference stretch the module's nodes span, its copies and what
        // lies between them. The module's total counts those copies, so only such a gene can be given it
        // -- copied whole (a single gene, or one with no private k-mers) or split between a near-identical
        // pair. A gene beside the copies, or one a short module repeats a piece of (a repeat inside a
        // gene), keeps only its own k-mer evidence. Measured on the reference stretch rather than the
        // nodes, since a gene whose exons differ between copies sits on nodes no path revisits.
        if (!dup_targets.empty()) {
            std::unordered_map<std::string, std::vector<std::pair<std::size_t, std::size_t>>> ref_iv;
            {
                std::size_t at = ref_meta.region_start_1based;
                for (const PathStep& st : ref_path->steps) {
                    const std::size_t len = graph.nodes.at(st.node_id).sequence.size();
                    if (len > 0) ref_iv[st.node_id].emplace_back(at, at + len - 1);
                    at += len;
                }
            }
            for (DupGeneTarget& t : dup_targets) {
                // A module's copies hold nodes no path revisits, so its stretch runs from its first node to its
                // last. A repeat unit's copies are its own nodes: a gene between two copies that are not
                // adjacent is not in them.
                std::vector<std::pair<std::size_t, std::size_t>> stretch;
                for (const std::string& id : t.module_nodes)
                    if (const auto it = ref_iv.find(id); it != ref_iv.end())
                        stretch.insert(stretch.end(), it->second.begin(), it->second.end());
                std::sort(stretch.begin(), stretch.end());
                std::vector<std::pair<std::size_t, std::size_t>> merged_iv;
                for (const auto& iv : stretch) {
                    if (!merged_iv.empty() && (t.module || iv.first <= merged_iv.back().second + 1))
                        merged_iv.back().second = std::max(merged_iv.back().second, iv.second);
                    else merged_iv.push_back(iv);
                }
                t.in_copies.assign(t.gene_idx.size(), 0);
                for (std::size_t k = 0; k < t.gene_idx.size(); ++k) {
                    const GeneFeature& g = genes[t.gene_idx[k]];
                    std::vector<std::pair<std::size_t, std::size_t>> cds = g.cds;
                    if (cds.empty()) cds.emplace_back(g.start_1based, g.end_1based);
                    std::size_t total = 0, on = 0;
                    for (const auto& [a, b] : cds) {
                        total += b - a + 1;
                        for (const auto& [lo, hi] : merged_iv) {
                            const std::size_t x = std::max(a, lo), y = std::min(b, hi);
                            if (x <= y) on += y - x + 1;
                        }
                    }
                    t.in_copies[k] = 2 * on > total ? 1 : 0;
                }
            }
            // A single gene beside the copies has no evidence of its own here: no row. And one row per
            // gene and haplotype: where a site has a MODULE record, a UNIT record of the same site (a
            // repeat inside or beside the module) does not write its genes again.
            std::unordered_set<std::size_t> module_bubbles;
            for (const DupGeneTarget& t : dup_targets) if (t.module) module_bubbles.insert(t.bubble_id);
            dup_targets.erase(std::remove_if(dup_targets.begin(), dup_targets.end(), [&](const DupGeneTarget& t) {
                return (t.gene_idx.size() == 1 && !t.in_copies[0]) || (!t.module && module_bubbles.count(t.bubble_id));
            }), dup_targets.end());
        }
        if (!dup_targets.empty()) {
            // Per-gene resolution (spelling every haplotype + building k-mer sets) is only meaningful when
            // a target folds >=2 genes; a single-gene DUP's copy number is just the module total.
            const bool any_multi = std::any_of(dup_targets.begin(), dup_targets.end(),
                [](const DupGeneTarget& t) { return t.gene_idx.size() >= 2; });

            // Module-total row (single gene, or a gene with no private k-mers). Evidence columns are "."
            // because no per-gene k-mer split was made.
            // A haplotype with no module total gets cn "." and is never marked reliable.
            auto emit_total = [&](const DupGeneTarget& t, const std::string& names, const std::string& ids,
                                  const char* reliable) {
                for (const PathRecord& p : graph.paths) {
                    const auto cit = t.sample_cn.find(p.name);
                    const bool known = cit != t.sample_cn.end();
                    const std::string total = known ? std::to_string(cit->second) : ".";
                    dup_gene_cn_rows.push_back(std::to_string(t.bubble_id) + '\t' + t.variant_id + '\t' +
                                               p.name + '\t' + names + '\t' + total + '\t' +
                                               (known ? reliable : "0") + "\t.\t.\t.\t.\t.\t" + ids);
                }
            };

            std::vector<std::string> hap_seq;                     // per graph.paths (spelled once)
            std::unordered_map<int, std::string> marker_cache;    // gene index -> discriminative marker
            std::vector<std::vector<std::string>> target_backgrounds;   // per target, per gene
            std::vector<std::vector<std::vector<std::string>>> target_exons;   // per target, per gene, per coding exon
            std::vector<std::vector<std::string>> target_spans;            // per target, per gene
            if (any_multi) {
                const std::string ref_full_seq = spell_path_steps_sequence(graph, ref_path->steps);
                const std::size_t region_start = ref_meta.region_start_1based;
                // The reference bases of [lo, hi] (1-based, inclusive), clipped to the region on both sides.
                auto ref_slice = [&](std::size_t lo, std::size_t hi) -> std::string {
                    const std::size_t region_end = region_start + ref_full_seq.size();   // exclusive
                    lo = std::max(lo, region_start);
                    hi = std::min(hi + 1, region_end);                                  // exclusive
                    if (lo >= hi) return std::string();
                    return ref_full_seq.substr(lo - region_start, hi - lo);
                };
                // Merged CDS, else the gene span (CDS-less genes). The CDS pieces are joined with an N,
                // which resets the k-mer window: joined end to end, every exon boundary would make k-1
                // k-mers that exist nowhere in the genome, private markers no haplotype can hit, and
                // dosage would read low by their share.
                auto marker_of = [&](int gi) -> std::string {
                    const GeneFeature& g = genes[gi];
                    if (g.cds.empty()) return ref_slice(g.start_1based, g.end_1based);
                    std::string out;
                    for (const auto& iv : g.cds) {
                        const std::string piece = ref_slice(iv.first, iv.second);
                        if (piece.empty()) continue;
                        if (!out.empty()) out += 'N';
                        out += piece;
                    }
                    return out;
                };
                for (const DupGeneTarget& t : dup_targets)
                    if (t.gene_idx.size() >= 2)
                        for (int gi : t.gene_idx)
                            if (!marker_cache.count(gi)) marker_cache[gi] = marker_of(gi);
                // Each gene's background: the full spans of its siblings in the target -- introns and
                // UTRs too, where the sibling screen above sees only their CDS -- minus the gene's own
                // span, which an annotated readthrough sibling can contain whole. Per target, since the
                // siblings differ between targets.
                auto background_of = [&](int gi, const std::vector<int>& siblings) -> std::string {
                    std::vector<char> keep(ref_full_seq.size(), 0);
                    auto mark = [&](const GeneFeature& g, char v) {
                        const std::size_t lo = std::max(g.start_1based, region_start);
                        const std::size_t hi = std::min(g.end_1based + 1, region_start + keep.size());
                        for (std::size_t x = lo; x < hi; ++x) keep[x - region_start] = v;
                    };
                    for (int o : siblings) if (o != gi) mark(genes[o], 1);
                    mark(genes[gi], 0);
                    std::string bg = ref_full_seq;
                    for (std::size_t x = 0; x < bg.size(); ++x) if (!keep[x]) bg[x] = 'N';
                    return bg;
                };
                target_backgrounds.assign(dup_targets.size(), {});
                for (std::size_t ti = 0; ti < dup_targets.size(); ++ti)
                    if (dup_targets[ti].gene_idx.size() >= 2)
                        for (int gi : dup_targets[ti].gene_idx)
                            target_backgrounds[ti].push_back(background_of(gi, dup_targets[ti].gene_idx));
                // Each gene's coding exons and whole span, for the hybrid reading. A gene beside the copies forms no
                // pair (no exons), but its span still keeps the k-mers it shares out of the pairs' private ones.
                target_exons.assign(dup_targets.size(), {});
                target_spans.assign(dup_targets.size(), {});
                for (std::size_t ti = 0; ti < dup_targets.size(); ++ti) {
                    const DupGeneTarget& t = dup_targets[ti];
                    if (t.gene_idx.size() < 2) continue;
                    target_exons[ti].resize(t.gene_idx.size());
                    target_spans[ti].resize(t.gene_idx.size());
                    for (std::size_t k = 0; k < t.gene_idx.size(); ++k) {
                        const GeneFeature& g = genes[t.gene_idx[k]];
                        // the canonical transcript's extent: a gene record grown by a readthrough holds its neighbour
                        target_spans[ti][k] = g.canonical_cds.empty()
                            ? ref_slice(g.start_1based, g.end_1based)
                            : ref_slice(g.canonical_cds.front().first, g.canonical_cds.back().second);
                        if (!t.in_copies[k]) continue;
                        for (const auto& iv : g.canonical_cds) {
                            std::string piece = ref_slice(iv.first, iv.second);
                            if (!piece.empty()) target_exons[ti][k].push_back(std::move(piece));
                        }
                    }
                }
                hap_seq.assign(graph.paths.size(), std::string());
                run_parallel(graph.paths.size(), options.threads, [&](std::size_t pi) {
                    hap_seq[pi] = spell_path_steps_sequence(graph, graph.paths[pi].steps);
                });
            }

            // Each haplotype's copies of each gene of a target, as reported (-1: none), for naming the genes of
            // events inside a module's copies below.
            std::vector<std::vector<std::vector<int>>> target_cn(dup_targets.size());
            // One row per gene and haplotype: a gene that several records of one site list (two repeats and
            // no module) is written by the first whose copies hold it, else by the first that lists it.
            std::map<std::pair<std::size_t, int>, std::size_t> gene_writer;
            for (int held_only = 1; held_only >= 0; --held_only)
                for (std::size_t ti = 0; ti < dup_targets.size(); ++ti)
                    for (std::size_t k = 0; k < dup_targets[ti].gene_idx.size(); ++k)
                        if (!held_only || dup_targets[ti].in_copies[k])
                            gene_writer.emplace(std::make_pair(dup_targets[ti].bubble_id, dup_targets[ti].gene_idx[k]), ti);
            const auto writes = [&](std::size_t ti, std::size_t k) {
                return gene_writer.at({dup_targets[ti].bubble_id, dup_targets[ti].gene_idx[k]}) == ti;
            };
            for (std::size_t ti = 0; ti < dup_targets.size(); ++ti) {
                const DupGeneTarget& t = dup_targets[ti];
                if (t.gene_idx.size() < 2) {   // single gene: the module total is its copy number
                    if (writes(ti, 0)) emit_total(t, genes[t.gene_idx[0]].gene_name, genes[t.gene_idx[0]].gene_id, "1");
                    continue;
                }
                std::vector<std::string> markers;
                markers.reserve(t.gene_idx.size());
                for (int gi : t.gene_idx) markers.push_back(marker_cache[gi]);
                const std::vector<std::string>& backgrounds = target_backgrounds[ti];
                // The record's copy number per haplotype (-1 where none), used to split a near-identical
                // paralog pair by its per-site allele fraction.
                std::vector<long> total(graph.paths.size(), -1);   // -1: no module total
                for (std::size_t pi = 0; pi < graph.paths.size(); ++pi) {
                    const auto cit = t.sample_cn.find(graph.paths[pi].name);
                    if (cit != t.sample_cn.end()) total[pi] = static_cast<long>(cit->second);
                }
                const std::vector<std::vector<GeneCnEvidence>> ev =
                    resolve_gene_cn(markers, hap_seq, total, 0.5, 31, &backgrounds, &t.in_copies, &target_spans[ti]);
                // Hybrid copies of paralogs in the copies, read exon by exon.
                const std::vector<std::vector<int>> hyb =
                    find_hybrid_copies(target_exons[ti], target_spans[ti], hap_seq);
                const auto hyb_of = [&](std::size_t pi, std::size_t k) {
                    return hyb[pi][k] < 0 ? std::string(".") : std::to_string(hyb[pi][k]);
                };
                target_cn[ti].assign(graph.paths.size(), std::vector<int>(t.gene_idx.size(), -1));
                for (std::size_t pi = 0; pi < graph.paths.size(); ++pi) {
                    const std::string& sname = graph.paths[pi].name;
                    const std::string prefix =
                        std::to_string(t.bubble_id) + '\t' + t.variant_id + '\t' + sname + '\t';
                    for (std::size_t k = 0; k < t.gene_idx.size(); ++k) {
                        GeneCnEvidence e = ev[pi][k];
                        // A hybrid copy is a run of each of its two genes: where the haplotype has one of the gene,
                        // runs do not raise its count.
                        if (hyb[pi][k] > 0) e.cn = e.rounded;
                        const std::string gname = genes[t.gene_idx[k]].gene_name;
                        const std::string& gid = genes[t.gene_idx[k]].gene_id;
                        // Its row is another record's. Only a repeat unit's target shares genes: a module's is the
                        // only target of its site.
                        if (!writes(ti, k)) continue;
                        if (e.separable && total[pi] < 0) {   // no module total: no copy number to report
                            dup_gene_cn_rows.push_back(prefix + gname + "\t.\t0\t.\t.\t.\t.\t.\t" + gid);
                        } else if (e.separable) {
                            char dbuf[32];
                            std::snprintf(dbuf, sizeof(dbuf), "%.2f", e.dosage);
                            // A gene in the copies has no more copies than the module: dosage above it counts
                            // copies elsewhere on the haplotype, so the total is reported and the row unreliable.
                            const bool over = t.in_copies[k] && e.cn > total[pi];
                            target_cn[ti][pi][k] = static_cast<int>(over ? total[pi] : e.cn);
                            dup_gene_cn_rows.push_back(
                                prefix + gname + '\t' + std::to_string(over ? total[pi] : e.cn) + (over ? "\t0\t" : "\t1\t") +
                                dbuf + '\t' + std::to_string(e.hits) + '\t' + std::to_string(e.priv_kmers) + '\t' +
                                (e.copies < 0 ? std::string(".") : std::to_string(e.copies)) + '\t' +
                                hyb_of(pi, k) + '\t' + gid);
                        } else {   // no private k-mers -> indistinguishable from a paralog; report the total
                            // when the gene is in the copies the total counts, else nothing
                            const auto cit = t.sample_cn.find(sname);
                            const std::string total = cit != t.sample_cn.end() && t.in_copies[k] ? std::to_string(cit->second) : ".";
                            if (total != ".") target_cn[ti][pi][k] = static_cast<int>(cit->second);
                            dup_gene_cn_rows.push_back(prefix + gname + '\t' + total + "\t0\t.\t.\t0\t.\t" + hyb_of(pi, k) +
                                                       '\t' + gid);
                        }
                    }
                }
            }

            // An event inside a gene of a module's copies names the genes its carriers' copies hold, not the
            // gene of the reference copy it happens to lie on: that gene gives way to those of its paralogs,
            // itself included, that some carrier has a copy of. An event reaching past the gene (deleting it,
            // say) keeps it: its carriers no longer hold that copy.
            for (std::size_t ti = 0; ti < dup_targets.size(); ++ti) {
                const DupGeneTarget& t = dup_targets[ti];
                if (!t.module || target_cn[ti].empty()) continue;
                std::vector<std::string> seqs;
                for (int gi : t.gene_idx) seqs.push_back(marker_cache[gi]);
                const std::vector<int> group = paralog_groups(seqs);
                std::unordered_map<std::string, std::size_t> copy_gene;   // gene name -> index, genes in the copies
                for (std::size_t k = 0; k < t.gene_idx.size(); ++k)
                    if (t.in_copies[k]) copy_gene.emplace(vcf_escape(genes[t.gene_idx[k]].gene_name), k);   // as GENES writes it
                for (OutRecord& rec : out_records) {
                    if (rec.bubble_id != t.bubble_id || rec.id == t.variant_id) continue;
                    std::vector<std::string> f;
                    {
                        std::string line = rec.line;
                        if (!line.empty() && line.back() == '\n') line.pop_back();
                        std::stringstream ls(line);
                        for (std::string x; std::getline(ls, x, '\t');) f.push_back(x);
                    }
                    if (f.size() < 9 + graph.paths.size() || f[7].find("SVTYPE=CNV") != std::string::npos) continue;
                    const std::size_t g0 = f[7].find("GENES=");
                    if (g0 == std::string::npos || (g0 > 0 && f[7][g0 - 1] != ';')) continue;
                    const std::size_t g1 = std::min(f[7].find(';', g0), f[7].size());
                    std::vector<std::string> names;
                    {
                        std::stringstream gs(f[7].substr(g0 + 6, g1 - g0 - 6));
                        for (std::string x; std::getline(gs, x, ',');) names.push_back(x);
                    }
                    std::size_t pos = 0, end = 0;
                    try {
                        pos = std::stoull(f[1]);
                        const std::size_t e0 = f[7].find("END=");
                        end = e0 != std::string::npos && (e0 == 0 || f[7][e0 - 1] == ';') ? std::stoull(f[7].substr(e0 + 4)) : pos;
                    } catch (...) { continue; }
                    const auto inside = [&](std::size_t k) {
                        const GeneFeature& gf = genes[t.gene_idx[k]];
                        return gf.start_1based <= pos + 1 && end <= gf.end_1based;
                    };
                    std::set<int> groups;
                    for (const std::string& nm : names)
                        if (const auto it = copy_gene.find(nm); it != copy_gene.end() && inside(it->second))
                            groups.insert(group[it->second]);
                    if (groups.empty()) continue;
                    std::set<std::string> held;   // the group's genes some carrier has a copy of
                    for (std::size_t pi = 0; pi < graph.paths.size(); ++pi) {
                        if (f[9 + pi].empty() || f[9 + pi][0] != '1') continue;
                        for (const auto& [nm, k] : copy_gene)
                            if (groups.count(group[k]) && target_cn[ti][pi][k] >= 1) held.insert(nm);
                    }
                    if (held.empty()) continue;   // no carrier's copies are known: the interval's genes stand
                    std::set<std::string> out;
                    for (const std::string& nm : names)
                        if (const auto it = copy_gene.find(nm); it == copy_gene.end() || !inside(it->second))
                            out.insert(nm);
                    out.insert(held.begin(), held.end());
                    std::string joined;
                    for (const std::string& nm : out) joined += (joined.empty() ? "" : ",") + nm;
                    f[7] = f[7].substr(0, g0 + 6) + joined + f[7].substr(g1);
                    std::string line;
                    for (std::size_t c = 0; c < f.size(); ++c) line += (c ? "\t" : "") + f[c];
                    rec.line = line + '\n';
                }
            }
        }

        const std::string dg_path = options.out_prefix + ".dup_gene_cn.tsv";
        std::ofstream dg_out(staged.stage(dg_path));
        if (!dg_out) {
            throw std::runtime_error("Failed to write per-gene DUP CN TSV: " + dg_path);
        }
        dg_out << "bubble_id\tvariant_id\tsample\tgenes\tcn\treliable\tdosage\thits\tpriv_kmers\tcopies\thybrid_copies\tgene_id\n";
        for (const std::string& r : dup_gene_cn_rows) dg_out << r << '\n';
        dg_out.close();
        if (!dg_out) throw std::runtime_error("Failed to finalize per-gene DUP CN TSV: " + dg_path);
    }

    // Written after the gene step, which can name the genes of events inside module copies.
    const std::string region_path = staged.stage(options.out_prefix + ".region.vcf");
    std::ofstream region_out(region_path);
    if (!region_out) {
        throw std::runtime_error("Failed to write region VCF: " + options.out_prefix + ".region.vcf");
    }
    write_vcf_header(region_out);

    for (const OutRecord& rec : out_records) emit(region_out, rec.line);
    region_out.close();
    if (!region_out) {
        throw std::runtime_error("Failed to finalize region VCF: " + options.out_prefix + ".region.vcf");
    }
    // Each site's file written whole and closed before the next is opened, so a run with many sites holds one
    // file open at a time, not one per site.
    if (options.write_per_bubble_vcf) {
        std::map<std::size_t, std::vector<const OutRecord*>> by_site;   // in coordinate order within a site
        for (const OutRecord& rec : out_records) if (rec.bubble_id != 0) by_site[rec.bubble_id].push_back(&rec);
        for (const auto& [bid, recs] : by_site) {
            const std::string path = options.out_prefix + ".bubble_" + std::to_string(bid) + ".vcf";
            std::ofstream f(staged.stage(path));
            if (!f) throw std::runtime_error("Failed to write bubble VCF: " + path + ": " + std::strerror(errno));
            write_vcf_header(f);
            for (const OutRecord* rec : recs) emit(f, rec->line);
            f.close();
            if (!f) throw std::runtime_error("Failed to finalize bubble VCF: " + path);
        }
        summary.per_bubble_vcfs = by_site.size();
    }

    // One commit for the whole family -- VCFs, variant nodes and the GTF sidecars together.
    staged.commit();
    summary.gene_tables_written = gtf_active;
    // Gene tables left at this prefix by an earlier --gtf run describe that run, not this one, when
    // this run wrote none (no --gtf, or no gene in the region). Dropped only after the commit.
    if (!gtf_active) {
        for (const char* suffix : {".node_genes.tsv", ".dup_gene_cn.tsv"}) {
            const std::string stale = options.out_prefix + suffix;
            if (is_input(stale)) continue;
            std::error_code ec;
            if (std::filesystem::remove(stale, ec) && !options.quiet)
                std::cerr << cli::log_prefix("call") << "removed " << stale << ", left by an earlier run at this prefix\n";
        }
    }
    // The same for the other optional outputs: variant nodes or carrier events from an earlier run with
    // other options would sit beside this run's region VCF and be read as current.
    {
        std::vector<std::string> stale;
        if (!options.write_variant_nodes) {
            stale.push_back(options.out_prefix + ".variant_nodes.tsv");
            stale.push_back(options.out_prefix + ".carrier_events.tsv");
        }
        for (const std::string& f : stale) {
            if (is_input(f)) continue;
            std::error_code ec;
            if (std::filesystem::remove(f, ec) && !options.quiet)
                std::cerr << cli::log_prefix("call") << "removed " << f << ", left by an earlier run at this prefix\n";
        }
    }
    // Only now are the superseded files safe to drop: until the new outputs are installed, removing
    // them would leave a failed run having destroyed the previous good result.
    if (!obsolete_bubble_vcfs.empty()) {
        std::size_t removed = 0;
        for (const std::filesystem::path& p : obsolete_bubble_vcfs) {
            if (is_input(p.string())) continue;
            std::error_code ec;
            if (std::filesystem::remove(p, ec) && !ec) ++removed;
            else if (ec) {
                std::cerr << cli::log_prefix("call") << "warning: could not remove stale per-bubble VCF " << p.string()
                          << ": " << ec.message() << " (it does not belong to this run)\n";
            }
        }
        if (removed > 0 && !options.quiet) {
            std::cerr << cli::log_prefix("call") << "removed " << removed
                      << " stale per-bubble VCF(s) from an earlier run at this prefix\n";
        }
    }

    if (summary_out) *summary_out = summary;
    // The command layer reports this summary via RunLog; no duplicate print here.
}

} // namespace panvar
