#include "panvar/transcripts.hpp"

#include "edlib.h"

#include "panvar/cli_utils.hpp"
#include "panvar/graph_utils.hpp"
#include "panvar/parallel.hpp"
#include "panvar/ref_path.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <set>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <utility>

namespace panvar {

const char* hst_status_name(HstStatus s) {
    switch (s) {
        case HstStatus::Ok: return "ok";
        case HstStatus::SpliceSiteChanged: return "splice_site_changed";
        case HstStatus::BoundaryMoved: return "boundary_moved";
        case HstStatus::ClippedAtEdge: return "clipped_at_edge";
        case HstStatus::IntronLost: return "intron_lost";
        case HstStatus::NotAsAnnotated: return "not_as_annotated";
        case HstStatus::MissingExon: return "missing_exon";
        case HstStatus::PartialCopy: return "partial_copy";
        case HstStatus::NoCopy: return "no_copy";
    }
    return "ok";
}

namespace {

// Chaining: how far back a link may reach, in anchors and in bases (or twice the gene's stretch; any
// length across a tandem stretch), and what an indel costs per base.
constexpr std::size_t kChainLookback = 200;
constexpr std::size_t kChainMaxGap = 100000;
constexpr double kChainIndelCost = 0.05;
// The largest gap between two chained nodes filled by an exact node alignment (cells); larger ones are
// filled from both ends only.
constexpr std::size_t kFillMaxCells = 4000000;
// Passes of two genes taking the same stretch of a walk: overlapping by more than this share of the shorter.
constexpr double kParalogOverlap = 0.5;
// A copy placed by sequence: the marker's exonic bases, aligned at the identity asked for, make a copy where
// they are covered to the share and a partial copy where covered to this much. The window reaches this share
// of the stretch past each of its ends, and is given up past this many times the stretch: the chain's
// anchors then disagree on where the gene lies.
constexpr double kPartialCoverage = 0.2;
constexpr double kWindowMargin = 0.1;
constexpr double kWindowMax = 3.0;
// An aligned base counts where the bases around it (this many on each side) match to this share, deletions of
// this length or more aside: unit-cost alignment sets a stretch the walk lacks at an end of the gene against
// whatever follows, base against base, with many short gaps.
constexpr std::size_t kLocalFlank = 15;
constexpr double kLocalIdentity = 0.7;
constexpr std::size_t kLongDeletion = 20;
// Fewer aligned bases than this between two gaps are part of the indel around them: an alignment may cut one
// indel in two around a base that happens to match.
constexpr std::size_t kBetweenGaps = 12;
// The largest alignments (stretch times window, in cells) tried: of the gene's span, within a band of a
// quarter of its bases, and of the marker's span, at any distance.
constexpr double kSequenceMaxCells = 2.5e9;
constexpr double kSequenceMaxCellsAny = 5e8;

struct NodeTable {
    std::vector<std::string> seq;
    std::unordered_map<std::string, int> id;
};

struct Walk {
    std::vector<int> node;
    std::vector<char> rev;
    std::vector<std::size_t> pref;   // bases before each step; size steps + 1
    std::size_t size() const { return node.size(); }
    std::size_t len() const { return pref.back(); }
};

Walk to_walk(const PathRecord& p, const NodeTable& t) {
    Walk w;
    w.node.reserve(p.steps.size());
    w.rev.reserve(p.steps.size());
    w.pref.reserve(p.steps.size() + 1);
    w.pref.push_back(0);
    for (const PathStep& s : p.steps) {
        const int v = t.id.at(s.node_id);
        w.node.push_back(v);
        w.rev.push_back(s.reverse ? 1 : 0);
        w.pref.push_back(w.pref.back() + t.seq[v].size());
    }
    return w;
}

// The same walk read from its end: steps in reverse order, each in the other orientation.
Walk reversed_walk(const Walk& w, const NodeTable& t) {
    Walk r;
    const std::size_t n = w.size();
    r.node.resize(n);
    r.rev.resize(n);
    r.pref.assign(n + 1, 0);
    for (std::size_t k = 0; k < n; ++k) {
        r.node[k] = w.node[n - 1 - k];
        r.rev[k] = w.rev[n - 1 - k] ? 0 : 1;
        r.pref[k + 1] = r.pref[k] + t.seq[r.node[k]].size();
    }
    return r;
}

std::string spell(const Walk& w, const NodeTable& t) {
    std::string s;
    s.reserve(w.len());
    for (std::size_t k = 0; k < w.size(); ++k) s += w.rev[k] ? reverse_complement(t.seq[w.node[k]]) : t.seq[w.node[k]];
    return s;
}

// Where each node occurs on a walk, as one sorted list per node.
struct Occurrences {
    std::vector<std::size_t> off;
    std::vector<std::size_t> pos;
    Occurrences(const Walk& w, std::size_t nodes) : off(nodes + 1, 0), pos(w.size()) {
        for (const int v : w.node) ++off[static_cast<std::size_t>(v) + 1];
        for (std::size_t v = 0; v < nodes; ++v) off[v + 1] += off[v];
        std::vector<std::size_t> fill(off.begin(), off.end() - 1);
        for (std::size_t k = 0; k < w.size(); ++k) pos[fill[static_cast<std::size_t>(w.node[k])]++] = k;
    }
};

using Cuts = std::vector<std::pair<std::size_t, std::size_t>>;   // 0-based, half-open, on the reference walk

struct Gene {
    std::vector<std::size_t> tx;
    std::size_t marker = 0;           // transcript whose exons decide copies
    std::size_t lo = 0, hi = 0;       // span on the reference walk
    std::size_t m_lo = 0, m_hi = 0;   // the marker's span: where the gene is, for telling copies apart
    std::size_t gs = 0, ge = 0;       // reference steps [gs, ge) holding the span, one more on each side
    std::vector<char> uniq;           // per step of [gs, ge): its node occurs once there
    std::size_t marker_bp = 0;
    std::vector<std::pair<std::size_t, std::size_t>> tandem_spans;   // its transcripts' tandem stretches
};

// One pass of a walk through a gene: per reference step of [gs, ge), the walk step it matches or -1.
struct Pass {
    std::size_t gene = 0;
    bool reverse = false;
    std::vector<long> match;
    std::size_t first = 0, last = 0;  // first and last matched steps, relative to gs
    std::size_t h_lo = 0, h_hi = 0;   // bases of the walk as written
    std::vector<std::size_t> walk_steps;   // the matched steps of the walk as written, ascending
    std::size_t marker_matched = 0;
    std::size_t matched = 0;          // reference bases of the gene's stretch matched, introns included
    double share = 0.0;
    bool at_walk_start = false, at_walk_end = false;   // the walk (oriented) starts or ends inside it
};

struct Anchor {
    std::size_t i, j;          // reference step, walk step
    std::size_t x, y, w;       // their bases, and the node's length
};

std::size_t overlap(std::size_t a0, std::size_t a1, std::size_t b0, std::size_t b1) {
    const std::size_t lo = std::max(a0, b0), hi = std::min(a1, b1);
    return hi > lo ? hi - lo : 0;
}

std::size_t marker_overlap(const Cuts& exons, std::size_t a, std::size_t b) {
    std::size_t s = 0;
    for (const auto& [lo, hi] : exons) s += overlap(a, b, lo, hi);
    return s;
}

// Weighted longest common subsequence of node steps (same node, same orientation; weight its length) between
// reference steps [i0, i1) and walk steps [j0, j1). With `late`, ties go to the latest walk steps.
void fill_gap(const Walk& R, std::size_t i0, std::size_t i1, const Walk& W, std::size_t j0, std::size_t j1,
              const NodeTable& t, std::vector<long>& match, std::size_t gs, bool late) {
    if (i1 <= i0 || j1 <= j0) return;
    const std::size_t na = i1 - i0, nb = j1 - j0;
    const auto eq = [&](std::size_t i, std::size_t j) { return R.node[i] == W.node[j] && R.rev[i] == W.rev[j]; };
    if ((na + 1) * (nb + 1) > kFillMaxCells) {   // too large: match equal steps inward from both ends
        std::size_t a = i0, b = j0;
        while (a < i1 && b < j1 && eq(a, b)) { match[a - gs] = static_cast<long>(b); ++a; ++b; }
        std::size_t c = i1, d = j1;
        while (c > a && d > b && eq(c - 1, d - 1)) { --c; --d; match[c - gs] = static_cast<long>(d); }
        return;
    }
    // late: align the reversed ranges, so the traceback (which favours the last cells) favours early ones.
    const auto ri = [&](std::size_t a) { return late ? i0 + a : i1 - 1 - a; };
    const auto rj = [&](std::size_t b) { return late ? j0 + b : j1 - 1 - b; };
    std::vector<std::size_t> dp((na + 1) * (nb + 1), 0);
    const auto at = [&](std::size_t a, std::size_t b) -> std::size_t& { return dp[a * (nb + 1) + b]; };
    for (std::size_t a = 1; a <= na; ++a)
        for (std::size_t b = 1; b <= nb; ++b) {
            std::size_t v = std::max(at(a - 1, b), at(a, b - 1));
            if (eq(ri(a - 1), rj(b - 1))) v = std::max(v, at(a - 1, b - 1) + t.seq[R.node[ri(a - 1)]].size());
            at(a, b) = v;
        }
    std::size_t a = na, b = nb;
    while (a > 0 && b > 0) {
        if (eq(ri(a - 1), rj(b - 1)) && at(a, b) == at(a - 1, b - 1) + t.seq[R.node[ri(a - 1)]].size()) {
            match[ri(a - 1) - gs] = static_cast<long>(rj(b - 1));
            --a; --b;
        } else if (at(a, b) == at(a, b - 1)) {
            --b;
        } else {
            --a;
        }
    }
}

// Colinear chains of the gene's once-occurring reference nodes on walk W, best first; each anchor used once.
std::vector<std::vector<Anchor>> chain_anchors(const Walk& R, const Gene& g, const Walk& W, const Occurrences& occ,
                                               bool reverse, std::size_t walk_steps) {
    std::vector<Anchor> an;
    for (std::size_t i = g.gs; i < g.ge; ++i) {
        if (!g.uniq[i - g.gs]) continue;
        const std::size_t v = static_cast<std::size_t>(R.node[i]);
        for (std::size_t q = occ.off[v]; q < occ.off[v + 1]; ++q) {
            const std::size_t p = occ.pos[q];
            const std::size_t j = reverse ? walk_steps - 1 - p : p;
            if (W.rev[j] != R.rev[i]) continue;
            an.push_back({i, j, R.pref[i], W.pref[j], R.pref[i + 1] - R.pref[i]});
        }
    }
    std::sort(an.begin(), an.end(), [](const Anchor& a, const Anchor& b) { return a.j < b.j; });
    const std::size_t n = an.size();
    std::vector<double> f(n);
    std::vector<long> prev(n, -1);
    // A gene with a tandem stretch links across any length: haplotypes carry any number of units.
    const std::size_t max_gap = g.tandem_spans.empty() ? std::max(kChainMaxGap, 2 * (R.pref[g.ge] - R.pref[g.gs]))
                                                       : SIZE_MAX;
    for (std::size_t a = 0; a < n; ++a) {
        f[a] = static_cast<double>(an[a].w);
        const std::size_t stop = a > kChainLookback ? a - kChainLookback : 0;
        for (std::size_t b = a; b-- > stop;) {
            if (an[b].i >= an[a].i) continue;
            const std::size_t dy = an[a].y - (an[b].y + an[b].w);
            if (dy > max_gap) break;
            const std::size_t dx = an[a].x - (an[b].x + an[b].w);
            if (dx > max_gap) continue;
            // A link across a tandem stretch pays nothing for the length the units add or remove.
            bool across = false;
            for (const auto& [t0, t1] : g.tandem_spans) across |= an[b].x + an[b].w < t1 && t0 < an[a].x;
            const double dd = across ? 0.0 : std::fabs(static_cast<double>(dx) - static_cast<double>(dy));
            const double s = f[b] + static_cast<double>(an[a].w) - kChainIndelCost * dd - (dd > 0 ? std::log2(dd) : 0.0);
            if (s > f[a]) { f[a] = s; prev[a] = static_cast<long>(b); }
        }
    }
    std::vector<std::size_t> order(n);
    std::iota(order.begin(), order.end(), std::size_t{0});
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return f[a] != f[b] ? f[a] > f[b] : a < b; });
    std::vector<char> used(n, 0);
    std::vector<std::vector<Anchor>> chains;
    for (const std::size_t end : order) {
        if (used[end]) continue;
        std::vector<Anchor> c;
        for (long k = static_cast<long>(end); k >= 0 && !used[static_cast<std::size_t>(k)]; k = prev[static_cast<std::size_t>(k)]) {
            used[static_cast<std::size_t>(k)] = 1;
            c.push_back(an[static_cast<std::size_t>(k)]);
        }
        std::reverse(c.begin(), c.end());
        chains.push_back(std::move(c));
    }
    return chains;
}

// A chain made into a pass: the gaps between its anchors and past its ends aligned node by node.
Pass make_pass(const Walk& R, const Gene& g, const Cuts& marker, const Walk& W, const std::vector<Anchor>& c,
               const NodeTable& t, std::size_t gene_index, bool reverse) {
    Pass p;
    p.gene = gene_index;
    p.reverse = reverse;
    p.match.assign(g.ge - g.gs, -1);
    for (std::size_t k = 0; k < c.size(); ++k) {
        p.match[c[k].i - g.gs] = static_cast<long>(c[k].j);
        if (k + 1 < c.size()) fill_gap(R, c[k].i + 1, c[k + 1].i, W, c[k].j + 1, c[k + 1].j, t, p.match, g.gs, true);
    }
    // Past the ends, the walk is searched over a little more than the reference stretch left to match.
    const std::size_t lead = R.pref[c.front().i] - R.pref[g.gs];
    std::size_t j0 = c.front().j;
    while (j0 > 0 && W.pref[c.front().j] - W.pref[j0 - 1] <= lead + lead / 2 + 1000) --j0;
    fill_gap(R, g.gs, c.front().i, W, j0, c.front().j, t, p.match, g.gs, true);
    const std::size_t trail = R.pref[g.ge] - R.pref[c.back().i + 1];
    std::size_t j1 = c.back().j + 1;
    while (j1 < W.size() && W.pref[j1 + 1] - W.pref[c.back().j + 1] <= trail + trail / 2 + 1000) ++j1;
    fill_gap(R, c.back().i + 1, g.ge, W, c.back().j + 1, j1, t, p.match, g.gs, false);

    bool any = false;
    for (std::size_t k = 0; k < p.match.size(); ++k) {
        if (p.match[k] < 0) continue;
        p.walk_steps.push_back(reverse ? W.size() - 1 - static_cast<std::size_t>(p.match[k]) : static_cast<std::size_t>(p.match[k]));
        if (!any) p.first = k;
        p.last = k;
        any = true;
        const std::size_t i = g.gs + k;
        p.marker_matched += marker_overlap(marker, R.pref[i], R.pref[i + 1]);
        p.matched += R.pref[i + 1] - R.pref[i];
    }
    std::sort(p.walk_steps.begin(), p.walk_steps.end());
    p.share = g.marker_bp ? static_cast<double>(p.marker_matched) / static_cast<double>(g.marker_bp) : 0.0;
    // The stretch of the walk the pass takes: its matched steps inside the marker's span, each cut to that
    // span. Flanking steps may be matched across other copies, transcripts reaching into a neighbour would
    // stretch it over that, and a long node can hold other genes too.
    std::size_t lo = SIZE_MAX, hi = 0;
    for (std::size_t k = p.first; k <= p.last; ++k) {
        const std::size_t i = g.gs + k;
        if (p.match[k] < 0 || R.pref[i + 1] <= g.m_lo || R.pref[i] >= g.m_hi) continue;
        const std::size_t w = W.pref[static_cast<std::size_t>(p.match[k])];
        lo = std::min(lo, w + std::max(R.pref[i], g.m_lo) - R.pref[i]);
        hi = std::max(hi, w + std::min(R.pref[i + 1], g.m_hi) - R.pref[i]);
    }
    if (lo >= hi) { lo = W.pref[static_cast<std::size_t>(p.match[p.first])]; hi = lo; }
    p.h_lo = reverse ? W.len() - hi : lo;
    p.h_hi = reverse ? W.len() - lo : hi;
    // The walk starts or ends inside the gene: less of it lies past the pass than the reference has there.
    const std::size_t w0 = W.pref[static_cast<std::size_t>(p.match[p.first])];
    const std::size_t w1 = W.pref[static_cast<std::size_t>(p.match[p.last]) + 1];
    p.at_walk_start = w0 < R.pref[g.gs + p.first] - R.pref[g.gs] + 1;
    p.at_walk_end = W.len() - w1 < R.pref[g.ge] - R.pref[g.gs + p.last + 1] + 1;
    return p;
}

struct CutOnWalk {
    bool ok = false;
    std::size_t pos = 0;   // on the oriented walk
    long shift = 0;        // reference bases the boundary moved by, positive when the exon grows
    bool in_variant = false;   // kept at its offset inside a substitution of the reference's length
};

// Carry reference cut c (between bases) to the walk through pass p. `start`: c opens an exon.
CutOnWalk map_cut(const Walk& R, const Gene& g, const Pass& p, const Walk& W, std::size_t c, bool start) {
    CutOnWalk out;
    // The step holding base c, or the step after the span when c ends the reference walk.
    const std::size_t i = static_cast<std::size_t>(std::upper_bound(R.pref.begin(), R.pref.end(), c) - R.pref.begin()) - 1;
    const std::size_t off = c - R.pref[i];
    const auto m = [&](std::size_t k) -> long { return k >= g.gs && k < g.ge ? p.match[k - g.gs] : -1; };
    if (off > 0 && m(i) >= 0) {
        out.ok = true;
        out.pos = W.pref[static_cast<std::size_t>(m(i))] + off;
        return out;
    }
    if (off == 0) {
        const long left = i > 0 ? m(i - 1) : -1, right = m(i);
        if (right >= 0 && (start || left < 0)) { out.ok = true; out.pos = W.pref[static_cast<std::size_t>(right)]; return out; }
        if (left >= 0) { out.ok = true; out.pos = W.pref[static_cast<std::size_t>(left) + 1]; return out; }
    }
    // Inside a variant: the nearest shared positions on either side.
    long kl = -1, kr = -1;
    for (std::size_t k = i; k-- > g.gs;) if (m(k) >= 0) { kl = static_cast<long>(k); break; }
    for (std::size_t k = i + 1; k < g.ge; ++k) if (m(k) >= 0) { kr = static_cast<long>(k); break; }
    if (kl < 0 && kr < 0) return out;
    const std::size_t xl = kl >= 0 ? R.pref[static_cast<std::size_t>(kl) + 1] : 0;
    const std::size_t xr = kr >= 0 ? R.pref[static_cast<std::size_t>(kr)] : 0;
    if (kl >= 0 && kr >= 0) {
        // The walk's unmatched stretch is as long as the reference's: a substitution (clustered SNPs the graph
        // folds into one node), where every base keeps its place and the boundary with it.
        const std::size_t wl = W.pref[static_cast<std::size_t>(m(static_cast<std::size_t>(kl))) + 1];
        const std::size_t wr = W.pref[static_cast<std::size_t>(m(static_cast<std::size_t>(kr)))];
        if (wl <= wr && wr - wl == xr - xl) {
            out.ok = true;
            out.in_variant = true;
            out.pos = wl + (c - xl);
            return out;
        }
    }
    bool go_left;
    if (kl < 0) go_left = false;
    else if (kr < 0) go_left = true;
    else if (c - xl != xr - c) go_left = c - xl < xr - c;
    else go_left = start;   // a tie keeps the variant inside the exon
    out.ok = true;
    if (go_left) {
        out.pos = W.pref[static_cast<std::size_t>(m(static_cast<std::size_t>(kl))) + 1];
        const long d = static_cast<long>(c - xl);
        out.shift = start ? d : -d;
    } else {
        out.pos = W.pref[static_cast<std::size_t>(m(static_cast<std::size_t>(kr)))];
        const long d = static_cast<long>(xr - c);
        out.shift = start ? -d : d;
    }
    return out;
}

bool exon_shared(const Walk& R, const Gene& g, const Pass& p, std::size_t cs, std::size_t ce) {
    for (std::size_t k = 0; k < p.match.size(); ++k) {
        if (p.match[k] < 0) continue;
        const std::size_t i = g.gs + k;
        if (R.pref[i] < ce && R.pref[i + 1] > cs) return true;
    }
    return false;
}

void add_detail(std::string& d, const std::string& note) {
    if (note.empty()) return;
    if (!d.empty()) d += ';';
    d += note;
}

std::string signed_text(long v) { return (v > 0 ? "+" : "") + std::to_string(v); }

// Exons of a transcript that repeat with a tandem unit: exons starting at the same node and offset of the
// reference walk, as a unit the reference walks several times holds the same exon in each pass. The stretch
// from the first to the last such exon is projected unit by unit, by sequence: the reference's window
// between the exons around the stretch is cut into a head, its units (each starting where the walk passes
// the stretch's first exon start, annotated or not) and a tail; on a haplotype the head and tail are
// aligned at the ends of its window and the units tile the rest, each aligned to the reference unit it
// resembles most, whose annotated exons (perhaps none) are lifted through that alignment.
struct Tandem {
    bool any = false;
    std::size_t lo = 0, hi = 0;   // genomic exon indices of the stretch, inclusive
    std::size_t rl = 0, rh = 0;   // the reference window: from the exon before the stretch to the one after
    std::size_t first = 0;        // where its first unit starts
    // where the stretch starts or ends the transcript: the window reaches the nearest node the gene's
    // stretch of the reference holds once (its step), which the haplotype's pass may share
    std::size_t lo_anchor = SIZE_MAX, hi_anchor = SIZE_MAX;
    std::string head, tail;       // the window's sequence before the first unit and after the last
    struct Unit {
        std::string seq;
        std::vector<std::array<std::size_t, 3>> exons;   // in the unit: start, end, the reference exon
        std::vector<std::uint64_t> kmers;               // distinct 15-mers, sorted
    };
    std::vector<Unit> units;
    std::size_t min_len = 0;
    std::vector<std::pair<std::string, bool>> probes;   // the units' distinct first and last 200 bases (true: first)
};

std::size_t step_at(const Walk& w, std::size_t c) {
    return static_cast<std::size_t>(std::upper_bound(w.pref.begin(), w.pref.end(), c) - w.pref.begin()) - 1;
}

constexpr std::size_t kUnitK = 15;

std::vector<std::uint64_t> unit_kmers(const std::string& s) {
    std::vector<std::uint64_t> v;
    std::uint64_t code = 0, valid = 0;
    const std::uint64_t mask = (1ULL << (2 * kUnitK)) - 1;
    for (const char c : s) {
        const int x = c == 'A' ? 0 : c == 'C' ? 1 : c == 'G' ? 2 : c == 'T' ? 3 : -1;
        if (x < 0) { valid = 0; continue; }
        code = ((code << 2) | static_cast<std::uint64_t>(x)) & mask;
        if (++valid >= kUnitK) v.push_back(code);
    }
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
    return v;
}

double jaccard(const std::vector<std::uint64_t>& a, const std::vector<std::uint64_t>& b) {
    std::size_t common = 0;
    for (std::size_t i = 0, j = 0; i < a.size() && j < b.size();) {
        if (a[i] < b[j]) ++i;
        else if (b[j] < a[i]) ++j;
        else { ++common; ++i; ++j; }
    }
    const std::size_t all = a.size() + b.size() - common;
    return all ? static_cast<double>(common) / static_cast<double>(all) : 0.0;
}

Tandem find_tandem(const Walk& R, const std::string& rseq, const Cuts& cuts, const Gene& g) {
    Tandem t;
    using Key = std::tuple<int, char, std::size_t>;
    std::vector<Key> start(cuts.size());
    std::map<Key, std::size_t> starts;
    for (std::size_t k = 0; k < cuts.size(); ++k) {
        const std::size_t i = step_at(R, cuts[k].first);
        start[k] = Key{R.node[i], R.rev[i], cuts[k].first - R.pref[i]};
        ++starts[start[k]];
    }
    for (std::size_t k = 0; k < cuts.size(); ++k) {
        if (starts[start[k]] < 2) continue;
        if (!t.any) t.lo = k;
        t.hi = k;
        t.any = true;
    }
    if (!t.any) return t;
    const std::size_t n = cuts.size();
    t.rl = t.lo > 0 ? cuts[t.lo - 1].second : cuts[t.lo].first;
    t.rh = t.hi + 1 < n ? cuts[t.hi + 1].first : cuts[t.hi].second;
    if (t.lo == 0)
        for (std::size_t i = step_at(R, t.rl); i-- > g.gs;)
            if (g.uniq[i - g.gs]) { t.lo_anchor = i; t.rl = R.pref[i + 1]; break; }
    if (t.hi + 1 == n)
        for (std::size_t i = step_at(R, t.rh - 1) + 1; i < g.ge; ++i)
            if (g.uniq[i - g.gs]) { t.hi_anchor = i; t.rh = R.pref[i]; break; }
    // the units: where the walk passes the first exon's start inside the window
    const auto [node, rev, off] = start[t.lo];
    std::vector<std::size_t> occ;
    for (std::size_t i = step_at(R, t.rl); i < R.size() && R.pref[i] < t.rh; ++i)
        if (R.node[i] == node && R.rev[i] == rev && R.pref[i] + off >= t.rl && R.pref[i] + off < t.rh) occ.push_back(R.pref[i] + off);
    std::vector<std::size_t> lens;
    for (std::size_t u = 0; u + 1 < occ.size(); ++u) lens.push_back(occ[u + 1] - occ[u]);
    std::sort(lens.begin(), lens.end());
    // the last unit runs a typical unit's length, and at least to the stretch's last exon
    const std::size_t typical = lens.empty() ? t.rh - occ.back() : lens[lens.size() / 2];
    const std::size_t last_end = std::min(t.rh, std::max(occ.back() + typical, cuts[t.hi].second));
    t.first = occ.front();
    t.head = rseq.substr(t.rl, occ.front() - t.rl);
    t.tail = rseq.substr(last_end, t.rh - last_end);
    for (std::size_t u = 0; u < occ.size(); ++u) {
        Tandem::Unit un;
        const std::size_t a = occ[u], b = u + 1 < occ.size() ? occ[u + 1] : last_end;
        un.seq = rseq.substr(a, b - a);
        for (std::size_t k = t.lo; k <= t.hi; ++k)
            if (cuts[k].first >= a && cuts[k].first < b)
                un.exons.push_back({cuts[k].first - a, std::min(cuts[k].second, b) - a, k});
        un.kmers = unit_kmers(un.seq);
        t.units.push_back(std::move(un));
    }
    t.min_len = SIZE_MAX;
    for (const auto& un : t.units) t.min_len = std::min(t.min_len, un.seq.size());
    for (const bool at_start : {true, false})
        for (const auto& un : t.units) {
            const std::size_t len = std::min<std::size_t>(200, un.seq.size());
            std::pair<std::string, bool> pr{at_start ? un.seq.substr(0, len) : un.seq.substr(un.seq.size() - len), at_start};
            if (std::find(t.probes.begin(), t.probes.end(), pr) == t.probes.end()) t.probes.push_back(std::move(pr));
        }
    return t;
}

struct Piece {
    std::size_t a = 0, b = 0;   // on the oriented walk
    std::size_t ref = 0;        // the reference exon it stands for
};

struct Placed {
    std::vector<Piece> pieces;
    std::size_t units = 0, resized = 0, unplaced = 0;
};

// edlib on two sequences: the edit distance and, with `path`, the alignment; -1 past `max_ed`.
EdlibAlignResult align_edlib(const std::string& q, const std::string& t, EdlibAlignMode mode, bool path, int max_ed) {
    return edlibAlign(q.data(), static_cast<int>(q.size()), t.data(), static_cast<int>(t.size()),
                      edlibNewAlignConfig(max_ed, mode, path ? EDLIB_TASK_PATH : EDLIB_TASK_LOC, nullptr, 0));
}

// Where query `q` ends when aligned from the start of `t` (its whole length, t's end free): the length of t
// it takes, or none past an edit distance of a third of q.
std::size_t prefix_end(const std::string& q, const std::string& t) {
    if (q.empty() || t.empty()) return 0;
    EdlibAlignResult r = align_edlib(q, t, EDLIB_MODE_SHW, false, static_cast<int>(std::max<std::size_t>(20, q.size() / 3)));
    const std::size_t end = r.editDistance >= 0 && r.numLocations > 0 ? static_cast<std::size_t>(r.endLocations[0]) + 1 : SIZE_MAX;
    edlibFreeAlignResult(r);
    return end;
}

// Where each base of `q` lies in `t`, from their alignment (global, or with `infix` q inside t with t's ends
// free): per base of q, the position in t before it and after it (equal where t lacks the base), and whether
// the alignment sets it against a base of t, equal or not. Unit-cost alignment places an indel's bases
// anywhere in a stretch of equal cost, so with `solid` the alignment is kept only in runs of at least 12
// matched bases; between two runs, the shorter side is set against the other's left or right end, whichever
// matches more bases, and the rest of the longer side is one indel. Without, every aligned base is kept
// where the alignment sets it: a stretch of clustered differences, with no run to keep, stays base for base.
// The alignment is tried within a band wide enough for the length difference and a quarter of the bases more,
// then, with `any_distance`, without one; else it fails past the band.
struct Lift {
    bool ok = false;
    std::vector<std::size_t> before, after;
    std::vector<char> aligned, matched;
    std::size_t t0 = 0, t1 = 0;   // the stretch of t the alignment covers
};

Lift lift(const std::string& q, const std::string& t, bool infix = false, bool solid = true, bool any_distance = true) {
    constexpr std::size_t kSolid = 12;
    const EdlibAlignMode mode = infix ? EDLIB_MODE_HW : EDLIB_MODE_NW;
    const std::size_t diff = infix ? 0 : q.size() > t.size() ? q.size() - t.size() : t.size() - q.size();
    EdlibAlignResult r = align_edlib(q, t, mode, true, static_cast<int>(diff + q.size() / 4 + 100));
    Lift out;
    if (r.editDistance < 0) {
        edlibFreeAlignResult(r);
        if (!any_distance) return out;
        r = align_edlib(q, t, mode, true, -1);
    }
    out.ok = true;
    out.aligned.assign(q.size(), 0);
    out.matched.assign(q.size(), 0);
    out.t0 = infix && r.numLocations > 0 ? static_cast<std::size_t>(r.startLocations[0]) : 0;
    std::vector<std::size_t>& before = out.before;
    std::vector<std::size_t>& after = out.after;
    before.resize(q.size());
    after.resize(q.size());
    std::vector<std::array<std::size_t, 3>> runs{{0, out.t0, 0}};   // query start, target start, length
    std::size_t qi = 0, ti = out.t0, rq = 0, rt = 0, len = 0;
    auto close = [&] {
        if (len >= kSolid) runs.push_back({rq, rt, len});
        len = 0;
    };
    for (int x = 0; x < r.alignmentLength; ++x) {
        const unsigned char op = r.alignment[x];
        if (op == EDLIB_EDOP_MATCH) {
            if (len == 0) { rq = qi; rt = ti; }
            out.aligned[qi] = out.matched[qi] = 1;
            before[qi] = ti; after[qi] = ti + 1;
            ++len; ++qi; ++ti;
            continue;
        }
        close();
        if (op == EDLIB_EDOP_MISMATCH) { out.aligned[qi] = 1; before[qi] = ti; after[qi] = ti + 1; }
        if (op == EDLIB_EDOP_INSERT) before[qi] = after[qi] = ti;   // a base t lacks
        if (op != EDLIB_EDOP_DELETE) ++qi;   // DELETE: a base of t only; INSERT: of q only
        if (op != EDLIB_EDOP_INSERT) ++ti;
    }
    close();
    edlibFreeAlignResult(r);
    out.t1 = infix ? ti : t.size();
    if (!solid) return out;
    runs.push_back({q.size(), out.t1, 0});
    for (std::size_t k = 0; k + 1 < runs.size(); ++k) {
        const auto [q0, t0, l0] = runs[k];
        for (std::size_t i = 0; i < l0; ++i) { before[q0 + i] = t0 + i; after[q0 + i] = t0 + i + 1; }
        const std::size_t ql = q0 + l0, tl = t0 + l0, qr = runs[k + 1][0], tr = runs[k + 1][1];
        const std::size_t lq = qr - ql, lt = tr - tl, m = std::min(lq, lt);
        std::size_t left = 0, right = 0;
        for (std::size_t i = 0; i < m; ++i) {
            left += q[ql + i] == t[tl + i];
            right += q[qr - m + i] == t[tr - m + i];
        }
        const bool at_left = left >= right;
        for (std::size_t i = 0; i < lq; ++i) {
            const std::size_t x = ql + i;
            if (at_left ? i < m : i >= lq - m) {
                const std::size_t y = at_left ? tl + i : tr - (lq - i);
                before[x] = y;
                after[x] = y + 1;
            } else {
                before[x] = after[x] = at_left ? tr : tl;   // a base t lacks
            }
        }
    }
    return out;
}

// A walk's units of one tandem stretch, kept for the haplotype's other transcripts with that stretch: each
// unit's bounds, the reference unit it resembles most and, once asked for, the alignment lifting its bases.
struct WalkUnit {
    std::size_t a = 0, b = 0, g = 0;
    bool lifted = false;
    std::vector<std::size_t> before, after;
};
using UnitCache = std::map<std::tuple<const char*, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t>, std::vector<WalkUnit>>;

std::vector<WalkUnit> walk_units(const std::string& ws, const Tandem& tm, std::size_t lo, std::size_t hi);

// The exons of a tandem stretch on walk W (its sequence `ws`) between bases lo and hi.
Placed tandem_exons(const std::string& ws, const Tandem& tm, std::size_t lo, std::size_t hi, UnitCache& cache) {
    Placed out;
    if (hi <= lo) return out;
    auto& units = cache.try_emplace({ws.data(), lo, hi, tm.rl, tm.rh, tm.first}).first->second;
    if (units.empty()) units = walk_units(ws, tm, lo, hi);
    out.units = units.size();
    for (WalkUnit& w : units) {
        const Tandem::Unit& un = tm.units[w.g];
        if (un.exons.empty()) continue;
        if (!w.lifted) {
            Lift l = lift(un.seq, ws.substr(w.a, w.b - w.a));
            w.before = std::move(l.before);
            w.after = std::move(l.after);
            w.lifted = true;
        }
        for (const auto& [a, b, ref] : un.exons) {
            const std::size_t x = w.before[a], y = w.after[b - 1];
            if (y <= x) { ++out.unplaced; continue; }
            out.resized += y - x != b - a;
            out.pieces.push_back({w.a + x, w.a + y, ref});
        }
    }
    return out;
}

// The units of a tandem stretch on walk W between bases lo and hi, each with its reference unit.
std::vector<WalkUnit> walk_units(const std::string& ws, const Tandem& tm, std::size_t lo, std::size_t hi) {
    std::vector<WalkUnit> out;
    std::size_t p = lo, q = hi;
    if (!tm.head.empty()) {   // the head, from the window's start
        const std::size_t e = prefix_end(tm.head, ws.substr(lo, std::min(hi - lo, tm.head.size() * 3 / 2 + 1000)));
        if (e != SIZE_MAX) p = lo + e;
    }
    if (!tm.tail.empty() && q > p) {   // the tail, back from the window's end
        const std::size_t span = std::min(q - p, tm.tail.size() * 3 / 2 + 1000);
        std::string rt(tm.tail.rbegin(), tm.tail.rend()), rw(ws.begin() + static_cast<long>(q - span), ws.begin() + static_cast<long>(q));
        std::reverse(rw.begin(), rw.end());
        const std::size_t e = prefix_end(rt, rw);
        if (e != SIZE_MAX) q -= e;
    }
    const std::size_t k = tm.units.size();
    // the reference unit sharing most 15-mers with a walk unit; ties go to the one nearest `want`
    auto best_unit = [&](const std::vector<std::uint64_t>& km, double want) {
        std::size_t best = 0;
        double share = -1.0;
        for (std::size_t g = 0; g < k; ++g) {
            const double j = jaccard(km, tm.units[g].kmers);
            const double d = std::fabs(static_cast<double>(g) - want), bd = std::fabs(static_cast<double>(best) - want);
            if (j > share || (j == share && d < bd)) { share = j; best = g; }
        }
        return std::make_pair(best, share);
    };
    // the units: each next one starts where the walk best holds a reference unit's first bases, or ends where
    // it best holds a unit's last bases, between half the shortest unit and the longest past the last start;
    // the nearer match to the reference decides, so an insertion or deletion near one end of a unit leaves
    // the boundary found from the other
    std::vector<std::pair<std::size_t, std::size_t>> seg;
    std::size_t max_len = 0;
    for (const auto& un : tm.units) max_len = std::max(max_len, un.seq.size());
    while (p < q) {
        if (!seg.empty() && q - p < tm.min_len / 2) {   // too little left for a unit: the last one's end
            seg.back().second = q;
            break;
        }
        const std::size_t from = p + tm.min_len / 2;
        std::size_t next = SIZE_MAX;
        int best = -1;
        if (from < q) {
            const std::string window = ws.substr(from, std::min(q - from, max_len + tm.min_len / 2));
            for (const auto& [probe, at_start] : tm.probes) {
                EdlibAlignResult r = align_edlib(probe, window, EDLIB_MODE_HW, false, static_cast<int>(probe.size() / 10));
                if (r.editDistance >= 0 && r.numLocations > 0 && (best < 0 || r.editDistance < best)) {
                    best = r.editDistance;
                    next = from + static_cast<std::size_t>(at_start ? r.startLocations[0] : r.endLocations[0] + 1);
                }
                edlibFreeAlignResult(r);
                if (best == 0) break;
            }
        }
        if (next == SIZE_MAX || next >= q) {
            // the rest is the last unit when it resembles one: past the stretch's end the window may run on
            if (seg.empty() || best_unit(unit_kmers(ws.substr(p, q - p)), 0).second >= 0.5) seg.emplace_back(p, q);
            break;
        }
        seg.emplace_back(p, next);
        p = next;
    }
    // each unit with the reference unit it resembles most (ties: the one nearest in relative position)
    const std::size_t n = seg.size();
    for (std::size_t j = 0; j < n; ++j) {
        const double want = n == k ? static_cast<double>(j) : std::floor(static_cast<double>(j * k) / static_cast<double>(n));
        WalkUnit w;
        w.a = seg[j].first;
        w.b = seg[j].second;
        w.g = best_unit(unit_kmers(ws.substr(w.a, w.b - w.a)), want).first;
        out.push_back(std::move(w));
    }
    return out;
}

// A copy placed by sequence: the reference stretch [q0, q0 + n) aligned inside a window of the oriented walk
// starting at base `a`.
struct SeqCopy {
    bool ok = false;
    bool reverse = false;
    bool marker_span = false;            // aligned over the marker's span, not the gene's
    std::size_t q0 = 0, a = 0;
    Lift map;
    std::size_t h_lo = 0, h_hi = 0;      // the aligned stretch, on the walk as written
    double coverage = 0.0, identity = 0.0;   // of the marker's exonic bases: aligned, and equal among those
};

// One copy of a gene on a walk: a pass (by nodes) or a sequence placement.
struct Copy {
    const Pass* pass = nullptr;
    const SeqCopy* seq = nullptr;
    bool reverse() const { return pass ? pass->reverse : seq->reverse; }
};

// Carry reference cut c to the walk through a sequence placement: exact where the base at the boundary is set
// against a base of the walk. A base the walk lacks sits, with its run of lacking bases, at one position of
// the walk: the boundary moves to the nearer end of the run, and a tie keeps the variant inside the exon.
CutOnWalk map_cut(const SeqCopy& sc, std::size_t c, bool start) {
    CutOnWalk out;
    const std::size_t n = sc.map.before.size();
    if (c < sc.q0 || c > sc.q0 + n) return out;
    const std::size_t i = c - sc.q0;
    const auto exact = [&](std::size_t k) { return sc.map.after[k] == sc.map.before[k] + 1; };
    if (start ? i < n && exact(i) : i > 0 && exact(i - 1)) {
        out.ok = true;
        out.pos = sc.a + (start ? sc.map.before[i] : sc.map.after[i - 1]);
        return out;
    }
    if (start ? i >= n : i == 0) return out;   // the cut lies at the stretch's end, outside it
    std::size_t ql = start ? i : i - 1, qr = ql;
    while (ql > 0 && !exact(ql - 1)) --ql;
    while (qr + 1 < n && !exact(qr + 1)) ++qr;
    ++qr;   // [ql, qr) lacking
    if (ql == 0 && qr == n) return out;
    out.ok = true;
    const long dl = static_cast<long>(c - (sc.q0 + ql)), dr = static_cast<long>(sc.q0 + qr - c);
    const bool left = ql == 0 ? false : qr == n ? true : dl != dr ? dl < dr : start;
    out.pos = sc.a + (left ? sc.map.after[ql - 1] : sc.map.before[qr]);
    const long d = left ? dl : dr;
    out.shift = left == start ? d : -d;
    return out;
}

CutOnWalk map_cut(const Walk& R, const Gene& g, const Copy& cp, const Walk& W, std::size_t c, bool start) {
    return cp.pass ? map_cut(R, g, *cp.pass, W, c, start) : map_cut(*cp.seq, c, start);
}

bool exon_shared(const Walk& R, const Gene& g, const Copy& cp, std::size_t cs, std::size_t ce) {
    if (cp.pass) return exon_shared(R, g, *cp.pass, cs, ce);
    const SeqCopy& sc = *cp.seq;
    const std::size_t n = sc.map.before.size();
    for (std::size_t x = std::max(cs, sc.q0); x < ce && x < sc.q0 + n; ++x)
        if (sc.map.after[x - sc.q0] == sc.map.before[x - sc.q0] + 1) return true;
    return false;
}

// The reference stretch [q_lo, q_hi) of a gene aligned inside the window of the walk that pass p, its best
// chain, points to: from the chain's anchors nearest the stretch's ends, with a margin for indels. Coverage
// and identity are measured on the marker's exonic bases.
SeqCopy align_stretch(const Walk& R, const std::string& rseq, const Gene& g, const Cuts& marker, const Pass& p, const Walk& W,
                      const std::string& wseq, std::size_t q_lo, std::size_t q_hi, bool any_distance) {
    SeqCopy sc;
    sc.reverse = p.reverse;
    const std::size_t len = q_hi - q_lo;
    // the walk positions of the stretch's two ends, each carried from the matched step nearest to it
    long a = 0, b = 0, da = -1, db = -1;
    for (std::size_t k = p.first; k <= p.last; ++k) {
        if (p.match[k] < 0) continue;
        const std::size_t i = g.gs + k, j = static_cast<std::size_t>(p.match[k]);
        const long x0 = static_cast<long>(R.pref[i]), x1 = static_cast<long>(R.pref[i + 1]);
        const long d0 = std::labs(x0 - static_cast<long>(q_lo)), d1 = std::labs(x1 - static_cast<long>(q_hi));
        if (da < 0 || d0 < da) { da = d0; a = static_cast<long>(W.pref[j]) + static_cast<long>(q_lo) - x0; }
        if (db < 0 || d1 < db) { db = d1; b = static_cast<long>(W.pref[j + 1]) + static_cast<long>(q_hi) - x1; }
    }
    const long margin = static_cast<long>(kWindowMargin * static_cast<double>(len)) + 100;
    a = std::max(0L, a - margin);
    b = std::min(static_cast<long>(W.len()), b + margin);
    if (da < 0 || b <= a || static_cast<double>(b - a) > kWindowMax * static_cast<double>(len) + 2 * static_cast<double>(margin)) return sc;
    sc.q0 = q_lo;
    sc.a = static_cast<std::size_t>(a);
    sc.map = lift(rseq.substr(q_lo, len), wseq.substr(sc.a, static_cast<std::size_t>(b - a)), true, false, any_distance);
    if (!sc.map.ok) return sc;
    {   // bases aligned where the alignment holds (kLocalIdentity); the others are bases the walk lacks
        Lift& l = sc.map;
        std::vector<std::size_t> match(len + 1, 0), count(len + 1, 0);
        for (std::size_t i = 0; i < len;) {
            std::size_t e = i;
            while (e < len && !l.aligned[e]) ++e;
            const bool long_deletion = e - i >= kLongDeletion;
            for (std::size_t x = i; x < std::max(e, i + 1); ++x) {
                match[x + 1] = match[x] + l.matched[x];
                count[x + 1] = count[x] + !(x < e && long_deletion);
            }
            i = std::max(e, i + 1);
        }
        std::vector<char> keep(len, 0);
        for (std::size_t i = 0; i < len; ++i) {
            const std::size_t lo = i > kLocalFlank ? i - kLocalFlank : 0, hi = std::min(len, i + kLocalFlank + 1);
            keep[i] = l.aligned[i] && static_cast<double>(match[hi] - match[lo]) >= kLocalIdentity * static_cast<double>(count[hi] - count[lo]);
        }
        for (std::size_t i = 0; i < len;) {   // runs of aligned bases without a gap, each between two gaps
            if (!l.aligned[i]) { ++i; continue; }
            std::size_t e = i + 1;
            while (e < len && l.aligned[e] && l.before[e] == l.after[e - 1]) ++e;
            if (i > 0 && e < len && e - i < kBetweenGaps)
                for (std::size_t x = i; x < e; ++x) keep[x] = 0;
            i = e;
        }
        for (std::size_t i = 0; i < len; ++i)
            if (!keep[i]) { l.aligned[i] = l.matched[i] = 0; l.after[i] = l.before[i]; }
    }
    std::size_t total = 0, aligned = 0, matched = 0;
    for (const auto& [s, e] : marker)
        for (std::size_t x = std::max(s, q_lo); x < std::min(e, q_hi); ++x) {
            ++total;
            aligned += sc.map.aligned[x - q_lo];
            matched += sc.map.matched[x - q_lo];
        }
    if (total == 0 || aligned == 0) return sc;
    sc.coverage = static_cast<double>(aligned) / static_cast<double>(total);
    sc.identity = static_cast<double>(matched) / static_cast<double>(aligned);
    // the stretch of the walk the copy takes: where the marker's span lies, as for a pass
    const std::size_t s0 = std::max(g.m_lo, q_lo) - q_lo, s1 = std::min(g.m_hi, q_hi) - q_lo;
    std::size_t lo = sc.a + sc.map.t0, hi = sc.a + sc.map.t1;
    if (s1 > s0 && sc.map.after[s1 - 1] > sc.map.before[s0]) { lo = sc.a + sc.map.before[s0]; hi = sc.a + sc.map.after[s1 - 1]; }
    sc.h_lo = p.reverse ? W.len() - hi : lo;
    sc.h_hi = p.reverse ? W.len() - lo : hi;
    sc.ok = true;
    return sc;
}

// A gene's copy found by sequence where its nodes fall short of the share. The gene's span is aligned first,
// within a band: a copy lies within it. Failing that, the marker's span is aligned at any distance, which
// is enough to tell a partial copy, or a copy of a short gene with a long isoform beside it; exons outside
// the marker's span are then not aligned.
SeqCopy place_by_sequence(const Walk& R, const std::string& rseq, const Gene& g, const Cuts& marker, const Pass& p,
                          const Walk& W, const std::string& wseq) {
    const auto cells = [](std::size_t len) {
        return static_cast<double>(len) * (static_cast<double>(len) * (1 + 2 * kWindowMargin) + 200);
    };
    SeqCopy sc;
    if (g.hi > g.lo && cells(g.hi - g.lo) <= kSequenceMaxCells) sc = align_stretch(R, rseq, g, marker, p, W, wseq, g.lo, g.hi, false);
    if (!sc.ok && g.m_hi > g.m_lo && cells(g.m_hi - g.m_lo) <= kSequenceMaxCellsAny) {
        sc = align_stretch(R, rseq, g, marker, p, W, wseq, g.m_lo, g.m_hi, true);
        sc.marker_span = g.m_lo != g.lo || g.m_hi != g.hi;
    }
    return sc;
}

// One transcript through one copy.
HstRecord project_one(const Walk& R, const std::string& rseq, const Gene& g, const Copy& cp, const Walk& W,
                      const std::string& wseq, const TranscriptFeature& tx, const Cuts& cuts,
                      const TranscriptsResult::EdgeClip& clip, const Tandem& tm, UnitCache& cache) {
    HstRecord r;
    r.reverse = cp.reverse();
    const std::size_t n = cuts.size();
    const bool minus = tx.strand == '-';
    // Exons and introns are numbered in transcript order over the whole annotated transcript, 1-based; a
    // transcript whose tandem exons follow the haplotype is numbered over its own exons.
    const auto exon_no = [&](std::size_t k) { return minus ? clip.exons - clip.first_exon - k : clip.first_exon + k + 1; };
    if (clip.bp_5p) add_detail(r.detail, "clipped:5':" + std::to_string(clip.bp_5p));
    if (clip.bp_3p) add_detail(r.detail, "clipped:3':" + std::to_string(clip.bp_3p));
    const auto in_tandem = [&](std::size_t k) { return tm.any && k >= tm.lo && k <= tm.hi; };
    std::vector<std::pair<std::size_t, std::size_t>> hx(n);
    bool missing = false;
    for (std::size_t k = 0; k < n; ++k) {
        if (in_tandem(k)) continue;
        const auto [cs, ce] = cuts[k];
        const std::string e = "e" + std::to_string(exon_no(k));
        const bool before = cp.pass && cs < R.pref[g.gs + cp.pass->first] && cp.pass->at_walk_start;
        const bool after = cp.pass && ce > R.pref[g.gs + cp.pass->last + 1] && cp.pass->at_walk_end;
        if (!exon_shared(R, g, cp, cs, ce)) {
            missing = true;
            const bool unaligned = cp.seq && (cs < cp.seq->q0 || ce > cp.seq->q0 + cp.seq->map.before.size());
            add_detail(r.detail, e + (before || after ? ":walk_end" : unaligned ? ":not_aligned" : ":missing"));
            continue;
        }
        const CutOnWalk a = map_cut(R, g, cp, W, cs, true), b = map_cut(R, g, cp, W, ce, false);
        if (!a.ok || !b.ok || b.pos <= a.pos) {
            missing = true;
            add_detail(r.detail, e + ":missing");
            continue;
        }
        hx[k] = {a.pos, b.pos};
        const auto note = [&](const CutOnWalk& x, bool start) {
            const std::string side = e + (start == minus ? ".end:" : ".start:");
            if (x.shift != 0) { ++r.moved; add_detail(r.detail, side + signed_text(x.shift)); }
            if (x.in_variant) add_detail(r.detail, side + "in_variant");
        };
        note(a, true);
        note(b, false);
    }
    std::vector<Piece> pieces;
    for (std::size_t k = 0; k < n && !missing; ++k) {
        if (!in_tandem(k)) {
            pieces.push_back({hx[k].first, hx[k].second, k});
        } else if (k == tm.lo) {
            if (!cp.pass) {   // a tandem stretch is placed through a pass's nodes
                missing = true;
                add_detail(r.detail, "tandem:missing");
                continue;
            }
            const Pass& p = *cp.pass;
            // The window between the exons around the stretch; where the stretch starts or ends the transcript,
            // the anchoring node the pass shares, else the pass's end (a cut carried into the repeat could land
            // in any unit).
            const auto shared = [&](std::size_t i) { return i != SIZE_MAX && p.match[i - g.gs] >= 0 ? p.match[i - g.gs] : -1L; };
            const long la = shared(tm.lo_anchor), ha = shared(tm.hi_anchor);
            const std::size_t lo = tm.lo > 0 ? hx[tm.lo - 1].second
                                 : la >= 0 ? W.pref[static_cast<std::size_t>(la) + 1] : W.pref[static_cast<std::size_t>(p.match[p.first])];
            const std::size_t hi = tm.hi + 1 < n ? hx[tm.hi + 1].first
                                 : ha >= 0 ? W.pref[static_cast<std::size_t>(ha)] : W.pref[static_cast<std::size_t>(p.match[p.last]) + 1];
            const Placed placed = tandem_exons(wseq, tm, lo, hi, cache);
            const std::vector<Piece>& found = placed.pieces;
            if (found.empty()) {
                missing = true;
                add_detail(r.detail, "tandem:missing");
            }
            if (found.size() != tm.hi - tm.lo + 1)
                add_detail(r.detail, "tandem:" + std::to_string(tm.hi - tm.lo + 1) + ">" + std::to_string(found.size()));
            if (placed.resized) add_detail(r.detail, "tandem_resized:" + std::to_string(placed.resized));
            if (placed.unplaced) add_detail(r.detail, "tandem_unplaced:" + std::to_string(placed.unplaced));
            pieces.insert(pieces.end(), found.begin(), found.end());
        }
    }
    const std::size_t m = pieces.size();
    const auto piece_no = [&](std::size_t q) { return tm.any ? (minus ? m - q : q + 1) : exon_no(pieces[q].ref); };
    for (std::size_t q = 1; q < m && !missing; ++q)
        if (pieces[q].a < pieces[q - 1].b) {
            missing = true;
            add_detail(r.detail, "e" + std::to_string(piece_no(q)) + ":overlaps_neighbour");
        }
    if (missing) {
        r.status = HstStatus::MissingExon;
        r.moved = 0;
        return r;
    }
    // Splice-site dinucleotides at each intron, compared with the reference's. An intron the haplotype
    // shortens below its two splice sites is lost: its exons are joined.
    bool lost = false;
    for (std::size_t q = 0; q + 1 < m; ++q) {
        const std::size_t intron = tm.any ? (minus ? m - 1 - q : q + 1)
                                 : minus ? clip.exons - 1 - clip.first_exon - pieces[q].ref : clip.first_exon + pieces[q].ref + 1;
        const std::size_t rl = cuts[pieces[q].ref].second, rr = cuts[pieces[q + 1].ref].first;
        const std::size_t hl = pieces[q].b, hr = pieces[q + 1].a;
        // The reference's own introns after the left exon and before the right one.
        if (pieces[q].ref + 1 >= n || cuts[pieces[q].ref + 1].first < rl + 4 || pieces[q + 1].ref == 0 ||
            rr < cuts[pieces[q + 1].ref - 1].second + 4)
            continue;
        if (hr < hl + 4) {
            lost = true;
            add_detail(r.detail, "i" + std::to_string(intron) + ":intron_lost");
            continue;
        }
        const std::string ref_left = rseq.substr(rl, 2), ref_right = rseq.substr(rr - 2, 2);
        const std::string hap_left = wseq.substr(hl, 2), hap_right = wseq.substr(hr - 2, 2);
        const auto note = [&](const char* site, const std::string& a, const std::string& b) {
            ++r.splice_changed;
            add_detail(r.detail, "i" + std::to_string(intron) + "." + site + ":" + (minus ? reverse_complement(a) : a) +
                                 ">" + (minus ? reverse_complement(b) : b));
        };
        if (ref_left != hap_left) note(minus ? "acceptor" : "donor", ref_left, hap_left);
        if (ref_right != hap_right) note(minus ? "donor" : "acceptor", ref_right, hap_right);
    }
    for (const Piece& x : pieces) r.sequence += wseq.substr(x.a, x.b - x.a);
    if (minus) r.sequence = reverse_complement(r.sequence);
    r.exons = m;
    for (std::size_t q = 0; q < m; ++q) {
        const Piece& x = pieces[minus ? m - 1 - q : q];
        r.exon_starts.push_back(r.reverse ? W.len() - x.b : x.a);
        r.exon_lengths.push_back(x.b - x.a);
    }
    const std::size_t lo = pieces.front().a, hi = pieces.back().b;
    r.path_start = r.reverse ? W.len() - hi : lo;
    r.path_end = r.reverse ? W.len() - lo : hi;
    r.status = lost ? HstStatus::IntronLost : clip.bp_5p + clip.bp_3p ? HstStatus::ClippedAtEdge
             : r.moved ? HstStatus::BoundaryMoved
             : r.splice_changed ? HstStatus::SpliceSiteChanged : HstStatus::Ok;
    return r;
}

bool apart_on_reference(const Gene& a, const Gene& b) { return a.m_hi <= b.m_lo || b.m_hi <= a.m_lo; }

// Two passes compete when they are of one gene and overlap on the walk, or of genes apart on the reference
// (paralogs folded onto the same nodes) and overlap by more than kParalogOverlap of the longer one.
bool compete(const Pass& x, const Pass& y, const std::vector<Gene>& genes) {
    const std::size_t ov = overlap(x.h_lo, x.h_hi, y.h_lo, y.h_hi);
    if (x.gene == y.gene) return ov > 0;
    if (!apart_on_reference(genes[x.gene], genes[y.gene])) return false;
    return static_cast<double>(ov) > kParalogOverlap * static_cast<double>(std::max(x.h_hi - x.h_lo, y.h_hi - y.h_lo));
}

// Bases of the walk (as written) in [lo, hi) that pass p matches to its gene's reference stretch.
std::size_t matched_bases(const Pass& p, const Walk& H, std::size_t lo, std::size_t hi) {
    std::size_t n = 0;
    auto it = std::upper_bound(H.pref.begin(), H.pref.end(), lo);
    for (auto k = std::lower_bound(p.walk_steps.begin(), p.walk_steps.end(),
                                   static_cast<std::size_t>(it - H.pref.begin()) - 1);
         k != p.walk_steps.end() && H.pref[*k] < hi; ++k)
        n += overlap(H.pref[*k], H.pref[*k + 1], lo, hi);
    return n;
}

// Pass x lies inside y, a longer pass of a gene apart on the reference, and y's alignment matches that
// stretch of the walk at least as well as x's does: the stretch is the part of y's copy that resembles x's gene
// (a pseudogene fragment, or a gene nested in the other copy of a duplicated gene), not a copy of x.
bool nested_in(const Pass& x, const Pass& y, const std::vector<Gene>& genes, const Walk& H) {
    if (x.gene == y.gene || !apart_on_reference(genes[x.gene], genes[y.gene])) return false;
    const double ov = static_cast<double>(overlap(x.h_lo, x.h_hi, y.h_lo, y.h_hi));
    if (!(ov > kParalogOverlap * static_cast<double>(x.h_hi - x.h_lo) &&
          ov <= kParalogOverlap * static_cast<double>(y.h_hi - y.h_lo)))
        return false;
    // Compared where y's alignment reaches, as y aligns only its own gene's stretch of the reference.
    const std::size_t lo = std::max(x.h_lo, H.pref[y.walk_steps.front()]);
    const std::size_t hi = std::min(x.h_hi, H.pref[y.walk_steps.back() + 1]);
    return matched_bases(y, H, lo, hi) >= matched_bases(x, H, lo, hi);
}

// Passes accepted best first, by the share of marker exons they walk, each against the ones accepted
// before it. A first round decides the copies of the longer genes; passes nested in those (nested_in) are
// dropped, and the rest compete again.
std::vector<std::size_t> accept_passes(const std::vector<Pass>& passes, const std::vector<Gene>& genes, const Walk& H,
                                       std::vector<std::string>* lost_to, const TranscriptsResult& res) {
    const auto name = [&](std::size_t gi) { return res.transcripts[genes[gi].tx.front()].gene_name; };
    std::vector<std::size_t> order(passes.size());
    std::iota(order.begin(), order.end(), std::size_t{0});
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
        const Pass& x = passes[a];
        const Pass& y = passes[b];
        if (x.share != y.share) return x.share > y.share;
        if (x.marker_matched != y.marker_matched) return x.marker_matched > y.marker_matched;
        if (x.matched != y.matched) return x.matched > y.matched;   // the whole stretch decides a tie
        return std::tie(x.gene, x.h_lo) < std::tie(y.gene, y.h_lo);
    });
    const auto greedy = [&](const std::vector<char>& dropped, bool note) {
        std::vector<std::size_t> kept;
        for (const std::size_t a : order) {
            if (dropped[a]) continue;
            bool ok = true;
            for (const std::size_t b : kept)
                if (compete(passes[a], passes[b], genes)) {
                    ok = false;
                    if (note && passes[a].gene != passes[b].gene && (*lost_to)[passes[a].gene].empty())
                        (*lost_to)[passes[a].gene] = "taken_by:" + name(passes[b].gene);
                    break;
                }
            if (ok) kept.push_back(a);
        }
        return kept;
    };
    std::vector<char> dropped(passes.size(), 0);
    const std::vector<std::size_t> first = greedy(dropped, false);
    for (std::size_t a = 0; a < passes.size(); ++a)
        for (const std::size_t b : first)
            if (nested_in(passes[a], passes[b], genes, H)) {
                dropped[a] = 1;
                if ((*lost_to)[passes[a].gene].empty()) (*lost_to)[passes[a].gene] = "inside:" + name(passes[b].gene);
                break;
            }
    return greedy(dropped, true);
}

}  // namespace

TranscriptsResult project_transcripts(const Graph& graph, const TranscriptsOptions& opt, HstSink& sink) {
    TranscriptsResult res;
    const bool project = !opt.gtf_path.empty();
    if (project) res.reference = resolve_reference_path_name(graph, opt.reference_path, "transcripts");
    const PathRecord* ref_path = nullptr;
    for (const PathRecord& p : graph.paths) {
        res.haplotypes.push_back(p.name);
        if (p.name == res.reference) ref_path = &p;
    }
    NodeTable t;
    {
        std::vector<std::string> ids;
        ids.reserve(graph.nodes.size());
        for (const auto& kv : graph.nodes) ids.push_back(kv.first);
        std::sort(ids.begin(), ids.end());
        for (const std::string& id : ids) {
            t.id.emplace(id, static_cast<int>(t.seq.size()));
            t.seq.push_back(graph.nodes.at(id).sequence);
        }
    }
    // Without a reference annotation nothing is projected: no genes, and the walk loop below adds the
    // annotated transcripts only.
    Walk R;
    std::string rseq;
    std::size_t base = 0, last = 0;
    std::vector<TranscriptFeature> partial, txs;
    if (project) {
        R = to_walk(*ref_path, t);
        rseq = spell(R, t);
        const ParsedReferencePath meta = parse_reference_path_label(ref_path->name, R.len(), ref_path->type == 'W');
        if (!meta.has_interval)
            throw std::runtime_error("transcripts: reference path '" + ref_path->name +
                                     "' carries no contig:start-end, so GTF coordinates cannot be placed on it");
        base = meta.region_start_1based;
        last = base + R.len() - 1;
        txs = parse_gtf_transcripts(opt.gtf_path, meta.chrom, base, last, &partial, "transcripts");
    }
    std::vector<TranscriptsResult::EdgeClip> tx_clip;
    for (const TranscriptFeature& t : txs) tx_clip.push_back({t.exons.size(), 0, 0, 0});
    // A transcript reaching past the path keeps its exons on the path, cut at the path's ends.
    for (TranscriptFeature& t : partial) {
        const auto cut = [&](const std::vector<std::pair<std::size_t, std::size_t>>& iv, std::size_t* left,
                             std::size_t* right, std::size_t* first) {
            std::vector<std::pair<std::size_t, std::size_t>> kept;
            for (std::size_t k = 0; k < iv.size(); ++k) {
                const auto [a, b] = iv[k];
                if (a < base) *left += std::min(b + 1, base) - a;
                if (b > last) *right += b - std::max(a, last + 1) + 1;
                if (b < base || a > last) continue;
                if (kept.empty() && first) *first = k;
                kept.emplace_back(std::max(a, base), std::min(b, last));
            }
            return kept;
        };
        TranscriptsResult::EdgeClip c;
        c.exons = t.exons.size();
        std::size_t left = 0, right = 0, unused = 0;
        t.exons = cut(t.exons, &left, &right, &c.first_exon);
        t.cds = cut(t.cds, &unused, &unused, nullptr);
        if (t.exons.empty()) {
            res.outside.push_back(std::move(t));
            continue;
        }
        c.bp_5p = t.strand == '-' ? right : left;
        c.bp_3p = t.strand == '-' ? left : right;
        txs.push_back(std::move(t));
        tx_clip.push_back(c);
    }

    // Genes, in reference order; transcripts regrouped to follow them.
    std::map<std::string, std::vector<std::size_t>> by_gene;
    for (std::size_t k = 0; k < txs.size(); ++k) by_gene[txs[k].gene_id].push_back(k);
    std::vector<std::pair<std::size_t, std::string>> gene_order;
    for (const auto& [gid, ks] : by_gene) {
        std::size_t lo = SIZE_MAX;
        for (const std::size_t k : ks) lo = std::min(lo, txs[k].exons.front().first);
        gene_order.emplace_back(lo, gid);
    }
    std::sort(gene_order.begin(), gene_order.end());
    std::vector<Gene> genes;
    std::vector<Cuts> cuts;
    for (const auto& [first_exon, gid] : gene_order) {
        Gene g;
        for (const std::size_t k : by_gene[gid]) {
            g.tx.push_back(res.transcripts.size());
            res.transcripts.push_back(std::move(txs[k]));
            res.clips.push_back(tx_clip[k]);
            const TranscriptFeature& tr = res.transcripts.back();
            Cuts c;
            for (const auto& [s, e] : tr.exons) c.emplace_back(s - base, e - base + 1);
            cuts.push_back(std::move(c));
        }
        // Marker: the canonical transcript, else the one with the most exonic bases (first such by id).
        std::size_t best_bp = 0;
        bool best_canon = false, have = false;
        g.lo = SIZE_MAX;
        for (const std::size_t k : g.tx) {
            std::size_t bp = 0;
            for (const auto& [a, b] : cuts[k]) bp += b - a;
            const bool canon = res.transcripts[k].canonical;
            if (!have || (canon && !best_canon) || (canon == best_canon && bp > best_bp)) {
                g.marker = k; best_bp = bp; best_canon = canon; have = true;
            }
            g.lo = std::min(g.lo, cuts[k].front().first);
            g.hi = std::max(g.hi, cuts[k].back().second);
        }
        g.marker_bp = best_bp;
        g.m_lo = cuts[g.marker].front().first;
        g.m_hi = cuts[g.marker].back().second;
        const auto step_of = [&](std::size_t c) {
            return static_cast<std::size_t>(std::upper_bound(R.pref.begin(), R.pref.end(), c) - R.pref.begin()) - 1;
        };
        const std::size_t s0 = step_of(g.lo), s1 = step_of(g.hi - 1);
        g.gs = s0 > 0 ? s0 - 1 : 0;
        g.ge = std::min(R.size(), s1 + 2);
        std::unordered_map<int, int> count;
        for (std::size_t i = g.gs; i < g.ge; ++i) ++count[R.node[i]];
        g.uniq.resize(g.ge - g.gs);
        for (std::size_t i = g.gs; i < g.ge; ++i) g.uniq[i - g.gs] = count[R.node[i]] == 1;
        genes.push_back(std::move(g));
    }
    for (std::size_t k = 0; k < res.transcripts.size(); ++k) {
        std::string s;
        for (const auto& [a, b] : cuts[k]) s += rseq.substr(a, b - a);
        res.reference_sequences.push_back(res.transcripts[k].strand == '-' ? reverse_complement(s) : s);
    }
    std::vector<Tandem> tandem;
    std::vector<std::size_t> gene_of(cuts.size(), 0);
    for (std::size_t gi = 0; gi < genes.size(); ++gi) for (const std::size_t k : genes[gi].tx) gene_of[k] = gi;
    for (std::size_t k = 0; k < cuts.size(); ++k) tandem.push_back(find_tandem(R, rseq, cuts[k], genes[gene_of[k]]));
    for (Gene& g : genes)
        for (const std::size_t k : g.tx)
            if (tandem[k].any) g.tandem_spans.emplace_back(cuts[k][tandem[k].lo].first, cuts[k][tandem[k].hi].second);
    // A gene's copies are told apart by its marker's coding exons where it has them: paralogs are compared
    // on what they encode, as untranslated ends differ more freely. A tandem stretch of the marker is left
    // out, as haplotypes carry its unit a different number of times.
    std::vector<Cuts> marker_cuts;
    for (Gene& g : genes) {
        const TranscriptFeature& m = res.transcripts[g.marker];
        Cuts c = cuts[g.marker];
        if (!m.cds.empty()) {
            c.clear();
            for (const auto& [s, e] : m.cds) c.emplace_back(s - base, e - base + 1);
        }
        const Tandem& tm = tandem[g.marker];
        if (tm.any) {
            const std::size_t t0 = cuts[g.marker][tm.lo].first, t1 = cuts[g.marker][tm.hi].second;
            Cuts kept;
            for (const auto& [a, b] : c) {
                if (a < t0) kept.emplace_back(a, std::min(b, t0));
                if (b > t1) kept.emplace_back(std::max(a, t1), b);
            }
            if (!kept.empty()) c = std::move(kept);
        }
        g.marker_bp = 0;
        for (const auto& [a, b] : c) g.marker_bp += b - a;
        marker_cuts.push_back(std::move(c));
    }

    const auto project_haplotype = [&](std::size_t h) {
        std::vector<HstRecord> out;
        const Walk H = to_walk(graph.paths[h], t);
        const Walk Hr = reversed_walk(H, t);
        const std::string hseq = spell(H, t);
        const std::string hrseq = reverse_complement(hseq);
        const Occurrences occ(H, t.seq.size());
        std::vector<Pass> passes;
        std::vector<double> best_share(genes.size(), 0.0);
        std::vector<Pass> weak(genes.size());   // per gene, its best pass below the share: where to look by sequence
        UnitCache unit_cache;
        for (std::size_t gi = 0; gi < genes.size(); ++gi) {
            const Gene& g = genes[gi];
            for (const bool rev : {false, true}) {
                const Walk& W = rev ? Hr : H;
                // Chains are completed only when their anchored bases inside the marker's span are a tenth
                // of the best chain's or more: the others are stray matches of single nodes, not passes.
                // Bases outside the span are left out, as long flanks would hide a short gene's other copies.
                const std::vector<std::vector<Anchor>> chains = chain_anchors(R, g, W, occ, rev, H.size());
                std::vector<std::size_t> bp(chains.size(), 0);
                std::size_t best = 0;
                for (std::size_t k = 0; k < chains.size(); ++k) {
                    for (const Anchor& a : chains[k]) bp[k] += overlap(a.x, a.x + a.w, g.m_lo, g.m_hi);
                    best = std::max(best, bp[k]);
                }
                for (std::size_t k = 0; k < chains.size(); ++k) {
                    if (10 * bp[k] < best || (bp[k] == 0 && best > 0)) continue;
                    Pass p = make_pass(R, g, marker_cuts[gi], W, chains[k], t, gi, rev);
                    best_share[gi] = std::max(best_share[gi], p.share);
                    if (p.share >= opt.min_exon_share) passes.push_back(std::move(p));
                    else if (p.marker_matched > weak[gi].marker_matched) weak[gi] = std::move(p);
                }
            }
        }
        std::vector<std::string> lost_to(genes.size());
        const std::vector<std::size_t> kept = accept_passes(passes, genes, H, &lost_to, res);
        std::vector<std::vector<const Pass*>> gene_passes(genes.size());
        std::vector<std::string> rival(passes.size());   // the best paralog a kept pass beat, and its share
        for (const std::size_t a : kept) {
            gene_passes[passes[a].gene].push_back(&passes[a]);
            const Pass* best = nullptr;
            for (const Pass& y : passes)
                if (y.gene != passes[a].gene && compete(passes[a], y, genes) && (!best || y.share > best->share)) best = &y;
            if (best) {
                char buf[64];
                std::snprintf(buf, sizeof buf, ":%.4f>%.4f", passes[a].share, best->share);
                rival[a] = "rival:" + res.transcripts[genes[best->gene].tx.front()].gene_name + buf;
            }
        }
        // A gene with no pass to the share, and not lost to another gene, is looked for by sequence where its
        // best chain lies. Placements are accepted best first (whole copies before partial ones, then by the
        // bases they match); one taking the stretch of an accepted pass or placement of a gene apart on the
        // reference is that gene's. Genes with a tandem stretch are placed by nodes only.
        std::vector<SeqCopy> placed(genes.size());
        {
            std::vector<SeqCopy> found(genes.size());
            std::vector<std::size_t> order;
            for (std::size_t gi = 0; gi < genes.size(); ++gi) {
                const Pass& wp = weak[gi];
                if (!gene_passes[gi].empty() || !lost_to[gi].empty() || !genes[gi].tandem_spans.empty() || wp.marker_matched == 0) continue;
                found[gi] = place_by_sequence(R, rseq, genes[gi], marker_cuts[gi], wp, wp.reverse ? Hr : H, wp.reverse ? hrseq : hseq);
                if (found[gi].ok && found[gi].identity >= opt.min_identity && found[gi].coverage >= kPartialCoverage) order.push_back(gi);
            }
            const auto rank = [&](std::size_t gi) {
                return std::make_tuple(found[gi].coverage >= opt.min_exon_share, found[gi].coverage * found[gi].identity);
            };
            std::stable_sort(order.begin(), order.end(), [&](std::size_t x, std::size_t y) { return rank(x) > rank(y); });
            for (const std::size_t gi : order) {
                const SeqCopy& sc = found[gi];
                const auto takes = [&](std::size_t other, std::size_t lo, std::size_t hi) {
                    if (other == gi || !apart_on_reference(genes[gi], genes[other])) return false;
                    const std::size_t ov = overlap(sc.h_lo, sc.h_hi, lo, hi);
                    return static_cast<double>(ov) > kParalogOverlap * static_cast<double>(std::max(sc.h_hi - sc.h_lo, hi - lo));
                };
                bool taken = false;
                for (const std::size_t a : kept) taken = taken || takes(passes[a].gene, passes[a].h_lo, passes[a].h_hi);
                for (std::size_t other = 0; other < genes.size(); ++other)
                    taken = taken || (placed[other].ok && takes(other, placed[other].h_lo, placed[other].h_hi));
                if (!taken) placed[gi] = sc;
            }
        }
        const auto fmt = [](const char* name, double v) {
            char buf[48];
            std::snprintf(buf, sizeof buf, "%s:%.3f", name, v);
            return std::string(buf);
        };
        for (std::size_t gi = 0; gi < genes.size(); ++gi) {
            auto& gp = gene_passes[gi];
            std::sort(gp.begin(), gp.end(), [](const Pass* a, const Pass* b) { return a->h_lo < b->h_lo; });
            const SeqCopy& sc = placed[gi];
            const bool by_sequence = sc.ok && sc.coverage >= opt.min_exon_share;
            for (const std::size_t k : genes[gi].tx) {
                if (gp.empty() && !by_sequence) {
                    HstRecord r;
                    r.transcript = k;
                    r.haplotype = h;
                    r.status = sc.ok ? HstStatus::PartialCopy : HstStatus::NoCopy;
                    r.detail = fmt("best_share", best_share[gi]);
                    if (sc.ok) { add_detail(r.detail, fmt("coverage", sc.coverage)); add_detail(r.detail, fmt("identity", sc.identity)); }
                    if (!lost_to[gi].empty()) add_detail(r.detail, lost_to[gi]);
                    out.push_back(std::move(r));
                    continue;
                }
                if (by_sequence) {
                    Copy cp;
                    cp.seq = &sc;
                    HstRecord r = project_one(R, rseq, genes[gi], cp, sc.reverse ? Hr : H, sc.reverse ? hrseq : hseq,
                                              res.transcripts[k], cuts[k], res.clips[k], tandem[k], unit_cache);
                    r.transcript = k;
                    r.haplotype = h;
                    r.copy = 1;
                    r.copies = 1;
                    add_detail(r.detail, sc.marker_span ? "placed_by_sequence:marker_span" : "placed_by_sequence");
                    add_detail(r.detail, fmt("best_share", best_share[gi]));
                    add_detail(r.detail, fmt("identity", sc.identity));
                    out.push_back(std::move(r));
                    continue;
                }
                for (std::size_t c = 0; c < gp.size(); ++c) {
                    const Pass& p = *gp[c];
                    Copy cp;
                    cp.pass = &p;
                    HstRecord r = project_one(R, rseq, genes[gi], cp, p.reverse ? Hr : H, p.reverse ? hrseq : hseq,
                                              res.transcripts[k], cuts[k], res.clips[k], tandem[k], unit_cache);
                    r.transcript = k;
                    r.haplotype = h;
                    r.copy = static_cast<int>(c + 1);
                    r.copies = static_cast<int>(gp.size());
                    const std::string& rv = rival[static_cast<std::size_t>(&p - passes.data())];
                    if (!rv.empty()) add_detail(r.detail, rv);
                    out.push_back(std::move(r));
                }
            }
        }
        return out;
    };

    res.has_annotation.assign(graph.paths.size(), 0);
    std::vector<std::vector<std::size_t>> annotated_on(graph.paths.size());
    if (!opt.hap_annotations.empty()) {
        res.annotated = load_hap_annotations(graph, opt.hap_annotations, opt.contig_alias, "transcripts");
        std::set<std::string> haps;
        std::ifstream table(opt.hap_annotations);
        for (std::string line; std::getline(table, line);) haps.insert(line.substr(0, line.find('\t')));
        for (std::size_t p = 0; p < graph.paths.size(); ++p) {
            const std::string& n = graph.paths[p].name;
            const std::size_t a = n.find('#'), b = a == std::string::npos ? a : n.find('#', a + 1);
            res.has_annotation[p] = b != std::string::npos && haps.count(n.substr(0, b)) != 0;
        }
        for (std::size_t k = 0; k < res.annotated.size(); ++k) annotated_on[res.annotated[k].path].push_back(k);
    }
    // The haplotype's own transcripts, cut from its walk; one with a projected record's id, copy and sequence
    // is that record, of both sources.
    const auto add_annotated = [&](std::size_t h, std::vector<HstRecord>& recs) {
        if (annotated_on[h].empty()) return;
        const std::string hseq = spell(to_walk(graph.paths[h], t), t);
        std::map<std::tuple<std::string, int>, std::size_t> projected;
        for (std::size_t k = 0; k < recs.size(); ++k)
            if (!recs[k].sequence.empty()) projected[{res.transcripts[recs[k].transcript].transcript_id, recs[k].copy}] = k;
        // The annotation's copies of each gene, found by gene id (a source gene), else by name.
        std::map<std::string, int> by_id, by_name;
        for (const std::size_t a : annotated_on[h]) {
            by_id[res.annotated[a].t.gene_id] = res.annotated[a].copies;
            by_name[res.annotated[a].t.gene_name] = res.annotated[a].copies;
        }
        for (HstRecord& r : recs) {
            const TranscriptFeature& t = res.transcripts[r.transcript];
            const auto i = by_id.find(t.gene_id), n = by_name.find(t.gene_name);
            r.annotated_copies = i != by_id.end() ? i->second : n != by_name.end() ? n->second : 0;
        }
        std::vector<HstRecord> extra;
        for (const std::size_t a : annotated_on[h]) {
            const AnnotatedTranscript& at = res.annotated[a];
            HstRecord r;
            r.source = kAnnotated;
            r.transcript = a;
            r.haplotype = h;
            r.copy = at.copy;
            r.copies = at.copies;
            r.annotated_copies = at.copies;
            r.record = at.record;
            for (const auto& [x, y] : at.t.exons) r.sequence += hseq.substr(x, y - x);
            if (at.t.strand == '-') r.sequence = reverse_complement(r.sequence);
            const auto it = projected.find({at.t.transcript_id, at.copy});
            if (it != projected.end() && recs[it->second].sequence == r.sequence) {
                recs[it->second].source |= kAnnotated;
                std::string& rec = recs[it->second].record;
                rec += (rec.empty() ? "" : ",") + at.record;
                add_detail(recs[it->second].detail, at.detail);
                continue;
            }
            const std::size_t n = at.t.exons.size();
            r.exons = n;
            for (std::size_t q = 0; q < n; ++q) {
                const auto [x, y] = at.t.exons[at.t.strand == '-' ? n - 1 - q : q];
                r.exon_starts.push_back(x);
                r.exon_lengths.push_back(y - x);
            }
            r.path_start = at.t.exons.front().first;
            r.path_end = at.t.exons.back().second;
            if (at.bp_5p) add_detail(r.detail, "clipped:5':" + std::to_string(at.bp_5p));
            if (at.bp_3p) add_detail(r.detail, "clipped:3':" + std::to_string(at.bp_3p));
            add_detail(r.detail, at.detail);
            r.status = at.bp_5p + at.bp_3p ? HstStatus::ClippedAtEdge : HstStatus::Ok;
            extra.push_back(std::move(r));
        }
        std::sort(extra.begin(), extra.end(), [&](const HstRecord& x, const HstRecord& y) {
            const TranscriptFeature& a = res.annotated[x.transcript].t;
            const TranscriptFeature& b = res.annotated[y.transcript].t;
            return std::tie(a.gene_name, a.transcript_id, x.copy, x.path_start) <
                   std::tie(b.gene_name, b.transcript_id, y.copy, y.path_start);
        });
        recs.insert(recs.end(), std::make_move_iterator(extra.begin()), std::make_move_iterator(extra.end()));
    };

    // Haplotypes in batches, so results reach the sink in path order without holding every sequence.
    const std::size_t threads = opt.threads ? opt.threads : std::max(1u, std::thread::hardware_concurrency());
    const std::size_t batch = threads * 4;
    for (std::size_t b0 = 0; b0 < graph.paths.size(); b0 += batch) {
        const std::size_t b1 = std::min(graph.paths.size(), b0 + batch);
        std::vector<std::vector<HstRecord>> got(b1 - b0);
        run_parallel(b1 - b0, threads, [&](std::size_t k) {
            if (project) got[k] = project_haplotype(b0 + k);
            add_annotated(b0 + k, got[k]);
        });
        // The reference walk must give back every transcript as annotated, in one of its copies.
        for (std::size_t k = 0; k < got.size(); ++k) {
            if (!project || &graph.paths[b0 + k] != ref_path) continue;
            std::vector<char> back(res.transcripts.size(), 0);
            for (const HstRecord& r : got[k])
                if ((r.source & kProjected) && r.sequence == res.reference_sequences[r.transcript]) back[r.transcript] = 1;
            for (std::size_t x = 0; x < back.size(); ++x)
                if (!back[x]) res.not_given_back.push_back(res.transcripts[x].transcript_id);
            for (HstRecord& r : got[k])
                if ((r.source & kProjected) && !back[r.transcript] && !r.sequence.empty()) r.status = HstStatus::NotAsAnnotated;
        }
        for (const auto& recs : got) sink.add(res, recs);
    }
    return res;
}

}  // namespace panvar
