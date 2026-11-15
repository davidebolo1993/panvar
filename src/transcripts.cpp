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
#include <unordered_set>
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
        case HstStatus::Discordant: return "discordant";
        case HstStatus::ParalogLift: return "paralog_lift";
        case HstStatus::MissingExon: return "missing_exon";
        case HstStatus::PartialCopy: return "partial_copy";
        case HstStatus::NoCopy: return "no_copy";
    }
    return "ok";
}

namespace {

// Chaining: how far back a link may reach, in anchors and in bases (or twice the gene's stretch; any
// length across a tandem stretch). Scores are integers, in 1/1280 of a base, so that equal chains tie
// exactly on every platform: a node scores its bases, and an indel costs 1/20 of a base per base and log2
// of its length.
constexpr std::size_t kChainLookback = 200;
constexpr std::size_t kChainMaxGap = 100000;
constexpr std::int64_t kScoreUnit = 1280;
constexpr std::int64_t kIndelPerBase = 64;
// The largest gap between two chained nodes filled by an exact node alignment (cells); larger ones are
// filled from both ends only.
constexpr std::size_t kFillMaxCells = 4000000;
// Passes of two genes taking the same stretch of a walk: overlapping by more than this share of the longer.
constexpr double kParalogOverlap = 0.5;
// A copy placed by sequence: the marker's exonic bases, aligned at the identity asked for, make a copy where
// they are covered to the share and a partial copy where covered to this much. The window reaches this share
// of the stretch past each of its ends, and is given up past this many times the stretch: the chain's
// anchors then disagree on where the gene lies.
constexpr double kPartialCoverage = 0.2;
constexpr double kWindowMargin = 0.1;
constexpr double kWindowMax = 3.0;
// An aligned base counts where the bases around it (this many on each side) match to this share of the
// matches, mismatches and gaps there, a gap of any length counting once: unit-cost alignment sets a stretch
// the walk lacks against whatever follows, base against base with many short gaps, while an indel between
// matching bases is one gap.
constexpr std::size_t kLocalFlank = 15;
constexpr double kLocalIdentity = 0.7;
// Gaps with fewer aligned bases than this between them are one indel, set again in one piece.
constexpr std::size_t kBetweenGaps = 12;
// The largest alignments (stretch times window, in cells) tried: of the gene's span, within a band of a
// quarter of its bases, and of the marker's span, at any distance.
constexpr double kSequenceMaxCells = 2.5e9;
constexpr double kSequenceMaxCellsAny = 5e8;
// Reference bases aligned with the stretch on each side of it: an indel across an end of the gene then lies
// inside the alignment, where it is a gap; at the alignment's end, the bases the walk lacks would be set
// against whatever follows.
constexpr std::size_t kEndFlank = 100;
// Bases of the reference walk on each side of a gene's span that its stretch takes in, at least one step.
constexpr std::size_t kGeneFlank = 100;

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
    std::size_t gs = 0, ge = 0;       // reference steps [gs, ge) holding the span and a flank on each side
    std::vector<char> uniq;           // per step of [gs, ge): its node occurs once there
    std::vector<char> own;            // per step of [gs, ge): the reference walks its node nowhere outside the span's steps
    std::size_t marker_bp = 0;
    std::vector<std::pair<std::size_t, std::size_t>> tandem_spans;   // its transcripts' tandem stretches
    Cuts exons;                       // every exon of its transcripts, merged, on the reference walk
    std::unordered_map<int, std::size_t> step_of;   // its once-occurring nodes: their reference step
};

// One pass of a walk through a gene: per reference step of [gs, ge), the walk step it matches or -1.
struct Pass {
    std::size_t gene = 0;
    bool reverse = false;
    std::vector<long> match;
    // where the walk's stretch is a second copy of a segment too (match_again): per reference step of that
    // copy, which `match` leaves out, its walk step or -1, and the copy's spans on the reference walk; empty
    // where there is none
    std::vector<long> again;
    std::vector<std::pair<std::size_t, std::size_t>> again_spans;
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

// What an indel of d bases costs a chain, in score units: per base, and log2(d), linear between powers of two.
std::int64_t indel_cost(std::size_t d) {
    if (d == 0) return 0;
    int msb = 0;
    while ((d >> (msb + 1)) != 0) ++msb;
    const std::uint64_t below = static_cast<std::uint64_t>(d) - (std::uint64_t{1} << msb);
    return static_cast<std::int64_t>(d) * kIndelPerBase + msb * kScoreUnit + static_cast<std::int64_t>((below * kScoreUnit) >> msb);
}

// Colinear chains of the gene's once-occurring reference nodes on walk W, best first; each anchor used once.
// Where a walk holds a stretch of the gene twice, moving from one copy to the next costs the same at every
// node, so equal chains differ only in where they move. A link with an indel inside an exon costs one unit
// more, so the move falls in an intron or outside the gene, never inside an exon; of equal links, the
// nearest on the walk wins.
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
    std::vector<std::int64_t> f(n);
    std::vector<long> prev(n, -1);
    // whether an exon of the gene starts before reference base b and ends after base a: a link from a node
    // ending at a to one starting at b moves inside it
    const auto in_exon = [&](std::size_t a, std::size_t b) {
        const auto it = std::partition_point(g.exons.begin(), g.exons.end(), [&](const auto& e) { return e.second <= a; });
        return it != g.exons.end() && it->first < b;
    };
    // A gene with a tandem stretch links across any length: haplotypes carry any number of units.
    const std::size_t max_gap = g.tandem_spans.empty() ? std::max(kChainMaxGap, 2 * (R.pref[g.ge] - R.pref[g.gs]))
                                                       : SIZE_MAX;
    for (std::size_t a = 0; a < n; ++a) {
        const std::int64_t w = static_cast<std::int64_t>(an[a].w) * kScoreUnit;
        f[a] = w;
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
            const std::size_t dd = across ? 0 : dx > dy ? dx - dy : dy - dx;
            const std::int64_t s = f[b] + w - indel_cost(dd) - (dd > 0 && in_exon(an[b].x + an[b].w, an[a].x) ? 1 : 0);
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
               const NodeTable& t, std::size_t gene_index, bool reverse, const std::vector<char>& taken) {
    Pass p;
    p.gene = gene_index;
    p.reverse = reverse;
    p.match.assign(g.ge - g.gs, -1);
    for (std::size_t k = 0; k < c.size(); ++k) {
        p.match[c[k].i - g.gs] = static_cast<long>(c[k].j);
        if (k + 1 < c.size()) fill_gap(R, c[k].i + 1, c[k + 1].i, W, c[k].j + 1, c[k + 1].j, t, p.match, g.gs, true);
    }
    // Past the ends, the walk is searched over a little more than the reference stretch left to match, up to
    // an anchor of a better chain (`taken`): that is another copy of the stretch. A gene with a tandem stretch
    // takes no such stop, as its units repeat on the walk and its pass must reach them all.
    const std::size_t lead = R.pref[c.front().i] - R.pref[g.gs];
    std::size_t j0 = c.front().j;
    while (j0 > 0 && W.pref[c.front().j] - W.pref[j0 - 1] <= lead + lead / 2 + 1000 && !taken[j0 - 1]) --j0;
    fill_gap(R, g.gs, c.front().i, W, j0, c.front().j, t, p.match, g.gs, true);
    const std::size_t trail = R.pref[g.ge] - R.pref[c.back().i + 1];
    std::size_t j1 = c.back().j + 1;
    while (j1 < W.size() && W.pref[j1 + 1] - W.pref[c.back().j + 1] <= trail + trail / 2 + 1000 && !taken[j1]) ++j1;
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

// Carry reference cut c (between bases) to the walk through a pass's match. `start`: c opens an exon.
CutOnWalk map_cut(const Walk& R, const Gene& g, const std::vector<long>& match, const Walk& W, std::size_t c, bool start) {
    CutOnWalk out;
    // The step holding base c, or the step after the span when c ends the reference walk.
    const std::size_t i = static_cast<std::size_t>(std::upper_bound(R.pref.begin(), R.pref.end(), c) - R.pref.begin()) - 1;
    const std::size_t off = c - R.pref[i];
    const auto m = [&](std::size_t k) -> long { return k >= g.gs && k < g.ge ? match[k - g.gs] : -1; };
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

bool exon_shared(const Walk& R, const Gene& g, const std::vector<long>& match, std::size_t cs, std::size_t ce) {
    for (std::size_t k = 0; k < match.size(); ++k) {
        if (match[k] < 0) continue;
        const std::size_t i = g.gs + k;
        if (R.pref[i] < ce && R.pref[i + 1] > cs) return true;
    }
    return false;
}

// Whether the matched steps over reference bases [cs, ce) lie further apart on the walk than on the reference
// by more than the stretch's length, with nodes of the gene between them: the match moves from one copy of a
// duplicated stretch to another inside it, and the bases between them are not the exon's.
bool exon_jumps(const Walk& R, const Gene& g, const std::vector<long>& match, const Walk& W, std::size_t cs, std::size_t ce) {
    const auto gene_between = [&](std::size_t j0, std::size_t j1) {
        for (std::size_t j = j0; j < j1; ++j)
            if (g.step_of.count(W.node[j])) return true;
        return false;
    };
    long pj = -1;
    std::size_t pi = 0;
    for (std::size_t k = 0; k < match.size(); ++k) {
        const std::size_t i = g.gs + k;
        if (match[k] < 0 || R.pref[i + 1] <= cs || R.pref[i] >= ce) continue;
        if (pj >= 0) {
            if (match[k] <= pj) return true;
            const std::size_t j0 = static_cast<std::size_t>(pj) + 1, j1 = static_cast<std::size_t>(match[k]);
            const std::size_t dy = W.pref[j1] - W.pref[j0], dx = R.pref[i] - R.pref[pi + 1];
            if (dy > dx + (ce - cs) && gene_between(j0, j1)) return true;
        }
        pj = match[k];
        pi = i;
    }
    return false;
}

// Bases of the reference stretch [cs, ce) on matched steps; with `own`, on those whose nodes are the gene's own.
std::size_t held_bases(const Walk& R, const Gene& g, const std::vector<long>& match, std::size_t cs, std::size_t ce, bool own) {
    std::size_t n = 0;
    for (std::size_t k = 0; k < match.size(); ++k)
        if (match[k] >= 0 && (!own || g.own[k])) n += overlap(R.pref[g.gs + k], R.pref[g.gs + k + 1], cs, ce);
    return n;
}

// A deletion from one copy of a segment into the next, where the graph folds the copies onto the same nodes,
// leaves the walk one copy with alleles of both. The pass follows one reference copy and leaves the other
// unmatched, though the walk passes nodes that only that copy has. A chain of such nodes (the gene's own, on
// reference steps the pass leaves out, inside the pass on the walk) is that other copy over the same
// stretch of the walk: the steps between its anchors are aligned node by node, so that an exon lying
// between them is found there too. Chains come best first; one over reference steps already taken is left.
void match_again(const Walk& R, const Gene& g, const Walk& W, const std::vector<std::vector<Anchor>>& chains,
                 const NodeTable& t, Pass& p) {
    const std::size_t w0 = static_cast<std::size_t>(p.match[p.first]), w1 = static_cast<std::size_t>(p.match[p.last]);
    for (const std::vector<Anchor>& c : chains) {
        std::vector<const Anchor*> left;
        for (const Anchor& a : c)
            if (p.match[a.i - g.gs] < 0 && a.j > w0 && a.j < w1) left.push_back(&a);
        if (left.size() < 2) continue;
        const std::size_t lo = R.pref[left.front()->i], hi = R.pref[left.back()->i + 1];
        bool taken = false;
        for (const auto& [s0, s1] : p.again_spans) taken = taken || (lo < s1 && s0 < hi);
        if (taken) continue;
        if (p.again.empty()) p.again.assign(p.match.size(), -1);
        for (std::size_t k = 0; k < left.size(); ++k) {
            p.again[left[k]->i - g.gs] = static_cast<long>(left[k]->j);
            if (k + 1 < left.size()) fill_gap(R, left[k]->i + 1, left[k + 1]->i, W, left[k]->j + 1, left[k + 1]->j, t, p.again, g.gs, true);
        }
        p.again_spans.emplace_back(lo, hi);
    }
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
    std::size_t last = 0;         // where its last unit ends
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
    t.last = last_end;
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
    // Units of one repeat are of comparable length, and each holds its exons whole. Equal starts that are not
    // (a node the reference walks twice for another reason) are no tandem stretch: the exons are left to
    // their nodes.
    bool units = lens.empty() || lens.back() <= 2 * lens.front();
    for (std::size_t u = 0; u < occ.size() && units; ++u) {
        const std::size_t a = occ[u], b = u + 1 < occ.size() ? occ[u + 1] : last_end;
        for (std::size_t k = t.lo; k <= t.hi; ++k)
            if (cuts[k].first >= a && cuts[k].first < b && cuts[k].second > b) units = false;
    }
    if (!units) return Tandem{};
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
using UnitCache = std::map<std::tuple<const char*, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t>,
                           std::vector<WalkUnit>>;

std::vector<WalkUnit> walk_units(const std::string& ws, const Tandem& tm, std::size_t lo, std::size_t hi);

// The exons of a tandem stretch on walk W (its sequence `ws`) between bases lo and hi.
Placed tandem_exons(const std::string& ws, const Tandem& tm, std::size_t lo, std::size_t hi, UnitCache& cache) {
    Placed out;
    if (hi <= lo) return out;
    auto& units = cache.try_emplace({ws.data(), lo, hi, tm.rl, tm.rh, tm.first, tm.last}).first->second;
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
    return cp.pass ? map_cut(R, g, cp.pass->match, W, c, start) : map_cut(*cp.seq, c, start);
}

bool exon_shared(const Walk& R, const Gene& g, const Copy& cp, std::size_t cs, std::size_t ce) {
    if (cp.pass) return exon_shared(R, g, cp.pass->match, cs, ce);
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
    {
        Lift& l = sc.map;
        const std::string q = rseq.substr(q_lo, len), t = wseq.substr(sc.a, static_cast<std::size_t>(b - a));
        // Gaps with fewer than kBetweenGaps aligned bases between them are one indel the alignment cut up: unit
        // costs let it scatter an indel's bases among others that happen to match. Their stretch is set again
        // as for a lift: the shorter side against the longer's left or right end, whichever matches more
        // bases, and the rest of the longer side one indel.
        const auto gap_at = [&](std::size_t i) { return !l.aligned[i] || (i > 0 && l.before[i] != l.after[i - 1]); };
        for (std::size_t i = 0; i < len;) {
            if (!gap_at(i)) { ++i; continue; }
            std::size_t g = i, gaps = 1, held = 0;   // the last gap of the cluster, its gaps, its aligned bases
            for (std::size_t k = i + 1; k < len && k - g <= kBetweenGaps; ++k)
                if (gap_at(k)) { gaps += !(k == g + 1 && !l.aligned[k] && !l.aligned[g]); g = k; }
            for (std::size_t k = i; k <= g; ++k) held += l.aligned[k];
            if (gaps > 1 && held > 0) {
                const std::size_t tl = i > 0 ? l.after[i - 1] : l.t0, tr = g + 1 < len ? l.before[g + 1] : l.t1;
                const std::size_t lq = g + 1 - i, lt = tr - tl, m = std::min(lq, lt);
                std::size_t left = 0, right = 0;
                for (std::size_t x = 0; x < m; ++x) {
                    left += q[i + x] == t[tl + x];
                    right += q[g + 1 - m + x] == t[tr - m + x];
                }
                const bool at_left = left >= right;
                for (std::size_t x = 0; x < lq; ++x) {
                    const bool set = at_left ? x < m : x >= lq - m;
                    const std::size_t y = !set ? (at_left ? tr : tl) : at_left ? tl + x : tr - (lq - x);
                    l.before[i + x] = y;
                    l.after[i + x] = y + set;
                    l.aligned[i + x] = set;
                    l.matched[i + x] = set && q[i + x] == t[y];
                }
            }
            i = g + 1;
        }
        // bases aligned where the alignment holds (kLocalIdentity); the others are bases the walk lacks
        std::vector<std::size_t> match(len + 1, 0), count(len + 1, 0);
        for (std::size_t i = 0; i < len; ++i) {
            // a gap opens at a base the walk lacks after one it has, and where the walk has bases between two
            const bool lacking = !l.aligned[i], opens = lacking && (i == 0 || l.aligned[i - 1]);
            const bool inserted = i > 0 && l.before[i] != l.after[i - 1];
            match[i + 1] = match[i] + l.matched[i];
            count[i + 1] = count[i] + (lacking ? 0 : 1) + opens + inserted;
        }
        for (std::size_t i = 0; i < len; ++i) {
            const std::size_t lo = i > kLocalFlank ? i - kLocalFlank : 0, hi = std::min(len, i + kLocalFlank + 1);
            if (l.aligned[i] && static_cast<double>(match[hi] - match[lo]) < kLocalIdentity * static_cast<double>(count[hi] - count[lo])) {
                l.aligned[i] = l.matched[i] = 0;
                l.after[i] = l.before[i];
            }
        }
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
// is enough to tell a partial copy, or a copy of a short gene with a long isoform beside it; exons reaching
// past the marker's span and its kEndFlank bases on each side are then not aligned.
SeqCopy place_by_sequence(const Walk& R, const std::string& rseq, const Gene& g, const Cuts& marker, const Pass& p,
                          const Walk& W, const std::string& wseq) {
    const auto cells = [](std::size_t len) {
        return static_cast<double>(len) * (static_cast<double>(len) * (1 + 2 * kWindowMargin) + 200);
    };
    const auto align = [&](std::size_t lo, std::size_t hi, bool any_distance) {
        return align_stretch(R, rseq, g, marker, p, W, wseq, lo > kEndFlank ? lo - kEndFlank : 0, std::min(rseq.size(), hi + kEndFlank),
                             any_distance);
    };
    SeqCopy sc;
    if (g.hi > g.lo && cells(g.hi - g.lo) <= kSequenceMaxCells) sc = align(g.lo, g.hi, false);
    if (!sc.ok && g.m_hi > g.m_lo && cells(g.m_hi - g.m_lo) <= kSequenceMaxCellsAny) {
        sc = align(g.m_lo, g.m_hi, true);
        sc.marker_span = g.m_lo != g.lo || g.m_hi != g.hi;
    }
    return sc;
}

// One transcript through one copy. With `own_share` (a pass below the share, which may lie on a paralog or
// on a diverged stretch), an exon counts only when that share of its bases is on nodes the reference walks
// nowhere outside the gene: the copy's test, exon by exon, tandem exons too.
HstRecord project_one(const Walk& R, const std::string& rseq, const Gene& g, const Copy& cp, const Walk& W,
                      const std::string& wseq, const TranscriptFeature& tx, const Cuts& cuts,
                      const TranscriptsResult::EdgeClip& clip, const Tandem& tm, UnitCache& cache, double own_share = 0.0) {
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
    std::vector<std::string> lacking(n), notes(n);   // per exon: why the copy lacks it, or its boundaries' notes
    // per exon the pass lacks: found where the walk's stretch is a second copy of its segment too (match_again)
    std::vector<char> again(n, 0);
    const auto in_second_copy = [&](std::size_t cs, std::size_t ce) {
        if (!cp.pass) return false;
        for (const auto& [s0, s1] : cp.pass->again_spans)
            if (s0 <= cs && ce <= s1) return true;
        return false;
    };
    const auto to_share = [&](std::size_t held, std::size_t cs, std::size_t ce) {
        return static_cast<double>(held) >= own_share * static_cast<double>(ce - cs);
    };
    for (std::size_t k = 0; k < n; ++k) {
        if (in_tandem(k)) continue;
        const auto [cs, ce] = cuts[k];
        const std::string e = "e" + std::to_string(exon_no(k));
        const bool before = cp.pass && cs < R.pref[g.gs + cp.pass->first] && cp.pass->at_walk_start;
        const bool after = cp.pass && ce > R.pref[g.gs + cp.pass->last + 1] && cp.pass->at_walk_end;
        CutOnWalk a, b;
        if (!exon_shared(R, g, cp, cs, ce)) {
            const bool unaligned = cp.seq && (cs < cp.seq->q0 || ce > cp.seq->q0 + cp.seq->map.before.size());
            lacking[k] = e + (before || after ? ":walk_end" : unaligned ? ":not_aligned" : ":missing");
        } else {
            a = map_cut(R, g, cp, W, cs, true);
            b = map_cut(R, g, cp, W, ce, false);
            if (!a.ok || !b.ok || b.pos <= a.pos) lacking[k] = e + ":missing";
            else if (cp.pass && exon_jumps(R, g, cp.pass->match, W, cs, ce)) lacking[k] = e + ":not_contiguous";
            else if (own_share > 0 && !to_share(held_bases(R, g, cp.pass->match, cs, ce, true), cs, ce))   // walked below the share, or on shared nodes
                lacking[k] = e + (to_share(held_bases(R, g, cp.pass->match, cs, ce, false), cs, ce) ? ":shared_nodes" : ":missing");
        }
        if (!lacking[k].empty() && in_second_copy(cs, ce) && exon_shared(R, g, cp.pass->again, cs, ce)) {
            a = map_cut(R, g, cp.pass->again, W, cs, true);
            b = map_cut(R, g, cp.pass->again, W, ce, false);
            again[k] = a.ok && b.ok && b.pos > a.pos && !exon_jumps(R, g, cp.pass->again, W, cs, ce) &&
                       (own_share <= 0 || to_share(held_bases(R, g, cp.pass->again, cs, ce, true), cs, ce));
        }
        if (!lacking[k].empty() && !again[k]) continue;
        hx[k] = {a.pos, b.pos};
        if (again[k]) notes[k] = e + ":merged_copy";
        const auto note = [&](const CutOnWalk& x, bool start) {
            const std::string side = e + (start == minus ? ".end:" : ".start:");
            if (x.shift != 0) { ++r.moved; add_detail(notes[k], side + signed_text(x.shift)); }
            if (x.in_variant) add_detail(notes[k], side + "in_variant");
        };
        note(a, true);
        note(b, false);
    }
    // An exon found on a second copy stays only in its place among the transcript's exons on the walk.
    for (std::size_t k = 0; k < n; ++k) {
        if (!again[k]) continue;
        const auto placed = [&](std::size_t q) { return !in_tandem(q) && (lacking[q].empty() || again[q]); };
        std::size_t lo = k, hi = k + 1;
        while (lo > 0 && !placed(lo - 1)) --lo;
        while (hi < n && !placed(hi)) ++hi;
        if ((lo > 0 && hx[lo - 1].second > hx[k].first) || (hi < n && hx[hi].first < hx[k].second)) again[k] = 0;
    }
    bool missing = false;
    for (std::size_t k = 0; k < n; ++k) {
        if (in_tandem(k)) continue;
        const bool lacks = !lacking[k].empty() && !again[k];
        missing = missing || lacks;
        add_detail(r.detail, lacks ? lacking[k] : notes[k]);
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
            // Below the share, each exon placed must lie on the gene's own nodes as the others must: the window
            // comes from a chain that may lie on a paralog.
            if (own_share > 0 && !missing) {
                std::set<int> mine;
                for (std::size_t q = 0; q < g.own.size(); ++q)
                    if (g.own[q]) mine.insert(R.node[g.gs + q]);
                for (const Piece& pc : found) {
                    std::size_t held = 0;
                    for (std::size_t j = step_at(W, pc.a); j < W.size() && W.pref[j] < pc.b; ++j)
                        if (mine.count(W.node[j])) held += overlap(W.pref[j], W.pref[j + 1], pc.a, pc.b);
                    if (!to_share(held, pc.a, pc.b)) {
                        missing = true;
                        add_detail(r.detail, "tandem:shared_nodes");
                        break;
                    }
                }
            }
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

// Pass x, over the stretch [x_lo, x_hi) of the walk, lies inside y, a longer pass of a gene apart on the
// reference, and y's alignment matches that stretch at least as well as x's does: the stretch is the part of
// y's copy that resembles x's gene (a pseudogene fragment, or a gene nested in the other copy of a duplicated
// gene), not a copy of x.
bool nested_in(const Pass& x, std::size_t x_lo, std::size_t x_hi, const Pass& y, const std::vector<Gene>& genes, const Walk& H) {
    if (x.gene == y.gene || !apart_on_reference(genes[x.gene], genes[y.gene])) return false;
    const double ov = static_cast<double>(overlap(x_lo, x_hi, y.h_lo, y.h_hi));
    if (!(ov > kParalogOverlap * static_cast<double>(x_hi - x_lo) &&
          ov <= kParalogOverlap * static_cast<double>(y.h_hi - y.h_lo)))
        return false;
    // Compared where y's alignment reaches, as y aligns only its own gene's stretch of the reference.
    const std::size_t lo = std::max(x_lo, H.pref[y.walk_steps.front()]);
    const std::size_t hi = std::min(x_hi, H.pref[y.walk_steps.back() + 1]);
    return matched_bases(y, H, lo, hi) >= matched_bases(x, H, lo, hi);
}

bool nested_in(const Pass& x, const Pass& y, const std::vector<Gene>& genes, const Walk& H) {
    return nested_in(x, x.h_lo, x.h_hi, y, genes, H);
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
        // then a fixed order: by place, then as the passes were found
        return std::tie(x.gene, x.h_lo, x.h_hi, x.reverse, a) < std::tie(y.gene, y.h_lo, y.h_hi, y.reverse, b);
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
    std::set<std::pair<std::string, std::string>> gtf_genes;   // every gene of the GTF, with annotations to judge
    if (project) {
        R = to_walk(*ref_path, t);
        rseq = spell(R, t);
        const ParsedReferencePath meta = parse_reference_path_label(ref_path->name, R.len(), ref_path->type == 'W');
        if (!meta.has_interval)
            throw std::runtime_error("transcripts: reference path '" + ref_path->name +
                                     "' carries no contig:start-end, so GTF coordinates cannot be placed on it");
        base = meta.region_start_1based;
        last = base + R.len() - 1;
        txs = parse_gtf_transcripts(opt.gtf_path, meta.chrom, base, last, &partial, "transcripts",
                                    opt.hap_annotations.empty() ? nullptr : &gtf_genes);
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
    std::vector<int> walked(t.seq.size(), 0);   // per node, the steps of the reference walk on it
    for (const int v : R.node) ++walked[static_cast<std::size_t>(v)];
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
        for (const std::size_t k : g.tx) g.exons.insert(g.exons.end(), cuts[k].begin(), cuts[k].end());
        std::sort(g.exons.begin(), g.exons.end());
        Cuts merged;
        for (const auto& e : g.exons) {
            if (!merged.empty() && e.first <= merged.back().second) merged.back().second = std::max(merged.back().second, e.second);
            else merged.push_back(e);
        }
        g.exons = std::move(merged);
        g.marker_bp = best_bp;
        g.m_lo = cuts[g.marker].front().first;
        g.m_hi = cuts[g.marker].back().second;
        const auto step_of = [&](std::size_t c) {
            return static_cast<std::size_t>(std::upper_bound(R.pref.begin(), R.pref.end(), c) - R.pref.begin()) - 1;
        };
        // The span with a flank of at least one step and kGeneFlank bases on each side: a boundary inside a
        // variant across the gene's first or last base needs the shared bases past the variant.
        const std::size_t s0 = step_of(g.lo), s1 = step_of(g.hi - 1);
        g.gs = s0;
        while (g.gs > 0 && (g.gs == s0 || R.pref[s0] - R.pref[g.gs] < kGeneFlank)) --g.gs;
        g.ge = s1 + 1;
        while (g.ge < R.size() && (g.ge == s1 + 1 || R.pref[g.ge] - R.pref[s1 + 1] < kGeneFlank)) ++g.ge;
        std::unordered_map<int, int> count, inside;   // steps per node over [gs, ge), and over the span alone
        for (std::size_t i = g.gs; i < g.ge; ++i) ++count[R.node[i]];
        for (std::size_t i = s0; i <= s1; ++i) ++inside[R.node[i]];
        g.uniq.resize(g.ge - g.gs);
        g.own.resize(g.ge - g.gs);
        for (std::size_t i = g.gs; i < g.ge; ++i) {
            g.uniq[i - g.gs] = count[R.node[i]] == 1;
            if (g.uniq[i - g.gs]) g.step_of.emplace(R.node[i], i);
            // the anchoring steps lie outside the gene: a neighbour's node there is not the gene's
            g.own[i - g.gs] = inside[R.node[i]] == walked[static_cast<std::size_t>(R.node[i])];
        }
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
        // per gene, its best pass below the share, where to look by sequence: by the marker's bases it matches,
        // then by all it matches, as a chain through introns and flanks alone still points to the gene
        std::vector<Pass> weak(genes.size());
        const auto reach = [](const Pass& p) { return std::make_pair(p.marker_matched, p.matched); };
        UnitCache unit_cache;
        std::vector<std::array<std::vector<std::vector<Anchor>>, 2>> chains_of(genes.size());   // per gene and orientation
        for (std::size_t gi = 0; gi < genes.size(); ++gi) {
            const Gene& g = genes[gi];
            for (const bool rev : {false, true}) {
                const Walk& W = rev ? Hr : H;
                // Chains are completed only when their anchored bases inside the marker's span are a tenth
                // of the best chain's or more: the others are stray matches of single nodes, not passes.
                // Bases outside the span are left out, as long flanks would hide a short gene's other copies.
                std::vector<std::vector<Anchor>>& chains = chains_of[gi][rev];
                chains = chain_anchors(R, g, W, occ, rev, H.size());
                std::vector<std::size_t> bp(chains.size(), 0);
                std::size_t best = 0;
                for (std::size_t k = 0; k < chains.size(); ++k) {
                    for (const Anchor& a : chains[k]) bp[k] += overlap(a.x, a.x + a.w, g.m_lo, g.m_hi);
                    best = std::max(best, bp[k]);
                }
                std::vector<char> taken(W.size(), 0);   // the walk steps of the chains before, for make_pass
                for (std::size_t k = 0; k < chains.size(); ++k) {
                    if (10 * bp[k] < best || (bp[k] == 0 && best > 0)) continue;
                    Pass p = make_pass(R, g, marker_cuts[gi], W, chains[k], t, gi, rev, taken);
                    if (g.tandem_spans.empty()) for (const Anchor& a : chains[k]) taken[a.j] = 1;
                    best_share[gi] = std::max(best_share[gi], p.share);
                    if (p.share >= opt.min_exon_share) passes.push_back(std::move(p));
                    else if (reach(p) > reach(weak[gi])) weak[gi] = std::move(p);
                }
            }
        }
        std::vector<std::string> lost_to(genes.size());
        const std::vector<std::size_t> kept = accept_passes(passes, genes, H, &lost_to, res);
        std::vector<std::vector<const Pass*>> gene_passes(genes.size());
        std::vector<std::string> rival(passes.size());   // the best paralog a kept pass beat, and its share
        for (const std::size_t a : kept) {
            match_again(R, genes[passes[a].gene], passes[a].reverse ? Hr : H, chains_of[passes[a].gene][passes[a].reverse], t, passes[a]);
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
        // reference, or nested in such a pass, is that gene's. Genes with a tandem stretch are placed by nodes
        // only.
        std::vector<SeqCopy> placed(genes.size());
        {
            std::vector<SeqCopy> found(genes.size());
            std::vector<std::size_t> order;
            for (std::size_t gi = 0; gi < genes.size(); ++gi) {
                const Pass& wp = weak[gi];
                if (!gene_passes[gi].empty() || !lost_to[gi].empty() || !genes[gi].tandem_spans.empty() || wp.matched == 0) continue;
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
                // taken by an accepted pass over the same stretch, or inside one that matches the stretch as well
                // as the chain the placement comes from: part of that gene's copy, as for a pass. The stretch is
                // the placement's: the chain may hold no node of the marker's span, and then has none.
                bool taken = false;
                for (const std::size_t a : kept)
                    taken = taken || takes(passes[a].gene, passes[a].h_lo, passes[a].h_hi) ||
                            nested_in(weak[gi], sc.h_lo, sc.h_hi, passes[a], genes, H);
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
            std::sort(gp.begin(), gp.end(), [](const Pass* a, const Pass* b) { return std::tie(a->h_lo, a->h_hi, a) < std::tie(b->h_lo, b->h_hi, b); });
            const SeqCopy& sc = placed[gi];
            const bool by_sequence = sc.ok && sc.coverage >= opt.min_exon_share;
            if (gp.empty() && !by_sequence) {
                // A gene without a copy, and not lost to another gene, may still hold whole transcripts where
                // its best chain lies: a deletion taking the marker's exons leaves those that do not use them.
                // Each transcript is projected through that chain; with one of them whole, the walk has a copy
                // below the share, whose other transcripts lack an exon.
                std::vector<HstRecord> whole;
                Pass& wp = weak[gi];
                if (lost_to[gi].empty() && wp.matched > 0) {
                    const Walk& W = wp.reverse ? Hr : H;
                    match_again(R, genes[gi], W, chains_of[gi][wp.reverse], t, wp);
                    Copy cp;
                    cp.pass = &wp;
                    bool any = false;
                    for (const std::size_t k : genes[gi].tx) {
                        whole.push_back(project_one(R, rseq, genes[gi], cp, W, wp.reverse ? hrseq : hseq, res.transcripts[k], cuts[k],
                                                    res.clips[k], tandem[k], unit_cache, opt.min_exon_share));
                        any = any || whole.back().status != HstStatus::MissingExon;
                    }
                    if (!any) whole.clear();
                }
                for (std::size_t x = 0; x < genes[gi].tx.size(); ++x) {
                    HstRecord r;
                    if (!whole.empty()) {
                        r = std::move(whole[x]);
                        r.copy = r.copies = 1;
                        r.copy_lo = weak[gi].h_lo;
                        r.copy_hi = weak[gi].h_hi;
                        add_detail(r.detail, "below_share");
                    } else {
                        r.status = sc.ok ? HstStatus::PartialCopy : HstStatus::NoCopy;
                        r.reverse = sc.ok && sc.reverse;
                    }
                    r.transcript = genes[gi].tx[x];
                    r.haplotype = h;
                    add_detail(r.detail, fmt("best_share", best_share[gi]));
                    if (sc.ok) { add_detail(r.detail, fmt("coverage", sc.coverage)); add_detail(r.detail, fmt("identity", sc.identity)); }
                    if (!lost_to[gi].empty()) add_detail(r.detail, lost_to[gi]);
                    out.push_back(std::move(r));
                }
                continue;
            }
            for (const std::size_t k : genes[gi].tx) {
                if (by_sequence) {
                    Copy cp;
                    cp.seq = &sc;
                    HstRecord r = project_one(R, rseq, genes[gi], cp, sc.reverse ? Hr : H, sc.reverse ? hrseq : hseq,
                                              res.transcripts[k], cuts[k], res.clips[k], tandem[k], unit_cache);
                    r.transcript = k;
                    r.haplotype = h;
                    r.copy = 1;
                    r.copies = 1;
                    r.copy_lo = sc.h_lo;
                    r.copy_hi = sc.h_hi;
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
                    r.copy_lo = p.h_lo;
                    r.copy_hi = p.h_hi;
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
    // Projected genes by id and by name, for the annotated transcripts' genes, and projected transcripts by id
    // (Ensembl ids without their version, as annotations and GTFs differ in it).
    std::unordered_map<std::string, std::size_t> gene_by_id, gene_by_name, tx_by_id;
    for (std::size_t gi = 0; gi < genes.size(); ++gi) {
        const TranscriptFeature& t0 = res.transcripts[genes[gi].tx.front()];
        gene_by_id.emplace(gene_id_key(t0.gene_id), gi);
        gene_by_name.emplace(t0.gene_name, gi);
    }
    for (std::size_t k = 0; k < res.transcripts.size(); ++k) tx_by_id.emplace(transcript_id_key(res.transcripts[k].transcript_id), k);
    // Genes the GTF has, on the path or not, by id and by name. Only these are judged as paralog lifts: the
    // place on the reference of a gene the GTF lacks is not known.
    std::unordered_set<std::string> gtf_ids, gtf_names;
    for (const auto& [id, name] : gtf_genes) {
        if (!id.empty()) gtf_ids.insert(gene_id_key(id));
        if (!name.empty()) gtf_names.insert(name);
    }
    // Each annotated transcript's gene: by its id where the GTF knows it, else by its name. `ann_gene` is the
    // projected gene (-1: none), `ann_known` whether the GTF has the gene at all.
    std::vector<long> ann_gene(res.annotated.size(), -1);
    std::vector<char> ann_known(res.annotated.size(), 0);
    for (std::size_t a = 0; a < res.annotated.size(); ++a) {
        AnnotatedTranscript& at = res.annotated[a];
        const std::string id = gene_id_key(at.t.gene_id);
        const auto gi = gene_by_id.find(id);
        if (gi != gene_by_id.end()) {
            ann_gene[a] = static_cast<long>(gi->second);
        } else if (gtf_ids.count(id) == 0) {
            const auto gn = gene_by_name.find(at.t.gene_name);
            if (gn != gene_by_name.end()) ann_gene[a] = static_cast<long>(gn->second);
        }
        ann_known[a] = ann_gene[a] >= 0 || gtf_ids.count(id) != 0 || gtf_names.count(at.t.gene_name) != 0;
        if (ann_gene[a] >= 0) at.projected_gene_id = res.transcripts[genes[static_cast<std::size_t>(ann_gene[a])].tx.front()].gene_id;
    }
    // A projected record the annotation may not overrule: its every boundary carried, or moved off a variant.
    const auto clean = [](HstStatus s) {
        return s == HstStatus::Ok || s == HstStatus::SpliceSiteChanged || s == HstStatus::BoundaryMoved || s == HstStatus::ClippedAtEdge;
    };
    // The haplotype's own transcripts, cut from its walk. Each annotated gene record is a copy, matched to the
    // projected copy of its gene it overlaps on the walk. One with a projected record's id, copy and sequence
    // is that record, of both sources; where they differ, a clean projected record wins. The annotation fills
    // a transcript or copy the projection fails on; a copy lying on another gene's copy is a paralog lift.
    const auto add_annotated = [&](std::size_t h, std::vector<HstRecord>& recs) {
        if (annotated_on[h].empty()) return;
        const std::string hseq = spell(to_walk(graph.paths[h], t), t);
        // The projected copies: the stretch each was found on, and that stretch widened to its written transcripts.
        struct ProjCopy {
            std::size_t gene = 0;
            int copy = 0;
            std::size_t s_lo = 0, s_hi = 0, lo = 0, hi = 0;
        };
        std::vector<ProjCopy> pcs;
        std::map<std::pair<std::size_t, int>, std::size_t> pc_of;
        std::map<std::pair<std::size_t, int>, std::size_t> rec_of;   // (transcript, copy) -> projected record
        std::vector<int> projected_copies(genes.size(), 0);
        for (std::size_t k = 0; k < recs.size(); ++k) {
            const HstRecord& r = recs[k];
            rec_of[{r.transcript, r.copy}] = k;
            if (r.copy == 0) continue;
            const std::size_t gi = gene_of[r.transcript];
            projected_copies[gi] = r.copies;
            const auto [it, fresh] = pc_of.emplace(std::make_pair(gi, r.copy), pcs.size());
            if (fresh) pcs.push_back({gi, r.copy, r.copy_lo, r.copy_hi, r.copy_lo, r.copy_hi});
            ProjCopy& pc = pcs[it->second];
            if (!r.sequence.empty()) {
                pc.lo = std::min(pc.lo, r.path_start);
                pc.hi = std::max(pc.hi, r.path_end);
            }
        }
        // The copy of another gene, apart from `gi` on the reference (any gene, for one off the reference path),
        // whose stretch holds over half of these exonic bases; -1 for none.
        const auto lifted_onto = [&](const Cuts& ex, long gi) -> long {
            std::size_t bp = 0, best_bp = 0;
            for (const auto& [a, b] : ex) bp += b - a;
            long best = -1;
            for (std::size_t c = 0; c < pcs.size(); ++c) {
                const ProjCopy& pc = pcs[c];
                if (gi >= 0 && (static_cast<std::size_t>(gi) == pc.gene || !apart_on_reference(genes[static_cast<std::size_t>(gi)], genes[pc.gene])))
                    continue;
                const std::size_t n = marker_overlap(ex, pc.s_lo, pc.s_hi);
                if (2 * n > bp && n > best_bp) { best = static_cast<long>(c); best_bp = n; }
            }
            return best;
        };
        // The annotation's copies: its gene records, with their gene among the projected ones (-1: none) and
        // their exonic bases, merged.
        struct AnnCopy {
            long gene = -1;
            bool placed = false;      // its gene is known on the reference: projected, or in the GTF off the path
            std::string key;          // the gene its copies are counted for, where not a projected one
            Cuts ex;
            std::vector<std::size_t> members;   // indices into res.annotated
            long match = -1, lift = -1;         // the projected copy of its gene it is, or of another gene it lies on
            bool kept = false;                  // holds a record that is not kept out
            int copy = 0;
        };
        std::vector<AnnCopy> acs;
        {
            std::map<std::pair<std::string, std::string>, std::size_t> of;
            for (const std::size_t a : annotated_on[h]) {
                const AnnotatedTranscript& at = res.annotated[a];
                const auto [it, fresh] = of.emplace(std::make_pair(at.gene_key, at.gene_record), acs.size());
                if (fresh) {
                    AnnCopy ac;
                    ac.gene = ann_gene[a];
                    ac.placed = ann_known[a] != 0;
                    ac.key = at.gene_key;
                    acs.push_back(std::move(ac));
                }
                AnnCopy& ac = acs[it->second];
                ac.members.push_back(a);
                ac.ex.insert(ac.ex.end(), at.t.exons.begin(), at.t.exons.end());
            }
        }
        for (AnnCopy& ac : acs) {
            std::sort(ac.ex.begin(), ac.ex.end());
            Cuts merged;
            for (const auto& iv : ac.ex) {
                if (!merged.empty() && iv.first <= merged.back().second) merged.back().second = std::max(merged.back().second, iv.second);
                else merged.push_back(iv);
            }
            ac.ex = std::move(merged);
            std::size_t best_bp = 0;
            for (std::size_t c = 0; ac.gene >= 0 && c < pcs.size(); ++c) {
                if (pcs[c].gene != static_cast<std::size_t>(ac.gene)) continue;
                const std::size_t n = marker_overlap(ac.ex, pcs[c].lo, pcs[c].hi);
                if (n > best_bp) { best_bp = n; ac.match = static_cast<long>(c); }
            }
            if (ac.match < 0 && ac.placed) ac.lift = lifted_onto(ac.ex, ac.gene);
        }

        // Each annotated transcript as a record; `copy_of` its annotated copy.
        std::vector<HstRecord> extra;
        std::vector<std::size_t> copy_of;
        std::set<std::size_t> joined;                 // projected records an annotated one is
        std::map<std::size_t, std::size_t> overruled; // projected records an annotated one differs from -> one such
        for (std::size_t ci = 0; ci < acs.size(); ++ci) {
            AnnCopy& ac = acs[ci];
            for (const std::size_t a : ac.members) {
                const AnnotatedTranscript& at = res.annotated[a];
                HstRecord r;
                r.source = kAnnotated;
                r.transcript = a;
                r.haplotype = h;
                r.annotated_copies = at.copies;
                r.record = at.record;
                for (const auto& [x, y] : at.t.exons) r.sequence += hseq.substr(x, y - x);
                if (at.t.strand == '-') r.sequence = reverse_complement(r.sequence);
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
                // its projected transcript, of the copy's gene
                long k = -1;
                if (ac.gene >= 0) {
                    const auto it = tx_by_id.find(transcript_id_key(at.t.transcript_id));
                    if (it != tx_by_id.end() && gene_of[it->second] == static_cast<std::size_t>(ac.gene)) k = static_cast<long>(it->second);
                }
                const auto lift = [&](long c) {
                    add_detail(r.detail, "on:" + res.transcripts[genes[pcs[static_cast<std::size_t>(c)].gene].tx.front()].gene_name);
                    if (!opt.prefer_annotation) r.status = HstStatus::ParalogLift;
                };
                // the projection fails here: the annotation fills it, unless it lies on another gene's copy or
                // has other exons than the transcript
                const auto fill = [&]() {
                    const long c = lifted_onto(at.t.exons, ac.gene);
                    if (c >= 0) return lift(c);
                    const std::size_t ref_exons = res.transcripts[static_cast<std::size_t>(k)].exons.size();
                    if (n != ref_exons && !opt.prefer_annotation) {
                        r.status = HstStatus::Discordant;
                        add_detail(r.detail, "exons:" + std::to_string(n) + "/" + std::to_string(ref_exons));
                        return;
                    }
                    add_detail(r.detail, "annotation_fill");
                };
                if (ac.lift >= 0) {
                    lift(ac.lift);
                } else if (ac.match >= 0 && k >= 0) {
                    const auto pit = rec_of.find({static_cast<std::size_t>(k), pcs[static_cast<std::size_t>(ac.match)].copy});
                    if (pit == rec_of.end()) {
                        fill();
                        ac.kept = ac.kept || !annotation_status(r.status);
                        extra.push_back(std::move(r));
                        copy_of.push_back(ci);
                        continue;
                    }
                    const std::size_t pi = pit->second;
                    HstRecord& pr = recs[pi];
                    if (!pr.sequence.empty() && pr.sequence == r.sequence) {
                        pr.source |= kAnnotated;
                        pr.record += (pr.record.empty() ? "" : ",") + at.record;
                        add_detail(pr.detail, at.detail);
                        joined.insert(pi);
                        continue;
                    }
                    if (pr.sequence.empty()) fill();
                    else if (opt.prefer_annotation) overruled.emplace(pi, extra.size());
                    else if (clean(pr.status)) r.status = HstStatus::Discordant;
                } else if (ac.match < 0 && k >= 0) {
                    fill();
                }
                ac.kept = ac.kept || !annotation_status(r.status);
                extra.push_back(std::move(r));
                copy_of.push_back(ci);
            }
        }
        // With the annotation preferred, a projected record it differs from, and none of its records joins,
        // is kept out.
        for (const auto& [pi, x] : overruled) {
            (void)x;
            if (joined.count(pi)) continue;
            recs[pi].status = HstStatus::Discordant;
            add_detail(recs[pi].detail, "annotated:" + extra[x].record);
        }
        // Copies: the projected ones, then the annotated ones the projection lacks, along the walk; a paralog
        // lift (unless the annotation is preferred) or a copy whose every record is kept out is none.
        std::vector<std::size_t> order(acs.size());
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](std::size_t x, std::size_t y) { return acs[x].ex.front() < acs[y].ex.front(); });
        std::vector<int> copies = projected_copies;
        std::map<std::string, int> other_copies;   // genes the projection lacks, by key
        for (const std::size_t ci : order) {
            AnnCopy& ac = acs[ci];
            if (ac.match >= 0) { ac.copy = pcs[static_cast<std::size_t>(ac.match)].copy; continue; }
            if (!ac.kept) continue;
            ac.copy = ac.gene >= 0 ? ++copies[static_cast<std::size_t>(ac.gene)] : ++other_copies[ac.key];
        }
        for (HstRecord& r : recs) r.copies = copies[gene_of[r.transcript]];
        for (std::size_t x = 0; x < extra.size(); ++x) {
            const AnnCopy& ac = acs[copy_of[x]];
            extra[x].copy = ac.copy;
            extra[x].copies = ac.gene >= 0 ? copies[static_cast<std::size_t>(ac.gene)] : other_copies[ac.key];
        }
        // The annotation's own copies of each projected gene.
        std::map<long, int> own_copies;
        for (const std::size_t a : annotated_on[h])
            if (ann_gene[a] >= 0) own_copies[ann_gene[a]] = res.annotated[a].copies;
        for (HstRecord& r : recs) {
            const auto c = own_copies.find(static_cast<long>(gene_of[r.transcript]));
            r.annotated_copies = c != own_copies.end() ? c->second : 0;
        }
        std::vector<std::size_t> idx(extra.size());
        std::iota(idx.begin(), idx.end(), 0);
        std::stable_sort(idx.begin(), idx.end(), [&](std::size_t x, std::size_t y) {
            const TranscriptFeature& a = res.annotated[extra[x].transcript].t;
            const TranscriptFeature& b = res.annotated[extra[y].transcript].t;
            return std::tie(a.gene_name, a.transcript_id, extra[x].copy, extra[x].path_start) <
                   std::tie(b.gene_name, b.transcript_id, extra[y].copy, extra[y].path_start);
        });
        for (const std::size_t x : idx) recs.push_back(std::move(extra[x]));
    };

    // Haplotypes in batches, so results reach the sink in path order without holding every sequence.
    const std::size_t threads = opt.threads ? opt.threads : std::max(1u, std::thread::hardware_concurrency());
    const std::size_t batch = threads * 4;
    for (std::size_t b0 = 0; b0 < graph.paths.size(); b0 += batch) {
        const std::size_t b1 = std::min(graph.paths.size(), b0 + batch);
        std::vector<std::vector<HstRecord>> got(b1 - b0);
        if (project) run_parallel(b1 - b0, threads, [&](std::size_t k) { got[k] = project_haplotype(b0 + k); });
        // The reference walk must give back every transcript as annotated, in one of its copies; this comes
        // before the haplotypes' own annotations, which judge projected records by their status.
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
        run_parallel(b1 - b0, threads, [&](std::size_t k) { add_annotated(b0 + k, got[k]); });
        for (const auto& recs : got) sink.add(res, recs);
    }
    return res;
}

}  // namespace panvar
