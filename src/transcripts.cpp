#include "panvar/transcripts.hpp"

#include "panvar/cli_utils.hpp"
#include "panvar/graph_utils.hpp"
#include "panvar/parallel.hpp"
#include "panvar/ref_path.hpp"

#include <algorithm>
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
        case HstStatus::MissingExon: return "missing_exon";
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
// from the first to the last such exon is cut into units, each starting where the stretch's first exon
// starts; a unit's layout is the kinds of exon it holds, in order. On a haplotype the same start cuts its
// walk into units, and each unit takes the layout of one reference unit, so the transcript follows the
// haplotype's own number of units and the reference walk gives back its annotation.
struct Tandem {
    bool any = false;
    std::size_t lo = 0, hi = 0;   // genomic exon indices of the stretch, inclusive
    struct Kind {
        int start_node = 0, end_node = 0;
        char start_rev = 0, end_rev = 0;
        std::size_t start_off = 0, end_off = 0;   // offsets of the exon's first and last base in their nodes
        std::size_t ref_exon = 0, len = 0;
    };
    std::vector<Kind> kinds;
    std::vector<std::vector<std::size_t>> units;   // per reference unit, its exons' kinds in order
    std::vector<std::vector<int>> unit_nodes;      // per reference unit, its nodes (sorted, distinct)
    std::vector<std::size_t> unit_len;             // per reference unit, its bases (the last: a whole unit's)
};

std::size_t step_at(const Walk& w, std::size_t c) {
    return static_cast<std::size_t>(std::upper_bound(w.pref.begin(), w.pref.end(), c) - w.pref.begin()) - 1;
}

Tandem find_tandem(const Walk& R, const Cuts& cuts) {
    Tandem t;
    std::vector<Tandem::Kind> kind(cuts.size());
    std::map<std::tuple<int, char, std::size_t>, std::size_t> starts;
    for (std::size_t k = 0; k < cuts.size(); ++k) {
        const std::size_t i = step_at(R, cuts[k].first), j = step_at(R, cuts[k].second - 1);
        kind[k] = {R.node[i], R.node[j], R.rev[i], R.rev[j], cuts[k].first - R.pref[i], cuts[k].second - 1 - R.pref[j],
                   k, cuts[k].second - cuts[k].first};
        ++starts[{R.node[i], R.rev[i], kind[k].start_off}];
    }
    for (std::size_t k = 0; k < cuts.size(); ++k) {
        if (starts[{kind[k].start_node, kind[k].start_rev, kind[k].start_off}] < 2) continue;
        if (!t.any) t.lo = k;
        t.hi = k;
        t.any = true;
    }
    if (!t.any) return t;
    const auto same_start = [&](const Tandem::Kind& x, const Tandem::Kind& y) {
        return std::tie(x.start_node, x.start_rev, x.start_off) == std::tie(y.start_node, y.start_rev, y.start_off);
    };
    std::vector<std::size_t> unit_start;
    for (std::size_t k = t.lo; k <= t.hi; ++k) {
        const Tandem::Kind& x = kind[k];
        std::size_t id = t.kinds.size();
        for (std::size_t q = 0; q < t.kinds.size(); ++q) {
            const Tandem::Kind& y = t.kinds[q];
            if (same_start(x, y) && std::tie(y.end_node, y.end_rev, y.end_off) == std::tie(x.end_node, x.end_rev, x.end_off))
                id = q;
        }
        if (id == t.kinds.size()) t.kinds.push_back(x);
        if (k == t.lo || same_start(x, kind[t.lo])) {
            t.units.emplace_back();
            unit_start.push_back(cuts[k].first);
        }
        t.units.back().push_back(id);
    }
    unit_start.push_back(cuts[t.hi].second);
    for (std::size_t u = 0; u + 1 < unit_start.size(); ++u) t.unit_len.push_back(unit_start[u + 1] - unit_start[u]);
    // The last unit ends at its last exon: it is compared over a whole unit, the median of the others.
    if (t.units.size() > 1) {
        std::vector<std::size_t> others(t.unit_len.begin(), t.unit_len.end() - 1);
        std::nth_element(others.begin(), others.begin() + others.size() / 2, others.end());
        t.unit_len.back() = std::max(t.unit_len.back(), others[others.size() / 2]);
    }
    for (std::size_t u = 0; u + 1 < unit_start.size(); ++u) {
        std::vector<int> nodes;
        const std::size_t end = std::min(unit_start[u] + t.unit_len[u], R.len());
        for (std::size_t i = step_at(R, unit_start[u]); i < R.size() && R.pref[i] < end; ++i) nodes.push_back(R.node[i]);
        std::sort(nodes.begin(), nodes.end());
        nodes.erase(std::unique(nodes.begin(), nodes.end()), nodes.end());
        t.unit_nodes.push_back(std::move(nodes));
    }
    return t;
}

struct Piece {
    std::size_t a = 0, b = 0;   // on the oriented walk
    std::size_t ref = 0;        // the reference exon it stands for
};

struct Placed {
    std::vector<Piece> pieces;
    std::size_t exact = 0, inferred = 0, resized = 0;
};

// One reference unit's layout placed on walk W between bases lo and hi, kind after kind: each exon where W
// holds its first and last nodes about its length apart. Where W lacks one of the two (a variant at the
// exon's first or last base), the exon is placed from the other at the reference's length (`inferred`).
// Where both are there further apart or closer than that (an indel inside the exon, up to half a unit), the
// exon spans them, with the haplotype's bases (`resized`).
Placed place_layout(const Walk& W, const Tandem& tm, const std::vector<std::size_t>& layout, std::size_t lo,
                    std::size_t hi) {
    std::vector<std::size_t> lens(tm.unit_len);
    std::nth_element(lens.begin(), lens.begin() + lens.size() / 2, lens.end());
    const std::size_t reach = lens[lens.size() / 2] / 2;
    Placed out;
    std::size_t pos = lo;
    const auto at = [&](std::size_t j, int node, char rev) { return W.node[j] == node && W.rev[j] == rev; };
    for (const std::size_t id : layout) {
        const Tandem::Kind& k = tm.kinds[id];
        const std::size_t slack = std::max<std::size_t>(30, k.len / 5);
        std::size_t a = SIZE_MAX, b = SIZE_MAX;
        bool exact = false, resized = false;
        for (std::size_t j = step_at(W, pos); j < W.size() && W.pref[j] < hi && a == SIZE_MAX; ++j) {
            if (!at(j, k.start_node, k.start_rev) || W.pref[j] + k.start_off < pos) continue;
            a = W.pref[j] + k.start_off;
            for (std::size_t e = j; e < W.size() && W.pref[e] <= a + k.len + std::max(slack, reach); ++e) {
                if (!at(e, k.end_node, k.end_rev)) continue;
                const std::size_t x = W.pref[e] + k.end_off + 1;
                if (x <= a || x > hi) continue;
                b = x;
                exact = x - a + slack >= k.len && x - a <= k.len + slack;
                resized = !exact;
                break;
            }
        }
        if (a != SIZE_MAX && b == SIZE_MAX && a + k.len <= hi) b = a + k.len;
        if (a == SIZE_MAX)   // no first node: the last one, at the reference's length before it
            for (std::size_t e = step_at(W, pos); e < W.size() && W.pref[e] < hi; ++e) {
                if (!at(e, k.end_node, k.end_rev)) continue;
                const std::size_t x = W.pref[e] + k.end_off + 1;
                if (x <= hi && x >= pos + k.len) { a = x - k.len; b = x; break; }
            }
        if (a == SIZE_MAX || b == SIZE_MAX) continue;
        out.pieces.push_back({a, b, k.ref_exon});
        (exact ? out.exact : resized ? out.resized : out.inferred) += 1;
        pos = b;
    }
    return out;
}

// The exons of a tandem stretch on walk W between bases lo and hi. The walk's units start where it holds the
// stretch's first exon start; a haplotype with the reference's number of units takes the reference's layouts
// in order, and otherwise each unit takes the layout that places most of its exons exactly, then the one
// of the reference unit sharing most of its nodes, then the one nearest in relative position.
Placed tandem_exons(const Walk& W, const Tandem& tm, std::size_t lo, std::size_t hi) {
    const Tandem::Kind& s = tm.kinds[tm.units[0][0]];
    std::vector<std::size_t> starts;
    for (std::size_t j = lo < W.len() ? step_at(W, lo) : W.size(); j < W.size() && W.pref[j] < hi; ++j)
        if (W.node[j] == s.start_node && W.rev[j] == s.start_rev && W.pref[j] + s.start_off >= lo &&
            W.pref[j] + s.start_off < hi)
            starts.push_back(W.pref[j] + s.start_off);
    Placed out;
    const std::size_t n = starts.size(), k = tm.units.size();
    for (std::size_t u = 0; u < n; ++u) {
        const std::size_t a = starts[u], b = u + 1 < n ? starts[u + 1] : hi;
        std::size_t pick = u;
        Placed best;
        if (n == k) {
            best = place_layout(W, tm, tm.units[u], a, b);
        } else {
            // The unit's nodes over at most a reference unit's length (the last unit runs to the stretch's end).
            std::vector<std::size_t> lens(tm.unit_len);
            std::nth_element(lens.begin(), lens.begin() + lens.size() / 2, lens.end());
            const std::size_t end = std::min(b, a + lens[lens.size() / 2]);
            std::vector<int> nodes;
            for (std::size_t j = step_at(W, a); j < W.size() && W.pref[j] < end; ++j) nodes.push_back(W.node[j]);
            std::sort(nodes.begin(), nodes.end());
            nodes.erase(std::unique(nodes.begin(), nodes.end()), nodes.end());
            double best_share = -1.0;
            for (std::size_t g = 0; g < k; ++g) {
                std::vector<int> common;
                std::set_intersection(nodes.begin(), nodes.end(), tm.unit_nodes[g].begin(), tm.unit_nodes[g].end(),
                                      std::back_inserter(common));
                const double share = static_cast<double>(common.size()) /
                                     static_cast<double>(nodes.size() + tm.unit_nodes[g].size() - common.size());
                const std::size_t near = u * k / n;
                Placed p = place_layout(W, tm, tm.units[g], a, b);
                const auto key = [&](const Placed& x, double sh, std::size_t gg) {
                    return std::make_tuple(x.exact, sh, -static_cast<long>(gg > near ? gg - near : near - gg));
                };
                if (g == 0 || key(p, share, g) > key(best, best_share, pick)) { best = std::move(p); best_share = share; pick = g; }
            }
        }
        // A walk unit missing the start (a variant there) runs over the next: its layout is placed again
        // while more than half a unit is left and the layout still places an exon.
        while (!best.pieces.empty()) {
            out.pieces.insert(out.pieces.end(), best.pieces.begin(), best.pieces.end());
            out.exact += best.exact;
            out.inferred += best.inferred;
            out.resized += best.resized;
            const std::size_t from = best.pieces.back().b;
            if (b < from || b - from <= tm.unit_len[pick] / 2) break;
            best = place_layout(W, tm, tm.units[pick], from, b);
        }
    }
    return out;
}

// One transcript through one pass.
HstRecord project_one(const Walk& R, const std::string& rseq, const Gene& g, const Pass& p, const Walk& W,
                      const std::string& wseq, const TranscriptFeature& tx, const Cuts& cuts,
                      const TranscriptsResult::EdgeClip& clip, const Tandem& tm) {
    HstRecord r;
    r.reverse = p.reverse;
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
        const bool before = cs < R.pref[g.gs + p.first], after = ce > R.pref[g.gs + p.last + 1];
        if (!exon_shared(R, g, p, cs, ce)) {
            missing = true;
            add_detail(r.detail, e + ((before && p.at_walk_start) || (after && p.at_walk_end) ? ":walk_end" : ":missing"));
            continue;
        }
        const CutOnWalk a = map_cut(R, g, p, W, cs, true), b = map_cut(R, g, p, W, ce, false);
        if (!a.ok || !b.ok || b.pos <= a.pos) {
            missing = true;
            add_detail(r.detail, e + ":missing");
            continue;
        }
        hx[k] = {a.pos, b.pos};
        if (a.shift != 0) { ++r.moved; add_detail(r.detail, e + (minus ? ".end:" : ".start:") + signed_text(a.shift)); }
        if (b.shift != 0) { ++r.moved; add_detail(r.detail, e + (minus ? ".start:" : ".end:") + signed_text(b.shift)); }
    }
    std::vector<Piece> pieces;
    for (std::size_t k = 0; k < n && !missing; ++k) {
        if (!in_tandem(k)) {
            pieces.push_back({hx[k].first, hx[k].second, k});
        } else if (k == tm.lo) {
            // Between the exons around the stretch, or the pass's ends where the stretch ends the transcript.
            const std::size_t lo = tm.lo > 0 ? hx[tm.lo - 1].second : W.pref[static_cast<std::size_t>(p.match[p.first])];
            const std::size_t hi = tm.hi + 1 < n ? hx[tm.hi + 1].first : W.pref[static_cast<std::size_t>(p.match[p.last]) + 1];
            const Placed placed = tandem_exons(W, tm, lo, hi);
            const std::vector<Piece>& found = placed.pieces;
            const std::size_t inferred = placed.inferred;
            if (found.empty()) {
                missing = true;
                add_detail(r.detail, "tandem:missing");
            }
            if (found.size() != tm.hi - tm.lo + 1)
                add_detail(r.detail, "tandem:" + std::to_string(tm.hi - tm.lo + 1) + ">" + std::to_string(found.size()));
            if (inferred) add_detail(r.detail, "tandem_inferred:" + std::to_string(inferred));
            if (placed.resized) add_detail(r.detail, "tandem_resized:" + std::to_string(placed.resized));
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
        r.exon_starts.push_back(p.reverse ? W.len() - x.b : x.a);
        r.exon_lengths.push_back(x.b - x.a);
    }
    const std::size_t lo = pieces.front().a, hi = pieces.back().b;
    r.path_start = p.reverse ? W.len() - hi : lo;
    r.path_end = p.reverse ? W.len() - lo : hi;
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
    for (const Cuts& c : cuts) tandem.push_back(find_tandem(R, c));
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
        for (std::size_t gi = 0; gi < genes.size(); ++gi) {
            auto& gp = gene_passes[gi];
            std::sort(gp.begin(), gp.end(), [](const Pass* a, const Pass* b) { return a->h_lo < b->h_lo; });
            for (const std::size_t k : genes[gi].tx) {
                if (gp.empty()) {
                    HstRecord r;
                    r.transcript = k;
                    r.haplotype = h;
                    r.status = HstStatus::NoCopy;
                    char buf[32];
                    std::snprintf(buf, sizeof buf, "%.3f", best_share[gi]);
                    r.detail = "best_share:" + std::string(buf);
                    if (!lost_to[gi].empty()) add_detail(r.detail, lost_to[gi]);
                    out.push_back(std::move(r));
                    continue;
                }
                for (std::size_t c = 0; c < gp.size(); ++c) {
                    const Pass& p = *gp[c];
                    HstRecord r = project_one(R, rseq, genes[gi], p, p.reverse ? Hr : H, p.reverse ? hrseq : hseq,
                                              res.transcripts[k], cuts[k], res.clips[k], tandem[k]);
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
        std::map<std::string, int> annotated_copies;   // per gene name
        for (const std::size_t a : annotated_on[h]) annotated_copies[res.annotated[a].t.gene_name] = res.annotated[a].copies;
        for (HstRecord& r : recs) {
            const auto it = annotated_copies.find(res.transcripts[r.transcript].gene_name);
            r.annotated_copies = it == annotated_copies.end() ? 0 : it->second;
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
            for (const auto& [x, y] : at.t.exons) r.sequence += hseq.substr(x, y - x);
            if (at.t.strand == '-') r.sequence = reverse_complement(r.sequence);
            const auto it = projected.find({at.t.transcript_id, at.copy});
            if (it != projected.end() && recs[it->second].sequence == r.sequence) {
                recs[it->second].source |= kAnnotated;
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
        }
        for (const auto& recs : got) sink.add(res, recs);
    }
    return res;
}

}  // namespace panvar
