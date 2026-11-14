#include "panvar/transcripts.hpp"

#include "panvar/cli_utils.hpp"
#include "panvar/graph_utils.hpp"
#include "panvar/parallel.hpp"
#include "panvar/ref_path.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
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
        case HstStatus::MissingExon: return "missing_exon";
        case HstStatus::NoCopy: return "no_copy";
    }
    return "ok";
}

namespace {

// Chaining: how far back a link may reach, in anchors and in bases, and what an indel costs per base.
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
    for (std::size_t a = 0; a < n; ++a) {
        f[a] = static_cast<double>(an[a].w);
        const std::size_t stop = a > kChainLookback ? a - kChainLookback : 0;
        for (std::size_t b = a; b-- > stop;) {
            if (an[b].i >= an[a].i) continue;
            const std::size_t dy = an[a].y - (an[b].y + an[b].w);
            if (dy > kChainMaxGap) break;
            const std::size_t dx = an[a].x - (an[b].x + an[b].w);
            if (dx > kChainMaxGap) continue;
            const double dd = std::fabs(static_cast<double>(dx) - static_cast<double>(dy));
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
    if (!d.empty()) d += ';';
    d += note;
}

std::string signed_text(long v) { return (v > 0 ? "+" : "") + std::to_string(v); }

// One transcript through one pass.
HstRecord project_one(const Walk& R, const std::string& rseq, const Gene& g, const Pass& p, const Walk& W,
                      const std::string& wseq, const TranscriptFeature& tx, const Cuts& cuts) {
    HstRecord r;
    r.reverse = p.reverse;
    const std::size_t n = cuts.size();
    const bool minus = tx.strand == '-';
    const auto exon_no = [&](std::size_t k) { return minus ? n - k : k + 1; };   // transcript order, 1-based
    std::vector<std::pair<std::size_t, std::size_t>> hx(n);
    bool missing = false;
    for (std::size_t k = 0; k < n; ++k) {
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
    for (std::size_t k = 1; k < n && !missing; ++k)
        if (hx[k].first < hx[k - 1].second) {
            missing = true;
            add_detail(r.detail, "e" + std::to_string(exon_no(k)) + ":overlaps_neighbour");
        }
    if (missing) {
        r.status = HstStatus::MissingExon;
        r.moved = 0;
        return r;
    }
    // Splice-site dinucleotides at each intron, compared with the reference's.
    for (std::size_t k = 0; k + 1 < n; ++k) {
        const std::size_t intron = minus ? n - 1 - k : k + 1;
        const std::size_t rl = cuts[k].second, rr = cuts[k + 1].first;
        const std::size_t hl = hx[k].second, hr = hx[k + 1].first;
        if (rr < rl + 4 || hr < hl + 4) continue;
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
    for (const auto& [a, b] : hx) r.sequence += wseq.substr(a, b - a);
    if (minus) r.sequence = reverse_complement(r.sequence);
    r.exons = n;
    for (std::size_t q = 0; q < n; ++q) {
        const auto [a, b] = hx[minus ? n - 1 - q : q];
        r.exon_starts.push_back(p.reverse ? W.len() - b : a);
        r.exon_lengths.push_back(b - a);
    }
    const std::size_t lo = hx.front().first, hi = hx.back().second;
    r.path_start = p.reverse ? W.len() - hi : lo;
    r.path_end = p.reverse ? W.len() - lo : hi;
    r.status = r.moved ? HstStatus::BoundaryMoved : r.splice_changed ? HstStatus::SpliceSiteChanged : HstStatus::Ok;
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
    res.reference = resolve_reference_path_name(graph, opt.reference_path, "transcripts");
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
    const Walk R = to_walk(*ref_path, t);
    const std::string rseq = spell(R, t);
    const ParsedReferencePath meta = parse_reference_path_label(ref_path->name, R.len(), ref_path->type == 'W');
    if (!meta.has_interval)
        throw std::runtime_error("transcripts: reference path '" + ref_path->name +
                                 "' carries no contig:start-end, so GTF coordinates cannot be placed on it");
    const std::size_t base = meta.region_start_1based;
    std::vector<TranscriptFeature> txs =
        parse_gtf_transcripts(opt.gtf_path, meta.chrom, base, base + R.len() - 1, &res.partial, "transcripts");

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
    // A gene's copies are told apart by its marker's coding exons where it has them: paralogs are compared
    // on what they encode, as untranslated ends differ more freely.
    std::vector<Cuts> marker_cuts;
    for (Gene& g : genes) {
        const TranscriptFeature& m = res.transcripts[g.marker];
        if (m.cds.empty()) { marker_cuts.push_back(cuts[g.marker]); continue; }
        Cuts c;
        for (const auto& [s, e] : m.cds) c.emplace_back(s - base, e - base + 1);
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
                // Chains are completed only when they hold a tenth of the best one's anchored bases: the
                // others are stray matches of single nodes, not passes.
                const std::vector<std::vector<Anchor>> chains = chain_anchors(R, g, W, occ, rev, H.size());
                std::vector<std::size_t> bp(chains.size(), 0);
                std::size_t best = 0;
                for (std::size_t k = 0; k < chains.size(); ++k) {
                    for (const Anchor& a : chains[k]) bp[k] += a.w;
                    best = std::max(best, bp[k]);
                }
                for (std::size_t k = 0; k < chains.size(); ++k) {
                    if (10 * bp[k] < best) continue;
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
                                              res.transcripts[k], cuts[k]);
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

    // Haplotypes in batches, so results reach the sink in path order without holding every sequence.
    const std::size_t threads = opt.threads ? opt.threads : std::max(1u, std::thread::hardware_concurrency());
    const std::size_t batch = threads * 4;
    for (std::size_t b0 = 0; b0 < graph.paths.size(); b0 += batch) {
        const std::size_t b1 = std::min(graph.paths.size(), b0 + batch);
        std::vector<std::vector<HstRecord>> got(b1 - b0);
        run_parallel(b1 - b0, threads, [&](std::size_t k) { got[k] = project_haplotype(b0 + k); });
        for (const auto& recs : got) sink.add(res, recs);
    }
    return res;
}

}  // namespace panvar
