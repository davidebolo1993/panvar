#include "panvar/copy_sites.hpp"

#include <algorithm>
#include <map>
#include <memory>
#include <set>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>

namespace panvar {
namespace {

// A crossing of a site, as bubble_crossings gives it: boundary steps left < right, and whether it runs source
// to sink in path order.
struct Cross {
    std::size_t site = 0;
    std::size_t left = 0, right = 0;
    bool source_to_sink = true;
};

// A crossing as a path is read, forwards or backwards (a haplotype stored on the other strand): the
// direction it reads the site in, and the steps it enters and leaves by.
struct Tok {
    const Cross* c = nullptr;
    bool fwd = true;
    std::size_t entry = 0, exit = 0;
};

std::vector<Tok> read_as(const std::vector<Cross>& cr, bool reversed) {
    std::vector<Tok> t;
    t.reserve(cr.size());
    for (const Cross& c : cr)
        t.push_back(reversed ? Tok{&c, !c.source_to_sink, c.right, c.left} : Tok{&c, c.source_to_sink, c.left, c.right});
    if (reversed) std::reverse(t.begin(), t.end());
    return t;
}

// How an alignment of a path's crossings to the reference's scores: matched crossings first, then fewest
// crossings left unmatched before the first match or after the last (a copy change has matched crossings on
// either side to be placed between), then how alike the matched copies are where the path holds fewer copies
// than the reference, then the earliest matches on the path, so the paths with one change agree on where it lies.
struct Score {
    bool valid = false;
    long long matches = 0, edge = 0;
    double sim = 0.0;
    long long at = 0;
};

bool better(const Score& a, const Score& b) {
    if (a.valid != b.valid) return a.valid;
    if (a.matches != b.matches) return a.matches > b.matches;
    if (a.edge != b.edge) return a.edge < b.edge;
    if (a.sim != b.sim) return a.sim > b.sim;
    return a.at < b.at;
}

struct Alignment {
    Score score;
    std::vector<std::pair<std::size_t, std::size_t>> pairs;   // (path token, reference token), increasing
};

// Monotone alignment of a path's crossings to the reference's: a pair is the same site read in the same
// direction. `sim(i, j)` scores a pair where it counts (see Score).
template <class Sim>
Alignment align(const std::vector<Tok>& h, const std::vector<Tok>& r, const Sim& sim) {
    const std::size_t nh = h.size(), nr = r.size();
    struct Cell { Score s; std::size_t pi = 0, pj = 0; };   // F: the best alignment ending in a pair here
    std::vector<Cell> f((nh + 1) * (nr + 1)), m((nh + 1) * (nr + 1));   // m: best F at or above-left, and where
    const auto at = [&](std::size_t i, std::size_t j) { return i * (nr + 1) + j; };
    for (std::size_t i = 1; i <= nh; ++i)
        for (std::size_t j = 1; j <= nr; ++j) {
            Cell& c = f[at(i, j)];
            if (h[i - 1].c->site == r[j - 1].c->site && h[i - 1].fwd == r[j - 1].fwd) {
                Score start{true, 0, static_cast<long long>(i - 1 + j - 1), 0.0, 0};
                const Cell& prev = m[at(i - 1, j - 1)];
                Score base = start;
                std::size_t pi = 0, pj = 0;
                if (prev.s.valid && better(prev.s, start)) { base = prev.s; pi = prev.pi; pj = prev.pj; }
                c.s = base;
                c.s.matches += 1;
                c.s.sim += sim(i - 1, j - 1);
                c.s.at += static_cast<long long>(i - 1);
                c.pi = pi; c.pj = pj;
            }
            Cell& best = m[at(i, j)];
            best = m[at(i - 1, j)];
            if (better(m[at(i, j - 1)].s, best.s)) best = m[at(i, j - 1)];
            if (better(c.s, best.s)) { best.s = c.s; best.pi = i; best.pj = j; }
        }
    Alignment out;
    std::size_t bi = 0, bj = 0;
    for (std::size_t i = 1; i <= nh; ++i)
        for (std::size_t j = 1; j <= nr; ++j) {
            Score s = f[at(i, j)].s;
            if (!s.valid) continue;
            s.edge += static_cast<long long>(nh - i + nr - j);
            if (better(s, out.score)) { out.score = s; bi = i; bj = j; }
        }
    while (bi > 0) {
        out.pairs.emplace_back(bi - 1, bj - 1);
        const Cell& c = f[at(bi, bj)];
        bi = c.pi; bj = c.pj;
    }
    std::reverse(out.pairs.begin(), out.pairs.end());
    return out;
}

// Length-weighted Jaccard of the interior nodes of two crossings.
double interior_similarity(const Graph& graph, const PathRecord& a, const Cross& x, const PathRecord& b, const Cross& y) {
    std::unordered_map<std::string, std::pair<std::size_t, std::size_t>> n;
    for (std::size_t k = x.left + 1; k < x.right; ++k) ++n[a.steps[k].node_id].first;
    for (std::size_t k = y.left + 1; k < y.right; ++k) ++n[b.steps[k].node_id].second;
    double lo = 0.0, hi = 0.0;
    for (const auto& [id, c] : n) {
        const auto it = graph.nodes.find(id);
        const double len = it == graph.nodes.end() ? 1.0 : static_cast<double>(it->second.sequence.size());
        lo += len * static_cast<double>(std::min(c.first, c.second));
        hi += len * static_cast<double>(std::max(c.first, c.second));
    }
    return hi > 0.0 ? lo / hi : 1.0;
}

PathWindow own_window(const PathRecord& p, std::size_t a, std::size_t b, bool forward) {
    return PathWindow{&p, std::min(a, b), std::max(a, b), forward};
}

} // namespace

std::vector<std::vector<Bubble>> copy_sites(const Graph& graph, const std::vector<BubblePathIndex>& path_indexes,
                                            const std::vector<Bubble>& bubbles, std::size_t ref_idx) {
    std::vector<std::vector<Bubble>> out(bubbles.size());
    if (ref_idx >= graph.paths.size()) return out;
    const std::size_t np = graph.paths.size();
    const PathRecord& ref = graph.paths[ref_idx];

    // Each path's crossings of the ordinary sites whose boundary orientations describe the walks, in path order.
    std::vector<std::vector<Cross>> cr(np);
    std::vector<char> multi(bubbles.size(), 0);
    for (std::size_t si = 0; si < bubbles.size(); ++si) {
        const Bubble& b = bubbles[si];
        if (b.inversion_parent || b.source == b.sink || b.windows ||
            !bubble_orientation_holds(graph.paths, path_indexes, b))
            continue;
        for (std::size_t pi = 0; pi < np; ++pi) {
            const auto c = bubble_crossings(graph.paths[pi], path_indexes[pi], b, false);
            if (c.size() > 1) multi[si] = 1;
            for (const BubblePathInterval& iv : c) cr[pi].push_back(Cross{si, iv.left, iv.right, iv.source_to_sink});
        }
    }
    if (std::none_of(multi.begin(), multi.end(), [](char c) { return c != 0; })) return out;
    // The reference's first and last nodes, where it walks each once, anchor a path that walks them once too: a
    // copy gained or lost at either end of the sites then still lies between two matched places.
    const std::size_t kStart = bubbles.size(), kEnd = bubbles.size() + 1;
    const auto is_multi = [&](std::size_t site) { return site < bubbles.size() && multi[site] != 0; };
    if (ref.steps.size() >= 2 && ref.steps.front().node_id != ref.steps.back().node_id) {
        const std::string& first = ref.steps.front().node_id;
        const std::string& last = ref.steps.back().node_id;
        for (std::size_t pi = 0; pi < np; ++pi) {
            const std::vector<PathStep>& st = graph.paths[pi].steps;
            std::size_t nf = 0, nl = 0, kf = 0, kl = 0;
            for (std::size_t k = 0; k < st.size(); ++k) {
                if (st[k].node_id == first) { ++nf; kf = k; }
                if (st[k].node_id == last) { ++nl; kl = k; }
            }
            if (pi == ref_idx && (nf != 1 || nl != 1)) {
                for (auto& v : cr) v.erase(std::remove_if(v.begin(), v.end(), [&](const Cross& c) { return c.site >= kStart; }), v.end());
                break;
            }
            if (nf == 1) cr[pi].push_back(Cross{kStart, kf, kf, st[kf].reverse == ref.steps.front().reverse});
            if (nl == 1) cr[pi].push_back(Cross{kEnd, kl, kl, st[kl].reverse == ref.steps.back().reverse});
        }
    }
    for (auto& v : cr)
        std::sort(v.begin(), v.end(), [](const Cross& a, const Cross& b) { return std::tie(a.left, a.right) < std::tie(b.left, b.right); });

    const std::vector<Tok> rt = read_as(cr[ref_idx], false);
    std::map<std::pair<std::size_t, bool>, std::size_t> ref_count;   // (site, direction) -> reference crossings
    for (const Tok& t : rt) ++ref_count[{t.c->site, t.fwd}];

    // Per path: its reading, and its tokens matched to each reference token.
    std::vector<std::vector<Tok>> reading(np);
    std::vector<Alignment> al(np);
    for (std::size_t pi = 0; pi < np; ++pi) {
        if (pi == ref_idx) {
            reading[pi] = rt;
            al[pi].score.valid = true;
            for (std::size_t j = 0; j < rt.size(); ++j) al[pi].pairs.emplace_back(j, j);
            continue;
        }
        if (cr[pi].empty()) continue;
        // A path walks the reference's first and last nodes the way it reads the locus: forwards where it walks
        // them as the reference does. Where it walks neither, or the two disagree, the reading that matches more
        // crossings is taken. An inversion between inverted copies reads almost as well backwards, the copies
        // then swapped, so the ends decide where they can.
        int ends = 0;   // +1 forwards, -1 backwards, 0 undecided
        for (const Cross& c : cr[pi])
            if (c.site >= kStart) {
                const int way = c.source_to_sink ? 1 : -1;
                ends = ends == 0 || ends == way ? way : 2;
            }
        for (const bool rev : {false, true}) {
            if ((ends == 1 && rev) || (ends == -1 && !rev)) continue;
            std::vector<Tok> h = read_as(cr[pi], rev);
            std::map<std::pair<std::size_t, bool>, std::size_t> hc;
            for (const Tok& t : h) ++hc[{t.c->site, t.fwd}];
            std::unordered_map<std::size_t, double> cache;
            const auto sim = [&](std::size_t i, std::size_t j) -> double {
                const auto key = std::make_pair(h[i].c->site, h[i].fwd);
                if (hc[key] >= ref_count[key]) return 0.0;
                const std::size_t k = i * rt.size() + j;
                const auto it = cache.find(k);
                if (it != cache.end()) return it->second;
                return cache[k] = interior_similarity(graph, graph.paths[pi], *h[i].c, ref, *rt[j].c);
            };
            Alignment a = align(h, rt, sim);
            if (!al[pi].score.valid || better(a.score, al[pi].score)) { al[pi] = std::move(a); reading[pi] = std::move(h); }
        }
    }

    // Copy sites: per reference crossing of a site crossed more than once, each path's window.
    using Windows = std::unordered_map<const PathRecord*, PathWindow>;
    std::vector<Windows> copy_win(rt.size());
    for (std::size_t j = 0; j < rt.size(); ++j)
        copy_win[j][&ref] = own_window(ref, rt[j].c->left, rt[j].c->right, rt[j].c->source_to_sink);
    // Gap sites: between two reference crossings, the paths whose copies between them differ in number.
    std::map<std::pair<std::size_t, std::size_t>, Windows> gap_win;
    std::vector<std::vector<std::size_t>> matched_of(np);   // per path, per reference token: its token, or none
    constexpr std::size_t kNone = static_cast<std::size_t>(-1);

    for (std::size_t pi = 0; pi < np; ++pi) {
        if (pi == ref_idx || !al[pi].score.valid) continue;
        const PathRecord& p = graph.paths[pi];
        const std::vector<Tok>& h = reading[pi];
        const auto& pairs = al[pi].pairs;
        matched_of[pi].assign(rt.size(), kNone);
        for (const auto& [i, j] : pairs) matched_of[pi][j] = i;
        // The path's window at each reference crossing it matches; widened below over crossings it cannot place.
        std::vector<std::pair<std::size_t, std::size_t>> win(rt.size(), {kNone, kNone});
        std::vector<std::size_t> held(rt.size(), 1);
        for (const auto& [i, j] : pairs) win[j] = {h[i].c->left, h[i].c->right};
        for (std::size_t k = 0; k <= pairs.size(); ++k) {
            const bool has_p = k > 0, has_n = k < pairs.size();
            const std::size_t i0 = has_p ? pairs[k - 1].first + 1 : 0, i1 = has_n ? pairs[k].first : h.size();
            const std::size_t j0 = has_p ? pairs[k - 1].second + 1 : 0, j1 = has_n ? pairs[k].second : rt.size();
            if (i0 == i1 && j0 == j1) {
                // Two crossings matched one after the other: where a site was read across copies, the stretch
                // between them lay inside its window. One that differs from the reference's is read as a gap.
                if (!has_p || !has_n || (!is_multi(rt[j0 - 1].c->site) && !is_multi(rt[j1].c->site)) ||
                    rt[j0 - 1].exit >= rt[j1].entry)
                    continue;
                const Tok& a = h[i0 - 1];
                const Tok& b = h[i1];
                const std::size_t lo = std::min(a.exit, b.entry), hi = std::max(a.exit, b.entry);
                bool same = hi - lo == rt[j1].entry - rt[j0 - 1].exit;
                for (std::size_t k = 1; same && lo + k < hi; ++k) {
                    const PathStep& x = a.exit < b.entry ? p.steps[lo + k] : p.steps[hi - k];
                    const PathStep& y = ref.steps[rt[j0 - 1].exit + k];
                    same = x.node_id == y.node_id && (x.reverse != (a.exit > b.entry)) == y.reverse;
                }
                if (!same) gap_win[{j0 - 1, j1}][&p] = own_window(p, a.exit, b.entry, a.exit < b.entry);
                continue;
            }
            bool involves = (has_p && is_multi(rt[j0 - 1].c->site)) || (has_n && is_multi(rt[j1].c->site));
            for (std::size_t i = i0; i < i1 && !involves; ++i) involves = is_multi(h[i].c->site);
            for (std::size_t j = j0; j < j1 && !involves; ++j) involves = is_multi(rt[j].c->site);
            if (!involves) continue;
            // Crossings the reference's read backwards, in reverse order, inside an inversion site: an inversion
            // there, which that site types from the walks. Read here too, it would be written twice.
            bool inverted = i1 - i0 == j1 - j0;
            for (std::size_t k = 0; inverted && i0 + k < i1; ++k) {
                const Tok& x = h[i0 + k];
                const Tok& y = rt[j1 - 1 - k];
                inverted = x.c->site == y.c->site && x.fwd != y.fwd && x.c->site < bubbles.size() &&
                           bubbles[x.c->site].parent_id != 0;
            }
            if (inverted) continue;
            // Between two matched crossings the change has a place: one gap site, read from the first's exit
            // to the second's entry on both walks. Where the two crossings meet at one boundary step on the
            // reference, nothing lies between: the gap starts at the first crossing's entry, and that crossing
            // is read there, not at its copy site.
            if (has_p && has_n) {
                const Tok& a = h[i0 - 1];
                const Tok& b = h[i1];
                const bool meet = rt[j0 - 1].exit >= rt[j1].entry;
                const std::size_t from = meet ? a.entry : a.exit;
                gap_win[{j0 - 1, j1}][&p] = own_window(p, from, b.entry, from < b.entry);
                if (meet) win[j0 - 1] = {kNone, kNone};
                continue;
            }
            // Otherwise a crossing left over goes into the window of a matched crossing of the same site and
            // direction, as one window from first to last crossing reads it.
            for (std::size_t i = i0; i < i1; ++i) {
                if (!is_multi(h[i].c->site)) continue;
                std::size_t best = kNone, dist = 0;
                for (std::size_t q = 0; q < pairs.size(); ++q) {
                    const std::size_t mi = pairs[q].first;
                    if (h[mi].c->site != h[i].c->site || h[mi].fwd != h[i].fwd) continue;
                    const std::size_t d = mi > i ? mi - i : i - mi;
                    if (best == kNone || d < dist) { best = q; dist = d; }
                }
                if (best == kNone) continue;
                auto& w = win[pairs[best].second];
                ++held[pairs[best].second];
                w.first = std::min(w.first, h[i].c->left);
                w.second = std::max(w.second, h[i].c->right);
            }
        }
        for (const auto& [i, j] : pairs) {
            if (!is_multi(rt[j].c->site) || win[j].first == kNone) continue;
            copy_win[j][&p] = own_window(p, win[j].first, win[j].second, h[i].c->source_to_sink);
            copy_win[j][&p].crossings = held[j];
        }
    }

    // A gap site's other paths, matching both its reference crossings, hold what the reference holds between them
    // (a path whose stretch differs is a carrier above) and read the reference's stretch: copies they hold
    // between are typed one at a time, not again here.
    for (auto& [key, w] : gap_win) {
        const auto [jp, jn] = key;
        const std::size_t from = rt[jp].exit >= rt[jn].entry ? rt[jp].entry : rt[jp].exit;
        const PathWindow rw = own_window(ref, from, rt[jn].entry, true);
        w[&ref] = rw;
        for (std::size_t pi = 0; pi < np; ++pi) {
            const PathRecord& p = graph.paths[pi];
            if (pi == ref_idx || w.count(&p) || matched_of[pi].empty()) continue;
            if (matched_of[pi][jp] != kNone && matched_of[pi][jn] != kNone) w[&p] = PathWindow{&ref, rw.first, rw.last, true};
        }
    }

    // The copy sites, per reference crossing in reference order, in place of their site.
    for (std::size_t j = 0; j < rt.size(); ++j) {
        const std::size_t si = rt[j].c->site;
        if (!is_multi(si)) continue;
        Bubble b = bubbles[si];
        b.windows = std::make_shared<const Windows>(std::move(copy_win[j]));
        out[si].push_back(std::move(b));
    }
    // The gap sites, each under the site of the reference crossing before it (after it, from the start). Its
    // boundaries are the steps the reference leaves that crossing by and enters the next by; its interior every
    // node a path reads between.
    for (auto& [key, w] : gap_win) {
        const auto [jp, jn] = key;
        const std::size_t si = rt[jp].c->site < bubbles.size() ? rt[jp].c->site : rt[jn].c->site;
        if (si >= bubbles.size()) continue;
        Bubble b;
        b.id = bubbles[si].id;
        const std::size_t from = rt[jp].exit >= rt[jn].entry ? rt[jp].entry : rt[jp].exit;
        b.source = ref.steps[from].node_id;
        b.source_reverse = ref.steps[from].reverse;
        b.sink = ref.steps[rt[jn].entry].node_id;
        b.sink_reverse = ref.steps[rt[jn].entry].reverse;
        std::set<std::string> inside;
        for (const auto& [path, pw] : w)
            for (std::size_t k = pw.first + 1; k < pw.last; ++k) inside.insert(pw.from->steps[k].node_id);
        inside.erase(b.source);
        inside.erase(b.sink);
        b.inside.assign(inside.begin(), inside.end());
        b.windows = std::make_shared<const Windows>(std::move(w));
        b.between_copies = true;
        if (out[si].empty() && !is_multi(si)) out[si].push_back(bubbles[si]);
        out[si].push_back(std::move(b));
    }
    return out;
}

} // namespace panvar
