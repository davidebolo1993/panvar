#include "panvar/module_fold.hpp"

#include "panvar/graph_utils.hpp"
#include "panvar/minimap2_align.hpp"
#include "panvar/parallel.hpp"

#include <edlib.h>

#include <algorithm>
#include <cctype>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <map>
#include <numeric>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <tuple>
#include <unordered_map>
#include <unordered_set>

namespace panvar {
namespace {

using Step = std::pair<std::string, bool>;   // node id, reverse
using Walk = std::vector<Step>;

struct FoldGraph {
    std::vector<std::string> order;                          // S lines in input order
    std::unordered_map<std::string, std::string> seq, tags;  // id -> sequence; id -> the rest of its S line
    std::vector<std::string> names;
    std::vector<Walk> walks;
};

std::vector<std::string> split_tabs(const std::string& line) {
    std::vector<std::string> f;
    std::size_t at = 0;
    while (true) {
        const std::size_t t = line.find('\t', at);
        f.push_back(line.substr(at, t == std::string::npos ? std::string::npos : t - at));
        if (t == std::string::npos) return f;
        at = t + 1;
    }
}

FoldGraph parse(const std::string& gfa) {
    FoldGraph g;
    std::istringstream in(gfa);
    std::string line;
    while (std::getline(in, line)) {
        if (line.size() < 2 || line[1] != '\t' || (line[0] != 'S' && line[0] != 'P')) continue;
        const std::vector<std::string> f = split_tabs(line);
        if (f.size() < 3) continue;
        if (line[0] == 'S') {
            g.order.push_back(f[1]);
            g.seq[f[1]] = f[2] == "*" ? std::string() : f[2];
            std::string rest;
            for (std::size_t k = 3; k < f.size(); ++k) rest += '\t' + f[k];
            g.tags[f[1]] = rest;
            continue;
        }
        Walk w;
        for (std::size_t s = 0; s < f[2].size();) {
            std::size_t e = f[2].find(',', s);
            if (e == std::string::npos) e = f[2].size();
            if (e - s >= 2) w.emplace_back(f[2].substr(s, e - s - 1), f[2][e - 1] == '-');
            s = e + 1;
        }
        g.names.push_back(f[1]);
        g.walks.push_back(std::move(w));
    }
    return g;
}

std::string spell(const FoldGraph& g, const Walk& w) {
    std::string out;
    for (const auto& [n, r] : w) out += r ? reverse_complement(g.seq.at(n)) : g.seq.at(n);
    return out;
}

Walk flip(const Walk& w) {
    Walk out(w.rbegin(), w.rend());
    for (auto& st : out) st.second = !st.second;
    return out;
}

bool near(std::size_t a, std::size_t b, double frac) {
    return static_cast<double>(a > b ? a - b : b - a) <= frac * static_cast<double>(b);
}

Minimap2Params params(const char* preset, std::size_t best_n, float pri, bool self = false) {
    Minimap2Params p;
    p.preset = preset; p.best_n = best_n; p.pri_ratio = pri; p.self = self;
    return p;
}

// ---- finding modules ----------------------------------------------------------------------------------

struct SelfHit { std::size_t qs, qe, ts, te, off; };

struct Array {
    std::size_t start = 0, end = 0, period = 0;
    std::string path, unit;
    double tiling = 0;   // the share of the array its unit's copies cover
};

// A sequence's arrays, from its same-strand self-alignments of min_copy_bp at min_identity (base-level: chains
// alone drift across a gene family). Hits whose spans overlap are one array, so a sequence can hold several.
// Candidate periods are the most recurring hit offsets and offset differences, the longest hit's offset, and
// fractions of each, since a multiple of the period tiles an array as well as the period does. The chosen
// period is the one whose unit, aligned back, covers the array best, then with more copies. The period is
// found before the size floor applies: an array whose period is under min_copy_bp is a short tandem repeat,
// left to panphorte, not folded at a multiple of its period. The unit is the stretch the hits share, at most
// one period.
std::vector<Array> find_arrays(const std::string& s, const std::string& name, const ModuleFoldOptions& o) {
    constexpr std::size_t kShortestPeriod = 100;   // fractions of a candidate tried down to this
    std::vector<SelfHit> all_hits;
    for (const Minimap2Hit& h : minimap2_hits(name, s, name, s, params("asm20", 64, 0.8f, true))) {
        if (h.reverse || h.query_start_bp >= h.target_start_bp) continue;
        if (h.query_end_bp - h.query_start_bp < o.min_copy_bp || h.gc_identity < o.min_identity) continue;
        all_hits.push_back({h.query_start_bp, h.query_end_bp, h.target_start_bp, h.target_end_bp,
                            h.target_start_bp - h.query_start_bp});
    }
    std::sort(all_hits.begin(), all_hits.end(), [](const SelfHit& x, const SelfHit& y) { return x.qs < y.qs; });
    std::vector<std::vector<SelfHit>> clusters;
    std::size_t reach = 0;
    for (const SelfHit& h : all_hits) {
        if (clusters.empty() || h.qs >= reach) clusters.emplace_back();
        clusters.back().push_back(h);
        reach = std::max(reach, h.te);
    }
    std::vector<Array> arrays;
    for (const std::vector<SelfHit>& hits : clusters) {
        std::vector<std::size_t> offs;
        for (const SelfHit& h : hits) offs.push_back(h.off);
        std::sort(offs.begin(), offs.end());
        offs.erase(std::unique(offs.begin(), offs.end()), offs.end());
        std::vector<std::size_t> vals = offs;
        for (std::size_t i = 0; i < offs.size(); ++i)
            for (std::size_t j = i + 1; j < offs.size(); ++j)
                if (offs[j] - offs[i] >= o.min_copy_bp) vals.push_back(offs[j] - offs[i]);
        auto support = [&](std::size_t p) {
            return std::count_if(vals.begin(), vals.end(), [&](std::size_t v) { return near(v, p, 0.02); });
        };
        auto hits_at = [&](std::size_t p) {
            std::vector<SelfHit> at;
            for (const SelfHit& h : hits) if (near(h.off, p, 0.02)) at.push_back(h);
            return at;
        };
        std::vector<std::size_t> ranked(vals.begin(), vals.end());
        std::sort(ranked.begin(), ranked.end());
        ranked.erase(std::unique(ranked.begin(), ranked.end()), ranked.end());
        std::stable_sort(ranked.begin(), ranked.end(), [&](std::size_t x, std::size_t y) { return support(x) > support(y); });
        std::vector<std::pair<std::size_t, std::vector<SelfHit>>> cands;
        auto taken = [&](std::size_t p) {
            return std::any_of(cands.begin(), cands.end(), [&](const auto& c) { return near(p, c.first, 0.02); });
        };
        for (std::size_t p : ranked) {
            if (cands.size() == 3) break;
            auto at = hits_at(p);
            if (!at.empty() && !taken(p)) cands.emplace_back(p, std::move(at));
        }
        std::vector<SelfHit> by_len = hits;
        std::sort(by_len.begin(), by_len.end(), [](const SelfHit& x, const SelfHit& y) { return x.qe - x.qs > y.qe - y.qs; });
        for (const SelfHit& h : by_len) {
            if (cands.size() >= 4) break;
            if (!taken(h.off)) cands.emplace_back(h.off, hits_at(h.off));
        }
        std::vector<std::pair<std::size_t, const std::vector<SelfHit>*>> trials;
        for (const auto& c : cands)
            for (std::size_t k = 1; k <= 6 && c.first / k >= kShortestPeriod; ++k) {
                const std::size_t q = c.first / k;
                if (std::none_of(trials.begin(), trials.end(), [&](const auto& t) { return near(q, t.first, 0.02); }))
                    trials.emplace_back(q, &c.second);
            }
        bool found = false;
        Array best;
        std::tuple<double, std::size_t, long long> best_key{-1.0, 0, 0};
        for (const auto& [p, atp] : trials) {
            std::size_t start = SIZE_MAX, end = 0, shared = 0;
            for (const SelfHit& h : *atp) {
                start = std::min(start, h.qs); end = std::max(end, h.te); shared = std::max(shared, h.qe - h.qs);
            }
            const std::size_t ulen = std::min(p, shared);
            if (start + ulen > s.size()) continue;
            const std::string unit = s.substr(start, ulen);
            std::vector<Minimap2Hit> all = minimap2_hits("u", unit, name, s, params("map-hifi", 200, 0.1f));
            std::sort(all.begin(), all.end(), [](const Minimap2Hit& x, const Minimap2Hit& y) {
                return x.target_end_bp - x.target_start_bp > y.target_end_bp - y.target_start_bp;
            });
            std::vector<const Minimap2Hit*> copies;
            for (const Minimap2Hit& h : all) {
                if (2 * (h.query_end_bp - h.query_start_bp) < ulen || h.gc_identity < o.min_identity) continue;
                if (std::none_of(copies.begin(), copies.end(), [&](const Minimap2Hit* t) {
                        return h.target_start_bp < t->target_end_bp && t->target_start_bp < h.target_end_bp; }))
                    copies.push_back(&h);
            }
            std::size_t lo = SIZE_MAX, hi = 0, covered = 0, n = 0;
            for (const Minimap2Hit* h : copies) {
                if (h->target_start_bp + ulen < start || h->target_end_bp > end + ulen) continue;
                lo = std::min(lo, h->target_start_bp); hi = std::max(hi, h->target_end_bp);
                covered += h->target_end_bp - h->target_start_bp; ++n;
            }
            if (n < 2) continue;
            const double tiling = std::round(100.0 * static_cast<double>(covered) / static_cast<double>(hi - lo)) / 100.0;
            const std::tuple<double, std::size_t, long long> key{tiling, n, -static_cast<long long>(p)};
            if (!found || key > best_key) {
                found = true; best_key = key;
                best.start = start; best.end = end; best.period = p; best.unit = unit; best.path = name; best.tiling = tiling;
            }
        }
        if (found && best.period >= o.min_copy_bp && best.unit.size() >= o.min_copy_bp) arrays.push_back(std::move(best));
    }
    return arrays;
}

using Module = ModuleUnit;

// Repeated sequence in s: a sample of its 31-mers (one in sixteen, by hash) that occur in it more than once.
std::unordered_set<std::uint64_t> repeated_kmers(const std::string& s) {
    constexpr std::size_t k = 31;
    std::unordered_map<std::uint64_t, int> seen;
    if (s.size() < k) return {};
    for (std::size_t i = 0; i + k <= s.size(); ++i) {
        const std::uint64_t h = std::hash<std::string_view>{}(std::string_view(s).substr(i, k));
        if ((h & 15) == 0) ++seen[h];
    }
    std::unordered_set<std::uint64_t> rep;
    for (const auto& [h, c] : seen) if (c > 1) rep.insert(h);
    return rep;
}

// Arrays from the reference, one haplotype per kilobase length class (a copy-number change changes a
// haplotype's length), and any haplotype repeating sequence those do not (an array can change without the
// length changing, when another changes the other way): half a module's worth of new repeated k-mers. Each
// sequence's arrays are placed on the reference by aligning their unit to it, at the copy identity floor, and
// grouped where their reference intervals overlap by half of the longer, transitively. A module's unit is the reference's own where the
// reference has an array there, else that of the period most arrays have, cut to a period another array of the
// module divides.
std::vector<Module> discover(const std::vector<std::string>& names, const std::vector<std::string>& seqs,
                             std::size_t ref, const ModuleFoldOptions& o) {
    std::vector<std::size_t> reps{ref};
    std::set<long long> classes{std::llround(static_cast<double>(seqs[ref].size()) / 1000.0)};
    for (std::size_t i = 0; i < seqs.size(); ++i)
        if (i != ref && classes.insert(std::llround(static_cast<double>(seqs[i].size()) / 1000.0)).second)
            reps.push_back(i);
    {
        std::vector<std::unordered_set<std::uint64_t>> rk(seqs.size());
        run_parallel(seqs.size(), o.threads, [&](std::size_t i) { rk[i] = repeated_kmers(seqs[i]); });
        std::unordered_set<std::uint64_t> covered;
        std::vector<char> is_rep(seqs.size(), 0);
        for (std::size_t r : reps) { is_rep[r] = 1; covered.insert(rk[r].begin(), rk[r].end()); }
        const std::size_t need = std::max<std::size_t>(1, o.min_copy_bp / 32);   // half a module, one k-mer in 16
        for (std::size_t i = 0; i < seqs.size(); ++i) {
            if (is_rep[i]) continue;
            std::size_t fresh = 0;
            for (const std::uint64_t h : rk[i]) fresh += covered.count(h) ? 0 : 1;
            if (fresh < need) continue;
            reps.push_back(i);
            covered.insert(rk[i].begin(), rk[i].end());
        }
    }
    std::vector<std::vector<Array>> found(reps.size());
    run_parallel(reps.size(), o.threads, [&](std::size_t k) { found[k] = find_arrays(seqs[reps[k]], names[reps[k]], o); });
    std::vector<Array> arrays;
    std::vector<std::size_t> from;   // the sampled sequence each array is from
    for (std::size_t k = 0; k < reps.size(); ++k)
        for (Array& a : found[k]) { arrays.push_back(std::move(a)); from.push_back(reps[k]); }
    std::vector<std::vector<Minimap2Hit>> placed(arrays.size());
    // Where each array's own flanks, anchor_bp of its sequence on either side, align once on the reference: the
    // reference positions of its start and end, SIZE_MAX where a flank does not place uniquely.
    std::vector<std::array<std::size_t, 2>> flank_at(arrays.size(), {SIZE_MAX, SIZE_MAX});
    run_parallel(arrays.size(), o.threads, [&](std::size_t i) {
        if (from[i] != ref) {
            // Both presets: map-hifi keeps a copy with a long indel in one alignment, asm20 reaches copies down
            // to the 90% floor; the span is taken over both.
            placed[i] = minimap2_hits("u", arrays[i].unit, names[ref], seqs[ref], params("map-hifi", 50, 0.1f));
            for (Minimap2Hit& h : minimap2_hits("u", arrays[i].unit, names[ref], seqs[ref], params("asm20", 50, 0.1f)))
                placed[i].push_back(std::move(h));
            const std::string& s = seqs[from[i]];
            const Array& a = arrays[i];
            auto place = [&](std::size_t b, std::size_t e, bool before) {
                if (e - b < o.anchor_bp / 2) return SIZE_MAX;
                std::vector<Minimap2Hit> one;
                for (Minimap2Hit& h : minimap2_hits("f", s.substr(b, e - b), names[ref], seqs[ref], params("map-hifi", 10, 0.5f)))
                    if (h.gc_identity >= o.anchor_identity && 5 * (h.query_end_bp - h.query_start_bp) >= 4 * (e - b))
                        one.push_back(std::move(h));
                if (one.size() != 1) return SIZE_MAX;
                // the side of the flank next to the array, on the reference's strand
                return before != one[0].reverse ? one[0].target_end_bp : one[0].target_start_bp;
            };
            flank_at[i] = {place(a.start > o.anchor_bp ? a.start - o.anchor_bp : 0, a.start, true),
                           place(a.end, std::min(s.size(), a.end + o.anchor_bp), false)};
        }
    });
    struct Cand { Array a; bool is_ref; std::size_t lo, hi; };
    std::vector<Cand> cands;
    for (std::size_t i = 0; i < arrays.size(); ++i) {
        const Array& a = arrays[i];
        if (from[i] == ref) { cands.push_back({a, true, a.start, a.end}); continue; }
        auto covering = [&](double c) {
            std::vector<const Minimap2Hit*> v;
            for (const Minimap2Hit& h : placed[i])
                if (static_cast<double>(h.query_end_bp - h.query_start_bp) >= c * static_cast<double>(a.unit.size()) &&
                    h.gc_identity >= o.min_copy_identity)
                    v.push_back(&h);
            return v;
        };
        if (covering(0.5).empty()) continue;
        // The span: runs of hits covering a quarter of the unit, each within a unit of the next, that hold a hit
        // covering half of it. The array's place on the reference is where its flanks align: the runs within a
        // unit of that place make the span (a spacer between two copies is in the flanks' interval), else the
        // nearest run. Where neither flank places, the run with the most matching bases. Another copy of the unit
        // elsewhere on the reference, full or partial, would otherwise stretch the span across every array between.
        // The array's ends are its self-alignments', which can run on through a lower-identity copy close beside
        // it: that copy is then one of the array's.
        std::vector<const Minimap2Hit*> q = covering(0.25);
        std::sort(q.begin(), q.end(), [](const Minimap2Hit* x, const Minimap2Hit* y) {
            return x->target_start_bp < y->target_start_bp; });
        const double half = 0.5 * static_cast<double>(a.unit.size());
        struct Run { std::size_t lo, hi, matches; };
        std::vector<Run> runs;
        for (std::size_t k = 0; k < q.size();) {
            Run r{q[k]->target_start_bp, q[k]->target_end_bp, 0};
            bool holds = false;
            for (; k < q.size() && q[k]->target_start_bp <= r.hi + a.unit.size(); ++k) {
                r.hi = std::max(r.hi, q[k]->target_end_bp);
                r.matches += q[k]->n_matches;
                holds = holds || static_cast<double>(q[k]->query_end_bp - q[k]->query_start_bp) >= half;
            }
            if (holds) runs.push_back(r);
        }
        std::size_t lo = SIZE_MAX, hi = 0;
        const auto [fs, fe] = flank_at[i];
        if (fs != SIZE_MAX || fe != SIZE_MAX) {
            const std::size_t at_lo = std::min(fs, fe), at_hi = fs == SIZE_MAX ? fe : fe == SIZE_MAX ? fs : std::max(fs, fe);
            auto dist = [&](const Run& r) {
                return r.hi < at_lo ? at_lo - r.hi : r.lo > at_hi ? r.lo - at_hi : std::size_t{0};
            };
            const Run* nearest = nullptr;
            for (const Run& r : runs) {
                if (dist(r) <= a.unit.size()) { lo = std::min(lo, r.lo); hi = std::max(hi, r.hi); }
                if (!nearest || dist(r) < dist(*nearest)) nearest = &r;
            }
            if (lo == SIZE_MAX && nearest) { lo = nearest->lo; hi = nearest->hi; }
        } else {
            std::size_t best = 0;
            for (const Run& r : runs)
                if (lo == SIZE_MAX || r.matches > best) { lo = r.lo; hi = r.hi; best = r.matches; }
        }
        cands.push_back({a, false, lo, hi});
    }
    std::vector<std::size_t> parent(cands.size());
    std::iota(parent.begin(), parent.end(), 0);
    std::function<std::size_t(std::size_t)> root = [&](std::size_t x) {
        while (parent[x] != x) x = parent[x] = parent[parent[x]];
        return x;
    };
    for (std::size_t i = 0; i < cands.size(); ++i)
        for (std::size_t j = 0; j < i; ++j) {
            const std::size_t lo = std::max(cands[i].lo, cands[j].lo), hi = std::min(cands[i].hi, cands[j].hi);
            const std::size_t longer = std::max(cands[i].hi - cands[i].lo, cands[j].hi - cands[j].lo);
            if (hi > lo && 2 * (hi - lo) >= longer) parent[root(i)] = root(j);
        }
    std::map<std::size_t, std::vector<std::size_t>> groups;
    for (std::size_t i = 0; i < cands.size(); ++i) groups[root(i)].push_back(i);
    std::vector<Module> modules;
    for (auto& [r, members] : groups) {
        std::sort(members.begin(), members.end(), [&](std::size_t a, std::size_t b) {
            return cands[a].a.period < cands[b].a.period; });
        std::vector<std::vector<std::size_t>> clusters;   // members by period, to within a tenth
        for (std::size_t m : members) {
            if (!clusters.empty() && cands[m].a.period - cands[clusters.back()[0]].a.period <=
                                         cands[clusters.back()[0]].a.period / 10)
                clusters.back().push_back(m);
            else clusters.push_back({m});
        }
        std::optional<std::size_t> own;
        for (std::size_t m : members) if (cands[m].is_ref) own = m;
        std::size_t rep = members[0];
        if (own) {
            rep = *own;
        } else {
            const std::vector<std::size_t>* best = &clusters[0];
            for (const auto& cl : clusters) if (cl.size() > best->size()) best = &cl;
            rep = (*best)[0];
            const auto copies_in = [&](std::size_t x) {
                return static_cast<double>(cands[x].a.end - cands[x].a.start) / static_cast<double>(cands[x].a.period);
            };
            for (std::size_t m : *best) if (copies_in(m) > copies_in(rep)) rep = m;
        }
        Module mod{cands[rep].a.unit, cands[rep].a.path, cands[rep].a.period, cands[rep].lo, cands[rep].hi};
        // A shorter period another array divides the module's by is its true unit only if that array's copies
        // tile it at least as well: one tiling half its array is a part the copies share, not their period.
        for (const auto& cl : clusters) {
            const std::size_t q = cands[cl[0]].a.period;
            if (q >= mod.period) continue;
            if (std::none_of(cl.begin(), cl.end(), [&](std::size_t m) { return cands[m].a.tiling >= cands[rep].a.tiling; }))
                continue;
            const long long k = std::llround(static_cast<double>(mod.period) / static_cast<double>(q));
            if (k >= 2 && std::fabs(static_cast<double>(mod.period) - static_cast<double>(k) * static_cast<double>(q)) <=
                              0.03 * static_cast<double>(mod.period)) {
                mod.unit = mod.unit.substr(0, std::min(q, mod.unit.size()));
                mod.period = q;
                break;
            }
        }
        modules.push_back(std::move(mod));
    }
    std::sort(modules.begin(), modules.end(), [](const Module& a, const Module& b) { return a.lo < b.lo; });
    return modules;
}

// ---- walk items -----------------------------------------------------------------------------------------

// U: template stretch [a,b); V: a difference replacing template [a,b) with `seq`; J: a jump within the
// template between two aligned pieces, with the sequence between them; D: the place of a long difference whose
// sequence follows as items on its own template; L: literal sequence. `t` is the template a U or V item lies
// on when it is not its piece's.
struct Item { char kind; std::size_t a, b; std::string seq; int t = -1; };

// A forward alignment's pieces against a template, as items. Each insertion or deletion is shifted left as
// far as the sequence allows (bounded by the previous difference), so one length change in a repeat is one
// node whichever offset the aligner chose.
std::vector<Item> aligned_items(const std::vector<const Minimap2Hit*>& pieces, const std::string& query,
                                const std::string& tmpl) {
    std::vector<Item> items;
    bool first = true;
    std::size_t prev_te = 0, prev_qe = 0;
    for (const Minimap2Hit* p : pieces) {
        if (!first) items.push_back({'J', prev_te, p->target_start_bp, query.substr(prev_qe, p->query_start_bp - prev_qe)});
        first = false;
        std::size_t u = p->target_start_bp, q = p->query_start_bp, run = u;
        for (const auto& [op, len_i] : p->cigar) {
            const std::size_t len = static_cast<std::size_t>(len_i);
            if (op == 7) { u += len; q += len; continue; }
            std::size_t at = u;
            std::string seq;
            if (op == 1) {          // insertion
                seq = query.substr(q, len);
                while (at > run && tmpl[at - 1] == seq.back()) { seq = tmpl[at - 1] + seq.substr(0, seq.size() - 1); --at; }
            } else if (op == 2) {   // deletion
                while (at > run && tmpl[at - 1] == tmpl[at + len - 1]) --at;
            }
            if (at > run) items.push_back({'U', run, at, {}});
            if (op == 8) {          // substitutions, one base each
                for (std::size_t k = 0; k < len; ++k) items.push_back({'V', u + k, u + k + 1, query.substr(q + k, 1)});
                u += len; q += len; run = u;
            } else if (op == 1) {
                items.push_back({'V', at, at, seq});
                q += len; run = at;
            } else if (op == 2) {
                items.push_back({'V', at, at + len, {}});
                u += len; run = at + len;
            }
        }
        if (u > run) items.push_back({'U', run, u, {}});
        if (u != p->target_end_bp || q != p->query_end_bp)
            throw std::runtime_error("an alignment's CIGAR does not cover it");
        prev_te = p->target_end_bp; prev_qe = p->query_end_bp;
    }
    return items;
}

struct Occurrence { std::vector<const Minimap2Hit*> pieces; std::size_t qs = 0, qe = 0; double cover = 0, gcid = 0; bool fwd = true; };

// A stretch's alignments to the unit, longest first where they do not overlap on the stretch, chained in
// stretch order into occurrences while they run on along the unit on one strand. A step back along the unit
// of more than `back` starts a new occurrence.
std::vector<Occurrence> occurrences(const std::vector<Minimap2Hit>& all, std::size_t unit_len, std::size_t back) {
    std::vector<const Minimap2Hit*> by_len;
    for (const Minimap2Hit& h : all) by_len.push_back(&h);
    std::sort(by_len.begin(), by_len.end(), [](const Minimap2Hit* a, const Minimap2Hit* b) {
        return a->query_end_bp - a->query_start_bp > b->query_end_bp - b->query_start_bp; });
    std::vector<const Minimap2Hit*> taken;
    for (const Minimap2Hit* h : by_len)
        if (std::none_of(taken.begin(), taken.end(), [&](const Minimap2Hit* t) {
                return h->query_start_bp < t->query_end_bp && t->query_start_bp < h->query_end_bp; }))
            taken.push_back(h);
    std::sort(taken.begin(), taken.end(), [](const Minimap2Hit* a, const Minimap2Hit* b) { return a->query_start_bp < b->query_start_bp; });
    std::vector<Occurrence> occ;
    for (const Minimap2Hit* h : taken) {
        const bool fwd = !h->reverse;
        if (!occ.empty() && occ.back().fwd == fwd) {
            const Minimap2Hit* last = occ.back().pieces.back();
            if (fwd ? h->target_start_bp + back >= last->target_end_bp : h->target_end_bp <= last->target_start_bp + back) {
                occ.back().pieces.push_back(h);
                continue;
            }
        }
        occ.emplace_back();
        occ.back().fwd = fwd;
        occ.back().pieces.push_back(h);
    }
    for (Occurrence& oc : occ) {
        std::vector<std::pair<std::size_t, std::size_t>> iv;
        double w = 0, wid = 0;
        oc.qs = SIZE_MAX;
        for (const Minimap2Hit* p : oc.pieces) {
            iv.emplace_back(p->target_start_bp, p->target_end_bp);
            const double len = static_cast<double>(p->query_end_bp - p->query_start_bp);
            w += len; wid += len * p->gc_identity;
            oc.qs = std::min(oc.qs, p->query_start_bp); oc.qe = std::max(oc.qe, p->query_end_bp);
        }
        std::sort(iv.begin(), iv.end());
        std::size_t tot = 0, end = 0;
        for (const auto& [a, b] : iv) if (b > end) { tot += b - std::max(a, end); end = b; }
        oc.cover = static_cast<double>(tot) / static_cast<double>(unit_len);
        oc.gcid = w > 0 ? wid / w : 0.0;
    }
    return occ;
}

// ---- folding one module -------------------------------------------------------------------------------

std::string fresh_id(const FoldGraph& g, std::uint64_t& next) {
    std::string id;
    do { id = std::to_string(next++); } while (g.seq.count(id));
    return id;
}

// A piece's tags: the `mn:Z:` name's range (`s7:100-400`, `fold_flank:1:0-3000`) is narrowed to the piece's
// bases [from, to) of the node, and a name without one (`s7`, `tip`) gets the piece's range, so the tag says
// which bases of the named sequence the piece holds. Other tags are kept as they are.
std::string piece_tags(const std::string& tags, std::size_t from, std::size_t to) {
    const std::size_t at = tags.find("mn:Z:");
    if (at == std::string::npos) return tags;
    const std::size_t name_end = std::min(tags.find('\t', at), tags.size());
    std::string mn = tags.substr(at + 5, name_end - at - 5);
    std::size_t base = 0;
    const std::size_t colon = mn.rfind(':'), dash = mn.rfind('-');
    auto digits = [&](std::size_t a, std::size_t b) {
        return a < b && std::all_of(mn.begin() + static_cast<long>(a), mn.begin() + static_cast<long>(b),
                                    [](char c) { return std::isdigit(static_cast<unsigned char>(c)) != 0; });
    };
    if (colon != std::string::npos && dash != std::string::npos && colon < dash && digits(colon + 1, dash) &&
        digits(dash + 1, mn.size())) {
        base = std::stoull(mn.substr(colon + 1, dash - colon - 1));
        mn.resize(colon);
    }
    return tags.substr(0, at + 5) + mn + ":" + std::to_string(base + from) + "-" + std::to_string(base + to) +
           tags.substr(name_end);
}

// Split node n after its k-th base, for every path walking it; returns the two pieces' ids.
std::pair<std::string, std::string> split_node(FoldGraph& g, const std::string& n, std::size_t k, std::uint64_t& next) {
    const std::string a = fresh_id(g, next);
    g.seq[a] = g.seq[n].substr(0, k);
    const std::string b = fresh_id(g, next);
    g.seq[b] = g.seq[n].substr(k);
    g.tags[a] = piece_tags(g.tags[n], 0, k);
    g.tags[b] = piece_tags(g.tags[n], k, g.seq[n].size());
    for (Walk& w : g.walks) {
        Walk out;
        for (const Step& st : w) {
            if (st.first != n) out.push_back(st);
            else if (st.second) { out.emplace_back(b, true); out.emplace_back(a, true); }
            else { out.emplace_back(a, false); out.emplace_back(b, false); }
        }
        w = std::move(out);
    }
    return {a, b};
}

// The cut nodes around [lo, hi) on the reference: backbone nodes, which the reference walks once forward, no
// path walks twice, and every path whose reference span covers them walks (a copy's own node, which haplotypes
// lacking the copy pass by, is not one; a contig ending before a node does not rule it out). A node of a module
// already folded (`copy_nodes`) is never one. One holding a boundary is split there (both boundaries in one
// node: split at hi, then at lo); otherwise the nearest such node outside is used. Returns false when a side
// has none.
bool cut_nodes(FoldGraph& g, std::size_t ref, std::size_t lo, std::size_t hi, std::uint64_t& next,
               const std::unordered_set<std::string>& copy_nodes, std::string& L, std::string& R) {
    auto candidates = [&]() {
        std::unordered_map<std::string, std::size_t> ref_count, at_of;   // the reference's nodes and their offsets
        {
            std::size_t at = 0;
            for (const Step& st : g.walks[ref]) { ++ref_count[st.first]; at_of.emplace(st.first, at); at += g.seq[st.first].size(); }
        }
        std::unordered_set<std::string> twice;
        std::vector<std::unordered_set<std::string>> walked(g.walks.size());
        std::vector<std::pair<std::size_t, std::size_t>> span(g.walks.size(), {SIZE_MAX, 0});   // reference offsets reached
        for (std::size_t p = 0; p < g.walks.size(); ++p) {
            for (const Step& st : g.walks[p]) {
                if (!walked[p].insert(st.first).second) twice.insert(st.first);
                if (ref_count[st.first] == 1) {
                    const std::size_t x = at_of[st.first];
                    span[p].first = std::min(span[p].first, x);
                    span[p].second = std::max(span[p].second, x + g.seq[st.first].size());
                }
            }
        }
        std::unordered_set<std::string> ok;
        for (const Step& st : g.walks[ref]) {
            if (st.second || ref_count[st.first] != 1 || twice.count(st.first) || copy_nodes.count(st.first)) continue;
            const std::size_t x = at_of[st.first], e = x + g.seq[st.first].size();
            bool backbone = true;
            for (std::size_t p = 0; p < g.walks.size() && backbone; ++p)
                if (span[p].first <= x && span[p].second >= e && !walked[p].count(st.first)) backbone = false;
            if (backbone) ok.insert(st.first);
        }
        return ok;
    };
    std::unordered_set<std::string> ok = candidates();
    std::optional<std::pair<std::string, std::size_t>> hold_lo, hold_hi;
    std::size_t at = 0;
    for (const Step& st : g.walks[ref]) {
        const std::size_t len = g.seq[st.first].size();
        if (ok.count(st.first) && at < lo && lo < at + len) hold_lo = std::make_pair(st.first, lo - at);
        if (ok.count(st.first) && at < hi && hi < at + len) hold_hi = std::make_pair(st.first, hi - at);
        at += len;
    }
    L.clear(); R.clear();
    std::string left_of_hi;
    if (hold_hi) std::tie(left_of_hi, R) = split_node(g, hold_hi->first, hold_hi->second, next);
    if (hold_lo) {
        const bool same = hold_hi && hold_hi->first == hold_lo->first;   // lo lies in the left piece now
        L = split_node(g, same ? left_of_hi : hold_lo->first, hold_lo->second, next).first;
    }
    if (L.empty() || R.empty()) {
        ok = candidates();
        at = 0;
        std::string left, right;
        for (const Step& st : g.walks[ref]) {
            const std::size_t len = g.seq[st.first].size();
            if (ok.count(st.first) && at + len <= lo) left = st.first;
            if (ok.count(st.first) && at >= hi && right.empty()) right = st.first;
            at += len;
        }
        if (L.empty()) L = left;
        if (R.empty()) R = right;
    }
    return !L.empty() && !R.empty();
}

struct Piece { int tmpl; std::vector<Item> items; std::string lit; };   // tmpl 0: a copy; 1-3: sequence around copies

// `copy_nodes` holds the copy nodes of the modules folded so far; this module's are added when it folds.
void fold_one(FoldGraph& folded, const std::vector<std::string>& own, std::size_t ref, Module mod,
              const ModuleFoldOptions& o, std::uint64_t& next, std::unordered_set<std::string>& copy_nodes,
              FoldedModule& rep) {
    FoldGraph g = folded;   // worked on a copy: a module not folded leaves no split node behind
    const std::string ref_seq = spell(g, g.walks[ref]);
    rep.period = mod.period; rep.unit_bp = mod.unit.size(); rep.unit_source = mod.source; rep.paths = g.walks.size();
    rep.ref_start = mod.lo; rep.ref_end = mod.hi;   // until cut, the candidate's span, so a refusal says where
    const std::size_t cand_lo = mod.lo, cand_hi = mod.hi;
    // An alignment's ends are uncertain by a few bases: a copy inserted right at the boundary must fall inside the
    // module, so the span is widened by 1% of the unit, at least 50 bp, before it is anchored and cut.
    {
        const std::size_t pad = std::max<std::size_t>(50, mod.unit.size() / 100);
        mod.lo = mod.lo > pad ? mod.lo - pad : 0;
        mod.hi = std::min(mod.hi + pad, ref_seq.size());
    }

    // Anchors unique on the reference, moved outward a kilobase at a time.
    auto unique = [&](std::size_t a, std::size_t b) {
        const std::string probe = ref_seq.substr(a, b - a);
        std::size_t n = 0;
        for (const Minimap2Hit& h : minimap2_hits("a", probe, "ref", ref_seq, params("map-hifi", 10, 0.5f)))
            if (h.gc_identity >= o.anchor_identity && 5 * (h.query_end_bp - h.query_start_bp) >= 4 * probe.size()) ++n;
        return n == 1;
    };
    const std::size_t lo0 = mod.lo, hi0 = mod.hi;
    while (mod.lo > o.anchor_bp && !unique(mod.lo - o.anchor_bp, mod.lo)) {
        if (lo0 - mod.lo >= o.max_anchor_shift) { rep.reason = "no anchor unique on the reference before the module"; return; }
        mod.lo -= std::min<std::size_t>(1000, mod.lo);
    }
    while (mod.hi + o.anchor_bp < ref_seq.size() && !unique(mod.hi, mod.hi + o.anchor_bp)) {
        if (mod.hi - hi0 >= o.max_anchor_shift) { rep.reason = "no anchor unique on the reference after the module"; return; }
        mod.hi += 1000;
    }

    std::string L, R;
    if (!cut_nodes(g, ref, mod.lo, mod.hi, next, copy_nodes, L, R)) {
        rep.reason = "no node the reference walks once, and no path twice, on one side of the module";
        return;
    }
    {
        std::size_t at = 0;
        for (const Step& st : g.walks[ref]) {
            if (st.first == L) rep.ref_start = at + g.seq[L].size();
            if (st.first == R) rep.ref_end = at;
            at += g.seq[st.first].size();
        }
    }

    // Spanning paths -- those walking both cut nodes, forward once read in the reference's orientation -- split
    // into the walk before, the stretch between, and the walk after. Others keep their walks.
    const std::size_t n = g.walks.size();
    std::vector<Walk> pre(n), post(n);
    std::vector<std::string> region(n), expect(n);
    std::vector<char> spans(n, 0), flipped(n, 0), exact(n, 0);
    for (std::size_t p = 0; p < n; ++p) {
        Walk w = g.walks[p];
        auto pos = [&](const std::string& id) {
            for (std::size_t k = 0; k < w.size(); ++k) if (w[k].first == id) return k;
            return w.size();
        };
        std::size_t i = pos(L), j = pos(R);
        if (i == w.size() || j == w.size()) continue;
        if (i > j) { w = flip(w); flipped[p] = 1; i = w.size() - 1 - i; j = w.size() - 1 - j; }
        if (w[i].second || w[j].second) continue;
        // A cut that takes in an earlier module's copies would rewrite them as this module's flank.
        for (std::size_t k = i + 1; k < j; ++k)
            if (copy_nodes.count(w[k].first)) {
                rep.ref_start = cand_lo; rep.ref_end = cand_hi;   // the candidate, as for any overlap
                rep.reason = "overlaps a module already folded";
                return;
            }
        spans[p] = 1;
        pre[p].assign(w.begin(), w.begin() + static_cast<long>(i) + 1);
        post[p].assign(w.begin() + static_cast<long>(j), w.end());
        region[p] = spell(g, Walk(w.begin() + static_cast<long>(i) + 1, w.begin() + static_cast<long>(j)));
    }
    if (!spans[ref]) { rep.reason = "the reference does not walk both cut nodes"; return; }

    // Each spanning path's stretch from its own sequence. The walk's last 2 kb before the stretch and first 2 kb
    // after it are aligned to the haplotype to find it and its strand; each junction is then placed exactly by an
    // edit-distance alignment of the walk's bases next to it (below), so a mismatch or an indel near the junction
    // neither moves it nor rules it out. Ambiguous or poor placements leave the path its rebuilt stretch.
    run_parallel(n, o.threads, [&](std::size_t p) {
        if (!spans[p] || own[p].empty()) return;
        const std::string a = spell(g, pre[p]), b = spell(g, post[p]);
        const std::string lp = a.substr(a.size() > 2000 ? a.size() - 2000 : 0), rp = b.substr(0, 2000);
        if (lp.size() < 200 || rp.size() < 200) return;   // shorter, the probe cannot be placed
        auto one = [&](const std::string& probe) -> std::optional<Minimap2Hit> {
            std::vector<Minimap2Hit> keep;
            for (const Minimap2Hit& h : minimap2_hits("p", probe, "h", own[p], params("asm5", 8, 0.1f)))
                if (10 * (h.query_end_bp - h.query_start_bp) >= 9 * probe.size()) keep.push_back(h);
            if (keep.size() != 1) return std::nullopt;
            return keep[0];
        };
        const std::optional<Minimap2Hit> l = one(lp), r = one(rp);
        if (!l || !r || l->reverse != r->reverse) return;
        const std::size_t len = own[p].size();
        const std::string oriented = l->reverse ? reverse_complement(own[p]) : own[p];
        // Estimates on the oriented sequence: where the walk's last base, and next first base, fall.
        const long long le_est = l->reverse ? static_cast<long long>(len - l->target_start_bp) + static_cast<long long>(lp.size() - l->query_end_bp)
                                            : static_cast<long long>(l->target_end_bp + (lp.size() - l->query_end_bp));
        const long long rs_est = r->reverse ? static_cast<long long>(len - r->target_end_bp) - static_cast<long long>(r->query_start_bp)
                                            : static_cast<long long>(r->target_start_bp) - static_cast<long long>(r->query_start_bp);
        // Place `q` in oriented[est - before, est + after): the end (want_end) or the start of its best alignment
        // nearest the estimate, or -1. `ambiguous` says whether it fits more than one place equally well.
        auto place = [&](const std::string& q, long long est, long long before, long long after, bool want_end,
                         bool& ambiguous) -> long long {
            ambiguous = false;
            const long long w0 = std::max<long long>(0, est - before), w1 = std::min<long long>(static_cast<long long>(len), est + after);
            if (w1 - w0 < static_cast<long long>(q.size())) return -1;
            const EdlibAlignResult res = edlibAlign(q.data(), static_cast<int>(q.size()), oriented.data() + w0,
                                                    static_cast<int>(w1 - w0),
                                                    edlibNewAlignConfig(static_cast<int>(q.size() / 10), EDLIB_MODE_HW, EDLIB_TASK_LOC, nullptr, 0));
            long long best = -1;
            if (res.status == EDLIB_STATUS_OK && res.editDistance >= 0)
                for (int i = 0; i < res.numLocations; ++i) {
                    const long long at = w0 + (want_end ? res.endLocations[i] + 1 : res.startLocations[i]);
                    if (best >= 0 && at != best) ambiguous = true;
                    if (best < 0 || std::llabs(at - est) < std::llabs(best - est)) best = at;
                }
            edlibFreeAlignResult(res);
            return best;
        };
        // The walk's bases beside the junction: 300, or the whole flank when shorter, doubled up to the flank while
        // they fit more than one place equally well (a repeat beside the junction); then the place nearest the
        // estimate.
        auto junction = [&](const std::string& flank, bool walk_end, long long est) -> long long {
            std::size_t n = std::min<std::size_t>(300, flank.size());
            for (;;) {
                const std::string q = walk_end ? flank.substr(flank.size() - n) : flank.substr(0, n);
                const long long w = static_cast<long long>(n);
                bool ambiguous = false;
                const long long at = place(q, est, walk_end ? 2 * w : w, walk_end ? w : 2 * w, walk_end, ambiguous);
                if (!ambiguous || n == flank.size()) return at;
                n = std::min(flank.size(), 2 * n);
            }
        };
        const long long le = junction(lp, true, le_est);
        const long long rs = junction(rp, false, rs_est);
        if (le < 0 || rs < 0 || le > rs) return;
        region[p] = oriented.substr(static_cast<std::size_t>(le), static_cast<std::size_t>(rs - le));
        exact[p] = 1;
    });
    for (std::size_t p = 0; p < n; ++p) {
        if (!spans[p]) continue;
        const std::string full = spell(g, pre[p]) + region[p] + spell(g, post[p]);
        expect[p] = flipped[p] ? reverse_complement(full) : full;
    }

    // Copies: forward occurrences covering min_unit_cover of the unit at the module's copy identity. That
    // identity comes from the reference's own copies when it holds two or more (their lowest, less two points),
    // within [min_copy_identity, min_identity]: paralogs the reference holds as copies count wherever they occur.
    std::vector<std::vector<Minimap2Hit>> hits(n);
    run_parallel(n, o.threads, [&](std::size_t p) {
        if (spans[p] && !region[p].empty()) hits[p] = minimap2_hits("r", region[p], "U", mod.unit, params("map-hifi", 50, 0.1f));
    });
    const std::size_t back = std::min<std::size_t>(200, mod.unit.size() / 20);
    double copy_identity = o.min_identity;
    {
        // Lowered only by copies the reference holds side by side (paralogs): one copy against a unit from another
        // haplotype says nothing about how far copies may diverge.
        double lowest = 1.0;
        std::size_t held = 0;
        for (const Occurrence& oc : occurrences(hits[ref], mod.unit.size(), back))
            if (oc.fwd && oc.cover >= o.min_unit_cover && oc.gcid >= o.min_copy_identity) { lowest = std::min(lowest, oc.gcid); ++held; }
        if (held >= 2) copy_identity = std::clamp(lowest - 0.02, o.min_copy_identity, o.min_identity);
    }
    rep.copy_identity = copy_identity;
    std::vector<std::vector<Piece>> plan(n);
    std::vector<std::size_t> copies(n, 0);
    for (std::size_t p = 0; p < n; ++p) {
        if (!spans[p]) continue;
        std::size_t at = 0;
        for (const Occurrence& oc : occurrences(hits[p], mod.unit.size(), back)) {
            if (!oc.fwd || oc.cover < o.min_unit_cover || oc.gcid < copy_identity || oc.qs < at) continue;
            if (oc.qs > at) plan[p].push_back({copies[p] == 0 ? 1 : 2, {}, region[p].substr(at, oc.qs - at)});
            plan[p].push_back({0, aligned_items(oc.pieces, region[p], mod.unit), {}});
            at = oc.qe;
            ++copies[p];
        }
        if (at < region[p].size()) plan[p].push_back({copies[p] == 0 ? 1 : 3, {}, region[p].substr(at)});
    }
    std::size_t max_copies = 0;
    for (std::size_t p = 0; p < n; ++p) if (spans[p]) max_copies = std::max(max_copies, copies[p]);
    if (max_copies < 2) { rep.reason = "no haplotype holds two copies"; return; }

    // Sequence around the copies, and long differences inside them, go on templates and share their segments:
    // as one node per distinct sequence they made the unit's ends hubs. Templates are made first-fit, the
    // reference's pieces first: a piece joins the first template of its kind it aligns to (forward, over half
    // its length, at 0.9 or more) or becomes one. Long differences get templates per position.
    std::vector<std::string> tmpl{mod.unit};
    std::vector<int> tmpl_kind{0};   // 0 unit; 1 before the copies; 2 between; 3 after; 4+ a difference's position
    std::map<std::pair<int, std::string>, std::vector<Item>> done;
    auto assign = [&](int kind, const std::string& q) -> std::vector<Item> {
        if (auto it = done.find({kind, q}); it != done.end()) return it->second;
        std::vector<Item> items;
        if (q.size() < o.template_min_bp) {
            items.push_back({'L', 0, 0, q});
            return done[{kind, q}] = items;
        }
        int chosen = -1;
        for (std::size_t t = 1; t < tmpl.size() && chosen < 0; ++t) {
            if (tmpl_kind[t] != kind) continue;
            if (tmpl[t] == q) { chosen = static_cast<int>(t); items.push_back({'U', 0, q.size(), {}}); break; }
            const std::vector<Minimap2Hit> hs = minimap2_hits("q", q, "t", tmpl[t], params("map-hifi", 1, 0.1f));
            const Minimap2Hit* best = nullptr;
            for (const Minimap2Hit& hh : hs)
                if (!hh.reverse && hh.gc_identity >= 0.9 && 2 * (hh.query_end_bp - hh.query_start_bp) >= q.size() &&
                    (best == nullptr || hh.query_end_bp - hh.query_start_bp > best->query_end_bp - best->query_start_bp))
                    best = &hh;
            if (best == nullptr) continue;
            chosen = static_cast<int>(t);
            if (best->query_start_bp > 0) items.push_back({'L', 0, 0, q.substr(0, best->query_start_bp)});
            for (Item& it : aligned_items({best}, q, tmpl[t])) items.push_back(std::move(it));
            if (best->query_end_bp < q.size()) items.push_back({'L', 0, 0, q.substr(best->query_end_bp)});
        }
        if (chosen < 0) {
            tmpl.push_back(q); tmpl_kind.push_back(kind);
            chosen = static_cast<int>(tmpl.size() - 1);
            items.push_back({'U', 0, q.size(), {}});
        }
        for (Item& it : items) if (it.kind != 'L') it.t = chosen;
        return done[{kind, q}] = items;
    };
    {
        std::vector<std::pair<std::size_t, std::size_t>> order;   // (path, piece), the reference's first
        for (std::size_t k = 0; k < plan[ref].size(); ++k) order.emplace_back(ref, k);
        for (std::size_t p = 0; p < n; ++p)
            if (p != ref) for (std::size_t k = 0; k < plan[p].size(); ++k) order.emplace_back(p, k);
        std::map<std::tuple<char, std::size_t, std::size_t>, int> site_kind;
        for (const auto& [p, k] : order) {
            Piece& pc = plan[p][k];
            if (pc.tmpl != 0) { pc.items = assign(pc.tmpl, pc.lit); pc.tmpl = 0; continue; }
            std::vector<Item> out;
            for (Item& it : pc.items) {
                if ((it.kind == 'V' || it.kind == 'J') && it.seq.size() >= o.template_min_bp) {
                    auto sk = site_kind.find({it.kind, it.a, it.b});
                    if (sk == site_kind.end())
                        sk = site_kind.emplace(std::make_tuple(it.kind, it.a, it.b), 4 + static_cast<int>(site_kind.size())).first;
                    out.push_back({'D', it.a, it.b, {}});
                    for (Item& sub : assign(sk->second, it.seq)) out.push_back(std::move(sub));
                } else {
                    out.push_back(std::move(it));
                }
            }
            pc.items = std::move(out);
        }
    }

    // Each template cut wherever an item on it starts, stops or differs; then every spanning path's walk.
    std::vector<std::set<std::size_t>> cuts(tmpl.size());
    for (std::size_t t = 0; t < tmpl.size(); ++t) cuts[t] = {0, tmpl[t].size()};
    for (const auto& pl : plan)
        for (const Piece& pc : pl)
            for (const Item& it : pc.items)
                if (it.kind != 'L') { auto& c = cuts[static_cast<std::size_t>(it.t >= 0 ? it.t : 0)]; c.insert(it.a); c.insert(it.b); }
    FoldGraph h = g;
    std::vector<std::map<std::size_t, std::string>> seg(tmpl.size());
    for (std::size_t t = 0; t < tmpl.size(); ++t) {
        const char* kind = tmpl_kind[t] == 0 ? "fold_unit" : tmpl_kind[t] == 2 ? "fold_spacer"
                         : tmpl_kind[t] <= 3 ? "fold_flank" : "fold_insertion";
        const std::vector<std::size_t> c(cuts[t].begin(), cuts[t].end());
        for (std::size_t k = 0; k + 1 < c.size(); ++k) {
            const std::string id = fresh_id(h, next);
            h.seq[id] = tmpl[t].substr(c[k], c[k + 1] - c[k]);
            h.tags[id] = std::string("\tmn:Z:") + kind + (t == 0 ? std::string() : ":" + std::to_string(t)) + ":" +
                         std::to_string(c[k]) + "-" + std::to_string(c[k + 1]);
            seg[t][c[k]] = id;
        }
    }
    std::map<std::tuple<int, char, std::size_t, std::size_t, std::string>, std::string> alt;
    std::map<std::string, std::string> lit;
    auto node_of = [&](auto& m, const auto& key, const std::string& s, const char* tag) {
        auto f = m.find(key);
        if (f == m.end()) {
            const std::string id = fresh_id(h, next);
            h.seq[id] = s; h.tags[id] = std::string("\tmn:Z:") + tag;
            f = m.emplace(key, id).first;
        }
        return f->second;
    };
    for (std::size_t p = 0; p < n; ++p) {
        if (!spans[p]) continue;
        Walk w = pre[p];
        for (const Piece& pc : plan[p])
            for (const Item& it : pc.items) {
                const std::size_t on = static_cast<std::size_t>(it.t >= 0 ? it.t : 0);
                if (it.kind == 'D') continue;
                if (it.kind == 'U') {
                    for (std::size_t x = it.a; x < it.b; x += h.seq[seg[on].at(x)].size()) w.emplace_back(seg[on].at(x), false);
                } else if (it.kind == 'L') {
                    if (!it.seq.empty()) w.emplace_back(node_of(lit, it.seq, it.seq, "fold_literal"), false);
                } else if (!it.seq.empty()) {
                    w.emplace_back(node_of(alt, std::make_tuple(static_cast<int>(on), it.kind, it.a, it.b, it.seq), it.seq,
                                           "fold_difference"), false);
                }
            }
        w.insert(w.end(), post[p].begin(), post[p].end());
        h.walks[p] = flipped[p] ? flip(w) : w;
    }

    // A path that would not spell as it should keeps its walk; the reference must fold.
    for (std::size_t p = 0; p < n; ++p) {
        if (!spans[p]) continue;
        if (spell(h, h.walks[p]) == expect[p]) continue;
        if (p == ref) { rep.reason = "the reference would not spell as it should"; return; }
        h.walks[p] = g.walks[p];
        spans[p] = 0;
    }
    rep.reference_copies = copies[ref];
    for (std::size_t p = 0; p < n; ++p) if (spans[p]) rep.folded_paths.push_back(g.names[p]);
    rep.min_copies = SIZE_MAX;
    for (std::size_t p = 0; p < n; ++p) {
        if (!spans[p]) continue;
        ++rep.paths_folded;
        rep.paths_exact += exact[p];
        rep.min_copies = std::min(rep.min_copies, copies[p]);
        rep.max_copies = std::max(rep.max_copies, copies[p]);
    }
    for (std::size_t p = 0; p < n; ++p) rep.paths_spanning += !pre[p].empty();
    for (std::size_t t = 1; t < tmpl.size(); ++t) (tmpl_kind[t] <= 3 ? rep.stretch_templates : rep.difference_templates)++;
    std::unordered_map<std::string, std::unordered_set<std::string>> nb;
    for (const Walk& w : h.walks)
        for (std::size_t k = 0; k + 1 < w.size(); ++k) { nb[w[k].first].insert(w[k + 1].first); nb[w[k + 1].first].insert(w[k].first); }
    for (const auto& [id, s] : nb) if (!g.seq.count(id) && s.size() >= o.hub_degree) ++rep.hubs;
    // This module's copy nodes: the unit's, the templates of long differences, and the difference nodes on
    // either. Those on the flank and spacer templates are not copies, and a later cut may take them in.
    for (std::size_t t = 0; t < tmpl.size(); ++t)
        if (t == 0 || tmpl_kind[t] >= 4) for (const auto& [x, id] : seg[t]) copy_nodes.insert(id);
    for (const auto& [key, id] : alt) {
        const std::size_t t = static_cast<std::size_t>(std::get<0>(key));
        if (t == 0 || tmpl_kind[t] >= 4) copy_nodes.insert(id);
    }
    folded = std::move(h);
    rep.folded = true;
}

} // namespace

std::vector<ModuleUnit> discover_modules(const std::vector<std::string>& names, const std::vector<std::string>& seqs,
                                         std::size_t ref, const ModuleFoldOptions& options) {
    return discover(names, seqs, ref, options);
}

ModuleFoldResult fold_modules(const std::string& gfa,
                              const std::vector<std::pair<std::string, std::string>>& haplotypes,
                              const std::string& reference,
                              const ModuleFoldOptions& options,
                              const std::vector<ModuleUnit>* known) {
    ModuleFoldResult result;
    result.gfa = gfa;
    FoldGraph g = parse(gfa);
    const auto rit = std::find(g.names.begin(), g.names.end(), reference);
    if (rit == g.names.end()) return result;
    const std::size_t ref = static_cast<std::size_t>(rit - g.names.begin());
    std::unordered_map<std::string, const std::string*> by_name;
    for (const auto& [nm, s] : haplotypes) by_name[nm] = &s;
    std::vector<std::string> own(g.names.size());
    for (std::size_t p = 0; p < g.names.size(); ++p)
        if (auto it = by_name.find(g.names[p]); it != by_name.end()) own[p] = *it->second;
    // Modules are found from the haplotypes' own sequences, which the rebuilt walks spell only within --min-var.
    std::vector<Module> modules;
    if (known != nullptr) {
        modules = *known;
    } else {
        std::vector<std::string> seqs(g.walks.size());
        for (std::size_t p = 0; p < g.walks.size(); ++p) seqs[p] = own[p].empty() ? spell(g, g.walks[p]) : own[p];
        modules = discover(g.names, seqs, ref, options);
    }
    std::uint64_t next = 1;
    for (const std::string& id : g.order) {
        char* end = nullptr;
        const unsigned long long v = std::strtoull(id.c_str(), &end, 10);
        if (end != nullptr && *end == '\0' && v >= next) next = v + 1;
    }
    bool any = false;
    // The copies' stretches of the modules already folded. A module's cut may run past its copies into a
    // neighbouring array's when no node between them can be cut; that array is still tried, and refused only
    // if its own cut would take in folded copies.
    std::vector<std::pair<std::size_t, std::size_t>> taken;
    std::unordered_set<std::string> copy_nodes;
    for (const Module& m : modules) {
        FoldedModule rep;
        if (std::any_of(taken.begin(), taken.end(), [&](const auto& t) { return m.lo < t.second && t.first < m.hi; })) {
            rep.ref_start = m.lo; rep.ref_end = m.hi;   // the candidate, so the line says which one was skipped
            rep.period = m.period; rep.unit_bp = m.unit.size(); rep.unit_source = m.source;
            rep.reason = "overlaps a module already folded";
        } else {
            try {
                fold_one(g, own, ref, m, options, next, copy_nodes, rep);
            } catch (const std::exception& e) {
                rep.folded = false;
                rep.reason = std::string("error: ") + e.what();
            }
        }
        if (rep.folded) {
            any = true;
            taken.emplace_back(m.lo, m.hi);
            // This cut may have taken in an earlier module's flank: that stretch now belongs to this module.
            for (FoldedModule& e : result.modules) {
                if (!e.folded || e.ref_end <= rep.ref_start || rep.ref_end <= e.ref_start) continue;
                if (e.ref_start < rep.ref_start) e.ref_end = rep.ref_start;
                else if (e.ref_end > rep.ref_end) e.ref_start = rep.ref_end;
            }
        }
        result.modules.push_back(std::move(rep));
    }
    if (!any) return result;

    // Emit the segments any walk uses, in input order then first use, and each walked link once.
    std::unordered_set<std::string> used;
    for (const Walk& w : g.walks) for (const Step& st : w) used.insert(st.first);
    std::vector<std::string> ids;
    std::unordered_set<std::string> listed;
    for (const std::string& id : g.order) if (used.count(id) && listed.insert(id).second) ids.push_back(id);
    for (const Walk& w : g.walks) for (const Step& st : w) if (listed.insert(st.first).second) ids.push_back(st.first);
    using Link = std::tuple<std::string, bool, std::string, bool>;
    std::set<Link> links;
    for (const Walk& w : g.walks)
        for (std::size_t k = 0; k + 1 < w.size(); ++k) {
            const Link fwd{w[k].first, w[k].second, w[k + 1].first, w[k + 1].second};
            const Link rev{w[k + 1].first, !w[k + 1].second, w[k].first, !w[k].second};
            links.insert(std::min(fwd, rev));
        }
    std::ostringstream out;
    out << "H\tVN:Z:1.0\n";
    for (const std::string& id : ids) out << "S\t" << id << '\t' << (g.seq[id].empty() ? "*" : g.seq[id]) << g.tags[id] << '\n';
    for (const auto& [a, ra, b, rb] : links) out << "L\t" << a << '\t' << (ra ? '-' : '+') << '\t' << b << '\t' << (rb ? '-' : '+') << "\t0M\n";
    for (std::size_t p = 0; p < g.walks.size(); ++p) {
        out << "P\t" << g.names[p] << '\t';
        for (std::size_t k = 0; k < g.walks[p].size(); ++k)
            out << (k ? "," : "") << g.walks[p][k].first << (g.walks[p][k].second ? '-' : '+');
        out << "\t*\n";
    }
    result.gfa = out.str();
    result.nodes = ids.size();
    result.links = links.size();
    // Degree per segment end, as the gate measures it. Self-loops count nodes with any link to themselves.
    std::unordered_map<std::string, std::array<std::unordered_set<std::string>, 2>> ends;
    std::unordered_set<std::string> looped;
    for (const auto& [a, ra, b, rb] : links) {
        if (a == b) { looped.insert(a); continue; }
        ends[a][ra ? 0 : 1].insert(b);
        ends[b][rb ? 1 : 0].insert(a);
    }
    for (const auto& [id, e] : ends) {
        const std::size_t d = std::max(e[0].size(), e[1].size());
        result.maxdeg = std::max(result.maxdeg, d);
        if (d >= options.gate_hub_degree) ++result.hubs;
    }
    result.selfloops = looped.size();
    for (std::size_t p = 0; p < g.walks.size(); ++p)
        result.spells_own.emplace_back(g.names[p], !own[p].empty() && spell(g, g.walks[p]) == own[p]);
    return result;
}

} // namespace panvar
