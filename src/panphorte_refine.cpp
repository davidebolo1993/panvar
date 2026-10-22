#include "panvar/panphorte_refine.hpp"

#include "panvar/cli_utils.hpp"
#include "panvar/gfa.hpp"
#include "panvar/gfa_io.hpp"
#include "panvar/graph_utils.hpp"
#include "panvar/poa.hpp"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <map>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace panvar {
namespace {

using Steps = std::vector<PathStep>;

std::string spell(const GfaModel& model, const Steps& steps) {
    std::string out;
    for (const auto& s : steps) {
        auto it = model.seq.find(s.node_id);
        if (it == model.seq.end()) continue;
        out += s.reverse ? reverse_complement(it->second) : it->second;
    }
    return out;
}

char orient(bool reverse) { return reverse ? '-' : '+'; }

// One arc read from either end (a+>b+ is b->a-), so a link stored from one side is found from the other.
std::string arc_key(const std::string& f, char fo, const std::string& t, char to) {
    const std::string fwd = f + fo + ">" + t + to;
    const std::string rev = t + (to == '+' ? '-' : '+') + ">" + f + (fo == '+' ? '-' : '+');
    return fwd <= rev ? fwd : rev;
}

// The interior a path carries between the two anchors, as in refine: every (a, b) occurrence pair is
// scored by how much of the site's interior it contains, the tightest span wins among equals, an exact
// tie is ambiguous, and a reverse crossing is canonicalized (lo/hi keep the original indices).
struct Interior {
    Steps steps;
    std::size_t lo = 0;
    std::size_t hi = 0;
    bool reversed = false;
};
std::optional<Interior> interior_between(const Steps& steps, const std::string& a, const std::string& b,
                                         const std::unordered_set<std::string>& region_interior,
                                         bool* ambiguous) {
    if (ambiguous != nullptr) *ambiguous = false;
    std::vector<std::size_t> at_a, at_b;
    for (std::size_t i = 0; i < steps.size(); ++i) {
        if (steps[i].node_id == a) at_a.push_back(i);
        if (steps[i].node_id == b) at_b.push_back(i);
    }
    if (at_a.empty() || at_b.empty()) return std::nullopt;
    struct Cand { std::size_t lo, hi, score, span; bool reversed; };
    std::vector<Cand> cands;
    const auto consider = [&](std::size_t lo, std::size_t hi, bool reversed) {
        if (hi <= lo) return;
        std::size_t score = 0;
        for (std::size_t k = lo + 1; k < hi; ++k)
            if (region_interior.count(steps[k].node_id) != 0) ++score;
        cands.push_back({lo, hi, score, hi - lo, reversed});
    };
    for (const std::size_t ia : at_a)
        for (const std::size_t ib : at_b) {
            if (ib > ia) consider(ia, ib, false);
            else if (ia > ib) consider(ib, ia, true);
        }
    if (cands.empty()) return std::nullopt;
    std::sort(cands.begin(), cands.end(), [](const Cand& x, const Cand& y) {
        if (x.score != y.score) return x.score > y.score;
        return x.span < y.span;
    });
    if (cands.size() > 1 && cands[0].score == cands[1].score && cands[0].span == cands[1].span) {
        if (ambiguous != nullptr) *ambiguous = true;
        return std::nullopt;
    }
    const Cand& best = cands.front();
    Interior out;
    out.lo = best.lo;
    out.hi = best.hi;
    out.reversed = best.reversed;
    out.steps.assign(steps.begin() + static_cast<long>(best.lo) + 1, steps.begin() + static_cast<long>(best.hi));
    if (best.reversed) {
        std::reverse(out.steps.begin(), out.steps.end());
        for (PathStep& st : out.steps) st.reverse = !st.reverse;
    }
    return out;
}

// A path's interior split around its ARRAY BLOCK: every step from the first REP step to the last -- the
// folded copies and whatever folding left between them. The block is kept verbatim and never realigned:
// a fragment repeated between copies (REP GT REP GT REP) belongs to the folded array, not to an unfolded
// duplication. What is realigned are the flanks before and after it. A path with no REP step is a
// non-carrier and is left as folded.
struct Shape {
    bool carrier = false;
    std::size_t block_lo = 0, block_hi = 0;   // [lo, hi) in the interior
    Steps pre, block, post;
};
Shape shape_of(const Steps& interior, const std::unordered_set<std::string>& reps) {
    Shape s;
    std::size_t first = interior.size(), last = 0;
    for (std::size_t i = 0; i < interior.size(); ++i)
        if (reps.count(interior[i].node_id) != 0) { first = std::min(first, i); last = i; }
    if (first == interior.size()) {
        s.block_lo = s.block_hi = interior.size();
        return s;
    }
    s.carrier = true;
    s.block_lo = first;
    s.block_hi = last + 1;
    s.pre.assign(interior.begin(), interior.begin() + static_cast<long>(first));
    s.block.assign(interior.begin() + static_cast<long>(first), interior.begin() + static_cast<long>(last) + 1);
    s.post.assign(interior.begin() + static_cast<long>(last) + 1, interior.end());
    return s;
}

// refine's MSA compaction: a run of columns traversed by exactly the same rows becomes one node.
struct InteriorGraph {
    std::vector<std::string> node_seq;               // local node index -> sequence
    std::vector<std::vector<std::size_t>> row_paths;  // per row: local node indices
};
InteriorGraph build_interior_graph(const std::vector<std::string>& rows) {
    using Cell = std::pair<int, char>;
    const std::size_t S = rows.size();
    const int L = rows.empty() ? 0 : static_cast<int>(rows[0].size());
    std::vector<std::vector<Cell>> row_cells(S);
    std::map<Cell, std::set<std::size_t>> users;
    for (std::size_t s = 0; s < S; ++s) {
        for (int c = 0; c < L; ++c) {
            const char ch = rows[s][static_cast<std::size_t>(c)];
            if (ch == '-') continue;
            Cell cell{c, ch};
            row_cells[s].push_back(cell);
            users[cell].insert(s);
        }
    }
    std::map<Cell, std::set<Cell>> out_edges, in_edges;
    std::map<std::pair<Cell, Cell>, std::set<std::size_t>> edge_users;
    for (std::size_t s = 0; s < S; ++s) {
        const auto& cells = row_cells[s];
        for (std::size_t k = 1; k < cells.size(); ++k) {
            out_edges[cells[k - 1]].insert(cells[k]);
            in_edges[cells[k]].insert(cells[k - 1]);
            edge_users[{cells[k - 1], cells[k]}].insert(s);
        }
    }
    std::map<Cell, Cell> merged;
    std::function<Cell(Cell)> find = [&](Cell x) -> Cell {
        auto it = merged.find(x);
        while (it != merged.end() && !(it->second == x)) {
            x = it->second;
            it = merged.find(x);
        }
        return x;
    };
    for (const auto& kv : users) {
        const Cell& u = kv.first;
        auto oi = out_edges.find(u);
        if (oi == out_edges.end() || oi->second.size() != 1) continue;
        const Cell v = *oi->second.begin();
        auto ii = in_edges.find(v);
        if (ii == in_edges.end() || ii->second.size() != 1) continue;
        if (users.at(u) == users.at(v) && users.at(u) == edge_users.at({u, v})) merged[v] = find(u);
    }
    std::map<Cell, std::vector<Cell>> group_cells;
    for (const auto& kv : users) group_cells[find(kv.first)].push_back(kv.first);
    InteriorGraph out;
    std::map<Cell, std::size_t> gid;
    for (auto& kv : group_cells) {
        auto& cells = kv.second;
        std::sort(cells.begin(), cells.end());
        std::string seq;
        for (const auto& c : cells) seq += c.second;
        gid[kv.first] = out.node_seq.size();
        out.node_seq.push_back(seq);
    }
    out.row_paths.resize(S);
    for (std::size_t s = 0; s < S; ++s) {
        std::vector<std::size_t> p;
        for (const auto& cell : row_cells[s]) {
            const std::size_t g = gid[find(cell)];
            if (p.empty() || p.back() != g) p.push_back(g);
        }
        out.row_paths[s] = std::move(p);
    }
    return out;
}

// One accepted site: the old flank nodes it retires, the rebuilt nodes (local index -> sequence), and
// each carrier's replacement interior.
struct SiteEdit {
    std::size_t row_index = 0;
    std::unordered_set<std::string> retire;
    std::vector<std::string> new_seq;
    // path -> (lo, hi, interior steps; rebuilt flank nodes carry "#<local index>" ids until assigned)
    std::unordered_map<std::size_t, std::tuple<std::size_t, std::size_t, Steps>> per_path;
    std::unordered_map<std::size_t, std::string> group_key;   // path -> identical-haplotype key
    std::unordered_set<std::size_t> changed;                  // paths with at least one flank step
};

using RepSignature = std::vector<std::tuple<std::string, char, std::size_t>>;
RepSignature rep_signature(const Steps& steps, const std::unordered_set<std::string>& reps) {
    RepSignature sig;
    for (const PathStep& st : steps) {
        if (reps.count(st.node_id) == 0) continue;
        if (!sig.empty() && std::get<0>(sig.back()) == st.node_id && std::get<1>(sig.back()) == orient(st.reverse))
            ++std::get<2>(sig.back());
        else
            sig.emplace_back(st.node_id, orient(st.reverse), 1);
    }
    return sig;
}

std::string steps_key(const Steps& steps) {
    std::string k;
    for (const PathStep& st : steps) k += st.node_id + orient(st.reverse) + ",";
    return k;
}

// Non-REP nodes and steps outside the array blocks, over every traversing path. `routes` replaces a
// carrier's flank steps by its rebuilt flank route when given.
struct Measure { std::size_t nodes = 0, steps = 0, single = 0; };
Measure measure(const std::unordered_map<std::size_t, Interior>& per_path,
                const std::unordered_map<std::size_t, Shape>& shapes,
                const std::unordered_map<std::size_t, std::vector<std::size_t>>* routes) {
    Measure m;
    std::unordered_map<std::string, std::size_t> users;
    for (const auto& [idx, in] : per_path) {
        std::unordered_set<std::string> mine;
        const std::vector<std::size_t>* route = nullptr;
        if (routes != nullptr) {
            const auto rit = routes->find(idx);
            if (rit != routes->end()) route = &rit->second;
        }
        if (route != nullptr) {
            for (const std::size_t n : *route) { ++m.steps; mine.insert("#" + std::to_string(n)); }
        } else {
            const Shape& s = shapes.at(idx);
            for (std::size_t k = 0; k < in.steps.size(); ++k)
                if (k < s.block_lo || k >= s.block_hi) { ++m.steps; mine.insert(in.steps[k].node_id); }
        }
        for (const auto& n : mine) ++users[n];
    }
    m.nodes = users.size();
    for (const auto& kv : users) if (kv.second == 1) ++m.single;
    return m;
}

std::optional<SiteEdit> process_site(const PanphorteRefineSite& site, const GfaModel& model,
                                     const std::vector<std::size_t>& path_order, std::size_t ref_idx,
                                     const std::unordered_set<std::string>& self_loop,
                                     const PanphorteRefineLimits& lim, PanphorteRefineRow& row) {
    const auto skip = [&](const std::string& reason, const std::string& detail) {
        row.status = "skipped:" + reason;
        row.detail = detail;
        return std::nullopt;
    };
    if (site.source.empty() || site.sink.empty() || site.source == site.sink ||
        site.interior.count(site.source) != 0 || site.interior.count(site.sink) != 0)
        return skip("not_a_clean_chain", "anchors do not reduce to one source and one sink");

    std::unordered_set<std::string> reps;
    for (const auto& n : site.interior)
        if (self_loop.count(n) != 0) reps.insert(n);

    std::unordered_map<std::size_t, Interior> per_path;
    std::unordered_map<std::size_t, Shape> shapes;
    std::size_t partial = 0;
    for (const std::size_t idx : path_order) {
        const Steps& steps = model.paths[idx].steps;
        bool ambiguous = false;
        auto got = interior_between(steps, site.source, site.sink, site.interior, &ambiguous);
        if (ambiguous)
            return skip("ambiguous_anchor_traversal",
                        "path " + gfa_path_name(model.paths[idx]) + " has two equally good traversals");
        // Interior nodes the chosen traversal does not cover: a partial or second traversal, which keeps
        // its old nodes and so keeps the old route walkable beside the rebuilt one.
        std::size_t outside = 0;
        for (std::size_t k = 0; k < steps.size(); ++k) {
            if (site.interior.count(steps[k].node_id) == 0) continue;
            if (got && k > got->lo && k < got->hi) continue;
            ++outside;
        }
        if (outside > 0) { ++partial; continue; }
        if (!got) continue;
        Shape shape = shape_of(got->steps, reps);
        // The duplication guard, on the flanks only: the array block is kept verbatim, so nothing inside
        // it is linearized. A flank node walked twice is an unfolded copy that realignment would flatten.
        std::unordered_map<std::string, int> c;
        for (std::size_t k = 0; k < got->steps.size(); ++k)
            if (k < shape.block_lo || k >= shape.block_hi) ++c[got->steps[k].node_id];
        for (const auto& kv : c)
            if (kv.second >= 2)
                return skip("unfolded_dup", "node " + kv.first + " x" + std::to_string(kv.second) + " on " +
                                                gfa_path_name(model.paths[idx]));
        shapes[idx] = std::move(shape);
        per_path[idx] = std::move(*got);
    }
    if (partial > 0)
        return skip("partial_paths", std::to_string(partial) +
                                         " path(s) touch the interior outside one anchor-to-anchor traversal");
    if (per_path.count(ref_idx) == 0) return skip("reference_not_traversing", "");
    row.paths = per_path.size();

    // Distinct flank sequences, reference first, one alignment for the flanks before the block and one
    // for those after it. Flanks are never aligned across the block: glued together they form a junction
    // no haplotype carries, and its repeat-like edges align ambiguously.
    std::vector<std::string> flank_seqs[2];
    std::unordered_map<std::string, std::size_t> flank_index[2];
    for (const std::size_t idx : path_order) {
        const auto it = shapes.find(idx);
        if (it == shapes.end() || !it->second.carrier) continue;
        ++row.carriers;
        const std::string seqs[2] = {spell(model, it->second.pre), spell(model, it->second.post)};
        for (int side = 0; side < 2; ++side)
            if (!seqs[side].empty() && flank_index[side].emplace(seqs[side], flank_seqs[side].size()).second)
                flank_seqs[side].push_back(seqs[side]);
    }
    row.distinct_pre = flank_seqs[0].size();
    row.distinct_post = flank_seqs[1].size();
    if (row.carriers == 0) return skip("no_carrier", "no traversing path walks a REP node");

    const Measure before = measure(per_path, shapes, nullptr);
    row.nodes_before = before.nodes;
    row.steps_before = before.steps;
    row.single_path_nodes_before = before.single;

    std::size_t longest = 0, total = 0;
    for (int side = 0; side < 2; ++side)
        for (const auto& r : flank_seqs[side]) { longest = std::max(longest, r.size()); total += r.size(); }
    const std::string sizes = "distinct flanks " + std::to_string(row.distinct_pre) + " before / " +
                              std::to_string(row.distinct_post) + " after, longest " + std::to_string(longest) +
                              " bp, total " + std::to_string(total) + " bp";
    if (longest > lim.max_poa_bp) return skip("segment_too_long", sizes);
    if (lim.max_poa_work > 0 && static_cast<double>(longest) * static_cast<double>(total) >
                                    static_cast<double>(lim.max_poa_work))
        return skip("poa_work", sizes);
    if (std::max(row.distinct_pre, row.distinct_post) > lim.max_walks) return skip("too_diverse", sizes);

    SiteEdit edit;
    std::vector<std::vector<std::size_t>> flank_routes[2];
    for (int side = 0; side < 2; ++side) {
        if (flank_seqs[side].empty()) continue;
        InteriorGraph ig = build_interior_graph(poa_msa(flank_seqs[side]));
        const std::size_t base = edit.new_seq.size();
        for (auto& sq : ig.node_seq) edit.new_seq.push_back(std::move(sq));
        for (auto& r : ig.row_paths) {
            for (std::size_t& n : r) n += base;
            flank_routes[side].push_back(std::move(r));
        }
        // The alignment must give every flank back exactly, or the site is left as folded.
        for (std::size_t r = 0; r < flank_seqs[side].size(); ++r) {
            std::string got;
            for (const std::size_t n : flank_routes[side][r]) got += edit.new_seq[n];
            if (got != flank_seqs[side][r])
                return skip("poa_not_lossless", std::string(side == 0 ? "before" : "after") + " flank " +
                                                    std::to_string(r) + " (" +
                                                    std::to_string(flank_seqs[side][r].size()) +
                                                    " bp) not reproduced");
        }
    }

    // each carrier's rebuilt flanks; non-carriers keep theirs
    std::unordered_map<std::size_t, std::vector<std::size_t>> pre_route, post_route, all_routes;
    for (const auto& [idx, shape] : shapes) {
        if (!shape.carrier) continue;
        const std::string seqs[2] = {spell(model, shape.pre), spell(model, shape.post)};
        auto& pr = pre_route[idx];
        auto& po = post_route[idx];
        if (!seqs[0].empty()) pr = flank_routes[0][flank_index[0].at(seqs[0])];
        if (!seqs[1].empty()) po = flank_routes[1][flank_index[1].at(seqs[1])];
        auto& all = all_routes[idx];
        all.insert(all.end(), pr.begin(), pr.end());
        all.insert(all.end(), po.begin(), po.end());
    }

    const Measure after = measure(per_path, shapes, &all_routes);
    row.nodes_after = after.nodes;
    row.steps_after = after.steps;
    row.single_path_nodes_after = after.single;
    row.measured = true;
    const bool simpler = row.nodes_after <= row.nodes_before && row.steps_after <= row.steps_before &&
                         (row.nodes_after < row.nodes_before || row.steps_after < row.steps_before);
    if (!simpler) {
        row.status = "kept_original:not_simpler";
        return std::nullopt;
    }

    for (const auto& [idx, shape] : shapes) {
        if (!shape.carrier) continue;
        const Interior& in = per_path.at(idx);
        Steps rebuilt;
        for (const std::size_t n : pre_route.at(idx)) rebuilt.push_back(PathStep{"#" + std::to_string(n), false});
        rebuilt.insert(rebuilt.end(), shape.block.begin(), shape.block.end());
        for (const std::size_t n : post_route.at(idx)) rebuilt.push_back(PathStep{"#" + std::to_string(n), false});
        if (in.reversed) {
            std::reverse(rebuilt.begin(), rebuilt.end());
            for (PathStep& st : rebuilt) st.reverse = !st.reverse;
        }
        edit.per_path[idx] = std::make_tuple(in.lo, in.hi, std::move(rebuilt));
        // identical haplotypes: the same flanks around the same array block
        edit.group_key[idx] = spell(model, shape.pre) + '\x01' + steps_key(shape.block) + '\x01' +
                              spell(model, shape.post);
        if (!shape.pre.empty() || !shape.post.empty()) edit.changed.insert(idx);
        for (const auto& st : shape.pre) edit.retire.insert(st.node_id);
        for (const auto& st : shape.post) edit.retire.insert(st.node_id);
    }
    row.status = "refined";
    return edit;
}

bool same_steps(const Steps& x, const Steps& y) {
    if (x.size() != y.size()) return false;
    for (std::size_t k = 0; k < x.size(); ++k)
        if (x[k].node_id != y[k].node_id || x[k].reverse != y[k].reverse) return false;
    return true;
}

}  // namespace

PanphorteRefineResult panphorte_refine_sites(GfaModel& model, const std::vector<PanphorteRefineSite>& sites,
                                             const std::string& resolved_reference,
                                             const PanphorteRefineLimits& limits) {
    PanphorteRefineResult result;
    std::unordered_set<std::string> self_loop;
    for (const auto& e : model.edges)
        if (e.from == e.to) self_loop.insert(e.from);

    std::size_t ref_idx = model.paths.size();
    for (std::size_t i = 0; i < model.paths.size(); ++i)
        if (gfa_path_name(model.paths[i]) == resolved_reference) { ref_idx = i; break; }
    if (ref_idx == model.paths.size())
        throw std::runtime_error("panphorte --refine: reference '" + resolved_reference + "' is not a path of the model");
    std::vector<std::size_t> path_order{ref_idx};
    for (std::size_t i = 0; i < model.paths.size(); ++i) if (i != ref_idx) path_order.push_back(i);

    std::vector<std::string> spelled_before;
    std::vector<RepSignature> reps_before;
    for (const GfaPath& p : model.paths) {
        spelled_before.push_back(spell(model, p.steps));
        reps_before.push_back(rep_signature(p.steps, self_loop));
    }

    std::vector<SiteEdit> edits;
    for (const PanphorteRefineSite& site : sites) {
        PanphorteRefineRow row;
        row.bubble_id = site.bubble_id;
        row.source = site.source;
        row.sink = site.sink;
        auto edit = process_site(site, model, path_order, ref_idx, self_loop, limits, row);
        if (edit) {
            edit->row_index = result.rows.size();
            edits.push_back(std::move(*edit));
            ++result.refined;
        } else if (row.status.rfind("kept_original", 0) == 0) {
            ++result.kept;
        } else {
            ++result.skipped;
        }
        result.rows.push_back(std::move(row));
    }
    if (edits.empty()) return result;

    // old interiors with their anchors, kept to prove afterwards that no old route is still walkable
    std::vector<std::tuple<std::size_t, std::size_t, Steps>> old_spans;   // path, edit, old [lo,hi] steps
    for (std::size_t e = 0; e < edits.size(); ++e)
        for (const auto& [idx, ed] : edits[e].per_path) {
            const Steps& st = model.paths[idx].steps;
            old_spans.emplace_back(idx, e, Steps(st.begin() + static_cast<long>(std::get<0>(ed)),
                                                 st.begin() + static_cast<long>(std::get<1>(ed)) + 1));
        }

    // real ids for the rebuilt nodes
    std::uint64_t max_id = 0;
    for (const std::string& id : model.node_order) {
        char* endp = nullptr;
        const std::uint64_t v = std::strtoull(id.c_str(), &endp, 10);
        if (endp != id.c_str() && *endp == '\0' && v > max_id) max_id = v;
    }
    for (SiteEdit& e : edits) {
        std::vector<std::string> ids;
        for (const std::string& s : e.new_seq) {
            std::string id;
            do id = std::to_string(++max_id); while (model.has_node(id));
            model.seq[id] = s;
            model.node_order.push_back(id);
            ids.push_back(id);
            ++result.nodes_added;
        }
        for (auto& [idx, ed] : e.per_path)
            for (PathStep& st : std::get<2>(ed))
                if (!st.node_id.empty() && st.node_id[0] == '#')
                    st.node_id = ids[std::stoul(st.node_id.substr(1))];
    }

    // splice, right to left per path, so earlier indices stay valid
    std::unordered_map<std::size_t, std::vector<std::tuple<std::size_t, std::size_t, const Steps*>>> by_path;
    for (const SiteEdit& e : edits)
        for (const auto& [idx, ed] : e.per_path)
            by_path[idx].emplace_back(std::get<0>(ed), std::get<1>(ed), &std::get<2>(ed));
    for (auto& [idx, list] : by_path) {
        std::sort(list.begin(), list.end(), [](const auto& x, const auto& y) { return std::get<0>(x) > std::get<0>(y); });
        Steps& st = model.paths[idx].steps;
        for (const auto& [lo, hi, ins] : list) {
            Steps next(st.begin(), st.begin() + static_cast<long>(lo) + 1);
            next.insert(next.end(), ins->begin(), ins->end());
            next.insert(next.end(), st.begin() + static_cast<long>(hi), st.end());
            st = std::move(next);
        }
    }

    // Retire old flank nodes no path still walks. A link an old interior walked and no final path walks
    // is dropped too: where old nodes survive (non-carriers still use them), the old route would otherwise
    // stay walkable through them. Then every adjacency the paths now need is added as a 0M link.
    std::unordered_set<std::string> retire;
    for (const SiteEdit& e : edits) retire.insert(e.retire.begin(), e.retire.end());
    for (const GfaPath& p : model.paths)
        for (const PathStep& s : p.steps) retire.erase(s.node_id);
    for (const auto& n : retire) { model.seq.erase(n); ++result.nodes_removed; }
    model.node_order.erase(std::remove_if(model.node_order.begin(), model.node_order.end(),
                                          [&](const std::string& n) { return retire.count(n) != 0; }),
                           model.node_order.end());
    std::unordered_set<std::string> old_arcs, walked;
    for (const auto& span : old_spans) {
        const Steps& old = std::get<2>(span);
        for (std::size_t k = 1; k < old.size(); ++k)
            old_arcs.insert(arc_key(old[k - 1].node_id, orient(old[k - 1].reverse), old[k].node_id, orient(old[k].reverse)));
    }
    for (const GfaPath& p : model.paths)
        for (std::size_t k = 1; k < p.steps.size(); ++k)
            walked.insert(arc_key(p.steps[k - 1].node_id, orient(p.steps[k - 1].reverse), p.steps[k].node_id,
                                  orient(p.steps[k].reverse)));
    std::vector<GfaEdge> kept;
    std::unordered_set<std::string> have;
    for (const GfaEdge& e : model.edges) {
        if (retire.count(e.from) != 0 || retire.count(e.to) != 0) continue;
        const std::string k = arc_key(e.from, e.from_orient, e.to, e.to_orient);
        if (old_arcs.count(k) != 0 && walked.count(k) == 0) { ++result.edges_removed; continue; }
        have.insert(k);
        kept.push_back(e);
    }
    for (const GfaPath& p : model.paths)
        for (std::size_t k = 1; k < p.steps.size(); ++k) {
            const PathStep& a = p.steps[k - 1];
            const PathStep& b = p.steps[k];
            if (have.insert(arc_key(a.node_id, orient(a.reverse), b.node_id, orient(b.reverse))).second) {
                kept.push_back(GfaEdge{a.node_id, orient(a.reverse), b.node_id, orient(b.reverse), "0M"});
                ++result.edges_added;
            }
        }
    model.edges = std::move(kept);

    // ACCEPTANCE: nothing is written unless all of these hold.
    // 1. Every path spells exactly what it spelled after folding.
    std::vector<std::string> broken;
    for (std::size_t i = 0; i < model.paths.size(); ++i)
        if (spell(model, model.paths[i].steps) != spelled_before[i] && broken.size() < 8)
            broken.push_back(gfa_path_name(model.paths[i]));
    if (!broken.empty())
        throw std::runtime_error("panphorte --refine: the rewrite changed the sequence of " +
                                 cli::join_with_comma(broken) + " relative to the folded graph");
    // 2. Every path walks the same REP nodes, in the same orientation, the same number of times.
    for (std::size_t i = 0; i < model.paths.size(); ++i)
        if (rep_signature(model.paths[i].steps, self_loop) != reps_before[i])
            throw std::runtime_error("panphorte --refine: the REP copies of " + gfa_path_name(model.paths[i]) +
                                     " changed (node, orientation or count)");
    // 3. Every consecutive step pair has a link.
    for (const GfaPath& p : model.paths)
        for (std::size_t k = 1; k < p.steps.size(); ++k) {
            const PathStep& a = p.steps[k - 1];
            const PathStep& b = p.steps[k];
            if (have.count(arc_key(a.node_id, orient(a.reverse), b.node_id, orient(b.reverse))) == 0)
                throw std::runtime_error("panphorte --refine: path " + gfa_path_name(p) +
                                         " steps across a missing link");
        }
    // 4. Carriers with the same flanks around the same array block walk the same route.
    for (const SiteEdit& e : edits) {
        std::unordered_map<std::string, Steps> route_of;
        for (const auto& [idx, key] : e.group_key) {
            const Steps& route = std::get<2>(e.per_path.at(idx));
            auto [it, fresh] = route_of.emplace(key, route);
            Steps flipped = route;   // a reverse-crossing path carries the same route reverse-complemented
            std::reverse(flipped.begin(), flipped.end());
            for (PathStep& st : flipped) st.reverse = !st.reverse;
            if (!fresh && !same_steps(it->second, route) && !same_steps(it->second, flipped))
                throw std::runtime_error("panphorte --refine: identical haplotypes walk different routes at bubble " +
                                         std::to_string(result.rows[e.row_index].bubble_id));
        }
    }
    // 5. No rewritten path's folded route is still walkable, or the site would be represented twice.
    {
        std::unordered_set<std::string> live_edges;
        for (const GfaEdge& e : model.edges) live_edges.insert(arc_key(e.from, e.from_orient, e.to, e.to_orient));
        for (const auto& [idx, e, old] : old_spans) {
            // Every flank step moved to a new node and the array block did not, so the old route differs
            // from the new one exactly when the path had a flank step.
            if (edits[e].changed.count(idx) == 0) continue;
            bool walkable = true;
            for (std::size_t k = 0; k < old.size() && walkable; ++k) {
                if (!model.has_node(old[k].node_id)) walkable = false;
                else if (k > 0 && live_edges.count(arc_key(old[k - 1].node_id, orient(old[k - 1].reverse),
                                                           old[k].node_id, orient(old[k].reverse))) == 0)
                    walkable = false;
            }
            if (walkable)
                throw std::runtime_error("panphorte --refine: the folded route of a rewritten path is still walkable "
                                         "at bubble " + std::to_string(result.rows[edits[e].row_index].bubble_id) +
                                         ", so the site would be represented twice");
        }
    }
    return result;
}

} // namespace panvar
