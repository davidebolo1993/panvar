#include "panvar/refine.hpp"

#include "panvar/bubbles.hpp"
#include "panvar/cli_utils.hpp"
#include "panvar/gfa.hpp"
#include "panvar/gfa_io.hpp"
#include "panvar/graph_sort.hpp"
#include "panvar/graph_utils.hpp"
#include "panvar/gtf.hpp"
#include "panvar/integrated_snarls.hpp"
#include "panvar/output.hpp"
#include "panvar/poa.hpp"

#include <algorithm>
#include <cctype>
#include <cstdio>
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

std::string to_lower(std::string s) {
    for (char& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}

std::string path_display_name(const GfaPath& p) {
    if (p.type == 'P') return p.name;
    std::string n = p.sample;
    if (!p.hap.empty()) n += "#" + p.hap;
    if (!p.seqid.empty()) n += "#" + p.seqid;
    return n;
}

std::string spell(const GfaModel& model, const Steps& steps) {
    std::string out;
    for (const auto& s : steps) {
        auto it = model.seq.find(s.node_id);
        if (it == model.seq.end()) continue;
        out += s.reverse ? reverse_complement(it->second) : it->second;
    }
    return out;
}

// Nodes carrying an L self-edge (panphorte's folded copy-number REP nodes).
std::unordered_set<std::string> self_loop_nodes(const GfaModel& model) {
    std::unordered_set<std::string> out;
    for (const auto& e : model.edges) {
        if (e.from == e.to) out.insert(e.from);
    }
    return out;
}

// oriented steps between anchors a,b, plus their [lo,hi] positions in the walk
struct Interior {
    Steps steps;
    std::size_t lo = 0;
    std::size_t hi = 0;
};
std::optional<Interior> interior_between(const Steps& steps, const std::string& a, const std::string& b) {
    std::size_t ia = steps.size(), ib = steps.size();
    for (std::size_t i = 0; i < steps.size(); ++i) {
        if (ia == steps.size() && steps[i].node_id == a) ia = i;
        if (ib == steps.size() && steps[i].node_id == b) ib = i;
    }
    if (ia == steps.size() || ib == steps.size()) return std::nullopt;
    const std::size_t lo = std::min(ia, ib), hi = std::max(ia, ib);
    Interior out;
    out.lo = lo;
    out.hi = hi;
    out.steps.assign(steps.begin() + static_cast<long>(lo) + 1, steps.begin() + static_cast<long>(hi));
    return out;
}

// an interior split at its REP blocks: the residual runs (segments, one more than there are REP
// blocks), the REP runs themselves, and the ordered REP-id skeleton every haplotype must share.
struct Split {
    std::vector<std::string> skeleton;
    std::vector<Steps> segments;
    std::vector<Steps> rep_blocks;
};
Split split_by_rep(const Steps& interior, const std::unordered_set<std::string>& rep_set) {
    Split sp;
    Steps cur;
    std::size_t i = 0;
    while (i < interior.size()) {
        const std::string& node = interior[i].node_id;
        if (rep_set.count(node)) {
            sp.segments.push_back(cur);
            cur.clear();
            std::size_t j = i;
            Steps block;
            while (j < interior.size() && interior[j].node_id == node) {
                block.push_back(interior[j]);
                ++j;
            }
            sp.skeleton.push_back(node);
            sp.rep_blocks.push_back(block);
            i = j;
        } else {
            cur.push_back(interior[i]);
            ++i;
        }
    }
    sp.segments.push_back(cur);
    return sp;
}

// collapse an MSA into a compact sub-graph: the new nodes (id,seq) and each row's node path.
// ids start at next_id, new nodes are all '+' oriented.
struct InteriorGraph {
    std::vector<std::pair<std::string, std::string>> added;
    std::vector<std::vector<std::string>> row_paths;
};
InteriorGraph build_interior_graph(const std::vector<std::string>& rows, long long& next_id) {
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
    // unitig compaction: chain u->v into one node when that edge is the only one either sees and
    // exactly the same rows traverse both (a stretch of columns shared by one set of haplotypes).
    std::map<Cell, Cell> merged;
    std::function<Cell(Cell)> find = [&](Cell x) -> Cell {
        auto it = merged.find(x);
        while (it != merged.end() && !(it->second == x)) {
            x = it->second;
            it = merged.find(x);
        }
        return x;
    };
    for (const auto& kv : users) {  // ordered by (col,base)
        const Cell& u = kv.first;
        auto oi = out_edges.find(u);
        if (oi == out_edges.end() || oi->second.size() != 1) continue;
        const Cell v = *oi->second.begin();
        auto ii = in_edges.find(v);
        if (ii == in_edges.end() || ii->second.size() != 1) continue;
        if (users.at(u) == users.at(v) && users.at(u) == edge_users.at({u, v})) {
            merged[v] = find(u);
        }
    }
    std::map<Cell, std::vector<Cell>> group_cells;
    for (const auto& kv : users) group_cells[find(kv.first)].push_back(kv.first);

    InteriorGraph out;
    std::map<Cell, std::string> gid;
    for (auto& kv : group_cells) {  // ordered by leader (col,base)
        auto& cells = kv.second;
        std::sort(cells.begin(), cells.end());
        std::string id = std::to_string(next_id++);
        std::string seq;
        for (const auto& c : cells) seq += c.second;
        gid[kv.first] = id;
        out.added.emplace_back(id, seq);
    }
    out.row_paths.resize(S);
    for (std::size_t s = 0; s < S; ++s) {
        std::vector<std::string> p;
        for (const auto& cell : row_cells[s]) {
            const std::string& g = gid[find(cell)];
            if (p.empty() || p.back() != g) p.push_back(g);
        }
        out.row_paths[s] = std::move(p);
    }
    return out;
}

// One rebuilt region: nodes to drop, nodes to add, and per (path index) the interior replacement.
struct RegionEdit {
    std::unordered_set<std::string> interior_rm;
    std::vector<std::pair<std::string, std::string>> added;
    std::unordered_map<std::size_t, std::tuple<std::size_t, std::size_t, Steps>> per_path;  // idx -> (lo,hi,new)
};

// Fuse bubbles that share a boundary node into region components (an over-split event spans snarls).
std::vector<std::vector<Bubble>> components(const std::vector<Bubble>& rows) {
    std::unordered_map<std::size_t, std::size_t> parent;
    for (const auto& b : rows) parent[b.id] = b.id;
    std::function<std::size_t(std::size_t)> find = [&](std::size_t x) {
        while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; }
        return x;
    };
    std::unordered_map<std::string, std::vector<std::size_t>> ep;
    for (const auto& b : rows) { ep[b.source].push_back(b.id); ep[b.sink].push_back(b.id); }
    for (const auto& kv : ep) {
        for (std::size_t k = 1; k < kv.second.size(); ++k) parent[find(kv.second[k])] = find(kv.second[0]);
    }
    std::unordered_map<std::size_t, std::vector<Bubble>> comps;
    for (const auto& b : rows) comps[find(b.id)].push_back(b);
    std::vector<std::vector<Bubble>> out;
    out.reserve(comps.size());
    for (auto& kv : comps) out.push_back(std::move(kv.second));
    return out;
}

// The two outer anchors (endpoints appearing once across pooled source/sink) + the interior node set.
struct Anchors {
    bool ok = false;
    std::string a, b;
    std::unordered_set<std::string> interior;
};
Anchors region_anchors(const std::vector<Bubble>& comp) {
    std::vector<std::string> endpoints;
    std::unordered_set<std::string> interior;
    for (const auto& r : comp) {
        endpoints.push_back(r.source);
        endpoints.push_back(r.sink);
        interior.insert(r.source);
        interior.insert(r.sink);
        for (const auto& n : r.inside) interior.insert(n);
    }
    std::unordered_map<std::string, int> cnt;
    for (const auto& e : endpoints) ++cnt[e];
    std::vector<std::string> once;
    for (const auto& e : endpoints) {
        if (cnt[e] == 1 && std::find(once.begin(), once.end(), e) == once.end()) once.push_back(e);
    }
    Anchors out;
    if (once.size() != 2) return out;
    out.ok = true;
    out.a = std::min(once[0], once[1]);
    out.b = std::max(once[0], once[1]);
    interior.erase(out.a);
    interior.erase(out.b);
    out.interior = std::move(interior);
    return out;
}

double median_len(std::vector<std::size_t> v) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    const std::size_t m = v.size() / 2;
    return (v.size() % 2) ? static_cast<double>(v[m])
                          : 0.5 * static_cast<double>(v[m - 1] + v[m]);
}

// Try to rebuild one region. Returns a RegionEdit (with note in *note) or nullopt if skipped.
std::optional<RegionEdit> process_region(const std::vector<Bubble>& comp, const GfaModel& model,
                                         const std::vector<std::size_t>& ref_and_paths, std::size_t ref_idx,
                                         const std::unordered_set<std::string>& self_loop,
                                         const RefineOptions& opt, long long& next_id, std::string& note) {
    std::string ids;
    for (const auto& b : comp) ids += (ids.empty() ? "" : ",") + std::to_string(b.id);
    const std::string tag = "bubbles [" + ids + "]";
    if (comp.size() < opt.min_bubbles) { note = tag + ": skip (single bubble; --min-bubbles)"; return std::nullopt; }

    const Anchors anch = region_anchors(comp);
    if (!anch.ok) { note = tag + ": skip (not a clean linear chain)"; return std::nullopt; }
    std::unordered_set<std::string> rep_set;
    for (const auto& n : anch.interior) if (self_loop.count(n)) rep_set.insert(n);

    // gather each traversing path's interior; skip on any unfolded (non-REP) revisit >=2x
    std::unordered_map<std::size_t, Interior> per_path;
    for (std::size_t idx : ref_and_paths) {
        auto got = interior_between(model.paths[idx].steps, anch.a, anch.b);
        if (!got) continue;
        std::unordered_map<std::string, int> c;
        for (const auto& s : got->steps) ++c[s.node_id];
        for (const auto& kv : c) {
            if (kv.second >= 2 && !rep_set.count(kv.first)) {
                note = tag + ": skip (unfolded dup: node " + kv.first + " x" + std::to_string(kv.second) + ")";
                return std::nullopt;
            }
        }
        per_path[idx] = std::move(*got);
    }
    if (!per_path.count(ref_idx)) { note = tag + ": skip (reference doesn't traverse region)"; return std::nullopt; }

    // split each interior at REP blocks; all traversing paths must share the ref REP skeleton
    std::unordered_map<std::size_t, Split> split;
    for (auto& kv : per_path) split[kv.first] = split_by_rep(kv.second.steps, rep_set);
    const std::vector<std::string>& ref_skel = split.at(ref_idx).skeleton;
    for (const auto& kv : split) {
        if (kv.second.skeleton != ref_skel) {
            note = tag + ": skip (REP skeleton mismatch)";
            return std::nullopt;
        }
    }
    const std::size_t nseg = ref_skel.size() + 1;

    // size / diversity guards per residual segment position
    for (std::size_t pos = 0; pos < nseg; ++pos) {
        std::vector<std::size_t> lens;
        std::set<std::string> distinct;
        for (const auto& kv : split) {
            std::string s = spell(model, kv.second.segments[pos]);
            lens.push_back(s.size());
            distinct.insert(s);
        }
        if (median_len(lens) > static_cast<double>(opt.max_poa_bp)) {
            note = tag + ": skip (segment " + std::to_string(pos) + " median > max-poa-bp)";
            return std::nullopt;
        }
        if (distinct.size() > opt.max_walks) {
            note = tag + ": skip (segment " + std::to_string(pos) + " too diverse, " +
                   std::to_string(distinct.size()) + " walks)";
            return std::nullopt;
        }
    }

    // rebuild each residual segment; map distinct seq -> rebuilt node-id path
    RegionEdit edit;
    std::vector<std::unordered_map<std::string, std::vector<std::string>>> seg_to_ids(nseg);
    for (std::size_t pos = 0; pos < nseg; ++pos) {
        std::vector<std::string> distinct;
        std::unordered_set<std::string> seen;
        for (std::size_t idx : ref_and_paths) {  // deterministic, reference-first order
            auto it = split.find(idx);
            if (it == split.end()) continue;
            std::string s = spell(model, it->second.segments[pos]);
            if (seen.insert(s).second) distinct.push_back(s);
        }
        std::vector<std::string> rows = poa_msa(distinct);
        InteriorGraph ig = build_interior_graph(rows, next_id);
        for (auto& a : ig.added) edit.added.push_back(a);
        for (std::size_t k = 0; k < distinct.size(); ++k) seg_to_ids[pos][distinct[k]] = ig.row_paths[k];
    }

    // reassemble each traversing path's interior: residual node ids + verbatim REP runs
    for (const auto& kv : per_path) {
        const std::size_t idx = kv.first;
        const Split& sp = split.at(idx);
        Steps rebuilt;
        for (std::size_t pos = 0; pos < nseg; ++pos) {
            const std::string s = spell(model, sp.segments[pos]);
            for (const auto& id : seg_to_ids[pos].at(s)) rebuilt.push_back(PathStep{id, false});
            if (pos < sp.rep_blocks.size()) {
                for (const auto& st : sp.rep_blocks[pos]) rebuilt.push_back(st);  // REP run verbatim
            }
        }
        edit.per_path[idx] = std::make_tuple(kv.second.lo, kv.second.hi, std::move(rebuilt));
    }

    for (const auto& n : anch.interior) if (!rep_set.count(n)) edit.interior_rm.insert(n);
    note = tag + ": REBUILT [" + (rep_set.empty() ? "plain" : "rep-residual") + "] " +
           std::to_string(ref_skel.size()) + " REP block(s), " + std::to_string(nseg) + " residual seg(s) -> +" +
           std::to_string(edit.added.size()) + " nodes";
    return edit;
}

}  // namespace

RefineSummary refine_graph(const RefineOptions& options) {
    RefineSummary summary;
    GfaModel model = read_gfa_model(options.gfa_path);
    std::vector<Bubble> bubbles = read_bubbles_csv(options.bubbles_csv_in);
    const std::unordered_set<std::string> self_loop = self_loop_nodes(model);

    // resolve reference path (name or case-insensitive substring) -> index
    const std::string ref_q = to_lower(options.reference_path);
    std::size_t ref_idx = model.paths.size();
    for (std::size_t i = 0; i < model.paths.size(); ++i) {
        if (to_lower(path_display_name(model.paths[i])).find(ref_q) != std::string::npos) { ref_idx = i; break; }
    }
    if (ref_idx == model.paths.size()) {
        throw std::runtime_error("refine: no path matching --reference-path '" + options.reference_path + "'");
    }

    // path indices to consider per region: reference first, then all others
    std::vector<std::size_t> path_order;
    path_order.push_back(ref_idx);
    for (std::size_t i = 0; i < model.paths.size(); ++i) if (i != ref_idx) path_order.push_back(i);

    // region components
    std::vector<std::vector<Bubble>> comps;
    if (options.only_bubble_ids.empty()) {
        comps = components(bubbles);
    } else {
        std::unordered_set<std::string> want(options.only_bubble_ids.begin(), options.only_bubble_ids.end());
        std::vector<Bubble> picked;
        for (const auto& b : bubbles) if (want.count(std::to_string(b.id))) picked.push_back(b);
        if (!picked.empty()) comps.push_back(picked);
    }

    long long next_id = 0;
    for (const auto& kv : model.seq) {
        char* end = nullptr;
        const long long v = std::strtoll(kv.first.c_str(), &end, 10);
        if (end != kv.first.c_str() && *end == '\0' && v > next_id) next_id = v;
    }
    ++next_id;

    std::vector<RegionEdit> edits;
    cli::ProgressBar progress(options.quiet ? "" : "Refining regions", comps.size());
    for (const auto& comp : comps) {
        std::string note;
        auto edit = process_region(comp, model, path_order, ref_idx, self_loop, options, next_id, note);
        if (edit) { edits.push_back(std::move(*edit)); ++summary.regions_rebuilt; }
        else ++summary.regions_skipped;
        progress.tick();
    }

    // splice every rebuilt interior back into the model
    std::unordered_set<std::string> removed;
    for (const auto& e : edits) removed.insert(e.interior_rm.begin(), e.interior_rm.end());
    for (const auto& e : edits) {
        for (const auto& a : e.added) {
            model.seq[a.first] = a.second;
            model.node_order.push_back(a.first);
            ++summary.nodes_added;
        }
    }
    // per-path edits, applied right-to-left so indices stay valid
    std::unordered_map<std::size_t, std::vector<std::tuple<std::size_t, std::size_t, Steps>>> path_edits;
    for (const auto& e : edits) {
        for (const auto& kv : e.per_path) path_edits[kv.first].push_back(kv.second);
    }
    for (auto& kv : path_edits) {
        std::sort(kv.second.begin(), kv.second.end(),
                  [](const auto& x, const auto& y) { return std::get<0>(x) > std::get<0>(y); });
        Steps& st = model.paths[kv.first].steps;
        for (const auto& ed : kv.second) {
            const std::size_t lo = std::get<0>(ed), hi = std::get<1>(ed);
            const Steps& ins = std::get<2>(ed);
            Steps next;
            next.insert(next.end(), st.begin(), st.begin() + static_cast<long>(lo) + 1);
            next.insert(next.end(), ins.begin(), ins.end());
            next.insert(next.end(), st.begin() + static_cast<long>(hi), st.end());
            st = std::move(next);
        }
    }
    // drop removed nodes -- but only those no path still walks.
    //
    // Only paths spanning BOTH anchors are rewritten above. A path that enters the interior without
    // crossing both anchors keeps its original steps, so deleting an interior node it references leaves
    // a path pointing at a node that no longer exists: a silently corrupt graph that every downstream
    // module then reads. Count references after the rewrite and keep anything still in use.
    if (!removed.empty()) {
        std::unordered_set<std::string> still_used;
        for (const auto& pth : model.paths) {
            for (const auto& st : pth.steps) {
                if (removed.count(st.node_id) != 0) still_used.insert(st.node_id);
            }
        }
        if (!still_used.empty()) {
            for (const auto& n : still_used) removed.erase(n);
            summary.nodes_retained_referenced = still_used.size();
        }
    }
    if (!removed.empty()) {
        for (const auto& n : removed) { model.seq.erase(n); ++summary.nodes_removed; }
        model.node_order.erase(std::remove_if(model.node_order.begin(), model.node_order.end(),
                                              [&](const std::string& n) { return removed.count(n) > 0; }),
                               model.node_order.end());
    }
    // rebuild edges: keep original edges not touching a removed node, then add any path adjacency the
    // graph now lacks (as 0M links)
    auto edge_key = [](const std::string& a, char oa, const std::string& b, char ob) {
        return a + oa + "\t" + b + ob;
    };
    std::vector<GfaEdge> kept;
    std::unordered_set<std::string> keys;
    for (const auto& e : model.edges) {
        if (removed.count(e.from) || removed.count(e.to)) continue;
        kept.push_back(e);
        keys.insert(edge_key(e.from, e.from_orient, e.to, e.to_orient));
    }
    for (const auto& p : model.paths) {
        for (std::size_t k = 1; k < p.steps.size(); ++k) {
            const PathStep& s1 = p.steps[k - 1];
            const PathStep& s2 = p.steps[k];
            const char o1 = s1.reverse ? '-' : '+';
            const char o2 = s2.reverse ? '-' : '+';
            const std::string key = edge_key(s1.node_id, o1, s2.node_id, o2);
            if (keys.insert(key).second) {
                kept.push_back(GfaEdge{s1.node_id, o1, s2.node_id, o2, "0M"});
            }
        }
    }
    model.edges = std::move(kept);

    // re-sort, re-snarl and write the output family, the same way panphorte does
    GraphSortOptions sort_opts;
    sort_opts.reference_path = options.reference_path;
    sort_opts.flip = !options.no_flip;
    sort_graph_reference(model, sort_opts);

    const std::string sorted_gfa = options.out_prefix + ".normalized.sorted.gfa";
    write_gfa_model(sorted_gfa, model);

    BubbleCallOptions bopts;
    bopts.reference_path = options.reference_path;
    bopts.snarl_pairs_override = find_top_level_snarls_cactus(snarl_input_from_model(model));
    bopts.quiet = options.quiet;
    ParseGfaOptions parse_options;
    parse_options.include_paths = true;
    parse_options.include_sequences = true;
    const Graph sorted_graph = parse_gfa(sorted_gfa, parse_options);
    const BubbleCallReport report = call_bubbles_report(sorted_graph, bopts);
    write_bubbles_csv(options.out_prefix + ".bubbles.csv", report.bubbles);
    write_bandage_node_colors_csv(options.out_prefix + ".bandage_nodes.csv", report.bubbles,
                                  report.non_snp_bubbles);
    if (!options.gtf_path.empty()) {
        emit_gene_annotation(sorted_graph, options.reference_path, options.gtf_path,
                             options.out_prefix + ".bandage_genes.csv");
    }
    summary.bubbles_after = report.bubbles.size();
    return summary;
}

}  // namespace panvar
