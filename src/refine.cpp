#include "panvar/refine.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>

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
#include <cstdint>

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

// Folded repeats, held fixed: node -> the folded site it belongs to. A run of steps on one site's nodes is
// one repeat block. The sites come from panphorte's unit graph table (<prefix>.panphorte.unit_graph.tsv),
// every node of every unit graph grouped by site. A self-loop alone does not mark a fold (a unit graph of
// several nodes has none), and realigning a fold as plain sequence would lay its copies out in a line.
// Nodes with an L self-edge the table does not name are each a group of their own.
using RepGroups = std::unordered_map<std::string, std::string>;

std::vector<std::string> split_tsv(const std::string& line) {
    std::vector<std::string> out;
    std::size_t start = 0;
    for (;;) {
        const std::size_t tab = line.find('\t', start);
        out.push_back(line.substr(start, tab == std::string::npos ? std::string::npos : tab - start));
        if (tab == std::string::npos) return out;
        start = tab + 1;
    }
}

std::vector<std::string> split_list(const std::string& s) {   // ';'-separated ids, '.' for none
    std::vector<std::string> out;
    if (s.empty() || s == ".") return out;
    std::size_t start = 0;
    for (;;) {
        const std::size_t semi = s.find(';', start);
        out.push_back(s.substr(start, semi == std::string::npos ? std::string::npos : semi - start));
        if (semi == std::string::npos) return out;
        start = semi + 1;
    }
}

// A panphorte table, read whole so it can be written back with refine's renumbered ids.
struct RepTable {
    std::string path;
    std::vector<std::string> header;
    std::vector<std::vector<std::string>> rows;
    std::size_t col(const std::string& name) const {
        const auto it = std::find(header.begin(), header.end(), name);
        return it == header.end() ? SIZE_MAX : static_cast<std::size_t>(it - header.begin());
    }
};

std::optional<RepTable> read_rep_table(const std::string& path) {
    std::ifstream in(path);
    if (!in) return std::nullopt;
    RepTable t;
    t.path = path;
    std::string line;
    if (!std::getline(in, line)) throw std::runtime_error("refine: " + path + " is empty");
    if (!line.empty() && line.back() == '\r') line.pop_back();
    t.header = split_tsv(line);
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        t.rows.push_back(split_tsv(line));
        if (t.rows.back().size() != t.header.size())
            throw std::runtime_error("refine: " + path + " has a row with " + std::to_string(t.rows.back().size()) +
                                     " columns, the header " + std::to_string(t.header.size()));
    }
    return t;
}

// The id columns of panphorte's unit graph table: one id, or a ';' list.
const std::set<std::string> kUnitIdColumns{"anchor_nodes", "nodes"};

// One crossing of a region: the steps strictly between an anchor walked into the region and the next
// step, which must be the other anchor walked out of it. `lo`/`hi` are the anchor indices in the walk.
// A crossing entered at the second anchor is canonicalized (steps reversed, orientation toggled) and
// `reversed` records it, so the rebuilt interior can be spliced back in the path's own direction.
struct Crossing {
    Steps steps;
    std::size_t lo = 0;
    std::size_t hi = 0;
    bool reversed = false;
};
struct PathCrossings {
    std::vector<Crossing> crossings;
    std::size_t outside = 0;   // interior steps no crossing covers
};
// Every crossing a path makes, in walk order. A haplotype crossing a region twice (the two copies of a
// duplication) has both rewritten, so no copy is left on the old nodes and the copies keep sharing nodes.
// A crossing never contains an anchor, so a duplication of the anchors is not realigned into an insertion.
//
// `a_in_rev`/`b_in_rev` are each anchor's orientation walking into the region; leaving through an anchor
// walks it the other way. Consecutive anchor hits that are not an entry followed by an exit (re-entering
// by the anchor just left, or a stretch between two exits outside the region) are not crossings. Interior
// steps no crossing covers are counted in `outside`.
PathCrossings crossings_of(const Steps& steps, const std::string& a, bool a_in_rev, const std::string& b,
                           bool b_in_rev, const std::unordered_set<std::string>& interior) {
    PathCrossings out;
    std::vector<std::size_t> hits;
    for (std::size_t k = 0; k < steps.size(); ++k)
        if (steps[k].node_id == a || steps[k].node_id == b) hits.push_back(k);
    std::vector<char> covered(steps.size(), 0);
    for (std::size_t t = 1; t < hits.size(); ++t) {
        const PathStep& x = steps[hits[t - 1]];
        const PathStep& y = steps[hits[t]];
        if (x.node_id == y.node_id) continue;
        const bool x_in = x.reverse == (x.node_id == a ? a_in_rev : b_in_rev);
        const bool y_out = y.reverse != (y.node_id == a ? a_in_rev : b_in_rev);
        if (!x_in || !y_out) continue;
        Crossing c;
        c.lo = hits[t - 1];
        c.hi = hits[t];
        c.reversed = x.node_id == b;
        c.steps.assign(steps.begin() + static_cast<long>(c.lo) + 1, steps.begin() + static_cast<long>(c.hi));
        if (c.reversed) {
            std::reverse(c.steps.begin(), c.steps.end());
            for (PathStep& st : c.steps) st.reverse = !st.reverse;
        }
        for (std::size_t k = c.lo + 1; k < c.hi; ++k) covered[k] = 1;
        out.crossings.push_back(std::move(c));
    }
    for (std::size_t k = 0; k < steps.size(); ++k)
        if (!covered[k] && interior.count(steps[k].node_id) != 0) ++out.outside;
    return out;
}

// An interior split at its unit blocks: the residual runs (segments, one more than there are unit
// blocks), the unit runs themselves, and the ordered skeleton of folded sites every haplotype must share.
// A block is a run of steps on one folded site's nodes: its unit graph's copies (one node looped, for a
// unit graph of one node).
struct Split {
    std::vector<std::string> skeleton;
    std::vector<Steps> segments;
    std::vector<Steps> rep_blocks;
};
Split split_by_rep(const Steps& interior, const std::unordered_set<std::string>& rep_set, const RepGroups& groups) {
    Split sp;
    Steps cur;
    std::size_t i = 0;
    while (i < interior.size()) {
        const std::string& node = interior[i].node_id;
        if (rep_set.count(node)) {
            sp.segments.push_back(cur);
            cur.clear();
            const std::string& site = groups.at(node);
            std::size_t j = i;
            Steps block;
            while (j < interior.size() && rep_set.count(interior[j].node_id) &&
                   groups.at(interior[j].node_id) == site) {
                block.push_back(interior[j]);
                ++j;
            }
            sp.skeleton.push_back(site);
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

// Collapse an MSA into a compact sub-graph: the new nodes (id,seq) and each row's node path.
// Ids start at next_id; new nodes are all '+' oriented.
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
    // Unitig compaction: chain u->v into one node when that edge is the only one either sees and
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
    // idx -> one (lo, hi, new interior) per crossing
    std::unordered_map<std::size_t, std::vector<std::tuple<std::size_t, std::size_t, Steps>>> per_path;
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

// Node ids as numbers (shorter first), so "999" sorts before "1000".
bool numeric_less(const std::string& a, const std::string& b) {
    if (a.size() != b.size()) return a.size() < b.size();
    return a < b;
}

// The two outer anchors (endpoints appearing once across pooled source/sink) + the interior node set.
struct Anchors {
    bool ok = false;
    std::string a, b;
    bool a_in_rev = false, b_in_rev = false;   // orientation walking from the anchor into the region
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
    // Read in the reference's direction: the anchor that is a site's source first (bubble orients every site
    // along the reference), else the smaller id as a number. The region's direction decides which way its
    // crossings are read for alignment, and so on which side of a repeat an indel's gap lands.
    const auto is_source = [&](const std::string& n) {
        for (const auto& r : comp) if (r.source == n) return true;
        return false;
    };
    const bool s0 = is_source(once[0]), s1 = is_source(once[1]);
    const bool first0 = s0 != s1 ? s0 : numeric_less(once[0], once[1]);
    out.a = first0 ? once[0] : once[1];
    out.b = first0 ? once[1] : once[0];
    // A source is walked into its bubble in the reference's orientation; a sink is entered from outside
    // the other way round.
    const auto in_rev = [&](const std::string& n) {
        for (const auto& r : comp) {
            if (r.source == n) return r.source_reverse;
            if (r.sink == n) return !r.sink_reverse;
        }
        return false;
    };
    out.a_in_rev = in_rev(out.a);
    out.b_in_rev = in_rev(out.b);
    interior.erase(out.a);
    interior.erase(out.b);
    out.interior = std::move(interior);
    return out;
}

// What process_region decided for one region: the edit when it was rebuilt, and the report fields.
struct RegionResult {
    std::optional<RegionEdit> edit;
    std::string reason;                    // rebuilt: plain | unit_residual; skipped: the guard that fired
    // What a skip points at: the haplotype and the node behind it ("." when the reason names neither).
    std::string blocking_path = ".", blocking_node = ".";
    std::size_t short_revisits = 0;        // rebuilt: crossings realigned under --min-dup-bp
    // A rebuilt region's work: crossings rewritten, bases realigned over all of them, and the distinct
    // sequences abPOA aligned with their bases.
    std::size_t crossings = 0, realigned_bp = 0, distinct_sequences = 0, aligned_bp = 0;
    std::string source = "-", sink = "-";  // the outer anchors, in the reference's direction
};

RegionResult process_region(const std::vector<Bubble>& comp, const GfaModel& model,
                            const std::vector<std::size_t>& ref_and_paths, std::size_t ref_idx,
                            const RepGroups& rep_groups, const RefineOptions& opt, long long& next_id) {
    RegionResult res;
    if (!comp.empty()) { res.source = comp.front().source; res.sink = comp.back().sink; }
    const auto skip = [&](const char* reason, std::string path = ".", std::string node = ".") {
        res.reason = reason;
        res.blocking_path = std::move(path);
        res.blocking_node = std::move(node);
        return res;
    };
    const Anchors anch = region_anchors(comp);
    if (anch.ok) {
        res.source = anch.a;
        res.sink = anch.b;
    }
    if (comp.size() < opt.min_bubbles)
        return skip("min_bubbles");
    if (!anch.ok) return skip("not_a_clean_chain");

    // Every crossing of every path. A path that touches the interior outside its crossings cannot be
    // rewritten there, and its old nodes would keep their old links beside the rebuilt ones, a route the
    // losslessness check cannot see; so the region is skipped. The interior checked is the CSV's nodes plus
    // every node a crossing walks, since the splice removes both.
    std::unordered_map<std::size_t, std::vector<Crossing>> per_path;
    std::unordered_set<std::string> region_nodes = anch.interior;
    for (const std::size_t idx : ref_and_paths) {
        PathCrossings pc = crossings_of(model.paths[idx].steps, anch.a, anch.a_in_rev, anch.b, anch.b_in_rev,
                                        anch.interior);
        for (const Crossing& c : pc.crossings)
            for (const PathStep& st : c.steps) region_nodes.insert(st.node_id);
        if (!pc.crossings.empty()) per_path[idx] = std::move(pc.crossings);
    }
    std::string partial_path;   // the first path touching the interior outside a crossing
    for (const std::size_t idx : ref_and_paths)
        if (partial_path.empty() &&
            crossings_of(model.paths[idx].steps, anch.a, anch.a_in_rev, anch.b, anch.b_in_rev, region_nodes)
                    .outside > 0)
            partial_path = gfa_path_name(model.paths[idx]);
    if (const auto it = per_path.find(ref_idx); it != per_path.end() && it->second.front().reversed)
        std::swap(res.source, res.sink);
    if (!partial_path.empty()) return skip("partial_paths", partial_path);
    if (!per_path.count(ref_idx)) return skip("reference_not_traversing");

    std::unordered_set<std::string> rep_set;
    for (const auto& n : anch.interior) if (rep_groups.count(n)) rep_set.insert(n);
    for (const auto& kv : per_path)
        for (const Crossing& c : kv.second)
            for (const PathStep& st : c.steps) if (rep_groups.count(st.node_id)) rep_set.insert(st.node_id);

    // An unfolded duplication (a node that is not a folded repeat, walked twice within one crossing) is
    // copy-number signal that realigning would linearize, so it skips the region for the whole panel. A
    // crossing whose revisits carry fewer than --min-dup-bp bases is realigned instead; the default 0 skips
    // any revisit. Every visit after a node's first counts its full length.
    std::size_t short_revisits = 0;
    for (const std::size_t idx : ref_and_paths) {
        const auto it = per_path.find(idx);
        if (it == per_path.end()) continue;
        for (const Crossing& c : it->second) {
            std::unordered_map<std::string, int> cnt;
            for (const PathStep& st : c.steps) ++cnt[st.node_id];
            std::size_t dup_bp = 0;
            std::string first;
            for (const PathStep& st : c.steps) {
                int& n = cnt[st.node_id];
                if (n < 2 || rep_set.count(st.node_id)) continue;
                if (first.empty()) first = st.node_id;
                dup_bp += model.seq.at(st.node_id).size() * static_cast<std::size_t>(n - 1);
                n = 0;   // each revisited node counted once, for all its visits after the first
            }
            if (first.empty()) continue;
            if (dup_bp >= opt.min_dup_bp)
                return skip("unfolded_dup", gfa_path_name(model.paths[idx]), first);
            ++short_revisits;
        }
    }

    // Split each crossing at unit blocks; every crossing must share the reference's unit skeleton.
    struct CrossingSplit { std::size_t idx; const Crossing* crossing; Split split; };
    std::vector<CrossingSplit> splits;   // reference first, then path order, crossings in walk order
    for (const std::size_t idx : ref_and_paths) {
        const auto it = per_path.find(idx);
        if (it == per_path.end()) continue;
        for (const Crossing& c : it->second) splits.push_back({idx, &c, split_by_rep(c.steps, rep_set, rep_groups)});
    }
    const std::vector<std::string>& ref_skel = splits.front().split.skeleton;
    for (const CrossingSplit& cs : splits)
        if (cs.split.skeleton != ref_skel)
            return skip("unit_skeleton_mismatch", gfa_path_name(model.paths[cs.idx]));
    const std::size_t nseg = ref_skel.size() + 1;

    // Size and diversity guards per residual segment position.
    for (std::size_t pos = 0; pos < nseg; ++pos) {
        // abPOA is given the distinct sequences, so every guard is measured over those: the decision
        // must not depend on how many haplotypes carry each sequence.
        std::set<std::string> distinct;
        for (const CrossingSplit& cs : splits) distinct.insert(spell(model, cs.split.segments[pos]));

        std::size_t longest = 0, distinct_total = 0;
        for (const std::string& d : distinct) {
            longest = std::max(longest, d.size());
            distinct_total += d.size();
        }
        // POA cost is driven by the longest sequence, so that is what --max-poa-bp bounds.
        if (longest > opt.max_poa_bp)
            return skip("segment_too_long");
        // Estimated work, a budget of its own: a product test derived from --max-poa-bp and --max-walks
        // would be redundant with them.
        if (opt.max_poa_work > 0 &&
            static_cast<double>(longest) * static_cast<double>(distinct_total) >
                static_cast<double>(opt.max_poa_work))
            return skip("poa_work");
        if (distinct.size() > opt.max_walks)
            return skip("too_diverse");
    }

    // Rebuild each residual segment; map distinct seq -> rebuilt node-id path. Identical crossings,
    // including the two copies of a duplication, get identical node paths.
    RegionEdit edit;
    std::vector<std::unordered_map<std::string, std::vector<std::string>>> seg_to_ids(nseg);
    for (std::size_t pos = 0; pos < nseg; ++pos) {
        std::vector<std::string> distinct;
        std::unordered_set<std::string> seen;
        for (const CrossingSplit& cs : splits) {  // deterministic, reference-first order
            std::string s = spell(model, cs.split.segments[pos]);
            if (seen.insert(s).second) distinct.push_back(s);
        }
        res.distinct_sequences += distinct.size();
        for (const std::string& d : distinct) res.aligned_bp += d.size();
        std::vector<std::string> rows = poa_msa(distinct);
        // abPOA reads anything but A, C, G and T as N and writes upper case, so it only places the gaps:
        // each row gets its own bytes back over its non-gap columns, and soft-masked and IUPAC bases come
        // through unchanged.
        for (std::size_t k = 0; k < distinct.size(); ++k) {
            std::size_t j = 0;
            for (char& c : rows[k]) {
                if (c == '-') continue;
                if (j == distinct[k].size()) break;
                c = distinct[k][j++];
            }
            if (j != distinct[k].size() ||
                static_cast<std::size_t>(std::count_if(rows[k].begin(), rows[k].end(), [](char c) { return c != '-'; })) != j)
                throw std::logic_error("refine: the POA alignment does not hold sequence " + std::to_string(k) +
                                       " of segment " + std::to_string(pos));
        }
        InteriorGraph ig = build_interior_graph(rows, next_id);
        for (auto& a : ig.added) edit.added.push_back(a);
        for (std::size_t k = 0; k < distinct.size(); ++k) seg_to_ids[pos][distinct[k]] = ig.row_paths[k];
    }

    // Reassemble each crossing: residual node ids + verbatim unit runs.
    for (const CrossingSplit& cs : splits) {
        const Split& sp = cs.split;
        Steps rebuilt;
        for (std::size_t pos = 0; pos < nseg; ++pos) {
            const std::string s = spell(model, sp.segments[pos]);
            for (const auto& id : seg_to_ids[pos].at(s)) rebuilt.push_back(PathStep{id, false});
            if (pos < sp.rep_blocks.size()) {
                for (const auto& st : sp.rep_blocks[pos]) rebuilt.push_back(st);  // unit run verbatim
            }
        }
        // The crossing was canonicalized for alignment; one entered at the second anchor gets the rebuilt
        // interior back in its own direction, or the splice would write it reverse-complemented.
        if (cs.crossing->reversed) {
            std::reverse(rebuilt.begin(), rebuilt.end());
            for (PathStep& st : rebuilt) st.reverse = !st.reverse;
        }
        edit.per_path[cs.idx].emplace_back(cs.crossing->lo, cs.crossing->hi, std::move(rebuilt));
    }

    // Drop the old interior: the region's nodes and whatever its crossings walked, unit nodes excepted.
    // A path still walking one after the splice is a fatal error in the caller.
    for (const auto& n : anch.interior) if (!rep_set.count(n)) edit.interior_rm.insert(n);
    for (const CrossingSplit& cs : splits)
        for (const PathStep& st : cs.crossing->steps)
            if (!rep_set.count(st.node_id)) edit.interior_rm.insert(st.node_id);

    res.reason = rep_set.empty() ? "plain" : "unit_residual";
    res.crossings = splits.size();
    for (const CrossingSplit& cs : splits)
        for (const Steps& seg : cs.split.segments) res.realigned_bp += spell(model, seg).size();
    res.short_revisits = short_revisits;
    res.edit = std::move(edit);
    return res;
}

// ---- twin merge ---------------------------------------------------------------------------------
//
// A twin is a node with the same sequence as another node in the same place: both entered from the
// same oriented neighbours and leaving to the same ones. Graph builders leave them where forward- and
// reverse-strand contigs carrying the same bases were attached to different nodes, so node identity,
// which call, describe and benchmark compare, reports a difference that is only assembly strand.
//
// Merging twins never changes what a path spells, so it needs no realignment and also runs in the
// duplicated modules POA refinement skips. It must never add a route through the graph, which the
// acceptance checks cannot see, so the rule is narrow and the result is checked:
//   * identical sequence, byte for byte, in the same orientation;
//   * identical oriented in- and out-neighbours, apart from each node's own self-links: one full key
//     per group, so a group is never chained together from partial matches;
//   * self-links are compared by kind (tandem n+->n+, right fold-back n+->n-, left fold-back n-->n+):
//     a member joins only if the keeper's kinds include its own, so the merged node carries no loop
//     combination that did not exist;
//   * a self-linked member joins only if some path walks it twice in a row: call reads every
//     self-linked node as a repeat unit, and an unwalked loop merged into a walked node would hand its
//     haplotypes a loop they never took;
//   * never two nodes linked to each other, which is how a homopolymer run is drawn (T -> T -> T). The
//     key already ensures it, since a link puts each node's handle among the other's neighbours; the
//     check below only guards that;
//   * never a node panphorte's repeat table names: the table points at those ids, and a unit graph's
//     copies are counted on them.
// After each pass, every route through a merged node (predecessor, node, successor) must have existed
// through one of its members; anything else is a bug and stops the run.
struct TwinMerge {
    std::string kept;
    std::string merged;
    std::size_t length_bp = 0;
    std::size_t steps_rewritten = 0;
};

struct TwinMergeResult {
    std::vector<TwinMerge> merges;
    std::size_t paths_walking_several = 0;   // paths that walked two members of one group
};

TwinMergeResult merge_twin_nodes(GfaModel& model, const std::unordered_set<std::string>& fixed) {
    TwinMergeResult result;
    auto flip = [](char o) { return o == '+' ? '-' : '+'; };
    std::set<std::size_t> walking_several;   // paths that walked two members of one group, in any pass
    for (;;) {
        // Oriented adjacency over handles ("id+" / "id-"), from both readings of every L line.
        std::map<std::string, std::set<std::string>> pred, succ;
        std::set<std::pair<std::string, std::string>> linked;   // unordered pairs of distinct nodes
        std::map<std::string, std::set<char>> loop_kinds;        // 't' tandem, 'r' n+->n-, 'l' n-->n+
        for (const GfaEdge& e : model.edges) {
            const std::string a = e.from + e.from_orient, b = e.to + e.to_orient;
            const std::string ra = e.from + flip(e.from_orient), rb = e.to + flip(e.to_orient);
            succ[a].insert(b); pred[b].insert(a);
            succ[rb].insert(ra); pred[ra].insert(rb);
            if (e.from == e.to) {
                loop_kinds[e.from].insert(e.from_orient == e.to_orient ? 't'
                                          : e.from_orient == '+' ? 'r' : 'l');
            } else {
                linked.insert(std::minmax(e.from, e.to));
            }
        }
        std::set<std::string> walked_loop;
        for (const GfaPath& p : model.paths)
            for (std::size_t k = 1; k < p.steps.size(); ++k)
                if (p.steps[k].node_id == p.steps[k - 1].node_id) walked_loop.insert(p.steps[k].node_id);
        auto without_self = [](const std::set<std::string>& hs, const std::string& id) {
            std::vector<std::string> v;
            for (const std::string& h : hs)
                if (h.compare(0, h.size() - 1, id) != 0) v.push_back(h);
            return v;
        };
        // One full key per node: sequence, then in- and out-neighbours of its forward handle.
        std::map<std::string, std::vector<std::string>> groups;
        for (const std::string& id : model.node_order) {
            const auto sit = model.seq.find(id);
            if (sit == model.seq.end() || sit->second.empty() || fixed.count(id)) continue;
            const auto ins = without_self(pred[id + '+'], id);
            const auto outs = without_self(succ[id + '+'], id);
            if (ins.empty() && outs.empty()) continue;   // no place to share
            std::string key = sit->second;
            key += '\x1f';
            for (const std::string& n : ins) { key += n; key += ','; }
            key += '\x1f';
            for (const std::string& n : outs) { key += n; key += ','; }
            groups[key].push_back(id);
        }
        std::unordered_map<std::string, std::string> rename;   // merged -> kept
        std::vector<std::vector<std::string>> classes;
        for (auto& kv : groups) {
            std::vector<std::string> ids = kv.second;
            if (ids.size() < 2) continue;
            std::sort(ids.begin(), ids.end(), numeric_less);
            // Looped members join only when their loop is walked.
            std::vector<std::string> cand;
            for (const std::string& id : ids)
                if (!loop_kinds.count(id) || walked_loop.count(id)) cand.push_back(id);
            if (cand.size() < 2) continue;
            // The keeper carries the widest set of loop kinds; a member whose kinds it does not include
            // would give the merged node a combination no node had, so that member stays separate.
            std::string keep = cand.front();
            for (const std::string& id : cand)
                if (loop_kinds[id].size() > loop_kinds[keep].size()) keep = id;
            const std::set<char>& kk = loop_kinds[keep];
            std::vector<std::string> members{keep};
            for (const std::string& id : cand) {
                if (id == keep) continue;
                const std::set<char>& k = loop_kinds[id];
                if (!std::includes(kk.begin(), kk.end(), k.begin(), k.end())) continue;
                bool independent = true;
                for (const std::string& m : members)
                    if (linked.count(std::minmax(id, m))) { independent = false; break; }
                if (independent) members.push_back(id);
            }
            if (members.size() < 2) continue;
            for (const std::string& m : members)
                if (m != keep) rename[m] = keep;
            classes.push_back(members);
        }
        if (rename.empty()) break;

        // A path that walked two members of one group visits the merged node twice from now on.
        std::unordered_set<std::string> keepers;
        for (const auto& r : rename) keepers.insert(r.second);
        for (std::size_t pi = 0; pi < model.paths.size(); ++pi) {
            std::map<std::string, std::set<std::string>> seen;   // keeper -> members walked
            for (const PathStep& st : model.paths[pi].steps) {
                const auto it = rename.find(st.node_id);
                if (it != rename.end()) seen[it->second].insert(st.node_id);
                else if (keepers.count(st.node_id)) seen[st.node_id].insert(st.node_id);
            }
            for (const auto& kv : seen)
                if (kv.second.size() >= 2) { walking_several.insert(pi); break; }
        }

        std::unordered_map<std::string, std::size_t> moved;
        for (GfaPath& p : model.paths)
            for (PathStep& st : p.steps) {
                const auto it = rename.find(st.node_id);
                if (it == rename.end()) continue;
                ++moved[st.node_id];
                st.node_id = it->second;
            }
        // Redirect the links of merged nodes and drop the duplicates that creates, an L line and its
        // reverse reading being one link. Links between untouched nodes are left exactly as they were.
        std::unordered_set<std::string> canon;
        auto key_of = [&](const GfaEdge& e) {
            const std::string fwd = e.from + e.from_orient + "\t" + e.to + e.to_orient;
            const std::string rev = e.to + flip(e.to_orient) + "\t" + e.from + flip(e.from_orient);
            return std::min(fwd, rev);
        };
        std::vector<GfaEdge> edges;
        for (const GfaEdge& e : model.edges)
            if (!rename.count(e.from) && !rename.count(e.to)) canon.insert(key_of(e));
        for (GfaEdge e : model.edges) {
            const bool touched = rename.count(e.from) || rename.count(e.to);
            if (touched) {
                if (const auto it = rename.find(e.from); it != rename.end()) e.from = it->second;
                if (const auto it = rename.find(e.to); it != rename.end()) e.to = it->second;
                if (!canon.insert(key_of(e)).second) continue;
            }
            edges.push_back(std::move(e));
        }
        model.edges = std::move(edges);

        // Route check: every predecessor -> merged node -> successor route must have existed through one
        // member. pred/succ still describe the graph before this pass; neighbours are compared after
        // renaming, so a member's own loop counts as a loop on the merged node.
        auto ren = [&](const std::string& h) {
            const auto it = rename.find(h.substr(0, h.size() - 1));
            return it == rename.end() ? h : it->second + h.back();
        };
        for (const std::vector<std::string>& cls : classes) {
            const std::string& keep = cls.front();
            for (const char o : {'+', '-'}) {
                std::set<std::string> new_in, new_out;
                for (const std::string& m : cls) {
                    for (const std::string& h : pred[m + o]) new_in.insert(ren(h));
                    for (const std::string& h : succ[m + o]) new_out.insert(ren(h));
                }
                for (const std::string& u : new_in)
                    for (const std::string& w : new_out) {
                        bool existed = false;
                        for (const std::string& m : cls) {
                            bool in_ok = false, out_ok = false;
                            for (const std::string& h : pred[m + o]) if (ren(h) == u) { in_ok = true; break; }
                            for (const std::string& h : succ[m + o]) if (ren(h) == w) { out_ok = true; break; }
                            if (in_ok && out_ok) { existed = true; break; }
                        }
                        if (!existed)
                            throw std::runtime_error("refine: merging twin nodes into " + keep +
                                                     " would create the route " + u + " -> " + keep + o +
                                                     " -> " + w + ", which no member had");
                    }
            }
        }

        std::vector<std::string> gone;
        for (const auto& kv : rename) gone.push_back(kv.first);
        std::sort(gone.begin(), gone.end(), numeric_less);
        for (const std::string& g : gone) {
            result.merges.push_back(TwinMerge{rename.at(g), g, model.seq.at(g).size(), moved[g]});
            model.seq.erase(g);
        }
        model.node_order.erase(std::remove_if(model.node_order.begin(), model.node_order.end(),
                                              [&](const std::string& n) { return rename.count(n) > 0; }),
                               model.node_order.end());
    }
    result.paths_walking_several = walking_several.size();
    return result;
}

}  // namespace

RefineSummary refine_graph(const RefineOptions& options) {
    RefineSummary summary;
    const auto say = [&](const std::string& message) {
        if (!options.quiet) std::cerr << cli::log_prefix("refine") << message << '\n';
    };
    // No output may name an input: committing the outputs would overwrite it.
    {
        const std::string outs[] = {
            options.out_prefix + ".normalized.sorted.gfa",
            options.out_prefix + ".bubbles.csv",
            options.out_prefix + ".bandage_nodes.csv",
            options.out_prefix + ".bandage_genes.csv",
            options.out_prefix + ".refine.report.tsv",
            options.out_prefix + ".panphorte.unit_graph.tsv",
        };
        for (const std::string& in : {options.gfa_path, options.bubbles_csv_in, options.gtf_path,
                                      options.panphorte_prefix.empty() ? std::string()
                                          : options.panphorte_prefix + ".panphorte.unit_graph.tsv"}) {
            if (in.empty()) continue;
            std::error_code ec;
            const auto ip = std::filesystem::weakly_canonical(in, ec);
            if (ec || ip.empty()) continue;
            for (const std::string& o : outs) {
                std::error_code e2;
                const auto op = std::filesystem::weakly_canonical(o, e2);
                if (!e2 && ip == op)
                    throw std::runtime_error("refine: output '" + o + "' is the same file as input '" +
                                             in + "'");
            }
        }
    }

    cli::require_regular_input("refine", "-i", options.gfa_path);   // read as a model and again as a graph
    GfaModel model = read_gfa_model(options.gfa_path);
    // Read once, inversion sites included, so --bubble-id can name one: -c may be a pipe.
    std::vector<Bubble> bubbles = read_bubbles_csv(options.bubbles_csv_in, true);
    std::unordered_set<std::string> inversion_ids;
    for (const Bubble& b : bubbles)
        if (b.inversion_parent) inversion_ids.insert(std::to_string(b.id));
    bubbles.erase(std::remove_if(bubbles.begin(), bubbles.end(), [](const Bubble& b) { return b.inversion_parent; }),
                  bubbles.end());
    std::string resolved_reference;

    // The graph contract every module applies: refine rewrites the graph, so a malformed input must be
    // refused rather than emitted looking repaired. The Graph view also resolves the reference, so that
    // rule is shared too.
    {
        ParseGfaOptions popts;
        popts.include_paths = true;
        popts.include_sequences = true;
        const Graph gview = parse_gfa(options.gfa_path, popts);
        validate_graph_paths(gview, "refine", true, true);
        resolved_reference = resolve_reference_path_name(gview, options.reference_path, "refine");
        say("input: haplotypes: " + std::to_string(model.paths.size()) + "; nodes: " +
            std::to_string(model.node_order.size()) + "; edges: " + std::to_string(model.edges.size()));
        say("reference: " + resolved_reference);

        // The CSV and the graph are separate inputs, so the CSV's ids must be unique and every node it
        // names must be in the graph.
        std::vector<std::string> missing;
        std::unordered_set<std::size_t> seen_ids;
        for (const Bubble& b : bubbles) {
            if (!seen_ids.insert(b.id).second)
                throw std::runtime_error("refine: duplicate bubble id in the CSV: " +
                                         std::to_string(b.id));
            const auto check = [&](const std::string& n) {
                if (gview.nodes.find(n) == gview.nodes.end() && missing.size() < 8)
                    missing.push_back("bubble " + std::to_string(b.id) + " node " + n);
            };
            check(b.source);
            check(b.sink);
            for (const std::string& n : b.inside) check(n);
        }
        if (!missing.empty())
            throw std::runtime_error(
                "refine: the bubbles CSV describes nodes the GFA does not contain (" +
                cli::join_with_comma(missing) + "); the CSV and the graph are not the same graph");
    }

    // The folded repeats to hold fixed (see RepGroups). Each table must describe this graph.
    RepGroups rep_groups;
    std::unordered_set<std::string> table_nodes;   // nodes panphorte's tables name: kept out of twin merge
    std::vector<RepTable> rep_tables;
    if (!options.panphorte_prefix.empty()) {
        // Node ids alone do not tie a table to a graph: sorted graphs reuse small ids, so a table from another
        // run can name nodes this graph has. Each row must agree with the graph: its node count is its list's,
        // its anchor is among its nodes, its unit_bp is its unit's length, a one-node unit graph spells its
        // unit (on either strand, as the sort may have flipped it), and no node belongs to two sites.
        const auto refuse = [&](const std::string& path, const std::string& why) {
            throw std::runtime_error("refine: " + path + " does not describe this graph: " + why +
                                     "; give the table panphorte wrote with this graph (--panphorte-prefix)");
        };
        std::unordered_map<std::string, std::string> site_of;   // node -> the site the table puts it in
        if (std::optional<RepTable> t = read_rep_table(options.panphorte_prefix + ".panphorte.unit_graph.tsv")) {
            const std::size_t site_col = t->col("input_bubble_id"), bp_col = t->col("unit_bp"), unit_col = t->col("unit");
            const std::size_t nodes_col = t->col("nodes"), n_col = t->col("n_nodes"), anchor_col = t->col("anchor_nodes");
            for (const std::size_t c : {site_col, bp_col, unit_col, nodes_col, n_col, anchor_col})
                if (c == SIZE_MAX)
                    throw std::runtime_error("refine: " + t->path + " is not a panphorte unit graph table (columns "
                                             "input_bubble_id, unit_bp, unit, anchor_nodes, n_nodes, nodes)");
            for (const auto& row : t->rows) {
                for (std::size_t c = 0; c < t->header.size(); ++c) {
                    if (!kUnitIdColumns.count(t->header[c])) continue;
                    for (const std::string& n : split_list(row[c])) {
                        if (model.seq.find(n) == model.seq.end())
                            throw std::runtime_error("refine: " + t->path + " names node " + n +
                                                     ", which the GFA does not contain; the table and the "
                                                     "graph are not the same graph");
                        const auto [it, fresh] = site_of.emplace(n, row[site_col]);
                        if (!fresh && it->second != row[site_col])
                            refuse(t->path, "node " + n + " is in site " + it->second + " and in site " + row[site_col]);
                        table_nodes.insert(n);
                        rep_groups.emplace(n, "site " + row[site_col]);
                    }
                }
                const std::vector<std::string> nodes = split_list(row[nodes_col]);
                if (row[n_col] != std::to_string(nodes.size()))
                    refuse(t->path, "site " + row[site_col] + " lists " + std::to_string(nodes.size()) +
                                        " nodes, its row says " + row[n_col]);
                for (const std::string& a : split_list(row[anchor_col]))
                    if (std::find(nodes.begin(), nodes.end(), a) == nodes.end())
                        refuse(t->path, "anchor node " + a + " of site " + row[site_col] + " is not among its nodes");
                if (row[bp_col] != std::to_string(row[unit_col].size()))
                    refuse(t->path, "site " + row[site_col] + "'s unit is " + std::to_string(row[unit_col].size()) +
                                        " bp, its row says " + row[bp_col]);
                if (nodes.size() == 1) {
                    std::string got = row[unit_col], want = model.seq.at(nodes.front());
                    for (char& ch : got) ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
                    for (char& ch : want) ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
                    if (got != want && got != reverse_complement(want))
                        refuse(t->path, "unit node " + nodes.front() + " (" + std::to_string(want.size()) +
                                            " bp) does not spell its site's unit (" + std::to_string(got.size()) +
                                            " bp), on either strand");
                }
            }
            rep_tables.push_back(std::move(*t));
        }
        summary.rep_tables_read = rep_tables.size();
    }
    if (!rep_tables.empty()) {
        std::set<std::string> sites;
        for (const auto& kv : rep_groups) sites.insert(kv.second);
        say("folded repeats: " + std::to_string(sites.size()) + " sites, " + std::to_string(table_nodes.size()) +
            " nodes (" + rep_tables.front().path + ")");
    } else {
        // Without the table a unit graph of several nodes is realigned as plain sequence, so this is shown
        // even with -q.
        if (options.panphorte_prefix.empty())
            std::cerr << cli::log_prefix("refine") << "warning: no panphorte unit graph table: the bubbles CSV "
                      << "name does not end in .bubbles.csv, so give its prefix with --panphorte-prefix; folded "
                      << "repeats are recognised by their self-loop only\n";
        else
            std::cerr << cli::log_prefix("refine") << "warning: no panphorte unit graph table at '"
                      << options.panphorte_prefix << ".panphorte.unit_graph.tsv'; folded repeats are recognised by "
                      << "their self-loop only\n";
    }
    for (const auto& e : model.edges)
        if (e.from == e.to) rep_groups.emplace(e.from, "node " + e.from);

    // The resolved name, matched exactly: the resolver already prefers an exact name and refuses an
    // ambiguous substring.
    std::size_t ref_idx = model.paths.size();
    for (std::size_t i = 0; i < model.paths.size(); ++i) {
        if (gfa_path_name(model.paths[i]) == resolved_reference) { ref_idx = i; break; }
    }
    if (ref_idx == model.paths.size()) {
        throw std::runtime_error("refine: resolved reference '" + resolved_reference +
                                 "' is not a path of the model");
    }

    // Every input path's spelled sequence, kept so the rewrite can be checked lossless. Indexed by path
    // position rather than name, so the check does not depend on path names round-tripping.
    std::vector<std::string> spelled_before;
    spelled_before.reserve(model.paths.size());
    for (const GfaPath& p : model.paths) spelled_before.push_back(spell(model, p.steps));

    // Path indices to consider per region: reference first, then all others.
    std::vector<std::size_t> path_order;
    path_order.push_back(ref_idx);
    for (std::size_t i = 0; i < model.paths.size(); ++i) if (i != ref_idx) path_order.push_back(i);

    // Region components.
    std::vector<std::vector<Bubble>> comps;
    if (options.only_bubble_ids.empty()) {
        comps = components(bubbles);
    } else {
        std::unordered_set<std::string> want(options.only_bubble_ids.begin(), options.only_bubble_ids.end());
        std::vector<Bubble> picked;
        std::unordered_set<std::string> found;
        for (const auto& b : bubbles) {
            const std::string id = std::to_string(b.id);
            if (want.count(id)) { picked.push_back(b); found.insert(id); }
        }
        // An id the CSV lacks is an error: an unchanged output family would read as "nothing needed refining".
        std::vector<std::string> absent;
        for (const std::string& w : want) if (!found.count(w)) absent.push_back(w);
        // An inversion site is in the CSV but not read: refine does not realign one, the sites inside it are
        // realigned on their own.
        if (!absent.empty()) {
            std::vector<std::string> missing, inversions;
            for (const std::string& w : absent) (inversion_ids.count(w) ? inversions : missing).push_back(w);
            std::sort(missing.begin(), missing.end());
            std::sort(inversions.begin(), inversions.end());
            std::string msg;
            if (!missing.empty()) msg = "--bubble-id not present in the bubbles CSV: " + cli::join_with_comma(missing);
            if (!inversions.empty())
                msg += (msg.empty() ? "" : "; ") + std::string("--bubble-id names an inversion site, which refine "
                       "does not realign (give the sites inside it): ") + cli::join_with_comma(inversions);
            throw std::runtime_error("refine: " + msg);
        }
        // Selected ids go through the same component grouping as auto mode, so disjoint bubbles stay
        // separate regions.
        comps = components(picked);
    }

    {
        std::size_t n_sites = 0;
        for (const auto& comp : comps) n_sites += comp.size();
        say("sites: " + std::to_string(n_sites) + "; regions: " + std::to_string(comps.size()) +
            " (touching sites joined)");
    }

    long long next_id = 0;
    for (const auto& kv : model.seq) {
        char* end = nullptr;
        const long long v = std::strtoll(kv.first.c_str(), &end, 10);
        if (end != kv.first.c_str() && *end == '\0' && v > next_id) next_id = v;
    }
    ++next_id;

    // Regions sorted by first anchor, as a number, so node numbering and output bytes do not depend on
    // hash order and the report reads in id order.
    std::sort(comps.begin(), comps.end(), [](const std::vector<Bubble>& a, const std::vector<Bubble>& b) {
        if (a.empty() || b.empty()) return b.empty() < a.empty();
        if (a.front().source != b.front().source) return numeric_less(a.front().source, b.front().source);
        return numeric_less(a.front().sink, b.front().sink);
    });

    // Sites of the CSV that overlap (one's node inside another region) would be rewritten twice over the
    // same steps, so both are skipped. Fused chains share only their anchors, and those are one region.
    // A region's nodes are the CSV's plus every node a crossing walks between its anchors, since the edit
    // replaces those steps whether or not the CSV lists them.
    std::vector<std::string> overlap_of(comps.size());
    {
        std::unordered_map<std::string, std::size_t> owner;
        for (std::size_t c = 0; c < comps.size(); ++c) {
            std::set<std::string> nodes;
            for (const Bubble& b : comps[c]) {
                nodes.insert(b.source);
                nodes.insert(b.sink);
                nodes.insert(b.inside.begin(), b.inside.end());
            }
            if (const Anchors anch = region_anchors(comps[c]); anch.ok)
                for (const std::size_t idx : path_order)
                    for (const Crossing& x : crossings_of(model.paths[idx].steps, anch.a, anch.a_in_rev, anch.b,
                                                          anch.b_in_rev, anch.interior).crossings)
                        for (const PathStep& st : x.steps) nodes.insert(st.node_id);
            for (const std::string& n : nodes) {
                const auto [it, fresh] = owner.emplace(n, c);
                if (fresh || it->second == c) continue;
                const std::size_t o = it->second;
                // The node both regions hold.
                if (overlap_of[c].empty()) overlap_of[c] = n;
                if (overlap_of[o].empty()) overlap_of[o] = n;
            }
        }
    }

    // The report gives every region's decision and reason. All outputs are staged and moved into place
    // only once the whole run has succeeded, so a failure at any step leaves no partial family.
    cli::StagedOutputs staged("refine");
    const std::string report_path = staged.stage(options.out_prefix + ".refine.report.tsv");
    std::ofstream rep(report_path);
    if (!rep) throw std::runtime_error("refine: cannot write " + options.out_prefix + ".refine.report.tsv");
    rep << "region\tn_bubbles\tsource\tsink\tdecision\treason\tcrossings\trealigned_bp\tdistinct_sequences\t"
           "aligned_bp\tnodes_removed\tnodes_added\tshort_revisits\tblocking_path\tblocking_node\n";

    struct RegionRow { std::size_t n_bubbles; RegionResult result; std::size_t edit_idx; std::size_t nodes_removed = 0, nodes_added = 0; };
    std::vector<RegionRow> rows;
    std::vector<RegionEdit> edits;
    cli::ProgressBar progress(options.quiet ? "" : cli::log_prefix("refine") + "refining regions", comps.size());
    for (std::size_t c = 0; c < comps.size(); ++c) {
        const auto& comp = comps[c];
        RegionResult result;
        if (!overlap_of[c].empty()) {
            if (const Anchors anch = region_anchors(comp); anch.ok) {
                result.source = anch.a;
                result.sink = anch.b;
            } else {
                result.source = comp.front().source;
                result.sink = comp.back().sink;
            }
            result.reason = "overlapping_sites";
            result.blocking_node = overlap_of[c];
        } else {
            result = process_region(comp, model, path_order, ref_idx, rep_groups, options, next_id);
        }
        std::size_t edit_idx = edits.size();
        std::size_t nodes_removed = 0, nodes_added = 0;
        if (result.edit) {
            nodes_removed = result.edit->interior_rm.size();
            nodes_added = result.edit->added.size();
            edits.push_back(std::move(*result.edit));
            result.edit.reset();
            ++summary.regions_rebuilt;
        } else {
            edit_idx = static_cast<std::size_t>(-1);
            ++summary.regions_skipped;
        }
        rows.push_back({comp.size(), std::move(result), edit_idx, nodes_removed, nodes_added});
        progress.tick();
    }

    {
        // Regions by outcome and reason, with the bases realigned over every crossing and the bases of
        // the distinct sequences abPOA aligned.
        std::map<std::string, std::size_t> rebuilt_by, skipped_by;
        std::size_t realigned_bp = 0, aligned_bp = 0;
        for (const RegionRow& r : rows) {
            const bool rebuilt = r.edit_idx != static_cast<std::size_t>(-1);
            ++(rebuilt ? rebuilt_by : skipped_by)[r.result.reason];
            realigned_bp += r.result.realigned_bp;
            aligned_bp += r.result.aligned_bp;
        }
        const auto line = [&](const char* what, const std::map<std::string, std::size_t>& by) {
            std::size_t n = 0;
            std::string reasons;
            for (const auto& [reason, k] : by) {
                n += k;
                reasons += (reasons.empty() ? "" : "; ") + reason + ": " + std::to_string(k);
            }
            return std::string(what) + ": " + std::to_string(n) + (reasons.empty() ? "" : " (" + reasons + ")");
        };
        say(line("rebuilt", rebuilt_by));
        say(line("skipped", skipped_by));
        say("realigned bp, all crossings: " + std::to_string(realigned_bp) + "; aligned bp, distinct sequences: " +
            std::to_string(aligned_bp));
    }

    // Splice every rebuilt interior back into the model.
    std::unordered_set<std::string> removed;
    for (const auto& e : edits) removed.insert(e.interior_rm.begin(), e.interior_rm.end());
    for (const auto& e : edits) {
        for (const auto& a : e.added) {
            model.seq[a.first] = a.second;
            model.node_order.push_back(a.first);
            ++summary.nodes_added;
        }
    }
    // Per-path edits, applied right-to-left so indices stay valid.
    std::unordered_map<std::size_t, std::vector<std::tuple<std::size_t, std::size_t, Steps>>> path_edits;
    for (const auto& e : edits) {
        for (const auto& kv : e.per_path)
            path_edits[kv.first].insert(path_edits[kv.first].end(), kv.second.begin(), kv.second.end());
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
    // An old interior node a path still walks would keep its old links, so the old route would stay
    // walkable beside the rebuilt one. Regions with such paths are skipped, so this is a bug if it fires.
    for (const auto& pth : model.paths)
        for (const auto& st : pth.steps)
            if (removed.count(st.node_id) != 0)
                throw std::runtime_error("refine: path " + gfa_path_name(pth) + " still walks old interior node " +
                                         st.node_id + " after the rebuild; nothing was written");
    {
        std::size_t region_no = 0;
        for (const RegionRow& r : rows) {
            const bool rebuilt = r.edit_idx != static_cast<std::size_t>(-1);
            const auto num = [&](std::size_t v) { return rebuilt ? std::to_string(v) : std::string("."); };
            rep << ++region_no << '\t' << r.n_bubbles << '\t' << r.result.source << '\t' << r.result.sink
                << '\t' << (rebuilt ? "rebuilt" : "skipped") << '\t' << r.result.reason << '\t'
                << num(r.result.crossings) << '\t' << num(r.result.realigned_bp) << '\t'
                << num(r.result.distinct_sequences) << '\t' << num(r.result.aligned_bp) << '\t'
                << num(r.nodes_removed) << '\t' << num(r.nodes_added) << '\t' << num(r.result.short_revisits)
                << '\t' << r.result.blocking_path << '\t' << r.result.blocking_node << '\n';
        }
    }
    if (!removed.empty()) {
        for (const auto& n : removed) { model.seq.erase(n); ++summary.nodes_removed; }
        model.node_order.erase(std::remove_if(model.node_order.begin(), model.node_order.end(),
                                              [&](const std::string& n) { return removed.count(n) > 0; }),
                               model.node_order.end());
    }
    // Rebuild edges: keep original edges not touching a removed node, then add any path adjacency the
    // graph now lacks (as 0M links). A link is one key read from either end: a+ -> b+ and b- -> a- are the
    // same link, and a path walking a stored link backwards must find it, not add its mirror image.
    auto edge_key = [](const std::string& a, char oa, const std::string& b, char ob) {
        const std::string fwd = a + oa + "\t" + b + ob;
        const std::string rev = b + static_cast<char>(ob == '+' ? '-' : '+') + "\t" + a +
                                static_cast<char>(oa == '+' ? '-' : '+');
        return std::min(fwd, rev);
    };
    std::vector<GfaEdge> kept;
    std::unordered_set<std::string> keys;
    for (const auto& e : model.edges) {
        if (removed.count(e.from) || removed.count(e.to)) continue;
        // The input's own duplicates too: a link written twice, or once in each reading, is one link.
        if (keys.insert(edge_key(e.from, e.from_orient, e.to, e.to_orient)).second) kept.push_back(e);
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

    // Twin merge, graph-wide, including the duplicated modules the POA rebuild skips, where strand-split
    // twins sit. It runs before the acceptance checks, so they cover it too.
    const TwinMergeResult tm = merge_twin_nodes(model, table_nodes);
    summary.twin_paths_walking_several = tm.paths_walking_several;
    summary.twins_merged = tm.merges.size();
    say("nodes: " + std::to_string(model.node_order.size()) + "; added: " +
        std::to_string(summary.nodes_added) + "; removed: " + std::to_string(summary.nodes_removed) +
        "; merged as twins: " + std::to_string(tm.merges.size()));
    say("twin nodes merged: " + std::to_string(tm.merges.size()) +
        (tm.paths_walking_several > 0 ? "; paths now walking a merged node twice: " +
                                            std::to_string(tm.paths_walking_several)
                                      : std::string()));

    // Acceptance: every path must still spell exactly what it spelled on the way in, and step only across
    // existing links. Checked before anything is written, so a violation cannot reach disk.
    {
        std::vector<std::string> broken;
        if (model.paths.size() != spelled_before.size()) {
            broken.push_back("path count " + std::to_string(spelled_before.size()) + " -> " +
                             std::to_string(model.paths.size()));
        } else {
            for (std::size_t i = 0; i < model.paths.size(); ++i) {
                const std::string now = spell(model, model.paths[i].steps);
                if (now != spelled_before[i] && broken.size() < 8) {
                    broken.push_back(gfa_path_name(model.paths[i]) + " (" +
                                     std::to_string(spelled_before[i].size()) + " bp -> " +
                                     std::to_string(now.size()) + " bp)");
                }
            }
        }
        if (!broken.empty())
            throw std::runtime_error("refine: the rewrite changed a haplotype's sequence, which it must "
                                     "never do: " + cli::join_with_comma(broken));

        // Every consecutive step pair must be joined by a link that exists, in the orientation walked.
        std::unordered_set<std::string> have;
        for (const GfaEdge& e : model.edges) {
            have.insert(e.from + e.from_orient + "\t" + e.to + e.to_orient);
            have.insert(e.to + static_cast<char>(e.to_orient == '+' ? '-' : '+') + "\t" + e.from +
                        static_cast<char>(e.from_orient == '+' ? '-' : '+'));
        }
        for (const GfaPath& p : model.paths) {
            for (std::size_t k = 1; k < p.steps.size(); ++k) {
                const std::string key = p.steps[k - 1].node_id +
                                        (p.steps[k - 1].reverse ? '-' : '+') + "\t" +
                                        p.steps[k].node_id + (p.steps[k].reverse ? '-' : '+');
                if (have.find(key) == have.end())
                    throw std::runtime_error("refine: path " + gfa_path_name(p) +
                                             " steps across a link the rebuilt graph does not have: " +
                                             key);
            }
        }
    }

    // Re-sort, re-snarl and write the output family, the same way panphorte does.
    GraphSortOptions sort_opts;
    sort_opts.reference_path = resolved_reference;
    sort_opts.flip = !options.no_flip;
    const GraphSortResult sort_result = sort_graph_reference(model, sort_opts);
    say("sorted along the reference: nodes turned around: " + std::to_string(sort_result.flipped_nodes) +
        "; loops placed as one block: " + std::to_string(sort_result.cycle_blocks));

    // Closed, then checked: closing surfaces a write that failed on the way to storage, and some
    // platforms refuse to rename an open file at commit.
    rep.close();
    if (!rep) throw std::runtime_error("refine: failed writing " + options.out_prefix +
                                       ".refine.report.tsv");
    const std::string sorted_gfa = staged.stage(options.out_prefix + ".normalized.sorted.gfa");
    write_gfa_model(sorted_gfa, model);

    BubbleCallOptions bopts;
    bopts.reference_path = sort_result.resolved_reference;
    // The caller's threshold, not BubbleCallOptions' default, so an input built with another one keeps
    // the bubbles refinement never touched.
    bopts.min_variant_bp = options.resnarl_min_variant_bp;
    bopts.merge_nearby_bp = options.resnarl_merge_nearby_bp;
    bopts.fuse_ambiguous = options.resnarl_fuse_ambiguous;
    bopts.snarl_pairs_override = find_top_level_snarls_cactus(snarl_input_from_model(model), &bopts.snarl_inward_override);
    bopts.snarl_source_supplied = true;
    bopts.quiet = options.quiet;
    bopts.log_module = "refine";
    ParseGfaOptions parse_options;
    parse_options.include_paths = true;
    parse_options.include_sequences = true;
    const Graph sorted_graph = parse_gfa(sorted_gfa, parse_options);
    const BubbleCallReport report = call_bubbles_report(sorted_graph, bopts);
    write_bubbles_csv(staged.stage(options.out_prefix + ".bubbles.csv"), report.bubbles);
    write_bandage_node_colors_csv(staged.stage(options.out_prefix + ".bandage_nodes.csv"),
                                  report.bubbles, report.non_snp_bubbles);
    if (!options.gtf_path.empty()) {
        if (emit_gene_annotation(sorted_graph, sort_result.resolved_reference, options.gtf_path,
                                 staged.stage(options.out_prefix + ".bandage_genes.csv"), "", "refine")) {
            summary.wrote_gene_annotation = true;   // listed among the outputs
        }
    }
    summary.bubbles_after = report.bubbles.size();
    say("sites after refining: " +
        std::to_string(std::count_if(report.bubbles.begin(), report.bubbles.end(),
                                     [](const Bubble& b) { return !b.inversion_parent; })) +
        " (ids reassigned)");

    // panphorte's repeat tables, carried forward under refine's prefix with the output graph's ids: the
    // re-sort renumbers every node, and the tables are what name a folded repeat's nodes (call counts a
    // unit graph's copies at the anchor its table names). The repeat nodes were held fixed and kept out
    // of twin merge, so each still exists (a one-node unit graph the sort flipped spells its unit
    // reverse-complemented, which the unit check reads on either strand).
    std::set<std::string> written;
    for (const RepTable& t : rep_tables) {
        const std::string suffix = t.path.substr(options.panphorte_prefix.size());
        const std::string out_path = options.out_prefix + suffix;
        std::ofstream out(staged.stage(out_path));
        if (!out) throw std::runtime_error("refine: cannot write " + out_path);
        for (std::size_t c = 0; c < t.header.size(); ++c) out << (c ? "\t" : "") << t.header[c];
        out << '\n';
        for (std::vector<std::string> row : t.rows) {
            for (std::size_t c = 0; c < t.header.size(); ++c) {
                if (!kUnitIdColumns.count(t.header[c]) || row[c] == ".") continue;
                std::string ids;
                for (const std::string& n : split_list(row[c])) {
                    const auto it = sort_result.id_remap.find(n);
                    const std::string id = it == sort_result.id_remap.end() ? n : it->second;
                    if (model.seq.find(id) == model.seq.end())
                        throw std::runtime_error("refine: repeat node " + n + " of " + t.path +
                                                 " is missing from the refined graph");
                    ids += (ids.empty() ? "" : ";") + id;
                }
                row[c] = ids;
            }
            for (std::size_t c = 0; c < row.size(); ++c) out << (c ? "\t" : "") << row[c];
            out << '\n';
        }
        out.close();
        if (!out) throw std::runtime_error("refine: failed writing " + out_path);
        written.insert(out_path);
        summary.rep_tables_written.push_back(out_path);
    }
    if (summary.wrote_gene_annotation) written.insert(options.out_prefix + ".bandage_genes.csv");

    // Everything succeeded: move the family into place.
    staged.commit();

    // Optional family files this run did not write (a gene table, repeat tables, or the obsolete REP
    // provenance table of an earlier run at this prefix) describe another graph, and `call -b <prefix>`
    // would read them beside this one. An input is never removed.
    // The merge count of twin nodes is in the log; a twin table is not part of the family.
    for (const char* suffix : {".bandage_genes.csv", ".panphorte.unit_graph.tsv", ".panphorte.rep_provenance.tsv",
                               ".refine.twins.tsv"}) {
        const std::string path = options.out_prefix + suffix;
        if (written.count(path) != 0) continue;
        std::error_code ec;
        bool is_input = false;
        for (const std::string& in : {options.gfa_path, options.bubbles_csv_in, options.gtf_path})
            is_input = is_input || (!in.empty() && std::filesystem::equivalent(path, in, ec));
        for (const RepTable& t : rep_tables) is_input = is_input || std::filesystem::equivalent(path, t.path, ec);
        if (!is_input && std::filesystem::remove(path, ec)) summary.stale_removed.push_back(path);
    }
    return summary;
}

}  // namespace panvar
