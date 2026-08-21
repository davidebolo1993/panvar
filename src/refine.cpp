#include "panvar/refine.hpp"

#include <filesystem>
#include <fstream>

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
    bool reversed = false;   // the path crosses sink->source; steps above are canonicalized
};
// The interior a path carries between two anchors.
//
// Taking the FIRST occurrence of each anchor independently is wrong wherever an anchor repeats: two
// first occurrences need not bound the same traversal, and on a reverse crossing they arrive in the
// opposite order, so the interior reads backwards. Every valid (a, b) pair is enumerated instead and
// scored by how much of the region's declared interior it contains. A reverse crossing is canonicalized
// by reversing the steps and toggling orientation, while `lo`/`hi` keep the ORIGINAL indices and
// `reversed` records the direction, so the edit can be spliced back where it came from.
//
// If two different pairs are equally good the traversal is genuinely ambiguous, and choosing one
// arbitrarily would rewrite a copy the caller did not mean. That is refused.
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
        for (std::size_t k = lo + 1; k < hi; ++k) {
            if (region_interior.empty() || region_interior.count(steps[k].node_id) != 0) ++score;
        }
        cands.push_back({lo, hi, score, hi - lo, reversed});
    };
    for (const std::size_t ia : at_a)
        for (const std::size_t ib : at_b) {
            if (ib > ia) consider(ia, ib, false);        // a ... b: forward
            else if (ia > ib) consider(ib, ia, true);    // b ... a: the path crosses in reverse
        }
    if (cands.empty()) return std::nullopt;

    // Most declared interior wins; among equals the tightest span, which is the enclosing traversal
    // rather than one spanning several copies.
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
    out.steps.assign(steps.begin() + static_cast<long>(best.lo) + 1,
                     steps.begin() + static_cast<long>(best.hi));
    if (best.reversed) {
        std::reverse(out.steps.begin(), out.steps.end());
        for (PathStep& st : out.steps) st.reverse = !st.reverse;
    }
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
    // A path that touches this region's interior without spanning both anchors is a PARTIAL traversal.
    // It cannot be rewritten, so its old nodes are retained -- and retaining them retains the old EDGES
    // between them, leaving the pre-refinement topology beside the refined one where a walk can still
    // take it. Sequence losslessness cannot see that: every path still spells the same bases. Skipping
    // the region is the conservative default; --partial-path-policy retain is experimental.
    std::size_t partial_paths = 0;
    {
        std::unordered_set<std::string> interior_set(anch.interior.begin(), anch.interior.end());
        for (std::size_t idx = 0; idx < model.paths.size(); ++idx) {
            if (interior_between(model.paths[idx].steps, anch.a, anch.b, interior_set, nullptr))
                continue;   // spans both anchors
            for (const PathStep& st : model.paths[idx].steps) {
                if (interior_set.count(st.node_id) != 0) { ++partial_paths; break; }
            }
        }
    }
    if (partial_paths > 0 && opt.partial_path_policy_skip) {
        note = tag + ": skip (" + std::to_string(partial_paths) +
               " path(s) traverse the interior without spanning both anchors; "
               "--partial-path-policy retain to rebuild anyway)";
        return std::nullopt;
    }

    const std::unordered_set<std::string> region_interior(anch.interior.begin(), anch.interior.end());
    for (std::size_t idx : ref_and_paths) {
        bool ambiguous = false;
        auto got = interior_between(model.paths[idx].steps, anch.a, anch.b, region_interior, &ambiguous);
        if (ambiguous) {
            note = tag + ": skip (ambiguous_anchor_traversal: path " + gfa_path_name(model.paths[idx]) +
                   " has two equally good traversals between the anchors)";
            return std::nullopt;
        }
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
        // abPOA is handed the DISTINCT sequences, so every guard is measured over those. Summing bases
        // across all carriers instead made the decision depend on cohort composition: duplicating an
        // identical haplotype could turn "rebuilt" into "skipped" without changing a single byte of
        // POA input.
        std::set<std::string> distinct;
        for (const auto& kv : split) distinct.insert(spell(model, kv.second.segments[pos]));
        const std::size_t carriers = split.size();

        std::size_t longest = 0, distinct_total = 0;
        for (const std::string& d : distinct) {
            longest = std::max(longest, d.size());
            distinct_total += d.size();
        }
        const std::string sizes = "carriers " + std::to_string(carriers) + ", distinct " +
                                  std::to_string(distinct.size()) + ", longest " +
                                  std::to_string(longest) + " bp, distinct total " +
                                  std::to_string(distinct_total) + " bp";

        // The median let a single very long outlier into abPOA no matter how large it was; POA cost is
        // driven by the longest sequence, and the median is blind to it.
        if (longest > opt.max_poa_bp) {
            note = tag + ": skip (segment " + std::to_string(pos) + " longest " +
                   std::to_string(longest) + " bp > max-poa-bp " + std::to_string(opt.max_poa_bp) +
                   "; " + sizes + ")";
            return std::nullopt;
        }
        // Estimated work, bounded independently. With longest <= max_poa_bp and distinct <= max_walks
        // a product test against max_poa_bp^2 * max_walks is mathematically redundant, so it is a real
        // separate budget or it is nothing.
        if (opt.max_poa_work > 0 &&
            static_cast<double>(longest) * static_cast<double>(distinct_total) >
                static_cast<double>(opt.max_poa_work)) {
            note = tag + ": skip (segment " + std::to_string(pos) + " estimated POA work " +
                   std::to_string(static_cast<unsigned long long>(longest) * distinct_total) +
                   " cells > max-poa-work " + std::to_string(opt.max_poa_work) + "; " + sizes + ")";
            return std::nullopt;
        }
        if (distinct.size() > opt.max_walks) {
            note = tag + ": skip (segment " + std::to_string(pos) + " too diverse; " + sizes + ")";
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
        // The interior was canonicalized to reference direction for alignment; a path that crosses the
        // region in reverse must get it back in ITS direction, or the splice writes a
        // reverse-complemented interior into a forward walk. Same length, different sequence -- which
        // is precisely what the losslessness check exists to catch, and did.
        if (kv.second.reversed) {
            std::reverse(rebuilt.begin(), rebuilt.end());
            for (PathStep& st : rebuilt) st.reverse = !st.reverse;
        }
        edit.per_path[idx] = std::make_tuple(kv.second.lo, kv.second.hi, std::move(rebuilt));
    }

    for (const auto& n : anch.interior) if (!rep_set.count(n)) edit.interior_rm.insert(n);
    note = tag + ": partial_paths=" + std::to_string(partial_paths) + " REBUILT [" + (rep_set.empty() ? "plain" : "rep-residual") + "] " +
           std::to_string(ref_skel.size()) + " REP block(s), " + std::to_string(nseg) + " residual seg(s) -> +" +
           std::to_string(edit.added.size()) + " nodes";
    return edit;
}

}  // namespace

RefineSummary refine_graph(const RefineOptions& options) {
    RefineSummary summary;
    // No output may name an input: the refined GFA is written before re-snarl, CSV, Bandage and GTF
    // work, so an aliased path destroys a file still being read.
    {
        const std::string outs[] = {
            options.out_prefix + ".normalized.sorted.gfa",
            options.out_prefix + ".bubbles.csv",
            options.out_prefix + ".bandage_nodes.csv",
            options.out_prefix + ".bandage_genes.csv",
            options.out_prefix + ".refine.report.tsv",
        };
        for (const std::string& in : {options.gfa_path, options.bubbles_csv_in, options.gtf_path}) {
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

    GfaModel model = read_gfa_model(options.gfa_path);
    std::vector<Bubble> bubbles = read_bubbles_csv(options.bubbles_csv_in);
    std::string resolved_reference;

    // The same graph contract every other module applies. refine REWRITES the graph, so accepting a
    // malformed one means emitting a repaired-looking graph that was never validated. The Graph view is
    // also what resolves the reference, so the rule is shared rather than reimplemented here.
    {
        ParseGfaOptions popts;
        popts.include_paths = true;
        popts.include_sequences = true;
        const Graph gview = parse_gfa(options.gfa_path, popts);
        validate_graph_paths(gview, "refine", true, true);
        if (gview.paths.empty())
            throw std::runtime_error("refine: input GFA has no P/W paths");
        resolved_reference = resolve_reference_path_name(gview, options.reference_path, "refine");

        // The CSV and the graph are separate inputs with nothing tying them together.
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

    const std::unordered_set<std::string> self_loop = self_loop_nodes(model);

    // The resolved name, matched exactly. Selecting the first case-insensitive SUBSTRING gave no
    // priority to an exact match and never rejected ambiguity, so `-r hap1` could silently pick
    // `hap10` if it came first in the file.
    std::size_t ref_idx = model.paths.size();
    for (std::size_t i = 0; i < model.paths.size(); ++i) {
        if (gfa_path_name(model.paths[i]) == resolved_reference) { ref_idx = i; break; }
    }
    if (ref_idx == model.paths.size()) {
        throw std::runtime_error("refine: resolved reference '" + resolved_reference +
                                 "' is not a path of the model");
    }

    // Every input path's spelled sequence, kept so the rewrite can be PROVEN lossless rather than
    // assumed to be. refine re-aligns interiors and rebuilds nodes; "sequence-preserving" was an
    // algorithmic argument with nothing checking it.
    // Indexed by path position, not keyed by name: a name-keyed map silently tolerates a naming rule
    // that does not round-trip, which is exactly the bug this replaced.
    std::vector<std::string> spelled_before;
    spelled_before.reserve(model.paths.size());
    for (const GfaPath& p : model.paths) spelled_before.push_back(spell(model, p.steps));

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
        std::unordered_set<std::string> found;
        for (const auto& b : bubbles) {
            const std::string id = std::to_string(b.id);
            if (want.count(id)) { picked.push_back(b); found.insert(id); }
        }
        // A --bubble-id naming nothing used to rebuild zero regions and still write a full, unchanged
        // output family with exit 0 -- indistinguishable from "nothing needed refining".
        std::vector<std::string> absent;
        for (const std::string& w : want) if (!found.count(w)) absent.push_back(w);
        if (!absent.empty()) {
            std::sort(absent.begin(), absent.end());
            throw std::runtime_error("refine: --bubble-id not present in the bubbles CSV: " +
                                     cli::join_with_comma(absent));
        }
        // Selected ids go through the SAME component grouping as auto mode. Forcing them into one
        // region assumed they were adjacent: two disjoint bubbles became a single region spanning
        // disconnected anchors, which then failed as one unit -- `--bubble-id 1,2` reported "rebuilt 0,
        // skipped 1" where auto mode rebuilt both.
        comps = components(picked);
    }

    long long next_id = 0;
    for (const auto& kv : model.seq) {
        char* end = nullptr;
        const long long v = std::strtoll(kv.first.c_str(), &end, 10);
        if (end != kv.first.c_str() && *end == '\0' && v > next_id) next_id = v;
    }
    ++next_id;

    // Component order was whatever the container produced, which makes the run's node numbering (and
    // so its output bytes) depend on hash iteration order. Sorted by first anchor, so two runs of the
    // same input agree byte for byte.
    std::sort(comps.begin(), comps.end(), [](const std::vector<Bubble>& a, const std::vector<Bubble>& b) {
        if (a.empty() || b.empty()) return b.empty() < a.empty();
        if (a.front().source != b.front().source) return a.front().source < b.front().source;
        return a.front().sink < b.front().sink;
    });

    // The reasons a region was skipped were composed and thrown away, so "skipped 3" was the whole
    // story a user got.
    // Nothing lands in its final location until the whole run has succeeded. The report used to be
    // opened before the losslessness check, so "a failed acceptance writes nothing" was false, and a
    // later re-snarl or GTF failure left a usable-looking subset of the family.
    cli::StagedOutputs staged("refine");
    const std::string report_path = staged.stage(options.out_prefix + ".refine.report.tsv");
    std::ofstream rep(report_path);
    if (!rep) throw std::runtime_error("refine: cannot write " + options.out_prefix + ".refine.report.tsv");
    rep << "region\tn_bubbles\tsource\tsink\tdecision\treason\n";

    std::vector<RegionEdit> edits;
    cli::ProgressBar progress(options.quiet ? "" : "Refining regions", comps.size());
    std::size_t region_no = 0;
    for (const auto& comp : comps) {
        std::string note;
        auto edit = process_region(comp, model, path_order, ref_idx, self_loop, options, next_id, note);
        ++region_no;
        rep << region_no << '\t' << comp.size() << '\t'
            << (comp.empty() ? std::string("-") : comp.front().source) << '\t'
            << (comp.empty() ? std::string("-") : comp.back().sink) << '\t'
            << (edit ? "rebuilt" : "skipped") << '\t' << (note.empty() ? "-" : note) << '\n';
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

    // ACCEPTANCE: every path must still spell exactly what it spelled on the way in. refine re-aligns
    // interiors and rebuilds nodes, and losslessness was an argument about the algorithm with nothing
    // testing it. Checked before anything is written, so a violation cannot reach disk.
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

    // re-sort, re-snarl and write the output family, the same way panphorte does
    GraphSortOptions sort_opts;
    sort_opts.reference_path = resolved_reference;
    sort_opts.flip = !options.no_flip;
    const GraphSortResult sort_result = sort_graph_reference(model, sort_opts);

    // Closed, then checked. flush() alone leaves the stream open across the commit below -- which
    // works on Unix and is not portable to a platform that refuses to rename an open file -- and it
    // does not turn a write that failed on the way to storage into an error, so a truncated report
    // could be committed as part of a successful family.
    rep.close();
    if (!rep) throw std::runtime_error("refine: failed writing " + options.out_prefix +
                                       ".refine.report.tsv");
    const std::string sorted_gfa = staged.stage(options.out_prefix + ".normalized.sorted.gfa");
    write_gfa_model(sorted_gfa, model);

    BubbleCallOptions bopts;
    bopts.reference_path = sort_result.resolved_reference;
    // The re-snarl silently used BubbleCallOptions' own 50 bp default, so an input built with a
    // different threshold could lose bubbles that refinement never touched.
    bopts.min_variant_bp = options.resnarl_min_variant_bp;
    bopts.snarl_pairs_override = find_top_level_snarls_cactus(snarl_input_from_model(model));
    bopts.snarl_source_supplied = true;
    bopts.quiet = options.quiet;
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
                                 staged.stage(options.out_prefix + ".bandage_genes.csv"))) {
            summary.wrote_gene_annotation = true;   // it was omitted from the reported output list
        }
    }
    summary.bubbles_after = report.bubbles.size();
    // Everything succeeded: move the family into place.
    staged.commit();
    return summary;
}

}  // namespace panvar
