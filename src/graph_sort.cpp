#include "panvar/graph_sort.hpp"

#include "panvar/gfa.hpp"
#include "panvar/graph_utils.hpp" // reverse_complement
#include "panvar/ref_path.hpp"

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <limits>
#include <queue>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace panvar {

std::string gfa_path_name(const GfaPath& p) {
    if (p.type == 'W') {
        return w_path_name(p.sample, p.hap, p.seqid, p.start, p.end);
    }
    return p.name;
}

namespace {

// Resolve the reference path as every module does (match_reference_name).
std::size_t resolve_reference(const GfaModel& model, const std::string& query) {
    std::vector<std::string> names;
    for (const GfaPath& p : model.paths) names.push_back(gfa_path_name(p));
    const std::vector<std::size_t> hits = match_reference_name(names, query);
    if (hits.size() == 1) return hits[0];
    if (hits.empty()) {
        // A W line named by its own (0-based) coordinates, which panvar names from start + 1.
        const std::string w = w_literal_as_path_name(query);
        for (std::size_t i = 0; !w.empty() && i < model.paths.size(); ++i)
            if (model.paths[i].type == 'W' && gfa_path_name(model.paths[i]) == w) return i;
        throw std::runtime_error("Reference path not found: " + query);
    }
    std::string msg = "Reference path ambiguous: " + query + " matches:";
    for (std::size_t i : hits) msg += "\n  " + gfa_path_name(model.paths[i]);
    throw std::runtime_error(msg);
}

char toggle(char o) { return o == '+' ? '-' : '+'; }

// Numeric value of an id for deterministic tie-breaking; non-numeric -> max.
std::uint64_t id_num(const std::string& id) {
    if (id.empty()) return std::numeric_limits<std::uint64_t>::max();
    for (char c : id) {
        if (!std::isdigit(static_cast<unsigned char>(c))) {
            return std::numeric_limits<std::uint64_t>::max();
        }
    }
    return std::strtoull(id.c_str(), nullptr, 10);
}

} // namespace

GraphSortResult sort_graph_reference(GfaModel& model, const GraphSortOptions& options) {
    GraphSortResult result;
    if (options.reference_path.empty()) {
        throw std::runtime_error("sort_graph_reference: reference_path is required");
    }
    const std::size_t ref_idx = resolve_reference(model, options.reference_path);
    const GfaPath& ref = model.paths[ref_idx];
    result.resolved_reference = gfa_path_name(ref);
    result.reference_steps = ref.steps.size();

    // Reference first-occurrence rank + orientation per node.
    std::unordered_map<std::string, std::size_t> ref_rank;
    std::unordered_map<std::string, bool> ref_reverse;
    for (std::size_t i = 0; i < ref.steps.size(); ++i) {
        const PathStep& s = ref.steps[i];
        if (ref_rank.emplace(s.node_id, i).second) {
            ref_reverse[s.node_id] = s.reverse;
        }
    }

    // ---- 1. flip nodes the reference traverses in reverse (preserves path spellings) ---
    if (options.flip) {
        std::unordered_set<std::string> flip_set;
        for (const auto& [id, rev] : ref_reverse) {
            if (rev) flip_set.insert(id);
        }
        if (!flip_set.empty()) {
            for (const std::string& id : flip_set) {
                auto it = model.seq.find(id);
                if (it != model.seq.end()) it->second = reverse_complement(it->second);
            }
            for (GfaEdge& e : model.edges) {
                if (flip_set.count(e.from)) e.from_orient = toggle(e.from_orient);
                if (flip_set.count(e.to)) e.to_orient = toggle(e.to_orient);
            }
            for (GfaPath& p : model.paths) {
                for (PathStep& s : p.steps) {
                    if (flip_set.count(s.node_id)) s.reverse = !s.reverse;
                }
            }
            result.flipped_nodes = flip_set.size();
        }
    }

    // ---- 2. reference-guided topological order -----------------------------------------
    const std::size_t INF = std::numeric_limits<std::size_t>::max();
    // Bucket each node by the rank of the nearest reference node (undirected BFS), so a
    // bubble interior shares its source's bucket and lands contiguously.
    std::unordered_map<std::string, std::size_t> bucket;
    {
        std::unordered_map<std::string, std::vector<std::string>> adj;
        adj.reserve(model.node_order.size() * 2);
        for (const GfaEdge& e : model.edges) {
            adj[e.from].push_back(e.to);
            adj[e.to].push_back(e.from);
        }
        std::vector<std::pair<std::size_t, std::string>> seeds;
        for (const auto& [id, r] : ref_rank) seeds.emplace_back(r, id);
        std::sort(seeds.begin(), seeds.end());
        std::queue<std::string> bfs;
        for (const auto& [r, id] : seeds) {
            if (bucket.emplace(id, r).second) bfs.push(id);
        }
        while (!bfs.empty()) {
            const std::string u = bfs.front();
            bfs.pop();
            const std::size_t bu = bucket[u];
            const auto ai = adj.find(u);
            if (ai == adj.end()) continue;
            for (const std::string& v : ai->second) {
                if (bucket.emplace(v, bu).second) bfs.push(v);
            }
        }
        for (const std::string& id : model.node_order) bucket.emplace(id, INF);
    }

    // Directed arcs from all path walks, in WALK direction. A step pair `a+ b+` says a precedes b, and
    // `b- a-` says the same thing walked backwards; reading the second as "b before a" would turn every
    // node a reverse-stored haplotype covers into a 2-cycle against the reference. A pair that switches
    // orientation takes its direction from whichever node is on the reference: after the flip above a
    // reference node reads forward along the reference, so its step orientation says which way the path
    // is moving there. That gives an arc to an off-reference allele node stored opposite to its flanks
    // (in `5- 1+ 2-`, 1 lies between 2 and 5), so it is numbered inside its own site.
    // Two reference nodes in opposite orientations are a genuine inversion and state no order.
    const std::size_t n = model.node_order.size();
    std::unordered_map<std::string, std::size_t> idx;
    idx.reserve(n * 2);
    for (std::size_t i = 0; i < n; ++i) idx.emplace(model.node_order[i], i);
    std::vector<std::vector<std::size_t>> succ(n);
    {
        std::unordered_set<std::uint64_t> seen_arc;
        for (const GfaPath& p : model.paths) {
            for (std::size_t i = 1; i < p.steps.size(); ++i) {
                const PathStep& s0 = p.steps[i - 1];
                const PathStep& s1 = p.steps[i];
                if (s0.node_id == s1.node_id) continue;          // self-arc (tandem repeat)
                const auto ia = idx.find(s0.node_id), ib = idx.find(s1.node_id);
                if (ia == idx.end() || ib == idx.end()) continue;
                // Orientation relative to the reference's reading of the node. After the flip the
                // reference reads every one of its nodes forward; with --no-flip it still reads some
                // reversed, which are numbered along the reference only when read this way.
                const bool rev0 = s0.reverse != (!options.flip && ref_reverse.count(s0.node_id) && ref_reverse.at(s0.node_id));
                const bool rev1 = s1.reverse != (!options.flip && ref_reverse.count(s1.node_id) && ref_reverse.at(s1.node_id));
                bool backwards = rev0;                           // walking against the reference
                if (rev0 != rev1) {
                    const bool r0 = ref_rank.count(s0.node_id) != 0, r1 = ref_rank.count(s1.node_id) != 0;
                    if (r0 == r1) continue;                      // inversion (or no anchor): no order
                    backwards = r0 ? rev0 : rev1;
                }
                const std::size_t from = backwards ? ib->second : ia->second;
                const std::size_t to = backwards ? ia->second : ib->second;
                if (seen_arc.insert((static_cast<std::uint64_t>(from) << 32) | to).second) succ[from].push_back(to);
            }
        }
    }

    struct Key {
        std::size_t bucket;
        std::uint64_t num;
        std::string id;
        bool operator>(const Key& o) const {
            if (bucket != o.bucket) return bucket > o.bucket;
            if (num != o.num) return num > o.num;
            return id > o.id;
        }
        bool operator<(const Key& o) const { return o > *this; }
    };
    std::vector<Key> key(n);
    for (std::size_t i = 0; i < n; ++i) key[i] = Key{bucket[model.node_order[i]], id_num(model.node_order[i]), model.node_order[i]};

    // Nodes on a cycle (a tandem repeat, a reference that revisits a region) cannot be ordered among
    // themselves by the walks, so each strongly connected component is placed as one block: blocks in
    // topological order, the smallest key first among those ready, and the nodes inside a block by key.
    // Breaking a cycle only when nothing else is ready would hold a stalled node back until every node
    // downstream of it is numbered, placing a repeat block after the unique flank that follows it.
    std::vector<std::size_t> comp(n, SIZE_MAX);
    std::size_t n_comp = 0;
    {   // Tarjan, iterative
        std::vector<std::size_t> index_of(n, SIZE_MAX), low(n, 0), stack;
        std::vector<char> on_stack(n, 0);
        std::size_t counter = 0;
        struct Frame { std::size_t v, next; };
        for (std::size_t root = 0; root < n; ++root) {
            if (index_of[root] != SIZE_MAX) continue;
            std::vector<Frame> call{{root, 0}};
            index_of[root] = low[root] = counter++;
            stack.push_back(root);
            on_stack[root] = 1;
            while (!call.empty()) {
                Frame& f = call.back();
                if (f.next < succ[f.v].size()) {
                    const std::size_t w = succ[f.v][f.next++];
                    if (index_of[w] == SIZE_MAX) {
                        index_of[w] = low[w] = counter++;
                        stack.push_back(w);
                        on_stack[w] = 1;
                        call.push_back({w, 0});
                    } else if (on_stack[w]) {
                        low[f.v] = std::min(low[f.v], index_of[w]);
                    }
                    continue;
                }
                const std::size_t v = f.v;
                if (low[v] == index_of[v]) {
                    std::size_t w;
                    do {
                        w = stack.back();
                        stack.pop_back();
                        on_stack[w] = 0;
                        comp[w] = n_comp;
                    } while (w != v);
                    ++n_comp;
                }
                call.pop_back();
                if (!call.empty()) low[call.back().v] = std::min(low[call.back().v], low[v]);
            }
        }
    }
    std::vector<std::vector<std::size_t>> members(n_comp);
    for (std::size_t i = 0; i < n; ++i) members[comp[i]].push_back(i);
    std::vector<Key> comp_key(n_comp);
    for (std::size_t c = 0; c < n_comp; ++c) {
        std::sort(members[c].begin(), members[c].end(), [&](std::size_t a, std::size_t b) { return key[a] < key[b]; });
        comp_key[c] = key[members[c].front()];
    }
    std::vector<std::vector<std::size_t>> comp_succ(n_comp);
    std::vector<std::size_t> comp_indeg(n_comp, 0);
    {
        std::unordered_set<std::uint64_t> seen;
        for (std::size_t v = 0; v < n; ++v)
            for (const std::size_t w : succ[v]) {
                const std::size_t a = comp[v], b = comp[w];
                if (a == b || !seen.insert((static_cast<std::uint64_t>(a) << 32) | b).second) continue;
                comp_succ[a].push_back(b);
                ++comp_indeg[b];
            }
    }
    using CompItem = std::pair<Key, std::size_t>;
    auto later = [](const CompItem& x, const CompItem& y) { return x.first > y.first; };
    std::priority_queue<CompItem, std::vector<CompItem>, decltype(later)> ready(later);
    for (std::size_t c = 0; c < n_comp; ++c) if (comp_indeg[c] == 0) ready.push({comp_key[c], c});
    std::vector<std::string> order;
    order.reserve(n);
    while (!ready.empty()) {
        const std::size_t c = ready.top().second;
        ready.pop();
        for (const std::size_t v : members[c]) order.push_back(model.node_order[v]);
        for (const std::size_t d : comp_succ[c]) if (--comp_indeg[d] == 0) ready.push({comp_key[d], d});
    }
    result.cycle_blocks = 0;
    for (std::size_t c = 0; c < n_comp; ++c) result.cycle_blocks += members[c].size() > 1;

    // ---- 3. renumber 1..N in sorted order ---------------------------------------------
    if (options.renumber) {
        std::unordered_map<std::string, std::string>& remap = result.id_remap;
        remap.reserve(order.size() * 2);
        for (std::size_t i = 0; i < order.size(); ++i) remap[order[i]] = std::to_string(i + 1);
        std::unordered_map<std::string, std::string> new_seq;
        new_seq.reserve(model.seq.size() * 2);
        for (const auto& [old_id, s] : model.seq) {
            const auto it = remap.find(old_id);
            new_seq[it != remap.end() ? it->second : old_id] = s;
        }
        model.seq = std::move(new_seq);
        model.node_order.clear();
        model.node_order.reserve(order.size());
        for (std::size_t i = 0; i < order.size(); ++i) model.node_order.push_back(std::to_string(i + 1));
        auto rm = [&](const std::string& id) -> std::string {
            const auto it = remap.find(id);
            return it != remap.end() ? it->second : id;
        };
        for (GfaEdge& e : model.edges) {
            e.from = rm(e.from);
            e.to = rm(e.to);
        }
        for (GfaPath& p : model.paths) {
            for (PathStep& s : p.steps) s.node_id = rm(s.node_id);
        }
    } else {
        model.node_order = order; // reorder only, keep ids
    }
    return result;
}

} // namespace panvar
