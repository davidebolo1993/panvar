#include "panvar/graph_sort.hpp"

#include "panvar/graph_utils.hpp" // reverse_complement

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <limits>
#include <queue>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace panvar {

std::string gfa_path_name(const GfaPath& p) {
    if (p.type == 'W') {
        return p.sample + "#" + p.hap + "#" + p.seqid + ":" + p.start + "-" + p.end;
    }
    return p.name;
}

namespace {

std::string to_lower(std::string s) {
    for (char& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}

// Resolve the reference path: exact name first, else a unique case-insensitive substring.
std::size_t resolve_reference(const GfaModel& model, const std::string& query) {
    for (std::size_t i = 0; i < model.paths.size(); ++i) {
        if (gfa_path_name(model.paths[i]) == query) return i;
    }
    const std::string q = to_lower(query);
    std::vector<std::size_t> hits;
    for (std::size_t i = 0; i < model.paths.size(); ++i) {
        if (to_lower(gfa_path_name(model.paths[i])).find(q) != std::string::npos) hits.push_back(i);
    }
    if (hits.size() == 1) return hits[0];
    if (hits.empty()) throw std::runtime_error("Reference path not found: " + query);
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

    // Directed arcs from all path walks (consecutive steps), for Kahn's topological order.
    std::unordered_map<std::string, std::unordered_set<std::string>> succ;
    std::unordered_map<std::string, std::size_t> indeg;
    indeg.reserve(model.node_order.size() * 2);
    for (const std::string& id : model.node_order) indeg[id] = 0;
    for (const GfaPath& p : model.paths) {
        for (std::size_t i = 1; i < p.steps.size(); ++i) {
            const std::string& a = p.steps[i - 1].node_id;
            const std::string& b = p.steps[i].node_id;
            if (a == b) continue; // ignore self-arcs (tandem repeats)
            if (!model.has_node(a) || !model.has_node(b)) continue;
            if (succ[a].insert(b).second) indeg[b] += 1;
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
    };
    auto make_key = [&](const std::string& id) { return Key{bucket[id], id_num(id), id}; };
    std::priority_queue<Key, std::vector<Key>, std::greater<Key>> ready;
    for (const std::string& id : model.node_order) {
        if (indeg[id] == 0) ready.push(make_key(id));
    }

    std::vector<std::string> order;
    order.reserve(model.node_order.size());
    std::unordered_set<std::string> emitted;
    emitted.reserve(model.node_order.size() * 2);
    auto emit = [&](const std::string& u) {
        order.push_back(u);
        emitted.insert(u);
        const auto si = succ.find(u);
        if (si == succ.end()) return;
        for (const std::string& v : si->second) {
            if (indeg[v] > 0 && --indeg[v] == 0) ready.push(make_key(v));
        }
    };
    while (order.size() < model.node_order.size()) {
        if (!ready.empty()) {
            const Key k = ready.top();
            ready.pop();
            if (emitted.count(k.id)) continue;
            emit(k.id);
        } else {
            // cycle: break by emitting the remaining node with the smallest key.
            std::string best;
            Key bk{INF, std::numeric_limits<std::uint64_t>::max(), ""};
            bool found = false;
            for (const std::string& id : model.node_order) {
                if (emitted.count(id)) continue;
                const Key k = make_key(id);
                if (!found || bk > k) {
                    bk = k;
                    best = id;
                    found = true;
                }
            }
            if (!found) break;
            emit(best);
        }
    }

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
