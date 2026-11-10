#include "panvar/background_graph.hpp"

#include "panvar/background.hpp"
#include "panvar/gz_reader.hpp"
#include "panvar/kmer.hpp"

#include <algorithm>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <unordered_map>

namespace panvar::bg::graph {

namespace {

std::uint64_t hash_line(const std::string& s, std::uint64_t h) {
    std::size_t i = 0;
    for (; i + 8 <= s.size(); i += 8) { std::uint64_t w; std::memcpy(&w, s.data() + i, 8); h = kmer::mix64(h ^ w); }
    std::uint64_t w = 0;
    std::memcpy(&w, s.data() + i, s.size() - i);
    return kmer::mix64(h ^ w ^ (static_cast<std::uint64_t>(s.size()) << 56));
}

char complement(char c) {
    switch (c) {
        case 'A': return 'T'; case 'C': return 'G'; case 'G': return 'C'; case 'T': return 'A';
        case 'a': return 't'; case 'c': return 'g'; case 'g': return 'c'; case 't': return 'a';
        default: return 'N';
    }
}

std::vector<std::string> split(const std::string& s, char sep) {
    std::vector<std::string> out;
    std::size_t b = 0;
    for (std::size_t i = 0; i <= s.size(); ++i)
        if (i == s.size() || s[i] == sep) { out.push_back(s.substr(b, i - b)); b = i + 1; }
    return out;
}

// The graph's nodes (sequence and name) and a reader of its haplotype walks.
struct Graph {
    std::unordered_map<std::string, std::uint32_t> id;
    std::vector<std::uint64_t> start;   // into seq; size nodes + 1
    std::string seq;
    std::uint64_t checksum = 0x5bd1e9955bd1e995ULL;

    std::uint64_t len(std::uint32_t v) const { return start[v + 1] - start[v]; }
    // Base i of node v read in orientation `rev`.
    char base(std::uint32_t v, bool rev, std::uint64_t i) const {
        return rev ? complement(seq[start[v] + len(v) - 1 - i]) : seq[start[v] + i];
    }
    void append(std::string& out, std::uint32_t v, bool rev, std::uint64_t from, std::uint64_t to) const {
        for (std::uint64_t i = from; i < to; ++i) out += base(v, rev, i);
    }
};

Graph load_nodes(const std::string& gfa) {
    Graph g;
    GzLineReader in(gfa);
    if (!in.ok()) throw std::runtime_error("background: cannot read " + gfa);
    std::string line;
    g.start.push_back(0);
    while (in.getline(line)) {
        g.checksum = hash_line(line, g.checksum);
        if (line.size() > 2 && line[0] == 'L' && line[1] == '\t') {   // haplotypes are spelled by concatenation
            const auto f = split(line, '\t');
            if (f.size() > 5 && f[5] != "*" && f[5] != "0M" && !f[5].empty())
                throw std::runtime_error("background: " + gfa + ": a link with overlap " + f[5] + "; only 0M or * links are supported");
            continue;
        }
        if (line.size() < 2 || line[0] != 'S' || line[1] != '\t') continue;
        const auto t1 = line.find('\t', 2);
        const std::string name = line.substr(2, t1 == std::string::npos ? std::string::npos : t1 - 2);
        if (t1 == std::string::npos) throw std::runtime_error("background: " + gfa + ": an S line without sequence");
        const auto t2 = line.find('\t', t1 + 1);
        std::string s = line.substr(t1 + 1, t2 == std::string::npos ? std::string::npos : t2 - t1 - 1);
        if (s == "*") throw std::runtime_error("background: " + gfa + ": node " + name + " has no sequence");
        if (!g.id.emplace(name, static_cast<std::uint32_t>(g.start.size() - 1)).second)
            throw std::runtime_error("background: " + gfa + ": node " + name + " defined twice");
        g.seq += s;
        g.start.push_back(g.seq.size());
    }
    return g;
}

struct Walk {
    std::string contig;
    std::uint64_t offset = 0;                        // coordinate of its first base on the contig
    std::vector<std::pair<std::uint32_t, bool>> steps;   // node, reversed
};

// Calls f(walk) for every P and W line, in file order.
template <class F>
void for_each_walk(const std::string& gfa, const Graph& g, F&& f) {
    GzLineReader in(gfa);
    if (!in.ok()) throw std::runtime_error("background: cannot read " + gfa);
    std::string line;
    Walk w;
    auto node = [&](const std::string& name) {
        const auto it = g.id.find(name);
        if (it == g.id.end()) throw std::runtime_error("background: " + gfa + ": a walk visits undefined node " + name);
        return it->second;
    };
    while (in.getline(line)) {
        if (line.size() < 2 || line[1] != '\t' || (line[0] != 'P' && line[0] != 'W')) continue;
        const auto t = split(line, '\t');
        w.steps.clear();
        w.offset = 0;
        if (line[0] == 'P') {
            if (t.size() < 3) throw std::runtime_error("background: " + gfa + ": a short P line");
            w.contig = t[1];
            if (t.size() > 3 && t[3] != "*" && !t[3].empty())
                for (const auto& ov : split(t[3], ','))
                    if (ov != "0M" && ov != "*") throw std::runtime_error("background: " + gfa + ": path " + t[1] + " has overlap " + ov + "; only 0M is supported");
            for (const auto& s : split(t[2], ',')) {
                if (s.size() < 2 || (s.back() != '+' && s.back() != '-')) throw std::runtime_error("background: " + gfa + ": bad P step " + s);
                w.steps.emplace_back(node(s.substr(0, s.size() - 1)), s.back() == '-');
            }
        } else {
            if (t.size() < 7) throw std::runtime_error("background: " + gfa + ": a short W line");
            w.contig = t[1] + "#" + t[2] + "#" + t[3];
            w.offset = t[4] == "*" ? 0 : std::stoull(t[4]);
            const std::string& ws = t[6];
            for (std::size_t i = 0; i < ws.size();) {
                if (ws[i] != '>' && ws[i] != '<') throw std::runtime_error("background: " + gfa + ": bad W walk");
                std::size_t e = i + 1;
                while (e < ws.size() && ws[e] != '>' && ws[e] != '<') ++e;
                w.steps.emplace_back(node(ws.substr(i + 1, e - i - 1)), ws[i] == '<');
                i = e;
            }
        }
        f(w);
    }
}

// The junction at step i: the node tuple whose spelled string holds the k-mers that start in step i's node
// and end after it, and that string.
struct Junction {
    std::string key;     // (node << 1 | reversed) as 4-byte words
    std::uint64_t s = 0; // starts in node i: min(len_i, k-1)
};
void junction_key(const Graph& g, const Walk& w, std::size_t i, Junction& j) {
    const std::uint64_t k1 = kK - 1;
    j.key.clear();
    auto put = [&](std::uint32_t v, bool rev) { const std::uint32_t x = (v << 1) | (rev ? 1u : 0u); j.key.append(reinterpret_cast<const char*>(&x), 4); };
    put(w.steps[i].first, w.steps[i].second);
    j.s = std::min<std::uint64_t>(g.len(w.steps[i].first), k1);
    std::uint64_t have = 0;
    for (std::size_t t = i + 1; t < w.steps.size() && have < k1; ++t) {
        put(w.steps[t].first, w.steps[t].second);
        have += g.len(w.steps[t].first);
    }
}
std::string junction_string(const Graph& g, const std::string& key) {
    const std::uint64_t k1 = kK - 1;
    std::string out;
    std::uint32_t x;
    std::memcpy(&x, key.data(), 4);
    const std::uint32_t v = x >> 1;
    const bool rev = x & 1;
    const std::uint64_t l = g.len(v), s = std::min<std::uint64_t>(l, k1);
    g.append(out, v, rev, l - s, l);
    std::uint64_t need = k1;
    for (std::size_t o = 4; o < key.size() && need > 0; o += 4) {
        std::memcpy(&x, key.data() + o, 4);
        const std::uint64_t take = std::min<std::uint64_t>(g.len(x >> 1), need);
        g.append(out, x >> 1, x & 1, 0, take);
        need -= take;
    }
    return out;
}

std::uint64_t acgt_of(const Graph& g, std::uint32_t v) {
    std::uint64_t n = 0;
    for (std::uint64_t i = g.start[v]; i < g.start[v + 1]; ++i) n += kmer::base_code(g.seq[i]) >= 0;
    return n;
}

}  // namespace

GraphCounts count_graph(const std::string& gfa, std::uint32_t, std::size_t, const std::function<bool(std::string, std::uint64_t)>& push) {
    const Graph g = load_nodes(gfa);
    GraphCounts out;
    out.checksum = g.checksum;
    const std::size_t nodes = g.start.size() - 1;
    std::vector<std::uint64_t> traversals(nodes, 0);
    std::unordered_map<std::string, std::uint64_t> junctions;
    std::size_t junction_bytes = 0;
    constexpr std::size_t kJunctionBudget = std::size_t{1} << 30;
    bool stopped = false;
    auto flush_junctions = [&] {
        for (auto& [key, n] : junctions) {
            if (stopped) break;
            std::string s = junction_string(g, key);
            if (s.size() >= kK && !push(std::move(s), n)) stopped = true;
        }
        junctions.clear();
        junction_bytes = 0;
    };
    Junction j;
    for_each_walk(gfa, g, [&](const Walk& w) {
        if (stopped) return;
        out.contigs.push_back(w.contig);
        std::uint64_t len = w.offset;
        for (const auto& st : w.steps) len += g.len(st.first);
        out.contig_lengths.push_back(len);
        ++out.haplotypes;
        for (std::size_t i = 0; i < w.steps.size(); ++i) {
            ++traversals[w.steps[i].first];
            if (i + 1 == w.steps.size()) continue;   // nothing follows: its k-mers are internal or none
            junction_key(g, w, i, j);
            auto [it, fresh] = junctions.emplace(j.key, 0);
            ++it->second;
            if (fresh) junction_bytes += j.key.size() + 64;
            if (junction_bytes > kJunctionBudget) flush_junctions();
        }
    });
    if (!stopped) flush_junctions();
    for (std::uint32_t v = 0; v < nodes && !stopped; ++v) {
        if (traversals[v] == 0) continue;
        out.bases += acgt_of(g, v) * traversals[v];
        if (g.len(v) >= kK && !push(g.seq.substr(g.start[v], g.len(v)), traversals[v])) stopped = true;
    }
    return out;
}

std::vector<anchors::Hit> graph_anchors(const std::string& gfa, const anchors::AnchorLookup& lookup, std::uint32_t genome, std::size_t) {
    const Graph g = load_nodes(gfa);
    const std::size_t nodes = g.start.size() - 1;
    // anchors inside each node, on its forward strand: (offset, id, forward)
    struct NodeAnchor { std::uint64_t off; std::uint32_t id; bool fwd; };
    std::vector<std::vector<NodeAnchor>> inside(nodes);
    for (std::uint32_t v = 0; v < nodes; ++v) {
        std::vector<anchors::Hit> h;
        anchors::scan_anchors(lookup, g.seq.data() + g.start[v], g.len(v), 0, 0, 0, h);
        for (const auto& x : h) inside[v].push_back({x.pos, x.anchor, x.forward});
    }
    std::unordered_map<std::string, std::vector<NodeAnchor>> junction_memo;
    std::vector<anchors::Hit> out;
    std::uint32_t contig = 0;
    Junction j;
    for_each_walk(gfa, g, [&](const Walk& w) {
        std::uint64_t pos = w.offset;
        for (std::size_t i = 0; i < w.steps.size(); ++i) {
            const auto [v, rev] = w.steps[i];
            const std::uint64_t l = g.len(v);
            if (pos + l >= (1ULL << 32)) throw std::runtime_error("background: " + gfa + ": a haplotype of 4 Gb or more; anchors keep 32-bit positions");
            for (const auto& a : inside[v]) {
                const std::uint64_t p = rev ? l - a.off - kK : a.off;
                out.push_back({a.id, genome, contig, static_cast<std::uint32_t>(pos + p), rev ? !a.fwd : a.fwd});
            }
            if (i + 1 < w.steps.size()) {
                junction_key(g, w, i, j);
                auto it = junction_memo.find(j.key);
                if (it == junction_memo.end()) {
                    const std::string s = junction_string(g, j.key);
                    std::vector<anchors::Hit> h;
                    anchors::scan_anchors(lookup, s.data(), s.size(), 0, 0, 0, h);
                    std::vector<NodeAnchor> v2;
                    for (const auto& x : h) v2.push_back({x.pos, x.anchor, x.forward});
                    it = junction_memo.emplace(j.key, std::move(v2)).first;
                }
                for (const auto& a : it->second)
                    out.push_back({a.id, genome, contig, static_cast<std::uint32_t>(pos + (l - j.s) + a.off), a.fwd});
            }
            pos += l;
        }
        ++contig;
    });
    return out;
}

std::vector<std::string> graph_contigs(const std::string& gfa, std::vector<std::uint64_t>* lengths) {
    const Graph g = load_nodes(gfa);
    std::vector<std::string> out;
    for_each_walk(gfa, g, [&](const Walk& w) {
        out.push_back(w.contig);
        if (lengths != nullptr) {
            std::uint64_t len = w.offset;
            for (const auto& st : w.steps) len += g.len(st.first);
            lengths->push_back(len);
        }
    });
    return out;
}

struct Haplotypes::Impl {
    Graph g;
    std::vector<Walk> walks;
    std::vector<std::vector<std::uint64_t>> starts;   // per walk, the contig coordinate of each step
};

Haplotypes::Haplotypes(const std::string& gfa) : impl_(std::make_unique<Impl>()) {
    impl_->g = load_nodes(gfa);
    for_each_walk(gfa, impl_->g, [&](const Walk& w) {
        impl_->walks.push_back(w);
        std::vector<std::uint64_t> st;
        std::uint64_t pos = w.offset;
        for (const auto& s : w.steps) { st.push_back(pos); pos += impl_->g.len(s.first); }
        st.push_back(pos);
        impl_->starts.push_back(std::move(st));
    });
}

Haplotypes::~Haplotypes() = default;

std::string Haplotypes::fetch(const std::string& contig, std::uint64_t start, std::uint64_t end) const {
    for (std::size_t wi = 0; wi < impl_->walks.size(); ++wi) {
        const Walk& w = impl_->walks[wi];
        const auto& st = impl_->starts[wi];
        if (w.contig != contig || start >= st.back() || end <= st.front()) continue;
        std::string out;
        std::size_t i = static_cast<std::size_t>(std::upper_bound(st.begin(), st.end(), std::max(start, st.front())) - st.begin()) - 1;
        for (; i + 1 < st.size() && st[i] < end; ++i) {
            const auto [v, rev] = w.steps[i];
            const std::uint64_t from = start > st[i] ? start - st[i] : 0;
            const std::uint64_t to = std::min<std::uint64_t>(impl_->g.len(v), end - st[i]);
            impl_->g.append(out, v, rev, from, to);
        }
        return out;
    }
    throw std::runtime_error("background: no haplotype " + contig + " covering " + std::to_string(start) + "-" + std::to_string(end));
}

void spell_graph(const std::string& gfa, const std::string& fasta_out) {
    const Graph g = load_nodes(gfa);
    std::ofstream out(fasta_out);
    if (!out) throw std::runtime_error("background: cannot write " + fasta_out);
    std::string s;
    for_each_walk(gfa, g, [&](const Walk& w) {
        s.clear();
        for (const auto& [v, rev] : w.steps) g.append(s, v, rev, 0, g.len(v));
        out << '>' << w.contig << '\n';
        for (std::size_t i = 0; i < s.size(); i += 80) out << s.substr(i, 80) << '\n';
    });
    if (!out) throw std::runtime_error("background: cannot write " + fasta_out);
}

}  // namespace panvar::bg::graph
