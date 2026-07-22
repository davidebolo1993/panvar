#include "panvar/rebuild.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_set>
#include <vector>

#include "panvar/gfa.hpp"
#include "panvar/parallel.hpp"

#include "gfa.h"
#include "minigraph.h"
#include "mgpriv.h" // mg_gchain_free: declared here rather than in the public header

namespace panvar {
namespace {

char comp(char c) {
    switch (c) {
        case 'A': return 'T'; case 'T': return 'A'; case 'C': return 'G'; case 'G': return 'C';
        case 'a': return 't'; case 't': return 'a'; case 'c': return 'g'; case 'g': return 'c';
        default: return c;
    }
}

std::string reverse_complement(const std::string& s) {
    std::string r(s.rbegin(), s.rend());
    for (char& c : r) c = comp(c);
    return r;
}

std::string spell(const Graph& g, const PathRecord& path) {
    std::string out;
    for (const PathStep& step : path.steps) {
        const auto it = g.nodes.find(step.node_id);
        if (it == g.nodes.end()) continue;
        out += step.reverse ? reverse_complement(it->second.sequence) : it->second.sequence;
    }
    return out;
}

void degree_stats(const Graph& g, std::size_t hub_degree, std::size_t& hubs, std::size_t& maxdeg) {
    hubs = 0;
    maxdeg = 0;
    for (const auto& kv : g.nodes) {
        const std::vector<std::string> nb = g.neighbors_of(kv.first);
        const std::unordered_set<std::string> distinct(nb.begin(), nb.end());
        if (distinct.size() > maxdeg) maxdeg = distinct.size();
        if (distinct.size() >= hub_degree) ++hubs;
    }
}

struct KmerStats {
    std::size_t distinct = 0; // how much DIFFERENT sequence the haplotype carries
    std::size_t total = 0;    // k-mer positions; total - distinct is repeat-copy redundancy
};

// Rolling 2-bit encoding: substr-per-position costs tens of millions of allocations at cohort scale.
// k <= 31 fits a uint64; above that fall back to substr.
KmerStats kmer_stats(const std::string& s, std::size_t k) {
    KmerStats st;
    if (k == 0 || s.size() < k) return st;
    st.total = s.size() - k + 1;
    if (k > 31) {
        std::unordered_set<std::string> set;
        set.reserve(s.size());
        for (std::size_t i = 0; i + k <= s.size(); ++i) set.insert(s.substr(i, k));
        st.distinct = set.size();
        return st;
    }
    std::unordered_set<std::uint64_t> set;
    set.reserve(s.size());
    const std::uint64_t mask = (k < 32) ? ((1ULL << (2 * k)) - 1) : ~0ULL;
    std::uint64_t code = 0;
    std::size_t valid = 0;
    for (std::size_t i = 0; i < s.size(); ++i) {
        int b;
        switch (s[i]) {
            case 'A': case 'a': b = 0; break;
            case 'C': case 'c': b = 1; break;
            case 'G': case 'g': b = 2; break;
            case 'T': case 't': b = 3; break;
            default: valid = 0; code = 0; continue; // ambiguity resets the roll
        }
        code = ((code << 2) | static_cast<std::uint64_t>(b)) & mask;
        if (++valid >= k) set.insert(code);
    }
    st.distinct = set.size();
    return st;
}

void copy_file(const std::string& from, const std::string& to) {
    std::ifstream in(from, std::ios::binary);
    if (!in) throw std::runtime_error("rebuild: cannot open input: " + from);
    std::ofstream out(to, std::ios::binary);
    if (!out) throw std::runtime_error("rebuild: cannot open output: " + to);
    out << in.rdbuf();
}

// One-record FASTAs, deleted on scope exit: mg_ggen consumes files, so haplotypes round-trip through
// disk. Named by rank so the progressive order stays visible when debugging a run.
class FastaScratch {
public:
    explicit FastaScratch(const std::string& dir) : dir_(dir) {
        std::string cmd = "mkdir -p '" + dir_ + "'";
        if (std::system(cmd.c_str()) != 0) throw std::runtime_error("rebuild: cannot create " + dir_);
    }
    ~FastaScratch() {
        const std::string cmd = "rm -rf '" + dir_ + "'";
        if (std::system(cmd.c_str()) != 0) {
            std::cerr << "[rebuild] warning: could not remove scratch dir " << dir_ << '\n';
        }
    }
    FastaScratch(const FastaScratch&) = delete;
    FastaScratch& operator=(const FastaScratch&) = delete;

    // `rank` is the position in richness order.
    std::string write(std::size_t rank, const std::string& name, const std::string& seq) {
        char buf[32];
        std::snprintf(buf, sizeof(buf), "/h%04zu.fa", rank);
        const std::string path = dir_ + buf;
        std::ofstream o(path);
        if (!o) throw std::runtime_error("rebuild: cannot write " + path);
        o << '>' << name << '\n';
        for (std::size_t i = 0; i < seq.size(); i += 60) o << seq.substr(i, 60) << '\n';
        return path;
    }

private:
    std::string dir_;
};

// Owns the gfa_t so every early return frees it.
class GfaHandle {
public:
    explicit GfaHandle(gfa_t* g) : g_(g) {}
    ~GfaHandle() { if (g_) gfa_destroy(g_); }
    GfaHandle(const GfaHandle&) = delete;
    GfaHandle& operator=(const GfaHandle&) = delete;
    gfa_t* get() const { return g_; }

private:
    gfa_t* g_;
};

} // namespace

RebuildSummary rebuild_graph(const RebuildOptions& options) {
    RebuildSummary sum;
    const auto t0 = std::chrono::steady_clock::now();
    auto secs = [&] {
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    };
    auto step = [&](const std::string& msg) {
        if (!options.quiet) std::cerr << "[rebuild " << static_cast<long>(secs()) << "s] " << msg << std::endl;
    };

    step("reading " + options.gfa_path);
    ParseGfaOptions po;
    po.include_paths = true;
    po.include_sequences = true;
    const Graph g = parse_gfa(options.gfa_path, po);

    sum.raw_nodes = g.nodes.size();
    sum.haplotypes = g.paths.size();
    degree_stats(g, options.hub_degree, sum.raw_hubs, sum.raw_maxdeg);
    step("ordering " + std::to_string(sum.haplotypes) + " haplotypes by k-mer richness");

    // ---- Criteria B: order haplotypes by k-mer richness ----
    // Lexicographic, diversity first: (distinct k-mers, total k-mers) descending. Abundance only
    // discriminates between haplotypes carrying the same amount of distinct sequence -- which is the
    // common case, not an edge case. Per-haplotype independent, so it runs across all workers.
    std::vector<std::string> hseq(g.paths.size());
    std::vector<std::pair<std::size_t, std::size_t>> score(g.paths.size());
    std::vector<std::size_t> idx(g.paths.size());
    run_parallel(g.paths.size(), options.threads, [&](std::size_t i) {
        hseq[i] = spell(g, g.paths[i]);
        const KmerStats st = kmer_stats(hseq[i], options.kmer);
        score[i] = {st.distinct, st.total};
    });
    std::size_t max_len = 1;
    for (std::size_t i = 0; i < g.paths.size(); ++i) {
        max_len = std::max(max_len, hseq[i].size());
        idx[i] = i;
    }
    std::sort(idx.begin(), idx.end(), [&](std::size_t a, std::size_t b) { return score[a] > score[b]; });
    if (!idx.empty()) sum.seed = g.paths[idx.front()].name;
    sum.raw_density = static_cast<double>(sum.raw_nodes) / (static_cast<double>(max_len) / 1000.0);

    // ---- Criteria A: the pathology gate ----
    sum.pathological = sum.raw_hubs >= options.min_hubs;
    if (!options.quiet) {
        std::cerr << "[rebuild " << static_cast<long>(secs()) << "s] gate: " << sum.raw_nodes
                  << " nodes, " << sum.haplotypes
                  << " haps; #deg>=" << options.hub_degree << "=" << sum.raw_hubs
                  << " maxdeg=" << sum.raw_maxdeg << " density=" << static_cast<long>(sum.raw_density)
                  << "/kb -> " << (sum.pathological ? "PATHOLOGICAL" : "healthy")
                  << "; seed=" << sum.seed << '\n';
    }
    if (!sum.pathological && !options.force) {
        if (!options.quiet) std::cerr << "[rebuild] healthy -> pass through unchanged\n";
        if (!options.out_path.empty()) copy_file(options.gfa_path, options.out_path);
        sum.out_nodes = sum.raw_nodes;
        return sum;
    }
    if (idx.empty()) throw std::runtime_error("rebuild: no haplotype paths in " + options.gfa_path);
    sum.ran = true;

    unsigned nthreads = options.threads > 0 ? static_cast<unsigned>(options.threads)
                                            : std::max(1u, std::thread::hardware_concurrency());
    const int n_threads = static_cast<int>(nthreads);

    // ---- progressive graph generation, in richness order ----
    // minigraph reads the first file as the seed graph and augments it with each subsequent one in turn,
    // so passing our richness order drives the progressive build directly.
    // Scratch FASTAs live beside --out by default, or under --tmp-dir in a dedicated subfolder (named
    // from the output basename) so cleanup removes only what we created, never the whole --tmp-dir.
    std::string scratch_dir = options.out_path + ".rebuild.tmp";
    if (!options.tmp_dir.empty()) {
        std::string base = options.out_path;
        const auto slash = base.find_last_of('/');
        if (slash != std::string::npos) base = base.substr(slash + 1);
        scratch_dir = options.tmp_dir + "/" + base + ".rebuild.tmp";
    }
    step("writing haplotype FASTAs to " + scratch_dir);
    FastaScratch scratch(scratch_dir);
    std::vector<std::string> fasta_paths;
    fasta_paths.reserve(idx.size());
    for (std::size_t rank = 0; rank < idx.size(); ++rank) {
        const std::size_t h = idx[rank];
        fasta_paths.push_back(scratch.write(rank, g.paths[h].name, hseq[h]));
    }

    mg_idxopt_t ipt;
    mg_mapopt_t opt;
    mg_ggopt_t gpt;
    mg_opt_set(nullptr, &ipt, &opt, &gpt);
    if (mg_opt_set("ggs", &ipt, &opt, &gpt) < 0) throw std::runtime_error("rebuild: mg_opt_set(ggs) failed");
    gpt.min_var_len = static_cast<int>(options.min_var);
    // Keep minigraph itself quiet (its per-sample logs are one block per haplotype, far too noisy at
    // cohort scale) and print our own throttled counter instead. mg_ggen augments one file per call
    // internally, so driving that loop here costs nothing extra and lets us report progress.
    mg_verbose = 1;

    GfaHandle handle(gfa_read(fasta_paths.front().c_str()));
    if (handle.get() == nullptr) {
        throw std::runtime_error("rebuild: minigraph could not seed from " + fasta_paths.front());
    }
    const std::size_t total = fasta_paths.size();
    step("generating graph: seed + " + std::to_string(total - 1) +
         " haplotypes in richness order (minigraph)");
    for (std::size_t i = 1; i < fasta_paths.size(); ++i) {
        const char* f = fasta_paths[i].c_str();
        if (mg_ggen(handle.get(), 1, &f, &ipt, &opt, &gpt, n_threads) != 0) {
            throw std::runtime_error("rebuild: mg_ggen failed on " + fasta_paths[i]);
        }
        if (i % 50 == 0 || i + 1 == total) {
            step("  added " + std::to_string(i + 1) + "/" + std::to_string(total) + " haplotypes; " +
                 std::to_string(handle.get()->n_seg) + " segments so far");
        }
    }
    gfa_t* out = handle.get();
    step("generated " + std::to_string(out->n_seg) + " segments, " + std::to_string(out->n_arc) +
         " arcs (min-var " + std::to_string(options.min_var) + ")");

    // ---- path recovery: map every haplotype back to the graph, in memory (no GAF round-trip) ----
    // Each mg_llchain_t.v is an ORIENTED vertex (seg<<1|strand), which is what lets inversion
    // traversals survive into the P lines.
    mg_idxopt_t ipt_map;
    mg_mapopt_t opt_map;
    mg_ggopt_t gpt_map;
    mg_opt_set(nullptr, &ipt_map, &opt_map, &gpt_map);
    if (mg_opt_set("asm", &ipt_map, &opt_map, &gpt_map) < 0) {
        throw std::runtime_error("rebuild: mg_opt_set(asm) failed");
    }
    mg_idx_t* gi = mg_index(out, &ipt_map, n_threads, &opt_map);
    if (gi == nullptr) throw std::runtime_error("rebuild: mg_index failed");
    step("recovering " + std::to_string(g.paths.size()) + " haplotype walks");

    // Parallel like minigraph's own kt_for: index and options shared read-only, one mg_tbuf_t per
    // worker (what makes mg_map re-entrant). Results go to distinct indices, so output is
    // thread-order independent.
    std::vector<std::vector<std::uint32_t>> walks(g.paths.size());
    std::vector<double> cover(g.paths.size(), 0.0);
    {
        std::atomic<std::size_t> next{0};
        std::vector<std::thread> pool;
        pool.reserve(nthreads);
        for (unsigned t = 0; t < nthreads; ++t) {
            pool.emplace_back([&] {
                mg_tbuf_t* tbuf = mg_tbuf_init();
                for (std::size_t h = next.fetch_add(1); h < g.paths.size(); h = next.fetch_add(1)) {
                    const std::string& q = hseq[h];
                    if (q.empty()) continue;
                    mg_gchains_t* gc = mg_map(gi, static_cast<int>(q.size()), q.c_str(), tbuf, &opt_map,
                                              g.paths[h].name.c_str());
                    if (gc != nullptr && gc->n_gc > 0) {
                        // Longest chain wins: concatenating chains would stitch an incoherent walk.
                        int best = 0;
                        for (int i = 1; i < gc->n_gc; ++i) {
                            if (gc->gc[i].qe - gc->gc[i].qs > gc->gc[best].qe - gc->gc[best].qs) best = i;
                        }
                        const mg_gchain_t& c = gc->gc[best];
                        walks[h].reserve(static_cast<std::size_t>(c.cnt));
                        for (int32_t j = 0; j < c.cnt; ++j) walks[h].push_back(gc->lc[c.off + j].v);
                        cover[h] = static_cast<double>(c.qe - c.qs) / static_cast<double>(q.size());
                    }
                    if (gc != nullptr) mg_gchain_free(gc);
                }
                mg_tbuf_destroy(tbuf);
            });
        }
        for (std::thread& th : pool) th.join();
    }
    mg_idx_destroy(gi);
    double cover_sum = 0.0;
    for (std::size_t h = 0; h < walks.size(); ++h) {
        if (!walks[h].empty()) ++sum.paths_recovered;
        cover_sum += cover[h];
    }
    sum.mean_query_cover =
        g.paths.empty() ? 0.0 : cover_sum / static_cast<double>(g.paths.size());

    // ---- emit plain GFA: S + oriented L + one P per haplotype ----
    // Orientation MUST be carried through: flattening links to +/+ silently destroys every inversion
    // bubble downstream.
    // Segments are renamed to consecutive integers: minigraph names them "s1", "s2", ... and tooling
    // that stores node ids as integers (odgi) rejects non-numeric names outright. The original name is
    // kept as an SN tag.
    if (!options.out_path.empty()) {
        std::ofstream o(options.out_path);
        if (!o) throw std::runtime_error("rebuild: cannot write " + options.out_path);
        auto seg_id = [](std::uint32_t seg) { return seg + 1; }; // 0-based index -> 1-based GFA id
        o << "H\tVN:Z:1.0\n";
        for (std::uint32_t i = 0; i < out->n_seg; ++i) {
            const gfa_seg_t& s = out->seg[i];
            o << "S\t" << seg_id(i) << '\t' << (s.seq ? s.seq : "*");
            if (s.name != nullptr) o << "\tSN:Z:" << s.name;
            o << '\n';
        }
        for (std::uint64_t k = 0; k < out->n_arc; ++k) {
            const gfa_arc_t& a = out->arc[k];
            if (a.del || a.comp) continue; // skip deleted and complement (dual) arcs, as gfa_print does
            o << "L\t" << seg_id(static_cast<std::uint32_t>(a.v_lv >> 33)) << '\t'
              << "+-"[(a.v_lv >> 32) & 1] << '\t' << seg_id(a.w >> 1) << '\t' << "+-"[a.w & 1]
              << "\t0M\n";
        }
        for (std::size_t h = 0; h < g.paths.size(); ++h) {
            if (walks[h].empty()) continue;
            o << "P\t" << g.paths[h].name << '\t';
            for (std::size_t j = 0; j < walks[h].size(); ++j) {
                if (j) o << ',';
                o << seg_id(walks[h][j] >> 1) << ((walks[h][j] & 1) ? '-' : '+');
            }
            o << "\t*\n";
        }
    }

    // ---- structural summary of what we emitted ----
    sum.out_nodes = out->n_seg;
    {
        std::vector<std::unordered_set<std::uint32_t>> nb(out->n_seg);
        std::size_t edges = 0;
        for (std::uint64_t k = 0; k < out->n_arc; ++k) {
            const gfa_arc_t& a = out->arc[k];
            if (a.del || a.comp) continue;
            ++edges;
            const std::uint32_t u = static_cast<std::uint32_t>(a.v_lv >> 33);
            const std::uint32_t v = a.w >> 1;
            nb[u].insert(v);
            nb[v].insert(u);
        }
        sum.out_edges = edges;
        for (const auto& s : nb) {
            sum.out_maxdeg = std::max(sum.out_maxdeg, s.size());
            if (s.size() >= options.hub_degree) ++sum.out_hubs;
        }
    }
    if (sum.paths_recovered < sum.haplotypes) {
        std::cerr << "[rebuild] warning: " << (sum.haplotypes - sum.paths_recovered) << " of "
                  << sum.haplotypes << " haplotypes could not be mapped back and have no P line\n";
    }
    if (!options.quiet) {
        std::cerr << "[rebuild " << static_cast<long>(secs()) << "s] emitted " << sum.out_nodes
                  << " nodes, " << sum.out_edges << " edges; maxdeg=" << sum.out_maxdeg << " #hub="
                  << sum.out_hubs << "; paths " << sum.paths_recovered << '/' << sum.haplotypes
                  << ", mean query cover " << sum.mean_query_cover << std::endl;
    }
    return sum;
}

} // namespace panvar
