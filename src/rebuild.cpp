#include "panvar/rebuild.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <atomic>
#include <cstdint>
#include <cstdio>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <random>
#include <thread>
#include <unordered_set>
#include <vector>

#include "panvar/align.hpp"
#include "panvar/gfa.hpp"
#include "panvar/graph_utils.hpp"
#include "panvar/parallel.hpp"

#include "gfa.h"
#include "minigraph.h"
#include "mgpriv.h" // mg_gchain_free: declared here rather than in the public header

namespace panvar {
namespace {

// Wall-clock HH:MM:SS, matching cli::RunLog's line prefix so rebuild's progress lines read the same as
// its RunLog summary lines and as every other module.
std::string hms() {
    const std::time_t t = std::time(nullptr);
    std::tm tm{};
#if defined(_WIN32)
    localtime_s(&tm, &t);
#else
    localtime_r(&t, &tm);
#endif
    char buf[16];
    std::strftime(buf, sizeof(buf), "%H:%M:%S", &tm);
    return std::string(buf);
}

// reverse_complement lives in graph_utils; a second copy here was one more thing that could drift.

std::string spell(const Graph& g, const PathRecord& path) {
    std::string out;
    for (const PathStep& step : path.steps) {
        const auto it = g.nodes.find(step.node_id);
        // validate_graph_paths() has already refused a graph where this could miss; the guard means a
        // future caller cannot silently spell a short sequence.
        if (it == g.nodes.end())
            throw std::runtime_error("rebuild: path " + path.name + " references missing node " + step.node_id);
        out += step.reverse ? reverse_complement(it->second.sequence) : it->second.sequence;
    }
    return out;
}

// Degree per HANDLE, not per node. A bidirected graph node has two ends, and pooling them conflates two
// different shapes: a node with 25 neighbours on each side is a clean two-sided branch, while one with
// 50 on a single side is the tangle the gate exists to find -- yet pooled they both read 50. The node's
// degree is therefore the larger of its two handle degrees. Self-loops are counted separately: they are
// how a tandem array appears after folding, and they are not pathology.
void degree_stats(const Graph& g, std::size_t hub_degree, std::size_t& hubs, std::size_t& maxdeg,
                  std::size_t* selfloops) {
    hubs = 0;
    maxdeg = 0;
    std::size_t loops = 0;
    for (const auto& kv : g.nodes) {
        std::unordered_set<std::string> left, right;
        bool loop = false;
        for (const Neighbor& n : kv.second.start) {
            if (n.node_id == kv.first) { loop = true; continue; }
            left.insert(n.node_id);
        }
        for (const Neighbor& n : kv.second.end) {
            if (n.node_id == kv.first) { loop = true; continue; }
            right.insert(n.node_id);
        }
        if (loop) ++loops;
        const std::size_t deg = std::max(left.size(), right.size());
        if (deg > maxdeg) maxdeg = deg;
        if (deg >= hub_degree) ++hubs;
    }
    if (selfloops != nullptr) *selfloops = loops;
}

struct KmerStats {
    std::size_t distinct = 0; // how much DIFFERENT sequence the haplotype carries
    std::size_t total = 0;    // k-mer positions; total - distinct is repeat-copy redundancy
};

// Rolling 2-bit encoding: substr-per-position costs tens of millions of allocations at cohort scale.
// k <= 31 fits a uint64; above that fall back to substr.
//
// Three properties this has to have, and previously did not:
//   CANONICAL   -- a k-mer and its reverse complement count as one. Without it a haplotype's richness
//                  depends on which strand the GFA happens to store it on, so the seed -- and hence the
//                  whole rebuild -- changes when an input is flipped.
//   CONSISTENT  -- `total` counts only the windows `distinct` counts. It was s.size()-k+1, including
//                  windows containing N, while `distinct` skipped them, so the redundancy figure
//                  (total - distinct) was wrong by however much ambiguity a haplotype carried.
//   UNIFORM     -- k > 31 behaves the same way. It kept every window verbatim, ambiguity and strand
//                  included, so crossing k=31 silently changed what the ranking meant.
KmerStats kmer_stats(const std::string& s, std::size_t k) {
    KmerStats st;
    if (k == 0 || s.size() < k) return st;
    if (k > 31) {
        std::unordered_set<std::string> set;
        set.reserve(s.size());
        for (std::size_t i = 0; i + k <= s.size(); ++i) {
            std::string up = s.substr(i, k);
            bool ok = true;
            for (char& c : up) {
                c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
                if (c != 'A' && c != 'C' && c != 'G' && c != 'T') { ok = false; break; }
            }
            if (!ok) continue;                        // ambiguity: not a countable window
            const std::string rc = reverse_complement(up);
            set.insert(up < rc ? up : rc);            // canonical
            ++st.total;
        }
        st.distinct = set.size();
        return st;
    }
    std::unordered_set<std::uint64_t> set;
    set.reserve(s.size());
    const std::uint64_t mask = (k < 32) ? ((1ULL << (2 * k)) - 1) : ~0ULL;
    const unsigned shift = static_cast<unsigned>(2 * (k - 1));
    std::uint64_t fwd = 0, rev = 0;
    std::size_t valid = 0;
    for (std::size_t i = 0; i < s.size(); ++i) {
        int b;
        switch (s[i]) {
            case 'A': case 'a': b = 0; break;
            case 'C': case 'c': b = 1; break;
            case 'G': case 'g': b = 2; break;
            case 'T': case 't': b = 3; break;
            default: valid = 0; fwd = 0; rev = 0; continue;  // ambiguity resets the roll
        }
        fwd = ((fwd << 2) | static_cast<std::uint64_t>(b)) & mask;
        rev = (rev >> 2) | (static_cast<std::uint64_t>(3 - b) << shift);
        if (++valid >= k) { set.insert(std::min(fwd, rev)); ++st.total; }
    }
    st.distinct = set.size();
    return st;
}

// A unique sibling of `final` -- outputs are staged here and renamed into place, so an interrupted or
// failed run never leaves a half-written graph where a complete one is expected.
std::filesystem::path staging_path(const std::string& final_path) {
    static std::atomic<unsigned> seq{0};
    std::random_device rd;
    const std::string tag = std::to_string(rd()) + "." + std::to_string(seq++);
    return std::filesystem::path(final_path + ".rebuild-tmp." + tag);
}

void commit_staged(const std::filesystem::path& staged, const std::string& final_path) {
    std::error_code ec;
    std::filesystem::rename(staged, final_path, ec);
    if (ec) {
        // Across filesystems rename fails; fall back to copy-then-remove, still not in place.
        std::filesystem::copy_file(staged, final_path,
                                   std::filesystem::copy_options::overwrite_existing, ec);
        std::filesystem::remove(staged);
        if (ec) throw std::runtime_error("rebuild: cannot move output into place: " + final_path);
    }
}

void copy_file(const std::string& from, const std::string& to) {
    // Staged then renamed: writing straight to `to` truncates it on open, which destroyed the input
    // outright when the same path was given for -i and -o. The guard in run_rebuild rejects that case,
    // but a pass-through should not depend on the caller having been checked.
    const std::filesystem::path staged = staging_path(to);
    {
        std::ifstream in(from, std::ios::binary);
        if (!in) throw std::runtime_error("rebuild: cannot open input: " + from);
        std::ofstream out(staged, std::ios::binary);
        if (!out) throw std::runtime_error("rebuild: cannot open output: " + staged.string());
        out << in.rdbuf();
        if (!out) throw std::runtime_error("rebuild: failed writing " + staged.string());
    }
    commit_staged(staged, to);
}

// One-record FASTAs, deleted on scope exit: mg_ggen consumes files, so haplotypes round-trip through
// disk. Named by rank so the progressive order stays visible when debugging a run.
class FastaScratch {
public:
    explicit FastaScratch(const std::string& dir) : dir_(dir) {
        // create_directory (not create_directories) so an EXISTING directory is an error rather than
        // something this object will delete on scope exit. The name carries a random tag, so two
        // concurrent runs cannot collide and neither can adopt a directory it did not make.
        std::error_code ec;
        std::filesystem::create_directories(std::filesystem::path(dir_).parent_path(), ec);
        if (!std::filesystem::create_directory(dir_, ec) || ec)
            throw std::runtime_error("rebuild: cannot create scratch directory " + dir_ +
                                     (ec ? " (" + ec.message() + ")" : " (already exists)"));
    }
    ~FastaScratch() {
        std::error_code ec;
        std::filesystem::remove_all(dir_, ec);
        if (ec) {
            // No shell fallback: `rm -rf` on a path this class did not necessarily create is a worse
            // outcome than a leftover directory, and the path is interpolated into a command line.
            std::cerr << "[rebuild " << hms() << "] warning: could not remove scratch dir " << dir_
                      << " (" << ec.message() << ")\n";
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

std::string fmt2(double v) {
    char b[32];
    std::snprintf(b, sizeof(b), "%.4f", v);
    return b;
}

// Enough digits to reproduce the decision. The four-decimal form printed a --min-recovered-identity of
// 0.999999 as "1.0000", so the audit could not be used to check why a path was rejected.
std::string fmt_exact(double v) {
    char b[40];
    std::snprintf(b, sizeof(b), "%.10g", v);
    return b;
}

RebuildSummary rebuild_graph(const RebuildOptions& options) {
    RebuildSummary sum;
    auto step = [&](const std::string& msg) {
        if (!options.quiet) std::cerr << "[rebuild " << hms() << "] " << msg << std::endl;
    };

    step("reading " + options.gfa_path);
    ParseGfaOptions po;
    po.include_paths = true;
    po.include_sequences = true;
    const Graph g = parse_gfa(options.gfa_path, po);
    validate_graph_paths(g, "rebuild", true, true);

    sum.raw_nodes = g.nodes.size();
    sum.haplotypes = g.paths.size();
    degree_stats(g, options.hub_degree, sum.raw_hubs, sum.raw_maxdeg, &sum.raw_selfloops);
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
    // Richness descending, then path NAME, then original index. std::sort is not stable, so equal
    // richness previously produced an implementation-defined order -- two runs on the same input could
    // seed from different haplotypes and build different graphs. Ties are the norm in a panel of
    // near-duplicate haplotypes, not an edge case.
    std::sort(idx.begin(), idx.end(), [&](std::size_t a, std::size_t b) {
        if (score[a] != score[b]) return score[a] > score[b];
        if (g.paths[a].name != g.paths[b].name) return g.paths[a].name < g.paths[b].name;
        return a < b;
    });
    if (!idx.empty()) sum.seed = g.paths[idx.front()].name;
    sum.raw_density = static_cast<double>(sum.raw_nodes) / (static_cast<double>(max_len) / 1000.0);

    // ---- Criteria A: the pathology gate ----
    sum.pathological = sum.raw_hubs >= options.min_hubs;
    if (!options.quiet) {
        std::cerr << "[rebuild " << hms() << "] gate: " << sum.raw_nodes
                  << " nodes, " << sum.haplotypes
                  << " haps; #deg>=" << options.hub_degree << "=" << sum.raw_hubs
                  << " maxdeg=" << sum.raw_maxdeg << " density=" << static_cast<long>(sum.raw_density)
                  << "/kb -> " << (sum.pathological ? "PATHOLOGICAL" : "healthy")
                  << "; seed=" << sum.seed << '\n';
    }
    if (!sum.pathological && !options.force) {
        step("healthy -> pass through unchanged");
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
    std::random_device scratch_rd;
    const std::string scratch_tag = ".rebuild.tmp." + std::to_string(scratch_rd());
    std::string scratch_dir = options.out_path + scratch_tag;
    if (!options.tmp_dir.empty()) {
        std::string base = options.out_path;
        const auto slash = base.find_last_of('/');
        if (slash != std::string::npos) base = base.substr(slash + 1);
        scratch_dir = options.tmp_dir + "/" + base + scratch_tag;
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
    // minigraph's length gates assume chromosome-scale input: ggsimple drops a chain wholesale when
    // its block length is under min_map_len (100 kb by default), before min_var_len is ever consulted.
    // A locus graph of a few tens of kb therefore augments nothing and collapses to the bare seed, so
    // scale both gates to the locus. Loci already above the defaults keep minigraph's own values.
    // The gate is derived from the whole length DISTRIBUTION, not from the seed. The seed is the richest
    // haplotype and can be far longer than the rest; a threshold scaled from it can exceed what a
    // shorter haplotype could ever produce, so that haplotype clears no chain and is dropped without
    // anything saying so. The median sets the scale and the SHORTEST haplotype caps it, so every
    // haplotype can in principle contribute a chain covering half of itself.
    std::vector<std::size_t> hlen;
    hlen.reserve(hseq.size());
    for (const std::string& h : hseq) if (!h.empty()) hlen.push_back(h.size());
    if (hlen.empty()) throw std::runtime_error("rebuild: every haplotype spelled an empty sequence");
    std::sort(hlen.begin(), hlen.end());
    const std::size_t min_len = hlen.front();
    const std::size_t med_len = hlen[hlen.size() / 2];
    const std::size_t seed_len = hseq[idx.front()].size();
    // Never above half the shortest haplotype, whatever else says.
    const int cap = static_cast<int>(std::max<std::size_t>(500, min_len / 2));
    if (options.min_align_len > 0) {
        gpt.min_map_len = static_cast<int>(options.min_align_len);
        if (gpt.min_map_len > cap) {
            step("warning: --min-align-len " + std::to_string(options.min_align_len) +
                 " exceeds half the shortest haplotype (" + std::to_string(min_len) +
                 " bp); haplotypes shorter than twice it cannot clear the gate");
        }
        if (static_cast<std::size_t>(gpt.min_depth_len) > options.min_align_len) {
            gpt.min_depth_len = static_cast<int>(std::max<std::size_t>(500, options.min_align_len / 5));
        }
        step("min alignment length " + std::to_string(gpt.min_map_len) + " (requested), min depth length " +
             std::to_string(gpt.min_depth_len));
    } else {
        const int base = (med_len < 2 * static_cast<std::size_t>(gpt.min_map_len))
                             ? static_cast<int>(std::max<std::size_t>(1000, med_len / 2))
                             : gpt.min_map_len;
        const int chosen = std::min(base, cap);
        if (chosen != gpt.min_map_len) {
            gpt.min_map_len = chosen;
            if (med_len < 2 * static_cast<std::size_t>(gpt.min_depth_len))
                gpt.min_depth_len = static_cast<int>(std::max<std::size_t>(500, med_len / 10));
            gpt.min_depth_len = std::min(gpt.min_depth_len, gpt.min_map_len);
            step("haplotype lengths " + std::to_string(min_len) + "-" + std::to_string(hlen.back()) +
                 " bp (median " + std::to_string(med_len) + ", seed " + std::to_string(seed_len) +
                 "): min alignment length " + std::to_string(gpt.min_map_len) +
                 (chosen == cap && cap < base ? " (capped by the shortest haplotype)" : "") +
                 ", min depth length " + std::to_string(gpt.min_depth_len));
        }
    }
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

    // A graph with no arcs carries no variation at all: every haplotype would spell the seed, which is
    // strictly worse than the input. Treat that as "not rebuildable" and keep the original graph.
    if (out->n_arc == 0 || out->n_seg <= 1) {
        step("degenerate rebuild (no variation recovered) -> pass through unchanged");
        if (!options.out_path.empty()) copy_file(options.gfa_path, options.out_path);
        sum.out_nodes = sum.raw_nodes;
        sum.ran = false;
        return sum;
    }

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
    std::vector<double> cover(g.paths.size(), 0.0);        // envelope (qe-qs)/len
    std::vector<double> matched(g.paths.size(), 0.0);      // matching bases / len
    std::vector<double> chain_id(g.paths.size(), 0.0);     // mlen / blen
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
                        // One chain wins -- concatenating chains would stitch an incoherent walk -- but
                        // "longest query SPAN" ranks by the outer envelope, so a chain reaching further
                        // through a large internal gap beats one that genuinely aligns more. Rank by
                        // matching bases, then identity within the aligned block, then mapping quality,
                        // with a deterministic tie-break so the choice is reproducible.
                        auto better = [](const mg_gchain_t& x, const mg_gchain_t& y) {
                            if (x.mlen != y.mlen) return x.mlen > y.mlen;
                            const double ix = x.blen > 0 ? static_cast<double>(x.mlen) / x.blen : 0.0;
                            const double iy = y.blen > 0 ? static_cast<double>(y.mlen) / y.blen : 0.0;
                            if (ix != iy) return ix > iy;
                            if (x.mapq != y.mapq) return x.mapq > y.mapq;
                            if (x.score != y.score) return x.score > y.score;
                            return x.qs < y.qs;
                        };
                        int best = 0;
                        for (int i = 1; i < gc->n_gc; ++i) if (better(gc->gc[i], gc->gc[best])) best = i;
                        const mg_gchain_t& c = gc->gc[best];
                        walks[h].reserve(static_cast<std::size_t>(c.cnt));
                        for (int32_t j = 0; j < c.cnt; ++j) walks[h].push_back(gc->lc[c.off + j].v);
                        // Three different questions, previously answered by one number:
                        //   cover   -- where the alignment begins and ends (the outer envelope)
                        //   matched -- how much of the query is genuinely aligned (envelope minus gaps)
                        //   chain_id-- how well it matches where it does align
                        // A chain spanning the whole query through a large internal gap scores high on
                        // the first and low on the second, which is exactly the case the envelope hid.
                        cover[h] = static_cast<double>(c.qe - c.qs) / static_cast<double>(q.size());
                        matched[h] = static_cast<double>(c.mlen) / static_cast<double>(q.size());
                        chain_id[h] = c.blen > 0 ? static_cast<double>(c.mlen) / static_cast<double>(c.blen)
                                                 : 0.0;
                    }
                    if (gc != nullptr) mg_gchain_free(gc);
                }
                mg_tbuf_destroy(tbuf);
            });
        }
        for (std::thread& th : pool) th.join();
    }
    mg_idx_destroy(gi);
    double cover_sum = 0.0, matched_sum = 0.0, chainid_sum = 0.0;
    double matched_min = 1.0;
    std::size_t n_recovered = 0;
    for (std::size_t h = 0; h < walks.size(); ++h) {
        if (!walks[h].empty()) {
            ++sum.paths_recovered;
            ++n_recovered;
            chainid_sum += chain_id[h];
        }
        cover_sum += cover[h];
        matched_sum += matched[h];
        matched_min = std::min(matched_min, matched[h]);
    }
    sum.mean_query_cover =
        g.paths.empty() ? 0.0 : cover_sum / static_cast<double>(g.paths.size());
    sum.mean_matched_cover =
        g.paths.empty() ? 0.0 : matched_sum / static_cast<double>(g.paths.size());
    sum.min_matched_cover = g.paths.empty() ? 0.0 : matched_min;
    sum.mean_chain_identity = n_recovered ? chainid_sum / static_cast<double>(n_recovered) : 0.0;
    step("recovered " + std::to_string(sum.paths_recovered) + "/" + std::to_string(g.paths.size()) +
         " walks; envelope cover " + fmt2(sum.mean_query_cover) + ", matched cover " +
         fmt2(sum.mean_matched_cover) + " (min " + fmt2(sum.min_matched_cover) + "), chain identity " +
         fmt2(sum.mean_chain_identity));

    // ---- emit plain GFA: S + oriented L + one P per haplotype ----
    // Orientation MUST be carried through: flattening links to +/+ silently destroys every inversion
    // bubble downstream.
    // Segments are renamed to consecutive integers: minigraph names them "s1", "s2", ... and tooling
    // that stores node ids as integers (odgi) rejects non-numeric names outright. The original name is
    // kept as an SN tag.
    std::filesystem::path staged_out;
    if (!options.out_path.empty()) {
        const std::filesystem::path staged = staging_path(options.out_path);
        std::ofstream o(staged);
        if (!o) throw std::runtime_error("rebuild: cannot write " + staged.string());
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
        o.flush();
        if (!o) throw std::runtime_error("rebuild: failed writing " + staged.string());
        o.close();
        staged_out = staged;   // committed only if the acceptance contract below is satisfied
    }

    std::vector<double> walk_id;
    // ---- recovered-walk identity: the end-to-end check ----
    // Re-spell each recovered walk from the REBUILT graph and compare it with the haplotype the input
    // spelled. Chain coverage and identity describe the alignment; this describes what a caller
    // actually gets back, so it catches an error anywhere in mapping, chain choice or emission. It is
    // also the quantity an acceptance threshold has to be stated in.
    {
        std::vector<std::string> seg_seq(out->n_seg);
        for (std::uint32_t i = 0; i < out->n_seg; ++i)
            seg_seq[i] = out->seg[i].seq ? out->seg[i].seq : std::string();
        walk_id.assign(g.paths.size(), -1.0);
        std::vector<double>& wid = walk_id;
        run_parallel(g.paths.size(), options.threads, [&](std::size_t h) {
            if (walks[h].empty() || hseq[h].empty()) return;
            std::string spelled;
            for (const std::uint32_t v : walks[h]) {
                const std::uint32_t seg = v >> 1;
                if (seg >= seg_seq.size()) return;
                const std::string& sq = seg_seq[seg];
                spelled += (v & 1) ? reverse_complement(sq) : sq;
            }
            if (spelled.empty()) return;
            const NwAlign nw = nw_edit_distance(spelled, hseq[h]);
            wid[h] = nw.aln_len ? std::max(0.0, 1.0 - static_cast<double>(nw.edits) /
                                                     static_cast<double>(nw.aln_len))
                                : 0.0;
        });
        double wsum = 0.0, wmin = 1.0;
        std::size_t wn = 0, wpoor = 0;
        for (const double v : wid) {
            if (v < 0.0) continue;
            wsum += v; wmin = std::min(wmin, v); ++wn;
            if (v < 0.99) ++wpoor;   // the count an acceptance threshold would act on
        }
        sum.walks_below_99 = wpoor;
        sum.walk_identity_checked = wn;
        sum.mean_walk_identity = wn ? wsum / static_cast<double>(wn) : 0.0;
        sum.min_walk_identity = wn ? wmin : 0.0;
        step("recovered-walk identity over " + std::to_string(wn) + " haplotypes: mean " +
             fmt2(sum.mean_walk_identity) + ", worst " + fmt2(sum.min_walk_identity) + "; " +
             std::to_string(wpoor) + " below 0.99");
    }

    // ---- acceptance contract, and the audit that explains its verdict ----
    //
    // minigraph augments variation ABOVE --min-var, so sub-threshold differences are collapsed by
    // construction and a recovered walk is never byte-identical to its haplotype. The contract is
    // therefore structural plus a threshold, not losslessness: every path must come back, its steps
    // must be connected by edges that exist, and it must spell what it spelled before to within the
    // declared bounds. Anything less and the rebuilt graph is discarded and the ORIGINAL passed
    // through, because a graph that silently drops or truncates a haplotype is worse than no rebuild:
    // everything downstream would agree with it.
    {
        // Every consecutive pair of steps on an emitted walk needs a link that actually exists, or the
        // P line describes a traversal the graph does not permit.
        // Both orientations of every link. A GFA stores each edge once and minigraph marks the dual as
        // `comp`, which the emitter skips -- so a walk traversing an edge the other way round finds no
        // arc and would be reported as dangling when the graph is perfectly well formed. For v -> w the
        // dual is (w^1) -> (v^1).
        std::unordered_set<std::uint64_t> arcs;
        arcs.reserve(static_cast<std::size_t>(out->n_arc) * 4);
        for (std::uint64_t k = 0; k < out->n_arc; ++k) {
            const gfa_arc_t& a = out->arc[k];
            if (a.del) continue;
            const std::uint32_t v = static_cast<std::uint32_t>(a.v_lv >> 32);
            arcs.insert((static_cast<std::uint64_t>(v) << 32) | a.w);
            arcs.insert((static_cast<std::uint64_t>(a.w ^ 1) << 32) | (v ^ 1));
        }
        std::size_t missing_paths = 0, dangling = 0, failing = 0;
        std::string first_bad;
        for (std::size_t h = 0; h < g.paths.size(); ++h) {
            if (walks[h].empty()) {
                ++missing_paths;
                if (first_bad.empty()) first_bad = g.paths[h].name + " (no walk recovered)";
                continue;
            }
            for (std::size_t j = 1; j < walks[h].size(); ++j) {
                const std::uint64_t key = (static_cast<std::uint64_t>(walks[h][j - 1]) << 32) | walks[h][j];
                if (!arcs.count(key)) {
                    ++dangling;
                    if (first_bad.empty()) first_bad = g.paths[h].name + " (a step pair has no edge)";
                    break;
                }
            }
            // An identity that could not be COMPUTED is not evidence of a good recovery; it is the
            // absence of evidence. Treating id < 0 as passing let exactly the unverifiable cases through
            // the one check meant to catch them.
            const double id = h < walk_id.size() ? walk_id[h] : -1.0;
            if (matched[h] < options.min_matched_cover || id < 0.0 ||
                id < options.min_recovered_identity) {
                ++failing;
                if (first_bad.empty())
                    first_bad = g.paths[h].name +
                                (id < 0.0 ? " (recovered-walk identity could not be computed)"
                                          : " (identity " + fmt2(id) + ", matched cover " +
                                                fmt2(matched[h]) + ")");
            }
        }
        sum.paths_failing = failing;
        sum.dangling_steps = dangling;

        // The reference, if named, is held to the same bounds and must be present: everything
        // downstream is reference-relative, so a reference that did not come back invalidates the lot.
        // Seeding stays richness-driven; this is about RECOVERY, not about which haplotype seeds.
        bool ref_ok = true;
        if (!options.reference_path.empty()) {
            // Exact match wins; otherwise a UNIQUE substring match; an ambiguous one is refused rather
            // than resolved by file order, which is not a property anybody intends to depend on.
            std::size_t ri = g.paths.size();
            for (std::size_t h = 0; h < g.paths.size(); ++h)
                if (g.paths[h].name == options.reference_path) { ri = h; break; }
            if (ri == g.paths.size()) {
                std::vector<std::size_t> hits;
                for (std::size_t h = 0; h < g.paths.size(); ++h)
                    if (g.paths[h].name.find(options.reference_path) != std::string::npos) hits.push_back(h);
                if (hits.size() == 1) ri = hits.front();
                else if (hits.size() > 1)
                    throw std::runtime_error("rebuild: --reference-path '" + options.reference_path +
                                             "' is ambiguous (" + std::to_string(hits.size()) +
                                             " matches, e.g. " + g.paths[hits[0]].name + " and " +
                                             g.paths[hits[1]].name + ")");
            }
            if (ri == g.paths.size()) {
                ref_ok = false;
                sum.reject_reason = "reference path not found in the input: " + options.reference_path;
            } else if (walks[ri].empty()) {
                ref_ok = false;
                sum.reject_reason = "reference path was not recovered: " + g.paths[ri].name;
            } else {
                const double rid = ri < walk_id.size() ? walk_id[ri] : -1.0;
                if (matched[ri] < options.min_matched_cover || rid < 0.0 ||
                    rid < options.min_recovered_identity) {
                    ref_ok = false;
                    sum.reject_reason = "reference path recovered below the contract: " +
                                        g.paths[ri].name + " (identity " + fmt2(rid < 0 ? 0.0 : rid) +
                                        ", matched cover " + fmt2(matched[ri]) + ")";
                }
            }
        }

        if (sum.reject_reason.empty()) {
            if (missing_paths > 0)
                sum.reject_reason = std::to_string(missing_paths) + " haplotype(s) not recovered, first: " + first_bad;
            else if (dangling > 0)
                sum.reject_reason = std::to_string(dangling) + " walk(s) contain a step pair with no edge, first: " + first_bad;
            else if (failing > 0)
                sum.reject_reason = std::to_string(failing) + " haplotype(s) below the contract, first: " + first_bad;
        }
        sum.accepted = ref_ok && sum.reject_reason.empty();

        // ---- audit sidecar: one row per path, so the verdict can be read rather than trusted ----
        const std::string audit = !options.audit_path.empty()
                                      ? options.audit_path
                                      : (options.out_path.empty() ? std::string()
                                                                  : options.out_path + ".rebuild_audit.tsv");
        if (!audit.empty()) {
            const std::filesystem::path audit_staged = staging_path(audit);
            std::ofstream a(audit_staged);
            if (!a) throw std::runtime_error("rebuild: cannot write audit " + audit_staged.string());
            {
                a << "path\toriginal_bp\trecovered_steps\tenvelope_cover\tmatched_cover\tchain_identity"
                     "\twalk_identity\tstatus\n";
                for (std::size_t h = 0; h < g.paths.size(); ++h) {
                    const double id = h < walk_id.size() ? walk_id[h] : -1.0;
                    std::string status = "ok";
                    if (walks[h].empty()) status = "not_recovered";
                    else if (matched[h] < options.min_matched_cover) status = "low_cover";
                    else if (id < 0.0) status = "identity_unavailable";
                    else if (id < options.min_recovered_identity) status = "low_identity";
                    a << g.paths[h].name << '\t' << hseq[h].size() << '\t' << walks[h].size() << '\t'
                      << fmt2(cover[h]) << '\t' << fmt2(matched[h]) << '\t' << fmt2(chain_id[h]) << '\t'
                      << (id < 0.0 ? std::string("NA") : fmt2(id)) << '\t' << status << '\n';
                }
                // The global verdict belongs beside the per-path rows, or a reader has to reconstruct
                // it from them and guess which bound applied.
                // What was DONE with the output, not just whether the contract held: --allow-loss emits
                // the rebuilt graph, and recording that as "rejected" described the opposite of what is
                // on disk.
                a << "#verdict\t"
                  << (sum.accepted ? "accepted"
                                   : (options.allow_loss ? "accepted_with_override" : "rejected"))
                  << '\n';
                a << "#reason\t" << (sum.reject_reason.empty() ? std::string("-") : sum.reject_reason)
                  << '\n';
                a << "#min_recovered_identity\t" << fmt_exact(options.min_recovered_identity) << '\n';
                a << "#min_matched_cover\t" << fmt_exact(options.min_matched_cover) << '\n';
                a << "#dangling_walks\t" << sum.dangling_steps << '\n';
                a.flush();
                if (!a) throw std::runtime_error("rebuild: failed writing audit " + audit_staged.string());
                a.close();
                commit_staged(audit_staged, audit);
                sum.audit_written = true;
            }
        }

        if (!staged_out.empty()) {
            if (sum.accepted || options.allow_loss) {
                if (!sum.accepted)
                    step("WARNING: accepting a rebuild that fails the contract (--allow-loss): " +
                         sum.reject_reason);
                commit_staged(staged_out, options.out_path);
            } else {
                std::error_code ec;
                std::filesystem::remove(staged_out, ec);
                step("REJECTED: " + sum.reject_reason + " -- passing the original graph through unchanged");
                copy_file(options.gfa_path, options.out_path);
            }
        }
    }

    // ---- structural summary of what we emitted ----
    sum.out_nodes = out->n_seg;
    {
        // Per HANDLE, exactly as the input gate measures it. Pooling both ends here while the input is
        // measured per end made the before/after numbers describe different quantities, so any claim
        // that rebuilding reduced tangling was comparing two different rulers.
        std::vector<std::array<std::unordered_set<std::uint32_t>, 2>> hs(out->n_seg);
        std::vector<char> has_loop(out->n_seg, 0);   // NODES with a self-loop, as the input side counts
        std::size_t edges = 0;
        for (std::uint64_t k = 0; k < out->n_arc; ++k) {
            const gfa_arc_t& a = out->arc[k];
            if (a.del || a.comp) continue;
            ++edges;
            const std::uint32_t v = static_cast<std::uint32_t>(a.v_lv >> 32);   // oriented vertex
            const std::uint32_t w = a.w;
            const std::uint32_t vs = v >> 1, ws = w >> 1;
            if (vs == ws) { has_loop[vs] = 1; continue; }
            // v is left by its end when forward, by its start when reverse; w is entered at its start
            // when forward, at its end when reverse.
            hs[vs][(v & 1) ? 0 : 1].insert(ws);
            hs[ws][(w & 1) ? 1 : 0].insert(vs);
        }
        sum.out_edges = edges;
        sum.out_selfloops = static_cast<std::size_t>(std::count(has_loop.begin(), has_loop.end(), 1));
        for (const auto& h : hs) {
            const std::size_t deg = std::max(h[0].size(), h[1].size());
            sum.out_maxdeg = std::max(sum.out_maxdeg, deg);
            if (deg >= options.hub_degree) ++sum.out_hubs;
        }
    }
    if (sum.paths_recovered < sum.haplotypes) {
        std::cerr << "[rebuild " << hms() << "] warning: " << (sum.haplotypes - sum.paths_recovered)
                  << " of " << sum.haplotypes << " haplotypes could not be mapped back and have no P line\n";
    }
    if (!options.quiet) {
        std::cerr << "[rebuild " << hms() << "] emitted " << sum.out_nodes
                  << " nodes, " << sum.out_edges << " edges; maxdeg=" << sum.out_maxdeg << " #hub="
                  << sum.out_hubs << "; paths " << sum.paths_recovered << '/' << sum.haplotypes
                  << ", mean query cover " << sum.mean_query_cover << std::endl;
    }
    return sum;
}

} // namespace panvar
