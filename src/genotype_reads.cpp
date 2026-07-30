#include "panvar/genotype_reads.hpp"

#include "panvar/syncmer.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <fstream>
#include <mutex>
#include <stdexcept>
#include <thread>
#include <unordered_map>

#include <zlib.h>
#include <kseq.h>

KSEQ_INIT(gzFile, gzread)

namespace panvar {
namespace {

std::uint64_t mix64(std::uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

// Must match genotype_markers.cpp exactly, or reads and panel disagree about what an adjacency is.
std::uint64_t adjacency_key(std::uint64_t a, std::uint64_t b, std::size_t gap) {
    const std::uint64_t lo = std::min(a, b);
    const std::uint64_t hi = std::max(a, b);
    const std::uint64_t g = static_cast<std::uint64_t>(std::min<std::size_t>(gap, 0xffff));
    return mix64(mix64(lo) ^ (mix64(hi) * 0x100000001b3ULL) ^ (g * 0x9e3779b97f4a7c15ULL));
}

double median_of(std::vector<double> v) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    return v[v.size() / 2];
}

} // namespace

ReadCounts count_reads(
    const std::vector<std::string>& read_paths,
    const ReadPanel& panel,
    std::size_t threads) {

    std::unordered_map<std::uint64_t, std::uint32_t> node_slot;
    std::unordered_map<std::uint64_t, std::uint32_t> edge_slot;
    node_slot.reserve(panel.node_codes.size() * 2);
    edge_slot.reserve(panel.edge_keys.size() * 2);
    for (std::uint32_t i = 0; i < panel.node_codes.size(); ++i) node_slot.emplace(panel.node_codes[i], i);
    for (std::uint32_t i = 0; i < panel.edge_keys.size(); ++i) edge_slot.emplace(panel.edge_keys[i], i);

    const std::size_t nthreads =
        std::max<std::size_t>(1, threads != 0 ? threads : std::thread::hardware_concurrency());

    ReadCounts total;
    total.node.assign(panel.node_codes.size(), 0);
    total.edge.assign(panel.edge_keys.size(), 0);

    const std::size_t k = panel.kmer_size;
    const std::size_t s = panel.syncmer_s != 0 ? panel.syncmer_s : default_syncmer_s(k);

    std::mutex merge_mu;
    for (const std::string& path : read_paths) {
        gzFile fp = gzopen(path.c_str(), "r");
        if (fp == nullptr) throw std::runtime_error("genotype: cannot open reads " + path);
        kseq_t* seq = kseq_init(fp);

        // One reader (kseq is not thread-safe) feeding batches to workers, each with its own
        // accumulator; `run_parallel` is not used here because it requires distinct output indices.
        std::mutex read_mu;
        std::vector<std::thread> pool;
        for (std::size_t t = 0; t < nthreads; ++t) {
            pool.emplace_back([&] {
                ReadCounts local;
                local.node.assign(panel.node_codes.size(), 0);
                local.edge.assign(panel.edge_keys.size(), 0);
                std::vector<std::string> batch;
                for (;;) {
                    batch.clear();
                    {
                        std::lock_guard<std::mutex> lock(read_mu);
                        while (batch.size() < 4096 && kseq_read(seq) >= 0) {
                            batch.emplace_back(seq->seq.s, seq->seq.l);
                        }
                    }
                    if (batch.empty()) break;
                    for (const std::string& r : batch) {
                        ++local.reads;
                        local.bases += r.size();
                        const std::vector<KmerOccurrence> sy = collect_syncmers(r, k, s);
                        local.syncmers += sy.size();
                        for (const KmerOccurrence& o : sy) {
                            const auto it = node_slot.find(o.code);
                            if (it != node_slot.end()) { ++local.node[it->second]; ++local.matched_syncmers; }
                        }
                        for (std::size_t i = 1; i < sy.size(); ++i) {
                            const std::uint64_t key =
                                adjacency_key(sy[i - 1].code, sy[i].code, sy[i].start - sy[i - 1].start);
                            const auto it = edge_slot.find(key);
                            if (it != edge_slot.end()) ++local.edge[it->second];
                            else ++local.novel_adjacencies;
                        }
                    }
                }
                std::lock_guard<std::mutex> lock(merge_mu);
                for (std::size_t i = 0; i < total.node.size(); ++i) total.node[i] += local.node[i];
                for (std::size_t i = 0; i < total.edge.size(); ++i) total.edge[i] += local.edge[i];
                total.reads += local.reads;
                total.bases += local.bases;
                total.syncmers += local.syncmers;
                total.matched_syncmers += local.matched_syncmers;
                total.novel_adjacencies += local.novel_adjacencies;
            });
        }
        for (std::thread& th : pool) th.join();
        kseq_destroy(seq);
        gzclose(fp);
    }
    return total;
}

std::vector<BlockDepth> estimate_depth(
    const ReadPanel& panel,
    const ReadCounts& counts,
    std::size_t min_anchors,
    double uneven_tolerance) {

    std::vector<BlockDepth> out(panel.anchor_slots.size());
    std::vector<double> all_anchor_counts;

    for (std::size_t bi = 0; bi < panel.anchor_slots.size(); ++bi) {
        BlockDepth& d = out[bi];
        d.block_index = bi;
        std::vector<double> v;
        v.reserve(panel.anchor_slots[bi].size());
        for (const std::uint32_t slot : panel.anchor_slots[bi]) {
            v.push_back(static_cast<double>(counts.node[slot]));
            all_anchor_counts.push_back(static_cast<double>(counts.node[slot]));
        }
        d.n_anchor = v.size();
        if (v.size() < min_anchors) continue;
        d.median = median_of(v);
        std::vector<double> dev;
        dev.reserve(v.size());
        for (const double x : v) dev.push_back(std::fabs(x - d.median));
        d.mad = median_of(dev);
        d.lambda_hap = d.median / 2.0;      // diploid: an invariant marker is carried by both copies
        d.usable = d.lambda_hap > 0.0;
        d.uneven = d.median > 0.0 && (d.mad / d.median) > uneven_tolerance;
    }

    // Region-wide fallback for blocks with too few anchors of their own.
    const double global_median = median_of(all_anchor_counts);
    for (BlockDepth& d : out) {
        if (!d.usable && global_median > 0.0) {
            d.median = global_median;
            d.lambda_hap = global_median / 2.0;
            d.usable = true;
        }
    }
    return out;
}

void write_read_audit(
    const std::string& out_prefix,
    const std::vector<Block>& chain,
    const ReadPanel& panel,
    const ReadCounts& counts,
    const std::vector<BlockDepth>& depth) {

    const std::string path = out_prefix + ".reads.depth.tsv";
    std::ofstream f(path);
    if (!f) throw std::runtime_error("genotype: cannot write " + path);
    f << "block_index\tblock_kind\tbubble_id\tn_alleles\tn_anchor\tanchor_median\tanchor_mad"
         "\tlambda_hap\tusable\tuneven\n";
    for (std::size_t bi = 0; bi < chain.size() && bi < depth.size(); ++bi) {
        const BlockDepth& d = depth[bi];
        f << chain[bi].index << '\t' << (chain[bi].kind == BlockKind::Bubble ? "bubble" : chain[bi].kind == BlockKind::Flank ? "flank" : "backbone")
          << '\t' << chain[bi].bubble_id << '\t'
          << (bi < panel.by_block.size() ? panel.by_block[bi].size() : 0) << '\t' << d.n_anchor << '\t'
          << d.median << '\t' << d.mad << '\t' << d.lambda_hap << '\t' << (d.usable ? 1 : 0) << '\t'
          << (d.uneven ? 1 : 0) << '\n';
    }
    (void)counts;
}

} // namespace panvar
