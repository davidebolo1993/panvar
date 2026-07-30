#include "panvar/genotype_blocks.hpp"

#include "panvar/bubble_alleles.hpp"
#include "panvar/graph_utils.hpp"
#include "panvar/parallel.hpp"
#include "panvar/syncmer.hpp"

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <stdexcept>

namespace panvar {
namespace {

// Group the per-haplotype spelled sequences of one block into alleles, preserving first-encounter
// order so allele indices are stable for a given graph.
BlockAlleles group_by_sequence(
    std::size_t block_index,
    const std::vector<std::string>& hap_names,
    const std::vector<std::string>& seqs,       // parallel to hap_names; empty = does not traverse
    std::size_t n_walk_alleles) {

    BlockAlleles out;
    out.block_index = block_index;
    out.n_walk_alleles = n_walk_alleles;
    std::unordered_map<std::string, std::size_t> seq_to_allele;
    for (std::size_t i = 0; i < hap_names.size(); ++i) {
        if (seqs[i].empty()) continue;
        auto it = seq_to_allele.find(seqs[i]);
        std::size_t ai;
        if (it == seq_to_allele.end()) {
            ai = out.allele_haplotypes.size();
            seq_to_allele.emplace(seqs[i], ai);
            out.allele_haplotypes.push_back(0);
            out.allele_bp.push_back(seqs[i].size());
            out.allele_seq.push_back(seqs[i]);
        } else {
            ai = it->second;
        }
        ++out.allele_haplotypes[ai];
        out.allele_of.emplace(hap_names[i], ai);
    }
    out.n_alleles = out.allele_haplotypes.size();
    return out;
}

} // namespace

std::vector<Block> build_block_chain(const std::vector<Bubble>& bubbles) {
    std::vector<Block> chain;
    if (bubbles.empty()) return chain;

    // Reference order is bubble id order; assert rather than assume, since it only holds for a graph
    // that went through the reference sort.
    for (std::size_t i = 1; i < bubbles.size(); ++i) {
        if (bubbles[i].id < bubbles[i - 1].id) {
            throw std::runtime_error(
                "genotype: bubbles are not in ascending id order; the block chain assumes reference "
                "order, which requires a reference-sorted graph");
        }
    }

    // A bubble's source/sink are not guaranteed to be in reference order -- a bubble whose paths
    // traverse the graph backwards is stored with them swapped (all of c4's and gstm1's are). Taking
    // them at face value makes the neighbouring flank and backbone intervals span straight across the
    // bubble's interior, so the same sequence lands in three blocks at once. Blocks then stop tiling
    // (measured: c4 1.56x, gstm1 1.36x, ankrd36c 1.93x of the haplotype length), and since confinement
    // drops any marker varying in more than one block, every bubble marker is destroyed by its own
    // neighbours. The graph is reference-sorted, so numeric node order IS reference order.
    auto node_pos = [](const std::string& n) {
        char* end = nullptr;
        const unsigned long long v = std::strtoull(n.c_str(), &end, 10);
        return (end != nullptr && *end == '\0') ? v : 0ULL;
    };
    auto left_of = [&](const Bubble& b) {
        return node_pos(b.source) <= node_pos(b.sink) ? b.source : b.sink;
    };
    auto right_of = [&](const Bubble& b) {
        return node_pos(b.source) <= node_pos(b.sink) ? b.sink : b.source;
    };

    // Leading flank: everything upstream of the first bubble. Skipping it leaves ~8.6% of each
    // haplotype (measured on cyp2d6) outside the chain, and reads landing there match no block --
    // which shows up as spurious "novel" adjacencies and would poison the off-panel detector.
    {
        Block b;
        b.index = chain.size();
        b.kind = BlockKind::Flank;
        b.sink = left_of(bubbles.front());   // empty source = from the start of the path
        chain.push_back(b);
    }
    for (std::size_t i = 0; i < bubbles.size(); ++i) {
        if (i > 0) {
            Block b;
            b.index = chain.size();
            b.kind = BlockKind::Backbone;
            b.source = right_of(bubbles[i - 1]);
            b.sink = left_of(bubbles[i]);
            if (b.source != b.sink) chain.push_back(b);   // adjacent bubbles share a boundary node
        }
        Block b;
        b.index = chain.size();
        b.kind = BlockKind::Bubble;
        b.bubble_id = bubbles[i].id;
        b.source = bubbles[i].source;
        b.sink = bubbles[i].sink;
        chain.push_back(b);
    }
    {
        Block b;                              // trailing flank; empty sink = to the end of the path
        b.index = chain.size();
        b.kind = BlockKind::Flank;
        b.source = right_of(bubbles.back());
        chain.push_back(b);
    }
    return chain;
}

BlockAlleles enumerate_block_alleles(
    const Graph& graph,
    const std::vector<BubblePathIndex>& path_indexes,
    const std::vector<Bubble>& bubbles,
    const Block& block,
    std::size_t threads) {

    std::vector<std::string> hap_names(graph.paths.size());
    std::vector<std::string> seqs(graph.paths.size());
    std::vector<std::string> signatures(graph.paths.size());

    const Bubble* bubble = nullptr;
    if (block.kind == BlockKind::Bubble) {
        for (const Bubble& b : bubbles) {
            if (b.id == block.bubble_id) { bubble = &b; break; }
        }
    }

    run_parallel(graph.paths.size(), threads, [&](std::size_t pi) {
        hap_names[pi] = graph.paths[pi].name;
        std::optional<std::vector<PathStep>> steps;
        if (bubble != nullptr) {
            steps = bubble_steps(graph.paths[pi], path_indexes[pi], *bubble);
        } else if (block.kind == BlockKind::Flank) {
            const bool leading = block.source.empty();
            steps = flank_steps(graph.paths[pi], path_indexes[pi],
                                leading ? block.sink : block.source, leading);
        } else {
            steps = interval_steps(graph.paths[pi], path_indexes[pi], block.source, block.sink);
        }
        if (!steps.has_value() || steps->empty()) return;
        seqs[pi] = spell_path_steps_sequence(graph, *steps);
        signatures[pi] = build_walk_signature(*steps);
    });

    std::unordered_set<std::string> walk_sigs;
    for (const std::string& s : signatures) {
        if (!s.empty()) walk_sigs.insert(s);
    }
    return group_by_sequence(block.index, hap_names, seqs, walk_sigs.size());
}

LinkageReport measure_linkage(
    const Graph& graph,
    const std::vector<BlockAlleles>& blocks) {

    LinkageReport rep;
    rep.n_haplotypes = graph.paths.size();
    rep.n_blocks = blocks.size();
    rep.collapse_curve.assign(blocks.size(), 0.0);
    rep.class_size.assign(graph.paths.size(), graph.paths.size());
    if (graph.paths.empty()) return rep;

    // For each haplotype, intersect the set of haplotypes still carrying the same allele, block by
    // block. Blocks a haplotype does not traverse are skipped for that haplotype rather than
    // treated as a mismatch.
    std::vector<std::vector<double>> per_block(blocks.size());
    for (auto& v : per_block) v.assign(graph.paths.size(), 0.0);

    run_parallel(graph.paths.size(), 0, [&](std::size_t hi) {
        const std::string& me = graph.paths[hi].name;
        std::vector<char> compat(graph.paths.size(), 1);
        std::size_t alive = graph.paths.size();
        for (std::size_t bi = 0; bi < blocks.size(); ++bi) {
            const BlockAlleles& B = blocks[bi];
            const auto mine = B.allele_of.find(me);
            if (mine != B.allele_of.end()) {
                alive = 0;
                for (std::size_t hj = 0; hj < graph.paths.size(); ++hj) {
                    if (!compat[hj]) continue;
                    const auto other = B.allele_of.find(graph.paths[hj].name);
                    if (other == B.allele_of.end() || other->second != mine->second) {
                        compat[hj] = 0;
                    } else {
                        ++alive;
                    }
                }
            }
            per_block[bi][hi] = static_cast<double>(alive);
        }
        rep.class_size[hi] = alive;
    });

    for (std::size_t bi = 0; bi < blocks.size(); ++bi) {
        double sum = 0.0;
        for (const double v : per_block[bi]) sum += v;
        rep.collapse_curve[bi] = sum / static_cast<double>(graph.paths.size());
    }
    rep.uniquely_identified =
        static_cast<std::size_t>(std::count(rep.class_size.begin(), rep.class_size.end(), std::size_t{1}));
    return rep;
}

void measure_leave_one_out(
    const Graph& graph,
    const std::vector<Block>& chain,
    const std::vector<BlockAlleles>& blocks,
    LinkageReport& rep) {

    const std::size_t n = graph.paths.size();
    rep.loo_bubble_agreement.assign(n, 0.0);
    if (n < 2) return;

    // Dense allele label per (block, haplotype); -1 when the haplotype does not traverse the block.
    std::vector<std::vector<long>> label(blocks.size(), std::vector<long>(n, -1));
    for (std::size_t bi = 0; bi < blocks.size(); ++bi) {
        for (std::size_t hi = 0; hi < n; ++hi) {
            const auto it = blocks[bi].allele_of.find(graph.paths[hi].name);
            if (it != blocks[bi].allele_of.end()) label[bi][hi] = static_cast<long>(it->second);
        }
    }
    std::vector<std::size_t> bubble_blocks;
    for (std::size_t bi = 0; bi < chain.size(); ++bi) {
        if (chain[bi].kind == BlockKind::Bubble) bubble_blocks.push_back(bi);
    }
    if (bubble_blocks.empty()) return;

    run_parallel(n, 0, [&](std::size_t hi) {
        // Nearest remaining haplotype by agreement over the whole chain.
        std::size_t best = n;
        std::size_t best_hits = 0;
        for (std::size_t hj = 0; hj < n; ++hj) {
            if (hj == hi) continue;
            std::size_t hits = 0;
            for (std::size_t bi = 0; bi < blocks.size(); ++bi) {
                const long a = label[bi][hi];
                if (a >= 0 && label[bi][hj] == a) ++hits;
            }
            if (best == n || hits > best_hits) { best = hj; best_hits = hits; }
        }
        if (best == n) return;
        std::size_t ok = 0;
        std::size_t scored = 0;
        for (const std::size_t bi : bubble_blocks) {
            const long a = label[bi][hi];
            if (a < 0) continue;
            ++scored;
            if (label[bi][best] == a) ++ok;
        }
        rep.loo_bubble_agreement[hi] =
            scored == 0 ? 0.0 : static_cast<double>(ok) / static_cast<double>(scored);
    });

    // Mosaic ceiling: per (haplotype, bubble block), can ANY other panel haplotype supply the same
    // allele? A mosaic model is free to switch between blocks, so this is its representational bound.
    std::vector<double> mosaic(n, 0.0);
    std::size_t singleton_cells = 0;
    std::size_t total_cells = 0;
    for (std::size_t hi = 0; hi < n; ++hi) {
        std::size_t ok = 0;
        std::size_t scored = 0;
        for (const std::size_t bi : bubble_blocks) {
            const long a = label[bi][hi];
            if (a < 0) continue;
            ++scored;
            ++total_cells;
            bool shared = false;
            for (std::size_t hj = 0; hj < n && !shared; ++hj) {
                if (hj != hi && label[bi][hj] == a) shared = true;
            }
            if (shared) ++ok; else ++singleton_cells;
        }
        mosaic[hi] = scored == 0 ? 0.0 : static_cast<double>(ok) / static_cast<double>(scored);
    }
    double msum = 0.0;
    for (const double v : mosaic) msum += v;
    rep.mosaic_ceiling = msum / static_cast<double>(n);
    rep.singleton_mass = total_cells == 0 ? 0.0
                                          : static_cast<double>(singleton_cells) / static_cast<double>(total_cells);

    double sum = 0.0;
    for (const double v : rep.loo_bubble_agreement) sum += v;
    rep.loo_mean_agreement = sum / static_cast<double>(n);
    rep.loo_perfect = static_cast<std::size_t>(
        std::count_if(rep.loo_bubble_agreement.begin(), rep.loo_bubble_agreement.end(),
                      [](double v) { return v >= 1.0; }));
}

void write_linkage_audit(
    const std::string& out_prefix,
    const Graph& graph,
    const std::vector<Block>& chain,
    const std::vector<BlockAlleles>& blocks,
    const LinkageReport& report) {

    const std::string bpath = out_prefix + ".audit.blocks.tsv";
    std::ofstream bf(bpath);
    if (!bf) throw std::runtime_error("genotype: cannot write " + bpath);
    bf << "block_index\tblock_kind\tbubble_id\tsource\tsink\tn_traversing\tn_alleles\tn_walk_alleles"
          "\tmax_allele_haplotypes\tsingleton_alleles\tmedian_allele_bp\tmean_compatible_after\n";
    for (std::size_t i = 0; i < chain.size(); ++i) {
        const BlockAlleles& B = blocks[i];
        std::size_t mx = 0;
        std::size_t singles = 0;
        for (const std::size_t n : B.allele_haplotypes) {
            mx = std::max(mx, n);
            if (n == 1) ++singles;
        }
        std::vector<std::size_t> bp = B.allele_bp;
        std::sort(bp.begin(), bp.end());
        const std::size_t med = bp.empty() ? 0 : bp[bp.size() / 2];
        bf << chain[i].index << '\t'
           << (chain[i].kind == BlockKind::Bubble ? "bubble" : chain[i].kind == BlockKind::Flank ? "flank" : "backbone") << '\t'
           << chain[i].bubble_id << '\t' << chain[i].source << '\t' << chain[i].sink << '\t'
           << B.allele_of.size() << '\t' << B.n_alleles << '\t' << B.n_walk_alleles << '\t'
           << mx << '\t' << singles << '\t' << med << '\t' << report.collapse_curve[i] << '\n';
    }

    const std::string lpath = out_prefix + ".audit.linkage.tsv";
    std::ofstream lf(lpath);
    if (!lf) throw std::runtime_error("genotype: cannot write " + lpath);
    lf << "haplotype\tlinkage_class_size\tuniquely_identified\tloo_bubble_agreement\n";
    for (std::size_t i = 0; i < graph.paths.size(); ++i) {
        lf << graph.paths[i].name << '\t' << report.class_size[i] << '\t'
           << (report.class_size[i] == 1 ? 1 : 0) << '\t'
           << (i < report.loo_bubble_agreement.size() ? report.loo_bubble_agreement[i] : 0.0) << '\n';
    }
}

NoveltyReport measure_novelty(
    const std::vector<Block>& chain,
    const std::vector<BlockAlleles>& blocks,
    std::size_t kmer_size,
    std::size_t syncmer_s) {

    NoveltyReport rep;
    double sum_s = 0.0;
    double sum_a = 0.0;
    const std::size_t s = syncmer_s != 0 ? syncmer_s : default_syncmer_s(kmer_size);

    for (std::size_t bi = 0; bi < chain.size(); ++bi) {
        if (chain[bi].kind != BlockKind::Bubble) continue;
        const BlockAlleles& B = blocks[bi];
        if (B.allele_seq.size() < 2) continue;

        // Inverted index: how many alleles carry each syncmer / adjacency. Counting once is O(total
        // markers); the naive alternative rescans every sibling per marker and dominated the runtime.
        std::vector<std::vector<std::uint64_t>> nodes(B.allele_seq.size());
        std::vector<std::vector<std::uint64_t>> edges(B.allele_seq.size());
        std::unordered_map<std::uint64_t, std::uint32_t> node_owners;
        std::unordered_map<std::uint64_t, std::uint32_t> edge_owners;
        for (std::size_t ai = 0; ai < B.allele_seq.size(); ++ai) {
            const std::vector<KmerOccurrence> sy = collect_syncmers(B.allele_seq[ai], kmer_size, s);
            std::unordered_set<std::uint64_t> n_seen;
            std::unordered_set<std::uint64_t> e_seen;
            for (std::size_t i = 0; i < sy.size(); ++i) {
                n_seen.insert(sy[i].code);
                if (i > 0) {
                    e_seen.insert(sy[i - 1].code ^ (sy[i].code * 0x9e3779b97f4a7c15ULL));
                }
            }
            nodes[ai].assign(n_seen.begin(), n_seen.end());
            edges[ai].assign(e_seen.begin(), e_seen.end());
            for (const std::uint64_t c : nodes[ai]) ++node_owners[c];
            for (const std::uint64_t c : edges[ai]) ++edge_owners[c];
        }
        for (std::size_t ai = 0; ai < B.allele_seq.size(); ++ai) {
            if (B.allele_haplotypes[ai] != 1) continue;
            ++rep.private_alleles;
            std::size_t seen_n = 0;
            std::size_t seen_e = 0;
            for (const std::uint64_t c : nodes[ai]) if (node_owners[c] > 1) ++seen_n;
            for (const std::uint64_t c : edges[ai]) if (edge_owners[c] > 1) ++seen_e;
            const double fn = nodes[ai].empty() ? 1.0
                                                : static_cast<double>(seen_n) / static_cast<double>(nodes[ai].size());
            const double fe = edges[ai].empty() ? 1.0
                                                : static_cast<double>(seen_e) / static_cast<double>(edges[ai].size());
            sum_s += fn;
            sum_a += fe;
            if (fn >= 1.0) ++rep.fully_reusable; else ++rep.with_novel_sequence;
        }
    }
    if (rep.private_alleles > 0) {
        rep.mean_syncmer_reuse = sum_s / static_cast<double>(rep.private_alleles);
        rep.mean_adjacency_reuse = sum_a / static_cast<double>(rep.private_alleles);
    }
    return rep;
}

} // namespace panvar
