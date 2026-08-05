#include "panvar/genotype_markers.hpp"

#include "panvar/graph_utils.hpp"
#include "panvar/parallel.hpp"
#include "panvar/syncmer.hpp"

#include <algorithm>
#include <chrono>
#include <fstream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <mutex>

namespace panvar {
namespace {

std::uint64_t mix64(std::uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

// Key for an adjacency between two consecutive syncmers. Order-normalized so a read and its
// reverse-complement agree, and the base gap is folded in: the same syncmer pair separated by a
// different distance is a different junction, which is exactly the signal an indel produces.
std::uint64_t adjacency_key(std::uint64_t a, std::uint64_t b, std::size_t gap) {
    const std::uint64_t lo = std::min(a, b);
    const std::uint64_t hi = std::max(a, b);
    const std::uint64_t g = static_cast<std::uint64_t>(std::min<std::size_t>(gap, 0xffff));
    return mix64(mix64(lo) ^ (mix64(hi) * 0x100000001b3ULL) ^ (g * 0x9e3779b97f4a7c15ULL));
}

using Counts = std::unordered_map<std::uint64_t, std::uint32_t>;

// Per-allele syncmer inventory for one bubble.
struct AlleleInventory {
    Counts nodes;
    Counts edges;
    std::size_t bp = 0;
    std::size_t syncmers_total = 0;
};

void collect_inventory(const std::string& seq, std::size_t k, std::size_t s, AlleleInventory& inv,
                       bool all_kmers = false) {
    const std::vector<KmerOccurrence> sy =
        all_kmers ? collect_canonical_kmer_occurrences(seq, k) : collect_syncmers(seq, k, s);
    inv.bp = seq.size();
    inv.syncmers_total = sy.size();
    for (const KmerOccurrence& o : sy) ++inv.nodes[o.code];
    for (std::size_t i = 1; i < sy.size(); ++i) {
        ++inv.edges[adjacency_key(sy[i - 1].code, sy[i].code, sy[i].start - sy[i - 1].start)];
    }
}

double jaccard(const Counts& a, const Counts& b) {
    if (a.empty() && b.empty()) return 1.0;
    std::size_t inter = 0;
    const Counts& small = a.size() <= b.size() ? a : b;
    const Counts& large = a.size() <= b.size() ? b : a;
    for (const auto& [c, n] : small) { (void)n; if (large.count(c) != 0) ++inter; }
    const std::size_t uni = a.size() + b.size() - inter;
    return uni == 0 ? 1.0 : static_cast<double>(inter) / static_cast<double>(uni);
}

// Separation between the alleles of one block, in one marker unit. A marker is informative when its
// multiplicity is not constant across alleles; presence/absence is the 0-vs-n special case, and equal
// presence at different copy number is the tandem case a private-set rule misses. Then
//   differ(a,b) = carried[a] + carried[b] - both(a,b) - equal(a,b)
// and an allele's usable separation is its minimum over siblings -- the hardest one to tell apart.
void compute_separation(
    const std::vector<AlleleInventory>& inv,
    bool edges,
    std::vector<std::size_t>& carried,
    std::vector<std::size_t>& min_sep,
    std::size_t& n_informative,
    std::size_t max_dense_alleles,
    std::size_t top_k,
    std::size_t sketch_size) {

    const std::size_t n = inv.size();
    // code -> the alleles carrying it, with multiplicity.
    std::unordered_map<std::uint64_t, std::vector<std::pair<std::uint32_t, std::uint32_t>>> present;
    for (std::size_t ai = 0; ai < n; ++ai) {
        const Counts& src = edges ? inv[ai].edges : inv[ai].nodes;
        for (const auto& [c, mult] : src) {
            present[c].emplace_back(static_cast<std::uint32_t>(ai), mult);
        }
    }

    carried.assign(n, 0);
    min_sep.assign(n, 0);
    n_informative = 0;

    // Informative = multiplicity not constant across alleles. Constant markers cannot separate
    // anything, whatever their count, so they are dropped before the pairwise work.
    std::vector<const std::vector<std::pair<std::uint32_t, std::uint32_t>>*> informative;
    std::vector<std::uint64_t> informative_codes;
    informative.reserve(present.size());
    informative_codes.reserve(present.size());
    for (const auto& [c, plist] : present) {
        bool inf = plist.size() < n;
        if (!inf) {
            for (std::size_t i = 1; i < plist.size(); ++i) {
                if (plist[i].second != plist[0].second) { inf = true; break; }
            }
        }
        if (!inf) continue;
        ++n_informative;
        informative.push_back(&plist);
        informative_codes.push_back(c);
        for (const auto& [ai, mult] : plist) { (void)mult; ++carried[ai]; }
    }
    if (n <= 1) {
        if (n == 1) min_sep[0] = carried[0];
        return;
    }

    // Per-allele marker lists (sparse path only), so each row can be accumulated independently.
    std::vector<std::vector<std::pair<const std::vector<std::pair<std::uint32_t, std::uint32_t>>*,
                                      std::uint32_t>>> by_allele(n);
    std::vector<std::vector<std::uint64_t>> codes_of(n);   // informative marker codes per allele
    for (std::size_t idx = 0; idx < informative.size(); ++idx) {
        const auto* plist = informative[idx];
        for (const auto& [ai, mult] : *plist) {
            by_allele[ai].emplace_back(plist, mult);
            codes_of[ai].push_back(informative_codes[idx]);
        }
    }

    // Alleles ordered by how many informative markers they carry: a sibling sharing nothing with `a`
    // gives differ = carried[a] + carried[b], so the best of those is the one carrying the fewest.
    std::vector<std::uint32_t> by_carried(n);
    for (std::size_t i = 0; i < n; ++i) by_carried[i] = static_cast<std::uint32_t>(i);
    std::sort(by_carried.begin(), by_carried.end(), [&](std::uint32_t x, std::uint32_t y) {
        if (carried[x] != carried[y]) return carried[x] < carried[y];
        return x < y;
    });

    if (top_k > 0 && n > top_k + 1) {
        // NOTE the error is one-directional and unsafe: min over a subset >= min over all, so this
        // path can only ever OVERSTATE how separable an allele is, never understate it. Measured on
        // the test loci at K=32 it overstates 0-8 block metrics per locus (worst deviation 11) and at
        // K=8 up to 16 (worst 79). Treat its output as optimistic and calibrate K against an exact
        // run before relying on it.
        // Bottom-k MinHash over the informative markers each allele carries. Marker codes are 2-bit
        // packed sequence, not hashes, so they must be mixed before taking minima or the sketch is
        // biased toward low-complexity k-mers.
        std::vector<std::vector<std::uint64_t>> sketch(n);
        for (std::size_t ai = 0; ai < n; ++ai) {
            std::vector<std::uint64_t> h;
            h.reserve(codes_of[ai].size());
            for (const std::uint64_t code : codes_of[ai]) h.push_back(mix64(code));
            std::sort(h.begin(), h.end());
            h.erase(std::unique(h.begin(), h.end()), h.end());
            if (h.size() > sketch_size) h.resize(sketch_size);
            sketch[ai] = std::move(h);
        }
        std::unordered_map<std::uint64_t, std::vector<std::uint32_t>> sketch_index;
        for (std::size_t ai = 0; ai < n; ++ai) {
            for (const std::uint64_t v : sketch[ai]) {
                sketch_index[v].push_back(static_cast<std::uint32_t>(ai));
            }
        }

        std::vector<std::uint32_t> hits(n, 0);
        std::vector<std::uint32_t> touched;
        std::vector<std::uint32_t> both(n, 0);
        std::vector<std::uint32_t> equal(n, 0);
        std::vector<std::uint32_t> touched2;
        for (std::size_t ai = 0; ai < n; ++ai) {
            for (const std::uint32_t bi : touched) hits[bi] = 0;
            touched.clear();
            for (const std::uint64_t v : sketch[ai]) {
                const auto it = sketch_index.find(v);
                if (it == sketch_index.end()) continue;
                for (const std::uint32_t bi : it->second) {
                    if (bi == ai) continue;
                    if (hits[bi] == 0) touched.push_back(bi);
                    ++hits[bi];
                }
            }
            std::vector<std::uint32_t> cand = touched;
            // Total order: ties on hit count must break on allele index, or `resize` below keeps an
            // arbitrary subset and the whole path stops being reproducible run to run.
            std::sort(cand.begin(), cand.end(), [&](std::uint32_t x, std::uint32_t y) {
                if (hits[x] != hits[y]) return hits[x] > hits[y];
                return x < y;
            });
            if (cand.size() > top_k) cand.resize(top_k);
            // differ(a,b) >= |carried[a]-carried[b]|, so the siblings that could still beat a sketch
            // miss are the ones closest in marker count. Expand outward from a's own rank.
            {
                const auto pos = std::lower_bound(by_carried.begin(), by_carried.end(), carried[ai],
                    [&](std::uint32_t x, std::size_t v) { return carried[x] < v; }) - by_carried.begin();
                for (std::size_t step = 0; step < top_k && cand.size() < top_k * 2; ++step) {
                    const long lo = static_cast<long>(pos) - static_cast<long>(step) - 1;
                    const std::size_t hi = pos + step;
                    if (lo >= 0 && by_carried[lo] != ai) cand.push_back(by_carried[lo]);
                    if (hi < by_carried.size() && by_carried[hi] != ai) cand.push_back(by_carried[hi]);
                }
            }
            std::sort(cand.begin(), cand.end());
            cand.erase(std::unique(cand.begin(), cand.end()), cand.end());

            for (const std::uint32_t bi : touched2) { both[bi] = 0; equal[bi] = 0; }
            touched2.clear();
            for (const auto& [plist, my_mult] : by_allele[ai]) {
                for (const auto& [bi, mult] : *plist) {
                    if (bi == ai) continue;
                    if (both[bi] == 0 && equal[bi] == 0) touched2.push_back(bi);
                    ++both[bi];
                    if (mult == my_mult) ++equal[bi];
                }
            }
            std::size_t best = static_cast<std::size_t>(-1);
            for (const std::uint32_t bi : cand) {
                best = std::min(best, carried[ai] + carried[bi] - both[bi] - equal[bi]);
            }
            min_sep[ai] = (best == static_cast<std::size_t>(-1)) ? carried[ai] : best;
        }
        return;
    }

    if (n <= max_dense_alleles) {
        // Dense pairwise scratch: fastest, but O(n^2) memory. Used whenever the block is small
        // enough that the matrix is cheap, which is every block in the current test loci.
        std::vector<std::uint32_t> both(n * n, 0);
        std::vector<std::uint32_t> equal(n * n, 0);
        for (const auto* plist : informative) {
            for (std::size_t i = 0; i < plist->size(); ++i) {
                for (std::size_t j = i + 1; j < plist->size(); ++j) {
                    const std::size_t x = (*plist)[i].first * n + (*plist)[j].first;
                    const std::size_t y = (*plist)[j].first * n + (*plist)[i].first;
                    ++both[x]; ++both[y];
                    if ((*plist)[i].second == (*plist)[j].second) { ++equal[x]; ++equal[y]; }
                }
            }
        }
        for (std::size_t ai = 0; ai < n; ++ai) {
            std::size_t worst = static_cast<std::size_t>(-1);
            for (std::size_t aj = 0; aj < n; ++aj) {
                if (aj == ai) continue;
                const std::size_t idx = ai * n + aj;
                worst = std::min(worst, carried[ai] + carried[aj] - both[idx] - equal[idx]);
            }
            min_sep[ai] = worst;
        }
        return;
    }

    // Sparse accumulator: one row at a time into a dense scratch vector, reset only where touched.
    // Memory is O(n_alleles) instead of the O(n_alleles^2) a full pairwise matrix would need, which
    // is what makes large panels feasible; results are identical to a dense pass. Kept serial on
    // purpose -- the caller already parallelises over blocks, and nesting the two oversubscribes.
    std::vector<std::uint32_t> both(n, 0);
    std::vector<std::uint32_t> equal(n, 0);
    std::vector<std::uint32_t> touched;
    touched.reserve(n);
    for (std::size_t ai = 0; ai < n; ++ai) {
        for (const std::uint32_t bi : touched) { both[bi] = 0; equal[bi] = 0; }
        touched.clear();
        for (const auto& [plist, my_mult] : by_allele[ai]) {
            for (const auto& [bi, mult] : *plist) {
                if (bi == ai) continue;
                if (both[bi] == 0 && equal[bi] == 0) touched.push_back(bi);
                ++both[bi];
                if (mult == my_mult) ++equal[bi];
            }
        }
        std::size_t best = static_cast<std::size_t>(-1);
        for (const std::uint32_t bi : touched) {
            best = std::min(best, carried[ai] + carried[bi] - both[bi] - equal[bi]);
        }
        for (const std::uint32_t bi : by_carried) {          // cheapest non-overlapping sibling
            if (bi == ai) continue;
            if (both[bi] != 0 || equal[bi] != 0) continue;
            best = std::min(best, carried[ai] + carried[bi]);
            break;
        }
        min_sep[ai] = (best == static_cast<std::size_t>(-1)) ? carried[ai] : best;
    }
}
} // namespace

MarkerPanel build_marker_panel(
    const Graph& graph,
    const std::vector<Bubble>& bubbles,
    const std::string& reference_path,
    const MarkerOptions& options) {

    const std::size_t k = options.kmer_size;
    const std::size_t s = options.syncmer_s != 0 ? options.syncmer_s : default_syncmer_s(k);

    MarkerPanel panel;
    panel.kmer_size = k;
    panel.syncmer_s = s;

    std::vector<BubblePathIndex> path_indexes(graph.paths.size());
    run_parallel(graph.paths.size(), options.threads, [&](std::size_t i) {
        path_indexes[i] = build_bubble_path_index(graph.paths[i]);
    });

    const std::unordered_set<std::string> selfloops = self_loop_nodes(graph);

    // ---- pass 1: per-bubble syncmer graphs, private-within-bubble candidates ----
    struct BubbleWork {
        BubbleMarkerReport report;
        std::vector<AlleleMarkers> markers;
        std::vector<Counts> cand_nodes;   // per allele: private node -> multiplicity in the walk
        std::vector<Counts> cand_edges;
    };
    std::vector<BubbleWork> work(bubbles.size());

    run_parallel(bubbles.size(), options.threads, [&](std::size_t bi) {
        const Bubble& bubble = bubbles[bi];
        BubbleWork& w = work[bi];
        BubbleAlleleSet set = enumerate_bubble_alleles(graph, path_indexes, bubble, reference_path);

        w.report.bubble_id = bubble.id;
        w.report.source = bubble.source;
        w.report.sink = bubble.sink;
        w.report.n_inside = bubble.inside.size();
        w.report.n_alleles = set.alleles.size();
        w.report.n_traversing = set.traversing.size();
        w.report.has_reference = set.has_reference;
        for (const std::string& n : bubble.inside) {
            if (selfloops.count(n) != 0) { w.report.folded = true; break; }
        }
        if (set.alleles.empty()) return;

        // Group walks by the sequence they spell. Two walks that spell the same bases are the same
        // allele to any read-based method, however differently the graph routes them, and keeping
        // them apart would leave both with no private markers.
        std::vector<std::string> seqs(set.alleles.size());
        for (std::size_t ai = 0; ai < set.alleles.size(); ++ai) {
            seqs[ai] = spell_path_steps_sequence(graph, set.alleles[ai].steps);
        }
        std::unordered_map<std::string, std::size_t> seq_to_group;
        std::vector<std::size_t> group_of(set.alleles.size(), 0);
        std::vector<std::size_t> group_first;              // representative walk index
        std::vector<std::size_t> group_haps;
        for (std::size_t ai = 0; ai < set.alleles.size(); ++ai) {
            auto it = seq_to_group.find(seqs[ai]);
            if (it == seq_to_group.end()) {
                const std::size_t gi = group_first.size();
                seq_to_group.emplace(seqs[ai], gi);
                group_of[ai] = gi;
                group_first.push_back(ai);
                group_haps.push_back(set.alleles[ai].members.size());
            } else {
                group_of[ai] = it->second;
                group_haps[it->second] += set.alleles[ai].members.size();
            }
        }
        w.report.n_alleles_walk = set.alleles.size();
        w.report.n_alleles = group_first.size();

        std::vector<AlleleInventory> inv(group_first.size());
        for (std::size_t gi = 0; gi < group_first.size(); ++gi) {
            collect_inventory(seqs[group_first[gi]], k, s, inv[gi]);
        }

        std::size_t ref_group = group_first.size();
        if (set.has_reference) {
            for (std::size_t ai = 0; ai < set.alleles.size(); ++ai) {
                if (set.alleles[ai].signature == set.reference_signature) { ref_group = group_of[ai]; break; }
            }
        }
        const std::size_t n_groups = group_first.size();

        // Informative markers: those whose multiplicity is NOT constant across the bubble's alleles.
        // Constant markers say nothing about which allele a read came from, whatever their count.
        // For each informative marker, `separating[a][b]` counts how many distinguish a from b:
        //   differ(a,b) = |A| + |B| - both(a,b) - equal(a,b)
        // where |A| is the informative markers a carries, both(a,b) those carried by each, and
        // equal(a,b) those carried by each at the SAME multiplicity.
        std::vector<std::size_t> carried_n, carried_e, minsep_n, minsep_e;
        std::size_t informative_n = 0;
        std::size_t informative_e = 0;
        compute_separation(inv, false, carried_n, minsep_n, informative_n, options.max_dense_alleles, options.sep_top_k, options.sketch_size);
        compute_separation(inv, true, carried_e, minsep_e, informative_e, options.max_dense_alleles, options.sep_top_k, options.sketch_size);

        std::size_t singleton_haps = 0;
        w.markers.resize(n_groups);
        w.cand_nodes.resize(n_groups);
        w.cand_edges.resize(n_groups);
        for (std::size_t ai = 0; ai < n_groups; ++ai) {
            AlleleMarkers& m = w.markers[ai];
            m.bubble_id = bubble.id;
            m.allele_id = ai;
            m.n_haplotypes = group_haps[ai];
            m.allele_bp = inv[ai].bp;
            m.is_reference = ai == ref_group;
            m.n_syncmers_total = inv[ai].syncmers_total;
            m.n_syncmers_distinct = inv[ai].nodes.size();
            m.n_informative_nodes = informative_n;
            m.n_informative_edges = informative_e;
            m.n_carried_nodes = carried_n[ai];
            m.n_carried_edges = carried_e[ai];
            m.min_separating_nodes = minsep_n[ai];
            m.min_separating_edges = minsep_e[ai];
            if (m.n_haplotypes == 1) { ++w.report.singleton_alleles; singleton_haps += 1; }

            for (const auto& [c, mult] : inv[ai].nodes) {
                if (mult > options.max_multiplicity) continue;
                w.cand_nodes[ai].emplace(c, mult);
            }
            for (const auto& [c, mult] : inv[ai].edges) {
                if (mult > options.max_multiplicity) continue;
                w.cand_edges[ai].emplace(c, mult);
            }

            double best = 0.0;
            std::size_t best_j = ai;
            for (std::size_t aj = 0; aj < n_groups; ++aj) {
                if (aj == ai) continue;
                const double j = jaccard(inv[ai].nodes, inv[aj].nodes);
                if (j > best) { best = j; best_j = aj; }
            }
            m.hardest_sibling = best_j;
            m.nearest_sibling_jaccard = best;
        }
        w.report.singleton_hap_frac =
            w.report.n_traversing == 0 ? 0.0
                                       : static_cast<double>(singleton_haps) /
                                             static_cast<double>(w.report.n_traversing);
    });

    // ---- pass 2: region uniqueness ----
    // A candidate is usable only if it occurs, across the whole panel, no more often than the
    // alleles of its own bubble account for. Anything extra is an occurrence somewhere else in the
    // region, where reads cannot be attributed to this allele.
    std::unordered_map<std::uint64_t, std::uint64_t> expected;   // candidate -> occurrences owed
    std::unordered_map<std::uint64_t, std::size_t> owner_bubble; // candidate -> bubble (0 = clash)
    if (options.require_region_unique) {
        for (std::size_t bi = 0; bi < work.size(); ++bi) {
            const BubbleWork& w = work[bi];
            for (std::size_t ai = 0; ai < w.markers.size(); ++ai) {
                const std::uint64_t haps = w.markers[ai].n_haplotypes;
                for (const auto& [c, mult] : w.cand_nodes[ai]) {
                    expected[c] += haps * mult;
                    auto it = owner_bubble.find(c);
                    if (it == owner_bubble.end()) owner_bubble.emplace(c, w.report.bubble_id);
                    else if (it->second != w.report.bubble_id) it->second = 0;
                }
                for (const auto& [c, mult] : w.cand_edges[ai]) {
                    expected[c] += haps * mult;
                    auto it = owner_bubble.find(c);
                    if (it == owner_bubble.end()) owner_bubble.emplace(c, w.report.bubble_id);
                    else if (it->second != w.report.bubble_id) it->second = 0;
                }
            }
        }

        // Dense slots so the per-path count pass needs no hashing of misses.
        std::unordered_map<std::uint64_t, std::uint32_t> slot;
        slot.reserve(expected.size() * 2);
        for (const auto& [c, n] : expected) { (void)n; const std::uint32_t id = static_cast<std::uint32_t>(slot.size()); slot.emplace(c, id); }

        const std::size_t nthreads =
            std::max<std::size_t>(1, options.threads != 0 ? options.threads
                                                          : std::thread::hardware_concurrency());
        std::vector<std::vector<std::uint64_t>> partial(nthreads,
                                                        std::vector<std::uint64_t>(slot.size(), 0));
        std::atomic<std::size_t> next_thread{0};
        std::vector<std::size_t> thread_of(nthreads);
        (void)thread_of;
        std::atomic<std::size_t> cursor{0};
        std::vector<std::thread> pool;
        for (std::size_t t = 0; t < nthreads; ++t) {
            pool.emplace_back([&, t] {
                std::vector<std::uint64_t>& acc = partial[t];
                for (;;) {
                    const std::size_t pi = cursor.fetch_add(1);
                    if (pi >= graph.paths.size()) break;
                    const std::string seq = spell_path_steps_sequence(graph, graph.paths[pi].steps);
                    const std::vector<KmerOccurrence> sy = collect_syncmers(seq, k, s);
                    for (const KmerOccurrence& o : sy) {
                        const auto it = slot.find(o.code);
                        if (it != slot.end()) ++acc[it->second];
                    }
                    for (std::size_t i = 1; i < sy.size(); ++i) {
                        const std::uint64_t key =
                            adjacency_key(sy[i - 1].code, sy[i].code, sy[i].start - sy[i - 1].start);
                        const auto it = slot.find(key);
                        if (it != slot.end()) ++acc[it->second];
                    }
                }
            });
        }
        for (std::thread& th : pool) th.join();
        (void)next_thread;

        std::vector<std::uint64_t> actual(slot.size(), 0);
        for (const auto& part : partial) {
            for (std::size_t i = 0; i < actual.size(); ++i) actual[i] += part[i];
        }

        for (BubbleWork& w : work) {
            for (std::size_t ai = 0; ai < w.markers.size(); ++ai) {
                AlleleMarkers& m = w.markers[ai];
                for (auto it = w.cand_nodes[ai].begin(); it != w.cand_nodes[ai].end();) {
                    const std::uint64_t c = it->first;
                    const bool clash = owner_bubble[c] == 0;
                    if (clash || actual[slot[c]] > expected[c]) {
                        ++m.n_nodes_lost_region;
                        it = w.cand_nodes[ai].erase(it);
                    } else {
                        ++it;
                    }
                }
                for (auto it = w.cand_edges[ai].begin(); it != w.cand_edges[ai].end();) {
                    const std::uint64_t c = it->first;
                    const bool clash = owner_bubble[c] == 0;
                    if (clash || actual[slot[c]] > expected[c]) {
                        ++m.n_edges_lost_region;
                        it = w.cand_edges[ai].erase(it);
                    } else {
                        ++it;
                    }
                }
            }
        }
    }

    // ---- finalize ----
    for (BubbleWork& w : work) {
        panel.bubbles.push_back(w.report);
        for (std::size_t ai = 0; ai < w.markers.size(); ++ai) {
            AlleleMarkers& m = w.markers[ai];
            m.node_codes.reserve(w.cand_nodes[ai].size());
            for (const auto& [c, n] : w.cand_nodes[ai]) { (void)n; m.node_codes.push_back(c); }
            m.edge_keys.reserve(w.cand_edges[ai].size());
            for (const auto& [c, n] : w.cand_edges[ai]) { (void)n; m.edge_keys.push_back(c); }
            std::sort(m.node_codes.begin(), m.node_codes.end());
            std::sort(m.edge_keys.begin(), m.edge_keys.end());
            panel.markers.push_back(std::move(m));
        }
    }
    return panel;
}

void write_marker_audit(const std::string& out_prefix, const MarkerPanel& panel) {
    const std::string bpath = out_prefix + ".audit.bubbles.tsv";
    std::ofstream bf(bpath);
    if (!bf) throw std::runtime_error("genotype: cannot write " + bpath);
    bf << "bubble_id\tsource\tsink\tn_inside\tn_alleles\tn_alleles_walk\tn_traversing\tsingleton_alleles"
          "\tsingleton_hap_frac\tfolded\thas_reference\n";
    for (const BubbleMarkerReport& b : panel.bubbles) {
        bf << b.bubble_id << '\t' << b.source << '\t' << b.sink << '\t' << b.n_inside << '\t'
           << b.n_alleles << '\t' << b.n_alleles_walk << '\t' << b.n_traversing << '\t' << b.singleton_alleles << '\t'
           << b.singleton_hap_frac << '\t' << (b.folded ? 1 : 0) << '\t' << (b.has_reference ? 1 : 0)
           << '\n';
    }

    const std::string apath = out_prefix + ".audit.alleles.tsv";
    std::ofstream af(apath);
    if (!af) throw std::runtime_error("genotype: cannot write " + apath);
    af << "bubble_id\tallele_id\tis_reference\tn_haplotypes\tallele_bp\tn_syncmers_total"
          "\tn_syncmers_distinct\tn_informative_nodes\tn_informative_edges\tn_carried_nodes"
          "\tn_carried_edges\tmin_separating_nodes\tmin_separating_edges\tn_nodes_lost_region"
          "\tn_edges_lost_region\thardest_sibling\tnearest_sibling_jaccard\n";
    for (const AlleleMarkers& m : panel.markers) {
        af << m.bubble_id << '\t' << m.allele_id << '\t' << (m.is_reference ? 1 : 0) << '\t'
           << m.n_haplotypes << '\t' << m.allele_bp << '\t' << m.n_syncmers_total << '\t'
           << m.n_syncmers_distinct << '\t' << m.n_informative_nodes << '\t'
           << m.n_informative_edges << '\t' << m.n_carried_nodes << '\t' << m.n_carried_edges
           << '\t' << m.min_separating_nodes << '\t' << m.min_separating_edges << '\t'
           << m.n_nodes_lost_region << '\t' << m.n_edges_lost_region << '\t' << m.hardest_sibling
           << '\t' << m.nearest_sibling_jaccard << '\n';
    }
}

std::vector<BlockMarkerStats> build_block_marker_panel(
    const std::vector<Block>& chain,
    const std::vector<BlockAlleles>& blocks,
    const MarkerOptions& options,
    ReadPanel* out_panel,
    const Graph* graph_for_region_uniqueness,
    bool want_separation_stats) {

    const std::size_t k = options.kmer_size;
    const std::size_t s = options.syncmer_s != 0 ? options.syncmer_s : default_syncmer_s(k);
    std::vector<BlockMarkerStats> out(chain.size());
    std::vector<std::vector<AlleleInventory>> kept_inv(out_panel ? chain.size() : 0);

    run_parallel(chain.size(), options.threads, [&](std::size_t bi) {
        const auto t0 = std::chrono::steady_clock::now();
        const BlockAlleles& B = blocks[bi];
        BlockMarkerStats& st = out[bi];
        st.block_index = chain[bi].index;
        st.is_bubble = chain[bi].kind == BlockKind::Bubble;
        st.bubble_id = chain[bi].bubble_id;
        st.n_alleles = B.allele_seq.size();
        st.n_haplotypes = B.allele_of.size();
        if (B.allele_seq.empty()) return;

        std::vector<std::size_t> bp = B.allele_bp;
        std::sort(bp.begin(), bp.end());
        st.block_bp = bp[bp.size() / 2];

        std::vector<AlleleInventory> inv(B.allele_seq.size());
        for (std::size_t ai = 0; ai < B.allele_seq.size(); ++ai) {
            collect_inventory(B.allele_seq[ai], k, s, inv[ai], options.all_kmers);
            for (const auto& [c, mult] : inv[ai].nodes) {
                (void)c;
                st.max_marker_multiplicity_seen = std::max<std::size_t>(st.max_marker_multiplicity_seen, mult);
                if (options.max_multiplicity != 0 && mult > options.max_multiplicity) ++st.markers_over_cap;
            }
            // A cap of 0 means "no cap": with multiplicity modelled explicitly, a marker repeated N
            // times is the copy-number signal rather than noise, so capping discards exactly the
            // evidence a tandem array carries.
            if (options.max_multiplicity != 0) {
                for (auto it = inv[ai].nodes.begin(); it != inv[ai].nodes.end();) {
                    it = (it->second > options.max_multiplicity) ? inv[ai].nodes.erase(it) : std::next(it);
                }
                for (auto it = inv[ai].edges.begin(); it != inv[ai].edges.end();) {
                    it = (it->second > options.max_multiplicity) ? inv[ai].edges.erase(it) : std::next(it);
                }
            }
        }

        // Pairwise separation is an audit statistic, not something the genotyper needs. It is the
        // single most expensive part of building the panel, so genotyping skips it.
        std::vector<std::size_t> carried_n, carried_e, sep_n, sep_e;
        if (want_separation_stats) {
            compute_separation(inv, false, carried_n, sep_n, st.n_informative_nodes, options.max_dense_alleles, options.sep_top_k, options.sketch_size);
            compute_separation(inv, true, carried_e, sep_e, st.n_informative_edges, options.max_dense_alleles, options.sep_top_k, options.sketch_size);
        } else {
            sep_n.assign(B.allele_seq.size(), 0);
            sep_e.assign(B.allele_seq.size(), 0);
        }

        std::vector<std::size_t> sn = sep_n;
        std::vector<std::size_t> se = sep_e;
        std::sort(sn.begin(), sn.end());
        std::sort(se.begin(), se.end());
        st.median_sep_nodes = sn.empty() ? 0 : sn[sn.size() / 2];
        st.median_sep_edges = se.empty() ? 0 : se[se.size() / 2];
        st.unseparable_alleles =
            static_cast<std::size_t>(std::count(sep_e.begin(), sep_e.end(), std::size_t{0}));

        std::size_t mass_n = 0;
        std::size_t mass_e = 0;
        std::size_t total = 0;
        for (std::size_t ai = 0; ai < B.allele_haplotypes.size(); ++ai) {
            total += B.allele_haplotypes[ai];
            if (sep_n[ai] >= options.min_markers) mass_n += B.allele_haplotypes[ai];
            if (sep_e[ai] >= options.min_markers) mass_e += B.allele_haplotypes[ai];
        }
        if (total > 0) {
            st.separable_mass_nodes = static_cast<double>(mass_n) / static_cast<double>(total);
            st.separable_mass_edges = static_cast<double>(mass_e) / static_cast<double>(total);
        }
        if (out_panel != nullptr) kept_inv[bi] = std::move(inv);
        st.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    });

    if (out_panel != nullptr) {
        out_panel->kmer_size = k;
        out_panel->syncmer_s = s;
        out_panel->all_kmers = options.all_kmers;
        out_panel->by_block.assign(chain.size(), {});
        out_panel->anchor_slots.assign(chain.size(), {});
        std::unordered_map<std::uint64_t, std::uint32_t> node_slot;
        std::unordered_map<std::uint64_t, std::uint32_t> edge_slot;
        auto slot_for = [](std::unordered_map<std::uint64_t, std::uint32_t>& m,
                           std::vector<std::uint64_t>& codes, std::uint64_t c) {
            auto it = m.find(c);
            if (it != m.end()) return it->second;
            const std::uint32_t id = static_cast<std::uint32_t>(codes.size());
            codes.push_back(c);
            m.emplace(c, id);
            return id;
        };
        // Confinement (below) must be judged on a property of the marker, not on which markers the
        // selected rule happened to retain. Counting "blocks that retained it" lets a permissive rule
        // inflate its own block counts and filter itself out, while a strict rule keeps markers that
        // are just as contaminated but invisible to the count. The property that actually matters is
        // where a marker's count VARIES with the genotype: a block in which every allele carries it at
        // the same multiplicity contributes a constant offset to every candidate pair alike, which is
        // harmless, whereas a block where it varies is a second signal mixed into the same count.
        // Which blocks each marker varies in, not merely how many. The list is what says where a
        // bubble's markers escape to, and a bubble that loses all of them to one distant block is a
        // different problem from one that loses them to many.
        std::unordered_map<std::uint64_t, std::vector<std::uint32_t>> vary_where;
        std::unordered_map<std::uint64_t, std::uint32_t> vary_blocks;
        // Blocks a marker OCCURS in at all, varying or not. Anchors are constant-multiplicity by
        // construction, so they never appear in vary_blocks and the vary-based confinement test can
        // never fire for them -- an anchor sitting in duplicated sequence keeps the read counts
        // contributed by every copy, and since a block may have only a few hundred anchors, a
        // duplication-biased subset drags the depth median with it.
        std::unordered_map<std::uint64_t, std::uint32_t> occ_blocks;
        std::unordered_map<std::uint64_t, std::uint32_t> evary_blocks;
        for (std::size_t bi = 0; bi < chain.size(); ++bi) {
            const std::vector<AlleleInventory>& inv = kept_inv[bi];
            out_panel->by_block[bi].resize(inv.size());
            // A marker carried by every allele at one and the same multiplicity separates nothing,
            // which is what makes it a clean depth anchor. Restrict to multiplicity 1 so the count
            // reads directly as per-copy coverage.
            std::unordered_map<std::uint64_t, std::uint32_t> seen;
            std::unordered_map<std::uint64_t, std::uint32_t> mult_of;
            for (const AlleleInventory& a : inv) {
                for (const auto& [c, mult] : a.nodes) {
                    ++seen[c];
                    auto m = mult_of.find(c);
                    if (m == mult_of.end()) mult_of.emplace(c, mult);
                    else if (m->second != mult) m->second = 0xffffffffu;   // not constant
                }
            }
            // Only informative markers reach the panel. A marker carried by every allele at the same
            // multiplicity contributes an identical term to every candidate pair, so it cancels in
            // the comparison -- while making any "fraction of my markers that were detected" score
            // degenerate, since it is detected for all alleles alike. Depth comes from the anchors
            // below instead, which is what those constant markers are actually good for.
            std::unordered_map<std::uint64_t, std::uint32_t> eseen;
            std::unordered_map<std::uint64_t, std::uint32_t> emult;
            for (const AlleleInventory& a : inv) {
                for (const auto& [c, mult] : a.edges) {
                    ++eseen[c];
                    auto m = emult.find(c);
                    if (m == emult.end()) emult.emplace(c, mult);
                    else if (m->second != mult) m->second = 0xffffffffu;
                }
            }
            // Rule-independent: where does this marker's count depend on the genotype at all?
            for (const auto& [c, n] : seen) {
                (void)n;
                ++occ_blocks[c];
            }
            for (const auto& [c, n] : seen) {
                if (n < inv.size() || mult_of[c] == 0xffffffffu) {
                    ++vary_blocks[c];
                    vary_where[c].push_back(static_cast<std::uint32_t>(bi));
                }
            }
            for (const auto& [c, n] : eseen) {
                if (n < inv.size() || emult[c] == 0xffffffffu) ++evary_blocks[c];
            }
            const bool presence_only =
                options.rule == MarkerRule::Mixed && chain[bi].kind != BlockKind::Bubble;
            auto informative_node = [&](std::uint64_t c) {
                if (presence_only) return seen[c] < inv.size();   // presence varies; ignore counts
                if (options.rule == MarkerRule::PanGenie) {
                    // unique to a single allele, and occurring exactly once within it
                    return seen[c] == 1 && mult_of[c] == 1;
                }
                if (options.rule == MarkerRule::Unique) {
                    return seen[c] == 1;          // unique to one allele, any copy number
                }
                return seen[c] < inv.size() || mult_of[c] == 0xffffffffu;
            };
            auto informative_edge = [&](std::uint64_t c) {
                return eseen[c] < inv.size() || emult[c] == 0xffffffffu;
            };
            for (std::size_t ai = 0; ai < inv.size(); ++ai) {
                for (const auto& [c, mult] : inv[ai].nodes) {
                    if (!informative_node(c)) continue;
                    // PanGenie stores presence, not count: an allele either has the k-mer or not.
                    const std::uint32_t stored =
                        (options.rule == MarkerRule::PanGenie || presence_only)
                            ? std::min<std::uint32_t>(mult, 1) : mult;
                    out_panel->by_block[bi][ai].nodes.emplace_back(
                        slot_for(node_slot, out_panel->node_codes, c), stored);
                }
                for (const auto& [c, mult] : inv[ai].edges) {
                    if (!informative_edge(c)) continue;
                    out_panel->by_block[bi][ai].edges.emplace_back(
                        slot_for(edge_slot, out_panel->edge_keys, c), mult);
                }
            }
            // A bypass allele spells nothing, so it carries no marker and would make "carried by every
            // allele" unsatisfiable -- the block would lose all its depth anchors. Anchors are about
            // coverage of the sequence that IS there, so measure them over the traversing alleles.
            const std::size_t n_real =
                inv.size() - (blocks[bi].bypass_allele >= 0 && inv.size() > 0 ? 1u : 0u);
            for (const auto& [c, n] : seen) {
                if (n == n_real && n_real > 0 && mult_of[c] == 1) {
                    out_panel->anchor_slots[bi].push_back(slot_for(node_slot, out_panel->node_codes, c));
                }
            }
        }

        // Region uniqueness. A syncmer occurring in more than one place accumulates read counts from
        // all of them, so its observed count no longer reflects the block it is being used to
        // genotype -- measured on cyp2d6, "private" markers of an absent allele carried counts of
        // 38-294 against a depth of 24, and the model duly preferred that absent allele.
        //
        // The comparison used for bubbles (panel-wide count vs what the bubble accounts for) does NOT
        // work here: blocks tile the haplotype, so the accounted-for total equals the panel-wide total
        // for every marker and nothing is ever dropped. What matters instead is CONFINEMENT -- a
        // marker is usable only if all of its occurrences fall inside the one block it is scoring.
        if (graph_for_region_uniqueness != nullptr && options.require_region_unique) {
            const Graph& graph = *graph_for_region_uniqueness;
            std::unordered_map<std::uint32_t, std::uint64_t> expected;    // node slot -> owed
            std::unordered_map<std::uint32_t, std::uint64_t> eexpected;   // edge slot -> owed
            for (std::size_t bi = 0; bi < chain.size(); ++bi) {
                for (std::size_t ai = 0; ai < out_panel->by_block[bi].size(); ++ai) {
                    const std::uint64_t haps = ai < blocks[bi].allele_haplotypes.size()
                                                   ? blocks[bi].allele_haplotypes[ai] : 0;
                    for (const auto& [slot, mult] : out_panel->by_block[bi][ai].nodes) {
                        expected[slot] += haps * mult;
                    }
                    for (const auto& [slot, mult] : out_panel->by_block[bi][ai].edges) {
                        eexpected[slot] += haps * mult;
                    }
                }
            }
            // How many distinct blocks carry each marker, and its total multiplicity within them.
            std::vector<std::uint32_t> blocks_with(out_panel->node_codes.size(), 0);
            std::vector<std::uint32_t> eblocks_with(out_panel->edge_keys.size(), 0);
            out_panel->vary_nodes = vary_blocks.size();
            out_panel->vary_edges = evary_blocks.size();
            for (const auto& [c, n] : vary_blocks) { (void)c; if (n <= 1) ++out_panel->confined_vary_nodes; }
            for (const auto& [c, n] : evary_blocks) { (void)c; if (n <= 1) ++out_panel->confined_vary_edges; }
            // Block-to-block overlap of informative markers, the shape of the confinement loss.
            out_panel->block_overlap.assign(chain.size(), {});
            {
                std::vector<std::vector<std::uint32_t>> tally(chain.size(),
                                                              std::vector<std::uint32_t>(chain.size(), 0));
                for (const auto& [c, where] : vary_where) {
                    (void)c;
                    if (where.size() < 2) continue;
                    for (const std::uint32_t x : where)
                        for (const std::uint32_t y : where)
                            if (x != y) ++tally[x][y];
                }
                for (std::size_t b = 0; b < chain.size(); ++b) {
                    for (std::size_t o = 0; o < chain.size(); ++o) {
                        if (tally[b][o] > 0) out_panel->block_overlap[b].emplace_back(o, tally[b][o]);
                    }
                    std::sort(out_panel->block_overlap[b].begin(), out_panel->block_overlap[b].end(),
                              [](const auto& x, const auto& y) { return x.second > y.second; });
                    if (out_panel->block_overlap[b].size() > 3) out_panel->block_overlap[b].resize(3);
                }
            }
            for (std::size_t sl = 0; sl < out_panel->node_codes.size(); ++sl) {
                const auto it = vary_blocks.find(out_panel->node_codes[sl]);
                blocks_with[sl] = it == vary_blocks.end() ? 0 : it->second;
            }
            for (std::size_t sl = 0; sl < out_panel->edge_keys.size(); ++sl) {
                const auto it = evary_blocks.find(out_panel->edge_keys[sl]);
                eblocks_with[sl] = it == evary_blocks.end() ? 0 : it->second;
            }
            std::vector<std::uint64_t> actual(out_panel->node_codes.size(), 0);
            std::vector<std::uint64_t> eactual(out_panel->edge_keys.size(), 0);
            std::mutex mu;
            run_parallel(graph.paths.size(), options.threads, [&](std::size_t pi) {
                const std::string seq = spell_path_steps_sequence(graph, graph.paths[pi].steps);
                const std::vector<KmerOccurrence> sy = options.all_kmers
                    ? collect_canonical_kmer_occurrences(seq, k) : collect_syncmers(seq, k, s);
                std::unordered_map<std::uint32_t, std::uint32_t> local;
                std::unordered_map<std::uint32_t, std::uint32_t> elocal;
                for (const KmerOccurrence& o : sy) {
                    const auto it = node_slot.find(o.code);
                    if (it != node_slot.end()) ++local[it->second];
                }
                for (std::size_t i = 1; i < sy.size(); ++i) {
                    const auto it = edge_slot.find(
                        adjacency_key(sy[i - 1].code, sy[i].code, sy[i].start - sy[i - 1].start));
                    if (it != edge_slot.end()) ++elocal[it->second];
                }
                std::lock_guard<std::mutex> lock(mu);
                for (const auto& [slot, n] : local) actual[slot] += n;
                for (const auto& [slot, n] : elocal) eactual[slot] += n;
            });
            // Diagnostic: when a marker is dropped for appearing in several blocks, are those blocks
            // ADJACENT (which would mean the shared boundary node between consecutive blocks is the
            // cause -- an artifact of how the chain is cut) or DISTANT (a genuine duplication
            // elsewhere in the region)? The answer decides whether the filter is correct.
            {
                std::vector<std::uint32_t> first_blk(out_panel->node_codes.size(), 0xffffffffu);
                std::vector<std::uint32_t> last_blk(out_panel->node_codes.size(), 0xffffffffu);
                for (std::size_t bi = 0; bi < chain.size(); ++bi) {
                    for (const auto& mset : out_panel->by_block[bi]) {
                        for (const auto& [slot, mult] : mset.nodes) {
                            (void)mult;
                            if (first_blk[slot] == 0xffffffffu) first_blk[slot] = static_cast<std::uint32_t>(bi);
                            last_blk[slot] = static_cast<std::uint32_t>(bi);
                        }
                    }
                }
                std::size_t adj = 0;
                std::size_t distant = 0;
                for (std::size_t sl = 0; sl < out_panel->node_codes.size(); ++sl) {
                    if (blocks_with[sl] <= 1) continue;
                    if (last_blk[sl] - first_blk[sl] <= 1) ++adj; else ++distant;
                }
                out_panel->dropped_adjacent_blocks = adj;
                out_panel->dropped_distant_blocks = distant;
            }
            std::size_t dropped = 0;
            for (std::size_t bi = 0; bi < chain.size(); ++bi) {
                for (const auto& mset : out_panel->by_block[bi]) {
                    out_panel->informative_before_filter += mset.nodes.size();
                    for (const auto& [slot, mult] : mset.nodes) {
                        (void)mult;
                        if (blocks_with[slot] > 1) ++out_panel->dropped_multi_block;
                        else if (actual[slot] > expected[slot]) ++out_panel->dropped_over_expected;
                    }
                }
            }
            // Snapshot before the erase, so a surviving marker can be interrogated afterwards.
            out_panel->dbg_vary.assign(out_panel->node_codes.size(), 0);
            out_panel->dbg_occ.assign(out_panel->node_codes.size(), 0);
            out_panel->dbg_actual.assign(out_panel->node_codes.size(), 0);
            out_panel->dbg_expected.assign(out_panel->node_codes.size(), 0);
            for (std::size_t sl = 0; sl < out_panel->node_codes.size(); ++sl) {
                out_panel->dbg_vary[sl] = blocks_with[sl];
                const auto io = occ_blocks.find(out_panel->node_codes[sl]);
                out_panel->dbg_occ[sl] = io == occ_blocks.end() ? 0 : io->second;
                out_panel->dbg_actual[sl] = actual[sl];
                const auto ie = expected.find(static_cast<std::uint32_t>(sl));
                out_panel->dbg_expected[sl] = ie == expected.end() ? 0 : ie->second;
            }
            std::vector<std::vector<std::vector<std::pair<std::uint32_t, std::uint32_t>>>> pre_filter;
            if (options.restore_stripped_alleles) {
                pre_filter.resize(chain.size());
                for (std::size_t bi = 0; bi < chain.size(); ++bi) {
                    pre_filter[bi].resize(out_panel->by_block[bi].size());
                    for (std::size_t ai = 0; ai < out_panel->by_block[bi].size(); ++ai) {
                        pre_filter[bi][ai] = out_panel->by_block[bi][ai].nodes;
                    }
                }
            }
            for (std::size_t bi = 0; bi < chain.size(); ++bi) {
                for (auto& mset : out_panel->by_block[bi]) {
                    const std::size_t before = mset.nodes.size();
                    mset.nodes.erase(
                        std::remove_if(mset.nodes.begin(), mset.nodes.end(),
                                       [&](const std::pair<std::uint32_t, std::uint32_t>& e) {
                                           // Ambiguous across blocks, or seen more often in the panel
                                           // than the blocks account for (sequence the tiling missed,
                                           // e.g. a repeat traversed more than once).
                                           return blocks_with[e.first] > 1 ||
                                                  actual[e.first] > expected[e.first];
                                       }),
                        mset.nodes.end());
                    dropped += before - mset.nodes.size();
                    const std::size_t ebefore = mset.edges.size();
                    mset.edges.erase(
                        std::remove_if(mset.edges.begin(), mset.edges.end(),
                                       [&](const std::pair<std::uint32_t, std::uint32_t>& e) {
                                           return eblocks_with[e.first] > 1 ||
                                                  eactual[e.first] > eexpected[e.first];
                                       }),
                        mset.edges.end());
                    dropped += ebefore - mset.edges.size();
                }
                // Confinement can strip an allele to nothing. When it does, the block loses the ability
                // to decide genotypes by COMPOSITION and falls back on the absolute count level of
                // whatever alleles still have markers -- the one scale-dependent decision this model
                // makes, and the one measured to be unreliable (a true heterozygote reads 1.78 x lambda
                // and a true homozygote 2.33 on the paralogous fixture, against 1.0 and 2.0).
                //
                // The trade this option offers: put the stripped allele's markers back, contamination
                // and all, and let the RATIO between alleles carry the decision. A paralogue inflates
                // both alleles' counts by a similar amount, so it corrupts absolute counts far more
                // than it corrupts their ratio -- which is why every SNV caller decides zygosity from
                // allele balance rather than from depth.
                if (options.restore_stripped_alleles) {
                    bool stripped = false;
                    for (std::size_t ai = 0; ai < out_panel->by_block[bi].size(); ++ai) {
                        if (blocks[bi].bypass_allele >= 0 &&
                            ai == static_cast<std::size_t>(blocks[bi].bypass_allele)) continue;
                        if (out_panel->by_block[bi][ai].nodes.empty() &&
                            !pre_filter[bi][ai].empty()) { stripped = true; break; }
                    }
                    if (stripped) {
                        for (std::size_t ai = 0; ai < out_panel->by_block[bi].size(); ++ai) {
                            out_panel->by_block[bi][ai].nodes = pre_filter[bi][ai];
                        }
                        ++out_panel->blocks_restored;
                    }
                }
                // Anchors need the OCCURRENCE-based test, not the vary-based one, and an expectation
                // of their own: they are absent from `by_block`, so `expected` is zero for them and the
                // old bound degenerated to "seen more than once per haplotype on average". An anchor is
                // carried by every allele of its block at multiplicity 1, so the panel should show it
                // exactly once per traversing haplotype.
                // Expect an anchor once per TRAVERSING haplotype; bypassers contribute no reads here.
                const std::uint64_t anchor_expected =
                    blocks[bi].n_traversing > 0 ? blocks[bi].n_traversing : blocks[bi].allele_of.size();
                auto& anch = out_panel->anchor_slots[bi];
                anch.erase(std::remove_if(anch.begin(), anch.end(),
                                          [&](std::uint32_t slot) {
                                              const auto it = occ_blocks.find(out_panel->node_codes[slot]);
                                              const std::uint32_t nb = it == occ_blocks.end() ? 0 : it->second;
                                              return nb > 1 || actual[slot] > anchor_expected;
                                          }),
                           anch.end());
            }
            out_panel->region_filtered_markers = dropped;
        }
    }
    return out;
}

void write_block_marker_audit(
    const std::string& out_prefix,
    const std::vector<BlockMarkerStats>& stats) {

    const std::string path = out_prefix + ".audit.blockmarkers.tsv";
    std::ofstream f(path);
    if (!f) throw std::runtime_error("genotype: cannot write " + path);
    f << "block_index\tblock_kind\tbubble_id\tn_alleles\tn_haplotypes\tblock_bp"
         "\tn_informative_nodes\tn_informative_edges\tmedian_sep_nodes\tmedian_sep_edges"
         "\tunseparable_alleles\tseparable_mass_nodes\tseparable_mass_edges"
         "\tmarkers_over_cap\tmax_multiplicity_seen\tseconds\n";
    for (const BlockMarkerStats& s : stats) {
        f << s.block_index << '\t' << (s.is_bubble ? "bubble" : "backbone") << '\t' << s.bubble_id
          << '\t' << s.n_alleles << '\t' << s.n_haplotypes << '\t' << s.block_bp << '\t'
          << s.n_informative_nodes << '\t' << s.n_informative_edges << '\t' << s.median_sep_nodes
          << '\t' << s.median_sep_edges << '\t' << s.unseparable_alleles << '\t'
          << s.separable_mass_nodes << '\t' << s.separable_mass_edges << '\t'
          << s.markers_over_cap << '\t' << s.max_marker_multiplicity_seen << '\t' << s.seconds << '\n';
    }
}

} // namespace panvar
