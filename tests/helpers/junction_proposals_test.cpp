// M9-G gates for the junction-aware proposals (JP1.1, junction_proposals.hpp) against an independent oracle: whole
// haplotypes are spelled and every inventory K-mer is counted in them and in the reads by direct string scanning, and the
// Poisson kernel is summed over the whole inventory. Score DIFFERENCES between candidates must agree (the fast path drops
// a candidate-independent constant).
#include "panvar/junction_proposals.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include <random>
#include <set>
#include <string>
#include <thread>
#include <vector>

using namespace panvar;

namespace {

int failures = 0;
void check(bool ok, const std::string& what) {
    std::cout << (ok ? "  PASS " : "  FAIL ") << what << '\n';
    if (!ok) ++failures;
}

std::mt19937 rng(20261020);
std::string random_seq(std::size_t n) {
    static const char* b = "ACGT";
    std::string s(n, 'A');
    for (char& c : s) c = b[rng() % 4];
    return s;
}
std::string revcomp(const std::string& s) {
    std::string r(s.rbegin(), s.rend());
    for (char& c : r) c = c == 'A' ? 'T' : c == 'C' ? 'G' : c == 'G' ? 'C' : 'A';
    return r;
}
std::string canon_str(const std::string& w) {
    const std::string r = revcomp(w);
    return w <= r ? w : r;
}

// reads tiling the given haplotypes
std::vector<ReadPair> reads_from(const std::vector<std::string>& haps, std::size_t len = 60, std::size_t step = 10) {
    std::vector<ReadPair> out;
    for (const std::string& h : haps) {
        for (std::size_t i = 0; i + len <= h.size(); i += step) {
            ReadPair p;
            p.first.bases = h.substr(i, len);
            p.first.qualities = std::string(len, 'I');
            p.second.bases = revcomp(h.substr(i, len));
            p.second.qualities = std::string(len, 'I');
            out.push_back(p);
        }
    }
    return out;
}

struct Fixture {
    std::vector<std::string> alleles;
    JunctionContext ctx;
    std::vector<ReadPair> reads;
};

// The oracle: the full Poisson score of candidate (x on hap 1, y on hap 2) over the inventory of the comparison.
double oracle(const Fixture& f, std::size_t x, std::size_t y, double rate, double bg) {
    const std::size_t k = kJunctionK;
    std::set<std::string> inventory;
    for (int h = 0; h < 2; ++h) {
        const std::string& pre = f.ctx.prefix[h];
        const std::string& suf = f.ctx.suffix[h];
        const std::string left = pre.substr(pre.size() > k - 1 ? pre.size() - (k - 1) : 0), right = suf.substr(0, std::min(k - 1, suf.size()));
        for (const std::string& a : f.alleles) {
            const std::string w = left + a + right;
            for (std::size_t i = 0; i + k <= w.size(); ++i) inventory.insert(canon_str(w.substr(i, k)));
        }
    }
    std::map<std::string, double> copies, observed;
    const std::string hap[2] = {f.ctx.prefix[0] + f.alleles[x] + f.ctx.suffix[0], f.ctx.prefix[1] + f.alleles[y] + f.ctx.suffix[1]};
    for (const std::string& h : hap) {
        for (std::size_t i = 0; i + k <= h.size(); ++i) {
            const std::string c = canon_str(h.substr(i, k));
            if (inventory.count(c)) copies[c] += 1;
        }
    }
    for (const ReadPair& p : f.reads) {
        for (const std::string* s : {&p.first.bases, &p.second.bases}) {
            for (std::size_t i = 0; i + k <= s->size(); ++i) {
                const std::string c = canon_str(s->substr(i, k));
                if (inventory.count(c)) observed[c] += 1;
            }
        }
    }
    double total = 0.0;
    for (const std::string& m : inventory) {
        const double lambda = rate * copies[m] + bg;
        total += observed[m] * std::log(lambda) - lambda;
    }
    return total;
}

Fixture make(std::vector<std::string> alleles, std::size_t inc1, std::size_t inc2, std::size_t truth1, std::size_t truth2,
             std::size_t flank = 120) {
    Fixture f;
    f.alleles = std::move(alleles);
    const std::string left = random_seq(flank), right = random_seq(flank);
    for (int h = 0; h < 2; ++h) {
        f.ctx.prefix[h] = left;
        f.ctx.suffix[h] = right;
    }
    f.ctx.incumbent[0] = inc1;
    f.ctx.incumbent[1] = inc2;
    f.reads = reads_from({left + f.alleles[truth1] + right, left + f.alleles[truth2] + right});
    return f;
}

const double kRate = 1.6, kBg = 0.012;

std::map<std::pair<std::size_t, std::size_t>, JunctionPairScore> by_pair(const std::vector<JunctionPairScore>& rows) {
    std::map<std::pair<std::size_t, std::size_t>, JunctionPairScore> m;
    for (const auto& r : rows) m[{r.x, r.y}] = r;
    return m;
}

double oracle_agreement(const Fixture& f) {
    const JunctionReadKmers t = JunctionReadKmers::from_fragments(f.reads);
    const auto rows = junction_rank_all(f.alleles, f.ctx, t, kRate, kBg);
    const auto m = by_pair(rows);
    const auto& base = m.at({0, 0});
    const double o0 = oracle(f, 0, 0, kRate, kBg);
    double worst = 0.0;
    for (const auto& [p, r] : m) {
        worst = std::max(worst, std::fabs((r.ab - base.ab) - (oracle(f, p.first, p.second, kRate, kBg) - o0)));
        if (p.first != p.second) worst = std::max(worst, std::fabs((r.ba - base.ab) - (oracle(f, p.second, p.first, kRate, kBg) - o0)));
    }
    return worst;
}

std::set<std::pair<std::string, std::string>> content_set(const std::vector<std::string>& alleles, const JunctionProposal& p) {
    std::set<std::pair<std::string, std::string>> s;
    for (const auto& [x, y] : p.pairs) s.insert(std::minmax(alleles[x], alleles[y]));
    return s;
}

}  // namespace

int main() {
    std::cout << "junction proposals (JP1.1):\n";
    // F1 short alleles only distinguishable through junctions (5 bp: no internal K-mer)
    {
        std::vector<std::string> al;
        for (int i = 0; i < 6; ++i) al.push_back(random_seq(5));
        const Fixture f = make(al, 0, 1, 3, 4);
        const auto rows = junction_rank_all(f.alleles, f.ctx, JunctionReadKmers::from_fragments(f.reads), kRate, kBg);
        check(rows[0].x == 3 && rows[0].y == 4 && rows[1].tie_class != rows[0].tie_class, "F1 a 5-bp allele pair ranks first through junctions, untied");
        check(oracle_agreement(f) < 1e-9, "F1 fast delta equals the whole-haplotype oracle");
    }
    // F2 empty allele
    {
        const Fixture f = make({"", random_seq(4), random_seq(9), random_seq(12)}, 1, 3, 0, 2);
        const auto rows = junction_rank_all(f.alleles, f.ctx, JunctionReadKmers::from_fragments(f.reads), kRate, kBg);
        check(rows[0].x == 0 && rows[0].y == 2, "F2 the empty allele is scored and ranks with the truth");
        check(oracle_agreement(f) < 1e-9, "F2 oracle");
    }
    // F3 K-mers spanning both boundaries of a very short allele, multiplicity kept
    {
        const Fixture f = make({"ACG", "TTA", "ACGACG"}, 0, 1, 0, 1);
        check(oracle_agreement(f) < 1e-9, "F3 both-boundary K-mers and multiplicities agree with the oracle");
    }
    // F4 repeated junction K-mer elsewhere in the locus (elsewhere multiplicity), and its refresh when elsewhere changes
    // under unchanged flanks
    {
        Fixture f = make({"A", "C", "G"}, 0, 1, 0, 1, 120);
        const std::string junction = f.ctx.prefix[0].substr(f.ctx.prefix[0].size() - (kJunctionK - 1)) + "A";
        for (int h = 0; h < 2; ++h) f.ctx.suffix[h] = f.ctx.suffix[h].substr(0, 60) + junction + junction + f.ctx.suffix[h].substr(60);
        check(oracle_agreement(f) < 1e-9, "F4 junction K-mer repeated elsewhere: oracle agreement");
        const JunctionReadKmers t = JunctionReadKmers::from_fragments(f.reads);
        const auto a = by_pair(junction_rank_all(f.alleles, f.ctx, t, kRate, kBg));
        Fixture g = f;   // same flanks (first/last 30 bases unchanged), elsewhere copies removed far away
        for (int h = 0; h < 2; ++h) g.ctx.suffix[h] = f.ctx.suffix[h].substr(0, 60) + random_seq(2 * kJunctionK) + f.ctx.suffix[h].substr(60 + 2 * kJunctionK);
        const auto b = by_pair(junction_rank_all(g.alleles, g.ctx, t, kRate, kBg));
        check(g.ctx.suffix[0].substr(0, 30) == f.ctx.suffix[0].substr(0, 30) && std::fabs(a.at({0, 0}).ab - b.at({0, 0}).ab) > 1e-6,
              "F4 same flanks, changed elsewhere multiplicity: the score changes");
        check(oracle_agreement(g) < 1e-9, "F4 refreshed ranking agrees with the oracle");
        // the complete cache key distinguishes the two incumbents (their whole label paths differ elsewhere)
        JunctionProposalCache cache(4);
        const auto k1 = JunctionProposalCache::key(1, {0, 0, 5}, {0, 1, 5});
        const auto k2 = JunctionProposalCache::key(1, {0, 0, 7}, {0, 1, 5});   // a far-block difference, same flanks
        cache.insert(k1, junction_local_pairs(f.alleles, f.ctx, t, kRate, kBg, 2));
        check(cache.find(k2) == nullptr && cache.find(k1) != nullptr, "F4 cache: a different elsewhere context is a miss, the same is a hit");
        check(JunctionProposalCache::key(1, {0, 1, 5}, {0, 0, 5}) == k1, "cache key is the unordered incumbent (homologue swap = same key)");
    }
    // F5 reverse-complement and F6 homologue-swap invariance, F7 allele-index permutation with sequences preserved
    {
        std::vector<std::string> al;
        for (int i = 0; i < 6; ++i) al.push_back(random_seq(static_cast<std::size_t>(rng() % 40)));
        const Fixture f = make(al, 0, 2, 1, 4);
        const JunctionReadKmers t = JunctionReadKmers::from_fragments(f.reads);
        const auto a = by_pair(junction_rank_all(f.alleles, f.ctx, t, kRate, kBg));
        Fixture r = f;
        for (auto& s : r.alleles) s = revcomp(s);
        for (int h = 0; h < 2; ++h) {
            r.ctx.prefix[h] = revcomp(f.ctx.suffix[h]);
            r.ctx.suffix[h] = revcomp(f.ctx.prefix[h]);
        }
        const auto b = by_pair(junction_rank_all(r.alleles, r.ctx, t, kRate, kBg));   // reads unchanged: canonical K-mers
        double worst = 0.0;
        for (const auto& [p, v] : a) worst = std::max(worst, std::fabs(v.best - b.at(p).best));
        check(worst < 1e-9, "F5 reverse-complement invariance");
        Fixture s = f;
        std::swap(s.ctx.prefix[0], s.ctx.prefix[1]);
        std::swap(s.ctx.suffix[0], s.ctx.suffix[1]);
        std::swap(s.ctx.incumbent[0], s.ctx.incumbent[1]);
        const auto c = by_pair(junction_rank_all(s.alleles, s.ctx, t, kRate, kBg));
        worst = 0.0;
        for (const auto& [p, v] : a) worst = std::max(worst, std::fabs(v.best - c.at(p).best));
        check(worst < 1e-9, "F6 homologue-swap invariance");
        const std::vector<std::size_t> perm = {3, 0, 5, 1, 4, 2};   // new index i holds old allele perm[i]
        Fixture q = f;
        for (std::size_t i = 0; i < perm.size(); ++i) q.alleles[i] = f.alleles[perm[i]];
        for (int h = 0; h < 2; ++h) q.ctx.incumbent[h] = static_cast<std::size_t>(std::find(perm.begin(), perm.end(), f.ctx.incumbent[h]) - perm.begin());
        const auto d = by_pair(junction_rank_all(q.alleles, q.ctx, t, kRate, kBg));
        worst = 0.0;
        for (const auto& [p, v] : d) worst = std::max(worst, std::fabs(v.best - a.at(std::minmax(perm[p.first], perm[p.second])).best));
        check(worst < 1e-9, "F7 allele-index permutation (sequences preserved)");
    }
    // F8 identical sequence under two labels: identical scores, one tie class
    {
        std::vector<std::string> al = {random_seq(10), random_seq(12), random_seq(8)};
        al.push_back(al[1]);
        const Fixture f = make(al, 0, 2, 1, 2);
        const auto m = by_pair(junction_rank_all(f.alleles, f.ctx, JunctionReadKmers::from_fragments(f.reads), kRate, kBg));
        check(std::fabs(m.at({1, 2}).best - m.at({2, 3}).best) < 1e-12 && m.at({1, 2}).tie_class == m.at({2, 3}).tie_class,
              "F8 identical sequences under different labels tie");
    }
    // F9 genuinely uninformative: equal-length alleles, no read touching the block
    {
        std::vector<std::string> al;
        for (int i = 0; i < 5; ++i) al.push_back(random_seq(20));
        Fixture f = make(al, 0, 2, 1, 2);
        f.reads = reads_from({random_seq(400)});
        const auto rows = junction_rank_all(f.alleles, f.ctx, JunctionReadKmers::from_fragments(f.reads), kRate, kBg);
        check(rows.front().tie_class == rows.back().tie_class, "F9 uninformative block stays one tie class");
    }
    // T1 a tie crossing the cutoff: content-based selection is invariant under allele-index permutation, repeats and threads
    {
        std::vector<std::string> al;
        for (int i = 0; i < 12; ++i) al.push_back(random_seq(20));
        Fixture f = make(al, 0, 1, 2, 3);
        f.reads = reads_from({random_seq(300)});   // no read touches the block: one tie class over all 78 pairs
        const JunctionReadKmers t = JunctionReadKmers::from_fragments(f.reads);
        const JunctionProposal p = junction_local_pairs(f.alleles, f.ctx, t, kRate, kBg, 32);
        check(p.cutoff_tie_crosses && p.pairs.size() == 32, "T1 the tie class crosses rank 32 and the shortlist is bounded at 32");
        const auto base = content_set(f.alleles, p);
        bool same = true;
        std::vector<std::size_t> perm(al.size());
        for (int trial = 0; trial < 5; ++trial) {
            for (std::size_t i = 0; i < perm.size(); ++i) perm[i] = i;
            std::shuffle(perm.begin(), perm.end(), rng);
            Fixture q = f;
            for (std::size_t i = 0; i < perm.size(); ++i) q.alleles[i] = f.alleles[perm[i]];
            for (int h = 0; h < 2; ++h) q.ctx.incumbent[h] = static_cast<std::size_t>(std::find(perm.begin(), perm.end(), f.ctx.incumbent[h]) - perm.begin());
            same = same && content_set(q.alleles, junction_local_pairs(q.alleles, q.ctx, t, kRate, kBg, 32)) == base;
        }
        check(same, "T1 the selected sequence-content set is invariant under 5 allele-index permutations");
        std::vector<std::set<std::pair<std::string, std::string>>> got(4);
        std::vector<std::thread> threads;
        for (int i = 0; i < 4; ++i) threads.emplace_back([&, i] { got[i] = content_set(f.alleles, junction_local_pairs(f.alleles, f.ctx, t, kRate, kBg, 32)); });
        for (auto& th : threads) th.join();
        check(std::all_of(got.begin(), got.end(), [&](const auto& g) { return g == base; }), "T1 identical selection across 4 concurrent invocations");
    }
    // C1 neighbouring-block move refreshes the context: a different left flank gives a different ranking and a different key
    {
        std::vector<std::string> al;
        for (int i = 0; i < 6; ++i) al.push_back(random_seq(6));
        const Fixture f = make(al, 0, 1, 2, 3);
        Fixture g = f;
        for (int h = 0; h < 2; ++h) g.ctx.prefix[h] = f.ctx.prefix[h].substr(0, f.ctx.prefix[h].size() - 10) + random_seq(10);
        const JunctionReadKmers t = JunctionReadKmers::from_fragments(f.reads);
        const auto a = by_pair(junction_rank_all(f.alleles, f.ctx, t, kRate, kBg));
        const auto b = by_pair(junction_rank_all(g.alleles, g.ctx, t, kRate, kBg));
        check(std::fabs(a.at({2, 3}).best - b.at({2, 3}).best) > 1e-6 && oracle_agreement(g) < 1e-9,
              "C1 a neighbouring-block change refreshes the ranking (oracle on the new context)");
        check(JunctionProposalCache::key(3, {1, 0, 2}, {1, 1, 2}) != JunctionProposalCache::key(3, {4, 0, 2}, {1, 1, 2}),
              "C1 the cache key changes with a neighbouring-block allele");
    }
    // C2 eligibility: all pair emissions tied AND more pairs than the limit
    {
        check(junction_block_eligible(std::vector<double>(40, 0.0), 32), "C2 tied block with 40 pairs is eligible");
        check(!junction_block_eligible(std::vector<double>(32, 0.0), 32), "C2 tied block with 32 pairs is not (every pair already generated)");
        std::vector<double> e(40, 0.0);
        e[7] = -1.0;
        check(!junction_block_eligible(e, 32), "C2 one differing emission: not eligible");
    }
    // C3 bounded cache eviction
    {
        JunctionProposalCache cache(2);
        JunctionProposal p;
        cache.insert(JunctionProposalCache::key(0, {1}, {1}), p);
        cache.insert(JunctionProposalCache::key(0, {2}, {2}), p);
        cache.insert(JunctionProposalCache::key(0, {3}, {3}), p);
        check(cache.size() == 2 && cache.find(JunctionProposalCache::key(0, {1}, {1})) == nullptr, "C3 the cache is bounded (oldest evicted)");
    }
    std::cout << (failures == 0 ? "junction_proposals: all passed\n" : "junction_proposals: FAILURES\n");
    return failures == 0 ? 0 : 1;
}
