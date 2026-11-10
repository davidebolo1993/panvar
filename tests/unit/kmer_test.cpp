// Unit tests of the shared k-mer primitives (include/panvar/kmer.hpp).
//   kmer_test [graph.gfa[.gz] ...]
// Each rolling primitive is checked against a brute-force version on random and low-complexity sequences,
// on both strands; with graphs, against collect_canonical_kmer_occurrences on their spelled paths.

#include "panvar/gfa.hpp"
#include "panvar/graph_utils.hpp"
#include "panvar/kmer.hpp"
#include "panvar/syncmer.hpp"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <random>
#include <set>
#include <string>
#include <tuple>
#include <vector>

using namespace panvar;

namespace {

int failures = 0;

void check(bool ok, const std::string& what) {
    std::printf("  %s %s\n", ok ? "ok  " : "FAIL", what.c_str());
    if (!ok) ++failures;
}

std::string random_sequence(std::mt19937_64& rng, std::size_t n, double n_rate = 0.0, const char* alphabet = "ACGT") {
    std::string s(n, 'A');
    std::uniform_real_distribution<double> u(0.0, 1.0);
    const std::size_t a = std::string(alphabet).size();
    for (auto& c : s) c = u(rng) < n_rate ? 'N' : alphabet[rng() % a];
    return s;
}

std::string revcomp(const std::string& s) {
    std::string r(s.rbegin(), s.rend());
    for (auto& c : r) c = c == 'A' ? 'T' : c == 'C' ? 'G' : c == 'G' ? 'C' : c == 'T' ? 'A' : 'N';
    return r;
}

// Brute force: every window of k ACGT bases, coded directly.
std::vector<std::tuple<std::uint64_t, std::size_t, bool>> brute_kmers(const std::string& s, std::size_t k) {
    std::vector<std::tuple<std::uint64_t, std::size_t, bool>> out;
    for (std::size_t i = 0; i + k <= s.size(); ++i) {
        std::uint64_t code = 0;
        bool clean = true;
        for (std::size_t j = 0; j < k; ++j) {
            const int b = kmer::base_code(s[i + j]);
            if (b < 0) { clean = false; break; }
            code = (code << 2) | static_cast<std::uint64_t>(b);
        }
        if (!clean) continue;
        const std::uint64_t rc = kmer::reverse_complement(code, k);
        out.emplace_back(code <= rc ? code : rc, i, code <= rc);
    }
    return out;
}

std::vector<std::tuple<std::uint64_t, std::size_t, bool>> rolling_kmers(const std::string& s, std::size_t k) {
    std::vector<std::tuple<std::uint64_t, std::size_t, bool>> out;
    kmer::for_each_canonical_kmer(s.data(), s.size(), k, [&](std::uint64_t c, std::size_t p, bool f) { out.emplace_back(c, p, f); });
    return out;
}

std::vector<std::pair<std::uint64_t, std::size_t>> rolling_syncmers(const std::string& s, std::size_t k, std::size_t sl) {
    std::vector<std::pair<std::uint64_t, std::size_t>> out;
    kmer::for_each_closed_syncmer(s.data(), s.size(), k, sl, [&](std::uint64_t c, std::size_t p, bool) { out.emplace_back(c, p); });
    return out;
}

std::vector<std::pair<std::uint64_t, std::size_t>> brute_syncmers(const std::string& s, std::size_t k, std::size_t sl) {
    std::vector<std::pair<std::uint64_t, std::size_t>> out;
    for (const auto& [c, p, f] : brute_kmers(s, k)) {
        (void)f;
        if (kmer::is_closed_syncmer(c, k, sl)) out.emplace_back(c, p);
    }
    return out;
}

double density(const std::string& s, std::size_t k, std::size_t sl) {
    const double kmers = static_cast<double>(brute_kmers(s, k).size());
    return kmers > 0 ? static_cast<double>(rolling_syncmers(s, k, sl).size()) / kmers : 0.0;
}

} // namespace

int main(int argc, char** argv) {
    std::mt19937_64 rng(20261110);
    const std::size_t ks[] = {15, 21, 25, 31};

    // Canonical k-mers: brute force, the existing collector, the strand flag.
    {
        bool same_brute = true, same_existing = true, flag_ok = true;
        for (int t = 0; t < 40; ++t) {
            const std::string s = random_sequence(rng, 50 + rng() % 4000, t % 3 == 0 ? 0.01 : 0.0);
            for (const std::size_t k : ks) {
                const auto r = rolling_kmers(s, k);
                same_brute = same_brute && r == brute_kmers(s, k);
                const auto e = collect_canonical_kmer_occurrences(s, k);
                same_existing = same_existing && e.size() == r.size();
                for (std::size_t i = 0; same_existing && i < e.size(); ++i)
                    same_existing = e[i].code == std::get<0>(r[i]) && e[i].start == std::get<1>(r[i]);
                for (const auto& [c, p, f] : r)
                    if (f != (kmer::decode(c, k) == std::string(s, p, k).substr(0, k))) { flag_ok = false; break; }
            }
        }
        check(same_brute, "canonical k-mers equal brute force (k 15/21/25/31, with and without N)");
        check(same_existing, "canonical k-mers equal collect_canonical_kmer_occurrences");
        check(flag_ok, "the forward flag says which strand the canonical code spells");
        check(kmer::reverse_complement(kmer::reverse_complement(0x1234567ULL, 31), 31) == 0x1234567ULL,
              "reverse complement is an involution");
    }

    // Closed syncmers: brute force, strand symmetry, density, composition.
    {
        bool same = true, symmetric = true;
        for (int t = 0; t < 40; ++t) {
            const std::string s = random_sequence(rng, 50 + rng() % 4000, t % 4 == 0 ? 0.005 : 0.0);
            for (const std::size_t k : ks) {
                const std::size_t sl = k == 15 ? 5 : 11;
                const auto r = rolling_syncmers(s, k, sl);
                same = same && r == brute_syncmers(s, k, sl);
                std::set<std::pair<std::uint64_t, std::size_t>> fw(r.begin(), r.end()), rv;
                for (const auto& [c, p] : rolling_syncmers(revcomp(s), k, sl)) rv.emplace(c, s.size() - k - p);
                symmetric = symmetric && fw == rv;
            }
        }
        check(same, "rolling closed syncmers equal the one-k-mer-at-a-time test");
        check(symmetric, "a sequence and its reverse complement have the same syncmers");
        // Low-complexity input must not break the queue (long runs of equal hashes).
        const std::string lowc = std::string(500, 'A') + std::string(300, 'N') + std::string(400, 'C') +
                                 random_sequence(rng, 600, 0.0, "AT");
        check(rolling_syncmers(lowc, 31, 11) == brute_syncmers(lowc, 31, 11), "homopolymers and N runs: rolling equals brute force");
        const double d = density(random_sequence(rng, 200000), 31, 11), expect = 2.0 / 21.0;
        char buf[160];
        std::snprintf(buf, sizeof buf, "density on random sequence is about 2/(k-s+1) (%.4f vs %.4f)", d, expect);
        check(std::fabs(d - expect) < 0.1 * expect, buf);
        const double at = density(random_sequence(rng, 100000, 0.0, "AAAATTTTCG"), 31, 11);
        const double gc = density(random_sequence(rng, 100000, 0.0, "GGGGCCCCAT"), 31, 11);
        std::snprintf(buf, sizeof buf, "density does not follow base composition (AT-rich %.4f, GC-rich %.4f)", at, gc);
        check(std::fabs(at - gc) < 0.15 * expect, buf);
    }

    // Low complexity.
    {
        auto code_of = [](const std::string& s) { std::uint64_t c = 0; for (char x : s) c = (c << 2) | kmer::base_code(x); return c; };
        const double homo = kmer::dust_score(code_of(std::string(31, 'A')), 31);
        std::string di, tri;
        for (int i = 0; i < 31; ++i) { di += "AT"[i % 2]; tri += "CAG"[i % 3]; }
        double rnd = 0.0;
        for (int t = 0; t < 200; ++t) rnd += kmer::dust_score(code_of(random_sequence(rng, 31)), 31) / 200.0;
        char buf[200];
        std::snprintf(buf, sizeof buf, "DUST: homopolymer 14.5 (%.2f), dinucleotide (%.2f) and trinucleotide (%.2f) repeats high, random low (%.2f)",
                      homo, kmer::dust_score(code_of(di), 31), kmer::dust_score(code_of(tri), 31), rnd);
        check(std::fabs(homo - 14.5) < 1e-9 && kmer::dust_score(code_of(di), 31) > 5 && kmer::dust_score(code_of(tri), 31) > 3 && rnd < 1.0, buf);
    }

    // Real graphs: every spelled path, against the existing collector.
    for (int a = 1; a < argc; ++a) {
        const Graph g = parse_gfa(argv[a]);
        std::size_t paths = 0, kmers = 0;
        bool same = true;
        for (const PathRecord& p : g.paths) {
            const std::string s = spell_path_steps_sequence(g, p.steps);
            const auto r = rolling_kmers(s, 31);
            const auto e = collect_canonical_kmer_occurrences(s, 31);
            same = same && e.size() == r.size();
            for (std::size_t i = 0; same && i < e.size(); ++i) same = e[i].code == std::get<0>(r[i]) && e[i].start == std::get<1>(r[i]);
            ++paths; kmers += r.size();
        }
        check(same, std::string(argv[a]) + ": " + std::to_string(paths) + " paths, " + std::to_string(kmers) +
                    " 31-mers equal collect_canonical_kmer_occurrences");
    }

    std::printf("\nkmer_test: %s\n", failures ? (std::to_string(failures) + " check(s) failed").c_str() : "all checks passed");
    return failures ? 1 : 0;
}
