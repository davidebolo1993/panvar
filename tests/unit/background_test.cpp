// Unit tests of the background store: the Elias-Fano sequence, the super-k-mer splitter, and a segment
// built from a synthetic genome, looked up k-mer by k-mer against a brute-force count.
//   background_test <scratch dir>

#include "panvar/background.hpp"
#include "panvar/elias_fano.hpp"
#include "panvar/kmer.hpp"

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <map>
#include <random>
#include <set>
#include <string>
#include <vector>

using namespace panvar;

namespace {

int failures = 0;

void check(bool ok, const std::string& what) {
    std::printf("  %s %s\n", ok ? "ok  " : "FAIL", what.c_str());
    if (!ok) ++failures;
}

std::string random_sequence(std::mt19937_64& rng, std::size_t n) {
    std::string s(n, 'A');
    for (auto& c : s) c = "ACGT"[rng() % 4];
    return s;
}

} // namespace

int main(int argc, char** argv) {
    if (argc < 2) { std::fprintf(stderr, "usage: background_test <scratch dir>\n"); return 2; }
    const std::filesystem::path dir = argv[1];
    std::filesystem::create_directories(dir);
    std::mt19937_64 rng(31);

    // Elias-Fano: membership and order on random sets of several sizes and universes.
    {
        bool ok = true;
        for (const std::size_t n : {0, 1, 2, 7, 100, 5000}) {
            for (const std::size_t ub : {8, 20, 62}) {
                std::set<std::uint64_t> s;
                const std::uint64_t umax = ub >= 64 ? ~0ULL : (1ULL << ub);
                while (s.size() < std::min<std::uint64_t>(n, umax)) s.insert(rng() % umax);
                const std::vector<std::uint64_t> v(s.begin(), s.end());
                const std::vector<std::uint64_t> w = ef::build(v, ub);
                const ef::View e(w.data());
                ok = ok && e.size() == v.size() && e.words() == w.size();
                for (std::size_t i = 0; i < v.size(); ++i) ok = ok && e.find(v[i]) == static_cast<std::int64_t>(i) && e.at(i) == v[i];
                std::vector<std::uint64_t> seen;
                e.for_each_value([&](std::uint64_t, std::uint64_t x) { seen.push_back(x); });
                ok = ok && seen == v;
                for (int t = 0; t < 2000; ++t) {
                    const std::uint64_t q = rng() % umax;
                    ok = ok && (e.find(q) >= 0) == (s.count(q) > 0);
                }
            }
        }
        check(ok, "Elias-Fano: find, at and iteration equal brute force (sizes 0 to 5000, universes 2^8 to 2^62)");
    }

    // Super-k-mers: they cover every k-mer exactly once, each in the partition partition_of gives it.
    {
        bool ok = true;
        for (int t = 0; t < 20; ++t) {
            std::string s = random_sequence(rng, 200 + rng() % 5000);
            if (t % 3 == 0) for (int i = 0; i < 5; ++i) s[rng() % s.size()] = 'N';
            std::vector<std::pair<std::size_t, std::uint32_t>> got;   // (k-mer start, partition)
            bg::for_each_super_kmer(s.data(), s.size(), 64, [&](std::size_t start, std::size_t len, std::uint32_t p) {
                for (std::size_t i = 0; i + bg::kK <= len; ++i) got.emplace_back(start + i, p);
            });
            std::vector<std::pair<std::size_t, std::uint32_t>> want;
            kmer::for_each_canonical_kmer(s.data(), s.size(), bg::kK, [&](std::uint64_t c, std::size_t p, bool) {
                want.emplace_back(p, bg::partition_of(c, 64));
            });
            ok = ok && got == want;
        }
        check(ok, "super-k-mers cover every k-mer once, in the partition its own minimizer gives");
    }

    // A segment from a synthetic genome with repeats: every k-mer's count, absent k-mers, both strands.
    {
        std::string genome = random_sequence(rng, 300000);
        const std::string unit = random_sequence(rng, 500);
        for (int i = 0; i < 400; ++i) genome += unit;          // a tandem repeat: counts above 255
        genome += random_sequence(rng, 100000) + std::string(50, 'N') + random_sequence(rng, 50000);
        const std::string fa = (dir / "g.fa").string();
        {
            std::ofstream f(fa);
            f << ">c1\n";
            for (std::size_t i = 0; i < genome.size(); i += 60) f << genome.substr(i, 60) << '\n';
            f << ">c2\n" << genome.substr(1000, 40000) << '\n';   // a second contig repeating a stretch
        }
        std::map<std::uint64_t, std::uint64_t> truth;
        for (const std::string& s : {genome, genome.substr(1000, 40000)})
            kmer::for_each_canonical_kmer(s.data(), s.size(), bg::kK, [&](std::uint64_t c, std::size_t, bool) { ++truth[c]; });
        bg::BuildOptions o;
        o.genomes = {{"g", fa}};
        o.segment_dir = (dir / "seg").string();
        o.tmp_dir = dir.string();
        o.partitions = 256;
        o.threads = 2;
        o.spill_files = 8;
        const bg::SegmentStats st = bg::build_segment(o);
        const bg::Segment seg(o.segment_dir);
        bool counts_ok = seg.distinct() == truth.size();
        std::uint64_t big = 0;
        for (const auto& [c, n] : truth) {
            counts_ok = counts_ok && seg.count(c) == n;
            big += n >= 255;
        }
        bool absent_ok = true;
        for (int t = 0; t < 20000; ++t) {
            const std::uint64_t c = kmer::canonical(rng() & kmer::mask_of(bg::kK), bg::kK);
            absent_ok = absent_ok && seg.count(c) == (truth.count(c) ? truth[c] : 0);
        }
        std::uint64_t dumped = 0;
        bool dump_ok = true;
        for (std::uint32_t p = 0; p < seg.partitions(); ++p)
            seg.for_each(p, [&](std::uint64_t c, std::uint64_t n) { ++dumped; dump_ok = dump_ok && truth[c] == n && bg::partition_of(c, 256) == p; });
        check(st.spill_files > 1, "the spill can be spread over several files (" + std::to_string(st.spill_files) + ")");
        check(counts_ok, "every k-mer's count equals brute force (" + std::to_string(truth.size()) + " k-mers, " +
                         std::to_string(big) + " counted 255 or more)");
        check(absent_ok, "random k-mers: present ones counted, absent ones 0");
        check(dump_ok && dumped == truth.size(), "iterating the partitions gives every k-mer once, in its own partition");

        // The same genome as two segments (split with a k-1 overlap): summed, and merged, equal to the whole.
        const std::filesystem::path store = dir / "store";
        std::filesystem::remove_all(store);
        const std::size_t cut = 200000;
        const std::string parts[2] = {genome.substr(0, cut), genome.substr(cut - (bg::kK - 1))};
        for (int i = 0; i < 2; ++i) {
            const std::string pfa = (dir / ("part" + std::to_string(i) + ".fa")).string();
            std::ofstream f(pfa);
            f << ">p\n" << parts[i] << '\n';
            if (i == 1) f << ">c2\n" << genome.substr(1000, 40000) << '\n';
            f.close();
            bg::BuildOptions po = o;
            po.genomes = {{"p" + std::to_string(i), pfa}};
            po.segment_dir = (store / ("seg-" + std::to_string(i))).string();
            bg::build_segment(po);
        }
        const bg::Store two(store.string());
        bool sum_ok = two.segments() == 2;
        for (const auto& [c, n] : truth) sum_ok = sum_ok && two.count(c) == n;
        bg::MergeOptions mo;
        mo.segment_dirs = two.segment_dirs();
        mo.out_segment_dir = (dir / "merged").string();
        mo.threads = 3;
        const bg::SegmentStats ms = bg::merge_segments(mo);
        const bg::Segment merged(mo.out_segment_dir);
        bool merge_ok = ms.distinct == truth.size();
        for (const auto& [c, n] : truth) merge_ok = merge_ok && merged.count(c) == n;
        check(sum_ok, "two segments of one genome: counts summed over them equal the whole genome's");
        check(merge_ok, "their merge is one segment with the same counts");

        // Compacted to SSHash: the same counts through the shared interface, 0 for absent k-mers.
        const std::filesystem::path one = dir / "one";
        std::filesystem::remove_all(one);
        std::filesystem::create_directories(one);
        std::filesystem::rename(mo.out_segment_dir, one / "seg-0");
        bg::CompactOptions co;
        co.store_dir = one.string();
        co.out_path = (dir / "one.sshash").string();
        co.tmp_dir = (dir / "compact.tmp").string();
        co.threads = 2;
        const bg::CompactStats cs = bg::compact_store(co);
        const auto compacted = bg::Background::open(co.out_path);
        const auto as_store = bg::Background::open(one.string());
        bool compact_ok = cs.kmers == truth.size() && compacted->distinct() == truth.size();
        for (const auto& [c, n] : truth) compact_ok = compact_ok && compacted->count(c) == n && as_store->count(c) == n;
        bool compact_absent = true;
        for (int t = 0; t < 20000; ++t) {
            const std::uint64_t c = kmer::canonical(rng() & kmer::mask_of(bg::kK), bg::kK);
            if (!truth.count(c)) compact_absent = compact_absent && compacted->count(c) == 0;
        }
        check(compact_ok, "compacted to SSHash (" + std::to_string(cs.strings) + " strings): every count unchanged");
        check(compact_absent, "compacted: absent k-mers give 0");
    }

    std::printf("\nbackground_test: %s\n", failures ? (std::to_string(failures) + " check(s) failed").c_str() : "all checks passed");
    return failures ? 1 : 0;
}
