// Unit tests of the background store: Elias-Fano, count encoding, radix sort, super-k-mers, segments
// counted under a memory bound (with and without folding identical super-k-mers) against brute force,
// stores of several segments, range merges, both string builders of compaction, graph counting and
// anchoring against the spelled haplotypes, incomplete files, and locating copies.
//   background_test <scratch dir>

#include "panvar/anchors.hpp"
#include "panvar/background.hpp"
#include "panvar/background_graph.hpp"
#include "panvar/elias_fano.hpp"
#include "panvar/kmer.hpp"

#include <algorithm>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <functional>
#include <map>
#include <random>
#include <set>
#include <string>
#include <tuple>
#include <vector>

using namespace panvar;
namespace fs = std::filesystem;

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

std::string rc(const std::string& x) {
    std::string r(x.rbegin(), x.rend());
    for (auto& c : r) c = c == 'A' ? 'T' : c == 'C' ? 'G' : c == 'G' ? 'C' : c == 'T' ? 'A' : 'N';
    return r;
}

void write_fasta(const std::string& path, const std::vector<std::pair<std::string, std::string>>& contigs) {
    std::ofstream f(path);
    for (const auto& [name, seq] : contigs) {
        f << '>' << name << '\n';
        for (std::size_t i = 0; i < seq.size(); i += 70) f << seq.substr(i, 70) << '\n';
    }
}

using Counts = std::map<std::uint64_t, std::uint64_t>;
void add_counts(Counts& c, const std::string& s) {
    kmer::for_each_canonical_kmer(s.data(), s.size(), bg::kK, [&](std::uint64_t x, std::size_t, bool) { ++c[x]; });
}

Counts dump(const bg::Segment& seg) {
    Counts out;
    for (std::uint32_t p = 0; p < seg.partitions(); ++p) seg.for_each(p, [&](std::uint64_t c, std::uint64_t n) { out[c] = n; });
    return out;
}

// A store folder holding the given segment folders (moved in), with a manifest listing them.
void make_store(const fs::path& store, std::uint32_t partitions, const std::vector<std::string>& segs) {
    fs::remove_all(store);
    fs::create_directories(store);
    bg::Manifest m;
    m.store_id = "test";
    m.partitions = partitions;
    for (const auto& s : segs) {
        const std::string dir = "seg-" + std::to_string(m.next_segment++);
        fs::rename(s, store / dir);
        const bg::SegmentInfo info = bg::read_segment_info((store / dir).string());
        m.segments.push_back({dir, info.distinct, info.occurrences, info.genomes});
    }
    bg::write_manifest(store.string(), m);
}

} // namespace

int main(int argc, char** argv) {
    if (argc < 2) { std::fprintf(stderr, "usage: background_test <scratch dir>\n"); return 2; }
    const fs::path dir = argv[1];
    fs::remove_all(dir);
    fs::create_directories(dir);
    std::mt19937_64 rng(31);

    // Elias-Fano: membership, order and the cursor on random sets of several sizes and universes.
    {
        bool ok = true;
        for (const std::size_t n : {0, 1, 2, 7, 100, 5000}) {
            for (const std::size_t ub : {8, 20, 62}) {
                std::set<std::uint64_t> s;
                const std::uint64_t umax = 1ULL << ub;
                while (s.size() < std::min<std::uint64_t>(n, umax)) s.insert(rng() % umax);
                const std::vector<std::uint64_t> v(s.begin(), s.end());
                const std::vector<std::uint64_t> w = ef::build(v, ub);
                const ef::View e(w.data());
                ok = ok && e.size() == v.size() && e.words() == w.size();
                for (std::size_t i = 0; i < v.size(); ++i) ok = ok && e.find(v[i]) == static_cast<std::int64_t>(i) && e.at(i) == v[i];
                std::vector<std::uint64_t> seen;
                ef::View::Cursor cur(e);
                for (std::uint64_t x; cur.next(x);) seen.push_back(x);
                ok = ok && seen == v;
                for (int t = 0; t < 2000; ++t) {
                    const std::uint64_t q = rng() % umax;
                    ok = ok && (e.find(q) >= 0) == (s.count(q) > 0);
                }
            }
        }
        bool refused = false;
        try { ef::build({1}, 64); } catch (const std::exception&) { refused = true; }
        check(ok && refused, "Elias-Fano: find, at and the cursor equal brute force; a 64-bit universe is refused");
    }

    // Counts: random, skewed and constant lists decode exactly, by index and by cursor.
    {
        bool ok = true;
        for (int t = 0; t < 60; ++t) {
            const std::size_t n = t == 0 ? 0 : 1 + rng() % 5000;
            std::vector<std::uint64_t> c(n);
            for (auto& x : c) {
                const int kind = t % 4;
                x = kind == 0 ? 1 + rng() % 3 : kind == 1 ? (rng() % 50 == 0 ? 1 + rng() % 100000 : 1 + rng() % 4)
                  : kind == 2 ? 466 : (rng() % 1000 == 0 ? (1ULL << 40) + rng() % 1000 : 400 + rng() % 120);
            }
            const auto w = bg::counts::encode(c);
            const bg::counts::View v(w.data(), n);
            ok = ok && v.words() == w.size();
            bg::counts::View::Cursor cur(v);
            for (std::size_t i = 0; i < n; ++i) ok = ok && v.get(i) == c[i] && cur.next() == c[i];
        }
        check(ok, "count encoding: exact by index and in sequence (small, skewed, constant and huge counts)");
        std::vector<std::uint64_t> c(10000, 466);
        check(bg::counts::View(bg::counts::encode(c).data(), c.size()).width() == 9, "a partition of counts near 466 costs 9 bits each");
    }

    // Radix sort equals std::sort.
    {
        bool ok = true;
        for (const std::size_t n : {0, 1, 10, 300, 100000}) {
            std::vector<std::uint64_t> v(n);
            for (auto& x : v) x = (rng() % 3 == 0) ? rng() % 50 : rng() & kmer::mask_of(31);
            auto w = v;
            std::sort(w.begin(), w.end());
            bg::radix_sort(v.data(), v.size());
            ok = ok && v == w;
        }
        check(ok, "radix sort equals std::sort");
    }

    // partition_of (reverse complement read once) equals the per-window definition.
    {
        bool ok = true;
        for (int t = 0; t < 20000; ++t) {
            const std::uint64_t code = rng() & kmer::mask_of(bg::kK);
            std::uint64_t best = ~0ULL;
            for (std::size_t j = 0; j + bg::kM <= bg::kK; ++j) {
                const std::uint64_t m = (code >> (2 * (bg::kK - bg::kM - j))) & kmer::mask_of(bg::kM);
                best = std::min(best, bg::minimizer_hash(kmer::canonical(m, bg::kM)));
            }
            ok = ok && bg::partition_of(code, 4096) == static_cast<std::uint32_t>(best & 4095) &&
                 bg::partition_of(code, 4096) == bg::partition_of(kmer::reverse_complement(code, bg::kK), 4096);
        }
        check(ok, "partition_of equals the minimizer definition, on both strands");
    }

    // Super-k-mers: they cover every k-mer exactly once, each in the partition partition_of gives it.
    {
        bool ok = true;
        for (int t = 0; t < 20; ++t) {
            std::string s = random_sequence(rng, 200 + rng() % 5000);
            if (t % 3 == 0) for (int i = 0; i < 5; ++i) s[rng() % s.size()] = 'N';
            std::vector<std::pair<std::size_t, std::uint32_t>> got, want;
            bg::for_each_super_kmer(s.data(), s.size(), 64, [&](std::size_t start, std::size_t len, std::uint32_t p) {
                for (std::size_t i = 0; i + bg::kK <= len; ++i) got.emplace_back(start + i, p);
            });
            kmer::for_each_canonical_kmer(s.data(), s.size(), bg::kK, [&](std::uint64_t c, std::size_t p, bool) {
                want.emplace_back(p, bg::partition_of(c, 64));
            });
            ok = ok && got == want;
        }
        check(ok, "super-k-mers cover every k-mer once, in the partition its own minimizer gives");
    }

    // Segments from a synthetic genome with repeats, counted several ways, all equal to brute force.
    std::string genome = random_sequence(rng, 300000);
    {
        const std::string unit = random_sequence(rng, 500);
        for (int i = 0; i < 400; ++i) genome += unit;          // a tandem repeat: counts above 255
        genome += random_sequence(rng, 100000) + std::string(50, 'N') + random_sequence(rng, 50000);
    }
    const std::string second = genome.substr(1000, 40000);
    const std::string fa = (dir / "g.fa").string();
    write_fasta(fa, {{"c1", genome}, {"c2", second}});
    {   // one very long line: chunks must stay bounded
        std::ofstream f(dir / "long.fa");
        f << ">c1\n" << genome << "\n>c2\n" << second << '\n';
    }
    Counts truth;
    add_counts(truth, genome);
    add_counts(truth, second);
    std::uint64_t big = 0;
    for (const auto& [c, n] : truth) big += n >= 255;
    bg::BuildOptions o;
    o.genomes = {{"g", fa, 0, ""}};
    o.tmp_dir = dir.string();
    o.partitions = 256;
    o.threads = 2;
    {
        struct Variant { std::string name; std::uint32_t spill; std::size_t memory; bool dedup; std::string file; };
        const std::vector<Variant> variants = {
            {"8 spill files, folding", 8, std::size_t{256} << 20, true, fa},
            {"1 spill file, no folding", 1, std::size_t{256} << 20, false, fa},
            {"a 2 MB bound (sorted runs on disk)", 4, std::size_t{2} << 20, true, fa},
            {"one line per contig", 8, std::size_t{256} << 20, true, (dir / "long.fa").string()},
        };
        for (std::size_t v = 0; v < variants.size(); ++v) {
            bg::BuildOptions vo = o;
            vo.segment_dir = (dir / ("seg" + std::to_string(v))).string();
            vo.spill_files = variants[v].spill;
            vo.memory_bytes = variants[v].memory;
            vo.dedup = variants[v].dedup;
            vo.genomes[0].path = variants[v].file;
            const bg::SegmentStats st = bg::build_segment(vo);
            const bg::Segment seg(vo.segment_dir);
            std::uint64_t occ = 0;
            for (const auto& [c, n] : truth) occ += n;
            bool ok = seg.distinct() == truth.size() && dump(seg) == truth && st.kmers == occ;
            for (const auto& [c, n] : truth) ok = ok && seg.count(c) == n;
            for (int t = 0; t < 5000; ++t) {
                const std::uint64_t c = kmer::canonical(rng() & kmer::mask_of(bg::kK), bg::kK);
                ok = ok && seg.count(c) == (truth.count(c) ? truth.at(c) : 0);
            }
            const bool runs = v != 2 || st.run_bytes > 0;
            check(ok && runs, "counts equal brute force with " + variants[v].name + " (" + std::to_string(truth.size()) + " k-mers, " +
                              std::to_string(big) + " at 255 or more" + (v == 2 ? ", " + std::to_string(st.run_bytes) + " bytes of runs" : "") + ")");
        }
    }

    // A genome longer than the reading chunk (1 Mb at a small memory bound), as one line and as 60 columns:
    // the same counts, every k-mer once, and the same anchors.
    {
        const std::string big = random_sequence(rng, 3200000);
        const std::string fa60 = (dir / "big60.fa").string(), fa1 = (dir / "big1.fa").string();
        write_fasta(fa60, {{"b", big}});
        { std::ofstream f(fa1); f << ">b\n" << big << '\n'; }
        const anchors::ReferenceAnchors ra = anchors::build_reference("b", fa60, 2);
        std::vector<std::vector<std::tuple<std::uint32_t, std::uint32_t, bool>>> places(2);
        std::vector<Counts> counts(2);
        std::vector<std::uint64_t> kmers(2);
        for (int v = 0; v < 2; ++v) {
            bg::BuildOptions bo = o;
            bo.genomes = {{"b", v == 0 ? fa60 : fa1, 0, ""}};
            bo.segment_dir = (dir / ("big-seg" + std::to_string(v))).string();
            bo.memory_bytes = std::size_t{16} << 20;
            bo.anchor_sets = {{"b", &ra}};
            kmers[v] = bg::build_segment(bo).kmers;
            counts[v] = dump(bg::Segment(bo.segment_dir));
            const anchors::SegmentAnchors sa((fs::path(bo.segment_dir) / "anchors-b.bin").string());
            std::vector<anchors::Place> ps;
            for (std::size_t a = 0; a < sa.anchors(); ++a) {
                sa.places(a, ps);
                for (const auto& p : ps) places[v].emplace_back(static_cast<std::uint32_t>(a), p.pos, p.repeated);
            }
        }
        bool at_ref = places[0].size() == ra.code.size();
        for (std::size_t i = 0; i < places[0].size() && at_ref; ++i) at_ref = std::get<1>(places[0][i]) == ra.at[std::get<0>(places[0][i])].pos;
        check(kmers[0] == big.size() - 30 && kmers[1] == kmers[0] && counts[0] == counts[1] && places[0] == places[1] && at_ref,
              "a 3.2 Mb genome read in 1 Mb chunks, one line or 60 columns: every k-mer once, every anchor where the reference has it");
    }

    // The same genome twice: folding stores identical super-k-mers once; counts double.
    {
        const std::string fa2 = (dir / "g2.fa").string();
        fs::copy_file(fa, fa2, fs::copy_options::overwrite_existing);
        bg::BuildOptions vo = o;
        vo.genomes = {{"a", fa, 0, ""}, {"b", fa2, 0, ""}};
        vo.segment_dir = (dir / "twice").string();
        vo.spill_files = 4;
        const bg::SegmentStats st = bg::build_segment(vo);
        vo.dedup = false;
        vo.segment_dir = (dir / "twice-nofold").string();
        const bg::SegmentStats st2 = bg::build_segment(vo);
        Counts doubled = truth;
        for (auto& [c, n] : doubled) n *= 2;
        const bg::Segment seg((dir / "twice").string());
        const bg::SegmentInfo info = bg::read_segment_info((dir / "twice").string());
        check(dump(seg) == doubled && st.spill_bytes * 3 < st2.spill_bytes * 2,
              "two identical genomes: counts double, folded spill " + std::to_string(st.spill_bytes) + " bytes against " + std::to_string(st2.spill_bytes));
        check(info.genomes.size() == 2 && info.genomes[0].checksum == info.genomes[1].checksum && info.genomes[0].bytes > 0,
              "identical files get the same sequence checksum and their size recorded");
    }

    // Stores: two segments summed; merged in three ranges into one segment of three parts files.
    {
        const std::size_t cut = 200000;
        const std::string parts[2] = {genome.substr(0, cut), genome.substr(cut - (bg::kK - 1))};
        std::vector<std::string> segs;
        for (int i = 0; i < 2; ++i) {
            const std::string pfa = (dir / ("part" + std::to_string(i) + ".fa")).string();
            if (i == 0) write_fasta(pfa, {{"p", parts[0]}});
            else write_fasta(pfa, {{"p", parts[1]}, {"c2", second}});
            bg::BuildOptions po = o;
            po.genomes = {{"p" + std::to_string(i), pfa, 0, ""}};
            po.segment_dir = (dir / ("built" + std::to_string(i))).string();
            bg::build_segment(po);
            segs.push_back(po.segment_dir);
        }
        const fs::path store = dir / "store";
        make_store(store, 256, segs);
        const bg::Store two(store.string());
        bool sum_ok = two.segments() == 2;
        for (const auto& [c, n] : truth) sum_ok = sum_ok && two.count(c) == n;
        std::vector<std::uint64_t> qs;
        for (const auto& [c, n] : truth) qs.push_back(c);
        for (int t = 0; t < 1000; ++t) qs.push_back(kmer::canonical(rng() & kmer::mask_of(bg::kK), bg::kK));
        const auto batch = two.count_batch(qs);
        bool batch_ok = true;
        for (std::size_t i = 0; i < qs.size(); ++i) batch_ok = batch_ok && batch[i] == (truth.count(qs[i]) ? truth.at(qs[i]) : 0);
        check(sum_ok, "two segments of one genome: counts summed over them equal the whole genome's");
        check(batch_ok, "batched lookups equal single lookups");
        const fs::path merged = dir / "merged";
        fs::create_directories(merged);
        bg::SegmentInfo info;
        info.partitions = 256;
        for (std::uint32_t i = 0; i < 3; ++i) {
            bg::MergeRangeOptions mo;
            mo.segment_dirs = {(store / "seg-0").string(), (store / "seg-1").string()};
            mo.first = 256 * i / 3;
            mo.end = 256 * (i + 1) / 3;
            mo.out_path = (merged / ("parts-" + std::to_string(i) + ".bin")).string();
            mo.threads = 3;
            info.distinct += bg::merge_range(mo).distinct;
            info.parts.push_back({"parts-" + std::to_string(i) + ".bin", mo.first, mo.end});
        }
        bg::write_segment_info(merged.string(), info);
        const bg::Segment m(merged.string());
        check(info.distinct == truth.size() && dump(m) == truth, "their merge, in three partition ranges, has the same counts");

        // Damaged files are refused, each on its own: no magic (an interrupted write), truncated, a corrupt
        // partition header.
        auto damaged = [&](const std::string& name, const std::function<void(const fs::path&)>& harm) {
            const fs::path b = dir / name;
            fs::remove_all(b);
            fs::create_directories(b);
            for (const char* f : {"parts-0.bin", "parts-1.bin", "parts-2.bin"}) fs::copy_file(merged / f, b / f);
            bg::write_segment_info(b.string(), info);
            harm(b / "parts-1.bin");
            try { bg::Segment x(b.string()); } catch (const std::exception&) { return true; }
            return false;
        };
        const bool no_magic = damaged("broken1", [](const fs::path& f) {
            std::fstream io(f, std::ios::in | std::ios::out | std::ios::binary);
            const std::uint64_t zero = 0;
            io.write(reinterpret_cast<const char*>(&zero), 8);
        });
        const bool truncated = damaged("broken2", [](const fs::path& f) { fs::resize_file(f, fs::file_size(f) - 16); });
        const bool bad_header = damaged("broken3", [](const fs::path& f) {
            // the first partition's Elias-Fano high-word count, made inconsistent
            std::fstream io(f, std::ios::in | std::ios::out | std::ios::binary);
            std::uint64_t first_off = 0;
            io.seekg(9 * 8);
            io.read(reinterpret_cast<char*>(&first_off), 8);
            io.seekp(static_cast<std::streamoff>((first_off + 3) * 8));
            const std::uint64_t bogus = 1ULL << 40;
            io.write(reinterpret_cast<const char*>(&bogus), 8);
        });
        check(no_magic && truncated && bad_header, "a parts file without its magic, truncated, or with a corrupt partition header is refused");

        // Compaction, both string builders: the same counts, 0 for absent k-mers, tied to the store.
        const fs::path one = dir / "one";
        make_store(one, 256, {merged.string()});
        for (const bool local : {false, true}) {
            bg::CompactOptions co;
            co.store_dir = one.string();
            co.out_path = (dir / "one.sshash").string();
            co.tmp_dir = (dir / "compact.tmp").string();
            co.threads = 2;
            co.local_strings = local;
            co.quiet = true;
            const bg::CompactStats cs = bg::compact_store(co);
            const auto compacted = bg::Background::open(co.out_path);
            bool ok = cs.kmers == truth.size() && compacted->distinct() == truth.size();
            for (const auto& [c, n] : truth) ok = ok && compacted->count(c) == n;
            for (int t = 0; t < 20000; ++t) {
                const std::uint64_t c = kmer::canonical(rng() & kmer::mask_of(bg::kK), bg::kK);
                if (!truth.count(c)) ok = ok && compacted->count(c) == 0;
            }
            check(ok, std::string("compacted (") + (local ? "partition-local" : "greedy") + " strings: " + std::to_string(cs.strings) +
                      " strings, " + std::to_string(cs.string_bases) + " bases): every count unchanged, absent k-mers 0");
        }
        {   // a store whose genomes changed makes its dictionary stale; a new generation alone (a reference
            // added, a merge) does not
            bg::Manifest mm = bg::read_manifest(one.string());
            ++mm.generation;
            bg::write_manifest(one.string(), mm);
            bool kept = true;
            try { bg::Background::open((dir / "one.sshash").string()); } catch (const std::exception&) { kept = false; }
            mm.segments[0].genomes.push_back({"another", "/elsewhere.fa", 1, "0000000000000000"});
            bg::write_manifest(one.string(), mm);
            bool stale = false;
            try { bg::Background::open((dir / "one.sshash").string()); } catch (const std::exception&) { stale = true; }
            check(kept && stale, "a dictionary is refused as stale when its store's genomes changed, not when only its generation did");
        }
        {   // an empty store compacts to an empty dictionary
            const std::string efa = (dir / "empty.fa").string();
            write_fasta(efa, {{"e", "NNNNACGT"}});
            bg::BuildOptions eo = o;
            eo.genomes = {{"e", efa, 0, ""}};
            eo.segment_dir = (dir / "empty-seg").string();
            bg::build_segment(eo);
            const fs::path es = dir / "empty-store";
            make_store(es, 256, {eo.segment_dir});
            bg::CompactOptions co;
            co.store_dir = es.string();
            co.out_path = (dir / "empty-store.sshash").string();
            co.tmp_dir = (dir / "compact.tmp").string();
            bg::compact_store(co);
            const auto e = bg::Background::open(co.out_path);
            check(e->distinct() == 0 && e->count(12345) == 0, "an empty store compacts to an empty background");
        }
    }

    // Graphs: counting the haplotypes of a GFA without spelling them equals counting them spelled; anchors too.
    {
        // a chain of short pieces (1 to 60 bp, so most k-mers cross several nodes) with alternatives
        std::vector<std::string> names, seqs;
        std::vector<std::vector<std::size_t>> options;   // per site, its node indexes
        for (int site = 0; site < 1500; ++site) {
            options.emplace_back();
            const int alts = rng() % 4 == 0 ? 2 : 1;
            for (int a = 0; a < alts; ++a) {
                names.push_back(std::to_string(names.size() + 1));
                seqs.push_back(random_sequence(rng, 1 + rng() % 60));
                if (rng() % 50 == 0) seqs.back()[0] = 'N';
                options.back().push_back(names.size() - 1);
            }
        }
        std::ofstream g(dir / "g.gfa");
        g << "H\tVN:Z:1.1\n";
        for (std::size_t i = 0; i < names.size(); ++i) g << "S\t" << names[i] << '\t' << seqs[i] << '\n';
        std::vector<std::vector<std::size_t>> haps;
        for (int h = 0; h < 6; ++h) {
            haps.emplace_back();
            for (const auto& opt : options) haps.back().push_back(opt[h == 0 ? 0 : rng() % opt.size()]);
        }
        for (int h = 0; h < 4; ++h) {   // P lines; hap 2 written in reverse
            g << "P\thap" << h << "\t";
            if (h == 2) for (std::size_t i = haps[h].size(); i-- > 0;) g << names[haps[h][i]] << '-' << (i ? "," : "");
            else for (std::size_t i = 0; i < haps[h].size(); ++i) g << names[haps[h][i]] << '+' << (i + 1 < haps[h].size() ? "," : "");
            g << "\t*\n";
        }
        for (int h = 4; h < 6; ++h) {   // W lines, the second starting at 1000 and visiting a stretch twice
            std::vector<std::size_t> walk = haps[h];
            if (h == 5) walk.insert(walk.begin() + 200, haps[h].begin() + 100, haps[h].begin() + 300);
            g << "W\tS" << h << "\t1\tchr1\t" << (h == 5 ? 1000 : 0) << "\t*\t";
            for (const std::size_t n : walk) g << '>' << names[n];
            g << '\n';
        }
        g.close();
        const std::string gfa = (dir / "g.gfa").string(), spelled = (dir / "g.spelled.fa").string();
        bg::graph::spell_graph(gfa, spelled);
        bg::BuildOptions go = o;
        go.genomes = {{"graph", gfa, 0, ""}};
        go.segment_dir = (dir / "graph-seg").string();
        bg::build_segment(go);
        go.genomes = {{"spelled", spelled, 0, ""}};
        go.segment_dir = (dir / "spelled-seg").string();
        bg::build_segment(go);
        const Counts from_graph = dump(bg::Segment((dir / "graph-seg").string()));
        const Counts from_fasta = dump(bg::Segment((dir / "spelled-seg").string()));
        check(!from_graph.empty() && from_graph == from_fasta, "a graph's haplotypes counted from nodes and junctions equal their spelled count (" +
                                                               std::to_string(from_graph.size()) + " k-mers)");
        // anchors: the first haplotype spelled is the reference
        std::string hap0;
        for (const std::size_t n : haps[0]) hap0 += seqs[n];
        const std::string rfa = (dir / "hap0.fa").string();
        write_fasta(rfa, {{"ref", hap0}});
        const anchors::ReferenceAnchors ra = anchors::build_reference("ref", rfa, 2);
        const anchors::AnchorLookup lookup(ra);
        auto key = [](const anchors::Hit& h) { return std::make_tuple(h.anchor, h.contig, h.pos, h.forward); };
        std::vector<std::tuple<std::uint32_t, std::uint32_t, std::uint32_t, bool>> a, b;
        for (const auto& h : bg::graph::graph_anchors(gfa, lookup, 0, 2)) a.push_back(key(h));
        std::vector<std::string> contigs;
        std::vector<std::uint64_t> lengths;
        for (auto h : anchors::map_fasta(lookup, spelled, 0, 2, contigs, lengths)) { if (h.contig == 5) h.pos += 1000; b.push_back(key(h)); }
        std::sort(a.begin(), a.end());
        std::sort(b.begin(), b.end());
        check(!a.empty() && a == b && contigs == bg::graph::graph_contigs(gfa),
              "anchors placed through the graph equal those found in the spelled haplotypes (" + std::to_string(a.size()) + ")");
        // through the counting pass, the anchor file and locate: every haplotype holds its copy
        bg::BuildOptions ga = o;
        ga.genomes = {{"graph", gfa, 0, ""}};
        ga.segment_dir = (dir / "graph-anchored").string();
        ga.anchor_sets = {{"ref", &ra}};
        bg::build_segment(ga);
        const anchors::SegmentAnchors gsa((dir / "graph-anchored" / "anchors-ref.bin").string());
        const auto gspans = anchors::locate(ra, {&gsa}, "ref", 10000, 30000);
        std::set<std::string> found;
        for (const auto& sp : gspans) found.insert(sp.contig);
        check(found.size() == 6, "located through a graph's anchors: a copy in each of its 6 haplotypes (" + std::to_string(found.size()) + ")");
    }

    // Locating: moved, reversed, split, tandem-duplicated and displaced copies; short intervals.
    {
        const std::string ref = random_sequence(rng, 2000000);
        const std::size_t ls = 900000, le = 1000000;
        const std::string b1 = ref.substr(0, 300000) + random_sequence(rng, 50000) + ref.substr(300000);
        const std::string b2 = random_sequence(rng, 200000) + rc(ref.substr(ls - 30000, (le - ls) + 60000)) + random_sequence(rng, 200000);
        const std::string b3a = random_sequence(rng, 100000) + ref.substr(ls - 20000, 70000);
        const std::string b3b = ref.substr(ls + 50000, 70000) + random_sequence(rng, 100000);
        // tandem: the locus twice in a row on one contig
        const std::string t1 = ref.substr(0, le) + ref.substr(ls, le - ls) + ref.substr(le);
        // displaced: the locus cut out and put back 1 Mb further on, the flanks joined where it was
        const std::string d1 = ref.substr(0, ls) + ref.substr(le, 900000) + ref.substr(ls, le - ls) + ref.substr(le + 900000);
        const std::string rfa = (dir / "ref.fa").string(), bfa = (dir / "b.fa").string(), tfa = (dir / "t.fa").string(), dfa = (dir / "d.fa").string();
        write_fasta(rfa, {{"chrR", ref}});
        write_fasta(bfa, {{"b1", b1}, {"b2", b2}, {"b3a", b3a}, {"b3b", b3b}});
        write_fasta(tfa, {{"t1", t1}});
        write_fasta(dfa, {{"d1", d1}});
        const anchors::ReferenceAnchors ra = anchors::build_reference("ref", rfa, 2);
        bg::BuildOptions ao = o;
        ao.genomes = {{"ref", rfa, 0, ""}, {"b", bfa, 0, ""}, {"t", tfa, 0, ""}, {"d", dfa, 0, ""}};
        ao.segment_dir = (dir / "anchor-seg").string();
        ao.anchor_sets = {{"ref", &ra}};
        bg::build_segment(ao);
        const anchors::SegmentAnchors sa((dir / "anchor-seg" / "anchors-ref.bin").string());
        const std::vector<std::string> want_genomes = {"ref", "b", "t", "d"};
        check(sa.fingerprint() == ra.fingerprint() && sa.genomes() == want_genomes, "the counting pass writes the segment's anchor file");
        std::uint64_t bytes = 0;
        const auto spans = anchors::locate(ra, {&sa}, "chrR", ls, le, {}, &bytes);
        auto of = [&](const std::string& g) { std::vector<anchors::Span> v; for (const auto& s : spans) if (s.genome == g) v.push_back(s); return v; };
        auto near = [](std::uint64_t a, std::uint64_t b) { return (a > b ? a - b : b - a) < 5000; };
        auto has = [&](const std::string& g, const std::string& contig, bool fwd, std::uint64_t s0, std::uint64_t e0) {
            for (const auto& sp : of(g)) if (sp.contig == contig && sp.forward == fwd && near(sp.start, s0) && near(sp.end, e0)) return true;
            return false;
        };
        check(ra.code.size() > 1500 && ra.code.size() < 2500, "anchors: about one per kilobase (" + std::to_string(ra.code.size()) + " on 2 Mb)");
        check(of("ref").size() == 1 && has("ref", "chrR", true, ls, le), "the reference locates its own interval");
        check(has("b", "b1", true, ls + 50000, le + 50000), "a moved copy is found where it moved");
        check(has("b", "b2", false, 230000, 230000 + (le - ls)), "a reversed copy elsewhere is found, on the other strand");
        check(has("b", "b3a", true, 120000, 170000) && has("b", "b3b", true, 0, 50000), "a copy split by a contig break is found as two partial copies");
        const auto tt = of("t");
        bool overlap = false;
        for (std::size_t i = 0; i < tt.size(); ++i)
            for (std::size_t j = i + 1; j < tt.size(); ++j)
                if (tt[i].contig == tt[j].contig && tt[i].start < tt[j].end && tt[j].start < tt[i].end) overlap = true;
        std::size_t chains = 0;
        for (const auto& s : tt) chains += s.chains;
        bool covers = false;
        for (const auto& s : tt) if (s.start <= ls + 5000 && s.end + 5000 >= le + (le - ls)) covers = true;
        bool both = has("t", "t1", true, ls, le) && has("t", "t1", true, le, le + (le - ls));
        check(!overlap && (both || (covers && chains == 2)), "a tandem duplicate gives two copies that do not overlap (" + std::to_string(tt.size()) + " span(s), " +
                                                              std::to_string(chains) + " chains)");
        check(has("d", "d1", true, ls + 900000, le + 900000), "a copy displaced from its flanks is still found (the flank-only chain does not stop the search)");
        check(bytes > 0 && bytes < 200000, "locating reads only the window's anchors (" + std::to_string(bytes) + " bytes)");
        // a short interval with fewer than three anchors inside
        std::size_t a0 = 0;
        while (ra.at[a0].pos < 500000) ++a0;
        const std::uint64_t s0 = ra.at[a0].pos, e0 = ra.at[a0].pos + 40;
        const auto small = anchors::locate(ra, {&sa}, "chrR", s0, e0);
        bool small_ok = false;
        for (const auto& sp : small) if (sp.genome == "ref" && sp.start <= s0 && sp.end >= e0 && sp.end - sp.start < 10000) small_ok = true;
        check(small_ok, "an interval with one anchor inside is still located");
        bool bad_refused = false;
        try { anchors::locate(ra, {&sa}, "chrR", le, ls); } catch (const std::exception&) { bad_refused = true; }
        try { anchors::locate(ra, {&sa}, "chrR", 0, 3000000); bad_refused = false; } catch (const std::exception&) {}
        check(bad_refused, "reversed or out-of-contig intervals are refused");

        // merging anchor files keeps every genome's places
        const std::string mpath = (dir / "anchors-merged.bin").string();
        anchors::merge_segment_anchors({(dir / "anchor-seg" / "anchors-ref.bin").string(), (dir / "anchor-seg" / "anchors-ref.bin").string()}, ra, mpath);
        const anchors::SegmentAnchors ms(mpath);
        const auto mspans = anchors::locate(ra, {&ms}, "chrR", ls, le);
        check(ms.genomes().size() == 8 && mspans.size() == 2 * spans.size(), "merged anchor files keep every genome's places");
    }

    // a haplotype contig held as several walks: an interval across them is cut from each, a gap refused
    {
        const std::string gfa = (dir / "split_walks.gfa").string();
        std::ofstream(gfa) << "H\tVN:Z:1.1\nS\ta\tACGTACGTAC\nS\tb\tGGGGGCCCCC\nS\tc\tTTTTTAAAAA\n"
                           << "W\ts\t1\tchrX\t0\t20\t>a>b\nW\ts\t1\tchrX\t20\t30\t>c\n"
                           << "W\ts\t2\tchrY\t0\t10\t>a\nW\ts\t2\tchrY\t15\t25\t<b\n";
        const bg::graph::Haplotypes h(gfa);
        const bool across = h.fetch("s#1#chrX", 5, 25) == "CGTACGGGGGCCCCCTTTTT" && h.fetch("s#1#chrX", 22, 30) == "TTTAAAAA";
        bool gap = false;
        try { h.fetch("s#2#chrY", 5, 20); } catch (const std::exception&) { gap = true; }
        check(across && gap, "graph intervals are cut across consecutive walks of a contig; a gap between walks is refused");
    }

    std::printf("\nbackground_test: %s\n", failures ? (std::to_string(failures) + " check(s) failed").c_str() : "all checks passed");
    return failures ? 1 : 0;
}
