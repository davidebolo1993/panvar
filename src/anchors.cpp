#include "panvar/anchors.hpp"

#include "panvar/gz_reader.hpp"
#include "panvar/kmer.hpp"

#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <cstdio>
#include <cstring>
#include <deque>
#include <map>
#include <mutex>
#include <stdexcept>
#include <thread>
#include <unordered_map>
#include <utility>

namespace panvar::anchors {

namespace {

constexpr std::uint64_t kSampleSeed = 0x9e3779b97f4a7c15ULL;
constexpr std::uint64_t kRefMagic = 0x3146455252434E41ULL;   // "ANCRREF1"
constexpr std::uint64_t kMapMagic = 0x3150414D52434E41ULL;   // "ANCRMAP1"

// Streams a FASTA's contigs to `threads` workers: f(contig index, sequence). Contig names are collected in
// file order.
void for_each_contig(const std::string& fasta, std::size_t threads, std::vector<std::string>& names,
                     const std::function<void(std::uint32_t, const std::string&)>& f) {
    GzLineReader in(fasta);
    if (!in.ok()) throw std::runtime_error("anchors: cannot read " + fasta);
    std::deque<std::pair<std::uint32_t, std::string>> queue;
    std::mutex m;
    std::condition_variable not_empty, not_full;
    bool closed = false;
    std::exception_ptr failure;
    std::vector<std::thread> pool;
    for (std::size_t t = 0; t < std::max<std::size_t>(1, threads); ++t)
        pool.emplace_back([&] {
            while (true) {
                std::pair<std::uint32_t, std::string> job;
                {
                    std::unique_lock<std::mutex> lk(m);
                    not_empty.wait(lk, [&] { return !queue.empty() || closed; });
                    if (queue.empty()) return;
                    job = std::move(queue.front());
                    queue.pop_front();
                    not_full.notify_one();
                }
                try {
                    f(job.first, job.second);
                } catch (...) {
                    std::lock_guard<std::mutex> lk(m);
                    if (!failure) failure = std::current_exception();
                }
            }
        });
    auto push = [&](std::uint32_t id, std::string seq) {
        std::unique_lock<std::mutex> lk(m);
        not_full.wait(lk, [&] { return queue.size() < std::max<std::size_t>(1, threads); });
        queue.emplace_back(id, std::move(seq));
        not_empty.notify_one();
    };
    std::string line, seq;
    bool have = false;
    try {
        while (in.getline(line)) {
            if (line.empty()) continue;
            if (line[0] == '>') {
                if (have) push(static_cast<std::uint32_t>(names.size() - 1), std::move(seq));
                seq.clear();
                const auto sp = line.find_first_of(" \t");
                names.push_back(line.substr(1, sp == std::string::npos ? std::string::npos : sp - 1));
                have = true;
                continue;
            }
            seq += line;
        }
        if (have) push(static_cast<std::uint32_t>(names.size() - 1), std::move(seq));
    } catch (...) {
        std::lock_guard<std::mutex> lk(m);
        if (!failure) failure = std::current_exception();
    }
    {
        std::lock_guard<std::mutex> lk(m);
        closed = true;
        not_empty.notify_all();
    }
    for (auto& t : pool) t.join();
    if (failure) std::rethrow_exception(failure);
}

struct Hit {
    std::uint64_t code;
    std::uint32_t contig, pos;
    bool forward;
};

void put64(std::FILE* f, std::uint64_t v) { if (std::fwrite(&v, 8, 1, f) != 1) throw std::runtime_error("anchors: write failed"); }
std::uint64_t get64(std::FILE* f) {
    std::uint64_t v = 0;
    if (std::fread(&v, 8, 1, f) != 1) throw std::runtime_error("anchors: truncated file");
    return v;
}
void put_names(std::FILE* f, const std::vector<std::string>& names) {
    put64(f, names.size());
    for (const auto& n : names) {
        put64(f, n.size());
        if (!n.empty() && std::fwrite(n.data(), 1, n.size(), f) != n.size()) throw std::runtime_error("anchors: write failed");
    }
}
std::vector<std::string> get_names(std::FILE* f) {
    std::vector<std::string> names(get64(f));
    for (auto& n : names) {
        n.resize(get64(f));
        if (!n.empty() && std::fread(&n[0], 1, n.size(), f) != n.size()) throw std::runtime_error("anchors: truncated file");
    }
    return names;
}
template <class T>
void put_vec(std::FILE* f, const std::vector<T>& v) {
    put64(f, v.size());
    if (!v.empty() && std::fwrite(v.data(), sizeof(T), v.size(), f) != v.size()) throw std::runtime_error("anchors: write failed");
}
template <class T>
std::vector<T> get_vec(std::FILE* f) {
    std::vector<T> v(get64(f));
    if (!v.empty() && std::fread(v.data(), sizeof(T), v.size(), f) != v.size()) throw std::runtime_error("anchors: truncated file");
    return v;
}

}  // namespace

bool sampled(std::uint64_t canonical_kmer) {
    return (kmer::mix64(canonical_kmer ^ kSampleSeed) >> 32) % kSampleOneIn == 0;
}

ReferenceAnchors build_reference(const std::string& fasta, std::size_t threads,
                                 const std::function<void(const std::string&)>& log) {
    std::vector<Hit> hits;
    std::mutex hm;
    ReferenceAnchors ref;
    for_each_contig(fasta, threads, ref.contigs, [&](std::uint32_t contig, const std::string& seq) {
        std::vector<Hit> local;
        kmer::for_each_closed_syncmer(seq.data(), seq.size(), kK, kS, [&](std::uint64_t c, std::size_t p, bool fwd) {
            if (sampled(c)) local.push_back({c, contig, static_cast<std::uint32_t>(p), fwd});
        });
        std::lock_guard<std::mutex> lk(hm);
        hits.insert(hits.end(), local.begin(), local.end());
    });
    std::sort(hits.begin(), hits.end(), [](const Hit& a, const Hit& b) { return a.code < b.code; });
    std::vector<Hit> unique;
    for (std::size_t i = 0; i < hits.size();) {
        std::size_t e = i + 1;
        while (e < hits.size() && hits[e].code == hits[i].code) ++e;
        if (e - i == 1) unique.push_back(hits[i]);
        i = e;
    }
    std::sort(unique.begin(), unique.end(), [](const Hit& a, const Hit& b) {
        return a.contig != b.contig ? a.contig < b.contig : a.pos < b.pos;
    });
    for (const Hit& h : unique) {
        ref.code.push_back(h.code);
        ref.at.push_back({h.contig, h.pos, h.forward});
    }
    if (log) log("anchors: " + std::to_string(ref.code.size()) + " unique of " + std::to_string(hits.size()) + " sampled syncmers");
    return ref;
}

GenomeMap map_genome(const ReferenceAnchors& ref, const std::string& fasta, std::size_t threads) {
    std::unordered_map<std::uint64_t, std::uint32_t> id;
    id.reserve(ref.code.size() * 2);
    for (std::size_t i = 0; i < ref.code.size(); ++i) id.emplace(ref.code[i], static_cast<std::uint32_t>(i));
    std::vector<std::pair<std::uint32_t, Placement>> hits;
    std::mutex hm;
    GenomeMap map;
    for_each_contig(fasta, threads, map.contigs, [&](std::uint32_t contig, const std::string& seq) {
        std::vector<std::pair<std::uint32_t, Placement>> local;
        kmer::for_each_closed_syncmer(seq.data(), seq.size(), kK, kS, [&](std::uint64_t c, std::size_t p, bool fwd) {
            if (!sampled(c)) return;
            const auto it = id.find(c);
            if (it != id.end()) local.push_back({it->second, {contig, static_cast<std::uint32_t>(p), fwd}});
        });
        std::lock_guard<std::mutex> lk(hm);
        hits.insert(hits.end(), local.begin(), local.end());
    });
    std::sort(hits.begin(), hits.end(), [](const auto& a, const auto& b) {
        if (a.first != b.first) return a.first < b.first;
        return a.second.contig != b.second.contig ? a.second.contig < b.second.contig : a.second.pos < b.second.pos;
    });
    const std::size_t n = ref.code.size();
    map.first.assign(n + 1, 0);
    map.occurrences.assign(n, 0);
    std::size_t h = 0;
    for (std::size_t a = 0; a < n; ++a) {
        map.first[a] = static_cast<std::uint32_t>(map.places.size());
        while (h < hits.size() && hits[h].first == a) {
            if (map.occurrences[a]++ < kMaxOccurrences) map.places.push_back(hits[h].second);
            ++h;
        }
    }
    map.first[n] = static_cast<std::uint32_t>(map.places.size());
    return map;
}

void save(const ReferenceAnchors& ref, const std::string& path) {
    std::FILE* f = std::fopen(path.c_str(), "wb");
    if (f == nullptr) throw std::runtime_error("anchors: cannot write " + path);
    put64(f, kRefMagic);
    put64(f, kK); put64(f, kS); put64(f, kSampleOneIn);
    put_names(f, ref.contigs);
    put_vec(f, ref.code);
    put_vec(f, ref.at);
    if (std::fclose(f) != 0) throw std::runtime_error("anchors: cannot write " + path);
}

void save(const GenomeMap& map, const std::string& path) {
    std::FILE* f = std::fopen(path.c_str(), "wb");
    if (f == nullptr) throw std::runtime_error("anchors: cannot write " + path);
    put64(f, kMapMagic);
    put_names(f, map.contigs);
    put_vec(f, map.first);
    put_vec(f, map.places);
    put_vec(f, map.occurrences);
    if (std::fclose(f) != 0) throw std::runtime_error("anchors: cannot write " + path);
}

ReferenceAnchors load_reference(const std::string& path) {
    std::FILE* f = std::fopen(path.c_str(), "rb");
    if (f == nullptr) throw std::runtime_error("anchors: cannot read " + path);
    ReferenceAnchors ref;
    try {
        if (get64(f) != kRefMagic || get64(f) != kK || get64(f) != kS || get64(f) != kSampleOneIn)
            throw std::runtime_error("anchors: " + path + " is not a reference anchor set of this version");
        ref.contigs = get_names(f);
        ref.code = get_vec<std::uint64_t>(f);
        ref.at = get_vec<Placement>(f);
    } catch (...) { std::fclose(f); throw; }
    std::fclose(f);
    return ref;
}

GenomeMap load_map(const std::string& path) {
    std::FILE* f = std::fopen(path.c_str(), "rb");
    if (f == nullptr) throw std::runtime_error("anchors: cannot read " + path);
    GenomeMap map;
    try {
        if (get64(f) != kMapMagic) throw std::runtime_error("anchors: " + path + " is not an anchor map");
        map.contigs = get_names(f);
        map.first = get_vec<std::uint32_t>(f);
        map.places = get_vec<Placement>(f);
        map.occurrences = get_vec<std::uint32_t>(f);
    } catch (...) { std::fclose(f); throw; }
    std::fclose(f);
    return map;
}

std::vector<Span> locate(const ReferenceAnchors& ref, const GenomeMap& map, const std::string& contig,
                         std::uint64_t start, std::uint64_t end, std::size_t flank, std::size_t min_inside,
                         std::uint64_t slack) {
    const auto cit = std::find(ref.contigs.begin(), ref.contigs.end(), contig);
    if (cit == ref.contigs.end()) throw std::runtime_error("anchors: the reference has no contig " + contig);
    const auto c = static_cast<std::uint32_t>(cit - ref.contigs.begin());
    // reference anchors of this contig, by position
    auto first_at_or_after = [&](std::uint64_t pos) {
        return static_cast<std::size_t>(std::lower_bound(ref.at.begin(), ref.at.end(), std::make_pair(c, pos),
            [](const Placement& a, const std::pair<std::uint32_t, std::uint64_t>& b) {
                return a.contig != b.first ? a.contig < b.first : a.pos < b.second;
            }) - ref.at.begin());
    };
    const std::size_t c_begin = first_at_or_after(0), c_end = first_at_or_after(~0ULL >> 1);
    const std::size_t in_begin = first_at_or_after(start);
    std::size_t in_end = in_begin;
    while (in_end < c_end && ref.at[in_end].pos + kK <= end) ++in_end;
    const std::size_t w_begin = in_begin >= c_begin + flank ? in_begin - flank : c_begin;
    const std::size_t w_end = std::min(c_end, in_end + flank);

    struct Point { std::uint64_t ref_pos, pos; std::uint32_t contig; bool same; bool inside; };
    std::map<std::pair<std::uint32_t, bool>, std::vector<Point>> groups;
    for (std::size_t a = w_begin; a < w_end; ++a)
        for (std::uint32_t i = map.first[a]; i < map.first[a + 1]; ++i) {
            const Placement& p = map.places[i];
            const bool same = p.forward == ref.at[a].forward;
            groups[{p.contig, same}].push_back({ref.at[a].pos, p.pos, p.contig, same, a >= in_begin && a < in_end});
        }
    std::vector<Span> out;
    for (auto& [key, pts] : groups) {
        std::sort(pts.begin(), pts.end(), [](const Point& a, const Point& b) { return a.ref_pos != b.ref_pos ? a.ref_pos < b.ref_pos : a.pos < b.pos; });
        std::vector<char> used(pts.size(), 0);
        while (true) {
            // longest chain over unused points: increasing reference position, genome position increasing (same
            // strand) or decreasing (opposite), each gap within the reference gap's 1.5 times plus slack
            const std::size_t n = pts.size();
            std::vector<std::size_t> len(n, 0), prev(n, SIZE_MAX);
            std::size_t best = SIZE_MAX;
            for (std::size_t j = 0; j < n; ++j) {
                if (used[j]) continue;
                len[j] = 1;
                for (std::size_t i = 0; i < j; ++i) {
                    if (used[i] || pts[i].ref_pos >= pts[j].ref_pos) continue;
                    const std::uint64_t dr = pts[j].ref_pos - pts[i].ref_pos;
                    const bool order = key.second ? pts[j].pos > pts[i].pos : pts[j].pos < pts[i].pos;
                    if (!order) continue;
                    const std::uint64_t dg = key.second ? pts[j].pos - pts[i].pos : pts[i].pos - pts[j].pos;
                    if (dg > dr + dr / 2 + slack || dr > dg + dg / 2 + slack) continue;
                    if (len[i] + 1 > len[j]) { len[j] = len[i] + 1; prev[j] = i; }
                }
                if (best == SIZE_MAX || len[j] > len[best]) best = j;
            }
            if (best == SIZE_MAX) break;
            std::vector<std::size_t> chain;
            for (std::size_t j = best; j != SIZE_MAX; j = prev[j]) chain.push_back(j);
            // The span: the inside anchors and the nearest chained flank anchor on each side (the flanks make the
            // chain sure; only the nearest one bounds the copy).
            std::size_t inside = 0, flanks = 0;
            std::uint64_t lo = ~0ULL, hi = 0;
            std::size_t before = SIZE_MAX, after = SIZE_MAX;
            for (std::size_t j : chain) {
                used[j] = 1;
                if (pts[j].inside) {
                    ++inside;
                    lo = std::min(lo, pts[j].pos);
                    hi = std::max(hi, pts[j].pos + kK);
                } else {
                    ++flanks;
                    if (pts[j].ref_pos < start) { if (before == SIZE_MAX || pts[j].ref_pos > pts[before].ref_pos) before = j; }
                    else if (after == SIZE_MAX || pts[j].ref_pos < pts[after].ref_pos) after = j;
                }
            }
            for (const std::size_t j : {before, after})
                if (j != SIZE_MAX) { lo = std::min(lo, pts[j].pos); hi = std::max(hi, pts[j].pos + kK); }
            if (inside < min_inside) break;
            out.push_back({map.contigs[key.first], lo, hi, key.second, inside, flanks});
        }
    }
    std::sort(out.begin(), out.end(), [](const Span& a, const Span& b) { return a.anchors_inside > b.anchors_inside; });
    return out;
}

}  // namespace panvar::anchors
