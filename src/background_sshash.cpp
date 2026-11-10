// The compacted background: a merged store rewritten as a string set and handed to SSHash, and the one
// lookup interface over both forms. The only file that includes SSHash's headers.

#include "panvar/background.hpp"

#include "include/dictionary_types.hpp"
// SSHash's templates are instantiated by including its sources, as its own tools do.
#include "src/builder/build.cpp"
#include "src/dictionary.cpp"
#include "src/info.cpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <mutex>
#include <numeric>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include <fcntl.h>
#include <unistd.h>
#include <zlib.h>

namespace panvar::bg {

namespace {

namespace fs = std::filesystem;

constexpr std::uint64_t kDictMagic = 0x3148535347425650ULL;   // "PVBGSSH1"

template <class F>
void parallel_partitions(std::uint32_t n, std::size_t threads, F&& f) {
    std::atomic<std::uint32_t> next{0};
    std::atomic<bool> failed{false};
    std::exception_ptr failure;
    std::mutex fm;
    std::vector<std::thread> pool;
    for (std::size_t t = 0; t < std::max<std::size_t>(1, threads); ++t)
        pool.emplace_back([&, t] {
            try {
                for (std::uint32_t i; !failed && (i = next++) < n;) f(i, t);
            } catch (...) {
                std::lock_guard<std::mutex> lk(fm);
                if (!failure) failure = std::current_exception();
                failed = true;
            }
        });
    for (auto& t : pool) t.join();
    if (failure) std::rethrow_exception(failure);
}

// panvar's header in front of the SSHash payload: magic, version, store id, generation, k-mers, payload start.
struct DictHeader {
    std::string store_id;
    std::uint64_t generation = 0, kmers = 0, payload = 0;
};

DictHeader read_header(std::istream& in, const std::string& path) {
    auto u64 = [&] {
        std::uint64_t v = 0;
        if (!in.read(reinterpret_cast<char*>(&v), 8)) throw std::runtime_error("background: " + path + " is truncated");
        return v;
    };
    if (u64() != kDictMagic) throw std::runtime_error("background: " + path + " is not a compacted background of this version (or was interrupted)");
    if (u64() != 1) throw std::runtime_error("background: " + path + ": unknown version");
    DictHeader h;
    const std::uint64_t idlen = u64();
    if (idlen > 256) throw std::runtime_error("background: " + path + " is corrupt");
    h.store_id.resize(idlen);
    in.read(&h.store_id[0], static_cast<std::streamsize>(idlen));
    in.seekg(static_cast<std::streamoff>((32 + idlen + 7) / 8 * 8));
    h.generation = u64();
    h.kmers = u64();
    h.payload = u64();
    return h;
}

void write_header(std::ostream& out, const DictHeader& h, std::uint64_t magic) {
    auto u64 = [&](std::uint64_t v) { out.write(reinterpret_cast<const char*>(&v), 8); };
    u64(magic);
    u64(1);
    u64(h.store_id.size());
    out.write(h.store_id.data(), static_cast<std::streamsize>(h.store_id.size()));
    for (std::size_t pad = (32 + h.store_id.size() + 7) / 8 * 8 - (24 + h.store_id.size()); pad > 0; --pad) out.put('\0');
    u64(h.generation);
    u64(h.kmers);
    u64(h.payload);
}

std::uint64_t header_size(const std::string& store_id) { return (32 + store_id.size() + 7) / 8 * 8 + 24; }

class StoreBackground final : public Background {
public:
    explicit StoreBackground(const std::string& dir) : store_(dir) {}
    std::uint64_t count(std::uint64_t c) const override { return store_.count(c); }
    std::vector<std::uint64_t> count_batch(const std::vector<std::uint64_t>& c) const override { return store_.count_batch(c); }
    std::uint64_t distinct() const override {
        if (store_.segments() == 1) return store_.segment(0).distinct();
        std::uint64_t n = 0;
        for (std::uint32_t p = 0; p < store_.partitions(); ++p) store_.for_each(p, [&](std::uint64_t, std::uint64_t) { ++n; });
        return n;
    }
    std::string kind() const override { return "store, " + std::to_string(store_.segments()) + " segment(s)"; }
    std::string store_id() const override { return store_.manifest().store_id; }
    std::uint64_t generation() const override { return store_.manifest().generation; }

private:
    Store store_;
};

// SSHash's k-mer: base i in bits 2i, with A, C, G, T as 0, 1, 3, 2.
sshash::dictionary_type::kmer_type to_sshash(std::uint64_t code) {
    std::uint64_t x = 0;
    for (std::size_t i = 0; i < kK; ++i) {
        const std::uint64_t c = (code >> (2 * (kK - 1 - i))) & 3;
        x |= (c ^ (c >> 1)) << (2 * i);
    }
    return sshash::dictionary_type::kmer_type(x);
}

class SshashBackground final : public Background {
public:
    explicit SshashBackground(const std::string& path) : map_(std::make_shared<Mapping>(path)) {
        std::ifstream in(path, std::ios::binary);
        if (!in) throw std::runtime_error("background: cannot read " + path);
        header_ = read_header(in, path);
        if (header_.kmers == 0) return;   // an empty background
        in.seekg(static_cast<std::streamoff>(header_.payload));
        essentials::generic_loader loader(in);
        std::shared_ptr<const void> owner(map_, map_->data());
        loader.set_mmap(static_cast<const std::uint8_t*>(map_->data()), map_->size(), owner);
        loader.visit(dict_);
        if (dict_.k() != kK) throw std::runtime_error("background: " + path + " is not a 31-mer dictionary");
        if (dict_.num_kmers() != header_.kmers) throw std::runtime_error("background: " + path + " is corrupt (k-mer count)");
    }
    std::uint64_t count(std::uint64_t c) const override {
        if (header_.kmers == 0) return 0;
        const auto r = dict_.lookup(to_sshash(c), true);
        return r.kmer_id == sshash::constants::invalid_uint64 ? 0 : dict_.weight(r.kmer_id);
    }
    std::uint64_t distinct() const override { return header_.kmers; }
    std::string kind() const override { return "compacted (SSHash)"; }
    std::string store_id() const override { return header_.store_id; }
    std::uint64_t generation() const override { return header_.generation; }

private:
    std::shared_ptr<Mapping> map_;
    DictHeader header_;
    sshash::dictionary_type dict_;
};

// Writes the weighted string set: ">id LN:i:L ab:Z:w1 w2 ..." then the sequence, gzip level 1.
class SpssWriter {
public:
    explicit SpssWriter(const std::string& path) : path_(path) {
        gz_ = gzopen(path.c_str(), "wb1");
        if (gz_ == nullptr) throw std::runtime_error("background: cannot write " + path);
    }
    ~SpssWriter() { if (gz_ != nullptr) gzclose(gz_); }
    void put(const std::string& seq, const std::vector<std::uint64_t>& w) {
        std::string& b = buf_;
        b += '>' + std::to_string(strings_) + " LN:i:" + std::to_string(seq.size()) + " ab:Z:";
        for (std::size_t i = 0; i < w.size(); ++i) { b += std::to_string(w[i]); if (i + 1 < w.size()) b += ' '; }
        b += '\n';
        b += seq;
        b += '\n';
        ++strings_;
        bases_ += seq.size();
        kmers_ += w.size();
        if (b.size() >= (std::size_t{4} << 20)) flush();
    }
    void flush() {
        if (!buf_.empty() && gzwrite(gz_, buf_.data(), static_cast<unsigned>(buf_.size())) != static_cast<int>(buf_.size()))
            throw std::runtime_error("background: cannot write " + path_);
        buf_.clear();
    }
    void close() {
        flush();
        const int rc = gzclose(gz_);
        gz_ = nullptr;
        if (rc != Z_OK) throw std::runtime_error("background: cannot write " + path_);
    }
    std::uint64_t strings_ = 0, bases_ = 0, kmers_ = 0;
private:
    std::string path_, buf_;
    gzFile gz_ = nullptr;
};

// The earlier builder: greedy simplitigs extended across the whole store (random access to every partition).
void greedy_strings(const Segment& seg, std::size_t threads, SpssWriter& out) {
    const std::uint32_t parts = seg.partitions();
    std::vector<std::uint64_t> base(parts + 1, 0);
    for (std::uint32_t p = 0; p < parts; ++p) base[p + 1] = base[p] + seg.partition_size(p);
    const std::uint64_t n = base[parts];
    std::vector<std::atomic<std::uint64_t>> visited((n + 63) / 64);
    for (auto& w : visited) w.store(0, std::memory_order_relaxed);
    auto claim = [&](std::uint64_t rank) {
        const std::uint64_t bit = 1ULL << (rank % 64);
        return (visited[rank / 64].fetch_or(bit, std::memory_order_relaxed) & bit) == 0;
    };
    auto find = [&](std::uint64_t code, std::uint64_t& rank, std::uint64_t& count) {
        const std::uint32_t p = partition_of(code, parts);
        const std::int64_t i = seg.find(code, p, &count);
        if (i < 0) return false;
        rank = base[p] + static_cast<std::uint64_t>(i);
        return true;
    };
    std::mutex om;
    const std::uint64_t kmask = kmer::mask_of(kK);
    parallel_partitions(parts, threads, [&](std::uint32_t p, std::size_t) {
        std::string seq, left;
        std::vector<std::uint64_t> weights, left_w, all_w;
        std::uint64_t i = 0;
        seg.for_each(p, [&](std::uint64_t code, std::uint64_t count) {
            const std::uint64_t rank = base[p] + i++;
            if (!claim(rank)) return;
            seq = kmer::decode(code, kK);
            weights.assign(1, count);
            std::uint64_t last = code, r = 0, c = 0;
            for (bool grew = true; grew;) {
                grew = false;
                for (std::uint64_t b = 0; b < 4 && !grew; ++b) {
                    const std::uint64_t nxt = ((last << 2) | b) & kmask;
                    if (find(kmer::canonical(nxt, kK), r, c) && claim(r)) { seq += "ACGT"[b]; weights.push_back(c); last = nxt; grew = true; }
                }
            }
            std::uint64_t first = code;
            left.clear();
            left_w.clear();
            for (bool grew = true; grew;) {
                grew = false;
                for (std::uint64_t b = 0; b < 4 && !grew; ++b) {
                    const std::uint64_t prv = (b << (2 * (kK - 1))) | (first >> 2);
                    if (find(kmer::canonical(prv, kK), r, c) && claim(r)) { left += "ACGT"[b]; left_w.push_back(c); first = prv; grew = true; }
                }
            }
            std::reverse(left.begin(), left.end());
            std::reverse(left_w.begin(), left_w.end());
            all_w = left_w;
            all_w.insert(all_w.end(), weights.begin(), weights.end());
            std::lock_guard<std::mutex> lk(om);
            out.put(left + seq, all_w);
        });
    });
}

// The partition-local builder: within each partition (in memory, no access to other partitions), greedy
// paths over its k-mers; then paths are glued end to end across partitions through a sorted table of their
// end k-mers. Equal-count neighbours are preferred, so counts change less often along a string.
void local_strings(const Segment& seg, std::size_t threads, SpssWriter& out, const std::function<void(const std::string&)>& log) {
    const std::uint32_t parts = seg.partitions();
    const std::uint64_t kmask = kmer::mask_of(kK);
    struct Frag {
        std::uint32_t arena;
        std::uint64_t offset;
        std::uint64_t left, right;       // first and last k-mer, read forward
        std::uint64_t left_w, right_w;   // their counts
    };
    std::vector<std::string> arenas(std::max<std::size_t>(1, threads));
    std::vector<std::vector<Frag>> frags(arenas.size());
    auto put_varint = [](std::string& b, std::uint64_t v) { while (v >= 0x80) { b += static_cast<char>((v & 0x7f) | 0x80); v >>= 7; } b += static_cast<char>(v); };

    parallel_partitions(parts, arenas.size(), [&](std::uint32_t p, std::size_t t) {
        std::vector<std::uint64_t> codes, cnt;
        seg.for_each(p, [&](std::uint64_t c, std::uint64_t n) { codes.push_back(c); cnt.push_back(n); });
        std::vector<char> used(codes.size(), 0);
        auto index = [&](std::uint64_t canonical) -> std::int64_t {
            const auto it = std::lower_bound(codes.begin(), codes.end(), canonical);
            return it != codes.end() && *it == canonical ? it - codes.begin() : -1;
        };
        std::vector<std::uint64_t> fwd_codes, fwd_w, back_codes, back_w;
        for (std::size_t s = 0; s < codes.size(); ++s) {
            if (used[s]) continue;
            used[s] = 1;
            // extend from the seed in one direction: `next` gives the oriented neighbour code for base b
            auto grow = [&](std::uint64_t from, std::uint64_t w, bool right, std::vector<std::uint64_t>& oc, std::vector<std::uint64_t>& ow) {
                oc.clear(); ow.clear();
                std::uint64_t cur = from, curw = w;
                while (true) {
                    std::int64_t pick = -1;
                    std::uint64_t pick_code = 0;
                    for (std::uint64_t b = 0; b < 4; ++b) {
                        const std::uint64_t nb = right ? ((cur << 2) | b) & kmask : (b << (2 * (kK - 1))) | (cur >> 2);
                        const std::int64_t i = index(kmer::canonical(nb, kK));
                        if (i < 0 || used[i]) continue;
                        if (pick < 0 || (cnt[i] == curw && cnt[pick] != curw)) { pick = i; pick_code = nb; }
                    }
                    if (pick < 0) break;
                    used[pick] = 1;
                    oc.push_back(pick_code);
                    ow.push_back(cnt[pick]);
                    cur = pick_code;
                    curw = cnt[pick];
                }
            };
            grow(codes[s], cnt[s], true, fwd_codes, fwd_w);
            grow(codes[s], cnt[s], false, back_codes, back_w);
            // the path read forward: back (reversed), seed, fwd
            std::vector<std::uint64_t> path, w;
            for (std::size_t i = back_codes.size(); i-- > 0;) { path.push_back(back_codes[i]); w.push_back(back_w[i]); }
            path.push_back(codes[s]); w.push_back(cnt[s]);
            path.insert(path.end(), fwd_codes.begin(), fwd_codes.end());
            w.insert(w.end(), fwd_w.begin(), fwd_w.end());
            std::string& a = arenas[t];
            Frag f{static_cast<std::uint32_t>(t), a.size(), path.front(), path.back(), w.front(), w.back()};
            put_varint(a, path.size());
            // bases: the first k-mer, then one base per further k-mer; 2 bits each
            const std::size_t nb = kK + path.size() - 1;
            const std::size_t at = a.size();
            a.resize(at + (nb + 3) / 4, '\0');
            auto set_base = [&](std::size_t i, std::uint64_t b) { a[at + i / 4] = static_cast<char>(static_cast<unsigned char>(a[at + i / 4]) | (b << (2 * (i % 4)))); };
            for (std::size_t i = 0; i < kK; ++i) set_base(i, (path[0] >> (2 * (kK - 1 - i))) & 3);
            for (std::size_t i = 1; i < path.size(); ++i) set_base(kK + i - 1, path[i] & 3);
            for (const std::uint64_t x : w) put_varint(a, x);
            frags[t].push_back(f);
        }
    });

    std::vector<Frag> all;
    for (auto& v : frags) { all.insert(all.end(), v.begin(), v.end()); std::vector<Frag>().swap(v); }
    const std::uint64_t nf = all.size();
    if (log) log("compact: " + std::to_string(nf) + " partition-local paths; gluing their ends");
    // end e = frag << 1 | side (0 left, 1 right); its outward code reads away from the path
    auto out_code = [&](std::uint64_t e) {
        const Frag& f = all[e >> 1];
        return (e & 1) ? f.right : kmer::reverse_complement(f.left, kK);
    };
    auto end_w = [&](std::uint64_t e) { const Frag& f = all[e >> 1]; return (e & 1) ? f.right_w : f.left_w; };
    std::vector<std::pair<std::uint64_t, std::uint64_t>> ends(2 * nf);
    for (std::uint64_t e = 0; e < 2 * nf; ++e) ends[e] = {kmer::canonical(out_code(e), kK), e};
    std::sort(ends.begin(), ends.end());
    std::vector<std::uint64_t> link(2 * nf, ~0ULL);
    std::vector<std::uint64_t> parent(nf);
    std::iota(parent.begin(), parent.end(), 0);
    auto root = [&](std::uint64_t x) { while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; } return x; };
    for (std::uint64_t e = 0; e < 2 * nf; ++e) {
        if (link[e] != ~0ULL) continue;
        const std::uint64_t o = out_code(e);
        std::uint64_t pick = ~0ULL;
        for (std::uint64_t b = 0; b < 4; ++b) {
            const std::uint64_t target = kmer::reverse_complement(((o << 2) | b) & kmask, kK);   // the other end's outward code
            const std::uint64_t key = kmer::canonical(target, kK);
            for (auto it = std::lower_bound(ends.begin(), ends.end(), std::make_pair(key, std::uint64_t{0})); it != ends.end() && it->first == key; ++it) {
                const std::uint64_t e2 = it->second;
                if (link[e2] != ~0ULL || out_code(e2) != target || root(e2 >> 1) == root(e >> 1)) continue;
                if (pick == ~0ULL || (end_w(e2) == end_w(e) && end_w(pick) != end_w(e))) pick = e2;
            }
        }
        if (pick == ~0ULL) continue;
        link[e] = pick;
        link[pick] = e;
        parent[root(e >> 1)] = root(pick >> 1);
    }
    std::vector<std::pair<std::uint64_t, std::uint64_t>>().swap(ends);

    // walk every chain from a free end, decoding paths in the orientation the chain reads them
    std::vector<char> done(nf, 0);
    auto get_varint = [](const unsigned char*& p) { std::uint64_t v = 0; for (int s = 0;; s += 7) { const unsigned char b = *p++; v |= static_cast<std::uint64_t>(b & 0x7f) << s; if (!(b & 0x80)) return v; } };
    std::string seq, piece;
    std::vector<std::uint64_t> w, pw;
    auto decode = [&](std::uint64_t fi, bool reverse) {
        const Frag& f = all[fi];
        const auto* p = reinterpret_cast<const unsigned char*>(arenas[f.arena].data() + f.offset);
        const std::uint64_t nk = get_varint(p);
        const std::size_t nb = kK + nk - 1;
        piece.resize(nb);
        for (std::size_t i = 0; i < nb; ++i) piece[i] = "ACGT"[(p[i / 4] >> (2 * (i % 4))) & 3];
        p += (nb + 3) / 4;
        pw.resize(nk);
        for (auto& x : pw) x = get_varint(p);
        if (reverse) {
            std::reverse(piece.begin(), piece.end());
            for (char& c : piece) c = c == 'A' ? 'T' : c == 'C' ? 'G' : c == 'G' ? 'C' : 'A';
            std::reverse(pw.begin(), pw.end());
        }
    };
    auto emit_chain = [&](std::uint64_t start_end) {
        seq.clear();
        w.clear();
        std::uint64_t e = start_end;   // the end we enter through
        while (true) {
            const std::uint64_t fi = e >> 1;
            done[fi] = 1;
            const bool reverse = (e & 1) == 1;   // entering through the right end reads it backwards
            decode(fi, reverse);
            if (seq.empty()) seq = piece;
            else seq.append(piece, kK - 1, std::string::npos);
            w.insert(w.end(), pw.begin(), pw.end());
            const std::uint64_t exit = (fi << 1) | (reverse ? 0 : 1);
            if (link[exit] == ~0ULL) break;
            e = link[exit];
        }
        out.put(seq, w);
    };
    for (std::uint64_t fi = 0; fi < nf; ++fi) {
        if (done[fi]) continue;
        if (link[fi << 1] == ~0ULL) emit_chain(fi << 1);
        else if (link[(fi << 1) | 1] == ~0ULL) emit_chain((fi << 1) | 1);
    }
    for (std::uint64_t fi = 0; fi < nf; ++fi)
        if (!done[fi]) throw std::runtime_error("background: compaction left a cycle of paths");
}

}  // namespace

std::unique_ptr<Background> Background::open(const std::string& path) {
    if (fs::is_directory(path)) return std::make_unique<StoreBackground>(path);
    if (!fs::is_regular_file(path)) throw std::runtime_error("background: no store or compacted background at " + path);
    auto dict = std::make_unique<SshashBackground>(path);
    // a dictionary next to its store must match it
    std::string base = path;
    if (base.size() > 7 && base.compare(base.size() - 7, 7, ".sshash") == 0) {
        base.resize(base.size() - 7);
        if (store_exists(base)) {
            const Manifest m = read_manifest(base);
            if (m.store_id != dict->store_id())
                throw std::runtime_error("background: " + path + " was compacted from another store than " + base);
            if (m.generation != dict->generation())
                throw std::runtime_error("background: " + path + " is stale: " + base + " changed since it was compacted (compact again)");
        }
    }
    return dict;
}

CompactStats compact_store(const CompactOptions& o) {
    auto log = [&](const std::string& m) { if (o.log) o.log(m); };
    const Manifest m = read_manifest(o.store_dir);
    if (m.segments.size() != 1)
        throw std::runtime_error("background: compacting needs a store of exactly one segment (run --merge first)");
    const Segment seg((fs::path(o.store_dir) / m.segments[0].dir).string());
    const std::uint64_t n = seg.distinct();
    CompactStats st;
    const std::string tmp_out = o.out_path + ".tmp";
    DictHeader h;
    h.store_id = m.store_id;
    h.generation = m.generation;
    h.kmers = n;
    h.payload = header_size(h.store_id);

    if (n == 0) {
        std::ofstream out(tmp_out, std::ios::binary);
        write_header(out, h, kDictMagic);
        out.flush();
        if (!out) throw std::runtime_error("background: cannot write " + tmp_out);
        out.close();
        fs::rename(tmp_out, o.out_path);
        st.bytes = fs::file_size(o.out_path);
        log("compact: the store holds no k-mer; wrote an empty background");
        return st;
    }

    fs::create_directories(o.tmp_dir);
    const std::string spss = (fs::path(o.tmp_dir) / "compact.spss.fa.gz").string();
    struct Cleanup { std::string a, b; ~Cleanup() { std::error_code ec; fs::remove(a, ec); fs::remove(b, ec); } } cleanup{spss, tmp_out};
    const auto t0 = std::chrono::steady_clock::now();
    log("compact: writing the string set of " + std::to_string(n) + " k-mers (" + (o.local_strings ? "partition-local paths, glued" : "greedy over the store") + ")");
    {
        SpssWriter w(spss);
        if (o.local_strings) local_strings(seg, o.threads, w, log);
        else greedy_strings(seg, o.threads, w);
        w.close();
        st.strings = w.strings_;
        st.string_bases = w.bases_;
        st.kmers = w.kmers_;
    }
    if (st.kmers != n)
        throw std::runtime_error("background: the string set holds " + std::to_string(st.kmers) + " k-mers, the store " + std::to_string(n));
    st.strings_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    log("compact: " + std::to_string(st.strings) + " strings, " + std::to_string(st.string_bases) + " bases (" +
        std::to_string(static_cast<double>(st.string_bases) / static_cast<double>(n)) + " bases per k-mer) in " +
        std::to_string(st.strings_seconds) + " s");

    st.minimizer = o.minimizer != 0 ? o.minimizer
        : std::clamp<std::size_t>(static_cast<std::size_t>(std::ceil(std::log(static_cast<double>(std::max<std::uint64_t>(n, 2))) / std::log(4.0))) + 4, 13, 21);
    sshash::build_configuration cfg;
    cfg.k = kK;
    cfg.m = st.minimizer;
    cfg.weighted = true;
    cfg.verbose = false;
    cfg.num_threads = std::max<std::size_t>(1, o.threads);
    cfg.ram_limit_in_GiB = std::max<std::size_t>(1, o.ram_gb);
    cfg.tmp_dirname = o.tmp_dir;
    log("compact: building the SSHash dictionary (minimizer length " + std::to_string(st.minimizer) + ")");
    const auto t1 = std::chrono::steady_clock::now();
    {
        sshash::dictionary_type dict;
        // SSHash prints a few lines to stdout: send them to stderr, or nowhere when quiet.
        std::cout.flush();
        std::fflush(stdout);
        const int saved = ::dup(1);
        const int sink = o.quiet ? ::open("/dev/null", O_WRONLY) : ::dup(2);
        if (saved >= 0 && sink >= 0) ::dup2(sink, 1);
        if (sink >= 0) ::close(sink);
        try {
            dict.build(spss, cfg);
        } catch (...) {
            std::cout.flush(); std::fflush(stdout);
            if (saved >= 0) { ::dup2(saved, 1); ::close(saved); }
            throw;
        }
        std::cout.flush(); std::fflush(stdout);
        if (saved >= 0) { ::dup2(saved, 1); ::close(saved); }
        if (dict.num_kmers() != n) throw std::runtime_error("background: SSHash holds " + std::to_string(dict.num_kmers()) + " k-mers, the store " + std::to_string(n));
        std::ofstream out(tmp_out, std::ios::binary);
        if (!out) throw std::runtime_error("background: cannot write " + tmp_out);
        write_header(out, h, 0);   // magic last
        essentials::generic_saver saver(out);
        saver.visit(dict);
        out.flush();
        if (!out) throw std::runtime_error("background: cannot write " + tmp_out);
        out.seekp(0);
        const std::uint64_t magic = kDictMagic;
        out.write(reinterpret_cast<const char*>(&magic), 8);
        out.flush();
        if (!out) throw std::runtime_error("background: cannot write " + tmp_out);
    }
    {
        const int fd = ::open(tmp_out.c_str(), O_RDONLY);
        if (fd >= 0) { ::fsync(fd); ::close(fd); }
    }
    fs::rename(tmp_out, o.out_path);
    st.sshash_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t1).count();
    st.bytes = fs::file_size(o.out_path);
    return st;
}

}  // namespace panvar::bg
