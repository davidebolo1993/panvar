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
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include <zlib.h>

namespace panvar::bg {

namespace {

namespace fs = std::filesystem;

template <class F>
void parallel_partitions(std::uint32_t n, std::size_t threads, F&& f) {
    std::atomic<std::uint32_t> next{0};
    std::exception_ptr failure;
    std::mutex fm;
    std::vector<std::thread> pool;
    for (std::size_t t = 0; t < std::max<std::size_t>(1, threads); ++t)
        pool.emplace_back([&] {
            try {
                for (std::uint32_t i; (i = next++) < n;) {
                    { std::lock_guard<std::mutex> lk(fm); if (failure) return; }
                    f(i);
                }
            } catch (...) {
                std::lock_guard<std::mutex> lk(fm);
                if (!failure) failure = std::current_exception();
            }
        });
    for (auto& t : pool) t.join();
    if (failure) std::rethrow_exception(failure);
}

class StoreBackground final : public Background {
public:
    explicit StoreBackground(const std::string& dir) : store_(dir) {
        std::uint64_t n = 0;
        for (const auto& d : store_.segment_dirs()) n += Segment(d).distinct();
        distinct_ = n;   // an upper bound when segments share k-mers
    }
    std::uint64_t count(std::uint64_t c) const override { return store_.count(c); }
    std::uint64_t distinct() const override { return distinct_; }
    std::string kind() const override { return "store, " + std::to_string(store_.segments()) + " segment(s)"; }

private:
    Store store_;
    std::uint64_t distinct_ = 0;
};

class SshashBackground final : public Background {
public:
    explicit SshashBackground(const std::string& path) {
        essentials::load(dict_, path.c_str());
        if (dict_.k() != kK) throw std::runtime_error("background: " + path + " is not a 31-mer dictionary");
    }
    std::uint64_t count(std::uint64_t c) const override {
        const std::string s = kmer::decode(c, kK);
        const auto r = dict_.lookup(s.c_str(), true);
        return r.kmer_id == sshash::constants::invalid_uint64 ? 0 : dict_.weight(r.kmer_id);
    }
    std::uint64_t distinct() const override { return dict_.num_kmers(); }
    std::string kind() const override { return "compacted (SSHash)"; }

private:
    sshash::dictionary_type dict_;
};

}  // namespace

std::unique_ptr<Background> Background::open(const std::string& path) {
    if (fs::is_directory(path)) return std::make_unique<StoreBackground>(path);
    if (fs::is_regular_file(path)) return std::make_unique<SshashBackground>(path);
    throw std::runtime_error("background: no store or compacted background at " + path);
}

CompactStats compact_store(const CompactOptions& o) {
    auto log = [&](const std::string& m) { if (o.log) o.log(m); };
    const auto dirs = segment_dirs_of(o.store_dir);
    if (dirs.size() != 1)
        throw std::runtime_error("background: compacting needs a store of one segment (run --merge first)");
    const Segment seg(dirs.front());
    const std::uint32_t parts = seg.partitions();
    std::vector<std::uint64_t> base(parts + 1, 0);
    for (std::uint32_t p = 0; p < parts; ++p) base[p + 1] = base[p] + seg.partition_size(p);
    const std::uint64_t n = base[parts];
    std::vector<std::atomic<std::uint64_t>> visited((n + 63) / 64);
    for (auto& w : visited) w.store(0, std::memory_order_relaxed);
    auto claim = [&](std::uint64_t rank) {   // true when this call is the one that visits it
        const std::uint64_t bit = 1ULL << (rank % 64);
        return (visited[rank / 64].fetch_or(bit, std::memory_order_relaxed) & bit) == 0;
    };
    // The rank and count of a canonical k-mer, or false when absent.
    auto find = [&](std::uint64_t code, std::uint64_t& rank, std::uint64_t& count) {
        const std::uint32_t p = partition_of(code, parts);
        const std::int64_t i = seg.find(code, p, &count);
        if (i < 0) return false;
        rank = base[p] + static_cast<std::uint64_t>(i);
        return true;
    };

    // The string set: greedy simplitigs, each k-mer once, with its count.
    fs::create_directories(o.tmp_dir);
    const std::string spss = (fs::path(o.tmp_dir) / "compact.spss.fa.gz").string();
    gzFile gz = gzopen(spss.c_str(), "wb1");
    if (gz == nullptr) throw std::runtime_error("background: cannot write " + spss);
    std::mutex gz_m;
    std::atomic<std::uint64_t> ids{0}, strings{0}, bases_written{0}, kmers_written{0};
    const std::uint64_t kmask = kmer::mask_of(kK);
    log("compact: writing the string set of " + std::to_string(n) + " k-mers");
    try {
        parallel_partitions(parts, o.threads, [&](std::uint32_t p) {
            std::string buf, seq, left;
            std::vector<std::uint64_t> weights, left_w;
            auto flush = [&] {
                std::lock_guard<std::mutex> lk(gz_m);
                if (!buf.empty() && gzwrite(gz, buf.data(), static_cast<unsigned>(buf.size())) != static_cast<int>(buf.size()))
                    throw std::runtime_error("background: cannot write " + spss);
                buf.clear();
            };
            std::uint64_t i = 0;
            seg.for_each(p, [&](std::uint64_t code, std::uint64_t count) {
                const std::uint64_t rank = base[p] + i++;
                if (!claim(rank)) return;
                seq = kmer::decode(code, kK);
                weights.assign(1, count);
                // right: the last k-mer's forward code, extended one base at a time
                std::uint64_t last = code, r = 0, c = 0;
                while (true) {
                    bool grew = false;
                    for (std::uint64_t b = 0; b < 4 && !grew; ++b) {
                        const std::uint64_t nxt = ((last << 2) | b) & kmask;
                        if (find(kmer::canonical(nxt, kK), r, c) && claim(r)) {
                            seq += "ACGT"[b];
                            weights.push_back(c);
                            last = nxt;
                            grew = true;
                        }
                    }
                    if (!grew) break;
                }
                // left: the first k-mer's forward code, extended backwards, collected reversed
                std::uint64_t first = code;
                left.clear();
                left_w.clear();
                while (true) {
                    bool grew = false;
                    for (std::uint64_t b = 0; b < 4 && !grew; ++b) {
                        const std::uint64_t prv = (b << (2 * (kK - 1))) | (first >> 2);
                        if (find(kmer::canonical(prv, kK), r, c) && claim(r)) {
                            left += "ACGT"[b];
                            left_w.push_back(c);
                            first = prv;
                            grew = true;
                        }
                    }
                    if (!grew) break;
                }
                std::reverse(left.begin(), left.end());
                std::reverse(left_w.begin(), left_w.end());
                const std::uint64_t id = ids++;
                buf += '>' + std::to_string(id) + " LN:i:" + std::to_string(left.size() + seq.size()) + " ab:Z:";
                for (std::size_t w = 0; w < left_w.size(); ++w) { buf += std::to_string(left_w[w]); buf += ' '; }
                for (std::size_t w = 0; w < weights.size(); ++w) { buf += std::to_string(weights[w]); if (w + 1 < weights.size()) buf += ' '; }
                buf += '\n';
                buf += left;
                buf += seq;
                buf += '\n';
                ++strings;
                bases_written += left.size() + seq.size();
                kmers_written += left_w.size() + weights.size();
                if (buf.size() >= (std::size_t{4} << 20)) flush();
            });
            flush();
        });
    } catch (...) {
        gzclose(gz);
        std::error_code ec;
        fs::remove(spss, ec);
        throw;
    }
    if (gzclose(gz) != Z_OK) throw std::runtime_error("background: cannot write " + spss);
    if (kmers_written != n) throw std::runtime_error("background: the string set holds " + std::to_string(kmers_written.load()) +
                                                     " k-mers, the store " + std::to_string(n));
    CompactStats st;
    st.kmers = kmers_written;
    st.strings = strings;
    st.string_bases = bases_written;
    log("compact: " + std::to_string(st.strings) + " strings, " + std::to_string(st.string_bases) + " bases (" +
        std::to_string(static_cast<double>(st.string_bases) / static_cast<double>(std::max<std::uint64_t>(1, n))) +
        " bases per k-mer)");

    // SSHash: its minimizer length grows with the k-mer count (about log4 of it, plus a margin).
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
    {
        sshash::dictionary_type dict;
        dict.build(spss, cfg);
        essentials::save(dict, o.out_path.c_str());
    }
    std::error_code ec;
    fs::remove(spss, ec);
    st.bytes = fs::file_size(o.out_path);
    return st;
}

}  // namespace panvar::bg
