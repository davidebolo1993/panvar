#include "panvar/background.hpp"

#include "panvar/gz_reader.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <unistd.h>

namespace panvar::bg {

namespace {

namespace fs = std::filesystem;

constexpr std::uint64_t kMagic = 0x3147455347425650ULL;   // "PVBGSEG1", little-endian
constexpr std::uint64_t kVersion = 1;
constexpr std::size_t kHeaderWords = 6;                    // magic, version, k, m, partitions, distinct
constexpr std::size_t kChunkBases = std::size_t{8} << 20;  // FASTA is cut into chunks of this many bases

std::uint64_t next_pow2(std::uint64_t x) {
    std::uint64_t p = 1;
    while (p < x) p <<= 1;
    return p;
}

std::uint32_t log2_exact(std::uint64_t x) {
    std::uint32_t l = 0;
    while ((std::uint64_t{1} << l) < x) ++l;
    return l;
}

// Every spill file stays open while the genomes are read.
void raise_open_file_limit(std::size_t files) {
    struct rlimit rl {};
    if (getrlimit(RLIMIT_NOFILE, &rl) != 0) return;
    const rlim_t want = static_cast<rlim_t>(files + 64);
    if (rl.rlim_cur >= want) return;
    rl.rlim_cur = rl.rlim_max == RLIM_INFINITY ? want : std::min<rlim_t>(want, rl.rlim_max);
    setrlimit(RLIMIT_NOFILE, &rl);
}

// The bases a genome holds, from its .fai when there is one, else estimated from the file size.
std::uint64_t estimated_bases(const std::string& path) {
    std::ifstream fai(path + ".fai");
    if (fai) {
        std::uint64_t total = 0;
        std::string line;
        while (std::getline(fai, line)) {
            const auto t1 = line.find('\t');
            if (t1 == std::string::npos) continue;
            total += std::strtoull(line.c_str() + t1 + 1, nullptr, 10);
        }
        if (total > 0) return total;
    }
    std::error_code ec;
    const std::uint64_t size = fs::file_size(path, ec);
    if (ec) return 0;
    std::ifstream in(path, std::ios::binary);
    unsigned char magic[2] = {0, 0};
    in.read(reinterpret_cast<char*>(magic), 2);
    return magic[0] == 0x1f && magic[1] == 0x8b ? size * 4 : size;   // gzip compresses DNA about 3.4-fold
}

// A spill record: partition (u16), length (u32), then the bases, 4 per byte, base i in bits 2*(i%4).
void append_record(std::string& buf, std::uint32_t part, const char* s, std::size_t len) {
    const std::uint16_t p16 = static_cast<std::uint16_t>(part);
    const std::uint32_t l32 = static_cast<std::uint32_t>(len);
    buf.append(reinterpret_cast<const char*>(&p16), 2);
    buf.append(reinterpret_cast<const char*>(&l32), 4);
    const std::size_t at = buf.size();
    buf.resize(at + (len + 3) / 4, '\0');
    for (std::size_t i = 0; i < len; ++i)
        buf[at + i / 4] = static_cast<char>(static_cast<unsigned char>(buf[at + i / 4]) |
                                            (kmer::base_code(s[i]) << (2 * (i % 4))));
}

// A bounded queue of FASTA chunks between the reader and the workers.
class ChunkQueue {
public:
    explicit ChunkQueue(std::size_t capacity) : capacity_(capacity) {}
    void push(std::string chunk) {
        std::unique_lock<std::mutex> lk(m_);
        not_full_.wait(lk, [&] { return q_.size() < capacity_; });
        q_.push_back(std::move(chunk));
        not_empty_.notify_one();
    }
    bool pop(std::string& chunk) {
        std::unique_lock<std::mutex> lk(m_);
        not_empty_.wait(lk, [&] { return !q_.empty() || closed_; });
        if (q_.empty()) return false;
        chunk = std::move(q_.front());
        q_.pop_front();
        not_full_.notify_one();
        return true;
    }
    void close() {
        std::lock_guard<std::mutex> lk(m_);
        closed_ = true;
        not_empty_.notify_all();
    }

private:
    std::size_t capacity_;
    std::deque<std::string> q_;
    bool closed_ = false;
    std::mutex m_;
    std::condition_variable not_empty_, not_full_;
};

// One partition's sorted k-mer codes and their counts, serialized: Elias-Fano codes, the counts as bytes
// (255 = see the overflow list), the overflow list (index << 32 | count).
std::vector<std::uint64_t> encode_partition(const std::vector<std::uint64_t>& codes, const std::vector<std::uint64_t>& counts) {
    std::vector<std::uint64_t> out = ef::build(codes, 2 * kK);
    const std::size_t n = codes.size(), count_words = (n + 7) / 8, at = out.size();
    out.resize(at + count_words, 0);
    auto* bytes = reinterpret_cast<std::uint8_t*>(out.data() + at);
    std::vector<std::uint64_t> overflow;
    for (std::size_t i = 0; i < n; ++i) {
        if (counts[i] >= 255) {
            bytes[i] = 255;
            if (counts[i] > 0xFFFFFFFFULL) throw std::runtime_error("background: a count above 2^32");
            overflow.push_back((static_cast<std::uint64_t>(i) << 32) | counts[i]);
        } else {
            bytes[i] = static_cast<std::uint8_t>(counts[i]);
        }
    }
    out.push_back(overflow.size());
    out.insert(out.end(), overflow.begin(), overflow.end());
    return out;
}

}  // namespace

SegmentStats build_segment(const BuildOptions& o) {
    if (o.partitions == 0 || o.partitions > 65536 || (o.partitions & (o.partitions - 1)) != 0)
        throw std::runtime_error("background: the partition count must be a power of two, at most 65536");
    if (o.genomes.empty()) throw std::runtime_error("background: no genome given");
    const std::size_t threads = std::max<std::size_t>(1, o.threads);
    auto log = [&](const std::string& m) { if (o.log) o.log(m); };

    // Spill files: as many as it takes for one to be counted in a thread's share of the memory (8 bytes per
    // k-mer occurrence, plus room for the vectors to grow).
    std::uint64_t bases = 0;
    for (const Genome& g : o.genomes) {
        const std::uint64_t b = estimated_bases(g.path);
        if (b == 0) throw std::runtime_error("background: cannot read " + g.path);
        bases += b;
    }
    const std::uint64_t per_thread = std::max<std::uint64_t>(std::uint64_t{64} << 20, o.memory_bytes / threads);
    if (o.spill_files != 0 && (o.spill_files & (o.spill_files - 1)) != 0)
        throw std::runtime_error("background: the spill file count must be a power of two");
    const std::uint32_t files = static_cast<std::uint32_t>(std::min<std::uint64_t>(
        o.partitions, o.spill_files != 0 ? o.spill_files : next_pow2((bases * 12 + per_thread - 1) / per_thread)));
    const std::uint32_t file_shift = log2_exact(files);
    raise_open_file_limit(files);

    const fs::path spill_dir = fs::path(o.tmp_dir.empty() ? fs::temp_directory_path() : fs::path(o.tmp_dir)) /
        ("panvar-background-" + std::to_string(::getpid()) + "-" +
         std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    fs::create_directories(spill_dir);
    struct SpillFile { std::FILE* fp = nullptr; std::mutex m; std::string path; };
    std::vector<SpillFile> spill(files);
    for (std::uint32_t b = 0; b < files; ++b) {
        spill[b].path = (spill_dir / ("s" + std::to_string(b) + ".spill")).string();
        spill[b].fp = std::fopen(spill[b].path.c_str(), "wb");
        if (spill[b].fp == nullptr) throw std::runtime_error("background: cannot write " + spill[b].path);
    }
    const std::size_t flush_at = std::clamp<std::size_t>((std::size_t{256} << 20) / (threads * files),
                                                         std::size_t{64} << 10, std::size_t{4} << 20);
    log("spill: " + std::to_string(files) + " files in " + spill_dir.string() + " (estimated " +
        std::to_string(bases) + " bases)");

    // Phase 1: the reader cuts the FASTA into chunks (each next chunk repeating the last k-1 bases of the
    // previous one), and workers cut chunks into super-k-mers appended to their partition's spill file.
    SegmentStats st;
    st.spill_files = files;
    std::atomic<std::uint64_t> acgt{0}, occurrences{0};
    std::exception_ptr failure;
    std::mutex failure_m;
    {
        ChunkQueue queue(threads * 2);
        std::vector<std::thread> workers;
        for (std::size_t t = 0; t < threads; ++t) {
            workers.emplace_back([&] {
                try {
                    std::vector<std::string> bufs(files);
                    auto flush = [&](std::uint32_t b) {
                        if (bufs[b].empty()) return;
                        std::lock_guard<std::mutex> lk(spill[b].m);
                        if (std::fwrite(bufs[b].data(), 1, bufs[b].size(), spill[b].fp) != bufs[b].size())
                            throw std::runtime_error("background: cannot write " + spill[b].path);
                        bufs[b].clear();
                    };
                    std::string chunk;
                    while (queue.pop(chunk)) {
                        std::uint64_t local_kmers = 0;
                        for_each_super_kmer(chunk.data(), chunk.size(), o.partitions,
                                            [&](std::size_t start, std::size_t len, std::uint32_t part) {
                            const std::uint32_t b = part & (files - 1);
                            append_record(bufs[b], part, chunk.data() + start, len);
                            local_kmers += len - kK + 1;
                            if (bufs[b].size() >= flush_at) flush(b);
                        });
                        occurrences += local_kmers;
                    }
                    for (std::uint32_t b = 0; b < files; ++b) flush(b);
                } catch (...) {
                    std::lock_guard<std::mutex> lk(failure_m);
                    if (!failure) failure = std::current_exception();
                }
            });
        }
        try {
            for (const Genome& g : o.genomes) {
                GzLineReader in(g.path);
                if (!in.ok()) throw std::runtime_error("background: cannot read " + g.path);
                std::string line, chunk;
                std::uint64_t genome_bases = 0;
                chunk.reserve(kChunkBases + kK);
                auto ship = [&](bool keep_tail) {
                    if (chunk.size() >= kK) {
                        std::string tail = keep_tail ? chunk.substr(chunk.size() - (kK - 1)) : std::string();
                        queue.push(std::move(chunk));
                        chunk = std::move(tail);
                        chunk.reserve(kChunkBases + kK);
                    } else if (!keep_tail) {
                        chunk.clear();
                    }
                };
                while (in.getline(line)) {
                    if (line.empty()) continue;
                    if (line[0] == '>') { ship(false); continue; }
                    for (char c : line) genome_bases += kmer::base_code(c) >= 0;
                    chunk += line;
                    if (chunk.size() >= kChunkBases) ship(true);
                }
                ship(false);
                acgt += genome_bases;
                log("read " + g.name + ": " + std::to_string(genome_bases) + " ACGT bases");
            }
        } catch (...) {
            std::lock_guard<std::mutex> lk(failure_m);
            if (!failure) failure = std::current_exception();
        }
        queue.close();
        for (auto& w : workers) w.join();
    }
    for (auto& s : spill) {
        if (s.fp != nullptr && std::fclose(s.fp) != 0 && !failure)
            failure = std::make_exception_ptr(std::runtime_error("background: cannot write " + s.path));
        s.fp = nullptr;
    }
    if (failure) { std::error_code ec; fs::remove_all(spill_dir, ec); std::rethrow_exception(failure); }
    st.bases = acgt;
    st.kmers = occurrences;
    for (const auto& s : spill) { std::error_code ec; st.spill_bytes += fs::file_size(s.path, ec); }
    log("spilled " + std::to_string(st.kmers) + " k-mer occurrences in " + std::to_string(st.spill_bytes) + " bytes");

    // Phase 2: workers take spill files one at a time, unpack, sort and count each partition in it, and
    // append its encoding to parts.bin.
    fs::create_directories(o.segment_dir);
    const std::string parts_path = (fs::path(o.segment_dir) / "parts.bin").string();
    std::FILE* out = std::fopen(parts_path.c_str(), "wb");
    if (out == nullptr) throw std::runtime_error("background: cannot write " + parts_path);
    std::vector<std::uint64_t> header(kHeaderWords + o.partitions, 0);
    header[0] = kMagic; header[1] = kVersion; header[2] = kK; header[3] = kM; header[4] = o.partitions;
    std::fwrite(header.data(), 8, header.size(), out);
    std::uint64_t out_words = header.size();
    std::vector<std::uint64_t> offsets(o.partitions, 0);
    std::mutex out_m;
    std::atomic<std::uint32_t> next_file{0};
    std::atomic<std::uint64_t> distinct{0};
    {
        std::vector<std::thread> workers;
        for (std::size_t t = 0; t < threads; ++t) {
            workers.emplace_back([&] {
                try {
                    const std::uint32_t per_file = o.partitions / files;
                    std::vector<std::vector<std::uint64_t>> codes(per_file);
                    std::string data, bases_s;
                    for (std::uint32_t b; (b = next_file++) < files;) {
                        {
                            std::ifstream in(spill[b].path, std::ios::binary);
                            data.assign(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
                        }
                        std::error_code ec;
                        fs::remove(spill[b].path, ec);
                        for (auto& v : codes) v.clear();
                        for (std::size_t at = 0; at + 6 <= data.size();) {
                            std::uint16_t p16; std::uint32_t len;
                            std::memcpy(&p16, data.data() + at, 2);
                            std::memcpy(&len, data.data() + at + 2, 4);
                            at += 6;
                            bases_s.resize(len);
                            for (std::uint32_t i = 0; i < len; ++i)
                                bases_s[i] = "ACGT"[(static_cast<unsigned char>(data[at + i / 4]) >> (2 * (i % 4))) & 3];
                            at += (len + 3) / 4;
                            auto& v = codes[p16 >> file_shift];
                            kmer::for_each_canonical_kmer(bases_s.data(), len, kK,
                                                          [&](std::uint64_t c, std::size_t, bool) { v.push_back(c); });
                        }
                        data.clear();
                        data.shrink_to_fit();
                        for (std::uint32_t j = 0; j < per_file; ++j) {
                            auto& v = codes[j];
                            std::sort(v.begin(), v.end());
                            std::vector<std::uint64_t> uniq, counts;
                            for (std::size_t i = 0; i < v.size();) {
                                std::size_t e = i + 1;
                                while (e < v.size() && v[e] == v[i]) ++e;
                                uniq.push_back(v[i]);
                                counts.push_back(e - i);
                                i = e;
                            }
                            const std::vector<std::uint64_t> blob = encode_partition(uniq, counts);
                            distinct += uniq.size();
                            const std::uint32_t p = (j << file_shift) | b;
                            std::lock_guard<std::mutex> lk(out_m);
                            offsets[p] = out_words;
                            if (std::fwrite(blob.data(), 8, blob.size(), out) != blob.size())
                                throw std::runtime_error("background: cannot write " + parts_path);
                            out_words += blob.size();
                            std::vector<std::uint64_t>().swap(v);
                        }
                    }
                } catch (...) {
                    std::lock_guard<std::mutex> lk(failure_m);
                    if (!failure) failure = std::current_exception();
                }
            });
        }
        for (auto& w : workers) w.join();
    }
    std::error_code ec;
    fs::remove_all(spill_dir, ec);
    if (failure) { std::fclose(out); std::rethrow_exception(failure); }
    header[5] = distinct;
    std::copy(offsets.begin(), offsets.end(), header.begin() + kHeaderWords);
    if (std::fseek(out, 0, SEEK_SET) != 0 || std::fwrite(header.data(), 8, header.size(), out) != header.size() ||
        std::fclose(out) != 0)
        throw std::runtime_error("background: cannot write " + parts_path);
    st.distinct = distinct;
    st.bytes = out_words * 8;
    return st;
}

Segment::Segment(const std::string& segment_dir) {
    const std::string path = (fs::path(segment_dir) / "parts.bin").string();
    const int fd = ::open(path.c_str(), O_RDONLY);
    if (fd < 0) throw std::runtime_error("background: cannot read " + path);
    struct stat sb {};
    if (::fstat(fd, &sb) != 0 || sb.st_size < static_cast<off_t>(kHeaderWords * 8)) {
        ::close(fd);
        throw std::runtime_error("background: " + path + " is not a background segment");
    }
    map_len_ = static_cast<std::size_t>(sb.st_size);
    map_ = ::mmap(nullptr, map_len_, PROT_READ, MAP_SHARED, fd, 0);
    ::close(fd);
    if (map_ == MAP_FAILED) { map_ = nullptr; throw std::runtime_error("background: cannot map " + path); }
    words_ = static_cast<const std::uint64_t*>(map_);
    if (words_[0] != kMagic || words_[1] != kVersion || words_[2] != kK || words_[3] != kM)
        throw std::runtime_error("background: " + path + " is not a background segment of this version");
    partitions_ = static_cast<std::uint32_t>(words_[4]);
    distinct_ = words_[5];
    offsets_ = words_ + kHeaderWords;
}

Segment::~Segment() {
    if (map_ != nullptr) ::munmap(map_, map_len_);
}

Segment::Part Segment::part(std::uint32_t p) const {
    Part pt;
    const std::uint64_t* w = words_ + offsets_[p];
    pt.codes = ef::View(w);
    const std::uint64_t n = pt.codes.size();
    pt.counts = reinterpret_cast<const std::uint8_t*>(w + pt.codes.words());
    const std::uint64_t* ov = w + pt.codes.words() + (n + 7) / 8;
    pt.overflow_n = ov[0];
    pt.overflow = ov + 1;
    return pt;
}

std::uint64_t Segment::count_at(const Part& pt, std::uint64_t index) const {
    const std::uint8_t c = pt.counts[index];
    if (c < 255) return c;
    std::uint64_t lo = 0, hi = pt.overflow_n;
    while (lo < hi) {
        const std::uint64_t mid = (lo + hi) / 2;
        if ((pt.overflow[mid] >> 32) < index) lo = mid + 1; else hi = mid;
    }
    if (lo < pt.overflow_n && (pt.overflow[lo] >> 32) == index) return pt.overflow[lo] & 0xFFFFFFFFULL;
    throw std::runtime_error("background: a count marked as overflow is missing from the overflow list");
}

std::uint64_t Segment::count(std::uint64_t canonical_kmer) const {
    const Part pt = part(partition_of(canonical_kmer, partitions_));
    const std::int64_t i = pt.codes.find(canonical_kmer);
    return i < 0 ? 0 : count_at(pt, static_cast<std::uint64_t>(i));
}

void Segment::for_each(std::uint32_t partition, const std::function<void(std::uint64_t, std::uint64_t)>& f) const {
    const Part pt = part(partition);
    pt.codes.for_each_value([&](std::uint64_t i, std::uint64_t code) { f(code, count_at(pt, i)); });
}

}  // namespace panvar::bg
