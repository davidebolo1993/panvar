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
#include <memory>
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

// Sorted codes with repeats -> distinct codes and their counts.
void collapse(const std::vector<std::uint64_t>& sorted, std::vector<std::uint64_t>& codes, std::vector<std::uint64_t>& counts) {
    codes.clear();
    counts.clear();
    for (std::size_t i = 0; i < sorted.size();) {
        std::size_t e = i + 1;
        while (e < sorted.size() && sorted[e] == sorted[i]) ++e;
        codes.push_back(sorted[i]);
        counts.push_back(e - i);
        i = e;
    }
}

// Writes a segment's parts.bin: a header with every partition's offset, then the partitions' encodings in
// whatever order threads finish them.
class PartsWriter {
public:
    PartsWriter(const std::string& segment_dir, std::uint32_t partitions)
        : path_((fs::path(segment_dir) / "parts.bin").string()), header_(kHeaderWords + partitions, 0), partitions_(partitions) {
        fs::create_directories(segment_dir);
        out_ = std::fopen(path_.c_str(), "wb");
        if (out_ == nullptr) throw std::runtime_error("background: cannot write " + path_);
        header_[0] = kMagic; header_[1] = kVersion; header_[2] = kK; header_[3] = kM; header_[4] = partitions;
        if (std::fwrite(header_.data(), 8, header_.size(), out_) != header_.size())
            throw std::runtime_error("background: cannot write " + path_);
        words_ = header_.size();
        written_.assign(partitions, 0);
    }
    ~PartsWriter() { if (out_ != nullptr) std::fclose(out_); }
    PartsWriter(const PartsWriter&) = delete;
    PartsWriter& operator=(const PartsWriter&) = delete;

    void put(std::uint32_t p, const std::vector<std::uint64_t>& codes, const std::vector<std::uint64_t>& counts) {
        const std::vector<std::uint64_t> blob = encode_partition(codes, counts);
        std::lock_guard<std::mutex> lk(m_);
        if (written_[p]) throw std::runtime_error("background: partition written twice");
        written_[p] = 1;
        header_[kHeaderWords + p] = words_;
        if (std::fwrite(blob.data(), 8, blob.size(), out_) != blob.size())
            throw std::runtime_error("background: cannot write " + path_);
        words_ += blob.size();
        distinct_ += codes.size();
    }
    // Writes the header; returns the file's size in bytes.
    std::uint64_t finish() {
        for (std::uint32_t p = 0; p < partitions_; ++p)
            if (!written_[p]) throw std::runtime_error("background: partition " + std::to_string(p) + " was not written");
        header_[5] = distinct_;
        if (std::fseek(out_, 0, SEEK_SET) != 0 || std::fwrite(header_.data(), 8, header_.size(), out_) != header_.size() ||
            std::fclose(out_) != 0) {
            out_ = nullptr;
            throw std::runtime_error("background: cannot write " + path_);
        }
        out_ = nullptr;
        return words_ * 8;
    }
    std::uint64_t distinct() const { return distinct_; }

private:
    std::string path_;
    std::vector<std::uint64_t> header_;
    std::uint32_t partitions_;
    std::FILE* out_ = nullptr;
    std::uint64_t words_ = 0, distinct_ = 0;
    std::vector<char> written_;
    std::mutex m_;
};

// Runs f(i) for i in [0, n) on `threads` threads; the first exception is rethrown after all stop.
template <class F>
void parallel_for(std::size_t n, std::size_t threads, F&& f) {
    std::atomic<std::size_t> next{0};
    std::exception_ptr failure;
    std::mutex fm;
    std::vector<std::thread> pool;
    for (std::size_t t = 0; t < std::max<std::size_t>(1, threads); ++t)
        pool.emplace_back([&] {
            try {
                for (std::size_t i; (i = next++) < n;) {
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
        o.partitions, o.spill_files != 0 ? o.spill_files : next_pow2((bases * 16 + per_thread - 1) / per_thread)));
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

    // Phase 2: threads take spill files one at a time, unpack, sort and count each partition in it, and
    // append its encoding to parts.bin.
    PartsWriter writer(o.segment_dir, o.partitions);
    try {
        parallel_for(files, threads, [&](std::size_t b) {
            const std::uint32_t per_file = o.partitions / files;
            std::vector<std::vector<std::uint64_t>> codes(per_file);
            std::string data, bases_s;
            {
                std::ifstream in(spill[b].path, std::ios::binary);
                data.assign(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
            }
            std::error_code ec;
            fs::remove(spill[b].path, ec);
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
                kmer::for_each_canonical_kmer(bases_s.data(), len, kK, [&](std::uint64_t c, std::size_t, bool) { v.push_back(c); });
            }
            std::string().swap(data);
            std::vector<std::uint64_t> uniq, counts;
            for (std::uint32_t j = 0; j < per_file; ++j) {
                auto& v = codes[j];
                std::sort(v.begin(), v.end());
                collapse(v, uniq, counts);
                std::vector<std::uint64_t>().swap(v);
                writer.put((j << file_shift) | static_cast<std::uint32_t>(b), uniq, counts);
            }
        });
    } catch (...) {
        std::error_code ec;
        fs::remove_all(spill_dir, ec);
        throw;
    }
    std::error_code ec;
    fs::remove_all(spill_dir, ec);
    st.bytes = writer.finish();
    st.distinct = writer.distinct();
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

namespace {

// One partition of several segments as distinct codes with summed counts, in code order.
void merged_partition(const std::vector<const Segment*>& segs, std::uint32_t p,
                      std::vector<std::uint64_t>& codes, std::vector<std::uint64_t>& counts) {
    std::vector<std::pair<std::uint64_t, std::uint64_t>> all;
    for (const Segment* s : segs) s->for_each(p, [&](std::uint64_t c, std::uint64_t n) { all.emplace_back(c, n); });
    if (segs.size() > 1) std::sort(all.begin(), all.end());
    codes.clear();
    counts.clear();
    for (const auto& [c, n] : all) {
        if (!codes.empty() && codes.back() == c) counts.back() += n;
        else { codes.push_back(c); counts.push_back(n); }
    }
}

}  // namespace

SegmentStats merge_segments(const MergeOptions& o) {
    if (o.segment_dirs.empty()) throw std::runtime_error("background: nothing to merge");
    std::vector<std::unique_ptr<Segment>> segs;
    std::vector<const Segment*> raw;
    for (const auto& d : o.segment_dirs) {
        segs.push_back(std::make_unique<Segment>(d));
        raw.push_back(segs.back().get());
        if (segs.back()->partitions() != segs.front()->partitions())
            throw std::runtime_error("background: segments with different partition counts cannot be merged");
    }
    const std::uint32_t partitions = segs.front()->partitions();
    PartsWriter writer(o.out_segment_dir, partitions);
    parallel_for(partitions, o.threads, [&](std::size_t p) {
        std::vector<std::uint64_t> codes, counts;
        merged_partition(raw, static_cast<std::uint32_t>(p), codes, counts);
        writer.put(static_cast<std::uint32_t>(p), codes, counts);
    });
    SegmentStats st;
    st.bytes = writer.finish();
    st.distinct = writer.distinct();
    return st;
}

std::vector<std::string> segment_dirs_of(const std::string& store_dir) {
    std::vector<std::pair<long, std::string>> found;
    std::error_code ec;
    for (const auto& e : fs::directory_iterator(store_dir, ec)) {
        const std::string name = e.path().filename().string();
        if (!e.is_directory() || name.rfind("seg-", 0) != 0) continue;
        const std::string num = name.substr(4);
        if (num.empty() || num.find_first_not_of("0123456789") != std::string::npos) continue;
        if (!fs::exists(e.path() / "parts.bin")) continue;
        found.emplace_back(std::stol(num), e.path().string());
    }
    std::sort(found.begin(), found.end());
    std::vector<std::string> out;
    for (auto& f : found) out.push_back(std::move(f.second));
    return out;
}

Store::Store(const std::string& store_dir) : dirs_(segment_dirs_of(store_dir)) {
    if (dirs_.empty()) throw std::runtime_error("background: " + store_dir + " holds no segment");
    for (const auto& d : dirs_) {
        segs_.push_back(std::make_unique<Segment>(d));
        if (segs_.back()->partitions() != segs_.front()->partitions())
            throw std::runtime_error("background: " + store_dir + " mixes segments with different partition counts");
    }
    partitions_ = segs_.front()->partitions();
}

std::uint64_t Store::count(std::uint64_t canonical_kmer) const {
    std::uint64_t n = 0;
    for (const auto& s : segs_) n += s->count(canonical_kmer);
    return n;
}

void Store::for_each(std::uint32_t partition, const std::function<void(std::uint64_t, std::uint64_t)>& f) const {
    std::vector<const Segment*> raw;
    for (const auto& s : segs_) raw.push_back(s.get());
    std::vector<std::uint64_t> codes, counts;
    merged_partition(raw, partition, codes, counts);
    for (std::size_t i = 0; i < codes.size(); ++i) f(codes[i], counts[i]);
}

}  // namespace panvar::bg
