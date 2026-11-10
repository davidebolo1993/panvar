#include "panvar/background.hpp"

#include "panvar/anchors.hpp"
#include "panvar/background_graph.hpp"
#include "panvar/gz_reader.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <memory>
#include <mutex>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <unistd.h>
#ifdef __linux__
#include <sched.h>
#endif

namespace panvar::bg {

namespace fs = std::filesystem;

// ---------------------------------------------------------------------------------------------------------
// Counts

namespace counts {

namespace {
std::uint32_t bit_length(std::uint64_t v) { return v == 0 ? 0 : 64 - static_cast<std::uint32_t>(__builtin_clzll(v)); }
std::uint64_t escape_of(std::uint32_t w) { return w >= 64 ? ~0ULL : (1ULL << w) - 1; }
}  // namespace

std::vector<std::uint64_t> encode(const std::vector<std::uint64_t>& c) {
    const std::uint64_t n = c.size();
    if (n == 0) return {0};
    std::uint64_t hist[65] = {}, ones[65] = {}, maxv = 0;
    for (const std::uint64_t v : c) {
        const std::uint32_t b = bit_length(v);
        ++hist[b];
        if ((v & (v + 1)) == 0) ++ones[b];   // 2^b - 1
        maxv = std::max(maxv, v);
    }
    const std::uint32_t full = std::max<std::uint32_t>(1, bit_length(maxv));
    std::uint32_t best_w = full;
    std::uint64_t best_cost = n * full;
    std::uint64_t above = 0;   // values with more than w bits
    for (std::uint32_t w = full - 1; w >= 1; --w) {
        above += hist[w + 1];
        const std::uint64_t exc = above + ones[w];
        const std::uint64_t cost = n * w + 128 * exc;
        if (cost < best_cost) { best_cost = cost; best_w = w; }
    }
    const std::uint32_t w = best_w;
    const bool exceptions = w < full;
    const std::uint64_t esc = escape_of(w);
    std::vector<std::uint64_t> idx, val;
    const std::uint64_t packed_words = (n * w + 63) / 64;
    std::vector<std::uint64_t> out(1 + packed_words, 0);
    std::uint64_t* packed = out.data() + 1;
    for (std::uint64_t i = 0; i < n; ++i) {
        std::uint64_t v = c[i];
        if (exceptions && v >= esc) { idx.push_back(i); val.push_back(v); v = esc; }
        const std::uint64_t bit = i * w;
        packed[bit / 64] |= v << (bit % 64);
        if (bit % 64 + w > 64) packed[bit / 64 + 1] |= v >> (64 - bit % 64);
    }
    out[0] = static_cast<std::uint64_t>(w) | (static_cast<std::uint64_t>(idx.size()) << 8);
    out.insert(out.end(), idx.begin(), idx.end());
    out.insert(out.end(), val.begin(), val.end());
    return out;
}

View::View(const std::uint64_t* w, std::uint64_t n) {
    width_ = static_cast<std::uint32_t>(w[0] & 0xff);
    n_exc_ = w[0] >> 8;
    if (n > 0 && (width_ == 0 || width_ > 64)) throw std::runtime_error("background: corrupt count block");
    packed_ = w + 1;
    const std::uint64_t pw = (n * width_ + 63) / 64;
    exc_index_ = packed_ + pw;
    exc_count_ = exc_index_ + n_exc_;
    words_ = 1 + pw + 2 * n_exc_;
    escape_ = escape_of(width_);
}

std::uint64_t View::raw(std::uint64_t i) const {
    const std::uint64_t bit = i * width_;
    std::uint64_t v = packed_[bit / 64] >> (bit % 64);
    if (bit % 64 + width_ > 64) v |= packed_[bit / 64 + 1] << (64 - bit % 64);
    return v & escape_;
}

std::uint64_t View::get(std::uint64_t i) const {
    const std::uint64_t v = raw(i);
    if (n_exc_ == 0 || v != escape_) return v;
    const std::uint64_t* it = std::lower_bound(exc_index_, exc_index_ + n_exc_, i);
    if (it == exc_index_ + n_exc_ || *it != i) throw std::runtime_error("background: a count marked as exception is missing");
    return exc_count_[it - exc_index_];
}

std::uint64_t View::Cursor::next() {
    const std::uint64_t v = v_->raw(i_++);
    if (v_->n_exc_ == 0 || v != v_->escape_) return v;
    return v_->exc_count_[e_++];
}

}  // namespace counts

// ---------------------------------------------------------------------------------------------------------
// Utilities

void radix_sort(std::uint64_t* a, std::size_t n) {
    struct Rec {
        static void sort(std::uint64_t* a, std::size_t n, int shift) {
            if (n < 256 || shift < 0) { std::sort(a, a + n); return; }
            std::size_t count[256] = {};
            for (std::size_t i = 0; i < n; ++i) ++count[(a[i] >> shift) & 255];
            for (std::size_t b = 0; b < 256; ++b)
                if (count[b] == n) { sort(a, n, shift - 8); return; }   // one bucket: next digit
            std::size_t head[256], tail[256], at = 0;
            for (std::size_t b = 0; b < 256; ++b) { head[b] = at; at += count[b]; tail[b] = at; }
            for (std::size_t b = 0; b < 256; ++b) {
                while (head[b] < tail[b]) {
                    std::uint64_t v = a[head[b]];
                    std::size_t d = (v >> shift) & 255;
                    while (d != b) {   // carry v to its bucket, take what was there
                        std::swap(v, a[head[d]++]);
                        d = (v >> shift) & 255;
                    }
                    a[head[b]++] = v;
                }
            }
            std::size_t start = 0;
            for (std::size_t b = 0; b < 256; ++b) {
                if (count[b] > 1) sort(a + start, count[b], shift - 8);
                start += count[b];
            }
        }
    };
    Rec::sort(a, n, 56);
}

std::size_t default_threads() {
    if (const char* s = std::getenv("SLURM_CPUS_PER_TASK")) {
        const long v = std::strtol(s, nullptr, 10);
        if (v > 0) return static_cast<std::size_t>(v);
    }
#ifdef __linux__
    cpu_set_t set;
    CPU_ZERO(&set);
    if (sched_getaffinity(0, sizeof set, &set) == 0) {
        const int c = CPU_COUNT(&set);
        if (c > 0) return static_cast<std::size_t>(c);
    }
#endif
    return std::max(1u, std::thread::hardware_concurrency());
}

namespace {

constexpr std::uint64_t kPartsMagic = 0x3254525047425650ULL;   // "PVBGPRT2", little-endian
constexpr std::uint64_t kVersion = 2;
constexpr std::size_t kHeaderWords = 9;   // magic, version, k, m, partitions, first, end, distinct, words
constexpr std::size_t kChunkBases = std::size_t{8} << 20;

bool little_endian() { const std::uint16_t x = 1; return *reinterpret_cast<const unsigned char*>(&x) == 1; }

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

void raise_open_file_limit(std::size_t files) {
    struct rlimit rl {};
    if (getrlimit(RLIMIT_NOFILE, &rl) != 0) return;
    const rlim_t want = static_cast<rlim_t>(files + 64);
    if (rl.rlim_cur >= want) return;
    rl.rlim_cur = rl.rlim_max == RLIM_INFINITY ? want : std::min<rlim_t>(want, rl.rlim_max);
    setrlimit(RLIMIT_NOFILE, &rl);
}

// A FILE* closed on scope exit.
struct File {
    std::FILE* fp = nullptr;
    File() = default;
    File(const std::string& path, const char* mode) : fp(std::fopen(path.c_str(), mode)) {}
    ~File() { if (fp != nullptr) std::fclose(fp); }
    File(const File&) = delete;
    File& operator=(const File&) = delete;
    bool close() { const bool ok = fp == nullptr || std::fclose(fp) == 0; fp = nullptr; return ok; }
};

void fsync_path(const std::string& path) {
    const int fd = ::open(path.c_str(), O_RDONLY);
    if (fd >= 0) { ::fsync(fd); ::close(fd); }
}

// The bases a genome holds, estimated from its .fai or its file size.
std::uint64_t estimated_bases(const std::string& path) {
    std::ifstream fai(path + ".fai");
    if (fai) {
        std::uint64_t total = 0;
        std::string line;
        while (std::getline(fai, line)) {
            const auto t1 = line.find('\t');
            if (t1 != std::string::npos) total += std::strtoull(line.c_str() + t1 + 1, nullptr, 10);
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

std::uint64_t hash_bytes(const char* s, std::size_t n, std::uint64_t h) {
    std::size_t i = 0;
    for (; i + 8 <= n; i += 8) { std::uint64_t w; std::memcpy(&w, s + i, 8); h = kmer::mix64(h ^ w); }
    std::uint64_t w = 0;
    std::memcpy(&w, s + i, n - i);
    return kmer::mix64(h ^ w ^ (static_cast<std::uint64_t>(n) << 56));
}

std::string hex64(std::uint64_t v) {
    char buf[17];
    std::snprintf(buf, sizeof buf, "%016llx", static_cast<unsigned long long>(v));
    return buf;
}

// Spill records: varint partition-in-file, varint length, varint multiplicity, then the bases, 4 per
// byte, base i in bits 2*(i%4).
void put_varint(std::string& buf, std::uint64_t v) {
    while (v >= 0x80) { buf += static_cast<char>((v & 0x7f) | 0x80); v >>= 7; }
    buf += static_cast<char>(v);
}
bool get_varint(const unsigned char*& p, const unsigned char* end, std::uint64_t& v) {
    v = 0;
    for (int shift = 0; shift < 64; shift += 7) {
        if (p >= end) return false;
        const unsigned char b = *p++;
        v |= static_cast<std::uint64_t>(b & 0x7f) << shift;
        if (!(b & 0x80)) return true;
    }
    throw std::runtime_error("background: corrupt spill record");
}
void pack_bases(std::string& buf, const char* s, std::size_t len) {
    const std::size_t at = buf.size();
    buf.resize(at + (len + 3) / 4, '\0');
    for (std::size_t i = 0; i < len; ++i)
        buf[at + i / 4] = static_cast<char>(static_cast<unsigned char>(buf[at + i / 4]) | (kmer::base_code(s[i]) << (2 * (i % 4))));
}

// Calls f(canonical code) for every 31-mer of packed bases.
template <class F>
void for_each_packed_kmer(const unsigned char* packed, std::size_t len, F&& f) {
    const std::uint64_t mask = kmer::mask_of(kK);
    const std::size_t rc_shift = 2 * (kK - 1);
    std::uint64_t fwd = 0, rev = 0;
    for (std::size_t i = 0; i < len; ++i) {
        const std::uint64_t b = (packed[i / 4] >> (2 * (i % 4))) & 3;
        fwd = ((fwd << 2) | b) & mask;
        rev = (rev >> 2) | ((3 - b) << rc_shift);
        if (i + 1 >= kK) f(fwd <= rev ? fwd : rev);
    }
}

// A unit of work for the counting threads: a stretch of sequence, how many times it counts, and where it
// sits (for anchors; `anchors` false for graph pieces, whose anchors are placed per haplotype elsewhere).
struct Item {
    std::string seq;
    std::uint64_t mult = 1;
    std::uint32_t genome = 0, contig = 0;
    std::uint64_t offset = 0;
    bool anchors = true;
};

// A bounded queue between the readers and the workers; abort() releases everyone when a thread fails.
class ItemQueue {
public:
    explicit ItemQueue(std::size_t capacity) : capacity_(capacity) {}
    bool push(Item item) {
        std::unique_lock<std::mutex> lk(m_);
        not_full_.wait(lk, [&] { return q_.size() < capacity_ || aborted_; });
        if (aborted_) return false;
        q_.push_back(std::move(item));
        not_empty_.notify_one();
        return true;
    }
    bool pop(Item& item) {
        std::unique_lock<std::mutex> lk(m_);
        not_empty_.wait(lk, [&] { return !q_.empty() || closed_ || aborted_; });
        if (aborted_ || q_.empty()) return false;
        item = std::move(q_.front());
        q_.pop_front();
        not_full_.notify_one();
        return true;
    }
    void close() { std::lock_guard<std::mutex> lk(m_); closed_ = true; not_empty_.notify_all(); }
    void abort() { std::lock_guard<std::mutex> lk(m_); aborted_ = true; not_empty_.notify_all(); not_full_.notify_all(); }
private:
    std::size_t capacity_;
    std::deque<Item> q_;
    bool closed_ = false, aborted_ = false;
    std::mutex m_;
    std::condition_variable not_empty_, not_full_;
};

// Identical super-k-mers of one spill file, stored once with their multiplicity, in a table of bounded
// size; when full, written out and emptied.
class DedupTable {
public:
    // The records take three quarters of the budget, reserved once (pages are touched only as they fill), the
    // slots the rest; nothing grows past that, so the table holds its share exactly.
    explicit DedupTable(std::size_t budget) : budget_(budget) {
        slots_.assign(1024, 0);
        arena_.reserve(budget * 3 / 4);
        max_slots_ = 1024;
        while (max_slots_ * 2 * 8 <= budget / 4) max_slots_ *= 2;
    }
    // `rec` is one spill record (as put_varint/pack_bases write it); appends to `out` when the table flushes.
    void add(std::uint64_t j, std::uint64_t len, std::uint64_t mult, const unsigned char* packed, std::string& out) {
        const std::size_t pbytes = (len + 3) / 4;
        std::uint64_t h = hash_bytes(reinterpret_cast<const char*>(packed), pbytes, j * 0x9e3779b97f4a7c15ULL ^ len);
        const std::uint64_t fp = (h >> 40) | 1;   // never 0
        std::uint64_t mask = slots_.size() - 1;
        for (std::uint64_t s = h & mask;; s = (s + 1) & mask) {
            const std::uint64_t slot = slots_[s];
            if (slot == 0) {
                const std::uint64_t off = arena_.size();
                if (off >= (1ULL << 40) || off + 12 + pbytes > arena_.capacity()) {
                    if (off == 0) { emit(out, j, len, mult, packed); return; }   // larger than the whole table
                    flush(out);
                    add(j, len, mult, packed, out);
                    return;
                }
                const std::uint32_t m32 = static_cast<std::uint32_t>(std::min<std::uint64_t>(mult, 0xFFFFFFFFULL));
                arena_.append(reinterpret_cast<const char*>(&m32), 4);
                const std::uint32_t j32 = static_cast<std::uint32_t>(j), l32 = static_cast<std::uint32_t>(len);
                arena_.append(reinterpret_cast<const char*>(&j32), 4);
                arena_.append(reinterpret_cast<const char*>(&l32), 4);
                arena_.append(reinterpret_cast<const char*>(packed), pbytes);
                slots_[s] = (fp << 40) | off;
                if (mult > 0xFFFFFFFFULL) big_.emplace_back(off, mult);
                if (++used_ * 10 > slots_.size() * 7) grow(out);
                return;
            }
            if ((slot >> 40) != fp) continue;
            const std::uint64_t off = slot & ((1ULL << 40) - 1);
            std::uint32_t ej, el;
            std::memcpy(&ej, arena_.data() + off + 4, 4);
            std::memcpy(&el, arena_.data() + off + 8, 4);
            if (ej != j || el != len || std::memcmp(arena_.data() + off + 12, packed, pbytes) != 0) continue;
            std::uint32_t m32;
            std::memcpy(&m32, arena_.data() + off, 4);
            if (static_cast<std::uint64_t>(m32) + mult > 0xFFFFFFFFULL) {   // keep exact: write it out separately
                emit(out, j, len, mult, packed);
                return;
            }
            m32 += static_cast<std::uint32_t>(mult);
            std::memcpy(&arena_[off], &m32, 4);
            ++hits_;
            return;
        }
    }
    void flush(std::string& out) {
        for (std::size_t off = 0; off < arena_.size();) {
            std::uint32_t m32, j32, l32;
            std::memcpy(&m32, arena_.data() + off, 4);
            std::memcpy(&j32, arena_.data() + off + 4, 4);
            std::memcpy(&l32, arena_.data() + off + 8, 4);
            std::uint64_t mult = m32;
            for (const auto& [bo, bm] : big_) if (bo == off) mult = bm;
            emit(out, j32, l32, mult, reinterpret_cast<const unsigned char*>(arena_.data() + off + 12));
            off += 12 + (l32 + 3) / 4;
        }
        arena_.clear();
        big_.clear();
        std::fill(slots_.begin(), slots_.end(), 0);
        used_ = 0;
    }
    std::uint64_t hits() const { return hits_; }
private:
    static void emit(std::string& out, std::uint64_t j, std::uint64_t len, std::uint64_t mult, const unsigned char* packed) {
        put_varint(out, j); put_varint(out, len); put_varint(out, mult);
        out.append(reinterpret_cast<const char*>(packed), (len + 3) / 4);
    }
    void grow(std::string& out) {
        if (slots_.size() >= max_slots_) { flush(out); return; }
        std::vector<std::uint64_t> old;
        old.swap(slots_);
        slots_.assign(old.size() * 2, 0);
        const std::uint64_t mask = slots_.size() - 1;
        for (const std::uint64_t slot : old) {
            if (slot == 0) continue;
            const std::uint64_t off = slot & ((1ULL << 40) - 1);
            std::uint32_t j32, l32;
            std::memcpy(&j32, arena_.data() + off + 4, 4);
            std::memcpy(&l32, arena_.data() + off + 8, 4);
            const std::uint64_t h = hash_bytes(arena_.data() + off + 12, (l32 + 3) / 4, j32 * 0x9e3779b97f4a7c15ULL ^ l32);
            std::uint64_t s = h & mask;
            while (slots_[s] != 0) s = (s + 1) & mask;
            slots_[s] = slot;
        }
    }
    std::size_t budget_;
    std::string arena_;
    std::vector<std::uint64_t> slots_;
    std::vector<std::pair<std::uint64_t, std::uint64_t>> big_;
    std::size_t used_ = 0, max_slots_ = 1024;
    std::uint64_t hits_ = 0;
};

// Writes one parts file: a header (magic written last), per partition its (offset, size) in words, then
// the partitions' encodings in whatever order threads finish them. Written to <path>.tmp and renamed.
class PartsWriter {
public:
    PartsWriter(const std::string& path, std::uint32_t partitions, std::uint32_t first, std::uint32_t end)
        : path_(path), tmp_(path + ".tmp"), first_(first), end_(end), table_(2 * (end - first), 0), written_(end - first, 0) {
        out_.fp = std::fopen(tmp_.c_str(), "wb");
        if (out_.fp == nullptr) throw std::runtime_error("background: cannot write " + tmp_);
        header_ = {0, kVersion, kK, kM, partitions, first, end, 0, 0};
        words_ = kHeaderWords + table_.size();
        std::vector<std::uint64_t> zeros(words_, 0);
        if (std::fwrite(zeros.data(), 8, zeros.size(), out_.fp) != zeros.size()) fail();
    }
    ~PartsWriter() { if (out_.fp != nullptr) { out_.close(); std::error_code ec; fs::remove(tmp_, ec); } }
    PartsWriter(const PartsWriter&) = delete;
    PartsWriter& operator=(const PartsWriter&) = delete;

    void put(std::uint32_t p, const std::vector<std::uint64_t>& codes, const std::vector<std::uint64_t>& cnt) {
        std::vector<std::uint64_t> blob = ef::build(codes, 2 * kK);
        const std::vector<std::uint64_t> cw = counts::encode(cnt);
        blob.insert(blob.end(), cw.begin(), cw.end());
        std::lock_guard<std::mutex> lk(m_);
        const std::uint32_t i = p - first_;
        if (p < first_ || p >= end_ || written_[i]) throw std::runtime_error("background: partition written twice or out of range");
        written_[i] = 1;
        table_[2 * i] = words_;
        table_[2 * i + 1] = blob.size();
        if (std::fwrite(blob.data(), 8, blob.size(), out_.fp) != blob.size()) fail();
        words_ += blob.size();
        distinct_ += codes.size();
    }
    // Completes the file (fsync, magic last) and renames it into place; returns its size in bytes.
    std::uint64_t finish() {
        for (std::uint32_t i = 0; i < end_ - first_; ++i)
            if (!written_[i]) throw std::runtime_error("background: partition " + std::to_string(first_ + i) + " was not written");
        header_[7] = distinct_;
        header_[8] = words_;
        if (std::fseek(out_.fp, 0, SEEK_SET) != 0 || std::fwrite(header_.data(), 8, header_.size(), out_.fp) != header_.size() ||
            std::fwrite(table_.data(), 8, table_.size(), out_.fp) != table_.size() || std::fflush(out_.fp) != 0 ||
            ::fsync(fileno(out_.fp)) != 0)
            fail();
        const std::uint64_t magic = kPartsMagic;
        if (std::fseek(out_.fp, 0, SEEK_SET) != 0 || std::fwrite(&magic, 8, 1, out_.fp) != 1 || std::fflush(out_.fp) != 0 ||
            ::fsync(fileno(out_.fp)) != 0 || !out_.close())
            fail();
        fs::rename(tmp_, path_);
        return words_ * 8;
    }
    std::uint64_t distinct() const { return distinct_; }

private:
    [[noreturn]] void fail() { throw std::runtime_error("background: cannot write " + tmp_); }
    std::string path_, tmp_;
    std::uint32_t first_, end_;
    std::vector<std::uint64_t> header_, table_;
    std::vector<char> written_;
    File out_;
    std::uint64_t words_ = 0, distinct_ = 0;
    std::mutex m_;
};

// Runs f(i) for i in [0, n) on `threads` threads; the first exception is rethrown after all stop.
template <class F>
void parallel_for(std::size_t n, std::size_t threads, F&& f) {
    std::atomic<std::size_t> next{0};
    std::atomic<bool> failed{false};
    std::exception_ptr failure;
    std::mutex fm;
    std::vector<std::thread> pool;
    for (std::size_t t = 0; t < std::max<std::size_t>(1, std::min(threads, n)); ++t)
        pool.emplace_back([&] {
            try {
                for (std::size_t i; !failed && (i = next++) < n;) f(i);
            } catch (...) {
                std::lock_guard<std::mutex> lk(fm);
                if (!failure) failure = std::current_exception();
                failed = true;
            }
        });
    for (auto& t : pool) t.join();
    if (failure) std::rethrow_exception(failure);
}

// Merges sorted (code, count) lists into one, summing equal codes.
void merge_runs(const std::vector<std::pair<const std::uint64_t*, const std::uint64_t*>>& runs, const std::vector<std::size_t>& lens,
                std::vector<std::uint64_t>& codes, std::vector<std::uint64_t>& cnt) {
    codes.clear();
    cnt.clear();
    using E = std::pair<std::uint64_t, std::size_t>;
    std::priority_queue<E, std::vector<E>, std::greater<E>> heap;
    std::vector<std::size_t> at(runs.size(), 0);
    for (std::size_t r = 0; r < runs.size(); ++r) if (lens[r] > 0) heap.emplace(runs[r].first[0], r);
    while (!heap.empty()) {
        const auto [code, r] = heap.top();
        heap.pop();
        const std::uint64_t c = runs[r].second[at[r]];
        if (!codes.empty() && codes.back() == code) cnt.back() += c;
        else { codes.push_back(code); cnt.push_back(c); }
        if (++at[r] < lens[r]) heap.emplace(runs[r].first[at[r]], r);
    }
}

// One spill file's partitions counted under a memory budget: occurrences gathered per partition, sorted
// and collapsed into a run when the budget fills, runs written to a run file when they alone exceed half
// of it, and every partition's runs merged at the end.
class PartitionCounter {
public:
    PartitionCounter(std::uint32_t per_file, std::size_t budget, std::string run_path)
        : acc_(per_file), budget_(budget), run_path_(std::move(run_path)) {}
    ~PartitionCounter() { if (runs_.fp != nullptr) { runs_.close(); } std::error_code ec; fs::remove(run_path_, ec); }

    void add(std::uint32_t j, std::uint64_t code, std::uint64_t mult) {
        Acc& a = acc_[j];
        // memory is counted by capacity: a vector that grows briefly holds its old and new buffers
        if (mult == 1) {
            const std::size_t before = a.raw.capacity();
            a.raw.push_back(code);
            if (a.raw.capacity() != before) used_ += (a.raw.capacity() + before) * 8;
        } else {
            const std::size_t before = a.wt.capacity();
            a.wt.emplace_back(code, mult);
            if (a.wt.capacity() != before) used_ += (a.wt.capacity() + before) * 16;
        }
        if (used_ > budget_) compact(true);
    }
    // The partition's distinct codes and counts; call once per partition after all adds.
    void finish(std::uint32_t j, std::vector<std::uint64_t>& codes, std::vector<std::uint64_t>& cnt) {
        if (!finished_) { compact(true); finished_ = true; if (runs_.fp != nullptr && std::fflush(runs_.fp) != 0) fail(); }
        Acc& a = acc_[j];
        if (a.disk.empty()) { codes.swap(a.rc); cnt.swap(a.rn); a = Acc(); return; }
        // the runs on disk are read in small chunks while they are merged, so a partition's runs need not fit
        struct Cursor {
            std::uint64_t off = 0, n = 0, at = 0, buf_at = 0;
            std::vector<std::uint64_t> c, k;
        };
        std::vector<Cursor> cur(a.disk.size());
        constexpr std::uint64_t kChunk = 1 << 15;
        auto refill = [&](Cursor& r) {
            const std::uint64_t take = std::min<std::uint64_t>(kChunk, r.n - r.at);
            r.c.resize(take);
            r.k.resize(take);
            if (std::fseek(runs_.fp, static_cast<long>(r.off + r.at * 8), SEEK_SET) != 0 || std::fread(r.c.data(), 8, take, runs_.fp) != take ||
                std::fseek(runs_.fp, static_cast<long>(r.off + (r.n + r.at) * 8), SEEK_SET) != 0 || std::fread(r.k.data(), 8, take, runs_.fp) != take)
                throw std::runtime_error("background: cannot read " + run_path_);
            r.buf_at = r.at;
        };
        using E = std::pair<std::uint64_t, std::size_t>;
        std::priority_queue<E, std::vector<E>, std::greater<E>> heap;
        for (std::size_t r = 0; r < cur.size(); ++r) {
            cur[r].off = a.disk[r].first;
            cur[r].n = a.disk[r].second;
            refill(cur[r]);
            heap.emplace(cur[r].c[0], r);
        }
        const std::size_t mem = cur.size();   // the run still in memory takes the last slot
        std::size_t mi = 0;
        if (!a.rc.empty()) heap.emplace(a.rc[0], mem);
        codes.clear();
        cnt.clear();
        while (!heap.empty()) {
            const auto [code, r] = heap.top();
            heap.pop();
            std::uint64_t c;
            if (r == mem) {
                c = a.rn[mi++];
                if (mi < a.rc.size()) heap.emplace(a.rc[mi], mem);
            } else {
                Cursor& x = cur[r];
                c = x.k[x.at - x.buf_at];
                if (++x.at < x.n) {
                    if (x.at - x.buf_at >= x.c.size()) refill(x);
                    heap.emplace(x.c[x.at - x.buf_at], r);
                }
            }
            if (!codes.empty() && codes.back() == code) cnt.back() += c;
            else { codes.push_back(code); cnt.push_back(c); }
        }
        a = Acc();
    }
    std::uint64_t run_bytes() const { return run_bytes_; }

private:
    struct Acc {
        std::vector<std::uint64_t> raw;
        std::vector<std::pair<std::uint64_t, std::uint64_t>> wt;
        std::vector<std::uint64_t> rc, rn;                         // collapsed run in memory
        std::vector<std::pair<std::uint64_t, std::uint64_t>> disk; // (offset, n) of runs on disk
    };
    [[noreturn]] void fail() { throw std::runtime_error("background: cannot write " + run_path_); }

    // Every partition's new occurrences sorted, collapsed and merged into its run; when what is held (runs
    // so far plus occurrences not yet collapsed) passes the budget, a run goes to the run file at once.
    void compact(bool may_spill) {
        std::size_t raw_left = 0, held = 0;
        for (const Acc& a : acc_) {
            raw_left += a.raw.capacity() * 8 + a.wt.capacity() * 16;
            held += (a.rc.capacity() + a.rn.capacity()) * 8;
        }
        std::vector<std::uint64_t> c1, n1, c2, n2;
        for (Acc& a : acc_) {
            if (!a.raw.empty() || !a.wt.empty()) {
                raw_left -= a.raw.capacity() * 8 + a.wt.capacity() * 16;
                held -= (a.rc.capacity() + a.rn.capacity()) * 8;
                radix_sort(a.raw.data(), a.raw.size());
                c1.clear(); n1.clear();
                for (std::size_t i = 0; i < a.raw.size();) {
                    std::size_t e = i + 1;
                    while (e < a.raw.size() && a.raw[e] == a.raw[i]) ++e;
                    c1.push_back(a.raw[i]); n1.push_back(e - i);
                    i = e;
                }
                std::vector<std::uint64_t>().swap(a.raw);
                std::sort(a.wt.begin(), a.wt.end());
                c2.clear(); n2.clear();
                for (const auto& [c, n] : a.wt) {
                    if (!c2.empty() && c2.back() == c) n2.back() += n;
                    else { c2.push_back(c); n2.push_back(n); }
                }
                std::vector<std::pair<std::uint64_t, std::uint64_t>>().swap(a.wt);
                std::vector<std::uint64_t> rc, rn;
                rc.reserve(c1.size() + c2.size() + a.rc.size());
                rn.reserve(rc.capacity());
                merge_runs({{c1.data(), n1.data()}, {c2.data(), n2.data()}, {a.rc.data(), a.rn.data()}}, {c1.size(), c2.size(), a.rc.size()}, rc, rn);
                rc.shrink_to_fit(); rn.shrink_to_fit();
                a.rc.swap(rc); a.rn.swap(rn);
                held += (a.rc.capacity() + a.rn.capacity()) * 8;
            }
            if (may_spill && !a.rc.empty() && held + raw_left > budget_ / 2) {
                held -= (a.rc.capacity() + a.rn.capacity()) * 8;
                spill(a);
            }
        }
        c1 = {}; n1 = {}; c2 = {}; n2 = {};
        used_ = held;
    }

    void spill(Acc& a) {
        if (runs_.fp == nullptr) {
            runs_.fp = std::fopen(run_path_.c_str(), "w+b");
            if (runs_.fp == nullptr) fail();
        }
        if (std::fseek(runs_.fp, 0, SEEK_END) != 0) fail();
        const long off = std::ftell(runs_.fp);
        if (std::fwrite(a.rc.data(), 8, a.rc.size(), runs_.fp) != a.rc.size() || std::fwrite(a.rn.data(), 8, a.rn.size(), runs_.fp) != a.rn.size())
            fail();
        a.disk.emplace_back(static_cast<std::uint64_t>(off), a.rc.size());
        run_bytes_ += a.rc.size() * 16;
        std::vector<std::uint64_t>().swap(a.rc);
        std::vector<std::uint64_t>().swap(a.rn);
    }

    std::vector<Acc> acc_;
    std::size_t budget_, used_ = 0;
    std::string run_path_;
    File runs_;
    std::uint64_t run_bytes_ = 0;
    bool finished_ = false;
};

enum class InputKind { Fasta, Gfa };

InputKind detect_input(const std::string& path) {
    GzLineReader in(path);
    if (!in.ok()) throw std::runtime_error("background: cannot read " + path);
    std::string line;
    while (in.getline(line)) {
        if (line.empty()) continue;
        if (line[0] == '>') return InputKind::Fasta;
        if (line[0] == '@') throw std::runtime_error("background: " + path + " looks like FASTQ; give genome assemblies (FASTA) or graphs (GFA)");
        if (line.size() > 1 && line[1] == '\t' && std::string("HSLPWCJ#").find(line[0]) != std::string::npos) return InputKind::Gfa;
        throw std::runtime_error("background: " + path + " is neither FASTA nor GFA");
    }
    return InputKind::Fasta;   // empty: a genome with no sequence
}

}  // namespace

// ---------------------------------------------------------------------------------------------------------
// Segment info

SegmentInfo read_segment_info(const std::string& dir) {
    const std::string path = (fs::path(dir) / "segment.tsv").string();
    std::ifstream f(path);
    if (!f) throw std::runtime_error("background: cannot read " + path);
    SegmentInfo info;
    std::string line;
    bool format_ok = false;
    while (std::getline(f, line)) {
        std::vector<std::string> t;
        std::istringstream ss(line);
        for (std::string x; std::getline(ss, x, '\t');) t.push_back(x);
        if (t.empty()) continue;
        if (t[0] == "format" && t.size() == 3 && t[1] == "panvar-background-segment" && t[2] == std::to_string(kVersion)) format_ok = true;
        else if (t[0] == "partitions" && t.size() == 2) info.partitions = static_cast<std::uint32_t>(std::stoul(t[1]));
        else if (t[0] == "parts" && t.size() == 4) info.parts.push_back({t[1], static_cast<std::uint32_t>(std::stoul(t[2])), static_cast<std::uint32_t>(std::stoul(t[3]))});
        else if (t[0] == "genome" && t.size() == 5) info.genomes.push_back({t[1], t[2], std::stoull(t[3]), t[4]});
        else if (t[0] == "distinct" && t.size() == 2) info.distinct = std::stoull(t[1]);
        else if (t[0] == "occurrences" && t.size() == 2) info.occurrences = std::stoull(t[1]);
        else if (t[0] == "anchors" && t.size() == 3) info.anchor_sets.emplace_back(t[1], std::stoull(t[2], nullptr, 16));
    }
    if (!format_ok) throw std::runtime_error("background: " + path + " is not a segment of this version");
    return info;
}

void write_segment_info(const std::string& dir, const SegmentInfo& info) {
    const std::string path = (fs::path(dir) / "segment.tsv").string(), tmp = path + ".tmp";
    std::ostringstream f;
    f << "format\tpanvar-background-segment\t" << kVersion << "\npartitions\t" << info.partitions << '\n';
    for (const auto& p : info.parts) f << "parts\t" << p.file << '\t' << p.first << '\t' << p.end << '\n';
    for (const auto& g : info.genomes) f << "genome\t" << g.name << '\t' << g.path << '\t' << g.bytes << '\t' << g.checksum << '\n';
    f << "distinct\t" << info.distinct << "\noccurrences\t" << info.occurrences << '\n';
    for (const auto& [name, fp] : info.anchor_sets) f << "anchors\t" << name << '\t' << hex64(fp) << '\n';
    write_file_durably(tmp, f.str());
    fs::rename(tmp, path);
}

// ---------------------------------------------------------------------------------------------------------
// Building a segment

SegmentStats build_segment(const BuildOptions& o) {
    if (!little_endian()) throw std::runtime_error("background: stores are little-endian; this machine is not");
    if (o.partitions == 0 || o.partitions > kMaxPartitions || (o.partitions & (o.partitions - 1)) != 0)
        throw std::runtime_error("background: the partition count must be a power of two, at most 65536");
    if (o.genomes.empty()) throw std::runtime_error("background: no genome given");
    if (o.genomes.size() > 65535) throw std::runtime_error("background: at most 65535 genomes in one segment");
    if (fs::exists(o.segment_dir)) {   // only a folder its caller created and claimed (holding job.lock) may exist
        for (const auto& e : fs::directory_iterator(o.segment_dir))
            if (e.path().filename() != "job.lock") throw std::runtime_error("background: " + o.segment_dir + " is not empty");
    }
    const std::size_t threads = std::max<std::size_t>(1, o.threads);
    auto log = [&](const std::string& m) { if (o.log) o.log(m); };

    std::vector<InputKind> kinds;
    std::uint64_t bases = 0;
    for (const Genome& g : o.genomes) {
        kinds.push_back(detect_input(g.path));
        bases += estimated_bases(g.path);
    }

    // Spill files: as many as it takes for one to be counted in a thread's share of the memory (8 bytes per
    // occurrence and room to sort; a file that still does not fit is counted in sorted runs on disk).
    const std::uint64_t per_thread = std::max<std::uint64_t>(std::uint64_t{64} << 20, o.memory_bytes / threads);
    // What the bound must also hold: the reader's chunks in flight, and the anchors found (20 bytes each; about
    // one per reference anchor and genome).
    const std::size_t chunk_bases = std::clamp<std::size_t>(o.memory_bytes / (16 * threads), std::size_t{1} << 20, kChunkBases);
    std::uint64_t anchor_bytes = 0;
    for (const AnchorSet& s : o.anchor_sets) anchor_bytes += s.ref->code.size() * 20ULL * o.genomes.size();
    const std::uint64_t fixed = chunk_bases * (3 * threads + 2) + anchor_bytes;
    const std::uint64_t left = o.memory_bytes > fixed ? o.memory_bytes - fixed : 0;
    // a partition is held whole while it is encoded (about 40 bytes per distinct k-mer)
    if (bases / o.partitions * 40 > o.memory_bytes / threads)
        log("warning: with " + std::to_string(o.partitions) + " partitions, one may need about " +
            std::to_string(bases / o.partitions * 40 >> 20) + " MB to encode, more than a thread's share of --memory; use more partitions");
    if (left < o.memory_bytes / 4)
        log("warning: --memory " + std::to_string(o.memory_bytes >> 20) + " MB leaves little for counting after the anchors (" +
            std::to_string(anchor_bytes >> 20) + " MB) and reading buffers");
    // a thread's share, less room for the copies made while runs are collapsed and merged
    const std::uint64_t count_budget = std::max<std::uint64_t>(std::uint64_t{1} << 20, std::max<std::uint64_t>(left, o.memory_bytes / 4) / threads * 5 / 8);
    if (o.spill_files != 0 && (o.spill_files & (o.spill_files - 1)) != 0)
        throw std::runtime_error("background: the spill file count must be a power of two");
    std::uint64_t want_files = o.spill_files != 0 ? o.spill_files : std::max<std::uint64_t>(1, next_pow2((bases * 12 + per_thread - 1) / per_thread));
    // every thread keeps a buffer of at least 4 KB per file: an eighth of the bound caps the file count
    if (o.spill_files == 0)
        while (want_files > 1 && threads * want_files * 4096 > o.memory_bytes / 8) want_files /= 2;
    const std::uint32_t files = static_cast<std::uint32_t>(std::min<std::uint64_t>(o.partitions, want_files));
    const std::uint32_t file_shift = log2_exact(files);
    const std::uint32_t per_file = o.partitions / files;
    raise_open_file_limit(files);

    const fs::path spill_dir = fs::path(o.tmp_dir.empty() ? fs::temp_directory_path() : fs::path(o.tmp_dir)) /
        ("panvar-background-" + std::to_string(::getpid()) + "-" +
         std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    fs::create_directories(spill_dir);
    struct SpillDirGuard {
        fs::path p;
        ~SpillDirGuard() { std::error_code ec; fs::remove_all(p, ec); }
    } spill_guard{spill_dir};
    JobLock spill_lock(spill_dir.string());

    struct SpillFile {
        File f;
        std::mutex m;
        std::string path;
        std::unique_ptr<DedupTable> dedup;
        std::string out;   // records the table flushed, to be written
    };
    std::vector<SpillFile> spill(files);
    // Memory in the spill phase: per-thread buffers up to an eighth of the budget, dedup tables up to three eighths.
    const std::size_t flush_at = std::clamp<std::size_t>(o.memory_bytes / 8 / (threads * files), std::size_t{4} << 10, std::size_t{4} << 20);
    const std::size_t table_budget = o.memory_bytes * 3 / 8 / files;
    for (std::uint32_t b = 0; b < files; ++b) {
        spill[b].path = (spill_dir / ("s" + std::to_string(b) + ".spill")).string();
        spill[b].f.fp = std::fopen(spill[b].path.c_str(), "wb");
        if (spill[b].f.fp == nullptr) throw std::runtime_error("background: cannot write " + spill[b].path);
        if (o.dedup && table_budget >= (std::size_t{1} << 20)) spill[b].dedup = std::make_unique<DedupTable>(table_budget);
    }
    log("spill: " + std::to_string(files) + " files in " + spill_dir.string() + " (estimated " + std::to_string(bases) +
        " bases" + (spill[0].dedup ? ", identical super-k-mers stored once" : "") + ")");

    // Anchor lookups, one per coordinate reference.
    std::vector<std::unique_ptr<anchors::AnchorLookup>> lookups;
    for (const AnchorSet& s : o.anchor_sets) lookups.push_back(std::make_unique<anchors::AnchorLookup>(*s.ref));

    // Phase 1: readers cut genomes into items (FASTA: chunks repeating the previous k-1 bases; graphs: node
    // and junction pieces with their traversal counts), workers cut items into super-k-mers appended to
    // their partition's spill file, and record anchors.
    SegmentStats st;
    st.spill_files = files;
    std::vector<Genome> genomes = o.genomes;
    std::vector<std::vector<std::string>> contig_names(genomes.size());
    std::vector<std::vector<std::uint64_t>> contig_lengths(genomes.size());
    std::vector<std::vector<std::vector<anchors::Hit>>> hits(o.anchor_sets.size());   // per set, per thread
    for (auto& h : hits) h.resize(threads + 1);   // the last for the reader (graph anchors)
    std::atomic<std::uint64_t> acgt{0}, occurrences{0};
    std::exception_ptr failure;
    std::mutex failure_m;
    ItemQueue queue(threads * 2);
    std::atomic<bool> aborted_run{false};
    auto record_failure = [&] {
        std::lock_guard<std::mutex> lk(failure_m);
        if (!failure) failure = std::current_exception();
        aborted_run = true;
        queue.abort();
    };
    {
        std::vector<std::thread> workers;
        for (std::size_t t = 0; t < threads; ++t) {
            workers.emplace_back([&, t] {
                try {
                    std::vector<std::string> bufs(files);
                    auto flush = [&](std::uint32_t b) {
                        if (bufs[b].empty()) return;
                        SpillFile& sf = spill[b];
                        std::lock_guard<std::mutex> lk(sf.m);
                        if (sf.dedup) {
                            const auto* p = reinterpret_cast<const unsigned char*>(bufs[b].data());
                            const auto* end = p + bufs[b].size();
                            while (p < end) {
                                std::uint64_t j, len, mult;
                                get_varint(p, end, j); get_varint(p, end, len); get_varint(p, end, mult);
                                sf.dedup->add(j, len, mult, p, sf.out);
                                p += (len + 3) / 4;
                            }
                            if (!sf.out.empty()) {
                                if (std::fwrite(sf.out.data(), 1, sf.out.size(), sf.f.fp) != sf.out.size())
                                    throw std::runtime_error("background: cannot write " + sf.path);
                                sf.out.clear();
                            }
                        } else if (std::fwrite(bufs[b].data(), 1, bufs[b].size(), sf.f.fp) != bufs[b].size()) {
                            throw std::runtime_error("background: cannot write " + sf.path);
                        }
                        bufs[b].clear();
                    };
                    Item item;
                    while (queue.pop(item)) {
                        std::uint64_t local = 0;
                        for_each_super_kmer(item.seq.data(), item.seq.size(), o.partitions, [&](std::size_t start, std::size_t len, std::uint32_t part) {
                            const std::uint32_t b = part & (files - 1);
                            put_varint(bufs[b], part >> file_shift);
                            put_varint(bufs[b], len);
                            put_varint(bufs[b], item.mult);
                            pack_bases(bufs[b], item.seq.data() + start, len);
                            local += len - kK + 1;
                            if (bufs[b].size() >= flush_at) flush(b);
                        });
                        occurrences += local * item.mult;
                        if (item.anchors)
                            for (std::size_t s = 0; s < lookups.size(); ++s)
                                anchors::scan_anchors(*lookups[s], item.seq.data(), item.seq.size(), item.offset, item.genome, item.contig, hits[s][t]);
                    }
                    for (std::uint32_t b = 0; b < files; ++b) flush(b);
                } catch (...) {
                    record_failure();
                }
            });
        }
        try {
            for (std::size_t gi = 0; gi < genomes.size(); ++gi) {
                Genome& g = genomes[gi];
                std::error_code ec;
                g.bytes = fs::file_size(g.path, ec);
                std::uint64_t checksum = 0x5bd1e9955bd1e995ULL, genome_bases = 0;
                if (kinds[gi] == InputKind::Gfa) {
                    graph::GraphCounts gc = graph::count_graph(g.path, static_cast<std::uint32_t>(gi), threads,
                        [&](std::string seq, std::uint64_t mult) {
                            Item it;
                            it.seq = std::move(seq);
                            it.mult = mult;
                            it.genome = static_cast<std::uint32_t>(gi);
                            it.anchors = false;
                            return queue.push(std::move(it));
                        });
                    if (aborted_run) break;
                    checksum = gc.checksum;
                    genome_bases = gc.bases;
                    contig_names[gi] = gc.contigs;
                    contig_lengths[gi] = gc.contig_lengths;
                    for (std::size_t s = 0; s < lookups.size() && !aborted_run; ++s) {
                        std::vector<anchors::Hit> gh = graph::graph_anchors(g.path, *lookups[s], static_cast<std::uint32_t>(gi), threads);
                        auto& dst = hits[s][threads];
                        if (dst.empty()) dst = std::move(gh);
                        else dst.insert(dst.end(), gh.begin(), gh.end());
                    }
                } else {
                    GzLineReader in(g.path);
                    if (!in.ok()) throw std::runtime_error("background: cannot read " + g.path);
                    std::string line, chunk;
                    std::uint32_t contig = 0;
                    std::uint64_t chunk_offset = 0, contig_pos = 0;
                    bool have = false, aborted = false;
                    chunk.reserve(chunk_bases + kK);
                    auto ship = [&](bool keep_tail) {
                        if (chunk.size() >= kK) {
                            std::string tail = keep_tail ? chunk.substr(chunk.size() - (kK - 1)) : std::string();
                            Item it;
                            it.seq = std::move(chunk);
                            it.genome = static_cast<std::uint32_t>(gi);
                            it.contig = contig;
                            it.offset = chunk_offset;
                            if (!queue.push(std::move(it))) aborted = true;
                            chunk = std::move(tail);
                            chunk_offset = contig_pos - chunk.size();
                            chunk.reserve(chunk_bases + kK);
                        } else if (!keep_tail) {
                            chunk.clear();
                            chunk_offset = contig_pos;
                        }
                    };
                    while (!aborted && in.getline(line)) {
                        checksum = hash_bytes(line.data(), line.size(), checksum);
                        if (line.empty()) continue;
                        if (line[0] == '>') {
                            if (have) { ship(false); contig_lengths[gi].push_back(contig_pos); }
                            const auto sp = line.find_first_of(" \t");
                            contig_names[gi].push_back(line.substr(1, sp == std::string::npos ? std::string::npos : sp - 1));
                            contig = static_cast<std::uint32_t>(contig_names[gi].size() - 1);
                            contig_pos = 0;
                            chunk.clear();
                            chunk_offset = 0;
                            have = true;
                            continue;
                        }
                        if (!have) throw std::runtime_error("background: " + g.path + ": sequence before the first '>' line");
                        for (std::size_t at = 0; at < line.size();) {   // long lines in pieces, so chunks stay bounded
                            const std::size_t take = std::min(line.size() - at, chunk_bases - chunk.size());
                            for (std::size_t i = at; i < at + take; ++i) genome_bases += kmer::base_code(line[i]) >= 0;
                            chunk.append(line, at, take);
                            contig_pos += take;
                            at += take;
                            if (contig_pos >= (1ULL << 32) && !lookups.empty())
                                throw std::runtime_error("background: " + g.path + ": a contig of 4 Gb or more; anchors keep 32-bit positions");
                            if (chunk.size() >= chunk_bases) ship(true);
                        }
                    }
                    if (have) { ship(false); contig_lengths[gi].push_back(contig_pos); }
                    if (aborted) break;
                }
                g.checksum = hex64(checksum);
                acgt += genome_bases;
                log("read " + g.name + ": " + std::to_string(genome_bases) + " ACGT bases" + (genome_bases == 0 ? " (no sequence)" : ""));
            }
        } catch (...) {
            record_failure();
        }
        queue.close();
        for (auto& w : workers) w.join();
    }
    std::uint64_t dedup_hits = 0;
    if (!failure) {
        try {
            for (auto& s : spill) {
                if (s.dedup) {
                    s.dedup->flush(s.out);
                    dedup_hits += s.dedup->hits();
                    if (!s.out.empty() && std::fwrite(s.out.data(), 1, s.out.size(), s.f.fp) != s.out.size())
                        throw std::runtime_error("background: cannot write " + s.path);
                    s.dedup.reset();
                    std::string().swap(s.out);
                }
            }
        } catch (...) {
            failure = std::current_exception();
        }
    }
    for (auto& s : spill)
        if (!s.f.close() && !failure) failure = std::make_exception_ptr(std::runtime_error("background: cannot write " + s.path));
    if (failure) std::rethrow_exception(failure);
    st.bases = acgt;
    st.kmers = occurrences;
    for (const auto& s : spill) { std::error_code ec; st.spill_bytes += fs::file_size(s.path, ec); }
    log("spilled " + std::to_string(st.kmers) + " k-mer occurrences in " + std::to_string(st.spill_bytes) + " bytes" +
        (dedup_hits > 0 ? " (" + std::to_string(dedup_hits) + " repeated super-k-mers folded)" : ""));

    // Phase 2: threads take spill files one at a time, count their partitions within the memory budget and
    // append the encodings to parts-0.bin.
    fs::create_directories(o.segment_dir);
    PartsWriter writer((fs::path(o.segment_dir) / "parts-0.bin").string(), o.partitions, 0, o.partitions);
    std::atomic<std::uint64_t> run_bytes{0};
    parallel_for(files, threads, [&](std::size_t b) {
        PartitionCounter counter(per_file, static_cast<std::size_t>(count_budget), spill[b].path + ".runs");
        {
            File in(spill[b].path, "rb");
            if (in.fp == nullptr) throw std::runtime_error("background: cannot read " + spill[b].path);
            std::vector<unsigned char> buf(std::size_t{4} << 20);
            std::size_t have = 0;
            while (true) {
                const std::size_t got = std::fread(buf.data() + have, 1, buf.size() - have, in.fp);
                if (got == 0 && std::ferror(in.fp)) throw std::runtime_error("background: cannot read " + spill[b].path);
                have += got;
                const unsigned char* p = buf.data();
                const unsigned char* end = buf.data() + have;
                while (p < end) {
                    const unsigned char* rec = p;
                    std::uint64_t j, len, mult;
                    if (!get_varint(p, end, j) || !get_varint(p, end, len) || !get_varint(p, end, mult) ||
                        static_cast<std::uint64_t>(end - p) < (len + 3) / 4) {
                        p = rec;
                        break;
                    }
                    if (j >= per_file || len < kK || mult == 0) throw std::runtime_error("background: corrupt spill file " + spill[b].path);
                    for_each_packed_kmer(p, len, [&](std::uint64_t code) { counter.add(static_cast<std::uint32_t>(j), code, mult); });
                    p += (len + 3) / 4;
                }
                const std::size_t consumed = static_cast<std::size_t>(p - buf.data()), left = have - consumed;
                if (got == 0) {
                    if (left != 0) throw std::runtime_error("background: truncated spill file " + spill[b].path);
                    break;
                }
                std::memmove(buf.data(), buf.data() + consumed, left);
                have = left;
                if (have == buf.size()) buf.resize(buf.size() * 2);   // one record larger than the buffer
            }
        }
        std::error_code ec;
        fs::remove(spill[b].path, ec);
        std::vector<std::uint64_t> codes, cnt;
        for (std::uint32_t j = 0; j < per_file; ++j) {
            counter.finish(j, codes, cnt);
            writer.put((j << file_shift) | static_cast<std::uint32_t>(b), codes, cnt);
        }
        run_bytes += counter.run_bytes();
    });
    st.run_bytes = run_bytes;
    st.bytes = writer.finish();
    st.distinct = writer.distinct();
    if (st.run_bytes > 0) log("counting wrote " + std::to_string(st.run_bytes) + " bytes of sorted runs (partitions above the memory share)");

    SegmentInfo info;
    info.partitions = o.partitions;
    info.parts.push_back({"parts-0.bin", 0, o.partitions});
    info.genomes = genomes;
    info.distinct = st.distinct;
    info.occurrences = st.kmers;
    for (std::size_t s = 0; s < o.anchor_sets.size(); ++s) {
        std::vector<anchors::Hit> all;
        for (auto& v : hits[s]) { all.insert(all.end(), v.begin(), v.end()); std::vector<anchors::Hit>().swap(v); }
        std::vector<std::string> names;
        for (const auto& g : genomes) names.push_back(g.name);
        anchors::write_segment_anchors((fs::path(o.segment_dir) / ("anchors-" + o.anchor_sets[s].name + ".bin")).string(),
                                       *o.anchor_sets[s].ref, names, contig_names, contig_lengths, all);
        info.anchor_sets.emplace_back(o.anchor_sets[s].name, o.anchor_sets[s].ref->fingerprint());
    }
    write_segment_info(o.segment_dir, info);
    fsync_path(o.segment_dir);
    return st;
}

// ---------------------------------------------------------------------------------------------------------
// Reading

Mapping::Mapping(const std::string& path) {
    const int fd = ::open(path.c_str(), O_RDONLY);
    if (fd < 0) throw std::runtime_error("background: cannot read " + path);
    struct stat sb {};
    if (::fstat(fd, &sb) != 0) { ::close(fd); throw std::runtime_error("background: cannot read " + path); }
    size_ = static_cast<std::size_t>(sb.st_size);
    if (size_ > 0) {
        data_ = ::mmap(nullptr, size_, PROT_READ, MAP_SHARED, fd, 0);
        if (data_ == MAP_FAILED) { data_ = nullptr; ::close(fd); throw std::runtime_error("background: cannot map " + path); }
    }
    ::close(fd);
}

Mapping::~Mapping() { if (data_ != nullptr) ::munmap(data_, size_); }

std::uint64_t parts_distinct(const std::string& path) {
    const Mapping m(path);
    const auto* w = static_cast<const std::uint64_t*>(m.data());
    if (m.size() < kHeaderWords * 8 || w[0] != kPartsMagic || w[1] != kVersion) throw std::runtime_error("background: " + path + " is incomplete");
    return w[7];
}

Segment::Segment(const std::string& segment_dir) : dir_(segment_dir), info_(read_segment_info(segment_dir)) {
    partitions_ = info_.partitions;
    if (partitions_ == 0 || partitions_ > kMaxPartitions || (partitions_ & (partitions_ - 1)) != 0)
        throw std::runtime_error("background: " + dir_ + ": bad partition count");
    where_.assign(partitions_, Where{});
    std::vector<char> seen(partitions_, 0);
    for (const PartsFile& pf : info_.parts) {
        const std::string path = (fs::path(dir_) / pf.file).string();
        maps_.push_back(std::make_unique<Mapping>(path));
        const Mapping& m = *maps_.back();
        const auto* w = static_cast<const std::uint64_t*>(m.data());
        const std::uint64_t total = m.size() / 8;
        auto bad = [&](const std::string& why) { return std::runtime_error("background: " + path + " " + why); };
        if (m.size() % 8 != 0 || total < kHeaderWords) throw bad("is not a complete parts file");
        if (w[0] != kPartsMagic) throw bad("is incomplete (interrupted write) or not a parts file of this version");
        if (w[1] != kVersion || w[2] != kK || w[3] != kM || w[4] != partitions_ || w[5] != pf.first || w[6] != pf.end || w[8] != total)
            throw bad("does not match its segment");
        if (pf.first >= pf.end || pf.end > partitions_ || total < kHeaderWords + 2ULL * (pf.end - pf.first)) throw bad("has a bad partition range");
        const std::uint64_t* table = w + kHeaderWords;
        const std::uint64_t data_start = kHeaderWords + 2ULL * (pf.end - pf.first);
        for (std::uint32_t p = pf.first; p < pf.end; ++p) {
            if (seen[p]++) throw bad("repeats partition " + std::to_string(p));
            const std::uint64_t off = table[2 * (p - pf.first)], size = table[2 * (p - pf.first) + 1];
            if (off < data_start || size < 7 || size > total || off > total - size) throw bad("has a partition outside the file");
            // the Elias-Fano header's fields must agree with each other and fit the partition's words
            const std::uint64_t* h = w + off;
            const std::uint64_t n = h[0], lb = h[1], max_high = h[2], hw = h[3], sc = h[4], lw = h[5];
            if (lb > 62 || n > size * 64 || max_high > (std::uint64_t{1} << 62) || hw != (n + max_high + 1 + 63) / 64 ||
                lw != (n * lb + 63) / 64 || sc == 0 || sc > hw * 64 / ef::kSelectSample + 1 || 6 + hw + sc + lw + 1 > size)
                throw bad("has a corrupt partition");
            const ef::View codes(h);
            const std::uint64_t ch = w[off + codes.words()];
            const std::uint64_t width = ch & 0xff, n_exc = ch >> 8;
            if ((n > 0 && (width == 0 || width > 64)) || n_exc > n || 1 + (n * width + 63) / 64 + 2 * n_exc + codes.words() != size)
                throw bad("has a corrupt partition");
            const counts::View cnt(w + off + codes.words(), codes.size());
            where_[p] = {w + off, size};
        }
    }
    for (std::uint32_t p = 0; p < partitions_; ++p)
        if (!seen[p]) throw std::runtime_error("background: " + dir_ + " lacks partition " + std::to_string(p));
}

Segment::Part Segment::part(std::uint32_t p) const {
    Part pt;
    pt.codes = ef::View(where_[p].words);
    pt.counts = counts::View(where_[p].words + pt.codes.words(), pt.codes.size());
    return pt;
}

std::uint64_t Segment::partition_size(std::uint32_t partition) const { return ef::View(where_[partition].words).size(); }

std::int64_t Segment::find(std::uint64_t canonical_kmer, std::uint32_t partition, std::uint64_t* count) const {
    const Part pt = part(partition);
    const std::int64_t i = pt.codes.find(canonical_kmer);
    if (i >= 0 && count != nullptr) *count = pt.counts.get(static_cast<std::uint64_t>(i));
    return i;
}

std::uint64_t Segment::count(std::uint64_t canonical_kmer, std::uint32_t partition) const {
    const Part pt = part(partition);
    const std::int64_t i = pt.codes.find(canonical_kmer);
    return i < 0 ? 0 : pt.counts.get(static_cast<std::uint64_t>(i));
}

void Segment::for_each(std::uint32_t partition, const std::function<void(std::uint64_t, std::uint64_t)>& f) const {
    const Part pt = part(partition);
    ef::View::Cursor cc(pt.codes);
    counts::View::Cursor nc(pt.counts);
    for (std::uint64_t code; cc.next(code);) f(code, nc.next());
}

void merge_partition(const std::vector<const Segment*>& segs, std::uint32_t p,
                     const std::function<void(std::uint64_t, std::uint64_t)>& f) {
    if (segs.size() == 1) { segs[0]->for_each(p, f); return; }
    std::vector<Segment::Part> parts;
    for (const Segment* s : segs) parts.push_back(s->part(p));
    std::vector<ef::View::Cursor> cc;
    std::vector<counts::View::Cursor> nc;
    for (const auto& pt : parts) { cc.emplace_back(pt.codes); nc.emplace_back(pt.counts); }
    using E = std::pair<std::uint64_t, std::size_t>;
    std::priority_queue<E, std::vector<E>, std::greater<E>> heap;
    for (std::size_t s = 0; s < parts.size(); ++s) { std::uint64_t c; if (cc[s].next(c)) heap.emplace(c, s); }
    bool have = false;
    std::uint64_t cur = 0, sum = 0;
    while (!heap.empty()) {
        const auto [code, s] = heap.top();
        heap.pop();
        const std::uint64_t n = nc[s].next();
        if (have && code == cur) sum += n;
        else { if (have) f(cur, sum); cur = code; sum = n; have = true; }
        std::uint64_t c;
        if (cc[s].next(c)) heap.emplace(c, s);
    }
    if (have) f(cur, sum);
}

SegmentStats merge_range(const MergeRangeOptions& o) {
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
    if (o.first >= o.end || o.end > partitions) throw std::runtime_error("background: bad merge range");
    PartsWriter writer(o.out_path, partitions, o.first, o.end);
    parallel_for(o.end - o.first, o.threads, [&](std::size_t i) {
        std::vector<std::uint64_t> codes, cnt;
        merge_partition(raw, o.first + static_cast<std::uint32_t>(i), [&](std::uint64_t c, std::uint64_t n) { codes.push_back(c); cnt.push_back(n); });
        writer.put(o.first + static_cast<std::uint32_t>(i), codes, cnt);
    });
    SegmentStats st;
    st.bytes = writer.finish();
    st.distinct = writer.distinct();
    return st;
}

}  // namespace panvar::bg
