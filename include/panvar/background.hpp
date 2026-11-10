#pragma once

// The background count store: every distinct canonical 31-mer of a set of genomes with its exact count,
// split into partitions by a hash of the k-mer's minimizer, so building, counting and merging run
// partition by partition on any number of threads. A batch of genomes is one segment; a store is one or
// more segments, and a k-mer's count is the sum over them.

#include "panvar/elias_fano.hpp"
#include "panvar/kmer.hpp"

#include <cstddef>
#include <cstdint>
#include <deque>
#include <functional>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace panvar::bg {

constexpr std::size_t kK = 31;             // k-mer length
constexpr std::size_t kM = 11;             // minimizer length
constexpr std::uint64_t kMinimizerSeed = 0x6a09e667f3bcc908ULL;   // keeps minimizer hashes apart from syncmers'

// The hash of a canonical m-mer as a minimizer.
inline std::uint64_t minimizer_hash(std::uint64_t canonical_mmer) { return kmer::mix64(canonical_mmer ^ kMinimizerSeed); }

// The partition of a canonical k-mer: its smallest canonical m-mer hash, modulo the partition count (a
// power of two). Both strands give the same partition.
inline std::uint32_t partition_of(std::uint64_t kmer_code, std::uint32_t partitions) {
    const std::uint64_t mmask = kmer::mask_of(kM);
    std::uint64_t best = ~0ULL;
    for (std::size_t j = 0; j + kM <= kK; ++j) {
        const std::uint64_t h = minimizer_hash(kmer::canonical((kmer_code >> (2 * (kK - kM - j))) & mmask, kM));
        if (h < best) best = h;
    }
    return static_cast<std::uint32_t>(best & (partitions - 1));
}

// Calls f(start, length, partition) for every super-k-mer of seq[0, n): a maximal run of consecutive
// 31-mers (no base other than ACGT inside) that share a partition, as the bases they span.
template <class F>
void for_each_super_kmer(const char* seq, std::size_t n, std::uint32_t partitions, F&& f) {
    if (n < kK) return;
    const std::uint64_t mmask = kmer::mask_of(kM);
    const std::size_t mrc = 2 * (kM - 1);
    std::uint64_t mf = 0, mr = 0;
    std::size_t filled = 0;
    std::deque<std::pair<std::uint64_t, std::size_t>> mins;   // (hash, m-mer start), increasing hashes
    std::size_t run_start = 0, run_end = 0;                     // bases of the current run
    std::uint32_t run_part = 0;
    bool in_run = false;
    auto close_run = [&] { if (in_run) f(run_start, run_end - run_start, run_part); in_run = false; };
    for (std::size_t i = 0; i < n; ++i) {
        const int b = kmer::base_code(seq[i]);
        if (b < 0) { close_run(); mf = mr = 0; filled = 0; mins.clear(); continue; }
        mf = ((mf << 2) | static_cast<std::uint64_t>(b)) & mmask;
        mr = (mr >> 2) | (static_cast<std::uint64_t>(3 - b) << mrc);
        ++filled;
        if (filled >= kM) {
            const std::uint64_t h = minimizer_hash(mf <= mr ? mf : mr);
            while (!mins.empty() && mins.back().first > h) mins.pop_back();
            mins.emplace_back(h, i + 1 - kM);
        }
        if (filled >= kK) {
            const std::size_t kstart = i + 1 - kK;
            while (mins.front().second < kstart) mins.pop_front();
            const auto part = static_cast<std::uint32_t>(mins.front().first & (partitions - 1));
            if (in_run && part == run_part) {
                run_end = i + 1;
            } else {
                close_run();
                in_run = true; run_part = part; run_start = kstart; run_end = i + 1;
            }
        }
    }
    close_run();
}

struct Genome {
    std::string name;
    std::string path;   // FASTA, plain or gzip/bgzip
};

struct BuildOptions {
    std::vector<Genome> genomes;
    std::string segment_dir;          // written: <segment_dir>/parts.bin
    std::string tmp_dir;              // spill files go in a fresh folder under it
    std::uint32_t partitions = 4096;  // a power of two, at most 65536
    std::size_t threads = 2;
    std::size_t memory_bytes = std::size_t{2} << 30;   // for the counting phase, across threads
    std::uint32_t spill_files = 0;    // 0: as many as the memory needs; else this many (a power of two)
    std::function<void(const std::string&)> log;       // progress lines
};

struct SegmentStats {
    std::uint64_t bases = 0;          // ACGT bases read
    std::uint64_t kmers = 0;          // 31-mer occurrences
    std::uint64_t distinct = 0;       // distinct canonical 31-mers
    std::uint64_t bytes = 0;          // parts.bin size
    std::uint64_t spill_bytes = 0;    // temporary files written
    std::uint32_t spill_files = 0;
};

SegmentStats build_segment(const BuildOptions& options);

struct MergeOptions {
    std::vector<std::string> segment_dirs;   // all with the same partition count
    std::string out_segment_dir;
    std::size_t threads = 2;
};

// Merges segments into one, partition by partition on threads: each k-mer once, its counts summed.
SegmentStats merge_segments(const MergeOptions& options);

// A segment, read in place.
class Segment {
public:
    explicit Segment(const std::string& segment_dir);
    ~Segment();
    Segment(const Segment&) = delete;
    Segment& operator=(const Segment&) = delete;

    std::uint32_t partitions() const { return partitions_; }
    std::uint64_t distinct() const { return distinct_; }
    // The count of a canonical 31-mer (0 when absent).
    std::uint64_t count(std::uint64_t canonical_kmer) const;
    // Every (k-mer, count) of one partition, in code order.
    void for_each(std::uint32_t partition, const std::function<void(std::uint64_t, std::uint64_t)>& f) const;
    // Distinct k-mers in a partition; a k-mer's index within its partition (-1 when absent) and its count.
    std::uint64_t partition_size(std::uint32_t partition) const;
    std::int64_t find(std::uint64_t canonical_kmer, std::uint32_t partition, std::uint64_t* count) const;

private:
    struct Part {
        ef::View codes;
        const std::uint8_t* counts = nullptr;
        const std::uint64_t* overflow = nullptr;   // (index << 32 | count), sorted by index
        std::uint64_t overflow_n = 0;
    };
    Part part(std::uint32_t p) const;
    std::uint64_t count_at(const Part& pt, std::uint64_t index) const;

    void* map_ = nullptr;
    std::size_t map_len_ = 0;
    const std::uint64_t* words_ = nullptr;
    const std::uint64_t* offsets_ = nullptr;   // per partition: word offset of its blob
    std::uint32_t partitions_ = 0;
    std::uint64_t distinct_ = 0;
};

// A store: every segment of it, read in place; a k-mer's count is the sum over the segments.
class Store {
public:
    explicit Store(const std::string& store_dir);
    std::uint32_t partitions() const { return partitions_; }
    std::size_t segments() const { return segs_.size(); }
    const std::vector<std::string>& segment_dirs() const { return dirs_; }
    std::uint64_t count(std::uint64_t canonical_kmer) const;
    // Every (k-mer, summed count) of one partition, in code order.
    void for_each(std::uint32_t partition, const std::function<void(std::uint64_t, std::uint64_t)>& f) const;

private:
    std::vector<std::string> dirs_;
    std::vector<std::unique_ptr<Segment>> segs_;
    std::uint32_t partitions_ = 0;
};

// The segment folders of a store (seg-<n>), in order.
std::vector<std::string> segment_dirs_of(const std::string& store_dir);

struct CompactOptions {
    std::string store_dir;            // a store of one segment (merge first)
    std::string out_path;             // the SSHash dictionary, <store>.sshash by convention
    std::string tmp_dir;              // the string set and SSHash's temporary files
    std::size_t threads = 2;
    std::size_t minimizer = 0;        // SSHash's minimizer length; 0 picks one from the k-mer count
    std::size_t ram_gb = 4;           // SSHash's construction memory limit
    std::function<void(const std::string&)> log;
};

struct CompactStats {
    std::uint64_t kmers = 0;          // distinct k-mers written
    std::uint64_t strings = 0;        // simplitigs
    std::uint64_t string_bases = 0;   // their total length
    std::uint64_t bytes = 0;          // the dictionary's size
    std::size_t minimizer = 0;
};

// Rewrites a merged store as a string set (each k-mer once, with its count) and builds an SSHash weighted
// dictionary from it.
CompactStats compact_store(const CompactOptions& options);

// Counts of canonical 31-mers from either form of a background: a store folder or a compacted dictionary.
class Background {
public:
    virtual ~Background() = default;
    virtual std::uint64_t count(std::uint64_t canonical_kmer) const = 0;
    virtual std::uint64_t distinct() const = 0;
    virtual std::string kind() const = 0;
    // A store folder, or a .sshash file written by compact_store.
    static std::unique_ptr<Background> open(const std::string& path);
};

} // namespace panvar::bg
