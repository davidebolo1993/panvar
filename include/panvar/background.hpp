#pragma once

// The background count store: every distinct canonical 31-mer of a set of genomes with its exact count,
// split into partitions by a hash of the k-mer's minimizer, so building, counting and merging run
// partition by partition on any number of threads or jobs. A batch of genomes is one segment; a store is
// the segments its manifest lists, and a k-mer's count is the sum over them.

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

namespace panvar::anchors { struct ReferenceAnchors; }

namespace panvar::bg {

constexpr std::size_t kK = 31;             // k-mer length
constexpr std::size_t kM = 11;             // minimizer length
constexpr std::uint64_t kMinimizerSeed = 0x6a09e667f3bcc908ULL;   // keeps minimizer hashes apart from syncmers'
constexpr std::uint32_t kMaxPartitions = 65536;

// The hash of a canonical m-mer as a minimizer.
inline std::uint64_t minimizer_hash(std::uint64_t canonical_mmer) { return kmer::mix64(canonical_mmer ^ kMinimizerSeed); }

// The partition of a k-mer (either strand): its smallest canonical m-mer hash, modulo the partition count
// (a power of two). The m-mers of the reverse complement are read off it once, not recomputed per window.
inline std::uint32_t partition_of(std::uint64_t kmer_code, std::uint32_t partitions) {
    const std::uint64_t mmask = kmer::mask_of(kM);
    const std::uint64_t rc = kmer::reverse_complement(kmer_code, kK);
    std::uint64_t best = ~0ULL;
    for (std::size_t j = 0; j + kM <= kK; ++j) {
        const std::uint64_t f = (kmer_code >> (2 * (kK - kM - j))) & mmask;
        const std::uint64_t r = (rc >> (2 * j)) & mmask;   // the reverse complement of f
        const std::uint64_t h = minimizer_hash(f <= r ? f : r);
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

// Counts of one partition, bit-packed at the width w minimizing n*w + 128*(counts that do not fit); a
// count equal to 2^w - 1 then means "see the exceptions" (index and count, sorted by index).
namespace counts {
std::vector<std::uint64_t> encode(const std::vector<std::uint64_t>& counts);
class View {
public:
    View() = default;
    View(const std::uint64_t* w, std::uint64_t n);
    std::uint64_t get(std::uint64_t i) const;
    std::uint64_t words() const { return words_; }
    std::uint32_t width() const { return width_; }
    // Sequential reading: the counts of indexes 0, 1, 2... in turn, exceptions found without search.
    class Cursor {
    public:
        explicit Cursor(const View& v) : v_(&v) {}
        std::uint64_t next();
    private:
        const View* v_;
        std::uint64_t i_ = 0, e_ = 0;
    };
private:
    friend class Cursor;
    std::uint64_t raw(std::uint64_t i) const;
    const std::uint64_t* packed_ = nullptr;
    const std::uint64_t* exc_index_ = nullptr;
    const std::uint64_t* exc_count_ = nullptr;
    std::uint64_t n_exc_ = 0, words_ = 0, escape_ = 0;
    std::uint32_t width_ = 0;
};
}  // namespace counts

struct Genome {
    std::string name;
    std::string path;            // FASTA (plain or gzip/bgzip) or GFA (P or W lines)
    std::uint64_t bytes = 0;     // the file's size when counted
    std::string checksum;        // of its sequence lines (and names), hex; filled by counting
};

// One parts file of a segment: partitions [first, end).
struct PartsFile {
    std::string file;
    std::uint32_t first = 0, end = 0;
};

// A segment's own record, segment.tsv in its folder.
struct SegmentInfo {
    std::uint32_t partitions = 0;
    std::vector<PartsFile> parts;
    std::vector<Genome> genomes;
    std::uint64_t distinct = 0, occurrences = 0;
    std::vector<std::pair<std::string, std::uint64_t>> anchor_sets;   // coordinate reference, its fingerprint
};
SegmentInfo read_segment_info(const std::string& segment_dir);
void write_segment_info(const std::string& segment_dir, const SegmentInfo& info);

// A coordinate reference's anchor set, as the counting pass maps it in every genome.
struct AnchorSet {
    std::string name;
    const anchors::ReferenceAnchors* ref = nullptr;
};

struct BuildOptions {
    std::vector<Genome> genomes;
    std::string segment_dir;          // written: parts-0.bin, segment.tsv, anchors-<ref>.bin
    std::string tmp_dir;              // spill files go in a fresh folder under it
    std::uint32_t partitions = 4096;  // a power of two, at most 65536
    std::size_t threads = 2;
    std::size_t memory_bytes = std::size_t{2} << 30;   // a bound for the whole run
    std::uint32_t spill_files = 0;    // 0: as many as the memory needs; else this many (a power of two)
    bool dedup = true;                // store identical super-k-mers once, with their multiplicity
    std::vector<AnchorSet> anchor_sets;   // mapped in every genome during the same pass
    std::function<void(const std::string&)> log;       // progress lines
    bool progress = false;            // a progress bar while partitions are counted (on a terminal)
};

struct SegmentStats {
    std::uint64_t bases = 0;          // ACGT bases read
    std::uint64_t kmers = 0;          // 31-mer occurrences
    std::uint64_t distinct = 0;       // distinct canonical 31-mers
    std::uint64_t bytes = 0;          // parts files' size
    std::uint64_t spill_bytes = 0;    // temporary files written
    std::uint64_t run_bytes = 0;      // sorted runs written when a spill file did not fit in memory
    std::uint32_t spill_files = 0;
};

// Counts the genomes into a new segment folder (which must not exist); fills the genomes' size and
// checksum in segment.tsv.
SegmentStats build_segment(const BuildOptions& options);

// Merges partitions [first, end) of several segments into one parts file (written to a temporary name,
// then renamed), streaming their sorted partitions: each k-mer once, its counts summed.
struct MergeRangeOptions {
    std::vector<std::string> segment_dirs;
    std::string out_path;
    std::uint32_t first = 0, end = 0;
    std::size_t threads = 2;
    bool progress = false;            // a progress bar over the partitions (on a terminal)
};
SegmentStats merge_range(const MergeRangeOptions& options);

// Distinct k-mers recorded in a parts file's header (the file must be complete).
std::uint64_t parts_distinct(const std::string& path);

// A memory-mapped read-only file.
class Mapping {
public:
    explicit Mapping(const std::string& path);
    ~Mapping();
    Mapping(const Mapping&) = delete;
    Mapping& operator=(const Mapping&) = delete;
    const void* data() const { return data_; }
    std::size_t size() const { return size_; }
private:
    void* data_ = nullptr;
    std::size_t size_ = 0;
};

// A segment, read in place.
class Segment {
public:
    explicit Segment(const std::string& segment_dir);
    std::uint32_t partitions() const { return partitions_; }
    std::uint64_t distinct() const { return info_.distinct; }
    const SegmentInfo& info() const { return info_; }
    const std::string& dir() const { return dir_; }
    // The count of a canonical 31-mer (0 when absent); `partition` as partition_of gives it.
    std::uint64_t count(std::uint64_t canonical_kmer) const { return count(canonical_kmer, partition_of(canonical_kmer, partitions_)); }
    std::uint64_t count(std::uint64_t canonical_kmer, std::uint32_t partition) const;
    // Every (k-mer, count) of one partition, in code order.
    void for_each(std::uint32_t partition, const std::function<void(std::uint64_t, std::uint64_t)>& f) const;
    std::uint64_t partition_size(std::uint32_t partition) const;
    // A k-mer's index within its partition (-1 when absent) and its count.
    std::int64_t find(std::uint64_t canonical_kmer, std::uint32_t partition, std::uint64_t* count) const;

    struct Part {
        ef::View codes;
        counts::View counts;
    };
    Part part(std::uint32_t p) const;

private:
    std::string dir_;
    SegmentInfo info_;
    std::vector<std::unique_ptr<Mapping>> maps_;
    struct Where { const std::uint64_t* words = nullptr; std::uint64_t size = 0; };
    std::vector<Where> where_;    // per partition
    std::uint32_t partitions_ = 0;
};

// The store's manifest: the only list of its segments.
struct ManifestSegment {
    std::string dir;              // seg-<n>
    std::uint64_t distinct = 0, occurrences = 0;
    std::vector<Genome> genomes;
};
struct Manifest {
    std::string store_id;         // random, fixed at creation
    std::uint64_t generation = 0; // +1 on every change
    std::uint32_t partitions = 0;
    std::uint64_t next_segment = 0;
    std::vector<std::pair<std::string, std::uint64_t>> references;   // coordinate references, fingerprints
    std::vector<ManifestSegment> segments;
};
bool store_exists(const std::string& store_dir);
Manifest read_manifest(const std::string& store_dir);
// Replaces manifest.json atomically (temporary file, fsync, rename, fsync of the folder).
void write_manifest(const std::string& store_dir, const Manifest& m);
// Writes a file and syncs it, checking every step (close included); removes it on failure.
void write_file_durably(const std::string& path, const std::string& content);
// What the store counts: a hash of its genomes' names and sequence checksums (not of how they are split
// into segments, nor of its anchors), so a merge or a new coordinate reference keeps it.
std::uint64_t content_id(const Manifest& m);

// An exclusive lock on a store's .lock file, held until destroyed; waits for other holders.
class StoreLock {
public:
    explicit StoreLock(const std::string& store_dir);
    ~StoreLock();
    StoreLock(const StoreLock&) = delete;
    StoreLock& operator=(const StoreLock&) = delete;
private:
    int fd_ = -1;
};

// A running job's claim on a folder (its job.lock, locked while the process lives), so --clean can tell
// a live job's files from a dead one's.
class JobLock {
public:
    explicit JobLock(const std::string& dir);
    ~JobLock();
    JobLock(const JobLock&) = delete;
    JobLock& operator=(const JobLock&) = delete;
private:
    int fd_ = -1;
};
bool job_alive(const std::string& dir);

// A store: the segments its manifest lists, read in place; a k-mer's count is the sum over them.
class Store {
public:
    explicit Store(const std::string& store_dir);
    std::uint32_t partitions() const { return partitions_; }
    std::size_t segments() const { return segs_.size(); }
    const Manifest& manifest() const { return manifest_; }
    const Segment& segment(std::size_t i) const { return *segs_[i]; }
    std::uint64_t count(std::uint64_t canonical_kmer) const;
    // Counts of many k-mers at once, grouped by partition so each partition of each segment is visited once.
    std::vector<std::uint64_t> count_batch(const std::vector<std::uint64_t>& canonical_kmers) const;
    // Every (k-mer, summed count) of one partition, in code order (a streaming k-way merge).
    void for_each(std::uint32_t partition, const std::function<void(std::uint64_t, std::uint64_t)>& f) const;

private:
    Manifest manifest_;
    std::vector<std::unique_ptr<Segment>> segs_;
    std::uint32_t partitions_ = 0;
};

// The sorted (k-mer, summed count) stream of one partition of several segments.
void merge_partition(const std::vector<const Segment*>& segs, std::uint32_t p,
                     const std::function<void(std::uint64_t, std::uint64_t)>& f);

struct CompactOptions {
    std::string store_dir;            // a store of one segment (merge first)
    std::string out_path;             // the compacted dictionary, <store>.sshash by convention
    std::string tmp_dir;              // the string set and SSHash's temporary files
    std::size_t threads = 2;
    std::size_t minimizer = 0;        // SSHash's minimizer length; 0 picks one from the k-mer count
    std::size_t ram_gb = 4;           // SSHash's construction memory limit
    bool local_strings = false;       // build strings within partitions and glue them: never reads across
                                      // partitions, but holds every path in memory; slower
    bool quiet = false;               // SSHash's own messages go to stderr unless quiet
    std::function<void(const std::string&)> log;
};

struct CompactStats {
    std::uint64_t kmers = 0;          // distinct k-mers written
    std::uint64_t strings = 0;        // strings in the set
    std::uint64_t string_bases = 0;   // their total length
    std::uint64_t bytes = 0;          // the dictionary's size
    std::size_t minimizer = 0;
    double strings_seconds = 0, sshash_seconds = 0;
};

// Rewrites a merged store as a string set (each k-mer once, with its count) and builds an SSHash weighted
// dictionary from it, written with a header that ties it to the store's id and content.
CompactStats compact_store(const CompactOptions& options);

// Counts of canonical 31-mers from either form of a background: a store folder or a compacted dictionary.
class Background {
public:
    virtual ~Background() = default;
    virtual std::uint64_t count(std::uint64_t canonical_kmer) const = 0;
    virtual std::vector<std::uint64_t> count_batch(const std::vector<std::uint64_t>& canonical_kmers) const;
    // Exact distinct k-mers (a store sums its segments' partitions by streaming when it has several).
    virtual std::uint64_t distinct() const = 0;
    virtual std::string kind() const = 0;
    virtual std::string store_id() const = 0;
    virtual std::uint64_t content() const = 0;   // the store's content_id when it was read or compacted
    // A store folder, or a compacted dictionary written by compact_store. A dictionary next to its store
    // (<store>.sshash beside <store>) must match the store's id and content.
    static std::unique_ptr<Background> open(const std::string& path);
};

// The process's peak memory so far, for logs: allocated memory (footprint) where the system reports it,
// and resident memory (which also counts pages of mapped files).
std::string memory_report();

// Default worker threads: Slurm's allocation, else the CPU affinity mask, else the hardware count.
std::size_t default_threads();

// In-place sort of 64-bit keys (MSD radix, std::sort for small buckets).
void radix_sort(std::uint64_t* keys, std::size_t n);

} // namespace panvar::bg
