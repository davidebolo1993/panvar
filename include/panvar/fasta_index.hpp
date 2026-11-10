#pragma once

// Reading intervals of a FASTA: through its .fai (and .gzi for a bgzip file) when it has them, else by
// one streaming pass over the file. Coordinates are 0-based, half-open.

#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace panvar::fasta {

struct Interval {
    std::string contig;
    std::uint64_t start = 0, end = 0;
};

class IndexedFasta {
public:
    // Opens `path` with its .fai; a gzip file must be bgzip with a .gzi. Throws when it cannot.
    explicit IndexedFasta(const std::string& path);
    ~IndexedFasta();
    IndexedFasta(const IndexedFasta&) = delete;
    IndexedFasta& operator=(const IndexedFasta&) = delete;
    static bool indexed(const std::string& path);
    bool has(const std::string& contig) const { return fai_.count(contig) != 0; }
    std::uint64_t length(const std::string& contig) const;
    // The bases of [start, end) of a contig, clamped to its length, newlines removed.
    std::string fetch(const std::string& contig, std::uint64_t start, std::uint64_t end);

private:
    struct Entry { std::uint64_t length, offset, line_bases, line_width; };
    void read_bytes(std::uint64_t offset, std::uint64_t n, std::string& out);
    std::string path_;
    std::map<std::string, Entry> fai_;
    bool bgzf_ = false;
    std::vector<std::pair<std::uint64_t, std::uint64_t>> gzi_;   // (compressed, uncompressed) block starts
    std::FILE* f_ = nullptr;
    std::uint64_t cached_coff_ = ~0ULL, cached_uoff_ = 0, cached_next_ = 0;   // the last block read, and the next one's offset
    std::string cached_;
};

// The sequences of many intervals of one FASTA, in order: indexed when possible, else one streaming pass.
std::vector<std::string> fetch_intervals(const std::string& path, const std::vector<Interval>& intervals);

}  // namespace panvar::fasta
