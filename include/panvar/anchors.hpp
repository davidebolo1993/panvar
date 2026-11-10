#pragma once

// Anchor maps: where reference-unique k-mers fall in each background genome, so that a locus given in
// reference coordinates can be found in any genome without aligning it. An anchor is a closed syncmer of
// the coordinate reference (k = 31, s = 11), kept by a hash at about one per kilobase, and only when it
// occurs exactly once in the reference. Each segment keeps, per coordinate reference, one anchor-major file:
// for every anchor, its places in every genome of the segment, read in place.

#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace panvar::bg { class Mapping; }

namespace panvar::anchors {

constexpr std::size_t kK = 31;
constexpr std::size_t kS = 11;
constexpr std::uint32_t kSampleOneIn = 100;   // of closed syncmers (about one per 10.5 bases): one per ~1 kb
constexpr std::size_t kMaxOccurrences = 4;    // places kept per anchor and genome; more marks it repeated there

// Whether a canonical k-mer is sampled as an anchor candidate (before the uniqueness test).
bool sampled(std::uint64_t canonical_kmer);

struct RefPlace {
    std::uint32_t contig = 0;
    std::uint32_t pos = 0;        // 0-based start of the k-mer
    bool forward = true;          // the reference reads as the anchor's canonical code
};

// The reference's anchors: canonical code and single place, ordered by contig then position.
struct ReferenceAnchors {
    std::string genome;                       // the coordinate reference's genome name
    std::vector<std::string> contigs;
    std::vector<std::uint64_t> contig_lengths;
    std::vector<std::uint64_t> code;
    std::vector<RefPlace> at;
    std::uint64_t fingerprint() const;        // of the codes, places and parameters
};

ReferenceAnchors build_reference(const std::string& genome, const std::string& fasta, std::size_t threads,
                                 const std::function<void(const std::string&)>& log = {});
void save(const ReferenceAnchors& ref, const std::string& path);
ReferenceAnchors load_reference(const std::string& path);

// One anchor found in a genome while it is read.
struct Hit {
    std::uint32_t anchor = 0;
    std::uint32_t genome = 0;     // index among the segment's genomes
    std::uint32_t contig = 0;     // index among that genome's contigs
    std::uint32_t pos = 0;
    bool forward = true;          // the genome reads as the anchor's canonical code
};

// Code -> anchor id, for scanning genomes.
class AnchorLookup {
public:
    explicit AnchorLookup(const ReferenceAnchors& ref);
    // -1 when the code is not an anchor.
    std::int64_t find(std::uint64_t canonical_kmer) const;
private:
    std::vector<std::uint64_t> slots_;   // open addressing: code << 1 | 1, 0 empty
    std::vector<std::uint32_t> ids_;
    std::uint64_t mask_ = 0;
};

// Calls f(hit) for every anchor in one stretch of a genome: seq[0, n) starting at `offset` of `contig`.
void scan_anchors(const AnchorLookup& lookup, const char* seq, std::size_t n, std::uint64_t offset,
                  std::uint32_t genome, std::uint32_t contig, std::vector<Hit>& out);

// Every anchor in a FASTA genome (contig names filled in file order).
std::vector<Hit> map_fasta(const AnchorLookup& lookup, const std::string& fasta, std::uint32_t genome, std::size_t threads,
                           std::vector<std::string>& contigs, std::vector<std::uint64_t>& lengths);

// Writes a segment's anchor file for one reference: per anchor, its places sorted by genome, contig and
// position; an anchor seen more than kMaxOccurrences times in a genome keeps one entry marked repeated.
void write_segment_anchors(const std::string& path, const ReferenceAnchors& ref,
                           const std::vector<std::string>& genome_names,
                           const std::vector<std::vector<std::string>>& contig_names,
                           const std::vector<std::vector<std::uint64_t>>& contig_lengths, std::vector<Hit>& hits);

// One place of an anchor in a genome of the segment.
struct Place {
    std::uint32_t genome = 0, contig = 0, pos = 0;
    bool forward = true;
    bool repeated = false;        // seen more than kMaxOccurrences times in this genome (pos not meaningful)
};

// A segment's anchor file, read in place.
class SegmentAnchors {
public:
    explicit SegmentAnchors(const std::string& path);
    ~SegmentAnchors();
    std::uint64_t fingerprint() const { return fingerprint_; }
    std::size_t anchors() const { return anchors_; }
    const std::vector<std::string>& genomes() const { return genomes_; }
    const std::vector<std::string>& contigs(std::size_t genome) const { return contigs_[genome]; }
    std::uint64_t contig_length(std::size_t genome, std::size_t contig) const { return lengths_[genome][contig]; }
    // Every place of anchor a, by genome.
    void places(std::size_t a, std::vector<Place>& out) const;
    std::uint64_t bytes_touched(std::size_t a_begin, std::size_t a_end) const;
private:
    std::unique_ptr<bg::Mapping> map_;
    std::uint64_t fingerprint_ = 0;
    std::size_t anchors_ = 0;
    std::vector<std::string> genomes_;
    std::vector<std::vector<std::string>> contigs_;
    std::vector<std::vector<std::uint64_t>> lengths_;
    const std::uint64_t* first_ = nullptr;    // anchors + 1 offsets into places
    const unsigned char* places_ = nullptr;   // 12 bytes each
};

// Joins segments' anchor files (genomes in segment order) into one, for a merged segment.
void merge_segment_anchors(const std::vector<std::string>& inputs, const ReferenceAnchors& ref, const std::string& out);

// A copy of a reference interval found in a genome.
struct Span {
    std::string genome;
    std::string contig;
    std::uint64_t start = 0, end = 0;   // 0-based, half-open, on the genome's contig
    bool forward = true;                 // the copy reads in the reference's orientation
    std::size_t anchors_inside = 0;      // chained anchors from inside the interval
    std::size_t anchors_flank = 0;       // chained anchors from the flanks
    std::size_t chains = 1;              // chains joined into this span because they overlapped
};

struct LocateOptions {
    std::size_t flank = 10;              // anchors taken on each side of the interval
    std::size_t min_inside = 3;          // inside anchors a chain needs (fewer when the interval has fewer)
    std::uint64_t slack = 200000;        // gap tolerance beyond 1.5 times the reference gap
    double gap_cost_per_kb = 0.05;       // chain score: 1 per anchor, minus this per kb of gap difference
};

// The copies of reference interval [start, end) of `contig` in every genome of the given anchor files:
// anchors inside the interval plus `flank` on each side are chained per genome (same contig and strand,
// reference order, gaps within bounds, best score first); each chain holding enough inside anchors is a
// copy, extended to the nearest chained flank anchor when that one sits where the reference puts it, and
// copies that overlap on a contig are joined.
std::vector<Span> locate(const ReferenceAnchors& ref, const std::vector<const SegmentAnchors*>& segments,
                         const std::string& contig, std::uint64_t start, std::uint64_t end,
                         const LocateOptions& options = {}, std::uint64_t* bytes_read = nullptr);

}  // namespace panvar::anchors
