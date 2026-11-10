#pragma once

// Anchor maps: where reference-unique k-mers fall in each background genome, so that a locus given in
// reference coordinates can be found in any genome without aligning it. An anchor is a closed syncmer of
// the coordinate reference (k = 31, s = 11), kept by a hash at about one per kilobase, and only when it
// occurs exactly once in the reference.

#include <cstddef>
#include <cstdint>
#include <functional>
#include <string>
#include <vector>

namespace panvar::anchors {

constexpr std::size_t kK = 31;
constexpr std::size_t kS = 11;
constexpr std::uint32_t kSampleOneIn = 100;   // of closed syncmers (about one per 10.5 bases): one per ~1 kb
constexpr std::size_t kMaxOccurrences = 4;    // places recorded per anchor in a genome

// A genome's contigs and, per anchor, where it occurs.
struct Placement {
    std::uint32_t contig = 0;
    std::uint32_t pos = 0;        // 0-based start of the k-mer
    bool forward = true;          // the k-mer reads as the anchor's canonical code on this strand
};

// The reference's anchors: canonical code, contig, position, strand; ordered by contig then position.
struct ReferenceAnchors {
    std::vector<std::string> contigs;
    std::vector<std::uint64_t> code;
    std::vector<Placement> at;    // its single place in the reference
};

// One genome's map: for each anchor (by index into the reference set), up to kMaxOccurrences places.
struct GenomeMap {
    std::vector<std::string> contigs;
    std::vector<std::uint32_t> first;       // per anchor, index into places; size anchors + 1
    std::vector<Placement> places;
    std::vector<std::uint32_t> occurrences; // per anchor, how many times it occurs (may exceed what is kept)
};

// Whether a canonical k-mer is sampled as an anchor candidate (before the uniqueness test).
bool sampled(std::uint64_t canonical_kmer);

ReferenceAnchors build_reference(const std::string& fasta, std::size_t threads,
                                 const std::function<void(const std::string&)>& log = {});
GenomeMap map_genome(const ReferenceAnchors& ref, const std::string& fasta, std::size_t threads);

void save(const ReferenceAnchors& ref, const std::string& path);
void save(const GenomeMap& map, const std::string& path);
ReferenceAnchors load_reference(const std::string& path);
GenomeMap load_map(const std::string& path);

// A copy of a reference interval found in a genome: its span there and the anchors that place it.
struct Span {
    std::string contig;
    std::uint64_t start = 0, end = 0;   // 0-based, half-open, on the genome's contig
    bool forward = true;                 // the copy reads in the reference's orientation
    std::size_t anchors_inside = 0;      // chained anchors from inside the interval
    std::size_t anchors_flank = 0;       // chained anchors from the flanks
};

// The copies of reference interval [start, end) of `contig` in a genome: anchors inside the interval
// plus `flank` on each side are chained (same contig, consistent order and strand, gaps within `slack`
// of the reference's gaps plus 50%); each chain holding at least `min_inside` interior anchors is a copy.
std::vector<Span> locate(const ReferenceAnchors& ref, const GenomeMap& map, const std::string& contig,
                         std::uint64_t start, std::uint64_t end, std::size_t flank = 10,
                         std::size_t min_inside = 3, std::uint64_t slack = 200000);

} // namespace panvar::anchors
