#pragma once

// Counting a pangenome graph's haplotypes without spelling them: every k-mer inside a node counts once per
// traversal of the node, and k-mers that cross node ends are read from the distinct junction strings
// (k-1 bases on each side) with how many haplotypes walk them. Haplotypes are the graph's P lines and W
// lines; anchors are placed per haplotype from per-node anchor lists.

#include "panvar/anchors.hpp"

#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace panvar::bg::graph {

struct GraphCounts {
    std::uint64_t checksum = 0;          // of the graph's lines
    std::uint64_t bases = 0;             // ACGT bases of all haplotypes, as if spelled
    std::uint64_t haplotypes = 0;        // walks
    std::vector<std::string> contigs;    // haplotype contigs: sample#hap#seq (PanSN) or path names
    std::vector<std::uint64_t> contig_lengths;   // the walk's end coordinate (its start offset plus its length)
};

// Streams the graph's k-mer material to `push(sequence, multiplicity)` (false stops); every haplotype's
// k-mers are produced exactly as often as in its spelled sequence. A node longer than `chunk_bases` goes in
// pieces overlapping by k - 1 bases; junction strings are flushed past `junction_budget` bytes. The node
// sequences themselves are held in memory (about a byte a base).
GraphCounts count_graph(const std::string& gfa, std::uint32_t genome, std::size_t threads,
                        const std::function<bool(std::string, std::uint64_t)>& push,
                        std::size_t chunk_bases = std::size_t{1} << 24, std::size_t junction_budget = std::size_t{1} << 30);

// The checksum counting records for a graph genome (GraphCounts::checksum).
std::uint64_t sequence_checksum(const std::string& gfa);

// Anchors placed on every haplotype contig (indexes as in GraphCounts::contigs).
std::vector<anchors::Hit> graph_anchors(const std::string& gfa, const anchors::AnchorLookup& lookup, std::uint32_t genome,
                                        std::size_t threads);

// The haplotype contig names and lengths, in the order counting and anchoring index them.
std::vector<std::string> graph_contigs(const std::string& gfa, std::vector<std::uint64_t>* lengths = nullptr);

// A graph's haplotypes, held to read stretches of them by contig coordinates.
class Haplotypes {
public:
    explicit Haplotypes(const std::string& gfa);
    ~Haplotypes();
    // The bases of [start, end) of a haplotype contig (W-line coordinates include the walk's start offset),
    // one piece per walk the interval crosses: counting never joined walks, so no k-mer spans two pieces.
    struct Piece {
        std::uint64_t start = 0;
        std::string seq;
    };
    std::vector<Piece> fetch(const std::string& contig, std::uint64_t start, std::uint64_t end) const;
private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

// The haplotypes spelled as FASTA (for checks).
void spell_graph(const std::string& gfa, const std::string& fasta_out);

}  // namespace panvar::bg::graph
