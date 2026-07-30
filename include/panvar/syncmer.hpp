#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace panvar {

// One canonical k-mer occurrence: the 2-bit code (min of forward and reverse-complement) and the
// 0-based start offset in the sequence it came from.
struct KmerOccurrence {
    std::uint64_t code = 0;
    std::size_t start = 0;
};

// A/C/G/T (either case) -> 0..3; anything else -> -1, which resets the rolling window.
int base_bits(char c);

// s for the closed-syncmer test when the caller does not pin one: max(1, min(11, (k+2)/3)).
std::size_t default_syncmer_s(std::size_t k);

// Every canonical k-mer of `sequence`, in order. k must be in [1,31] (2-bit packed into one
// uint64_t); anything else yields an empty result. Non-ACGT bases break the window.
std::vector<KmerOccurrence> collect_canonical_kmer_occurrences(
    const std::string& sequence,
    std::size_t k);

// Closed syncmer: keep a k-mer iff its minimal s-mer sits at the first or last position.
// Edgar 2021, https://doi.org/10.7717/peerj.10805. Operates on the already-canonical code.
bool is_closed_syncmer(std::uint64_t code, std::size_t k, std::size_t s);

// The closed syncmers of `sequence`, in order — the sparse anchor set a sequence is projected onto.
std::vector<KmerOccurrence> collect_syncmers(
    const std::string& sequence,
    std::size_t k,
    std::size_t s);

// The ACGT string a 2-bit code spells.
std::string decode_kmer(std::uint64_t code, std::size_t k);

} // namespace panvar
