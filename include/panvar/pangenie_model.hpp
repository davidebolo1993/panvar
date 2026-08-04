#pragma once

// A faithful reimplementation of PanGenie's genotyping model, run on panvar's own panel and read
// counts so the two can be compared on identical inputs.
//
// Why reimplement rather than run the tool: PanGenie does not build here (jellyfish 2.3.1 uses
// std::get_temporary_buffer, removed from the standard library, while PanGenie requires C++20 -- no
// setting compiles both), and jellyfish is only its k-mer counter, which we do not need because we
// already have counts. Everything below is transcribed from the source at github.com/eblerjana/pangenie
// with the file and reasoning noted at each step, so a reader can check it against the original.
//
// This is NOT an approximation of the marker rule inside our model, which is what `--marker-rule
// pangenie` does. It is their emission, their transition and their forward-backward.

#include "panvar/genotype.hpp"
#include "panvar/genotype_blocks.hpp"
#include "panvar/genotype_reads.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace panvar {

struct PanGenieOptions {
    // src/transitionprobabilitycomputer.cpp: distance = (to - from) * 0.000004 * recombrate * Ne
    double recombrate = 1.26;
    double effective_N = 25000.0;
    // src/probabilitytable.cpp: added to each copy-number probability before normalising, so a k-mer
    // whose count is impossible under every copy number cannot zero the whole product.
    double regularization = 0.001;
    std::size_t threads = 0;
};

// One genotype per block, in the same shape panvar's own caller reports, so the two are directly
// comparable. `undefined` marks a block PanGenie's marker rule left with no usable k-mer at all --
// their emission then returns a uniform 1.0 for every allele pair, which is not the same thing as a
// confident call and must not be scored as one.
struct PanGenieCall {
    std::size_t block_index = 0;
    std::size_t allele1 = 0;
    std::size_t allele2 = 0;
    double gq = 0.0;
    std::size_t n_markers = 0;
    bool undefined = false;
};

// Their marker rule, from src/uniquekmercomputer.cpp:
//   occurences[kmer] gets the allele only when `entry.second == 1`  (occurs exactly once in it)
//   `(genomic_count - local_count) != 0` -> skip                    (not unique genome-wide)
//   `local_count > 1` -> skip                                       (carried by more than one allele)
// and UniqueKmers stores presence, not count (`bool kmer_on_allele`), so an expected count is
// bool + bool and cannot exceed 2. Genome-wide uniqueness has no analogue here -- we only have the
// locus -- so it is approximated by panvar's own region confinement, which is noted where it matters.
std::vector<PanGenieCall> genotype_pangenie(
    const std::vector<Block>& chain,
    const std::vector<BlockAlleles>& blocks,
    const ReadPanel& panel,
    const ReadCounts& counts,
    const std::vector<BlockDepth>& depth,
    const std::vector<std::string>& haplotype_names,
    const PanGenieOptions& options);

} // namespace panvar
