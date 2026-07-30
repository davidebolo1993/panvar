#pragma once

#include "panvar/genotype_blocks.hpp"
#include "panvar/genotype_reads.hpp"

#include <string>
#include <vector>

namespace panvar {

// Everything needed to genotype a sample that depends only on the PANEL, not on the reads: the block
// chain, each block's alleles and which haplotype carries which, and the marker panel with its depth
// anchors. Building it is the expensive part of a run (spelling every path, extracting syncmers,
// enforcing region uniqueness) and it is identical for every sample, so a cohort should build it once
// and reuse it.
//
// Deliberately excluded: `BlockAlleles::allele_seq`. The spelled sequence of every allele is needed
// while markers are being extracted and to resolve leave-one-out truth, but not to genotype -- so it
// is neither stored nor reloaded.
struct GenotypeIndex {
    std::vector<Block> chain;
    std::vector<BlockAlleles> blocks;          // without allele_seq
    ReadPanel panel;
    std::vector<std::string> haplotype_names;
    std::size_t kmer_size = 31;
    std::size_t syncmer_s = 0;
};

// Binary, versioned. Not portable across endianness -- an index is a cache, rebuild it if in doubt.
void write_genotype_index(const std::string& path, const GenotypeIndex& index);
GenotypeIndex read_genotype_index(const std::string& path);

} // namespace panvar
