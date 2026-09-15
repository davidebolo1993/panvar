#pragma once

#include <string>
#include <vector>

#include "panvar/genotype_model.hpp"

namespace panvar {

// Exact Li-Stephens prior of one fixed allele mosaic after marginalising over every compatible panel
// template history. `switch_probability[b]` is the probability at the boundary from block b to b+1
// in T = (1-r)I + (r/n)11^T. Complexity is O(blocks * templates).
double log_haploid_mosaic_prior(
    const PanelAlleleMatrix& panel,
    const MosaicPath& mosaic,
    const std::vector<double>& switch_probability);

// Prior of an unordered diploid mosaic. The two homologues evolve independently; two distinct paths
// therefore aggregate both orders and receive the exact log(2) multiplicity.
//
// This is correct only for SEQUENCE-CANONICAL paths. Two different allele paths can spell identical
// complete sequences, and no read-based method can tell them apart; for those, use the sequence
// priors below, which aggregate every representation and decide the log(2) on bytes, not on paths.
double log_unordered_diplotype_prior(
    const PanelAlleleMatrix& panel,
    const MosaicDiplotype& diplotype,
    const std::vector<double>& switch_probability);

// Exact prior of one complete haplotype sequence: the sum, over EVERY allele path that spells it,
// of that path's exact Li-Stephens prior. Computed as a forward recursion over (block, offset in the
// sequence, template), so no path is enumerated and the value does not depend on which candidates a
// proposal happened to retain. -inf when no allele path spells the sequence.
double log_haploid_sequence_prior(
    const PreparedGenotypeModel& model,
    const std::string& sequence,
    const std::vector<double>& switch_probability);

// Prior of an unordered pair of complete haplotype sequences. Both orders of two distinct sequences
// are the same unordered state, so they receive log(2); two identical sequences do not, however
// different the allele paths that produced them.
double log_unordered_sequence_pair_prior(
    const PreparedGenotypeModel& model,
    const std::string& first,
    const std::string& second,
    const std::vector<double>& switch_probability);

} // namespace panvar
