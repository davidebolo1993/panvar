#pragma once

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
double log_unordered_diplotype_prior(
    const PanelAlleleMatrix& panel,
    const MosaicDiplotype& diplotype,
    const std::vector<double>& switch_probability);

} // namespace panvar
