#pragma once

#include <string>
#include <vector>

namespace panvar {

// `panvar associate`: GWAS on a BIMBAM dosage genotype matrix (from `describe`) against a
// phenotype/covariate table. Per-feature GLM (linear or logistic), MAF filter on the final
// genotypes, and region-appropriate multiple-testing correction (Bonferroni + Benjamini-Hochberg).
int run_associate_command(const std::vector<std::string>& args);

} // namespace panvar
