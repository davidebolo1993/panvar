#pragma once

// M3-E (EXPERIMENTAL, default off; profile key block_support_rule full_data_likelihood_v1_experimental): block-genotype
// support under the FULL DATA LIKELIHOOD (fragment + exposure, the copying prior excluded), kept separate from phase.
// One shared implementation for the pooled and the unpooled summary paths. Contract: tests/analysis/m3e/DESIGN.md.
//
// Candidate unit: a sequence-level group (members spelling one unordered sequence pair). A group's genotype at a block is
// the SET of its members' unordered allele pairs there; a group whose members disagree at a block is projection-ambiguous
// there and is never resolved by choosing a representative. Per block b, with the called group w:
//   w projection-ambiguous at b                                   -> AMBIGUOUS_PROJECTION
//   g = w's genotype at b; same = groups all of whose members carry g; disagreeing = groups none of whose members do;
//   mixed = the rest. s / d = the same / disagreeing group with the largest likelihood L (ties: lowest index).
//   no disagreeing group                                          -> NO_EVALUATED_ALTERNATIVE (vacuous)
//   L_s - u_s > L_d + u_d and every mixed group certifiably below s -> DATA_SUPPORTED
//   L_d - u_d > L_s + u_s                                         -> PRIOR_DEPENDENT (the data prefer another genotype)
//   a mixed group not certifiably below s (d not above)           -> AMBIGUOUS_COMPETITOR
//   otherwise                                                     -> DATA_TIED
// Support conditional on the evaluated candidate set: not a calibrated probability, not a confidence, not phasing.
// Phase is NOT relabelled: the existing fragment contenders are only classified (genotype-changing / -preserving / mixed).

#include <cstddef>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace panvar {

using BlockGenotypeRow = std::vector<std::pair<std::size_t, std::size_t>>;   // unordered allele pair per block

struct SupportGroup {
    std::vector<BlockGenotypeRow> rows;   // distinct unordered genotype rows of the group's members
    double fragment = 0.0, exposure = 0.0, prior = 0.0, bound = 0.0;   // bound: the certification numeric_error_bound
};

struct BlockDataSupport {
    std::string label = "NOT_EVALUATED";
    double margin = std::numeric_limits<double>::quiet_NaN();   // L_s - L_d
    bool certified = false;
    long competitor = -1;                                       // d, an index into the groups; -1 when none
    double delta_fragment = std::numeric_limits<double>::quiet_NaN(), delta_exposure = std::numeric_limits<double>::quiet_NaN(),
           delta_likelihood = std::numeric_limits<double>::quiet_NaN(), delta_prior = std::numeric_limits<double>::quiet_NaN();
    std::size_t disagreeing = 0, mixed = 0;
};

struct PhaseCompetitorClasses {
    std::size_t genotype_changing = 0, genotype_preserving = 0, mixed = 0;
};

struct DataSupportResult {
    std::vector<BlockDataSupport> blocks;
    PhaseCompetitorClasses phase;
};

DataSupportResult block_data_support_v1(const std::vector<SupportGroup>& groups, std::size_t winner, std::size_t block_count);

// The genotype-support term the experimental quality filter adds for a label: "" (none) or its reason.
std::string data_support_quality_reason(const std::string& label);

}  // namespace panvar
