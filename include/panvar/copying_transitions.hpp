#pragma once

#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

#include "panvar/genotype_locus.hpp"

namespace panvar {

// M4-B: the ONE construction of the copying prior's per-boundary redraw probabilities, shared by every supported
// stage (proposals, ranking, refinement, projection, dosage competitors and audits) and by the diagnostic score-set.
//
//   T_b = (1 - r_b) I + r_b U,   U = uniform over the n template rows (the current one included; duplicate rows count)
//
//   constant_v1                         r_b = constant_r for every boundary (the default; unchanged behaviour)
//   distance_ref_core_v1_experimental   r_b = 1 - exp(-d_b / lambda), d_b = distance between consecutive reference
//                                       core midpoints (LocusGeometry, computed against the frame), lambda = a declared
//                                       copying length in bp. Computed as -expm1(-d/lambda): exact for small d,
//                                       saturating to 1 (T -> U) for d >> lambda. d_b = 0 gives r_b = 0 (T = I): two
//                                       informative states at one reference position cannot recombine under this model.
//
// r_b is a copying-model donor-redraw probability, not a biological crossover probability. The sequence-level
// recursions (log_haploid_sequence_prior and friends) consume the vector unchanged: every allele path spelling the
// bytes is still summed.
struct CopyingTransitionSpec {
    std::string mode = "constant_v1";   // constant_v1 | distance_ref_core_v1_experimental
    double constant_r = 0.01;
    double copying_length_bp = 0.0;     // distance mode only
};

// A requested specification the inputs cannot honour (missing/unavailable geometry, invalid parameter, size mismatch).
// Never caught to fall back to constant transitions.
class UnsupportedCopyingTransition : public std::invalid_argument {
public:
    using std::invalid_argument::invalid_argument;
};

struct CopyingTransitions {
    std::string mode;
    std::vector<double> r;              // one per block boundary
    std::vector<double> distance_bp;    // distance mode only
    double constant_r = 0.0;
    double copying_length_bp = 0.0;
    std::string geometry_digest = "none";
    std::string vector_digest;          // SHA-256 over the exact bits of r (and the mode)
    std::size_t zero_distance_boundaries = 0;              // distance mode: d_b == 0
    std::size_t zero_distance_informative_boundaries = 0;  // ... between two informative (polymorphic-frame) states
    bool distance_mode() const noexcept { return mode == "distance_ref_core_v1_experimental"; }
};

CopyingTransitions build_copying_transitions(std::size_t blocks, const LocusGeometry* geometry, const CopyingTransitionSpec& spec);

// The hazard-matched DIAGNOSTIC control length: the lambda giving the same total hazard as constant r over this
// geometry's boundaries, sum_b d_b / ((B-1) * -log(1-r)). Depends on block count by construction; a control only.
double hazard_matched_copying_length(const LocusGeometry& geometry, double constant_r);

}  // namespace panvar
