#pragma once

// M9-J (opt-in; profile keys compound_annotation off|contract_v2 and compound_interval_declarations <file>): compound
// interval annotations beside sequence_blocks, implementing the M9-I.1 contract v2. An ANNOTATION only: predictions,
// candidates, scores, priors, search, block rows, phase assertions, GQ, filters, FASTA policy and dosage are unchanged.
//
// Target of an interval [first..last] of consecutive blocks: the UNORDERED pair of concatenated interval sequences of a
// labelling -- exact sequence (no SNP or indel tolerance), including within-interval phase, invariant to homologue exchange
// and to synonymous block labellings; not structural equivalence.
//
// Over the explored pool (every sequence group, every recorded labelling of it, each group's mass counted once):
//   INVALID      a non-finite or negative mass: no bounds
//   UNASSESSED   no group, or total mass 0
//   masses are normalised by their total (compensated summation, group order); the discrepancy total - 1 is reported
//   support_lower  groups whose projections all equal the target
//   support_upper  groups with any projection equal to the target, plus groups without recorded labellings (missing)
//   error_lower    groups with labellings and no projection equal to the target
//   error_upper    groups with any projection unequal to the target, plus missing groups
//   every bound is a direct compensated sum (never 1 - something); zero-mass groups are kept and counted
//   competitor (from projections, not mass): NO_EVALUATED_ALTERNATIVE | EVALUATED_ALTERNATIVE |
//              EVALUATED_ALTERNATIVE_ZERO_MASS | UNKNOWN (missing groups, no recorded difference); never evidence of
//              correctness
//   target resolution (M9-J.1: decided independently of the mass bounds, over the decision-set groups only):
//              NO_DECISION_SET          no group is flagged as a decision-set member
//              INCOMPLETE_DECISION_SET  some decision-set group has no recorded labelling (its projection is unknown)
//              UNRESOLVED_DECISION_SET  every decision-set group is projected and the projections differ
//              RESOLVED                 every decision-set group is projected and all its labellings give one target
//              only RESOLVED has mass_reference the_resolved_target; otherwise the bounds refer to the representative's
//              projection only and mass_reference says why. Missing NON-decision groups never affect resolution.
// Scope: explored candidates and recorded labellings only; pool-conditional and uncalibrated; nothing about unsearched
// alternatives. Intervals are conservative declarations, not minimal or exhaustive.
#include <cstddef>
#include <functional>
#include <string>
#include <utility>
#include <vector>

namespace panvar {

// ---- declarations: a model-bound artifact generated from the model inventory
struct CompoundDeclaration {
    std::string id;
    std::vector<std::string> block_ids;   // consecutive in model order
    std::string source;
};
struct CompoundDeclarations {
    std::string path = "none", sha256 = "NA", model_fingerprint;
    std::vector<std::string> model_block_ids;   // the declaring model's ordered block ids
    std::vector<CompoundDeclaration> intervals;
};
// Throws std::runtime_error on a malformed file.
CompoundDeclarations read_compound_declarations(const std::string& path);

struct CompoundSpan {
    std::size_t first = 0, last = 0;
    std::vector<std::string> sources;   // declaration ids and/or "representation_interval"
};
// Validates the declarations against the active model (fingerprint, ordered block ids, consecutive interval blocks) and
// coalesces them with the run's representation intervals (overlapping spans merge; adjacent spans do not).
// declaration_status: NONE | VALID | INVALID_FINGERPRINT | INVALID_BLOCK_IDS; an invalid file contributes no span.
std::vector<CompoundSpan> resolve_compound_spans(const CompoundDeclarations* declarations, const std::string& model_fingerprint,
                                                 const std::vector<std::string>& block_ids,
                                                 const std::vector<std::pair<std::size_t, std::size_t>>& representation_intervals,
                                                 std::string* declaration_status);

// ---- assessment
using CompoundAlleleBytes = std::function<const std::string&(std::size_t block, std::size_t allele)>;
struct CompoundGroup {
    double mass = 0.0;
    std::vector<std::pair<const std::vector<std::size_t>*, const std::vector<std::size_t>*>> members;   // every recorded labelling
    bool decision_set = false;
};
struct CompoundAnnotationRow {
    std::size_t first = 0, last = 0;
    std::string blocks, sources, status = "UNASSESSED", target_resolution = "NA", mass_reference = "NA", competitor = "NA";
    std::size_t decision_set_projections = 0;
    std::size_t decision_set_groups = 0, decision_set_missing_groups = 0;   // M9-J.1
    std::string decision_set_projection_digests = "NA";
    std::string target_a_sha256 = "NA", target_b_sha256 = "NA";
    std::size_t target_a_bp = 0, target_b_bp = 0;
    double support_lower = 0.0, support_upper = 0.0, error_lower = 0.0, error_upper = 0.0;
    double original_total = 0.0, normalisation_discrepancy = 0.0;
    std::size_t groups = 0, missing_groups = 0, ambiguous_groups = 0, zero_mass_groups = 0;
    bool bounds = false;   // false: no bounds are reported (INVALID, UNASSESSED, not assessed)
};
CompoundAnnotationRow compound_assess(const CompoundAlleleBytes& bytes, std::size_t first, std::size_t last,
                                      const std::vector<CompoundGroup>& groups, const std::vector<std::size_t>& representative_first,
                                      const std::vector<std::size_t>& representative_second);

// A run's annotation table (one row per span). `status_override` rows (provisional / no prediction) carry no bounds.
struct CompoundAnnotation {
    bool enabled = false;
    std::string declaration_path = "none", declaration_sha256 = "NA", declaration_status = "NONE";
    std::vector<CompoundAnnotationRow> rows;
    double seconds = 0.0;
};
std::string render_compound_annotation(const CompoundAnnotation& annotation, const std::vector<std::string>& block_ids);

}  // namespace panvar
