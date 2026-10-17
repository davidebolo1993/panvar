#pragma once
// Sequence-first reporting (output_schema sequence_first_v1; opt-in, experimental). A REPORTING layer over the pooled
// decision the caller has already made: it never changes a score, prior, candidate, ranking, search step, cap, legacy
// output or dosage/phase decision. Contract: tests/genotype_sequence_label_reporting_design.md sections 6-7 (S1-A),
// as corrected for S1-B (tests/analysis/s1b/README.md).
//
// Terms:
//   decision set   -- the numerical top set of the explored, exactly scored pool (PooledDecision::top_groups);
//   group          -- one unordered sequence pair; its members are the explored allele labellings that spell it;
//   RESOLVED       -- one group in the decision set: resolution WITHIN the declared decision set, never a claim of
//                     high confidence or biological correctness;
//   EXPLORED       -- block projections are those of the explored members only. Uniqueness and error-mass bounds are
//                     conditional on them; an unexplored labelling of the same sequence may exist.
#include "panvar/genotype_repair.hpp"

#include <cstddef>
#include <functional>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace panvar {

// One block of the sequence-first block table.
struct SequenceBlockRow {
    std::string block_id;
    // IDENTIFIED (one unordered genotype over every explored labelling of every decision-set group),
    // REPRESENTATION_AMBIGUOUS (one decision-set group, several explored labellings disagree here),
    // GENOTYPE_UNRESOLVED (several decision-set groups, genotypes disagree here),
    // PROVISIONAL_LABEL (a provisional incumbent's own labelling; no decision set, no confidence).
    std::string status;
    std::pair<std::size_t, std::size_t> genotype{0, 0};   // IDENTIFIED / PROVISIONAL_LABEL only
    // Observed genotypes over the decision set with the groups carrying each ("a/b:g1,g2;c/d:g3"); "NA" when unused.
    std::string observed_genotypes = "NA";
    // IDENTIFIED only: conditional error-mass bounds over the explored pool, each group counted once with its
    // candidate-set probability. lower: groups whose every explored labelling has another genotype here; upper: groups
    // with any explored labelling that has another genotype here. NaN otherwise.
    double error_mass_lower = std::numeric_limits<double>::quiet_NaN();
    double error_mass_upper = std::numeric_limits<double>::quiet_NaN();
    std::string interval = "none";   // the representation interval of the representative group holding this block
    // LEGACY DIAGNOSTIC: the existing representative-based block quality (one labelling per group), not the
    // projection-aware quantity above. NaN for a provisional report.
    double legacy_representative_gq = std::numeric_limits<double>::quiet_NaN();
    bool legacy_numerically_separated = false;
};

// A representation interval of one decision-set group: a conservative ENCLOSING block range (not claimed minimal)
// that covers every block where two explored labellings of the group differ, including unchanged blocks between them
// (coupled alternatives). Tuples are the whole observed labellings restricted to the range, per haplotype (A spells the
// group's first sequence by canonical digest order), never Cartesian combinations of per-block alternatives.
struct RepresentationInterval {
    std::size_t group = 0;
    std::size_t first_block = 0;
    std::size_t last_block = 0;
    std::vector<std::pair<std::vector<std::size_t>, std::vector<std::size_t>>> tuples;   // sorted, distinct
};

// One row of the alternatives table (the decision set, then the next groups by rank up to a fixed count). Sequences
// are identified by the caller's internal 128-bit digest and length; no sequence is copied.
struct SequenceAlternative {
    std::size_t rank = 0;   // 1-based rank in the exact ranking
    std::size_t group = 0;
    bool in_decision_set = false;
    bool representative = false;
    double score = 0.0, fragment_term = 0.0, exposure_term = 0.0, prior = 0.0, numeric_error_bound = 0.0;
    double candidate_set_probability = 0.0;
    std::size_t explored_labellings = 0;
    bool projection_ambiguous = false;
    SequenceDigest first, second;
    std::vector<std::size_t> blocks_with_other_genotypes;   // vs the representative group, over explored labellings
};

struct SequenceReport {
    bool built = false;
    std::string prediction = "NO_PREDICTION";   // PREDICTION | PROVISIONAL | NO_PREDICTION
    // RESOLVED | GENOTYPE_RESOLVED_SEQUENCE_SET | SEQUENCE_TIED (PREDICTION); NOT_DETERMINED (PROVISIONAL); NA
    std::string sequence_status = "NA";
    std::string completion = "NA";              // completed_ranking | round_cap | cumulative_cap | the stop reason
    std::string completion_note;
    std::vector<std::string> block_ids;
    std::size_t representative_group = 0;
    std::size_t representative_member = 0;      // index into PooledDecision::candidates (PREDICTION only)
    std::vector<std::size_t> representative_first, representative_second;   // one explored labelling, oriented A/B
    std::string hap1_bases, hap2_bases;         // the representative sequence pair, ordered by SHA-256
    std::string hap1_sha256, hap2_sha256;
    std::string representative_tie_break = "highest_score";   // or canonical_sha256 (exact score ties)
    std::size_t exact_score_ties = 1;
    std::size_t decision_set_groups = 0, pool_groups = 0, explored_labellings = 0;
    double representative_probability = std::numeric_limits<double>::quiet_NaN();   // candidate-conditional
    double decision_set_probability = std::numeric_limits<double>::quiet_NaN();     // candidate-conditional
    double representative_score = std::numeric_limits<double>::quiet_NaN();
    std::string legacy_phase_support = "NA";
    std::vector<SequenceBlockRow> blocks;
    std::vector<RepresentationInterval> intervals;
    std::vector<SequenceAlternative> alternatives;
};

// Spells one haplotype from its allele labels (the caller's own spelling function).
using HaplotypeSpeller = std::function<std::string(const std::vector<std::size_t>&)>;

// After a completed, valid ranking. Deterministic and order-independent: the representative is the highest-scoring
// decision-set group; exact score ties are broken by the sorted pair of SHA-256 digests of the spelled sequences.
SequenceReport build_sequence_report(const PooledDecision& decided, const std::vector<std::string>& block_ids,
                                     const HaplotypeSpeller& spell, std::size_t alternative_rows = 5);

// A declared-cap stop with a last committed, fully scored incumbent: PROVISIONAL, decision set NOT_DETERMINED, no
// confidence, no intervals, no alternatives.
SequenceReport provisional_sequence_report(const std::vector<std::size_t>& first, const std::vector<std::size_t>& second,
                                           double total_score, const std::string& stop, const std::string& note,
                                           const std::vector<std::string>& block_ids, const HaplotypeSpeller& spell);

// Outputs (written only under sequence_first_v1):
//   <prefix>.sequence_locus.tsv          key/value locus record (always; NO_PREDICTION records the reason)
//   <prefix>.sequence_blocks.tsv         one row per block (prediction only)
//   <prefix>.sequence_haplotypes.fa      the representative sequence pair, explicitly marked (prediction only)
//   <prefix>.representation_intervals.tsv, <prefix>.sequence_alternatives.tsv (PREDICTION only; header-only otherwise)
struct SequenceLegacy {
    std::string run_status = "NA", call_class = "NA", no_call_reason = "NA";
};
void write_sequence_report(const SequenceReport& report, const SequenceLegacy& legacy, const std::string& prefix,
                           const std::string& no_prediction_reason);
// The status keys the caller adds under sequence_first_v1.
std::string sequence_report_status_keys(const SequenceReport& report);
// Exit semantics under sequence_first_v1: 0 PREDICTION, 3 PROVISIONAL, 1 NO_PREDICTION.
int sequence_report_exit_code(const SequenceReport& report);

}  // namespace panvar
