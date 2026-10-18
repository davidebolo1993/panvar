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
    std::string hap_a_sha256, hap_b_sha256;   // A: the group's sequence with the smaller SHA-256 (FASTA hap1 if representative)
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

// S1-D (sequence_first_v2): a predicted haplotype's membership in the ACTIVE (curated) panel, by exact sequence identity
// in the locus frame. A panel haplotype is identified by every active template that spells it (sequence-equivalent
// aliases share the match), each with the orientation of its source path relative to the locus chain.
struct HaplotypeOrigin {
    bool evaluated = false;
    bool panel = false;
    std::vector<std::string> paths;          // sorted template names spelling exactly this sequence
    std::vector<std::string> orientations;   // per path: forward | reverse_complement (source path vs locus chain)
};

struct SequenceReport {
    bool built = false;
    std::string schema = "sequence_first_v1";   // sequence_first_v1 | sequence_first_v2 (S1-D)
    std::string prediction = "NO_PREDICTION";   // PREDICTION | PROVISIONAL | NO_PREDICTION
    // RESOLVED | GENOTYPE_RESOLVED_SEQUENCE_SET | SEQUENCE_TIED (PREDICTION); NOT_DETERMINED (PROVISIONAL); NA
    std::string sequence_status = "NA";
    std::string completion = "NA";              // completed_ranking | round_cap | cumulative_cap | the stop reason
    std::string completion_note;
    std::vector<std::string> block_ids;
    std::size_t representative_group = 0;
    std::size_t representative_member = 0;      // index into PooledDecision::candidates (PREDICTION only)
    // one explored labelling: representative_first spells hap1_bases, representative_second spells hap2_bases (S1-C)
    std::vector<std::size_t> representative_first, representative_second;
    std::string hap1_bases, hap2_bases;         // the representative sequence pair, ordered by SHA-256
    std::string hap1_sha256, hap2_sha256;
    // canonical_sha256 (several decision-set groups, all numerically tied) or highest_score (one group)
    std::string representative_tie_break = "highest_score";
    std::size_t exact_score_ties = 1;   // decision-set groups with a bitwise-equal best score (information only)
    std::size_t decision_set_groups = 0, pool_groups = 0, explored_labellings = 0;
    double representative_probability = std::numeric_limits<double>::quiet_NaN();   // candidate-conditional
    double decision_set_probability = std::numeric_limits<double>::quiet_NaN();     // candidate-conditional
    double representative_score = std::numeric_limits<double>::quiet_NaN();
    std::string legacy_phase_support = "NA";
    std::vector<SequenceBlockRow> blocks;
    std::vector<RepresentationInterval> intervals;
    std::vector<SequenceAlternative> alternatives;
    bool publication_failed = false;   // S1-C: writing the report failed; nothing partial is left
    std::string publication_error;
    // S1-D (sequence_first_v2 only)
    HaplotypeOrigin origin[2];                  // hap1, hap2
    std::string graph_sha256 = "NA", model_fingerprint = "NA";
};

// S1-D: classifies hap1/hap2 against the active panel (template names, their spelled locus-frame sequences, and
// whether each source path runs opposite to the locus chain). Lookup by SHA-256, membership by exact string equality.
void annotate_panel_membership(SequenceReport& report, const std::vector<std::string>& template_names,
                               const std::vector<std::string>& template_sequences,
                               const std::vector<bool>& template_reverse);
// S1-D: whether sequence_first_v2 writes the FASTA, and why: OFF_PANEL_PAIR (written), PANEL_REPRESENTED_PAIR or
// NO_PREDICTION (not written). sequence_first_v1 always writes it for a prediction.
std::pair<bool, std::string> sequence_fasta_decision(const SequenceReport& report);

// Spells one haplotype from its allele labels (the caller's own spelling function).
using HaplotypeSpeller = std::function<std::string(const std::vector<std::size_t>&)>;

// After a completed, valid ranking. Deterministic and order-independent: the representative is the canonical member of
// the numerically inseparable decision set (every decision-set group is tied with the best under the caller's
// certified bounds): the smallest sorted pair of SHA-256 digests of the spelled sequences. It is a representative, not
// necessarily the strictly highest floating-point score. Scores are written with 12
// significant digits (as mosaic_scores.tsv): a prior's last bits depend on the panel's sample order.
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
    std::string dosage_mode = "NA";   // the run's --dosage-mode (the report is the caller decision BEFORE that stage)
    std::string dosage_mode_source = "NA";   // S1-D: flag | contract default (sequence_first_v2)
};
void write_sequence_report(const SequenceReport& report, const SequenceLegacy& legacy, const std::string& prefix,
                           const std::string& no_prediction_reason);
// Publishes the report (S1-C): content files first, the locus commit record last. On any failure every owned file of
// this prefix is removed, a NO_PREDICTION record naming the failure is attempted, `report` is replaced by a failed
// report (exit 2) and false is returned.
bool publish_sequence_report(SequenceReport& report, const SequenceLegacy& legacy, const std::string& prefix,
                             const std::string& no_prediction_reason);
// The exact suffixes a run owns under its prefix for this report (files and their atomic_write staging names).
const std::vector<std::string>& sequence_report_outputs();
// The status keys the caller adds under sequence_first_v1.
std::string sequence_report_status_keys(const SequenceReport& report);
// Exit semantics under sequence_first_v1: 0 PREDICTION, 3 PROVISIONAL, 1 NO_PREDICTION, 2 publication failed.
int sequence_report_exit_code(const SequenceReport& report);

}  // namespace panvar
