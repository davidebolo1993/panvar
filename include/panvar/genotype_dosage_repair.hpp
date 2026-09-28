#pragma once

// Dosage repair (docs/algorithms/genotype_dosage.md, section 8): for an AUTO_REPAIR block, dosage-admissible
// candidates are constructed (pair route: panel allele pairs; repeated-unit route: unit multisets screened by the
// declared count model), spliced into the reported homologues, scored exactly by the production likelihood, and
// chosen by the top-genotype rule and, within one content multiset, the panel-structure allocation rule. Shadow only:
// nothing here changes the reported genotype.

#include "panvar/fragment_likelihood.hpp"
#include "panvar/genotype_dosage.hpp"

#include <cstddef>
#include <cstdint>
#include <map>
#include <optional>
#include <string>
#include <vector>

namespace panvar::dosage {

struct RepairSettings {
    std::size_t pair_cap = 4096;
    std::size_t unit_pool = 36;
    std::size_t cap_multisets = 64, cap_haploid = 256, cap_diploid = 4096;
    std::size_t em_iterations = 3000;
    double em_tolerance = 1e-10;
    std::size_t large_optional_run_bp = 1000;
    double tie_tolerance_nats = 1e-9;       // CertificationOptions::minimum_tie_tolerance_nats
    double guard_tolerance = 1e-6;
    std::size_t chunk_haplotypes = 128;
};

// ---------------------------------------------------------------------------------------------- unit library

struct Unit {
    std::size_t id = 0;
    std::string sequence;
    char state = 'O';                       // 'I': carries an optional run >= large_optional_run_bp
    bool internal_capable = false;          // last node -> entry node is a panel edge
    bool terminal_capable = false;          // last node -> right flank is a panel edge
    bool anchor_compatible = false;         // leaves the bases trimmed into the left anchor unchanged
    std::vector<std::vector<std::string>> walks;             // distinct unit walks (steps as "<node><+|->")
    std::vector<std::array<std::size_t, 3>> carriers;        // (allele id, position 1-based, allele unit count)
};

struct UnitLibrary {
    std::string entry_step;
    std::string flank_left, flank_right;    // flank_left includes the pre-entry segment
    std::size_t trim_left = 0, trim_right = 0;
    std::vector<Unit> units;
    std::vector<std::vector<std::size_t>> panel;   // per allele id: its unit tuple
    std::size_t kmax = 0;
    std::string block_sequence(const std::vector<std::size_t>& tuple) const;
};

UnitLibrary build_unit_library(const PreparedGenotypeBlock& block, const NodeSequences& nodes, const RepairSettings& settings);

// ---------------------------------------------------------------------------------------------- scoring

struct PersistedScore {
    double total = 0.0, fragment = 0.0, exposure = 0.0;
};

struct ScoringInputs {
    const PreparedGenotypeModel* model = nullptr;
    MosaicPath first, second;               // the reported call's allele paths
    const ReportedCall* call = nullptr;
    const std::vector<ReadPair>* reads = nullptr;
    const LibraryModel* library = nullptr;
    FragmentScoringOptions options;
    std::vector<double> switch_probability;
    PersistedScore persisted;
};

struct BlockCandidate {
    std::string label;
    std::string orientation = "AB";
    std::string ctx1_sequence, ctx2_sequence;          // block sequences placed in the first / second homologue
    std::string ctx1_id, ctx2_id;                      // panel allele id, or the full SHA-256 of a generated sequence
    std::size_t rank = 0;                              // construction order (structure rank or pair order)
    // pair route
    std::size_t allele_a = 0, allele_b = 0;
    // repeated-unit route
    std::vector<std::size_t> units_a, units_b, multiset;
    std::size_t multiset_rank = 0;
    double proposal_loglik = 0.0;
    std::string states;
    double panel_log_prior = 0.0;
    // exact scores
    double fragment_term = 0.0, exposure_term = 0.0, likelihood = 0.0;
    bool in_top_set = false, chosen_content = false, chosen = false, allocation_tied = false;
};

struct RepairResult {
    std::string route;                      // pair / repeated_unit
    std::string outcome = "REPAIRED";       // REPAIRED or TRIGGERED_UNREPAIRED
    std::string reason = "none";
    std::string scope = "COMPLETE";
    std::string scope_reason = "none";
    std::string cap_reached = "none";
    std::vector<BlockCandidate> candidates;
    long chosen = -1;
    std::size_t top_set_size = 0;
    bool orientation_tie = false;
    std::size_t allocation_tie = 0;         // distinct genotypes in the allocation-equivalence set (prior, then likelihood within tolerance)
    bool representative_only = false;      // the chosen arrangement only serialises an unresolved allocation-equivalence set
    std::size_t multisets_evaluated = 0, multisets_retained = 0, haploid = 0, diploid = 0;
    std::vector<int> m_range;
    std::size_t internal_capable = 0, terminal_capable = 0;
    double beta = 0.0, phi = 0.0;
    std::size_t em_iterations = 0;
    double original_likelihood = 0.0, original_fragment = 0.0, original_exposure = 0.0;
    double repair_separation = 0.0;
    double conditional_repair_gap = 0.0;
    bool gap_defined = false;
    double seconds_construct = 0.0, seconds_score = 0.0;
};

RepairResult repair_block(const BlockAudit& audit, const CensusBlock& census, const BlockFeatureModel& features,
                          const NodeSequences& nodes, const ScoringInputs& scoring, const DosageSettings& settings,
                          const RepairSettings& repair_settings);

// The choice among exactly scored candidates (spec section 8): the top-genotype rule, orientation ties, and on the
// repeated-unit route the panel-structure allocation within one content multiset. Sets outcome/reason on a tie.
void choose_repair(RepairResult& result, CensusRoute route, double tie_tolerance_nats);

std::string repair_candidates_table(const std::vector<std::pair<std::string, RepairResult>>& repairs);
std::string repair_blocks_table(const std::vector<std::pair<std::string, RepairResult>>& repairs);

std::string sha256_hex(const std::string& bytes);

// ---------------------------------------------------------------------------------------------- provenance

// A panel sequence the provenance may compare against: only the panel supplied to the caller.
struct PanelSequence {
    std::string id;          // template name (full haplotypes) or block allele id
    std::string sequence;
};

struct NearestPanel {
    std::string id;
    std::string sha256;
    long distance = 0;                  // global edit distance, unit costs
    std::size_t alignment_columns = 0;  // '=' + 'X' + 'I' + 'D' columns of the reported alignment
    double normalized = 0.0;            // distance / alignment_columns
    std::string cigar;                  // extended CIGAR, query = the repaired sequence, target = the panel sequence
};

// Every panel sequence at the minimal global edit distance from `query`, ordered by (SHA-256 of the sequence, id):
// invariant to the panel's order. Ties are all reported.
std::vector<NearestPanel> nearest_panel_sequences(const std::string& query, const std::vector<PanelSequence>& panel);

// Registers a generated sequence's display id (the first `length` hex digits of its full SHA-256). False on a
// collision: two distinct full ids sharing one display id, which fails the run closed.
bool register_display_id(std::map<std::string, std::string>& shown, const std::string& full_id, std::size_t length = 12);

// SHA-256 of the sorted (id, sequence SHA-256) list: what was searched, independent of its order.
std::string panel_digest(const std::vector<PanelSequence>& panel);

// ---------------------------------------------------------------------------------------------- allocation test

// Truth-free allocation observability (spec section 8.3; prototype sv7_selfreview/alloc_test.py): every alternative
// allocation of the reported unit copies is rebuilt as a junction-preserving transposition and scored exactly; an
// identical per-fragment contribution vector (1e-9-nat grid) and exposure (1e-9) means no fragment separates the two
// allocations. An identifiability annotation only: it never changes a call or triggers repair.
struct AllocationRow {
    std::string arrangement;     // "<hap1 copy indices>|<hap2 copy indices>"
    int flank_order = 0;
    std::string shape;
    std::uint64_t fragment_fnv64_1e9 = 0;
    std::size_t fragments_differing = 0;
    double exposure_delta = 0.0;
    double likelihood_delta = 0.0;
    bool identical = false;
};
struct AllocationTestResult {
    std::string status = "NOT_APPLICABLE";   // RESOLVED, UNRESOLVED, INCONCLUSIVE, NOT_APPLICABLE
    long switch_coordinate = -1;
    std::size_t copies = 0;
    std::uint64_t reference_fnv64_1e9 = 0;
    double reference_likelihood = 0.0;
    std::vector<AllocationRow> rows;
    double seconds = 0.0;
};
AllocationTestResult allocation_observability_test(const UnitLibrary& lib, const std::vector<std::size_t>& hap1_units,
                                                   const std::vector<std::size_t>& hap2_units, std::size_t block,
                                                   const ScoringInputs& scoring, const RepairSettings& settings);

// ---------------------------------------------------------------------------------------------- shadow audit

// The call a finished run reports: the equivalence set's rank-1 representative, with its persisted scores.
struct PersistedCall {
    MosaicPath first, second;
    std::string hap1, hap2;
    PersistedScore score;
    double prior = 0.0;
};
std::optional<PersistedCall> read_persisted_call(const std::string& run_prefix);

struct DosageShadowOptions {
    bool repair = true;
    bool provenance = true;
    bool allocation_test = false;   // shadow-only, off until its measured cost is reviewed
    std::size_t threads = 1;
    std::size_t chunk_haplotypes = 128;
};

struct DosageShadowSummary {
    std::string audit_status = "COMPLETE";   // COMPLETE, NOT_RUN (no persisted call), REFUSED_LIBRARY
    std::string note;
    std::size_t molecule_span_bp = 0, m = 0, sensitive_blocks = 0, repairs = 0, repairs_conditional = 0, allocation_blocks_tested = 0;
    std::map<std::string, std::size_t> status_counts, reason_counts;
    double seconds_census = 0.0, seconds_features = 0.0, seconds_estimate = 0.0, seconds_repair = 0.0, seconds_provenance = 0.0, seconds_allocation = 0.0,
           seconds_total = 0.0;
};

// Runs census, features, estimate, classification and (optionally) repair of the call persisted at run_prefix, and
// writes out_prefix.dosage_{census,features,blocks,repairs,repair_candidates,status}.tsv. Writes no other file and
// changes no call. Truth-free.
DosageShadowSummary run_dosage_shadow(const PreparedGenotypeModel& model, const NodeSequences& nodes,
                                      const std::vector<ReadPair>& reads, const LibraryModel& library,
                                      const FragmentScoringOptions& scoring, const std::vector<double>& switch_probability,
                                      const std::string& run_prefix, const std::string& out_prefix,
                                      const DosageShadowOptions& options,
                                      const std::optional<PersistedCall>& call_override = std::nullopt);

} // namespace panvar::dosage
