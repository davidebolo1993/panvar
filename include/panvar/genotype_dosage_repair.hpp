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
    std::size_t allocation_tie = 0;         // distinct genotypes tied on (panel prior, likelihood) at the pick
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

std::string repair_candidates_table(const std::vector<std::pair<std::string, RepairResult>>& repairs);
std::string repair_blocks_table(const std::vector<std::pair<std::string, RepairResult>>& repairs);

std::string sha256_hex(const std::string& bytes);

} // namespace panvar::dosage
