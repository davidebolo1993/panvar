#pragma once

// panvar genotype --dosage-mode off|audit|repair (docs/algorithms/genotype_dosage.md sections 20 and 22) and the
// primary locus-level output, <prefix>.locus_genotype.tsv. Everything here runs after the caller's outputs are
// committed and reads them back from disk, so the same code integrates a fresh run and replays a sealed one.

#include "panvar/genotype_dosage_repair.hpp"
#include "panvar/marker_model.hpp"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <optional>
#include <string>
#include <vector>

namespace panvar::dosage {

// The caller's call-bearing outputs. In repair mode each one that exists is copied byte for byte to
// <prefix>.original<suffix> before anything is rewritten, so <prefix>.original is a complete legacy prefix the
// evaluator scores as it scores any run.
const std::vector<std::string>& legacy_call_outputs();

// Every file the dosage stage or the locus output can leave under a prefix: a run removes them before it starts.
const std::vector<std::string>& dosage_stage_outputs();

struct DosageIntegration {
    std::string mode = "off";                  // off | audit | repair
    std::string audit_status = "NOT_RUN";       // the stage's own status (COMPLETE, NOT_RUN, REFUSED_LIBRARY, FAILED)
    std::string audit_note = "none";
    const DosageShadowSummary* summary = nullptr;
    // repair mode only: the COMPLETE_DIPLOTYPE_REPAIR blocks that replaced the reported genotype, and those that did
    // not, with the fail-closed reason
    std::vector<std::string> applied;
    std::map<std::string, std::string> not_applied;
    // audit mode only, a record for the headline filters: each COMPLETE_DIPLOTYPE_REPAIR block with the guard reason
    // repair mode would refuse it for, or headline::kAuditModeNotApplied when repair mode would apply it
    std::map<std::string, std::string> audit_not_applied;
    std::optional<MosaicPath> reported_first, reported_second;   // set only when a repair was applied
    double reported_fragment = 0.0, reported_exposure = 0.0, reported_prior = 0.0, reported_separation = 0.0;
    bool originals_written = false;
    // basis of the reported filters: caller, recomputed_for_reported_genotype, or not_evaluated_after_dosage_repair
    std::string filter_basis = "caller";
    std::string filter_note = "none";
};

// What re-annotating a reported repaired genotype needs (section 23): the reads, the library and its declared marker
// background, the marker settings, and the scoring settings with the repair caps. Null: the reported genotype's
// filters are NOT_EVALUATED_AFTER_DOSAGE_REPAIR.
struct ReportedFilterInputs {
    const std::vector<ReadPair>* reads = nullptr;
    const LibraryModel* library = nullptr;
    std::optional<double> declared_marker_background;
    std::size_t kmer_size = 31, syncmer_s = 15;
    FragmentScoringOptions scoring;
    std::uint64_t repair_max_physical_origins = 0, repair_max_state_classes = 0, repair_max_piece_hits = 0, repair_max_latent_pairs = 0;
    double switch_probability = 0.01;
    bool representation_guard = true;
};

// One block's representation_conflict_rule_v1 result, as the caller computes it for its own call.
struct GuardRow {
    std::string filter = "PASS", reason = "none", competitor = "NA";
    double marker_margin = 0.0;
    double fragment_diff = std::numeric_limits<double>::quiet_NaN(), exposure_diff = std::numeric_limits<double>::quiet_NaN();
    double likelihood_margin = std::numeric_limits<double>::quiet_NaN(), statistic = std::numeric_limits<double>::quiet_NaN();
};
// The caller's marker emissions, rebuilt from the model, the reads and the declared background exactly as panvar
// genotype builds them (marker depth fit, proposal parameters, emissions).
PreparedProposalEmissions rebuild_marker_emissions(const PreparedGenotypeModel& model, const ReportedFilterInputs& in);
// The guard over every block of the diplotype (first, second) whose exact score is `called` (prior `called_prior`):
// the marker's best pair, its margin, and where decisive the substituted diplotype scored exactly in both assignments.
std::vector<GuardRow> representation_guard_rows(const PreparedGenotypeModel& model, const PreparedProposalEmissions& emissions,
                                                const ReportedFilterInputs& in, const MosaicPath& first, const MosaicPath& second,
                                                const DiplotypeScore& called, double called_prior);
// The caller's per-block read-support rule (genotype_repair.cpp): a block is READS-supported unless some rival whose
// fragment term is not below the reported diplotype's, within both numeric error bounds, carries another genotype there.
struct Rival {
    double fragment_term = 0.0, error_bound = 0.0;
    std::vector<bool> differs;   // per block: its unordered genotype differs from the reported one
};
std::vector<bool> read_support(double fragment_term, double error_bound, const std::vector<Rival>& rivals, std::size_t blocks);
// The caller's ranked sequence candidates (mosaic_scores.tsv with their representatives' paths from
// mosaic_candidates.tsv), in ranking order: fragment term, bound and unordered genotype per block.
struct RankedCandidate {
    double fragment_term = 0.0, error_bound = 0.0;
    bool winner = false;
    std::vector<std::pair<std::size_t, std::size_t>> genotype;
};
std::vector<RankedCandidate> read_caller_ranking(const std::string& prefix);
// The caller's filter columns for one block, formatted as the caller writes them (12 significant digits, NA for NaN;
// the marker margin NA only when the guard is off).
std::vector<std::string> guard_cells(const GuardRow& row);
// The caller's locus-level filter and quality summary (status keys, in the caller's order) from per-block rows, and the
// library-sensitive blocks of the reported genotype (library_sensitivity_audit_v1), which add their locus flag.
std::vector<std::pair<std::string, std::string>> locus_filter_summary(const std::vector<std::string>& block_ids,
                                                                      const std::vector<std::string>& filters,
                                                                      const std::vector<std::string>& supports,
                                                                      const std::vector<std::string>& qualities,
                                                                      const std::vector<std::string>& library_sensitive = {},
                                                                      bool library_not_evaluated = false);

// off / audit: copies the caller's outputs from run_prefix to out_prefix when they differ, and changes nothing.
// repair: additionally writes <out>.original<suffix>, then applies every eligible COMPLETE_DIPLOTYPE_REPAIR to the
// block tables and the haplotype FASTA, appends the dosage columns to both block tables and the dosage keys to the
// status. Every new file content is built before the first rewrite is written: a refusal leaves no partial repair.
DosageIntegration apply_dosage_mode(const PreparedGenotypeModel* model, const std::string& mode,
                                    const std::string& run_prefix, const std::string& out_prefix,
                                    const DosageShadowSummary* summary, const std::string& audit_status,
                                    const std::string& audit_note, const ReportedFilterInputs* filters = nullptr);

// <out>.locus_genotype.tsv: one header row and one row for the locus (schema in section 22.3). Written in every
// mode; it reads the committed outputs under out_prefix. model may be null (switches are then NA).
void write_locus_genotype(const PreparedGenotypeModel* model, const std::string& out_prefix, const DosageIntegration& d);

// The minimum number of panel-template switches along one haplotype path, and the blocks where a new template
// segment starts under the canonical left-greedy decomposition (each segment extended while any template carries
// every allele so far). nullopt when some allele is carried by no template.
struct TemplateSwitches {
    std::size_t count = 0;
    std::vector<std::size_t> blocks;
};
std::optional<TemplateSwitches> template_switches(const PreparedGenotypeModel& model, const MosaicPath& path);

// The locus table's column names, in order (the schema gate reads this).
const std::vector<std::string>& locus_genotype_columns();
// The trailing columns repair mode appends to genotypes.tsv and block_genotypes.tsv, in order.
const std::vector<std::string>& dosage_block_columns();

} // namespace panvar::dosage
