#pragma once

// panvar genotype --dosage-mode off|audit|repair (docs/algorithms/genotype_dosage.md sections 20 and 22) and the
// primary locus-level output, <prefix>.locus_genotype.tsv. Everything here runs after the caller's outputs are
// committed and reads them back from disk, so the same code integrates a fresh run and replays a sealed one.

#include "panvar/genotype_dosage_repair.hpp"

#include <cstddef>
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
    std::optional<MosaicPath> reported_first, reported_second;   // set only when a repair was applied
    double reported_fragment = 0.0, reported_exposure = 0.0, reported_prior = 0.0, reported_separation = 0.0;
    bool originals_written = false;
};

// off / audit: copies the caller's outputs from run_prefix to out_prefix when they differ, and changes nothing.
// repair: additionally writes <out>.original<suffix>, then applies every eligible COMPLETE_DIPLOTYPE_REPAIR to the
// block tables and the haplotype FASTA, appends the dosage columns to both block tables and the dosage keys to the
// status. Every new file content is built before the first rewrite is written: a refusal leaves no partial repair.
DosageIntegration apply_dosage_mode(const PreparedGenotypeModel* model, const std::string& mode,
                                    const std::string& run_prefix, const std::string& out_prefix,
                                    const DosageShadowSummary* summary, const std::string& audit_status,
                                    const std::string& audit_note);

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
