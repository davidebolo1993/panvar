#pragma once

// The three public headline filters (docs/algorithms/genotype_headline_filters.md; the frozen design is SPEC_v2 of the
// C4 headline-filter audit): GENOTYPE_FILTER, DOSAGE_FILTER and PHASE_FILTER, with deterministic reasons and
// CALL_SUMMARY. They are derived summaries of the fields the run already reports -- the representation filter,
// library sensitivity, the dosage stage's status and repair fields, zygosity and phase_support -- and replace none of
// them. There is no global PASS/FAIL value.
//
// Values, in fail-closed precedence (the headline is the most severe applying status; every reason is kept):
//     NO_CALL > NOT_EVALUATED > LOW_QUALITY > PASS > NOT_APPLICABLE
// A locus takes the most severe of its units (blocks for genotype and dosage, adjacent-block boundaries for phase),
// NOT_APPLICABLE only when every unit is, and lists each reason with its unit IDs (REASON@id1,id2). Per-status unit
// counts are reported beside it, so a fail-closed locus with one non-PASS block is not mistaken for a discarded one.

#include <array>
#include <cstddef>
#include <map>
#include <string>
#include <vector>

namespace panvar::headline {

constexpr const char* kHeadlineRule = "headline_filters_v2";
extern const std::array<const char*, 5> kStatuses;   // NO_CALL, NOT_EVALUATED, LOW_QUALITY, PASS, NOT_APPLICABLE
int severity(const std::string& status);

struct Verdict {
    std::string status = "PASS";
    std::vector<std::string> reasons;   // in the rule's order, without duplicates; empty = "none"
    std::string reason_text() const;
};
// The most severe status of the parts, with every part's reasons (duplicates removed, first occurrence kept).
Verdict combine(const std::vector<Verdict>& parts);

// What the run reported, locus-wide.
struct LocusFacts {
    std::string call_class = "no_call";          // strict_call, genotype_equivalent_call, no_call
    std::string no_call_reason = "NA";
    std::string phase_support = "NA";
    std::string library_sensitivity = "NOT_APPLICABLE";   // NOT_APPLICABLE, STABLE, SENSITIVE, NOT_EVALUATED
    std::string dosage_mode = "off";             // off, audit, repair
    std::string dosage_audit_status = "NA";      // the dosage stage's status (COMPLETE when it ran to completion)
};
// What the run reported for one called block.
struct BlockFacts {
    std::string block_id;
    std::string zygosity;                        // HET, HOM
    std::string representation_filter;           // PASS, LOWQUAL_REPRESENTATION_CONFLICT, NOT_EVALUATED...
    bool library_sensitive = false;              // listed in library_sensitive_blocks
    bool library_not_evaluated = false;          // listed in library_sensitivity_not_evaluated_blocks
    std::string dosage_status = "NA";            // NA when the stage did not report the block
    std::string repair_class = "NA";
    std::string dosage_content_scope = "NA";
    bool dosage_repair_applied = false;
};

Verdict genotype_verdict(const LocusFacts& locus, const BlockFacts& block);
Verdict dosage_verdict(const LocusFacts& locus, const BlockFacts& block);
Verdict phase_boundary_verdict(const LocusFacts& locus, const BlockFacts& left, const BlockFacts& right);

struct Aggregate {
    std::string status = "NOT_APPLICABLE";
    std::string reasons = "none";                // REASON@ids;REASON@ids
    std::map<std::string, std::size_t> counts;   // per status
    std::size_t units = 0, applicable = 0;
};
Aggregate aggregate(const std::vector<std::pair<std::string, Verdict>>& units);

struct LocusHeadline {
    Aggregate genotype, dosage, phase;
    std::string call_summary;                    // GENOTYPE_<s>;DOSAGE_<s>;PHASE_<s>
};

// Computes the headline for a finished run prefix from its own outputs (<prefix>.mosaic_status.tsv,
// .block_genotypes.tsv, .genotypes.tsv, .dosage_blocks.tsv), appends the four block columns (genotype_filter,
// genotype_filter_reasons, dosage_filter, dosage_filter_reasons) to both block tables, replacing them if present, and
// returns the locus values. `block_ids` names the model's blocks for a run with no block table (a no-call).
LocusHeadline apply_headline(const std::string& out_prefix, const std::vector<std::string>& block_ids,
                             const std::string& dosage_mode, const std::string& dosage_audit_status);

// The locus-table columns the headline adds, in order, and their values.
const std::vector<std::string>& locus_columns();
std::map<std::string, std::string> locus_values(const LocusHeadline& headline);

} // namespace panvar::headline
