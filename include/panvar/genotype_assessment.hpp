#pragma once

// assessment_v1 (M15-G): the primary v2 assessment of a published prediction, at locus and block level. The contract is
// tests/analysis/m15g/CONTRACT.md and its normative reference implementation tests/analysis/m15g/assessment_spec.py;
// this is a line-for-line port (the strings must agree byte for byte).
//
// It reads ONLY the run prefix's own outputs (.sequence_locus.tsv, .sequence_blocks.tsv, .mosaic_status.tsv,
// .block_genotypes.tsv, .dosage_status.tsv, .dosage_blocks.tsv) and uses a legacy block fact only in verified context:
// a completed v2 PREDICTION whose legacy decision exists and whose legacy reported alleles (sorted) equal the v2
// published genotype at that block. Everything else is NOT_ASSESSED with a reason; nothing is scored or recomputed.
// The old headline fields (headline_filters_v2) are left exactly as they are.

#include <string>
#include <utility>
#include <vector>

namespace panvar::assessment {

inline constexpr const char* kRule = "assessment_v1";

const std::vector<std::string>& locus_keys();     // the locus rows, in order
const std::vector<std::string>& block_columns();  // the block columns, in order

struct View {
    std::vector<std::pair<std::string, std::string>> locus;                       // locus_keys() order
    std::vector<std::pair<std::string, std::vector<std::string>>> blocks;         // block_id -> block_columns() order
};

// The assessment of a v2 run prefix, from its files.
View assess_prefix(const std::string& prefix);
// The exact text assessment_spec.py's render() prints (tests).
std::string render(const View& view);
// Writes the view: rows into .sequence_locus.tsv and columns into .sequence_blocks.tsv (replacing earlier assessment rows
// and columns, never duplicating them), mirrored into .locus_genotype.tsv (with headline_filters_status) and the legacy
// block tables when they exist. Throws on an unwritable file.
void write_view(const std::string& prefix, const View& view);
// assess_prefix + write_view, for a sequence_first_v2 prefix.
void write_assessment(const std::string& prefix);

}  // namespace panvar::assessment
