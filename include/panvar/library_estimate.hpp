#pragma once

// Automatic library estimation, DIAGNOSTIC ONLY (docs/algorithms/library_estimate.md). Nothing here is used by a call:
// the estimate is written beside a run, or by genotype_probe library-estimate, and compared with the declared library.
//
// Read-level quantities come from the reads directly. The insert (fragment) distribution comes from pairs placed
// decisively against the model's panel templates: both mates hit anchor k-mers (present exactly once, on one strand,
// in every template), the pair is FR-concordant, and the implied insert is identical on every template, so it depends
// on no candidate genotype. Truth-free: it reads the model graph and the reads only.

#include "panvar/genotype_model.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace panvar {

struct LibraryEstimateOptions {
    std::size_t k = 31;
    std::size_t anchor_sampling = 4;          // keep anchor k-mers whose hash is 0 modulo this (index size)
    std::size_t max_pairs = 20000;            // deterministic sample: the pairs with the smallest salted name hashes
    std::size_t min_accepted_pairs = 500;     // fewer accepted placements: the insert falls back (measured, library_estimate.md)
    std::size_t error_pairs = 2000;           // accepted pairs used for the invariant-site mismatch rate
    double outlier_sds = 6.0;                 // |insert - median| beyond this many robust sd: outlier
    double max_outlier_fraction = 0.05;       // above: contamination warning (--library auto refuses)
    double max_discordant_fraction = 0.05;    // orientation-discordant among anchored pairs above: contamination (auto refuses)
    double multimodal_separation_sds = 4.0;   // two-group split separated by this many within-group sd ...
    double multimodal_min_fraction = 0.05;    // ... with each group at least this fraction: multimodal
    std::size_t marker_kmer = 31, marker_syncmer = 15;   // the caller's marker settings, for the depth-fit route
    std::size_t threads = 1;
    std::string salt = "panvar-library-estimate-v1";
};

struct LibraryEstimateRow {
    std::string parameter, value, method;
    std::size_t pairs_examined = 0, pairs_accepted = 0;
    std::string spread = "NA", status = "ESTIMATED", reason = "none", source = "reads";
};

struct LibraryEstimate {
    std::vector<LibraryEstimateRow> rows;
    std::vector<long> accepted_inserts;       // in sample order
    std::size_t pairs_total = 0, pairs_sampled = 0, both_anchored = 0, concordant = 0, orientation_discordant = 0,
                mate_inconsistent = 0, not_decisive = 0, accepted = 0, outliers = 0, anchors = 0, templates = 0;
    double insert_median = 0.0, insert_robust_sd = 0.0, insert_mean = 0.0, insert_sd = 0.0;
    double read1_length_median = 0.0, read2_length_median = 0.0;
    double invariant_mismatch_rate = 0.0, phred_error_mean = 0.0, anchor_depth_per_copy = 0.0, fragment_start_intensity = 0.0,
           marker_fit_intensity = 0.0;
    bool multimodal = false;
    double seconds_index = 0.0, seconds_placement = 0.0, seconds_total = 0.0;
};

// declared: an explicit library to compare with (and the fallback source when the insert cannot be estimated); may be
// null, in which case a fallback is reported as missing.
LibraryEstimate estimate_library(const PreparedGenotypeModel& model, const std::vector<ReadPair>& reads,
                                 const LibraryEstimateOptions& options, const LibraryModel* declared);

// <out>.library_estimate.tsv: parameter, value, method, pairs_examined, pairs_accepted, spread, status, reason, source,
// declared_value.
std::string library_estimate_table(const LibraryEstimate& estimate, const LibraryModel* declared);

// --library auto (docs/algorithms/library_estimate.md section 5). The insert centre and spread come from the estimate
// (mean and sd of the accepted inserts within the outlier bound); every other library value comes from `preset`. The
// estimate is refused -- never used, never replaced by a guess -- when accepted pairs are fewer than min_accepted_pairs,
// the insert distribution is multimodal, or the outlier or orientation-discordant fraction exceeds its bound; the
// refusal names --library FILE.
struct AutoLibrary {
    bool accepted = false;
    LibraryModel library;
    std::string refusal;                                              // empty when accepted
    std::vector<std::pair<std::string, std::string>> sources;        // parameter -> reads | placement_estimate | short_read_preset | explicit_override
};
AutoLibrary build_auto_library(const LibraryEstimate& estimate, const LibraryModel& preset, const LibraryEstimateOptions& options);

// The rows appended to <out>.library_estimate.tsv for the library the call used: one per parameter, with its value and
// source (parameter "library.<name>", status USED).
std::string library_used_rows(const LibraryModel& used, const std::vector<std::pair<std::string, std::string>>& sources);

} // namespace panvar
