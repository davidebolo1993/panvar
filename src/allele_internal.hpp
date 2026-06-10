#pragma once

// Internal interfaces shared between the allele clustering core (allele.cpp) and
// the per-bubble similarity diagnostics (allele_similarity.cpp). Not part of the
// public panvar API.

#include <cstddef>
#include <string>
#include <vector>

#include "call_internal.hpp"
#include "panvar/bubbles.hpp"
#include "panvar/gfa.hpp"

namespace panvar {

// Absolute and length-normalized pairwise distance matrices for the unique
// alleles of one bubble.
struct DistanceMatrices {
    std::vector<std::vector<int>> abs;
    std::vector<std::vector<double>> norm;
};

// Output paths for the per-bubble similarity diagnostics bundle.
struct SimilarityReportPaths {
    std::string bubble_dir;
    std::string alleles_tsv;
    std::string matrix_norm_tsv;
    std::string matrix_abs_tsv;
    std::string stats_tsv;
    std::string cluster_signatures_dir;
};

inline double safe_div(double num, double den) {
    return den == 0.0 ? 0.0 : (num / den);
}

SimilarityReportPaths similarity_report_paths_for_bubble(
    const std::string& output_dir,
    std::size_t bubble_id);

void write_similarity_reports_for_bubble(
    const SimilarityReportPaths& paths,
    const Bubble& bubble,
    const Graph& graph,
    const std::vector<UniqueAllele>& unique_alleles,
    const std::vector<Cluster>& clusters,
    const std::vector<std::size_t>& cluster_of_unique,
    const DistanceMatrices& dists,
    double min_similarity);

} // namespace panvar
