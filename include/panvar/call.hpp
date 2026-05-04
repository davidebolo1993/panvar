#pragma once

#include <string>

#include "panvar/allele.hpp"

namespace panvar {

void call_variants_from_precomputed_alleles(
    const Graph& graph,
    const AlleleCallOptions& options,
    const std::string& clusters_csv_in_path,
    const std::string& assignments_csv_in_path,
    AlleleCallSummary* summary_out = nullptr);

} // namespace panvar
