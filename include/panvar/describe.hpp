#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace panvar {

struct DescribeOptions {
    std::string vcf_in_path;
    std::string out_dir = "describe_out";
    std::string gtf_path;
    std::vector<std::string> gene_match_patterns;
    std::vector<std::size_t> size_bins = {100, 1000};
    bool quiet = false;
};

struct DescribeSummary {
    std::size_t bubbles = 0;
    std::size_t events = 0;
    std::size_t haplotype_rows = 0;
    std::size_t files_written = 0;
};

void describe_from_region_vcf(
    const DescribeOptions& options,
    DescribeSummary* summary_out = nullptr);

} // namespace panvar
