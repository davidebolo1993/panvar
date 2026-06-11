#pragma once

#include <cstddef>
#include <string>

namespace panvar {

// Parsed genomic coordinates from a reference path name such as
// "GRCh38#0#chr6:31891045-32123783" -> chrom="chr6", region_start_1based=31891045.
// When the name carries no parseable interval, has_interval is false and chrom
// falls back to a sensible token (after the last '#', or from "chr...").
struct ParsedReferencePath {
    std::string chrom;
    std::size_t region_start_1based = 1;
    bool has_interval = false;
};

ParsedReferencePath parse_reference_path_label(const std::string& path_name);

} // namespace panvar
