#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace panvar {

// Parsed genomic coordinates from a reference path name such as
// "GRCh38#0#chr6:31891045-32123783" -> chrom="chr6", region_start_1based=31891046.
// The interval's convention is read from the path's length (path_length_bp, 0 when unknown):
// end - start bases is half-open with a 0-based start, as pggb, odgi and impg name subgraph
// paths, so the first base is start + 1; end - start + 1 bases is 1-based inclusive, as panvar
// names W-line paths. A length matching neither, or no length, is read as half-open and sets
// length_mismatch when a length was given -- except for a W-line path (one_based), whose name panvar
// built 1-based from the W line's 0-based start, so no guess is made: reading a mismatching one as
// 0-based would shift it by one base twice; such a name may end "-*" when the W line gives a start and no end.
// When the name carries no parseable interval,
// has_interval is false and chrom falls back to a sensible token (after the last '#', or from
// "chr...").
struct ParsedReferencePath {
    std::string chrom;
    std::size_t region_start_1based = 1;
    bool has_interval = false;
    bool half_open = false;        // the interval was read as 0-based half-open
    bool length_mismatch = false;  // the interval's length is neither end - start nor end - start + 1
};

// Which of `names` a -r query names: the exact name; else the one name containing the query as written;
// else, only when no name contains it as written, the one containing it in any case. So "GRCh38" picks
// "GRCh38#0#chr6" beside "grch38#1#chr6". The names it matched are in `matches`: none or several mean no
// answer.
std::vector<std::size_t> match_reference_name(const std::vector<std::string>& names, const std::string& query);

ParsedReferencePath parse_reference_path_label(const std::string& path_name, std::size_t path_length_bp = 0,
                                               bool one_based = false);

} // namespace panvar
