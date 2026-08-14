#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace panvar {

struct PanphorteOptions {
    std::string gfa_path;
    std::string bubbles_csv_in;     // module-1 bubbles CSV
    std::string out_prefix;         // writes <prefix>.normalized.gfa + <prefix>.panphorte.report.tsv
    std::size_t min_unit_bp = 50;   // minimum repeat-unit span to normalize
    // Interior-span filter applied when the normalized graph is re-snarled under --reference-path.
    // This used to be BubbleCallOptions' own default, applied silently to a CSV that may have been
    // produced with a different one.
    std::size_t resnarl_min_variant_bp = 50;
    std::size_t min_copies = 2;     // minimum tandem copies to normalize
    // Minimum fraction of bubble-traversing haplotypes carrying a >=min_copies array for the bubble to
    // be normalized. Separates a population VNTR (folded) from a rare private duplication of a gene or
    // segmental module (left for `call`, since folding would collapse paralogs). Default 0.5.
    double min_array_prevalence = 0.5;
    double max_interruption_frac = 0.25; // tolerance for interrupting bases within an array
    // Minimum identity to treat a block as a copy of the repeat unit. 1.0 = exact
    // (current behavior); < 1.0 enables approximate, similarity-based detection.
    double min_similarity = 1.0;
    std::size_t threads = 0;        // 0 = hardware concurrency (approximate detection)
    std::vector<std::size_t> bubble_ids; // if non-empty, restrict to these bubbles
    // When set, after normalization the graph is internally sorted+flipped along this
    // reference and re-snarled (cactus), writing <prefix>.normalized.sorted.gfa and
    // <prefix>.bubbles.csv so `call` can run with no external tools.
    std::string reference_path;
    bool no_flip = false;
    bool quiet = false;
};

struct PanphorteSummary {
    std::size_t bubbles_seen = 0;
    std::size_t bubbles_normalized = 0;
    std::size_t paths_rewritten = 0;
    std::size_t nodes_removed = 0;
    std::size_t nodes_added = 0;
    std::size_t edges_added = 0;
    // Set when --reference-path triggers the internal sort + re-snarl.
    bool sorted = false;
    std::size_t resnarled_bubbles = 0;
};

void panphorte_normalize(const PanphorteOptions& options, PanphorteSummary* summary_out = nullptr);

} // namespace panvar
