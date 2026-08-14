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
    // A copy whose boundary falls inside a node cannot be folded without splitting that node. Declining
    // only that copy leaves the site half REP and half literal, and `call` counts REP occurrences -- so
    // a haplotype carrying one copy is reported CN 0. The site is therefore refused as a whole.
    // Setting this true restores per-copy refusal: more sites fold, at the cost of that false CN on the
    // affected haplotypes. Off by default, and it warns with the counts.
    bool allow_partial_boundary = false;
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
    // Reference-coordinate GTF projected onto the sorted normalized graph. Handled here rather than by
    // the caller so the gene annotation is part of the same staged transaction as everything else: run
    // afterwards it could fail with the normalized family already on disk.
    std::string gtf_path;
    bool quiet = false;
};

struct PanphorteSummary {
    std::size_t bubbles_seen = 0;
    std::size_t bubbles_normalized = 0;
    std::size_t paths_rewritten = 0;
    std::size_t nodes_removed = 0;
    std::size_t nodes_added = 0;
    // Of nodes_added, the per-occurrence fragments emitted around an accepted copy (the bases inside
    // the replaced step range that fall outside the copy). One array with thousands of copies can add
    // thousands of these, so the growth they account for is reported on its own.
    std::size_t fragment_nodes_added = 0;
    std::size_t edges_added = 0;
    // Oriented links that only the pre-rewrite paths traversed, dropped so the replaced branch does not
    // survive beside the normalized one.
    std::size_t edges_removed = 0;
    // Set when --reference-path triggers the internal sort + re-snarl.
    bool sorted = false;
    std::size_t resnarled_bubbles = 0;
    bool genes_written = false;
};

void panphorte_normalize(const PanphorteOptions& options, PanphorteSummary* summary_out = nullptr);

} // namespace panvar
