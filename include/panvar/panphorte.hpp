#pragma once

#include <cstddef>
#include <map>
#include <string>
#include <vector>

namespace panvar {

struct PanphorteOptions {
    std::string gfa_path;
    std::string bubbles_csv_in;     // module-1 bubbles CSV
    std::string out_prefix;         // writes <prefix>.normalized.gfa + <prefix>.panphorte.report.tsv
    std::size_t min_unit_bp = 50;   // minimum repeat-unit span to normalize
    // Interior-span filter of the re-snarl under --reference-path. Set here, not taken from
    // BubbleCallOptions, whose default may differ from the one that produced the input CSV.
    std::size_t resnarl_min_variant_bp = 50;
    // The re-snarl fuses sites whose variants could be placed in each other, as bubble does, so the output
    // CSV keeps the sites bubble fused. Proximity fusion (bubble's --merge-nearby-bp) is off by default:
    // after folding it can join a folded module to a neighbouring repeat site, and a copy-number record's
    // CNBP, measured across its whole bubble, would then take in the other site's indels.
    std::size_t resnarl_merge_nearby_bp = 0;
    bool resnarl_fuse_ambiguous = true;
    std::size_t min_copies = 2;     // minimum tandem copies to normalize
    // Minimum fraction of bubble-traversing haplotypes carrying a >=min_copies array for the bubble to
    // be normalized. Separates a population VNTR (folded) from a rare private duplication of a gene or
    // segmental module (left for `call`, since folding would collapse paralogs).
    double min_array_prevalence = 0.5;
    double max_interruption_frac = 0.25; // tolerance for interrupting bases within an array
    // Minimum identity for a block to count as a copy of the repeat unit: 1.0 is exact; < 1.0 enables
    // approximate, alignment-based detection.
    double min_similarity = 1.0;
    // With min_similarity < 1, copies are detected approximately but folded onto a unit graph (shared
    // segments, each copy's differences as nodes of their own), so every path keeps its sequence. False
    // (--lossy-units) builds it from the reference's copies only: the reference keeps its sequence, every
    // other haplotype walks the plain unit once per copy.
    bool lossless_units = true;
    std::size_t threads = 0;        // 0 = hardware concurrency (approximate detection)
    std::vector<std::size_t> bubble_ids; // if non-empty, restrict to these bubbles
    // When set, after normalization the graph is internally sorted+flipped along this
    // reference and re-snarled (cactus), writing <prefix>.normalized.sorted.gfa and
    // <prefix>.bubbles.csv so `call` can run with no external tools.
    std::string reference_path;
    bool no_flip = false;
    // Reference-coordinate GTF projected onto the sorted normalized graph. Handled here so the annotation
    // is part of the same staged transaction: run afterwards, it could fail with the outputs already on disk.
    std::string gtf_path;
    bool quiet = false;
    // --refine: after folding, realign the flanks around the folded array at every folded site (one
    // alignment before the array, one after, over the haplotypes that walk the unit node), keeping the
    // result only where it is simpler. The array itself and haplotypes without a unit node are untouched.
    // Off: output unchanged. Needs --reference-path. Limits as refine's --max-poa-bp / --max-walks /
    // --max-poa-work, applied per flank alignment.
    bool refine = false;
    std::size_t refine_max_poa_bp = 5000;
    std::size_t refine_max_walks = 500;
    std::size_t refine_max_poa_work = 0;
};

struct PanphorteSummary {
    std::size_t bubbles_seen = 0;
    std::size_t bubbles_normalized = 0;
    std::size_t paths_rewritten = 0;
    std::size_t nodes_removed = 0;
    std::size_t nodes_added = 0;
    // Of nodes_added, the fragment nodes: bases inside a replaced step range that fall outside the copies.
    // A large array can add many of these, so they are reported on their own.
    std::size_t fragment_nodes_added = 0;
    std::size_t edges_added = 0;
    // Oriented links that only the pre-rewrite paths traversed, dropped so the replaced branch does not
    // survive beside the normalized one.
    std::size_t edges_removed = 0;
    // Rewritten paths whose original route is still walkable end to end, so their site would be
    // represented twice. Those sites are declined (status declined_shared_route) and the run repeated, so in
    // a completed run this is 0.
    std::size_t routes_surviving = 0;
    // Set when --reference-path triggers the internal sort + re-snarl.
    bool sorted = false;
    std::size_t resnarled_bubbles = 0;
    bool genes_written = false;
    // --refine: folded sites rebuilt, kept as folded because the rebuild was not simpler, and skipped.
    std::size_t refine_refined = 0;
    std::size_t refine_kept = 0;
    std::size_t refine_skipped = 0;
    std::size_t refine_nodes_added = 0;
    std::size_t refine_nodes_removed = 0;
    // Sites by report status (normalized, below_prevalence, ...).
    std::map<std::string, std::size_t> status_counts;
};

void panphorte_normalize(const PanphorteOptions& options, PanphorteSummary* summary_out = nullptr);

} // namespace panvar
