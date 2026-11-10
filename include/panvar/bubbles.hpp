#pragma once

#include <cstddef>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include "panvar/gfa.hpp"
#include "panvar/integrated_snarls.hpp"

namespace panvar {

struct Bubble {
    std::size_t id = 0;
    std::string source;
    std::string sink;
    // Orientation of each boundary as the reference traverses it. A boundary is a handle, not a node: a
    // node can bound a snarl on either side, and merging has to know which sides face each other.
    // Forward when the reference does not traverse the boundary.
    bool source_reverse = false;
    bool sink_reverse = false;
    // False when read from a site table with no orientation columns: the two above are then defaults.
    bool orientation_known = true;
    std::vector<std::string> inside;
    // Traversal support: paths that cross the bubble at all. Crossing a site is not evidence for any
    // particular allele, so on a full panel this is close to the panel size.
    std::size_t path_support = 0;
    // Allele support: distinct source->sink walks, paths taking the reference's walk, and the support of
    // the best- and least-supported non-reference walks.
    std::size_t distinct_alleles = 0;
    std::size_t ref_allele_support = 0;
    std::size_t alt_allele_support_max = 0;
    std::size_t alt_allele_support_min = 0;
    // Interior span in bp, summed over the interior nodes a path traverses: not the size of the difference
    // between alleles (a 1 bp substitution inside a 1 kb allele has a 1 kb span).
    std::size_t min_inside_bp = 0;
    std::size_t max_inside_bp = 0;
    // Supporting paths with inside span >= --min-variant-bp (drives the size filter).
    std::size_t long_path_support = 0;
    // True when an internal node is seen in both orientations across supporting paths.
    bool inversion_signal = false;
    // An inversion site: the stretch between the two flanks of an inversion read from the haplotypes'
    // walks (path_inversion.hpp); its breakpoints lie in inverted repeats the graph collapses, so no snarl
    // holds it. Its interior is every node a path walks between the flanks, and the ordinary sites inside
    // it name it in `parent_id`: the one kind of site that overlaps others.
    bool inversion_parent = false;
    std::size_t parent_id = 0;   // the inversion site enclosing this one (0: none)
    // Boundaries of the sites inside an inversion site; set by view_inversion_sites_on_own_nodes, not
    // stored.
    std::vector<std::pair<std::string, std::string>> inner_bounds;
};

struct BubbleCallOptions {
    // Optional JSON-lines from `vg view -R -j <snarls.pb>`, whose top-level pairs replace the internal
    // finder (integrated_snarls.hpp).
    std::string snarls_input_path;
    // Reference path (name or unique substring) used to order the internal snarl finder.
    std::string reference_path;
    // Precomputed top-level (source, sink) pairs; when non-empty they override both --snarls and the
    // internal finder.
    std::vector<std::pair<std::string, std::string>> snarl_pairs_override;
    // Inward sides of each override pair, parallel to `snarl_pairs_override` (--snarls reads them from
    // the `backward` flags). Empty: inferred from how the paths walk each pair.
    std::vector<SnarlInward> snarl_inward_override;
    // Whether a snarl source was supplied, whether or not it found anything: a graph with no snarl (a
    // linear path) then reports an empty table rather than an error.
    bool snarl_source_supplied = false;
    bool collect_snarl_debug = false;
    // Keep bubbles only if at least one supporting path has inside length >= this threshold
    // (unless inversion_signal is detected). Set to 0 to disable size filtering.
    std::size_t min_variant_bp = 50;
    // Drop bubbles whose longest supporting path exceeds this (0: no cap). Bounds the per-bubble walk
    // cost on hypervariable regions, at the price of SVs larger than the cap.
    std::size_t max_variant_bp = 0;
    std::size_t min_path_support = 0;
    // Minimum support for the best-supported non-reference allele: unlike min_path_support, this asks
    // whether an alternate has been seen more than once.
    std::size_t min_alt_support = 0;
    // Merge bubbles within this many reference bp between their interiors (facing boundaries included; a
    // shared boundary counts its own length). 0 disables merging.
    std::size_t merge_nearby_bp = 0;
    // Fuse sites whose variants could be placed in each other: every allele's change against the
    // reference allele, slid as far left and right as the reference allows, covers a stretch, and sites
    // whose stretches overlap or touch are one locus the graph split. A group holding a retained site is
    // fused, sites under the size floor included.
    bool fuse_ambiguous = false;
    // Suppress the stderr progress bar (output files are unaffected).
    bool quiet = false;
    // Module named in log lines: the module running the snarl pass (panphorte and refine re-snarl).
    std::string log_module = "bubble";
};

struct SnarlDebugEntry {
    std::size_t candidate_id = 0;
    std::string source;
    std::string sink;
    std::size_t inside_node_count = 0;
    std::size_t n_paths = 0;
    std::size_t min_inside_bp = 0;
    std::size_t long_path_support = 0;
    bool inversion_signal = false;
    bool accepted = false;
};

// What each step did, for the run's log.
struct BubbleRunCounts {
    std::size_t snarls = 0;            // boundary pairs found or read
    std::size_t candidates = 0;        // sites scored
    std::size_t widened = 0;           // sites inside a loop widened to the sequence around it
    std::size_t passed = 0;            // sites passing the filters
    std::size_t merged = 0;            // sites made by merging nearby ones
    std::size_t fused = 0;             // sites made by fusing over placements
    std::size_t overlaps_dropped = 0;  // sites dropped for overlapping a smaller one
    std::size_t rescued = 0;           // sites joined from pieces under the size floor
    std::size_t inversion_sites = 0;
};

struct BubbleCallReport {
    std::vector<Bubble> bubbles;
    BubbleRunCounts counts;
    // Candidate non-SNP bubbles before final filters (used for overview coloring).
    std::vector<Bubble> non_snp_bubbles;
    std::vector<SnarlDebugEntry> snarl_debug;
};

BubbleCallReport call_bubbles_report(const Graph& graph, const BubbleCallOptions& options);

// An inversion site viewed on its own nodes: its interior less the interiors of the sites inside it, whose
// boundaries it keeps -- the stretch the inversion turns, walked forward or reversed. A step from one inner
// site's boundary to its other is that site's (its deletion), not the inversion site's: see
// inner_site_step.
void view_inversion_sites_on_own_nodes(std::vector<Bubble>& bubbles);
bool inner_site_step(const Bubble& bubble, const std::string& from, const std::string& to);

} // namespace panvar
