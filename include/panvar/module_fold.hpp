#pragma once

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace panvar {

// Copy-number modules on a rebuilt graph. The generator keeps a duplication's copies apart, so a copy gain is a
// long insertion and nothing counts copies. A module is a stretch the haplotypes repeat at --min-identity or
// more (gap-compressed), bounded by anchors unique on the reference. Folding it makes every copy walk one unit
// subgraph -- the unit cut at every difference a copy has from it, each distinct difference its own node -- so
// a copy-number change is a change in passes through the unit's nodes. Each path's module stretch is taken from
// the haplotype's own sequence, since the rebuilt walk spells differences under --min-var as the backbone's.
// A path is folded only if it then spells exactly what it should; one that would not keeps its walk.
struct ModuleFoldOptions {
    std::size_t min_copy_bp = 3000;   // shortest repeated stretch that can make a module
    double min_identity = 0.95;       // identity of the repeat that makes a module; the most a copy is held to
    double min_copy_identity = 0.90;  // the least a copy is held to (see FoldedModule::copy_identity)
    double anchor_identity = 0.95;    // a hit this close makes a flank or an anchor not unique on the reference
    double min_unit_cover = 0.5;      // share of the unit a copy covers
    std::size_t anchor_bp = 2000;     // length of the reference anchors bounding a module
    std::size_t max_anchor_shift = 50000;   // how far an anchor may move to become unique
    std::size_t template_min_bp = 200;      // sequence around copies, or a difference, aligned to templates
    std::size_t hub_degree = 50;            // distinct neighbours that make a module segment a hub (rebuild passes --hub-degree)
    std::size_t gate_hub_degree = 50;       // the same for the folded graph's summary (rebuild's gate)
    std::size_t threads = 1;
};

struct FoldedModule {
    std::size_t ref_start = 0, ref_end = 0;   // the folded stretch on the reference path, 0-based, end exclusive
    std::size_t period = 0;                   // one copy's step along the array
    std::size_t unit_bp = 0;                  // the unit every copy walks
    std::string unit_source;                  // the path the unit was taken from
    double copy_identity = 0;                 // a copy's bar: the reference copies' lowest identity less 0.02,
                                              // within [min_copy_identity, min_identity]
    std::size_t reference_copies = 0;
    std::size_t min_copies = 0, max_copies = 0;   // over the folded paths
    std::size_t paths = 0;                    // all paths
    std::size_t paths_spanning = 0;           // paths walking both cut nodes (the others keep their walks)
    std::size_t paths_folded = 0;             // spanning paths whose folded walk spells as it should
    std::size_t paths_exact = 0;              // folded paths whose module stretch is their own sequence
    std::vector<std::string> folded_paths;    // the paths folded through this module, by name
    std::size_t stretch_templates = 0;        // templates for the sequence around the copies
    std::size_t difference_templates = 0;     // templates for long differences inside the copies
    std::size_t hubs = 0;                     // new segments with hub_degree or more distinct neighbours
    bool folded = false;
    std::string reason;                       // why not, when not folded
};

struct ModuleFoldResult {
    std::string gfa;                          // the folded graph (the input unchanged when nothing folded)
    std::vector<FoldedModule> modules;
    // Per emitted path, by name: whether its path spells the given haplotype sequence exactly.
    std::vector<std::pair<std::string, bool>> spells_own;
    std::size_t nodes = 0, links = 0, maxdeg = 0, hubs = 0, selfloops = 0;   // of the folded graph; selfloops: nodes linked to themselves
};

// A copy-number module found before folding: its unit, the path the unit was taken from, one copy's step along
// the array, and its stretch [lo, hi) on the reference.
struct ModuleUnit {
    std::string unit, source;
    std::size_t period = 0, lo = 0, hi = 0;
};

// The modules in a locus, from the sequences given (`ref` indexes the reference among them).
std::vector<ModuleUnit> discover_modules(const std::vector<std::string>& names, const std::vector<std::string>& seqs,
                                         std::size_t ref, const ModuleFoldOptions& options);

// `gfa` is the rebuilt graph (S and P lines; links are re-derived from the walks). `haplotypes` gives each
// emitted path's own sequence, by P-line name; `reference` names the reference path. `known` gives the modules
// when they were already found from these sequences; otherwise they are found here.
ModuleFoldResult fold_modules(const std::string& gfa,
                              const std::vector<std::pair<std::string, std::string>>& haplotypes,
                              const std::string& reference,
                              const ModuleFoldOptions& options,
                              const std::vector<ModuleUnit>* known = nullptr);

} // namespace panvar
