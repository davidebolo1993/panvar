#include "panvar/genotype_locus.hpp"

#include "panvar/bubble_path.hpp"
#include "panvar/graph_utils.hpp"

#include <algorithm>
#include <cstdlib>
#include <map>
#include <optional>
#include <set>
#include <stdexcept>
#include <unordered_set>

namespace panvar {
namespace {

enum class ChainBlockKind { Flank, Bubble, Backbone };

struct ChainBlock {
    ChainBlockKind kind = ChainBlockKind::Bubble;
    std::size_t bubble_id = 0;
    const Bubble* bubble = nullptr;
    std::string source;      // empty on a leading flank
    std::string sink;        // empty on a trailing flank
    // Two bubbles can abut, sharing one boundary node with no backbone between them. Both would then
    // spell that node and the chain would count the same bases twice. The later bubble gives it up.
    // Named rather than front/back, because a reverse-oriented bubble carries its left endpoint at
    // the END of its canonical source->sink steps.
    std::string trim_node;
};

std::uint64_t mix64(std::uint64_t value) {
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31);
}

class Digest {
public:
    void add(const std::string& value) {
        for (const unsigned char byte : value) {
            state_ = mix64(state_ ^ byte);
        }
        add(value.size());
    }
    void add(std::size_t value) { state_ = mix64(state_ ^ static_cast<std::uint64_t>(value)); }

    std::string hex() const {
        static const char* digits = "0123456789abcdef";
        std::string out(16, '0');
        for (std::size_t i = 0; i < 16; ++i) {
            out[15 - i] = digits[(state_ >> (4 * i)) & 0xfULL];
        }
        return out;
    }

private:
    std::uint64_t state_ = 0xcbf29ce484222325ULL;
};

// The graph is reference-sorted, so numeric node order is reference order. A bubble whose paths run
// backwards is stored with source and sink swapped; taking them at face value makes the neighbouring
// backbone intervals span straight across the bubble's interior, putting the same sequence in
// several blocks at once and destroying the markers of every block that touches it.
unsigned long long node_position(const std::string& node) {
    char* end = nullptr;
    const unsigned long long value = std::strtoull(node.c_str(), &end, 10);
    return (end != nullptr && *end == '\0') ? value : 0ULL;
}

const std::string& left_boundary(const Bubble& bubble) {
    return node_position(bubble.source) <= node_position(bubble.sink) ? bubble.source : bubble.sink;
}

const std::string& right_boundary(const Bubble& bubble) {
    return node_position(bubble.source) <= node_position(bubble.sink) ? bubble.sink : bubble.source;
}

// A chain is a partition. `bubble` may legitimately report a container and its children together,
// but keeping both would place the contained sequence in two blocks. Keep the leaves: the finer
// decomposition is what a panel-based genotyper needs, and a container spanning tens of kb is
// near-singleton across the panel and cannot represent a held-out sample.
std::vector<const Bubble*> drop_contained(
    const std::vector<Bubble>& bubbles,
    std::size_t& dropped) {
    std::vector<std::unordered_set<std::string>> inside(bubbles.size());
    for (std::size_t i = 0; i < bubbles.size(); ++i) {
        inside[i].insert(bubbles[i].inside.begin(), bubbles[i].inside.end());
    }
    std::vector<const Bubble*> kept;
    for (std::size_t i = 0; i < bubbles.size(); ++i) {
        bool contains_another = false;
        for (std::size_t j = 0; j < bubbles.size() && !contains_another; ++j) {
            if (i == j || inside[j].empty() || inside[j].size() >= inside[i].size()) continue;
            contains_another = true;
            for (const std::string& node : inside[j]) {
                if (inside[i].count(node) == 0) {
                    contains_another = false;
                    break;
                }
            }
        }
        if (contains_another) {
            ++dropped;
        } else {
            kept.push_back(&bubbles[i]);
        }
    }
    if (kept.empty()) {
        dropped = 0;
        for (const Bubble& bubble : bubbles) kept.push_back(&bubble);
    }
    return kept;
}

std::vector<ChainBlock> build_chain(
    const std::vector<const Bubble*>& bubbles,
    std::size_t& trimmed) {
    std::vector<ChainBlock> chain;

    // Leading flank. Skipping it would leave several percent of every template outside the chain,
    // and the round-trip gate would fail on the first template.
    ChainBlock leading;
    leading.kind = ChainBlockKind::Flank;
    leading.sink = left_boundary(*bubbles.front());
    chain.push_back(leading);

    for (std::size_t i = 0; i < bubbles.size(); ++i) {
        std::string abutting;
        if (i > 0) {
            ChainBlock backbone;
            backbone.kind = ChainBlockKind::Backbone;
            backbone.source = right_boundary(*bubbles[i - 1]);
            backbone.sink = left_boundary(*bubbles[i]);
            if (backbone.source != backbone.sink) {
                chain.push_back(backbone);
            } else {
                abutting = backbone.source;
                ++trimmed;
            }
        }
        ChainBlock block;
        block.kind = ChainBlockKind::Bubble;
        block.bubble_id = bubbles[i]->id;
        block.bubble = bubbles[i];
        block.source = bubbles[i]->source;
        block.sink = bubbles[i]->sink;
        block.trim_node = abutting;
        chain.push_back(block);
    }

    ChainBlock trailing;
    trailing.kind = ChainBlockKind::Flank;
    trailing.source = right_boundary(*bubbles.back());
    chain.push_back(trailing);
    return chain;
}

// One template's walk across one chain block, already canonicalised to chain orientation by the
// path primitives. `traverses` separates a template that walks the block and finds nothing between
// the boundaries (an empty allele) from one that does not reach the block at all (a bypass, which is
// usually a deletion spanning the boundaries). Both spell no bases, so they collapse to one
// canonical allele, but they are distinct source records and stay separately auditable.
struct BlockWalk {
    bool traverses = false;
    std::string signature;
    std::string sequence;
};

BlockWalk walk_block(
    const Graph& graph,
    const PathRecord& path,
    const BubblePathIndex& index,
    const ChainBlock& block) {
    BlockWalk walk;
    std::optional<std::vector<PathStep>> steps;
    if (block.kind == ChainBlockKind::Bubble) {
        steps = bubble_steps(path, index, *block.bubble);
        if (steps.has_value() && !block.trim_node.empty() && steps->size() > 1) {
            if (steps->front().node_id == block.trim_node) {
                steps->erase(steps->begin());
            } else if (steps->back().node_id == block.trim_node) {
                steps->pop_back();
            }
        }
    } else if (block.kind == ChainBlockKind::Flank) {
        const bool leading = block.source.empty();
        const std::string& boundary = leading ? block.sink : block.source;
        steps = flank_steps(path, index, boundary, leading);
        if (!steps.has_value() && index.positions.find(boundary) != index.positions.end()) {
            // The boundary is on the path, but nothing lies beyond it: the flank is EMPTY, not
            // absent. A bubble owns its boundary nodes, so a template whose path starts exactly at
            // the first bubble has a zero-length leading flank and has still observed the block.
            // Only a boundary the path never reaches is missing coverage.
            steps.emplace();
        }
    } else {
        // The boundary nodes belong to the bubbles either side, so the backbone takes only what lies
        // strictly between them. `interval_steps` decides traversal; `interval_interior_steps`
        // returning nothing means an empty interior, not an absent template.
        if (interval_steps(path, index, block.source, block.sink).has_value()) {
            walk.traverses = true;
            steps = interval_interior_steps(path, index, block.source, block.sink);
            if (!steps.has_value()) steps.emplace();
        }
    }
    if (!steps.has_value()) return walk;
    walk.traverses = true;
    walk.signature = build_walk_signature(*steps);
    walk.sequence = spell_path_steps_sequence(graph, *steps);
    return walk;
}

// Orientation is derived structurally, from the order in which a template visits the chain's ordered
// boundary nodes — never by comparing spelled bytes against the authoritative sequence. Deriving it
// from the bytes would make the round-trip gate tautological: it would pass by construction on any
// chain, including one whose blocks do not tile.
std::optional<ChainOrientation> template_orientation(
    const PathRecord& path,
    const std::vector<std::string>& ordered_boundaries) {
    std::map<std::string, std::size_t> first_occurrence;
    for (std::size_t at = 0; at < path.steps.size(); ++at) {
        first_occurrence.emplace(path.steps[at].node_id, at);
    }
    std::vector<std::size_t> positions;
    for (const std::string& node : ordered_boundaries) {
        const auto found = first_occurrence.find(node);
        if (found != first_occurrence.end()) positions.push_back(found->second);
    }
    if (positions.size() < 2) return std::nullopt;

    const bool ascending = std::is_sorted(positions.begin(), positions.end());
    const bool descending = std::is_sorted(positions.rbegin(), positions.rend());
    if (ascending == descending) return std::nullopt;   // neither, or a single repeated position
    return ascending ? ChainOrientation::Forward : ChainOrientation::Reverse;
}

// A bubble's source and sink are traversed by every path through it, so their bases sit in the
// variant block for every allele. Where the graph cuts those boundary nodes is representation, not
// model: splitting one moves bases between a variant block and its context neighbour while leaving
// the concatenated locus identical. Trimming each variant block down to the span over which its
// alleles actually differ removes that freedom — the maximal trim is the same whichever way the
// boundary node was cut — and it is also the honest definition of an SV allele.
//
// The shared flanks move outward into the context blocks, which are never trimmed: their boundaries
// are then defined by where the variant alleles diverge, which is representation-independent too.
// A block carrying a bypass allele has no common flank to trim, so it is left alone.
struct BlockSequences {
    std::vector<std::string> source_ids;
    std::vector<std::string> sequences;   // parallel to source_ids
    bool has_bypass = false;
};

std::size_t common_prefix(const std::vector<std::string>& values) {
    if (values.empty()) return 0;
    std::size_t length = values.front().size();
    for (std::size_t i = 1; i < values.size(); ++i) {
        const std::string& other = values[i];
        length = std::min(length, other.size());
        std::size_t at = 0;
        while (at < length && values.front()[at] == other[at]) ++at;
        length = at;
    }
    return length;
}

std::size_t common_suffix(const std::vector<std::string>& values) {
    if (values.empty()) return 0;
    std::size_t length = values.front().size();
    for (std::size_t i = 1; i < values.size(); ++i) {
        const std::string& other = values[i];
        length = std::min(length, other.size());
        std::size_t at = 0;
        while (at < length &&
               values.front()[values.front().size() - 1 - at] == other[other.size() - 1 - at]) {
            ++at;
        }
        length = at;
    }
    return length;
}

struct TrimOutcome {
    std::size_t prefix = 0;
    std::size_t suffix = 0;
    bool representation_dependent = false;
};

TrimOutcome trim_variant_block_to_divergence(
    BlockSequences& before,
    BlockSequences& block,
    BlockSequences& after,
    std::size_t& moved_bases) {
    TrimOutcome outcome;
    if (block.sequences.size() < 2) return outcome;
    if (block.has_bypass) {
        // The bypassing template has none of the boundary bytes the traversing alleles share, so they
        // cannot be moved out uniformly. If any exist, the block boundary is where the graph cut a
        // node, and a node split would change the model.
        std::vector<std::string> traversing;
        for (std::size_t i = 0; i < block.sequences.size(); ++i) {
            if (block.source_ids[i] != "bypass") traversing.push_back(block.sequences[i]);
        }
        if (!traversing.empty()) {
            std::size_t shortest = traversing.front().size();
            for (const std::string& value : traversing) shortest = std::min(shortest, value.size());
            const std::size_t prefix = std::min(common_prefix(traversing), shortest);
            const std::size_t suffix = std::min(common_suffix(traversing), shortest - prefix);
            outcome.representation_dependent = prefix + suffix > 0;
        }
        return outcome;
    }
    std::size_t shortest = block.sequences.front().size();
    for (const std::string& value : block.sequences) shortest = std::min(shortest, value.size());
    std::size_t prefix = std::min(common_prefix(block.sequences), shortest);
    std::size_t suffix = std::min(common_suffix(block.sequences), shortest - prefix);
    if (prefix == 0 && suffix == 0) return outcome;

    const std::string prefix_bases = block.sequences.front().substr(0, prefix);
    const std::string suffix_bases =
        block.sequences.front().substr(block.sequences.front().size() - suffix);
    for (std::string& value : block.sequences) {
        value = value.substr(prefix, value.size() - prefix - suffix);
    }
    // The trimmed flanks are identical for every allele, so every template gains exactly the same
    // bases in the neighbouring context block and no template's spelled locus changes.
    for (std::string& value : before.sequences) value += prefix_bases;
    for (std::string& value : after.sequences) value = suffix_bases + value;
    moved_bases += prefix + suffix;
    outcome.prefix = prefix;
    outcome.suffix = suffix;
    return outcome;
}

} // namespace

LocusInputResult build_locus_input(
    const Graph& graph,
    const std::vector<Bubble>& bubbles,
    const LocusInputOptions& options) {
    if (options.locus_id.empty()) {
        throw std::invalid_argument("genotype locus ID is empty");
    }
    if (graph.paths.empty()) {
        throw std::invalid_argument("genotype locus graph has no panel paths");
    }
    if (bubbles.empty()) {
        throw std::invalid_argument("genotype locus has no bubbles to build a block chain from");
    }
    for (std::size_t i = 1; i < bubbles.size(); ++i) {
        if (bubbles[i].id < bubbles[i - 1].id) {
            throw std::invalid_argument(
                "genotype locus bubbles are not in ascending id order; the block chain assumes "
                "reference order, which requires a reference-sorted graph");
        }
    }

    LocusInputResult result;
    LocusInputReport& report = result.report;

    const std::vector<const Bubble*> kept =
        drop_contained(bubbles, report.contained_bubbles_dropped);
    const std::vector<ChainBlock> chain =
        build_chain(kept, report.abutting_boundaries_trimmed);

    // Block identifiers and per-template refusal reasons are needed before any allele is emitted:
    // the frame decision below names the blocks a template failed to cover.
    std::vector<std::string> block_ids(chain.size());
    for (std::size_t bi = 0; bi < chain.size(); ++bi) {
        block_ids[bi] = (chain[bi].kind == ChainBlockKind::Bubble ? "sv" : "ctx") +
                        std::string(":") + std::to_string(bi);
    }
    std::vector<std::string> partial_reason_by_template(graph.paths.size());

    std::vector<std::string> ordered_boundaries;
    for (const Bubble* bubble : kept) {
        ordered_boundaries.push_back(left_boundary(*bubble));
        ordered_boundaries.push_back(right_boundary(*bubble));
    }

    std::vector<BubblePathIndex> indexes;
    indexes.reserve(graph.paths.size());
    for (const PathRecord& path : graph.paths) {
        indexes.push_back(build_bubble_path_index(path));
    }

    // Decide every template's frame before any allele is emitted, so a refused template contributes
    // no allele record and cannot widen a block it was never part of.
    std::vector<std::optional<ChainOrientation>> orientations(graph.paths.size());
    for (std::size_t pi = 0; pi < graph.paths.size(); ++pi) {
        orientations[pi] = template_orientation(graph.paths[pi], ordered_boundaries);
        if (!orientations[pi].has_value()) {
            report.refused_templates.push_back(graph.paths[pi].name);
            partial_reason_by_template[pi] = "no consistent chain frame";
            PartialTemplateReport partial;
            partial.name = graph.paths[pi].name;
            partial.reason = partial_reason_by_template[pi];
            report.partial_templates.push_back(std::move(partial));
        } else if (*orientations[pi] == ChainOrientation::Forward) {
            ++report.forward_templates;
        } else {
            ++report.reverse_templates;
        }
    }

    result.input.locus_id = options.locus_id;
    result.input.chain_orientation = ChainOrientation::Forward;
    result.input.blocks.reserve(chain.size());
    // Every span of the locus is a block, so nothing is left for a shared invariant context. The
    // frame is kept because a locus whose between-bubble spans really are invariant across the whole
    // panel would collapse those blocks to one allele each, which is the same model.
    result.input.invariant_segments.assign(chain.size() + 1, std::string());

    std::vector<std::vector<std::string>> allele_source_by_template(
        graph.paths.size(), std::vector<std::string>());

    // Pass one: walk every block once and record what each template did there, emitting no allele
    // yet. Whether a template may contribute at all cannot be decided until the whole chain has been
    // walked, and the divergence trim below moves bases between neighbours, so nothing can be
    // finalised earlier either.
    std::vector<std::map<std::string, std::string>> sequence_by_source(chain.size());
    for (std::size_t bi = 0; bi < chain.size(); ++bi) {
        for (std::size_t pi = 0; pi < graph.paths.size(); ++pi) {
            if (!orientations[pi].has_value()) continue;
            const BlockWalk walk = walk_block(graph, graph.paths[pi], indexes[pi], chain[bi]);
            std::string source_id;
            if (!walk.traverses) {
                source_id = "bypass";
            } else if (walk.signature.empty()) {
                source_id = "empty";
            } else {
                source_id = "walk:" + walk.signature;
            }
            sequence_by_source[bi].emplace(source_id, walk.sequence);
            allele_source_by_template[pi].push_back(source_id);
        }
    }

    // Pass two: decide each template's frame. A block a template does not walk is an admissible
    // deletion only when the template is observed on BOTH sides of it -- the path brackets the
    // missing region, so that region really is absent from this haplotype. With no observation on one
    // side the path simply stopped: the block is uncovered, the template never said anything about
    // it, and the bases the path does carry there belong to no block at all. Uncovered sequence must
    // never become a zero-length deletion allele, so such a template's frame is partial and it
    // contributes nothing -- no allele, no bypass count, and so no representation-dependence flag.
    std::vector<bool> complete(graph.paths.size(), false);
    for (std::size_t pi = 0; pi < graph.paths.size(); ++pi) {
        if (!orientations[pi].has_value()) continue;
        const std::vector<std::string>& sources = allele_source_by_template[pi];
        std::vector<std::string> uncovered;
        bool missing_left = false;
        bool missing_right = false;
        for (std::size_t bi = 0; bi < sources.size(); ++bi) {
            if (sources[bi] != "bypass") continue;
            bool observed_left = false;
            bool observed_right = false;
            for (std::size_t j = 0; j < bi; ++j) {
                if (sources[j] != "bypass") { observed_left = true; break; }
            }
            for (std::size_t j = bi + 1; j < sources.size(); ++j) {
                if (sources[j] != "bypass") { observed_right = true; break; }
            }
            if (observed_left && observed_right) continue;   // bracketed on both sides: a deletion
            uncovered.push_back(block_ids[bi]);
            if (!observed_left) missing_left = true;
            if (!observed_right) missing_right = true;
        }
        complete[pi] = uncovered.empty();
        if (complete[pi]) continue;

        std::string names;
        for (const std::string& block : uncovered) {
            if (!names.empty()) names += ", ";
            names += block;
        }
        PartialTemplateReport partial;
        partial.name = graph.paths[pi].name;
        partial.uncovered_blocks = uncovered;
        partial.missing_terminal_anchor =
            missing_left && missing_right ? "both" : (missing_left ? "left" : "right");
        partial.reason = "locus_coverage_partial: missing " + partial.missing_terminal_anchor +
                         "-terminal coverage at " + names;
        partial_reason_by_template[pi] = partial.reason;
        report.partial_templates.push_back(std::move(partial));
        report.refused_templates.push_back(graph.paths[pi].name);
    }

    // Pass three: emit alleles, from templates with a complete frame only.
    std::vector<BlockSequences> block_sequences(chain.size());
    std::vector<std::map<std::string, std::size_t>> slot_by_source(chain.size());
    std::vector<std::size_t> bypass_counts(chain.size(), 0);
    for (std::size_t bi = 0; bi < chain.size(); ++bi) {
        BlockSequences& out = block_sequences[bi];
        for (std::size_t pi = 0; pi < graph.paths.size(); ++pi) {
            if (!complete[pi]) continue;
            const std::string& source_id = allele_source_by_template[pi][bi];
            if (source_id == "bypass") {
                // A bypass that survives the frame decision is bracketed on both sides: a genuine
                // deletion, which is a well-defined state and must stay representable.
                out.has_bypass = true;
                ++bypass_counts[bi];
            }
            const auto found = slot_by_source[bi].find(source_id);
            if (found == slot_by_source[bi].end()) {
                slot_by_source[bi].emplace(source_id, out.source_ids.size());
                out.source_ids.push_back(source_id);
                out.sequences.push_back(sequence_by_source[bi].at(source_id));
            }
        }
    }

    std::size_t moved_bases = 0;
    std::vector<std::size_t> trimmed_prefix(chain.size(), 0);
    std::vector<std::size_t> trimmed_suffix(chain.size(), 0);
    std::vector<bool> dependent(chain.size(), false);
    for (std::size_t bi = 0; bi < chain.size(); ++bi) {
        if (chain[bi].kind != ChainBlockKind::Bubble) continue;
        if (bi == 0 || bi + 1 >= chain.size()) continue;
        const TrimOutcome outcome = trim_variant_block_to_divergence(
            block_sequences[bi - 1], block_sequences[bi], block_sequences[bi + 1], moved_bases);
        trimmed_prefix[bi] = outcome.prefix;
        trimmed_suffix[bi] = outcome.suffix;
        if (outcome.representation_dependent) {
            // Both neighbours share this block's node-defined boundary.
            dependent[bi - 1] = true;
            dependent[bi] = true;
            dependent[bi + 1] = true;
        }
    }
    report.boundary_bases_moved_to_context = moved_bases;

    for (std::size_t bi = 0; bi < chain.size(); ++bi) {
        const bool is_variant = chain[bi].kind == ChainBlockKind::Bubble;
        const BlockSequences& sequences = block_sequences[bi];

        GenotypeBlockInput out;
        out.kind = is_variant ? GenotypeBlockKind::Variant : GenotypeBlockKind::Context;
        out.block_id = block_ids[bi];

        LocusBlockReport block_report;
        block_report.block_id = out.block_id;
        block_report.kind = out.kind;
        block_report.bubble_id = is_variant ? chain[bi].bubble_id : 0;
        block_report.source_node = chain[bi].source;
        block_report.sink_node = chain[bi].sink;
        block_report.bypass_templates = bypass_counts[bi];
        block_report.trimmed_prefix_bp = trimmed_prefix[bi];
        block_report.trimmed_suffix_bp = trimmed_suffix[bi];
        block_report.representation_dependent = dependent[bi];
        if (dependent[bi]) report.representation_dependent_blocks.push_back(block_report.block_id);

        std::set<std::string> distinct_sequences;
        for (std::size_t a = 0; a < sequences.source_ids.size(); ++a) {
            out.alleles.push_back({sequences.source_ids[a], sequences.sequences[a]});
            distinct_sequences.insert(sequences.sequences[a]);
            const std::size_t size = sequences.sequences[a].size();
            if (a == 0 || size < block_report.min_bp) block_report.min_bp = size;
            block_report.max_bp = std::max(block_report.max_bp, size);
        }
        block_report.source_alleles = out.alleles.size();
        block_report.canonical_alleles = distinct_sequences.size();

        report.total_source_alleles += block_report.source_alleles;
        report.total_canonical_alleles += block_report.canonical_alleles;
        if (is_variant) {
            ++report.variant_blocks;
        } else {
            ++report.context_blocks;
        }
        report.blocks.push_back(std::move(block_report));
        result.input.blocks.push_back(std::move(out));
    }

    // Classify every block a template does not simply traverse, now that the whole chain is known.
    // A block the template does not walk is an admissible deletion only when the template is observed
    // on BOTH sides of it: the path brackets the missing region, so the region really is absent from
    // that haplotype. With no observation on one side the path stops short, and the block is
    // uncovered -- the template never said anything about it. Reporting only, for now.
    for (std::size_t pi = 0; pi < graph.paths.size(); ++pi) {
        if (!orientations[pi].has_value()) continue;
        const std::vector<std::string>& sources = allele_source_by_template[pi];
        for (std::size_t bi = 0; bi < sources.size() && bi < report.blocks.size(); ++bi) {
            if (sources[bi] != "bypass" && sources[bi] != "empty") continue;
            TemplateBlockCoverage row;
            row.template_name = graph.paths[pi].name;
            row.block_id = report.blocks[bi].block_id;
            if (sources[bi] == "empty") {
                row.coverage = BlockCoverage::EmptyInterior;
            } else {
                bool observed_left = false;
                bool observed_right = false;
                for (std::size_t j = 0; j < bi; ++j) {
                    if (sources[j] != "bypass") { observed_left = true; break; }
                }
                for (std::size_t j = bi + 1; j < sources.size(); ++j) {
                    if (sources[j] != "bypass") { observed_right = true; break; }
                }
                if (observed_left && observed_right) {
                    row.coverage = BlockCoverage::BridgedBypass;
                } else {
                    row.coverage = BlockCoverage::Uncovered;
                    row.uncovered_on_left = !observed_left;
                }
            }
            report.non_traversing_coverage.push_back(std::move(row));
        }
    }

    // Curation is applied at the one boundary where graph paths become panel templates. A curated
    // template contributes neither a panel row nor an authoritative sequence, so the round trip
    // verifies exactly the panel the model scores. It cannot be done later: a partial frame fails
    // model construction closed, which is the behaviour that must stay intact.
    const std::set<std::string> excluded(options.excluded_templates.begin(),
                                         options.excluded_templates.end());
    std::set<std::string> excluded_matched;

    for (std::size_t pi = 0; pi < graph.paths.size(); ++pi) {
        PanelTemplateInput out;
        out.name = graph.paths[pi].name;
        if (excluded.count(out.name) != 0) {
            excluded_matched.insert(out.name);
            LocusInputReport::CuratedExclusion record;
            record.template_name = out.name;
            record.reason = complete[pi] ? std::string("complete frame, excluded by curation")
                                         : partial_reason_by_template[pi];
            report.curated_exclusions.push_back(std::move(record));
            // The frame pass above already recorded this template as refused, because it runs before
            // curation does. Leaving both records reports one template twice under two headings -- one
            // refused AND one curated -- from which a reader concludes the panel lost two rows, and
            // worse, that a partial frame survived into a model that in fact never saw it.
            report.refused_templates.erase(
                std::remove(report.refused_templates.begin(), report.refused_templates.end(), out.name),
                report.refused_templates.end());
            report.partial_templates.erase(
                std::remove_if(report.partial_templates.begin(), report.partial_templates.end(),
                               [&out](const PartialTemplateReport& partial) {
                                   return partial.name == out.name;
                               }),
                report.partial_templates.end());
            continue;
        }
        if (!complete[pi]) {
            out.frame_status = TemplateFrameStatus::Partial;
            out.partial_reason = partial_reason_by_template[pi];
            result.input.templates.push_back(std::move(out));
            continue;
        }
        out.frame_status = TemplateFrameStatus::Complete;
        out.source_orientation = *orientations[pi];
        out.allele_source_ids = std::move(allele_source_by_template[pi]);
        result.input.templates.push_back(std::move(out));

        AuthoritativePanelSequence authoritative;
        authoritative.template_name = graph.paths[pi].name;
        authoritative.source_frame_bases = spell_path_steps_sequence(graph, graph.paths[pi].steps);
        result.authoritative.push_back(std::move(authoritative));
    }

    // A manifest naming a template this graph does not have is refused, never ignored. Silently
    // curating nothing is the failure mode that would let a stale manifest travel with a changed graph
    // and quietly stop excluding the template it was written for.
    if (excluded_matched.size() != excluded.size()) {
        std::string missing;
        for (const std::string& name : excluded) {
            if (excluded_matched.count(name) != 0) continue;
            if (!missing.empty()) missing += ", ";
            missing += name;
        }
        throw std::invalid_argument(
            "curation names templates that are not paths of this graph: " + missing +
            "; a manifest that has drifted from its graph must be corrected, not ignored");
    }
    if (!report.curated_exclusions.empty()) {
        // Ordered by the set, so the digest depends on which templates were curated and not on the
        // order they were listed in.
        Digest digest;
        digest.add(excluded.size());
        for (const std::string& name : excluded) digest.add(name);
        report.curation_digest = digest.hex();
    }
    return result;
}

const char* block_coverage_label(BlockCoverage coverage, bool uncovered_on_left) {
    switch (coverage) {
        case BlockCoverage::Traversing: return "traversing";
        case BlockCoverage::EmptyInterior: return "empty_between_boundaries";
        case BlockCoverage::BridgedBypass: return "bridged_bypass";
        case BlockCoverage::Uncovered: return uncovered_on_left ? "uncovered_left" : "uncovered_right";
    }
    return "unknown";
}

void require_representation_independent(const LocusInputReport& report) {
    if (report.representation_dependent_blocks.empty()) return;
    std::string names;
    for (const std::string& block : report.representation_dependent_blocks) {
        if (!names.empty()) names += ", ";
        names += block;
    }
    throw std::invalid_argument(
        "block boundaries depend on graph representation at " + names +
        ": a bypass allele prevents the divergence trim, so a node split would change the model; "
        "version one refuses such a locus");
}

std::optional<PathProjection> project_path_onto_locus(
    const Graph& graph,
    const std::vector<Bubble>& bubbles,
    const LocusInputReport& report,
    const std::string& path_name) {
    const PathRecord* path = nullptr;
    for (const PathRecord& candidate : graph.paths) {
        if (candidate.name == path_name) {
            path = &candidate;
            break;
        }
    }
    if (path == nullptr) {
        throw std::invalid_argument("projection path is not in the graph: " + path_name);
    }
    if (bubbles.empty()) {
        throw std::invalid_argument("projection needs the bubbles the model was built from");
    }
    std::size_t dropped = 0;
    std::size_t trimmed = 0;
    const std::vector<const Bubble*> kept = drop_contained(bubbles, dropped);
    const std::vector<ChainBlock> chain = build_chain(kept, trimmed);
    if (chain.size() != report.blocks.size()) {
        throw std::invalid_argument(
            "projection chain and model report disagree on block count; they were not built from the "
            "same bubbles");
    }
    std::vector<std::string> ordered_boundaries;
    for (const Bubble* bubble : kept) {
        ordered_boundaries.push_back(left_boundary(*bubble));
        ordered_boundaries.push_back(right_boundary(*bubble));
    }
    const std::optional<ChainOrientation> orientation = template_orientation(*path, ordered_boundaries);
    if (!orientation.has_value()) return std::nullopt;

    const BubblePathIndex index = build_bubble_path_index(*path);
    std::vector<BlockWalk> walks;
    walks.reserve(chain.size());
    for (const ChainBlock& block : chain) walks.push_back(walk_block(graph, *path, index, block));

    PathProjection projection;
    projection.orientation = *orientation;
    projection.blocks.resize(chain.size());
    for (std::size_t bi = 0; bi < chain.size(); ++bi) projection.blocks[bi] = walks[bi].sequence;

    for (std::size_t bi = 0; bi < chain.size(); ++bi) {
        const std::size_t prefix = report.blocks[bi].trimmed_prefix_bp;
        const std::size_t suffix = report.blocks[bi].trimmed_suffix_bp;
        if (prefix + suffix == 0 || !walks[bi].traverses || !projection.blocks[bi].has_value()) continue;
        std::string& own = *projection.blocks[bi];
        if (own.size() < prefix + suffix) {
            projection.blocks[bi].reset();
            continue;
        }
        const std::string moved_prefix = own.substr(0, prefix);
        const std::string moved_suffix = own.substr(own.size() - suffix);
        own = own.substr(prefix, own.size() - prefix - suffix);
        if (bi > 0 && projection.blocks[bi - 1].has_value()) *projection.blocks[bi - 1] += moved_prefix;
        if (bi + 1 < chain.size() && projection.blocks[bi + 1].has_value()) {
            *projection.blocks[bi + 1] = moved_suffix + *projection.blocks[bi + 1];
        }
    }
    return projection;
}

std::string model_fingerprint(const PreparedGenotypeModel& model) {
    Digest digest;
    digest.add(model.locus_id());
    digest.add(static_cast<std::size_t>(model.chain_orientation()));
    digest.add(model.blocks().size());
    for (const PreparedGenotypeBlock& block : model.blocks()) {
        // Block IDs carry chain positions, not node identifiers, so they survive a node split. Allele
        // sequences rather than source aliases are hashed: a split changes every walk signature while
        // leaving the statistical model identical.
        digest.add(block.block_id);
        digest.add(static_cast<std::size_t>(block.kind));
        digest.add(block.alleles.size());
        for (const CanonicalBlockAllele& allele : block.alleles) {
            digest.add(allele.sequence);
        }
    }
    for (const std::string& context : model.invariant_segments()) {
        digest.add(context);
    }
    digest.add(model.panel().template_count());
    for (std::size_t t = 0; t < model.panel().template_count(); ++t) {
        digest.add(model.panel().template_names[t]);
        digest.add(static_cast<std::size_t>(model.template_source_orientations()[t]));
        for (const std::size_t allele_id : model.panel().template_alleles[t]) {
            digest.add(allele_id);
        }
    }
    return digest.hex();
}

} // namespace panvar
