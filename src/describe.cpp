#include "panvar/describe.hpp"

#include "panvar/syncmer.hpp"

#include "panvar/bubble_path.hpp"
#include "panvar/variant_nodes.hpp"
#include "panvar/cli_utils.hpp"
#include "panvar/cosigt.hpp"
#include "panvar/gfa.hpp"
#include "panvar/graph_utils.hpp"
#include "panvar/ref_path.hpp"
#include "panvar/output.hpp"
#include "panvar/parallel.hpp"

#include <algorithm>
#include <functional>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <system_error>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <mutex>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include <zlib.h>
#include <cstdio>
#include <string_view>

namespace panvar {
namespace {

struct KmerStats {
    std::size_t present_paths = 0;
    std::uint64_t total_count = 0;
    std::uint32_t min_nonzero_count = std::numeric_limits<std::uint32_t>::max();
    std::uint32_t max_count = 0;
};

struct PathMeta {
    std::size_t path_index = 0;
    std::string path_name;
    std::string sample;
    std::string haplotype;
    std::size_t length_bp = 0;
};

struct SequenceSegment {
    std::size_t start = 0;
    std::size_t end = 0;
    std::string node_id;
};

struct BubblePathSequence {
    std::string sequence;
    std::vector<SequenceSegment> segments;
    bool complete = false;
};

struct KmerPathCount {
    std::uint32_t count = 0;
    std::unordered_set<std::string> nodes;
};

struct SparseKmerEntry {
    std::size_t feature_id = 0;
    std::uint32_t count = 0;
    std::vector<std::string> nodes;
};

struct BubbleDescribeResult {
    std::size_t paths = 0;
    std::size_t features = 0;            // k-mers kept
    std::size_t features_total = 0;      // k-mer candidates before filter
    bool matrix_written = false;
    std::string status = "ok";
    std::string matrix_reason = ".";
    std::string feature_map_path = ".";
    std::string matrix_path = ".";
    std::string counts_jsonl_path = ".";
    std::size_t node_features = 0;       // nodes kept
    std::size_t node_features_total = 0; // node candidates before filter
    std::size_t edge_features = 0;       // edges kept
    std::size_t edge_features_total = 0; // edge candidates before filter
    bool graph_matrix_written = false;
    std::string graph_feature_map_path = ".";
    std::string graph_matrix_path = ".";
    std::vector<std::uint32_t> traversing_paths;  // paths (graph index) crossing the site
    std::unordered_map<std::uint64_t, std::vector<std::string>> kmer_nodes;  // kept k-mer code -> its node set
};

class GzipWriter {
public:
    explicit GzipWriter(const std::string& path)
        : path_(path), file_(gzopen(path.c_str(), "wb6")) {
        if (file_ == nullptr) {
            throw std::runtime_error("Failed to write gzip output: " + path);
        }
    }

    GzipWriter(const GzipWriter&) = delete;
    GzipWriter& operator=(const GzipWriter&) = delete;

    ~GzipWriter() {
        if (file_ != nullptr) {
            gzclose(file_);
        }
    }

    void write(const std::string& text) {
        if (text.empty()) {
            return;
        }
        const int written = gzwrite(file_, text.data(), static_cast<unsigned int>(text.size()));
        if (written != static_cast<int>(text.size())) {
            throw std::runtime_error("Failed to write gzip output: " + path_);
        }
    }

    void close() {
        if (file_ != nullptr) {
            if (gzclose(file_) != Z_OK) {
                file_ = nullptr;
                throw std::runtime_error("Failed to close gzip output: " + path_);
            }
            file_ = nullptr;
        }
    }

private:
    std::string path_;
    gzFile file_ = nullptr;
};

std::string tsv_sanitize(std::string value) {
    for (char& c : value) {
        if (c == '\t' || c == '\n' || c == '\r') {
            c = ' ';
        }
    }
    return value;
}

std::string json_escape(const std::string& value) {
    std::string out;
    out.reserve(value.size() + 8);
    for (const unsigned char c : value) {
        switch (c) {
            case '"':
                out += "\\\"";
                break;
            case '\\':
                out += "\\\\";
                break;
            case '\b':
                out += "\\b";
                break;
            case '\f':
                out += "\\f";
                break;
            case '\n':
                out += "\\n";
                break;
            case '\r':
                out += "\\r";
                break;
            case '\t':
                out += "\\t";
                break;
            default:
                if (c < 0x20) {
                    constexpr char hex[] = "0123456789abcdef";
                    out += "\\u00";
                    out.push_back(hex[(c >> 4) & 0xF]);
                    out.push_back(hex[c & 0xF]);
                } else {
                    out.push_back(static_cast<char>(c));
                }
                break;
        }
    }
    return out;
}

std::pair<std::string, std::string> parse_sample_haplotype(const std::string& path_name) {
    const std::size_t first_hash = path_name.find('#');
    if (first_hash == std::string::npos) {
        return {path_name, "."};
    }
    const std::string sample = path_name.substr(0, first_hash);
    const std::size_t second_hash = path_name.find('#', first_hash + 1);
    if (second_hash == std::string::npos) {
        return {sample, path_name.substr(first_hash + 1)};
    }
    return {sample, path_name.substr(first_hash + 1, second_hash - first_hash - 1)};
}

std::size_t effective_syncmer_s(const DescribeOptions& options) {
    if (options.syncmer_s != 0) {
        return options.syncmer_s;
    }
    return default_syncmer_s(options.kmer_size);
}

std::vector<std::size_t> select_feature_occurrence_indices(
    const std::vector<KmerOccurrence>& occurrences,
    const DescribeOptions& options) {

    std::vector<std::size_t> out;
    if (occurrences.empty()) {
        return out;
    }

    if (options.feature_mode == DescribeFeatureMode::AllKmers) {
        out.reserve(occurrences.size());
        for (std::size_t i = 0; i < occurrences.size(); ++i) {
            out.push_back(i);
        }
        return out;
    }

    // Syncmer (default): keep closed syncmers only.
    const std::size_t s = effective_syncmer_s(options);
    out.reserve(occurrences.size() / 8 + 1);
    for (std::size_t i = 0; i < occurrences.size(); ++i) {
        if (is_closed_syncmer(occurrences[i].code, options.kmer_size, s)) {
            out.push_back(i);
        }
    }
    return out;
}

void count_canonical_kmers(
    const std::string& sequence,
    std::size_t k,
    std::unordered_map<std::uint64_t, std::uint32_t>& counts) {

    counts.clear();
    if (k == 0 || k > 31 || sequence.size() < k) {
        return;
    }

    const std::size_t possible = sequence.size() - k + 1;
    counts.reserve(possible);

    const std::uint64_t mask = (1ULL << (2 * k)) - 1ULL;
    const std::size_t rc_shift = 2 * (k - 1);
    std::uint64_t fwd = 0;
    std::uint64_t rev = 0;
    std::size_t filled = 0;

    for (const char c : sequence) {
        const int b = base_bits(c);
        if (b < 0) {
            fwd = 0;
            rev = 0;
            filled = 0;
            continue;
        }
        const std::uint64_t ub = static_cast<std::uint64_t>(b);
        const std::uint64_t cb = static_cast<std::uint64_t>(3 - b);
        fwd = ((fwd << 2) | ub) & mask;
        rev = (rev >> 2) | (cb << rc_shift);
        if (filled < k) {
            ++filled;
        }
        if (filled >= k) {
            const std::uint64_t canonical = std::min(fwd, rev);
            auto& count = counts[canonical];
            if (count != std::numeric_limits<std::uint32_t>::max()) {
                ++count;
            }
        }
    }
}

void count_selected_features(
    const std::string& sequence,
    const DescribeOptions& options,
    std::unordered_map<std::uint64_t, std::uint32_t>& counts) {

    counts.clear();
    if (options.feature_mode == DescribeFeatureMode::AllKmers) {
        count_canonical_kmers(sequence, options.kmer_size, counts);
        return;
    }

    const std::vector<KmerOccurrence> occurrences =
        collect_canonical_kmer_occurrences(sequence, options.kmer_size);
    const std::vector<std::size_t> selected = select_feature_occurrence_indices(occurrences, options);
    counts.reserve(selected.size());
    for (const std::size_t idx : selected) {
        auto& count = counts[occurrences[idx].code];
        if (count != std::numeric_limits<std::uint32_t>::max()) {
            ++count;
        }
    }
}

void add_window_nodes(
    const std::vector<SequenceSegment>& segments,
    std::size_t window_start,
    std::size_t window_end,
    std::unordered_set<std::string>& nodes,
    std::size_t& segment_hint) {

    while (segment_hint < segments.size() && segments[segment_hint].end <= window_start) {
        ++segment_hint;
    }
    for (std::size_t i = segment_hint; i < segments.size() && segments[i].start < window_end; ++i) {
        nodes.insert(segments[i].node_id);
    }
}

void count_canonical_kmers_with_nodes(
    const BubblePathSequence& path_sequence,
    std::size_t k,
    std::unordered_map<std::uint64_t, KmerPathCount>& counts) {

    counts.clear();
    const std::string& sequence = path_sequence.sequence;
    if (k == 0 || k > 31 || sequence.size() < k) {
        return;
    }

    const std::size_t possible = sequence.size() - k + 1;
    counts.reserve(possible);

    const std::uint64_t mask = (1ULL << (2 * k)) - 1ULL;
    const std::size_t rc_shift = 2 * (k - 1);
    std::uint64_t fwd = 0;
    std::uint64_t rev = 0;
    std::size_t filled = 0;
    std::size_t segment_hint = 0;

    for (std::size_t pos = 0; pos < sequence.size(); ++pos) {
        const int b = base_bits(sequence[pos]);
        if (b < 0) {
            fwd = 0;
            rev = 0;
            filled = 0;
            continue;
        }
        const std::uint64_t ub = static_cast<std::uint64_t>(b);
        const std::uint64_t cb = static_cast<std::uint64_t>(3 - b);
        fwd = ((fwd << 2) | ub) & mask;
        rev = (rev >> 2) | (cb << rc_shift);
        if (filled < k) {
            ++filled;
        }
        if (filled >= k) {
            const std::uint64_t canonical = std::min(fwd, rev);
            auto& entry = counts[canonical];
            if (entry.count != std::numeric_limits<std::uint32_t>::max()) {
                ++entry.count;
            }
            const std::size_t window_start = pos + 1 - k;
            const std::size_t window_end = pos + 1;
            add_window_nodes(path_sequence.segments, window_start, window_end, entry.nodes, segment_hint);
        }
    }
}

void count_selected_features_with_nodes(
    const BubblePathSequence& path_sequence,
    const DescribeOptions& options,
    std::unordered_map<std::uint64_t, KmerPathCount>& counts) {

    counts.clear();
    if (options.feature_mode == DescribeFeatureMode::AllKmers) {
        count_canonical_kmers_with_nodes(path_sequence, options.kmer_size, counts);
        return;
    }

    const std::vector<KmerOccurrence> occurrences =
        collect_canonical_kmer_occurrences(path_sequence.sequence, options.kmer_size);
    const std::vector<std::size_t> selected = select_feature_occurrence_indices(occurrences, options);
    counts.reserve(selected.size());

    std::size_t segment_hint = 0;
    for (const std::size_t idx : selected) {
        const KmerOccurrence& occurrence = occurrences[idx];
        auto& entry = counts[occurrence.code];
        if (entry.count != std::numeric_limits<std::uint32_t>::max()) {
            ++entry.count;
        }
        add_window_nodes(
            path_sequence.segments,
            occurrence.start,
            occurrence.start + options.kmer_size,
            entry.nodes,
            segment_hint);
    }
}

void update_kmer_stats(
    const std::unordered_map<std::uint64_t, std::uint32_t>& counts,
    std::unordered_map<std::uint64_t, KmerStats>& stats) {

    for (const auto& [code, count] : counts) {
        auto& s = stats[code];
        ++s.present_paths;
        s.total_count += count;
        s.min_nonzero_count = std::min(s.min_nonzero_count, count);
        s.max_count = std::max(s.max_count, count);
    }
}

// Whether a feature is informative: one whose count varies among its carriers always is (copy number);
// otherwise it needs min(present, absent) > min_paths. A feature no path observes is not.
bool feature_passes_filter(
    std::size_t present_paths,
    std::uint32_t min_nonzero_count,
    std::uint32_t max_count,
    std::size_t path_count,
    std::size_t min_paths) {

    if (path_count == 0) {
        return false;
    }
    if (min_nonzero_count != max_count) {
        return true;
    }
    const std::size_t absent = path_count - present_paths;
    const std::size_t minor = std::min(present_paths, absent);
    return minor > min_paths;
}

std::vector<std::uint64_t> select_discriminative_features(
    const std::unordered_map<std::uint64_t, KmerStats>& stats,
    std::size_t path_count,
    std::size_t min_paths) {

    std::vector<std::uint64_t> features;
    features.reserve(stats.size());
    for (const auto& [code, s] : stats) {
        if (feature_passes_filter(s.present_paths, s.min_nonzero_count, s.max_count, path_count, min_paths)) {
            features.push_back(code);
        }
    }
    std::sort(features.begin(), features.end());
    return features;
}

std::string feature_name(std::size_t feature_id) {
    return "K" + std::to_string(feature_id);
}

std::unordered_map<std::uint64_t, std::size_t> make_feature_id_map(
    const std::vector<std::uint64_t>& features) {

    std::unordered_map<std::uint64_t, std::size_t> out;
    out.reserve(features.size() * 2 + 1);
    for (std::size_t i = 0; i < features.size(); ++i) {
        out[features[i]] = i + 1;
    }
    return out;
}

// Node and edge dosage: per path, how often it walks each node and each oriented edge of the site. A
// count, not a copy-number call: a node revisited elsewhere also counts more than once.

struct GraphFeatureStat {
    std::size_t present_paths = 0;
    std::uint64_t total_count = 0;
    std::uint32_t min_nonzero_count = std::numeric_limits<std::uint32_t>::max();
    std::uint32_t max_count = 0;
};

// Ordered maps give a deterministic feature order (stable N*/E* ids).
struct GraphDosage {
    std::map<std::string, std::uint32_t> node_counts;
    std::map<std::string, std::uint32_t> edge_counts;
};

// One path's counts at a site, taken once and read by every output of the site.
struct PathWork {
    std::vector<std::pair<std::uint64_t, std::uint32_t>> kmers;   // by code
    GraphDosage dosage;
};

// The count of k-mer `code` in a path's counts, 0 when absent.
std::uint32_t kmer_count(const PathWork& w, std::uint64_t code) {
    const auto it = std::lower_bound(w.kmers.begin(), w.kmers.end(), std::make_pair(code, std::uint32_t{0}));
    return it != w.kmers.end() && it->first == code ? it->second : 0;
}

std::string node_feature_name(std::size_t feature_id) {
    return "N" + std::to_string(feature_id);
}

std::string edge_feature_name(std::size_t feature_id) {
    return "E" + std::to_string(feature_id);
}

// A path's walk across a bubble: every crossing, on the bubble's own nodes only. `run_start[i]` marks a
// step the walk did not reach from the step before it -- a crossing's source, or a node after a stretch
// outside the bubble: no edge is counted into it and no k-mer spans it. `crossing_start[i]` marks the
// first step (the source) of each crossing, whose last step is its sink.
struct BubbleWalk {
    std::vector<PathStep> steps;
    std::vector<char> run_start;
    std::vector<char> crossing_start;
    std::vector<std::size_t> path_pos;   // each step's index in the path
    std::vector<char> forward;           // whether its crossing reads the path forward
    bool empty() const { return steps.empty(); }
};

// How much of one step's sequence and node a site's scope keeps.
struct KeepSpan {
    bool full = false;        // a variant node: all of it
    bool boundary = false;    // the site's source or sink: kept on every crossing
    std::size_t prefix = 0;   // bases kept at the node's start (a variant node follows downstream)
    std::size_t suffix = 0;   // bases kept at the node's end (a variant node precedes upstream)
    bool no_node = false;     // a visit of a boundary another described site counts: no node dosage here
    bool any() const { return full || boundary || prefix != 0 || suffix != 0; }
};

// What a site keeps of one node under --variant-nodes, over all its walks, in the node's own orientation:
// all of it, its first `head` and last `tail` bases, and whether it is a boundary.
struct NodeKeep {
    bool full = false;
    bool boundary = false;
    std::size_t head = 0;
    std::size_t tail = 0;
};
using SiteKeep = std::unordered_map<std::string, NodeKeep>;

// Per path, the sorted path positions at which a site's crossings visit the boundaries it shares and counts.
using SharedVisits = std::vector<std::vector<std::size_t>>;

// What of a site's walks is described: with --variant-nodes, its variant nodes and the flank; and a
// boundary shared with another described site (one's sink is the next one's source) is counted at one
// of them, its owner, on every visit the owner's crossings make. `shared_nodes` and `shared_kmers` map
// such a boundary to its owner's visits: on those, the node loses its dosage here (`shared_nodes`), or
// keeps only the k-1 bases next to this site's interior (`shared_kmers`), so the k-mers across the
// junction stay and those inside the node do not. A visit only this site's crossings make counts here.
struct SiteScope {
    const std::unordered_set<std::string>* variant_nodes = nullptr;
    // With --variant-nodes, the bases kept on every walk: the union of each walk's own flank, so a node
    // or k-mer is masked alike on every walk and its count stays a visit count. A boundary shared for
    // k-mers keeps the union over the sites sharing it, so both spell it alike.
    const SiteKeep* site_keep = nullptr;
    std::size_t flank_bp = 0;
    std::unordered_map<std::string, const SharedVisits*> shared_nodes;
    std::unordered_map<std::string, const SharedVisits*> shared_kmers;
    // This site's shared boundaries, owned or not, for nodes and for k-mers, whose visits are split
    // between the two sites (some only the other site's crossings make). Only the sum is a visit count,
    // so the node's row and the k-mers inside it are pooled here whether or not they pass the filter.
    std::unordered_set<std::string> pooled_nodes;
    std::unordered_set<std::string> pooled_kmer_nodes;
    // Shared boundaries counted here: one shorter than k-1 is spelled with the bases beyond it, so the
    // k-mers from this site's interior across it into the next site's are counted, here.
    std::unordered_set<std::string> extend_kmers;
    std::size_t kmer_size = 0;
    bool masked() const { return variant_nodes != nullptr || !shared_nodes.empty() || !shared_kmers.empty(); }
};

// Oriented steps (node id, then '+' or '-') that some walk of a site takes just before a variant node
// (`before`) or just after one (`after`). A walk stepping straight from one of the first to one of the
// second skips a record's nodes, as a deletion's carriers or an insertion's non-carriers do, and the
// flank is kept around that junction too.
struct VariantJunctions {
    std::unordered_set<std::string> before;
    std::unordered_set<std::string> after;
};

// Per-step keep spec for --variant-nodes and the flank (defined below). Node and edge dosage use it per
// node, counted or not; k-mers use the base offsets.
std::vector<KeepSpan> variant_keep_mask(
    const Graph& graph,
    const BubbleWalk& walk,
    const std::unordered_set<std::string>& variant_nodes,
    std::size_t flank_bp,
    const VariantJunctions* junctions);

// Keep the site's source and sink on every crossing (defined below).
void keep_bubble_boundaries(const Graph& graph, const BubbleWalk& walk,
                            std::size_t flank_bp, std::vector<KeepSpan>& keep);

// The keep spec of one walk under a site's scope (defined below).
std::vector<KeepSpan> walk_keep(const Graph& graph, const BubbleWalk& walk, const SiteScope& scope,
                                std::size_t path_index);

// The bubble's source, sink and interior. Every path of a bubble asks for them and a worker handles one
// bubble at a time, so the set is built once per bubble and thread rather than once per path.
const std::unordered_set<std::string>& bubble_members(const Bubble& bubble) {
    thread_local std::size_t cached_id = 0;
    thread_local std::string cached_source, cached_sink;
    thread_local std::size_t cached_inside = 0;
    thread_local bool cached = false;
    thread_local std::unordered_set<std::string> members;
    if (!cached || cached_id != bubble.id || cached_source != bubble.source || cached_sink != bubble.sink ||
        cached_inside != bubble.inside.size()) {
        members.clear();
        members.insert(bubble.source);
        members.insert(bubble.sink);
        members.insert(bubble.inside.begin(), bubble.inside.end());
        cached = true;
        cached_id = bubble.id;
        cached_source = bubble.source;
        cached_sink = bubble.sink;
        cached_inside = bubble.inside.size();
    }
    return members;
}

// Whether this bubble's crossings are read with bubble_crossings' fallback: only when no path crosses
// it through the boundary sides its table gives (bubble_orientation_holds), as inspect and call decide.
// Every path of a bubble asks, and a worker handles one bubble at a time, so it is kept per thread.
bool crossing_fallback(const Graph& graph, const std::vector<BubblePathIndex>& path_indexes, const Bubble& bubble) {
    thread_local const Graph* cached_graph = nullptr;
    thread_local std::size_t cached_id = 0;
    thread_local std::string cached_source, cached_sink;
    thread_local bool cached_value = false;
    if (cached_graph != &graph || cached_id != bubble.id || cached_source != bubble.source || cached_sink != bubble.sink) {
        cached_value = !bubble_orientation_holds(graph.paths, path_indexes, bubble);
        cached_graph = &graph;
        cached_id = bubble.id;
        cached_source = bubble.source;
        cached_sink = bubble.sink;
    }
    return cached_value;
}

BubbleWalk canonical_steps_for_bubble_path(
    const Graph& graph,
    const std::vector<BubblePathIndex>& path_indexes,
    const Bubble& bubble,
    std::size_t path_index) {

    // Every crossing on its own, read source to sink, as inspect reads them; nothing between two
    // crossings belongs to the site. A crossing with no interior node, a deletion, counts like any other.
    const PathRecord& path = graph.paths[path_index];
    BubbleWalk walk;
    const std::unordered_set<std::string>& members = bubble_members(bubble);
    for (const BubblePathInterval& iv :
         bubble_crossings(path, path_indexes[path_index], bubble, crossing_fallback(graph, path_indexes, bubble))) {
        std::vector<std::size_t> pos;
        const std::vector<PathStep> steps = canonical_bubble_path_steps(path, bubble, iv, &pos);
        if (steps.empty()) continue;
        // Only the bubble's nodes: a site table can leave a crossing's node out of its interior, and
        // that node is then not a feature of this bubble.
        bool left = true;
        bool first = true;
        for (std::size_t j = 0; j < steps.size(); ++j) {
            const PathStep& step = steps[j];
            if (!members.count(step.node_id)) {
                left = true;
                continue;
            }
            // An inversion site's step across a site inside it belongs to that site.
            if (!first && bubble.inversion_parent && inner_site_step(bubble, walk.steps.back().node_id, step.node_id))
                left = true;
            walk.run_start.push_back(left && !walk.steps.empty() ? 1 : 0);
            walk.crossing_start.push_back(first ? 1 : 0);
            walk.path_pos.push_back(pos[j]);
            walk.forward.push_back(iv.source_to_sink ? 1 : 0);
            walk.steps.push_back(step);
            left = false;
            first = false;
        }
    }
    return walk;
}

// Node ids and oriented-edge keys share one pool, so a node id containing the edge separator would be
// summed with an edge; it is refused.
void reject_node_edge_key_collision(const std::string& node_id) {
    if (node_id.find('>') != std::string::npos)
        throw std::runtime_error("describe: node id '" + node_id + "' contains '>', which is the "
                                 "oriented-edge key separator; node and edge features share one "
                                 "namespace and would be pooled together");
}

// Node and edge dosage over a path's walk. With `keep`, a node counts only on kept steps, and an edge only
// when both its ends are kept.
GraphDosage count_graph_dosage(const BubbleWalk& walk,
                               const std::vector<KeepSpan>* keep = nullptr) {
    const std::vector<PathStep>& steps = walk.steps;
    // Counted on views of the walk's own ids, and keyed by string once per distinct node and edge: a long
    // walk visits the same few nodes many times.
    std::unordered_map<std::string_view, std::uint32_t> nodes;
    struct EdgeRef {
        std::string_view from, to;
        bool from_rev, to_rev;
        bool operator==(const EdgeRef& o) const {
            return from == o.from && to == o.to && from_rev == o.from_rev && to_rev == o.to_rev;
        }
    };
    struct EdgeHash {
        std::size_t operator()(const EdgeRef& e) const {
            const std::size_t h = std::hash<std::string_view>{}(e.from) * 1000003u ^ std::hash<std::string_view>{}(e.to);
            return h * 4u + (e.from_rev ? 2u : 0u) + (e.to_rev ? 1u : 0u);
        }
    };
    std::unordered_map<EdgeRef, std::uint32_t, EdgeHash> edges;
    // A node the flank only partly reaches is counted whole: a node dosage has no partial count.
    for (std::size_t i = 0; i < steps.size(); ++i) {
        if (keep != nullptr && (!(*keep)[i].any() || (*keep)[i].no_node)) continue;
        ++nodes[steps[i].node_id];
    }
    for (std::size_t i = 0; i + 1 < steps.size(); ++i) {
        if (walk.run_start[i + 1]) continue;
        if (keep != nullptr && (!(*keep)[i].any() || !(*keep)[i + 1].any())) continue;
        ++edges[EdgeRef{steps[i].node_id, steps[i + 1].node_id, steps[i].reverse, steps[i + 1].reverse}];
    }
    GraphDosage dosage;
    for (const auto& [id, c] : nodes) {
        reject_node_edge_key_collision(std::string(id));
        dosage.node_counts.emplace(std::string(id), c);
    }
    for (const auto& [e, c] : edges)
        dosage.edge_counts.emplace(std::string(e.from) + (e.from_rev ? '-' : '+') + ">" + std::string(e.to) +
                                       (e.to_rev ? '-' : '+'),
                                   c);
    return dosage;
}

void update_graph_stats(
    const std::map<std::string, std::uint32_t>& counts,
    std::map<std::string, GraphFeatureStat>& stats) {

    for (const auto& [key, count] : counts) {
        auto& s = stats[key];
        ++s.present_paths;
        s.total_count += count;
        s.min_nonzero_count = std::min(s.min_nonzero_count, count);
        s.max_count = std::max(s.max_count, count);
    }
}

// The informative features, in key order.
std::vector<std::string> select_discriminative_graph_features(
    const std::map<std::string, GraphFeatureStat>& stats,
    std::size_t path_count,
    std::size_t min_paths) {

    std::vector<std::string> features;
    features.reserve(stats.size());
    for (const auto& [key, s] : stats) {
        if (feature_passes_filter(s.present_paths, s.min_nonzero_count, s.max_count, path_count, min_paths)) {
            features.push_back(key);
        }
    }
    return features;
}

void accumulate_graph_stats(
    const std::vector<PathWork>& work,
    std::map<std::string, GraphFeatureStat>& node_stats,
    std::map<std::string, GraphFeatureStat>& edge_stats) {

    for (const PathWork& w : work) {
        update_graph_stats(w.dosage.node_counts, node_stats);
        update_graph_stats(w.dosage.edge_counts, edge_stats);
    }
}

void write_graph_feature_map(
    const std::string& path,
    const std::vector<std::string>& node_features,
    const std::vector<std::string>& edge_features,
    const std::map<std::string, GraphFeatureStat>& node_stats,
    const std::map<std::string, GraphFeatureStat>& edge_stats,
    std::size_t path_count) {

    GzipWriter out(path);
    out.write("feature_id\tfeature_name\tfeature_type\tlabel\tpaths_present\tmin_count\tmax_count\ttotal_count\n");

    auto write_block = [&](const std::vector<std::string>& features,
                           const std::map<std::string, GraphFeatureStat>& stats,
                           const char* type,
                           std::string (*namer)(std::size_t)) {
        for (std::size_t i = 0; i < features.size(); ++i) {
            const auto it = stats.find(features[i]);
            if (it == stats.end()) {
                continue;
            }
            const GraphFeatureStat& s = it->second;
            const std::uint32_t min_count =
                (path_count > 0 && s.present_paths == path_count) ? s.min_nonzero_count : 0;
            out.write(std::to_string(i + 1));
            out.write("\t");
            out.write(namer(i + 1));
            out.write("\t");
            out.write(type);
            out.write("\t");
            out.write(tsv_sanitize(features[i]));
            out.write("\t");
            out.write(std::to_string(s.present_paths));
            out.write("\t");
            out.write(std::to_string(min_count));
            out.write("\t");
            out.write(std::to_string(s.max_count));
            out.write("\t");
            out.write(std::to_string(s.total_count));
            out.write("\n");
        }
    };

    write_block(node_features, node_stats, "node", node_feature_name);
    write_block(edge_features, edge_stats, "edge", edge_feature_name);
    out.close();
}

void write_graph_matrix(
    const std::string& path,
    const Bubble& bubble,
    const std::vector<PathMeta>& paths,
    const std::vector<PathWork>& work,
    const std::vector<std::string>& node_features,
    const std::vector<std::string>& edge_features) {

    GzipWriter out(path);
    std::string header = "bubble_id\tsample\thaplotype\tpath_name";
    for (std::size_t i = 0; i < node_features.size(); ++i) {
        header.push_back('\t');
        header += node_feature_name(i + 1);
    }
    for (std::size_t i = 0; i < edge_features.size(); ++i) {
        header.push_back('\t');
        header += edge_feature_name(i + 1);
    }
    header.push_back('\n');
    out.write(header);

    for (std::size_t pi = 0; pi < paths.size(); ++pi) {
        const PathMeta& meta = paths[pi];
        const GraphDosage& dosage = work[pi].dosage;

        std::string line;
        line.reserve(128 + (node_features.size() + edge_features.size()) * 3);
        line += std::to_string(bubble.id);
        line.push_back('\t');
        line += tsv_sanitize(meta.sample);
        line.push_back('\t');
        line += tsv_sanitize(meta.haplotype);
        line.push_back('\t');
        line += tsv_sanitize(meta.path_name);
        for (const std::string& key : node_features) {
            line.push_back('\t');
            const auto it = dosage.node_counts.find(key);
            line += (it == dosage.node_counts.end()) ? "0" : std::to_string(it->second);
        }
        for (const std::string& key : edge_features) {
            line.push_back('\t');
            const auto it = dosage.edge_counts.find(key);
            line += (it == dosage.edge_counts.end()) ? "0" : std::to_string(it->second);
        }
        line.push_back('\n');
        out.write(line);
    }
    out.close();
}

std::vector<std::string> sorted_nodes(const std::unordered_set<std::string>& nodes) {
    std::vector<std::string> out(nodes.begin(), nodes.end());
    std::sort(out.begin(), out.end());
    return out;
}

std::string join_nodes(const std::vector<std::string>& nodes) {
    std::string out;
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        if (i > 0) {
            out.push_back(';');
        }
        out += tsv_sanitize(nodes[i]);
    }
    return out;
}

std::string oriented_key(const PathStep& step) { return step.node_id + (step.reverse ? '-' : '+'); }

// Per-step keep spec for --variant-nodes: a variant node whole, and the flank_bp bases nearest it, or
// nearest a junction that skips variant nodes, along the walk, at the end of each node that faces it.
std::vector<KeepSpan> variant_keep_mask(
    const Graph& graph,
    const BubbleWalk& walk,
    const std::unordered_set<std::string>& variant_nodes,
    std::size_t flank_bp,
    const VariantJunctions* junctions) {

    const std::vector<PathStep>& steps = walk.steps;
    const std::size_t n = steps.size();
    std::vector<char> is_var(n, 0);
    std::vector<KeepSpan> keep(n);
    std::vector<std::size_t> len(n, 0);
    for (std::size_t i = 0; i < n; ++i) {
        const auto it = graph.nodes.find(steps[i].node_id);
        len[i] = (it == graph.nodes.end()) ? 0 : it->second.sequence.size();
        is_var[i] = variant_nodes.count(steps[i].node_id) ? 1 : 0;
        keep[i].full = is_var[i] != 0;
    }
    if (flank_bp == 0) {
        return keep;
    }
    // junction[i]: the step from i to i + 1 skips variant nodes another walk takes there.
    std::vector<char> junction(n, 0);
    if (junctions != nullptr)
        for (std::size_t i = 0; i + 1 < n; ++i)
            junction[i] = !walk.run_start[i + 1] && !is_var[i] && !is_var[i + 1] &&
                          junctions->before.count(oriented_key(steps[i])) &&
                          junctions->after.count(oriented_key(steps[i + 1]));
    // Forward sweep: `gap` is the distance from the last variant node or junction to this node's start, so
    // its first flank_bp - gap bases are kept; the backward sweep does the same for node ends. SIZE_MAX
    // means none yet. A flank does not reach across a re-entry, whose bases are not on the walk.
    std::size_t gap = std::numeric_limits<std::size_t>::max();
    for (std::size_t i = 0; i < n; ++i) {
        if (walk.run_start[i]) gap = std::numeric_limits<std::size_t>::max();
        if (is_var[i]) { gap = 0; continue; }
        if (gap != std::numeric_limits<std::size_t>::max()) {
            if (gap < flank_bp) keep[i].prefix = std::min(len[i], flank_bp - gap);
            gap += len[i];
        }
        if (junction[i]) gap = 0;
    }
    gap = std::numeric_limits<std::size_t>::max();
    for (std::size_t i = n; i-- > 0;) {
        if (is_var[i]) { gap = 0; }
        else {
            if (junction[i]) gap = 0;
            if (gap != std::numeric_limits<std::size_t>::max()) {
                if (gap < flank_bp) keep[i].suffix = std::min(len[i], flank_bp - gap);
                gap += len[i];
            }
        }
        if (walk.run_start[i]) gap = std::numeric_limits<std::size_t>::max();
    }
    return keep;
}

// Keep each crossing's source and sink, and the flank of each facing the interior: a deletion of the whole
// interior leaves its carriers only the boundaries. Boundary nodes themselves are constant across
// crossers, and the filter drops them.
void keep_bubble_boundaries(const Graph& graph, const BubbleWalk& walk,
                            std::size_t flank_bp, std::vector<KeepSpan>& keep) {
    const std::vector<PathStep>& steps = walk.steps;
    if (steps.empty() || keep.size() != steps.size()) return;
    auto node_len = [&](std::size_t i) {
        const auto it = graph.nodes.find(steps[i].node_id);
        return it == graph.nodes.end() ? std::size_t{0} : it->second.sequence.size();
    };
    for (std::size_t i = 0; i < steps.size(); ++i) {
        if (walk.crossing_start[i]) {   // a source: its downstream end faces the interior
            keep[i].boundary = true;
            keep[i].suffix = std::max(keep[i].suffix, std::min(node_len(i), flank_bp));
        }
        if (i + 1 == steps.size() || walk.crossing_start[i + 1]) {   // a sink: its upstream end does
            keep[i].boundary = true;
            keep[i].prefix = std::max(keep[i].prefix, std::min(node_len(i), flank_bp));
        }
    }
}

// One walk's own keep spec under --variant-nodes: its variant nodes and junctions, their flank and its
// boundaries.
std::vector<KeepSpan> walk_variant_keep(const Graph& graph, const BubbleWalk& walk, const SiteScope& scope,
                                        const VariantJunctions* junctions) {
    std::vector<KeepSpan> keep = variant_keep_mask(graph, walk, *scope.variant_nodes, scope.flank_bp, junctions);
    keep_bubble_boundaries(graph, walk, scope.flank_bp, keep);
    return keep;
}

// Adds one walk's keep spec to the site's, per node in the node's own orientation.
void merge_site_keep(const BubbleWalk& walk, const std::vector<KeepSpan>& keep, SiteKeep& site) {
    for (std::size_t i = 0; i < walk.steps.size(); ++i) {
        const KeepSpan& k = keep[i];
        if (!k.any()) continue;
        NodeKeep& n = site[walk.steps[i].node_id];
        n.full = n.full || k.full;
        n.boundary = n.boundary || k.boundary;
        const bool rev = walk.steps[i].reverse;
        n.head = std::max(n.head, rev ? k.suffix : k.prefix);
        n.tail = std::max(n.tail, rev ? k.prefix : k.suffix);
    }
}

// Under --variant-nodes, every walk keeps the bases any walk keeps: a node beside a deletion is then
// read on its carriers too, and so is the junction they walk in its place.
SiteKeep build_site_keep(const Graph& graph, const std::vector<BubblePathIndex>& path_indexes,
                         const Bubble& bubble, const SiteScope& scope) {
    std::vector<BubbleWalk> walks;
    for (std::size_t i = 0; i < graph.paths.size(); ++i) {
        BubbleWalk walk = canonical_steps_for_bubble_path(graph, path_indexes, bubble, i);
        if (!walk.empty()) walks.push_back(std::move(walk));
    }
    VariantJunctions junctions;
    for (const BubbleWalk& walk : walks)
        for (std::size_t i = 0; i + 1 < walk.steps.size(); ++i) {
            if (walk.run_start[i + 1]) continue;
            const bool var_here = scope.variant_nodes->count(walk.steps[i].node_id) != 0;
            const bool var_next = scope.variant_nodes->count(walk.steps[i + 1].node_id) != 0;
            if (!var_here && var_next) junctions.before.insert(oriented_key(walk.steps[i]));
            if (var_here && !var_next) junctions.after.insert(oriented_key(walk.steps[i + 1]));
        }
    SiteKeep site_keep;
    for (const BubbleWalk& walk : walks)
        merge_site_keep(walk, walk_variant_keep(graph, walk, scope, &junctions), site_keep);
    return site_keep;
}

// Whether `path_index` visits position `pos` on one of the owner's crossings.
bool owner_visits(const SharedVisits& visits, std::size_t path_index, std::size_t pos) {
    if (path_index >= visits.size()) return false;
    const std::vector<std::size_t>& v = visits[path_index];
    return std::binary_search(v.begin(), v.end(), pos);
}

// The keep spec of one walk under `scope`; empty when every base and node counts.
std::vector<KeepSpan> walk_keep(const Graph& graph, const BubbleWalk& walk, const SiteScope& scope,
                                std::size_t path_index) {
    std::vector<KeepSpan> keep;
    if (!scope.masked()) return keep;
    if (scope.variant_nodes != nullptr && scope.site_keep != nullptr) {
        keep.assign(walk.steps.size(), KeepSpan{});
        for (std::size_t i = 0; i < walk.steps.size(); ++i) {
            const auto it = scope.site_keep->find(walk.steps[i].node_id);
            if (it == scope.site_keep->end()) continue;
            const NodeKeep& n = it->second;
            const bool rev = walk.steps[i].reverse;
            keep[i].full = n.full;
            keep[i].boundary = n.boundary;
            keep[i].prefix = rev ? n.tail : n.head;
            keep[i].suffix = rev ? n.head : n.tail;
        }
    } else if (scope.variant_nodes != nullptr) {
        keep = walk_variant_keep(graph, walk, scope, nullptr);
    } else {
        keep.assign(walk.steps.size(), KeepSpan{});
        for (KeepSpan& k : keep) k.full = true;
    }
    const std::size_t junction = scope.kmer_size > 0 ? scope.kmer_size - 1 : 0;
    for (std::size_t i = 0; i < walk.steps.size(); ++i) {
        const bool source = walk.crossing_start[i] != 0;
        const bool sink = i + 1 == walk.steps.size() || walk.crossing_start[i + 1] != 0;
        if (!source && !sink) continue;
        const std::string& id = walk.steps[i].node_id;
        KeepSpan& k = keep[i];
        const auto sn = scope.shared_nodes.find(id);
        if (sn != scope.shared_nodes.end() && owner_visits(*sn->second, path_index, walk.path_pos[i]))
            k.no_node = true;
        const auto sk = scope.shared_kmers.find(id);
        if (sk != scope.shared_kmers.end() && owner_visits(*sk->second, path_index, walk.path_pos[i])) {
            const auto it = graph.nodes.find(id);
            const std::size_t len = it == graph.nodes.end() ? 0 : it->second.sequence.size();
            // Only the end facing the interior: a source's last bases, a sink's first ones.
            const std::size_t inner = std::min(len, junction);
            // The bases kept next to the interior, the whole node when its two kept ends meet.
            const std::size_t had = k.full || k.prefix + k.suffix >= len ? len : (source ? k.suffix : k.prefix);
            k.full = false;
            k.boundary = true;
            if (source) { k.suffix = std::min(had, inner); k.prefix = 0; }
            else { k.prefix = std::min(had, inner); k.suffix = 0; }
        }
    }
    return keep;
}

// The kept sequence of one step, the rest replaced by 'N' of the same length: the k-mer scanner resets
// on 'N', and offsets and node provenance are unchanged.
std::string masked_step_sequence(const std::string& node_seq, bool reverse, const KeepSpan& k) {
    const std::string oriented = reverse ? reverse_complement(node_seq) : node_seq;
    if (k.full) return oriented;
    if (!k.any()) return std::string(oriented.size(), 'N');
    std::string out(oriented.size(), 'N');
    // prefix and suffix are in the walk's reading direction, as `oriented` is.
    const std::size_t pre = std::min(k.prefix, oriented.size());
    for (std::size_t i = 0; i < pre; ++i) out[i] = oriented[i];
    const std::size_t suf = std::min(k.suffix, oriented.size());
    for (std::size_t i = 0; i < suf; ++i) out[oriented.size() - 1 - i] = oriented[oriented.size() - 1 - i];
    return out;
}

// The bases beyond step `pos` of `path`, in the reading direction of a crossing (`forward`): the `need` bases
// before it or after it, with the node each came from.
std::string bases_beyond(const Graph& graph, const PathRecord& path, std::size_t pos, bool forward, bool before,
                         std::size_t need, std::vector<SequenceSegment>* segments) {
    std::vector<std::pair<std::string, std::string>> pieces;   // node id, its bases in reading order
    std::size_t have = 0;
    const bool up = forward != before;   // walking the path towards higher indices
    for (std::size_t j = pos; have < need;) {
        if (up ? j + 1 >= path.steps.size() : j == 0) break;
        j = up ? j + 1 : j - 1;
        const PathStep& step = path.steps[j];
        const auto it = graph.nodes.find(step.node_id);
        if (it == graph.nodes.end() || it->second.sequence.empty()) break;
        const bool rev = forward ? step.reverse : !step.reverse;
        pieces.emplace_back(step.node_id, rev ? reverse_complement(it->second.sequence) : it->second.sequence);
        have += pieces.back().second.size();
    }
    if (before) std::reverse(pieces.begin(), pieces.end());
    std::string out;
    for (const auto& pc : pieces) out += pc.second;
    const std::size_t drop = out.size() > need ? out.size() - need : 0;
    const std::size_t from = before ? drop : 0;   // the bases next to the step
    if (segments != nullptr) {
        std::size_t at = 0;
        for (const auto& pc : pieces) {
            const std::size_t lo = std::max(at, from), hi = std::min(at + pc.second.size(), from + out.size() - drop);
            if (lo < hi) segments->push_back(SequenceSegment{lo - from, hi - from, pc.first});
            at += pc.second.size();
        }
    }
    return out.substr(from, out.size() - drop);
}

// A path's walk across the site spelled, bases outside the scope masked to 'N'; with `segments`, the node
// each base came from. Not complete, and empty, when a node of the walk has no sequence.
BubblePathSequence spell_walk(const Graph& graph, const BubbleWalk& walk, const std::vector<KeepSpan>& keep,
                              const SiteScope& scope, std::size_t path_index, bool segments) {
    BubblePathSequence out;
    if (walk.empty()) return out;
    const std::vector<PathStep>& steps = walk.steps;
    std::size_t total_len = 0;
    std::vector<const std::string*> seqs(steps.size());
    for (std::size_t si = 0; si < steps.size(); ++si) {
        const auto it = graph.nodes.find(steps[si].node_id);
        if (it == graph.nodes.end() || it->second.sequence.empty()) return out;
        seqs[si] = &it->second.sequence;
        total_len += it->second.sequence.size();
    }
    out.sequence.reserve(total_len + steps.size());
    if (segments) out.segments.reserve(steps.size());
    const auto add = [&](const std::string& piece, const std::string& node) {
        const std::size_t start = out.sequence.size();
        out.sequence += piece;
        if (segments) out.segments.push_back(SequenceSegment{start, out.sequence.size(), node});
    };
    const auto add_beyond = [&](std::size_t si, bool before, std::size_t need) {
        std::vector<SequenceSegment> ext;
        const std::size_t start = out.sequence.size();
        out.sequence += bases_beyond(graph, graph.paths[path_index], walk.path_pos[si], walk.forward[si] != 0,
                                     before, need, segments ? &ext : nullptr);
        for (const SequenceSegment& e : ext) out.segments.push_back(SequenceSegment{start + e.start, start + e.end, e.node_id});
    };
    for (std::size_t si = 0; si < steps.size(); ++si) {
        const PathStep& step = steps[si];
        const std::string& node_seq = *seqs[si];
        // A re-entry is spelled as one 'N', which resets the k-mer scanner: no k-mer spans it.
        if (walk.run_start[si]) out.sequence.push_back('N');
        const std::string piece = !keep.empty() ? masked_step_sequence(node_seq, step.reverse, keep[si])
                                  : step.reverse ? reverse_complement(node_seq) : node_seq;
        const bool source = walk.crossing_start[si] != 0;
        const bool sink = si + 1 == steps.size() || walk.crossing_start[si + 1] != 0;
        const bool extend = (source || sink) && scope.extend_kmers.count(step.node_id) &&
                            piece.size() + 1 < scope.kmer_size && piece.find('N') == std::string::npos;
        const std::size_t need = extend ? scope.kmer_size - 1 - piece.size() : 0;
        if (extend && source) add_beyond(si, true, need);
        add(piece, step.node_id);
        if (extend && sink) add_beyond(si, false, need);
    }
    out.complete = true;
    return out;
}

BubblePathSequence path_sequence_with_segments_for_bubble(
    const Graph& graph,
    const Bubble& bubble,
    const std::vector<BubblePathIndex>& path_indexes,
    std::size_t path_index,
    const SiteScope& scope) {
    const BubbleWalk walk = canonical_steps_for_bubble_path(graph, path_indexes, bubble, path_index);
    return spell_walk(graph, walk, walk_keep(graph, walk, scope, path_index), scope, path_index, true);
}

// The paths crossing the site, each walked once: its k-mer counts (into `stats` too) and its node and edge
// counts, in `work`, one entry per path returned.
std::vector<PathMeta> collect_path_metadata_and_stats(
    const Graph& graph,
    const Bubble& bubble,
    const std::vector<BubblePathIndex>& path_indexes,
    const DescribeOptions& options,
    std::unordered_map<std::uint64_t, KmerStats>& stats,
    const SiteScope& scope,
    std::vector<PathWork>& work) {

    std::vector<PathMeta> paths;
    paths.reserve(graph.paths.size());
    work.clear();
    std::unordered_map<std::uint64_t, std::uint32_t> counts;

    for (std::size_t i = 0; i < graph.paths.size(); ++i) {
        const BubbleWalk walk = canonical_steps_for_bubble_path(graph, path_indexes, bubble, i);
        if (walk.empty()) continue;
        const std::vector<KeepSpan> keep = walk_keep(graph, walk, scope, i);
        PathWork w;
        std::size_t length = 0;
        if (options.emit_kmers) {
            const BubblePathSequence spelled = spell_walk(graph, walk, keep, scope, i, false);
            if (!spelled.complete) continue;
            length = spelled.sequence.size();
            count_selected_features(spelled.sequence, options, counts);
            update_kmer_stats(counts, stats);
            w.kmers.assign(counts.begin(), counts.end());
            std::sort(w.kmers.begin(), w.kmers.end());
        } else {
            // Not spelled: the walk is complete when every node has a sequence.
            bool complete = true;
            for (std::size_t si = 0; si < walk.steps.size() && complete; ++si) {
                const auto it = graph.nodes.find(walk.steps[si].node_id);
                complete = it != graph.nodes.end() && !it->second.sequence.empty();
                if (complete) length += it->second.sequence.size() + (walk.run_start[si] ? 1 : 0);
            }
            if (!complete) continue;
        }
        if (options.emit_graph) w.dosage = count_graph_dosage(walk, scope.masked() ? &keep : nullptr);

        const auto [sample, haplotype] = parse_sample_haplotype(graph.paths[i].name);
        PathMeta meta;
        meta.path_index = i;
        meta.path_name = graph.paths[i].name;
        meta.sample = sample;
        meta.haplotype = haplotype;
        meta.length_bp = length;
        paths.push_back(std::move(meta));
        work.push_back(std::move(w));
    }

    return paths;
}

void write_feature_map(
    const std::string& path,
    const std::vector<std::uint64_t>& features,
    const std::unordered_map<std::uint64_t, KmerStats>& stats,
    const std::vector<std::unordered_set<std::string>>& feature_nodes,
    std::size_t k,
    std::size_t path_count) {

    GzipWriter out(path);
    out.write("feature_id\tfeature_name\tencoded_kmer\tkmer\tpaths_present\tmin_count\tmax_count\ttotal_count\tnode_count\tnodes\n");
    for (std::size_t i = 0; i < features.size(); ++i) {
        const std::uint64_t code = features[i];
        const auto it = stats.find(code);
        if (it == stats.end()) {
            continue;
        }
        const KmerStats& s = it->second;
        const std::uint32_t min_count = (s.present_paths == path_count) ? s.min_nonzero_count : 0;
        out.write(std::to_string(i + 1));
        out.write("\t");
        out.write(feature_name(i + 1));
        out.write("\t");
        out.write(std::to_string(code));
        out.write("\t");
        out.write(decode_kmer(code, k));
        out.write("\t");
        out.write(std::to_string(s.present_paths));
        out.write("\t");
        out.write(std::to_string(min_count));
        out.write("\t");
        out.write(std::to_string(s.max_count));
        out.write("\t");
        out.write(std::to_string(s.total_count));
        const std::size_t feature_id = i + 1;
        const std::vector<std::string> nodes =
            feature_id < feature_nodes.size() ? sorted_nodes(feature_nodes[feature_id]) : std::vector<std::string>{};
        out.write("\t");
        out.write(std::to_string(nodes.size()));
        out.write("\t");
        out.write(join_nodes(nodes));
        out.write("\n");
    }
    out.close();
}

std::vector<SparseKmerEntry> filtered_feature_counts_with_nodes(
    const std::unordered_map<std::uint64_t, KmerPathCount>& counts,
    const std::unordered_map<std::uint64_t, std::size_t>& feature_id_by_code) {

    std::vector<SparseKmerEntry> out;
    out.reserve(std::min(counts.size(), feature_id_by_code.size()));
    for (const auto& [code, path_count] : counts) {
        const auto id_it = feature_id_by_code.find(code);
        if (id_it == feature_id_by_code.end()) {
            continue;
        }
        SparseKmerEntry entry;
        entry.feature_id = id_it->second;
        entry.count = path_count.count;
        entry.nodes = sorted_nodes(path_count.nodes);
        out.push_back(std::move(entry));
    }
    std::sort(out.begin(), out.end(), [](const SparseKmerEntry& a, const SparseKmerEntry& b) {
        return a.feature_id < b.feature_id;
    });
    return out;
}

void write_sparse_jsonl(
    const std::string& path,
    const Graph& graph,
    const Bubble& bubble,
    const std::vector<BubblePathIndex>& path_indexes,
    const std::vector<PathMeta>& paths,
    const std::unordered_map<std::uint64_t, std::size_t>& feature_id_by_code,
    const DescribeOptions& options,
    std::vector<std::unordered_set<std::string>>& feature_nodes,
    const SiteScope& scope) {

    GzipWriter out(path);
    std::unordered_map<std::uint64_t, KmerPathCount> counts;

    for (const PathMeta& meta : paths) {
        const BubblePathSequence path_sequence = path_sequence_with_segments_for_bubble(
            graph,
            bubble,
            path_indexes,
            meta.path_index,
            scope);
        if (!path_sequence.complete) {
            continue;
        }
        count_selected_features_with_nodes(path_sequence, options, counts);
        const auto sparse = filtered_feature_counts_with_nodes(counts, feature_id_by_code);

        out.write("{\"bubble_id\":");
        out.write(std::to_string(bubble.id));
        out.write(",\"sample\":\"");
        out.write(json_escape(meta.sample));
        out.write("\",\"haplotype\":\"");
        out.write(json_escape(meta.haplotype));
        out.write("\",\"path_name\":\"");
        out.write(json_escape(meta.path_name));
        out.write("\",\"path_length_bp\":");
        out.write(std::to_string(meta.length_bp));
        out.write(",\"kmers\":[");
        for (std::size_t i = 0; i < sparse.size(); ++i) {
            if (i > 0) {
                out.write(",");
            }
            if (sparse[i].feature_id < feature_nodes.size()) {
                for (const std::string& node_id : sparse[i].nodes) {
                    feature_nodes[sparse[i].feature_id].insert(node_id);
                }
            }
            out.write("[");
            out.write(std::to_string(sparse[i].feature_id));
            out.write(",");
            out.write(std::to_string(sparse[i].count));
            out.write("]");
        }
        out.write("]}\n");
    }
    out.close();
}

void write_wide_matrix(
    const std::string& path,
    const Bubble& bubble,
    const std::vector<PathMeta>& paths,
    const std::vector<PathWork>& work,
    const std::vector<std::uint64_t>& features) {

    GzipWriter out(path);
    std::string header = "bubble_id\tsample\thaplotype\tpath_name";
    for (std::size_t i = 0; i < features.size(); ++i) {
        header.push_back('\t');
        header += feature_name(i + 1);
    }
    header.push_back('\n');
    out.write(header);

    for (std::size_t pi = 0; pi < paths.size(); ++pi) {
        const PathMeta& meta = paths[pi];
        std::string line;
        line.reserve(128 + features.size() * 3);
        line += std::to_string(bubble.id);
        line.push_back('\t');
        line += tsv_sanitize(meta.sample);
        line.push_back('\t');
        line += tsv_sanitize(meta.haplotype);
        line.push_back('\t');
        line += tsv_sanitize(meta.path_name);
        for (const std::uint64_t code : features) {
            line.push_back('\t');
            line += std::to_string(kmer_count(work[pi], code));
        }
        line.push_back('\n');
        out.write(line);
    }
    out.close();
}

// A pooled feature's counts by column -- a path's index in the graph, or a sample's in the --samples
// order -- sorted by column. Indices rather than names: a pool holds one entry per carrier cell.
using Carriers = std::vector<std::pair<std::uint32_t, std::uint32_t>>;

// Add `src` into `dst`, both sorted by column.
void add_carriers(Carriers& dst, Carriers&& src) {
    if (dst.empty()) { dst = std::move(src); return; }
    Carriers out;
    out.reserve(dst.size() + src.size());
    std::size_t i = 0, j = 0;
    while (i < dst.size() || j < src.size()) {
        if (j == src.size() || (i < dst.size() && dst[i].first < src[j].first)) out.push_back(dst[i++]);
        else if (i == dst.size() || src[j].first < dst[i].first) out.push_back(src[j++]);
        else {
            const std::uint64_t sum = std::uint64_t{dst[i].second} + src[j].second;
            out.emplace_back(dst[i].first, static_cast<std::uint32_t>(
                std::min<std::uint64_t>(sum, std::numeric_limits<std::uint32_t>::max())));
            ++i; ++j;
        }
    }
    dst.swap(out);
}

// K-mer cohort pool: kmer code -> its carriers, over the features passed in, summed across
// bubbles (the same canonical k-mer in two bubbles sums). Feeds the BIMBAM export.
using KmerPool = std::unordered_map<std::uint64_t, Carriers>;

void accumulate_kmer_counts(
    const std::vector<PathMeta>& paths,
    const std::vector<PathWork>& work,
    const std::vector<std::uint64_t>& features,
    KmerPool& pool) {

    const std::unordered_set<std::uint64_t> keep(features.begin(), features.end());
    if (keep.empty()) {
        return;
    }
    for (std::size_t pi = 0; pi < paths.size(); ++pi) {
        for (const auto& [code, c] : work[pi].kmers) {
            if (c > 0 && keep.count(code) != 0) {
                pool[code].emplace_back(static_cast<std::uint32_t>(paths[pi].path_index), c);
            }
        }
    }
}

// Node and edge pool, keyed by node id or "from>to" edge key: the graph analogue of KmerPool.
using GraphPool = std::map<std::string, Carriers>;

void accumulate_graph_counts(
    const std::vector<PathMeta>& paths,
    const std::vector<PathWork>& work,
    const std::vector<std::string>& node_features,
    const std::vector<std::string>& edge_features,
    GraphPool& pool) {

    const std::unordered_set<std::string> keep_nodes(node_features.begin(), node_features.end());
    const std::unordered_set<std::string> keep_edges(edge_features.begin(), edge_features.end());
    if (keep_nodes.empty() && keep_edges.empty()) {
        return;
    }
    for (std::size_t pi = 0; pi < paths.size(); ++pi) {
        const GraphDosage& dosage = work[pi].dosage;
        const auto col = static_cast<std::uint32_t>(paths[pi].path_index);
        for (const auto& [key, c] : dosage.node_counts) {
            if (c > 0 && keep_nodes.count(key) != 0) pool[key].emplace_back(col, c);
        }
        for (const auto& [key, c] : dosage.edge_counts) {
            if (c > 0 && keep_edges.count(key) != 0) pool[key].emplace_back(col, c);
        }
    }
}


std::string format_dosage(double d) {
    if (d == std::floor(d) && std::fabs(d) < 1e15) {
        return std::to_string(static_cast<long long>(d));
    }
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.6g", d);
    return buf;
}

// --scale-dosage: a feature's observed dosages mapped linearly onto 0..2. A feature with one dosage
// throughout has no range to map; it reads 0, so the scaled matrix holds 0..2 only.
double scaled_dosage(double v, double lo, double hi) {
    return hi > lo ? (v - lo) / (hi - lo) * 2.0 : 0.0;
}

// Reference coordinates for the sidecars: each node the reference walks spans [first, last] base (1-based, on
// the interval in the reference's name), and each site its boundaries' span. A row lies over its nodes on the
// reference; a node only other haplotypes walk takes its site's span, which is where its variation is.
struct RefSpans {
    bool known = false;
    std::string chrom;
    std::unordered_map<std::string, std::pair<std::size_t, std::size_t>> node;
    std::unordered_map<std::size_t, std::pair<std::size_t, std::size_t>> site;
};

// "\tchrom\tpos\tend" of a row over `nodes` (ids separated by ';' or, in an edge's key, '>', each perhaps ending
// in its orientation) kept at `bubbles`; "\t.\t.\t." when it has no place on the reference.
std::string span_columns(const RefSpans& spans, const std::string& nodes, const std::vector<std::size_t>& bubbles) {
    if (!spans.known) return "\t.\t.\t.";
    std::size_t lo = SIZE_MAX, hi = 0;
    bool off_reference = false;
    const bool edge = nodes.find('>') != std::string::npos;
    std::size_t start = 0;
    for (std::size_t i = 0; i <= nodes.size(); ++i) {
        if (i < nodes.size() && nodes[i] != ';' && nodes[i] != '>') continue;
        std::string id = nodes.substr(start, i - start);
        start = i + 1;
        // an edge's key gives each node its orientation; a plain id keeps whatever it ends in
        if (edge && !id.empty() && (id.back() == '+' || id.back() == '-')) id.pop_back();
        if (id.empty() || id == ".") continue;
        const auto it = spans.node.find(id);
        if (it == spans.node.end()) { off_reference = true; continue; }
        lo = std::min(lo, it->second.first); hi = std::max(hi, it->second.second);
    }
    if (off_reference || lo == SIZE_MAX)
        for (std::size_t b : bubbles)
            if (const auto it = spans.site.find(b); it != spans.site.end()) {
                lo = std::min(lo, it->second.first); hi = std::max(hi, it->second.second);
            }
    if (lo == SIZE_MAX) return "\t.\t.\t.";
    return "\t" + spans.chrom + "\t" + std::to_string(lo) + "\t" + std::to_string(hi);
}

// One BIMBAM row: a carrier holds its dosage, a column that crosses the feature's sites without carrying
// it reads 0, and one that does not cross them reads NA.
struct BimbamRow {
    std::string id;                                   // k-mer sequence, or node/edge id
    std::string nodes;                                // the nodes it lies on
    std::vector<std::size_t> bubbles;                 // the sites it was kept at
    const Carriers* carriers = nullptr;               // column -> dosage
    std::string encoding;                             // syncmer/all, or node/edge
};

// Rows by feature id, so the output does not depend on hash order.
inline void sort_rows(std::vector<BimbamRow>& rows) {
    std::sort(rows.begin(), rows.end(),
              [](const BimbamRow& a, const BimbamRow& b) { return a.id < b.id; });
}

// <out>/{haplotype|sample}/{kmers|graph|variant}/: one associate input each (matrix, sidecar, columns).
std::filesystem::path substrate_dir(const std::filesystem::path& out_dir, const char* agg,
                                    const char* sub) {
    const std::filesystem::path d = out_dir / agg / sub;
    std::filesystem::create_directories(d);
    return d;
}

void write_bimbam_rows(
    const std::string& bimbam_path,
    GzipWriter& annot,
    const std::string& layer,
    const RefSpans& spans,
    const std::vector<BimbamRow>& rows,
    const std::vector<std::string>& sample_order,
    const std::unordered_map<std::size_t, std::vector<char>>& bubble_traversers,
    bool scale_dosage,
    const std::size_t* min_paths = nullptr,
    std::unordered_set<std::string>* dropped = nullptr,
    const std::unordered_set<std::string>* skip = nullptr) {

    GzipWriter geno(bimbam_path);
    for (const BimbamRow& row : rows) {
        if (skip != nullptr && skip->count(row.id)) continue;
        std::vector<double> vals(sample_order.size(), std::numeric_limits<double>::quiet_NaN());
        std::vector<double> carried(sample_order.size(), 0.0);
        std::vector<char> carries(sample_order.size(), 0);
        for (const auto& [col, c] : *row.carriers) { carried[col] = c; carries[col] = 1; }
        double lo = std::numeric_limits<double>::infinity(), hi = -std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i < sample_order.size(); ++i) {
            // A pooled dosage is a sum over the feature's sites, so it is observed only where all of them
            // are; a partial sum would read as a low dosage.
            bool all_observed = true;
            for (std::size_t b : row.bubbles) {
                const auto bt = bubble_traversers.find(b);
                if (bt == bubble_traversers.end() || !bt->second[i]) {
                    all_observed = false;
                    break;
                }
            }
            if (all_observed) vals[i] = carried[i];
            if (std::isfinite(vals[i])) { lo = std::min(lo, vals[i]); hi = std::max(hi, vals[i]); }
        }
        // A row pooled over several sites was informative at each, but its sum need not be (an allele some
        // carry at one site and the others at the next sums to 1 for all), so the filter is applied again
        // over the columns observed at all its sites. Summed over samples it can be constant too.
        if (min_paths != nullptr) {
            std::size_t observed = 0, present = 0;
            double min_nz = std::numeric_limits<double>::infinity(), max_v = 0.0;
            for (double v : vals) {
                if (!std::isfinite(v)) continue;
                ++observed;
                if (v > 0) { ++present; min_nz = std::min(min_nz, v); max_v = std::max(max_v, v); }
            }
            const auto as_count = [](double v) {
                return static_cast<std::uint32_t>(std::min(v, static_cast<double>(std::numeric_limits<std::uint32_t>::max())));
            };
            if (!feature_passes_filter(present, present ? as_count(min_nz) : 0, as_count(max_v), observed, *min_paths)) {
                if (dropped != nullptr) dropped->insert(row.id);
                continue;
            }
        }
        // annotation sidecar: feature_id, layer, encoding, bubbles, nodes, chrom, pos, end
        std::string a = tsv_sanitize(row.id);
        a += '\t'; a += layer;
        a += '\t'; a += row.encoding.empty() ? std::string("count") : row.encoding;
        a += '\t';
        for (std::size_t i = 0; i < row.bubbles.size(); ++i) { if (i) a += ';'; a += std::to_string(row.bubbles[i]); }
        a += '\t'; a += row.nodes; a += span_columns(spans, row.nodes, row.bubbles); a += '\n';
        annot.write(a);
        // BIMBAM mean-genotype line: id, allele1, allele2, dosages... (--scale-dosage maps to 0..2)
        std::string line = tsv_sanitize(row.id);
        line += ", A, B";
        for (double v : vals) {
            line += ", ";
            if (!std::isfinite(v)) line += "NA";
            else line += format_dosage(scale_dosage ? scaled_dosage(v, lo, hi) : v);
        }
        line += '\n';
        geno.write(line);
    }
    geno.close();
}

BubbleDescribeResult describe_one_bubble(
    const DescribeOptions& options,
    const Graph& graph,
    const std::vector<BubblePathIndex>& path_indexes,
    const Bubble& bubble,
    const SiteScope& site_scope,
    KmerPool* kmer_pool = nullptr,
    GraphPool* graph_pool = nullptr) {

    // A site sharing a boundary for k-mers gets its keep spec from the caller, merged with its neighbour's.
    SiteKeep site_keep;
    SiteScope scope = site_scope;
    if (scope.variant_nodes != nullptr && scope.site_keep == nullptr) {
        site_keep = build_site_keep(graph, path_indexes, bubble, scope);
        scope.site_keep = &site_keep;
    }

    BubbleDescribeResult result;
    const std::filesystem::path out_dir(options.out_dir);
    const std::filesystem::path bubble_dir = out_dir / ("bubble_" + std::to_string(bubble.id));
    result.feature_map_path = (bubble_dir / "kmer_features.tsv.gz").string();
    result.matrix_path = (bubble_dir / "kmer_matrix.tsv.gz").string();
    result.counts_jsonl_path = (bubble_dir / "kmer_counts.jsonl.gz").string();
    result.graph_feature_map_path = (bubble_dir / "graph_features.tsv.gz").string();
    result.graph_matrix_path = (bubble_dir / "graph_matrix.tsv.gz").string();

    std::unordered_map<std::uint64_t, KmerStats> stats;
    stats.reserve(4096);
    std::vector<PathWork> work;
    const std::vector<PathMeta> paths = collect_path_metadata_and_stats(
        graph,
        bubble,
        path_indexes,
        options,
        stats,
        scope,
        work);

    result.paths = paths.size();
    if (paths.empty()) {
        result.status = "skipped:no-paths";
        result.feature_map_path = ".";
        result.matrix_path = ".";
        result.counts_jsonl_path = ".";
        result.matrix_reason = "skipped:no-paths";
        result.graph_feature_map_path = ".";
        result.graph_matrix_path = ".";
        return result;
    }
    std::filesystem::create_directories(bubble_dir);
    result.traversing_paths.reserve(paths.size());
    for (const PathMeta& m : paths) result.traversing_paths.push_back(static_cast<std::uint32_t>(m.path_index));

    // ---- k-mer substrate (skipped under --only-graph) ----
    if (options.emit_kmers) {
        const std::vector<std::uint64_t> features =
            select_discriminative_features(stats, paths.size(), options.min_feature_paths);
        result.features = features.size();
        result.features_total = stats.size();
        const auto feature_id_by_code = make_feature_id_map(features);

        std::vector<std::unordered_set<std::string>> feature_nodes(features.size() + 1);
        write_sparse_jsonl(
            result.counts_jsonl_path,
            graph,
            bubble,
            path_indexes,
            paths,
            feature_id_by_code,
            options,
            feature_nodes,
            scope);
        write_feature_map(result.feature_map_path, features, stats, feature_nodes, options.kmer_size, paths.size());

        // Each kept k-mer's nodes, for the pooled sidecar.
        if (kmer_pool != nullptr) {
            for (const std::uint64_t code : features) {
                const auto it = feature_id_by_code.find(code);
                if (it != feature_id_by_code.end() && it->second < feature_nodes.size())
                    result.kmer_nodes.emplace(code, sorted_nodes(feature_nodes[it->second]));
            }
        }

        if (kmer_pool != nullptr) {
            // The k-mers inside a shared boundary are pooled whatever the filter says here.
            std::vector<std::uint64_t> pooled = features;
            std::map<std::uint64_t, std::vector<std::string>> boundary_nodes;
            std::unordered_map<std::uint64_t, std::uint32_t> inside;
            for (const std::string& n : scope.pooled_kmer_nodes) {
                const auto nit = graph.nodes.find(n);
                if (nit == graph.nodes.end()) continue;
                count_canonical_kmers(nit->second.sequence, options.kmer_size, inside);
                for (const auto& [code, c] : inside)
                    if (stats.count(code) && !std::binary_search(features.begin(), features.end(), code))
                        boundary_nodes[code].push_back(n);
            }
            for (auto& [code, ns] : boundary_nodes) {
                pooled.push_back(code);
                result.kmer_nodes.emplace(code, sorted_nodes(std::unordered_set<std::string>(ns.begin(), ns.end())));
            }
            accumulate_kmer_counts(paths, work, pooled, *kmer_pool);
        }

        const bool wide_allowed_by_cap =
            options.max_wide_features == 0 || features.size() <= options.max_wide_features;
        if (options.write_wide_matrix && (options.force_wide_matrix || wide_allowed_by_cap)) {
            write_wide_matrix(result.matrix_path, bubble, paths, work, features);
            result.matrix_written = true;
            result.matrix_reason = "written";
        } else if (!options.write_wide_matrix) {
            result.matrix_path = ".";
            result.matrix_reason = "skipped:disabled";
        } else {
            result.matrix_path = ".";
            result.matrix_reason = "skipped:feature-cap";
        }
    } else {
        result.feature_map_path = ".";
        result.matrix_path = ".";
        result.counts_jsonl_path = ".";
        result.matrix_reason = "skipped:only-graph";
    }

    // ---- node + edge dosage substrate (skipped under --only-kmers) ----
    if (options.emit_graph) {
        std::map<std::string, GraphFeatureStat> node_stats;
        std::map<std::string, GraphFeatureStat> edge_stats;
        accumulate_graph_stats(work, node_stats, edge_stats);
        const std::vector<std::string> node_features =
            select_discriminative_graph_features(node_stats, paths.size(), options.min_feature_paths);
        const std::vector<std::string> edge_features =
            select_discriminative_graph_features(edge_stats, paths.size(), options.min_feature_paths);
        result.node_features = node_features.size();
        result.node_features_total = node_stats.size();
        result.edge_features = edge_features.size();
        result.edge_features_total = edge_stats.size();
        if (graph_pool != nullptr) {
            // A shared boundary's node row is pooled whatever the filter says here.
            std::vector<std::string> pooled = node_features;
            for (const std::string& n : scope.pooled_nodes)
                if (node_stats.count(n)) pooled.push_back(n);
            accumulate_graph_counts(paths, work, pooled, edge_features, *graph_pool);
        }
        write_graph_feature_map(
            result.graph_feature_map_path, node_features, edge_features, node_stats, edge_stats, paths.size());
        const bool graph_wide_allowed =
            options.max_wide_features == 0 ||
            node_features.size() + edge_features.size() <= options.max_wide_features;
        if (options.write_wide_matrix && (options.force_wide_matrix || graph_wide_allowed)) {
            write_graph_matrix(result.graph_matrix_path, bubble, paths, work, node_features, edge_features);
            result.graph_matrix_written = true;
        } else {
            result.graph_matrix_path = ".";
        }
    } else {
        result.graph_feature_map_path = ".";
        result.graph_matrix_path = ".";
    }

    return result;
}

// Split on a single delimiter (keeps empty fields).
std::vector<std::string> split_on(const std::string& s, char delim) {
    std::vector<std::string> out;
    std::size_t start = 0;
    for (std::size_t p = 0; p <= s.size(); ++p) {
        if (p == s.size() || s[p] == delim) { out.push_back(s.substr(start, p - start)); start = p + 1; }
    }
    return out;
}

// Value of an INFO key ("KEY=val;KEY2=val2" -> val; "" if absent or flag-only).
std::string info_value(const std::string& info, const std::string& key) {
    const std::string needle = key + "=";
    std::size_t p = 0;
    while (p < info.size()) {
        std::size_t e = info.find(';', p);
        if (e == std::string::npos) e = info.size();
        if (info.compare(p, needle.size(), needle) == 0)
            return info.substr(p + needle.size(), e - (p + needle.size()));
        p = e + 1;
    }
    return std::string();
}

// A GENES value as call writes it: names joined by ',', each percent-encoded (%3B %2C %3D %3A %25 %20
// %09 %0D %0A). Decoded and joined by ';', as associate's --node-genes column is; a decoded tab or line
// break would split the sidecar's row, so it reads as a space.
std::string vcf_unescape_list(const std::string& value) {
    if (value.empty() || value == ".") return value;
    std::string out;
    for (std::size_t i = 0; i < value.size(); ++i) {
        char c = value[i];
        if (c == ',') { out += ';'; continue; }
        if (c == '%' && i + 2 < value.size() && std::isxdigit(static_cast<unsigned char>(value[i + 1])) &&
            std::isxdigit(static_cast<unsigned char>(value[i + 2]))) {
            c = static_cast<char>(std::stoi(value.substr(i + 1, 2), nullptr, 16));
            i += 2;
        }
        out += (c == '\t' || c == '\r' || c == '\n') ? ' ' : c;
    }
    return out;
}

// One variant row read from call's VCF: CN on a copy-number record, an allele indicator per ALT on a
// multiallelic one, GT otherwise; GT '.' reads NA.
struct VariantBimbam {
    // af and an are the VCF's, over the panel's haplotypes; the sample level computes its own.
    std::string id, svtype, bubble, nodes, gene = ".", af = ".", an = ".";
    std::string chrom = ".", pos = ".", end = ".";   // the record's CHROM, POS and END
    std::unordered_map<std::string, double> dose;  // haplotype -> dosage; absent reads NA
    // Whether the haplotype carries this row's ALT, from GT: a copy number is above 0 on the reference too.
    std::unordered_map<std::string, char> is_alt;
};

void emit_variant_substrate(const DescribeOptions& options, DescribeSummary& summary) {
    std::ifstream vin(options.variant_vcf_path);
    if (!vin) throw std::runtime_error("describe --variant-vcf: cannot open " + options.variant_vcf_path);
    const std::filesystem::path out_dir(options.out_dir);
    std::filesystem::create_directories(out_dir);

    std::vector<std::string> hap_order;     // VCF sample columns (haplotype paths) = BIMBAM column order
    std::vector<VariantBimbam> rows;
    // A record's nodes come from variant_nodes.tsv, which lists every node its carriers walk (EVENT_NODES
    // of a copy-number record is only its boundary). The two files must come from one call run: every
    // record of the described sites has its row, with the same site and type, and every row its record.
    struct NodeRow { std::string bubble, svtype, nodes; bool matched = false; };
    std::unordered_map<std::string, NodeRow> record_nodes;
    if (!options.variant_nodes_path.empty()) {
        std::ifstream vn(options.variant_nodes_path);
        if (!vn) throw std::runtime_error("describe: failed to open variant nodes: " + options.variant_nodes_path);
        std::string l;
        bool head = true;
        while (std::getline(vn, l)) {
            if (!l.empty() && l.back() == '\r') l.pop_back();
            if (head) {
                // Without call's header the first record would be read as one, and its absence reported as
                // a mismatch between the two files.
                if (l.empty()) continue;
                const std::vector<std::string> h = split_on(l, '\t');
                if (h.size() < 4 || h[0] != "variant_id" || h[1] != "bubble_id" || h[2] != "svtype" || h[3] != "node_ids")
                    throw std::runtime_error("describe: " + options.variant_nodes_path + ": the header is not call's "
                                             "(variant_id, bubble_id, svtype, node_ids); this is not a "
                                             "variant_nodes.tsv written by `call`");
                head = false;
                continue;
            }
            const std::vector<std::string> f = split_on(l, '\t');
            if (f.size() >= 4 && !f[0].empty() && !f[3].empty()) record_nodes[f[0]] = NodeRow{f[1], f[2], f[3]};
        }
    }
    const auto mismatch = [&](const std::string& what) {
        return std::runtime_error("describe: --variant-vcf " + options.variant_vcf_path + " and --variant-nodes " +
                                  options.variant_nodes_path + " are not from one call run: " + what);
    };
    std::string line;
    bool saw_header = false;
    std::size_t lineno = 0;
    std::unordered_set<std::string> seen_ids;     // record IDs as they appear in the VCF
    std::unordered_set<std::string> emitted_ids;  // feature ids actually written, generated included
    while (std::getline(vin, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        if (line[0] == '#') {
            if (line.rfind("#CHROM", 0) == 0) {
                std::vector<std::string> f = split_on(line, '\t');
                for (std::size_t i = 9; i < f.size(); ++i) hap_order.push_back(f[i]);
                // A repeated column would make a dosage's haplotype ambiguous.
                std::unordered_set<std::string> seen_h;
                for (const std::string& h : hap_order)
                    if (!seen_h.insert(h).second)
                        throw std::runtime_error("describe --variant-vcf: duplicate sample column '" + h +
                                                 "' in " + options.variant_vcf_path);
                saw_header = true;
            }
            continue;
        }
        ++lineno;
        if (!saw_header)
            throw std::runtime_error("describe --variant-vcf: a data line appears before the #CHROM "
                                     "header in " + options.variant_vcf_path);
        std::vector<std::string> f = split_on(line, '\t');
        // A line of the wrong width would put dosages under the wrong haplotypes.
        if (f.size() != 9 + hap_order.size())
            throw std::runtime_error("describe --variant-vcf: " + options.variant_vcf_path + " line " +
                                     std::to_string(lineno) + " has " + std::to_string(f.size()) +
                                     " fields, the #CHROM header declares " +
                                     std::to_string(9 + hap_order.size()));
        const std::string& id = f[2];
        if (id == "." || id.empty())
            throw std::runtime_error("describe --variant-vcf: " + options.variant_vcf_path + " line " +
                                     std::to_string(lineno) + " has no record ID; each row is named by its "
                                     "record's ID");
        if (!seen_ids.insert(id).second)
            throw std::runtime_error("describe --variant-vcf: duplicate record ID '" + id + "' in " +
                                     options.variant_vcf_path);
        const std::string& info = f[7];
        const std::string svtype = info_value(info, "SVTYPE");
        const std::string bubble = info_value(info, "BUBBLE_ID");
        if (bubble.empty())
            throw std::runtime_error("describe --variant-vcf: record '" + id + "' carries no BUBBLE_ID");
        // --bubble-id restricts this substrate too.
        if (!options.bubble_ids.empty() &&
            std::find(options.bubble_ids.begin(), options.bubble_ids.end(),
                      static_cast<std::size_t>(std::strtoull(bubble.c_str(), nullptr, 10))) == options.bubble_ids.end())
            continue;
        std::string nodes = info_value(info, "EVENT_NODES");
        if (!options.variant_nodes_path.empty()) {
            const auto rn = record_nodes.find(id);
            if (rn == record_nodes.end()) throw mismatch("record '" + id + "' has no row in the node table");
            if (rn->second.bubble != bubble || rn->second.svtype != svtype)
                throw mismatch("record '" + id + "' is " + svtype + " at site " + bubble + " in the VCF but " +
                               rn->second.svtype + " at site " + rn->second.bubble + " in the node table");
            rn->second.matched = true;
            nodes = rn->second.nodes;
        }
        const std::string gene = vcf_unescape_list(info_value(info, "GENES"));
        // AF is Number=A, one value per ALT; a multiallelic record's rows each take their own.
        const std::string af = info_value(info, "AF");
        const std::vector<std::string> af_parts = split_on(af, ',');
        const std::string an = info_value(info, "AN");
        // The ALT column sets the rows and the GT indices they test: REF plus one allele per ALT.
        const std::size_t n_alt = (f[4] == "." || f[4].empty()) ? 0 : split_on(f[4], ',').size();
        const int nalleles = static_cast<int>(n_alt) + 1;
        // An optional NALLELES (REF plus the ALTs) must agree with it.
        const std::string nalleles_s = info_value(info, "NALLELES");
        if (!nalleles_s.empty()) {
            std::size_t used = 0;
            int declared = 0;
            try { declared = std::stoi(nalleles_s, &used); } catch (const std::exception&) { used = 0; }
            if (used != nalleles_s.size() || declared < 1)
                throw std::runtime_error("describe --variant-vcf: record '" + id + "' has NALLELES='" +
                                         nalleles_s + "', which is not a positive integer");
            if (declared != nalleles)
                throw std::runtime_error("describe --variant-vcf: record '" + id + "' declares NALLELES=" +
                                         nalleles_s + " but its ALT column carries " +
                                         std::to_string(n_alt) + " allele(s); the two must agree "
                                         "(NALLELES counts REF plus the ALTs)");
        }
        // A list of the wrong length would label an allele with another allele's frequency.
        if (!af.empty() && af != "." &&
            af_parts.size() != static_cast<std::size_t>(nalleles > 1 ? nalleles - 1 : 1))
            throw std::runtime_error("describe --variant-vcf: record '" + id + "' has AF with " +
                                     std::to_string(af_parts.size()) + " value(s), which is not one "
                                     "per ALT (Number=A)");
        const std::vector<std::string> fmt_keys = split_on(f[8], ':');
        int gt_i = -1, cn_i = -1;
        for (std::size_t i = 0; i < fmt_keys.size(); ++i) {
            if (fmt_keys[i] == "GT") gt_i = static_cast<int>(i);
            else if (fmt_keys[i] == "CN") cn_i = static_cast<int>(i);
        }
        // A copy-number record (CNV; DUP in older VCFs) is one row, its CN, however many ALTs it lists.
        // Any other record has one row per ALT.
        const bool is_dup = (svtype == "CNV" || svtype == "DUP");
        const bool is_multi = n_alt > 1 && !is_dup;
        const int n_rows = is_multi ? nalleles - 1 : 1;
        for (int a = 0; a < n_rows; ++a) {
            VariantBimbam v;
            v.id = is_multi ? (id + "_a" + std::to_string(a + 1)) : id;
            // Record ids and the per-ALT ids made from them share one namespace.
            if (!emitted_ids.insert(v.id).second)
                throw std::runtime_error("describe --variant-vcf: feature id '" + v.id +
                                         "' is emitted twice; per-ALT ids and record IDs share one "
                                         "namespace and must be unique across both");
            v.svtype = svtype.empty() ? "." : svtype;
            v.bubble = bubble; v.nodes = nodes.empty() ? "." : nodes;
            v.chrom = f[0]; v.pos = f[1];
            {   // END, else the last base REF spans
                const std::string e = info_value(info, "END");
                v.end = !e.empty() ? e : std::to_string(std::strtoull(f[1].c_str(), nullptr, 10) + (f[3].empty() ? 0 : f[3].size() - 1));
            }
            v.gene = gene.empty() ? "." : gene;
            if (af.empty()) v.af = ".";
            else if (!is_multi) v.af = af_parts.size() == 1 ? af : ".";
            else v.af = (static_cast<std::size_t>(a) < af_parts.size()) ? af_parts[static_cast<std::size_t>(a)] : ".";
            v.an = an.empty() ? "." : an;
            for (std::size_t s = 0; s < hap_order.size() && 9 + s < f.size(); ++s) {
                const std::vector<std::string> sub = split_on(f[9 + s], ':');
                if (gt_i < 0)
                    throw std::runtime_error("describe --variant-vcf: record '" + id +
                                             "' has no GT in FORMAT; the first field is not a genotype");
                const std::string gt = gt_i < static_cast<int>(sub.size()) ? sub[gt_i] : std::string(".");
                if (gt == "." || gt.empty()) continue;  // not crossing the site: NA
                // Each column is one haplotype, so a genotype is one allele index.
                if (gt.find('/') != std::string::npos || gt.find('|') != std::string::npos)
                    throw std::runtime_error("describe --variant-vcf: sample '" + hap_order[s] +
                                             "' has the diploid genotype '" + gt + "' at record '" + id +
                                             "'; describe reads call's per-haplotype VCF (call without "
                                             "--diploid), and --samples gives the per-sample level");
                int gti = 0;
                {
                    std::size_t used = 0;
                    try { gti = std::stoi(gt, &used); } catch (const std::exception&) { used = 0; }
                    if (used != gt.size() || gti < 0)
                        throw std::runtime_error("describe --variant-vcf: sample '" + hap_order[s] +
                                                 "' has the genotype '" + gt + "' at record '" + id +
                                                 "', which is not an allele index");
                }
                const int max_allele = nalleles > 1 ? nalleles - 1 : 1;
                if (gti > max_allele)
                    throw std::runtime_error("describe --variant-vcf: record '" + id + "' gives sample '" +
                                             hap_order[s] + "' allele index " + std::to_string(gti) +
                                             ", but the record has only " + std::to_string(max_allele) +
                                             " ALT allele(s)");
                double d;
                if (is_dup) {
                    // A copy number has no GT stand-in: presence is a different quantity.
                    if (cn_i < 0 || cn_i >= static_cast<int>(sub.size()) || sub[cn_i] == ".")
                        throw std::runtime_error("describe --variant-vcf: copy-number record '" + id +
                                                 "' has no usable CN for sample '" + hap_order[s] +
                                                 "'; a copy-number feature cannot fall back to GT presence");
                    std::size_t cused = 0;
                    double cn = 0.0;
                    try { cn = std::stod(sub[cn_i], &cused); } catch (const std::exception&) { cused = 0; }
                    if (cused != sub[cn_i].size() || !std::isfinite(cn) || cn < 0.0)
                        throw std::runtime_error("describe --variant-vcf: copy-number record '" + id +
                                                 "' has CN '" + sub[cn_i] + "' for sample '" + hap_order[s] +
                                                 "', which is not a finite non-negative number");
                    d = cn;
                } else if (is_multi) {
                    d = (gti == a + 1) ? 1.0 : 0.0;
                } else {
                    d = static_cast<double>(gti);
                }
                v.dose[hap_order[s]] = d;
                v.is_alt[hap_order[s]] = static_cast<char>((is_multi ? (gti == a + 1) : (gti >= 1)) ? 1 : 0);
            }
            rows.push_back(std::move(v));
        }
    }

    if (!saw_header)
        throw std::runtime_error("describe --variant-vcf: " + options.variant_vcf_path +
                                 " has no #CHROM header line");
    for (const auto& [rid, row] : record_nodes) {
        if (row.matched) continue;
        const std::size_t b = static_cast<std::size_t>(std::strtoull(row.bubble.c_str(), nullptr, 10));
        if (options.bubble_ids.empty() ||
            std::find(options.bubble_ids.begin(), options.bubble_ids.end(), b) != options.bubble_ids.end())
            throw mismatch("node-table row '" + rid + "' has no record in the VCF");
    }
    if (rows.empty())
        throw std::runtime_error("describe --variant-vcf: " + options.variant_vcf_path +
                                 " yielded no variant rows");

    // af and an are the VCF's unless given: the sample level passes its own, over its samples.
    auto write_annot = [&](GzipWriter& annot, const VariantBimbam& v,
                           const std::string* af_override = nullptr,
                           const std::string* an_override = nullptr) {
        std::string a = tsv_sanitize(v.id);
        a += "\tvariant\tdosage\t"; a += v.bubble; a += '\t'; a += v.nodes;
        a += '\t'; a += v.svtype; a += '\t'; a += v.gene;
        a += '\t'; a += (af_override ? *af_override : v.af);
        a += '\t'; a += (an_override ? *an_override : v.an);
        a += '\t'; a += v.chrom; a += '\t'; a += v.pos; a += '\t'; a += v.end;
        a += '\n';
        annot.write(a);
    };
    const std::string annot_header = "feature_id\tlayer\tencoding\tbubbles\tnodes\tsvtype\tgene\tAF\tAN\tchrom\tpos\tend\n";

    // Haplotype level, in haplotype/variant/.
    {
        const std::filesystem::path dir = substrate_dir(out_dir, "haplotype", "variant");
        GzipWriter sout((dir / "samples.txt.gz").string());
        for (const std::string& s : hap_order) { sout.write(s); sout.write("\n"); }
        sout.close();
        GzipWriter annot((dir / "feature_annot.variant.tsv.gz").string());
        annot.write(annot_header);
        GzipWriter geno((dir / "bimbam_variant.bimbam.gz").string());
        for (const VariantBimbam& v : rows) {
            write_annot(annot, v);
            double lo = std::numeric_limits<double>::infinity(), hi = -lo;
            for (const auto& kv : v.dose) { lo = std::min(lo, kv.second); hi = std::max(hi, kv.second); }
            std::string g = tsv_sanitize(v.id); g += ", A, B";
            for (const std::string& s : hap_order) {
                g += ", ";
                const auto it = v.dose.find(s);
                if (it == v.dose.end()) g += "NA";
                else g += format_dosage(options.scale_dosage ? scaled_dosage(it->second, lo, hi) : it->second);
            }
            g += '\n';
            geno.write(g);
        }
        geno.close(); annot.close();
        summary.files_written += 3;
    }

    // Sample level: the sum of a sample's haplotype dosages, observed only when each haplotype is.
    if (!options.samples_path.empty()) {
        // Every haplotype in the table must be a VCF column.
        std::unordered_set<std::string> vcf_haps(hap_order.begin(), hap_order.end());
        // The graph pass, when it runs, has noted the table's ploidy already.
        const CosigtTable cosigt = read_cosigt_table(options.samples_path, "describe", &vcf_haps,
                                                     "a haplotype column in --variant-vcf",
                                                     options.quiet || options.emit_kmers || options.emit_graph);
        const auto& sample_order = cosigt.sample_order;
        const std::filesystem::path dir = substrate_dir(out_dir, "sample", "variant");
        GzipWriter sout((dir / "samples.txt.gz").string());
        for (const std::string& s : sample_order) { sout.write(s); sout.write("\n"); }
        sout.close();
        GzipWriter annot((dir / "feature_annot.variant.tsv.gz").string());
        annot.write(annot_header);
        GzipWriter geno((dir / "bimbam_variant.bimbam.gz").string());
        for (const VariantBimbam& v : rows) {
            std::unordered_map<std::string, double> samp;
            std::unordered_set<std::string> covered;
            // AN and AC as the VCF defines them, over this matrix's samples: one allele per assigned
            // haplotype observed here, and those carrying this row's ALT.
            std::size_t an_s = 0, ac_s = 0;
            const auto observed = [&](const std::string& c) { return v.dose.count(c) != 0; };
            for (const auto& [sm, slots] : cosigt.slots) {
                std::vector<const std::string*> pick;
                for (const auto& slot : slots) {
                    const std::string* c = slot_contig(slot, observed);
                    if (c == nullptr) break;
                    pick.push_back(c);
                }
                if (pick.size() != slots.size()) continue;   // NA
                covered.insert(sm);
                for (const std::string* h : pick) {
                    samp[sm] += v.dose.at(*h);
                    const auto ai = v.is_alt.find(*h);
                    ++an_s;
                    if (ai != v.is_alt.end() && ai->second != 0) ++ac_s;
                }
            }
            const std::string an_str = std::to_string(an_s);
            const std::string af_str =
                an_s ? format_dosage(static_cast<double>(ac_s) / static_cast<double>(an_s)) : std::string(".");
            write_annot(annot, v, &af_str, &an_str);
            double lo = std::numeric_limits<double>::infinity(), hi = -lo;
            for (const std::string& s : covered) { lo = std::min(lo, samp[s]); hi = std::max(hi, samp[s]); }
            std::string g = tsv_sanitize(v.id); g += ", A, B";
            for (const std::string& s : sample_order) {
                g += ", ";
                if (!covered.count(s)) g += "NA";
                else g += format_dosage(options.scale_dosage ? scaled_dosage(samp[s], lo, hi) : samp[s]);
            }
            g += '\n';
            geno.write(g);
        }
        geno.close(); annot.close();
        summary.files_written += 3;
    }

    summary.variant_rows = rows.size();
}

} // namespace

void describe_kmers_from_graph(
    const DescribeOptions& options,
    DescribeSummary* summary_out) {

    if (options.gfa_path.empty()) {
        throw std::runtime_error("Describe requires gfa_path");
    }
    if (options.bubbles_csv_in.empty()) {
        throw std::runtime_error("Describe requires bubbles_csv_in");
    }
    if (options.kmer_size == 0 || options.kmer_size > 31) {
        throw std::runtime_error("Describe k-mer size must be in [1,31]");
    }
    if (options.feature_mode == DescribeFeatureMode::Syncmer) {
        const std::size_t s = effective_syncmer_s(options);
        if (s == 0 || s >= options.kmer_size) {
            throw std::runtime_error("Describe syncmer s-mer size must be > 0 and < k-mer size");
        }
    }

    // Inputs are validated before anything is written.
    ParseGfaOptions parse_options;
    parse_options.include_paths = true;
    parse_options.include_sequences = true;
    const Graph graph = parse_gfa(options.gfa_path, parse_options);
    require_graph_paths(graph.nodes.size(), graph.paths.size(), "describe");
    // Walks are spelled by concatenating whole segments: steps must name nodes and links, path names be
    // unique, and overlaps zero.
    validate_graph_paths(graph, "describe", true, true);

    std::vector<Bubble> all_bubbles = read_bubbles_csv(options.bubbles_csv_in, /*with_inversion_sites=*/true);
    if (all_bubbles.empty())
        throw std::runtime_error("describe: " + options.bubbles_csv_in + " lists no bubbles");
    // Without the boundaries' orientation a path crossing a site backwards, more than once, is not read
    // crossing by crossing.
    if (!all_bubbles.front().orientation_known)
        throw std::runtime_error("describe: " + options.bubbles_csv_in + " has no source_orient/sink_orient "
                                 "columns; rerun bubble (or panphorte/refine) to write a current site table");
    {
        std::unordered_set<std::size_t> seen;
        for (const Bubble& b : all_bubbles) {
            if (!seen.insert(b.id).second)
                throw std::runtime_error("describe: bubbles CSV has a duplicate bubble id " +
                                         std::to_string(b.id) + "; which site a feature belongs to "
                                         "would be undefined");
            auto need = [&](const std::string& n) {
                if (!graph.nodes.count(n))
                    throw std::runtime_error("describe: bubble " + std::to_string(b.id) + " names node '" +
                                             n + "', which is not in " + options.gfa_path +
                                             " -- the CSV and the graph do not describe the same data");
            };
            need(b.source);
            need(b.sink);
            for (const std::string& n : b.inside) need(n);
        }
        for (const std::size_t want : options.bubble_ids)
            if (!seen.count(want))
                throw std::runtime_error("describe: --bubble-id " + std::to_string(want) +
                                         " is not in " + options.bubbles_csv_in);
    }
    // Reference coordinates for the sidecars (-r): each node's span on the reference, and each site's.
    RefSpans ref_spans;
    if (!options.reference_path.empty()) {
        const std::string ref_name = resolve_reference_path_name(graph, options.reference_path, "describe");
        const PathRecord* ref = nullptr;
        for (const PathRecord& p : graph.paths) if (p.name == ref_name) ref = &p;
        std::size_t total = 0;
        for (const PathStep& st : ref->steps) total += graph.nodes.at(st.node_id).sequence.size();
        const ParsedReferencePath meta = parse_reference_path_label(ref_name, total, ref->type == 'W');
        ref_spans.known = true;
        ref_spans.chrom = meta.chrom;
        std::size_t at = meta.region_start_1based;
        for (const PathStep& st : ref->steps) {
            const std::size_t len = graph.nodes.at(st.node_id).sequence.size();
            if (len == 0) continue;
            auto [it, fresh] = ref_spans.node.try_emplace(st.node_id, at, at + len - 1);
            if (!fresh) { it->second.first = std::min(it->second.first, at); it->second.second = std::max(it->second.second, at + len - 1); }
            at += len;
        }
        for (const Bubble& b : all_bubbles) {
            const auto s0 = ref_spans.node.find(b.source), s1 = ref_spans.node.find(b.sink);
            if (s0 == ref_spans.node.end() || s1 == ref_spans.node.end()) continue;
            ref_spans.site[b.id] = {std::min(s0->second.first, s1->second.first), std::max(s0->second.second, s1->second.second)};
        }
    }
    // An inversion site is described on its own nodes, the inner sites' interiors left to them. Variant
    // rows are checked against its whole interior.
    const auto variant_members = bubble_member_nodes(all_bubbles);
    view_inversion_sites_on_own_nodes(all_bubbles);
    std::filesystem::create_directories(options.out_dir);
    std::unordered_set<std::size_t> bubble_filter;
    bubble_filter.reserve(options.bubble_ids.size() * 2 + 1);
    for (const std::size_t id : options.bubble_ids) {
        bubble_filter.insert(id);
    }

    // --variant-nodes: site -> its records' nodes; only those sites are described.
    std::unordered_map<std::size_t, std::unordered_set<std::string>> variant_nodes_by_bubble;
    const bool variant_mode = !options.variant_nodes_path.empty();
    if (variant_mode) {
        const VariantNodes vn = load_variant_nodes(options.variant_nodes_path, variant_members, "describe");
        for (const auto& [bid, recs] : vn.records) {
            auto& set = variant_nodes_by_bubble[bid];
            for (const CalledRecord& rec : recs) set.insert(rec.nodes.begin(), rec.nodes.end());
        }
        // A table with no record (call writes one for a run with none) would describe nothing.
        if (variant_nodes_by_bubble.empty())
            throw std::runtime_error("describe: --variant-nodes " + options.variant_nodes_path + " lists no "
                                     "variant, so nothing would be described; drop --variant-nodes to describe "
                                     "every bubble");
        // So would a selected site with no record.
        for (const std::size_t want : options.bubble_ids)
            if (!variant_nodes_by_bubble.count(want))
                throw std::runtime_error("describe: --bubble-id " + std::to_string(want) + " has no variant in " +
                                         options.variant_nodes_path + ", so --variant-nodes leaves nothing of it "
                                         "to describe");
    }

    const std::filesystem::path out_dir(options.out_dir);
    const std::string index_path = (out_dir / "describe.index.tsv").string();

    std::ofstream index_out(index_path);
    if (!index_out) {
        throw std::runtime_error("Failed to write describe index TSV: " + index_path);
    }
    index_out
        << "bubble_id\tstatus\tpaths\t"
        << "kmer_candidates\tkmer_kept\tkmer_discarded\t"
        << "feature_map_tsv_gz\tmatrix_tsv_gz\tcounts_jsonl_gz\tmatrix_written\tmatrix_reason\t"
        << "node_candidates\tnode_kept\tnode_discarded\t"
        << "edge_candidates\tedge_kept\tedge_discarded\t"
        << "graph_features_tsv_gz\tgraph_matrix_tsv_gz\tgraph_matrix_written\n";

    DescribeSummary summary;
    summary.files_written = 1; // the index
    summary.haplotypes = graph.paths.size();
    summary.nodes = graph.nodes.size();

    // The sites described: --bubble-id, and with --variant-nodes those with a record.
    std::vector<const Bubble*> to_process;
    for (const Bubble& bubble : all_bubbles) {
        if (!bubble_filter.empty() && bubble_filter.find(bubble.id) == bubble_filter.end()) continue;
        if (variant_mode && variant_nodes_by_bubble.find(bubble.id) == variant_nodes_by_bubble.end()) continue;
        to_process.push_back(&bubble);
    }
    // Nested or overlapping sites would count the same interior twice. Checked among the sites described,
    // so selecting one of such a pair works.
    {
        std::unordered_map<std::string, std::size_t> owner;
        for (const Bubble* b : to_process)
            for (const std::string& n : b->inside) {
                const auto [it, fresh] = owner.emplace(n, b->id);
                if (!fresh)
                    throw std::runtime_error("describe: bubbles " + std::to_string(it->second) + " and " +
                                             std::to_string(b->id) + " both claim interior node " + n +
                                             "; overlapping or nested sites would count its features twice "
                                             "-- describe a table without them, or select one with --bubble-id");
            }
    }
    // Only the described sites' nodes are indexed, as in call and inspect.
    std::vector<BubblePathIndex> path_indexes;
    {
        std::unordered_set<std::string> site_nodes;
        for (const Bubble* b : to_process) {
            site_nodes.insert(b->source);
            site_nodes.insert(b->sink);
            site_nodes.insert(b->inside.begin(), b->inside.end());
        }
        path_indexes.reserve(graph.paths.size());
        for (const auto& path : graph.paths) path_indexes.push_back(build_bubble_path_index(path, site_nodes));
    }

    // A node two described sites both walk -- one's sink and the next one's source, or a site's boundary
    // inside an inversion site -- is counted at one of them, its owner: the site holding it in its
    // interior, else the site with the smallest id. The owner counts every visit its crossings make; a
    // visit only the other site's crossings make, as when a path repeats that site alone, counts there.
    // K-mers follow the same rule among the sites that count k-mers, which an inversion site does not.
    std::vector<SiteScope> scopes(to_process.size());
    std::unordered_map<std::size_t, std::size_t> site_index;   // bubble id -> index in to_process
    for (std::size_t bi = 0; bi < to_process.size(); ++bi) site_index.emplace(to_process[bi]->id, bi);
    std::vector<std::unordered_set<std::string>> owned_shared(to_process.size());
    std::vector<SharedVisits> owner_visits_by_site(to_process.size());
    // Per site, the boundaries it shares for k-mers, owned or not; and each (site, owner, node, k-mers?)
    // where the site is not the owner.
    std::vector<std::unordered_set<std::string>> kmer_shared(to_process.size());
    struct SharedPair { std::size_t site, owner; std::string node; bool kmers; };
    std::vector<SharedPair> shared_pairs;
    for (const bool kmers : {false, true}) {
        std::unordered_map<std::string, std::size_t> interior_owner, boundary_owner, boundary_sites;
        for (const Bubble* b : to_process) {
            if (kmers && b->inversion_parent) continue;
            for (const std::string& n : b->inside) interior_owner[n] = b->id;
            for (const std::string* n : {&b->source, &b->sink}) {
                const auto [it, fresh] = boundary_owner.emplace(*n, b->id);
                if (!fresh) it->second = std::min(it->second, b->id);
                ++boundary_sites[*n];
            }
        }
        for (std::size_t bi = 0; bi < to_process.size(); ++bi) {
            const Bubble* b = to_process[bi];
            if (kmers && b->inversion_parent) continue;
            for (const std::string* n : {&b->source, &b->sink}) {
                const auto io = interior_owner.find(*n);
                const std::size_t owner = io != interior_owner.end() ? io->second : boundary_owner.at(*n);
                if (owner != b->id) {
                    const std::size_t oi = site_index.at(owner);
                    owned_shared[oi].insert(*n);
                    (kmers ? scopes[bi].shared_kmers : scopes[bi].shared_nodes)[*n] = &owner_visits_by_site[oi];
                    shared_pairs.push_back(SharedPair{bi, oi, *n, kmers});
                    if (kmers) { kmer_shared[bi].insert(*n); kmer_shared[oi].insert(*n); }
                } else if (kmers && boundary_sites.at(*n) > 1) {
                    scopes[bi].extend_kmers.insert(*n);
                }
            }
        }
    }
    {
        std::vector<std::size_t> owners;
        for (std::size_t bi = 0; bi < to_process.size(); ++bi)
            if (!owned_shared[bi].empty()) owners.push_back(bi);
        run_parallel(owners.size(), options.threads, [&](std::size_t k) {
            const std::size_t oi = owners[k];
            SharedVisits& visits = owner_visits_by_site[oi];
            visits.assign(graph.paths.size(), {});
            for (std::size_t p = 0; p < graph.paths.size(); ++p) {
                const BubbleWalk walk = canonical_steps_for_bubble_path(graph, path_indexes, *to_process[oi], p);
                for (std::size_t i = 0; i < walk.steps.size(); ++i)
                    if (owned_shared[oi].count(walk.steps[i].node_id)) visits[p].push_back(walk.path_pos[i]);
                std::sort(visits[p].begin(), visits[p].end());
            }
        });
    }
    // A shared boundary's visits are split when some path visits it on a crossing only the other site
    // makes; only then are its rows pooled whatever the filter says (see SiteScope), since otherwise the
    // owner counts every visit and the extra rows would only cost memory.
    {
        std::vector<std::size_t> sites;
        for (const SharedPair& sp : shared_pairs) sites.push_back(sp.site);
        std::sort(sites.begin(), sites.end());
        sites.erase(std::unique(sites.begin(), sites.end()), sites.end());
        std::vector<std::unordered_set<std::string>> split(to_process.size());
        run_parallel(sites.size(), options.threads, [&](std::size_t k) {
            const std::size_t bi = sites[k];
            const SiteScope& sc = scopes[bi];
            for (std::size_t p = 0; p < graph.paths.size(); ++p) {
                const BubbleWalk walk = canonical_steps_for_bubble_path(graph, path_indexes, *to_process[bi], p);
                for (std::size_t i = 0; i < walk.steps.size(); ++i) {
                    const std::string& id = walk.steps[i].node_id;
                    for (const auto* m : {&sc.shared_nodes, &sc.shared_kmers}) {
                        const auto it = m->find(id);
                        if (it != m->end() && !owner_visits(*it->second, p, walk.path_pos[i])) split[bi].insert(id);
                    }
                }
            }
        });
        for (const SharedPair& sp : shared_pairs) {
            if (!split[sp.site].count(sp.node)) continue;
            for (const std::size_t si : {sp.site, sp.owner})
                (sp.kmers ? scopes[si].pooled_kmer_nodes : scopes[si].pooled_nodes).insert(sp.node);
        }
    }
    for (std::size_t bi = 0; bi < to_process.size(); ++bi) {
        if (variant_mode) scopes[bi].variant_nodes = &variant_nodes_by_bubble.at(to_process[bi]->id);
        scopes[bi].flank_bp = options.variant_flank_bp;
        scopes[bi].kmer_size = options.kmer_size;
    }
    // Under --variant-nodes, the sites sharing a boundary for k-mers keep the union of their kept bases
    // on it: the owner then spells on its visits what the other site would, and the other site spells
    // what the owner would on the visits only it makes, so each k-mer inside the node counts once per visit.
    std::vector<SiteKeep> shared_site_keep(to_process.size());
    if (variant_mode) {
        std::vector<std::size_t> sharing;
        for (std::size_t bi = 0; bi < to_process.size(); ++bi)
            if (!kmer_shared[bi].empty()) sharing.push_back(bi);
        run_parallel(sharing.size(), options.threads, [&](std::size_t k) {
            const std::size_t bi = sharing[k];
            shared_site_keep[bi] = build_site_keep(graph, path_indexes, *to_process[bi], scopes[bi]);
        });
        std::unordered_map<std::string, NodeKeep> merged;
        for (const std::size_t bi : sharing)
            for (const std::string& n : kmer_shared[bi]) {
                const auto it = shared_site_keep[bi].find(n);
                if (it == shared_site_keep[bi].end()) continue;
                NodeKeep& m = merged[n];
                m.full = m.full || it->second.full;
                m.boundary = m.boundary || it->second.boundary;
                m.head = std::max(m.head, it->second.head);
                m.tail = std::max(m.tail, it->second.tail);
            }
        for (const std::size_t bi : sharing) {
            for (const std::string& n : kmer_shared[bi]) {
                const auto it = merged.find(n);
                if (it != merged.end()) shared_site_keep[bi][n] = it->second;
            }
            scopes[bi].site_keep = &shared_site_keep[bi];
        }
    }

    KmerPool kmer_pool;
    // Per-haplotype node/edge dosage, pooled across bubbles for BIMBAM and --samples.
    GraphPool graph_pool;
    // Which pools to fill.
    const bool want_kmer_pool = options.emit_kmers && options.bimbam;
    const bool want_graph_pool =
        options.emit_graph && (options.bimbam || !options.samples_path.empty());
    const bool want_bimbam = options.bimbam && (options.emit_kmers || options.emit_graph);
    // For NA: the sites each feature was kept at, and the paths crossing each site.
    std::unordered_map<std::uint64_t, std::vector<std::size_t>> feature_bubbles_k;
    std::unordered_map<std::string, std::vector<std::size_t>> feature_bubbles_g;
    std::unordered_map<std::uint64_t, std::set<std::string>> feature_nodes_k;  // pooled k-mer code -> node set
    std::unordered_map<std::size_t, std::vector<char>> bubble_traversers;   // bubble -> path columns crossing it
    cli::ProgressBar progress(options.quiet ? "" : "Describing bubbles", to_process.size());

    // Sites run in parallel, each writing its own files and filling its own pools; pools and index rows
    // are merged in site order, so the output does not depend on the thread count.
    std::vector<BubbleDescribeResult> results(to_process.size());
    std::vector<KmerPool> task_kmer(want_kmer_pool ? to_process.size() : 0);
    std::vector<GraphPool> task_graph(want_graph_pool ? to_process.size() : 0);
    // An inversion site has graph features only: between inverted repeats, the k-mers across each
    // breakpoint are the reverse complements of the reference's, the same canonical k-mers.
    DescribeOptions inversion_options = options;
    inversion_options.emit_kmers = false;
    // A site's pools and index row are merged once every earlier site has been, in site order, so the
    // output does not depend on the thread count and a finished site's pools are freed early.
    const auto merge_site = [&](std::size_t bi) {
        const Bubble& bubble = *to_process[bi];
        progress.tick();
        if (want_bimbam && !results[bi].traversing_paths.empty()) {
            auto& tset = bubble_traversers[bubble.id];
            tset.assign(graph.paths.size(), 0);
            for (const std::uint32_t p : results[bi].traversing_paths) tset[p] = 1;
        }
        if (want_kmer_pool) {
            for (auto& [code, carriers] : task_kmer[bi]) {
                add_carriers(kmer_pool[code], std::move(carriers));
                if (want_bimbam) feature_bubbles_k[code].push_back(bubble.id);
            }
            KmerPool().swap(task_kmer[bi]);
            if (want_bimbam) {
                for (auto& [code, ns] : results[bi].kmer_nodes) {
                    auto& dst = feature_nodes_k[code];
                    for (const std::string& n : ns) dst.insert(n);
                }
            }
        }
        if (want_graph_pool) {
            for (auto& [feature, carriers] : task_graph[bi]) {
                add_carriers(graph_pool[feature], std::move(carriers));
                if (want_bimbam) feature_bubbles_g[feature].push_back(bubble.id);
            }
            GraphPool().swap(task_graph[bi]);
        }
        const BubbleDescribeResult& result = results[bi];
        ++summary.bubbles_processed;
        if (result.paths > 0) {
            ++summary.bubbles_with_paths;
            summary.paths_written += result.paths;
            summary.features_written += result.features;
            summary.features_candidates += result.features_total;
            if (result.feature_map_path != ".") {
                summary.jsonl_files_written += 1;
                summary.files_written += 2; // feature map + sparse JSONL
                if (result.matrix_written) {
                    summary.matrix_files_written += 1;
                    summary.files_written += 1;
                }
            }
            summary.node_edge_features_written += result.node_features + result.edge_features;
            summary.node_edge_candidates += result.node_features_total + result.edge_features_total;
            if (options.emit_graph) {
                summary.files_written += 1; // graph feature map
                if (result.graph_matrix_written) {
                    summary.graph_matrix_files_written += 1;
                    summary.files_written += 1;
                }
            }
        }

        // Relative to the output directory, which is staged and renamed into place.
        const auto rel = [&](const std::string& p) {
            return p == "." ? p : std::filesystem::path(p).lexically_relative(out_dir).generic_string();
        };
        index_out
            << bubble.id << '\t'
            << result.status << '\t'
            << result.paths << '\t'
            << result.features_total << '\t'
            << result.features << '\t'
            << (result.features_total - result.features) << '\t'
            << rel(result.feature_map_path) << '\t'
            << rel(result.matrix_path) << '\t'
            << rel(result.counts_jsonl_path) << '\t'
            << (result.matrix_written ? "1" : "0") << '\t'
            << result.matrix_reason << '\t'
            << result.node_features_total << '\t'
            << result.node_features << '\t'
            << (result.node_features_total - result.node_features) << '\t'
            << result.edge_features_total << '\t'
            << result.edge_features << '\t'
            << (result.edge_features_total - result.edge_features) << '\t'
            << rel(result.graph_feature_map_path) << '\t'
            << rel(result.graph_matrix_path) << '\t'
            << (result.graph_matrix_written ? "1" : "0") << '\n';
        BubbleDescribeResult().kmer_nodes.swap(results[bi].kmer_nodes);
    };
    std::mutex merge_mutex;
    std::vector<char> done(to_process.size(), 0);
    std::size_t next_merge = 0;
    run_parallel(to_process.size(), options.threads, [&](std::size_t bi) {
        const Bubble& bubble = *to_process[bi];
        results[bi] = describe_one_bubble(bubble.inversion_parent ? inversion_options : options, graph,
                                          path_indexes, bubble, scopes[bi],
                                          want_kmer_pool ? &task_kmer[bi] : nullptr,
                                          want_graph_pool ? &task_graph[bi] : nullptr);
        if (bubble.inversion_parent && options.emit_kmers && results[bi].status == "ok")
            results[bi].matrix_reason = "skipped:inversion-site";
        const std::lock_guard<std::mutex> lock(merge_mutex);
        done[bi] = 1;
        while (next_merge < to_process.size() && done[next_merge]) merge_site(next_merge++);
    });
    if (next_merge != to_process.size()) throw std::runtime_error("describe: a site was not merged");
    cli::close_output(index_out, index_path);

    // Pooled matrices under haplotype/{kmers,graph}/. Rows the filter drops here are left out of the
    // sample level too, which also drops the rows constant across samples.
    std::unordered_set<std::string> pooled_dropped_k, pooled_dropped_g;
    if (want_bimbam) {
        std::vector<std::string> sample_order;
        sample_order.reserve(graph.paths.size());
        for (const auto& p : graph.paths) sample_order.push_back(p.name);
        const char* pool_header = "feature_id\tlayer\tencoding\tbubbles\tnodes\tchrom\tpos\tend\n";
        auto write_samples = [&](const std::filesystem::path& dir) {
            GzipWriter sout((dir / "samples.txt.gz").string());
            for (const std::string& s : sample_order) { sout.write(s); sout.write("\n"); }
            sout.close();
        };
        if (options.emit_kmers) {
            const std::filesystem::path dir = substrate_dir(out_dir, "haplotype", "kmers");
            write_samples(dir);
            GzipWriter annot((dir / "feature_annot.kmers.tsv.gz").string());
            annot.write(pool_header);
            std::vector<BimbamRow> rows;
            rows.reserve(kmer_pool.size());
            for (const auto& [code, carriers] : kmer_pool) {
                BimbamRow r;
                r.id = decode_kmer(code, options.kmer_size);
                r.encoding = options.feature_mode == DescribeFeatureMode::Syncmer ? "syncmer" : "all";
                const auto nit = feature_nodes_k.find(code);
                r.nodes = (nit != feature_nodes_k.end() && !nit->second.empty())
                    ? join_nodes(std::vector<std::string>(nit->second.begin(), nit->second.end())) : ".";
                const auto fb = feature_bubbles_k.find(code);
                if (fb != feature_bubbles_k.end()) r.bubbles = fb->second;
                r.carriers = &carriers;
                rows.push_back(std::move(r));
            }
            sort_rows(rows);
            write_bimbam_rows((dir / "bimbam_kmers.bimbam.gz").string(), annot, "kmer", ref_spans,
                              rows, sample_order, bubble_traversers, options.scale_dosage,
                              &options.min_feature_paths, &pooled_dropped_k);
            annot.close();
            summary.files_written += 3;
        }
        if (options.emit_graph) {
            const std::filesystem::path dir = substrate_dir(out_dir, "haplotype", "graph");
            write_samples(dir);
            GzipWriter annot((dir / "feature_annot.graph.tsv.gz").string());
            annot.write(pool_header);
            std::vector<BimbamRow> rows;
            rows.reserve(graph_pool.size());
            for (const auto& [feature, carriers] : graph_pool) {
                BimbamRow r;
                r.id = feature;
                r.nodes = feature;
                r.encoding = feature.find('>') != std::string::npos ? "edge" : "node";
                const auto fb = feature_bubbles_g.find(feature);
                if (fb != feature_bubbles_g.end()) r.bubbles = fb->second;
                r.carriers = &carriers;
                rows.push_back(std::move(r));
            }
            sort_rows(rows);
            write_bimbam_rows((dir / "bimbam_graph.bimbam.gz").string(), annot, "graph", ref_spans,
                              rows, sample_order, bubble_traversers, options.scale_dosage,
                              &options.min_feature_paths, &pooled_dropped_g);
            annot.close();
            summary.files_written += 3;
        }
    }

    // --samples: each pooled row summed over each sample's haplotypes (one listed twice counts twice).
    if (!options.samples_path.empty()) {
        std::unordered_set<std::string> graph_haps;
        for (const auto& pr : graph.paths) graph_haps.insert(pr.name);
        const CosigtTable cosigt = read_cosigt_table(options.samples_path, "describe", &graph_haps,
                                                     "a path in the input graph", options.quiet);
        const std::vector<std::string>& sample_order_s = cosigt.sample_order;
        // Each path's sample columns, once per slot naming it (a haplotype listed twice counts twice).
        std::unordered_map<std::string, std::uint32_t> path_col, sample_col;
        for (std::size_t i = 0; i < graph.paths.size(); ++i) path_col.emplace(graph.paths[i].name, static_cast<std::uint32_t>(i));
        for (std::size_t i = 0; i < sample_order_s.size(); ++i) sample_col.emplace(sample_order_s[i], static_cast<std::uint32_t>(i));
        std::vector<std::vector<std::uint32_t>> path_samples(graph.paths.size());
        for (const auto& [path, samples] : cosigt.path_to_samples)
            for (const std::string& sm : samples) path_samples[path_col.at(path)].push_back(sample_col.at(sm));
        // A pool's carriers summed into sample columns.
        std::vector<std::uint64_t> sum(sample_order_s.size(), 0);
        auto to_samples = [&](const Carriers& carriers) {
            std::vector<std::uint32_t> touched;
            for (const auto& [col, c] : carriers)
                for (const std::uint32_t sc : path_samples[col]) {
                    if (sum[sc] == 0) touched.push_back(sc);
                    sum[sc] += c;
                }
            std::sort(touched.begin(), touched.end());
            Carriers out;
            out.reserve(touched.size());
            for (const std::uint32_t sc : touched) {
                out.emplace_back(sc, static_cast<std::uint32_t>(
                    std::min<std::uint64_t>(sum[sc], std::numeric_limits<std::uint32_t>::max())));
                sum[sc] = 0;
            }
            return out;
        };

        if (want_bimbam) {
            // A sample's dosage is the sum over its haplotypes, so it is complete at a bubble only when
            // each of them is observed there: exactly one of a haplotype's contigs traverses it. The pools
            // sum every contig's counts, which at a bubble where the sample is known are its one
            // traversing contig's.
            std::unordered_map<std::size_t, std::vector<char>> bubble_traversers_s;
            for (const auto& [b, crossing] : bubble_traversers) {
                auto& dst = bubble_traversers_s[b];
                dst.assign(sample_order_s.size(), 0);
                const std::vector<char>& tv = crossing;
                const auto observed = [&](const std::string& c) {
                    const auto it = path_col.find(c);
                    return it != path_col.end() && tv[it->second] != 0;
                };
                for (const auto& [sm, slots] : cosigt.slots) {
                    bool all = !slots.empty();
                    for (const auto& slot : slots)
                        if (slot_contig(slot, observed) == nullptr) { all = false; break; }
                    if (all) dst[sample_col.at(sm)] = 1;
                }
            }
            const char* pool_header = "feature_id\tlayer\tencoding\tbubbles\tnodes\tchrom\tpos\tend\n";
            // A row varying across haplotypes can be constant across samples (carried only by paths in no
            // sample, or by every sample alike): such a row is dropped, whatever --min-paths.
            const std::size_t constant_only = 0;
            auto write_samples_s = [&](const std::filesystem::path& dir) {
                GzipWriter sout((dir / "samples.txt.gz").string());
                for (const std::string& s : sample_order_s) { sout.write(s); sout.write("\n"); }
                sout.close();
            };
            if (options.emit_kmers) {
                const std::filesystem::path dir = substrate_dir(out_dir, "sample", "kmers");
                write_samples_s(dir);
                GzipWriter annot((dir / "feature_annot.kmers.tsv.gz").string());
                annot.write(pool_header);
                std::vector<Carriers> sample_carriers;
                sample_carriers.reserve(kmer_pool.size());
                std::vector<BimbamRow> rows; rows.reserve(kmer_pool.size());
                for (const auto& [code, carriers] : kmer_pool) {
                    sample_carriers.push_back(to_samples(carriers));
                    BimbamRow r; r.id = decode_kmer(code, options.kmer_size);
                    r.encoding = options.feature_mode == DescribeFeatureMode::Syncmer ? "syncmer" : "all";
                    const auto nit = feature_nodes_k.find(code);
                    r.nodes = (nit != feature_nodes_k.end() && !nit->second.empty())
                        ? join_nodes(std::vector<std::string>(nit->second.begin(), nit->second.end())) : ".";
                    const auto fb = feature_bubbles_k.find(code);
                    if (fb != feature_bubbles_k.end()) r.bubbles = fb->second;
                    r.carriers = &sample_carriers.back();
                    rows.push_back(std::move(r));
                }
                sort_rows(rows);
                write_bimbam_rows((dir / "bimbam_kmers.bimbam.gz").string(), annot, "kmer", ref_spans,
                                  rows, sample_order_s, bubble_traversers_s, options.scale_dosage,
                                  &constant_only, nullptr, &pooled_dropped_k);
                annot.close();
                summary.files_written += 3;
            }
            if (options.emit_graph) {
                const std::filesystem::path dir = substrate_dir(out_dir, "sample", "graph");
                write_samples_s(dir);
                GzipWriter annot((dir / "feature_annot.graph.tsv.gz").string());
                annot.write(pool_header);
                std::vector<Carriers> sample_carriers;
                sample_carriers.reserve(graph_pool.size());
                std::vector<BimbamRow> rows; rows.reserve(graph_pool.size());
                for (const auto& [feature, carriers] : graph_pool) {
                    sample_carriers.push_back(to_samples(carriers));
                    BimbamRow r; r.id = feature; r.nodes = feature;
                    r.encoding = feature.find('>') != std::string::npos ? "edge" : "node";
                    const auto fb = feature_bubbles_g.find(feature);
                    if (fb != feature_bubbles_g.end()) r.bubbles = fb->second;
                    r.carriers = &sample_carriers.back();
                    rows.push_back(std::move(r));
                }
                sort_rows(rows);
                write_bimbam_rows((dir / "bimbam_graph.bimbam.gz").string(), annot, "graph", ref_spans,
                                  rows, sample_order_s, bubble_traversers_s, options.scale_dosage,
                                  &constant_only, nullptr, &pooled_dropped_g);
                annot.close();
                summary.files_written += 3;
            }
            summary.samples = sample_order_s.size();
        }
    }

    if (summary_out != nullptr) {
        *summary_out = summary;
    }
}

// The variant substrate, read from the VCF alone, so --only-variant needs no GFA.
void describe_variant_from_vcf(const DescribeOptions& options, DescribeSummary* summary_out) {
    DescribeSummary local;
    DescribeSummary& summary = summary_out ? *summary_out : local;
    if (!options.emit_kmers && !options.emit_graph) std::filesystem::create_directories(options.out_dir);
    emit_variant_substrate(options, summary);
}

} // namespace panvar
