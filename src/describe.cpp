#include "panvar/describe.hpp"

#include "panvar/syncmer.hpp"

#include "panvar/bubble_path.hpp"
#include "panvar/cli_utils.hpp"
#include "panvar/gfa.hpp"
#include "panvar/graph_utils.hpp"
#include "panvar/output.hpp"
#include "panvar/parallel.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include <zlib.h>

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
    std::size_t features = 0;            // discriminative k-mers kept
    std::size_t features_total = 0;      // k-mer candidates before filter
    bool matrix_written = false;
    std::string status = "ok";
    std::string matrix_reason = ".";
    std::string feature_map_path = ".";
    std::string matrix_path = ".";
    std::string counts_jsonl_path = ".";
    std::size_t node_features = 0;       // discriminative nodes kept
    std::size_t node_features_total = 0; // node candidates before filter
    std::size_t edge_features = 0;       // discriminative edges kept
    std::size_t edge_features_total = 0; // edge candidates before filter
    bool graph_matrix_written = false;
    std::string graph_feature_map_path = ".";
    std::string graph_matrix_path = ".";
    std::vector<std::string> traversing_paths;  // path names that traverse this bubble (for BIMBAM NA)
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

std::string feature_mode_name(DescribeFeatureMode mode) {
    switch (mode) {
        case DescribeFeatureMode::AllKmers:
            return "all";
        case DescribeFeatureMode::Syncmer:
            return "syncmer";
    }
    return "all";
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

// Keep a feature if it's informative: copy-number features (count varies) always pass; otherwise it
// needs a symmetric minor-presence cut min(present, absent) > min_paths. min_paths==0 keeps near-all.
bool feature_passes_filter(
    std::size_t present_paths,
    std::uint32_t min_nonzero_count,
    std::uint32_t max_count,
    std::size_t path_count,
    std::size_t min_paths) {

    if (path_count == 0) {
        return true;
    }
    if (min_nonzero_count != max_count) {
        return true; // copy-number signal: informative regardless of presence breadth
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

// Node/edge dosage: the graph-coordinate substrate paralleling the k-mers - per-path traversal counts
// per inside node and per oriented edge. This is a descriptive count, not a CN call (a node revisited
// elsewhere also counts >1); the edge layer carries the real adjacency/tandem signal. CN calling lives
// in variant_call.cpp.

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

std::string node_feature_name(std::size_t feature_id) {
    return "N" + std::to_string(feature_id);
}

std::string edge_feature_name(std::size_t feature_id) {
    return "E" + std::to_string(feature_id);
}

std::string oriented_node_token(const PathStep& step) {
    return step.node_id + (step.reverse ? '-' : '+');
}

std::string edge_key(const PathStep& from, const PathStep& to) {
    return oriented_node_token(from) + ">" + oriented_node_token(to);
}

// Forward decl (defined below): per-step keep mask for --variant-nodes (+flank). Used by the graph
// dosage counters here so they honour the variant-restricted scope, same as the k-mer substrate.
std::vector<char> variant_keep_mask(
    const Graph& graph,
    const std::vector<PathStep>& steps,
    const std::unordered_set<std::string>& variant_nodes,
    std::size_t flank_bp);

std::vector<PathStep> canonical_steps_for_bubble_path(
    const Bubble& bubble,
    const PathRecord& path,
    const BubblePathIndex& index) {

    const auto interval = find_best_bubble_path_interval(index, bubble);
    if (!interval.has_value()) {
        return {};
    }
    return canonical_bubble_path_steps(path, bubble, *interval);
}

// Count node/edge dosage over a path's bubble walk. When `keep` is non-null (variant-restricted mode)
// a node is counted only on kept steps, and an edge only when BOTH of its endpoints are kept, so the
// graph substrate honours --variant-nodes (+flank) exactly like the k-mer substrate.
GraphDosage count_graph_dosage(const std::vector<PathStep>& steps,
                               const std::vector<char>* keep = nullptr) {
    GraphDosage dosage;
    for (std::size_t i = 0; i < steps.size(); ++i) {
        if (keep != nullptr && !(*keep)[i]) continue;
        ++dosage.node_counts[steps[i].node_id];
    }
    for (std::size_t i = 0; i + 1 < steps.size(); ++i) {
        if (keep != nullptr && (!(*keep)[i] || !(*keep)[i + 1])) continue;
        ++dosage.edge_counts[edge_key(steps[i], steps[i + 1])];
    }
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

// Drop features present in every path with identical dosage (no variance ->
// useless for association). Input map iterates sorted, so output is sorted too.
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
    const Graph& graph,
    const Bubble& bubble,
    const std::vector<BubblePathIndex>& path_indexes,
    const std::vector<PathMeta>& paths,
    std::map<std::string, GraphFeatureStat>& node_stats,
    std::map<std::string, GraphFeatureStat>& edge_stats,
    const std::unordered_set<std::string>* variant_nodes,
    std::size_t flank_bp) {

    for (const PathMeta& meta : paths) {
        const std::vector<PathStep> steps = canonical_steps_for_bubble_path(
            bubble, graph.paths[meta.path_index], path_indexes[meta.path_index]);
        if (steps.empty()) {
            continue;
        }
        const std::vector<char> keep =
            variant_nodes != nullptr ? variant_keep_mask(graph, steps, *variant_nodes, flank_bp)
                                     : std::vector<char>();
        const GraphDosage dosage = count_graph_dosage(steps, variant_nodes != nullptr ? &keep : nullptr);
        update_graph_stats(dosage.node_counts, node_stats);
        update_graph_stats(dosage.edge_counts, edge_stats);
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
    const Graph& graph,
    const Bubble& bubble,
    const std::vector<BubblePathIndex>& path_indexes,
    const std::vector<PathMeta>& paths,
    const std::vector<std::string>& node_features,
    const std::vector<std::string>& edge_features,
    const std::unordered_set<std::string>* variant_nodes,
    std::size_t flank_bp) {

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

    for (const PathMeta& meta : paths) {
        const std::vector<PathStep> steps = canonical_steps_for_bubble_path(
            bubble, graph.paths[meta.path_index], path_indexes[meta.path_index]);
        if (steps.empty()) {
            continue;
        }
        const std::vector<char> keep =
            variant_nodes != nullptr ? variant_keep_mask(graph, steps, *variant_nodes, flank_bp)
                                     : std::vector<char>();
        const GraphDosage dosage = count_graph_dosage(steps, variant_nodes != nullptr ? &keep : nullptr);

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

// Per-step keep mask for `--variant-nodes`. Keep a step (emit its real sequence) when it is a variant
// node, or -- with flank_bp > 0 -- within flank_bp along the path of one, so k-mers spanning into the
// flanking sequence survive. flank_bp == 0 is the strict variant-only mask. Path-specific, so per path.
std::vector<char> variant_keep_mask(
    const Graph& graph,
    const std::vector<PathStep>& steps,
    const std::unordered_set<std::string>& variant_nodes,
    std::size_t flank_bp) {

    const std::size_t n = steps.size();
    std::vector<char> is_var(n, 0), keep(n, 0);
    std::vector<std::size_t> len(n, 0);
    for (std::size_t i = 0; i < n; ++i) {
        const auto it = graph.nodes.find(steps[i].node_id);
        len[i] = (it == graph.nodes.end()) ? 0 : it->second.sequence.size();
        is_var[i] = variant_nodes.count(steps[i].node_id) ? 1 : 0;
        keep[i] = is_var[i];
    }
    if (flank_bp == 0) {
        return keep;
    }
    // Forward sweep: `gap` is the bp from the right edge of the last variant node to the left edge
    // of the current node; keep while still within flank_bp. Backward sweep mirrors it for the left
    // flank. SIZE_MAX marks "no variant node seen yet" so leading/trailing runs are not kept.
    std::size_t gap = std::numeric_limits<std::size_t>::max();
    for (std::size_t i = 0; i < n; ++i) {
        if (is_var[i]) { gap = 0; continue; }
        if (gap != std::numeric_limits<std::size_t>::max()) {
            if (gap < flank_bp) keep[i] = 1;
            gap += len[i];
        }
    }
    gap = std::numeric_limits<std::size_t>::max();
    for (std::size_t i = n; i-- > 0;) {
        if (is_var[i]) { gap = 0; continue; }
        if (gap != std::numeric_limits<std::size_t>::max()) {
            if (gap < flank_bp) keep[i] = 1;
            gap += len[i];
        }
    }
    return keep;
}

// When `variant_nodes` is non-null, bases spelled from nodes outside the keep mask are replaced
// by 'N' (same length). Because the k-mer scanner resets its window on any non-ACGT base, this
// confines all k-mer/syncmer generation to the kept runs — without disturbing segment offsets or
// node provenance. (`--variant-nodes` restriction, optionally widened by `flank_bp`.)
std::string path_sequence_for_bubble(
    const Graph& graph,
    const Bubble& bubble,
    const PathRecord& path,
    const BubblePathIndex& index,
    bool* complete_out,
    const std::unordered_set<std::string>* variant_nodes = nullptr,
    std::size_t flank_bp = 0) {

    if (complete_out != nullptr) {
        *complete_out = false;
    }
    const auto interval = find_best_bubble_path_interval(index, bubble);
    if (!interval.has_value()) {
        return {};
    }
    const auto steps = canonical_bubble_path_steps(path, bubble, *interval);
    if (steps.empty()) {
        return {};
    }
    if (variant_nodes == nullptr) {
        bool complete = false;
        std::string sequence = spell_path_steps_sequence(graph, steps, &complete);
        if (complete_out != nullptr) {
            *complete_out = complete;
        }
        return complete ? sequence : std::string{};
    }
    const std::vector<char> keep = variant_keep_mask(graph, steps, *variant_nodes, flank_bp);
    std::string sequence;
    for (std::size_t si = 0; si < steps.size(); ++si) {
        const PathStep& step = steps[si];
        const auto node_it = graph.nodes.find(step.node_id);
        if (node_it == graph.nodes.end() || node_it->second.sequence.empty()) {
            return {};
        }
        const auto& node_seq = node_it->second.sequence;
        if (keep[si]) {
            sequence += step.reverse ? reverse_complement(node_seq) : node_seq;
        } else {
            sequence.append(node_seq.size(), 'N');
        }
    }
    if (complete_out != nullptr) {
        *complete_out = true;
    }
    return sequence;
}

BubblePathSequence path_sequence_with_segments_for_bubble(
    const Graph& graph,
    const Bubble& bubble,
    const PathRecord& path,
    const BubblePathIndex& index,
    const std::unordered_set<std::string>* variant_nodes = nullptr,
    std::size_t flank_bp = 0) {

    BubblePathSequence out;
    const auto interval = find_best_bubble_path_interval(index, bubble);
    if (!interval.has_value()) {
        return out;
    }
    const auto steps = canonical_bubble_path_steps(path, bubble, *interval);
    if (steps.empty()) {
        return out;
    }

    std::size_t total_len = 0;
    for (const PathStep& step : steps) {
        const auto node_it = graph.nodes.find(step.node_id);
        if (node_it == graph.nodes.end() || node_it->second.sequence.empty()) {
            return out;
        }
        total_len += node_it->second.sequence.size();
    }

    const std::vector<char> keep =
        variant_nodes != nullptr ? variant_keep_mask(graph, steps, *variant_nodes, flank_bp)
                                 : std::vector<char>{};
    out.sequence.reserve(total_len);
    out.segments.reserve(steps.size());
    for (std::size_t si = 0; si < steps.size(); ++si) {
        const PathStep& step = steps[si];
        const auto& node_seq = graph.nodes.at(step.node_id).sequence;
        const std::size_t start = out.sequence.size();
        if (variant_nodes != nullptr && !keep[si]) {
            out.sequence.append(node_seq.size(), 'N');     // masked: no k-mers from non-kept nodes
        } else if (step.reverse) {
            out.sequence += reverse_complement(node_seq);
        } else {
            out.sequence += node_seq;
        }
        out.segments.push_back(SequenceSegment{start, out.sequence.size(), step.node_id});
    }
    out.complete = true;
    return out;
}

std::vector<PathMeta> collect_path_metadata_and_stats(
    const Graph& graph,
    const Bubble& bubble,
    const std::vector<BubblePathIndex>& path_indexes,
    const DescribeOptions& options,
    std::unordered_map<std::uint64_t, KmerStats>& stats,
    const std::unordered_set<std::string>* variant_nodes = nullptr) {

    std::vector<PathMeta> paths;
    paths.reserve(graph.paths.size());
    std::unordered_map<std::uint64_t, std::uint32_t> counts;

    for (std::size_t i = 0; i < graph.paths.size(); ++i) {
        bool complete = false;
        const std::string sequence = path_sequence_for_bubble(graph, bubble, graph.paths[i], path_indexes[i], &complete, variant_nodes, options.variant_flank_bp);
        if (!complete) {
            continue;
        }

        const auto [sample, haplotype] = parse_sample_haplotype(graph.paths[i].name);
        PathMeta meta;
        meta.path_index = i;
        meta.path_name = graph.paths[i].name;
        meta.sample = sample;
        meta.haplotype = haplotype;
        meta.length_bp = sequence.size();
        paths.push_back(std::move(meta));

        count_selected_features(sequence, options, counts);
        update_kmer_stats(counts, stats);
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
    const std::unordered_set<std::string>* variant_nodes = nullptr) {

    GzipWriter out(path);
    std::unordered_map<std::uint64_t, KmerPathCount> counts;

    for (const PathMeta& meta : paths) {
        const BubblePathSequence path_sequence = path_sequence_with_segments_for_bubble(
            graph,
            bubble,
            graph.paths[meta.path_index],
            path_indexes[meta.path_index],
            variant_nodes,
            options.variant_flank_bp);
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
    const Graph& graph,
    const Bubble& bubble,
    const std::vector<BubblePathIndex>& path_indexes,
    const std::vector<PathMeta>& paths,
    const std::vector<std::uint64_t>& features,
    const DescribeOptions& options,
    const std::unordered_set<std::string>* variant_nodes = nullptr) {

    GzipWriter out(path);
    std::string header = "bubble_id\tsample\thaplotype\tpath_name";
    for (std::size_t i = 0; i < features.size(); ++i) {
        header.push_back('\t');
        header += feature_name(i + 1);
    }
    header.push_back('\n');
    out.write(header);

    std::unordered_map<std::uint64_t, std::uint32_t> counts;
    for (const PathMeta& meta : paths) {
        bool complete = false;
        const std::string sequence = path_sequence_for_bubble(
            graph,
            bubble,
            graph.paths[meta.path_index],
            path_indexes[meta.path_index],
            &complete,
            variant_nodes,
            options.variant_flank_bp);
        if (!complete) {
            continue;
        }
        count_selected_features(sequence, options, counts);

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
            const auto it = counts.find(code);
            if (it == counts.end()) {
                line.push_back('0');
            } else {
                line += std::to_string(it->second);
            }
        }
        line.push_back('\n');
        out.write(line);
    }
    out.close();
}

void write_params_json(const DescribeOptions& options, const std::string& path) {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Failed to write describe params JSON: " + path);
    }
    out << "{\n"
        << "  \"mode\": \"canonical_kmer_counts\",\n"
        << "  \"gfa\": \"" << json_escape(options.gfa_path) << "\",\n"
        << "  \"bubbles_csv\": \"" << json_escape(options.bubbles_csv_in) << "\",\n"
        << "  \"kmer_size\": " << options.kmer_size << ",\n"
        << "  \"feature_mode\": \"" << feature_mode_name(options.feature_mode) << "\",\n"
        << "  \"syncmer_s\": " << effective_syncmer_s(options) << ",\n"
        << "  \"canonical_reverse_complement\": true,\n"
        << "  \"non_discriminative_kmers_removed\": true,\n"
        << "  \"min_feature_paths\": " << options.min_feature_paths << ",\n"
        << "  \"copy_number_features_exempt_from_min_paths\": true,\n"
        << "  \"node_edge_dosage_tables\": true,\n"
        << "  \"node_provenance\": \"feature_map\",\n"
        << "  \"sparse_jsonl_tuple\": \"[feature_id,count]\",\n"
        << "  \"max_wide_features\": " << options.max_wide_features << ",\n"
        << "  \"wide_matrix_requested\": " << (options.write_wide_matrix ? "true" : "false") << ",\n"
        << "  \"force_wide_matrix\": " << (options.force_wide_matrix ? "true" : "false") << "\n"
        << "}\n";
}

// K-mer cohort pool: kmer code -> (path name -> summed count), over discriminative features only,
// aggregated across bubbles (the same canonical k-mer in two bubbles sums). Feeds the BIMBAM export.
using KmerPool = std::unordered_map<std::uint64_t, std::unordered_map<std::string, std::uint64_t>>;

void accumulate_kmer_counts(
    const Graph& graph,
    const Bubble& bubble,
    const std::vector<BubblePathIndex>& path_indexes,
    const std::vector<PathMeta>& paths,
    const std::vector<std::uint64_t>& features,
    const DescribeOptions& options,
    const std::unordered_set<std::string>* variant_nodes,
    KmerPool& pool) {

    const std::unordered_set<std::uint64_t> keep(features.begin(), features.end());
    if (keep.empty()) {
        return;
    }
    std::unordered_map<std::uint64_t, std::uint32_t> counts;
    for (const PathMeta& meta : paths) {
        bool complete = false;
        const std::string sequence = path_sequence_for_bubble(
            graph, bubble, graph.paths[meta.path_index], path_indexes[meta.path_index], &complete, variant_nodes, options.variant_flank_bp);
        if (!complete) {
            continue;
        }
        count_selected_features(sequence, options, counts);
        for (const auto& [code, c] : counts) {
            if (c > 0 && keep.count(code) != 0) {
                pool[code][meta.path_name] += c;
            }
        }
    }
}

// Per-feature, per-haplotype node/edge dosage pool: the graph-substrate analogue of KmerPool,
// keyed by the real graph coordinate (a node id, or a "from>to" edge key) so it stays traceable
// (node features join to call's variant_nodes.tsv by node id). Only discriminative features are pooled.
using GraphPool = std::map<std::string, std::map<std::string, std::uint32_t>>;

void accumulate_graph_counts(
    const Graph& graph,
    const Bubble& bubble,
    const std::vector<BubblePathIndex>& path_indexes,
    const std::vector<PathMeta>& paths,
    const std::vector<std::string>& node_features,
    const std::vector<std::string>& edge_features,
    GraphPool& pool,
    const std::unordered_set<std::string>* variant_nodes,
    std::size_t flank_bp) {

    const std::unordered_set<std::string> keep_nodes(node_features.begin(), node_features.end());
    const std::unordered_set<std::string> keep_edges(edge_features.begin(), edge_features.end());
    if (keep_nodes.empty() && keep_edges.empty()) {
        return;
    }
    for (const PathMeta& meta : paths) {
        const std::vector<PathStep> steps = canonical_steps_for_bubble_path(
            bubble, graph.paths[meta.path_index], path_indexes[meta.path_index]);
        if (steps.empty()) {
            continue;
        }
        const std::vector<char> keep =
            variant_nodes != nullptr ? variant_keep_mask(graph, steps, *variant_nodes, flank_bp)
                                     : std::vector<char>();
        const GraphDosage dosage = count_graph_dosage(steps, variant_nodes != nullptr ? &keep : nullptr);
        for (const auto& [key, c] : dosage.node_counts) {
            if (c > 0 && keep_nodes.count(key) != 0) pool[key][meta.path_name] += c;
        }
        for (const auto& [key, c] : dosage.edge_counts) {
            if (c > 0 && keep_edges.count(key) != 0) pool[key][meta.path_name] += c;
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

// BIMBAM mean-genotype dosage export (GEMMA / associate). One normalized row per feature: a carrier
// holds its dosage, a sample that traverses the feature's bubble(s) but doesn't carry it gets 0, and
// one that traverses none gets NA (missing, distinct from absent).
struct BimbamRow {
    std::string id;                                   // decoded k-mer or node/edge id (globally unique)
    std::string nodes;                                // provenance ("." for k-mers; the id itself for graph)
    std::vector<std::size_t> bubbles;                 // bubble(s) the feature was discriminative in
    std::unordered_map<std::string, double> carriers; // path -> dosage
};

// Per-substrate output directory: <out>/{haplotype|sample}/{kmers|graph|variant}/, created on demand.
// Each holds a self-contained associate input -- the BIMBAM matrix, its feature_annot, and the column
// order (samples.txt.gz) -- so the substrates and the aggregated-vs-per-haplotype split never co-mingle.
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
    const std::vector<BimbamRow>& rows,
    const std::vector<std::string>& sample_order,
    const std::unordered_map<std::size_t, std::unordered_set<std::string>>& bubble_traversers,
    bool scale_dosage) {

    GzipWriter geno(bimbam_path);
    for (const BimbamRow& row : rows) {
        // annotation sidecar: feature_id, layer, encoding, bubbles, nodes
        std::string a = tsv_sanitize(row.id);
        a += '\t'; a += layer;
        a += "\tcount\t";
        for (std::size_t i = 0; i < row.bubbles.size(); ++i) { if (i) a += ';'; a += std::to_string(row.bubbles[i]); }
        a += '\t'; a += row.nodes; a += '\n';
        annot.write(a);

        // Per-sample dosage (NaN = missing/NA): a carrier's count, 0 if it traverses the bubble but
        // does not carry the feature, else NA.
        std::vector<double> vals(sample_order.size(), std::numeric_limits<double>::quiet_NaN());
        double lo = std::numeric_limits<double>::infinity(), hi = -std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i < sample_order.size(); ++i) {
            const auto cit = row.carriers.find(sample_order[i]);
            if (cit != row.carriers.end()) {
                vals[i] = cit->second;
            } else {
                for (std::size_t b : row.bubbles) {
                    const auto bt = bubble_traversers.find(b);
                    if (bt != bubble_traversers.end() && bt->second.count(sample_order[i])) { vals[i] = 0.0; break; }
                }
            }
            if (std::isfinite(vals[i])) { lo = std::min(lo, vals[i]); hi = std::max(hi, vals[i]); }
        }
        // BIMBAM mean-genotype line: id, allele1, allele2, dosages... (--scale-dosage maps to 0..2)
        const bool do_scale = scale_dosage && std::isfinite(lo) && hi > lo;
        std::string line = tsv_sanitize(row.id);
        line += ", A, B";
        for (double v : vals) {
            line += ", ";
            if (!std::isfinite(v)) line += "NA";
            else line += format_dosage(do_scale ? (v - lo) / (hi - lo) * 2.0 : v);
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
    const std::unordered_set<std::string>* variant_nodes = nullptr,
    KmerPool* kmer_pool = nullptr,
    GraphPool* graph_pool = nullptr) {

    BubbleDescribeResult result;
    const std::filesystem::path out_dir(options.out_dir);
    const std::filesystem::path bubble_dir = out_dir / ("bubble_" + std::to_string(bubble.id));
    std::filesystem::create_directories(bubble_dir);
    result.feature_map_path = (bubble_dir / "kmer_features.tsv.gz").string();
    result.matrix_path = (bubble_dir / "kmer_matrix.tsv.gz").string();
    result.counts_jsonl_path = (bubble_dir / "kmer_counts.jsonl.gz").string();
    result.graph_feature_map_path = (bubble_dir / "graph_features.tsv.gz").string();
    result.graph_matrix_path = (bubble_dir / "graph_matrix.tsv.gz").string();

    std::unordered_map<std::uint64_t, KmerStats> stats;
    stats.reserve(4096);
    const std::vector<PathMeta> paths = collect_path_metadata_and_stats(
        graph,
        bubble,
        path_indexes,
        options,
        stats,
        variant_nodes);

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
    result.traversing_paths.reserve(paths.size());
    for (const PathMeta& m : paths) result.traversing_paths.push_back(m.path_name);

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
            variant_nodes);
        write_feature_map(result.feature_map_path, features, stats, feature_nodes, options.kmer_size, paths.size());

        // Carry each kept k-mer's node provenance up to the pooled BIMBAM/feature_annot (so a pooled
        // k-mer keeps a node link for traceback and gene annotation; per-bubble detail stays in the map).
        if (kmer_pool != nullptr) {
            for (const std::uint64_t code : features) {
                const auto it = feature_id_by_code.find(code);
                if (it != feature_id_by_code.end() && it->second < feature_nodes.size())
                    result.kmer_nodes.emplace(code, sorted_nodes(feature_nodes[it->second]));
            }
        }

        if (kmer_pool != nullptr) {
            accumulate_kmer_counts(graph, bubble, path_indexes, paths, features, options, variant_nodes, *kmer_pool);
        }

        const bool wide_allowed_by_cap =
            options.max_wide_features == 0 || features.size() <= options.max_wide_features;
        if (options.write_wide_matrix && (options.force_wide_matrix || wide_allowed_by_cap)) {
            write_wide_matrix(result.matrix_path, graph, bubble, path_indexes, paths, features, options, variant_nodes);
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
        accumulate_graph_stats(graph, bubble, path_indexes, paths, node_stats, edge_stats,
                               variant_nodes, options.variant_flank_bp);
        const std::vector<std::string> node_features =
            select_discriminative_graph_features(node_stats, paths.size(), options.min_feature_paths);
        const std::vector<std::string> edge_features =
            select_discriminative_graph_features(edge_stats, paths.size(), options.min_feature_paths);
        result.node_features = node_features.size();
        result.node_features_total = node_stats.size();
        result.edge_features = edge_features.size();
        result.edge_features_total = edge_stats.size();
        if (graph_pool != nullptr) {
            accumulate_graph_counts(graph, bubble, path_indexes, paths, node_features, edge_features,
                                    *graph_pool, variant_nodes, options.variant_flank_bp);
        }
        write_graph_feature_map(
            result.graph_feature_map_path, node_features, edge_features, node_stats, edge_stats, paths.size());
        if (options.write_wide_matrix) {
            write_graph_matrix(
                result.graph_matrix_path, graph, bubble, path_indexes, paths, node_features, edge_features,
                variant_nodes, options.variant_flank_bp);
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

// Parse a cosigt sample->haplotype table: "sample<TAB>hap1,hap2<TAB>..." (header row tolerated).
// Returns (haplotype-path -> samples carrying it, sample column order). Shared by the k-mer/graph
// SAMPLE-level export and the variant-level export so the diploid summation is defined in one place.
std::pair<std::unordered_map<std::string, std::vector<std::string>>, std::vector<std::string>>
read_cosigt_table(const std::string& path) {
    std::ifstream sin(path);
    if (!sin) throw std::runtime_error("Failed to open --samples file: " + path);
    std::unordered_map<std::string, std::vector<std::string>> path_to_samples;
    std::vector<std::string> sample_order;
    std::string line;
    bool first = true;
    while (std::getline(sin, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        std::vector<std::string> fields = split_on(line, '\t');
        if (fields.empty()) continue;
        if (first) {
            first = false;
            std::string h = fields[0];
            for (char& c : h) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
            if (h == "sample" || h == "sample.id" || h == "sample_id" || (!h.empty() && h[0] == '#')) continue;
        }
        const std::string& sample = fields[0];
        if (sample.empty()) continue;
        sample_order.push_back(sample);
        for (std::size_t fi = 1; fi < fields.size(); ++fi)
            for (const std::string& hap : split_on(fields[fi], ','))
                if (!hap.empty()) path_to_samples[hap].push_back(sample);
    }
    return {std::move(path_to_samples), std::move(sample_order)};
}

// Variant-level BIMBAM: the SV calls as the GWAS unit. Reads call's VCF (FORMAT GT[:CN]), one dosage
// row per variant (per ALT for multiallelic): CN (DUP) / allele indicator (multiallelic) / GT 0|1 else;
// GT=. -> NA. Same format/NA semantics as the k-mer/graph BIMBAM, but the unit is the variant -- an
// honest testing denominator (features within a variant are correlated).
struct VariantBimbam {
    std::string id, svtype, bubble, nodes, gene = ".", af = ".", an = ".";
    std::unordered_map<std::string, double> dose;  // haplotype -> dosage (present only; absent => NA)
};

void emit_variant_substrate(const DescribeOptions& options, DescribeSummary& summary) {
    std::ifstream vin(options.variant_vcf_path);
    if (!vin) throw std::runtime_error("describe --variant-vcf: cannot open " + options.variant_vcf_path);
    const std::filesystem::path out_dir(options.out_dir);
    std::filesystem::create_directories(out_dir);

    std::vector<std::string> hap_order;     // VCF sample columns (haplotype paths) = BIMBAM column order
    std::vector<VariantBimbam> rows;
    std::string line;
    while (std::getline(vin, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        if (line[0] == '#') {
            if (line.rfind("#CHROM", 0) == 0) {
                std::vector<std::string> f = split_on(line, '\t');
                for (std::size_t i = 9; i < f.size(); ++i) hap_order.push_back(f[i]);
            }
            continue;
        }
        std::vector<std::string> f = split_on(line, '\t');
        if (f.size() < 10) continue;
        const std::string& id = f[2];
        const std::string& info = f[7];
        const std::string svtype = info_value(info, "SVTYPE");
        const std::string bubble = info_value(info, "BUBBLE_ID");
        const std::string nodes = info_value(info, "EVENT_NODES");
        const std::string gene = info_value(info, "GENES");
        const std::string af = info_value(info, "AF");
        const std::string an = info_value(info, "AN");
        const std::string nalleles_s = info_value(info, "NALLELES");
        const int nalleles = nalleles_s.empty() ? 1 : std::max(1, std::atoi(nalleles_s.c_str()));
        const std::vector<std::string> fmt_keys = split_on(f[8], ':');
        int gt_i = -1, cn_i = -1;
        for (std::size_t i = 0; i < fmt_keys.size(); ++i) {
            if (fmt_keys[i] == "GT") gt_i = static_cast<int>(i);
            else if (fmt_keys[i] == "CN") cn_i = static_cast<int>(i);
        }
        const bool is_dup = (svtype == "DUP");
        const bool is_multi = nalleles > 1;
        // NALLELES counts REF plus the ALTs (variant_call.cpp writes alt_seqs.size() + 1, and the
        // header says so), while a dosage row is emitted per ALT. Using it directly produced one extra
        // row per multiallelic site, testing gt == NALLELES, which no sample can carry -- an all-zero
        // column for an allele that does not exist.
        const int n_rows = is_multi ? nalleles - 1 : 1;
        for (int a = 0; a < n_rows; ++a) {
            VariantBimbam v;
            v.id = is_multi ? (id + "_a" + std::to_string(a + 1)) : id;
            v.svtype = svtype.empty() ? "." : svtype;
            v.bubble = bubble; v.nodes = nodes.empty() ? "." : nodes;
            v.gene = gene.empty() ? "." : gene;
            v.af = af.empty() ? "." : af; v.an = an.empty() ? "." : an;
            for (std::size_t s = 0; s < hap_order.size() && 9 + s < f.size(); ++s) {
                const std::vector<std::string> sub = split_on(f[9 + s], ':');
                const std::string gt = (gt_i >= 0 && gt_i < static_cast<int>(sub.size()))
                                           ? sub[gt_i] : (sub.empty() ? std::string(".") : sub[0]);
                if (gt == "." || gt.empty()) continue;  // bubble not traversed -> NA (absent)
                double d;
                if (is_dup && cn_i >= 0 && cn_i < static_cast<int>(sub.size()) && sub[cn_i] != ".")
                    d = std::atof(sub[cn_i].c_str());
                else if (is_multi)
                    d = (std::atoi(gt.c_str()) == a + 1) ? 1.0 : 0.0;
                else
                    d = std::atof(gt.c_str());  // 0/1 presence
                v.dose[hap_order[s]] = d;
            }
            rows.push_back(std::move(v));
        }
    }

    auto write_annot = [&](GzipWriter& annot, const VariantBimbam& v) {
        std::string a = tsv_sanitize(v.id);
        a += "\tvariant\tdosage\t"; a += v.bubble; a += '\t'; a += v.nodes;
        a += '\t'; a += v.svtype; a += '\t'; a += v.gene; a += '\t'; a += v.af; a += '\t'; a += v.an; a += '\n';
        annot.write(a);
    };
    const std::string annot_header = "feature_id\tlayer\tencoding\tbubbles\tnodes\tsvtype\tgene\tAF\tAN\n";

    // Haplotype-level BIMBAM + sidecar + column order, in haplotype/variant/.
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
            const bool do_scale = options.scale_dosage && hi > lo;
            std::string g = tsv_sanitize(v.id); g += ", A, B";
            for (const std::string& s : hap_order) {
                g += ", ";
                const auto it = v.dose.find(s);
                if (it == v.dose.end()) g += "NA";
                else g += format_dosage(do_scale ? (it->second - lo) / (hi - lo) * 2.0 : it->second);
            }
            g += '\n';
            geno.write(g);
        }
        geno.close(); annot.close();
        summary.files_written += 3;
    }

    // Sample-level (diploid) BIMBAM: sum a sample's haplotype dosages; NA only if no haplotype traverses.
    // In sample/variant/, self-contained (matrix + feature_annot + column order).
    if (!options.samples_path.empty()) {
        const auto cosigt = read_cosigt_table(options.samples_path);
        const auto& path_to_samples = cosigt.first;
        const auto& sample_order = cosigt.second;
        const std::filesystem::path dir = substrate_dir(out_dir, "sample", "variant");
        GzipWriter sout((dir / "samples.txt.gz").string());
        for (const std::string& s : sample_order) { sout.write(s); sout.write("\n"); }
        sout.close();
        GzipWriter annot((dir / "feature_annot.variant.tsv.gz").string());
        annot.write(annot_header);
        GzipWriter geno((dir / "bimbam_variant.bimbam.gz").string());
        for (const VariantBimbam& v : rows) {
            write_annot(annot, v);
            std::unordered_map<std::string, double> samp;
            std::unordered_set<std::string> covered;
            for (const auto& [hap, d] : v.dose) {
                const auto it = path_to_samples.find(hap);
                if (it == path_to_samples.end()) continue;
                for (const std::string& s : it->second) { samp[s] += d; covered.insert(s); }
            }
            double lo = std::numeric_limits<double>::infinity(), hi = -lo;
            for (const std::string& s : covered) { lo = std::min(lo, samp[s]); hi = std::max(hi, samp[s]); }
            const bool do_scale = options.scale_dosage && hi > lo;
            std::string g = tsv_sanitize(v.id); g += ", A, B";
            for (const std::string& s : sample_order) {
                g += ", ";
                if (!covered.count(s)) g += "NA";
                else g += format_dosage(do_scale ? (samp[s] - lo) / (hi - lo) * 2.0 : samp[s]);
            }
            g += '\n';
            geno.write(g);
        }
        geno.close(); annot.close();
        summary.files_written += 3;
    }

    if (!options.quiet)
        std::cerr << "[describe] variant substrate: " << rows.size() << " variant rows from "
                  << options.variant_vcf_path << "\n";
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

    std::filesystem::create_directories(options.out_dir);

    ParseGfaOptions parse_options;
    parse_options.include_paths = true;
    parse_options.include_sequences = true;
    const Graph graph = parse_gfa(options.gfa_path, parse_options);
    if (graph.paths.empty()) {
        throw std::runtime_error("Input GFA has no P/W paths; describe requires paths");
    }

    const std::vector<Bubble> all_bubbles = read_bubbles_csv(options.bubbles_csv_in);
    std::unordered_set<std::size_t> bubble_filter;
    bubble_filter.reserve(options.bubble_ids.size() * 2 + 1);
    for (const std::size_t id : options.bubble_ids) {
        bubble_filter.insert(id);
    }

    // --variant-nodes: restrict to called-variation nodes. Parse <prefix>.variant_nodes.tsv
    // (variant_id, bubble_id, svtype, node_ids) into bubble_id -> {node ids}; only those
    // bubbles are processed and only their variant nodes contribute k-mers.
    std::unordered_map<std::size_t, std::unordered_set<std::string>> variant_nodes_by_bubble;
    const bool variant_mode = !options.variant_nodes_path.empty();
    if (variant_mode) {
        std::ifstream vn(options.variant_nodes_path);
        if (!vn) {
            throw std::runtime_error("Failed to open --variant-nodes file: " + options.variant_nodes_path);
        }
        std::string line;
        bool header = true;
        while (std::getline(vn, line)) {
            if (header) { header = false; continue; }  // skip "variant_id\tbubble_id\tsvtype\tnode_ids"
            if (line.empty()) continue;
            std::vector<std::string> f;
            std::size_t start = 0;
            for (std::size_t p = 0; p <= line.size(); ++p) {
                if (p == line.size() || line[p] == '\t') { f.push_back(line.substr(start, p - start)); start = p + 1; }
            }
            if (f.size() < 4) continue;
            const std::size_t bid = static_cast<std::size_t>(std::strtoull(f[1].c_str(), nullptr, 10));
            auto& set = variant_nodes_by_bubble[bid];
            std::size_t s = 0;
            const std::string& nodes = f[3];
            for (std::size_t p = 0; p <= nodes.size(); ++p) {
                if (p == nodes.size() || nodes[p] == ',') { if (p > s) set.insert(nodes.substr(s, p - s)); s = p + 1; }
            }
        }
    }

    std::vector<BubblePathIndex> path_indexes;
    path_indexes.reserve(graph.paths.size());
    for (const auto& path : graph.paths) {
        path_indexes.push_back(build_bubble_path_index(path));
    }

    const std::filesystem::path out_dir(options.out_dir);
    const std::string index_path = (out_dir / "describe.index.tsv").string();
    const std::string params_path = (out_dir / "describe.params.json").string();
    write_params_json(options, params_path);

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
    summary.files_written = 2; // index + params

    // bubbles to process (variant mode restricts to bubbles present in variant_nodes.tsv)
    std::vector<const Bubble*> to_process;
    for (const Bubble& bubble : all_bubbles) {
        if (!bubble_filter.empty() && bubble_filter.find(bubble.id) == bubble_filter.end()) continue;
        if (variant_mode && variant_nodes_by_bubble.find(bubble.id) == variant_nodes_by_bubble.end()) continue;
        to_process.push_back(&bubble);
    }

    KmerPool kmer_pool;
    // Per-haplotype node/edge dosage, pooled across bubbles for BIMBAM and --samples.
    GraphPool graph_pool;
    // Which pooled substrates to accumulate (the BIMBAM genotype export needs the cohort pools).
    const bool want_kmer_pool = options.emit_kmers && options.bimbam;
    const bool want_graph_pool =
        options.emit_graph && (options.bimbam || !options.samples_path.empty());
    const bool want_bimbam = options.bimbam && (options.emit_kmers || options.emit_graph);
    // For BIMBAM NA: which bubble(s) a feature is discriminative in, and which paths traverse each bubble.
    std::unordered_map<std::uint64_t, std::vector<std::size_t>> feature_bubbles_k;
    std::unordered_map<std::string, std::vector<std::size_t>> feature_bubbles_g;
    std::unordered_map<std::uint64_t, std::set<std::string>> feature_nodes_k;  // pooled k-mer code -> node set
    std::unordered_map<std::size_t, std::unordered_set<std::string>> bubble_traversers;
    cli::ProgressBar progress(options.quiet ? "" : "Describing bubbles", to_process.size());

    // Process bubbles in parallel: describe_one_bubble writes its own per-bubble files (independent
    // paths) and accumulates into THREAD-LOCAL pools here; the pools and the index rows are merged
    // back in bubble order on the main thread, so the output is byte-identical to a single-thread run.
    std::vector<BubbleDescribeResult> results(to_process.size());
    std::vector<KmerPool> task_kmer(want_kmer_pool ? to_process.size() : 0);
    std::vector<GraphPool> task_graph(want_graph_pool ? to_process.size() : 0);
    run_parallel(to_process.size(), options.threads, [&](std::size_t bi) {
        const Bubble& bubble = *to_process[bi];
        const std::unordered_set<std::string>* vnodes = nullptr;
        if (variant_mode) {
            const auto it = variant_nodes_by_bubble.find(bubble.id);
            if (it != variant_nodes_by_bubble.end()) vnodes = &it->second;
        }
        results[bi] = describe_one_bubble(options, graph, path_indexes, bubble, vnodes,
                                          want_kmer_pool ? &task_kmer[bi] : nullptr,
                                          want_graph_pool ? &task_graph[bi] : nullptr);
    });

    for (std::size_t bi = 0; bi < to_process.size(); ++bi) {
        const Bubble& bubble = *to_process[bi];
        progress.tick();
        if (want_bimbam && !results[bi].traversing_paths.empty()) {
            auto& tset = bubble_traversers[bubble.id];
            for (const std::string& p : results[bi].traversing_paths) tset.insert(p);
        }
        if (want_kmer_pool) {
            for (auto& [code, carriers] : task_kmer[bi]) {
                auto& dst = kmer_pool[code];
                for (const auto& [path, c] : carriers) dst[path] += c;
                if (want_bimbam) feature_bubbles_k[code].push_back(bubble.id);
            }
            if (want_bimbam) {
                for (auto& [code, ns] : results[bi].kmer_nodes) {
                    auto& dst = feature_nodes_k[code];
                    for (const std::string& n : ns) dst.insert(n);
                }
            }
        }
        if (want_graph_pool) {
            for (auto& [feature, carriers] : task_graph[bi]) {
                auto& dst = graph_pool[feature];
                for (const auto& [path, c] : carriers) dst[path] += c;
                if (want_bimbam) feature_bubbles_g[feature].push_back(bubble.id);
            }
        }
        const BubbleDescribeResult& result = results[bi];
        ++summary.bubbles_processed;
        if (result.paths > 0) {
            ++summary.bubbles_with_paths;
            summary.paths_written += result.paths;
            summary.features_written += result.features;
            summary.features_candidates += result.features_total;
            summary.jsonl_files_written += 1;
            summary.files_written += 2; // feature map + sparse JSONL
            if (result.matrix_written) {
                summary.matrix_files_written += 1;
                summary.files_written += 1;
            }
            summary.node_edge_features_written += result.node_features + result.edge_features;
            summary.node_edge_candidates += result.node_features_total + result.edge_features_total;
            summary.files_written += 1; // graph feature map
            if (result.graph_matrix_written) {
                summary.graph_matrix_files_written += 1;
                summary.files_written += 1;
            }
        }

        index_out
            << bubble.id << '\t'
            << result.status << '\t'
            << result.paths << '\t'
            << result.features_total << '\t'
            << result.features << '\t'
            << (result.features_total - result.features) << '\t'
            << result.feature_map_path << '\t'
            << result.matrix_path << '\t'
            << result.counts_jsonl_path << '\t'
            << (result.matrix_written ? "1" : "0") << '\t'
            << result.matrix_reason << '\t'
            << result.node_features_total << '\t'
            << result.node_features << '\t'
            << (result.node_features_total - result.node_features) << '\t'
            << result.edge_features_total << '\t'
            << result.edge_features << '\t'
            << (result.edge_features_total - result.edge_features) << '\t'
            << result.graph_feature_map_path << '\t'
            << result.graph_matrix_path << '\t'
            << (result.graph_matrix_written ? "1" : "0") << '\n';
    }

    // Pooled BIMBAM mean-genotype dosage, per substrate under haplotype/{kmers,graph}/: each folder gets
    // its own matrix, feature_annot, and column order (samples.txt.gz).
    if (want_bimbam) {
        std::vector<std::string> sample_order;
        sample_order.reserve(graph.paths.size());
        for (const auto& p : graph.paths) sample_order.push_back(p.name);
        const char* pool_header = "feature_id\tlayer\tencoding\tbubbles\tnodes\n";
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
                const auto nit = feature_nodes_k.find(code);  // aggregated node provenance (across bubbles)
                r.nodes = (nit != feature_nodes_k.end() && !nit->second.empty())
                    ? join_nodes(std::vector<std::string>(nit->second.begin(), nit->second.end())) : ".";
                const auto fb = feature_bubbles_k.find(code);
                if (fb != feature_bubbles_k.end()) r.bubbles = fb->second;
                for (const auto& [path, c] : carriers) r.carriers[path] = static_cast<double>(c);
                rows.push_back(std::move(r));
            }
            write_bimbam_rows((dir / "bimbam_kmers.bimbam.gz").string(), annot, "kmer",
                              rows, sample_order, bubble_traversers, options.scale_dosage);
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
                r.nodes = feature;  // node/edge id is its own provenance
                const auto fb = feature_bubbles_g.find(feature);
                if (fb != feature_bubbles_g.end()) r.bubbles = fb->second;
                for (const auto& [path, c] : carriers) r.carriers[path] = static_cast<double>(c);
                rows.push_back(std::move(r));
            }
            write_bimbam_rows((dir / "bimbam_graph.bimbam.gz").string(), annot, "graph",
                              rows, sample_order, bubble_traversers, options.scale_dosage);
            annot.close();
            summary.files_written += 3;
        }
    }

    // --samples (cosigt): aggregate any per-haplotype pool into per-SAMPLE dosage by summing over
    // the sample's assigned haplotype paths (a haplotype listed twice -> doubled, i.e. a homozygous
    // diploid genotype). Parsed once and reused for the k-mer and node/edge sample files.
    if (!options.samples_path.empty()) {
        const auto cosigt = read_cosigt_table(options.samples_path);
        const std::unordered_map<std::string, std::vector<std::string>>& path_to_samples = cosigt.first;
        const std::vector<std::string>& sample_order_s = cosigt.second;  // sample (column) order

        // Node/edge substrate (complementary graph-local layer): diploid summation, keyed by node id /
        // edge so it joins to call's variant_nodes.tsv. Built once here, reused for the sample-level BIMBAM.
        GraphPool graph_sample_pool;
        if (want_bimbam && options.emit_graph) {
            for (const auto& [feature, carriers] : graph_pool) {
                auto& out = graph_sample_pool[feature];
                for (const auto& [path, count] : carriers) {
                    const auto it = path_to_samples.find(path);
                    if (it == path_to_samples.end()) continue;
                    for (const std::string& s : it->second) out[s] += count;
                }
            }
        }

        // Sample-level BIMBAM (diploid, summed dosage) -- the realistic per-sample genotype input to
        // `panvar associate`. Coverage (for NA) is lifted from per-haplotype traversers to samples.
        if (want_bimbam) {
            std::unordered_map<std::size_t, std::unordered_set<std::string>> bubble_traversers_s;
            for (const auto& [b, paths] : bubble_traversers) {
                auto& dst = bubble_traversers_s[b];
                for (const std::string& path : paths) {
                    const auto it = path_to_samples.find(path);
                    if (it == path_to_samples.end()) continue;
                    for (const std::string& s : it->second) dst.insert(s);
                }
            }
            const char* pool_header = "feature_id\tlayer\tencoding\tbubbles\tnodes\n";
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
                KmerPool sample_pool;
                for (const auto& [code, carriers] : kmer_pool) {
                    auto& out = sample_pool[code];
                    for (const auto& [path, count] : carriers) {
                        const auto it = path_to_samples.find(path);
                        if (it == path_to_samples.end()) continue;
                        for (const std::string& s : it->second) out[s] += count;
                    }
                }
                std::vector<BimbamRow> rows; rows.reserve(sample_pool.size());
                for (const auto& [code, carriers] : sample_pool) {
                    BimbamRow r; r.id = decode_kmer(code, options.kmer_size);
                    const auto nit = feature_nodes_k.find(code);
                    r.nodes = (nit != feature_nodes_k.end() && !nit->second.empty())
                        ? join_nodes(std::vector<std::string>(nit->second.begin(), nit->second.end())) : ".";
                    const auto fb = feature_bubbles_k.find(code);
                    if (fb != feature_bubbles_k.end()) r.bubbles = fb->second;
                    for (const auto& [s, c] : carriers) r.carriers[s] = static_cast<double>(c);
                    rows.push_back(std::move(r));
                }
                write_bimbam_rows((dir / "bimbam_kmers.bimbam.gz").string(), annot, "kmer",
                                  rows, sample_order_s, bubble_traversers_s, options.scale_dosage);
                annot.close();
                summary.files_written += 3;
            }
            if (options.emit_graph) {
                const std::filesystem::path dir = substrate_dir(out_dir, "sample", "graph");
                write_samples_s(dir);
                GzipWriter annot((dir / "feature_annot.graph.tsv.gz").string());
                annot.write(pool_header);
                std::vector<BimbamRow> rows; rows.reserve(graph_sample_pool.size());
                for (const auto& [feature, carriers] : graph_sample_pool) {
                    BimbamRow r; r.id = feature; r.nodes = feature;
                    const auto fb = feature_bubbles_g.find(feature);
                    if (fb != feature_bubbles_g.end()) r.bubbles = fb->second;
                    for (const auto& [s, c] : carriers) r.carriers[s] = static_cast<double>(c);
                    rows.push_back(std::move(r));
                }
                write_bimbam_rows((dir / "bimbam_graph.bimbam.gz").string(), annot, "graph",
                                  rows, sample_order_s, bubble_traversers_s, options.scale_dosage);
                annot.close();
                summary.files_written += 3;
            }
            if (!options.quiet) {
                std::cerr << "[describe] sample-level BIMBAM: " << sample_order_s.size()
                          << " samples -> " << (out_dir / "sample").string() << "/{kmers,graph}/\n";
            }
        }
    }

    if (summary_out != nullptr) {
        *summary_out = summary;
    }
}

// Public entry: variant substrate on its own. Independent of describe_kmers_from_graph so --only-variant
// emits just the SV-call BIMBAM without loading the GFA. emit_variant_substrate holds the logic.
void describe_variant_from_vcf(const DescribeOptions& options, DescribeSummary* summary_out) {
    DescribeSummary local;
    DescribeSummary& summary = summary_out ? *summary_out : local;
    emit_variant_substrate(options, summary);
}

} // namespace panvar
