#include "panvar/describe.hpp"

#include "panvar/bubble_path.hpp"
#include "panvar/cli_utils.hpp"
#include "panvar/gfa.hpp"
#include "panvar/graph_utils.hpp"
#include "panvar/output.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
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

struct KmerOccurrence {
    std::uint64_t code = 0;
    std::size_t start = 0;
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

int base_bits(char c) {
    switch (c) {
        case 'A':
        case 'a':
            return 0;
        case 'C':
        case 'c':
            return 1;
        case 'G':
        case 'g':
            return 2;
        case 'T':
        case 't':
            return 3;
        default:
            return -1;
    }
}

std::string feature_mode_name(DescribeFeatureMode mode) {
    switch (mode) {
        case DescribeFeatureMode::AllKmers:
            return "all";
        case DescribeFeatureMode::Minimizer:
            return "minimizer";
        case DescribeFeatureMode::Syncmer:
            return "syncmer";
    }
    return "all";
}

std::size_t effective_syncmer_s(const DescribeOptions& options) {
    if (options.syncmer_s != 0) {
        return options.syncmer_s;
    }
    return std::max<std::size_t>(1, std::min<std::size_t>(11, (options.kmer_size + 2) / 3));
}

std::vector<KmerOccurrence> collect_canonical_kmer_occurrences(
    const std::string& sequence,
    std::size_t k) {

    std::vector<KmerOccurrence> out;
    if (k == 0 || k > 31 || sequence.size() < k) {
        return out;
    }

    const std::size_t possible = sequence.size() - k + 1;
    out.reserve(possible);

    const std::uint64_t mask = (1ULL << (2 * k)) - 1ULL;
    const std::size_t rc_shift = 2 * (k - 1);
    std::uint64_t fwd = 0;
    std::uint64_t rev = 0;
    std::size_t filled = 0;

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
            out.push_back(KmerOccurrence{std::min(fwd, rev), pos + 1 - k});
        }
    }

    return out;
}

bool is_closed_syncmer(std::uint64_t code, std::size_t k, std::size_t s) {
    if (s == 0 || s > k) {
        return false;
    }
    if (s == k) {
        return true;
    }

    const std::uint64_t mask = (1ULL << (2 * s)) - 1ULL;
    const std::size_t last_offset = k - s;
    std::uint64_t best = std::numeric_limits<std::uint64_t>::max();
    bool best_at_end = false;

    for (std::size_t offset = 0; offset <= last_offset; ++offset) {
        const std::size_t shift = 2 * (k - s - offset);
        const std::uint64_t sub = (code >> shift) & mask;
        const bool at_end = offset == 0 || offset == last_offset;
        if (sub < best) {
            best = sub;
            best_at_end = at_end;
        } else if (sub == best && at_end) {
            best_at_end = true;
        }
    }

    return best_at_end;
}

void select_minimizer_indices_for_run(
    const std::vector<KmerOccurrence>& occurrences,
    std::size_t begin,
    std::size_t end,
    std::size_t window,
    std::vector<bool>& selected) {

    const std::size_t len = end - begin;
    if (len == 0) {
        return;
    }
    if (len <= window) {
        std::size_t best = begin;
        for (std::size_t i = begin + 1; i < end; ++i) {
            if (occurrences[i].code < occurrences[best].code) {
                best = i;
            }
        }
        selected[best] = true;
        return;
    }

    std::vector<std::size_t> deque;
    deque.reserve(window + 1);
    std::size_t front = 0;

    for (std::size_t i = begin; i < end; ++i) {
        while (front < deque.size() && deque[front] + window <= i) {
            ++front;
        }
        while (deque.size() > front && occurrences[deque.back()].code > occurrences[i].code) {
            deque.pop_back();
        }
        deque.push_back(i);
        if (i + 1 - begin >= window) {
            selected[deque[front]] = true;
        }
        if (front > 64 && front * 2 > deque.size()) {
            deque.erase(deque.begin(), deque.begin() + static_cast<std::ptrdiff_t>(front));
            front = 0;
        }
    }
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

    if (options.feature_mode == DescribeFeatureMode::Syncmer) {
        const std::size_t s = effective_syncmer_s(options);
        out.reserve(occurrences.size() / 8 + 1);
        for (std::size_t i = 0; i < occurrences.size(); ++i) {
            if (is_closed_syncmer(occurrences[i].code, options.kmer_size, s)) {
                out.push_back(i);
            }
        }
        return out;
    }

    std::vector<bool> selected(occurrences.size(), false);
    const std::size_t window = std::max<std::size_t>(1, options.minimizer_window);
    std::size_t run_begin = 0;
    while (run_begin < occurrences.size()) {
        std::size_t run_end = run_begin + 1;
        while (run_end < occurrences.size() &&
               occurrences[run_end].start == occurrences[run_end - 1].start + 1) {
            ++run_end;
        }
        select_minimizer_indices_for_run(occurrences, run_begin, run_end, window, selected);
        run_begin = run_end;
    }

    out.reserve(occurrences.size() / window + 1);
    for (std::size_t i = 0; i < selected.size(); ++i) {
        if (selected[i]) {
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

std::string decode_kmer(std::uint64_t code, std::size_t k) {
    std::string out(k, 'A');
    for (std::size_t i = 0; i < k; ++i) {
        const std::size_t shift = 2 * (k - 1 - i);
        const std::uint64_t b = (code >> shift) & 0x3ULL;
        switch (b) {
            case 0:
                out[i] = 'A';
                break;
            case 1:
                out[i] = 'C';
                break;
            case 2:
                out[i] = 'G';
                break;
            default:
                out[i] = 'T';
                break;
        }
    }
    return out;
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

// Shared keep rule for both feature layers. A feature is kept when it is
// informative for association:
//  - copy-number features (count varies across carrying paths) are always kept;
//  - otherwise it must pass a symmetric minor-presence (MAF-style) cut:
//    min(present, absent) > min_paths.
// min_paths == 0 reproduces the legacy rule (drop only features present in every
// path with one identical count).
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

// --- Node and edge dosage features ----------------------------------------
// Per-bubble association substrate that parallels the k-mer features but is
// keyed by graph coordinates: how many times each path traverses each inside
// node (node dosage) and each oriented step-to-step transition (edge dosage).
// Shares node IDs with the k-mer feature map and the future graph-native calls.
//
// Node dosage is a DESCRIPTIVE traversal count, not a copy-number call: a count
// > 1 may be a tandem duplication or just the same node revisited elsewhere in
// the walk. Adjacency (the real tandem signal) is captured by the edge layer -
// a tandem block repeats an edge, scattered reuse does not. True CN/duplication
// (adjacent-block detection with motif/region size thresholds) is left to the
// graph-native variant-calling step.

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

GraphDosage count_graph_dosage(const std::vector<PathStep>& steps) {
    GraphDosage dosage;
    for (const PathStep& step : steps) {
        ++dosage.node_counts[step.node_id];
    }
    for (std::size_t i = 0; i + 1 < steps.size(); ++i) {
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
    std::map<std::string, GraphFeatureStat>& edge_stats) {

    for (const PathMeta& meta : paths) {
        const std::vector<PathStep> steps = canonical_steps_for_bubble_path(
            bubble, graph.paths[meta.path_index], path_indexes[meta.path_index]);
        if (steps.empty()) {
            continue;
        }
        const GraphDosage dosage = count_graph_dosage(steps);
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

    for (const PathMeta& meta : paths) {
        const std::vector<PathStep> steps = canonical_steps_for_bubble_path(
            bubble, graph.paths[meta.path_index], path_indexes[meta.path_index]);
        if (steps.empty()) {
            continue;
        }
        const GraphDosage dosage = count_graph_dosage(steps);

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

std::string path_sequence_for_bubble(
    const Graph& graph,
    const Bubble& bubble,
    const PathRecord& path,
    const BubblePathIndex& index,
    bool* complete_out) {

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
    bool complete = false;
    std::string sequence = spell_path_steps_sequence(graph, steps, &complete);
    if (complete_out != nullptr) {
        *complete_out = complete;
    }
    return complete ? sequence : std::string{};
}

BubblePathSequence path_sequence_with_segments_for_bubble(
    const Graph& graph,
    const Bubble& bubble,
    const PathRecord& path,
    const BubblePathIndex& index) {

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

    out.sequence.reserve(total_len);
    out.segments.reserve(steps.size());
    for (const PathStep& step : steps) {
        const auto& node_seq = graph.nodes.at(step.node_id).sequence;
        const std::size_t start = out.sequence.size();
        if (step.reverse) {
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
    std::unordered_map<std::uint64_t, KmerStats>& stats) {

    std::vector<PathMeta> paths;
    paths.reserve(graph.paths.size());
    std::unordered_map<std::uint64_t, std::uint32_t> counts;

    for (std::size_t i = 0; i < graph.paths.size(); ++i) {
        bool complete = false;
        const std::string sequence = path_sequence_for_bubble(graph, bubble, graph.paths[i], path_indexes[i], &complete);
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
    std::vector<std::unordered_set<std::string>>& feature_nodes) {

    GzipWriter out(path);
    std::unordered_map<std::uint64_t, KmerPathCount> counts;

    for (const PathMeta& meta : paths) {
        const BubblePathSequence path_sequence = path_sequence_with_segments_for_bubble(
            graph,
            bubble,
            graph.paths[meta.path_index],
            path_indexes[meta.path_index]);
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
    const DescribeOptions& options) {

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
            &complete);
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
        << "  \"minimizer_window\": " << options.minimizer_window << ",\n"
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

BubbleDescribeResult describe_one_bubble(
    const DescribeOptions& options,
    const Graph& graph,
    const std::vector<BubblePathIndex>& path_indexes,
    const Bubble& bubble) {

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
        stats);

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
        feature_nodes);
    write_feature_map(result.feature_map_path, features, stats, feature_nodes, options.kmer_size, paths.size());

    const bool wide_allowed_by_cap =
        options.max_wide_features == 0 || features.size() <= options.max_wide_features;
    if (options.write_wide_matrix && (options.force_wide_matrix || wide_allowed_by_cap)) {
        write_wide_matrix(result.matrix_path, graph, bubble, path_indexes, paths, features, options);
        result.matrix_written = true;
        result.matrix_reason = "written";
    } else if (!options.write_wide_matrix) {
        result.matrix_path = ".";
        result.matrix_reason = "skipped:disabled";
    } else {
        result.matrix_path = ".";
        result.matrix_reason = "skipped:feature-cap";
    }

    // Node + edge dosage features (association substrate), keyed by the same
    // graph coordinates as the k-mer node provenance.
    std::map<std::string, GraphFeatureStat> node_stats;
    std::map<std::string, GraphFeatureStat> edge_stats;
    accumulate_graph_stats(graph, bubble, path_indexes, paths, node_stats, edge_stats);
    const std::vector<std::string> node_features =
        select_discriminative_graph_features(node_stats, paths.size(), options.min_feature_paths);
    const std::vector<std::string> edge_features =
        select_discriminative_graph_features(edge_stats, paths.size(), options.min_feature_paths);
    result.node_features = node_features.size();
    result.node_features_total = node_stats.size();
    result.edge_features = edge_features.size();
    result.edge_features_total = edge_stats.size();
    write_graph_feature_map(
        result.graph_feature_map_path, node_features, edge_features, node_stats, edge_stats, paths.size());
    if (options.write_wide_matrix) {
        write_graph_matrix(
            result.graph_matrix_path, graph, bubble, path_indexes, paths, node_features, edge_features);
        result.graph_matrix_written = true;
    } else {
        result.graph_matrix_path = ".";
    }

    return result;
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
    if (options.feature_mode == DescribeFeatureMode::Minimizer && options.minimizer_window == 0) {
        throw std::runtime_error("Describe minimizer window must be > 0");
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

    for (const Bubble& bubble : all_bubbles) {
        if (!bubble_filter.empty() && bubble_filter.find(bubble.id) == bubble_filter.end()) {
            continue;
        }
        if (!options.quiet) {
            std::cerr << "[describe] bubble " << bubble.id << "\n";
        }
        BubbleDescribeResult result = describe_one_bubble(options, graph, path_indexes, bubble);
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

    if (summary_out != nullptr) {
        *summary_out = summary;
    }
}

} // namespace panvar
