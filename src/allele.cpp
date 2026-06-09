#include "panvar/allele.hpp"

#include "panvar/bubble_path.hpp"
#include <algorithm>
#include <atomic>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <queue>
#include <regex>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include <zlib.h>

#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wc99-extensions"
#endif
#include <minimap.h>
#if defined(__clang__)
#pragma clang diagnostic pop
#endif

#include "panvar/bubbles.hpp"
#include "panvar/graph_utils.hpp"
#include "panvar/output.hpp"
#include "call_internal.hpp"

namespace panvar {
namespace {

struct OdgiPathEntry {
    std::string path_name;
    std::size_t cluster_id = 0;
    std::size_t interval_start = 0;
    std::size_t interval_end = 0;
    bool source_to_sink = true;
};

struct ReferenceWindow {
    std::size_t range_start_bp = 0;
    std::size_t range_end_bp = 0;
};

struct DistanceMatrices {
    std::vector<std::vector<int>> abs;
    std::vector<std::vector<double>> norm;
};

struct SimilarityReportPaths {
    std::string bubble_dir;
    std::string alleles_tsv;
    std::string matrix_norm_tsv;
    std::string matrix_abs_tsv;
    std::string stats_tsv;
    std::string cluster_signatures_dir;
};

constexpr const char* kClusterPalette[] = {
    "#E41A1C",
    "#377EB8",
    "#4DAF4A",
    "#FF7F00",
    "#A65628",
    "#F781BF",
    "#999999",
    "#66C2A5",
    "#FC8D62",
    "#8DA0CB",
    "#E78AC3",
    "#A6D854",
};

constexpr std::size_t kClusterPaletteSize = sizeof(kClusterPalette) / sizeof(kClusterPalette[0]);
constexpr std::size_t kSequenceSketchKmer = 15;
constexpr std::size_t kSequenceSketchSize = 1024;
constexpr std::size_t kWalkSketchShingle = 3;
constexpr std::size_t kWalkSketchSize = 512;
// Auto switch heuristic for very large sequence bubbles in --distance-mode auto:
// if pair_count * max_token_len exceeds this threshold, prefer approximate
// sequence distances and greedy threshold clustering.
constexpr long double kAutoSequenceFastWorkThreshold = 2.5e8L;

double sketch_reject_margin(double min_similarity) {
    // Lower estimated identities than this margin below threshold are safely
    // considered dissimilar without DP.
    // Keep margins wider for permissive thresholds and tighter for strict ones.
    const double margin = 0.015 + 0.10 * (1.0 - min_similarity);
    return std::clamp(margin, 0.015, 0.03);
}

double sketch_accept_margin(double min_similarity) {
    // Require estimated identity to exceed threshold by a small margin before
    // accepting similarity without DP.
    // Keep this margin small enough that high thresholds (e.g. 0.95) still
    // have a practical fast-accept region.
    const double margin = 0.002 + 0.05 * (1.0 - min_similarity);
    return std::clamp(margin, 0.002, 0.01);
}

int bounded_levenshtein_distance(
    const std::string& a,
    const std::string& b,
    int max_dist);

int bounded_weighted_levenshtein_distance(
    const std::vector<std::uint64_t>& a,
    const std::vector<int>& a_weights,
    const std::vector<std::uint64_t>& b,
    const std::vector<int>& b_weights,
    int max_dist);

std::string csv_escape(const std::string& v) {
    if (v.find_first_of(",\"\n\r") == std::string::npos) {
        return v;
    }
    std::string out;
    out.reserve(v.size() + 2);
    out.push_back('"');
    for (const char c : v) {
        if (c == '"') {
            out.push_back('"');
            out.push_back('"');
        } else {
            out.push_back(c);
        }
    }
    out.push_back('"');
    return out;
}

std::string join_size_t(const std::vector<std::size_t>& values, char delim) {
    if (values.empty()) {
        return "";
    }
    std::string out;
    out.reserve(values.size() * 4);
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (i > 0) {
            out.push_back(delim);
        }
        out += std::to_string(values[i]);
    }
    return out;
}

std::string shell_quote(const std::string& in) {
    std::string out;
    out.reserve(in.size() + 8);
    out.push_back('\'');
    for (const char c : in) {
        if (c == '\'') {
            out += "'\\''";
        } else {
            out.push_back(c);
        }
    }
    out.push_back('\'');
    return out;
}

std::optional<ReferenceWindow> compute_reference_window(
    const PathRecord& reference_path,
    std::size_t interval_start,
    std::size_t interval_end,
    const std::unordered_map<std::string, Node>& nodes,
    std::size_t flank_nodes) {

    if (reference_path.steps.empty()) {
        return std::nullopt;
    }

    const std::size_t n = reference_path.steps.size();
    const std::size_t start = std::min(interval_start, interval_end);
    const std::size_t end = std::max(interval_start, interval_end);
    if (start >= n || end >= n) {
        return std::nullopt;
    }

    const std::size_t lo = (start > flank_nodes) ? (start - flank_nodes) : 0;
    const std::size_t hi = std::min(n - 1, end + flank_nodes);
    const auto pref = path_prefix_bp(reference_path, nodes);

    ReferenceWindow out;
    out.range_start_bp = pref[lo];
    out.range_end_bp = pref[hi + 1];
    return out;
}

std::unordered_map<std::size_t, std::string> cluster_palette(const std::set<std::size_t>& cluster_ids) {
    std::unordered_map<std::size_t, std::string> out;
    out.reserve(cluster_ids.size() * 2);
    std::size_t i = 0;
    for (const std::size_t cluster_id : cluster_ids) {
        out[cluster_id] = kClusterPalette[i % kClusterPaletteSize];
        ++i;
    }
    return out;
}

std::string html_escape(const std::string& in) {
    std::string out;
    out.reserve(in.size() + 16);
    for (const char c : in) {
        switch (c) {
            case '&':
                out += "&amp;";
                break;
            case '<':
                out += "&lt;";
                break;
            case '>':
                out += "&gt;";
                break;
            case '"':
                out += "&quot;";
                break;
            case '\'':
                out += "&#39;";
                break;
            default:
                out.push_back(c);
                break;
        }
    }
    return out;
}

double safe_div(double num, double den) {
    return den == 0.0 ? 0.0 : (num / den);
}

std::uint64_t hash_step_token(const PathStep& step) {
    // 64-bit FNV-1a over node id + orientation.
    std::uint64_t h = 1469598103934665603ULL;
    for (const unsigned char c : step.node_id) {
        h ^= static_cast<std::uint64_t>(c);
        h *= 1099511628211ULL;
    }
    h ^= (step.reverse ? 0xF0ULL : 0x0FULL);
    h *= 1099511628211ULL;
    return h;
}

std::uint64_t splitmix64(std::uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

int dna_code(char c) {
    switch (std::toupper(static_cast<unsigned char>(c))) {
        case 'A':
            return 0;
        case 'C':
            return 1;
        case 'G':
            return 2;
        case 'T':
            return 3;
        default:
            return -1;
    }
}

std::vector<std::uint64_t> build_sequence_minhash_sketch(
    const std::string& sequence,
    std::size_t kmer_size,
    std::size_t sketch_size) {

    if (kmer_size == 0 || sketch_size == 0 || kmer_size > 31 || sequence.size() < kmer_size) {
        return {};
    }

    std::priority_queue<std::uint64_t> top_hashes;
    std::unordered_set<std::uint64_t> seen;
    seen.reserve(std::min<std::size_t>(sequence.size(), 200000));

    const std::uint64_t mask = (kmer_size == 32)
                                   ? std::numeric_limits<std::uint64_t>::max()
                                   : ((1ULL << (2 * kmer_size)) - 1ULL);
    std::uint64_t forward = 0;
    std::uint64_t reverse = 0;
    std::size_t valid_len = 0;

    for (const char c : sequence) {
        const int base = dna_code(c);
        if (base < 0) {
            forward = 0;
            reverse = 0;
            valid_len = 0;
            continue;
        }

        forward = ((forward << 2) | static_cast<std::uint64_t>(base)) & mask;
        reverse = (reverse >> 2) |
                  (static_cast<std::uint64_t>(3 - base) << (2 * (kmer_size - 1)));
        ++valid_len;
        if (valid_len < kmer_size) {
            continue;
        }

        const std::uint64_t canonical = std::min(forward, reverse);
        const std::uint64_t h = splitmix64(canonical);
        if (!seen.insert(h).second) {
            continue;
        }
        if (top_hashes.size() < sketch_size) {
            top_hashes.push(h);
        } else if (h < top_hashes.top()) {
            top_hashes.pop();
            top_hashes.push(h);
        }
    }

    std::vector<std::uint64_t> sketch(top_hashes.size());
    for (std::size_t i = sketch.size(); i > 0; --i) {
        sketch[i - 1] = top_hashes.top();
        top_hashes.pop();
    }
    std::sort(sketch.begin(), sketch.end());
    return sketch;
}

std::vector<std::uint64_t> build_walk_minhash_sketch(
    const std::vector<std::uint64_t>& tokens,
    std::size_t shingle_size,
    std::size_t sketch_size) {

    if (shingle_size == 0 || sketch_size == 0 || tokens.size() < shingle_size) {
        return {};
    }

    std::priority_queue<std::uint64_t> top_hashes;
    std::unordered_set<std::uint64_t> seen;
    seen.reserve(std::min<std::size_t>(tokens.size(), 200000));

    for (std::size_t i = 0; i + shingle_size <= tokens.size(); ++i) {
        std::uint64_t h = 1469598103934665603ULL;
        for (std::size_t k = 0; k < shingle_size; ++k) {
            std::uint64_t v = splitmix64(tokens[i + k] + (0x9e3779b97f4a7c15ULL * (k + 1)));
            h ^= v;
            h *= 1099511628211ULL;
        }
        h = splitmix64(h);
        if (!seen.insert(h).second) {
            continue;
        }
        if (top_hashes.size() < sketch_size) {
            top_hashes.push(h);
        } else if (h < top_hashes.top()) {
            top_hashes.pop();
            top_hashes.push(h);
        }
    }

    std::vector<std::uint64_t> sketch(top_hashes.size());
    for (std::size_t i = sketch.size(); i > 0; --i) {
        sketch[i - 1] = top_hashes.top();
        top_hashes.pop();
    }
    std::sort(sketch.begin(), sketch.end());
    return sketch;
}

double sketch_jaccard(
    const std::vector<std::uint64_t>& a,
    const std::vector<std::uint64_t>& b) {

    if (a.empty() || b.empty()) {
        return 0.0;
    }
    std::size_t i = 0;
    std::size_t j = 0;
    std::size_t inter = 0;
    while (i < a.size() && j < b.size()) {
        if (a[i] == b[j]) {
            ++inter;
            ++i;
            ++j;
        } else if (a[i] < b[j]) {
            ++i;
        } else {
            ++j;
        }
    }
    const std::size_t uni = a.size() + b.size() - inter;
    if (uni == 0) {
        return 0.0;
    }
    return static_cast<double>(inter) / static_cast<double>(uni);
}

double estimate_identity_from_jaccard(double jaccard) {
    // Approximate identity from k-mer Jaccard relation:
    // identity ~= 2J / (1 + J)
    if (jaccard <= 0.0) {
        return 0.0;
    }
    const double id = (2.0 * jaccard) / (1.0 + jaccard);
    return std::clamp(id, 0.0, 1.0);
}

std::vector<std::uint64_t> build_walk_tokens(const std::vector<PathStep>& steps) {
    std::vector<std::uint64_t> tokens;
    tokens.reserve(steps.size());
    for (const auto& step : steps) {
        tokens.push_back(hash_step_token(step));
    }
    return tokens;
}

std::vector<int> build_walk_step_weights(
    const std::vector<PathStep>& steps,
    const std::unordered_map<std::string, Node>& nodes) {

    std::vector<int> weights;
    weights.reserve(steps.size());
    for (const auto& step : steps) {
        std::size_t len = 1;
        const auto node_it = nodes.find(step.node_id);
        if (node_it != nodes.end()) {
            len = std::max<std::size_t>(1, node_it->second.sequence.size());
        }
        const int w = static_cast<int>(std::min<std::size_t>(
            len,
            static_cast<std::size_t>(std::numeric_limits<int>::max() / 4)));
        weights.push_back(std::max(1, w));
    }
    return weights;
}

std::size_t sum_step_weights(const std::vector<int>& weights) {
    std::size_t sum = 0;
    for (const int w : weights) {
        sum += static_cast<std::size_t>(std::max(1, w));
    }
    return sum;
}

std::size_t token_length_for_distance(const UniqueAllele& allele) {
    if (allele.uses_sequence_similarity) {
        return allele.compare_token.size();
    }
    if (allele.compare_steps_weight_sum > 0) {
        return allele.compare_steps_weight_sum;
    }
    return allele.compare_steps.size();
}

bool is_sequence_cluster_mode(ClusterMode mode) {
    return mode == ClusterMode::Sequence || mode == ClusterMode::SequenceFast;
}

double elapsed_seconds(const std::chrono::steady_clock::time_point& start) {
    const auto now = std::chrono::steady_clock::now();
    const auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(now - start).count();
    return static_cast<double>(ms) / 1000.0;
}

int allowed_max_edits(double min_similarity, std::size_t max_len) {
    if (max_len == 0) {
        return 0;
    }
    const double raw = (1.0 - min_similarity) * static_cast<double>(max_len);
    // Keep threshold deterministic around values like 0.1 * 50 = 5, where
    // binary floating-point can produce 4.999999999...
    const double eps = 1e-9;
    int allowed = static_cast<int>(std::floor(raw + eps));
    if (allowed < 0) {
        allowed = 0;
    }
    const int max_len_i = static_cast<int>(std::min<std::size_t>(
        max_len,
        static_cast<std::size_t>(std::numeric_limits<int>::max())));
    if (allowed > max_len_i) {
        allowed = max_len_i;
    }
    return allowed;
}

int exact_levenshtein_distance(const std::string& a, const std::string& b) {
    const std::size_t max_len = std::max(a.size(), b.size());
    if (max_len > static_cast<std::size_t>(std::numeric_limits<int>::max() - 1)) {
        throw std::runtime_error("Allele token too long for exact distance computation");
    }
    return bounded_levenshtein_distance(a, b, static_cast<int>(max_len));
}

int exact_levenshtein_distance(
    const std::vector<std::uint64_t>& a,
    const std::vector<int>& a_weights,
    const std::vector<std::uint64_t>& b,
    const std::vector<int>& b_weights) {

    if (a.size() != a_weights.size() || b.size() != b_weights.size()) {
        throw std::runtime_error("Weighted walk token/weight length mismatch");
    }
    std::size_t sum_a = 0;
    for (const int w : a_weights) {
        sum_a += static_cast<std::size_t>(std::max(1, w));
    }
    std::size_t sum_b = 0;
    for (const int w : b_weights) {
        sum_b += static_cast<std::size_t>(std::max(1, w));
    }
    const std::size_t upper = sum_a + sum_b;
    if (upper > static_cast<std::size_t>(std::numeric_limits<int>::max() - 2)) {
        throw std::runtime_error("Weighted walk token lists too long for exact distance computation");
    }
    return bounded_weighted_levenshtein_distance(
        a,
        a_weights,
        b,
        b_weights,
        static_cast<int>(upper));
}

std::size_t resolve_thread_count(std::size_t requested) {
    if (requested > 0) {
        return requested;
    }
    const unsigned int hw = std::thread::hardware_concurrency();
    if (hw == 0) {
        return 1;
    }
    return static_cast<std::size_t>(hw);
}

DistanceMatrices build_distance_matrices(
    const std::vector<UniqueAllele>& unique_alleles,
    double min_similarity,
    bool fast_distance,
    std::size_t exact_distance_max_bp,
    std::size_t requested_threads,
    bool show_progress,
    std::size_t bubble_id,
    bool approximate_walk_for_large_bubble,
    bool approximate_sequence_for_fast_mode) {

    const std::size_t n = unique_alleles.size();
    DistanceMatrices mats;
    mats.abs.assign(n, std::vector<int>(n, 0));
    mats.norm.assign(n, std::vector<double>(n, 0.0));
    if (n <= 1) {
        return mats;
    }

    std::vector<std::vector<std::uint64_t>> sequence_sketches;
    std::vector<std::vector<std::uint64_t>> walk_sketches;
    sequence_sketches.resize(n);
    walk_sketches.resize(n);
    if (fast_distance) {
        for (std::size_t i = 0; i < n; ++i) {
            const auto& allele = unique_alleles[i];
            if (allele.uses_sequence_similarity) {
                if (!allele.compare_token.empty()) {
                    sequence_sketches[i] =
                        build_sequence_minhash_sketch(allele.compare_token, kSequenceSketchKmer, kSequenceSketchSize);
                }
            } else {
                if (!allele.compare_steps.empty()) {
                    walk_sketches[i] =
                        build_walk_minhash_sketch(allele.compare_steps, kWalkSketchShingle, kWalkSketchSize);
                }
            }
        }
    }

    const std::size_t threads = std::max<std::size_t>(1, std::min<std::size_t>(resolve_thread_count(requested_threads), n));
    const std::uint64_t total_pairs = static_cast<std::uint64_t>(n) * static_cast<std::uint64_t>(n - 1) / 2ULL;
    std::atomic<std::size_t> next_i{0};
    std::atomic<std::uint64_t> done_pairs{0};
    std::atomic<bool> progress_done{false};
    std::thread progress_thread;

    if (show_progress && total_pairs >= 50000) {
        progress_thread = std::thread([&]() {
            using clock = std::chrono::steady_clock;
            auto last = clock::now();
            std::uint64_t last_done = 0;
            while (!progress_done.load(std::memory_order_acquire)) {
                std::this_thread::sleep_for(std::chrono::seconds(2));
                const auto now = clock::now();
                const std::uint64_t done = done_pairs.load(std::memory_order_relaxed);
                const double pct = (total_pairs == 0)
                                       ? 100.0
                                       : (100.0 * static_cast<double>(done) / static_cast<double>(total_pairs));
                const double dt = std::chrono::duration_cast<std::chrono::milliseconds>(now - last).count() / 1000.0;
                const double rate = (dt > 0.0)
                                        ? (static_cast<double>(done - last_done) / dt)
                                        : 0.0;
                std::cerr
                    << "[allele] bubble " << bubble_id
                    << " distance progress: " << done << "/" << total_pairs
                    << " (" << std::fixed << std::setprecision(1) << pct << "%), "
                    << std::setprecision(0) << rate << " pairs/s\n";
                last = now;
                last_done = done;
                if (done >= total_pairs) {
                    break;
                }
            }
        });
    }

    auto compute_pair = [&](std::size_t i, std::size_t j) {
        const auto& lhs = unique_alleles[i];
        const auto& rhs = unique_alleles[j];
        const bool use_sequence = lhs.uses_sequence_similarity && rhs.uses_sequence_similarity;
        const std::size_t lhs_len = token_length_for_distance(lhs);
        const std::size_t rhs_len = token_length_for_distance(rhs);
        const std::size_t max_len = std::max(lhs_len, rhs_len);

        int d_abs = 0;
        bool decided_without_dp = false;
        const int allowed = allowed_max_edits(min_similarity, max_len);
        const int length_gap = static_cast<int>(
            (lhs_len > rhs_len) ? (lhs_len - rhs_len) : (rhs_len - lhs_len));
        if (length_gap > allowed) {
            d_abs = allowed + 1;
            decided_without_dp = true;
        }

        if (use_sequence && fast_distance) {
            const auto& s1 = sequence_sketches[i];
            const auto& s2 = sequence_sketches[j];
            if (approximate_sequence_for_fast_mode) {
                if (!s1.empty() && !s2.empty()) {
                    const double jacc = sketch_jaccard(s1, s2);
                    const double est_id = estimate_identity_from_jaccard(jacc);
                    const int max_len_i = static_cast<int>(std::min<std::size_t>(
                        max_len,
                        static_cast<std::size_t>(std::numeric_limits<int>::max())));
                    int est_d = static_cast<int>(std::llround(
                        (1.0 - est_id) * static_cast<double>(max_len)));
                    est_d = std::clamp(est_d, 0, max_len_i);
                    d_abs = est_d;
                } else {
                    // Fallback when sequence sketches are unavailable.
                    d_abs = (length_gap > allowed) ? (allowed + 1) : allowed;
                }
                decided_without_dp = true;
            } else if (!s1.empty() && !s2.empty()) {
                const double jacc = sketch_jaccard(s1, s2);
                const double est_id = estimate_identity_from_jaccard(jacc);
                const double reject_margin = sketch_reject_margin(min_similarity);
                const double accept_margin = sketch_accept_margin(min_similarity);
                if (est_id + reject_margin < min_similarity) {
                    d_abs = allowed + 1;
                    decided_without_dp = true;
                } else if (est_id > min_similarity + accept_margin) {
                    const int est_d = static_cast<int>(std::llround((1.0 - est_id) * static_cast<double>(max_len)));
                    d_abs = std::clamp(est_d, 0, allowed);
                    decided_without_dp = true;
                }
            }
        } else if (!use_sequence && fast_distance) {
            const auto& s1 = walk_sketches[i];
            const auto& s2 = walk_sketches[j];
            if (!s1.empty() && !s2.empty()) {
                const double jacc = sketch_jaccard(s1, s2);
                const double est_id = estimate_identity_from_jaccard(jacc);
                const double reject_margin = sketch_reject_margin(min_similarity);
                if (est_id + reject_margin < min_similarity) {
                    d_abs = allowed + 1;
                    decided_without_dp = true;
                } else if (approximate_walk_for_large_bubble) {
                    // For very large bubbles in walk mode, avoid DP and use
                    // sketch-estimated distance directly.
                    const int max_len_i = static_cast<int>(std::min<std::size_t>(
                        max_len,
                        static_cast<std::size_t>(std::numeric_limits<int>::max())));
                    const int est_d = static_cast<int>(std::llround((1.0 - est_id) * static_cast<double>(max_len)));
                    d_abs = std::clamp(est_d, 0, max_len_i);
                    if (d_abs > allowed) {
                        d_abs = allowed + 1;
                    }
                    decided_without_dp = true;
                }
            } else if (approximate_walk_for_large_bubble) {
                // Fallback when sketches are unavailable (very short walks):
                // use a conservative threshold-oriented estimate.
                d_abs = (length_gap > allowed) ? (allowed + 1) : allowed;
                decided_without_dp = true;
            }
        }

        if (!decided_without_dp) {
            if (!fast_distance) {
                if (use_sequence) {
                    d_abs = exact_levenshtein_distance(lhs.compare_token, rhs.compare_token);
                } else {
                    d_abs = exact_levenshtein_distance(
                        lhs.compare_steps,
                        lhs.compare_step_weights,
                        rhs.compare_steps,
                        rhs.compare_step_weights);
                }
            } else {
                // In auto mode, use threshold-bounded DP for all uncertain pairs.
                // This is exact when distance <= allowed, and early-outs when > allowed.
                if (use_sequence) {
                    d_abs = bounded_levenshtein_distance(lhs.compare_token, rhs.compare_token, allowed);
                } else {
                    d_abs = bounded_weighted_levenshtein_distance(
                        lhs.compare_steps,
                        lhs.compare_step_weights,
                        rhs.compare_steps,
                        rhs.compare_step_weights,
                        allowed);
                }
                if (d_abs > allowed) {
                    d_abs = allowed + 1;
                }
            }
        }

        const double denom = static_cast<double>(std::max<std::size_t>(1, max_len));
        const double d_norm = std::min(1.0, static_cast<double>(d_abs) / denom);
        mats.abs[i][j] = d_abs;
        mats.abs[j][i] = d_abs;
        mats.norm[i][j] = d_norm;
        mats.norm[j][i] = d_norm;
        done_pairs.fetch_add(1, std::memory_order_relaxed);
    };

    (void)exact_distance_max_bp;

    auto worker = [&]() {
        while (true) {
            const std::size_t i = next_i.fetch_add(1, std::memory_order_relaxed);
            if (i >= n) {
                break;
            }
            for (std::size_t j = i + 1; j < n; ++j) {
                compute_pair(i, j);
            }
        }
    };

    if (threads == 1) {
        worker();
    } else {
        std::vector<std::thread> pool;
        pool.reserve(threads);
        for (std::size_t t = 0; t < threads; ++t) {
            pool.emplace_back(worker);
        }
        for (auto& th : pool) {
            th.join();
        }
    }

    progress_done.store(true, std::memory_order_release);
    if (progress_thread.joinable()) {
        progress_thread.join();
    }

    return mats;
}

void write_similarity_alleles_tsv(
    const std::string& output_path,
    const Bubble& bubble,
    const std::vector<UniqueAllele>& unique_alleles,
    const std::vector<std::size_t>& cluster_of_unique) {

    std::ofstream out(output_path);
    if (!out) {
        throw std::runtime_error("Failed to write similarity allele table: " + output_path);
    }

    out << "bubble_id\tallele_id\tcluster_id\tpath_support\tsequence_length\tmode\tsignature\n";
    for (std::size_t i = 0; i < unique_alleles.size(); ++i) {
        const auto& allele = unique_alleles[i];
        out << bubble.id << '\t'
            << allele.allele_id << '\t'
            << cluster_of_unique[i] << '\t'
            << allele.path_support << '\t'
            << allele.sequence_length << '\t'
            << (allele.uses_sequence_similarity ? "sequence" : "walk") << '\t'
            << allele.signature << '\n';
    }
}

void write_similarity_distance_matrix_tsv(
    const std::string& output_path,
    const std::vector<UniqueAllele>& unique_alleles,
    const std::vector<std::vector<double>>& dist) {

    std::ofstream out(output_path);
    if (!out) {
        throw std::runtime_error("Failed to write similarity distance matrix: " + output_path);
    }

    out << std::fixed << std::setprecision(6);
    out << "allele_id";
    for (const auto& allele : unique_alleles) {
        out << '\t' << allele.allele_id;
    }
    out << '\n';

    for (std::size_t i = 0; i < unique_alleles.size(); ++i) {
        out << unique_alleles[i].allele_id;
        for (std::size_t j = 0; j < unique_alleles.size(); ++j) {
            out << '\t' << dist[i][j];
        }
        out << '\n';
    }
}

void write_similarity_distance_matrix_abs_tsv(
    const std::string& output_path,
    const std::vector<UniqueAllele>& unique_alleles,
    const std::vector<std::vector<int>>& dist_abs) {

    std::ofstream out(output_path);
    if (!out) {
        throw std::runtime_error("Failed to write absolute distance matrix: " + output_path);
    }

    out << "allele_id";
    for (const auto& allele : unique_alleles) {
        out << '\t' << allele.allele_id;
    }
    out << '\n';

    for (std::size_t i = 0; i < unique_alleles.size(); ++i) {
        out << unique_alleles[i].allele_id;
        for (std::size_t j = 0; j < unique_alleles.size(); ++j) {
            out << '\t' << dist_abs[i][j];
        }
        out << '\n';
    }
}

std::vector<double> compute_silhouette_scores(
    const std::vector<std::size_t>& cluster_of_unique,
    const std::vector<std::vector<double>>& dist_norm) {

    const std::size_t n = cluster_of_unique.size();
    std::vector<double> scores(n, 0.0);
    if (n <= 1 || dist_norm.size() != n) {
        return scores;
    }

    std::unordered_map<std::size_t, std::vector<std::size_t>> members;
    members.reserve(n * 2);
    for (std::size_t i = 0; i < n; ++i) {
        members[cluster_of_unique[i]].push_back(i);
    }
    if (members.size() <= 1) {
        return scores;
    }

    for (std::size_t i = 0; i < n; ++i) {
        const auto& own = members.at(cluster_of_unique[i]);
        double a = 0.0;
        if (own.size() > 1) {
            for (const std::size_t j : own) {
                if (j != i) {
                    a += dist_norm[i][j];
                }
            }
            a /= static_cast<double>(own.size() - 1);
        }

        double b = std::numeric_limits<double>::infinity();
        for (const auto& [cluster_id, idxs] : members) {
            if (cluster_id == cluster_of_unique[i] || idxs.empty()) {
                continue;
            }
            double mean = 0.0;
            for (const std::size_t j : idxs) {
                mean += dist_norm[i][j];
            }
            mean /= static_cast<double>(idxs.size());
            b = std::min(b, mean);
        }
        if (!std::isfinite(b)) {
            scores[i] = 0.0;
            continue;
        }
        const double denom = std::max(a, b);
        scores[i] = (denom <= 0.0) ? 0.0 : ((b - a) / denom);
    }
    return scores;
}

void write_similarity_cluster_stats_tsv(
    const std::string& output_path,
    const std::vector<UniqueAllele>& unique_alleles,
    const std::vector<std::size_t>& cluster_of_unique,
    const std::vector<std::vector<double>>& dist_norm,
    const std::vector<double>& silhouette_scores,
    const std::vector<std::vector<int>>& dist_abs,
    double min_similarity) {

    std::unordered_map<std::size_t, std::vector<std::size_t>> members;
    for (std::size_t i = 0; i < cluster_of_unique.size(); ++i) {
        members[cluster_of_unique[i]].push_back(i);
    }

    std::vector<std::size_t> cluster_ids;
    cluster_ids.reserve(members.size());
    for (const auto& [cluster_id, _idxs] : members) {
        cluster_ids.push_back(cluster_id);
    }
    std::sort(cluster_ids.begin(), cluster_ids.end());

    std::ofstream out(output_path);
    if (!out) {
        throw std::runtime_error("Failed to write similarity cluster stats: " + output_path);
    }
    out << std::fixed << std::setprecision(6);
    out << "cluster_id\tmember_alleles\tpath_support_sum\t"
           "mean_silhouette\tpath_weighted_mean_silhouette\t"
           "mean_intra_distance_norm\tmean_intra_similarity\tmean_intra_edit_distance\t"
           "max_intra_distance_norm\tmax_intra_edit_distance\t"
           "mean_to_other_clusters_distance_norm\tmean_to_other_clusters_similarity\t"
           "mean_to_other_clusters_edit_distance\t"
           "min_to_other_clusters_distance_norm\tmax_to_other_clusters_distance_norm\t"
           "min_to_other_clusters_edit_distance\tmax_to_other_clusters_edit_distance\t"
           "nearest_other_cluster_id\tnearest_other_min_distance_norm\tnearest_other_min_similarity\t"
           "separation_margin_norm\tthreshold_margin_norm\tseparation_quality\n";

    for (const std::size_t cluster_id : cluster_ids) {
        const auto& idxs = members[cluster_id];
        double intra_sum_norm = 0.0;
        double intra_sum_abs = 0.0;
        std::size_t intra_pairs = 0;
        double intra_max_norm = 0.0;
        double intra_max_abs = 0.0;
        for (std::size_t a = 0; a < idxs.size(); ++a) {
            for (std::size_t b = a + 1; b < idxs.size(); ++b) {
                intra_sum_norm += dist_norm[idxs[a]][idxs[b]];
                intra_sum_abs += static_cast<double>(dist_abs[idxs[a]][idxs[b]]);
                intra_max_norm = std::max(intra_max_norm, dist_norm[idxs[a]][idxs[b]]);
                intra_max_abs = std::max(intra_max_abs, static_cast<double>(dist_abs[idxs[a]][idxs[b]]));
                ++intra_pairs;
            }
        }
        const double intra_mean_norm = safe_div(intra_sum_norm, static_cast<double>(intra_pairs));
        const double intra_mean_abs = safe_div(intra_sum_abs, static_cast<double>(intra_pairs));

        double between_cluster_sum_norm = 0.0;
        double between_cluster_sum_abs = 0.0;
        std::size_t between_cluster_count = 0;
        double between_cluster_min_norm = std::numeric_limits<double>::infinity();
        double between_cluster_max_norm = 0.0;
        double between_cluster_min_abs = std::numeric_limits<double>::infinity();
        double between_cluster_max_abs = 0.0;
        std::size_t nearest_other_cluster_id = 0;
        double nearest_other_min_norm = std::numeric_limits<double>::infinity();
        for (const std::size_t other_cluster : cluster_ids) {
            if (other_cluster == cluster_id) {
                continue;
            }
            const auto& other = members[other_cluster];
            double sum_norm = 0.0;
            double sum_abs = 0.0;
            std::size_t pairs = 0;
            double pair_min_norm = std::numeric_limits<double>::infinity();
            for (const std::size_t i : idxs) {
                for (const std::size_t j : other) {
                    sum_norm += dist_norm[i][j];
                    sum_abs += static_cast<double>(dist_abs[i][j]);
                    pair_min_norm = std::min(pair_min_norm, dist_norm[i][j]);
                    ++pairs;
                }
            }
            if (pairs == 0) {
                continue;
            }
            const double m_norm = sum_norm / static_cast<double>(pairs);
            const double m_abs = sum_abs / static_cast<double>(pairs);
            between_cluster_sum_norm += m_norm;
            between_cluster_sum_abs += m_abs;
            ++between_cluster_count;
            between_cluster_min_norm = std::min(between_cluster_min_norm, m_norm);
            between_cluster_max_norm = std::max(between_cluster_max_norm, m_norm);
            between_cluster_min_abs = std::min(between_cluster_min_abs, m_abs);
            between_cluster_max_abs = std::max(between_cluster_max_abs, m_abs);
            if (pair_min_norm < nearest_other_min_norm) {
                nearest_other_min_norm = pair_min_norm;
                nearest_other_cluster_id = other_cluster;
            }
        }
        const double between_cluster_mean_norm = safe_div(
            between_cluster_sum_norm,
            static_cast<double>(between_cluster_count));
        const double between_cluster_mean_abs = safe_div(
            between_cluster_sum_abs,
            static_cast<double>(between_cluster_count));
        if (!std::isfinite(between_cluster_min_norm)) {
            between_cluster_min_norm = 0.0;
        }
        if (!std::isfinite(between_cluster_min_abs)) {
            between_cluster_min_abs = 0.0;
        }
        if (!std::isfinite(nearest_other_min_norm)) {
            nearest_other_min_norm = 0.0;
            nearest_other_cluster_id = 0;
        }

        const double max_norm_dist = std::clamp(1.0 - min_similarity, 0.0, 1.0);
        const double separation_margin_norm = nearest_other_min_norm - intra_max_norm;
        const double threshold_margin_norm = nearest_other_min_norm - max_norm_dist;
        std::string separation_quality = "single_cluster";
        if (cluster_ids.size() > 1) {
            if (threshold_margin_norm < 0.0 || separation_margin_norm < 0.0) {
                separation_quality = "poor";
            } else if (threshold_margin_norm < 0.02 || separation_margin_norm < 0.02) {
                separation_quality = "borderline";
            } else {
                separation_quality = "good";
            }
        }

        std::size_t path_support_sum = 0;
        double silhouette_sum = 0.0;
        double weighted_silhouette_sum = 0.0;
        for (const std::size_t idx : idxs) {
            path_support_sum += unique_alleles[idx].path_support;
            silhouette_sum += silhouette_scores[idx];
            weighted_silhouette_sum +=
                silhouette_scores[idx] * static_cast<double>(unique_alleles[idx].path_support);
        }

        out << cluster_id << '\t'
            << idxs.size() << '\t'
            << path_support_sum << '\t'
            << safe_div(silhouette_sum, static_cast<double>(idxs.size())) << '\t'
            << safe_div(weighted_silhouette_sum, static_cast<double>(path_support_sum)) << '\t'
            << intra_mean_norm << '\t'
            << (1.0 - intra_mean_norm) << '\t'
            << intra_mean_abs << '\t'
            << intra_max_norm << '\t'
            << intra_max_abs << '\t'
            << between_cluster_mean_norm << '\t'
            << (1.0 - between_cluster_mean_norm) << '\t'
            << between_cluster_mean_abs << '\t'
            << between_cluster_min_norm << '\t'
            << between_cluster_max_norm << '\t'
            << between_cluster_min_abs << '\t'
            << between_cluster_max_abs << '\t'
            << nearest_other_cluster_id << '\t'
            << nearest_other_min_norm << '\t'
            << (1.0 - nearest_other_min_norm) << '\t'
            << separation_margin_norm << '\t'
            << threshold_margin_norm << '\t'
            << separation_quality << '\n';
    }
}

std::size_t node_length_bp(
    const std::unordered_map<std::string, Node>& nodes,
    const std::string& node_id) {

    const auto it = nodes.find(node_id);
    if (it == nodes.end() || it->second.sequence.empty()) {
        return 1;
    }
    return std::max<std::size_t>(1, it->second.sequence.size());
}

std::string allele_node_orientation_cell(
    const UniqueAllele& allele,
    const std::string& node_id) {

    std::string cell;
    for (const PathStep& step : allele.steps) {
        if (step.node_id == node_id) {
            cell.push_back(step.reverse ? '-' : '+');
        }
    }
    return cell.empty() ? std::string("|") : cell;
}

std::vector<std::string> cluster_signature_node_order(
    const Cluster& cluster,
    const std::vector<UniqueAllele>& unique_alleles) {

    std::vector<std::string> order;
    std::unordered_set<std::string> seen;

    auto add_nodes = [&](std::size_t allele_idx) {
        if (allele_idx >= unique_alleles.size()) {
            return;
        }
        for (const PathStep& step : unique_alleles[allele_idx].steps) {
            if (seen.insert(step.node_id).second) {
                order.push_back(step.node_id);
            }
        }
    };

    add_nodes(cluster.representative_idx);
    for (const std::size_t idx : cluster.member_unique_idxs) {
        add_nodes(idx);
    }
    return order;
}

void write_similarity_cluster_signature_tables(
    const std::string& output_dir,
    const std::unordered_map<std::string, Node>& nodes,
    const std::vector<UniqueAllele>& unique_alleles,
    const std::vector<Cluster>& clusters) {

    std::filesystem::create_directories(output_dir);
    for (const Cluster& cluster : clusters) {
        const std::string output_path = output_dir + "/cluster_" + std::to_string(cluster.cluster_id) + ".tsv";
        std::ofstream out(output_path);
        if (!out) {
            throw std::runtime_error("Failed to write cluster signature table: " + output_path);
        }

        out << "node_id";
        for (const std::size_t idx : cluster.member_unique_idxs) {
            out << "\tallele_" << unique_alleles[idx].allele_id;
        }
        out << "\tnode_length_bp\n";

        const std::vector<std::string> node_order = cluster_signature_node_order(cluster, unique_alleles);
        for (const std::string& node_id : node_order) {
            out << node_id;
            for (const std::size_t idx : cluster.member_unique_idxs) {
                out << '\t' << allele_node_orientation_cell(unique_alleles[idx], node_id);
            }
            out << '\t' << node_length_bp(nodes, node_id) << '\n';
        }
    }
}

struct UPGMANode {
    bool is_leaf = false;
    std::size_t allele_idx = 0;
    int left = -1;
    int right = -1;
    double height = 0.0;
    std::vector<std::size_t> leaves;
};

struct UPGMATree {
    std::vector<UPGMANode> nodes;
    int root = -1;
};

UPGMATree build_upgma_tree(const std::vector<std::vector<double>>& dist) {
    UPGMATree tree;
    const std::size_t n = dist.size();
    tree.nodes.reserve(n * 2);
    if (n == 0) {
        return tree;
    }

    std::vector<int> active_nodes;
    active_nodes.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        UPGMANode leaf;
        leaf.is_leaf = true;
        leaf.allele_idx = i;
        leaf.height = 0.0;
        leaf.leaves = {i};
        tree.nodes.push_back(std::move(leaf));
        active_nodes.push_back(static_cast<int>(i));
    }

    while (active_nodes.size() > 1) {
        double best = std::numeric_limits<double>::infinity();
        std::size_t best_i = 0;
        std::size_t best_j = 1;

        for (std::size_t i = 0; i < active_nodes.size(); ++i) {
            for (std::size_t j = i + 1; j < active_nodes.size(); ++j) {
                const auto& a = tree.nodes[static_cast<std::size_t>(active_nodes[i])].leaves;
                const auto& b = tree.nodes[static_cast<std::size_t>(active_nodes[j])].leaves;
                double sum = 0.0;
                std::size_t pairs = 0;
                for (const auto ai : a) {
                    for (const auto bi : b) {
                        sum += dist[ai][bi];
                        ++pairs;
                    }
                }
                if (pairs == 0) {
                    continue;
                }
                const double mean = sum / static_cast<double>(pairs);
                if (mean < best) {
                    best = mean;
                    best_i = i;
                    best_j = j;
                }
            }
        }

        const int left_id = active_nodes[best_i];
        const int right_id = active_nodes[best_j];
        UPGMANode parent;
        parent.is_leaf = false;
        parent.left = left_id;
        parent.right = right_id;
        parent.height = best / 2.0;
        parent.leaves = tree.nodes[static_cast<std::size_t>(left_id)].leaves;
        const auto& right_leaves = tree.nodes[static_cast<std::size_t>(right_id)].leaves;
        parent.leaves.insert(parent.leaves.end(), right_leaves.begin(), right_leaves.end());

        const int parent_id = static_cast<int>(tree.nodes.size());
        tree.nodes.push_back(std::move(parent));

        active_nodes[best_i] = parent_id;
        active_nodes.erase(active_nodes.begin() + static_cast<std::ptrdiff_t>(best_j));
    }

    tree.root = active_nodes.front();
    return tree;
}

SimilarityReportPaths similarity_report_paths_for_bubble(
    const std::string& output_dir,
    std::size_t bubble_id) {

    SimilarityReportPaths out;
    out.bubble_dir = output_dir + "/bubble_" + std::to_string(bubble_id);
    out.alleles_tsv = out.bubble_dir + "/alleles.tsv";
    out.matrix_norm_tsv = out.bubble_dir + "/distance_matrix_norm.tsv";
    out.matrix_abs_tsv = out.bubble_dir + "/distance_matrix_abs.tsv";
    out.stats_tsv = out.bubble_dir + "/cluster_stats.tsv";
    out.cluster_signatures_dir = out.bubble_dir + "/cluster_signatures";
    return out;
}

void write_similarity_reports_for_bubble(
    const SimilarityReportPaths& paths,
    const Bubble& bubble,
    const Graph& graph,
    const std::vector<UniqueAllele>& unique_alleles,
    const std::vector<Cluster>& clusters,
    const std::vector<std::size_t>& cluster_of_unique,
    const DistanceMatrices& dists,
    double min_similarity) {

    if (unique_alleles.empty()) {
        return;
    }

    std::filesystem::create_directories(paths.bubble_dir);
    const std::vector<double> silhouette_scores =
        compute_silhouette_scores(cluster_of_unique, dists.norm);

    write_similarity_alleles_tsv(paths.alleles_tsv, bubble, unique_alleles, cluster_of_unique);
    write_similarity_distance_matrix_tsv(paths.matrix_norm_tsv, unique_alleles, dists.norm);
    write_similarity_distance_matrix_abs_tsv(paths.matrix_abs_tsv, unique_alleles, dists.abs);
    write_similarity_cluster_stats_tsv(
        paths.stats_tsv,
        unique_alleles,
        cluster_of_unique,
        dists.norm,
        silhouette_scores,
        dists.abs,
        min_similarity);
    write_similarity_cluster_signature_tables(
        paths.cluster_signatures_dir,
        graph.nodes,
        unique_alleles,
        clusters);
}

std::string build_walk_signature(const std::vector<PathStep>& steps) {
    std::string sig;
    sig.reserve(steps.size() * 8);
    for (std::size_t i = 0; i < steps.size(); ++i) {
        if (i > 0) {
            sig.push_back(',');
        }
        sig += steps[i].node_id;
        sig.push_back(steps[i].reverse ? '-' : '+');
    }
    return sig;
}

int bounded_levenshtein_distance(
    const std::string& a,
    const std::string& b,
    int max_dist) {

    std::size_t a_begin = 0;
    std::size_t b_begin = 0;
    std::size_t a_end = a.size();
    std::size_t b_end = b.size();

    while (a_begin < a_end && b_begin < b_end && a[a_begin] == b[b_begin]) {
        ++a_begin;
        ++b_begin;
    }
    while (a_end > a_begin && b_end > b_begin && a[a_end - 1] == b[b_end - 1]) {
        --a_end;
        --b_end;
    }

    const int n = static_cast<int>(a_end - a_begin);
    const int m = static_cast<int>(b_end - b_begin);
    if (std::abs(n - m) > max_dist) {
        return max_dist + 1;
    }
    if (n == 0) {
        return m;
    }
    if (m == 0) {
        return n;
    }

    const int inf = max_dist + 1;

    int prev_start = 0;
    int prev_end = std::min(m, max_dist);
    std::vector<int> prev(static_cast<std::size_t>(prev_end - prev_start + 1), inf);
    for (int j = prev_start; j <= prev_end; ++j) {
        prev[static_cast<std::size_t>(j - prev_start)] = j;
    }

    for (int i = 1; i <= n; ++i) {
        const int curr_start = std::max(0, i - max_dist);
        const int curr_end = std::min(m, i + max_dist);
        std::vector<int> curr(static_cast<std::size_t>(curr_end - curr_start + 1), inf);

        int row_best = inf;
        for (int j = curr_start; j <= curr_end; ++j) {
            int best = inf;
            if (j == 0) {
                best = i;
            } else {
                int deletion = inf;
                int insertion = inf;
                int substitution = inf;

                if (j >= prev_start && j <= prev_end) {
                    deletion = prev[static_cast<std::size_t>(j - prev_start)] + 1;
                }
                if (j - 1 >= curr_start) {
                    insertion = curr[static_cast<std::size_t>((j - 1) - curr_start)] + 1;
                }
                if (j - 1 >= prev_start && j - 1 <= prev_end) {
                    const int cost =
                        (a[a_begin + static_cast<std::size_t>(i - 1)] ==
                         b[b_begin + static_cast<std::size_t>(j - 1)])
                            ? 0
                            : 1;
                    substitution = prev[static_cast<std::size_t>((j - 1) - prev_start)] + cost;
                }
                best = std::min({deletion, insertion, substitution});
            }

            curr[static_cast<std::size_t>(j - curr_start)] = best;
            row_best = std::min(row_best, best);
        }

        if (row_best > max_dist) {
            return max_dist + 1;
        }

        prev.swap(curr);
        prev_start = curr_start;
        prev_end = curr_end;
    }

    if (m < prev_start || m > prev_end) {
        return max_dist + 1;
    }
    return prev[static_cast<std::size_t>(m - prev_start)];
}

int bounded_weighted_levenshtein_distance(
    const std::vector<std::uint64_t>& a,
    const std::vector<int>& a_weights,
    const std::vector<std::uint64_t>& b,
    const std::vector<int>& b_weights,
    int max_dist) {

    if (a.size() != a_weights.size() || b.size() != b_weights.size()) {
        throw std::runtime_error("Weighted walk token/weight length mismatch");
    }

    const int inf = std::max(1, max_dist + 1);
    auto saturating_add = [&](int x, int y) -> int {
        if (x >= inf || y >= inf) {
            return inf;
        }
        if (x > inf - y) {
            return inf;
        }
        return x + y;
    };

    const std::size_t n = a.size();
    const std::size_t m = b.size();
    if (n == 0) {
        int total = 0;
        for (const int w : b_weights) {
            total = saturating_add(total, std::max(1, w));
            if (total > max_dist) {
                return max_dist + 1;
            }
        }
        return total;
    }
    if (m == 0) {
        int total = 0;
        for (const int w : a_weights) {
            total = saturating_add(total, std::max(1, w));
            if (total > max_dist) {
                return max_dist + 1;
            }
        }
        return total;
    }

    int total_a = 0;
    for (const int w : a_weights) {
        total_a = saturating_add(total_a, std::max(1, w));
    }
    int total_b = 0;
    for (const int w : b_weights) {
        total_b = saturating_add(total_b, std::max(1, w));
    }
    if (std::abs(total_a - total_b) > max_dist) {
        return max_dist + 1;
    }

    std::vector<int> prev(m + 1, inf);
    prev[0] = 0;
    for (std::size_t j = 1; j <= m; ++j) {
        prev[j] = saturating_add(prev[j - 1], std::max(1, b_weights[j - 1]));
    }
    if (*std::min_element(prev.begin(), prev.end()) > max_dist) {
        return max_dist + 1;
    }

    std::vector<int> curr(m + 1, inf);
    for (std::size_t i = 1; i <= n; ++i) {
        const int del_cost = std::max(1, a_weights[i - 1]);
        curr[0] = saturating_add(prev[0], del_cost);
        int row_best = curr[0];

        for (std::size_t j = 1; j <= m; ++j) {
            const int ins_cost = std::max(1, b_weights[j - 1]);
            const int sub_cost = (a[i - 1] == b[j - 1]) ? 0 : std::max(del_cost, ins_cost);

            const int deletion = saturating_add(prev[j], del_cost);
            const int insertion = saturating_add(curr[j - 1], ins_cost);
            const int substitution = saturating_add(prev[j - 1], sub_cost);
            const int best = std::min({deletion, insertion, substitution});
            curr[j] = best;
            row_best = std::min(row_best, best);
        }

        if (row_best > max_dist) {
            return max_dist + 1;
        }
        prev.swap(curr);
    }

    if (prev[m] > max_dist) {
        return max_dist + 1;
    }
    return prev[m];
}

void collect_subtree_leaves(
    const UPGMATree& tree,
    int node_id,
    std::vector<std::size_t>& out) {

    if (node_id < 0) {
        return;
    }
    const auto& node = tree.nodes[static_cast<std::size_t>(node_id)];
    if (node.is_leaf) {
        out.push_back(node.allele_idx);
        return;
    }
    collect_subtree_leaves(tree, node.left, out);
    collect_subtree_leaves(tree, node.right, out);
}

std::vector<Cluster> cluster_unique_alleles_from_tree(
    const std::vector<UniqueAllele>& unique_alleles,
    const std::vector<std::vector<double>>& dist_norm,
    const UPGMATree& tree,
    double min_similarity) {

    if (unique_alleles.empty() || tree.root < 0) {
        return {};
    }

    const double max_norm_dist = std::clamp(1.0 - min_similarity, 0.0, 1.0);
    // UPGMA node height stores half of the merge distance.
    const double cut_height = max_norm_dist / 2.0;
    const double eps = 1e-12;

    std::vector<std::vector<std::size_t>> member_sets;
    member_sets.reserve(unique_alleles.size());

    std::function<void(int)> cut = [&](int node_id) {
        const auto& node = tree.nodes[static_cast<std::size_t>(node_id)];
        if (node.is_leaf) {
            member_sets.push_back({node.allele_idx});
            return;
        }

        if (node.height <= cut_height + eps) {
            std::vector<std::size_t> leaves;
            leaves.reserve(node.leaves.size());
            collect_subtree_leaves(tree, node_id, leaves);
            member_sets.push_back(std::move(leaves));
            return;
        }

        cut(node.left);
        cut(node.right);
    };
    cut(tree.root);

    std::vector<Cluster> clusters;
    clusters.reserve(member_sets.size());

    for (auto& members : member_sets) {
        if (members.empty()) {
            continue;
        }

        Cluster cluster;
        cluster.cluster_id = clusters.size() + 1;
        cluster.member_unique_idxs = std::move(members);
        cluster.total_path_support = 0;

        std::size_t rep_idx = cluster.member_unique_idxs.front();
        double best_max_dist = std::numeric_limits<double>::infinity();
        double best_mean_dist = std::numeric_limits<double>::infinity();
        for (const auto idx : cluster.member_unique_idxs) {
            cluster.total_path_support += unique_alleles[idx].path_support;

            double sum = 0.0;
            double max_d = 0.0;
            std::size_t count = 0;
            for (const auto jdx : cluster.member_unique_idxs) {
                if (idx == jdx) {
                    continue;
                }
                const double d = dist_norm[idx][jdx];
                sum += d;
                max_d = std::max(max_d, d);
                ++count;
            }
            const double mean_d = safe_div(sum, static_cast<double>(count));

            bool choose = false;
            if (max_d + 1e-12 < best_max_dist) {
                choose = true;
            } else if (std::abs(max_d - best_max_dist) <= 1e-12) {
                if (mean_d + 1e-12 < best_mean_dist) {
                    choose = true;
                } else if (std::abs(mean_d - best_mean_dist) <= 1e-12) {
                    const auto& cand = unique_alleles[idx];
                    const auto& curr = unique_alleles[rep_idx];
                    if (cand.path_support > curr.path_support) {
                        choose = true;
                    } else if (cand.path_support == curr.path_support &&
                               cand.allele_id < curr.allele_id) {
                        choose = true;
                    }
                }
            }

            if (choose) {
                rep_idx = idx;
                best_max_dist = max_d;
                best_mean_dist = mean_d;
            }
        }
        cluster.representative_idx = rep_idx;

        std::sort(cluster.member_unique_idxs.begin(), cluster.member_unique_idxs.end(), [&](std::size_t lhs, std::size_t rhs) {
            if (unique_alleles[lhs].path_support != unique_alleles[rhs].path_support) {
                return unique_alleles[lhs].path_support > unique_alleles[rhs].path_support;
            }
            return unique_alleles[lhs].allele_id < unique_alleles[rhs].allele_id;
        });
        clusters.push_back(std::move(cluster));
    }

    return clusters;
}

struct DisjointSet {
    explicit DisjointSet(std::size_t n) : parent(n), rank(n, 0) {
        std::iota(parent.begin(), parent.end(), 0);
    }

    std::size_t find(std::size_t x) {
        while (parent[x] != x) {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        return x;
    }

    void unite(std::size_t a, std::size_t b) {
        a = find(a);
        b = find(b);
        if (a == b) {
            return;
        }
        if (rank[a] < rank[b]) {
            std::swap(a, b);
        }
        parent[b] = a;
        if (rank[a] == rank[b]) {
            ++rank[a];
        }
    }

    std::vector<std::size_t> parent;
    std::vector<unsigned char> rank;
};

std::vector<Cluster> cluster_unique_alleles_threshold_graph(
    const std::vector<UniqueAllele>& unique_alleles,
    const std::vector<std::vector<double>>& dist_norm,
    double min_similarity) {

    const std::size_t n = unique_alleles.size();
    if (n == 0) {
        return {};
    }

    const double max_norm_dist = std::clamp(1.0 - min_similarity, 0.0, 1.0);
    const double eps = 1e-12;
    DisjointSet dsu(n);

    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i + 1; j < n; ++j) {
            if (dist_norm[i][j] <= max_norm_dist + eps) {
                dsu.unite(i, j);
            }
        }
    }

    std::unordered_map<std::size_t, std::vector<std::size_t>> members_by_root;
    members_by_root.reserve(n * 2);
    for (std::size_t i = 0; i < n; ++i) {
        members_by_root[dsu.find(i)].push_back(i);
    }

    std::vector<std::vector<std::size_t>> member_sets;
    member_sets.reserve(members_by_root.size());
    for (auto& kv : members_by_root) {
        member_sets.push_back(std::move(kv.second));
    }
    std::sort(member_sets.begin(), member_sets.end(), [&](const auto& a, const auto& b) {
        const auto min_a = std::min_element(a.begin(), a.end(), [&](std::size_t lhs, std::size_t rhs) {
            return unique_alleles[lhs].allele_id < unique_alleles[rhs].allele_id;
        });
        const auto min_b = std::min_element(b.begin(), b.end(), [&](std::size_t lhs, std::size_t rhs) {
            return unique_alleles[lhs].allele_id < unique_alleles[rhs].allele_id;
        });
        return unique_alleles[*min_a].allele_id < unique_alleles[*min_b].allele_id;
    });

    std::vector<Cluster> clusters;
    clusters.reserve(member_sets.size());
    for (auto& members : member_sets) {
        if (members.empty()) {
            continue;
        }

        Cluster cluster;
        cluster.cluster_id = clusters.size() + 1;
        cluster.member_unique_idxs = std::move(members);
        cluster.total_path_support = 0;

        std::size_t rep_idx = cluster.member_unique_idxs.front();
        double best_max_dist = std::numeric_limits<double>::infinity();
        double best_mean_dist = std::numeric_limits<double>::infinity();
        for (const auto idx : cluster.member_unique_idxs) {
            cluster.total_path_support += unique_alleles[idx].path_support;

            double sum = 0.0;
            double max_d = 0.0;
            std::size_t count = 0;
            for (const auto jdx : cluster.member_unique_idxs) {
                if (idx == jdx) {
                    continue;
                }
                const double d = dist_norm[idx][jdx];
                sum += d;
                max_d = std::max(max_d, d);
                ++count;
            }
            const double mean_d = safe_div(sum, static_cast<double>(count));

            bool choose = false;
            if (max_d + 1e-12 < best_max_dist) {
                choose = true;
            } else if (std::abs(max_d - best_max_dist) <= 1e-12) {
                if (mean_d + 1e-12 < best_mean_dist) {
                    choose = true;
                } else if (std::abs(mean_d - best_mean_dist) <= 1e-12) {
                    const auto& cand = unique_alleles[idx];
                    const auto& curr = unique_alleles[rep_idx];
                    if (cand.path_support > curr.path_support) {
                        choose = true;
                    } else if (cand.path_support == curr.path_support &&
                               cand.allele_id < curr.allele_id) {
                        choose = true;
                    }
                }
            }
            if (choose) {
                rep_idx = idx;
                best_max_dist = max_d;
                best_mean_dist = mean_d;
            }
        }
        cluster.representative_idx = rep_idx;

        std::sort(cluster.member_unique_idxs.begin(), cluster.member_unique_idxs.end(), [&](std::size_t lhs, std::size_t rhs) {
            if (unique_alleles[lhs].path_support != unique_alleles[rhs].path_support) {
                return unique_alleles[lhs].path_support > unique_alleles[rhs].path_support;
            }
            return unique_alleles[lhs].allele_id < unique_alleles[rhs].allele_id;
        });
        clusters.push_back(std::move(cluster));
    }

    return clusters;
}

std::vector<Cluster> cluster_unique_alleles_greedy_threshold(
    const std::vector<UniqueAllele>& unique_alleles,
    const std::vector<std::vector<double>>& dist_norm,
    double min_similarity) {

    const std::size_t n = unique_alleles.size();
    if (n == 0) {
        return {};
    }
    const double max_norm_dist = std::clamp(1.0 - min_similarity, 0.0, 1.0);
    const double eps = 1e-12;

    std::vector<std::size_t> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](std::size_t lhs, std::size_t rhs) {
        if (unique_alleles[lhs].path_support != unique_alleles[rhs].path_support) {
            return unique_alleles[lhs].path_support > unique_alleles[rhs].path_support;
        }
        return unique_alleles[lhs].allele_id < unique_alleles[rhs].allele_id;
    });

    std::vector<Cluster> clusters;
    clusters.reserve(n);

    for (const auto idx : order) {
        std::size_t best_cluster = std::numeric_limits<std::size_t>::max();
        double best_dist = std::numeric_limits<double>::infinity();

        for (std::size_t c = 0; c < clusters.size(); ++c) {
            const double d = dist_norm[idx][clusters[c].representative_idx];
            if (d <= max_norm_dist + eps) {
                if (d + eps < best_dist) {
                    best_dist = d;
                    best_cluster = c;
                } else if (std::abs(d - best_dist) <= eps) {
                    const auto curr_rep = clusters[c].representative_idx;
                    const auto best_rep = clusters[best_cluster].representative_idx;
                    if (unique_alleles[curr_rep].path_support > unique_alleles[best_rep].path_support ||
                        (unique_alleles[curr_rep].path_support == unique_alleles[best_rep].path_support &&
                         unique_alleles[curr_rep].allele_id < unique_alleles[best_rep].allele_id)) {
                        best_cluster = c;
                    }
                }
            }
        }

        if (best_cluster == std::numeric_limits<std::size_t>::max()) {
            Cluster cluster;
            cluster.cluster_id = clusters.size() + 1;
            cluster.representative_idx = idx;
            cluster.member_unique_idxs.push_back(idx);
            cluster.total_path_support = unique_alleles[idx].path_support;
            clusters.push_back(std::move(cluster));
        } else {
            auto& cluster = clusters[best_cluster];
            cluster.member_unique_idxs.push_back(idx);
            cluster.total_path_support += unique_alleles[idx].path_support;
        }
    }

    for (auto& cluster : clusters) {
        std::size_t rep_idx = cluster.representative_idx;
        double best_max_dist = std::numeric_limits<double>::infinity();
        double best_mean_dist = std::numeric_limits<double>::infinity();
        for (const auto idx : cluster.member_unique_idxs) {
            double sum = 0.0;
            double max_d = 0.0;
            std::size_t count = 0;
            for (const auto jdx : cluster.member_unique_idxs) {
                if (idx == jdx) {
                    continue;
                }
                const double d = dist_norm[idx][jdx];
                sum += d;
                max_d = std::max(max_d, d);
                ++count;
            }
            const double mean_d = safe_div(sum, static_cast<double>(count));

            bool choose = false;
            if (max_d + eps < best_max_dist) {
                choose = true;
            } else if (std::abs(max_d - best_max_dist) <= eps) {
                if (mean_d + eps < best_mean_dist) {
                    choose = true;
                } else if (std::abs(mean_d - best_mean_dist) <= eps) {
                    const auto& cand = unique_alleles[idx];
                    const auto& curr = unique_alleles[rep_idx];
                    if (cand.path_support > curr.path_support ||
                        (cand.path_support == curr.path_support &&
                         cand.allele_id < curr.allele_id)) {
                        choose = true;
                    }
                }
            }
            if (choose) {
                rep_idx = idx;
                best_max_dist = max_d;
                best_mean_dist = mean_d;
            }
        }
        cluster.representative_idx = rep_idx;

        std::sort(cluster.member_unique_idxs.begin(), cluster.member_unique_idxs.end(), [&](std::size_t lhs, std::size_t rhs) {
            if (unique_alleles[lhs].path_support != unique_alleles[rhs].path_support) {
                return unique_alleles[lhs].path_support > unique_alleles[rhs].path_support;
            }
            return unique_alleles[lhs].allele_id < unique_alleles[rhs].allele_id;
        });
    }

    return clusters;
}

void skip_json_whitespace(const std::string& text, std::size_t& pos) {
    while (pos < text.size() &&
           std::isspace(static_cast<unsigned char>(text[pos])) != 0) {
        ++pos;
    }
}

std::string parse_json_quoted_string(const std::string& text, std::size_t& pos) {
    if (pos >= text.size() || text[pos] != '"') {
        throw std::runtime_error("Invalid predefined cluster JSON: expected quoted string");
    }
    ++pos;
    std::string out;
    out.reserve(32);
    while (pos < text.size()) {
        const char c = text[pos++];
        if (c == '"') {
            return out;
        }
        if (c != '\\') {
            out.push_back(c);
            continue;
        }
        if (pos >= text.size()) {
            throw std::runtime_error("Invalid predefined cluster JSON: dangling escape");
        }
        const char esc = text[pos++];
        switch (esc) {
            case '"':
            case '\\':
            case '/':
                out.push_back(esc);
                break;
            case 'b':
                out.push_back('\b');
                break;
            case 'f':
                out.push_back('\f');
                break;
            case 'n':
                out.push_back('\n');
                break;
            case 'r':
                out.push_back('\r');
                break;
            case 't':
                out.push_back('\t');
                break;
            default:
                throw std::runtime_error(
                    "Invalid predefined cluster JSON: unsupported escape sequence");
        }
    }
    throw std::runtime_error("Invalid predefined cluster JSON: unterminated string");
}

std::unordered_map<std::string, std::string> load_predefined_path_cluster_map_json(
    const std::string& json_path) {

    std::ifstream in(json_path);
    if (!in) {
        throw std::runtime_error("Failed to read predefined clusters JSON: " + json_path);
    }
    std::string text(
        (std::istreambuf_iterator<char>(in)),
        std::istreambuf_iterator<char>());
    std::size_t pos = 0;
    skip_json_whitespace(text, pos);
    if (pos >= text.size() || text[pos] != '{') {
        throw std::runtime_error("Invalid predefined clusters JSON: expected top-level object");
    }
    ++pos;

    std::unordered_map<std::string, std::string> out;
    out.reserve(1024);

    while (true) {
        skip_json_whitespace(text, pos);
        if (pos >= text.size()) {
            throw std::runtime_error("Invalid predefined clusters JSON: unexpected EOF");
        }
        if (text[pos] == '}') {
            ++pos;
            break;
        }

        const std::string key = parse_json_quoted_string(text, pos);
        skip_json_whitespace(text, pos);
        if (pos >= text.size() || text[pos] != ':') {
            throw std::runtime_error("Invalid predefined clusters JSON: expected ':'");
        }
        ++pos;
        skip_json_whitespace(text, pos);
        if (pos >= text.size()) {
            throw std::runtime_error("Invalid predefined clusters JSON: missing value");
        }

        std::string value;
        if (text[pos] == '"') {
            value = parse_json_quoted_string(text, pos);
        } else {
            const std::size_t value_start = pos;
            while (pos < text.size() && text[pos] != ',' && text[pos] != '}') {
                ++pos;
            }
            value = text.substr(value_start, pos - value_start);
            std::size_t lo = 0;
            while (lo < value.size() &&
                   std::isspace(static_cast<unsigned char>(value[lo])) != 0) {
                ++lo;
            }
            std::size_t hi = value.size();
            while (hi > lo &&
                   std::isspace(static_cast<unsigned char>(value[hi - 1])) != 0) {
                --hi;
            }
            value = value.substr(lo, hi - lo);
        }
        if (value.empty()) {
            value = "__EMPTY__";
        }
        out[key] = value;

        skip_json_whitespace(text, pos);
        if (pos >= text.size()) {
            throw std::runtime_error("Invalid predefined clusters JSON: unexpected EOF");
        }
        if (text[pos] == ',') {
            ++pos;
            continue;
        }
        if (text[pos] == '}') {
            ++pos;
            break;
        }
        throw std::runtime_error("Invalid predefined clusters JSON: expected ',' or '}'");
    }

    skip_json_whitespace(text, pos);
    if (pos != text.size()) {
        throw std::runtime_error("Invalid predefined clusters JSON: trailing content");
    }
    return out;
}

std::vector<Cluster> cluster_unique_alleles_from_predefined_path_labels(
    const std::vector<UniqueAllele>& unique_alleles,
    const std::vector<PathAssignment>& assignments,
    const std::unordered_map<std::string, std::string>& label_by_path,
    std::size_t* missing_assignment_count_out,
    std::size_t* missing_path_count_out) {

    if (missing_assignment_count_out != nullptr) {
        *missing_assignment_count_out = 0;
    }
    if (missing_path_count_out != nullptr) {
        *missing_path_count_out = 0;
    }
    if (unique_alleles.empty() || assignments.empty()) {
        return {};
    }

    std::unordered_map<std::string, std::vector<std::size_t>> unique_members_by_label;
    std::unordered_map<std::string, std::size_t> support_by_label;
    std::unordered_set<std::string> missing_paths;
    unique_members_by_label.reserve(assignments.size() * 2);
    support_by_label.reserve(assignments.size() * 2);

    std::size_t missing_assignment_count = 0;
    for (const auto& assignment : assignments) {
        auto it = label_by_path.find(assignment.path_name);
        bool missing = (it == label_by_path.end());
        std::string label =
            missing ? ("__UNMAPPED__:" + assignment.path_name) : it->second;
        if (label.empty()) {
            label = "__EMPTY__";
        }
        if (missing) {
            ++missing_assignment_count;
            missing_paths.insert(assignment.path_name);
        }
        unique_members_by_label[label].push_back(assignment.unique_idx);
        support_by_label[label] += 1;
    }
    if (missing_assignment_count_out != nullptr) {
        *missing_assignment_count_out = missing_assignment_count;
    }
    if (missing_path_count_out != nullptr) {
        *missing_path_count_out = missing_paths.size();
    }

    std::vector<std::string> labels;
    labels.reserve(unique_members_by_label.size());
    for (const auto& kv : unique_members_by_label) {
        labels.push_back(kv.first);
    }
    std::sort(labels.begin(), labels.end());

    std::vector<Cluster> clusters;
    clusters.reserve(labels.size());
    for (std::size_t i = 0; i < labels.size(); ++i) {
        const std::string& label = labels[i];
        auto members = unique_members_by_label[label];
        std::sort(members.begin(), members.end());
        members.erase(std::unique(members.begin(), members.end()), members.end());
        if (members.empty()) {
            continue;
        }

        Cluster cluster;
        cluster.cluster_id = i + 1;
        cluster.member_unique_idxs = std::move(members);
        cluster.total_path_support = support_by_label[label];
        if (cluster.total_path_support == 0) {
            for (const std::size_t unique_idx : cluster.member_unique_idxs) {
                cluster.total_path_support += unique_alleles[unique_idx].path_support;
            }
        }

        std::size_t rep_idx = cluster.member_unique_idxs.front();
        for (const std::size_t unique_idx : cluster.member_unique_idxs) {
            const auto& cand = unique_alleles[unique_idx];
            const auto& curr = unique_alleles[rep_idx];
            if (cand.path_support > curr.path_support ||
                (cand.path_support == curr.path_support && cand.allele_id < curr.allele_id)) {
                rep_idx = unique_idx;
            }
        }
        cluster.representative_idx = rep_idx;

        std::sort(cluster.member_unique_idxs.begin(), cluster.member_unique_idxs.end(), [&](std::size_t lhs, std::size_t rhs) {
            if (unique_alleles[lhs].path_support != unique_alleles[rhs].path_support) {
                return unique_alleles[lhs].path_support > unique_alleles[rhs].path_support;
            }
            return unique_alleles[lhs].allele_id < unique_alleles[rhs].allele_id;
        });

        clusters.push_back(std::move(cluster));
    }

    return clusters;
}

struct DotplotPoint {
    std::size_t ref_pos = 0;
    std::size_t query_pos = 0;
    char strand = 'F';
};

struct AtomicVariantEvent {
    std::size_t ref_offset_start_bp = 0;
    std::size_t ref_offset_end_bp = 0;
    std::size_t alt_offset_start_bp = 0;
    std::size_t alt_offset_end_bp = 0;
    std::string event_type;
    std::string event_subtype = ".";
    std::size_t event_bp = 0;
    long long svlen = 0;
    std::size_t inserted_bp = 0;
    std::string inserted_seq;
    bool has_dup_evidence = false;
    double dup_best_similarity = 0.0;
    std::size_t dup_ref_start_bp = 0;
    std::size_t dup_ref_end_bp = 0;
    char dup_orientation = '+';
    std::size_t dup_unit_bp = 0;
    std::size_t dup_ref_copy_number = 0;
    std::size_t dup_added_copies = 0;
    std::size_t dup_alt_copy_number = 0;
    double dup_copy_ratio = 0.0;
    int cn_delta = 0;
    bool preserve_ins_svlen = false;
    std::string evidence_source = "unknown";
};

struct CigarRun {
    char op = 'M';
    std::size_t len = 0;
};

struct MinimapPrimaryAlignment {
    bool ok = false;
    bool reverse = false;
    int edit_distance = -1;
    std::size_t query_start_bp = 0;
    std::size_t query_end_bp = 0;
    std::size_t target_start_bp = 0;
    std::size_t target_end_bp = 0;
    std::string cigar_extended;
};

std::size_t longest_common_prefix_bp(const std::string& a, const std::string& b) {
    const std::size_t max_i = std::min(a.size(), b.size());
    std::size_t i = 0;
    while (i < max_i && a[i] == b[i]) {
        ++i;
    }
    return i;
}

std::size_t longest_common_suffix_bp(
    const std::string& a,
    const std::string& b,
    std::size_t skip_prefix) {

    const std::size_t max_sfx = std::min(a.size(), b.size());
    std::size_t sfx = 0;
    while (sfx < max_sfx &&
           (skip_prefix + sfx) < a.size() &&
           (skip_prefix + sfx) < b.size() &&
           a[a.size() - 1 - sfx] == b[b.size() - 1 - sfx]) {
        ++sfx;
    }
    return sfx;
}

ParsedReferencePath parse_reference_path_label(const std::string& path_name) {
    ParsedReferencePath out;
    out.chrom = path_name;

    auto pick_chrom = [&](const std::string& prefix) {
        if (prefix.empty()) {
            return std::string(path_name);
        }
        const std::size_t hash = prefix.rfind('#');
        if (hash != std::string::npos && hash + 1 < prefix.size()) {
            return prefix.substr(hash + 1);
        }
        const std::size_t chr_pos = prefix.find("chr");
        if (chr_pos != std::string::npos) {
            return prefix.substr(chr_pos);
        }
        return prefix;
    };

    const std::size_t colon = path_name.rfind(':');
    const std::size_t dash = path_name.rfind('-');
    if (colon != std::string::npos && dash != std::string::npos && dash > colon + 1 && dash + 1 < path_name.size()) {
        const std::string start_str = path_name.substr(colon + 1, dash - colon - 1);
        const std::string end_str = path_name.substr(dash + 1);
        bool digits_only = !start_str.empty() && !end_str.empty();
        for (const char c : start_str) {
            if (!std::isdigit(static_cast<unsigned char>(c))) {
                digits_only = false;
                break;
            }
        }
        if (digits_only) {
            for (const char c : end_str) {
                if (!std::isdigit(static_cast<unsigned char>(c))) {
                    digits_only = false;
                    break;
                }
            }
        }
        if (digits_only) {
            try {
                const std::size_t parsed_start = static_cast<std::size_t>(std::stoull(start_str));
                const std::size_t parsed_end = static_cast<std::size_t>(std::stoull(end_str));
                if (parsed_start > 0 && parsed_end >= parsed_start) {
                    out.has_interval = true;
                    out.region_start_1based = parsed_start;
                    out.chrom = pick_chrom(path_name.substr(0, colon));
                    return out;
                }
            } catch (const std::exception&) {
            }
        }
    }

    out.chrom = pick_chrom(path_name);
    return out;
}

char vcf_base(char c) {
    const char up = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    switch (up) {
        case 'A':
        case 'C':
        case 'G':
        case 'T':
        case 'N':
            return up;
        default:
            return 'N';
    }
}

std::vector<CigarRun> parse_extended_cigar(const std::string& cigar) {
    std::vector<CigarRun> runs;
    if (cigar.empty()) {
        return runs;
    }

    std::size_t i = 0;
    while (i < cigar.size()) {
        if (!std::isdigit(static_cast<unsigned char>(cigar[i]))) {
            throw std::runtime_error("Malformed CIGAR: missing length before op");
        }
        std::size_t len = 0;
        while (i < cigar.size() && std::isdigit(static_cast<unsigned char>(cigar[i]))) {
            len = (len * 10) + static_cast<std::size_t>(cigar[i] - '0');
            ++i;
        }
        if (i >= cigar.size()) {
            throw std::runtime_error("Malformed CIGAR: missing op");
        }
        const char op = cigar[i++];
        if (len == 0) {
            continue;
        }
        runs.push_back({op, len});
    }
    return runs;
}

struct CigarTraceRow {
    std::size_t op_index = 0;
    char op = 'M';
    std::size_t len = 0;
    std::size_t ref_start_bp = 0;
    std::size_t ref_end_bp = 0;
    std::size_t alt_start_bp = 0;
    std::size_t alt_end_bp = 0;
    std::string interpreted_event = ".";
    bool sv_sized = false;
};

[[maybe_unused]] std::string interpreted_event_from_cigar_op(char op) {
    if (op == '=' || op == 'M') {
        return "MATCH";
    }
    if (op == 'X') {
        return "MISMATCH";
    }
    if (op == 'I') {
        // query=ALT, target=reference.
        return "INS(ref->alt)";
    }
    if (op == 'D') {
        // query=ALT, target=reference.
        return "DEL(ref->alt)";
    }
    return "UNSUPPORTED";
}

[[maybe_unused]] std::vector<CigarTraceRow> build_cigar_trace_rows(
    const std::vector<CigarRun>& runs,
    std::size_t min_sv_bp) {

    std::vector<CigarTraceRow> rows;
    rows.reserve(runs.size());
    std::size_t ref_pos = 0;
    std::size_t alt_pos = 0;
    for (std::size_t i = 0; i < runs.size(); ++i) {
        const auto& run = runs[i];
        CigarTraceRow row;
        row.op_index = i + 1;
        row.op = run.op;
        row.len = run.len;
        row.ref_start_bp = ref_pos;
        row.alt_start_bp = alt_pos;
        row.interpreted_event = interpreted_event_from_cigar_op(run.op);
        row.sv_sized = (run.op == 'I' || run.op == 'D') && run.len >= min_sv_bp;

        if (run.op == '=' || run.op == 'M' || run.op == 'X') {
            ref_pos += run.len;
            alt_pos += run.len;
        } else if (run.op == 'I') {
            alt_pos += run.len;
        } else if (run.op == 'D') {
            ref_pos += run.len;
        } else {
            throw std::runtime_error(std::string("Unsupported CIGAR op: ") + run.op);
        }

        row.ref_end_bp = ref_pos;
        row.alt_end_bp = alt_pos;
        rows.push_back(std::move(row));
    }
    return rows;
}

[[maybe_unused]] void write_cigar_trace_tsv(const std::string& output_path, const std::vector<CigarTraceRow>& rows) {
    std::ofstream out(output_path);
    if (!out) {
        throw std::runtime_error("Failed to write CIGAR trace TSV: " + output_path);
    }
    out << "op_index\top\tlen\tref_start_bp\tref_end_bp\talt_start_bp\talt_end_bp\tinterpreted_event\tsv_sized\n";
    for (const auto& row : rows) {
        out << row.op_index << '\t'
            << row.op << '\t'
            << row.len << '\t'
            << row.ref_start_bp << '\t'
            << row.ref_end_bp << '\t'
            << row.alt_start_bp << '\t'
            << row.alt_end_bp << '\t'
            << row.interpreted_event << '\t'
            << (row.sv_sized ? 1 : 0) << '\n';
    }
}

std::string preview_text(const std::string& text, std::size_t max_chars) {
    if (text.size() <= max_chars) {
        return text;
    }
    if (max_chars < 8) {
        return text.substr(0, max_chars);
    }
    const std::size_t keep = (max_chars - 3) / 2;
    return text.substr(0, keep) + "..." + text.substr(text.size() - keep);
}

std::vector<AtomicVariantEvent> extract_atomic_events_from_cigar(
    const std::vector<CigarRun>& runs,
    const std::string& alt_seq,
    std::size_t min_sv_bp,
    std::size_t ref_start_bp = 0,
    std::size_t alt_start_bp = 0) {

    std::vector<AtomicVariantEvent> events;
    std::size_t ref_pos = ref_start_bp;
    std::size_t alt_pos = alt_start_bp;
    for (const auto& run : runs) {
        if (run.op == '=' || run.op == 'M' || run.op == 'X') {
            ref_pos += run.len;
            alt_pos += run.len;
            continue;
        }
        // minimap2 CIGAR semantics:
        //   I = insertion in query (ALT) relative to reference
        //   D = deletion in query (ALT) relative to reference
        // Here query=ALT and target=reference.
        if (run.op == 'I') {
            if (run.len >= min_sv_bp) {
                AtomicVariantEvent ev;
                ev.event_type = "INS";
                ev.ref_offset_start_bp = ref_pos;
                ev.ref_offset_end_bp = ref_pos;
                ev.alt_offset_start_bp = alt_pos;
                ev.alt_offset_end_bp = alt_pos + run.len;
                ev.inserted_bp = run.len;
                ev.event_bp = run.len;
                ev.svlen = static_cast<long long>(run.len);
                ev.inserted_seq = alt_seq.substr(alt_pos, run.len);
                ev.evidence_source = "cigar";
                events.push_back(std::move(ev));
            }
            alt_pos += run.len;
            continue;
        }
        if (run.op == 'D') {
            if (run.len >= min_sv_bp) {
                AtomicVariantEvent ev;
                ev.event_type = "DEL";
                ev.ref_offset_start_bp = ref_pos;
                ev.ref_offset_end_bp = ref_pos + run.len;
                ev.alt_offset_start_bp = alt_pos;
                ev.alt_offset_end_bp = alt_pos;
                ev.event_bp = run.len;
                ev.svlen = -static_cast<long long>(run.len);
                ev.cn_delta = -1;
                ev.evidence_source = "cigar";
                events.push_back(std::move(ev));
            }
            ref_pos += run.len;
            continue;
        }
        throw std::runtime_error(std::string("Unsupported CIGAR op: ") + run.op);
    }
    return events;
}

struct PafRecord {
    std::string query_name;
    std::size_t query_len = 0;
    std::size_t query_start = 0;
    std::size_t query_end = 0;
    char strand = '+';
    std::string target_name;
    std::size_t target_len = 0;
    std::size_t target_start = 0;
    std::size_t target_end = 0;
    std::size_t n_matches = 0;
    std::size_t aln_block_len = 0;
    int mapq = 60;
    bool primary = true;
    int cm = -1;
    int s1 = std::numeric_limits<int>::min();
    int s2 = std::numeric_limits<int>::min();
    double dv = std::numeric_limits<double>::quiet_NaN();
    std::string cigar;
};

struct SplitAlignmentSegment {
    std::size_t query_start_bp = 0;
    std::size_t query_end_bp = 0;
    std::size_t ref_start_bp = 0;
    std::size_t ref_end_bp = 0;
    bool reverse = false;
};

std::vector<AtomicVariantEvent> extract_split_events_from_records(
    const std::vector<PafRecord>& records,
    char required_strand,
    std::size_t min_sv_bp,
    std::size_t max_sv_bp,
    std::size_t query_len,
    bool split_ins_use_geometric_svlen) {

    std::vector<SplitAlignmentSegment> segments;
    segments.reserve(records.size());
    for (const auto& rec : records) {
        if (required_strand == '+' || required_strand == '-') {
            if (rec.strand != required_strand) {
                continue;
            }
        }
        if (rec.query_end <= rec.query_start || rec.target_end <= rec.target_start) {
            continue;
        }

        std::size_t q_start = rec.query_start;
        std::size_t q_end = rec.query_end;
        const bool reverse = (rec.strand == '-');
        if (reverse) {
            if (query_len < rec.query_end || query_len < rec.query_start) {
                continue;
            }
            q_start = query_len - rec.query_end;
            q_end = query_len - rec.query_start;
        }
        if (q_end <= q_start) {
            continue;
        }
        segments.push_back({
            q_start,
            q_end,
            rec.target_start,
            rec.target_end,
            reverse});
    }

    if (segments.size() < 2) {
        return {};
    }

    std::sort(segments.begin(), segments.end(), [](const SplitAlignmentSegment& lhs, const SplitAlignmentSegment& rhs) {
        if (lhs.query_start_bp != rhs.query_start_bp) {
            return lhs.query_start_bp < rhs.query_start_bp;
        }
        if (lhs.query_end_bp != rhs.query_end_bp) {
            return lhs.query_end_bp < rhs.query_end_bp;
        }
        if (lhs.ref_start_bp != rhs.ref_start_bp) {
            return lhs.ref_start_bp < rhs.ref_start_bp;
        }
        return lhs.ref_end_bp < rhs.ref_end_bp;
    });

    constexpr long long kQueryGapTolerance = 50;
    constexpr long long kQueryOverlapTolerance = 50;
    constexpr long long kRefGapTolerance = 50;
    constexpr long long kRefOverlapTolerance = 50;

    std::vector<AtomicVariantEvent> out;
    out.reserve(segments.size());
    for (std::size_t i = 0; i + 1 < segments.size(); ++i) {
        const auto& curr = segments[i];
        const auto& next = segments[i + 1];
        if (curr.reverse != next.reverse) {
            continue;
        }

        const long long distance_on_read =
            static_cast<long long>(next.query_start_bp) - static_cast<long long>(curr.query_end_bp);
        if (distance_on_read < -kQueryOverlapTolerance) {
            continue;
        }

        long long distance_on_ref = 0;
        if (curr.reverse) {
            distance_on_ref =
                static_cast<long long>(curr.ref_start_bp) - static_cast<long long>(next.ref_end_bp);
        } else {
            distance_on_ref =
                static_cast<long long>(next.ref_start_bp) - static_cast<long long>(curr.ref_end_bp);
        }

        if (distance_on_ref >= -kRefOverlapTolerance) {
            const long long deviation = distance_on_read - distance_on_ref;
            if (deviation >= static_cast<long long>(min_sv_bp) &&
                distance_on_ref <= kRefGapTolerance) {
                AtomicVariantEvent ev;
                ev.event_type = "INS";
                const std::size_t ref_anchor = curr.reverse ? curr.ref_start_bp : curr.ref_end_bp;
                ev.ref_offset_start_bp = ref_anchor;
                ev.ref_offset_end_bp = ref_anchor;
                ev.alt_offset_start_bp = std::min<std::size_t>(curr.query_end_bp, query_len);
                ev.alt_offset_end_bp = std::min<std::size_t>(query_len, ev.alt_offset_start_bp + static_cast<std::size_t>(deviation));
                ev.event_bp = static_cast<std::size_t>(deviation);
                ev.inserted_bp = ev.event_bp;
                ev.svlen = static_cast<long long>(ev.event_bp);
                ev.preserve_ins_svlen = split_ins_use_geometric_svlen;
                ev.evidence_source = "split";
                out.push_back(std::move(ev));
            } else if (deviation <= -static_cast<long long>(min_sv_bp) &&
                       deviation >= -static_cast<long long>(max_sv_bp) &&
                       distance_on_read <= kQueryGapTolerance) {
                AtomicVariantEvent ev;
                ev.event_type = "DEL";
                const std::size_t del_len = static_cast<std::size_t>(-deviation);
                const std::size_t ref_anchor = curr.reverse ? next.ref_end_bp : curr.ref_end_bp;
                ev.ref_offset_start_bp = ref_anchor;
                ev.ref_offset_end_bp = ref_anchor + del_len;
                ev.alt_offset_start_bp = curr.query_end_bp;
                ev.alt_offset_end_bp = curr.query_end_bp;
                ev.event_bp = del_len;
                ev.svlen = -static_cast<long long>(del_len);
                ev.cn_delta = -1;
                ev.evidence_source = "split";
                out.push_back(std::move(ev));
            }
        } else {
            if (distance_on_read > kQueryGapTolerance) {
                continue;
            }
            const long long deviation = distance_on_read - distance_on_ref;
            if (deviation < static_cast<long long>(min_sv_bp)) {
                continue;
            }
            AtomicVariantEvent ev;
            ev.event_type = "INS";
            const std::size_t ref_anchor = curr.reverse ? curr.ref_start_bp : next.ref_start_bp;
            ev.ref_offset_start_bp = ref_anchor;
            ev.ref_offset_end_bp = ref_anchor;
            ev.alt_offset_start_bp = std::min<std::size_t>(curr.query_end_bp, query_len);
            ev.alt_offset_end_bp = std::min<std::size_t>(query_len, ev.alt_offset_start_bp + static_cast<std::size_t>(deviation));
            ev.event_bp = static_cast<std::size_t>(deviation);
            ev.inserted_bp = ev.event_bp;
            ev.svlen = static_cast<long long>(ev.event_bp);
            ev.preserve_ins_svlen = split_ins_use_geometric_svlen;
            ev.evidence_source = "split";
            out.push_back(std::move(ev));
        }
    }
    return out;
}

struct DupSearchResult {
    bool found = false;
    double best_similarity = 0.0;
    int best_edit_distance = std::numeric_limits<int>::max();
    std::size_t ref_start_bp = 0;
    std::size_t ref_end_bp = 0;
    bool reverse = false;
};

std::vector<std::size_t> candidate_dup_starts_from_kmers(
    const std::string& query,
    const std::string& reference,
    std::size_t k,
    std::size_t max_candidates) {

    if (query.size() < k || reference.size() < k || k == 0) {
        return {};
    }

    std::unordered_map<std::string, std::vector<std::size_t>> ref_hits;
    ref_hits.reserve(reference.size() * 2);
    constexpr std::size_t kMaxRefHitsPerKmer = 96;
    for (std::size_t i = 0; i + k <= reference.size(); ++i) {
        const std::string key = reference.substr(i, k);
        auto& hits = ref_hits[key];
        if (hits.size() < kMaxRefHitsPerKmer) {
            hits.push_back(i);
        }
    }

    std::unordered_map<std::size_t, std::size_t> votes;
    votes.reserve(query.size() * 4);
    const std::size_t q_step = std::max<std::size_t>(1, query.size() / 64);
    for (std::size_t q = 0; q + k <= query.size(); q += q_step) {
        const std::string key = query.substr(q, k);
        const auto it = ref_hits.find(key);
        if (it == ref_hits.end()) {
            continue;
        }
        for (const auto ref_pos : it->second) {
            const std::size_t raw_start = (ref_pos >= q) ? (ref_pos - q) : 0;
            const std::size_t max_start = reference.size() > query.size() ? (reference.size() - query.size()) : 0;
            const std::size_t start = std::min(raw_start, max_start);
            votes[start] += 1;
        }
    }

    std::vector<std::pair<std::size_t, std::size_t>> ranked(votes.begin(), votes.end());
    std::sort(ranked.begin(), ranked.end(), [](const auto& lhs, const auto& rhs) {
        if (lhs.second != rhs.second) {
            return lhs.second > rhs.second;
        }
        return lhs.first < rhs.first;
    });

    std::vector<std::size_t> starts;
    starts.reserve(std::min(max_candidates, ranked.size()));
    for (std::size_t i = 0; i < ranked.size() && i < max_candidates; ++i) {
        starts.push_back(ranked[i].first);
    }
    return starts;
}

std::string canonical_kmer_window(
    const std::string& sequence,
    std::size_t start,
    std::size_t k) {

    if (k == 0 || start + k > sequence.size()) {
        return {};
    }

    std::string kmer;
    kmer.reserve(k);
    for (std::size_t i = 0; i < k; ++i) {
        const char up = static_cast<char>(std::toupper(static_cast<unsigned char>(sequence[start + i])));
        if (up != 'A' && up != 'C' && up != 'G' && up != 'T') {
            return {};
        }
        kmer.push_back(up);
    }
    const std::string rc = reverse_complement(kmer);
    return (rc < kmer) ? rc : kmer;
}

std::size_t count_exact_non_overlapping_matches(
    const std::string& query,
    const std::string& target) {

    if (query.empty() || target.empty() || query.size() > target.size()) {
        return 0;
    }

    std::vector<std::pair<std::size_t, std::size_t>> intervals;
    intervals.reserve(16);
    auto collect_orientation = [&](const std::string& pattern) {
        std::size_t pos = target.find(pattern);
        while (pos != std::string::npos) {
            intervals.push_back({pos, pos + pattern.size()});
            pos = target.find(pattern, pos + 1);
        }
    };

    collect_orientation(query);
    const std::string query_rc = reverse_complement(query);
    if (query_rc != query) {
        collect_orientation(query_rc);
    }
    if (intervals.empty()) {
        return 0;
    }

    std::sort(intervals.begin(), intervals.end(), [](const auto& lhs, const auto& rhs) {
        if (lhs.first != rhs.first) {
            return lhs.first < rhs.first;
        }
        return lhs.second < rhs.second;
    });

    std::size_t count = 0;
    std::size_t current_end = 0;
    bool have_group = false;
    for (const auto& [start, end] : intervals) {
        if (!have_group || start >= current_end) {
            ++count;
            current_end = end;
            have_group = true;
        } else if (end > current_end) {
            current_end = end;
        }
    }
    return count;
}

struct DupCopyCountEstimate {
    bool valid = false;
    std::size_t copy_count = 0;
    std::size_t exact_hits = 0;
    std::size_t informative_kmers = 0;
    double nonzero_fraction = 0.0;
};

struct DupRefinementDebug {
    bool attempted = false;
    bool dup_match_found = false;
    bool similarity_pass = false;
    double min_similarity = 0.0;
    double best_similarity = 0.0;
    std::size_t ref_start_bp = 0;
    std::size_t ref_end_bp = 0;
    char orientation = '+';
    std::size_t dup_unit_bp = 0;
    bool has_valid_window = false;
    DupCopyCountEstimate ref_est;
    DupCopyCountEstimate alt_est;
    std::size_t added_by_length = 0;
    std::size_t raw_added = 0;
    std::string decision = "skip";
};

[[maybe_unused]] DupCopyCountEstimate estimate_dup_copy_count(
    const std::string& unit_seq,
    const std::string& target_seq) {

    DupCopyCountEstimate out;
    out.exact_hits = count_exact_non_overlapping_matches(unit_seq, target_seq);

    if (unit_seq.empty() || target_seq.empty()) {
        out.valid = true;
        out.copy_count = 0;
        return out;
    }

    std::size_t copy_count = out.exact_hits;
    const std::size_t k = std::clamp<std::size_t>(unit_seq.size() / 8, 9, 21);
    if (unit_seq.size() >= k && target_seq.size() >= k) {
        std::unordered_map<std::string, std::size_t> target_kmer_counts;
        target_kmer_counts.reserve(target_seq.size());
        for (std::size_t i = 0; i + k <= target_seq.size(); ++i) {
            const std::string key = canonical_kmer_window(target_seq, i, k);
            if (key.empty()) {
                continue;
            }
            target_kmer_counts[key] += 1;
        }

        std::unordered_set<std::string> seen_unit_kmers;
        seen_unit_kmers.reserve(unit_seq.size());
        std::vector<std::size_t> all_counts;
        std::vector<std::size_t> nonzero_counts;
        all_counts.reserve(unit_seq.size());
        nonzero_counts.reserve(unit_seq.size());
        constexpr std::size_t kMaxInformativeKmerCount = 64;

        for (std::size_t i = 0; i + k <= unit_seq.size(); ++i) {
            const std::string key = canonical_kmer_window(unit_seq, i, k);
            if (key.empty() || !seen_unit_kmers.insert(key).second) {
                continue;
            }
            const auto it = target_kmer_counts.find(key);
            const std::size_t count = (it == target_kmer_counts.end()) ? 0 : it->second;
            if (count > kMaxInformativeKmerCount) {
                continue;
            }
            all_counts.push_back(count);
            if (count > 0) {
                nonzero_counts.push_back(count);
            }
        }

        out.informative_kmers = all_counts.size();
        if (!all_counts.empty()) {
            out.nonzero_fraction =
                static_cast<double>(nonzero_counts.size()) /
                static_cast<double>(all_counts.size());
        }

        if (all_counts.size() >= 8) {
            if (nonzero_counts.empty() || (out.nonzero_fraction < 0.15 && out.exact_hits == 0)) {
                copy_count = 0;
            } else {
                std::sort(nonzero_counts.begin(), nonzero_counts.end());
                const std::size_t q_idx = static_cast<std::size_t>(
                    std::floor(0.80 * static_cast<double>(nonzero_counts.size() - 1)));
                const std::size_t q80 = nonzero_counts[std::min(q_idx, nonzero_counts.size() - 1)];
                copy_count = std::max(copy_count, static_cast<std::size_t>(std::llround(static_cast<double>(q80))));
            }
        }
    }

    out.valid = true;
    out.copy_count = copy_count;
    return out;
}

[[maybe_unused]] DupSearchResult classify_dup_like_insert(
    const std::string& inserted_seq,
    const std::string& reference_seq,
    double min_similarity) {

    DupSearchResult out;
    if (inserted_seq.empty() || reference_seq.empty()) {
        return out;
    }

    auto update_best = [&](double sim, int d, bool reverse, std::size_t start, std::size_t win_len) {
        if (!out.found ||
            sim > out.best_similarity + 1e-12 ||
            (std::abs(sim - out.best_similarity) <= 1e-12 && d < out.best_edit_distance)) {
            out.found = true;
            out.best_similarity = sim;
            out.best_edit_distance = d;
            out.ref_start_bp = start;
            out.ref_end_bp = start + win_len;
            out.reverse = reverse;
        }
    };

    auto evaluate = [&](const std::string& query,
                        const std::vector<std::uint64_t>& query_sketch,
                        bool use_sketch_eval,
                        bool reverse,
                        std::size_t start,
                        std::size_t win_len,
                        double approx_min_similarity) {
        if (win_len == 0 || start >= reference_seq.size() || start + win_len > reference_seq.size()) {
            return;
        }

        const std::size_t max_len = std::max(query.size(), win_len);
        if (use_sketch_eval) {
            if (query_sketch.empty()) {
                return;
            }
            constexpr std::size_t kSketchK = 15;
            constexpr std::size_t kSketchSize = 256;
            const std::string window = reference_seq.substr(start, win_len);
            const auto window_sketch = build_sequence_minhash_sketch(window, kSketchK, kSketchSize);
            if (window_sketch.empty()) {
                return;
            }
            const double jacc = sketch_jaccard(query_sketch, window_sketch);
            const double sim = estimate_identity_from_jaccard(jacc);
            if (sim + 1e-12 < approx_min_similarity) {
                return;
            }
            const int d_est = static_cast<int>(std::llround(
                (1.0 - sim) * static_cast<double>(std::max<std::size_t>(1, max_len))));
            update_best(sim, d_est, reverse, start, win_len);
            return;
        }

        const std::string window = reference_seq.substr(start, win_len);
        const int allowed = allowed_max_edits(min_similarity, max_len);
        const int d = bounded_levenshtein_distance(query, window, allowed);
        if (d > allowed) {
            return;
        }
        const double sim = 1.0 - (static_cast<double>(d) / static_cast<double>(std::max<std::size_t>(1, max_len)));
        update_best(sim, d, reverse, start, win_len);
    };

    auto search_orientation = [&](const std::string& query, bool reverse) {
        const std::size_t qlen = query.size();
        constexpr std::size_t kMaxDpQueryBp = 5000;
        constexpr std::size_t kSmallSeedCandidates = 128;
        constexpr std::size_t kLargeSeedCandidates = 24;
        constexpr double kLargeApproxSlack = 0.05;

        const bool use_sketch_eval =
            (qlen > kMaxDpQueryBp || qlen > (reference_seq.size() / 2));
        const std::size_t max_candidates = use_sketch_eval ? kLargeSeedCandidates : kSmallSeedCandidates;
        const double approx_min_similarity = std::clamp(min_similarity - kLargeApproxSlack, 0.0, 1.0);

        std::vector<std::uint64_t> query_sketch;
        if (use_sketch_eval) {
            constexpr std::size_t kSketchK = 15;
            constexpr std::size_t kSketchSize = 256;
            query_sketch = build_sequence_minhash_sketch(query, kSketchK, kSketchSize);
        }

        const std::size_t k = std::clamp<std::size_t>(qlen / 10, 9, 17);
        std::vector<std::size_t> starts = candidate_dup_starts_from_kmers(query, reference_seq, k, max_candidates);

        if (use_sketch_eval && starts.empty()) {
            if (reference_seq.size() > qlen) {
                const std::size_t max_start = reference_seq.size() - qlen;
                const std::size_t coarse = 12;
                const std::size_t step = std::max<std::size_t>(1, max_start / coarse);
                for (std::size_t s = 0; s <= max_start && starts.size() < coarse; s += step) {
                    starts.push_back(s);
                }
                if (starts.empty() || starts.back() != max_start) {
                    starts.push_back(max_start);
                }
            } else {
                starts.push_back(0);
            }
        }

        const std::vector<double> len_scales = use_sketch_eval
            ? std::vector<double>{0.9, 1.0, 1.1}
            : std::vector<double>{0.8, 0.9, 1.0, 1.1, 1.2};

        const std::size_t jitter = use_sketch_eval
            ? std::max<std::size_t>(8, qlen / 40)
            : std::max<std::size_t>(4, qlen / 20);
        const std::vector<long long> shifts = use_sketch_eval
            ? std::vector<long long>{
                -static_cast<long long>(jitter),
                0,
                static_cast<long long>(jitter)}
            : std::vector<long long>{
                -static_cast<long long>(2 * jitter),
                -static_cast<long long>(jitter),
                0,
                static_cast<long long>(jitter),
                static_cast<long long>(2 * jitter)};

        for (const std::size_t base_start : starts) {
            for (const auto scale : len_scales) {
                std::size_t win_len = static_cast<std::size_t>(std::llround(scale * static_cast<double>(qlen)));
                win_len = std::max<std::size_t>(1, win_len);
                if (win_len > reference_seq.size()) {
                    continue;
                }
                const std::size_t max_start = reference_seq.size() - win_len;
                for (const auto shift : shifts) {
                    const long long shifted = static_cast<long long>(base_start) + shift;
                    const std::size_t start = (shifted < 0)
                                                  ? 0
                                                  : std::min<std::size_t>(static_cast<std::size_t>(shifted), max_start);
                    evaluate(query, query_sketch, use_sketch_eval, reverse, start, win_len, approx_min_similarity);
                }
            }
        }

        // Fallback if k-mer seeds miss highly diverged duplicates.
        if (!out.found && qlen <= reference_seq.size()) {
            const std::size_t step = use_sketch_eval
                ? std::max<std::size_t>(1, (reference_seq.size() - qlen) / 12)
                : std::max<std::size_t>(1, qlen / 10);
            const std::size_t win_len = qlen;
            for (std::size_t start = 0; start + win_len <= reference_seq.size(); start += step) {
                evaluate(query, query_sketch, use_sketch_eval, reverse, start, win_len, approx_min_similarity);
            }
            const std::size_t last_start = reference_seq.size() - win_len;
            evaluate(query, query_sketch, use_sketch_eval, reverse, last_start, win_len, approx_min_similarity);
        }
    };

    // Keep fast exact fallback first.
    const std::size_t exact_fwd = reference_seq.find(inserted_seq);
    if (exact_fwd != std::string::npos) {
        out.found = true;
        out.best_similarity = 1.0;
        out.best_edit_distance = 0;
        out.ref_start_bp = exact_fwd;
        out.ref_end_bp = exact_fwd + inserted_seq.size();
        out.reverse = false;
        return out;
    }
    const std::string inserted_rc = reverse_complement(inserted_seq);
    const std::size_t exact_rev = reference_seq.find(inserted_rc);
    if (exact_rev != std::string::npos) {
        out.found = true;
        out.best_similarity = 1.0;
        out.best_edit_distance = 0;
        out.ref_start_bp = exact_rev;
        out.ref_end_bp = exact_rev + inserted_seq.size();
        out.reverse = true;
        return out;
    }

    search_orientation(inserted_seq, false);
    search_orientation(inserted_rc, true);
    return out;
}
std::string inserted_sequence_for_event(
    const AtomicVariantEvent& ev,
    const std::string& alt_seq) {

    if (!ev.inserted_seq.empty()) {
        return ev.inserted_seq;
    }
    if (ev.alt_offset_end_bp > ev.alt_offset_start_bp &&
        ev.alt_offset_end_bp <= alt_seq.size()) {
        return alt_seq.substr(ev.alt_offset_start_bp, ev.alt_offset_end_bp - ev.alt_offset_start_bp);
    }
    return {};
}

void classify_insertion_events(
    std::vector<AtomicVariantEvent>& events,
    const std::string& ref_seq,
    const std::string& alt_seq,
    bool classify_ins) {

    constexpr double kDupMinSimilarity = 0.75;
    constexpr std::size_t kTandemWindowBp = 1000;
    for (auto& ev : events) {
        if (ev.event_type != "INS") {
            continue;
        }
        if (!classify_ins) {
            continue;
        }
        ev.event_subtype = "NOVEL";

        const std::string inserted_seq = inserted_sequence_for_event(ev, alt_seq);
        if (inserted_seq.empty()) {
            continue;
        }
        ev.inserted_bp = std::max(ev.inserted_bp, inserted_seq.size());

        const DupSearchResult dup = classify_dup_like_insert(
            inserted_seq,
            ref_seq,
            kDupMinSimilarity);
        if (!dup.found || dup.ref_end_bp <= dup.ref_start_bp) {
            continue;
        }

        ev.has_dup_evidence = true;
        ev.dup_best_similarity = dup.best_similarity;
        ev.dup_ref_start_bp = dup.ref_start_bp;
        ev.dup_ref_end_bp = dup.ref_end_bp;
        ev.dup_orientation = dup.reverse ? '-' : '+';
        ev.dup_unit_bp = std::max<std::size_t>(1, dup.ref_end_bp - dup.ref_start_bp);

        const std::size_t anchor = std::min(ev.ref_offset_start_bp, ref_seq.size());
        const std::size_t local_lo = (anchor > kTandemWindowBp) ? (anchor - kTandemWindowBp) : 0;
        const std::size_t local_hi = std::min<std::size_t>(ref_seq.size(), anchor + kTandemWindowBp);
        const bool local_dup =
            (ev.dup_ref_start_bp < local_hi) && (ev.dup_ref_end_bp > local_lo);
        if (local_dup) {
            ev.event_subtype = dup.reverse ? "DUP_TANDEM_INV" : "DUP_TANDEM";
        } else {
            ev.event_subtype = dup.reverse ? "DUP_INTERSPERSED_INV" : "DUP_INTERSPERSED";
        }

        const std::size_t unit_start = std::min(ev.dup_ref_start_bp, ref_seq.size());
        const std::size_t unit_end = std::min(ev.dup_ref_end_bp, ref_seq.size());
        std::string unit_seq;
        if (unit_end > unit_start) {
            unit_seq = ref_seq.substr(unit_start, unit_end - unit_start);
        }

        std::size_t ref_cn = 1;
        std::size_t alt_cn_obs = 0;
        if (!unit_seq.empty()) {
            const DupCopyCountEstimate ref_est = estimate_dup_copy_count(unit_seq, ref_seq);
            const DupCopyCountEstimate alt_est = estimate_dup_copy_count(unit_seq, alt_seq);
            if (ref_est.valid && ref_est.copy_count > 0) {
                ref_cn = ref_est.copy_count;
            }
            if (alt_est.valid) {
                alt_cn_obs = alt_est.copy_count;
            }
        }

        const std::size_t added_by_length = std::max<std::size_t>(
            1,
            static_cast<std::size_t>(std::llround(
                static_cast<double>(inserted_seq.size()) /
                static_cast<double>(std::max<std::size_t>(1, ev.dup_unit_bp)))));
        std::size_t alt_cn = (alt_cn_obs >= ref_cn) ? alt_cn_obs : (ref_cn + added_by_length);
        alt_cn = std::max<std::size_t>(alt_cn, ref_cn + added_by_length);

        ev.dup_ref_copy_number = ref_cn;
        ev.dup_alt_copy_number = alt_cn;
        ev.dup_added_copies = (alt_cn > ref_cn) ? (alt_cn - ref_cn) : 0;
        ev.dup_copy_ratio = (ref_cn > 0)
            ? (static_cast<double>(alt_cn) / static_cast<double>(ref_cn))
            : 0.0;
        ev.cn_delta = static_cast<int>(ev.dup_added_copies);
    }
}

[[maybe_unused]] std::vector<DotplotPoint> collect_dotplot_points(
    const std::string& ref,
    const std::string& query,
    std::size_t k,
    std::size_t max_points) {

    if (k == 0 || ref.size() < k || query.size() < k || max_points == 0) {
        return {};
    }

    std::unordered_map<std::string, std::vector<std::size_t>> ref_kmers;
    ref_kmers.reserve(ref.size());
    constexpr std::size_t kMaxRefHitsPerKmer = 256;
    for (std::size_t i = 0; i + k <= ref.size(); ++i) {
        const std::string key = ref.substr(i, k);
        auto& hits = ref_kmers[key];
        if (hits.size() < kMaxRefHitsPerKmer) {
            hits.push_back(i);
        }
    }

    std::vector<DotplotPoint> points;
    points.reserve(max_points);
    const std::size_t cap_fwd = (max_points + 1) / 2;
    const std::size_t cap_rev = max_points - cap_fwd;

    std::size_t fwd_count = 0;
    for (std::size_t q = 0; q + k <= query.size() && fwd_count < cap_fwd; ++q) {
        const std::string key = query.substr(q, k);
        const auto it = ref_kmers.find(key);
        if (it == ref_kmers.end()) {
            continue;
        }
        for (const auto r : it->second) {
            if (fwd_count >= cap_fwd) {
                break;
            }
            points.push_back({r, q, 'F'});
            ++fwd_count;
        }
    }

    const std::string query_rc = reverse_complement(query);
    std::size_t rev_count = 0;
    for (std::size_t q = 0; q + k <= query_rc.size() && rev_count < cap_rev; ++q) {
        const std::string key = query_rc.substr(q, k);
        const auto it = ref_kmers.find(key);
        if (it == ref_kmers.end()) {
            continue;
        }
        const std::size_t mapped_q = query.size() - q - k;
        for (const auto r : it->second) {
            if (rev_count >= cap_rev) {
                break;
            }
            points.push_back({r, mapped_q, 'R'});
            ++rev_count;
        }
    }

    return points;
}

[[maybe_unused]] void write_dotplot_matches_tsv(
    const std::string& output_path,
    const std::vector<DotplotPoint>& points) {

    std::ofstream out(output_path);
    if (!out) {
        throw std::runtime_error("Failed to write dotplot matches TSV: " + output_path);
    }
    out << "ref_pos\tquery_pos\tstrand\n";
    for (const auto& p : points) {
        out << p.ref_pos << '\t' << p.query_pos << '\t' << p.strand << '\n';
    }
}

[[maybe_unused]] void write_dotplot_svg(
    const std::string& output_path,
    std::size_t bubble_id,
    std::size_t cluster_id,
    const std::string& event_type,
    std::size_t k,
    const std::string& ref_seq,
    const std::string& query_seq,
    const std::vector<DotplotPoint>& points) {

    const double width = 940.0;
    const double height = 940.0;
    const double margin = 80.0;
    const double plot_w = width - (2.0 * margin);
    const double plot_h = height - (2.0 * margin);
    const double ref_den = std::max(1.0, static_cast<double>(ref_seq.size() - (ref_seq.empty() ? 0 : 1)));
    const double qry_den = std::max(1.0, static_cast<double>(query_seq.size() - (query_seq.empty() ? 0 : 1)));

    std::size_t forward_points = 0;
    std::size_t reverse_points = 0;
    for (const auto& p : points) {
        if (p.strand == 'R') {
            ++reverse_points;
        } else {
            ++forward_points;
        }
    }

    std::ofstream out(output_path);
    if (!out) {
        throw std::runtime_error("Failed to write dotplot SVG: " + output_path);
    }

    out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    out << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << width
        << "\" height=\"" << height << "\" viewBox=\"0 0 " << width << " " << height << "\">\n";
    out << "<rect x=\"0\" y=\"0\" width=\"" << width << "\" height=\"" << height << "\" fill=\"white\"/>\n";
    out << "<text x=\"" << margin << "\" y=\"28\" font-family=\"Arial, sans-serif\" font-size=\"18\" font-weight=\"bold\">"
        << "Bubble " << bubble_id << " Cluster " << cluster_id << " Dotplot</text>\n";
    out << "<text x=\"" << margin << "\" y=\"48\" font-family=\"Arial, sans-serif\" font-size=\"12\" fill=\"#555\">"
        << "event=" << html_escape(event_type) << ", k=" << k
        << ", forward points=" << forward_points << ", reverse points=" << reverse_points
        << "</text>\n";

    out << "<rect x=\"" << margin << "\" y=\"" << margin << "\" width=\"" << plot_w << "\" height=\"" << plot_h
        << "\" fill=\"none\" stroke=\"#DDDDDD\"/>\n";
    out << "<line x1=\"" << margin << "\" y1=\"" << margin << "\" x2=\"" << (margin + plot_w)
        << "\" y2=\"" << (margin + plot_h) << "\" stroke=\"#CCCCCC\" stroke-dasharray=\"4,4\"/>\n";
    out << "<line x1=\"" << (margin + plot_w) << "\" y1=\"" << margin << "\" x2=\"" << margin
        << "\" y2=\"" << (margin + plot_h) << "\" stroke=\"#CCCCCC\" stroke-dasharray=\"4,4\"/>\n";

    for (const auto& p : points) {
        const double x = margin + (static_cast<double>(p.ref_pos) / ref_den) * plot_w;
        const double y = margin + (static_cast<double>(p.query_pos) / qry_den) * plot_h;
        const char* color = (p.strand == 'R') ? "#E41A1C" : "#377EB8";
        out << "<circle cx=\"" << x << "\" cy=\"" << y << "\" r=\"1.2\" fill=\"" << color
            << "\" fill-opacity=\"0.75\"/>\n";
    }

    out << "<text x=\"" << (margin + (plot_w / 2.0)) << "\" y=\"" << (height - 24.0)
        << "\" text-anchor=\"middle\" font-family=\"Arial, sans-serif\" font-size=\"12\">Reference sequence position (bp)</text>\n";
    out << "<text x=\"24\" y=\"" << (margin + (plot_h / 2.0))
        << "\" transform=\"rotate(-90,24," << (margin + (plot_h / 2.0))
        << ")\" text-anchor=\"middle\" font-family=\"Arial, sans-serif\" font-size=\"12\">Cluster sequence position (bp)</text>\n";

    out << "<rect x=\"" << (width - 250.0) << "\" y=\"58\" width=\"10\" height=\"10\" fill=\"#377EB8\"/>\n";
    out << "<text x=\"" << (width - 236.0) << "\" y=\"67\" font-family=\"Arial, sans-serif\" font-size=\"10\">forward matches</text>\n";
    out << "<rect x=\"" << (width - 250.0) << "\" y=\"74\" width=\"10\" height=\"10\" fill=\"#E41A1C\"/>\n";
    out << "<text x=\"" << (width - 236.0) << "\" y=\"83\" font-family=\"Arial, sans-serif\" font-size=\"10\">reverse-complement matches</text>\n";

    out << "</svg>\n";
}

struct DotplotGeneBand {
    std::size_t ref_start_bp = 0;
    std::size_t ref_end_bp = 0;
    std::string name;
};

struct DotplotMatchSegment {
    std::size_t ref_start_bp = 0;
    std::size_t ref_end_bp = 0;
    double query_start_bp = 0.0;
    double query_end_bp = 0.0;
};

ParsedReferencePath parse_reference_path_label(const std::string& path_name);
std::vector<CigarRun> parse_extended_cigar(const std::string& cigar);

std::vector<std::string> split_tab_fields(const std::string& line) {
    std::vector<std::string> out;
    std::string field;
    std::istringstream iss(line);
    while (std::getline(iss, field, '\t')) {
        out.push_back(field);
    }
    return out;
}

std::string trim_ascii_whitespace(const std::string& text) {
    std::size_t lo = 0;
    while (lo < text.size() && std::isspace(static_cast<unsigned char>(text[lo]))) {
        ++lo;
    }
    std::size_t hi = text.size();
    while (hi > lo && std::isspace(static_cast<unsigned char>(text[hi - 1]))) {
        --hi;
    }
    return text.substr(lo, hi - lo);
}

std::string gtf_attribute_value(const std::string& attributes, const std::string& key) {
    const std::string needle = key + " ";
    std::size_t pos = attributes.find(needle);
    if (pos == std::string::npos) {
        return {};
    }
    pos += needle.size();
    while (pos < attributes.size() && std::isspace(static_cast<unsigned char>(attributes[pos]))) {
        ++pos;
    }
    if (pos >= attributes.size()) {
        return {};
    }
    if (attributes[pos] == '"') {
        const std::size_t end_quote = attributes.find('"', pos + 1);
        if (end_quote == std::string::npos) {
            return {};
        }
        return attributes.substr(pos + 1, end_quote - pos - 1);
    }
    const std::size_t end = attributes.find(';', pos);
    return trim_ascii_whitespace(attributes.substr(pos, end == std::string::npos ? std::string::npos : (end - pos)));
}

std::string fallback_gene_name(const std::string& raw) {
    if (raw.empty()) {
        return ".";
    }
    return raw;
}

std::string normalize_chrom_for_lookup(const std::string& chrom) {
    std::string out = trim_ascii_whitespace(chrom);
    if (out.empty()) {
        return out;
    }
    if (out.rfind("chr", 0) == 0 || out.rfind("CHR", 0) == 0) {
        out[0] = 'c';
        out[1] = 'h';
        out[2] = 'r';
    }
    return out;
}

const std::vector<GeneAnnotation>* find_gene_annotations(
    const GeneAnnotationIndex& index,
    const std::string& chrom) {

    auto find_exact = [&](const std::string& key) -> const std::vector<GeneAnnotation>* {
        const auto it = index.find(key);
        if (it == index.end()) {
            return nullptr;
        }
        return &it->second;
    };

    const std::string norm = normalize_chrom_for_lookup(chrom);
    if (const auto* hit = find_exact(norm)) {
        return hit;
    }
    if (norm.rfind("chr", 0) == 0 && norm.size() > 3) {
        if (const auto* hit = find_exact(norm.substr(3))) {
            return hit;
        }
    } else if (!norm.empty()) {
        if (const auto* hit = find_exact("chr" + norm)) {
            return hit;
        }
    }
    return nullptr;
}

GeneAnnotationIndex load_gene_annotations_from_gtf(const std::string& gtf_path) {
    GeneAnnotationIndex index;
    if (gtf_path.empty()) {
        return index;
    }

    auto ingest_line = [&](const std::string& line) {
        if (line.empty() || line[0] == '#') {
            return;
        }
        const auto fields = split_tab_fields(line);
        if (fields.size() < 9) {
            return;
        }
        const std::string feature = fields[2];
        if (feature != "gene") {
            return;
        }

        std::size_t start1 = 0;
        std::size_t end1 = 0;
        try {
            start1 = static_cast<std::size_t>(std::stoull(fields[3]));
            end1 = static_cast<std::size_t>(std::stoull(fields[4]));
        } catch (const std::exception&) {
            return;
        }
        if (start1 == 0 || end1 < start1) {
            return;
        }

        std::string gene_name = gtf_attribute_value(fields[8], "gene_name");
        if (gene_name.empty()) {
            gene_name = gtf_attribute_value(fields[8], "gene_id");
        }
        gene_name = fallback_gene_name(gene_name);

        const std::string chrom = normalize_chrom_for_lookup(fields[0]);
        GeneAnnotation item;
        item.start0 = start1 - 1;
        item.end0 = end1;
        item.name = std::move(gene_name);
        index[chrom].push_back(std::move(item));
    };

    auto has_gzip_suffix = [](const std::string& path) {
        if (path.size() < 3) {
            return false;
        }
        const std::string suffix = path.substr(path.size() - 3);
        return suffix == ".gz" || suffix == ".GZ";
    };

    if (has_gzip_suffix(gtf_path)) {
        gzFile gz = gzopen(gtf_path.c_str(), "rb");
        if (gz == nullptr) {
            throw std::runtime_error("Failed to read gzipped GTF: " + gtf_path);
        }

        std::string chunk(1 << 15, '\0');
        std::string line;
        while (gzgets(gz, chunk.data(), static_cast<int>(chunk.size())) != nullptr) {
            line.assign(chunk.c_str());
            while (!line.empty() && line.back() != '\n' && !gzeof(gz)) {
                if (gzgets(gz, chunk.data(), static_cast<int>(chunk.size())) == nullptr) {
                    break;
                }
                line.append(chunk.c_str());
            }
            if (!line.empty() && line.back() == '\n') {
                line.pop_back();
            }
            ingest_line(line);
        }

        int gz_err = Z_OK;
        const char* gz_msg = gzerror(gz, &gz_err);
        gzclose(gz);
        if (gz_err != Z_OK && gz_err != Z_STREAM_END) {
            const std::string err_text = (gz_msg == nullptr) ? "unknown zlib error" : std::string(gz_msg);
            throw std::runtime_error("Failed while reading gzipped GTF '" + gtf_path + "': " + err_text);
        }
    } else {
        std::ifstream in(gtf_path);
        if (!in) {
            throw std::runtime_error("Failed to read GTF: " + gtf_path);
        }
        std::string line;
        while (std::getline(in, line)) {
            ingest_line(line);
        }
    }

    for (auto& [chrom, genes] : index) {
        std::sort(genes.begin(), genes.end(), [](const GeneAnnotation& lhs, const GeneAnnotation& rhs) {
            if (lhs.start0 != rhs.start0) {
                return lhs.start0 < rhs.start0;
            }
            if (lhs.end0 != rhs.end0) {
                return lhs.end0 < rhs.end0;
            }
            return lhs.name < rhs.name;
        });
        genes.erase(
            std::unique(genes.begin(), genes.end(), [](const GeneAnnotation& lhs, const GeneAnnotation& rhs) {
                return lhs.start0 == rhs.start0 && lhs.end0 == rhs.end0 && lhs.name == rhs.name;
            }),
            genes.end());
    }
    return index;
}

std::size_t query_forward_coord_from_record(const PafRecord& rec, std::size_t q_oriented) {
    if (rec.strand == '-') {
        if (rec.query_end >= q_oriented) {
            return rec.query_end - q_oriented;
        }
        return rec.query_start;
    }
    return rec.query_start + q_oriented;
}

[[maybe_unused]] std::vector<DotplotMatchSegment> collect_match_segments_from_records(const std::vector<PafRecord>& records) {
    std::vector<DotplotMatchSegment> segments;
    segments.reserve(records.size() * 16);
    for (const auto& rec : records) {
        bool emitted = false;
        if (!rec.cigar.empty()) {
            const auto runs = parse_extended_cigar(rec.cigar);
            std::size_t q_oriented = 0;
            std::size_t t_oriented = 0;
            const std::size_t q_span = rec.query_end >= rec.query_start ? (rec.query_end - rec.query_start) : 0;
            const std::size_t t_span = rec.target_end >= rec.target_start ? (rec.target_end - rec.target_start) : 0;
            for (const auto& run : runs) {
                if (run.len == 0) {
                    continue;
                }
                if (run.op == '=' || run.op == 'M' || run.op == 'X') {
                    if (q_oriented + run.len > q_span || t_oriented + run.len > t_span) {
                        break;
                    }
                    const std::size_t q0 = query_forward_coord_from_record(rec, q_oriented);
                    const std::size_t q1 = query_forward_coord_from_record(rec, q_oriented + run.len);
                    const std::size_t t0 = rec.target_start + t_oriented;
                    const std::size_t t1 = rec.target_start + t_oriented + run.len;
                    segments.push_back({
                        t0,
                        t1,
                        static_cast<double>(q0),
                        static_cast<double>(q1)});
                    emitted = true;
                    q_oriented += run.len;
                    t_oriented += run.len;
                    continue;
                }
                if (run.op == 'I') {
                    if (q_oriented + run.len > q_span) {
                        break;
                    }
                    q_oriented += run.len;
                    continue;
                }
                if (run.op == 'D') {
                    if (t_oriented + run.len > t_span) {
                        break;
                    }
                    t_oriented += run.len;
                    continue;
                }
            }
        }

        if (!emitted && rec.target_end > rec.target_start) {
            const double q0 = static_cast<double>(rec.strand == '-' ? rec.query_end : rec.query_start);
            const double q1 = static_cast<double>(rec.strand == '-' ? rec.query_start : rec.query_end);
            segments.push_back({
                rec.target_start,
                rec.target_end,
                q0,
                q1});
        }
    }
    return segments;
}

std::vector<DotplotGeneBand> collect_dotplot_gene_bands(
    const GeneAnnotationIndex* gene_index,
    const std::string& reference_label,
    std::size_t ref_len,
    const std::vector<std::string>& gene_name_matches) {

    std::vector<DotplotGeneBand> out;
    if (gene_index == nullptr || gene_index->empty() || ref_len == 0 || gene_name_matches.empty()) {
        return out;
    }

    struct GeneNameMatcher {
        std::optional<std::regex> regex;
        std::string lower;
    };
    std::vector<GeneNameMatcher> matchers;
    matchers.reserve(gene_name_matches.size());
    for (const auto& raw_match : gene_name_matches) {
        if (raw_match.empty()) {
            continue;
        }
        GeneNameMatcher matcher;
        try {
            matcher.regex.emplace(raw_match, std::regex_constants::icase);
        } catch (const std::regex_error&) {
        }
        matcher.lower = raw_match;
        std::transform(
            matcher.lower.begin(),
            matcher.lower.end(),
            matcher.lower.begin(),
            [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
        matchers.push_back(std::move(matcher));
    }
    if (matchers.empty()) {
        return out;
    }

    const ParsedReferencePath ref_meta = parse_reference_path_label(reference_label);
    if (!ref_meta.has_interval) {
        return out;
    }
    const auto* genes = find_gene_annotations(*gene_index, ref_meta.chrom);
    if (genes == nullptr || genes->empty()) {
        return out;
    }

    const std::size_t region_start0 = ref_meta.region_start_1based - 1;
    const std::size_t region_end0 = region_start0 + ref_len;
    for (const auto& gene : *genes) {
        bool gene_name_ok = false;
        std::string gene_name_lc;
        bool have_gene_name_lc = false;
        for (const auto& matcher : matchers) {
            if (matcher.regex.has_value()) {
                if (std::regex_search(gene.name, *matcher.regex)) {
                    gene_name_ok = true;
                    break;
                }
                continue;
            }
            if (!have_gene_name_lc) {
                gene_name_lc = gene.name;
                std::transform(
                    gene_name_lc.begin(),
                    gene_name_lc.end(),
                    gene_name_lc.begin(),
                    [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
                have_gene_name_lc = true;
            }
            if (gene_name_lc.find(matcher.lower) != std::string::npos) {
                gene_name_ok = true;
                break;
            }
        }
        if (!gene_name_ok) {
            continue;
        }
        if (gene.end0 <= region_start0 || gene.start0 >= region_end0) {
            continue;
        }
        const std::size_t clipped_start = std::max(gene.start0, region_start0);
        const std::size_t clipped_end = std::min(gene.end0, region_end0);
        if (clipped_end <= clipped_start) {
            continue;
        }
        out.push_back({
            clipped_start - region_start0,
            clipped_end - region_start0,
            gene.name});
    }

    std::sort(out.begin(), out.end(), [](const DotplotGeneBand& lhs, const DotplotGeneBand& rhs) {
        if (lhs.ref_start_bp != rhs.ref_start_bp) {
            return lhs.ref_start_bp < rhs.ref_start_bp;
        }
        if (lhs.ref_end_bp != rhs.ref_end_bp) {
            return lhs.ref_end_bp < rhs.ref_end_bp;
        }
        return lhs.name < rhs.name;
    });
    return out;
}

[[maybe_unused]] std::vector<std::pair<std::size_t, std::size_t>> project_gene_band_to_query(
    const DotplotGeneBand& band,
    const std::vector<DotplotMatchSegment>& segments,
    std::size_t query_len) {

    std::vector<std::pair<std::size_t, std::size_t>> spans;
    for (const auto& seg : segments) {
        if (seg.ref_end_bp <= band.ref_start_bp || seg.ref_start_bp >= band.ref_end_bp) {
            continue;
        }
        const std::size_t overlap_start = std::max(seg.ref_start_bp, band.ref_start_bp);
        const std::size_t overlap_end = std::min(seg.ref_end_bp, band.ref_end_bp);
        if (overlap_end <= overlap_start || seg.ref_end_bp <= seg.ref_start_bp) {
            continue;
        }
        const double seg_len = static_cast<double>(seg.ref_end_bp - seg.ref_start_bp);
        const double frac_lo = static_cast<double>(overlap_start - seg.ref_start_bp) / seg_len;
        const double frac_hi = static_cast<double>(overlap_end - seg.ref_start_bp) / seg_len;
        const double q0 = seg.query_start_bp + frac_lo * (seg.query_end_bp - seg.query_start_bp);
        const double q1 = seg.query_start_bp + frac_hi * (seg.query_end_bp - seg.query_start_bp);
        const double q_lo = std::max(0.0, std::min(q0, q1));
        const double q_hi = std::min(static_cast<double>(query_len), std::max(q0, q1));
        if (q_hi <= q_lo) {
            continue;
        }
        spans.emplace_back(
            static_cast<std::size_t>(std::floor(q_lo)),
            static_cast<std::size_t>(std::ceil(q_hi)));
    }

    if (spans.empty()) {
        return spans;
    }
    std::sort(spans.begin(), spans.end());
    std::vector<std::pair<std::size_t, std::size_t>> merged;
    merged.reserve(spans.size());
    for (const auto& span : spans) {
        if (merged.empty() || span.first > merged.back().second + 3) {
            merged.push_back(span);
            continue;
        }
        merged.back().second = std::max(merged.back().second, span.second);
    }
    return merged;
}

std::string cigar_from_runs(const std::vector<CigarRun>& runs) {
    std::ostringstream oss;
    for (const auto& run : runs) {
        if (run.len == 0) {
            continue;
        }
        oss << run.len << run.op;
    }
    return oss.str();
}

struct Minimap2MappingResult {
    MinimapPrimaryAlignment best;
    std::vector<PafRecord> records;
};

std::size_t clamp_i32_to_size(std::int32_t value, std::size_t cap) {
    if (value <= 0) {
        return 0;
    }
    return std::min<std::size_t>(static_cast<std::size_t>(value), cap);
}

std::vector<CigarRun> cigar_runs_from_minimap2(const mm_reg1_t& reg) {
    std::vector<CigarRun> runs;
    if (reg.p == nullptr || reg.p->n_cigar == 0) {
        return runs;
    }
    runs.reserve(reg.p->n_cigar);
    for (std::size_t i = 0; i < reg.p->n_cigar; ++i) {
        const std::uint32_t packed = reg.p->cigar[i];
        const std::size_t len = static_cast<std::size_t>(packed >> 4);
        const std::size_t op_idx = static_cast<std::size_t>(packed & 0x0fu);
        if (len == 0 || op_idx >= 10) {
            continue;
        }
        runs.push_back({MM_CIGAR_STR[op_idx], len});
    }
    return runs;
}

int edit_distance_from_runs(const std::vector<CigarRun>& runs) {
    std::size_t edit_bp = 0;
    for (const auto& run : runs) {
        if (run.op == 'X' || run.op == 'I' || run.op == 'D') {
            edit_bp += run.len;
        }
    }
    if (edit_bp > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        return std::numeric_limits<int>::max();
    }
    return static_cast<int>(edit_bp);
}

Minimap2MappingResult map_query_to_reference_minimap2(
    const std::string& query_name,
    const std::string& query_seq,
    const std::string& target_name,
    const std::string& target_seq,
    const std::string& minimap_preset,
    std::size_t minimap_best_n,
    bool minimap_emit_secondary) {

    Minimap2MappingResult out;
    if (query_seq.empty() || target_seq.empty()) {
        return out;
    }
    if (query_seq.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
        target_seq.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        return out;
    }

    mm_idxopt_t idx_opt;
    mm_mapopt_t map_opt;
    if (mm_set_opt(nullptr, &idx_opt, &map_opt) < 0) {
        return out;
    }
    const std::string preset = minimap_preset.empty() ? std::string("asm20") : minimap_preset;
    if (mm_set_opt(preset.c_str(), &idx_opt, &map_opt) < 0) {
        return out;
    }
    map_opt.flag |= (MM_F_CIGAR | MM_F_EQX | MM_F_OUT_CG);
    map_opt.flag &= ~MM_F_OUT_SAM;
    if (!minimap_emit_secondary) {
        map_opt.flag |= MM_F_NO_PRINT_2ND;
    }
    const int requested_best_n = static_cast<int>(
        std::min<std::size_t>(
            static_cast<std::size_t>(std::numeric_limits<int>::max()),
            std::max<std::size_t>(static_cast<std::size_t>(1), minimap_best_n)));
    map_opt.best_n = std::max(map_opt.best_n, requested_best_n);

    const char* seq_ptrs[1] = {target_seq.c_str()};
    const char* name_ptrs[1] = {target_name.c_str()};
    mm_idx_t* idx = mm_idx_str(
        idx_opt.w,
        idx_opt.k,
        (idx_opt.flag & MM_I_HPC) != 0 ? 1 : 0,
        idx_opt.bucket_bits,
        1,
        seq_ptrs,
        name_ptrs);
    if (idx == nullptr) {
        return out;
    }
    mm_mapopt_update(&map_opt, idx);

    mm_tbuf_t* tbuf = mm_tbuf_init();
    if (tbuf == nullptr) {
        mm_idx_destroy(idx);
        return out;
    }

    int n_regs = 0;
    mm_reg1_t* regs = mm_map(
        idx,
        static_cast<int>(query_seq.size()),
        query_seq.c_str(),
        &n_regs,
        tbuf,
        &map_opt,
        query_name.c_str());
    if (regs == nullptr || n_regs <= 0) {
        mm_tbuf_destroy(tbuf);
        mm_idx_destroy(idx);
        return out;
    }

    int best_idx = -1;
    for (int i = 0; i < n_regs; ++i) {
        if (best_idx < 0) {
            best_idx = i;
            continue;
        }
        const bool curr_primary = regs[i].parent == regs[i].id;
        const bool best_primary = regs[best_idx].parent == regs[best_idx].id;
        if (curr_primary != best_primary) {
            if (curr_primary) {
                best_idx = i;
            }
            continue;
        }
        if (regs[i].score > regs[best_idx].score ||
            (regs[i].score == regs[best_idx].score && regs[i].blen > regs[best_idx].blen) ||
            (regs[i].score == regs[best_idx].score &&
             regs[i].blen == regs[best_idx].blen &&
             regs[i].mlen > regs[best_idx].mlen)) {
            best_idx = i;
        }
    }

    out.records.reserve(static_cast<std::size_t>(n_regs));
    for (int i = 0; i < n_regs; ++i) {
        const mm_reg1_t& reg = regs[i];
        PafRecord rec;
        rec.query_name = query_name;
        rec.query_len = query_seq.size();
        rec.query_start = clamp_i32_to_size(reg.qs, query_seq.size());
        rec.query_end = clamp_i32_to_size(reg.qe, query_seq.size());
        if (rec.query_end < rec.query_start) {
            std::swap(rec.query_start, rec.query_end);
        }
        rec.strand = reg.rev ? '-' : '+';
        rec.target_name = target_name;
        rec.target_len = target_seq.size();
        rec.target_start = clamp_i32_to_size(reg.rs, target_seq.size());
        rec.target_end = clamp_i32_to_size(reg.re, target_seq.size());
        if (rec.target_end < rec.target_start) {
            std::swap(rec.target_start, rec.target_end);
        }
        rec.n_matches = reg.mlen > 0 ? static_cast<std::size_t>(reg.mlen) : 0;
        rec.aln_block_len = reg.blen > 0 ? static_cast<std::size_t>(reg.blen) : 0;
        rec.mapq = reg.mapq;
        rec.primary = (i == best_idx);
        rec.cm = reg.cnt;
        rec.s1 = reg.score;
        rec.s2 = reg.subsc;
        if (std::isfinite(reg.div) && reg.div >= 0.0f) {
            rec.dv = static_cast<double>(reg.div);
        }
        const auto runs = cigar_runs_from_minimap2(reg);
        rec.cigar = cigar_from_runs(runs);
        out.records.push_back(std::move(rec));

        if (i == best_idx) {
            out.best.ok = true;
            out.best.reverse = reg.rev;
            out.best.query_start_bp = out.records.back().query_start;
            out.best.query_end_bp = out.records.back().query_end;
            out.best.target_start_bp = out.records.back().target_start;
            out.best.target_end_bp = out.records.back().target_end;
            out.best.cigar_extended = out.records.back().cigar;
            out.best.edit_distance = edit_distance_from_runs(runs);
        }
    }

    for (int i = 0; i < n_regs; ++i) {
        std::free(regs[i].p);
    }
    std::free(regs);
    mm_tbuf_destroy(tbuf);
    mm_idx_destroy(idx);
    return out;
}

void write_paf_records(const std::string& output_path, const std::vector<PafRecord>& records) {
    std::ofstream out(output_path);
    if (!out) {
        throw std::runtime_error("Failed to write PAF: " + output_path);
    }
    for (const auto& rec : records) {
        out << rec.query_name << '\t'
            << rec.query_len << '\t'
            << rec.query_start << '\t'
            << rec.query_end << '\t'
            << rec.strand << '\t'
            << rec.target_name << '\t'
            << rec.target_len << '\t'
            << rec.target_start << '\t'
            << rec.target_end << '\t'
            << rec.n_matches << '\t'
            << rec.aln_block_len << '\t'
            << rec.mapq
            << "\ttp:A:" << (rec.primary ? "P" : "S");
        if (rec.cm >= 0) {
            out << "\tcm:i:" << rec.cm;
        }
        if (rec.s1 != std::numeric_limits<int>::min()) {
            out << "\ts1:i:" << rec.s1;
        }
        if (rec.s2 != std::numeric_limits<int>::min()) {
            out << "\ts2:i:" << rec.s2;
        }
        if (std::isfinite(rec.dv)) {
            std::ostringstream dv_ss;
            dv_ss << std::fixed << std::setprecision(4) << rec.dv;
            out << "\tdv:f:" << dv_ss.str();
        }
        if (!rec.cigar.empty()) {
            out << "\tcg:Z:" << rec.cigar;
        }
        out << '\n';
    }
}

void write_dotplot_svg_from_paf(
    const std::string& output_path,
    std::size_t bubble_id,
    std::size_t cluster_id,
    const std::string& ref_name,
    const std::string& query_name,
    std::size_t ref_len,
    std::size_t query_len,
    const std::vector<PafRecord>& records,
    const std::vector<AtomicVariantEvent>& called_events,
    const GeneAnnotationIndex* gene_index,
    const std::string& gtf_reference_label,
    const std::vector<std::string>& gene_name_matches) {

    const double width = 940.0;
    const double height = 1040.0;
    const double margin = 80.0;
    const double legend_h = 110.0;
    const double plot_w = width - (2.0 * margin);
    const double plot_h = height - (2.0 * margin) - legend_h;
    const double ref_den = std::max(1.0, static_cast<double>(std::max<std::size_t>(1, ref_len)));
    const double qry_den = std::max(1.0, static_cast<double>(std::max<std::size_t>(1, query_len)));
    const std::string ref_for_gtf = gtf_reference_label.empty() ? ref_name : gtf_reference_label;
    const auto gene_bands = collect_dotplot_gene_bands(gene_index, ref_for_gtf, ref_len, gene_name_matches);
    std::string gene_match_label;
    for (const auto& pattern : gene_name_matches) {
        if (pattern.empty()) {
            continue;
        }
        if (!gene_match_label.empty()) {
            gene_match_label += " | ";
        }
        gene_match_label += pattern;
    }

    constexpr const char* kGeneColors[] = {
        "#E31A1C",
        "#33A02C",
        "#1F78B4",
        "#FF7F00",
        "#6A3D9A",
        "#B15928",
        "#17A2B8",
        "#D95F02",
        "#7570B3",
        "#1B9E77",
    };
    constexpr std::size_t kGeneColorN = sizeof(kGeneColors) / sizeof(kGeneColors[0]);

    std::unordered_map<std::string, std::string> gene_color_by_name;
    std::vector<std::string> gene_legend_order;
    gene_color_by_name.reserve(gene_bands.size() * 2);
    gene_legend_order.reserve(gene_bands.size());
    for (const auto& band : gene_bands) {
        if (gene_color_by_name.find(band.name) != gene_color_by_name.end()) {
            continue;
        }
        const std::string color = kGeneColors[gene_legend_order.size() % kGeneColorN];
        gene_color_by_name.emplace(band.name, color);
        gene_legend_order.push_back(band.name);
    }

    std::ofstream out(output_path);
    if (!out) {
        throw std::runtime_error("Failed to write PAF dotplot SVG: " + output_path);
    }
    out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    out << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << width
        << "\" height=\"" << height << "\" viewBox=\"0 0 " << width << " " << height << "\">\n";
    out << "<rect x=\"0\" y=\"0\" width=\"" << width << "\" height=\"" << height << "\" fill=\"white\"/>\n";
    out << "<text x=\"" << margin << "\" y=\"28\" font-family=\"Arial, sans-serif\" font-size=\"18\" font-weight=\"bold\">"
        << "Bubble " << bubble_id << " Cluster " << cluster_id << " Dotplot (PAF)</text>\n";

    out << "<rect x=\"" << margin << "\" y=\"" << margin << "\" width=\"" << plot_w << "\" height=\"" << plot_h
        << "\" fill=\"none\" stroke=\"#DDDDDD\"/>\n";
    out << "<line x1=\"" << margin << "\" y1=\"" << (margin + plot_h) << "\" x2=\"" << (margin + plot_w)
        << "\" y2=\"" << margin << "\" stroke=\"#DDDDDD\" stroke-dasharray=\"4,4\"/>\n";

    auto plot_x = [&](std::size_t ref_pos) -> double {
        return margin + (static_cast<double>(ref_pos) / ref_den) * plot_w;
    };
    auto plot_xd = [&](double ref_pos) -> double {
        return margin + (ref_pos / ref_den) * plot_w;
    };
    auto plot_y = [&](std::size_t query_pos) -> double {
        return margin + plot_h - (static_cast<double>(query_pos) / qry_den) * plot_h;
    };
    auto plot_yd = [&](double query_pos) -> double {
        return margin + plot_h - (query_pos / qry_den) * plot_h;
    };
    auto format_bp = [](std::size_t v) -> std::string {
        const std::string s = std::to_string(v);
        std::string out;
        out.reserve(s.size() + (s.size() / 3));
        for (std::size_t i = 0; i < s.size(); ++i) {
            const std::size_t from_end = s.size() - i;
            if (i > 0 && (from_end % 3) == 0) {
                out.push_back(',');
            }
            out.push_back(s[i]);
        }
        return out;
    };
    auto clamp_ref_bp = [&](std::size_t pos) -> std::size_t {
        return std::min(pos, ref_len);
    };
    auto clamp_qry_bp = [&](std::size_t pos) -> std::size_t {
        return std::min(pos, query_len);
    };

    auto draw_gene_overlays = [&](double t0, double t1, double q0, double q1, const char* dash, double stroke_w) {
        if (gene_bands.empty() || t1 <= t0) {
            return;
        }
        const double ref_span = t1 - t0;
        for (const auto& band : gene_bands) {
            const double band_start = static_cast<double>(band.ref_start_bp);
            const double band_end = static_cast<double>(band.ref_end_bp);
            if (band_end <= t0 || band_start >= t1) {
                continue;
            }
            const double overlap_start = std::max(t0, band_start);
            const double overlap_end = std::min(t1, band_end);
            if (overlap_end <= overlap_start) {
                continue;
            }
            const auto color_it = gene_color_by_name.find(band.name);
            if (color_it == gene_color_by_name.end()) {
                continue;
            }
            const double frac0 = (overlap_start - t0) / ref_span;
            const double frac1 = (overlap_end - t0) / ref_span;
            const double q_overlap_start = q0 + frac0 * (q1 - q0);
            const double q_overlap_end = q0 + frac1 * (q1 - q0);
            out << "<line x1=\"" << plot_xd(overlap_start)
                << "\" y1=\"" << plot_yd(q_overlap_start)
                << "\" x2=\"" << plot_xd(overlap_end)
                << "\" y2=\"" << plot_yd(q_overlap_end)
                << "\" stroke=\"" << color_it->second
                << "\" stroke-width=\"" << (stroke_w + 1.00)
                << "\" stroke-opacity=\"0.95\" stroke-linecap=\"round\""
                << dash << "/>\n";
        }
    };

    for (const auto& rec : records) {
        const char* color = "#1F77B4";
        const char* dash = rec.primary ? "" : " stroke-dasharray=\"3,2\"";
        const double stroke_w = rec.primary ? 1.35 : 1.05;
        const double stroke_alpha = rec.primary ? 0.90 : 0.75;

        bool emitted = false;
        if (!rec.cigar.empty()) {
            const auto runs = parse_extended_cigar(rec.cigar);
            std::size_t q_oriented = 0;
            std::size_t t_oriented = 0;
            const std::size_t q_span = rec.query_end >= rec.query_start ? (rec.query_end - rec.query_start) : 0;
            const std::size_t t_span = rec.target_end >= rec.target_start ? (rec.target_end - rec.target_start) : 0;
            for (const auto& run : runs) {
                if (run.len == 0) {
                    continue;
                }
                if (run.op == '=' || run.op == 'M' || run.op == 'X') {
                    if (q_oriented + run.len > q_span || t_oriented + run.len > t_span) {
                        break;
                    }
                    const std::size_t q0 = query_forward_coord_from_record(rec, q_oriented);
                    const std::size_t q1 = query_forward_coord_from_record(rec, q_oriented + run.len);
                    const std::size_t t0 = rec.target_start + t_oriented;
                    const std::size_t t1 = rec.target_start + t_oriented + run.len;
                    const double x1 = plot_x(t0);
                    const double y1 = plot_y(q0);
                    const double x2 = plot_x(t1);
                    const double y2 = plot_y(q1);
                    out << "<line x1=\"" << x1 << "\" y1=\"" << y1
                        << "\" x2=\"" << x2 << "\" y2=\"" << y2
                        << "\" stroke=\"" << color << "\" stroke-width=\"" << stroke_w
                        << "\" stroke-opacity=\"" << stroke_alpha << "\" stroke-linecap=\"round\""
                        << dash << "/>\n";
                    draw_gene_overlays(
                        static_cast<double>(t0),
                        static_cast<double>(t1),
                        static_cast<double>(q0),
                        static_cast<double>(q1),
                        dash,
                        stroke_w);
                    q_oriented += run.len;
                    t_oriented += run.len;
                    emitted = true;
                    continue;
                }
                if (run.op == 'I') {
                    if (q_oriented + run.len > q_span) {
                        break;
                    }
                    q_oriented += run.len;
                    continue;
                }
                if (run.op == 'D') {
                    if (t_oriented + run.len > t_span) {
                        break;
                    }
                    t_oriented += run.len;
                    continue;
                }
            }
        }

        if (!emitted) {
            const double x1 = plot_x(rec.target_start);
            const double x2 = plot_x(rec.target_end);
            const double y1 = plot_y(rec.strand == '-' ? rec.query_end : rec.query_start);
            const double y2 = plot_y(rec.strand == '-' ? rec.query_start : rec.query_end);
            out << "<line x1=\"" << x1 << "\" y1=\"" << y1
                << "\" x2=\"" << x2 << "\" y2=\"" << y2
                << "\" stroke=\"" << color << "\" stroke-width=\"" << stroke_w
                << "\" stroke-opacity=\"" << stroke_alpha << "\"" << dash << "/>\n";
            draw_gene_overlays(
                static_cast<double>(rec.target_start),
                static_cast<double>(rec.target_end),
                static_cast<double>(rec.strand == '-' ? rec.query_end : rec.query_start),
                static_cast<double>(rec.strand == '-' ? rec.query_start : rec.query_end),
                dash,
                stroke_w);
        }
    }

    // Variant breakpoint guides:
    // DEL -> vertical boundaries on reference axis
    // INS -> single vertical insertion breakpoint on reference axis
    // INV -> both reference and assembly boundaries
    const char* bp_color = "#8A8A8A";
    const char* bp_dash = "5,4";
    const double bp_width = 1.15;
    for (const auto& ev : called_events) {
        const std::size_t ref_start = clamp_ref_bp(std::min(ev.ref_offset_start_bp, ev.ref_offset_end_bp));
        const std::size_t ref_end = clamp_ref_bp(std::max(ev.ref_offset_start_bp, ev.ref_offset_end_bp));
        const std::size_t qry_start = clamp_qry_bp(std::min(ev.alt_offset_start_bp, ev.alt_offset_end_bp));
        const std::size_t qry_end = clamp_qry_bp(std::max(ev.alt_offset_start_bp, ev.alt_offset_end_bp));

        auto draw_vertical = [&](std::size_t ref_pos) {
            const double x = plot_x(ref_pos);
            out << "<line x1=\"" << x << "\" y1=\"" << margin
                << "\" x2=\"" << x << "\" y2=\"" << (margin + plot_h)
                << "\" stroke=\"" << bp_color << "\" stroke-width=\"" << bp_width
                << "\" stroke-opacity=\"0.90\" stroke-dasharray=\"" << bp_dash << "\"/>\n";
        };
        auto draw_horizontal = [&](std::size_t qry_pos) {
            const double y = plot_y(qry_pos);
            out << "<line x1=\"" << margin << "\" y1=\"" << y
                << "\" x2=\"" << (margin + plot_w) << "\" y2=\"" << y
                << "\" stroke=\"" << bp_color << "\" stroke-width=\"" << bp_width
                << "\" stroke-opacity=\"0.90\" stroke-dasharray=\"" << bp_dash << "\"/>\n";
        };

        if (ev.event_type == "DEL") {
            draw_vertical(ref_start);
            if (ref_end != ref_start) {
                draw_vertical(ref_end);
            }
            continue;
        }
        if (ev.event_type == "INS") {
            const std::size_t ref_anchor = ref_start;
            draw_vertical(ref_anchor);
            continue;
        }
        if (ev.event_type == "INV") {
            draw_vertical(ref_start);
            if (ref_end != ref_start) {
                draw_vertical(ref_end);
            }
            draw_horizontal(qry_start);
            if (qry_end != qry_start) {
                draw_horizontal(qry_end);
            }
        }
    }

    // Reference coordinate ticks (equally spaced and intentionally sparse).
    const ParsedReferencePath ref_meta = parse_reference_path_label(ref_for_gtf);
    const bool has_ref_interval = ref_meta.has_interval;
    constexpr std::size_t tick_intervals = 5;
    const double axis_y = margin + plot_h;
    for (std::size_t i = 0; i <= tick_intervals; ++i) {
        const std::size_t local_bp = (ref_len * i) / tick_intervals;
        const double x = plot_x(local_bp);
        out << "<line x1=\"" << x << "\" y1=\"" << axis_y
            << "\" x2=\"" << x << "\" y2=\"" << (axis_y + 6.0)
            << "\" stroke=\"#B5B5B5\" stroke-width=\"1\"/>\n";
        std::size_t label_bp = local_bp;
        if (has_ref_interval) {
            label_bp = ref_meta.region_start_1based + local_bp;
        }
        out << "<text x=\"" << x << "\" y=\"" << (axis_y + 18.0)
            << "\" text-anchor=\"middle\" font-family=\"Arial, sans-serif\" font-size=\"10\" fill=\"#666666\">"
            << format_bp(label_bp) << "</text>\n";
    }

    out << "<text x=\"" << (margin + (plot_w / 2.0)) << "\" y=\"" << (margin + plot_h + 30.0)
        << "\" text-anchor=\"middle\" font-family=\"Arial, sans-serif\" font-size=\"12\">"
        << html_escape(ref_name) << "</text>\n";
    out << "<text x=\"22\" y=\"" << (margin + (plot_h / 2.0))
        << "\" transform=\"rotate(-90,22," << (margin + (plot_h / 2.0))
        << ")\" text-anchor=\"middle\" font-family=\"Arial, sans-serif\" font-size=\"12\">"
        << html_escape(query_name) << "</text>\n";

    if (!gene_legend_order.empty()) {
        const double legend_title_y = margin + plot_h + 56.0;
        out << "<text x=\"" << margin << "\" y=\"" << legend_title_y
            << "\" font-family=\"Arial, sans-serif\" font-size=\"12\" font-weight=\"bold\">"
            << "Gene highlights (" << html_escape(gene_match_label) << ")</text>\n";

        double cursor_x = margin;
        double cursor_y = legend_title_y + 22.0;
        const double x_limit = width - margin;
        const double y_limit = height - 12.0;
        bool truncated = false;
        for (const auto& gene_name : gene_legend_order) {
            const auto color_it = gene_color_by_name.find(gene_name);
            if (color_it == gene_color_by_name.end()) {
                continue;
            }
            const double entry_w = 36.0 + std::min<double>(300.0, static_cast<double>(gene_name.size()) * 6.8);
            if (cursor_x + entry_w > x_limit) {
                cursor_x = margin;
                cursor_y += 20.0;
            }
            if (cursor_y > y_limit) {
                truncated = true;
                break;
            }
            out << "<line x1=\"" << cursor_x << "\" y1=\"" << cursor_y
                << "\" x2=\"" << (cursor_x + 18.0) << "\" y2=\"" << cursor_y
                << "\" stroke=\"" << color_it->second
                << "\" stroke-width=\"3\" stroke-linecap=\"round\"/>\n";
            out << "<text x=\"" << (cursor_x + 24.0) << "\" y=\"" << (cursor_y + 4.0)
                << "\" font-family=\"Arial, sans-serif\" font-size=\"11\">"
                << html_escape(gene_name) << "</text>\n";
            cursor_x += entry_w;
        }
        if (truncated) {
            out << "<text x=\"" << margin << "\" y=\"" << std::min(height - 8.0, cursor_y + 14.0)
                << "\" font-family=\"Arial, sans-serif\" font-size=\"10\" fill=\"#666666\">"
                << "... more genes omitted from legend"
                << "</text>\n";
        }
    }
    out << "</svg>\n";
}

[[maybe_unused]] void normalize_cluster_cn_deltas(std::vector<AtomicVariantEvent>& events, std::size_t min_sv_bp) {
    if (events.empty()) {
        return;
    }

    struct CnEventAnchor {
        std::size_t idx = 0;
        std::size_t start = 0;
        std::size_t end = 0;
        int delta = 0;
        std::size_t event_bp = 0;
    };

    auto event_span = [](const AtomicVariantEvent& ev) -> std::pair<std::size_t, std::size_t> {
        std::size_t start = std::min(ev.ref_offset_start_bp, ev.alt_offset_start_bp);
        std::size_t end = std::max(ev.ref_offset_end_bp, ev.alt_offset_end_bp);
        if (end <= start) {
            const std::size_t fallback = std::max<std::size_t>(1, ev.event_bp);
            end = start + fallback;
        }
        return {start, end};
    };

    std::vector<CnEventAnchor> anchors;
    anchors.reserve(events.size());
    for (std::size_t i = 0; i < events.size(); ++i) {
        const auto& ev = events[i];
        if (ev.cn_delta == 0) {
            continue;
        }
        const auto [start, end] = event_span(ev);
        anchors.push_back({i, start, end, ev.cn_delta, ev.event_bp});
    }
    if (anchors.empty()) {
        return;
    }

    for (const auto& anchor : anchors) {
        events[anchor.idx].cn_delta = 0;
    }

    // Preserve multiple independent CN events in one cluster while collapsing
    // close-by fragments from alignment decomposition into a single CN carrier.
    const std::size_t merge_gap_bp = std::max<std::size_t>(10, min_sv_bp);
    const std::size_t no_index = std::numeric_limits<std::size_t>::max();
    for (int sign : {-1, 1}) {
        std::vector<CnEventAnchor> local;
        local.reserve(anchors.size());
        for (const auto& anchor : anchors) {
            if ((sign < 0 && anchor.delta < 0) || (sign > 0 && anchor.delta > 0)) {
                local.push_back(anchor);
            }
        }
        if (local.empty()) {
            continue;
        }
        std::sort(local.begin(), local.end(), [](const CnEventAnchor& lhs, const CnEventAnchor& rhs) {
            if (lhs.start != rhs.start) {
                return lhs.start < rhs.start;
            }
            if (lhs.end != rhs.end) {
                return lhs.end < rhs.end;
            }
            return lhs.idx < rhs.idx;
        });

        std::size_t group_end = local.front().end;
        std::size_t representative_idx = local.front().idx;
        std::size_t representative_bp = local.front().event_bp;
        int representative_delta = local.front().delta;

        auto flush_group = [&]() {
            if (representative_idx == no_index) {
                return;
            }
            if (sign < 0) {
                events[representative_idx].cn_delta = -std::max(1, std::abs(representative_delta));
            } else {
                events[representative_idx].cn_delta = std::max(1, representative_delta);
            }
        };

        auto consider_representative = [&](const CnEventAnchor& anchor) {
            if (sign < 0) {
                const int current_abs = std::abs(representative_delta);
                const int candidate_abs = std::abs(anchor.delta);
                if (candidate_abs > current_abs ||
                    (candidate_abs == current_abs && anchor.event_bp > representative_bp) ||
                    (candidate_abs == current_abs && anchor.event_bp == representative_bp && anchor.idx < representative_idx)) {
                    representative_idx = anchor.idx;
                    representative_bp = anchor.event_bp;
                    representative_delta = anchor.delta;
                }
                return;
            }
            if (anchor.delta > representative_delta ||
                (anchor.delta == representative_delta && anchor.event_bp > representative_bp) ||
                (anchor.delta == representative_delta && anchor.event_bp == representative_bp && anchor.idx < representative_idx)) {
                representative_idx = anchor.idx;
                representative_bp = anchor.event_bp;
                representative_delta = anchor.delta;
            }
        };

        for (std::size_t i = 1; i < local.size(); ++i) {
            const CnEventAnchor& anchor = local[i];
            if (anchor.start <= group_end + merge_gap_bp) {
                group_end = std::max(group_end, anchor.end);
                consider_representative(anchor);
                continue;
            }

            flush_group();
            group_end = anchor.end;
            representative_idx = anchor.idx;
            representative_bp = anchor.event_bp;
            representative_delta = anchor.delta;
        }
        flush_group();
    }
}

std::vector<AtomicVariantEvent> merge_nearby_events(
    const std::vector<AtomicVariantEvent>& events,
    const std::string& ref_seq,
    const std::string& alt_seq,
    std::size_t merge_gap_bp,
    std::size_t min_sv_bp) {

    auto event_ref_span_bp = [](const AtomicVariantEvent& ev) -> std::size_t {
        if (ev.ref_offset_end_bp <= ev.ref_offset_start_bp) {
            return 0;
        }
        return ev.ref_offset_end_bp - ev.ref_offset_start_bp;
    };
    auto event_alt_span_bp = [](const AtomicVariantEvent& ev) -> std::size_t {
        if (ev.alt_offset_end_bp <= ev.alt_offset_start_bp) {
            return 0;
        }
        return ev.alt_offset_end_bp - ev.alt_offset_start_bp;
    };
    auto effective_span_bp = [&](const AtomicVariantEvent& ev) -> std::size_t {
        if (ev.event_type == "DEL") {
            const std::size_t ref_span = event_ref_span_bp(ev);
            return ref_span > 0 ? ref_span : ev.event_bp;
        }
        if (ev.event_type == "INS") {
            const std::size_t alt_span = event_alt_span_bp(ev);
            return alt_span > 0 ? alt_span : ev.event_bp;
        }
        return ev.event_bp;
    };
    auto cross_source_size_similarity = [&](const AtomicVariantEvent& lhs, const AtomicVariantEvent& rhs) -> double {
        const std::size_t lhs_bp = std::max<std::size_t>(1, effective_span_bp(lhs));
        const std::size_t rhs_bp = std::max<std::size_t>(1, effective_span_bp(rhs));
        const std::size_t small = std::min(lhs_bp, rhs_bp);
        const std::size_t large = std::max(lhs_bp, rhs_bp);
        return static_cast<double>(small) / static_cast<double>(large);
    };
    auto can_cross_source_merge = [&](const AtomicVariantEvent& lhs, const AtomicVariantEvent& rhs) -> bool {
        if (lhs.evidence_source == rhs.evidence_source) {
            return true;
        }
        if (lhs.evidence_source == "unknown" || rhs.evidence_source == "unknown") {
            return true;
        }
        // Prevent over-collapsing distinct signals (e.g. a long split-chain DEL
        // and a shorter CIGAR DEL).
        constexpr double kMinCrossSourceSizeSimilarity = 0.70;
        return cross_source_size_similarity(lhs, rhs) >= kMinCrossSourceSizeSimilarity;
    };

    std::vector<AtomicVariantEvent> sorted = events;
    std::sort(sorted.begin(), sorted.end(), [](const AtomicVariantEvent& lhs, const AtomicVariantEvent& rhs) {
        if (lhs.event_type != rhs.event_type) {
            return lhs.event_type < rhs.event_type;
        }
        if (lhs.ref_offset_start_bp != rhs.ref_offset_start_bp) {
            return lhs.ref_offset_start_bp < rhs.ref_offset_start_bp;
        }
        if (lhs.alt_offset_start_bp != rhs.alt_offset_start_bp) {
            return lhs.alt_offset_start_bp < rhs.alt_offset_start_bp;
        }
        return lhs.event_bp < rhs.event_bp;
    });

    std::vector<AtomicVariantEvent> merged;
    merged.reserve(sorted.size());
    for (const auto& ev : sorted) {
        if (!(ev.event_type == "INS" || ev.event_type == "DEL")) {
            if (ev.event_bp >= min_sv_bp) {
                merged.push_back(ev);
            }
            continue;
        }
        if (merged.empty() || merged.back().event_type != ev.event_type) {
            merged.push_back(ev);
            continue;
        }

        AtomicVariantEvent& curr = merged.back();
        const bool ref_near = ev.ref_offset_start_bp <= (curr.ref_offset_end_bp + merge_gap_bp);
        const bool alt_near = ev.alt_offset_start_bp <= (curr.alt_offset_end_bp + merge_gap_bp);
        if (!ref_near || !alt_near) {
            merged.push_back(ev);
            continue;
        }
        if (!can_cross_source_merge(curr, ev)) {
            merged.push_back(ev);
            continue;
        }

        curr.ref_offset_end_bp = std::max(curr.ref_offset_end_bp, ev.ref_offset_end_bp);
        curr.alt_offset_end_bp = std::max(curr.alt_offset_end_bp, ev.alt_offset_end_bp);
        curr.ref_offset_start_bp = std::min(curr.ref_offset_start_bp, ev.ref_offset_start_bp);
        curr.alt_offset_start_bp = std::min(curr.alt_offset_start_bp, ev.alt_offset_start_bp);
        curr.preserve_ins_svlen = curr.preserve_ins_svlen || ev.preserve_ins_svlen;

        if (curr.event_type == "DEL") {
            if (curr.ref_offset_end_bp > curr.ref_offset_start_bp &&
                curr.ref_offset_end_bp <= ref_seq.size()) {
                curr.event_bp = curr.ref_offset_end_bp - curr.ref_offset_start_bp;
                curr.svlen = -static_cast<long long>(curr.event_bp);
                curr.inserted_bp = 0;
                curr.inserted_seq.clear();
            }
        } else if (curr.event_type == "INS") {
            std::size_t alt_span_bp = 0;
            if (curr.alt_offset_end_bp > curr.alt_offset_start_bp &&
                curr.alt_offset_end_bp <= alt_seq.size()) {
                alt_span_bp = curr.alt_offset_end_bp - curr.alt_offset_start_bp;
            }
            std::size_t merged_bp = alt_span_bp;
            if (curr.preserve_ins_svlen) {
                const std::size_t curr_svlen_bp = curr.svlen > 0 ? static_cast<std::size_t>(curr.svlen) : 0;
                const std::size_t ev_svlen_bp = ev.svlen > 0 ? static_cast<std::size_t>(ev.svlen) : 0;
                merged_bp = std::max(merged_bp, curr.event_bp);
                merged_bp = std::max(merged_bp, ev.event_bp);
                merged_bp = std::max(merged_bp, curr_svlen_bp);
                merged_bp = std::max(merged_bp, ev_svlen_bp);
            }
            curr.event_bp = merged_bp;
            curr.svlen = static_cast<long long>(merged_bp);
            curr.inserted_bp = merged_bp;
            if (alt_span_bp > 0 && curr.alt_offset_start_bp + alt_span_bp <= alt_seq.size()) {
                curr.inserted_seq = alt_seq.substr(curr.alt_offset_start_bp, alt_span_bp);
            } else {
                curr.inserted_seq.clear();
            }
        }
    }

    std::vector<AtomicVariantEvent> filtered;
    filtered.reserve(merged.size());
    for (auto& ev : merged) {
        if (ev.event_bp >= min_sv_bp || ev.event_type == "INV") {
            filtered.push_back(std::move(ev));
        }
    }
    return filtered;
}

std::string event_sequence_for_vcf(
    const AtomicVariantEvent& ev,
    const std::string& ref_seq,
    const std::string& alt_seq) {

    if (ev.event_type == "INS") {
        if (!ev.inserted_seq.empty()) {
            return ev.inserted_seq;
        }
        if (ev.alt_offset_end_bp > ev.alt_offset_start_bp &&
            ev.alt_offset_end_bp <= alt_seq.size()) {
            return alt_seq.substr(ev.alt_offset_start_bp, ev.alt_offset_end_bp - ev.alt_offset_start_bp);
        }
        return {};
    }
    if (ev.event_type == "DEL" || ev.event_type == "INV") {
        if (ev.ref_offset_end_bp > ev.ref_offset_start_bp &&
            ev.ref_offset_end_bp <= ref_seq.size()) {
            return ref_seq.substr(ev.ref_offset_start_bp, ev.ref_offset_end_bp - ev.ref_offset_start_bp);
        }
        return {};
    }
    return {};
}

std::string atomic_event_label(const AtomicVariantEvent& ev) {
    std::ostringstream oss;
    oss << ev.event_type;
    if (ev.event_subtype != ".") {
        oss << ":" << ev.event_subtype;
    }
    oss << " ref[" << ev.ref_offset_start_bp << "," << ev.ref_offset_end_bp << ")"
        << " alt[" << ev.alt_offset_start_bp << "," << ev.alt_offset_end_bp << ")"
        << " bp=" << ev.event_bp
        << " svlen=" << ev.svlen
        << " cn_delta=" << ev.cn_delta
        << " src=" << ev.evidence_source;
    if (ev.inserted_bp > 0) {
        oss << " inserted_bp=" << ev.inserted_bp;
    }
    if (ev.has_dup_evidence) {
        oss << " dup_sim=" << std::fixed << std::setprecision(4) << ev.dup_best_similarity
            << " dup_unit_bp=" << ev.dup_unit_bp
            << " ref_cn=" << ev.dup_ref_copy_number
            << " alt_cn=" << ev.dup_alt_copy_number
            << " added=" << ev.dup_added_copies;
    }
    return oss.str();
}

std::pair<std::string, std::string> vcf_alt_and_type(
    const std::string& event_type,
    const std::string& event_subtype);

void write_cluster_pairwise_vcf(
    const std::string& output_path,
    const ParsedReferencePath& reference_meta,
    const std::string& reference_label,
    std::size_t bubble_id,
    std::size_t cluster_id,
    std::size_t reference_cluster_id,
    const std::string& orientation,
    double best_norm,
    const std::string& ref_seq,
    const std::string& alt_seq,
    const std::vector<AtomicVariantEvent>& events) {

    std::ofstream out(output_path);
    if (!out) {
        throw std::runtime_error("Failed to write cluster pairwise VCF: " + output_path);
    }

    out << "##fileformat=VCFv4.2\n";
    out << "##source=panvar\n";
    out << "##reference=" << reference_label << "\n";
    out << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant\">\n";
    out << "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Structural variant type\">\n";
    out << "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length difference ALT-REF\">\n";
    out << "##INFO=<ID=BUBBLE_ID,Number=1,Type=Integer,Description=\"panvar bubble identifier\">\n";
    out << "##INFO=<ID=CLUSTER_ID,Number=1,Type=Integer,Description=\"Cluster identifier within bubble\">\n";
    out << "##INFO=<ID=REF_CLUSTER_ID,Number=1,Type=Integer,Description=\"Reference cluster identifier within bubble\">\n";
    out << "##INFO=<ID=EVENT,Number=1,Type=String,Description=\"panvar event type\">\n";
    out << "##INFO=<ID=INS_SUBTYPE,Number=1,Type=String,Description=\"INS subtype (NOVEL or DUP-like subtype)\">\n";
    out << "##INFO=<ID=ORIENT,Number=1,Type=String,Description=\"Best orientation of cluster vs reference\">\n";
    out << "##INFO=<ID=BEST_NORM_ED,Number=1,Type=Float,Description=\"Best normalized edit distance to reference cluster\">\n";
    out << "##INFO=<ID=DUP_SIM,Number=1,Type=Float,Description=\"Best similarity of inserted sequence to duplicated source candidate\">\n";
    out << "##INFO=<ID=DUP_REF_START,Number=1,Type=Integer,Description=\"1-based start of duplicated source interval on reference\">\n";
    out << "##INFO=<ID=DUP_REF_END,Number=1,Type=Integer,Description=\"1-based end of duplicated source interval on reference\">\n";
    out << "##INFO=<ID=DUP_ORIENT,Number=1,Type=String,Description=\"Orientation of duplicated source match (+/-)\">\n";
    out << "##INFO=<ID=DUP_UNIT_BP,Number=1,Type=Integer,Description=\"Duplicated source unit length in bp\">\n";
    out << "##INFO=<ID=DUP_REF_CN,Number=1,Type=Integer,Description=\"Estimated reference copy number of duplication unit\">\n";
    out << "##INFO=<ID=DUP_ALT_CN,Number=1,Type=Integer,Description=\"Estimated ALT copy number of duplication unit\">\n";
    out << "##INFO=<ID=DUP_ADDED,Number=1,Type=Integer,Description=\"Estimated number of added copies in ALT\">\n";
    out << "##INFO=<ID=DUP_COPY_RATIO,Number=1,Type=Float,Description=\"Estimated ALT/REF copy ratio for duplication unit\">\n";
    out << "##INFO=<ID=INSSEQ,Number=1,Type=String,Description=\"Inserted sequence\">\n";
    out << "##INFO=<ID=DELSEQ,Number=1,Type=String,Description=\"Deleted sequence from reference\">\n";
    out << "##INFO=<ID=INVSEQ,Number=1,Type=String,Description=\"Inverted reference sequence\">\n";
    out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

    const std::size_t region_start_1based = reference_meta.has_interval ? reference_meta.region_start_1based : 1;
    const std::string chrom = reference_meta.chrom.empty() ? std::string("ref") : reference_meta.chrom;
    for (std::size_t i = 0; i < events.size(); ++i) {
        const auto& ev = events[i];
        const auto [alt, svtype] = vcf_alt_and_type(ev.event_type, ev.event_subtype);
        const std::size_t start_off = std::min(ev.ref_offset_start_bp, ev.ref_offset_end_bp);
        const std::size_t end_off = std::max(ev.ref_offset_start_bp, ev.ref_offset_end_bp);

        std::size_t pos_1based = region_start_1based + start_off;
        if (pos_1based == 0) {
            pos_1based = 1;
        }
        std::size_t info_end = region_start_1based + end_off;
        if (svtype == "INS" || info_end < pos_1based) {
            info_end = pos_1based;
        }

        const std::size_t ref_base_off = std::min(start_off, ref_seq.empty() ? 0 : (ref_seq.size() - 1));
        const std::string ref_base =
            ref_seq.empty() ? std::string("N") : std::string(1, vcf_base(ref_seq[ref_base_off]));

        const std::string event_seq = event_sequence_for_vcf(ev, ref_seq, alt_seq);
        long long svlen = ev.svlen;
        if (!event_seq.empty()) {
            if (svtype == "DEL") {
                svlen = -static_cast<long long>(event_seq.size());
            } else if (svtype == "INS") {
                if (!ev.preserve_ins_svlen) {
                    svlen = static_cast<long long>(event_seq.size());
                } else if (svlen <= 0) {
                    svlen = static_cast<long long>(event_seq.size());
                }
            } else if (svtype == "INV") {
                svlen = static_cast<long long>(event_seq.size());
            }
        }

        out << chrom << '\t'
            << pos_1based << '\t'
            << "B" << bubble_id << "_C" << cluster_id << "_E" << (i + 1) << '\t'
            << ref_base << '\t'
            << alt << '\t'
            << ".\tPASS\t";
        out << "END=" << info_end
            << ";SVTYPE=" << svtype
            << ";SVLEN=" << svlen
            << ";BUBBLE_ID=" << bubble_id
            << ";CLUSTER_ID=" << cluster_id
            << ";REF_CLUSTER_ID=" << reference_cluster_id
            << ";EVENT=" << ev.event_type
            << ";INS_SUBTYPE=" << ev.event_subtype
            << ";ORIENT=" << orientation
            << ";BEST_NORM_ED=" << std::fixed << std::setprecision(6) << best_norm;
        if (ev.has_dup_evidence) {
            const std::size_t dup_start_1based = region_start_1based + ev.dup_ref_start_bp;
            const std::size_t dup_end_1based = region_start_1based + ev.dup_ref_end_bp;
            out << ";DUP_SIM=" << std::fixed << std::setprecision(6) << ev.dup_best_similarity
                << ";DUP_REF_START=" << dup_start_1based
                << ";DUP_REF_END=" << dup_end_1based
                << ";DUP_ORIENT=" << ev.dup_orientation
                << ";DUP_UNIT_BP=" << ev.dup_unit_bp
                << ";DUP_REF_CN=" << ev.dup_ref_copy_number
                << ";DUP_ALT_CN=" << ev.dup_alt_copy_number
                << ";DUP_ADDED=" << ev.dup_added_copies
                << ";DUP_COPY_RATIO=" << std::fixed << std::setprecision(6) << ev.dup_copy_ratio;
        }
        if (!event_seq.empty()) {
            if (svtype == "INS") {
                out << ";INSSEQ=" << event_seq;
            } else if (svtype == "DEL") {
                out << ";DELSEQ=" << event_seq;
            } else if (svtype == "INV") {
                out << ";INVSEQ=" << event_seq;
            }
        }
        out << '\n';
    }
}

[[maybe_unused]] void write_cluster_debug_trace(
    const std::string& output_path,
    const Bubble& bubble,
    const std::string& reference_path,
    std::size_t reference_cluster_id,
    std::size_t cluster_id,
    std::size_t representative_allele_id,
    std::size_t min_sv_bp,
    double dup_min_similarity,
    std::size_t ref_len,
    std::size_t alt_len,
    long long len_delta,
    int orientation_margin,
    const std::string& orientation,
    int dist_fwd,
    int dist_rev,
    int best_dist,
    double best_norm,
    const std::string& best_cigar,
    const std::string& cigar_tsv_path,
    const std::vector<CigarTraceRow>& cigar_rows,
    const std::vector<AtomicVariantEvent>& events_before_cn_norm,
    const std::vector<AtomicVariantEvent>& events_after_cn_norm) {

    std::ofstream out(output_path);
    if (!out) {
        throw std::runtime_error("Failed to write call debug trace: " + output_path);
    }
    (void)dup_min_similarity;

    std::size_t sv_sized_indels = 0;
    std::size_t ins_ops = 0;
    std::size_t del_ops = 0;
    for (const auto& row : cigar_rows) {
        if (row.op == 'D') {
            ++ins_ops;
        } else if (row.op == 'I') {
            ++del_ops;
        }
        if (row.sv_sized) {
            ++sv_sized_indels;
        }
    }

    out << "bubble_id: " << bubble.id << "\n";
    out << "source: " << bubble.source << "\n";
    out << "sink: " << bubble.sink << "\n";
    out << "reference_path: " << reference_path << "\n";
    out << "reference_cluster_id: " << reference_cluster_id << "\n";
    out << "cluster_id: " << cluster_id << "\n";
    out << "representative_allele_id: " << representative_allele_id << "\n\n";

    out << "STEP 1 - WE_HAVE\n";
    out << "ref_len_bp=" << ref_len
        << ", alt_len_bp=" << alt_len
        << ", len_delta_bp=" << len_delta
        << ", min_sv_bp=" << min_sv_bp
        << ", ins_typing=disabled"
        << "\n\n";

    out << "STEP 2 - WE_LOOK_FOR\n";
    out << "1) Best orientation from forward vs reverse-complement alignment\n";
    out << "2) SV-sized INS/DEL operations in CIGAR (>= min_sv_bp)\n";
    out << "3) Split-read junction deviations across PAF chains (>= min_sv_bp)\n";
    out << "4) If no SV-sized indel but large net length shift: fallback net INS/DEL by prefix/suffix\n";
    out << "5) Merge nearby concordant events (small match gaps collapsed)\n";
    out << "6) Emit only INV/DEL/INS in this simplified mode (no INS subtype typing)\n\n";

    out << "STEP 3 - WE_FOUND\n";
    out << "alignment_forward_edit_distance=" << dist_fwd << "\n";
    out << "alignment_reverse_edit_distance=" << dist_rev << "\n";
    out << "orientation_margin_bp=" << orientation_margin << "\n";
    out << "chosen_orientation=" << orientation << "\n";
    out << "best_edit_distance=" << best_dist << "\n";
    out << "best_edit_distance_norm=" << std::fixed << std::setprecision(6) << best_norm << "\n";
    out << "best_cigar_preview=" << preview_text(best_cigar, 160) << "\n";
    out << "cigar_trace_tsv=" << cigar_tsv_path << "\n";
    out << "cigar_ops=" << cigar_rows.size()
        << ", sv_sized_indel_ops=" << sv_sized_indels
        << ", ins_ops(raw D)=" << ins_ops
        << ", del_ops(raw I)=" << del_ops
        << "\n";
    if (events_before_cn_norm.empty()) {
        out << "atomic_events_before_cn_norm=none\n";
    } else {
        out << "atomic_events_before_cn_norm:\n";
        for (std::size_t i = 0; i < events_before_cn_norm.size(); ++i) {
            out << "  event_" << (i + 1) << ": " << atomic_event_label(events_before_cn_norm[i]) << "\n";
        }
    }
    out << "\n";

    out << "STEP 4 - DECISION\n";
    if (events_after_cn_norm.empty()) {
        out << "final_call=MATCH (no qualifying SV event from this cluster)\n";
    } else {
        out << "final_atomic_events:\n";
        for (std::size_t i = 0; i < events_after_cn_norm.size(); ++i) {
            out << "  event_" << (i + 1) << ": " << atomic_event_label(events_after_cn_norm[i]);
            out << "\n";
        }
    }
}

VariantBubbleReport write_variant_reports_for_bubble(
    const Bubble& bubble,
    const std::vector<UniqueAllele>& unique_alleles,
    const std::vector<Cluster>& clusters,
    const std::vector<std::size_t>& cluster_of_unique,
    const std::vector<PathAssignment>& assignments,
    const std::string& reference_path,
    const std::string& reference_sequence_global,
    std::size_t min_sv_bp,
    bool write_debug_files,
    const std::string& debug_output_dir,
    const GeneAnnotationIndex* dotplot_gene_index,
    const std::vector<std::string>& dotplot_gene_matches,
    const std::string& minimap_preset,
    std::size_t minimap_best_n,
    bool minimap_emit_secondary,
    bool split_ins_use_geometric_svlen,
    bool classify_ins) {

    VariantBubbleReport report;
    report.clusters_total = clusters.size();

    std::string bubble_debug_dir;
    if (write_debug_files) {
        if (debug_output_dir.empty()) {
            throw std::runtime_error("--debug-out-dir cannot be empty when call debug output is enabled");
        }
        bubble_debug_dir = debug_output_dir + "/bubble_" + std::to_string(bubble.id);
        std::filesystem::create_directories(bubble_debug_dir);
    }

    auto write_cluster_status_text = [&](const std::string& cluster_debug_dir,
                                         const VariantBubbleReport::ClusterDebugStatus& status) {
        if (!write_debug_files || cluster_debug_dir.empty()) {
            return;
        }
        const std::string status_path = cluster_debug_dir + "/status.txt";
        std::ofstream out(status_path);
        if (!out) {
            throw std::runtime_error("Failed to write cluster status text: " + status_path);
        }
        out << "bubble_id=" << bubble.id << "\n";
        out << "cluster_id=" << status.cluster_id << "\n";
        out << "status=" << status.status << "\n";
        out << "is_reference_cluster=" << (status.is_reference_cluster ? 1 : 0) << "\n";
        out << "representative_allele_id=" << status.representative_allele_id << "\n";
        out << "representative_haplotype=" << status.representative_haplotype << "\n";
        out << "reference_length_bp=" << status.reference_length_bp << "\n";
        out << "cluster_length_bp=" << status.cluster_length_bp << "\n";
        out << "minimap_best_ok=" << (status.minimap_best_ok ? 1 : 0) << "\n";
        out << "orientation=" << status.orientation << "\n";
        out << std::fixed << std::setprecision(6);
        out << "best_norm_ed=" << status.best_edit_distance_norm << "\n";
        out << "paf_records=" << status.paf_records << "\n";
        out << "event_count=" << status.event_count << "\n";
        out << "dotplot_written=" << (status.dotplot_written ? 1 : 0) << "\n";
        out << "pairwise_vcf_written=" << (status.pairwise_vcf_written ? 1 : 0) << "\n";
    };

    auto write_bubble_debug_status = [&]() {
        if (!write_debug_files || bubble_debug_dir.empty()) {
            return;
        }
        const std::string bubble_status_path = bubble_debug_dir + "/bubble_status.tsv";
        std::ofstream bubble_out(bubble_status_path);
        if (!bubble_out) {
            throw std::runtime_error("Failed to write bubble debug status TSV: " + bubble_status_path);
        }
        std::size_t ok_clusters = 0;
        std::size_t clusters_with_events = 0;
        std::size_t dotplots = 0;
        std::size_t cluster_vcfs = 0;
        for (const auto& cluster_status : report.cluster_statuses) {
            if (cluster_status.minimap_best_ok) {
                report.clusters_with_minimap_hit += 1;
            }
            if (cluster_status.event_count > 0) {
                clusters_with_events += 1;
            }
            if (cluster_status.dotplot_written) {
                dotplots += 1;
            }
            if (cluster_status.pairwise_vcf_written) {
                cluster_vcfs += 1;
            }
            if (cluster_status.status.rfind("ok", 0) == 0) {
                ok_clusters += 1;
            }
        }
        report.clusters_with_events = clusters_with_events;

        bubble_out
            << "bubble_id\tstatus\thas_reference_assignment\treference_cluster_id\treference_interval_start\t"
            << "reference_interval_end\tunique_alleles\tclusters_total\tclusters_observed\tclusters_ok\t"
            << "clusters_with_minimap_hit\tclusters_with_events\tdotplots_written\tcluster_pairwise_vcfs\t"
            << "region_vcf_rows\n";
        bubble_out
            << bubble.id << '\t'
            << report.status << '\t'
            << (report.has_reference_assignment ? 1 : 0) << '\t'
            << report.reference_cluster_id << '\t'
            << report.reference_interval_start << '\t'
            << report.reference_interval_end << '\t'
            << unique_alleles.size() << '\t'
            << report.clusters_total << '\t'
            << report.cluster_statuses.size() << '\t'
            << ok_clusters << '\t'
            << report.clusters_with_minimap_hit << '\t'
            << report.clusters_with_events << '\t'
            << dotplots << '\t'
            << cluster_vcfs << '\t'
            << report.vcf_rows.size() << '\n';

        const std::string cluster_status_path = bubble_debug_dir + "/cluster_status.tsv";
        std::ofstream cluster_out(cluster_status_path);
        if (!cluster_out) {
            throw std::runtime_error("Failed to write cluster debug status TSV: " + cluster_status_path);
        }
        cluster_out
            << "bubble_id\tcluster_id\tstatus\tis_reference_cluster\trepresentative_allele_id\t"
            << "representative_haplotype\treference_length_bp\tcluster_length_bp\tpaf_records\t"
            << "minimap_best_ok\torientation\tbest_norm_ed\tevent_count\tdotplot_written\t"
            << "pairwise_vcf_written\n";
        for (const auto& cluster_status : report.cluster_statuses) {
            cluster_out
                << bubble.id << '\t'
                << cluster_status.cluster_id << '\t'
                << cluster_status.status << '\t'
                << (cluster_status.is_reference_cluster ? 1 : 0) << '\t'
                << cluster_status.representative_allele_id << '\t'
                << cluster_status.representative_haplotype << '\t'
                << cluster_status.reference_length_bp << '\t'
                << cluster_status.cluster_length_bp << '\t'
                << cluster_status.paf_records << '\t'
                << (cluster_status.minimap_best_ok ? 1 : 0) << '\t'
                << cluster_status.orientation << '\t'
                << std::fixed << std::setprecision(6) << cluster_status.best_edit_distance_norm << '\t'
                << cluster_status.event_count << '\t'
                << (cluster_status.dotplot_written ? 1 : 0) << '\t'
                << (cluster_status.pairwise_vcf_written ? 1 : 0) << '\n';
        }
    };

    auto return_with_status = [&](const std::string& status) -> VariantBubbleReport {
        report.status = status;
        write_bubble_debug_status();
        return report;
    };

    if (unique_alleles.empty() || clusters.empty()) {
        return return_with_status("skipped:no-clusters");
    }

    struct BestPathCluster {
        std::size_t cluster_id = 0;
        bool source_to_sink = false;
        std::size_t allele_length = 0;
    };
    std::unordered_map<std::string, BestPathCluster> best_cluster_by_path;
    best_cluster_by_path.reserve(assignments.size() * 2);
    for (const auto& assignment : assignments) {
        if (assignment.unique_idx >= cluster_of_unique.size()) {
            continue;
        }
        const std::size_t cluster_id = cluster_of_unique[assignment.unique_idx];
        if (cluster_id == 0) {
            continue;
        }
        const std::size_t allele_length = unique_alleles[assignment.unique_idx].sequence_length;
        const auto it = best_cluster_by_path.find(assignment.path_name);
        if (it == best_cluster_by_path.end() || (!it->second.source_to_sink && assignment.source_to_sink)) {
            best_cluster_by_path[assignment.path_name] = {cluster_id, assignment.source_to_sink, allele_length};
        }
    }
    report.cluster_by_path.reserve(best_cluster_by_path.size());
    report.allele_length_by_path.reserve(best_cluster_by_path.size());
    for (const auto& [path_name, best] : best_cluster_by_path) {
        report.cluster_by_path[path_name] = best.cluster_id;
        report.allele_length_by_path[path_name] = best.allele_length;
    }

    std::optional<PathAssignment> reference_assignment;
    for (const auto& a : assignments) {
        if (a.path_name != reference_path) {
            continue;
        }
        if (!reference_assignment.has_value() || (!reference_assignment->source_to_sink && a.source_to_sink)) {
            reference_assignment = a;
        }
    }
    if (!reference_assignment.has_value()) {
        return return_with_status("skipped:no-reference-assignment");
    }
    report.has_reference_assignment = true;
    report.reference_interval_start = reference_assignment->interval_start;
    report.reference_interval_end = reference_assignment->interval_end;

    if (reference_assignment->unique_idx >= cluster_of_unique.size()) {
        return return_with_status("skipped:bad-reference-assignment");
    }

    const std::size_t ref_cluster_id = cluster_of_unique[reference_assignment->unique_idx];
    report.reference_cluster_id = ref_cluster_id;
    if (ref_cluster_id == 0) {
        return return_with_status("skipped:no-reference-cluster");
    }

    const std::string& ref_seq = unique_alleles[reference_assignment->unique_idx].sequence;
    if (ref_seq.empty()) {
        return return_with_status("skipped:no-reference-sequence");
    }

    std::string gtf_reference_label = reference_path;
    const ParsedReferencePath gtf_ref_meta = parse_reference_path_label(reference_path);
    if (!reference_sequence_global.empty() &&
        !ref_seq.empty() &&
        gtf_ref_meta.has_interval) {
        const std::size_t found = reference_sequence_global.find(ref_seq);
        if (found != std::string::npos) {
            const std::size_t start_1based = gtf_ref_meta.region_start_1based + found;
            const std::size_t end_1based = start_1based + ref_seq.size() - 1;
            gtf_reference_label = gtf_ref_meta.chrom + ":" +
                                  std::to_string(start_1based) + "-" +
                                  std::to_string(end_1based);
        }
    }

    std::unordered_map<std::size_t, std::pair<std::string, bool>> representative_path_by_unique;
    representative_path_by_unique.reserve(unique_alleles.size() * 2);
    for (const auto& assignment : assignments) {
        auto it = representative_path_by_unique.find(assignment.unique_idx);
        if (it == representative_path_by_unique.end() ||
            (!it->second.second && assignment.source_to_sink)) {
            representative_path_by_unique[assignment.unique_idx] =
                {assignment.path_name, assignment.source_to_sink};
        }
    }

    for (const auto& cluster : clusters) {
        std::vector<std::size_t> member_alleles;
        member_alleles.reserve(cluster.member_unique_idxs.size());
        for (const auto idx : cluster.member_unique_idxs) {
            member_alleles.push_back(unique_alleles[idx].allele_id);
        }
        const std::string member_alleles_text = join_size_t(member_alleles, ';');

        std::size_t representative_unique_idx = cluster.representative_idx;
        std::size_t rep_allele_id = unique_alleles[representative_unique_idx].allele_id;
        std::string representative_haplotype = ".";
        std::string cluster_seq = unique_alleles[representative_unique_idx].sequence;
        if (cluster.cluster_id == ref_cluster_id) {
            representative_unique_idx = reference_assignment->unique_idx;
            rep_allele_id = unique_alleles[reference_assignment->unique_idx].allele_id;
            representative_haplotype = reference_path;
            cluster_seq = ref_seq;
        } else {
            const auto path_it = representative_path_by_unique.find(representative_unique_idx);
            if (path_it != representative_path_by_unique.end()) {
                representative_haplotype = path_it->second.first;
            }
        }

        VariantBubbleReport::ClusterDebugStatus cluster_status;
        cluster_status.cluster_id = cluster.cluster_id;
        cluster_status.is_reference_cluster = (cluster.cluster_id == ref_cluster_id);
        cluster_status.representative_allele_id = rep_allele_id;
        cluster_status.representative_haplotype = representative_haplotype;
        cluster_status.reference_length_bp = ref_seq.size();
        cluster_status.cluster_length_bp = cluster_seq.size();

        std::string cluster_debug_dir;
        if (write_debug_files) {
            cluster_debug_dir = bubble_debug_dir + "/cluster_" + std::to_string(cluster.cluster_id);
            std::filesystem::create_directories(cluster_debug_dir);
        }

        if (cluster_seq.empty()) {
            cluster_status.status = "skipped:empty-cluster-sequence";
            write_cluster_status_text(cluster_debug_dir, cluster_status);
            report.cluster_statuses.push_back(std::move(cluster_status));
            continue;
        }

        const std::string query_name =
            representative_haplotype + "|bubble=" + std::to_string(bubble.id) +
            "|cluster=" + std::to_string(cluster.cluster_id) +
            "|allele=" + std::to_string(rep_allele_id);
        const std::size_t ref_allele_id = unique_alleles[reference_assignment->unique_idx].allele_id;
        const std::string reference_seq_id =
            reference_path + "|bubble=" + std::to_string(bubble.id) +
            "|cluster=" + std::to_string(ref_cluster_id) +
            "|allele=" + std::to_string(ref_allele_id);
        const std::string target_name = reference_seq_id;

        const Minimap2MappingResult map_result = map_query_to_reference_minimap2(
            query_name,
            cluster_seq,
            target_name,
            ref_seq,
            minimap_preset,
            minimap_best_n,
            minimap_emit_secondary);
        if (!map_result.best.ok) {
            cluster_status.status = "skipped:no-minimap-best-hit";
            write_cluster_status_text(cluster_debug_dir, cluster_status);
            report.cluster_statuses.push_back(std::move(cluster_status));
            continue;
        }
        cluster_status.minimap_best_ok = true;

        const bool rev_better = map_result.best.reverse;
        const std::string orientation = rev_better ? "-" : "+";
        cluster_status.orientation = orientation;
        const std::string cluster_seq_rc = reverse_complement(cluster_seq);
        const std::string& cmp_seq = rev_better ? cluster_seq_rc : cluster_seq;
        const int best_dist = std::max(0, map_result.best.edit_distance);
        const std::size_t max_len = std::max<std::size_t>(1, std::max(ref_seq.size(), cmp_seq.size()));
        const double best_norm = static_cast<double>(best_dist) / static_cast<double>(max_len);
        cluster_status.best_edit_distance_norm = best_norm;
        const long long len_delta = static_cast<long long>(cmp_seq.size()) - static_cast<long long>(ref_seq.size());
        const std::size_t lcp_bp = longest_common_prefix_bp(ref_seq, cmp_seq);
        const std::size_t lcs_bp = longest_common_suffix_bp(ref_seq, cmp_seq, lcp_bp);
        std::vector<CigarRun> best_runs;
        if (!map_result.best.cigar_extended.empty()) {
            best_runs = parse_extended_cigar(map_result.best.cigar_extended);
        }
        const std::size_t best_ref_start = std::min(map_result.best.target_start_bp, ref_seq.size());
        std::size_t best_alt_start = std::min(map_result.best.query_start_bp, cmp_seq.size());
        if (rev_better) {
            if (map_result.best.query_end_bp <= cluster_seq.size()) {
                best_alt_start = cluster_seq.size() - map_result.best.query_end_bp;
            } else {
                best_alt_start = 0;
            }
        }
        if (best_alt_start > cmp_seq.size()) {
            best_alt_start = 0;
        }

        std::vector<PafRecord> paf_records;
        if (!map_result.records.empty()) {
            paf_records = map_result.records;
        } else {
            PafRecord coarse;
            coarse.query_name = query_name;
            coarse.query_len = cluster_seq.size();
            coarse.query_start = 0;
            coarse.query_end = cluster_seq.size();
            coarse.strand = rev_better ? '-' : '+';
            coarse.target_name = target_name;
            coarse.target_len = ref_seq.size();
            coarse.target_start = 0;
            coarse.target_end = ref_seq.size();
            const std::size_t block = std::max(cluster_seq.size(), ref_seq.size());
            coarse.n_matches = static_cast<std::size_t>(std::max(0LL, static_cast<long long>(block) - best_dist));
            coarse.aln_block_len = block;
            coarse.mapq = 60;
            coarse.primary = true;
            paf_records.push_back(std::move(coarse));
        }
        cluster_status.paf_records = paf_records.size();

        if (write_debug_files) {
            const std::string reference_fasta = cluster_debug_dir + "/reference.fa";
            const std::string representative_fasta = cluster_debug_dir + "/representative.fa";
            const std::string paf_path = cluster_debug_dir + "/alignment.paf";

            {
                std::ofstream out(reference_fasta);
                if (!out) {
                    throw std::runtime_error("Failed to write reference FASTA: " + reference_fasta);
                }
                out << ">" << reference_seq_id << "\n" << ref_seq << "\n";
            }
            {
                std::ofstream out(representative_fasta);
                if (!out) {
                    throw std::runtime_error("Failed to write representative FASTA: " + representative_fasta);
                }
                out << ">" << query_name << "\n" << cluster_seq << "\n";
            }
            write_paf_records(paf_path, paf_records);
        }

        std::vector<AtomicVariantEvent> atomic_events;
        if (cluster.cluster_id != ref_cluster_id && map_result.best.ok) {
            atomic_events = extract_atomic_events_from_cigar(best_runs, cmp_seq, 1, best_ref_start, best_alt_start);
            const char orientation_char = rev_better ? '-' : '+';
            std::vector<AtomicVariantEvent> split_events = extract_split_events_from_records(
                map_result.records,
                orientation_char,
                min_sv_bp,
                100000,
                cmp_seq.size(),
                split_ins_use_geometric_svlen);
            atomic_events.insert(
                atomic_events.end(),
                std::make_move_iterator(split_events.begin()),
                std::make_move_iterator(split_events.end()));

            bool has_sv_indel = false;
            for (const auto& ev : atomic_events) {
                if (ev.event_type == "INS" || ev.event_type == "DEL") {
                    has_sv_indel = true;
                    break;
                }
            }
            if (!has_sv_indel) {
                if (len_delta >= static_cast<long long>(min_sv_bp)) {
                    const std::size_t end_keep = std::min<std::size_t>(lcs_bp, cmp_seq.size() - lcp_bp);
                    const std::size_t span = cmp_seq.size() - lcp_bp - end_keep;
                    if (span >= min_sv_bp && lcp_bp + span <= cmp_seq.size()) {
                        AtomicVariantEvent ev;
                        ev.event_type = "INS";
                        ev.ref_offset_start_bp = lcp_bp;
                        ev.ref_offset_end_bp = lcp_bp;
                        ev.alt_offset_start_bp = lcp_bp;
                        ev.alt_offset_end_bp = lcp_bp + span;
                        ev.event_bp = span;
                        ev.inserted_bp = span;
                        ev.inserted_seq = cmp_seq.substr(lcp_bp, span);
                        ev.svlen = static_cast<long long>(span);
                        ev.evidence_source = "fallback";
                        atomic_events.push_back(std::move(ev));
                    }
                } else if ((-len_delta) >= static_cast<long long>(min_sv_bp)) {
                    const std::size_t ref_end_keep = std::min<std::size_t>(lcs_bp, ref_seq.size() - lcp_bp);
                    const std::size_t span = ref_seq.size() - lcp_bp - ref_end_keep;
                    if (span >= min_sv_bp && lcp_bp + span <= ref_seq.size()) {
                        AtomicVariantEvent ev;
                        ev.event_type = "DEL";
                        ev.ref_offset_start_bp = lcp_bp;
                        ev.ref_offset_end_bp = lcp_bp + span;
                        ev.alt_offset_start_bp = lcp_bp;
                        ev.alt_offset_end_bp = lcp_bp;
                        ev.event_bp = span;
                        ev.svlen = -static_cast<long long>(span);
                        ev.cn_delta = -1;
                        ev.evidence_source = "fallback";
                        atomic_events.push_back(std::move(ev));
                    }
                }
            }

            // Keep explicit inversion context as dedicated event when reverse orientation is clearly preferred.
            if (rev_better && best_norm <= 0.25) {
                AtomicVariantEvent inv_event;
                inv_event.event_type = "INV";
                inv_event.event_subtype = ".";
                inv_event.ref_offset_start_bp = 0;
                inv_event.ref_offset_end_bp = ref_seq.size();
                inv_event.alt_offset_start_bp = 0;
                inv_event.alt_offset_end_bp = cmp_seq.size();
                inv_event.event_bp = std::max(ref_seq.size(), cmp_seq.size());
                inv_event.svlen = len_delta;
                inv_event.evidence_source = "orientation";
                atomic_events.insert(atomic_events.begin(), std::move(inv_event));
            }
        }
        atomic_events = merge_nearby_events(
            atomic_events,
            ref_seq,
            cmp_seq,
            std::max<std::size_t>(1, min_sv_bp / 5),
            min_sv_bp);
        for (auto& ev : atomic_events) {
            ev.event_subtype = ".";
            if (ev.event_type == "DEL") {
                ev.cn_delta = -1;
            } else if (ev.event_type == "INS") {
                ev.cn_delta = 1;
            } else {
                ev.cn_delta = 0;
            }
        }
        if (classify_ins) {
            classify_insertion_events(atomic_events, ref_seq, cmp_seq, true);
        }

        if (write_debug_files && !cluster_debug_dir.empty()) {
            const std::string dot_svg = cluster_debug_dir + "/dotplot.svg";
            write_dotplot_svg_from_paf(
                dot_svg,
                bubble.id,
                cluster.cluster_id,
                target_name,
                query_name,
                ref_seq.size(),
                cluster_seq.size(),
                paf_records,
                atomic_events,
                dotplot_gene_index,
                gtf_reference_label,
                dotplot_gene_matches);

            ++report.dotplots_written;
            ++report.debug_reports_written;
            cluster_status.dotplot_written = true;

            const std::string cluster_vcf_path = cluster_debug_dir + "/cluster_vs_reference.vcf";
            const ParsedReferencePath debug_ref_meta = parse_reference_path_label(gtf_reference_label);
            write_cluster_pairwise_vcf(
                cluster_vcf_path,
                debug_ref_meta,
                gtf_reference_label,
                bubble.id,
                cluster.cluster_id,
                ref_cluster_id,
                orientation,
                best_norm,
                ref_seq,
                cmp_seq,
                atomic_events);
            cluster_status.pairwise_vcf_written = true;
        }

        if (!atomic_events.empty()) {
            std::size_t event_index = 0;
            for (const auto& ev : atomic_events) {
                ++event_index;
                VariantBubbleReport::VariantRecordForVcf row;
                row.cluster_id = cluster.cluster_id;
                row.event_index = event_index;
                row.representative_allele_id = rep_allele_id;
                row.cluster_path_support = cluster.total_path_support;
                row.event_type = ev.event_type;
                row.event_subtype = ev.event_subtype;
                row.orientation = orientation;
                row.reference_length = ref_seq.size();
                row.cluster_length = cmp_seq.size();
                row.length_delta = ev.svlen;
                row.best_edit_distance_norm = best_norm;
                row.inserted_bp = ev.inserted_bp;
                row.member_alleles = member_alleles_text;
                row.ref_offset_start_bp = ev.ref_offset_start_bp;
                row.ref_offset_end_bp = ev.ref_offset_end_bp;
                row.alt_offset_start_bp = ev.alt_offset_start_bp;
                row.alt_offset_end_bp = ev.alt_offset_end_bp;
                row.event_bp = ev.event_bp;
                row.event_sequence = event_sequence_for_vcf(ev, ref_seq, cmp_seq);
                row.representative_haplotype = representative_haplotype;
                row.cn_delta = ev.cn_delta;
                row.preserve_ins_svlen = ev.preserve_ins_svlen;
                row.has_dup_evidence = ev.has_dup_evidence;
                row.dup_best_similarity = ev.dup_best_similarity;
                row.dup_ref_start_bp = ev.dup_ref_start_bp;
                row.dup_ref_end_bp = ev.dup_ref_end_bp;
                row.dup_orientation = ev.dup_orientation;
                row.dup_unit_bp = ev.dup_unit_bp;
                row.dup_ref_copy_number = ev.dup_ref_copy_number;
                row.dup_added_copies = ev.dup_added_copies;
                row.dup_alt_copy_number = ev.dup_alt_copy_number;
                row.dup_copy_ratio = ev.dup_copy_ratio;
                report.vcf_rows.push_back(std::move(row));
            }
        }
        cluster_status.event_count = atomic_events.size();
        cluster_status.status = atomic_events.empty() ? "ok:no-sv-events" : "ok:sv-events";
        write_cluster_status_text(cluster_debug_dir, cluster_status);
        report.cluster_statuses.push_back(std::move(cluster_status));
    }
    write_bubble_debug_status();
    return report;
}

std::pair<std::string, std::string> vcf_alt_and_type(
    const std::string& event_type,
    const std::string& event_subtype) {
    (void)event_subtype;
    if (event_type == "DEL") {
        return {"<DEL>", "DEL"};
    }
    if (event_type == "INS") {
        return {"<INS>", "INS"};
    }
    if (event_type == "INV") {
        return {"<INV>", "INV"};
    }
    return {"<INS>", "INS"};
}

void write_region_level_vcf(
    const std::string& output_path,
    const std::string& reference_path,
    const ParsedReferencePath& reference_meta,
    const std::string& reference_sequence,
    const std::vector<std::size_t>& reference_prefix_bp,
    const std::vector<RegionVcfBubble>& bubbles,
    const std::vector<std::string>& sample_names,
    std::size_t region_start_1based,
    std::size_t merge_window_bp,
    const std::string& merge_mode,
    std::size_t lenient_window_bp,
    double lenient_min_ref_jaccard,
    double min_seq_similarity,
    double max_seq_edit_fraction,
    std::size_t* records_written_out) {

    struct PendingVcfRow {
        std::size_t bubble_id = 0;
        std::string source;
        std::string sink;
        std::size_t cluster_id = 0;
        std::size_t event_index = 0;
        std::size_t reference_cluster_id = 0;
        std::string event_type;
        std::string event_subtype = ".";
        std::string orientation;
        long long length_delta = 0;
        double best_norm = 0.0;
        std::size_t inserted_bp = 0;
        std::size_t support = 0;
        std::string member_alleles;
        std::size_t pos_1based = 1;
        std::size_t end_1based = 1;
        const std::unordered_map<std::string, std::size_t>* cluster_by_path = nullptr;
        std::string ref_base = "N";
        std::string event_sequence;
        bool preserve_ins_svlen = false;
        bool has_dup_evidence = false;
        double dup_best_similarity = 0.0;
        std::size_t dup_ref_start_bp = 0;
        std::size_t dup_ref_end_bp = 0;
        char dup_orientation = '+';
        std::size_t dup_unit_bp = 0;
        std::size_t dup_ref_copy_number = 0;
        std::size_t dup_added_copies = 0;
        std::size_t dup_alt_copy_number = 0;
        double dup_copy_ratio = 0.0;
        std::string pangene_cn_delta;
        std::string pangene_gain_genes;
        std::string pangene_loss_genes;
        std::size_t pangene_gain_copies = 0;
        std::size_t pangene_loss_copies = 0;
    };

    std::string merge_mode_lower = merge_mode;
    std::transform(
        merge_mode_lower.begin(),
        merge_mode_lower.end(),
        merge_mode_lower.begin(),
        [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    const bool lenient_merge = (merge_mode_lower == "lenient");
    const std::size_t lenient_window =
        std::max<std::size_t>(merge_window_bp, lenient_window_bp);
    const double bounded_edit_cap =
        std::max(0.0, std::min(max_seq_edit_fraction, 1.0 - min_seq_similarity));

    auto sequence_similarity = [&](const std::string& lhs, const std::string& rhs) -> double {
        if (lhs.empty() && rhs.empty()) {
            return 1.0;
        }
        if (lhs.empty() || rhs.empty()) {
            return 0.0;
        }
        const std::size_t max_len = std::max(lhs.size(), rhs.size());
        const int max_edits = static_cast<int>(std::llround(bounded_edit_cap * static_cast<double>(max_len)));
        const int d = bounded_levenshtein_distance(lhs, rhs, std::max(0, max_edits));
        if (d > max_edits) {
            return 0.0;
        }
        return 1.0 - (static_cast<double>(d) / static_cast<double>(std::max<std::size_t>(1, max_len)));
    };
    auto ref_jaccard = [](const PendingVcfRow& lhs, const PendingVcfRow& rhs) -> double {
        const std::size_t a0 = std::min(lhs.pos_1based, lhs.end_1based);
        const std::size_t a1 = std::max(lhs.pos_1based, lhs.end_1based);
        const std::size_t b0 = std::min(rhs.pos_1based, rhs.end_1based);
        const std::size_t b1 = std::max(rhs.pos_1based, rhs.end_1based);
        if (a1 <= a0 || b1 <= b0) {
            return 0.0;
        }
        const std::size_t lo = std::max(a0, b0);
        const std::size_t hi = std::min(a1, b1);
        if (hi <= lo) {
            return 0.0;
        }
        const double inter = static_cast<double>(hi - lo);
        const double uni = static_cast<double>((a1 - a0) + (b1 - b0) - (hi - lo));
        if (uni <= 0.0) {
            return 0.0;
        }
        return inter / uni;
    };

    std::vector<PendingVcfRow> rows;
    rows.reserve(bubbles.size() * 2);
    for (const auto& bubble : bubbles) {
        const std::size_t lo_step = std::min(bubble.reference_interval_start, bubble.reference_interval_end);
        const std::size_t hi_step = std::max(bubble.reference_interval_start, bubble.reference_interval_end);
        if (lo_step + 1 >= reference_prefix_bp.size() || hi_step + 1 >= reference_prefix_bp.size()) {
            continue;
        }

        const std::size_t local_start_bp = reference_prefix_bp[lo_step];
        const std::size_t local_end_bp_exclusive = reference_prefix_bp[hi_step + 1];
        if (local_end_bp_exclusive < local_start_bp) {
            continue;
        }

        for (const auto& row : bubble.variant_rows) {
            PendingVcfRow pending;
            pending.bubble_id = bubble.bubble_id;
            pending.source = bubble.source;
            pending.sink = bubble.sink;
            pending.cluster_id = row.cluster_id;
            pending.event_index = row.event_index;
            pending.reference_cluster_id = bubble.reference_cluster_id;
            pending.event_type = row.event_type;
            pending.event_subtype = row.event_subtype;
            pending.orientation = row.orientation;
            pending.length_delta = row.length_delta;
            pending.best_norm = row.best_edit_distance_norm;
            pending.inserted_bp = row.inserted_bp;
            pending.support = row.cluster_path_support;
            pending.member_alleles = row.member_alleles;
            pending.preserve_ins_svlen = row.preserve_ins_svlen;
            pending.has_dup_evidence = row.has_dup_evidence;
            pending.dup_best_similarity = row.dup_best_similarity;
            pending.dup_ref_start_bp = row.dup_ref_start_bp;
            pending.dup_ref_end_bp = row.dup_ref_end_bp;
            pending.dup_orientation = row.dup_orientation;
            pending.dup_unit_bp = row.dup_unit_bp;
            pending.dup_ref_copy_number = row.dup_ref_copy_number;
            pending.dup_added_copies = row.dup_added_copies;
            pending.dup_alt_copy_number = row.dup_alt_copy_number;
            pending.dup_copy_ratio = row.dup_copy_ratio;
            pending.pangene_cn_delta = row.pangene_cn_delta;
            pending.pangene_gain_genes = row.pangene_gain_genes;
            pending.pangene_loss_genes = row.pangene_loss_genes;
            pending.pangene_gain_copies = row.pangene_gain_copies;
            pending.pangene_loss_copies = row.pangene_loss_copies;

            const std::size_t start_off = std::min(row.ref_offset_start_bp, row.ref_offset_end_bp);
            const std::size_t end_off = std::max(row.ref_offset_start_bp, row.ref_offset_end_bp);
            const std::size_t bounded_start_off = std::min(start_off, local_end_bp_exclusive - local_start_bp);
            const std::size_t bounded_end_off = std::min(end_off, local_end_bp_exclusive - local_start_bp);

            std::size_t pos_1based = region_start_1based + local_start_bp + bounded_start_off;
            if (pos_1based == 0) {
                pos_1based = 1;
            }
            std::size_t end_1based = region_start_1based + local_start_bp + bounded_end_off;
            if (end_1based < pos_1based) {
                end_1based = pos_1based;
            }
            pending.pos_1based = pos_1based;
            pending.end_1based = end_1based;
            pending.cluster_by_path = &bubble.cluster_by_path;
            pending.ref_base = "N";
            const std::size_t ref_base_off = local_start_bp + bounded_start_off;
            if (!reference_sequence.empty() && ref_base_off < reference_sequence.size()) {
                pending.ref_base = std::string(1, vcf_base(reference_sequence[ref_base_off]));
            }
            pending.event_sequence = row.event_sequence;
            if (pending.event_sequence.empty()) {
                if ((row.event_type == "DEL" || row.event_type == "INV") &&
                    row.ref_offset_end_bp > row.ref_offset_start_bp &&
                    (local_start_bp + row.ref_offset_end_bp) <= reference_sequence.size()) {
                    pending.event_sequence = reference_sequence.substr(
                        local_start_bp + row.ref_offset_start_bp,
                        row.ref_offset_end_bp - row.ref_offset_start_bp);
                }
            }
            rows.push_back(std::move(pending));
        }
    }

    std::sort(rows.begin(), rows.end(), [](const PendingVcfRow& lhs, const PendingVcfRow& rhs) {
        if (lhs.bubble_id != rhs.bubble_id) {
            return lhs.bubble_id < rhs.bubble_id;
        }
        if (lhs.pos_1based != rhs.pos_1based) {
            return lhs.pos_1based < rhs.pos_1based;
        }
        if (lhs.end_1based != rhs.end_1based) {
            return lhs.end_1based < rhs.end_1based;
        }
        if (lhs.cluster_id != rhs.cluster_id) {
            return lhs.cluster_id < rhs.cluster_id;
        }
        return lhs.event_index < rhs.event_index;
    });

    struct MergedVcfGroup {
        PendingVcfRow seed;
        std::vector<const PendingVcfRow*> members;
        std::unordered_set<std::size_t> carrier_clusters;
        std::size_t support_total = 0;
        std::vector<std::string> member_ids;
    };

    auto mergeable = [&](const MergedVcfGroup& group, const PendingVcfRow& row) {
        const auto& seed = group.seed;
        if (seed.bubble_id != row.bubble_id) {
            return false;
        }
        if (seed.event_type != row.event_type) {
            return false;
        }
        if (seed.event_subtype != row.event_subtype &&
            (seed.event_subtype != "." || row.event_subtype != ".")) {
            return false;
        }
        const std::size_t pos_gap = (seed.pos_1based > row.pos_1based)
                                        ? (seed.pos_1based - row.pos_1based)
                                        : (row.pos_1based - seed.pos_1based);
        const bool strict_pos_ok = (pos_gap <= merge_window_bp);
        if (seed.event_type != "INS") {
            const std::size_t end_gap = (seed.end_1based > row.end_1based)
                                            ? (seed.end_1based - row.end_1based)
                                            : (row.end_1based - seed.end_1based);
            const bool strict_end_ok = (end_gap <= merge_window_bp);
            if (!(strict_pos_ok && strict_end_ok)) {
                if (!lenient_merge) {
                    return false;
                }
                if (pos_gap > lenient_window) {
                    return false;
                }
                if (ref_jaccard(seed, row) < lenient_min_ref_jaccard) {
                    return false;
                }
            }
        } else if (!strict_pos_ok) {
            if (!lenient_merge || pos_gap > lenient_window) {
                return false;
            }
        }
        if (!seed.event_sequence.empty() && !row.event_sequence.empty()) {
            const double sim = sequence_similarity(seed.event_sequence, row.event_sequence);
            if (sim < min_seq_similarity) {
                return false;
            }
        }
        return true;
    };

    std::vector<MergedVcfGroup> merged;
    merged.reserve(rows.size());
    for (const auto& row : rows) {
        std::size_t idx = merged.size();
        for (std::size_t i = 0; i < merged.size(); ++i) {
            if (mergeable(merged[i], row)) {
                idx = i;
                break;
            }
        }
        if (idx == merged.size()) {
            MergedVcfGroup group;
            group.seed = row;
            group.members.push_back(&row);
            group.carrier_clusters.insert(row.cluster_id);
            group.support_total = row.support;
            group.member_ids.push_back("C" + std::to_string(row.cluster_id) + "E" + std::to_string(row.event_index));
            merged.push_back(std::move(group));
            continue;
        }

        auto& group = merged[idx];
        group.members.push_back(&row);
        group.carrier_clusters.insert(row.cluster_id);
        group.support_total += row.support;
        group.member_ids.push_back("C" + std::to_string(row.cluster_id) + "E" + std::to_string(row.event_index));
        group.seed.preserve_ins_svlen = group.seed.preserve_ins_svlen || row.preserve_ins_svlen;
        if (row.pos_1based < group.seed.pos_1based) {
            group.seed.pos_1based = row.pos_1based;
            group.seed.ref_base = row.ref_base;
        }
        group.seed.end_1based = std::max(group.seed.end_1based, row.end_1based);
        group.seed.best_norm = std::min(group.seed.best_norm, row.best_norm);
        if ((!row.pangene_cn_delta.empty() && group.seed.pangene_cn_delta.empty()) ||
            (row.pangene_gain_copies + row.pangene_loss_copies >
             group.seed.pangene_gain_copies + group.seed.pangene_loss_copies)) {
            group.seed.pangene_cn_delta = row.pangene_cn_delta;
            group.seed.pangene_gain_genes = row.pangene_gain_genes;
            group.seed.pangene_loss_genes = row.pangene_loss_genes;
            group.seed.pangene_gain_copies = row.pangene_gain_copies;
            group.seed.pangene_loss_copies = row.pangene_loss_copies;
        }
        if ((!group.seed.has_dup_evidence && row.has_dup_evidence) ||
            (group.seed.has_dup_evidence && row.has_dup_evidence &&
             row.dup_best_similarity > group.seed.dup_best_similarity)) {
            group.seed.event_subtype = row.event_subtype;
            group.seed.has_dup_evidence = row.has_dup_evidence;
            group.seed.dup_best_similarity = row.dup_best_similarity;
            group.seed.dup_ref_start_bp = row.dup_ref_start_bp;
            group.seed.dup_ref_end_bp = row.dup_ref_end_bp;
            group.seed.dup_orientation = row.dup_orientation;
            group.seed.dup_unit_bp = row.dup_unit_bp;
            group.seed.dup_ref_copy_number = row.dup_ref_copy_number;
            group.seed.dup_added_copies = row.dup_added_copies;
            group.seed.dup_alt_copy_number = row.dup_alt_copy_number;
            group.seed.dup_copy_ratio = row.dup_copy_ratio;
        }
        if (row.event_type == "INS" &&
            row.preserve_ins_svlen &&
            row.length_delta > group.seed.length_delta) {
            group.seed.length_delta = row.length_delta;
            group.seed.inserted_bp = row.inserted_bp;
        }
        if (row.event_sequence.size() > group.seed.event_sequence.size()) {
            group.seed.event_sequence = row.event_sequence;
            group.seed.length_delta = row.length_delta;
            group.seed.inserted_bp = row.inserted_bp;
            group.seed.preserve_ins_svlen = row.preserve_ins_svlen;
            group.seed.event_subtype = row.event_subtype;
            group.seed.has_dup_evidence = row.has_dup_evidence;
            group.seed.dup_best_similarity = row.dup_best_similarity;
            group.seed.dup_ref_start_bp = row.dup_ref_start_bp;
            group.seed.dup_ref_end_bp = row.dup_ref_end_bp;
            group.seed.dup_orientation = row.dup_orientation;
            group.seed.dup_unit_bp = row.dup_unit_bp;
            group.seed.dup_ref_copy_number = row.dup_ref_copy_number;
            group.seed.dup_added_copies = row.dup_added_copies;
            group.seed.dup_alt_copy_number = row.dup_alt_copy_number;
            group.seed.dup_copy_ratio = row.dup_copy_ratio;
            group.seed.pangene_cn_delta = row.pangene_cn_delta;
            group.seed.pangene_gain_genes = row.pangene_gain_genes;
            group.seed.pangene_loss_genes = row.pangene_loss_genes;
            group.seed.pangene_gain_copies = row.pangene_gain_copies;
            group.seed.pangene_loss_copies = row.pangene_loss_copies;
        }
    }

    std::sort(merged.begin(), merged.end(), [](const MergedVcfGroup& lhs, const MergedVcfGroup& rhs) {
        if (lhs.seed.bubble_id != rhs.seed.bubble_id) {
            return lhs.seed.bubble_id < rhs.seed.bubble_id;
        }
        if (lhs.seed.pos_1based != rhs.seed.pos_1based) {
            return lhs.seed.pos_1based < rhs.seed.pos_1based;
        }
        return lhs.seed.end_1based < rhs.seed.end_1based;
    });

    std::ofstream out(output_path);
    if (!out) {
        throw std::runtime_error("Failed to write region-level VCF: " + output_path);
    }

    out << "##fileformat=VCFv4.2\n";
    out << "##source=panvar\n";
    out << "##reference=" << reference_path << "\n";
    out << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant\">\n";
    out << "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Structural variant type\">\n";
    out << "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length difference ALT-REF\">\n";
    out << "##INFO=<ID=BUBBLE_ID,Number=1,Type=Integer,Description=\"panvar bubble identifier\">\n";
    out << "##INFO=<ID=REF_CLUSTER_ID,Number=1,Type=Integer,Description=\"Reference cluster identifier within bubble\">\n";
    out << "##INFO=<ID=EVENT,Number=1,Type=String,Description=\"panvar event type\">\n";
    out << "##INFO=<ID=INS_SUBTYPE,Number=1,Type=String,Description=\"INS subtype (NOVEL or DUP-like subtype)\">\n";
    out << "##INFO=<ID=ORIENT,Number=1,Type=String,Description=\"Best orientation of cluster vs reference\">\n";
    out << "##INFO=<ID=BEST_NORM_ED,Number=1,Type=Float,Description=\"Best normalized edit distance to reference cluster\">\n";
    out << "##INFO=<ID=INS_BP,Number=1,Type=Integer,Description=\"Inserted bp size\">\n";
    out << "##INFO=<ID=DUP_SIM,Number=1,Type=Float,Description=\"Best similarity of inserted sequence to duplicated source candidate\">\n";
    out << "##INFO=<ID=DUP_REF_START,Number=1,Type=Integer,Description=\"1-based start of duplicated source interval on reference\">\n";
    out << "##INFO=<ID=DUP_REF_END,Number=1,Type=Integer,Description=\"1-based end of duplicated source interval on reference\">\n";
    out << "##INFO=<ID=DUP_ORIENT,Number=1,Type=String,Description=\"Orientation of duplicated source match (+/-)\">\n";
    out << "##INFO=<ID=DUP_UNIT_BP,Number=1,Type=Integer,Description=\"Duplicated source unit length in bp\">\n";
    out << "##INFO=<ID=DUP_REF_CN,Number=1,Type=Integer,Description=\"Estimated reference copy number of duplication unit\">\n";
    out << "##INFO=<ID=DUP_ALT_CN,Number=1,Type=Integer,Description=\"Estimated ALT copy number of duplication unit\">\n";
    out << "##INFO=<ID=DUP_ADDED,Number=1,Type=Integer,Description=\"Estimated number of added copies in ALT\">\n";
    out << "##INFO=<ID=DUP_COPY_RATIO,Number=1,Type=Float,Description=\"Estimated ALT/REF copy ratio for duplication unit\">\n";
    out << "##INFO=<ID=PANGENE_CN_DELTA,Number=1,Type=String,Description=\"pangene gene copy deltas for this event cluster (gene:ref>alt, comma-separated)\">\n";
    out << "##INFO=<ID=PANGENE_GAIN_GENES,Number=1,Type=String,Description=\"pangene genes with copy gains in event cluster vs reference (gene(+delta), comma-separated)\">\n";
    out << "##INFO=<ID=PANGENE_LOSS_GENES,Number=1,Type=String,Description=\"pangene genes with copy losses in event cluster vs reference (gene(-delta), comma-separated)\">\n";
    out << "##INFO=<ID=PANGENE_GAIN_COPIES,Number=1,Type=Integer,Description=\"Total gained gene copies from pangene comparison\">\n";
    out << "##INFO=<ID=PANGENE_LOSS_COPIES,Number=1,Type=Integer,Description=\"Total lost gene copies from pangene comparison\">\n";
    out << "##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description=\"Total path support for merged event\">\n";
    out << "##INFO=<ID=CLUSTERS,Number=.,Type=Integer,Description=\"Cluster IDs carrying this merged event\">\n";
    out << "##INFO=<ID=MERGED_EVENTS,Number=.,Type=String,Description=\"Merged member event IDs (C<cluster>E<event>)\">\n";
    out << "##INFO=<ID=INSSEQ,Number=1,Type=String,Description=\"Inserted sequence\">\n";
    out << "##INFO=<ID=DELSEQ,Number=1,Type=String,Description=\"Deleted sequence from reference\">\n";
    out << "##INFO=<ID=INVSEQ,Number=1,Type=String,Description=\"Inverted reference sequence\">\n";
    out << "##INFO=<ID=BUBBLE_SOURCE,Number=1,Type=String,Description=\"Bubble source node id\">\n";
    out << "##INFO=<ID=BUBBLE_SINK,Number=1,Type=String,Description=\"Bubble sink node id\">\n";
    out << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Haplotype genotype (0=reference-like,1=event carrier,.=missing)\">\n";
    out << "##FORMAT=<ID=BC,Number=1,Type=Integer,Description=\"Assigned bubble cluster id for this sample at this bubble (missing if unassigned)\">\n";
    out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    for (const auto& sample : sample_names) {
        out << '\t' << sample;
    }
    out << '\n';

    std::size_t records_written = 0;
    for (std::size_t i = 0; i < merged.size(); ++i) {
        const auto& group = merged[i];
        const auto [alt, svtype] = vcf_alt_and_type(group.seed.event_type, ".");
        std::size_t info_end = group.seed.end_1based;
        if (svtype == "INS") {
            info_end = group.seed.pos_1based;
        }

        long long svlen = group.seed.length_delta;
        if (!group.seed.event_sequence.empty()) {
            if (svtype == "DEL") {
                svlen = -static_cast<long long>(group.seed.event_sequence.size());
            } else if (svtype == "INS") {
                if (!group.seed.preserve_ins_svlen) {
                    svlen = static_cast<long long>(group.seed.event_sequence.size());
                } else if (svlen <= 0) {
                    svlen = static_cast<long long>(group.seed.event_sequence.size());
                }
            } else if (svtype == "INV") {
                svlen = static_cast<long long>(group.seed.event_sequence.size());
            }
        }

        std::vector<std::size_t> cluster_list(group.carrier_clusters.begin(), group.carrier_clusters.end());
        std::sort(cluster_list.begin(), cluster_list.end());
        std::ostringstream clusters_csv;
        for (std::size_t c = 0; c < cluster_list.size(); ++c) {
            if (c > 0) {
                clusters_csv << ',';
            }
            clusters_csv << cluster_list[c];
        }
        std::ostringstream merged_events_csv;
        for (std::size_t m = 0; m < group.member_ids.size(); ++m) {
            if (m > 0) {
                merged_events_csv << ',';
            }
            merged_events_csv << group.member_ids[m];
        }

        const std::string record_id =
            "B" + std::to_string(group.seed.bubble_id) +
            "_M" + std::to_string(i + 1);

        out << reference_meta.chrom << '\t'
            << group.seed.pos_1based << '\t'
            << record_id << '\t'
            << group.seed.ref_base << '\t'
            << alt << '\t'
            << ".\tPASS\t";
        out << "END=" << info_end
            << ";SVTYPE=" << svtype
            << ";SVLEN=" << svlen
            << ";BUBBLE_ID=" << group.seed.bubble_id
            << ";REF_CLUSTER_ID=" << group.seed.reference_cluster_id
            << ";EVENT=" << group.seed.event_type
            << ";INS_SUBTYPE=" << group.seed.event_subtype
            << ";ORIENT=" << group.seed.orientation
            << ";BEST_NORM_ED=" << std::fixed << std::setprecision(6) << group.seed.best_norm
            << ";INS_BP=" << group.seed.inserted_bp
            << ";SUPPORT=" << group.support_total
            << ";CLUSTERS=" << clusters_csv.str()
            << ";MERGED_EVENTS=" << merged_events_csv.str();
        if (group.seed.has_dup_evidence) {
            const std::size_t dup_start_1based = region_start_1based + group.seed.dup_ref_start_bp;
            const std::size_t dup_end_1based = region_start_1based + group.seed.dup_ref_end_bp;
            out << ";DUP_SIM=" << std::fixed << std::setprecision(6) << group.seed.dup_best_similarity
                << ";DUP_REF_START=" << dup_start_1based
                << ";DUP_REF_END=" << dup_end_1based
                << ";DUP_ORIENT=" << group.seed.dup_orientation
                << ";DUP_UNIT_BP=" << group.seed.dup_unit_bp
                << ";DUP_REF_CN=" << group.seed.dup_ref_copy_number
                << ";DUP_ALT_CN=" << group.seed.dup_alt_copy_number
                << ";DUP_ADDED=" << group.seed.dup_added_copies
                << ";DUP_COPY_RATIO=" << std::fixed << std::setprecision(6) << group.seed.dup_copy_ratio;
        }
        if (!group.seed.pangene_cn_delta.empty()) {
            out << ";PANGENE_CN_DELTA=" << group.seed.pangene_cn_delta
                << ";PANGENE_GAIN_COPIES=" << group.seed.pangene_gain_copies
                << ";PANGENE_LOSS_COPIES=" << group.seed.pangene_loss_copies;
            if (!group.seed.pangene_gain_genes.empty()) {
                out << ";PANGENE_GAIN_GENES=" << group.seed.pangene_gain_genes;
            }
            if (!group.seed.pangene_loss_genes.empty()) {
                out << ";PANGENE_LOSS_GENES=" << group.seed.pangene_loss_genes;
            }
        }
        if (!group.seed.event_sequence.empty()) {
            if (svtype == "INS") {
                out << ";INSSEQ=" << group.seed.event_sequence;
            } else if (svtype == "DEL") {
                out << ";DELSEQ=" << group.seed.event_sequence;
            } else if (svtype == "INV") {
                out << ";INVSEQ=" << group.seed.event_sequence;
            }
        }
        out << ";BUBBLE_SOURCE=" << group.seed.source
            << ";BUBBLE_SINK=" << group.seed.sink;
        out << "\tGT:BC";

        for (const auto& sample : sample_names) {
            if (group.seed.cluster_by_path == nullptr) {
                out << "\t.:.";
                continue;
            }
            const auto it = group.seed.cluster_by_path->find(sample);
            if (it == group.seed.cluster_by_path->end()) {
                out << "\t.:.";
                continue;
            }
            if (group.carrier_clusters.count(it->second) > 0) {
                out << "\t1:" << it->second;
            } else {
                out << "\t0:" << it->second;
            }
        }
        out << '\n';
        ++records_written;
    }

    if (records_written_out != nullptr) {
        *records_written_out = records_written;
    }
}

struct UnsupervisedLengthCnModel {
    bool valid = false;
    std::size_t bubble_id = 0;
    std::vector<std::size_t> thresholds;
    std::vector<int> bin_labels;
    std::size_t reference_bin = 0;
    double fit_r2 = 0.0;
    double fit_score = 0.0;
    std::size_t train_samples = 0;
    std::size_t range_bp = 0;
};

std::size_t classify_length_bin(std::size_t allele_length, const std::vector<std::size_t>& thresholds) {
    std::size_t idx = 0;
    while (idx < thresholds.size() && allele_length >= thresholds[idx]) {
        ++idx;
    }
    return idx;
}

std::vector<std::size_t> collect_large_gap_thresholds(
    const std::vector<std::size_t>& unique_lengths,
    std::size_t min_sv_bp,
    std::size_t max_thresholds,
    std::size_t* range_bp_out) {

    if (range_bp_out != nullptr) {
        *range_bp_out = 0;
    }
    if (unique_lengths.size() < 2) {
        return {};
    }
    const std::size_t range_bp = unique_lengths.back() - unique_lengths.front();
    if (range_bp_out != nullptr) {
        *range_bp_out = range_bp;
    }
    if (range_bp == 0) {
        return {};
    }

    const double gap_cutoff = std::max<double>(
        static_cast<double>(std::max<std::size_t>(1, min_sv_bp)),
        0.10 * static_cast<double>(range_bp));

    std::vector<std::pair<std::size_t, std::size_t>> candidates;
    candidates.reserve(unique_lengths.size());
    for (std::size_t i = 0; i + 1 < unique_lengths.size(); ++i) {
        const std::size_t lo = unique_lengths[i];
        const std::size_t hi = unique_lengths[i + 1];
        if (hi <= lo) {
            continue;
        }
        const std::size_t gap = hi - lo;
        if (static_cast<double>(gap) < gap_cutoff) {
            continue;
        }
        candidates.push_back({gap, lo + (gap / 2)});
    }
    if (candidates.empty()) {
        return {};
    }

    if (candidates.size() > max_thresholds) {
        std::sort(candidates.begin(), candidates.end(), [](const auto& lhs, const auto& rhs) {
            if (lhs.first != rhs.first) {
                return lhs.first > rhs.first;
            }
            return lhs.second < rhs.second;
        });
        candidates.resize(max_thresholds);
    }

    std::vector<std::size_t> thresholds;
    thresholds.reserve(candidates.size());
    for (const auto& item : candidates) {
        thresholds.push_back(item.second);
    }
    std::sort(thresholds.begin(), thresholds.end());
    thresholds.erase(std::unique(thresholds.begin(), thresholds.end()), thresholds.end());
    return thresholds;
}

double unsupervised_partition_r2(
    const std::vector<std::size_t>& lengths,
    const std::vector<std::size_t>& thresholds,
    std::vector<std::size_t>* bin_counts_out) {

    if (lengths.empty()) {
        if (bin_counts_out != nullptr) {
            bin_counts_out->clear();
        }
        return 0.0;
    }

    const std::size_t bin_count = thresholds.size() + 1;
    std::vector<std::size_t> counts(bin_count, 0);
    std::vector<long double> sums(bin_count, 0.0L);
    std::vector<long double> sums_sq(bin_count, 0.0L);

    long double total_sum = 0.0L;
    long double total_sq = 0.0L;
    for (const auto len : lengths) {
        const std::size_t bin = classify_length_bin(len, thresholds);
        const long double value = static_cast<long double>(len);
        counts[bin] += 1;
        sums[bin] += value;
        sums_sq[bin] += value * value;
        total_sum += value;
        total_sq += value * value;
    }

    const long double n = static_cast<long double>(lengths.size());
    const long double overall_mean = total_sum / std::max<long double>(1.0L, n);
    const long double total_sst = std::max<long double>(0.0L, total_sq - n * overall_mean * overall_mean);
    if (total_sst <= 0.0L) {
        if (bin_counts_out != nullptr) {
            *bin_counts_out = std::move(counts);
        }
        return 0.0;
    }

    long double sse = 0.0L;
    for (std::size_t b = 0; b < bin_count; ++b) {
        if (counts[b] == 0) {
            continue;
        }
        const long double n_b = static_cast<long double>(counts[b]);
        const long double mean_b = sums[b] / n_b;
        const long double sse_b = std::max<long double>(0.0L, sums_sq[b] - n_b * mean_b * mean_b);
        sse += sse_b;
    }

    if (bin_counts_out != nullptr) {
        *bin_counts_out = std::move(counts);
    }
    const long double r2 = 1.0L - (sse / total_sst);
    if (!std::isfinite(static_cast<double>(r2))) {
        return 0.0;
    }
    return std::clamp(static_cast<double>(r2), 0.0, 1.0);
}

UnsupervisedLengthCnModel train_unsupervised_length_cn_model(
    const std::vector<RegionVcfBubble>& bubbles,
    const std::string& reference_path,
    int cn_baseline,
    std::size_t min_sv_bp) {

    UnsupervisedLengthCnModel best;
    const std::size_t min_samples = 12;
    const std::size_t max_thresholds = 8;
    const std::size_t min_model_range = std::max<std::size_t>(200, min_sv_bp * 4);

    for (const auto& bubble : bubbles) {
        const auto ref_it = bubble.allele_length_by_path.find(reference_path);
        if (ref_it == bubble.allele_length_by_path.end()) {
            continue;
        }

        std::vector<std::size_t> lengths;
        lengths.reserve(bubble.allele_length_by_path.size());
        std::vector<std::size_t> unique_lengths;
        unique_lengths.reserve(bubble.allele_length_by_path.size());
        for (const auto& [sample, len_bp] : bubble.allele_length_by_path) {
            (void)sample;
            lengths.push_back(len_bp);
            unique_lengths.push_back(len_bp);
        }
        if (lengths.size() < min_samples) {
            continue;
        }
        std::sort(unique_lengths.begin(), unique_lengths.end());
        unique_lengths.erase(std::unique(unique_lengths.begin(), unique_lengths.end()), unique_lengths.end());
        if (unique_lengths.size() < 2) {
            continue;
        }

        std::size_t range_bp = 0;
        std::vector<std::size_t> thresholds =
            collect_large_gap_thresholds(unique_lengths, min_sv_bp, max_thresholds, &range_bp);
        if (thresholds.empty() || range_bp < min_model_range) {
            continue;
        }

        std::vector<std::size_t> bin_counts;
        const double fit_r2 = unsupervised_partition_r2(lengths, thresholds, &bin_counts);
        if (!std::isfinite(fit_r2)) {
            continue;
        }
        const std::size_t min_major_bin_samples = std::max<std::size_t>(3, lengths.size() / 50);
        std::size_t major_bins = 0;
        std::size_t sparse_edge_bins = 0;
        std::size_t sparse_internal_bins = 0;
        std::size_t major_samples = 0;
        for (std::size_t b = 0; b < bin_counts.size(); ++b) {
            const std::size_t count = bin_counts[b];
            if (count >= min_major_bin_samples) {
                ++major_bins;
                major_samples += count;
                continue;
            }
            if (b == 0 || b + 1 == bin_counts.size()) {
                ++sparse_edge_bins;
            } else {
                ++sparse_internal_bins;
            }
        }
        if (fit_r2 < 0.60 || major_bins < 2 || sparse_internal_bins > 0) {
            continue;
        }

        const std::size_t reference_bin = classify_length_bin(ref_it->second, thresholds);
        if (reference_bin >= bin_counts.size()) {
            continue;
        }

        const double major_fraction =
            static_cast<double>(major_samples) /
            static_cast<double>(std::max<std::size_t>(1, lengths.size()));
        std::vector<int> bin_labels(bin_counts.size(), cn_baseline);
        for (std::size_t b = 0; b < bin_labels.size(); ++b) {
            const int delta = static_cast<int>(b) - static_cast<int>(reference_bin);
            bin_labels[b] = std::max(0, cn_baseline + delta);
        }

        const double fit_score =
            (2.20 * fit_r2) +
            (0.30 * static_cast<double>(major_bins - 1)) +
            (0.25 * std::log10(static_cast<double>(range_bp) + 1.0)) +
            (0.10 * std::log10(static_cast<double>(lengths.size()) + 1.0)) +
            (0.20 * major_fraction) -
            (0.05 * static_cast<double>(sparse_edge_bins));

        UnsupervisedLengthCnModel local;
        local.valid = true;
        local.bubble_id = bubble.bubble_id;
        local.thresholds = std::move(thresholds);
        local.bin_labels = std::move(bin_labels);
        local.reference_bin = reference_bin;
        local.fit_r2 = fit_r2;
        local.fit_score = fit_score;
        local.train_samples = lengths.size();
        local.range_bp = range_bp;

        if (!best.valid ||
            local.fit_score > best.fit_score + 1e-12 ||
            (std::abs(local.fit_score - best.fit_score) <= 1e-12 &&
             local.fit_r2 > best.fit_r2 + 1e-12) ||
            (std::abs(local.fit_score - best.fit_score) <= 1e-12 &&
             std::abs(local.fit_r2 - best.fit_r2) <= 1e-12 &&
             local.train_samples > best.train_samples) ||
            (std::abs(local.fit_score - best.fit_score) <= 1e-12 &&
             std::abs(local.fit_r2 - best.fit_r2) <= 1e-12 &&
             local.train_samples == best.train_samples &&
             local.bubble_id < best.bubble_id)) {
            best = std::move(local);
        }
    }

    return best;
}

[[maybe_unused]] void write_cn_reports(
    const std::string& cn_by_sample_path,
    const std::string& cn_by_event_path,
    int cn_baseline,
    std::size_t min_sv_bp,
    const std::string& reference_path,
    const std::vector<std::string>& sample_names,
    const std::vector<RegionVcfBubble>& bubbles,
    std::size_t* cn_sample_rows_out,
    std::size_t* cn_event_rows_out) {

    std::ofstream sample_out(cn_by_sample_path);
    if (!sample_out) {
        throw std::runtime_error("Failed to write CN sample report TSV: " + cn_by_sample_path);
    }
    std::ofstream event_out(cn_by_event_path);
    if (!event_out) {
        throw std::runtime_error("Failed to write CN event report TSV: " + cn_by_event_path);
    }

    std::unordered_map<std::string, int> cn_by_sample;
    std::unordered_map<std::string, int> delta_by_sample;
    std::unordered_map<std::string, std::size_t> events_by_sample;
    cn_by_sample.reserve(sample_names.size() * 2);
    delta_by_sample.reserve(sample_names.size() * 2);
    events_by_sample.reserve(sample_names.size() * 2);
    for (const auto& sample : sample_names) {
        cn_by_sample[sample] = cn_baseline;
        delta_by_sample[sample] = 0;
        events_by_sample[sample] = 0;
    }

    event_out
        << "bubble_id\tsource\tsink\tcluster_id\tevent_index\tevent_type\tevent_subtype\t"
        << "cn_delta\tcluster_path_support\tsamples_with_event\ttotal_cn_delta\n";

    std::size_t event_rows = 0;
    for (const auto& bubble : bubbles) {
        for (const auto& row : bubble.variant_rows) {
            if (row.cn_delta == 0) {
                continue;
            }
            std::size_t samples_with_event = 0;
            long long total_cn_delta = 0;
            for (const auto& sample : sample_names) {
                const auto it = bubble.cluster_by_path.find(sample);
                if (it == bubble.cluster_by_path.end()) {
                    continue;
                }
                if (it->second != row.cluster_id) {
                    continue;
                }
                cn_by_sample[sample] += row.cn_delta;
                delta_by_sample[sample] += row.cn_delta;
                events_by_sample[sample] += 1;
                ++samples_with_event;
                total_cn_delta += row.cn_delta;
            }

            event_out
                << bubble.bubble_id << '\t'
                << bubble.source << '\t'
                << bubble.sink << '\t'
                << row.cluster_id << '\t'
                << row.event_index << '\t'
                << row.event_type << '\t'
                << row.event_subtype << '\t'
                << row.cn_delta << '\t'
                << row.cluster_path_support << '\t'
                << samples_with_event << '\t'
                << total_cn_delta << '\n';
            ++event_rows;
        }
    }

    UnsupervisedLengthCnModel cn_model_unsupervised;
    const RegionVcfBubble* model_bubble = nullptr;
    if (!reference_path.empty()) {
        cn_model_unsupervised =
            train_unsupervised_length_cn_model(bubbles, reference_path, cn_baseline, min_sv_bp);
        if (cn_model_unsupervised.valid) {
            for (const auto& bubble : bubbles) {
                if (bubble.bubble_id == cn_model_unsupervised.bubble_id) {
                    model_bubble = &bubble;
                    break;
                }
            }
        }
    }

    sample_out
        << "sample_name\tbaseline_cn\tpredicted_cn\ttotal_delta\tevents_applied\t"
        << "event_based_cn\tcalibrated_cn\tcn_model\tcn_model_bubble\tcn_model_quality\n";
    std::size_t sample_rows = 0;
    for (const auto& sample : sample_names) {
        const int event_based_cn = cn_by_sample[sample];
        std::optional<int> calibrated_cn;
        std::string model_name = "event";
        if (cn_model_unsupervised.valid && model_bubble != nullptr) {
            const auto len_it = model_bubble->allele_length_by_path.find(sample);
            if (len_it != model_bubble->allele_length_by_path.end()) {
                const std::size_t bin = classify_length_bin(len_it->second, cn_model_unsupervised.thresholds);
                if (bin < cn_model_unsupervised.bin_labels.size()) {
                    calibrated_cn = cn_model_unsupervised.bin_labels[bin];
                    model_name = "length-unsupervised";
                }
            } else {
                // Missing assignment at the calibration bubble is treated as one
                // bin below the smallest observed class (e.g. CN0 when ref CN=2).
                const int inferred =
                    cn_baseline - static_cast<int>(cn_model_unsupervised.reference_bin) - 1;
                calibrated_cn = std::max(0, inferred);
                model_name = "length-unsupervised";
            }
        }
        const int final_cn = calibrated_cn.has_value() ? *calibrated_cn : event_based_cn;

        sample_out
            << sample << '\t'
            << cn_baseline << '\t'
            << final_cn << '\t'
            << delta_by_sample[sample] << '\t'
            << events_by_sample[sample] << '\t'
            << event_based_cn << '\t';
        if (calibrated_cn.has_value()) {
            sample_out << *calibrated_cn;
        } else {
            sample_out << '.';
        }
        sample_out << '\t' << model_name << '\t';
        if (cn_model_unsupervised.valid) {
            sample_out << cn_model_unsupervised.bubble_id << '\t'
                       << std::fixed << std::setprecision(6) << cn_model_unsupervised.fit_score;
        } else {
            sample_out << ".\t.";
        }
        sample_out << '\n';
        ++sample_rows;
    }

    if (cn_sample_rows_out != nullptr) {
        *cn_sample_rows_out = sample_rows;
    }
    if (cn_event_rows_out != nullptr) {
        *cn_event_rows_out = event_rows;
    }
}

void write_headers(
    std::ofstream& clusters_out,
    std::ofstream& assignments_out,
    std::ofstream* cluster_sequences_out) {
    clusters_out
        << "bubble_id,source,sink,cluster_id,representative_allele_id,total_path_support,"
        << "member_allele_count,member_alleles\n";

    assignments_out
        << "bubble_id,source,sink,path_name,cluster_id,allele_id,allele_length,interval_start,interval_end,"
        << "source_to_sink\n";

    if (cluster_sequences_out != nullptr) {
        (*cluster_sequences_out)
            << "bubble_id,source,sink,cluster_id,representative_allele_id,"
            << "total_path_support,sequence_length,sequence\n";
    }
}

auto* kParseReferencePathLabelImpl = &parse_reference_path_label;
auto* kLoadGeneAnnotationsFromGtfImpl = &load_gene_annotations_from_gtf;
auto* kWriteVariantReportsForBubbleImpl = &write_variant_reports_for_bubble;
auto* kWriteRegionLevelVcfImpl = &write_region_level_vcf;

} // namespace

ParsedReferencePath parse_reference_path_label(const std::string& path_name) {
    return (*kParseReferencePathLabelImpl)(path_name);
}

GeneAnnotationIndex load_gene_annotations_from_gtf(const std::string& gtf_path) {
    return (*kLoadGeneAnnotationsFromGtfImpl)(gtf_path);
}

VariantBubbleReport write_variant_reports_for_bubble(
    const Bubble& bubble,
    const std::vector<UniqueAllele>& unique_alleles,
    const std::vector<Cluster>& clusters,
    const std::vector<std::size_t>& cluster_of_unique,
    const std::vector<PathAssignment>& assignments,
    const std::string& reference_path,
    const std::string& reference_sequence_global,
    std::size_t min_sv_bp,
    bool write_debug_files,
    const std::string& debug_output_dir,
    const GeneAnnotationIndex* dotplot_gene_index,
    const std::vector<std::string>& dotplot_gene_matches,
    const std::string& minimap_preset,
    std::size_t minimap_best_n,
    bool minimap_emit_secondary,
    bool split_ins_use_geometric_svlen,
    bool classify_ins) {
    return (*kWriteVariantReportsForBubbleImpl)(
        bubble,
        unique_alleles,
        clusters,
        cluster_of_unique,
        assignments,
        reference_path,
        reference_sequence_global,
        min_sv_bp,
        write_debug_files,
        debug_output_dir,
        dotplot_gene_index,
        dotplot_gene_matches,
        minimap_preset,
        minimap_best_n,
        minimap_emit_secondary,
        split_ins_use_geometric_svlen,
        classify_ins);
}

void write_region_level_vcf(
    const std::string& output_path,
    const std::string& reference_path,
    const ParsedReferencePath& reference_meta,
    const std::string& reference_sequence,
    const std::vector<std::size_t>& reference_prefix_bp,
    const std::vector<RegionVcfBubble>& bubbles,
    const std::vector<std::string>& sample_names,
    std::size_t region_start_1based,
    std::size_t merge_window_bp,
    const std::string& merge_mode,
    std::size_t lenient_window_bp,
    double lenient_min_ref_jaccard,
    double min_seq_similarity,
    double max_seq_edit_fraction,
    std::size_t* records_written_out) {
    (*kWriteRegionLevelVcfImpl)(
        output_path,
        reference_path,
        reference_meta,
        reference_sequence,
        reference_prefix_bp,
        bubbles,
        sample_names,
        region_start_1based,
        merge_window_bp,
        merge_mode,
        lenient_window_bp,
        lenient_min_ref_jaccard,
        min_seq_similarity,
        max_seq_edit_fraction,
        records_written_out);
}

void validate_call_merge_options(const AlleleCallOptions& options) {
    std::string merge_mode = options.vcf_merge_mode;
    std::transform(
        merge_mode.begin(),
        merge_mode.end(),
        merge_mode.begin(),
        [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    if (!(merge_mode == "strict" || merge_mode == "lenient")) {
        throw std::runtime_error("vcf_merge_mode must be one of: strict, lenient");
    }
    if (options.vcf_merge_lenient_window_bp == 0) {
        throw std::runtime_error("vcf_merge_lenient_window_bp must be >= 1");
    }
    if (!(options.vcf_merge_lenient_min_ref_jaccard >= 0.0 &&
          options.vcf_merge_lenient_min_ref_jaccard <= 1.0)) {
        throw std::runtime_error("vcf_merge_lenient_min_ref_jaccard must be in [0,1]");
    }
    if (!(options.vcf_merge_min_seq_similarity >= 0.0 &&
          options.vcf_merge_min_seq_similarity <= 1.0)) {
        throw std::runtime_error("vcf_merge_min_seq_similarity must be in [0,1]");
    }
    if (!(options.vcf_merge_max_seq_edit_fraction >= 0.0 &&
          options.vcf_merge_max_seq_edit_fraction <= 1.0)) {
        throw std::runtime_error("vcf_merge_max_seq_edit_fraction must be in [0,1]");
    }
}

void call_alleles_to_csv(
    const Graph& graph,
    const AlleleCallOptions& options,
    const std::string& clusters_csv_path,
    const std::string& assignments_csv_path,
    AlleleCallSummary* summary_out) {

    if (!(options.min_similarity > 0.0 && options.min_similarity <= 1.0)) {
        throw std::runtime_error("min_similarity must be in (0,1]");
    }
    validate_call_merge_options(options);

    std::ofstream clusters_out(clusters_csv_path);
    if (!clusters_out) {
        throw std::runtime_error("Failed to write allele clusters CSV: " + clusters_csv_path);
    }
    std::ofstream assignments_out(assignments_csv_path);
    if (!assignments_out) {
        throw std::runtime_error("Failed to write allele assignments CSV: " + assignments_csv_path);
    }
    std::ofstream cluster_sequences_out;
    if (options.write_cluster_sequences) {
        if (options.cluster_sequences_csv_path.empty()) {
            throw std::runtime_error(
                "--cluster-sequences-csv cannot be empty when representative sequence export is enabled");
        }
        cluster_sequences_out.open(options.cluster_sequences_csv_path);
        if (!cluster_sequences_out) {
            throw std::runtime_error(
                "Failed to write cluster sequence CSV: " + options.cluster_sequences_csv_path);
        }
    }
    write_headers(
        clusters_out,
        assignments_out,
        options.write_cluster_sequences ? &cluster_sequences_out : nullptr);

    if (options.bubbles_csv_in.empty()) {
        throw std::runtime_error(
            "Module 'allele' requires module-1 bubbles input (use --bubble-prefix-in in CLI).");
    }
    const std::vector<Bubble> bubbles = read_bubbles_csv(options.bubbles_csv_in);
    std::unordered_map<std::string, std::string> predefined_cluster_labels_by_path;
    const bool use_predefined_clusters = !options.predefined_clusters_json_path.empty();
    if (use_predefined_clusters) {
        predefined_cluster_labels_by_path =
            load_predefined_path_cluster_map_json(options.predefined_clusters_json_path);
        if (predefined_cluster_labels_by_path.empty()) {
            throw std::runtime_error(
                "Predefined clusters JSON is empty: " + options.predefined_clusters_json_path);
        }
    }

    std::vector<BubblePathIndex> path_indexes;
    path_indexes.reserve(graph.paths.size());
    for (const auto& path : graph.paths) {
        path_indexes.push_back(build_bubble_path_index(path));
    }

    AlleleCallSummary summary;
    summary.bubbles_processed = bubbles.size();
    const auto run_start = std::chrono::steady_clock::now();
    const bool need_allele_sequence =
        is_sequence_cluster_mode(options.cluster_mode) ||
        options.write_cluster_sequences ||
        options.write_region_vcf ||
        options.write_debug_reports;

    std::ofstream odgi_manifest;
    std::string odgi_manifest_path;
    const std::string odgi_input = options.odgi_input_path;
    const auto graph_paths = path_records_by_name(graph);
    const PathRecord* reference_path_record = nullptr;
    const bool need_reference_path =
        options.write_odgi_viz_inputs || options.write_region_vcf || options.write_debug_reports ||
        options.skip_bubbles_without_reference;
    std::vector<std::size_t> reference_prefix_bp;
    std::string reference_sequence;
    ParsedReferencePath reference_meta;

    if (need_reference_path) {
        if (options.reference_path.empty()) {
            throw std::runtime_error("--reference-path is required for reference-aware exports");
        }
        const auto ref_it = graph_paths.find(options.reference_path);
        if (ref_it == graph_paths.end()) {
            throw std::runtime_error("Reference path not found in graph: " + options.reference_path);
        }
        reference_path_record = ref_it->second;
        reference_prefix_bp = path_prefix_bp(*reference_path_record, graph.nodes);
        bool complete_ref_sequence = false;
        reference_sequence = spell_path_steps_sequence(graph, reference_path_record->steps, &complete_ref_sequence);
        if (!complete_ref_sequence) {
            reference_sequence.clear();
        }
        reference_meta = ::panvar::parse_reference_path_label(options.reference_path);
    }

    GeneAnnotationIndex dotplot_gene_index;
    const GeneAnnotationIndex* dotplot_gene_ptr = nullptr;
    if (options.write_debug_reports && !options.dotplot_gtf_path.empty()) {
        dotplot_gene_index = ::panvar::load_gene_annotations_from_gtf(options.dotplot_gtf_path);
        dotplot_gene_ptr = &dotplot_gene_index;
    }

    std::vector<RegionVcfBubble> region_vcf_bubbles;
    std::vector<std::string> vcf_sample_names;
    if (options.write_region_vcf) {
        region_vcf_bubbles.reserve(bubbles.size());
        std::unordered_set<std::string> seen_samples;
        seen_samples.reserve(graph.paths.size() * 2);
        for (const auto& path : graph.paths) {
            if (seen_samples.insert(path.name).second) {
                vcf_sample_names.push_back(path.name);
            }
        }
    }

    if (options.write_odgi_viz_inputs) {
        if (options.odgi_viz_out_dir.empty()) {
            throw std::runtime_error("--odgi-viz-out-dir is required when ODGI viz export is enabled");
        }
        if (odgi_input.empty()) {
            throw std::runtime_error("--odgi-input cannot be empty when ODGI viz export is enabled");
        }

        std::filesystem::create_directories(options.odgi_viz_out_dir);
        odgi_manifest_path = options.odgi_viz_out_dir + "/manifest.tsv";

        odgi_manifest.open(odgi_manifest_path);
        if (!odgi_manifest) {
            throw std::runtime_error("Failed to write ODGI manifest: " + odgi_manifest_path);
        }
        odgi_manifest
            << "bubble_id\tsource\tsink\tbubble_dir\trows\tclusters\trange_start_bp\trange_end_bp\t"
            << "paths_file\tpath_colors_file\tviz_script\tpng_file\todgi_rc\tstatus\n";
    }

    std::ofstream similarity_manifest;
    if (options.write_similarity_reports) {
        if (options.similarity_out_dir.empty()) {
            throw std::runtime_error("--similarity-out-dir is required when similarity reports are enabled");
        }
        std::filesystem::create_directories(options.similarity_out_dir);
        const std::string manifest_path = options.similarity_out_dir + "/manifest.tsv";
        similarity_manifest.open(manifest_path);
        if (!similarity_manifest) {
            throw std::runtime_error("Failed to write similarity manifest: " + manifest_path);
        }
        similarity_manifest
            << "bubble_id\tsource\tsink\tstatus\tunique_alleles\tclusters\tbubble_dir\t"
            << "alleles_tsv\tdistance_matrix_norm_tsv\tdistance_matrix_abs_tsv\t"
            << "cluster_stats_tsv\tcluster_signatures_dir\n";
    }

    if (options.write_debug_reports) {
        if (options.debug_out_dir.empty()) {
            throw std::runtime_error("--debug-out-dir cannot be empty when call debug output is enabled");
        }
        std::filesystem::create_directories(options.debug_out_dir);
    }

    auto emit_similarity_manifest_skip = [&](const Bubble& bubble, const std::string& status) {
        if (!options.write_similarity_reports) {
            return;
        }
        const SimilarityReportPaths sim_paths =
            similarity_report_paths_for_bubble(options.similarity_out_dir, bubble.id);
        similarity_manifest << bubble.id << '\t'
                            << bubble.source << '\t'
                            << bubble.sink << '\t'
                            << status << '\t'
                            << 0 << '\t'
                            << 0 << '\t'
                            << sim_paths.bubble_dir << '\t'
                            << "." << '\t'
                            << "." << '\t'
                            << "." << '\t'
                            << "." << '\t'
                            << "." << '\n';
    };

    auto emit_odgi_manifest_skip = [&](const Bubble& bubble, const std::string& status) {
        if (!options.write_odgi_viz_inputs) {
            return;
        }
        const std::string bubble_dir = options.odgi_viz_out_dir + "/bubble_" + std::to_string(bubble.id);
        odgi_manifest << bubble.id << '\t'
                      << bubble.source << '\t'
                      << bubble.sink << '\t'
                      << bubble_dir << '\t'
                      << 0 << '\t'
                      << 0 << '\t'
                      << 0 << '\t'
                      << 0 << '\t'
                      << "." << '\t'
                      << "." << '\t'
                      << "." << '\t'
                      << "." << '\t'
                      << -1 << '\t'
                      << status << '\n';
    };

    std::size_t predefined_missing_assignments_total = 0;
    std::size_t predefined_missing_paths_total = 0;

    for (std::size_t bubble_idx = 0; bubble_idx < bubbles.size(); ++bubble_idx) {
        const auto& bubble = bubbles[bubble_idx];
        const auto bubble_start = std::chrono::steady_clock::now();
        if (options.show_progress) {
            std::cerr
                << "[allele] bubble " << (bubble_idx + 1) << "/" << bubbles.size()
                << " id=" << bubble.id
                << " inside_nodes=" << bubble.inside.size()
                << "\n";
        }

        std::unordered_map<std::string, std::size_t> unique_by_signature;
        std::vector<UniqueAllele> unique_alleles;
        std::vector<PathAssignment> assignments;
        const auto extract_start = std::chrono::steady_clock::now();

        unique_by_signature.reserve(graph.paths.size() * 2);
        unique_alleles.reserve(graph.paths.size());
        assignments.reserve(graph.paths.size());
        bool has_reference_assignment = false;

        for (std::size_t path_idx = 0; path_idx < graph.paths.size(); ++path_idx) {
            const auto& path = graph.paths[path_idx];
            const auto& index = path_indexes[path_idx];

            const auto best_interval = find_best_bubble_path_interval(index, bubble);
            if (!best_interval.has_value()) {
                continue;
            }

            std::vector<PathStep> allele_steps = canonical_bubble_path_steps(path, bubble, *best_interval);
            if (allele_steps.size() < 2) {
                continue;
            }
            if (allele_steps.front().node_id != bubble.source || allele_steps.back().node_id != bubble.sink) {
                continue;
            }

            std::string signature = build_walk_signature(allele_steps);
            std::size_t unique_idx = 0;

            const auto found = unique_by_signature.find(signature);
            if (found == unique_by_signature.end()) {
                UniqueAllele allele;
                allele.steps = std::move(allele_steps);
                allele.signature = std::move(signature);
                allele.path_support = 0;

                bool has_sequence = false;
                if (need_allele_sequence) {
                    allele.sequence = spell_path_steps_sequence(graph, allele.steps, &has_sequence);
                }
                allele.compare_steps = build_walk_tokens(allele.steps);
                allele.compare_step_weights = build_walk_step_weights(allele.steps, graph.nodes);
                allele.compare_steps_weight_sum = sum_step_weights(allele.compare_step_weights);
                if (is_sequence_cluster_mode(options.cluster_mode)) {
                    if (has_sequence) {
                        allele.compare_token = allele.sequence;
                        allele.sequence_length = allele.sequence.size();
                        allele.uses_sequence_similarity = true;
                    } else {
                        allele.compare_token = allele.signature;
                        allele.sequence_length =
                            (allele.compare_steps_weight_sum > 0)
                                ? allele.compare_steps_weight_sum
                                : allele.steps.size();
                        allele.uses_sequence_similarity = false;
                    }
                } else {
                    allele.compare_token = allele.signature;
                    allele.sequence_length =
                        has_sequence
                            ? allele.sequence.size()
                            : ((allele.compare_steps_weight_sum > 0)
                                   ? allele.compare_steps_weight_sum
                                   : allele.steps.size());
                    allele.uses_sequence_similarity = false;
                }

                unique_alleles.push_back(std::move(allele));
                unique_idx = unique_alleles.size() - 1;
                unique_by_signature.emplace(unique_alleles.back().signature, unique_idx);
            } else {
                unique_idx = found->second;
            }

            unique_alleles[unique_idx].path_support += 1;

            PathAssignment assignment;
            assignment.unique_idx = unique_idx;
            assignment.path_name = path.name;
            assignment.interval_start = best_interval->left;
            assignment.interval_end = best_interval->right;
            assignment.source_to_sink = best_interval->source_to_sink;
            assignments.push_back(std::move(assignment));
            if (!options.reference_path.empty() && path.name == options.reference_path) {
                has_reference_assignment = true;
            }
        }

        if (assignments.empty()) {
            emit_similarity_manifest_skip(bubble, "skipped:no-path-assignments");
            emit_odgi_manifest_skip(bubble, "skipped:no-path-assignments");
            if (options.show_progress) {
                std::cerr
                    << "[allele] bubble " << bubble.id
                    << " skipped (no source->sink path assignments), elapsed="
                    << std::fixed << std::setprecision(2) << elapsed_seconds(bubble_start) << "s\n";
            }
            continue;
        }

        if (options.skip_bubbles_without_reference && !has_reference_assignment) {
            summary.bubbles_skipped_no_reference += 1;
            emit_similarity_manifest_skip(bubble, "skipped:no-reference-assignment-filter");
            emit_odgi_manifest_skip(bubble, "skipped:no-reference-assignment-filter");
            if (options.show_progress) {
                std::cerr
                    << "[allele] bubble " << bubble.id
                    << " skipped by --skip-no-reference-bubbles, elapsed="
                    << std::fixed << std::setprecision(2) << elapsed_seconds(bubble_start) << "s\n";
            }
            continue;
        }

        if (options.show_progress) {
            std::cerr
                << "[allele] bubble " << bubble.id
                << " extracted paths=" << assignments.size()
                << ", unique_alleles=" << unique_alleles.size()
                << ", elapsed=" << std::fixed << std::setprecision(2)
                << elapsed_seconds(extract_start) << "s\n";
        }

        summary.bubbles_with_assignments += 1;

        for (std::size_t i = 0; i < unique_alleles.size(); ++i) {
            unique_alleles[i].allele_id = i + 1;
        }

        const std::uint64_t pair_count =
            static_cast<std::uint64_t>(unique_alleles.size()) *
            static_cast<std::uint64_t>(unique_alleles.size() - 1) / 2ULL;
        bool all_sequence_tokens = true;
        std::size_t max_sequence_token_len = 0;
        for (const auto& allele : unique_alleles) {
            if (!allele.uses_sequence_similarity) {
                all_sequence_tokens = false;
                break;
            }
            max_sequence_token_len = std::max(max_sequence_token_len, allele.compare_token.size());
        }
        const long double sequence_work =
            static_cast<long double>(pair_count) *
            static_cast<long double>(std::max<std::size_t>(1, max_sequence_token_len));

        const bool by_allele_cap =
            options.max_upgma_alleles > 0 &&
            unique_alleles.size() > options.max_upgma_alleles;
        const bool by_sequence_work = sequence_work >= kAutoSequenceFastWorkThreshold;
        const bool auto_sequence_fast =
            !use_predefined_clusters &&
            options.cluster_mode == ClusterMode::Sequence &&
            options.fast_distance &&
            all_sequence_tokens &&
            (by_allele_cap || by_sequence_work);

        const bool use_sequence_fast_clustering =
            !use_predefined_clusters &&
            (options.cluster_mode == ClusterMode::SequenceFast || auto_sequence_fast);
        const bool use_threshold_graph_clustering =
            !use_predefined_clusters &&
            !use_sequence_fast_clustering &&
            options.max_upgma_alleles > 0 &&
            unique_alleles.size() > options.max_upgma_alleles;

        const bool need_distance_matrix =
            !use_predefined_clusters || options.write_similarity_reports;
        DistanceMatrices dists;
        if (need_distance_matrix) {
            const auto dist_start = std::chrono::steady_clock::now();
            dists = build_distance_matrices(
                unique_alleles,
                options.min_similarity,
                options.fast_distance,
                options.exact_distance_max_bp,
                options.threads,
                options.show_progress,
                bubble.id,
                use_threshold_graph_clustering,
                use_sequence_fast_clustering);
            if (options.show_progress) {
                std::cerr
                    << "[allele] bubble " << bubble.id
                    << " distance matrix computed, elapsed="
                    << std::fixed << std::setprecision(2) << elapsed_seconds(dist_start) << "s\n";
            }
        }

        UPGMATree tree;
        std::vector<Cluster> clusters;
        if (use_predefined_clusters) {
            std::size_t missing_assignments = 0;
            std::size_t missing_paths = 0;
            clusters = cluster_unique_alleles_from_predefined_path_labels(
                unique_alleles,
                assignments,
                predefined_cluster_labels_by_path,
                &missing_assignments,
                &missing_paths);
            predefined_missing_assignments_total += missing_assignments;
            predefined_missing_paths_total += missing_paths;
            if (options.show_progress) {
                std::cerr
                    << "[allele] bubble " << bubble.id
                    << " predefined clustering from JSON labels"
                    << " (clusters=" << clusters.size()
                    << ", unmapped_assignments=" << missing_assignments
                    << ", unmapped_paths=" << missing_paths
                    << ")\n";
            }
        } else if (use_sequence_fast_clustering) {
            if (options.show_progress) {
                std::cerr
                    << "[allele] bubble " << bubble.id
                    << " sequence fast-path clustering (sketch-estimated distances, greedy threshold assignment)";
                if (auto_sequence_fast) {
                    std::cerr
                        << " [auto-triggered: pair_count=" << pair_count
                        << ", max_token_len=" << max_sequence_token_len
                        << ", workload=" << std::fixed << std::setprecision(0)
                        << static_cast<double>(sequence_work)
                        << ", threshold=" << static_cast<double>(kAutoSequenceFastWorkThreshold)
                        << "]";
                }
                std::cerr << "\n";
            }
            clusters = cluster_unique_alleles_greedy_threshold(
                unique_alleles,
                dists.norm,
                options.min_similarity);
        } else if (use_threshold_graph_clustering) {
            if (options.show_progress) {
                std::cerr
                    << "[allele] bubble " << bubble.id
                    << " switching to threshold-graph clustering"
                    << " (unique_alleles=" << unique_alleles.size()
                    << " > max_upgma_alleles=" << options.max_upgma_alleles << ")\n";
            }
            clusters = cluster_unique_alleles_threshold_graph(
                unique_alleles,
                dists.norm,
                options.min_similarity);
        } else {
            tree = build_upgma_tree(dists.norm);
            clusters = cluster_unique_alleles_from_tree(
                unique_alleles,
                dists.norm,
                tree,
                options.min_similarity);
        }
        if (clusters.empty()) {
            emit_similarity_manifest_skip(bubble, "skipped:no-clusters");
            emit_odgi_manifest_skip(bubble, "skipped:no-clusters");
            if (options.show_progress) {
                std::cerr
                    << "[allele] bubble " << bubble.id
                    << " skipped (0 clusters), elapsed="
                    << std::fixed << std::setprecision(2) << elapsed_seconds(bubble_start) << "s\n";
            }
            continue;
        }

        summary.unique_alleles += unique_alleles.size();
        summary.clusters += clusters.size();
        summary.assignments += assignments.size();

        std::vector<std::size_t> cluster_of_unique(unique_alleles.size(), 0);
        for (const auto& cluster : clusters) {
            for (const std::size_t unique_idx : cluster.member_unique_idxs) {
                cluster_of_unique[unique_idx] = cluster.cluster_id;
            }
        }

        for (const auto& cluster : clusters) {
            const UniqueAllele& rep = unique_alleles[cluster.representative_idx];
            std::vector<std::size_t> member_alleles;
            member_alleles.reserve(cluster.member_unique_idxs.size());
            for (const std::size_t unique_idx : cluster.member_unique_idxs) {
                member_alleles.push_back(unique_alleles[unique_idx].allele_id);
            }

            clusters_out
                << bubble.id << ','
                << csv_escape(bubble.source) << ','
                << csv_escape(bubble.sink) << ','
                << cluster.cluster_id << ','
                << rep.allele_id << ','
                << cluster.total_path_support << ','
                << member_alleles.size() << ','
                << csv_escape(join_size_t(member_alleles, ';'))
                << '\n';

            if (options.write_cluster_sequences) {
                cluster_sequences_out
                    << bubble.id << ','
                    << csv_escape(bubble.source) << ','
                    << csv_escape(bubble.sink) << ','
                    << cluster.cluster_id << ','
                    << rep.allele_id << ','
                    << cluster.total_path_support << ','
                    << rep.sequence_length << ','
                    << csv_escape(rep.sequence)
                    << '\n';
                summary.cluster_sequences_written += 1;
            }
        }

        for (const auto& assignment : assignments) {
            const UniqueAllele& allele = unique_alleles[assignment.unique_idx];
            const std::size_t cluster_id = cluster_of_unique[assignment.unique_idx];
            assignments_out
                << bubble.id << ','
                << csv_escape(bubble.source) << ','
                << csv_escape(bubble.sink) << ','
                << csv_escape(assignment.path_name) << ','
                << cluster_id << ','
                << allele.allele_id << ','
                << allele.sequence_length << ','
                << assignment.interval_start << ','
                << assignment.interval_end << ','
                << (assignment.source_to_sink ? 1 : 0)
                << '\n';
        }
        clusters_out.flush();
        assignments_out.flush();
        if (options.write_cluster_sequences) {
            cluster_sequences_out.flush();
        }

        if (options.write_similarity_reports) {
            const SimilarityReportPaths sim_paths =
                similarity_report_paths_for_bubble(options.similarity_out_dir, bubble.id);

            write_similarity_reports_for_bubble(
                sim_paths,
                bubble,
                graph,
                unique_alleles,
                clusters,
                cluster_of_unique,
                dists,
                options.min_similarity);

            similarity_manifest << bubble.id << '\t'
                                << bubble.source << '\t'
                                << bubble.sink << '\t'
                                << "ok" << '\t'
                                << unique_alleles.size() << '\t'
                                << clusters.size() << '\t'
                                << sim_paths.bubble_dir << '\t'
                                << sim_paths.alleles_tsv << '\t'
                                << sim_paths.matrix_norm_tsv << '\t'
                                << sim_paths.matrix_abs_tsv << '\t'
                                << sim_paths.stats_tsv << '\t'
                                << sim_paths.cluster_signatures_dir << '\n';
            summary.similarity_reports_written += 1;
        }

        const bool need_variant_logic = options.write_region_vcf || options.write_debug_reports;
        if (need_variant_logic) {
            VariantBubbleReport report = ::panvar::write_variant_reports_for_bubble(
                bubble,
                unique_alleles,
                clusters,
                cluster_of_unique,
                assignments,
                options.reference_path,
                reference_sequence,
                options.min_sv_bp,
                options.write_debug_reports,
                options.debug_out_dir,
                dotplot_gene_ptr,
                options.dotplot_gene_matches,
                options.minimap_preset,
                options.minimap_best_n,
                options.minimap_emit_secondary,
                options.split_ins_use_geometric_svlen,
                options.classify_ins);

            summary.debug_reports_written += report.debug_reports_written;
            summary.dotplots_written += report.dotplots_written;

            if (options.write_region_vcf &&
                report.has_reference_assignment && !report.vcf_rows.empty()) {
                RegionVcfBubble item;
                item.bubble_id = bubble.id;
                item.source = bubble.source;
                item.sink = bubble.sink;
                item.reference_cluster_id = report.reference_cluster_id;
                item.reference_interval_start = report.reference_interval_start;
                item.reference_interval_end = report.reference_interval_end;
                item.cluster_by_path = std::move(report.cluster_by_path);
                item.allele_length_by_path = std::move(report.allele_length_by_path);
                item.variant_rows = std::move(report.vcf_rows);
                region_vcf_bubbles.push_back(std::move(item));
            }
        }

        if (options.write_odgi_viz_inputs) {
            std::unordered_map<std::string, OdgiPathEntry> best_by_path;
            best_by_path.reserve(assignments.size() * 2);
            for (const auto& assignment : assignments) {
                OdgiPathEntry entry;
                entry.path_name = assignment.path_name;
                entry.cluster_id = cluster_of_unique[assignment.unique_idx];
                entry.interval_start = assignment.interval_start;
                entry.interval_end = assignment.interval_end;
                entry.source_to_sink = assignment.source_to_sink;

                const auto found = best_by_path.find(entry.path_name);
                if (found == best_by_path.end() || (!found->second.source_to_sink && entry.source_to_sink)) {
                    best_by_path[entry.path_name] = std::move(entry);
                }
            }

            std::string paths_file = ".";
            std::string colors_file = ".";
            std::string viz_script = ".";
            std::string png_file = ".";
            int odgi_rc = -1;
            std::size_t range_start_bp = 0;
            std::size_t range_end_bp = 0;
            std::string status = "skipped:no-reference-assignment";
            std::size_t rows_exported = 0;
            std::size_t cluster_count = 0;
            const std::string bubble_dir = options.odgi_viz_out_dir + "/bubble_" + std::to_string(bubble.id);

            const auto ref_row_it = best_by_path.find(options.reference_path);
            if (ref_row_it != best_by_path.end()) {
                auto ref_window = compute_reference_window(
                    *reference_path_record,
                    ref_row_it->second.interval_start,
                    ref_row_it->second.interval_end,
                    graph.nodes,
                    options.flank_nodes);

                if (!ref_window.has_value()) {
                    status = "skipped:invalid-reference-window";
                } else {
                    range_start_bp = ref_window->range_start_bp;
                    range_end_bp = ref_window->range_end_bp;
                    status = "ok";

                    std::vector<OdgiPathEntry> ordered;
                    ordered.reserve(best_by_path.size());
                    for (auto& [path_name, entry] : best_by_path) {
                        (void)path_name;
                        ordered.push_back(std::move(entry));
                    }
                    std::sort(ordered.begin(), ordered.end(), [](const OdgiPathEntry& lhs, const OdgiPathEntry& rhs) {
                        if (lhs.cluster_id != rhs.cluster_id) {
                            return lhs.cluster_id < rhs.cluster_id;
                        }
                        return lhs.path_name < rhs.path_name;
                    });

                    if (options.max_paths_per_bubble > 0 && ordered.size() > options.max_paths_per_bubble) {
                        bool ref_in_prefix = false;
                        std::size_t ref_index = ordered.size();
                        for (std::size_t i = 0; i < ordered.size(); ++i) {
                            if (ordered[i].path_name == options.reference_path) {
                                ref_index = i;
                                ref_in_prefix = i < options.max_paths_per_bubble;
                                break;
                            }
                        }

                        ordered.resize(options.max_paths_per_bubble);
                        if (!ref_in_prefix && ref_index < best_by_path.size() && !ordered.empty()) {
                            ordered.back() = ref_row_it->second;
                        }
                    }

                    rows_exported = ordered.size();
                    std::set<std::size_t> cluster_ids;
                    for (const auto& row : ordered) {
                        cluster_ids.insert(row.cluster_id);
                    }
                    cluster_count = cluster_ids.size();
                    const auto color_by_cluster = cluster_palette(cluster_ids);

                    std::filesystem::create_directories(bubble_dir);
                    paths_file = bubble_dir + "/paths.txt";
                    colors_file = bubble_dir + "/path_colors.tsv";
                    viz_script = bubble_dir + "/viz.sh";
                    png_file = bubble_dir + "/plot.png";

                    std::ofstream paths_out(paths_file);
                    std::ofstream colors_out(colors_file);
                    if (!paths_out || !colors_out) {
                        throw std::runtime_error(
                            "Failed to write ODGI input files for bubble " + std::to_string(bubble.id));
                    }
                    for (const auto& row : ordered) {
                        paths_out << row.path_name << '\n';
                        const auto color_it = color_by_cluster.find(row.cluster_id);
                        const std::string color = (color_it == color_by_cluster.end()) ? "#444444" : color_it->second;
                        colors_out << row.path_name << '\t' << color << '\n';
                    }
                    paths_out.close();
                    colors_out.close();

                    const std::string range =
                        options.reference_path + ":" + std::to_string(range_start_bp) + "-" +
                        std::to_string(range_end_bp);
                    const std::size_t image_height = std::max<std::size_t>(
                        400,
                        120 + rows_exported * (options.odgi_path_height + 2));
                    const std::string command = options.odgi_bin + " viz"
                                                " -i " + shell_quote(odgi_input) +
                                                " -o " + shell_quote(png_file) +
                                                " -r " + shell_quote(range) +
                                                " -p " + shell_quote(paths_file) +
                                                " -F " + shell_quote(colors_file) +
                                                " -x " + std::to_string(options.odgi_width) +
                                                " -a " + std::to_string(options.odgi_path_height) +
                                                " -y " + std::to_string(image_height);

                    std::ofstream script_out(viz_script);
                    if (!script_out) {
                        throw std::runtime_error("Failed to write ODGI script: " + viz_script);
                    }
                    script_out << "#!/usr/bin/env bash\nset -euo pipefail\n" << command << '\n';
                    script_out.close();
                    std::filesystem::permissions(
                        viz_script,
                        std::filesystem::perms::owner_exec | std::filesystem::perms::group_exec |
                            std::filesystem::perms::others_exec,
                        std::filesystem::perm_options::add);

                    summary.bubbles_with_odgi_exports += 1;

                    if (options.run_odgi) {
                        odgi_rc = std::system(command.c_str());
                        if (odgi_rc == 0) {
                            summary.odgi_runs_ok += 1;
                        } else {
                            summary.odgi_runs_failed += 1;
                            status = "run-failed";
                        }
                    }
                }
            }

            odgi_manifest << bubble.id << '\t'
                          << bubble.source << '\t'
                          << bubble.sink << '\t'
                          << bubble_dir << '\t'
                          << rows_exported << '\t'
                          << cluster_count << '\t'
                          << range_start_bp << '\t'
                          << range_end_bp << '\t'
                          << paths_file << '\t'
                          << colors_file << '\t'
                          << viz_script << '\t'
                          << png_file << '\t'
                          << odgi_rc << '\t'
                          << status << '\n';
        }

        if (options.show_progress) {
            std::cerr
                << "[allele] bubble " << bubble.id
                << " done: unique=" << unique_alleles.size()
                << ", clusters=" << clusters.size()
                << ", paths=" << assignments.size()
                << ", elapsed=" << std::fixed << std::setprecision(2)
                << elapsed_seconds(bubble_start) << "s"
                << ", total=" << elapsed_seconds(run_start) << "s\n";
        }
    }

    if (options.write_region_vcf) {
        if (options.region_vcf_path.empty()) {
            throw std::runtime_error("--vcf-out cannot be empty when region-level VCF export is enabled");
        }
        if (reference_prefix_bp.empty()) {
            throw std::runtime_error("Reference path prefix lengths are unavailable for region-level VCF export");
        }
        const std::size_t region_start_1based =
            reference_meta.has_interval ? reference_meta.region_start_1based : 1;
        ::panvar::write_region_level_vcf(
            options.region_vcf_path,
            options.reference_path,
            reference_meta,
            reference_sequence,
            reference_prefix_bp,
            region_vcf_bubbles,
            vcf_sample_names,
            region_start_1based,
            options.vcf_merge_window_bp,
            options.vcf_merge_mode,
            options.vcf_merge_lenient_window_bp,
            options.vcf_merge_lenient_min_ref_jaccard,
            options.vcf_merge_min_seq_similarity,
            options.vcf_merge_max_seq_edit_fraction,
            &summary.region_vcf_records);
    }

    if (summary_out != nullptr) {
        *summary_out = summary;
    }

    if (options.show_progress) {
        if (use_predefined_clusters && predefined_missing_assignments_total > 0) {
            std::cerr
                << "[allele] warning: " << predefined_missing_assignments_total
                << " path assignments (across " << predefined_missing_paths_total
                << " unique paths) were not found in --clusters-json and were placed in synthetic singleton labels\n";
        }
        std::cerr
            << "[allele] completed " << bubbles.size() << " bubbles in "
            << std::fixed << std::setprecision(2) << elapsed_seconds(run_start) << "s\n";
    }
}

} // namespace panvar
