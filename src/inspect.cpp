#include "panvar/inspect.hpp"

#include "panvar/bubble_path.hpp"
#include "panvar/bubbles.hpp"
#include "panvar/cli_utils.hpp"
#include "panvar/gfa.hpp"
#include "panvar/graph_utils.hpp"
#include "panvar/output.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <queue>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <zlib.h>

namespace panvar {
namespace {

void print_inspect_help() {
    std::cout
        << "Usage:\n"
        << "  panvar inspect -i <graph.gfa> (-b <prefix> | -c <bubbles.csv>) [--bubble-id <N>] [options]\n\n"
        << "Options:\n"
        << "  -i, --gfa <path>                 Input GFA file (required)\n"
        << "  -c, --bubbles-csv <path>         Module-1 bubbles CSV (required if no prefix)\n"
        << "  -b, --bubble-prefix-in <prefix>  Module-1 output prefix from 'panvar bubble'\n"
        << "                                   (auto-uses <prefix>.bubbles.csv)\n"
        << "      --bubble-id <N>              Bubble ID to inspect (default: inspect all bubbles)\n"
        << "  -o, --out-prefix <prefix>        Output prefix (default: inspect)\n"
        << "      --fasta-out <path>           Explicit FASTA.GZ path (requires --bubble-id)\n"
        << "      --table-out <path>           Explicit node-count TSV path (requires --bubble-id)\n"
        << "      --edge-table-out <path>      Explicit edge-count TSV path (requires --bubble-id)\n"
        << "      --cluster                    Group paths by source->sink walk; write a\n"
        << "                                   <prefix>.bubble_<N>.clusters.tsv per bubble\n"
        << "      --cluster-similarity <f>     Walk similarity threshold for --cluster (default: 0.90)\n"
        << "  -q, --quiet                      Disable the progress bar\n"
        << "  -h, --help                       Show this help\n";
}

struct InspectNodeCount {
    std::size_t total = 0;
    std::size_t forward = 0;
    std::size_t reverse = 0;
};

std::string tsv_sanitize(std::string value) {
    for (char& c : value) {
        if (c == '\t' || c == '\n' || c == '\r') {
            c = ' ';
        }
    }
    return value;
}

void gz_write_string(gzFile file, const std::string& text, const std::string& output_path) {
    if (text.empty()) {
        return;
    }
    const int written = gzwrite(file, text.data(), static_cast<unsigned int>(text.size()));
    if (written != static_cast<int>(text.size())) {
        throw std::runtime_error("Failed to write compressed FASTA: " + output_path);
    }
}

void gz_write_fasta_record(
    gzFile file,
    const std::string& output_path,
    const std::string& name,
    const std::string& sequence) {

    gz_write_string(file, ">" + name + "\n", output_path);
    constexpr std::size_t kLineWidth = 80;
    for (std::size_t i = 0; i < sequence.size(); i += kLineWidth) {
        gz_write_string(file, sequence.substr(i, kLineWidth) + "\n", output_path);
    }
}

struct InspectBubbleResult {
    std::size_t paths_written = 0;
    std::size_t clusters_written = 0;
    std::string fasta_out_path;
    std::string table_out_path;
    std::string edge_table_out_path;
    std::string node_lengths_out_path;
    std::string clusters_out_path;
};

// Optional sidecar outputs requested for a bubble.
struct InspectEmit {
    std::string node_lengths_out_path;  // always written
    std::string clusters_out_path;      // written only when cluster is true
    bool cluster = false;
    double cluster_similarity = 0.90;
};

// Orientation-aware edge key for a step pair, e.g. "12+>13-". Adjacency-aware so a
// tandem self-loop (the same edge traversed repeatedly) shows up as a high count.
std::string edge_key(const PathStep& a, const PathStep& b) {
    return a.node_id + (a.reverse ? '-' : '+') + ">" + b.node_id + (b.reverse ? '-' : '+');
}

// One path's traversal of a bubble: node counts and edge counts.
struct InspectPathRow {
    std::string name;
    std::size_t sequence_length = 0;
    std::unordered_map<std::string, InspectNodeCount> node_counts;
    std::unordered_map<std::string, std::size_t> edge_counts;
};

// --- Path clustering by source->sink walk -----------------------------------
//
// A walk is summarized as an oriented, bp-weighted token multiset: each step
// contributes its node's bp length to the token "<node_id><strand>". Two walks
// are compared with weighted Jaccard = sum(min) / sum(max) over the token union,
// which captures inversions (strand in the token) and copy number (repeats add
// weight) while staying order-insensitive. std::map keeps tokens sorted so the
// Jaccard merge is a linear two-pointer pass.
using TokenWeights = std::map<std::string, std::size_t>;

TokenWeights build_token_weights(const Graph& graph, const std::vector<PathStep>& steps) {
    TokenWeights weights;
    for (const auto& step : steps) {
        const auto node_it = graph.nodes.find(step.node_id);
        const std::size_t len =
            node_it == graph.nodes.end() ? 1 : std::max<std::size_t>(1, node_it->second.sequence.size());
        std::string token = step.node_id;
        token.push_back(step.reverse ? '-' : '+');
        weights[token] += len;
    }
    return weights;
}

double weighted_jaccard_tokens(const TokenWeights& a, const TokenWeights& b) {
    std::size_t inter = 0;
    std::size_t uni = 0;
    auto ia = a.begin();
    auto ib = b.begin();
    while (ia != a.end() && ib != b.end()) {
        if (ia->first == ib->first) {
            inter += std::min(ia->second, ib->second);
            uni += std::max(ia->second, ib->second);
            ++ia;
            ++ib;
        } else if (ia->first < ib->first) {
            uni += ia->second;
            ++ia;
        } else {
            uni += ib->second;
            ++ib;
        }
    }
    for (; ia != a.end(); ++ia) {
        uni += ia->second;
    }
    for (; ib != b.end(); ++ib) {
        uni += ib->second;
    }
    return uni == 0 ? 0.0 : static_cast<double>(inter) / static_cast<double>(uni);
}

// --- Fast non-greedy clustering: walk MinHash sketch + threshold-graph CC ------
//
// A fast, order-independent clusterer. Each walk is summarized by a bottom-k MinHash
// sketch over oriented node-step shingles; sketch Jaccard estimates identity; walks
// within `--cluster-similarity` are united with a disjoint-set forest, so clusters are
// the connected components (transitive, order-independent).
//
// The sketch is **multiplicity-aware**: a shingle seen k times contributes k distinct
// sketch elements (its occurrence index is folded into the hash), so the sketch Jaccard
// approximates the shingle *multiset* Jaccard rather than the plain set Jaccard. This
// matters for tandem repeats (e.g. LPA KIV-2): two haplotypes that share the same repeat
// unit but at different copy numbers have nearly identical shingle *sets* but very
// different *multisets*, so a set sketch would merge all copy numbers into one cluster
// while the multiset sketch separates them by copy number — consistent with the
// copy-number-aware `weighted_jaccard_tokens` fallback used for un-shingleable walks.
std::uint64_t splitmix64(std::uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

// hash_step_token (oriented per-step 64-bit hash) comes from graph_utils.hpp.

constexpr std::size_t kWalkSketchShingle = 3;
constexpr std::size_t kWalkSketchSize = 512;

std::vector<std::uint64_t> build_walk_minhash_sketch(
    const std::vector<std::uint64_t>& tokens,
    std::size_t shingle_size,
    std::size_t sketch_size) {

    if (shingle_size == 0 || sketch_size == 0 || tokens.size() < shingle_size) {
        return {};
    }
    std::priority_queue<std::uint64_t> top_hashes;  // max-heap -> keep the bottom k
    // Per-shingle occurrence counter: the n-th occurrence of a shingle is salted with n,
    // so repeats yield distinct sketch elements (multiset rather than set semantics).
    std::unordered_map<std::uint64_t, std::uint32_t> occ;
    for (std::size_t i = 0; i + shingle_size <= tokens.size(); ++i) {
        std::uint64_t h = 1469598103934665603ULL;
        for (std::size_t k = 0; k < shingle_size; ++k) {
            const std::uint64_t v = splitmix64(tokens[i + k] + (0x9e3779b97f4a7c15ULL * (k + 1)));
            h ^= v;
            h *= 1099511628211ULL;
        }
        h = splitmix64(h);
        const std::uint32_t n = occ[h]++;  // 0 for the first occurrence, 1 for the next, ...
        const std::uint64_t e = splitmix64(h + 0x9e3779b97f4a7c15ULL * (n + 1));
        if (top_hashes.size() < sketch_size) {
            top_hashes.push(e);
        } else if (e < top_hashes.top()) {
            top_hashes.pop();
            top_hashes.push(e);
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

double sketch_jaccard(const std::vector<std::uint64_t>& a, const std::vector<std::uint64_t>& b) {
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
    return uni == 0 ? 0.0 : static_cast<double>(inter) / static_cast<double>(uni);
}

// identity ~= 2J / (1 + J) from k-mer/shingle Jaccard.
double estimate_identity_from_jaccard(double jaccard) {
    if (jaccard <= 0.0) {
        return 0.0;
    }
    return std::clamp((2.0 * jaccard) / (1.0 + jaccard), 0.0, 1.0);
}

// DisjointSet (union-find) is shared from graph_utils.hpp.

// One distinct canonical walk and the paths realizing it.
struct UniqueWalk {
    std::string signature;
    TokenWeights weights;
    std::vector<std::uint64_t> sketch;  // walk MinHash sketch (for CC clustering)
    std::vector<std::string> members;   // path names, in encounter order
};

struct ClusterOut {
    std::size_t cluster_id = 0;
    std::size_t n_paths = 0;
    std::string representative;
    std::vector<std::string> members;
};

// Connected-components clustering: edge iff sketch-estimated identity >=
// threshold; clusters = components (transitive). Representative = member minimizing
// max-then-mean intra-cluster distance, tie-broken by higher support then signature.
std::vector<ClusterOut> cluster_unique_walks_cc(
    const std::vector<UniqueWalk>& uniques, double threshold) {

    const std::size_t n = uniques.size();
    if (n == 0) {
        return {};
    }
    std::vector<std::vector<double>> dist(n, std::vector<double>(n, 0.0));
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i + 1; j < n; ++j) {
            // Walks too short to shingle have empty sketches; fall back to the exact
            // bp-weighted Jaccard so short bubbles still cluster instead of all-singleton.
            const double id = (uniques[i].sketch.empty() || uniques[j].sketch.empty())
                ? weighted_jaccard_tokens(uniques[i].weights, uniques[j].weights)
                : estimate_identity_from_jaccard(sketch_jaccard(uniques[i].sketch, uniques[j].sketch));
            dist[i][j] = dist[j][i] = 1.0 - id;
        }
    }

    const double max_norm_dist = std::clamp(1.0 - threshold, 0.0, 1.0);
    const double eps = 1e-12;
    DisjointSet dsu(n);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i + 1; j < n; ++j) {
            if (dist[i][j] <= max_norm_dist + eps) {
                dsu.unite(i, j);
            }
        }
    }

    std::unordered_map<std::size_t, std::vector<std::size_t>> by_root;
    for (std::size_t i = 0; i < n; ++i) {
        by_root[dsu.find(i)].push_back(i);
    }

    std::vector<ClusterOut> out;
    out.reserve(by_root.size());
    for (auto& kv : by_root) {
        const std::vector<std::size_t>& members = kv.second;
        std::size_t rep = members.front();
        double best_max = std::numeric_limits<double>::infinity();
        double best_mean = std::numeric_limits<double>::infinity();
        std::size_t best_support = 0;
        for (const std::size_t a : members) {
            double sum = 0.0;
            double mx = 0.0;
            std::size_t cnt = 0;
            for (const std::size_t b : members) {
                if (a == b) {
                    continue;
                }
                sum += dist[a][b];
                mx = std::max(mx, dist[a][b]);
                ++cnt;
            }
            const double mean = cnt == 0 ? 0.0 : sum / static_cast<double>(cnt);
            const std::size_t support = uniques[a].members.size();
            const bool better =
                mx + eps < best_max ||
                (std::abs(mx - best_max) <= eps && mean + eps < best_mean) ||
                (std::abs(mx - best_max) <= eps && std::abs(mean - best_mean) <= eps &&
                 (support > best_support ||
                  (support == best_support && uniques[a].signature < uniques[rep].signature)));
            if (better) {
                rep = a;
                best_max = mx;
                best_mean = mean;
                best_support = support;
            }
        }
        ClusterOut co;
        for (const std::size_t m : members) {
            co.n_paths += uniques[m].members.size();
            for (const auto& name : uniques[m].members) {
                co.members.push_back(name);
            }
        }
        co.representative = uniques[rep].members.front();
        std::sort(co.members.begin(), co.members.end());
        out.push_back(std::move(co));
    }
    std::sort(out.begin(), out.end(), [](const ClusterOut& a, const ClusterOut& b) {
        if (a.n_paths != b.n_paths) {
            return a.n_paths > b.n_paths;
        }
        return a.representative < b.representative;
    });
    for (std::size_t i = 0; i < out.size(); ++i) {
        out[i].cluster_id = i;
    }
    return out;
}

InspectBubbleResult write_inspect_outputs_for_bubble(
    const Graph& graph,
    const Bubble& bubble,
    const std::string& fasta_out_path,
    const std::string& table_out_path,
    const std::string& edge_table_out_path,
    const InspectEmit& emit) {

    cli::ensure_parent_dir_for_file(fasta_out_path);
    cli::ensure_parent_dir_for_file(table_out_path);
    cli::ensure_parent_dir_for_file(edge_table_out_path);

    gzFile fasta = gzopen(fasta_out_path.c_str(), "wb");
    if (fasta == nullptr) {
        throw std::runtime_error("Failed to write compressed FASTA: " + fasta_out_path);
    }

    std::unordered_set<std::string> inside_nodes;
    inside_nodes.reserve(bubble.inside.size() * 2);
    inside_nodes.insert(bubble.inside.begin(), bubble.inside.end());

    // Pass 1: collect per-path node + edge counts and stream the FASTA. Edge columns
    // are not known up front, so the matrices are written after the union is built.
    std::vector<InspectPathRow> rows;
    std::vector<std::string> edge_order;          // edge columns, first-seen order
    std::unordered_set<std::string> edge_seen;
    // Distinct canonical walks, collected only when clustering is requested.
    std::unordered_map<std::string, std::size_t> sig_to_unique;
    std::vector<UniqueWalk> uniques;
    try {
        for (const auto& path : graph.paths) {
            const BubblePathIndex index = build_bubble_path_index(path);
            const auto interval = find_best_bubble_path_interval(index, bubble);
            if (!interval.has_value()) {
                continue;
            }

            const std::vector<PathStep> steps = canonical_bubble_path_steps(path, bubble, *interval);
            if (steps.empty()) {
                continue;
            }

            const std::string sequence = spell_path_steps_sequence(graph, steps);
            InspectPathRow row;
            row.name = tsv_sanitize(path.name);
            row.sequence_length = sequence.size();
            row.node_counts.reserve(bubble.inside.size());
            for (const auto& step : steps) {
                if (inside_nodes.find(step.node_id) == inside_nodes.end()) {
                    continue;
                }
                auto& c = row.node_counts[step.node_id];
                c.total += 1;
                if (step.reverse) {
                    c.reverse += 1;
                } else {
                    c.forward += 1;
                }
            }
            for (std::size_t i = 1; i < steps.size(); ++i) {
                const std::string key = edge_key(steps[i - 1], steps[i]);
                if (edge_seen.insert(key).second) {
                    edge_order.push_back(key);
                }
                row.edge_counts[key] += 1;
            }

            const std::string fasta_name =
                row.name +
                " bubble=" + std::to_string(bubble.id) +
                " source=" + bubble.source +
                " sink=" + bubble.sink +
                " length_bp=" + std::to_string(row.sequence_length) +
                " source_to_sink=" + (interval->source_to_sink ? "1" : "0") +
                " interval=" + std::to_string(interval->left) + "-" + std::to_string(interval->right);
            gz_write_fasta_record(fasta, fasta_out_path, fasta_name, sequence);

            if (emit.cluster) {
                const std::string sig = build_walk_signature(steps);
                const auto it = sig_to_unique.find(sig);
                if (it == sig_to_unique.end()) {
                    sig_to_unique.emplace(sig, uniques.size());
                    UniqueWalk uw;
                    uw.signature = sig;
                    uw.weights = build_token_weights(graph, steps);
                    std::vector<std::uint64_t> tokens;
                    tokens.reserve(steps.size());
                    for (const auto& step : steps) {
                        tokens.push_back(hash_step_token(step));
                    }
                    uw.sketch = build_walk_minhash_sketch(tokens, kWalkSketchShingle, kWalkSketchSize);
                    uw.members.push_back(row.name);
                    uniques.push_back(std::move(uw));
                } else {
                    uniques[it->second].members.push_back(row.name);
                }
            }

            rows.push_back(std::move(row));
        }
    } catch (...) {
        gzclose(fasta);
        throw;
    }
    if (gzclose(fasta) != Z_OK) {
        throw std::runtime_error("Failed to close compressed FASTA: " + fasta_out_path);
    }

    // Stable, deterministic edge column order.
    std::sort(edge_order.begin(), edge_order.end());

    std::ofstream table_out(table_out_path);
    if (!table_out) {
        throw std::runtime_error("Failed to write inspect table: " + table_out_path);
    }
    std::ofstream edge_out(edge_table_out_path);
    if (!edge_out) {
        throw std::runtime_error("Failed to write inspect edge table: " + edge_table_out_path);
    }

    table_out << "path_name\tpath_length_bp";
    for (const auto& node_id : bubble.inside) {
        table_out << "\tnode." << tsv_sanitize(node_id);
    }
    table_out << "\n";
    edge_out << "path_name\tpath_length_bp";
    for (const auto& key : edge_order) {
        edge_out << "\tedge." << key;
    }
    edge_out << "\n";

    InspectBubbleResult result;
    result.fasta_out_path = fasta_out_path;
    result.table_out_path = table_out_path;
    result.edge_table_out_path = edge_table_out_path;
    result.node_lengths_out_path = emit.node_lengths_out_path;
    for (const InspectPathRow& row : rows) {
        table_out << row.name << '\t' << row.sequence_length;
        for (const auto& node_id : bubble.inside) {
            const auto it = row.node_counts.find(node_id);
            if (it == row.node_counts.end()) {
                table_out << "\t0:0:0";
            } else {
                table_out << '\t' << it->second.total << ':' << it->second.forward << ':'
                          << it->second.reverse;
            }
        }
        table_out << '\n';

        edge_out << row.name << '\t' << row.sequence_length;
        for (const auto& key : edge_order) {
            const auto it = row.edge_counts.find(key);
            edge_out << '\t' << (it == row.edge_counts.end() ? 0 : it->second);
        }
        edge_out << '\n';
        ++result.paths_written;
    }

    // Per-node bp lengths, in the same order as the node.* columns above, so the
    // node coverage heatmap can length-scale its x-axis without reordering.
    cli::ensure_parent_dir_for_file(emit.node_lengths_out_path);
    std::ofstream node_lengths_out(emit.node_lengths_out_path);
    if (!node_lengths_out) {
        throw std::runtime_error("Failed to write node lengths TSV: " + emit.node_lengths_out_path);
    }
    node_lengths_out << "node_id\tlength_bp\n";
    for (const auto& node_id : bubble.inside) {
        const auto it = graph.nodes.find(node_id);
        const std::size_t len = it == graph.nodes.end() ? 0 : it->second.sequence.size();
        node_lengths_out << tsv_sanitize(node_id) << '\t' << len << '\n';
    }

    if (emit.cluster) {
        const std::vector<ClusterOut> clusters =
            cluster_unique_walks_cc(uniques, emit.cluster_similarity);
        cli::ensure_parent_dir_for_file(emit.clusters_out_path);
        std::ofstream clusters_out(emit.clusters_out_path);
        if (!clusters_out) {
            throw std::runtime_error("Failed to write clusters TSV: " + emit.clusters_out_path);
        }
        clusters_out << "cluster_id\tn_paths\trepresentative_path\tmembers\n";
        for (const auto& c : clusters) {
            clusters_out << c.cluster_id << '\t' << c.n_paths << '\t' << c.representative << '\t';
            for (std::size_t i = 0; i < c.members.size(); ++i) {
                if (i > 0) {
                    clusters_out << ';';
                }
                clusters_out << c.members[i];
            }
            clusters_out << '\n';
        }
        result.clusters_written = clusters.size();
        result.clusters_out_path = emit.clusters_out_path;
    }

    return result;
}

} // namespace

int run_inspect_command(const std::vector<std::string>& args) {
    if (args.empty()) {
        print_inspect_help();
        return 0;
    }

    std::string gfa_path;
    std::string bubbles_csv_path;
    std::string bubble_prefix_in;
    std::string out_prefix = "inspect";
    std::string fasta_out_path;
    std::string table_out_path;
    std::string edge_table_out_path;
    std::size_t bubble_id = 0;
    bool cluster = false;
    double cluster_similarity = 0.90;
    bool quiet = false;

    for (std::size_t i = 0; i < args.size(); ++i) {
        const std::string& arg = args[i];
        auto require_value = [&](const std::string& flag) -> const std::string& {
            if (i + 1 >= args.size()) {
                throw std::runtime_error("Missing value after " + flag);
            }
            return args[++i];
        };

        if (arg == "-h" || arg == "--help") {
            print_inspect_help();
            return 0;
        }
        if (arg == "-i" || arg == "--gfa") {
            gfa_path = require_value(arg);
            continue;
        }
        if (arg == "-c" || arg == "--bubbles-csv") {
            bubbles_csv_path = require_value(arg);
            continue;
        }
        if (arg == "-b" || arg == "--bubble-prefix-in") {
            bubble_prefix_in = require_value(arg);
            continue;
        }
        if (arg == "--bubble-id") {
            bubble_id = cli::parse_size_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "-o" || arg == "--out-prefix") {
            out_prefix = require_value(arg);
            continue;
        }
        if (arg == "--fasta-out") {
            fasta_out_path = require_value(arg);
            continue;
        }
        if (arg == "--table-out") {
            table_out_path = require_value(arg);
            continue;
        }
        if (arg == "--edge-table-out") {
            edge_table_out_path = require_value(arg);
            continue;
        }
        if (arg == "--cluster") {
            cluster = true;
            continue;
        }
        if (arg == "--cluster-similarity") {
            cluster_similarity = cli::parse_similarity_arg(arg, require_value(arg));
            continue;
        }
        if (arg == "-q" || arg == "--quiet") {
            quiet = true;
            continue;
        }
        throw std::runtime_error("Unknown option: " + arg);
    }

    if (gfa_path.empty()) {
        throw std::runtime_error("Missing required input: --gfa <path>");
    }
    if (!bubble_prefix_in.empty()) {
        const std::string derived = bubble_prefix_in + ".bubbles.csv";
        if (bubbles_csv_path.empty()) {
            bubbles_csv_path = derived;
        } else if (bubbles_csv_path != derived) {
            throw std::runtime_error(
                "Conflicting bubble inputs: --bubble-prefix-in resolves to '" +
                derived + "' but --bubbles-csv is '" + bubbles_csv_path + "'");
        }
    }
    if (bubbles_csv_path.empty()) {
        throw std::runtime_error("Missing required input: --bubbles-csv <path> or --bubble-prefix-in <prefix>");
    }
    if (bubble_id == 0 && (!fasta_out_path.empty() || !table_out_path.empty() || !edge_table_out_path.empty())) {
        throw std::runtime_error("--fasta-out/--table-out/--edge-table-out require --bubble-id; use --out-prefix when inspecting all bubbles");
    }

    ParseGfaOptions parse_options;
    parse_options.include_paths = true;
    parse_options.include_sequences = true;
    const Graph graph = parse_gfa(gfa_path, parse_options);
    if (graph.paths.empty()) {
        throw std::runtime_error("Input GFA has no P/W paths; inspect requires paths");
    }

    const std::vector<Bubble> bubbles = read_bubbles_csv(bubbles_csv_path);
    if (bubbles.empty()) {
        throw std::runtime_error("Bubbles CSV has no bubbles: " + bubbles_csv_path);
    }

    std::vector<const Bubble*> selected_bubbles;
    if (bubble_id == 0) {
        selected_bubbles.reserve(bubbles.size());
        for (const auto& bubble : bubbles) {
            selected_bubbles.push_back(&bubble);
        }
    } else {
        const auto bubble_it = std::find_if(
            bubbles.begin(),
            bubbles.end(),
            [&](const Bubble& b) { return b.id == bubble_id; });
        if (bubble_it == bubbles.end()) {
            throw std::runtime_error("Bubble ID not found in bubbles CSV: " + std::to_string(bubble_id));
        }
        selected_bubbles.push_back(&(*bubble_it));
    }

    // Verbose per-bubble block is useful for a single bubble; for an all-bubbles run
    // it just floods stdout, so we show a progress bar on stderr instead.
    const bool single_bubble = bubble_id != 0;

    std::size_t total_paths_written = 0;
    std::cout
        << "Input graph: " << gfa_path << "\n"
        << "Bubble source: " << bubbles_csv_path << "\n"
        << "Bubbles inspected: " << selected_bubbles.size() << "\n";

    cli::ProgressBar progress((single_bubble || quiet) ? "" : "Inspecting bubbles",
                              single_bubble ? 0 : selected_bubbles.size());

    for (const Bubble* bubble_ptr : selected_bubbles) {
        const Bubble& bubble = *bubble_ptr;
        std::string bubble_fasta_out_path = fasta_out_path;
        std::string bubble_table_out_path = table_out_path;
        std::string bubble_edge_table_out_path = edge_table_out_path;
        if (bubble_fasta_out_path.empty()) {
            bubble_fasta_out_path = out_prefix + ".bubble_" + std::to_string(bubble.id) + ".paths.fa.gz";
        }
        if (bubble_table_out_path.empty()) {
            bubble_table_out_path = out_prefix + ".bubble_" + std::to_string(bubble.id) + ".node_counts.tsv";
        }
        if (bubble_edge_table_out_path.empty()) {
            bubble_edge_table_out_path = out_prefix + ".bubble_" + std::to_string(bubble.id) + ".edge_counts.tsv";
        }

        InspectEmit emit;
        emit.node_lengths_out_path =
            out_prefix + ".bubble_" + std::to_string(bubble.id) + ".node_lengths.tsv";
        emit.clusters_out_path =
            out_prefix + ".bubble_" + std::to_string(bubble.id) + ".clusters.tsv";
        emit.cluster = cluster;
        emit.cluster_similarity = cluster_similarity;

        const InspectBubbleResult result = write_inspect_outputs_for_bubble(
            graph,
            bubble,
            bubble_fasta_out_path,
            bubble_table_out_path,
            bubble_edge_table_out_path,
            emit);
        total_paths_written += result.paths_written;
        progress.tick();

        if (single_bubble) {
            std::cout
                << "Bubble ID: " << bubble.id << "\n"
                << "Source/sink: " << bubble.source << "/" << bubble.sink << "\n"
                << "Inside nodes: " << bubble.inside.size() << "\n"
                << "Paths written: " << result.paths_written << "\n"
                << "Wrote: " << result.fasta_out_path << "\n"
                << "Wrote: " << result.table_out_path << "\n"
                << "Wrote: " << result.edge_table_out_path << "\n"
                << "Wrote: " << result.node_lengths_out_path << "\n";
            if (emit.cluster) {
                std::cout << "Clusters: " << result.clusters_written << "\n"
                          << "Wrote: " << result.clusters_out_path << "\n";
            }
        }
    }
    progress.done();

    std::cout << "Total paths written: " << total_paths_written << "\n";

    return 0;
}

} // namespace panvar
