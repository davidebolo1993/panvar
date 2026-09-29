#include "panvar/inspect.hpp"

#include "panvar/bubble_path.hpp"
#include "panvar/bubbles.hpp"
#include "panvar/cli_utils.hpp"
#include "panvar/cosigt.hpp"
#include "panvar/gfa.hpp"
#include "panvar/graph_utils.hpp"
#include "panvar/output.hpp"
#include "panvar/ref_path.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <filesystem>
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
        << "      --genotype-table             Write <prefix>.bubble_<N>.genotypes.tsv.gz (one row per\n"
        << "                                   haplotype, one column per allele, 0/1) and the matching\n"
        << "                                   <prefix>.bubble_<N>.alleles.fa.gz; identical sequences share\n"
        << "                                   an allele. A haplotype not crossing the bubble reads NA\n"
        << "      --samples <tsv>              cosigt sample->haplotype table: genotype rows become samples\n"
        << "                                   (allele copy counts, NA if any assigned haplotype does not\n"
        << "                                   cross the bubble). Implies --genotype-table\n"
        << "  -r, --reference-path <name>      Orient --genotype-table alleles to this path's strand\n"
        << "                                   (default: bubble source->sink)\n"
        << "  -q, --quiet                     Disable the progress bar\n"
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
        throw std::runtime_error("Failed to write compressed output: " + output_path);
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
    std::size_t alleles_written = 0;
    bool reference_seen = false;
};

// Optional sidecar outputs requested for a bubble.
struct InspectEmit {
    std::string node_lengths_out_path;  // always written
    std::string walks_out_path;         // always written: the paths FASTA's records as node walks
    std::string clusters_out_path;      // written only when cluster is true
    bool cluster = false;
    double cluster_similarity = 0.90;
    // --genotype-table
    bool genotype_table = false;
    std::string genotypes_out_path;
    std::string alleles_out_path;
    std::string reference_path;  // resolved path name; empty keeps source->sink
    std::string allele_walks_out_path;  // --genotype-table: the allele FASTA's walks
    std::string nodes_json_out_path;    // --genotype-table: the bubble's nodes, sequences and positions
    const std::vector<std::string>* samples = nullptr;  // cosigt row order; null = haplotype rows
    const std::unordered_map<std::string, std::vector<std::string>>* sample_haplotypes = nullptr;
};

// Owns a gzFile so an exception cannot leak it; close() reports the failure a destructor would swallow.
class GzOut {
public:
    explicit GzOut(const std::string& path) : path_(path), file_(gzopen(path.c_str(), "wb")) {
        if (file_ == nullptr) throw std::runtime_error("Failed to write compressed output: " + path);
    }
    GzOut(const GzOut&) = delete;
    GzOut& operator=(const GzOut&) = delete;
    ~GzOut() {
        if (file_ != nullptr) gzclose(file_);
    }
    gzFile get() const { return file_; }
    const std::string& path() const { return path_; }
    void close() {
        const int rc = gzclose(file_);
        file_ = nullptr;
        if (rc != Z_OK) throw std::runtime_error("Failed to close compressed output: " + path_);
    }

private:
    std::string path_;
    gzFile file_ = nullptr;
};

// Orientation-aware edge key for a step pair, e.g. "12+>13-". Adjacency-aware so a
// tandem self-loop (the same edge traversed repeatedly) shows up as a high count.
std::string edge_key(const PathStep& a, const PathStep& b) {
    return a.node_id + (a.reverse ? '-' : '+') + ">" + b.node_id + (b.reverse ? '-' : '+');
}

// A walk as oriented node ids, "1+,7-,8+", in the direction it is spelled.
std::string walk_text(const std::vector<PathStep>& steps) {
    std::string out;
    for (const PathStep& step : steps) {
        if (!out.empty()) out.push_back(',');
        out += step.node_id;
        out.push_back(step.reverse ? '-' : '+');
    }
    return out;
}

// The same walk read on the other strand: steps reversed, each orientation flipped. It spells the
// reverse complement of what the original spells.
std::vector<PathStep> reverse_walk(const std::vector<PathStep>& steps) {
    std::vector<PathStep> out(steps.rbegin(), steps.rend());
    for (PathStep& step : out) step.reverse = !step.reverse;
    return out;
}

// One path's traversal of a bubble: node counts and edge counts.
struct InspectPathRow {
    std::string name;
    std::size_t sequence_length = 0;
    std::unordered_map<std::string, InspectNodeCount> node_counts;
    std::unordered_map<std::string, std::size_t> edge_counts;
};

// Path clustering by source->sink walk: summarize a walk as an oriented, bp-weighted token multiset
// (token "<node><strand>" weighted by node bp), compared with weighted Jaccard sum(min)/sum(max).
// Order-insensitive; captures inversion (strand) and copy number (repeat weight). Sorted map so the
// merge is a linear two-pointer pass.
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

// Fast clustering: each walk -> bottom-k MinHash sketch over oriented node-step shingles; sketch
// Jaccard estimates identity; walks within --cluster-similarity are unioned (disjoint-set), clusters
// = connected components. The sketch is multiplicity-aware (a shingle's occurrence index folds into
// the hash) so it approximates the multiset Jaccard - needed to separate tandem-repeat haplotypes
// sharing a unit at different copy numbers. MinHash: Broder 1997; Mash: Ondov et al. 2016.
std::uint64_t splitmix64(std::uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

// hash_step_token (oriented per-step 64-bit hash) comes from graph_utils.hpp.

constexpr std::size_t kWalkSketchShingle = 3;
constexpr std::size_t kWalkSketchSize = 512;

// `out_complete`, when given, reports whether the sketch holds EVERY shingle occurrence rather than
// the bottom k of them. A complete sketch is the walk's full shingle multiset, so two complete
// sketches can be compared exactly instead of estimated.
std::vector<std::uint64_t> build_walk_minhash_sketch(
    const std::vector<std::uint64_t>& tokens,
    std::size_t shingle_size,
    std::size_t sketch_size,
    bool* out_complete = nullptr) {

    if (out_complete != nullptr) *out_complete = false;
    if (shingle_size == 0 || sketch_size == 0 || tokens.size() < shingle_size) {
        return {};
    }
    if (out_complete != nullptr) {
        *out_complete = (tokens.size() - shingle_size + 1) <= sketch_size;
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

// Bottom-k Jaccard over the UNION of the two sketches.
//
// The obvious form -- intersection over union of the two stored sketches -- is biased whenever either
// sketch is truncated, i.e. a walk longer than kWalkSketchSize + shingle - 1 steps. Each sketch then
// holds the smallest hashes of its OWN shingle set, so a shared shingle can sit inside one sketch and
// below the other's cutoff. The bias grows with the size difference, which is exactly the tandem-array
// case: two walks of the same repeat unit at different copy numbers.
//
// The standard unbiased estimator takes the k smallest distinct values of the union, k = min(|a|,|b|),
// and asks how many are in both. Every value considered is at or below both cutoffs, so membership is
// decidable from the sketches alone.
double sketch_jaccard(const std::vector<std::uint64_t>& a, const std::vector<std::uint64_t>& b) {
    if (a.empty() || b.empty()) {
        return 0.0;
    }
    const std::size_t k = std::min(a.size(), b.size());
    std::size_t i = 0;
    std::size_t j = 0;
    std::size_t taken = 0;
    std::size_t both = 0;
    while (taken < k && (i < a.size() || j < b.size())) {
        const bool take_a = (j >= b.size()) || (i < a.size() && a[i] < b[j]);
        const bool take_b = (i >= a.size()) || (j < b.size() && b[j] < a[i]);
        if (take_a) {
            ++i;
        } else if (take_b) {
            ++j;
        } else {                      // equal: present in both
            ++i;
            ++j;
            ++both;
        }
        ++taken;
    }
    return taken == 0 ? 0.0 : static_cast<double>(both) / static_cast<double>(taken);
}

// Exact multiset Jaccard of two COMPLETE sketches.
//
// The bottom-k estimator above deliberately subsamples: it looks at k = min(|a|,|b|) values even when
// both sketches already hold every shingle the walks have. On the worked example that estimates
// J = 0.5 where the exact value sitting in memory is 0.667, and near a clustering threshold that
// difference splits or joins a pair for no reason. When neither sketch was truncated there is nothing
// to estimate -- the sketches ARE the multisets -- so intersection over union is computed directly.
double exact_multiset_jaccard(const std::vector<std::uint64_t>& a,
                              const std::vector<std::uint64_t>& b) {
    if (a.empty() || b.empty()) {
        return 0.0;
    }
    std::size_t i = 0, j = 0, both = 0;
    while (i < a.size() && j < b.size()) {
        if (a[i] < b[j]) ++i;
        else if (b[j] < a[i]) ++j;
        else { ++both; ++i; ++j; }
    }
    const std::size_t uni = a.size() + b.size() - both;
    return uni == 0 ? 0.0 : static_cast<double>(both) / static_cast<double>(uni);
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
    bool sketch_complete = false;       // sketch holds every shingle, so it can be compared exactly
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
    // This is an all-pairs O(U^2) comparison over a dense U x U matrix of doubles. That is fine for
    // the reviewed panels (hundreds of distinct walks) and becomes the dominant cost on a large
    // cohort. No LSH or banding candidate stage exists, so the limit is stated rather than hidden.
    constexpr std::size_t kWarnUniqueWalks = 2000;
    constexpr std::size_t kMaxUniqueWalks = 25000;      // ~5 GB of matrix alone
    if (n > kMaxUniqueWalks) {
        throw std::runtime_error(
            "inspect: " + std::to_string(n) + " distinct walks exceeds the clustering limit of " +
            std::to_string(kMaxUniqueWalks) + ". Clustering compares every pair over a dense " +
            "matrix, so this run would need about " +
            std::to_string((n * n * sizeof(double)) / (1024ULL * 1024ULL * 1024ULL)) +
            " GB. Restrict to one bubble with --bubble-id, or drop --cluster.");
    }
    if (n > kWarnUniqueWalks) {
        std::cerr << "[inspect] WARNING: clustering " << n << " distinct walks compares every pair "
                  << "over a dense matrix (about "
                  << (n * n * sizeof(double)) / (1024ULL * 1024ULL) << " MB and "
                  << (n * (n - 1) / 2) << " comparisons). This implementation is not cohort-scale.\n";
    }

    std::vector<std::vector<double>> dist(n, std::vector<double>(n, 0.0));
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i + 1; j < n; ++j) {
            // Three cases, in order of how much is actually known. Walks too short to shingle have
            // empty sketches and fall back to the exact bp-weighted Jaccard, so short bubbles still
            // cluster instead of coming out all-singleton. Two COMPLETE sketches are the full shingle
            // multisets, so they are compared exactly rather than estimated. Only when at least one
            // was truncated is the bottom-k estimator needed.
            double id;
            if (uniques[i].sketch.empty() || uniques[j].sketch.empty()) {
                id = weighted_jaccard_tokens(uniques[i].weights, uniques[j].weights);
            } else if (uniques[i].sketch_complete && uniques[j].sketch_complete) {
                id = estimate_identity_from_jaccard(
                    exact_multiset_jaccard(uniques[i].sketch, uniques[j].sketch));
            } else {
                id = estimate_identity_from_jaccard(
                    sketch_jaccard(uniques[i].sketch, uniques[j].sketch));
            }
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
        // The lexicographically smallest member of the medoid walk, not the first one encountered:
        // identical walks are pooled in GFA record order, so `members.front()` made the reported
        // representative -- and therefore the sort key of the whole clusters TSV -- depend on how the
        // P/W records happened to be ordered in the file.
        co.representative = *std::min_element(uniques[rep].members.begin(),
                                              uniques[rep].members.end());
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

// --genotype-table state for one bubble: each distinct spelled sequence is one allele, and every
// crossing path points at the allele it spells.
struct AlleleAssignment {
    std::unordered_map<std::string, std::size_t> index_of_sequence;
    std::vector<const std::string*> sequences;  // first-seen order; keys of index_of_sequence
    std::unordered_map<std::string, std::size_t> allele_of_path;
    bool reference_seen = false;
    bool reference_reversed = false;
    // Per allele (first-seen index): the distinct walks that spell it, each with the paths taking it,
    // in GFA path order. Almost always one walk; two walks can spell the same bases.
    std::vector<std::vector<std::pair<std::vector<PathStep>, std::vector<std::string>>>> walks_of_allele;
    // Every crossing path's walk, source->sink, in GFA path order, for the node dictionary.
    std::vector<std::pair<std::string, std::vector<PathStep>>> walks;
    // The reference path's crossing, as step indices on that path: [reference_lo, reference_hi].
    const PathRecord* reference = nullptr;
    std::size_t reference_lo = 0;
    std::size_t reference_hi = 0;
};

// Writes the allele FASTA (">alleleK hapA,hapB,...") and the row x allele table whose columns are the
// FASTA ids. Alleles are numbered by descending haplotype count, ties by first appearance in GFA path
// order, so the numbering never depends on hashing. A row is all NA when its haplotype -- or any of a
// sample's haplotypes -- does not cross the bubble: a partial sum would read as a real low count.
std::size_t write_genotype_outputs(const Graph& graph, const AlleleAssignment& alleles,
                                   const InspectEmit& emit) {
    const std::size_t n = alleles.sequences.size();
    std::vector<std::size_t> support(n, 0);
    for (const auto& entry : alleles.allele_of_path) ++support[entry.second];
    std::vector<std::size_t> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(),
                     [&](std::size_t a, std::size_t b) { return support[a] > support[b]; });
    std::vector<std::size_t> rank(n);
    for (std::size_t r = 0; r < n; ++r) rank[order[r]] = r;

    std::vector<std::string> members(n);
    for (const auto& path : graph.paths) {
        const auto it = alleles.allele_of_path.find(path.name);
        if (it == alleles.allele_of_path.end()) continue;
        std::string& list = members[rank[it->second]];
        if (!list.empty()) list.push_back(',');
        list += path.name;
    }

    GzOut fasta(emit.alleles_out_path);
    for (std::size_t r = 0; r < n; ++r) {
        const std::string& sequence = *alleles.sequences[order[r]];
        gz_write_fasta_record(fasta.get(), fasta.path(), "allele" + std::to_string(r + 1) + " " + members[r],
                              alleles.reference_reversed ? reverse_complement(sequence) : sequence);
    }
    fasta.close();

    // The allele FASTA's walks: one record per allele under the allele FASTA's own header, or, when
    // several walks spell one allele, one record per walk marked walk=i/n and listing only the
    // haplotypes taking it. Oriented exactly as the allele sequence is, so each walk spells its allele.
    GzOut walks(emit.allele_walks_out_path);
    for (std::size_t r = 0; r < n; ++r) {
        const auto& entries = alleles.walks_of_allele[order[r]];
        for (std::size_t w = 0; w < entries.size(); ++w) {
            std::string header = "allele" + std::to_string(r + 1) + " ";
            if (entries.size() == 1) {
                header += members[r];
            } else {
                std::string takers;
                for (const std::string& name : entries[w].second) takers += (takers.empty() ? "" : ",") + name;
                header += takers + " walk=" + std::to_string(w + 1) + "/" + std::to_string(entries.size());
            }
            const std::vector<PathStep> walk =
                alleles.reference_reversed ? reverse_walk(entries[w].first) : entries[w].first;
            gz_write_string(walks.get(), ">" + header + "\n" + walk_text(walk) + "\n", walks.path());
        }
    }
    walks.close();

    GzOut table(emit.genotypes_out_path);
    std::string line = emit.samples == nullptr ? "haplotype" : "sample";
    for (std::size_t r = 0; r < n; ++r) line += "\tallele" + std::to_string(r + 1);
    line += '\n';
    gz_write_string(table.get(), line, table.path());
    const auto write_row = [&](const std::string& label, const std::vector<std::string>& haplotypes) {
        std::vector<std::size_t> counts(n, 0);
        bool missing = false;
        for (const std::string& hap : haplotypes) {
            const auto it = alleles.allele_of_path.find(hap);
            if (it == alleles.allele_of_path.end()) {
                missing = true;
                break;
            }
            ++counts[rank[it->second]];
        }
        line = tsv_sanitize(label);
        for (std::size_t r = 0; r < n; ++r) {
            line += '\t';
            line += missing ? "NA" : std::to_string(counts[r]);
        }
        line += '\n';
        gz_write_string(table.get(), line, table.path());
    };
    if (emit.samples == nullptr) {
        for (const auto& path : graph.paths) write_row(path.name, {path.name});
    } else {
        for (const std::string& sample : *emit.samples) write_row(sample, emit.sample_haplotypes->at(sample));
    }
    table.close();
    return n;
}

std::string json_string(const std::string& text) {
    std::string out = "\"";
    for (const char c : text) {
        switch (c) {
            case '"': out += "\\\""; break;
            case '\\': out += "\\\\"; break;
            case '\n': out += "\\n"; break;
            case '\t': out += "\\t"; break;
            case '\r': out += "\\r"; break;
            default:
                if (static_cast<unsigned char>(c) < 0x20) {
                    char buffer[8];
                    std::snprintf(buffer, sizeof(buffer), "\\u%04x", static_cast<unsigned char>(c));
                    out += buffer;
                } else {
                    out.push_back(c);
                }
        }
    }
    return out + "\"";
}

// The bubble's node dictionary: every node on any crossing walk, with its forward-strand sequence,
// length, how many haplotypes visit it, and where it sits.
//
// Position is anchored to the REFERENCE path's own walk through the bubble (-r). A node the reference
// visits gets ref_offset, the 0-based bp offset of its first base from the first base of that walk
// (the source node included), along the reference path's direction; when the path name carries an
// interval (chr:start-end) it also gets ref_start/ref_end, 1-based on that interval. A node the
// reference does not visit has no reference coordinate: it gets the two reference nodes that flank it
// on the first walk visiting it, and the reference offsets it lies between. `order` ranks every node
// along the reference -- reference nodes at their offset, others just after their left flank -- and
// falls back to first appearance across the walks when no reference crosses the bubble.
void write_nodes_json(const Graph& graph, const Bubble& bubble, const AlleleAssignment& alleles,
                      const InspectEmit& emit) {
    std::vector<std::string> nodes;                        // first appearance across the walks
    std::unordered_map<std::string, std::size_t> first_seen;
    std::unordered_map<std::string, std::size_t> haplotypes;
    for (const auto& walk : alleles.walks) {
        std::unordered_set<std::string> here;
        for (const PathStep& step : walk.second) {
            if (first_seen.emplace(step.node_id, nodes.size()).second) nodes.push_back(step.node_id);
            if (here.insert(step.node_id).second) ++haplotypes[step.node_id];
        }
    }
    const auto length_of = [&](const std::string& id) -> std::size_t {
        const auto it = graph.nodes.find(id);
        return it == graph.nodes.end() ? 0 : it->second.sequence.size();
    };

    // The reference walk, in the reference path's own direction, and each node's first offset on it.
    std::unordered_map<std::string, std::size_t> ref_offset;
    std::unordered_map<std::string, std::size_t> ref_occurrences;
    bool has_interval = false;
    std::size_t absolute_first = 0;   // 1-based position of the walk's first base on the interval
    if (alleles.reference != nullptr) {
        const auto& steps = alleles.reference->steps;
        const ParsedReferencePath label = parse_reference_path_label(alleles.reference->name);
        has_interval = label.has_interval;
        std::size_t before = 0;
        for (std::size_t i = 0; i < alleles.reference_lo; ++i) before += length_of(steps[i].node_id);
        absolute_first = label.region_start_1based + before;
        std::size_t offset = 0;
        for (std::size_t i = alleles.reference_lo; i <= alleles.reference_hi; ++i) {
            ref_offset.emplace(steps[i].node_id, offset);
            ++ref_occurrences[steps[i].node_id];
            offset += length_of(steps[i].node_id);
        }
    }
    const bool referenced = alleles.reference != nullptr;
    const auto offset_end = [&](const std::string& id) {
        const std::size_t length = length_of(id);
        return ref_offset.at(id) + (length == 0 ? 0 : length - 1);
    };

    // Flanking reference nodes of a node the reference does not visit, on the first walk visiting it.
    struct Flank {
        std::string left;
        std::string right;
    };
    std::unordered_map<std::string, Flank> flank;
    if (referenced) {
        for (const auto& walk : alleles.walks) {
            const auto& steps = walk.second;
            for (std::size_t i = 0; i < steps.size(); ++i) {
                const std::string& id = steps[i].node_id;
                if (ref_offset.count(id) != 0 || flank.count(id) != 0) continue;
                Flank f;
                for (std::size_t j = i; j-- > 0;) {
                    if (ref_offset.count(steps[j].node_id) != 0) { f.left = steps[j].node_id; break; }
                }
                for (std::size_t j = i + 1; j < steps.size(); ++j) {
                    if (ref_offset.count(steps[j].node_id) != 0) { f.right = steps[j].node_id; break; }
                }
                // Along the reference, not along this walk: the flank with the smaller offset comes first.
                if (!f.left.empty() && !f.right.empty() && ref_offset.at(f.right) < ref_offset.at(f.left)) {
                    std::swap(f.left, f.right);
                }
                flank.emplace(id, f);
            }
        }
    }

    // `order`: along the reference when there is one, else first appearance.
    std::vector<std::string> ordered = nodes;
    if (referenced) {
        const auto key = [&](const std::string& id) {
            if (ref_offset.count(id) != 0) return std::make_tuple(ref_offset.at(id), 0, first_seen.at(id));
            const Flank& f = flank.at(id);
            const std::size_t after = !f.left.empty() ? offset_end(f.left) : 0;
            return std::make_tuple(after, 1, first_seen.at(id));
        };
        std::stable_sort(ordered.begin(), ordered.end(),
                         [&](const std::string& a, const std::string& b) { return key(a) < key(b); });
    }

    GzOut out(emit.nodes_json_out_path);
    std::string text = "{\n  \"bubble\": " + std::to_string(bubble.id) + ",\n  \"source\": " +
                       json_string(bubble.source) + ",\n  \"sink\": " + json_string(bubble.sink) +
                       ",\n  \"reference_path\": " +
                       (emit.reference_path.empty() ? std::string("null") : json_string(emit.reference_path)) +
                       ",\n  \"reference_crosses_bubble\": " + (referenced ? "true" : "false") +
                       ",\n  \"reference_strand\": " +
                       (referenced ? json_string(alleles.reference_reversed ? "sink_to_source" : "source_to_sink")
                                   : std::string("null")) +
                       ",\n  \"haplotypes_crossing\": " + std::to_string(alleles.walks.size()) + ",\n  \"nodes\": {";
    gz_write_string(out.get(), text, out.path());
    for (std::size_t r = 0; r < ordered.size(); ++r) {
        const std::string& id = ordered[r];
        const auto node = graph.nodes.find(id);
        const std::string sequence = node == graph.nodes.end() ? std::string() : node->second.sequence;
        text = std::string(r == 0 ? "\n" : ",\n") + "    " + json_string(id) + ": {\"seq\": " + json_string(sequence) +
               ", \"length\": " + std::to_string(sequence.size()) + ", \"order\": " + std::to_string(r) +
               ", \"haplotypes\": " + std::to_string(haplotypes.at(id)) + ", \"on_reference\": ";
        if (!referenced) {
            text += "null";
        } else if (ref_offset.count(id) != 0) {
            const std::size_t offset = ref_offset.at(id);
            text += "true, \"ref_offset\": " + std::to_string(offset);
            if (ref_occurrences.at(id) > 1) text += ", \"ref_occurrences\": " + std::to_string(ref_occurrences.at(id));
            if (has_interval) {
                text += ", \"ref_start\": " + std::to_string(absolute_first + offset) +
                        ", \"ref_end\": " + std::to_string(absolute_first + offset_end(id));
            }
        } else {
            const Flank& f = flank.at(id);
            const auto side = [&](const std::string& anchor, bool left) -> std::string {
                if (anchor.empty()) return "null";
                return std::to_string(left ? offset_end(anchor) : ref_offset.at(anchor));
            };
            text += "false, \"between\": [" + (f.left.empty() ? std::string("null") : json_string(f.left)) + ", " +
                    (f.right.empty() ? std::string("null") : json_string(f.right)) + "], \"between_ref_offsets\": [" +
                    side(f.left, true) + ", " + side(f.right, false) + "]";
            if (has_interval && !f.left.empty() && !f.right.empty()) {
                text += ", \"between_ref_positions\": [" + std::to_string(absolute_first + offset_end(f.left)) + ", " +
                        std::to_string(absolute_first + ref_offset.at(f.right)) + "]";
            }
        }
        text += "}";
        gz_write_string(out.get(), text, out.path());
    }
    gz_write_string(out.get(), "\n  }\n}\n", out.path());
    out.close();
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

    GzOut walks_file(emit.walks_out_path);

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
    AlleleAssignment alleles;
    try {
        for (const auto& path : graph.paths) {
            // bubble_steps, not the inside-node-only interval finder: a haplotype that crosses the
            // bubble with NO interior node -- a pure deletion, or the short side of an insertion -- has
            // no interval and was dropped, so inspect emitted one allele where `bubble` reported two.
            // It is a real allele and it is often the interesting one.
            const BubblePathIndex index = build_bubble_path_index(path);
            BubblePathInterval used{};
            std::vector<std::size_t> step_indices;
            const auto steps_opt = bubble_steps(path, index, bubble, &used,
                                                emit.genotype_table ? &step_indices : nullptr);
            if (!steps_opt.has_value() || steps_opt->empty()) {
                continue;
            }
            const std::vector<PathStep>& steps = *steps_opt;
            const BubblePathInterval* interval = &used;

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
            // The same record as a walk: its header, then the oriented nodes that spell exactly the
            // sequence above (source and sink included, as the sequence includes them), on one line.
            gz_write_string(walks_file.get(), ">" + fasta_name + "\n" + walk_text(steps) + "\n", walks_file.path());

            if (emit.genotype_table) {
                const auto slot = alleles.index_of_sequence.try_emplace(sequence, alleles.sequences.size());
                if (slot.second) alleles.sequences.push_back(&slot.first->first);
                alleles.allele_of_path.emplace(path.name, slot.first->second);
                if (alleles.walks_of_allele.size() < alleles.sequences.size()) {
                    alleles.walks_of_allele.resize(alleles.sequences.size());
                }
                auto& walks_here = alleles.walks_of_allele[slot.first->second];
                const auto same_walk = std::find_if(walks_here.begin(), walks_here.end(), [&](const auto& entry) {
                    return walk_text(entry.first) == walk_text(steps);
                });
                if (same_walk == walks_here.end()) {
                    walks_here.push_back({steps, {path.name}});
                } else {
                    same_walk->second.push_back(path.name);
                }
                alleles.walks.emplace_back(path.name, steps);
                if (path.name == emit.reference_path && !step_indices.empty()) {
                    // Every walk is spelled source->sink; the reference's own step indices run
                    // backwards exactly when it crosses sink->source, i.e. when its strand is the
                    // reverse complement of that spelling.
                    alleles.reference_seen = true;
                    alleles.reference_reversed = step_indices.front() > step_indices.back();
                    alleles.reference = &path;
                    alleles.reference_lo = *std::min_element(step_indices.begin(), step_indices.end());
                    alleles.reference_hi = *std::max_element(step_indices.begin(), step_indices.end());
                }
            }

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
                    uw.sketch = build_walk_minhash_sketch(tokens, kWalkSketchShingle, kWalkSketchSize,
                                                          &uw.sketch_complete);
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
    walks_file.close();

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
    // Closed and CHECKED, not left to the destructor. Destruction closes the file but has nowhere to
    // report a write that failed on its way to storage, so a truncated table would be committed as
    // part of a successful family. The same applies to every stream below.
    node_lengths_out.close();
    if (!node_lengths_out)
        throw std::runtime_error("Failed writing node lengths TSV: " + emit.node_lengths_out_path);

    if (emit.genotype_table) {
        cli::ensure_parent_dir_for_file(emit.genotypes_out_path);
        cli::ensure_parent_dir_for_file(emit.alleles_out_path);
        result.alleles_written = write_genotype_outputs(graph, alleles, emit);
        cli::ensure_parent_dir_for_file(emit.nodes_json_out_path);
        write_nodes_json(graph, bubble, alleles, emit);
        result.reference_seen = alleles.reference_seen;
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
        clusters_out.close();
        if (!clusters_out)
            throw std::runtime_error("Failed writing clusters TSV: " + emit.clusters_out_path);
        result.clusters_written = clusters.size();
        result.clusters_out_path = emit.clusters_out_path;
    }

    // table_out and edge_out are still open here; both are checked before returning so a late write
    // failure on either becomes an error rather than a short file inside a committed family.
    table_out.close();
    if (!table_out) throw std::runtime_error("Failed writing node counts TSV: " + table_out_path);
    edge_out.close();
    if (!edge_out) throw std::runtime_error("Failed writing edge counts TSV: " + edge_table_out_path);

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
    bool genotype_table = false;
    std::string samples_path;
    std::string reference_query;

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
        if (arg == "--genotype-table") {
            genotype_table = true;
            continue;
        }
        if (arg == "--samples") {
            samples_path = require_value(arg);
            genotype_table = true;
            continue;
        }
        if (arg == "-r" || arg == "--reference-path") {
            reference_query = require_value(arg);
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
    if (!reference_query.empty() && !genotype_table) {
        throw std::runtime_error("--reference-path orients only the --genotype-table output; pass --genotype-table "
                                 "(or --samples) with it");
    }

    ParseGfaOptions parse_options;
    parse_options.include_paths = true;
    parse_options.include_sequences = true;
    const Graph graph = parse_gfa(gfa_path, parse_options);
    // The same contract every other module applies. Without it a duplicate path name produced two
    // rows labelled identically -- nothing downstream could tell which haplotype was which -- and a
    // non-zero overlap silently inflated every path_length_bp, because inspect spells by
    // concatenating whole nodes and an overlap means those nodes share bases.
    validate_graph_paths(graph, "inspect", true, true);
    if (graph.paths.empty()) {
        throw std::runtime_error("Input GFA has no P/W paths; inspect requires paths");
    }

    std::string reference_path;
    std::vector<std::string> sample_order;
    std::unordered_map<std::string, std::vector<std::string>> sample_haplotypes;
    if (genotype_table) {
        // The allele FASTA lists haplotypes as a comma-separated description after the id, so a name
        // carrying whitespace or a comma would split into names that exist nowhere. PanSN never does;
        // refuse rather than rename, because the names are the join key back to the cosigt table.
        for (const auto& path : graph.paths) {
            if (path.name.find_first_of(", \t") != std::string::npos) {
                throw std::runtime_error("inspect --genotype-table: path name '" + path.name +
                                         "' contains whitespace or a comma and cannot be listed in the "
                                         "alleles FASTA header");
            }
        }
        if (!reference_query.empty()) {
            reference_path = resolve_reference_path_name(graph, reference_query, "inspect");
        }
        if (!samples_path.empty()) {
            std::unordered_set<std::string> graph_haps;
            for (const auto& path : graph.paths) graph_haps.insert(path.name);
            auto cosigt = read_cosigt_table(samples_path, "inspect", &graph_haps, "a path in the input graph");
            for (const auto& [hap, samples] : cosigt.first) {
                for (const std::string& sample : samples) sample_haplotypes[sample].push_back(hap);
            }
            sample_order = std::move(cosigt.second);
        }
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

    // The bubbles CSV and the GFA are separate inputs and nothing ties them together, so a CSV from
    // another graph -- or from the same graph before a rewrite renumbered it -- used to run to
    // completion: node lengths came back as 0, the count matrix had a column per phantom node and no
    // rows, and the run exited 0. Every selected bubble must name nodes this graph actually has.
    {
        std::vector<std::string> missing;
        for (const Bubble* b : selected_bubbles) {
            const auto check = [&](const std::string& n) {
                if (graph.nodes.find(n) == graph.nodes.end() && missing.size() < 8) {
                    missing.push_back("bubble " + std::to_string(b->id) + " node " + n);
                }
            };
            check(b->source);
            check(b->sink);
            for (const std::string& n : b->inside) check(n);
        }
        if (!missing.empty()) {
            throw std::runtime_error(
                "inspect: the bubbles CSV describes nodes the GFA does not contain (" +
                cli::join_with_comma(missing) +
                "); the CSV and the graph are not the same graph");
        }
    }

    // Verbose per-bubble block is useful for a single bubble; for an all-bubbles run
    // it just floods stdout, so we show a progress bar on stderr instead.
    const bool single_bubble = bubble_id != 0;

    std::size_t total_paths_written = 0;
    cli::RunLog log("inspect", quiet);
    log.info("input " + gfa_path + " (" + std::to_string(selected_bubbles.size()) +
             (selected_bubbles.size() == 1 ? " bubble)" : " bubbles)"));

    std::vector<std::string> single_outputs;
    std::string single_info;

    // Every output path this run will write, explicit and derived alike, resolved before anything is
    // opened. Checking only the three explicit flags missed the larger half of the family: an explicit
    // --table-out naming the derived <prefix>.bubble_N.node_lengths.tsv collided silently and one
    // clobbered the other, which is exactly the case a preflight exists to catch.
    {
        std::vector<std::string> finals;
        for (const Bubble* b : selected_bubbles) {
            const std::string stem = out_prefix + ".bubble_" + std::to_string(b->id);
            const std::string final_fasta =
                fasta_out_path.empty() ? stem + ".paths.fa.gz" : fasta_out_path;
            finals.push_back(final_fasta);
            std::string final_walks = stem + ".walks.fa.gz";
            if (!fasta_out_path.empty()) {
                std::string walk_stem = fasta_out_path;
                for (const char* suffix : {".gz", ".fa", ".fasta"}) {
                    const std::string sfx(suffix);
                    if (walk_stem.size() > sfx.size() &&
                        walk_stem.compare(walk_stem.size() - sfx.size(), sfx.size(), sfx) == 0) {
                        walk_stem.resize(walk_stem.size() - sfx.size());
                    }
                }
                final_walks = walk_stem + ".walks.fa.gz";
            }
            finals.push_back(final_walks);
            finals.push_back(table_out_path.empty() ? stem + ".node_counts.tsv" : table_out_path);
            finals.push_back(edge_table_out_path.empty() ? stem + ".edge_counts.tsv"
                                                         : edge_table_out_path);
            finals.push_back(stem + ".node_lengths.tsv");
            if (cluster) finals.push_back(stem + ".clusters.tsv");
            if (genotype_table) {
                finals.push_back(stem + ".genotypes.tsv.gz");
                finals.push_back(stem + ".alleles.fa.gz");
                finals.push_back(stem + ".alleles.walks.fa.gz");
                finals.push_back(stem + ".nodes.json.gz");
            }
        }
        std::vector<std::string> inputs = {gfa_path, bubbles_csv_path};
        if (!samples_path.empty()) inputs.push_back(samples_path);
        cli::reject_output_collisions("inspect", finals, inputs);
    }

    cli::ProgressBar progress((single_bubble || quiet) ? "" : "Inspecting bubbles",
                              single_bubble ? 0 : selected_bubbles.size());

    // An all-bubbles run writes five files per bubble. Failing at bubble 400 of 500 used to leave 399
    // complete-looking bubbles on disk beside a non-zero exit, which is the shape of result a later
    // command consumes without noticing.
    cli::StagedOutputs staged("inspect");

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
        // The walks companion of the paths FASTA: beside an explicit --fasta-out, else under the prefix.
        std::string final_walks = out_prefix + ".bubble_" + std::to_string(bubble.id) + ".walks.fa.gz";
        if (!fasta_out_path.empty()) {
            std::string stem = fasta_out_path;
            for (const char* suffix : {".gz", ".fa", ".fasta"}) {
                const std::string sfx(suffix);
                if (stem.size() > sfx.size() && stem.compare(stem.size() - sfx.size(), sfx.size(), sfx) == 0) {
                    stem.resize(stem.size() - sfx.size());
                }
            }
            final_walks = stem + ".walks.fa.gz";
        }
        emit.walks_out_path = staged.stage(final_walks);
        emit.clusters_out_path =
            out_prefix + ".bubble_" + std::to_string(bubble.id) + ".clusters.tsv";
        emit.cluster = cluster;
        emit.cluster_similarity = cluster_similarity;

        // The writer receives staged paths; everything reported below is the final destination.
        const std::string final_fasta = bubble_fasta_out_path;
        const std::string final_table = bubble_table_out_path;
        const std::string final_edges = bubble_edge_table_out_path;
        const std::string final_lengths = emit.node_lengths_out_path;
        const std::string final_clusters = emit.clusters_out_path;
        bubble_fasta_out_path = staged.stage(final_fasta);
        bubble_table_out_path = staged.stage(final_table);
        bubble_edge_table_out_path = staged.stage(final_edges);
        emit.node_lengths_out_path = staged.stage(final_lengths);
        emit.clusters_out_path = staged.stage(final_clusters);
        const std::string final_genotypes =
            out_prefix + ".bubble_" + std::to_string(bubble.id) + ".genotypes.tsv.gz";
        const std::string final_alleles = out_prefix + ".bubble_" + std::to_string(bubble.id) + ".alleles.fa.gz";
        const std::string final_allele_walks =
            out_prefix + ".bubble_" + std::to_string(bubble.id) + ".alleles.walks.fa.gz";
        const std::string final_nodes_json = out_prefix + ".bubble_" + std::to_string(bubble.id) + ".nodes.json.gz";
        if (genotype_table) {
            emit.genotype_table = true;
            emit.genotypes_out_path = staged.stage(final_genotypes);
            emit.alleles_out_path = staged.stage(final_alleles);
            emit.allele_walks_out_path = staged.stage(final_allele_walks);
            emit.nodes_json_out_path = staged.stage(final_nodes_json);
            emit.reference_path = reference_path;
            if (!samples_path.empty()) {
                emit.samples = &sample_order;
                emit.sample_haplotypes = &sample_haplotypes;
            }
        }

        InspectBubbleResult result = write_inspect_outputs_for_bubble(
            graph,
            bubble,
            bubble_fasta_out_path,
            bubble_table_out_path,
            bubble_edge_table_out_path,
            emit);
        result.fasta_out_path = final_fasta;
        result.table_out_path = final_table;
        result.edge_table_out_path = final_edges;
        result.node_lengths_out_path = final_lengths;
        result.clusters_out_path = final_clusters;
        if (result.paths_written == 0) {
            throw std::runtime_error(
                "inspect: no path crosses bubble " + std::to_string(bubble.id) + " (" + bubble.source +
                ".." + bubble.sink + "); the bubbles CSV does not describe this graph");
        }
        // Requiring at least one crossing closes the dangerous case -- a CSV from a different graph --
        // but not the quiet one: a STALE graph with the same node ids and only some of the original
        // paths still produces crossings, just fewer of them. The CSV records how many `bubble` saw,
        // so the two can simply be compared. A mismatch is not an error (a legitimate path-dropping
        // transform produces one), but it must not pass unremarked, because every count downstream is
        // then taken over a different panel than the one the CSV describes.
        if (bubble.path_support > 0 && result.paths_written != bubble.path_support) {
            std::cerr << "[inspect] WARNING: bubble " << bubble.id << " is crossed by "
                      << result.paths_written << " path(s), but the bubbles CSV records path_support="
                      << bubble.path_support
                      << ". The graph and the CSV describe different panels; counts here are over the "
                         "graph's paths.\n";
        }
        if (genotype_table && !reference_path.empty() && !result.reference_seen) {
            std::cerr << "[inspect] WARNING: reference path " << reference_path << " does not cross bubble "
                      << bubble.id << "; its alleles are left in source->sink orientation.\n";
        }
        total_paths_written += result.paths_written;
        progress.tick();

        if (single_bubble) {
            single_info = "bubble " + std::to_string(bubble.id) + " (source/sink " + bubble.source +
                          "/" + bubble.sink + ", " + std::to_string(bubble.inside.size()) +
                          " inside nodes, " + std::to_string(result.paths_written) + " paths" +
                          (emit.cluster ? ", " + std::to_string(result.clusters_written) + " clusters" : "") +
                          ")";
            single_outputs = {result.fasta_out_path, final_walks, result.table_out_path,
                              result.edge_table_out_path, result.node_lengths_out_path};
            if (emit.cluster) {
                single_outputs.push_back(result.clusters_out_path);
            }
            if (genotype_table) {
                single_outputs.push_back(final_genotypes);
                single_outputs.push_back(final_alleles);
                single_outputs.push_back(final_allele_walks);
                single_outputs.push_back(final_nodes_json);
            }
        }
    }
    progress.done();
    staged.commit();

    if (single_bubble) {
        log.info(single_info);
        log.wrote(single_outputs);
    } else {
        log.info("inspected " + std::to_string(selected_bubbles.size()) + " bubbles; " +
                 std::to_string(total_paths_written) + " paths written");
        log.info("wrote per-bubble outputs under " + out_prefix + ".bubble_<id>.*");
    }
    log.done();

    return 0;
}

} // namespace panvar
