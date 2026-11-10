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
#include <cstdio>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <zlib.h>
#if !defined(_WIN32)
#include <sys/stat.h>
#endif

namespace panvar {
namespace {

void print_inspect_help() {
    std::cout
        << "Usage:\n"
        << "  panvar inspect -i <graph.gfa> (-b <prefix> | -c <bubbles.csv>) [options]\n\n"
        << "Write, for each site, every haplotype's allele (sequence and node walk) and how many times\n"
        << "it walks each node and link, for checking a site or drawing it as a heatmap.\n\n"
        << "Options:\n"
        << "  -i, --gfa <path>                  Sorted GFA from bubble, panphorte or refine (required)\n"
        << "  -b, --bubble-prefix <prefix>      Sites from <prefix>.bubbles.csv (required unless -c)\n"
        << "  -c, --bubbles-csv <path>          Site table, by path\n"
        << "  -o, --out-prefix <prefix>         Output folder: <prefix>/bubble_<N>/ per site, files named\n"
        << "                                    <name>.bubble_<N>.* (<name> the folder's name; default: inspect)\n"
        << "      --bubble-id <N>               Only this site; repeatable (default: all)\n"
        << "      --cluster                     Also group crossings with similar walks,\n"
        << "                                    <name>.bubble_<N>.clusters.tsv\n"
        << "      --cluster-similarity <f>      Share of bases a walk has in common with its group's\n"
        << "                                    representative, for --cluster; above 1 a percentage (default 0.90)\n"
        << "      --genotype-table              Also write an allele table (a row per haplotype, a column per\n"
        << "                                    distinct allele), the alleles, a node dictionary and\n"
        << "                                    <prefix>/<name>.sites.tsv, who crosses each site\n"
        << "      --samples <cosigt.tsv>        Allele table rows are samples (sample<TAB>hap1<TAB>hap2...,\n"
        << "                                    as for describe); implies --genotype-table\n"
        << "  -r, --reference-path <name>       With --genotype-table: alleles on this path's strand, node\n"
        << "                                    positions on it\n"
        << "      --fasta-out <path>            Crossings FASTA path, walks beside it; one --bubble-id only\n"
        << "      --table-out <path>            Node-count table path; one --bubble-id only\n"
        << "      --edge-table-out <path>       Link-count table path; one --bubble-id only\n"
        << "                                    (each gzip-compressed when its name ends in .gz)\n"
        << "  -q, --quiet                       No progress logs\n"
        << "  -h, --help                        Show this help\n";
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

bool ends_with_gz(const std::string& path) {   // in any case: out.FA.GZ is gzip too
    if (path.size() < 3) return false;
    const std::string end = path.substr(path.size() - 3);
    return end[0] == '.' && (end[1] == 'g' || end[1] == 'G') && (end[2] == 'z' || end[2] == 'Z');
}

// One output file, gzip-compressed or plain. Callers decide by the final name (ends in .gz), since the
// path written to is a staged one.
class TextOut {
public:
    TextOut(const std::string& path, bool gz) : path_(path), gz_(gz) {
        if (gz_) {
            file_ = gzopen(path.c_str(), "wb");
            if (file_ == nullptr) throw std::runtime_error("inspect: cannot write " + path);
        } else {
            plain_.open(path);
            if (!plain_) throw std::runtime_error("inspect: cannot write " + path);
        }
    }
    ~TextOut() { if (file_ != nullptr) gzclose(file_); }
    TextOut(const TextOut&) = delete;
    TextOut& operator=(const TextOut&) = delete;
    void write(const std::string& text) {
        if (text.empty()) return;
        if (gz_) {
            if (gzwrite(file_, text.data(), static_cast<unsigned int>(text.size())) != static_cast<int>(text.size()))
                throw std::runtime_error("inspect: failed writing " + path_);
        } else {
            plain_ << text;
            if (!plain_) throw std::runtime_error("inspect: failed writing " + path_);
        }
    }
    // Closed and checked: a write that failed on its way to storage must not be committed as complete.
    void close() {
        if (gz_) {
            gzFile f = file_;
            file_ = nullptr;
            if (gzclose(f) != Z_OK) throw std::runtime_error("inspect: failed writing " + path_);
        } else {
            plain_.close();
            if (!plain_) throw std::runtime_error("inspect: failed writing " + path_);
        }
    }

private:
    std::string path_;
    bool gz_;
    gzFile file_ = nullptr;
    std::ofstream plain_;
};

struct InspectBubbleResult {
    std::size_t paths_written = 0;      // rows: one per crossing
    std::size_t paths_crossing = 0;     // distinct paths with at least one crossing
    std::size_t paths_touching_only = 0;  // paths on both boundaries with no crossing in the table's orientation
    std::size_t clusters_written = 0;
    std::size_t alleles_written = 0;
    std::size_t samples_na = 0;          // --samples: sample rows that are NA
    bool reference_crosses = false;      // -r: the reference has a crossing
    bool read_without_orientations = false;   // crossings read by the fallback (see bubble_crossings)
    std::string fasta_out_path;
    std::string table_out_path;
    std::string edge_table_out_path;
    std::string node_lengths_out_path;
    std::string clusters_out_path;
};

// Optional sidecar outputs requested for a bubble, and how each output is compressed (by its final
// name, since the writer is handed staged paths).
struct InspectEmit {
    std::string node_lengths_out_path;  // always written
    std::string walks_out_path;         // always written: the paths FASTA's records as node walks
    std::string clusters_out_path;      // written only when cluster is true
    bool cluster = false;
    double cluster_similarity = 0.90;
    bool fasta_gz = true, table_gz = false, edge_table_gz = false;
    // --genotype-table
    bool genotype_table = false;
    std::string alleles_out_path, allele_walks_out_path, nodes_json_out_path, genotypes_out_path;
    const PathRecord* reference = nullptr;            // -r: orients the alleles and places the nodes
    std::size_t reference_index = 0;                  // its index in graph.paths
    const ParsedReferencePath* reference_meta = nullptr;
    const std::vector<std::size_t>* reference_prefix = nullptr;   // bases before each reference step
    const CosigtTable* samples = nullptr;             // --samples: rows are samples
    const std::vector<char>* path_reversed = nullptr; // per path: stored on the reverse strand
};

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

// The same walk on the other strand: steps reversed, each orientation flipped. It spells the reverse
// complement of what the original spells.
std::vector<PathStep> reverse_walk(const std::vector<PathStep>& steps) {
    std::vector<PathStep> out(steps.rbegin(), steps.rend());
    for (PathStep& step : out) step.reverse = !step.reverse;
    return out;
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

// Node ids in graph order: a sorted graph numbers its nodes along the reference, so numeric ids compare as
// numbers (text order would put 10 before 9); other ids compare as text, after the numeric ones.
bool node_id_less(const std::string& a, const std::string& b) {
    const bool na = !a.empty() && a.find_first_not_of("0123456789") == std::string::npos;
    const bool nb = !b.empty() && b.find_first_not_of("0123456789") == std::string::npos;
    if (na != nb) return na;
    if (na && a.size() != b.size()) return a.size() < b.size();
    return a < b;
}

// Edge columns in the order of their first node, then their second, by node_id_less.
bool edge_key_less(const std::string& a, const std::string& b) {
    const auto parts = [](const std::string& k) {
        const std::size_t gt = k.find('>');
        return std::make_pair(k.substr(0, gt - 1), k.substr(gt + 1, k.size() - gt - 2));
    };
    const auto [a1, a2] = parts(a);
    const auto [b1, b2] = parts(b);
    if (a1 != b1) return node_id_less(a1, b1);
    if (a2 != b2) return node_id_less(a2, b2);
    return a < b;
}

// Oriented edge key for a step pair, e.g. "12+>13-". Counted per adjacent pair, so a tandem self-loop
// (the same edge walked repeatedly) shows up as a high count.
std::string edge_key(const PathStep& a, const PathStep& b) {
    return a.node_id + (a.reverse ? '-' : '+') + ">" + b.node_id + (b.reverse ? '-' : '+');
}

bool looks_like_crossing_label(const std::string& name) {
    const std::size_t at = name.rfind('@');
    return at != std::string::npos && at + 1 < name.size() &&
           name.find_first_not_of("0123456789", at + 1) == std::string::npos;
}

// A row's label, used by rows, FASTA records and cluster members: the path name, with "@k" appended when
// the path crosses the bubble more than once (k the crossing's 1-based rank in path order). A name that
// itself ends in "@<number>" always takes the suffix ("x@1" crossing once is "x@1@1"), so every label reads
// back to one (path, crossing) and "x@1" cannot collide with the first crossing of a path "x".
std::string crossing_label(const std::string& name, std::size_t index, std::size_t count) {
    return count > 1 || looks_like_crossing_label(name) ? name + "@" + std::to_string(index + 1) : name;
}

// One crossing of a bubble by a path: node counts and edge counts.
struct InspectPathRow {
    std::string name;                   // crossing_label
    std::size_t crossing = 1;           // 1-based, in path order
    std::size_t crossings = 1;
    std::size_t sequence_length = 0;
    std::unordered_map<std::string, InspectNodeCount> node_counts;
    std::unordered_map<std::string, std::size_t> edge_counts;
};

// Walk clustering: two walks compare by the bases they share between the site's boundaries, node by node:
// each oriented node a walk steps on is an element weighing its length, its n-th visit a distinct element
// (so a repeat unit at different copy numbers differs). The similarity is twice the shared bases over the
// bases of both walks (cluster_unique_walks), so a difference counts by its length: a long insertion
// outweighs any number of single-base differences, and an extra copy weighs the copy. The crossing's first
// and last steps, the boundaries every crossing shares, are left out, or long boundaries would make every
// pair alike.
std::uint64_t splitmix64(std::uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

// hash_step_token (oriented per-step 64-bit hash) comes from graph_utils.hpp.

// A walk's elements, (oriented node visit, bases), sorted, and their total bases.
struct WalkBases {
    std::vector<std::pair<std::uint64_t, std::uint64_t>> elements;
    std::uint64_t bases = 0;
};

WalkBases walk_bases(const Graph& graph, const std::vector<PathStep>& steps) {
    WalkBases out;
    out.elements.reserve(steps.size());
    std::unordered_map<std::uint64_t, std::uint32_t> visits;
    for (std::size_t i = 1; i + 1 < steps.size(); ++i) {
        const PathStep& step = steps[i];
        const std::uint64_t t = hash_step_token(step);
        const std::uint32_t n = visits[t]++;
        const auto node = graph.nodes.find(step.node_id);
        const std::uint64_t bases = std::max<std::uint64_t>(1, node == graph.nodes.end() ? 0 : node->second.sequence.size());
        out.elements.emplace_back(splitmix64(t + 0x9e3779b97f4a7c15ULL * (n + 1)), bases);
        out.bases += bases;
    }
    std::sort(out.elements.begin(), out.elements.end());
    return out;
}

// Twice the bases on elements both walks hold, over the bases of both walks.
double shared_bases(const WalkBases& a, const WalkBases& b) {
    std::uint64_t shared = 0;
    for (std::size_t i = 0, j = 0; i < a.elements.size() && j < b.elements.size();) {
        if (a.elements[i].first < b.elements[j].first) ++i;
        else if (b.elements[j].first < a.elements[i].first) ++j;
        else { shared += a.elements[i].second; ++i; ++j; }
    }
    const std::uint64_t both = a.bases + b.bases;
    return both == 0 ? 1.0 : 2.0 * static_cast<double>(shared) / static_cast<double>(both);   // two empty interiors
}

// One distinct canonical walk and the paths realizing it.
struct UniqueWalk {
    std::string signature;
    WalkBases bases;                    // for clustering
    std::vector<std::string> members;   // row labels (see crossing_label), in encounter order
};

struct ClusterOut {
    std::size_t cluster_id = 0;
    std::size_t n_paths = 0;
    std::string representative;
    std::vector<std::string> members;
};

// Representative clustering. Walks are taken most-carried first (ties by signature); each joins the
// cluster whose representative it is most similar to, at or above the threshold, or else founds one.
// Every member is within the threshold of its representative, and membership is not transitive, so a
// copy-number ladder of near-identical walks does not chain into one cluster.
std::vector<ClusterOut> cluster_unique_walks(
    const std::vector<UniqueWalk>& uniques, double threshold) {

    const std::size_t n = uniques.size();
    if (n == 0) {
        return {};
    }
    const auto identity = [&](std::size_t i, std::size_t j) {
        return shared_bases(uniques[i].bases, uniques[j].bases);
    };
    // Each walk is compared with the representatives only: U x C comparisons, no pairwise matrix.
    constexpr std::size_t kWarnComparisons = 50000000;
    std::vector<std::size_t> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
        if (uniques[a].members.size() != uniques[b].members.size())
            return uniques[a].members.size() > uniques[b].members.size();
        return uniques[a].signature < uniques[b].signature;
    });
    const double eps = 1e-12;
    std::vector<std::size_t> reps;                    // representative walk of each cluster
    std::vector<std::vector<std::size_t>> groups;     // walks of each cluster
    std::size_t comparisons = 0;
    bool warned = false;
    for (const std::size_t w : order) {
        std::size_t best = reps.size();
        double best_id = -1.0;
        for (std::size_t c = 0; c < reps.size(); ++c) {
            const double id = identity(reps[c], w);
            ++comparisons;
            if (id + eps >= threshold && id > best_id + eps) { best = c; best_id = id; }
        }
        if (!warned && comparisons > kWarnComparisons) {
            std::cerr << cli::log_prefix("inspect") << "warning: clustering " << n << " distinct walks needs over "
                      << kWarnComparisons << " walk comparisons; restrict with --bubble-id or drop --cluster\n";
            warned = true;
        }
        if (best == reps.size()) {
            reps.push_back(w);
            groups.push_back({w});
        } else {
            groups[best].push_back(w);
        }
    }

    std::vector<ClusterOut> out;
    out.reserve(groups.size());
    for (std::size_t c = 0; c < groups.size(); ++c) {
        ClusterOut co;
        for (const std::size_t m : groups[c]) {
            co.n_paths += uniques[m].members.size();
            for (const auto& name : uniques[m].members) co.members.push_back(name);
        }
        // The representative walk's smallest member label, so the name reported does not depend on the
        // order of the GFA's path records.
        co.representative = *std::min_element(uniques[reps[c]].members.begin(), uniques[reps[c]].members.end());
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

// --genotype-table state for one bubble. A haplotype's allele is what it spells across the site: its
// crossings' sequences in forward-strand order (a reverse-stored path's are turned back), usually one,
// several where the path crosses more than once. Each distinct allele is one column.
struct AlleleAssignment {
    std::unordered_map<std::string, std::size_t> index_of_key;       // allele key -> first-seen index
    std::vector<std::vector<std::string>> sequences;                 // per allele: its crossings' sequences
    std::unordered_map<std::string, std::size_t> allele_of_path;     // path name -> allele
    // Per allele: the distinct walk tuples that spell it, each with the paths taking it, in GFA path order.
    std::vector<std::vector<std::pair<std::vector<std::vector<PathStep>>, std::vector<std::string>>>> walks_of_allele;
    // Every crossing path's walks, for the node dictionary.
    std::vector<std::pair<std::string, std::vector<std::vector<PathStep>>>> walks;
    bool reference_seen = false;
    bool reference_reversed = false;   // the reference crosses the site sink -> source
    std::size_t reference_lo = 0, reference_hi = 0;   // the reference's (first) crossing, as its step indices
};

// The allele FASTA (">alleleK hap1,hap2,...", a crossing per record when an allele has several), its walks,
// and the row x allele table whose columns are the FASTA ids. Alleles are numbered by descending haplotype
// count, ties by first appearance in GFA path order, so the numbering never depends on hashing. A row is NA
// when its haplotype -- or any of a sample's haplotypes -- does not cross the site: a partial sum would read
// as a real low count. Under -r, alleles and walks are on the reference's strand.
std::size_t write_genotype_outputs(const Graph& graph, const AlleleAssignment& alleles, const InspectEmit& emit,
                                   std::size_t& samples_na) {
    const std::size_t n = alleles.sequences.size();
    std::vector<std::size_t> support(n, 0);
    for (const auto& entry : alleles.allele_of_path) ++support[entry.second];
    std::vector<std::size_t> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return support[a] > support[b]; });
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
    // On the reference's strand, each crossing reads as the reference crosses the site, and a tuple's
    // crossings, kept in forward-strand order, come in the order the reference meets them: reversed when the
    // reference is stored on the reverse strand. The two differ where the site runs against the forward strand.
    const bool rc = alleles.reference_reversed;
    const bool flip = alleles.reference_seen && emit.path_reversed != nullptr &&
                      (*emit.path_reversed)[emit.reference_index];
    const auto oriented = [&](std::vector<std::string> seqs) {
        if (flip) std::reverse(seqs.begin(), seqs.end());
        if (rc) for (std::string& x : seqs) x = reverse_complement(x);
        return seqs;
    };
    const auto oriented_walks = [&](std::vector<std::vector<PathStep>> ws) {
        if (flip) std::reverse(ws.begin(), ws.end());
        if (rc) for (auto& w : ws) w = reverse_walk(w);
        return ws;
    };
    TextOut fasta(emit.alleles_out_path, true);
    for (std::size_t r = 0; r < n; ++r) {
        const std::vector<std::string> seqs = oriented(alleles.sequences[order[r]]);
        for (std::size_t c = 0; c < seqs.size(); ++c) {
            std::string header = ">allele" + std::to_string(r + 1) + " " + members[r];
            if (seqs.size() > 1) header += " crossing=" + std::to_string(c + 1) + "/" + std::to_string(seqs.size());
            fasta.write(header + "\n");
            for (std::size_t i = 0; i < seqs[c].size(); i += 80) fasta.write(seqs[c].substr(i, 80) + "\n");
        }
    }
    fasta.close();
    // One record per allele under the allele FASTA's own header -- or, when several walks spell one allele,
    // one per walk, marked walk=i/n and listing only the haplotypes taking it -- oriented as the allele is,
    // so each walk spells its allele. Several crossings: a record per crossing, marked as in the allele FASTA.
    TextOut walks(emit.allele_walks_out_path, true);
    for (std::size_t r = 0; r < n; ++r) {
        const auto& entries = alleles.walks_of_allele[order[r]];
        for (std::size_t w = 0; w < entries.size(); ++w) {
            std::string header = ">allele" + std::to_string(r + 1) + " ";
            if (entries.size() == 1) header += members[r];
            else {
                std::string takers;
                for (const std::string& name : entries[w].second) takers += (takers.empty() ? "" : ",") + name;
                header += takers + " walk=" + std::to_string(w + 1) + "/" + std::to_string(entries.size());
            }
            const std::vector<std::vector<PathStep>> cws = oriented_walks(entries[w].first);
            for (std::size_t c = 0; c < cws.size(); ++c) {
                std::string h = header;
                if (cws.size() > 1) h += " crossing=" + std::to_string(c + 1) + "/" + std::to_string(cws.size());
                walks.write(h + "\n" + walk_text(cws[c]) + "\n");
            }
        }
    }
    walks.close();
    TextOut table(emit.genotypes_out_path, true);
    std::string line = emit.samples == nullptr ? "haplotype" : "sample";
    for (std::size_t r = 0; r < n; ++r) line += "\tallele" + std::to_string(r + 1);
    table.write(line + "\n");
    const auto write_row = [&](const std::string& label, const std::vector<const std::string*>& haps, bool missing) {
        std::vector<std::size_t> counts(n, 0);
        for (const std::string* hap : haps) {
            const auto it = alleles.allele_of_path.find(*hap);
            if (it == alleles.allele_of_path.end()) { missing = true; break; }
            ++counts[rank[it->second]];
        }
        std::string row = tsv_sanitize(label);
        for (std::size_t r = 0; r < n; ++r) row += '\t' + (missing ? std::string("NA") : std::to_string(counts[r]));
        table.write(row + "\n");
    };
    if (emit.samples == nullptr) {
        for (const auto& path : graph.paths) write_row(path.name, {&path.name}, false);
    } else {
        // Each of a sample's haplotypes is the one contig of its slot that crosses the site; none, or two,
        // leaves the sample unknown here -- the rule describe --samples uses.
        const auto crosses = [&](const std::string& c) { return alleles.allele_of_path.count(c) != 0; };
        for (const std::string& sample : emit.samples->sample_order) {
            std::vector<const std::string*> haps;
            bool missing = false;
            for (const std::vector<std::string>& slot : emit.samples->slots.at(sample)) {
                const std::string* c = slot_contig(slot, crosses);
                if (c == nullptr) { missing = true; break; }
                haps.push_back(c);
            }
            samples_na += missing;
            write_row(sample, haps, missing);
        }
    }
    table.close();
    return n;
}

// The site's node dictionary: every node on any crossing walk, with its forward-strand sequence, length,
// how many haplotypes visit it, and where it sits. A node on the -r path's (first) crossing has ref_offset,
// the 0-based bp offset of its first base from the crossing's first base (its first boundary node
// included), along the reference path's direction, plus ref_start/ref_end when the path name carries an
// interval (read as call reads it). Any other node has the two reference nodes flanking it on the first
// walk visiting it, and the offsets between them. `order` ranks nodes along the reference (others just
// after their left flank), or by first appearance when no reference crosses the site.
void write_nodes_json(const Graph& graph, const Bubble& bubble, const AlleleAssignment& alleles,
                      const InspectEmit& emit) {
    std::vector<std::string> nodes;
    std::unordered_map<std::string, std::size_t> first_seen, haplotypes;
    for (const auto& walk : alleles.walks) {
        std::unordered_set<std::string> here;
        for (const auto& cw : walk.second)
            for (const PathStep& step : cw) {
                if (first_seen.emplace(step.node_id, nodes.size()).second) nodes.push_back(step.node_id);
                if (here.insert(step.node_id).second) ++haplotypes[step.node_id];
            }
    }
    const auto length_of = [&](const std::string& id) -> std::size_t {
        const auto it = graph.nodes.find(id);
        return it == graph.nodes.end() ? 0 : it->second.sequence.size();
    };
    std::unordered_map<std::string, std::size_t> ref_offset, ref_occurrences;
    const bool referenced = emit.reference != nullptr && alleles.reference_seen;
    const bool has_interval = referenced && emit.reference_meta->has_interval;
    std::size_t absolute_first = 0;   // genomic position of the crossing's first base
    if (referenced) {
        const auto& steps = emit.reference->steps;
        const std::size_t before = (*emit.reference_prefix)[alleles.reference_lo];
        absolute_first = emit.reference_meta->region_start_1based + before;
        std::size_t offset = 0;
        for (std::size_t i = alleles.reference_lo; i <= alleles.reference_hi; ++i) {
            ref_offset.emplace(steps[i].node_id, offset);
            ++ref_occurrences[steps[i].node_id];
            offset += length_of(steps[i].node_id);
        }
    }
    const auto offset_end = [&](const std::string& id) {
        const std::size_t length = length_of(id);
        return ref_offset.at(id) + (length == 0 ? 0 : length - 1);
    };
    struct Flank { std::string left, right; };
    std::unordered_map<std::string, Flank> flank;
    if (referenced) {
        for (const auto& walk : alleles.walks)
            for (const auto& steps : walk.second)
                for (std::size_t i = 0; i < steps.size(); ++i) {
                    const std::string& id = steps[i].node_id;
                    if (ref_offset.count(id) != 0 || flank.count(id) != 0) continue;
                    Flank f;
                    for (std::size_t j = i; j-- > 0;)
                        if (ref_offset.count(steps[j].node_id) != 0) { f.left = steps[j].node_id; break; }
                    for (std::size_t j = i + 1; j < steps.size(); ++j)
                        if (ref_offset.count(steps[j].node_id) != 0) { f.right = steps[j].node_id; break; }
                    // Along the reference, not along this walk: the flank with the smaller offset first.
                    if (!f.left.empty() && !f.right.empty() && ref_offset.at(f.right) < ref_offset.at(f.left))
                        std::swap(f.left, f.right);
                    flank.emplace(id, f);
                }
    }
    std::vector<std::string> ordered = nodes;
    if (referenced) {
        const auto key = [&](const std::string& id) {
            if (ref_offset.count(id) != 0) return std::make_tuple(ref_offset.at(id), 0, first_seen.at(id));
            const Flank& f = flank.at(id);
            return std::make_tuple(!f.left.empty() ? offset_end(f.left) : std::size_t{0}, 1, first_seen.at(id));
        };
        std::stable_sort(ordered.begin(), ordered.end(),
                         [&](const std::string& a, const std::string& b) { return key(a) < key(b); });
    }
    TextOut out(emit.nodes_json_out_path, true);
    out.write("{\n  \"bubble\": " + std::to_string(bubble.id) + ",\n  \"source\": " + json_string(bubble.source) +
              ",\n  \"sink\": " + json_string(bubble.sink) + ",\n  \"reference_path\": " +
              (emit.reference == nullptr ? std::string("null") : json_string(emit.reference->name)) +
              ",\n  \"reference_crosses_bubble\": " + (referenced ? "true" : "false") + ",\n  \"reference_strand\": " +
              (referenced ? json_string(alleles.reference_reversed ? "sink_to_source" : "source_to_sink")
                          : std::string("null")) +
              ",\n  \"haplotypes_crossing\": " + std::to_string(alleles.walks.size()) + ",\n  \"nodes\": {");
    for (std::size_t r = 0; r < ordered.size(); ++r) {
        const std::string& id = ordered[r];
        const auto node = graph.nodes.find(id);
        const std::string sequence = node == graph.nodes.end() ? std::string() : node->second.sequence;
        std::string text = std::string(r == 0 ? "\n" : ",\n") + "    " + json_string(id) + ": {\"seq\": " +
                           json_string(sequence) + ", \"length\": " + std::to_string(sequence.size()) +
                           ", \"order\": " + std::to_string(r) + ", \"haplotypes\": " +
                           std::to_string(haplotypes.at(id)) + ", \"on_reference\": ";
        if (!referenced) {
            text += "null";
        } else if (ref_offset.count(id) != 0) {
            const std::size_t offset = ref_offset.at(id);
            text += "true, \"ref_offset\": " + std::to_string(offset);
            if (ref_occurrences.at(id) > 1) text += ", \"ref_occurrences\": " + std::to_string(ref_occurrences.at(id));
            if (has_interval)
                text += ", \"ref_start\": " + std::to_string(absolute_first + offset) +
                        ", \"ref_end\": " + std::to_string(absolute_first + offset_end(id));
        } else {
            const Flank& f = flank.at(id);
            const auto side = [&](const std::string& anchor, bool left) -> std::string {
                if (anchor.empty()) return "null";
                return std::to_string(left ? offset_end(anchor) : ref_offset.at(anchor));
            };
            text += "false, \"between\": [" + (f.left.empty() ? std::string("null") : json_string(f.left)) + ", " +
                    (f.right.empty() ? std::string("null") : json_string(f.right)) + "], \"between_ref_offsets\": [" +
                    side(f.left, true) + ", " + side(f.right, false) + "]";
            if (has_interval && !f.left.empty() && !f.right.empty())
                text += ", \"between_ref_positions\": [" + std::to_string(absolute_first + offset_end(f.left)) + ", " +
                        std::to_string(absolute_first + ref_offset.at(f.right)) + "]";
        }
        out.write(text + "}");
    }
    out.write("\n  }\n}\n");
    out.close();
}

InspectBubbleResult write_inspect_outputs_for_bubble(
    const Graph& graph,
    const std::vector<BubblePathIndex>& path_indexes,
    const Bubble& bubble,
    const std::string& fasta_out_path,
    const std::string& table_out_path,
    const std::string& edge_table_out_path,
    const InspectEmit& emit) {

    cli::ensure_parent_dir_for_file(fasta_out_path);
    cli::ensure_parent_dir_for_file(table_out_path);
    cli::ensure_parent_dir_for_file(edge_table_out_path);

    TextOut fasta(fasta_out_path, emit.fasta_gz);
    TextOut walks_file(emit.walks_out_path, true);
    AlleleAssignment alleles;

    std::unordered_set<std::string> inside_nodes;
    inside_nodes.reserve(bubble.inside.size() * 2);
    inside_nodes.insert(bubble.inside.begin(), bubble.inside.end());
    // Node columns in graph order, so a heatmap reads along the site.
    std::vector<std::string> node_columns = bubble.inside;
    std::sort(node_columns.begin(), node_columns.end(), node_id_less);

    // Pass 1: collect per-crossing node and edge counts and stream the FASTA. Edge columns are not known
    // up front, so the tables are written after the union is built.
    std::vector<InspectPathRow> rows;
    std::vector<std::string> edge_order;          // edge columns; sorted into graph order before writing
    std::unordered_set<std::string> edge_seen;
    InspectBubbleResult result;
    std::unordered_set<std::string> labels_seen;
    // Distinct canonical walks, collected only when clustering is requested.
    std::unordered_map<std::string, std::size_t> sig_to_unique;
    std::vector<UniqueWalk> uniques;
    // When the table's orientations hold, a path that only touches the boundaries from outside is no
    // allele (see bubble_crossings); such paths are counted and reported rather than left as missing rows.
    const bool crossing_fallback = !bubble_orientation_holds(graph.paths, path_indexes, bubble);
    result.read_without_orientations = crossing_fallback;
    {
        for (std::size_t pi = 0; pi < graph.paths.size(); ++pi) {
            const PathRecord& path = graph.paths[pi];
            // Each crossing on its own (bubble_crossings), not one window from the first boundary to the
            // last: a path crossing twice gets a row per crossing, and the stretch between crossings is
            // part of neither. A crossing with no interior node (a pure deletion, the short side of an
            // insertion) is a crossing like any other.
            std::vector<std::pair<BubblePathInterval, std::vector<PathStep>>> crossings;
            std::vector<std::vector<char>> joined;   // per crossing, per step: reached from the step before
            for (const BubblePathInterval& iv : bubble_crossings(path, path_indexes[pi], bubble, crossing_fallback)) {
                std::vector<PathStep> steps = canonical_bubble_path_steps(path, bubble, iv);
                std::vector<char> j(steps.size(), 1);
                if (bubble.inversion_parent) {
                    // An inversion site on its own nodes: the sites inside it are left out. Its other
                    // variants stay, so a crossing's own-node walk is still close to one per haplotype.
                    // A step that follows a left-out stretch has `joined` cleared, so no edge is counted
                    // across the gap.
                    std::vector<PathStep> own;
                    std::vector<char> oj;
                    bool gap = false;
                    for (std::size_t k = 0; k < steps.size(); ++k) {
                        const bool keep = k == 0 || k + 1 == steps.size() || inside_nodes.count(steps[k].node_id);
                        if (!keep) { gap = true; continue; }
                        if (!own.empty() && inner_site_step(bubble, own.back().node_id, steps[k].node_id)) gap = true;
                        own.push_back(steps[k]);
                        oj.push_back(own.size() == 1 || !gap ? 1 : 0);
                        gap = false;
                    }
                    steps = std::move(own);
                    j = std::move(oj);
                }
                if (!steps.empty()) { crossings.emplace_back(iv, std::move(steps)); joined.push_back(std::move(j)); }
            }
            if (crossings.empty()) {
                if (path_indexes[pi].positions.count(bubble.source) && path_indexes[pi].positions.count(bubble.sink))
                    ++result.paths_touching_only;
                continue;
            }
            ++result.paths_crossing;
            if (emit.genotype_table) {
                std::vector<std::string> seqs;
                std::vector<std::vector<PathStep>> ws;
                for (const auto& [iv, steps] : crossings) {
                    seqs.push_back(spell_path_steps_sequence(graph, steps));
                    ws.push_back(steps);
                }
                // A path stored on the reverse strand meets the crossings in reverse order, so its order is
                // turned back: an allele's crossings run along the forward strand whatever the path's.
                if (emit.path_reversed != nullptr && (*emit.path_reversed)[pi]) {
                    std::reverse(seqs.begin(), seqs.end());
                    std::reverse(ws.begin(), ws.end());
                }
                std::string key;
                for (const std::string& x : seqs) { key += x; key.push_back('\n'); }
                const auto slot = alleles.index_of_key.try_emplace(key, alleles.sequences.size());
                if (slot.second) {
                    alleles.sequences.push_back(seqs);
                    alleles.walks_of_allele.emplace_back();
                }
                const std::size_t a = slot.first->second;
                alleles.allele_of_path[path.name] = a;
                auto& entries = alleles.walks_of_allele[a];
                auto e = std::find_if(entries.begin(), entries.end(), [&](const auto& x) {
                    if (x.first.size() != ws.size()) return false;
                    for (std::size_t k = 0; k < ws.size(); ++k)
                        if (walk_text(x.first[k]) != walk_text(ws[k])) return false;
                    return true;
                });
                if (e == entries.end()) entries.push_back({ws, {path.name}});
                else e->second.push_back(path.name);
                alleles.walks.emplace_back(path.name, std::move(ws));
                if (emit.reference != nullptr && pi == emit.reference_index && !alleles.reference_seen) {
                    alleles.reference_seen = true;
                    alleles.reference_reversed = !crossings.front().first.source_to_sink;
                    alleles.reference_lo = std::min(crossings.front().first.left, crossings.front().first.right);
                    alleles.reference_hi = std::max(crossings.front().first.left, crossings.front().first.right);
                }
            }
            for (std::size_t ci = 0; ci < crossings.size(); ++ci) {
                const BubblePathInterval* interval = &crossings[ci].first;
                const std::vector<PathStep>& steps = crossings[ci].second;
                const std::string label = crossing_label(tsv_sanitize(path.name), ci, crossings.size());
                // Labels are unique by construction (crossing_label) once path names are; the check
                // guards that.
                if (!labels_seen.insert(label).second) {
                    throw std::runtime_error("inspect: two rows of bubble " + std::to_string(bubble.id) +
                                             " would both be labelled '" + label + "'");
                }
                // clusters.tsv lists members separated by ';', so a name holding one cannot be listed.
                if (emit.cluster && label.find(';') != std::string::npos) {
                    throw std::runtime_error("inspect: path '" + path.name + "' has a ';' in its name, which "
                                             "clusters.tsv uses to separate members; rename it or drop --cluster");
                }

                const std::string sequence = spell_path_steps_sequence(graph, steps);
                InspectPathRow row;
                row.name = label;
                row.crossing = ci + 1;
                row.crossings = crossings.size();
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
                    if (!joined[ci][i]) continue;   // a site inside an inversion site lay between them
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
                    " interval=" + std::to_string(interval->left) + "-" + std::to_string(interval->right) +
                    " crossing=" + std::to_string(row.crossing) + "/" + std::to_string(row.crossings);
                fasta.write(">" + fasta_name + "\n");
                // The same record as a walk: its header, then the oriented nodes that spell exactly the
                // sequence (source and sink included), on one line.
                walks_file.write(">" + fasta_name + "\n" + walk_text(steps) + "\n");
                constexpr std::size_t kLineWidth = 80;
                for (std::size_t i = 0; i < sequence.size(); i += kLineWidth)
                    fasta.write(sequence.substr(i, kLineWidth) + "\n");

                if (emit.cluster) {
                    const std::string sig = build_walk_signature(steps);
                    const auto it = sig_to_unique.find(sig);
                    if (it == sig_to_unique.end()) {
                        sig_to_unique.emplace(sig, uniques.size());
                        UniqueWalk uw;
                        uw.signature = sig;
                        uw.bases = walk_bases(graph, steps);
                        uw.members.push_back(row.name);
                        uniques.push_back(std::move(uw));
                    } else {
                        uniques[it->second].members.push_back(row.name);
                    }
                }

                rows.push_back(std::move(row));
            }
        }
    }
    fasta.close();
    walks_file.close();
    if (emit.genotype_table) {
        result.alleles_written = write_genotype_outputs(graph, alleles, emit, result.samples_na);
        result.reference_crosses = alleles.reference_seen;
        write_nodes_json(graph, bubble, alleles, emit);
    }

    std::sort(edge_order.begin(), edge_order.end(), edge_key_less);

    TextOut table_file(table_out_path, emit.table_gz);
    TextOut edge_file(edge_table_out_path, emit.edge_table_gz);
    std::ostringstream table_out, edge_out;
    table_out << "path_name\tpath_length_bp\tcrossing\tcrossings";
    for (const auto& node_id : node_columns) {
        table_out << "\tnode." << tsv_sanitize(node_id);
    }
    table_out << "\n";
    edge_out << "path_name\tpath_length_bp\tcrossing\tcrossings";
    for (const auto& key : edge_order) {
        edge_out << "\tedge." << key;
    }
    edge_out << "\n";

    result.fasta_out_path = fasta_out_path;
    result.table_out_path = table_out_path;
    result.edge_table_out_path = edge_table_out_path;
    result.node_lengths_out_path = emit.node_lengths_out_path;
    for (const InspectPathRow& row : rows) {
        table_out << row.name << '\t' << row.sequence_length << '\t' << row.crossing << '\t'
                  << row.crossings;
        for (const auto& node_id : node_columns) {
            const auto it = row.node_counts.find(node_id);
            if (it == row.node_counts.end()) {
                table_out << "\t0:0:0";
            } else {
                table_out << '\t' << it->second.total << ':' << it->second.forward << ':'
                          << it->second.reverse;
            }
        }
        table_out << '\n';
        table_file.write(table_out.str());
        table_out.str(std::string());

        edge_out << row.name << '\t' << row.sequence_length << '\t' << row.crossing << '\t'
                  << row.crossings;
        for (const auto& key : edge_order) {
            const auto it = row.edge_counts.find(key);
            edge_out << '\t' << (it == row.edge_counts.end() ? 0 : it->second);
        }
        edge_out << '\n';
        edge_file.write(edge_out.str());
        edge_out.str(std::string());
        ++result.paths_written;
    }
    table_file.write(table_out.str());   // the header, when there was no row
    edge_file.write(edge_out.str());

    // Per-node bp lengths, in the same order as the node.* columns above, so the
    // node coverage heatmap can length-scale its x-axis without reordering.
    cli::ensure_parent_dir_for_file(emit.node_lengths_out_path);
    std::ofstream node_lengths_out(emit.node_lengths_out_path);
    if (!node_lengths_out) {
        throw std::runtime_error("Failed to write node lengths TSV: " + emit.node_lengths_out_path);
    }
    node_lengths_out << "node_id\tlength_bp\n";
    for (const auto& node_id : node_columns) {
        const auto it = graph.nodes.find(node_id);
        const std::size_t len = it == graph.nodes.end() ? 0 : it->second.sequence.size();
        node_lengths_out << tsv_sanitize(node_id) << '\t' << len << '\n';
    }
    // Closed and checked here, as is every stream below: the destructor cannot report a failed write,
    // and a truncated table must not be committed.
    node_lengths_out.close();
    if (!node_lengths_out)
        throw std::runtime_error("Failed writing node lengths TSV: " + emit.node_lengths_out_path);

    if (emit.cluster) {
        const std::vector<ClusterOut> clusters =
            cluster_unique_walks(uniques, emit.cluster_similarity);
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

    // Both tables are closed and checked last, so a late write failure is an error, not a short file.
    table_file.close();
    edge_file.close();

    return result;
}

// Where a bubble's walks FASTA goes: beside the paths FASTA. With an explicit --fasta-out X.fa, X.fasta or
// X.fna (optionally .gz) it is X.walks.fa.gz; otherwise <prefix>.bubble_<N>.walks.fa.gz.
std::string walks_path_for(const std::string& explicit_fasta, const std::string& stem) {
    if (explicit_fasta.empty()) return stem + ".walks.fa.gz";
    std::string base = explicit_fasta;
    for (const char* sfx : {".gz", ".fa", ".fasta", ".fna"})
        if (base.size() > std::string(sfx).size() && base.compare(base.size() - std::string(sfx).size(), std::string::npos, sfx) == 0)
            base.erase(base.size() - std::string(sfx).size());
    return base + ".walks.fa.gz";
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
    std::vector<std::size_t> bubble_ids;   // --bubble-id, repeatable; empty = all bubbles
    bool cluster = false;
    double cluster_similarity = 0.90;
    bool cluster_similarity_given = false;
    bool quiet = false;
    bool genotype_table = false;
    std::string samples_path, reference_query;

    for (std::size_t i = 0; i < args.size(); ++i) {
        const std::string arg = cli::renamed_input_flag(args[i], "inspect");
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
        if (arg == "-b" || arg == "--bubble-prefix") {
            bubble_prefix_in = require_value(arg);
            continue;
        }
        if (arg == "--bubble-id") {
            const std::size_t id = cli::parse_size_arg(arg, require_value(arg));
            // Bubble ids start at 1; no --bubble-id at all means every bubble.
            if (id == 0) throw std::runtime_error("--bubble-id must be 1 or more; omit it to inspect all bubbles");
            if (std::find(bubble_ids.begin(), bubble_ids.end(), id) == bubble_ids.end()) bubble_ids.push_back(id);
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
            cluster_similarity_given = true;
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
            if (reference_query.empty()) throw std::runtime_error("-r/--reference-path needs a path name");
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
    if (cluster_similarity_given && !cluster) {
        throw std::runtime_error("--cluster-similarity applies only with --cluster");
    }
    if (!reference_query.empty() && !genotype_table)
        throw std::runtime_error("-r/--reference-path orients the --genotype-table outputs; it does nothing without them");
    if (!bubble_prefix_in.empty()) {
        const std::string derived = bubble_prefix_in + ".bubbles.csv";
        if (bubbles_csv_path.empty()) {
            bubbles_csv_path = derived;
        } else if (bubbles_csv_path != derived) {
            throw std::runtime_error(
                "Conflicting bubble inputs: --bubble-prefix resolves to '" +
                derived + "' but --bubbles-csv is '" + bubbles_csv_path + "'");
        }
    }
    if (bubbles_csv_path.empty()) {
        throw std::runtime_error("Missing required input: --bubbles-csv <path> or --bubble-prefix <prefix>");
    }
    if (bubble_ids.size() != 1 && (!fasta_out_path.empty() || !table_out_path.empty() || !edge_table_out_path.empty())) {
        throw std::runtime_error("--fasta-out/--table-out/--edge-table-out name one bubble's files and need exactly one "
                                 "--bubble-id; use --out-prefix otherwise");
    }

    cli::RunLog log("inspect", quiet);
    log.info("reading: " + gfa_path);
    ParseGfaOptions parse_options;
    parse_options.include_paths = true;
    parse_options.include_sequences = true;
    const Graph graph = parse_gfa(gfa_path, parse_options);
    log.info("input: haplotypes: " + std::to_string(graph.paths.size()) + "; nodes: " +
             std::to_string(graph.nodes.size()));
    // The path contract every module applies: unique path names, so each row names one haplotype, and
    // zero overlaps, since inspect spells by concatenating whole nodes.
    validate_graph_paths(graph, "inspect", true, true);

    // A FASTA id ends at the first whitespace, so a path name holding one would lose its end in every
    // record; every run refuses it.
    for (const PathRecord& p : graph.paths)
        if (p.name.find_first_of(" \t") != std::string::npos)
            throw std::runtime_error("inspect: path '" + p.name + "' has whitespace in its name, which would cut "
                                     "its FASTA ids");
    const PathRecord* reference = nullptr;
    std::size_t reference_index = 0;
    ParsedReferencePath reference_meta;
    std::size_t reference_total_bp = 0;
    std::vector<std::size_t> reference_prefix;
    std::optional<CosigtTable> samples;
    if (genotype_table) {
        // The allele FASTA lists haplotypes comma-separated, so --genotype-table also refuses commas.
        for (const PathRecord& p : graph.paths)
            if (p.name.find(',') != std::string::npos)
                throw std::runtime_error("inspect --genotype-table: path '" + p.name + "' has a comma in its "
                                         "name, which the allele FASTA uses to list haplotypes");
        if (!reference_query.empty()) {
            const std::string name = resolve_reference_path_name(graph, reference_query, "inspect");
            for (std::size_t i = 0; i < graph.paths.size(); ++i)
                if (graph.paths[i].name == name) { reference = &graph.paths[i]; reference_index = i; }
            // Bases before each step, computed once for every site's positions.
            reference_prefix.assign(1, 0);
            for (const PathStep& st : reference->steps) {
                reference_total_bp += graph.nodes.at(st.node_id).sequence.size();
                reference_prefix.push_back(reference_total_bp);
            }
            // The interval read as call reads it: by the path's length, a W line's name 1-based.
            reference_meta = parse_reference_path_label(reference->name, reference_total_bp, reference->type == 'W');
        }
        if (!samples_path.empty()) {
            std::unordered_set<std::string> names;
            for (const PathRecord& p : graph.paths) names.insert(p.name);
            samples = read_cosigt_table(samples_path, "inspect", &names, "a path in the input graph", quiet);
        }
    }

    std::vector<Bubble> bubbles = read_bubbles_csv(bubbles_csv_path, /*with_inversion_sites=*/true);
    // An inversion site encloses the sites inside it and is shown on its own nodes: the stretch the
    // inversion turns.
    view_inversion_sites_on_own_nodes(bubbles);
    // A site table with no site (bubble found none of --min-variant-bp) is a result, not an error: nothing is
    // inspected, and with --genotype-table sites.tsv is written with its header only.
    // Ids name the output files, so each must be unique (the reader refuses id 0).
    {
        std::unordered_set<std::size_t> seen;
        for (const Bubble& b : bubbles) {
            if (!seen.insert(b.id).second)
                throw std::runtime_error("inspect: the bubbles CSV has bubble id " + std::to_string(b.id) + " twice");
        }
    }

    std::vector<const Bubble*> selected_bubbles;
    if (bubble_ids.empty()) {
        selected_bubbles.reserve(bubbles.size());
        for (const auto& bubble : bubbles) {
            selected_bubbles.push_back(&bubble);
        }
    } else {
        for (const std::size_t id : bubble_ids) {
            const auto bubble_it = std::find_if(bubbles.begin(), bubbles.end(), [&](const Bubble& b) { return b.id == id; });
            if (bubble_it == bubbles.end()) {
                throw std::runtime_error("Bubble ID not found in bubbles CSV: " + std::to_string(id));
            }
            selected_bubbles.push_back(&(*bubble_it));
        }
    }

    // The bubbles CSV and the GFA are separate inputs, so every selected bubble must name nodes this graph
    // has: a CSV from another graph, or from a renumbered one, is refused rather than written as empty
    // tables.
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

    // One --bubble-id gets a per-bubble summary; any other run gets a progress bar instead, since a
    // summary per bubble would flood the log.
    const bool single_bubble = selected_bubbles.size() == 1 && !bubble_ids.empty();

    std::size_t total_paths_written = 0;
    std::size_t total_paths_crossing = 0;
    std::size_t total_alleles = 0;
    std::size_t total_not_crossing = 0;   // (haplotype, site) pairs with no crossing: NA rows
    log.info("sites: " + std::to_string(bubbles.size()) + " (" + bubbles_csv_path + "); inspected: " +
             std::to_string(selected_bubbles.size()));
    log.info("reference: " + (reference != nullptr ? reference->name
                                                   : std::string(genotype_table ? "none (no -r)" : "none")));

    std::vector<std::string> single_outputs;
    std::string single_info;
    std::set<std::string> run_outputs;   // every final path this run writes

    // The prefix is the run's folder: <prefix>/<name>.sites.tsv and a folder per site, <prefix>/bubble_<N>/,
    // whose files keep their full names (<name>.bubble_<N>.<suffix>, <name> the folder's name), so a file
    // copied out still says which site it is. Site files left beside the prefix (<prefix>.bubble_<N>.<suffix>)
    // are swept below too. A trailing '/' is dropped, and a prefix such as '.' takes the name of the folder it
    // points to, so no file name starts with a dot.
    while (out_prefix.size() > 1 && out_prefix.back() == '/') out_prefix.pop_back();
    std::string prefix_name = std::filesystem::path(out_prefix).filename().string();
    if (prefix_name.empty() || prefix_name == "." || prefix_name == "..") {
        std::error_code ec;
        std::filesystem::path folder = std::filesystem::absolute(out_prefix, ec).lexically_normal();
        if (folder.filename().empty()) folder = folder.parent_path();
        prefix_name = ec ? std::string() : folder.filename().string();
    }
    if (out_prefix.empty() || prefix_name.empty() || prefix_name == "." || prefix_name == "..")
        throw std::runtime_error("inspect: --out-prefix '" + out_prefix + "' names no folder");
    {
        std::error_code ec;
        if (std::filesystem::exists(out_prefix, ec) && !std::filesystem::is_directory(out_prefix, ec))
            throw std::runtime_error("inspect: --out-prefix " + out_prefix + " is a file; inspect writes a folder there");
    }
    const auto site_dir = [&](std::size_t id) { return out_prefix + "/bubble_" + std::to_string(id); };
    const auto site_stem = [&](std::size_t id) {
        return site_dir(id) + "/" + prefix_name + ".bubble_" + std::to_string(id);
    };
    // Folders this run makes, removed again when it fails (when empty), the sites' before the run's.
    struct MadeDirs {
        std::vector<std::string> dirs;
        bool keep = false;
        ~MadeDirs() {
            if (keep) return;
            std::error_code ec;
            for (const std::string& d : dirs) std::filesystem::remove(d, ec);
        }
    } made_dirs;
    for (const Bubble* b : selected_bubbles) {
        std::error_code ec;
        if (!std::filesystem::exists(site_dir(b->id), ec)) made_dirs.dirs.push_back(site_dir(b->id));
    }
    {
        std::error_code ec;
        if (!std::filesystem::exists(out_prefix, ec)) made_dirs.dirs.push_back(out_prefix);
    }

    // Every output path this run will write, explicit and derived alike, is resolved and checked for
    // collisions with each other and with the inputs before anything is opened.
    {
        std::vector<std::string> finals;
        for (const Bubble* b : selected_bubbles) {
            const std::string stem = site_stem(b->id);
            finals.push_back(fasta_out_path.empty() ? stem + ".paths.fa.gz" : fasta_out_path);
            finals.push_back(table_out_path.empty() ? stem + ".node_counts.tsv" : table_out_path);
            finals.push_back(edge_table_out_path.empty() ? stem + ".edge_counts.tsv"
                                                         : edge_table_out_path);
            finals.push_back(stem + ".node_lengths.tsv");
            finals.push_back(walks_path_for(fasta_out_path, stem));
            if (cluster) finals.push_back(stem + ".clusters.tsv");
            if (genotype_table)
                for (const char* sfx : {".alleles.fa.gz", ".alleles.walks.fa.gz", ".nodes.json.gz", ".genotypes.tsv.gz"})
                    finals.push_back(stem + sfx);
        }
        if (genotype_table) finals.push_back(out_prefix + "/" + prefix_name + ".sites.tsv");
        std::vector<std::string> inputs = {gfa_path, bubbles_csv_path};
        if (!samples_path.empty()) inputs.push_back(samples_path);
        cli::reject_output_collisions("inspect", finals, inputs);
        run_outputs.insert(finals.begin(), finals.end());
    }

    // Path indexes cover only the nodes of the sites being read (boundaries and interiors, all a crossing
    // reads). Sites are taken in batches of about an eighth of the graph's nodes, each indexed on its own, so
    // memory stays bounded when the sites cover most of the graph, at the cost of one pass over the paths per
    // batch.
    std::vector<std::size_t> batch_start;   // index into selected_bubbles where each batch begins
    {
        const std::size_t budget = std::max<std::size_t>(20000, graph.nodes.size() / 8);
        std::size_t in_batch = 0;
        for (std::size_t i = 0; i < selected_bubbles.size(); ++i) {
            const std::size_t n = selected_bubbles[i]->inside.size() + 2;
            if (i == 0 || in_batch + n > budget) { batch_start.push_back(i); in_batch = 0; }
            in_batch += n;
        }
        batch_start.push_back(selected_bubbles.size());
    }
    std::vector<BubblePathIndex> path_indexes;
    const auto index_batch = [&](std::size_t from, std::size_t to) {
        std::unordered_set<std::string> site_nodes;
        for (std::size_t i = from; i < to; ++i) {
            const Bubble* b = selected_bubbles[i];
            site_nodes.insert(b->source);
            site_nodes.insert(b->sink);
            site_nodes.insert(b->inside.begin(), b->inside.end());
        }
        path_indexes.clear();
        path_indexes.shrink_to_fit();
        path_indexes.reserve(graph.paths.size());
        for (const PathRecord& path : graph.paths) path_indexes.push_back(build_bubble_path_index(path, site_nodes));
    };
    // A path's strand: a sorted graph numbers its nodes along the reference, so a path stored
    // reverse-complemented ends on a smaller id than it starts on. Its ends are used rather than its
    // reverse steps, which a large inversion can outweigh on a forward path.
    std::vector<char> path_reversed;
    if (genotype_table) {
        path_reversed.reserve(graph.paths.size());
        for (const PathRecord& path : graph.paths)
            path_reversed.push_back(!path.steps.empty() &&
                                    node_id_less(path.steps.back().node_id, path.steps.front().node_id));
    }

    cli::ProgressBar progress((single_bubble || quiet) ? "" : cli::log_prefix("inspect") + "inspecting sites",
                              single_bubble ? 0 : selected_bubbles.size());

    // Every output is staged and committed only after the whole run succeeds, so a failed run leaves no
    // complete-looking subset of bubbles on disk.
    cli::StagedOutputs staged("inspect");

    // --genotype-table: a row per site of who has an allele there, so coverage reads without counting NA
    // rows. A sample row is NA when any of its haplotypes is.
    const std::string sites_path = out_prefix + "/" + prefix_name + ".sites.tsv";
    const std::string sites_header = "bubble_id\tsource\tsink\thaplotypes\tcrossing\tna\ttouching_only"
                                     "\treference_crosses\talleles\tsamples\tsamples_na";
    std::vector<std::pair<std::size_t, std::string>> sites_rows;   // (bubble id, row), written in id order
    std::vector<std::size_t> no_crossing, support_differs, without_orientations, reference_misses;
    std::size_t touching_only = 0;

    std::size_t next_batch = 0;
    for (std::size_t bi = 0; bi < selected_bubbles.size(); ++bi) {
        if (bi == batch_start[next_batch]) {
            index_batch(bi, batch_start[next_batch + 1]);
            ++next_batch;
        }
        const Bubble& bubble = *selected_bubbles[bi];
        std::string bubble_fasta_out_path = fasta_out_path;
        std::string bubble_table_out_path = table_out_path;
        std::string bubble_edge_table_out_path = edge_table_out_path;
        const std::string stem = site_stem(bubble.id);
        if (bubble_fasta_out_path.empty()) bubble_fasta_out_path = stem + ".paths.fa.gz";
        if (bubble_table_out_path.empty()) bubble_table_out_path = stem + ".node_counts.tsv";
        if (bubble_edge_table_out_path.empty()) bubble_edge_table_out_path = stem + ".edge_counts.tsv";

        InspectEmit emit;
        emit.node_lengths_out_path = stem + ".node_lengths.tsv";
        emit.clusters_out_path = stem + ".clusters.tsv";
        emit.cluster = cluster;
        emit.cluster_similarity = cluster_similarity;
        const std::string final_walks = walks_path_for(fasta_out_path, stem);
        emit.walks_out_path = staged.stage(final_walks);
        emit.genotype_table = genotype_table;
        if (genotype_table) {
            emit.alleles_out_path = staged.stage(stem + ".alleles.fa.gz");
            emit.allele_walks_out_path = staged.stage(stem + ".alleles.walks.fa.gz");
            emit.nodes_json_out_path = staged.stage(stem + ".nodes.json.gz");
            emit.genotypes_out_path = staged.stage(stem + ".genotypes.tsv.gz");
            emit.reference = reference;
            emit.reference_index = reference_index;
            emit.reference_meta = &reference_meta;
            emit.reference_prefix = &reference_prefix;
            emit.samples = samples ? &*samples : nullptr;
            emit.path_reversed = &path_reversed;
        }
        emit.fasta_gz = ends_with_gz(bubble_fasta_out_path);
        emit.table_gz = ends_with_gz(bubble_table_out_path);
        emit.edge_table_gz = ends_with_gz(bubble_edge_table_out_path);

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
        if (cluster) emit.clusters_out_path = staged.stage(final_clusters);

        InspectBubbleResult result = write_inspect_outputs_for_bubble(
            graph,
            path_indexes,
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
        // Collected per site, reported once per run (below). A site no path crosses has empty tables; the
        // paths crossing are compared with the table's path_support (paths, not rows), since a graph with the
        // same node ids but fewer paths still has crossings.
        const auto note = [](std::vector<std::size_t>& ids, std::size_t id) { ids.push_back(id); };
        if (result.paths_written == 0) note(no_crossing, bubble.id);
        if (bubble.path_support > 0 && result.paths_crossing != bubble.path_support) note(support_differs, bubble.id);
        if (result.read_without_orientations) note(without_orientations, bubble.id);
        if (genotype_table && reference != nullptr && !result.reference_crosses) note(reference_misses, bubble.id);
        touching_only += result.paths_touching_only;
        if (genotype_table) {
            const std::size_t n = graph.paths.size();
            sites_rows.emplace_back(
                bubble.id, std::to_string(bubble.id) + '\t' + tsv_sanitize(bubble.source) + '\t' +
                               tsv_sanitize(bubble.sink) + '\t' + std::to_string(n) + '\t' +
                               std::to_string(result.paths_crossing) + '\t' +
                               std::to_string(n - result.paths_crossing) + '\t' +
                               std::to_string(result.paths_touching_only) + '\t' +
                               (reference == nullptr ? "NA" : result.reference_crosses ? "1" : "0") + '\t' +
                               std::to_string(result.alleles_written) + '\t' +
                               (samples ? std::to_string(samples->sample_order.size()) : "NA") + '\t' +
                               (samples ? std::to_string(result.samples_na) : "NA"));
        }
        total_paths_written += result.paths_written;
        total_paths_crossing += result.paths_crossing;
        total_alleles += result.alleles_written;
        total_not_crossing += graph.paths.size() - std::min(graph.paths.size(), result.paths_crossing);
        progress.tick();

        if (single_bubble) {
            single_info = "site: " + std::to_string(bubble.id) + "; boundaries: " + bubble.source + ".." +
                          bubble.sink + "; interior nodes: " + std::to_string(bubble.inside.size()) +
                          "; haplotypes crossing: " + std::to_string(result.paths_crossing) +
                          (result.paths_written != result.paths_crossing
                               ? "; crossings: " + std::to_string(result.paths_written)
                               : std::string()) +
                          (emit.cluster ? "; clusters: " + std::to_string(result.clusters_written) : "");
            single_outputs = {result.fasta_out_path, result.table_out_path,
                              result.edge_table_out_path, result.node_lengths_out_path, final_walks};
            if (genotype_table)
                for (const char* sfx : {".alleles.fa.gz", ".alleles.walks.fa.gz", ".nodes.json.gz", ".genotypes.tsv.gz"})
                    single_outputs.push_back(stem + sfx);
            if (genotype_table) single_outputs.push_back(sites_path);
            if (emit.cluster) {
                single_outputs.push_back(result.clusters_out_path);
            }
        }
    }
    progress.done();
    // No selected site crossed by any path: the table describes another graph, and nothing is written.
    if (!no_crossing.empty() && no_crossing.size() == selected_bubbles.size())
        throw std::runtime_error("inspect: no path crosses any selected site (first: " +
                                 std::to_string(no_crossing.front()) + "); the site table does not describe this graph");
    // One line per condition, however many sites it affects; shown with -q too.
    const auto warn_sites = [&](const std::vector<std::size_t>& ids, const std::string& what) {
        if (ids.empty()) return;
        std::string first;
        for (std::size_t k = 0; k < ids.size() && k < 5; ++k) first += (k ? "," : "") + std::to_string(ids[k]);
        std::cerr << cli::log_prefix("inspect") << "warning: " << what << ": " << ids.size() << " (first: " << first
                  << (ids.size() > 5 ? ",..." : "") << ")\n";
    };
    warn_sites(no_crossing, "sites no haplotype crosses, written with empty tables");
    warn_sites(support_differs, "sites whose crossing haplotypes differ from the table's path_support (the graph "
                                "and the table describe different panels; counts here are the graph's)");
    warn_sites(without_orientations, "sites read without the table's boundary orientations (none in the table, or "
                                     "no haplotype crosses as they say): each haplotype read from the stretch "
                                     "between the boundaries holding the site's nodes");
    warn_sites(reference_misses, "sites the reference does not cross (alleles stay source to sink, nodes have no "
                                 "reference position)");
    if (touching_only > 0)
        log.info("haplotypes touching a site's boundaries without crossing it (no row): " + std::to_string(touching_only));
    if (genotype_table) {
        std::sort(sites_rows.begin(), sites_rows.end());
        std::string text = sites_header + "\n";
        for (const auto& row : sites_rows) text += row.second + "\n";
        TextOut sites(staged.stage(sites_path), false);
        sites.write(text);
        sites.close();
    }
    staged.commit();
    made_dirs.keep = true;
    // Files an earlier run left under this prefix would read as this run's, so they are removed after the
    // commit (a failed run removes nothing). Only inspect's own names (<prefix>.bubble_<id>.<any of its
    // suffixes>, <prefix>.sites.tsv) are touched, never a file this run wrote, the GFA or the bubbles CSV.
    // A run over every bubble sweeps the prefix's whole family; a --bubble-id run only its bubbles' files and
    // a sites.tsv it did not write.
    {
        static const char* const kSuffixes[] = {".paths.fa.gz", ".node_counts.tsv", ".edge_counts.tsv",
                                                ".node_lengths.tsv", ".clusters.tsv", ".walks.fa.gz",
                                                ".alleles.fa.gz", ".alleles.walks.fa.gz", ".nodes.json.gz",
                                                ".genotypes.tsv.gz"};
        std::set<std::string> candidates;
        std::set<std::string> folders;   // site folders, removed when the sweep leaves them empty
        const std::filesystem::path pfx(out_prefix);
        const std::filesystem::path beside = pfx.has_parent_path() ? pfx.parent_path() : std::filesystem::path(".");
        const auto near = [&](const std::string& name) { return pfx.has_parent_path() ? (beside / name).string() : name; };
        const std::string file_stem = prefix_name + ".bubble_";
        // An id as inspect writes it after `lead`: digits, no leading zero (ids start at 1); bubble_007 is not
        // one. The position after the id, or npos.
        const auto id_end = [](const std::string& name, const std::string& lead) -> std::size_t {
            if (name.compare(0, lead.size(), lead) != 0 || name.size() == lead.size() || name[lead.size()] == '0')
                return std::string::npos;
            const std::size_t e = name.find_first_not_of("0123456789", lead.size());
            return e == lead.size() ? std::string::npos : (e == std::string::npos ? name.size() : e);
        };
        // A site's files: in its folder, and flat beside the prefix or in a folder <prefix>.bubble_<N>/ beside it.
        const auto add_site = [&](const std::string& id) {
            const std::string name = file_stem + id;
            for (const std::string& folder : {out_prefix + "/bubble_" + id, near(name)}) {
                folders.insert(folder);
                for (const char* sfx : kSuffixes) candidates.insert(folder + "/" + name + sfx);
            }
            for (const char* sfx : kSuffixes) candidates.insert(near(name + sfx));
        };
        if (bubble_ids.empty()) {
            std::error_code ec;
            for (std::filesystem::directory_iterator it(out_prefix, ec), end; !ec && it != end; it.increment(ec)) {
                const std::string name = it->path().filename().string();
                if (id_end(name, "bubble_") == name.size()) add_site(name.substr(7));
            }
            for (std::filesystem::directory_iterator it(beside, ec), end; !ec && it != end; it.increment(ec)) {
                const std::string name = it->path().filename().string();
                const std::size_t e = id_end(name, file_stem);
                if (e != std::string::npos && (e == name.size() || name[e] == '.'))
                    add_site(name.substr(file_stem.size(), e - file_stem.size()));
            }
        } else {
            for (const Bubble* b : selected_bubbles) add_site(std::to_string(b->id));
        }
        // sites.tsv is a common name, so only a file starting with inspect's header is taken for one.
        bool stale_sites = false;
        if (!genotype_table) {
            std::ifstream in(sites_path);
            std::string first;
            stale_sites = in && std::getline(in, first) && first == sites_header;
        }
        if (stale_sites) candidates.insert(sites_path);
        {   // a sites.tsv beside the prefix
            std::ifstream in(out_prefix + ".sites.tsv");
            std::string first;
            if (in && std::getline(in, first) && first == sites_header) candidates.insert(out_prefix + ".sites.tsv");
        }
        std::size_t removed = 0;
        bool removed_sites = false;
        // A leftover is kept when it is a file this run wrote or read, under any name. Each file's identity
        // (device, inode) is read once, so the check stays linear in the files.
        std::vector<std::string> kept_files(run_outputs.begin(), run_outputs.end());
        for (const std::string& in : {gfa_path, bubbles_csv_path, samples_path})
            if (!in.empty()) kept_files.push_back(in);
#if !defined(_WIN32)
        const auto identity = [](const std::string& path) -> std::optional<std::pair<std::uint64_t, std::uint64_t>> {
            struct stat st {};
            if (::stat(path.c_str(), &st) != 0) return std::nullopt;
            return std::make_pair(static_cast<std::uint64_t>(st.st_dev), static_cast<std::uint64_t>(st.st_ino));
        };
        std::set<std::pair<std::uint64_t, std::uint64_t>> kept_ids;
        for (const std::string& f : kept_files)
            if (const auto id = identity(f)) kept_ids.insert(*id);
#endif
        for (const std::string& path : candidates) {
            std::error_code ec;
            bool keep = run_outputs.count(path) != 0;
#if !defined(_WIN32)
            if (!keep)
                if (const auto id = identity(path)) keep = kept_ids.count(*id) != 0;
#else
            for (std::size_t k = 0; k < kept_files.size() && !keep; ++k)
                keep = std::filesystem::equivalent(path, kept_files[k], ec);
#endif
            if (!keep && std::filesystem::remove(path, ec)) {
                if (path == sites_path) removed_sites = true;
                else ++removed;
            }
        }
        for (const std::string& f : folders) {   // a folder this run did not fill, now empty
            std::error_code ec;
            if (std::filesystem::is_directory(f, ec) && std::filesystem::is_empty(f, ec)) std::filesystem::remove(f, ec);
        }
        if (removed > 0)
            log.info("removed " + std::to_string(removed) + " file(s) an earlier run left at " + out_prefix);
        if (removed_sites) log.info("removed " + sites_path + ", which an earlier run left");
    }

    if (single_bubble) {
        log.info(single_info);
        // Under the run's folder (the default names), the files are listed relative to it, as with many sites.
        const std::string under = out_prefix + "/";
        bool all_under = true;
        std::string names;
        for (const std::string& f : single_outputs) {
            all_under = all_under && f.compare(0, under.size(), under) == 0;
            if (all_under) names += (names.empty() ? "" : ", ") + f.substr(under.size());
        }
        if (all_under)
            log.info("wrote " + std::to_string(single_outputs.size()) + " files in " + under + ": " + names);
        else
            log.wrote(single_outputs);
    } else {
        log.info("crossings: (haplotype, site) pairs: " + std::to_string(total_paths_crossing) +
                 "; rows: " + std::to_string(total_paths_written) + " (a haplotype crossing a site twice has two)");
        log.info(genotype_table ? "genotype tables: alleles: " + std::to_string(total_alleles) +
                                      "; (haplotype, site) pairs with no crossing (NA): " +
                                      std::to_string(total_not_crossing)
                                : std::string("genotype tables: not asked (--genotype-table)"));
        const std::size_t per_site = (run_outputs.size() - (genotype_table ? 1 : 0)) /
                                     std::max<std::size_t>(1, selected_bubbles.size());
        log.info("wrote " + std::to_string(run_outputs.size()) + " files in " + out_prefix + "/: bubble_<id>/ (" +
                 std::to_string(per_site) + " per site)" + (genotype_table ? ", " + prefix_name + ".sites.tsv" : std::string()));
    }
    log.done();

    return 0;
}

} // namespace panvar
