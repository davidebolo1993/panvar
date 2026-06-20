#include "panvar/gtf.hpp"

#include "panvar/cli_utils.hpp"
#include "panvar/graph_utils.hpp"
#include "panvar/gz_reader.hpp"
#include "panvar/minimap2_align.hpp"

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <functional>
#include <iostream>
#include <stdexcept>

namespace panvar {
namespace {

std::string to_lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
}

// Strip a leading "chr" (case-insensitive) and lowercase, so "chr6" and "6" compare equal.
std::string norm_chrom(const std::string& c) {
    std::string s = to_lower(c);
    if (s.size() > 3 && s.compare(0, 3, "chr") == 0) s.erase(0, 3);
    return s;
}

// Extract a double-quoted GTF attribute value: key "value";  -> value (empty if absent).
std::string gtf_attr(const std::string& attrs, const std::string& key) {
    const std::string needle = key + " \"";
    std::size_t p = attrs.find(needle);
    if (p == std::string::npos) return std::string();
    p += needle.size();
    const std::size_t q = attrs.find('"', p);
    if (q == std::string::npos) return std::string();
    return attrs.substr(p, q - p);
}

std::vector<std::string> split_tab(const std::string& line) {
    std::vector<std::string> out;
    std::size_t start = 0;
    while (true) {
        const std::size_t tab = line.find('\t', start);
        if (tab == std::string::npos) { out.push_back(line.substr(start)); break; }
        out.push_back(line.substr(start, tab - start));
        start = tab + 1;
    }
    return out;
}

// Numeric-aware node-id ordering (ids are typically integers), lexicographic fallback.
bool node_id_less(const std::string& a, const std::string& b) {
    const bool da = !a.empty() && std::all_of(a.begin(), a.end(), [](unsigned char c){ return std::isdigit(c); });
    const bool db = !b.empty() && std::all_of(b.begin(), b.end(), [](unsigned char c){ return std::isdigit(c); });
    if (da && db) {
        if (a.size() != b.size()) return a.size() < b.size();
        return a < b;
    }
    return a < b;
}

// A small, visually distinct qualitative palette (ColorBrewer Set1/Set2/Dark2 mix).
const char* gene_color(const std::string& gene_name) {
    static const char* kPalette[] = {
        "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628",
        "#F781BF", "#1B9E77", "#D95F02", "#7570B3", "#66A61E", "#E6AB02"};
    std::uint64_t h = 1469598103934665603ULL; // FNV-1a
    for (const unsigned char c : gene_name) { h ^= c; h *= 1099511628211ULL; }
    return kPalette[h % (sizeof(kPalette) / sizeof(kPalette[0]))];
}

} // namespace

bool is_pansn(const std::string& name) {
    const std::size_t first = name.find('#');
    if (first == std::string::npos) return false;
    if (name.find('#', first + 1) == std::string::npos) return false;
    return parse_reference_path_label(name).has_interval;
}

const PathRecord* resolve_reference_path(const Graph& graph, const std::string& query) {
    for (const PathRecord& p : graph.paths) {
        if (p.name == query) return &p;
    }
    const std::string needle = to_lower(query);
    const PathRecord* hit = nullptr;
    int n = 0;
    for (const PathRecord& p : graph.paths) {
        if (to_lower(p.name).find(needle) != std::string::npos) { hit = &p; ++n; }
    }
    return n == 1 ? hit : nullptr;
}

std::vector<GeneFeature> parse_gtf(const std::string& gtf_path,
                                   const std::string& ref_chrom,
                                   std::size_t lo_1based,
                                   std::size_t hi_1based) {
    GzLineReader reader(gtf_path);
    if (!reader.ok()) {
        throw std::runtime_error("Failed to open GTF: " + gtf_path);
    }
    const std::string want = norm_chrom(ref_chrom);
    std::vector<GeneFeature> genes;
    bool chrom_seen = false;
    std::string line;
    while (reader.getline(line)) {
        if (line.empty() || line[0] == '#') continue;
        const std::vector<std::string> f = split_tab(line);
        if (f.size() < 9) continue;
        if (f[2] != "gene") continue;
        if (norm_chrom(f[0]) != want) continue;
        chrom_seen = true;
        std::size_t start = 0, end = 0;
        try {
            start = static_cast<std::size_t>(std::stoull(f[3]));
            end = static_cast<std::size_t>(std::stoull(f[4]));
        } catch (const std::exception&) { continue; }
        if (end < lo_1based || start > hi_1based) continue; // no overlap with the region
        const std::string biotype = gtf_attr(f[8], "gene_biotype");
        if (biotype == "lncRNA") continue;  // long ncRNAs span whole loci and blur paralog annotation
        GeneFeature g;
        g.chrom = f[0];
        g.start_1based = start;
        g.end_1based = end;
        g.gene_id = gtf_attr(f[8], "gene_id");
        g.gene_name = gtf_attr(f[8], "gene_name");
        if (g.gene_name.empty()) g.gene_name = g.gene_id;
        g.biotype = biotype;
        genes.push_back(std::move(g));
    }
    if (!chrom_seen) {
        std::cerr << "warning: gtf chrom '" << ref_chrom << "' not found in " << gtf_path
                  << "; continuing without gene annotation\n";
    }
    return genes;
}

std::unordered_map<std::string, std::vector<int>> project_genes_to_nodes(
    const Graph& graph,
    const PathRecord& ref_path,
    const ParsedReferencePath& ref_meta,
    const std::vector<GeneFeature>& genes) {

    std::unordered_map<std::string, std::vector<int>> node_genes;
    if (genes.empty()) return node_genes;
    const std::vector<std::size_t> pref = path_prefix_bp(ref_path, graph.nodes);
    const std::size_t base = ref_meta.region_start_1based;
    for (std::size_t k = 0; k < ref_path.steps.size(); ++k) {
        const std::string& id = ref_path.steps[k].node_id;
        const std::size_t node_lo = base + pref[k];          // 1-based start of this node
        const std::size_t node_hi = base + pref[k + 1] - 1;  // inclusive end
        std::vector<int>& hits = node_genes[id];
        for (int gi = 0; gi < static_cast<int>(genes.size()); ++gi) {
            const GeneFeature& g = genes[gi];
            if (g.end_1based < node_lo || g.start_1based > node_hi) continue;
            if (std::find(hits.begin(), hits.end(), gi) == hits.end()) hits.push_back(gi);
        }
        if (hits.empty()) node_genes.erase(id);
    }
    return node_genes;
}

namespace {
std::vector<std::string> sorted_annotated_nodes(
    const std::unordered_map<std::string, std::vector<int>>& node_genes) {
    std::vector<std::string> ids;
    ids.reserve(node_genes.size());
    for (const auto& kv : node_genes) ids.push_back(kv.first);
    std::sort(ids.begin(), ids.end(), node_id_less);
    return ids;
}
} // namespace

void write_node_genes_tsv(const std::string& output_path,
                          const std::unordered_map<std::string, std::vector<int>>& node_genes,
                          const std::vector<GeneFeature>& genes) {
    cli::ensure_parent_dir_for_file(output_path);
    std::ofstream out(output_path);
    if (!out) throw std::runtime_error("Failed to write node_genes TSV: " + output_path);
    out << "node_id\tgenes\n";
    for (const std::string& id : sorted_annotated_nodes(node_genes)) {
        out << id << '\t';
        const std::vector<int>& gi = node_genes.at(id);
        for (std::size_t i = 0; i < gi.size(); ++i) {
            if (i) out << ';';
            out << genes[gi[i]].gene_name;
        }
        out << '\n';
    }
}

void write_bandage_gene_colors_csv(const std::string& output_path,
                                   const std::unordered_map<std::string, std::vector<int>>& node_genes,
                                   const std::vector<GeneFeature>& genes) {
    cli::ensure_parent_dir_for_file(output_path);
    std::ofstream out(output_path);
    if (!out) throw std::runtime_error("Failed to write Bandage gene CSV: " + output_path);
    out << "Name,Colour,Gene\n";
    for (const std::string& id : sorted_annotated_nodes(node_genes)) {
        const std::vector<int>& gi = node_genes.at(id);
        std::string label;
        for (std::size_t i = 0; i < gi.size(); ++i) {
            if (i) label += ';';
            label += genes[gi[i]].gene_name;
        }
        out << id << ',' << gene_color(genes[gi.front()].gene_name) << ',' << label << '\n';
    }
}

namespace {

// Total length covered by a set of [lo,hi) intervals after merging overlaps.
std::size_t union_len(std::vector<std::pair<std::size_t, std::size_t>> iv) {
    if (iv.empty()) return 0;
    std::sort(iv.begin(), iv.end());
    std::size_t total = 0, cur_lo = iv[0].first, cur_hi = iv[0].second;
    for (std::size_t i = 1; i < iv.size(); ++i) {
        if (iv[i].first <= cur_hi) { cur_hi = std::max(cur_hi, iv[i].second); }
        else { total += cur_hi - cur_lo; cur_lo = iv[i].first; cur_hi = iv[i].second; }
    }
    return total + (cur_hi - cur_lo);
}

double interval_overlap_frac(std::size_t a0, std::size_t a1,
                             const std::vector<std::pair<std::size_t, std::size_t>>& iv) {
    if (a1 <= a0) return 0.0;
    std::size_t ov = 0;
    for (const auto& b : iv) {
        const std::size_t lo = std::max(a0, b.first), hi = std::min(a1, b.second);
        if (hi > lo) ov += hi - lo;
    }
    return static_cast<double>(ov) / static_cast<double>(a1 - a0);
}

// One chained copy locus of a single gene on the haplotype.
struct CopyLocus { std::size_t t0, t1; double identity; double qcov; };

// Chain a gene's minimap2 hits into copy loci: hits that are target-adjacent but cover DISJOINT
// query regions belong to the same copy (a short copy split around a missing insertion); hits that
// reuse the same query region are SEPARATE copies (e.g. tandem-adjacent full copies).
std::vector<CopyLocus> chain_gene_hits(const std::vector<Minimap2Hit>& hits, std::size_t gene_len,
                                       double min_identity) {
    // Detect copies with gap-compressed identity (so a short form whose large indel deflates block
    // identity still passes), but carry BLOCK identity for the cross-gene competition (so near-identical
    // paralogs do not steal each other's loci just because indels were discounted).
    struct H { std::size_t qs, qe, ts, te; double id; };  // id = block identity, for competition
    std::vector<H> hs;
    for (const Minimap2Hit& h : hits) {
        if (!h.ok || h.aln_block_len == 0 || h.gc_identity < min_identity) continue;
        hs.push_back({h.query_start_bp, h.query_end_bp, h.target_start_bp, h.target_end_bp, h.identity()});
    }
    std::sort(hs.begin(), hs.end(), [](const H& a, const H& b) { return a.ts < b.ts; });
    struct Cluster { std::size_t t0, t1; std::vector<std::pair<std::size_t, std::size_t>> q; double idsum, idw; };
    std::vector<Cluster> clusters;
    const std::size_t max_gap = gene_len;  // target adjacency window; query-overlap is the real guard
    for (const H& h : hs) {
        Cluster* dst = nullptr;
        for (Cluster& c : clusters) {
            const bool target_adjacent = h.ts <= c.t1 + max_gap;
            const bool query_disjoint = interval_overlap_frac(h.qs, h.qe, c.q) < 0.5;
            if (target_adjacent && query_disjoint) { dst = &c; break; }
        }
        const double w = static_cast<double>(h.qe - h.qs);
        if (dst == nullptr) {
            clusters.push_back({h.ts, h.te, {{h.qs, h.qe}}, h.id * w, w});
        } else {
            dst->t1 = std::max(dst->t1, h.te);
            dst->t0 = std::min(dst->t0, h.ts);
            dst->q.emplace_back(h.qs, h.qe);
            dst->idsum += h.id * w; dst->idw += w;
        }
    }
    std::vector<CopyLocus> loci;
    for (const Cluster& c : clusters) {
        const double qcov = gene_len ? static_cast<double>(union_len(c.q)) / static_cast<double>(gene_len) : 0.0;
        loci.push_back({c.t0, c.t1, c.idw > 0 ? c.idsum / c.idw : 0.0, qcov});
    }
    return loci;
}

} // namespace

std::vector<int> assign_gene_copies(const std::vector<std::string>& gene_seqs,
                                    const std::string& hap_seq,
                                    const std::string& preset,
                                    double min_identity,
                                    double min_qcov) {
    std::vector<int> counts(gene_seqs.size(), 0);
    if (hap_seq.empty()) return counts;

    struct Cand { int gene; std::size_t t0, t1; double identity; };
    std::vector<Cand> cands;
    for (std::size_t gi = 0; gi < gene_seqs.size(); ++gi) {
        const std::string& gs = gene_seqs[gi];
        if (gs.empty()) continue;
        const std::vector<Minimap2Hit> hits = minimap2_hits("gene", gs, "hap", hap_seq, preset);
        for (const CopyLocus& L : chain_gene_hits(hits, gs.size(), min_identity)) {
            if (L.qcov < min_qcov) continue;
            cands.push_back({static_cast<int>(gi), L.t0, L.t1, L.identity});
        }
    }
    // Competitive assignment: claim loci in descending (block) identity, skipping any that overlaps an
    // already-claimed locus by >50% of the shorter interval; the winning gene gets the copy.
    std::sort(cands.begin(), cands.end(), [](const Cand& a, const Cand& b) { return a.identity > b.identity; });
    std::vector<std::pair<std::size_t, std::size_t>> claimed;
    for (const Cand& c : cands) {
        const std::size_t clen = c.t1 > c.t0 ? c.t1 - c.t0 : 1;
        bool overlaps = false;
        for (const auto& cl : claimed) {
            const std::size_t ov_lo = std::max(c.t0, cl.first), ov_hi = std::min(c.t1, cl.second);
            if (ov_hi <= ov_lo) continue;
            const std::size_t llen = cl.second > cl.first ? cl.second - cl.first : 1;
            if (static_cast<double>(ov_hi - ov_lo) / static_cast<double>(std::min(clen, llen)) > 0.5) {
                overlaps = true; break;
            }
        }
        if (overlaps) continue;
        claimed.emplace_back(c.t0, c.t1);
        ++counts[c.gene];
    }
    return counts;
}

std::vector<int> gene_collapse_groups(const std::vector<std::string>& gene_seqs,
                                      const std::string& preset,
                                      double max_identity) {
    const std::size_t n = gene_seqs.size();
    std::vector<int> group(n);
    for (std::size_t i = 0; i < n; ++i) group[i] = static_cast<int>(i);
    std::function<int(int)> find = [&](int x) { while (group[x] != x) { group[x] = group[group[x]]; x = group[x]; } return x; };
    for (std::size_t i = 0; i < n; ++i) {
        if (gene_seqs[i].empty()) continue;
        for (std::size_t j = i + 1; j < n; ++j) {
            if (gene_seqs[j].empty()) continue;
            // Near-identical over most of the shorter gene (block identity) => same collapse group.
            const Minimap2Hit h = minimap2_best_hit("gj", gene_seqs[j], "gi", gene_seqs[i], preset);
            if (!h.ok) continue;
            const double qcov = static_cast<double>(h.query_end_bp - h.query_start_bp) /
                                static_cast<double>(gene_seqs[j].size());
            if (h.identity() > max_identity && qcov > 0.5) group[find(static_cast<int>(i))] = find(static_cast<int>(j));
        }
    }
    for (std::size_t i = 0; i < n; ++i) group[i] = find(static_cast<int>(i));  // canonical root per gene
    return group;
}

bool emit_gene_annotation(const Graph& graph,
                          const std::string& ref_query,
                          const std::string& gtf_path,
                          const std::string& bandage_csv_out,
                          const std::string& node_genes_out) {
    const PathRecord* ref_path = resolve_reference_path(graph, ref_query);
    if (ref_path == nullptr) {
        std::cerr << "warning: --gtf given but reference path '" << ref_query
                  << "' not found/ambiguous; skipping gene annotation\n";
        return false;
    }
    if (!is_pansn(ref_path->name)) {
        std::cerr << "warning: --gtf given but reference path '" << ref_path->name
                  << "' is not PanSN (sample#hap#contig:start-end); skipping gene annotation\n";
        return false;
    }
    const ParsedReferencePath meta = parse_reference_path_label(ref_path->name);
    const std::vector<std::size_t> pref = path_prefix_bp(*ref_path, graph.nodes);
    const std::size_t lo = meta.region_start_1based;
    const std::size_t hi = lo + (pref.empty() ? 0 : pref.back()) - 1;
    const std::vector<GeneFeature> genes = parse_gtf(gtf_path, meta.chrom, lo, hi);
    const auto node_genes = project_genes_to_nodes(graph, *ref_path, meta, genes);
    write_bandage_gene_colors_csv(bandage_csv_out, node_genes, genes);
    if (!node_genes_out.empty()) write_node_genes_tsv(node_genes_out, node_genes, genes);
    return true;
}

} // namespace panvar
