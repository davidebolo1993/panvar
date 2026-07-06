#include "panvar/gtf.hpp"

#include "panvar/cli_utils.hpp"
#include "panvar/graph_utils.hpp"
#include "panvar/gz_reader.hpp"

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <fstream>
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
    std::unordered_map<std::string, std::size_t> gene_by_id;                 // gene_id -> index in genes
    std::unordered_map<std::string, std::vector<std::pair<std::size_t, std::size_t>>> cds_by_id;
    bool chrom_seen = false;
    std::string line;
    while (reader.getline(line)) {
        if (line.empty() || line[0] == '#') continue;
        const std::vector<std::string> f = split_tab(line);
        if (f.size() < 9) continue;
        const bool is_gene = (f[2] == "gene"), is_cds = (f[2] == "CDS");
        if (!is_gene && !is_cds) continue;
        if (norm_chrom(f[0]) != want) continue;
        std::size_t start = 0, end = 0;
        try {
            start = static_cast<std::size_t>(std::stoull(f[3]));
            end = static_cast<std::size_t>(std::stoull(f[4]));
        } catch (const std::exception&) { continue; }
        if (is_cds) {
            // Collect coding-exon intervals keyed by gene_id; attached to their gene after the pass.
            const std::string gid = gtf_attr(f[8], "gene_id");
            if (!gid.empty()) cds_by_id[gid].emplace_back(start, end);
            continue;
        }
        chrom_seen = true;
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
        if (!g.gene_id.empty()) gene_by_id[g.gene_id] = genes.size();
        genes.push_back(std::move(g));
    }
    // Attach merged (sorted, overlap-unioned) CDS intervals to each kept gene.
    for (auto& kv : cds_by_id) {
        const auto git = gene_by_id.find(kv.first);
        if (git == gene_by_id.end()) continue;   // CDS of a gene outside the region / filtered out
        std::vector<std::pair<std::size_t, std::size_t>>& iv = kv.second;
        std::sort(iv.begin(), iv.end());
        std::vector<std::pair<std::size_t, std::size_t>> merged;
        for (const auto& p : iv) {
            if (!merged.empty() && p.first <= merged.back().second + 1)
                merged.back().second = std::max(merged.back().second, p.second);
            else
                merged.push_back(p);
        }
        genes[git->second].cds = std::move(merged);
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
