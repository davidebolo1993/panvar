#pragma once

#include <cstddef>
#include <string>
#include <unordered_map>
#include <vector>

#include "panvar/gfa.hpp"
#include "panvar/ref_path.hpp"

namespace panvar {

// One gene-level GTF feature (col3 == "gene"), in reference coordinates.
struct GeneFeature {
    std::string gene_id;
    std::string gene_name;   // falls back to gene_id when the GTF row has no gene_name
    std::string biotype;
    std::string chrom;       // as written in the GTF (e.g. "6")
    std::size_t start_1based = 0;
    std::size_t end_1based = 0;
};

// True iff `name` looks like PanSN `sample#hap#contig:start-end`: at least two '#'
// separators and a parseable ":start-end" interval. Annotation is gated on this.
bool is_pansn(const std::string& name);

// Resolve a reference path in `graph` by exact name, else unique case-insensitive
// substring. Returns nullptr when not found or ambiguous (callers decide how to react).
const PathRecord* resolve_reference_path(const Graph& graph, const std::string& query);

// Parse gene features from a (optionally gzipped) Ensembl/GENCODE GTF, keeping only
// col3=="gene" rows on `ref_chrom` (lenient: a leading "chr" is stripped on both sides,
// case-insensitive, so "chr6" matches "6") whose interval overlaps [lo_1based, hi_1based].
// If no row ever matches `ref_chrom`, a warning is logged and an empty vector returned --
// never a hard error.
std::vector<GeneFeature> parse_gtf(const std::string& gtf_path,
                                   const std::string& ref_chrom,
                                   std::size_t lo_1based,
                                   std::size_t hi_1based);

// Project genes onto the reference walk's nodes: each reference node spans
// [start + prefix[k], start + prefix[k+1]) in reference bp; map node_id -> the indices
// (into `genes`) of every gene overlapping that node. Only reference nodes appear.
std::unordered_map<std::string, std::vector<int>> project_genes_to_nodes(
    const Graph& graph,
    const PathRecord& ref_path,
    const ParsedReferencePath& ref_meta,
    const std::vector<GeneFeature>& genes);

// node_id <tab> gene names (';'-joined), one row per annotated node, node-id sorted.
void write_node_genes_tsv(const std::string& output_path,
                          const std::unordered_map<std::string, std::vector<int>>& node_genes,
                          const std::vector<GeneFeature>& genes);

// Bandage-loadable CSV: Name,Colour,Gene (one row per annotated reference node, coloured
// by a deterministic per-gene palette, labelled with the gene name(s)).
void write_bandage_gene_colors_csv(const std::string& output_path,
                                   const std::unordered_map<std::string, std::vector<int>>& node_genes,
                                   const std::vector<GeneFeature>& genes);

// Per-gene copy number for a haplotype by competitive realignment: align each gene's reference
// sequence (gene_seqs, parallel to the gene set) against the haplotype sequence and assign each
// distinct target locus to the gene it aligns BEST to (highest identity). Returns per-gene copy
// counts (parallel to gene_seqs) = number of loci assigned to each gene. This separates collapsed
// paralogs (CYP2D6 vs 2D7) that graph multiplicity cannot. min_identity / min_qcov filter weak hits.
std::vector<int> assign_gene_copies(const std::vector<std::string>& gene_seqs,
                                    const std::string& hap_seq,
                                    const std::string& preset = "asm20",
                                    double min_identity = 0.90,
                                    double min_qcov = 0.50);

// One-shot annotation for bubble/panphorte: resolve the reference path in `graph` by
// `ref_query`, gate on PanSN (warn+return false otherwise), parse genes over the reference
// region, project them onto reference nodes, and write the Bandage gene CSV (and, if
// `node_genes_out` is non-empty, a node_genes TSV). Returns true when a CSV was written.
bool emit_gene_annotation(const Graph& graph,
                          const std::string& ref_query,
                          const std::string& gtf_path,
                          const std::string& bandage_csv_out,
                          const std::string& node_genes_out = "");

} // namespace panvar
