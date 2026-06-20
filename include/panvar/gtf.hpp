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

// Per-gene copy number for a haplotype by competitive realignment. Each gene's reference sequence
// (gene_seqs, parallel to the gene set) is aligned to the haplotype; same-gene hits that are
// target-adjacent but query-DISJOINT are chained into ONE copy (so a copy split around a missing
// insertion still counts once), while query-OVERLAPPING hits are distinct copies. Detection uses a
// gap-compressed identity floor (so a structurally truncated copy is not rejected for the indel);
// each copy is assigned to the gene it aligns BEST to by BLOCK identity (so near-identical paralogs do
// not steal each other's copies), deduplicated across genes by >50% target overlap. Returns per-gene
// copy counts. Meaningful only for genes in a reliable (singleton) collapse group.
std::vector<int> assign_gene_copies(const std::vector<std::string>& gene_seqs,
                                    const std::string& hap_seq,
                                    const std::string& preset = "asm20",
                                    double min_identity = 0.90,
                                    double min_qcov = 0.50);

// Cluster genes into collapse groups: two genes are in the same group when one's reference aligns to
// the other at > max_identity (block) over most of its length (near-identical paralogs the realignment
// cannot tell apart — C4A/C4B ~99.9%, GSTM1/GSTM2 ~100%). Returns a group id per gene. A gene alone in
// its group is RELIABLE (its per-gene copy number is trustworthy); a group of >1 is UNRELIABLE and
// should be reported as a combined total over the collapsing genes, not split.
std::vector<int> gene_collapse_groups(const std::vector<std::string>& gene_seqs,
                                      const std::string& preset = "asm20",
                                      double max_identity = 0.98);

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
