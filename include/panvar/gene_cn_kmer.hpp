#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace panvar {

// Per (haplotype, gene) evidence from the k-mer copy-number resolver. Written into dup_gene_cn.tsv so a
// call is human-auditable: "1 copy of gene X because its `hits` private k-mers over a set of `priv_kmers`
// give dosage ~= 1".
struct GeneCnEvidence {
    long priv_kmers = 0;   // |private k-mer set| for this gene (denominator; constant across haplotypes)
    long hits = 0;         // private-k-mer occurrences in this haplotype (with multiplicity)
    double dosage = 0.0;   // hits / priv_kmers, ~= number of copies of this gene
    int cn = 0;            // rounded copy-number estimate (round(dosage))
    bool separable = false; // false when the gene has no private k-mers (indistinguishable from a paralog)
};

// Pooled private-k-mer dosage resolver for a folded paralog cluster. `gene_markers[g]` is gene g's
// discriminative reference sequence (merged CDS, or the gene span when a gene has no CDS). Builds each
// gene's PRIVATE canonical k-mer set (k-mers unique to it vs all other genes in the group), then for each
// haplotype in `hap_seqs` counts private-k-mer occurrences -> dosage -> rounded CN. A gene with no private
// k-mers is marked not separable (its dosage is meaningless; the caller reports the module total instead).
// Returns evidence indexed [haplotype][gene]. `k` must be <= 31.
std::vector<std::vector<GeneCnEvidence>> resolve_gene_cn_kmer(
    const std::vector<std::string>& gene_markers,
    const std::vector<std::string>& hap_seqs,
    std::size_t k = 31);

// Routed per-gene copy-number resolver. Pooled private-k-mer dosage (above) for divergent paralogs, but
// NEAR-IDENTICAL pairs (whose k-mer Jaccard exceeds `near_identical_jaccard`) are resolved by PER-SITE
// CONSENSUS instead: the pooled dosage blurs on such pairs because gene-conversion mosaics smear the count,
// whereas per-site splits the known module `total[hap]` by the median across the pair's divergent sites of
// each site's A-vs-B allele-k-mer fraction (each site independently sees the whole split, so a converted
// site is out-voted). Only pairs (near-identical components of size exactly 2) are routed; larger clusters
// stay pooled. Returns evidence indexed [haplotype][gene]. `total[hap]` is the module copy number used to
// split a near-identical pair (from the coverage route).
std::vector<std::vector<GeneCnEvidence>> resolve_gene_cn(
    const std::vector<std::string>& gene_markers,
    const std::vector<std::string>& hap_seqs,
    const std::vector<long>& total,
    double near_identical_jaccard = 0.5,
    std::size_t k = 31);

} // namespace panvar
