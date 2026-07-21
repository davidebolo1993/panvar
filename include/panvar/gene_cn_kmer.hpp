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

// Pooled private-k-mer dosage for a folded paralog cluster. `gene_markers[g]` is gene g's discriminative
// reference (merged CDS, else gene span). Builds each gene's private canonical k-mer set (unique against
// the others), then counts occurrences per haplotype -> dosage -> rounded CN. A gene with no private
// k-mers is not separable and the caller reports the module total. Evidence indexed [haplotype][gene].
std::vector<std::vector<GeneCnEvidence>> resolve_gene_cn_kmer(
    const std::vector<std::string>& gene_markers,
    const std::vector<std::string>& hap_seqs,
    std::size_t k = 31);

// Routed per-gene copy-number resolver: pooled dosage (above) for divergent paralogs, but near-identical
// pairs (k-mer Jaccard over `near_identical_jaccard`) go to per-site consensus instead, since gene
// conversion smears the pooled count. Per-site splits the known module `total[hap]` by the median
// A-vs-B allele-k-mer fraction across divergent sites, so a converted site is out-voted. Only size-2
// components are routed; larger clusters stay pooled.
std::vector<std::vector<GeneCnEvidence>> resolve_gene_cn(
    const std::vector<std::string>& gene_markers,
    const std::vector<std::string>& hap_seqs,
    const std::vector<long>& total,
    double near_identical_jaccard = 0.5,
    std::size_t k = 31);

} // namespace panvar
