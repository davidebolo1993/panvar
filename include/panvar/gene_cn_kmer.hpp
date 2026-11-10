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
    int cn = 0;            // copy number: round(dosage) raised to `copies`, or a near-identical pair's split
    int copies = -1;       // runs of hits read as copies of the gene; -1 where none could raise the count
    int rounded = 0;       // the count from pooled hits (or the pair split) alone, before runs raise it
    bool separable = false; // false when the gene has no private k-mers (indistinguishable from a paralog)
};

// Pooled private-k-mer dosage for a folded paralog cluster. `gene_markers[g]` is gene g's discriminative
// reference (merged CDS, else gene span). Builds each gene's private canonical k-mer set, then counts
// occurrences per haplotype -> dosage -> rounded CN. A gene with no private k-mers is not separable and
// the caller reports the module total. Evidence indexed [haplotype][gene].
//
// Privacy is judged against the sibling genes and, through `backgrounds`, the rest of the locus.
//
// With `spans` (each gene's reference span, first to last coding base), the hits are also read as copies
// (`copies`), for the genes `read_runs` marks (all when not given): hits further apart than both 10 kb and the
// gene's length start a new run, so copies closer than that read as one run (and keep the pooled count). A run
// holding a fifth of the gene's private k-mers or more is a candidate, aligned only when the candidates outnumber
// the rounded count. It is one copy when the gene aligns there at 98% identity or more, the stretch it matched
// holds half the run and no copy already found, and that stretch is closer to the gene than to the other genes
// (`backgrounds`). Versions of a gene can keep very different shares of the private k-mers, so two copies of a
// low-share version add up to under 1.5 copies' worth of hits.
//
// `backgrounds`, when given, holds one sequence per gene: its siblings' full spans (introns and UTRs
// included) with the gene's own span masked. A marker k-mer occurring there is dropped as well, since
// a haplotype carrying the sibling carries it whether or not it carries the gene.
std::vector<std::vector<GeneCnEvidence>> resolve_gene_cn_kmer(
    const std::vector<std::string>& gene_markers,
    const std::vector<std::string>& hap_seqs,
    std::size_t k = 31,
    const std::vector<std::string>* backgrounds = nullptr,
    const std::vector<std::string>* spans = nullptr,
    const std::vector<char>* read_runs = nullptr);

// Routed per-gene copy-number resolver: pooled dosage (above) for divergent paralogs, but near-identical
// pairs (k-mer Jaccard over `near_identical_jaccard`) go to per-site consensus instead, since gene
// conversion smears the pooled count. Per-site splits the known module `total[hap]` by the median
// A-vs-B allele-k-mer fraction across divergent sites, so a converted site is out-voted. Only size-2
// components are routed, and only when the pooled dosages show one gene of the pair per module copy;
// larger clusters, and pairs whose two genes share every copy, stay pooled. `total[hap]` < 0 means the
// haplotype has no module total, and such a haplotype is left unresolved. `in_copies`, when given, says
// per gene whether it lies in the module's copies: the total counts those copies only, so a pair is split
// by it only when both of its genes do. A pooled count is raised to the gene's runs of hits (`copies`) where
// those are more, past the module total too: the caller reports such a gene as unreliable.
std::vector<std::vector<GeneCnEvidence>> resolve_gene_cn(
    const std::vector<std::string>& gene_markers,
    const std::vector<std::string>& hap_seqs,
    const std::vector<long>& total,
    double near_identical_jaccard = 0.5,
    std::size_t k = 31,
    const std::vector<std::string>* backgrounds = nullptr,
    const std::vector<char>* in_copies = nullptr,
    const std::vector<std::string>* spans = nullptr);

// Hybrid copies of paralog pairs, read copy by copy. For two genes of a target whose coding exons are
// homologous (`gene_exons[g]`, one sequence per merged coding exon in gene order; `gene_spans[g]` the whole
// gene), each exon of one is matched to the other's by shared 15-mers, and each exon's k-mers absent from
// every other gene's span are its private ones. On each haplotype the hits of both genes' private k-mers
// form copies, each read on its own strand (a new copy where an exon comes again, where the exon order runs
// back against the strand by more than a stray hit, or where hits lie further apart than both 10 kb and the
// longer gene), and each copy's exon reads one gene where at least half of that gene's private
// k-mers for the exon are present and under half of the other's. Only exons that tell the panel's copies
// apart count: the copies whose majority is a gene must show that gene there in 80% of cases. A pair
// with fewer than `min_exons` such exons is not assessed. A copy is a hybrid when the other gene reads on
// one run of exons at either end of it, a single crossover; other genes' exons inside a copy are
// conversion and leave it alone.
// Returns [haplotype][gene]: the haplotype's hybrid copies of the gene with a paralog, -1 where no pair of
// the gene was assessed. A gene in two pairs counts its hybrid copies of both, overlapping ones once.
std::vector<std::vector<int>> find_hybrid_copies(
    const std::vector<std::vector<std::string>>& gene_exons,
    const std::vector<std::string>& gene_spans,
    const std::vector<std::string>& hap_seqs,
    std::size_t min_exons = 5,
    std::size_t k = 31);

// Paralog groups of genes: genes whose sequences (one per gene, e.g. merged coding sequence) share at least
// `min_shared` of the smaller one's canonical k-mers are in one group, linked transitively. Returns each gene's
// group id; a gene with an empty sequence is a group of its own.
std::vector<int> paralog_groups(const std::vector<std::string>& seqs, std::size_t k = 15, double min_shared = 0.5);

} // namespace panvar
