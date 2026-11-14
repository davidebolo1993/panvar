#!/usr/bin/env python3
"""A locus whose gene shares its first three exons with a paralog elsewhere, for the transcripts in the locus
index. The gene has five exons (the second and fourth 40 bases); the paralog, on the other strand, has the
same first three and two of its own, with other introns. Exon k-mers of the shared exons are shared; the
transcript's junctions exist in no genome, so only the transcripts tell which are the paralog's too.

    junction_fixture.py <out_dir>
writes
  ref.fa, ga.fa        the reference (the locus and the paralog) and a genome with small changes in the locus
  gp.fa                a genome with a processed pseudogene (the spliced transcript and a poly-A tail)
  gn.fa                an unrelated genome
  locus.gfa, genes.gtf the locus (the reference's and ga's copies, one node each) and its gene GJ
  locus_p.gfa          a locus elsewhere, one of whose two haplotypes carries the processed pseudogene
  locus_b.gfa, genes_b.gtf  the paralog's locus over several nodes, ga's walk written reversed and carrying a
                       SNP four bases into the fourth exon; its gene GP, on the minus strand
  tx_own.fa            a transcriptome holding the gene's transcript and unrelated ones
  tx_paralog.fa        the same with the paralog's transcript
  tx_similar.fa        the gene's record named by its symbol alone, and the paralog's transcript under a symbol
                       that differs by a suffix and under the gene's symbol with another gene id
  junction_kmers.txt   the transcript's k-mers no locus copy holds: all, those of the first two junctions, and
                       those each record of splice_changed.fa adds
  splice_changed.fa, .tsv  two transcript records with a changed splice site, and their table rows: one with
                       status splice_site_changed (the first exon joined to the third), one whose status
                       boundary_moved covers the change (the first exon joined to the fourth)
  junction_nodes_b.tsv per k-mer of GP's transcripts that no copy of its locus holds, the node of its first base
  reads.fq             120-base reads across the first two junctions: reads of the paralog's mRNA as well
  reads_own.fq         120-base reads across the last two junctions: the gene's own
  reads_other.fq       reads of the processed pseudogene running into its flank and its tail (pp*), and reads
                       across an exon skip no transcript has (skip*)
"""
import os
import random
import sys

out = sys.argv[1]
os.makedirs(out, exist_ok=True)
rng = random.Random(5)
COMP = str.maketrans("ACGT", "TGCA")


def seq(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def rc(s):
    return s.translate(COMP)[::-1]


def mutate(s, rate):
    b = list(s)
    for i in range(len(b)):
        if rng.random() < rate:
            b[i] = rng.choice([c for c in "ACGT" if c != b[i]])
    return "".join(b)


def write(name, contigs):
    with open(os.path.join(out, name + ".fa"), "w") as f:
        for cname, s in contigs:
            f.write(f">{cname}\n")
            for i in range(0, len(s), 60):
                f.write(s[i:i + 60] + "\n")


def intron(n):
    return "GT" + seq(n - 4) + "AG"


def gene_of(exons, introns):
    """The gene's sequence and its exons' intervals in it."""
    s, iv = "", []
    for k, e in enumerate(exons):
        iv.append((len(s), len(s) + len(e)))
        s += e + (introns[k] if k < len(introns) else "")
    return s, iv


canon = lambda x: min(x, rc(x))
kmers = lambda s: {canon(s[i:i + 31]) for i in range(len(s) - 30)}

e1, e2, e3, e4, e5 = seq(300), seq(40), seq(300), seq(40), seq(300)
p4, p5 = seq(40), seq(300)                               # the paralog's own last exons
ls, le = 120000, 140000                                  # the locus
gene_at = 5000
gene, gene_exons = gene_of([e1, e2, e3, e4, e5], [intron(2000), intron(1500), intron(1200), intron(900)])
paralog, par_exons = gene_of([e1, e2, e3, p4, p5], [intron(3000), intron(1000), intron(1400), intron(700)])
bs, be = 225000, 245000                                  # the paralog's locus; the paralog on the minus strand
par_at = 5000
ref = list(seq(300000))
ref[ls + gene_at:ls + gene_at + len(gene)] = list(gene)
ref[bs + par_at:bs + par_at + len(paralog)] = list(rc(paralog))
ref = "".join(ref)
hap_a = ref[ls:ls + 1000] + mutate(ref[ls + 1000:ls + gene_at], 0.003) + ref[ls + gene_at:le]   # changes outside the gene
tx = e1 + e2 + e3 + e4 + e5
ptx = e1 + e2 + e3 + p4 + p5
write("ref", [("chr1", ref)])
flank_l, flank_r = seq(100000), seq(100000)
write("gp", [("chr1", flank_l + tx + "A" * 30 + flank_r)])
write("gn", [("chr1", seq(50000))])

# the locus: one node per copy
with open(os.path.join(out, "locus.gfa"), "w") as f:
    f.write("H\tVN:Z:1.0\nS\t1\t" + ref[ls:le] + "\nS\t2\t" + hap_a + "\n")
    f.write(f"P\tref#1#chr1:{ls}-{le}\t1+\t*\nP\tga#1#chr1:{ls}-{le}\t2+\t*\n")
attr = 'gene_id "ENSG00000000001"; gene_name "GJ"; transcript_id "TJ1"; transcript_name "TJ1"; tag "Ensembl_canonical";'
with open(os.path.join(out, "genes.gtf"), "w") as f:
    at = ls + gene_at
    f.write(f"chr1\ttest\tgene\t{at + 1}\t{at + len(gene)}\t.\t+\t.\t{attr}\n")
    for a, b in gene_exons:
        f.write(f"chr1\ttest\texon\t{at + a + 1}\t{at + b}\t.\t+\t.\t{attr}\n")

# a locus elsewhere whose second haplotype holds the processed pseudogene, between two shared nodes
ps, pm, pe = 60000, 70000, 80000
with open(os.path.join(out, "locus_p.gfa"), "w") as f:
    f.write(f"H\tVN:Z:1.0\nS\t1\t{ref[ps:pm]}\nS\t2\t{ref[pm:pe]}\nS\t3\t{tx + 'A' * 30}\n")
    f.write("L\t1\t+\t2\t+\t0M\nL\t1\t+\t3\t+\t0M\nL\t3\t+\t2\t+\t0M\n")
    f.write(f"P\tref#1#chr1:{ps}-{pe}\t1+,2+\t*\nP\tgq#1#chr1:{ps}-{pe}\t1+,3+,2+\t*\n")

# the paralog's locus: nodes cut inside exons and introns, a SNP of ga four bases into the fourth exon (so
# the junction k-mers over it are ga's own), ga's walk written from its end
gb = ref[bs:be]
L = len(paralog)
gpos = lambda o: par_at + L - 1 - o                     # the locus position of the paralog's base o (its own strand)
snp = gpos(par_exons[3][0] + 4)
cuts = sorted({0, len(gb), 2000, gpos(par_exons[4][0] + 100), snp, snp + 1, gpos(par_exons[3][0] - 300),
               gpos(par_exons[2][0] + 150) + 1, gpos(par_exons[1][0] + 20), gpos(par_exons[0][0] + 120), 15000})
nodes = [(str(k + 1), a, b) for k, (a, b) in enumerate(zip(cuts, cuts[1:]))]
alt_id = str(len(nodes) + 1)
alt_base = {"A": "C", "C": "G", "G": "T", "T": "A"}[gb[snp]]
gb_a = gb[:snp] + alt_base + gb[snp + 1:]
write("ga", [("chr1", ref[:ls] + hap_a + ref[le:bs] + gb_a + ref[be:])])
ref_walk = [n for n, _, _ in nodes]
ga_walk = [alt_id if a == snp else n for n, a, _ in nodes]
with open(os.path.join(out, "locus_b.gfa"), "w") as f:
    f.write("H\tVN:Z:1.0\n")
    for n, a, b in nodes:
        f.write(f"S\t{n}\t{gb[a:b]}\n")
    f.write(f"S\t{alt_id}\t{alt_base}\n")
    walks = [[(n, "+") for n in ref_walk], [(n, "-") for n in reversed(ga_walk)]]
    for x, ox, y, oy in sorted({(x, ox, y, oy) for w in walks for (x, ox), (y, oy) in zip(w, w[1:])}):
        f.write(f"L\t{x}\t{ox}\t{y}\t{oy}\t0M\n")
    for name, w in zip(("ref", "ga"), walks):
        f.write(f"P\t{name}#1#chr1:{bs}-{be}\t{','.join(n + o for n, o in w)}\t*\n")
attr_b = 'gene_id "ENSG00000000002"; gene_name "GP"; transcript_id "TP1"; transcript_name "TP1"; tag "Ensembl_canonical";'
with open(os.path.join(out, "genes_b.gtf"), "w") as f:
    at = bs + par_at
    f.write(f"chr1\ttest\tgene\t{at + 1}\t{at + L}\t.\t-\t.\t{attr_b}\n")
    for a, b in sorted((L - b, L - a) for a, b in par_exons):
        f.write(f"chr1\ttest\texon\t{at + a + 1}\t{at + b}\t.\t-\t.\t{attr_b}\n")

# per k-mer of the paralog's transcripts that no copy of its locus holds, the node of its first base: the
# reference's transcript first, then ga's for the k-mers only it has
node_at = lambda p, walk: next(walk[k] for k, (_, a, b) in enumerate(nodes) if a <= p < b)
tx_pos = [gpos(o) for a, b in par_exons for o in range(a, b)]   # the locus position of each transcript base
genomic = kmers(gb) | kmers(gb_a)
ptx_a = "".join(rc(gb_a)[len(gb) - 1 - p] for p in tx_pos)
expected = {}
for s, walk in ((ptx, ref_walk), (ptx_a, ga_walk)):
    for t in range(len(s) - 30):
        k = canon(s[t:t + 31])
        if k not in genomic and k not in expected:
            expected[k] = node_at(tx_pos[t], walk)
with open(os.path.join(out, "junction_nodes_b.tsv"), "w") as f:
    for k in sorted(expected):
        f.write(f"{k}\t{expected[k]}\n")

# transcriptomes, named transcript|gene id|gene name: the gene's own transcript (its id with a version) and
# unrelated ones; then the paralog's too. tx_similar.fa: the gene's record by its symbol alone, and the
# paralog's transcript twice, under a symbol that differs by a suffix and under the gene's symbol with
# another gene id
others = [(f"TU{k}|ENSG000000001{k}0.1|GU{k}", seq(900)) for k in range(5)]
with open(os.path.join(out, "tx_own.fa"), "w") as f:
    for name, s in [("TJ1|ENSG00000000001.7|GJ", tx)] + others:
        f.write(f">{name}\n{s}\n")
with open(os.path.join(out, "tx_paralog.fa"), "w") as f:
    for name, s in [("TJ1|ENSG00000000001.7|GJ", tx), ("TP1|ENSG00000000002.3|GP", ptx)] + others:
        f.write(f">{name}\n{s}\n")
with open(os.path.join(out, "tx_similar.fa"), "w") as f:
    for name, s in [("TJ1 gene=GJ", tx), ("TP1 gene=GJ.2", ptx), ("TP9|ENSG00000000002.3|GJ", ptx)] + others:
        f.write(f">{name}\n{s}\n")

own = kmers(ref[ls:le]) | kmers(hap_a)
with open(os.path.join(out, "junction_kmers.txt"), "w") as f:
    f.write(f"{len(kmers(tx) - own)} {len(kmers(e1 + e2 + e3) - own)} {len(kmers(e1 + e3) - own - kmers(tx))} "
            f"{len(kmers(e1 + e4) - own - kmers(tx) - kmers(e1 + e3))}\n")


def fastq(name, reads):
    with open(os.path.join(out, name), "w") as f:
        for rname, r in reads:
            f.write(f"@{rname}\n{r}\n+\n{'I' * len(r)}\n")


fastq("reads.fq", [(f"j{k}", tx[s:s + 120]) for k, s in ((k, 300 - 70 + rng.randrange(0, 40)) for k in range(200))])
j3 = len(e1 + e2 + e3)
fastq("reads_own.fq", [(f"o{k}", tx[s:s + 120]) for k, s in ((k, j3 - 70 + rng.randrange(0, 40)) for k in range(200))])
pseudo = flank_l[-200:] + tx + "A" * 30 + flank_r[:200]
skip = e1 + e3
masked = e1 + e4
with open(os.path.join(out, "splice_changed.fa"), "w") as f:
    f.write(f">TJ9|GJ|ref#1#chr1:{ls}-{le}|1 source=projected status=splice_site_changed\n{skip}\n")
    f.write(f">TJ8|GJ|ref#1#chr1:{ls}-{le}|1 source=projected status=boundary_moved\n{masked}\n")
with open(os.path.join(out, "splice_changed.tsv"), "w") as f:   # rows of the transcripts' table, without its header
    for tid, status, moved, third, detail in (("TJ9", "splice_site_changed", 0, 2, "i1.donor:GT>GA"),
                                              ("TJ8", "boundary_moved", 1, 3, "e1.end:+1;i1.donor:GT>GA")):
        a, b = gene_exons[third]
        f.write("\t".join(map(str, [tid, tid, "ENSG00000000001", "GJ", ".", "+", f"ref#1#chr1:{ls}-{le}", 1, 1, "projected", ".",
                                     status, "+", gene_at, gene_at + b, 300 + b - a, 2, ".", ".", ".", moved, 1,
                                     f"{gene_at},{gene_at + a}", f"300,{b - a}", detail])) + "\n")
fastq("reads_other.fq",
      [(f"pp_head{k}", pseudo[140 + 2 * k:260 + 2 * k]) for k in range(20)] +                       # flank into exon 1
      [(f"pp_tail{k}", pseudo[len(pseudo) - 330 + 2 * k:len(pseudo) - 210 + 2 * k]) for k in range(20)] +   # exon 5 into the tail
      [(f"skip{k}", skip[300 - 70 + 2 * k:300 + 50 + 2 * k]) for k in range(20)])
