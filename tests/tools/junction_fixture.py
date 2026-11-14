#!/usr/bin/env python3
"""A locus whose gene has a paralog elsewhere with the same exons and other introns, for the transcripts
in the locus index: the exons' k-mers are shared with the paralog, so only the transcript's junctions are
private. A processed pseudogene (the spliced transcript) in another genome holds the junctions too.

    junction_fixture.py <out_dir>
writes ref.fa (the reference: the locus and the paralog), ga.fa (the locus with small changes), gp.fa (a
genome with the processed pseudogene), locus.gfa (the reference's and ga's copies), genes.gtf (the gene,
three exons, the middle one 40 bases), junction_kmers.txt (how many of the transcript's k-mers no locus
copy holds), and reads.fq (120-base reads across both junctions of the
transcript, so each holds two private sites only with the transcripts in the index).
"""
import os
import random
import sys

out = sys.argv[1]
os.makedirs(out, exist_ok=True)
rng = random.Random(5)


def seq(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


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


e1, e2, e3 = seq(300), seq(40), seq(300)
ls, le = 120000, 140000
gene_at = 5000                                           # in the locus
intron1, intron2 = "GT" + seq(1996) + "AG", "GT" + seq(1496) + "AG"
gene = e1 + intron1 + e2 + intron2 + e3
paralog = e1 + "GT" + seq(2996) + "AG" + e2 + "GT" + seq(996) + "AG" + e3   # same exons, other introns
ref = list(seq(300000))
ref[ls + gene_at:ls + gene_at + len(gene)] = list(gene)
ref[230000:230000 + len(paralog)] = list(paralog)
ref = "".join(ref)
hap_a = ref[ls:ls + 1000] + mutate(ref[ls + 1000:ls + gene_at], 0.003) + ref[ls + gene_at:le]   # changes outside the gene
write("ref", [("chr1", ref)])
write("ga", [("chr1", ref[:ls] + hap_a + ref[le:])])
write("gp", [("chr1", seq(100000) + e1 + e2 + e3 + "A" * 30 + seq(100000))])
with open(os.path.join(out, "locus.gfa"), "w") as f:
    f.write("H\tVN:Z:1.0\nS\t1\t" + ref[ls:le] + "\nS\t2\t" + hap_a + "\n")
    f.write(f"P\tref#1#chr1:{ls}-{le}\t1+\t*\nP\tga#1#chr1:{ls}-{le}\t2+\t*\n")
exons, at = [], ls + gene_at
for e, intron in ((e1, intron1), (e2, intron2), (e3, "")):
    exons.append((at + 1, at + len(e)))
    at += len(e) + len(intron)
attr = 'gene_id "GJ"; gene_name "GJ"; transcript_id "TJ1"; transcript_name "TJ1"; tag "Ensembl_canonical";'
with open(os.path.join(out, "genes.gtf"), "w") as f:
    f.write(f"chr1\ttest\tgene\t{exons[0][0]}\t{exons[-1][1]}\t.\t+\t.\t{attr}\n")
    for a, b in exons:
        f.write(f"chr1\ttest\texon\t{a}\t{b}\t.\t+\t.\t{attr}\n")
tx = e1 + e2 + e3
comp = str.maketrans("ACGT", "TGCA")
canon = lambda x: min(x, x.translate(comp)[::-1])
kmers = lambda s: {canon(s[i:i + 31]) for i in range(len(s) - 30)}
with open(os.path.join(out, "junction_kmers.txt"), "w") as f:   # the transcript's k-mers no locus copy holds
    f.write(f"{len(kmers(tx) - kmers(ref[ls:le]) - kmers(hap_a))}\n")
with open(os.path.join(out, "reads.fq"), "w") as f:
    for k in range(200):
        s = 300 - 70 + rng.randrange(0, 40)                # across both junctions
        r = tx[s:s + 120]
        f.write(f"@j{k}\n{r}\n+\n{'I' * len(r)}\n")
