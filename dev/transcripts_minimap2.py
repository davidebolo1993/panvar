#!/usr/bin/env python3
"""`panvar transcripts` against splice alignment: each reference transcript (spliced here from the
reference path and the GTF) is aligned to each haplotype's spelled path with minimap2 -x splice:hq, and
the sequence the alignment gives on the haplotype is compared with the HST of that transcript, haplotype
and copy (the alignment overlapping the HST most).

usage: transcripts_minimap2.py <graph.gfa[.gz]> <reference path name> <genes.gtf.gz> <transcripts prefix> \
           <minimap2> <work dir> [threads]

Each difference is classed:
  boundary_moved        the HST's own status: a boundary inside a variant, moved
  alignment_clipped     minimap2 left transcript bases unaligned (short end exons); not compared further
  alignment_other_copy  an aligned exon lies over 1 kb outside the HST: taken from another copy
  alignment_exon_count  the alignment has other exons, and the HST is no further from the reference
  junction_choice       a junction placed elsewhere, and the HST is no further from the reference
  repeat_shift          the bases between the two boundary choices repeat beside them (a tenth may differ)
  projection_worse      none of these: the HST is further from the reference than the alignment
Transcripts reaching past the path, and HSTs of a tandem stretch, are not compared.
"""
import collections
import gzip
import os
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from transcripts_common import edit_distances, rc, walks  # noqa: E402

gfa, refname, gtf, prefix, mm2, work = sys.argv[1:7]
threads = int(sys.argv[7]) if len(sys.argv) > 7 else 2
os.makedirs(work, exist_ok=True)
seqs = walks(gfa)

# reference transcripts; the path's interval is half-open (pggb) or inclusive, told apart by its length
m = re.search(r"#([^#]+):(\d+)-(\d+)$", refname)
chrom, start, end = m.group(1).replace("chr", ""), int(m.group(2)), int(m.group(3))
rseq = seqs[refname]
first = start + 1 if len(rseq) == end - start else start   # 1-based position of the path's first base
exons, strand = collections.defaultdict(list), {}
with gzip.open(gtf, "rt") as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "exon" or f[0].replace("chr", "") != chrom:
            continue
        t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
        exons[t].append((int(f[3]), int(f[4])))
        strand[t] = f[6]
reftx = {}
for t, ex in exons.items():
    ex.sort()
    if ex[0][0] >= first and ex[-1][1] < first + len(rseq):
        s = "".join(rseq[a - first:b - first + 1] for a, b in ex)
        reftx[t] = rc(s) if strand[t] == "-" else s
with open(f"{work}/reftx.fa", "w") as out:
    out.writelines(f">{t}\n{s}\n" for t, s in reftx.items())

hst, rows = {}, []
with open(f"{prefix}.hst.fa") as fh:
    for line in fh:
        if line.startswith(">"):
            key = line[1:].split()[0]
        else:
            hst[key] = line.strip()
with open(f"{prefix}.hst.tsv") as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    for line in fh:
        r = dict(zip(hdr, line.rstrip("\n").split("\t")))
        if r["exon_starts"] != "." and r["transcript_id"] in reftx and "tandem:" not in r["detail"]:
            rows.append(r)


def align(i_h):
    i, h = i_h
    fa, paf = f"{work}/h{i}.fa", f"{work}/h{i}.paf"
    if not os.path.exists(paf):
        with open(fa, "w") as out:
            out.write(f">{i}\n{seqs[h]}\n")
        with open(paf + ".tmp", "w") as out:
            subprocess.run([mm2, "-x", "splice:hq", "-c", "--cs=long", "-N", "20", "-p", "0.8", "-t", "1", fa,
                            f"{work}/reftx.fa"], stdout=out, stderr=subprocess.DEVNULL, check=True)
        os.rename(paf + ".tmp", paf)
        os.remove(fa)
    return h, paf


def from_cs(cs, ts):
    """The haplotype's bases of the aligned exons (target strand) and the exons' intervals."""
    out, blocks, pos, bs = [], [], ts, ts
    for o, v in re.findall(r"([=:*+\-~])([A-Za-z0-9]+)", cs):
        if o in "=-":
            out.append(v.upper())
            pos += len(v)
        elif o == "*":
            out.append(v[0].upper())
            pos += 1
        elif o == "~":
            blocks.append((bs, pos))
            pos += int(re.match(r"[a-z]{2}(\d+)", v).group(1))
            bs = pos
    blocks.append((bs, pos))
    return "".join(out), blocks


aln = collections.defaultdict(list)
with ThreadPoolExecutor(threads) as ex:
    for h, paf in ex.map(align, enumerate(sorted({r["haplotype"] for r in rows}))):
        for line in open(paf):
            f = line.rstrip("\n").split("\t")
            s, blocks = from_cs(next(x[5:] for x in f[12:] if x.startswith("cs:Z:")), int(f[7]))
            aln[(f[0], h)].append(dict(qlen=int(f[1]), qs=int(f[2]), qe=int(f[3]), ts=int(f[7]), te=int(f[8]),
                                       seq=rc(s) if f[4] == "-" else s, blocks=blocks))

cls, per_gene, examples, pending = collections.Counter(), collections.defaultdict(collections.Counter), collections.defaultdict(list), []


def note(r, c, *extra):
    cls[c] += 1
    per_gene[r["gene_name"]][c] += 1
    if c != "identical" and len(examples[c]) < 8:
        examples[c].append((f"{r['transcript_id']}|{r['gene_name']}|{r['haplotype']}|{r['copy']}", r["status"], *extra))


for r in rows:
    ours = hst[f"{r['transcript_id']}|{r['gene_name']}|{r['haplotype']}|{r['copy']}"]
    lo, hi = int(r["path_start"]), int(r["path_end"])
    cands = [a for a in aln.get((r["transcript_id"], r["haplotype"]), []) if a["te"] > lo and a["ts"] < hi]
    if not cands:
        note(r, "no_alignment")
        continue
    a = max(cands, key=lambda a: min(a["te"], hi) - max(a["ts"], lo))
    mine = sorted((int(x), int(x) + int(y)) for x, y in zip(r["exon_starts"].split(","), r["exon_lengths"].split(",")))
    if a["seq"] == ours:
        note(r, "identical")
    elif r["status"] == "boundary_moved":
        note(r, "boundary_moved")
    elif a["qs"] > 0 or a["qe"] < a["qlen"]:
        note(r, "alignment_clipped")
    else:
        pending.append((r, ours, a, mine))

d = edit_distances([(f"o{k}", reftx[r["transcript_id"]], ours) for k, (r, ours, a, mine) in enumerate(pending)] +
                   [(f"a{k}", reftx[r["transcript_id"]], a["seq"]) for k, (r, ours, a, mine) in enumerate(pending)])
for k, (r, ours, a, mine) in enumerate(pending):
    eo, ea = d[f"o{k}"], d[f"a{k}"]
    lo, hi = int(r["path_start"]), int(r["path_end"])
    if any(y1 < lo - 1000 or y0 > hi + 1000 for y0, y1 in a["blocks"]):
        c = "alignment_other_copy"
    elif len(a["blocks"]) != len(mine):
        c = "alignment_exon_count" if eo <= ea else "projection_worse"
    elif eo <= ea:
        c = "junction_choice"
    else:
        hs = seqs[r["haplotype"]]

        def periodic(u, v):
            u, v = min(u, v), max(u, v)
            n = v - u
            near = lambda x, y: len(x) == len(y) and sum(p != q for p, q in zip(x, y)) <= max(1, n // 10)
            return u == v or near(hs[u:v], hs[u - n:u]) or near(hs[u:v], hs[v:v + n])

        c = "repeat_shift" if all(periodic(x0, y0) and periodic(x1, y1)
                                  for (x0, x1), (y0, y1) in zip(mine, a["blocks"])) else "projection_worse"
    note(r, c, f"edits: ours {eo}, alignment {ea}")

total = sum(cls.values())
print(f"HSTs compared: {total}")
for c, n in cls.most_common():
    print(f"  {c}\t{n}\t{100 * n / total:.2f}%")
print("per gene:")
for g in sorted(per_gene):
    print(f"  {g}\t" + " ".join(f"{c}={n}" for c, n in per_gene[g].most_common()))
for c, ex in examples.items():
    print(f"examples {c}:")
    for e in ex:
        print("  ", *e)
