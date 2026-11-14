#!/usr/bin/env python3
"""A reference transcriptome for `panvar background --add --transcriptome`: every transcript of a GTF, spliced
from the reference FASTA, named `transcript_id|gene_id|gene_name` so recruit can tell a locus's own genes.

    reference_transcripts.py <reference.fa[.gz]> <genes.gtf[.gz]> <out.fa>

Contigs are matched with or without a `chr` prefix; transcripts on contigs the FASTA lacks are left out and
counted. GENCODE's transcripts FASTA and `gffread -w` output (headers with `gene=<name>`) work as they are.
"""
import collections
import gzip
import re
import sys

fa, gtf, out = sys.argv[1:4]
opener = lambda p: gzip.open(p, "rt") if p.endswith(".gz") else open(p)
COMP = str.maketrans("ACGTacgt", "TGCAtgca")
bare = lambda c: c[3:] if c.startswith("chr") else c

exons, meta = collections.defaultdict(list), {}
with opener(gtf) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "exon":
            continue
        attr = dict(re.findall(r'(\S+) "([^"]*)"', f[8]))
        t = attr["transcript_id"]
        exons[t].append((int(f[3]), int(f[4])))
        meta[t] = (bare(f[0]), f[6], attr.get("gene_id", "."), attr.get("gene_name", attr.get("gene_id", ".")))
wanted = {m[0] for m in meta.values()}
by_contig = collections.defaultdict(list)
for t, m in meta.items():
    by_contig[m[0]].append(t)

written, seen = 0, set()
with opener(fa) as fh, open(out, "w") as o:
    def flush(name, chunks):
        global written
        if name is None or bare(name) not in wanted:
            return
        seq = "".join(chunks)
        seen.add(bare(name))
        for t in by_contig[bare(name)]:
            _, strand, gid, gname = meta[t]
            s = "".join(seq[a - 1:b] for a, b in sorted(exons[t])).upper()
            if strand == "-":
                s = s.translate(COMP)[::-1]
            o.write(f">{t}|{gid}|{gname}\n{s}\n")
            written += 1
    name, chunks = None, []
    for line in fh:
        if line.startswith(">"):
            flush(name, chunks)
            name, chunks = line[1:].split()[0], []
        elif bare(name) in wanted:
            chunks.append(line.strip())
    flush(name, chunks)
missing = sum(len(by_contig[c]) for c in wanted - seen)
print(f"transcripts written: {written}; on contigs the FASTA lacks: {missing}", file=sys.stderr)
