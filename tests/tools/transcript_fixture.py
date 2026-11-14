#!/usr/bin/env python3
"""Synthetic locus for `panvar transcripts`: a reference with three genes and haplotypes carrying an exon
SNP, an exon deletion, a splice-site SNP, a deletion across an exon boundary, a duplicated gene (walked
twice over the same nodes), a deleted gene, and a walk written in the reverse orientation. The expected
transcripts are spliced from each haplotype's own sequence, independently of panvar.

  transcript_fixture.py write <dir>              writes locus.gfa, genes.gtf, expected.tsv
  transcript_fixture.py check <dir> <prefix>     compares <prefix>.hst.tsv/.hst.fa with expected.tsv
"""
import random
import sys

COMP = str.maketrans("ACGT", "TGCA")


def rc(s):
    return s.translate(COMP)[::-1]


L = 12000
OFFSET = 1000                       # the reference path is chr1:1000-13000, half-open
GENES = {                           # gene -> strand, transcripts -> exons (0-based, half-open on the path)
    "GA": ("+", {"TA1": [(1000, 1300), (2000, 2060), (3000, 3300)], "TA2": [(1000, 1300), (3000, 3300)]}),
    "GB": ("-", {"TB1": [(5000, 5100), (5600, 5800)]}),
    "GC": ("+", {"TC1": [(8000, 8200), (8500, 8700)]}),
}
MODULE = (7900, 8800)               # the stretch holding GC, duplicated or deleted
SNP_A, SNP_B, SPLICE = 1100, 5050, 1300
DEL_EXON = (1900, 2250)             # exon 2 of GA
DEL_EDGE = (2995, 3003)             # across the start of GA's exon 3
# haplotype -> variants
HAPS = {
    "h1": [("snp", SNP_A), ("snp", SNP_B)],
    "h2": [("del", DEL_EXON)],
    "h3": [("snp", SPLICE)],
    "h4": [("del", DEL_EDGE)],
    "h5": [("dup", MODULE)],
    "h6": [("del", MODULE)],
}
REVERSED = {"h7": "h1", "h8": "h5"}  # written from the end: h7 is h1, h8 is h5


def reference():
    rnd = random.Random(7)
    s = [rnd.choice("ACGT") for _ in range(L)]
    for _, (strand, txs) in GENES.items():   # canonical splice sites on every intron
        for exons in txs.values():
            for (a, b), (c, d) in zip(exons, exons[1:]):
                donor, acceptor = ("GT", "AG") if strand == "+" else ("CT", "AC")
                s[b:b + 2] = donor
                s[c - 2:c] = acceptor
    return "".join(s)


def alt(base):
    return {"A": "C", "C": "G", "G": "T", "T": "A"}[base]


def build():
    ref = reference()
    cuts = {0, L, *MODULE, SNP_A, SNP_A + 1, SNP_B, SNP_B + 1, SPLICE, SPLICE + 1, *DEL_EXON, *DEL_EDGE}
    cuts = sorted(cuts)
    nodes, ref_nodes = {}, []
    for a, b in zip(cuts, cuts[1:]):
        nid = str(len(nodes) + 1)
        nodes[nid] = ref[a:b]
        ref_nodes.append((nid, a, b))
    alt_node = {}
    for p in (SNP_A, SNP_B, SPLICE):
        nid = str(len(nodes) + 1)
        nodes[nid] = alt(ref[p])
        alt_node[p] = nid
    walks = {"ref": [(n, "+") for n, _, _ in ref_nodes]}
    seqs = {"ref": ref}
    for h, variants in HAPS.items():
        steps = []
        for nid, a, b in ref_nodes:
            take = nid
            for kind, v in variants:
                if kind == "snp" and a == v:
                    take = alt_node[v]
                if kind == "del" and v[0] <= a and b <= v[1]:
                    take = None
            if take:
                steps.append((take, "+"))
        for kind, v in variants:   # a duplication walks the module's nodes a second time
            if kind == "dup":
                idx = [i for i, (nid, a, b) in enumerate(ref_nodes) if v[0] <= a and b <= v[1]]
                mod = [(ref_nodes[i][0], "+") for i in idx]
                at = [s[0] for s in steps].index(mod[-1][0]) + 1
                steps[at:at] = mod
        walks[h] = steps
        seqs[h] = "".join(nodes[n] for n, _ in steps)
    for h, src in REVERSED.items():
        walks[h] = [(n, "-" if o == "+" else "+") for n, o in reversed(walks[src])]
        seqs[h] = rc(seqs[src])
    return ref, nodes, walks, seqs


def hap_map(variants):
    """Reference position -> haplotype position, per base kept (dup: the copies' offsets)."""
    m, pos = {}, 0
    deleted = [v for k, v in variants if k == "del"]
    dup = [v for k, v in variants if k == "dup"]
    copies = []
    for p in range(L):
        if any(a <= p < b for a, b in deleted):
            continue
        m[p] = pos
        pos += 1
        if dup and p == dup[0][1] - 1:   # end of the module: the second copy follows
            copies.append(dict((q, m[q] + (pos - m[dup[0][0]])) for q in range(*dup[0])))
            pos += dup[0][1] - dup[0][0]
    return m, copies


def expected(ref, seqs):
    rows = []
    for h in ["ref", *HAPS, *REVERSED]:
        src = REVERSED.get(h, h)
        variants = HAPS.get(src, [])
        hseq = seqs[src]
        m, dup_copies = hap_map(variants)
        for g, (strand, txs) in GENES.items():
            maps = [m]
            if g == "GC" and dup_copies:
                maps = [{q: m[q] for q in range(*MODULE)}, dup_copies[0]]
            for t, exons in txs.items():
                for c, mp in enumerate(maps, 1):
                    if g == "GC" and ("del", MODULE) in variants:
                        rows.append((t, g, h, 0, "no_copy", "", ""))
                        continue
                    status, detail, parts = "ok", [], []
                    for k, (a, b) in enumerate(exons):
                        n = k + 1 if strand == "+" else len(exons) - k
                        if not any(p in mp for p in range(a, b)):
                            status = "missing_exon"
                            detail.append(f"e{n}:missing")
                            continue
                        lo = a
                        if ("del", DEL_EDGE) in variants and a == DEL_EDGE[1] - 3:
                            lo = DEL_EDGE[1]   # moved to the nearer shared position, 3 bp into the exon
                            status, detail = "boundary_moved", detail + [f"e{n}.start:-3"]
                        parts.append(hseq[mp[lo]:mp[b - 1] + 1])
                    if status == "missing_exon":
                        rows.append((t, g, h, c, status, ";".join(detail), ""))
                        continue
                    if ("snp", SPLICE) in variants and g == "GA":
                        status, detail = "splice_site_changed", [f"i1.donor:GT>{alt('G')}T"]
                    s = "".join(parts)
                    rows.append((t, g, h, c, status, ";".join(detail), rc(s) if strand == "-" else s))
    return rows


def name(h):
    return "ref#0#chr1:1000-13000" if h == "ref" else f"{h}#1#ctg{h[1:]}:0-{L}"


def write(d):
    ref, nodes, walks, seqs = build()
    edges = set()
    for steps in walks.values():
        for (a, oa), (b, ob) in zip(steps, steps[1:]):
            edges.add((a, oa, b, ob))
    with open(f"{d}/locus.gfa", "w") as out:
        out.write("H\tVN:Z:1.0\n")
        for n, s in nodes.items():
            out.write(f"S\t{n}\t{s}\n")
        for a, oa, b, ob in sorted(edges):
            out.write(f"L\t{a}\t{oa}\t{b}\t{ob}\t0M\n")
        for h, steps in walks.items():
            out.write(f"P\t{name(h)}\t{','.join(n + o for n, o in steps)}\t*\n")
    with open(f"{d}/genes.gtf", "w") as out:
        for g, (strand, txs) in GENES.items():
            lo = min(a for ex in txs.values() for a, _ in ex) + OFFSET + 1
            hi = max(b for ex in txs.values() for _, b in ex) + OFFSET
            attr = f'gene_id "{g}"; gene_name "{g}"; gene_biotype "protein_coding";'
            out.write(f"chr1\ttest\tgene\t{lo}\t{hi}\t.\t{strand}\t.\t{attr}\n")
            for t, exons in txs.items():
                ta = attr + f' transcript_id "{t}"; transcript_name "{t}"; transcript_biotype "protein_coding";'
                if t in ("TA1", "TB1", "TC1"):
                    ta += ' tag "Ensembl_canonical";'
                out.write(f"chr1\ttest\ttranscript\t{exons[0][0] + OFFSET + 1}\t{exons[-1][1] + OFFSET}\t.\t{strand}\t.\t{ta}\n")
                for a, b in exons:
                    out.write(f"chr1\ttest\texon\t{a + OFFSET + 1}\t{b + OFFSET}\t.\t{strand}\t.\t{ta}\n")
    with open(f"{d}/expected.tsv", "w") as out:
        out.write("transcript_id\tgene_name\thaplotype\tcopy\tstatus\tdetail\tsequence\n")
        for t, g, h, c, st, det, s in expected(ref, seqs):
            out.write(f"{t}\t{g}\t{name(h)}\t{c}\t{st}\t{det or '.'}\t{s or '.'}\n")


def check(d, prefix):
    exp = {}
    with open(f"{d}/expected.tsv") as fh:
        fh.readline()
        for line in fh:
            t, g, h, c, st, det, s = line.rstrip("\n").split("\t")
            exp[(t, h, c)] = (st, det, s)
    got, fa = {}, {}
    with open(f"{prefix}.hst.fa") as fh:
        for line in fh:
            if line.startswith(">"):
                key = line[1:].split()[0]
            else:
                fa[key] = line.strip()
    with open(f"{prefix}.hst.tsv") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(hdr, line.rstrip("\n").split("\t")))
            key = f"{r['transcript_id']}|{r['gene_name']}|{r['haplotype']}|{r['copy']}"
            got[(r["transcript_id"], r["haplotype"], r["copy"])] = (r, fa.get(key, "."))
    bad = 0
    for k in sorted(set(exp) | set(got)):
        if k not in got or k not in exp:
            print(f"  row {'missing' if k not in got else 'unexpected'}: {k}")
            bad += 1
            continue
        st, det, s = exp[k]
        r, seq = got[k]
        problems = []
        if r["status"] != st:
            problems.append(f"status {r['status']} != {st}")
        if det != "." and not all(x in r["detail"].split(";") for x in det.split(";")):
            problems.append(f"detail {r['detail']} lacks {det}")
        if seq != s:
            problems.append("sequence differs")
        if problems:
            print(f"  {k}: " + "; ".join(problems))
            bad += 1
    print(f"rows checked: {len(exp)}; wrong: {bad}")
    return 1 if bad else 0


if __name__ == "__main__":
    if sys.argv[1] == "write":
        write(sys.argv[2])
    else:
        sys.exit(check(sys.argv[2], sys.argv[3]))
