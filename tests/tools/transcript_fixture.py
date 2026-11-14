#!/usr/bin/env python3
"""Synthetic loci for `panvar transcripts`, each a reference with genes and haplotypes carrying variants
relative to it. The expected transcripts are spliced from each haplotype's own sequence, independently of
panvar.

  transcript_fixture.py cases                          the case names
  transcript_fixture.py write <case> <dir>             writes locus.gfa, genes.gtf, expected.tsv
  transcript_fixture.py check <dir> <prefix>           compares <prefix>.hst.tsv/.hst.fa with expected.tsv

Variants on the reference (0-based, half-open): ("snp", p), ("del", (a, b)) and ("dup", (a, b)), a stretch
walked twice over the same nodes. A haplotype in "reversed" is another one written from its end.
"""
import random
import sys

COMP = str.maketrans("ACGT", "TGCA")
OFFSET = 1000                       # every reference path is chr1:1000-<1000+L>, half-open


def rc(s):
    return s.translate(COMP)[::-1]


def alt(base):
    return {"A": "C", "C": "G", "G": "T", "T": "A"}[base]


# gene -> (strand, {transcript: exons}); exons on the reference path, possibly reaching past its ends
CASES = {
    # an exon SNP, an exon deletion, a splice-site SNP, a deletion across an exon boundary, a duplicated and
    # a deleted gene, a minus-strand gene, walks written reversed
    "main": dict(
        L=12000,
        genes={"GA": ("+", {"TA1": [(1000, 1300), (2000, 2060), (3000, 3300)], "TA2": [(1000, 1300), (3000, 3300)]}),
               "GB": ("-", {"TB1": [(5000, 5100), (5600, 5800)]}),
               "GC": ("+", {"TC1": [(8000, 8200), (8500, 8700)]})},
        haps={"h1": [("snp", 1100), ("snp", 5050)], "h2": [("del", (1900, 2250))], "h3": [("snp", 1300)],
              "h4": [("del", (2995, 3003))], "h5": [("dup", (7900, 8800))], "h6": [("del", (7900, 8800))]},
        reversed={"h7": "h1", "h8": "h5"}),
    # two genes apart inside one long node
    "bignode": dict(
        L=3000,
        genes={"GE": ("+", {"TE1": [(1000, 1050), (1100, 1150)]}), "GF": ("-", {"TF1": [(1700, 1750), (1800, 1850)]})},
        haps={"h1": [("snp", 2500)]},
        reversed={"h2": "h1"}),
    # a short gene duplicated between long flanks
    "flank": dict(
        L=20200,
        genes={"GG": ("+", {"TG1": [(10020, 10080)]})},
        haps={"h1": [("dup", (10000, 10100))]},
        reversed={"h2": "h1"}),
}


def reference(case):
    rnd = random.Random(7)
    s = [rnd.choice("ACGT") for _ in range(case["L"])]
    for strand, txs in case["genes"].values():   # canonical splice sites on every intron inside the path
        for exons in txs.values():
            for (a, b), (c, d) in zip(exons, exons[1:]):
                donor, acceptor = ("GT", "AG") if strand == "+" else ("CT", "AC")
                if 0 <= b and b + 2 <= case["L"]:
                    s[b:b + 2] = donor
                if 0 <= c - 2 and c <= case["L"]:
                    s[c - 2:c] = acceptor
    return "".join(s)


def build(case):
    L = case["L"]
    ref = reference(case)
    cuts = {0, L}
    for vs in case["haps"].values():
        for kind, v in vs:
            cuts.update((v, v + 1) if kind == "snp" else v)
    cuts = sorted(cuts)
    nodes, ref_nodes = {}, []
    for a, b in zip(cuts, cuts[1:]):
        nid = str(len(nodes) + 1)
        nodes[nid] = ref[a:b]
        ref_nodes.append((nid, a, b))
    alt_node = {}
    for vs in case["haps"].values():
        for kind, v in vs:
            if kind == "snp" and v not in alt_node:
                nid = str(len(nodes) + 1)
                nodes[nid] = alt(ref[v])
                alt_node[v] = nid
    walks, seqs, maps = {"ref": [(n, "+") for n, _, _ in ref_nodes]}, {"ref": ref}, {"ref": [dict((p, p) for p in range(L))]}
    for h, variants in case["haps"].items():
        steps, pos, m, second = [], 0, {}, None
        for nid, a, b in ref_nodes:
            take = nid
            for kind, v in variants:
                if kind == "snp" and a == v:
                    take = alt_node[v]
                if kind == "del" and v[0] <= a and b <= v[1]:
                    take = None
            if take:
                steps.append((take, "+"))
                for p in range(a, b):
                    m[p] = pos + p - a
                pos += b - a
            for kind, v in variants:   # a duplication walks its nodes again after its last one
                if kind == "dup" and b == v[1]:
                    dup = [(n, "+") for n, x, y in ref_nodes if v[0] <= x and y <= v[1]]
                    steps += dup
                    second = dict((p, pos + p - v[0]) for p in range(*v))
                    pos += v[1] - v[0]
        walks[h] = steps
        seqs[h] = "".join(nodes[n] for n, _ in steps)
        maps[h] = [m] + ([second] if second else [])
    for h, src in case.get("reversed", {}).items():
        walks[h] = [(n, "-" if o == "+" else "+") for n, o in reversed(walks[src])]
        seqs[h] = rc(seqs[src])
        maps[h] = maps[src]
    return ref, nodes, walks, seqs, maps


def project(case, ref, hseq, m, strand, exons):
    """One transcript through one copy's map: (status, notes, sequence)."""
    L = case["L"]
    n_all = len(exons)
    inside = [(max(a, 0), min(b, L)) for a, b in exons if b > 0 and a < L]
    num = {k: (k + 1 if strand == "+" else n_all - k) for k in range(n_all)}
    first_in = next(k for k, (a, b) in enumerate(exons) if b > 0 and a < L)
    notes, missing, moved = [], False, False

    def cut(c, start, e):
        """The haplotype position of reference cut c (between bases c - 1 and c)."""
        nonlocal moved
        if c in m:
            return m[c]
        if c - 1 in m:
            return m[c - 1] + 1
        xl = max(p for p in m if p < c) + 1   # inside a deletion: the nearer shared position
        xr = min(p for p in m if p > c)
        dl, dr = c - xl, xr - c
        left = dl < dr or (dl == dr and start)   # a tie keeps the variant inside the exon
        moved = True
        shift = (dl if start else -dl) if left else (-dr if start else dr)
        notes.append(f"e{e}{'.start:' if start == (strand == '+') else '.end:'}{shift:+d}")
        return m[xl - 1] + 1 if left else m[xr]

    hx = []
    for k, (a, b) in enumerate(inside):
        e = num[first_in + k]
        if not any(p in m for p in range(a, b)):
            missing = True
            notes.append(f"e{e}:missing")
            continue
        hx.append((cut(a, True, e), cut(b, False, e)))
    if missing:
        return "missing_exon", notes, ""
    lost = changed = False
    for k in range(len(inside) - 1):
        intron = (n_all - (first_in + k) - 1) if strand == "-" else first_in + k + 1
        hl, hr = hx[k][1], hx[k + 1][0]
        if hr < hl + 4:
            lost = True
            notes.append(f"i{intron}:intron_lost")
            continue
        rl, rr = inside[k][1], inside[k + 1][0]
        for site, r, h in ((("donor", "acceptor")[strand == "-"], ref[rl:rl + 2], hseq[hl:hl + 2]),
                           (("acceptor", "donor")[strand == "-"], ref[rr - 2:rr], hseq[hr - 2:hr])):
            if r != h:
                changed = True
                notes.append(f"i{intron}.{site}:{rc(r) if strand == '-' else r}>{rc(h) if strand == '-' else h}")
    clipped = len(inside) < n_all or inside[0][0] != exons[first_in][0] or inside[-1][1] != exons[first_in + len(inside) - 1][1]
    s = "".join(hseq[a:b] for a, b in hx)
    status = ("intron_lost" if lost else "clipped_at_edge" if clipped else "boundary_moved" if moved
              else "splice_site_changed" if changed else "ok")
    return status, notes, rc(s) if strand == "-" else s


def expected(case, ref, seqs, maps):
    rows = []
    for h in ["ref", *case["haps"], *case.get("reversed", {})]:
        for g, (strand, txs) in case["genes"].items():
            lo = min(max(a, 0) for ex in txs.values() for a, b in ex if b > 0 and a < case["L"])
            hi = max(min(b, case["L"]) for ex in txs.values() for a, b in ex if b > 0 and a < case["L"])
            # copies: the maps holding the gene's bases
            copies = [m for m in maps[h] if any(p in m for p in range(lo, hi))]
            for t, exons in txs.items():
                if not any(b > 0 and a < case["L"] for a, b in exons):
                    continue
                if not copies:
                    rows.append((t, g, h, 0, "no_copy", [], ""))
                for c, m in enumerate(copies, 1):
                    st, notes, s = project(case, ref, seqs[case.get("reversed", {}).get(h, h)], m, strand, exons)
                    rows.append((t, g, h, c, st, notes, s))
    return rows


def name(h, L):
    return f"ref#0#chr1:{OFFSET}-{OFFSET + L}" if h == "ref" else f"{h}#1#ctg{h[1:]}:0-{L}"


def write(case_name, d):
    case = CASES[case_name]
    ref, nodes, walks, seqs, maps = build(case)
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
            out.write(f"P\t{name(h, case['L'])}\t{','.join(n + o for n, o in steps)}\t*\n")
    with open(f"{d}/genes.gtf", "w") as out:
        for g, (strand, txs) in case["genes"].items():
            lo = min(a for ex in txs.values() for a, _ in ex) + OFFSET + 1
            hi = max(b for ex in txs.values() for _, b in ex) + OFFSET
            attr = f'gene_id "{g}"; gene_name "{g}"; gene_biotype "protein_coding";'
            out.write(f"chr1\ttest\tgene\t{lo}\t{hi}\t.\t{strand}\t.\t{attr}\n")
            for k, (t, exons) in enumerate(txs.items()):
                ta = attr + f' transcript_id "{t}"; transcript_name "{t}"; transcript_biotype "protein_coding";'
                if k == 0:
                    ta += ' tag "Ensembl_canonical";'
                out.write(f"chr1\ttest\ttranscript\t{exons[0][0] + OFFSET + 1}\t{exons[-1][1] + OFFSET}\t.\t{strand}\t.\t{ta}\n")
                for a, b in exons:
                    out.write(f"chr1\ttest\texon\t{a + OFFSET + 1}\t{b + OFFSET}\t.\t{strand}\t.\t{ta}\n")
    with open(f"{d}/expected.tsv", "w") as out:
        out.write("transcript_id\tgene_name\thaplotype\tcopy\tstatus\tdetail\tsequence\n")
        for t, g, h, c, st, notes, s in expected(case, ref, seqs, maps):
            out.write(f"{t}\t{g}\t{name(h, case['L'])}\t{c}\t{st}\t{';'.join(notes) or '.'}\t{s or '.'}\n")


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
            k = (r["transcript_id"], r["haplotype"], r["copy"])
            key = f"{r['transcript_id']}|{r['gene_name']}|{r['haplotype']}|{r['copy']}"
            got[k] = (r, fa.get(key) or ".")
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
    if sys.argv[1] == "cases":
        print(" ".join(CASES))
    elif sys.argv[1] == "write":
        write(sys.argv[2], sys.argv[3])
    else:
        sys.exit(check(sys.argv[2], sys.argv[3]))
