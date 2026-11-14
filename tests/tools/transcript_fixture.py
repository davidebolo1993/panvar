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
    # an intron deleted whole, one cut to its splice sites, one lost with the exon ends around it
    "intron": dict(
        L=5000,
        genes={"GI": ("+", {"TI1": [(1000, 1200), (2000, 2200), (3000, 3200)]}),
               "GJ": ("-", {"TJ1": [(4000, 4100), (4300, 4400)]})},
        haps={"h1": [("del", (1200, 2000)), ("del", (4100, 4300))], "h2": [("del", (1202, 1998))],
              "h3": [("del", (2150, 3050))]},
        reversed={"h4": "h1"}),
    # transcripts reaching past either end of the reference path, and one with no exon on it
    "edge": dict(
        L=5000,
        genes={"GD": ("+", {"TD1": [(4500, 4700), (5100, 5300)]}), "GK": ("+", {"TK1": [(3800, 3900), (4900, 5100)]}),
               "GH": ("-", {"TH1": [(-300, -100), (200, 400)]}), "GL": ("+", {"TL1": [(-500, -400), (5400, 5500)]})},
        haps={"h1": [("snp", 4550)]},
        reversed={"h2": "h1"}),
    # exons inside a repeat unit the reference walks three times and haplotypes 1, 3 or 5 times
    "tandem": dict(units={"ref": 3, "h1": 1, "h2": 5, "h3": 3}, reversed={"h4": "h2"}),
    # reference units whose exons differ in length (100, 150, 100 bases from one start)
    "tandem_layout": dict(units={"ref": 3, "h1": 3, "h2": 1, "h3": 5}, layout=[400, 450, 400], reversed={"h4": "h1"}),
    # the main case with haplotypes' own annotations (write_annotated)
    "annotated": dict(),
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
    before = sum(min(b, 0) - a for a, b in exons if a < 0)   # bases past the path's start and end
    after = sum(b - max(a, L) for a, b in exons if b > L)
    clipped = before + after > 0
    for end, bp in (("5'", before if strand == "+" else after), ("3'", after if strand == "+" else before)):
        if bp:
            notes.append(f"clipped:{end}:{bp}")
    s = "".join(hseq[a:b] for a, b in hx)
    status = ("intron_lost" if lost else "clipped_at_edge" if clipped else "boundary_moved" if moved
              else "splice_site_changed" if changed else "ok")
    return status, notes, rc(s) if strand == "-" else s


def expected(case, ref, seqs, maps):
    rows = []
    for h in ["ref", *case["haps"], *case.get("reversed", {})]:
        for g, (strand, txs) in case["genes"].items():
            if not any(b > 0 and a < case["L"] for ex in txs.values() for a, b in ex):
                continue   # no exon on the path: not projected
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


def write_tandem(case, d):
    """Flank A (2 kb), a 1 kb unit walked n times over the same nodes, flank B (2 kb). Gene GT has an exon in
    A, one in each unit and one in B, so each haplotype's transcript has n unit exons. The reference's unit
    exons start at one node and end at base 400 or 450 of the unit (`layout`, per reference unit); a
    haplotype's unit j takes the layout of reference unit j when the counts agree, else of unit j * k / n.
    `insert` puts 150 bases inside the exon of one unit of a haplotype."""
    rnd = random.Random(11)
    seq = lambda n: "".join(rnd.choice("ACGT") for _ in range(n))
    A, U, B, I = list(seq(2000)), list(seq(1000)), list(seq(2000)), seq(150)
    A[1600:1602], U[298:300], U[400:402], U[450:452], B[498:500] = "GT", "AG", "GT", "GT", "AG"
    A, U, B = "".join(A), "".join(U), "".join(B)
    nodes = {"1": A, "2": U[:300], "3": U[300:350], "4": U[350:400], "5": U[400:450], "6": U[450:], "7": B, "8": I}
    units, k = case["units"], case["units"]["ref"]
    layout = case.get("layout", [400] * k)
    insert = case.get("insert", {})

    def walk(h, n):
        w = ["1"]
        for j in range(n):
            w += ["2", "3", "8", "4", "5", "6"] if insert.get(h) == j else ["2", "3", "4", "5", "6"]
        return w + ["7"]

    walks = {h: [(x, "+") for x in walk(h, n)] for h, n in units.items()}
    for h, src in case.get("reversed", {}).items():
        walks[h] = [(x, "-") for x in reversed(walk(src, units[src]))]
    L = 2000 + 1000 * k + 2000
    exons = ([(1500, 1600)] + [(2000 + 1000 * u + 300, 2000 + 1000 * u + layout[u]) for u in range(k)] +
             [(2000 + 1000 * k + 500, 2000 + 1000 * k + 600)])
    with open(f"{d}/locus.gfa", "w") as out:
        out.write("H\tVN:Z:1.0\n")
        for n, s in nodes.items():
            out.write(f"S\t{n}\t{s}\n")
        for e in ("1 2", "2 3", "3 4", "3 8", "8 4", "4 5", "5 6", "6 2", "6 7"):
            a, b = e.split()
            out.write(f"L\t{a}\t+\t{b}\t+\t0M\n")
        for h, steps in walks.items():
            out.write(f"P\t{name(h, L)}\t{','.join(n + o for n, o in steps)}\t*\n")
    attr = 'gene_id "GT"; gene_name "GT"; gene_biotype "protein_coding"; transcript_id "TT1"; transcript_name "TT1"; tag "Ensembl_canonical";'
    with open(f"{d}/genes.gtf", "w") as out:
        out.write(f"chr1\ttest\ttranscript\t{exons[0][0] + OFFSET + 1}\t{exons[-1][1] + OFFSET}\t.\t+\t.\t{attr}\n")
        for a, b in exons:
            out.write(f"chr1\ttest\texon\t{a + OFFSET + 1}\t{b + OFFSET}\t.\t+\t.\t{attr}\n")
    with open(f"{d}/expected.tsv", "w") as out:
        out.write("transcript_id\tgene_name\thaplotype\tcopy\tstatus\tdetail\tsequence\n")
        for h in walks:
            src = case.get("reversed", {}).get(h, h)
            n = units[src]
            notes = [f"tandem:{k}>{n}"] if n != k else []
            body = ""
            for j in range(n):
                end = layout[j if n == k else j * k // n]
                body += U[300:350] + (I if insert.get(src) == j else "") + U[350:end]
            if src in insert:
                notes.append("tandem_resized:1")
            out.write(f"TT1\tGT\t{name(h, L)}\t1\tok\t{';'.join(notes) or '.'}\t{A[1500:1600] + body + B[500:600]}\n")


def write_annotated(d):
    """The main case, with annotations of h1 (GFF3: the reference transcript TA1, a novel isoform, and a
    second copy of GB the projection does not see), h2 (GFF3: a transcript past the path's end), h3 (GFF3 on
    an aliased contig name) and h5 (GTF: both copies of GC, as two gene records)."""
    write("main", d)
    case = CASES["main"]
    L = case["L"]
    ref, nodes, walks, seqs, maps = build(case)
    m1 = maps["h1"][0]
    pos = lambda m, p: m[p] + 1   # 1-based on the haplotype's contig (its path starts at 0)
    gff = lambda ctg, kind, a, b, strand, attrs: f"{ctg}\tCAT\t{kind}\t{a}\t{b}\t.\t{strand}\t.\t{attrs}\n"
    rows = []   # expected annotated-only rows: transcript, gene, haplotype, copy, status, detail, sequence
    h1 = [gff("ctg1", "gene", pos(m1, 1000), pos(m1, 3299), "+", "ID=g1;gene_name=GA"),
          gff("ctg1", "transcript", pos(m1, 1000), pos(m1, 3299), "+", "ID=TA1;Parent=g1;transcript_id=TA1")]
    h1 += [gff("ctg1", "exon", pos(m1, a), pos(m1, b - 1), "+", "Parent=TA1") for a, b in case["genes"]["GA"][1]["TA1"]]
    novel = [(1000, 1300), (3000, 3150)]   # a novel isoform: exon 3 ends early
    h1 += [gff("ctg1", "transcript", pos(m1, 1000), pos(m1, 3149), "+", "ID=nov1;Parent=g1;transcript_class=putative_novel_isoform")]
    h1 += [gff("ctg1", "exon", pos(m1, a), pos(m1, b - 1), "+", "Parent=nov1") for a, b in novel]
    rows.append(("nov1", "GA", "h1", 1, "ok", "transcript_class:putative_novel_isoform", "".join(seqs["h1"][m1[a]:m1[b - 1] + 1] for a, b in novel)))
    gb2 = [(6500, 6600), (6800, 6900)]       # a second GB copy, minus strand, where the projection sees none
    h1 += [gff("ctg1", "gene", pos(m1, 5000), pos(m1, 5799), "-", "ID=g2;gene_name=GB"),
           gff("ctg1", "transcript", pos(m1, 5000), pos(m1, 5799), "-", "ID=TB1;Parent=g2;transcript_id=TB1")]
    h1 += [gff("ctg1", "exon", pos(m1, a), pos(m1, b - 1), "-", "Parent=TB1") for a, b in case["genes"]["GB"][1]["TB1"]]
    h1 += [gff("ctg1", "gene", pos(m1, 6500), pos(m1, 6899), "-", "ID=g3;gene_name=GB"),
           gff("ctg1", "transcript", pos(m1, 6500), pos(m1, 6899), "-", "ID=TB1x;Parent=g3;transcript_id=TB1")]
    h1 += [gff("ctg1", "exon", pos(m1, a), pos(m1, b - 1), "-", "Parent=TB1x") for a, b in gb2]
    rows.append(("TB1", "GB", "h1", 2, "ok", ".", rc("".join(seqs["h1"][m1[a]:m1[b - 1] + 1] for a, b in gb2))))
    m2 = maps["h2"][0]
    edge = [(m2[11500], m2[11700 - 1] + 1), (len(seqs["h2"]) + 100, len(seqs["h2"]) + 200)]   # path coordinates
    h2 = [gff("ctg2", "transcript", edge[0][0] + 1, edge[1][1], "+", "ID=edge1;gene_name=GZ;transcript_id=edge1")]
    h2 += [gff("ctg2", "exon", a + 1, b, "+", "Parent=edge1") for a, b in edge]
    rows.append(("edge1", "GZ", "h2", 1, "clipped_at_edge", "clipped:3':100", seqs["h2"][edge[0][0]:edge[0][1]]))
    m3 = maps["h3"][0]
    alias_ex = [(4000, 4200)]
    h3 = [gff("CTG3ALT", "transcript", pos(m3, 4000), pos(m3, 4199), "+", "ID=al1;gene_name=GY;transcript_id=al1")]
    h3 += [gff("CTG3ALT", "exon", pos(m3, a), pos(m3, b - 1), "+", "Parent=al1") for a, b in alias_ex]
    rows.append(("al1", "GY", "h3", 1, "ok", ".", seqs["h3"][m3[4000]:m3[4199] + 1]))
    h5 = []
    for k, m in enumerate(maps["h5"], 1):   # GTF: both GC copies, as projected
        attr = f'gene_id "gc{k}"; gene_name "GC"; transcript_id "TC1";'
        h5 += [f'ctg5\tCAT\texon\t{pos(m, a)}\t{pos(m, b - 1)}\t.\t+\t.\t{attr}\n' for a, b in case["genes"]["GC"][1]["TC1"]]
    for h, lines, ext in (("h1", h1, "gff3"), ("h2", h2, "gff3"), ("h3", h3, "gff3"), ("h5", h5, "gtf")):
        with open(f"{d}/{h}.{ext}", "w") as out:
            out.write("##gff-version 3\n" if ext == "gff3" else "")
            out.writelines(lines)
    with open(f"{d}/hap_annotations.tsv", "w") as out:
        out.write("haplotype\tpath\n")
        for h, ext in (("h1", "gff3"), ("h2", "gff3"), ("h3", "gff3"), ("h5", "gtf")):
            out.write(f"{h}#1\t{d}/{h}.{ext}\n")
    with open(f"{d}/alias.tsv", "w") as out:
        out.write("# names of one contig per row\nctg3\tCTG3ALT\n")
    # TA1 on h1 and both TC1 copies on h5 are the projected records, of both sources
    lines = open(f"{d}/expected.tsv").read().splitlines()
    with open(f"{d}/expected.tsv", "w") as out:
        out.write(lines[0] + "\tsource\n")
        for line in lines[1:]:
            f = line.split("\t")
            both = (f[0], f[2].split("#")[0]) in (("TA1", "h1"), ("TC1", "h5"), ("TB1", "h1"))
            out.write(line + ("\tprojected+annotated" if both else "\tprojected") + "\n")
        for t, g, h, c, st, det, sq in rows:
            out.write(f"{t}\t{g}\t{name(h, L)}\t{c}\t{st}\t{det}\t{sq}\tannotated\n")


def write(case_name, d):
    case = CASES[case_name]
    if case_name.startswith("tandem"):
        return write_tandem(case, d)
    if case_name == "annotated":
        return write_annotated(d)
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
            if not any(b > 0 and a < case["L"] for ex in txs.values() for a, b in ex):
                continue   # no exon on the path: not projected
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
            t, g, h, c, st, det, s, *src = line.rstrip("\n").split("\t")
            exp[(t, h, c, src[0] if src else "projected")] = (st, det, s)
    got, fa = {}, {}
    with open(f"{prefix}.hst.fa") as fh:
        for line in fh:
            if line.startswith(">"):
                key = line[1:].split()[0]
                if key in fa:
                    sys.exit(f"duplicate FASTA name: {key}")
                fa[key] = None
            else:
                fa[key] = line.strip()
    with open(f"{prefix}.hst.tsv") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(hdr, line.rstrip("\n").split("\t")))
            k = (r["transcript_id"], r["haplotype"], r["copy"], r["source"])
            if k in got:
                sys.exit(f"duplicate row: {k}")
            key = f"{r['transcript_id']}|{r['gene_name']}|{r['haplotype']}|{r['copy']}"
            key += "|annotated" if r["source"] == "annotated" else ""
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
    return 1 if bad or not exp else 0


if __name__ == "__main__":
    if sys.argv[1] == "cases":
        print(" ".join(CASES))
    elif sys.argv[1] == "write":
        write(sys.argv[2], sys.argv[3])
    else:
        sys.exit(check(sys.argv[2], sys.argv[3]))
