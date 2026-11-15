#!/usr/bin/env python3
"""Synthetic loci for `panvar transcripts`, each a reference with genes and haplotypes carrying variants
relative to it. The expected transcripts are spliced from each haplotype's own sequence, independently of
panvar.

  transcript_fixture.py cases                          the case names
  transcript_fixture.py write <case> <dir>             writes locus.gfa, genes.gtf, expected.tsv
  transcript_fixture.py check <dir> <prefix> [<exp>]   compares <prefix>.hst.tsv/.hst.fa with expected.tsv (or <exp>)

Variants on the reference (0-based, half-open): ("snp", p), ("del", (a, b)), ("sub", (a, b)), a stretch of the
same length with every third base changed, or every `sub_step`-th (clustered SNPs folded into one node),
and ("dup", (a, b)), a stretch walked twice over the same nodes. A haplotype in "reversed" is another one written from its end.
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
    "tandem_layout": dict(units={"ref": 3, "h1": 3, "h2": 1, "h3": 5}, layout=[[(300, 400)], [(300, 450)], [(300, 400)]],
                          reversed={"h4": "h1"}),
    # a haplotype with 50 units between 2 kb flanks
    "tandem_long": dict(units={"ref": 3, "h1": 50}),
    # 150 bases inserted inside one unit's exon
    "tandem_insert": dict(units={"ref": 3, "h1": 3, "h2": 4}, insert={"h1": (1, 1, 150), "h2": (3, 1, 150)}),
    # a unit whose exon's first node differs (the first unit of one haplotype, a middle one of another)
    "tandem_nostart": dict(units={"ref": 3, "h1": 3, "h2": 4}, alt={"h1": (0, "b"), "h2": (2, "b")}),
    # a node walked twice inside a unit, the exon ending at its second visit
    "tandem_revisit": dict(units={"ref": 3, "h1": 3, "h2": 2}, pieces=[("S", 50), ("E", 50), ("M", 400), ("T", 450)],
                           order=["S", "E", "M", "E", "T"], layout=[[(0, 550)]] * 3),
    # a reference unit without exons between two with
    "tandem_empty": dict(units={"ref": 3, "h1": 3, "h2": 4, "h3": 2}, layout=[[(300, 400)], [], [(300, 400)]]),
    # an exon that passes its own first node again: the starts there are no units, and the exons are left to
    # the nodes
    "tandem_overrun": dict(units={"ref": 3, "h1": 3}, pieces=[("E", 50), ("M", 400), ("T", 500)],
                           order=["E", "M", "E", "T"], layout=[[(0, 500)]] * 3, reversed={"h2": "h1"}),
    # whole genes duplicated in tandem over several nodes, the first copy with a deletion or an SNP, the second
    # the reference's: a chain moving from one copy to the other inside an exon would join them
    "dup_gene": dict(
        L=9000,
        genes={"G0": ("+", {"T0a": [(600, 842), (949, 1230), (1616, 1838), (2100, 2370)], "T0b": [(600, 842), (1616, 1838), (2100, 2370)]}),
               "G1": ("+", {"T1a": [(3195, 3386), (3855, 4046)]}),
               "G2": ("-", {"T2a": [(4778, 4982), (5189, 5451)]})},
        haps={"h1": [("snp", 3184)], "h2": [("dup", (4710, 5562)), ("del", (4052, 4055))],
              "h3": [("dup", (519, 2420)), ("del", (1840, 1843)), ("del", (5704, 5712))],
              "h4": [("dup", (546, 2445)), ("snp", 3859), ("snp", 5439), ("del", (1214, 1219))], "h5": [("snp", 4973)]},
        reversed={"h6": "h1", "h7": "h3"}),
    # the same, the first copy's deletion across the end of the gene's last exon
    "dup_gene_end": dict(
        L=9000,
        genes={"G0": ("-", {"T0a": [(600, 813), (1297, 1593)]}), "G1": ("+", {"T1a": [(2595, 2813), (3177, 3346)]}),
               "G2": ("+", {"T2a": [(4028, 4213), (4303, 4577)]})},
        haps={"h1": [("dup", (3998, 4608)), ("del", (2599, 2611))], "h2": [("del", (2794, 2802))],
              "h3": [("del", (6133, 6134)), ("del", (4200, 4208))], "h4": [("dup", (3922, 4652)), ("del", (4575, 4583)), ("snp", 4585)],
              "h5": [("snp", 4569)]},
        reversed={"h6": "h1", "h7": "h4"}),
    # stretches duplicated in tandem inside a gene, over several nodes: one holding an exon, one starting
    # inside an exon, one ending inside an exon, one holding two exons. Each haplotype has one copy of the
    # gene, its exons whole (`intragenic`: its transcripts are the reference's).
    "dup_intragenic": dict(
        L=6000,
        genes={"GI": ("+", {"TI1": [(1000, 1200), (2000, 2200), (3000, 3200), (4000, 4200)],
                            "TI2": [(1000, 1200), (2000, 2200), (4000, 4200)]})},
        haps={"h1": [("dup", (1800, 2600))], "h2": [("dup", (2100, 2900))], "h3": [("snp", 1900), ("snp", 2400), ("snp", 2700), ("snp", 3500)],
              "h4": [("dup", (1500, 3500))], "h5": [("dup", (2900, 4100))]},
        intragenic={"h1", "h2", "h4", "h5"},
        reversed={"h6": "h2", "h7": "h5"}),
    # exon boundaries inside a substitution block: one that keeps the splice site, one that changes it
    "substitution": dict(
        L=5000,
        genes={"GS": ("+", {"TS1": [(1000, 1300), (2000, 2100), (3000, 3300)]})},
        haps={"h1": [("sub", (1296, 1306))], "h2": [("sub", (1996, 2006))], "h3": [("snp", 1300), ("snp", 2999)]},
        reversed={"h4": "h1"}),
    # substitution blocks across a gene's first and last bases, cut by other haplotypes' SNPs into several
    # nodes: the boundary keeps its offset, as inside the gene
    "substitution_ends": dict(
        L=5000,
        genes={"GS": ("+", {"TS1": [(1000, 1300), (2000, 2300)]}), "GT": ("-", {"TT1": [(3000, 3200), (3500, 3800)]})},
        haps={"h1": [("sub", (985, 1012)), ("sub", (3790, 3816))], "h2": [("snp", 990), ("snp", 995), ("snp", 3805), ("snp", 3810)]},
        reversed={"h3": "h1"}),
    # a gene whose exons lie mostly on allele nodes (placed by sequence, with a boundary inside a deletion on
    # one haplotype), and a gene most of whose exons a deletion removes (a partial copy)
    "divergent": dict(
        L=6000, sub_step=5,
        genes={"GH": ("+", {"TH1": [(1000, 1300), (2000, 2300), (3000, 3300)], "TH2": [(1000, 1300), (3000, 3300)]}),
               "GP": ("-", {"TP1": [(4500, 4700), (5000, 5300)]})},
        haps={"h1": [("sub", (1000, 1120)), ("sub", (2000, 2105))], "h2": [("del", (4950, 5350))],
              "h3": [("sub", (1000, 1120)), ("sub", (2000, 2105)), ("del", (2297, 2305))]},
        reversed={"h4": "h1"}),
    # a gene whose every exon lies on an allele node, 95% identical: no exonic node is shared, and the chain
    # through its introns and flanks points to it
    "allele_exons": dict(
        L=5000, sub_step=20,
        genes={"GQ": ("+", {"TQ1": [(1000, 1300), (2000, 2300), (3000, 3300)]})},
        haps={"h1": [("sub", (1000, 1300)), ("sub", (2000, 2300)), ("sub", (3000, 3300))]},
        reversed={"h2": "h1"}),
    # a gene placed by sequence with a deletion across the start of its first exon and one across the end of
    # its last: the ends of the gene, where an alignment has nothing past them to hold on to
    "end_indel": dict(
        L=6000, sub_step=5,
        genes={"GH": ("+", {"TH1": [(1000, 1300), (2000, 2300), (3000, 3300)]})},
        haps={"h1": [("del", (997, 1002)), ("sub", (1010, 1130)), ("sub", (2000, 2105)), ("del", (3296, 3303))]},
        reversed={"h2": "h1"}),
    # a gene whose marker loses an exon to a deletion that leaves other transcripts whole: a copy below the
    # share, with those transcripts written
    "partial": dict(
        L=8000,
        genes={"GW": ("-", {"TW1": [(1000, 1100), (3000, 3600)], "TW2": [(1000, 1100), (1400, 1500)], "TW3": [(3200, 3600)],
                            "TW4": [(1000, 1100), (1400, 1500), (4500, 4600)]})},
        haps={"h1": [("del", (2500, 4000))], "h2": [("snp", 1450)]},
        reversed={"h3": "h1"}),
    # a deleted gene whose exons a paralog elsewhere shares: on the paralog's copy, its transcripts stay
    # without a copy (write_paralog)
    "partial_paralog": dict(),
    # the same, with GX's other transcript on exons in tandem (write_tandem_paralog)
    "paralog_tandem": dict(),
    # tandem exons whole on the gene's own nodes, the marker's exon deleted: a copy below the share
    "partial_tandem": dict(),
    # a gene deleted beside a paralog that follows it on the same node (write_adjacent_paralog)
    "paralog_adjacent": dict(),
    # a short gene on the nodes of a stretch inside a longer gene, and a haplotype without the short gene
    # whose copy of the longer one carries an allele there (write_nested)
    "nested": dict(),
    # the same with the short gene's exon on the middle node alone: its chain on that haplotype holds only the
    # nodes beside the exon
    "nested_flank": dict(),
    # a segment the reference walks twice inside a gene; haplotypes holding it once, with the second copy's
    # alleles or the first's, and one lacking an exon of its second copy (write_merged)
    "merged": dict(),
    # the main case with haplotypes' own annotations (write_annotated)
    "annotated": dict(),
    # the main case with annotations the projection contradicts, fills or that lie on another gene's copy
    # (write_annotation_rules)
    "annotation_rules": dict(),
}
SHARE, IDENTITY, PARTIAL = 0.8, 0.93, 0.2   # panvar's --min-exon-share and --min-identity, and its partial copy


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
            if kind == "sub" and v not in alt_node:
                nid = str(len(nodes) + 1)
                nodes[nid] = "".join(alt(c) if (p - v[0]) % case.get("sub_step", 3) == 0 else c
                                     for p, c in enumerate(ref[v[0]:v[1]], v[0]))
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
                if kind == "sub" and v[0] <= a and b <= v[1]:   # the block's first node stands for the whole
                    take = None
                    if a == v[0]:
                        steps.append((alt_node[v], "+"))
                        for p in range(*v):
                            m[p] = pos + p - v[0]
                        pos += v[1] - v[0]
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


def project(case, ref, hseq, m, strand, exons, subs=()):
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
        if any(a < c < b for a, b in subs):   # kept at its offset inside a same-length substitution
            notes.append(f"e{e}{'.start:' if start == (strand == '+') else '.end:'}in_variant")
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
            src = case.get("reversed", {}).get(h, h)
            if src != h:   # copies are numbered along the walk as written
                copies.reverse()
            subs = [v for kind, v in case["haps"].get(src, []) if kind == "sub"]
            snps = {v for kind, v in case["haps"].get(src, []) if kind == "snp"}
            # a copy below the share still holds the transcripts with every exon on it: with one of them whole,
            # the haplotype has a copy, and its other transcripts lack an exon
            whole = len(copies) == 1 and any(project(case, ref, seqs[src], copies[0], strand, ex, subs)[0] != "missing_exon"
                                             for ex in txs.values() if any(b > 0 and a < case["L"] for a, b in ex))

            def below(t, exons, status):
                if not whole:
                    return (t, g, h, 0, status, [], "")
                st, notes, s = project(case, ref, seqs[src], copies[0], strand, exons, subs)
                return (t, g, h, 1, st, notes + ["below_share"], s)

            for t, exons in txs.items():
                if not any(b > 0 and a < case["L"] for a, b in exons):
                    continue
                if src in case.get("intragenic", ()):
                    st, notes, s = project(case, ref, ref, maps["ref"][0], strand, exons)
                    rows.append((t, g, h, 1, st, notes, s))
                    continue
                if not copies:
                    rows.append((t, g, h, 0, "no_copy", [], ""))
                by_sequence = False
                if len(copies) == 1:
                    # the marker's exonic bases the copy shares by nodes; below the share, the copy is placed by
                    # sequence when they are covered and alike enough, and is partial when covered less
                    marker = [p for a, b in next(iter(txs.values())) for p in range(max(a, 0), min(b, case["L"]))]
                    held = [p for p in marker if p in copies[0]]
                    changed = [p for p in held if p in snps or
                               any(a <= p < b and (p - a) % case.get("sub_step", 3) == 0 for a, b in subs)]
                    on_alleles = [p for p in held if p in snps or any(a <= p < b for a, b in subs)]
                    if (len(held) - len(on_alleles)) / len(marker) < SHARE:
                        alike = held and 1 - len(changed) / len(held) >= IDENTITY
                        if not alike or len(held) / len(marker) < PARTIAL:
                            rows.append(below(t, exons, "no_copy"))
                            continue
                        if len(held) / len(marker) < SHARE:
                            rows.append(below(t, exons, "partial_copy"))
                            continue
                        by_sequence = True
                for c, m in enumerate(copies, 1):
                    # a duplication's second copy walks the reference's nodes: no substitution there
                    on_subs = () if by_sequence or m is not maps[src][0] else subs
                    st, notes, s = project(case, ref, seqs[src], m, strand, exons, on_subs)
                    rows.append((t, g, h, c, st, notes + ["placed_by_sequence"] * by_sequence, s))
    return rows


def name(h, L):
    return f"ref#0#chr1:{OFFSET}-{OFFSET + L}" if h == "ref" else f"{h}#1#ctg{h[1:]}:0-{L}"


PIECES = [("a", 300), ("b", 50), ("c", 50), ("d", 50), ("e", 550)]   # the default unit, 1 kb


def write_tandem(case, d):
    """Flank A (2 kb), a 1 kb unit walked n times over the same nodes, flank B (2 kb). Gene GT has an exon in
    A, the exons of each reference unit's layout (offsets in the unit; a unit may have none), and one in B.
    A unit is a walk over pieces (`order`; a piece may come twice). A haplotype's unit j takes the layout of
    reference unit j when the counts agree. `alt` puts a changed copy of one piece in
    one unit of a haplotype; `insert` puts bases after a piece of one unit. Units are told apart by sequence:
    with identical units, only the last (which runs into flank B) differs, so a haplotype's last unit takes
    the reference's last layout and the others the layout nearest their relative position."""
    rnd = random.Random(11)
    seq = lambda n: "".join(rnd.choice("ACGT") for _ in range(n))
    A, B = list(seq(2000)), list(seq(2000))
    A[1600:1602], B[498:500] = "GT", "AG"
    A, B = "".join(A), "".join(B)
    pieces = case.get("pieces", PIECES)
    piece = {n: seq(length) for n, length in pieces}
    order = case.get("order", [n for n, _ in pieces])
    units, k = case["units"], case["units"]["ref"]
    layout = case.get("layout", [[(300, 400)]] * k)
    alt, insert = case.get("alt", {}), case.get("insert", {})
    nodes = {"A": A, "B": B, **piece}
    for h, (j, n) in alt.items():
        s = list(piece[n])
        s[len(s) // 2] = {"A": "C", "C": "G", "G": "T", "T": "A"}[s[len(s) // 2]]
        nodes[f"{n}x"] = "".join(s)
    for h, (j, i, length) in insert.items():
        nodes["I"] = seq(length) if "I" not in nodes else nodes["I"]
    ids = {n: str(k + 1) for k, n in enumerate(nodes)}

    def unit(h, j):   # (node names, sequence, insertion offset and length)
        names = list(order)
        if h in alt and alt[h][0] == j:
            names[names.index(alt[h][1])] = alt[h][1] + "x"
        ins = None
        if h in insert and insert[h][0] == j:
            at = insert[h][1] + 1
            ins = (sum(len(nodes[x]) for x in names[:at]), len(nodes["I"]))
            names.insert(at, "I")
        return names, "".join(nodes[x] for x in names), ins

    def walk(h, n):
        w = ["A"]
        for j in range(n):
            w += unit(h, j)[0]
        return w + ["B"]

    walks = {h: [(ids[x], "+") for x in walk(h, n)] for h, n in units.items()}
    for h, src in case.get("reversed", {}).items():
        walks[h] = [(ids[x], "-") for x in reversed(walk(src, units[src]))]
    edges = set()
    for steps in walks.values():
        for (x, ox), (y, oy) in zip(steps, steps[1:]):
            edges.add((x, ox, y, oy))
    L = 2000 + 1000 * k + 2000
    exons = ([(1500, 1600)] + [(2000 + 1000 * u + a, 2000 + 1000 * u + b) for u in range(k) for a, b in layout[u]] +
             [(2000 + 1000 * k + 500, 2000 + 1000 * k + 600)])
    with open(f"{d}/locus.gfa", "w") as out:
        out.write("H\tVN:Z:1.0\n")
        for n, s in nodes.items():
            out.write(f"S\t{ids[n]}\t{s}\n")
        for x, ox, y, oy in sorted(edges):
            out.write(f"L\t{x}\t{ox}\t{y}\t{oy}\t0M\n")
        for h, steps in walks.items():
            out.write(f"P\t{name(h, L)}\t{','.join(n + o for n, o in steps)}\t*\n")
    attr = 'gene_id "GT"; gene_name "GT"; gene_biotype "protein_coding"; transcript_id "TT1"; transcript_name "TT1"; tag "Ensembl_canonical";'
    with open(f"{d}/genes.gtf", "w") as out:
        out.write(f"chr1\ttest\ttranscript\t{exons[0][0] + OFFSET + 1}\t{exons[-1][1] + OFFSET}\t.\t+\t.\t{attr}\n")
        for a, b in exons:
            out.write(f"chr1\ttest\texon\t{a + OFFSET + 1}\t{b + OFFSET}\t.\t+\t.\t{attr}\n")
    ref_count = sum(len(x) for x in layout)
    with open(f"{d}/expected.tsv", "w") as out:
        out.write("transcript_id\tgene_name\thaplotype\tcopy\tstatus\tdetail\tsequence\n")
        for h in walks:
            src = case.get("reversed", {}).get(h, h)
            n = units[src]
            body, count, resized = "", 0, 0
            for j in range(n):
                _, useq, ins = unit(src, j)
                g = j if n == k else k - 1 if j == n - 1 else min(j * k // n, k - 2)
                for a, b in layout[g]:
                    if ins and ins[0] <= a:
                        a, b = a + ins[1], b + ins[1]
                    elif ins and a < ins[0] < b:
                        b += ins[1]
                        resized += 1
                    body += useq[a:b]
                    count += 1
            notes = ([f"tandem:{ref_count}>{count}"] if count != ref_count else []) + ([f"tandem_resized:{resized}"] if resized else [])
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
    rows.append(("nov1", "GA", "h1", 1, "ok", "transcript_class:putative_novel_isoform", "".join(seqs["h1"][m1[a]:m1[b - 1] + 1] for a, b in novel), "nov1"))
    # two records of one source transcript, with different exons, under one gene
    h1 += [gff("ctg1", "gene", pos(m1, 9000), pos(m1, 9149), "+", "ID=g4;gene_name=GQ")]
    for rec, end in (("cat1", 9100), ("cat2", 9150)):
        h1 += [gff("ctg1", "transcript", pos(m1, 9000), pos(m1, end - 1), "+", f"ID={rec};Parent=g4;source_transcript=ENSTQ1.1"),
               gff("ctg1", "exon", pos(m1, 9000), pos(m1, end - 1), "+", f"Parent={rec}")]
        rows.append(("ENSTQ1", "GQ", "h1", 1, "ok", "source_transcript:ENSTQ1.1", seqs["h1"][m1[9000]:m1[end - 1] + 1], rec))
    gb2 = [(6500, 6600), (6800, 6900)]       # a second GB copy, minus strand, where the projection sees none
    h1 += [gff("ctg1", "gene", pos(m1, 5000), pos(m1, 5799), "-", "ID=g2;gene_name=GB"),
           gff("ctg1", "transcript", pos(m1, 5000), pos(m1, 5799), "-", "ID=TB1;Parent=g2;transcript_id=TB1")]
    h1 += [gff("ctg1", "exon", pos(m1, a), pos(m1, b - 1), "-", "Parent=TB1") for a, b in case["genes"]["GB"][1]["TB1"]]
    h1 += [gff("ctg1", "gene", pos(m1, 6500), pos(m1, 6899), "-", "ID=g3;gene_name=GB"),
           gff("ctg1", "transcript", pos(m1, 6500), pos(m1, 6899), "-", "ID=TB1x;Parent=g3;transcript_id=TB1")]
    h1 += [gff("ctg1", "exon", pos(m1, a), pos(m1, b - 1), "-", "Parent=TB1x") for a, b in gb2]
    rows.append(("TB1", "GB", "h1", 2, "ok", ".", rc("".join(seqs["h1"][m1[a]:m1[b - 1] + 1] for a, b in gb2)), "TB1x"))
    m2 = maps["h2"][0]
    edge = [(m2[11500], m2[11700 - 1] + 1), (len(seqs["h2"]) + 100, len(seqs["h2"]) + 200)]   # path coordinates
    h2 = [gff("ctg2", "transcript", edge[0][0] + 1, edge[1][1], "+", "ID=edge1;gene_name=GZ;transcript_id=edge1")]
    h2 += [gff("ctg2", "exon", a + 1, b, "+", "Parent=edge1") for a, b in edge]
    rows.append(("edge1", "GZ", "h2", 1, "clipped_at_edge", "clipped:3':100", seqs["h2"][edge[0][0]:edge[0][1]], "edge1"))
    # a model with CDS rows only, and an exon whose transcript record is missing
    cds = [(m2[9000], m2[9099] + 1), (m2[9300], m2[9399] + 1)]
    h2 += [gff("ctg2", "mRNA", cds[0][0] + 1, cds[1][1], "-", "ID=cds1;gene_name=GW;transcript_id=cds1")]
    h2 += [gff("ctg2", "CDS", a + 1, b, "-", "Parent=cds1") for a, b in cds]
    h2 += [gff("ctg2", "exon", m2[9500] + 1, m2[9599] + 1, "+", "Parent=nowhere")]
    rows.append(("cds1", "GW", "h2", 1, "ok", ".", rc("".join(seqs["h2"][a:b] for a, b in cds)), "cds1"))
    m3 = maps["h3"][0]
    alias_ex = [(4000, 4200)]
    h3 = [gff("CTG3ALT", "transcript", pos(m3, 4000), pos(m3, 4199), "+", "ID=al1;gene_name=GY;transcript_id=al1")]
    h3 += [gff("CTG3ALT", "exon", pos(m3, a), pos(m3, b - 1), "+", "Parent=al1") for a, b in alias_ex]
    rows.append(("al1", "GY", "h3", 1, "ok", ".", seqs["h3"][m3[4000]:m3[4199] + 1], "al1"))
    # GC under another name, with its source gene: one gene for the copy count
    h3 += [gff("CTG3ALT", "gene", pos(m3, 8000), pos(m3, 8699), "+", "ID=g6;gene_name=GCALT;source_gene=GC.1"),
           gff("CTG3ALT", "transcript", pos(m3, 8000), pos(m3, 8699), "+", "ID=c6;Parent=g6;transcript_id=TC1")]
    h3 += [gff("CTG3ALT", "exon", pos(m3, a), pos(m3, b - 1), "+", "Parent=c6") for a, b in case["genes"]["GC"][1]["TC1"]]
    h5 = []
    for k, m in enumerate(maps["h5"], 1):   # GTF: both GC copies, as projected
        # the two records share their source gene under different names, the second with CAT's copy suffix:
        # one gene, two copies
        attr = f'gene_id "gc{k}"; gene_name "{"GC" if k == 1 else "GCALT"}"; transcript_id "TC1"; source_gene "GC.1{"" if k == 1 else "_2"}";'
        h5 += [f'ctg5\tCAT\texon\t{pos(m, a)}\t{pos(m, b - 1)}\t.\t+\t.\t{attr}\n' for a, b in case["genes"]["GC"][1]["TC1"]]
    # a GTF transcript whose CDS row comes before its exon row with other attributes: the exon row's count
    m5 = maps["h5"][0]
    h5 += [f'ctg5\tCAT\tCDS\t{pos(m5, 10020)}\t{pos(m5, 10079)}\t.\t+\t0\tgene_id "gv"; gene_name "WRONG"; transcript_id "TV1"; transcript_name "WRONG";\n',
           f'ctg5\tCAT\texon\t{pos(m5, 10000)}\t{pos(m5, 10099)}\t.\t+\t.\tgene_id "gv"; gene_name "GV"; transcript_id "TV1"; transcript_name "TV1";\n']
    rows.append(("TV1", "GV", "h5", 1, "ok", ".", seqs["h5"][m5[10000]:m5[10099] + 1], "TV1"))
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
        out.write(lines[0] + "\tsource\trecord\n")
        for line in lines[1:]:
            f = line.split("\t")
            both = (f[0], f[2].split("#")[0]) in (("TA1", "h1"), ("TC1", "h5"), ("TB1", "h1"), ("TC1", "h3"))
            out.write(line + ("\tprojected+annotated" if both else "\tprojected") + "\t.\n")
        for t, g, h, c, st, det, sq, rec in rows:
            out.write(f"{t}\t{g}\t{name(h, L)}\t{c}\t{st}\t{det}\t{sq}\tannotated\t{rec}\n")


def write_annotation_rules(d):
    """The main case, with annotations: h1 TA1 twice, as projected and with exon 3 ending early, and a novel
    isoform; h2 TA1 (its exon 2 deleted) with another exon 2, and without it; h4 TA1 with exon 3 starting
    later than projected, on GA's copy a gene the GTF lacks, genes the GTF has off the path (GO; one by its
    versioned id alone; one sourced by id that shares GA's name) and a novel gene (CAT's source_gene=None); h6 (no GC) GC on GB's copy, with a novel isoform, and GC on its own further
    on. expected.tsv is the default run's, expected_prefer.tsv the run with --prefer-annotation's."""
    write("main", d)
    case = CASES["main"]
    L = case["L"]
    ref, nodes, walks, seqs, maps = build(case)
    gff = lambda ctg, kind, a, b, strand, attrs: f"{ctg}\tCAT\t{kind}\t{a}\t{b}\t.\t{strand}\t.\t{attrs}\n"
    lines = {h: [] for h in ("h1", "h2", "h4", "h6")}
    rows = []   # transcript, gene, haplotype, record, exons, default (copy, status, detail), preferred (copy, status, detail)

    def gene(h, gid, attrs, lo, hi):
        m = maps[h][0]
        lines[h].append(gff(f"ctg{h[1:]}", "gene", m[lo] + 1, m[hi - 1] + 1, "+", f"ID={gid};{attrs}"))

    def tx(h, gid, rec, tid, gname, ex, default, preferred, attrs=""):
        m = maps[h][0]
        lines[h].append(gff(f"ctg{h[1:]}", "transcript", m[ex[0][0]] + 1, m[ex[-1][1] - 1] + 1, "+",
                            f"ID={rec};Parent={gid};transcript_id={tid}{attrs}"))
        lines[h] += [gff(f"ctg{h[1:]}", "exon", m[a] + 1, m[b - 1] + 1, "+", f"Parent={rec}") for a, b in ex]
        seq = "".join(seqs[h][m[a]:m[b - 1] + 1] for a, b in ex)
        rows.append((tid, gname, h, rec, seq, default, preferred))

    ta1 = case["genes"]["GA"][1]["TA1"]
    gene("h1", "a1", "gene_name=GA;source_gene=GA.1", 1000, 3300)
    tx("h1", "a1", "same1", "TA1", "GA", ta1, None, None)   # joins the projected record
    tx("h1", "a1", "short1", "TA1", "GA", ta1[:2] + [(3000, 3250)], (1, "discordant", "."), (1, "ok", "."))
    tx("h1", "a1", "nov1", "nov1", "GA", [(1000, 1300), (3000, 3150)], (1, "ok", "."), (1, "ok", "."))
    gene("h2", "a2", "gene_name=GA;source_gene=GA.1", 1000, 3300)
    tx("h2", "a2", "fill2", "TA1", "GA", [(1000, 1300), (2300, 2360), (3000, 3300)],
       (1, "ok", "annotation_fill"), (1, "ok", "annotation_fill"))
    tx("h2", "a2", "trunc2", "TA1", "GA", [(1000, 1300), (3000, 3300)], (1, "discordant", "exons:2/3"), (1, "ok", "annotation_fill"))
    gene("h4", "a4", "gene_name=GA;source_gene=GA.1", 1000, 3300)
    tx("h4", "a4", "late4", "TA1", "GA", ta1[:2] + [(3010, 3300)], (1, "discordant", "."), (1, "ok", "."))
    gene("h4", "x4", "gene_name=XL;source_gene=ENSGX9.1", 1100, 2600)
    tx("h4", "x4", "xl4", "xl4", "XL", [(1100, 1250), (2500, 2600)], (1, "ok", "."), (1, "ok", "."))
    gene("h4", "o4", "gene_name=GO", 1000, 3150)
    tx("h4", "o4", "go4", "go4", "GO", [(1000, 1200), (3050, 3150)], (0, "paralog_lift", "on:GA"), (1, "ok", "on:GA"))
    gene("h4", "e4", "gene_id=ENSG00000999901.1", 1000, 3200)
    tx("h4", "e4", "ev4", "ev4", "ENSG00000999901.1", [(1000, 1200), (3100, 3200)], (0, "paralog_lift", "on:GA"),
       (1, "ok", "on:GA"))
    gene("h4", "b4", "gene_name=GA;source_gene=ENSG00000999902.1", 1000, 3300)
    tx("h4", "b4", "gb4", "gb4", "GA", [(1000, 1300), (3050, 3300)], (0, "paralog_lift", "on:GA"), (1, "ok", "on:GA"))
    gene("h4", "n4", "gene_name=NV;source_gene=None", 1500, 1600)
    tx("h4", "n4", "nv4", "nv4", "NV", [(1500, 1600)], (1, "ok", "."), (1, "ok", "."))
    gene("h6", "c6", "gene_name=GC;source_gene=GC.1", 5000, 5700)
    tx("h6", "c6", "lift6", "TC1", "GC", [(5000, 5200), (5500, 5700)], (0, "paralog_lift", "on:GB"), (1, "ok", "on:GB"))
    tx("h6", "c6", "novl6", "novl6", "GC", [(5000, 5200)], (0, "paralog_lift", "on:GB"), (1, "ok", "on:GB"))
    gene("h6", "f6", "gene_name=GC;source_gene=GC.1_2", 10000, 10700)
    tx("h6", "f6", "fill6", "TC1", "GC", [(10000, 10200), (10500, 10700)], (1, "ok", "annotation_fill"), (2, "ok", "annotation_fill"))
    for h, ls in lines.items():
        with open(f"{d}/{h}.gff3", "w") as out:
            out.write("##gff-version 3\n")
            out.writelines(ls)
    with open(f"{d}/hap_annotations.tsv", "w") as out:
        out.write("haplotype\tpath\n")
        for h in lines:
            out.write(f"{h}#1\t{d}/{h}.gff3\n")
    with open(f"{d}/alias.tsv", "w") as out:
        out.write("# names of one contig per row\n")
    with open(f"{d}/genes.gtf", "a") as out:   # GO, past the reference path
        a = 'gene_id "GO"; gene_name "GO"; transcript_id "TO1";'
        out.write(f"chr1\ttest\texon\t{OFFSET + L + 40001}\t{OFFSET + L + 40200}\t.\t+\t.\t{a}\n")
        # off the path too: a versioned id without a name, and a gene sharing GA's name
        for k, a in ((41001, 'gene_id "ENSG00000999901.1"; transcript_id "TE1";'),
                     (42001, 'gene_id "ENSG00000999902"; gene_name "GA"; transcript_id "TB9";')):
            out.write(f"chr1\ttest\texon\t{OFFSET + L + k}\t{OFFSET + L + k + 199}\t.\t+\t.\t{a}\n")
    proj = open(f"{d}/expected.tsv").read().splitlines()
    for name_, k in (("expected.tsv", 5), ("expected_prefer.tsv", 6)):
        with open(f"{d}/{name_}", "w") as out:
            out.write(proj[0] + "\tsource\trecord\n")
            for line in proj[1:]:
                f = line.split("\t")
                h = f[2].split("#")[0]
                if f[0] == "TA1" and h == "h4" and k == 6:   # the annotation overrules the projection
                    f[4], f[5] = "discordant", "annotated:late4"
                both = (f[0], h) == ("TA1", "h1")
                out.write("\t".join(f) + ("\tprojected+annotated" if both else "\tprojected") + "\t.\n")
            for t, g, h, rec, seq, *want in rows:
                if want[k - 5] is None:
                    continue
                c, st, det = want[k - 5]
                out.write(f"{t}\t{g}\t{name(h, L)}\t{c}\t{st}\t{det}\t{seq}\tannotated\t{rec}\n")


def write_nested(d, flank=False):
    """Gene GL: exons L1 and L2 around an intron U (three nodes). Gene GS, further on: one exon, the same three
    nodes walked again, or with `flank` the middle one alone. A haplotype without GS, whose copy of GL has an
    allele of U's middle node, 96% identical: the chain of GS's nodes there lies inside GL's copy, which
    matches the stretch as well, so the haplotype has no copy of GS. The same haplotype reversed."""
    rnd = random.Random(13)
    seq = lambda n: "".join(rnd.choice("ACGT") for _ in range(n))
    size = dict(A=1000, L1=400, Ua=40, Ub=220, Uc=40, L2=400, C=1500, D=1000, T=500)
    nodes = {n: seq(k) for n, k in size.items()}
    nodes["Ubx"] = "".join(alt(c) if i % 25 == 0 else c for i, c in enumerate(nodes["Ub"]))
    ids = {n: str(k + 1) for k, n in enumerate(nodes)}
    walk = {"ref": ["A", "L1", "Ua", "Ub", "Uc", "L2", "C", "Ua", "Ub", "Uc", "D", "T"],
            "h1": ["A", "L1", "Ua", "Ubx", "Uc", "L2", "T"]}
    steps = {h: [(ids[n], "+") for n in w] for h, w in walk.items()}
    steps["h2"] = [(n, "-") for n, _ in reversed(steps["h1"])]
    L = sum(size[n] for n in walk["ref"])
    at, pos = {}, 0
    for k, n in enumerate(walk["ref"]):
        at[k] = pos
        pos += size[n]
    exons = {"GL": [(at[1], at[1] + 400), (at[5], at[5] + 400)], "GS": [(at[8], at[8] + 220)] if flank else [(at[7], at[7] + 300)]}
    with open(f"{d}/locus.gfa", "w") as out:
        out.write("H\tVN:Z:1.0\n")
        for n, s in nodes.items():
            out.write(f"S\t{ids[n]}\t{s}\n")
        for x, ox, y, oy in sorted({(x, ox, y, oy) for w in steps.values() for (x, ox), (y, oy) in zip(w, w[1:])}):
            out.write(f"L\t{x}\t{ox}\t{y}\t{oy}\t0M\n")
        for h, w in steps.items():
            out.write(f"P\t{name(h, L)}\t{','.join(n + o for n, o in w)}\t*\n")
    with open(f"{d}/genes.gtf", "w") as out:
        for g, ex in exons.items():
            attr = (f'gene_id "{g}"; gene_name "{g}"; gene_biotype "protein_coding"; transcript_id "T{g}"; '
                    f'transcript_name "T{g}"; tag "Ensembl_canonical";')
            out.write(f"chr1\ttest\ttranscript\t{ex[0][0] + OFFSET + 1}\t{ex[-1][1] + OFFSET}\t.\t+\t.\t{attr}\n")
            for a, b in ex:
                out.write(f"chr1\ttest\texon\t{a + OFFSET + 1}\t{b + OFFSET}\t.\t+\t.\t{attr}\n")
    gl, gs = nodes["L1"] + nodes["L2"], nodes["Ub"] if flank else nodes["Ua"] + nodes["Ub"] + nodes["Uc"]
    with open(f"{d}/expected.tsv", "w") as out:
        out.write("transcript_id\tgene_name\thaplotype\tcopy\tstatus\tdetail\tsequence\n")
        out.write(f"TGL\tGL\t{name('ref', L)}\t1\tok\t.\t{gl}\nTGS\tGS\t{name('ref', L)}\t1\tok\t.\t{gs}\n")
        for h in ("h1", "h2"):
            out.write(f"TGL\tGL\t{name(h, L)}\t1\tok\t.\t{gl}\nTGS\tGS\t{name(h, L)}\t0\tno_copy\t.\t.\n")


def write_merged(d):
    """Gene GM over P, a segment S, U, and S again: the reference walks S twice, over the same nodes but for one
    allele node on each side of S's middle node, the first ten bases before it (x1, y1 in the first copy; x2,
    y2 in the second). TM1 has an exon in P, one in the first S and one on the middle node of the second; TM2 one in P
    and one on the middle node of the first. h1 lacks the stretch from one copy into the other: it holds S once, with the second
    copy's alleles, so that copy is the second too and TM1's last exon is found on it (merged_copy), where
    TM2's is. h2 holds both copies and lacks the middle node in the second, and h4 holds S once with the
    first copy's alleles: the exon is missing, not taken from the first copy. h3 is h1 reversed."""
    rnd = random.Random(17)
    seq = lambda n: "".join(rnd.choice("ACGT") for _ in range(n))
    size = dict(P=1000, Sa=200, x1=1, x2=1, Sl=10, Sb=100, y1=1, y2=1, Sc=200, U=1500, Q=1000)
    nodes = {n: seq(k) for n, k in size.items()}
    nodes["x2"], nodes["y2"] = alt(nodes["x1"]), alt(nodes["y1"])
    ids = {n: str(k + 1) for k, n in enumerate(nodes)}
    S1, S2 = ["Sa", "x1", "Sl", "Sb", "y1", "Sc"], ["Sa", "x2", "Sl", "Sb", "y2", "Sc"]
    walk = {"ref": ["P", *S1, "U", *S2, "U", "Q"], "h1": ["P", *S2, "U", "Q"],
            "h2": ["P", *S1, "U", "Sa", "x2", "Sl", "y2", "Sc", "U", "Q"], "h4": ["P", *S1, "U", "Q"]}
    steps = {h: [(ids[n], "+") for n in w] for h, w in walk.items()}
    steps["h3"] = [(n, "-") for n, _ in reversed(steps["h1"])]
    L = sum(size[n] for n in walk["ref"])
    at, pos = [], 0
    for n in walk["ref"]:
        at.append(pos)
        pos += size[n]
    exons = {"TM1": [(300, 900), (at[1] + 50, at[1] + 150), (at[11], at[11] + 100)], "TM2": [(300, 900), (at[4], at[4] + 100)]}
    with open(f"{d}/locus.gfa", "w") as out:
        out.write("H\tVN:Z:1.0\n")
        for n, s in nodes.items():
            out.write(f"S\t{ids[n]}\t{s}\n")
        for x, ox, y, oy in sorted({(x, ox, y, oy) for w in steps.values() for (x, ox), (y, oy) in zip(w, w[1:])}):
            out.write(f"L\t{x}\t{ox}\t{y}\t{oy}\t0M\n")
        for h, w in steps.items():
            out.write(f"P\t{name(h, L)}\t{','.join(n + o for n, o in w)}\t*\n")
    with open(f"{d}/genes.gtf", "w") as out:
        for t, ex in exons.items():
            attr = (f'gene_id "GM"; gene_name "GM"; gene_biotype "protein_coding"; transcript_id "{t}"; '
                    f'transcript_name "{t}";' + (' tag "Ensembl_canonical";' if t == "TM1" else ""))
            out.write(f"chr1\ttest\ttranscript\t{ex[0][0] + OFFSET + 1}\t{ex[-1][1] + OFFSET}\t.\t+\t.\t{attr}\n")
            for a, b in ex:
                out.write(f"chr1\ttest\texon\t{a + OFFSET + 1}\t{b + OFFSET}\t.\t+\t.\t{attr}\n")
    tm1 = nodes["P"][300:900] + nodes["Sa"][50:150] + nodes["Sb"]
    tm2 = nodes["P"][300:900] + nodes["Sb"]
    with open(f"{d}/expected.tsv", "w") as out:
        out.write("transcript_id\tgene_name\thaplotype\tcopy\tstatus\tdetail\tsequence\n")
        for h in ("ref", "h1", "h2", "h3", "h4"):
            if h in ("h2", "h4"):
                out.write(f"TM1\tGM\t{name(h, L)}\t1\tmissing_exon\te3:missing\t.\n")
            else:
                out.write(f"TM1\tGM\t{name(h, L)}\t1\tok\t{'.' if h == 'ref' else 'e3:merged_copy'}\t{tm1}\n")
            out.write(f"TM2\tGM\t{name(h, L)}\t1\tok\t.\t{tm2}\n")


def node_locus(d, seed, size, walk, exons):
    """Writes locus.gfa and genes.gtf of a locus drawn as nodes: `size` per node, `walk` per haplotype as node
    names (an "r" + name entry is that haplotype written reversed), `exons` per gene and transcript as
    (step, start, end) on the reference walk, the first transcript the canonical one. Returns the nodes'
    sequences, the reference steps' positions and the path's length."""
    rnd = random.Random(seed)
    nodes = {n: "".join(rnd.choice("ACGT") for _ in range(k)) for n, k in size.items()}
    ids = {n: str(k + 1) for k, n in enumerate(nodes)}
    steps = {}
    for h, w in walk.items():
        steps[h] = [(ids[n], "+") for n in w] if isinstance(w, list) else [(n, "-") for n, _ in reversed(steps[w])]
    at, pos = [], 0
    for n in walk["ref"]:
        at.append(pos)
        pos += size[n]
    L = pos
    with open(f"{d}/locus.gfa", "w") as out:
        out.write("H\tVN:Z:1.0\n")
        for n, s in nodes.items():
            out.write(f"S\t{ids[n]}\t{s}\n")
        for x, ox, y, oy in sorted({(x, ox, y, oy) for w in steps.values() for (x, ox), (y, oy) in zip(w, w[1:])}):
            out.write(f"L\t{x}\t{ox}\t{y}\t{oy}\t0M\n")
        for h, w in steps.items():
            out.write(f"P\t{name(h, L)}\t{','.join(n + o for n, o in w)}\t*\n")
    with open(f"{d}/genes.gtf", "w") as out:
        for g, txs in exons.items():
            for k, (t, ex) in enumerate(txs.items()):
                ex = [(at[i] + a, at[i] + b) for i, a, b in ex]
                attr = (f'gene_id "{g}"; gene_name "{g}"; gene_biotype "protein_coding"; transcript_id "{t}"; '
                        f'transcript_name "{t}";' + (' tag "Ensembl_canonical";' if k == 0 else ""))
                out.write(f"chr1\ttest\ttranscript\t{ex[0][0] + OFFSET + 1}\t{ex[-1][1] + OFFSET}\t.\t+\t.\t{attr}\n")
                for a, b in ex:
                    out.write(f"chr1\ttest\texon\t{a + OFFSET + 1}\t{b + OFFSET}\t.\t+\t.\t{attr}\n")
    return nodes, at, L


def write_rows(d, L, rows):
    """expected.tsv from (transcript, gene, haplotype, copy, status, sequence[, detail]) rows."""
    with open(f"{d}/expected.tsv", "w") as out:
        out.write("transcript_id\tgene_name\thaplotype\tcopy\tstatus\tdetail\tsequence\n")
        for t, g, h, c, st, s, *det in rows:
            out.write(f"{t}\t{g}\t{name(h, L)}\t{c}\t{st}\t{det[0] if det else '.'}\t{s or '.'}\n")


def write_paralog(d):
    """Gene GX: exons X1, X2 and a long one in B. Gene GY, further on: X1 and X2 walked again around another
    intron, and an exon of its own. A haplotype lacks GX with its flank B. GX's nodes X1 and X2 are on it, in
    GY's copy, and TX2 has no other exon: but the reference walks them outside GX too, so they are not GX's
    there, and GX has no copy. The same haplotype reversed."""
    size = dict(A=1000, X1=300, I1=500, X2=300, B=3000, S=200, I2=500, C=500, Y3=300, T=500)
    walk = {"ref": ["A", "X1", "I1", "X2", "B", "S", "X1", "I2", "X2", "C", "Y3", "T"],
            "h1": ["A", "X1", "I2", "X2", "C", "Y3", "T"], "h2": "h1"}
    exons = {"GX": {"TX1": [(1, 0, 300), (3, 0, 300), (4, 200, 2800)], "TX2": [(1, 0, 300), (3, 0, 300)]},
             "GY": {"TY1": [(6, 0, 300), (8, 0, 300), (10, 0, 300)]}}
    nodes, _, L = node_locus(d, 19, size, walk, exons)
    tx2 = nodes["X1"] + nodes["X2"]
    spliced = {"TX1": tx2 + nodes["B"][200:2800], "TX2": tx2, "TY1": tx2 + nodes["Y3"]}
    write_rows(d, L, [(t, g, h, *((0, "no_copy", "") if g == "GX" and h != "ref" else (1, "ok", spliced[t])))
                      for h in ("ref", "h1", "h2") for g, txs in exons.items() for t in txs])


def write_tandem_paralog(d, paralog):
    """Gene GX: TX1 on its own node M, TX2 on X twice around I, between B and C. With `paralog`, gene GY,
    further on, walks B, X, I, X and C again between exons of its own, and a haplotype lacks GX: GX's best
    chain, through B, I and C, lies on GY's copy, where TX2's exons would be placed unit by unit, but the
    reference walks those nodes outside GX too, so GX has no copy. Without, a haplotype lacks M alone: TX2 is
    whole on GX's own nodes, a copy below the share. The same haplotype reversed."""
    size = dict(A=1000, M=300, B=400, X=200, I=300, C=400, S=1500, Y1=300, Y2=300, T=500)
    gx = ["B", "X", "I", "X", "C"]
    walk = ({"ref": ["A", "M", *gx, "S", "Y1", *gx, "Y2", "T"], "h1": ["A", "S", "Y1", *gx, "Y2", "T"]} if paralog
            else {"ref": ["A", "M", *gx, "T"], "h1": ["A", *gx, "T"]})
    walk["h2"] = "h1"
    exons = {"GX": {"TX1": [(1, 0, 300)], "TX2": [(3, 0, 150), (5, 0, 150)]}}
    if paralog:
        exons["GY"] = {"TY1": [(8, 0, 300), (14, 0, 300)]}
    nodes, _, L = node_locus(d, 23, size, walk, exons)
    spliced = {"TX1": nodes["M"], "TX2": nodes["X"][:150] * 2, "TY1": nodes["Y1"] + nodes["Y2"]}

    def row(t, g, h):
        if g == "GY" or h == "ref":
            return (t, g, h, 1, "ok", spliced[t])
        if paralog:
            return (t, g, h, 0, "no_copy", "")
        return (t, g, h, 1, "ok", spliced[t], "below_share") if t == "TX2" else (t, g, h, 1, "missing_exon", "", "e1:missing;below_share")
    write_rows(d, L, [row(t, g, h) for h in ("ref", "h1", "h2") for g, txs in exons.items() for t in txs])


def write_adjacent_paralog(d):
    """Gene GX: TX1 on M and X, TX2 on X alone. Gene GY, right after it: X walked again, and Y. A haplotype
    lacks M and the first X; the X left is too little of TX1 for a partial copy. GX's stretch of the reference reaches one step past its last exon, onto GY's X,
    but X is not GX's own: GX has no copy, and GY keeps its X. The same haplotype reversed."""
    size = dict(A=1000, M=1200, X=300, Y=300, Q=1000)
    walk = {"ref": ["A", "M", "X", "X", "Y", "Q"], "h1": ["A", "X", "Y", "Q"], "h2": "h1"}
    exons = {"GX": {"TX1": [(1, 0, 1000), (2, 50, 250)], "TX2": [(2, 50, 250)]},
             "GY": {"TY1": [(3, 50, 250), (4, 50, 250)]}}
    nodes, _, L = node_locus(d, 29, size, walk, exons)
    x = nodes["X"][50:250]
    spliced = {"TX1": nodes["M"][:1000] + x, "TX2": x, "TY1": x + nodes["Y"][50:250]}
    write_rows(d, L, [(t, g, h, *((0, "no_copy", "") if g == "GX" and h != "ref" else (1, "ok", spliced[t])))
                      for h in ("ref", "h1", "h2") for g, txs in exons.items() for t in txs])


def write(case_name, d):
    case = CASES[case_name]
    if case_name.startswith("tandem"):
        return write_tandem(case, d)
    if case_name.startswith("nested"):
        return write_nested(d, case_name == "nested_flank")
    if case_name == "annotated":
        return write_annotated(d)
    if case_name == "annotation_rules":
        return write_annotation_rules(d)
    if case_name == "merged":
        return write_merged(d)
    if case_name == "partial_paralog":
        return write_paralog(d)
    if case_name in ("paralog_tandem", "partial_tandem"):
        return write_tandem_paralog(d, case_name == "paralog_tandem")
    if case_name == "paralog_adjacent":
        return write_adjacent_paralog(d)
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


def check(d, prefix, expected_name="expected.tsv"):
    exp = {}
    with open(f"{d}/{expected_name}") as fh:
        fh.readline()
        for line in fh:
            t, g, h, c, st, det, s, *src = line.rstrip("\n").split("\t")
            source = src[0] if src else "projected"
            exp[(t, h, c, source, src[1] if source == "annotated" else ".")] = (st, det, s, g)
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
            k = (r["transcript_id"], r["haplotype"], r["copy"], r["source"],
                 r.get("record") if r["source"] == "annotated" else ".")
            if k in got:
                sys.exit(f"duplicate row: {k}")
            key = f"{r['transcript_id']}|{r['gene_name']}|{r['haplotype']}|{r['copy']}"
            key += f"|annotated:{r.get('record')}" if r["source"] == "annotated" else ""
            got[k] = (r, fa.get(key) or ".")
    bad = 0
    for k in sorted(set(exp) | set(got)):
        if k not in got or k not in exp:
            print(f"  row {'missing' if k not in got else 'unexpected'}: {k}")
            bad += 1
            continue
        st, det, s, g = exp[k]
        r, seq = got[k]
        problems = []
        if r["gene_name"] != g:
            problems.append(f"gene {r['gene_name']} != {g}")
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
        sys.exit(check(*sys.argv[2:5]))
