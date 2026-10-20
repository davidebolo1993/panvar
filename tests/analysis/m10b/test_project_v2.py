#!/usr/bin/env python3
"""M10-B.1 fixtures for project_v2.py. Toy classes and toy class-defining sites are FIXTURE ASSUMPTIONS. The adversarial
'A' fixtures EXPOSE limitations of the 50-bp convention: each records whether the intended distinction is preserved, and
a limitation is expected and reported, never tuned away."""
import hashlib
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import project  # noqa: E402  (version 1, only to document its missing-walk defect)
import project_v2 as p2  # noqa: E402

fails, limits = 0, []


def check(cond, what):
    global fails
    print(("  PASS " if cond else "  FAIL ") + what)
    fails += 0 if cond else 1


def expose(name, preserved, observed):
    limits.append((name, preserved, observed))
    print(f"  {'PRESERVED ' if preserved else 'LIMITATION'} {name}: {observed}")


def rnd(tag, n):
    out, i = "", 0
    while len(out) < n:
        out += "".join("ACGT"[b % 4] for b in hashlib.sha256(f"{tag}:{i}".encode()).digest())
        i += 1
    return out[:n]


def write_model(d, name, nodes, blocks, kinds=None, trims=None, sources=None, seq_override=None, walks_recorded=None):
    """blocks: {b: [walk or None, ...]}; sequences follow the implementation's boundary transformation."""
    kinds = kinds or {b: ("SV" if b % 2 else "context") for b in blocks}
    trims = trims or {}
    sp = lambda w: "".join(nodes[s[:-1]] if s[-1] == "+" else p2.rc(nodes[s[:-1]]) for s in w)  # noqa: E731
    seqs = {}
    for b, ws in blocks.items():
        for a, w in enumerate(ws):
            if seq_override and (b, a) in seq_override:
                seqs[(b, a)] = seq_override[(b, a)]
                continue
            s = sp(w or [])
            if kinds[b] == "SV":
                P, S = trims.get(b, (0, 0))
                s = s[P:len(s) - S]
            else:
                if b - 1 in blocks and kinds[b - 1] == "SV" and trims.get(b - 1, (0, 0))[1]:
                    s = sp(blocks[b - 1][0])[-trims[b - 1][1]:] + s
                if b + 1 in blocks and kinds[b + 1] == "SV" and trims.get(b + 1, (0, 0))[0]:
                    s = s + sp(blocks[b + 1][0])[:trims[b + 1][0]]
            seqs[(b, a)] = s
    p = os.path.join(d, name)
    with open(p + ".block_alleles.tsv", "w") as f:
        f.write("block_index\tblock_id\tkind\tallele_index\tlength\tsequence\n")
        for (b, a), s in sorted(seqs.items()):
            f.write(f"{b}\t{('sv:' if kinds[b] == 'SV' else 'ctx:')}{b}\t{kinds[b]}\t{a}\t{len(s)}\t{s}\n")
    with open(p + ".walks.tsv", "w") as f:
        f.write("block_index\tallele_index\twalk_index\twalk\n")
        for b, ws in blocks.items():
            for a, w in enumerate(ws):
                rec = (walks_recorded or {}).get((b, a), [w] if w is not None else [])
                for i, x in enumerate(rec):
                    f.write(f"{b}\t{a}\t{i}\t{','.join(x)}\n")
    with open(p + ".sources.tsv", "w") as f:
        f.write("block_index\tallele_index\twalk_sources\tempty_sources\tbypass_sources\tother_sources\n")
        for b, ws in blocks.items():
            for a, w in enumerate(ws):
                src = (sources or {}).get((b, a), (1 if w is not None else 0, 0, 0, 0))
                f.write(f"{b}\t{a}\t" + "\t".join(map(str, src)) + "\n")
    with open(p + ".blocks.tsv", "w") as f:
        f.write("block_index\tmodel_block_id\treport_block_id\tkind\ttrimmed_prefix_bp\ttrimmed_suffix_bp\trepresentation_dependent\t"
                "bypass_templates\tsource_node\tsink_node\n")
        for b in blocks:
            P, S = trims.get(b, (0, 0))
            f.write(f"{b}\tb{b}\tb{b}\t{kinds[b]}\t{P}\t{S}\t0\t0\t-\t-\n")
    with open(p + ".nodes.tsv", "w") as f:
        f.write("node\tlength\tsequence\n")
        for k, v in nodes.items():
            f.write(f"{k}\t{len(v)}\t{v}\n")
    open(p + ".panel.tsv", "w").write("template\talleles\nT1\t" + ",".join("0" for _ in blocks) + "\n")
    open(p + ".summary.tsv", "w").write("key\tvalue\nmodel_fingerprint\ttoy\n")
    return p


N = {k: rnd(k, 60) for k in ("X", "Y", "Z", "U", "V", "J1", "J2", "B", "L")}
N.update({"s0": "A", "s1": "C", "c0": "GAT", "c1": "GTT", "Ydup": N["Y"], "g1": "G"})
S1, S2 = rnd("S1", 60), rnd("S2", 60)
N.update({"S1": S1, "S2": S2, "M1": S1[:-5], "M2": S1[-5:] + S2})
N.update({f"u{i}": N["U"][6 * i:6 * i + 6] for i in range(10)})          # U as ten 6-bp nodes
N.update({f"i{i}": rnd("INS", 55)[5 * i:5 * i + 5] for i in range(11)})   # a 55-bp insertion as eleven 5-bp nodes
N["Usnp"] = N["U"][:30] + ("A" if N["U"][30] != "A" else "C") + N["U"][31:]   # U with one base changed, a long node
N["B1"], N["B2"] = rnd("BB", 30) + rnd("B1", 30), rnd("BB", 30) + rnd("B2", 30)   # share 30 bytes

BLOCKS = {0: [["X+"]],
          1: [["U+"], ["U+", "U+"], ["U+", "U+", "U+"], ["U+", "s0+", "V+"], ["U+", "s1+", "V+"], ["V+", "s0+", "U+"],
              ["U-", "s0+", "V+"], ["U+", "J1+", "V+"], ["U+", "J2+", "V+"], ["Y+"], [], ["U+", "c0+", "V+"], ["U+", "c1+", "V+"],
              ["S1+"], ["M1+"],
              ["L+"] + [f"u{i}+" for i in range(10)] + ["s0+"],   # a15 U spelled by short nodes, plus a neutral SNP
              ["L+", "U+"],                                        # a16 U as one long node
              ["U+", "V+"],                                        # a17 U V
              ["U+"] + [f"i{i}+" for i in range(11)] + ["V+"],     # a18 U, 55-bp insertion of short nodes, V
              ["Usnp+", "s0+", "V+"],                              # a19 the SNP inside a long node
              None,                                                # a20 walk metadata absent (non-empty sequence)
              None,                                                # a21 another, unrelated, absent walk
              ["J1+", "s1+"],                                      # a22 filler (a10 is the valid empty record)
              ["U+", "V+", "s1+"]],                                # a23 partially missing: 2 walk sources, 1 recorded
          2: [["Z+"], ["Ydup+", "Z+"], ["S2+"], ["M2+"], ["S2+", "s1+"]]}
SEQ_OVERRIDE = {(1, 20): rnd("missing-a", 70), (1, 21): rnd("missing-b", 70)}
# a10: a valid empty record (an 'empty' source, no walk); the model keeps one allele per sequence, so only one empty allele
SOURCES = {(1, 10): (0, 1, 0, 0), (1, 20): (0, 0, 0, 0), (1, 21): (0, 0, 0, 0), (1, 23): (2, 0, 0, 0)}
RECORDED = {(1, 10): []}

with tempfile.TemporaryDirectory() as d:
    prefix = write_model(d, "toy", N, BLOCKS, seq_override=SEQ_OVERRIDE, sources=SOURCES, walks_recorded=RECORDED)
    m = p2.Model(prefix)
    pj = lambda t, mm=m: mm.project(mm.spell(0, t), 0, 2)  # noqa: E731
    dip = lambda a, b, mm=m: p2.compare([pj(a[0], mm), pj(a[1], mm)], [pj(b[0], mm), pj(b[1], mm)])[0]  # noqa: E731
    R = (0, 0, 0)
    # --- ported version-1 fixtures (verdict names of version 2)
    check(pj((0, 9, 0))[1] == pj((0, 10, 1))[1] and pj((0, 9, 0))[0] == "RESOLVED", "F1 synonymous split + duplicate node: identical")
    check(dip(((0, 3, 0), R), ((0, 4, 0), R)) == "UNIQUE_AGREEMENT", "F2 neutral short SNP: agreement (declared convention)")
    check(dip(((0, 1, 0), R), ((0, 2, 0), R)) == "DISAGREEMENT", "F4 copy number 2 vs 3: disagreement")
    v = p2.compare([pj((0, 2, 0)), pj(R)], [pj((0, 1, 0)), pj((0, 1, 0))])
    check(v == ("DISAGREEMENT", "DOSAGE_EQUAL_ALLOCATION_DIFFERS"), f"F5 3+1 vs 2+2 at equal dosage: {v}")
    check(dip(((0, 2, 0), R), (R, (0, 2, 0))) == "UNIQUE_AGREEMENT", "F5b homologue exchange: agreement")
    check(dip(((0, 3, 0), R), ((0, 5, 0), R)) == "DISAGREEMENT" and dip(((0, 3, 0), R), ((0, 6, 0), R)) == "DISAGREEMENT",
          "F6 order and orientation: disagreement")
    check(dip(((0, 7, 0), R), ((0, 8, 0), R)) == "DISAGREEMENT", "F7 junction J1 vs J2: disagreement")
    check(pj((0, 13, 2))[0] == "AMBIGUOUS" and dip(((0, 13, 2), R), ((0, 14, 3), R)) == "SET_IDENTICAL", "F9 one sequence, two parses: SET_IDENTICAL")
    check(dip(((0, 13, 2), R), ((0, 13, 4), R)) == "COMPATIBLE", "F9c intersecting sets: COMPATIBLE")
    off = m.spell(0, (0, 3, 0)).replace(N["V"][:10], "N" * 10)
    po = m.project(off, 0, 2)
    check(po[0] == "UNKNOWN" and p2.compare([po, pj(R)], [po, pj(R)])[0] == "UNSUPPORTED", "F10 off-model: UNKNOWN; UNKNOWN vs UNKNOWN unsupported")
    # --- missing metadata (item 2)
    check(pj((0, 20, 0))[0] == "INCOMPLETE", "M1 absent walk, non-empty allele: INCOMPLETE")
    check(dip(((0, 20, 0), R), ((0, 21, 0), R)) == "UNSUPPORTED", "M1b two unrelated absent-walk sequences: UNSUPPORTED, never agreement")
    check(pj((0, 10, 0))[0] == "RESOLVED" and "valid empty" in pj((0, 10, 0))[3], "M2 valid empty record: RESOLVED (empty contribution)")
    check(pj((0, 23, 0))[0] == "INCOMPLETE", "M3 partially missing alternative walks (2 sources, 1 recorded): INCOMPLETE")
    v1 = project.Model(prefix)
    a1, a2 = v1.project(v1.spell(0, (0, 20, 0)), 0, 2), v1.project(v1.spell(0, (0, 21, 0)), 0, 2)
    b0 = v1.project(v1.spell(0, R), 0, 2)
    check(a1[0] == "RESOLVED" and project.compare([a1, b0], [a2, b0])[0] == "EQUAL",
          "M4 version-1 defect reproduced: absent walks RESOLVED and two unrelated sequences EQUAL (version 1 kept unchanged)")
    # --- adversarial: the long-node convention (item 3)
    expose("A1 same unit as one long node vs ten short nodes (+ neutral SNP)", dip(((0, 15, 0), R), ((0, 16, 0), R)) == "UNIQUE_AGREEMENT",
           dip(((0, 15, 0), R), ((0, 16, 0), R)))
    expose("A2 a 55-bp insertion made only of short nodes", dip(((0, 18, 0), R), ((0, 17, 0), R)) == "DISAGREEMENT",
           dip(((0, 18, 0), R), ((0, 17, 0), R)))
    expose("A3 neutral SNP inside a long node vs as a separate short node", dip(((0, 19, 0), R), ((0, 3, 0), R)) == "UNIQUE_AGREEMENT",
           dip(((0, 19, 0), R), ((0, 3, 0), R)))
    expose("A4 shifted block boundaries spelling the same compound sequence", pj((0, 9, 0))[1] == pj((0, 10, 1))[1], "identical projections")
    c0, c1 = N["c0"], N["c1"]
    md = p2.Model(prefix, class_defining_sequences=[c0, c1])
    undecl = dip(((0, 11, 0), R), ((0, 12, 0), R))
    decl = p2.compare([pj((0, 11, 0), md), pj(R, md)], [pj((0, 12, 0), md), pj(R, md)])[0]
    neutral_decl = p2.compare([pj((0, 3, 0), md), pj(R, md)], [pj((0, 4, 0), md), pj(R, md)])[0]
    expose("A5 class-defining short change vs neutral short change", undecl == "UNIQUE_AGREEMENT" and decl == "DISAGREEMENT" and
           neutral_decl == "UNIQUE_AGREEMENT", f"undeclared {undecl}; declared class site {decl}; neutral under declaration {neutral_decl}")
    check(undecl == "UNIQUE_AGREEMENT" and decl == "DISAGREEMENT" and neutral_decl == "UNIQUE_AGREEMENT",
          "A5 is distinguished only through the explicit declaration")
    mg = p2.Model(prefix, class_defining_sequences=[N["g1"]])   # declare the 1-bp site 'G' by sequence class
    captured = p2.compare([pj((0, 3, 0), mg), pj(R, mg)], [pj((0, 4, 0), mg), pj(R, mg)])[0]
    expose("A5b a 1-bp class site declared by sequence class (G) leaves the neutral C change ignored",
           captured == "UNIQUE_AGREEMENT", f"neutral A-vs-C change under the G declaration: {captured} "
           "(G and C share one canonical class, so a sequence-class declaration cannot name a short site)")
    # --- byte-verified class identity: forced collision
    mc = p2.Model(prefix, hasher=lambda b: "0000000000")
    same = [dip(x, y, mc) == dip(x, y) for x, y in [(((0, 1, 0), R), ((0, 2, 0), R)), (((0, 3, 0), R), ((0, 5, 0), R)),
                                                    (((0, 7, 0), R), ((0, 8, 0), R)), (((0, 3, 0), R), ((0, 4, 0), R))]]
    check(all(same) and len(set(mc.registry.values())) == len(mc.registry) and any("~" in v for v in mc.registry.values()),
          "C1 forced hash collision: classes stay distinct by bytes (ids disambiguated), every verdict unchanged")
    # --- scope and boundary transformation
    tb = {0: [["X+"]], 1: [["B+", "U+"], ["B+", "V+"]], 2: [["Z+"]]}
    pt = p2.Model(write_model(d, "trim", N, tb, trims={1: (60, 0)}))
    au = pt.audit()
    s_sv = pt.project(pt.seq[(1, 0)], 1, 1)
    check(all(r["verified"] for r in au.values()) and s_sv[0] == "RESOLVED" and len(next(iter(s_sv[1]))) == 1,
          "S1 trimmed variant block: relation verified; clip on a node boundary drops the moved boundary node")
    check(pt.project(pt.seq[(0, 0)], 0, 0)[0] == "UNSUPPORTED", "S3 context scope receiving bytes from a trimmed neighbour: UNSUPPORTED")
    check(pt.project(pt.spell(0, (0, 0, 0)), 0, 2)[0] == "RESOLVED", "S0 whole haplotype over trimmed blocks: RESOLVED")
    tc = {0: [["X+"]], 1: [["B1+", "U+"], ["B2+", "V+"]], 2: [["Z+"]]}
    pc = p2.Model(write_model(d, "cut", N, tc, trims={1: (30, 0)}))
    check(pc.project(pc.seq[(1, 0)], 1, 1)[0] == "UNSUPPORTED", "S2 the clip cuts a long node: UNSUPPORTED")
    pw = p2.Model(write_model(d, "wrong", N, tb, trims={1: (60, 0)}))
    pw.trim[1] = (59, 0)   # the report disagrees with the bytes
    check(not pw.audit()[1]["verified"] and pw.project(pw.seq[(1, 0)], 1, 1)[0] == "UNSUPPORTED", "S4 trim not reproducible: audit fails, UNSUPPORTED")

print("limitations exposed: " + "; ".join(f"{n} -> {'preserved' if ok else 'NOT preserved'}" for n, ok, _ in limits))
print("m10b.1 fixtures: all passed" if not fails else "m10b.1 fixtures: FAILURES")
sys.exit(1 if fails else 0)
