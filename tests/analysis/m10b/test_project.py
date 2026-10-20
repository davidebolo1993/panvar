#!/usr/bin/env python3
"""M10-B frozen fixtures for project.py (toy model in the m10b_walk_dump format). Toy classes and the toy class-defining
site are FIXTURE ASSUMPTIONS, not biological rules."""
import hashlib
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import project  # noqa: E402

fails = 0


def check(cond, what):
    global fails
    print(("  PASS " if cond else "  FAIL ") + what)
    fails += 0 if cond else 1


def rnd(tag, n):   # deterministic pseudo-random DNA
    out, i = "", 0
    while len(out) < n:
        out += "".join("ACGT"[b % 4] for b in hashlib.sha256(f"{tag}:{i}".encode()).digest())
        i += 1
    return out[:n]


# nodes: long (60 bp) units and junctions, a duplicate-sequence node, short SNP nodes, split nodes for the ambiguity case
N = {"X": rnd("X", 60), "Y": rnd("Y", 60), "Z": rnd("Z", 60), "U": rnd("U", 60), "V": rnd("V", 60), "J1": rnd("J1", 60),
     "J2": rnd("J2", 60), "s0": "A", "s1": "C", "c0": "G", "c1": "T"}
N["Ydup"] = N["Y"]                                  # a second node with Y's sequence (graph redundancy)
S1, S2 = rnd("S1", 60), rnd("S2", 60)
N["S1"], N["S2"], N["M1"], N["M2"] = S1, S2, S1[:-5], S1[-5:] + S2   # one sequence, two long-node segmentations


def sp(walk):
    return "".join(N[s[:-1]] if s[-1] == "+" else project.rc(N[s[:-1]]) for s in walk)


# blocks: 0 ctx (flank), 1 sv, 2 ctx; alleles listed as walks (sequence = spelled walk). Two context-2 alleles give the
# synonymous split: sv allele [X,Y] + ctx [Z]  ==  sv allele [X] + ctx [Ydup,Z].
W = {0: [["X+"]],
     1: [["U+"],                                  # a0  one unit
         ["U+", "U+"],                            # a1  two units (CN change)
         ["U+", "U+", "U+"],                      # a2  three units
         ["U+", "s0+", "V+"],                     # a3  U V with SNP s0
         ["U+", "s1+", "V+"],                     # a4  U V with SNP s1 (neutral small difference)
         ["V+", "s0+", "U+"],                     # a5  order change
         ["U-", "s0+", "V+"],                     # a6  orientation change
         ["U+", "J1+", "V+"],                     # a7  junction J1
         ["U+", "J2+", "V+"],                     # a8  junction J2
         ["Y+"],                                  # a9  for the synonymous split
         [],                                      # a10 empty allele
         ["U+", "c0+", "V+"],                     # a11 toy class site c0
         ["U+", "c1+", "V+"],                     # a12 toy class site c1
         ["S1+"],                                 # a13 ambiguity parse A (sv part)
         ["M1+"]],                                # a14 ambiguity parse B (sv part)
     2: [["Z+"], ["Ydup+", "Z+"], ["S2+"], ["M2+"], ["S2+", "s1+"]]}


def write_model(d, class_defining_site=False):
    p = os.path.join(d, "toy")
    with open(p + ".block_alleles.tsv", "w") as f:
        f.write("block_index\tblock_id\tkind\tallele_index\tlength\tsequence\n")
        for b, alleles in W.items():
            for a, w in enumerate(alleles):
                s = sp(w)
                f.write(f"{b}\t{['ctx:0', 'sv:1', 'ctx:2'][b]}\t{'SV' if b == 1 else 'context'}\t{a}\t{len(s)}\t{s}\n")
    with open(p + ".walks.tsv", "w") as f:
        f.write("block_index\tallele_index\twalk_index\twalk\n")
        for b, alleles in W.items():
            for a, w in enumerate(alleles):
                f.write(f"{b}\t{a}\t0\t{','.join(w)}\n")
    with open(p + ".nodes.tsv", "w") as f:
        f.write("node\tlength\tsequence\n")
        for k, v in N.items():
            f.write(f"{k}\t{len(v)}\t{v}\n")
    with open(p + ".panel.tsv", "w") as f:   # intact templates: the mosaic (0,2,0)+... below is deliberately absent
        f.write("template\talleles\nT1\t0,0,0\nT2\t0,3,0\nT3\t0,1,0\n")
    with open(p + ".summary.tsv", "w") as f:
        f.write("key\tvalue\nmodel_fingerprint\ttoy\n")
    return p


with tempfile.TemporaryDirectory() as d:
    m = project.Model(write_model(d))
    pj = lambda t: m.project(m.spell(0, t), 0, 2)  # noqa: E731
    dip = lambda a, b: project.compare([pj(a[0]), pj(a[1])], [pj(b[0]), pj(b[1])])  # noqa: E731
    H = lambda *t: tuple(t)  # noqa: E731

    # F1 synonymous block split (and a duplicate-sequence node): same sequence, different labelling -> identical projection
    t_a, t_b = H(0, 9, 0), H(0, 10, 1)
    check(m.spell(0, t_a) == m.spell(0, t_b) and pj(t_a)[1] == pj(t_b)[1] and pj(t_a)[0] == "RESOLVED",
          "F1 synonymous split + duplicate node: identical, RESOLVED projection")
    # F2 neutral small difference (SNP s0 vs s1, below the convention, not declared): EQUAL for a declared reason
    check(dip((H(0, 3, 0), H(0, 0, 0)), (H(0, 4, 0), H(0, 0, 0)))[0] == "EQUAL", "F2 neutral SNP ignored (below ANCHOR_BP, undeclared)")
    # F3 declared class-defining short difference: c0 vs c1 differ only when declared
    check(dip((H(0, 11, 0), H(0, 0, 0)), (H(0, 12, 0), H(0, 0, 0)))[0] == "EQUAL", "F3a toy class site undeclared: ignored")
    cdef = {"U" + hashlib.sha256(min(N[k], project.rc(N[k])).encode()).hexdigest()[:10] for k in ("c0", "c1")}
    md = project.Model(m.prefix, class_defining=cdef)
    pjd = lambda t: md.project(md.spell(0, t), 0, 2)  # noqa: E731
    check(project.compare([pjd(H(0, 11, 0)), pjd(H(0, 0, 0))], [pjd(H(0, 12, 0)), pjd(H(0, 0, 0))])[0] == "DISJOINT",
          "F3b toy class site declared (explicit rule): DISJOINT")
    # F4 copy-number change
    check(dip((H(0, 1, 0), H(0, 0, 0)), (H(0, 2, 0), H(0, 0, 0)))[0] == "DISJOINT", "F4 copy number 2 vs 3: DISJOINT")
    # F5 3+1 vs 2+2 at equal total dosage
    v = dip((H(0, 2, 0), H(0, 0, 0)), (H(0, 1, 0), H(0, 1, 0)))
    check(v == ("DISJOINT", "DOSAGE_EQUAL_ALLOCATION_DIFFERS"), f"F5 3+1 vs 2+2: DISJOINT, dosage equal flagged ({v})")
    check(dip((H(0, 2, 0), H(0, 0, 0)), (H(0, 0, 0), H(0, 2, 0)))[0] == "EQUAL", "F5b homologue exchange: EQUAL")
    # F6 order and orientation
    check(dip((H(0, 3, 0), H(0, 0, 0)), (H(0, 5, 0), H(0, 0, 0)))[0] == "DISJOINT", "F6a order U,V vs V,U: DISJOINT")
    check(dip((H(0, 3, 0), H(0, 0, 0)), (H(0, 6, 0), H(0, 0, 0)))[0] == "DISJOINT", "F6b orientation U+ vs U-: DISJOINT")
    # F7 junction change (long junction nodes)
    check(dip((H(0, 7, 0), H(0, 0, 0)), (H(0, 8, 0), H(0, 0, 0)))[0] == "DISJOINT", "F7 junction J1 vs J2: DISJOINT")
    # F8 allowed mosaic absent from the intact templates: projected like any spelled sequence
    mosaic = H(0, 2, 0)
    check(list(mosaic) not in m.panel.values() and pj(mosaic)[0] == "RESOLVED" and
          dip((mosaic, H(0, 0, 0)), (mosaic, H(0, 0, 0)))[0] == "EQUAL", "F8 allowed mosaic not in the panel: RESOLVED, EQUAL to itself")
    # F9 ambiguous parsing: S1+S2 spelled with long nodes S1|S2 or M1|M2 -> set-valued; same sequence -> EQUAL_SET_VALUED
    amb_a, amb_b = H(0, 13, 2), H(0, 14, 3)
    pa = pj(amb_a)
    check(m.spell(0, amb_a) == m.spell(0, amb_b) and pa[0] == "AMBIGUOUS" and len(pa[1]) == 2 and pa[1] == pj(amb_b)[1],
          "F9a same sequence, two structural parses: AMBIGUOUS (2 skeletons), labelling-invariant")
    check(dip((amb_a, H(0, 0, 0)), (amb_b, H(0, 0, 0)))[0] == "EQUAL_SET_VALUED", "F9b identical set-valued projections: EQUAL_SET_VALUED")
    near = H(0, 13, 4)   # S1 + S2 + SNP: only the S1|S2 parse exists
    check(pj(near)[0] == "RESOLVED" and dip((amb_a, H(0, 0, 0)), (near, H(0, 0, 0)))[0] == "COMPATIBLE",
          "F9c set intersecting a resolved projection: COMPATIBLE (not EQUAL)")
    # F10 missing structural information: an off-model sequence -> UNKNOWN; two UNKNOWNs are never equated
    off = m.spell(0, H(0, 3, 0)).replace(N["V"][:10], "N" * 10)
    po = m.project(off, 0, 2)
    check(po[0] == "UNKNOWN", "F10a off-model sequence: UNKNOWN")
    check(project.compare([po, pj(H(0, 0, 0))], [po, pj(H(0, 0, 0))])[0] == "UNKNOWN", "F10b UNKNOWN vs the same UNKNOWN: UNKNOWN, never EQUAL")
    # F11 matching dosage alone is not a structural match (order differs at equal class counts)
    check(dip((H(0, 3, 0), H(0, 0, 0)), (H(0, 5, 0), H(0, 0, 0)))[0] != "EQUAL", "F11 equal class counts, different order: not EQUAL")
print("m10b fixtures: all passed" if not fails else "m10b fixtures: FAILURES")
sys.exit(1 if fails else 0)
