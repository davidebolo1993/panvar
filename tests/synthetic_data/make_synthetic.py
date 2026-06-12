#!/usr/bin/env python3
"""Generate a tiny synthetic pangenome exercising every call event type + confounders.

One reference path + sample haplotypes, each carrying ONE isolated event so the
expected VCF is unambiguous. Loci (between single-copy anchors A0..A6); each indel
bubble keeps a "stay" node so neither allele has an empty interior (which is how real
bubbles look, and which the caller's interval finder requires):
  DEL      : A0-S1-DEL1-A1; s_del drops DEL1 (keeps S1).
  INS      : A1-S2-A2;      s_ins inserts INS1 after S2 (novel sequence).
  INV      : reference has INV1 forward; s_inv traverses it reverse.
  SELFDUP  : R has a self-loop; reference R x1, s_selfdup R x3 (self-loop DUP).
  PEAKDUP  : folded G via G-C-G cycle (no self-loop). reference G x2 (baseline),
             s_peakdup1/2 G x3 (peak-multiplicity DUP, +1), s_null G x1 (a deletion,
             must NOT fire as a DUP).
  CONFOUND : B folded x2 in ALL paths (invariant cluster background) + node D whose
             deletion (s_confdel) makes the bubble exist. B must NOT yield a DUP.
"""
import random
random.seed(7)

names = {}            # name -> node id (str)
seqs = {}             # id -> sequence
order = []            # node ids in creation order
def mk(name, length):
    if name in names:
        return names[name]
    nid = str(len(order) + 1)
    names[name] = nid
    seqs[nid] = "".join(random.choice("ACGT") for _ in range(length))
    order.append(nid)
    return nid

# nodes
A = [mk(f"A{i}", 40) for i in range(7)]
S1 = mk("S1", 50)   # "stay" node so the DEL bubble interior is never empty
S2 = mk("S2", 50)   # "stay" node so the INS bubble interior is never empty
DEL1 = mk("DEL1", 60)
INS1 = mk("INS1", 70)
INV1 = mk("INV1", 80)
R    = mk("R", 60)
C    = mk("C", 15)
G    = mk("G", 60)
B    = mk("B", 60)
Cb   = mk("Cb", 15)
D    = mk("D", 55)

def step(name, rev=False):
    return (names[name], rev)

def variant(**kw):
    """Build a path = reference with one locus replaced."""
    p = []
    # DEL locus: A0 - S1 - [DEL1] - A1
    p += [step("A0"), step("S1")] + kw.get("del_block", [step("DEL1")]) + [step("A1")]
    # INS locus: S2 - [INS1?] - A2
    p += [step("S2")] + kw.get("ins_block", []) + [step("A2")]
    # INV locus: [INV1] - A3
    p += kw.get("inv_block", [step("INV1")]) + [step("A3")]
    # SELFDUP locus: [R ...] - A4
    p += kw.get("self_block", [step("R")]) + [step("A4")]
    # PEAKDUP locus: [G C G ...] - A5
    p += kw.get("peak_block", [step("G"), step("C"), step("G")]) + [step("A5")]
    # CONFOUND locus: B Cb B - [D] - A6
    p += [step("B"), step("Cb"), step("B")] + kw.get("conf_block", [step("D")]) + [step("A6")]
    return p

paths = {}
paths["synref#0#chrS:1-2000"]     = variant()                                      # reference
paths["s_del#0#chrS:1-2000"]      = variant(del_block=[])                          # drop DEL1
paths["s_ins#0#chrS:1-2000"]      = variant(ins_block=[step("INS1")])              # insert INS1
paths["s_inv#0#chrS:1-2000"]      = variant(inv_block=[step("INV1", True)])        # INV1 reverse
paths["s_selfdup#0#chrS:1-2000"]  = variant(self_block=[step("R"), step("R"), step("R")])  # R x3
paths["s_peakdup1#0#chrS:1-2000"] = variant(peak_block=[step("G"), step("C"), step("G"), step("C"), step("G")])  # G x3
paths["s_peakdup2#0#chrS:1-2000"] = variant(peak_block=[step("G"), step("C"), step("G"), step("C"), step("G")])  # G x3
paths["s_null#0#chrS:1-2000"]     = variant(peak_block=[step("G")])                # G x1 (deletion, not a DUP)
paths["s_confdel#0#chrS:1-2000"]  = variant(conf_block=[])                         # drop D
paths["s_ref2#0#chrS:1-2000"]     = variant()                                      # identical to reference

# edges: all oriented adjacencies across paths (dedup), orientation-aware
edges = set()
for p in paths.values():
    for (na, ra), (nb, rb) in zip(p, p[1:]):
        edges.add((na, "-" if ra else "+", nb, "-" if rb else "+"))

with open("tests/synthetic_data/syn.gfa", "w") as f:
    f.write("H\tVN:Z:1.0\n")
    for nid in order:
        f.write(f"S\t{nid}\t{seqs[nid]}\n")
    for a, ao, b, bo in sorted(edges):
        f.write(f"L\t{a}\t{ao}\t{b}\t{bo}\t0M\n")
    for name, p in paths.items():
        walk = ",".join(f"{nid}{'-' if rev else '+'}" for nid, rev in p)
        f.write(f"P\t{name}\t{walk}\t*\n")

# expected-truth file for the smoke assertions
with open("tests/synthetic_data/expected.tsv", "w") as f:
    f.write("svtype\tcarriers\tnote\n")
    f.write("DEL\ts_del\tdeletion of DEL1 (60bp)\n")
    f.write("INS\ts_ins\tnovel insertion INS1 (70bp)\n")
    f.write("INV\ts_inv\tinversion of INV1 (80bp)\n")
    f.write("DUP\ts_selfdup\tR x3 vs ref x1 (self-loop), REF_CN=1\n")
    f.write("DUP\ts_peakdup1,s_peakdup2\tG x3 vs ref x2 (peak multiplicity, +1), REF_CN=2\n")
    f.write("DEL\ts_null\tG x1 vs ref x2 = deletion, must NOT be a DUP\n")
    f.write("DEL\ts_confdel\tdeletion of D; B folded x2 invariant must NOT yield a DUP\n")
print("wrote tests/synthetic_data/syn.gfa")
print(f"nodes={len(order)} edges={len(edges)} paths={len(paths)}")
