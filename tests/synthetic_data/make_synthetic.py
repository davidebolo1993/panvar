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
             s_peakdup1/2 G x3 (+1), s_null G x1 (-1); with --cn all three are carriers
             of one module CNV, the loss included.
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

print("wrote tests/synthetic_data/syn.gfa")
print(f"nodes={len(order)} edges={len(edges)} paths={len(paths)}")


# ---------------------------------------------------------------------------
# Second graph: W-lines (instead of P-lines) + a multi-copy locus whose copies
# appear in MIXED orientation (forward/reverse). Exercises W-line parsing and
# orientation handling end to end.
# ---------------------------------------------------------------------------
wnames = {}; wseqs = {}; worder = []
def wmk(name, length):
    if name in wnames: return wnames[name]
    nid = str(len(worder) + 1)
    wnames[name] = nid; wseqs[nid] = "".join(random.choice("ACGT") for _ in range(length))
    worder.append(nid); return nid

W = [wmk(f"W{i}", 40) for i in range(3)]
IV = wmk("IV", 80)   # inversion locus
U  = wmk("U", 60)    # repeat unit, folded via U-Cc-U (no self-loop)
Cc = wmk("Cc", 15)
def wstep(name, rev=False): return (wnames[name], rev)

def wvariant(**kw):
    p = []
    p += [wstep("W0")] + kw.get("inv_block", [wstep("IV")]) + [wstep("W1")]
    p += kw.get("dup_block", [wstep("U"), wstep("Cc"), wstep("U")]) + [wstep("W2")]
    return p

# walks: reference U x2 (forward); s_winv inverts IV; s_wmixdup has U x3 as +,-,+
wpaths = {}
wpaths[("synref", 0)]   = wvariant()
wpaths[("s_winv", 0)]   = wvariant(inv_block=[wstep("IV", True)])
wpaths[("s_wmixdup", 0)]= wvariant(dup_block=[wstep("U"), wstep("Cc"), wstep("U", True), wstep("Cc"), wstep("U")])  # +,-,+
wpaths[("s_wref2", 0)]  = wvariant()

wedges = set()
for p in wpaths.values():
    for (na, ra), (nb, rb) in zip(p, p[1:]):
        wedges.add((na, "-" if ra else "+", nb, "-" if rb else "+"))

with open("tests/synthetic_data/syn_w.gfa", "w") as f:
    f.write("H\tVN:Z:1.1\n")
    for nid in worder:
        f.write(f"S\t{nid}\t{wseqs[nid]}\n")
    for a, ao, b, bo in sorted(wedges):
        f.write(f"L\t{a}\t{ao}\t{b}\t{bo}\t0M\n")
    for (sample, hap), p in wpaths.items():
        span = sum(len(wseqs[nid]) for nid, _ in p)
        walk = "".join(("<" if rev else ">") + nid for nid, rev in p)
        # W  sample  hap  seqid  start  end  walk
        f.write(f"W\t{sample}\t{hap}\tchrW\t0\t{span}\t{walk}\n")

print("wrote tests/synthetic_data/syn_w.gfa")
print(f"w-nodes={len(worder)} w-edges={len(wedges)} w-walks={len(wpaths)}")
