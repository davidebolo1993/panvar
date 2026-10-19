#!/usr/bin/env python3
"""M6-E analysis (panel only; frozen with DESIGN.md).

Inputs: the prepare unit's tables, the alignment events, and M1-D's panel component FASTA (allowlisted; nothing else).

Correspondence (per component type; reference coordinates of the frozen reference component, W = 10 bp):
  SIMPLE       every component aligns a base here and no component has an indel (or a large-deletion edge) within W
  CONDITIONAL  absent (inside a deletion >= 50 bp) in some components, simple in every component that has it
  UNCERTAIN    inside a small deletion or within W of any indel / large-deletion edge in some component
Sites: reference positions with >= 1 SNV in some component; alleles counted by physical copy, haplotype and donor.
Co-occurrence:
  (a) within-copy linkage among informative SIMPLE sites (minor allele >= 3 physical copies; non-reference vs reference):
      r^2 across physical copies; null = each site permuted independently across copies (20 permutations, seed 20261019)
  (b) within-haplotype copy similarity on sites with simple correspondence in BOTH period and tail coordinates (mapped
      through the single tail-reference vs period-reference alignment, excluding its indel neighbourhoods): Hamming
      distance of copy pairs in the same haplotype against all copy pairs from different donors with the same type
      composition (median percentile; 0.5 = no excess similarity); length-class concordance likewise
  (c) unique-sequence structure: pairwise Hamming between unique components over SIMPLE SNV sites; single-linkage
      cluster counts at 0, 5, 20 and 50 differences

  analyse.py <prepare dir> <events tsv> <out prefix>
"""
import builtins
import csv
import os
import sys
from collections import Counter, defaultdict

import numpy as np

prep, events_path, prefix = sys.argv[1:4]
M1D_FA = "/Users/davide.bolognini/panvar/tests/results/m1d_unit_reachability_20261019/runs/inventory.a1/inventory.fa"
ALLOW = {os.path.abspath(p) for p in [M1D_FA, events_path, f"{prep}/components.tsv", f"{prep}/copies.tsv", f"{prep}/prepare.tsv", f"{prep}/pairs.tsv"]}
_open = builtins.open


def guarded(path, mode="r", *a, **k):
    p = os.path.abspath(path)
    if "w" in mode or "a" in mode:
        assert p.startswith(os.path.dirname(os.path.abspath(prefix))), f"write outside the unit directory: {p}"
    else:
        assert p in ALLOW, f"not an allowlisted panel-only input: {p}"
    return _open(path, mode, *a, **k)


builtins.open = guarded
W, LARGE, MINOR, PERMS = 10, 50, 3, 20
T = lambda p: list(csv.DictReader(open(p), delimiter="\t"))


def fasta(path):
    recs, name = {}, None
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            name = line[1:].split()[0]
            recs[name] = []
        elif name:
            recs[name].append(line)
    return {k: "".join(v) for k, v in recs.items()}


comp = fasta(M1D_FA)
comps = {r["unique_record"]: r for r in T(f"{prep}/components.tsv")}
copies = T(f"{prep}/copies.tsv")
refs = {l.split("\t")[1]: l.split("\t")[2].split()[0] for l in open(f"{prep}/prepare.tsv") if l.startswith("reference\t")}
ev = defaultdict(list)
summ = {}
for r in T(events_path):
    if r["event"] == "SUMMARY":
        summ[r["label"]] = r
    else:
        ev[r["label"]].append(r)
S = []   # summary rows


def say(k, v):
    S.append((k, v))


class Aligned:
    """One component aligned to a reference: base at each reference position (None where absent)."""

    def __init__(self, label, L):
        self.snv = {}
        self.absent = np.zeros(L, dtype=bool)
        self.large = np.zeros(L, dtype=bool)
        self.unc = np.zeros(L, dtype=bool)
        self.end_gap = False
        for e in ev[label]:
            p, n = int(e["ref_pos"]), int(e["ref_len"])
            if e["event"] == "SNV":
                self.snv[p] = e["query_bases"]
            elif e["event"] == "DEL":
                self.absent[p:p + n] = True
                if n >= LARGE:
                    self.large[p:p + n] = True
                    self.unc[max(0, p - W):min(L, p + W)] = True
                    self.unc[max(0, p + n - W):min(L, p + n + W)] = True
                else:
                    self.unc[max(0, p - W):min(L, p + n + W)] = True
                self.end_gap |= p == 0 or p + n == L
            else:   # INS before p
                self.unc[max(0, p - W):min(L, p + W)] = True
                self.end_gap |= p == 0 or p == L
        self.unc |= self.absent & ~self.large


def analyse_type(kind, label_prefix, ref_name, members):
    R = comp[ref_name]
    L = len(R)
    al = {k: Aligned(f"{label_prefix}:{k}", L) for k in members}
    absent_large = np.sum([a.large for a in al.values()], axis=0) if al else np.zeros(L)
    unc = np.any([a.unc for a in al.values()], axis=0)
    cls = np.where(unc, 2, np.where(absent_large > 0, 1, 0))   # 0 SIMPLE, 1 CONDITIONAL, 2 UNCERTAIN
    out = dict(R=R, L=L, al=al, cls=cls)
    say(f"{label_prefix}.reference", f"{ref_name} ({L} bp)")
    for c, nm in enumerate(["SIMPLE", "CONDITIONAL", "UNCERTAIN"]):
        say(f"{label_prefix}.positions_{nm}", f"{int(np.sum(cls == c))} ({np.mean(cls == c):.4f})")
    say(f"{label_prefix}.components_with_end_gaps", sum(a.end_gap for a in al.values()))
    say(f"{label_prefix}.distance_to_reference", "min %s median %s max %s" % tuple(
        int(x) for x in np.percentile([int(summ[f'{label_prefix}:{k}']['distance']) for k in members], [0, 50, 100], method="nearest")))
    # merged non-simple intervals
    iv = []
    for c in (1, 2):
        pos = np.flatnonzero(cls == c)
        if len(pos):
            br = np.flatnonzero(np.diff(pos) > 1)
            for s_, e_ in zip(np.r_[pos[0], pos[br + 1]], np.r_[pos[br], pos[-1]]):
                n_comp = sum(1 for a in al.values() if (a.unc[s_:e_ + 1].any() if c == 2 else a.large[s_:e_ + 1].any()))
                iv.append((["SIMPLE", "CONDITIONAL", "UNCERTAIN"][c], int(s_), int(e_) + 1, int(e_ - s_ + 1), n_comp))
    out["intervals"] = iv
    say(f"{label_prefix}.uncertain_intervals", f"{sum(1 for i in iv if i[0] == 'UNCERTAIN')} intervals, {sum(i[3] for i in iv if i[0] == 'UNCERTAIN')} bp; largest "
        + ",".join(f"{i[1]}-{i[2]}({i[3]} bp,{i[4]} comps)" for i in sorted((i for i in iv if i[0] == 'UNCERTAIN'), key=lambda i: -i[3])[:3]))
    say(f"{label_prefix}.conditional_intervals", ",".join(f"{i[1]}-{i[2]}({i[3]} bp, absent in {i[4]} comps)" for i in iv if i[0] == "CONDITIONAL") or "none")
    # clustered SNVs (another SNV within W in the same component)
    clustered = set()
    for a in al.values():
        ps = sorted(a.snv)
        for x, y in zip(ps, ps[1:]):
            if y - x <= W:
                clustered.update((x, y))
    out["clustered"] = clustered
    return out


def copy_alleles(kind, A, positions):
    """Matrix sites x physical copies of the given type: allele characters ('' where absent)."""
    cps = [c for c in copies if c["type"] == kind]
    M = []
    for x in positions:
        row = []
        for c in cps:
            a = A["al"][c["unique_record"]]
            row.append("" if a.absent[x] else a.snv.get(x, A["R"][x]))
        M.append(row)
    return cps, M


def site_table(kind, A, label_prefix):
    cps = [c for c in copies if c["type"] == kind]
    positions = sorted({p for a in A["al"].values() for p in a.snv})
    _, M = copy_alleles(kind, A, positions)
    rows = []
    for x, row in zip(positions, M):
        cnt = Counter(v for v in row if v)
        hap, don = defaultdict(set), defaultdict(set)
        for v, c in zip(row, cps):
            if v:
                hap[v].add(c["template"])
                don[v].add(c["donor"])
        alle = [v for v, _ in cnt.most_common()]
        minor = cnt[alle[1]] if len(alle) > 1 else 0
        rows.append(dict(type=label_prefix, ref_pos=x, cls=["SIMPLE", "CONDITIONAL", "UNCERTAIN"][A["cls"][x]], ref_base=A["R"][x], alleles=len(alle),
                         clustered=int(x in A["clustered"]), present_copies=sum(cnt.values()), absent_copies=sum(1 for v in row if not v),
                         copies=",".join(f"{v}:{cnt[v]}" for v in alle), haplotypes=",".join(f"{v}:{len(hap[v])}" for v in alle),
                         donors=",".join(f"{v}:{len(don[v])}" for v in alle), minor_copies=minor,
                         minor_haplotypes=len(hap[alle[1]]) if len(alle) > 1 else 0, minor_donors=len(don[alle[1]]) if len(alle) > 1 else 0))
    for c in ("SIMPLE", "CONDITIONAL", "UNCERTAIN"):
        rs = [r for r in rows if r["cls"] == c]
        say(f"{label_prefix}.snv_sites_{c}", f"{len(rs)} (multiallelic {sum(1 for r in rs if r['alleles'] > 2)}, clustered {sum(r['clustered'] for r in rs)})")
    for unit in ("copies", "haplotypes", "donors"):
        b = Counter()
        for r in rows:
            if r["cls"] != "UNCERTAIN" and r[f"minor_{unit}"] > 0:
                m = r[f"minor_{unit}"]
                b["1" if m == 1 else "2" if m == 2 else "3-5" if m <= 5 else "6-20" if m <= 20 else ">20"] += 1
        say(f"{label_prefix}.minor_allele_spectrum_by_{unit}", ", ".join(f"{k}:{b[k]}" for k in ["1", "2", "3-5", "6-20", ">20"]))
    return rows


def ld(kind, A, rows, label_prefix, rng):
    inf = [r for r in rows if r["cls"] == "SIMPLE" and r["minor_copies"] >= MINOR]
    say(f"{label_prefix}.informative_sites", len(inf))
    if len(inf) < 2:
        return
    pos = [r["ref_pos"] for r in inf]
    cps, M = copy_alleles(kind, A, pos)
    X = np.array([[1.0 if v != A["R"][x] else 0.0 for v in row] for x, row in zip(pos, M)])
    iu = np.triu_indices(len(pos), 1)
    r2 = np.corrcoef(X)[iu] ** 2
    d = np.abs(np.subtract.outer(pos, pos))[iu]
    say(f"{label_prefix}.ld_pairs", len(r2))
    say(f"{label_prefix}.ld_r2", f"median {np.median(r2):.3f}; >=0.8 {np.mean(r2 >= 0.8):.4f}; =1 {np.mean(r2 >= 1 - 1e-9):.4f}")
    for lo, hi in [(0, 1000), (1000, 5000), (5000, 10**9)]:
        m = (d >= lo) & (d < hi)
        if m.any():
            say(f"{label_prefix}.ld_r2_distance_{lo}-{hi if hi < 10**9 else 'inf'}", f"pairs {int(m.sum())}; median {np.median(r2[m]):.3f}; >=0.8 {np.mean(r2[m] >= 0.8):.4f}")
    nul8, nul1 = [], []
    for _ in range(PERMS):
        Y = np.array([rng.permutation(x) for x in X])
        q = np.corrcoef(Y)[iu] ** 2
        nul8.append(np.mean(q >= 0.8))
        nul1.append(np.mean(q >= 1 - 1e-9))
    say(f"{label_prefix}.ld_null_independent_sites", f">=0.8 mean {np.mean(nul8):.4f} (max {np.max(nul8):.4f}); =1 mean {np.mean(nul1):.4f}")
    # perfectly linked groups and distinct unit patterns
    parent = list(range(len(pos)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    for k in np.flatnonzero(r2 >= 1 - 1e-9):
        a, b = find(iu[0][k]), find(iu[1][k])
        parent[a] = b
    groups = Counter(find(i) for i in range(len(pos)))
    say(f"{label_prefix}.perfectly_linked_groups", f"{len(groups)} groups for {len(pos)} sites; largest {max(groups.values())}")
    say(f"{label_prefix}.distinct_copy_patterns", f"{len({tuple(col) for col in X.T})} patterns over {X.shape[1]} physical copies")


def unique_structure(kind, A, rows, label_prefix):
    pos = [r["ref_pos"] for r in rows if r["cls"] == "SIMPLE"]
    names = sorted(A["al"])
    X = np.array([[A["al"][k].snv.get(x, A["R"][x]) for x in pos] for k in names]) if pos else np.zeros((len(names), 0))
    D = np.array([[int(np.sum(X[i] != X[j])) for j in range(len(names))] for i in range(len(names))])
    iu = np.triu_indices(len(names), 1)
    say(f"{label_prefix}.unique_pairwise_hamming_simple_sites", "min %d median %d max %d" % (D[iu].min(), np.median(D[iu]), D[iu].max()))
    for th in (0, 5, 20, 50):
        parent = list(range(len(names)))

        def find(i):
            while parent[i] != i:
                parent[i] = parent[parent[i]]
                i = parent[i]
            return i
        for i, j in zip(*iu):
            if D[i, j] <= th:
                parent[find(i)] = find(j)
        say(f"{label_prefix}.single_linkage_clusters_le_{th}", len({find(i) for i in range(len(names))}))


rng = np.random.default_rng(20261019)
members = lambda kind: sorted((k for k in comps if comps[k]["type"] == kind), key=lambda k: int(k.split("_u")[1]))
A = {}
site_rows = []
for kind, lab, ref in [("prefix", "prefix", refs["prefix"]), ("period", "period", refs["period"]), ("tail", "tail", refs["tail"])]:
    A[lab] = analyse_type(kind, lab, ref, members(kind))
    rows = site_table(kind, A[lab], lab)
    site_rows += rows
    if kind != "prefix":
        ld(kind, A[lab], rows, lab, rng)
        unique_structure(kind, A[lab], rows, lab)
# reference sensitivity (periods against the most frequent short period): correspondence classes only
As = analyse_type("period", "period_short", refs["period_short"], members("period"))
site_table("period", As, "period_short")

# cross mapping: tail reference (query) aligned to the period reference
P, Tl = A["period"], A["tail"]
cross = Aligned("cross:tailref_vs_periodref", P["L"])
p2t = np.full(P["L"], -1)
qi, ri = 0, 0
evs = sorted(ev["cross:tailref_vs_periodref"], key=lambda e: (int(e["ref_pos"]), 0 if e["event"] == "INS" else 1))
k = 0
while ri < P["L"]:
    while k < len(evs) and int(evs[k]["ref_pos"]) == ri and evs[k]["event"] == "INS":
        qi += int(evs[k]["query_len"])
        k += 1
    if k < len(evs) and int(evs[k]["ref_pos"]) == ri and evs[k]["event"] == "DEL":
        ri += int(evs[k]["ref_len"])
        k += 1
        continue
    if k < len(evs) and int(evs[k]["ref_pos"]) == ri and evs[k]["event"] == "SNV":
        k += 1
    p2t[ri] = qi
    qi += 1
    ri += 1
cs = summ["cross:tailref_vs_periodref"]
say("cross.tailref_vs_periodref", f"distance {cs['distance']}; SNV {cs['ref_len']}; INS {cs['query_pos']} ({cs['query_len']} bp); DEL {cs['ref_bases']} ({cs['query_bases']} bp)")
mapped = (p2t >= 0) & ~cross.unc
say("cross.period_positions_mapped_simple", f"{int(mapped.sum())} of {P['L']} ({mapped.mean():.4f})")

# within-haplotype copy similarity on jointly simple, mapped sites
joint = [s for s in range(P["L"]) if P["cls"][s] == 0 and mapped[s] and Tl["cls"][p2t[s]] == 0]
cps = [c for c in copies if c["type"] in ("period", "tail")]


def allele(c, s):
    if c["type"] == "period":
        a = P["al"][c["unique_record"]]
        return a.snv.get(s, P["R"][s])
    t = p2t[s]
    a = Tl["al"][c["unique_record"]]
    return a.snv.get(t, Tl["R"][t])


vary = []
for s in joint:
    col = [allele(c, s) for c in cps]
    cnt = Counter(col).most_common()
    if len(cnt) > 1 and cnt[1][1] >= MINOR:
        vary.append(col)
say("joint.sites_simple_in_both_and_mapped", len(joint))
say("joint.informative_sites", len(vary))
X = np.array(vary).T if vary else np.zeros((len(cps), 0))
lclass = [comps[c["unique_record"]]["length_class"] for c in cps]
within = defaultdict(list)
between = defaultdict(list)
cls_w, cls_b = defaultdict(list), defaultdict(list)
for i in range(len(cps)):
    for j in range(i + 1, len(cps)):
        key = "-".join(sorted((cps[i]["type"], cps[j]["type"])))
        dist = int(np.sum(X[i] != X[j]))
        same = lclass[i] == lclass[j]
        if cps[i]["template"] == cps[j]["template"]:
            within[key].append(dist)
            cls_w[key].append(same)
        elif cps[i]["donor"] != cps[j]["donor"]:
            between[key].append(dist)
            cls_b[key].append(same)
for key in sorted(within):
    w, b = np.array(within[key]), np.sort(np.array(between[key]))
    pct = [(np.searchsorted(b, x, "left") + np.searchsorted(b, x, "right")) / 2 / len(b) for x in w]
    say(f"joint.within_haplotype_{key}", f"pairs {len(w)}; mean Hamming {w.mean():.2f} vs between-donor {b.mean():.2f} ({len(b)} pairs); "
        f"median percentile {np.median(pct):.3f}; identical {np.mean(w == 0):.3f} vs {np.mean(b == 0):.3f}")
    say(f"joint.length_class_concordance_{key}", f"same class within haplotype {np.mean(cls_w[key]):.3f} vs between donors {np.mean(cls_b[key]):.3f}")

# outputs
with open(f"{prefix}.sites.tsv", "w") as f:
    w = csv.DictWriter(f, fieldnames=list(site_rows[0]), delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(site_rows)
with open(f"{prefix}.intervals.tsv", "w") as f:
    f.write("type\tclass\tbegin\tend\tbp\tcomponents\n")
    for lab in ("prefix", "period", "tail", "period_short"):
        src = As if lab == "period_short" else A[lab]
        for i in src["intervals"]:
            f.write(f"{lab}\t{i[0]}\t{i[1]}\t{i[2]}\t{i[3]}\t{i[4]}\n")
with open(f"{prefix}.summary.tsv", "w") as f:
    for k_, v in S:
        f.write(f"{k_}\t{v}\n")
print(open(f"{prefix}.summary.tsv").read())
