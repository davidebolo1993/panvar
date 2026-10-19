#!/usr/bin/env python3
"""M3-H read-level decomposition of the matched pair (diagnostic; truth-guided only through the frozen pair; reads the
score-set --decompose outputs; changes nothing).

A = R_t0 (truth ctx:0 pair (14,56)), B = R_c0 (called (4,52)); delta = contribution_A - contribution_B per read pair (each
pair once, the scorer's own masses). Reported, at 2x and at 5x:
  * counts favouring A, favouring B, exactly neutral (|delta| <= 1e-9) and effectively neutral (|delta| < 0.01 nats);
  * concentration: signed and absolute delta carried by the top 1/5/10/20/50 pairs; delta by placement class;
  * nesting: 5x pairs matched to 2x by name (verified identical bytes upstream); the 5x delta split into shared vs added;
    shared pairs' 2x vs 5x deltas compared (the libraries differ in fragment_start_intensity only);
  * discriminating content at ctx:0: 31-mers of one candidate's block allele absent from the other's, per homologue, and
    their clustered positions (block coordinates);
  * strongest contributors (top 10 by |delta| at 2x; top 10 added at 5x): per mate, gapless placements (25-mer seeds,
    both strands) on the four haplotypes -- best mismatch count, offsets at the best count (placement ambiguity), overlap
    with ctx:0 and with discriminating k-mers -- the FR insert of the best same-haplotype mate placement, base qualities,
    and the scorer's origin counts and class. Gapless placement cannot see indels: a read with no placement within 8
    mismatches is reported unresolved, not inferred.
Flag (frozen): a contributor with |delta| >= 0.5 whose best proper-geometry placement (insert within mean +- 4 sd)
has strictly fewer total mismatches on the candidate it scores AGAINST is a CANDIDATE_INCONSISTENCY for targeted
inspection -- reported, not interpreted as a bug.

  read_detail.py <score 2x dir> <score 5x dir> <reads 2x dir> <reads 5x dir> <spell prefix> <out prefix>
"""
import csv
import sys
from collections import defaultdict

s2, s5, r2, r5, spell, prefix = sys.argv[1:7]
K, SEED, MAXMM = 31, 25, 8
MEAN, SD = 350.0, 40.0
COMP = str.maketrans("ACGTN", "TGCAN")
rc = lambda s: s.translate(COMP)[::-1]


def table(p):
    return list(csv.DictReader(open(p), delimiter="\t"))


def fq(p):
    L = open(p).read().split("\n")
    return [(L[i][1:].split()[0].rsplit("/", 1)[0], L[i + 1], L[i + 3]) for i in range(0, len(L) - 3, 4)]


def kv(p):
    return {l.split("\t")[0]: l.rstrip("\n").split("\t")[1] for l in open(p) if "\t" in l}


fa, name = {}, None
for line in open(spell + ".fa"):
    line = line.strip()
    if line.startswith(">"):
        name = line[1:]
        fa[name] = []
    else:
        fa[name].append(line)
fa = {k: "".join(v) for k, v in fa.items()}
coords = {r["record"]: (int(r["block_begin"]), int(r["block_end"])) for r in table(spell + ".spelled.tsv")}
HAPS = {"A1": "R_t0.1", "A2": "R_t0.2", "B1": "R_c0.1", "B2": "R_c0.2"}
out = []
# identity: truth bytes
if "truth.1" in fa:
    tr = sorted([fa["truth.1"], fa["truth.2"]])
    out.append(f"sequence identity: R_t0 spells the truth pair {sorted([fa['R_t0.1'], fa['R_t0.2']]) == tr}; R_c0 differs from it "
               f"{sorted([fa['R_c0.1'], fa['R_c0.2']]) != tr}")
for k in (1, 2):
    a, b = fa[f"R_t0.{k}"], fa[f"R_c0.{k}"]
    ca, cb = coords[f"R_t0.{k}"], coords[f"R_c0.{k}"]
    out.append(f"homologue {k}: outside ctx:0 identical {a[:ca[0]] == b[:cb[0]] and a[ca[1]:] == b[cb[1]:]}; ctx:0 bp A {ca[1] - ca[0]} B {cb[1] - cb[0]}")
# discriminating k-mers per homologue (block coordinates)
disc = {}
for k in (1, 2):
    for x, y in (("A", "B"), ("B", "A")):
        hx, hy = HAPS[f"{x}{k}"], HAPS[f"{y}{k}"]
        bx, by = fa[hx][coords[hx][0]:coords[hx][1]], fa[hy][coords[hy][0]:coords[hy][1]]
        ky = {by[i:i + K] for i in range(len(by) - K + 1)}
        pos = [i for i in range(len(bx) - K + 1) if bx[i:i + K] not in ky]
        clusters = []
        for p in pos:
            if clusters and p - clusters[-1][1] <= K:
                clusters[-1][1] = p
            else:
                clusters.append([p, p])
        disc[f"{x}{k}"] = set(pos)
        out.append(f"discriminating {K}-mers in {x}{k} ({hx}) absent from {y}{k}: {len(pos)} in {len(clusters)} cluster(s) at block offsets "
                   + ", ".join(f"{c[0]}-{c[1] + K}" for c in clusters[:12]) + (" ..." if len(clusters) > 12 else ""))


def frag_table(d):
    stem = [p for p in __import__("glob").glob(d + "/score.decomp.*.fragments.tsv")][0]
    return {int(r["fragment"]): (float(r["delta"]), r["class"]) for r in table(stem)}, stem


def summary(d):
    return kv([p for p in __import__("glob").glob(d + "/score.decomp.*.summary.tsv")][0])


reads2, reads5 = fq(r2 + "/reads_1.fq"), fq(r5 + "/reads_1.fq")
mates2, mates5 = fq(r2 + "/reads_2.fq"), fq(r5 + "/reads_2.fq")
d2, stem2 = frag_table(s2)
d5, stem5 = frag_table(s5)
idx2 = {reads2[i][0]: i for i in range(len(reads2))}
for label, d, sm in (("2x", d2, summary(s2)), ("5x", d5, summary(s5))):
    vals = [v[0] for v in d.values()]
    order = sorted(vals, key=lambda x: -abs(x))
    tot, absum = sum(vals), sum(abs(x) for x in vals)
    out.append(f"[{label}] pairs {len(vals)}; delta(A-B) {tot:.6f} (scorer {float(sm['fragment_term_A_scorer']) - float(sm['fragment_term_B_scorer']):.6f}; "
               f"canonical-minus-scorer A {sm['fragment_term_A_from_canonical_minus_scorer']} B {sm['fragment_term_B_from_canonical_minus_scorer']}; "
               f"placement_evidence_interpretable {sm['placement_evidence_interpretable']}; exposure A-B "
               f"{float(sm['exposure_term_A']) - float(sm['exposure_term_B']):.3e})")
    out.append(f"[{label}] favour A (truth) {sum(x > 1e-9 for x in vals)}, favour B (called) {sum(x < -1e-9 for x in vals)}, exactly neutral "
               f"{sum(abs(x) <= 1e-9 for x in vals)}; effectively neutral (|d|<0.01) {sum(abs(x) < 0.01 for x in vals)}; non-neutral A "
               f"{sum(x >= 0.01 for x in vals)} B {sum(x <= -0.01 for x in vals)}; sum|d| {absum:.4f}")
    out.append(f"[{label}] concentration (signed / share of sum|d|): " + "; ".join(
        f"top{t} {sum(order[:t]):.4f} / {sum(abs(x) for x in order[:t]) / absum:.3f}" for t in (1, 5, 10, 20, 50)))
    by = defaultdict(lambda: [0, 0.0])
    for v, c in d.values():
        by[c][0] += 1
        by[c][1] += v
    out.append(f"[{label}] by class: " + "; ".join(f"{c} {n} pairs {s:.4f}" for c, (n, s) in sorted(by.items())))
# nesting
shared, added = 0.0, 0.0
diffs = []
names5 = [r[0] for r in reads5]
for i, nm in enumerate(names5):
    if nm in idx2:
        shared += d5[i][0]
        diffs.append((abs(d5[i][0] - d2[idx2[nm]][0]), nm, d2[idx2[nm]][0], d5[i][0]))
    else:
        added += d5[i][0]
diffs.sort(reverse=True)
sign_flips = sum(1 for _, _, a, b in diffs if abs(a) >= 0.01 and abs(b) >= 0.01 and (a > 0) != (b > 0))
out.append(f"[nesting] 5x pairs {len(names5)}: shared with 2x {len(diffs)} (delta at 5x {shared:.6f}; same pairs at 2x {sum(v[0] for v in d2.values()):.6f}), "
           f"added {len(names5) - len(diffs)} (delta {added:.6f}); shared-pair |d5 - d2| max {diffs[0][0]:.4f}, sign flips among non-neutral {sign_flips}")

# read placement for strongest contributors
index = defaultdict(list)
for h, rec in HAPS.items():
    s = fa[rec]
    for i in range(0, len(s) - SEED + 1):
        index[s[i:i + SEED]].append((h, i))


def place(read):
    best = {}
    for strand, seq in (("+", read), ("-", rc(read))):
        cand = set()
        for j in range(0, len(seq) - SEED + 1, 8):
            for h, i in index.get(seq[j:j + SEED], ()):
                cand.add((h, i - j))
        for h, off in cand:
            s = fa[HAPS[h]]
            if off < 0 or off + len(seq) > len(s):
                continue
            mm = sum(1 for a, b in zip(seq, s[off:off + len(seq)]) if a != b)
            if mm <= MAXMM:
                best.setdefault(h, []).append((mm, off, strand))
    return {h: sorted(v) for h, v in best.items()}


def describe(fid, d, reads, mates, label):
    n1, s1, q1 = reads[fid]
    _, s2_, q2 = mates[fid]
    p1, p2 = place(s1), place(s2_)
    rows = []
    totals = {}
    for h, rec in HAPS.items():
        a, b = p1.get(h, []), p2.get(h, [])
        bm1 = a[0][0] if a else None
        bm2 = b[0][0] if b else None
        amb1 = sum(1 for x in a if x[0] == bm1) if a else 0
        amb2 = sum(1 for x in b if x[0] == bm2) if b else 0
        geo = None
        for x in a:
            for y in b:
                if x[2] != y[2]:
                    left, right = min(x[1], y[1]), max(x[1] + len(s1), y[1] + len(s2_))
                    ins = right - left
                    if abs(ins - MEAN) <= 4 * SD and (geo is None or x[0] + y[0] < geo[0]):
                        geo = (x[0] + y[0], ins, left)
        cb = coords[rec]
        inblk = geo is not None and geo[2] < cb[1] and geo[2] + geo[1] > cb[0]
        dk = disc.get(h, set())
        over = geo is not None and any(cb[0] + p >= geo[2] - K and cb[0] + p < geo[2] + geo[1] for p in dk)
        totals[h] = geo[0] if geo else None
        rows.append(f"{h}:mm {bm1}/{bm2} amb {amb1}/{amb2} proper {('mm ' + str(geo[0]) + ' insert ' + str(geo[1])) if geo else 'none'}"
                    f"{' ctx0' if inblk else ''}{' disc' if over else ''}")
    tA = min([x for x in (totals['A1'], totals['A2']) if x is not None], default=None)
    tB = min([x for x in (totals['B1'], totals['B2']) if x is not None], default=None)
    flag = ""
    if abs(d[fid][0]) >= 0.5 and tA is not None and tB is not None:
        if (d[fid][0] < 0 and tA < tB) or (d[fid][0] > 0 and tB < tA):
            flag = "CANDIDATE_INCONSISTENCY"
    quals = sorted(set(q1 + q2))
    return (f"[{label}] pair {n1} (fragment {fid}) delta {d[fid][0]:.4f} class {d[fid][1]} | " + " | ".join(rows) +
            f" | best proper mm A {tA} B {tB} | quals {''.join(quals)} {flag}"), flag


flags = []
top2 = sorted(d2, key=lambda f: -abs(d2[f][0]))[:10]
for f in top2:
    text, fl = describe(f, d2, reads2, mates2, "2x top")
    out.append(text)
    flags += [(f, d2[f][0])] if fl else []
added_ids = [i for i, nm in enumerate(names5) if nm not in idx2]
top5 = sorted(added_ids, key=lambda f: -abs(d5[f][0]))[:10]
for f in top5:
    text, fl = describe(f, d5, reads5, mates5, "5x added top")
    out.append(text)
# scorer origin counts for the strongest contributors (placement ambiguity as the scorer saw it)
for label, d in (("2x", s2), ("5x", s5)):
    tops = table([p for p in __import__("glob").glob(d + "/score.decomp.*.top_fragments.tsv")][0])[:10]
    out.append(f"[{label}] scorer origins (A1/A2/B1/B2 total, in-block) for top 10: " + "; ".join(
        f"f{r['fragment']} {r['A1_origins']}/{r['A2_origins']}/{r['B1_origins']}/{r['B2_origins']} in {r['A1_in_block']}/{r['A2_in_block']}/{r['B1_in_block']}/{r['B2_in_block']}"
        for r in tops))
flag_sum = sum(v for _, v in flags)
out.append(f"FLAGS {len(flags)} among the 2x top 10 (combined delta {flag_sum:.4f}; net 2x deficit {sum(v[0] for v in d2.values()):.4f})")
open(prefix + ".txt", "w").write("\n".join(out) + "\n")
print("\n".join(out))
