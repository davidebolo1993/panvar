#!/usr/bin/env python3
"""M3-H prior attribution (diagnostic; no prior change; the current copying model only, r = 0.01).

The copying prior (src/ls_hmm.cpp): start uniform over the n curated panel templates; between blocks T = (1-r)I + (r/n)11^T,
where r is a donor REDRAW probability (a copying-model parameter, not a biological recombination rate); a template emits
only the allele it carries. log P(haplotype) sums every compatible template history (forward recursion).

GATE (before anything is reported): for each candidate, log P(h1) + log P(h2) + log 2 (distinct haplotypes) from this
independent recursion must equal the caller-code sequence-level pair prior (score-set's `prior`, which also sums every
allele path spelling the same bytes) within 1e-9. Equality means no alternative spelling contributes for these
haplotypes, so the per-history attribution below decomposes the sequence-level prior itself.

Attribution, per homologue (the two candidates differ at block 0 only):
  P(h) = (1/n) sum_{t in C(a0)} beta0(t),   beta0(t) = (1-r) e1(t) beta1(t)  +  (r/n) S1,   S1 = sum_u e1(u) beta1(u)
  C(a0): templates carrying the block-0 allele (donor compatibility);
  'continue' = histories whose block-0 donor also copies block 1 (no redraw at the ctx:0/sv:1 boundary);
  'redraw'   = histories that redraw at that boundary (|C| r/n S1; identical per-template, so it scales with |C|).
Also: posterior expected number of donor changes (redraw to a different template) over the whole haplotype, and the
block-0 donor posterior.

  prior_attrib.py <panel-alleles table> <score-set scores table> <out prefix>
"""
import csv
import math
import sys

panel_path, scores_path, prefix = sys.argv[1], sys.argv[2], sys.argv[3]
R = 0.01
rows = []
for line in open(panel_path):
    if line.startswith("#") or line.startswith("template\t"):
        continue
    name, alleles = line.rstrip("\n").split("\t")
    if name.startswith("MOSAIC#"):   # curated out of the panel (mosaic_exclude.tsv)
        continue
    rows.append((name, list(map(int, alleles.split(",")))))
n = len(rows)
B = len(rows[0][1])
L = lambda xs: max(xs) + math.log(sum(math.exp(x - max(xs)) for x in xs)) if xs and max(xs) > -math.inf else -math.inf
ls, lw = math.log1p(-R), math.log(R) - math.log(n)


def emit(b, a):
    return [0.0 if rows[t][1][b] == a else -math.inf for t in range(n)]


def forward(path):
    al = [[(-math.log(n) + e) for e in emit(0, path[0])]]
    for b in range(1, B):
        prev = al[-1]
        tot = L(prev)
        e = emit(b, path[b])
        al.append([(L([prev[u] + ls, tot + lw]) + e[u]) if e[u] == 0.0 else -math.inf for u in range(n)])
    return al


def backward(path):
    be = [[0.0] * n for _ in range(B)]
    for b in range(B - 2, -1, -1):
        e = emit(b + 1, path[b + 1])
        nxt = [be[b + 1][u] + e[u] for u in range(n)]
        s = L(nxt)
        be[b] = [L([ls + nxt[t], lw + s]) for t in range(n)]
    return be


def analyse(path):
    al, be = forward(path), backward(path)
    lp = L(al[-1])
    # expected donor changes: sum_b sum_{t != u} alpha_b(t) (r/n) e(u) beta(u) / P
    changes = 0.0
    for b in range(B - 1):
        e = emit(b + 1, path[b + 1])
        tot = L(al[b])
        for u in range(n):
            if e[u] != 0.0:
                continue
            into = L([al[b][t] for t in range(n) if t != u])
            if into > -math.inf:
                changes += math.exp(into + lw + be[b + 1][u] - lp)
    # boundary-0 split
    e1 = emit(1, path[1])
    S1 = L([e1[u] + be[1][u] for u in range(n)])
    comp = [t for t in range(n) if rows[t][1][0] == path[0]]
    cont = L([-math.log(n) + ls + e1[t] + be[1][t] for t in comp])
    redraw = (-math.log(n) + math.log(len(comp)) + lw + S1) if comp else -math.inf
    post0 = sorted(((math.exp(al[0][t] + be[0][t] - lp), rows[t][0]) for t in comp), reverse=True)
    cont_names = [rows[t][0] for t in comp if e1[t] == 0.0]
    return dict(lp=lp, changes=changes, comp=len(comp), comp_names=[rows[t][0] for t in comp], cont=cont, redraw=redraw,
                cont_names=cont_names, post0=post0, check=L([cont, redraw]))


sc = {r["label"]: r for r in csv.DictReader(open(scores_path), delimiter="\t")}
out, res = [], {}
gate_ok = True
for lab in ("R_c0", "R_t0"):
    h1, h2 = list(map(int, sc[lab]["hap1_alleles"].split(","))), list(map(int, sc[lab]["hap2_alleles"].split(",")))
    a, b = analyse(h1), analyse(h2)
    pair = a["lp"] + b["lp"] + (math.log(2.0) if h1 != h2 else 0.0)
    diff = pair - float(sc[lab]["prior"])
    gate_ok &= abs(diff) <= 1e-9 and abs(a["check"] - a["lp"]) <= 1e-9 and abs(b["check"] - b["lp"]) <= 1e-9
    res[lab] = (a, b)
    out.append(f"GATE {lab}: independent pair prior {pair:.12f} vs score-set sequence-level prior {sc[lab]['prior']} (diff {diff:.2e}); "
               f"boundary split reconstructs each haplotype: {abs(a['check'] - a['lp']) < 1e-9 and abs(b['check'] - b['lp']) < 1e-9}")
out.append(f"GATE {'PASS' if gate_ok else 'FAIL'} (n = {n} curated templates, {B} blocks, r = {R})")
lines = ["candidate\thomologue\tblock0_allele\tlog_prior\tcompatible_templates\tcontinue_into_block1_templates\tlog_continue\tlog_redraw\t"
         "continue_share\texpected_donor_changes\tblock0_donor_posterior_top"]
for lab, (a, b) in res.items():
    for k, h in ((1, a), (2, b)):
        allele = sc[lab][f"hap{k}_alleles"].split(",")[0]
        share = math.exp(h["cont"] - h["lp"] - math.log(1) + math.log(1)) if h["cont"] > -math.inf else 0.0
        share = math.exp(h["cont"] - L([h["cont"], h["redraw"]])) if h["cont"] > -math.inf else 0.0
        lines.append("\t".join(map(str, [lab, k, allele, f"{h['lp']:.6f}", f"{h['comp']} ({','.join(h['comp_names'])})" if h["comp"] <= 6 else h["comp"],
                                         f"{len(h['cont_names'])} ({','.join(h['cont_names'][:6])}{',...' if len(h['cont_names']) > 6 else ''})",
                                         f"{h['cont']:.6f}", f"{h['redraw']:.6f}", f"{share:.4f}", f"{h['changes']:.4f}",
                                         "; ".join(f"{nm} {p:.3f}" for p, nm in h["post0"][:3])])))
if gate_ok:
    for k in (0, 1):
        dc = res["R_t0"][k]["lp"] - res["R_c0"][k]["lp"]
        out.append(f"homologue {k + 1}: truth allele {sc['R_t0'][f'hap{k + 1}_alleles'].split(',')[0]} vs called "
                   f"{sc['R_c0'][f'hap{k + 1}_alleles'].split(',')[0]}: delta log prior {dc:.6f}; compatible templates "
                   f"{res['R_t0'][k]['comp']} vs {res['R_c0'][k]['comp']}; continuing into block 1 {len(res['R_t0'][k]['cont_names'])} vs "
                   f"{len(res['R_c0'][k]['cont_names'])}; expected donor changes {res['R_t0'][k]['changes']:.3f} vs {res['R_c0'][k]['changes']:.3f}")
    out.append(f"pair delta prior (truth - called) {float(sc['R_t0']['prior']) - float(sc['R_c0']['prior']):.6f} = homologue sum "
               f"{sum(res['R_t0'][k]['lp'] - res['R_c0'][k]['lp'] for k in (0, 1)):.6f}")
open(prefix + ".haplotypes.tsv", "w").write("\n".join(lines) + "\n")
open(prefix + ".summary.txt", "w").write("\n".join(out) + "\n")
print("\n".join(out))
print("\n".join(lines))
sys.exit(0 if gate_ok else 1)
