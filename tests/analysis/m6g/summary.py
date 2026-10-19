#!/usr/bin/env python3
"""M6-G summary (frozen with DESIGN.md). Structural posterior = uniform weights over the six allocations (explicit,
identical for RES and PANEL; induced totals 2/3/4 = 1/6, 2/6, 3/6) x the marginal likelihood.

Per fixture and method: allocation posterior, total-copy posterior, and for the truth against every other allocation
d log m = d E_q[log L] - d KL (read evidence vs prior/complexity). Negative control: the per-read independently
marginalised score's top allocation, reported beside the correct one.

Frozen criteria:
  A_missing_variation (truth 1+2, total 3; novel alleles at both admitted sites in all 3 copies)
     PANEL premise met: PANEL's P(total 3) < 0.5
     RES: TOTAL_CORRECT P(total 3) >= 0.95 | MANUFACTURED_WRONG_TOTAL some other total >= 0.95 | TOTAL_UNRESOLVED
  B_unbalanced_allocation (truth 3+1): CORRECT P(3+1) >= 0.95 | BALANCED_PREFERENCE 2+2 is the top allocation | OTHER
  C_uninformative_allocation (truth 2+2): RETAINED max conditional allocation posterior within total 4 < 0.95 | OVERCONFIDENT

  summary.py <runs dir> <out prefix>
"""
import csv
import glob
import math
import os
import sys
from collections import defaultdict

runs, prefix = sys.argv[1:3]
T = lambda p: list(csv.DictReader(open(p), delimiter="\t"))
KV = lambda p: {l.split("\t")[0]: l.rstrip("\n").split("\t")[1] for l in open(p) if "\t" in l}
TRUTH = {"A_missing_variation": "1+2", "B_unbalanced_allocation": "3+1", "C_uninformative_allocation": "2+2"}
tot = lambda s: sum(int(x) for x in s.split("+"))


def lse(v):
    m = max(v)
    return m + math.log(sum(math.exp(x - m) for x in v))


rows, verdicts, contrasts = [], [], []
for d in sorted(glob.glob(f"{runs}/run.*")):
    if not os.path.exists(f"{d}/g.run.tsv"):
        continue
    info = KV(f"{d}/g.run.tsv")
    fx = info["fixture"]
    st = T(f"{d}/g.structures.tsv")
    post = {}
    for method in ("RES", "PANEL"):
        r = {x["structure"]: x for x in st if x["method"] == method}
        lm = {k: float(v["log_marginal"]) for k, v in r.items()}
        z = lse(list(lm.values()))
        p = {k: math.exp(v - z) for k, v in lm.items()}
        nc = {k: float(v["negative_control_per_read"]) for k, v in r.items()}
        pt = defaultdict(float)
        for k, v in p.items():
            pt[tot(k)] += v
        post[method] = (p, pt)
        truth = TRUTH[fx]
        for k in sorted(r, key=lambda k: (tot(k), k)):
            rows.append(dict(fixture=fx, method=method, structure=k, total=tot(k), posterior=f"{p[k]:.6g}", total_posterior=f"{pt[tot(k)]:.6g}",
                             log_marginal=r[k]["log_marginal"], posterior_mean_log_likelihood=r[k]["posterior_mean_log_likelihood"], kl=r[k]["kl"],
                             posterior_mean_novel_site_copies=r[k]["posterior_mean_novel_site_copies"], negative_control=r[k]["negative_control_per_read"]))
            if k != truth:
                contrasts.append(dict(fixture=fx, method=method, truth=truth, versus=k,
                                      d_log_m=f"{lm[truth] - lm[k]:.6f}",
                                      d_ElogL=f"{float(r[truth]['posterior_mean_log_likelihood']) - float(r[k]['posterior_mean_log_likelihood']):.6f}",
                                      d_KL=f"{float(r[truth]['kl']) - float(r[k]['kl']):.6f}",
                                      d_negative_control=f"{nc[truth] - nc[k]:.6f}"))
        top = max(p, key=p.get)
        nctop = max(nc, key=nc.get)
        if fx == "A_missing_variation":
            if method == "PANEL":
                v = "PREMISE_MET" if pt[3] < 0.5 else "PREMISE_NOT_MET"
            else:
                v = "TOTAL_CORRECT" if pt[3] >= 0.95 else ("MANUFACTURED_WRONG_TOTAL" if any(pt[t] >= 0.95 for t in pt if t != 3) else "TOTAL_UNRESOLVED")
        elif fx == "B_unbalanced_allocation":
            v = "CORRECT" if p["3+1"] >= 0.95 else ("BALANCED_PREFERENCE" if top == "2+2" else "OTHER")
        else:
            within = [p[k] / pt[4] for k in p if tot(k) == 4] if pt[4] > 0 else [1.0]
            v = "RETAINED" if max(within) < 0.95 else "OVERCONFIDENT"
        verdicts.append(dict(fixture=fx, method=method, verdict=v, top_allocation=top, P_top=f"{p[top]:.6g}", P_truth=f"{p[truth]:.6g}",
                             P_total2=f"{pt[2]:.6g}", P_total3=f"{pt[3]:.6g}", P_total4=f"{pt[4]:.6g}", negative_control_top=nctop,
                             truth_cumulative_residual_prior_cost=info["cumulative_residual_prior_cost"], truth_novel_site_copies=info["truth_novel_site_copies"],
                             seconds=info["seconds"], peak_rss_gib=info["peak_rss_gib"], candidates=info["candidates_scored"], unique_haplotypes=info["unique_haplotypes"],
                             fragments=info["fragments"]))


def write(name, rs):
    if rs:
        with open(f"{prefix}.{name}.tsv", "w") as f:
            w = csv.DictWriter(f, fieldnames=list(rs[0]), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(rs)


write("structures", rows)
write("contrasts", contrasts)
write("verdicts", verdicts)
for v in verdicts:
    print("\t".join(str(x) for x in v.values()))
