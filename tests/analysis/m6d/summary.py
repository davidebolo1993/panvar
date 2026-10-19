#!/usr/bin/env python3
"""M6-D summary (synthetic; truth is the declared fixture truth; frozen with DESIGN.md).

From each run's compact state table (absolute alleles, likelihood, both log priors) this script independently recomputes,
per scaffold and prior, log marginal = log sum_theta P(theta|G) L(theta), the posterior-mean log likelihood and
KL(posterior || prior), and checks:
  C1  the recomputation equals the executable's scaffold table (|d| <= 1e-6)
  C2  the original-prior marginals reproduce M6-C's recorded marginals for the same fixture/seed (|d| <= 1e-6)
  C3  log m = E_q log L - KL holds for every scaffold and prior (|d| <= 1e-9)
Per run and prior (uniform weights over the nine allocations, unchanged from M6-C; they induce total-copy prior masses
2/9, 3/9, 4/9 for totals 3, 4, 5):
  allocation posterior; total-copy posterior; per-total evidence = log mean of the allocation marginals within a total
  (the posterior odds between totals divided by the induced prior odds);
  allocation class as M6-C: CORRECT (P_truth >= 0.95), UNCERTAIN (truth in a tied top set, ties on log marginals within
  1e-6), WRONG (another allocation >= 0.95), SPLIT;
  dosage class: DOSAGE_CORRECT (truth total >= 0.95), DOSAGE_WRONG (another total >= 0.95), DOSAGE_SPLIT;
  allocation within the truth's total: conditional posterior of the truth allocation.
Attribution (F2 and F4): for the truth against every other allocation, the log-marginal difference under each prior,
decomposed exactly into d(E_q log L) - d(KL); the centring effect = d_common - d_original.
Pooled-content equality is never used to declare two allocations equivalent; ties are read from the log marginals only.

  summary.py <runs dir> <m6c runs dir> <ledger> <out prefix>
"""
import csv
import glob
import math
import os
import sys
from collections import defaultdict

runs, m6c, ledger, prefix = sys.argv[1:5]
T = lambda p: list(csv.DictReader(open(p), delimiter="\t"))
KV = lambda p: {l.split("\t")[0]: l.rstrip("\n").split("\t")[1] for l in open(p) if "\t" in l}
PRIORS = ["original", "common_content"]
PRIOR_COL = {"original": "log_prior_original", "common_content": "log_prior_common_content"}
TOTAL_PRIOR = {3: 2 / 9, 4: 3 / 9, 5: 4 / 9}
total = lambda lab: sum(int(x) for x in lab.split("+"))


def lse(v):
    m = max(v)
    return m if m == -math.inf else m + math.log(sum(math.exp(x - m) for x in v))


checks, srows, rrows, arows = [], [], [], []
for d in sorted(glob.glob(f"{runs}/run.*")):
    if not os.path.exists(f"{d}/d.run.tsv"):
        continue
    info = KV(f"{d}/d.run.tsv")
    fx, seed, truth = info["fixture"], info["seed"], info["truth"]
    states = defaultdict(list)
    for r in T(f"{d}/d.states.tsv"):
        states[r["scaffold"]].append(r)
    exe = {(r["scaffold"], r["prior"]): r for r in T(f"{d}/d.scaffolds.tsv")}
    old = {r["scaffold"]: float(r["log_marginal"]) for r in T(f"{m6c}/run.{fx}.{seed}/n.scaffolds.tsv")}
    res = {}
    for lab, rows in states.items():
        ll = [float(r["log_likelihood"]) for r in rows]
        for pr in PRIORS:
            lp = [float(r[PRIOR_COL[pr]]) for r in rows]
            terms = [a + b for a, b in zip(lp, ll)]
            # Amendment 1: explicit normalisation and centred log q (not through the rounding of log m)
            top = max(terms)
            w = [math.exp(t - top) for t in terms]
            W = math.fsum(w)
            m = top + math.log(W)
            q = [wi / W for wi in w]
            el = math.fsum(qi * li for qi, li in zip(q, ll))
            kl = math.fsum(qi * ((t - top) - math.log(W) - p) for qi, t, p in zip(q, terms, lp) if qi > 0)
            alt_mean = sum(qi * int(r["alt_count"]) for qi, r in zip(q, rows))
            mass = math.fsum(math.exp(p) for p in lp)
            res[(lab, pr)] = dict(m=m, el=el, kl=kl, alt=alt_mean, states=len(rows), base=exe[(lab, pr)]["base_allele"])
            checks.append((fx, seed, lab, pr, "C1_recompute", m - float(exe[(lab, pr)]["log_marginal"]), 1e-6))
            checks.append((fx, seed, lab, pr, "C3_kl_identity", m - (el - kl), 1e-9))
            checks.append((fx, seed, lab, pr, "C0_prior_mass", mass - 1.0, 1e-12))
            if pr == "original":
                checks.append((fx, seed, lab, pr, "C2_reproduces_m6c", m - old[lab], 1e-6))
    labs = sorted(states, key=lambda s: (total(s), s))
    for pr in PRIORS:
        logs = {lab: res[(lab, pr)]["m"] for lab in labs}
        z = lse(list(logs.values()))
        post = {lab: math.exp(logs[lab] - z) for lab in labs}
        tot_post = defaultdict(float)
        for lab in labs:
            tot_post[total(lab)] += post[lab]
        evid = {t: lse([logs[l] for l in labs if total(l) == t]) - math.log(sum(1 for l in labs if total(l) == t)) for t in TOTAL_PRIOR}
        top = max(logs, key=logs.get)
        tied = [k for k in labs if abs(logs[k] - logs[top]) <= 1e-6]
        if post[truth] >= 0.95:
            cls = "CORRECT"
        elif truth in tied and len(tied) > 1:
            cls = "UNCERTAIN"
        elif post[top] >= 0.95:
            cls = "WRONG"
        else:
            cls = "SPLIT"
        tt = total(truth)
        top_total = max(tot_post, key=tot_post.get)
        dcls = "DOSAGE_CORRECT" if tot_post[tt] >= 0.95 else ("DOSAGE_WRONG" if tot_post[top_total] >= 0.95 else "DOSAGE_SPLIT")
        within = post[truth] / tot_post[tt] if tot_post[tt] > 0 else float("nan")
        same_total = [logs[l] for l in labs if total(l) == tt]
        rrows.append(dict(fixture=fx, seed=seed, truth=truth, prior=pr, allocation_class=cls, P_truth=f"{post[truth]:.6g}", top=top,
                          P_top=f"{post[top]:.6g}", tied_top="/".join(tied), dosage_class=dcls,
                          P_total3=f"{tot_post[3]:.6g}", P_total4=f"{tot_post[4]:.6g}", P_total5=f"{tot_post[5]:.6g}",
                          evidence_total3=f"{evid[3]:.6f}", evidence_total4=f"{evid[4]:.6f}", evidence_total5=f"{evid[5]:.6f}",
                          P_truth_within_truth_total=f"{within:.6g}", spread_log_m_within_truth_total=f"{max(same_total) - min(same_total):.6g}"))
        for lab in labs:
            r = res[(lab, pr)]
            srows.append(dict(fixture=fx, seed=seed, scaffold=lab, base_allele=r["base"], prior=pr, nuisance_states=r["states"],
                              log_marginal=f"{r['m']:.10f}", posterior_mean_log_likelihood=f"{r['el']:.10f}", kl=f"{r['kl']:.10f}",
                              posterior_mean_alt_count=f"{r['alt']:.6f}", posterior=f"{post[lab]:.6g}", log_m_minus_best=f"{logs[lab] - logs[top]:.6f}"))
    for lab in labs:
        if lab == truth:
            continue
        d = {}
        for pr in PRIORS:
            a, b = res[(truth, pr)], res[(lab, pr)]
            d[pr] = (a["m"] - b["m"], a["el"] - b["el"], a["kl"] - b["kl"])
        arows.append(dict(fixture=fx, seed=seed, truth=truth, versus=lab, versus_base=res[(lab, "original")]["base"],
                          d_log_m_original=f"{d['original'][0]:.6f}", d_ElogL_original=f"{d['original'][1]:.6f}", d_KL_original=f"{d['original'][2]:.6f}",
                          d_log_m_common=f"{d['common_content'][0]:.6f}", d_ElogL_common=f"{d['common_content'][1]:.6f}", d_KL_common=f"{d['common_content'][2]:.6f}",
                          centring_effect=f"{d['common_content'][0] - d['original'][0]:.6f}"))


def write(name, rows):
    if not rows:
        return
    with open(f"{prefix}.{name}.tsv", "w") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


write("scaffolds", srows)
write("runs", rrows)
write("attribution", arows)
failed = [c for c in checks if not abs(c[5]) <= c[6]]
with open(f"{prefix}.checks.tsv", "w") as f:
    f.write("check\tn\tmax_abs\ttolerance\tresult\n")
    for name in sorted({c[4] for c in checks}):
        cs = [c for c in checks if c[4] == name]
        mx = max(abs(c[5]) for c in cs)
        f.write(f"{name}\t{len(cs)}\t{mx:.3g}\t{cs[0][6]:g}\t{'PASS' if mx <= cs[0][6] else 'FAIL'}\n")
led = T(ledger)
with open(f"{prefix}.cost.tsv", "w") as f:
    f.write("unit\twall_seconds\tpeak_rss_gib\texit\n")
    for r in led:
        f.write(f"{r['unit']}\t{r['wall_seconds']}\t{r['peak_rss_gib']}\t{r['exit']}\n")
print(open(f"{prefix}.checks.tsv").read())
for r in rrows:
    print("\t".join(str(v) for v in r.values()))
sys.exit(1 if failed else 0)
