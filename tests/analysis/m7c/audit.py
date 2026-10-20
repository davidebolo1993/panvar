#!/usr/bin/env python3
"""M7-C read-only audit of the M7-B.1 fixed-pool comparison (no genotyping; M7-B / M7-B.1 packages are read only).

  A  attempts: decode every concatenated Brotli stream of sol_ext..csv.br (large window; complete consumption); gate: every
     declared genotype has exactly attempts 1..20 at stage 2; reconstruct res.json lik_mean and lik_sd; weighted score check;
     negative controls (truncated stream must fail to decode; a missing or duplicated attempt must fail the gate)
  B  GC: reproduce Locityper's cumulative-tail boundaries (bg/depth.rs blur_boundary_values) from the per-GC background
     window counts and compare with the logged values; classify target windows (weight >= 0.001) as in sparse bins (< 100
     windows), outside the cumulative-tail boundaries (blurred), or outside the observed background GC range
  C  score terms: mean aln_lik / depth_lik / lik over the 20 stage-2 attempts; frozen contrasts top vs best correct-class
     candidate and top vs the exact caller candidate, weighted by the actual coefficients (1 - s, 1 + s), s = 0.85
  D  variance semantics: res.json lik_sd versus the recovered attempt variance

  audit.py <out dir>
"""
import csv
import gzip
import json
import math
import os
import re
import statistics
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from brotli_members import decode_members  # noqa: E402

out = sys.argv[1]
B1 = "/Users/davide.bolognini/panvar/tests/results/m7b1_fixed_pool_20261020/runs"
B0 = "/Users/davide.bolognini/panvar/tests/results/m7b_fixed_pool_20261020/runs"
SKEW = 0.85
WA, WD = 1 - SKEW, 1 + SKEW
LN10 = math.log(10)
T = lambda p: list(csv.DictReader(open(p), delimiter="\t"))
canon = lambda g: ",".join(sorted(g.split(",")))
LOCUS = {"c4": "C4", "cyp": "CYP2D6"}
TRUTH = {"c4_loo_na19650_10": "2+2", "cyp_loo_na19650_5": "2+3", "c4_lzo_hg00171_5": "2+2", "cyp_lzo_hg02769_5": "2+2"}
mapping = T(f"{B0}/pools/mapping.tsv")


def write(name, rows):
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(f"{out}/{name}.tsv", "w") as f:
        w = csv.DictWriter(f, fieldnames=keys, delimiter="\t", lineterminator="\n", restval="")
        w.writeheader()
        w.writerows(rows)


def gate(rows, declared):
    att = defaultdict(list)
    for r in rows:
        if r["stage"] == "2":
            att[canon(r["genotype"])].append(int(r["attempt"]))
    return set(att) == declared and all(sorted(a) == list(range(1, 21)) for a in att.values())


units = sorted(u for u in os.listdir(B1) if u.startswith("loc.") and os.path.isdir(f"{B1}/{u}/gt"))
A, Crows, Drows, attempts_by_run = [], [], [], {}
for u in units:
    cell, pool = u.split(".")[1], u.split(".")[2]
    gdir = f"{B1}/{u}/gt/loci/{LOCUS[cell.split('_')[0]]}"
    data = open(f"{gdir}/sol_ext..csv.br", "rb").read()
    txt, members = decode_members(data)
    rows = list(csv.DictReader(txt.splitlines(), delimiter="\t"))
    declared = {l.split("\t")[1] for l in open(f"{B0}/pools/{cell}.priors.{pool}.tsv")}
    res = json.load(gzip.open(f"{gdir}/res.json.gz", "rt"))
    opts = {canon(o["genotype"]): o for o in res["options"]}
    st2 = defaultdict(list)
    for r in rows:
        if r["stage"] == "2":
            st2[canon(r["genotype"])].append(r)
    attempts_by_run[u] = st2
    lik = {g: [float(r["lik"]) for r in rs] for g, rs in st2.items()}
    mean_err = max(abs(statistics.mean(v) - opts[g]["lik_mean"]) for g, v in lik.items())
    var_s_err = max(abs(statistics.variance(v) * LN10 - opts[g]["lik_sd"]) for g, v in lik.items())
    var_p_err = max(abs(statistics.pvariance(v) * LN10 - opts[g]["lik_sd"]) for g, v in lik.items())
    w_err = max(abs(WA * float(r["aln_lik"]) + WD * float(r["depth_lik"]) - float(r["lik"])) for rs in st2.values() for r in rs)
    # negative controls
    try:
        decode_members(data[:-7])
        trunc = "NOT DETECTED"
    except Exception:
        trunc = "detected"
    try:
        cut = members[0][0] + members[1][0] // 2
        decode_members(data[:cut])
        trunc_mid = "NOT DETECTED"
    except Exception:
        trunc_mid = "detected"
    g0 = sorted(st2)[0]
    missing = [r for r in rows if not (r["stage"] == "2" and canon(r["genotype"]) == g0 and r["attempt"] == "20")]
    dup = rows + [next(r for r in rows if r["stage"] == "2" and canon(r["genotype"]) == g0)]
    A.append(dict(run=u, streams=len(members), bytes=len(data), all_bytes_consumed=sum(n for n, _ in members) == len(data),
                  stage2_genotypes=len(st2), declared=len(declared), attempt_gate="PASS" if gate(rows, declared) else "FAIL",
                  neg_truncated_tail=trunc, neg_truncated_mid_member=trunc_mid,
                  neg_missing_attempt="gate FAIL (detected)" if not gate(missing, declared) else "NOT DETECTED",
                  neg_duplicate_attempt="gate FAIL (detected)" if not gate(dup, declared) else "NOT DETECTED",
                  max_err_mean_vs_res=f"{mean_err:.2e}", max_err_lik_sd_vs_sample_var_x_ln10=f"{var_s_err:.2e}",
                  max_err_lik_sd_vs_population_var_x_ln10=f"{var_p_err:.2e}", max_err_weighted_score=f"{w_err:.2e}"))
    # D: variance semantics per genotype (top 3 by mean)
    for g in sorted(lik, key=lambda g: -statistics.mean(lik[g]))[:3]:
        v = statistics.variance(lik[g])
        Drows.append(dict(run=u, genotype=g, res_lik_mean=opts[g]["lik_mean"], recovered_mean=statistics.mean(lik[g]), res_lik_sd=opts[g]["lik_sd"],
                          recovered_sample_variance=v, variance_times_ln10=v * LN10, attempt_sd=math.sqrt(v), standard_error_of_mean=math.sqrt(v / 20),
                          log10_of_variance_for_contrast=math.log10(v)))
    # C: score-term contrasts on D pools
    if pool != "D":
        continue
    dmap = defaultdict(list)
    for m in mapping:
        if m["cell"] == cell and m["tier"] == "D":
            dmap[m["canonical_diplotype"]].append(m)
    mean = {g: {k: statistics.mean(float(r[k]) for r in rs) for k in ("aln_lik", "depth_lik", "lik")} for g, rs in st2.items()}
    top = max(mean, key=lambda g: mean[g]["lik"])
    correct = [g for g in mean if dmap[g][0]["unit_class"] == TRUTH[cell]]
    best_correct = max(correct, key=lambda g: mean[g]["lik"])
    caller = next(g for g in dmap if any("reported" in m["labels"].split("|") for m in dmap[g]))
    lab = lambda g: "||".join(m["labels"] for m in dmap[g])
    rank = {g: i + 1 for i, g in enumerate(sorted(mean, key=lambda g: -mean[g]["lik"]))}
    for name, g in (("best correct-class candidate", best_correct), ("exact caller candidate", caller)):
        d = {k: mean[g][k] - mean[top][k] for k in mean[top]}
        Crows.append(dict(run=u, contrast=f"{name} minus top", candidate=lab(g)[:60], candidate_class=dmap[g][0]["unit_class"], candidate_rank=rank[g],
                          top=lab(top)[:60], top_class=dmap[top][0]["unit_class"], d_total=d["lik"], d_alignment_weighted=WA * d["aln_lik"],
                          d_depth_weighted=WD * d["depth_lik"], d_alignment_raw=d["aln_lik"], d_depth_raw=d["depth_lik"],
                          reconstruction_error=WA * d["aln_lik"] + WD * d["depth_lik"] - d["lik"]))

# B: GC
Brows = []
for cell in TRUTH:
    sup = {int(r["gc"]): int(r["nwindows"]) for r in T(f"{B1}/preproc.{cell}/gc_support.tsv")}
    n = 101
    counts = [sup.get(g, 0) for g in range(n)]
    ends = [sum(counts[:g + 1]) for g in range(n)]
    starts = [sum(counts[:g]) for g in range(n)]
    m = ends[-1]
    left_ix = next((i for i in range(n) if ends[i] >= 100), n)
    p = next((k for k, i in enumerate(reversed(range(n))) if m - starts[i] >= 100), n)
    right_ix = n - p
    logged = re.findall(r"GC-content < (\d+) or > (\d+)", re.sub(r"\x1b\[[0-9;]*m", "", open(f"{B1}/preproc.{cell}/commands.log").read()))
    observed = [g for g in range(n) if counts[g] > 0]
    for u in [u for u in units if u.split(".")[1] == cell]:
        gdir = f"{B1}/{u}/gt/loci/{LOCUS[cell.split('_')[0]]}"
        win = [l.rstrip("\n").split("\t") for l in gzip.open(f"{gdir}/windows.bed.gz", "rt") if not l.startswith("#")]
        used = [int(w[3]) for w in win if float(w[7]) >= 0.001]
        N = len(used)
        sparse = sum(1 for g in used if counts[g] < 100)
        blurred = sum(1 for g in used if g < left_ix or g > right_ix)
        outside = sum(1 for g in used if g < min(observed) or g > max(observed))
        Brows.append(dict(run=u, background_windows=m, logged_boundaries=f"<{logged[0][0]} or >{logged[0][1]}" if logged else "not logged",
                          reproduced_boundaries=f"<{left_ix} or >{right_ix}", reproduced_equals_logged=bool(logged) and (str(left_ix), str(right_ix)) == logged[0],
                          observed_background_gc=f"{min(observed)}-{max(observed)}", target_windows=N,
                          sparse_bin_windows=sparse, sparse_bin_fraction=f"{sparse / N:.4f}",
                          outside_cumulative_tail_boundaries=blurred, outside_boundaries_fraction=f"{blurred / N:.4f}",
                          outside_observed_range=outside, outside_observed_fraction=f"{outside / N:.4f}"))
os.makedirs(out, exist_ok=True)
write("attempts", A)
write("gc", Brows)
write("score_terms", Crows)
write("variance", Drows)
for name in ("attempts", "gc", "score_terms"):
    print(open(f"{out}/{name}.tsv").read())
ok = all(r["attempt_gate"] == "PASS" and r["all_bytes_consumed"] and "NOT" not in (r["neg_truncated_tail"] + r["neg_truncated_mid_member"] + r["neg_missing_attempt"] + r["neg_duplicate_attempt"])
         for r in A) and all(r["reproduced_equals_logged"] for r in Brows)
sys.exit(0 if ok else 1)
