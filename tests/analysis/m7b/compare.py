#!/usr/bin/env python3
"""M7-B comparison (frozen with DESIGN.md). Rankings and gaps are compared WITHIN each method only.

Truth labels (evaluation only; from records): C4 NA19650 truth 2+2 (M1-B alts header), CYP NA19650 truth 2+3 (M1-B alts
header), controls 2+2 with the caller correct (M6-A summary, caller_structure CORRECT, caller units 2+2).

  compare.py <unit dir>
"""
import math
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from m7b_common import *  # noqa: E402,F401,F403

d = sys.argv[1]
RUNS = f"{OUT}/runs"
TRUTH = {"c4_loo_na19650_10": "2+2", "cyp_loo_na19650_5": "2+3", "c4_lzo_hg00171_5": "2+2", "cyp_lzo_hg02769_5": "2+2"}
LN10 = math.log(10)
mapping = table(f"{RUNS}/pools/mapping.tsv")
td = {(r["cell"], r["canonical_diplotype"]): r for r in table(f"{RUNS}/truthdist/truthdist.tsv")} if os.path.exists(f"{RUNS}/truthdist/truthdist.tsv") else {}
rank_rows, summary, seeds, membership = [], [], [], []


def loc_opts(unit):
    p = f"{RUNS}/{unit}/options.tsv"
    return {r["genotype"]: r for r in table(p)} if os.path.exists(p) and os.path.exists(f"{RUNS}/{unit}/checks.tsv") and all(
        r["ok"] == "True" for r in table(f"{RUNS}/{unit}/checks.tsv")) else None


def ranking(scores):
    order = sorted(scores, key=lambda k: -scores[k])
    return order, {k: i + 1 for i, k in enumerate(order)}


for cell, c in CELLS.items():
    rows = [m for m in mapping if m["cell"] == cell and m["tier"] == "D"]
    dips = OrderedDict()
    for m in rows:
        dips.setdefault(m["canonical_diplotype"], []).append(m)
    pv = {}
    fs = f"{RUNS}/panvar.D.{cell}/fs.scores_truth_free.tsv"
    if os.path.exists(fs):
        for r in table(fs):
            pv[r["label"]] = r
    lik, tot = {}, {}
    for dip, ms in dips.items():
        vals = [pv[f"row{m['row']}"] for m in ms if f"row{m['row']}" in pv]
        if vals:
            lik[dip] = float(vals[0]["fragment_term"]) + float(vals[0]["exposure_term"])
            tot[dip] = max(float(v["total"]) for v in vals)
    loc = loc_opts(f"loc.{cell}.D.s1")
    lmean = {g: float(r["lik_mean_log10"]) for g, r in loc.items()} if loc else {}
    methods = {"panvar_F+E": lik, "panvar_total": tot, "locityper_s1": lmean}
    cls = {dip: ms[0]["unit_class"] for dip, ms in dips.items()}
    totc = {dip: int(ms[0]["unit_total"]) if ms[0]["unit_total"] not in ("", "NA") else None for dip, ms in dips.items()}
    caller = next(dip for dip, ms in dips.items() if any("reported" in m["labels"].split("|") for m in ms))
    truth = TRUTH[cell]
    tt = sum(int(x) for x in truth.split("+"))
    for meth, sc in methods.items():
        if not sc:
            summary.append(dict(cell=cell, method=meth, status="NOT_AVAILABLE"))
            continue
        order, rk = ranking(sc)
        unit = "nats" if meth.startswith("panvar") else "log10"
        for dip in order:
            ex = td.get((cell, dip), {})
            rank_rows.append(dict(cell=cell, method=meth, rank=rk[dip], canonical_diplotype=dip, labels="||".join(m["labels"] for m in dips[dip]),
                                  unit_class=cls[dip], unit_total=totc[dip], score=sc[dip], score_unit=unit, gap_to_top=sc[order[0]] - sc[dip],
                                  block_excess=ex.get("block_excess", ""), is_caller=dip == caller, truth_class=cls[dip] == truth))
        top = order[0]
        right_total = [k for k in order if totc[k] == tt]
        wrong_total = [k for k in order if totc[k] is not None and totc[k] != tt]
        right_cls = [k for k in right_total if cls[k] == truth]
        wrong_alloc = [k for k in right_total if cls[k] != truth]
        g = lambda a, b: (sc[a] - sc[b]) if a and b else ""
        summary.append(dict(
            cell=cell, method=meth, status="OK", score_unit=unit, top=top, top_class=cls[top], top_total=totc[top],
            top_block_excess=td.get((cell, top), {}).get("block_excess", ""),
            total_copy_number=f"best true-total rank {rk[right_total[0]] if right_total else 'none'} vs best wrong-total rank {rk[wrong_total[0]] if wrong_total else 'none'}; "
                              f"gap {g(right_total[0] if right_total else None, wrong_total[0] if wrong_total else None)}",
            allocation=f"best true-class rank {rk[right_cls[0]] if right_cls else 'none'} vs best wrong-allocation (true total) rank "
                       f"{rk[wrong_alloc[0]] if wrong_alloc else 'none'}; gap {g(right_cls[0] if right_cls else None, wrong_alloc[0] if wrong_alloc else None)}",
            excess_within_true_class=";".join(f"r{rk[k]}:{td.get((cell, k), {}).get('block_excess', '?')}" for k in right_cls),
            caller_rank=rk.get(caller), caller_gap_to_top=sc[top] - sc[caller] if caller in sc else "",
            control_deterioration=("YES" if rk.get(caller) != 1 else "no") if c["control"] else "n/a"))
    # seed stability
    for a, b in (("s1", "s1r"), ("s1", "s2")):
        A, B = loc_opts(f"loc.{cell}.D.{a}"), loc_opts(f"loc.{cell}.D.{b}")
        if not (A and B):
            continue
        sa = {g: float(r["lik_mean_log10"]) for g, r in A.items()}
        sb = {g: float(r["lik_mean_log10"]) for g, r in B.items()}
        oa, ra = ranking(sa)
        ob, rb = ranking(sb)
        diffs = {g: sb[g] - sa[g] for g in sa if g in sb}
        margin = sa[oa[0]] - sa[oa[1]] if len(oa) > 1 else float("inf")
        noise = max(abs(diffs[oa[0]]), abs(diffs[oa[1]])) if len(oa) > 1 else 0
        seeds.append(dict(cell=cell, comparison=f"{a}_vs_{b}", top_a=oa[0], top_b=ob[0], same_top=oa[0] == ob[0], rank_changes=sum(1 for g in ra if ra[g] != rb.get(g)),
                          max_abs_diff_log10=max(abs(v) for v in diffs.values()), mean_abs_diff_log10=sum(abs(v) for v in diffs.values()) / len(diffs),
                          winning_margin_a_log10=margin, top_two_seed_diff_log10=noise, verdict="INCONCLUSIVE" if margin <= noise else "STABLE"))
    # pool membership
    for pool in ("DX", "DO"):
        P_ = loc_opts(f"loc.{cell}.{pool}.s1")
        if not (P_ and loc):
            continue
        diffs = {g: float(P_[g]["lik_mean_log10"]) - float(loc[g]["lik_mean_log10"]) for g in loc if g in P_}
        allp = {g: float(r["lik_mean_log10"]) for g, r in P_.items()}
        o, rk = ranking(allp)
        membership.append(dict(cell=cell, pool=pool, d_genotypes=len(diffs), max_abs_shift_log10=max(abs(v) for v in diffs.values()),
                               d_top_in_D_run=max(lmean, key=lmean.get), d_rank_order_preserved=sorted(diffs, key=lambda g: -float(loc[g]["lik_mean_log10"])) ==
                               sorted(diffs, key=lambda g: -float(P_[g]["lik_mean_log10"])), top_in_pool_run=o[0],
                               top_in_pool_run_tier="D" if o[0] in loc else "non-D (diagnostic only)"))
write_tsv(f"{d}/rankings.tsv", rank_rows)
write_tsv(f"{d}/summary.tsv", summary)
if seeds:
    write_tsv(f"{d}/seeds.tsv", seeds)
if membership:
    write_tsv(f"{d}/membership.tsv", membership)
for r in summary:
    print("\t".join(f"{k}={v}" for k, v in r.items()))
for r in seeds + membership:
    print(r)
