#!/usr/bin/env python3
"""Reconcile the generic report's data/ tables against the sealed C4 tables they must reproduce:
  1 the 484-job campaign package (call classes)
  2 the 649-cell final_26401b4 package (call classes, per-block evaluator values, layer A, mosaic strata)
  3 the headline_filters_v2 audit (locus headline values; policies.tsv counts, every stratum)
  4 the phase evaluation (summaries by experiment, switch class x coverage, truth switch category, per-call status,
    LOO projection)
Usage: reconcile.py --report <report dir> [--results-root DIR] [--out FILE]. Exit 1 on any mismatch."""
import os
import sys

import numpy as np
import pandas as pd

import c4_common as C

p = C.parser(__doc__.split("\n")[0], work=False)
p.add_argument("--report", required=True, help="the report directory (its data/ is read)")
p.add_argument("--out", help="reconciliation table (default: <report>/../reconciliation.tsv)")
args = p.parse_args()
_, RES = C.roots(args)
REP = os.path.abspath(args.report)
D = os.path.join(REP, "data")
F = os.path.join(RES, C.FINAL, "package")
H = os.path.join(RES, C.HEADLINE)
PH = os.path.join(RES, C.PHASE, "tables")
CAMPAIGN_RUNS = os.path.join(RES, C.CAMPAIGN, "package", "runs.tsv")
LABEL = {"lzo": "LZO", "loo": "LOO", "mosaic": "mosaic", "fresh_lzo": "fresh_LZO", "fresh_loo": "fresh_LOO"}
rd = lambda p: pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False)
rows = []


def check(section, what, expected, got):
    ok = str(expected) == str(got)
    rows.append({"section": section, "check": what, "expected": expected, "report": got, "match": ok})
    if not ok:
        print(f"MISMATCH {section}: {what}: expected {expected}, report {got}")


runs = rd(f"{D}/runs.tsv")
cls = lambda frame: dict(frame.call_class.value_counts().sort_index())

# ---- 1 the 484-job campaign
camp = rd(CAMPAIGN_RUNS)
check("1 campaign 484", "jobs", len(camp), int(runs.run_id.str.startswith("campaign/").sum()))
check("1 campaign 484", "call classes", cls(camp), cls(runs[runs.run_id.str.startswith("campaign/")]))
cmap = dict(zip("campaign/" + camp.job, camp.call_class))
check("1 campaign 484", "per-job call class", "identical",
      "identical" if all(cmap[r] == c for r, c in zip(runs.run_id, runs.call_class) if r in cmap) else "differs")
call = rd(f"{D}/callability.tsv")
for e in ("lzo", "loo", "mosaic"):
    sub = camp[camp.experiment == e]
    for cc, n in sorted(sub.call_class.value_counts().items()):
        got = pd.to_numeric(call[(call.experiment == e) & (call.call_class == cc)].runs).sum()
        check("1 campaign 484", f"{e} {cc} runs (callability.tsv)", n, got)

# ---- 2 the 649-cell package
pr = rd(f"{F}/runs.tsv")
pr["run_id"] = pr.dataset + "/" + pr.job
check("2 replay 649", "runs", len(pr), len(runs))
check("2 replay 649", "call classes", cls(pr), cls(runs))
nc = call[call.call_class == "no_call"]
check("2 replay 649", "no-call reasons", "tie_branching_incomplete:2",
      ";".join(f"{r}:{pd.to_numeric(g.runs).sum()}" for r, g in nc.groupby("no_call_reason")))
pb = rd(f"{F}/blocks.tsv")
pb["run_id"] = pb.dataset + "/" + pb.job
ev = rd(f"{D}/block_evaluation.tsv")
m = pb.merge(ev, on=["run_id", "block_id"], suffixes=("_pkg", "_rep"))
check("2 replay 649", "reported blocks", len(pb), len(ev))
check("2 replay 649", "blocks joined", len(pb), len(m))
num = lambda s: pd.to_numeric(s.replace({"NA": None, "": None, "nan": None}), errors="coerce")
for col in ("distance", "floor", "excess", "normalized_excess"):
    a, b = num(m[col + "_pkg"]), num(m[col + "_rep"])
    same = ((a - b).abs() <= 1e-12 * np.maximum(1, a.abs())) | (a.isna() & b.isna())
    check("2 replay 649", f"per-block {col} (7117 rows)", "all equal", "all equal" if same.all() else f"{int((~same).sum())} differ")
check("2 replay 649", "per-block truth_defined", "all equal",
      "all equal" if (m.truth_defined_pkg == m.truth_defined_rep).all() else "differ")
la = rd(f"{F}/layer_A_primary.tsv")
defined = ev[ev.truth_defined == "True"]
a0 = la[la.experiment == "all"].iloc[0]
for k, v in (("blocks", len(ev)), ("truth_defined", len(defined)), ("truth_undefined", len(ev) - len(defined)),
             ("exact", int((num(defined.distance) == 0).sum())), ("floor_attained", int((num(defined.excess) == 0).sum())),
             ("excess_total", float(num(defined.excess).sum())), ("hard_no_calls", int((runs.call_class == "no_call").sum()))):
    check("2 replay 649", f"layer A all: {k}", a0[k] if k != "excess_total" else float(a0[k]), v)
for r in la[la.experiment != "all"].itertuples():
    g = defined[(defined.experiment == r.experiment) & (defined.coverage == r.coverage)]
    check("2 replay 649", f"layer A {r.experiment} {r.coverage}x floor attained / defined",
          f"{r.floor_attained}/{r.truth_defined}", f"{int((num(g.excess) == 0).sum())}/{len(g)}")
ms = rd(f"{F}/mosaic_strata.tsv")
ms = ms[(ms.stratum == "switch_class") & (ms.quality == "all_quality")]
mx = rd(f"{D}/mosaic_excess_summary.tsv")
for r in ms.itertuples():
    g = mx[(mx.switch_class == r.value) & (mx.coverage == r.coverage)]
    check("2 replay 649", f"mosaic {r.value} {r.coverage}x floor attained / blocks", f"{r.floor_attained}/{r.blocks}",
          f"{g.floor_attained.iloc[0]}/{g.denominator.iloc[0]}" if len(g) else "absent")

# ---- 3 the headline audit
dl = rd(f"{H}/derived_v2/loci.tsv")
hl = rd(f"{D}/headline_locus_status.tsv")
for f in ("GENOTYPE_FILTER", "DOSAGE_FILTER", "PHASE_FILTER"):
    exp = {k: int(v) for k, v in dl[f].value_counts().sort_index().items()}
    got = {s: int(pd.to_numeric(g.runs).sum()) for s, g in hl[hl["filter"] == f].groupby("status")}
    check("3 headline audit", f"locus {f} counts", exp, dict(sorted(got.items())))
pol = rd(f"{H}/eval_v2/policies.tsv")
fp = rd(f"{D}/filter_performance.tsv")
name = {"GENOTYPE_FILTER==PASS (reference)": "GENOTYPE_FILTER == PASS", "EXACT_GENOTYPE_USABLE": "EXACT_GENOTYPE_USABLE",
        "DOSAGE_CONTENT_USABLE": "DOSAGE_CONTENT_USABLE"}
cols = [("correct_retained", "retained_correct"), ("incorrect_retained", "retained_incorrect"),
        ("correct_excluded", "filtered_correct"), ("incorrect_excluded", "filtered_incorrect")]
n_pol = 0
for r in pol[pol.policy.isin(name)].itertuples():
    sel = (fp.policy == name[r.policy]) & (fp.stratum == r.stratum)
    if r.experiment:
        sel &= fp.experiment == {v: k for k, v in LABEL.items()}[r.experiment]
    if r.coverage:
        sel &= fp.coverage == r.coverage
    if r.block_id:
        sel &= fp.block_id == r.block_id
    g = fp[sel]
    what = f"{r.policy} | {r.stratum} {r.experiment} {r.coverage} {r.block_id}".strip()
    exp = [r.blocks] + [getattr(r, a) for a, _ in cols]
    got = ([g.units.iloc[0]] + [g[b].iloc[0] for _, b in cols]) if len(g) == 1 else f"{len(g)} rows"
    if r.policy != "DOSAGE_CONTENT_USABLE":
        exp += [r.catastrophic_retained, r.catastrophic_excluded]
        got = got + [g.catastrophic_retained.iloc[0], g.catastrophic_filtered.iloc[0]] if len(g) == 1 else got
    check("3 headline audit", what, exp, got)
    n_pol += 1
print(f"headline policy rows reconciled: {n_pol}")

# ---- 4 phase
pa = rd(f"{D}/phase_accuracy.tsv")
for c in ("assessable", "correct", "switch_errors"):
    pa[c] = pd.to_numeric(pa[c])
se = rd(f"{PH}/summary_by_experiment.tsv")
se = se[se.phase_support == "ALL"]
for r in se.itertuples():
    e = {v: k for k, v in LABEL.items()}[r.experiment]
    g = pa[pa.experiment == e]
    check("4 phase", f"{r.experiment} assessable / correct / switch errors",
          f"{r.assessable_boundaries}/{r.correct_boundaries}/{r.switch_errors}",
          f"{g.assessable.sum()}/{g.correct.sum()}/{g.switch_errors.sum()}")
pls = rd(f"{D}/phase_locus_status.tsv")
for r in se.itertuples():
    e = {v: k for k, v in LABEL.items()}[r.experiment]
    g = pls[pls.experiment == e]
    check("4 phase", f"{r.experiment} strict / GE / exact phased locus",
          f"{r.strict_calls}/{r.genotype_equivalent_unresolved}/{r.exact_phased_locus_calls}",
          f"{pd.to_numeric(g.strict_calls).sum()}/{pd.to_numeric(g.genotype_equivalent_unphased).sum()}/"
          f"{pd.to_numeric(g.exact_phased_locus).sum()}")
sc = rd(f"{PH}/summary_mosaic_by_switch_class_coverage.tsv")
pc = rd(f"{D}/phase_by_switch_class.tsv")
for r in sc[sc.phase_support == "ALL"].itertuples():
    g = pc[(pc.switch_class == r.mosaic_class) & (pc.coverage == r.coverage)]
    check("4 phase", f"mosaic {r.mosaic_class} {r.coverage}x assessable / switch errors",
          f"{r.assessable_boundaries}/{r.switch_errors}", f"{g.assessable.iloc[0]}/{g.switch_errors.iloc[0]}" if len(g) else "absent")
tc = rd(f"{PH}/summary_by_experiment_truth_switch_category.tsv")
pt = rd(f"{D}/phase_by_boundary_type.tsv")
for r in tc[tc.phase_support == "ALL"].itertuples():
    e = {v: k for k, v in LABEL.items()}[r.experiment]
    cat = r.truth_switch_category if e == "mosaic" else "no switch plan"
    g = pt[(pt.experiment == e) & (pt.truth_switch_category == cat)]
    check("4 phase", f"{r.experiment} truth switch {r.truth_switch_category} assessable / correct",
          f"{r.assessable_boundaries}/{r.correct_boundaries}", f"{g.assessable.iloc[0]}/{g.correct.iloc[0]}" if len(g) else "absent")
lp = rd(f"{PH}/loo_sequence_projected_summary.tsv")
pp = rd(f"{D}/phase_projected.tsv")
for r in lp[(lp.coverage == "") & (lp.phase_support == "")].itertuples():
    e = {v: k for k, v in LABEL.items()}[r.experiment]
    g = pp[pp.experiment == e]
    check("4 phase", f"{r.experiment} projected evaluated / concordant", f"{r.evaluated}/{r.concordant}",
          f"{pd.to_numeric(g.evaluated).sum()}/{pd.to_numeric(g.concordant).sum()}")

out = pd.DataFrame(rows)
out.to_csv(args.out or os.path.join(os.path.dirname(REP), "reconciliation.tsv"), sep="\t", index=False, lineterminator="\n")
bad = int((~out.match).sum())
print(out.groupby("section").match.agg(["sum", "count"]).to_string())
print(f"{len(out) - bad} of {len(out)} checks match")
sys.exit(1 if bad else 0)
