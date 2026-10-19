#!/usr/bin/env python3
"""M6-F evaluation (frozen with DESIGN.md). For each evaluation family F and stratum: fit on the training families only
(m6f_core.fit_fold: inner leave-one-family-out selection of each model's grid point), then score every copy of F at the
selected parameters on one shared observation. Writes per-copy, per-family-stratum, per-family, per-stratum tables and
the decision under the frozen rule:

  per comparison (C1-C0, C2-C0, C2-C1), pooled per family over strata:
    meaningful   mean paired difference per evaluated copy >= 2.0 nats
    consistent   positive in >= 75 % of evaluated families AND the largest family's share of the summed positive
                 differences <= 0.5
  decision: C2_IMPROVES_CONSISTENTLY if C2-C1 is meaningful and consistent; else C1_SUFFICES if C1-C0 is meaningful and
            consistent; else PARK. (Strata get the same flags, descriptively.)

  evaluate.py <prepare dir> <events tsv> <out prefix>
"""
import builtins
import csv
import os
import sys
import time
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
prep, events_path, prefix = sys.argv[1:4]
import m6f_core as core  # noqa: E402

ALLOW = {os.path.abspath(p) for p in [core.M1D_FA, events_path, f"{prep}/copies.tsv", f"{prep}/prepare.tsv", f"{prep}/families.tsv"]}
_open = builtins.open
OUTDIR = os.path.dirname(os.path.abspath(prefix))


def guarded(path, mode="r", *a, **k):
    p = os.path.abspath(path)
    if "w" in mode or "a" in mode:
        assert p.startswith(OUTDIR), f"write outside the unit directory: {p}"
    else:
        assert p in ALLOW or p.startswith(OUTDIR), f"not an allowlisted input: {p}"
    return _open(path, mode, *a, **k)


builtins.open = guarded
copies, refs, comps, families, evaluation = core.load(prep, events_path)
strata = sorted(refs)
ref_len = {s: len(next(iter(comps[s].values())).base) for s in strata}
out, fits, timing = [], [], []
for fam in evaluation:
    for s in strata:
        mine = [c for c in copies if c["stratum"] == s and c["family"] == fam]
        if not mine:
            continue
        t0 = time.time()
        fit = core.fit_fold(copies, comps[s], s, fam)
        t_fit = time.time() - t0
        if fit is None:
            for c in mine:
                out.append(dict(family=fam, donor=c["donor"], copy=c["copy"], stratum=s, status="UNASSESSABLE_NO_TRAINING_LIBRARY"))
            continue
        sel = fit["selected"]
        fits.append(dict(family=fam, stratum=s, library_components=len(fit["library"]), library_weight=sum(fit["library"].values()),
                         mask_positions=int(fit["mask"].sum()), reference_bp=ref_len[s], inner_targets=fit["inner_targets"],
                         C0_eps=sel["C0"][0], C1_eps=sel["C1"][0], C2_eps=sel["C2"][0], C2_rho=sel["C2"][1], fit_seconds=round(t_fit, 3)))
        for c in mine:
            obs = core.observation(fit["library"], comps[s], fit["mask"], c["unique_record"])
            if obs.n_scored == 0:
                out.append(dict(family=fam, donor=c["donor"], copy=c["copy"], stratum=s, status="UNASSESSABLE_NO_SCORED_POSITION"))
                continue
            sc = core.score_selected(obs, sel)
            n = obs.n_scored
            out.append(dict(family=fam, donor=c["donor"], copy=c["copy"], stratum=s, status="SCORED", unique_record=c["unique_record"],
                            reference_bp=ref_len[s], excluded_training_mask=obs.excluded_training_mask, excluded_heldout_alignment=obs.excluded_heldout_alignment,
                            scored=n, coverage=round(n / ref_len[s], 4), training_variable_scored=len(obs.var_pos), novel_at_training_monomorphic=obs.novel_at_monomorphic,
                            logC0=sc["C0"], logC1=sc["C1"], logC2=sc["C2"], d_C1_C0=sc["C1"] - sc["C0"], d_C2_C0=sc["C2"] - sc["C0"], d_C2_C1=sc["C2"] - sc["C1"],
                            per_base_C0=sc["C0"] / n, per_base_C1=sc["C1"] / n, per_base_C2=sc["C2"] / n, observation_digest=sc["digest"][:16]))
        timing.append((fam, s, time.time() - t0))


def write(name, rows):
    if not rows:
        return
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(f"{prefix}.{name}.tsv", "w") as f:
        w = csv.DictWriter(f, fieldnames=keys, delimiter="\t", lineterminator="\n", restval="")
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.6f}" if isinstance(v, float) else v) for k, v in r.items()})


write("copies", out)
write("fits", fits)
scored = [r for r in out if r["status"] == "SCORED"]
DIFFS = ["d_C1_C0", "d_C2_C0", "d_C2_C1"]


def flags(per_family, n_copies):
    res = {}
    for d in DIFFS:
        vals = [v[d] for v in per_family.values()]
        tot = sum(vals)
        pos = [v for v in vals if v > 0]
        share = max(pos) / sum(pos) if pos else float("nan")
        frac = len(pos) / len(vals) if vals else float("nan")
        mean = tot / n_copies if n_copies else float("nan")
        res[d] = dict(sum=tot, mean_per_copy=mean, families_positive=f"{len(pos)}/{len(vals)}", largest_share=share,
                      meaningful=mean >= 2.0, consistent=bool(vals) and frac >= 0.75 and share <= 0.5)
    return res


fam_rows, strat_rows = [], []
for scope in strata + ["ALL"]:
    rows = [r for r in scored if scope == "ALL" or r["stratum"] == scope]
    if not rows:
        continue
    per_family = defaultdict(lambda: defaultdict(float))
    for r in rows:
        pf = per_family[r["family"]]
        for k in ["logC0", "logC1", "logC2"] + DIFFS + ["scored"]:
            pf[k] += r[k]
        pf["copies"] += 1
    for f_, v in sorted(per_family.items()):
        fam_rows.append(dict(scope=scope, family=f_, copies=int(v["copies"]), scored_positions=int(v["scored"]), logC0=v["logC0"], logC1=v["logC1"], logC2=v["logC2"],
                             d_C1_C0=v["d_C1_C0"], d_C2_C0=v["d_C2_C0"], d_C2_C1=v["d_C2_C1"],
                             per_base_d_C2_C1=v["d_C2_C1"] / v["scored"], per_base_d_C1_C0=v["d_C1_C0"] / v["scored"]))
    fl = flags(per_family, len(rows))
    for d, v in fl.items():
        strat_rows.append(dict(scope=scope, comparison=d, families=len(per_family), copies=len(rows), sum=v["sum"], mean_per_copy=v["mean_per_copy"],
                               mean_per_base=v["sum"] / sum(r["scored"] for r in rows), families_positive=v["families_positive"],
                               largest_family_share=v["largest_share"], meaningful=v["meaningful"], consistent=v["consistent"]))
write("families", fam_rows)
write("strata", strat_rows)
allf = {r["comparison"]: r for r in strat_rows if r["scope"] == "ALL"}
if allf and allf["d_C2_C1"]["meaningful"] and allf["d_C2_C1"]["consistent"]:
    decision = "C2_IMPROVES_CONSISTENTLY"
elif allf and allf["d_C1_C0"]["meaningful"] and allf["d_C1_C0"]["consistent"]:
    decision = "C1_SUFFICES"
else:
    decision = "PARK"
with open(f"{prefix}.decision.tsv", "w") as f:
    f.write(f"decision\t{decision}\n")
    f.write(f"copies_scored\t{len(scored)}\ncopies_unassessable\t{sum(1 for r in out if r['status'] != 'SCORED')}\n")
    f.write(f"seconds_total\t{sum(t for *_, t in timing):.2f}\n")
    for fam, s, t in timing:
        f.write(f"seconds\t{fam}\t{s}\t{t:.2f}\n")
print(open(f"{prefix}.decision.tsv").read())
for r in strat_rows:
    print("\t".join(f"{v:.4f}" if isinstance(v, float) else str(v) for v in r.values()))
