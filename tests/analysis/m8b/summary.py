#!/usr/bin/env python3
"""M8-B summary (frozen with DESIGN.md).

Per locus: complete training objective profiles per mode; selection (max; ties within 1e-9 -> smaller r / larger lambda);
BOUNDARY if the selection is a grid edge; flatness (best - worst, best - legacy); family sensitivity (selection after
removing each training family's contribution); isolation (training outputs identical with scrambled eval rows); outer
evaluation of every eval family (all-masked included) at legacy r = 0.01, fitted constant, fitted distance: paired
per-family deltas, meaningful (mean >= 1 nat per haplotype) and consistent (>= 75 % families positive, largest positive
share <= 0.5); whole-sequence support counts.

  summary.py <runs dir> <out prefix>
"""
import csv
import filecmp
import os
import sys
from collections import defaultdict

runs, prefix = sys.argv[1:3]
T = lambda p: list(csv.DictReader(open(p), delimiter="\t"))
LEGACY = 0.01
out_rows, prof_rows, eval_rows, fam_rows, sup_rows = [], [], [], [], []


def select(profile, mode):
    pts = sorted(((float(r["param"]), float(r["training_objective"])) for r in profile if r["mode"] == mode), key=lambda x: x[0])
    best = max(v for _, v in pts)
    tied = [p for p, v in pts if v >= best - 1e-9]
    pick = min(tied) if mode == "constant" else max(tied)
    return pick, pts, pick in (pts[0][0], pts[-1][0])


for locus in ("C4", "CYP2D6", "ACOT"):
    d = f"{runs}/fit.{locus}"
    if not os.path.exists(f"{d}/h.profile.tsv"):
        continue
    prof = T(f"{d}/h.profile.tsv")
    famp = T(f"{d}/h.family_profile.tsv")
    info = {r["key"]: r["value"] for r in T(f"{d}/h.info.tsv")}
    iso = all(filecmp.cmp(f"{d}/h.{x}.tsv", f"{runs}/fit_scrambled.{locus}/h.{x}.tsv", shallow=False) for x in ("profile", "family_profile")) \
        if os.path.exists(f"{runs}/fit_scrambled.{locus}/h.profile.tsv") else None
    sel = {}
    for mode in ("constant", "distance"):
        pick, pts, edge = select(prof, mode)
        vals = [v for _, v in pts]
        legacy = next(float(r["training_objective"]) for r in prof if r["mode"] == "constant" and abs(float(r["param"]) - LEGACY) < 1e-12)
        changed = 0
        fams = sorted({r["family"] for r in famp})
        for f in fams:
            sub = defaultdict(float)
            for r in prof:
                if r["mode"] == mode:
                    sub[float(r["param"])] = float(r["training_objective"])
            for r in famp:
                if r["family"] == f and r["mode"] == mode:
                    sub[float(r["param"])] -= float(r["contribution"])
            best = max(sub.values())
            tied = [p for p, v in sub.items() if v >= best - 1e-9]
            changed += (min(tied) if mode == "constant" else max(tied)) != pick
        sel[mode] = pick
        for p, v in pts:
            prof_rows.append(dict(locus=locus, mode=mode, param=p, training_objective=v, minus_best=v - max(vals)))
        out_rows.append(dict(locus=locus, mode=mode, selected=pick, boundary="BOUNDARY" if edge else "interior", profile_range=max(vals) - min(vals),
                             best_minus_legacy=max(vals) - legacy, families_whose_removal_changes_selection=f"{changed}/{len(fams)}",
                             training_haplotypes=prof[0]["training_haplotypes"], training_masked_blocks=prof[0]["training_masked_blocks"],
                             minus_infinite=max(int(r["minus_infinite"]) for r in prof if r["mode"] == mode), isolation=("PASS" if iso else "FAIL") if iso is not None else "not run",
                             eval_families=info.get("eval_families"), lambda_hm_bp=info.get("lambda_hm_bp")))
    ev = T(f"{d}/h.eval.tsv")
    if ev:
        B = int(info["blocks"])
        per = defaultdict(lambda: defaultdict(float))
        hap = defaultdict(set)
        masked = defaultdict(int)
        for r in ev:
            p = float(r["param"])
            key = ("legacy" if r["mode"] == "constant" and abs(p - LEGACY) < 1e-12 else None,
                   "fitted_constant" if r["mode"] == "constant" and abs(p - sel["constant"]) < 1e-12 * max(1, p) else None,
                   "fitted_distance" if r["mode"] == "distance" and abs(p - sel["distance"]) <= 1e-9 * p else None)
            for k in key:
                if k:
                    per[r["family"]][k] += float(r["loglik"])
            hap[r["family"]].add(r["template"])
            if r["mode"] == "constant" and abs(p - LEGACY) < 1e-12:
                masked[r["family"]] += int(r["masked_blocks"])
        fams = sorted(per)
        for f in fams:
            fam_rows.append(dict(locus=locus, family=f, haplotypes=len(hap[f]), mask_fraction=masked[f] / (B * len(hap[f])),
                                 legacy=per[f]["legacy"], fitted_constant=per[f]["fitted_constant"], fitted_distance=per[f]["fitted_distance"],
                                 d_constant_vs_legacy=per[f]["fitted_constant"] - per[f]["legacy"], d_distance_vs_constant=per[f]["fitted_distance"] - per[f]["fitted_constant"]))
        nh = sum(len(hap[f]) for f in fams)
        for comp in ("d_constant_vs_legacy", "d_distance_vs_constant"):
            vals = [r[comp] for r in fam_rows if r["locus"] == locus]
            pos = [v for v in vals if v > 1e-9]
            share = max(pos) / sum(pos) if pos else float("nan")
            mean = sum(vals) / nh
            eval_rows.append(dict(locus=locus, comparison=comp, families=len(vals), haplotypes=nh, total=sum(vals), mean_per_haplotype=mean,
                                  families_positive=f"{len(pos)}/{len(vals)}", largest_positive_share=share, meaningful=mean >= 1.0,
                                  consistent=bool(vals) and len(pos) / len(vals) >= 0.75 and (share <= 0.5 if pos else False)))
    sp = f"{runs}/support.{locus}/h.support.tsv"
    if os.path.exists(sp):
        s = T(sp)
        for role in ("eval", "train"):
            rs = [r for r in s if r["role"] == role]
            if rs:
                sup_rows.append(dict(locus=locus, role=role, haplotypes=len(rs), tested=sum(1 for r in rs if not r["status"].startswith("UNTESTED")),
                                     whole_sequence_supported=sum(1 for r in rs if r["status"] == "WHOLE_SEQUENCE_SUPPORTED"),
                                     whole_sequence_unsupported=sum(1 for r in rs if r["status"] == "WHOLE_SEQUENCE_UNSUPPORTED"),
                                     untested=sum(1 for r in rs if r["status"].startswith("UNTESTED")),
                                     supported_despite_masked_blocks=sum(1 for r in rs if r["status"] == "WHOLE_SEQUENCE_SUPPORTED" and r["bank_masked_blocks"] not in ("", "0"))))


def write(name, rows):
    if rows:
        with open(f"{prefix}.{name}.tsv", "w") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(rows)


write("selection", out_rows)
write("profiles", prof_rows)
write("evaluation", eval_rows)
write("families", fam_rows)
write("support", sup_rows)
for name in ("selection", "evaluation", "support"):
    if os.path.exists(f"{prefix}.{name}.tsv"):
        print(open(f"{prefix}.{name}.tsv").read())
