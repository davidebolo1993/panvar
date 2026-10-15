#!/usr/bin/env python3
"""Summarise bounded cache measurements: stage costs, work counters, completion and scientific identity.

Usage: perf_summary.py <out dir>   (reads <out dir>/runs/<cell>.<mode>/call.*; writes stage_costs.tsv, identity.tsv)
Stages (seconds, from the caller's own status): preparation = model + markers; proposal = proposal + spelling;
placement = union fragment preparation + placement and mass (union cache fills are inside it and reported separately as
"of which"); combination = diploid combination;
repair = adaptive repair (its scoring included); ties = tie branching; dosage = the dosage stage total; other = total -
the stages named (library sensitivity, annotation, output). An incomplete run reports what it completed and is never
turned into a speed-up. Complete: exit 0, or exit 1 with a scientific no-call; a run stopped by a declared resource cap
(status note names the cap) or by a limit is incomplete.
"""
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import perf_compare  # noqa: E402

OUT = sys.argv[1]


def kv(path):
    out = {}
    if os.path.exists(path):
        for line in open(path).read().splitlines():
            p = line.split("\t")
            if len(p) >= 2:
                out.setdefault(p[0], p[1])
    return out


def f(st, *keys):
    total = 0.0
    for k in keys:
        try:
            total += float(st.get(k, 0) or 0)
        except ValueError:
            pass
    return total


rows = []
runs = os.path.join(OUT, "runs")
for d in sorted(os.listdir(runs)):
    cell, mode = d.rsplit(".", 1)
    p = os.path.join(runs, d, "call")
    if not os.path.exists(p + ".run.tsv"):
        continue
    rec = dict(zip(*[l.split("\t") for l in open(p + ".run.tsv").read().splitlines()[:2]]))
    st, ds = kv(p + ".mosaic_status.tsv"), kv(p + ".dosage_status.tsv")
    total = f(st, "seconds_total")
    stages = {"preparation": f(st, "seconds_model", "seconds_markers"),
              "proposal": f(st, "seconds_proposal", "seconds_spelling"),
              "placement": f(st, "seconds_rescoring_fragment_preparation", "seconds_rescoring_placement_and_mass"),
              "combination": f(st, "seconds_rescoring_diploid_combination"),
              "repair": f(st, "seconds_repair"), "ties": f(st, "seconds_tie_branching"),
              "dosage": f(ds, "seconds_dosage_total")}
    named = sum(v for k, v in stages.items() if k != "dosage")
    capped = "resource cap" in st.get("note", "") or rec["exit"] not in ("0", "1")
    complete = "no" if capped else "yes"
    rows.append({"cell": cell, "mode": mode, "exit": rec["exit"], "complete": complete,
                 "incomplete_reason": st.get("note", "")[:110] if capped else "none", "status": st.get("status", "NA"),
                 "call_class": st.get("call_class", "NA"), "wall_s": rec["wall_seconds"], "caller_total_s": f"{total:.1f}",
                 **{f"{k}_s": f"{v:.1f}" for k, v in stages.items()},
                 "of_which_cache_fill_s": f"{f(st, 'seconds_union_placement_cache_fill', 'seconds_union_origin_cache_fill'):.1f}",
                 "other_s": f"{max(0.0, total - named):.1f}",
                 "peak_rss_GiB": f"{int(rec['peak_rss_bytes']) / 2 ** 30:.2f}" if rec["peak_rss_bytes"].isdigit() else "NA",
                 "proposed_starts": st.get("proposed_starts", "NA"), "full_read_verifications": st.get("full_read_verifications", "NA"),
                 "dp_evaluations": st.get("dp_evaluations", "NA"), "physical_origins": st.get("physical_origins", "NA"),
                 "statistical_state_classes": st.get("statistical_state_classes", "NA"),
                 "union_placement_cache_hits": st.get("union_placement_cache_hits", "NA"),
                 "union_origin_cache_hits": st.get("union_origin_cache_hits", "NA"),
                 "emission_cache_dp_hits": st.get("emission_cache_dp_hits", "NA"),
                 "repair_rounds_run": st.get("repair_rounds_run", "NA")})
with open(os.path.join(OUT, "stage_costs.tsv"), "w") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
ident = []
for cell in sorted({r["cell"] for r in rows}):
    a, b = os.path.join(runs, f"{cell}.default"), os.path.join(runs, f"{cell}.c4v6")
    ra = next((r for r in rows if r["cell"] == cell and r["mode"] == "default"), None)
    rb = next((r for r in rows if r["cell"] == cell and r["mode"] == "c4v6"), None)
    if not (ra and rb):
        ident.append({"cell": cell, "both_complete": "no", "identical": "NA", "detail": "one arm missing", "speedup_total": "NA"})
        continue
    if ra["complete"] != "yes" or rb["complete"] != "yes":
        ident.append({"cell": cell, "both_complete": "no", "identical": "NA",
                      "detail": f"default {ra['exit']}/{ra['status']}, c4v6 {rb['exit']}/{rb['status']}", "speedup_total": "NA"})
        continue
    ok, problems, diag, n = perf_compare.compare(a, b)
    ident.append({"cell": cell, "both_complete": "yes", "identical": "yes" if ok else "NO",
                  "detail": f"{n} files compared; {len(diag)} diagnostic keys changed" + ("" if ok else "; " + " | ".join(problems[:5])),
                  "speedup_total": f"{float(ra['caller_total_s']) / max(1e-9, float(rb['caller_total_s'])):.2f}"})
with open(os.path.join(OUT, "identity.tsv"), "w") as fh:
    w = csv.DictWriter(fh, fieldnames=list(ident[0]), delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(ident)
for r in ident:
    print(r)
