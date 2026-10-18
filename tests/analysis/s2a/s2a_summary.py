#!/usr/bin/env python3
"""S2-A real-data summary: per run, time to the first evaluated prediction (checkpoint) apart from the total runtime,
peak memory, the early evaluation's work beside the union pass's, and the final completion status; off/on completed
outputs compared (paths normalised; status timing keys, the declared early-incumbent rows and repair_memory excluded).
No accuracy claim.  s2a_summary.py <results dir>
"""
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s2a_fixtures as fx  # noqa: E402

KEYS = ["status", "sequence_prediction", "sequence_completion", "sequence_exit_code", "early_incumbent_status",
        "early_incumbent_union_check", "checkpoint_commits", "checkpoint_last_stage",
        "checkpoint_seconds_first_evaluated_prediction", "early_incumbent_seconds_preparation",
        "early_incumbent_seconds_selection", "early_incumbent_seconds_evaluation", "early_incumbent_pairs_examined",
        "early_incumbent_candidate_fragment_evaluations", "early_incumbent_physical_origins",
        "early_incumbent_footprint_after_bytes", "candidate_fragment_evaluations", "physical_origins",
        "seconds_total", "elapsed_seconds", "phys_footprint_peak_bytes", "checkpoint_retained", "phase", "note"]


def main(out):
    ledger = {r["unit"]: r for r in csv.DictReader(open(os.path.join(out, "ledger.tsv")), delimiter="\t")}
    rows = []
    for unit in ("c4.off", "c4.on", "cyp.off", "cyp.on", "acot30.on"):
        d = os.path.join(out, "runs", unit)
        st = fx.kv(os.path.join(d, "call.mosaic_status.tsv"))
        loc = fx.kv(os.path.join(d, "call.sequence_locus.tsv"))
        row = {"unit": unit, "exit": ledger.get(unit, {}).get("exit", "NA"),
               "wall_seconds": ledger.get(unit, {}).get("wall_seconds", "NA"),
               "peak_rss_gib": ledger.get(unit, {}).get("peak_rss_gib", "NA"),
               "locus_prediction": loc.get("prediction", "NA"), "locus_completion": loc.get("completion", "NA"),
               "checkpoint_stage": loc.get("checkpoint_stage", "NA")}
        for k in KEYS:
            row[k] = st.get(k, "NA")[:160]
        rows.append(row)
    for cell in ("c4", "cyp"):
        a, b = os.path.join(out, "runs", cell + ".off"), os.path.join(out, "runs", cell + ".on")
        if os.path.exists(a) and os.path.exists(b):
            diffs = fx.compare(a, b)
            rows.append({"unit": f"{cell}: off vs on", "exit": "-", "note": "identical" if not diffs else " | ".join(diffs[:4])})
    keys = list(rows[0])
    with open(os.path.join(out, "summary.tsv"), "w") as f:
        w = csv.DictWriter(f, fieldnames=keys, delimiter="\t", lineterminator="\n", restval="-")
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        print("\t".join(f"{k}={r.get(k, '-')}" for k in keys if r.get(k, "-") not in ("-", "NA")))
    return 0


if __name__ == "__main__":
    sys.exit(main(os.path.abspath(sys.argv[1])))
