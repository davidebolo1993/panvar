#!/usr/bin/env python3
"""The ten non-tie no-calls of C4 production-v1, from the frozen package only (no genotyping, no truth).

    scripts/c4_non_tie_no_call_inventory.py --campaign <campaign dir> --out <table.tsv>

Seven union-scoring-cap refusals (the ledger's rescore_incomplete: the union pass reached its physical-origin
cap before any repair) and three repair round caps. For each: job, experiment, coverage, status, the cap
and how far past or inside it the run was, the repair's last-round improvement for a round cap, runtime and
footprint. The cheapest representative of each mechanism is the one with the smallest recorded wall time.
"""
import argparse
import csv


def key_values(path):
    out = {}
    with open(path) as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t", 1)
            if len(parts) == 2 and parts[0] not in out:
                out[parts[0]] = parts[1]
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--campaign", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    jobs = {r["job"]: r for r in csv.DictReader(open(f"{args.campaign}/jobs.tsv"), delimiter="\t")}
    ledger = [r for r in csv.DictReader(open(f"{args.campaign}/completeness.tsv"), delimiter="\t")
              if r["outcome"] == "no_call" and r["detail"] != "improving_move_tie"]
    rows = []
    for entry in sorted(ledger, key=lambda r: (r["detail"], float(r["wall_seconds"]))):
        prefix = jobs[entry["job"]]["out_prefix"]
        status = key_values(prefix + ".mosaic_status.tsv")
        mechanism = "union_scoring_cap" if entry["detail"] == "rescore_incomplete" else entry["detail"]
        row = {"job": entry["job"], "experiment": entry["campaign"], "subject": entry["subject"],
               "coverage": entry["coverage"], "mechanism": mechanism, "status": status.get("status", "NA"),
               "ledger_detail": entry["detail"]}
        if mechanism == "union_scoring_cap":
            cap = int(status["max_physical_origins"])
            seen = int(status.get("observed_physical_origins", "0"))
            row.update({"cap": "max_physical_origins", "cap_value": cap, "observed": seen, "margin": cap - seen,
                        "observed_piece_hits": status.get("observed_piece_hits", "NA"),
                        "repair_rounds_run": "NA", "repair_moves": "NA", "last_round_improvement": "NA"})
        else:
            with open(prefix + ".repair_rounds.tsv") as handle:
                last = list(csv.DictReader(handle, delimiter="\t"))[-1]
            row.update({"cap": "repair_max_rounds", "cap_value": status["repair_max_rounds"],
                        "observed": status["repair_rounds_run"], "margin": int(status["repair_max_rounds"]) - int(status["repair_rounds_run"]),
                        "observed_piece_hits": "NA", "repair_rounds_run": status["repair_rounds_run"],
                        "repair_moves": status["repair_moves"], "last_round_improvement": last["delta_total"]})
        row.update({"repair_margin_physical_origins": status.get("repair_margin_physical_origins", "NA"),
                    "wall_seconds": entry["wall_seconds"], "peak_footprint_bytes": entry["peak_footprint_bytes"]})
        rows.append(row)
    cheapest = {}
    for row in rows:
        if row["mechanism"] not in cheapest or float(row["wall_seconds"]) < float(cheapest[row["mechanism"]]["wall_seconds"]):
            cheapest[row["mechanism"]] = row
    for row in rows:
        row["cheapest_representative"] = "true" if cheapest[row["mechanism"]] is row else "false"
    with open(args.out, "w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    for row in rows:
        print("\t".join(str(row[k]) for k in ("job", "mechanism", "coverage", "status", "cap_value", "observed", "margin",
                                              "last_round_improvement", "wall_seconds", "cheapest_representative")))


if __name__ == "__main__":
    main()
