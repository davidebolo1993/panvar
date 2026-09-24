#!/usr/bin/env python3
"""DIAGNOSTIC (v2/post-hoc): one row per synthetic-mosaic no-call cell of C4 production-v1.

    scripts/c4_mosaic_tie_summary.py --campaign <campaign dir> --out <diagnostic dir>

Reads the frozen campaign (completeness.tsv, each job's production status) and the diagnostic outputs of
tests/genotype_c4_mosaic_tie_diagnostic.sh, and writes <out>/mosaic_no_call_diagnostic.tsv: the original
refusal, the tied neighbours and their distinct unordered genotypes, branches, the pooled top set, the
diagnostic outcome, distance/floor/excess where the outcome can be evaluated, and the added cost. Every
diagnostic outcome is conditional on a pruned search and is not a call.
"""
import argparse
import csv
import pathlib


def kv(path):
    out = {}
    try:
        for line in open(path):
            parts = line.rstrip("\n").split("\t", 1)
            if len(parts) == 2 and parts[0] not in out:
                out[parts[0]] = parts[1]
    except FileNotFoundError:
        pass
    return out


def time_of(path):
    wall = footprint = "NA"
    try:
        for line in open(path):
            parts = line.split()
            if len(parts) >= 2 and parts[1] == "real":
                wall = parts[0]
            if line.strip().endswith("peak memory footprint"):
                footprint = parts[0]
    except FileNotFoundError:
        pass
    return wall, footprint


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--campaign", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    campaign, out = pathlib.Path(args.campaign), pathlib.Path(args.out)
    ledger = list(csv.DictReader(open(campaign / "completeness.tsv"), delimiter="\t"))
    rows = []
    for entry in sorted(ledger, key=lambda r: r["job"]):
        if not (entry["campaign"] == "mosaic" and entry["outcome"] == "no_call"):
            continue
        d = out / entry["job"]
        diag = kv(d / "call.diag.status.tsv")
        ties = list(csv.DictReader(open(d / "call.diag.terminal_ties.tsv"), delimiter="\t")) \
            if (d / "call.diag.terminal_ties.tsv").exists() else []
        first = [t for t in ties if t["tie"] == "0"]
        branches = list(csv.DictReader(open(d / "call.diag.branches.tsv"), delimiter="\t")) \
            if (d / "call.diag.branches.tsv").exists() else []
        evaluation = kv(d / "call.diag_eval.summary.tsv")
        if not evaluation and entry["detail"] == "union_scoring_cap":
            evaluation = kv(d / "call.eval.summary.tsv")   # the rerun's own call, when it made one
        rerun = kv(d / "call.mosaic_status.tsv")
        wall, footprint = time_of(d / "time")
        rows.append({
            "job": entry["job"], "plan": entry["subject"], "coverage": entry["coverage"],
            "original_refusal": entry["detail"],
            "tied_neighbours": len(first) if first else "NA",
            "tied_distinct_genotypes": len({t["unordered_genotype"] for t in first}) if first else "NA",
            "tied_repaired_blocks": ",".join(sorted({t["repaired_block"] for t in first})) if first else "NA",
            "nested_ties": max(len({t["tie"] for t in ties}) - 1, 0) if ties else "NA",
            "branches_run": diag.get("branches_run", "NA"),
            "branches_converged": sum(b["outcome"] == "converged" for b in branches) if branches else "NA",
            "branch_terminal_genotypes": diag.get("distinct_branch_terminal_genotypes", "NA"),
            "pooled_top_size": diag.get("pooled_top_size", "NA"),
            "pooled_top_genotypes": diag.get("pooled_top_genotypes", "NA"),
            "pooled_top_from": "NA",
            "diagnostic_outcome": (f"union_cap_rerun:{rerun.get('status', 'NA')}" if entry["detail"] == "union_scoring_cap" else "")
                                  + ((";" if entry["detail"] == "union_scoring_cap" else "") + diag["decision"] if "decision" in diag else ""),
            "incomplete_reason": diag.get("incomplete_reason", ""),
            "distance": evaluation.get("unordered_distance_total", "NA"),
            "floor": evaluation.get("availability_floor_total", "NA"),
            "excess": evaluation.get("excess_over_floor_total", "NA"),
            "evaluated_as": evaluation.get("evaluated_diplotype", "NA"),
            "diagnostic_seconds": diag.get("seconds", "NA"),
            "run_wall_seconds": wall, "run_peak_footprint_bytes": footprint,
            "production_wall_seconds": entry.get("wall_seconds", "NA"),
            "union_cap_rerun_physical_origins": rerun.get("observed_physical_origins", rerun.get("physical_origins", "NA"))
            if entry["detail"] == "union_scoring_cap" else "",
            "union_cap_rerun_overrides": rerun.get("settings_overridden_after_profile", "") if entry["detail"] == "union_scoring_cap" else "",
        })
        top = d / "call.diag.pooled_top.tsv"
        if top.exists():
            rows[-1]["pooled_top_from"] = ",".join(sorted({r["scored_in"].split("_")[0] for r in csv.DictReader(open(top), delimiter="\t")}))
    with open(out / "mosaic_no_call_diagnostic.tsv", "w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    for row in rows:
        print("\t".join(str(row[k]) for k in ("job", "original_refusal", "tied_neighbours", "tied_distinct_genotypes",
                                                "branches_run", "pooled_top_size", "pooled_top_genotypes",
                                                "diagnostic_outcome", "distance", "floor", "excess")))


if __name__ == "__main__":
    main()
