#!/usr/bin/env python3
"""M1-B: historical vs current outcome for each audit case (saved outputs only).

Target-block verdict (declared before execution):
  PERSISTS    the fresh call carries the same unordered allele pair at the target block as the historical call
  DISAPPEARS  the fresh target-block excess is <= 100 bp
  CHANGES     otherwise (a different allele pair whose excess is still > 100 bp)
  NOT_ASSESSED the fresh caller or its evaluation did not complete
Also reported: model fingerprint equality (allele indices are comparable only when equal), status, prediction, phase support,
whole-locus distance/floor/excess, every other SV block whose called pair or excess changed, and the fresh effective settings
against the M1-A frozen baseline settings for the matching baseline cell (each difference listed).

  outcome.py AUDIT_SET.tsv RUNS_DIR BASELINE_SETTINGS.tsv OUT_PREFIX
"""
import csv
import os
import sys

BASELINE_CELL = {"c4_10": "c4_30", "c4_30": "c4_30", "cyp_5": "cyp_5", "cyp_30": "cyp_30"}


def kv(path):
    out = {}
    if os.path.exists(path):
        for line in open(path, errors="replace"):
            k, _, v = line.rstrip("\n").partition("\t")
            out[k] = v
    return out


def blocks(path):
    return {r["block_id"]: r for r in csv.DictReader(open(path), delimiter="\t")} if os.path.exists(path) else {}


def totals(bl):
    d = f = 0
    for r in bl.values():
        if r["unordered_distance"] not in ("NA", ""):
            d += int(float(r["unordered_distance"]))
            f += int(float(r["availability_floor"]))
    return d, f, d - f


def main(audit, runs, baseline, out):
    cases = list(csv.DictReader(open(audit), delimiter="\t"))
    base = {}
    for r in csv.DictReader(open(baseline), delimiter="\t"):
        base.setdefault(r["cell"], {})[r["setting"]] = r["value"]
    rows, diffs = [], []
    for c in cases:
        name, blk = c["case"], c["target_block"]
        hs, hb = kv(c["historical_run"] + ".mosaic_status.tsv"), blocks(c["historical_eval"] + ".blocks.tsv")
        fp = f"{runs}/call.{name}/call"
        fs, fl = kv(fp + ".mosaic_status.tsv"), kv(fp + ".sequence_locus.tsv")
        fb = blocks(f"{runs}/eval.{name}/eval.blocks.tsv")
        h, f = hb.get(blk), fb.get(blk)
        if not f or not fs:
            verdict = "NOT_ASSESSED"
        elif sorted([h["called_allele1"], h["called_allele2"]]) == sorted([f["called_allele1"], f["called_allele2"]]):
            verdict = "PERSISTS"
        elif float(f["excess"]) <= 100:
            verdict = "DISAPPEARS"
        else:
            verdict = "CHANGES"
        other = []
        for b, hr in hb.items():
            fr = fb.get(b)
            if b == blk or fr is None or hr["kind"] != "SV":
                continue
            if sorted([hr["called_allele1"], hr["called_allele2"]]) != sorted([fr["called_allele1"], fr["called_allele2"]]) \
                    or hr["excess"] != fr["excess"]:
                other.append(f"{b}:{hr['called_allele1']}/{hr['called_allele2']}(x{hr['excess']})->"
                             f"{fr['called_allele1']}/{fr['called_allele2']}(x{fr['excess']})")
        ht, ft = totals(hb), totals(fb)
        rows.append({
            "case": name, "target_block": blk,
            "fingerprint_equal": str(hs.get("model_fingerprint") == fs.get("model_fingerprint")) if fs else "NA",
            "historical_status": hs.get("status", "NA"), "fresh_status": fs.get("status", "NA"),
            "fresh_prediction": fl.get("prediction", "NA"), "fresh_completion": fl.get("completion", "NA"),
            "historical_phase": hs.get("phase_support", "NA"), "fresh_phase": fs.get("phase_support", "NA"),
            "historical_target": f"{h['called_allele1']}/{h['called_allele2']} {h['called1_bp']}+{h['called2_bp']} bp excess {h['excess']}" if h else "NA",
            "truth_target": f"{h['truth1_bp']}+{h['truth2_bp']} bp floor {h['availability_floor']}" if h else "NA",
            "fresh_target": f"{f['called_allele1']}/{f['called_allele2']} {f['called1_bp']}+{f['called2_bp']} bp excess {f['excess']} support {f['genotype_support']}" if f else "NA",
            "verdict": verdict,
            "historical_locus_distance_floor_excess": "/".join(map(str, ht)),
            "fresh_locus_distance_floor_excess": "/".join(map(str, ft)) if fb else "NA",
            "other_sv_block_changes": ";".join(other) or "none",
            "fresh_seconds_total": fs.get("seconds_total", "NA"), "fresh_repair_rounds": fs.get("repair_rounds_run", "NA"),
            "fresh_tie_branches": fs.get("tie_branches_run", "NA"),
        })
        st = {k: v.split("\t")[0] for k, v in kv(fp + ".settings.tsv").items()}
        b = base.get(BASELINE_CELL[name], {})
        for k in sorted(set(st) | set(b)):
            if k in ("setting", "key"):
                continue
            if st.get(k) != b.get(k):
                diffs.append((name, k, b.get(k, "ABSENT"), st.get(k, "ABSENT")))
    keys = list(rows[0].keys())
    with open(out + ".outcome.tsv", "w") as fh:
        fh.write("\t".join(keys) + "\n")
        for r in rows:
            fh.write("\t".join(r[k] for k in keys) + "\n")
    with open(out + ".settings_vs_baseline.tsv", "w") as fh:
        fh.write("case\tsetting\tbaseline\tfresh\n")
        for d in diffs:
            fh.write("\t".join(d) + "\n")
    for r in rows:
        print(f"{r['case']}: {r['verdict']} | fresh {r['fresh_status']} | {r['fresh_target']} | hist {r['historical_target']}")
    print(f"settings differences vs baseline: {len(diffs)} (listed in {out}.settings_vs_baseline.tsv)")


if __name__ == "__main__":
    main(*sys.argv[1:5])
