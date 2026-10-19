#!/usr/bin/env python3
"""M3-F accuracy per arm (TRUTH; runs only after the predictions are frozen).

Per cell and arm, from the evaluator (eval.<cell>.<arm>/eval.summary.tsv, eval.blocks.tsv):
  * sequence: excess_over_floor_total, availability floor, unordered distance;
  * structure, kept separate: blocks_exact, sv_blocks_exact, SV vs context unordered distance;
  * phase, kept separate: phased distance and phase-switch blocks (NA for an unphased call), the call's phase_support;
  * genotype-data support (support v2), separate from phase: per block, genotype_support_v2x / quality_filter_v2x against
    whether that block is exact (per-block excess 0).
Outcome of on vs off per cell (frozen):
  REFUSAL_OR_NO_CALL  either arm wrote no evaluable call (kept in every denominator; named by arm)
  RECOVERED           on excess 0 < off excess
  IMPROVED            0 < on excess < off excess, and no structural loss
  HARMFUL             on excess > off excess, or fewer exact blocks, or fewer exact SV blocks
  MIXED               lower excess but fewer exact blocks or SV blocks (reported as harmful in structure)
  CHANGED_EQUAL       same excess and exactness, different selected pair or block projection
  UNCHANGED           same excess and exactness and the same selection
STATUS_DOWNGRADE is flagged separately (COMPLETE -> genotype-equivalent / incomplete / no-call).
Control stop rule (target 'control'): exit 1 if any control is HARMFUL, MIXED, REFUSAL_OR_NO_CALL in the on arm only, or
STATUS_DOWNGRADE.

  accuracy_m3f.py <runs dir> control|all <out prefix> cell...
"""
import csv
import os
import sys
from collections import Counter

runs, target, prefix, cells = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:]
KEYS = ["status", "run_no_call", "excess_over_floor_total", "availability_floor_total", "unordered_distance_total", "blocks", "blocks_exact",
        "sv_blocks", "sv_blocks_exact", "sv_unordered_distance_total", "context_unordered_distance_total", "phased_distance_total",
        "phase_switch_blocks", "phase_support", "pass_blocks", "pass_blocks_exact", "floor_attained_in_candidates",
        "closest_candidate_total_score_gap_to_winner", "closest_candidate_fragment_gap_to_winner", "truth_candidate_source"]
RANK = {"COMPLETE": 0, "CONDITIONAL_GENOTYPE_EQUIVALENT": 1}


def kv(path):
    if not os.path.exists(path):
        return {}
    out = {}
    for line in open(path):
        k, _, v = line.rstrip("\n").partition("\t")
        out.setdefault(k, v)
    return out


def table(path):
    return list(csv.DictReader(open(path), delimiter="\t")) if os.path.exists(path) else []


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


rows, support, outcomes = [], Counter(), []
for cell in cells:
    arms = {}
    for arm in ("off", "on"):
        ev = kv(f"{runs}/eval.{cell}.{arm}/eval.summary.tsv")
        st = kv(f"{runs}/run.{cell}.{arm}/call.mosaic_status.tsv")
        arms[arm] = (ev, st)
        for k in KEYS:
            rows.append((cell, arm, k, ev.get(k, "-")))
        rows.append((cell, arm, "caller_status", st.get("status", "-")))
        calls = {r["block_id"]: r for r in table(f"{runs}/run.{cell}.{arm}/call.block_genotypes.tsv")}
        # Amendment 1: a failed evaluation (no summary) can leave a partial block table; it contributes no blocks
        for b in (table(f"{runs}/eval.{cell}.{arm}/eval.blocks.tsv") if ev else []):
            c = calls.get(b["block_id"], {})
            exact = "exact" if b.get("excess") == "0" else "wrong"
            support[(cell, arm, c.get("genotype_support_v2x", "-"), c.get("quality_filter_v2x", "-"), b.get("kind", "-"), exact)] += 1
    (eo, so), (en, sn) = arms["off"], arms["on"]
    xo, xn = num(eo.get("excess_over_floor_total")), num(en.get("excess_over_floor_total"))
    bo, bn = num(eo.get("blocks_exact")), num(en.get("blocks_exact"))
    vo, vn = num(eo.get("sv_blocks_exact")), num(en.get("sv_blocks_exact"))
    if None in (xo, xn, bo, bn, vo, vn):
        outcome = "REFUSAL_OR_NO_CALL(" + ",".join(a for a, x in (("off", xo), ("on", xn)) if x is None) + ")"
    elif xn > xo or bn < bo or vn < vo:
        outcome = "MIXED" if xn < xo else "HARMFUL"
    elif xn < xo:
        outcome = "RECOVERED" if xn == 0 else "IMPROVED"
    else:
        mt = {}
        cmp_dir = f"{runs}/compare/m3f.matched.tsv" if os.path.exists(f"{runs}/compare/m3f.matched.tsv") else f"{runs}/control.check/m3f.matched.tsv"
        for r in table(cmp_dir):
            if r["cell"] == cell:
                mt[r["key"]] = r["changed"]
        changed = any(mt.get(k) == "1" for k in ("winner_unordered_rows_digest", "called_sequence_digest", "block_projection_digest"))
        outcome = "CHANGED_EQUAL" if changed else "UNCHANGED"
    ro, rn = RANK.get(so.get("status"), 2), RANK.get(sn.get("status"), 2)
    downgrade = rn > ro
    outcomes.append((cell, outcome, int(downgrade), so.get("status", "-"), sn.get("status", "-"), eo.get("excess_over_floor_total", "-"),
                     en.get("excess_over_floor_total", "-"), eo.get("blocks_exact", "-"), en.get("blocks_exact", "-"), eo.get("sv_blocks_exact", "-"),
                     en.get("sv_blocks_exact", "-")))
with open(f"{prefix}.accuracy.tsv", "w") as f:
    f.write("cell\tarm\tkey\tvalue\n")
    for r in rows:
        f.write("\t".join(str(x) for x in r) + "\n")
with open(f"{prefix}.support_by_exactness.tsv", "w") as f:
    f.write("cell\tarm\tgenotype_support_v2x\tquality_filter_v2x\tblock_kind\tblock_exact\tblocks\n")
    for k, v in sorted(support.items()):
        f.write("\t".join(k) + f"\t{v}\n")
with open(f"{prefix}.outcome.tsv", "w") as f:
    f.write("cell\toutcome\tstatus_downgrade\tstatus_off\tstatus_on\texcess_off\texcess_on\tblocks_exact_off\tblocks_exact_on\tsv_exact_off\tsv_exact_on\n")
    for r in outcomes:
        f.write("\t".join(str(x) for x in r) + "\n")
print(open(f"{prefix}.outcome.tsv").read(), end="")
if target == "control":
    bad = [r for r in outcomes if r[1] in ("HARMFUL", "MIXED") or r[2] == 1 or r[1] == "REFUSAL_OR_NO_CALL(on)"]
    print("control stop rule:", "STOP (a control worsened)" if bad else "PASS (no control worsened)")
    sys.exit(1 if bad else 0)
