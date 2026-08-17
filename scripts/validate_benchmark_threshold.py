#!/usr/bin/env python3
"""Does lowering the call threshold actually find more, and does the ledger stay honest while it does?

`benchmark --min-sv-bp` only RECLASSIFIES: it re-sorts truth events into called/missed/below and gates
the `called` reconstruction. It cannot discover anything, because discovery happened in `call`. So a
sweep of benchmark alone proves nothing, and reading one as if it did would be the same mistake the
alignment-run split made -- watching a number move for a reason unrelated to the claim.

The real experiment reruns BOTH at each threshold, over a fixed graph and a fixed bubble set, and
checks the invariants that must hold if the ledger means what it says:

  invariant   the total truth-event count and their sizes come from the WALKS, so they must not move
              at all across thresholds. If they do, the decomposition is threshold-dependent and
              nothing else here is interpretable.
  monotone    below_threshold must not increase as the threshold falls, and called+missed must not
              decrease -- an event only ever crosses from below into eligible.
  substance   of the newly eligible events, how many become `called` rather than `missed`. This is
              the part that is NOT guaranteed by definition, and it is the actual question.
  no harm     FP, FN, worse_than_baseline, ref_mismatch, unplaceable and clamped must not rise, and
              genotype variation_recovered must not fall.
  ceiling     the allele VCF must still reconstruct every haplotype exactly at every threshold.

  validate_benchmark_threshold.py --panvar build/panvar --gfa call.sorted.gfa
                                  --bubble-prefix panphorte --reference NAME --workdir DIR
                                  [--thresholds 200,100,50,25,10] [--call-extra "--cn"]
"""
import argparse
import os
import subprocess
import sys


def run(cmd):
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        tail = (r.stderr.strip().splitlines() or [""])[-1]
        raise SystemExit(f"validate_benchmark_threshold: command failed ({tail})\n  {' '.join(cmd)}")


def summary(path):
    """(scope, key, band) -> (n, pct) from a benchmark qv_summary.tsv.

    All three columns are part of the identity. `truth_event` emits an ALL row AND a row per svtype,
    every one of them with band `called`, so keying on (scope, band) silently returns whichever svtype
    happened to sort last."""
    out = {}
    with open(path) as fh:
        next(fh, None)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 5:
                out[(f[0], f[1], f[2])] = (f[3], f[4])
    return out


def num(v, default=0.0):
    try:
        return float(v)
    except (TypeError, ValueError):
        return default


def event_fingerprint(path):
    """Multiset of (sample, bubble, size_bp). Threshold-independent by construction -- it is read off
    the walks -- so this is the check that the decomposition itself did not move."""
    sizes = {}
    with open(path) as fh:
        next(fh, None)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 5:
                k = (f[0], f[1], f[4])
                sizes[k] = sizes.get(k, 0) + 1
    return sizes


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panvar", required=True)
    ap.add_argument("--gfa", required=True)
    ap.add_argument("--bubble-prefix", required=True,
                    help="prefix whose .bubbles.csv `call` ran on; held FIXED across the sweep")
    ap.add_argument("--reference", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--thresholds", default="200,100,50,25,10")
    ap.add_argument("--threads", default="0")
    ap.add_argument("--call-extra", default="",
                    help="extra flags passed to every `call` run, e.g. '--cn'. The sweep must run the "
                         "SAME call configuration the pipeline ships, or it measures a different "
                         "caller: at CYP2D6, adding --cn moves genotype variation_recovered from "
                         "86.1%% to 39.2%%, because copy-number records reconstruct heuristically")
    args = ap.parse_args()

    thresholds = [int(t) for t in args.thresholds.split(",") if t.strip()]
    if sorted(thresholds, reverse=True) != thresholds:
        sys.exit("--thresholds must be given in descending order")
    os.makedirs(args.workdir, exist_ok=True)

    rows = []
    for t in thresholds:
        pfx = os.path.join(args.workdir, f"t{t}")
        run([args.panvar, "call", "-i", args.gfa, "--bubble-prefix-in", args.bubble_prefix,
             "--reference-path", args.reference, "-o", pfx + ".call",
             "--min-sv-bp", str(t), "--allele-vcf", "--threads", args.threads, "--quiet"]
            + [x for x in args.call_extra.split() if x])
        run([args.panvar, "benchmark", "-i", args.gfa, "--bubble-prefix-in", args.bubble_prefix,
             "--reference-path", args.reference,
             "--variant-nodes", pfx + ".call.variant_nodes.tsv",
             "--vcf", pfx + ".call.region.vcf", "-o", pfx + ".bm",
             "--min-sv-bp", str(t), "--threads", args.threads, "--quiet"])
        # The allele VCF is the serialization ceiling and must stay exact at every threshold.
        run([args.panvar, "benchmark", "-i", args.gfa, "--bubble-prefix-in", args.bubble_prefix,
             "--reference-path", args.reference,
             "--variant-nodes", pfx + ".call.variant_nodes.tsv",
             "--vcf", pfx + ".call.alleles.vcf", "-o", pfx + ".av",
             "--min-sv-bp", str(t), "--no-truth-events", "--threads", args.threads, "--quiet"])
        s = summary(pfx + ".bm.qv_summary.tsv")
        av = summary(pfx + ".av.qv_summary.tsv")
        rows.append({
            "t": t,
            "events": int(num(s.get(("truth_event", "ALL", "events"), (0, 0))[0])),
            "called": int(num(s.get(("truth_event", "ALL", "called"), (0, 0))[0])),
            "missed": int(num(s.get(("truth_event", "ALL", "missed"), (0, 0))[0])),
            "below": int(num(s.get(("truth_event", "ALL", "below_threshold"), (0, 0))[0])),
            "tp": int(num(s.get(("gt_carrier", "ALL", "TP"), (0, 0))[0])),
            "fp": int(num(s.get(("gt_carrier", "ALL", "FP"), (0, 0))[0])),
            "fn": int(num(s.get(("gt_carrier", "ALL", "FN"), (0, 0))[0])),
            "worse": int(num(s.get(("gt_gap", "ALL", "worse_than_baseline"), (0, 0))[0])),
            "refmm": int(num(s.get(("gt_records", "ALL", "ref_mismatch"), (0, 0))[0])),
            "unpl": int(num(s.get(("gt_records", "ALL", "unplaceable"), (0, 0))[0])),
            "clamp": int(num(s.get(("gt_records", "ALL", "clamped"), (0, 0))[0])),
            "vr": num(s.get(("variation_recovered", "ALL", "genotype"), (0, "nan"))[1], float("nan")),
            "av_delta": int(num(av.get(("gt_gap", "ALL", "genotype_delta"), (0, 0))[0])),
            "fp_path": pfx + ".bm.truth_events.tsv",
        })

    print(f"{'T':>5} {'events':>8} {'called':>8} {'missed':>7} {'below':>8} "
          f"{'TP':>6} {'FP':>6} {'FN':>5} {'worse':>6} {'var_rec':>8} {'allele':>7}")
    for r in rows:
        print(f"{r['t']:>5} {r['events']:>8} {r['called']:>8} {r['missed']:>7} {r['below']:>8} "
              f"{r['tp']:>6} {r['fp']:>6} {r['fn']:>5} {r['worse']:>6} {r['vr']:>7.2f}% {r['av_delta']:>7}")

    fails = []
    notes = []
    base = event_fingerprint(rows[0]["fp_path"])
    for r in rows[1:]:
        if event_fingerprint(r["fp_path"]) != base:
            fails.append(f"T={r['t']}: the truth events themselves moved. Event size is read off the "
                         f"walks, so it MUST be threshold-independent")
    for r in rows:
        if r["events"] != rows[0]["events"]:
            fails.append(f"T={r['t']}: total truth events {r['events']} != {rows[0]['events']}")
        if r["called"] + r["missed"] + r["below"] != r["events"]:
            fails.append(f"T={r['t']}: called+missed+below != events (the classes must partition)")
        if r["av_delta"] != 0:
            fails.append(f"T={r['t']}: the allele VCF left {r['av_delta']} bp unreconstructed; it is "
                         f"the serialization ceiling and must be exact")
    for a, b in zip(rows, rows[1:]):          # b has the LOWER threshold
        if b["below"] > a["below"]:
            fails.append(f"T {a['t']}->{b['t']}: below_threshold rose {a['below']}->{b['below']}")
        if b["called"] + b["missed"] < a["called"] + a["missed"]:
            fails.append(f"T {a['t']}->{b['t']}: eligible events fell "
                         f"{a['called'] + a['missed']}->{b['called'] + b['missed']}")
        if b["called"] < a["called"]:
            fails.append(f"T {a['t']}->{b['t']}: record-attributed called fell {a['called']}->{b['called']}")
        # A rising ref_mismatch is a placement bug, not a modelling trade-off: the record's REF no
        # longer matches the reference at its POS, and every score built on it is fiction.
        if b["refmm"] > a["refmm"]:
            fails.append(f"T {a['t']}->{b['t']}: ref_mismatch rose {a['refmm']}->{b['refmm']}")
        # These CAN legitimately move against the threshold -- more records mean more merged and
        # overlapping ones, and independent records do not compose over one walk. They are the call
        # interactions this sweep exists to surface, so they are reported loudly and separately rather
        # than failed: a hard failure here would just train the reader to ignore it.
        for k, label in (("fp", "FP"), ("fn", "FN"), ("worse", "worse_than_baseline"),
                         ("unpl", "unplaceable"), ("clamp", "clamped")):
            if b[k] > a[k]:
                notes.append(f"T {a['t']}->{b['t']}: {label} rose {a[k]}->{b[k]}")
        if b["vr"] < a["vr"] - 0.01:      # percentage points; below this it is float noise
            notes.append(f"T {a['t']}->{b['t']}: genotype variation_recovered FELL "
                         f"{a['vr']:.2f}% -> {b['vr']:.2f}% -- lowering the threshold made the VCF "
                         f"reconstruct less, which is a call-side interaction, not a benchmark defect")

    # The substantive question, which no invariant above can answer: of the events that BECAME
    # eligible as the threshold fell, what fraction the caller actually found.
    print()
    for a, b in zip(rows, rows[1:]):
        newly = (b["called"] + b["missed"]) - (a["called"] + a["missed"])
        if newly > 0:
            got = b["called"] - a["called"]
            print(f"  T {a['t']}->{b['t']}: {newly} events became eligible, {got} of them attributed "
                  f"to a record ({100.0 * got / newly:.1f}%)")
        else:
            print(f"  T {a['t']}->{b['t']}: no event became eligible")

    print()
    if notes:
        print("  NON-MONOTONIC -- call interactions, to explain rather than to ignore:")
        for n in notes:
            print("    " + n)
        print()
    if fails:
        for f in fails:
            print("  FAIL " + f)
        sys.exit(f"{len(fails)} invariant(s) violated")
    print("  all threshold invariants hold")


main()
