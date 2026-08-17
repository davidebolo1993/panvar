#!/usr/bin/env python3
"""Is the panel-spacing copy-number step a property of the LOCUS or of the COHORT?

`call --cn-unit-spacing` takes the per-copy step from the gaps between the panel's own copy-state
clusters. That is measurably better calibrated than the reference-derived unit, but it makes the
answer depend on who is in the panel -- and unlike the reference ratio, nothing about it is fixed by
the reference. This resamples the cohort and re-measures.

Two experiments, both driven from outside the caller by rewriting the panel rather than by
instrumenting the estimator:

  subsample     keep a fraction of the non-reference haplotypes, R times. A step that is a property
                of the locus should barely move; one that is a property of the cohort will.
  reference     re-run with a different path as the reference. CN is anchored on REF_CN, so this
                changes what "reference-like" means and is the other half of the stability question.

Reports the step's spread, and -- more importantly -- how many per-haplotype CN calls change, since
that is what a downstream consumer actually sees.

  validate_spacing_stability.py --panvar build/panvar --gfa norm.sorted.gfa --bubbles b.bubbles.csv
                                --reference NAME [--replicates 12] [--fraction 0.7] [--workdir DIR]
"""
import argparse
import collections
import os
import random
import statistics
import subprocess
import sys


def read_paths(gfa):
    names = []
    for line in open(gfa):
        if line.startswith("P\t") or line.startswith("W\t"):
            names.append(line.split("\t")[1])
    return names


def write_subset(gfa, out, keep):
    """Copy the GFA keeping only P/W lines whose name is in `keep`. S/L lines are untouched, so the
    graph is unchanged and only the panel differs -- which is the variable under test."""
    kept = 0
    with open(out, "w") as fh:
        for line in open(gfa):
            if line.startswith(("P\t", "W\t")):
                if line.split("\t")[1] in keep:
                    fh.write(line); kept += 1
            else:
                fh.write(line)
    return kept


def call_once(panvar, gfa, bubbles, ref, prefix, spacing=True):
    cmd = [panvar, "call", "-i", gfa, "-c", bubbles, "-r", ref, "-o", prefix, "--cn", "-q"]
    if spacing:
        cmd.append("--cn-unit-spacing")
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        return None, (r.stderr.strip().splitlines() or [""])[-1]
    return prefix + ".region.vcf", None


def module_records(vcf):
    """record id -> (step_bp, {sample: cn})"""
    out, hdr = {}, None
    for line in open(vcf):
        if line.startswith("#CHROM"):
            hdr = line.rstrip("\n").split("\t"); continue
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        inf = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
        if inf.get("CN_METHOD") != "MODULE_BP":
            continue
        keys = f[8].split(":")
        if "CN" not in keys:
            continue
        ci = keys.index("CN")
        cn = {}
        for i in range(9, len(hdr)):
            v = f[i].split(":")
            if len(v) > ci and v[ci] != ".":
                cn[hdr[i]] = int(v[ci])
        out[f[2]] = (inf.get("CN_STEP_BP"), inf.get("CN_STEP_SUPPORT"), cn)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panvar", required=True)
    ap.add_argument("--gfa", required=True)
    ap.add_argument("--bubbles", required=True)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--replicates", type=int, default=12)
    ap.add_argument("--fraction", type=float, default=0.7)
    ap.add_argument("--workdir", default=None)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    work = args.workdir or os.path.join(os.path.dirname(args.gfa) or ".", "spacing_stability")
    os.makedirs(work, exist_ok=True)
    paths = read_paths(args.gfa)
    others = [p for p in paths if p != args.reference]
    if not others:
        sys.exit("validate_spacing_stability: the GFA has no non-reference paths")

    full_vcf, err = call_once(args.panvar, args.gfa, args.bubbles, args.reference,
                              os.path.join(work, "full"))
    if full_vcf is None:
        sys.exit(f"validate_spacing_stability: the full panel does not support spacing mode: {err}")
    base = module_records(full_vcf)
    if not base:
        sys.exit("validate_spacing_stability: no MODULE_BP record in the full-panel run")

    print(f"### full panel   {len(paths)} haplotypes")
    for rid, (step, support, cn) in sorted(base.items()):
        print(f"  {rid:24s} CN_STEP_BP={step}  support(clusters,gaps,dropped)={support}  n_cn={len(cn)}")

    rng = random.Random(args.seed)
    k = max(4, int(args.fraction * len(others)))
    steps = collections.defaultdict(list)
    flips = collections.Counter()
    compared = collections.Counter()
    refused = 0
    for r in range(args.replicates):
        keep = set(rng.sample(others, k)) | {args.reference}
        sub = os.path.join(work, f"sub{r}.gfa")
        write_subset(args.gfa, sub, keep)
        vcf, err = call_once(args.panvar, sub, args.bubbles, args.reference,
                             os.path.join(work, f"sub{r}"))
        if vcf is None:
            refused += 1          # a thinner panel can legitimately stop supporting the model
            continue
        for rid, (step, _sup, cn) in module_records(vcf).items():
            if step is not None:
                steps[rid].append(int(step))
            if rid in base:
                for s, c in cn.items():
                    if s in base[rid][2]:
                        compared[rid] += 1
                        if base[rid][2][s] != c:
                            flips[rid] += 1

    print(f"\n### {args.replicates} replicates at {args.fraction:.0%} of the panel "
          f"({k} of {len(others)} non-reference haplotypes)")
    if refused:
        print(f"  {refused} replicate(s) REFUSED spacing mode outright -- the thinner panel could not "
              f"support an estimate, which is the model declining rather than guessing")
    for rid in sorted(base):
        vals = steps.get(rid, [])
        if not vals:
            print(f"  {rid:24s} no replicate produced a step")
            continue
        lo, hi = min(vals), max(vals)
        full_step = int(base[rid][0]) if base[rid][0] else 0
        spread = 100.0 * (hi - lo) / full_step if full_step else float("nan")
        n = compared[rid]
        line = f"  {rid:24s} step {lo}..{hi} (full {full_step}, spread {spread:.1f}%)"
        if n:
            line += (f"   CN changed on {flips[rid]}/{n} haplotype-replicate pairs "
                     f"({100.0 * flips[rid] / n:.2f}%)")
        print(line)

    print("\n### reference switch")
    alt = None
    for p in others:
        if module_records(full_vcf):
            alt = p
            break
    if alt is None:
        print("  no alternative reference available")
        return
    vcf, err = call_once(args.panvar, args.gfa, args.bubbles, alt, os.path.join(work, "altref"))
    if vcf is None:
        print(f"  switching the reference to {alt} makes spacing unsupported: {err}")
        return
    alt_rec = module_records(vcf)
    print(f"  reference {args.reference}  ->  {alt}")
    print(f"  MODULE_BP records {len(base)} -> {len(alt_rec)}"
          f"   (record ids are reference-relative, so a direct per-record match is not expected)")
    for rid, (step, support, _cn) in sorted(alt_rec.items()):
        print(f"    {rid:24s} CN_STEP_BP={step}  support={support}")


main()
