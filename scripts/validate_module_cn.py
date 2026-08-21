#!/usr/bin/env python3
"""Test the MODULE_BP linear model against INDEPENDENT copy-number truth.

`call` defines CN as round(hbp/u), so asking whether hbp is close to u*CN only asks whether the
measurements form an integer-like lattice -- it cannot tell you the lattice is in the right place. This
compares the raw dosage against an external truth instead, and never refits u: recalibrating the slope
against truth and then reporting the residual would hide exactly the bad unit we are looking for.

  u    = CN_UNIT_BP, the reference-calibrated unit (INFO)
  x_h  = FORMAT:CNR_RAW, this haplotype's hbp/u before rounding
  C_h  = truth copy number, from a pangene BED (optionally shifted by a constant offset)

Reported: exact integer agreement, raw error and its spread, the slope of x on C through the origin
(1.0 means the unit is right; below 1 means it is too large), a normalized model residual, and the same
broken out by true CN and by direction, because a unit that is wrong by a scale factor looks fine at the
CN it was calibrated on and drifts either side of it.

  validate_module_cn.py --vcf call.region.vcf --truth gstm.bed --genes GSTM1,GSTM2,GSTM4,GSTM5
                        [--offset auto|N] [--label NAME]
"""
import argparse
import collections
import math
import statistics
import sys


def read_truth(path, genes):
    """molecule -> number of annotated rows whose gene is in `genes`."""
    want = {g for g in genes.split(",") if g}
    per = collections.Counter()
    seen = set()
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 2:
                continue
            seen.add(f[0])
            if f[1] in want:
                per[f[0]] += 1
    return {m: per.get(m, 0) for m in seen}


def strip_hap(name):
    return name.split(":")[0] if ":" in name else name


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--genes", required=True)
    ap.add_argument("--offset", default="auto",
                    help="constant shift between module CN and the summed gene count "
                         "('auto' takes the mode, or give an integer)")
    ap.add_argument("--label", default="module")
    args = ap.parse_args()

    truth = read_truth(args.truth, args.genes)

    hdr, rec = None, None
    with open(args.vcf) as fh:
        for line in fh:
            if line.startswith("#CHROM"):
                hdr = line.rstrip("\n").split("\t")
                continue
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            inf = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
            if inf.get("CN_METHOD") != "MODULE_BP" or "CNR_RAW" not in f[8]:
                continue
            span = abs(int(inf.get("END", f[1])) - int(f[1]))
            if rec is None or span > rec[0]:
                rec = (span, f, inf)
    if rec is None:
        sys.exit("validate_module_cn: no MODULE_BP record with FORMAT:CNR_RAW in " + args.vcf)
    _span, f, inf = rec
    keys = f[8].split(":")
    ri, ci = keys.index("CNR_RAW"), keys.index("CN")
    unit = float(inf.get("CN_UNIT_BP", 0)) or float("nan")

    pairs = []   # (truth C, raw dosage x, integer CN)
    for i in range(9, len(hdr)):
        v = f[i].split(":")
        if len(v) <= ri or v[ri] == ".":
            continue
        key = strip_hap(hdr[i])
        # truth files key on the molecule name; try the full path name then the stripped form
        C = truth.get(hdr[i], truth.get(key))
        if C is None:
            continue
        pairs.append((C, float(v[ri]), int(v[ci])))
    if not pairs:
        sys.exit("validate_module_cn: no haplotype matched the truth file")

    if args.offset == "auto":
        off = collections.Counter(cn - C for C, _x, cn in pairs).most_common(1)[0][0]
    else:
        off = int(args.offset)
    # truth on the MODULE's scale
    data = [(C + off, x, cn) for C, x, cn in pairs]

    # An auto offset is taken from the CALLED CN, so it absorbs any constant error in REF_CN. What
    # remains is a test of relative dosage, not of the absolute anchor. Report the unshifted numbers
    # alongside so the two questions stay separate; --offset N makes the shift a stated assumption.
    if args.offset == "auto":
        exact0 = sum(1 for C, _x, cn in pairs if cn == C)
        print(f"### {args.label}   OFFSET IS AUTO ({off:+d}), taken from the called CN itself")
        print(f"  Against UNSHIFTED truth: integer agreement {exact0}/{len(pairs)} "
              f"({100.0*exact0/len(pairs):.1f}%)")
        if off != 0:
            print(f"  The {off:+d} shift below is therefore NOT validated -- it is assumed. Everything "
                  f"after this line tests relative dosage only.")
            print(f"  Re-run with --offset 0, or with a biologically justified constant, to test the "
                  f"absolute REF_CN anchor.")
        print()

    exact = sum(1 for C, _x, cn in data if cn == C)
    raw_err = [x - C for C, x, _cn in data]
    slope = (sum(x * C for C, x, _ in data) / sum(C * C for C, _x, _ in data)) if data else float("nan")
    num = sum((x * unit - unit * C) ** 2 for C, x, _ in data)
    den = sum((unit * C) ** 2 for C, _x, _ in data)
    nres = math.sqrt(num / den) if den else float("nan")

    print(f"### {args.label}   n={len(data)}   unit={unit:.0f}   truth offset {off:+d}")
    print(f"  integer agreement   {exact}/{len(data)} ({100.0*exact/len(data):.1f}%)")
    print(f"  raw error (x - C)   mean {statistics.mean(raw_err):+.3f}   "
          f"median {statistics.median(raw_err):+.3f}   "
          f"max|.| {max(abs(e) for e in raw_err):.3f}")
    print(f"  slope through 0     {slope:.4f}   "
          f"({'unit too large' if slope < 0.98 else 'unit too small' if slope > 1.02 else 'unit consistent'})")
    print(f"  normalized residual {nres:.4f}")
    by = collections.defaultdict(list)
    for C, x, cn in data:
        by[C].append((x, cn))
    print("  by true CN:")
    for C in sorted(by):
        xs = [x for x, _ in by[C]]
        ok = sum(1 for _x, cn in by[C] if cn == C)
        print(f"    C={C}  n={len(by[C]):4d}  dosage {min(xs):.3f}..{max(xs):.3f}  "
              f"mean error {statistics.mean(x - C for x in xs):+.3f}  integer ok {ok}/{len(by[C])}")


main()
