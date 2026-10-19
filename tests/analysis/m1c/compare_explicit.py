#!/usr/bin/env python3
"""M1-C gate: explicit sequences spelled from panel paths must reproduce the path-based scorer exactly.

Compared, between the path-based run (PATH_PREFIX) and the explicit run (EXPLICIT_PREFIX), both score-set outputs:
  * per candidate label present in both: fragment and exposure terms (hexfloat, bit for bit), prior and total (17 digits);
  * per unique haplotype of the explicit run, matched to the path run by (sequence FNV-64, length): exposure (hexfloat),
    the fragment mass column (FNV-64 of its raw doubles) and fragments with mass;
  * controls: ctl_swap and ctl_dup_label equal `reported` bit for bit; ctl_homozygous has a defined prior and a finite
    likelihood; ctl_truth has prior_support none and no total.
Exit 0 only if every check passes.

  compare_explicit.py PATH_PREFIX EXPLICIT_PREFIX OUT.tsv
"""
import csv
import sys


def rows(path):
    return {r["label"]: r for r in csv.DictReader(open(path), delimiter="\t")}


def main(path_prefix, explicit_prefix, out):
    p, e = rows(path_prefix + ".scores_truth_free.tsv"), rows(explicit_prefix + ".scores_truth_free.tsv")
    hp = {(r["sequence_fnv64"], r["bp"]): r for r in csv.DictReader(open(path_prefix + ".haplotypes.tsv"), delimiter="\t")}
    he = list(csv.DictReader(open(explicit_prefix + ".haplotypes.tsv"), delimiter="\t"))
    checks = []
    for label in sorted(set(p) & set(e)):
        for key in ("fragment_term_hex", "exposure_term_hex", "prior", "total"):
            checks.append((f"candidate {label} {key}", p[label][key] == e[label][key], f"{p[label][key]} | {e[label][key]}"))
    matched = 0
    for r in he:
        k = (r["sequence_fnv64"], r["bp"])
        if k not in hp:
            continue
        matched += 1
        for key in ("exposure_hex", "mass_column_fnv64", "fragments_with_mass"):
            checks.append((f"haplotype {k[0]} {key}", hp[k][key] == r[key], f"{hp[k][key]} | {r[key]}"))
    # every explicit haplotype except the two off-panel truth haplotypes must have its path-run counterpart
    checks.append(("explicit haplotypes matched to path haplotypes", matched == len(he) - 2, f"{matched} of {len(he)}"))
    rep = e["reported"]
    for ctl in ("ctl_swap", "ctl_dup_label"):
        for key in ("fragment_term_hex", "exposure_term_hex", "prior", "total"):
            checks.append((f"{ctl} {key} equals reported", e[ctl][key] == rep[key], f"{e[ctl][key]} | {rep[key]}"))
    h = e["ctl_homozygous"]
    checks.append(("ctl_homozygous prior defined", h["prior_support"] == "defined", h["prior"]))
    checks.append(("ctl_homozygous likelihood finite", h["likelihood"] not in ("-inf", "nan", "inf"), h["likelihood"]))
    t = e["ctl_truth"]
    checks.append(("ctl_truth prior_support none", t["prior_support"] == "none", t["prior"]))
    checks.append(("ctl_truth total NA", t["total"] == "NA", t["total"]))
    with open(out, "w") as f:
        f.write("check\tpass\tdetail\n")
        for name, ok, detail in checks:
            f.write(f"{name}\t{ok}\t{detail}\n")
    failed = [c for c in checks if not c[1]]
    print(f"{len(checks)} checks, {len(failed)} failed")
    for c in failed:
        print("FAILED:", c[0], c[2])
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:4]))
