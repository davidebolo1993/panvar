#!/usr/bin/env python3
"""M1-C: the explicit-sequence validation set. Every M1-B C4 10x path candidate, re-declared as two explicit sequences
spelled from its own allele paths (genotype_probe spell), with the same expected totals; plus controls:
  ctl_swap        the reported pair with homologues exchanged                (must equal `reported` bit for bit)
  ctl_dup_label   the reported pair again under another label                (must equal `reported` bit for bit)
  ctl_homozygous  the reported 3-unit haplotype twice (duplicate sequence)  (a homozygote: prior defined, no log 2)
  ctl_truth       the truth pair (off-panel)                                 (prior_support none; no total)

  validation_set.py --candidates M1B_CANDIDATES --spell PREFIX --out OUT_PREFIX
"""
import argparse


def main():
    ap = argparse.ArgumentParser()
    for x in ("candidates", "spell", "out"):
        ap.add_argument("--" + x, required=True)
    a = ap.parse_args()
    iv = {}
    for line in open(a.spell + ".spelled.tsv"):
        f = line.rstrip("\n").split("\t")
        if f[0] != "record":
            iv[f[0]] = (f[3], f[4])
    fa = a.spell + ".fa"
    rows = []
    for line in open(a.candidates):
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        label, expected = f[0], (f[3] if len(f) > 3 else "NA")
        rows.append((label, f"{label}.1", f"{label}.2", expected))
    rows += [("ctl_swap", "reported.2", "reported.1", "NA"), ("ctl_dup_label", "reported.1", "reported.2", "NA"),
             ("ctl_homozygous", "reported.1", "reported.1", "NA"), ("ctl_truth", "truth.1", "truth.2", "NA")]
    with open(a.out + ".explicit.tsv", "w") as out:
        out.write("# label\tfasta\trecord1\trecord2\tbegin1\tend1\tbegin2\tend2\texpected_total\n")
        for label, r1, r2, expected in rows:
            out.write(f"{label}\t{fa}\t{r1}\t{r2}\t{iv[r1][0]}\t{iv[r1][1]}\t{iv[r2][0]}\t{iv[r2][1]}\t{expected}\n")
    print(open(a.out + ".explicit.tsv").read(), end="")


if __name__ == "__main__":
    main()
