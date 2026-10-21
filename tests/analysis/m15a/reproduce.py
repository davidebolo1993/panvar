#!/usr/bin/env python3
"""M15-A archived-evidence reproduction (no caller run): the M14-A full-panel runs (built from the current sources; panvar
relinked at M14-A's build unit) against M11-D's archived runs (M11-C build), for the 15 shared cells: rank-1 score, block
paths, haplotype SHA-256, status, model fingerprint, and the release-candidate output files' presence."""
import csv
import os
import sys

R = "/Users/davide.bolognini/panvar/tests/results"
CELLS = ["cyp.HG02769", "cyp.HG02953", "cyp.HG03456", "cyp.NA19650", "cyp.NA19705", "cyp.NA19983", "c4.NA19705", "c4.NA19317",
         "c4.HG03065", "c4.HG02018", "c4.HG00513", "c4.NA12329", "acot.HG03209", "acot.HG01975", "acot.NA19909"]
FILES = ("locus_genotype.tsv", "block_genotypes.tsv", "sequence_locus.tsv", "sequence_blocks.tsv", "sequence_haplotypes.tsv",
         "representation_intervals.tsv", "sequence_alternatives.tsv", "sequence_compound.tsv", "dosage_status.tsv", "settings.tsv")


def key(p):
    r = next(x for x in csv.DictReader(open(p + ".equivalence_set.tsv"), delimiter="\t") if x["rank"] == "1" and x["representative"] == "true")
    h = {x["haplotype"]: x["sha256"] for x in csv.DictReader(open(p + ".sequence_haplotypes.tsv"), delimiter="\t")}
    s = {l.split("\t")[0]: l.rstrip("\n").split("\t")[1] for l in open(p + ".mosaic_status.tsv") if "\t" in l}
    return (r["score"], r["hap1_alleles"], r["hap2_alleles"], tuple(sorted(h.items())), s.get("status"), s.get("model_fingerprint"))


def main(out):
    rows, ok = [], 0
    for c in CELLS:
        a, b = f"{R}/m11d_dev_pilot_20261020/runs/call.{c}/call", f"{R}/m14a_representative_panel_20261021/runs/call.{c}.full/call"
        same = key(a) == key(b)
        files = all(os.path.exists(f"{b}.{f}") for f in FILES)
        ok += same and files
        rows.append(f"{c}\t{same}\t{files}\n")
    open(f"{out}/reproduction.tsv", "w").writelines(["cell\tidentical_to_m11d\trelease_outputs_present\n"] + rows)
    print("".join(rows) + f"{ok} of {len(CELLS)} identical with every release-candidate output present")
    return 0 if ok == len(CELLS) else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
