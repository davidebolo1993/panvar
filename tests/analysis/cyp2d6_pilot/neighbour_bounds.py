#!/usr/bin/env python3
"""CYP2D6/D7 pilot, AMENDMENT_3: the repair neighbourhood bound of every frozen campaign model, from the model alone.

The caller's bound (src/genotype_repair.cpp, repair_neighbour_bound) is
    sum over blocks b of  2 * min(repair_local_pairs, a_b * (a_b + 1) / 2),   a_b = the model's alleles in block b.
It reads no reads and no truth. a_b is taken here from `genotype_probe panel-alleles` on a curation-equivalent copy of
each model (its NA18989#1 P line removed; panel-alleles takes no --panel-curation), after checking that the copy's
preflight fingerprint equals the curated model's. a_b is the number of distinct allele indices the active templates
carry in block b; the indices must be exactly 0..a_b-1 and their sum over blocks must equal the curated preflight's
canonical_alleles, so no model allele is unaccounted for.

Usage: neighbour_bounds.py <pilot package> <genotype_probe> <scratch dir> <repair_local_pairs>
Writes <pilot package>/model_neighbour_bounds.tsv (one row per model and block, with the model total on every row)
and prints the maximum.
"""
import csv
import os
import subprocess
import sys

REF = "grch38#1#chr22:42031864-42245566"
CURATED = "NA18989#1#haplotype1-0000031:36617797-36825011"


def kv(path):
    out = {}
    for line in open(path):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 2:
            out[p[0]] = p[1]
    return out


def main(pkg, probe, scratch, local_pairs):
    os.makedirs(scratch, exist_ok=True)
    with open(os.path.join(pkg, "data", "models.tsv")) as f:
        models = list(csv.DictReader(f, delimiter="\t"))
    rows, totals = [], {}
    for m in models:
        src = os.path.join(pkg, m["gfa"])
        copy = os.path.join(scratch, f"{m['model']}.curated_equivalent.gfa")
        removed = 0
        with open(src) as fin, open(copy, "w") as fout:
            for line in fin:
                if line.startswith("P\t") and line.split("\t")[1] == CURATED:
                    removed += 1
                    continue
                fout.write(line)
        if removed != 1:
            sys.exit(f"{m['model']}: expected to remove the one curated path, removed {removed}")
        pf = copy + ".preflight.tsv"
        subprocess.run([probe, "preflight", "--gfa", copy, "--reference", REF, "--out", pf], check=True,
                       capture_output=True)
        p = kv(pf)
        curated_pf = kv(os.path.join(pkg, "data", "graphs", f"{m['model']}.preflight.tsv"))
        if p["model_fingerprint"] != m["model_fingerprint"] or p["model_fingerprint"] != curated_pf["model_fingerprint"]:
            sys.exit(f"{m['model']}: the curation-equivalent copy has fingerprint {p['model_fingerprint']}, "
                     f"not the curated model's {m['model_fingerprint']}")
        pa = copy + ".panel_alleles.tsv"
        subprocess.run([probe, "panel-alleles", "--gfa", copy, "--reference", REF, "--out", pa], check=True,
                       capture_output=True)
        alleles = [l.rstrip("\n").split("\t")[1].split(",") for l in open(pa)
                   if l.strip() and not l.startswith(("#", "template"))]
        nblocks = len(alleles[0])
        if nblocks != int(curated_pf["blocks"]) or len(alleles) != int(curated_pf["templates"]):
            sys.exit(f"{m['model']}: panel-alleles shape {len(alleles)}x{nblocks} differs from the curated preflight")
        counts = []
        for b in range(nblocks):
            used = {int(r[b]) for r in alleles}
            if used != set(range(len(used))):
                sys.exit(f"{m['model']}: block {b} allele indices are not 0..a-1: {sorted(used)}")
            counts.append(len(used))
        if sum(counts) != int(curated_pf["canonical_alleles"]):
            sys.exit(f"{m['model']}: sum of block alleles {sum(counts)} != canonical_alleles "
                     f"{curated_pf['canonical_alleles']}")
        contrib = [2 * min(local_pairs, a * (a + 1) // 2) for a in counts]
        total = sum(contrib)
        totals[m["model"]] = total
        for b, (a, c) in enumerate(zip(counts, contrib)):
            rows.append({"model": m["model"], "model_fingerprint": m["model_fingerprint"], "blocks": nblocks,
                         "block_index": b, "block_alleles": a, "unordered_allele_pairs": a * (a + 1) // 2,
                         "repair_local_pairs": local_pairs, "block_contribution": c, "total_bound": total})
    with open(os.path.join(pkg, "model_neighbour_bounds.tsv"), "w") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    for k, v in totals.items():
        print(f"{k}\t{v}")
    print(f"max\t{max(totals.values())}")


if __name__ == "__main__":
    main(os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2]), os.path.abspath(sys.argv[3]), int(sys.argv[4]))
