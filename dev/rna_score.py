#!/usr/bin/env python3
"""Scores recruit tables of simulated RNA reads (dev/rna_sim.py) against their truth.

    rna_score.py --locus <name> --truth <truth.tsv> --out <scores.tsv> <arm>:<class>:<table> [...]
Each table is recruit's output for one read class (short, kinnex, ont) in one arm. A fragment missing from
the table was not screened: unrecruited.

Rows (`group`), each with fragments, certain, ambiguous, certain/fragments and (certain+ambiguous)/fragments:
    locus:<perturbation>  the held-out haplotypes' molecules, by perturbation, and locus:all
    junction              unperturbed molecules whose reads cross a splice junction
    exonic                unperturbed molecules whose reads cross none
    whole                 long reads of whole unperturbed molecules (for ONT those not cut at the 5' end)
    paralog, pseudogene, transcriptome   negatives: certain is false certain
    precision             certain locus fragments over all certain fragments, in this mix of reads
A wide table of certain/fragments per arm goes to standard output.
"""
import argparse
import collections
import sys

PERTURB = ["annotated", "skip", "retained_intron", "pre_mrna", "alt_end"]
NEG = ["paralog", "pseudogene", "transcriptome"]


def fates(path, locus):
    out = {}
    with open(path) as f:
        head = f.readline().rstrip("\n").split("\t")
        i_frag, i_loc, i_fate = head.index("fragment"), head.index("locus"), head.index("fate")
        for line in f:
            r = line.rstrip("\n").split("\t")
            if r[i_loc] == locus:
                out[r[i_frag]] = r[i_fate]
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--locus", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("sets", nargs="+")
    a = ap.parse_args()

    truth = collections.defaultdict(list)   # class -> rows
    with open(a.truth) as f:
        head = f.readline().rstrip("\n").split("\t")
        for line in f:
            r = dict(zip(head, line.rstrip("\n").split("\t")))
            truth[r["class"]].append(r)

    rows = []
    wide = collections.OrderedDict()
    arms = []
    for spec in a.sets:
        arm, cls, table = spec.split(":", 2)
        if arm not in arms:
            arms.append(arm)
        fate = fates(table, a.locus)
        groups = collections.defaultdict(lambda: [0, 0, 0])   # fragments, certain, ambiguous

        def add(g, r):
            v = groups[g]
            v[0] += 1
            x = fate.get(r["read"], "unrecruited")
            v[1] += x == "certain"
            v[2] += x == "ambiguous"
        for r in truth[cls]:
            src = r["source"]
            if src == "locus":
                add("locus:" + r["perturbation"], r)
                add("locus:all", r)
                if r["perturbation"] == "annotated":
                    add("junction" if r["junction_read"] == "1" else "exonic", r)
                    if cls != "short" and r["full_length"] == "1":
                        add("whole", r)
            else:
                add(src, r)
        tp = groups["locus:all"][1]
        fp = sum(groups[n][1] for n in NEG)
        groups["precision"] = [tp + fp, tp, 0]
        order = ["locus:" + p for p in PERTURB] + ["locus:all", "junction", "exonic", "whole"] + NEG + ["precision"]
        for g in order:
            if g not in groups:
                continue
            n, c, amb = groups[g]
            rows.append([a.locus, cls, arm, g, n, c, amb, f"{c / n:.4f}" if n else ".",
                         f"{(c + amb) / n:.4f}" if n else "."])
            wide.setdefault((cls, g), {})[arm] = (c, n)

    with open(a.out, "w") as f:
        f.write("locus\tclass\tarm\tgroup\tfragments\tcertain\tambiguous\tcertain_share\trecruited_share\n")
        for r in rows:
            f.write("\t".join(map(str, r)) + "\n")
    w = sys.stdout
    w.write("class\tgroup\tn\t" + "\t".join(arms) + "\n")
    for (cls, g), v in wide.items():
        n = max(x[1] for x in v.values())
        cells = []
        for arm in arms:
            if arm not in v:
                cells.append(".")
            elif g in NEG:
                cells.append(str(v[arm][0]))   # false certain fragments
            else:
                cells.append(f"{v[arm][0] / v[arm][1]:.3f}" if v[arm][1] else ".")
        w.write(f"{cls}\t{g}\t{n}\t" + "\t".join(cells) + "\n")


if __name__ == "__main__":
    main()
