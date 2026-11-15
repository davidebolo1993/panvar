#!/usr/bin/env python3
"""Scores recruit tables of simulated RNA reads (dev/rna_sim.py) against their truth.

    rna_score.py score --locus <name> --truth <truth.tsv> --out <scores.tsv> <arm>:<class>:<table> [...]
    rna_score.py merge --out <merged.tsv> <scores.tsv> [...]
`score`: each table is recruit's output (plain or gzipped) for one read class (short, kinnex, ont) in one
arm. A fragment missing from the table was not screened: unrecruited. Rows (`group`), each with
fragments, certain, ambiguous, certain/fragments and (certain+ambiguous)/fragments:
    locus:<perturbation>  the held-out haplotypes' projected transcripts and the molecules made from them,
                          by perturbation, and locus:all
    junction, exonic      unperturbed molecules whose reads cross a splice junction, or none
    whole                 long reads of whole unperturbed molecules (for ONT those not cut at the 5' end)
    cat:all, cat:junction, cat:whole   the same for transcripts of the haplotypes' own annotations, cut
                          by their coordinates (a second truth, kept out of locus:all and precision)
    paralog, paralog_hap, pseudogene, transcriptome   negatives: certain is false certain
    precision             certain locus fragments over all certain fragments, in this mix of reads
    gene:<name>           locus reads of one gene; cat_gene:<name> and neg:<source>:<name> likewise
`merge`: per-seed score tables pooled (fragments summed) with the mean, lowest and highest share over
the seeds.
    rna_score.py table <merged.tsv> [...]
`table`: the main groups of merged tables, one row per locus, class and group, one column per arm: the
pooled share with the seeds' range, and for negatives the false certain fragments.
"""
import argparse
import collections
import gzip
import sys

PERTURB = ["annotated", "skip", "retained_intron", "pre_mrna", "alt_end"]
NEG = ["paralog", "paralog_hap", "pseudogene", "transcriptome"]


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def fates(path, locus):
    out = {}
    with opener(path) as f:
        head = f.readline().rstrip("\n").split("\t")
        i_frag, i_loc, i_fate = head.index("fragment"), head.index("locus"), head.index("fate")
        for line in f:
            r = line.rstrip("\n").split("\t")
            if r[i_loc] == locus:
                out[r[i_frag]] = r[i_fate]
    return out


def score(a):
    truth = collections.defaultdict(list)   # class -> rows
    with opener(a.truth) as f:
        head = f.readline().rstrip("\n").split("\t")
        for line in f:
            r = dict(zip(head, line.rstrip("\n").split("\t")))
            truth[r["class"]].append(r)

    rows, wide, arms = [], collections.OrderedDict(), []
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
                add("gene:" + r["gene"], r)
                if r["perturbation"] == "annotated":
                    add("junction" if r["junction_read"] == "1" else "exonic", r)
                    if cls != "short" and r["full_length"] == "1":
                        add("whole", r)
            elif src == "cat":
                add("cat:all", r)
                add("cat_gene:" + r["gene"], r)
                if r["junction_read"] == "1":
                    add("cat:junction", r)
                if cls != "short" and r["full_length"] == "1":
                    add("cat:whole", r)
            else:
                add(src, r)
                add(f"neg:{src}:{r['gene']}", r)
        tp = groups["locus:all"][1]
        fp = sum(groups[n][1] for n in NEG)
        groups["precision"] = [tp + fp, tp, 0]
        main = (["locus:" + p for p in PERTURB] + ["locus:all", "junction", "exonic", "whole",
                                                    "cat:all", "cat:junction", "cat:whole"] + NEG + ["precision"])
        order = [g for g in main if g in groups] + sorted(g for g in groups if g not in main)
        for g in order:
            n, c, amb = groups[g]
            rows.append([a.locus, cls, arm, g, n, c, amb, f"{c / n:.4f}" if n else ".",
                         f"{(c + amb) / n:.4f}" if n else "."])
            if g in main:
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


def merge(a):
    acc = collections.OrderedDict()
    for path in a.tables:
        with open(path) as f:
            head = f.readline().rstrip("\n").split("\t")
            for line in f:
                r = dict(zip(head, line.rstrip("\n").split("\t")))
                key = (r["locus"], r["class"], r["arm"], r["group"])
                v = acc.setdefault(key, [0, 0, 0, []])
                n, c, amb = int(r["fragments"]), int(r["certain"]), int(r["ambiguous"])
                v[0] += n
                v[1] += c
                v[2] += amb
                if n:
                    v[3].append(c / n)
    with open(a.out, "w") as f:
        f.write("locus\tclass\tarm\tgroup\tfragments\tcertain\tambiguous\tcertain_share\tseeds\tseed_mean\t"
                "seed_min\tseed_max\n")
        for (locus, cls, arm, g), (n, c, amb, shares) in acc.items():
            m = f"{sum(shares) / len(shares):.4f}" if shares else "."
            lo = f"{min(shares):.4f}" if shares else "."
            hi = f"{max(shares):.4f}" if shares else "."
            share = f"{c / n:.4f}" if n else "."
            f.write(f"{locus}\t{cls}\t{arm}\t{g}\t{n}\t{c}\t{amb}\t{share}\t{len(shares)}\t{m}\t{lo}\t{hi}\n")


MAIN = (["locus:" + p for p in PERTURB] + ["locus:all", "junction", "exonic", "whole", "cat:all",
                                           "cat:junction", "cat:whole"] + NEG + ["precision"])


def table(a):
    cells, arms = collections.OrderedDict(), []
    for path in a.tables:
        with open(path) as f:
            head = f.readline().rstrip("\n").split("\t")
            for line in f:
                r = dict(zip(head, line.rstrip("\n").split("\t")))
                if r["group"] not in MAIN:
                    continue
                if r["arm"] not in arms:
                    arms.append(r["arm"])
                key = (r["locus"], r["class"], r["group"])
                n, c = int(r["fragments"]), int(r["certain"])
                if r["group"] in NEG:
                    v = f"{c}/{n}"
                elif n:
                    v = f"{c / n:.3f} ({float(r['seed_min']):.3f}-{float(r['seed_max']):.3f})"
                else:
                    v = "."
                cells.setdefault(key, {})[r["arm"]] = v
    order = {g: i for i, g in enumerate(MAIN)}
    print("locus\tclass\tgroup\t" + "\t".join(arms))
    for key in sorted(cells, key=lambda k: (k[0], ["short", "kinnex", "ont"].index(k[1]), order[k[2]])):
        print("\t".join(key) + "\t" + "\t".join(cells[key].get(x, ".") for x in arms))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("score")
    p.add_argument("--locus", required=True)
    p.add_argument("--truth", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("sets", nargs="+")
    p = sub.add_parser("merge")
    p.add_argument("--out", required=True)
    p.add_argument("tables", nargs="+")
    p = sub.add_parser("table")
    p.add_argument("tables", nargs="+")
    a = ap.parse_args()
    {"score": score, "merge": merge, "table": table}[a.cmd](a)


if __name__ == "__main__":
    main()
