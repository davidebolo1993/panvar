#!/usr/bin/env python3
"""Scores dev/recruit_seeds.sh: per read set and screen, the fragments certain and recruited (certain or
ambiguous) for the set's locus (the genome set: for any of the five), and the screens' seconds per
technology.

    recruit_seeds_score.py <seeds dir> <sim dir>
"""
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from locityper_score import LOCI, TECHS, reads_path  # noqa: E402
from recruit_score import read_fastq  # noqa: E402

SCREENS = ["syncmer", "minimizer", "all"]


def counts(path, locus):
    certain, recruited = set(), set()
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if locus is not None and r["locus"] != locus:
                continue
            if r["fate"] == "certain":
                certain.add(r["fragment"])
            if r["fate"] in ("certain", "ambiguous"):
                recruited.add(r["fragment"])
    return len(certain), len(recruited)


def main():
    out, sim = sys.argv[1:3]
    cols = ["tech", "locus", "source", "fragments"] + [f"{s}_{x}" for s in SCREENS for x in ("certain", "recruited")]
    print("\t".join(cols))
    for tech in TECHS:
        for locus in LOCI:
            for src in ["locus", "near", "decoy", "boundary"]:
                n = len(list(read_fastq(reads_path(sim, locus, tech, src))))
                row = [tech, locus, src, str(n)]
                for s in SCREENS:
                    path = f"{out}/{s}.{locus}.{src}.{tech}.tsv"
                    if not os.path.exists(path):
                        sys.exit(f"recruit_seeds_score: {path} is missing")
                    row += [str(x) for x in counts(path, locus)]
                print("\t".join(row))
        n = len(list(read_fastq(reads_path(sim, LOCI[0], tech, "new"))))
        row = [tech, "any", "new", str(n)]
        for s in SCREENS:
            row += [str(x) for x in counts(f"{out}/{s}.new.{tech}.tsv", None)]
        print("\t".join(row))
    seconds = {}
    with open(f"{out}/times.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            seconds[(r["screen"], r["tech"])] = seconds.get((r["screen"], r["tech"]), 0) + float(r["seconds"])
    for tech in TECHS:
        print("\t".join(["seconds", tech, "all sets", "."] + [f"{seconds.get((s, tech), 0):.1f}\t." for s in SCREENS]))


if __name__ == "__main__":
    main()
