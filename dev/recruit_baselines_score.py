#!/usr/bin/env python3
"""Scores dev/recruit_baselines.sh: per read set, the fragments each tool takes, beside panvar's certain
and recruited (certain or ambiguous) ones and, given a dev/locityper_compare.sh run, Locityper's. Locus
sets: the locus's fragments taken; near, decoy and the CHM13 genome: fragments from elsewhere taken;
boundary: fragments straddling an edge, and wholly outside, taken. On the genome set BBDuk and
mirabait bait the five loci at once, so their count there is of fragments taken to any of them.

    recruit_baselines_score.py <baselines dir> <sim dir> [locityper_compare dir]
"""
import glob
import gzip
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from locityper_score import LOCI, TECHS, locityper_names, name_of, panvar_sets, reads_path  # noqa: E402
from recruit_score import locus_interval, read_fastq, truth_span  # noqa: E402


def fastq_names(paths):
    names = set()
    for path in paths:
        opener = gzip.open if path.endswith(".gz") else open
        with opener(path, "rt") as f:
            for i, line in enumerate(f):
                if i % 4 == 0:
                    names.add(name_of(line))
    return names


def main():
    out, sim = sys.argv[1:3]
    lt = sys.argv[3] if len(sys.argv) > 3 else None
    mirabait = bool(glob.glob(f"{out}/mirabait.*"))
    cols = ["tech", "locus", "source", "fragments", "bbduk"] + (["mirabait"] if mirabait else []) 
    cols += (["locityper"] if lt else []) + ["panvar_certain", "panvar_recruited"]
    print("\t".join(cols))
    cache = {}
    for tech in TECHS:
        for locus in LOCI:
            for src in ["locus", "near", "decoy", "boundary", "new"]:
                reads = reads_path(sim, locus, tech, src)
                if reads not in cache:
                    cache[reads] = list(read_fastq(reads))
                frags = cache[reads]
                label = f"all.{tech}.new" if src == "new" else f"{locus}.{tech}.{src}"
                sets = [fastq_names([f"{out}/bbduk.{label}.fq.gz"])]
                if mirabait:
                    sets.append(fastq_names(glob.glob(f"{out}/mirabait.{label}*.f*q*")))
                if lt:
                    sets.append(locityper_names(f"{lt}/lt.{locus}.{tech}.{src}.fq"))
                table = (f"{sim}/benchmark/g38.{locus}.{tech}.locus.tsv" if src == "locus" or not lt
                         else f"{lt}/pv.all.{tech}.new.tsv" if src == "new" else f"{lt}/pv.{locus}.{tech}.{src}.tsv")
                if src != "locus" and not lt:
                    table = f"{sim}/benchmark/g38.{locus}.{tech}.{src}.tsv"
                sets.extend(panvar_sets(table, locus) if os.path.exists(table) else (set(), set()))
                if src == "boundary":
                    ls, le = locus_interval(sim, locus)
                    kinds = {}
                    for name, length in frags:
                        x, y = truth_span(tech, name, length)
                        inside = max(0, min(y, le) - max(x, ls))
                        kinds[name] = "in" if inside >= y - x else "out" if inside == 0 else "straddle"
                    for kind in ("straddle", "out"):
                        group = {f for f, k in kinds.items() if k == kind}
                        print("\t".join(str(x) for x in [tech, locus, f"boundary_{kind}", len(group)] + [len(s & group) for s in sets]))
                    continue
                names = {f for f, _ in frags}
                print("\t".join(str(x) for x in [tech, locus, src, len(frags)] + [len(s & names) for s in sets]))


if __name__ == "__main__":
    main()
