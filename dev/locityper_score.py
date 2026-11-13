#!/usr/bin/env python3
"""Scores dev/locityper_compare.sh: per read set, the fragments Locityper recruits beside panvar's certain
and recruited (certain or ambiguous) ones. Locus sets: the share of the locus's fragments; near, decoy and
the CHM13 genome: fragments from elsewhere taken; boundary: fragments straddling an edge taken, and
fragments wholly outside taken.

    locityper_score.py <compare dir> <sim dir>
"""
import csv
import gzip
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from recruit_score import locus_interval, read_fastq, truth_span  # noqa: E402

LOCI = ["c4", "cyp2d", "gstm", "acot", "lpa"]
TECHS = ["short", "illumina", "ont", "hifi"]


def name_of(header):
    n = header[1:].split()[0]
    return n[:-2] if n[-2:] in ("/1", "/2") else n


def locityper_names(path):
    names = set()
    if not os.path.exists(path):   # Locityper writes a file even when it recruits nothing
        sys.exit(f"locityper_score: {path} is missing: that run did not complete")
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        for i, line in enumerate(f):
            if i % 4 == 0:
                names.add(name_of(line))
    return names


def panvar_sets(path, locus):
    certain, recruited = set(), set()
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["locus"] != locus:
                continue
            if r["fate"] == "certain":
                certain.add(r["fragment"])
            if r["fate"] in ("certain", "ambiguous"):
                recruited.add(r["fragment"])
    return certain, recruited


def reads_path(sim, locus, tech, src):
    pre = f"{sim}/new/{tech}_new" if src == "new" else f"{sim}/{locus}/{tech}_{src}"
    return pre + ("_1.fq.gz" if tech in ("short", "illumina") else ".fq.gz")


def main():
    out, sim = sys.argv[1:3]
    cols = ["tech", "locus", "source", "fragments", "locityper", "panvar_certain", "panvar_recruited"]
    print("\t".join(cols))
    cache = {}
    for tech in TECHS:
        for locus in LOCI:
            for src in ["locus", "near", "decoy", "boundary", "new"]:
                reads = reads_path(sim, locus, tech, src)
                if reads not in cache:
                    cache[reads] = list(read_fastq(reads))
                frags = cache[reads]
                lt = locityper_names(f"{out}/lt.{locus}.{tech}.{src}.fq")
                table = (f"{sim}/benchmark/g38.{locus}.{tech}.locus.tsv" if src == "locus"
                         else f"{out}/pv.all.{tech}.new.tsv" if src == "new" else f"{out}/pv.{locus}.{tech}.{src}.tsv")
                pc, pr = panvar_sets(table, locus)
                n = len(frags)
                if src == "boundary":
                    ls, le = locus_interval(sim, locus)
                    kinds = {}
                    for name, length in frags:
                        x, y = truth_span(tech, name, length)
                        inside = max(0, min(y, le) - max(x, ls))
                        kinds[name] = "in" if inside >= y - x else "out" if inside == 0 else "straddle"
                    for kind in ("straddle", "out"):
                        group = {f for f, k in kinds.items() if k == kind}
                        row = [tech, locus, f"boundary_{kind}", len(group), len(lt & group), len(pc & group), len(pr & group)]
                        print("\t".join(str(x) for x in row))
                    continue
                row = [tech, locus, src, n, len(lt), len(pc), len(pr)]
                print("\t".join(str(x) for x in row))


if __name__ == "__main__":
    main()
