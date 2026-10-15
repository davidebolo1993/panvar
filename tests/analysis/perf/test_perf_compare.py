#!/usr/bin/env python3
"""Negative controls for perf_compare.py on a real pair of runs (copied; the originals are only read).

Usage: test_perf_compare.py <run dir A> <run dir B>   (two scientifically identical runs of one cell)
Requires: the pair passes unaltered and with a timing or cache-counter edit; and each scientific edit to B fails: a
candidate score, a block genotype, a genotype filter, a FASTA base, the call class, a repair-trajectory total and a
non-budget setting.
"""
import os
import shutil
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import perf_compare  # noqa: E402


def edit_table(path, column, change):
    lines = open(path).read().split("\n")
    head = lines[0].split("\t")
    row = lines[1].split("\t")
    row[head.index(column)] = change(row[head.index(column)])
    lines[1] = "\t".join(row)
    open(path, "w").write("\n".join(lines))


def edit_key(path, key, change):
    lines = open(path).read().split("\n")
    for n, line in enumerate(lines):
        p = line.split("\t")
        if p[0] == key:
            p[1] = change(p[1])
            lines[n] = "\t".join(p)
            open(path, "w").write("\n".join(lines))
            return
    raise KeyError(key)


def edit_fasta(path):
    lines = open(path).read().split("\n")
    n = next(k for k, l in enumerate(lines) if l and not l.startswith(">"))
    k = len(lines[n]) // 2
    lines[n] = lines[n][:k] + ("C" if lines[n][k] != "C" else "G") + lines[n][k + 1:]
    open(path, "w").write("\n".join(lines))


def bump(v):
    return v[:-1] + ("1" if v[-1] != "1" else "2")


CASES = [
    ("unaltered", None, True),
    ("timing key", lambda d: edit_key(f"{d}/call.mosaic_status.tsv", "seconds_total", bump), True),
    ("cache counter", lambda d: edit_key(f"{d}/call.mosaic_status.tsv", "dp_evaluations", bump), True),
    ("candidate score", lambda d: edit_table(f"{d}/call.mosaic_scores.tsv", "score", bump), False),
    ("block genotype", lambda d: edit_table(f"{d}/call.block_genotypes.tsv", "allele1", bump), False),
    ("genotype filter", lambda d: edit_table(f"{d}/call.block_genotypes.tsv", "genotype_filter",
                                             lambda v: "LOW_QUALITY" if v == "PASS" else "PASS"), False),
    ("FASTA base", lambda d: edit_fasta(f"{d}/call.mosaic_haplotypes.fa"), False),
    ("call class", lambda d: edit_key(f"{d}/call.mosaic_status.tsv", "call_class", lambda v: "no_call"), False),
    ("repair total", lambda d: edit_table(f"{d}/call.repair_rounds.tsv", "total", bump), False),
    ("non-budget setting", lambda d: edit_key(f"{d}/call.settings.tsv", "kmer_size", lambda v: "29"), False),
]


def main():
    src_a, src_b = sys.argv[1:3]
    failures = 0
    for name, edit, should_pass in CASES:
        tmp = tempfile.mkdtemp(prefix="perf_compare_")
        try:
            da, db = os.path.join(tmp, "a"), os.path.join(tmp, "b")
            shutil.copytree(src_a, da)
            shutil.copytree(src_b, db)
            if edit:
                edit(db)
            ok, problems, _, _ = perf_compare.compare(da, db)
        finally:
            shutil.rmtree(tmp)
        good = ok == should_pass
        failures += not good
        print(("PASS " if good else "FAIL ") + f"{name}: {'identical' if ok else 'different'} "
              f"(expected {'identical' if should_pass else 'different'}) {problems[:1]}")
    print(f"{len(CASES)} cases, {failures} failed")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
