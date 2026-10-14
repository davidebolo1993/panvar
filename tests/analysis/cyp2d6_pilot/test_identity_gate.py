#!/usr/bin/env python3
"""Negative controls for identity_gate.py on a real cell (its outputs are copied; the originals are only read).

Usage: test_identity_gate.py <pilot package> <job> <run dir A> <profile A> <run dir B> <profile B> <declared keys>
Requires: the unaltered pair passes; a timing-only edit passes; and each of these edits to run B fails: a block genotype,
a candidate score, a genotype filter, a FASTA base, a scientific status counter, a repair-trajectory value, and a
declared setting that disagrees with its profile.
"""
import os
import shutil
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import identity_gate  # noqa: E402


def edit_table(path, column, change):
    lines = open(path).read().split("\n")
    head = lines[0].split("\t")
    i = head.index(column)
    row = lines[1].split("\t")
    row[i] = change(row[i])
    lines[1] = "\t".join(row)
    open(path, "w").write("\n".join(lines))


def edit_key(path, key, change):
    lines = open(path).read().split("\n")
    for n, line in enumerate(lines):
        p = line.split("\t")
        if p[0] == key:
            p[1] = change(p[1])
            lines[n] = "\t".join(p)
            break
    else:
        raise KeyError(key)
    open(path, "w").write("\n".join(lines))


def edit_fasta(path):
    lines = open(path).read().split("\n")
    n = next(k for k, l in enumerate(lines) if l and not l.startswith(">"))
    k = len(lines[n]) // 2
    lines[n] = lines[n][:k] + ("C" if lines[n][k] != "C" else "G") + lines[n][k + 1:]
    open(path, "w").write("\n".join(lines))


def bump(v):
    return v[:-1] + ("1" if v[-1] != "1" else "2")


def main():
    pkg, job, run_a, prof_a, run_b, prof_b, declared = sys.argv[1:8]
    declared = set(declared.split(","))
    prof_a, prof_b = os.path.join(pkg, prof_a), os.path.join(pkg, prof_b)
    cases = [
        ("unaltered", None, True),
        ("timing field only (seconds_total)", lambda d: edit_key(f"{d}/{job}.mosaic_status.tsv", "seconds_total", bump), True),
        ("block genotype (allele1)", lambda d: edit_table(f"{d}/{job}.block_genotypes.tsv", "allele1", bump), False),
        ("candidate score", lambda d: edit_table(f"{d}/{job}.mosaic_scores.tsv", "score", bump), False),
        ("genotype filter", lambda d: edit_table(f"{d}/{job}.block_genotypes.tsv", "genotype_filter",
                                                 lambda v: "LOW_QUALITY" if v == "PASS" else "PASS"), False),
        ("FASTA base", lambda d: edit_fasta(f"{d}/{job}.mosaic_haplotypes.fa"), False),
        ("scientific counter (statistical_state_classes)",
         lambda d: edit_key(f"{d}/{job}.mosaic_status.tsv", "statistical_state_classes", bump), False),
        ("repair trajectory (repair_rounds total)", lambda d: edit_table(f"{d}/{job}.repair_rounds.tsv", "total", bump), False),
        ("declared setting disagrees with its profile",
         lambda d: edit_key(f"{d}/{job}.mosaic_status.tsv", "max_state_classes", lambda v: "500000000"), False),
    ]
    failures = 0
    for name, edit, should_pass in cases:
        tmp = tempfile.mkdtemp(prefix="identity_gate_")
        try:
            da, db = os.path.join(tmp, os.path.basename(run_a)), os.path.join(tmp, os.path.basename(run_b))
            os.makedirs(da)
            os.makedirs(db)
            for src, dst in ((os.path.join(pkg, "data", run_a), da), (os.path.join(pkg, "data", run_b), db)):
                for f in os.listdir(src):
                    if f.startswith(job + "."):
                        shutil.copy(os.path.join(src, f), dst)
            if edit:
                edit(db)
            ok, detail = identity_gate.compare(job, da, prof_a, db, prof_b, declared)
        finally:
            shutil.rmtree(tmp)
        good = ok == should_pass
        failures += not good
        print(("PASS " if good else "FAIL ") + f"{name}: gate {'passes' if ok else 'fails'} "
              f"(expected {'pass' if should_pass else 'fail'}) -- {detail[:160]}")
    print(f"{len(cases)} cases, {failures} failed")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
