#!/usr/bin/env python3
"""P1.1: negative and positive controls for the P1 wrapper (p1_compare.py) on a real pair of completed runs.

The pair is copied; the originals are only read. Each case edits run B's copy and runs p1_compare.py on (A, B'):
  must FAIL (exit 1): a candidate score, a score component, a block genotype, a block filter, phase support, a FASTA
      base, a repair-trajectory total, a dosage table field, and a completed call turned into a scientific no-call;
  must be UNRESOLVED (exit 2): a completion status turned into an interrupted one (TIME_LIMIT, RESOURCE_REFUSAL, a cap
      stop) on either side, and a missing status;
  may PASS (exit 0): unaltered; declared timing (seconds_total, a trajectory's seconds column); declared provenance
      (prepared_locus_* rows added); permitted executed work (indexed_union_executed_dps).
  test_p1_compare.py <run dir A> <run dir B>    (two scientifically identical completed runs, prefix "call")
"""
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
failures = []


def edit_table(path, column, change, row=1):
    lines = open(path).read().split("\n")
    head = lines[0].split("\t")
    cells = lines[row].split("\t")
    cells[head.index(column)] = change(cells[head.index(column)])
    lines[row] = "\t".join(cells)
    open(path, "w").write("\n".join(lines))


def edit_status(path, changes):
    out, seen = [], set()
    for line in (l for l in open(path).read().split("\n") if l != ""):
        key = line.split("\t", 1)[0]
        if key in changes:
            out.append(key + "\t" + changes[key])
            seen.add(key)
        else:
            out.append(line)
    out += [k + "\t" + v for k, v in changes.items() if k not in seen]
    open(path, "w").write("\n".join(out) + "\n")


def flip_fasta(path):
    lines = open(path).read().split("\n")
    seq = lines[1]   # the FASTA is wrapped (60 columns): flip a base inside the first sequence line
    at = len(seq) // 2
    lines[1] = seq[:at] + ("C" if seq[at] != "C" else "G") + seq[at + 1:]
    open(path, "w").write("\n".join(lines))


def run_case(name, src_a, src_b, edit, expect, on_a=False):
    tmp = tempfile.mkdtemp(prefix="p11_wrap_")
    a, b = os.path.join(tmp, "A"), os.path.join(tmp, "B")
    shutil.copytree(src_a, a)
    shutil.copytree(src_b, b)
    if edit:
        edit(a if on_a else b)
    r = subprocess.run([sys.executable, os.path.join(HERE, "p1_compare.py"), "control", a, b, "call"],
                       capture_output=True, text=True)
    got = {0: "PASS", 1: "FAIL", 2: "UNRESOLVED"}.get(r.returncode, f"exit {r.returncode}")
    ok = got == expect
    print(f"{'ok' if ok else 'FAILED'}: {name}: {got} (expected {expect})")
    if not ok:
        failures.append(name)
        print(r.stdout[:600], r.stderr[-600:])
    shutil.rmtree(tmp)


def main(src_a, src_b):
    T = lambda f: lambda d: os.path.join(d, "call." + f)  # noqa: E731
    bump = lambda x: repr(float(x) + 0.5)  # noqa: E731
    cases = [
        ("unaltered pair", None, "PASS"),
        ("timing: seconds_total", lambda d: edit_status(T("mosaic_status.tsv")(d), {"seconds_total": "999.0"}), "PASS"),
        ("timing: trajectory seconds", lambda d: edit_table(T("repair_rounds.tsv")(d), "seconds", lambda x: "9.99"), "PASS"),
        ("provenance: prepared_locus rows added",
         lambda d: edit_status(T("mosaic_status.tsv")(d), {"prepared_locus": "1", "prepared_locus_origin": "loaded",
                                                          "prepared_locus_panel_source": "loaded"}), "PASS"),
        ("executed work: indexed_union_executed_dps",
         lambda d: edit_status(T("mosaic_status.tsv")(d), {"indexed_union_executed_dps": "1"}), "PASS"),
        ("candidate score", lambda d: edit_table(T("mosaic_scores.tsv")(d), "score", bump), "FAIL"),
        ("score component (fragment term)", lambda d: edit_table(T("mosaic_scores.tsv")(d), "fragment_term", bump), "FAIL"),
        ("block genotype", lambda d: edit_table(T("block_genotypes.tsv")(d), "allele2", lambda x: x + "_edited"), "FAIL"),
        ("block filter", lambda d: edit_table(T("block_genotypes.tsv")(d), "filter", lambda x: "EDITED"), "FAIL"),
        ("phase support", lambda d: edit_table(T("block_genotypes.tsv")(d), "phase_support", lambda x: "EDITED"), "FAIL"),
        ("FASTA base", lambda d: flip_fasta(T("mosaic_haplotypes.fa")(d)), "FAIL"),
        ("repair trajectory total", lambda d: edit_table(T("repair_rounds.tsv")(d), "total", bump), "FAIL"),
        ("call turned into a scientific no-call",
         lambda d: edit_status(T("mosaic_status.tsv")(d), {"status": "RESCORE_INCOMPLETE", "call_class": "no_call",
                                                          "no_call_reason": "distinct_tie"}), "FAIL"),
        ("completion: TIME_LIMIT", lambda d: edit_status(T("mosaic_status.tsv")(d), {"status": "TIME_LIMIT"}), "UNRESOLVED"),
        ("completion: RESOURCE_REFUSAL", lambda d: edit_status(T("mosaic_status.tsv")(d), {"status": "RESOURCE_REFUSAL"}),
         "UNRESOLVED"),
        ("completion: cap stop (RESCORE_INCOMPLETE, no scientific reason)",
         lambda d: edit_status(T("mosaic_status.tsv")(d), {"status": "RESCORE_INCOMPLETE", "no_call_reason": "none"}),
         "UNRESOLVED"),
        ("completion: missing status", lambda d: os.remove(T("mosaic_status.tsv")(d)), "UNRESOLVED"),
    ]
    if os.path.exists(os.path.join(src_b, "call.dosage_blocks.tsv")):
        head = open(os.path.join(src_b, "call.dosage_blocks.tsv")).readline().rstrip("\n").split("\t")
        col = "dosage_status" if "dosage_status" in head else head[-1]
        cases.append((f"dosage table field ({col})",
                      lambda d: edit_table(T("dosage_blocks.tsv")(d), col, lambda x: x + "_edited"), "FAIL"))
    for name, edit, expect in cases:
        run_case(name, src_a, src_b, edit, expect)
    run_case("completion: TIME_LIMIT on side A", src_a, src_b,
             lambda d: edit_status(T("mosaic_status.tsv")(d), {"status": "TIME_LIMIT"}), "UNRESOLVED", on_a=True)
    print("ALL PASSED" if not failures else "FAILURES: " + ", ".join(failures))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1], sys.argv[2]))
