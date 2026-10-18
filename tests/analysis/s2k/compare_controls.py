#!/usr/bin/env python3
"""S2-K: recompare S2-J's saved C4 and ACOT 5x outputs with comparator v2, then negative and positive controls on COPIES
of the saved ACOT 5x reuse-only output (scratch directory; the saved results are never modified). Every negative control
must FAIL (a scientific change or an undeclared status change); every positive control must PASS (a declared allowance).
  compare_controls.py <out dir> <scratch dir>
"""
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
S2J = "/Users/davide.bolognini/panvar/tests/results/s2j_origin_admission_20261018/runs"
V2 = os.path.join(HERE, "s2j_compare_v2.py")


def compare(a, b, out):
    return subprocess.run([sys.executable, V2, a, b, out], capture_output=True, text=True).returncode


def edit_line(path, predicate, change):
    lines = open(path).read().split("\n")
    for i, line in enumerate(lines):
        if predicate(i, line):
            lines[i] = change(line)
            break
    else:
        raise SystemExit(f"control could not be applied to {path}")
    open(path, "w").write("\n".join(lines))


def set_status(key, value):
    return lambda path: edit_line(path, lambda i, l: l.startswith(key + "\t"), lambda l: f"{key}\t{value}")


def bump_cell(column):
    def change(line):
        cells = line.split("\t")
        cells[column] = cells[column] + "1" if cells[column] not in ("", "NA") else "0"
        return "\t".join(cells)
    return change


CONTROLS = [  # (name, file, mutation, must_pass)
    ("score changed", "call.mosaic_scores.tsv", lambda p: edit_line(p, lambda i, l: i == 1, bump_cell(3)), False),
    ("candidate removed", "call.mosaic_candidates.tsv",
     lambda p: open(p, "w").write("\n".join(open(p).read().split("\n")[:-2] + [""])), False),
    ("trajectory changed", "call.repair_rounds.tsv", lambda p: edit_line(p, lambda i, l: i == 1, bump_cell(5)), False),
    ("filter changed", "call.block_genotypes.tsv", lambda p: edit_line(p, lambda i, l: i == 1, bump_cell(-2)), False),
    ("prediction changed", "call.sequence_locus.tsv", lambda p: edit_line(p, lambda i, l: l.startswith("prediction\t"),
                                                                           lambda l: "prediction\tPROVISIONAL"), False),
    ("undeclared status value", "call.mosaic_status.tsv", set_status("winner_set_size", "2"), False),
    ("undeclared status key added", "call.mosaic_status.tsv", lambda p: open(p, "a").write("unexpected_new_key\t1\n"), False),
    ("declared: origin cache end bytes", "call.mosaic_status.tsv", set_status("repair_origin_cache_bytes_end", "1"), True),
    ("declared: prefilter misses", "call.mosaic_status.tsv", set_status("union_prefilter_cache_misses", "1"), True),
]


def main(out, scratch):
    os.makedirs(out, exist_ok=True)
    rows, ok = [], True
    for tag in ("c4", "acot5"):
        code = compare(f"{S2J}/caller.{tag}.existing", f"{S2J}/caller.{tag}.reuse", f"{out}/recompare_{tag}.tsv")
        rows.append(f"recompare_{tag}\t{'PASS' if code == 0 else 'FAIL'}")
        ok &= code == 0
    for name, file, mutate, must_pass in CONTROLS:
        d = os.path.join(scratch, "control_" + name.replace(" ", "_").replace(":", ""))
        if os.path.exists(d):
            shutil.rmtree(d)
        shutil.copytree(f"{S2J}/caller.acot5.reuse", d)
        mutate(os.path.join(d, file))
        code = compare(f"{S2J}/caller.acot5.existing", d, os.path.join(d, "compare.tsv"))
        good = (code == 0) == must_pass
        ok &= good
        rows.append(f"control\t{name}\texpected_{'pass' if must_pass else 'fail'}\tobserved_{'pass' if code == 0 else 'fail'}\t"
                    f"{'OK' if good else 'WRONG'}")
    open(f"{out}/controls.tsv", "w").write("\n".join(rows) + "\n")
    print("\n".join(rows))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1], sys.argv[2]))
