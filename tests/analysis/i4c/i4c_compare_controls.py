#!/usr/bin/env python3
"""I4-C: negative and positive controls for i4c_compare's typed tie-branch comparison.

Copies the I4-B NA19650 views run (read only) into a temporary directory, applies one change at a time and checks the
qualified verdict: every scientific column, row order, row identity and file presence change must FAIL; a change to an
allowed column alone must PASS (QUALIFIED). The recorded runs are never modified.
Usage: i4c_compare_controls.py <control run dir> <views run dir>
"""
import os
import shutil
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import i4c_compare  # noqa: E402


def mutate(path, column=None, row=1, value=None, op=None):
    lines = open(path).read().rstrip("\n").split("\n")
    head = lines[0].split("\t")
    if op == "delete_row":
        del lines[row]
    elif op == "swap_rows":
        lines[1], lines[2] = lines[2], lines[1]
    elif op == "add_column":
        lines = [l + "\textra" if i == 0 else l + "\t0" for i, l in enumerate(lines)]
    else:
        f = lines[row].split("\t")
        c = head.index(column)
        f[c] = value(f[c]) if callable(value) else value
        lines[row] = "\t".join(f)
    open(path, "w").write("\n".join(lines) + "\n")


def bump(x):
    try:
        return str(int(x) + 1)
    except ValueError:
        try:
            return repr(float(x) + 1.0)
        except ValueError:
            return x + "_changed"


def main(control, views):
    cases_fail = [
        ("tie_branches.tsv", dict(column="branch", value=bump)),
        ("tie_branches.tsv", dict(column="parent_tie", value=bump)),
        ("tie_branches.tsv", dict(column="start_hap1", value=bump)),
        ("tie_branches.tsv", dict(column="start_hap2", value=bump)),
        ("tie_branches.tsv", dict(column="terminal_hap1", value=bump)),
        ("tie_branches.tsv", dict(column="terminal_hap2", value=bump)),
        ("tie_branches.tsv", dict(column="start_total", value=bump)),
        ("tie_branches.tsv", dict(column="terminal_total", value=bump)),
        ("tie_branches.tsv", dict(column="terminal_genotype", value=bump)),
        ("tie_branches.tsv", dict(column="outcome", value="changed_outcome")),
        ("tie_branches.tsv", dict(column="stop", value="changed_stop")),
        ("tie_branches.tsv", dict(column="rounds", value=bump)),
        ("tie_branches.tsv", dict(column="executed_physical_origins", value=bump)),
        ("tie_branches.tsv", dict(op="delete_row", row=2)),
        ("tie_branches.tsv", dict(op="swap_rows")),
        ("tie_branches.tsv", dict(op="add_column")),
        ("tie_branch_memory.tsv", dict(column="event", value="changed_event")),
        ("tie_branch_memory.tsv", dict(column="branch", value=bump)),
        ("tie_branch_memory.tsv", dict(column="cumulative_physical_origins", value=bump)),
        ("tie_branch_memory.tsv", dict(column="slots_live", value=bump)),
        ("tie_branch_memory.tsv", dict(op="delete_row", row=3)),
        ("tie_branch_memory.tsv", dict(op="swap_rows")),
        ("tie_branch_memory.tsv", "remove_file"),
    ]
    cases_pass = [
        ("tie_branches.tsv", dict(column="seconds", value="999.0")),
        ("tie_branch_memory.tsv", dict(column="phys_footprint_peak", value=bump)),
        ("tie_branch_memory.tsv", dict(column="logical_total", value=bump)),
    ]
    failures = []
    with tempfile.TemporaryDirectory() as tmp:
        for expect_ok, cases in ((False, cases_fail), (True, cases_pass)):
            for i, (table, change) in enumerate(cases):
                d = os.path.join(tmp, f"case_{int(expect_ok)}_{i}")
                shutil.copytree(views, d)
                target = os.path.join(d, f"call.{table}")
                if change == "remove_file":
                    os.remove(target)   # one named file in this temporary copy
                else:
                    mutate(target, **change)
                _, qok, problems, _, _, _ = i4c_compare.compare(control, d, "call", "i4b")
                label = f"{table} {change if isinstance(change, str) else change.get('column') or change.get('op')}"
                if qok != expect_ok:
                    failures.append(f"{label}: expected {'pass' if expect_ok else 'fail'}, got {'pass' if qok else 'fail'}")
                else:
                    print(f"  {'pass' if expect_ok else 'fail'} as required: {label}" + ("" if expect_ok else f"  ({problems[0][:80]})"))
    if failures:
        print("CONTROLS FAILED:\n  " + "\n  ".join(failures))
        return 1
    print(f"controls: {len(cases_fail)} changes refused, {len(cases_pass)} allowed-column changes accepted")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1], sys.argv[2]))
