#!/usr/bin/env python3
"""S2-L: comparator v3 on saved outputs only (no genotyping rerun).
  1. Recompare S2-K's 6 GiB ACOT 30x run against S2-J's completed run (the S2-K compare_6gib failure), declared setting
     footprint_limit_bytes; categories operational, profile, checkpoints, envelope.
  2. Controls on RELOCATED copies of the 6 GiB run (path-bearing values renamed to the copy's own directory): every
     negative control must fail FOR ITS OWN REASON, every positive control must pass.
  compare_controls_v3.py <out dir> <scratch dir>
"""
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
V3 = os.path.join(HERE, "s2l_compare.py")
S2J = "/Users/davide.bolognini/panvar/tests/results/s2j_origin_admission_20261018/runs/caller.acot30.reuse"
S2K = "/Users/davide.bolognini/panvar/tests/results/s2k_post_refinement_20261019/runs/caller.acot30.6gib"
ARGS = ["--declared-settings", "footprint_limit_bytes", "--categories", "operational,profile,checkpoints,envelope"]


def compare(a, b, out):
    code = subprocess.run([sys.executable, V3, a, b, out, *ARGS], capture_output=True, text=True).returncode
    reason = "-"
    for line in open(out):
        if line.startswith("science_identical\t"):
            reason = line.rstrip("\n").split("\t", 1)[1]
    return code, reason


def lines_of(path):
    return open(path).read().split("\n")


def write_lines(path, lines):
    open(path, "w").write("\n".join(lines))


def set_key(key, value):
    def mutate(path):
        lines = lines_of(path)
        hit = [i for i, l in enumerate(lines) if l.startswith(key + "\t")]
        if not hit:
            raise SystemExit(f"{key} not in {path}")
        cells = lines[hit[0]].split("\t")
        cells[1] = value
        lines[hit[0]] = "\t".join(cells)
        write_lines(path, lines)
    return mutate


def drop_key(key):
    return lambda p: write_lines(p, [l for l in lines_of(p) if not l.startswith(key + "\t")])


def bump(row, column):
    def mutate(path):
        lines = lines_of(path)
        cells = lines[row].split("\t")
        cells[column] = cells[column] + "1"
        lines[row] = "\t".join(cells)
        write_lines(path, lines)
    return mutate


# (label, file, mutation, must_pass, substring the failure reason must contain)
CONTROLS = [
    ("score changed", "call.mosaic_scores.tsv", bump(1, 3), False, "call.mosaic_scores.tsv differs"),
    ("candidate removed", "call.mosaic_candidates.tsv", lambda p: write_lines(p, lines_of(p)[:-2] + [""]), False,
     "call.mosaic_candidates.tsv differs"),
    ("trajectory changed", "call.repair_rounds.tsv", bump(1, 5), False, "trajectory differs"),
    ("filter changed", "call.block_genotypes.tsv", bump(1, -2), False, "call.block_genotypes.tsv differs"),
    ("prediction changed", "call.sequence_locus.tsv", set_key("prediction", "PROVISIONAL"), False, "call.sequence_locus.tsv differs"),
    ("effective scientific setting changed", "call.settings.tsv", set_key("mismatch_band", "7"), False, "mismatch_band"),
    ("resource cap changed (settings)", "call.settings.tsv", set_key("max_physical_origins", "1"), False, "max_physical_origins"),
    ("resource cap changed (status)", "call.mosaic_status.tsv", set_key("max_union_cached_prefilters_consumed", "1"), False,
     "max_union_cached_prefilters_consumed"),
    ("prefilter configuration changed", "call.mosaic_status.tsv", set_key("union_prefilter_cache_budget_bytes", "1"), False,
     "union_prefilter_cache_budget_bytes"),
    ("undeclared setting added", "call.settings.tsv", lambda p: open(p, "a").write("extra_setting\t1\ttest\n"), False, "extra_setting"),
    ("undeclared status value", "call.mosaic_status.tsv", set_key("winner_set_size", "2"), False, "winner_set_size"),
    ("status field missing", "call.mosaic_status.tsv", drop_key("numeric_separation_gap"), False, "numeric_separation_gap (only in A)"),
    ("status field added", "call.mosaic_status.tsv", lambda p: open(p, "a").write("unexpected_new_key\t1\n"), False,
     "unexpected_new_key (only in B)"),
    ("declared: audited prefilter counter", "call.mosaic_status.tsv", set_key("union_prefilter_cache_decisions_held", "1"), True, ""),
    ("declared: timing", "call.mosaic_status.tsv", set_key("seconds_total", "1.0"), True, ""),
]


def main(out, scratch):
    os.makedirs(out, exist_ok=True)
    os.makedirs(scratch, exist_ok=True)
    rows, ok = [], True
    code, reason = compare(S2J, S2K, f"{out}/recompare_acot30_6gib_vs_s2j.tsv")
    rows.append(f"recompare\tacot30_6gib_vs_s2j\t{'PASS' if code == 0 else 'FAIL'}\t{reason}")
    ok &= code == 0
    for label, file, mutate, must_pass, expect in CONTROLS:
        d = os.path.join(scratch, "control_" + "".join(ch if ch.isalnum() else "_" for ch in label))
        if os.path.exists(d):
            shutil.rmtree(d)
        shutil.copytree(S2K, d)
        for entry in os.listdir(d):   # relocation: path-bearing values name the copy's own directory
            p = os.path.join(d, entry)
            text = open(p, errors="surrogateescape").read()
            if S2K in text:
                open(p, "w", errors="surrogateescape").write(text.replace(S2K, d))
        mutate(os.path.join(d, file))
        code, reason = compare(S2J, d, os.path.join(d, "compare_v3.tsv"))
        good = (code == 0) if must_pass else (code != 0 and expect in reason)
        ok &= good
        rows.append(f"control\t{label}\texpected_{'pass' if must_pass else 'fail'}\tobserved_{'pass' if code == 0 else 'fail'}\t"
                    f"{'OK' if good else 'WRONG'}\t{reason}")
    open(f"{out}/controls_v3.tsv", "w").write("\n".join(rows) + "\n")
    print("\n".join(rows))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1], sys.argv[2]))
