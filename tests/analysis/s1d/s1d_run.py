#!/usr/bin/env python3
"""S1-D validation runner (frozen manifest; one process at a time; <= 4 threads; <= 300 s per unit; 1200 s budget,
counted as time spent in this batch plus the declared timeout of the next unit; a unit that does not fit is NOT_RUN).

Order: the sequence-report unit test (v2 membership, references, FASTA rule, lifecycle) and the CLI contract (v2 dosage
policy) -> caller fixtures (s1d_fixtures.py: panel pair, homozygous, aliases, two off-panel, provisional, one panel + one
off-panel via curation, reverse orientation, legacy identity, loaded index, same-prefix lifecycle) -> the end-to-end
suite (+ test_sequence_first_v2) -> a real extraction example: sequence_first_v2 on the existing C4 5x cell (P1
development overlay; default dosage audit) and extraction of every panel haplotype from the retained C4 graph -> every
other enabled test, heaviest last.

  s1d_run.py manifest | next | <unit>
"""
import csv
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "s1b"))
sys.path.insert(0, os.path.join(HERE, "..", "s1c"))
sys.path.insert(0, os.path.join(HERE, "..", "p1"))
import s1b_run  # noqa: E402  (test lists, test-name rule)
import p1_run  # noqa: E402  (run_timed)

WT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
BUILD = f"{WT}/build"
OUT = "/Users/davide.bolognini/panvar/tests/results/s1d_panel_references_20261018"
BUDGET = 1200
tn = s1b_run.test_name
MANIFEST = [("t.genotype_sequence_report", "ctest", "genotype_sequence_report", 30),
            ("t.genotype_cli_contract", "ctest", "genotype_cli_contract", 30),
            ("fixtures.A", "fixtures", "F3,F9,F1,F7full,F7", 240),
            ("fixtures.B", "fixtures", "F11u,F11c,F12", 240),
            ("t.genotype_command_e2e", "ctest", "genotype_command_e2e", 300),
            ("c4.v2", "c4_v2", "c4_5", 120),
            ("c4.extract", "c4_extract", "c4_5", 60)]
MANIFEST += [(f"t.{tn(t)}", "ctest", tn(t), 60) for t in s1b_run.LIGHT + s1b_run.REST_SMALL + s1b_run.OTHER]
MANIFEST += [(f"t.{t}", "ctest", t, s) for t, s in s1b_run.HEAVY]
# Amendment 1 (after t.genotype_command_e2e failed on a test defect, before further execution; it stays recorded): the
# S1-D e2e test compared the legacy block table byte for byte between audit and off, but that table's dosage_filter /
# dosage_filter_reasons columns record the dosage mode by design (off: DOSAGE_NOT_EVALUATED:dosage_mode_off). The test
# now requires every other column identical and the off record to say not evaluated; rerun as .a1.
AMENDMENT_1 = [("t.genotype_command_e2e.a1", "ctest", "genotype_command_e2e", 300)]
ORDER = [m[0] for m in MANIFEST[:5]] + [m[0] for m in AMENDMENT_1] + [m[0] for m in MANIFEST[5:]]
LEDGER = "unit\trevision\tdirty_tree\tdeclared_timeout_seconds\twall_seconds\tpeak_rss_gib\texit\ttimed_out\tnote\n"


def spent():
    p = f"{OUT}/ledger.tsv"
    return sum(float(r["wall_seconds"]) for r in csv.DictReader(open(p), delimiter="\t")) if os.path.exists(p) else 0.0


def done():
    p = f"{OUT}/ledger.tsv"
    return {r["unit"] for r in csv.DictReader(open(p), delimiter="\t")} if os.path.exists(p) else set()


def state():
    rev = subprocess.run(["git", "-C", WT, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "-C", WT, "status", "--porcelain"], capture_output=True, text=True).stdout.strip()
    return rev, len([l for l in dirty.split("\n") if l.strip()])


def write_manifest():
    os.makedirs(OUT, exist_ok=True)
    text = "unit\tkind\ttarget\tdeclared_timeout_seconds\n" + "".join("\t".join(str(x) for x in m) + "\n" for m in MANIFEST)
    p = f"{OUT}/MANIFEST.tsv"
    if os.path.exists(p):
        assert open(p).read() == text, "a different manifest is already frozen"
    else:
        open(p, "w").write(text)
    print(text, end="")


def record(unit, rev, dirty, timeout, wall, rss, code, timed_out, note):
    new = not os.path.exists(f"{OUT}/ledger.tsv")
    with open(f"{OUT}/ledger.tsv", "a") as f:
        if new:
            f.write(LEDGER)
        f.write(f"{unit}\t{rev}\t{dirty}\t{timeout}\t{wall:.2f}\t{rss:.3f}\t{code}\t{int(bool(timed_out))}\t{note}\n")


def run_unit(unit):
    u, kind, target, timeout = next(m for m in MANIFEST + AMENDMENT_1 if m[0] == unit)
    if any(m[0] == unit for m in AMENDMENT_1):
        assert os.path.exists(f"{OUT}/AMENDMENT_1.tsv"), "freeze amendment 1 first"
    assert open(f"{OUT}/MANIFEST.tsv").read() and timeout <= 300 and unit not in done()
    rev, dirty = state()
    assert dirty == 0, "the tree must be committed and clean"
    used = spent()
    if used + timeout > BUDGET:
        record(unit, rev, dirty, timeout, 0.0, 0.0, "not_run", 0, f"NOT RUN: {used:.1f} s + {timeout} s > {BUDGET} s")
        print(f"{unit}: NOT RUN")
        return 3
    d = f"{OUT}/runs/{unit}"
    assert not os.path.exists(d), f"{d} exists; a run is never repeated into the same place"
    os.makedirs(d)
    me = os.path.abspath(__file__)
    if kind == "ctest":
        cmd = ["nice", "ctest", "--test-dir", BUILD, "-R", f"^panvar_{target}$", "--output-on-failure", "-j1"]
    elif kind == "fixtures":
        cmd = ["nice", sys.executable, os.path.join(HERE, "s1d_fixtures.py"), f"{OUT}/fixtures", f"{BUILD}/panvar",
               f"{BUILD}/genotype_prepared_locus_driver", target]
    elif kind == "c4_v2":
        cell = p1_run.cells()[target]
        prof = p1_run.profile_for(cell, "opt", d, timeout, target)
        args = p1_run.genotype_args(cell, prof, f"{d}/call")
        i = args.index("--dosage-mode")
        del args[i:i + 2]   # sequence_first_v2 decides dosage: audit by default
        cmd = ["nice", f"{BUILD}/panvar", "genotype", *args, "--output-schema", "sequence_first_v2"]
    else:
        cmd = [sys.executable, os.path.join(HERE, "extract_panel_haplotype.py"), "--check", p1_run.cells()[target]["graph"],
               f"{OUT}/runs/c4.v2/call.sequence_haplotypes.tsv"]
    cmd = ["/usr/bin/time", "-l", *cmd]
    t = time.time()
    with open(f"{d}/run.log", "w") as log:
        log.write(f"$ (declared timeout {timeout} s) " + " ".join(cmd) + "\n")
        log.flush()
        code, timed_out = p1_run.run_timed(cmd, log, timeout)
    wall = time.time() - t
    rss = 0
    text = open(f"{d}/run.log").read()
    for line in text.splitlines():
        if "maximum resident set size" in line:
            rss = int(line.split()[0])
    note = "-"
    if kind == "ctest":
        m = re.search(r"Test +#\d+: \S+ \.+\s*\**\s*(\w[\w ]*?)\s+([\d.]+) sec", text)
        note = f"ctest: {m.group(1)} {m.group(2)} s" if m else "ctest: no result line"
    record(unit, rev, dirty, timeout, wall, rss / 2**30, code, timed_out, note)
    print(f"{unit}: exit {code}{' TIMEOUT' if timed_out else ''}, {wall:.1f} s, {note} (spent {spent():.1f}/{BUDGET})")
    return 0


if __name__ == "__main__":
    arg = sys.argv[1]
    if arg == "manifest":
        write_manifest()
    elif arg == "amendment1":
        text = "".join("\t".join(str(x) for x in m) + "\n" for m in AMENDMENT_1)
        p = f"{OUT}/AMENDMENT_1.tsv"
        assert not os.path.exists(p) or open(p).read() == text
        open(p, "w").write(text)
        print(text, end="")
    elif arg == "next":
        todo = [u for u in ORDER if u not in done()]
        sys.exit(run_unit(todo[0]) if todo else print("all units recorded") or 0)
    else:
        sys.exit(run_unit(arg))
