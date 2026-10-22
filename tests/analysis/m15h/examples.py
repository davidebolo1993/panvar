#!/usr/bin/env python3
"""M15-H units: the final help texts (captured), and preserved outputs (help/docs batch: no scientific or assessment change).
Usage: examples.py help_text|preserved_c4 <unit dir>"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "m15d"))
import examples as E  # noqa: E402

R, WT = E.R, E.WT
M15G = f"{R}/m15g_assessment_20261022/runs/a3.ex.complete_c4/c4"
CHECKS = []


def check(ok, what):
    CHECKS.append((ok, what))
    print(("PASS " if ok else "FAIL ") + what)


def ex_help_text(d):
    for flag, name in (("--help", "help.txt"), ("--help-advanced", "help_advanced.txt")):
        p = subprocess.run([E.PANVAR, "genotype", flag], capture_output=True, text=True)
        open(f"{d}/{name}", "w").write(p.stdout)
        check(p.returncode == 0 and len(p.stdout) > 1000, f"{flag}: exit 0, captured to {name} ({len(p.stdout.splitlines())} lines)")
    return 0


def ex_preserved_c4(d):
    """The same command as M15-G's a3.ex.complete_c4: every output identical apart from the run's own name, timing and
    resource rows (the help batch changed no scientific, assessment or settings code path)."""
    out = f"{d}/c4"
    p = subprocess.run(E.cmd("c4", out), capture_output=True, text=True)
    check(p.returncode == 0, "c4: exit 0")
    new = sorted(f[2:] for f in os.listdir(d) if f.startswith("c4."))   # A1: the suffix keeps its leading dot
    old = sorted(f[2:] for f in os.listdir(os.path.dirname(M15G)) if f.startswith("c4."))
    check(new == old, f"c4: the same output files ({len(new)})")
    volatile = ("seconds", "_bytes", "rss_", "footprint", "elapsed", "threads_peak", "_peak",
                # A2: the exact caches' served-work counters, scheduling-dependent at 4 threads (the e2e S2-G allowance)
                "emission_cache_", "union_prefilter_cache_", "_executed_dps", "cached_prefilters_consumed")
    for suffix in new:
        a, b = open(out + suffix).read(), open(M15G + suffix).read()
        a = a.replace(out, "<RUN>").replace("c4_seq", "<RUN>_seq")
        b = b.replace(M15G, "<RUN>").replace("c4_seq", "<RUN>_seq")
        if suffix in (".mosaic_status.tsv", ".dosage_status.tsv", ".repair_memory.tsv", ".repair_rounds.tsv", ".library_estimate.tsv",
                      ".library_sensitivity.tsv"):
            keep = lambda t: [l for l in t.splitlines() if not any(v in l.split("\t")[0] for v in volatile)]
            if suffix in (".repair_memory.tsv", ".repair_rounds.tsv"):
                continue   # timing and memory logs
            check(keep(a) == keep(b), f"c4: {suffix} identical apart from timing and resource rows")
        else:
            check(a == b, f"c4: {suffix} identical")
    return 0


if __name__ == "__main__":
    name, d = sys.argv[1], sys.argv[2]
    os.makedirs(d, exist_ok=True)
    globals()["ex_" + name](d)
    bad = [w for ok, w in CHECKS if not ok]
    print(f"{len(CHECKS) - len(bad)} of {len(CHECKS)} checks passed" + (f"; FAILED: {bad}" if bad else ""))
    sys.exit(1 if bad else 0)
