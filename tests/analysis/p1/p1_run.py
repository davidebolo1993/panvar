#!/usr/bin/env python3
"""P1: validation of the consolidated development baseline (frozen manifest; one process at a time).

MANIFEST (below; written to OUT/MANIFEST.tsv and committed before any execution) lists every unit in execution order
with its DECLARED process-group timeout (<= 300 s). A unit starts only if the experimental time already spent (ledger.tsv,
failures and comparisons included) plus its declared timeout fits BUDGET = 900 s; otherwise it is recorded as not run.
Builds are not experimental units (their times are reported separately).

Builds:  p1      this worktree (the consolidated baseline)
         stack   i4-c-seed-reuse 4cad112 (the reviewed research stack: the optimized reference)
         prod    s1-a-sequence-reporting 5a4285f (src/include identical to genotype-mosaic cee58de: the production caller)
Arms (one composed profile per run, tests/analysis/p1/compose_profile.py):
         opt     the cell's base profile + operational overrides + the full development overlay
         ref     the cell's base profile + operational overrides + the reference route (production caches and
                 consumption caps only; every scoring backend at its default: the I-series "c4v6" arm)
Operational overrides (resource guards only, declared): time_limit_seconds = declared timeout - 20 s (the caller
records TIME_LIMIT before the process-group kill); footprint_limit_bytes = 4 GiB (the local test envelope, not a
production requirement). Scientific settings are each cell's own (perf cells: base profile, --preset short-read,
--dosage-mode repair, --threads 4; CYP2D6 with its panel curation), unchanged from every I-series run.

  p1_run.py manifest | selfcheck | <unit> | next      (`next` runs the next unit not yet in the ledger)
"""
import csv
import os
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "cyp2d6_pilot"))
sys.path.insert(0, os.path.join(HERE, "..", "perf"))
from pilot import run_timed  # noqa: E402
import compose_profile  # noqa: E402

WT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
WORKTREES = "/Users/davide.bolognini/panvar-worktrees"
BUILDS = {"p1": f"{WT}/build", "stack": f"{WORKTREES}/i4-c-seed-reuse/build",
          "prod": f"{WORKTREES}/s1-a-sequence-reporting/build"}
BUILD_REVS = {"p1": None, "stack": "4cad112e240ff925417eac6484bfa31e699ee106",
              "prod": "5a4285f"}
PERF = "/Users/davide.bolognini/panvar/tests/results/genotype_perf_bounded_20261015"
I3D = "/Users/davide.bolognini/panvar/tests/results/i3d_panel_extension_20261016/index"
C4READS = "/Users/davide.bolognini/panvar/tests/results/genotype_c4_production_11fe681/reads"
OUT = "/Users/davide.bolognini/panvar/tests/results/p1_consolidation_20261017"
BUDGET = 900
LEDGER_HEADER = ("unit\tbuild\trevision\tdirty_paths\tdeclared_timeout_seconds\twall_seconds\tpeak_rss_gib\texit\t"
                 "timed_out\tstatus\tseconds_total\tphys_footprint_peak_bytes\tnote\n")
GIB4 = str(4 << 30)


def cells():
    c = {r["id"]: r for r in csv.DictReader(open(f"{PERF}/cells.tsv"), delimiter="\t")}
    b = dict(c["c4_5"])   # a second C4 sample for sequential reuse: same graph, reference, base profile; HG00171 at 5x
    b.update(id="c4B_5", reads1=f"{C4READS}/HG00171/cov5/reads_1.fq", reads2=f"{C4READS}/HG00171/cov5/reads_2.fq",
             library=f"{C4READS}/HG00171/cov5/library.txt")
    c["c4B_5"] = b
    return c


# unit, kind, cell(s), build, arm, declared timeout (s), extra
MANIFEST = [
    ("selfcheck", "selfcheck", "", "", "", 10, ""),
    ("o.c4_5", "cli", "c4_5", "p1", "opt", 45, ""),
    ("x.c4_5", "cli", "c4_5", "stack", "opt", 45, ""),
    ("r.c4_5", "cli", "c4_5", "p1", "ref", 60, ""),
    ("p.c4_5", "cli", "c4_5", "prod", "ref", 60, ""),
    ("o.cyp_5", "cli", "cyp_5", "p1", "opt", 45, ""),
    ("x.cyp_5", "cli", "cyp_5", "stack", "opt", 45, ""),
    ("r.cyp_5", "cli", "cyp_5", "p1", "ref", 60, ""),
    ("p.cyp_5", "cli", "cyp_5", "prod", "ref", 60, ""),
    ("o.acot_5", "cli", "acot_5", "p1", "opt", 90, ""),
    ("x.acot_5", "cli", "acot_5", "stack", "opt", 90, ""),
    ("r.acot_5", "cli", "acot_5", "p1", "ref", 120, ""),
    ("w.c4_5", "write-panel", "c4_5", "p1", "opt", 45, "index/c4_5.panel.pvpl"),
    ("o.c4B_5", "cli", "c4B_5", "p1", "opt", 45, ""),
    ("l.c4_5.seq", "loaded", "c4_5,c4B_5", "p1", "opt", 90, "index/c4_5.panel.pvpl"),
    ("w.cyp_5", "write-panel", "cyp_5", "p1", "opt", 45, "index/cyp_5.panel.pvpl"),
    ("l.cyp_5", "loaded", "cyp_5", "p1", "opt", 45, "index/cyp_5.panel.pvpl"),
    ("l.acot_5.existing", "loaded", "acot_5", "p1", "opt", 90, f"{I3D}/acot_5.panel.pvpl"),
    ("n.c4_index_cyp_args", "refuse", "cyp_5", "p1", "opt", 30, "index/c4_5.panel.pvpl"),
    ("n.cyp_index_uncurated", "refuse", "cyp_5:uncurated", "p1", "opt", 30, "index/cyp_5.panel.pvpl"),
    ("n.c4_index_corrupt", "refuse", "c4_5", "p1", "opt", 30, "index/c4_5.corrupt.pvpl"),
    ("compare", "compare", "", "", "", 60, ""),
    ("t.comparators", "pytests", "", "", "", 30, ""),
    ("t.e2e", "ctest", "^panvar_genotype_command_e2e$", "p1", "", 240, ""),
    ("t.rest", "ctest", "-E:^panvar_genotype_(command_e2e|block_delta|branch_lifetime)$", "p1", "", 300, ""),
    ("t.heavy", "ctest", "^panvar_genotype_(block_delta|branch_lifetime)$", "p1", "", 240, ""),
    # optional 30x sentinels (C4 and CYP2D6): run only if the remaining budget admits them
    ("o.c4_30", "cli", "c4_30", "p1", "opt", 120, ""),
    ("x.c4_30", "cli", "c4_30", "stack", "opt", 120, ""),
    ("o.cyp_30", "cli", "cyp_30", "p1", "opt", 120, ""),
    ("x.cyp_30", "cli", "cyp_30", "stack", "opt", 120, ""),
]
# comparisons, declared before execution: (kind, run A, run B, prefix A, prefix B)
PAIRS = [
    ("stack", "o.c4_5", "x.c4_5"), ("stack", "o.cyp_5", "x.cyp_5"), ("stack", "o.acot_5", "x.acot_5"),
    ("backend", "o.c4_5", "r.c4_5"), ("backend", "o.cyp_5", "r.cyp_5"), ("backend", "o.acot_5", "r.acot_5"),
    ("production", "r.c4_5", "p.c4_5"), ("production", "r.cyp_5", "p.cyp_5"),
    ("index", "o.c4_5", "l.c4_5.seq:c4_5"), ("index", "o.c4B_5", "l.c4_5.seq:c4B_5"),
    ("index", "o.cyp_5", "l.cyp_5:cyp_5"), ("index", "o.acot_5", "l.acot_5.existing:acot_5"),
    ("stack", "o.c4_30", "x.c4_30"), ("stack", "o.cyp_30", "x.cyp_30"),
]


def spent():
    p = f"{OUT}/ledger.tsv"
    return sum(float(r["wall_seconds"]) for r in csv.DictReader(open(p), delimiter="\t")) if os.path.exists(p) else 0.0


def done():
    p = f"{OUT}/ledger.tsv"
    return {r["unit"] for r in csv.DictReader(open(p), delimiter="\t")} if os.path.exists(p) else set()


def write_manifest():
    os.makedirs(OUT, exist_ok=True)
    text = "unit\tkind\tcells\tbuild\tarm\tdeclared_timeout_seconds\textra\n" + "".join(
        "\t".join(str(x) for x in m) + "\n" for m in MANIFEST)
    p = f"{OUT}/MANIFEST.tsv"
    if os.path.exists(p):
        assert open(p).read() == text, "a different manifest is already frozen"
    else:
        open(p, "w").write(text)
    print(text, end="")


def profile_for(cell, arm, d, timeout, tag):
    mode = "full" if arm == "opt" else "reference"
    out = f"{d}/profile.{tag}.tsv"
    return compose_profile.compose(out, f"p1-{arm}-{tag}", cell["base_profile"], mode, [
        f"time_limit_seconds={timeout - 20}:declared unit timeout minus 20 s (records TIME_LIMIT before the kill)",
        f"footprint_limit_bytes={GIB4}:4 GiB local test envelope (not a production requirement)"])


def genotype_args(cell, profile, prefix, curated=True):
    args = ["--preset", "short-read", "--profile", profile, "--threads", "4", "--graph", cell["graph"],
            "--reference-path", cell["reference"], "--reads1", cell["reads1"], "--reads2", cell["reads2"],
            "--library", cell["library"], "--dosage-mode", "repair", "--out", prefix, "--quiet"]
    if cell["curation"] != "none" and curated:
        args[args.index("--reads1"):args.index("--reads1")] = ["--panel-curation", cell["curation"]]
    return args


def git_state(build):
    src = os.path.dirname(BUILDS[build])
    rev = subprocess.run(["git", "-C", src, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    dirty = len([l for l in subprocess.run(["git", "-C", src, "status", "--porcelain", "--untracked-files=no"],
                                           capture_output=True, text=True).stdout.split("\n") if l.strip()])
    return rev, dirty


def status_of(prefix):
    p = prefix + ".mosaic_status.tsv"
    if not os.path.exists(p):
        return {}
    return dict(l.split("\t", 1) for l in open(p).read().split("\n") if "\t" in l)


def selfcheck():
    import perf_bounded
    rows = compose_profile.overlay_rows()
    for k, v in perf_bounded.C4V6_CACHE.items():
        assert rows[k] == str(v), (k, rows[k], v)
    assert set(compose_profile.REFERENCE_ROUTE) == set(perf_bounded.C4V6_CACHE), "reference route = the c4v6 arm"
    for build, rev in BUILD_REVS.items():
        if rev:
            got, dirty = git_state(build)
            assert got.startswith(rev) and dirty == 0, (build, got, dirty)
    print("selfcheck: overlay caches equal perf_bounded.C4V6_CACHE; reference route = c4v6 arm; builds at their revisions")
    return 0


def run_unit(unit):
    u, kind, cellspec, build, arm, timeout, extra = next(m for m in MANIFEST if m[0] == unit)
    assert open(f"{OUT}/MANIFEST.tsv").read(), "freeze the manifest first"
    assert timeout <= 300 and unit not in done()
    used = spent()
    if used + timeout > BUDGET:
        note = f"NOT RUN: {used:.1f} s spent + declared {timeout} s > {BUDGET} s"
        new = not os.path.exists(f"{OUT}/ledger.tsv")
        with open(f"{OUT}/ledger.tsv", "a") as f:
            if new:
                f.write(LEDGER_HEADER)
            f.write(f"{unit}\t{build or '-'}\t-\t-\t{timeout}\t0.00\t0\tnot_run\t0\tNOT_RUN\tna\tna\t{note}\n")
        print(f"{unit}: {note}")
        return 3
    d = f"{OUT}/runs/{unit}"
    assert not os.path.exists(d), f"{d} exists; a run is never repeated into the same place"
    os.makedirs(d)
    c = cells()
    note = ""
    if kind == "selfcheck":
        cmd = [sys.executable, os.path.abspath(__file__), "_selfcheck"]
    elif kind == "compare":
        cmd = [sys.executable, os.path.abspath(__file__), "_compare"]
    elif kind == "pytests":
        cmd = [sys.executable, os.path.abspath(__file__), "_pytests"]
    elif kind == "ctest":
        sel = ["-E", cellspec[3:]] if cellspec.startswith("-E:") else ["-R", cellspec]
        cmd = ["nice", "ctest", "--test-dir", BUILDS[build], *sel, "--output-on-failure", "-j1"]
    elif kind == "cli":
        cell = c[cellspec]
        prof = profile_for(cell, arm, d, timeout, cellspec)
        cmd = ["nice", f"{BUILDS[build]}/panvar", "genotype", *genotype_args(cell, prof, f"{d}/call")]
    else:   # driver units
        index = extra if extra.startswith("/") else f"{OUT}/{extra}"
        files = []
        for spec in cellspec.split(","):
            name, _, flag = spec.partition(":")
            cell = c[name]
            prof = profile_for(cell, arm, d, timeout, name)
            os.makedirs(f"{d}/{name}", exist_ok=True)
            args = genotype_args(cell, prof, f"{d}/{name}/call", curated=(flag != "uncurated"))
            files.append(f"{d}/{name}/args.txt")
            open(files[-1], "w").write("".join(a + "\n" for a in args))
        driver = f"{BUILDS[build]}/genotype_prepared_locus_driver"
        if kind == "write-panel":
            assert not os.path.exists(index), "an index is never overwritten"
            os.makedirs(os.path.dirname(index), exist_ok=True)
            cmd = ["nice", driver, "write-panel", GIB4, "0", f"{d}/driver", index, *files]
        else:
            if unit == "n.c4_index_corrupt":
                src = f"{OUT}/index/c4_5.panel.pvpl"
                assert os.path.exists(src) and not os.path.exists(index)
                data = bytearray(open(src, "rb").read())
                at = len(data) // 2
                data[at] ^= 0x01   # one payload bit flipped in a COPY; the fresh index itself is untouched
                open(index, "wb").write(bytes(data))
                note = f"bit 0 of byte {at} of {len(data)} flipped in a copy of index/c4_5.panel.pvpl"
            assert os.path.exists(index), index
            cmd = ["nice", driver, "loaded", GIB4, "0", f"{d}/driver", index, *files]
    cmd = ["/usr/bin/time", "-l", *cmd]
    rev, dirty = git_state(build) if build else git_state("p1")
    t = time.time()
    with open(f"{d}/run.log", "w") as log:
        log.write(f"$ (declared timeout {timeout} s) " + " ".join(cmd) + "\n")
        log.flush()
        code, timed_out = run_timed(cmd, log, timeout)
    wall = time.time() - t
    rss = 0
    for line in open(f"{d}/run.log"):
        if "maximum resident set size" in line:
            rss = int(line.split()[0])
    prefixes = []
    if kind == "cli":
        prefixes = [f"{d}/call"]
    elif kind in ("loaded", "refuse"):
        prefixes = [f"{d}/{s.partition(':')[0]}/call" for s in cellspec.split(",")]
    st = [status_of(p) for p in prefixes]
    status = ",".join(s.get("status", "none") for s in st) or "na"
    secs = ",".join(s.get("seconds_total", "na") for s in st) or "na"
    peak = ",".join(s.get("resources_phys_footprint_peak_bytes", "na") for s in st) or "na"
    if kind == "refuse":
        refused = code != 0 and not timed_out and all(not s for s in st)
        err = "none"
        if os.path.exists(f"{d}/driver.batch_resource.tsv"):
            err = dict(l.split("\t", 1) for l in open(f"{d}/driver.batch_resource.tsv").read().split("\n")
                       if "\t" in l).get("error", "none")
        note = (note + "; " if note else "") + ("REFUSED before any output" if refused else "NOT REFUSED") + \
            f"; error: {err[:240]}"
    new = not os.path.exists(f"{OUT}/ledger.tsv")
    with open(f"{OUT}/ledger.tsv", "a") as f:
        if new:
            f.write(LEDGER_HEADER)
        f.write(f"{unit}\t{build or '-'}\t{rev}\t{dirty}\t{timeout}\t{wall:.2f}\t{rss / 2**30:.3f}\t{code}\t"
                f"{int(bool(timed_out))}\t{status}\t{secs}\t{peak}\t{note or '-'}\n")
    print(f"{unit}: exit {code}{' TIMEOUT' if timed_out else ''}, {wall:.1f} s, rss {rss / 2**30:.2f} GiB, "
          f"status {status} {note} (spent {spent():.1f}/{BUDGET})")
    return 0


def run_dir(spec):
    unit, _, sample = spec.partition(":")
    return (f"{OUT}/runs/{unit}/{sample}" if sample else f"{OUT}/runs/{unit}"), "call"


def compare_all():
    import p1_compare
    rows = ["\t".join(p1_compare.HEADER)]
    failed = 0
    for kind, a, b in PAIRS:
        da, na = run_dir(a)
        db, nb = run_dir(b)
        if not (os.path.exists(f"{da}/{na}.mosaic_status.tsv") and os.path.exists(f"{db}/{nb}.mosaic_status.tsv")):
            rows.append("\t".join([kind, a, b, "NOT_AVAILABLE"] + ["-"] * 7))
            continue
        r = subprocess.run([sys.executable, os.path.join(HERE, "p1_compare.py"), kind, da, db, na],
                           capture_output=True, text=True)
        rows.append(r.stdout.strip() or "\t".join([kind, a, b, "ERROR", r.stderr.strip()[-200:]]))
        failed += r.returncode != 0
    open(f"{OUT}/compare.tsv", "w").write("\n".join(rows) + "\n")
    print("\n".join(rows))
    failed += leakage()
    return 1 if failed else 0


# outputs that name panel sources or spell candidates (a curated-out template must never appear in them)
LEAK_SCOPE = (".block_genotypes.tsv", ".genotypes.tsv", ".locus_genotype.tsv", ".mosaic_candidates.tsv",
              ".mosaic_scores.tsv", ".equivalence_set.tsv", ".mosaic_haplotypes.fa", ".repair_rounds.tsv")


def leakage():
    """The curated CYP2D6 panel: no template the curation removes may appear in any scientific output of the loaded
    (index) run or of its cold run. Other files that merely LIST the curated names (e.g. a curation record) are
    counted separately, never as leakage."""
    cur = cells()["cyp_5"]["curation"]
    names = [l.split("#")[0].strip() for l in open(cur) if l.split("#")[0].strip()]
    rows = ["run\tcurated_names\tscientific_files_scanned\tleaks\tother_files_naming_curated"]
    bad = 0
    for spec in ("o.cyp_5", "l.cyp_5:cyp_5"):
        d, name = run_dir(spec)
        if not os.path.isdir(d):
            rows.append(f"{spec}\t{len(names)}\t0\tNOT_AVAILABLE\t-")
            continue
        files = sorted(f for f in os.listdir(d) if f.startswith(name + "."))
        sci = [f for f in files if f.endswith(LEAK_SCOPE)]
        leaks = [f"{f}:{n}" for f in sci for n in names if n in open(os.path.join(d, f), errors="replace").read()]
        other = [f for f in files if f not in sci and any(n in open(os.path.join(d, f), errors="replace").read()
                                                            for n in names)]
        bad += bool(leaks)
        rows.append(f"{spec}\t{len(names)}\t{len(sci)}\t{';'.join(leaks) or 'none'}\t{';'.join(other) or 'none'}")
    open(f"{OUT}/leakage.tsv", "w").write("\n".join(rows) + "\n")
    print("\n".join(rows))
    return bad


def pytests():
    failed = 0
    tests = sorted(os.path.join(dp, f) for dp, _, fs in os.walk(os.path.join(HERE, "..")) for f in fs
                   if f.startswith("test_") and f.endswith("_compare.py"))
    for t in tests:
        r = subprocess.run([sys.executable, t], capture_output=True, text=True)
        print(f"{os.path.relpath(t, HERE)}: {'ok' if r.returncode == 0 else 'FAILED'}")
        if r.returncode != 0:
            print(r.stdout[-800:], r.stderr[-800:])
        failed += r.returncode != 0
    return 1 if failed else 0


if __name__ == "__main__":
    arg = sys.argv[1]
    if arg == "manifest":
        write_manifest()
    elif arg == "_selfcheck":
        sys.exit(selfcheck())
    elif arg == "_compare":
        sys.exit(compare_all())
    elif arg == "_pytests":
        sys.exit(pytests())
    elif arg == "next":
        todo = [m[0] for m in MANIFEST if m[0] not in done()]
        sys.exit(run_unit(todo[0]) if todo else print("all units recorded") or 0)
    else:
        sys.exit(run_unit(arg))
