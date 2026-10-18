#!/usr/bin/env python3
"""S2-G runner (frozen manifests; one process at a time; <= 4 threads; <= 300 s per unit; 1200 s budget counted as time
spent in this batch -- development runs included -- plus the declared timeout of the next unit; a unit that does not fit
is NOT_RUN).

Part A (frozen before any S2-G execution): one complete ACOT 30x reference -- ordinary blocked (S2-E route) with the
two declared validation-only resource extensions, physical-origin cap 1,000,000,000 and executed state-class cap
1,000,000,000, every other setting unchanged -- then a bit-for-bit comparison of all 108,939 candidates' score
components and the numerical top set against the S2-F combined run (blocked_fine.acot30 at 3ab2a11; the kernel source is
unchanged since). If another limit binds or the comparison fails: stop.
Part B (frozen after the integration gates were developed, before its execution): the gates (caller blocked-route gate,
the S2-F kernel gates, the full e2e) -> C4 and CYP2D6 5x full-caller control/integrated runs and their explicit-allowance
comparison (s2g_compare.py) -> ONE integrated ACOT 30x full-caller attempt with the original main-arm caps and unchanged
refinement (prioritised over ACOT 5x) -> ACOT 5x control/integrated/compare if budget permits. Caller arms: the P1 'opt'
profile (indexed union backend, fine views, origin reuse; time limit = declared timeout - 20 s; 4 GiB footprint guard),
sequence_first_v2, early_incumbent experimental_v1; integrated adds union_fragment_blocks fine_core_4096_experimental.

  s2g_run.py manifest A|B | next A|B | <unit> | dev <name> <note> -- <command...>
"""
import csv
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "p1"))
import p1_run  # noqa: E402

WT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
BUILD = f"{WT}/build"
OUT = "/Users/davide.bolognini/panvar/tests/results/s2g_blocked_integration_20261018"
S2C = "/Users/davide.bolognini/panvar/tests/results/s2c_demand_driven_20261018/runs"
S2F = "/Users/davide.bolognini/panvar/tests/results/s2f_blocked_fine_aggregates_20261018/runs"
BUDGET = 1200
EXTENDED_ORIGIN_CAP = 1_000_000_000        # declared by Sol (S2-G 1) before execution; validation only
EXTENDED_STATE_CLASS_CAP = 1_000_000_000   # declared by Sol (S2-G 1) before execution; validation only
MANIFEST_A = [("ref.acot30.ordinary_blocked_ext", "reference", "acot_30", 300),
              ("ref.acot30.compare_vs_combined", "ref_compare", "acot_30", 30)]
MANIFEST_B = [("t.genotype_blocked_union_route", "ctest", "genotype_blocked_union_route", 60),
              ("t.genotype_blocked_fine", "ctest", "genotype_blocked_fine", 120),
              ("t.genotype_blocked_union", "ctest", "genotype_blocked_union", 60),
              ("t.genotype_command_e2e", "ctest", "genotype_command_e2e", 300)]
for _cell, _tag, _t in (("c4_5", "c4", 90), ("cyp_5", "cyp", 90)):
    MANIFEST_B += [(f"caller.{_tag}.control", "caller", f"{_cell}:control", _t),
                   (f"caller.{_tag}.integrated", "caller", f"{_cell}:integrated", _t),
                   (f"compare.{_tag}", "compare", _tag, 30)]
MANIFEST_B += [("caller.acot30.integrated", "caller", "acot_30:integrated", 300),
               ("caller.acot5.control", "caller", "acot_5:control", 150),
               ("caller.acot5.integrated", "caller", "acot_5:integrated", 150),
               ("compare.acot5", "compare", "acot5", 30)]
TAGS = {"c4_5": "c4", "cyp_5": "cyp", "acot_5": "acot5", "acot_30": "acot30"}
MANIFESTS = {"A": MANIFEST_A, "B": MANIFEST_B}


def ledger_rows():
    p = f"{OUT}/ledger.tsv"
    return list(csv.DictReader(open(p), delimiter="\t")) if os.path.exists(p) else []


def spent():
    return sum(float(r["wall_seconds"]) for r in ledger_rows())


def done():
    return {r["unit"] for r in ledger_rows()}


def state():
    rev = subprocess.run(["git", "-C", WT, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "-C", WT, "status", "--porcelain"], capture_output=True, text=True).stdout.strip()
    return rev, len([line for line in dirty.split("\n") if line.strip()])


def write_manifest(part):
    os.makedirs(OUT, exist_ok=True)
    rows = MANIFESTS[part]
    assert rows, f"manifest {part} is empty"
    text = "unit\tkind\ttarget\tdeclared_timeout_seconds\n" + "".join("\t".join(str(x) for x in m) + "\n" for m in rows)
    p = f"{OUT}/MANIFEST_{part}.tsv"
    if os.path.exists(p):
        assert open(p).read() == text, "a different manifest is already frozen"
    else:
        open(p, "w").write(text)
    print(text, end="")


def record(unit, rev, dirty, timeout, wall, rss, code, timed_out, note):
    os.makedirs(OUT, exist_ok=True)
    p = f"{OUT}/ledger.tsv"
    if not os.path.exists(p):
        open(p, "w").write("unit\trevision\tdirty_tree\tdeclared_timeout_seconds\twall_seconds\tpeak_rss_gib\texit\ttimed_out\tnote\n")
    with open(p, "a") as f:
        f.write(f"{unit}\t{rev}\t{dirty}\t{timeout}\t{wall:.2f}\t{rss:.3f}\t{code}\t{int(bool(timed_out))}\t{note}\n")


def timed(cmd, d, timeout):
    cmd = ["/usr/bin/time", "-l", *cmd]
    t = time.time()
    with open(f"{d}/run.log", "w") as log:
        log.write(f"$ (declared timeout {timeout} s) " + " ".join(cmd) + "\n")
        log.flush()
        code, timed_out = p1_run.run_timed(cmd, log, timeout)
    wall = time.time() - t
    text = open(f"{d}/run.log").read()
    rss = 0
    for line in text.splitlines():
        if "maximum resident set size" in line:
            rss = int(line.split()[0])
    return code, timed_out, wall, rss, text


def tsv(path):
    out = {}
    if os.path.exists(path):
        for line in open(path, errors="replace"):
            k, _, v = line.rstrip("\n").partition("\t")
            out[k] = v
    return out


def caller_command(cell_id, arm, timeout, d):
    cell = p1_run.cells()[cell_id]
    prof = p1_run.profile_for(cell, "opt", d, timeout, cell_id)
    args = p1_run.genotype_args(cell, prof, f"{d}/call")
    i = args.index("--dosage-mode")
    del args[i:i + 2]   # sequence_first_v2: dosage is audit-only by contract
    extra = f"{d}/profile.s2g.tsv"
    text = "profile\ts2g-caller\nearly_incumbent\texperimental_v1\n"
    if arm == "integrated":
        text += "union_fragment_blocks\tfine_core_4096_experimental\n"
    open(extra, "w").write(text)
    args[args.index("--profile") + 2:args.index("--profile") + 2] = ["--profile", extra]
    return ["nice", f"{BUILD}/panvar", "genotype", *args, "--output-schema", "sequence_first_v2"]


def command(unit, kind, target, timeout, d):
    harness = f"{BUILD}/s2c_progressive_harness"
    if kind == "ctest":
        return ["nice", "ctest", "--test-dir", BUILD, "-R", f"^panvar_{target}$", "--output-on-failure", "-j1"], "-"
    if kind == "caller":
        cell_id, arm = target.split(":")
        return caller_command(cell_id, arm, timeout, d), f"caller {arm}"
    if kind == "compare":
        a, b = f"{OUT}/runs/caller.{target}.control", f"{OUT}/runs/caller.{target}.integrated"
        if not (os.path.exists(f"{a}/call.mosaic_status.tsv") and os.path.exists(f"{b}/call.mosaic_status.tsv")):
            return None, "NOT RUN (condition): a caller arm did not run"
        return [sys.executable, os.path.join(HERE, "s2g_compare.py"), a, b, f"{d}/compare.tsv"], "-"
    if kind == "reference":
        cell = p1_run.cells()[target]
        inputs = f"{S2C}/dump.acot30/call.union_kernel_inputs.bin"
        return ["nice", harness, "blocked", inputs, cell["reads1"], cell["reads2"], f"{d}/blocked", "4096", str(4 << 30),
                str(EXTENDED_ORIGIN_CAP), str(EXTENDED_STATE_CLASS_CAP)], \
            (f"VALIDATION-ONLY REFERENCE: physical-origin cap {EXTENDED_ORIGIN_CAP}, state-class cap {EXTENDED_STATE_CLASS_CAP} "
             "(declared before execution)")
    if kind == "ref_compare":
        ref = f"{OUT}/runs/ref.acot30.ordinary_blocked_ext"
        if tsv(f"{ref}/blocked.tsv").get("completed") != "yes":
            return None, "NOT RUN (condition): the reference did not complete"
        inputs = f"{S2C}/dump.acot30/call.union_kernel_inputs.bin"
        return [harness, "compare", inputs, f"{ref}/blocked", f"{S2F}/blocked_fine.acot30/blocked_fine", f"{d}/compare.tsv"], "-"
    raise ValueError(kind)


def run_unit(unit):
    part, (u, kind, target, timeout) = next((p, m) for p, rows in MANIFESTS.items() for m in rows if m[0] == unit)
    assert os.path.exists(f"{OUT}/MANIFEST_{part}.tsv") and timeout <= 300 and unit not in done()
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
    cmd, note = command(unit, kind, target, timeout, d)
    if cmd is None:
        record(unit, rev, dirty, timeout, 0.0, 0.0, "not_run", 0, note)
        print(f"{unit}: {note}")
        return 3
    code, timed_out, wall, rss, text = timed(cmd, d, timeout)
    if kind == "ctest":
        m = re.search(r"Test +#\d+: \S+ \.+\s*\**\s*(\w[\w ]*?)\s+([\d.]+) sec", text)
        note = f"ctest: {m.group(1)} {m.group(2)} s" if m else "ctest: no result line"
    record(unit, rev, dirty, timeout, wall, rss / 2**30, code, timed_out, note)
    print(f"{unit}: exit {code}{' TIMEOUT' if timed_out else ''}, {wall:.1f} s, {note} (spent {spent():.1f}/{BUDGET})")
    return 0


def dev(name, note, cmd):
    """A development run: counted against the budget, recorded, never a result."""
    used = spent()
    assert used + 300 <= BUDGET, f"development run refused: {used:.1f} s + 300 s > {BUDGET} s"
    d = f"{OUT}/dev/{name}"
    assert not os.path.exists(d), f"{d} exists"
    os.makedirs(d)
    code, timed_out, wall, rss, text = timed(cmd, d, 300)
    record(name, "development", "-", "-", wall, rss / 2**30, code, timed_out, f"development run (counted): {note}")
    print(text[-4000:])
    print(f"{name}: exit {code}, {wall:.1f} s (spent {spent():.1f}/{BUDGET})")
    return code


if __name__ == "__main__":
    arg = sys.argv[1]
    if arg == "manifest":
        write_manifest(sys.argv[2])
    elif arg == "next":
        rows = MANIFESTS[sys.argv[2]]
        todo = [m[0] for m in rows if m[0] not in done()]
        sys.exit(run_unit(todo[0]) if todo else print("all units recorded") or 0)
    elif arg == "devcaller":   # devcaller <name> <cell> <arm> <timeout>: a counted development caller run
        name, cell_id, arm, timeout = sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5])
        d = f"{OUT}/dev/{name}"
        os.makedirs(d)
        sys.exit(dev(name + ".run", f"caller {cell_id} {arm} (calibration)", caller_command(cell_id, arm, timeout, d)))
    elif arg == "devcompare":
        sys.exit(dev(sys.argv[2], "comparator calibration", [sys.executable, os.path.join(HERE, "s2g_compare.py"),
                                                              f"{OUT}/dev/{sys.argv[3]}", f"{OUT}/dev/{sys.argv[4]}",
                                                              f"{OUT}/dev/{sys.argv[2]}.tsv"]))
    elif arg == "dev":
        sep = sys.argv.index("--")
        sys.exit(dev(sys.argv[2], sys.argv[3], sys.argv[sep + 1:]))
    else:
        sys.exit(run_unit(arg))
