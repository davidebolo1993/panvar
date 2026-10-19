#!/usr/bin/env python3
"""S2-L runner -- bounded cross-locus validation (frozen manifest; one process at a time; 4 threads; <= 300 s per unit;
1200 s budget counted as time spent in this batch -- development runs included -- plus the declared timeout of the next
unit; a unit that does not fit is NOT_RUN and the matrix is reported incomplete).

Arms (same revision, same inputs, same cell profile = the P1 'opt' development profile with its scientific settings and
every work/consumption cap; sequence_first_v2; S2-A early_incumbent experimental_v1; 6 GiB envelope
footprint_limit_bytes 6,442,450,944; time limit = declared timeout - 20 s; parked options off):
  control   ordinary (unblocked) indexed union scoring; existing origin-cache admission
  combined  + union_fragment_blocks fine_core_4096_experimental + repair_round1_origin_admission reuse_only_experimental
The ONLY declared differences are those two settings. Pairs are compared with s2l_compare.py (v3), categories
operational, profile, route, admission, checkpoints. A control refusal at a work cap is a control resource limitation,
never equivalence and never a reason to raise a cap. A genuine tie is a completed outcome.

Order: gates -> comparator v3 controls (saved outputs only) -> C4, CYP2D6, ACOT at 5x, then C4, CYP2D6 at 30x, then
ACOT 30x (combined first) -> one repeat of the combined ACOT 30x cell if the budget permits -> the summary.

  s2l_run.py manifest B | next B | <unit> | dev <name> <note> -- <command...>
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
OUT = "/Users/davide.bolognini/panvar/tests/results/s2l_cross_locus_20261019"
S2J = "/Users/davide.bolognini/panvar/tests/results/s2j_origin_admission_20261018"
S2H = "/Users/davide.bolognini/panvar/tests/results/s2h_refinement_memory_20261018"
BUDGET = 1200
CELLS = [("c4_5", "c4_5", 90), ("cyp_5", "cyp_5", 90), ("acot_5", "acot_5", 150),
         ("c4_30", "c4_30", 150), ("cyp_30", "cyp_30", 150)]
MANIFEST_B = [("t.genotype_blocked_union_route", "ctest", "genotype_blocked_union_route", 60),
              ("t.genotype_round1_admission", "ctest", "genotype_round1_admission", 30),
              ("compare.v3_controls", "controls", "-", 30)]
for _cell, _tag, _t in CELLS:
    MANIFEST_B += [(f"control.{_tag}", "caller", f"{_cell}:control", _t),
                   (f"combined.{_tag}", "caller", f"{_cell}:combined", _t),
                   (f"compare.{_tag}", "pair", f"control.{_tag}:combined.{_tag}", 30)]
MANIFEST_B += [("combined.acot_30", "caller", "acot_30:combined", 300),
               ("control.acot_30", "caller", "acot_30:control", 300),
               ("compare.acot_30", "pair", "control.acot_30:combined.acot_30", 30),
               ("combined.acot_30.repeat", "caller", "acot_30:combined", 300),
               ("compare.acot_30.repeat", "pair", "combined.acot_30:combined.acot_30.repeat", 30),
               ("summary", "summary", "-", 30)]
GIB6 = 6442450944
MANIFESTS = {"B": MANIFEST_B}


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
    extra = f"{d}/profile.s2l.tsv"
    text = f"profile\ts2l-{arm}\nearly_incumbent\texperimental_v1\nfootprint_limit_bytes\t{GIB6}\n"
    if arm == "combined":
        text += "union_fragment_blocks\tfine_core_4096_experimental\nrepair_round1_origin_admission\treuse_only_experimental\n"
    open(extra, "w").write(text)
    args[args.index("--profile") + 2:args.index("--profile") + 2] = ["--profile", extra]
    return ["nice", f"{BUILD}/panvar", "genotype", *args, "--output-schema", "sequence_first_v2"]


def command(unit, kind, target, timeout, d):
    if kind == "ctest":
        return ["nice", "ctest", "--test-dir", BUILD, "-R", f"^panvar_{target}$", "--output-on-failure", "-j1"], "-"
    if kind == "controls":
        return [sys.executable, os.path.join(HERE, "compare_controls_v3.py"), d, f"{d}/scratch"], "-"
    if kind == "summary":
        return [sys.executable, os.path.join(HERE, "s2l_summary.py"), f"{OUT}/runs", f"{d}/summary.tsv"], "-"
    if kind == "pair":
        x, y = (f"{OUT}/runs/{u}" for u in target.split(":"))
        if not (os.path.exists(f"{x}/call.mosaic_status.tsv") and os.path.exists(f"{y}/call.mosaic_status.tsv")):
            return None, "NOT RUN (condition): an arm did not run"
        sx, sy = tsv(f"{x}/call.mosaic_status.tsv"), tsv(f"{y}/call.mosaic_status.tsv")
        done = ("COMPLETE", "CONDITIONAL_GENOTYPE_EQUIVALENT", "BLOCK_CONSENSUS")
        if sx.get("status") not in done or sy.get("status") not in done:
            return None, (f"NOT RUN (condition): not both completed ({sx.get('status')} / {sy.get('status')}); "
                          "a refusal is a resource limitation, not an equivalence")
        cats = "operational,profile,checkpoints" if target.endswith(".repeat") else "operational,profile,route,admission,checkpoints"
        declared = "" if target.endswith(".repeat") else "union_fragment_blocks,repair_round1_origin_admission"
        return [sys.executable, os.path.join(HERE, "s2l_compare.py"), x, y, f"{d}/compare.tsv",
                "--declared-settings", declared, "--categories", cats], "-"
    cell_id, arm = target.split(":")
    return caller_command(cell_id, arm, timeout, d), f"{arm} arm"


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
    elif arg == "dev":
        sep = sys.argv.index("--")
        sys.exit(dev(sys.argv[2], sys.argv[3], sys.argv[sep + 1:]))
    else:
        sys.exit(run_unit(arg))
