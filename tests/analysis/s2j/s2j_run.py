#!/usr/bin/env python3
"""S2-J runner (frozen manifest; one process at a time; <= 4 threads; <= 300 s per unit; 1200 s budget counted as time
spent in this batch -- development runs included -- plus the declared timeout of the next unit; a unit that does not fit
is NOT_RUN).

Order: the gates (round-1 admission, scoring context, caller diagnostics + admission) -> ONE ACOT 30x first-round replay
with repair_round1_origin_admission reuse_only_experimental (S2-H's recorded anchor; the declared policy override;
original caps, 4 threads, 4 GiB guard) -> C4 and ACOT 5x existing vs reuse-only full callers (S2-G integrated arm) and
their explicit-allowance comparison -> CONDITIONAL: if the replay committed round 1 and the gates passed, ONE integrated
ACOT 30x full caller with the policy (unchanged settings; live refinement memory on for round boundaries) -> CONDITIONAL
(Sol, declared before execution): if that caller stops FIRST on the 4 GiB footprint guard, ONE separately recorded
follow-up with a 6 GiB ceiling (footprint_limit_bytes 6,442,450,944), everything else unchanged -- a 6 GiB-envelope run,
never a 4 GiB result. No unit is repeated; no timeout or cap is raised.

  s2j_run.py manifest B | next B | <unit> | dev <name> <note> -- <command...>
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
OUT = "/Users/davide.bolognini/panvar/tests/results/s2j_origin_admission_20261018"
S2H = "/Users/davide.bolognini/panvar/tests/results/s2h_refinement_memory_20261018"
BUDGET = 1200
MANIFEST_B = [("t.genotype_round1_admission", "ctest", "genotype_round1_admission", 30),
              ("t.genotype_scoring_context", "ctest", "genotype_scoring_context", 30),
              ("t.genotype_blocked_union_route", "ctest", "genotype_blocked_union_route", 60),
              ("replay.acot30.round1.reuse", "caller", "acot_30:replay", 300),
              ("caller.c4.existing", "caller", "c4_5:existing", 90),
              ("caller.c4.reuse", "caller", "c4_5:reuse", 90),
              ("compare.c4", "compare", "c4", 30),
              ("caller.acot5.existing", "caller", "acot_5:existing", 150),
              ("caller.acot5.reuse", "caller", "acot_5:reuse", 150),
              ("compare.acot5", "compare", "acot5", 30),
              ("caller.acot30.reuse", "caller", "acot_30:reuse_full", 300),
              ("caller.acot30.reuse_6gib", "caller", "acot_30:reuse_full_6gib", 300)]
GIB4, GIB6 = 4294967296, 6442450944
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
    del args[i:i + 2]
    extra = f"{d}/profile.s2j.tsv"
    text = "profile\ts2j-caller\nunion_fragment_blocks\tfine_core_4096_experimental\n"
    schema = "sequence_first_v2"
    if arm != "existing":
        text += "repair_round1_origin_admission\treuse_only_experimental\n"
    if arm == "replay":   # S2-H's recorded anchor; no early pass, no sequence report
        text += (f"diagnostic_refinement_memory\t1\nearly_incumbent\toff\ndiagnostic_refinement_replay\t"
                 f"{S2H}/runs/caller.acot30.instrumented/call.refinement_replay_inputs.tsv\n")
        schema = "legacy"
    else:
        text += "early_incumbent\texperimental_v1\n"
        if arm.startswith("reuse_full"):
            text += "diagnostic_refinement_memory\t1\n"   # round-boundary memory (measurement only; gated unchanged)
    if arm == "reuse_full_6gib":
        text += f"footprint_limit_bytes\t{GIB6}\n"   # Sol's declared follow-up: the only change
    open(extra, "w").write(text)
    args[args.index("--profile") + 2:args.index("--profile") + 2] = ["--profile", extra]
    return ["nice", f"{BUILD}/panvar", "genotype", *args, "--output-schema", schema]


def status_of(unit):
    return tsv(f"{OUT}/runs/{unit}/call.mosaic_status.tsv")


def command(unit, kind, target, timeout, d):
    if kind == "ctest":
        return ["nice", "ctest", "--test-dir", BUILD, "-R", f"^panvar_{target}$", "--output-on-failure", "-j1"], "-"
    if kind == "compare":
        a, b = f"{OUT}/runs/caller.{target}.existing", f"{OUT}/runs/caller.{target}.reuse"
        if not (os.path.exists(f"{a}/call.mosaic_status.tsv") and os.path.exists(f"{b}/call.mosaic_status.tsv")):
            return None, "NOT RUN (condition): a caller arm did not run"
        return [sys.executable, os.path.join(HERE, "s2j_compare.py"), a, b, f"{d}/compare.tsv"], "-"
    cell_id, arm = target.split(":")
    if arm == "reuse_full":
        replay = tsv(f"{OUT}/runs/replay.acot30.round1.reuse/call.refinement_replay.tsv")
        gates = [r for r in ledger_rows() if r["unit"].startswith("t.")]
        if replay.get("rounds_completed") != "1" or any(r["exit"] != "0" for r in gates) or len(gates) < 3:
            return None, "NOT RUN (condition): the replay did not commit round 1 or a gate failed"
    if arm == "reuse_full_6gib":
        st = status_of("caller.acot30.reuse")
        if not (st.get("status") == "RESOURCE_REFUSAL" and f"exceeded the declared {GIB4}" in st.get("note", "")):
            return None, f"NOT RUN (condition): the 4 GiB caller did not stop first on the footprint guard ({st.get('status', 'absent')})"
        return caller_command(cell_id, arm, timeout, d), "6 GiB-ENVELOPE FOLLOW-UP (declared by Sol before execution): footprint_limit_bytes 6442450944, everything else unchanged"
    return caller_command(cell_id, arm, timeout, d), f"caller {arm}"


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
