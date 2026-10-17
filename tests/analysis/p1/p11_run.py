#!/usr/bin/env python3
"""P1.1: validation closure of the consolidated baseline (frozen manifest; one process at a time; 4 threads).

The C++ implementation under test stays 392fc0b: every unit first checks that src/, include/, CMakeLists.txt and
tests/*.cpp are unchanged since 392fc0b, and the build directory is the P1 build (not rebuilt). A unit starts only if the
time already spent in THIS batch (ledger.tsv) plus its declared timeout fits BUDGET = 1200 s; otherwise it is recorded
NOT_RUN. P1's runs and results are read only (the P1 runs are the saved consolidated/stack/reference/index runs).

Order: validation-tool controls (composer, policy audit of every P1 profile, wrapper controls) -> every enabled CTest
test not yet passed on this implementation, ONE TEST PER UNIT, block_delta and branch_lifetime first -> one more ACOT 5x
reference-route run with a longer declared operational time limit (unit ceiling 300 s -> caller limit 280 s; scientific
settings, work and memory caps unchanged; the first timed-out run r.acot_5 is kept) -> mandatory comparisons with
outcome bookkeeping -> coverage inventory -> optional C4/CYP2D6 30x consolidated-vs-stack pairs.

  p11_run.py manifest | next | <unit>
"""
import csv
import json
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import p1_run  # noqa: E402  (cells, profile composition, arguments, builds, status reader, run_timed)

WT = p1_run.WT
P1 = p1_run.OUT
OUT = "/Users/davide.bolognini/panvar/tests/results/p11_validation_closure_20261017"
BUDGET = 1200
IMPL = "392fc0b"
CT = "ctest"
# enabled tests not yet passed on 392fc0b (P1 logs: 12 passed, seeding_exactness interrupted); declared timeouts from
# earlier cost data (block_delta ~77 s, branch_lifetime ~92 s, report ~97 s) with margin; 60 s otherwise
TESTS = [("panvar_genotype_block_delta", 200), ("panvar_genotype_branch_lifetime", 220),
         ("panvar_genotype_seeding_exactness", 60)] + [(t, 240 if t == "panvar_genotype_report" else 60) for t in [
    "panvar_genotype_headline", "panvar_genotype_cli_contract", "panvar_genotype_curation_contract",
    "panvar_library_sensitivity", "panvar_genotype_block_quality", "panvar_genotype_mosaic_proposal",
    "panvar_genotype_grouped_beam", "panvar_genotype_source_tree_clean", "panvar_genotype_mosaic_binding",
    "panvar_genotype_heldout_mosaic_manifest", "panvar_genotype_mosaic_truth_orientation", "panvar_genotype_heldout_binding",
    "panvar_genotype_lzo_consolidation", "panvar_genotype_subsample_nested", "panvar_genotype_ladder_selection",
    "panvar_genotype_report", "panvar_genotype_loo_development_gates", "panvar_genotype_pair_override",
    "panvar_genotype_production_campaign_plan", "panvar_genotype_downstream_schema", "panvar_genotype_v5_profile",
    "panvar_genotype_v6_profile", "panvar_genotype_repair_pilot_invariance", "panvar_smoke_synthetic",
    "panvar_associate_stats", "panvar_bubble_stats", "panvar_inspect_stats", "panvar_call_stats",
    "panvar_benchmark_stats", "panvar_describe_stats", "panvar_panphorte_stats", "panvar_refine_stats",
    "panvar_smoke_real_c4"]]
MANIFEST = [("composer.tests", "py", "test_compose_profile.py", "", "", 30, ""),
            ("policy.audit", "audit", "", "", "", 30, ""),
            ("wrapper.controls", "py", "test_p1_compare.py", "", "", 120,
             f"{P1}/runs/o.c4_5|{P1}/runs/o.c4_5.r2")]
MANIFEST += [(f"t.{t.replace('panvar_', '')}", "ctest", t, "p1", "", s, "") for t, s in TESTS]
MANIFEST += [("r2.acot_5", "cli", "acot_5", "p1", "ref", 300, ""),
             ("compare", "compare", "", "", "", 60, ""),
             ("coverage", "coverage", "", "", "", 30, ""),
             ("o.c4_30", "cli", "c4_30", "p1", "opt", 150, ""),
             ("x.c4_30", "cli", "c4_30", "stack", "opt", 150, ""),
             ("o.cyp_30", "cli", "cyp_30", "p1", "opt", 150, ""),
             ("x.cyp_30", "cli", "cyp_30", "stack", "opt", 150, ""),
             ("compare.30x", "compare", "", "", "", 60, "")]
# mandatory pairs (a validation-complete verdict needs every one to PASS); run dirs: P1 saved runs or this batch's
MANDATORY = [
    ("stack", "P1:o.c4_5", "P1:x.c4_5"), ("stack", "P1:o.cyp_5", "P1:x.cyp_5"), ("stack", "P1:o.acot_5", "P1:x.acot_5"),
    ("backend", "P1:o.c4_5", "P1:r.c4_5"), ("backend", "P1:o.cyp_5", "P1:r.cyp_5"), ("backend", "P1:o.acot_5", "P11:r2.acot_5"),
    ("production", "P1:r.c4_5", "P1:p.c4_5"), ("production", "P1:r.cyp_5", "P1:p.cyp_5"),
    ("index", "P1:o.c4_5", "P1:l.c4_5.seq:c4_5"), ("index", "P1:o.c4B_5", "P1:l.c4_5.seq:c4B_5"),
    ("index", "P1:o.cyp_5", "P1:l.cyp_5:cyp_5"), ("index", "P1:o.acot_5", "P1:l.acot_5.existing:acot_5"),
]
INFORMATIONAL = [("repeat", "P1:o.c4_5", "P1:o.c4_5.r2"), ("backend_first_attempt", "P1:o.acot_5", "P1:r.acot_5")]
HIGH_COVERAGE = [("stack_30x", "P11:o.c4_30", "P11:x.c4_30"), ("stack_30x", "P11:o.cyp_30", "P11:x.cyp_30")]
LEDGER_HEADER = ("unit\tbuild\trevision\tcpp_unchanged_since_392fc0b\tdeclared_timeout_seconds\twall_seconds\tpeak_rss_gib\t"
                 "exit\ttimed_out\tstatus\tseconds_total\tphys_footprint_peak_bytes\tnote\n")


def spent():
    p = f"{OUT}/ledger.tsv"
    return sum(float(r["wall_seconds"]) for r in csv.DictReader(open(p), delimiter="\t")) if os.path.exists(p) else 0.0


def done():
    p = f"{OUT}/ledger.tsv"
    return {r["unit"] for r in csv.DictReader(open(p), delimiter="\t")} if os.path.exists(p) else set()


def cpp_unchanged():
    return subprocess.run(["git", "-C", WT, "diff", "--quiet", IMPL, "HEAD", "--", "src", "include", "CMakeLists.txt",
                           "tests/*.cpp"]).returncode == 0 and subprocess.run(
        ["git", "-C", WT, "diff", "--quiet", "--", "src", "include", "CMakeLists.txt", "tests/*.cpp"]).returncode == 0


def write_manifest():
    os.makedirs(OUT, exist_ok=True)
    text = "unit\tkind\ttarget\tbuild\tarm\tdeclared_timeout_seconds\textra\n" + "".join(
        "\t".join(str(x) for x in m) + "\n" for m in MANIFEST)
    text += "".join(f"#pair\t{c}\t{k}\t{a}\t{b}\n" for c, ps in (("mandatory", MANDATORY), ("informational", INFORMATIONAL),
                                                                  ("high_coverage", HIGH_COVERAGE)) for k, a, b in ps)
    p = f"{OUT}/MANIFEST.tsv"
    if os.path.exists(p):
        assert open(p).read() == text, "a different manifest is already frozen"
    else:
        open(p, "w").write(text)
    print(text, end="")


def run_dir(spec):
    batch, rest = spec.split(":", 1)
    unit, _, sample = rest.partition(":")
    root = P1 if batch == "P1" else OUT
    return (f"{root}/runs/{unit}/{sample}" if sample else f"{root}/runs/{unit}"), "call"


def compare(pairs_by_class, out_name):
    import p1_compare
    rows = ["class\t" + "\t".join(p1_compare.HEADER)]
    verdicts = {}
    for cls, pairs in pairs_by_class:
        for kind, a, b in pairs:
            da, name = run_dir(a)
            db, _ = run_dir(b)
            r = subprocess.run([sys.executable, os.path.join(HERE, "p1_compare.py"), kind, da, db, name],
                               capture_output=True, text=True)
            line = r.stdout.strip() or "\t".join([kind, a, b, "ERROR", r.stderr.strip()[-300:]])
            verdict = line.split("\t")[3]
            line = "\t".join([kind, a, b] + line.split("\t")[3:])
            rows.append(f"{cls}\t{line}")
            verdicts.setdefault(cls, []).append((kind, a, b, verdict))
    mand = verdicts.get("mandatory", [])
    if mand:
        bad = [f"{k} {a} vs {b}: {v}" for k, a, b, v in mand if v != "PASS"]
        overall = ("VALIDATION_COMPLETE: every mandatory comparison PASSED" if not bad else
                   "VALIDATION_INCOMPLETE: " + "; ".join(bad))
        rows.append(f"#verdict\t{overall}")
    open(f"{OUT}/{out_name}", "w").write("\n".join(rows) + "\n")
    print("\n".join(rows))
    return 0 if all(v == "PASS" for c in verdicts for *_, v in verdicts[c] if c == "mandatory") else 1


def policy_audit():
    import compose_profile as cp
    found = sorted(os.path.join(dp, f) for dp, _, fs in os.walk(f"{P1}/runs") for f in fs
                   if f.startswith("profile.") and f.endswith(".tsv"))
    bad = 0
    for p in found:
        try:
            keys = cp.check_composed(p)
            print(f"ok: {os.path.relpath(p, P1)}: overrides {keys}")
        except ValueError as e:
            bad += 1
            print(f"VIOLATION: {os.path.relpath(p, P1)}: {e}")
    print(f"{len(found)} recorded P1 profiles audited; {bad} violations")
    return 1 if bad or not found else 0


def coverage():
    """Per-test inventory: CTest's list (with DISABLED) against every saved log of this implementation."""
    j = json.loads(subprocess.run([CT, "--test-dir", f"{WT}/build", "--show-only=json-v1"], capture_output=True,
                                  text=True).stdout)
    logs = [(f"P1:{u}", f"{P1}/runs/{u}/run.log") for u in ("t.e2e", "t.focused")]
    logs += [(f"P1.1:{m[0]}", f"{OUT}/runs/{m[0]}/run.log") for m in MANIFEST if m[1] == "ctest"]
    result = {}
    for unit, log in logs:
        if not os.path.exists(log):
            continue
        text = open(log).read()
        for m in re.finditer(r"Start +\d+: (\S+)", text):
            result.setdefault(m.group(1), []).append((unit, "INTERRUPTED"))
        for m in re.finditer(r"Test +#\d+: (\S+) \.+\s*\**\s*(Passed|Failed|Not Run|Timeout|Exception)[^\n]*?([\d.]+) sec", text):
            entries = result.setdefault(m.group(1), [])
            entries[:] = [e for e in entries if not (e[0] == unit and e[1] == "INTERRUPTED")]
            entries.append((unit, f"{m.group(2).upper()} {m.group(3)} s"))
    rows = ["test\tdisabled\tverdict\tevidence"]
    counts = {}
    for t in j["tests"]:
        props = {p["name"]: p["value"] for p in t.get("properties", [])}
        ev = result.get(t["name"], [])
        if props.get("DISABLED"):
            verdict = "DISABLED (preserved)"
        elif any(v.startswith("PASSED") for _, v in ev):
            verdict = "PASSED"
        elif any(v.startswith(("FAILED", "TIMEOUT", "EXCEPTION")) for _, v in ev):
            verdict = "FAILED"
        elif any(v == "INTERRUPTED" for _, v in ev):
            verdict = "INTERRUPTED"
        else:
            verdict = "NEVER_RUN"
        counts[verdict] = counts.get(verdict, 0) + 1
        rows.append(f"{t['name']}\t{'yes' if props.get('DISABLED') else 'no'}\t{verdict}\t"
                    f"{'; '.join(f'{u}: {v}' for u, v in ev) or '-'}")
    enabled = sum(v for k, v in counts.items() if not k.startswith("DISABLED"))
    rows.append(f"#summary\t{json.dumps(counts)}\tenabled={enabled}\t"
                f"{'ENABLED_SUITE_COMPLETE' if counts.get('PASSED', 0) == enabled else 'ENABLED_SUITE_INCOMPLETE'}")
    open(f"{OUT}/coverage.tsv", "w").write("\n".join(rows) + "\n")
    print("\n".join(rows))
    return 0 if counts.get("PASSED", 0) == enabled else 1


def run_unit(unit):
    u, kind, target, build, arm, timeout, extra = next(m for m in MANIFEST if m[0] == unit)
    assert open(f"{OUT}/MANIFEST.tsv").read() and timeout <= 300 and unit not in done()
    os.makedirs(OUT, exist_ok=True)
    unchanged = cpp_unchanged()
    assert unchanged, "the C++ implementation changed since 392fc0b"
    used = spent()
    if used + timeout > BUDGET:
        note = f"NOT RUN: {used:.1f} s spent + declared {timeout} s > {BUDGET} s"
        new = not os.path.exists(f"{OUT}/ledger.tsv")
        with open(f"{OUT}/ledger.tsv", "a") as f:
            if new:
                f.write(LEDGER_HEADER)
            f.write(f"{unit}\t{build or '-'}\t-\t{int(unchanged)}\t{timeout}\t0.00\t0\tnot_run\t0\tNOT_RUN\tna\tna\t{note}\n")
        print(f"{unit}: {note}")
        return 3
    d = f"{OUT}/runs/{unit}"
    assert not os.path.exists(d), f"{d} exists; a run is never repeated into the same place"
    os.makedirs(d)
    me = os.path.abspath(__file__)
    if kind == "py":
        cmd = [sys.executable, os.path.join(HERE, target), *(extra.split("|") if extra else [])]
    elif kind in ("audit", "compare", "coverage"):
        cmd = [sys.executable, me, f"_{kind}", unit]
    elif kind == "ctest":
        cmd = ["nice", CT, "--test-dir", p1_run.BUILDS[build], "-R", f"^{target}$", "--output-on-failure", "-j1"]
    else:
        cell = p1_run.cells()[target]
        prof = p1_run.profile_for(cell, arm, d, timeout, target)
        cmd = ["nice", f"{p1_run.BUILDS[build]}/panvar", "genotype", *p1_run.genotype_args(cell, prof, f"{d}/call")]
    cmd = ["/usr/bin/time", "-l", *cmd]
    rev, _ = p1_run.git_state(build) if build else p1_run.git_state("p1")
    t = time.time()
    with open(f"{d}/run.log", "w") as log:
        log.write(f"$ (declared timeout {timeout} s) " + " ".join(cmd) + "\n")
        log.flush()
        code, timed_out = p1_run.run_timed(cmd, log, timeout)
    wall = time.time() - t
    rss = 0
    for line in open(f"{d}/run.log"):
        if "maximum resident set size" in line:
            rss = int(line.split()[0])
    st = p1_run.status_of(f"{d}/call") if kind == "cli" else {}
    note = "-"
    if kind == "ctest":
        m = re.search(r"Test +#\d+: \S+ \.+\s*\**\s*(\w[\w ]*?)\s+([\d.]+) sec", open(f"{d}/run.log").read())
        note = f"ctest: {m.group(1)} {m.group(2)} s" if m else "ctest: no result line"
    new = not os.path.exists(f"{OUT}/ledger.tsv")
    with open(f"{OUT}/ledger.tsv", "a") as f:
        if new:
            f.write(LEDGER_HEADER)
        f.write(f"{unit}\t{build or '-'}\t{rev}\t{int(unchanged)}\t{timeout}\t{wall:.2f}\t{rss / 2**30:.3f}\t{code}\t"
                f"{int(bool(timed_out))}\t{st.get('status', 'na')}\t{st.get('seconds_total', 'na')}\t"
                f"{st.get('resources_phys_footprint_peak_bytes', 'na')}\t{note}\n")
    print(f"{unit}: exit {code}{' TIMEOUT' if timed_out else ''}, {wall:.1f} s, status {st.get('status', 'na')} {note} "
          f"(spent {spent():.1f}/{BUDGET})")
    return 0


if __name__ == "__main__":
    arg = sys.argv[1]
    if arg == "manifest":
        write_manifest()
    elif arg == "_audit":
        sys.exit(policy_audit())
    elif arg == "_coverage":
        sys.exit(coverage())
    elif arg == "_compare":
        if sys.argv[2] == "compare.30x":
            sys.exit(compare([("high_coverage", HIGH_COVERAGE)], "compare_30x.tsv"))
        sys.exit(compare([("mandatory", MANDATORY), ("informational", INFORMATIONAL)], "compare.tsv"))
    elif arg == "next":
        todo = [m[0] for m in MANIFEST if m[0] not in done()]
        sys.exit(run_unit(todo[0]) if todo else print("all units recorded") or 0)
    else:
        sys.exit(run_unit(arg))
