#!/usr/bin/env python3
"""S1-B validation runner (frozen manifest; one process at a time; <= 4 threads; <= 300 s per unit; 1200 s budget,
counted as time spent in this batch plus the declared timeout of the next unit; a unit that does not fit is recorded
NOT_RUN). The implementation under test is this worktree's build; every unit records the revision and requires a clean
tree for src/, include/, CMakeLists.txt and tests/*.cpp.

Order: the new unit test and the CLI contract -> fixture gates (s1b_fixtures.py, two parts) -> the end-to-end suite ->
legacy-facing focused suites -> one CYP2D6 30x development sentinel under sequence_first_v1 (the existing perf cell,
P1's development overlay and frozen settings; no truth is read) and its legacy-identity comparison against P1.1's saved
consolidated run of the same cell -> every other enabled CTest test, one per unit, heaviest last.

  s1b_run.py manifest | next | <unit>
"""
import csv
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "p1"))
sys.path.insert(0, HERE)
import p1_run  # noqa: E402  (cells, profile composition, arguments, status reader, run_timed)

WT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
BUILD = f"{WT}/build"
OUT = "/Users/davide.bolognini/panvar/tests/results/s1b_sequence_first_20261018"
P11 = "/Users/davide.bolognini/panvar/tests/results/p11_validation_closure_20261017/runs/o.cyp_30"
BUDGET = 1200
LIGHT = ["review_regressions", "dosage", "headline", "downstream_schema", "curation_contract", "prepared_locus_file",
         "batch_resource", "source_tree_clean", "library_sensitivity", "block_quality"]
REST_SMALL = ["mosaic_core", "locus_adapter", "fragment_oracle", "fragment_scorer", "scoring_context", "mosaic_proposal",
              "grouped_beam", "mosaic_binding", "heldout_mosaic_manifest", "mosaic_truth_orientation", "heldout_binding",
              "lzo_consolidation", "subsample_nested", "ladder_selection", "loo_development_gates", "pair_override",
              "production_campaign_plan", "v5_profile", "v6_profile", "repair_pilot_invariance"]
OTHER = ["library_estimate", "smoke_synthetic", "associate_stats", "bubble_stats", "inspect_stats", "call_stats",
         "benchmark_stats", "describe_stats", "panphorte_stats", "refine_stats", "smoke_real_c4"]
HEAVY = [("genotype_seeding_exactness", 60), ("genotype_union_placement_cache", 180), ("genotype_block_delta", 240),
         ("genotype_branch_lifetime", 280), ("genotype_report", 280)]


def test_name(short):
    return short if short.startswith(("genotype_", "library_", "smoke_")) or short.endswith("_stats") else "genotype_" + short


MANIFEST = [("t.genotype_sequence_report", "ctest", "genotype_sequence_report", 30),
            ("t.genotype_cli_contract", "ctest", "genotype_cli_contract", 30),
            ("fixtures.A", "fixtures", "F1,F2,F3,F4", 240),
            ("fixtures.B", "fixtures", "F5,F6,F7,F7full,F9", 240),
            ("t.genotype_command_e2e", "ctest", "genotype_command_e2e", 300)]
MANIFEST += [(f"t.{test_name(t)}", "ctest", test_name(t), 60) for t in LIGHT]
MANIFEST += [("cyp30.v1", "sentinel", "cyp_30", 150), ("cyp30.compare", "compare", "cyp_30", 30)]
MANIFEST += [(f"t.{test_name(t)}", "ctest", test_name(t), 60) for t in REST_SMALL + OTHER]
MANIFEST += [(f"t.{t}", "ctest", t, s) for t, s in HEAVY]
# Amendment 1 (committed after fixtures.A, before any further execution; the units already run stay recorded):
#   fixture variants found the caller's sequence-pair prior differing in its last bits with sample order, which made the
#   17-digit alternatives table differ and could make a bitwise "exact tie" depend on panel order. The implementation
#   now takes the canonical SHA-256 minimum over the whole (numerically tied) decision set and writes score components
#   with 12 significant digits; the oracle follows the same representative rule and orients labellings by sequence
#   (its first version compared two distinct sequences' labellings across haplotypes); the harness declares the
#   status rows for the dosage-stage copy call.original.mosaic_status.tsv too. The C++ changed, so the unit test and the
#   CLI contract rerun; fixtures A and B run as .a1 units into fixtures_a1 (fixtures.B is superseded, never run).
AMENDMENT_1 = [("t.genotype_sequence_report.a1", "ctest", "genotype_sequence_report", 30),
               ("t.genotype_cli_contract.a1", "ctest", "genotype_cli_contract", 30),
               ("fixtures.A.a1", "fixtures", "F1,F2,F3,F4", 240),
               ("fixtures.B.a1", "fixtures", "F5,F6,F7,F7full,F9", 240)]
SUPERSEDED = {"fixtures.B"}
ORDER = [m[0] for m in MANIFEST[:3]] + [m[0] for m in AMENDMENT_1] + \
    [m[0] for m in MANIFEST[3:] if m[0] not in SUPERSEDED]
LEDGER = ("unit\trevision\tdirty_cpp\tdeclared_timeout_seconds\twall_seconds\tpeak_rss_gib\texit\ttimed_out\tnote\n")


def spent():
    p = f"{OUT}/ledger.tsv"
    return sum(float(r["wall_seconds"]) for r in csv.DictReader(open(p), delimiter="\t")) if os.path.exists(p) else 0.0


def done():
    p = f"{OUT}/ledger.tsv"
    return {r["unit"] for r in csv.DictReader(open(p), delimiter="\t")} if os.path.exists(p) else set()


def state():
    rev = subprocess.run(["git", "-C", WT, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "-C", WT, "status", "--porcelain", "--", "src", "include", "CMakeLists.txt", "tests"],
                           capture_output=True, text=True).stdout.strip()
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


def cyp_compare():
    import s1b_fixtures
    d = f"{OUT}/runs/cyp30.v1"
    problems = s1b_fixtures.legacy_identity(P11, d)
    loc = s1b_fixtures.kv(f"{d}/call.sequence_locus.tsv")
    print("legacy identity against P1.1 o.cyp_30 (392fc0b build, legacy schema):",
          "IDENTICAL (declared additions only)" if not problems else "DIFFERENT: " + " | ".join(problems[:8]))
    for k in ("prediction", "sequence_status", "representation_intervals", "decision_set_groups", "exact_score_ties",
              "representative_tie_break", "identified_blocks", "representation_ambiguous_blocks",
              "genotype_unresolved_blocks", "representative_probability", "decision_set_probability", "legacy_run_status",
              "legacy_no_call_reason", "exit_code"):
        print(f"{k}\t{loc.get(k)}")
    return 1 if problems else 0


def run_unit(unit):
    u, kind, target, timeout = next(m for m in MANIFEST + AMENDMENT_1 if m[0] == unit)
    assert open(f"{OUT}/MANIFEST.tsv").read() and timeout <= 300 and unit not in done() and unit not in SUPERSEDED
    if any(m[0] == unit for m in AMENDMENT_1):
        assert os.path.exists(f"{OUT}/AMENDMENT_1.tsv"), "freeze amendment 1 first"
    rev, dirty = state()
    assert dirty == 0, "the implementation tree must be committed and clean"
    used = spent()
    if used + timeout > BUDGET:
        record(unit, rev, dirty, timeout, 0.0, 0.0, "not_run", 0, f"NOT RUN: {used:.1f} s + {timeout} s > {BUDGET} s")
        print(f"{unit}: NOT RUN ({used:.1f} + {timeout} > {BUDGET})")
        return 3
    d = f"{OUT}/runs/{unit}"
    assert not os.path.exists(d), f"{d} exists; a run is never repeated into the same place"
    os.makedirs(d)
    if kind == "ctest":
        cmd = ["nice", "ctest", "--test-dir", BUILD, "-R", f"^panvar_{target}$", "--output-on-failure", "-j1"]
    elif kind == "fixtures":
        fixture_dir = f"{OUT}/fixtures_a1" if unit.endswith(".a1") else f"{OUT}/fixtures"
        cmd = ["nice", sys.executable, os.path.join(HERE, "s1b_fixtures.py"), fixture_dir, f"{BUILD}/panvar",
               f"{BUILD}/cyp2d6_block_alleles", f"{BUILD}/genotype_prepared_locus_driver", target]
    elif kind == "compare":
        cmd = [sys.executable, os.path.abspath(__file__), "_compare"]
    else:
        cell = p1_run.cells()[target]
        prof = p1_run.profile_for(cell, "opt", d, timeout, target)
        cmd = ["nice", f"{BUILD}/panvar", "genotype", *p1_run.genotype_args(cell, prof, f"{d}/call"),
               "--output-schema", "sequence_first_v1"]
    cmd = ["/usr/bin/time", "-l", *cmd]
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
    note = "-"
    text = open(f"{d}/run.log").read()
    if kind == "ctest":
        m = re.search(r"Test +#\d+: \S+ \.+\s*\**\s*(\w[\w ]*?)\s+([\d.]+) sec", text)
        note = f"ctest: {m.group(1)} {m.group(2)} s" if m else "ctest: no result line"
    elif kind == "sentinel":
        st = p1_run.status_of(f"{d}/call")
        note = f"legacy status {st.get('status')}/{st.get('no_call_reason')}; sequence {st.get('sequence_prediction')}/" \
               f"{st.get('sequence_status')}; seconds_total {st.get('seconds_total')}"
    record(unit, rev, dirty, timeout, wall, rss / 2**30, code, timed_out, note)
    print(f"{unit}: exit {code}{' TIMEOUT' if timed_out else ''}, {wall:.1f} s, {note} (spent {spent():.1f}/{BUDGET})")
    return 0


if __name__ == "__main__":
    arg = sys.argv[1]
    if arg == "manifest":
        write_manifest()
    elif arg == "_compare":
        sys.exit(cyp_compare())
    elif arg == "amendment1":
        text = "".join("\t".join(str(x) for x in m) + "\n" for m in AMENDMENT_1) + "fixtures.B\tSUPERSEDED\n"
        p = f"{OUT}/AMENDMENT_1.tsv"
        assert not os.path.exists(p) or open(p).read() == text
        open(p, "w").write(text)
        print(text, end="")
    elif arg == "next":
        todo = [u for u in ORDER if u not in done()]
        sys.exit(run_unit(todo[0]) if todo else print("all units recorded") or 0)
    else:
        sys.exit(run_unit(arg))
