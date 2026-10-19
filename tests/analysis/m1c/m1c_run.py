#!/usr/bin/env python3
"""M1-C runner -- isolate C4 allocation from sequence-content differences (frozen manifest; one process at a time; <= 4
threads; <= 300 s per unit; 1200 s budget counted as time spent in this batch -- development runs included -- plus the
declared timeout of the next unit; 6 GiB envelope). Diagnostic evaluator only: no caller run, no production change.

Order: mass-reconciliation gates (two synthetic fixtures, then the M1-B C4 10x decomposition rerun) -> spelling (C4 10x,
30x: the M1-B candidates and the truth pair) -> explicit-sequence gate (the M1-B 10x candidates re-declared as explicit
sequences plus controls, compared bit for bit with the path-based run) -> construction of the frozen allocation
alternatives -> explicit scoring with decompositions (only if the explicit gate passed) -> content -> summary.

  m1c_run.py manifest B | next B | <unit> | dev <name> <note> -- <command...>
"""
import csv
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "p1"))
import p1_run  # noqa: E402

WT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
BUILD = f"{WT}/build"
OUT = "/Users/davide.bolognini/panvar/tests/results/m1c_allocation_20261019"
M1B = "/Users/davide.bolognini/panvar/tests/results/m1b_structural_20261019"
BUDGET = 1200
CASES = {r["case"]: r for r in csv.DictReader(open(os.path.join(HERE, "..", "m1b", "AUDIT_SET.tsv")), delimiter="\t")}
MANIFEST_B = [("gate.mass_fixture.s1", "fixture", "20261019", 30), ("gate.mass_fixture.s2", "fixture", "7", 30),
              ("gate.mass_c4_10", "massc4", "c4_10", 120),
              ("spell.c4_10", "spell", "c4_10", 60), ("spell.c4_30", "spell", "c4_30", 60),
              ("gate.explicit_set", "valset", "c4_10", 20), ("gate.explicit_c4_10", "valscore", "c4_10", 120),
              ("gate.explicit_compare", "valcompare", "c4_10", 20),
              ("construct.c4_10", "construct", "c4_10", 60), ("construct.c4_30", "construct", "c4_30", 60),
              ("score.c4_10", "score", "c4_10", 150), ("score.c4_30", "score", "c4_30", 300),
              ("content.c4_10", "content", "c4_10", 90), ("content.c4_30", "content", "c4_30", 90),
              ("summary", "summary", "-", 30)]
# Amendment 1 (declared before further execution): both fixture units failed before scoring on a harness defect -- the
# fixture's LibraryModel had no provenance string, which validate_library_model refuses (recorded, unchanged). The fixture
# now declares one; nothing else changed. New units rerun the two fixtures; the rest of manifest B continues unchanged.
MANIFEST_A1 = [("gate.mass_fixture.s1.a1", "fixture", "20261019", 30), ("gate.mass_fixture.s2.a1", "fixture", "7", 30)]
# Amendment 2 (declared before further execution): score.c4_30 would record every origin of 22 unique haplotypes at 30x;
# M1-B's 16-haplotype 30x decomposition peaked at 4.0 GB, so this projects to about 5.5 GB, too close to the 6 GiB
# envelope. score.c4_30 is WITHDRAWN before execution and replaced by two scorer runs over disjoint alternative subsets
# that share `reported` (a candidate's score does not depend on the others scored beside it -- validated bit for bit by
# the explicit gate); `reported` must be bit-equal across the two parts. Content and summary units are unchanged in
# meaning; the summary reads both parts.
PARTS = {"1": ["reported", "realloc_A_p1", "realloc_A_p2", "realloc_B_p1", "realloc_B_p2"],
         "2": ["reported", "saved_2+2", "truth", "truthblock_in_reported_a", "truthblock_in_reported_b",
               "reportedblock_in_truth_a", "reportedblock_in_truth_b"]}
MANIFEST_A2 = [("score.c4_30.part1", "scorepart", "c4_30:1", 150), ("score.c4_30.part2", "scorepart", "c4_30:2", 200)]
WITHDRAWN = {"score.c4_30": "WITHDRAWN by Amendment 2 before execution (memory envelope); replaced by score.c4_30.part1/part2"}
MANIFESTS = {"B": MANIFEST_B, "A1": MANIFEST_A1, "A2": MANIFEST_A2}


def ledger_rows():
    p = f"{OUT}/ledger.tsv"
    return list(csv.DictReader(open(p), delimiter="\t")) if os.path.exists(p) else []


def spent():
    return sum(float(r["wall_seconds"]) for r in ledger_rows())


def done():
    return {r["unit"] for r in ledger_rows()}


def exit_of(unit):
    return next((r["exit"] for r in ledger_rows() if r["unit"] == unit), None)


def state():
    rev = subprocess.run(["git", "-C", WT, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "-C", WT, "status", "--porcelain"], capture_output=True, text=True).stdout.strip()
    return rev, len([line for line in dirty.split("\n") if line.strip()])


def write_manifest(part):
    os.makedirs(OUT, exist_ok=True)
    rows = MANIFESTS[part]
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


def reported_total(case):
    rows = list(csv.DictReader(open(f"{M1B}/runs/call.{case}/call.equivalence_set.tsv"), delimiter="\t"))
    return next(r for r in rows if r["rank"] == "1")["score"]


def comparator_total(case):
    for line in open(f"{M1B}/runs/alts.{case}/{case}.candidates.tsv"):
        f = line.rstrip("\n").split("\t")
        if f[0] == "saved_2+2":
            return f[3]
    raise KeyError("saved_2+2")


def scorer(c, case, d, extra):
    return ["nice", f"{BUILD}/genotype_probe", "score-set", "--out", f"{d}/score", "--reference", c["reference"],
            "--donor", c["sample"], "--run", f"{M1B}/runs/call.{case}/call", "--model-gfa", c["model_graph"],
            "--library", c["library"], "--reads1", c["reads1"], "--reads2", c["reads2"], "--threads", "4", *extra]


def scorepart(target, d):
    if exit_of("gate.explicit_compare") != "0":
        return None, "NOT RUN (condition): the explicit-sequence gate did not pass"
    case, part = target.split(":")
    c = CASES[case]
    keep = PARTS[part]
    src = f"{OUT}/runs/construct.{case}/alternatives.explicit.tsv"
    rows = [l for l in open(src) if not l.startswith("#") and l.strip() and l.split("\t")[0] in keep]
    open(f"{d}/alternatives.part{part}.explicit.tsv", "w").write("".join(rows))
    labels = [l.split("\t")[0] for l in rows]
    pairs = [f"reported,{x}" for x in labels if x.startswith("realloc_")]
    pairs += [f"reported,{x}" for x in ("truth", "saved_2+2", "truthblock_in_reported_a") if x in labels]
    if "reportedblock_in_truth_a" in labels:
        pairs.append("truth,reportedblock_in_truth_a")
    return scorer(c, case, d, ["--explicit", f"{d}/alternatives.part{part}.explicit.tsv", "--decompose", ";".join(pairs),
                               "--decompose-block", c["target_block"]]), f"explicit allocation alternatives, part {part}"


def command(unit, kind, target, timeout, d):
    if kind == "fixture":
        return ["nice", f"{BUILD}/genotype_probe", "mass-reconcile", "--out", f"{d}/fixture", "--seed", target], "synthetic fixture"
    if kind == "summary":
        return [sys.executable, os.path.join(HERE, "summary.py"), f"{OUT}/runs", f"{d}/m1c"], "-"
    if unit in WITHDRAWN:
        return None, WITHDRAWN[unit]
    if kind == "scorepart":
        return scorepart(target, d)
    c = CASES[target]
    cands = f"{M1B}/runs/alts.{target}/{target}.candidates.tsv"
    if kind == "massc4":
        return scorer(c, target, d, ["--candidates", cands, "--decompose", "reported,saved_2+2", "--decompose-block",
                                     c["target_block"]]), "M1-B decomposition rerun with mass reconciliation"
    if kind == "spell":
        return ["nice", f"{BUILD}/genotype_probe", "spell", "--out", f"{d}/spell", "--reference", c["reference"],
                "--model-gfa", c["model_graph"], "--candidates", cands, "--block", c["target_block"], "--truth-gfa",
                c["truth_graph"], "--donor", c["sample"]], "diagnostic (truth)"
    if kind == "valset":
        return [sys.executable, os.path.join(HERE, "validation_set.py"), "--candidates", cands, "--spell",
                f"{OUT}/runs/spell.{target}/spell", "--out", f"{d}/validation"], "-"
    if kind == "valscore":
        return scorer(c, target, d, ["--explicit", f"{OUT}/runs/gate.explicit_set/validation.explicit.tsv"]), "explicit gate"
    if kind == "valcompare":
        return [sys.executable, os.path.join(HERE, "compare_explicit.py"), f"{OUT}/runs/gate.mass_c4_10/score",
                f"{OUT}/runs/gate.explicit_c4_10/score", f"{d}/compare.tsv"], "explicit gate comparison"
    if kind == "construct":
        return [sys.executable, os.path.join(HERE, "construct.py"), "--spell", f"{OUT}/runs/spell.{target}/spell", "--out",
                f"{d}/alternatives", "--expected-reported", reported_total(target), "--expected-comparator",
                comparator_total(target)], "diagnostic construction (truth)"
    if kind == "score":
        if exit_of("gate.explicit_compare") != "0":
            return None, "NOT RUN (condition): the explicit-sequence gate did not pass"
        alts = f"{OUT}/runs/construct.{target}/alternatives.explicit.tsv"
        labels = [l.split("\t")[0] for l in open(alts) if not l.startswith("#") and l.strip()]
        pairs = [f"reported,{x}" for x in labels if x.startswith("realloc_")]
        pairs += [f"reported,{x}" for x in ("truth", "saved_2+2", "truthblock_in_reported_a") if x in labels]
        if "reportedblock_in_truth_a" in labels:
            pairs.append("truth,reportedblock_in_truth_a")
        return scorer(c, target, d, ["--explicit", alts, "--decompose", ";".join(pairs), "--decompose-block",
                                     c["target_block"]]), "explicit allocation alternatives"
    if kind == "content":
        return [sys.executable, os.path.join(HERE, "content.py"), "--explicit",
                f"{OUT}/runs/construct.{target}/alternatives.explicit.tsv", "--out", f"{d}/{target}"], "diagnostic (truth)"
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
