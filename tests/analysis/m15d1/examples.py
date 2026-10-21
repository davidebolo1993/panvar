#!/usr/bin/env python3
"""M15-D.1: binary-level reproductions of the basic-command corrections, on the development C4 and CYP2D6 cells (no
reserved family). Each example runs `panvar genotype` and asserts the frozen expectations (DESIGN.md section 2).
Usage: examples.py <example> <unit dir>"""
import os
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "m15d"))
import examples as E  # noqa: E402  (the M15-D cells, readers and prediction key)

R = E.R
M15D = f"{R}/m15d_basic_command_20261021/runs"
CHECKS = []


def check(ok, what):
    CHECKS.append((ok, what))
    print(("PASS " if ok else "FAIL ") + what)


def run(a):
    print("$", " ".join(a).replace(R, "<results>").replace(E.WT, "<worktree>"))
    p = subprocess.run(a, capture_output=True, text=True)
    print(f"exit {p.returncode}")
    if p.stderr.strip():
        print("stderr:", p.stderr.strip()[-700:])
    return p


def profile(d, name, rows):
    f = f"{d}/{name}.tsv"
    open(f, "w").write(f"profile\tm15d1-{name}\n" + "".join(f"{k}\t{v}\n" for k, v in rows))
    return f


def exits_agree(prefix, code, prediction, what, status_exit=True):
    """The process exit, .sequence_locus.tsv exit_code, the status's sequence_exit_code and (stops, failures) exit_code."""
    sl, st = E.KV(prefix + ".sequence_locus.tsv"), E.KV(prefix + ".mosaic_status.tsv")
    c = str(code)
    print(f"{what}: locus exit_code {sl.get('exit_code')}, status sequence_exit_code {st.get('sequence_exit_code')}, "
          f"status exit_code {st.get('exit_code')}, prediction {sl.get('prediction')}/{st.get('sequence_prediction')}, "
          f"status {st.get('status')}")
    check(sl.get("exit_code") == c and st.get("sequence_exit_code") == c and (not status_exit or st.get("exit_code") == c),
          f"{what}: every exit field says {c}")
    check(sl.get("prediction") == prediction and st.get("sequence_prediction") == prediction, f"{what}: {prediction} in report and status")
    return sl, st


def ex_library_preflight(d):
    """A missing, unreadable or invalid declared library is refused (exit 2) with the prefix untouched."""
    good = open(f"{E.CELLS['c4']['reads']}/library.txt").read()
    cases = [("missing", None, "cannot open library model file"), ("unreadable", good, "cannot open library model file"),
             ("repeated", good + "insert_sd_bp 40\n", "insert_sd_bp"), ("unknown", good + "gap_open 4\n", "gap_open")]
    sentinels = [".mosaic_status.tsv", ".sequence_locus.tsv", ".settings.tsv", ".genotypes.tsv", ".sequence_haplotypes.fa"]
    for name, content, needle in cases:
        sub = f"{d}/{name}"
        os.makedirs(sub, exist_ok=True)
        lib = f"{sub}/library.txt"
        if content is not None:
            open(lib, "w").write(content)
            if name == "unreadable":
                os.chmod(lib, 0)
        for s in sentinels:
            open(f"{sub}/c4{s}", "w").write(f"sentinel {s}\n")
        p = run(E.cmd("c4", f"{sub}/c4", lib))
        if name == "unreadable":
            os.chmod(lib, 0o600)
        kept = all(open(f"{sub}/c4{s}").read() == f"sentinel {s}\n" for s in sentinels)
        files = sorted(os.listdir(sub))
        want = sorted([f"c4{s}" for s in sentinels] + ([] if content is None else ["library.txt"]))
        check(p.returncode == 2 and needle in p.stderr, f"{name} library: exit 2, the message names '{needle}'")
        check(kept and files == want, f"{name} library: every file at the prefix byte-identical, nothing added ({files})")
    return 0


def stop(d, name, point, code, prediction):
    out = f"{d}/c4"
    p = run(E.cmd("c4", out, "declared", ["--profile", profile(d, name, [("test_interrupt_point", point)])]))
    check(p.returncode == code, f"{name}: process exit {code} (got {p.returncode})")
    return exits_agree(out, code, prediction, name)


def ex_term_time(d):
    sl, st = stop(d, "time_stop", "before_first_evaluation", 124, "NO_PREDICTION")
    check(st.get("status") == "TIME_LIMIT" and not os.path.exists(f"{d}/c4.sequence_haplotypes.tsv"), "time stop: TIME_LIMIT, no sequences")
    return 0


def ex_term_footprint(d):
    sl, st = stop(d, "footprint_stop", "footprint_before_first_evaluation", 125, "NO_PREDICTION")
    check(st.get("status") == "RESOURCE_REFUSAL" and st.get("no_call_reason") == "footprint_limit" and
          not os.path.exists(f"{d}/c4.sequence_haplotypes.tsv"), "footprint stop: RESOURCE_REFUSAL (footprint_limit), no sequences")
    return 0


def ex_term_checkpoint(d):
    sl, st = stop(d, "checkpoint_stop", "after_checkpoint", 3, "PROVISIONAL")
    check(sl.get("completion") == "budget_stop" and os.path.exists(f"{d}/c4.sequence_haplotypes.tsv"),
          "checkpoint stop: budget_stop, the checkpoint's sequences kept")
    return 0


def ex_term_nocall(d):
    """An ordinary no-prediction outcome: a declared neighbourhood bound of 1 with no early checkpoint to keep."""
    out = f"{d}/c4"
    p = run(E.cmd("c4", out, "declared", ["--profile", profile(d, "no_early", [("early_incumbent", "off")]),
                                         "--repair-max-neighbours", "1"]))
    check(p.returncode == 1, f"ordinary no prediction: process exit 1 (got {p.returncode})")
    sl, st = exits_agree(out, 1, "NO_PREDICTION", "ordinary no prediction", status_exit=False)
    print(f"no_call_reason {st.get('no_call_reason')}")
    return 0


def ex_term_locus_failure(d):
    """A reporting failure after the sequence report was published (the locus table cannot be written): 2."""
    out = f"{d}/c4"
    os.makedirs(f"{out}.locus_genotype.tsv.writing/blocker", exist_ok=True)
    p = run(E.cmd("c4", out, "declared"))
    check(p.returncode == 2, f"locus output failure: process exit 2 (got {p.returncode})")
    sl, st = E.KV(out + ".sequence_locus.tsv"), E.KV(out + ".mosaic_status.tsv")
    print(f"locus exit_code {sl.get('exit_code')}, status sequence_exit_code {st.get('sequence_exit_code')}, exit_code "
          f"{st.get('exit_code')}, locus_output {sl.get('locus_output', '')[:120]}")
    check(sl.get("exit_code") == "2" and st.get("sequence_exit_code") == "2" and st.get("exit_code") == "2",
          "locus output failure: every exit field says 2")
    check(sl.get("locus_output", "").startswith("FAILED: ") and st.get("locus_output", "").startswith("FAILED: "),
          "locus output failure: named in the report and the status")
    os.rmdir(f"{out}.locus_genotype.tsv.writing/blocker")
    os.rmdir(f"{out}.locus_genotype.tsv.writing")
    return 0


def ex_cyp_audit(d):
    """The known CYP2D6 audit limitation, reported honestly (its unit-library model is not changed here)."""
    out = f"{d}/cyp2d6"
    p = run(E.cmd("cyp2d6", out, "declared", ["--dosage-check"]))
    check(p.returncode == 0, f"cyp2d6 --dosage-check: exit 0 (got {p.returncode})")
    archived, absent = E.CELLS["cyp2d6"]["archived"], f"{M15D}/a2.ex.dosage_absent_cyp2d6/cyp2d6"
    check(E.prediction_key(out) == E.prediction_key(archived), "cyp2d6: prediction identical to the archived M11-D run")
    check(E.prediction_key(out) == E.prediction_key(absent), "cyp2d6: prediction identical to M15-D's run without the audit")
    ds = E.KV(out + ".dosage_status.tsv")
    print(f"dosage_audit_status {ds.get('dosage_audit_status')}: {ds.get('dosage_audit_note')}")
    check(ds.get("dosage_audit_status") == "FAILED" and "alleles disagree on their flanks or anchor trims" in ds.get("dosage_audit_note", ""),
          "cyp2d6: the audit reports FAILED with its concrete reason")
    loc = (E.T(out + ".locus_genotype.tsv") or [{}])[0]
    print(f"locus dosage_filter {loc.get('dosage_filter')} ({loc.get('dosage_filter_reasons', '')[:160]}); summary {loc.get('call_summary')}")
    blocks = E.T(out + ".block_genotypes.tsv")
    block_filters = sorted({r.get("dosage_filter") for r in blocks})
    print("block dosage_filter values:", block_filters)
    check(loc.get("dosage_filter") == "NOT_EVALUATED" and "audit_FAILED" in loc.get("dosage_filter_reasons", "") and
          "CONSISTENT" not in loc.get("dosage_filter_reasons", ""), "cyp2d6: locus DOSAGE_FILTER NOT_EVALUATED (audit_FAILED), never consistent")
    check(bool(blocks) and all(f in ("NOT_EVALUATED", "NO_CALL") for f in block_filters),
          f"cyp2d6: no block reports a dosage verdict ({block_filters})")
    check("WARNING: the requested dosage audit did not complete (FAILED" in p.stderr,
          "cyp2d6: the failure is visible on stderr when --dosage-check was requested")
    return 0


if __name__ == "__main__":
    name, d = sys.argv[1], sys.argv[2]
    os.makedirs(d, exist_ok=True)
    globals()["ex_" + name](d)
    bad = [w for ok, w in CHECKS if not ok]
    print(f"{len(CHECKS) - len(bad)} of {len(CHECKS)} checks passed" + (f"; FAILED: {bad}" if bad else ""))
    sys.exit(1 if bad else 0)
