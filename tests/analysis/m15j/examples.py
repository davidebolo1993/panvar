#!/usr/bin/env python3
"""M15-J units (DESIGN.md section 4): ordinary runs and the dosage check unchanged against M15-G, removed options and forbidden
profile settings refused before the output prefix is touched, the panel index unchanged against M15-E, and the unrepaired
view unchanged against the M15-I binary. Usage: examples.py <example> <unit dir>"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "m15d"))
import examples as E  # noqa: E402
sys.path.insert(0, os.path.join(HERE, "..", "m15e"))
import importlib.util  # noqa: E402
_spec = importlib.util.spec_from_file_location("m15e_examples", os.path.join(HERE, "..", "m15e", "examples.py"))
M15E_X = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(M15E_X)

R, WT = E.R, E.WT
M15G = f"{R}/m15g_assessment_20261022/runs"
M15E = f"{R}/m15e_panel_index_20261021/runs"
OUT = f"{R}/m15j_interface_20261022"
BASELINE = f"{OUT}/baseline/panvar.m15i"   # the M15-I candidate binary (SHA-256 3f7e03a4...), copied before any M15-J change
CHECKS = []
# The M15-H rule: timing, memory and the exact caches' scheduling-dependent served-work counters are volatile
VOLATILE = ("seconds", "_bytes", "rss_", "footprint", "elapsed", "threads_peak", "_peak", "emission_cache_",
            "union_prefilter_cache_", "_executed_dps", "cached_prefilters_consumed")
KEYED = (".mosaic_status.tsv", ".dosage_status.tsv", ".library_estimate.tsv", ".library_sensitivity.tsv")
LOGS = (".repair_memory.tsv", ".repair_rounds.tsv")


def check(ok, what):
    CHECKS.append((ok, what))
    print(("PASS " if ok else "FAIL ") + what)
    return ok


def run(a):
    return E.run(a)


def norm(text, prefix):
    return text.replace(prefix, "<RUN>").replace(os.path.basename(prefix) + "_seq", "<RUN>_seq")


def same_outputs(new, old, what):
    """Every output of `new` equals `old`'s, apart from the run's own name, timing, memory and volatile counters."""
    nd, od = os.path.dirname(new), os.path.dirname(old)
    nb, ob = os.path.basename(new), os.path.basename(old)
    a = sorted(f[len(nb):] for f in os.listdir(nd) if f.startswith(nb + "."))
    b = sorted(f[len(ob):] for f in os.listdir(od) if f.startswith(ob + "."))
    check(a == b, f"{what}: the same output files ({len(a)})")
    keep = lambda t: [l for l in t.splitlines() if not any(v in l.split("\t")[0] for v in VOLATILE)]
    for suffix in a:
        if suffix in LOGS or suffix not in b:
            continue
        x, y = norm(open(new + suffix, errors="replace").read(), new), norm(open(old + suffix, errors="replace").read(), old)
        if suffix in KEYED:
            check(keep(x) == keep(y), f"{what}: {suffix} identical apart from timing and resource rows")
        else:
            check(x == y, f"{what}: {suffix} identical")


def ordinary(cell, d):
    out = f"{d}/{cell}"
    p = run(E.cmd(cell, out))
    check(p.returncode == 0, f"{cell}: exit 0 with automatic library estimation")
    st = E.KV(out + ".mosaic_status.tsv")
    check(st.get("library_provenance", "").startswith("auto:"), f"{cell}: the library was estimated ({st.get('library_provenance', '')[:60]})")
    same_outputs(out, f"{M15G}/a3.ex.complete_{cell}/{cell}", f"{cell} vs M15-G")
    s = E.settings(out)
    check(s.get("output_schema", ("", ""))[0] == "sequence_first_v2", f"{cell}: settings record output_schema sequence_first_v2 "
          f"(source {s.get('output_schema', ('', ''))[1]})")
    return 0


def ex_ordinary_c4(d):
    return ordinary("c4", d)


def ex_ordinary_cyp2d6(d):
    return ordinary("cyp2d6", d)


def ex_ordinary_acot(d):
    return ordinary("acot", d)


def ex_dosage_check_cyp2d6(d):
    """The M15-G ex.dosage_cyp2d6 command (declared library, --dosage-check): outputs unchanged; and its prediction equal to the
    same command without --dosage-check."""
    out, plain = f"{d}/cyp2d6", f"{d}/cyp2d6_plain"
    p = run(E.cmd("cyp2d6", out, "declared", ["--dosage-check"]))
    check(p.returncode == 0, "cyp2d6 --dosage-check: exit 0")
    same_outputs(out, f"{M15G}/ex.dosage_cyp2d6/cyp2d6", "cyp2d6 --dosage-check vs M15-G")
    q = run(E.cmd("cyp2d6", plain, "declared"))
    check(q.returncode == 0, "cyp2d6 without --dosage-check: exit 0")
    unchanged(out, plain, "cyp2d6")
    return 0


def unchanged(out, plain, cell):
    check(E.prediction_key(out) == E.prediction_key(plain), f"{cell}: --dosage-check leaves the prediction unchanged")
    for suffix in (".sequence_haplotypes.tsv", ".mosaic_scores.tsv", ".mosaic_candidates.tsv", ".equivalence_set.tsv",
                   ".sequence_alternatives.tsv"):
        a, b = M15E_X.named(out, suffix), M15E_X.named(plain, suffix)
        check(a is not None and a == b, f"{cell}: {suffix} unchanged by --dosage-check")
    dc = E.KV(out + ".sequence_locus.tsv").get("dosage_check")
    dp = E.KV(plain + ".sequence_locus.tsv").get("dosage_check")
    check(dp == "NOT_REQUESTED" and dc not in (None, "NOT_REQUESTED"), f"{cell}: dosage_check {dc} with the flag, {dp} without")


def ex_dosage_check_c4(d):
    out = f"{d}/c4"
    p = run(E.cmd("c4", out, extra=["--dosage-check"]))
    check(p.returncode == 0, "c4 --dosage-check: exit 0")
    unchanged(out, f"{OUT}/runs/ex.ordinary_c4/c4", "c4")
    return 0


# Every removed spelling, forbidden profile setting and repeated profile: exit 2, the message, nothing touched.
REFUSALS = [
    (["--output-schema", "legacy"], "--output-schema was removed"),
    (["--output-schema", "sequence_first_v2"], "--output-schema was removed"),
    (["--dosage-mode", "repair"], "--dosage-mode was removed"),
    (["--dosage-mode", "audit"], "--dosage-mode was removed"),
    (["--dosage-audit", "shadow"], "--dosage-audit was removed"),
    (["--phase-audit-candidates", "8"], "set the profile key phase_audit_candidates"),
    (["--repair-origin-cache-bytes", "1048576"], "set the profile key repair_origin_cache_bytes"),
    (["--repair-tie-seed-cache-bytes", "1048576"], "set the profile key repair_tie_seed_cache_bytes"),
    (["--dosage-allocation-test"], "--dosage-allocation-test was removed"),
    (["--library-diagnostic"], "--library-diagnostic was removed"),
    (["--max-candidates", "4"], "--max-candidates is now --max-mosaic-candidates"),
    ("profile:output_schema\tlegacy", "'output_schema' cannot select a reporting contract"),
    ("profile:output_schema\tsequence_first_v1", "'output_schema' cannot select a reporting contract"),
    ("profile:dosage_mode\trepair", "'dosage_mode' is not a profile setting"),
    ("profile:dosage_audit\tshadow", "'dosage_audit' is not a profile setting"),
    ("twice", "--profile was given 2 times"),
]


def ex_refusals(d):
    sentinel = "M15-J sentinel: a file a refused run must leave untouched\n"
    for i, (extra, needle) in enumerate(REFUSALS):
        case = f"{d}/case{i:02d}"
        os.makedirs(case, exist_ok=True)
        prefix = f"{case}/c4"
        for suffix in (".locus_genotype.tsv", ".mosaic_status.tsv", ".sequence_locus.tsv"):
            open(prefix + suffix, "w").write(sentinel)
        before = sorted(os.listdir(case))
        if isinstance(extra, str) and extra.startswith("profile:"):
            prof = f"{case}/profile.tsv"
            open(prof, "w").write("profile\tm15j_forbidden\n" + extra[len("profile:"):] + "\n")
            args, label = ["--profile", prof], extra[len("profile:"):].replace("\t", " ")
            before = sorted(os.listdir(case))
        elif extra == "twice":
            prof = f"{case}/profile.tsv"
            open(prof, "w").write("profile\tm15j_plain\n")
            args, label = ["--profile", prof, "--profile", prof], "--profile twice"
            before = sorted(os.listdir(case))
        else:
            args, label = list(extra), " ".join(extra)
        p = run(E.cmd("c4", prefix, extra=args))
        check(p.returncode == 2 and needle in p.stderr, f"refused, exit 2, with the message: {label} ({p.returncode})")
        untouched = sorted(os.listdir(case)) == before and all(
            open(prefix + s).read() == sentinel for s in (".locus_genotype.tsv", ".mosaic_status.tsv", ".sequence_locus.tsv"))
        check(untouched, f"the output prefix untouched: {label}")
    return 0


def index_equal_but_provenance(new, old):
    """A1: the index equals M15-E's apart from the builder's recorded source revision (40 hex characters in the header) and the
    32-byte header digest that covers it. Every other byte, the payload included, must be identical."""
    a, b = open(new, "rb").read(), open(old, "rb").read()
    check(len(a) == len(b), f"c4: the index has M15-E's size ({len(a)} bytes)")
    diff = [i for i in range(min(len(a), len(b))) if a[i] != b[i]]
    hexd = set(b"0123456789abcdef")
    def run_of_hex(x, i):
        lo = i
        while lo > 0 and x[lo - 1] in hexd:
            lo -= 1
        hi = i
        while hi < len(x) and x[hi] in hexd:
            hi += 1
        return lo, hi
    rev = set()
    for i in diff:
        if 24 <= i < 56:
            continue
        lo, hi = run_of_hex(a, i)
        if (lo, hi) == run_of_hex(b, i) and hi - lo == 40:
            rev.add((lo, hi))
    allowed = lambda i: 24 <= i < 56 or any(lo <= i < hi for lo, hi in rev)
    check(diff and all(allowed(i) for i in diff) and len(rev) == 1,
          f"c4: the index equals M15-E's apart from the builder revision field {sorted(rev)} and the header digest "
          f"[24, 56) ({len(diff)} bytes differ)")
    if rev:
        lo, hi = next(iter(rev))
        print(f"builder revision recorded: M15-J {a[lo:hi].decode()} vs M15-E {b[lo:hi].decode()}")


def ex_index_c4(d):
    """The M15-E index build and load on C4 with the M15-J binary: the index bytes equal M15-E's index; the indexed run equals
    the cold run (the M15-E checks) and both equal M15-E's runs of the same commands."""
    idx_dir = f"{d}/index"
    os.makedirs(idx_dir, exist_ok=True)
    idx = f"{idx_dir}/c4.pvpl"
    p = run(M15E_X.locus_cmd("c4", ["--build-panel-index", idx]))
    check(p.returncode == 0 and os.path.exists(idx), "c4: index built")
    old_idx = f"{M15E}/ex.build_c4/index/c4.pvpl"
    index_equal_but_provenance(idx, old_idx)
    cold, ind = f"{d}/c4.cold", f"{d}/c4.index"
    p1 = run(E.cmd("c4", cold))
    p2 = run(E.cmd("c4", ind, extra=["--panel-index", idx]))
    check(p1.returncode == p2.returncode == 0, f"c4: cold and indexed exit 0 ({p1.returncode}, {p2.returncode})")
    M15E_X.CHECKS = CHECKS
    M15E_X.check = check
    M15E_X.compare("c4", cold, ind, "c4 cold vs indexed (M15-J)")
    for mine, theirs, what in ((cold, f"{M15E}/ex.index_c4/c4.cold", "cold"), (ind, f"{M15E}/ex.index_c4/c4.index", "indexed")):
        for suffix in M15E_X.SCIENCE:
            a, b = M15E_X.named(mine, suffix), M15E_X.named(theirs, suffix)
            if suffix in (".sequence_blocks.tsv", ".block_genotypes.tsv", ".genotypes.tsv"):
                continue   # M15-G added the assessment columns after M15-E; compared against M15-G in ex.ordinary_c4
            check(a == b, f"c4 {what}: {suffix} identical to M15-E")
        check(E.prediction_key(mine) == E.prediction_key(theirs), f"c4 {what}: prediction identical to M15-E")
    return 0


def ex_unrefined(d):
    """--adaptive-repair off on C4: the M15-J binary and the M15-I binary give the same outputs (behaviour unchanged pending
    decision 2.2); the settings record the preset's v2 as inactive with the reason."""
    new, old = f"{d}/new/c4", f"{d}/m15i/c4"
    os.makedirs(os.path.dirname(new), exist_ok=True)
    os.makedirs(os.path.dirname(old), exist_ok=True)
    a = E.cmd("c4", new, extra=["--adaptive-repair", "off"])
    b = E.cmd("c4", old, extra=["--adaptive-repair", "off"])
    b[0] = BASELINE
    p, q = run(a), run(b)
    check(p.returncode == q.returncode, f"c4 unrepaired: the same exit from both binaries ({p.returncode}, {q.returncode})")
    same_outputs(new, old, "c4 unrepaired (M15-J vs M15-I binary)")
    s = E.settings(new)
    schema = s.get("output_schema", ("", ""))
    check(schema[0] == "legacy" and "(inactive: the unrepaired view has no pooled decision to report)" in schema[1],
          f"c4 unrepaired: output_schema legacy, recorded inactive with the reason ({schema[1]})")
    st = E.KV(new + ".mosaic_status.tsv")
    print(f"unrepaired c4: status {st.get('status')}, call_class {st.get('call_class')}, phase_support {st.get('phase_support')}, "
          f"sequence report present: {os.path.exists(new + '.sequence_locus.tsv')}")
    return 0


if __name__ == "__main__":
    name, d = sys.argv[1], sys.argv[2]
    os.makedirs(d, exist_ok=True)
    globals()["ex_" + name](d)
    bad = [w for ok, w in CHECKS if not ok]
    print(f"{len(CHECKS) - len(bad)} of {len(CHECKS)} checks passed" + (f"; FAILED: {bad}" if bad else ""))
    sys.exit(1 if bad else 0)
