#!/usr/bin/env python3
"""M15-L units (DESIGN.md section 5): each command runs with the M15-L binary and with the M15-K binary (baseline/panvar.m15k).
The M15-L run must keep exactly the primary reports beside the prefix and every other report in <prefix>.diagnostics/, and every
file must equal the M15-K file of the same name (the M15-H rule: names, timing, scheduling counters), apart from the declared
reference change (the locus table's main_fasta). Usage: examples.py <example> <unit dir>"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "m15d"))
import examples as E  # noqa: E402

R, WT = E.R, E.WT
OUT = f"{R}/m15l_diagnostics_dir_20261022"
BASELINE = f"{OUT}/baseline/panvar.m15k"   # the M15-K binary (5c4707a5...), copied before any M15-L change
M15K_INDEX = f"{R}/m15k_unrefined_v2_20261022/runs/a2.ex.index_c4/index/c4.pvpl"
PRIMARY = {".sequence_locus.tsv", ".sequence_blocks.tsv", ".sequence_haplotypes.tsv", ".sequence_haplotypes.fa"}
VOLATILE = ("seconds", "_bytes", "rss_", "footprint", "elapsed", "threads_peak", "_peak", "emission_cache_",
            "union_prefilter_cache_", "_executed_dps", "cached_prefilters_consumed")
KEYED = (".mosaic_status.tsv", ".dosage_status.tsv", ".library_estimate.tsv", ".library_sensitivity.tsv")
LOGS = (".repair_memory.tsv", ".repair_rounds.tsv")
CHECKS = []


def check(ok, what):
    CHECKS.append((ok, what))
    print(("PASS " if ok else "FAIL ") + what)
    return ok


def both(name, d, extra=(), library=None, cell="c4", setup=None):
    """The same command with both binaries, each in its own directory; returns (new prefix, old prefix, new exit, old exit)."""
    new, old = f"{d}/new/{cell}", f"{d}/m15k/{cell}"
    for x in (new, old):
        os.makedirs(os.path.dirname(x), exist_ok=True)
        if setup:
            setup(x)
    a = E.cmd(cell, new, library, list(extra))
    b = E.cmd(cell, old, library, list(extra))
    b[0] = BASELINE
    p, q = E.run(a), E.run(b)
    check(p.returncode == q.returncode, f"{name}: the same exit from both binaries ({p.returncode}, {q.returncode})")
    return new, old, p.returncode, q.returncode


def layout(name, new, old):
    """The M15-L layout of `new`, and every file equal to the M15-K file of the same name."""
    base = os.path.basename(new)
    beside = sorted(f[len(base):] for f in os.listdir(os.path.dirname(new)) if f.startswith(base + ".") and f != base + ".diagnostics")
    diag_dir = new + ".diagnostics"
    diag = sorted(f[len(base):] for f in os.listdir(diag_dir) if f.startswith(base + ".")) if os.path.isdir(diag_dir) else []
    old_files = sorted(f[len(base):] for f in os.listdir(os.path.dirname(old)) if f.startswith(base + "."))
    check(set(beside) <= PRIMARY, f"{name}: only primary reports beside the prefix ({beside})")
    check(not (set(diag) & PRIMARY), f"{name}: no primary report in the diagnostics directory")
    check(sorted(beside + diag) == old_files, f"{name}: the same {len(old_files)} reports as M15-K, by name "
          f"({len(beside)} beside the prefix, {len(diag)} in {base}.diagnostics/)")
    for suffix in old_files:
        if suffix in LOGS:
            continue
        mine = new + suffix if suffix in PRIMARY else f"{diag_dir}/{base}{suffix}"
        if not os.path.exists(mine):
            continue
        x = open(mine, errors="replace").read()
        y = open(old + suffix, errors="replace").read()
        # the run's own paths: the diagnostics prefix first, then the primary prefix, then M15-K's prefix
        x = x.replace(f"{diag_dir}/{base}", "<RUN>").replace(new, "<RUN>")
        y = y.replace(old, "<RUN>")
        if suffix == ".locus_genotype.tsv":   # the declared reference change (DESIGN.md section 3)
            x = x.replace(f"../{base}.sequence_haplotypes.fa", f"{base}.sequence_haplotypes.fa")
        if suffix in KEYED:
            keep = lambda t: [l for l in t.splitlines() if not any(v in l.split("\t")[0] for v in VOLATILE)]
            check(keep(x) == keep(y), f"{name}: {suffix} equal to M15-K (timing and resource rows aside)")
        else:
            check(x == y, f"{name}: {suffix} equal to M15-K")
    return beside, diag


def cell_unit(cell, d, extra=(), library=None):
    new, old, rc, _ = both(cell, d, extra, library, cell)
    beside, _ = layout(cell, new, old)
    check(".sequence_locus.tsv" in beside, f"{cell}: the outcome record beside the prefix")
    print(f"{cell}: exit {rc}; {E.KV(new + '.sequence_locus.tsv').get('assessment_summary')}")
    return 0


def ex_layout_c4(d):
    return cell_unit("c4", d)


def ex_layout_cyp2d6(d):
    return cell_unit("cyp2d6", d)


def ex_layout_acot(d):
    return cell_unit("acot", d)


def ex_dosage_check_cyp2d6(d):
    return cell_unit("cyp2d6", d, ["--dosage-check"], "declared")


def ex_indexed_c4(d):
    """--panel-index with the M15-K index, so both binaries load one file."""
    check(os.path.exists(M15K_INDEX), "the M15-K index exists")
    new, old, rc, _ = both("indexed c4", d, ["--panel-index", M15K_INDEX])
    layout("indexed c4", new, old)
    st = E.KV(f"{new}.diagnostics/c4.mosaic_status.tsv")
    check(st.get("prepared_locus_origin") == "loaded" and st.get("panel_index") == M15K_INDEX,
          "indexed c4: the status (in the diagnostics directory) records the loaded index")
    return 0


def stop_profile(d, name, point):
    prof = f"{d}/{name}.tsv"
    open(prof, "w").write(f"profile\tm15l-{name}\ntest_interrupt_point\t{point}\n")
    return prof


def ex_provisional(d):
    prof = stop_profile(d, "provisional", "after_checkpoint")
    new, old, rc, _ = both("provisional", d, ["--profile", prof], "declared")
    check(rc == 3, f"provisional: exit 3 ({rc})")
    beside, _ = layout("provisional", new, old)
    sl = E.KV(new + ".sequence_locus.tsv")
    check(sl.get("prediction") == "PROVISIONAL" and sl.get("exit_code") == "3",
          "provisional: the PROVISIONAL record beside the prefix, exit_code 3")
    return 0


def ex_failed(d):
    # (1) a watchdog stop before any evaluated candidate: exit 124, the stop record beside the prefix
    prof = stop_profile(d, "absent", "before_first_evaluation")
    new, old, rc, _ = both("absent", f"{d}/absent", ["--profile", prof], "declared")
    check(rc == 124, f"absent: exit 124 ({rc})")
    layout("absent", new, old)
    sl = E.KV(new + ".sequence_locus.tsv")
    check(sl.get("prediction") == "NO_PREDICTION" and sl.get("exit_code") == "124",
          "absent: the NO_PREDICTION stop record beside the prefix, exit_code 124")
    # (2) a reporting failure after publication: exit 2, recorded in the primary record (the locus table's write is blocked)
    def block_locus_table(prefix):
        base = os.path.basename(prefix)
        target = f"{prefix}.diagnostics/{base}.locus_genotype.tsv.writing/blocker" if "/new/" in prefix else \
                 f"{prefix}.locus_genotype.tsv.writing/blocker"
        os.makedirs(target, exist_ok=True)
    new, old, rc, _ = both("reporting failure", f"{d}/failure", (), "declared", setup=block_locus_table)
    sl = E.KV(new + ".sequence_locus.tsv")
    check(rc == 2 and sl.get("exit_code") == "2" and sl.get("locus_output", "").startswith("FAILED: "),
          "reporting failure: exit 2, recorded in the primary record beside the prefix")
    so = E.KV(old + ".sequence_locus.tsv")
    strip = lambda m: {k: v.replace(old, "<RUN>").replace(new, "<RUN>").replace(f"{new}.diagnostics/c4", "<RUN>")
                       for k, v in m.items() if "seconds" not in k}
    check(strip(sl) == strip(so), "reporting failure: the primary record equals M15-K's (its paths aside)")
    # (3) a refused run (a declared library that does not exist): exit 2, nothing created
    ref = f"{d}/refused/c4"
    os.makedirs(os.path.dirname(ref), exist_ok=True)
    p = E.run(E.cmd("c4", ref, f"{d}/no_such_library.txt"))
    check(p.returncode == 2 and os.listdir(os.path.dirname(ref)) == [],
          f"refused: exit 2 and nothing created, not even the diagnostics directory ({os.listdir(os.path.dirname(ref))})")
    return 0


def ex_stale_flat(d):
    """A prefix holding an earlier run's flat reports: those report names are removed; an unrelated file is kept."""
    out = f"{d}/c4"
    os.makedirs(d, exist_ok=True)
    stale = [".mosaic_status.tsv", ".block_genotypes.tsv", ".locus_genotype.tsv", ".settings.tsv", ".mosaic_haplotypes.fa"]
    for s in stale:
        open(out + s, "w").write("an earlier flat run\n")
    open(out + ".notes.txt", "w").write("not a run report\n")
    open(f"{d}/unrelated.tsv", "w").write("not a run report\n")
    p = E.run(E.cmd("c4", out))
    check(p.returncode == 0, "stale flat: exit 0")
    check(not any(os.path.exists(out + s) for s in stale), "stale flat: the earlier flat reports are removed")
    check(os.path.exists(out + ".notes.txt") and os.path.exists(f"{d}/unrelated.tsv"), "stale flat: unrelated files are kept")
    check(os.path.exists(f"{out}.diagnostics/c4.mosaic_status.tsv") and os.path.exists(out + ".sequence_locus.tsv"),
          "stale flat: the run's own reports are in the new layout")
    return 0


if __name__ == "__main__":
    name, d = sys.argv[1], sys.argv[2]
    os.makedirs(d, exist_ok=True)
    globals()["ex_" + name](d)
    bad = [w for ok, w in CHECKS if not ok]
    print(f"{len(CHECKS) - len(bad)} of {len(CHECKS)} checks passed" + (f"; FAILED: {bad}" if bad else ""))
    sys.exit(1 if bad else 0)
