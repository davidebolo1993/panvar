#!/usr/bin/env python3
"""Build and validate the C4 genotype report from inputs made by build_inputs.py (no genotyping). Writes, under
--work-dir: report/ (the report), the gate runs report_rerun/, report_redraw/, report_modality/ and report_qc/,
reconciliation.tsv, run_record.tsv (provenance) and validation.txt. Every output directory must be new.

  V1 inputs: every seal verifies before and after; every input file the report reads has the same SHA-256 afterwards
  V2 determinism: a second run gives byte-identical data/*.tsv, figures and report.md; the manifests differ only in the
     command line (the --out-dir)
  V3 reconciliation: reconcile.py (campaign 484, replay 649, headline audit, phase evaluation)
  V4 every plotted value is recoverable from data/: every figure redrawn from the TSVs read back from disk (the
     in-memory tables are discarded) is byte-identical to the original
  V5 modality: run metadata with modality=short-read for every run gives the same aggregates (plus the modality column)
     and byte-identical figures
  V6 truth-free mode reads no truth input, draws no accuracy figure and shows no block-length information
  V7 (with --reference-report) the same scientific values as an earlier report: every data/ table of the reference
     equal, after applying --rename OLD=NEW column renames to it; the only tables it lacks are the --new-table ones
  V8 block_lengths.tsv: recomputed independently from the evaluator tables the manifest lists (every truth1_bp and
     truth2_bp of a row with both defined), and its minimum and maximum equal the C4 ranges in C4_TRUTH_ALLELE_RANGES
Exit 1 on any failure.
"""
import os
import subprocess
import sys

import numpy as np
import pandas as pd
import pandas.api.types as ptypes

import c4_common as C

lines = []

# The observed C4 truth-allele length ranges (bp) of the evaluated dataset, per block: the expected values the report's
# block labels must show. C4-specific, so they live here and never in the generic tool.
C4_TRUTH_ALLELE_RANGES = {"ctx:0": (12306, 12402), "sv:1": (0, 7670), "ctx:2": (20576, 20927), "sv:3": (117, 420),
                          "ctx:4": (325, 375), "sv:5": (0, 117), "ctx:6": (47951, 48003), "sv:7": (26346, 98187),
                          "ctx:8": (14267, 14431), "sv:9": (0, 69), "ctx:10": (63396, 63467)}


def report(ok, what):
    lines.append(("PASS " if ok else "FAIL ") + what)
    print(lines[-1], flush=True)
    return ok


def tree(d, sub):
    base = os.path.join(d, sub)
    return {f: C.sha256(os.path.join(base, f)) for f in sorted(os.listdir(base))}


def check_block_lengths(out, inputs):
    """V8: the table against the evaluator tables themselves (read as text, parsed with float), and against the
    expected C4 ranges."""
    man = pd.read_csv(os.path.join(inputs, "evaluation_manifest.tsv"), sep="\t", dtype=str, keep_default_na=False)
    vals = {}
    for p in man.path:
        e = pd.read_csv(os.path.join(inputs, p), sep="\t", dtype=str, keep_default_na=False)
        for r in e.itertuples():
            if r.truth1_bp not in ("NA", "") and r.truth2_bp not in ("NA", ""):
                vals.setdefault(r.block_id, []).extend([float(r.truth1_bp), float(r.truth2_bp)])
    t = pd.read_csv(os.path.join(out, "data", "block_lengths.tsv"), sep="\t", dtype=str, keep_default_na=False).set_index("block_id")
    same = sorted(t.index) == sorted(vals) and all(
        int(t.loc[b, "truth_allele_observations"]) == len(v) and float(t.loc[b, "truth_allele_bp_min"]) == min(v)
        and float(t.loc[b, "truth_allele_bp_max"]) == max(v) and float(t.loc[b, "truth_allele_bp_median"]) == float(np.median(v))
        for b, v in vals.items())
    ok = report(same, f"V8 block_lengths.tsv ({len(t)} blocks, {sum(len(v) for v in vals.values())} truth-allele "
                      f"observations) equals an independent recomputation from the {len(man)} evaluator tables")
    got = {b: (float(t.loc[b, "truth_allele_bp_min"]), float(t.loc[b, "truth_allele_bp_max"])) for b in t.index}
    bad = [b for b, r in C4_TRUTH_ALLELE_RANGES.items() if got.get(b) != r]
    ok &= report(not bad and len(got) == len(C4_TRUTH_ALLELE_RANGES),
                 f"V8 the minimum and maximum of all {len(C4_TRUTH_ALLELE_RANGES)} blocks equal the expected C4 ranges; "
                 f"differing: {bad or 'none'}; labels: " + ", ".join(f"{b} {t.loc[b, 'displayed_length_range']}" for b in t.index))
    return ok


def main():
    p = C.parser(__doc__.split("\n")[0])
    p.add_argument("--reference-report", help="an earlier report directory whose data/ values must be reproduced")
    p.add_argument("--rename", action="append", default=[], help="OLD=NEW: a column renamed since the reference report")
    p.add_argument("--new-table", action="append", default=[], help="NAME: a data table added since the reference report")
    args = p.parse_args()
    repo, results = C.roots(args)
    sys.path.insert(0, os.path.join(repo, "scripts"))
    import plot_genotype_results as cli
    from genotype_report import pipeline

    work = os.path.abspath(args.work_dir)
    inputs = os.path.join(work, "inputs")
    if not os.path.isdir(inputs):
        sys.exit(f"{inputs} not found: run build_inputs.py first")
    outs = {k: os.path.join(work, k) for k in ("report", "report_rerun", "report_redraw", "report_modality", "report_qc")}
    for p in list(outs.values()) + [os.path.join(work, "inputs_modality")]:
        if os.path.exists(p):
            sys.exit(f"{p} exists: every output directory must be new")

    def run(a, cls=pipeline.Report):
        r = cls(cli.parse(a), ["plot_genotype_results.py"] + a)
        r.load()
        os.makedirs(r.out, exist_ok=True)
        r.aggregate()
        r.draw()
        r.write_report()
        r.write_manifest()
        return r

    class FromDisk(pipeline.Report):
        """Aggregates as usual (writing data/), then discards every in-memory table and draws from the TSVs read back
        from disk, typed only by each column's schema (bool, integer, nullable integer, float, text)."""
        def aggregate(self):
            super().aggregate()
            for name, mem in list(self.tables.items()):
                t = pd.read_csv(os.path.join(self.out, "data", name + ".tsv"), sep="\t", dtype=str, keep_default_na=False)
                for c in t.columns:
                    kind = mem[c].dtype
                    if ptypes.is_bool_dtype(kind):
                        t[c] = t[c].map({"True": True, "False": False})
                    elif str(kind) == "Int64":
                        t[c] = t[c].map(lambda v: None if v == "NA" else int(v)).astype("Int64")
                    elif ptypes.is_integer_dtype(kind):
                        t[c] = t[c].astype(kind)
                    elif ptypes.is_float_dtype(kind):
                        t[c] = t[c].map(lambda v: float("nan") if v == "NA" else float(v)).astype(float)   # exact
                    else:
                        t[c] = t[c].astype(kind)
                        if mem[c].isna().any():
                            t[c] = t[c].where(t[c] != "NA", np.nan)
                        if len(mem[c].dropna()) and set(mem[c].dropna().map(type)) == {bool}:
                            t[c] = t[c].map(lambda v: {"True": True, "False": False}.get(v, v) if isinstance(v, str) else v)
                self.tables[name] = t

    C.write_record(os.path.join(work, "run_record.tsv"), C.provenance(repo, results))
    ok = report(C.verify_seals(results), "V1 every seal verifies before the runs")
    a = run(C.tool_argv(results, inputs, outs["report"]))
    inputs_before = {p: C.sha256(p) for p in sorted(set(os.path.abspath(s) for s in a.sources))}
    b = run(C.tool_argv(results, inputs, outs["report_rerun"]))
    for sub in ("data", "figures"):
        ta, tb = tree(a.out, sub), tree(b.out, sub)
        ok &= report(ta == tb, f"V2 {sub}/: {len(ta)} files byte-identical across two runs")
    ok &= report(C.sha256(os.path.join(a.out, "report.md")) == C.sha256(os.path.join(b.out, "report.md")), "V2 report.md byte-identical")
    ma = pd.read_csv(os.path.join(a.out, "manifest.tsv"), sep="\t", dtype=str, keep_default_na=False)
    mb = pd.read_csv(os.path.join(b.out, "manifest.tsv"), sep="\t", dtype=str, keep_default_na=False)
    diff = [r.name for r, s in zip(ma.itertuples(index=False), mb.itertuples(index=False)) if tuple(r) != tuple(s)]
    ok &= report(len(ma) == len(mb) and diff == ["command"], f"V2 manifest.tsv ({len(ma)} rows) differs only in: {diff}")
    rc = subprocess.run([sys.executable, os.path.join(C.HERE, "reconcile.py"), "--report", a.out, "--repo-root", repo,
                         "--results-root", results, "--out", os.path.join(work, "reconciliation.tsv")],
                        capture_output=True, text=True)
    ok &= report(rc.returncode == 0, "V3 reconciliation: " + (rc.stdout.strip().splitlines() or ["no output: " + rc.stderr[-300:]])[-1])
    c = run(C.tool_argv(results, inputs, outs["report_redraw"]), FromDisk)
    ta = tree(a.out, "figures")
    ok &= report(ta == tree(c.out, "figures"), f"V4 all {len(ta)} figure files redrawn from data/*.tsv read back from disk are byte-identical")
    data_backed = all(os.path.isfile(os.path.join(a.out, "data", d + ".tsv")) for _, _, _, ds in a.figures for d in ds)
    ok &= report(data_backed and len(a.figures) > 0, f"V4 every one of {len(a.figures)} figures names backing TSVs that exist in data/")
    os.makedirs(os.path.join(work, "inputs_modality"))
    meta = pd.read_csv(os.path.join(inputs, "run_metadata.tsv"), sep="\t", dtype=str, keep_default_na=False)
    mpath = os.path.join(work, "inputs_modality", "run_metadata.tsv")
    meta.assign(modality="short-read").to_csv(mpath, sep="\t", index=False, lineterminator="\n")
    m = run(C.tool_argv(results, inputs, outs["report_modality"], mpath))
    same, carried = True, True
    for f in sorted(os.listdir(os.path.join(a.out, "data"))):
        x = pd.read_csv(os.path.join(a.out, "data", f), sep="\t", dtype=str, keep_default_na=False)
        y = pd.read_csv(os.path.join(m.out, "data", f), sep="\t", dtype=str, keep_default_na=False)
        if "modality" not in y.columns and f != "switch_map.tsv":
            carried = False
            print("no modality column in", f)
        y = y.drop(columns=[c for c in ("modality",) if c in y.columns])
        if not x.equals(y):
            same = False
            print("differs apart from modality:", f)
    ok &= report(carried, "V5 every aggregate carries the modality column (switch_map.tsv is the supplied plan description, "
                          "not an aggregate of runs)")
    ok &= report(same, "V5 every aggregate equals the no-modality run once the modality column is dropped")
    ok &= report(tree(a.out, "figures") == tree(m.out, "figures"), "V5 every figure byte-identical with a single modality")
    q = run(C.tool_argv(results, inputs, outs["report_qc"], truth=False))
    names = [n for _, n, _, _ in q.figures]
    acc = [n for n in names if n.startswith(("eval_", "mosaic_block", "mosaic_excess", "dosage_content", "phase_accuracy"))]
    qin = pd.read_csv(os.path.join(q.out, "manifest.tsv"), sep="\t", dtype=str, keep_default_na=False)
    truth_read = [v for k, v in zip(qin.kind, qin.value) if k == "input" and (C.HEADLINE in v or C.PHASE in v or "eval" in os.path.basename(v)
                                                                               or "switch" in os.path.basename(v))]
    with open(os.path.join(q.out, "report.md")) as fh:
        truth_free = "truth-free QC report" in fh.read()
    ok &= report(not acc and not truth_read and truth_free and not os.path.exists(os.path.join(q.out, "data", "block_evaluation.tsv")),
                 f"V6 truth-free run: {len(names)} figures, no accuracy figure, no truth input read "
                 f"({len(qin[qin.kind == 'input'])} inputs, all package or run metadata), {len(q.skipped)} skipped with reasons")
    with open(os.path.join(q.out, "report.md")) as fh:
        qtext = fh.read()
    ok &= report(not os.path.exists(os.path.join(q.out, "data", "block_lengths.tsv")) and "truth-allele" not in qtext
                 and not any("block_lengths" in ds for _, _, _, ds in q.figures),
                 "V6 truth-free run: no block_lengths.tsv, no figure drawn from it, no block-length text in report.md")
    if args.reference_report:
        ref = os.path.join(os.path.abspath(args.reference_report), "data")
        renames = dict(r.split("=", 1) for r in args.rename)
        new = sorted(n + ".tsv" for n in args.new_table)
        fr, fa = sorted(os.listdir(ref)), sorted(os.listdir(os.path.join(a.out, "data")))
        differ = []
        for f in fr:
            x = pd.read_csv(os.path.join(ref, f), sep="\t", dtype=str, keep_default_na=False).rename(columns=renames)
            y = pd.read_csv(os.path.join(a.out, "data", f), sep="\t", dtype=str, keep_default_na=False) if f in fa else None
            if y is None or not x.equals(y):
                differ.append(f)
        added = sorted(set(fa) - set(fr))
        ok &= report(not differ and added == new, f"V7 all {len(fr)} data tables of the reference report "
                     f"{args.reference_report} are present and equal after renames {renames or 'none'} (differing: "
                     f"{differ or 'none'}); tables added: {added or 'none'}, declared new: {new or 'none'}")
    ok &= check_block_lengths(a.out, inputs)
    after = {p: C.sha256(p) for p in inputs_before}
    ok &= report(after == inputs_before, f"V1 all {len(after)} input files unchanged after five runs")
    ok &= report(C.verify_seals(results), "V1 every seal verifies after the runs")
    with open(os.path.join(work, "validation.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
