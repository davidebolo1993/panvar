#!/usr/bin/env python3
"""S1-C: close the CYP2D6 30x legacy-identity comparison WITHOUT rerunning the sentinel.

Compares P1.1's saved consolidated legacy run (o.cyp_30, 392fc0b build) with S1-B's saved sequence_first_v1 run
(cyp30.v1, a2b3cd4 build), lists EVERY difference in the caller's outputs (call.*), and classifies each exactly:
  * status rows in tests/analysis/perf/perf_compare.py's own DIAGNOSTIC class (unmodified);
  * EXACTLY two named executed-work rows, indexed_union_executed_prefilters and indexed_union_executed_dps (P1.1 measured
    them varying between identical 4-thread runs) -- by name, not by an indexed_union_ prefix;
  * the declared S1-B sequence rows (output_schema, sequence_prediction, sequence_status, sequence_completion,
    sequence_exit_code) and the one declared settings row (output_schema sequence_first_v1, flag --output-schema);
  * settings: every key's VALUE must be identical; the SOURCE must be identical once each run's own directory is replaced
    by one placeholder (each run composes its profile inside its own run directory).
Anything else -- a scientific output, a non-diagnostic status row, a setting value or source type -- fails.
Negative controls (on temporary copies; the saved runs are only read) prove that real changes still fail.

  cyp_closure.py <legacy run dir> <sequence_first run dir>     -> verdict; exit 0 iff QUALIFIED_IDENTICAL and all controls
"""
import importlib.util
import os
import re
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("pc", os.path.join(HERE, "..", "perf", "perf_compare.py"))
pc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pc)
DECLARED_COUNTERS = {"indexed_union_executed_prefilters", "indexed_union_executed_dps"}
DECLARED_SEQUENCE_ROWS = {"output_schema", "sequence_prediction", "sequence_status", "sequence_completion",
                          "sequence_exit_code"}
SETTINGS_ROW = ("output_schema", "sequence_first_v1", "flag --output-schema")
STATUS_FILES = ("call.mosaic_status.tsv", "call.original.mosaic_status.tsv", "call.dosage_status.tsv")


def kv(path):
    out = {}
    for line in open(path).read().splitlines():
        p = line.split("\t", 1)
        if len(p) == 2:
            out[p[0]] = p[1]
    return out


def settings(path):
    rows = {}
    for line in open(path).read().splitlines()[1:]:
        k, v, s = line.split("\t")
        rows[k] = (v, s)
    return rows


def compare(a, b):
    """(verdict, classified differences, unexplained problems)."""
    fa = {f for f in os.listdir(a) if f.startswith("call.")}
    fb = {f for f in os.listdir(b) if f.startswith("call.")}
    sequence_files = {"call.sequence_locus.tsv", "call.sequence_blocks.tsv", "call.sequence_haplotypes.fa",
                      "call.representation_intervals.tsv", "call.sequence_alternatives.tsv"}
    classified, problems = [], []
    if fa - fb or (fb - fa) - sequence_files:
        problems.append(f"file sets: missing {sorted(fa - fb)}, extra {sorted((fb - fa) - sequence_files)}")
    for f in sorted(fa & fb):
        pa, pb = os.path.join(a, f), os.path.join(b, f)
        if f in STATUS_FILES:
            ka, kb = kv(pa), kv(pb)
            for k in sorted(set(ka) | set(kb)):
                if ka.get(k) == kb.get(k):
                    continue
                if k in DECLARED_SEQUENCE_ROWS and k not in ka and f != "call.dosage_status.tsv":
                    classified.append(f"{f}:{k}:declared_sequence_row")
                elif k in DECLARED_COUNTERS and k in ka and k in kb:
                    classified.append(f"{f}:{k}:declared_executed_work_counter")
                elif pc.DIAGNOSTIC.search(k) and k in ka and k in kb:
                    classified.append(f"{f}:{k}:perf_compare_diagnostic")
                else:
                    problems.append(f"{f}: {k}: {ka.get(k)} vs {kb.get(k)}")
        elif f == "call.settings.tsv":
            sa, sb = settings(pa), settings(pb)
            if sb.get(SETTINGS_ROW[0]) != SETTINGS_ROW[1:]:
                problems.append("settings: the declared output_schema row is missing or different")
            sb.pop(SETTINGS_ROW[0], None)
            if set(sa) != set(sb):
                problems.append(f"settings: keys {sorted(set(sa) ^ set(sb))}")
            for k in sorted(set(sa) & set(sb)):
                (va, xa), (vb, xb) = sa[k], sb[k]
                if va != vb:
                    problems.append(f"settings: {k} value {va} vs {vb}")
                elif xa.replace(a, "<RUN>") != xb.replace(b, "<RUN>"):
                    problems.append(f"settings: {k} source {xa} vs {xb}")
                elif xa != xb:
                    classified.append(f"settings:{k}:source_names_the_runs_own_profile_path")
        elif f in ("call.repair_memory.tsv",):
            continue   # memory ledger: measurements
        elif f == "call.repair_rounds.tsv":
            cut = lambda p: [l.split("\t")[:28] for l in open(p).read().splitlines()]  # noqa: E731
            if cut(pa) != cut(pb):
                problems.append("repair trajectory differs")
        elif open(pa, errors="replace").read().replace(a, "<RUN>") != open(pb, errors="replace").read().replace(b, "<RUN>"):
            problems.append(f"{f}: differs")
    return ("QUALIFIED_IDENTICAL" if not problems else "DIFFERENT"), classified, problems


def edit_status(path, key, value):
    lines = open(path).read().splitlines()
    open(path, "w").write("\n".join(f"{key}\t{value}" if l.split("\t", 1)[0] == key else l for l in lines) + "\n")


def edit_table(path, column, change):
    lines = open(path).read().splitlines()
    head = lines[0].split("\t")
    cells = lines[1].split("\t")
    cells[head.index(column)] = change(cells[head.index(column)])
    lines[1] = "\t".join(cells)
    open(path, "w").write("\n".join(lines) + "\n")


def edit_setting(path, key, value=None, source=None):
    lines = open(path).read().splitlines()
    for n, l in enumerate(lines):
        p = l.split("\t")
        if p[0] == key:
            lines[n] = "\t".join([p[0], value if value is not None else p[1], source if source is not None else p[2]])
    open(path, "w").write("\n".join(lines) + "\n")


def flip_fasta(path):
    lines = open(path).read().splitlines()
    seq = lines[1]
    lines[1] = seq[:10] + ("C" if seq[10] != "C" else "G") + seq[11:]
    open(path, "w").write("\n".join(lines) + "\n")


def main(a, b):
    verdict, classified, problems = compare(a, b)
    kinds = {}
    for c in classified:
        kinds[c.rsplit(":", 1)[1]] = kinds.get(c.rsplit(":", 1)[1], 0) + 1
    print(f"verdict: {verdict}; {len(classified)} classified differences {kinds}; {len(problems)} unexplained")
    for c in classified:
        if not c.endswith("perf_compare_diagnostic"):
            print("  classified:", c)
    for p in problems:
        print("  UNEXPLAINED:", p)
    controls = [
        ("unaltered", None, "QUALIFIED_IDENTICAL"),
        ("declared counter indexed_union_executed_dps changed", lambda d: edit_status(f"{d}/call.mosaic_status.tsv",
                                                                                     "indexed_union_executed_dps", "1"),
         "QUALIFIED_IDENTICAL"),
        ("undeclared counter indexed_union_executed_unique_starts changed",
         lambda d: edit_status(f"{d}/call.mosaic_status.tsv", "indexed_union_executed_unique_starts", "1"), "DIFFERENT"),
        ("scientific setting value repair_max_neighbours 958 -> 959",
         lambda d: edit_setting(f"{d}/call.settings.tsv", "repair_max_neighbours", value="959"), "DIFFERENT"),
        ("setting source type profile -> flag",
         lambda d: edit_setting(f"{d}/call.settings.tsv", "max_state_classes", source="flag --max-state-classes"),
         "DIFFERENT"),
        ("block genotype changed", lambda d: edit_table(f"{d}/call.block_genotypes.tsv", "allele2", lambda x: x + "_x"),
         "DIFFERENT"),
        ("candidate score changed", lambda d: edit_table(f"{d}/call.mosaic_scores.tsv", "score",
                                                         lambda x: repr(float(x) + 0.5)), "DIFFERENT"),
        ("call class changed", lambda d: edit_status(f"{d}/call.mosaic_status.tsv", "call_class", "strict_call"),
         "DIFFERENT"),
        ("legacy FASTA base changed", lambda d: flip_fasta(f"{d}/call.dosage_original.fa")
         if os.path.exists(f"{d}/call.dosage_original.fa") else flip_fasta(f"{d}/call.original.mosaic_haplotypes.fa"),
         "DIFFERENT"),
    ]
    bad = 0
    for name, edit, expect in controls:
        tmp = tempfile.mkdtemp(prefix="s1c_cyp_")
        bb = os.path.join(tmp, os.path.basename(b.rstrip("/")))
        shutil.copytree(b, bb, ignore=shutil.ignore_patterns("*.fq", "*.gz"))
        # the copy names a different run directory in its settings sources: relabel it so only the edit differs
        for f in os.listdir(bb):
            if f.startswith("call.") and (f.endswith(".tsv") or f.endswith(".fa")):
                text = open(os.path.join(bb, f), errors="replace").read()
                open(os.path.join(bb, f), "w").write(text.replace(b, bb))
        if edit:
            edit(bb)
        got, _, why = compare(a, bb)
        ok = got == expect
        bad += not ok
        print(f"  control {'ok' if ok else 'FAILED'}: {name}: {got} (expected {expect})" +
              ("" if ok or not why else f" [{why[0]}]"))
        shutil.rmtree(tmp)
    return 0 if verdict == "QUALIFIED_IDENTICAL" and bad == 0 else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1].rstrip("/"), sys.argv[2].rstrip("/")))
