#!/usr/bin/env python3
"""S1-D.1 fixture gates: the sequence_first_v2 FASTA policy across the whole output prefix, and the dosage contract at
the CLI (the direct API and prepared routes are gated in genotype_command_e2e, test_sequence_first_v2_contract).

Reuses the S1-D fixtures (s1d_fixtures.py builders; F3, F9, F1 panel/panel, F11c panel/off-panel, F7full off/off, F7
provisional, F11u, F12 reverse). Per v2 case (CLI, no dosage flag):
  - membership (independent, from the GFA), extraction --check and the FASTA decision, as in S1-D;
  - the COMPLETE FASTA inventory of the prefix: exactly {call.sequence_haplotypes.fa} for an off-panel pair, nothing for
    a panel pair; a non-empty dosage FASTA would be an audit diagnostic and is listed separately (none expected here);
  - no output references a missing prefix file (call.<name>.tsv|fa, also inside paths and before '#');
  - the locus table links only the authoritative FASTA (main_fasta) and names no suppressed record (hap ids NA);
  - dosage: default audit with source contract_default_sequence_first_v2; the audit status is recorded as found
    (NOT_RUN on a legacy no-call is the disclosed limitation, never counted as a pass).
F11c also runs the CLI dosage arms: --dosage-mode audit / off, --dosage-audit off, and --dosage-mode repair (refused,
exit 2, a sentinel prefix left byte-identical). F11u: legacy identity (unset schema, --dosage-mode audit) with only the
declared differences (no sequence reports; v2 lacks the legacy FASTA and the empty dosage FASTAs; the locus table
differs only in main_fasta, hap1_id, hap2_id, dosage_representative_fasta), and the loaded persisted index (write-panel, then loaded) equal to cold.
Lifecycle on one prefix seeded with stale run-owned FASTAs: off-panel -> panel-only -> off-panel -> watchdog ->
off-panel -> refused repair. Cross-library non-regression: earlier legacy and sequence_first_v1 runs (S1-C library
c744035: F1 base and F7 base, with --dosage-mode repair; S1-D library c947bfc: F11u legacy, audit) rerun with this
library and compared byte for byte (paths normalised; status timing keys and repair_memory excluded).
  s1d1_fixtures.py <out dir> <panvar> <genotype_prepared_locus_driver> [section,...]
  sections: cases, arms, legacy, lifecycle, crosslib
  s1d1_fixtures.py --check-prefix <run dir> <expected fasta_reason> <graph.gfa>   (the real-data check)
"""
import csv
import os
import re
import shlex
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "../s1d"))
import s1d_fixtures as s1d  # noqa: E402  (also puts fixtures, pilot, s1a_fixtures, extract_panel_haplotype on sys.path)

fixtures, s1a_fixtures = s1d.fixtures, s1d.s1a_fixtures
REF = s1d.REF
RESULTS = "/Users/davide.bolognini/panvar/tests/results"
SUPPRESSED = {"call.mosaic_haplotypes.fa", "call.dosage_original.fa", "call.dosage_repaired.fa",
              "call.dosage_repaired_block_alleles.fa"}
REF_PATTERN = re.compile(r"call\.[A-Za-z0-9_.]+?\.(?:tsv|fa|fasta)\b(?!\.writing)")   # a staging name is no reference
VOLATILE = ("seconds", "footprint", "rss", "peak", "_bytes", "worker", "threads_", "profile_sha256", "timing")


def fasta_inventory(d):
    core, diagnostic = set(), set()
    for f in os.listdir(d):
        if f.startswith("call.") and (f.endswith(".fa") or f.endswith(".fasta")):
            if f.startswith("call.dosage_") and os.path.getsize(os.path.join(d, f)) > 0:
                diagnostic.add(f)
            else:
                core.add(f)
    return core, diagnostic


def dangling(d):
    out = []
    for f in sorted(os.listdir(d)):
        p = os.path.join(d, f)
        if not f.startswith("call.") or not os.path.isfile(p):
            continue
        for name in REF_PATTERN.findall(open(p, errors="replace").read()):
            if not os.path.exists(os.path.join(d, name)):
                out.append(f"{f}->{name}")
    return sorted(set(out))


def locus_row(path):
    rows = list(csv.DictReader(open(path), delimiter="\t")) if os.path.exists(path) else []
    return rows[0] if rows else {}


def case_args(out, kind, profile, curation):
    """Builds the S1-D fixture graph and reads under <out>/<kind> and returns (gfa, base args, reads dir)."""
    d = os.path.join(out, kind)
    os.makedirs(d, exist_ok=True)
    gfa = os.path.join(d, "graph.gfa")
    if kind in ("F1", "F3"):
        nodes, rows = fixtures.build(kind)
        fixtures.write_gfa(gfa, nodes, rows)
        reads = os.path.join(s1d.OLD, kind, "reads")
    elif kind == "F9":
        import s1b_fixtures
        nodes, rows, donor, gap, _ = s1b_fixtures.build("F9")
        fixtures.write_gfa(gfa, nodes, rows)
        reads = os.path.join(d, "reads")
        s1a_fixtures.tiled_reads(nodes, donor, gap, reads)
    elif kind.startswith("F7"):
        nodes, rows, donor, gap = s1a_fixtures.build("F7")
        fixtures.write_gfa(gfa, nodes, rows)
        reads = os.path.join(d, "reads")
        s1a_fixtures.tiled_reads(nodes, donor, gap, reads)
    else:
        nodes, rows, donor, reversed_paths = s1d.f11(reverse_s21=(kind == "F12"))
        s1d.write_gfa(gfa, nodes, rows, reversed_paths)
        reads = os.path.join(d, "reads")
        s1a_fixtures.tiled_reads(nodes, donor, None, reads)
    args = ["--preset", "short-read", "--profile", os.path.join(out, f"profile_{profile}.tsv"), "--threads", "2",
            "--graph", gfa, "--reference-path", REF]
    if curation:
        args += ["--panel-curation", os.path.join(out, curation)]
    args += ["--reads1", os.path.join(reads, "reads_1.fq"), "--reads2", os.path.join(reads, "reads_2.fq"), "--library",
             os.path.join(reads, "library.txt"), "--quiet"]
    return d, gfa, args


def v2_checks(prefix_dir, expect_reason):
    """The whole-prefix checks of one v2 run; returns (fields, problems)."""
    loc = s1d.kv(os.path.join(prefix_dir, "call.sequence_locus.tsv"))
    core, diagnostic = fasta_inventory(prefix_dir)
    off_panel = expect_reason == "OFF_PANEL_PAIR"
    want = {"call.sequence_haplotypes.fa"} if off_panel else set()
    problems = []
    if loc.get("fasta_reason") != expect_reason:
        problems.append(f"fasta_reason {loc.get('fasta_reason')}")
    if core != want:
        problems.append(f"inventory {sorted(core)} != {sorted(want)}")
    refs = dangling(prefix_dir)
    if refs:
        problems.append("dangling " + ",".join(refs[:3]))
    row = locus_row(os.path.join(prefix_dir, "call.locus_genotype.tsv"))
    if row.get("main_fasta") != ("call.sequence_haplotypes.fa" if off_panel else "NA") or row.get("hap1_id") != "NA" or \
            row.get("hap2_id") != "NA" or row.get("original_fasta") != "NA":
        problems.append(f"locus links main_fasta={row.get('main_fasta')} hap1_id={row.get('hap1_id')}")
    status = s1d.kv(os.path.join(prefix_dir, "call.dosage_status.tsv"))
    fields = {"fasta_inventory": ",".join(sorted(core)) or "none",
              "dosage_diagnostic_fasta": ",".join(sorted(diagnostic)) or "none",
              "dangling_references": len(refs), "main_fasta": row.get("main_fasta"),
              "dosage_mode": loc.get("dosage_mode"), "dosage_mode_source": loc.get("dosage_mode_source"),
              "dosage_audit_status": status.get("dosage_audit_status", "absent"),
              "locus_dosage_audit_status": row.get("dosage_audit_status")}
    return loc, fields, problems


def section_cases(out, panvar, driver, rows):
    plan = [("F3", "default", None, "PANEL_REPRESENTED_PAIR"), ("F9", "k0", None, "PANEL_REPRESENTED_PAIR"),
            ("F1", "default", None, "PANEL_REPRESENTED_PAIR"), ("F7full", "k0", None, "OFF_PANEL_PAIR"),
            ("F7", "round1", None, "OFF_PANEL_PAIR"), ("F11u", "k0", None, "PANEL_REPRESENTED_PAIR"),
            ("F11c", "k0", "curation_F11.txt", "OFF_PANEL_PAIR"), ("F12", "k0", None, "PANEL_REPRESENTED_PAIR")]
    for kind, profile, curation, expect in plan:
        d, gfa, args = case_args(out, kind, profile, curation)
        v2 = os.path.join(d, "v2")
        os.makedirs(v2, exist_ok=True)
        code = s1d.run([panvar, "genotype", *args, "--output-schema", "sequence_first_v2", "--out", f"{v2}/call"],
                       os.path.join(d, "v2.log"))
        loc, fields, problems = v2_checks(v2, expect)
        _, member = s1d.membership_check(gfa, f"{v2}/call.sequence_haplotypes.tsv", {"S2#1#c"} if curation else set())
        check = subprocess.run([sys.executable, os.path.join(HERE, "../s1d/extract_panel_haplotype.py"), "--check", gfa,
                                f"{v2}/call.sequence_haplotypes.tsv"], capture_output=True, text=True)
        extraction = check.stdout.strip().splitlines()[-1] if check.stdout.strip() else "no output"
        problems += member
        if check.returncode != 0:
            problems.append("extraction FAILED")
        if fields["dosage_mode"] != "audit" or fields["dosage_mode_source"] != "contract_default_sequence_first_v2":
            problems.append("dosage default not audit")
        legacy_class = loc.get("legacy_call_class")
        audit_note = ("audit ran" if fields["dosage_audit_status"] == "COMPLETE" else
                      "NOT_RUN on a legacy no-call (disclosed limitation; not validated)" if legacy_class == "no_call" and
                      fields["dosage_audit_status"] == "NOT_RUN" else "UNEXPECTED")
        if audit_note == "UNEXPECTED":
            problems.append(f"audit status {fields['dosage_audit_status']} with legacy class {legacy_class}")
        rows.append({"section": "cases", "case": kind, "exit": code, "prediction": loc.get("prediction"),
                     "hap_origins": f"{loc.get('hap1_origin')}/{loc.get('hap2_origin')}", "fasta_reason": loc.get("fasta_reason"),
                     **fields, "legacy_call_class": legacy_class, "audit_note": audit_note, "extraction": extraction,
                     "result": "PASS" if not problems else "FAIL: " + " | ".join(problems)})
        print("\t".join(str(v) for v in rows[-1].values()), flush=True)


def section_arms(out, panvar, rows):
    """CLI dosage arms on F11c (an off-panel pair: the sequence FASTA is the only FASTA in every accepted arm)."""
    d, _, args = case_args(out, "F11c", "k0", "curation_F11.txt")
    arms = [("audit_flag", ["--dosage-mode", "audit"], "audit", "flag --dosage-mode", True),
            ("off_flag", ["--dosage-mode", "off"], "off", "flag --dosage-mode", False),
            ("dosage_audit_off", ["--dosage-audit", "off"], "off", "flag --dosage-audit", False)]
    for label, extra, mode, source, audited in arms:
        p = os.path.join(d, "arm_" + label)
        os.makedirs(p, exist_ok=True)
        code = s1d.run([panvar, "genotype", *args, *extra, "--output-schema", "sequence_first_v2", "--out", f"{p}/call"],
                       os.path.join(d, f"arm_{label}.log"))
        loc, fields, problems = v2_checks(p, "OFF_PANEL_PAIR")
        if fields["dosage_mode"] != mode or fields["dosage_mode_source"] != source:
            problems.append(f"recorded {fields['dosage_mode']}/{fields['dosage_mode_source']}")
        has_audit = os.path.exists(f"{p}/call.dosage_status.tsv") and os.path.exists(f"{p}/call.dosage_blocks.tsv")
        if has_audit != audited or (audited and fields["dosage_audit_status"] != "COMPLETE"):
            problems.append(f"audit artefacts present={has_audit}")
        row = locus_row(f"{p}/call.locus_genotype.tsv")
        if row.get("dosage_mode") != mode or (fields["locus_dosage_audit_status"] == "COMPLETE") != audited:
            problems.append(f"locus dosage {row.get('dosage_mode')}/{fields['locus_dosage_audit_status']}")
        rows.append({"section": "arms", "case": "F11c:" + label, "exit": code, "prediction": loc.get("prediction"),
                     "hap_origins": f"{loc.get('hap1_origin')}/{loc.get('hap2_origin')}", "fasta_reason": loc.get("fasta_reason"),
                     **fields, "legacy_call_class": loc.get("legacy_call_class"),
                     "audit_note": "audit ran" if audited else "not evaluated (declared off)", "extraction": "-",
                     "result": "PASS" if not problems else "FAIL: " + " | ".join(problems)})
        print("\t".join(str(v) for v in rows[-1].values()), flush=True)
    # repair: refused at argument resolution, the prefix untouched
    p = os.path.join(d, "arm_repair")
    os.makedirs(p, exist_ok=True)
    sentinel = {f"{p}/call.mosaic_status.tsv": "sentinel status\n", f"{p}/call.sequence_haplotypes.fa": ">stale\nACGT\n"}
    for f, text in sentinel.items():
        open(f, "w").write(text)
    code = s1d.run([panvar, "genotype", *args, "--dosage-mode", "repair", "--output-schema", "sequence_first_v2", "--out",
                    f"{p}/call"], os.path.join(d, "arm_repair.log"))
    log = open(os.path.join(d, "arm_repair.log")).read()
    untouched = sorted(os.listdir(p)) == sorted(os.path.basename(f) for f in sentinel) and \
        all(open(f).read() == t for f, t in sentinel.items())
    ok = code == 2 and untouched and "legacy contract" in log
    rows.append({"section": "arms", "case": "F11c:repair_flag", "exit": code, "prediction": "-", "hap_origins": "-",
                 "fasta_reason": "-", "fasta_inventory": "sentinel only", "dosage_diagnostic_fasta": "-",
                 "dangling_references": "-", "main_fasta": "-", "dosage_mode": "-", "dosage_mode_source": "-",
                 "dosage_audit_status": "-", "locus_dosage_audit_status": "-", "legacy_call_class": "-",
                 "audit_note": "refused before the prefix", "extraction": "-",
                 "result": "PASS" if ok else f"FAIL: exit {code} untouched={untouched}"})
    print("\t".join(str(v) for v in rows[-1].values()), flush=True)


def compare_dirs(a, b, allowed_missing=(), locus_columns=()):
    """Byte comparison of two run directories with paths normalised; status timing keys and repair_memory excluded."""
    diffs = []
    fa = {f for f in os.listdir(a) if f.startswith("call.")}
    fb = {f for f in os.listdir(b) if f.startswith("call.")}
    if (fa - set(allowed_missing)) != fb:
        diffs.append(f"file sets {sorted((fa - set(allowed_missing)) ^ fb)}")
    for f in sorted(fa & fb):
        if f == "call.repair_memory.tsv":
            continue
        ta = open(os.path.join(a, f), errors="replace").read().replace(a, "<D>")
        tb = open(os.path.join(b, f), errors="replace").read().replace(b, "<D>")
        if f.endswith("_status.tsv"):
            ka = {k: v for k, v in s1d.kv(os.path.join(a, f)).items() if not any(x in k for x in VOLATILE)}
            kb = {k: v for k, v in s1d.kv(os.path.join(b, f)).items() if not any(x in k for x in VOLATILE)}
            ka = {k: v.replace(a, "<D>") for k, v in ka.items()}
            kb = {k: v.replace(b, "<D>") for k, v in kb.items()}
            if ka != kb:
                diffs.append(f"{f}: {sorted(set(ka.items()) ^ set(kb.items()))[:2]}")
        elif f == "call.repair_rounds.tsv":
            cut = lambda t: [line.split("\t")[:28] for line in t.splitlines()]  # noqa: E731
            if cut(ta) != cut(tb):
                diffs.append("trajectory")
        elif f == "call.locus_genotype.tsv" and locus_columns:
            ra, rb = locus_row(os.path.join(a, f)), locus_row(os.path.join(b, f))
            ra = {k: v.replace(os.path.basename(a), "") for k, v in ra.items()}
            rb = {k: v.replace(os.path.basename(b), "") for k, v in rb.items()}
            changed = sorted(k for k in set(ra) | set(rb) if ra.get(k) != rb.get(k))
            if set(changed) - set(locus_columns):
                diffs.append(f"locus columns {changed}")
        elif ta != tb:
            diffs.append(f"{f} differs")
    return diffs


def section_legacy(out, panvar, driver, rows):
    d, _, args = case_args(out, "F11u", "k0", None)
    v2 = os.path.join(d, "v2")
    if not os.path.exists(f"{v2}/call.sequence_locus.tsv"):
        os.makedirs(v2, exist_ok=True)
        s1d.run([panvar, "genotype", *args, "--output-schema", "sequence_first_v2", "--out", f"{v2}/call"],
                os.path.join(d, "v2.log"))
    leg = os.path.join(d, "legacy")
    os.makedirs(leg, exist_ok=True)
    s1d.run([panvar, "genotype", *args, "--dosage-mode", "audit", "--out", f"{leg}/call"], os.path.join(d, "legacy.log"))
    # legacy-only files allowed missing from v2: the suppressed legacy FASTA and the EMPTY dosage FASTAs
    suppressed = {f for f in SUPPRESSED if os.path.exists(f"{leg}/{f}") and
                  (f == "call.mosaic_haplotypes.fa" or os.path.getsize(f"{leg}/{f}") == 0)}
    v2_view = os.path.join(d, "v2_legacy_view")   # v2 without its sequence reports, for the comparison
    os.makedirs(v2_view, exist_ok=True)
    for f in os.listdir(v2_view):
        os.remove(os.path.join(v2_view, f))
    for f in os.listdir(v2):
        if f.startswith("call.") and ".sequence_" not in f and ".representation_intervals" not in f:
            text = open(os.path.join(v2, f), errors="replace").read().replace(v2, v2_view)
            if f == "call.settings.tsv":
                text = text.replace("output_schema\tsequence_first_v2\tflag --output-schema\n", "", 1)
            if f == "call.mosaic_status.tsv":
                keep = [line for line in text.splitlines(True) if line.split("\t", 1)[0] not in
                        {"output_schema", "sequence_prediction", "sequence_status", "sequence_completion",
                         "sequence_exit_code", "sequence_publication", "sequence_fasta", "sequence_fasta_reason"}]
                text = "".join(keep)
            open(os.path.join(v2_view, f), "w").write(text)
    diffs = compare_dirs(leg, v2_view, allowed_missing=suppressed, locus_columns=("main_fasta", "hap1_id", "hap2_id",
                                                                                 "dosage_representative_fasta"))
    leg_core, _ = fasta_inventory(leg)
    rows.append({"section": "legacy", "case": "F11u:legacy_vs_v2", "exit": "-", "prediction": "-", "hap_origins": "-",
                 "fasta_reason": "-", "fasta_inventory": "legacy: " + ",".join(sorted(leg_core)),
                 "dosage_diagnostic_fasta": "-", "dangling_references": len(dangling(leg)), "main_fasta": "-",
                 "dosage_mode": "audit", "dosage_mode_source": "-", "dosage_audit_status": "-",
                 "locus_dosage_audit_status": "-", "legacy_call_class": "-",
                 "audit_note": "declared differences only: " + ",".join(sorted(suppressed)) +
                 "; locus main_fasta, hap ids, dosage_representative_fasta",
                 "extraction": "-", "result": "PASS" if not diffs and not dangling(leg) else "FAIL: " + " | ".join(diffs)})
    print("\t".join(str(v) for v in rows[-1].values()), flush=True)
    # cold vs loaded persisted index (the prepared-locus command route: default audit, the same FASTA policy)
    idx = os.path.join(d, "F11u.panel.pvpl")
    if os.path.exists(idx):
        os.remove(idx)
    ld = os.path.join(d, "v2.loaded")
    os.makedirs(ld, exist_ok=True)
    argfile = os.path.join(d, "v2.loaded.args.txt")
    open(argfile, "w").write("".join(a + "\n" for a in args + ["--output-schema", "sequence_first_v2", "--out", f"{ld}/call"]))
    w = s1d.run([driver, "write-panel", str(4 << 30), "0", os.path.join(d, "drv_write"), idx, argfile],
                os.path.join(d, "index_write.log"))
    lc = s1d.run([driver, "loaded", str(4 << 30), "0", os.path.join(d, "drv_load"), idx, argfile],
                 os.path.join(d, "index_load.log"))
    files = ("call.sequence_locus.tsv", "call.sequence_haplotypes.tsv", "call.sequence_blocks.tsv",
             "call.representation_intervals.tsv", "call.sequence_alternatives.tsv", "call.locus_genotype.tsv",
             "call.dosage_blocks.tsv")
    same = [f for f in files if open(f"{v2}/{f}").read().replace(v2, "<D>") != open(f"{ld}/{f}").read().replace(ld, "<D>")]
    ds_a = {k: v for k, v in s1d.kv(f"{v2}/call.dosage_status.tsv").items() if not any(x in k for x in VOLATILE)}
    ds_b = {k: v for k, v in s1d.kv(f"{ld}/call.dosage_status.tsv").items() if not any(x in k for x in VOLATILE)}
    if ds_a != ds_b:
        same.append("call.dosage_status.tsv")
    inv_same = fasta_inventory(v2) == fasta_inventory(ld)
    loc, fields, problems = v2_checks(ld, "PANEL_REPRESENTED_PAIR")
    ok = w == 0 and lc == 0 and not same and inv_same and not problems
    rows.append({"section": "legacy", "case": "F11u:loaded_index", "exit": f"write {w} load {lc}", "prediction": loc.get("prediction"),
                 "hap_origins": f"{loc.get('hap1_origin')}/{loc.get('hap2_origin')}", "fasta_reason": loc.get("fasta_reason"),
                 **fields, "legacy_call_class": loc.get("legacy_call_class"), "audit_note": "loaded == cold",
                 "extraction": "-", "result": "PASS" if ok else f"FAIL: differ {same} inventory_same={inv_same} {problems}"})
    print("\t".join(str(v) for v in rows[-1].values()), flush=True)


def section_lifecycle(out, panvar, rows):
    nodes, frows, donor, _ = s1d.f11()
    d = os.path.join(out, "lifecycle")
    os.makedirs(d, exist_ok=True)
    gfa = os.path.join(d, "graph.gfa")
    s1d.write_gfa(gfa, nodes, frows)
    reads = os.path.join(d, "reads")
    s1a_fixtures.tiled_reads(nodes, donor, None, reads)
    prefix = os.path.join(d, "call")
    for f in os.listdir(d):
        if f.startswith("call."):
            os.remove(os.path.join(d, f))
    # stale run-owned FASTAs a legacy repair run could have left (non-empty): the first v2 run must remove them
    for f in ("call.mosaic_haplotypes.fa", "call.original.mosaic_haplotypes.fa", "call.dosage_repaired.fa",
              "call.dosage_original.fa", "call.dosage_repaired_block_alleles.fa", "call.sequence_haplotypes.fa"):
        open(os.path.join(d, f), "w").write(">stale\nACGT\n")
    curation = ["--panel-curation", os.path.join(out, "curation_F11.txt")]
    steps = [("off_panel", curation, "k0", [], 0, {"call.sequence_haplotypes.fa"}),
             ("panel_only", [], "k0", [], 0, set()),
             ("off_panel_again", curation, "k0", [], 0, {"call.sequence_haplotypes.fa"}),
             ("watchdog", [], "watchdog", [], 124, set()),
             ("off_panel_after_failure", curation, "k0", [], 0, {"call.sequence_haplotypes.fa"}),
             ("refused_repair", curation, "k0", ["--dosage-mode", "repair"], 2, None)]
    results, ok = [], True
    for label, extra, profile, flags, want_exit, want_inventory in steps:
        before = {f: open(os.path.join(d, f), "rb").read() for f in os.listdir(d) if f.startswith("call.")}
        args = ["--preset", "short-read", "--profile", os.path.join(out, f"profile_{profile}.tsv"), "--threads", "2",
                "--graph", gfa, "--reference-path", REF, *extra, "--reads1", os.path.join(reads, "reads_1.fq"),
                "--reads2", os.path.join(reads, "reads_2.fq"), "--library", os.path.join(reads, "library.txt"), *flags,
                "--output-schema", "sequence_first_v2", "--out", prefix, "--quiet"]
        code = s1d.run([panvar, "genotype", *args], os.path.join(d, f"{label}.log"))
        core, diagnostic = fasta_inventory(d)
        refs = dangling(d)
        if want_inventory is None:   # refused: byte-identical prefix
            after = {f: open(os.path.join(d, f), "rb").read() for f in os.listdir(d) if f.startswith("call.")}
            good = code == want_exit and after == before
            results.append(f"{label}: exit {code}, prefix {'unchanged' if after == before else 'CHANGED'}")
        else:
            seq = sorted(f for f in os.listdir(d) if f.startswith("call.sequence") or f.startswith("call.representation"))
            good = code == want_exit and core == want_inventory and not diagnostic and not refs
            if label == "watchdog":
                good = good and not seq and s1d.kv(prefix + ".mosaic_status.tsv").get("status") == "TIME_LIMIT"
            results.append(f"{label}: exit {code}, FASTA {','.join(sorted(core)) or 'none'}, dangling {len(refs)}" +
                           (f", sequence files {len(seq)}" if label == "watchdog" else ""))
        ok = ok and good
        if not good:
            results[-1] += " <- FAIL"
    rows.append({"section": "lifecycle", "case": "F11 one prefix (stale seeded)", "exit": "-", "prediction": "-",
                 "hap_origins": "-", "fasta_reason": "-", "fasta_inventory": "-", "dosage_diagnostic_fasta": "-",
                 "dangling_references": "-", "main_fasta": "-", "dosage_mode": "-", "dosage_mode_source": "-",
                 "dosage_audit_status": "-", "locus_dosage_audit_status": "-", "legacy_call_class": "-",
                 "audit_note": " | ".join(results), "extraction": "-", "result": "PASS" if ok else "FAIL"})
    print("\t".join(str(v) for v in rows[-1].values()), flush=True)


def section_crosslib(out, panvar, rows):
    """Earlier-library runs rerun with this library: legacy and sequence_first_v1 byte-identical (non-regression)."""
    s1c = os.path.join(RESULTS, "s1c_reporting_integrity_20261018/fixtures")
    s1dd = os.path.join(RESULTS, "s1d_panel_references_20261018/fixtures_a2")
    cases = [("S1-C F1 base legacy", os.path.join(s1c, "F1/base.legacy.log"), os.path.join(s1c, "F1/base.legacy")),
             ("S1-C F1 base v1", os.path.join(s1c, "F1/base.v1.log"), os.path.join(s1c, "F1/base.v1")),
             ("S1-C F7 base legacy", os.path.join(s1c, "F7/base.legacy.log"), os.path.join(s1c, "F7/base.legacy")),
             ("S1-C F7 base v1", os.path.join(s1c, "F7/base.v1.log"), os.path.join(s1c, "F7/base.v1")),
             ("S1-D F11u legacy", os.path.join(s1dd, "F11u/legacy.log"), os.path.join(s1dd, "F11u/legacy"))]
    d = os.path.join(out, "crosslib")
    os.makedirs(d, exist_ok=True)
    for label, log, old in cases:
        cmd = shlex.split(open(log).readline()[2:].strip())
        new = os.path.join(d, label.replace(" ", "_"))
        os.makedirs(new, exist_ok=True)
        for f in os.listdir(new):
            os.remove(os.path.join(new, f))
        cmd[0] = panvar
        i = cmd.index("--out")
        cmd[i + 1] = os.path.join(new, "call")
        code = s1d.run(cmd, new + ".log")
        diffs = compare_dirs(old, new)
        core, _ = fasta_inventory(new)
        refs = dangling(new)
        rows.append({"section": "crosslib", "case": label, "exit": code, "prediction": "-", "hap_origins": "-",
                     "fasta_reason": "-", "fasta_inventory": ",".join(sorted(core)) or "none", "dosage_diagnostic_fasta": "-",
                     "dangling_references": len(refs), "main_fasta": "-",
                     "dosage_mode": cmd[cmd.index("--dosage-mode") + 1] if "--dosage-mode" in cmd else "default",
                     "dosage_mode_source": "-", "dosage_audit_status": "-", "locus_dosage_audit_status": "-",
                     "legacy_call_class": "-", "audit_note": "byte comparison vs " + os.path.relpath(old, RESULTS),
                     "extraction": "-", "result": "PASS" if not diffs else "FAIL: " + " | ".join(diffs[:4])})
        print("\t".join(str(v) for v in rows[-1].values()), flush=True)


def main(out, panvar, driver, only=None):
    os.makedirs(out, exist_ok=True)
    base = "profile\tfixture-bounds\ntime_limit_seconds\t300\nfootprint_limit_bytes\t4294967296\n"
    profiles = {"default": base, "round1": base + "repair_rounds\t1\nmax_mosaic_candidates\t0\n",
                "k0": base + "max_mosaic_candidates\t0\n",
                "watchdog": "profile\ts1d-watchdog\ntime_limit_seconds\t0.000001\nfootprint_limit_bytes\t4294967296\n"}
    for name, text in profiles.items():
        open(os.path.join(out, f"profile_{name}.tsv"), "w").write(text)
    open(os.path.join(out, "curation_F11.txt"), "w").write("# S1-D: the only (2,1,2) template\nS2#1#c\n")
    sections = only or ["cases", "arms", "legacy", "lifecycle", "crosslib"]
    rows = []
    if "cases" in sections:
        section_cases(out, panvar, driver, rows)
    if "arms" in sections:
        section_arms(out, panvar, rows)
    if "legacy" in sections:
        section_legacy(out, panvar, driver, rows)
    if "lifecycle" in sections:
        section_lifecycle(out, panvar, rows)
    if "crosslib" in sections:
        section_crosslib(out, panvar, rows)
    path = os.path.join(out, "fixtures.tsv")
    kept = [r for r in csv.DictReader(open(path), delimiter="\t")] if os.path.exists(path) else []
    kept = [r for r in kept if r["section"] not in sections]
    with open(path, "w") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(kept + rows)
    return 0 if all(r["result"] == "PASS" for r in rows) else 1


def check_prefix(d, expect_reason, gfa):
    """The real-data check: whole-prefix policy, dosage default and extraction of every panel haplotype."""
    loc, fields, problems = v2_checks(d, expect_reason)
    check = subprocess.run([sys.executable, os.path.join(HERE, "../s1d/extract_panel_haplotype.py"), "--check", gfa,
                            f"{d}/call.sequence_haplotypes.tsv"], capture_output=True, text=True)
    if check.returncode != 0:
        problems.append("extraction FAILED")
    if fields["dosage_mode"] != "audit" or fields["dosage_mode_source"] != "contract_default_sequence_first_v2":
        problems.append("dosage default not audit")
    if fields["dosage_audit_status"] != "COMPLETE" or not os.path.exists(f"{d}/call.dosage_blocks.tsv"):
        problems.append(f"audit {fields['dosage_audit_status']}")
    print(check.stdout.strip())
    for k in ("prediction", "sequence_status", "hap1_origin", "hap1_panel_paths", "hap2_origin", "hap2_panel_paths",
              "fasta_emitted", "fasta_reason", "legacy_fasta", "graph_sha256", "model_fingerprint"):
        print(f"{k}\t{loc.get(k)}")
    for k, v in fields.items():
        print(f"{k}\t{v}")
    print("files\t" + ",".join(sorted(f for f in os.listdir(d) if f.startswith("call."))))
    print("result\t" + ("PASS" if not problems else "FAIL: " + " | ".join(problems)))
    return 0 if not problems else 1


if __name__ == "__main__":
    if sys.argv[1] == "--check-prefix":
        sys.exit(check_prefix(os.path.abspath(sys.argv[2]), sys.argv[3], os.path.abspath(sys.argv[4])))
    sys.exit(main(os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2]), os.path.abspath(sys.argv[3]),
                  sys.argv[4].split(",") if len(sys.argv) > 4 else None))
