#!/usr/bin/env python3
"""M15-D: the basic short-read command, run and checked on development controls only (no reserved family). Each example
runs its command(s), asserts the frozen expectations (DESIGN.md section 3) and exits 0 only if every expectation holds.
Observed values (estimated library, outcome, runtime) are printed for the record and never checked against truth.
Usage: examples.py <example> <unit dir>"""
import csv
import hashlib
import os
import shutil
import subprocess
import sys
import time

csv.field_size_limit(sys.maxsize)
WT = "/Users/davide.bolognini/panvar-worktrees/s1-d-panel-references"
R = "/Users/davide.bolognini/panvar/tests/results"
PANVAR = f"{WT}/build/panvar"
M11D = f"{R}/m11d_dev_pilot_20261020/runs"
CELLS = {
    "c4": dict(graph=f"{R}/genotype_c4_production_11fe681/graphs/NA19705.loo.gfa", ref="HG00096#1#haplotype1-0000024:31848049-32074186",
               cur=None, reads=f"{R}/genotype_c4_production_11fe681/reads/NA19705/cov5", archived=f"{M11D}/call.c4.NA19705/call"),
    "cyp2d6": dict(graph=f"{R}/cyp2d6_shortread_pilot_v1/data/graphs/loo.HG02769.gfa", ref="grch38#1#chr22:42031864-42245566",
                   cur=f"{WT}/tests/genotype_cyp2d6_curation.tsv", reads=f"{R}/cyp2d6_shortread_pilot_v1/data/reads/HG02769/cov5",
                   archived=f"{M11D}/call.cyp.HG02769/call"),
    "acot": dict(graph=f"{R}/genotype_perf_bounded_20261015/inputs/acot/acot.gfa", ref="GRCh38#0#chr14:73416151-73668669",
                 cur=f"{R}/m9h_jp11_validation_20261020/runs/prepare/curation.LOO.txt",
                 reads=f"{R}/m8c_copying_challenge_20261020/runs/prepare.a1/reads/HG03209/cov5", archived=f"{M11D}/call.acot.HG03209/call"),
}
# The accelerated configuration the preset must carry, each recorded with the preset as its source (DESIGN.md section 2).
ACCELERATED = {"union_scoring_backend": "indexed_fine_experimental", "adaptive_window_discovery": "fine_experimental",
               "indexed_union_store_layout": "exact", "indexed_union_evidence_layout": "fine_views_experimental",
               "adaptive_fused_seed": "1", "combination_prepared_library_terms": "1",
               "gap_recursion_backend": "fixed_geometry_experimental", "dosage_census_backend": "interned_experimental",
               "union_fragment_blocks": "fine_core_4096_experimental", "repair_round1_origin_admission": "reuse_only_experimental",
               "repair_junction_proposals": "experimental_v1", "early_incumbent": "experimental_v1", "output_schema": "sequence_first_v2",
               "compound_annotation": "contract_v2", "footprint_limit_bytes": "6442450944", "max_state_classes": "600000000",
               "repair_max_neighbours": "model", "emission_cache_bytes": "536870912", "indexed_union_store_bytes": "2147483648"}
# The archived M11-D runs (the previous accelerated candidate, by profile) differ from a no-profile declared-library run in
# exactly these settings: their declared watchdog guard, diagnostics, locus-specific neighbourhood cap and an explicit
# repair_joint_search off (the built-in value, absent from a default table).
DECLARED_DIFF = {"time_limit_seconds", "block_diagnostics", "repair_max_neighbours", "repair_joint_search"}


def T(p):
    return list(csv.DictReader((l for l in open(p) if not l.startswith("#")), delimiter="\t")) if os.path.exists(p) else []


def KV(p):
    return dict(l.rstrip("\n").split("\t", 1) for l in open(p) if "\t" in l and not l.startswith("#")) if os.path.exists(p) else {}


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest() if os.path.exists(p) else None


def fasta_key(p):
    """A2: a FASTA's content with the record names' prefix-derived part removed (records are named <prefix>_seq_hapN)."""
    if not os.path.exists(p):
        return None
    lines = []
    for line in open(p):
        if line.startswith(">"):
            name, _, rest = line[1:].partition(" ")
            line = ">" + name[name.find("_seq_"):] + " " + rest
        lines.append(line)
    return hashlib.sha256("".join(lines).encode()).hexdigest()


def untimed(p):
    """A2: a key/value or table file without its timing rows (keys starting seconds_ or containing _seconds)."""
    if not os.path.exists(p):
        return None
    keep = [l for l in open(p) if not (l.split("\t", 1)[0].startswith("seconds") or "_seconds" in l.split("\t", 1)[0])]
    return hashlib.sha256("".join(keep).encode()).hexdigest()


def cmd(cell, out, library=None, extra=(), reference=True, reads=None):
    c = CELLS[cell]
    rd = reads or c["reads"]
    a = [PANVAR, "genotype", "--graph", c["graph"]]
    if reference:
        a += ["--reference-path", c["ref"]]
    if c["cur"]:
        a += ["--panel-curation", c["cur"]]
    a += ["--reads1", f"{rd}/reads_1.fq", "--reads2", f"{rd}/reads_2.fq", "--threads", "4", "--out", out, "--quiet"]
    if library == "declared":
        a += ["--library", f"{c['reads']}/library.txt"]
    elif library:
        a += ["--library", library]
    return a + list(extra)


def run(a):
    print("$", " ".join(a).replace(R, "<results>").replace(WT, "<worktree>"))
    t = time.time()
    p = subprocess.run(a, capture_output=True, text=True)
    p.wall = time.time() - t
    print(f"exit {p.returncode}  wall {p.wall:.1f} s")
    if p.stderr.strip():
        print("stderr:", p.stderr.strip()[-600:])
    return p


def prediction_key(prefix):
    win = next((r for r in T(prefix + ".equivalence_set.tsv") if r["rank"] == "1" and r["representative"] == "true"), {})
    hs = {r["haplotype"]: (r["sha256"], r["origin"]) for r in T(prefix + ".sequence_haplotypes.tsv")}
    st = KV(prefix + ".mosaic_status.tsv")
    return (win.get("score"), win.get("hap1_alleles"), win.get("hap2_alleles"), tuple(sorted(hs.items())), st.get("status"),
            st.get("model_fingerprint"))


def settings(prefix):
    return {r["key"]: (r["value"], r["source"]) for r in T(prefix + ".settings.tsv")}


CHECKS = []


def check(ok, what):
    CHECKS.append((ok, what))
    print(("PASS " if ok else "FAIL ") + what)


def fasta_rule(cell, out):
    """The sequence-first v2 FASTA rule: both haplotypes on the panel -> no FASTA and extraction references; either
    off-panel -> one FASTA with both records, matching the table's SHA-256; never the legacy FASTA."""
    hs = T(out + ".sequence_haplotypes.tsv")
    fa = out + ".sequence_haplotypes.fa"
    off = any(r["origin"] != "PANEL" for r in hs)
    check(len(hs) == 2, f"{cell}: sequence_haplotypes.tsv lists two haplotypes")
    check(os.path.exists(fa) == off, f"{cell}: FASTA present iff a haplotype is off-panel (origins {[r['origin'] for r in hs]})")
    if off:
        recs, name = {}, None
        for line in open(fa):
            line = line.rstrip("\n")
            if line.startswith(">"):
                name = line[1:].split()[0]
                recs[name] = ""
            elif name:
                recs[name] += line
        check(len(recs) == 2 and {hashlib.sha256(v.encode()).hexdigest() for v in recs.values()} <= {r["sha256"] for r in hs},
              f"{cell}: both FASTA records match the table's SHA-256")
    else:
        p = run([sys.executable, f"{WT}/tests/analysis/s1d/extract_panel_haplotype.py", "--check", CELLS[cell]["graph"], out + ".sequence_haplotypes.tsv"])
        check(p.returncode == 0, f"{cell}: panel haplotypes recovered by extraction reference (hash and length checked)")
    check(not os.path.exists(out + ".mosaic_haplotypes.fa"), f"{cell}: no legacy FASTA under v2")


def outcome(prefix, p):
    sl, st = KV(prefix + ".sequence_locus.tsv"), KV(prefix + ".mosaic_status.tsv")
    pred = sl.get("prediction")
    exits = {"PREDICTION": 0, "PROVISIONAL": 3, "NO_PREDICTION": 1}
    check(pred in exits and p.returncode == exits[pred], f"exit {p.returncode} matches the reported outcome {pred} ({sl.get('completion')})")
    print(f"outcome: prediction {pred}, completion {sl.get('completion')}, sequence_status {sl.get('sequence_status')}, "
          f"status {st.get('status')}, no_call_reason {st.get('no_call_reason', 'NA')}")
    return pred


def report_library(prefix):
    for suffix in (".library_estimate.tsv", ".library_sensitivity.tsv"):
        f = prefix + suffix
        print(f"--- {suffix}:", "".join(open(f).readlines()[:40]) if os.path.exists(f) else "absent")
    print("library_provenance:", KV(prefix + ".mosaic_status.tsv").get("library_provenance"))
    return os.path.exists(prefix + ".library_estimate.tsv")


def resources(prefix):
    st = KV(prefix + ".mosaic_status.tsv")
    print(f"resources: seconds_total {st.get('seconds_total')}, phys_footprint_peak {st.get('resources_phys_footprint_peak_bytes')}, "
          f"repair_neighbour_bound {st.get('repair_neighbour_bound')}, repair_max_neighbours {st.get('repair_max_neighbours')}")


def noprofile(cell, d):
    """The ordinary command: no profile, no --library (auto), no dosage flag."""
    out = f"{d}/{cell}"
    p = run(cmd(cell, out))
    pred = outcome(out, p)
    s = settings(out)
    bad = {k: s.get(k) for k, v in ACCELERATED.items() if s.get(k) != (v, "preset short-read (default)")}
    check(not bad, f"{cell}: the accelerated configuration is selected from the preset ({bad or 'all keys'})")
    check("time_limit_seconds" not in s or s["time_limit_seconds"][0] in ("0", "0.0", "0.000000"), f"{cell}: no time limit")
    check(not any(v[1].startswith("profile") for v in s.values()), f"{cell}: no setting comes from a profile")
    st = KV(out + ".mosaic_status.tsv")
    check(st.get("repair_max_neighbours") == "model" and st.get("repair_neighbour_bound", "").isdigit(),
          f"{cell}: model neighbourhood bound recorded ({st.get('repair_neighbour_bound')})")
    est = report_library(out)
    check(est and st.get("library_provenance", "").startswith("auto: insert from"), f"{cell}: library estimated from the reads")
    loc = (T(out + ".locus_genotype.tsv") or [{}])[0]
    print(f"locus: genotype_filter {loc.get('genotype_filter')} ({loc.get('genotype_filter_reasons')}), dosage_filter "
          f"{loc.get('dosage_filter')} ({loc.get('dosage_filter_reasons')}), phase_filter {loc.get('phase_filter')}")
    check(loc.get("dosage_mode") == "off" and (loc.get("dosage_filter") == "NO_CALL" or
          (loc.get("dosage_filter") == "NOT_EVALUATED" and "dosage_mode_off" in loc.get("dosage_filter_reasons", ""))),
          f"{cell}: dosage not run, reported not evaluated")
    check(KV(out + ".sequence_locus.tsv").get("dosage_mode") in ("off", None), f"{cell}: the sequence report records dosage off")
    if pred in ("PREDICTION", "PROVISIONAL"):
        fasta_rule(cell, out)
    resources(out)
    return 0


def ex_noprofile_c4(d):
    return noprofile("c4", d)


def ex_noprofile_cyp2d6(d):
    return noprofile("cyp2d6", d)


def ex_noprofile_acot(d):
    return noprofile("acot", d)


def ex_auto_equivalence_c4(d):
    """--library auto equals an omitted --library (the no-profile C4 run)."""
    ref = os.path.join(os.path.dirname(d), "ex.noprofile_c4", "c4")
    out = f"{d}/c4.auto"
    p = run(cmd("c4", out, library="auto"))
    check(prediction_key(out) == prediction_key(ref), "omitted vs --library auto: identical prediction")
    for suffix in (".library_estimate.tsv", ".library_sensitivity.tsv", ".library_sensitivity_blocks.tsv", ".sequence_haplotypes.tsv",
                   ".sequence_blocks.tsv", ".block_genotypes.tsv", ".genotypes.tsv", ".equivalence_set.tsv"):
        check(untimed(out + suffix) == untimed(ref + suffix), f"omitted vs --library auto: {suffix} identical apart from timing rows "
              f"({'present' if sha(ref + suffix) else 'absent in both' if not sha(out + suffix) else 'absent in one'})")
    check(fasta_key(out + ".sequence_haplotypes.fa") == fasta_key(ref + ".sequence_haplotypes.fa"),
          "omitted vs --library auto: .sequence_haplotypes.fa identical apart from the prefix-derived record names")
    a, b = settings(out), settings(ref)
    check({k: v[0] for k, v in a.items()} == {k: v[0] for k, v in b.items()}, "omitted vs --library auto: every effective setting equal")
    x, y = KV(out + ".mosaic_status.tsv"), KV(ref + ".mosaic_status.tsv")
    lib = {k for k in set(x) | set(y) if "library" in k and x.get(k) != y.get(k)}
    print("library status rows that differ:", {k: (x.get(k), y.get(k)) for k in sorted(lib)})
    check(all("seconds" in k for k in lib - {"library_invocation"}) and x.get("library_invocation") == "explicit_auto"
          and y.get("library_invocation") == "omitted_means_auto",
          "omitted vs --library auto: library status rows equal apart from timing and the invocation provenance row")
    print(f"exit {p.returncode}")
    return 0


def declared(cell, d):
    """No profile, the declared library and --dosage-check: the previous accelerated candidate's core results (M11-D)."""
    out = f"{d}/{cell}"
    p = run(cmd(cell, out, library="declared", extra=["--dosage-check"]))
    check(p.returncode == 0, f"{cell}: exit 0")
    a = CELLS[cell]["archived"]
    check(prediction_key(out) == prediction_key(a), f"{cell}: prediction identical to the archived M11-D accelerated run")
    for suffix in (".sequence_haplotypes.tsv", ".equivalence_set.tsv"):
        check(sha(out + suffix) == sha(a + suffix), f"{cell}: {suffix} identical to the archived run")
    check(fasta_key(out + ".sequence_haplotypes.fa") == fasta_key(a + ".sequence_haplotypes.fa"),
          f"{cell}: .sequence_haplotypes.fa identical to the archived run apart from the prefix-derived record names")
    new, old = settings(out), settings(a)
    diff = sorted(k for k in set(new) | set(old) if new.get(k, (None,))[0] != old.get(k, (None,))[0])
    check(set(diff) <= DECLARED_DIFF, f"{cell}: settings differ from the archived run only in {sorted(DECLARED_DIFF)} (got {diff})")
    check(new.get("repair_joint_search") is None and old.get("repair_joint_search", ("off",))[0] == "off",
          f"{cell}: repair_joint_search is the built-in off in both")
    print("differing settings:", {k: (new.get(k, ('absent',))[0], old.get(k, ('absent',))[0]) for k in diff})
    ds, dsa = KV(out + ".dosage_status.tsv"), KV(a + ".dosage_status.tsv")
    src = KV(out + ".sequence_locus.tsv").get("dosage_mode_source")
    print(f"dosage: new {ds.get('dosage_audit_status')} (source {src}); archived {dsa.get('dosage_audit_status')} "
          f"(source {KV(a + '.sequence_locus.tsv').get('dosage_mode_source')})")
    strip = lambda m: {k: v for k, v in m.items() if not k.startswith("seconds_")}
    check(src == "flag --dosage-check" and strip(ds) == strip(dsa),
          f"{cell}: --dosage-check runs the archived audit (dosage status equal apart from timing)")
    st, sa = KV(out + ".mosaic_status.tsv"), KV(a + ".mosaic_status.tsv")
    print(f"runtime: new {st.get('seconds_total')} s vs archived {sa.get('seconds_total')} s; peak footprint new "
          f"{st.get('resources_phys_footprint_peak_bytes')} vs archived {sa.get('resources_phys_footprint_peak_bytes')}")
    return out


def ex_declared_c4(d):
    declared("c4", d)
    return 0


def ex_declared_cyp2d6(d):
    declared("cyp2d6", d)
    return 0


def ex_declared_acot(d):
    declared("acot", d)
    return 0


def dosage_absent(cell, d):
    """The same declared-library run without --dosage-check: the prediction is unchanged; dosage diagnostics and status differ."""
    pre = os.path.basename(d).split("ex.")[0]   # amendment units ("a2.ex.*") read their own declared run
    ref = os.path.join(os.path.dirname(d), f"{pre}ex.declared_{cell}", cell)
    out = f"{d}/{cell}"
    p = run(cmd(cell, out, library="declared"))
    check(p.returncode == 0, f"{cell}: exit 0")
    check(prediction_key(out) == prediction_key(ref), f"{cell}: absent vs --dosage-check: identical prediction")
    for suffix in (".sequence_haplotypes.tsv", ".sequence_haplotypes.fa", ".equivalence_set.tsv", ".mosaic_scores.tsv", ".mosaic_candidates.tsv"):
        check(sha(out + suffix) == sha(ref + suffix), f"{cell}: absent vs --dosage-check: {suffix} identical")
    la, lr = (T(out + ".locus_genotype.tsv") or [{}])[0], (T(ref + ".locus_genotype.tsv") or [{}])[0]
    for name, l in (("absent", la), ("--dosage-check", lr)):
        print(f"{name}: dosage_mode {l.get('dosage_mode')}, audit {l.get('dosage_audit_status')}, GENOTYPE {l.get('genotype_filter')}, "
              f"DOSAGE {l.get('dosage_filter')} ({l.get('dosage_filter_reasons')}), summary {l.get('call_summary')}")
    check(la.get("dosage_mode") == "off" and la.get("dosage_filter") == "NOT_EVALUATED" and "dosage_mode_off" in la.get("dosage_filter_reasons", ""),
          f"{cell}: without the flag dosage is off and DOSAGE_FILTER is NOT_EVALUATED (dosage_mode_off)")
    check(lr.get("dosage_mode") == "audit" and lr.get("dosage_audit_status") == "COMPLETE" and "dosage_mode_off" not in lr.get("dosage_filter_reasons", ""),
          f"{cell}: with --dosage-check the audit completes and the dosage filter is evaluated from it")
    check(la.get("genotype_filter") == lr.get("genotype_filter") and la.get("genotype_filter_reasons") == lr.get("genotype_filter_reasons"),
          f"{cell}: GENOTYPE_FILTER identical with and without the audit")
    check(KV(out + ".sequence_locus.tsv").get("dosage_mode") == "off" and KV(ref + ".sequence_locus.tsv").get("dosage_mode_source") == "flag --dosage-check",
          f"{cell}: the sequence report records off vs flag --dosage-check")
    only_ref = sorted(f[len(cell):] for f in os.listdir(os.path.dirname(ref)) if f.startswith(cell + ".")
                      and not os.path.exists(f"{d}/{f}"))
    only_out = sorted(f[len(cell):] for f in os.listdir(d) if f.startswith(cell + ".") and not os.path.exists(f"{os.path.dirname(ref)}/{f}"))
    print("files only with --dosage-check:", only_ref, "| only without:", only_out)
    check(all(f.startswith(".dosage") for f in only_ref) and not only_out, f"{cell}: the audit adds only dosage diagnostics")
    check(not any(f.endswith(".fa") and "dosage" not in f for f in only_ref), f"{cell}: the audit adds no prediction FASTA")
    st, sr = KV(out + ".mosaic_status.tsv"), KV(ref + ".mosaic_status.tsv")
    print(f"runtime: without audit {p.wall:.1f} s (seconds_total {st.get('seconds_total')}); with --dosage-check seconds_total "
          f"{sr.get('seconds_total')}")
    return 0


def ex_dosage_absent_c4(d):
    return dosage_absent("c4", d)


def ex_dosage_absent_cyp2d6(d):
    return dosage_absent("cyp2d6", d)


def ex_dosage_absent_acot(d):
    return dosage_absent("acot", d)


def ex_invalid(d):
    """Contradictory or incomplete commands fail with exit 2 and a message naming the problem, before touching outputs: a
    sentinel file already at the prefix is left byte-identical and nothing else is written."""
    def refused(label, a, needle, prefix, untouched=True):
        sentinel = prefix + ".settings.tsv"
        open(sentinel, "w").write("sentinel\n")
        p = run(a)
        msg = p.stderr + p.stdout
        left = sorted(f for f in os.listdir(d) if f.startswith(os.path.basename(prefix) + "."))
        kept = left == [os.path.basename(sentinel)] and open(sentinel).read() == "sentinel\n"
        if untouched:
            check(p.returncode == 2 and needle in msg and kept, f"refused ({label}): exit 2, message names '{needle}', the prefix untouched ({left})")
        else:   # an input error found when the run reads it: exit 2 and the message checked; the prefix state recorded only
            check(p.returncode == 2 and needle in msg, f"refused ({label}): exit 2, message names '{needle}'")
            print(f"recorded: prefix after the refusal {left}; sentinel kept {kept}")

    refused("--dosage-check with --dosage-mode off", cmd("c4", f"{d}/r1", "declared", ["--dosage-check", "--dosage-mode", "off"]),
            "--dosage-check contradicts", f"{d}/r1")
    refused("--dosage-mode off then --dosage-check", cmd("c4", f"{d}/r2", "declared", ["--dosage-mode", "off", "--dosage-check"]),
            "--dosage-check contradicts", f"{d}/r2")
    refused("--dosage-check with --dosage-audit off", cmd("c4", f"{d}/r3", "declared", ["--dosage-check", "--dosage-audit", "off"]),
            "--dosage-check contradicts", f"{d}/r3")
    refused("--dosage-mode repair under the default v2 report", cmd("c4", f"{d}/r4", "declared", ["--dosage-mode", "repair"]),
            "--output-schema legacy --dosage-mode repair", f"{d}/r4")
    refused("--dosage-check with the legacy repair", cmd("c4", f"{d}/r5", "declared", ["--dosage-check", "--output-schema", "legacy", "--dosage-mode", "repair"]),
            "--dosage-check contradicts", f"{d}/r5")
    refused("explicit sequence_first_v2 with the unrepaired view", cmd("c4", f"{d}/r6", "declared", ["--output-schema", "sequence_first_v2", "--adaptive-repair", "off"]),
            "needs the adaptive repair", f"{d}/r6")
    refused("no --reference-path on a graph without RS:Z", cmd("c4", f"{d}/r7", "declared", reference=False), "--reference-path", f"{d}/r7")
    refused("a missing library file", cmd("c4", f"{d}/r8", f"{d}/no_such_library.txt"), "no_such_library.txt", f"{d}/r8", untouched=False)
    refused("--preset long-read", cmd("c4", f"{d}/r9", "declared", ["--preset", "long-read"]), "not implemented", f"{d}/r9")
    return 0


def ex_unrepaired_view(d):
    """--adaptive-repair off without a profile: the preset's repair-only features switch off and are recorded inactive."""
    out = f"{d}/c4.off"
    p = run(cmd("c4", out, "declared", ["--adaptive-repair", "off"]))
    check(p.returncode in (0, 1), f"unrepaired view: runs under the legacy schema (exit {p.returncode})")
    s = settings(out)
    for k, v in (("output_schema", "legacy"), ("early_incumbent", "off"), ("compound_annotation", "off"), ("repair_junction_proposals", "off")):
        check(s.get(k, ("", ""))[0] == v and "(inactive: " in s.get(k, ("", ""))[1], f"unrepaired view: {k} {v}, recorded inactive ({s.get(k)})")
    check(s.get("union_scoring_backend", ("", ""))[0] == "indexed_fine_experimental", "unrepaired view: exact backends stay selected")
    check(not os.path.exists(out + ".sequence_locus.tsv"), "unrepaired view: no sequence-first report")
    return 0


def ex_auto_refused(d):
    """The automatic library refuses (LIBRARY_REFUSED) on too few pairs; no fallback, no prediction."""
    reads = f"{d}/reads"
    os.makedirs(reads, exist_ok=True)
    for m in ("1", "2"):
        with open(f"{CELLS['c4']['reads']}/reads_{m}.fq") as src, open(f"{reads}/reads_{m}.fq", "w") as dst:
            for i, line in enumerate(src):
                if i >= 300 * 4:
                    break
                dst.write(line)
    out = f"{d}/c4.few"
    p = run(cmd("c4", out, reads=reads))
    st = KV(out + ".mosaic_status.tsv")
    est = report_library(out)
    print(f"status {st.get('status')}; status rows naming the refusal:", {k: v for k, v in st.items() if "refus" in k or "library" in k})
    check(st.get("status") == "LIBRARY_REFUSED", "too few pairs: status LIBRARY_REFUSED")
    check(p.returncode == 1, f"too few pairs: exit 1 (NO_PREDICTION), got {p.returncode}")
    sl = KV(out + ".sequence_locus.tsv")
    check(sl.get("prediction", "NO_PREDICTION") == "NO_PREDICTION", f"too few pairs: no prediction ({sl.get('prediction', 'no sequence report')})")
    check(not os.path.exists(out + ".sequence_haplotypes.fa") and not os.path.exists(out + ".mosaic_haplotypes.fa"), "too few pairs: no FASTA")
    check(bool(est), "too few pairs: the estimate and its refusal are written")
    return 0


def ex_resource_stop(d):
    """An expert override (one key) binds a time limit on the ordinary ACOT command: never reported as a completed prediction."""
    prof = f"{d}/time8.tsv"
    open(prof, "w").write("profile\tm15d-time-limit-8s\ntime_limit_seconds\t8\n")
    out = f"{d}/acot.tight"
    p = run(cmd("acot", out, "declared", ["--profile", prof]))
    sl, st = KV(out + ".sequence_locus.tsv"), KV(out + ".mosaic_status.tsv")
    print(f"prediction {sl.get('prediction')} completion {sl.get('completion')} status {st.get('status')}")
    check(p.returncode in (3, 124), f"binding time limit: exit {p.returncode} (3 PROVISIONAL after a checkpoint, or 124 before one)")
    check((p.returncode == 3 and sl.get("prediction") == "PROVISIONAL") or (p.returncode == 124 and sl.get("prediction") != "PREDICTION"),
          "binding time limit: never reported as a completed prediction")
    s = settings(out)
    if os.path.exists(out + ".settings.tsv"):
        check(s.get("time_limit_seconds", ("", ""))[1].startswith("profile") and s.get("union_scoring_backend", ("", ""))[1].startswith("preset"),
              "binding time limit: the override is recorded from the profile, everything else from the preset")
    else:   # A2: a stop before any checkpoint writes the minimal status only (documented); it names the profile and the limit
        check(st.get("profile") == prof and st.get("time_limit_seconds") == "8" and st.get("status") == "TIME_LIMIT",
              "binding time limit: the minimal status names the override profile and its 8 s limit")
    print(f"recorded: sequence_locus exit_code {sl.get('exit_code')} vs process exit {p.returncode}")
    return 0


if __name__ == "__main__":
    name, d = sys.argv[1], sys.argv[2]
    os.makedirs(d, exist_ok=True)
    globals()["ex_" + name](d)
    bad = [w for ok, w in CHECKS if not ok]
    print(f"{len(CHECKS) - len(bad)} of {len(CHECKS)} checks passed" + (f"; FAILED: {bad}" if bad else ""))
    sys.exit(1 if bad else 0)
