#!/usr/bin/env python3
"""M15-G verification units (frozen before execution; development controls only, no reserved family).
Usage: examples.py <unit> <unit dir>"""
import hashlib
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "m15d"))
import assessment_spec as S  # noqa: E402
import examples as E  # noqa: E402  (the M15-D development cells and readers)

R = E.R
WT = E.WT
ASSESS = f"{WT}/build/genotype_assessment_test"
M15E = f"{R}/m15e_panel_index_20261021/runs"
CHECKS = []
READ_SUFFIXES = [".sequence_locus.tsv", ".sequence_blocks.tsv", ".mosaic_status.tsv", ".block_genotypes.tsv", ".dosage_status.tsv",
                 ".dosage_blocks.tsv"]


def check(ok, what):
    CHECKS.append((ok, what))
    print(("PASS " if ok else "FAIL ") + what)


def run(a):
    print("$", " ".join(a).replace(R, "<results>").replace(WT, "<worktree>"))
    p = subprocess.run(a, capture_output=True, text=True)
    print(f"exit {p.returncode}")
    if p.stderr.strip():
        print("stderr:", p.stderr.strip()[-500:])
    return p


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest() if os.path.exists(p) else None


def cpp_matches_spec(prefix, what):
    want = S.render(prefix)
    got = subprocess.run([ASSESS, "show", prefix], capture_output=True, text=True).stdout
    check(got == want, f"{what}: the C++ view equals the specification's")
    return got == want


def written_view(prefix):
    """The view as WRITTEN by the run (locus rows in .sequence_locus.tsv, block columns in .sequence_blocks.tsv)."""
    sl = E.KV(prefix + ".sequence_locus.tsv")
    blocks = E.T(prefix + ".sequence_blocks.tsv")
    return {k: sl.get(k) for k in S.LOCUS_KEYS}, [(b["block_id"], {c: b.get(c) for c in S.BLOCK_COLUMNS}) for b in blocks]


def written_equals_spec(prefix, what):
    """What the run wrote equals the specification evaluated on the run's own (final) files."""
    L, Bk = S.assess(S.read(prefix))
    wl, wb = written_view(prefix)
    check(wl == L, f"{what}: the written locus rows equal the specification ({wl.get('assessment_summary')})")
    check(wb == Bk, f"{what}: the written block columns equal the specification ({len(wb)} blocks)")
    print(f"{what}: summary {L['assessment_summary']}\n    findings={L['findings']}\n    not_assessed={L['not_assessed']}\n"
          f"    alternatives={L['dosage_alternatives']}")
    return L


def without(path, drop_cols):
    """A TSV's text with the named columns removed (to compare the pre-existing columns only)."""
    if not os.path.exists(path):
        return None
    lines = open(path).read().rstrip("\n").split("\n")
    head = lines[0].split("\t")
    keep = [i for i, c in enumerate(head) if c not in drop_cols]
    return "\n".join("\t".join(r.split("\t")[i] if i < len(r.split("\t")) else "" for i in keep) for r in lines)


def named(text, *names):
    for n in names:
        text = text.replace(n, "<RUN>")
    return text


# ---------------------------------------------------------------- units ------------------------------------------------
def ex_replay_saved(d):
    """Every saved development v2 prefix (reserved excluded): the C++ view equals the specification's, read-only."""
    prefixes = []
    for root, dirs, files in os.walk(R):
        if "m15b_reserved" in root or "/index" in root or root.startswith(f"{R}/m15g"):
            dirs[:] = []
            continue
        for f in files:
            if f.endswith(".sequence_locus.tsv"):
                p = os.path.join(root, f[: -len(".sequence_locus.tsv")])
                if E.KV(p + ".sequence_locus.tsv").get("schema_version") == "sequence_first_v2":
                    prefixes.append(p)
    prefixes.sort()
    before = {p + s: sha(p + s) for p in prefixes for s in READ_SUFFIXES}
    bad, kinds = [], {}
    for p in prefixes:
        want, got = S.render(p), subprocess.run([ASSESS, "show", p], capture_output=True, text=True).stdout
        if got != want:
            bad.append(p)
        L, _ = S.assess(S.read(p))
        kinds[L["assessment_summary"].split(";")[0]] = kinds.get(L["assessment_summary"].split(";")[0], 0) + 1
    after = {k: sha(k) for k in before}
    print(f"replayed {len(prefixes)} saved v2 prefixes; outcomes {kinds}; differing: {[b.replace(R + '/', '') for b in bad[:5]]}")
    check(not bad, f"replay: the C++ view equals the specification on {len(prefixes) - len(bad)} of {len(prefixes)} saved prefixes")
    check(before == after, f"replay: {len(before)} source files byte-identical before and after (read-only)")
    return 0


def complete(cell, d):
    """A completed no-profile run: the primary view written and equal to the specification; predictions, candidate
    scores, sequences, phase inference and every pre-existing column (old headline fields included) unchanged from M15-E's
    run of the same command."""
    out = f"{d}/{cell}"
    p = run(E.cmd(cell, out))
    old = f"{M15E}/ex.index_{cell}/{cell}.cold"
    check(p.returncode == 0, f"{cell}: exit 0")
    written_equals_spec(out, cell)
    cpp_matches_spec(out, cell)
    check(E.prediction_key(out) == E.prediction_key(old), f"{cell}: prediction identical to M15-E (score, paths, haplotypes, status, fingerprint)")
    for suffix in (".sequence_haplotypes.tsv", ".mosaic_scores.tsv", ".mosaic_candidates.tsv", ".equivalence_set.tsv",
                   ".representation_intervals.tsv", ".sequence_alternatives.tsv", ".sequence_compound.tsv"):
        a, b = open(out + suffix).read() if os.path.exists(out + suffix) else None, open(old + suffix).read() if os.path.exists(old + suffix) else None
        check(a is not None and named(a, out, os.path.basename(out)) == named(b, old, os.path.basename(old)), f"{cell}: {suffix} unchanged")
    check(E.fasta_key(out + ".sequence_haplotypes.fa") == E.fasta_key(old + ".sequence_haplotypes.fa"), f"{cell}: FASTA unchanged")
    for suffix in (".sequence_blocks.tsv", ".block_genotypes.tsv", ".genotypes.tsv"):
        a, b = without(out + suffix, S.BLOCK_COLUMNS), without(old + suffix, S.BLOCK_COLUMNS)
        check(a is not None and a == b, f"{cell}: {suffix} pre-existing columns unchanged (old headline columns included)")
    a = without(out + ".locus_genotype.tsv", S.LOCUS_KEYS + ["headline_filters_status"])
    b = without(old + ".locus_genotype.tsv", [])
    check(named(a, os.path.basename(out)) == named(b, os.path.basename(old)), f"{cell}: locus table pre-existing columns unchanged")
    check(E.T(out + ".locus_genotype.tsv")[0].get("headline_filters_status") == "superseded_by_assessment_v1",
          f"{cell}: the old headline fields are labelled superseded")
    so, sn = E.KV(old + ".mosaic_status.tsv"), E.KV(out + ".mosaic_status.tsv")
    for k in ("status", "call_class", "phase_support", "no_call_reason", "library_sensitivity_status", "library_sensitive_blocks"):
        check(sn.get(k) == so.get(k), f"{cell}: status {k} unchanged ({sn.get(k)})")
    lo = {k: v for k, v in E.KV(old + ".sequence_locus.tsv").items()}
    ln = {k: v for k, v in E.KV(out + ".sequence_locus.tsv").items() if k not in S.LOCUS_KEYS}
    diff = sorted(k for k in set(lo) | set(ln) if named(str(lo.get(k)), old, os.path.basename(old)) != named(str(ln.get(k)), out, os.path.basename(out))
                  and "seconds" not in k)
    check(not diff, f"{cell}: the sequence record's pre-existing rows unchanged ({diff or 'none'})")
    return 0


def ex_complete_c4(d):
    return complete("c4", d)


def ex_complete_cyp2d6(d):
    return complete("cyp2d6", d)


def ex_complete_acot(d):
    return complete("acot", d)


def ex_dosage_cyp2d6(d):
    """The known CYP2D6 audit failure: dosage UNAVAILABLE with its reason, never agreement; the prediction unchanged."""
    out = f"{d}/cyp2d6"
    p = run(E.cmd("cyp2d6", out, "declared", ["--dosage-check"]))
    check(p.returncode == 0, "cyp2d6 --dosage-check: exit 0")
    L = written_equals_spec(out, "cyp2d6 --dosage-check")
    cpp_matches_spec(out, "cyp2d6 --dosage-check")
    check(L["dosage_check"] == "UNAVAILABLE" and "DOSAGE:AUDIT_FAILED(unit library: alleles disagree on their flanks or anchor trims in sv:9)@locus" in L["not_assessed"],
          "cyp2d6: dosage UNAVAILABLE, the audit failure named once at the locus")
    check(E.prediction_key(out) == E.prediction_key(f"{R}/m15d1_correctness_20261021/runs/ex.cyp_audit/cyp2d6"),
          "cyp2d6: prediction identical to M15-D.1's run of the same command")
    return 0


def ex_legacy_tie(d):
    """A v2 PREDICTION beside a legacy distinct_tie (no legacy block table): primary rows exist, checks NOT_ASSESSED."""
    F = f"{R}/s1b_sequence_first_20261018/fixtures"
    out = f"{d}/call"
    p = run([E.PANVAR, "genotype", "--profile", f"{F}/profile_default.tsv", "--threads", "2", "--graph", f"{F}/F1/base.gfa",
             "--reference-path", "ref#0#chain", "--reads1", f"{R}/representation_fixtures_20261015/F1/reads/reads_1.fq",
             "--reads2", f"{R}/representation_fixtures_20261015/F1/reads/reads_2.fq", "--library",
             f"{R}/representation_fixtures_20261015/F1/reads/library.txt", "--quiet", "--out", out])
    sl = E.KV(out + ".sequence_locus.tsv")
    check(p.returncode == 0 and sl.get("prediction") == "PREDICTION" and sl.get("legacy_no_call_reason") == "distinct_tie",
          f"legacy tie: v2 PREDICTION beside legacy distinct_tie (exit {p.returncode})")
    check(not os.path.exists(out + ".block_genotypes.tsv"), "legacy tie: no legacy block table exists")
    L = written_equals_spec(out, "legacy tie")
    cpp_matches_spec(out, "legacy tie")
    check(L["sequence_resolution"] == "RESOLVED" and L["genotype_check"] == "NOT_ASSESSED" and L["read_support"] == "NOT_ASSESSED"
          and "LABELS_AMBIGUOUS" in L["findings"] and L["phase"] == "NOT_ASSESSED", "legacy tie: resolution kept, unverifiable checks NOT_ASSESSED")
    blocks = E.T(out + ".sequence_blocks.tsv")
    check(bool(blocks) and all(b.get("genotype_check") == "NOT_ASSESSED" for b in blocks), "legacy tie: primary block rows exist, every check NOT_ASSESSED")
    return 0


def stop(d, name, point, code, prediction):
    out = f"{d}/c4"
    prof = f"{d}/{name}.tsv"
    open(prof, "w").write(f"profile\tm15g-{name}\ntest_interrupt_point\t{point}\n")
    p = run(E.cmd("c4", out, "declared", ["--profile", prof]))
    check(p.returncode == code, f"{name}: exit {code} (got {p.returncode})")
    L = written_equals_spec(out, name)
    cpp_matches_spec(out, name)
    check(E.KV(out + ".sequence_locus.tsv").get("prediction") == prediction, f"{name}: {prediction}")
    return L


def ex_provisional(d):
    L = stop(d, "provisional", "after_checkpoint", 3, "PROVISIONAL")
    check(L["sequence_resolution"] == "NOT_DETERMINED" and L["findings"] == "none" and L["read_support"] == "NOT_ASSESSED"
          and L["phase"] == "NOT_ASSESSED", "provisional: no decision set and no completed support manufactured")
    blocks = E.T(f"{d}/c4.sequence_blocks.tsv")
    check(bool(blocks) and all(b.get("block_resolution") == "PROVISIONAL_LABEL" and b.get("read_support") == "NOT_ASSESSED" for b in blocks),
          "provisional: block rows labelled PROVISIONAL_LABEL, support NOT_ASSESSED")
    return 0


def ex_absent(d):
    L = stop(d, "absent", "before_first_evaluation", 124, "NO_PREDICTION")
    check(L["genotype_check"] == "NA" and L["phase"] == "NA" and L["assessment_summary"].startswith("NO_PREDICTION;"),
          "absent prediction: every assessment NA")
    return 0


def ex_reporting_failure(d):
    out = f"{d}/c4"
    os.makedirs(f"{out}.locus_genotype.tsv.writing/blocker", exist_ok=True)
    p = run(E.cmd("c4", out, "declared"))
    sl = E.KV(out + ".sequence_locus.tsv")
    check(p.returncode == 2 and sl.get("exit_code") == "2" and sl.get("locus_output", "").startswith("FAILED: "),
          "reporting failure: exit 2, recorded in the sequence record")
    L = written_equals_spec(out, "reporting failure")
    cpp_matches_spec(out, "reporting failure")
    check(L["assessment_rule"] == "assessment_v1", "reporting failure: the primary view of the published prediction is still written")
    os.rmdir(f"{out}.locus_genotype.tsv.writing/blocker")
    os.rmdir(f"{out}.locus_genotype.tsv.writing")
    return 0


if __name__ == "__main__":
    name, d = sys.argv[1], sys.argv[2]
    os.makedirs(d, exist_ok=True)
    globals()["ex_" + name](d)
    bad = [w for ok, w in CHECKS if not ok]
    print(f"{len(CHECKS) - len(bad)} of {len(CHECKS)} checks passed" + (f"; FAILED: {bad}" if bad else ""))
    sys.exit(1 if bad else 0)
