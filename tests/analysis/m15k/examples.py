#!/usr/bin/env python3
"""M15-K units (DESIGN.md section 4): refined runs unchanged against M15-J; unrefined runs against the M15-J binary, every
difference classified as expected (DESIGN.md section 1) or unexpected; the index build/load with producer provenance; a rebuild
after a new commit without reconfiguring. Usage: examples.py <example> <unit dir>"""
import csv
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "m15d"))
import examples as E  # noqa: E402
import importlib.util  # noqa: E402


def _load(name, rel):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, "..", rel, "examples.py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


M15E_X = _load("m15e_examples", "m15e")
M15J_X = _load("m15j_examples", "m15j")

R, WT = E.R, E.WT
M15J = f"{R}/m15j_interface_20261022/runs"
OUT = f"{R}/m15k_unrefined_v2_20261022"
BASELINE = f"{OUT}/baseline/panvar.m15j"   # the M15-J binary (57d1ffb2...), copied before any M15-K change
CHECKS = []
VOLATILE = M15J_X.VOLATILE
ASSESS_BLOCK = ["block_resolution", "genotype_check", "genotype_check_reasons", "dosage_check", "dosage_check_reason",
                "dosage_alternative", "read_support", "read_support_reason"]
ASSESS_LOCUS = ["assessment_rule", "sequence_resolution", "genotype_check", "genotype_check_coverage", "dosage_check",
                "dosage_check_coverage", "read_support", "read_support_coverage", "phase", "phase_reason", "findings",
                "not_assessed", "dosage_alternatives", "assessment_summary", "headline_filters_status"]
V2_FILES = {".sequence_locus.tsv", ".sequence_blocks.tsv", ".sequence_haplotypes.tsv", ".sequence_haplotypes.fa",
            ".sequence_alternatives.tsv", ".representation_intervals.tsv", ".sequence_compound.tsv"}
LOGS = (".repair_memory.tsv", ".repair_rounds.tsv")


def check(ok, what):
    CHECKS.append((ok, what))
    print(("PASS " if ok else "FAIL ") + what)
    return ok


M15J_X.CHECKS = CHECKS
M15J_X.check = check


# ---- refined runs: unchanged against M15-J ------------------------------------------------------------------------------
def refined(cell, d):
    out = f"{d}/{cell}"
    p = E.run(E.cmd(cell, out))
    check(p.returncode == 0, f"{cell}: exit 0 (refined, automatic library)")
    M15J_X.same_outputs(out, f"{M15J}/ex.ordinary_{cell}/{cell}", f"{cell} refined vs M15-J")
    return 0


def ex_refined_c4(d):
    return refined("c4", d)


def ex_refined_cyp2d6(d):
    return refined("cyp2d6", d)


def ex_refined_acot(d):
    return refined("acot", d)


# ---- unrefined runs: M15-K against the M15-J binary, every difference classified -------------------------------------------
def rows(text):
    return list(csv.DictReader(text.splitlines(), delimiter="\t"))


def kv(text):
    return dict(l.split("\t", 1) for l in text.splitlines() if "\t" in l)


def unrefined(cell, d):
    new, old = f"{d}/new/{cell}", f"{d}/m15j/{cell}"
    os.makedirs(os.path.dirname(new), exist_ok=True)
    os.makedirs(os.path.dirname(old), exist_ok=True)
    a = E.cmd(cell, new, extra=["--adaptive-repair", "off"])
    b = E.cmd(cell, old, extra=["--adaptive-repair", "off"])
    b[0] = BASELINE
    p, q = E.run(a), E.run(b)
    expected, unexpected = [], []
    nb, ob = os.path.basename(new), os.path.basename(old)
    nf = {f[len(nb):] for f in os.listdir(os.path.dirname(new)) if f.startswith(nb + ".")}
    of = {f[len(ob):] for f in os.listdir(os.path.dirname(old)) if f.startswith(ob + ".")}
    for f in sorted(nf - of):
        (expected if f in V2_FILES else unexpected).append(f"{f}: new file")
    for f in sorted(of - nf):
        (expected if f == ".mosaic_haplotypes.fa" else unexpected).append(f"{f}: no longer written")
    sl = E.KV(new + ".sequence_locus.tsv")
    st_new, st_old = E.KV(new + ".mosaic_status.tsv"), E.KV(old + ".mosaic_status.tsv")
    phase_changed = st_new.get("phase_support") != st_old.get("phase_support")
    for f in sorted(nf & of):
        if f in LOGS:
            continue
        x = open(new + f, errors="replace").read().replace(new, "<RUN>").replace(nb + "_seq", "<RUN>_seq")
        y = open(old + f, errors="replace").read().replace(old, "<RUN>").replace(ob + "_seq", "<RUN>_seq")
        if f in (".mosaic_candidates.tsv", ".mosaic_scores.tsv"):
            check(x == y, f"{cell}: {f} identical (candidates, objective components, union ranking)")
            continue
        if x == y:
            continue
        if f in (".settings.tsv",):
            a_, b_ = {r["key"]: r for r in rows(x)}, {r["key"]: r for r in rows(y)}
            for k in sorted(set(a_) | set(b_)):
                if a_.get(k) == b_.get(k):
                    continue
                ok = (k == "output_schema" and a_[k]["value"] == "sequence_first_v2") or \
                     (k == "compound_annotation" and a_[k]["value"] == "contract_v2") or \
                     (k == "early_incumbent" and a_[k]["value"] == b_[k]["value"] == "off")
                (expected if ok else unexpected).append(f"{f}: {k} {b_.get(k, {}).get('value')} -> {a_.get(k, {}).get('value')}")
            continue
        if f in (".mosaic_status.tsv", ".dosage_status.tsv", ".library_estimate.tsv", ".library_sensitivity.tsv"):
            if f == ".library_estimate.tsv":   # a table: timing rows aside
                keep = lambda t: [l for l in t.splitlines() if "seconds" not in l.split("\t")[0]]
                (expected if keep(x) == keep(y) else unexpected).append(f"{f}: " + ("timing rows only" if keep(x) == keep(y) else "content"))
                continue
            a_, b_ = kv(x), kv(y)
            for k in sorted(set(a_) | set(b_)):
                if a_.get(k) == b_.get(k) or any(v in k for v in VOLATILE):
                    continue
                ok = k.startswith("sequence_") or k.startswith("compound_") or k == "output_schema" or \
                     (k == "phase_support" and a_.get(k) == sl.get("legacy_phase_support"))
                (expected if ok else unexpected).append(f"{f}: {k} {b_.get(k)} -> {a_.get(k)}")
            continue
        if f in (".block_genotypes.tsv", ".genotypes.tsv", ".locus_genotype.tsv"):
            ra, rb = rows(x), rows(y)
            added = [c for c in (ra[0].keys() if ra else []) if rb and c not in rb[0]]
            assess = ASSESS_LOCUS if f == ".locus_genotype.tsv" else ASSESS_BLOCK
            for c in added:
                (expected if c in assess else unexpected).append(f"{f}: column {c} added")
            if len(ra) != len(rb):
                unexpected.append(f"{f}: {len(rb)} -> {len(ra)} rows")
                continue
            changed = sorted({c for r1, r2 in zip(ra, rb) for c in r2 if c in r1 and r1[c] != r2[c]})
            for c in changed:
                ok = (f == ".locus_genotype.tsv" and c in ("main_fasta", "prefix")) or (phase_changed and "phase" in c)
                (expected if ok else unexpected).append(f"{f}: column {c} values changed")
            continue
        unexpected.append(f"{f}: content differs")
    print(f"{cell}: exits M15-J {q.returncode} -> M15-K {p.returncode}")
    for e in expected:
        print(f"  expected   {e}")
    for u in unexpected:
        print(f"  UNEXPECTED {u}")
    check(not unexpected, f"{cell}: every difference from M15-J is one DESIGN.md section 1 expects ({len(expected)} expected, "
                          f"{len(unexpected)} unexpected)")
    pred = sl.get("prediction")
    check(pred in ("PREDICTION", "PROVISIONAL", "NO_PREDICTION") and
          p.returncode == (0 if pred == "PREDICTION" else 3 if pred == "PROVISIONAL" else 1),
          f"{cell}: v2 {pred} {sl.get('sequence_status')}, exit {p.returncode} follows the v2 contract")
    if pred == "PREDICTION":
        scores = rows(open(new + ".mosaic_scores.tsv").read())
        best = max(float(r["score"]) for r in scores)
        rep = float(sl.get("representative_score", "nan"))
        check(abs(rep - best) <= 1e-6 * max(1.0, abs(best)), f"{cell}: the representative's score {rep} is the best union score {best}")
    check(st_new.get("phase_support") == sl.get("legacy_phase_support"),
          f"{cell}: the status phase is the pooled decision's ({st_new.get('phase_support')}; M15-J {st_old.get('phase_support')})")
    print(f"{cell}: assessment {sl.get('assessment_summary')}")
    return 0


def ex_unrefined_c4(d):
    return unrefined("c4", d)


def ex_unrefined_cyp2d6(d):
    return unrefined("cyp2d6", d)


def ex_unrefined_acot(d):
    return unrefined("acot", d)


# ---- the index: producer provenance apart from compatibility ----------------------------------------------------------------
def head_commit():
    return subprocess.run(["git", "-C", WT, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()


def build_record(text):
    return dict(l.split("\t", 1) for l in text.strip().splitlines()[1:] if "\t" in l)


def payload_equal_but_provenance(new, old, new_commit, new_tree, old_commit, old_tree):
    """The index equals `old` apart from the header digest [24, 56) and the two provenance strings (whose lengths may differ):
    every byte before the commit (outside the digest) and every byte after the tree string identical."""
    a, b = open(new, "rb").read(), open(old, "rb").read()
    ia, ib = a.find(new_commit.encode()), b.find(old_commit.encode())
    if ia < 0 or ib < 0 or ia != ib:
        return False, f"commit offsets {ia}, {ib}"
    before_ok = all(a[i] == b[i] for i in range(ia) if not 24 <= i < 56)
    ta = a.find(new_tree.encode(), ia + len(new_commit))
    tb = b.find(old_tree.encode(), ib + len(old_commit))
    if ta < 0 or tb < 0 or ta - ia > 64 or tb - ib > 64:
        return False, "tree strings not found next to the commit"
    tail_a, tail_b = a[ta + len(new_tree):], b[tb + len(old_tree):]
    return before_ok and tail_a == tail_b, f"{len(a)} vs {len(b)} bytes; length prefixes between the strings may differ"


def ex_index_c4(d):
    idx_dir = f"{d}/index"
    os.makedirs(idx_dir, exist_ok=True)
    idx = f"{idx_dir}/c4.pvpl"
    p = E.run(M15E_X.locus_cmd("c4", ["--build-panel-index", idx]))
    rec = build_record(p.stdout)
    open(f"{idx_dir}/c4.build.tsv", "w").write(p.stdout)
    check(p.returncode == 0 and os.path.exists(idx), "c4: index built")
    head = head_commit()
    check(rec.get("producer_source_commit") == head and rec.get("producer_source_tree") == "clean",
          f"c4: producer provenance is the source built: {rec.get('producer_source_commit')} {rec.get('producer_source_tree')} "
          f"(HEAD {head})")
    old = f"{M15J}/a1.ex.index_c4/index/c4.pvpl"
    ok, how = payload_equal_but_provenance(idx, old, rec.get("producer_source_commit", ""), rec.get("producer_source_tree", ""),
                                           "e1b900df838221bdde8edf85412310bd0fe4ae76", "0")
    check(ok, f"c4: the index equals M15-J's apart from the provenance strings and the header digest ({how})")
    for k in ("schema_version", "preparation_algorithm_version", "panel_extension_version", "marker_panel_algorithm_version",
              "identity_gfa_sha256", "model_fingerprint"):
        old_rec = build_record(open(f"{R}/m15e_panel_index_20261021/runs/ex.build_c4/index/c4.build.tsv").read())
        check(rec.get(k) == old_rec.get(k), f"c4: compatibility identity {k} unchanged ({rec.get(k)})")
    cold, ind = f"{d}/c4.cold", f"{d}/c4.index"
    p1, p2 = E.run(E.cmd("c4", cold)), E.run(E.cmd("c4", ind, extra=["--panel-index", idx]))
    check(p1.returncode == p2.returncode == 0, f"c4: cold and indexed exit 0 ({p1.returncode}, {p2.returncode})")
    M15E_X.CHECKS = CHECKS
    M15E_X.check = check
    M15E_X.compare("c4", cold, ind, "c4 cold vs indexed (M15-K)")
    M15J_X.same_outputs(ind, f"{M15J}/a1.ex.index_c4/c4.index", "c4 indexed vs M15-J indexed")
    return 0


def ex_provenance_rebuild(d):
    """A new commit, then `cmake --build` WITHOUT reconfiguring: the next index records the new HEAD."""
    before = head_commit()
    c = subprocess.run(["git", "-C", WT, "commit", "--allow-empty", "-q", "-m",
                        "M15-K ex.provenance_rebuild: an empty commit to show the index records the source revision built "
                        "(no source change)"], capture_output=True, text=True)
    after = head_commit()
    check(c.returncode == 0 and after != before, f"a new commit: {before[:12]} -> {after[:12]}")
    b = subprocess.run(["cmake", "--build", f"{WT}/build", "-j", "4", "--target", "panvar"], capture_output=True, text=True)
    check(b.returncode == 0, "rebuilt without reconfiguring")
    hdr = open(f"{WT}/build/generated/panvar_source_provenance.h").read()
    check(after in hdr and '"clean"' in hdr, f"the generated header names the new HEAD, clean")
    idx = f"{d}/c4.pvpl"
    p = E.run(M15E_X.locus_cmd("c4", ["--build-panel-index", idx]))
    rec = build_record(p.stdout)
    check(p.returncode == 0 and rec.get("producer_source_commit") == after and rec.get("producer_source_tree") == "clean",
          f"the index built after the rebuild records the new HEAD: {rec.get('producer_source_commit')} {rec.get('producer_source_tree')}")
    return 0


if __name__ == "__main__":
    name, d = sys.argv[1], sys.argv[2]
    os.makedirs(d, exist_ok=True)
    globals()["ex_" + name](d)
    bad = [w for ok, w in CHECKS if not ok]
    print(f"{len(CHECKS) - len(bad)} of {len(CHECKS)} checks passed" + (f"; FAILED: {bad}" if bad else ""))
    sys.exit(1 if bad else 0)
