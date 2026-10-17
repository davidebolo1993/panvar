#!/usr/bin/env python3
"""S1-B fixture gates: the implemented sequence_first_v1 outputs against the fixture oracle (s1b_oracle.py), legacy
identity, permutation invariance and exit semantics, on S1-A's fixtures plus one new case.

Fixtures (reads and graphs exactly as S1-A built them; F1-F4 reads reused read-only):
  F1-F4 representation fixtures; F5 partial resolution / genuine tie; F6 lower-ranked projection-ambiguous group;
  F7 provisional (round cap 1), F7full its completed control;
  F9 (new) an UNOBSERVED projection: F1's graph with the (-,R) labelling absent from the panel, donor (R,-)/(R,-), no
     mosaics: the caller never explores (-,R), so it reports sv:1/sv:2 IDENTIFIED with projection_coverage EXPLORED,
     while the complete oracle finds them representation-ambiguous -- why EXPLORED is not completeness.
Per fixture: base/rev/ren/sam variants under sequence_first_v1 (outputs must be byte-identical across variants), plus
the base variant with the schema unset (every legacy output must equal the sequence_first_v1 run's, except the
declared status rows, the one settings row and the added sequence files). F2 base is also run through a persisted
prepared-locus index (development driver: write-panel then loaded) and must equal its cold sequence outputs.
Oracle agreement is required exactly when the oracle finds no unexplored projection; otherwise the difference is
recorded and the caller's report must be a restriction of the oracle's (EXPLORED subset of COMPLETE).

  s1b_fixtures.py <out dir> <panvar> <cyp2d6_block_alleles> <genotype_prepared_locus_driver> [F1,F2,...]
"""
import csv
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
for sub in ("../representation_fixtures", "../cyp2d6_pilot", "../s1a", "."):
    sys.path.insert(0, os.path.join(HERE, sub))
import fixtures  # noqa: E402
import pilot  # noqa: E402
import random  # noqa: E402
import s1a_fixtures  # noqa: E402
import s1b_oracle  # noqa: E402

REF = fixtures.REF
OLD = "/Users/davide.bolognini/panvar/tests/results/representation_fixtures_20261015"
# declared before execution: (legacy exit, v1 exit, v1 prediction, v1 sequence_status), from S1-A's sealed findings
EXPECT = {"F1": (1, 0, "PREDICTION", "RESOLVED"), "F2": (1, 0, "PREDICTION", "RESOLVED"),
          "F3": (0, 0, "PREDICTION", "GENOTYPE_RESOLVED_SEQUENCE_SET"), "F4": (1, 0, "PREDICTION", "RESOLVED"),
          "F5": (1, 0, "PREDICTION", "SEQUENCE_TIED"), "F6": (0, 0, "PREDICTION", "RESOLVED"),
          "F7": (1, 3, "PROVISIONAL", "NOT_DETERMINED"), "F7full": (0, 0, "PREDICTION", "RESOLVED"),
          "F9": (0, 0, "PREDICTION", "RESOLVED")}
SEQUENCE_FILES = ["call.sequence_locus.tsv", "call.sequence_blocks.tsv", "call.representation_intervals.tsv",
                  "call.sequence_alternatives.tsv", "call.sequence_haplotypes.fa"]
# measurements, plus the indexed-union executed-work counters that vary between identical runs at >1 thread (P1/P1.1)
VOLATILE = re.compile(r"seconds|footprint|rss|peak|_bytes|timing|worker|threads_|profile_sha256|indexed_union_executed")


def build(kind):
    if kind in ("F1", "F2", "F3", "F4"):
        nodes, rows_ = fixtures.build(kind)
        return nodes, rows_, None, None, os.path.join(OLD, kind, "reads")
    if kind == "F9":
        nodes, rows_ = fixtures.build("F1")
        named = dict(rows_)
        rm, rr, none = named["S1#1#c"], named["S2#2#c"], named[REF]   # (R,-), (R,R), (-,-)
        rows9 = [(REF, none), ("S1#1#c", rm), ("S1#2#c", rr), ("S2#1#c", rm), ("S2#2#c", none), ("S3#1#c", rr)]
        return nodes, rows9, [rm, rm], None, None
    nodes, rows_, donor, gap = s1a_fixtures.build("F7" if kind.startswith("F7") else kind)
    return nodes, rows_, donor, gap, None


def kv(path):
    out = {}
    if os.path.exists(path):
        for line in open(path).read().splitlines():
            p = line.split("\t", 1)
            if len(p) == 2:
                out[p[0]] = p[1]
    return out


def tsv(path):
    return list(csv.DictReader(open(path), delimiter="\t")) if os.path.exists(path) else []


def norm(text, d):
    return text.replace(d, "<D>")


def legacy_identity(d_legacy, d_v1):
    """Problems ([] when the v1 run's legacy outputs equal the unset run's, with only the declared additions)."""
    problems = []
    a = {f for f in os.listdir(d_legacy) if f.startswith("call.")}   # the caller's outputs (not runner logs/profiles)
    b = {f for f in os.listdir(d_v1) if f.startswith("call.")}
    extra = b - a
    if extra - set(SEQUENCE_FILES) or a - b:
        problems.append(f"file sets: v1 extra {sorted(extra - set(SEQUENCE_FILES))}, missing {sorted(a - b)}")
    for f in sorted(a & b):
        pa, pb = os.path.join(d_legacy, f), os.path.join(d_v1, f)
        if f.endswith("_status.tsv"):
            ka, kb = kv(pa), kv(pb)
            if f == "call.mosaic_status.tsv":
                for k in ("output_schema", "sequence_prediction", "sequence_status", "sequence_completion",
                          "sequence_exit_code"):
                    if k not in kb:
                        problems.append(f"{f}: declared row {k} missing")
                    kb.pop(k, None)
            if set(ka) != set(kb):
                problems.append(f"{f}: rows added/removed {sorted(set(ka) ^ set(kb))}")
            for k, v in ka.items():
                if not VOLATILE.search(k) and norm(v, d_legacy) != norm(kb.get(k, ""), d_v1):
                    problems.append(f"{f}: {k}: {v} vs {kb.get(k)}")
        elif f == "call.settings.tsv":
            tb = open(pb).read()
            row = "output_schema\tsequence_first_v1\tflag --output-schema\n"
            if tb.count(row) != 1:
                problems.append("settings: the one output_schema row")
            if tb.replace(row, "", 1) != open(pa).read():
                problems.append("settings: other rows changed")
        elif f in ("call.repair_memory.tsv", "call.log"):
            continue
        elif f == "call.repair_rounds.tsv":
            cut = lambda p: [l.split("\t")[:28] for l in open(p).read().splitlines()]  # noqa: E731
            if cut(pa) != cut(pb):
                problems.append("repair trajectory differs")
        elif norm(open(pa, errors="replace").read(), d_legacy) != norm(open(pb, errors="replace").read(), d_v1):
            problems.append(f"{f}: differs")
    return problems


def close(x, y):
    return abs(x - y) <= 1e-5 * max(abs(x), abs(y)) + 1e-12


def oracle_check(d, dump):
    """(agreement verdict, notes). Exact agreement when the oracle has no unexplored projection."""
    o = s1b_oracle.oracle(os.path.join(d, "call"), dump)
    loc = kv(os.path.join(d, "call.sequence_locus.tsv"))
    blocks = tsv(os.path.join(d, "call.sequence_blocks.tsv"))
    ivs = [r for r in tsv(os.path.join(d, "call.representation_intervals.tsv")) if r["representative_group"] == "yes"]
    notes = []
    if loc.get("prediction") != o["prediction"]:
        notes.append(f"prediction {loc.get('prediction')} vs oracle {o['prediction']}")
    if o["prediction"] == "NO_PREDICTION":
        return ("AGREE" if not notes else "DISAGREE"), notes, o
    if loc.get("sequence_status") != o["sequence_status"]:
        notes.append(f"sequence_status {loc.get('sequence_status')} vs {o['sequence_status']}")
    if (loc.get("representative_hap1_sha256"), loc.get("representative_hap2_sha256")) != o["rep"]:
        notes.append("representative differs")
    for row, ob in zip(blocks, o["blocks"]):
        if (row["block_id"], row["block_status"]) != ob[:2] or (ob[2] != "NA" and row["genotype"] != ob[2]):
            notes.append(f"block {row['block_id']}: {row['block_status']} {row['genotype']} vs {ob[1]} {ob[2]}")
        elif ob[1] == "IDENTIFIED" and not (close(float(row["error_mass_lower"]), ob[3]) and
                                            close(float(row["error_mass_upper"]), ob[4])):
            notes.append(f"block {row['block_id']}: bounds [{row['error_mass_lower']}, {row['error_mass_upper']}] vs "
                         f"[{ob[3]:.6g}, {ob[4]:.6g}]")
    if len(blocks) != len(o["blocks"]):
        notes.append("block count")
    mine = []
    for r in ivs:
        tuples = sorted(tuple(sorted((tuple(int(x) for x in a.split(",")), tuple(int(x) for x in b.split(",")))))
                        for a, b in (t.split("/") for t in r["tuples_hapA_hapB"].split(";")))
        mine.append((r["interval"], tuples))
    theirs = [(i, [tuple(t) for t in ts]) for i, ts in o["intervals"]]
    if mine != theirs:
        notes.append(f"intervals {[(i, len(t)) for i, t in mine]} vs oracle {[(i, len(t)) for i, t in theirs]}")
    if not notes:
        return "AGREE", notes, o
    return ("EXPLORED_RESTRICTION" if o["projections_unexplored"] > 0 else "DISAGREE"), notes, o


def run(cmd, log, timeout=300):
    with open(log, "w") as f:
        f.write("$ " + " ".join(cmd) + "\n")
        f.flush()
        code, _ = pilot.run_timed(cmd, f, timeout)
    return code


def main(out, panvar, dumper, driver, only=None):
    os.makedirs(out, exist_ok=True)
    base = "profile\tfixture-bounds\ntime_limit_seconds\t300\nfootprint_limit_bytes\t4294967296\n"
    profiles = {"default": base, "round1": base + "repair_rounds\t1\nmax_mosaic_candidates\t0\n",
                "k0": base + "max_mosaic_candidates\t0\n"}
    for name, text in profiles.items():
        open(os.path.join(out, f"profile_{name}.tsv"), "w").write(text)
    plan = [("F1", "default"), ("F2", "default"), ("F3", "default"), ("F4", "default"), ("F5", "default"),
            ("F6", "default"), ("F7", "round1"), ("F7full", "k0"), ("F9", "k0")]
    if only:
        plan = [p for p in plan if p[0] in only]
    summary = []
    for kind, profile in plan:
        d = os.path.join(out, kind)
        os.makedirs(d, exist_ok=True)
        nodes, rows_, donor, gap, reads = build(kind)
        if reads is None:
            reads = os.path.join(d, "reads")
            s1a_fixtures.tiled_reads(nodes, donor, gap, reads)
        perm = list(nodes)
        random.Random(7).shuffle(perm)
        variants = {"base": {}, "rev": {"reverse": True}, "ren": {"rename": {n: str(100 + i) for i, n in enumerate(perm)}},
                    "sam": {"samples": {"S1": "Z1", "S2": "A9", "S3": "M5", "S4": "B2"}}}
        outputs, exits = {}, {}
        for name, opts in variants.items():
            gfa = os.path.join(d, f"{name}.gfa")
            fixtures.write_gfa(gfa, nodes, rows_, **opts)
            fixtures.run([dumper, gfa, REF, "none", os.path.join(d, f"{name}.model")], os.path.join(d, f"{name}.model.log"))
            arms = ["v1", "legacy"] if name == "base" else ["v1"]
            for arm in arms:
                rd = os.path.join(d, f"{name}.{arm}")
                os.makedirs(rd, exist_ok=True)
                args = ["--preset", "short-read", "--profile", os.path.join(out, f"profile_{profile}.tsv"), "--threads", "2",
                        "--graph", gfa, "--reference-path", REF, "--reads1", os.path.join(reads, "reads_1.fq"), "--reads2",
                        os.path.join(reads, "reads_2.fq"), "--library", os.path.join(reads, "library.txt"),
                        "--dosage-mode", "repair", "--out", os.path.join(rd, "call"), "--quiet"]
                if arm == "v1":
                    args += ["--output-schema", "sequence_first_v1"]
                exits[(name, arm)] = run([panvar, "genotype", *args], os.path.join(d, f"{name}.{arm}.log"))
                if arm == "v1":
                    open(os.path.join(d, f"{name}.v1.args.txt"), "w").write("".join(a + "\n" for a in args))
            outputs[name] = {f: norm(open(os.path.join(d, f"{name}.v1", f)).read(), os.path.join(d, f"{name}.v1"))
                             for f in SEQUENCE_FILES if os.path.exists(os.path.join(d, f"{name}.v1", f))}
        verdict, notes, o = oracle_check(os.path.join(d, "base.v1"), os.path.join(d, "base.model"))
        legacy_problems = legacy_identity(os.path.join(d, "base.legacy"), os.path.join(d, "base.v1"))
        invariant = all(outputs[v] == outputs["base"] for v in outputs)
        loc = kv(os.path.join(d, "base.v1", "call.sequence_locus.tsv"))
        want = EXPECT[kind]
        got = (exits[("base", "legacy")], exits[("base", "v1")], loc.get("prediction"), loc.get("sequence_status"))
        index_check = "not_run"
        if kind == "F2":   # persisted prepared-locus index: write-panel, then loaded, against the cold v1 outputs
            idx = os.path.join(d, "F2.panel.pvpl")
            if os.path.exists(idx):
                os.remove(idx)
            ld = os.path.join(d, "base.v1.loaded")
            os.makedirs(ld, exist_ok=True)
            args = open(os.path.join(d, "base.v1.args.txt")).read().replace(os.path.join(d, "base.v1", "call"),
                                                                             os.path.join(ld, "call"))
            open(os.path.join(d, "base.v1.loaded.args.txt"), "w").write(args)
            w = run([driver, "write-panel", str(4 << 30), "0", os.path.join(d, "drv_write"), idx,
                     os.path.join(d, "base.v1.loaded.args.txt")], os.path.join(d, "index_write.log"))
            lcode = run([driver, "loaded", str(4 << 30), "0", os.path.join(d, "drv_load"), idx,
                         os.path.join(d, "base.v1.loaded.args.txt")], os.path.join(d, "index_load.log"))
            same = all(norm(open(os.path.join(ld, f)).read(), ld) == outputs["base"][f] for f in outputs["base"])
            index_check = f"write {w}, load {lcode}, sequence outputs {'identical' if same else 'DIFFERENT'}"
        row = {"fixture": kind, "legacy_exit": got[0], "v1_exit": got[1], "prediction": got[2], "sequence_status": got[3],
               "expected": "/".join(str(x) for x in want), "expectation_met": "yes" if got == want else "NO",
               "oracle": verdict, "oracle_notes": " | ".join(notes) or "none",
               "projections_unexplored": o["projections_unexplored"],
               "legacy_identity": "identical" if not legacy_problems else " | ".join(legacy_problems[:6]),
               "variants_identical": "yes" if invariant else "NO",
               "intervals": loc.get("representation_intervals", "NA"), "decision_set_groups": loc.get("decision_set_groups", "NA"),
               "representative_probability": loc.get("representative_probability", "NA"), "index_route": index_check}
        summary.append(row)
        print("\t".join(str(v) for v in row.values()), flush=True)
    path = os.path.join(out, "fixtures.tsv")
    kept = [r for r in tsv(path) if r["fixture"] not in {k for k, _ in plan}]
    with open(path, "w") as f:
        w = csv.DictWriter(f, fieldnames=list(summary[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(kept + summary)
    bad = [r for r in summary if r["expectation_met"] != "yes" or r["oracle"] == "DISAGREE" or
           r["legacy_identity"] != "identical" or r["variants_identical"] != "yes" or
           (r["index_route"] != "not_run" and "identical" not in r["index_route"].split("outputs ")[-1])]
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2]), os.path.abspath(sys.argv[3]),
                  os.path.abspath(sys.argv[4]), sys.argv[5].split(",") if len(sys.argv) > 5 else None))
