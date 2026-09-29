#!/usr/bin/env python3
"""Unit tests for scripts/plot_genotype_results.py (the panvar genotype reporting tool), on small synthetic tables:
aggregation counts and denominators computed by hand, the normalised-excess formula, truth-undefined exclusion, every
refusal, the modality contract (absent: no column; one: figures identical to absent; several: faceted), byte-identical
reruns, and directory-mode loading (including pre-headline outputs). Usage: genotype_report_test.py [scripts dir]."""
import hashlib
import os
import shutil
import sys
import tempfile
import unittest

import pandas as pd

SCRIPTS = os.path.abspath(sys.argv.pop(1)) if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts")
sys.path.insert(0, SCRIPTS)
import plot_genotype_results as cli  # noqa: E402
from genotype_report import aggregate as A  # noqa: E402
from genotype_report import load as L  # noqa: E402
from genotype_report.common import Refusal, wilson  # noqa: E402

BLOCKS = ["ctx:0", "sv:1", "ctx:2"]


def write(path, rows, header):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join(str(r[h]) for h in header) + "\n")


def synthetic(root):
    """Six runs: experiments A and B at 2x and 10x, A-2x has a no-call and a genotype-equivalent call. Returns the
    package, run metadata and evaluation paths."""
    runs, blocks, meta, ev = [], [], [], []
    spec = [("r1", "A", "2", "strict_call"), ("r2", "A", "2", "no_call"), ("r3", "A", "10", "genotype_equivalent_call"),
            ("r4", "B", "2", "strict_call"), ("r5", "B", "10", "strict_call"), ("r6", "B", "10", "strict_call")]
    for rid, exp, cov, cc in spec:
        runs.append({"run_id": rid, "run_status": "COMPLETE" if cc != "no_call" else "REPAIR_INCOMPLETE", "call_class": cc,
                     "no_call_reason": "none" if cc != "no_call" else "tie_branching_incomplete", "phase_support": "PRIOR_ONLY",
                     "genotype_filter": "PASS", "dosage_filter": "NOT_EVALUATED", "phase_filter": "LOW_QUALITY",
                     "seconds_total": "10.5", "threads_requested": "8", "settings_digest": "settings:abc", "preset": "short-read"})
        meta.append({"run_id": rid, "experiment": exp, "coverage": cov, "subject": "s" + rid})
        if cc == "no_call":
            continue
        for i, b in enumerate(BLOCKS):
            gf = "LOW_QUALITY" if (rid, b) in (("r1", "sv:1"), ("r5", "sv:1")) else "PASS"
            df = "PASS" if b == "sv:1" else "NOT_APPLICABLE"
            blocks.append({"run_id": rid, "block_index": i, "block_id": b, "zygosity": "HET", "conditional_GQ": str(10 * (i + 1)),
                           "genotype_filter": gf, "genotype_filter_reasons": "none", "dosage_filter": df,
                           "dosage_filter_reasons": "none", "dosage_status": "CONSISTENT" if b == "sv:1" else "NOT_SENSITIVE"})
            # truth: r4 ctx:2 is undefined; excess 0 except r1 sv:1 (2000, catastrophic), r5 sv:1 (5) and r6 ctx:0 (3)
            excess = {("r1", "sv:1"): 2000, ("r5", "sv:1"): 5, ("r6", "ctx:0"): 3}.get((rid, b), 0)
            undef = (rid, b) == ("r4", "ctx:2")
            ev.append({"run_id": rid, "block_id": b, "unordered_distance": "NA" if undef else excess + 1,
                       "availability_floor": "NA" if undef else 1, "excess": "NA" if undef else excess,
                       "truth1_bp": "NA" if undef else 100, "truth2_bp": "NA" if undef else 120,
                       "called1_bp": 90, "called2_bp": 150})
    pkg = os.path.join(root, "pkg")
    write(os.path.join(pkg, "runs.tsv"), runs, list(runs[0]))
    write(os.path.join(pkg, "blocks.tsv"), blocks, list(blocks[0]))
    write(os.path.join(root, "meta", "meta.tsv"), meta, list(meta[0]))
    write(os.path.join(root, "meta", "eval.tsv"), ev, list(ev[0]))
    return pkg, os.path.join(root, "meta", "meta.tsv"), os.path.join(root, "meta", "eval.tsv")


def run_cli(argv):
    args = cli.parse(argv)
    from genotype_report.pipeline import run
    return run(args, ["plot_genotype_results.py"] + argv)


def sha(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def digest_tree(d, sub):
    return {f: sha(os.path.join(d, sub, f)) for f in sorted(os.listdir(os.path.join(d, sub)))}


class Report(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="genotype_report_test_")
        self.pkg, self.meta, self.ev = synthetic(self.tmp)

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def out(self, name):
        return os.path.join(self.tmp, "out", name)

    def full(self, name, meta=None, extra=()):
        return run_cli(["--input", self.pkg, "--out-dir", self.out(name), "--run-metadata", meta or self.meta,
                        "--evaluation-blocks", self.ev] + list(extra))

    def data(self, name, table):
        return pd.read_csv(os.path.join(self.out(name), "data", table + ".tsv"), sep="\t", dtype=str, keep_default_na=False)

    # ---------------------------------------------------------------------------------------------- aggregation
    def test_callability_counts_and_denominators(self):
        self.full("a")
        c = self.data("a", "callability")
        a2 = c[(c.experiment == "A") & (c.coverage == "2")].set_index("call_class")
        self.assertEqual(a2.loc["no_call", "no_call_reason"], "tie_branching_incomplete")
        self.assertEqual((a2.loc["no_call", "runs"], a2.loc["no_call", "denominator"]), ("1", "2"))
        self.assertEqual(c[(c.experiment == "B") & (c.coverage == "10")].runs.tolist(), ["2"])

    def test_block_accuracy_excludes_undefined_truth(self):
        self.full("a")
        acc = self.data("a", "block_accuracy")
        # B, ctx:2, 2x: the only run (r4) has undefined truth -> denominator 0, counted as truth_undefined
        r = acc[(acc.panel == "B") & (acc.block_id == "ctx:2") & (acc.coverage == "2") & (acc["set"] == "all called blocks")].iloc[0]
        self.assertEqual((r.denominator, r.truth_undefined), ("0", "1"))
        # B, sv:1, 10x: r5 excess 5, r6 excess 0 -> 1/2
        r = acc[(acc.panel == "B") & (acc.block_id == "sv:1") & (acc.coverage == "10") & (acc["set"] == "all called blocks")].iloc[0]
        self.assertEqual((r.floor_attained, r.denominator), ("1", "2"))
        lo, hi = wilson(1, 2)
        self.assertAlmostEqual(float(r.ci_low), lo)
        self.assertAlmostEqual(float(r.ci_high), hi)

    def test_normalized_excess_formula(self):
        self.full("a")
        e = self.data("a", "block_evaluation")
        r = e[(e.run_id == "r5") & (e.block_id == "sv:1")].iloc[0]
        # excess 5 / max(100 + 120, 90 + 150) = 5 / 240
        self.assertEqual(r.alignment_denominator_bp, "240.0")
        self.assertEqual(float(r.normalized_excess), 5 / 240)   # written and read back exactly
        u = e[(e.run_id == "r4") & (e.block_id == "ctx:2")].iloc[0]
        self.assertEqual((u.truth_defined, u.normalized_excess, u.excess), ("False", "NA", "NA"))

    def test_filter_performance_by_hand(self):
        self.full("a")
        fp = self.data("a", "filter_performance")
        g = fp[(fp.policy == "GENOTYPE_FILTER == PASS") & (fp.stratum == "all")].iloc[0]
        # 14 truth-defined blocks (15 reported, r4 ctx:2 undefined); incorrect: r1 sv:1, r5 sv:1 (both filtered LOW_QUALITY),
        # r6 ctx:0 (retained); catastrophic: r1 sv:1 (filtered)
        self.assertEqual([g.units, g.retained_correct, g.retained_incorrect, g.filtered_correct, g.filtered_incorrect],
                         ["14", "11", "1", "0", "2"])
        self.assertEqual([g.catastrophic_retained, g.catastrophic_filtered], ["0", "1"])
        self.assertAlmostEqual(float(g.bad_call_detection_sensitivity), 2 / 3)
        e = fp[(fp.policy == "EXACT_GENOTYPE_USABLE") & (fp.stratum == "all")].iloc[0]
        self.assertEqual(e.units, "14")

    def test_exact_genotype_usable_rule(self):
        self.full("a")
        e = self.data("a", "block_evaluation")
        rule = (e.genotype_filter == "PASS") & e.dosage_filter.isin(["PASS", "NOT_APPLICABLE"])
        self.assertTrue(((e.exact_genotype_usable == "True") == rule).all())

    # ---------------------------------------------------------------------------------------------- refusals
    def test_accuracy_refused_without_truth(self):
        r = run_cli(["--input", self.pkg, "--out-dir", self.out("qc"), "--run-metadata", self.meta])
        self.assertFalse(any(n.startswith("eval_") for _, n, _, _ in r.figures))
        self.assertTrue(any("truth is never inferred" in why for _, why in r.skipped))
        self.assertFalse(os.path.exists(os.path.join(self.out("qc"), "data", "block_evaluation.tsv")))
        text = open(os.path.join(self.out("qc"), "report.md")).read()
        self.assertIn("truth-free QC report", text)

    def test_refusals(self):
        bad_eval = os.path.join(self.tmp, "meta", "bad_eval.tsv")
        write(bad_eval, [{"run_id": "r9", "block_id": "sv:1", "unordered_distance": 0, "availability_floor": 0, "excess": 0}],
              ["run_id", "block_id", "unordered_distance", "availability_floor", "excess"])
        with self.assertRaisesRegex(Refusal, "runs not in the input"):
            run_cli(["--input", self.pkg, "--out-dir", self.out("x1"), "--evaluation-blocks", bad_eval])
        with self.assertRaisesRegex(Refusal, "overlaps the input"):
            run_cli(["--input", self.pkg, "--out-dir", os.path.join(self.pkg, "report")])
        with self.assertRaisesRegex(Refusal, "exactly once"):
            run_cli(["--input", self.pkg, "--out-dir", self.out("x2"), "--run-metadata", self.meta, "--panel-groups", "G=A"])
        # mixed block models
        b = pd.read_csv(os.path.join(self.pkg, "blocks.tsv"), sep="\t", dtype=str)
        b.loc[b.run_id == "r5", "block_id"] = b.loc[b.run_id == "r5", "block_id"].str.replace("ctx:2", "ctx:9")
        b.to_csv(os.path.join(self.pkg, "blocks.tsv"), sep="\t", index=False)
        with self.assertRaisesRegex(Refusal, "disagree on the block model"):
            run_cli(["--input", self.pkg, "--out-dir", self.out("x3")])

    def test_malformed_schema(self):
        r = pd.read_csv(os.path.join(self.pkg, "runs.tsv"), sep="\t", dtype=str)
        r.drop(columns=["call_class"]).to_csv(os.path.join(self.pkg, "runs.tsv"), sep="\t", index=False)
        with self.assertRaisesRegex(Refusal, "missing required column"):
            run_cli(["--input", self.pkg, "--out-dir", self.out("x4")])

    def test_normalization_refused_without_lengths(self):
        e = pd.read_csv(self.ev, sep="\t", dtype=str)
        short = os.path.join(self.tmp, "meta", "short_eval.tsv")
        e.drop(columns=["truth1_bp", "truth2_bp", "called1_bp", "called2_bp"]).to_csv(short, sep="\t", index=False)
        r = run_cli(["--input", self.pkg, "--out-dir", self.out("n"), "--run-metadata", self.meta, "--evaluation-blocks", short])
        self.assertTrue(any("normalisation denominator is unavailable" in why for _, why in r.skipped))
        self.assertFalse(any("normalized_excess" in n for _, n, _, _ in r.figures))

    # ---------------------------------------------------------------------------------------------- modality
    def test_modality_contract(self):
        self.full("none")
        self.assertNotIn("modality", self.data("none", "callability").columns)
        m = pd.read_csv(self.meta, sep="\t", dtype=str)
        one = os.path.join(self.tmp, "meta", "one.tsv")
        m.assign(modality="short-read").to_csv(one, sep="\t", index=False)
        self.full("one", one)
        # every aggregate carries the modality; dropping it gives exactly the tables of the run without it
        for f in sorted(os.listdir(os.path.join(self.out("none"), "data"))):
            a = self.data("none", f[:-4])
            b = self.data("one", f[:-4])
            if "modality" in b.columns:
                self.assertTrue((b.modality == "short-read").all(), f)
                b = b.drop(columns=["modality"])
            self.assertTrue(a.equals(b), f"{f} differs apart from the modality column")
        # figures: with one modality every figure is byte-identical to the run without the column
        self.assertEqual(digest_tree(self.out("none"), "figures"), digest_tree(self.out("one"), "figures"))
        # two modalities: every figure is drawn once per modality
        two = os.path.join(self.tmp, "meta", "two.tsv")
        m.assign(modality=["short-read", "short-read", "short-read", "long-read", "long-read", "long-read"]).to_csv(two, sep="\t", index=False)
        r = self.full("two", two)
        names = [n for _, n, _, _ in r.figures]
        self.assertIn("qc_callability.short-read", names)
        self.assertIn("qc_callability.long-read", names)
        self.assertNotIn("qc_callability", names)
        # an empty modality is refused, never assumed
        empty = os.path.join(self.tmp, "meta", "empty.tsv")
        m.assign(modality=["short-read", "", "short-read", "short-read", "short-read", "short-read"]).to_csv(empty, sep="\t", index=False)
        with self.assertRaisesRegex(Refusal, "empty modality"):
            self.full("empty", empty)

    # ---------------------------------------------------------------------------------------------- determinism
    def test_rerun_is_byte_identical(self):
        self.full("r1")
        self.full("r2")
        for sub in ("data", "figures"):
            self.assertEqual(digest_tree(self.out("r1"), sub), digest_tree(self.out("r2"), sub), sub)
        self.assertEqual(sha(os.path.join(self.out("r1"), "report.md")), sha(os.path.join(self.out("r2"), "report.md")))

    def test_inputs_are_not_modified(self):
        before = {p: sha(p) for p in (self.meta, self.ev, os.path.join(self.pkg, "runs.tsv"), os.path.join(self.pkg, "blocks.tsv"))}
        self.full("i")
        after = {p: sha(p) for p in before}
        self.assertEqual(before, after)
        self.assertEqual(sorted(os.listdir(self.pkg)), ["blocks.tsv", "runs.tsv"])

    # ---------------------------------------------------------------------------------------------- directory mode
    def test_directory_mode_and_legacy_columns(self):
        d = os.path.join(self.tmp, "dir")
        for rid, legacy in (("x/a", False), ("x/b", True)):
            p = os.path.join(d, rid)
            write(p + ".locus_genotype.tsv", [{"run_status": "COMPLETE", "call_class": "strict_call", "no_call_reason": "none"}],
                  ["run_status", "call_class", "no_call_reason"])
            with open(p + ".mosaic_status.tsv", "w") as f:
                f.write("key\tvalue\nstatus\tCOMPLETE\nseconds_total\t3.5\nthreads_requested\t8\n")
            rows = [{"block_id": b, "zygosity": "HET", "filter": "PASS", "dosage_filter": "PASS"} for b in BLOCKS]
            write(p + ".block_genotypes.tsv", rows, ["block_id", "zygosity", "filter", "dosage_filter"])
            # an original prefix beside it is not a run
            write(p + ".original.locus_genotype.tsv", [{"run_status": "X", "call_class": "strict_call"}], ["run_status", "call_class"])
        runs, blocks, _ = L.load_directory(d)
        self.assertEqual(runs.run_id.tolist(), ["x/a", "x/b"])
        # pre-headline outputs: dosage_filter is the dosage-stage component, never read as the headline filter
        self.assertTrue((blocks.dosage_filter == "NA").all())
        self.assertTrue((blocks.genotype_filter == "NA").all())
        r = run_cli(["--input", d, "--out-dir", self.out("dir")])
        self.assertTrue(any("headline_filters_v2" in why for _, why in r.skipped))

    def test_switch_metadata_and_order(self):
        m = pd.read_csv(self.meta, sep="\t", dtype=str)
        sw = os.path.join(self.tmp, "meta", "switch.tsv")
        write(sw, [{"subject": "sr4", "switch_class": "one_switch", "switch_before_block_hap1": "2", "switch_before_block_hap2": "none"},
                   {"subject": "sr5", "switch_class": "control", "switch_before_block_hap1": "none", "switch_before_block_hap2": "none"}],
              ["subject", "switch_class", "switch_before_block_hap1", "switch_before_block_hap2"])
        r = self.full("sw", extra=["--switch-metadata", sw])
        mp = self.data("sw", "switch_map")
        self.assertEqual(mp.switch_class.iloc[0], "control")   # fewest switches first
        row = mp[(mp.plan == "sr4") & (mp.hap1_switch == "True")]
        self.assertEqual(row.boundary.tolist(), ["sv:1>ctx:2"])   # a switch before block 2 sits between blocks 1 and 2
        bad = os.path.join(self.tmp, "meta", "switch_bad.tsv")
        write(bad, [{"subject": "sr4", "switch_class": "x", "switch_before_block_hap1": "3", "switch_before_block_hap2": "none"}],
              ["subject", "switch_class", "switch_before_block_hap1", "switch_before_block_hap2"])
        with self.assertRaisesRegex(Refusal, "outside the 3-block model"):
            self.full("swb", extra=["--switch-metadata", bad])


if __name__ == "__main__":
    unittest.main(verbosity=1)
