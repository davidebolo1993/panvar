#!/usr/bin/env python3
"""Unit tests for scripts/plot_genotype_results.py (the panvar genotype reporting tool), on small synthetic tables:
aggregation counts and denominators computed by hand, the normalised-excess formula, truth-undefined exclusion, every
refusal, the modality contract (absent: no column; one: figures identical to absent; several: faceted), byte-identical
reruns, directory-mode loading (including pre-headline outputs), the per-coverage n / zero block labels (read from
block_metric_summary, laid out without overlap for 1, 2 and 4 coverages) and the observed truth-allele length ranges
under the block IDs (evaluation only). Usage: genotype_report_test.py [scripts dir]."""
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
from genotype_report import plots as P  # noqa: E402
from genotype_report.common import Refusal, format_length_range, wilson  # noqa: E402

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
        self.assertEqual(r.diploid_length_denominator_bp, "240.0")
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
        with open(os.path.join(self.out("qc"), "report.md")) as f:
            text = f.read()
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
        self.assertTrue(any("diploid length denominator is unavailable" in why for _, why in r.skipped))
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

    # ---------------------------------------------------------------------------------------------- block labels
    def counts(self, r, metric, set_name, row_key="panel", rows=None):
        """{(row, block, coverage): (text, x, colour)} of the n / zero labels the ordered-block figure draws."""
        rows = rows or [(p, p) for p, _ in r.ctx.panels]
        pts = r.tables["block_evaluation"]
        pts = pts if set_name == "all called blocks" else pts[pts.exact_genotype_usable == True]  # noqa: E712
        fig = P.fig_ordered_blocks(pts, r.tables["block_metric_summary"], r.info, metric, dict(A.METRICS)[metric], set_name, rows, row_key,
                                   lengths=r.tables.get("block_lengths"))
        ccol = P.coverage_colors(r.info.coverages)
        w = 0.8 / len(r.info.coverages)
        out = {}
        for ax, (_, rval) in zip(fig.axes, rows):
            for t in [c for c in ax.get_children() if c.get_gid() == "block_count"]:
                x = t.xy[0]
                bi = int(round(x))
                ci = int(round((x - bi) / w + (len(r.info.coverages) - 1) / 2))
                c = r.info.coverages[ci]
                self.assertEqual(t.get_color(), ccol[c])   # drawn in its coverage's colour, at its dodged position
                out[(rval, r.info.block_order[bi], c)] = t.get_text()
        return fig, out

    def test_per_coverage_counts_by_hand(self):
        r = self.full("c")
        _, ex = self.counts(r, "excess", "all called blocks")
        # B: r4 at 2x; r5, r6 at 10x. r6 ctx:0 excess 3, r5 sv:1 excess 5, r4 ctx:2 truth-undefined
        self.assertEqual(ex[("B", "ctx:0", "2")], "n=1(1)")
        self.assertEqual(ex[("B", "ctx:0", "10")], "n=2(1)")
        self.assertEqual(ex[("B", "sv:1", "10")], "n=2(1)")
        self.assertEqual(ex[("B", "ctx:2", "2")], "n=0(0,1)")
        self.assertEqual(ex[("B", "ctx:2", "10")], "n=2(2)")
        # A: r1 at 2x (r2 is a no-call and contributes nothing), r3 at 10x; r1 sv:1 excess 2000
        self.assertEqual(ex[("A", "sv:1", "2")], "n=1(0)")
        self.assertEqual(ex[("A", "sv:1", "10")], "n=1(1)")
        self.assertEqual(len(ex), 2 * 3 * 2)   # rows x blocks x coverages
        # raw distance is excess + 1 here, so no exact sequence distance zero anywhere
        _, di = self.counts(r, "distance", "all called blocks")
        self.assertEqual(di[("B", "ctx:0", "10")], "n=2(0)")
        self.assertTrue(all(v.split("(")[1].split(",")[0].rstrip(")") == "0" for v in di.values()))
        # EXACT_GENOTYPE_USABLE drops r1 sv:1 and r5 sv:1 (GENOTYPE_FILTER LOW_QUALITY)
        _, us = self.counts(r, "excess", "EXACT_GENOTYPE_USABLE")
        self.assertEqual(us[("B", "sv:1", "10")], "n=1(1)")
        self.assertEqual(us[("A", "sv:1", "2")], "n=0(0)")
        self.assertEqual(us[("B", "sv:1", "2")], ex[("B", "sv:1", "2")])
        self.assertEqual(us[("B", "ctx:0", "10")], ex[("B", "ctx:0", "10")])
        # row headings: R = called run instances plotted in the row (A: r1, r3; r2 is a no-call; B: r4, r5, r6)
        fig, _ = self.counts(r, "excess", "all called blocks")
        self.assertEqual([ax.get_title(loc="left") for ax in fig.axes if ax.get_title(loc="left")], ["A (R=2)", "B (R=3)"])
        self.assertEqual(fig._suptitle.get_text(), "Excess over the availability floor (bp) by block \u2014 all called blocks")
        P.plt.close(fig)
        # the labels are the summary table's values, as written to data/
        sm = self.data("c", "block_metric_summary")
        g = sm[(sm.panel == "B") & (sm.block_id == "sv:1") & (sm.coverage == "10") & (sm.metric == "excess")].set_index("set")
        self.assertEqual((g.loc["all called blocks", "n"], g.loc["EXACT_GENOTYPE_USABLE", "n"]), ("2", "1"))

    def test_titles_state_only_the_comparison(self):
        """Every figure title names the metric or comparison and the plotted set; axis, filter and denominator
        explanations live in the captions."""
        titles, save = [], P.save

        def record(fig, stem):
            titles.extend([fig._suptitle.get_text()] if fig._suptitle else [])
            titles.extend(ax.get_title(loc=k) for ax in fig.axes for k in ("left", "center") if ax.get_title(loc=k))
            save(fig, stem)
        P.save = record
        try:
            sw = os.path.join(self.tmp, "meta", "switch.tsv")
            write(sw, [{"subject": "sr4", "switch_class": "one_switch", "switch_before_block_hap1": "2", "switch_before_block_hap2": "none"}],
                  ["subject", "switch_class", "switch_before_block_hap1", "switch_before_block_hap2"])
            r = self.full("t", extra=["--switch-metadata", sw])
        finally:
            P.save = save
        self.assertGreater(len(titles), 20)
        # (policy names such as "GENOTYPE_FILTER == PASS" are the plotted set and may appear)
        banned = ("symlog", "linear in", "Wilson", "USABLE (", "denominator", "cells:", "labels:", "(solid)", "never pooled",
                  "fail-closed", "may carry", "conditional on", "called runs", "units per block", "[", "\n")
        self.assertEqual([t for t in titles if any(b in t for b in banned)], [])
        # the explanations moved to the captions
        caps = " ".join(cap for _, _, cap, _ in r.figures)
        for phrase in ("symlog axis: linear in [0, 1] bp", "95% Wilson", "worst of its units", "may carry several",
                       "never pooled", "labels are retained / all", "n = truth-defined observations plotted; parentheses "
                       "contain the number at zero and, when present, the number of additional truth-undefined observations.",
                       "not the number of unique biological donors"):
            self.assertIn(phrase, caps)

    def test_length_range_format(self):
        cases = {(117, 117): "117 bp", (1, 69): "1\u201369 bp", (12306, 12402): "12.3\u201312.4 kb",
                 (26346, 98187): "26.3\u201398.2 kb", (0, 7670): "0\u20137.7 kb", (0, 69): "0\u201369 bp",
                 (12350, 12350): "12.4 kb", (999, 1000): "0.999\u20131.000 kb",
                 (47951, 48003): "47.95\u201348.00 kb", (117.5, 420): "117.5\u2013420 bp"}
        for (lo, hi), want in cases.items():
            self.assertEqual(format_length_range(lo, hi), want, (lo, hi))

    def test_block_lengths_variable_alleles(self):
        e = pd.read_csv(self.ev, sep="\t", dtype=str, keep_default_na=False)
        lengths = {"ctx:0": {"r1": (1000, 1100), "r3": (1000, 1250), "r4": (990, 1000), "r5": (1000, 1000), "r6": (1100, 1300)},
                   "sv:1": {"r1": (0, 7), "r3": (69, 1), "r4": (7, 7), "r5": (3, 0), "r6": (4, 50)},
                   "ctx:2": {"r1": (12306, 12402), "r3": (12350, 12350), "r5": (12400, 12306), "r6": (12402, 12380)}}
        for i, row in e.iterrows():
            t = lengths[row.block_id].get(row.run_id)
            if t:   # r4 ctx:2 stays truth-undefined (NA) and is excluded
                e.loc[i, ["truth1_bp", "truth2_bp"]] = [str(t[0]), str(t[1])]
        var = os.path.join(self.tmp, "meta", "var_eval.tsv")
        e.to_csv(var, sep="\t", index=False)
        r = run_cli(["--input", self.pkg, "--out-dir", self.out("len"), "--run-metadata", self.meta, "--evaluation-blocks", var])
        bl = self.data("len", "block_lengths").set_index("block_id")
        # ctx:0: 990, 1000 x5, 1100 x2, 1250, 1300 -> median (1000 + 1000) / 2
        self.assertEqual(bl.loc["ctx:0"].tolist(), ["0", "10", "990.0", "1000.0", "1300.0", "1.0\u20131.3 kb"])
        # sv:1: 0 0 1 3 4 7 7 7 50 69 -> median (4 + 7) / 2 = 5.5
        self.assertEqual(bl.loc["sv:1"].tolist(), ["1", "10", "0.0", "5.5", "69.0", "0\u201369 bp"])
        # ctx:2: 8 observations (r4 undefined): 12306 12306 12350 12350 12380 12400 12402 12402 -> median 12365
        self.assertEqual(bl.loc["ctx:2"].tolist(), ["2", "8", "12306.0", "12365.0", "12402.0", "12.3\u201312.4 kb"])
        self.assertEqual(list(self.data("len", "block_lengths").columns),
                         ["block_index", "block_id", "truth_allele_observations", "truth_allele_bp_min", "truth_allele_bp_median",
                          "truth_allele_bp_max", "displayed_length_range"])
        # every model-order block figure names the table and shows two-line labels
        for name in ("eval_block_excess_all", "eval_block_accuracy_heatmap", "eval_block_accuracy_intervals", "eval_filter_by_block"):
            self.assertIn("block_lengths", [d for _, n, _, ds in r.figures if n == name for d in ds], name)
        self.assertEqual(P.block_labels(r.info.block_order, r.tables["block_lengths"]),
                         ["ctx:0\n1.0\u20131.3 kb", "sv:1\n0\u201369 bp", "ctx:2\n12.3\u201312.4 kb"])
        with open(os.path.join(self.out("len"), "report.md")) as f:
            self.assertIn("It is not a reference span", f.read())

    def test_truth_free_has_no_block_lengths(self):
        r = run_cli(["--input", self.pkg, "--out-dir", self.out("qc"), "--run-metadata", self.meta])
        self.assertFalse(os.path.exists(os.path.join(self.out("qc"), "data", "block_lengths.tsv")))
        self.assertNotIn("block_lengths", [d for _, _, _, ds in r.figures for d in ds])
        with open(os.path.join(self.out("qc"), "report.md")) as f:
            text = f.read()
        self.assertNotIn("truth-allele", text)
        self.assertNotIn("block_lengths", text)
        self.assertEqual(P.block_labels(r.info.block_order, r.tables.get("block_lengths")), BLOCKS)
        # an evaluation without truth lengths has none either
        e = pd.read_csv(self.ev, sep="\t", dtype=str)
        short = os.path.join(self.tmp, "meta", "short_eval.tsv")
        e.drop(columns=["truth1_bp", "truth2_bp", "called1_bp", "called2_bp"]).to_csv(short, sep="\t", index=False)
        run_cli(["--input", self.pkg, "--out-dir", self.out("nolen"), "--run-metadata", self.meta, "--evaluation-blocks", short])
        self.assertFalse(os.path.exists(os.path.join(self.out("nolen"), "data", "block_lengths.tsv")))

    def assert_no_overlap(self, fig):
        fig.canvas.draw()
        rd = fig.canvas.get_renderer()
        boxes = []
        for ax in fig.axes:
            labels = [c.get_window_extent(rd) for c in ax.get_children() if c.get_gid() == "block_count"]
            others = [ax.title.get_window_extent(rd)] + ([fig._suptitle.get_window_extent(rd)] if fig._suptitle else [])
            for i, a in enumerate(labels):
                self.assertTrue(fig.bbox.contains(a.x0, a.y0) and fig.bbox.contains(a.x1, a.y1), "label outside the figure")
                for b in labels[i + 1:] + others:
                    self.assertFalse(a.overlaps(b), "overlapping labels")
            boxes += labels
            ticks = [t.get_window_extent(rd) for t in ax.get_xticklabels() if t.get_text()]
            for i, a in enumerate(ticks):
                for b in ticks[i + 1:]:
                    self.assertFalse(a.overlaps(b), "overlapping block labels")
        return len(boxes)

    def test_count_label_layout_1_2_4_coverages(self):
        m = pd.read_csv(self.meta, sep="\t", dtype=str)
        for covs in (["30"] * 6, None, ["1", "2", "5", "30", "5", "1"]):
            path = self.meta
            if covs:
                path = os.path.join(self.tmp, "meta", f"cov{len(set(covs))}.tsv")
                m.assign(coverage=covs).to_csv(path, sep="\t", index=False)
            name = f"lay{len(set(covs)) if covs else 2}"
            r = self.full(name, path)
            fig, ex = self.counts(r, "excess", "all called blocks")
            present = {(e.panel, e.coverage) for e in r.tables["block_metric_summary"].itertuples()}
            self.assertEqual(len(ex), len(present) * len(BLOCKS), name)
            self.assertEqual(self.assert_no_overlap(fig), len(ex), name)
            P.plt.close(fig)
            # the widest counts a cohort could have: every label at five digits, with an undefined line
            sm = r.tables["block_metric_summary"].copy()
            sm[["n", "zero", "truth_undefined"]] = 88888
            fig = P.fig_ordered_blocks(r.tables["block_evaluation"], sm, r.info, "excess", "excess", "all called blocks",
                                       [(p, p) for p, _ in r.ctx.panels], "panel", lengths=r.tables["block_lengths"])
            self.assertEqual(self.assert_no_overlap(fig), len(ex), name)
            P.plt.close(fig)
        # four coverages in the data: the per-coverage labels alternate between two heights
        self.assertEqual(P.count_layout(["n=88888(88888,88888)"], 4)[0], 2)
        self.assertEqual(P.count_layout(["n=1(1)"], 1)[0], 1)


if __name__ == "__main__":
    unittest.main(verbosity=1)
