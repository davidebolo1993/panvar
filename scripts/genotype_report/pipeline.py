"""Orchestration: load and validate every input, compute every aggregate (written to data/ before any figure is drawn),
draw the figures whose inputs are present (recording a reason for every skipped one), and write report.md and
manifest.tsv. The input directories are only read.
"""
import os
import platform
import re
import subprocess
import sys
from dataclasses import dataclass

import matplotlib
import numpy as np
import pandas as pd

from . import aggregate as A
from . import load as L
from . import plots as P
from .common import (BLOCK_LENGTH_CAPTION, CALLED, DOSAGE_CONTENT_USABLE_RULE, EXACT_GENOTYPE_USABLE_RULE,
                     NORMALIZED_EXCESS_FORMULA, PHASED_GENOTYPE_USABLE_RULE, Refusal, sha256_file, write_tsv)

TOOL = "panvar genotype report"
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


@dataclass
class Info:
    experiments: list
    coverages: list
    panels: list
    block_order: list


def _inside(a, b):
    a, b = os.path.realpath(a), os.path.realpath(b)
    return a == b or a.startswith(b.rstrip(os.sep) + os.sep)


def _slug(text):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", text).strip("_") or "x"


def count_caption(metric, table):
    zero = "exact sequence distance zero" if metric == "distance" else "availability floor attained"
    return (f"Above each block, one label per coverage at that coverage's dodged position and in its colour: n = the "
            "observations in that row, block and coverage that are called, evaluated, truth-defined, retained by the "
            f"plotted set and non-missing for the metric; zero = those among n whose plotted value is exactly 0 ({zero}). "
            "Truth-undefined observations are excluded and, when present, counted as undefined=k for that coverage; "
            f"no-calls contribute no block observations. n and zero are coverage-specific, read from data/{table}.tsv.")


class Report:
    def __init__(self, args, argv):
        self.args, self.argv = args, argv
        self.out = os.path.abspath(args.out_dir)
        self.tables, self.figures, self.skipped, self.sources = {}, [], [], []

    # ---------------------------------------------------------------------------------------------- inputs
    def load(self):
        a = self.args
        inputs = [a.input] + [p for p in (a.run_metadata, a.evaluation_blocks, a.phase_boundaries, a.phase_calls,
                                          a.phase_projected, a.switch_metadata, a.dosage_content_evaluation,
                                          a.representative_diagnostics) if p]
        for p in inputs:
            if not os.path.exists(p):
                raise Refusal(f"input not found: {p}")
            d = p if os.path.isdir(p) else os.path.dirname(os.path.abspath(p))
            if _inside(self.out, d) or _inside(d, self.out):
                raise Refusal(f"--out-dir {self.out} overlaps the input location {d}: inputs are never written to")
        if os.path.exists(self.out) and os.listdir(self.out):
            raise Refusal(f"--out-dir {self.out} exists and is not empty")
        kind = a.input_kind
        if kind == "auto":
            kind = "package" if all(os.path.isfile(os.path.join(a.input, f)) for f in ("runs.tsv", "blocks.tsv")) else "directory"
        self.input_kind = kind
        runs, blocks, src = (L.load_package if kind == "package" else L.load_directory)(a.input)
        self.sources += src
        blocks, order = L.validate(runs, blocks)
        runs, meta = L.load_run_metadata(a.run_metadata, runs)
        if a.run_metadata:
            self.sources.append(a.run_metadata)
        ctx = A.build_context(runs, blocks, order, a.panel_groups, a.catastrophic_bp)
        if a.evaluation_blocks:
            ctx.evaluation, ctx.has_normalization, src = L.load_evaluation(a.evaluation_blocks, ctx.runs, ctx.blocks)
            self.sources += src
        if a.dosage_content_evaluation:
            if not a.evaluation_blocks:
                raise Refusal("--dosage-content-evaluation needs --evaluation-blocks (the content is judged on evaluated blocks)")
            ctx.content = L.load_dosage_content(a.dosage_content_evaluation, ctx.blocks)
            self.sources.append(a.dosage_content_evaluation)
        if a.representative_diagnostics:
            ctx.representatives = L.load_representatives(a.representative_diagnostics, ctx.blocks)
            self.sources.append(a.representative_diagnostics)
        if a.switch_metadata:
            ctx.switch_plans, ctx.switch_key = L.load_switch_metadata(a.switch_metadata, ctx.runs, len(order))
            self.sources.append(a.switch_metadata)
        if a.phase_boundaries:
            ctx.phase_boundaries = L.load_phase_boundaries(a.phase_boundaries, ctx.runs)
            self.sources.append(a.phase_boundaries)
        if a.phase_calls:
            ctx.phase_calls = L.load_phase_calls(a.phase_calls, ctx.runs)
            self.sources.append(a.phase_calls)
        if a.phase_projected:
            ctx.phase_projected = L.load_phase_projected(a.phase_projected, ctx.runs)
            self.sources.append(a.phase_projected)
        self.ctx = ctx
        self.info = Info(ctx.experiments, ctx.coverages, ctx.panels, ctx.block_order)

    # ---------------------------------------------------------------------------------------------- aggregates
    def table(self, name, frame):
        self.tables[name] = frame
        write_tsv(frame, os.path.join(self.out, "data", name + ".tsv"))
        return frame

    def skip(self, figure, reason):
        self.skipped.append((figure, reason))

    def aggregate(self):
        c = self.ctx
        self.table("runs", c.runs.drop(columns=[x for x in ("_meta_order",) if x in c.runs.columns]))
        self.table("blocks", A.sort_frame(c.blocks.drop(columns=["conditional_GQ_value"]), c, ["run_id"]))
        self.table("callability", A.callability(c))
        if c.has_headline:
            self.table("headline_block_status", A.headline_block_status(c))
            self.table("headline_locus_status", A.headline_locus_status(c))
            self.table("headline_block_reasons", A.headline_block_reasons(c))
        self.table("gq_values", A.gq_values(c))
        self.table("gq_summary", A.gq_summary(c))
        self.table("genotype_support", A.genotype_support(c))
        self.table("phase_support", A.phase_support(c))
        self.table("dosage_states", A.dosage_states(c))
        self.table("repair_rounds", A.repair_rounds(c))
        self.table("runtime", A.runtime_points(c))
        self.table("runtime_summary", A.runtime_summary(c))
        self.table("runtime_groups", A.runtime_groups(c))
        self.table("library", A.library(c))
        self.table("provenance", A.provenance(c))
        if c.evaluation is not None:
            pts = self.table("block_evaluation", A.evaluated_blocks(c))
            lengths = A.block_lengths(c)
            if lengths is not None:
                self.table("block_lengths", lengths)
            self.table("block_metric_summary", A.block_metric_summary(c, pts))
            self.table("block_accuracy", A.block_accuracy(c, pts))
            self.table("gq_diagnostic", A.gq_diagnostic(c, pts))
            if c.has_headline:
                self.table("filter_performance", A.filter_performance(c, pts))
            if c.switch_plans is not None:
                self.table("mosaic_block_metric_summary", A.block_metric_summary(c, pts[pts.switch_class.fillna("") != ""], by="switch_class"))
                self.table("mosaic_excess_summary", A.mosaic_excess_summary(c, pts))
            if c.content is not None and c.has_headline:
                self.table("dosage_content_usable", A.dosage_content_usable(c, pts))
                self.table("dosage_content_vs_exact", A.dosage_content_vs_exact(c, pts))
        if c.switch_plans is not None:
            self.table("switch_map", A.switch_map(c))
            self.table("mosaic_callability", A.mosaic_callability(c))
        if c.representatives is not None:
            self.table("dosage_representatives", A.representatives(c))
        if c.phase_boundaries is not None:
            self.table("phase_boundary_status", A.phase_boundary_status(c))
            self.table("phase_accuracy", A.phase_accuracy(c))
            self.table("phase_by_boundary_type", A.phase_by_boundary_type(c))
            if c.switch_plans is not None:
                self.table("phase_by_switch_class", A.phase_by_switch_class(c))
        if c.phase_calls is not None:
            self.table("phase_locus_status", A.phase_locus_status(c))
        if c.phase_projected is not None:
            self.table("phase_projected", A.phase_projected(c))
        self.table("strict_phase_support", A.strict_phase_support(c))

    # ---------------------------------------------------------------------------------------------- figures
    def figure(self, name, section, data, draw, caption):
        """Draw once, or once per modality when more than one is present (the figure name gains .<modality>)."""
        c = self.ctx
        facets = [(None, None)] if len(c.modalities) <= 1 else [(m, "." + _slug(m)) for m in c.modalities]
        for modality, suffix in facets:
            frames = []
            for d in data:
                t = self.tables[d]
                if modality is not None and "modality" in t.columns:
                    t = t[t.modality == modality]
                frames.append(t)
            fig = draw(*frames)
            if modality is not None:
                fig.suptitle((fig._suptitle.get_text() + " | " if fig._suptitle else "") + f"modality: {modality}", fontsize=9)
            fname = name + (suffix or "")
            P.save(fig, os.path.join(self.out, "figures", fname))
            self.figures.append((section, fname, caption + (f" Modality: {modality}." if modality else ""), data))

    def _lengths(self):
        """The block_lengths table as extra figure data, and the caption sentence, when it exists (evaluation with
        truth lengths); nothing otherwise, so a truth-free report shows no block-length information."""
        return (["block_lengths"], " " + BLOCK_LENGTH_CAPTION) if "block_lengths" in self.tables else ([], "")

    def draw(self):
        c, i = self.ctx, self.info
        T = self.tables
        os.makedirs(os.path.join(self.out, "figures"), exist_ok=True)
        # ---- QC (truth-free)
        self.figure("qc_callability", "QC", ["callability"], lambda t: P.fig_callability(t, i),
                    "Runs by call class, no-calls split by their recorded reason; n is the runs in the cell. Truth-free.")
        if c.has_headline:
            self.figure("qc_headline_block_status", "QC", ["headline_block_status"],
                        lambda t: P.fig_status_panels(t, i, ["GENOTYPE_FILTER", "DOSAGE_FILTER"], "blocks",
                                                      "Headline filters per reported block"),
                        "Block-level GENOTYPE_FILTER and DOSAGE_FILTER of called runs (no-call runs have no reported blocks; "
                        "see qc_callability). A status is the filter's verdict on the evidence, not an accuracy label.")
            self.figure("qc_headline_locus_status", "QC", ["headline_locus_status"],
                        lambda t: P.fig_status_panels(t, i, ["GENOTYPE_FILTER", "DOSAGE_FILTER", "PHASE_FILTER"], "runs",
                                                      "Headline filters per locus (fail-closed: the worst unit)"),
                        "Locus-level headline filters: each is the worst of its units, so one non-PASS block makes the locus "
                        "non-PASS. Truth-free.")
            self.figure("qc_headline_reasons", "QC", ["headline_block_reasons"], lambda t: P.fig_reasons(t, i),
                        "Reasons attached to block-level headline verdicts, by experiment. Truth-free.")
        else:
            for f in ("qc_headline_block_status", "qc_headline_locus_status", "qc_headline_reasons"):
                self.skip(f, "the outputs carry no headline filters (produced before headline_filters_v2)")
        self.figure("qc_gq_distribution", "QC", ["gq_values"], lambda t: P.fig_gq(t, i),
                    "Conditional GQ of reported blocks. GQ is conditional on the scored candidate set and is not globally "
                    "calibrated. Truth-free: no accuracy is implied.")
        self.figure("qc_support", "QC", ["genotype_support", "phase_support"], lambda g, p: P.fig_support(g, p, i),
                    "Genotype support of reported blocks and phase support of runs by call class. Truth-free.")
        self.figure("qc_dosage_states", "QC", ["dosage_states"], lambda t: P.fig_dosage_states(t, i, A.DOSAGE_OUTCOMES),
                    "Dosage-stage state of reported blocks: complete repair (applied, available in audit mode, or refused by a "
                    "guard), content-only correction, triggered unrepaired, SUSPECT, UNDERPOWERED, CONSISTENT, not "
                    "dosage-sensitive, or not evaluated (dosage mode off). Truth-free.")
        if T["repair_rounds"].repair_rounds_run.notna().any():
            self.figure("qc_repair_rounds", "QC", ["repair_rounds"], lambda t: P.fig_repair_rounds(t, i),
                        "Adaptive repair rounds run per run, coloured by call class. Truth-free.")
        else:
            self.skip("qc_repair_rounds", "no run records repair_rounds_run")
        if T["runtime"].seconds_total.notna().any():
            groups = "; ".join(f"{r.compat_group} ({r.runs} runs, {r.distinct_graphs} graph digest(s))"
                               for r in T["runtime_groups"].itertuples())
            self.figure("qc_runtime", "QC", ["runtime"], lambda t: P.fig_runtime(t, i),
                        "Wall time, peak physical footprint, repair rounds and dosage-stage time per run. Colour is the "
                        "compatibility group (settings digest, threads, execution group, and modality when supplied); groups "
                        "are never pooled. The genotyper records no host or binary identity: execution_group, from the run "
                        "metadata, must distinguish machine, environment and comparable workload conditions whenever runtime "
                        f"is compared. Groups: {groups}.")
        else:
            self.skip("qc_runtime", "no run records seconds_total")
        self.figure("qc_library", "QC", ["library"], lambda t: P.fig_library(t, i),
                    "Library source (declared file or auto estimate) and library-sensitivity status per run. Truth-free.")
        self.figure("qc_provenance", "QC", ["provenance"], lambda t: P.fig_provenance(t, i),
                    "Runs per preset (and whether it was the default), profile and settings digest. Truth-free.")
        # ---- evaluation
        if c.evaluation is None:
            for f in ("eval_block_*", "eval_block_accuracy_*", "eval_filter_*", "eval_gq_diagnostic", "mosaic_block_*",
                      "mosaic_excess_by_class", "dosage_content_*"):
                self.skip(f, "no --evaluation-blocks: accuracy figures are refused (truth is never inferred)")
        else:
            self._draw_evaluation()
        self._draw_switch_and_phase()

    def _draw_evaluation(self):
        c, i, T = self.ctx, self.info, self.tables
        pts = T["block_evaluation"]
        rows = [(p, p) for p, _ in c.panels]
        und = int((pts.evaluated & ~pts.truth_defined).sum())
        ne = int((~pts.evaluated).sum())
        denom = (f"Denominators: evaluated truth-defined blocks of called runs ({und} truth-undefined and {ne} unevaluated "
                 "blocks are excluded and counted; no-call runs contribute no blocks).")
        ld, lcap = self._lengths()
        sets = ["all called blocks"] + (["EXACT_GENOTYPE_USABLE"] if c.has_headline else [])
        for metric, label in A.METRICS:
            if metric == "normalized_excess" and not c.has_normalization:
                self.skip("eval_block_normalized_excess_*", "the evaluation lacks truth1_bp, truth2_bp, called1_bp, "
                          "called2_bp: the diploid length denominator is unavailable, so length-normalised excess is refused")
                continue
            for s in sets:
                sub = pts if s == "all called blocks" else pts[pts.exact_genotype_usable == True]  # noqa: E712
                slug = "all" if s == "all called blocks" else "exact_genotype_usable"
                rule = "" if s == "all called blocks" else f" EXACT_GENOTYPE_USABLE is exactly: {EXACT_GENOTYPE_USABLE_RULE}."
                formula = (f" Diploid-length-normalised excess = {NORMALIZED_EXCESS_FORMULA}, a deterministic length scale, not an "
                           "observed alignment length." if metric == "normalized_excess" else "")
                main = " For leave-one-out rows excess is the primary accuracy metric; raw distance is shown separately." if metric != "distance" else \
                    " Raw distance includes what the panel cannot spell; for leave-one-out rows read the excess figures."
                self.figure(f"eval_block_{metric}_{slug}", "Ordered block performance", ["block_evaluation", "block_metric_summary"] + ld,
                            lambda t, sm, *l, m=metric, lab=label, sn=s: P.fig_ordered_blocks(
                                t if sn == "all called blocks" else t[t.exact_genotype_usable == True], sm, i, m, lab, sn,  # noqa: E712
                                rows, "panel", lengths=l[0] if l else None),
                            f"{label} per reported block, every block in model order, one row per experiment group "
                            f"({'; '.join(f'{p}: {", ".join(m)}' for p, m in c.panels)}); coverages dodged, medians black."
                            f"{rule}{formula}{main} {count_caption(metric, 'block_metric_summary')} {denom}{lcap}")
        if c.has_headline:
            pass
        else:
            self.skip("eval_block_*_exact_genotype_usable", "no headline filters: EXACT_GENOTYPE_USABLE cannot be formed")
        self.figure("eval_block_accuracy_heatmap", "Block accuracy", ["block_accuracy"] + ld,
                    lambda t, *l: P.fig_accuracy_heatmap(t, i, sets, l[0] if l else None),
                    "Floor attainment (excess 0; for leave-zero-out and mosaic cells the floor is 0 wherever the truth is "
                    "spellable, so this is exactness) by block and coverage, all called blocks versus EXACT_GENOTYPE_USABLE "
                    f"({EXACT_GENOTYPE_USABLE_RULE}). Cells: attained / denominator. {denom}{lcap}")
        self.figure("eval_block_accuracy_intervals", "Block accuracy", ["block_accuracy"] + ld,
                    lambda t, *l: P.fig_accuracy_intervals(t, i, sets, l[0] if l else None),
                    "The same rates as points with 95% Wilson intervals; numerators and denominators are in "
                    "data/block_accuracy.tsv." + lcap)
        if c.has_headline:
            self.figure("eval_filter_counts", "Headline-filter performance", ["filter_performance"],
                        lambda t: P.fig_filter_counts(t, i),
                        "Retained and filtered units, split by correctness, for each headline filter and usage policy. "
                        "Sequence policies: correct = excess 0 over evaluated truth-defined blocks. DOSAGE_CONTENT_USABLE "
                        f"({DOSAGE_CONTENT_USABLE_RULE}) is judged on the supplied dosage-content truth over dosage-sensitive "
                        "blocks. PHASE_FILTER is a separate use dimension, judged only on phase truth (data/filter_performance.tsv, "
                        "run level) and never folded into the exact-genotype policies.")
            self.figure("eval_filter_rates", "Headline-filter performance", ["filter_performance"],
                        lambda t: P.fig_filter_rates(t, i),
                        "Retention (retained / units), precision among retained, and bad-call detection sensitivity "
                        "(filtered incorrect / incorrect) by experiment and coverage.")
            self.figure("eval_filter_catastrophic", "Headline-filter performance", ["filter_performance"],
                        lambda t: P.fig_filter_catastrophic(t, i, c.catastrophic_bp),
                        f"Catastrophic errors (excess >= {c.catastrophic_bp:g} bp) retained and filtered per sequence policy.")
            self.figure("eval_filter_by_block", "Headline-filter performance", ["filter_performance"] + ld,
                        lambda t, *l: P.fig_filter_by_block(t, i, l[0] if l else None),
                        f"EXACT_GENOTYPE_USABLE ({EXACT_GENOTYPE_USABLE_RULE}) retention and precision by block and experiment."
                        + lcap)
        else:
            self.skip("eval_filter_*", "no headline filters in the outputs")
        self.figure("eval_gq_diagnostic", "GQ diagnostic", ["gq_diagnostic"],
                    lambda t: P.fig_gq_diagnostic(t, i, A.GQ_LABELS),
                    "DIAGNOSTIC ONLY. Floor attainment by conditional-GQ bin with counts and 95% Wilson intervals. GQ is "
                    "conditional on the scored candidate set and is not yet globally calibrated; no GQ threshold is "
                    "recommended.")
        if c.content is not None and c.has_headline:
            self.figure("dosage_content_usable", "Dosage", ["dosage_content_usable"], lambda t: P.fig_dosage_content(t, i),
                        f"DOSAGE_CONTENT_USABLE ({DOSAGE_CONTENT_USABLE_RULE}) over dosage-sensitive blocks (DOSAGE_FILTER not "
                        "NOT_APPLICABLE or NO_CALL): retention, and correctness of the reported dosage content against the "
                        "supplied content truth. Content is kept separate from exact sequence.")
            self.figure("dosage_content_vs_exact", "Dosage", ["dosage_content_vs_exact"], lambda t: P.fig_content_vs_exact(t, i),
                        "Correct dosage content versus exact-diplotype resolution (excess 0), for blocks where both truths "
                        "are defined: a correct content does not imply an exact diplotype.")
        else:
            self.skip("dosage_content_usable, dosage_content_vs_exact",
                      "no --dosage-content-evaluation (content truth is never inferred)" if c.content is None else "no headline filters")

    def _draw_switch_and_phase(self):
        c, i, T = self.ctx, self.info, self.tables
        self.figure("dosage_outcomes", "Dosage", ["dosage_states"], lambda t: P.fig_dosage_states(t, i, A.DOSAGE_OUTCOMES),
                    "Complete repair, content-only correction, SUSPECT, UNDERPOWERED and triggered-unrepaired counts per "
                    "experiment and coverage (the same table as qc_dosage_states). Truth-free.")
        if c.representatives is not None:
            self.figure("dosage_representatives", "Dosage", ["dosage_representatives"], lambda t: P.fig_representatives(t, i),
                        "NON-AUTHORITATIVE diagnostic: the sequence distance of a conditional representative (a serialisation "
                        "of an unresolved set) against the reported genotype's. Never a call and never an accuracy claim.")
        else:
            self.skip("dosage_representatives", "no --representative-diagnostics")
        if c.switch_plans is None:
            self.skip("mosaic_*", "no --switch-metadata (switches are never inferred from calls)")
        else:
            classes = A.switch_class_order(c)
            note = ("A plan alternates each homologue between its two selected source templates (A and B) at the listed "
                    "boundaries; it does not draw a new donor at each switch.")
            ld, lcap = self._lengths()
            self.figure("mosaic_switch_map", "Mosaic switches", ["switch_map"] + ld,
                        lambda t, *l: P.fig_switch_map(t, i, classes, l[0] if l else None),
                        "Predefined truth switch boundaries per plan, from the supplied switch metadata. The symbols are "
                        "truth-boundary annotations, not error observations. " + note + lcap)
            self.figure("mosaic_callability", "Mosaic switches", ["mosaic_callability"],
                        lambda t: P.fig_mosaic_callability(t, i, classes), "Call class by switch class and coverage. " + note)
            if c.evaluation is not None:
                pts = T["block_evaluation"]
                self.figure("mosaic_excess_by_class", "Mosaic switches", ["block_evaluation"],
                            lambda t: P.fig_mosaic_excess(t[t.switch_class.fillna("") != ""], i, classes),
                            "Excess and diploid-length-normalised excess per block by switch class and coverage (summaries in "
                            "data/mosaic_excess_summary.tsv). " + note)
                plans = T["switch_map"]
                markers = {}
                for k in classes:
                    markers[k] = [(p, dict(zip(g.boundary_index, g.category))) for p, g in plans[plans.switch_class == k].groupby("plan", sort=True)]
                for metric, label in A.METRICS:
                    if metric == "normalized_excess" and not c.has_normalization:
                        continue
                    self.figure(f"mosaic_block_{metric}", "Mosaic switches", ["block_evaluation", "mosaic_block_metric_summary"] + ld,
                                lambda t, sm, *l, m=metric, lab=label: P.fig_ordered_blocks(
                                    t[t.switch_class.fillna("") != ""], sm, i, m, lab, "all called blocks",
                                    [(k, k) for k in classes], "switch_class", markers, l[0] if l else None),
                                f"{label} by block, one row per switch class. The switch symbols inside the top of each "
                                "panel are truth-boundary annotations from the switch metadata (triangle up: haplotype 1, "
                                "down: haplotype 2, diamond: both), one per plan; they are not error observations. "
                                + count_caption(metric, "mosaic_block_metric_summary") + " " + note + lcap)
        # ---- phase
        sps = T["strict_phase_support"]
        self.figure("phase_callability", "Phase", ["phase_support"], lambda t: P.fig_phase_callability(t, i),
                    "Runs by call class and phase_support. Truth-free. " + self._phase_support_sentence(sps))
        if c.phase_boundaries is None:
            self.skip("phase_boundary_status, phase_accuracy, phase_by_*", "no --phase-boundaries")
        else:
            acc = T["phase_accuracy"]
            n, err = int(acc.assessable.sum()), int(acc.switch_errors.sum())
            zero = ("Zero observed switch errors occurred among assessable primary boundaries " if err == 0 else
                    f"{err} switch error(s) occurred among ") + f"({err}/{n} assessable boundaries)."
            ge = ("Genotype-equivalent calls are unphased and are excluded as not assessable, never counted correct.")
            support = self._phase_support_sentence(sps)
            self.figure("phase_boundary_status", "Phase", ["phase_boundary_status"], lambda t: P.fig_phase_boundary_status(t, i),
                        f"Adjacent-block boundaries by assessment status and reason. {ge} {support}")
            self.figure("phase_accuracy", "Phase", ["phase_accuracy"], lambda t: P.fig_phase_accuracy(t, i),
                        f"Correct boundaries and switch errors among assessable primary boundaries, with the switch-error "
                        f"rate (95% Wilson). {zero} {support} {ge}")
            self.figure("phase_by_boundary_type", "Phase", ["phase_by_boundary_type"], lambda t: P.fig_phase_by_boundary_type(t, i),
                        "Boundaries by truth switch type (from the switch metadata, not from predictions). " + zero)
            if c.switch_plans is not None:
                classes = A.switch_class_order(c)
                self.figure("phase_by_switch_class", "Phase", ["phase_by_switch_class"],
                            lambda t: P.fig_phase_by_class(t, i, classes),
                            f"Switch-error rate by mosaic switch class and coverage. {zero} {ge}")
            else:
                self.skip("phase_by_switch_class", "no --switch-metadata")
        if c.phase_calls is not None:
            self.figure("phase_locus_status", "Phase", ["phase_locus_status"], lambda t: P.fig_phase_locus(t, i),
                        "Full-locus phased-call status from the phase evaluation's per-call table: strict calls that are an "
                        "exact phased locus, other strict calls, genotype-equivalent calls (unphased, excluded rather than "
                        "counted correct) and no-calls.")
        else:
            self.skip("phase_locus_status", "no --phase-calls")
        if c.phase_projected is not None:
            self.figure("phase_projected_secondary", "Phase", ["phase_projected"], lambda t: P.fig_phase_projected(t, i),
                        "SECONDARY: the sequence-projected phase concordance is a concordance diagnostic, not primary "
                        "phase accuracy.")
        else:
            self.skip("phase_projected_secondary", "no --phase-projected")

    def _phase_support_sentence(self, sps):
        parts = []
        for e in self.ctx.experiments:
            s = sps[sps.experiment == e]
            if not len(s):
                continue
            total = int(s.strict_calls.sum())
            if s.phase_support.nunique() == 1:
                parts.append(f"{e}: all {total} strict calls {s.phase_support.iloc[0]}")
            else:
                parts.append(f"{e}: " + ", ".join(f"{int(r.strict_calls)} {r.phase_support}" for r in s.itertuples()))
        return ("Strict-call phase support - " + "; ".join(parts) + ".") if parts else ""

    # ---------------------------------------------------------------------------------------------- report and manifest
    def write_report(self):
        c, T, a = self.ctx, self.tables, self.args
        runs = c.runs
        called = runs.call_class.isin(CALLED)
        lines = [f"# {a.title}", "",
                 f"Generated by `{TOOL}` from a {self.input_kind} input. Inputs, their SHA-256 and every output's hash are in "
                 "[manifest.tsv](manifest.tsv); every plotted value is in a TSV under [data/](data/).", "",
                 "## Inputs and denominators", "",
                 f"- **Runs:** {len(runs)} ({int((runs.call_class == 'strict_call').sum())} strict, "
                 f"{int((runs.call_class == 'genotype_equivalent_call').sum())} genotype-equivalent, "
                 f"{int((runs.call_class == 'no_call').sum())} no-call, {int((~runs.call_class.isin(CALLED + ('no_call',))).sum())} "
                 "without a recorded call class).",
                 f"- **Reported blocks:** {len(c.blocks)} from {int(called.sum())} called runs; block model "
                 f"({len(c.block_order)} blocks, model order): {', '.join(c.block_order)}.",
                 f"- **Experiments:** {', '.join(c.experiments)}; coverages: {', '.join(c.coverages)}"
                 + ("" if a.run_metadata else " (no --run-metadata: every run is in one group)") + ".",
                 f"- **Modality:** " + (", ".join(c.modalities) + (" (figures drawn once per modality)" if len(c.modalities) > 1 else "")
                                        if c.modalities else "not supplied (none assumed; aggregates carry no modality column)") + ".",
                 f"- **Headline filters:** {'present' if c.has_headline else 'absent (outputs predate headline_filters_v2)'}.",
                 ]
        if c.evaluation is not None:
            pts = T["block_evaluation"]
            lines += [f"- **Evaluation:** {int(pts.evaluated.sum())} evaluated blocks; {int((pts.evaluated & ~pts.truth_defined).sum())} "
                      f"truth-undefined (excluded from every accuracy denominator and counted); {int((~pts.evaluated).sum())} "
                      "reported blocks without an evaluation row (excluded and counted)."]
        else:
            lines += ["- **Evaluation:** none supplied. This is a truth-free QC report: nothing in it is labelled accurate or "
                      "erroneous."]
        lines += ["", "## Definitions", "",
                  "- **Call class:** strict_call (a phased diplotype), genotype_equivalent_call (an unordered genotype, "
                  "unphased), no_call (no genotype; the reason is recorded).",
                  "- **Headline filters:** GENOTYPE_FILTER, DOSAGE_FILTER and PHASE_FILTER, each PASS, LOW_QUALITY, "
                  "NOT_EVALUATED, NO_CALL or NOT_APPLICABLE (genotype_headline_filters.md). A locus value is the worst of "
                  "its units.",
                  f"- **EXACT_GENOTYPE_USABLE:** exactly `{EXACT_GENOTYPE_USABLE_RULE}`. It is never called \"high quality\".",
                  f"- **DOSAGE_CONTENT_USABLE:** `{DOSAGE_CONTENT_USABLE_RULE}`, judged on dosage content, not sequence.",
                  f"- **PHASED_GENOTYPE_USABLE:** `{PHASED_GENOTYPE_USABLE_RULE}`; PHASE_FILTER is an independent use "
                  "dimension and is never folded into exact-genotype filtering.",
                  "- **Evaluator metrics (the evaluator's own fields, never recomputed):** raw diploid distance = "
                  "`unordered_distance`; availability floor = `availability_floor`; excess = `excess` (the leave-one-out "
                  "metric). Floor attained = excess 0. A block is truth-defined when the evaluator records both truth "
                  "allele lengths.",
                  f"- **Diploid-length-normalised excess:** `{NORMALIZED_EXCESS_FORMULA}` (field `normalized_excess`; its "
                  "denominator is `diploid_length_denominator_bp`), from the evaluator's recorded allele lengths. The "
                  "denominator is a deterministic diploid length scale, not an observed alignment length. Refused when "
                  "those lengths are absent.",
                  f"- **Catastrophic:** excess >= {c.catastrophic_bp:g} bp.",
                  "- **Intervals:** 95% Wilson score intervals.",
                  "- **GQ:** conditional on the scored candidate set and not globally calibrated; no threshold is recommended.",
                  "- **Axes:** distance and excess use a symlog axis (linear in [0, 1] bp, logarithmic above) so zero-valued "
                  "calls and catastrophic outliers both stay visible; length-normalised excess uses linear [0, 1e-4], logarithmic above.",
                  "- **Runtime pooling:** runs are pooled only within a compatibility group (settings digest, threads, "
                  "execution group" + (", modality" if c.modalities else "") + "); the genotyper records no host or "
                  "binary identity, so `execution_group` must distinguish machine, environment and comparable workload "
                  "conditions whenever runtime is compared.", ""]
        if c.evaluation is not None:
            lines.insert(len(lines) - 1,   # the definitions end with a blank line
                         "- **Block counts (ordered block figures):** n and zero are coverage-specific: one label per "
                         "coverage above each block, from `block_metric_summary.tsv` (n: called, evaluated, truth-defined, "
                         "retained by the plotted set, non-missing; zero: among n, value exactly 0 - availability floor "
                         "attained for excess and normalised excess, exact sequence distance zero for raw distance).")
        if "block_lengths" in T:
            lines.insert(len(lines) - 1, "- **Block length labels:** " + BLOCK_LENGTH_CAPTION)
        lines += ["## Figures", ""]
        section = None
        for sec, name, caption, data in self.figures:
            if sec != section:
                lines += [f"### {sec}", ""]
                section = sec
            lines += [f"- **[{name}](figures/{name}.png)** ([pdf](figures/{name}.pdf)): {caption} "
                      f"Data: {', '.join(f'[{d}.tsv](data/{d}.tsv)' for d in data)}."]
        lines += ["", "## Skipped figures", ""]
        lines += [f"- `{f}`: {why}." for f, why in self.skipped] or ["- none"]
        lines += ["", "## Data tables", ""]
        lines += [f"- [{n}.tsv](data/{n}.tsv): {len(t)} rows" for n, t in self.tables.items()]
        with open(os.path.join(self.out, "report.md"), "w") as f:
            f.write("\n".join(lines) + "\n")

    def write_manifest(self):
        rows = [("meta", "command", " ".join(self.argv), "", ""), ("meta", "tool", TOOL, "", ""),
                ("meta", "python", platform.python_version(), "", ""), ("meta", "pandas", pd.__version__, "", ""),
                ("meta", "numpy", np.__version__, "", ""), ("meta", "matplotlib", matplotlib.__version__, "", "")]
        repo = os.path.dirname(os.path.dirname(SCRIPT_DIR))
        try:
            commit = subprocess.run(["git", "-C", repo, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
            dirty = subprocess.run(["git", "-C", repo, "status", "--porcelain", "--", "scripts/genotype_report",
                                    "scripts/plot_genotype_results.py"], capture_output=True, text=True).stdout.strip()
            rows.append(("meta", "script_commit", (commit or "unknown") + (" (script files modified)" if dirty else ""), "", ""))
        except OSError:
            rows.append(("meta", "script_commit", "unknown", "", ""))
        scripts = sorted(os.path.join(SCRIPT_DIR, f) for f in os.listdir(SCRIPT_DIR) if f.endswith(".py"))
        scripts.append(os.path.join(repo, "scripts", "plot_genotype_results.py"))
        for p in scripts:
            if os.path.isfile(p):
                rows.append(("script", os.path.relpath(p, repo), p, sha256_file(p), str(os.path.getsize(p))))
        for p in sorted(set(os.path.abspath(s) for s in self.sources)):
            rows.append(("input", os.path.basename(p), p, sha256_file(p), str(os.path.getsize(p))))
        outs = []
        for d, _, files in os.walk(self.out):
            for f in files:
                p = os.path.join(d, f)
                if os.path.relpath(p, self.out) != "manifest.tsv":
                    outs.append(p)
        for p in sorted(outs):
            rows.append(("output", os.path.relpath(p, self.out), os.path.relpath(p, self.out), sha256_file(p), str(os.path.getsize(p))))
        write_tsv(pd.DataFrame(rows, columns=["kind", "name", "value", "sha256", "bytes"]), os.path.join(self.out, "manifest.tsv"))


def run(args, argv):
    r = Report(args, argv)
    r.load()
    os.makedirs(r.out, exist_ok=True)
    r.aggregate()
    r.draw()
    r.write_report()
    r.write_manifest()
    return r
