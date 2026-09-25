#!/usr/bin/env python3
"""Consolidate the final C4 production campaigns into one package, and draw the figure.

    scripts/consolidate_c4_production.py --out <campaign dir>

Reads the campaign's jobs.tsv, completeness.tsv and every job's status, evaluation and timing, and writes
under <campaign dir>/package/:

    runs.tsv              one row per intended job: outcome, call class, no-call reason, phase support,
                          distance, floor, excess, exact blocks, runtime and memory
    blocks.tsv            one row per job and block of every CALL (strict or genotype-equivalent), every
                          quality -- the same rows as blocks_all_quality.tsv
    blocks_all_quality.tsv   every called block, with the caller's per-block filter (PASS / LOWQUAL_* /
                          NOT_EVALUATED / NA)
    blocks_pass_only.tsv  representation-PASS: the called blocks whose filter is exactly PASS (fail closed:
                          NOT_EVALUATED and NA are not PASS)
    blocks_high_quality.tsv   the caller's combined block quality: quality_filter == PASS (block_quality_rule_v1 --
                          representation PASS and genotype_support READS), read, never rebuilt here; NA (a run
                          that predates the column) is not PASS. Not calibrated accuracy: a PASS block can still
                          be wrong where both channels agree with the call (HG00512.2 sv:1). This is the
                          BLOCK-level field: a call's locus_filter (the representation guard's locus summary) can
                          be PASS while most of its blocks are LOW_QUALITY (no_read_support), and it never
                          decides this table
    blocks_undefined_truth.tsv   the called blocks whose truth is undefined (a truth haplotype does not spell the
                          block). They stay in the three tables above, flagged truth_defined false, and are
                          counted apart (undefined_truth_blocks): no block-accuracy denominator includes them
    summary.tsv           per experiment and coverage: strict / genotype-equivalent / no-call rates, no-call
                          reasons, exact genotype rate over all intended samples and given a call, exact blocks
                          and SV blocks, total and normalized edit distance, distance above the availability
                          floor, phase-support counts, runtime and memory distributions
    fallback_diagnostic.tsv   the initial-union genotype of every improving-move-tie no-call, evaluated as a
                          LABELLED DIAGNOSTIC FALLBACK. It is never a call and never enters any rate above.
    c4_production_block_distance.png / .pdf   the figure

Called-only accuracy uses calls only: a no-call contributes no block row and enters no block denominator. Over-all-
intended rates count every no-call as not exact. Block-accuracy denominators count truth-defined blocks only.

A call is decided by the final status and call class (and search_conditional names a tie-branched one), never
by repair_converged: a tie-branched COMPLETE call keeps the production search's repair_converged false. The
quality filter never turns a call into a no-call: callability counts every call whatever its filter, and the
LowQual rate is reported beside it, separately, over calls only, and so is the NOT_EVALUATED rate; each of
the three block tables reports its denominator and the fraction of called blocks it retains.
"""
import argparse
import csv
import os
import pathlib

import numpy as np
import pandas as pd

EXPERIMENTS = [("lzo", "Ordinary LZO"), ("mosaic", "Synthetic-mosaic LZO"), ("loo", "LOO")]
COVERAGES = [2, 5, 10, 30]


def key_values(path):
    out = {}
    try:
        with open(path) as handle:
            for line in handle:
                parts = line.rstrip("\n").split("\t", 1)
                if len(parts) == 2 and parts[0] not in out:
                    out[parts[0]] = parts[1]
    except FileNotFoundError:
        return {}
    return out


def number(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return np.nan


def distribution(series, prefix):
    values = pd.to_numeric(series, errors="coerce").dropna()
    if values.empty:
        return {f"{prefix}_median": "NA", f"{prefix}_p90": "NA", f"{prefix}_max": "NA"}
    return {f"{prefix}_median": f"{values.median():.6g}", f"{prefix}_p90": f"{values.quantile(0.9):.6g}",
            f"{prefix}_max": f"{values.max():.6g}"}


def consolidate(out, figure=True):
    package = out / "package"
    package.mkdir(exist_ok=True)
    jobs = pd.read_csv(out / "jobs.tsv", sep="\t", dtype=str)
    ledger = pd.read_csv(out / "completeness.tsv", sep="\t", dtype=str).set_index("job")
    runs, blocks, fallback = [], [], []
    for job in jobs.itertuples(index=False):
        prefix = job.out_prefix
        status = key_values(prefix + ".mosaic_status.tsv")
        summary = key_values(prefix + ".eval.summary.tsv")
        entry = ledger.loc[job.job] if job.job in ledger.index else {}
        outcome = entry.get("outcome", "pending") if len(entry) else "pending"
        called = outcome in ("strict_call", "genotype_equivalent_call")
        truth_bp = np.nan
        block_rows = []
        if os.path.exists(prefix + ".eval.blocks.tsv"):
            frame = pd.read_csv(prefix + ".eval.blocks.tsv", sep="\t", dtype=str)
            truth_bp = pd.to_numeric(frame.truth1_bp, errors="coerce").sum() + \
                pd.to_numeric(frame.truth2_bp, errors="coerce").sum()
            if called:
                for row in frame.itertuples(index=False):
                    d = number(row.unordered_distance)
                    t = number(row.truth1_bp) + number(row.truth2_bp)
                    defined = not (np.isnan(number(row.truth1_bp)) or np.isnan(number(row.truth2_bp)))
                    block_rows.append({"job": job.job, "experiment": job.campaign, "subject": job.subject,
                                       "coverage": int(job.coverage), "call_class": outcome, "block_id": row.block_id,
                                       "kind": row.kind, "truth_defined": "true" if defined else "false",
                                       "unordered_distance": row.unordered_distance,
                                       "normalized_distance": d / t if t and not np.isnan(d) else np.nan,
                                       "availability_floor": row.availability_floor, "excess": row.excess,
                                       "genotype_support": row.genotype_support,
                                       "phase_support": status.get("phase_support", "NA"),
                                       "filter": getattr(row, "filter", "NA"),
                                       "quality_filter": getattr(row, "quality_filter", "NA"),
                                       "quality_reason": getattr(row, "quality_reason", "NA")})
        blocks.extend(block_rows)
        distance = number(summary.get("unordered_distance_total")) if called else np.nan
        reason = status.get("no_call_reason", "NA")
        if outcome == "no_call" and reason in ("NA", ""):
            # Stops before the repair write no no_call_reason; name them from the status itself.
            state = status.get("status", "")
            note = status.get("note", "")
            reason = ("union_scoring_cap" if state == "RESCORE_INCOMPLETE" and "resource cap" in note
                      else state.lower() if state else "unknown")
        runs.append({
            "job": job.job, "experiment": job.campaign, "subject": job.subject, "coverage": int(job.coverage),
            "outcome": outcome, "status": status.get("status", "NA"),
            "call_class": status.get("call_class", "NA"), "no_call_reason": reason if outcome == "no_call" else "none",
            "phase_support": status.get("phase_support", "NA") if called else "NA",
            "search_conditional": status.get("search_conditional", "none"),
            "locus_filter": status.get("locus_filter", "NA") if called else "NA",
            "lowqual_blocks": status.get("lowqual_blocks", "NA") if called else "NA",
            "not_evaluated_blocks": status.get("not_evaluated_blocks", "NA") if called else "NA",
            "representation_guard_complete": status.get("representation_guard_complete", "NA") if called else "NA",
            "pass_distance": number(summary.get("pass_unordered_distance_total")) if called else np.nan,
            "pass_excess": number(summary.get("pass_excess_over_floor_total")) if called else np.nan,
            "high_quality_excess": number(summary.get("high_quality_excess_over_floor_total")) if called else np.nan,
            "final_decision": status.get("final_decision", "NA"), "top_set_size": status.get("top_set_size", "NA"),
            "distance": distance, "floor": number(summary.get("availability_floor_total")),
            "excess": number(summary.get("excess_over_floor_total")) if called else np.nan,
            "normalized_distance": distance / truth_bp if called and truth_bp else np.nan,
            "blocks": number(summary.get("blocks")), "sv_blocks": number(summary.get("sv_blocks")),
            "undefined_truth_blocks": number(summary.get("undefined_truth_blocks")) if called else np.nan,
            "blocks_truth_defined": sum(r["truth_defined"] == "true" for r in block_rows) if called else np.nan,
            "sv_blocks_truth_defined": sum(r["truth_defined"] == "true" and r["kind"] == "SV" for r in block_rows)
                                       if called else np.nan,
            "quality_lowqual_blocks": status.get("quality_lowqual_blocks", "NA") if called else "NA",
            "quality_lowqual_block_ids": status.get("quality_lowqual_block_ids", "NA") if called else "NA",
            "blocks_exact": number(summary.get("blocks_exact")) if called else np.nan,
            "sv_blocks_exact": number(summary.get("sv_blocks_exact")) if called else np.nan,
            "exact": called and distance == 0,
            "repair_rounds": status.get("repair_rounds_run", "NA"), "repair_moves": status.get("repair_moves", "NA"),
            "attempts": entry.get("attempts", "NA") if len(entry) else "NA",
            "wall_seconds": number(entry.get("wall_seconds")) if len(entry) else np.nan,
            "cpu_seconds": number(entry.get("user_seconds")) + number(entry.get("sys_seconds")) if len(entry) else np.nan,
            "peak_footprint_bytes": number(entry.get("peak_footprint_bytes")) if len(entry) else np.nan,
            "peak_rss_bytes": number(entry.get("peak_rss_bytes")) if len(entry) else np.nan,
        })
        if outcome == "no_call" and status.get("no_call_reason") == "improving_move_tie":
            diag = key_values(prefix + ".fallback.summary.tsv")
            fallback.append({"job": job.job, "experiment": job.campaign, "subject": job.subject,
                             "coverage": int(job.coverage), "label": "DIAGNOSTIC FALLBACK - NOT A CALL",
                             "evaluated_diplotype": diag.get("evaluated_diplotype", "NA"),
                             "union_ranking_status": status.get("union_ranking_status", "NA"),
                             "union_top_hap1": status.get("union_top_hap1", "NA"),
                             "union_top_hap2": status.get("union_top_hap2", "NA"),
                             "distance": diag.get("unordered_distance_total", "NA"),
                             "floor": diag.get("availability_floor_total", "NA"),
                             "excess": diag.get("excess_over_floor_total", "NA"),
                             "blocks_exact": diag.get("blocks_exact", "NA"),
                             "sv_blocks_exact": diag.get("sv_blocks_exact", "NA")})
    runs = pd.DataFrame(runs)
    blocks = pd.DataFrame(blocks)
    runs.to_csv(package / "runs.tsv", sep="\t", index=False)
    blocks.to_csv(package / "blocks.tsv", sep="\t", index=False)
    blocks.to_csv(package / "blocks_all_quality.tsv", sep="\t", index=False)
    passing = (blocks["filter"].astype(str) == "PASS") if not blocks.empty else pd.Series(dtype=bool)
    quality = (blocks["quality_filter"].astype(str) == "PASS") if not blocks.empty else pd.Series(dtype=bool)
    (blocks[passing] if not blocks.empty else blocks).to_csv(package / "blocks_pass_only.tsv", sep="\t", index=False)
    (blocks[quality] if not blocks.empty else blocks).to_csv(package / "blocks_high_quality.tsv", sep="\t",
                                                             index=False)
    (blocks[blocks["truth_defined"] == "false"] if not blocks.empty else blocks).to_csv(
        package / "blocks_undefined_truth.tsv", sep="\t", index=False)
    pd.DataFrame(fallback).to_csv(package / "fallback_diagnostic.tsv", sep="\t", index=False)

    rows = []
    reasons = sorted(set(runs.no_call_reason) - {"none"})
    for experiment, _ in EXPERIMENTS:
        for coverage in COVERAGES:
            part = runs[(runs.experiment == experiment) & (runs.coverage == coverage)]
            if part.empty:
                continue
            calls = part[part.outcome.isin(["strict_call", "genotype_equivalent_call"])]
            n = len(part)
            row = {"experiment": experiment, "coverage": coverage, "intended": n,
                   "accounted": int(part.outcome.isin(["strict_call", "genotype_equivalent_call", "no_call"]).sum()),
                   "strict_calls": int((part.outcome == "strict_call").sum()),
                   "genotype_equivalent_calls": int((part.outcome == "genotype_equivalent_call").sum()),
                   "no_calls": int((part.outcome == "no_call").sum()),
                   "failed_or_pending": int((~part.outcome.isin(["strict_call", "genotype_equivalent_call",
                                                                   "no_call"])).sum())}
            row["strict_call_rate"] = row["strict_calls"] / n
            row["genotype_equivalent_call_rate"] = row["genotype_equivalent_calls"] / n
            row["no_call_rate"] = row["no_calls"] / n
            for reason in reasons:
                row[f"no_call_{reason}"] = int((part.no_call_reason == reason).sum())
            row["exact_genotype_rate_all_intended"] = int(part.exact.sum()) / n
            row["exact_genotype_rate_given_call"] = int(calls.exact.sum()) / len(calls) if len(calls) else np.nan
            row["floor_attained_rate_given_call"] = float((calls.excess == 0).mean()) if len(calls) else np.nan
            # Exact blocks over the truth-defined blocks of calls (an undefined block can be neither exact nor wrong).
            row["exact_blocks_given_call"] = (calls.blocks_exact.sum() / calls.blocks_truth_defined.sum()
                                              if len(calls) and calls.blocks_truth_defined.sum() else np.nan)
            row["exact_sv_blocks_given_call"] = (calls.sv_blocks_exact.sum() / calls.sv_blocks_truth_defined.sum()
                                                 if len(calls) and calls.sv_blocks_truth_defined.sum() else np.nan)
            row["distance_total_given_call"] = calls.distance.sum() if len(calls) else np.nan
            row["distance_mean_given_call"] = calls.distance.mean() if len(calls) else np.nan
            row["normalized_distance_mean_given_call"] = calls.normalized_distance.mean() if len(calls) else np.nan
            row["excess_total_given_call"] = calls.excess.sum() if len(calls) else np.nan
            row["excess_mean_given_call"] = calls.excess.mean() if len(calls) else np.nan
            # Quality, reported beside callability and never folded into it. lowqual_calls is the REPRESENTATION
            # guard's locus filter; the combined block quality is quality_filter, counted separately below.
            flagged = calls.locus_filter.astype(str).str.startswith("LOWQUAL")
            row["tie_branched_calls"] = int((calls.search_conditional == "tie_branched").sum())
            row["lowqual_calls"] = int(flagged.sum())
            row["lowqual_call_rate_given_call"] = float(flagged.mean()) if len(calls) else np.nan
            here = blocks[(blocks.experiment == experiment) & (blocks.coverage == coverage)] if not blocks.empty else blocks
            flagged_blocks = here["filter"].astype(str).str.startswith("LOWQUAL") if not here.empty else pd.Series(dtype=bool)
            unevaluated = calls.locus_filter.astype(str) == "NOT_EVALUATED"
            row["not_evaluated_calls"] = int(unevaluated.sum())
            row["not_evaluated_call_rate_given_call"] = float(unevaluated.mean()) if len(calls) else np.nan
            row["called_blocks"] = len(here)
            row["lowqual_blocks"] = int(flagged_blocks.sum()) if len(here) else 0
            row["block_filter_rate_given_call"] = float(flagged_blocks.mean()) if len(here) else np.nan
            filters = here["filter"].astype(str) if len(here) else pd.Series(dtype=str)
            not_evaluated = ~filters.isin(["PASS"]) & ~filters.str.startswith("LOWQUAL") if len(here) else pd.Series(dtype=bool)
            row["not_evaluated_blocks"] = int(not_evaluated.sum()) if len(here) else 0
            row["not_evaluated_block_rate_given_call"] = float(not_evaluated.mean()) if len(here) else np.nan
            # The three block tables, each with its accuracy denominator (truth-defined called blocks) and retained
            # fraction; undefined-truth blocks are counted apart.
            defined = (here["truth_defined"] == "true") if len(here) else pd.Series(dtype=bool)
            pass_here = (filters == "PASS") & defined if len(here) else pd.Series(dtype=bool)
            high_here = (here["quality_filter"].astype(str) == "PASS") & defined if len(here) else pd.Series(dtype=bool)
            row["undefined_truth_blocks"] = int((~defined).sum()) if len(here) else 0
            row["all_quality_blocks"] = int(defined.sum()) if len(here) else 0
            row["representation_pass_blocks"] = int(pass_here.sum()) if len(here) else 0
            row["representation_pass_fraction"] = float(pass_here.sum() / defined.sum()) if len(here) and defined.sum() else np.nan
            row["high_quality_blocks"] = int(high_here.sum()) if len(here) else 0
            row["high_quality_fraction"] = float(high_here.sum() / defined.sum()) if len(here) and defined.sum() else np.nan
            # The combined block quality (quality_filter) over every called block, and the calls carrying any
            # LOW_QUALITY block -- not the locus filter.
            qualities = here["quality_filter"].astype(str) if len(here) else pd.Series(dtype=str)
            row["quality_pass_blocks"] = int((qualities == "PASS").sum()) if len(here) else 0
            row["quality_low_quality_blocks"] = int((qualities == "LOW_QUALITY").sum()) if len(here) else 0
            row["quality_unrecorded_blocks"] = int((~qualities.isin(["PASS", "LOW_QUALITY"])).sum()) if len(here) else 0
            low_jobs = set(here.job[qualities == "LOW_QUALITY"]) if len(here) else set()
            row["calls_with_low_quality_blocks"] = int(calls.job.isin(low_jobs).sum())
            row["excess_total_pass_only_given_call"] = calls.pass_excess.sum() if len(calls) else np.nan
            row["excess_total_high_quality_given_call"] = calls.high_quality_excess.sum() if len(calls) else np.nan
            for support in ("READS_WITHIN_CANDIDATES", "PRIOR_ONLY", "UNRESOLVED"):
                row[f"phase_{support}"] = int((calls.phase_support == support).sum())
            row.update(distribution(part.wall_seconds, "wall_seconds"))
            row.update(distribution(part.peak_footprint_bytes, "peak_footprint_bytes"))
            rows.append(row)
    pd.DataFrame(rows).to_csv(package / "summary.tsv", sep="\t", index=False)
    if figure:
        plot(runs, blocks, package)
    print(f"package written to {package}")


def plot(runs, blocks, package):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colour = {"SV": "#B2182B", "context": "#2166AC"}
    marker = {"strict_call": "o", "genotype_equivalent_call": "^"}
    reason_colour = {}
    palette = ["#4D4D4D", "#E08214", "#8073AC", "#1B7837", "#D6604D", "#878787"]
    fig = plt.figure(figsize=(15, 7.2))
    grid = fig.add_gridspec(2, 3, height_ratios=[4, 1], hspace=0.08, wspace=0.18)
    rng = np.random.default_rng(20260923)
    for column, (experiment, title) in enumerate(EXPERIMENTS):
        axis = fig.add_subplot(grid[0, column])
        lower = fig.add_subplot(grid[1, column], sharex=axis)
        part = blocks[blocks.experiment == experiment] if not blocks.empty else blocks
        for i, coverage in enumerate(COVERAGES):
            for kind, offset in (("SV", -0.18), ("context", 0.18)):
                here = part[(part.coverage == coverage) & (part.kind == kind)] if not part.empty else part
                values = pd.to_numeric(here.unordered_distance, errors="coerce").dropna() if not here.empty else []
                if len(values) == 0:
                    continue
                shown = np.log10(np.asarray(values, dtype=float) + 1.0)
                if np.ptp(shown) > 0:
                    parts = axis.violinplot([shown], positions=[i + offset], widths=0.32, showextrema=False)
                    for body in parts["bodies"]:
                        body.set_facecolor(colour[kind])
                        body.set_alpha(0.28)
                        body.set_edgecolor("none")
                for call_class, symbol in marker.items():
                    mask = (here.call_class == call_class).to_numpy()
                    if not mask.any():
                        continue
                    y = np.log10(pd.to_numeric(here.unordered_distance[mask], errors="coerce").to_numpy() + 1.0)
                    x = i + offset + rng.uniform(-0.1, 0.1, size=len(y))
                    axis.scatter(x, y, s=7, marker=symbol, color=colour[kind], alpha=0.55, linewidths=0)
        axis.set_title(title)
        # The axis reaches the largest distance in ANY panel, so no block falls off the top unseen and the
        # three experiments share one scale.
        largest = pd.to_numeric(blocks.unordered_distance, errors="coerce").max() if not blocks.empty else 0
        top = max(20000.0, float(largest) if largest == largest else 0.0)
        ticks = [t for t in (0, 1, 10, 100, 1000, 10000, 100000, 1000000) if t <= top * 10]
        axis.set_yticks(np.log10(np.array(ticks) + 1.0))
        axis.set_yticklabels([f"{t:,}" for t in ticks])
        axis.set_ylim(-0.15, np.log10(top + 1.0) + 0.12)
        if column == 0:
            axis.set_ylabel("phase-invariant block distance from truth (bp, log)")
        plt.setp(axis.get_xticklabels(), visible=False)
        axis.grid(axis="y", color="#DDDDDD", linewidth=0.6)
        # No-calls are runs, not blocks, and get no distance: their counts by reason, drawn apart.
        sub = runs[runs.experiment == experiment]
        bottoms = np.zeros(len(COVERAGES))
        for reason in sorted(set(sub.no_call_reason) - {"none"}):
            if reason not in reason_colour:
                reason_colour[reason] = palette[len(reason_colour) % len(palette)]
            counts = np.array([((sub.coverage == c) & (sub.no_call_reason == reason)).sum() for c in COVERAGES])
            lower.bar(range(len(COVERAGES)), counts, bottom=bottoms, width=0.5, color=reason_colour[reason])
            bottoms += counts
        totals = [(sub.coverage == c).sum() for c in COVERAGES]
        for i, (made, total) in enumerate(zip(bottoms, totals)):
            lower.text(i, made + 0.2, f"{int(made)}/{int(total)}", ha="center", va="bottom", fontsize=8)
        lower.set_xticks(range(len(COVERAGES)))
        lower.set_xticklabels([f"{c}x" for c in COVERAGES])
        lower.set_ylim(0, max(max(totals) if totals else 1, 1) * 1.25)
        if column == 0:
            lower.set_ylabel("no-calls")
    handles = [plt.Line2D([], [], color=colour[k], marker="s", linestyle="", markersize=8, alpha=0.5, label=f"{k} block")
               for k in colour]
    handles += [plt.Line2D([], [], color="#555555", marker=m, linestyle="", markersize=6,
                           label="strict call" if c == "strict_call" else "genotype-equivalent call")
                for c, m in marker.items()]
    handles += [plt.Rectangle((0, 0), 1, 1, color=v, label=f"no-call: {k.replace('_', ' ')}")
                for k, v in reason_colour.items()]
    fig.legend(handles=handles, loc="lower center", ncol=min(len(handles), 6), frameon=False, fontsize=9)
    fig.suptitle("C4 production caller (panel-anchored adaptive completion): per-block distance from truth")
    fig.subplots_adjust(bottom=0.14, top=0.9)
    for suffix in ("png", "pdf"):
        fig.savefig(package / f"c4_production_block_distance.{suffix}", dpi=160)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", required=True)
    parser.add_argument("--no-figure", action="store_true", help="write the tables only")
    args = parser.parse_args()
    consolidate(pathlib.Path(args.out).resolve(), figure=not args.no_figure)


if __name__ == "__main__":
    main()
