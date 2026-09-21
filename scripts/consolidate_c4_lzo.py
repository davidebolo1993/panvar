#!/usr/bin/env python3
"""Consolidate the C4 leave-zero-out results into one reproducible package.

    scripts/consolidate_c4_lzo.py --heldout <run dir> --manifest <manifest.tsv>
                                  --ordinary <cohort.tsv> --out <dir>

Read-only. It runs no caller, simulates nothing and changes no call: it turns finished outputs into
tables and plots, and writes the raw TSV behind every figure.

Two experiments are reported side by side and NEVER pooled, because their truths are different objects:

  ordinary LZO          truth is two complete panel haplotypes (the held-out donor's own rows).
  synthetic mosaic LZO  truth is assembled block by block from panel haplotypes and is generally not a
                        panel row at all.

Within the mosaic experiment the two arms are also never pooled. K=0 conditions on every panel
diplotype; K=128 conditions on those plus up to 128 proposed mosaics. Their conditional quantities --
GQ above all -- answer different questions, so a pooled GQ or a pooled accuracy would be meaningless.

Three more separations are enforced throughout:
  * calls and no-calls (a refused run has no genotype, and NA is never written as zero);
  * strict calls and CONDITIONAL_CONSENSUS_DIAGNOSTIC rows (the latter are not calls and enter no
    strict rate);
  * SV and context blocks wherever a stratified number is reported.
"""

import argparse
import csv
import os
import statistics

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

COVERAGES = ["30", "10", "5", "2"]
ARMS = ["K0", "K128"]
ARM_K = {"K0": "0", "K128": "128"}
# One fixed colour per arm, used in every figure so the arms are never visually confused.
ARM_COLOUR = {"K0": "#2166AC", "K128": "#B2182B"}
NA = "NA"


def read_kv(path):
    out = {}
    if not os.path.exists(path):
        return out
    for line in open(path):
        parts = line.rstrip("\n").split("\t")
        if len(parts) == 2:
            out[parts[0]] = parts[1]
    return out


def read_tsv(path):
    if not os.path.exists(path):
        return []
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, header, rows):
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def number(value):
    """A value that is genuinely numeric, or None. NA and empty are None, never zero."""
    if value is None or value in ("", NA, "absent"):
        return None
    try:
        return float(value)
    except ValueError:
        return None


def cell(value):
    return NA if value is None else value


def plans_of(manifest):
    rows = []
    for line in open(manifest):
        parts = line.rstrip("\n").split("\t")
        if len(parts) >= 9 and parts[0] == "plan" and parts[1] != "name":
            rows.append((parts[1], parts[2]))
    return rows


# --------------------------------------------------------------------------------------------------
# Run-level table
# --------------------------------------------------------------------------------------------------

RUN_HEADER = [
    "plan", "class", "coverage", "arm", "status", "candidate_scope", "final_candidate_scope",
    "search_exhaustive", "tied_winners", "fragments_total", "wall_s", "max_rss_bytes",
    "phase_support", "phased_haplotype_emitted", "genotype_support", "min_block_gq",
    "separation_gap", "exact_blocks", "total_distance", "max_block_distance",
    "exact_sv_blocks", "exact_context_blocks",
]


def build_runs(heldout, plans):
    diagnostics = {}
    for row in read_tsv(os.path.join(heldout, "diagnostics.tsv")):
        diagnostics[(row["plan"], row["coverage"], row["K"])] = row
    runtime = {}
    for row in read_tsv(os.path.join(heldout, "runtime.tsv")):
        runtime[(row["plan"], row["coverage"])] = row

    rows = []
    for plan, klass in plans:
        for coverage in COVERAGES:
            for arm in ARMS:
                k = ARM_K[arm]
                prefix = os.path.join(heldout, plan, f"call.{coverage}.K{k}")
                status = read_kv(prefix + ".mosaic_status.tsv")
                summary = read_kv(os.path.join(heldout, plan, f"eval.{coverage}.K{k}.summary.tsv"))
                diag = diagnostics.get((plan, coverage, k), {})
                time_row = runtime.get((plan, coverage), {})
                complete = status.get("status") == "COMPLETE"

                # Per-block maxima come from the evaluator's own table, and only for a real call.
                worst, exact_context = None, None
                if complete:
                    blocks = read_tsv(os.path.join(heldout, plan, f"eval.{coverage}.K{k}.blocks.tsv"))
                    distances = [number(b["unordered_distance"]) for b in blocks]
                    distances = [d for d in distances if d is not None]
                    worst = int(max(distances)) if distances else None
                    exact_context = sum(
                        1 for b in blocks
                        if b["kind"] == "context" and number(b["unordered_distance"]) == 0
                    )

                rows.append([
                    plan, klass, coverage, arm,
                    status.get("status", NA),
                    status.get("candidate_scope", NA),
                    status.get("final_candidate_scope", NA),
                    status.get("proposal_search_exhaustive", NA),
                    status.get("winner_set_size", NA),
                    cell(diag.get("fragments_total")),
                    cell(time_row.get("wall_s")),
                    cell(time_row.get("max_rss_bytes")),
                    status.get("phase_support", NA),
                    "yes" if os.path.exists(prefix + ".mosaic_haplotypes.fa") else "no",
                    cell(diag.get("genotype_support")),
                    cell(diag.get("min_gq")),
                    cell(diag.get("separation_gap")),
                    cell(summary.get("blocks_exact")) if complete else NA,
                    cell(summary.get("unordered_distance_total")) if complete else NA,
                    cell(worst) if complete else NA,
                    cell(summary.get("sv_blocks_exact")) if complete else NA,
                    cell(exact_context) if complete else NA,
                ])
    return rows


# --------------------------------------------------------------------------------------------------
# Per-block table: strict calls and conditional consensus, labelled and never mixed
# --------------------------------------------------------------------------------------------------

BLOCK_HEADER = [
    "plan", "class", "coverage", "arm", "source", "block_id", "block_index", "kind",
    "distance", "exact", "truth1_bp", "truth2_bp", "normalized_distance",
    "genotype_support", "conditional_gq", "decided_by",
]


def build_blocks(heldout, plans):
    rows = []
    for plan, klass in plans:
        for coverage in COVERAGES:
            for arm in ARMS:
                k = ARM_K[arm]
                status = read_kv(os.path.join(heldout, plan, f"call.{coverage}.K{k}.mosaic_status.tsv"))
                complete = status.get("status") == "COMPLETE"
                if complete:
                    table = read_tsv(os.path.join(heldout, plan, f"eval.{coverage}.K{k}.blocks.tsv"))
                    genotypes = {
                        row["block_id"]: row
                        for row in read_tsv(
                            os.path.join(heldout, plan, f"call.{coverage}.K{k}.block_genotypes.tsv"))
                    }
                    source = "strict_call"
                else:
                    table = read_tsv(os.path.join(heldout, plan, f"tied.{coverage}.K{k}.blocks.tsv"))
                    genotypes = {}
                    source = "conditional_consensus"

                for index, row in enumerate(table):
                    kind = row.get("kind", NA)
                    distance = number(row.get("unordered_distance", row.get("distance")))
                    t1 = number(row.get("truth1_bp"))
                    t2 = number(row.get("truth2_bp"))
                    span = (t1 or 0) + (t2 or 0)
                    normalized = distance / span if (distance is not None and span > 0) else None
                    genotype_row = genotypes.get(row.get("block_id"), {})
                    support = genotype_row.get("genotype_support", row.get("genotype_support", NA))
                    rows.append([
                        plan, klass, coverage, arm, source, row.get("block_id", NA), index, kind,
                        cell(None if distance is None else int(distance)),
                        NA if distance is None else ("true" if distance == 0 else "false"),
                        cell(None if t1 is None else int(t1)),
                        cell(None if t2 is None else int(t2)),
                        NA if normalized is None else f"{normalized:.8f}",
                        support if support else NA,
                        genotype_row.get("conditional_GQ", NA),
                        # PRIOR_ONLY is not read support: it means the reads did not decide this block.
                        support if support in ("READS", "PRIOR_ONLY") else NA,
                    ])
    return rows


# --------------------------------------------------------------------------------------------------
# Tied-winner diagnostics -- NOT CALLS
# --------------------------------------------------------------------------------------------------

TIED_HEADER = [
    "plan", "class", "coverage", "arm", "label", "winner_set_size", "true_pair_in_tied_set",
    "any_tied_candidate_exact", "min_distance", "median_distance", "max_distance",
    "unanimous_blocks", "blocks", "exact_consensus_blocks", "consensus_total_distance",
    "consensus_max_block_distance", "consensus_sv_subtotal", "consensus_context_subtotal",
]


def build_tied(heldout, plans):
    rows = []
    for plan, klass in plans:
        for coverage in COVERAGES:
            for arm in ARMS:
                k = ARM_K[arm]
                summary = read_kv(os.path.join(heldout, plan, f"tied.{coverage}.K{k}.summary.tsv"))
                if not summary:
                    continue
                rows.append([
                    plan, klass, coverage, arm, summary.get("label", NA),
                    summary.get("winner_set_size", NA), summary.get("true_pair_in_tied_set", NA),
                    summary.get("any_tied_candidate_exact", NA),
                    summary.get("min_whole_locus_distance", NA),
                    summary.get("median_whole_locus_distance", NA),
                    summary.get("max_whole_locus_distance", NA),
                    summary.get("unanimous_blocks", NA), summary.get("blocks", NA),
                    summary.get("exact_consensus_blocks", NA),
                    summary.get("consensus_total_distance", NA),
                    summary.get("consensus_max_block_distance", NA),
                    summary.get("consensus_sv_subtotal", NA),
                    summary.get("consensus_context_subtotal", NA),
                ])
    return rows


# --------------------------------------------------------------------------------------------------
# Phase, denominators, ordinary LZO
# --------------------------------------------------------------------------------------------------

PHASE_CATEGORIES = ["READS_AUDITED", "READS_WITHIN_CANDIDATES", "PRIOR_ONLY", "UNRESOLVED"]


def build_phase(runs):
    rows = []
    for arm in ARMS:
        arm_runs = [r for r in runs if r[3] == arm]
        counts = {c: 0 for c in PHASE_CATEGORIES}
        other = 0
        for run in arm_runs:
            if run[12] in counts:
                counts[run[12]] += 1
            else:
                other += 1
        emitted = sum(1 for r in arm_runs if r[13] == "yes")
        for category in PHASE_CATEGORIES:
            rows.append([arm, category, counts[category], len(arm_runs)])
        if other:
            rows.append([arm, "OTHER", other, len(arm_runs)])
        rows.append([arm, "phased_haplotype_emitted", emitted, len(arm_runs)])
    return rows


def build_denominators(runs, tied, plans):
    controls = {plan for plan, klass in plans if klass == "panel_control"}
    rows = []
    for include_controls in (True, False):
        label = "all_plans" if include_controls else "excluding_panel_controls"
        for arm in ARMS:
            subset = [r for r in runs if r[3] == arm and (include_controls or r[0] not in controls)]
            intended = len(subset)
            complete = [r for r in subset if r[4] == "COMPLETE"]
            refused = [r for r in subset if r[4] != "COMPLETE"]
            exact = [r for r in complete if number(r[18]) == 0]
            blocks_exact = [number(r[17]) for r in complete if number(r[17]) is not None]
            distances = [number(r[18]) for r in complete if number(r[18]) is not None]
            tied_subset = [t for t in tied if t[3] == arm and (include_controls or t[0] not in controls)]
            consensus_exact = [t for t in tied_subset if number(t[14]) == 0]

            rows.append([label, arm, "intended_runs", intended, ""])
            rows.append([label, arm, "completed_runs", len(complete),
                         f"{100 * len(complete) / intended:.1f}%" if intended else NA])
            rows.append([label, arm, "no_call_runs", len(refused),
                         f"{100 * len(refused) / intended:.1f}%" if intended else NA])
            rows.append([label, arm, "exact_calls_over_intended", len(exact),
                         f"{100 * len(exact) / intended:.1f}%" if intended else NA])
            rows.append([label, arm, "exact_calls_given_a_call", len(exact),
                         f"{100 * len(exact) / len(complete):.1f}%" if complete else NA])
            rows.append([label, arm, "mean_exact_blocks_given_a_call",
                         f"{statistics.mean(blocks_exact):.2f}" if blocks_exact else NA,
                         f"over {len(blocks_exact)} calls x 11 blocks"])
            rows.append([label, arm, "total_distance_given_a_call",
                         int(sum(distances)) if distances else NA,
                         f"over {len(distances)} calls, {11 * len(distances)} blocks"])
            rows.append([label, arm, "refused_runs_with_exact_conditional_consensus",
                         len(consensus_exact),
                         f"of {len(tied_subset)} refused (DIAGNOSTIC, not calls)"])
    return rows


def build_ordinary(path):
    rows = read_tsv(path)
    if not rows:
        return [], []
    complete = [r for r in rows if r["status"] == "COMPLETE"]
    exact = [r for r in complete if number(r["unordered_distance"]) == 0]
    blocks_exact = [number(r["blocks_exact"]) for r in complete]
    sv_exact = [number(r["sv_blocks_exact"]) for r in complete]
    phase = {}
    for r in complete:
        phase[r["phase_support"]] = phase.get(r["phase_support"], 0) + 1
    summary = [
        ["donors_intended", len(rows), ""],
        ["completed_runs", len(complete), f"{100 * len(complete) / len(rows):.1f}%"],
        ["no_call_runs", len(rows) - len(complete), f"{100 * (len(rows) - len(complete)) / len(rows):.1f}%"],
        ["exact_calls_over_intended", len(exact), f"{100 * len(exact) / len(rows):.1f}%"],
        ["exact_calls_given_a_call", len(exact), f"{100 * len(exact) / len(complete):.1f}%"],
        ["mean_exact_blocks_given_a_call", f"{statistics.mean(blocks_exact):.2f}", "of 11"],
        ["mean_exact_sv_blocks_given_a_call", f"{statistics.mean(sv_exact):.2f}", "of 5"],
        ["total_distance_given_a_call",
         int(sum(number(r["unordered_distance"]) for r in complete)),
         f"over {len(complete)} calls"],
        ["total_sv_distance", int(sum(number(r["sv_distance"]) for r in complete)), ""],
        ["total_context_distance", int(sum(number(r["context_distance"]) for r in complete)), ""],
    ]
    for key in sorted(phase):
        summary.append([f"phase_{key}", phase[key], f"of {len(complete)} calls"])
    return summary, rows


# --------------------------------------------------------------------------------------------------
# Plots. Every figure writes the table it was drawn from.
# --------------------------------------------------------------------------------------------------

def block_order(blocks):
    seen = []
    for row in blocks:
        key = (int(row[6]), row[5], row[7])
        if key not in seen:
            seen.append(key)
    return sorted(seen)


def plot_block_distance(blocks, out, normalized=False):
    column = 12 if normalized else 8
    name = "plot_block_distance_normalized" if normalized else "plot_block_distance"
    title = "length-normalized phase-invariant diploid distance" if normalized else \
            "phase-invariant diploid edit distance"
    order = block_order(blocks)
    header = ["arm", "source", "block_index", "block_id", "kind", "value"]
    rows = []
    for row in blocks:
        value = number(row[column])
        if value is None:
            continue
        rows.append([row[3], row[4], row[6], row[5], row[7], value])
    write_tsv(os.path.join(out, name + ".tsv"), header, rows)

    figure, axes = plt.subplots(len(ARMS), 1, figsize=(12, 7), sharex=True)
    for axis, arm in zip(axes, ARMS):
        for source, filled in (("strict_call", True), ("conditional_consensus", False)):
            xs, ys = [], []
            for index, block_id, kind in order:
                values = [r[5] for r in rows if r[0] == arm and r[1] == source and r[2] == index]
                # Violins only where there are enough points for a distribution to mean anything.
                if len(values) >= 8 and not normalized:
                    body = axis.violinplot([values], positions=[index], widths=0.7,
                                           showextrema=False, showmedians=True)
                    for part in body["bodies"]:
                        part.set_facecolor(ARM_COLOUR[arm])
                        part.set_alpha(0.18 if filled else 0.08)
                xs.extend([index + (0.0 if filled else 0.22)] * len(values))
                ys.extend(values)
            axis.scatter(xs, ys, s=26,
                         facecolors=ARM_COLOUR[arm] if filled else "none",
                         edgecolors=ARM_COLOUR[arm], linewidths=1.2, zorder=3,
                         label="strict call" if filled else "conditional consensus (not a call)")
        axis.set_yscale("symlog", linthresh=1)
        axis.set_ylabel(f"{arm}\n{'normalized' if normalized else 'distance'}")
        axis.grid(alpha=0.25, axis="y")
        axis.legend(loc="upper left", fontsize=8, framealpha=0.9)
    labels = [f"{block_id}\n({kind})" for _, block_id, kind in order]
    axes[-1].set_xticks([index for index, _, _ in order])
    axes[-1].set_xticklabels(labels, fontsize=8)
    axes[0].set_title(f"C4 synthetic-mosaic LZO: {title}, by block (arms never pooled)")
    figure.tight_layout()
    figure.savefig(os.path.join(out, name + ".png"), dpi=150)
    plt.close(figure)


def plot_class_coverage(runs, tied, out):
    tied_by_key = {(t[0], t[2], t[3]): t for t in tied}
    counts = {}
    for run in runs:
        plan, klass, coverage, arm, status = run[0], run[1], run[2], run[3], run[4]
        if status == "COMPLETE":
            outcome = "exact strict call" if number(run[18]) == 0 else "inexact strict call"
        else:
            diagnostic = tied_by_key.get((plan, coverage, arm))
            exact = diagnostic is not None and number(diagnostic[14]) == 0
            outcome = ("no-call, exact conditional consensus" if exact
                       else "no-call, inexact conditional consensus")
        counts[(klass, coverage, arm, outcome)] = counts.get((klass, coverage, arm, outcome), 0) + 1

    outcomes = ["exact strict call", "inexact strict call",
                "no-call, exact conditional consensus", "no-call, inexact conditional consensus"]
    colours = ["#1A9850", "#FDAE61", "#74ADD1", "#D73027"]
    classes = sorted({r[1] for r in runs})
    header = ["class", "coverage", "arm", "outcome", "runs"]
    rows = [[k, c, a, o, n] for (k, c, a, o), n in sorted(counts.items())]
    write_tsv(os.path.join(out, "plot_class_coverage.tsv"), header, rows)

    figure, axes = plt.subplots(1, len(ARMS), figsize=(13, 5), sharey=True)
    for axis, arm in zip(axes, ARMS):
        labels, bottoms = [], []
        for klass in classes:
            for coverage in COVERAGES:
                labels.append(f"{klass}\n{coverage}x")
                bottoms.append(0)
        for outcome, colour in zip(outcomes, colours):
            heights, position = [], 0
            for klass in classes:
                for coverage in COVERAGES:
                    heights.append(counts.get((klass, coverage, arm, outcome), 0))
                    position += 1
            axis.bar(range(len(labels)), heights, bottom=bottoms, color=colour, label=outcome)
            bottoms = [b + h for b, h in zip(bottoms, heights)]
        axis.set_xticks(range(len(labels)))
        axis.set_xticklabels(labels, fontsize=6, rotation=90)
        axis.set_title(arm)
        axis.grid(alpha=0.25, axis="y")
    axes[0].set_ylabel("runs")
    axes[-1].legend(fontsize=7, loc="upper right")
    figure.suptitle("C4 mosaic LZO outcomes by class and coverage (calls and no-calls kept apart)")
    figure.tight_layout()
    figure.savefig(os.path.join(out, "plot_class_coverage.png"), dpi=150)
    plt.close(figure)


def plot_paired(runs, tied, out):
    by_key = {(r[0], r[2], r[3]): r for r in runs}
    tied_by_key = {(t[0], t[2], t[3]): t for t in tied}
    header = ["plan", "class", "coverage", "k0_distance", "k128_distance", "k128_outcome",
              "k128_consensus_distance"]
    rows = []
    for (plan, coverage, arm), run in sorted(by_key.items()):
        if arm != "K0":
            continue
        other = by_key.get((plan, coverage, "K128"))
        k0 = number(run[18])
        if other is not None and other[4] == "COMPLETE":
            rows.append([plan, run[1], coverage, cell(None if k0 is None else int(k0)),
                         cell(int(number(other[18]))), "call", NA])
        else:
            diagnostic = tied_by_key.get((plan, coverage, "K128"))
            consensus = number(diagnostic[14]) if diagnostic else None
            rows.append([plan, run[1], coverage, cell(None if k0 is None else int(k0)), NA,
                         "no_call", cell(None if consensus is None else int(consensus))])
    write_tsv(os.path.join(out, "plot_paired_k.tsv"), header, rows)

    figure, axis = plt.subplots(figsize=(9, 6))
    refused_y = 0.4   # a distinct band for refusals, so they cannot vanish from the figure
    for row in rows:
        k0 = number(row[3])
        if k0 is None:
            continue
        start = max(k0, refused_y)
        if row[5] == "call":
            end = max(number(row[4]), refused_y)
            axis.plot([0, 1], [start, end], color="#4D4D4D", alpha=0.5, zorder=1)
            axis.scatter([1], [end], color=ARM_COLOUR["K128"], s=32, zorder=3)
        else:
            axis.plot([0, 1], [start, refused_y / 2], color="#B2182B", alpha=0.5,
                      linestyle="--", zorder=1)
            axis.scatter([1], [refused_y / 2], facecolors="none", edgecolors="#B2182B", s=46, zorder=3)
        axis.scatter([0], [start], color=ARM_COLOUR["K0"], s=32, zorder=3)
    axis.set_yscale("log")
    axis.set_xticks([0, 1])
    axis.set_xticklabels(["K=0 (panel arm)", "K=128 (mosaic arm)"])
    axis.axhline(refused_y, color="#999999", linewidth=0.8)
    axis.text(1.02, refused_y / 2, "NO CALL\n(refused)", fontsize=8, color="#B2182B", va="center")
    axis.set_ylabel("whole-locus phase-invariant distance (log; refusals in their own band)")
    axis.set_title("Paired K=0 to K=128, same plan and coverage")
    axis.grid(alpha=0.25, axis="y")
    figure.tight_layout()
    figure.savefig(os.path.join(out, "plot_paired_k.png"), dpi=150)
    plt.close(figure)


def plot_tied(tied, out):
    header = ["plan", "coverage", "winner_set_size", "min", "median", "max", "consensus"]
    rows = [[t[0], t[2], int(t[5]), int(t[8]), int(t[9]), int(t[10]), int(t[14])] for t in tied]
    rows.sort()
    write_tsv(os.path.join(out, "plot_tied_winner.tsv"), header, rows)
    if not rows:
        return
    figure, axis = plt.subplots(figsize=(11, 5))
    positions = range(len(rows))
    for position, row in zip(positions, rows):
        axis.vlines(position, max(row[3], 0.4), max(row[5], 0.4), color="#999999", linewidth=2)
        axis.scatter([position], [max(row[4], 0.4)], color="#4D4D4D", s=28, label="median" if position == 0 else "")
        axis.scatter([position], [max(row[6], 0.4)], facecolors="none", edgecolors="#B2182B", s=52,
                     label="conditional consensus" if position == 0 else "")
    axis.set_xticks(list(positions))
    axis.set_xticklabels([f"{r[0]}\n{r[1]}x  n={r[2]}" for r in rows], fontsize=6, rotation=90)
    axis.set_yscale("symlog", linthresh=1)
    axis.set_ylabel("distance to truth")
    axis.set_title("Refused runs: tied-candidate distance range and the consensus "
                   "(CONDITIONAL_CONSENSUS_DIAGNOSTIC, not calls)")
    axis.legend(fontsize=8)
    axis.grid(alpha=0.25, axis="y")
    figure.tight_layout()
    figure.savefig(os.path.join(out, "plot_tied_winner.png"), dpi=150)
    plt.close(figure)


def plot_gq(runs, out):
    header = ["arm", "outcome", "coverage", "plan", "min_block_gq"]
    rows = []
    for run in runs:
        gq = number(run[16])
        if gq is None:
            continue
        outcome = "exact call" if (run[4] == "COMPLETE" and number(run[18]) == 0) else (
            "inexact call" if run[4] == "COMPLETE" else "no call")
        rows.append([run[3], outcome, run[2], run[0], gq])
    write_tsv(os.path.join(out, "plot_gq.tsv"), header, rows)

    figure, axes = plt.subplots(1, len(ARMS), figsize=(11, 5), sharey=True)
    outcomes = ["exact call", "inexact call", "no call"]
    for axis, arm in zip(axes, ARMS):
        for index, outcome in enumerate(outcomes):
            values = [r[4] for r in rows if r[0] == arm and r[1] == outcome]
            if not values:
                continue
            axis.scatter([index] * len(values), values, s=30, alpha=0.75,
                         color=ARM_COLOUR[arm])
        axis.set_xticks(range(len(outcomes)))
        axis.set_xticklabels(outcomes, fontsize=8)
        axis.set_yscale("log")
        axis.set_title(arm)
        axis.grid(alpha=0.25, axis="y")
    axes[0].set_ylabel("minimum block conditional GQ (log)")
    figure.suptitle("Conditional, UNCALIBRATED GQ by arm and outcome -- arms never pooled")
    figure.tight_layout()
    figure.savefig(os.path.join(out, "plot_gq.png"), dpi=150)
    plt.close(figure)


def plot_runtime(heldout, out):
    rows = []
    for row in read_tsv(os.path.join(heldout, "runtime.tsv")):
        rows.append([row["plan"], row["coverage"], float(row["wall_s"]),
                     int(row["max_rss_bytes"]) / 1073741824.0])
    write_tsv(os.path.join(out, "plot_runtime.tsv"),
              ["plan", "coverage", "wall_s", "peak_rss_gb"], rows)
    figure, axes = plt.subplots(1, 2, figsize=(11, 4))
    for axis, index, label in ((axes[0], 2, "wall seconds"), (axes[1], 3, "peak RSS (GB)")):
        for position, coverage in enumerate(COVERAGES):
            values = [r[index] for r in rows if r[1] == coverage]
            axis.scatter([position] * len(values), values, s=30, color="#2166AC", alpha=0.75)
        axis.set_xticks(range(len(COVERAGES)))
        axis.set_xticklabels([f"{c}x" for c in COVERAGES])
        axis.set_ylabel(label)
        axis.grid(alpha=0.25, axis="y")
    figure.suptitle("Shared scoring pass per plan and coverage (one run serves both arms)")
    figure.tight_layout()
    figure.savefig(os.path.join(out, "plot_runtime.png"), dpi=150)
    plt.close(figure)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--heldout", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--ordinary", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    os.makedirs(args.out, exist_ok=True)

    plans = plans_of(args.manifest)
    runs = build_runs(args.heldout, plans)
    blocks = build_blocks(args.heldout, plans)
    tied = build_tied(args.heldout, plans)

    write_tsv(os.path.join(args.out, "runs.tsv"), RUN_HEADER, runs)
    write_tsv(os.path.join(args.out, "blocks.tsv"), BLOCK_HEADER, blocks)
    write_tsv(os.path.join(args.out, "tied.tsv"), TIED_HEADER, tied)
    write_tsv(os.path.join(args.out, "phase.tsv"), ["arm", "category", "runs", "of"],
              build_phase(runs))
    write_tsv(os.path.join(args.out, "denominators.tsv"),
              ["population", "arm", "quantity", "value", "note"],
              build_denominators(runs, tied, plans))
    ordinary_summary, ordinary_rows = build_ordinary(args.ordinary)
    write_tsv(os.path.join(args.out, "ordinary_lzo.tsv"), ["quantity", "value", "note"],
              ordinary_summary)

    plot_block_distance(blocks, args.out, normalized=False)
    plot_block_distance(blocks, args.out, normalized=True)
    plot_class_coverage(runs, tied, args.out)
    plot_paired(runs, tied, args.out)
    plot_tied(tied, args.out)
    plot_gq(runs, args.out)
    plot_runtime(args.heldout, args.out)

    print(f"runs={len(runs)} blocks={len(blocks)} tied_diagnostics={len(tied)} "
          f"ordinary_donors={len(ordinary_rows)}")
    print(f"written to {args.out}")


if __name__ == "__main__":
    main()
