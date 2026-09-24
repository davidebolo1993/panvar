#!/usr/bin/env python3
"""C4 production-v1: result tables and figures by mosaic switch class and by block, from the frozen package.

    scripts/c4_production_v1_analysis.py --package <campaign>/package --out <campaign>/analysis_v1

Reads only the completed production-v1 package (runs.tsv, blocks.tsv, fallback_diagnostic.tsv, summary.tsv)
and the frozen held-out mosaic manifest; genotypes nothing and never writes into the package. Every
aggregate that summary.tsv also carries is recomputed here from runs.tsv and blocks.tsv and compared; any
mismatch is a hard failure (crosscheck.tsv records every comparison).

Conventions, kept everywhere:
  * a no-call is not an observation of distance: it is counted, by reason, in its own strip or column and
    never enters a called-distance distribution or a distance of zero;
  * distance is the phase-invariant unordered genotype edit distance; excess is distance above the
    availability floor (the LOO quantity to read, since a raw distance of zero is often impossible);
  * a block whose truth is undefined (not spellable in the model's frame) is neither exact nor inexact: it
    is counted as truth_undefined and excluded from block rates.

Writes mosaic_by_class.tsv, block_accuracy.tsv, crosscheck.tsv and the figures
c4_mosaic_by_switch_class_callability, c4_mosaic_by_switch_class_distance, c4_production_block_exact_heatmap,
c4_production_block_distance_by_id and c4_mosaic_block_exact_by_class (each .png and .pdf).
"""
import argparse
import pathlib
import sys

import numpy as np
import pandas as pd

REPO = pathlib.Path(__file__).resolve().parent.parent
MANIFEST = REPO / "tests" / "genotype_c4_heldout_mosaic_manifest.tsv"
COVERAGES = [2, 5, 10, 30]
EXPERIMENTS = [("lzo", "Ordinary LZO"), ("mosaic", "Synthetic-mosaic LZO"), ("loo", "LOO")]
CLASSES = ["panel_control", "one_homologue_once", "both_homologues", "three_switches", "dense_switching"]
BLOCKS = ["ctx:0", "sv:1", "ctx:2", "sv:3", "ctx:4", "sv:5", "ctx:6", "sv:7", "ctx:8", "sv:9", "ctx:10"]
CALLS = ("strict_call", "genotype_equivalent_call")
REASON_COLOUR = {"improving_move_tie": "#E08214", "union_scoring_cap": "#4D4D4D", "round_cap": "#8073AC"}
CALL_COLOUR = {"strict_call": "#2166AC", "genotype_equivalent_call": "#92C5DE"}


class Mismatch(RuntimeError):
    pass


def plan_classes():
    rows = [line.rstrip("\n").split("\t") for line in open(MANIFEST) if line.startswith("plan\t")]
    header = rows[0]
    table = pd.DataFrame(rows[1:], columns=header)
    return dict(zip(table["name"], table["class"]))


def load(package):
    runs = pd.read_csv(package / "runs.tsv", sep="\t")
    blocks = pd.read_csv(package / "blocks.tsv", sep="\t")
    fallback = pd.read_csv(package / "fallback_diagnostic.tsv", sep="\t")
    summary = pd.read_csv(package / "summary.tsv", sep="\t")
    classes = plan_classes()
    for frame in (runs, blocks, fallback):
        frame["mosaic_class"] = np.where(frame.experiment == "mosaic", frame.subject.map(classes), "")
    if runs[(runs.experiment == "mosaic") & runs.mosaic_class.isna()].shape[0]:
        raise Mismatch("a mosaic plan has no class in the frozen manifest")
    runs["called"] = runs.outcome.isin(CALLS)
    return runs, blocks, fallback, summary


# ---- the cross-check against summary.tsv --------------------------------------------------------

def crosscheck(runs, blocks, summary, out):
    records = []

    def compare(experiment, coverage, key, mine, theirs):
        both_nan = (isinstance(mine, float) and np.isnan(mine)) and (isinstance(theirs, float) and np.isnan(theirs))
        ok = both_nan or (np.isclose(float(mine), float(theirs), rtol=1e-12, atol=1e-9))
        records.append({"experiment": experiment, "coverage": coverage, "quantity": key,
                        "recomputed": mine, "summary": theirs, "agrees": ok})

    # Block-level totals must add up to the run-level totals the summary was built from.
    for job, rows in blocks.groupby("job"):
        run = runs.loc[runs.job == job].iloc[0]
        d = pd.to_numeric(rows.unordered_distance, errors="coerce")
        if not np.isclose(d.sum(), run.distance) or (d == 0).sum() != run.blocks_exact:
            raise Mismatch(f"{job}: block distances or exact blocks do not add up to the run's")
        sv = rows.kind == "SV"
        if (d[sv] == 0).sum() != run.sv_blocks_exact:
            raise Mismatch(f"{job}: exact SV blocks do not add up")
    for row in summary.itertuples(index=False):
        part = runs[(runs.experiment == row.experiment) & (runs.coverage == row.coverage)]
        calls = part[part.called]
        cb = blocks[(blocks.experiment == row.experiment) & (blocks.coverage == row.coverage)]
        d = pd.to_numeric(cb.unordered_distance, errors="coerce")
        n = len(part)
        compare(row.experiment, row.coverage, "intended", n, row.intended)
        compare(row.experiment, row.coverage, "strict_calls", (part.outcome == "strict_call").sum(), row.strict_calls)
        compare(row.experiment, row.coverage, "genotype_equivalent_calls",
                (part.outcome == "genotype_equivalent_call").sum(), row.genotype_equivalent_calls)
        compare(row.experiment, row.coverage, "no_calls", (part.outcome == "no_call").sum(), row.no_calls)
        for reason in ("improving_move_tie", "round_cap", "union_scoring_cap"):
            compare(row.experiment, row.coverage, f"no_call_{reason}", (part.no_call_reason == reason).sum(),
                    getattr(row, f"no_call_{reason}"))
        compare(row.experiment, row.coverage, "exact_genotype_rate_all_intended", part.exact.sum() / n,
                row.exact_genotype_rate_all_intended)
        compare(row.experiment, row.coverage, "exact_genotype_rate_given_call",
                calls.exact.sum() / len(calls) if len(calls) else np.nan, row.exact_genotype_rate_given_call)
        compare(row.experiment, row.coverage, "floor_attained_rate_given_call",
                float((calls.excess == 0).mean()) if len(calls) else np.nan, row.floor_attained_rate_given_call)
        # From BLOCKS, not runs: the block table must reproduce the run-level block rates.
        compare(row.experiment, row.coverage, "exact_blocks_given_call",
                (d == 0).sum() / (len(calls) * len(BLOCKS)) if len(calls) else np.nan, row.exact_blocks_given_call)
        sv = cb.kind == "SV"
        compare(row.experiment, row.coverage, "exact_sv_blocks_given_call",
                (d[sv] == 0).sum() / (len(calls) * 5) if len(calls) else np.nan, row.exact_sv_blocks_given_call)
        compare(row.experiment, row.coverage, "distance_total_given_call", d.sum() if len(calls) else np.nan,
                row.distance_total_given_call)
        compare(row.experiment, row.coverage, "excess_total_given_call",
                pd.to_numeric(cb.excess, errors="coerce").sum() if len(calls) else np.nan, row.excess_total_given_call)
        for support in ("READS_WITHIN_CANDIDATES", "PRIOR_ONLY", "UNRESOLVED"):
            compare(row.experiment, row.coverage, f"phase_{support}", (calls.phase_support == support).sum(),
                    getattr(row, f"phase_{support}"))
    table = pd.DataFrame(records)
    table.to_csv(out / "crosscheck.tsv", sep="\t", index=False)
    bad = table[~table.agrees]
    if len(bad):
        raise Mismatch("recomputed aggregates disagree with summary.tsv:\n" + bad.to_string(index=False))
    return len(table)


# ---- tables ------------------------------------------------------------------------------------

def mosaic_by_class(runs, out):
    rows = []
    mosaic = runs[runs.experiment == "mosaic"]
    for cls in CLASSES:
        for coverage in COVERAGES:
            part = mosaic[(mosaic.mosaic_class == cls) & (mosaic.coverage == coverage)]
            calls = part[part.called]
            row = {"mosaic_class": cls, "coverage": coverage, "intended_plans": len(part),
                   "plans": ",".join(sorted(part.subject)),
                   "strict_calls": int((part.outcome == "strict_call").sum()),
                   "genotype_equivalent_calls": int((part.outcome == "genotype_equivalent_call").sum()),
                   "no_calls": int((part.outcome == "no_call").sum())}
            for reason in REASON_COLOUR:
                row[f"no_call_{reason}"] = int((part.no_call_reason == reason).sum())
            row["exact_calls"] = int(calls.exact.sum())
            row["exact_over_intended"] = f"{int(calls.exact.sum())}/{len(part)}"
            row["exact_given_call"] = f"{int(calls.exact.sum())}/{len(calls)}" if len(calls) else "0/0"
            row["floor_attained_given_call"] = f"{int((calls.excess == 0).sum())}/{len(calls)}" if len(calls) else "0/0"
            for quantity in ("distance", "excess"):
                values = calls[quantity]
                row[f"{quantity}_median_given_call"] = values.median() if len(values) else np.nan
                row[f"{quantity}_max_given_call"] = values.max() if len(values) else np.nan
            rows.append(row)
    table = pd.DataFrame(rows)
    table.to_csv(out / "mosaic_by_class.tsv", sep="\t", index=False)
    return table


def block_accuracy(runs, blocks, out):
    rows = []
    groups = [(e, "", c) for e, _ in EXPERIMENTS for c in COVERAGES]
    groups += [("mosaic", cls, c) for cls in CLASSES for c in COVERAGES]
    for experiment, cls, coverage in groups:
        part = runs[(runs.experiment == experiment) & (runs.coverage == coverage)]
        cb = blocks[(blocks.experiment == experiment) & (blocks.coverage == coverage)]
        if cls:
            part = part[part.mosaic_class == cls]
            cb = cb[cb.mosaic_class == cls]
        called = int(part.called.sum())
        for block in BLOCKS:
            here = cb[cb.block_id == block]
            d = pd.to_numeric(here.unordered_distance, errors="coerce")
            x = pd.to_numeric(here.excess, errors="coerce")
            defined = d.dropna()
            exact = int((defined == 0).sum())
            rows.append({
                "experiment": experiment, "coverage": coverage, "mosaic_class": cls or "all",
                "block_id": block, "kind": "SV" if block.startswith("sv") else "context",
                "intended_jobs": len(part), "called_jobs": called, "truth_undefined": int(d.isna().sum()),
                "exact_called_blocks": exact,
                "exact_rate_given_call": exact / len(defined) if len(defined) else np.nan,
                # No-calls are unavailable, not zero-distance: an exact block needs a call.
                "exact_rate_over_intended": exact / len(part) if len(part) else np.nan,
                "distance_median_given_call": defined.median() if len(defined) else np.nan,
                "distance_mean_given_call": defined.mean() if len(defined) else np.nan,
                "distance_max_given_call": defined.max() if len(defined) else np.nan,
                "excess_median_given_call": x.dropna().median() if len(defined) else np.nan,
                "excess_mean_given_call": x.dropna().mean() if len(defined) else np.nan,
                "excess_max_given_call": x.dropna().max() if len(defined) else np.nan,
            })
    table = pd.DataFrame(rows)
    table.to_csv(out / "block_accuracy.tsv", sep="\t", index=False)
    return table


# ---- figures -----------------------------------------------------------------------------------

def save(fig, out, name):
    for suffix in ("png", "pdf"):
        fig.savefig(out / f"{name}.{suffix}", dpi=160)


def symlog(axis, top):
    axis.set_yscale("symlog", linthresh=1.0, linscale=0.6)
    axis.set_ylim(-0.3, max(10.0, top) * 2.0)
    axis.axhline(0, color="#BBBBBB", linewidth=0.6, zorder=0)


def mosaic_callability(runs, table, out, plt):
    fig, axes = plt.subplots(1, len(CLASSES), figsize=(17, 4.6), sharey=True)
    for axis, cls in zip(axes, CLASSES):
        rows = table[table.mosaic_class == cls].set_index("coverage")
        x = np.arange(len(COVERAGES))
        bottom = np.zeros(len(COVERAGES))
        for key, colour, label in (("strict_calls", CALL_COLOUR["strict_call"], "strict call"),
                                   ("genotype_equivalent_calls", CALL_COLOUR["genotype_equivalent_call"],
                                    "genotype-equivalent call")):
            counts = np.array([rows.loc[c, key] for c in COVERAGES])
            axis.bar(x, counts, bottom=bottom, width=0.6, color=colour, label=label)
            bottom += counts
        for reason, colour in REASON_COLOUR.items():
            counts = np.array([rows.loc[c, f"no_call_{reason}"] for c in COVERAGES])
            axis.bar(x, counts, bottom=bottom, width=0.6, color=colour, hatch="//", edgecolor="white",
                     label=f"no-call: {reason.replace('_', ' ')}")
            bottom += counts
        for i, c in enumerate(COVERAGES):
            axis.text(i, 2.12, f"exact {rows.loc[c, 'exact_over_intended']}\nper call {rows.loc[c, 'exact_given_call']}",
                      ha="center", va="bottom", fontsize=7.5)
        # The two plans, named: with two per class a distribution would be pretence.
        plans = runs[(runs.experiment == "mosaic") & (runs.mosaic_class == cls)]
        names = sorted(plans.subject.unique())
        for j, name in enumerate(names):
            for i, c in enumerate(COVERAGES):
                outcome = plans[(plans.subject == name) & (plans.coverage == c)].outcome.iloc[0]
                symbol = {"strict_call": "o", "genotype_equivalent_call": "^", "no_call": "x"}[outcome]
                axis.scatter(i - 0.12 + 0.24 * j, j + 0.5, marker=symbol, color="black", s=18, zorder=5)
        axis.set_title(cls.replace("_", " "))
        axis.set_xticks(x)
        axis.set_xticklabels([f"{c}x" for c in COVERAGES])
        axis.set_ylim(0, 2.9)
        axis.set_yticks([0, 1, 2])
    axes[0].set_ylabel("plans (two per class)")
    handles, labels = axes[0].get_legend_handles_labels()
    handles += [plt.Line2D([], [], color="black", marker=m, linestyle="", label=l)
                for m, l in (("o", "plan: strict call"), ("^", "plan: genotype-equivalent"), ("x", "plan: no-call"))]
    fig.legend(handles=handles, loc="lower center", ncol=5, frameon=False, fontsize=8)
    fig.suptitle("C4 production-v1, synthetic mosaics: callability by switch class (left marker = plan _1, right = plan _2)")
    fig.subplots_adjust(bottom=0.24, top=0.86, wspace=0.08)
    save(fig, out, "c4_mosaic_by_switch_class_callability")
    plt.close(fig)


def mosaic_distance(runs, out, plt):
    mosaic = runs[runs.experiment == "mosaic"]
    top = float(mosaic[mosaic.called][["distance", "excess"]].max().max() or 1)
    fig = plt.figure(figsize=(17, 6.2))
    grid = fig.add_gridspec(3, len(CLASSES), height_ratios=[3, 3, 1.1], hspace=0.12, wspace=0.08)
    for k, cls in enumerate(CLASSES):
        plans = mosaic[mosaic.mosaic_class == cls]
        names = sorted(plans.subject.unique())
        for r, quantity in enumerate(("distance", "excess")):
            axis = fig.add_subplot(grid[r, k])
            for j, name in enumerate(names):
                for i, c in enumerate(COVERAGES):
                    run = plans[(plans.subject == name) & (plans.coverage == c)].iloc[0]
                    if not run.called:
                        continue
                    axis.scatter(i - 0.12 + 0.24 * j, run[quantity],
                                 marker="o" if run.outcome == "strict_call" else "^",
                                 color=["#1B7837", "#762A83"][j], s=26, zorder=5)
            symlog(axis, top)
            axis.set_xlim(-0.5, len(COVERAGES) - 0.5)
            axis.set_xticks(range(len(COVERAGES)))
            axis.set_xticklabels([])
            if k == 0:
                axis.set_ylabel(f"{quantity} (bp, symlog)")
            else:
                axis.set_yticklabels([])
            if r == 0:
                axis.set_title(cls.replace("_", " "))
        strip = fig.add_subplot(grid[2, k])
        for j, name in enumerate(names):
            for i, c in enumerate(COVERAGES):
                run = plans[(plans.subject == name) & (plans.coverage == c)].iloc[0]
                if run.called:
                    continue
                strip.scatter(i - 0.12 + 0.24 * j, 0.5, marker="s", s=40,
                              color=REASON_COLOUR.get(run.no_call_reason, "#999999"))
        strip.set_ylim(0, 1)
        strip.set_xlim(-0.5, len(COVERAGES) - 0.5)
        strip.set_yticks([])
        strip.set_xticks(range(len(COVERAGES)))
        strip.set_xticklabels([f"{c}x" for c in COVERAGES])
        if k == 0:
            strip.set_ylabel("no-calls", rotation=0, ha="right", va="center")
    handles = [plt.Line2D([], [], color=c, marker="o", linestyle="", label=f"plan _{j + 1}")
               for j, c in enumerate(("#1B7837", "#762A83"))]
    handles += [plt.Line2D([], [], color="#555555", marker=m, linestyle="", label=l)
                for m, l in (("o", "strict call"), ("^", "genotype-equivalent call"))]
    handles += [plt.Rectangle((0, 0), 1, 1, color=c, label=f"no-call: {r.replace('_', ' ')}") for r, c in REASON_COLOUR.items()]
    fig.legend(handles=handles, loc="lower center", ncol=7, frameon=False, fontsize=8)
    fig.suptitle("C4 production-v1, synthetic mosaics: distance and excess over the availability floor, calls only")
    fig.subplots_adjust(bottom=0.13, top=0.9)
    save(fig, out, "c4_mosaic_by_switch_class_distance")
    plt.close(fig)


def exact_heatmap(table, out, plt):
    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list("exact", ["#B2182B", "#FDDBC7", "#F7F7F7", "#D1E5F0", "#2166AC"])
    fig, axes = plt.subplots(len(EXPERIMENTS), 1, figsize=(15, 9.5))
    for axis, (experiment, title) in zip(axes, EXPERIMENTS):
        part = table[(table.experiment == experiment) & (table.mosaic_class == "all")]
        grid = np.full((len(COVERAGES), len(BLOCKS)), np.nan)
        for i, c in enumerate(COVERAGES):
            for j, b in enumerate(BLOCKS):
                row = part[(part.coverage == c) & (part.block_id == b)].iloc[0]
                grid[i, j] = row.exact_rate_given_call
                label = ("no calls" if row.called_jobs == 0 else f"{row.exact_rate_given_call:.0%}") + \
                        f"\n{row.called_jobs}/{row.intended_jobs}"
                if row.truth_undefined:
                    label += f"\n({row.truth_undefined} undef.)"
                if row.called_jobs == 0:
                    axis.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, facecolor="#DDDDDD", hatch="xx",
                                                 edgecolor="#999999"))
                axis.text(j, i, label, ha="center", va="center", fontsize=7)
        image = axis.imshow(np.ma.masked_invalid(grid), cmap=cmap, vmin=0, vmax=1, aspect="auto")
        axis.set_xticks(range(len(BLOCKS)))
        axis.set_xticklabels(BLOCKS)
        axis.set_yticks(range(len(COVERAGES)))
        axis.set_yticklabels([f"{c}x" for c in COVERAGES])
        axis.set_title(f"{title}: exact-block rate given a call (cell: rate, called/intended)", fontsize=10)
    fig.colorbar(image, ax=axes, fraction=0.02, pad=0.01, label="exact-block rate given a call")
    save(fig, out, "c4_production_block_exact_heatmap")
    plt.close(fig)


def distance_by_id(runs, blocks, out, plt):
    top = float(pd.to_numeric(blocks.unordered_distance, errors="coerce").max())
    fig = plt.figure(figsize=(17, 11))
    grid = fig.add_gridspec(len(EXPERIMENTS) * 2, 1, height_ratios=[4, 0.7] * len(EXPERIMENTS), hspace=0.35)
    rng = np.random.default_rng(20260924)
    offsets = {c: (k - 1.5) * 0.19 for k, c in enumerate(COVERAGES)}
    colours = {2: "#FDB863", 5: "#E08214", 10: "#8073AC", 30: "#2D004B"}
    for e, (experiment, title) in enumerate(EXPERIMENTS):
        axis = fig.add_subplot(grid[2 * e])
        part = blocks[blocks.experiment == experiment]
        for j, b in enumerate(BLOCKS):
            if b.startswith("sv"):
                axis.axvspan(j - 0.5, j + 0.5, color="#F4E1E1", zorder=0)
            for c in COVERAGES:
                d = pd.to_numeric(part[(part.block_id == b) & (part.coverage == c)].unordered_distance,
                                  errors="coerce").dropna().to_numpy()
                x = j + offsets[c] + rng.uniform(-0.06, 0.06, size=len(d))
                axis.scatter(x, d, s=5, color=colours[c], alpha=0.6, linewidths=0, zorder=3)
        symlog(axis, top)
        axis.set_xlim(-0.5, len(BLOCKS) - 0.5)
        axis.set_xticks(range(len(BLOCKS)))
        axis.set_xticklabels(BLOCKS)
        axis.set_ylabel("block distance (bp, symlog)")
        axis.set_title(f"{title}: every called block, by block ID (SV blocks shaded)", fontsize=10)
        strip = fig.add_subplot(grid[2 * e + 1])
        sub = runs[runs.experiment == experiment]
        for k, c in enumerate(COVERAGES):
            made = 0
            for reason, colour in REASON_COLOUR.items():
                count = int(((sub.coverage == c) & (sub.no_call_reason == reason)).sum())
                strip.bar(k, count, bottom=made, width=0.5, color=colour)
                made += count
            strip.text(k, made + 0.3, f"{made}/{int((sub.coverage == c).sum())} no-call", ha="center", fontsize=7)
        strip.set_xticks(range(len(COVERAGES)))
        strip.set_xticklabels([f"{c}x" for c in COVERAGES])
        strip.set_ylim(0, max(4, int((sub.outcome == 'no_call').groupby(sub.coverage).sum().max() or 0) * 1.6))
        strip.set_ylabel("no-calls", fontsize=8)
    handles = [plt.Line2D([], [], color=colours[c], marker="o", linestyle="", label=f"{c}x") for c in COVERAGES]
    handles += [plt.Rectangle((0, 0), 1, 1, color=c, label=f"no-call: {r.replace('_', ' ')}") for r, c in REASON_COLOUR.items()]
    fig.legend(handles=handles, loc="lower center", ncol=7, frameon=False, fontsize=8)
    fig.suptitle("C4 production-v1: phase-invariant block distance from truth by block ID (calls only; no-calls below)")
    fig.subplots_adjust(bottom=0.05, top=0.95)
    save(fig, out, "c4_production_block_distance_by_id")
    plt.close(fig)


def mosaic_block_heatmap(table, out, plt):
    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list("exact", ["#B2182B", "#FDDBC7", "#F7F7F7", "#D1E5F0", "#2166AC"])
    part = table[(table.experiment == "mosaic") & (table.mosaic_class != "all")]
    labels = [f"{cls.replace('_', ' ')} {c}x" for cls in CLASSES for c in COVERAGES]
    grid = np.full((len(labels), len(BLOCKS)), np.nan)
    fig, axis = plt.subplots(figsize=(13, 10))
    for i, (cls, c) in enumerate([(cls, c) for cls in CLASSES for c in COVERAGES]):
        for j, b in enumerate(BLOCKS):
            row = part[(part.mosaic_class == cls) & (part.coverage == c) & (part.block_id == b)].iloc[0]
            if row.called_jobs == 0:
                axis.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, facecolor="#DDDDDD", hatch="xx",
                                             edgecolor="#999999"))
                axis.text(j, i, f"0/{row.intended_jobs}", ha="center", va="center", fontsize=7)
                continue
            grid[i, j] = row.exact_rate_given_call
            axis.text(j, i, f"{row.exact_called_blocks}/{row.called_jobs} of {row.intended_jobs}",
                      ha="center", va="center", fontsize=6.5)
    axis.imshow(np.ma.masked_invalid(grid), cmap=cmap, vmin=0, vmax=1, aspect="auto")
    axis.set_xticks(range(len(BLOCKS)))
    axis.set_xticklabels(BLOCKS)
    axis.set_yticks(range(len(labels)))
    axis.set_yticklabels(labels, fontsize=8)
    for boundary in range(len(COVERAGES), len(labels), len(COVERAGES)):
        axis.axhline(boundary - 0.5, color="black", linewidth=0.8)
    axis.set_title("Synthetic mosaics by switch class: exact blocks / called plans of intended plans\n"
                   "(two plans per class: raw counts, not rates to generalise from; hatched = no call)", fontsize=10)
    fig.tight_layout()
    save(fig, out, "c4_mosaic_block_exact_by_class")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--package", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    package = pathlib.Path(args.package).resolve()
    out = pathlib.Path(args.out).resolve()
    if out == package or package in out.parents and out.name == package.name:
        sys.exit("refusing to write into the frozen package")
    out.mkdir(parents=True, exist_ok=True)
    runs, blocks, fallback, summary = load(package)
    try:
        checked = crosscheck(runs, blocks, summary, out)
    except Mismatch as mismatch:
        print(f"HARD FAILURE: {mismatch}")
        return 1
    classes = mosaic_by_class(runs, out)
    accuracy = block_accuracy(runs, blocks, out)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    mosaic_callability(runs, classes, out, plt)
    mosaic_distance(runs, out, plt)
    exact_heatmap(accuracy, out, plt)
    distance_by_id(runs, blocks, out, plt)
    mosaic_block_heatmap(accuracy, out, plt)
    print(f"{checked} aggregates agree with summary.tsv; tables and figures written to {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
