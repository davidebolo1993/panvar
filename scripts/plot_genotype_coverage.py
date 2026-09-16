#!/usr/bin/env python3
"""Per-block true edit distance by coverage, one row per arm (review round 8).

    scripts/plot_genotype_coverage.py --blocks <blocks.tsv> --runs <runs.tsv> --out <dir> [--dpi 150]

Input is the table written by tests/genotype_coverage_c4.sh: one row per donor, arm, coverage and block, with
the unordered diploid edlib distance to truth, whether it is exact, the conditional GQ, and the genotype and
phase support.

The figure is one violin per block and coverage, panel arm above and mosaic arm below, of the true unordered
diploid edit distance. At 30x the violins collapse to zero, which is the milestone result; the lower coverages
are where the arms can separate. Underneath, the no-call rate (runs that produced no genotype) and the exact
rate are drawn per coverage, per arm, because a distance violin cannot show a run that never called.

Arms are never pooled: the panel arm conditions on every panel diplotype, the mosaic arm on those plus its
mosaics, so their conditional quantities answer different questions.
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ARMS = ["panel", "mosaic"]
COLOUR = {"panel": "#2166AC", "mosaic": "#B2182B"}


def violins(axis, frame, blocks, coverages, colour):
    positions = []
    data = []
    for block_index, block in enumerate(blocks):
        for coverage_index, coverage in enumerate(coverages):
            values = frame[(frame.block == block) & (frame.coverage == coverage)].unordered_diploid_edlib
            values = pd.to_numeric(values, errors="coerce").dropna().to_numpy()
            if values.size == 0:
                continue
            positions.append(block_index + (coverage_index - (len(coverages) - 1) / 2) * 0.13)
            data.append(values if values.size > 1 else np.repeat(values, 2))
    if data:
        parts = axis.violinplot(data, positions=positions, widths=0.12, showextrema=False, showmedians=True)
        for body in parts["bodies"]:
            body.set_facecolor(colour)
            body.set_alpha(0.65)
            body.set_edgecolor("none")
        parts["cmedians"].set_color("#14201f")
        parts["cmedians"].set_linewidth(0.8)
    axis.set_xticks(range(len(blocks)))
    axis.set_xticklabels(blocks, rotation=45, ha="right", fontsize=8)
    axis.set_ylabel("unordered diploid\nedit distance to truth")
    axis.grid(axis="y", color="#e6ecec", linewidth=0.8)
    axis.set_axisbelow(True)


def rates(axis, runs, coverages):
    width = 0.35
    for index, arm in enumerate(ARMS):
        mine = runs[runs.arm == arm]
        no_call = [float((mine[mine.coverage == c].run_no_call == "true").mean() or 0.0) for c in coverages]
        called = mine[mine.run_no_call == "false"]
        exact = []
        for coverage in coverages:
            rows = called[called.coverage == coverage]
            exact.append(float((pd.to_numeric(rows.unordered_distance, errors="coerce") == 0).mean()) if len(rows) else 0.0)
        offset = (index - 0.5) * width
        axis.bar(np.arange(len(coverages)) + offset, exact, width, color=COLOUR[arm], alpha=0.8,
                 label="%s: donors fully exact" % arm)
        axis.bar(np.arange(len(coverages)) + offset, no_call, width, bottom=exact, color=COLOUR[arm], alpha=0.35,
                 hatch="//", label="%s: no-call" % arm)
    axis.set_xticks(range(len(coverages)))
    axis.set_xticklabels(["%gx" % c for c in coverages])
    axis.set_ylim(0, 1.05)
    axis.set_ylabel("rate over donors")
    axis.set_xlabel("coverage")
    axis.legend(fontsize=8, frameon=False, ncol=2)
    axis.grid(axis="y", color="#e6ecec", linewidth=0.8)
    axis.set_axisbelow(True)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--blocks", required=True)
    parser.add_argument("--runs", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--dpi", type=int, default=150)
    args = parser.parse_args()
    os.makedirs(args.out, exist_ok=True)

    blocks = pd.read_csv(args.blocks, sep="\t")
    runs = pd.read_csv(args.runs, sep="\t")
    coverages = sorted(blocks.coverage.unique())
    block_ids = sorted(blocks.block.unique(), key=lambda b: (b.split(":")[0], int(b.split(":")[1])))
    role = blocks.role.iloc[0] if len(blocks) else "unknown"

    figure, axes = plt.subplots(3, 1, figsize=(max(8, 0.8 * len(block_ids)), 10), constrained_layout=True,
                                gridspec_kw={"height_ratios": [3, 3, 2]})
    for axis, arm in zip(axes[:2], ARMS):
        violins(axis, blocks[blocks.arm == arm], block_ids, coverages, COLOUR[arm])
        axis.set_title("%s arm, %s donors: %s" % (arm, role, ", ".join("%gx" % c for c in coverages)), fontsize=10)
    rates(axes[2], runs, coverages)
    figure.suptitle("C4 leave-zero-out by nested coverage", fontsize=12)
    path = os.path.join(args.out, "coverage_blocks_%s.png" % role)
    figure.savefig(path, dpi=args.dpi)
    print("wrote %s" % path)

    # A compact table beside the figure: nothing in the plot that is not also readable as numbers.
    summary = []
    for arm in ARMS:
        for coverage in coverages:
            mine = runs[(runs.arm == arm) & (runs.coverage == coverage)]
            called = mine[mine.run_no_call == "false"]
            distances = pd.to_numeric(called.unordered_distance, errors="coerce")
            rows = blocks[(blocks.arm == arm) & (blocks.coverage == coverage)]
            exact_blocks = (rows.exact == "true").sum()
            summary.append({"arm": arm, "coverage": coverage, "runs": len(mine),
                            "no_calls": int((mine.run_no_call == "true").sum()),
                            "donors_fully_exact": int((distances == 0).sum()),
                            "blocks_exact": int(exact_blocks), "blocks": int(len(rows)),
                            "distance_total": int(distances.fillna(0).sum())})
    table = pd.DataFrame(summary)
    table.to_csv(os.path.join(args.out, "coverage_summary_%s.tsv" % role), sep="\t", index=False)
    print(table.to_string(index=False))


if __name__ == "__main__":
    main()
