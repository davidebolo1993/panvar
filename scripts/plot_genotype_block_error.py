#!/usr/bin/env python3
"""Ordered-block raincloud figures of block-call error for `panvar genotype` evaluations (plotting only).

Each figure is one row of the --rows table: two vertically aligned panels sharing the ordered block axis, excess over
the availability floor (bp, symlog: linear in [0, 1] bp, logarithmic above) and normalised excess (linear). In every
block and coverage, a half violin and the same observations as deterministic jittered points, with Q1, median and Q3
marks, under an n=(C,Z,U) label. Also writes the counts behind every label, every drawn point, the locus companion table
and the captions. Nothing is filtered, capped or recomputed from sequences; the values are the supplied ones.

Inputs (TSV):
  --points        run_id, row, coverage, block_index, block_id, truth_defined, availability_floor, excess,
                  normalized_excess: one row per called block observation (no-call runs contribute none); the metrics
                  are NA exactly where truth is undefined
  --block-lengths block_id, displayed_length_range (shown under each block ID)
  --rows          row, title, stem, locus_mode (sequence_correct | floor_attaining)
  --run-counts    row, coverage, called_runs, no_call_runs
  --expected-called (optional) row, coverage, called_runs: refuse unless the points match
  --phase-calls   (optional) run_id, exact_phased_locus
  --phase-projected (optional) run_id, projected_status, projected_concordant
  --layer-note, --length-note, --normalization-note, --xlabel: caption and axis text describing the inputs

Writes <out-dir>/figures/<stem>.png|pdf, block_counts.tsv, plotted_points.tsv, violins.tsv, <locus-table> and
captions.md. An --out-dir that exists and is not empty is refused. Exit 0 on success, 2 on a refusal.
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import pandas as pd  # noqa: E402

from genotype_report import aggregate as A  # noqa: E402
from genotype_report import plots as P  # noqa: E402
from genotype_report.common import Refusal, coverage_order, read_tsv, require_columns, to_bool, to_number, write_tsv  # noqa: E402

EXCESS_LABEL = "Excess over the availability floor (bp)"
NORMALIZED_LABEL = "Diploid-length-normalized excess"   # excess / max(truth1+truth2, called1+called2) bp
XLABEL = "block (genomic order)"


def parse(argv):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    for a in ("--points", "--block-lengths", "--rows", "--run-counts", "--out-dir"):
        p.add_argument(a, required=True)
    for a in ("--expected-called", "--phase-calls", "--phase-projected"):
        p.add_argument(a)
    p.add_argument("--locus-table", default="locus_accuracy.tsv")
    p.add_argument("--layer-note", default="")
    p.add_argument("--length-note", default="Below each block ID: the range shown in the supplied block-length table.")
    p.add_argument("--normalization-note", default="")
    p.add_argument("--excess-label", default=EXCESS_LABEL)
    p.add_argument("--normalized-label", default=NORMALIZED_LABEL)
    p.add_argument("--xlabel", default=XLABEL, help="x-axis label; say what the range under each block ID is")
    p.add_argument("--slide", action="store_true",
                   help="16:9 presentation layout (13.333 x 7.5 in): no count labels on the axes (denominators stay in "
                        "block_counts.tsv and the caption), slide-size type, the stems take --stem-suffix")
    p.add_argument("--stem-suffix", default="_slide", help="appended to each stem with --slide")
    p.add_argument("--footer", default="Not plotted: truth-undefined observations and no-call runs. "
                                       "Denominators n=(C,Z,U) per block and coverage: block_counts.tsv.",
                   help="footer line of the --slide figures")
    return p.parse_args(argv)


def load(a):
    pts = read_tsv(a.points, "block points")
    require_columns(pts, A.POINT_COLUMNS, "block points")
    for c in ("coverage", "block_index"):
        pts[c] = pts[c].astype(int)
    td = to_bool(pts.truth_defined)
    if td.isna().any():
        raise Refusal("block points: truth_defined must be true or false on every row")
    pts["truth_defined"] = td.astype(bool)
    for c in ("availability_floor", "excess", "normalized_excess"):
        pts[c] = to_number(pts[c])
    rc = read_tsv(a.run_counts, "run counts")
    for c in ("coverage", "called_runs", "no_call_runs"):
        rc[c] = rc[c].astype(int)
    ex = None
    if a.expected_called:
        ex = read_tsv(a.expected_called, "expected called runs")
        for c in ("coverage", "called_runs"):
            ex[c] = ex[c].astype(int)
    return pts, rc, ex


def caption(row, counts, rc, a):
    r = rc[rc.row == row["row"]].sort_values("coverage")
    nocalls = ", ".join(f"{c}x: {n}" for c, n in zip(r.coverage, r.no_call_runs))
    called = ", ".join(f"{c}x: {n}" for c, n in zip(r.coverage, r.called_runs))
    k = counts[counts.row == row["row"]]
    return "\n".join([
        f"### {row['title']}  (`{row['stem']}`)", "",
        f"{a.layer_note}".strip(), "",
        f"Top: {a.excess_label.lower()}; symlog axis, linear from 0 to 1 bp and logarithmic above. Bottom: "
        f"{a.normalized_label.lower()}, linear axis. {a.normalization_note}".strip(), "",
        "Each block and coverage shows one half violin and the same observations as points jittered on its open side "
        "(deterministic offsets from a hash of run and block). Short black marks across the violin are Q1 and Q3; the "
        "thicker mark is the median (linear interpolation between order statistics). All violins share one maximum "
        f"width; sample size is given only by {'the companion table' if a.slide else 'the label'}. Densities are Gaussian kernel estimates (Scott's bandwidth) in "
        "the displayed axis coordinate, evaluated only over the observed range, so they never extend below zero or beyond "
        "the largest observation; a distribution with zero variance is drawn as a flat bar at its single value. No "
        "observation is capped: each axis reaches its largest value.", "",
        ("Denominators: no count labels are drawn on this slide; block_counts.tsv gives n=(C,Z,U) for every block and "
         "coverage. " if a.slide else "Label n=(C,Z,U) above each block and coverage: ") +
        "C = called block observations for that group, coverage and block, "
        "including observations whose benchmark truth is undefined; Z = truth-defined observations with excess exactly 0; "
        "U = observations whose benchmark truth is undefined, for which no distance is plotted. Plotted points per "
        "violin = C − U. U counts undefined truth, not missing genotype calls.", "",
        f"Called runs (each contributes one observation per block): {called}. No-call runs contribute no block "
        f"observation and are excluded from C: {nocalls}. "
        f"Total U in this figure: {int(k.U.sum())} observation(s)." + (exceptions(k) if a.slide else ""), "",
        a.length_note.strip(), ""])


def exceptions(k):
    """The cells whose truth is undefined, named once (slide captions carry no per-cell labels)."""
    u = k[k.U > 0]
    if not len(u):
        return " No block has a truth-undefined observation."
    parts = []
    for b, g in u.groupby("block_id", sort=False):
        parts.append(f"{b} (" + ", ".join(f"{c}x: {n}" for c, n in zip(g.coverage, g.U)) + ")")
    return " Truth-undefined observations, not plotted: " + "; ".join(parts) + "."


def run(argv):
    a = parse(argv)
    if os.path.exists(a.out_dir) and os.listdir(a.out_dir):
        raise Refusal(f"--out-dir {a.out_dir} exists and is not empty")
    pts, rc, ex = load(a)
    A.validate_block_points(pts, rc, ex)
    rows = read_tsv(a.rows, "rows")
    lengths = read_tsv(a.block_lengths, "block lengths")
    order = list(pts.sort_values("block_index").drop_duplicates("block_id").block_id)
    ticks = P.block_labels(order, lengths)
    counts = A.block_call_counts(pts)
    os.makedirs(os.path.join(a.out_dir, "figures"), exist_ok=True)
    drawn_pts, drawn_v, caps, layout = [], [], [f"# Captions\n"], []
    for _, row in rows.iterrows():
        d = pts[pts.row == row["row"]]
        covs = coverage_order(d.coverage)
        cnt = {(r.block_id, r.coverage): (r.C, r.Z, r.U) for r in counts[counts.row == row["row"]].itertuples()}
        plotted = d[d.truth_defined]
        if a.slide:
            row = row.copy()
            row["stem"] = row["stem"] + a.stem_suffix
        fig, drawn = P.fig_block_raincloud(plotted, order, covs,
                                           [("excess", a.excess_label, "symlog_bp"), ("normalized_excess", a.normalized_label, "linear")],
                                           cnt, row["title"], ticks, xlabel=a.xlabel, slide=a.slide,
                                           footer=a.footer if a.slide else None)
        if a.slide:
            clipped, overlapping = P.text_layout_problems(fig)
            layout.append({"figure": row["stem"], "width_in": fig.get_figwidth(), "height_in": fig.get_figheight(),
                           "clipped_texts": len(clipped), "overlapping_text_pairs": len(overlapping),
                           "detail": "; ".join(clipped + overlapping) or "none"})
        P.save(fig, os.path.join(a.out_dir, "figures", row["stem"]), dpi=P.SLIDE["dpi"] if a.slide else None)
        drawn_pts += [(row["row"],) + p for p in drawn["points"]]
        drawn_v += [(row["row"],) + v for v in drawn["violins"]]
        caps.append(caption(row, counts, rc, a))
    write_tsv(counts, os.path.join(a.out_dir, "block_counts.tsv"))
    if a.slide:
        write_tsv(pd.DataFrame(layout), os.path.join(a.out_dir, "layout_check.tsv"))
    write_tsv(pd.DataFrame(drawn_pts, columns=["row", "metric", "run_id", "block_id", "coverage", "x", "y"]),
              os.path.join(a.out_dir, "plotted_points.tsv"))
    write_tsv(pd.DataFrame(drawn_v, columns=["row", "metric", "block_id", "coverage", "kind", "n", "q1", "median", "q3"]),
              os.path.join(a.out_dir, "violins.tsv"))
    pc = pp = None
    if a.phase_calls:
        pc = read_tsv(a.phase_calls, "phase calls")
        pc["exact_phased_locus"] = to_bool(pc.exact_phased_locus)
    if a.phase_projected:
        pp = read_tsv(a.phase_projected, "projected phase")
        pp["projected_concordant"] = to_bool(pp.projected_concordant).fillna(False)
    modes = dict(zip(rows.row, rows.locus_mode))
    write_tsv(A.locus_accuracy(pts, rc, modes, pc, pp), os.path.join(a.out_dir, a.locus_table))
    with open(os.path.join(a.out_dir, "captions.md"), "w") as f:
        f.write("\n".join(caps))
    return a.out_dir


def main(argv=None):
    try:
        out = run(sys.argv[1:] if argv is None else argv)
    except Refusal as e:
        print(f"REFUSED: {e}", file=sys.stderr)
        return 2
    print(f"figures -> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
