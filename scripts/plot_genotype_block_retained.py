#!/usr/bin/env python3
"""Retained-call raincloud slide for `panvar genotype` evaluations (plotting only): the block calls passing a
truth-free retention rule (--retain-column == --retain-value), one symlog-bp raincloud panel per row of --rows on fixed
reference y-limits, half violins and quartiles recomputed from the retained calls only, points at their deterministic
raincloud positions. Excluded calls are not drawn; they are listed with their fields.

Inputs (TSV):
  --points      run_id, row, coverage, block_index, block_id, excess and the --retain-column: one row per truth-defined,
                evaluated observation (the whole population, retained and excluded); other columns are carried
  --rows        row, label
  --block-lengths block_id, displayed_length_range (row order = block axis; block_index must match)
  --panel-axes  row, ymin, ymax: the fixed y-limits per panel (a point outside them is refused)

Truth (excess) is never used to decide membership: the rule reads only --retain-column.

Writes <out-dir>/figures/<stem>.png|pdf (16:9 slide), retained_points.tsv, retained_violins.tsv, retention_counts.tsv
(row, coverage, block: total, retained, excluded), excluded_calls.tsv, panel_axes.tsv, layout_check.tsv, captions.md.
Exit 0, or 2 on a refusal.
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import pandas as pd  # noqa: E402

from genotype_report import plots as P  # noqa: E402
from genotype_report.common import Refusal, read_tsv, require_columns, to_number, write_tsv  # noqa: E402

REQUIRED = ["run_id", "row", "coverage", "block_index", "block_id", "excess"]


def parse(argv):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    for a in ("--points", "--rows", "--block-lengths", "--panel-axes", "--retain-column", "--retain-value", "--out-dir"):
        p.add_argument(a, required=True)
    p.add_argument("--coverages", required=True, help="comma list: the coverage dodge order")
    p.add_argument("--stem", default="retained_block_excess_slide")
    p.add_argument("--title", default="Block excess among retained calls")
    p.add_argument("--xlabel", default="block (genomic order)")
    p.add_argument("--footer", default="")
    p.add_argument("--caption-note", default="")
    return p.parse_args(argv)


def run(argv):
    a = parse(argv)
    if os.path.exists(a.out_dir) and os.listdir(a.out_dir):
        raise Refusal(f"--out-dir {a.out_dir} exists and is not empty")
    pts = read_tsv(a.points, "points")
    require_columns(pts, REQUIRED + [a.retain_column], "points")
    if pts.duplicated(["run_id", "block_id"]).any():
        raise Refusal("a (run, block) observation appears twice")
    for c in ("coverage", "block_index"):
        pts[c] = pts[c].astype(int)
    pts["excess"] = to_number(pts.excess)
    if pts.excess.isna().any() or (pts.excess < 0).any():
        raise Refusal("every observation needs a non-negative excess (truth-defined, evaluated calls only)")
    retained_mask = pts[a.retain_column] == a.retain_value
    rtab = read_tsv(a.rows, "rows")
    rows = list(zip(rtab.row, rtab.label))
    if set(pts.row) - set(rtab.row):
        raise Refusal("points in rows not listed in --rows")
    lengths = read_tsv(a.block_lengths, "block lengths")
    order = list(lengths.block_id)
    pos = {b: i for i, b in enumerate(order)}
    if any(b not in pos or pos[b] != i for b, i in zip(pts.block_id, pts.block_index)):
        raise Refusal("a point's block is missing from, or at another position in, the block-length table")
    covs = [int(c) for c in a.coverages.split(",")]
    if set(pts.coverage) - set(covs):
        raise Refusal("points at a coverage not listed in --coverages")
    axes = read_tsv(a.panel_axes, "panel axes")
    ylims = {r.row: (float(r.ymin), float(r.ymax)) for r in axes.itertuples()}
    if set(rtab.row) - set(ylims):
        raise Refusal("--panel-axes lacks a row of --rows")
    kept = pts[retained_mask]
    try:
        fig, drawn = P.fig_block_raincloud_rows(kept, order, rows, covs, P.block_labels(order, lengths), a.title, ylims,
                                                xlabel=a.xlabel, footer=a.footer or None)
    except ValueError as e:
        raise Refusal(str(e))
    clipped, overlapping = P.text_layout_problems(fig)
    size = (fig.get_figwidth(), fig.get_figheight())
    axes_tab = pd.DataFrame(fig.panel_axes, columns=["row", "ymin", "ymax", "yscale"])
    os.makedirs(os.path.join(a.out_dir, "figures"), exist_ok=True)
    P.save(fig, os.path.join(a.out_dir, "figures", a.stem), dpi=P.SLIDE["dpi"])
    d = pd.DataFrame(drawn["points"], columns=["row", "metric", "run_id", "block_id", "coverage", "x", "excess"]).drop(columns="metric")
    carry = [c for c in pts.columns if c not in ("row", "coverage", "excess")]
    d = d.merge(pts[carry], on=["run_id", "block_id"], how="left", validate="one_to_one")
    if len(d) != len(kept):
        raise Refusal(f"drew {len(d)} points for {len(kept)} retained observations")
    write_tsv(d, os.path.join(a.out_dir, "retained_points.tsv"))
    write_tsv(pd.DataFrame(drawn["violins"], columns=["row", "metric", "block_id", "coverage", "kind", "n", "q1", "median", "q3"]),
              os.path.join(a.out_dir, "retained_violins.tsv"))
    t = pts.assign(retained=retained_mask)
    counts = (t.groupby(["row", "coverage", "block_index", "block_id"], sort=True)
              .agg(total=("run_id", "size"), retained=("retained", "sum")).reset_index())
    counts["excluded"] = counts.total - counts.retained
    write_tsv(counts, os.path.join(a.out_dir, "retention_counts.tsv"))
    write_tsv(pts[~retained_mask].sort_values(["row", "block_index", "coverage", "run_id"], kind="stable"),
              os.path.join(a.out_dir, "excluded_calls.tsv"))
    write_tsv(axes_tab, os.path.join(a.out_dir, "panel_axes.tsv"))
    write_tsv(pd.DataFrame([{"figure": a.stem, "width_in": size[0], "height_in": size[1], "clipped_texts": len(clipped),
                             "overlapping_text_pairs": len(overlapping), "detail": "; ".join(clipped + overlapping) or "none"}]),
              os.path.join(a.out_dir, "layout_check.tsv"))
    lines = [f"# Caption\n", f"### {a.title}  (`{a.stem}`)", "", f"Retention rule: {a.retain_column} == {a.retain_value}.", "",
             a.caption_note.strip(), "",
             "Per panel (retained / all): " + "; ".join(f"{lab}: {int((kept.row == rv).sum()):,} / {int((pts.row == rv).sum()):,}"
                                                         for rv, lab in rows) + "."]
    if a.footer:
        lines += ["", a.footer]
    with open(os.path.join(a.out_dir, "captions.md"), "w") as f:
        f.write("\n".join(lines) + "\n")
    return a.out_dir


def main(argv=None):
    try:
        out = run(sys.argv[1:] if argv is None else argv)
    except Refusal as e:
        print(f"REFUSED: {e}", file=sys.stderr)
        return 2
    print(f"figure -> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
