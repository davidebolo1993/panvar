#!/usr/bin/env python3
"""Slide figure of block calls coloured by a categorical outcome (plotting only): vertically stacked panels, one per
row of --rows, ordered blocks on x, excess over the availability floor on the symlog bp axis with the same limits,
coverage dodge, marker shapes and deterministic jitter as the conditional-GQ figure; every point is drawn and
coloured by the presentation category its dosage status (and repair flag) maps to.

Inputs (TSV):
  --points      run_id, row, coverage, block_index, block_id, excess and the --category-column (default dosage_status);
                dosage_repair_applied when a rule conditions on it: one row per plotted observation (excess present on
                every row); other columns are carried into the outputs
  --rules       category, label, values (comma list of --category-column values; the column may also be named
                dosage_status), repair_applied (any | true | false), color, size, edge (a colour or none): in drawing
                order; every point must match exactly one rule
  --category-column  the points column the rules classify (default dosage_status)
  --rows        row, label
  --block-lengths block_id, displayed_length_range
  --unplotted   (optional) observations of the same layer that cannot be plotted (no defined error value), with
                run_id, block_id, row and a reason column; counted in the outputs (and the footer, unless
                --unplotted-note is empty), never drawn
  --panel-axes  (optional) row, ymin, ymax: fixed y-limits per panel (e.g. a reference figure's panel_axes.tsv), so a
                subset of points keeps that figure's axes; a point outside them is refused
  --coverages   (optional) comma list: the coverage dodge order, so a subset keeps the reference positions
  --blocks-from-lengths  take the x-axis blocks, in order, from --block-lengths (a subset may miss some blocks); the
                points' block_index must equal the position there
  --footer      (optional) a footer line; joined with the unplotted note when both are given

Writes <out-dir>/figures/<stem>.png|pdf (16:9 slide), plotted_points.tsv, category_counts.tsv (row, coverage, block,
category), panel_axes.tsv, layout_check.tsv, unplotted.tsv (when given) and captions.md. Exit 0, or 2 on a refusal.
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import pandas as pd  # noqa: E402

from genotype_report import plots as P  # noqa: E402
from genotype_report.common import Refusal, coverage_order, read_tsv, require_columns, to_number, write_tsv  # noqa: E402

REQUIRED = ["run_id", "row", "coverage", "block_index", "block_id", "excess"]
RULE_COLUMNS = ["category", "label", "repair_applied", "color", "size", "edge"]


def parse(argv):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    for a in ("--points", "--rules", "--rows", "--block-lengths", "--out-dir"):
        p.add_argument(a, required=True)
    p.add_argument("--unplotted")
    p.add_argument("--unplotted-note", default="{n} observation(s) have no defined error value and are not plotted.",
                   help="footer sentence for --unplotted; {n} is replaced by their count")
    p.add_argument("--stem", default="block_categories_slide")
    p.add_argument("--title", default="Block calls by outcome")
    p.add_argument("--xlabel", default="block (genomic order)")
    p.add_argument("--caption-note", default="")
    p.add_argument("--panel-axes")
    p.add_argument("--coverages")
    p.add_argument("--blocks-from-lengths", action="store_true")
    p.add_argument("--footer", default="")
    p.add_argument("--category-column", default="dosage_status")
    return p.parse_args(argv)


def categorize(pts, rules, column="dosage_status"):
    """The category of every point: the one rule whose values contain the point's `column` value and whose repair
    condition holds. A point matching no rule, or several, is refused."""
    hits = []
    vcol = "values" if "values" in rules.columns else "dosage_status"
    needs_flag = (rules.repair_applied != "any").any()
    if needs_flag and "dosage_repair_applied" not in pts.columns:
        raise Refusal("a rule conditions on dosage_repair_applied, which the points lack")
    applied = (pts.dosage_repair_applied.astype(str).str.lower() == "true") if needs_flag else pd.Series(False, index=pts.index)
    for r in rules.itertuples():
        statuses = {s.strip() for s in getattr(r, vcol).split(",")}
        cond = {"any": pd.Series(True, index=pts.index), "true": applied, "false": ~applied}.get(r.repair_applied)
        if cond is None:
            raise Refusal(f"rule {r.category}: repair_applied must be any, true or false")
        hits.append(pts[column].isin(statuses) & cond)
    m = pd.concat(hits, axis=1)
    m.columns = list(rules.category)
    n = m.sum(axis=1)
    if (n == 0).any():
        raise Refusal(f"no rule covers {column} value(s) {', '.join(sorted(set(pts[column][n == 0])))}")
    if (n > 1).any():
        raise Refusal("a point matches more than one category rule")
    return m.idxmax(axis=1)


def run(argv):
    a = parse(argv)
    if os.path.exists(a.out_dir) and os.listdir(a.out_dir):
        raise Refusal(f"--out-dir {a.out_dir} exists and is not empty")
    pts = read_tsv(a.points, "points")
    require_columns(pts, REQUIRED + [a.category_column], "points")
    if pts.duplicated(["run_id", "block_id"]).any():
        raise Refusal("a (run, block) observation appears twice")
    for c in ("coverage", "block_index"):
        pts[c] = pts[c].astype(int)
    pts["excess"] = to_number(pts.excess)
    if pts.excess.isna().any() or (pts.excess < 0).any():
        raise Refusal("every plotted observation needs a non-negative excess")
    rules = read_tsv(a.rules, "category rules")
    require_columns(rules, RULE_COLUMNS, "category rules")
    if "values" not in rules.columns and "dosage_status" not in rules.columns:
        raise Refusal("category rules need a values (or dosage_status) column")
    pts["category"] = categorize(pts, rules, a.category_column)
    rtab = read_tsv(a.rows, "rows")
    rows = list(zip(rtab.row, rtab.label))
    if set(pts.row) - set(rtab.row):
        raise Refusal("points in rows not listed in --rows")
    lengths = read_tsv(a.block_lengths, "block lengths")
    if a.blocks_from_lengths:
        order = list(lengths.block_id)
        pos = {b: i for i, b in enumerate(order)}
        if any(b not in pos or pos[b] != i for b, i in zip(pts.block_id, pts.block_index)):
            raise Refusal("a point's block is missing from, or at another position in, the block-length table")
    else:
        order = list(pts.sort_values("block_index").drop_duplicates("block_id").block_id)
    ticks = P.block_labels(order, lengths)
    covs = coverage_order(pts.coverage)
    if a.coverages:
        covs = [int(c) for c in a.coverages.split(",")]
        if set(pts.coverage) - set(covs):
            raise Refusal("points at a coverage not listed in --coverages")
    ylims = None
    if a.panel_axes:
        ax = read_tsv(a.panel_axes, "panel axes")
        ylims = {r.row: (float(r.ymin), float(r.ymax)) for r in ax.itertuples()}
        if set(rtab.row) - set(ylims):
            raise Refusal("--panel-axes lacks a row of --rows")
    unplotted = None
    footer = a.footer or None
    if a.unplotted:
        unplotted = read_tsv(a.unplotted, "unplotted observations")
        require_columns(unplotted, ["run_id", "block_id", "row", "reason"], "unplotted observations")
        clash = set(zip(unplotted.run_id, unplotted.block_id)) & set(zip(pts.run_id, pts.block_id))
        if clash:
            raise Refusal(f"{len(clash)} unplotted observation(s) are also in --points")
        note = a.unplotted_note.replace("{n}", str(len(unplotted)))
        footer = " ".join(x for x in (a.footer, note) if x) or None
    cats = [(r.category, r.label, r.color, float(r.size), None if r.edge == "none" else r.edge) for r in rules.itertuples()]
    try:
        fig, drawn = P.fig_block_categories(pts, order, rows, covs, ticks, a.title, cats, xlabel=a.xlabel, footer=footer, ylims=ylims)
    except ValueError as e:
        raise Refusal(str(e))
    clipped, overlapping = P.text_layout_problems(fig)
    os.makedirs(os.path.join(a.out_dir, "figures"), exist_ok=True)
    axes_tab = pd.DataFrame(fig.panel_axes, columns=["row", "ymin", "ymax", "yscale"])
    size = (fig.get_figwidth(), fig.get_figheight())
    P.save(fig, os.path.join(a.out_dir, "figures", a.stem), dpi=P.SLIDE["dpi"])
    d = pd.DataFrame(drawn, columns=["row", "run_id", "block_id", "coverage", "x", "excess", "category"])
    extra = [c for c in pts.columns if c not in ("row", "coverage", "excess", "category")]
    d = d.merge(pts[["run_id", "block_id"] + [c for c in extra if c not in ("run_id", "block_id")]], on=["run_id", "block_id"],
                how="left", validate="one_to_one")
    if len(d) != len(pts):
        raise Refusal(f"drew {len(d)} points for {len(pts)} observations")
    write_tsv(d, os.path.join(a.out_dir, "plotted_points.tsv"))
    counts = d.groupby(["row", "coverage", "block_index", "block_id", "category"], sort=True).size().rename("points").reset_index()
    write_tsv(counts, os.path.join(a.out_dir, "category_counts.tsv"))
    write_tsv(axes_tab, os.path.join(a.out_dir, "panel_axes.tsv"))
    write_tsv(pd.DataFrame([{"figure": a.stem, "width_in": size[0], "height_in": size[1], "clipped_texts": len(clipped),
                             "overlapping_text_pairs": len(overlapping), "detail": "; ".join(clipped + overlapping) or "none"}]),
              os.path.join(a.out_dir, "layout_check.tsv"))
    if unplotted is not None:
        write_tsv(unplotted, os.path.join(a.out_dir, "unplotted.tsv"))
    tot = d.category.value_counts()
    lines = [f"# Caption\n", f"### {a.title}  (`{a.stem}`)", "", a.caption_note.strip(), "",
             "Categories (points): " + "; ".join(f"{r.label}: {int(tot.get(r.category, 0)):,}" for r in rules.itertuples()) + ".", "",
             "Per panel: " + "; ".join(f"{lab}: {int((d.row == rv).sum()):,}" for rv, lab in rows) + "."]
    if footer:
        lines += ["", footer]
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
