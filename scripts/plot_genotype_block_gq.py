#!/usr/bin/env python3
"""Conditional-GQ diagnostic for `panvar genotype` evaluations (plotting only): one figure of vertically stacked panels,
one per row of --rows, ordered blocks on x, excess over the availability floor on the same symlog bp axis as the
ordered-block raincloud figures, coverage dodged and encoded by marker shape, every point coloured continuously by
conditional GQ through log10(1 + GQ), with the colour bar labelled in GQ values.

No point is removed or capped for its GQ: a missing GQ is drawn in neutral grey and an infinite one in the top colour
with a black outline; both are counted in the outputs.

Inputs (TSV):
  --points        run_id, row, coverage, block_index, block_id, excess, conditional_GQ (one row per plotted observation:
                  only evaluated, truth-defined observations; excess must be present); optional experiment, subject,
                  layer columns are carried into the point table
  --block-lengths block_id, displayed_length_range (shown under each block ID)
  --rows          row, label (panel order and titles)
  --caption-note  caption text describing the inputs

Writes <out-dir>/figures/<stem>.png|pdf, plotted_points_gq.tsv, gq_counts.tsv, panel_axes.tsv and captions.md. An --out-dir that
exists and is not empty is refused. Exit 0 on success, 2 on a refusal.
"""
import argparse
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import pandas as pd  # noqa: E402

from genotype_report import plots as P  # noqa: E402
from genotype_report.common import Refusal, coverage_order, read_tsv, require_columns, to_number, write_tsv  # noqa: E402

REQUIRED = ["run_id", "row", "coverage", "block_index", "block_id", "excess", "conditional_GQ"]


def parse(argv):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    for a in ("--points", "--block-lengths", "--rows", "--out-dir"):
        p.add_argument(a, required=True)
    p.add_argument("--stem", default="conditional_gq_block_excess")
    p.add_argument("--title", default="Conditional genotype quality by block")
    p.add_argument("--ylabel", default="Excess over the availability floor (bp)")
    p.add_argument("--xlabel", default="block (genomic order)")
    p.add_argument("--colorbar-label", default="Conditional GQ")
    p.add_argument("--caption-note", default="")
    p.add_argument("--slide", action="store_true", help="16:9 presentation layout (13.333 x 7.5 in), slide-size type")
    p.add_argument("--footer", default="", help="footer line of the --slide figure")
    return p.parse_args(argv)


def gq_number(s):
    """Conditional GQ from its text: NA/empty -> NaN (missing), inf/+inf/Infinity -> +inf; anything else must parse."""
    out = []
    for v in s:
        t = str(v).strip()
        if t in ("", "NA", "nan", "NaN"):
            out.append(math.nan)
            continue
        try:
            x = float(t)
        except ValueError:
            raise Refusal(f"conditional_GQ is not a number: {t}")
        if x < 0 or x == -math.inf:
            raise Refusal(f"negative conditional_GQ: {t}")
        out.append(x)
    return pd.Series(out, index=s.index, dtype=float)


def caption(a, counts, rows):
    lines = [f"### {a.title}  (`{a.stem}`)", ""]
    if a.caption_note.strip():
        lines += [a.caption_note.strip(), ""]
    lines += [
        "Each point is one truth-defined, evaluated block call. Vertical position: excess over the availability floor "
        "(bp; symlog axis, linear from 0 to 1 bp and logarithmic above). Colour: conditional GQ, mapped through "
        "log10(1 + GQ) and labelled in GQ values on the colour bar. Marker shape and horizontal position within each "
        "block encode coverage; the horizontal offset within a coverage is a deterministic jitter (a hash of run and "
        "block) carrying no information.", "",
        "Conditional GQ measures how far the called genotype is separated from the alternative genotypes among the "
        "candidates actually evaluated. It is not a calibrated probability that the call is correct in the real "
        "world, and a high GQ cannot detect a missing truth allele, an incomplete candidate space or any other model "
        "failure. No GQ threshold or filter is applied: every observation is drawn.", "",
        "Missing GQ would be drawn in neutral grey and an infinite GQ in the top colour with a black outline; in this "
        f"figure there are {int(counts.missing.sum())} missing and {int(counts.infinite.sum())} infinite values.", ""]
    per = ", ".join(f"{lab}: {int(counts[counts.row == r].points.sum())}" for r, lab in rows)
    lines += [f"Points per panel: {per}.", ""]
    return "\n".join(lines)


def run(argv):
    a = parse(argv)
    if os.path.exists(a.out_dir) and os.listdir(a.out_dir):
        raise Refusal(f"--out-dir {a.out_dir} exists and is not empty")
    pts = read_tsv(a.points, "GQ points")
    require_columns(pts, REQUIRED, "GQ points")
    if pts.duplicated(["run_id", "block_id"]).any():
        raise Refusal("a (run, block) observation appears twice")
    for c in ("coverage", "block_index"):
        pts[c] = pts[c].astype(int)
    pts["excess"] = to_number(pts.excess)
    if pts.excess.isna().any():
        raise Refusal("a plotted observation has no excess (only evaluated, truth-defined observations belong here)")
    if (pts.excess < 0).any():
        raise Refusal("negative excess")
    pts["conditional_GQ_text"] = pts.conditional_GQ
    pts["conditional_GQ"] = gq_number(pts.conditional_GQ)
    rtab = read_tsv(a.rows, "rows")
    rows = list(zip(rtab.row, rtab.label))
    unknown = sorted(set(pts.row) - set(rtab.row))
    if unknown:
        raise Refusal(f"points in rows not listed in --rows: {', '.join(unknown)}")
    lengths = read_tsv(a.block_lengths, "block lengths")
    order = list(pts.sort_values("block_index").drop_duplicates("block_id").block_id)
    ticks = P.block_labels(order, lengths)
    covs = coverage_order(pts.coverage)
    fig, drawn = P.fig_block_gq(pts, order, rows, covs, ticks, a.title, ylabel=a.ylabel, xlabel=a.xlabel,
                                cbar_label=a.colorbar_label, slide=a.slide, footer=a.footer or None)
    os.makedirs(os.path.join(a.out_dir, "figures"), exist_ok=True)
    axes_tab = pd.DataFrame(fig.panel_axes, columns=["row", "ymin", "ymax", "yscale"])
    if a.slide:
        clipped, overlapping = P.text_layout_problems(fig)
        write_tsv(pd.DataFrame([{"figure": a.stem, "width_in": fig.get_figwidth(), "height_in": fig.get_figheight(),
                                 "clipped_texts": len(clipped), "overlapping_text_pairs": len(overlapping),
                                 "detail": "; ".join(clipped + overlapping) or "none"}]),
                  os.path.join(a.out_dir, "layout_check.tsv"))
    P.save(fig, os.path.join(a.out_dir, "figures", a.stem), dpi=P.SLIDE["dpi"] if a.slide else None)
    write_tsv(axes_tab, os.path.join(a.out_dir, "panel_axes.tsv"))
    d = pd.DataFrame(drawn, columns=["row", "run_id", "block_id", "coverage", "x", "excess", "conditional_GQ", "color_value_log10_1p_gq",
                                     "color_hex", "gq_class"])
    extra = [c for c in pts.columns if c not in REQUIRED + ["conditional_GQ_text"]]
    d = d.merge(pts[["run_id", "block_id", "conditional_GQ_text"] + extra], on=["run_id", "block_id"], how="left", validate="one_to_one")
    if len(d) != len(pts):
        raise Refusal(f"drew {len(d)} points for {len(pts)} observations")
    write_tsv(d, os.path.join(a.out_dir, "plotted_points_gq.tsv"))
    counts = (d.assign(missing=d.gq_class == "missing", infinite=d.gq_class == "+inf")
              .groupby(["row", "coverage", "block_id"], sort=True)
              .agg(points=("run_id", "size"), missing=("missing", "sum"), infinite=("infinite", "sum")).reset_index())
    write_tsv(counts, os.path.join(a.out_dir, "gq_counts.tsv"))
    with open(os.path.join(a.out_dir, "captions.md"), "w") as f:
        f.write("# Caption\n\n" + caption(a, counts, rows))
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
