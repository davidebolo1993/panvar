#!/usr/bin/env python3
"""Slide figures of a synthetic mosaic benchmark design, from its switch map (plotting only; nothing is simulated).

Two 16:9 figures (13.333 x 7.5 in; PNG 3200 x 1800 px, vector PDF):
  <design-stem>  one example plan per class: H1 and H2 tracks over the ordered blocks, coloured by source A or B
                 (both homologues start on A and toggle at every switch), switches drawn on block boundaries;
  <map-stem>     every plan's switch boundaries, H1-only / H2-only / both distinguished.

Inputs (TSV):
  --switch-map  switch_class, plan, source_a, source_b, boundary_index, boundary ('left>right'), hap1_switch, hap2_switch
  --classes     class, label, example_plan (row order of the design figure; the example plan is labelled as such)
  --note        a line set under the design figure (repeatable)
  --caption-note caption text describing the experiment

Writes <out-dir>/figures/<stem>.png|pdf, mosaic_design_tracks.tsv (every plan, homologue and block; example rows
marked), switch_plans.tsv (donors, switch positions and counts per plan), layout_check.tsv and captions.md. Refuses
a switch map whose plans do not share one chained block order. Exit 0 on success, 2 on a refusal.
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import pandas as pd  # noqa: E402

from genotype_report import aggregate as A  # noqa: E402
from genotype_report import plots as P  # noqa: E402
from genotype_report.common import Refusal, read_tsv, write_tsv  # noqa: E402


def parse(argv):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    for a in ("--switch-map", "--classes", "--out-dir"):
        p.add_argument(a, required=True)
    p.add_argument("--title", default="Synthetic mosaic benchmark design")
    p.add_argument("--map-title", default="Mosaic switch plans")
    p.add_argument("--design-stem", default="mosaic_benchmark_design_slide")
    p.add_argument("--map-stem", default="mosaic_switch_map_full_slide")
    p.add_argument("--note", action="append", default=[])
    p.add_argument("--caption-note", default="")
    return p.parse_args(argv)


def run(argv):
    a = parse(argv)
    if os.path.exists(a.out_dir) and os.listdir(a.out_dir):
        raise Refusal(f"--out-dir {a.out_dir} exists and is not empty")
    t = read_tsv(a.switch_map, "switch map")
    order, tracks = A.switch_source_tracks(t)
    cls = read_tsv(a.classes, "classes")
    for c in cls.itertuples():
        if c.example_plan not in set(tracks.plan) or tracks[tracks.plan == c.example_plan].switch_class.iloc[0] != c[1]:
            raise Refusal(f"example plan {c.example_plan} is not a plan of class {c[1]}")
    unknown = sorted(set(tracks.switch_class) - set(cls["class"]))
    if unknown:
        raise Refusal(f"classes in the switch map but not in --classes: {', '.join(unknown)}")
    tracks["example_plan"] = tracks.plan.isin(set(cls.example_plan))
    summary = A.switch_plan_summary(tracks)
    os.makedirs(os.path.join(a.out_dir, "figures"), exist_ok=True)
    layout = []
    rows = [(c.example_plan, f"{c.label}\nexample: {c.example_plan}") for c in cls.itertuples()]
    fig = P.fig_mosaic_design(tracks, order, rows, a.title, notes=a.note)
    plans = [p for c in cls["class"] for p in sorted(summary[summary.switch_class == c].plan)]
    fig2 = P.fig_switch_map_slide(tracks, order, plans, a.map_title)
    for f, stem in ((fig, a.design_stem), (fig2, a.map_stem)):
        clipped, overlapping = P.text_layout_problems(f)
        layout.append({"figure": stem, "width_in": f.get_figwidth(), "height_in": f.get_figheight(), "clipped_texts": len(clipped),
                       "overlapping_text_pairs": len(overlapping), "detail": "; ".join(clipped + overlapping) or "none"})
        P.save(f, os.path.join(a.out_dir, "figures", stem), dpi=P.SLIDE["dpi"])
    write_tsv(tracks[["switch_class", "plan", "example_plan", "homologue", "block_index", "block_id", "source", "switch_before_block",
                      "source_a", "source_b"]], os.path.join(a.out_dir, "mosaic_design_tracks.tsv"))
    write_tsv(summary, os.path.join(a.out_dir, "switch_plans.tsv"))
    write_tsv(pd.DataFrame(layout), os.path.join(a.out_dir, "layout_check.tsv"))
    per = summary.groupby("switch_class", sort=False).plan.size()
    controls = summary[(summary.h1_switches == 0) & (summary.h2_switches == 0)]
    donors = "; ".join(f"{r.plan}: A = {r.source_a}, B = {r.source_b}" for r in summary.itertuples())
    cap = [f"# Captions\n", f"### {a.title}  (`{a.design_stem}`)", "",
           f"{len(summary)} plans in {len(per)} classes ({', '.join(f'{k}: {v}' for k, v in per.items())}); "
           f"{len(summary) - len(controls)} recombine their sources and {len(controls)} have no switch (controls). "
           "Each row shows one example plan of its class; the other plan of the class switches at its own boundaries "
           f"(see `{a.map_stem}`). Both homologues start on source A and toggle at every switch: a later switch returns to "
           "the other source of the same pair, never to a third donor.", "",
           a.caption_note.strip(), "",
           f"### {a.map_title}  (`{a.map_stem}`)", "",
           "Every plan's switch boundaries; markers distinguish a switch on H1 only, on H2 only, or on both. Sources of "
           f"each plan: {donors}.", ""]
    with open(os.path.join(a.out_dir, "captions.md"), "w") as f:
        f.write("\n".join(cap))
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
