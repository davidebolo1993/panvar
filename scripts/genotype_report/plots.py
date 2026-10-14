"""Figures. Every function draws from aggregate tables only (the TSVs in data/), never from runs, so figures can be
restyled without recomputing anything. Colours, orders, jitter and axes are fixed; related figures share them.
"""
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.backends.backend_agg import RendererAgg  # noqa: E402
from matplotlib.font_manager import FontProperties  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

from .common import STATUSES, stable_jitter  # noqa: E402

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
                     "legend.fontsize": 7, "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.spines.top": False,
                     "axes.spines.right": False, "figure.dpi": 100, "savefig.dpi": 150, "pdf.fonttype": 42,
                     "svg.hashsalt": "panvar-genotype-report"})

STATUS_COLORS = {"PASS": "#2f7d4a", "LOW_QUALITY": "#e0a030", "NOT_EVALUATED": "#7c8aa5", "NO_CALL": "#b53b2e",
                 "NOT_APPLICABLE": "#d8d8d8", "NA": "#f2f2f2"}
CALL_COLORS = {"strict_call": "#1d5a85", "genotype_equivalent_call": "#79aed6", "no_call": "#b53b2e", "NA": "#b0b0b0"}
OUTCOME_COLORS = {"retained_correct": "#2f7d4a", "retained_incorrect": "#b53b2e", "filtered_correct": "#9fd0ad",
                  "filtered_incorrect": "#eba59d"}
COVERAGE_COLORS = ["#9ecae1", "#4a98c9", "#1f64a8", "#08306b", "#6a51a3", "#3f007d", "#bcbddc", "#fdae6b"]
CATEGORICAL = ["#1d5a85", "#e0a030", "#2f7d4a", "#b53b2e", "#6a51a3", "#8c564b", "#17becf", "#7f7f7f", "#bcbd22",
               "#e377c2", "#393b79", "#637939", "#8c6d31", "#843c39", "#7b4173"]
SWITCH_MARKERS = {"hap1_only": ("^", "#1d5a85", "haplotype 1 switches"), "hap2_only": ("v", "#e0a030", "haplotype 2 switches"),
                  "both": ("D", "#b53b2e", "both haplotypes switch")}
SYMLOG_NOTE = {"bp": "symlog axis: linear in [0, 1] bp, logarithmic above", "norm": "symlog axis: linear in [0, 1e-4], logarithmic above"}


def save(fig, path_stem, dpi=None):
    fig.savefig(path_stem + ".png", metadata={"Software": None}, **({"dpi": dpi} if dpi else {}))
    fig.savefig(path_stem + ".pdf", metadata={"CreationDate": None, "ModDate": None, "Producer": None, "Creator": None})
    plt.close(fig)


def _legend(ax, **kwargs):
    """A legend only when the axes has labelled artists (an empty legend call warns)."""
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(**kwargs)


def _err(value, low, high):
    """Error-bar lengths from an interval, clipped at 0 against floating-point rounding."""
    return [[max(0.0, value - low)], [max(0.0, high - value)]]


_MEASURE = RendererAgg(10, 10, 720)


def text_width_in(text, size, weight="normal"):
    """The rendered width in inches of the widest line of `text` in the report font (deterministic: the Agg renderer)."""
    prop = FontProperties(family=plt.rcParams["font.family"], size=size, weight=weight)
    return max(_MEASURE.get_text_width_height_descent(line, prop, ismath=False)[0] for line in text.split("\n")) / 720


def block_labels(block_order, lengths=None):
    """Block tick labels in model order: the block ID, and below it the observed truth-allele length range when a
    block_lengths table is supplied (evaluation mode only) and the block has truth-defined observations."""
    if lengths is None or not len(lengths):
        return list(block_order)
    shown = {b: r for b, r in zip(lengths.block_id, lengths.displayed_length_range) if isinstance(r, str) and r != "NA"}
    return [f"{b}\n{shown[b]}" if b in shown else b for b in block_order]


def block_unit_in(labels, size, minimum):
    """Inches per block so horizontal block labels never touch their neighbours."""
    return max(minimum, max(text_width_in(t, size) for t in labels) + 0.12)


def categorical_colors(values):
    return {v: CATEGORICAL[i % len(CATEGORICAL)] for i, v in enumerate(values)}


def coverage_colors(coverages):
    return {c: COVERAGE_COLORS[i % len(COVERAGE_COLORS)] for i, c in enumerate(coverages)}


def cell_labels(table, info, keys=("experiment", "coverage")):
    cells = []
    for e in info.experiments:
        for c in info.coverages:
            if ((table.experiment == e) & (table.coverage == c)).any():
                cells.append((e, c))
    return cells


def _cell_text(e, c):
    return f"{e}\n{c}x" if c != "all" else e


def stacked_status(ax, table, cells, category, order, colors, value, title, denominator="denominator", legend=True):
    """Stacked bars of `value` per (experiment, coverage) cell, stacked by `category` in `order`; the cell's
    denominator is printed above its bar."""
    x = np.arange(len(cells))
    bottom = np.zeros(len(cells))
    for cat in order:
        h = []
        for e, c in cells:
            s = table[(table.experiment == e) & (table.coverage == c) & (table[category] == cat)]
            h.append(s[value].sum())
        h = np.array(h, dtype=float)
        if h.sum() == 0:
            continue
        ax.bar(x, h, bottom=bottom, color=colors.get(cat, "#999999"), width=0.8, label=cat, edgecolor="white", linewidth=0.3)
        bottom += h
    for i, (e, c) in enumerate(cells):
        s = table[(table.experiment == e) & (table.coverage == c)]
        n = int(s[denominator].iloc[0]) if len(s) else 0
        ax.text(i, bottom[i], f"n={n}", ha="center", va="bottom", fontsize=6)
    ax.set_xticks(x, [_cell_text(e, c) for e, c in cells], fontsize=6)
    ax.set_ylabel(value)
    ax.set_title(title)
    ax.set_ylim(0, max(bottom.max() * 1.12, 1))
    if legend:
        _legend(ax, loc="upper left", bbox_to_anchor=(1.0, 1.0), frameon=False)


# ------------------------------------------------------------------------------------------------ QC
def fig_callability(t, info):
    t = t.copy()
    t["category"] = np.where(t.call_class == "no_call", "no_call: " + t.no_call_reason, t.call_class)
    order = [c for c in ["strict_call", "genotype_equivalent_call", "NA"] if c in set(t.category)]
    order += sorted(c for c in set(t.category) if c.startswith("no_call"))
    colors = dict(CALL_COLORS)
    reds = ["#b53b2e", "#7a1f16", "#d9776b", "#e8a79f"]
    for i, c in enumerate(o for o in order if o.startswith("no_call")):
        colors[c] = reds[i % len(reds)]
    cells = cell_labels(t, info)
    fig, ax = plt.subplots(figsize=(max(7, 0.45 * len(cells) + 5), 3.4))
    stacked_status(ax, t, cells, "category", order, colors, "runs", "Call class by experiment and coverage (runs)")
    fig.tight_layout()
    return fig


def fig_status_panels(t, info, filters, unit, title):
    cells = cell_labels(t, info)
    fig, axes = plt.subplots(len(filters), 1, figsize=(max(7, 0.45 * len(cells) + 5), 2.6 * len(filters)), squeeze=False)
    for ax, f in zip(axes[:, 0], filters):
        stacked_status(ax, t[t["filter"] == f], cells, "status", STATUSES + ["NA"], STATUS_COLORS, unit, f"{f} ({unit})",
                       legend=False)
    present = [st for st in STATUSES + ["NA"] if (t.status == st).any()]
    axes[0, 0].legend(handles=[Patch(color=STATUS_COLORS[st], label=st) for st in present], frameon=False,
                      loc="upper left", bbox_to_anchor=(1.0, 1.0))
    fig.suptitle(title, fontsize=9)
    fig.tight_layout()
    return fig


def fig_reasons(t, info):
    filters = ["GENOTYPE_FILTER", "DOSAGE_FILTER"]
    exps = info.experiments
    fig, axes = plt.subplots(1, 2, figsize=(11, 0.5 + 0.22 * max(4, t.groupby("filter").reason.nunique().max() * 1.0) * 1.6))
    for ax, f in zip(axes, filters):
        s = t[(t["filter"] == f) & (t.reason != "none")]
        reasons = sorted(s.reason.unique())
        y = np.arange(len(reasons))
        h = 0.8 / max(1, len(exps))
        cols = categorical_colors(exps)
        for i, e in enumerate(exps):
            v = [s[(s.experiment == e) & (s.reason == r)].blocks.sum() for r in reasons]
            ax.barh(y + (i - (len(exps) - 1) / 2) * h, v, height=h, color=cols[e], label=e)
        ax.set_yticks(y, reasons, fontsize=6)
        ax.invert_yaxis()
        ax.set_xlabel("blocks carrying the reason")
        ax.set_title(f"{f} reasons")
        if f == filters[0]:
            _legend(ax, frameon=False)
    fig.tight_layout()
    return fig


def fig_gq(points, info):
    cells = cell_labels(points, info)
    fig, ax = plt.subplots(figsize=(max(7, 0.45 * len(cells) + 5), 3.4))
    for i, (e, c) in enumerate(cells):
        s = points[(points.experiment == e) & (points.coverage == c)]
        v = s.conditional_GQ.dropna()
        if len(v):
            ax.boxplot([v], positions=[i], widths=0.6, showfliers=False, medianprops={"color": "black"},
                       boxprops={"color": "#1d5a85"}, whiskerprops={"color": "#1d5a85"}, capprops={"color": "#1d5a85"})
            xs = [i + stable_jitter(f"{r}|{b}", 0.5) for r, b in zip(s.run_id[s.conditional_GQ.notna()], s.block_id[s.conditional_GQ.notna()])]
            ax.scatter(xs, v, s=2, color="#1d5a85", alpha=0.25, linewidths=0)
        ax.text(i, 1.02, f"n={len(v)}", transform=ax.get_xaxis_transform(), ha="center", fontsize=6)
    ax.set_xticks(range(len(cells)), [_cell_text(e, c) for e, c in cells], fontsize=6)
    ax.set_ylabel("conditional GQ")
    ax.set_title("Conditional GQ of reported blocks", pad=14)
    fig.tight_layout()
    return fig


def fig_support(gs, ps, info):
    cells_g, cells_p = cell_labels(gs, info), cell_labels(ps, info)
    fig, axes = plt.subplots(2, 1, figsize=(max(7, 0.45 * max(len(cells_g), len(cells_p)) + 6), 5.4))
    g_order = sorted(gs.genotype_support.unique())
    stacked_status(axes[0], gs, cells_g, "genotype_support", g_order, categorical_colors(g_order), "blocks",
                   "Genotype support (blocks)")
    ps = ps.copy()
    ps["category"] = ps.call_class + ": " + ps.phase_support
    p_order = sorted(ps.category.unique())
    stacked_status(axes[1], ps, cells_p, "category", p_order, categorical_colors(p_order), "runs",
                   "Phase support by call class (runs)")
    fig.tight_layout()
    return fig


def fig_dosage_states(t, info, outcomes):
    cells = cell_labels(t, info)
    fig, ax = plt.subplots(figsize=(max(7, 0.45 * len(cells) + 5), 3.6))
    order = [o for o in outcomes if o in set(t.dosage_outcome)]
    stacked_status(ax, t, cells, "dosage_outcome", order, categorical_colors(outcomes), "blocks",
                   "Dosage-stage state of reported blocks")
    fig.tight_layout()
    return fig


def fig_repair_rounds(t, info):
    cells = cell_labels(t, info)
    fig, ax = plt.subplots(figsize=(max(7, 0.45 * len(cells) + 5), 3.2))
    for i, (e, c) in enumerate(cells):
        s = t[(t.experiment == e) & (t.coverage == c)]
        for cc, g in s.groupby("call_class", sort=True):
            v = g.repair_rounds_run
            ok = v.notna()
            xs = [i + stable_jitter(r, 0.6) for r in g.run_id[ok]]
            ax.scatter(xs, v[ok], s=8, color=CALL_COLORS.get(cc, "#999999"), label=cc, linewidths=0)
        ax.text(i, 1.02, f"n={len(s)}", transform=ax.get_xaxis_transform(), ha="center", fontsize=6)
    handles = [Line2D([], [], marker="o", ls="", color=CALL_COLORS[c], label=c) for c in CALL_COLORS if c in set(t.call_class)]
    ax.legend(handles=handles, frameon=False, loc="upper left", bbox_to_anchor=(1, 1))
    ax.set_xticks(range(len(cells)), [_cell_text(e, c) for e, c in cells], fontsize=6)
    ax.set_ylabel("repair rounds run")
    ax.set_title("Adaptive repair rounds per run", pad=14)
    fig.tight_layout()
    return fig


def fig_runtime(t, info):
    cells = cell_labels(t, info)
    groups = sorted(t.compat_group.unique())
    gcol = categorical_colors(groups)
    markers = {"strict_call": "o", "genotype_equivalent_call": "s", "no_call": "x", "NA": "^"}
    panels = [("seconds_total", "wall time (s, log)", True), ("footprint_peak_gib", "peak physical footprint (GiB)", False),
              ("repair_rounds_run", "repair rounds run", False), ("seconds_dosage_total", "dosage stage (s, log)", True)]
    panels = [p for p in panels if t[p[0]].notna().any()]
    fig, axes = plt.subplots(len(panels), 1, figsize=(max(9, 0.45 * len(cells) + 7), 2.3 * len(panels) + 0.4), squeeze=False)
    for ax, (col, label, logy) in zip(axes[:, 0], panels):
        for i, (e, c) in enumerate(cells):
            s = t[(t.experiment == e) & (t.coverage == c) & t[col].notna()]
            ng = s.compat_group.nunique()
            for gi, (g, gg) in enumerate(s.groupby("compat_group", sort=True)):
                off = (gi - (ng - 1) / 2) * (0.7 / max(ng, 1))
                for cc, h in gg.groupby("call_class", sort=True):
                    xs = [i + off + stable_jitter(r, 0.5 / max(ng, 1)) for r in h.run_id]
                    ax.scatter(xs, h[col], s=7, color=gcol[g], marker=markers.get(cc, "o"), linewidths=0.6 if cc == "no_call" else 0)
            ax.text(i, 1.02, f"n={len(s)}", transform=ax.get_xaxis_transform(), ha="center", fontsize=6)
        if logy:
            ax.set_yscale("log")
        ax.set_ylabel(label)
        ax.set_xticks(range(len(cells)), [_cell_text(e, c) for e, c in cells], fontsize=6)
    handles = [Patch(color=gcol[g], label=g) for g in groups]
    handles += [Line2D([], [], marker=markers[c], ls="", color="#444444", label=c) for c in markers if c in set(t.call_class)]
    axes[0, 0].legend(handles=handles, frameon=False, loc="upper left", bbox_to_anchor=(1, 1))
    fig.suptitle("Runtime and resources per run", fontsize=9)
    fig.tight_layout()
    return fig


def fig_library(t, info):
    t = t.copy()
    t["category"] = t.library_mode + ": " + t.library_sensitivity_status
    cells = cell_labels(t, info)
    order = sorted(t.category.unique())
    fig, ax = plt.subplots(figsize=(max(7, 0.45 * len(cells) + 5), 3.2))
    stacked_status(ax, t, cells, "category", order, categorical_colors(order), "runs",
                   "Library source and library-sensitivity status (runs)")
    fig.tight_layout()
    return fig


def fig_provenance(t, info):
    t = t.copy()
    t["label"] = ("preset " + t.preset + " (" + t.preset_source + ") | profile " + t.profile_name + " | " + t.settings_digest
                  + " | threads " + t.threads_requested)
    labels = sorted(t.label.unique())
    exps = info.experiments
    fig, ax = plt.subplots(figsize=(10, 0.8 + 0.35 * len(labels)))
    y = np.arange(len(labels))
    left = np.zeros(len(labels))
    cols = categorical_colors(exps)
    for e in exps:
        v = np.array([t[(t.label == l) & (t.experiment == e)].runs.sum() for l in labels], dtype=float)
        ax.barh(y, v, left=left, color=cols[e], label=e)
        left += v
    ax.set_yticks(y, labels, fontsize=6)
    ax.invert_yaxis()
    ax.set_xlabel("runs")
    ax.set_title("Settings provenance: preset, profile and settings digest (runs)")
    _legend(ax, frameon=False, loc="lower right")
    fig.tight_layout()
    return fig


# ------------------------------------------------------------------------------------------------ ordered blocks
def _symlog(ax, metric, top=None):
    """Symlog with the floor slightly below zero, so zero-valued calls sit visibly above the axis line rather than on it;
    `top` shares one scale across rows."""
    if metric == "normalized_excess":
        ax.set_yscale("symlog", linthresh=1e-4, linscale=0.5)
        ax.set_ylim(-4e-5, max(top or 0, 1e-3) * 1.8)
    else:
        ax.set_yscale("symlog", linthresh=1, linscale=0.5)
        ax.set_ylim(-0.4, max(top or 0, 10) * 1.8)
    ax.axhline(0, color="#bbbbbb", lw=0.5, zorder=0)


COUNT_FONT = 4.5        # the per-coverage count labels only; every other text keeps its size
COUNT_GAP_IN = 0.03
BLOCK_DODGE = 0.9       # the share of a block the dodged coverages span in the ordered-block figures


def _count_labels(summary, set_name, metric, row_key, rows, blocks, covs):
    """{(row, block, coverage): 'n=<n>(<zero>)', or 'n=<n>(<zero>,<undefined>)' when truth-undefined observations are
    present} from block_metric_summary, for the coverages a row has in that table (a block with no row for such a
    coverage has n=0). n, zero and truth_undefined are read, never recomputed."""
    s = summary[(summary["set"] == set_name) & (summary.metric == metric)]
    out = {}
    for _, rval in rows:
        r = s[s[row_key] == rval]
        present = [c for c in covs if (r.coverage == c).any()]
        cells = {(b, c): g for (b, c), g in r.groupby(["block_id", "coverage"], sort=False)}
        for b in blocks:
            for c in present:
                g = cells.get((b, c))
                n, z, u = (0, 0, 0) if g is None else (int(g.n.iloc[0]), int(g.zero.iloc[0]), int(g.truth_undefined.iloc[0]))
                out[(rval, b, c)] = f"n={n}({z},{u})" if u else f"n={n}({z})"
    return out


def count_layout(labels, ncov):
    """Inches per block such that the count labels, all on one baseline, never overlap: the closest label centres
    (adjacent coverages, or the last coverage of a block and the first of the next) are one label width plus a gap
    apart."""
    width = max([text_width_in(t, COUNT_FONT) for t in labels] or [0]) + COUNT_GAP_IN
    w = BLOCK_DODGE / max(1, ncov)
    closest = min(w, 1 - BLOCK_DODGE + w) if ncov > 1 else 1.0
    return width / closest


def fig_ordered_blocks(points, summary, info, metric, label, set_name, rows, row_key, markers=None, lengths=None):
    """rows: [(row label, filter value)]; points: evaluation point table; summary: its block_metric_summary (the source
    of every n / zero / undefined label). Every block in model order, coverages dodged, each call a point, the median a
    black bar; above each block one label per coverage at its dodged position, in its colour, all on one baseline."""
    covs = info.coverages
    ccol = coverage_colors(covs)
    nb = len(info.block_order)
    w = BLOCK_DODGE / max(1, len(covs))
    counts = _count_labels(summary, set_name, metric, row_key, rows, info.block_order, covs)
    ticks = block_labels(info.block_order, lengths)
    unit = max(count_layout(list(counts.values()), len(covs)), block_unit_in(ticks, 7, 0.9))
    band_pt = 2 + COUNT_FONT * 1.25 + 1.5
    fig, axes = plt.subplots(len(rows), 1, figsize=(unit * nb + 2.6, (2.4 + band_pt / 72) * len(rows) + 0.6),
                             squeeze=False, sharex=True)
    shown = points[points[row_key].isin([v for _, v in rows]) & points.truth_defined & points.evaluated]
    top = shown[metric].max() if shown[metric].notna().any() else 0
    for ax, (rlabel, rval) in zip(axes[:, 0], rows):
        s = points[points[row_key] == rval]
        d = s[s.truth_defined & s.evaluated & s[metric].notna()]
        for bi, b in enumerate(info.block_order):
            for ci, c in enumerate(covs):
                g = d[(d.block_id == b) & (d.coverage == c)]
                x0 = bi + (ci - (len(covs) - 1) / 2) * w
                if len(g):
                    xs = [x0 + stable_jitter(f"{r}|{b}", w * 0.7) for r in g.run_id]
                    ax.scatter(xs, g[metric], s=5, color=ccol[c], linewidths=0, alpha=0.8)
                    ax.plot([x0 - w * 0.4, x0 + w * 0.4], [g[metric].median()] * 2, color="black", lw=1)
                if (rval, b, c) in counts:
                    ax.annotate(counts[(rval, b, c)], (x0, 1.0), xycoords=ax.get_xaxis_transform(),
                                xytext=(0, 2), textcoords="offset points", ha="center", va="baseline",
                                fontsize=COUNT_FONT, color=ccol[c], annotation_clip=False, gid="block_count")
        _symlog(ax, metric, top)
        ax.set_xlim(-0.5, nb - 0.5)
        ax.set_ylabel(label, fontsize=7)
        ax.set_title(f"{rlabel} (R={d.run_id.nunique()})", loc="left", fontsize=8, pad=band_pt + 3)   # plotted run IDs
        for bi in range(nb - 1):
            ax.axvline(bi + 0.5, color="#e5e5e5", lw=0.5, zorder=0)
        if markers is not None and rval in markers:
            _switch_markers(ax, markers[rval])
    axes[-1, 0].set_xticks(range(nb), ticks, rotation=0)
    axes[-1, 0].set_xlabel("block (model order)")
    handles = [Line2D([], [], marker="o", ls="", color=ccol[c], label=f"{c}x" if c != "all" else "all") for c in covs]
    handles.append(Line2D([], [], color="black", lw=1, label="median"))
    if markers is not None:
        handles += [Line2D([], [], marker=m, ls="", color=col, label=lab) for m, col, lab in SWITCH_MARKERS.values()]
    axes[0, 0].legend(handles=handles, frameon=False, loc="upper left", bbox_to_anchor=(1.0, 1.0))
    fig.suptitle(f"{_cap(label)} by block \u2014 {set_name}", fontsize=9)
    fig.tight_layout()
    return fig


def _cap(text):
    return text[:1].upper() + text[1:]


def _switch_markers(ax, plan_cats):
    """plan_cats: [(plan, {boundary_index: category})]. One marker per plan at the boundary (between blocks i and
    i+1), stacked downward from the top."""
    trans = ax.get_xaxis_transform()
    for k, (plan, cats) in enumerate(plan_cats):
        for bi, cat in cats.items():
            if cat in SWITCH_MARKERS:
                m, col, _ = SWITCH_MARKERS[cat]
                ax.scatter([bi + 0.5], [0.95 - 0.07 * k], marker=m, color=col, s=14, transform=trans, clip_on=False, zorder=5)


# ------------------------------------------------------------------------------------------------ block rainclouds
# Presentation layout: a fixed 16:9 canvas (13.333 x 7.5 in; PNG 3200 x 1800 px at 240 dpi, PDF 960 x 540 pt), slide-size
# type, no per-cell count labels (their denominators go to a companion table). Values, scales and marks are unchanged.
SLIDE = {"figsize": (40 / 3, 7.5), "dpi": 240, "title": 13, "label": 10, "tick": 9, "block": 8, "legend": 9, "footer": 8.5,
         "raincloud_point": 7, "gq_point": 14}
RAINCLOUD_GRID = 256          # density evaluation points per distribution
RAINCLOUD_HALF = 0.42         # the half violin's maximum width, as a share of one coverage slot
RAINCLOUD_GAP = 0.06          # space between the violin's flat edge and the first jitter column, same unit
RAINCLOUD_JITTER = 0.34       # width of the jittered point band on the open side, same unit
SYMLOG_BP = {"linthresh": 1.0, "linscale": 0.5, "base": 10}   # the accepted bp treatment: linear in [0, 1], log above


def display_transform(scale):
    """(forward, inverse) maps between data values and the axis's displayed coordinate: identity for 'linear', the
    matplotlib symlog transform for 'symlog_bp'. Densities are estimated in the displayed coordinate so their shape is
    the one the reader sees."""
    if scale == "linear":
        return (lambda v: np.asarray(v, dtype=float)), (lambda t: np.asarray(t, dtype=float))
    if scale == "symlog_bp":
        from matplotlib.scale import SymmetricalLogTransform
        fwd = SymmetricalLogTransform(SYMLOG_BP["base"], SYMLOG_BP["linthresh"], SYMLOG_BP["linscale"])
        inv = fwd.inverted()
        return (lambda v: fwd.transform_non_affine(np.asarray(v, dtype=float))), \
               (lambda t: inv.transform_non_affine(np.asarray(t, dtype=float)))
    raise ValueError(f"unknown scale {scale}")


def half_violin(values, scale, grid=RAINCLOUD_GRID):
    """The density outline of one distribution: (y data values, width in [0, 1]) with the maximum width 1, so every
    violin has the same maximum width whatever its sample size.

    A Gaussian kernel density with Scott's bandwidth (sd * n^(-1/5), sample sd) is estimated in the displayed
    coordinate and evaluated only over the observed range [min, max] (no tails beyond the data, so never below zero for
    non-negative values). Fewer than two observations, or zero variance in the displayed coordinate, return the explicit
    degenerate outline: a single flat bar at the common value (kind 'point_mass'), never a kernel spread around it."""
    v = np.asarray(values, dtype=float)
    if not len(v):
        return {"kind": "empty", "y": np.array([]), "width": np.array([])}
    fwd, inv = display_transform(scale)
    t = fwd(v)
    lo, hi = float(t.min()), float(t.max())
    sd = float(t.std(ddof=1)) if len(t) > 1 else 0.0
    if len(t) < 2 or sd == 0.0 or hi == lo:
        y = float(inv(np.array([lo]))[0])
        return {"kind": "point_mass", "y": np.array([y, y]), "width": np.array([1.0, 1.0])}
    h = sd * len(t) ** (-1 / 5)
    g = np.linspace(lo, hi, grid)
    dens = np.exp(-0.5 * ((g[:, None] - t[None, :]) / h) ** 2).sum(axis=1)
    return {"kind": "kde", "y": inv(g), "width": dens / dens.max()}


def quartiles(values):
    """Q1, median, Q3 by linear interpolation between order statistics (numpy's default percentile method)."""
    q1, med, q3 = np.percentile(np.asarray(values, dtype=float), [25, 50, 75])
    return float(q1), float(med), float(q3)


def count_triple_label(c, z, u):
    """The per-block, per-coverage label: always all three values, U included when 0."""
    return f"n=({c},{z},{u})"


def fig_block_raincloud(points, block_order, coverages, panels, counts, title, ticks, xlabel="block (genomic order)", slide=False,
                        footer=None):
    """One figure of vertically aligned panels sharing the ordered block axis; in every block one raincloud per
    coverage: a half violin (left), the same observations as deterministic jittered points on the open (right) side,
    and Q1 / median / Q3 marks across the violin (the median thicker).

    points: one row per plotted observation with run_id, block_id, coverage and one column per panel metric;
    panels: [(metric, y label, 'symlog_bp' | 'linear')]; counts: {(block, coverage): (C, Z, U)} printed above the top
    panel. Every value is plotted (no capping): the axis always reaches the largest observation.
    Returns (figure, drawn) where drawn lists every point as (metric, run_id, block_id, coverage, x, y) and every
    violin as (metric, block_id, coverage, kind, n, q1, median, q3).
    slide=True: the SLIDE canvas and type, no count labels above the axes, the coverage legend above the panels and an
    optional footer line; every value, density, quartile, jitter and axis limit is the same as in the audit layout."""
    covs = list(coverages)
    ccol = coverage_colors(covs)
    nb = len(block_order)
    w = BLOCK_DODGE / max(1, len(covs))
    if slide:
        fig, axes = plt.subplots(len(panels), 1, figsize=SLIDE["figsize"], squeeze=False, sharex=True, layout="constrained")
        counts = {}
    else:
        labels = [count_triple_label(*counts[(b, c)]) for b in block_order for c in covs if (b, c) in counts]
        unit = max(count_layout(labels, len(covs)), block_unit_in(ticks, 7, 0.9))
        band_pt = 2 + COUNT_FONT * 1.25 + 1.5
        fig, axes = plt.subplots(len(panels), 1, figsize=(unit * nb + 2.4, 2.5 * len(panels) + band_pt / 72 + 0.9),
                                 squeeze=False, sharex=True)
    drawn_points, drawn_violins = [], []
    for ax, (metric, ylabel, scale) in zip(axes[:, 0], panels):
        for bi, b in enumerate(block_order):
            for ci, c in enumerate(covs):
                x0 = bi + (ci - (len(covs) - 1) / 2) * w
                g = points[(points.block_id == b) & (points.coverage == c) & points[metric].notna()]
                g = g.sort_values("run_id", kind="stable")
                if len(g):
                    vals = g[metric].to_numpy(dtype=float)
                    hv = half_violin(vals, scale)
                    half = RAINCLOUD_HALF * w
                    if hv["kind"] == "kde":
                        ax.fill_betweenx(hv["y"], x0 - half * hv["width"], x0, color=ccol[c], alpha=0.45, lw=0.4,
                                         edgecolor=ccol[c], zorder=2)
                    else:   # zero variance: a flat bar at the common value, the full violin width
                        ax.plot([x0 - half, x0], hv["y"], color=ccol[c], lw=2.2, solid_capstyle="butt", zorder=2)
                    q1, med, q3 = quartiles(vals)
                    for q, lw in ((q1, 0.7), (q3, 0.7), (med, 1.8)):
                        ax.plot([x0 - half, x0], [q, q], color="#222222", lw=lw, solid_capstyle="butt", zorder=4)
                    band = RAINCLOUD_JITTER * w
                    xs = [x0 + RAINCLOUD_GAP * w + band / 2 + stable_jitter(f"{r}|{b}", band) for r in g.run_id]
                    ax.scatter(xs, vals, s=SLIDE["raincloud_point"] if slide else 4, color=ccol[c], linewidths=0, alpha=0.85, zorder=3)
                    drawn_points += [(metric, r, b, c, x, y) for r, x, y in zip(g.run_id, xs, vals)]
                    drawn_violins.append((metric, b, c, hv["kind"], len(vals), q1, med, q3))
                if ax is axes[0, 0] and (b, c) in counts:
                    ax.annotate(count_triple_label(*counts[(b, c)]), (x0, 1.0), xycoords=ax.get_xaxis_transform(),
                                xytext=(0, 2), textcoords="offset points", ha="center", va="baseline",
                                fontsize=COUNT_FONT, color=ccol[c], annotation_clip=False, gid="block_count")
        top = points[metric].max() if points[metric].notna().any() else 0.0
        if scale == "symlog_bp":
            _symlog(ax, "excess", top)
        else:
            span = max(float(top), 1e-2)
            ax.set_ylim(-0.03 * span, span * 1.08)
            ax.axhline(0, color="#bbbbbb", lw=0.5, zorder=0)
        ax.set_xlim(-0.5, nb - 0.5)
        ax.set_ylabel(ylabel, fontsize=SLIDE["label"] if slide else 7)
        if slide:
            ax.tick_params(axis="y", labelsize=SLIDE["tick"])
        for bi in range(nb - 1):
            ax.axvline(bi + 0.5, color="#e5e5e5", lw=0.5, zorder=0)
    handles = [Patch(facecolor=ccol[c], edgecolor=ccol[c], alpha=0.6, label=f"{c}x" if c != "all" else "all") for c in covs]
    handles += [Line2D([], [], color="#222222", lw=1.8, label="median"), Line2D([], [], color="#222222", lw=0.7, label="Q1, Q3")]
    if slide:
        axes[-1, 0].set_xticks(range(nb), ticks, rotation=0, fontsize=SLIDE["block"])
        axes[-1, 0].set_xlabel(xlabel, fontsize=SLIDE["label"])
        _slide_frame(fig, title, handles, "coverage", footer)
        return fig, {"points": drawn_points, "violins": drawn_violins}
    axes[0, 0].set_title(" ", pad=band_pt + 2)   # room for the count labels above the top panel
    axes[-1, 0].set_xticks(range(nb), ticks, rotation=0)
    axes[-1, 0].set_xlabel(xlabel)
    axes[0, 0].legend(handles=handles, frameon=False, loc="upper left", bbox_to_anchor=(1.0, 1.0), title="coverage")
    fig.suptitle(title, fontsize=9)
    fig.tight_layout()
    return fig, {"points": drawn_points, "violins": drawn_violins}


def _slide_frame(fig, title, handles, legend_title, footer):
    """Slide title, one legend row above the panels (its title as a leading label, so it never stacks into the figure
    title) and an optional footer line, all inside the constrained layout."""
    fig.suptitle(title, fontsize=SLIDE["title"])
    lead = Line2D([], [], ls="", marker="", label=f"{legend_title}:")
    top = fig.axes[0]
    top.legend(handles=[lead] + list(handles), loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=len(handles) + 1, frameon=False,
               fontsize=SLIDE["legend"], handlelength=1.4, columnspacing=1.2, borderaxespad=0.2)
    if footer:
        fig.get_layout_engine().set(rect=(0, 0.035, 1, 0.965))
        fig.text(0.005, 0.006, footer, ha="left", va="bottom", fontsize=SLIDE["footer"], color="#333333", gid="slide_footer")


def text_layout_problems(fig):
    """Every visible text of a drawn figure (titles, axis labels, tick labels, legend entries, colour-bar text, footer):
    the ones outside the canvas, and the pairs of distinct texts whose rendered boxes overlap. Returns (clipped,
    overlapping) as lists of strings; both empty means nothing is clipped or overlapping."""
    fig.canvas.draw()
    rd = fig.canvas.get_renderer()
    undrawn = set()   # tick labels outside their axis's view interval exist but are never drawn
    for ax in fig.axes:
        for axis in (ax.xaxis, ax.yaxis):
            lo, hi = sorted(axis.get_view_interval())
            for tick in axis.get_major_ticks() + axis.get_minor_ticks():
                if not (lo - 1e-12 * max(1.0, abs(lo)) <= tick.get_loc() <= hi + 1e-12 * max(1.0, abs(hi))):
                    undrawn.update((id(tick.label1), id(tick.label2)))
    texts = [t for t in fig.findobj(matplotlib.text.Text) if t.get_visible() and t.get_text().strip() and id(t) not in undrawn]
    boxes = [(t, t.get_window_extent(rd)) for t in texts]
    W, H = fig.bbox.width, fig.bbox.height
    clipped = [t.get_text() for t, b in boxes if b.x0 < -0.5 or b.y0 < -0.5 or b.x1 > W + 0.5 or b.y1 > H + 0.5]
    overlapping = []
    for i, (t, a) in enumerate(boxes):
        for u, b in boxes[i + 1:]:
            if a.overlaps(b) and a.width > 0 and b.width > 0:
                overlapping.append(f"{t.get_text()!r} / {u.get_text()!r}")
    return clipped, overlapping


# ------------------------------------------------------------------------------------------------ block GQ points
COVERAGE_MARKERS = ["o", "s", "^", "D", "v", "P", "X", "h"]
GQ_CMAP = "viridis"
GQ_MISSING_COLOR = "#9a9a9a"
GQ_POINT_JITTER = 0.7   # share of one coverage slot the jittered points span, centred on the slot


def gq_color_value(gq):
    """The colour coordinate log10(1 + GQ): NaN for a missing GQ, +inf for an infinite one."""
    g = np.asarray(gq, dtype=float)
    with np.errstate(invalid="ignore"):
        return np.where(np.isnan(g), np.nan, np.where(np.isposinf(g), np.inf, np.log10(1.0 + g)))


def gq_color_scale(gq):
    """(norm, colormap, colorbar ticks as actual GQ values). The scale spans log10(1 + GQ) from 0 to the largest finite
    GQ; ticks are 0 and the powers of ten 1, 10, 100, ... up to the largest finite GQ."""
    v = gq_color_value(gq)
    finite = v[np.isfinite(v)]
    top = float(finite.max()) if len(finite) else 1.0
    top = max(top, np.log10(2.0))
    norm = matplotlib.colors.Normalize(vmin=0.0, vmax=top)
    ticks = [0.0] + [10.0 ** k for k in range(0, 12) if np.log10(1 + 10.0 ** k) <= top + 1e-12]
    return norm, matplotlib.colormaps[GQ_CMAP], ticks


def gq_point_colors(gq, norm, cmap):
    """Hex RGBA per point: missing -> the neutral grey, +inf -> the top of the scale, finite -> cmap(norm(log10(1+GQ)))."""
    v = gq_color_value(gq)
    out = []
    for x in v:
        if np.isnan(x):
            out.append(GQ_MISSING_COLOR)
        else:
            out.append(matplotlib.colors.to_hex(cmap(norm(norm.vmax if np.isinf(x) else x)), keep_alpha=False))
    return out


def fig_block_gq(points, block_order, rows, coverages, ticks, title, metric="excess", ylabel="excess (bp)",
                 xlabel="block (genomic order)", cbar_label="Conditional GQ", slide=False, footer=None):
    """Vertically stacked panels (rows: [(row value, panel label)]), ordered blocks on x, `metric` on the symlog bp axis
    with each panel's limits set exactly as in fig_block_raincloud's top panel for the same row, coverage dodged and
    encoded by marker shape, every point coloured continuously by conditional GQ through log10(1 + GQ) with the colour
    bar labelled in GQ. Missing GQ: grey; +inf: top colour with a black outline. No point is dropped for its GQ.
    Returns (figure, drawn): one (row, run_id, block_id, coverage, x, y, gq, colour value, hex colour, gq_class) per
    point; figure.panel_axes lists (row, ymin, ymax, yscale) per panel.
    slide=True: the SLIDE canvas and type, the coverage legend above the panels, larger markers and an optional footer;
    every value, x position, colour and axis limit is the same as in the audit layout."""
    covs = list(coverages)
    marker = {c: COVERAGE_MARKERS[i % len(COVERAGE_MARKERS)] for i, c in enumerate(covs)}
    nb = len(block_order)
    w = BLOCK_DODGE / max(1, len(covs))
    norm, cmap, cticks = gq_color_scale(points.conditional_GQ)
    if slide:
        fig, axes = plt.subplots(len(rows), 1, figsize=SLIDE["figsize"], squeeze=False, sharex=True, layout="constrained")
    else:
        unit = block_unit_in(ticks, 7, 0.9)
        fig, axes = plt.subplots(len(rows), 1, figsize=(max(unit * nb, 6) + 2.6, 2.4 * len(rows) + 0.9), squeeze=False, sharex=True)
    drawn = []
    any_missing = any_inf = False
    for ax, (rval, rlabel) in zip(axes[:, 0], rows):
        d = points[points.row == rval]
        for bi, b in enumerate(block_order):
            for ci, c in enumerate(covs):
                g = d[(d.block_id == b) & (d.coverage == c)].sort_values("run_id", kind="stable")
                if not len(g):
                    continue
                x0 = bi + (ci - (len(covs) - 1) / 2) * w
                xs = [x0 + stable_jitter(f"{r}|{b}", w * GQ_POINT_JITTER) for r in g.run_id]
                gq = g.conditional_GQ.to_numpy(dtype=float)
                cols = gq_point_colors(gq, norm, cmap)
                cls = ["missing" if np.isnan(x) else "+inf" if np.isposinf(x) else "finite" for x in gq]
                edge = ["black" if k == "+inf" else "none" for k in cls]
                ax.scatter(xs, g[metric], s=SLIDE["gq_point"] if slide else 9, c=cols, marker=marker[c], edgecolors=edge, linewidths=[0.6 if e == "black" else 0 for e in edge],
                           zorder=3)
                any_missing |= "missing" in cls
                any_inf |= "+inf" in cls
                vals = gq_color_value(gq)
                drawn += [(rval, r, b, c, x, y, q, v, col, k) for r, x, y, q, v, col, k in
                          zip(g.run_id, xs, g[metric], gq, vals, cols, cls)]
        top = d[metric].max() if d[metric].notna().any() else 0.0
        _symlog(ax, "excess", top)
        ax.set_xlim(-0.5, nb - 0.5)
        if not slide:
            ax.set_ylabel(ylabel, fontsize=7)
        ax.set_title(rlabel, loc="left", fontsize=SLIDE["label"] if slide else 8)
        if slide:
            ax.tick_params(axis="y", labelsize=SLIDE["tick"])
        for bi in range(nb - 1):
            ax.axvline(bi + 0.5, color="#e5e5e5", lw=0.5, zorder=0)
    handles = [Line2D([], [], marker=marker[c], ls="", color="#555555", label=f"{c}x" if c != "all" else "all") for c in covs]
    if any_missing:
        handles.append(Line2D([], [], marker="o", ls="", color=GQ_MISSING_COLOR, label="GQ missing"))
    if any_inf:
        handles.append(Line2D([], [], marker="o", ls="", color=cmap(1.0), markeredgecolor="black", label="GQ = +∞"))
    if slide:
        axes[-1, 0].set_xticks(range(nb), ticks, rotation=0, fontsize=SLIDE["block"])
        axes[-1, 0].set_xlabel(xlabel, fontsize=SLIDE["label"])
        fig.supylabel(ylabel, fontsize=SLIDE["label"])   # one label for the stacked panels
        _slide_frame(fig, title, handles, "coverage", footer)
        sm = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
        cb = fig.colorbar(sm, ax=list(axes[:, 0]), fraction=0.02, pad=0.01, aspect=35)
        cb.set_ticks([float(np.log10(1 + t)) for t in cticks], labels=[f"{t:g}" for t in cticks])
        cb.ax.tick_params(labelsize=SLIDE["tick"])
        cb.set_label(cbar_label, fontsize=SLIDE["label"])
        fig.panel_axes = [(rval, *ax.get_ylim(), ax.get_yscale()) for ax, (rval, _) in zip(axes[:, 0], rows)]
        return fig, drawn
    axes[-1, 0].set_xticks(range(nb), ticks, rotation=0)
    axes[-1, 0].set_xlabel(xlabel)
    axes[0, 0].legend(handles=handles, frameon=False, loc="upper left", bbox_to_anchor=(1.0, 1.0), title="coverage")
    fig.suptitle(title, fontsize=9)
    fig.tight_layout()
    sm = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
    cb = fig.colorbar(sm, ax=list(axes[:, 0]), fraction=0.015, pad=0.09, aspect=40)
    cb.set_ticks([float(np.log10(1 + t)) for t in cticks], labels=[f"{t:g}" for t in cticks])
    cb.set_label(cbar_label, fontsize=7)
    fig.panel_axes = [(rval, *ax.get_ylim(), ax.get_yscale()) for ax, (rval, _) in zip(axes[:, 0], rows)]
    return fig, drawn


# ------------------------------------------------------------------------------------------------ block accuracy
def fig_accuracy_heatmap(t, info, sets, lengths=None):
    rows = [p for p, _ in info.panels]
    ticks = block_labels(info.block_order, lengths)
    bottom_in = max(0.8, text_width_in("\n".join(ticks), 6) + 0.25)
    fig, axes = plt.subplots(len(rows), len(sets), figsize=(0.55 * len(info.block_order) * len(sets) + 2.5,
                                                             bottom_in + 0.6 + (0.38 * len(info.coverages) + 0.7) * len(rows)),
                             squeeze=False, sharex=True)
    cmap = matplotlib.colormaps["YlGnBu"].with_extremes(bad="#eeeeee")
    for ri, p in enumerate(rows):
        for si, sname in enumerate(sets):
            ax = axes[ri, si]
            s = t[(t.panel == p) & (t["set"] == sname)]
            M = np.full((len(info.coverages), len(info.block_order)), np.nan)
            for ci, c in enumerate(info.coverages):
                for bi, b in enumerate(info.block_order):
                    g = s[(s.coverage == c) & (s.block_id == b)]
                    if len(g) and g.denominator.iloc[0] > 0:
                        M[ci, bi] = g.rate.iloc[0]
                        ax.text(bi, ci, f"{int(g.floor_attained.iloc[0])}\n/{int(g.denominator.iloc[0])}", ha="center",
                                va="center", fontsize=5, color="black" if g.rate.iloc[0] < 0.7 else "white")
            im = ax.imshow(np.ma.masked_invalid(M), cmap=cmap, vmin=0, vmax=1, aspect="auto")
            ax.set_yticks(range(len(info.coverages)), [f"{c}x" if c != "all" else "all" for c in info.coverages])
            ax.set_xticks(range(len(info.block_order)), ticks, rotation=90, fontsize=6)
            ax.set_title(f"{p} \u2014 {sname}", fontsize=8)
    fig.subplots_adjust(left=0.07, right=0.86, top=1 - 0.9 / fig.get_figheight(), bottom=bottom_in / fig.get_figheight(), hspace=0.35, wspace=0.12)
    fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.6, label="floor-attainment rate (excess 0)")
    fig.suptitle("Floor attainment by block and coverage", fontsize=9)
    return fig


def fig_accuracy_intervals(t, info, sets, lengths=None):
    rows = [p for p, _ in info.panels]
    ticks = block_labels(info.block_order, lengths)
    covs = info.coverages
    ccol = coverage_colors(covs)
    w = 0.8 / max(1, len(covs))
    fig, axes = plt.subplots(len(rows), len(sets), figsize=(5.5 * len(sets) + 1.5, 2.2 * len(rows) + 0.6), squeeze=False, sharey=True)
    for ri, p in enumerate(rows):
        for si, sname in enumerate(sets):
            ax = axes[ri, si]
            s = t[(t.panel == p) & (t["set"] == sname)]
            for ci, c in enumerate(covs):
                g = s[s.coverage == c].set_index("block_id")
                xs, ys, lo, hi = [], [], [], []
                for bi, b in enumerate(info.block_order):
                    if b in g.index and g.loc[b, "denominator"] > 0:
                        xs.append(bi + (ci - (len(covs) - 1) / 2) * w)
                        ys.append(g.loc[b, "rate"])
                        lo.append(max(0.0, g.loc[b, "rate"] - g.loc[b, "ci_low"]))
                        hi.append(max(0.0, g.loc[b, "ci_high"] - g.loc[b, "rate"]))
                if xs:
                    ax.errorbar(xs, ys, yerr=[lo, hi], fmt="o", ms=2.5, color=ccol[c], elinewidth=0.7, label=f"{c}x")
            ax.set_ylim(-0.03, 1.03)
            ax.set_xticks(range(len(info.block_order)), ticks, rotation=90, fontsize=6)
            ax.set_title(f"{p} \u2014 {sname}", fontsize=8)
            if si == 0:
                ax.set_ylabel("floor-attainment rate")
    _legend(axes[0, -1], frameon=False, loc="upper left", bbox_to_anchor=(1, 1))
    fig.suptitle("Floor-attainment rate by block and coverage", fontsize=9)
    fig.tight_layout()
    return fig


# ------------------------------------------------------------------------------------------------ filter performance
def fig_filter_counts(t, info):
    t = t[t.stratum.isin(["all", "experiment"])]
    policies = [p for p in t.policy.unique() if "locus" not in p]
    cols = ["all"] + info.experiments
    fig, axes = plt.subplots(len(policies), 1, figsize=(max(7, 1.0 * len(cols) + 3), 2.2 * len(policies)), squeeze=False)
    for ax, p in zip(axes[:, 0], policies):
        s = t[t.policy == p]
        x = np.arange(len(cols))
        bottom = np.zeros(len(cols))
        for k in ("retained_correct", "retained_incorrect", "filtered_correct", "filtered_incorrect"):
            v = np.array([s[(s.stratum == "all") if c == "all" else (s.experiment == c)][k].sum() for c in cols], dtype=float)
            ax.bar(x, v, bottom=bottom, color=OUTCOME_COLORS[k], label=k.replace("_", " "), width=0.75)
            bottom += v
        for i, c in enumerate(cols):
            ax.text(i, bottom[i], f"n={int(bottom[i])}", ha="center", va="bottom", fontsize=6)
        ax.set_xticks(x, cols, fontsize=7)
        ax.set_title(p, fontsize=8, loc="left")
        ax.set_ylabel("units")
        ax.set_ylim(0, max(bottom.max() * 1.15, 1))
    _legend(axes[0, 0], frameon=False, loc="upper left", bbox_to_anchor=(1, 1))
    fig.tight_layout()
    return fig


def fig_filter_rates(t, info):
    t = t[t.stratum == "experiment_coverage"]
    policies = [p for p in t.policy.unique() if "locus" not in p]
    cells = [(e, c) for e in info.experiments for c in info.coverages if ((t.experiment == e) & (t.coverage == c)).any()]
    metrics = [("retention", "#1d5a85", "o"), ("precision_retained", "#2f7d4a", "s"), ("bad_call_detection_sensitivity", "#b53b2e", "^")]
    fig, axes = plt.subplots(len(policies), 1, figsize=(max(7, 0.45 * len(cells) + 3), 2.2 * len(policies)), squeeze=False, sharex=True)
    for ax, p in zip(axes[:, 0], policies):
        s = t[t.policy == p]
        for m, col, mk in metrics:
            v = [s[(s.experiment == e) & (s.coverage == c)][m].sum() if ((s.experiment == e) & (s.coverage == c)).any() else np.nan
                 for e, c in cells]
            ax.plot(range(len(cells)), v, marker=mk, ls="", color=col, ms=4, label=m.replace("_", " "))
        ax.set_ylim(-0.03, 1.05)
        ax.set_title(p, fontsize=8, loc="left")
    axes[-1, 0].set_xticks(range(len(cells)), [_cell_text(e, c) for e, c in cells], fontsize=6)
    _legend(axes[0, 0], frameon=False, loc="upper left", bbox_to_anchor=(1, 1))
    fig.suptitle("Retention, precision among retained, and bad-call detection sensitivity by experiment and coverage", fontsize=9)
    fig.tight_layout()
    return fig


def fig_filter_catastrophic(t, info, catastrophic_bp):
    t = t[(t.stratum.isin(["all", "experiment"])) & t.catastrophic_retained.notna()]
    policies = list(t.policy.unique())
    cols = ["all"] + info.experiments
    fig, ax = plt.subplots(figsize=(max(7, 0.9 * len(cols) * len(policies) / 2 + 3), 3.2))
    w = 0.8 / max(1, len(policies))
    pc = categorical_colors(policies)
    for pi, p in enumerate(policies):
        s = t[t.policy == p]
        for ci, c in enumerate(cols):
            r = s[(s.stratum == "all") if c == "all" else (s.experiment == c)]
            kept, filt = r.catastrophic_retained.sum(), r.catastrophic_filtered.sum()
            x = ci + (pi - (len(policies) - 1) / 2) * w
            ax.bar(x, kept, width=w * 0.9, color=pc[p], label=p if ci == 0 else None)
            ax.bar(x, filt, bottom=kept, width=w * 0.9, color=pc[p], alpha=0.3)
            ax.text(x, kept + filt, f"{int(kept)}/{int(kept + filt)}", ha="center", va="bottom", fontsize=5)
    ax.set_xticks(range(len(cols)), cols)
    ax.set_ylabel(f"catastrophic blocks (excess >= {catastrophic_bp:g} bp)")
    ax.set_title("Catastrophic errors retained and filtered by policy", fontsize=8)
    _legend(ax, frameon=False, loc="upper left", bbox_to_anchor=(1, 1))
    fig.tight_layout()
    return fig


def fig_filter_by_block(t, info, lengths=None):
    t = t[(t.stratum == "experiment_block") & (t.policy == "EXACT_GENOTYPE_USABLE")]
    exps = [e for e in info.experiments if (t.experiment == e).any()]
    ticks = block_labels(info.block_order, lengths)
    unit = block_unit_in(ticks, 7, 0.8)
    fig, axes = plt.subplots(len(exps), 1, figsize=(max(8, unit * len(info.block_order) + 2), 2.0 * len(exps) + 0.5), squeeze=False, sharex=True)
    for ax, e in zip(axes[:, 0], exps):
        s = t[t.experiment == e].set_index("block_id")
        x = np.arange(len(info.block_order))
        ret = [s.loc[b, "retention"] if b in s.index else np.nan for b in info.block_order]
        pre = [s.loc[b, "precision_retained"] if b in s.index else np.nan for b in info.block_order]
        ax.plot(x - 0.1, ret, "o", color="#1d5a85", ms=4, label="retention")
        ax.plot(x + 0.1, pre, "s", color="#2f7d4a", ms=4, label="precision among retained")
        for i, b in enumerate(info.block_order):
            if b in s.index:
                ax.text(i, 1.06, f"{int(s.loc[b, 'units'])}", ha="center", fontsize=5.5)
        ax.set_ylim(-0.03, 1.15)
        ax.set_title(e, fontsize=8, loc="left")
    axes[-1, 0].set_xticks(range(len(info.block_order)), ticks)
    _legend(axes[0, 0], frameon=False, loc="upper left", bbox_to_anchor=(1, 1))
    fig.suptitle("EXACT_GENOTYPE_USABLE retention and precision by block", fontsize=9)
    fig.tight_layout()
    return fig


# ------------------------------------------------------------------------------------------------ mosaic
def fig_mosaic_callability(t, info, classes):
    t = t.copy()
    t["category"] = np.where(t.call_class == "no_call", "no_call: " + t.no_call_reason, t.call_class)
    cells = [(k, c) for k in classes for c in info.coverages if ((t.switch_class == k) & (t.coverage == c)).any()]
    t = t.rename(columns={"switch_class": "experiment"})
    order = [c for c in ["strict_call", "genotype_equivalent_call", "NA"] if c in set(t.category)]
    order += sorted(c for c in set(t.category) if c.startswith("no_call"))
    fig, ax = plt.subplots(figsize=(max(7, 0.45 * len(cells) + 5), 3.4))
    stacked_status(ax, t, cells, "category", order, {**CALL_COLORS, **{o: "#b53b2e" for o in order if o.startswith("no_call")}},
                   "runs", "Call class by mosaic switch class and coverage (runs)")
    fig.tight_layout()
    return fig


def fig_mosaic_excess(points, info, classes):
    covs = info.coverages
    ccol = coverage_colors(covs)
    metrics = [("excess", "excess (bp)")] + ([("normalized_excess", "diploid-length-normalised excess")] if points.normalized_excess.notna().any() else [])
    fig, axes = plt.subplots(len(metrics), 1, figsize=(max(7, 1.2 * len(classes) + 3), 2.8 * len(metrics)), squeeze=False)
    w = 0.8 / max(1, len(covs))
    for ax, (m, lab) in zip(axes[:, 0], metrics):
        d = points[points.truth_defined & points.evaluated & points[m].notna()]
        for ki, k in enumerate(classes):
            for ci, c in enumerate(covs):
                g = d[(d.switch_class == k) & (d.coverage == c)]
                x0 = ki + (ci - (len(covs) - 1) / 2) * w
                if len(g):
                    xs = [x0 + stable_jitter(f"{r}|{b}", w * 0.7) for r, b in zip(g.run_id, g.block_id)]
                    ax.scatter(xs, g[m], s=4, color=ccol[c], linewidths=0, alpha=0.8)
                    ax.plot([x0 - w * 0.4, x0 + w * 0.4], [g[m].median()] * 2, color="black", lw=1)
            ax.text(ki, 1.01, f"n={len(d[d.switch_class == k])}", transform=ax.get_xaxis_transform(), ha="center", fontsize=6)
        _symlog(ax, m, d[m].max() if len(d) else 0)
        ax.set_xticks(range(len(classes)), classes, fontsize=7)
        ax.set_ylabel(lab + " per block")
        ax.set_title(f"{_cap(lab)} by switch class and coverage",
                     fontsize=8, loc="left", pad=12)
    handles = [Line2D([], [], marker="o", ls="", color=ccol[c], label=f"{c}x") for c in covs]
    axes[0, 0].legend(handles=handles, frameon=False, loc="upper left", bbox_to_anchor=(1, 1))
    fig.tight_layout()
    return fig


def fig_switch_map(t, info, classes, lengths=None):
    plans = []
    for k in classes:
        plans += [(k, p) for p in sorted(t[t.switch_class == k].plan.unique())]
    nb = len(info.block_order)
    ticks = block_labels(info.block_order, lengths)
    unit = block_unit_in(ticks, 7, 0.8)
    fig, ax = plt.subplots(figsize=(max(8, unit * nb + 3), 0.6 + 0.3 * len(plans) + (0.15 if lengths is not None else 0)))
    for yi, (k, p) in enumerate(plans):
        s = t[t.plan == p]
        for r in s.itertuples():
            if r.category in SWITCH_MARKERS:
                m, col, _ = SWITCH_MARKERS[r.category]
                ax.scatter([r.boundary_index + 0.5], [yi], marker=m, color=col, s=30)
        ax.axhline(yi, color="#eeeeee", lw=0.5, zorder=0)
    src = {p: (t[t.plan == p].source_a.iloc[0], t[t.plan == p].source_b.iloc[0]) for _, p in plans}
    ax.set_yticks(range(len(plans)), [f"{p} (A {src[p][0]} / B {src[p][1]})" for k, p in plans], fontsize=6)
    ax.invert_yaxis()
    ax.set_xticks(range(nb), ticks)
    ax.set_xlim(-0.5, nb - 0.5)
    for bi in range(nb - 1):
        ax.axvline(bi + 0.5, color="#f0f0f0", lw=0.5, zorder=0)
    ax.set_xlabel("block (model order); markers sit on the boundary where the homologue changes source template")
    handles = [Line2D([], [], marker=m, ls="", color=col, label=lab) for m, col, lab in SWITCH_MARKERS.values()]
    ax.legend(handles=handles, frameon=False, loc="upper left", bbox_to_anchor=(1, 1))
    ax.set_title("Truth switch map by plan", fontsize=8)
    fig.tight_layout()
    return fig


# ------------------------------------------------------------------------------------------------ dosage
def fig_dosage_content(t, info):
    cells = cell_labels(t, info)
    fig, axes = plt.subplots(2, 1, figsize=(max(7, 0.45 * len(cells) + 5), 5.2))
    x = np.arange(len(cells))
    get = lambda col: np.array([t[(t.experiment == e) & (t.coverage == c)][col].sum() for e, c in cells], dtype=float)
    sens, ret = get("dosage_sensitive_blocks"), get("retained")
    axes[0].bar(x, sens, color="#d8d8d8", label="dosage-sensitive blocks")
    axes[0].bar(x, ret, color="#1d5a85", label="retained (DOSAGE_FILTER == PASS)")
    for i in range(len(cells)):
        axes[0].text(i, sens[i], f"{int(ret[i])}/{int(sens[i])}", ha="center", va="bottom", fontsize=6)
    axes[0].set_title("DOSAGE_CONTENT_USABLE retention", fontsize=8, loc="left")
    _legend(axes[0], frameon=False, loc="upper left", bbox_to_anchor=(1, 1))
    bottom = np.zeros(len(cells))
    for k, col in (("retained_content_correct", OUTCOME_COLORS["retained_correct"]),
                   ("retained_content_incorrect", OUTCOME_COLORS["retained_incorrect"]),
                   ("filtered_content_correct", OUTCOME_COLORS["filtered_correct"]),
                   ("filtered_content_incorrect", OUTCOME_COLORS["filtered_incorrect"])):
        v = get(k)
        axes[1].bar(x, v, bottom=bottom, color=col, label=k.replace("_", " "))
        bottom += v
    und = get("content_truth_undefined")
    for i in range(len(cells)):
        axes[1].text(i, bottom[i], f"n={int(bottom[i])}" + (f"\nundef={int(und[i])}" if und[i] else ""), ha="center", va="bottom", fontsize=5.5)
    axes[1].set_title("Dosage content correctness", fontsize=8, loc="left")
    _legend(axes[1], frameon=False, loc="upper left", bbox_to_anchor=(1, 1))
    for ax in axes:
        ax.set_xticks(x, [_cell_text(e, c) for e, c in cells], fontsize=6)
        ax.set_ylim(0, max(ax.get_ylim()[1] * 1.12, 1))
    fig.tight_layout()
    return fig


def fig_content_vs_exact(t, info):
    panels = [p for p, _ in info.panels if (t.panel == p).any()]
    fig, axes = plt.subplots(1, len(panels), figsize=(3.2 * len(panels) + 1, 3.0), squeeze=False)
    for ax, p in zip(axes[0], panels):
        s = t[t.panel == p]
        M = np.zeros((2, 2))
        for i, dc in enumerate(["content correct", "content incorrect"]):
            for j, ex in enumerate(["floor attained", "floor not attained"]):
                M[i, j] = s[(s.dosage_content == dc) & (s.exact_diplotype == ex)].blocks.sum()
        ax.imshow(M, cmap="Blues", vmin=0, vmax=max(M.max(), 1))
        for i in range(2):
            for j in range(2):
                ax.text(j, i, f"{int(M[i, j])}", ha="center", va="center", fontsize=8,
                        color="white" if M[i, j] > M.max() * 0.6 else "black")
        ax.set_xticks([0, 1], ["exact diplotype\n(floor attained)", "not exact"], fontsize=6)
        ax.set_yticks([0, 1], ["content correct", "content incorrect"], fontsize=6)
        ax.set_title(p, fontsize=8)
    fig.suptitle("Dosage content versus exact-diplotype resolution", fontsize=9)
    fig.tight_layout()
    return fig


def fig_representatives(t, info):
    fig, ax = plt.subplots(figsize=(5.5, 4))
    roles = sorted(t.role.unique())
    rc = categorical_colors(roles)
    for r in roles:
        s = t[t.role == r]
        ax.scatter(s.original_distance + 1, s.representative_phased_distance + 1, s=12, color=rc[r], label=r)
    lim = max(t.original_distance.max(), t.representative_phased_distance.max()) + 10
    ax.plot([1, lim], [1, lim], color="#999999", lw=0.7, ls="--")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("reported-genotype distance + 1 (bp, log)")
    ax.set_ylabel("representative phased distance + 1 (bp, log)")
    ax.set_title("NON-AUTHORITATIVE diagnostic: conditional representatives' sequence distance", fontsize=8)
    _legend(ax, frameon=False, fontsize=6)
    fig.tight_layout()
    return fig


# ------------------------------------------------------------------------------------------------ phase
def fig_phase_callability(t, info):
    t = t.copy()
    t["category"] = t.call_class + ": " + t.phase_support
    cells = cell_labels(t, info)
    order = sorted(t.category.unique())
    fig, ax = plt.subplots(figsize=(max(7, 0.45 * len(cells) + 5), 3.4))
    stacked_status(ax, t, cells, "category", order, categorical_colors(order), "runs",
                   "Phase callability: call class and phase_support (runs)")
    fig.tight_layout()
    return fig


def fig_phase_boundary_status(t, info):
    t = t.copy()
    t["category"] = t.status + ": " + t.reason
    cells = cell_labels(t, info)
    order = sorted(t.category.unique())
    fig, ax = plt.subplots(figsize=(max(6, 0.45 * len(cells) + 3.5), 3.6))
    stacked_status(ax, t, cells, "category", order, categorical_colors(order), "boundaries",
                   "Adjacent-block boundaries by assessment status")
    fig.tight_layout()
    return fig


def fig_phase_accuracy(t, info):
    cells = cell_labels(t, info)
    fig, axes = plt.subplots(2, 1, figsize=(max(7, 0.45 * len(cells) + 5), 5.0), sharex=True)
    x = np.arange(len(cells))
    get = lambda col: np.array([t[(t.experiment == e) & (t.coverage == c)][col].sum() for e, c in cells], dtype=float)
    ok, err = get("correct"), get("switch_errors")
    axes[0].bar(x, ok, color="#2f7d4a", label="correct boundaries")
    axes[0].bar(x, err, bottom=ok, color="#b53b2e", label="switch errors")
    for i in range(len(cells)):
        axes[0].text(i, ok[i] + err[i], f"{int(err[i])}/{int(ok[i] + err[i])}", ha="center", va="bottom", fontsize=6)
    axes[0].set_ylabel("assessable boundaries")
    axes[0].set_ylim(0, max((ok + err).max() * 1.15, 1))
    _legend(axes[0], frameon=False, loc="upper left", bbox_to_anchor=(1, 1))
    axes[0].set_title("Switch errors among assessable boundaries", fontsize=8, loc="left")
    for i, (e, c) in enumerate(cells):
        r = t[(t.experiment == e) & (t.coverage == c)].iloc[0]
        if r.assessable > 0:
            axes[1].errorbar([i], [r.switch_error_rate], yerr=_err(r.switch_error_rate, r.ci_low, r.ci_high),
                             fmt="o", color="#b53b2e", ms=3)
    axes[1].set_ylim(-0.02, max(0.2, axes[1].get_ylim()[1]))
    axes[1].set_ylabel("switch-error rate (95% Wilson)")
    axes[1].set_xticks(x, [_cell_text(e, c) for e, c in cells], fontsize=6)
    fig.tight_layout()
    return fig


def fig_phase_by_class(t, info, classes):
    cells = [(k, c) for k in classes for c in info.coverages if ((t.switch_class == k) & (t.coverage == c)).any()]
    fig, ax = plt.subplots(figsize=(max(7, 0.45 * len(cells) + 5), 3.4))
    for i, (k, c) in enumerate(cells):
        r = t[(t.switch_class == k) & (t.coverage == c)].iloc[0]
        if r.assessable > 0:
            ax.errorbar([i], [r.switch_error_rate], yerr=_err(r.switch_error_rate, r.ci_low, r.ci_high),
                        fmt="o", color="#b53b2e", ms=3)
        ax.text(i, 1.01, f"{int(r.switch_errors)}/{int(r.assessable)}", transform=ax.get_xaxis_transform(), ha="center", fontsize=5.5)
    ax.set_xticks(range(len(cells)), [f"{k}\n{c}x" for k, c in cells], fontsize=5.5, rotation=90)
    ax.set_ylim(-0.02, 1.0)
    ax.set_ylabel("switch-error rate (95% Wilson)")
    ax.set_title("Switch-error rate by mosaic switch class and coverage", fontsize=8, pad=14)
    fig.tight_layout()
    return fig


def fig_phase_by_boundary_type(t, info):
    exps = [e for e in info.experiments if (t.experiment == e).any()]
    cats = ["none", "hap1_only", "hap2_only", "both", "no switch plan"]
    fig, ax = plt.subplots(figsize=(max(7, 1.2 * len(exps) * 2 + 2), 3.4))
    w = 0.8 / len(cats)
    cc = categorical_colors(cats)
    for ei, e in enumerate(exps):
        for ki, k in enumerate(cats):
            r = t[(t.experiment == e) & (t.truth_switch_category == k)]
            if not len(r):
                continue
            r = r.iloc[0]
            x = ei + (ki - (len(cats) - 1) / 2) * w
            ax.bar(x, r.correct, width=w * 0.9, color=cc[k], label=k if ei == 0 else None)
            ax.bar(x, r.switch_errors, bottom=r.correct, width=w * 0.9, color="#b53b2e")
            ax.text(x, r.correct + r.switch_errors, f"{int(r.switch_errors)}/{int(r.assessable)}", ha="center", va="bottom", fontsize=5)
    ax.set_xticks(range(len(exps)), exps)
    ax.set_ylabel("assessable boundaries (red: switch errors)")
    ax.set_title("Phase at boundaries by truth switch type", fontsize=8)
    _legend(ax, frameon=False, loc="upper left", bbox_to_anchor=(1, 1))
    fig.tight_layout()
    return fig


def fig_phase_locus(t, info):
    cells = cell_labels(t, info)
    fig, ax = plt.subplots(figsize=(max(7, 0.45 * len(cells) + 5), 3.4))
    x = np.arange(len(cells))
    get = lambda col: np.array([t[(t.experiment == e) & (t.coverage == c)][col].sum() for e, c in cells], dtype=float)
    strict, ge, nc, exact = get("strict_calls"), get("genotype_equivalent_unphased"), get("no_calls"), get("exact_phased_locus")
    ax.bar(x, exact, color="#2f7d4a", label="strict: exact phased locus")
    ax.bar(x, strict - exact, bottom=exact, color="#79aed6", label="strict: not exact phased locus")
    ax.bar(x, ge, bottom=strict, color="#d8d8d8", label="genotype-equivalent (unphased, excluded)")
    ax.bar(x, nc, bottom=strict + ge, color="#b53b2e", label="no-call")
    for i in range(len(cells)):
        ax.text(i, strict[i] + ge[i] + nc[i], f"n={int(strict[i] + ge[i] + nc[i])}", ha="center", va="bottom", fontsize=6)
    ax.set_xticks(x, [_cell_text(e, c) for e, c in cells], fontsize=6)
    ax.set_ylabel("runs")
    ax.set_title("Full-locus phased-call status", fontsize=8)
    _legend(ax, frameon=False, loc="upper left", bbox_to_anchor=(1, 1))
    fig.tight_layout()
    return fig


def fig_phase_projected(t, info):
    cells = cell_labels(t, info)
    fig, ax = plt.subplots(figsize=(max(7, 0.45 * len(cells) + 5), 3.2))
    for i, (e, c) in enumerate(cells):
        r = t[(t.experiment == e) & (t.coverage == c)].iloc[0]
        if r.evaluated > 0:
            ax.errorbar([i], [r.concordance], yerr=_err(r.concordance, r.ci_low, r.ci_high), fmt="s", color="#6a51a3", ms=3)
        ax.text(i, 1.01, f"{int(r.concordant)}/{int(r.evaluated)}", transform=ax.get_xaxis_transform(), ha="center", fontsize=5.5)
    ax.set_xticks(range(len(cells)), [_cell_text(e, c) for e, c in cells], fontsize=6)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("sequence-projected concordance")
    ax.set_title("SECONDARY diagnostic: sequence-projected phase concordance", fontsize=8, pad=14)
    fig.tight_layout()
    return fig


# ------------------------------------------------------------------------------------------------ GQ diagnostic
def fig_gq_diagnostic(t, info, labels):
    panels = [p for p, _ in info.panels if (t.panel == p).any()]
    fig, axes = plt.subplots(1, len(panels), figsize=(3.6 * len(panels) + 1, 3.4), squeeze=False, sharey=True)
    for ax, p in zip(axes[0], panels):
        s = t[t.panel == p].set_index("gq_bin")
        bins = [b for b in labels + ["GQ not reported"] if b in s.index]
        for i, b in enumerate(bins):
            r = s.loc[b]
            ax.errorbar([i], [r.rate], yerr=_err(r.rate, r.ci_low, r.ci_high), fmt="o", color="#1d5a85", ms=3)
            ax.text(i, 1.02, f"{int(r.floor_attained)}/{int(r.denominator)}", transform=ax.get_xaxis_transform(), ha="center",
                    fontsize=5, rotation=90, va="bottom")
        ax.set_xticks(range(len(bins)), bins, rotation=90, fontsize=6)
        ax.set_ylim(-0.03, 1.03)
        ax.set_title(p, fontsize=8, pad=26)
    axes[0, 0].set_ylabel("floor-attainment rate (95% Wilson)")
    fig.suptitle("DIAGNOSTIC ONLY: floor attainment by conditional-GQ bin", fontsize=9)
    fig.tight_layout()
    return fig
