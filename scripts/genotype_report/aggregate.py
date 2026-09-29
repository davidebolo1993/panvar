"""Aggregation: pure functions from the normalised tables to the TSVs every figure is drawn from. No plotting here, so a
figure can be restyled from data/ without recomputing anything.

Grouping keys are `strata + [experiment, coverage]`, where strata is ["modality"] when the run metadata supplies a
modality and [] otherwise (a modality is never assumed). Every aggregate keeps its denominator beside its numerator.
"""
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from .common import CALL_CLASSES, CALLED, STATUSES, Refusal, coverage_order, ordered_unique, to_number, wilson

GQ_BINS = [0, 10, 20, 30, 40, 60, 80, 99, np.inf]
GQ_LABELS = ["[0,10)", "[10,20)", "[20,30)", "[30,40)", "[40,60)", "[60,80)", "[80,99)", ">=99"]
DOSAGE_OUTCOMES = ["complete repair applied", "complete repair available, not applied (audit)",
                   "complete repair not applied (guard)", "content-only correction", "triggered, unrepaired", "SUSPECT",
                   "UNDERPOWERED", "CONSISTENT", "not dosage-sensitive", "not evaluated"]


@dataclass
class Context:
    runs: pd.DataFrame
    blocks: pd.DataFrame
    block_order: list
    experiments: list
    coverages: list
    panels: list                     # [(panel name, [experiments])]
    strata: list                     # ["modality"] or []
    modalities: list
    has_headline: bool
    catastrophic_bp: float
    evaluation: pd.DataFrame = None
    has_normalization: bool = False
    content: pd.DataFrame = None
    representatives: pd.DataFrame = None
    switch_plans: dict = None
    switch_key: str = None
    phase_boundaries: pd.DataFrame = None
    phase_calls: pd.DataFrame = None
    phase_projected: pd.DataFrame = None
    notes: list = field(default_factory=list)

    @property
    def G(self):
        return self.strata + ["experiment", "coverage"]


def parse_panel_groups(text, experiments):
    """'NAME=exp1,exp2;NAME2=exp3': the rows of the ordered-block figures. Every experiment must be placed exactly
    once; without the option each experiment is its own row."""
    if not text:
        return [(e, [e]) for e in experiments]
    panels, placed = [], []
    for part in [p for p in text.split(";") if p.strip()]:
        if "=" not in part:
            raise Refusal(f"--panel-groups: {part!r} is not NAME=experiment[,experiment]")
        name, members = part.split("=", 1)
        members = [m.strip() for m in members.split(",") if m.strip()]
        unknown = [m for m in members if m not in experiments]
        if unknown:
            raise Refusal(f"--panel-groups: unknown experiment(s) {', '.join(unknown)} (known: {', '.join(experiments)})")
        panels.append((name.strip(), members))
        placed += members
    dup = sorted({m for m in placed if placed.count(m) > 1})
    left = [e for e in experiments if e not in placed]
    if dup or left:
        raise Refusal("--panel-groups must place every experiment exactly once"
                      + (f"; repeated: {', '.join(dup)}" if dup else "") + (f"; missing: {', '.join(left)}" if left else ""))
    return panels


def build_context(runs, blocks, block_order, panel_groups=None, catastrophic_bp=1000.0):
    runs = runs.copy()
    if "_meta_order" in runs.columns:
        runs = runs.sort_values(["_meta_order"], kind="stable")
    experiments = ordered_unique(runs.experiment)
    coverages = coverage_order(runs.coverage)
    strata = ["modality"] if "modality" in runs.columns else []
    modalities = sorted(runs.modality.unique()) if strata else []
    panels = parse_panel_groups(panel_groups, experiments)
    panel_of = {e: p for p, members in panels for e in members}
    runs["panel"] = runs.experiment.map(panel_of)
    if "execution_group" not in runs.columns:
        runs["execution_group"] = "not recorded"
    runs["compat_group"] = ("settings " + runs.settings_digest.str.replace("settings:", "", regex=False).str[:12]
                            + " | threads " + runs.threads_requested + " | " + runs.execution_group)
    if len(modalities) > 1:   # one modality is constant: the label, and so every figure, stays as without it
        runs["compat_group"] = runs.modality + " | " + runs.compat_group
    carry = ["run_id", "experiment", "coverage", "panel", "call_class"] + strata + [
        c for c in ("subject", "execution_group") if c in runs.columns]
    b = blocks.merge(runs[carry], on="run_id", how="left")
    b["conditional_GQ_value"] = to_number(b.conditional_GQ)
    has_headline = bool((b.genotype_filter != "NA").any()) if len(b) else bool((runs.genotype_filter != "NA").any())
    if has_headline and len(b) and (b.genotype_filter == "NA").any():
        raise Refusal("mixed runs: some block tables carry the headline filters and some do not (outputs from before and "
                      "after headline_filters_v2); report them separately")
    if has_headline:
        b["exact_genotype_usable"] = (b.genotype_filter == "PASS") & b.dosage_filter.isin(["PASS", "NOT_APPLICABLE"])
    runs = runs.reset_index(drop=True)
    return Context(runs=runs, blocks=b, block_order=block_order, experiments=experiments, coverages=coverages,
                   panels=panels, strata=strata, modalities=modalities, has_headline=has_headline,
                   catastrophic_bp=catastrophic_bp)


# ------------------------------------------------------------------------------------------------ deterministic order
def sort_frame(frame, ctx, extra=None):
    """Sort by every known key in a fixed precedence and rank; unknown keys alphabetically. Rows never depend on the
    input order."""
    ranks = {"modality": {m: i for i, m in enumerate(ctx.modalities)},
             "panel": {p: i for i, (p, _) in enumerate(ctx.panels)},
             "experiment": {e: i for i, e in enumerate(ctx.experiments)},
             "coverage": {c: i for i, c in enumerate(ctx.coverages)},
             "call_class": {c: i for i, c in enumerate(CALL_CLASSES)},
             "status": {s: i for i, s in enumerate(STATUSES)},
             "set": {"all called blocks": 0, "EXACT_GENOTYPE_USABLE": 1},
             "dosage_outcome": {s: i for i, s in enumerate(DOSAGE_OUTCOMES)},
             "gq_bin": {s: i for i, s in enumerate(GQ_LABELS)}}
    precedence = ["modality", "policy", "stratum", "set", "metric", "filter", "panel", "switch_class", "experiment", "coverage",
                  "compat_group", "block_index", "boundary_index", "call_class", "status", "dosage_outcome", "gq_bin"]
    precedence += list(extra or [])
    keys = [k for k in precedence if k in frame.columns]
    keys += [c for c in frame.columns if c not in keys and frame[c].dtype == object]
    if frame.empty:
        return frame.reset_index(drop=True)
    tmp = frame.copy()
    sort_cols = []
    for k in keys:
        col = "__k_" + k
        if k in ranks:
            tmp[col] = tmp[k].map(lambda v, r=ranks[k]: (r.get(v, len(r)), str(v)))
        elif tmp[k].dtype == object:
            tmp[col] = tmp[k].astype(str)
        else:
            tmp[col] = tmp[k]
        sort_cols.append(col)
    tmp = tmp.sort_values(sort_cols, kind="stable")
    return tmp.drop(columns=sort_cols).reset_index(drop=True)


def _count(frame, keys, name, ctx, denominator_keys=None):
    if frame.empty:
        return pd.DataFrame(columns=keys + [name, "denominator"])
    out = frame.groupby(keys, dropna=False, sort=False).size().reset_index(name=name)
    dk = denominator_keys if denominator_keys is not None else ctx.G
    den = frame.groupby(dk, dropna=False, sort=False).size().reset_index(name="denominator")
    out = out.merge(den, on=dk, how="left")
    return sort_frame(out, ctx)


# ------------------------------------------------------------------------------------------------ truth-free QC
def callability(ctx):
    r = ctx.runs.copy()
    r["no_call_reason"] = np.where(r.call_class.isin(CALLED), "none", r.no_call_reason)
    return _count(r, ctx.G + ["call_class", "no_call_reason"], "runs", ctx)


def headline_block_status(ctx):
    rows = []
    for f in ("genotype", "dosage"):
        t = ctx.blocks.copy()
        t["filter"] = f.upper() + "_FILTER"
        t["status"] = t[f + "_filter"]
        rows.append(t)
    t = pd.concat(rows, ignore_index=True)
    return _count(t, ctx.G + ["filter", "status"], "blocks", ctx, ctx.G + ["filter"])


def headline_locus_status(ctx):
    rows = []
    for f in ("genotype", "dosage", "phase"):
        t = ctx.runs.copy()
        t["filter"] = f.upper() + "_FILTER"
        t["status"] = t[f + "_filter"]
        rows.append(t)
    t = pd.concat(rows, ignore_index=True)
    return _count(t, ctx.G + ["filter", "status"], "runs", ctx, ctx.G + ["filter"])


def headline_block_reasons(ctx):
    rows = []
    for f in ("genotype", "dosage"):
        t = ctx.blocks[ctx.G + ["run_id", "block_id", f + "_filter", f + "_filter_reasons"]].copy()
        t.columns = ctx.G + ["run_id", "block_id", "status", "reasons"]
        t["filter"] = f.upper() + "_FILTER"
        t["reason"] = t.reasons.str.split(",")
        rows.append(t.explode("reason"))
    t = pd.concat(rows, ignore_index=True)
    out = t.groupby(ctx.strata + ["experiment", "filter", "status", "reason"], sort=False).size().reset_index(name="blocks")
    den = ctx.blocks.groupby(ctx.strata + ["experiment"], sort=False).size().reset_index(name="denominator_blocks")
    return sort_frame(out.merge(den, on=ctx.strata + ["experiment"]), ctx)


def gq_values(ctx):
    t = ctx.blocks[ctx.G + ["panel", "run_id", "block_index", "block_id", "conditional_GQ"]].copy()
    t["conditional_GQ"] = to_number(t.conditional_GQ)
    return sort_frame(t, ctx, ["run_id"])


def gq_summary(ctx):
    t = gq_values(ctx)
    rows = []
    for keys, g in t.groupby(ctx.G, sort=False):
        v = g.conditional_GQ.dropna()
        keys = keys if isinstance(keys, tuple) else (keys,)
        rows.append({**dict(zip(ctx.G, keys)), "blocks": len(g), "gq_missing": int(g.conditional_GQ.isna().sum()),
                     "median": v.median() if len(v) else np.nan, "q25": v.quantile(0.25) if len(v) else np.nan,
                     "q75": v.quantile(0.75) if len(v) else np.nan, "min": v.min() if len(v) else np.nan,
                     "max": v.max() if len(v) else np.nan})
    return sort_frame(pd.DataFrame(rows), ctx)


def genotype_support(ctx):
    return _count(ctx.blocks, ctx.G + ["genotype_support"], "blocks", ctx)


def phase_support(ctx):
    return _count(ctx.runs, ctx.G + ["call_class", "phase_support"], "runs", ctx)


def dosage_outcome_of(frame, dosage_mode):
    s, applied = frame.dosage_status, frame.dosage_repair_applied == "true"
    reason = frame.dosage_filter_reasons if "dosage_filter_reasons" in frame.columns else pd.Series("NA", index=frame.index)
    out = np.select(
        [s.eq("REPAIRED") & applied, s.eq("REPAIRED") & reason.str.contains("COMPLETE_REPAIR_AVAILABLE_NOT_APPLIED"),
         s.eq("REPAIRED"), s.eq("CONTENT_REPAIRED"), s.isin(["TRIGGERED_UNREPAIRED", "AUTO_REPAIR"]), s.eq("SUSPECT"),
         s.eq("UNDERPOWERED"), s.eq("CONSISTENT"), s.eq("NOT_SENSITIVE")],
        DOSAGE_OUTCOMES[:9], default="not evaluated")
    return pd.Series(out, index=frame.index)


def dosage_states(ctx):
    t = ctx.blocks.merge(ctx.runs[["run_id", "dosage_mode"]], on="run_id", how="left")
    t["dosage_outcome"] = dosage_outcome_of(t, t.dosage_mode)
    return _count(t, ctx.G + ["dosage_mode", "dosage_outcome", "dosage_status", "repair_class"], "blocks", ctx)


def repair_rounds(ctx):
    t = ctx.runs[ctx.G + ["run_id", "call_class", "repair_rounds_run"]].copy()
    t["repair_rounds_run"] = to_number(t.repair_rounds_run)
    return sort_frame(t, ctx, ["run_id"])


def runtime_points(ctx):
    t = ctx.runs[ctx.G + ["compat_group", "run_id", "call_class", "seconds_total", "footprint_peak_bytes",
                          "repair_rounds_run", "seconds_dosage_total"]].copy()
    for c in ("seconds_total", "footprint_peak_bytes", "repair_rounds_run", "seconds_dosage_total"):
        t[c] = to_number(t[c])
    t["footprint_peak_gib"] = t.footprint_peak_bytes / 2 ** 30
    return sort_frame(t, ctx, ["run_id"])


def runtime_groups(ctx):
    """One row per compatibility group: what its runs share. Pooling across groups never happens; a figure cell that
    holds more than one group shows them separately."""
    rows = []
    for g, r in ctx.runs.groupby("compat_group", sort=True):
        rows.append({**({"modality": ";".join(sorted(r.modality.unique()))} if ctx.strata else {}), "compat_group": g, "runs": len(r), "settings_digest": ";".join(sorted(r.settings_digest.unique())),
                     "profile_name": ";".join(sorted(r.profile_name.unique())), "preset": ";".join(sorted(r.preset.unique())),
                     "threads_requested": ";".join(sorted(r.threads_requested.unique())),
                     "execution_group": ";".join(sorted(r.execution_group.unique())),
                     "distinct_graphs": r.graph_sha256.nunique(),
                     "experiments": ";".join(ordered_unique(r.experiment))})
    return pd.DataFrame(rows)


def runtime_summary(ctx):
    t = runtime_points(ctx)
    rows = []
    for keys, g in t.groupby(ctx.G + ["compat_group"], sort=False):
        row = dict(zip(ctx.G + ["compat_group"], keys))
        row["runs"] = len(g)
        for c in ("seconds_total", "footprint_peak_gib", "repair_rounds_run", "seconds_dosage_total"):
            v = g[c].dropna()
            row[c + "_n"] = len(v)
            row[c + "_median"] = v.median() if len(v) else np.nan
            row[c + "_max"] = v.max() if len(v) else np.nan
        rows.append(row)
    return sort_frame(pd.DataFrame(rows), ctx)


def library(ctx):
    return _count(ctx.runs, ctx.G + ["library_mode", "library_sensitivity_status"], "runs", ctx)


def provenance(ctx):
    t = ctx.runs.copy()
    t["profile_sha256"] = t.profile_sha256.str[:16]
    keys = ctx.strata + ["experiment", "preset", "preset_source", "profile_name", "profile_sha256", "settings_digest",
                         "threads_requested", "execution_group"]
    return _count(t, keys, "runs", ctx, ctx.strata + ["experiment"])


# ------------------------------------------------------------------------------------------------ evaluation
def evaluated_blocks(ctx):
    """Every reported block with its evaluation (point level). evaluated = an evaluator row exists; accuracy
    denominators use truth-defined evaluated blocks only."""
    e = ctx.evaluation
    b = ctx.blocks.merge(e, on=["run_id", "block_id"], how="left", indicator=True)
    b["evaluated"] = b["_merge"] == "both"
    b = b.drop(columns=["_merge"])
    b["truth_defined"] = b.truth_defined.fillna(False).astype(bool)
    b["floor_attained"] = np.where(b.truth_defined, b.excess == 0, np.nan)
    b["exact"] = np.where(b.truth_defined, b.distance == 0, np.nan)
    b["catastrophic"] = np.where(b.truth_defined, b.excess >= ctx.catastrophic_bp, np.nan)
    if ctx.switch_plans is not None:
        b["switch_class"] = _switch_class(ctx, b)
    keep = ctx.G + ["panel", "run_id", "call_class", "block_index", "block_id", "kind", "evaluated", "truth_defined", "distance",
                    "floor", "excess", "normalized_excess", "diploid_length_denominator_bp", "floor_attained", "exact", "catastrophic",
                    "conditional_GQ_value", "genotype_filter", "dosage_filter", "exact_genotype_usable", "switch_class"]
    for c in keep:
        if c not in b.columns:
            b[c] = np.nan
    return sort_frame(b[keep], ctx, ["run_id"])


def _sets(frame, ctx):
    yield "all called blocks", frame
    if ctx.has_headline:
        yield "EXACT_GENOTYPE_USABLE", frame[frame.exact_genotype_usable == True]  # noqa: E712


METRICS = [("distance", "raw diploid distance (bp)"), ("excess", "excess over the availability floor (bp)"),
           ("normalized_excess", "diploid-length-normalised excess")]


def block_metric_summary(ctx, points, by="panel"):
    rows = []
    for metric, _ in METRICS:
        if metric == "normalized_excess" and not ctx.has_normalization:
            continue
        for set_name, s in _sets(points, ctx):
            for keys, g in s.groupby(ctx.strata + [by, "block_index", "block_id", "coverage"], sort=False):
                keys = keys if isinstance(keys, tuple) else (keys,)
                d = g[g.truth_defined & g.evaluated]
                v = d[metric].dropna()
                rows.append({**dict(zip(ctx.strata + [by, "block_index", "block_id", "coverage"], keys)), "set": set_name,
                             "metric": metric, "blocks_reported": len(g), "not_evaluated": int((~g.evaluated).sum()),
                             "truth_undefined": int((g.evaluated & ~g.truth_defined).sum()), "n": len(v),
                             "zero": int((v == 0).sum()), "median": v.median() if len(v) else np.nan,
                             "mean": v.mean() if len(v) else np.nan, "max": v.max() if len(v) else np.nan})
    return sort_frame(pd.DataFrame(rows), ctx)


def block_accuracy(ctx, points, by="panel"):
    rows = []
    for set_name, s in _sets(points, ctx):
        for keys, g in s.groupby(ctx.strata + [by, "block_index", "block_id", "coverage"], sort=False):
            keys = keys if isinstance(keys, tuple) else (keys,)
            d = g[g.truth_defined & g.evaluated]
            k, n = int((d.excess == 0).sum()), len(d)
            lo, hi = wilson(k, n)
            rows.append({**dict(zip(ctx.strata + [by, "block_index", "block_id", "coverage"], keys)), "set": set_name,
                         "blocks_reported": len(g), "not_evaluated": int((~g.evaluated).sum()),
                         "truth_undefined": int((g.evaluated & ~g.truth_defined).sum()), "denominator": n,
                         "floor_attained": k, "rate": k / n if n else np.nan, "ci_low": lo, "ci_high": hi,
                         "exact_distance_zero": int((d.distance == 0).sum())})
    return sort_frame(pd.DataFrame(rows), ctx)


def _performance(units, retained, correct, catastrophic=None):
    """units: boolean index frame restricted to truth-defined units."""
    inc = ~correct
    row = {"units": int(len(units)), "retained_correct": int((retained & correct).sum()),
           "retained_incorrect": int((retained & inc).sum()), "filtered_correct": int((~retained & correct).sum()),
           "filtered_incorrect": int((~retained & inc).sum())}
    n_ret = row["retained_correct"] + row["retained_incorrect"]
    n_inc = row["retained_incorrect"] + row["filtered_incorrect"]
    row["retention"] = n_ret / len(units) if len(units) else np.nan
    row["precision_retained"] = row["retained_correct"] / n_ret if n_ret else np.nan
    row["bad_call_detection_sensitivity"] = row["filtered_incorrect"] / n_inc if n_inc else np.nan
    if catastrophic is not None:
        row["catastrophic_retained"] = int((retained & catastrophic).sum())
        row["catastrophic_filtered"] = int((~retained & catastrophic).sum())
    else:
        row["catastrophic_retained"] = row["catastrophic_filtered"] = np.nan
    return row


def filter_performance(ctx, points):
    """Headline filters and usage policies against the supplied truth. Sequence basis: correct = excess 0 (floor
    attained), catastrophic = excess >= catastrophic_bp, over evaluated truth-defined blocks. Content basis
    (DOSAGE_CONTENT_USABLE): the supplied content correctness, over dosage-sensitive blocks with defined content truth.
    Phase is judged only on phase truth (run level), never folded into the exact-genotype policies."""
    if not ctx.has_headline:
        raise Refusal("the runs carry no headline filters")
    d = points[points.truth_defined & points.evaluated].copy()
    d["correct"] = d.excess == 0
    d["cat"] = d.excess >= ctx.catastrophic_bp
    policies = [("GENOTYPE_FILTER == PASS", "sequence: excess 0", d.genotype_filter == "PASS"),
                ("DOSAGE_FILTER in {PASS, NOT_APPLICABLE}", "sequence: excess 0", d.dosage_filter.isin(["PASS", "NOT_APPLICABLE"])),
                ("EXACT_GENOTYPE_USABLE", "sequence: excess 0", d.exact_genotype_usable == True)]  # noqa: E712
    strata = [("all", []), ("experiment", ["experiment"]), ("experiment_coverage", ["experiment", "coverage"]),
              ("experiment_block", ["experiment", "block_index", "block_id"]), ("block", ["block_index", "block_id"])]
    rows = []
    for name, basis, keep in policies:
        d["_ret"] = keep
        for sname, keys in strata:
            groups = [((), d)] if not keys and not ctx.strata else d.groupby(ctx.strata + keys, sort=False)
            for k, g in groups:
                k = k if isinstance(k, tuple) else (k,)
                rows.append({"policy": name, "correctness_basis": basis, "stratum": sname,
                             **dict(zip(ctx.strata + keys, k)), **_performance(g, g._ret, g.correct, g.cat)})
    if ctx.content is not None:
        c = points.merge(ctx.content, on=["run_id", "block_id"], how="inner")
        sens = c[~c.dosage_filter.isin(["NOT_APPLICABLE", "NO_CALL"]) & c.content_correct.notna()].copy()
        sens["_ret"] = sens.dosage_filter == "PASS"
        sens["_ok"] = sens.content_correct == True  # noqa: E712
        for sname, keys in strata:
            groups = [((), sens)] if not keys and not ctx.strata else sens.groupby(ctx.strata + keys, sort=False)
            for k, g in groups:
                k = k if isinstance(k, tuple) else (k,)
                rows.append({"policy": "DOSAGE_CONTENT_USABLE", "correctness_basis": "dosage content (supplied)",
                             "stratum": sname, **dict(zip(ctx.strata + keys, k)), **_performance(g, g._ret, g._ok)})
    if ctx.phase_calls is not None:
        pc = ctx.phase_calls.merge(ctx.runs[ctx.G + ["run_id", "call_class", "genotype_filter", "dosage_filter", "phase_filter"]],
                                   on="run_id")
        pc = pc[(pc.call_class == "strict_call") & (pc.assessed_boundaries > 0)].copy()
        run_policies = [("PHASE_FILTER == PASS (locus)", "phase: every assessed boundary correct", pc.phase_filter == "PASS",
                         pc.all_assessable_correct == True),  # noqa: E712
                        ("PHASED_GENOTYPE_USABLE (locus)", "phase and genotype: exact phased locus",
                         (pc.genotype_filter == "PASS") & pc.dosage_filter.isin(["PASS", "NOT_APPLICABLE"]) & (pc.phase_filter == "PASS"),
                         pc.exact_phased_locus == True)]  # noqa: E712
        for name, basis, keep, ok in run_policies:
            pc["_ret"], pc["_ok"] = keep, ok
            for sname, keys in strata[:3]:
                groups = [((), pc)] if not keys and not ctx.strata else pc.groupby(ctx.strata + keys, sort=False)
                for k, g in groups:
                    k = k if isinstance(k, tuple) else (k,)
                    rows.append({"policy": name, "correctness_basis": basis, "stratum": sname + " (strict-call runs)",
                                 **dict(zip(ctx.strata + keys, k)), **_performance(g, g._ret, g._ok)})
    out = pd.DataFrame(rows)
    for c in ("experiment", "coverage", "block_index", "block_id"):
        if c not in out.columns:
            out[c] = ""
    out[["experiment", "coverage", "block_id"]] = out[["experiment", "coverage", "block_id"]].fillna("")
    out["block_index"] = out.block_index.fillna(-1).astype(int)
    order = {p: i for i, p in enumerate(ordered_unique(out.policy))}
    sorder = {s: i for i, s in enumerate(ordered_unique(out.stratum))}
    out["_p"], out["_s"] = out.policy.map(order), out.stratum.map(sorder)
    exp_r = {e: i for i, e in enumerate([""] + ctx.experiments)}
    cov_r = {c: i for i, c in enumerate([""] + ctx.coverages)}
    out["_e"], out["_c"] = out.experiment.map(exp_r), out.coverage.map(cov_r)
    lead = ["modality"] if ctx.strata else []
    out = out.sort_values(lead + ["_p", "_s", "_e", "_c", "block_index"], kind="stable").drop(columns=["_p", "_s", "_e", "_c"])
    cols = ["policy", "correctness_basis", "stratum"] + lead + ["experiment", "coverage", "block_index", "block_id", "units",
            "retained_correct", "retained_incorrect", "filtered_correct", "filtered_incorrect", "retention",
            "precision_retained", "bad_call_detection_sensitivity", "catastrophic_retained", "catastrophic_filtered"]
    out = out[cols].reset_index(drop=True)
    out["block_index"] = out.block_index.where(out.block_index >= 0).astype("Int64")   # NA outside block strata
    for c in ("catastrophic_retained", "catastrophic_filtered"):
        out[c] = out[c].astype("Int64")   # counts; NA where the basis has no catastrophic flag
    return out


# ------------------------------------------------------------------------------------------------ mosaic switches
def _switch_class(ctx, frame):
    key = frame[ctx.switch_key] if ctx.switch_key in frame.columns else pd.Series("NA", index=frame.index)
    return key.map(lambda k: ctx.switch_plans[k]["switch_class"] if k in ctx.switch_plans else "")


def switch_class_order(ctx):
    """Classes by their fewest total switches, then name: a panel control (0) first."""
    best = {}
    for p in ctx.switch_plans.values():
        n = len(p["hap1"]) + len(p["hap2"])
        best[p["switch_class"]] = min(n, best.get(p["switch_class"], n))
    return sorted(best, key=lambda c: (best[c], c))


def switch_map(ctx):
    rows = []
    for key in sorted(ctx.switch_plans):
        p = ctx.switch_plans[key]
        for i in range(1, len(ctx.block_order)):
            h1, h2 = i in p["hap1"], i in p["hap2"]
            rows.append({"switch_class": p["switch_class"], "plan": key, "source_a": p["source_a"], "source_b": p["source_b"],
                         "boundary_index": i - 1, "boundary": f"{ctx.block_order[i - 1]}>{ctx.block_order[i]}",
                         "hap1_switch": h1, "hap2_switch": h2,
                         "category": "both" if h1 and h2 else "hap1_only" if h1 else "hap2_only" if h2 else "none"})
    out = pd.DataFrame(rows)
    order = {c: i for i, c in enumerate(switch_class_order(ctx))}
    out["_c"] = out.switch_class.map(order)
    return out.sort_values(["_c", "plan", "boundary_index"], kind="stable").drop(columns=["_c"]).reset_index(drop=True)


def mosaic_runs(ctx):
    r = ctx.runs.copy()
    r["switch_class"] = _switch_class(ctx, r)
    return r[r.switch_class != ""]


def mosaic_callability(ctx):
    r = mosaic_runs(ctx)
    r["no_call_reason"] = np.where(r.call_class.isin(CALLED), "none", r.no_call_reason)
    out = _count(r, ctx.strata + ["switch_class", "coverage", "call_class", "no_call_reason"], "runs", ctx,
                 ctx.strata + ["switch_class", "coverage"])
    order = {c: i for i, c in enumerate(switch_class_order(ctx))}
    out["_c"] = out.switch_class.map(order)
    return out.sort_values(["_c"], kind="stable").drop(columns=["_c"]).reset_index(drop=True)


def mosaic_excess_summary(ctx, points):
    m = points[points.switch_class.fillna("") != ""]
    rows = []
    for keys, g in m.groupby(ctx.strata + ["switch_class", "coverage"], sort=False):
        keys = keys if isinstance(keys, tuple) else (keys,)
        d = g[g.truth_defined & g.evaluated]
        row = {**dict(zip(ctx.strata + ["switch_class", "coverage"], keys)), "runs": g.run_id.nunique(), "blocks": len(g),
               "truth_undefined": int((g.evaluated & ~g.truth_defined).sum()), "denominator": len(d),
               "floor_attained": int((d.excess == 0).sum())}
        for metric in ("excess", "normalized_excess"):
            v = d[metric].dropna()
            row[metric + "_median"] = v.median() if len(v) else np.nan
            row[metric + "_mean"] = v.mean() if len(v) else np.nan
            row[metric + "_max"] = v.max() if len(v) else np.nan
        rows.append(row)
    out = pd.DataFrame(rows)
    order = {c: i for i, c in enumerate(switch_class_order(ctx))}
    cov = {c: i for i, c in enumerate(ctx.coverages)}
    out["_c"], out["_v"] = out.switch_class.map(order), out.coverage.map(cov)
    return out.sort_values(ctx.strata + ["_c", "_v"], kind="stable").drop(columns=["_c", "_v"]).reset_index(drop=True)


# ------------------------------------------------------------------------------------------------ dosage
def dosage_content_usable(ctx, points):
    c = points.merge(ctx.content, on=["run_id", "block_id"], how="left")
    sens = c[~c.dosage_filter.isin(["NOT_APPLICABLE", "NO_CALL"])].copy()
    rows = []
    for keys, g in sens.groupby(ctx.G, sort=False):
        keys = keys if isinstance(keys, tuple) else (keys,)
        ret = g.dosage_filter == "PASS"
        defined = g.content_correct.notna()
        ok = g.content_correct == True  # noqa: E712
        n_ret_def = int((ret & defined).sum())
        rows.append({**dict(zip(ctx.G, keys)), "dosage_sensitive_blocks": len(g), "retained": int(ret.sum()),
                     "retention": ret.mean() if len(g) else np.nan, "content_truth_defined": int(defined.sum()),
                     "content_truth_undefined": int((~defined).sum()), "retained_content_correct": int((ret & ok).sum()),
                     "retained_content_incorrect": int((ret & defined & ~ok).sum()),
                     "filtered_content_correct": int((~ret & ok).sum()), "filtered_content_incorrect": int((~ret & defined & ~ok).sum()),
                     "content_precision_retained": int((ret & ok).sum()) / n_ret_def if n_ret_def else np.nan})
    return sort_frame(pd.DataFrame(rows), ctx)


def dosage_content_vs_exact(ctx, points):
    c = points.merge(ctx.content, on=["run_id", "block_id"], how="inner")
    c = c[c.content_correct.notna() & c.truth_defined & c.evaluated].copy()
    c["dosage_content"] = np.where(c.content_correct == True, "content correct", "content incorrect")  # noqa: E712
    c["exact_diplotype"] = np.where(c.excess == 0, "floor attained", "floor not attained")
    c["dosage_filter_pass"] = c.dosage_filter == "PASS"
    return _count(c, ctx.strata + ["panel", "dosage_filter_pass", "dosage_content", "exact_diplotype"], "blocks", ctx,
                  ctx.strata + ["panel", "dosage_filter_pass"])


def representatives(ctx):
    r = ctx.representatives.merge(ctx.blocks[ctx.G + ["panel", "run_id", "block_id", "block_index"]], on=["run_id", "block_id"])
    return sort_frame(r, ctx, ["run_id"])


# ------------------------------------------------------------------------------------------------ phase
def _boundaries(ctx):
    b = ctx.phase_boundaries.merge(ctx.runs[ctx.G + ["run_id", "call_class", "phase_support"] +
                                           ([ctx.switch_key] if ctx.switch_key and ctx.switch_key in ctx.runs.columns and ctx.switch_key != "run_id" else [])],
                                   on="run_id", how="left")
    if ctx.switch_plans is not None:
        b["switch_class"] = _switch_class(ctx, b)
        key = b[ctx.switch_key] if ctx.switch_key in b.columns else b.run_id
        def cat(k, i):
            p = ctx.switch_plans.get(k)
            if p is None:
                return "no switch plan"
            h1, h2 = (i + 1) in p["hap1"], (i + 1) in p["hap2"]
            return "both" if h1 and h2 else "hap1_only" if h1 else "hap2_only" if h2 else "none"
        b["truth_switch_category"] = [cat(k, i) for k, i in zip(key, b.boundary_index)]
        if "table_truth_switch_category" in b.columns:
            plan = b.truth_switch_category != "no switch plan"
            clash = plan & (b.table_truth_switch_category != b.truth_switch_category)
            if clash.any():
                raise Refusal(f"the phase table's truth_switch_category disagrees with the switch metadata at {int(clash.sum())} "
                              "boundaries")
    else:
        b["switch_class"] = ""
        b["truth_switch_category"] = "no switch plan"
    return b


def phase_boundary_status(ctx):
    b = _boundaries(ctx)
    return _count(b, ctx.G + ["status", "reason"], "boundaries", ctx)


def _phase_rows(frame, keys, ctx):
    rows = []
    groups = frame.groupby(keys, sort=False) if keys else [((), frame)]
    for k, g in groups:
        k = k if isinstance(k, tuple) else (k,)
        a = g[g.status == "ASSESSED"]
        n, ok, err = len(a), int((a.phase_correct == True).sum()), int((a.switch_error == True).sum())  # noqa: E712
        lo, hi = wilson(err, n)
        rows.append({**dict(zip(keys, k)), "runs": g.run_id.nunique(), "boundaries": len(g), "assessable": n,
                     "correct": ok, "switch_errors": err, "switch_error_rate": err / n if n else np.nan,
                     "ci_low": lo, "ci_high": hi})
    return rows


def phase_accuracy(ctx):
    return sort_frame(pd.DataFrame(_phase_rows(_boundaries(ctx), ctx.G, ctx)), ctx)


def phase_by_switch_class(ctx):
    b = _boundaries(ctx)
    b = b[b.switch_class != ""]
    out = pd.DataFrame(_phase_rows(b, ctx.strata + ["switch_class", "coverage"], ctx))
    order = {c: i for i, c in enumerate(switch_class_order(ctx))}
    cov = {c: i for i, c in enumerate(ctx.coverages)}
    out["_c"], out["_v"] = out.switch_class.map(order), out.coverage.map(cov)
    return out.sort_values(ctx.strata + ["_c", "_v"], kind="stable").drop(columns=["_c", "_v"]).reset_index(drop=True)


def phase_by_boundary_type(ctx):
    b = _boundaries(ctx)
    out = pd.DataFrame(_phase_rows(b, ctx.strata + ["experiment", "truth_switch_category"], ctx))
    cat = {c: i for i, c in enumerate(["none", "hap1_only", "hap2_only", "both", "no switch plan"])}
    out["_t"] = out.truth_switch_category.map(cat)
    exp = {e: i for i, e in enumerate(ctx.experiments)}
    out["_e"] = out.experiment.map(exp)
    return out.sort_values(ctx.strata + ["_e", "_t"], kind="stable").drop(columns=["_e", "_t"]).reset_index(drop=True)


def phase_locus_status(ctx):
    c = ctx.phase_calls.merge(ctx.runs[ctx.G + ["run_id", "call_class"]], on="run_id", how="left")
    rows = []
    for keys, g in c.groupby(ctx.G, sort=False):
        keys = keys if isinstance(keys, tuple) else (keys,)
        strict = g[g.call_class == "strict_call"]
        rows.append({**dict(zip(ctx.G, keys)), "runs": len(g), "strict_calls": len(strict),
                     "genotype_equivalent_unphased": int((g.call_class == "genotype_equivalent_call").sum()),
                     "no_calls": int((g.call_class == "no_call").sum()),
                     "exact_phased_locus": int((g.exact_phased_locus == True).sum()),  # noqa: E712
                     "all_assessable_correct": int((g.all_assessable_correct == True).sum()),  # noqa: E712
                     "any_switch_error": int((g.any_switch_error == True).sum())})  # noqa: E712
    return sort_frame(pd.DataFrame(rows), ctx)


def phase_projected(ctx):
    p = ctx.phase_projected.merge(ctx.runs[ctx.G + ["run_id"]], on="run_id", how="left")
    rows = []
    for keys, g in p.groupby(ctx.G, sort=False):
        keys = keys if isinstance(keys, tuple) else (keys,)
        e = g[g.projected_status == "EVALUATED"]
        k, n = int((e.projected_concordant == True).sum()), len(e)  # noqa: E712
        lo, hi = wilson(k, n)
        rows.append({**dict(zip(ctx.G, keys)), "boundaries": len(g), "evaluated": n, "concordant": k,
                     "concordance": k / n if n else np.nan, "ci_low": lo, "ci_high": hi})
    return sort_frame(pd.DataFrame(rows), ctx)


def strict_phase_support(ctx):
    """For captions: the phase support of every strict call, per experiment."""
    s = ctx.runs[ctx.runs.call_class == "strict_call"]
    return _count(s, ctx.strata + ["experiment", "phase_support"], "strict_calls", ctx, ctx.strata + ["experiment"])


# ------------------------------------------------------------------------------------------------ GQ diagnostic
def gq_diagnostic(ctx, points):
    d = points[points.truth_defined & points.evaluated].copy()
    d["gq_bin"] = pd.cut(d.conditional_GQ_value, GQ_BINS, right=False, labels=GQ_LABELS).astype(object).fillna("GQ not reported")
    rows = []
    for keys, g in d.groupby(ctx.strata + ["panel", "gq_bin"], sort=False):
        keys = keys if isinstance(keys, tuple) else (keys,)
        k, n = int((g.excess == 0).sum()), len(g)
        lo, hi = wilson(k, n)
        rows.append({**dict(zip(ctx.strata + ["panel", "gq_bin"], keys)), "denominator": n, "floor_attained": k,
                     "rate": k / n if n else np.nan, "ci_low": lo, "ci_high": hi})
    return sort_frame(pd.DataFrame(rows), ctx)
