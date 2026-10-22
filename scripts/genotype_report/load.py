"""Loading and schema validation. Two input forms normalise to the same two tables:

  runs    one row per genotyper run (run_id, call class, headline locus filters, provenance, resources)
  blocks  one row per reported block of a called run, in model order (block_index)

and every optional input (run metadata, evaluation, dosage content, switch plans, phase tables) is validated and keyed
by run_id before anything is aggregated. Nothing here writes a file.
"""
import hashlib
import os

import pandas as pd

from .common import (CALLED, NORMALIZATION_COLUMNS, Refusal, attach_block_id, attach_run_id, read_tsv, require_columns,
                     to_bool, to_number)

RUN_COLUMNS = ["run_id", "run_status", "call_class", "no_call_reason", "phase_support", "dosage_mode", "dosage_audit_status",
               "locus_repair_scope", "genotype_filter", "genotype_filter_reasons", "dosage_filter", "dosage_filter_reasons",
               "phase_filter", "phase_filter_reasons", "call_summary", "seconds_total", "footprint_peak_bytes",
               "repair_rounds_run", "seconds_dosage_total", "threads_requested", "profile_name", "profile_sha256", "preset",
               "preset_source", "settings_digest", "graph_sha256", "library_mode", "library_sensitivity_status"]
BLOCK_COLUMNS = ["run_id", "block_index", "block_id", "kind", "zygosity", "genotype_support", "conditional_GQ",
                 "representation_filter", "quality_filter", "dosage_status", "repair_class", "dosage_repair_applied",
                 "dosage_content_scope", "genotype_filter", "genotype_filter_reasons", "dosage_filter", "dosage_filter_reasons"]
HEADLINE_BLOCK = ["genotype_filter", "genotype_filter_reasons", "dosage_filter", "dosage_filter_reasons"]
HEADLINE_LOCUS = ["genotype_filter", "genotype_filter_reasons", "dosage_filter", "dosage_filter_reasons", "phase_filter",
                  "phase_filter_reasons", "call_summary"]
METADATA_KNOWN = ["experiment", "coverage", "subject", "modality", "execution_group"]


def _kv(path):
    out = {}
    if os.path.isfile(path):
        with open(path) as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 2 and parts[0] not in out:
                    out[parts[0]] = parts[1]
    return out


def _settings_digest(prefix, status):
    """A digest of every effective setting except the worker count, from <prefix>.settings.tsv; for outputs that
    predate it, the profile and preset identities the status records."""
    path = prefix + ".settings.tsv"
    if os.path.isfile(path):
        t = read_tsv(path, "settings table")
        require_columns(t, ["key", "value"], path)
        text = "".join(f"{k}\t{v}\n" for k, v in zip(t.key, t.value) if k != "threads")
        return "settings:" + hashlib.sha256(text.encode()).hexdigest()[:16]
    return "profile:" + status.get("profile_sha256", "NA")[:16] + "|preset:" + status.get("preset", "none")


def primary_prefix(prefix):
    """M15-L: a run's diagnostics live in <P>.diagnostics/<name>.*; the run is identified by its primary prefix P. A flat
    prefix (runs before M15-L, or the internal API) is its own primary prefix."""
    d, name = os.path.split(prefix)
    if os.path.basename(d) == name + ".diagnostics":
        return os.path.join(os.path.dirname(d), name)
    return prefix


def _read_prefix(root, prefix, sources):
    run_id = os.path.relpath(primary_prefix(prefix), root).replace(os.sep, "/")
    L = read_tsv(prefix + ".locus_genotype.tsv", "locus table")
    if len(L) != 1:
        raise Refusal(f"{prefix}.locus_genotype.tsv: expected one row, found {len(L)}")
    require_columns(L, ["run_status", "call_class"], prefix + ".locus_genotype.tsv")
    loc = L.iloc[0]
    sources.append(prefix + ".locus_genotype.tsv")
    status = _kv(prefix + ".mosaic_status.tsv")
    if status:
        sources.append(prefix + ".mosaic_status.tsv")
    dstatus = _kv(prefix + ".dosage_status.tsv")
    if dstatus:
        sources.append(prefix + ".dosage_status.tsv")
    lib = _kv(prefix + ".library_sensitivity.tsv")
    if os.path.isfile(prefix + ".settings.tsv"):
        sources.append(prefix + ".settings.tsv")
    provenance = status.get("library_provenance", "")
    library_mode = lib.get("library_mode") or status.get("library_mode") or ("file" if provenance.startswith("declared") else "NA")
    if library_mode == "auto" or status.get("library_mode") == "auto":
        library_mode = "auto"
    run = {"run_id": run_id, "run_status": loc.run_status, "call_class": loc.call_class,
           "no_call_reason": loc.get("no_call_reason", "NA"), "phase_support": loc.get("phase_support", status.get("phase_support", "NA")),
           "dosage_mode": loc.get("dosage_mode", status.get("dosage_mode", "off")),
           "dosage_audit_status": loc.get("dosage_audit_status", "NA"), "locus_repair_scope": loc.get("locus_repair_scope", "NA"),
           "seconds_total": status.get("seconds_total", "NA"),
           "footprint_peak_bytes": status.get("resources_phys_footprint_peak_bytes", "NA"),
           "repair_rounds_run": status.get("repair_rounds_run", "NA"),
           "seconds_dosage_total": dstatus.get("seconds_dosage_total", "NA"),
           "threads_requested": status.get("threads_requested", "NA"), "profile_name": status.get("profile_name", "NA"),
           "profile_sha256": status.get("profile_sha256", "NA"), "preset": status.get("preset", "none"),
           "preset_source": status.get("preset_source", "NA"), "settings_digest": _settings_digest(prefix, status),
           "graph_sha256": status.get("sha256_gfa", "NA"), "library_mode": library_mode,
           "library_sensitivity_status": lib.get("library_sensitivity_status", status.get("library_sensitivity_status", "NA"))}
    has_headline = "genotype_filter" in L.columns
    for c in HEADLINE_LOCUS:
        run[c] = loc[c] if has_headline else "NA"
    blocks = []
    table = prefix + ".block_genotypes.tsv"
    if os.path.isfile(table):
        sources.append(table)
        B = read_tsv(table, "block table")
        require_columns(B, ["block_id"], table)
        if "genotype_filter" not in B.columns and "dosage_filter" in B.columns:
            # outputs before headline_filters_v2: `dosage_filter` is the dosage stage's component, since renamed
            B = B.rename(columns={"dosage_filter": "dosage_stage_filter", "dosage_filter_reason": "dosage_stage_filter_reason"})
        dosage = {}
        if os.path.isfile(prefix + ".dosage_blocks.tsv"):
            sources.append(prefix + ".dosage_blocks.tsv")
            D = read_tsv(prefix + ".dosage_blocks.tsv", "dosage block table")
            dosage = {r.block_id: r for r in D.itertuples(index=False)}
        for i, r in enumerate(B.to_dict("records")):
            d = dosage.get(r["block_id"])
            pick = lambda col: r[col] if col in r and r[col] != "" else (getattr(d, col, "NA") if d is not None else "NA")
            blocks.append({"run_id": run_id, "block_index": i, "block_id": r["block_id"], "kind": r.get("kind", "NA"),
                           "zygosity": r.get("zygosity", "NA"), "genotype_support": r.get("genotype_support", "NA"),
                           "conditional_GQ": r.get("conditional_GQ", "NA"), "representation_filter": r.get("filter", "NA"),
                           "quality_filter": r.get("quality_filter", "NA"), "dosage_status": pick("dosage_status"),
                           "repair_class": pick("repair_class"), "dosage_repair_applied": r.get("dosage_repair_applied", "NA"),
                           "dosage_content_scope": pick("dosage_content_scope"),
                           **{c: (r[c] if c in r else "NA") for c in HEADLINE_BLOCK}})
    elif loc.call_class in CALLED:
        raise Refusal(f"{prefix}: a {loc.call_class} run without .block_genotypes.tsv")
    return run, blocks


def load_directory(root):
    """Every <prefix>.locus_genotype.tsv under root (recursively; <prefix>.original.* dosage-repair originals are not
    runs). run_id is the run's primary prefix path relative to root (M15-L: the locus table is in <P>.diagnostics/, and the
    run is P)."""
    root = os.path.abspath(root)
    prefixes = []
    for d, _, files in os.walk(root):
        for f in files:
            if f.endswith(".locus_genotype.tsv"):
                p = os.path.join(d, f[: -len(".locus_genotype.tsv")])
                if not p.endswith(".original"):
                    prefixes.append(p)
    if not prefixes:
        raise Refusal(f"no panvar genotype output prefix (<prefix>.locus_genotype.tsv) under {root}")
    runs, blocks, sources = [], [], []
    for p in sorted(prefixes):
        r, b = _read_prefix(root, p, sources)
        runs.append(r)
        blocks.extend(b)
    return pd.DataFrame(runs, columns=RUN_COLUMNS), pd.DataFrame(blocks, columns=BLOCK_COLUMNS), sources


def load_package(root):
    """A consolidated package: runs.tsv and blocks.tsv in the normalised schema (run_id or dataset/job keys). Optional
    columns absent from the package read as NA."""
    rp, bp = os.path.join(root, "runs.tsv"), os.path.join(root, "blocks.tsv")
    runs = attach_run_id(read_tsv(rp, "package runs.tsv"), rp)
    blocks = attach_block_id(attach_run_id(read_tsv(bp, "package blocks.tsv"), bp), bp)
    require_columns(runs, ["run_id", "run_status", "call_class"], rp)
    if "dosage_filter" in blocks.columns and "genotype_filter" not in blocks.columns:
        raise Refusal(f"{bp}: a dosage_filter column without genotype_filter is the pre-headline_filters_v2 dosage "
                      "component (now dosage_stage_filter); rename it or supply the headline columns")
    if "block_index" not in blocks.columns:
        blocks["block_index"] = blocks.groupby("run_id", sort=False).cumcount().astype(str)
    for c in RUN_COLUMNS:
        if c not in runs.columns:
            runs[c] = "NA"
    for c in BLOCK_COLUMNS:
        if c not in blocks.columns:
            blocks[c] = "NA"
    extra_run = [c for c in runs.columns if c not in RUN_COLUMNS]
    return runs[RUN_COLUMNS + extra_run], blocks[BLOCK_COLUMNS], [rp, bp]


def validate(runs, blocks):
    """Refuse duplicate runs, blocks of unknown or uncalled runs, and runs whose block sequences disagree: one report
    covers one locus model."""
    if runs.run_id.duplicated().any():
        raise Refusal("duplicate run_id: " + ", ".join(runs.run_id[runs.run_id.duplicated()].unique()[:5]))
    unknown = sorted(set(blocks.run_id) - set(runs.run_id))
    if unknown:
        raise Refusal("blocks for unknown runs: " + ", ".join(unknown[:5]))
    blocks = blocks.copy()
    blocks["block_index"] = pd.to_numeric(blocks.block_index, errors="coerce")
    if blocks.block_index.isna().any():
        raise Refusal("block_index is not an integer")
    blocks = blocks.sort_values(["run_id", "block_index"], kind="stable")
    if blocks.duplicated(["run_id", "block_id"]).any():
        raise Refusal("a block_id appears twice in one run")
    seqs = blocks.groupby("run_id", sort=True).block_id.apply(tuple)
    distinct = seqs.value_counts()
    if len(distinct) > 1:
        a, b = distinct.index[0], distinct.index[1]
        raise Refusal(f"runs disagree on the block model ({len(distinct)} distinct block sequences, e.g. {','.join(a)} vs "
                      f"{','.join(b)}): one report covers one locus model")
    order = list(distinct.index[0]) if len(distinct) else []
    called = set(runs.run_id[runs.call_class.isin(CALLED)])
    blockless = sorted(called - set(blocks.run_id))
    if blockless:
        raise Refusal("called runs without blocks: " + ", ".join(blockless[:5]))
    return blocks.reset_index(drop=True), order


def load_run_metadata(path, runs):
    """Run metadata: a run key and any of experiment, coverage, subject, modality, execution_group (other columns are
    carried). It never overrides a value the runs already record differently. Without it every run is in one group."""
    runs = runs.copy()
    if path is None:
        for c in ("experiment", "coverage"):
            if c not in runs.columns or (runs[c] == "NA").all():
                runs[c] = "all"
        return runs, None
    meta = attach_run_id(read_tsv(path, "run metadata"), path)
    if meta.run_id.duplicated().any():
        raise Refusal(f"{path}: duplicate run keys")
    missing = sorted(set(runs.run_id) - set(meta.run_id))
    if missing:
        raise Refusal(f"{path}: no metadata for {len(missing)} run(s), e.g. {', '.join(missing[:3])}")
    extra = sorted(set(meta.run_id) - set(runs.run_id))
    if extra:
        raise Refusal(f"{path}: metadata for {len(extra)} run(s) not in the input, e.g. {', '.join(extra[:3])}")
    if "modality" in meta.columns and (meta.modality.isin(["", "NA"])).any():
        raise Refusal(f"{path}: an empty modality; give every run its modality or drop the column (none is assumed)")
    keep = [c for c in meta.columns if c not in ("dataset", "job", "cell")]
    merged = runs.merge(meta[keep], on="run_id", how="left", suffixes=("", "__meta"))
    for c in keep:
        if c + "__meta" in merged.columns:
            clash = (merged[c] != "NA") & (merged[c] != merged[c + "__meta"])
            if clash.any():
                raise Refusal(f"{path}: {c} disagrees with the runs for {int(clash.sum())} run(s)")
            merged[c] = merged[c + "__meta"]
            merged = merged.drop(columns=[c + "__meta"])
    for c in ("experiment", "coverage"):
        if c not in merged.columns:
            merged[c] = "all"
    order = list(meta.run_id)
    merged["_meta_order"] = merged.run_id.map({r: i for i, r in enumerate(order)})
    return merged, meta


def load_evaluation(path, runs, blocks):
    """The evaluator's per-block fields: a consolidated TSV (run key, block_id, unordered_distance, availability_floor,
    excess and, for diploid-length-normalised excess, truth1_bp, truth2_bp, called1_bp, called2_bp), or a manifest (run key, path) of
    per-run evaluator block tables. Every row must match a reported block."""
    t = read_tsv(path, "evaluation blocks")
    sources = [path]
    if "path" in t.columns and "unordered_distance" not in t.columns:
        t = attach_run_id(t, path)
        parts = []
        base = os.path.dirname(os.path.abspath(path))
        for rid, p in zip(t.run_id, t.path):
            full = p if os.path.isabs(p) else os.path.join(base, p)
            e = read_tsv(full, "evaluator block table")
            e.insert(0, "run_id", rid)
            parts.append(e)
            sources.append(full)
        t = pd.concat(parts, ignore_index=True) if parts else t
    else:
        t = attach_run_id(t, path)
    t = attach_block_id(t, path)
    require_columns(t, ["run_id", "block_id", "unordered_distance", "availability_floor", "excess"], path)
    if t.duplicated(["run_id", "block_id"]).any():
        raise Refusal(f"{path}: an evaluated block appears twice")
    unknown = sorted(set(t.run_id) - set(runs.run_id))
    if unknown:
        raise Refusal(f"{path}: evaluation for runs not in the input, e.g. {', '.join(unknown[:3])}")
    uncalled = sorted(set(t.run_id) & set(runs.run_id[~runs.call_class.isin(CALLED)]))
    if uncalled:
        raise Refusal(f"{path}: evaluation rows for runs that emitted no genotype, e.g. {', '.join(uncalled[:3])}")
    known = set(zip(blocks.run_id, blocks.block_id))
    stray = [k for k in zip(t.run_id, t.block_id) if k not in known]
    if stray:
        raise Refusal(f"{path}: {len(stray)} evaluated block(s) the runs did not report, e.g. {stray[0]}")
    has_norm = all(c in t.columns for c in NORMALIZATION_COLUMNS)
    out = pd.DataFrame({"run_id": t.run_id, "block_id": t.block_id,
                        "distance": to_number(t.unordered_distance), "floor": to_number(t.availability_floor),
                        "excess": to_number(t.excess)})
    if has_norm:
        for c in NORMALIZATION_COLUMNS:
            out[c] = to_number(t[c])
        out["truth_defined"] = out.truth1_bp.notna() & out.truth2_bp.notna()
        cols = pd.concat([out.truth1_bp + out.truth2_bp, out.called1_bp + out.called2_bp], axis=1).max(axis=1)
        out["diploid_length_denominator_bp"] = cols.where(out.truth_defined)
        out["normalized_excess"] = (out.excess / out.diploid_length_denominator_bp).where(out.truth_defined & (cols > 0))
    else:
        out["truth_defined"] = out.distance.notna()
        for c in ("truth1_bp", "truth2_bp"):   # the observed truth-allele lengths (block_lengths.tsv) need only these
            if c in t.columns:
                out[c] = to_number(t[c])
    # a truth-undefined block has no distance: it is neither correct nor incorrect
    for c in ("distance", "floor", "excess"):
        out[c] = out[c].where(out.truth_defined)
    return out, has_norm, sources


def load_dosage_content(path, blocks):
    t = attach_block_id(attach_run_id(read_tsv(path, "dosage content evaluation"), path), path)
    require_columns(t, ["run_id", "block_id", "content_correct"], path)
    if t.duplicated(["run_id", "block_id"]).any():
        raise Refusal(f"{path}: a block appears twice")
    known = set(zip(blocks.run_id, blocks.block_id))
    stray = [k for k in zip(t.run_id, t.block_id) if k not in known]
    if stray:
        raise Refusal(f"{path}: {len(stray)} block(s) the runs did not report, e.g. {stray[0]}")
    out = t[["run_id", "block_id"]].copy()
    out["content_correct"] = to_bool(t.content_correct)
    out["content_truth_source"] = t["content_truth_source"] if "content_truth_source" in t.columns else "NA"
    return out


def load_representatives(path, blocks):
    t = attach_block_id(attach_run_id(read_tsv(path, "representative diagnostics"), path), path)
    require_columns(t, ["run_id", "block_id", "original_distance", "representative_phased_distance"], path)
    known = set(zip(blocks.run_id, blocks.block_id))
    stray = [k for k in zip(t.run_id, t.block_id) if k not in known]
    if stray:
        raise Refusal(f"{path}: {len(stray)} block(s) the runs did not report, e.g. {stray[0]}")
    out = t[["run_id", "block_id"]].copy()
    for c in ("original_distance", "representative_phased_distance"):
        out[c] = to_number(t[c])
    for c in ("repair_class", "role", "truth_class"):
        out[c] = t[c] if c in t.columns else "NA"
    return out


def _indices(text):
    text = text.strip()
    if text in ("", "none", "NA"):
        return ()
    try:
        return tuple(sorted(int(x) for x in text.split(",") if x.strip()))
    except ValueError:
        raise Refusal(f"switch positions must be block indices: {text!r}")


def load_switch_metadata(path, runs, n_blocks):
    """Mosaic switch plans: subject (joined to the run metadata's subject) or a run key; switch_class;
    switch_before_block_hap1 and switch_before_block_hap2, comma-separated 0-based block indices i where that
    homologue switches source template between block i-1 and block i ('none' for no switch); optional source_a,
    source_b (the plan's two source templates). Switches are never inferred from calls."""
    t = read_tsv(path, "switch metadata")
    require_columns(t, ["switch_class", "switch_before_block_hap1", "switch_before_block_hap2"], path)
    key = "subject" if "subject" in t.columns else None
    if key is None:
        t = attach_run_id(t, path)
        key = "run_id"
    elif "subject" not in runs.columns:
        raise Refusal(f"{path} is keyed by subject but the run metadata has no subject column")
    if t[key].duplicated().any():
        raise Refusal(f"{path}: duplicate {key}")
    plans = {}
    for r in t.to_dict("records"):
        h1, h2 = _indices(r["switch_before_block_hap1"]), _indices(r["switch_before_block_hap2"])
        for i in h1 + h2:
            if not 1 <= i < n_blocks:
                raise Refusal(f"{path}: switch before block {i} is outside the {n_blocks}-block model")
        plans[r[key]] = {"switch_class": r["switch_class"], "hap1": h1, "hap2": h2,
                         "source_a": r.get("source_a", "NA"), "source_b": r.get("source_b", "NA")}
    return plans, key


def load_phase_boundaries(path, runs):
    t = attach_run_id(read_tsv(path, "phase boundaries"), path)
    require_columns(t, ["run_id", "boundary_index", "boundary", "status", "phase_correct", "switch_error"], path)
    unknown = sorted(set(t.run_id) - set(runs.run_id))
    if unknown:
        raise Refusal(f"{path}: boundaries for runs not in the input, e.g. {', '.join(unknown[:3])}")
    if t.duplicated(["run_id", "boundary_index"]).any():
        raise Refusal(f"{path}: a boundary appears twice")
    out = t[["run_id", "boundary_index", "boundary", "status"]].copy()
    out["boundary_index"] = pd.to_numeric(out.boundary_index)
    out["reason"] = t["reason"] if "reason" in t.columns else "NA"
    out["phase_correct"] = to_bool(t.phase_correct)
    out["switch_error"] = to_bool(t.switch_error)
    if "truth_switch_category" in t.columns:
        out["table_truth_switch_category"] = t.truth_switch_category
    return out


def load_phase_calls(path, runs):
    t = attach_run_id(read_tsv(path, "phase calls"), path)
    require_columns(t, ["run_id", "exact_phased_locus", "assessed_boundaries", "correct_boundaries", "switch_errors",
                        "all_assessable_correct", "any_switch_error"], path)
    unknown = sorted(set(t.run_id) - set(runs.run_id))
    if unknown:
        raise Refusal(f"{path}: phase calls for runs not in the input, e.g. {', '.join(unknown[:3])}")
    out = t[["run_id"]].copy()
    for c in ("exact_phased_locus", "all_assessable_correct", "any_switch_error"):
        out[c] = to_bool(t[c])
    for c in ("assessed_boundaries", "correct_boundaries", "switch_errors"):
        out[c] = pd.to_numeric(t[c])
    return out


def load_phase_projected(path, runs):
    t = attach_run_id(read_tsv(path, "sequence-projected phase"), path)
    require_columns(t, ["run_id", "boundary_index", "projected_status", "projected_concordant"], path)
    unknown = sorted(set(t.run_id) - set(runs.run_id))
    if unknown:
        raise Refusal(f"{path}: rows for runs not in the input, e.g. {', '.join(unknown[:3])}")
    out = t[["run_id", "boundary_index", "projected_status"]].copy()
    out["projected_concordant"] = to_bool(t.projected_concordant)
    return out
