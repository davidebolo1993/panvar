#!/usr/bin/env python3
"""P1: compose ONE genotype profile from (the analysis's own base profile) + (declared operational overrides) + (the
development performance overlay, profiles/genotype_dev_performance_v1.tsv), so that a run records a single profile
path and SHA-256 that contains every setting.

Rules (P1.1: every violation raises ValueError, so `python -O` cannot strip a check):
  * mode is exactly one of MODES (full, reference, none);
  * the overlay may contain only the keys in OVERLAY_SECTIONS (backend choices, cache/workspace budgets, consumption
    guards, reporting); a scientific or model key there is refused;
  * an operational override is KEY=VALUE:REASON with a non-empty value and reason; KEY must be in
    OPERATIONAL_OVERRIDES (resource guards only: wall-clock and physical-footprint limits) and must replace a key the
    base already sets. Scientific/model settings (k-mer, copying prior, emission model, ...) and work caps are refused;
  * no key may come from two sources except a declared override;
  * `reference` mode keeps only REFERENCE_ROUTE: the caches and consumption caps that exist in the production caller,
    with the scoring backends at their defaults (the I-series' "c4v6" reference arm).
The composed file lists its inputs and their SHA-256 in comments.

  compose_profile.py OUT NAME BASE_PROFILE {full|reference|none} [KEY=VALUE:REASON ...]
  compose_profile.py --check COMPOSED_PROFILE ...      (audit already composed profiles against the policy)
"""
import hashlib
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
OVERLAY = os.path.abspath(os.path.join(HERE, "..", "..", "..", "profiles", "genotype_dev_performance_v1.tsv"))
OVERLAY_SECTIONS = {
    "performance backend choices": ["union_scoring_backend", "adaptive_window_discovery", "indexed_union_store_layout",
                                    "indexed_union_evidence_layout", "adaptive_fused_seed",
                                    "marker_proposal_reuse_table_bytes", "combination_prepared_library_terms",
                                    "gap_recursion_backend", "dosage_census_backend"],
    "cache and workspace budgets": ["indexed_union_store_bytes", "adaptive_fine_store_bytes", "emission_cache_bytes",
                                    "repair_origin_cache_bytes", "union_placement_cache_bytes",
                                    "union_origin_cache_bytes", "union_prefilter_cache_bytes",
                                    "start_sort_scratch_bytes"],
    "resource capacities and guards": ["adaptive_fine_max_served_placements", "max_union_cached_placements_consumed",
                                       "max_union_cached_origins_consumed", "max_union_cached_states_consumed",
                                       "max_union_cached_prefilters_consumed"],
    "reporting": ["timing_detail"],
}
MODES = ("full", "reference", "none")
# P1.1: the only keys an operational override may change (resource guards that stop a run; they never change a
# completed result). Everything else -- scientific/model settings and work caps -- comes from the base profile as is.
OPERATIONAL_OVERRIDES = {
    "time_limit_seconds": "wall-clock guard (the watchdog records TIME_LIMIT)",
    "footprint_limit_bytes": "physical-footprint guard (the watchdog records RESOURCE_REFUSAL)",
}
REFERENCE_ROUTE = ["emission_cache_bytes", "repair_origin_cache_bytes", "union_placement_cache_bytes",
                   "max_union_cached_placements_consumed", "union_origin_cache_bytes",
                   "max_union_cached_origins_consumed", "max_union_cached_states_consumed",
                   "union_prefilter_cache_bytes", "max_union_cached_prefilters_consumed", "start_sort_scratch_bytes"]


def sha256(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def refuse(message):
    raise ValueError(message)


def read_profile(path):
    rows = []
    for line in open(path).read().split("\n"):
        line = line.rstrip("\r")
        if not line or line.startswith("#"):
            continue
        if "\t" not in line:
            refuse(f"{path}: expected key<TAB>value: {line}")
        key, value = line.split("\t", 1)
        rows.append((key, value))
    keys = [k for k, _ in rows]
    if len(keys) != len(set(keys)):
        refuse(f"{path}: a key is set twice")
    return rows


def parse_override(item):
    """KEY=VALUE:REASON -> (key, value, reason); refuses anything else."""
    if "=" not in item or ":" not in item.split("=", 1)[1]:
        refuse(f"override must be KEY=VALUE:REASON: {item!r}")
    key, rest = item.split("=", 1)
    value, reason = rest.split(":", 1)
    if not key or not value.strip() or not reason.strip():
        refuse(f"override needs a key, a value and a reason: {item!r}")
    if key not in OPERATIONAL_OVERRIDES:
        refuse(f"override of {key!r} refused: only operational guards {sorted(OPERATIONAL_OVERRIDES)} may be overridden "
               "(scientific/model settings and work caps come from the base profile)")
    return key, value, reason


def overlay_rows():
    rows = dict(read_profile(OVERLAY))
    if rows.pop("profile", None) != "genotype-dev-performance-v1":
        refuse("the overlay must name itself genotype-dev-performance-v1")
    allowed = [k for keys in OVERLAY_SECTIONS.values() for k in keys]
    unknown = set(rows) - set(allowed)
    if unknown:
        refuse(f"overlay holds keys outside its declared sections (scientific keys are refused): {unknown}")
    if set(rows) != set(allowed):
        refuse(f"overlay is missing declared keys: {set(allowed) - set(rows)}")
    return rows


def compose(out, name, base, mode, overrides):
    if mode not in MODES:
        refuse(f"mode must be one of {MODES}: {mode!r}")
    base_rows = [(k, v) for k, v in read_profile(base) if k != "profile"]
    base_keys = {k for k, _ in base_rows}
    over = {}
    for item in overrides:
        key, value, reason = parse_override(item)
        if key not in base_keys:
            refuse(f"override {key} must replace a key the base profile sets")
        if key in over:
            refuse(f"override {key} twice")
        over[key] = (value, reason)
    lines = ["# composed by tests/analysis/p1/compose_profile.py (P1); do not edit",
             f"# base profile (analysis settings): {base} sha256 {sha256(base)}"]
    if mode != "none":
        lines.append(f"# overlay ({mode}): {OVERLAY} sha256 {sha256(OVERLAY)}")
    lines += [f"profile\t{name}", "## analysis base profile (scientific/model settings and its own guards)"]
    lines += [f"{k}\t{v}" for k, v in base_rows if k not in over]
    if over:
        lines.append("## operational overrides (declared; resource guards only)")
        for k, (v, reason) in over.items():
            lines += [f"# {k}: {reason}", f"{k}\t{v}"]
    if mode != "none":
        rows = overlay_rows()
        for section, keys in OVERLAY_SECTIONS.items():
            chosen = [k for k in keys if mode == "full" or k in REFERENCE_ROUTE]
            if not chosen:
                continue
            lines.append(f"## overlay: {section}")
            for k in chosen:
                if k in base_keys:
                    refuse(f"{k} is set by both the base profile and the overlay")
                lines.append(f"{k}\t{rows[k]}")
    text = "\n".join(lines) + "\n"
    open(out, "w").write(text)
    read_profile(out)   # duplicate check on the result
    check_composed(out)
    return out


def check_composed(path):
    """Audits a composed profile: every row in its operational-override section is an allowed guard, preceded by its
    reason; no key repeats. Returns the overridden keys."""
    read_profile(path)
    keys, section, reason_seen = [], None, False
    for line in open(path).read().split("\n"):
        if line.startswith("## "):
            section = line
            continue
        if section and section.startswith("## operational overrides"):
            if line.startswith("# ") and ":" in line:
                reason_seen = bool(line.split(":", 1)[1].strip())
                continue
            if "\t" in line:
                key = line.split("\t", 1)[0]
                if key not in OPERATIONAL_OVERRIDES:
                    refuse(f"{path}: override of {key!r} is not an operational guard")
                if not reason_seen:
                    refuse(f"{path}: override of {key!r} has no reason")
                keys.append(key)
                reason_seen = False
    return keys


if __name__ == "__main__":
    if sys.argv[1] == "--check":
        for p in sys.argv[2:]:
            print(p, check_composed(p))
    else:
        compose(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5:])
