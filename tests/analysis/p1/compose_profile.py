#!/usr/bin/env python3
"""P1: compose ONE genotype profile from (the analysis's own base profile) + (declared operational overrides) + (the
development performance overlay, profiles/genotype_dev_performance_v1.tsv), so that a run records a single profile
path and SHA-256 that contains every setting.

Rules:
  * the overlay may contain only the keys in OVERLAY_SECTIONS (backend choices, cache/workspace budgets, consumption
    guards, reporting); a scientific or model key there is refused;
  * no key may come from two sources, except an operational override, which must replace a key the base already sets
    (it is declared with its reason and written in its own section);
  * `reference` mode keeps only REFERENCE_ROUTE: the caches and consumption caps that exist in the production caller,
    with the scoring backends at their defaults (the I-series' "c4v6" reference arm).
The composed file lists its inputs and their SHA-256 in comments.

  compose_profile.py OUT NAME BASE_PROFILE {full|reference|none} [KEY=VALUE:REASON ...]
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
REFERENCE_ROUTE = ["emission_cache_bytes", "repair_origin_cache_bytes", "union_placement_cache_bytes",
                   "max_union_cached_placements_consumed", "union_origin_cache_bytes",
                   "max_union_cached_origins_consumed", "max_union_cached_states_consumed",
                   "union_prefilter_cache_bytes", "max_union_cached_prefilters_consumed", "start_sort_scratch_bytes"]


def sha256(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def read_profile(path):
    rows = []
    for line in open(path).read().split("\n"):
        line = line.rstrip("\r")
        if not line or line.startswith("#"):
            continue
        key, value = line.split("\t", 1)
        rows.append((key, value))
    keys = [k for k, _ in rows]
    assert len(keys) == len(set(keys)), f"{path}: a key is set twice"
    return rows


def overlay_rows():
    rows = dict(read_profile(OVERLAY))
    assert rows.pop("profile") == "genotype-dev-performance-v1"
    allowed = [k for keys in OVERLAY_SECTIONS.values() for k in keys]
    unknown = set(rows) - set(allowed)
    assert not unknown, f"overlay holds keys outside its declared sections (scientific keys are refused): {unknown}"
    assert set(rows) == set(allowed), f"overlay is missing declared keys: {set(allowed) - set(rows)}"
    return rows


def compose(out, name, base, mode, overrides):
    base_rows = [(k, v) for k, v in read_profile(base) if k != "profile"]
    base_keys = {k for k, _ in base_rows}
    over = {}
    for item in overrides:
        kv, reason = item.split(":", 1)
        key, value = kv.split("=", 1)
        assert key in base_keys, f"override {key} must replace a key the base profile sets"
        assert key not in over, f"override {key} twice"
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
                assert k not in base_keys, f"{k} is set by both the base profile and the overlay"
                lines.append(f"{k}\t{rows[k]}")
    open(out, "w").write("\n".join(lines) + "\n")
    read_profile(out)   # duplicate check on the result
    return out


if __name__ == "__main__":
    compose(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5:])
