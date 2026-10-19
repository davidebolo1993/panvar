#!/usr/bin/env python3
"""S2-L: whole-process summary of every control/combined run (saved outputs only).
  s2l_summary.py <runs dir> <out.tsv>
Columns: outcome and first binding limit (refusals: the status note); prediction, sequence status, completion and phase
support; elapsed total (status) and wall (/usr/bin/time); union (seconds_rescoring), refinement (seconds_repair) and
post-search (seconds_detail_adaptive_post_search) seconds; process peak footprint (/usr/bin/time, else the status);
rounds, moves, tie branching, library-sensitivity status; union executed/served origins and the union origin-cap margin;
refinement executed origins and margin; round-1 admission counts; identities (model fingerprint, input SHA-256 digest,
effective-settings digest over setting VALUES, prepared locus).
"""
import hashlib
import os
import re
import sys

COLUMNS = ["unit", "status", "first_binding_limit", "prediction", "sequence_status", "completion", "phase_support",
           "seconds_total", "wall_seconds", "union_seconds", "refinement_seconds", "post_search_seconds", "peak_footprint_bytes",
           "rounds", "moves", "tie_branched", "library_sensitivity", "union_executed_origins", "union_served_origins",
           "union_origin_cap_margin", "refinement_executed_origins", "refinement_origin_margin", "admission_singletons_not_retained",
           "admission_retained", "model_fingerprint", "inputs_sha256_digest", "settings_values_sha256", "prepared_locus"]


def kv(path):
    out = {}
    if os.path.exists(path):
        for line in open(path, errors="replace"):
            k, _, v = line.rstrip("\n").partition("\t")
            out[k] = v
    return out


def main(runs, out):
    rows = []
    for unit in sorted(os.listdir(runs)):
        if not (unit.startswith("control.") or unit.startswith("combined.")):
            continue
        d = os.path.join(runs, unit)
        st, locus = kv(f"{d}/call.mosaic_status.tsv"), kv(f"{d}/call.sequence_locus.tsv")
        log = open(f"{d}/run.log", errors="replace").read() if os.path.exists(f"{d}/run.log") else ""
        wall = re.search(r"([\d.]+) real", log)
        peak = re.search(r"(\d+)\s+peak memory footprint", log)
        settings = kv(f"{d}/call.settings.tsv")
        values = "\n".join(f"{k}={v.split(chr(9))[0]}" for k, v in sorted(settings.items()) if k != "key")
        inputs = "\n".join(f"{k}={v}" for k, v in sorted(st.items()) if k.startswith("sha256_"))
        refused = st.get("status") not in ("COMPLETE", "CONDITIONAL_GENOTYPE_EQUIVALENT", "BLOCK_CONSENSUS")

        def num(k):
            try:
                return int(st[k])
            except (KeyError, ValueError):
                return None
        cap, used = num("max_physical_origins"), num("physical_origins")
        row = {
            "unit": unit, "status": st.get("status", "NA"),
            "first_binding_limit": st.get("note", "NA") if refused else "none",
            "prediction": locus.get("prediction", "NA"), "sequence_status": locus.get("sequence_status", "NA"),
            "completion": locus.get("completion", "NA"), "phase_support": st.get("phase_support", "NA"),
            "seconds_total": st.get("seconds_total", st.get("elapsed_seconds", "NA")), "wall_seconds": wall.group(1) if wall else "NA",
            "union_seconds": st.get("seconds_rescoring", "NA"), "refinement_seconds": st.get("seconds_repair", "NA"),
            "post_search_seconds": st.get("seconds_detail_adaptive_post_search", "NA"),
            "peak_footprint_bytes": peak.group(1) if peak else st.get("resources_phys_footprint_peak_bytes", st.get("phys_footprint_peak_bytes", "NA")),
            "rounds": st.get("repair_rounds_run", "NA"), "moves": st.get("repair_moves", "NA"), "tie_branched": st.get("tie_branched", "NA"),
            "library_sensitivity": st.get("library_sensitivity", "NA"), "union_executed_origins": st.get("physical_origins", "NA"),
            "union_served_origins": st.get("reused_physical_origins", "NA"),
            "union_origin_cap_margin": str(cap - used) if cap is not None and used is not None else "NA",
            "refinement_executed_origins": st.get("repair_cumulative_physical_origins", "NA"),
            "refinement_origin_margin": st.get("repair_margin_physical_origins", "NA"),
            "admission_singletons_not_retained": st.get("repair_admission_singleton_misses_not_retained", "NA"),
            "admission_retained": st.get("repair_admission_retained_while_restricted", "NA"),
            "model_fingerprint": st.get("model_fingerprint", "NA"),
            "inputs_sha256_digest": hashlib.sha256(inputs.encode()).hexdigest()[:16] if inputs else "NA",
            "settings_values_sha256": hashlib.sha256(values.encode()).hexdigest()[:16] if values else "NA",
            "prepared_locus": st.get("prepared_locus", "none"),
        }
        rows.append(row)
    with open(out, "w") as o:
        o.write("\t".join(COLUMNS) + "\n")
        for r in rows:
            o.write("\t".join(str(r[c]).replace("\t", " ") for c in COLUMNS) + "\n")
    print(open(out).read(), end="")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1], sys.argv[2]))
