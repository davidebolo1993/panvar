#!/usr/bin/env python3
"""CYP2D6/D7 pilot, AMENDMENT_4: the state-class capacity sentinel gate (call outputs only; no truth is read).

Usage: sentinel_gate_a4.py <pilot package> <state-class cap> <time ceiling seconds>
Cells: the 30x sentinels lzo.HG02769.30, mosaic.mosaic_one_1.30, loo.NA19983.30 and their 5x cells, all in
data/runs_a4/; the 5x cells are compared with the amendment-3 run (data/runs_a3/, cap 200 M).

30x, per cell (reported in sentinel_gate_a4_counters.tsv):
  S1 completed (exit 0 or 1) and no declared resource cap fired: the status note names no cap and every union-scoring
     counter is within its cap (statistical_state_classes, physical_origins, piece_hits) and every repair margin is
     non-negative;
  S2 statistical_state_classes < the state-class cap (consumed and margin reported);
  S3 peak RSS < 4 GiB and wall < the time ceiling;
  LZO the ordinary LZO 30x sentinel emits an authoritative call: exit 0, call_class strict_call or
      genotype_equivalent_call.
5x, per cell:
  S4 every scientific output is byte-identical to the amendment-3 run. Permitted differences, and only these: fields
     recording the declared cap or the profile that declares it (max_state_classes, time_limit_seconds, profile,
     profile_name, profile_sha256, their settings.tsv rows), time and memory measurements, and the run directory name.
Writes sentinel_gate_a4.tsv and sentinel_gate_a4_counters.tsv into the package; exit 0 only if every check passes.
"""
import os
import re
import sys

TIMING = re.compile(r"second|wall|rss|footprint|phys|_gib|timing|_bytes_peak|peak_bytes|elapsed|throughput|per_second")
DECLARED = re.compile(r"^(max_state_classes|time_limit_seconds|profile|profile_name|profile_sha256)$")
GIB4 = 4 * 2 ** 30
J30 = ["lzo.HG02769.30", "mosaic.mosaic_one_1.30", "loo.NA19983.30"]
J5 = ["lzo.HG02769.5", "mosaic.mosaic_one_1.5", "loo.NA19983.5"]
COUNTERS = ["status", "call_class", "no_call_reason", "note", "statistical_state_classes", "merged_state_classes",
            "physical_origins", "piece_hits", "reused_latent_alignment_pairs", "max_state_classes",
            "max_physical_origins", "max_piece_hits", "max_mass_cells", "scorer_mass_matrix_bytes",
            "scorer_placement_rows_peak", "scorer_proposed_starts_peak", "repair_rounds_run",
            "repair_cumulative_physical_origins", "repair_margin_physical_origins", "repair_cumulative_state_classes",
            "repair_margin_state_classes", "repair_cumulative_piece_hits", "repair_margin_piece_hits",
            "repair_cumulative_latent_pairs", "repair_margin_latent_pairs", "tie_branches_run",
            "seconds_total", "resources_rss_peak_bytes", "resources_phys_footprint_peak_bytes"]


def kv(path):
    out = {}
    for l in open(path).read().splitlines():
        p = l.split("\t")
        if len(p) >= 2:
            out[p[0]] = p[1]
    return out


def run_rec(path):
    lines = open(path).read().splitlines()
    return dict(zip(lines[0].split("\t"), lines[1].split("\t")))


def normalised(path, a, b):
    lines = open(path, errors="replace").read().replace(b, a).splitlines()
    if not lines or "\t" not in lines[0]:
        return "\n".join(lines)
    head = lines[0].split("\t")
    if head[:2] == ["key", "value"] or os.path.basename(path).endswith(".settings.tsv"):
        return "\n".join(l for l in lines
                         if not (TIMING.search(l.split("\t")[0]) or DECLARED.match(l.split("\t")[0])))
    keep = [i for i, h in enumerate(head) if not (TIMING.search(h) or DECLARED.match(h))]
    return "\n".join("\t".join(p[i] for i in keep if i < len(p)) for p in (l.split("\t") for l in lines))


def main(pkg, cap, ceiling):
    rows, counters = [], []

    def check(job, gate, ok, detail):
        rows.append((job, gate, "PASS" if ok else "FAIL", detail))

    for job in J30:
        prefix = os.path.join(pkg, "data", "runs_a4", job)
        rec, st = run_rec(prefix + ".run.tsv"), kv(prefix + ".mosaic_status.tsv")
        counters.append({"job": job, "exit": rec["exit"], "wall_seconds": rec["wall_seconds"],
                         "peak_rss_bytes": rec["peak_rss_bytes"], **{k: st.get(k, "NA") for k in COUNTERS},
                         "state_class_margin": cap - int(st["statistical_state_classes"])
                         if st.get("statistical_state_classes", "NA").isdigit() else "NA"})
        note = st.get("note", "")
        within = []
        for used, limit in (("statistical_state_classes", "max_state_classes"), ("physical_origins", "max_physical_origins"),
                            ("piece_hits", "max_piece_hits")):
            u, l_ = st.get(used, "NA"), st.get(limit, "NA")
            within.append(u.isdigit() and l_.isdigit() and int(u) <= int(l_))
        margins = [int(v) >= 0 for k, v in st.items() if k.startswith("repair_margin_") and v.lstrip("-").isdigit()]
        fired = (rec["exit"] not in ("0", "1") or re.search(r"\bcap\b|limit", note) is not None
                 or not all(within) or not all(margins))
        check(job, "S1", not fired, f"exit {rec['exit']}; status {st.get('status')}; call_class {st.get('call_class')}; "
                                    f"no_call_reason {st.get('no_call_reason')}; note: {note[:120]}")
        sc = st.get("statistical_state_classes", "NA")
        check(job, "S2", sc.isdigit() and int(sc) < cap, f"state classes {sc} of {cap}")
        rss, wall = int(rec["peak_rss_bytes"]), float(rec["wall_seconds"])
        check(job, "S3", 0 < rss < GIB4 and wall < ceiling, f"peak RSS {rss / 2 ** 30:.2f} GiB, wall {wall:.0f} s")
        if job.startswith("lzo."):
            check(job, "LZO", rec["exit"] == "0" and st.get("call_class") in ("strict_call", "genotype_equivalent_call"),
                  f"exit {rec['exit']}, call_class {st.get('call_class')}")
    for job in J5:
        d3, d4 = os.path.join(pkg, "data", "runs_a3"), os.path.join(pkg, "data", "runs_a4")
        skip = (job + ".call.log", job + ".run.tsv")
        f3 = sorted(f for f in os.listdir(d3) if f.startswith(job + ".") and f not in skip)
        f4 = sorted(f for f in os.listdir(d4) if f.startswith(job + ".") and f not in skip)
        identical = [f for f in f4 if f in f3 and open(os.path.join(d3, f), "rb").read() == open(os.path.join(d4, f), "rb").read()]
        differ = [f for f in f4 if f in f3 and f not in identical and
                  normalised(os.path.join(d3, f), "runs_a3", "runs_a4") != normalised(os.path.join(d4, f), "runs_a3", "runs_a4")]
        check(job, "S4", f3 == f4 and not differ,
              f"{len(f4)} files; {len(identical)} byte-identical; the rest differ only in declared-cap, profile, "
              f"time/memory or path fields" if f3 == f4 and not differ else
              f"file sets equal {f3 == f4}; scientific differences in {differ}")
    with open(os.path.join(pkg, "sentinel_gate_a4.tsv"), "w") as f:
        f.write("job\tgate\tstatus\tdetail\n" + "".join("\t".join(r) + "\n" for r in rows))
    with open(os.path.join(pkg, "sentinel_gate_a4_counters.tsv"), "w") as f:
        keys = list(counters[0])
        f.write("\t".join(keys) + "\n" + "".join("\t".join(str(c[k]) for k in keys) + "\n" for c in counters))
    for r in rows:
        print("\t".join(r))
    failed = [r for r in rows if r[2] != "PASS"]
    print(f"{len(rows)} checks, {len(failed)} failed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(os.path.abspath(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3])))
