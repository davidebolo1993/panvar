#!/usr/bin/env python3
"""CYP2D6/D7 pilot, AMENDMENT_3: the amended sentinel gate (call outputs only; no truth is read).

Usage: sentinel_gate_a3.py <pilot package> <tag> <repeat tag> <declared cap> <job>...
For each job, from data/runs_<tag>/<job>.* and its repeat in data/runs_<repeat tag>/:
  G1 the cell completed (exit 0 call or 1 no-call); a no-call's reason is recorded, and it is never neighbourhood_bound;
  G2 the run's declared repair_max_neighbours equals the declared cap and is >= its repair_neighbour_bound;
  G3 every round of every repair-round table (.repair_rounds.tsv, and its .original copy) has
     neighbours_generated <= the declared cap, and generated = scored + duplicate + suppressed_visited + self_skipped
     + rejected_infinite_prior (no neighbour is dropped for any other reason, so none is truncated by the cap);
  G4 peak RSS < 4 GiB and wall < 1200 s;
  G5 deterministic: every output file of the run exists in the repeat and is byte-identical, except fields that measure
     time or memory (names matching TIMING) and the run's own output directory, which are compared after removal.
Also checks that the ordinary LZO sentinel (lzo.*) produced a call (exit 0).
Writes <pilot package>/sentinel_gate_a3.tsv; exit 0 only when every check passes.
"""
import os
import re
import sys

TIMING = re.compile(r"second|wall|rss|footprint|phys|_gib|timing|_bytes_peak|peak_bytes|elapsed|throughput|per_second")
GIB4 = 4 * 2 ** 30


def table(path):
    lines = [l for l in open(path).read().splitlines() if l and not l.startswith("#")]
    head = lines[0].split("\t")
    return head, [dict(zip(head, l.split("\t"))) for l in lines[1:]]


def kv(path):
    out = {}
    for l in open(path).read().splitlines():
        p = l.split("\t")
        if len(p) >= 2:
            out[p[0]] = p[1]
    return out


def normalised(path, run_dir, repeat_dir):
    """File content with the run directory name unified and time/memory fields removed (key/value or column tables)."""
    text = open(path, errors="replace").read().replace(repeat_dir, run_dir)
    lines = text.splitlines()
    if not lines or "\t" not in lines[0]:
        return text
    head = lines[0].split("\t")
    if head[:2] == ["key", "value"]:
        return "\n".join(l for l in lines if not TIMING.search(l.split("\t")[0]))
    keep = [i for i, h in enumerate(head) if not TIMING.search(h)]
    return "\n".join("\t".join(p[i] for i in keep if i < len(p)) for p in (l.split("\t") for l in lines))


def main(pkg, tag, rtag, cap, jobs):
    rows = []

    def check(job, gate, ok, detail):
        rows.append((job, gate, "PASS" if ok else "FAIL", detail))

    run_dir, repeat_dir = f"runs_{tag}", f"runs_{rtag}"
    for job in jobs:
        prefix = os.path.join(pkg, "data", run_dir, job)
        rprefix = os.path.join(pkg, "data", repeat_dir, job)
        rec = table(prefix + ".run.tsv")[1][0]
        st = kv(prefix + ".mosaic_status.tsv")
        reason = st.get("no_call_reason", "NA")
        check(job, "G1", rec["exit"] in ("0", "1") and reason != "neighbourhood_bound",
              f"exit {rec['exit']}; status {st.get('status')}; no_call_reason {reason}")
        declared, bound = int(st.get("repair_max_neighbours", -1)), int(st.get("repair_neighbour_bound", 10 ** 12))
        check(job, "G2", declared == cap and declared >= bound, f"declared {declared}, model bound {bound}")
        tables = [p for p in (prefix + ".repair_rounds.tsv", prefix + ".original.repair_rounds.tsv") if os.path.exists(p)]
        rounds, problems = 0, []
        for t in tables:
            for r in table(t)[1]:
                rounds += 1
                g = int(r["neighbours_generated"])
                parts = sum(int(r[k]) for k in ("neighbours_scored", "duplicate", "suppressed_visited", "self_skipped",
                                                  "rejected_infinite_prior"))
                if g > cap or g != parts:
                    problems.append(f"{os.path.basename(t)} round {r['round']}: generated {g}, accounted {parts}")
        check(job, "G3", not problems and (rounds > 0 or rec["exit"] == "1"),
              f"{rounds} rounds in {len(tables)} table(s); max generated "
              f"{max([int(r['neighbours_generated']) for t in tables for r in table(t)[1]] or [0])}; "
              + ("; ".join(problems) or "every neighbour accounted"))
        rss, wall = int(rec["peak_rss_bytes"]), float(rec["wall_seconds"])
        check(job, "G4", 0 < rss < GIB4 and wall < 1200, f"peak RSS {rss / 2 ** 30:.2f} GiB, wall {wall:.0f} s")
        d, rd = os.path.dirname(prefix), os.path.dirname(rprefix)
        mine = sorted(f for f in os.listdir(d) if f.startswith(job + ".") and f not in (job + ".call.log", job + ".run.tsv"))
        theirs = sorted(f for f in os.listdir(rd) if f.startswith(job + ".") and f not in (job + ".call.log", job + ".run.tsv"))
        differ = [f for f in mine if f in theirs and
                  normalised(os.path.join(d, f), run_dir, repeat_dir) != normalised(os.path.join(rd, f), run_dir, repeat_dir)]
        identical = [f for f in mine if f in theirs and open(os.path.join(d, f), "rb").read() == open(os.path.join(rd, f), "rb").read()]
        check(job, "G5", mine == theirs and not differ,
              f"{len(mine)} files; {len(identical)} byte-identical; differing beyond time/memory/path: {differ or 'none'}"
              + ("" if mine == theirs else f"; file sets differ: {sorted(set(mine) ^ set(theirs))}"))
        if job.startswith("lzo."):
            check(job, "LZO", rec["exit"] == "0", f"ordinary LZO sentinel exit {rec['exit']}")
    with open(os.path.join(pkg, "sentinel_gate_a3.tsv"), "w") as f:
        f.write("job\tgate\tstatus\tdetail\n" + "".join("\t".join(r) + "\n" for r in rows))
    for r in rows:
        print("\t".join(r))
    failed = [r for r in rows if r[2] != "PASS"]
    print(f"{len(rows)} checks, {len(failed)} failed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(os.path.abspath(sys.argv[1]), sys.argv[2], sys.argv[3], int(sys.argv[4]), sys.argv[5:]))
