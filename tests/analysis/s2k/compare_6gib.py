#!/usr/bin/env python3
"""S2-K: the 6 GiB-envelope ACOT 30x confirmation vs S2-J's completed run (caller.acot30.reuse), same code policy.
Every scientific output byte-identical after path normalisation (genotypes, block genotypes, candidates, scores,
equivalence set, locus genotype, sequence report files, dosage files, library sensitivity), the rounds table's scientific
columns identical, and every status key identical EXCEPT, explicitly:
  * the declared envelope: footprint_limit_bytes (status and settings: value and source);
  * S2-K's checkpoint rows: footprint_checkpoints and footprint_checkpoint<N>_* (S2-K adds post-refinement stages and the
    allocator field, so indices shift);
  * profile identity (profile_sha256, settings sources naming the arm's profile file);
  * timing, process resources and scheduling-variable counters exactly as comparator v2 declares them.
  compare_6gib.py <s2j run dir> <6 GiB run dir> <out.tsv>
"""
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s2j_compare_v2 as v2  # noqa: E402

CHECKPOINT_RE = re.compile(r"^footprint_checkpoint([0-9]+_.*|s)$")
ENVELOPE = {"footprint_limit_bytes", "profile_sha256"}


def allowed(k):
    return (k in ENVELOPE or CHECKPOINT_RE.match(k) or k in v2.RESOURCES or k in v2.SCHEDULING or k in v2.V2_ALLOWED or
            v2.TIMING_RE.search(k) or v2.RESOURCE_RE.match(k))


def main(a, b, out):
    skip = re.compile(r"^(profile\..*|run\.log|call\.refinement_.*|call\.repair_memory\.tsv)$")
    fa = {f for f in os.listdir(a) if not skip.match(f)}
    fb = {f for f in os.listdir(b) if not skip.match(f)}
    rows, bad = [], []
    if fa != fb:
        bad.append(f"file sets differ: {sorted(fa ^ fb)}")
    for f in sorted(fa & fb):
        pa, pb = os.path.join(a, f), os.path.join(b, f)
        if f == "call.mosaic_status.tsv":
            x, y = v2.kv(pa, a), v2.kv(pb, b)
            unexpected = [k for k in sorted(set(x) | set(y)) if x.get(k) != y.get(k) and not allowed(k)]
            rows.append(f"{f}\tunexpected_differing_keys\t{','.join(unexpected) or 'none'}")
            if unexpected:
                bad.append("status keys outside the allowances: " + ",".join(unexpected))
            continue
        if f == "call.settings.tsv":
            x = {k: v.split("\t")[0] for k, v in v2.kv(pa, a).items()}
            y = {k: v.split("\t")[0] for k, v in v2.kv(pb, b).items()}
            x.pop("footprint_limit_bytes", None)
            y.pop("footprint_limit_bytes", None)
            if x != y:
                bad.append("settings values differ beyond the declared envelope")
            continue
        if f == "call.repair_rounds.tsv":
            cut = lambda p: [line.split("\t")[:28] for line in open(p).read().splitlines()]  # noqa: E731
            if cut(pa) != cut(pb):
                bad.append("trajectory differs")
            continue
        if f == "call.dosage_status.tsv":
            x = {k: v for k, v in v2.kv(pa, a).items() if not k.startswith("seconds_")}
            y = {k: v for k, v in v2.kv(pb, b).items() if not k.startswith("seconds_")}
            if x != y:
                bad.append("dosage status differs")
            continue
        if open(pa, errors="replace").read().replace(a, "<D>") != open(pb, errors="replace").read().replace(b, "<D>"):
            bad.append(f"{f} differs")
    with open(out, "w") as o:
        o.write("key\tvalue\ncomparator\tcompare_6gib\n")
        o.write(f"science_identical\t{'yes' if not bad else 'NO: ' + ' | '.join(bad)}\n")
        for r in rows:
            o.write(r + "\n")
    print(open(out).read(), end="")
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main(os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2]), sys.argv[3]))
