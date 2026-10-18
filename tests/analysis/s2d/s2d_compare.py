#!/usr/bin/env python3
"""S2-D: complete-run identity between the control (early_incumbent, bulk indexed union) and the progressive caller run.
Every output byte-identical (paths normalised) except: the status and settings (declared rows, timing/resource keys
listed), repair_memory (timing) and dosage_status timing keys. Prints the differing status keys.
  s2d_compare.py <control dir> <progressive dir> <out.tsv>
"""
import os
import sys

DECLARED = ("union_progressive", "checkpoint_")   # progressive-only rows and the extra prefix checkpoint commit
VOLATILE = ("seconds", "footprint", "rss", "peak", "_bytes", "worker", "threads_", "profile", "timing", "elapsed")


def kv(p):
    out = {}
    for line in open(p, errors="replace").read().splitlines():
        a = line.split("\t", 1)
        if len(a) == 2:
            out[a[0]] = a[1]
    return out


def main(a, b, out):
    fa = {f for f in os.listdir(a) if f.startswith("call.")}
    fb = {f for f in os.listdir(b) if f.startswith("call.")}
    rows, bad = [], []
    if fa != fb:
        bad.append(f"file sets differ: {sorted(fa ^ fb)}")
    for f in sorted(fa & fb):
        if f in ("call.repair_memory.tsv",):
            continue
        ta = open(os.path.join(a, f), errors="replace").read().replace(a, "<D>")
        tb = open(os.path.join(b, f), errors="replace").read().replace(b, "<D>")
        if f in ("call.mosaic_status.tsv", "call.dosage_status.tsv", "call.settings.tsv"):
            ka, kb = kv(os.path.join(a, f)), kv(os.path.join(b, f))
            keys = sorted(k for k in set(ka) | set(kb)
                          if not any(x in k for x in VOLATILE) and not k.startswith(DECLARED) and
                          ka.get(k, "").replace(a, "<D>") != kb.get(k, "").replace(b, "<D>"))
            if f == "call.settings.tsv":
                keys = [k for k in keys if k not in ("union_progressive_schedule",)]
            rows.append(f"{f}\tdiffering_keys\t{','.join(keys) or 'none'}")
        elif f == "call.repair_rounds.tsv":
            cut = lambda t: [line.split("\t")[:28] for line in t.splitlines()]  # noqa: E731
            if cut(ta) != cut(tb):
                bad.append("trajectory differs")
        elif ta != tb:
            bad.append(f"{f} differs")
    with open(out, "w") as o:
        o.write("key\tvalue\n")
        o.write(f"outputs_identical\t{'yes' if not bad else 'NO: ' + ' | '.join(bad)}\n")
        for r in rows:
            o.write(r + "\n")
    print(open(out).read(), end="")
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main(os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2]), sys.argv[3]))
