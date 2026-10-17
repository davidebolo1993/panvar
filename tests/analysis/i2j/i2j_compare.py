#!/usr/bin/env python3
"""I2-J: is the interned dosage census scientifically invisible on one cell?

i2i_compare.compare (perf_compare's byte-identical scientific files -- dosage tables including dosage_census.tsv,
candidate tables and scores -- status key by key, settings apart from budget keys, the complete adaptive trajectory,
with the I2-A..I2-I provenance/work/diagnostic additions allowed) plus exactly the census backend's provenance and work:
dosage_census_backend (dosage status and settings row) and dosage_census_work_* (work counters). Its coarse timers,
seconds_dosage_census_*, are already timing keys.
Usage: i2j_compare.py <run dir A> <run dir B> [prefix]
"""
import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "i2i"))
import i2i_compare  # noqa: E402

perf_compare = i2i_compare.perf_compare
perf_compare.DIAGNOSTIC = re.compile(perf_compare.DIAGNOSTIC.pattern + r"|^dosage_census_backend$|^dosage_census_work_")
compare = i2i_compare.compare

if __name__ == "__main__":
    ok, problems, diagnostics, compared = compare(sys.argv[1], sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else "call")
    print(("IDENTICAL" if ok else "DIFFERENT") + f": {compared} files compared; allowed keys that changed: "
          f"{len(diagnostics)}" + ("" if ok else "; " + " | ".join(problems[:12])))
    sys.exit(0 if ok else 1)
