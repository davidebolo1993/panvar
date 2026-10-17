#!/usr/bin/env python3
"""I2-K: is the union-to-adaptive emission handoff scientifically invisible on one cell?

i2j_compare.compare (perf_compare's byte-identical scientific files -- dosage and census tables, candidate tables and
scores -- status key by key, settings apart from budget keys, the complete adaptive trajectory's scientific columns,
with the I2-A..I2-J provenance/work/diagnostic additions allowed) plus exactly the handoff's provenance and work block:
adaptive_emission_handoff and adaptive_emission_handoff_* (status and settings row). The step-1 storage keys
(indexed_union_emission_cache_*, union_end_phys_footprint_bytes) are already diagnostic.
Usage: i2k_compare.py <run dir A> <run dir B> [prefix]
"""
import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "i2j"))
import i2j_compare  # noqa: E402

perf_compare = i2j_compare.perf_compare
perf_compare.DIAGNOSTIC = re.compile(perf_compare.DIAGNOSTIC.pattern + r"|^adaptive_emission_handoff$|^adaptive_emission_handoff_")
compare = i2j_compare.compare

if __name__ == "__main__":
    ok, problems, diagnostics, compared = compare(sys.argv[1], sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else "call")
    print(("IDENTICAL" if ok else "DIFFERENT") + f": {compared} files compared; allowed keys that changed: "
          f"{len(diagnostics)}" + ("" if ok else "; " + " | ".join(problems[:12])))
    sys.exit(0 if ok else 1)
