#!/usr/bin/env python3
"""I3-A: is analysing a sample against a prepared locus scientifically invisible?

i2k_compare.compare (perf_compare's byte-identical scientific files -- dosage and census tables, candidate tables and
scores -- status key by key, settings apart from budget keys, the complete adaptive trajectory's scientific columns,
with the I2-A..I2-K provenance/work/diagnostic additions allowed) plus exactly the prepared-locus provenance keys:
prepared_locus, prepared_locus_id, prepared_locus_key_sha256 (status). Timings (seconds_detail_model_*) are already
timing keys.
Usage: i3a_compare.py <run dir A> <run dir B> [prefix]
"""
import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "i2k"))
import i2k_compare  # noqa: E402

perf_compare = i2k_compare.perf_compare
perf_compare.DIAGNOSTIC = re.compile(perf_compare.DIAGNOSTIC.pattern + r"|^prepared_locus$|^prepared_locus_id$|^prepared_locus_key_sha256$")
compare = i2k_compare.compare

if __name__ == "__main__":
    ok, problems, diagnostics, compared = compare(sys.argv[1], sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else "call")
    print(("IDENTICAL" if ok else "DIFFERENT") + f": {compared} files compared; allowed keys that changed: "
          f"{len(diagnostics)}" + ("" if ok else "; " + " | ".join(problems[:12])))
    sys.exit(0 if ok else 1)
