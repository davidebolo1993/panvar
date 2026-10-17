#!/usr/bin/env python3
"""I2-B: are the control (indexed union, reference adaptive search) and the fine arm (indexed union, fine adaptive window
discovery) scientifically identical on one cell?

i2a_compare.compare -- perf_compare's byte-identical scientific files, key-by-key status, settings apart from budget
keys, and the complete adaptive trajectory's scientific columns (.repair_rounds.tsv) -- with exactly one further
extension of the allowed class: the adaptive window block adaptive_window_*, adaptive_fine_* and
worker_seconds_adaptive_* (status), and the three experimental profile rows adaptive_window_discovery /
adaptive_fine_store_bytes / adaptive_fine_max_served_placements (settings). Nothing else may differ.
Usage: i2b_compare.py <run dir A> <run dir B> [prefix, default call]
"""
import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "i2a"))
import i2a_compare  # noqa: E402

perf_compare = i2a_compare.perf_compare
ADAPTIVE = r"|^adaptive_window_|^adaptive_fine_|^worker_seconds_adaptive_"
perf_compare.DIAGNOSTIC = re.compile(perf_compare.DIAGNOSTIC.pattern + ADAPTIVE)
compare = i2a_compare.compare

if __name__ == "__main__":
    ok, problems, diagnostics, compared = compare(sys.argv[1], sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else "call")
    print(("IDENTICAL" if ok else "DIFFERENT") + f": {compared} files compared; allowed keys that changed: "
          f"{len(diagnostics)}" + ("" if ok else "; " + " | ".join(problems[:12])))
    sys.exit(0 if ok else 1)
