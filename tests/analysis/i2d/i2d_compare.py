#!/usr/bin/env python3
"""I2-D: is fused anchor seeding scientifically invisible on one cell?

i2b_compare.compare (perf_compare's byte-identical scientific files, key-by-key status, settings apart from budget keys,
and the complete adaptive trajectory's scientific columns, with the I2-A/I2-B provenance blocks allowed), with exactly two
further allowed keys: adaptive_fused_seed (provenance; also its settings row) and adaptive_full_placements (the work
counter the change reduces). Nothing else may differ. Usage: i2d_compare.py <run dir A> <run dir B> [prefix]
"""
import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "i2b"))
import i2b_compare  # noqa: E402

perf_compare = i2b_compare.perf_compare
perf_compare.DIAGNOSTIC = re.compile(perf_compare.DIAGNOSTIC.pattern + r"|^adaptive_fused_seed$|^adaptive_full_placements$")
compare = i2b_compare.compare

if __name__ == "__main__":
    ok, problems, diagnostics, compared = compare(sys.argv[1], sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else "call")
    print(("IDENTICAL" if ok else "DIFFERENT") + f": {compared} files compared; allowed keys that changed: "
          f"{len(diagnostics)}" + ("" if ok else "; " + " | ".join(problems[:12])))
    sys.exit(0 if ok else 1)
