#!/usr/bin/env python3
"""I2-H: is the fixed-geometry emitting recursion scientifically invisible on one cell?

i2f_compare.compare (perf_compare's byte-identical scientific files -- dosage outputs, candidate tables and scores included --
status key by key, settings apart from budget keys, the complete adaptive trajectory, with the I2-A..I2-F provenance/work
keys allowed) plus exactly two allowed keys: gap_recursion_backend and gap_recursion_backend_scope (provenance; the first
also as its settings row). No status key began with gap_recursion_ before this batch.
Usage: i2h_compare.py <run dir A> <run dir B> [prefix]
"""
import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "i2f"))
import i2f_compare  # noqa: E402

perf_compare = i2f_compare.perf_compare
perf_compare.DIAGNOSTIC = re.compile(perf_compare.DIAGNOSTIC.pattern + r"|^gap_recursion_backend$|^gap_recursion_backend_scope$")
compare = i2f_compare.compare

if __name__ == "__main__":
    ok, problems, diagnostics, compared = compare(sys.argv[1], sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else "call")
    print(("IDENTICAL" if ok else "DIFFERENT") + f": {compared} files compared; allowed keys that changed: "
          f"{len(diagnostics)}" + ("" if ok else "; " + " | ".join(problems[:12])))
    sys.exit(0 if ok else 1)
