#!/usr/bin/env python3
"""I3-C: is analysing a sample against a LOADED index scientifically invisible?

i3a_compare.compare (byte-identical scientific files, status key by key, settings apart from budget keys, the complete
adaptive trajectory's scientific columns, the I2-A..I3-A diagnostic additions) plus exactly the I3-C runtime keys:
prepared_locus_origin and prepared_locus_file_sha256 (status).
Usage: i3c_compare.py <run dir A> <run dir B> [prefix]
"""
import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "i3a"))
import i3a_compare  # noqa: E402

perf_compare = i3a_compare.perf_compare
perf_compare.DIAGNOSTIC = re.compile(perf_compare.DIAGNOSTIC.pattern + r"|^prepared_locus_origin$|^prepared_locus_file_sha256$")
compare = i3a_compare.compare

if __name__ == "__main__":
    ok, problems, diagnostics, compared = compare(sys.argv[1], sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else "call")
    print(("IDENTICAL" if ok else "DIFFERENT") + f": {compared} files compared; allowed keys that changed: "
          f"{len(diagnostics)}" + ("" if ok else "; " + " | ".join(problems[:12])))
    sys.exit(0 if ok else 1)
