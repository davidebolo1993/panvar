#!/usr/bin/env python3
"""I3-D: is analysing a sample against a panel-loaded (or v1-loaded) index scientifically invisible?

i3c_compare.compare plus exactly the I3-D runtime key prepared_locus_panel_source (status).
Usage: i3d_compare.py <run dir A> <run dir B> [prefix]
"""
import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "i3c"))
import i3c_compare  # noqa: E402

perf_compare = i3c_compare.perf_compare
perf_compare.DIAGNOSTIC = re.compile(perf_compare.DIAGNOSTIC.pattern + r"|^prepared_locus_panel_source$")
compare = i3c_compare.compare

if __name__ == "__main__":
    ok, problems, diagnostics, compared = compare(sys.argv[1], sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else "call")
    print(("IDENTICAL" if ok else "DIFFERENT") + f": {compared} files compared; allowed keys that changed: "
          f"{len(diagnostics)}" + ("" if ok else "; " + " | ".join(problems[:12])))
    sys.exit(0 if ok else 1)
