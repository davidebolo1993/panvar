#!/usr/bin/env python3
"""I2-I: does the verification census leave every scientific output unchanged?

i2h_compare.compare (perf_compare's byte-identical scientific files -- dosage outputs, candidate tables and scores included --
status key by key, settings apart from budget keys, the complete adaptive trajectory, with the I2-A..I2-H provenance/work
keys allowed) plus exactly two diagnostic additions: the census's own output file (<prefix>.verification_census.tsv, which
only a census run writes) and its settings row (verification_census_period). Nothing else may differ.
Usage: i2i_compare.py <run dir A> <run dir B> [prefix]
"""
import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "i2h"))
import i2h_compare  # noqa: E402

perf_compare = i2h_compare.perf_compare
perf_compare.DIAGNOSTIC = re.compile(perf_compare.DIAGNOSTIC.pattern + r"|^verification_census_period$")
CENSUS_SUFFIX = ".verification_census.tsv"
_listdir = os.listdir


def compare(dir_a, dir_b, prefix="call"):
    census = prefix + CENSUS_SUFFIX
    os.listdir = lambda d: [f for f in _listdir(d) if f != census]   # the census file is the census arm's own output
    try:
        return i2h_compare.compare(dir_a, dir_b, prefix)
    finally:
        os.listdir = _listdir


if __name__ == "__main__":
    ok, problems, diagnostics, compared = compare(sys.argv[1], sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else "call")
    print(("IDENTICAL" if ok else "DIFFERENT") + f": {compared} files compared; allowed keys that changed: "
          f"{len(diagnostics)}" + ("" if ok else "; " + " | ".join(problems[:12])))
    sys.exit(0 if ok else 1)
