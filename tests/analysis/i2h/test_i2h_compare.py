#!/usr/bin/env python3
"""Negative controls for i2h_compare's two allowed keys (and the inherited controls)."""
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
for d in (HERE, os.path.join(HERE, "..", "i2f"), os.path.join(HERE, "..", "i2e"), os.path.join(HERE, "..", "i2d"),
          os.path.join(HERE, "..", "i2b"), os.path.join(HERE, "..", "i2a")):
    sys.path.insert(0, d)
import i2h_compare  # noqa: E402
import test_i2a_compare as base  # noqa: E402
import test_i2f_compare as previous  # noqa: E402


def main():
    root = tempfile.mkdtemp()
    try:
        ref = base.run(root, "ref")
        cases = [
            ("backend provenance", True, dict(status=base.STATUS + "gap_recursion_backend\tfixed_geometry_experimental\n"
                                                       "gap_recursion_backend_scope\tunion_adaptive_dosage\n")),
            ("backend settings row", True, dict(settings=base.SETTINGS + "gap_recursion_backend\tfixed_geometry_experimental\tprofile x\n")),
            ("a near-miss key", False, dict(status=base.STATUS + "gap_recursion_backend_rows\t1\n")),
            ("a scientific key", False, dict(status=base.STATUS.replace("COMPLETE", "TIED"))),
        ]
        for i, (label, expect, kw) in enumerate(cases):
            ok = i2h_compare.compare(ref, base.run(root, f"c{i}", **kw))[0]
            assert ok == expect, f"{label}: compare returned {ok}, expected {expect}"
        print(f"i2h_compare: {len(cases)} controls as expected")
    finally:
        shutil.rmtree(root)
    previous.main()


if __name__ == "__main__":
    main()
