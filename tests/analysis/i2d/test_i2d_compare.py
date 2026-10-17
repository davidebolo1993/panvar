#!/usr/bin/env python3
"""Negative controls for i2d_compare's two allowed keys (and the inherited I2-B/I2-A controls)."""
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
for d in (HERE, os.path.join(HERE, "..", "i2b"), os.path.join(HERE, "..", "i2a")):
    sys.path.insert(0, d)
import i2d_compare  # noqa: E402
import test_i2a_compare as base  # noqa: E402
import test_i2b_compare as previous  # noqa: E402


def main():
    root = tempfile.mkdtemp()
    try:
        ref = base.run(root, "ref")
        cases = [
            ("fused provenance and placement counter", True, dict(status=base.STATUS + "adaptive_fused_seed\t1\nadaptive_full_placements\t2\n")),
            ("fused settings row", True, dict(settings=base.SETTINGS + "adaptive_fused_seed\t1\tprofile x\n")),
            ("a near-miss key", False, dict(status=base.STATUS + "adaptive_full_placements_extra\tx\n")),
            ("a scientific key", False, dict(status=base.STATUS.replace("COMPLETE", "REPAIR_INCOMPLETE"))),
        ]
        for i, (label, expect, kw) in enumerate(cases):
            ok = i2d_compare.compare(ref, base.run(root, f"c{i}", **kw))[0]
            assert ok == expect, f"{label}: compare returned {ok}, expected {expect}"
        print(f"i2d_compare: {len(cases)} controls as expected")
    finally:
        shutil.rmtree(root)
    previous.main()


if __name__ == "__main__":
    main()
