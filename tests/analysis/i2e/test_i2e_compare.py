#!/usr/bin/env python3
"""Negative controls for i2e_compare's one allowed family (and the inherited controls)."""
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
for d in (HERE, os.path.join(HERE, "..", "i2d"), os.path.join(HERE, "..", "i2b"), os.path.join(HERE, "..", "i2a")):
    sys.path.insert(0, d)
import i2e_compare  # noqa: E402
import test_i2a_compare as base  # noqa: E402
import test_i2d_compare as previous  # noqa: E402


def main():
    root = tempfile.mkdtemp()
    try:
        ref = base.run(root, "ref")
        cases = [
            ("reuse block", True, dict(status=base.STATUS + "marker_proposal_reuse\ton\nmarker_proposal_contributions_reused\t9\n")),
            ("reuse settings row", True, dict(settings=base.SETTINGS + "marker_proposal_reuse_table_bytes\t4\tprofile x\n")),
            ("a candidate table", False, dict(genotypes="g\t9\n")),
            ("a scientific key", False, dict(status=base.STATUS.replace("strict_call", "no_call"))),
        ]
        for i, (label, expect, kw) in enumerate(cases):
            ok = i2e_compare.compare(ref, base.run(root, f"c{i}", **kw))[0]
            assert ok == expect, f"{label}: compare returned {ok}, expected {expect}"
        print(f"i2e_compare: {len(cases)} controls as expected")
    finally:
        shutil.rmtree(root)
    previous.main()


if __name__ == "__main__":
    main()
