#!/usr/bin/env python3
"""Negative controls for i2i_compare's two diagnostic additions (and the inherited controls)."""
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
for d in (HERE, os.path.join(HERE, "..", "i2h"), os.path.join(HERE, "..", "i2f"), os.path.join(HERE, "..", "i2e"),
          os.path.join(HERE, "..", "i2d"), os.path.join(HERE, "..", "i2b"), os.path.join(HERE, "..", "i2a")):
    sys.path.insert(0, d)
import i2i_compare  # noqa: E402
import test_i2a_compare as base  # noqa: E402
import test_i2h_compare as previous  # noqa: E402


def main():
    root = tempfile.mkdtemp()
    try:
        ref = base.run(root, "ref")
        cases = [
            ("the census file on one side", True, dict(extra="call.verification_census.tsv")),
            ("the census settings row", True, dict(settings=base.SETTINGS + "verification_census_period\t128\tprofile x\n")),
            ("another output on one side", False, dict(extra="call.verification_census_rows.tsv")),
            ("a census status key", False, dict(status=base.STATUS + "verification_census_requests\t9\n")),
            ("a scientific key", False, dict(status=base.STATUS.replace("COMPLETE", "TIED"))),
        ]
        for i, (label, expect, kw) in enumerate(cases):
            ok = i2i_compare.compare(ref, base.run(root, f"c{i}", **kw))[0]
            assert ok == expect, f"{label}: compare returned {ok}, expected {expect}"
        print(f"i2i_compare: {len(cases)} controls as expected")
    finally:
        shutil.rmtree(root)
    previous.main()


if __name__ == "__main__":
    main()
