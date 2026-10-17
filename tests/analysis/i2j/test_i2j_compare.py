#!/usr/bin/env python3
"""Negative controls for i2j_compare's census-backend additions (and the inherited controls)."""
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
for d in (HERE, os.path.join(HERE, "..", "i2i"), os.path.join(HERE, "..", "i2h"), os.path.join(HERE, "..", "i2f"),
          os.path.join(HERE, "..", "i2e"), os.path.join(HERE, "..", "i2d"), os.path.join(HERE, "..", "i2b"),
          os.path.join(HERE, "..", "i2a")):
    sys.path.insert(0, d)
import i2j_compare  # noqa: E402
import test_i2a_compare as base  # noqa: E402
import test_i2i_compare as previous  # noqa: E402


def main():
    root = tempfile.mkdtemp()
    try:
        ref = base.run(root, "ref")
        cases = [
            ("backend provenance and work", True, dict(status=base.STATUS + "dosage_census_backend\tinterned_experimental\n"
                                                                "dosage_census_work_pairs_compared\t9\nseconds_dosage_census_pairs\t1\n")),
            ("backend settings row", True, dict(settings=base.SETTINGS + "dosage_census_backend\tinterned_experimental\tprofile x\n")),
            ("a near-miss key", False, dict(status=base.STATUS + "dosage_census_backends\t1\n")),
            ("a census result key", False, dict(status=base.STATUS + "dosage_census_sensitive_blocks\t3\n")),
            ("a scientific key", False, dict(status=base.STATUS.replace("COMPLETE", "TIED"))),
        ]
        for i, (label, expect, kw) in enumerate(cases):
            ok = i2j_compare.compare(ref, base.run(root, f"c{i}", **kw))[0]
            assert ok == expect, f"{label}: compare returned {ok}, expected {expect}"
        print(f"i2j_compare: {len(cases)} controls as expected")
    finally:
        shutil.rmtree(root)
    previous.main()


if __name__ == "__main__":
    main()
