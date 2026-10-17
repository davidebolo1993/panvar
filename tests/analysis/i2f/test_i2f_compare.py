#!/usr/bin/env python3
"""Negative controls for i2f_compare's one allowed family (and the inherited controls)."""
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
for d in (HERE, os.path.join(HERE, "..", "i2e"), os.path.join(HERE, "..", "i2d"), os.path.join(HERE, "..", "i2b"),
          os.path.join(HERE, "..", "i2a")):
    sys.path.insert(0, d)
import i2f_compare  # noqa: E402
import test_i2a_compare as base  # noqa: E402
import test_i2e_compare as previous  # noqa: E402


def main():
    root = tempfile.mkdtemp()
    try:
        ref = base.run(root, "ref")
        cases = [
            ("combination block", True, dict(status=base.STATUS + "combination_prepared_library_terms\t1\n"
                                                      "combination_invariant_transcendentals_avoided\t9\n")),
            ("combination settings row", True, dict(settings=base.SETTINGS + "combination_prepared_library_terms\t1\tprofile x\n")),
            ("a score key outside the family", False, dict(status=base.STATUS + "prepared_combination_score\t1\n")),
            ("a scientific key", False, dict(status=base.STATUS.replace("COMPLETE", "TIED"))),
        ]
        for i, (label, expect, kw) in enumerate(cases):
            ok = i2f_compare.compare(ref, base.run(root, f"c{i}", **kw))[0]
            assert ok == expect, f"{label}: compare returned {ok}, expected {expect}"
        print(f"i2f_compare: {len(cases)} controls as expected")
    finally:
        shutil.rmtree(root)
    previous.main()


if __name__ == "__main__":
    main()
