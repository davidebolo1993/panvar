#!/usr/bin/env python3
"""Negative controls for i3a_compare's prepared-locus provenance keys (and the inherited controls)."""
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
for d in (HERE, os.path.join(HERE, "..", "i2k"), os.path.join(HERE, "..", "i2j"), os.path.join(HERE, "..", "i2i"),
          os.path.join(HERE, "..", "i2h"), os.path.join(HERE, "..", "i2f"), os.path.join(HERE, "..", "i2e"),
          os.path.join(HERE, "..", "i2d"), os.path.join(HERE, "..", "i2b"), os.path.join(HERE, "..", "i2a")):
    sys.path.insert(0, d)
import i3a_compare  # noqa: E402
import test_i2a_compare as base  # noqa: E402
import test_i2k_compare as previous  # noqa: E402


def main():
    root = tempfile.mkdtemp()
    try:
        ref = base.run(root, "ref")
        cases = [
            ("prepared-locus provenance", True, dict(status=base.STATUS + "prepared_locus\t1\nprepared_locus_id\t3\n"
                                                              "prepared_locus_key_sha256\tabc\n")),
            ("a near-miss key", False, dict(status=base.STATUS + "prepared_locus_digest\tabc\n")),
            ("a model key", False, dict(status=base.STATUS + "model_fingerprint\tchanged\n")),
            ("a scientific key", False, dict(status=base.STATUS.replace("COMPLETE", "TIED"))),
        ]
        for i, (label, expect, kw) in enumerate(cases):
            ok = i3a_compare.compare(ref, base.run(root, f"c{i}", **kw))[0]
            assert ok == expect, f"{label}: compare returned {ok}, expected {expect}"
        print(f"i3a_compare: {len(cases)} controls as expected")
    finally:
        shutil.rmtree(root)
    previous.main()


if __name__ == "__main__":
    main()
