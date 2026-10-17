#!/usr/bin/env python3
"""Negative controls for i2b_compare's extension (and the I2-A controls it inherits)."""
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "i2a"))
import i2b_compare  # noqa: E402
import test_i2a_compare as base  # noqa: E402

ROUNDS = "round\tstart_total\tmoved\tseconds\tpiece_hits\n1\t-10.5\tfalse\t3.2\t100\n"


def main():
    root = tempfile.mkdtemp()
    try:
        def run(name, rounds=ROUNDS, **kw):
            d = base.run(root, name, **kw)
            open(os.path.join(d, "call.repair_rounds.tsv"), "w").write(rounds)
            return d
        ref = run("ref")
        cases = [
            ("adaptive window block", True, dict(status=base.STATUS + "adaptive_window_discovery\tfine_experimental\n"
                                                  "adaptive_fine_units_served\t9\nworker_seconds_adaptive_window_rows\t1.5\n")),
            ("adaptive settings rows", True, dict(settings=base.SETTINGS + "adaptive_window_discovery\tfine_experimental\tprofile x\n"
                                                  "adaptive_fine_store_bytes\t1\tprofile x\n")),
            ("trajectory timing and executed work", True, dict(rounds=ROUNDS.replace("3.2\t100", "1.1\t7"))),
            ("trajectory scientific column", False, dict(rounds=ROUNDS.replace("-10.5", "-10.6"))),
            ("a scientific key outside the block", False, dict(status=base.STATUS + "adaptive_call\tx\n")),
            ("a genotype table", False, dict(genotypes="g\t3\n")),
        ]
        for i, (label, expect, kw) in enumerate(cases):
            ok = i2b_compare.compare(ref, run(f"c{i}", **kw))[0]
            assert ok == expect, f"{label}: compare returned {ok}, expected {expect}"
        print(f"i2b_compare: {len(cases)} controls as expected")
    finally:
        shutil.rmtree(root)
    base.main()


if __name__ == "__main__":
    main()
