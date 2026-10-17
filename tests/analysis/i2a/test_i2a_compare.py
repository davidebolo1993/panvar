#!/usr/bin/env python3
"""Negative controls for i2a_compare's single extension: allowed keys pass, everything else still fails."""
import os
import shutil
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import i2a_compare  # noqa: E402

STATUS = "key\tvalue\nstatus\tCOMPLETE\ncall_class\tstrict_call\nphysical_origins\t10\naccepted_mate_placements\t5\n"
SETTINGS = "key\tvalue\tsource\nkmer_size\t21\tdefault\n"


def run(root, name, status=STATUS, settings=SETTINGS, genotypes="g\t1\n", extra=None):
    d = os.path.join(root, name)
    os.makedirs(d)
    if extra:
        open(os.path.join(d, extra), "w").write("x\n")
    for suffix, text in ((".mosaic_status.tsv", status), (".settings.tsv", settings), (".genotypes.tsv", genotypes)):
        open(os.path.join(d, "call" + suffix), "w").write(text)
    return d


def main():
    root = tempfile.mkdtemp()
    try:
        base = run(root, "base")
        cases = [
            ("indexed provenance keys", True, dict(status=STATUS + "union_scoring_backend\tindexed_fine_experimental\n"
                                                    "indexed_union_store_charged_bytes\t9\nthreads_indexed_union_pairing\t4\n")),
            ("executed accepted placements", True, dict(status=STATUS.replace("accepted_mate_placements\t5", "accepted_mate_placements\t7"))),
            ("experimental settings rows", True, dict(settings=SETTINGS + "union_scoring_backend\tindexed_fine_experimental\tprofile x\n"
                                                      "indexed_union_store_bytes\t1\tprofile x\n")),
            ("a scientific status key", False, dict(status=STATUS.replace("strict_call", "genotype_equivalent"))),
            ("a scientific key the extension does not name", False, dict(status=STATUS + "indexed_call\tx\n")),
            ("a genotype table", False, dict(genotypes="g\t2\n")),
            ("a non-budget setting", False, dict(settings=SETTINGS.replace("21", "23"))),
            ("a runner record on one side only", True, dict(extra="call.run.tsv")),
            ("an output on one side only", False, dict(extra="call.equivalence_set.tsv")),
        ]
        for i, (label, expect, kw) in enumerate(cases):
            ok = i2a_compare.compare(base, run(root, f"c{i}", **kw))[0]
            assert ok == expect, f"{label}: compare returned {ok}, expected {expect}"
        print(f"i2a_compare: {len(cases)} controls as expected")
    finally:
        shutil.rmtree(root)


if __name__ == "__main__":
    main()
