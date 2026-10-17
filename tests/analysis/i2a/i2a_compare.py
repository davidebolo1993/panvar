#!/usr/bin/env python3
"""I2-A: are the reference and the experimental indexed union backend scientifically identical on one cell?

perf_compare.compare (byte-identical genotype/locus/block tables, candidates and scores, equivalence sets, FASTA, dosage
tables and library sensitivity; status and dosage status key by key; settings apart from budget keys; the repair
trajectory's scientific columns), with exactly one extension of its allowed (DIAGNOSTIC) class:
  * the indexed backend's provenance block: indexed_union_*, threads_indexed_union_*, union_scoring_backend* (status) and
    the two experimental profile rows union_scoring_backend / indexed_union_store_bytes (settings);
  * accepted_mate_placements: executed work (the indexed route counts placements accepted per discovered unit, the
    reference per haplotype), like the other executed-work counters perf_compare already allows.
A file present in one run only fails, unless perf_compare never compares it (NOT_COMPARED, e.g. a runner's own
.run.tsv record). Nothing else is allowed to differ. Usage: i2a_compare.py <run dir A> <run dir B> [prefix, default call]
"""
import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "perf"))
import perf_compare  # noqa: E402

INDEXED = r"|^indexed_union_|^threads_indexed_union_|^union_scoring_backend|^accepted_mate_placements$"
perf_compare.DIAGNOSTIC = re.compile(perf_compare.DIAGNOSTIC.pattern + INDEXED)


def compare(dir_a, dir_b, name="call"):
    """perf_compare.compare, except that a file present on one side only is not a problem when perf_compare itself
    never compares it (NOT_COMPARED: a runner's own record such as .run.tsv, the call log, the repair memory ledger)."""
    ok, problems, diagnostics, compared = perf_compare.compare(dir_a, dir_b, name)
    kept = []
    for p in problems:
        if p.startswith("file sets differ: "):
            files = re.findall(r"'([^']+)'", p)
            if files and all(f.endswith(perf_compare.NOT_COMPARED) for f in files):
                continue
        kept.append(p)
    return not kept, kept, diagnostics, compared

if __name__ == "__main__":
    ok, problems, diagnostics, compared = compare(sys.argv[1], sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else "call")
    print(("IDENTICAL" if ok else "DIFFERENT") + f": {compared} files compared; allowed keys that changed: "
          f"{len(diagnostics)}" + ("" if ok else "; " + " | ".join(problems[:12])))
    sys.exit(0 if ok else 1)
