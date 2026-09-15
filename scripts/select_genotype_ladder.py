#!/usr/bin/env python3
"""Select the mosaic genotyper's (state limit, K) from a development LZO ladder.

    scripts/select_genotype_ladder.py <ladder-dir>

Reads <ladder-dir>/ladder.tsv, written by tests/genotype_ladder_c4.sh, and each run's block calls.
The rule is the review-round-5 one, and it never looks at truth distance:

* For each donor the reference is the largest setting: the largest state limit, and the largest K.
* A setting (S, K) is stable for a donor when every setting (S', K') with S' >= S and K' >= K completed
  with the same winner haplotypes and the same unordered block calls as the reference.
* The selected setting is the smallest one -- smallest state limit, then smallest K -- that is stable
  for every donor. Settings are compared only within the grid that every donor ran.
* A donor whose two largest state limits disagree at the largest K is listed for extension to the next
  state limit; no setting is selected while any donor is unresolved.
"""

import csv
import hashlib
import os
import sys


def block_call_signature(run_prefix):
    path = run_prefix + ".block_genotypes.tsv"
    if not os.path.exists(path):
        return "NA"
    calls = []
    with open(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            pair = sorted((int(row["allele1"]), int(row["allele2"])))
            calls.append("%s:%d/%d" % (row["block_id"], pair[0], pair[1]))
    return hashlib.sha256("\n".join(calls).encode()).hexdigest()[:12]


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    ladder_dir = sys.argv[1]
    runs = {}
    with open(os.path.join(ladder_dir, "ladder.tsv")) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = (row["donor"], int(row["state_limit"]), int(row["K"]))
            prefix = os.path.join(ladder_dir, row["donor"], "lzo", "S%s.K%s" % (row["state_limit"], row["K"]))
            runs[key] = {
                "status": row["status"],
                "winner": row["winner_hash"],
                "blocks": block_call_signature(prefix),
                "distance": row["unordered_distance"],
                "excess": row["excess"],
                "phase": row["phase_support"],
            }

    donors = sorted({d for d, _, _ in runs})
    limits = sorted({s for _, s, _ in runs})
    ks = sorted({k for _, _, k in runs})
    grid = [(s, k) for s in limits for k in ks]
    missing = [(d, s, k) for d in donors for s, k in grid if (d, s, k) not in runs]
    if missing:
        sys.exit("ladder is incomplete; missing runs: %s" % missing[:10])

    def agrees(donor, a, b):
        x, y = runs[(donor,) + a], runs[(donor,) + b]
        return (x["status"] == "COMPLETE" and y["status"] == "COMPLETE"
                and x["winner"] == y["winner"] and x["blocks"] == y["blocks"])

    reference = (limits[-1], ks[-1])
    stable = {}
    unresolved = []
    print("donor\tstate_limit\tK\tstatus\twinner\tblock_calls\tphase\tdistance\texcess\tstable")
    for donor in donors:
        if len(limits) >= 2 and not agrees(donor, (limits[-2], ks[-1]), reference):
            unresolved.append(donor)
        for s, k in grid:
            above = [(s2, k2) for s2, k2 in grid if s2 >= s and k2 >= k]
            stable[(donor, s, k)] = all(agrees(donor, setting, reference) for setting in above)
            run = runs[(donor, s, k)]
            print("\t".join(str(v) for v in (donor, s, k, run["status"], run["winner"], run["blocks"],
                                             run["phase"], run["distance"], run["excess"],
                                             "true" if stable[(donor, s, k)] else "false")))

    print()
    if unresolved:
        print("extend_to_next_state_limit\t%s" % ",".join(unresolved))
        print("selected\tnone: the two largest state limits disagree for these donors")
        return
    for s, k in grid:
        if all(stable[(donor, s, k)] for donor in donors):
            print("selected\tstate_limit=%d\tK=%d" % (s, k))
            return
    print("selected\tnone: no setting is stable for every donor")


if __name__ == "__main__":
    main()
