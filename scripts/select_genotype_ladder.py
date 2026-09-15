#!/usr/bin/env python3
"""Select the mosaic genotyper's frozen arms from a development LZO ladder.

    scripts/select_genotype_ladder.py <ladder-dir>

Reads <ladder-dir>/ladder.tsv, written by tests/genotype_ladder_c4.sh, and each run's block calls. The rule is the
review-round-7 one; it never reads truth distance.

Every panel diplotype is always a candidate, and K counts only the off-panel mosaics added to it.

* Run eligibility: status COMPLETE, depth fit converged, panel baseline complete.
* For each donor the reference is the largest setting: the largest state limit at the largest K.
* The mosaic arm is the smallest setting (state limit first, then K) with K > 0 that satisfies, for every
  development donor, all three frozen-setting requirements (review round 7):
    - primary calls stable: every setting at or above it in both dimensions (K > 0) is eligible and gives the
      reference winner haplotypes and unordered block calls;
    - primary K reached;
    - the declared phase audit reachable (ranked, not MOSAIC_TARGET_UNREACHED) at that state limit.
  If no setting satisfies all three, nothing is selected. Audit reachability is required so that a frozen
  setting can distinguish a negative phase result (PRIOR_ONLY) from an undecided one (UNRESOLVED).
* The panel arm is K = 0 at the mosaic arm's state limit. Whether it gives the reference call is reported per
  donor; it is not a selection criterion.
* A donor whose two largest state limits disagree at the largest K is listed for extension to the next state
  limit, and no arm is selected while any donor is unresolved.
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
                "eligible": row["status"] == "COMPLETE" and row["marker_converged"] == "true"
                            and row["panel_baseline_complete"] == "true",
                "status": row["status"],
                "winner": row["winner_hash"],
                "blocks": block_call_signature(prefix),
                "reached": row["target"] == "true",
                "phase": row["phase_support"],
                "audit_reachable": row.get("phase_audit_status", "NA") in ("COMPLETE", "RESCORE_INCOMPLETE"),
                "distance": row["unordered_distance"],
                "excess": row["excess"],
            }

    donors = sorted({d for d, _, _ in runs})
    limits = sorted({s for _, s, _ in runs})
    ks = sorted({k for _, _, k in runs})
    positive = [k for k in ks if k > 0]
    grid = [(s, k) for s in limits for k in ks]
    missing = [(d, s, k) for d in donors for s, k in grid if (d, s, k) not in runs]
    if missing:
        sys.exit("ladder is incomplete; missing runs: %s" % missing[:10])
    if 0 not in ks or not positive:
        sys.exit("the ladder must contain K = 0 and at least one positive K")

    reference = (limits[-1], ks[-1])

    def same_call(donor, setting):
        x, y = runs[(donor,) + setting], runs[(donor,) + reference]
        return x["eligible"] and y["eligible"] and x["winner"] == y["winner"] and x["blocks"] == y["blocks"]

    print("donor\tstate_limit\tK\tstatus\teligible\twinner\tblock_calls\tK_reached\taudit_reachable\tphase\t"
          "distance\texcess\tsame_as_reference\tmosaic_stable")
    stable = {}
    unresolved = []
    for donor in donors:
        if len(limits) >= 2 and not same_call(donor, (limits[-2], ks[-1])):
            unresolved.append(donor)
        for s, k in grid:
            if k > 0:
                above = [(s2, k2) for s2, k2 in grid if s2 >= s and k2 >= k and k2 > 0]
                stable[(donor, s, k)] = all(same_call(donor, setting) for setting in above)
            run = runs[(donor, s, k)]
            print("\t".join(str(v) for v in (
                donor, s, k, run["status"], str(run["eligible"]).lower(), run["winner"], run["blocks"],
                str(run["reached"]).lower(), str(run["audit_reachable"]).lower(), run["phase"], run["distance"],
                run["excess"],
                str(same_call(donor, (s, k))).lower(),
                "NA" if k == 0 else str(stable[(donor, s, k)]).lower())))

    print()
    if unresolved:
        print("extend_to_next_state_limit\t%s" % ",".join(unresolved))
        print("selected\tnone: the two largest state limits disagree for these donors")
        return
    stable_settings = [(s, k) for s in limits for k in positive if all(stable[(d, s, k)] for d in donors)]
    qualifying = [c for c in stable_settings
                  if all(runs[(d,) + c]["reached"] and runs[(d,) + c]["audit_reachable"] for d in donors)]
    for s, k in stable_settings:
        if (s, k) not in qualifying:
            print("rejected\tstate_limit=%d\tK=%d\tK_reached_for_every_donor=%s\taudit_reachable_for_every_donor=%s"
                  % (s, k, str(all(runs[(d, s, k)]["reached"] for d in donors)).lower(),
                     str(all(runs[(d, s, k)]["audit_reachable"] for d in donors)).lower()))
    if not qualifying:
        print("selected\tnone: no positive-K setting is stable, reaches K and reaches the audit for every donor")
        return
    mosaic = qualifying[0]
    print("mosaic_arm\tstate_limit=%d\tK=%d\tK_reached_for_every_donor=true\taudit_reachable_for_every_donor=true"
          % mosaic)
    panel = (mosaic[0], 0)
    agreeing = [d for d in donors if same_call(d, panel)]
    print("panel_arm\tstate_limit=%d\tK=0\tsame_call_as_reference=%d/%d%s"
          % (panel[0], len(agreeing), len(donors),
             "" if len(agreeing) == len(donors) else "\tdiffers=" + ",".join(d for d in donors if d not in agreeing)))


if __name__ == "__main__":
    main()
