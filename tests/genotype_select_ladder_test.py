#!/usr/bin/env python3
"""Gate for scripts/select_genotype_ladder.py (review round 7).

    tests/genotype_select_ladder_test.py <path to select_genotype_ladder.py>

Two synthetic development ladders, written to a temporary directory, with identical exact calls everywhere:

* the declared phase audit is unreachable at state limit 256 and reachable above it: the frozen setting must be
  state limit 1024, K = 8 -- genotype stability alone would have chosen 256;
* the same ladder with the audit reachable at 256: the selection must fall to 256, so the audit clause, not the
  calls, is what moved it.

A third ladder where no setting reaches its audit must select nothing.
"""

import os
import subprocess
import sys
import tempfile

DONORS = ["D1", "D2"]
LIMITS = [256, 1024, 4096]
KS = [0, 8, 16]
COLUMNS = ["donor", "state_limit", "K", "status", "marker_converged", "panel_baseline_complete", "winner_hash",
           "target", "phase_support", "phase_audit_status", "unordered_distance", "excess"]


def write_ladder(root, audit_reachable_at):
    with open(os.path.join(root, "ladder.tsv"), "w") as ladder:
        ladder.write("\t".join(COLUMNS) + "\n")
        for donor in DONORS:
            os.makedirs(os.path.join(root, donor, "lzo"), exist_ok=True)
            for limit in LIMITS:
                for k in KS:
                    audit = "COMPLETE" if limit in audit_reachable_at else "MOSAIC_TARGET_UNREACHED"
                    phase = "PRIOR_ONLY" if limit in audit_reachable_at else "UNRESOLVED"
                    row = [donor, str(limit), str(k), "COMPLETE", "true", "true", "winner-" + donor, "true",
                           phase, audit, "0", "0"]
                    ladder.write("\t".join(row) + "\n")
                    with open(os.path.join(root, donor, "lzo", "S%d.K%d.block_genotypes.tsv" % (limit, k)), "w") as calls:
                        calls.write("block_id\tkind\treportable\tallele1\tallele2\n")
                        calls.write("sv:1\tSV\ttrue\t3\t1\nctx:2\tcontext\tfalse\t0\t0\n")


def select(selector, audit_reachable_at):
    with tempfile.TemporaryDirectory() as root:
        write_ladder(root, audit_reachable_at)
        result = subprocess.run([sys.executable, selector, root], capture_output=True, text=True, check=True)
    return result.stdout


def arm(output, name):
    for line in output.splitlines():
        if line.startswith(name + "\t"):
            return line
    return None


def require(condition, message):
    if not condition:
        sys.exit("genotype ladder selection: FAIL: " + message)


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    selector = sys.argv[1]

    audit_above_256 = select(selector, {1024, 4096})
    mosaic = arm(audit_above_256, "mosaic_arm")
    require(mosaic is not None and "state_limit=1024\tK=8\t" in mosaic,
            "a setting stable at 256 but audit-incomplete there must freeze 1024, got: %s" % mosaic)
    require("rejected\tstate_limit=256\tK=8\t" in audit_above_256 and "audit_reachable_for_every_donor=false" in audit_above_256,
            "the audit-incomplete 256 setting must be reported as rejected")
    panel = arm(audit_above_256, "panel_arm")
    require(panel is not None and "state_limit=1024\tK=0" in panel, "the panel arm must share the mosaic arm's limit")

    audit_everywhere = select(selector, {256, 1024, 4096})
    mosaic = arm(audit_everywhere, "mosaic_arm")
    require(mosaic is not None and "state_limit=256\tK=8\t" in mosaic,
            "with the audit reachable at 256 the same calls must select 256, got: %s" % mosaic)

    audit_nowhere = select(selector, set())
    require(arm(audit_nowhere, "mosaic_arm") is None and arm(audit_nowhere, "selected") is not None,
            "a ladder whose audit is never reachable must select nothing")
    print("genotype ladder selection: PASS (audit-incomplete 256 -> 1024; audit at 256 -> 256; no audit -> none)")


if __name__ == "__main__":
    main()
