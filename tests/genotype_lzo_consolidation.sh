#!/usr/bin/env bash
# Gate: the C4 LZO consolidation's two load-bearing computations.
#
#   tests/genotype_lzo_consolidation.sh <scratch-dir> [consolidated-dir]
#
# 1. The phase-invariant diploid distance itself. The truth pair is UNORDERED, so a predicted pair must
#    score the same however it is written down: D = min(e(T1,P1)+e(T2,P2), e(T1,P2)+e(T2,P1)). If that
#    symmetry ever broke, every per-block distance and every "exact" verdict in the report would depend
#    on which haplotype the caller happened to emit first. Tested on known answers, not on the data.
#
# 2. The tied-winner consensus invariants, asserted on the real CONDITIONAL_CONSENSUS_DIAGNOSTIC tables
#    produced by `genotype_probe tied-consensus`. These are arithmetic identities: if they fail, the
#    diagnostic is inconsistent with itself and nothing derived from it can be reported.
#
# When a consolidated directory is given, the gate also checks that the plot source tables reproduce the
# run-level table's own counts -- a plot that disagrees with its raw TSV is the failure this catches.
#
# Fast and read-only: it runs no caller, simulates nothing, and needs no graph.
set -uo pipefail

WORK="${1:?usage: genotype_lzo_consolidation.sh <scratch-dir> [consolidated-dir]}"
CONSOLIDATED="${2:-}"
rm -rf "$WORK"; mkdir -p "$WORK"
fail=0

echo "=== phase-invariant diploid distance ==="
python3 - <<'PY' || fail=$((fail + 1))
import sys

def lev(a, b):
    if a == b:
        return 0
    previous = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        current = [i]
        for j, cb in enumerate(b, 1):
            current.append(min(previous[j] + 1, current[j - 1] + 1, previous[j - 1] + (ca != cb)))
        previous = current
    return previous[-1]

def diploid(t1, t2, p1, p2):
    return min(lev(t1, p1) + lev(t2, p2), lev(t1, p2) + lev(t2, p1))

bad = 0
cases = [
    # truth1, truth2, pred1, pred2, expected
    ("ACGT", "ACGT", "ACGT", "ACGT", 0),            # identical homozygote
    ("ACGT", "TTTT", "TTTT", "ACGT", 0),            # correct but written in the other order
    ("ACGT", "TTTT", "ACGT", "TTTT", 0),            # correct in this order too
    ("ACGT", "ACGT", "ACGA", "ACGT", 1),            # one substitution on one homologue
    ("ACGT", "ACGT", "ACGA", "ACGA", 2),            # one on each
    ("AAAA", "CCCC", "AAA", "CCCC", 1),             # a deletion is one edit, not four
    ("AAAA", "CCCC", "CCCC", "AAAA", 0),            # swap invariance with distinct homologues
    ("", "ACGT", "ACGT", "", 0),                    # empty allele, swapped
]
for t1, t2, p1, p2, want in cases:
    got = diploid(t1, t2, p1, p2)
    if got != want:
        print(f"  FAIL  D({t1!r},{t2!r};{p1!r},{p2!r}) = {got}, want {want}")
        bad += 1

# Swap invariance must hold for arbitrary inputs, not only the cases above.
import itertools
alphabet = "AC"
for t1 in ("".join(x) for x in itertools.product(alphabet, repeat=3)):
    for p1 in ("".join(x) for x in itertools.product(alphabet, repeat=3)):
        t2, p2 = "AAC", "CCA"
        if diploid(t1, t2, p1, p2) != diploid(t2, t1, p1, p2):
            print(f"  FAIL  truth order changed the distance at {t1}/{p1}")
            bad += 1
        if diploid(t1, t2, p1, p2) != diploid(t1, t2, p2, p1):
            print(f"  FAIL  predicted order changed the distance at {t1}/{p1}")
            bad += 1

# The diploid distance can never exceed scoring each homologue in the written order.
for t1, t2, p1, p2, _ in cases:
    if diploid(t1, t2, p1, p2) > lev(t1, p1) + lev(t2, p2):
        print("  FAIL  min() returned more than the direct assignment")
        bad += 1

if bad:
    print(f"  {bad} phase-invariance failure(s)")
    sys.exit(1)
print(f"  ok    {len(cases)} known answers, swap invariance over 64 pairs, min() bound")
PY

echo "=== tied-winner consensus invariants (CONDITIONAL_CONSENSUS_DIAGNOSTIC tables) ==="
HELDOUT="${HELDOUT:-}"
if [ -z "$HELDOUT" ] || [ ! -d "$HELDOUT" ]; then
  echo "  skipped: set HELDOUT=<run dir> to assert against real diagnostic tables"
else
  python3 - "$HELDOUT" <<'PY' || fail=$((fail + 1))
import glob, os, sys

root = sys.argv[1]
summaries = sorted(glob.glob(os.path.join(root, "*", "tied.*.summary.tsv")))
if not summaries:
    print("  FAIL  no tied-consensus summaries found"); sys.exit(1)

bad = 0
for path in summaries:
    kv = {}
    for line in open(path):
        parts = line.rstrip("\n").split("\t")
        if len(parts) == 2:
            kv[parts[0]] = parts[1]
    name = os.path.relpath(path, root)
    num = lambda k: int(kv[k])

    if kv.get("label") != "CONDITIONAL_CONSENSUS_DIAGNOSTIC":
        print(f"  FAIL  {name}: not labelled as a diagnostic"); bad += 1
    # Subtotals must reconstruct the total exactly.
    if num("consensus_sv_subtotal") + num("consensus_context_subtotal") != num("consensus_total_distance"):
        print(f"  FAIL  {name}: SV + context != total"); bad += 1
    # Counts are bounded by the chain.
    if not (0 <= num("exact_consensus_blocks") <= num("unanimous_blocks") <= num("blocks")):
        print(f"  FAIL  {name}: block counts out of order"); bad += 1
    # Order statistics.
    lo, mid, hi = num("min_whole_locus_distance"), num("median_whole_locus_distance"), num("max_whole_locus_distance")
    if not (lo <= mid <= hi):
        print(f"  FAIL  {name}: min/median/max not ordered"); bad += 1
    # "Exactly correct candidate exists" must mean exactly "minimum distance is zero".
    exact_flag = kv["any_tied_candidate_exact"] == "true"
    if exact_flag != (lo == 0):
        print(f"  FAIL  {name}: any_tied_candidate_exact={exact_flag} but min distance={lo}"); bad += 1
    if (kv["true_pair_in_tied_set"] == "true") != exact_flag:
        print(f"  FAIL  {name}: true-pair-present disagrees with exact-candidate"); bad += 1
    # Full unanimity means every tied winner has the same block genotypes, so they cannot differ in
    # distance, and the consensus IS that shared genotype.
    if num("unanimous_blocks") == num("blocks"):
        if lo != hi:
            print(f"  FAIL  {name}: unanimous on every block but candidate distances differ"); bad += 1
        if num("consensus_total_distance") != lo:
            print(f"  FAIL  {name}: consensus distance != the unanimous candidates' distance"); bad += 1
    if num("consensus_max_block_distance") > num("consensus_total_distance"):
        print(f"  FAIL  {name}: max block distance exceeds the total"); bad += 1

    # The per-candidate table must agree with the summary's order statistics.
    candidates = path.replace(".summary.tsv", ".candidates.tsv")
    rows = [line.rstrip("\n").split("\t") for line in open(candidates)][1:]
    totals = sorted(int(r[1]) for r in rows)
    if len(rows) != num("winner_set_size"):
        print(f"  FAIL  {name}: candidate rows {len(rows)} != winner_set_size {num('winner_set_size')}"); bad += 1
    if totals and (totals[0] != lo or totals[-1] != hi):
        print(f"  FAIL  {name}: candidate table min/max disagree with the summary"); bad += 1
    for row in rows:
        if int(row[3]) + int(row[4]) != int(row[1]):
            print(f"  FAIL  {name}: candidate {row[0]} SV+context != total"); bad += 1
        if (row[5] == "true") != (int(row[1]) == 0):
            print(f"  FAIL  {name}: candidate {row[0]} exact flag disagrees with its distance"); bad += 1

if bad:
    print(f"  {bad} invariant failure(s) over {len(summaries)} diagnostics"); sys.exit(1)
print(f"  ok    {len(summaries)} diagnostics: subtotals, bounds, order statistics, unanimity, candidate tables")
PY
fi

echo "=== plot tables reproduce the run-level table ==="
if [ -z "$CONSOLIDATED" ] || [ ! -d "$CONSOLIDATED" ]; then
  echo "  skipped: pass a consolidated directory to cross-check the plot sources"
else
  python3 - "$CONSOLIDATED" <<'PY' || fail=$((fail + 1))
import csv, os, sys

root = sys.argv[1]
def rows(name):
    path = os.path.join(root, name)
    if not os.path.exists(path):
        print(f"  FAIL  missing {name}"); sys.exit(1)
    return list(csv.DictReader(open(path), delimiter="\t"))

runs = rows("runs.tsv")
bad = 0
if len(runs) != 80:
    print(f"  FAIL  runs.tsv has {len(runs)} rows, expected 80"); bad += 1
pairs = {(r["plan"], r["coverage"]) for r in runs}
if len(pairs) != 40:
    print(f"  FAIL  runs.tsv covers {len(pairs)} plan/coverage pairs, expected 40"); bad += 1

# Outcome plot source must partition the runs exactly.
outcomes = rows("plot_class_coverage.tsv")
total = sum(int(r["runs"]) for r in outcomes)
if total != len(runs):
    print(f"  FAIL  outcome table counts {total} runs, runs.tsv has {len(runs)}"); bad += 1

# Paired plot must carry one row per plan/coverage pair, refusals included.
paired = rows("plot_paired_k.tsv")
if len(paired) != 40:
    print(f"  FAIL  paired table has {len(paired)} rows, expected 40"); bad += 1
refused = sum(1 for r in paired if r["k128_outcome"] == "no_call")
declared = sum(1 for r in runs if r["arm"] == "K128" and r["status"] != "COMPLETE")
if refused != declared:
    print(f"  FAIL  paired table shows {refused} refusals, runs.tsv shows {declared}"); bad += 1

# Per-block table: strict calls only, 11 blocks each.
blocks = rows("blocks.tsv")
called = {(r["plan"], r["coverage"], r["arm"]) for r in blocks if r["source"] == "strict_call"}
complete = {(r["plan"], r["coverage"], r["arm"]) for r in runs if r["status"] == "COMPLETE"}
if called != complete:
    print(f"  FAIL  per-block strict rows cover {len(called)} arms, {len(complete)} are COMPLETE"); bad += 1
for key in called:
    n = sum(1 for r in blocks if r["source"] == "strict_call" and (r["plan"], r["coverage"], r["arm"]) == key)
    if n != 11:
        print(f"  FAIL  {key} has {n} block rows, expected 11"); bad += 1
        break

if bad:
    print(f"  {bad} cross-check failure(s)"); sys.exit(1)
print(f"  ok    80 run rows, 40 pairs, outcome/paired/per-block tables agree with runs.tsv")
PY
fi

if [ "$fail" -ne 0 ]; then echo "genotype_lzo_consolidation: $fail gate(s) FAILED"; exit 1; fi
echo "genotype_lzo_consolidation: all gates passed"
