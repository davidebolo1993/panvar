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

echo "=== the ordinary-LZO source is the canonical one, and reproduces its expected counts ==="
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ORDINARY="${ORDINARY:-$REPO/tests/results/genotype_cohort_c4_r6_5937384}"
ORDINARY_COVERAGE="${ORDINARY_COVERAGE:-$REPO/tests/results/genotype_coverage_c4_held_out_abb35a0}"
SUPERSEDED_COHORT="$REPO/tests/results/genotype_cohort_c4_b1684ae"
if [ ! -d "$ORDINARY" ]; then
  echo "  skipped: $ORDINARY is not present (results directories are local artifacts)"
else
  python3 - "$ORDINARY" "$ORDINARY_COVERAGE" "$SUPERSEDED_COHORT" "${CONSOLIDATED:-}" <<'PY' || fail=$((fail + 1))
import csv, os, sys

ordinary, coverage_dir, superseded, consolidated = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
bad = 0

def kv(path):
    out = {}
    if os.path.exists(path):
        for line in open(path):
            parts = line.rstrip("\n").split("\t")
            if len(parts) == 2:
                out[parts[0]] = parts[1]
    return out

def tsv(path):
    if not os.path.exists(path):
        return []
    return list(csv.DictReader(open(path), delimiter="\t"))

# 1. The 30x ordinary source must be round 6 (5937384), never the superseded round-5 cohort.
commit = kv(os.path.join(ordinary, "signature.tsv")).get("git_commit", "")
if not commit.startswith("5937384"):
    print(f"  FAIL  ordinary source commit is {commit or 'missing'}, expected 5937384..."); bad += 1
if commit.startswith("b1684ae"):
    print("  FAIL  the SUPERSEDED round-5 cohort is being used as the ordinary source"); bad += 1

# 2. Its canonical counts, per arm, recomputed from the per-donor rows.
rows = tsv(os.path.join(ordinary, "cohort.tsv"))
if not rows:
    print("  FAIL  ordinary cohort.tsv is missing or empty"); bad += 1
arms = sorted({r["arm"] for r in rows})
if arms != ["mosaic", "panel"]:
    print(f"  FAIL  ordinary arms are {arms}, expected ['mosaic', 'panel']"); bad += 1
for arm in arms:
    subset = [r for r in rows if r["arm"] == arm]
    complete = [r for r in subset if r["status"] == "COMPLETE"]
    no_call = [r for r in subset if r["run_no_call"] == "true"]
    fully = [r for r in complete
             if int(r["unordered_distance"]) == 0 and int(r["blocks_exact"]) == int(r["blocks"])]
    got = (len(subset), len(complete), len(no_call), len(fully),
           sum(int(r["blocks_exact"]) for r in subset), sum(int(r["blocks"]) for r in subset),
           sum(int(r["sv_blocks_exact"]) for r in subset), sum(int(r["sv_blocks"]) for r in subset),
           sum(int(r["unordered_distance"]) for r in subset))
    want = (56, 56, 0, 56, 616, 616, 280, 280, 0)
    if got != want:
        print(f"  FAIL  ordinary {arm} arm is {got}, expected {want}"); bad += 1
    else:
        print(f"  ok    ordinary {arm}: 56/56 called, 56/56 fully exact, 616/616 blocks, "
              f"280/280 SV, distance 0")

# 3. The nested-coverage ordinary result, recomputed per arm and coverage.
expected = {"2": (55, 615), "5": (55, 615), "10": (56, 616),
            "15": (56, 616), "20": (56, 616), "30": (56, 616)}
runs = tsv(os.path.join(coverage_dir, "runs.tsv"))
if not runs:
    print(f"  FAIL  {coverage_dir}/runs.tsv missing"); bad += 1
else:
    cover_commit = kv(os.path.join(coverage_dir, "signature.tsv")).get("git_commit", "")
    if not cover_commit.startswith("abb35a0"):
        print(f"  FAIL  nested-coverage commit is {cover_commit or 'missing'}, expected abb35a0..."); bad += 1
    for arm in sorted({r["arm"] for r in runs}):
        for cov, (want_donors, want_blocks) in expected.items():
            subset = [r for r in runs if r["arm"] == arm and r["coverage"] == cov]
            if len(subset) != 56:
                print(f"  FAIL  nested {arm} {cov}x has {len(subset)} runs, expected 56"); bad += 1
                continue
            no_calls = sum(1 for r in subset if r["run_no_call"] == "true")
            fully = sum(1 for r in subset
                        if int(r["unordered_distance"]) == 0 and int(r["blocks_exact"]) == 11)
            blocks_exact = sum(int(r["blocks_exact"]) for r in subset)
            if (no_calls, fully, blocks_exact) != (0, want_donors, want_blocks):
                print(f"  FAIL  nested {arm} {cov}x: no_calls={no_calls} donors={fully} "
                      f"blocks={blocks_exact}, expected 0/{want_donors}/{want_blocks}"); bad += 1
    if not bad:
        print("  ok    nested coverage: 55/56 and 615/616 at 2x and 5x, 56/56 and 616/616 at "
              "10x-30x, zero no-calls, both arms")

# 4. The round-5 cohort must be marked superseded wherever it still exists.
if os.path.isdir(superseded) and not os.path.exists(os.path.join(superseded, "SUPERSEDED.txt")):
    print("  FAIL  the round-5 cohort is present but not marked SUPERSEDED"); bad += 1

# 5. A consolidated package must record the canonical source, and must not cite the superseded one.
if consolidated and os.path.isdir(consolidated):
    ordinary_rows = tsv(os.path.join(consolidated, "ordinary_lzo.tsv"))
    commits = {r["value"] for r in ordinary_rows if r["quantity"] == "source_commit"}
    if not commits or any(not c.startswith("5937384") for c in commits):
        print(f"  FAIL  consolidated ordinary_lzo.tsv cites {commits or 'no'} source commit"); bad += 1
    if not os.path.exists(os.path.join(consolidated, "ordinary_lzo_coverage.tsv")):
        print("  FAIL  consolidated package has no ordinary_lzo_coverage.tsv"); bad += 1
    for name in ("ordinary_lzo.tsv", "provenance.tsv"):
        path = os.path.join(consolidated, name)
        if os.path.exists(path) and "b1684ae" in open(path).read():
            print(f"  FAIL  {name} still cites the superseded b1684ae cohort"); bad += 1

if bad:
    print(f"  {bad} ordinary-LZO provenance failure(s)"); sys.exit(1)
print("  ok    canonical ordinary-LZO provenance and counts")
PY
fi

if [ "$fail" -ne 0 ]; then echo "genotype_lzo_consolidation: $fail gate(s) FAILED"; exit 1; fi
echo "genotype_lzo_consolidation: all gates passed"
