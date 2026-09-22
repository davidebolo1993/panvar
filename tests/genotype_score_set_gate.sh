#!/usr/bin/env bash
# Gate: score-set reproduces declared expectations, and REFUSES when one does not.
#
#   tests/genotype_score_set_gate.sh <build dir> <scratch dir>
#
# score-set exists to score a handful of declared diplotypes so a candidate the search never built
# can be compared against candidates it did build. That comparison is only worth anything if the
# already-scored candidates come back with the scores they had before -- otherwise the new number is
# being compared against a differently-configured pass.
#
# So the guard is the whole point, and a guard that has never been seen to refuse is not a guard.
# This checks both directions:
#
#   accept  -- the seed and the two recorded winners, declared with their persisted totals, reproduce
#              and the command reports every one of them.
#   refuse  -- the same run with ONE expectation perturbed by 0.01 nats exits 4 and writes no table.
#
# 0.01 nats is far below any difference this diagnostic is used to interpret (the smallest is 22.8),
# so a guard that catches it catches everything that matters.
#
# tests/results is gitignored, so a clean checkout has nothing to run against. The gate says so.
set -uo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="${1:?usage: genotype_score_set_gate.sh <build dir> <scratch dir>}"
WORK="${2:?usage: genotype_score_set_gate.sh <build dir> <scratch dir>}"
PROBE="$BUILD/genotype_probe"
PKG="$REPO/tests/results/genotype_loo_c4_development_7e37bb4"
REF='HG00096#1#haplotype1-0000024:31848049-32074186'
rm -rf "$WORK"; mkdir -p "$WORK"

fail=0
same() { if [ "$2" = "$3" ]; then echo "  ok: $1 = $2"; else echo "  GATE FAILED: $1 is '$2', expected '$3'"; fail=1; fi; }

[ -x "$PROBE" ] || { echo "genotype_probe not built at $PROBE" >&2; exit 1; }

D="$PKG/inputs/NA20509"
if [ ! -f "$D/c4.loo.gfa.gz" ] || [ ! -f "$PKG/NA20509/call.10.K128.mosaic_status.tsv" ]; then
  echo "skipped: the campaign artifacts under tests/results are absent (gitignored)."
  exit 0
fi

G="$WORK/c4.gfa"; L="$WORK/c4.loo.gfa"; R1="$WORK/r1.fq"; R2="$WORK/r2.fq"
gzip -dc "$PKG/inputs/c4.gfa.gz"  > "$G"  2>/dev/null || { echo "skipped: no full graph"; exit 0; }
gzip -dc "$D/c4.loo.gfa.gz"       > "$L"  2>/dev/null || { echo "skipped: no LOO graph"; exit 0; }
gzip -dc "$D/cov10/reads_1.fq.gz" > "$R1" 2>/dev/null || { echo "skipped: no reads"; exit 0; }
gzip -dc "$D/cov10/reads_2.fq.gz" > "$R2" 2>/dev/null || { echo "skipped: no reads"; exit 0; }

# The three candidates the search really scored, with the totals it really recorded.
cat > "$WORK/good.tsv" <<'EOF'
seed_329	15,10,41,12,3,6,58,11,51,1,99	48,45,6,14,6,0,1,112,29,1,119	-152228.321247
R16_winner_sv7	15,10,41,12,3,6,58,11,51,1,99	48,45,6,14,6,0,1,108,29,1,119	-151954.429577
R32_winner_ctx2	15,10,87,12,3,6,58,11,51,1,99	48,45,37,14,6,0,1,112,29,1,119	-151931.605119
EOF
# Identical, except the R=16 winner is declared 0.01 nats away from what it scores.
sed 's/-151954\.429577/-151954.439577/' "$WORK/good.tsv" > "$WORK/bad.tsv"

run() {
  "$PROBE" score-set --reference "$REF" --donor NA20509 \
    --model-gfa "$L" --truth-gfa "$G" --run "$PKG/NA20509/call.10.K128" \
    --candidates "$1" --reads1 "$R1" --reads2 "$R2" --library "$D/cov10/library.txt" \
    --threads 8 --out "$2" > "$2.log" 2>&1
  echo "$?"
}

echo "== declared expectations that hold are accepted =="
code="$(run "$WORK/good.tsv" "$WORK/accept")"
same "exit code" "$code" "0"
same "expectations reproduced" "$(sed -n 's/.*, \([0-9]*\) expectations reproduced.*/\1/p' "$WORK/accept.log")" "3"
if [ -s "$WORK/accept.scores.tsv" ]; then echo "  ok: a score table was written"; else
  echo "  GATE FAILED: no score table"; fail=1; fi

echo "== one expectation off by 0.01 nats refuses, and writes nothing =="
code="$(run "$WORK/bad.tsv" "$WORK/refuse")"
same "exit code" "$code" "4"
if grep -q "REFUSED" "$WORK/refuse.log"; then echo "  ok: the refusal is reported"; else
  echo "  GATE FAILED: no refusal reported"; fail=1; fi
if [ -f "$WORK/refuse.scores.tsv" ]; then
  echo "  GATE FAILED: a score table was written despite refusing"; fail=1
else echo "  ok: no score table was written"; fi

echo
if [ "$fail" -ne 0 ]; then echo "score-set gate FAILED"; exit 1; fi
echo "score-set gate passed"
