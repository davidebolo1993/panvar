#!/usr/bin/env bash
# Gate: repair-pilot's candidate construction and certification are invariant to candidate ORDER and
# to exchanging a seed's two homologues.
#
#   tests/genotype_repair_pilot_invariance.sh <build dir> <scratch dir>
#
# Both properties are checked INSIDE the pilot, from the scores its single pass already produced, so
# neither costs an extra likelihood call:
#
#   order  -- the same prefix certified back to front must return the same winner, the same candidate
#             count and the same winning score. A ranking that depended on presentation order would
#             make every result an artefact of how the union happened to be assembled.
#   swap   -- rebuilding the repair construction with each seed's two haplotypes exchanged must yield
#             the identical candidate key set, because both homologue assignments are generated for
#             every local pair. If it did not, "repair block b of this seed" would mean two different
#             things depending on which haplotype was written first.
#
# Run at --local-pairs 1, the cheapest configuration that still exercises both paths: one local pair
# per block per seed rather than 32. The scoring cost is dominated by the original candidates, which
# are scored either way, so this is as cheap as the property can honestly be tested on a real model.
#
# tests/results is gitignored, so a clean checkout has no campaign outputs to run against. That is not
# a failure: the gate says what it skipped.
set -uo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="${1:?usage: genotype_repair_pilot_invariance.sh <build dir> <scratch dir>}"
WORK="${2:?usage: genotype_repair_pilot_invariance.sh <build dir> <scratch dir>}"
PROBE="$BUILD/genotype_probe"
PKG="$REPO/tests/results/genotype_loo_c4_development_7e37bb4"
REF='HG00096#1#haplotype1-0000024:31848049-32074186'
rm -rf "$WORK"; mkdir -p "$WORK"

fail=0
same() { if [ "$2" = "$3" ]; then echo "  ok: $1 = $2"; else echo "  GATE FAILED: $1 is '$2', expected '$3'"; fail=1; fi; }

[ -x "$PROBE" ] || { echo "genotype_probe not built at $PROBE" >&2; exit 1; }

D="$PKG/inputs/NA20509"
if [ ! -f "$D/c4.loo.gfa.gz" ] || [ ! -f "$PKG/NA20509/call.10.K128.mosaic_scores.tsv" ]; then
  echo "skipped: the campaign artifacts under tests/results are absent (gitignored)."
  echo "         Invariance is checked where those outputs exist."
  exit 0
fi

G="$WORK/c4.gfa"; L="$WORK/c4.loo.gfa"; R1="$WORK/r1.fq"; R2="$WORK/r2.fq"
gzip -dc "$PKG/inputs/c4.gfa.gz"        > "$G"  2>/dev/null || { echo "skipped: no full graph"; exit 0; }
gzip -dc "$D/c4.loo.gfa.gz"             > "$L"  2>/dev/null || { echo "skipped: no LOO graph"; exit 0; }
gzip -dc "$D/cov10/reads_1.fq.gz"       > "$R1" 2>/dev/null || { echo "skipped: no reads"; exit 0; }
gzip -dc "$D/cov10/reads_2.fq.gz"       > "$R2" 2>/dev/null || { echo "skipped: no reads"; exit 0; }

o="$WORK/invariance"
if ! "$PROBE" repair-pilot --reference "$REF" --donor NA20509 \
      --model-gfa "$L" --truth-gfa "$G" --run "$PKG/NA20509/call.10.K128" \
      --reads1 "$R1" --reads2 "$R2" --library "$D/cov10/library.txt" \
      --local-pairs 1 --prefixes 1 --threads 8 --out "$o" > "$o.log" 2>&1; then
  echo "  GATE FAILED: the pilot did not complete"; head -3 "$o.log" | sed 's/^/      /'; exit 1
fi

field() { awk -F'\t' -v k="$1" '$1==k{print $2}' "$o.summary.tsv"; }

echo "== candidate-order and homologue-swap invariance =="
same "order_invariant" "$(field order_invariant)" "true"
same "swap_invariant"  "$(field swap_invariant)"  "true"

echo "== the reproduction guard still holds at this configuration =="
same "originals_reproduce" "$(field originals_reproduce)" "true"
same "originals_unkeyed"   "$(field originals_unkeyed)"   "0"
same "originals_mismatched" "$(field originals_mismatched)" "0"

echo "== repairs were actually constructed (a vacuous pass proves nothing) =="
kept="$(field repairs_kept)"
if [ "${kept:-0}" -gt 0 ]; then echo "  ok: repairs_kept = $kept"; else
  echo "  GATE FAILED: no repairs were built, so neither property was exercised"; fail=1
fi

echo
if [ "$fail" -ne 0 ]; then echo "repair-pilot invariance gate FAILED"; exit 1; fi
echo "repair-pilot invariance gate passed"
