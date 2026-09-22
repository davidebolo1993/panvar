#!/usr/bin/env bash
# Gate: `genotype_probe compare --alleles-a/--alleles-b` isolates exactly one block.
#
#   tests/genotype_pair_override.sh <build dir> <scratch dir>
#
# Comparing two runs' winners does NOT isolate a block. The HG04036 sv:7 comparison that way differed
# in seven of eleven blocks and reported a 147.58-nat advantage for the longer allele; the same
# comparison with only sv:7 replaced reports 19.44. Attributing a score difference to one substitution
# therefore requires the override, and this gate holds it to that: one differing block, a refusal when
# an allele is out of range, and an allele scan whose floor reproduces the recorded availability floor.
#
# The C4 campaign artifacts live under tests/results, which is gitignored, so a clean checkout has
# nothing to run against. That is not a failure: the gate checks what it can and says what it skipped.
set -uo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="${1:?usage: genotype_pair_override.sh <build dir> <scratch dir>}"
WORK="${2:?usage: genotype_pair_override.sh <build dir> <scratch dir>}"
PROBE="$BUILD/genotype_probe"
PKG="$REPO/tests/results/genotype_loo_c4_development_7e37bb4"
REF='HG00096#1#haplotype1-0000024:31848049-32074186'
rm -rf "$WORK"; mkdir -p "$WORK"

fail=0
ok()   { echo "  ok: $1"; }
bad()  { echo "  GATE FAILED: $1"; fail=1; }
same() { if [ "$2" = "$3" ]; then ok "$1 = $2"; else bad "$1 is '$2', expected '$3'"; fi; }

[ -x "$PROBE" ] || { echo "genotype_probe not built at $PROBE" >&2; exit 1; }

if [ ! -d "$PKG/inputs" ]; then
  echo "skipped: $PKG/inputs is absent (tests/results is gitignored)."
  echo "         The override's argument handling is exercised only where the campaign artifacts exist."
  exit 0
fi

# Inputs: the preserved LOO graph, full graph and 2x reads -- the smallest read set, because this gate
# checks the override's semantics rather than any scientific quantity.
D=NA20509
G="$WORK/c4.gfa"; L="$WORK/c4.loo.gfa"; R1="$WORK/r1.fq"; R2="$WORK/r2.fq"
gzip -dc "$PKG/inputs/c4.gfa.gz"            > "$G"  2>/dev/null || { echo "skipped: no full graph"; exit 0; }
gzip -dc "$PKG/inputs/$D/c4.loo.gfa.gz"     > "$L"  2>/dev/null || { echo "skipped: no LOO graph"; exit 0; }
gzip -dc "$PKG/inputs/$D/cov2/reads_1.fq.gz" > "$R1" 2>/dev/null || { echo "skipped: no reads"; exit 0; }
gzip -dc "$PKG/inputs/$D/cov2/reads_2.fq.gz" > "$R2" 2>/dev/null || { echo "skipped: no reads"; exit 0; }
LIB="$PKG/inputs/$D/cov2/library.txt"
RUN="$PKG/$D/call.2.K0"
[ -f "$RUN.mosaic_scores.tsv" ] || { echo "skipped: no preserved run tables"; exit 0; }

run_compare() { # <out> <extra args...>
  local out="$1"; shift
  "$PROBE" compare --reference "$REF" --donor "$D" --block sv:7 \
    --model-gfa "$L" --truth-gfa "$G" --run-a "$RUN" \
    --reads1 "$R1" --reads2 "$R2" --library "$LIB" --depth 2 \
    --out "$out" "$@" > "$out.log" 2>&1
}

echo "== an out-of-range allele is refused =="
if run_compare "$WORK/bad" --alleles-a 999,0 --alleles-b 0,0; then
  bad "an allele outside the block was accepted"
else
  if grep -q "names an allele outside block sv:7" "$WORK/bad.log"; then
    ok "refused by name, and said which block"
  else
    bad "refused, but not with the expected message"; head -2 "$WORK/bad.log" | sed 's/^/      /'
  fi
fi

echo "== a malformed specification is refused =="
if run_compare "$WORK/malformed" --alleles-a 12; then
  bad "a specification without a comma was accepted"
else
  grep -q "expects <allele1>,<allele2>" "$WORK/malformed.log" && ok "refused a missing second allele" \
    || bad "refused, but not with the expected message"
fi

echo "== overriding one block leaves exactly one differing block =="
if run_compare "$WORK/single" --alleles-a 12,11 --alleles-b 111,11; then
  same "differing_blocks" "$(awk -F'\t' '$1=="differing_blocks"{print $2}' "$WORK/single.summary.tsv")" "sv:7"
  same "alleles_a_override" "$(awk -F'\t' '$1=="alleles_a_override"{print $2}' "$WORK/single.summary.tsv")" "12,11"
  same "run_b (optional, omitted here)" "$(awk -F'\t' '$1=="run_b"{print $2}' "$WORK/single.summary.tsv")" "none"
  # The complete likelihood must be reported, not the fragment term alone: a fragment-only comparison
  # inverted this model's preference once already.
  for key in log_likelihood_a log_likelihood_b total_gap_b_minus_a exposure_gap_b_minus_a; do
    v="$(awk -F'\t' -v k="$key" '$1==k{print $2}' "$WORK/single.summary.tsv")"
    [ -n "$v" ] && ok "$key reported ($v)" || bad "$key is missing from the summary"
  done
  # total = fragment + exposure, to the tolerance the decomposition already declares.
  awk -F'\t' '$1=="fragment_term_gap_b_minus_a"{f=$2} $1=="exposure_gap_b_minus_a"{e=$2}
              $1=="total_gap_b_minus_a"{t=$2}
    END{ d=t-(f+e); if (d<0) d=-d;
         if (d < 1e-6) printf "  ok: total = fragment + exposure (residual %.3g)\n", d;
         else { printf "  GATE FAILED: total - (fragment + exposure) = %.6g\n", d; exit 1 } }' \
    "$WORK/single.summary.tsv" || fail=1
else
  bad "the single-block override run failed"; head -3 "$WORK/single.log" | sed 's/^/      /'
fi

echo "== the allele scan reproduces the recorded availability floor =="
if [ -f "$WORK/single.alleles.tsv" ]; then
  f1="$(awk -F'\t' '$1=="floor_distance_truth1"{print $2}' "$WORK/single.summary.tsv")"
  f2="$(awk -F'\t' '$1=="floor_distance_truth2"{print $2}' "$WORK/single.summary.tsv")"
  same "floor total at sv:7 (recorded availability floor is 25)" "$((f1 + f2))" "25"
  n="$(awk -F'\t' 'NR>1' "$WORK/single.alleles.tsv" | wc -l | tr -d ' ')"
  same "alleles scanned" "$n" "$(awk -F'\t' '$1=="block_alleles"{print $2}' "$WORK/single.summary.tsv")"
else
  bad "no allele scan was written"
fi

echo
if [ "$fail" -ne 0 ]; then echo "pair-override gate FAILED"; exit 1; fi
echo "pair-override gate passed"
