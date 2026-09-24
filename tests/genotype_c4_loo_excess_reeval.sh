#!/usr/bin/env bash
# Re-evaluates every C4 production-v1 LOO CALL with the evaluator that reports floor alleles, into a separate
# tree (the frozen package is never touched). Evaluation only: truth is read, no call is made or changed.
#
#   tests/genotype_c4_loo_excess_reeval.sh <evaluator (genotype_probe)> <campaign dir> <out dir>
set -uo pipefail
PROBE="${1:?}"; CAMPAIGN="${2:?}"; OUT="${3:?}"
REF='HG00096#1#haplotype1-0000024:31848049-32074186'
mkdir -p "$OUT"
shasum -a 256 "$PROBE" | cut -d' ' -f1 > "$OUT/evaluator.sha256"
awk -F'\t' 'NR>1 && $2=="loo" && ($5=="strict_call" || $5=="genotype_equivalent_call"){print $1}' "$CAMPAIGN/completeness.tsv" |
while read -r job; do
  row="$(awk -F'\t' -v j="$job" '$1==j' "$CAMPAIGN/jobs.tsv")"
  IFS=$'\t' read -r _ _ _ _ scheme label model _ _ _ _ truth _ _ _ _ _ prefix <<< "$row"
  [ -f "$OUT/$job.summary.tsv" ] && continue
  "$PROBE" evaluate --scheme "$scheme" --donor "$label" --reference "$REF" --run-prefix "$prefix" --model-gfa "$model" \
    --truth-gfa "$truth" --out "$OUT/$job" > "$OUT/$job.log" 2>&1 || echo "evaluation failed: $job"
  rm -f "$OUT/$job".*.sorted.gfa
done
echo "re-evaluated $(ls "$OUT"/*.summary.tsv | wc -l | tr -d ' ') LOO calls"
