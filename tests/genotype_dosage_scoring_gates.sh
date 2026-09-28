#!/usr/bin/env bash
# Repair-scoring lifecycle gates (Sol, 2026-09-28, part C) on one persisted campaign call with an AUTO repair:
#   1. batch-size invariance: --repair-batch 1 (a release after every candidate, maximal slot reuse) and a batch larger
#      than the candidate set (no release until the end) give byte-identical candidate and repair tables;
#   2. exact-at-cap: the measured cumulative physical-origin total as the cap succeeds with identical tables; one less
#      refuses (TRIGGERED_UNREPAIRED / resource_refusal) with no candidate row;
#   3. the reproduction guard is first: a perturbed persisted total refuses (winner_not_reproduced) with no candidate
#      row and no repaired haplotype placed.
# Usage: <genotype_probe> <jobs.tsv> <cell> <scratch dir>
set -uo pipefail
PROBE=$1; JOBS=$2; CELL=$3; OUT=$4
mkdir -p "$OUT"
IFS=$'\t' read -r job camp donor cov scheme lab model r1 r2 lib cur truth x1 x2 x3 x4 salt prefix <<< "$(awk -F'\t' -v j="$CELL" '$1==j' "$JOBS")"
run() { local name=$1; shift
  "$PROBE" dosage-replay --model-gfa "$model" --reference 'HG00096#1#haplotype1-0000024:31848049-32074186' --panel-curation "$cur" \
    --run "$prefix" --reads1 "$r1" --reads2 "$r2" --library "$lib" --threads 8 --repair 1 --provenance 0 --out "$OUT/$name" "$@" > "$OUT/$name.log" 2>&1; }
fail=0
same() { cmp -s "$OUT/$1.dosage_repair_candidates.tsv" "$OUT/$2.dosage_repair_candidates.tsv" && cmp -s "$OUT/$1.dosage_repairs.tsv" "$OUT/$2.dosage_repairs.tsv"; }
run ref
run b1 --repair-batch 1
run bmax --repair-batch 1000000
same ref b1 && same ref bmax && echo "PASS batch invariance (batch 1 releases $(awk -F'\t' 'NR==2{print $11}' $OUT/b1.dosage_repair_scoring.tsv), slots $(awk -F'\t' 'NR==2{print $12"/"$13}' $OUT/b1.dosage_repair_scoring.tsv) used/created)" || { echo "FAIL batch invariance"; fail=1; }
P=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)c[$i]=i} NR==2{print $c["physical_origins"]}' "$OUT/ref.dosage_repair_scoring.tsv")
run cap --repair-max-physical-origins "$P"
run capm1 --repair-max-physical-origins $((P - 1))
same ref cap && echo "PASS exact-at-cap ($P physical origins) succeeds, identical" || { echo "FAIL exact-at-cap"; fail=1; }
st=$(awk -F'\t' 'NR>1 && $4 ~ /TRIGGERED_UNREPAIRED/ {print $4"/"$5}' "$OUT/capm1.dosage_blocks.tsv")
rows=$(($(wc -l < "$OUT/capm1.dosage_repair_candidates.tsv") - 1))
[ "$st" = "TRIGGERED_UNREPAIRED/resource_refusal" ] && [ "$rows" -eq 0 ] && echo "PASS cap-1 refuses ($st), $rows candidate rows" || { echo "FAIL cap-1: '$st', $rows rows"; fail=1; }
run guard --perturb-persisted-total 0.001
gst=$(awk -F'\t' 'NR>1 && $4 ~ /TRIGGERED_UNREPAIRED/ {print $4"/"$5}' "$OUT/guard.dosage_blocks.tsv")
grows=$(($(wc -l < "$OUT/guard.dosage_repair_candidates.tsv") - 1))
[ "$gst" = "TRIGGERED_UNREPAIRED/winner_not_reproduced" ] && [ "$grows" -eq 0 ] && echo "PASS guard first ($gst), $grows candidate rows" || { echo "FAIL guard: '$gst', $grows rows"; fail=1; }
[ $fail -eq 0 ] && echo "PASS $CELL" || echo "FAILED $CELL"
exit $fail
