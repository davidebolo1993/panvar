#!/usr/bin/env bash
# Known-truth COMPLETE_DIPLOTYPE_REPAIR gate (docs/algorithms/genotype_dosage.md section 22.6). The four C4 pair-route
# repairs have no defined truth at sv:1, so this gate makes one: a leave-zero-out call that is exact at <block> gets a
# dosage error injected there (dosage-synthetic, scored exactly as the repair guard scores it), and repair mode must
# find the true pair again. Checks: the stage classes the block COMPLETE_DIPLOTYPE_REPAIR and applies it; the repaired
# pair and both repaired haplotypes equal the real call's; the evaluator scores the injected original at a positive
# distance and the repaired prefix exactly as it scores the real call.
# Usage: genotype_dosage_synthetic_gate.sh <genotype_probe> <jobs.tsv> <lzo cell> <block> <hap1 allele|hap2 allele> <scratch dir>
set -uo pipefail
PROBE=$1; JOBS=$2; CELL=$3; BLOCK=$4; INJECT=$5; OUT=$6
REF='HG00096#1#haplotype1-0000024:31848049-32074186'
mkdir -p "$OUT"
IFS=$'\t' read -r job camp donor cov scheme lab model r1 r2 lib cur truth x1 x2 x3 x4 salt prefix <<< "$(awk -F'\t' -v j="$CELL" '$1==j' "$JOBS")"
[ -n "${prefix:-}" ] || { echo "no job $CELL"; exit 2; }
fail=0
"$PROBE" dosage-synthetic --model-gfa "$model" --reference "$REF" --panel-curation "$cur" --run "$prefix" --reads1 "$r1" --reads2 "$r2" \
  --library "$lib" --block "$BLOCK" --alleles "$INJECT" --threads 8 --out "$OUT/injected" > "$OUT/synthetic.log" 2>&1 || { echo "FAIL dosage-synthetic"; cat "$OUT/synthetic.log"; exit 1; }
"$PROBE" dosage-replay --model-gfa "$model" --reference "$REF" --panel-curation "$cur" --run "$OUT/injected" --reads1 "$r1" --reads2 "$r2" \
  --library "$lib" --threads 8 --mode repair --out "$OUT/repaired" > "$OUT/replay.log" 2>&1 || { echo "FAIL dosage-replay"; tail -3 "$OUT/replay.log"; exit 1; }
cls=$(awk -F'\t' -v b="$BLOCK" 'NR==1{for(i=1;i<=NF;i++)h[$i]=i} $1==b{print $h["repair_class"]}' "$OUT/repaired.dosage_blocks.tsv")
[ "$cls" = COMPLETE_DIPLOTYPE_REPAIR ] || { echo "FAIL $BLOCK classed $cls"; fail=1; }
applied=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i} NR==2{print $h["dosage_complete_applied_block_ids"]}' "$OUT/repaired.locus_genotype.tsv")
[ "$applied" = "$BLOCK" ] || { echo "FAIL applied '$applied'"; fail=1; }
real=$(awk -F'\t' -v b="$BLOCK" '$1==b{print $4"/"$5}' "$prefix.block_genotypes.tsv")
got=$(awk -F'\t' -v b="$BLOCK" '$1==b{print $4"/"$5}' "$OUT/repaired.block_genotypes.tsv")
inj=$(awk -F'\t' -v b="$BLOCK" '$1==b{print $4"/"$5}' "$OUT/repaired.original.block_genotypes.tsv")
echo "$BLOCK: real $real, injected $inj, repaired $got"
[ "$got" = "$real" ] || { echo "FAIL the repair did not restore the real pair"; fail=1; }
cmp -s <(grep -v '>' "$prefix.mosaic_haplotypes.fa") <(grep -v '>' "$OUT/repaired.mosaic_haplotypes.fa") || { echo "FAIL repaired haplotypes differ from the real call's"; fail=1; }
for v in repaired repaired.original; do
  "$PROBE" evaluate --scheme "$scheme" --donor "$donor" --reference "$REF" --run-prefix "$OUT/$v" --model-gfa "$model" --truth-gfa "$truth" \
    --out "$OUT/$v.eval" > "$OUT/$v.evaluate.log" 2>&1 || { echo "FAIL evaluate $v"; fail=1; }
done
d_orig=$(awk -F'\t' -v b="$BLOCK" '$4==b{print $10}' "$OUT/repaired.original.eval.blocks.tsv")
d_rep=$(awk -F'\t' -v b="$BLOCK" '$4==b{print $10}' "$OUT/repaired.eval.blocks.tsv")
echo "evaluator distance at $BLOCK: injected $d_orig, repaired $d_rep"
[ "${d_orig:-0}" -gt 0 ] 2>/dev/null || { echo "FAIL the injected call is not wrong at $BLOCK"; fail=1; }
cmp -s "$prefix.eval.blocks.tsv" "$OUT/repaired.eval.blocks.tsv" || { echo "FAIL the repaired prefix does not evaluate as the real call"; fail=1; }
grep -q "^call_source	locus_genotype_dosage_repaired" "$OUT/repaired.eval.summary.tsv" || { echo "FAIL the evaluator did not read the repaired call"; fail=1; }
[ $fail -eq 0 ] && echo "PASS $CELL $BLOCK $INJECT" || echo "FAILED $CELL"
exit $fail
