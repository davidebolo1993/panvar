#!/usr/bin/env bash
# Thread identity gate (docs/algorithms/genotype_dosage.md section 14.6): the shadow audit of one persisted campaign call,
# with repair, provenance and the allocation test, at 1, 2, 4 and 8 threads. Every dosage table must be byte-identical;
# only dosage_status.tsv and dosage_repair_scoring.tsv (wall times, scheduling-dependent cache hits) are exempt. Usage: <genotype_probe> <jobs.tsv> <cell> <scratch dir>
set -uo pipefail
PROBE=$1; JOBS=$2; CELL=$3; OUT=$4
mkdir -p "$OUT"
IFS=$'\t' read -r job camp donor cov scheme lab model r1 r2 lib cur truth x1 x2 x3 x4 salt prefix <<< "$(awk -F'\t' -v j="$CELL" '$1==j' "$JOBS")"
for t in 1 2 4 8; do
  "$PROBE" dosage-replay --model-gfa "$model" --reference 'HG00096#1#haplotype1-0000024:31848049-32074186' --panel-curation "$cur" \
    --run "$prefix" --reads1 "$r1" --reads2 "$r2" --library "$lib" --threads $t --repair 1 --allocation-test 1 --out "$OUT/t$t" > "$OUT/t$t.log" 2>&1 \
    || { echo "FAIL run at $t threads"; exit 1; }
done
fail=0
for f in "$OUT"/t1.dosage_*; do
  suf=${f#"$OUT/t1"}
  case "$suf" in .dosage_status.tsv|.dosage_repair_scoring.tsv) continue;; esac
  # output-prefix paths inside FASTA headers (the equivalence-set link) are replaced, as the shadow identity gate does
  for t in 2 4 8; do cmp -s <(sed "s#$OUT/t1#PREFIX#g" "$f") <(sed "s#$OUT/t$t#PREFIX#g" "$OUT/t$t$suf") || { echo "FAIL $suf differs at $t threads"; fail=1; }; done
done
n=$(ls "$OUT"/t1.dosage_* | grep -vc status)
[ $fail -eq 0 ] && echo "PASS $CELL: $n dosage tables byte-identical at 1/2/4/8 threads" || echo "FAILED $CELL"
exit $fail
