#!/usr/bin/env bash
# The representation-conflict guard (representation_conflict_rule_v1) on three frozen C4 cells, through the
# production path (v2 profile).
#
#   tests/genotype_representation_guard_gate.sh <build dir> <campaign dir> <diagnostics_v2 dir> <out dir>
#
#   1  annotation only: guard on and off give identical pre-filter outputs -- every table apart from the
#      appended filter columns, the status apart from the guard's own keys, the FASTA apart from its
#      locus_filter field; and a v1-call cell's pre-filter tables equal the frozen campaign job's
#   2  the flags are what the frozen truth-free diagnostic table says: HG01114.2 and HG01114.5 sv:7
#      LOWQUAL_REPRESENTATION_CONFLICT with the same evidence, every other called block PASS, and the LZO
#      control PASS at every block
#   3  1, 4 and 8 threads give identical filter columns
#   4  no truth input: the only inputs the run hashed (and so read) are the graph, the reads and the library
set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="${1:?usage: genotype_representation_guard_gate.sh <build> <campaign> <diagnostics_v2> <out>}"
CAMPAIGN="${2:?}"; DIAG="${3:?}"; OUT="${4:?}"
PROBE="$BUILD/genotype_probe"
PROFILE="$REPO/profiles/genotype_c4_v2_campaign_local.tsv"
REF='HG00096#1#haplotype1-0000024:31848049-32074186'
FROZEN="$(dirname "$DIAG")/v2_representation_guard/truth_free_flags.tsv"
mkdir -p "$OUT"
fail=0
same() { if [ "$2" = "$3" ]; then echo "  ok: $1"; else echo "  GATE FAILED: $1 ($2 vs $3)"; fail=1; fi; }
key() { awk -F'\t' -v k="$1" '$1==k{print $2; exit}' "$2"; }
digest() { shasum -a 256 | cut -d' ' -f1; }
run() {   # run <name> <job> <threads> [flags]
  local name="$1" job="$2" threads="$3"; shift 3
  IFS=$'\t' read -r _ _ _ _ _ _ model r1 r2 lib cur _ <<< "$(awk -F'\t' -v j="$job" '$1==j' "$CAMPAIGN/jobs.tsv")"
  local CUR=(); [ "$cur" != none ] && CUR=(--panel-curation "$cur")
  mkdir -p "$OUT/$name"
  "$PROBE" tie-branch --diagnostic 0 --profile "$PROFILE" --threads "$threads" --gfa "$model" --reference "$REF" \
    --reads1 "$r1" --reads2 "$r2" --library "$lib" ${CUR[@]+"${CUR[@]}"} --out-prefix "$OUT/$name/call" "$@" \
    > "$OUT/$name/log" 2>&1
  rm -f "$OUT/$name"/*.sorted.gfa
}
prefilter() {   # every output apart from the guard's own columns, keys and FASTA field
  local p="$1"
  for t in genotypes block_genotypes; do awk -F'\t' 'BEGIN{OFS="\t"}{NF=NF-8; print}' "$p.$t.tsv"; done
  cat "$p.mosaic_candidates.tsv" "$p.mosaic_scores.tsv" "$p.equivalence_set.tsv"
  cut -f1-28 "$p.repair_rounds.tsv"
  [ -f "$p.mosaic_haplotypes.fa" ] && sed -e 's/ locus_filter=[A-Z_]*//' -e 's/^>[^ ]*_hap/>hap/' "$p.mosaic_haplotypes.fa"
  [ -f "$p.tie_branches.tsv" ] && cut -f1-13 "$p.tie_branches.tsv"
  grep -vE '^(seconds_|resources_|threads_|phys_footprint|rss_|emission_cache_dp_hits|emission_cache_dp_misses|emission_cache_executed_dps|scorer_|representation_|locus_filter|lowqual_blocks|settings_overridden|sha256_)' "$p.mosaic_status.tsv"
}
filters() { awk -F'\t' 'NR>1{print $1"\t"$(NF-7)"\t"$(NF-6)"\t"$(NF-5)"\t"$(NF-4)"\t"$(NF-3)"\t"$(NF-2)"\t"$(NF-1)"\t"$NF}' "$1.block_genotypes.tsv"; }
evidence() {   # <run prefix> <block>: filter, marker margin, fragment diff, exposure diff, to one decimal
  awk -F'\t' -v b="$2" 'NR>1 && $1==b{printf "%s %.1f %.1f %.1f", $(NF-7), $(NF-5), $(NF-3), $(NF-2)}' "$1.block_genotypes.tsv"
}
frozen_evidence() {   # <cell> <block> from the frozen truth-free table
  awk -F'\t' -v c="$1" -v b="$2" 'NR==1{for(i=1;i<=NF;i++) col[$i]=i; next} $col["cell"]==c && $col["block"]==b{printf "%s %.1f %.1f %.1f", $col["filter"], $col["marker_margin"], $col["fragment_diff"], $col["exposure_diff"]}' "$FROZEN"
}

echo "== runs"
run h2_on8 loo.HG01114.2 8
run h2_off8 loo.HG01114.2 8 --representation-guard 0
run h2_on4 loo.HG01114.2 4
run h2_on1 loo.HG01114.2 1
run lzo_on8 lzo.HG03732.5 8
run lzo_off8 lzo.HG03732.5 8 --representation-guard 0
run h5_on8 loo.HG01114.5 8
run h5_off8 loo.HG01114.5 8 --representation-guard 0

echo "== 1. annotation only"
for pair in h2 lzo h5; do
  same "$pair: pre-filter outputs, guard on = off" "$(prefilter "$OUT/${pair}_on8/call" | digest)" "$(prefilter "$OUT/${pair}_off8/call" | digest)"
done
for pair in h2:loo.HG01114.2 lzo:lzo.HG03732.5; do
  name="${pair%%:*}"; job="${pair#*:}"
  frozen="$(awk -F'\t' -v j="$job" '$1==j{print $18}' "$CAMPAIGN/jobs.tsv")"
  same "$name: pre-filter tables = the frozen campaign job's" \
    "$( { for t in genotypes block_genotypes; do awk -F'\t' 'BEGIN{OFS="\t"}{NF=NF-8; print}' "$OUT/${name}_on8/call.$t.tsv"; done; cat "$OUT/${name}_on8/call."{mosaic_candidates,mosaic_scores,equivalence_set}.tsv; } | digest)" \
    "$( { cat "$frozen."{genotypes,block_genotypes,mosaic_candidates,mosaic_scores,equivalence_set}.tsv; } | digest)"
done
same "guard off writes NOT_EVALUATED" "$(filters "$OUT/h2_off8/call" | cut -f2 | sort -u | tr '\n' ' ')" "NOT_EVALUATED "

echo "== 2. the flags"
same "HG01114.2 lowqual blocks" "$(key lowqual_blocks "$OUT/h2_on8/call.mosaic_status.tsv")" "sv:7"
same "HG01114.2 locus filter" "$(key locus_filter "$OUT/h2_on8/call.mosaic_status.tsv")" "LOWQUAL_REPRESENTATION_CONFLICT"
same "HG01114.2 sv:7 evidence = the frozen diagnostic's" "$(evidence "$OUT/h2_on8/call" sv:7)" "$(frozen_evidence loo.HG01114.2 sv:7)"
same "HG01114.5 (tie-branched) lowqual blocks" "$(key lowqual_blocks "$OUT/h5_on8/call.mosaic_status.tsv")" "sv:7"
same "HG01114.5 sv:7 evidence = the frozen diagnostic's" "$(evidence "$OUT/h5_on8/call" sv:7)" "$(frozen_evidence loo.HG01114.5 sv:7)"
same "HG01114.5 is still a tie-branched call" "$(key search_conditional "$OUT/h5_on8/call.mosaic_status.tsv")/$(key call_class "$OUT/h5_on8/call.mosaic_status.tsv")" "tie_branched/genotype_equivalent_call"
same "LZO control: every block PASS" "$(filters "$OUT/lzo_on8/call" | cut -f2 | sort -u | tr '\n' ' ')" "PASS "
same "LZO control locus filter" "$(key locus_filter "$OUT/lzo_on8/call.mosaic_status.tsv")" "PASS"
same "the FASTA records the locus filter" "$(grep -c 'locus_filter=PASS' "$OUT/lzo_on8/call.mosaic_haplotypes.fa")" "2"

echo "== 3. thread invariance"
same "filter columns: 8 = 4 threads" "$(filters "$OUT/h2_on8/call" | digest)" "$(filters "$OUT/h2_on4/call" | digest)"
same "filter columns: 8 = 1 thread" "$(filters "$OUT/h2_on8/call" | digest)" "$(filters "$OUT/h2_on1/call" | digest)"

echo "== 4. no truth input"
for name in h2_on8 lzo_on8 h5_on8; do
  same "$name: inputs hashed" "$(grep '^sha256_' "$OUT/$name/call.mosaic_status.tsv" | cut -f1 | sort | tr '\n' ' ')" "sha256_gfa sha256_library sha256_reads1 sha256_reads2 "
done
[ "$fail" -eq 0 ] && echo "representation guard gate: PASS" || { echo "representation guard gate: FAILED"; exit 1; }
