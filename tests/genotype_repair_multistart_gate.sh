#!/usr/bin/env bash
# Gates for the diagnostic multi-start repair search, on the frozen NA20509 10x LOO case.
#
#   tests/genotype_repair_multistart_gate.sh <build dir> <scratch dir>
#
#   1  seed construction is invariant to the order of the input candidate table
#   2  exactly at max_repair_starts proceeds; one below refuses before repair, with no final output
#   3  seed-order reversal: identical per-seed trajectories, pooled candidates and final ranking
#   4  shared-cache and isolated-context executions: identical scientific results
#   5  1 thread and 8 threads: identical
#   6  truth present and truth absent: identical search tables
#   7  a one-seed multi-start reproduces the single-start search bit for bit
#   8  per-seed and cumulative cap-minus-one refusals leave no final output
#
# "Scientific" means the search itself: trajectories (rounds columns 1-28 and the delta's work), seeds,
# optima, pooled identities and the exact rankings. Per-seed work totals, emission-cache hit counts and
# timings legitimately depend on cache warmth and processing order and are excluded where noted.
set -uo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="${1:?usage: genotype_repair_multistart_gate.sh <build dir> <scratch dir>}"
WORK="${2:?usage: genotype_repair_multistart_gate.sh <build dir> <scratch dir>}"
PROBE="$BUILD/genotype_probe"
PKG="$REPO/tests/results/genotype_loo_c4_development_7e37bb4"
REF='HG00096#1#haplotype1-0000024:31848049-32074186'
D="$PKG/inputs/NA20509"
[ -x "$PROBE" ] || { echo "genotype_probe not built at $PROBE" >&2; exit 1; }
if [ ! -f "$D/c4.loo.gfa.gz" ] || [ ! -f "$PKG/NA20509/call.10.K128.mosaic_candidates.tsv" ]; then
  echo "skipped: the frozen campaign artifacts under tests/results are absent (gitignored)."; exit 0
fi
rm -rf "$WORK"; mkdir -p "$WORK"
G="$WORK/c4.gfa"; L="$WORK/c4.loo.gfa"; R1="$WORK/r1.fq"; R2="$WORK/r2.fq"
gzip -dc "$PKG/inputs/c4.gfa.gz" > "$G"; gzip -dc "$D/c4.loo.gfa.gz" > "$L"
gzip -dc "$D/cov10/reads_1.fq.gz" > "$R1"; gzip -dc "$D/cov10/reads_2.fq.gz" > "$R2"

fail=0
same() { if [ "$2" = "$3" ]; then echo "  ok: $1 = $2"; else echo "  GATE FAILED: $1 is '$2', expected '$3'"; fail=1; fi; }
exists() { [ -e "$1" ] && echo present || echo absent; }
key() { awk -F'\t' -v k="$2" '$1==k{print $2}' "$1"; }
col() { awk -F'\t' -v r="$2" -v c="$3" 'NR==1{for(i=1;i<=NF;i++)h[$i]=i} $1==r{print $h[c]}' "$1"; }
pick() {   # pick <table> <columns...>; a missing table can never compare equal to anything
  local t="$1"; shift
  if [ ! -f "$t" ]; then echo "MISSING TABLE $t $RANDOM$RANDOM"; return; fi
  awk -F'\t' -v want="$*" 'BEGIN{n=split(want,w," ")} NR==1{for(i=1;i<=NF;i++)h[$i]=i}
    {line=""; for(k=1;k<=n;k++){line=line (k>1?"\t":"") $h[w[k]]} print line}' "$t"
}
search() {   # search <out> <extra args...>
  local o="$1"; shift
  mkdir -p "$(dirname "$o")"
  "$PROBE" repair-search --reference "$REF" --donor NA20509 --model-gfa "$L" --run "$PKG/NA20509/call.10.K128" \
    --reads1 "$R1" --reads2 "$R2" --library "$D/cov10/library.txt" \
    --delta 1 --origin-cache 1 --origin-cache-bytes 402653184 --origin-cache-policy planned \
    --max-physical-origins 19600000000 --max-state-classes 16800000000 --max-piece-hits 56000000000 \
    --max-latent-pairs 56000000000 --local-pairs 32 --rounds 11 --max-neighbors 670 "$@" --out "$o" > "$o.log" 2>&1
  echo "$?"
}
MS=(--multi-start 1 --max-repair-starts 16 --seed-max-physical-origins 4900000000 --seed-max-state-classes 4200000000
    --seed-max-piece-hits 14000000000 --seed-max-latent-pairs 14000000000)
SEARCH_COLS="round start_total moved repaired_block local_pair_rank assignment allele_x allele_y hap1_alleles hap2_alleles fragment_term exposure_term prior total delta_total status unique_best separation neighbours_generated neighbours_scored duplicate suppressed_visited self_skipped rejected_infinite_prior scored_pairs order_invariant swap_invariant origins_reused origins_recomputed fragments_reduced fragments_copied"
SEED_COLS="seed candidate_id selected_as group_rank group_sources union_source start_score unordered_genotype hap1_alleles hap2_alleles rounds moves stop_reason optimum_total optimum_hap1 optimum_hap2 optimum_genotype"
RANK_COLS="rank score fragment_term exposure_term prior numeric_error_bound in_winner_set members union_members rep_hap1 rep_hap2 unordered_genotype"
SUMMARY_KEYS="seeds pool_identities pool_score_mismatches initial_class initial_winners final_class final_winners final_best_score final_certified_gap call_status distinct_optima trajectory_merges"
same_search() {   # same_search <label> <reference prefix> <other prefix>
  local label="$1" a="$2" b="$3" ok=1 n
  if [ ! -f "$a.summary.tsv" ] || [ ! -f "$b.summary.tsv" ]; then
    echo "  GATE FAILED: $label: a run to compare produced no summary"; fail=1; return
  fi
  n="$(key "$a.summary.tsv" seeds)"
  for ((k = 0; k < n; k++)); do
    diff -q <(pick "$a.seed$k.rounds.tsv" $SEARCH_COLS) <(pick "$b.seed$k.rounds.tsv" $SEARCH_COLS) > /dev/null || { echo "  GATE FAILED: $label: seed $k trajectory differs"; ok=0; }
  done
  diff -q <(pick "$a.seeds.tsv" $SEED_COLS) <(pick "$b.seeds.tsv" $SEED_COLS) > /dev/null || { echo "  GATE FAILED: $label: seed outcomes differ"; ok=0; }
  for t in initial_ranking final_ranking; do
    diff -q <(pick "$a.$t.tsv" $RANK_COLS) <(pick "$b.$t.tsv" $RANK_COLS) > /dev/null || { echo "  GATE FAILED: $label: $t differs"; ok=0; }
  done
  for k2 in $SUMMARY_KEYS; do
    [ "$(key "$a.summary.tsv" "$k2")" = "$(key "$b.summary.tsv" "$k2")" ] || { echo "  GATE FAILED: $label: summary $k2 differs"; ok=0; }
  done
  if [ "$ok" -eq 1 ]; then echo "  ok: $label: identical trajectories, seeds, pooled identities and rankings"; else fail=1; fi
}

echo "== reference multi-start (8 threads, shared context, canonical seed order) =="
same "reference exit" "$(search "$WORK/base/ms" --threads 8 "${MS[@]}" --truth-gfa "$G")" "0"
same "reference status" "$(key "$WORK/base/ms.run_status.tsv" status)" "COMPLETE"
same "reference pool mismatches" "$(key "$WORK/base/ms.summary.tsv" pool_score_mismatches)" "0"
SEEDS="$(key "$WORK/base/ms.summary.tsv" seeds)"
echo "  seeds: $SEEDS; call status $(key "$WORK/base/ms.summary.tsv" call_status)"

echo "== 1. seed construction is invariant to input candidate order =="
for s in 7 1234567; do
  same "shuffled $s exit" "$(search "$WORK/shuffle$s/ms" --threads 8 "${MS[@]}" --seeds-only 1 --shuffle-union-input $s)" "0"
  if diff -q <(pick "$WORK/base/ms.seed_set.tsv" seed candidate_id selected_as group_rank group_sources score hap1_alleles hap2_alleles) \
             <(pick "$WORK/shuffle$s/ms.seed_set.tsv" seed candidate_id selected_as group_rank group_sources score hap1_alleles hap2_alleles) > /dev/null; then
    echo "  ok: input order $s gives the identical seed set"
  else echo "  GATE FAILED: input order $s changed the seed set"; fail=1; fi
done

echo "== 2. max_repair_starts: exactly at the bound proceeds, one below refuses without truncation =="
same "at the bound" "$(search "$WORK/atcap/ms" --threads 8 "${MS[@]}" --seeds-only 1 --max-repair-starts "$SEEDS")" "0"
same "one below" "$(search "$WORK/under/ms" --threads 8 "${MS[@]}" --seeds-only 1 --max-repair-starts "$((SEEDS - 1))")" "5"
same "one below status" "$(key "$WORK/under/ms.run_status.tsv" status)" "RESOURCE_REFUSAL"
for f in summary seeds final_ranking seed_set; do same "one below wrote no $f" "$(exists "$WORK/under/ms.$f.tsv")" "absent"; done

echo "== 3. seed-order reversal =="
same "reversed exit" "$(search "$WORK/reversed/ms" --threads 8 "${MS[@]}" --reverse-seeds 1 --truth-gfa "$G")" "0"
same_search "reversed seed order" "$WORK/base/ms" "$WORK/reversed/ms"

echo "== 4. shared cache versus isolated contexts =="
same "isolated exit" "$(search "$WORK/isolated/ms" --threads 8 "${MS[@]}" --isolate-seeds 1 --truth-gfa "$G")" "0"
same_search "isolated contexts" "$WORK/base/ms" "$WORK/isolated/ms"

echo "== 5. 1 thread versus 8 threads =="
same "1-thread exit" "$(search "$WORK/one_thread/ms" --threads 1 "${MS[@]}" --truth-gfa "$G")" "0"
same_search "1 thread" "$WORK/base/ms" "$WORK/one_thread/ms"

echo "== 6. truth present versus truth absent =="
same "truth-absent exit" "$(search "$WORK/no_truth/ms" --threads 8 "${MS[@]}")" "0"
same_search "truth absent" "$WORK/base/ms" "$WORK/no_truth/ms"

echo "== 7. one-seed multi-start reproduces the single-start search =="
same "single-start exit" "$(search "$WORK/single/rs" --threads 8 --truth-gfa "$G")" "0"
same "one-seed exit" "$(search "$WORK/oneseed/ms" --threads 8 "${MS[@]}" --seed-candidates 329 --truth-gfa "$G")" "0"
same "one-seed seeds" "$(key "$WORK/oneseed/ms.summary.tsv" seeds)" "1"
WORK_COLS="$SEARCH_COLS unique_haplotypes haplotypes_new haplotypes_reused bases_scanned piece_hits physical_origins state_classes window_bases oc_hits oc_misses oc_stored origins_constructed origins_loaded"
if diff -q <(pick "$WORK/single/rs.rounds.tsv" $WORK_COLS) <(pick "$WORK/oneseed/ms.seed0.rounds.tsv" $WORK_COLS) > /dev/null; then
  echo "  ok: the one-seed trajectory is identical to single-start in every search and work column"
else echo "  GATE FAILED: the one-seed trajectory differs from single-start"; fail=1; fi

echo "== 8. per-seed and cumulative cap-minus-one refusals leave no final output =="
PER_SEED="$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next}{if($h["origins_used"]>m)m=$h["origins_used"]}END{print m}' "$WORK/base/ms.seeds.tsv")"
TOTAL="$(key "$WORK/base/ms.summary.tsv" cumulative_physical_origins)"
same "per-seed bound minus one exit" "$(search "$WORK/seedcap/ms" --threads 8 "${MS[@]}" --seed-max-physical-origins "$((PER_SEED - 1))" --truth-gfa "$G")" "5"
same "cumulative bound minus one exit" "$(search "$WORK/processcap/ms" --threads 8 "${MS[@]}" --max-physical-origins "$((TOTAL - 1))" --truth-gfa "$G")" "5"
for o in seedcap processcap; do
  same "$o status" "$(key "$WORK/$o/ms.run_status.tsv" status)" "RESOURCE_REFUSAL"
  for f in summary seeds final_ranking initial_ranking; do same "$o wrote no $f" "$(exists "$WORK/$o/ms.$f.tsv")" "absent"; done
  same "$o staged graph removed" "$(exists "$WORK/$o/ms.model.sorted.gfa")" "absent"
done

echo
if [ "$fail" -ne 0 ]; then echo "multi-start gate FAILED"; exit 1; fi
echo "multi-start gate passed"
