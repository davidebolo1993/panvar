#!/usr/bin/env bash
# Gates for the DIAGNOSTIC (v2/post-hoc) improving-move-tie branching, on one frozen C4 production-v1 tie.
#
#   tests/genotype_tie_branch_gate.sh <build dir> <work dir> <campaign dir> [job]
#
#   1  production outputs under the diagnostic equal the frozen campaign job's (status apart from wall,
#      memory and worker gauges; the repair trajectory's search columns)
#   2  diagnostic on and off: identical production outputs; at 1 thread every counter, the emission
#      cache's included, since its hit/miss split depends on scheduling under 8 threads (the total of
#      requests does not) -- which is also why the 8-thread comparisons set those four counters aside
#   3  reversed branch order: identical decision, terminal ties, pooled top set and branch endpoints
#   4  4 threads against 8: identical diagnostic tables apart from wall and memory columns
#   5  truth-free diagnostic tables are byte-identical before and after truth is evaluated
set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="${1:?usage: genotype_tie_branch_gate.sh <build> <work> <campaign> [job]}"
WORK="${2:?}"
CAMPAIGN="${3:?}"
JOB="${4:-mosaic.both_homologues_2.2}"
PROBE="$BUILD/genotype_probe"
REF='HG00096#1#haplotype1-0000024:31848049-32074186'
PROFILE="$REPO/profiles/genotype_c4_campaign_local.tsv"
rm -rf "$WORK"; mkdir -p "$WORK"
row="$(awk -F'\t' -v j="$JOB" '$1==j' "$CAMPAIGN/jobs.tsv")"
[ -n "$row" ] || { echo "no job $JOB in $CAMPAIGN/jobs.tsv" >&2; exit 1; }
IFS=$'\t' read -r _ campaign subject cov scheme label model r1 r2 lib cur truth _ _ _ _ _ prefix <<< "$row"
CUR=(); [ "$cur" != none ] && CUR=(--panel-curation "$cur")
fail=0
same() { if [ "$2" = "$3" ]; then echo "  ok: $1"; else echo "  GATE FAILED: $1 ($2 vs $3)"; fail=1; fi; }
run() {   # run <name> <threads> [extra...]
  local name="$1" threads="$2"; shift 2
  mkdir -p "$WORK/$name"
  "$PROBE" tie-branch --profile "$PROFILE" --threads "$threads" --gfa "$model" --reference "$REF" --reads1 "$r1" \
    --reads2 "$r2" --library "$lib" "${CUR[@]}" --out-prefix "$WORK/$name/call" "$@" > "$WORK/$name/log" 2>&1
}
# The status apart from what a rerun legitimately changes: wall time, process memory and worker gauges.
status_science() {   # 8 threads: the cache's hit/miss split is scheduling-dependent
  grep -vE '^(seconds_|resources_|threads_|phys_footprint|rss_|emission_cache_dp_hits|emission_cache_dp_misses|emission_cache_executed_dps|scorer_emission_cache_arena_bytes)' "$1"
}
status_all() { grep -vE '^(seconds_|resources_|threads_|phys_footprint|rss_)' "$1"; }
rounds_science() { cut -f1-28 "$1"; }
diag_science() {   # the diagnostic's tables apart from wall and memory
  local p="$1"
  grep -vE '^(seconds|phys_footprint_max_bytes|branch_order)\b' "$p.diag.status.tsv"
  cat "$p.diag.terminal_ties.tsv" "$p.diag.pooled_top.tsv" 2>/dev/null
  awk -F'\t' 'BEGIN{OFS="\t"} {$15=""; $17=""; print}' "$p.diag.branches.tsv"
}
digest() { shasum -a 256 "$@" | cut -d' ' -f1 | tr '\n' ' '; }

echo "== $JOB =="
run on8 8
run off8 8 --diagnostic 0
run rev8 8 --reverse-branches 1
run on4 4
run on1 1
run off1 1 --diagnostic 0
same "diagnostic ran and finished" "$(awk -F'\t' '$1=="state"{print $2}' "$WORK/on8/call.diag.status.tsv")" "COMPLETE"
same "no diagnostic files when off" "$(ls "$WORK/off8" | grep -c '\.diag\.')" "0"
same "1. production status equals the frozen campaign's" \
  "$(diff <(status_science "$prefix.mosaic_status.tsv" | grep -v '^sha256_\|^profile\|^settings_overridden') <(status_science "$WORK/on8/call.mosaic_status.tsv" | grep -v '^sha256_\|^profile\|^settings_overridden') > /dev/null && echo identical || echo DIFFERS)" identical
same "1. production trajectory equals the frozen campaign's" \
  "$(diff <(rounds_science "$prefix.repair_rounds.tsv") <(rounds_science "$WORK/on8/call.repair_rounds.tsv") > /dev/null && echo identical || echo DIFFERS)" identical
same "2. production status: diagnostic on = off" \
  "$(diff <(status_science "$WORK/on8/call.mosaic_status.tsv") <(status_science "$WORK/off8/call.mosaic_status.tsv") > /dev/null && echo identical || echo DIFFERS)" identical
same "2. production trajectory: diagnostic on = off" \
  "$(diff <(rounds_science "$WORK/on8/call.repair_rounds.tsv") <(rounds_science "$WORK/off8/call.repair_rounds.tsv") > /dev/null && echo identical || echo DIFFERS)" identical
same "2. 1 thread: every production status counter, the cache's included, on = off" \
  "$(diff <(status_all "$WORK/on1/call.mosaic_status.tsv" | grep -v '^settings_overridden') <(status_all "$WORK/off1/call.mosaic_status.tsv" | grep -v '^settings_overridden') > /dev/null && echo identical || echo DIFFERS)" identical
same "2. no production genotype either way" \
  "$(ls "$WORK/on8" "$WORK/off8" | grep -cE '\.(genotypes|block_genotypes)\.tsv$|\.fa$')" "0"
same "3. reversed branch order: identical decision, ties, pooled top set" \
  "$(diff <(grep -vE '^(seconds|phys_footprint_max_bytes|branch_order)\b' "$WORK/on8/call.diag.status.tsv"; cat "$WORK/on8/call.diag.terminal_ties.tsv" "$WORK/on8/call.diag.pooled_top.tsv") \
          <(grep -vE '^(seconds|phys_footprint_max_bytes|branch_order)\b' "$WORK/rev8/call.diag.status.tsv"; cat "$WORK/rev8/call.diag.terminal_ties.tsv" "$WORK/rev8/call.diag.pooled_top.tsv") > /dev/null && echo identical || echo DIFFERS)" identical
same "3. reversed branch order: identical branch endpoints" \
  "$(diff <(cut -f3-14 "$WORK/on8/call.diag.branches.tsv" | sort) <(cut -f3-14 "$WORK/rev8/call.diag.branches.tsv" | sort) > /dev/null && echo identical || echo DIFFERS)" identical
same "4. 4 threads = 8 threads" \
  "$(diff <(diag_science "$WORK/on8/call") <(diag_science "$WORK/on4/call") > /dev/null && echo identical || echo DIFFERS)" identical
before="$(digest "$WORK/on8/call.diag."*.tsv)"
rep="$(awk -F'\t' '$1=="representative_hap1"{a=$2} $1=="representative_hap2"{b=$2} END{print a"|"b}' "$WORK/on8/call.diag.status.tsv")"
"$PROBE" evaluate --scheme "$scheme" --donor "$label" --reference "$REF" --run-prefix "$WORK/on8/call" --model-gfa "$model" \
  --truth-gfa "$truth" "${CUR[@]}" --fallback-alleles "$rep" --evaluated-label diagnostic_tie_branching_v2_posthoc_not_a_call \
  --out "$WORK/on8/call.diag_eval" > "$WORK/on8/eval.log" 2>&1
same "5. evaluation ran" "$?" "0"
same "5. truth-free tables unchanged by evaluation" "$(digest "$WORK/on8/call.diag."*.tsv)" "$before"
[ "$fail" -eq 0 ] && echo "tie-branch gate: PASS" || { echo "tie-branch gate: FAILED"; exit 1; }
