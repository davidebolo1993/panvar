#!/usr/bin/env bash
# Developer harness, not a regression test: does running donors concurrently cost more than it saves?
#
#   BUILD=<build dir> READS=<existing coverage tree> tests/genotype_jobs_benchmark.sh <out-dir>
#
# Round 8 saw four C4 runs take 42-61 minutes of wall clock on 57-146 seconds of CPU, and reruns in
# isolation finished in 21-59 seconds with identical calls. That is machine contention, not a slow
# caller, so the question is what concurrency actually buys. This measures total cohort makespan -- the
# number that matters to someone waiting for a cohort -- at JOBS=1 and JOBS=2.
#
# Reads are NOT simulated here. They are taken from an existing coverage tree, so both settings score
# byte-identical input and nothing but the scheduling differs. The calls of the two settings must match
# exactly; the script fails if they do not.
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="${BUILD:-$REPO/build}"
READS="${READS:?READS must point at a coverage tree holding cached reads}"
OUT="${1:?usage: BUILD=<build> READS=<coverage tree> tests/genotype_jobs_benchmark.sh <out-dir>}"
COVERAGE="${COVERAGE:-30}"
FROZEN="$REPO/tests/genotype_coverage_settings.tsv"
setting() { awk -F'\t' -v key="$1" '$1==key{print $2}' "$FROZEN"; }
STATE_LIMIT="$(setting state_limit)"
PANEL_K="$(setting panel_arm_mosaic_candidates)"
MOSAIC_K="$(setting mosaic_arm_mosaic_candidates)"
AUDIT="$(setting phase_audit_candidates)"
[ -n "$STATE_LIMIT" ] && [ -n "$PANEL_K" ] && [ -n "$MOSAIC_K" ] && [ -n "$AUDIT" ] ||
  { echo "incomplete $FROZEN" >&2; exit 1; }
REF="HG00096#1#haplotype1-0000024:31848049-32074186"
PROBE="$BUILD/genotype_probe"

# A fixed donor subset, taken in sorted order so the two settings score the same donors in the same
# order. Enough donors to keep two workers busy, few enough that the benchmark is minutes not hours.
DONOR_COUNT="${DONOR_COUNT:-8}"
DONORS="${DONORS:-$(ls -1 "$READS" | while read -r d; do
  [ -f "$READS/$d/cov$COVERAGE/reads_1.fq" ] && echo "$d"; done | sort | head -n "$DONOR_COUNT")}"
[ -n "$DONORS" ] || { echo "no donors with cached cov$COVERAGE reads under $READS" >&2; exit 1; }

mkdir -p "$OUT"
[ -z "$(ls -A "$OUT" 2>/dev/null)" ] || { echo "$OUT is not empty" >&2; exit 1; }
gzip -dc "$REPO/tests/real_data/c4.gfa.gz" > "$OUT/c4.gfa"

{
  printf 'key\tvalue\n'
  printf 'git_commit\t%s\n' "$(git -C "$REPO" rev-parse HEAD)"
  printf 'git_tree_dirty\t%s\n' "$(git -C "$REPO" status --porcelain | grep -q . && echo true || echo false)"
  printf 'probe_sha256\t%s\n' "$(shasum -a 256 "$PROBE" | cut -d' ' -f1)"
  printf 'script_sha256\t%s\n' "$(shasum -a 256 "${BASH_SOURCE[0]}" | cut -d' ' -f1)"
  printf 'graph_sha256\t%s\n' "$(shasum -a 256 "$OUT/c4.gfa" | cut -d' ' -f1)"
  printf 'reads_tree\t%s\ncoverage\t%s\ndonors\t%s\n' "$READS" "$COVERAGE" "$(echo $DONORS)"
  printf 'state_limit\t%s\npanel_arm_K\t%s\nmosaic_arm_K\t%s\nphase_audit\t%s\n' \
    "$STATE_LIMIT" "$PANEL_K" "$MOSAIC_K" "$AUDIT"
  printf 'nproc\t%s\n' "$(sysctl -n hw.ncpu 2>/dev/null || nproc)"
} > "$OUT/signature.tsv"

# Scores one donor from cached reads. No simulation and no subsampling: the input is already on disk and
# identical between settings, so only the scheduling differs.
run_donor() {
  local donor="$1" dir="$OUT/$JOBS_TAG/$1"
  mkdir -p "$dir"
  "$PROBE" genotype --gfa "$OUT/c4.gfa" --reference "$REF" \
    --reads1 "$READS/$donor/cov$COVERAGE/reads_1.fq" --reads2 "$READS/$donor/cov$COVERAGE/reads_2.fq" \
    --library "$READS/$donor/cov$COVERAGE/library.txt" \
    --out-prefix "$dir/lzo" --state-limit "$STATE_LIMIT" --prefixes "$PANEL_K,$MOSAIC_K" \
    --phase-audit "$AUDIT" > "$dir/lzo.log" 2>&1
  echo "done $donor"
}
export -f run_donor
export OUT PROBE REF READS COVERAGE STATE_LIMIT PANEL_K MOSAIC_K AUDIT

# The call hash ignores everything that is legitimately allowed to differ between two schedulings:
# timings, resource counters and the work counters that record how much of a shared pass each arm saw.
call_hash() {
  local tag="$1"
  find "$OUT/$tag" -name "lzo.K*.mosaic_status.tsv" | sort | while read -r f; do
    printf '%s\n' "${f#$OUT/$tag/}"
    grep -vE '^(seconds_|peak_rss|wall_|threads_|.*_evaluations|.*_work|beam_extensions)' "$f"
  done | shasum -a 256 | cut -d' ' -f1
}

for JOBS_TAG in jobs1 jobs2; do
  JOBS=${JOBS_TAG#jobs}
  export JOBS_TAG
  mkdir -p "$OUT/$JOBS_TAG"
  start=$(date +%s)
  printf '%s\n' $DONORS | xargs -P "$JOBS" -I{} bash -c 'run_donor "$@"' _ {} > "$OUT/$JOBS_TAG.log" 2>&1
  end=$(date +%s)
  printf '%s\t%s\n' "$JOBS_TAG" "$((end - start))" >> "$OUT/makespan.tsv"
  echo "$JOBS_TAG makespan $((end - start))s"
done

H1="$(call_hash jobs1)"
H2="$(call_hash jobs2)"
{
  printf 'setting\tmakespan_seconds\n'
  cat "$OUT/makespan.tsv"
  printf 'jobs1_call_hash\t%s\njobs2_call_hash\t%s\ncalls_identical\t%s\n' \
    "$H1" "$H2" "$([ "$H1" = "$H2" ] && echo true || echo false)"
} > "$OUT/summary.tsv"
cat "$OUT/summary.tsv"
[ "$H1" = "$H2" ] || { echo "FAIL: concurrency changed the calls" >&2; exit 1; }
echo "jobs benchmark: calls identical at JOBS=1 and JOBS=2"
