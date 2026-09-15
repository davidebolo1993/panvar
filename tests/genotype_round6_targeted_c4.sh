#!/usr/bin/env bash
# Developer harness, not a regression test: the review-round-6 targeted rerun on C4 (steps 9-10).
#
#   BUILD=<clean build dir> tests/genotype_round6_targeted_c4.sh <frozen-cohort-dir> <out-dir>
#
# Reruns only the donors the round-5 held-out cohort got wrong or refused -- NA20355 (false depth refusal),
# HG01573, HG02492, HG04217 and NA19331 (truth never proposed) -- and HG02011, an exact control, on the SAME
# simulated reads as the frozen round-5 cohort, at state limit 1024 and K = 8 with a declared K = 256 phase
# audit. The frozen cohort directory is read, never written. Run it from a clean worktree.
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="${BUILD:-$REPO/build}"
FROZEN="${1:?usage: tests/genotype_round6_targeted_c4.sh <frozen-cohort-dir> <out-dir>}"
OUT="${2:?usage: tests/genotype_round6_targeted_c4.sh <frozen-cohort-dir> <out-dir>}"
DONORS="${DONORS:-NA20355 HG01573 HG02492 HG04217 NA19331 HG02011}"
STATE_LIMIT="${STATE_LIMIT:-1024}"
K="${K:-8}"
AUDIT="${AUDIT:-256}"
JOBS="${JOBS:-3}"
REF="HG00096#1#haplotype1-0000024:31848049-32074186"
PROBE="$BUILD/genotype_probe"

mkdir -p "$OUT"
if /usr/bin/time -l true >/dev/null 2>&1; then TIME_FLAG="-l"; else TIME_FLAG="-v"; fi
rss_of() {
  if [ "$TIME_FLAG" = "-l" ]; then awk '/maximum resident set size/{print $1}' "$1"
  else awk -F: '/Maximum resident set size/{gsub(/ /,"",$2); print $2*1024}' "$1"; fi
}
wall_of() {
  if [ "$TIME_FLAG" = "-l" ]; then awk '/ real /{print $1; exit}' "$1"
  else awk -F': ' '/Elapsed \(wall clock\)/{n=split($2,p,":"); s=0; for(i=1;i<=n;i++) s=s*60+p[i]; print s}' "$1"; fi
}

{
  printf 'key\tvalue\n'
  printf 'git_commit\t%s\n' "$(git -C "$REPO" rev-parse HEAD)"
  printf 'git_tree_dirty\t%s\n' "$(git -C "$REPO" status --porcelain | grep -q . && echo true || echo false)"
  printf 'probe_sha256\t%s\n' "$(shasum -a 256 "$PROBE" | cut -d' ' -f1)"
  printf 'script_sha256\t%s\n' "$(shasum -a 256 "${BASH_SOURCE[0]}" | cut -d' ' -f1)"
  printf 'frozen_cohort\t%s\nfrozen_cohort_commit\t%s\n' "$FROZEN" "$(awk -F'\t' '$1=="git_commit"{print $2}' "$FROZEN/signature.tsv")"
  printf 'graph_sha256\t%s\n' "$(shasum -a 256 "$FROZEN/c4.gfa" | cut -d' ' -f1)"
  printf 'state_limit\t%s\nmax_candidates\t%s\nphase_audit_candidates\t%s\ndonors\t%s\n' "$STATE_LIMIT" "$K" "$AUDIT" "$(echo $DONORS)"
} > "$OUT/signature.tsv"

run_donor() {
  local donor="$1" prefix="$OUT/$1"
  /usr/bin/time "$TIME_FLAG" "$PROBE" genotype --gfa "$FROZEN/c4.gfa" --reference "$REF" \
    --reads1 "$FROZEN/$donor/sim/reads_1.fq" --reads2 "$FROZEN/$donor/sim/reads_2.fq" \
    --library "$FROZEN/$donor/sim/library.txt" --out-prefix "$prefix" --state-limit "$STATE_LIMIT" \
    --prefixes "$K" --phase-audit "$AUDIT" > "$prefix.log" 2> "$prefix.time"
  printf 'peak_rss_bytes\t%s\nwall_seconds\t%s\n' "$(rss_of "$prefix.time")" "$(wall_of "$prefix.time")" > "$prefix.resources.tsv"
  if [ -f "$prefix.K$K.mosaic_status.tsv" ]; then
    "$PROBE" evaluate --scheme lzo --donor "$donor" --reference "$REF" --run-prefix "$prefix.K$K" \
      --model-gfa "$FROZEN/c4.gfa" --truth-gfa "$FROZEN/c4.gfa" --out "$prefix.K$K.eval" >> "$prefix.log" 2>&1
  fi
  echo "done $donor"
}
export -f run_donor rss_of wall_of
export OUT FROZEN PROBE REF STATE_LIMIT K AUDIT TIME_FLAG
printf '%s\n' $DONORS | xargs -P "$JOBS" -I{} bash -c 'run_donor "$@"' _ {}

{
  printf 'donor\tstatus\tmarker_termination\tmarker_converged\tpanel_baseline_complete\tpanel_diplotypes\tmosaics\ttruth_in_candidates\ttruth_source\tblocks_exact\tsv_blocks_exact\tunordered_distance\texcess\tphase_within_candidates\tphase_support\tphase_audit_status\tmin_block_GQ\twall_seconds\tpeak_rss_bytes\n'
  for donor in $DONORS; do
    p="$OUT/$donor.K$K"; st="$p.mosaic_status.tsv"; s="$p.eval.summary.tsv"
    [ -f "$st" ] || st="$OUT/$donor.mosaic_status.tsv"
    w() { { [ -f "$st" ] && awk -F'\t' -v key="$1" '$1==key{print $2}' "$st"; } || true; }
    v() { { [ -f "$s" ] && awk -F'\t' -v key="$1" '$1==key{print $2}' "$s"; } || true; }
    gq="$( { [ -f "$p.block_genotypes.tsv" ] && awk -F'\t' 'NR>1{if(min==""||$13+0<min+0)min=$13} END{print min}' "$p.block_genotypes.tsv"; } || true)"
    row=()
    for value in "$donor" "$(w status)" "$(w marker_depth_termination)" "$(w marker_depth_converged)" \
      "$(w panel_baseline_complete)" "$(w panel_diplotypes)" "$(w mosaic_candidates)" "$(v truth_in_candidates)" \
      "$(v truth_candidate_source)" "$(v blocks_exact)" "$(v sv_blocks_exact)" "$(v unordered_distance_total)" \
      "$(v excess_over_floor_total)" "$(w phase_support_within_candidates)" "$(w phase_support)" \
      "$(w phase_audit_status)" "$gq" "$(awk -F'\t' '$1=="wall_seconds"{print $2}' "$OUT/$donor.resources.tsv")" \
      "$(awk -F'\t' '$1=="peak_rss_bytes"{print $2}' "$OUT/$donor.resources.tsv")"; do
      row+=("${value:-NA}")
    done
    (IFS=$'\t'; echo "${row[*]}")
  done
} > "$OUT/targeted.tsv"
column -t -s $'\t' "$OUT/targeted.tsv"
