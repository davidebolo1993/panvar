#!/usr/bin/env bash
# Developer harness, not a regression test: the leave-zero-out state-limit x K ladder on the frozen C4
# development donors (review round 5; K = 0 panel-only control and a K = 256 phase audit added in round 6).
#
#   BUILD=<clean build dir> tests/genotype_ladder_c4.sh <out-dir>
#
# For each development donor in tests/genotype_development_donors.tsv, simulates reads once, then for each
# merged-state limit runs the proposal once with K = the largest prefix and ranks every K prefix of it
# (genotype_probe genotype --prefixes), each identical to a standalone run with that K. Every prefix is
# evaluated against the donor's truth. Held-out donors are never touched. Outputs land under the given
# directory; run it from a clean worktree so the signature names one commit.
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="${BUILD:-$REPO/build}"
OUT="${1:?usage: tests/genotype_ladder_c4.sh <out-dir>}"
MANIFEST="$REPO/tests/genotype_development_donors.tsv"
DONORS="${DONORS:-$(awk -F'\t' '$2=="development"{print $1}' "$MANIFEST")}"
STATE_LIMITS="${STATE_LIMITS:-256 1024 4096}"
PREFIXES="${PREFIXES:-0,8,16,32,64,128,256}"
PHASE_AUDIT="${PHASE_AUDIT:-256}"
SEED="${SEED:-7}"
DEPTH="${DEPTH:-30}"
JOBS="${JOBS:-4}"
REF="HG00096#1#haplotype1-0000024:31848049-32074186"
PROBE="$BUILD/genotype_probe"

mkdir -p "$OUT"
gzip -dc "$REPO/tests/real_data/c4.gfa.gz" > "$OUT/c4.gfa"

if /usr/bin/time -l true >/dev/null 2>&1; then TIME_FLAG="-l"; else TIME_FLAG="-v"; fi
rss_of() {   # bytes, from either BSD or GNU time output
  if [ "$TIME_FLAG" = "-l" ]; then awk '/maximum resident set size/{print $1}' "$1"
  else awk -F: '/Maximum resident set size/{gsub(/ /,"",$2); print $2*1024}' "$1"; fi
}
wall_of() {  # seconds
  if [ "$TIME_FLAG" = "-l" ]; then awk '/ real /{print $1; exit}' "$1"
  else awk -F': ' '/Elapsed \(wall clock\)/{n=split($2,p,":"); s=0; for(i=1;i<=n;i++) s=s*60+p[i]; print s}' "$1"; fi
}

{
  printf 'key\tvalue\n'
  printf 'git_commit\t%s\n' "$(git -C "$REPO" rev-parse HEAD)"
  printf 'git_tree_dirty\t%s\n' "$(git -C "$REPO" status --porcelain | grep -q . && echo true || echo false)"
  printf 'probe_sha256\t%s\n' "$(shasum -a 256 "$PROBE" | cut -d' ' -f1)"
  printf 'script_sha256\t%s\n' "$(shasum -a 256 "${BASH_SOURCE[0]}" | cut -d' ' -f1)"
  printf 'manifest_sha256\t%s\n' "$(shasum -a 256 "$MANIFEST" | cut -d' ' -f1)"
  printf 'graph_sha256\t%s\n' "$(shasum -a 256 "$OUT/c4.gfa" | cut -d' ' -f1)"
  printf 'donors\t%s\nstate_limits\t%s\nprefixes\t%s\nphase_audit\t%s\nseed\t%s\ndepth\t%s\nreference\t%s\n' \
    "$(echo $DONORS)" "$STATE_LIMITS" "$PREFIXES" "$PHASE_AUDIT" "$SEED" "$DEPTH" "$REF"
} > "$OUT/signature.tsv"

run_donor() {
  local donor="$1" dir="$OUT/$1"
  mkdir -p "$dir/sim" "$dir/lzo"
  "$PROBE" simulate --gfa "$OUT/c4.gfa" --reference "$REF" --donor "$donor" --depth "$DEPTH" \
    --seed "$SEED" --out "$dir/sim" > "$dir/sim.log" 2>&1
  for limit in $STATE_LIMITS; do
    local prefix="$dir/lzo/S$limit"
    /usr/bin/time "$TIME_FLAG" "$PROBE" genotype --gfa "$OUT/c4.gfa" --reference "$REF" \
      --reads1 "$dir/sim/reads_1.fq" --reads2 "$dir/sim/reads_2.fq" --library "$dir/sim/library.txt" \
      --out-prefix "$prefix" --state-limit "$limit" --prefixes "$PREFIXES" --phase-audit "$PHASE_AUDIT" > "$prefix.log" 2> "$prefix.time"
    printf 'peak_rss_bytes\t%s\nwall_seconds\t%s\n' "$(rss_of "$prefix.time")" "$(wall_of "$prefix.time")" \
      > "$prefix.resources.tsv"
    for k in ${PREFIXES//,/ }; do
      "$PROBE" evaluate --scheme lzo --donor "$donor" --reference "$REF" --run-prefix "$prefix.K$k" \
        --model-gfa "$OUT/c4.gfa" --truth-gfa "$OUT/c4.gfa" --out "$prefix.K$k.eval" >> "$prefix.log" 2>&1
    done
  done
  echo "done $donor"
}
export -f run_donor rss_of wall_of
export OUT PROBE REF SEED DEPTH STATE_LIMITS PREFIXES PHASE_AUDIT TIME_FLAG

printf '%s\n' $DONORS | xargs -P "$JOBS" -I{} bash -c 'run_donor "$@"' _ {}

# ---- collect ------------------------------------------------------------------------------------------
{
  printf 'donor\tstate_limit\tK\tstatus\tunique\tfound\tpruned\texhaustive\ttarget\tconditional\tcandidate_scope\tglobal_certified\tpanel_baseline_complete\tmarker_converged\ttruth_in_candidates\ttruth_proposal_rank\ttruth_score_rank\tphase_support\tphase_within_candidates\tphase_audit_status\tmin_block_GQ\tcontenders\tread_supported_blocks\tblocks_exact\tsv_blocks_exact\tread_supported_blocks_exact\tunordered_distance\tsv_distance\tcontext_distance\tavailability_floor\texcess\tnumeric_separation_gap\tmarker_termination\twinner_hash\tladder_wall_seconds\tladder_peak_rss_bytes\n'
  for donor in $DONORS; do
    for limit in $STATE_LIMITS; do
      prefix="$OUT/$donor/lzo/S$limit"
      for k in ${PREFIXES//,/ }; do
        p="$prefix.K$k"; s="$p.eval.summary.tsv"; st="$p.mosaic_status.tsv"
        v() { [ -f "$s" ] && awk -F'\t' -v key="$1" '$1==key{print $2}' "$s" || echo NA; }
        w() { [ -f "$st" ] && awk -F'\t' -v key="$1" '$1==key{print $2}' "$st" || echo NA; }
        gq() { { [ -f "$p.block_genotypes.tsv" ] && awk -F'\t' 'NR>1{if(min==""||$13+0<min+0)min=$13} END{print min}' "$p.block_genotypes.tsv"; } || echo NA; }
        hash="NA"; [ -f "$p.mosaic_haplotypes.fa" ] && hash="$(grep -v '^>' "$p.mosaic_haplotypes.fa" | shasum -a 256 | cut -c1-12)"
        printf '%s\t' "$donor" "$limit" "$k" "$(v status)" "$(w mosaic_candidates)" "$(w mosaic_candidates_found)" \
          "$(w proposal_states_pruned)" "$(w proposal_search_exhaustive)" "$(w mosaic_target_reached)" \
          "$(w candidate_set_conditional)" "$(w candidate_scope)" "$(w global_certified)" "$(w panel_baseline_complete)" "$(w marker_depth_converged)" \
          "$(v truth_in_candidates)" "$(v truth_proposal_rank)" "$(v truth_fragment_score_rank)" \
          "$(w phase_support)" "$(w phase_support_within_candidates)" "$(w phase_audit_status)" "$(gq)" "$(v fragment_contenders)" "$(v read_supported_blocks)" "$(v blocks_exact)" \
          "$(v sv_blocks_exact)" "$(v read_supported_blocks_exact)" "$(v unordered_distance_total)" \
          "$(v sv_unordered_distance_total)" "$(v context_unordered_distance_total)" \
          "$(v availability_floor_total)" "$(v excess_over_floor_total)" "$(w numeric_separation_gap)$(w certified_gap)" \
          "$(w marker_depth_termination)" "$hash" \
          "$(awk -F'\t' '$1=="wall_seconds"{print $2}' "$prefix.resources.tsv")"
        printf '%s\n' "$(awk -F'\t' '$1=="peak_rss_bytes"{print $2}' "$prefix.resources.tsv")"
      done
    done
  done
} > "$OUT/ladder.tsv"
column -t -s $'\t' "$OUT/ladder.tsv"
