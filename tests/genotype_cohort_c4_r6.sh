#!/usr/bin/env bash
# Developer harness, not a regression test: the corrected held-out leave-zero-out C4 cohort (review round 6).
#
#   BUILD=<clean build dir> tests/genotype_cohort_c4_r6.sh <out-dir>
#
# Settings come only from tests/genotype_frozen_settings_r6.tsv, committed after the round-6 development ladder
# and before this runs. Two arms are ranked from ONE proposal and ONE scoring pass per donor: the panel arm
# (K = 0, the complete panel baseline alone) and the mosaic arm (the frozen positive K), each with the declared
# phase audit. Donors are the manifest's held_out rows. The frozen round-5 cohort and the round-6 targeted outputs
# are never written. Run it from a clean worktree so the signature names one commit.
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="${BUILD:-$REPO/build}"
OUT="${1:?usage: tests/genotype_cohort_c4_r6.sh <out-dir>}"
MANIFEST="$REPO/tests/genotype_development_donors.tsv"
FROZEN="$REPO/tests/genotype_frozen_settings_r6.tsv"
setting() { awk -F'\t' -v key="$1" '$1==key{print $2}' "$FROZEN"; }
STATE_LIMIT="$(setting state_limit)"
PANEL_K="$(setting panel_arm_mosaic_candidates)"
MOSAIC_K="$(setting mosaic_arm_mosaic_candidates)"
AUDIT="$(setting phase_audit_candidates)"
SEED="$(setting seed)"
DEPTH="$(setting depth)"
[ -n "$STATE_LIMIT" ] && [ "$PANEL_K" = 0 ] && [ -n "$MOSAIC_K" ] && [ "$MOSAIC_K" -gt 0 ] && [ -n "$AUDIT" ] &&
  [ "$AUDIT" -gt "$MOSAIC_K" ] && [ -n "$SEED" ] && [ -n "$DEPTH" ] || { echo "incomplete $FROZEN" >&2; exit 1; }
DONORS="${DONORS:-$(awk -F'\t' '$2=="held_out"{print $1}' "$MANIFEST")}"
JOBS="${JOBS:-4}"
REF="HG00096#1#haplotype1-0000024:31848049-32074186"
PROBE="$BUILD/genotype_probe"

mkdir -p "$OUT"
[ -z "$(ls -A "$OUT" 2>/dev/null)" ] || { echo "$OUT is not empty; the corrected cohort writes a new directory" >&2; exit 1; }
gzip -dc "$REPO/tests/real_data/c4.gfa.gz" > "$OUT/c4.gfa"

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
  printf 'manifest_sha256\t%s\n' "$(shasum -a 256 "$MANIFEST" | cut -d' ' -f1)"
  printf 'frozen_settings_sha256\t%s\n' "$(shasum -a 256 "$FROZEN" | cut -d' ' -f1)"
  printf 'graph_sha256\t%s\n' "$(shasum -a 256 "$OUT/c4.gfa" | cut -d' ' -f1)"
  printf 'state_limit\t%s\npanel_arm_K\t%s\nmosaic_arm_K\t%s\nphase_audit\t%s\nseed\t%s\ndepth\t%s\nreference\t%s\ndonors\t%s\n' \
    "$STATE_LIMIT" "$PANEL_K" "$MOSAIC_K" "$AUDIT" "$SEED" "$DEPTH" "$REF" "$(echo $DONORS)"
} > "$OUT/signature.tsv"

run_donor() {
  local donor="$1" dir="$OUT/$1"
  mkdir -p "$dir/sim"
  "$PROBE" simulate --gfa "$OUT/c4.gfa" --reference "$REF" --donor "$donor" --depth "$DEPTH" \
    --seed "$SEED" --out "$dir/sim" > "$dir/sim.log" 2>&1
  local prefix="$dir/lzo"
  /usr/bin/time "$TIME_FLAG" "$PROBE" genotype --gfa "$OUT/c4.gfa" --reference "$REF" \
    --reads1 "$dir/sim/reads_1.fq" --reads2 "$dir/sim/reads_2.fq" --library "$dir/sim/library.txt" \
    --out-prefix "$prefix" --state-limit "$STATE_LIMIT" --prefixes "$PANEL_K,$MOSAIC_K" --phase-audit "$AUDIT" \
    > "$prefix.log" 2> "$prefix.time"
  printf 'peak_rss_bytes\t%s\nwall_seconds\t%s\n' "$(rss_of "$prefix.time")" "$(wall_of "$prefix.time")" > "$prefix.resources.tsv"
  for k in "$PANEL_K" "$MOSAIC_K"; do
    if [ -f "$prefix.K$k.mosaic_status.tsv" ]; then
      "$PROBE" evaluate --scheme lzo --donor "$donor" --reference "$REF" --run-prefix "$prefix.K$k" \
        --model-gfa "$OUT/c4.gfa" --truth-gfa "$OUT/c4.gfa" --out "$prefix.K$k.eval" >> "$prefix.log" 2>&1
    fi
  done
  echo "done $donor"
}
export -f run_donor rss_of wall_of
export OUT PROBE REF SEED DEPTH STATE_LIMIT PANEL_K MOSAIC_K AUDIT TIME_FLAG
printf '%s\n' $DONORS | xargs -P "$JOBS" -I{} bash -c 'run_donor "$@"' _ {}

# ---- collect: one row per donor and arm ---------------------------------------------------------------
{
  printf 'donor\tarm\tK\tstatus\trun_no_call\tmarker_termination\tpanel_baseline_complete\tmosaic_candidates\tmosaic_target_reached\tcandidate_scope\ttruth_in_candidates\ttruth_source\tblocks\tblocks_exact\tsv_blocks\tsv_blocks_exact\tread_supported_blocks\tread_supported_blocks_exact\tunordered_distance\tsv_distance\tcontext_distance\texcess\tphase_within_candidates\tphase_support\tphase_audit_status\tmin_block_GQ\tnumeric_separation_gap\twinner_hash\twall_seconds\tpeak_rss_bytes\n'
  for donor in $DONORS; do
    for arm in panel mosaic; do
      k="$PANEL_K"; [ "$arm" = mosaic ] && k="$MOSAIC_K"
      p="$OUT/$donor/lzo.K$k"; st="$p.mosaic_status.tsv"; s="$p.eval.summary.tsv"; r="$OUT/$donor/lzo.resources.tsv"
      [ -f "$st" ] || st="$OUT/$donor/lzo.mosaic_status.tsv"
      w() { { [ -f "$st" ] && awk -F'\t' -v key="$1" '$1==key{print $2}' "$st" | head -1; } || true; }
      v() { { [ -f "$s" ] && awk -F'\t' -v key="$1" '$1==key{print $2}' "$s" | head -1; } || true; }
      x() { { [ -f "$r" ] && awk -F'\t' -v key="$1" '$1==key{print $2}' "$r"; } || true; }
      gq="$( { [ -f "$p.block_genotypes.tsv" ] && awk -F'\t' 'NR>1{if(min==""||$13+0<min+0)min=$13} END{print min}' "$p.block_genotypes.tsv"; } || true)"
      hash="$( { [ -f "$p.mosaic_haplotypes.fa" ] && grep -v '^>' "$p.mosaic_haplotypes.fa" | shasum -a 256 | cut -c1-16; } || true)"
      no_call=true; [ "$(w status)" = COMPLETE ] && no_call=false
      row=()
      for value in "$donor" "$arm" "$k" "$(w status)" "$no_call" "$(w marker_depth_termination)" \
        "$(w panel_baseline_complete)" "$(w mosaic_candidates)" "$(w mosaic_target_reached)" "$(w candidate_scope)" \
        "$(v truth_in_candidates)" "$(v truth_candidate_source)" "$(v blocks)" "$(v blocks_exact)" "$(v sv_blocks)" \
        "$(v sv_blocks_exact)" "$(v read_supported_blocks)" "$(v read_supported_blocks_exact)" \
        "$(v unordered_distance_total)" "$(v sv_unordered_distance_total)" "$(v context_unordered_distance_total)" \
        "$(v excess_over_floor_total)" "$(w phase_support_within_candidates)" "$(w phase_support)" \
        "$(w phase_audit_status)" "$gq" "$(w numeric_separation_gap)" "$hash" "$(x wall_seconds)" "$(x peak_rss_bytes)"; do
        row+=("${value:-NA}")
      done
      (IFS=$'\t'; echo "${row[*]}")
    done
  done
} > "$OUT/cohort.tsv"

# ---- aggregate per arm, with block accuracy stratified by conditional GQ ------------------------------------
{
  printf 'arm\tquantity\tvalue\n'
  for arm in panel mosaic; do
    awk -F'\t' -v arm="$arm" 'NR==1{for(i=1;i<=NF;i++)c[$i]=i; next} $c["arm"]==arm {
        donors++
        if ($c["run_no_call"]=="true") { nocall++; next }
        called++; blocks+=$c["blocks"]; exact+=$c["blocks_exact"]; sv+=$c["sv_blocks"]; svexact+=$c["sv_blocks_exact"]
        rs+=$c["read_supported_blocks"]; rsexact+=$c["read_supported_blocks_exact"]
        dist+=$c["unordered_distance"]; exc+=$c["excess"]; if ($c["excess"]+0>0) withexcess++
        if ($c["excess"]+0==0) donorexact++
        if ($c["truth_in_candidates"]=="true") truthin++
        if ($c["mosaic_target_reached"]=="false") deficient++
        if ($c["phase_support"]=="READS_AUDITED") audited++
        if ($c["phase_support"]=="READS_WITHIN_CANDIDATES") within++ }
      END {
        printf "%s\tdonors\t%d\n%s\trun_no_call\t%d\n%s\tcalled\t%d\n", arm, donors, arm, nocall, arm, called
        printf "%s\tdonors_fully_exact\t%d\n%s\tblocks_emitted\t%d\n%s\tblocks_exact\t%d\n", arm, donorexact, arm, blocks, arm, exact
        printf "%s\tintended_blocks\t%d\n", arm, blocks + nocall * 11
        printf "%s\tsv_blocks_exact\t%d/%d\n%s\tread_supported_blocks_exact\t%d/%d\n", arm, svexact, sv, arm, rsexact, rs
        printf "%s\tunordered_distance_total\t%d\n%s\texcess_total\t%d\n%s\tdonors_with_excess\t%d\n", arm, dist, arm, exc, arm, withexcess
        printf "%s\ttruth_in_candidates\t%d\n%s\tmosaic_target_unreached\t%d\n", arm, truthin, arm, deficient
        printf "%s\tphase_reads_audited\t%d\n%s\tphase_reads_within_candidates\t%d\n", arm, audited, arm, within }' "$OUT/cohort.tsv"
  done
  # GQ strata over every emitted block: conditional GQ is not a calibrated error rate at 30x LZO; reported, not thresholded.
  for arm in panel mosaic; do
    for donor in $DONORS; do
      k="$PANEL_K"; [ "$arm" = mosaic ] && k="$MOSAIC_K"
      e="$OUT/$donor/lzo.K$k.eval.blocks.tsv"; b="$OUT/$donor/lzo.K$k.block_genotypes.tsv"
      [ -f "$e" ] && [ -f "$b" ] || continue
      paste <(awk -F'\t' 'NR>1{print $4"\t"$10}' "$e") <(awk -F'\t' 'NR>1{print $1"\t"$13}' "$b") | awk -v arm="$arm" -F'\t' '{print arm"\t"$1"\t"$2"\t"$3"\t"$4}'
    done
  done | awk -F'\t' '{ gq=$5; s = (gq=="inf" || gq+0>=100) ? "GQ>=100" : (gq+0>=50 ? "GQ50-100" : (gq+0>=20 ? "GQ20-50" : "GQ<20"))
      if ($2!=$4) { bad++; next } if ($3=="NA") next; n[$1"\t"s]++; if ($3+0==0) ok[$1"\t"s]++ }
    END { for (k in n) printf "%s\tblocks_exact_%s\t%d/%d\n", substr(k,1,index(k,"\t")-1), substr(k,index(k,"\t")+1), ok[k], n[k]
          if (bad) printf "all\tblock_join_mismatches\t%d\n", bad }'
} > "$OUT/cohort.summary.tsv"
column -t -s $'\t' "$OUT/cohort.summary.tsv"
