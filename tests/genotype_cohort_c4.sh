#!/usr/bin/env bash
# Developer harness, not a regression test: the held-out leave-zero-out C4 cohort at frozen settings
# (review round 5).
#
#   BUILD=<clean build dir> tests/genotype_cohort_c4.sh <out-dir>
#
# Settings come only from tests/genotype_frozen_settings.tsv, committed after the development ladder and
# before this runs. Donors are the manifest's held_out rows; development donors are never run here. Each
# donor is simulated once, genotyped once at the frozen (state limit, K) through the development prefix
# harness with a single prefix -- identical to a standalone run with that K -- and evaluated against its
# truth. Run it from a clean worktree so the signature names one commit.
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="${BUILD:-$REPO/build}"
OUT="${1:?usage: tests/genotype_cohort_c4.sh <out-dir>}"
MANIFEST="$REPO/tests/genotype_development_donors.tsv"
FROZEN="$REPO/tests/genotype_frozen_settings.tsv"
setting() { awk -F'\t' -v key="$1" '$1==key{print $2}' "$FROZEN"; }
STATE_LIMIT="$(setting state_limit)"
K="$(setting max_candidates)"
SEED="$(setting seed)"
DEPTH="$(setting depth)"
[ -n "$STATE_LIMIT" ] && [ -n "$K" ] && [ -n "$SEED" ] && [ -n "$DEPTH" ] || { echo "incomplete $FROZEN" >&2; exit 1; }
DONORS="${DONORS:-$(awk -F'\t' '$2=="held_out"{print $1}' "$MANIFEST")}"
JOBS="${JOBS:-4}"
REF="HG00096#1#haplotype1-0000024:31848049-32074186"
PROBE="$BUILD/genotype_probe"

mkdir -p "$OUT"
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
  printf 'state_limit\t%s\nmax_candidates\t%s\nseed\t%s\ndepth\t%s\nreference\t%s\ndonors\t%s\n' \
    "$STATE_LIMIT" "$K" "$SEED" "$DEPTH" "$REF" "$(echo $DONORS)"
} > "$OUT/signature.tsv"

run_donor() {
  local donor="$1" dir="$OUT/$1"
  mkdir -p "$dir/sim"
  "$PROBE" simulate --gfa "$OUT/c4.gfa" --reference "$REF" --donor "$donor" --depth "$DEPTH" \
    --seed "$SEED" --out "$dir/sim" > "$dir/sim.log" 2>&1
  local prefix="$dir/lzo"
  /usr/bin/time "$TIME_FLAG" "$PROBE" genotype --gfa "$OUT/c4.gfa" --reference "$REF" \
    --reads1 "$dir/sim/reads_1.fq" --reads2 "$dir/sim/reads_2.fq" --library "$dir/sim/library.txt" \
    --out-prefix "$prefix" --state-limit "$STATE_LIMIT" --prefixes "$K" > "$prefix.log" 2> "$prefix.time"
  printf 'peak_rss_bytes\t%s\nwall_seconds\t%s\n' "$(rss_of "$prefix.time")" "$(wall_of "$prefix.time")" \
    > "$prefix.resources.tsv"
  "$PROBE" evaluate --scheme lzo --donor "$donor" --reference "$REF" --run-prefix "$prefix.K$K" \
    --model-gfa "$OUT/c4.gfa" --truth-gfa "$OUT/c4.gfa" --out "$prefix.K$K.eval" >> "$prefix.log" 2>&1
  echo "done $donor"
}
export -f run_donor rss_of wall_of
export OUT PROBE REF SEED DEPTH STATE_LIMIT K TIME_FLAG

printf '%s\n' $DONORS | xargs -P "$JOBS" -I{} bash -c 'run_donor "$@"' _ {}

# ---- collect ------------------------------------------------------------------------------------------
COLUMNS_OUT="donor status note run_no_call unique found exhaustive target conditional candidate_scope global_certified truth_in_candidates truth_proposal_rank truth_score_rank phase_support phase_unsupported contenders blocks blocks_exact sv_blocks sv_blocks_exact read_supported_blocks read_supported_blocks_exact block_unsupported unordered_distance sv_distance context_distance availability_floor excess read_supported_distance read_supported_excess certified_gap marker_termination wall_seconds peak_rss_bytes"
{
  echo "$COLUMNS_OUT" | tr ' ' '\t'
  for donor in $DONORS; do
    p="$OUT/$donor/lzo.K$K"; s="$p.eval.summary.tsv"; st="$p.mosaic_status.tsv"; r="$OUT/$donor/lzo.resources.tsv"
    v() { [ -f "$s" ] && awk -F'\t' -v key="$1" '$1==key{print $2}' "$s" | head -1 || true; }
    w() { [ -f "$st" ] && awk -F'\t' -v key="$1" '$1==key{print $2}' "$st" | head -1 || true; }
    x() { [ -f "$r" ] && awk -F'\t' -v key="$1" '$1==key{print $2}' "$r" || true; }
    na() { [ -n "$1" ] && printf '%s' "$1" || printf 'NA'; }
    row=()
    for value in "$donor" "$(w status)" "$(w note | tr '\t' ' ' | cut -c1-80)" "$(v run_no_call)" \
      "$(w unique_candidates)" "$(w unique_candidates_found)" "$(w proposal_search_exhaustive)" \
      "$(w proposal_target_reached)" "$(w candidate_set_conditional)" "$(w candidate_scope)" \
      "$(w global_certified)" "$(v truth_in_candidates)" "$(v truth_proposal_rank)" \
      "$(v truth_fragment_score_rank)" "$(v phase_support)" "$(v phase_unsupported)" \
      "$(v fragment_contenders)" "$(v blocks)" "$(v blocks_exact)" "$(v sv_blocks)" "$(v sv_blocks_exact)" \
      "$(v read_supported_blocks)" "$(v read_supported_blocks_exact)" "$(v block_unsupported)" \
      "$(v unordered_distance_total)" "$(v sv_unordered_distance_total)" "$(v context_unordered_distance_total)" \
      "$(v availability_floor_total)" "$(v excess_over_floor_total)" \
      "$(v read_supported_unordered_distance_total)" "$(v read_supported_excess_over_floor_total)" \
      "$(w certified_gap)" "$(w marker_depth_termination)" "$(x wall_seconds)" "$(x peak_rss_bytes)"; do
      row+=("$(na "$value")")
    done
    (IFS=$'\t'; echo "${row[*]}")
  done
} > "$OUT/cohort.tsv"

# Aggregate: run no-calls are counted over every donor; block accuracy only over emitted genotypes, and
# again over read-supported genotypes. Unsupported blocks and prior-only phase are not no-calls.
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) c[$i]=i; next}
  { donors++
    if ($c["run_no_call"]!="false") { nocall++; next }
    called++; blocks+=$c["blocks"]; exact+=$c["blocks_exact"]; sv+=$c["sv_blocks"]; svexact+=$c["sv_blocks_exact"]
    rs+=$c["read_supported_blocks"]; rsexact+=$c["read_supported_blocks_exact"]; unsupported+=$c["block_unsupported"]
    dist+=$c["unordered_distance"]; svd+=$c["sv_distance"]; ctxd+=$c["context_distance"]
    flo+=$c["availability_floor"]; exc+=$c["excess"]; if ($c["excess"]>0) withexcess++
    if ($c["phase_unsupported"]=="true") phaseu++
    if ($c["truth_in_candidates"]=="true") truthin++
    if ($c["global_certified"]=="true") global++ }
  END {
    printf "key\tvalue\ndonors\t%d\nrun_no_call\t%d\ncalled\t%d\n", donors, nocall, called
    printf "blocks_emitted\t%d\nblocks_exact\t%d\nsv_blocks\t%d\nsv_blocks_exact\t%d\n", blocks, exact, sv, svexact
    printf "read_supported_blocks\t%d\nread_supported_blocks_exact\t%d\nblock_unsupported\t%d\n", rs, rsexact, unsupported
    printf "unordered_distance_total\t%d\nsv_distance_total\t%d\ncontext_distance_total\t%d\n", dist, svd, ctxd
    printf "availability_floor_total\t%d\nexcess_total\t%d\ndonors_with_excess\t%d\n", flo, exc, withexcess
    printf "phase_unsupported\t%d\ntruth_in_candidates\t%d\nglobal_certified\t%d\n", phaseu, truthin, global }' \
  "$OUT/cohort.tsv" > "$OUT/cohort.summary.tsv"
column -t -s $'\t' "$OUT/cohort.summary.tsv"
