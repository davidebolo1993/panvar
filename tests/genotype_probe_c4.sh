#!/usr/bin/env bash
# Developer probe, not a regression test: the Phase 8 C4 protocol on one frozen development donor.
#
#   tests/genotype_probe_c4.sh [out-dir]
#
# Runs the leave-zero-out K ladder and a leave-one-out probe for the same donor and the same reads,
# then evaluates every run against the donor's truth. Outputs land under tests/results/ (gitignored).
# NA19240 is one of the four donors frozen as development/regression donors in
# archive/genotype-factor-2026-09-11:tests/genotype_baseline_cohort.sh, so the held-out cohort stays
# untouched.
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="${BUILD:-$REPO/build}"
OUT="${1:-$REPO/tests/results/genotype_probe_c4}"
DONOR="${DONOR:-NA19240}"
SEED="${SEED:-7}"
DEPTH="${DEPTH:-30}"
LADDER="${LADDER:-8 16 32 64 128 256}"
LOO_K="${LOO_K:-64 256}"
REF="HG00096#1#haplotype1-0000024:31848049-32074186"
PANVAR="$BUILD/panvar"
PROBE="$BUILD/genotype_probe"

mkdir -p "$OUT/sim" "$OUT/lzo" "$OUT/loo"
gzip -dc "$REPO/tests/real_data/c4.gfa.gz" > "$OUT/c4.gfa"

if /usr/bin/time -l true >/dev/null 2>&1; then TIME_FLAG="-l"; else TIME_FLAG="-v"; fi
rss_of() {   # bytes, from either BSD or GNU time output
  if [ "$TIME_FLAG" = "-l" ]; then awk '/maximum resident set size/{print $1}' "$1"
  else awk -F: '/Maximum resident set size/{gsub(/ /,"",$2); print $2*1024}' "$1"; fi
}

# ---- frozen checkpoint signature ------------------------------------------------------------------
{
  printf 'key\tvalue\n'
  printf 'git_commit\t%s\n' "$(git -C "$REPO" rev-parse HEAD)"
  printf 'git_tree_dirty\t%s\n' "$(git -C "$REPO" status --porcelain | grep -q . && echo true || echo false)"
  printf 'panvar_sha256\t%s\n' "$(shasum -a 256 "$PANVAR" | cut -d' ' -f1)"
  printf 'probe_sha256\t%s\n' "$(shasum -a 256 "$PROBE" | cut -d' ' -f1)"
  printf 'script_sha256\t%s\n' "$(shasum -a 256 "${BASH_SOURCE[0]}" | cut -d' ' -f1)"
  printf 'graph_sha256\t%s\n' "$(shasum -a 256 "$OUT/c4.gfa" | cut -d' ' -f1)"
  printf 'donor\t%s\nseed\t%s\ndepth\t%s\nladder\t%s\nloo_k\t%s\nreference\t%s\n' \
    "$DONOR" "$SEED" "$DEPTH" "$LADDER" "$LOO_K" "$REF"
} > "$OUT/signature.tsv"

# ---- simulation and the leave-one-out graph -------------------------------------------------------
"$PROBE" simulate --gfa "$OUT/c4.gfa" --reference "$REF" --donor "$DONOR" --depth "$DEPTH" \
  --seed "$SEED" --out "$OUT/sim"
"$PROBE" loo-gfa --gfa "$OUT/c4.gfa" --donor "$DONOR" --out "$OUT/c4.loo.gfa"

run_and_evaluate() {   # scheme graph K
  local scheme="$1" graph="$2" k="$3" prefix="$OUT/$1/K$3"
  set +e
  # the unrepaired caller this ladder was built on: explicit since the short-read preset became the default
  /usr/bin/time "$TIME_FLAG" "$PANVAR" genotype --adaptive-repair off --gfa "$graph" --reference-path "$REF" \
    --reads1 "$OUT/sim/reads_1.fq" --reads2 "$OUT/sim/reads_2.fq" --library "$OUT/sim/library.txt" \
    --out-prefix "$prefix" --max-mosaic-candidates "$k" -q 2> "$prefix.time"
  set -e
  printf 'peak_rss_bytes\t%s\n' "$(rss_of "$prefix.time")" > "$prefix.resources.tsv"
  "$PROBE" evaluate --scheme "$scheme" --donor "$DONOR" --reference "$REF" --run-prefix "$prefix" \
    --model-gfa "$graph" --truth-gfa "$OUT/c4.gfa" --out "$prefix.eval"
}

for k in $LADDER; do run_and_evaluate lzo "$OUT/c4.gfa" "$k"; done
for k in $LOO_K; do run_and_evaluate loo "$OUT/c4.loo.gfa" "$k"; done

# ---- collect ------------------------------------------------------------------------------------------
{
  printf 'scheme\tK\tstatus\tstate_limit\tunique\tfound\tpruned\texhaustive\ttarget\tconditional\ttruth_in_candidates\ttruth_proposal_rank\ttruth_score_rank\tphase_support\tcontenders\tgenotype_read_supported\tnumeric_separation_gap\tfragment_gap\tblocks_exact\tsv_blocks_exact\tunordered_distance\tavailability_floor\texcess\tphase_switch_blocks\ttruth_spellable\tproposal_seconds\tseconds\tpeak_rss_bytes\twinner_hash\n'
  for scheme in lzo loo; do
    ks="$LADDER"; [ "$scheme" = loo ] && ks="$LOO_K"
    for k in $ks; do
      p="$OUT/$scheme/K$k"; s="$p.eval.summary.tsv"
      v() { awk -F'\t' -v key="$1" '$1==key{print $2}' "$s"; }
      hash="NA"; [ -f "$p.mosaic_haplotypes.fa" ] && hash="$(grep -v '^>' "$p.mosaic_haplotypes.fa" | shasum -a 256 | cut -c1-12)"
      printf '%s\t' "$scheme" "$k" "$(v status)" "$(v proposal_state_limit)" "$(v mosaic_candidates)" \
        "$(v mosaic_candidates_found)" "$(v proposal_states_pruned)" "$(v proposal_search_exhaustive)" \
        "$(v mosaic_target_reached)" "$(v candidate_set_conditional)" "$(v truth_in_candidates)" \
        "$(v truth_proposal_rank)" "$(v truth_fragment_score_rank)" "$(v phase_support)" \
        "$(v fragment_contenders)" "$(v blocks_genotype_read_supported)" "$(v numeric_separation_gap)" \
        "$(v fragment_term_gap_to_runner_up)" "$(v blocks_exact)" "$(v sv_blocks_exact)" \
        "$(v unordered_distance_total)" "$(v availability_floor_total)" "$(v excess_over_floor_total)" \
        "$(v phase_switch_blocks)" "$(v truth_hap1_spellable)/$(v truth_hap2_spellable)" \
        "$(v seconds_proposal)" "$(v seconds_total)" "$(awk -F'\t' '{print $2}' "$p.resources.tsv")"
      printf '%s\n' "$hash"
    done
  done
} > "$OUT/ladder.tsv"
column -t -s $'\t' "$OUT/ladder.tsv"
