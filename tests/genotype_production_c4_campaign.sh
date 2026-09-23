#!/usr/bin/env bash
# The final C4 production campaigns, PREPARED BUT NOT LAUNCHED: held-out LZO, held-out synthetic
# mosaic and eligible-donor LOO, each at 2x, 5x, 10x and 30x, every call made by the one production
# command -- `panvar genotype --profile profiles/genotype_c4_campaign_local.tsv` (the science of
# profiles/genotype_c4.tsv with the local 3.5 GiB process limit) -- which is never told what kind of
# sample it is given. The campaign kind decides only the MODEL graph the call receives and the truth
# the evaluation compares against afterwards.
#
#   tests/genotype_production_c4_campaign.sh plan <out dir>
#       Writes <out>/jobs.tsv, <out>/settings.tsv and <out>/estimate.tsv, and checks every source input
#       it will need exists and is bound by content. Generates no reads, calls nothing, evaluates nothing.
#   JOBS=<j> THREADS=<t> tests/genotype_production_c4_campaign.sh launch <build dir> <out dir>
#       Runs the planned jobs through scripts/c4_production_campaign.py, JOBS at a time. REFUSES unless
#       GENOTYPE_C4_CAMPAIGN_APPROVED names the reviewed plan's jobs.tsv sha256: a campaign that can be
#       launched by accident is not frozen.
#   tests/genotype_production_c4_campaign.sh summarize <out dir>
#       One row per job from whatever has finished: status, decision, phase support, distance, floor,
#       excess, wall and peak memory.
#
# COHORTS, all from committed manifests, none re-selected here:
#   lzo     the 56 held-out donors of tests/genotype_development_donors.tsv, model = the full C4 graph,
#           reads = the frozen held-out 30x sources (tests/results/genotype_coverage_c4_held_out_abb35a0).
#           Their 2x/5x/10x subsets were never frozen as reads, and the selector that drew them at abb35a0
#           has since been retired, so they are regenerated -- genuinely missing -- by the current nested
#           selector UNSALTED (fnv1a64_fmix64_top53_v1), the development LOO campaign's rule for sources
#           drawn from the shared unsalted stream, and bound by content on creation. The package's own
#           retained counts are kept in package_kept_fragments as a record, not a gate. LZO and LOO of one
#           donor score identical subsets.
#   loo     the 55 eligible held-out donors of tests/genotype_loo_heldout_eligibility.tsv (HG00096, the
#           coordinate-reference donor, is ineligible), model = the C4 graph with the donor's two paths
#           removed (genotype_probe loo-gfa), the same reads as lzo, truth = the full C4 graph.
#   mosaic  the ten held-out plans of tests/genotype_c4_heldout_mosaic_manifest.tsv, model = the plan's
#           staged graph with the two synthetic haplotypes curated out of the panel, reads simulated and
#           nested by tests/genotype_c4_heldout_mosaic.sh (stages manifest stage verify reads nested),
#           truth = the staged graph's MOSAIC paths.
#
# TRUTH is never a call input: the call receives the model graph, reads, library and (mosaic only) the
# curation manifest. `launch` refuses any job whose call command names its truth file, except LZO,
# where the full graph is by definition both the model and the truth, and mosaic, where the synthetic
# paths are removed from the panel by curation and remain only for projection.
set -uo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MODE="${1:?usage: genotype_production_c4_campaign.sh plan <out> | launch <build> <out> | summarize <out>}"
shift
PROFILE="$REPO/profiles/genotype_c4_campaign_local.tsv"
SCIENCE_PROFILE="$REPO/profiles/genotype_c4.tsv"
GRAPH_GZ="$REPO/tests/real_data/c4.gfa.gz"
GRAPH_GZ_SHA=bf67be79db3ca05c2544162be0e18e4f0e3325872dd0d872407b577275870aa0
GRAPH_SHA=b13af390c28dd8b43a913c43f123bb1fcf639e0569f99783fa12a778820c8733
DONORS="$REPO/tests/genotype_development_donors.tsv"
ELIGIBILITY="$REPO/tests/genotype_loo_heldout_eligibility.tsv"
MOSAIC_MANIFEST="$REPO/tests/genotype_c4_heldout_mosaic_manifest.tsv"
HELDOUT_READS="$REPO/tests/results/genotype_coverage_c4_held_out_abb35a0"
TIMING="$REPO/tests/genotype_production_c4_timing.tsv"
REF='HG00096#1#haplotype1-0000024:31848049-32074186'
COVERAGES="2 5 10 30"
SOURCE_COVERAGE=30
sha() { if command -v sha256sum > /dev/null; then sha256sum "$1" | cut -d' ' -f1; else shasum -a 256 "$1" | cut -d' ' -f1; fi; }
key() { awk -F'\t' -v k="$2" '$1==k{print $2; exit}' "$1"; }

plan() {
  local out="$1" bad=0
  mkdir -p "$out"
  [ "$(sha "$GRAPH_GZ")" = "$GRAPH_GZ_SHA" ] || { echo "REFUSED: $GRAPH_GZ moved" >&2; exit 1; }
  local jobs="$out/jobs.tsv"
  printf "job\tcampaign\tsubject\tcoverage\tscheme\tdonor_label\tmodel_gfa\treads1\treads2\tlibrary\tcuration\ttruth_gfa\tsource_reads1\tsource_reads2\tsource_library\tpackage_kept_fragments\tsubsample_salt\tout_prefix\n" > "$jobs"
  local n=0
  # ---- lzo and loo: held-out donors, the frozen held-out read sources
  local donor cov src kept
  while read -r donor; do
    src="$HELDOUT_READS/$donor/source"
    if [ ! -f "$src/reads_1.fq.gz" ] || [ ! -f "$src/reads_2.fq.gz" ] || [ ! -f "$src/library.txt" ]; then
      echo "MISSING held-out source reads for $donor under $src" >&2; bad=1; continue
    fi
    for cov in $COVERAGES; do
      kept="$(key "$HELDOUT_READS/$donor/cov$cov/subsample.tsv" kept_fragments 2>/dev/null)"
      if [ "$cov" != "$SOURCE_COVERAGE" ] && [ -z "$kept" ]; then
        echo "MISSING recorded subsample for $donor at ${cov}x" >&2; bad=1
      fi
      for campaign in lzo loo; do
        if [ "$campaign" = loo ] && ! grep -q "^eligible"$'\t'"$donor\$" "$ELIGIBILITY"; then continue; fi
        local base="$out/$campaign/$donor"
        local model="$out/graphs/c4.gfa"
        [ "$campaign" = loo ] && model="$out/graphs/$donor.loo.gfa"
        n=$((n + 1))
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
          "$campaign.$donor.$cov" "$campaign" "$donor" "$cov" "$campaign" "$donor" "$model" \
          "$out/reads/$donor/cov$cov/reads_1.fq" "$out/reads/$donor/cov$cov/reads_2.fq" \
          "$out/reads/$donor/cov$cov/library.txt" "none" "$out/graphs/c4.gfa" "$src/reads_1.fq.gz" \
          "$src/reads_2.fq.gz" "$src/library.txt" "${kept:-all}" "none" "$base/call.$cov" >> "$jobs"
      done
    done
  done < <(awk -F'\t' '!/^#/ && $2=="held_out"{print $1}' "$DONORS")
  # ---- mosaic: the ten held-out plans, prepared by the frozen held-out mosaic harness
  local name
  while read -r name; do
    for cov in $COVERAGES; do
      n=$((n + 1))
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "mosaic.$name.$cov" mosaic "$name" "$cov" mosaic MOSAIC "$out/mosaic_inputs/$name/staged.gfa" \
        "$out/mosaic_inputs/$name/cov$cov/reads_1.fq" "$out/mosaic_inputs/$name/cov$cov/reads_2.fq" \
        "$out/mosaic_inputs/$name/cov$cov/library.txt" "$out/mosaic_inputs/mosaic_exclude.tsv" \
        "$out/mosaic_inputs/$name/staged.gfa" "harness" "harness" "harness" "harness" "$name" \
        "$out/mosaic/$name/call.$cov" >> "$jobs"
    done
  done < <(awk -F'\t' '$1=="plan" && $2!="name"{print $2}' "$MOSAIC_MANIFEST")

  # ---- settings, bound by content
  {
    printf "key\tvalue\n"
    printf "campaign\tc4-production-final\nlaunched\tfalse\n"
    printf "git_commit\t%s\ngit_tree_dirty\t%s\n" "$(git -C "$REPO" rev-parse HEAD)" \
      "$([ -z "$(git -C "$REPO" status --porcelain -- src include profiles tests/genotype_production_c4_campaign.sh)" ] && echo false || echo true)"
    printf "profile\t%s\nprofile_sha256\t%s\n" "profiles/genotype_c4_campaign_local.tsv" "$(sha "$PROFILE")"
    printf "science_profile\t%s\nscience_profile_sha256\t%s\n" "profiles/genotype_c4.tsv" "$(sha "$SCIENCE_PROFILE")"
    printf "graph\ttests/real_data/c4.gfa.gz\ngraph_gz_sha256\t%s\ngraph_sha256\t%s\n" "$GRAPH_GZ_SHA" "$GRAPH_SHA"
    printf "donor_manifest_sha256\t%s\neligibility_sha256\t%s\nmosaic_manifest_sha256\t%s\n" \
      "$(sha "$DONORS")" "$(sha "$ELIGIBILITY")" "$(sha "$MOSAIC_MANIFEST")"
    printf "reference\t%s\ncoverages\t%s\nsource_coverage\t%s\n" "$REF" "${COVERAGES// /,}" "$SOURCE_COVERAGE"
    printf "jobs\t%s\njobs_lzo\t%s\njobs_loo\t%s\njobs_mosaic\t%s\n" "$n" \
      "$(awk -F'\t' '$2=="lzo"' "$jobs" | wc -l | tr -d ' ')" "$(awk -F'\t' '$2=="loo"' "$jobs" | wc -l | tr -d ' ')" \
      "$(awk -F'\t' '$2=="mosaic"' "$jobs" | wc -l | tr -d ' ')"
    printf "jobs_sha256\t%s\n" "$(sha "$jobs")"
    printf "execution\tserial; one call at a time, each with the profile's 8 threads\n"
    printf "truth_policy\ttruth is an evaluation input only; no call command names a truth file except LZO (model = truth by definition) and mosaic (synthetic paths curated out of the panel)\n"
  } > "$out/settings.tsv"

  # ---- runtime and disk, from the production timings measured on this commit's regressions
  if [ -f "$TIMING" ]; then
    awk -F'\t' -v jobs="$jobs" '
      FNR==NR { if ($1 != "coverage" && $1 !~ /^#/) { wall[$1]=$2; disk[$1]=$3; fp[$1]=$4 } next }
      FNR>1 { count[$4]++ }
      END {
        for (c in count) { printf "%s\t%d\t%s\t%s\t%s\t%.2f\t%.0f\n", c, count[c], wall[c], disk[c], fp[c], count[c]*wall[c]/3600, count[c]*disk[c]; h+=count[c]*wall[c]; b+=count[c]*disk[c] }
        printf "all\t%d\tNA\tNA\tNA\t%.2f\t%.0f\n", FNR-1, h/3600, b
      }' "$TIMING" "$jobs" | sort -t$'\t' -k1,1n > "$out/estimate.body"
    { printf "coverage\tjobs\tseconds_per_job\tbytes_per_job\tpeak_footprint_bytes_per_job\tserial_hours\tbytes_total\n"
      grep -v '^all' "$out/estimate.body"; grep '^all' "$out/estimate.body"; } > "$out/estimate.tsv"
    rm -f "$out/estimate.body"
  fi
  echo "planned $n jobs under $out (not launched)"
  [ -f "$out/estimate.tsv" ] && column -t -s$'\t' "$out/estimate.tsv"
  [ "$bad" -eq 0 ] || { echo "plan INCOMPLETE: missing inputs above" >&2; exit 1; }
}

launch() {
  local build="$1" out="$2"
  local jobs="$out/jobs.tsv"
  [ -f "$jobs" ] || { echo "no plan at $out; run plan first" >&2; exit 1; }
  if [ "${GENOTYPE_C4_CAMPAIGN_APPROVED:-}" != "$(sha "$jobs")" ]; then
    echo "REFUSED: launching requires review of $jobs and GENOTYPE_C4_CAMPAIGN_APPROVED=<its sha256>." >&2
    exit 2
  fi
  : "${JOBS:?JOBS (concurrent calls) must be set by the scheduler calibration}"
  : "${THREADS:?THREADS (per call) must be set by the scheduler calibration}"
  # The resumable runner: classification, retries of external failures only, completeness after every job,
  # then evaluation, consolidation and the figure (scripts/c4_production_campaign.py).
  exec python3 "$REPO/scripts/c4_production_campaign.py" launch --build "$build" --out "$out" \
    --jobs "$JOBS" --threads "$THREADS"
}

summarize() {
  local out="$1" jobs="$1/jobs.tsv"
  printf "job\tstatus\tfinal_decision\tphase_support\tdistance\tfloor\texcess\twall_seconds\tpeak_footprint_bytes\n"
  tail -n +2 "$jobs" | while IFS=$'\t' read -r job _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ prefix; do
    local st="$prefix.mosaic_status.tsv" ev="$prefix.eval.summary.tsv"
    [ -f "$st" ] || continue
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$job" "$(key "$st" status)" "$(key "$st" final_decision)" \
      "$(key "$st" phase_support)" "$(key "$ev" unordered_distance_total 2>/dev/null)" \
      "$(key "$ev" availability_floor_total 2>/dev/null)" "$(key "$ev" excess_over_floor_total 2>/dev/null)" \
      "$(awk '/ real /{print $1}' "$prefix.time" 2>/dev/null)" "$(awk '/peak memory footprint/{print $1}' "$prefix.time" 2>/dev/null)"
  done
}

case "$MODE" in
  plan) plan "${1:?plan <out dir>}" ;;
  launch) launch "${1:?launch <build dir> <out dir>}" "${2:?launch <build dir> <out dir>}" ;;
  summarize) summarize "${1:?summarize <out dir>}" ;;
  *) echo "unknown mode $MODE" >&2; exit 2 ;;
esac
