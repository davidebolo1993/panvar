#!/usr/bin/env bash
# The final C4 production campaigns, PREPARED BUT NOT LAUNCHED: held-out LZO, held-out synthetic
# mosaic and eligible-donor LOO, each at 2x, 5x, 10x and 30x, every call made by the one production
# command -- `panvar genotype --profile profiles/genotype_c4_v6_campaign_local.tsv` -- which is never told
# what kind of sample it is given. The campaign kind decides only the MODEL graph the call receives and the
# truth the evaluation compares against afterwards.
#
# PROFILES. The production profile is the v6 campaign-local profile. Its science is profiles/genotype_c4_v4.tsv
# (tie branching, cap 16; the representation guard; 32 repair rounds): v4 campaign-local is v4 key for key
# apart from the process settings (name, threads, footprint and time limits), v5 adds exactly the reviewed
# campaign-local resource and exact-acceleration settings (tests/genotype_v5_profile_gate.sh), and v6 exactly
# the 768 MiB repair origin cache (tests/genotype_v6_profile_gate.sh). The plan refuses unless that whole chain
# holds. v6 sets no tie-branch seed cache (repair_tie_seed_cache_bytes stays 0: measured, not adopted). The v1 profiles (genotype_c4.tsv, genotype_c4_campaign_local.tsv) are the frozen bc5dd97 campaign's and
# are never used here.
#
#   tests/genotype_production_c4_campaign.sh plan <out dir> <production build dir>
#       Writes <out>/jobs.tsv, <out>/settings.tsv, <out>/estimate.tsv, <out>/source_inputs.tsv (every frozen held-out
#       source read file and library the plan will read: absolute path, sha256, bytes) and <out>/plan_manifest.tsv
#       (the sha256 and size of jobs.tsv, settings.tsv and source_inputs.tsv). The manifest's own sha256 is the
#       plan's approval token. Generates no reads, calls nothing, evaluates nothing.
#       <out> must be new (never a directory a campaign has run in) and named genotype_c4_production_<build
#       commit, 7 hex>. The build must be a clean tree at PRODUCTION_COMMIT, or at a descendant whose production
#       sources (src, include, CMakeLists.txt, external) are identical to it; its commit and binary hashes are
#       recorded. CONCURRENT_JOBS (default 1) and THREADS_PER_CALL (default 8) are bound into the plan; more
#       than one concurrent job needs SCHEDULER_CALIBRATION=<calibration.tsv> selecting it.
#   tests/genotype_production_c4_campaign.sh launch <build dir> <out dir>
#       Runs the planned jobs through scripts/c4_production_campaign.py with the plan's scheduler. REFUSES
#       unless GENOTYPE_C4_CAMPAIGN_APPROVED names the reviewed plan's manifest sha256 (a campaign that can be
#       launched by accident is not frozen). The runner then refuses any manifest component that changed, a binary,
#       profile, scheduler or campaign file that is not the plan's, and any frozen source that changed -- all before
#       anything is staged or run.
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
MODE="${1:?usage: genotype_production_c4_campaign.sh plan <out> <build> | launch <build> <out> | summarize <out>}"
shift
PROFILE_REL=profiles/genotype_c4_v6_campaign_local.tsv
SCIENCE_PROFILE_REL=profiles/genotype_c4_v4.tsv
PROFILE="$REPO/$PROFILE_REL"
SCIENCE_PROFILE="$REPO/$SCIENCE_PROFILE_REL"
# The reviewed production build: census-free tie-branch seeds (fefa6c0, accepted with the v6 sentinels), plus the
# exact early branch-cap stop and the bounded seed cache (off: no campaign profile sets it).
PRODUCTION_COMMIT=11fe68189e9c7175b3f426b6a9f6c0948ecd387b
PRODUCTION_SOURCES="src include CMakeLists.txt external"
OLD_CAMPAIGN="$REPO/tests/results/genotype_c4_production_bc5dd97/campaign"
GRAPH_GZ="$REPO/tests/real_data/c4.gfa.gz"
GRAPH_GZ_SHA=bf67be79db3ca05c2544162be0e18e4f0e3325872dd0d872407b577275870aa0
GRAPH_SHA=b13af390c28dd8b43a913c43f123bb1fcf639e0569f99783fa12a778820c8733
DONORS="$REPO/tests/genotype_development_donors.tsv"
ELIGIBILITY="$REPO/tests/genotype_loo_heldout_eligibility.tsv"
MOSAIC_MANIFEST="$REPO/tests/genotype_c4_heldout_mosaic_manifest.tsv"
HELDOUT_READS_CANONICAL="$REPO/tests/results/genotype_coverage_c4_held_out_abb35a0"
# A gate may plan against a mirror of the package (to mutate a source without touching the frozen one); the plan
# records which root it used, and source_inputs.tsv binds every file it names.
HELDOUT_READS="${GENOTYPE_C4_HELDOUT_READS:-$HELDOUT_READS_CANONICAL}"
TIMING="$REPO/tests/genotype_production_c4_timing.tsv"
REF='HG00096#1#haplotype1-0000024:31848049-32074186'
COVERAGES="2 5 10 30"
SOURCE_COVERAGE=30
sha() { if command -v sha256sum > /dev/null; then sha256sum "$1" | cut -d' ' -f1; else shasum -a 256 "$1" | cut -d' ' -f1; fi; }
key() { awk -F'\t' -v k="$2" '$1==k{print $2; exit}' "$1"; }

# The profile chain v4 -> v4 campaign-local -> v5 -> v6, each step exactly its reviewed differences.
profile_chain() {
  science() { grep -v '^#' "$1" | grep -v '^$' | grep -v -E '^(profile|threads|footprint_limit_bytes|time_limit_seconds)'$'\t' | sort; }
  diff <(science "$SCIENCE_PROFILE") <(science "$REPO/profiles/genotype_c4_v4_campaign_local.tsv") > /dev/null \
    || { echo "REFUSED: v4 campaign-local does not carry the science of $SCIENCE_PROFILE_REL" >&2; return 1; }
  bash "$REPO/tests/genotype_v5_profile_gate.sh" > /dev/null || { echo "REFUSED: the v5 profile gate fails" >&2; return 1; }
  bash "$REPO/tests/genotype_v6_profile_gate.sh" > /dev/null || { echo "REFUSED: the v6 profile gate fails" >&2; return 1; }
  [ -z "$(key "$PROFILE" repair_tie_seed_cache_bytes)" ] || { echo "REFUSED: $PROFILE_REL sets a seed cache" >&2; return 1; }
  local k
  for k in "repair_tie_branching 1" "repair_tie_max_branches 16" "representation_guard 1" "repair_rounds 32"; do
    [ "$(key "$SCIENCE_PROFILE" "${k% *}")" = "${k#* }" ] || { echo "REFUSED: $SCIENCE_PROFILE_REL ${k% *} is not ${k#* }" >&2; return 1; }
    [ "$(key "$PROFILE" "${k% *}")" = "${k#* }" ] || { echo "REFUSED: $PROFILE_REL ${k% *} is not ${k#* }" >&2; return 1; }
  done
}

# The production build: a clean tree whose production sources are PRODUCTION_COMMIT's. Prints its commit.
production_build() {
  local build="$1" tree head
  [ -x "$build/panvar" ] && [ -x "$build/genotype_probe" ] || { echo "REFUSED: $build has no panvar and genotype_probe" >&2; return 1; }
  tree="$(git -C "$build/.." rev-parse --show-toplevel 2>/dev/null)" || { echo "REFUSED: $build is not in a git tree" >&2; return 1; }
  [ -z "$(git -C "$tree" status --porcelain)" ] || { echo "REFUSED: the build tree $tree is not clean" >&2; return 1; }
  head="$(git -C "$tree" rev-parse HEAD)"
  if [ "$head" != "$PRODUCTION_COMMIT" ]; then
    git -C "$tree" merge-base --is-ancestor "$PRODUCTION_COMMIT" "$head" 2> /dev/null \
      || { echo "REFUSED: the build tree is at $head, not $PRODUCTION_COMMIT or a descendant" >&2; return 1; }
    # shellcheck disable=SC2086
    git -C "$tree" diff --quiet "$PRODUCTION_COMMIT" "$head" -- $PRODUCTION_SOURCES \
      || { echo "REFUSED: the build tree's production sources differ from $PRODUCTION_COMMIT" >&2; return 1; }
  fi
  echo "$head"
}

plan() {
  local out="$1" build="$2" bad=0 build_commit
  [ "$(sha "$GRAPH_GZ")" = "$GRAPH_GZ_SHA" ] || { echo "REFUSED: $GRAPH_GZ moved" >&2; exit 1; }
  profile_chain || exit 1
  build_commit="$(production_build "$build")" || exit 1
  # A fresh directory, named from the build: never the frozen v1 campaign, never one a campaign has run in.
  mkdir -p "$(dirname "$out")"
  local resolved; resolved="$(cd "$(dirname "$out")" && pwd)/$(basename "$out")"
  [ "$resolved" != "$OLD_CAMPAIGN" ] || { echo "REFUSED: $out is the frozen bc5dd97 campaign" >&2; exit 1; }
  out="$resolved"   # every path in jobs.tsv is absolute, so the runner may be started from anywhere
  [ "$(basename "$out")" = "genotype_c4_production_${build_commit:0:7}" ] \
    || { echo "REFUSED: the plan directory must be named genotype_c4_production_${build_commit:0:7}" >&2; exit 1; }
  local used
  for used in launch.tsv completeness.tsv inputs.tsv package lzo loo mosaic reads graphs mosaic_inputs; do
    [ ! -e "$out/$used" ] || { echo "REFUSED: $out already holds $used; a campaign directory is never reused" >&2; exit 1; }
  done
  local concurrent="${CONCURRENT_JOBS:-1}" threads="${THREADS_PER_CALL:-8}" calibration="none"
  if [ "$concurrent" != 1 ]; then
    [ -n "${SCHEDULER_CALIBRATION:-}" ] && [ "$(key "$SCHEDULER_CALIBRATION" selected)" = "JOBS=$concurrent THREADS=$threads" ] \
      || { echo "REFUSED: $concurrent concurrent jobs need a scheduler calibration selecting JOBS=$concurrent THREADS=$threads" >&2; exit 1; }
    calibration="$SCHEDULER_CALIBRATION $(sha "$SCHEDULER_CALIBRATION")"
  fi
  mkdir -p "$out"
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
    # Dirtiness of what the campaign itself uses: its profile, harness, runner, consolidation and manifests.
    # The binary comes from a separate clean build, whose commit the runner records at launch.
    printf "git_commit\t%s\ncampaign_files_dirty\t%s\n" "$(git -C "$REPO" rev-parse HEAD)" \
      "$([ -z "$(git -C "$REPO" status --porcelain -- profiles tests/genotype_production_c4_campaign.sh tests/genotype_v5_profile_gate.sh tests/genotype_v6_profile_gate.sh scripts/c4_production_campaign.py scripts/consolidate_c4_production.py tests/genotype_c4_heldout_mosaic.sh "$DONORS" "$ELIGIBILITY" "$MOSAIC_MANIFEST")" ] && echo false || echo true)"
    printf "production_profile\t%s\nproduction_profile_sha256\t%s\n" "$PROFILE_REL" "$(sha "$PROFILE")"
    printf "production_profile_name\t%s\n" "$(key "$PROFILE" profile)"
    printf "science_profile\t%s\nscience_profile_sha256\t%s\n" "$SCIENCE_PROFILE_REL" "$(sha "$SCIENCE_PROFILE")"
    printf "profile_chain\tv4 science = v4 campaign-local science; v5 gate; v6 gate; no seed cache: verified\n"
    local f
    for f in profiles/genotype_c4_v4_campaign_local.tsv profiles/genotype_c4_v5_campaign_local.tsv \
             tests/genotype_v5_profile_gate.sh tests/genotype_v6_profile_gate.sh; do
      printf "%s_sha256\t%s\n" "$f" "$(sha "$REPO/$f")"
    done
    printf "production_commit\t%s\nbuild_dir\t%s\nbuild_commit\t%s\nbuild_tree_clean\ttrue\n" \
      "$PRODUCTION_COMMIT" "$(cd "$build" && pwd)" "$build_commit"
    printf "build_production_sources\t%s identical to production_commit\n" "$PRODUCTION_SOURCES"
    printf "production_binary_sha256\t%s\nevaluator_binary_sha256\t%s\n" "$(sha "$build/panvar")" "$(sha "$build/genotype_probe")"
    for f in tests/genotype_production_c4_campaign.sh scripts/c4_production_campaign.py \
             scripts/consolidate_c4_production.py tests/genotype_c4_heldout_mosaic.sh; do
      printf "%s_sha256\t%s\n" "$f" "$(sha "$REPO/$f")"
    done
    printf "graph\ttests/real_data/c4.gfa.gz\ngraph_gz_sha256\t%s\ngraph_sha256\t%s\n" "$GRAPH_GZ_SHA" "$GRAPH_SHA"
    printf "donor_manifest_sha256\t%s\neligibility_sha256\t%s\nmosaic_manifest_sha256\t%s\n" \
      "$(sha "$DONORS")" "$(sha "$ELIGIBILITY")" "$(sha "$MOSAIC_MANIFEST")"
    printf "reference\t%s\ncoverages\t%s\nsource_coverage\t%s\n" "$REF" "${COVERAGES// /,}" "$SOURCE_COVERAGE"
    printf "jobs\t%s\njobs_lzo\t%s\njobs_loo\t%s\njobs_mosaic\t%s\n" "$n" \
      "$(awk -F'\t' '$2=="lzo"' "$jobs" | wc -l | tr -d ' ')" "$(awk -F'\t' '$2=="loo"' "$jobs" | wc -l | tr -d ' ')" \
      "$(awk -F'\t' '$2=="mosaic"' "$jobs" | wc -l | tr -d ' ')"
    printf "jobs_sha256\t%s\n" "$(sha "$jobs")"
    printf "heldout_reads_root\t%s\nheldout_reads_root_canonical\t%s\n" "$HELDOUT_READS" \
      "$([ "$HELDOUT_READS" = "$HELDOUT_READS_CANONICAL" ] && echo true || echo false)"
    printf "concurrent_jobs\t%s\nthreads_per_call\t%s\nscheduler_calibration\t%s\n" "$concurrent" "$threads" "$calibration"
    printf "execution\t%s concurrent call(s), %s threads each; each call may approach 3.5 GiB of footprint under a 4 GiB limit\n" "$concurrent" "$threads"
    printf "outcomes\tcompleted: strict_call, genotype_equivalent_call, no_call (tie_branching_incomplete, divergent_tie and every other declared no_call_reason); a no-call is never retried and never stops the campaign\n"
    printf "estimate_provenance\tHISTORICAL: tests/genotype_production_c4_timing.tsv was measured with the v1 profile (genotype_c4.tsv) at b754647 on an Apple M4; it is not a prediction for the production profile\n"
    printf "truth_policy\ttruth is an evaluation input only; no call command names a truth file except LZO (model = truth by definition) and mosaic (synthetic paths curated out of the panel)\n"
  } > "$out/settings.tsv"

  # ---- runtime and disk: HISTORICAL v1 timings (see estimate_provenance), labelled as such in the table
  if [ -f "$TIMING" ]; then
    awk -F'\t' -v jobs="$jobs" '
      FNR==NR { if ($1 != "coverage" && $1 !~ /^#/) { wall[$1]=$2; disk[$1]=$3; fp[$1]=$4 } next }
      FNR>1 { count[$4]++ }
      END {
        for (c in count) { printf "%s\t%d\t%s\t%s\t%s\t%.2f\t%.0f\n", c, count[c], wall[c], disk[c], fp[c], count[c]*wall[c]/3600, count[c]*disk[c]; h+=count[c]*wall[c]; b+=count[c]*disk[c] }
        printf "all\t%d\tNA\tNA\tNA\t%.2f\t%.0f\n", FNR-1, h/3600, b
      }' "$TIMING" "$jobs" | sort -t$'\t' -k1,1n > "$out/estimate.body"
    { printf "# HISTORICAL, NOT A PREDICTION FOR THE PRODUCTION PROFILE: per-job seconds, bytes and footprints measured with the v1 profile\n"
      printf "# (profiles/genotype_c4.tsv, 11 rounds, no tie branching, no representation guard) at b754647 on an Apple M4.\n"
      printf "coverage\tjobs\tseconds_per_job\tbytes_per_job\tpeak_footprint_bytes_per_job\tserial_hours\tbytes_total\n"
      grep -v '^all' "$out/estimate.body"; grep '^all' "$out/estimate.body"; } > "$out/estimate.tsv"
    rm -f "$out/estimate.body"
  fi
  # ---- the frozen sources, bound by content now; then the manifest over the plan's components
  {
    printf "path\tsha256\tbytes\trole\n"
    awk -F'\t' 'NR>1 && $13!="harness"{print $13"\treads_1"; print $14"\treads_2"; print $15"\tlibrary"}' "$jobs" | sort -u |
      while IFS=$'\t' read -r path role; do
        [ -f "$path" ] || { echo "MISSING source $path" >&2; continue; }
        printf "%s\t%s\t%s\t%s\n" "$path" "$(sha "$path")" "$(wc -c < "$path" | tr -d ' ')" "$role"
      done
  } > "$out/source_inputs.tsv"
  {
    printf "component\tsha256\tbytes\n"
    local c
    for c in jobs.tsv settings.tsv source_inputs.tsv; do
      printf "%s\t%s\t%s\n" "$c" "$(sha "$out/$c")" "$(wc -c < "$out/$c" | tr -d ' ')"
    done
  } > "$out/plan_manifest.tsv"
  echo "planned $n jobs under $out (not launched)"
  echo "source inputs: $(($(wc -l < "$out/source_inputs.tsv") - 1)); plan manifest sha256 (the approval token): $(sha "$out/plan_manifest.tsv")"
  [ -f "$out/estimate.tsv" ] && grep -v '^#' "$out/estimate.tsv" | column -t -s$'\t' && echo "(historical v1 timings; not a prediction for the production profile)"
  [ "$bad" -eq 0 ] || { echo "plan INCOMPLETE: missing inputs above" >&2; exit 1; }
}

launch() {
  local build="$1" out="$2"
  local jobs="$out/jobs.tsv"
  [ -f "$jobs" ] || { echo "no plan at $out; run plan first" >&2; exit 1; }
  [ -f "$out/plan_manifest.tsv" ] || { echo "REFUSED: no plan manifest at $out" >&2; exit 2; }
  if [ "${GENOTYPE_C4_CAMPAIGN_APPROVED:-}" != "$(sha "$out/plan_manifest.tsv")" ]; then
    echo "REFUSED: launching requires review of $out/plan_manifest.tsv and GENOTYPE_C4_CAMPAIGN_APPROVED=<its sha256>." >&2
    exit 2
  fi
  # The scheduler is the plan's; the runner refuses a binary or profile that is not.
  local jobs_n threads_n
  jobs_n="$(key "$out/settings.tsv" concurrent_jobs)"; threads_n="$(key "$out/settings.tsv" threads_per_call)"
  [ -n "$jobs_n" ] && [ -n "$threads_n" ] || { echo "REFUSED: $out/settings.tsv binds no scheduler" >&2; exit 2; }
  # The resumable runner: classification, retries of external failures only, completeness after every job,
  # then evaluation, consolidation and the figure (scripts/c4_production_campaign.py).
  exec python3 "$REPO/scripts/c4_production_campaign.py" launch --build "$build" --out "$out" \
    --jobs "$jobs_n" --threads "$threads_n" ${EVALUATOR:+--evaluator "$EVALUATOR"}
}

summarize() {
  local out="$1" jobs="$1/jobs.tsv"
  printf "job\tstatus\tcall_class\tno_call_reason\tfinal_decision\tphase_support\tdistance\tfloor\texcess\twall_seconds\tpeak_footprint_bytes\n"
  tail -n +2 "$jobs" | while IFS=$'\t' read -r job _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ prefix; do
    local st="$prefix.mosaic_status.tsv" ev="$prefix.eval.summary.tsv"
    [ -f "$st" ] || continue
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$job" "$(key "$st" status)" "$(key "$st" call_class)" \
      "$(key "$st" no_call_reason)" "$(key "$st" final_decision)" \
      "$(key "$st" phase_support)" "$(key "$ev" unordered_distance_total 2>/dev/null)" \
      "$(key "$ev" availability_floor_total 2>/dev/null)" "$(key "$ev" excess_over_floor_total 2>/dev/null)" \
      "$(awk '/ real /{print $1}' "$prefix.time" 2>/dev/null)" "$(awk '/peak memory footprint/{print $1}' "$prefix.time" 2>/dev/null)"
  done
}

case "$MODE" in
  plan) plan "${1:?plan <out dir> <production build dir>}" "${2:?plan <out dir> <production build dir>}" ;;
  launch) launch "${1:?launch <build dir> <out dir>}" "${2:?launch <build dir> <out dir>}" ;;
  summarize) summarize "${1:?summarize <out dir>}" ;;
  *) echo "unknown mode $MODE" >&2; exit 2 ;;
esac
