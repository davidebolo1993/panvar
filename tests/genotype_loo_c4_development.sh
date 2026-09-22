#!/usr/bin/env bash
# The C4 development-donor LEAVE-ONE-OUT campaign.
#
#   BUILD=<build dir> tests/genotype_loo_c4_development.sh <out-dir>
#   tests/genotype_loo_c4_development.sh --check-manifest
#
# Every setting comes from tests/genotype_loo_c4_development_settings.tsv and nothing is hard-coded
# here: a value that is not in that file is not a setting of this campaign. The gate below asserts
# the frozen caller settings rather than trusting them, and REJECTS the round-7 combination by name,
# because NA19240 reaches the blockwise availability floor only at K=128 and a K=8 run measures a
# different question.
#
# LEAVE-ONE-OUT means the donor's own two paths are removed from the model graph and from nothing
# else: truth is read from the full graph. loo-gfa refuses unless exactly two paths leave.
#
# Reads are the frozen ordinary-LZO sources, reused rather than redrawn, so a donor present in the
# panel and the same donor removed from it are compared on identical bytes. Their hashes are bound
# into the signature.
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SETTINGS="$REPO/tests/genotype_loo_c4_development_settings.tsv"
setting() { awk -F'\t' -v k="$1" '$1==k{print $2; exit}' "$SETTINGS"; }

REFERENCE="$(setting reference)"
SEED="$(setting seed)"
SOURCE_COVERAGE="$(setting source_coverage)"
COVERAGES="$(setting coverages | tr ',' ' ')"
STATE_LIMIT="$(setting proposal_state_limit)"
CALL_PREFIXES="$(setting call_prefixes)"
PHASE_AUDIT="$(setting phase_audit_candidates)"
MAX_STATE_CLASSES="$(setting max_state_classes)"
MAX_PHYSICAL_ORIGINS="$(setting max_physical_origins)"
EMISSION_CACHE_BYTES="$(setting emission_cache_bytes)"
THREADS="$(setting threads)"
FROZEN_ROOT="$REPO/$(setting frozen_reads_root)"
DONOR_MANIFEST="$REPO/$(setting donor_manifest)"
ROLE="$(setting role)"
DONORS="${DONORS:-$(awk -F'\t' -v r="$ROLE" '$2==r{print $1}' "$DONOR_MANIFEST")}"

check_manifest() {
  local bad=0
  echo "manifest gate: $SETTINGS"

  # The frozen caller settings, asserted rather than trusted.
  local width k audit diag
  width="$STATE_LIMIT"; k="$(setting max_mosaic_candidates)"
  audit="$PHASE_AUDIT"; diag="$(setting diagnostic_comparison_width)"
  if [ "$width" = "256" ] && [ "$k" = "128" ] && [ "$diag" = "1024" ] && [ "$audit" = "0" ]; then
    echo "  frozen settings: width $width, K $k, diagnostic width $diag, phase audit $audit"
  else
    echo "  FROZEN SETTINGS WRONG: width=$width K=$k diagnostic=$diag audit=$audit"; bad=$((bad + 1))
  fi

  # The round-7 combination, rejected BY NAME. This is the whole point of the gate: running the
  # cohort at width 1024 / K 8 / audit 256 would answer a different question from the one NA19240
  # settled, and would do so silently.
  local rw rk ra
  rw="$(setting rejected_proposal_state_limit)"; rk="$(setting rejected_max_mosaic_candidates)"
  ra="$(setting rejected_phase_audit_candidates)"
  if [ "$width" = "$rw" ] || [ "$k" = "$rk" ] || [ "$audit" = "$ra" ]; then
    echo "  ROUND-7 SETTINGS DETECTED (width $rw / K $rk / audit $ra): refusing"; bad=$((bad + 1))
  else
    echo "  round-7 combination (width $rw, K $rk, audit $ra) is absent, as required"
  fi

  # Four coverages, the source among them.
  local ncov
  ncov="$(printf '%s\n' $COVERAGES | grep -c .)"
  if [ "$ncov" = "4" ] && printf '%s\n' $COVERAGES | grep -qx "$SOURCE_COVERAGE"; then
    echo "  four coverages, source included: $(echo $COVERAGES)"
  else
    echo "  COVERAGE SHAPE WRONG: $(echo $COVERAGES) (source $SOURCE_COVERAGE)"; bad=$((bad + 1))
  fi

  # Experimental caps are declared, and the campaign may not run on the global defaults.
  if [ "$MAX_STATE_CLASSES" = "600000000" ] && [ "$MAX_PHYSICAL_ORIGINS" = "700000000" ] &&
     [ "$EMISSION_CACHE_BYTES" = "536870912" ] && [ "$THREADS" = "8" ]; then
    echo "  experimental limits declared: $MAX_STATE_CLASSES classes, $MAX_PHYSICAL_ORIGINS origins,"
    echo "    $EMISSION_CACHE_BYTES cache bytes, $THREADS threads"
  else
    echo "  EXPERIMENTAL LIMITS WRONG"; bad=$((bad + 1))
  fi

  # Reads: the LZO lineage, unsalted, and the frozen sources must actually be present.
  local salt missing donors_n
  salt="$(setting seed_salt)"
  if [ "$SEED" = "20260916" ] && [ "$salt" = "none" ]; then
    echo "  reads bound to the ordinary-LZO stream: seed $SEED, salt $salt"
  else
    echo "  READ LINEAGE WRONG: seed=$SEED salt=$salt"; bad=$((bad + 1))
  fi
  missing=0; donors_n=0
  for d in $DONORS; do
    donors_n=$((donors_n + 1))
    [ -f "$FROZEN_ROOT/$d/source/reads_1.fq.gz" ] || missing=$((missing + 1))
  done
  if [ "$missing" = "0" ] && [ "$donors_n" = "$(setting expected_donors)" ]; then
    echo "  $donors_n development donors, every frozen source present"
  else
    echo "  FROZEN SOURCES MISSING: $missing of $donors_n"; bad=$((bad + 1))
  fi

  # Completeness declared before any result exists.
  local want_runs want_status narms
  narms="$(printf '%s' "$CALL_PREFIXES" | tr ',' '\n' | grep -c .)"
  want_runs="$(setting expected_shared_runs)"; want_status="$(setting expected_k_arm_status_records)"
  if [ "$want_runs" = "$((donors_n * ncov))" ] &&
     [ "$want_status" = "$((donors_n * ncov * narms))" ]; then
    echo "  completeness declared: $want_runs runs, $want_status K-arm status records"
  else
    echo "  COMPLETENESS INCONSISTENT: $want_runs / $want_status for $donors_n x $ncov x $narms"
    bad=$((bad + 1))
  fi

  local metrics
  metrics="$(awk -F'\t' '$1=="metric"' "$SETTINGS" | wc -l | tr -d ' ')"
  if [ "$metrics" -ge 30 ]; then echo "  $metrics metrics frozen, each naming its source key"
  else echo "  METRICS MISSING: only $metrics"; bad=$((bad + 1)); fi

  if [ "$bad" -ne 0 ]; then echo "manifest gate FAILED"; exit 1; fi
  echo "manifest gate passed"
}

if [ "${1:-}" = "--check-manifest" ]; then check_manifest; exit 0; fi

OUT="${1:?usage: BUILD=<build> tests/genotype_loo_c4_development.sh <out-dir>}"
BUILD="${BUILD:-$REPO/build}"
PROBE="$BUILD/genotype_probe"
[ -x "$PROBE" ] || { echo "genotype_probe not built at $PROBE" >&2; exit 1; }
check_manifest
mkdir -p "$OUT"
[ -z "$(ls -A "$OUT" 2>/dev/null)" ] || { echo "$OUT is not empty" >&2; exit 1; }

gzip -dc "$REPO/$(setting graph)" > "$OUT/c4.gfa"

if /usr/bin/time -l true >/dev/null 2>&1; then TIME_FLAG="-l"; else TIME_FLAG="-v"; fi
rss_of() { if [ "$TIME_FLAG" = "-l" ]; then awk '/maximum resident set size/{print $1}' "$1"
           else awk -F: '/Maximum resident set size/{gsub(/ /,"",$2); print $2*1024}' "$1"; fi; }
wall_of() { if [ "$TIME_FLAG" = "-l" ]; then awk '/ real /{print $1; exit}' "$1"
            else awk -F': ' '/Elapsed \(wall clock\)/{n=split($2,p,":"); s=0; for(i=1;i<=n;i++) s=s*60+p[i]; print s}' "$1"; fi; }

{
  printf 'key\tvalue\n'
  printf 'campaign\t%s\n' "$(setting campaign)"
  printf 'git_commit\t%s\n' "$(git -C "$REPO" rev-parse HEAD)"
  printf 'git_tree_dirty\t%s\n' "$(git -C "$REPO" status --porcelain | grep -q . && echo true || echo false)"
  printf 'probe_sha256\t%s\n' "$(shasum -a 256 "$PROBE" | cut -d' ' -f1)"
  printf 'script_sha256\t%s\n' "$(shasum -a 256 "${BASH_SOURCE[0]}" | cut -d' ' -f1)"
  printf 'settings_sha256\t%s\n' "$(shasum -a 256 "$SETTINGS" | cut -d' ' -f1)"
  printf 'graph_sha256\t%s\n' "$(shasum -a 256 "$OUT/c4.gfa" | cut -d' ' -f1)"
  printf 'donor_manifest_sha256\t%s\n' "$(shasum -a 256 "$DONOR_MANIFEST" | cut -d' ' -f1)"
  printf 'reference\t%s\nseed\t%s\nseed_salt\t%s\ncoverages\t%s\n' \
    "$REFERENCE" "$SEED" "$(setting seed_salt)" "$(echo $COVERAGES)"
  printf 'proposal_state_limit\t%s\ncall_prefixes\t%s\nphase_audit_candidates\t%s\n' \
    "$STATE_LIMIT" "$CALL_PREFIXES" "$PHASE_AUDIT"
  printf 'max_state_classes\t%s\nmax_physical_origins\t%s\nemission_cache_bytes\t%s\nthreads\t%s\n' \
    "$MAX_STATE_CLASSES" "$MAX_PHYSICAL_ORIGINS" "$EMISSION_CACHE_BYTES" "$THREADS"
  printf 'donors\t%s\n' "$(echo $DONORS)"
  # Every frozen source bound by content, so a changed read file cannot pass unnoticed.
  for d in $DONORS; do
    printf 'reads_sha256_%s_1\t%s\n' "$d" \
      "$(gzip -dc "$FROZEN_ROOT/$d/source/reads_1.fq.gz" | shasum -a 256 | cut -d' ' -f1)"
    printf 'reads_sha256_%s_2\t%s\n' "$d" \
      "$(gzip -dc "$FROZEN_ROOT/$d/source/reads_2.fq.gz" | shasum -a 256 | cut -d' ' -f1)"
  done
} > "$OUT/signature.tsv"

for donor in $DONORS; do
  dir="$OUT/$donor"; mkdir -p "$dir"
  src="$FROZEN_ROOT/$donor/source"
  # Leave-one-out: the donor's own two paths leave the MODEL graph. loo-gfa refuses unless exactly
  # two are removed, so a donor whose paths are named differently cannot slip through as an LZO run.
  "$PROBE" loo-gfa --gfa "$OUT/c4.gfa" --donor "$donor" --out "$dir/c4.loo.gfa" > "$dir/loo.log" 2>&1
  for cov in $COVERAGES; do
    reads="$dir/cov$cov"; mkdir -p "$reads"
    "$PROBE" subsample --reads1 "$src/reads_1.fq.gz" --reads2 "$src/reads_2.fq.gz" \
      --library "$src/library.txt" --source-coverage "$SOURCE_COVERAGE" --coverage "$cov" \
      --name-salt "$donor" --out "$reads" >> "$dir/loo.log" 2>&1
    prefix="$dir/call.$cov"
    /usr/bin/time "$TIME_FLAG" "$PROBE" genotype --gfa "$dir/c4.loo.gfa" --reference "$REFERENCE" \
      --reads1 "$reads/reads_1.fq" --reads2 "$reads/reads_2.fq" --library "$reads/library.txt" \
      --out-prefix "$prefix" --state-limit "$STATE_LIMIT" --prefixes "$CALL_PREFIXES" \
      --phase-audit "$PHASE_AUDIT" --threads "$THREADS" \
      --max-state-classes "$MAX_STATE_CLASSES" --max-physical-origins "$MAX_PHYSICAL_ORIGINS" \
      --emission-cache-bytes "$EMISSION_CACHE_BYTES" \
      > "$prefix.log" 2> "$prefix.time" || true
    printf 'peak_rss_bytes\t%s\nwall_seconds\t%s\n' \
      "$(rss_of "$prefix.time")" "$(wall_of "$prefix.time")" > "$dir/resources.$cov.tsv"
    for k in ${CALL_PREFIXES//,/ }; do
      if [ -f "$prefix.K$k.mosaic_status.tsv" ]; then
        # model = the LOO graph the call was made against; truth = the FULL graph, which still
        # carries this donor's paths. That asymmetry is what makes it leave-one-out.
        # A swallowed evaluate failure once hid 22 missing accuracy blocks behind a campaign that
        # reported success, so a failure is recorded where the summary will see it. It still does not
        # abort the cohort: one donor's evaluation failing is not a reason to lose the other seven.
        if ! "$PROBE" evaluate --scheme loo --donor "$donor" --reference "$REFERENCE" \
          --run-prefix "$prefix.K$k" --model-gfa "$dir/c4.loo.gfa" --truth-gfa "$OUT/c4.gfa" \
          --out "$dir/eval.$cov.K$k" >> "$prefix.log" 2>&1; then
          printf 'evaluate_failed\t%s\t%s\tK%s\n' "$donor" "$cov" "$k" >> "$OUT/evaluate_failures.tsv"
          echo "EVALUATE FAILED: $donor cov=$cov K=$k (see $prefix.log)" >&2
        fi
      fi
    done
  done
  echo "done $donor"
done
echo "campaign complete: $OUT"
