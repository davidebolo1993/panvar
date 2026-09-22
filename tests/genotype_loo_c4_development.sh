#!/usr/bin/env bash
# The C4 development-donor LEAVE-ONE-OUT campaign.
#
#   BUILD=<build dir> tests/genotype_loo_c4_development.sh <out-dir>
#   tests/genotype_loo_c4_development.sh --check-manifest
#
# Every setting comes from tests/genotype_loo_c4_development_settings.tsv and nothing is hard-coded
# here: a value that is not in that file is not a setting of this campaign. The gate below asserts the
# frozen caller settings rather than trusting them, and REJECTS the round-7 combination by name,
# because NA19240 reaches the blockwise availability floor only at K=128 and a K=8 run measures a
# different question.
#
# LEAVE-ONE-OUT means the donor's own two paths are removed from the model graph and from nothing
# else: truth is read from the full graph. loo-gfa refuses unless exactly two paths leave.
#
# THINNING IS UNSALTED. The frozen sources are one realization per donor from a shared, unsalted
# stream, so thinning them with a donor-salted selector would choose a different subset from the one
# ordinary LZO chose and destroy the paired comparison. subsample is therefore called WITHOUT
# --name-salt, every regenerated subset is gated against a frozen retained-fragment count, and at 30x
# the regenerated FASTQ must be byte-identical to the frozen source.
#
# NOTHING IS SWALLOWED. A failed call, a failed evaluation, a missing status and a legitimate no-call
# are four different things, recorded separately and reconciled at the end; the campaign exits nonzero
# unless every declared run is accounted for and every acceptance value is met. An earlier revision
# used `|| true` on evaluate and reported success while 22 of 38 accuracy blocks were missing.
#
# The functions below are sourced by tests/genotype_loo_c4_development_gates.sh. A gate that
# reimplements the thing it guards proves only that the copy agrees with itself.
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
# tests/results is gitignored (245 MB of frozen FASTQs), so a clean checkout -- a detached worktree
# for a reviewer-facing run, or CI -- has no read sources. The root is therefore overridable, the way
# the sibling harnesses take BUILD, DONORS and JOBS, so a run from a clean worktree can point at the
# frozen sources without the campaign's code being dirty. Whatever it points at is gated by content:
# every source hash is checked against the settings before a read is scored, and bound into the
# signature afterwards, so an override changes where the reads are read from and not which reads they are.
FROZEN_ROOT="${FROZEN_ROOT:-$REPO/$(setting frozen_reads_root)}"
DONOR_MANIFEST="$REPO/$(setting donor_manifest)"
ROLE="$(setting role)"
DONORS="${DONORS:-$(awk -F'\t' -v r="$ROLE" '$2==r{print $1}' "$DONOR_MANIFEST")}"

# ---- content binding ---------------------------------------------------------------------------
content_sha() { shasum -a 256 "$1" | cut -d' ' -f1; }
gz_content_sha() { gzip -dc "$1" | shasum -a 256 | cut -d' ' -f1; }

# The frozen expectation rows, addressed by name rather than by line number.
expected_source() { # <donor> <field>
  awk -F'\t' -v d="$2" -v f="$3" '$1=="source" && $2==d && $3==f{print $4; exit}' "$1"
}
expected_subset() { # <settings> <donor> <coverage>
  awk -F'\t' -v d="$2" -v c="$3" '$1=="subset" && $2==d && $3==c{print $4; exit}' "$1"
}
expected_lineage() { # <settings> <path>
  awk -F'\t' -v p="$2" '$1=="lineage_sha256" && $2==p{print $3; exit}' "$1"
}

# Item 3, half one: the frozen inputs are what they claim to be, before a single read is scored.
# Prints one FAIL line per mismatch and returns nonzero.
check_frozen_inputs() { # <settings> <frozen root> <donors...>
  local settings="$1" root="$2"; shift 2
  local bad=0 d got want
  for d in "$@"; do
    for pair in "reads_1.fq.gz reads1_sha256" "reads_2.fq.gz reads2_sha256"; do
      set -- $pair
      got="$(gz_content_sha "$root/$d/source/$1")"
      want="$(expected_source "$settings" "$d" "$2")"
      if [ -z "$want" ]; then echo "FAIL no frozen expectation for $d $2"; bad=1
      elif [ "$got" != "$want" ]; then echo "FAIL $d $1 content moved: $got != $want"; bad=1; fi
    done
    got="$(content_sha "$root/$d/source/library.txt")"
    want="$(expected_source "$settings" "$d" "library_sha256")"
    if [ "$got" != "$want" ]; then echo "FAIL $d library.txt moved: $got != $want"; bad=1; fi
  done
  return $bad
}

# Item 3, half two: the regenerated subset is the one the settings froze. The retained count pins the
# selector; the mate counts pin pairing; at the source coverage the bytes themselves must match, which
# is the strongest available statement that this campaign scores the reads ordinary LZO scored.
check_subset() { # <settings> <donor> <coverage> <subset dir> <frozen source dir>
  local settings="$1" d="$2" cov="$3" dir="$4" src="$5"
  local bad=0 kept want salt n1 n2
  kept="$(awk -F'\t' '$1=="kept_fragments"{print $2; exit}' "$dir/subsample.tsv")"
  want="$(expected_subset "$settings" "$d" "$cov")"
  salt="$(awk -F'\t' '$1=="name_salt"{print $2; exit}' "$dir/subsample.tsv")"
  if [ -z "$want" ]; then echo "FAIL no frozen subset count for $d ${cov}x"; bad=1
  elif [ "$kept" != "$want" ]; then echo "FAIL $d ${cov}x retained $kept fragments, frozen $want"; bad=1; fi
  # A salted thinning is the specific mistake this campaign was corrected for; it is named, not implied.
  if [ "$salt" != "none" ]; then
    echo "FAIL $d ${cov}x thinned with name_salt '$salt'; this campaign requires the unsalted selector"
    bad=1
  fi
  n1=$(( $(wc -l < "$dir/reads_1.fq") / 4 )); n2=$(( $(wc -l < "$dir/reads_2.fq") / 4 ))
  if [ "$n1" != "$n2" ]; then echo "FAIL $d ${cov}x mates out of step: R1=$n1 R2=$n2"; bad=1
  elif [ "$n1" != "$kept" ]; then echo "FAIL $d ${cov}x wrote $n1 fragments, recorded $kept"; bad=1; fi
  if [ "$cov" = "$SOURCE_COVERAGE" ]; then
    if [ "$(content_sha "$dir/reads_1.fq")" != "$(gz_content_sha "$src/reads_1.fq.gz")" ] ||
       [ "$(content_sha "$dir/reads_2.fq")" != "$(gz_content_sha "$src/reads_2.fq.gz")" ]; then
      echo "FAIL $d ${cov}x is not byte-identical to the frozen source"; bad=1
    fi
  fi
  return $bad
}

# Item 7. Six run statuses exist -- COMPLETE, MODEL_INCOMPLETE, CANDIDATE_INCOMPLETE,
# SPELLING_INCOMPLETE, RESCORE_INCOMPLETE, BLOCK_CONSENSUS -- and only two of them are a result the
# reads are entitled to produce. A run that never wrote a K-arm status is NOT a no-call; it is missing,
# and its cause is usually named in the unsuffixed status file, which is read rather than ignored.
classify_run() { # <call prefix> <K> <failures file>
  local prefix="$1" k="$2" failures="$3" status
  if [ -f "$failures" ] && grep -q "	call	$(basename "$prefix")	K$k\$" "$failures" 2>/dev/null; then
    echo "FAILED"; return
  fi
  if [ ! -f "$prefix.K$k.mosaic_status.tsv" ]; then
    if [ -f "$prefix.mosaic_status.tsv" ]; then
      status="$(awk -F'\t' '$1=="status"{print $2; exit}' "$prefix.mosaic_status.tsv")"
      echo "MISSING_ARM:$status"
    else
      echo "MISSING"
    fi
    return
  fi
  status="$(awk -F'\t' '$1=="status"{print $2; exit}' "$prefix.K$k.mosaic_status.tsv")"
  case "$status" in
    COMPLETE) echo "COMPLETE" ;;
    RESCORE_INCOMPLETE|BLOCK_CONSENSUS) echo "NO_CALL:$status" ;;
    *) echo "INCOMPLETE:$status" ;;
  esac
}

run_note() { # <call prefix> -- the cause a degraded run recorded, if any
  local prefix="$1"
  [ -f "$prefix.mosaic_status.tsv" ] || return 0
  awk -F'\t' '$1=="note"{print $2; exit}' "$prefix.mosaic_status.tsv"
}

# Item 6. Declared completeness, evaluation completeness, the NA19240 acceptance values and the
# performance envelope, checked against what is on disk. Returns nonzero if anything is unmet.
completeness_verdict() { # <out dir> <settings> <donors...>
  local out="$1" settings="$2"; shift 2
  local donors=("$@")
  local bad=0 runs=0 arms=0 complete=0 nocall=0 incomplete=0 missing=0 failed=0
  local d cov k verdict prefix summary
  local prefixes; prefixes="$(setting call_prefixes | tr ',' ' ')"
  for d in "${donors[@]}"; do
    for cov in $COVERAGES; do
      prefix="$out/$d/call.$cov"
      if [ -f "$out/$d/resources.$cov.tsv" ]; then runs=$((runs + 1)); fi
      for k in $prefixes; do
        arms=$((arms + 1))
        verdict="$(classify_run "$prefix" "$k" "$out/failures.tsv")"
        case "$verdict" in
          COMPLETE)
            complete=$((complete + 1))
            summary="$out/$d/eval.$cov.K$k.summary.tsv"
            if [ ! -f "$summary" ]; then
              echo "FAIL $d ${cov}x K$k is COMPLETE but was never evaluated"; bad=1
            elif ! grep -q '^unordered_distance_total	' "$summary"; then
              echo "FAIL $d ${cov}x K$k is COMPLETE but its summary carries no accuracy block"; bad=1
            fi ;;
          NO_CALL:*)
            nocall=$((nocall + 1))
            [ -f "$out/$d/eval.$cov.K$k.summary.tsv" ] || {
              echo "FAIL $d ${cov}x K$k no-called (${verdict#NO_CALL:}) but was never evaluated"; bad=1; } ;;
          INCOMPLETE:*)
            incomplete=$((incomplete + 1))
            echo "FAIL $d ${cov}x K$k incomplete: ${verdict#INCOMPLETE:} $(run_note "$prefix")"; bad=1 ;;
          MISSING_ARM:*)
            missing=$((missing + 1))
            echo "FAIL $d ${cov}x K$k wrote no arm status; run says ${verdict#MISSING_ARM:}: $(run_note "$prefix")"
            bad=1 ;;
          MISSING)
            missing=$((missing + 1)); echo "FAIL $d ${cov}x K$k produced no status at all"; bad=1 ;;
          FAILED)
            failed=$((failed + 1)); echo "FAIL $d ${cov}x K$k the call itself failed"; bad=1 ;;
        esac
      done
    done
  done

  local want_runs want_arms
  want_runs="$(setting expected_shared_runs)"; want_arms="$(setting expected_k_arm_status_records)"
  [ "$runs" = "$want_runs" ] || { echo "FAIL $runs shared runs, $want_runs declared"; bad=1; }
  [ "$arms" = "$want_arms" ] || { echo "FAIL $arms arm records, $want_arms declared"; bad=1; }
  # -s, not -f: the run creates this file empty, and "zero recorded failures" is the good case.
  if [ -s "$out/failures.tsv" ]; then
    echo "FAIL $(wc -l < "$out/failures.tsv" | tr -d ' ') recorded failures:"; sed 's/^/       /' "$out/failures.tsv"
    bad=1
  fi

  # The acceptance values, at the source coverage, for the donor the settings name. Read from a
  # setting rather than hard-coded, so the same verdict serves the held-out cohort and can be gated.
  local adonor want_k0 want_k128 want_floor got want
  adonor="$(setting acceptance_donor)"
  want_k0="$(awk -F'\t' -v d="$adonor" '$1=="acceptance" && $2==d "_K0_distance"{print $3; exit}' "$settings")"
  want_k128="$(awk -F'\t' -v d="$adonor" '$1=="acceptance" && $2==d "_K128_distance"{print $3; exit}' "$settings")"
  want_floor="$(awk -F'\t' -v d="$adonor" '$1=="acceptance" && $2==d "_K128_attains_availability_floor"{print $3; exit}' "$settings")"
  for k in 0 128; do
    if [ "$k" = "0" ]; then want="$want_k0"; else want="$want_k128"; fi
    summary="$out/$adonor/eval.$SOURCE_COVERAGE.K$k.summary.tsv"
    got=""
    if [ -f "$summary" ]; then
      got="$(awk -F'\t' '$1=="unordered_distance_total"{print $2; exit}' "$summary")"
    fi
    if [ "$got" != "$want" ]; then
      echo "FAIL $adonor K$k distance '${got:-none}', acceptance $want"; bad=1
    fi
  done
  summary="$out/$adonor/eval.$SOURCE_COVERAGE.K128.summary.tsv"
  got=""
  if [ -f "$summary" ]; then
    got="$(awk -F'\t' '$1=="floor_attained_in_candidates"{print $2; exit}' "$summary")"
  fi
  if [ "$got" != "$want_floor" ]; then
    echo "FAIL $adonor K128 floor attained '${got:-none}', acceptance $want_floor"; bad=1
  fi

  # The performance envelope at the source coverage.
  local max_s max_rss secs rss res
  max_s="$(awk -F'\t' '$1=="acceptance" && $2=="max_30x_seconds"{print $3; exit}' "$settings")"
  max_rss="$(awk -F'\t' '$1=="acceptance" && $2=="max_rss_bytes"{print $3; exit}' "$settings")"
  for d in "${donors[@]}"; do
    res="$out/$d/resources.$SOURCE_COVERAGE.tsv"
    if [ ! -f "$res" ]; then continue; fi
    secs="$(awk -F'\t' '$1=="wall_seconds"{print $2; exit}' "$res")"
    rss="$(awk -F'\t' '$1=="peak_rss_bytes"{print $2; exit}' "$res")"
    if awk -v s="${secs:-0}" -v m="$max_s" 'BEGIN{exit !(s > m)}'; then
      echo "FAIL $d ${SOURCE_COVERAGE}x took ${secs}s, limit ${max_s}s"; bad=1
    fi
    if [ "${rss:-0}" -gt "$max_rss" ]; then
      echo "FAIL $d ${SOURCE_COVERAGE}x peak RSS $rss, limit $max_rss"; bad=1
    fi
  done

  printf 'complete\t%s\nlegitimate_no_call\t%s\nincomplete\t%s\nmissing\t%s\nfailed\t%s\n' \
    "$complete" "$nocall" "$incomplete" "$missing" "$failed" > "$out/accounting.tsv"
  echo "accounting: $complete complete, $nocall legitimate no-calls, $incomplete incomplete, $missing missing, $failed failed"
  return $bad
}

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

  # Reads: the LZO lineage, unsalted seed AND unsalted thinning.
  local salt thinning
  salt="$(setting seed_salt)"; thinning="$(setting thinning_selector)"
  if [ "$SEED" = "20260916" ] && [ "$salt" = "none" ] && [ "$thinning" = "unsalted" ] &&
     [ "$(setting thinning_name_salt)" = "none" ]; then
    echo "  reads bound to the ordinary-LZO stream: seed $SEED, salt $salt, thinning $thinning"
  else
    echo "  READ LINEAGE WRONG: seed=$SEED salt=$salt thinning=$thinning"; bad=$((bad + 1))
  fi

  # Frozen expectations must exist for every donor and coverage before any of them is produced.
  local missing=0 donors_n=0 d cov
  for d in $DONORS; do
    donors_n=$((donors_n + 1))
    [ -f "$FROZEN_ROOT/$d/source/reads_1.fq.gz" ] || missing=$((missing + 1))
    for field in reads1_sha256 reads2_sha256 library_sha256 fragments; do
      [ -n "$(expected_source "$SETTINGS" "$d" "$field")" ] || {
        echo "  NO FROZEN EXPECTATION: $d $field"; bad=$((bad + 1)); }
    done
    for cov in $COVERAGES; do
      [ -n "$(expected_subset "$SETTINGS" "$d" "$cov")" ] || {
        echo "  NO FROZEN SUBSET COUNT: $d ${cov}x"; bad=$((bad + 1)); }
    done
  done
  # Two different questions, and conflating them made this gate fail on every clean checkout: whether
  # the SETTINGS are consistent is a property of the manifest, while whether the READS are on this
  # machine is a property of the machine. The run path demands both; a settings check demands only the
  # first. A root that exists but is missing some donors is a real fault either way.
  if [ "$donors_n" != "$(setting expected_donors)" ]; then
    echo "  DONOR COUNT WRONG: $donors_n, expected $(setting expected_donors)"; bad=$((bad + 1))
  elif [ ! -d "$FROZEN_ROOT" ]; then
    if [ "${REQUIRE_FROZEN_SOURCES:-0}" = "1" ]; then
      echo "  FROZEN SOURCES ABSENT: $FROZEN_ROOT"; bad=$((bad + 1))
    else
      echo "  $donors_n donors declared; frozen sources not on this machine, settings checked alone"
    fi
  elif [ "$missing" != "0" ]; then
    echo "  FROZEN SOURCES MISSING: $missing of $donors_n under $FROZEN_ROOT"; bad=$((bad + 1))
  else
    echo "  $donors_n development donors, every frozen source present and expected"
  fi

  # Both source lineages bound by content.
  local lp lw lg
  for lp in "$(setting reads_lineage)" "$(setting caller_settings_lineage)"; do
    lw="$(expected_lineage "$SETTINGS" "$lp")"
    lg="$( [ -f "$REPO/$lp" ] && content_sha "$REPO/$lp" )"
    if [ -n "$lw" ] && [ "$lw" = "$lg" ]; then echo "  lineage bound: $lp"
    else echo "  LINEAGE MOVED OR UNBOUND: $lp ($lg vs $lw)"; bad=$((bad + 1)); fi
  done

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

  # The unresolved held-out design case must still be recorded as unresolved.
  if [ "$(awk -F'\t' '$1=="design_case" && $2=="HG00096" && $3=="resolution"{print $4; exit}' "$SETTINGS")" \
       = "UNRESOLVED_PENDING_REVIEW" ]; then
    echo "  HG00096 reference-donor case recorded as unresolved; held-out cohort stays closed"
  else
    echo "  HG00096 DESIGN CASE NOT RECORDED"; bad=$((bad + 1))
  fi

  local metrics
  metrics="$(awk -F'\t' '$1=="metric"' "$SETTINGS" | wc -l | tr -d ' ')"
  if [ "$metrics" -ge 30 ]; then echo "  $metrics metrics frozen, each naming its source key"
  else echo "  METRICS MISSING: only $metrics"; bad=$((bad + 1)); fi

  if [ "$bad" -ne 0 ]; then echo "manifest gate FAILED"; return 1; fi
  echo "manifest gate passed"
  return 0
}

if [ "${1:-}" = "--check-manifest" ]; then check_manifest; exit $?; fi
# Sourced by the gate: define the functions, run nothing.
if [ "${1:-}" = "--source-only" ]; then return 0 2>/dev/null || exit 0; fi

OUT="${1:?usage: BUILD=<build> tests/genotype_loo_c4_development.sh <out-dir>}"
BUILD="${BUILD:-$REPO/build}"
PROBE="$BUILD/genotype_probe"
[ -x "$PROBE" ] || { echo "genotype_probe not built at $PROBE" >&2; exit 1; }
# About to score reads, so their absence is fatal here even though it is not when checking settings.
REQUIRE_FROZEN_SOURCES=1
check_manifest
mkdir -p "$OUT"
[ -z "$(ls -A "$OUT" 2>/dev/null)" ] || { echo "$OUT is not empty" >&2; exit 1; }

gzip -dc "$REPO/$(setting graph)" > "$OUT/c4.gfa"

# Item 3: the frozen inputs, before anything is scored.
if ! check_frozen_inputs "$SETTINGS" "$FROZEN_ROOT" $DONORS; then
  echo "frozen inputs do not match the settings; refusing" >&2; exit 1
fi
echo "frozen inputs match the settings"

# The reference donor cannot be left out of its own coordinate frame: removing it leaves a graph with
# no reference path, and the run degrades to MODEL_INCOMPLETE instead of refusing. Caught here, by
# name, rather than discovered as a missing status file. See the HG00096 design case in the settings.
for donor in $DONORS; do
  if printf '%s' "$REFERENCE" | grep -q "^$donor#"; then
    echo "$donor is the reference donor; leave-one-out would remove the reference path" >&2
    exit 1
  fi
done

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
  printf 'probe_sha256\t%s\n' "$(content_sha "$PROBE")"
  printf 'script_sha256\t%s\n' "$(content_sha "${BASH_SOURCE[0]}")"
  printf 'settings_sha256\t%s\n' "$(content_sha "$SETTINGS")"
  printf 'graph_sha256\t%s\n' "$(content_sha "$OUT/c4.gfa")"
  printf 'donor_manifest_sha256\t%s\n' "$(content_sha "$DONOR_MANIFEST")"
  # Item 4: both source lineages, bound by content alongside everything else.
  printf 'reads_lineage_sha256\t%s\n' "$(content_sha "$REPO/$(setting reads_lineage)")"
  printf 'caller_settings_lineage_sha256\t%s\n' "$(content_sha "$REPO/$(setting caller_settings_lineage)")"
  printf 'reference\t%s\nseed\t%s\nseed_salt\t%s\nthinning_selector\t%s\ncoverages\t%s\n' \
    "$REFERENCE" "$SEED" "$(setting seed_salt)" "$(setting thinning_selector)" "$(echo $COVERAGES)"
  printf 'proposal_state_limit\t%s\ncall_prefixes\t%s\nphase_audit_candidates\t%s\n' \
    "$STATE_LIMIT" "$CALL_PREFIXES" "$PHASE_AUDIT"
  printf 'max_state_classes\t%s\nmax_physical_origins\t%s\nemission_cache_bytes\t%s\nthreads\t%s\n' \
    "$MAX_STATE_CLASSES" "$MAX_PHYSICAL_ORIGINS" "$EMISSION_CACHE_BYTES" "$THREADS"
  printf 'donors\t%s\n' "$(echo $DONORS)"
  for d in $DONORS; do
    printf 'reads_sha256_%s_1\t%s\n' "$d" "$(gz_content_sha "$FROZEN_ROOT/$d/source/reads_1.fq.gz")"
    printf 'reads_sha256_%s_2\t%s\n' "$d" "$(gz_content_sha "$FROZEN_ROOT/$d/source/reads_2.fq.gz")"
    printf 'library_sha256_%s\t%s\n' "$d" "$(content_sha "$FROZEN_ROOT/$d/source/library.txt")"
  done
} > "$OUT/signature.tsv"

# Item 4, continued: derived inputs are bound as they are produced -- each donor's LOO graph, and each
# subset FASTQ with its library. A subset that is regenerated or edited later cannot pass unnoticed.
printf 'kind\tdonor\tcoverage\tfile\tsha256\n' > "$OUT/inputs.tsv"
: > "$OUT/failures.tsv"

for donor in $DONORS; do
  dir="$OUT/$donor"; mkdir -p "$dir"
  src="$FROZEN_ROOT/$donor/source"
  # Leave-one-out: the donor's own two paths leave the MODEL graph. loo-gfa refuses unless exactly
  # two are removed, so a donor whose paths are named differently cannot slip through as an LZO run.
  if ! "$PROBE" loo-gfa --gfa "$OUT/c4.gfa" --donor "$donor" --out "$dir/c4.loo.gfa" \
       > "$dir/loo.log" 2>&1; then
    printf '%s\tloo_gfa\t%s\t-\n' "$donor" "$donor" >> "$OUT/failures.tsv"
    echo "LOO GRAPH FAILED: $donor" >&2
    continue
  fi
  printf 'loo_graph\t%s\t-\tc4.loo.gfa\t%s\n' "$donor" "$(content_sha "$dir/c4.loo.gfa")" >> "$OUT/inputs.tsv"

  for cov in $COVERAGES; do
    reads="$dir/cov$cov"; mkdir -p "$reads"
    # UNSALTED: no --name-salt. The frozen sources are unsalted, so the subsets must be too.
    if ! "$PROBE" subsample --reads1 "$src/reads_1.fq.gz" --reads2 "$src/reads_2.fq.gz" \
         --library "$src/library.txt" --source-coverage "$SOURCE_COVERAGE" --coverage "$cov" \
         --out "$reads" >> "$dir/loo.log" 2>&1; then
      printf '%s\tsubsample\t%s\t%s\n' "$donor" "$cov" "-" >> "$OUT/failures.tsv"
      echo "SUBSAMPLE FAILED: $donor ${cov}x" >&2
      continue
    fi
    if ! check_subset "$SETTINGS" "$donor" "$cov" "$reads" "$src"; then
      printf '%s\tsubset_gate\t%s\t-\n' "$donor" "$cov" >> "$OUT/failures.tsv"
      echo "SUBSET GATE FAILED: $donor ${cov}x" >&2
      continue
    fi
    for f in reads_1.fq reads_2.fq library.txt; do
      printf 'subset\t%s\t%s\t%s\t%s\n' "$donor" "$cov" "$f" "$(content_sha "$reads/$f")" >> "$OUT/inputs.tsv"
    done

    prefix="$dir/call.$cov"
    if ! /usr/bin/time "$TIME_FLAG" "$PROBE" genotype --gfa "$dir/c4.loo.gfa" --reference "$REFERENCE" \
      --reads1 "$reads/reads_1.fq" --reads2 "$reads/reads_2.fq" --library "$reads/library.txt" \
      --out-prefix "$prefix" --state-limit "$STATE_LIMIT" --prefixes "$CALL_PREFIXES" \
      --phase-audit "$PHASE_AUDIT" --threads "$THREADS" \
      --max-state-classes "$MAX_STATE_CLASSES" --max-physical-origins "$MAX_PHYSICAL_ORIGINS" \
      --emission-cache-bytes "$EMISSION_CACHE_BYTES" \
      > "$prefix.log" 2> "$prefix.time"; then
      # Recorded per arm so the classifier can tell a failed call from a missing one.
      for k in ${CALL_PREFIXES//,/ }; do
        printf '%s\tcall\t%s\tK%s\n' "$donor" "$(basename "$prefix")" "$k" >> "$OUT/failures.tsv"
      done
      echo "CALL FAILED: $donor ${cov}x (see $prefix.log)" >&2
    fi
    printf 'peak_rss_bytes\t%s\nwall_seconds\t%s\n' \
      "$(rss_of "$prefix.time")" "$(wall_of "$prefix.time")" > "$dir/resources.$cov.tsv"

    for k in ${CALL_PREFIXES//,/ }; do
      if [ -f "$prefix.K$k.mosaic_status.tsv" ]; then
        # model = the LOO graph the call was made against; truth = the FULL graph, which still
        # carries this donor's paths. That asymmetry is what makes it leave-one-out.
        if ! "$PROBE" evaluate --scheme loo --donor "$donor" --reference "$REFERENCE" \
          --run-prefix "$prefix.K$k" --model-gfa "$dir/c4.loo.gfa" --truth-gfa "$OUT/c4.gfa" \
          --out "$dir/eval.$cov.K$k" >> "$prefix.log" 2>&1; then
          printf '%s\tevaluate\t%s\tK%s\n' "$donor" "$cov" "$k" >> "$OUT/failures.tsv"
          echo "EVALUATE FAILED: $donor cov=$cov K=$k (see $prefix.log)" >&2
        fi
      fi
    done
  done
  echo "done $donor"
done

echo
# Item 5 and 6: the campaign's verdict is computed from what is on disk, and it is the exit status.
if ! completeness_verdict "$OUT" "$SETTINGS" $DONORS; then
  echo "campaign INCOMPLETE: $OUT" >&2
  exit 1
fi
echo "campaign complete: $OUT"
