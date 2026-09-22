#!/usr/bin/env bash
# Gate: the C4 development LOO campaign refuses inputs and results that have moved underneath it.
#
#   tests/genotype_loo_c4_development_gates.sh <scratch-dir>
#
# The campaign's first run reported success while 22 of 38 accuracy blocks were missing, because a
# swallowed evaluate failure looked exactly like a donor with nothing to say. It also thinned unsalted
# frozen sources with a donor-salted selector, scoring a different subset than the one it claimed to
# pair against. Both are silent failures, so each mutation below must REFUSE rather than report.
#
# The harness's own functions are sourced rather than reimplemented: a gate that reimplements the
# thing it guards proves only that the copy agrees with itself.
#
# Runs on stub files in about a second. It does not need the C4 graph and does not call the genotyper.
set -uo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
HARNESS="$REPO/tests/genotype_loo_c4_development.sh"
WORK="${1:?usage: genotype_loo_c4_development_gates.sh <scratch-dir>}"
rm -rf "$WORK"; mkdir -p "$WORK"

# shellcheck source=/dev/null
source "$HARNESS" --source-only
set +e  # the harness brings `set -e`; this gate deliberately runs things that must fail

# The stub world below overwrites the harness's globals, and check_manifest reads fourteen of them.
# The real values are captured here and restored wholesale before the committed manifest is checked:
# restoring only the one that happened to break is how a gate starts passing for the wrong reason.
REAL_SETTINGS="$SETTINGS"; REAL_COVERAGES="$COVERAGES"; REAL_SOURCE_COVERAGE="$SOURCE_COVERAGE"
REAL_STATE_LIMIT="$STATE_LIMIT"; REAL_PHASE_AUDIT="$PHASE_AUDIT"; REAL_CALL_PREFIXES="$CALL_PREFIXES"
restore_real() {
  SETTINGS="$REAL_SETTINGS"; COVERAGES="$REAL_COVERAGES"; SOURCE_COVERAGE="$REAL_SOURCE_COVERAGE"
  STATE_LIMIT="$REAL_STATE_LIMIT"; PHASE_AUDIT="$REAL_PHASE_AUDIT"; CALL_PREFIXES="$REAL_CALL_PREFIXES"
}

fail=0
refused() { # <label> <status>
  if [ "$2" -eq 0 ]; then echo "  GATE FAILED: $1 was ACCEPTED"; fail=1; else echo "  refused: $1"; fi
}
accepted() { # <label> <status>
  if [ "$2" -eq 0 ]; then echo "  accepted: $1"; else echo "  GATE FAILED: $1 was REFUSED"; fail=1; fi
}
same() { # <label> <got> <want>
  if [ "$2" = "$3" ]; then echo "  $1: $2"; else echo "  GATE FAILED: $1 is '$2', expected '$3'"; fail=1; fi
}

# ---- a stub world ------------------------------------------------------------------------------
D=DSTUB
FROZEN="$WORK/frozen/$D/source"; mkdir -p "$FROZEN"
OUT="$WORK/out"; mkdir -p "$OUT/$D/cov30"
STUB="$WORK/settings.tsv"

make_reads() { # <dir> <n> [name prefix]
  local dir="$1" n="$2" p="${3:-frag}" i
  : > "$dir/reads_1.fq"; : > "$dir/reads_2.fq"
  for ((i = 1; i <= n; i++)); do
    printf '@%s_%d/1\nACGT\n+\nIIII\n' "$p" "$i" >> "$dir/reads_1.fq"
    printf '@%s_%d/2\nTGCA\n+\nIIII\n' "$p" "$i" >> "$dir/reads_2.fq"
  done
}
make_reads "$FROZEN" 10
printf 'insert_mean_bp 350\nfragment_start_intensity 0.075\n' > "$FROZEN/library.txt"
gzip -kf "$FROZEN/reads_1.fq" "$FROZEN/reads_2.fq"

R1="$(gz_content_sha "$FROZEN/reads_1.fq.gz")"
R2="$(gz_content_sha "$FROZEN/reads_2.fq.gz")"
LIB="$(content_sha "$FROZEN/library.txt")"

{
  printf 'key\tvalue\n'
  printf 'acceptance_donor\t%s\n' "$D"
  printf 'call_prefixes\t0,128\n'
  printf 'expected_shared_runs\t1\nexpected_k_arm_status_records\t2\n'
  printf 'source\t%s\treads1_sha256\t%s\n' "$D" "$R1"
  printf 'source\t%s\treads2_sha256\t%s\n' "$D" "$R2"
  printf 'source\t%s\tlibrary_sha256\t%s\n' "$D" "$LIB"
  printf 'source\t%s\tfragments\t10\n' "$D"
  printf 'subset\t%s\t30\t10\n' "$D"
  printf 'acceptance\t%s_K0_distance\t6\n' "$D"
  printf 'acceptance\t%s_K128_distance\t4\n' "$D"
  printf 'acceptance\t%s_K128_attains_availability_floor\ttrue\n' "$D"
  printf 'acceptance\tmax_30x_seconds\t360\n'
  printf 'acceptance\tmax_rss_bytes\t2147483648\n'
} > "$STUB"

# Point the sourced functions at the stub world.
SETTINGS="$STUB"; COVERAGES="30"; SOURCE_COVERAGE=30

subset_dir() { # rebuild a clean, correct 30x subset
  local dir="$OUT/$D/cov30"
  rm -rf "$dir"; mkdir -p "$dir"
  make_reads "$dir" 10
  cp "$FROZEN/library.txt" "$dir/library.txt"
  printf 'key\tvalue\nkept_fragments\t10\nname_salt\tnone\n' > "$dir/subsample.tsv"
}
arm() { # <K> <status> -- write a K-arm status file
  printf 'key\tvalue\nstatus\t%s\n' "$2" > "$OUT/$D/call.30.K$1.mosaic_status.tsv"
}
summary() { # <K> <distance> [floor_attained]
  { printf 'key\tvalue\n'
    [ -n "${2:-}" ] && printf 'unordered_distance_total\t%s\n' "$2"
    printf 'floor_attained_in_candidates\t%s\n' "${3:-true}"
  } > "$OUT/$D/eval.30.K$1.summary.tsv"
}
healthy() { # the state a good run leaves behind
  rm -f "$OUT/$D"/call.30.K*.mosaic_status.tsv "$OUT/$D"/eval.30.K*.summary.tsv \
        "$OUT/$D/call.30.mosaic_status.tsv"
  : > "$OUT/failures.tsv"
  printf 'wall_seconds\t100\npeak_rss_bytes\t1000000\n' > "$OUT/$D/resources.30.tsv"
  arm 0 COMPLETE; arm 128 COMPLETE
  summary 0 6 true; summary 128 4 true
  subset_dir
}

echo "== positive control: a healthy campaign is accepted =="
healthy
check_frozen_inputs "$STUB" "$WORK/frozen" "$D" > /dev/null; accepted "frozen inputs" $?
check_subset "$STUB" "$D" 30 "$OUT/$D/cov30" "$FROZEN" > /dev/null; accepted "30x subset" $?
completeness_verdict "$OUT" "$STUB" "$D" > "$WORK/verdict.log" 2>&1; accepted "completeness verdict" $?

echo "== mutation 1: changed read content =="
healthy
printf '@frag_1/1\nTTTT\n+\nIIII\n' > "$FROZEN/reads_1.fq"; gzip -kf "$FROZEN/reads_1.fq"
check_frozen_inputs "$STUB" "$WORK/frozen" "$D" > /dev/null; refused "a changed frozen read" $?
make_reads "$FROZEN" 10; gzip -kf "$FROZEN/reads_1.fq" "$FROZEN/reads_2.fq"

echo "== mutation 2: salted thinning =="
healthy
printf 'key\tvalue\nkept_fragments\t10\nname_salt\t%s\n' "$D" > "$OUT/$D/cov30/subsample.tsv"
check_subset "$STUB" "$D" 30 "$OUT/$D/cov30" "$FROZEN" > /dev/null; refused "a donor-salted thinning" $?

echo "== mutation 3: a missing mate =="
healthy
head -n 36 "$OUT/$D/cov30/reads_2.fq" > "$WORK/short.fq"; mv "$WORK/short.fq" "$OUT/$D/cov30/reads_2.fq"
check_subset "$STUB" "$D" 30 "$OUT/$D/cov30" "$FROZEN" > /dev/null; refused "mates out of step" $?

echo "== mutation 3b: a subset that is not the frozen one =="
healthy
make_reads "$OUT/$D/cov30" 9
printf 'key\tvalue\nkept_fragments\t9\nname_salt\tnone\n' > "$OUT/$D/cov30/subsample.tsv"
check_subset "$STUB" "$D" 30 "$OUT/$D/cov30" "$FROZEN" > /dev/null; refused "a retained count off by one" $?

echo "== mutation 4: a failed call =="
healthy
printf '%s\tcall\tcall.30\tK0\n' "$D" > "$OUT/failures.tsv"
same "classify_run on a failed call" "$(classify_run "$OUT/$D/call.30" 0 "$OUT/failures.tsv")" "FAILED"
completeness_verdict "$OUT" "$STUB" "$D" > /dev/null 2>&1; refused "a failed call" $?

echo "== mutation 5: a failed evaluation =="
healthy
rm -f "$OUT/$D/eval.30.K0.summary.tsv"
completeness_verdict "$OUT" "$STUB" "$D" > /dev/null 2>&1; refused "a COMPLETE run that was never evaluated" $?
healthy
summary 0 "" true   # a summary with no accuracy block
completeness_verdict "$OUT" "$STUB" "$D" > /dev/null 2>&1; refused "an evaluation with no accuracy block" $?

echo "== mutation 6: a missing status =="
healthy
rm -f "$OUT/$D/call.30.K0.mosaic_status.tsv"
same "classify_run with no status" "$(classify_run "$OUT/$D/call.30" 0 "$OUT/failures.tsv")" "MISSING"
completeness_verdict "$OUT" "$STUB" "$D" > /dev/null 2>&1; refused "a missing arm status" $?
# The HG00096 shape: no arm status, but the unsuffixed status names the cause.
printf 'key\tvalue\nstatus\tMODEL_INCOMPLETE\nnote\tReference path not found\n' \
  > "$OUT/$D/call.30.mosaic_status.tsv"
same "classify_run with only an unsuffixed status" \
  "$(classify_run "$OUT/$D/call.30" 0 "$OUT/failures.tsv")" "MISSING_ARM:MODEL_INCOMPLETE"
same "the cause is reported, not lost" "$(run_note "$OUT/$D/call.30")" "Reference path not found"

echo "== item 7: a legitimate no-call is not a failure =="
healthy
arm 128 RESCORE_INCOMPLETE
same "classify_run on an equivalence set" \
  "$(classify_run "$OUT/$D/call.30" 128 "$OUT/failures.tsv")" "NO_CALL:RESCORE_INCOMPLETE"
arm 128 BLOCK_CONSENSUS
same "classify_run on a block-consensus run" \
  "$(classify_run "$OUT/$D/call.30" 128 "$OUT/failures.tsv")" "NO_CALL:BLOCK_CONSENSUS"
# ... while a machinery failure with the same "not COMPLETE" shape is NOT a no-call.
arm 128 SPELLING_INCOMPLETE
same "classify_run on an incomplete run" \
  "$(classify_run "$OUT/$D/call.30" 128 "$OUT/failures.tsv")" "INCOMPLETE:SPELLING_INCOMPLETE"

echo "== the manifest gate still rejects the round-7 combination =="
restore_real
check_manifest > "$WORK/manifest.log" 2>&1; accepted "the committed manifest" $?
STUB7="$WORK/round7.tsv"
sed -e 's/^proposal_state_limit\t256/proposal_state_limit\t1024/' \
    -e 's/^max_mosaic_candidates\t128/max_mosaic_candidates\t8/' \
    -e 's/^phase_audit_candidates\t0/phase_audit_candidates\t256/' \
    "$SETTINGS" > "$STUB7"
SETTINGS="$STUB7"; STATE_LIMIT=1024; PHASE_AUDIT=256
check_manifest > "$WORK/round7.log" 2>&1; refused "the round-7 combination" $?
grep -q "ROUND-7 SETTINGS DETECTED" "$WORK/round7.log"
same "round-7 is refused by name" "$?" "0"

echo "== the frozen-source check separates a clean checkout from a broken one =="
restore_real
SAVED_FROZEN="$FROZEN_ROOT"
# tests/results is gitignored, so a clean checkout legitimately has no reads. An earlier revision of
# this gate conflated that with a fault and failed on every fresh checkout and in CI.
FROZEN_ROOT="$WORK/no-such-root"
REQUIRE_FROZEN_SOURCES=0 check_manifest > "$WORK/absent.log" 2>&1
accepted "settings alone, with no reads on the machine" $?
# The run path is about to score reads, so there the same absence is fatal.
REQUIRE_FROZEN_SOURCES=1 check_manifest > "$WORK/absent_strict.log" 2>&1
refused "a campaign run with no reads on the machine" $?
# A root that exists but is missing a donor is a fault in either mode: that is a broken checkout, not
# an absent one, and it is how a campaign silently scores seven donors while claiming eight.
PARTIAL="$WORK/partial"; mkdir -p "$PARTIAL"
skipped=""
for d in $DONORS; do
  if [ -z "$skipped" ]; then skipped="$d"; continue; fi
  mkdir -p "$PARTIAL/$d/source"; : > "$PARTIAL/$d/source/reads_1.fq.gz"
done
FROZEN_ROOT="$PARTIAL"
REQUIRE_FROZEN_SOURCES=0 check_manifest > "$WORK/partial.log" 2>&1
refused "a frozen root missing one donor" $?
grep -q "FROZEN SOURCES MISSING: 1 of 8" "$WORK/partial.log"
same "the missing donor is counted, not just noticed" "$?" "0"
FROZEN_ROOT="$SAVED_FROZEN"

echo
if [ "$fail" -ne 0 ]; then echo "LOO development gates FAILED"; exit 1; fi
echo "LOO development gates passed"
