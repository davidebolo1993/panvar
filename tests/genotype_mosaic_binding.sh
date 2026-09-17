#!/usr/bin/env bash
# Gate: the mosaic experiment harness refuses to interpret outputs whose inputs have moved underneath
# them.  Round 14 drew a false conclusion because a staged graph was regenerated while its FASTQs were
# reused, so every later call scored reads from one pair of haplotypes against a truth projected from
# another.  These gates exercise the signature machinery that now prevents that -- cheaply, on stub
# files, without running the five-plan experiment, which takes minutes and needs the real C4 graph.
#
#   tests/genotype_mosaic_binding.sh <scratch-dir>
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
HARNESS="$REPO/tests/genotype_c4_mosaic_experiment.sh"
WORK="${1:?usage: genotype_mosaic_binding.sh <scratch-dir>}"
rm -rf "$WORK"; mkdir -p "$WORK"

# The harness's own signature function, sourced rather than reimplemented: a gate that reimplements the
# thing it guards proves only that the copy agrees with itself.
PROBE="$WORK/fake_probe"; printf 'stub\n' > "$PROBE"
SEED=20260916 DEPTH=30 READ_LENGTH=100 INSERT_MEAN=350 INSERT_SD=40 SIGMAS=4
SUBSTITUTION_ERROR=0.002 OUTLIER_WEIGHT=0.001
eval "$(awk '/^input_bundle_signature\(\) \{/,/^\}/' "$HARNESS")"

plan() {
  local d="$WORK/$1"; mkdir -p "$d/source"
  printf 'graph\n'   > "$d/staged.gfa"
  printf 'reads one\n' > "$d/source/reads_1.fq"
  printf 'reads two\n' > "$d/source/reads_2.fq"
  printf 'library\n'  > "$d/source/library.txt"
  echo "$d"
}
fail=0
check() {  # check <name> <expected: same|differ> <sig-a> <sig-b>
  local name="$1" want="$2" a="$3" b="$4" got
  if [ "$a" = "$b" ]; then got=same; else got=differ; fi
  if [ "$got" = "$want" ]; then
    printf '  ok    %s (%s)\n' "$name" "$got"
  else
    printf '  FAIL  %s: wanted %s, got %s\n' "$name" "$want" "$got"; fail=$((fail+1))
  fi
}

d="$(plan one)"
base="$(input_bundle_signature "$d" one)"

# 1. a modified reads_2.fq must be detected -- the file the old read-only signature never covered.
printf 'reads two modified\n' > "$d/source/reads_2.fq"
check "modified reads_2.fq" differ "$base" "$(input_bundle_signature "$d" one)"
printf 'reads two\n' > "$d/source/reads_2.fq"

# 2. a modified library.txt must be detected: the declared insert model is part of what was scored.
printf 'library modified\n' > "$d/source/library.txt"
check "modified library.txt" differ "$base" "$(input_bundle_signature "$d" one)"
printf 'library\n' > "$d/source/library.txt"

# 3. a missing input is a refusal, never compatibility.
mv "$d/source/reads_2.fq" "$d/source/hidden"
sig_missing="$(input_bundle_signature "$d" one)"
case "$sig_missing" in
  MISSING:*) printf '  ok    missing input refuses (%s)\n' "$sig_missing" ;;
  *) printf '  FAIL  missing input produced a usable signature: %s\n' "$sig_missing"; fail=$((fail+1)) ;;
esac
mv "$d/source/hidden" "$d/source/reads_2.fq"

# 4. a changed simulation parameter with an UNCHANGED graph must still change the signature, or a
#    reparameterised draw would silently reuse the previous reads.
old_depth="$DEPTH"; DEPTH=10
check "changed simulation parameter" differ "$base" "$(input_bundle_signature "$d" one)"
DEPTH="$old_depth"

# 5. the plan salt is part of the bundle: two plans sharing a graph must not share a signature.
check "plan salt" differ "$base" "$(input_bundle_signature "$d" two)"

# 6. identical content in a different directory stays valid: the signature must be content-addressed,
#    not path-addressed, or copying a tree looks like tampering.
cp -R "$d" "$WORK/copied"
check "copied directory" same "$base" "$(input_bundle_signature "$WORK/copied" one)"

# 7. restoring the original content restores the original signature.
check "restored content" same "$base" "$(input_bundle_signature "$d" one)"

[ "$fail" -eq 0 ] || { echo "binding gates failed: $fail"; exit 1; }
echo "mosaic input binding: all gates pass"
