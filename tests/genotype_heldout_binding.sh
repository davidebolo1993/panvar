#!/usr/bin/env bash
# Gate: the held-out mosaic harness refuses to interpret results whose inputs or outputs have moved
# underneath them.
#
#   tests/genotype_heldout_binding.sh <scratch-dir>
#
# The signature machinery is the only thing standing between "these are held-out results" and "these are
# some outputs that happen to be in the directory". Round 14 drew a false conclusion from exactly that
# gap: a staged graph was regenerated while its FASTQs were reused, and every later call scored reads
# from one pair of haplotypes against a truth projected from another. Nothing announced it.
#
# So each mutation below is a way that could happen again, and each must REFUSE rather than report.
# Runs on stub files in seconds: it does not need the real C4 graph and does not run the ten-plan
# experiment.
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
HARNESS="$REPO/tests/genotype_c4_heldout_mosaic.sh"
WORK="${1:?usage: genotype_heldout_binding.sh <scratch-dir>}"
rm -rf "$WORK"; mkdir -p "$WORK"

# ---- a stub experiment -------------------------------------------------------------------------------
OUT="$WORK/out"
MANIFEST="$WORK/manifest.tsv"
DRAW="$WORK/draw.py"
ORIENTATION="$WORK/orientation.tsv"
PANEL_ALLELES="$WORK/panel_alleles.tsv"
GRAPH_GZ="$WORK/graph.gfa.gz"
CURATION="$OUT/mosaic_exclude.tsv"
PROBE="$WORK/fake_probe"
mkdir -p "$OUT"

printf 'plan\tname\tclass\tdonor_a\tdonor_b\tswap_hap_b\tswitches_hap1\tswitches_hap2\trequire_off_panel\n' > "$MANIFEST"
printf 'plan\tp1\tboth_homologues\tDONOR_A\tDONOR_B\t0\t3\t7\tboth\n' >> "$MANIFEST"
printf 'draw\n' > "$DRAW"
printf 'orientation\n' > "$ORIENTATION"
printf 'alleles\n' > "$PANEL_ALLELES"
printf 'graph\n' > "$GRAPH_GZ"
printf 'probe\n' > "$PROBE"
printf 'MOSAIC#1#synthetic:0-0\nMOSAIC#2#synthetic:0-0\n' > "$CURATION"

# Caller settings, as the harness holds them.
SOURCE_COVERAGE=30 SEED=20260921 READ_LENGTH=100 INSERT_MEAN=350 INSERT_SD=40 SIGMAS=4
SUBSTITUTION_ERROR=0.002 OUTLIER_WEIGHT=0.001 STATE_LIMIT=256 CALL_PREFIXES=0,128
PHASE_AUDIT=0 MIN_VARIANT_BP=50
COVERAGES="30 10 5 2"

plan_rows() { awk -F'\t' '$1=="plan" && $2!="name"' "$MANIFEST"; }

# The harness's own functions, sourced rather than reimplemented: a gate that reimplements the thing it
# guards proves only that the copy agrees with itself.
eval "$(awk '/^content_sha\(\) \{/,/^\}/' "$HARNESS")"
eval "$(awk '/^caller_settings_blob\(\) \{/,/^\}/' "$HARNESS")"
eval "$(awk '/^experiment_signature\(\) \{/,/^\}/' "$HARNESS")"
eval "$(awk '/^experiment_signature_id\(\) \{/' "$HARNESS")"
eval "$(awk '/^signature_incomplete\(\) \{/,/^\}/' "$HARNESS")"
eval "$(awk '/^call_signature\(\) \{/,/^\}/' "$HARNESS")"
eval "$(awk '/^require_call_signature\(\) \{/,/^\}/' "$HARNESS")"

# The post-subsampling gate, lifted out of the harness so the same program is exercised here.
awk '/<<.PY.$/{inside=1; body=""; next}
     inside && /^PY$/{ if (body ~ /NESTING BROKEN/) printf "%s", body; inside=0; next }
     inside{ body = body $0 "\n" }' "$HARNESS" > "$WORK/nested_gate.py"
[ -s "$WORK/nested_gate.py" ] || { echo "could not extract the nesting gate from the harness"; exit 1; }

plan_dir() {
  local d="$OUT/p1"
  mkdir -p "$d"
  printf 'staged graph\n' > "$d/staged.gfa"
  printf 'truth\n' > "$d/truth.tsv"
  local cov
  for cov in $COVERAGES; do
    mkdir -p "$d/cov$cov"
    printf 'reads one\n' > "$d/cov$cov/reads_1.fq"
    printf 'reads two\n' > "$d/cov$cov/reads_2.fq"
    printf 'library\n'   > "$d/cov$cov/library.txt"
  done
  printf '%s' "$d"
}
D="$(plan_dir)"

fail=0
check() {  # check <name> <expected: same|differ> <a> <b>
  local name="$1" want="$2" a="$3" b="$4" got
  if [ "$a" = "$b" ]; then got=same; else got=differ; fi
  if [ "$got" = "$want" ]; then printf '  ok    %s (%s)\n' "$name" "$got"
  else printf '  FAIL  %s: wanted %s, got %s\n' "$name" "$want" "$got"; fail=$((fail + 1)); fi
}
refuses() {  # refuses <name> <command...>
  local name="$1"; shift
  if "$@" > /dev/null 2>&1; then
    printf '  FAIL  %s: interpretation was ALLOWED\n' "$name"; fail=$((fail + 1))
  else
    printf '  ok    %s refused\n' "$name"
  fi
}
accepts() {
  local name="$1"; shift
  if "$@" > /dev/null 2>&1; then printf '  ok    %s accepted\n' "$name"
  else printf '  FAIL  %s: refused when it should have been accepted\n' "$name"; fail=$((fail + 1)); fi
}

echo "=== the experiment signature binds every input ==="
base="$(experiment_signature_id p1 30)"

printf 'reads one MODIFIED\n' > "$D/cov30/reads_1.fq"
check "stale reads" differ "$base" "$(experiment_signature_id p1 30)"
printf 'reads one\n' > "$D/cov30/reads_1.fq"

printf 'library MODIFIED\n' > "$D/cov30/library.txt"
check "changed library" differ "$base" "$(experiment_signature_id p1 30)"
printf 'library\n' > "$D/cov30/library.txt"

printf 'MOSAIC#1#synthetic:0-0\n' > "$CURATION"
check "changed curation file" differ "$base" "$(experiment_signature_id p1 30)"
printf 'MOSAIC#1#synthetic:0-0\nMOSAIC#2#synthetic:0-0\n' > "$CURATION"

old_width="$STATE_LIMIT"; STATE_LIMIT=1024
check "changed caller setting" differ "$base" "$(experiment_signature_id p1 30)"
STATE_LIMIT="$old_width"

printf 'probe REBUILT\n' > "$PROBE"
check "changed binary" differ "$base" "$(experiment_signature_id p1 30)"
printf 'probe\n' > "$PROBE"

printf 'plan\tp1\tboth_homologues\tDONOR_A\tDONOR_B\t0\t4\t7\tboth\n' > "$WORK/row"
{ head -1 "$MANIFEST"; cat "$WORK/row"; } > "$MANIFEST.new" && mv "$MANIFEST.new" "$MANIFEST"
check "changed plan row" differ "$base" "$(experiment_signature_id p1 30)"
{ head -1 "$MANIFEST"; printf 'plan\tp1\tboth_homologues\tDONOR_A\tDONOR_B\t0\t3\t7\tboth\n'; } > "$MANIFEST.new"
mv "$MANIFEST.new" "$MANIFEST"
check "plan row restored" same "$base" "$(experiment_signature_id p1 30)"

mv "$D/cov30/reads_2.fq" "$D/cov30/hidden"
if signature_incomplete p1 30; then printf '  ok    missing input refuses (MISSING component)\n'
else printf '  FAIL  missing input produced a complete signature\n'; fail=$((fail + 1)); fi
mv "$D/cov30/hidden" "$D/cov30/reads_2.fq"

echo "=== the call signature is verified before anything is interpreted ==="
for k in ${CALL_PREFIXES//,/ }; do
  printf 'status\n'     > "$D/call.30.K$k.mosaic_status.tsv"
  printf 'scores\n'     > "$D/call.30.K$k.mosaic_scores.tsv"
  printf 'candidates\n' > "$D/call.30.K$k.mosaic_candidates.tsv"
  printf 'blocks\n'     > "$D/call.30.K$k.block_genotypes.tsv"
done
call_signature p1 30 > "$D/call.30.signature"
accepts "unchanged run" require_call_signature p1 30

printf 'status TAMPERED\n' > "$D/call.30.K0.mosaic_status.tsv"
refuses "stale call output" require_call_signature p1 30
printf 'status\n' > "$D/call.30.K0.mosaic_status.tsv"
accepts "call output restored" require_call_signature p1 30

printf 'reads one MODIFIED\n' > "$D/cov30/reads_1.fq"
refuses "reads changed after calling" require_call_signature p1 30
printf 'reads one\n' > "$D/cov30/reads_1.fq"

STATE_LIMIT=1024
refuses "caller setting changed after calling" require_call_signature p1 30
STATE_LIMIT="$old_width"

printf 'MOSAIC#1#synthetic:0-0\n' > "$CURATION"
refuses "curation changed after calling" require_call_signature p1 30
printf 'MOSAIC#1#synthetic:0-0\nMOSAIC#2#synthetic:0-0\n' > "$CURATION"

printf 'probe REBUILT\n' > "$PROBE"
refuses "binary changed after calling" require_call_signature p1 30
printf 'probe\n' > "$PROBE"

mv "$D/call.30.K128.block_genotypes.tsv" "$D/call.30.K128.hidden"
refuses "missing arm output" require_call_signature p1 30
mv "$D/call.30.K128.hidden" "$D/call.30.K128.block_genotypes.tsv"

rm -f "$D/call.30.signature"
refuses "no call signature at all" require_call_signature p1 30

echo "=== the post-subsampling gate catches a broken nested subset ==="
nest="$WORK/nest"
write_reads() {  # write_reads <dir> <count>
  local dir="$1" n="$2" i
  mkdir -p "$dir"
  : > "$dir/reads_1.fq"; : > "$dir/reads_2.fq"
  for i in $(seq 1 "$n"); do
    printf '@frag_%d/1\nACGT\n+\nIIII\n' "$i" >> "$dir/reads_1.fq"
    printf '@frag_%d/2\nACGT\n+\nIIII\n' "$i" >> "$dir/reads_2.fq"
  done
}
# Properly nested: 2x subset of 5x subset of 10x subset of 30x.
write_reads "$nest/cov30" 40
write_reads "$nest/cov10" 20
write_reads "$nest/cov5" 10
write_reads "$nest/cov2" 4
accepts "properly nested subsets" python3 "$WORK/nested_gate.py" "$nest" p1 "$COVERAGES"

# A fragment present at 5x but absent at 10x: the subsets are no longer nested, so a coverage
# comparison would be comparing two different draws.
printf '@frag_99/1\nACGT\n+\nIIII\n' >> "$nest/cov5/reads_1.fq"
printf '@frag_99/2\nACGT\n+\nIIII\n' >> "$nest/cov5/reads_2.fq"
refuses "broken nested subset" python3 "$WORK/nested_gate.py" "$nest" p1 "$COVERAGES"
write_reads "$nest/cov5" 10

# Mates out of step: reads_2 short by one record. Counted rather than `head -n -4`, which is a GNU
# extension and truncates nothing on the BSD head this runs against.
records="$(wc -l < "$nest/cov10/reads_2.fq")"
head -n "$((records - 4))" "$nest/cov10/reads_2.fq" > "$nest/cov10/trimmed"
mv "$nest/cov10/trimmed" "$nest/cov10/reads_2.fq"
refuses "mates out of step" python3 "$WORK/nested_gate.py" "$nest" p1 "$COVERAGES"
write_reads "$nest/cov10" 20

# Mate names disagreeing at one record, with the counts still equal.
python3 - "$nest/cov10/reads_2.fq" <<'PY'
import sys
path = sys.argv[1]
lines = open(path).read().splitlines()
lines[4] = "@frag_OTHER/2"
open(path, "w").write("\n".join(lines) + "\n")
PY
refuses "mate names disagree" python3 "$WORK/nested_gate.py" "$nest" p1 "$COVERAGES"
write_reads "$nest/cov10" 20
accepts "nesting restored" python3 "$WORK/nested_gate.py" "$nest" p1 "$COVERAGES"

echo "=== the completeness gate follows the declared output contract ==="
eval "$(awk '/^arm_output_problems\(\) \{/,/^\}/' "$HARNESS")"
arm="$WORK/arm"; mkdir -p "$arm"

mk_arm() {  # mk_arm <status|-> [output suffixes...]
  rm -f "$arm"/a.*
  if [ "$1" != "-" ]; then printf 'key\tvalue\nstatus\t%s\n' "$1" > "$arm/a.mosaic_status.tsv"; fi
  shift
  local f
  for f in "$@"; do printf 'x\n' > "$arm/a.$f"; done
}
check_arm() {  # check_arm <name> <ok|fail>
  local name="$1" want="$2" n
  n="$(arm_output_problems "$arm/a" | grep -c . || true)"
  if { [ "$want" = ok ] && [ "$n" -eq 0 ]; } || { [ "$want" = fail ] && [ "$n" -gt 0 ]; }; then
    printf '  ok    %s (%s problem(s))\n' "$name" "$n"
  else
    printf '  FAIL  %s: wanted %s, got %s problem(s)\n' "$name" "$want" "$n"
    arm_output_problems "$arm/a" | sed 's/^/          /'
    fail=$((fail + 1))
  fi
}

COMPLETE_SET="mosaic_scores.tsv mosaic_candidates.tsv block_genotypes.tsv genotypes.tsv mosaic_haplotypes.fa"

mk_arm COMPLETE $COMPLETE_SET
check_arm "COMPLETE with every output" ok

# 1. a complete arm missing a block output must fail.
mk_arm COMPLETE mosaic_scores.tsv mosaic_candidates.tsv genotypes.tsv mosaic_haplotypes.fa
check_arm "COMPLETE missing block_genotypes" fail

# 2. an incomplete arm without block output must PASS: status only is the contract, and reading that
#    absence as a missing output is what stopped an experiment in which nothing was missing.
mk_arm RESCORE_INCOMPLETE mosaic_scores.tsv mosaic_candidates.tsv
check_arm "RESCORE_INCOMPLETE with status and diagnostics only" ok
mk_arm CANDIDATE_INCOMPLETE
check_arm "CANDIDATE_INCOMPLETE with status alone" ok

# 3. an incomplete arm CARRYING a block output must fail: a refused run emitted calls.
mk_arm RESCORE_INCOMPLETE mosaic_scores.tsv block_genotypes.tsv
check_arm "RESCORE_INCOMPLETE carrying block_genotypes" fail
mk_arm MODEL_INCOMPLETE genotypes.tsv
check_arm "MODEL_INCOMPLETE carrying genotypes" fail

# 4. block consensus carrying a phased FASTA must fail: phase is unresolved under consensus.
mk_arm BLOCK_CONSENSUS mosaic_scores.tsv mosaic_candidates.tsv block_genotypes.tsv genotypes.tsv
check_arm "BLOCK_CONSENSUS without FASTA" ok
mk_arm BLOCK_CONSENSUS $COMPLETE_SET
check_arm "BLOCK_CONSENSUS carrying phased FASTA" fail

# 5. a missing status must fail, whatever else is present.
mk_arm - $COMPLETE_SET
check_arm "missing status file" fail
printf 'key\tvalue\nnote\tno status key\n' > "$arm/a.mosaic_status.tsv"
check_arm "status file declaring no status" fail

if [ "$fail" -ne 0 ]; then echo "genotype_heldout_binding: $fail gate(s) FAILED"; exit 1; fi
echo "genotype_heldout_binding: all gates passed"
