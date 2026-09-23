#!/usr/bin/env bash
# Gates for the experimental greedy coordinate-repair search.
#
#   tests/genotype_repair_search_gate.sh <build dir> <scratch dir> [main run prefix]
#
# Six properties, in the order they matter:
#
#   1,2  the search REPRODUCES the campaign's own result and then goes past it -- round 1 must land on
#        the persisted R=32 ctx:2 winner, round 2 must reach the separately scored ctx:2 + sv:7
#        compound. These two read the main run's rounds table (third argument); they are skipped, not
#        faked, when it is absent, because re-running the 32-pair search inside a gate is not cheap.
#   3    order and homologue-swap invariance, reported per round by the search itself.
#   4    cycle prevention is real and deterministic: pre-occupying round 1's winner with --visited must
#        make the search pick something else, and two identical invocations must agree byte for byte.
#        Counting suppressed neighbours would not prove this -- a count of zero would pass silently.
#   5    exact-at-cap proceeds, cap-minus-one refuses before allocating.
#   6    construction never reads truth: the rounds table must be byte-identical with and without
#        --truth-gfa. If any selection consulted truth, dropping it would change the table.
#   7    resuming mid-search continues exactly as running straight through would.
#   8    an interrupted run can never be read as a completed one: a stop through the time-limit path
#        leaves only <out>.*.partial.tsv rows identical to the complete run's, a TIME_LIMIT status and
#        no summary; a cap refusal leaves RESOURCE_REFUSAL and nothing final; a complete run's
#        origin-cache demand reconciles with the work its rounds report.
#
# Gates 3-6 run at --local-pairs 1, where a round is small, since none of them depends on the
# neighbourhood being wide.
#
# tests/results is gitignored, so a clean checkout has nothing to run against. The gate says so.
set -uo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="${1:?usage: genotype_repair_search_gate.sh <build dir> <scratch dir> [main run prefix]}"
WORK="${2:?usage: genotype_repair_search_gate.sh <build dir> <scratch dir> [main run prefix]}"
MAIN="${3:-}"
PROBE="$BUILD/genotype_probe"
PKG="$REPO/tests/results/genotype_loo_c4_development_7e37bb4"
REF='HG00096#1#haplotype1-0000024:31848049-32074186'
rm -rf "$WORK"; mkdir -p "$WORK"

fail=0
same() { if [ "$2" = "$3" ]; then echo "  ok: $1 = $2"; else echo "  GATE FAILED: $1 is '$2', expected '$3'"; fail=1; fi; }
# An absent or empty value is not evidence of a difference: it is an absent measurement.
differs() {
  if [ -z "${2//\//}" ]; then echo "  GATE FAILED: $1 was not measured at all"; fail=1
  elif [ "$2" != "$3" ]; then echo "  ok: $1 ($2) differs from $3"
  else echo "  GATE FAILED: $1 is still '$2'"; fail=1; fi
}

[ -x "$PROBE" ] || { echo "genotype_probe not built at $PROBE" >&2; exit 1; }

D="$PKG/inputs/NA20509"
if [ ! -f "$D/c4.loo.gfa.gz" ] || [ ! -f "$PKG/NA20509/call.10.K128.mosaic_status.tsv" ]; then
  echo "skipped: the campaign artifacts under tests/results are absent (gitignored)."
  exit 0
fi

G="$WORK/c4.gfa"; L="$WORK/c4.loo.gfa"; R1="$WORK/r1.fq"; R2="$WORK/r2.fq"
gzip -dc "$PKG/inputs/c4.gfa.gz"  > "$G"  2>/dev/null || { echo "skipped: no full graph"; exit 0; }
gzip -dc "$D/c4.loo.gfa.gz"       > "$L"  2>/dev/null || { echo "skipped: no LOO graph"; exit 0; }
gzip -dc "$D/cov10/reads_1.fq.gz" > "$R1" 2>/dev/null || { echo "skipped: no reads"; exit 0; }
gzip -dc "$D/cov10/reads_2.fq.gz" > "$R2" 2>/dev/null || { echo "skipped: no reads"; exit 0; }

search() {   # search <out> <extra args...>
  local o="$1"; shift
  "$PROBE" repair-search --reference "$REF" --donor NA20509 \
    --model-gfa "$L" --run "$PKG/NA20509/call.10.K128" \
    --reads1 "$R1" --reads2 "$R2" --library "$D/cov10/library.txt" \
    --threads 8 --out "$o" "$@" > "$o.log" 2>&1
  echo "$?"
}
col() { awk -F'\t' -v r="$2" -v c="$3" 'NR==1{for(i=1;i<=NF;i++)h[$i]=i} $1==r{print $h[c]}' "$1"; }
key() { awk -F'\t' -v k="$2" '$1==k{print $2}' "$1"; }

echo "== 1. round 1 reproduces the persisted R=32 ctx:2 winner =="
if [ -n "$MAIN" ] && [ -f "$MAIN.rounds.tsv" ]; then
  same "round 1 block"  "$(col "$MAIN.rounds.tsv" 1 repaired_block)" "ctx:2"
  same "round 1 total"  "$(col "$MAIN.rounds.tsv" 1 total)" "-151931.605119"
  same "round 1 hap1"   "$(col "$MAIN.rounds.tsv" 1 hap1_alleles)" "15,10,87,12,3,6,58,11,51,1,99"
  same "round 1 hap2"   "$(col "$MAIN.rounds.tsv" 1 hap2_alleles)" "48,45,37,14,6,0,1,112,29,1,119"
  echo "== 2. round 2 produces and can select the scored ctx:2 + sv:7 compound =="
  # The compound is -151657.713214, the value score-set produced for it independently. Greedy ascent
  # does NOT select it: a ctx:6 repair outscores it. So the property to check is reachability, not
  # selection-by-default -- asserting the compound wins would assert that the search is worse than it
  # is. Two things establish it: round 2's certified separation equals winner minus compound, which
  # places the compound exactly second; and suppressing the winner makes round 2 select it outright.
  COMPOUND=-151657.713214
  won="$(col "$MAIN.rounds.tsv" 2 total)"
  if awk -v w="$won" -v c="$COMPOUND" 'BEGIN{exit !(w>c)}'; then
    echo "  ok: round 2's winner ($won) outscores the compound ($COMPOUND)"
  else
    echo "  GATE FAILED: round 2's winner $won does not outscore the compound"; fail=1
  fi
  gap="$(col "$MAIN.rounds.tsv" 2 separation)"
  if awk -v g="$gap" -v w="$won" -v c="$COMPOUND" 'BEGIN{exit !((g-(w-c))^2 < 1e-10)}'; then
    echo "  ok: the certified separation ($gap) is winner minus compound, so the compound ranks second"
  else
    echo "  GATE FAILED: separation $gap is not winner minus compound; the compound is not the runner-up"
    fail=1
  fi
  b1="$(col "$MAIN.rounds.tsv" 2 hap1_alleles)"
  b2="$(col "$MAIN.rounds.tsv" 2 hap2_alleles)"
  code="$(search "$WORK/compound" --local-pairs 32 --rounds 2 --visited "$b1|$b2")"
  same "suppressed re-run exit" "$code" "0"
  same "round 2 then repairs"   "$(col "$WORK/compound.rounds.tsv" 2 repaired_block)" "sv:7"
  same "round 2 then scores"    "$(col "$WORK/compound.rounds.tsv" 2 total)" "$COMPOUND"
  same "round 2 then moves"     "$(col "$WORK/compound.rounds.tsv" 2 moved)" "true"
  echo "== 3. invariance, every round of the main run =="
  same "order_invariant" "$(key "$MAIN.summary.tsv" order_invariant)" "true"
  same "swap_invariant"  "$(key "$MAIN.summary.tsv" swap_invariant)"  "true"
else
  echo "  skipped: no main run prefix given, so gates 1-3 have nothing to read."
fi

echo "== 4. cycle prevention is real, and the search is deterministic =="
base="$(search "$WORK/plain" --local-pairs 1 --rounds 1)"
same "plain exit" "$base" "0"
won1="$(col "$WORK/plain.rounds.tsv" 1 hap1_alleles)"
won2="$(col "$WORK/plain.rounds.tsv" 1 hap2_alleles)"
echo "  round 1 without suppression: $won1 / $won2"
code="$(search "$WORK/blocked" --local-pairs 1 --rounds 1 --visited "$won1|$won2")"
same "blocked exit" "$code" "0"
same "the pre-occupied candidate was suppressed" \
     "$(awk -F'\t' 'NR==2{print ($22>0)?"true":"false"}' "$WORK/blocked.rounds.tsv")" "true"
differs "round 1 winner under suppression" \
        "$(col "$WORK/blocked.rounds.tsv" 1 hap1_alleles)/$(col "$WORK/blocked.rounds.tsv" 1 hap2_alleles)" \
        "$won1/$won2"
code="$(search "$WORK/repeat" --local-pairs 1 --rounds 1)"
same "repeat exit" "$code" "0"
if diff -q <(cut -f1-28 "$WORK/plain.rounds.tsv") <(cut -f1-28 "$WORK/repeat.rounds.tsv") > /dev/null; then
  echo "  ok: two identical invocations agree exactly"
else
  echo "  GATE FAILED: the search is not deterministic"; fail=1
fi

echo "== 5. exact-at-cap proceeds, cap-minus-one refuses =="
cap="$(key "$WORK/plain.summary.tsv" neighbor_cap)"
echo "  computed neighbourhood bound: $cap"
code="$(search "$WORK/atcap" --local-pairs 1 --rounds 1 --max-neighbors "$cap")"
same "exactly at the cap" "$code" "0"
code="$(search "$WORK/under" --local-pairs 1 --rounds 1 --max-neighbors "$((cap - 1))")"
same "one below the cap" "$code" "3"
if [ -f "$WORK/under.rounds.tsv" ]; then
  echo "  GATE FAILED: a refused run still wrote a rounds table"; fail=1
else echo "  ok: the refusal wrote no rounds table"; fi

echo "== 6. construction never reads truth =="
code="$(search "$WORK/notruth" --local-pairs 1 --rounds 2)"
same "without truth" "$code" "0"
code="$(search "$WORK/withtruth" --local-pairs 1 --rounds 2 --truth-gfa "$G")"
same "with truth" "$code" "0"
same "truth was evaluated" "$(key "$WORK/withtruth.summary.tsv" truth_evaluated)" "true"
# Columns 29-30 are seconds and peak RSS, which differ between any two runs.
if diff -q <(cut -f1-28 "$WORK/notruth.rounds.tsv") <(cut -f1-28 "$WORK/withtruth.rounds.tsv") > /dev/null; then
  echo "  ok: the rounds table is identical with and without truth"
else
  echo "  GATE FAILED: the rounds table depends on truth being present"; fail=1
fi

echo "== 7. resuming mid-search continues exactly as running straight through would =="
# Resumption is only safe if it is indistinguishable from not having stopped. Run four rounds, then
# resume from round 2's winner with the first two rounds' candidates restored as visited, and require
# rounds 3 and 4 to come out identical. Comparing only the final score would not catch a resumed leg
# that took a different route to a similar place.
code="$(search "$WORK/straight" --local-pairs 1 --rounds 4)"
same "straight-through exit" "$code" "0"
seed1="$(key "$WORK/straight.summary.tsv" start_hap1)"
seed2="$(key "$WORK/straight.summary.tsv" start_hap2)"
r1a="$(col "$WORK/straight.rounds.tsv" 1 hap1_alleles)"; r1b="$(col "$WORK/straight.rounds.tsv" 1 hap2_alleles)"
r2a="$(col "$WORK/straight.rounds.tsv" 2 hap1_alleles)"; r2b="$(col "$WORK/straight.rounds.tsv" 2 hap2_alleles)"
r2t="$(col "$WORK/straight.rounds.tsv" 2 total)"
code="$(search "$WORK/resumed" --local-pairs 1 --rounds 2 --first-round 3           --start "$r2a|$r2b" --start-score "$r2t" --visited "$seed1|$seed2;$r1a|$r1b;$r2a|$r2b")"
same "resumed exit" "$code" "0"
same "the resumed leg declares it resumed" "$(key "$WORK/resumed.summary.tsv" resumed)" "true"
if diff -q <(awk -F'\t' 'NR==1 || $1>=3' "$WORK/straight.rounds.tsv" | cut -f1-28) \
           <(cut -f1-28 "$WORK/resumed.rounds.tsv") > /dev/null; then
  echo "  ok: rounds 3-4 are identical whether resumed or run straight through"
else
  echo "  GATE FAILED: the resumed leg diverges from the straight-through search"; fail=1
fi
# And the guard really is a guard: a start-score one nat off must refuse.
code="$(search "$WORK/badresume" --local-pairs 1 --rounds 1 --first-round 3 \
          --start "$r2a|$r2b" --start-score "$(awk -v t="$r2t" 'BEGIN{printf "%.6f", t+1}')")"
same "a wrong declared start score refuses" "$code" "4"


echo "== 8. an interrupted run can never be read as a completed one =="
exists() { [ -e "$1" ] && echo present || echo absent; }
DELTA=(--local-pairs 1 --delta 1 --origin-cache 1 --origin-cache-bytes 268435456)
code="$(search "$WORK/whole" "${DELTA[@]}" --rounds 3 --truth-gfa "$G")"
same "complete run exit" "$code" "0"
same "complete run status" "$(key "$WORK/whole.run_status.tsv" status)" "COMPLETE"
same "complete run leaves no partial rounds" "$(exists "$WORK/whole.rounds.partial.tsv")" "absent"
same "complete run demand reconciles" "$(key "$WORK/whole.run_status.tsv" reconciled)" "true"
code="$(search "$WORK/cut" "${DELTA[@]}" --rounds 3 --truth-gfa "$G" --interrupt-after-rounds 2)"
same "interrupted exit" "$code" "124"
same "interrupted status" "$(key "$WORK/cut.run_status.tsv" status)" "TIME_LIMIT"
same "interrupted is not complete" "$(key "$WORK/cut.run_status.tsv" complete)" "false"
same "interrupted rounds committed" "$(key "$WORK/cut.run_status.tsv" rounds_committed)" "2"
for f in rounds summary distance memory origin_demand; do
  same "interrupted run wrote no final $f table" "$(exists "$WORK/cut.$f.tsv")" "absent"
done
if diff -q <(head -3 "$WORK/whole.rounds.tsv" | cut -f1-28) <(cut -f1-28 "$WORK/cut.rounds.partial.tsv") > /dev/null; then
  echo "  ok: the two committed rounds are identical to the complete run's"
else
  echo "  GATE FAILED: committed rounds differ from the complete run's"; fail=1
fi
if diff -q <(head -4 "$WORK/whole.distance.tsv") "$WORK/cut.distance.partial.tsv" > /dev/null; then
  echo "  ok: the committed distances are identical to the complete run's"
else
  echo "  GATE FAILED: committed distances differ from the complete run's"; fail=1
fi
code="$(search "$WORK/capped" "${DELTA[@]}" --rounds 2 --max-physical-origins 1000)"
same "cap refusal exit" "$code" "5"
same "cap refusal status" "$(key "$WORK/capped.run_status.tsv" status)" "RESOURCE_REFUSAL"
same "cap refusal is not complete" "$(key "$WORK/capped.run_status.tsv" complete)" "false"
same "cap refusal wrote no rounds table" "$(exists "$WORK/capped.rounds.tsv")" "absent"
same "cap refusal wrote no summary" "$(exists "$WORK/capped.summary.tsv")" "absent"
# Every exit path removes the staged graphs and keeps the diagnostic record. The interrupted run exits
# through _Exit, which runs no destructor, so its cleanup is explicit; the others unwind.
code="$(search "$WORK/footprint" "${DELTA[@]}" --rounds 2 --footprint-limit-bytes 1)"
same "footprint refusal exit" "$code" "125"
same "footprint refusal status" "$(key "$WORK/footprint.run_status.tsv" status)" "RESOURCE_REFUSAL"
code="$(search "$WORK/broken" "${DELTA[@]}" --rounds 1 --start "1|2" --start-score 0)"
differs "exception exit" "$code" "0"
same "exception status" "$(key "$WORK/broken.run_status.tsv" status)" "FAILED"
for o in cut capped footprint broken; do
  same "$o: staged model graph removed" "$(exists "$WORK/$o.model.sorted.gfa")" "absent"
  same "$o: staged truth graph removed" "$(exists "$WORK/$o.truth.sorted.gfa")" "absent"
  same "$o: run status kept" "$(exists "$WORK/$o.run_status.tsv")" "present"
done
for f in rounds memory origin_demand distance round_truth; do
  same "interrupted run kept its partial $f table" "$(exists "$WORK/cut.$f.partial.tsv")" "present"
done

echo
if [ "$fail" -ne 0 ]; then echo "repair-search gate FAILED"; exit 1; fi
echo "repair-search gate passed"
