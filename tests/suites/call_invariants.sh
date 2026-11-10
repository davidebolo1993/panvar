#!/usr/bin/env bash
# call_invariants.sh - properties `panvar call` must hold whatever the input's incidental form.
#
# call_stats.sh pins coordinates on forward-chain fixtures. The real loci showed that several things
# the output must NOT depend on are not tested anywhere, and most of them never fire on the six
# committed loci (0 of 55 bubbles run against the reference there), so only a synthetic fixture
# reaches them:
#
#   1. bubble orientation   the same bubble given as (source,sink) and as (sink,source)
#   2. P-line order         the same graph with its paths reversed or shuffled
#   3. node numbering       the same graph renumbered, its VCF mapped back
#   4. chaining             A~B~C merge into one record although A and C do not match directly
#   5. substitution         a substitution and a pure deletion of the same block
#   6. short UNIT            a 9 bp self-loop walked ten times gives one UNIT DUP, with and without --cn
#
# A case whose fix is a later plan item is an expected failure, tagged with that item's id: it
# prints 'xfail' while it still fails and 'XPASS' once it passes, and neither fails the run. When
# an XPASS appears, the fix has landed and the case should become a hard assertion.
#
# Fixtures come from tools/inv_fixtures.py; the record comparisons from tools/invcheck.py, which
# pairs records with vcfdiff.py's passes (content, id, site, node set, carriers) so a renamed record is still
# compared with its counterpart.
#
#   call_invariants.sh <panvar-binary> <out-dir>
set -uo pipefail

BIN="${1:?usage: call_invariants.sh <panvar> <outdir>}"
OUT="${2:?}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TOOLS="$HERE/tools"
mkdir -p "$OUT"
OUT="$OUT/run.$$.$(date +%s)"
rm -rf "$OUT"; mkdir -p "$OUT"
fails=0; passes=0; xfails=0; xpasses=0
ok()  { printf "  ok    %s\n" "$1"; passes=$((passes + 1)); }
bad() { printf "  FAIL  %s\n" "$1"; fails=$((fails + 1)); }
IC() { python3 -I "$TOOLS/invcheck.py" "$@"; }

# check <label> <command...>: a hard assertion; the command's own lines say what differed
check() {
  local label="$1" out; shift
  out=$("$@" 2>&1)
  if [ $? -eq 0 ]; then ok "$label"; else bad "$label"; [ -n "$out" ] && printf '%s\n' "$out" | head -8; fi
}
# xcheck <plan-id> <label> <command...>: an expected failure until <plan-id> lands
xcheck() {
  local id="$1" label="$2" out; shift 2
  out=$("$@" 2>&1)
  if [ $? -eq 0 ]; then
    printf "  XPASS %s [%s]: passes now; make it a hard assertion\n" "$label" "$id"; xpasses=$((xpasses + 1))
  else
    printf "  xfail %s [%s]\n" "$label" "$id"; xfails=$((xfails + 1))
    [ -n "$out" ] && printf '%s\n' "$out" | head -4 | sed 's/^/  /'
  fi
}
# eq <label> <got> <want>
eq() { [ "$2" = "$3" ] && ok "$1 = $3" || bad "$1: expected $3, got ${2:-<missing>}"; }

# call_ <gfa> <csv> <ref> <out-prefix> [extra args...]; a failed run leaves no VCF behind
call_() {
  local g=$1 c=$2 r=$3 o=$4; shift 4
  rm -f "$o.region.vcf"
  "$BIN" call -i "$g" -c "$c" -r "$r" -o "$o" --no-per-bubble-vcf -q "$@" > "$o.log" 2>&1 \
    || printf "  (call exited non-zero for %s: %s)\n" "$(basename "$o")" "$(head -1 "$o.log")"
}

python3 -I "$TOOLS/inv_fixtures.py" "$OUT" || { echo "fixture generation failed"; exit 1; }

# ---------------------------------------------------------------- 1. bubble orientation
# ref 1+2+4+5; node 3 inserted after 2, node 2 deleted, node 2 inverted, node 2 replaced by node 6.
# The same bubble as the row (1,4) and as (4,1), i.e. a CSV that disagrees with -r about direction.
# A reverse row is unreachable through bubble/panphorte (they orient to the reference) but not
# through call's own -r, a hand-made CSV or a third-party one. Being fixed under C2/C3/C4, so this
# is a hard assertion.
echo "-- 1. a bubble given against the reference direction (C2/C3/C4)"
D="$OUT/rev"; R='ref#0#chr1:1001-1165'
for t in fwd rev; do
  call_ "$D/g.gfa" "$D/$t.bubbles.csv" "$R" "$D/$t"
done
check "the same records and carriers either way" IC carriers "$D/fwd.region.vcf" "$D/rev.region.vcf"
check "INSSEQ/DELSEQ/INVSEQ are spelled on the reference strand either way (C2)" \
  IC same "$D/fwd.region.vcf" "$D/rev.region.vcf" --fields INSSEQ,DELSEQ,INVSEQ
check "POS/END of the INS, DEL and INV do not move (C4)" \
  IC same "$D/fwd.region.vcf" "$D/rev.region.vcf" --fields POS,END
# A substitution is one record. As a DEL+INS pair its INS arm anchored on the walk-order preceding
# node, which a reverse bubble put at the far end of the replaced block. In the fixture above node 6
# matches node 2 here and there by chance, so the change measures under --min-sv-bp and there is no
# SUB record to compare. Node 2 on A and C only and node 6 on G and T only share no base: hapSUB is
# then one 60 bp replacement with a net change of 5 bp, a SUB record.
D2="$OUT/revsub"; mkdir -p "$D2"
awk -F'\t' -v OFS='\t' '$1 == "S" && $2 == "2" { gsub(/G/, "A", $3); gsub(/T/, "C", $3) }
                        $1 == "S" && $2 == "6" { gsub(/A/, "G", $3); gsub(/C/, "T", $3) } { print }' \
  "$D/g.gfa" > "$D2/g.gfa"
for t in fwd rev; do
  call_ "$D2/g.gfa" "$D/$t.bubbles.csv" "$R" "$D2/$t"
done
nsub() { awk -F'\t' '!/^#/ && $8 ~ /SVTYPE=SUB;/' "$1" | wc -l | tr -d ' '; }
eq "SUB records to compare (fwd rev)" "$(nsub "$D2/fwd.region.vcf") $(nsub "$D2/rev.region.vcf")" "1 1"
check "the substitution keeps its POS, END and sequences either way" \
  IC same "$D2/fwd.region.vcf" "$D2/rev.region.vcf" --fields POS,END,REF,ALT,SVLEN --only SVTYPE=SUB

# ---------------------------------------------------------------- 2. P-line order
# tie: two 60 bp deletions of identical sequence 30 bp apart merge into one record whose
# representative is a size tie. It used to go to the first member in path order, so POS moved; the
# tie now goes to the most-carried exact member, then to content (ORD).
# Carriers were invariant in every permuted real-data run, so those are hard. The shuffle seeds are
# ones that put hapB before hapA; a seed that keeps their order would XPASS by chance.
echo "-- 2. P-line order"
for fx in tie rev; do
  D="$OUT/$fx"; G="$D/g.gfa"; C="$D/g.bubbles.csv"; R='ref#0#chr1:1-170'
  [ "$fx" = rev ] && { C="$D/fwd.bubbles.csv"; R='ref#0#chr1:1001-1165'; }
  call_ "$G" "$C" "$R" "$D/ord0"
  for mode in rev 4 12; do
    python3 -I "$TOOLS/permute_paths.py" "$G" "$D/perm_$mode.gfa" "$mode"
    call_ "$D/perm_$mode.gfa" "$C" "$R" "$D/ord_$mode"
    check "$fx, paths $mode: the same carrier sets" IC carriers "$D/ord0.region.vcf" "$D/ord_$mode.region.vcf"
    if [ "$fx" = tie ]; then
      check "$fx, paths $mode: identical records (representative of a size tie)" \
        IC same "$D/ord0.region.vcf" "$D/ord_$mode.region.vcf" --fields ALL
    else
      check "$fx, paths $mode: identical records" \
        IC same "$D/ord0.region.vcf" "$D/ord_$mode.region.vcf" --fields ALL
    fi
  done
done

# ---------------------------------------------------------------- 3. node renumbering
# Every node and bubble id renumbered through a random bijection, the VCF mapped back with
# unmap_vcf.py. Anything still different depends on the numbering itself.
echo "-- 3. node and bubble renumbering"
for fx in rev tie chain sub; do
  D="$OUT/$fx"; G="$D/g.gfa"; C="$D/g.bubbles.csv"; R='ref#0#chr1:1-170'
  [ "$fx" = rev ] && { C="$D/fwd.bubbles.csv"; R='ref#0#chr1:1001-1165'; }
  [ "$fx" = chain ] && R='ref#0#chr1:1-100'
  [ "$fx" = sub ] && { G="$D/sub_delfirst.gfa"; R='ref#0#chr1:1-140'; }   # the SUB maps back too
  call_ "$G" "$C" "$R" "$D/num0"
  for seed in 3 4; do
    python3 -I "$TOOLS/renumber.py" "$G" "$C" "$D/num_$seed.gfa" "$D/num_$seed.bubbles.csv" "$seed"
    call_ "$D/num_$seed.gfa" "$D/num_$seed.bubbles.csv" "$R" "$D/num_$seed"
    for x in region; do
      python3 -I "$TOOLS/unmap_vcf.py" "$D/num_$seed.$x.vcf" "$D/num_$seed.gfa.map" "$D/num_${seed}_back.$x.vcf"
      check "$fx, seed $seed: the $x VCF maps back unchanged" \
        IC same "$D/num0.$x.vcf" "$D/num_${seed}_back.$x.vcf" --fields ALL
    done
  done
done

# ---------------------------------------------------------------- 4. single-linkage chaining
# Three 100 bp insertions on distinct nodes: A and C each differ from B in 15 bases at opposite ends,
# so A~B and B~C pass --merge-seq-identity 0.80 and A~C (0.70) does not. One record, and
# MERGE_DIAMETER must be there to show how far the chain reaches. The control drops B. Single linkage
# is the node comparison's merge (--no-align-alleles); events found by aligning alleles merge around a
# representative instead, and each must match it.
echo "-- 4. chaining A~B~C"
D="$OUT/chain"; R='ref#0#chr1:1-100'
call_ "$D/g.gfa" "$D/g.bubbles.csv" "$R" "$D/ch" --no-align-alleles
call_ "$D/g_ac.gfa" "$D/g_ac.bubbles.csv" "$R" "$D/ch_ac" --no-align-alleles
eq "control: A and C alone stay separate INS records" "$(IC count "$D/ch_ac.region.vcf" SVTYPE=INS)" "2"
eq "with B they chain into one INS record" "$(IC count "$D/ch.region.vcf" SVTYPE=INS)" "1"
eq "that record merges all three" "$(IC get "$D/ch.region.vcf" SVTYPE=INS NMERGED)" "3"
eq "and is carried by all three" "$(IC get "$D/ch.region.vcf" SVTYPE=INS AC)" "3"
md=$(IC get "$D/ch.region.vcf" SVTYPE=INS MERGE_DIAMETER)
[ -n "$md" ] && [ "$md" != "<absent>" ] && ok "MERGE_DIAMETER is reported ($md)" \
                                        || bad "the chained record has no MERGE_DIAMETER"

# ---------------------------------------------------------------- 5. substitution
# hapSUB replaces node 2 (60 bp) by node 6 (58 bp), hapDEL deletes node 2. The substitution is one
# record: its net change of -2 bp is under --min-sv-bp, so it is a SUB spelling both sequences, and the
# deletion stays a DEL of its own. As a DEL+INS pair the two DELs merged and the record reported a
# 60 bp loss for hapSUB, whose sequence is 2 bp shorter.
echo "-- 5. a substitution and a pure deletion of the same block"
D="$OUT/sub"; R='ref#0#chr1:1-140'
for o in delfirst subfirst; do call_ "$D/sub_$o.gfa" "$D/g.bubbles.csv" "$R" "$D/$o"; done
for o in delfirst subfirst; do
  eq "$o: the DEL is hapDEL's alone (AC:SVLEN:GT hapDEL:GT hapSUB)" \
     "$(IC get "$D/$o.region.vcf" SVTYPE=DEL AC):$(IC get "$D/$o.region.vcf" SVTYPE=DEL SVLEN):$(IC get "$D/$o.region.vcf" SVTYPE=DEL FMT:GT:hapDEL#1#c1):$(IC get "$D/$o.region.vcf" SVTYPE=DEL FMT:GT:hapSUB#1#c1)" "1:-60:1:0"
  eq "$o: the substitution is one SUB record carried by hapSUB (SVLEN:GT)" \
     "$(IC get "$D/$o.region.vcf" SVTYPE=SUB SVLEN):$(IC get "$D/$o.region.vcf" SVTYPE=SUB FMT:GT:hapSUB#1#c1)" "-2:1"
done
check "the records and carriers do not depend on P-line order" \
  IC same "$D/delfirst.region.vcf" "$D/subfirst.region.vcf" --fields POS,REF,ALT,SVLEN,GT

# ---------------------------------------------------------------- 6. short self-loop UNIT
# A 9 bp self-loop walked ten times by one haplotype: +81 bp, an SV-sized copy change on a unit far
# below --min-sv-bp. The UNIT route counts it; the original binary skipped units shorter than
# --min-sv-bp in default mode (and --cn fell back to PEAK).
echo "-- 6. a short self-loop expanded tenfold"
D="$OUT/loop"; R='ref#0#chr1:1-129'
for mode in default cn; do
  x=""; [ "$mode" = cn ] && x="--cn"
  # shellcheck disable=SC2086
  call_ "$D/g.gfa" "$D/g.bubbles.csv" "$R" "$D/$mode" $x
  V="$D/$mode.region.vcf"
  eq "$mode: one DUP record" "$(IC count "$V" SVTYPE=CNV)" "1"
  eq "$mode: it is a UNIT call" "$(IC get "$V" SVTYPE=CNV CN_METHOD)" "UNIT"
  eq "$mode: CN of the tenfold haplotype" "$(IC get "$V" SVTYPE=CNV FMT:CN:hapx10#1#c1)" "10"
  eq "$mode: CN of a reference-like haplotype" "$(IC get "$V" SVTYPE=CNV FMT:CN:hapR#1#c1)" "1"
  eq "$mode: SVLEN is the 81 bp gain" "$(IC get "$V" SVTYPE=CNV SVLEN)" "81"
done

echo
printf "call_invariants: %d passed, %d FAILED, %d xfail, %d XPASS\n" "$passes" "$fails" "$xfails" "$xpasses"
[ "$fails" -eq 0 ]
