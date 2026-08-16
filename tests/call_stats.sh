#!/usr/bin/env bash
# call_stats.sh - coordinate contract assertions for `panvar call`.
#
# real_smoke.sh checks that call RUNS and that its files exist; it asserts nothing about where a
# record lands. The coordinates are the part of a VCF every downstream consumer trusts silently, so
# they need a fixture whose every expected number is derivable by hand rather than read back from
# the caller's own output.
#
# The reference chain below starts at genomic 1, so a genomic coordinate IS an offset along the path:
#
#   n1 100 bp    1..100     flank
#   n2  70 bp  101..170     deleted by hapdel
#   n3 100 bp  171..270     flank
#   n4  30 bp  271..300     flank
#   n5  80 bp  301..380     inverted by hapinv
#   n6 100 bp  381..480     flank
#   n7  60 bp  481..540     flank
#   nI  60 bp     --        inserted by hapins, between n6 and n7
#
# A symbolic DEL/INV is anchored on the base PRECEDING the event, so REF carries one real base and
# the event occupies POS+1 onwards. That makes the expected records:
#
#   DEL  POS 100  END 170  SVLEN -70    REF = last base of n1
#   INV  POS 300  END 380  SVLEN  80    REF = last base of n4
#   INS  POS 480  END 480  SVLEN  60    REF = last base of n6
#
#   call_stats.sh <panvar-binary> <out-dir>
set -uo pipefail

BIN="${1:?usage: call_stats.sh <panvar> <outdir>}"
OUT="${2:?}"
mkdir -p "$OUT"
OUT="$OUT/run.$$.$(date +%s)"
rm -rf "$OUT"; mkdir -p "$OUT"
fails=0
ok()  { printf "  ok   %s\n" "$1"; }
bad() { printf "  FAIL %s\n" "$1"; fails=$((fails + 1)); }

# Deterministic sequence, distinct per node: repeat a fixed 40-mer rotated by $2 so no two nodes
# share content (identical flanks would give the snarl finder a different graph than intended).
blk() { local n=$1 r=${2:-0} u="ACGTTGCAATTCCGGATCAGGTCAAGCTTGACCTAGGACT" out=""
        u="${u:$r}${u:0:$r}"
        while [ ${#out} -lt "$n" ]; do out="$out$u"; done; printf '%s' "${out:0:$n}"; }

N1=$(blk 100 0); N2=$(blk 70 3);  N3=$(blk 100 7); N4=$(blk 30 11)
N5=$(blk 80 17); N6=$(blk 100 23); N7=$(blk 60 29); NI=$(blk 60 31)
REF_NAME='ref#0#chr1:1-540'

{ printf 'H\tVN:Z:1.0\n'
  printf "S\tn1\t%s\nS\tn2\t%s\nS\tn3\t%s\nS\tn4\t%s\n" "$N1" "$N2" "$N3" "$N4"
  printf "S\tn5\t%s\nS\tn6\t%s\nS\tn7\t%s\nS\tnI\t%s\n" "$N5" "$N6" "$N7" "$NI"
  # forward links along the reference chain, plus the deletion bypass and the insertion detour
  for e in "n1 + n2 +" "n2 + n3 +" "n1 + n3 +" "n3 + n4 +" "n4 + n5 +" "n5 + n6 +" \
           "n6 + n7 +" "n6 + nI +" "nI + n7 +" "n4 + n5 -" "n5 - n6 +"; do
    set -- $e; printf "L\t%s\t%s\t%s\t%s\t0M\n" "$1" "$2" "$3" "$4"
  done
  printf 'P\t%s\tn1+,n2+,n3+,n4+,n5+,n6+,n7+\t*\n' "$REF_NAME"
  printf 'P\thapref#1#chr1\tn1+,n2+,n3+,n4+,n5+,n6+,n7+\t*\n'
  printf 'P\thapdel#1#chr1\tn1+,n3+,n4+,n5+,n6+,n7+\t*\n'
  printf 'P\thapinv#1#chr1\tn1+,n2+,n3+,n4+,n5-,n6+,n7+\t*\n'
  printf 'P\thapins#1#chr1\tn1+,n2+,n3+,n4+,n5+,n6+,nI+,n7+\t*\n'; } > "$OUT/coord.gfa"

"$BIN" bubble -i "$OUT/coord.gfa" -r "$REF_NAME" -o "$OUT/cb" -q >/dev/null 2>&1
nb=$(tail -n +2 "$OUT/cb.bubbles.csv" 2>/dev/null | wc -l | tr -d ' ')
[ "$nb" = "3" ] && ok "the fixture presents three bubbles (deletion, inversion, insertion)" \
                || bad "expected 3 bubbles in the fixture, got $nb"

"$BIN" call -i "$OUT/cb.sorted.gfa" -c "$OUT/cb.bubbles.csv" -r "$REF_NAME" \
      -o "$OUT/cc" -q >/dev/null 2>&1
VCF="$OUT/cc.region.vcf"
[ -s "$VCF" ] || { bad "call wrote no region VCF"; printf "%d assertion(s) failed\n" "$((fails + 1))"; exit 1; }

# field <svtype> <column>  -- POS=2 REF=4; END/SVLEN come out of INFO
field() { awk -F'\t' -v t="$1" -v c="$2" '$0!~/^#/ && $8~("SVTYPE="t";") {print $c; exit}' "$VCF"; }
info()  { awk -F'\t' -v t="$1" -v k="$2" '$0!~/^#/ && $8~("SVTYPE="t";"){
            n=split($8,a,";"); for(i=1;i<=n;i++){split(a[i],kv,"="); if(kv[1]==k){print kv[2]; exit}}}' "$VCF"; }

check() { # check <label> <got> <want>
  [ "$2" = "$3" ] && ok "$1 = $3" || bad "$1: expected $3, got ${2:-<missing>}"; }

# The event coordinates, each derivable from the node lengths at the top of this file.
check "DEL POS (base before the deleted span)"  "$(field DEL 2)"      "100"
check "DEL END (last deleted base)"             "$(info  DEL END)"    "170"
check "DEL SVLEN"                               "$(info  DEL SVLEN)"  "-70"
check "INV POS (base before the inverted span)" "$(field INV 2)"      "300"
check "INV END (last inverted base)"            "$(info  INV END)"    "380"
check "INV SVLEN"                               "$(info  INV SVLEN)"  "80"
check "INS POS (base before the insertion)"     "$(field INS 2)"      "480"
check "INS END (INS does not span reference)"   "$(info  INS END)"    "480"
check "INS SVLEN"                               "$(info  INS SVLEN)"  "60"

# REF must be the reference base AT POS, which is the last base of the preceding flank. Anchoring on
# the first affected base instead would put the first base of n2/n5 here and still look well formed,
# so this is what separates the two conventions.
check "DEL REF is the last base of n1" "$(field DEL 4)" "${N1: -1}"
check "INV REF is the last base of n4" "$(field INV 4)" "${N4: -1}"
check "INS REF is the last base of n6" "$(field INS 4)" "${N6: -1}"

# END-POS must equal the deleted length, i.e. the record's own interval closes on the event.
dp=$(field DEL 2); de=$(info DEL END)
[ -n "$dp" ] && [ -n "$de" ] && [ "$((de - dp))" = "70" ] \
  && ok "DEL END-POS equals the deleted length (70)" \
  || bad "DEL END-POS: expected 70, got $((${de:-0} - ${dp:-0}))"
ip=$(field INV 2); ie=$(info INV END)
[ -n "$ip" ] && [ -n "$ie" ] && [ "$((ie - ip))" = "80" ] \
  && ok "INV END-POS equals the inverted length (80)" \
  || bad "INV END-POS: expected 80, got $((${ie:-0} - ${ip:-0}))"

# DELSEQ/INVSEQ are the sequences the record says it changes; they must be what the reference
# actually holds at POS+1..END. Independent of merging, and it pins the offset directly.
[ "$(info DEL DELSEQ)" = "$N2" ] && ok "DELSEQ is exactly n2, the deleted node" \
                                 || bad "DELSEQ does not equal n2"
[ "$(info INS INSSEQ)" = "$NI" ] && ok "INSSEQ is exactly nI, the inserted node" \
                                 || bad "INSSEQ does not equal nI"

printf "\n"
if [ "$fails" -eq 0 ]; then printf "call_stats: all assertions passed\n"; exit 0; fi
printf "call_stats: %d assertion(s) failed\n" "$fails"; exit 1
