#!/usr/bin/env bash
# refine_stats.sh - contract assertions for `panvar refine`.
#
# real_smoke.sh only checks that refine's files exist and that `call` runs afterwards, which passes
# whether or not the rewrite was correct.
#
# Fixtures carry a backbone longer than their alleles: the cactus decomposition roots at the longest
# path, so on a graph where the bubble IS most of the graph the site is absorbed into the root and no
# snarl is reported. That is a property of the decomposition, not of refine.
#
#   refine_stats.sh <panvar-binary> <out-dir>
set -uo pipefail

BIN="${1:?usage: refine_stats.sh <panvar> <outdir>}"
OUT="${2:?}"
mkdir -p "$OUT"
OUT="$OUT/run.$$.$(date +%s)"
rm -rf "$OUT"; mkdir -p "$OUT"
fails=0
ok()  { printf "  ok   %s\n" "$1"; }
bad() { printf "  FAIL %s\n" "$1"; fails=$((fails + 1)); }
rebuilt() { sed -n 's/.*rebuilt \([0-9]*\) region(s), skipped \([0-9]*\).*/\1 \2/p' "$1"; }

# Deterministic pseudo-random sequence without python: repeat a fixed 40-mer.
blk() { local n=$1 out="" u="ACGTTGCAATTCCGGATCAGGTCAAGCTTGACCTAGGACT"
        while [ ${#out} -lt "$n" ]; do out="$out$u"; done; printf '%s' "${out:0:$n}"; }

BB1=$(blk 3000); BB2=$(blk 3000); BB3=$(blk 3000)
A1=$(blk 120);   A2=$(blk 118);   A3=$(blk 117)

# Two DISJOINT bubbles separated by a backbone node, so they are two components.
{ printf 'H\tVN:Z:1.0\n'
  printf "S\t1\t%s\nS\t2\t%s\nS\t3\t%s\nS\t4\t%s\nS\t5\t%s\nS\t6\t%s\nS\t7\t%s\nS\t8\t%s\n" \
         "$BB1" "$A1" "$A2" "$BB2" "$BB3" "$A1" "$A3" "$BB1"
  for e in "1 2" "1 3" "2 4" "3 4" "4 5" "5 6" "5 7" "6 8" "7 8"; do
    set -- $e; printf "L\t%s\t+\t%s\t+\t0M\n" "$1" "$2"
  done
  printf 'P\tref\t1+,2+,4+,5+,6+,8+\t*\n'
  printf 'P\talt\t1+,3+,4+,5+,7+,8+\t*\n'
  printf 'P\tmix\t1+,2+,4+,5+,7+,8+\t*\n'; } > "$OUT/dis.gfa"
"$BIN" bubble -i "$OUT/dis.gfa" -r ref -o "$OUT/disb" --min-variant-bp 0 -q >/dev/null 2>&1
[ "$(tail -n +2 "$OUT/disb.bubbles.csv" | wc -l | tr -d ' ')" = "2" ] \
  && ok "the fixture presents two disjoint bubbles" \
  || bad "expected 2 bubbles in the fixture, got $(tail -n +2 "$OUT/disb.bubbles.csv" | wc -l | tr -d ' ')"

# ---------------------------------------------------------------- targeted multi-bubble selection
# Selected ids must go through the same component grouping as auto mode. Forcing them into one region
# assumed they were adjacent, so two disjoint bubbles became one region spanning disconnected anchors
# and failed as a unit: "rebuilt 0, skipped 1" where auto rebuilt both.
"$BIN" refine -i "$OUT/disb.sorted.gfa" --bubbles-csv-in "$OUT/disb.bubbles.csv" -r ref \
       -o "$OUT/auto" > "$OUT/auto.log" 2>&1
"$BIN" refine -i "$OUT/disb.sorted.gfa" --bubbles-csv-in "$OUT/disb.bubbles.csv" -r ref \
       -o "$OUT/sel" --bubble-id 1,2 > "$OUT/sel.log" 2>&1
[ "$(rebuilt "$OUT/auto.log")" = "2 0" ] \
  && ok "auto mode rebuilds both disjoint regions" \
  || bad "auto mode gave '$(rebuilt "$OUT/auto.log")', expected '2 0'"
[ "$(rebuilt "$OUT/sel.log")" = "2 0" ] \
  && ok "--bubble-id 1,2 splits the selection into components, as auto mode does" \
  || bad "--bubble-id 1,2 gave '$(rebuilt "$OUT/sel.log")', expected '2 0'"

# ---------------------------------------------------------------- a requested id that is not there
rm -f "$OUT/miss".*
"$BIN" refine -i "$OUT/disb.sorted.gfa" --bubbles-csv-in "$OUT/disb.bubbles.csv" -r ref \
       -o "$OUT/miss" --bubble-id 999 -q >/dev/null 2>&1
[ "$?" -ne 0 ] && ok "--bubble-id naming nothing is refused" \
               || bad "--bubble-id 999 succeeded, rebuilding nothing and writing a full output family"
[ "$(ls "$OUT/miss".* 2>/dev/null | wc -l | tr -d ' ')" = "0" ] \
  && ok "a refused run writes no output family" \
  || bad "the refused run left $(ls "$OUT/miss".* 2>/dev/null | wc -l | tr -d ' ') file(s)"

# ---------------------------------------------------------------- sequence losslessness
# refine re-aligns interiors and rebuilds nodes. Every haplotype must still spell exactly what it did.
spell() {
  awk -v want="$2" '
    $1=="S"{seq[$2]=$3}
    $1=="P" && $2==want{
      n=split($3,a,","); out=""
      for(i=1;i<=n;i++){ s=a[i]; id=substr(s,1,length(s)-1); o=substr(s,length(s)); t=seq[id]
        if(o=="-"){ r=""; for(j=length(t);j>0;j--){ c=substr(t,j,1)
            r=r (c=="A"?"T":c=="T"?"A":c=="C"?"G":c=="G"?"C":c) } t=r }
        out=out t }
      print out }' "$1"
}
allsame=1
for h in ref alt mix; do
  [ "$(spell "$OUT/disb.sorted.gfa" "$h")" = "$(spell "$OUT/auto.normalized.sorted.gfa" "$h")" ] || allsame=0
done
[ "$allsame" = "1" ] \
  && ok "every haplotype spells exactly the same sequence after refinement" \
  || bad "refinement changed a haplotype's sequence"

# ---------------------------------------------------------------- determinism
"$BIN" refine -i "$OUT/disb.sorted.gfa" --bubbles-csv-in "$OUT/disb.bubbles.csv" -r ref \
       -o "$OUT/d1" -q >/dev/null 2>&1
"$BIN" refine -i "$OUT/disb.sorted.gfa" --bubbles-csv-in "$OUT/disb.bubbles.csv" -r ref \
       -o "$OUT/d2" -q >/dev/null 2>&1
cmp -s "$OUT/d1.normalized.sorted.gfa" "$OUT/d2.normalized.sorted.gfa" \
  && ok "two runs of the same input are byte-identical" \
  || bad "repeated runs differ; component iteration order is leaking into the output"

# ---------------------------------------------------------------- reference spelling
"$BIN" refine -i "$OUT/disb.sorted.gfa" --bubbles-csv-in "$OUT/disb.bubbles.csv" -r REF \
       -o "$OUT/up" -q >/dev/null 2>&1
cmp -s "$OUT/d1.bubbles.csv" "$OUT/up.bubbles.csv" \
  && ok "a case-insensitive reference alias gives the same result as the exact name" \
  || bad "-r REF and -r ref produced different output"

# An ambiguous reference must be an error rather than the first file-order match.
{ printf 'H\tVN:Z:1.0\n'
  printf "S\t1\t%s\nS\t2\t%s\nS\t3\t%s\nS\t4\t%s\n" "$BB1" "$A1" "$A2" "$BB2"
  for e in "1 2" "1 3" "2 4" "3 4"; do set -- $e; printf "L\t%s\t+\t%s\t+\t0M\n" "$1" "$2"; done
  printf 'P\thapA\t1+,2+,4+\t*\nP\thapB\t1+,3+,4+\t*\n'; } > "$OUT/amb.gfa"
"$BIN" bubble -i "$OUT/amb.gfa" -r hapA -o "$OUT/ambb" --min-variant-bp 0 -q >/dev/null 2>&1
"$BIN" refine -i "$OUT/ambb.sorted.gfa" --bubbles-csv-in "$OUT/ambb.bubbles.csv" -r hap \
       -o "$OUT/ambo" -q >/dev/null 2>&1
[ "$?" -ne 0 ] && ok "an ambiguous reference is refused rather than resolved by file order" \
               || bad "-r hap matched two paths and was accepted"

# ---------------------------------------------------------------- the graph and CSV contracts
{ printf 'H\tVN:Z:1.0\n'
  printf "S\t1\t%s\nS\t2\t%s\nS\t3\t%s\nS\t4\t%s\n" "$BB1" "$A1" "$A2" "$BB2"
  for e in "1 2" "1 3" "3 4"; do set -- $e; printf "L\t%s\t+\t%s\t+\t0M\n" "$1" "$2"; done
  printf 'P\tref\t1+,2+,4+\t*\nP\talt\t1+,3+,4+\t*\n'; } > "$OUT/nolink.gfa"
"$BIN" refine -i "$OUT/nolink.gfa" --bubbles-csv-in "$OUT/disb.bubbles.csv" -r ref \
       -o "$OUT/nl" -q >/dev/null 2>&1
[ "$?" -ne 0 ] && ok "a path step with no link behind it is refused" \
               || bad "a malformed graph was accepted"

printf 'bubble_id,source,sink,path_support,min_inside_bp,max_inside_bp,inside_nodes\n9,900,902,2,10,10,"901"\n' \
  > "$OUT/ghost.csv"
"$BIN" refine -i "$OUT/disb.sorted.gfa" --bubbles-csv-in "$OUT/ghost.csv" -r ref -o "$OUT/gh" \
       -q >/dev/null 2>&1
[ "$?" -ne 0 ] && ok "a bubbles CSV naming nodes absent from the GFA is refused" \
               || bad "a mismatched CSV ran to completion"

"$BIN" refine -i "$OUT/disb.sorted.gfa" --bubbles-csv-in "$OUT/disb.bubbles.csv" -r ref \
       -o "$OUT/disb" -q >/dev/null 2>&1
[ "$?" -ne 0 ] && ok "an output that would overwrite an input is refused" \
               || bad "an output was allowed to overwrite an input"

# ---------------------------------------------------------------- the POA cost guard bounds the worst case
# A guard on the MEDIAN let one very long outlier into abPOA regardless of its size. --max-poa-bp 100
# must skip a region whose longest interior is 120 bp even though most are shorter.
"$BIN" refine -i "$OUT/disb.sorted.gfa" --bubbles-csv-in "$OUT/disb.bubbles.csv" -r ref \
       -o "$OUT/guard" --max-poa-bp 100 > "$OUT/guard.log" 2>&1
[ "$(rebuilt "$OUT/guard.log")" = "0 2" ] \
  && ok "--max-poa-bp bounds the longest sequence, not the median" \
  || bad "guard gave '$(rebuilt "$OUT/guard.log")', expected '0 2'"
grep -q "longest" "$OUT/guard.refine.report.tsv" \
  && ok "the report names the guard that fired" \
  || bad "the report does not explain the skip: $(tail -1 "$OUT/guard.refine.report.tsv")"

# ---------------------------------------------------------------- the decision report
[ -s "$OUT/auto.refine.report.tsv" ] \
  && ok "a decision report is written" \
  || bad "no refine.report.tsv was written"
[ "$(tail -n +2 "$OUT/auto.refine.report.tsv" | wc -l | tr -d ' ')" = "2" ] \
  && ok "the report carries one row per region" \
  || bad "the report has $(tail -n +2 "$OUT/auto.refine.report.tsv" | wc -l | tr -d ' ') rows, expected 2"

echo
if [ "$fails" -eq 0 ]; then echo "refine_stats: all assertions passed"; exit 0; fi
echo "refine_stats: $fails assertion(s) FAILED"; exit 1
