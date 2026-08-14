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

# ---------------------------------------------------------------- a folded REP is held fixed
# panphorte's self-loop REP node carries copy number. refine re-aligns the residual flanks around it
# and must leave the REP run itself verbatim: same node, same count, same orientation.
UNIT=$(blk 64)
{ printf 'H\tVN:Z:1.0\n'
  printf "S\t1\t%s\nS\t2\t%s\nS\t3\t%s\nS\t4\t%s\nS\t5\t%s\n" "$BB1" "$UNIT" "$(blk 40)" "$(blk 38)" "$BB2"
  for e in "1 2" "2 2" "2 3" "2 4" "3 5" "4 5"; do
    set -- $e; printf "L\t%s\t+\t%s\t+\t0M\n" "$1" "$2"
  done
  printf 'P\tref\t1+,2+,2+,2+,3+,5+\t*\n'
  printf 'P\talt\t1+,2+,2+,2+,4+,5+\t*\n'
  printf 'P\tthird\t1+,2+,2+,4+,5+\t*\n'; } > "$OUT/rep.gfa"
"$BIN" bubble -i "$OUT/rep.gfa" -r ref -o "$OUT/repb" --min-variant-bp 0 -q >/dev/null 2>&1
before_ref=$(spell "$OUT/repb.sorted.gfa" ref)
"$BIN" refine -i "$OUT/repb.sorted.gfa" --bubbles-csv-in "$OUT/repb.bubbles.csv" -r ref \
       -o "$OUT/repr" -q >/dev/null 2>&1
if [ -s "$OUT/repr.normalized.sorted.gfa" ]; then
  [ "$(spell "$OUT/repr.normalized.sorted.gfa" ref)" = "$before_ref" ] \
    && ok "a folded REP region preserves the reference's spelled sequence" \
    || bad "the REP region changed the reference's sequence"
  selfloop=$(awk '$1=="L" && $2==$4 {print $2}' "$OUT/repr.normalized.sorted.gfa" | head -1)
  [ -n "$selfloop" ] \
    && ok "the self-loop REP node survives refinement" \
    || bad "the self-loop edge was lost, so copy number is no longer recoverable"
  n=$(awk -v r="$selfloop" '$1=="P" && $2=="ref"{n=split($3,a,","); c=0
        for(i=1;i<=n;i++){s=a[i]; if(substr(s,1,length(s)-1)==r) c++} print c}' \
        "$OUT/repr.normalized.sorted.gfa")
  [ "$n" = "3" ] \
    && ok "the REP run keeps its copy count (3x)" \
    || bad "the REP is traversed ${n}x after refinement, expected 3"
else
  bad "the REP fixture produced no output"
fi

# ---------------------------------------------------------------- an unfolded non-REP revisit is skipped
# A plain interior node visited twice is an unfolded copy-number signal. POA would linearize it and
# destroy what `call` reconstructs, so the region must be skipped rather than rebuilt.
# Node 2 is revisited via 3 and carries NO self-loop, so panphorte never folded it: it is an unfolded
# copy-number signal, not a REP. (A self-loop would make it a REP and a different guard would fire.)
{ printf 'H\tVN:Z:1.0\n'
  printf "S\t1\t%s\nS\t2\t%s\nS\t3\t%s\nS\t4\t%s\nS\t5\t%s\n" \
         "$BB1" "$(blk 120)" "$(blk 30)" "$BB2" "$(blk 118)"
  for e in "1 2" "2 3" "3 2" "2 4" "1 5" "5 4"; do
    set -- $e; printf "L\t%s\t+\t%s\t+\t0M\n" "$1" "$2"
  done
  printf 'P\tref\t1+,2+,3+,2+,4+\t*\nP\talt\t1+,5+,4+\t*\n'; } > "$OUT/unf.gfa"
"$BIN" bubble -i "$OUT/unf.gfa" -r ref -o "$OUT/unfb" --min-variant-bp 0 -q >/dev/null 2>&1
"$BIN" refine -i "$OUT/unfb.sorted.gfa" --bubbles-csv-in "$OUT/unfb.bubbles.csv" -r ref \
       -o "$OUT/unfr" > "$OUT/unfr.log" 2>&1
grep -q "unfolded dup" "$OUT/unfr.refine.report.tsv" 2>/dev/null \
  && ok "an unfolded non-REP revisit is skipped, with the reason recorded" \
  || bad "the unfolded duplication was not skipped: $(tail -1 "$OUT/unfr.refine.report.tsv" 2>/dev/null)"

# ---------------------------------------------------------------- partial traversal
# A path entering the interior without spanning both anchors cannot be rewritten. Retaining its nodes
# also retains the old edges, so refined and unrefined topology coexist -- skip is the default.
{ printf 'H\tVN:Z:1.0\n'
  printf "S\t1\t%s\nS\t2\t%s\nS\t3\t%s\nS\t4\t%s\nS\t5\t%s\n" \
         "$BB1" "$(blk 120)" "$(blk 118)" "$BB2" "$(blk 90)"
  for e in "1 2" "1 3" "2 4" "3 4" "5 2"; do set -- $e; printf "L\t%s\t+\t%s\t+\t0M\n" "$1" "$2"; done
  printf 'P\tref\t1+,2+,4+\t*\nP\talt\t1+,3+,4+\t*\nP\tpartial\t5+,2+,4+\t*\n'; } > "$OUT/part.gfa"
"$BIN" bubble -i "$OUT/part.gfa" -r ref -o "$OUT/partb" --min-variant-bp 0 -q >/dev/null 2>&1
"$BIN" refine -i "$OUT/partb.sorted.gfa" --bubbles-csv-in "$OUT/partb.bubbles.csv" -r ref \
       -o "$OUT/ps" > "$OUT/ps.log" 2>&1
[ "$(rebuilt "$OUT/ps.log")" = "0 1" ] \
  && ok "a region with a partial traverser is skipped by default" \
  || bad "partial traversal gave '$(rebuilt "$OUT/ps.log")', expected '0 1'"
"$BIN" refine -i "$OUT/partb.sorted.gfa" --bubbles-csv-in "$OUT/partb.bubbles.csv" -r ref \
       -o "$OUT/pr" --partial-path-policy retain > "$OUT/pr.log" 2>&1
[ "$(rebuilt "$OUT/pr.log")" = "1 0" ] \
  && ok "--partial-path-policy retain rebuilds it instead" \
  || bad "retain gave '$(rebuilt "$OUT/pr.log")', expected '1 0'"

# ---------------------------------------------------------------- W-line graphs
# parse_gfa names a W path sample#hap#seqid:start-end. refine kept a private naming rule that dropped
# the :start-end suffix, so the resolved reference matched no path and every W graph failed.
W="$(cd "$(dirname "$0")" && pwd)/synthetic_data/syn_w.gfa"
if [ -f "$W" ]; then
  wref=$(awk '$1=="W"{print $2"#"$3"#"$4":"$5"-"$6; exit}' "$W")
  "$BIN" bubble -i "$W" -r "$wref" -o "$OUT/wb" --min-variant-bp 0 -q >/dev/null 2>&1
  "$BIN" refine -i "$OUT/wb.sorted.gfa" --bubbles-csv-in "$OUT/wb.bubbles.csv" -r "$wref" \
         -o "$OUT/wr" -q >/dev/null 2>&1
  [ "$?" -eq 0 ] && ok "a W-line graph resolves its reference and refines" \
                 || bad "refine failed on a W-line graph"
else
  ok "W-line fixture not present, skipped"
fi

# ---------------------------------------------------------------- duplicate carriers must not decide
# abPOA is handed DISTINCT sequences, so adding an identical haplotype changes no POA input. Summing
# bases across carriers made cohort composition decide whether a region was rebuilt.
{ printf 'H\tVN:Z:1.0\n'
  printf "S\t1\t%s\nS\t2\t%s\nS\t3\t%s\nS\t4\t%s\n" "$BB1" "$(blk 120)" "$(blk 118)" "$BB2"
  for e in "1 2" "1 3" "2 4" "3 4"; do set -- $e; printf "L\t%s\t+\t%s\t+\t0M\n" "$1" "$2"; done
  printf 'P\tref\t1+,2+,4+\t*\nP\talt\t1+,3+,4+\t*\n'; } > "$OUT/dup1.gfa"
{ cat "$OUT/dup1.gfa"
  for i in 1 2 3 4 5 6 7 8; do printf "P\tcopy%s\t1+,3+,4+\t*\n" "$i"; done; } > "$OUT/dup2.gfa"
for g in dup1 dup2; do
  "$BIN" bubble -i "$OUT/$g.gfa" -r ref -o "$OUT/${g}b" --min-variant-bp 0 -q >/dev/null 2>&1
  "$BIN" refine -i "$OUT/${g}b.sorted.gfa" --bubbles-csv-in "$OUT/${g}b.bubbles.csv" -r ref \
         -o "$OUT/${g}r" --max-poa-work 200000 > "$OUT/${g}r.log" 2>&1
done
[ "$(rebuilt "$OUT/dup1r.log")" = "$(rebuilt "$OUT/dup2r.log")" ] \
  && ok "replicating identical carriers does not change the POA decision" \
  || bad "decision changed with duplicate carriers: '$(rebuilt "$OUT/dup1r.log")' vs '$(rebuilt "$OUT/dup2r.log")'"

# ---------------------------------------------------------------- shuffled bubble rows
# Component order must come from the graph, not from CSV row order.
head -1 "$OUT/disb.bubbles.csv" > "$OUT/shuf.csv"
tail -n +2 "$OUT/disb.bubbles.csv" | tail -r >> "$OUT/shuf.csv" 2>/dev/null || \
  { tail -n +2 "$OUT/disb.bubbles.csv" | tac >> "$OUT/shuf.csv"; }
"$BIN" refine -i "$OUT/disb.sorted.gfa" --bubbles-csv-in "$OUT/shuf.csv" -r ref -o "$OUT/sh" \
       -q >/dev/null 2>&1
cmp -s "$OUT/d1.normalized.sorted.gfa" "$OUT/sh.normalized.sorted.gfa" \
  && ok "reordering the bubbles CSV gives byte-identical output" \
  || bad "output depends on the order of rows in the bubbles CSV"

# ---------------------------------------------------------------- the re-snarl threshold is usable
"$BIN" refine -i "$OUT/disb.sorted.gfa" --bubbles-csv-in "$OUT/disb.bubbles.csv" -r ref \
       -o "$OUT/rs" --resnarl-min-variant-bp 0 -q >/dev/null 2>&1
[ "$?" -eq 0 ] && [ -s "$OUT/rs.bubbles.csv" ] \
  && ok "--resnarl-min-variant-bp is accepted and produces a CSV" \
  || bad "--resnarl-min-variant-bp is advertised but not parsed"

# ---------------------------------------------------------------- an output may not name the GTF
printf '#gtf\n' > "$OUT/genes.gtf"
"$BIN" refine -i "$OUT/disb.sorted.gfa" --bubbles-csv-in "$OUT/disb.bubbles.csv" -r ref \
       -o "$OUT/genes" --gtf "$OUT/genes.gtf" -q >/dev/null 2>&1
# <prefix>.bandage_genes.csv vs --gtf genes.gtf do not collide; the direct collision is the check below
cp "$OUT/disb.bubbles.csv" "$OUT/alias.bandage_genes.csv"
"$BIN" refine -i "$OUT/disb.sorted.gfa" --bubbles-csv-in "$OUT/disb.bubbles.csv" -r ref \
       -o "$OUT/alias" --gtf "$OUT/alias.bandage_genes.csv" -q >/dev/null 2>&1
[ "$?" -ne 0 ] && ok "an output that would overwrite the --gtf input is refused" \
               || bad "the --gtf input was left open to being overwritten"

echo
if [ "$fails" -eq 0 ]; then echo "refine_stats: all assertions passed"; exit 0; fi
echo "refine_stats: $fails assertion(s) FAILED"; exit 1
