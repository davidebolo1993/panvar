#!/usr/bin/env bash
# refine_stats.sh - contract assertions for `panvar refine`.
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
rebuilt() { echo "$(sed -n 's/.*\] rebuilt: \([0-9]*\).*/\1/p' "$1") $(sed -n 's/.*\] skipped: \([0-9]*\).*/\1/p' "$1")"; }

# Deterministic pseudo-random sequence without python (a Park-Miller generator, exact in any awk). Not a
# repeated motif: that is a tandem repeat, which the re-snarl extends sites to.
blk() { awk -v n="$1" 'BEGIN{ b="ACGT"; x=1; for(i=0;i<n;i++){ x=(x*16807)%2147483647
          printf "%s", substr(b, int(x*4/2147483647)+1, 1) } }'; }

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
"$BIN" bubble -i "$OUT/dis.gfa" -r ref -o "$OUT/disb" --min-variant-bp 0 --merge-nearby-bp 0 -q >/dev/null 2>&1
[ "$(tail -n +2 "$OUT/disb.bubbles.csv" | wc -l | tr -d ' ')" = "2" ] \
  && ok "the fixture presents two disjoint bubbles" \
  || bad "expected 2 bubbles in the fixture, got $(tail -n +2 "$OUT/disb.bubbles.csv" | wc -l | tr -d ' ')"

# ---------------------------------------------------------------- targeted multi-bubble selection
# Selected ids go through the same component grouping as auto mode: forced into one region, two disjoint
# bubbles would become one region spanning disconnected anchors and fail as a unit.
"$BIN" refine -i "$OUT/disb.sorted.gfa" --bubbles-csv "$OUT/disb.bubbles.csv" -r ref \
       -o "$OUT/auto" > "$OUT/auto.log" 2>&1
"$BIN" refine -i "$OUT/disb.sorted.gfa" --bubbles-csv "$OUT/disb.bubbles.csv" -r ref \
       -o "$OUT/sel" --bubble-id 1,2 > "$OUT/sel.log" 2>&1
[ "$(rebuilt "$OUT/auto.log")" = "2 0" ] \
  && ok "auto mode rebuilds both disjoint regions" \
  || bad "auto mode gave '$(rebuilt "$OUT/auto.log")', expected '2 0'"
[ "$(rebuilt "$OUT/sel.log")" = "2 0" ] \
  && ok "--bubble-id 1,2 splits the selection into components, as auto mode does" \
  || bad "--bubble-id 1,2 gave '$(rebuilt "$OUT/sel.log")', expected '2 0'"

# ---------------------------------------------------------------- a requested id that is not there
rm -f "$OUT/miss".*
"$BIN" refine -i "$OUT/disb.sorted.gfa" --bubbles-csv "$OUT/disb.bubbles.csv" -r ref \
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
"$BIN" refine -i "$OUT/disb.sorted.gfa" --bubbles-csv "$OUT/disb.bubbles.csv" -r ref \
       -o "$OUT/d1" -q >/dev/null 2>&1
"$BIN" refine -i "$OUT/disb.sorted.gfa" --bubbles-csv "$OUT/disb.bubbles.csv" -r ref \
       -o "$OUT/d2" -q >/dev/null 2>&1
cmp -s "$OUT/d1.normalized.sorted.gfa" "$OUT/d2.normalized.sorted.gfa" \
  && ok "two runs of the same input are byte-identical" \
  || bad "repeated runs differ; component iteration order is leaking into the output"

# ---------------------------------------------------------------- reference spelling
"$BIN" refine -i "$OUT/disb.sorted.gfa" --bubbles-csv "$OUT/disb.bubbles.csv" -r REF \
       -o "$OUT/up" -q >/dev/null 2>&1
cmp -s "$OUT/d1.bubbles.csv" "$OUT/up.bubbles.csv" \
  && ok "a case-insensitive reference alias gives the same result as the exact name" \
  || bad "-r REF and -r ref produced different output"

# An ambiguous reference must be an error rather than the first file-order match.
{ printf 'H\tVN:Z:1.0\n'
  printf "S\t1\t%s\nS\t2\t%s\nS\t3\t%s\nS\t4\t%s\n" "$BB1" "$A1" "$A2" "$BB2"
  for e in "1 2" "1 3" "2 4" "3 4"; do set -- $e; printf "L\t%s\t+\t%s\t+\t0M\n" "$1" "$2"; done
  printf 'P\thapA\t1+,2+,4+\t*\nP\thapB\t1+,3+,4+\t*\n'; } > "$OUT/amb.gfa"
"$BIN" bubble -i "$OUT/amb.gfa" -r hapA -o "$OUT/ambb" --min-variant-bp 0 --merge-nearby-bp 0 -q >/dev/null 2>&1
"$BIN" refine -i "$OUT/ambb.sorted.gfa" --bubbles-csv "$OUT/ambb.bubbles.csv" -r hap \
       -o "$OUT/ambo" -q >/dev/null 2>&1
[ "$?" -ne 0 ] && ok "an ambiguous reference is refused rather than resolved by file order" \
               || bad "-r hap matched two paths and was accepted"

# ---------------------------------------------------------------- the graph and CSV contracts
{ printf 'H\tVN:Z:1.0\n'
  printf "S\t1\t%s\nS\t2\t%s\nS\t3\t%s\nS\t4\t%s\n" "$BB1" "$A1" "$A2" "$BB2"
  for e in "1 2" "1 3" "3 4"; do set -- $e; printf "L\t%s\t+\t%s\t+\t0M\n" "$1" "$2"; done
  printf 'P\tref\t1+,2+,4+\t*\nP\talt\t1+,3+,4+\t*\n'; } > "$OUT/nolink.gfa"
# Only the site on nodes 1-4, so the missing link is the one thing wrong.
head -2 "$OUT/disb.bubbles.csv" > "$OUT/nolink.bubbles.csv"
"$BIN" refine -i "$OUT/nolink.gfa" --bubbles-csv "$OUT/nolink.bubbles.csv" -r ref \
       -o "$OUT/nl" -q >/dev/null 2>"$OUT/nl.err"
grep -q "steps from 2+ to 4+ but the graph has no such link" "$OUT/nl.err" \
  && ok "a path step with no link behind it is refused" \
  || bad "a path step with no link: $(head -1 "$OUT/nl.err")"

printf 'bubble_id,source,sink,path_support,min_inside_bp,max_inside_bp,inside_nodes\n9,900,902,2,10,10,"901"\n' \
  > "$OUT/ghost.csv"
"$BIN" refine -i "$OUT/disb.sorted.gfa" --bubbles-csv "$OUT/ghost.csv" -r ref -o "$OUT/gh" \
       -q >/dev/null 2>&1
[ "$?" -ne 0 ] && ok "a bubbles CSV naming nodes absent from the GFA is refused" \
               || bad "a mismatched CSV ran to completion"

"$BIN" refine -i "$OUT/disb.sorted.gfa" --bubbles-csv "$OUT/disb.bubbles.csv" -r ref \
       -o "$OUT/disb" -q >/dev/null 2>&1
[ "$?" -ne 0 ] && ok "an output that would overwrite an input is refused" \
               || bad "an output was allowed to overwrite an input"

# ---------------------------------------------------------------- the POA cost guard bounds the worst case
# The guard is on the longest sequence, not the median, so one very long outlier cannot reach abPOA.
# The interiors are 118 and 120 bp in one region, 117 and 120 bp in the other: --max-poa-bp 119 skips
# both on their longest, where a limit on the median or the shortest would rebuild at least one.
"$BIN" refine -i "$OUT/disb.sorted.gfa" --bubbles-csv "$OUT/disb.bubbles.csv" -r ref \
       -o "$OUT/guard" --max-poa-bp 119 > "$OUT/guard.log" 2>&1
[ "$(rebuilt "$OUT/guard.log")" = "0 2" ] \
  && ok "--max-poa-bp bounds the longest sequence, not the median" \
  || bad "guard gave '$(rebuilt "$OUT/guard.log")', expected '0 2'"
[ "$(awk -F'\t' 'NR>1{print $6}' "$OUT/guard.refine.report.tsv" | sort -u)" = "segment_too_long" ] \
  && ok "the report names the guard that fired" \
  || bad "the report does not explain the skip: $(tail -1 "$OUT/guard.refine.report.tsv")"

# ---------------------------------------------------------------- the decision report
[ -s "$OUT/auto.refine.report.tsv" ] \
  && ok "a decision report is written" \
  || bad "no refine.report.tsv was written"
[ "$(tail -n +2 "$OUT/auto.refine.report.tsv" | wc -l | tr -d ' ')" = "2" ] \
  && ok "the report carries one row per region" \
  || bad "the report has $(tail -n +2 "$OUT/auto.refine.report.tsv" | wc -l | tr -d ' ') rows, expected 2"

# ---------------------------------------------------------------- a folded unit is held fixed
# panphorte's self-loop unit node carries copy number. refine re-aligns the residual flanks around it
# and must leave the unit run itself verbatim: same node, same count, same orientation.
UNIT=$(blk 64)
{ printf 'H\tVN:Z:1.0\n'
  printf "S\t1\t%s\nS\t2\t%s\nS\t3\t%s\nS\t4\t%s\nS\t5\t%s\n" "$BB1" "$UNIT" "$(blk 40)" "$(blk 38)" "$BB2"
  for e in "1 2" "2 2" "2 3" "2 4" "3 5" "4 5"; do
    set -- $e; printf "L\t%s\t+\t%s\t+\t0M\n" "$1" "$2"
  done
  printf 'P\tref\t1+,2+,2+,2+,3+,5+\t*\n'
  printf 'P\talt\t1+,2+,2+,2+,4+,5+\t*\n'
  printf 'P\tthird\t1+,2+,2+,4+,5+\t*\n'; } > "$OUT/rep.gfa"
"$BIN" bubble -i "$OUT/rep.gfa" -r ref -o "$OUT/repb" --min-variant-bp 0 --merge-nearby-bp 0 -q >/dev/null 2>&1
before_ref=$(spell "$OUT/repb.sorted.gfa" ref)
"$BIN" refine -i "$OUT/repb.sorted.gfa" --bubbles-csv "$OUT/repb.bubbles.csv" -r ref \
       -o "$OUT/repr" -q >/dev/null 2>&1
if [ -s "$OUT/repr.normalized.sorted.gfa" ]; then
  [ "$(spell "$OUT/repr.normalized.sorted.gfa" ref)" = "$before_ref" ] \
    && ok "a folded unit region preserves the reference's spelled sequence" \
    || bad "the unit region changed the reference's sequence"
  selfloop=$(awk '$1=="L" && $2==$4 {print $2}' "$OUT/repr.normalized.sorted.gfa" | head -1)
  [ -n "$selfloop" ] \
    && ok "the self-loop unit node survives refinement" \
    || bad "the self-loop edge was lost, so copy number is no longer recoverable"
  n=$(awk -v r="$selfloop" '$1=="P" && $2=="ref"{n=split($3,a,","); c=0
        for(i=1;i<=n;i++){s=a[i]; if(substr(s,1,length(s)-1)==r) c++} print c}' \
        "$OUT/repr.normalized.sorted.gfa")
  [ "$n" = "3" ] \
    && ok "the unit run keeps its copy count (3x)" \
    || bad "the unit is traversed ${n}x after refinement, expected 3"
else
  bad "the unit fixture produced no output"
fi

# ---------------------------------------------------------------- an unfolded non-unit revisit is skipped
# A plain interior node visited twice is an unfolded copy-number signal. POA would linearize it and
# destroy what `call` reconstructs, so the region must be skipped rather than rebuilt.
# Node 2 is revisited via 3 and carries NO self-loop, so panphorte never folded it: it is an unfolded
# copy-number signal, not a unit. (A self-loop would make it a unit and a different guard would fire.)
{ printf 'H\tVN:Z:1.0\n'
  printf "S\t1\t%s\nS\t2\t%s\nS\t3\t%s\nS\t4\t%s\nS\t5\t%s\n" \
         "$BB1" "$(blk 120)" "$(blk 30)" "$BB2" "$(blk 118)"
  for e in "1 2" "2 3" "3 2" "2 4" "1 5" "5 4"; do
    set -- $e; printf "L\t%s\t+\t%s\t+\t0M\n" "$1" "$2"
  done
  printf 'P\tref\t1+,2+,3+,2+,4+\t*\nP\talt\t1+,5+,4+\t*\n'; } > "$OUT/unf.gfa"
"$BIN" bubble -i "$OUT/unf.gfa" -r ref -o "$OUT/unfb" --min-variant-bp 0 --merge-nearby-bp 0 -q >/dev/null 2>&1
"$BIN" refine -i "$OUT/unfb.sorted.gfa" --bubbles-csv "$OUT/unfb.bubbles.csv" -r ref \
       -o "$OUT/unfr" > "$OUT/unfr.log" 2>&1
# report columns by name: <file> <column>, first data row
rrow() { awk -F'\t' -v w="$2" 'NR==1{for(i=1;i<=NF;i++)c[$i]=i;next} NR==2{print $(c[w])}' "$1" 2>/dev/null; }
[ "$(rrow "$OUT/unfr.refine.report.tsv" decision):$(rrow "$OUT/unfr.refine.report.tsv" reason)" = "skipped:unfolded_dup" ] \
  && ok "an unfolded non-unit revisit is skipped, with the reason recorded" \
  || bad "the unfolded duplication was not skipped: $(tail -1 "$OUT/unfr.refine.report.tsv" 2>/dev/null)"
# --min-dup-bp lets a revisit carrying fewer bases through: node 2 is 120 bp, so 200 realigns it
"$BIN" refine -i "$OUT/unfb.sorted.gfa" --bubbles-csv "$OUT/unfb.bubbles.csv" -r ref \
       -o "$OUT/unfm" --min-dup-bp 200 > "$OUT/unfm.log" 2>&1
[ "$(rrow "$OUT/unfm.refine.report.tsv" decision)" = "rebuilt" ] \
  && [ "$(rrow "$OUT/unfm.refine.report.tsv" short_revisits)" -ge 1 ] \
  && ok "a revisit shorter than --min-dup-bp is realigned, and the report says so" \
  || bad "--min-dup-bp 200 on a 120 bp revisit: $(tail -1 "$OUT/unfm.refine.report.tsv" 2>/dev/null)"

# ---------------------------------------------------------------- a haplotype crossing twice
# dup walks the region 1..4 twice, carrying both copies of the anchors. Each crossing is rebuilt on its
# own, so dup still walks both anchors twice and no old node survives: neither crossing is left on the old
# nodes, and the second copy of the anchors is not realigned into an insertion.
"${PYTHON:-python3}" - "$OUT" <<'PYEOF'
import random, sys
out = sys.argv[1]
r = random.Random(5)
rnd = lambda n: ''.join(r.choice('ACGT') for _ in range(n))
X, Y = rnd(40), rnd(40)
nodes = {'0': rnd(300), '1': rnd(60), '2': X + 'A' + Y, '3': X + 'C' + Y, '4': rnd(60), '9': rnd(300)}
paths = {'ref': '0+,1+,2+,4+,9+', 'h2': '0+,1+,3+,4+,9+', 'h3': '0+,1+,3+,4+,9+',
         'dup': '0+,1+,2+,4+,1+,3+,4+,9+'}
edges = set()
for w in paths.values():
    st = w.split(',')
    edges.update((a[:-1], a[-1], b[:-1], b[-1]) for a, b in zip(st, st[1:]))
with open(f'{out}/twice.gfa', 'w') as o:
    o.write('H\tVN:Z:1.0\n')
    for k, v in nodes.items(): o.write(f'S\t{k}\t{v}\n')
    for a, ao, b, bo in sorted(edges): o.write(f'L\t{a}\t{ao}\t{b}\t{bo}\t0M\n')
    for k, v in paths.items(): o.write(f'P\t{k}\t{v}\t*\n')
with open(f'{out}/twice_in.bubbles.csv', 'w') as o:
    o.write('bubble_id,source,source_orient,sink,sink_orient,inside_node_count,total_node_count,path_support,'
            'distinct_alleles,ref_allele_support,alt_allele_support_max,alt_allele_support_min,min_inside_bp,'
            'max_inside_bp,inside_nodes\n1,1,+,4,+,2,4,4,2,1,1,1,81,81,"2;3"\n')
PYEOF
"$BIN" refine -i "$OUT/twice.gfa" --bubbles-csv "$OUT/twice_in.bubbles.csv" -r ref -o "$OUT/twice" \
       > "$OUT/twice.log" 2>&1
# anchors 1 and 4 keep their sequences; count dup's visits to the nodes spelling them
tw_visits=$(awk -F'\t' -v a="$(awk '$1=="S" && $2=="1"{print $3}' "$OUT/twice.gfa")" \
                       -v b="$(awk '$1=="S" && $2=="4"{print $3}' "$OUT/twice.gfa")" '
  $1=="S" && ($3==a || $3==b) {anc[$2]=1}
  $1=="P" && $2=="dup" {n=split($3,s,","); for(i=1;i<=n;i++) if (substr(s[i],1,length(s[i])-1) in anc) c++; print c}' \
  "$OUT/twice.normalized.sorted.gfa" 2>/dev/null)
[ "$(rrow "$OUT/twice.refine.report.tsv" decision)" = "rebuilt" ] && [ "$tw_visits" = "4" ] \
  && [ "$(rrow "$OUT/twice.refine.report.tsv" crossings)" = "5" ] \
  && ok "both crossings of a twice-crossing haplotype are rebuilt and its anchors stay duplicated" \
  || bad "twice-crossing haplotype: $(tail -1 "$OUT/twice.refine.report.tsv" 2>/dev/null), anchor visits '$tw_visits' (expected 4)"

# ---------------------------------------------------------------- partial traversal
# A path entering the interior without spanning both anchors cannot be rewritten. Retaining its nodes
# also retains the old edges, so refined and unrefined topology coexist -- the region is skipped.
{ printf 'H\tVN:Z:1.0\n'
  printf "S\t1\t%s\nS\t2\t%s\nS\t3\t%s\nS\t4\t%s\nS\t5\t%s\n" \
         "$BB1" "$(blk 120)" "$(blk 118)" "$BB2" "$(blk 90)"
  for e in "1 2" "1 3" "2 4" "3 4" "5 2"; do set -- $e; printf "L\t%s\t+\t%s\t+\t0M\n" "$1" "$2"; done
  printf 'P\tref\t1+,2+,4+\t*\nP\talt\t1+,3+,4+\t*\nP\tpartial\t5+,2+,4+\t*\n'; } > "$OUT/part.gfa"
"$BIN" bubble -i "$OUT/part.gfa" -r ref -o "$OUT/partb" --min-variant-bp 0 --merge-nearby-bp 0 -q >/dev/null 2>&1
"$BIN" refine -i "$OUT/partb.sorted.gfa" --bubbles-csv "$OUT/partb.bubbles.csv" -r ref \
       -o "$OUT/ps" > "$OUT/ps.log" 2>&1
[ "$(rebuilt "$OUT/ps.log")" = "0 1" ] \
  && ok "a region with a partial traverser is skipped" \
  || bad "partial traversal gave '$(rebuilt "$OUT/ps.log")', expected '0 1'"

# ---------------------------------------------------------------- W-line graphs
# parse_gfa names a W path sample#hap#seqid:start-end, the 0-based W start written 1-based. refine must
# name it the same way, or the resolved reference matches no path.
W="$(cd "$(dirname "$0")/.." && pwd)/synthetic_data/syn_w.gfa"
if [ -f "$W" ]; then
  wref=$(awk '$1=="W"{print $2"#"$3"#"$4":"$5+1"-"$6; exit}' "$W")
  "$BIN" bubble -i "$W" -r "$wref" -o "$OUT/wb" --min-variant-bp 0 --merge-nearby-bp 0 -q >/dev/null 2>&1
  "$BIN" refine -i "$OUT/wb.sorted.gfa" --bubbles-csv "$OUT/wb.bubbles.csv" -r "$wref" \
         -o "$OUT/wr" -q >/dev/null 2>&1
  [ "$?" -eq 0 ] && ok "a W-line graph resolves its reference and refines" \
                 || bad "refine failed on a W-line graph"
else
  ok "W-line fixture not present, skipped"
fi

# ---------------------------------------------------------------- duplicate carriers must not decide
# abPOA is handed DISTINCT sequences, so adding an identical haplotype changes no POA input, and the
# size limits must not change either. The work over the distinct sequences is 120 x 238 = 28,560; summed
# over the carriers it would be 120 x 1,182 = 141,840 in dup2, so a budget of 100,000 tells them apart.
{ printf 'H\tVN:Z:1.0\n'
  printf "S\t1\t%s\nS\t2\t%s\nS\t3\t%s\nS\t4\t%s\n" "$BB1" "$(blk 120)" "$(blk 118)" "$BB2"
  for e in "1 2" "1 3" "2 4" "3 4"; do set -- $e; printf "L\t%s\t+\t%s\t+\t0M\n" "$1" "$2"; done
  printf 'P\tref\t1+,2+,4+\t*\nP\talt\t1+,3+,4+\t*\n'; } > "$OUT/dup1.gfa"
{ cat "$OUT/dup1.gfa"
  for i in 1 2 3 4 5 6 7 8; do printf "P\tcopy%s\t1+,3+,4+\t*\n" "$i"; done; } > "$OUT/dup2.gfa"
for g in dup1 dup2; do
  "$BIN" bubble -i "$OUT/$g.gfa" -r ref -o "$OUT/${g}b" --min-variant-bp 0 --merge-nearby-bp 0 -q >/dev/null 2>&1
  "$BIN" refine -i "$OUT/${g}b.sorted.gfa" --bubbles-csv "$OUT/${g}b.bubbles.csv" -r ref \
         -o "$OUT/${g}r" --max-poa-work 100000 > "$OUT/${g}r.log" 2>&1
done
[ "$(rebuilt "$OUT/dup1r.log")" = "1 0" ] && [ "$(rebuilt "$OUT/dup2r.log")" = "1 0" ] \
  && ok "replicating identical carriers does not change the POA decision" \
  || bad "decision changed with duplicate carriers: '$(rebuilt "$OUT/dup1r.log")' vs '$(rebuilt "$OUT/dup2r.log")'"

# ---------------------------------------------------------------- shuffled bubble rows
# Component order must come from the graph, not from CSV row order.
head -1 "$OUT/disb.bubbles.csv" > "$OUT/shuf.csv"
tail -n +2 "$OUT/disb.bubbles.csv" | tail -r >> "$OUT/shuf.csv" 2>/dev/null || \
  { tail -n +2 "$OUT/disb.bubbles.csv" | tac >> "$OUT/shuf.csv"; }
"$BIN" refine -i "$OUT/disb.sorted.gfa" --bubbles-csv "$OUT/shuf.csv" -r ref -o "$OUT/sh" \
       -q >/dev/null 2>&1
cmp -s "$OUT/d1.normalized.sorted.gfa" "$OUT/sh.normalized.sorted.gfa" \
  && ok "reordering the bubbles CSV gives byte-identical output" \
  || bad "output depends on the order of rows in the bubbles CSV"
# The CSV is read sorted by id, so reordering rows alone hands refine the same list; swapping the two
# sites' ids reverses that list.
{ head -1 "$OUT/disb.bubbles.csv"
  tail -n +2 "$OUT/disb.bubbles.csv" | awk -F, -v OFS=, '{$1 = 3 - $1; print}'; } > "$OUT/swap.csv"
"$BIN" refine -i "$OUT/disb.sorted.gfa" --bubbles-csv "$OUT/swap.csv" -r ref -o "$OUT/sw" \
       -q >/dev/null 2>&1
cmp -s "$OUT/d1.normalized.sorted.gfa" "$OUT/sw.normalized.sorted.gfa" \
  && cmp -s "$OUT/d1.refine.report.tsv" "$OUT/sw.refine.report.tsv" \
  && ok "swapping the sites' ids gives byte-identical output and report" \
  || bad "output depends on the sites' ids"

# ---------------------------------------------------------------- the re-snarl threshold is usable
"$BIN" refine -i "$OUT/disb.sorted.gfa" --bubbles-csv "$OUT/disb.bubbles.csv" -r ref \
       -o "$OUT/rs" --resnarl-min-variant-bp 0 -q >/dev/null 2>&1
[ "$?" -eq 0 ] && [ -s "$OUT/rs.bubbles.csv" ] \
  && ok "--resnarl-min-variant-bp is accepted and produces a CSV" \
  || bad "--resnarl-min-variant-bp is advertised but not parsed"

# ---------------------------------------------------------------- an output may not name the GTF
printf '#gtf\n' > "$OUT/genes.gtf"
"$BIN" refine -i "$OUT/disb.sorted.gfa" --bubbles-csv "$OUT/disb.bubbles.csv" -r ref \
       -o "$OUT/genes" --gtf "$OUT/genes.gtf" -q >/dev/null 2>&1
# <prefix>.bandage_genes.csv vs --gtf genes.gtf do not collide; the direct collision is the check below
cp "$OUT/disb.bubbles.csv" "$OUT/alias.bandage_genes.csv"
"$BIN" refine -i "$OUT/disb.sorted.gfa" --bubbles-csv "$OUT/disb.bubbles.csv" -r ref \
       -o "$OUT/alias" --gtf "$OUT/alias.bandage_genes.csv" -q >/dev/null 2>&1
[ "$?" -ne 0 ] && ok "an output that would overwrite the --gtf input is refused" \
               || bad "the --gtf input was left open to being overwritten"

# ---------------------------------------------------------------- twin merge
# A twin is the same sequence in the same place on a second node -- here the graph builder put the
# forward contigs' copy of a 22 bp element on node 3 and the reverse contig's on its looped twin 4. Node
# identity then reads assembly strand as a difference. refine merges twins
# without realigning anything; the controls must stay: node 6 has the same sequence and predecessor
# but another successor, and 8 -> 9 is a homopolymer run (T, then T again), not two copies.
TW=$(blk 22)
{ printf 'H\tVN:Z:1.0\n'
  printf "S\t1\t%s\nS\t2\t%s\nS\t3\t%s\nS\t4\t%s\nS\t5\t%s\nS\t6\t%s\nS\t7\t%s\n" \
         "$BB1" "$(blk 40)" "$TW" "$TW" "$BB2" "$TW" "$(blk 37)"
  printf "S\t8\tT\nS\t9\tT\nS\t10\t%s\n" "$BB3"
  for e in "1 + 2 +" "2 + 3 +" "3 + 5 +" "2 + 4 +" "4 + 4 +" "4 + 5 +" "2 + 6 +" "6 + 7 +" "7 + 5 +" \
           "5 + 8 +" "8 + 9 +" "8 + 10 +" "9 + 10 +"; do
    set -- $e; printf "L\t%s\t%s\t%s\t%s\t0M\n" "$1" "$2" "$3" "$4"
  done
  printf 'P\tref\t1+,2+,3+,5+,8+,10+\t*\n'
  printf 'P\tfwd\t1+,2+,3+,5+,8+,9+,10+\t*\n'
  printf 'P\trev\t10-,8-,5-,4-,4-,2-,1-\t*\n'
  printf 'P\tctl\t1+,2+,6+,7+,5+,8+,10+\t*\n'; } > "$OUT/twin.gfa"
"$BIN" bubble -i "$OUT/twin.gfa" -r ref -o "$OUT/twinb" --min-variant-bp 0 --merge-nearby-bp 0 -q >/dev/null 2>&1
# No bubble is listed, so the POA rebuild never runs and only the twin merge can change the graph:
# what is tested here is the merge rule itself.
head -1 "$OUT/twinb.bubbles.csv" > "$OUT/twinb.none.csv"
spell_paths() { awk -F'\t' '$1=="S"{s[$2]=$3} $1=="P"{p[$2]=$3} END{for(n in p){m=split(p[n],st,","); o="";
  for(i=1;i<=m;i++){id=substr(st[i],1,length(st[i])-1); q=s[id]; if(substr(st[i],length(st[i]))=="-"){r="";
  for(j=length(q);j>0;j--){c=substr(q,j,1); r=r (c=="A"?"T":c=="C"?"G":c=="G"?"C":"A")} q=r} o=o q} print n"\t"o}}' "$1" | sort; }
"$BIN" refine -i "$OUT/twinb.sorted.gfa" --bubbles-csv "$OUT/twinb.none.csv" -r ref -o "$OUT/twinr" \
  > "$OUT/twinr.log" 2>&1 || bad "refine failed on the twin fixture: $(tail -1 "$OUT/twinr.log")"
seq_id() { awk -F'\t' -v s="$2" '$1=="S" && $3==s{print $2}' "$1" | sort | tr '\n' ' '; }
twins_now=$(seq_id "$OUT/twinr.normalized.sorted.gfa" "$TW")
[ "$(echo $twins_now | wc -w | tr -d ' ')" = "2" ] \
  && ok "the strand-split twin is merged; the twin with another successor stays a separate node" \
  || bad "nodes spelling the 22 bp element after refine: '$twins_now' (expected the merged twin and the control)"
grep -q "twin nodes merged: 1" "$OUT/twinr.log" \
  && ok "the log counts the merge" || bad "the log does not count one twin merge: $(grep twin "$OUT/twinr.log")"
# The node counts add up: input nodes + added - removed - merged = output nodes.
awk '/input: .*nodes: /{sub(/.*nodes: /, ""); split($0, a, ";"); n0 = a[1]}
     /\] nodes: /{sub(/.*\] nodes: /, ""); gsub(/[a-z ]+: /, ""); split($0, a, ";"); n1 = a[1]; d = a[2] - a[3] - a[4]}
     END{exit !(n0 + d == n1 && n0 != "")}' "$OUT/twinr.log" \
  && ok "the node-count line adds up with the merged twins" \
  || bad "node counts do not add up: $(grep -E 'input:|\] nodes:' "$OUT/twinr.log" | tr '\n' ' ')"
[ ! -e "$OUT/twinr.refine.twins.tsv" ] && ok "no twin table is written" || bad "a twin table was written"
awk -F'\t' '$1=="L" && $2==$4 && $3==$5' "$OUT/twinr.normalized.sorted.gfa" | grep -q . \
  && ok "the merged node keeps the self-loop" || bad "the merged twin lost its self-loop"
[ "$(awk -F'\t' '$1=="S" && $3=="T"' "$OUT/twinr.normalized.sorted.gfa" | wc -l | tr -d ' ')" = "2" ] \
  && ok "a homopolymer run (T then T) is not merged into one node" \
  || bad "the homopolymer run was collapsed"
cmp -s <(spell_paths "$OUT/twin.gfa") <(spell_paths "$OUT/twinr.normalized.sorted.gfa") \
  && ok "every path spells exactly what it did before the merge" \
  || bad "a path's sequence changed across the twin merge"
# A looped twin whose loop no path walks must stay separate: call reads every self-linked node as a
# repeat unit, so merging it into the reference's node would hand every haplotype a loop it never took.
awk -F'\t' 'BEGIN{OFS="\t"} $1=="P" && $2=="rev"{n=split($3,s,","); o=""; for(i=1;i<=n;i++){ if(i>1 && s[i]==s[i-1]) continue; o=o (o==""?"":",") s[i]} $3=o} {print}' \
  "$OUT/twinb.sorted.gfa" > "$OUT/twinu.gfa"
"$BIN" refine -i "$OUT/twinu.gfa" --bubbles-csv "$OUT/twinb.none.csv" -r ref -o "$OUT/twinu" > "$OUT/twinu.log" 2>&1 \
  || bad "refine failed on the unwalked-loop twin fixture: $(tail -1 "$OUT/twinu.log")"
grep -q "twin nodes merged: 0" "$OUT/twinu.log" \
  && ok "a looped twin whose loop no path walks is not merged" \
  || bad "an unwalked looped twin was merged: $(grep 'twin nodes' "$OUT/twinu.log")"
# Twins with different KINDS of self-link -- a tandem loop on one, a fold-back on the other, both
# walked -- must stay apart: the merged node would carry both, and a walk through tandem then
# fold-back would become possible that no haplotype has.
{ printf 'H\tVN:Z:1.0\n'
  printf "S\t1\t%s\nS\t2\tGATTACAGG\nS\t3\tGATTACAGG\nS\t4\t%s\n" "$(blk 60)" "$(blk 61)"
  for e in "1 + 2 +" "1 + 3 +" "2 + 4 +" "3 + 4 +" "2 + 2 +" "3 + 3 -"; do
    set -- $e; printf "L\t%s\t%s\t%s\t%s\t0M\n" "$1" "$2" "$3" "$4"; done
  printf 'P\tref\t1+,2+,4+\t*\nP\th1\t1+,2+,2+,4+\t*\nP\th2\t1+,3+,3-,1-\t*\nP\th3\t1+,3+,4+\t*\n'; } > "$OUT/twinmix.gfa"
"$BIN" refine -i "$OUT/twinmix.gfa" --bubbles-csv "$OUT/twinb.none.csv" -r ref -o "$OUT/twinmix" > "$OUT/twinmix.log" 2>&1 \
  || bad "refine failed on the mixed self-link twin fixture: $(tail -1 "$OUT/twinmix.log")"
grep -q "twin nodes merged: 0" "$OUT/twinmix.log" \
  && ok "twins with different kinds of self-link (tandem vs fold-back) are not merged" \
  || bad "twins with a tandem loop and a fold-back were merged into one node: $(grep 'twin nodes' "$OUT/twinmix.log")"

# ---------------------------------------------------------------- the re-snarl fuses by placement
# A 72 bp (AC)n microsatellite holding a 2 bp deletion (under the size floor) and an expansion. With POA
# kept off the region (--max-poa-bp 1), refine only re-snarls it: the output CSV keeps the site bubble
# fused over the repeat (2..7), and --no-fuse-ambiguous re-snarls it to the expansion's own boundaries.
ac() { local o="" i; for ((i = 0; i < $1; i++)); do o="${o}AC"; done; printf '%s' "$o"; }
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\t1\t%s\nS\t9\t%s\n' "$(blk 300)" "$(blk 300 | rev)"
  printf 'S\t2\t%s\nS\t3\tAC\nS\t5\t%s\nS\t6\t%s\nS\t7\t%s\nS\t8\t%s\n' "$(ac 10)" "$(ac 10)" "$(ac 5)" "$(ac 40)" "$(ac 10)"
  for e in "1 2" "2 3" "3 5" "2 5" "5 6" "5 7" "6 8" "7 8" "8 9"; do set -- $e; printf 'L\t%s\t+\t%s\t+\t0M\n' "$1" "$2"; done
  printf 'P\tref\t1+,2+,3+,5+,6+,8+,9+\t*\nP\thdel\t1+,2+,5+,6+,8+,9+\t*\nP\thexp\t1+,2+,3+,5+,7+,8+,9+\t*\n'; } > "$OUT/tr.gfa"
"$BIN" bubble -i "$OUT/tr.gfa" -r ref -o "$OUT/trb" --merge-nearby-bp 0 -q >/dev/null 2>&1
"$BIN" refine -i "$OUT/trb.sorted.gfa" --bubble-prefix "$OUT/trb" -r ref -o "$OUT/trr" --max-poa-bp 1 -q >/dev/null 2>&1
"$BIN" refine -i "$OUT/trb.sorted.gfa" --bubble-prefix "$OUT/trb" -r ref -o "$OUT/trrn" --max-poa-bp 1 --no-fuse-ambiguous -q >/dev/null 2>&1
csv_sites() { tail -n +2 "$1" 2>/dev/null | cut -d, -f2,4 | paste -sd' ' -; }
[ "$(csv_sites "$OUT/trr.bubbles.csv"):$(csv_sites "$OUT/trrn.bubbles.csv")" = "2,7:4,7" ] \
  && ok "the re-snarl fuses a repeat's sites by placement (2,7), and not with --no-fuse-ambiguous (4,7)" \
  || bad "refine re-snarl: $(csv_sites "$OUT/trr.bubbles.csv") / $(csv_sites "$OUT/trrn.bubbles.csv"), expected 2,7 / 4,7"


# ---------------------------------------------------------------- folded repeats from panphorte's tables
# A unit graph of several nodes has no self-loop, nor does a one-node one no haplotype walks twice in a row;
# read by its self-loop alone, either would be realigned as plain sequence. refine reads
# <prefix>.panphorte.unit_graph.tsv, holds those nodes fixed and writes the table forward with its own ids.
PC="$(cd "$(dirname "$0")/.." && pwd)/synthetic_data/panphorte_cases"
walk() { awk -F'\t' -v p="$2" '$1=="P" && $2==p{print $3}' "$1"; }
ugcol() { awk -F'\t' -v w="$2" 'NR==1{for(i=1;i<=NF;i++)c[$i]=i;next} {print $(c[w])}' "$1" 2>/dev/null; }
# two phase-rotated one-node folds of one site, no self-loop on the second: hX walks two copies
"$BIN" bubble -i "$PC/refine_rep_no_loop.gfa" -r ref -o "$OUT/rnb" --min-variant-bp 0 --merge-nearby-bp 0 -q >/dev/null 2>&1
"$BIN" panphorte -i "$OUT/rnb.sorted.gfa" -b "$OUT/rnb" -o "$OUT/rnp" -r ref --threads 1 -q >/dev/null 2>&1
"$BIN" refine -i "$OUT/rnp.normalized.sorted.gfa" -b "$OUT/rnp" -r ref -o "$OUT/rnr" -q >/dev/null 2>&1
[ -n "$(walk "$OUT/rnp.normalized.sorted.gfa" hX)" ] \
  && [ "$(walk "$OUT/rnr.normalized.sorted.gfa" hX)" = "$(walk "$OUT/rnp.normalized.sorted.gfa" hX)" ] \
  && [ "$(ugcol "$OUT/rnr.panphorte.unit_graph.tsv" nodes | paste -sd' ' -)" = "5 7" ] \
  && ok "a one-node fold without a self-loop is held fixed (hX keeps both copies) and its table carried forward" \
  || bad "one-node fold without a self-loop: hX $(walk "$OUT/rnp.normalized.sorted.gfa" hX) -> $(walk "$OUT/rnr.normalized.sorted.gfa" hX), table '$(ugcol "$OUT/rnr.panphorte.unit_graph.tsv" nodes | paste -sd' ' -)'"
# a two-node unit graph nobody revisits (the copies differ only in letter case)
"$BIN" bubble -i "$PC/case_period.gfa" -r arr -o "$OUT/cpb" --min-variant-bp 0 --merge-nearby-bp 0 -q >/dev/null 2>&1
"$BIN" panphorte -i "$OUT/cpb.sorted.gfa" -b "$OUT/cpb" -o "$OUT/cpp" -r arr --min-similarity 0.95 --threads 1 -q >/dev/null 2>&1
"$BIN" refine -i "$OUT/cpp.normalized.sorted.gfa" -b "$OUT/cpp" -r arr -o "$OUT/cpr" -q >/dev/null 2>&1
[ "$(walk "$OUT/cpr.normalized.sorted.gfa" arr)" = "1+,2+,3+,4+" ] \
  && [ "$(cut -f7,9 "$OUT/cpr.panphorte.unit_graph.tsv" 2>/dev/null | tail -1)" = "$(printf '2;3\t2;3')" ] \
  && ok "a unit graph without a self-loop is held fixed, its table carried forward" \
  || bad "unit graph without a self-loop: arr walks $(walk "$OUT/cpr.normalized.sorted.gfa" arr), table '$(tail -1 "$OUT/cpr.panphorte.unit_graph.tsv" 2>/dev/null)'"
# A unit graph walked 1-3 times, beside a flank the builder split (h2/h4 carry a 60 bp deletion as two
# nodes): the flank is realigned around the repeat block, every haplotype keeps its copies -- also with a
# --min-dup-bp above the copies' bases, which must not realign them as a plain duplication.
UU=$(blk 300 | rev); UX=$(blk 150 | tr ACGT TGCA)
snv() { local c=${1:$2:1}; [ "$c" = C ] && c=G || c=C; printf '%s' "${1:0:$2}$c${1:$(($2 + 1))}"; }
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\t1\t%s\nS\t2\t%s\nS\t3\t%s\nS\t5\t%s\nS\t4\t%s\nS\t6\t%s\nS\t7\t%s\nS\t8\t%s\n' "$BB1" "$UU" "$(snv "$UU" 30)" \
    "$(snv "$UU" 40)" "$BB2" "$UX" "${UX:0:40}" "${UX:100}"
  for e in "1 6" "6 2" "1 7" "7 8" "8 2" "2 2" "2 3" "3 4" "2 4" "2 5" "5 2"; do set -- $e; printf 'L\t%s\t+\t%s\t+\t0M\n' "$1" "$2"; done
  printf 'P\tref\t1+,6+,2+,3+,4+\t*\nP\th1\t1+,6+,2+,2+,4+\t*\nP\th2\t1+,7+,8+,2+,5+,2+,4+\t*\n'
  printf 'P\th3\t1+,6+,2+,4+\t*\nP\th4\t1+,7+,8+,2+,2+,2+,4+\t*\n'; } > "$OUT/ugf.gfa"
"$BIN" bubble -i "$OUT/ugf.gfa" -r ref -o "$OUT/ugfb" --min-variant-bp 0 --merge-nearby-bp 0 -q >/dev/null 2>&1
"$BIN" panphorte -i "$OUT/ugfb.sorted.gfa" -b "$OUT/ugfb" -o "$OUT/ugfp" -r ref --min-similarity 0.95 --threads 1 -q >/dev/null 2>&1
"$BIN" refine -i "$OUT/ugfp.normalized.sorted.gfa" -b "$OUT/ugfp" -r ref -o "$OUT/ugfr" -q >/dev/null 2>&1
"$BIN" refine -i "$OUT/ugfp.normalized.sorted.gfa" -b "$OUT/ugfp" -r ref -o "$OUT/ugfd" --min-dup-bp 100000 -q >/dev/null 2>&1
passes() { awk -F'\t' 'FNR==1 && FNR==NR{for(i=1;i<=NF;i++)col[$i]=i} FNR==NR{v=$(col["anchor_nodes"]); if(FNR>1 && v!="."){k=split(v,x,";"); for(j=1;j<=k;j++) an[x[j]]=1}; next}
                       $1=="P"{n=split($3,a,","); c=0; for(i=1;i<=n;i++) if (substr(a[i],1,length(a[i])-1) in an) c++; printf "%s %d ", $2, c}' \
                  "$1.panphorte.unit_graph.tsv" "$1.normalized.sorted.gfa" 2>/dev/null; }
for r in ugfr ugfd; do
  [ "$(passes "$OUT/ugfp")" = "ref 2 h1 2 h2 3 h3 1 h4 3 " ] && [ "$(passes "$OUT/$r")" = "$(passes "$OUT/ugfp")" ] \
    && [ "$(awk -F'\t' 'NR==2{print $5, $6}' "$OUT/$r.refine.report.tsv")" = "rebuilt unit_residual" ] \
    && ok "unit graph beside a split flank ($r): flank rebuilt around the repeat, copies 2/2/3/1/3 kept" \
    || bad "unit graph beside a split flank ($r): $(awk -F'\t' 'NR==2{print $5, $6, $7}' "$OUT/$r.refine.report.tsv"), passes '$(passes "$OUT/$r")' (panphorte '$(passes "$OUT/ugfp")')"
done

# ---------------------------------------------------------------- letters POA does not keep
# abPOA reads soft-masked bases as upper case and anything else as N; decoded from its codes, the rows
# would change the haplotypes. Each row keeps its own bytes.
csvhdr='bubble_id,source,source_orient,sink,sink_orient,inside_node_count,total_node_count,path_support,distinct_alleles,ref_allele_support,alt_allele_support_max,alt_allele_support_min,min_inside_bp,max_inside_bp,inside_nodes'
{ printf 'H\tVN:Z:1.0\nS\t1\tGG\nS\t2\tacgt\nS\t3\tacgR\nS\t4\tTT\n'
  for e in "1 2" "2 4" "1 3" "3 4"; do set -- $e; printf 'L\t%s\t+\t%s\t+\t0M\n' "$1" "$2"; done
  printf 'P\tref\t1+,2+,4+\t*\nP\talt\t1+,3+,4+\t*\n'; } > "$OUT/lc.gfa"
printf '%s\n1,1,+,4,+,2,4,2,2,1,1,1,4,4,"2;3"\n' "$csvhdr" > "$OUT/lc.bubbles.csv"
"$BIN" refine -i "$OUT/lc.gfa" --bubbles-csv "$OUT/lc.bubbles.csv" -r ref -o "$OUT/lcr" -q > "$OUT/lcr.log" 2>&1
[ "$(awk -F'\t' 'NR==2{print $5}' "$OUT/lcr.refine.report.tsv" 2>/dev/null)" = "rebuilt" ] \
  && [ "$(grep -c '^S.*acg[tR]*$' "$OUT/lcr.normalized.sorted.gfa")" -ge 1 ] && grep -q $'^S\t[0-9]*\tR$' "$OUT/lcr.normalized.sorted.gfa" \
  && ok "lower-case and IUPAC bases come through POA unchanged" \
  || bad "lower-case/IUPAC site: $(head -c 200 "$OUT/lcr.log")"

# ---------------------------------------------------------------- one L line per link
# A link stored one way and walked the other (L 2 - 1 -, path 1+,2+) is written once.
{ printf 'H\tVN:Z:1.0\nS\t1\t%s\nS\t2\t%s\nL\t2\t-\t1\t-\t0M\nP\tref\t1+,2+\t*\nP\thap\t1+,2+\t*\n' "$(blk 40)" "$(blk 41)"; } > "$OUT/rl.gfa"
printf '%s\n' "$csvhdr" > "$OUT/rl.bubbles.csv"
"$BIN" refine -i "$OUT/rl.gfa" --bubbles-csv "$OUT/rl.bubbles.csv" -r ref -o "$OUT/rlr" -q >/dev/null 2>&1
[ "$(grep -c '^L' "$OUT/rlr.normalized.sorted.gfa" 2>/dev/null)" = "1" ] \
  && ok "a link walked against its stored direction is not written twice" \
  || bad "reverse-stored link: $(grep -c '^L' "$OUT/rlr.normalized.sorted.gfa" 2>/dev/null) L lines"

# ---------------------------------------------------------------- overlapping sites, --bubble-id, stale files
# A site inside another (a CSV not from bubble) would be rewritten twice over the same steps. Both are
# skipped with a reason.
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\t1\t%s\nS\t2\t%s\nS\t3\t%s\nS\t4\t%s\nS\t5\t%s\nS\t6\t%s\nS\t7\t%s\n' "$BB1" "$(blk 56)" "$A1" "$A2" "$(blk 57)" "$BB2" "$A3"
  for e in "1 2" "2 3" "2 4" "3 5" "4 5" "5 6" "1 7" "7 6"; do set -- $e; printf 'L\t%s\t+\t%s\t+\t0M\n' "$1" "$2"; done
  printf 'P\tref\t1+,2+,3+,5+,6+\t*\nP\th1\t1+,2+,4+,5+,6+\t*\nP\th2\t1+,7+,6+\t*\n'; } > "$OUT/nest.gfa"
printf '%s\n1,1,+,6,+,5,7,3,3,1,1,1,100,300,"2;3;4;5;7"\n2,2,+,5,+,2,4,2,2,1,1,1,117,120,"3;4"\n' "$csvhdr" > "$OUT/nest.bubbles.csv"
"$BIN" refine -i "$OUT/nest.gfa" --bubbles-csv "$OUT/nest.bubbles.csv" -r ref -o "$OUT/nestr" -q > "$OUT/nestr.log" 2>&1
[ "$(cut -f6 "$OUT/nestr.refine.report.tsv" 2>/dev/null | tail -n +2 | paste -sd' ' -)" = "overlapping_sites overlapping_sites" ] \
  && ok "overlapping sites are both skipped (overlapping_sites)" \
  || bad "overlapping sites: $(cut -f5-7 "$OUT/nestr.refine.report.tsv" 2>/dev/null | tail -n +2 | paste -sd'|' -) $(head -c 200 "$OUT/nestr.log")"
"$BIN" refine -i "$OUT/disb.sorted.gfa" --bubbles-csv "$OUT/disb.bubbles.csv" -r ref \
       -o "$OUT/rep2" --bubble-id 1 --bubble-id 2 > "$OUT/rep2.log" 2>&1
[ "$(rebuilt "$OUT/rep2.log")" = "2 0" ] \
  && ok "--bubble-id is repeatable" || bad "--bubble-id 1 --bubble-id 2 gave '$(rebuilt "$OUT/rep2.log")', expected '2 0'"
empty_ok=1
for v in "" "," "1," "x"; do
  "$BIN" refine -i "$OUT/disb.sorted.gfa" --bubbles-csv "$OUT/disb.bubbles.csv" -r ref -o "$OUT/emp" \
         --bubble-id "$v" -q >/dev/null 2>&1 && empty_ok=0
done
[ "$empty_ok" = 1 ] && ok "an empty or malformed --bubble-id is refused, not read as every site" \
  || bad "an empty or malformed --bubble-id was accepted"
echo stale > "$OUT/stl.bandage_genes.csv"; echo stale > "$OUT/stl.panphorte.unit_graph.tsv"
"$BIN" refine -i "$OUT/disb.sorted.gfa" --bubbles-csv "$OUT/disb.bubbles.csv" -r ref -o "$OUT/stl" -q >/dev/null 2>&1
[ ! -e "$OUT/stl.bandage_genes.csv" ] && [ ! -e "$OUT/stl.panphorte.unit_graph.tsv" ] && [ -s "$OUT/stl.normalized.sorted.gfa" ] \
  && ok "gene and repeat tables of an earlier run at the prefix are removed" \
  || bad "stale optional outputs survived a run that did not write them"
# ---- fixtures in synthetic_data/refine_cases ------------------------------------------------------
RC="$(cd "$(dirname "$0")/.." && pwd)/synthetic_data/refine_cases"
# A table is checked against the graph, not only for its node ids: one from another run can name nodes of
# this graph (sorted graphs reuse small ids), and must not be carried forward with this graph's sequence.
"$BIN" refine -i "$RC/table_conflict.gfa" --bubbles-csv "$RC/table_conflict.bubbles.csv" \
  --panphorte-prefix "$RC/stale_phase" -r ref -o "$OUT/r3fo" -q >/dev/null 2>"$OUT/r3fo.err"; rc=$?
[ "$rc" -ne 0 ] && grep -q "does not describe this graph" "$OUT/r3fo.err" && [ ! -e "$OUT/r3fo.normalized.sorted.gfa" ] \
  && ok "a unit graph table whose one-node unit is not its node's sequence is refused" \
  || bad "foreign unit graph table: exit $rc, $(head -c 120 "$OUT/r3fo.err")"
"$BIN" refine -i "$RC/table_conflict.gfa" -b "$RC/table_conflict" -r ref -o "$OUT/r3tc" -q >/dev/null 2>"$OUT/r3tc.err"; rc=$?
[ "$rc" -ne 0 ] && grep -q "is in site 1 and in site 999" "$OUT/r3tc.err" \
  && ok "a table putting one node in two sites is refused" \
  || bad "contradictory unit graph table: exit $rc, $(head -c 120 "$OUT/r3tc.err")"
# -b A beside --bubbles-csv B would read B's sites with A's tables.
"$BIN" refine -i "$RC/table_conflict.gfa" -b "$RC/table_conflict" --bubbles-csv "$RC/unit_only.bubbles.csv" \
  -r ref -o "$OUT/r3mx" -q >/dev/null 2>"$OUT/r3mx.err"; rc=$?
[ "$rc" -ne 0 ] && grep -q "name different runs" "$OUT/r3mx.err" \
  && ok "-b and --bubbles-csv naming different runs is refused" \
  || bad "mixed families: exit $rc, $(head -c 120 "$OUT/r3mx.err")"
# Two rows whose listed nodes are disjoint still overlap when one's crossings walk the other's nodes
# (alt_outer walks 3, rows list 2;7 and 4;8): both are skipped, and the run stays lossless.
"$BIN" refine -i "$RC/hidden_overlap.gfa" --bubbles-csv "$RC/hidden_overlap.bubbles.csv" -r ref -o "$OUT/r3ho" \
  -q >/dev/null 2>"$OUT/r3ho.err"; rc=$?
[ "$rc" -eq 0 ] && [ "$(awk -F'\t' 'NR>1{print $6}' "$OUT/r3ho.refine.report.tsv" | sort -u)" = "overlapping_sites" ] \
  && ok "rows overlapping through walked (unlisted) nodes are both skipped" \
  || bad "hidden overlap: exit $rc, reasons $(awk -F'\t' 'NR>1{print $6}' "$OUT/r3ho.refine.report.tsv" 2>/dev/null | paste -sd, -)"
# A link the input writes in both readings is one link.
"$BIN" refine -i "$RC/dup_links.gfa" --bubbles-csv "$RC/empty.bubbles.csv" -r ref -o "$OUT/r3dl" -q >/dev/null 2>&1
[ "$(awk '$1=="L"' "$OUT/r3dl.normalized.sorted.gfa" | wc -l | tr -d ' ')" = "1" ] \
  && ok "an input link written in both readings is emitted once" \
  || bad "duplicate input links: $(awk '$1=="L"' "$OUT/r3dl.normalized.sorted.gfa" | wc -l | tr -d ' ') L lines"

# ---------------------------------------------------------------- spec fixtures (tests/tools/spec_gfa.py)
SPEC_GFA="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/tools/spec_gfa.py"
spec() { mkdir -p "$OUT/$1"; cat > "$OUT/$1/spec"; python3 "$SPEC_GFA" "$OUT/$1/spec" "$OUT/$1"; }
# The CSV lists only node 2, but alt's crossing walks 3 and frag walks 3 outside any crossing. Rebuilding
# would remove 3 from alt and keep it, with its old links, for frag: the old route beside the new one.
spec unlisted <<'EOF2'
N 1 300 1
N 2 120 2
N 3 118 3
N 4 300 4
N 5 50 5
P ref 1+,2+,4+
P alt 1+,3+,4+
P frag 5+,3+,4+
B 1 4 2
EOF2
"$BIN" refine -i "$OUT/unlisted/g.gfa" -b "$OUT/unlisted/g" -r ref -o "$OUT/unl" -q >/dev/null 2>&1; rc=$?
[ "$rc" -eq 0 ] && [ "$(rrow "$OUT/unl.refine.report.tsv" reason)" = "partial_paths" ] \
  && ok "a path walking a crossing's unlisted node outside any crossing skips the region" \
  || bad "unlisted crossing node: exit $rc, reason '$(rrow "$OUT/unl.refine.report.tsv" reason)'"
# A chain whose site ids run against the reference: a skipped region still reports its outer anchors.
spec revchain <<'EOF2'
N 1 300 1
N a1 40 2
N a2 35 3
N 2 60 4
N b1 40 5
N b2 30 6
N 3 300 7
P ref 1+,a1+,2+,b1+,3+
P h1 1+,a2+,2+,b2+,3+
P h2 1+,a1+,2+,b2+,3+
B 2 3 b1;b2
B 1 2 a1;a2
EOF2
"$BIN" refine -i "$OUT/revchain/g.gfa" -b "$OUT/revchain/g" -r ref -o "$OUT/rvc" --min-bubbles 3 -q >/dev/null 2>&1
[ "$(rrow "$OUT/rvc.refine.report.tsv" reason) $(rrow "$OUT/rvc.refine.report.tsv" source) $(rrow "$OUT/rvc.refine.report.tsv" sink)" = "min_bubbles 1 3" ] \
  && ok "a region skipped before its crossings are read reports its outer anchors (1..3)" \
  || bad "skipped region anchors: $(cut -f3-6 "$OUT/rvc.refine.report.tsv" | tail -1)"
# A unit graph table named explicitly must exist; one implied by -b is warned about, even with -q.
"$BIN" refine -i "$OUT/revchain/g.gfa" -b "$OUT/revchain/g" -r ref -o "$OUT/rvm" -q \
  --panphorte-prefix "$OUT/no_such_run" >/dev/null 2>"$OUT/rvm.err"; rc=$?
[ "$rc" -ne 0 ] && grep -q "no_such_run.panphorte.unit_graph.tsv" "$OUT/rvm.err" && [ ! -e "$OUT/rvm.normalized.sorted.gfa" ] \
  && ok "a missing table named by --panphorte-prefix is an error" \
  || bad "missing explicit table: exit $rc, $(head -1 "$OUT/rvm.err")"
"$BIN" refine -i "$OUT/revchain/g.gfa" -b "$OUT/revchain/g" -r ref -o "$OUT/rvq" -q >/dev/null 2>"$OUT/rvq.err"
grep -q "warning: no panphorte unit graph table" "$OUT/rvq.err" \
  && ok "a missing implied table is warned about with -q" \
  || bad "no warning for a missing implied table under -q"
# A CSV whose name does not end in .bubbles.csv implies no table: the warning says to name one.
cp "$OUT/disb.bubbles.csv" "$OUT/sites.csv"
"$BIN" refine -i "$OUT/disb.sorted.gfa" -c "$OUT/sites.csv" -r ref -o "$OUT/rsc" -q >/dev/null 2>"$OUT/rsc.err"
grep -q "warning: no panphorte unit graph table.*--panphorte-prefix" "$OUT/rsc.err" && ! grep -q "'.panphorte" "$OUT/rsc.err" \
  && ok "a CSV with no table prefix is warned about with --panphorte-prefix named" \
  || bad "warning for a CSV with no table prefix: $(head -1 "$OUT/rsc.err")"
# An inversion site is in the CSV, but refine does not realign it: --bubble-id naming it says so.
{ cat "$OUT/disb.bubbles.csv"; printf '3,1,+,8,+,6,8,3,2,2,1,1,470,476,"2;3;4;5;6;7",inversion,\n'; } > "$OUT/inv.bubbles.csv"
"$BIN" refine -i "$OUT/disb.sorted.gfa" -c "$OUT/inv.bubbles.csv" -r ref -o "$OUT/rinv" --bubble-id 3 -q \
  >/dev/null 2>"$OUT/rinv.err"; rc=$?
[ "$rc" -ne 0 ] && grep -q "names an inversion site" "$OUT/rinv.err" && ! grep -q "not present" "$OUT/rinv.err" \
  && ok "--bubble-id naming an inversion site is refused as one, not as absent" \
  || bad "--bubble-id on an inversion site: exit $rc, $(head -1 "$OUT/rinv.err")"
# The CSV is read once, so the same holds when it comes from a pipe.
"$BIN" refine -i "$OUT/disb.sorted.gfa" -c <(cat "$OUT/inv.bubbles.csv") -r ref -o "$OUT/rinvp" --bubble-id 3 -q \
  >/dev/null 2>"$OUT/rinvp.err"
"$BIN" refine -i "$OUT/disb.sorted.gfa" -c <(cat "$OUT/disb.bubbles.csv") -r ref -o "$OUT/rinvp" --bubble-id 999 -q \
  >/dev/null 2>>"$OUT/rinvp.err"
grep -q "names an inversion site" "$OUT/rinvp.err" && grep -q "not present in the bubbles CSV: 999" "$OUT/rinvp.err" \
  && ok "--bubble-id errors name the inversion site or the absent id with a piped CSV" \
  || bad "--bubble-id with a piped CSV: $(cat "$OUT/rinvp.err" | tr '\n' ' ')"

echo
if [ "$fails" -eq 0 ]; then echo "refine_stats: all assertions passed"; exit 0; fi
echo "refine_stats: $fails assertion(s) FAILED"; exit 1
