#!/usr/bin/env bash
# panphorte_stats.sh - contract assertions for `panvar panphorte`.
#
# The registered C4 smoke test passes with five normalized=no rows, so it never exercises folding at
# all. Every fixture here is small enough that the expected outcome is derivable by hand.
#
# One structural note that shapes the fixtures: the cactus decomposition roots at the longest path, so
# on a graph where the bubble IS most of the graph the site is absorbed into the root and no snarl is
# reported. Every fixture therefore carries a backbone longer than its alleles. This is a property of
# the decomposition, not of panphorte.
#
#   panphorte_stats.sh <panvar-binary> <out-dir>
set -uo pipefail

BIN="${1:?usage: panphorte_stats.sh <panvar> <outdir>}"
OUT="${2:?}"
mkdir -p "$OUT"
OUT="$OUT/run.$$.$(date +%s)"
rm -rf "$OUT"; mkdir -p "$OUT"
fails=0
ok()  { printf "  ok   %s\n" "$1"; }
bad() { printf "  FAIL %s\n" "$1"; fails=$((fails + 1)); }
# value of <column> in the first data row of the report
rep() { awk -F'\t' -v w="$2" 'NR==1{for(i=1;i<=NF;i++)c[$i]=i;next} NR==2{print $(c[w])}' "$1"; }
pline() { awk -v p="$2" '$1=="P" && $2==p {print $3}' "$1"; }

UNIT="ACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCA"   # 64 bp
BB1="$(head -c 2000 /dev/zero | tr '\0' 'A')"
BB2="$(head -c 2000 /dev/zero | tr '\0' 'G')"

# ---------------------------------------------------------------- prevalence counts every haplotype
# One two-copy carrier and nine haplotypes that cross the site directly. True prevalence is 1/10, so a
# 50% gate must refuse. Counting only haplotypes with an interior node made the denominator the set of
# carriers, and this read 1/1.
{ printf 'H\tVN:Z:1.0\n'
  printf "S\t1\t%s\nS\t2\t%s\nS\t3\t%s\n" "$BB1" "$UNIT" "$BB2"
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t1\t+\t3\t+\t0M\n'
  printf 'P\tarr\t1+,2+,2+,3+\t*\n'
  for i in 1 2 3 4 5 6 7 8 9; do printf "P\td%s\t1+,3+\t*\n" "$i"; done; } > "$OUT/prev.gfa"
"$BIN" bubble -i "$OUT/prev.gfa" -r arr -o "$OUT/prevb" --min-variant-bp 0 -q >/dev/null 2>&1
"$BIN" panphorte -i "$OUT/prevb.sorted.gfa" -c "$OUT/prevb.bubbles.csv" -o "$OUT/prevp" -q >/dev/null 2>&1
[ "$(rep "$OUT/prevp.panphorte.report.tsv" normalized)" = "no" ] \
  && ok "one carrier among nine deletion alleles does not clear a 50% gate" \
  || bad "prevalence 1/10 was normalized at the default gate"
[ "$(rep "$OUT/prevp.panphorte.report.tsv" n_traversing)" = "10" ] \
  && ok "the prevalence denominator counts direct source->sink crossings" \
  || bad "n_traversing is $(rep "$OUT/prevp.panphorte.report.tsv" n_traversing), expected 10"
"$BIN" panphorte -i "$OUT/prevb.sorted.gfa" -c "$OUT/prevb.bubbles.csv" -o "$OUT/prevp2" \
       --min-array-prevalence 0.05 -q >/dev/null 2>&1
[ "$(rep "$OUT/prevp2.panphorte.report.tsv" normalized)" = "yes" ] \
  && ok "the same site folds once the gate is lowered below its true prevalence" \
  || bad "prevalence 1/10 did not fold at a 0.05 gate"

# ---------------------------------------------------------------- CN 0/1/2/3 through one site REP
# --min-copies confirms the SITE; it is not a per-haplotype threshold. Every haplotype carrying at
# least one copy must fold through the same node, so copy number read from REP multiplicity is 3/2/1/0.
{ printf 'H\tVN:Z:1.0\n'
  printf "S\t1\t%s\nS\t2\t%s\nS\t3\t%s\n" "$BB1" "$UNIT" "$BB2"
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t1\t+\t3\t+\t0M\n'
  printf 'P\tcn3\t1+,2+,2+,2+,3+\t*\nP\tcn2\t1+,2+,2+,3+\t*\n'
  printf 'P\tcn1\t1+,2+,3+\t*\nP\tcn0\t1+,3+\t*\n'; } > "$OUT/cn.gfa"
"$BIN" bubble -i "$OUT/cn.gfa" -r cn3 -o "$OUT/cnb" --min-variant-bp 0 -q >/dev/null 2>&1
"$BIN" panphorte -i "$OUT/cnb.sorted.gfa" -c "$OUT/cnb.bubbles.csv" -o "$OUT/cnp" -q >/dev/null 2>&1
if [ "$(rep "$OUT/cnp.panphorte.report.tsv" normalized)" = "yes" ]; then
  ok "a site carried at 3/2/1/0 copies is folded"
  rep_node=$(pline "$OUT/cnp.normalized.gfa" cn3 | tr ',' '\n' | sed -n '2p' | tr -d '+-')
  for want in "cn3 3" "cn2 2" "cn1 1"; do
    set -- $want
    n=$(pline "$OUT/cnp.normalized.gfa" "$1" | tr ',' '\n' | grep -c "^${rep_node}[+-]$")
    [ "$n" = "$2" ] && ok "$1 traverses the site REP $2x" \
                    || bad "$1 traverses the REP ${n}x, expected $2"
  done
  n0=$(pline "$OUT/cnp.normalized.gfa" cn0 | tr ',' '\n' | grep -c "^${rep_node}[+-]$")
  [ "$n0" = "0" ] && ok "cn0 does not traverse the REP at all" \
                  || bad "cn0 traverses the REP ${n0}x, expected 0"
  [ "$(rep "$OUT/cnp.panphorte.report.tsv" min_copies)" = "1" ] \
    && ok "the report records a one-copy haplotype (min_copies=1)" \
    || bad "min_copies is $(rep "$OUT/cnp.panphorte.report.tsv" min_copies), expected 1"
else
  bad "the 3/2/1/0 site was not folded at all"
fi

# ---------------------------------------------------------------- one motif per site
# Two haplotypes carrying unrelated repeats. Prevalence must be asked of each motif on its own; a
# poly-C array lends no support to a poly-T one.
CU="$(head -c 64 /dev/zero | tr '\0' 'C')"
TU="$(head -c 64 /dev/zero | tr '\0' 'T')"
{ printf 'H\tVN:Z:1.0\n'
  printf "S\t1\t%s\nS\t2\t%s\nS\t3\t%s\nS\t4\t%s\n" "$BB1" "$CU" "$TU" "$BB2"
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t2\t+\t0M\nL\t2\t+\t4\t+\t0M\n'
  printf 'L\t1\t+\t3\t+\t0M\nL\t3\t+\t3\t+\t0M\nL\t3\t+\t4\t+\t0M\n'
  printf 'P\tpC\t1+,2+,2+,4+\t*\nP\tpT\t1+,3+,3+,4+\t*\n'; } > "$OUT/motif.gfa"
"$BIN" bubble -i "$OUT/motif.gfa" -r pC -o "$OUT/motifb" --min-variant-bp 0 -q >/dev/null 2>&1
"$BIN" panphorte -i "$OUT/motifb.sorted.gfa" -c "$OUT/motifb.bubbles.csv" -o "$OUT/motifp" \
       --min-array-prevalence 0.75 -q >/dev/null 2>&1
[ "$(rep "$OUT/motifp.panphorte.report.tsv" n_motifs)" = "2" ] \
  && ok "two unrelated repeats in one bubble are counted as two motifs" \
  || bad "n_motifs is $(rep "$OUT/motifp.panphorte.report.tsv" n_motifs), expected 2"
[ "$(rep "$OUT/motifp.panphorte.report.tsv" normalized)" = "no" ] \
  && ok "neither motif clears a 0.75 gate on its own half of the cohort" \
  || bad "a 0.75 gate was cleared by pooling two unrelated motifs"

# ---------------------------------------------------------------- REP nodes are site-local
# Two separate sites carrying the SAME repeat unit must not share a node: that would join two loci in
# the graph, so a walk could leave one and arrive in the other.
{ printf 'H\tVN:Z:1.0\n'
  printf "S\t1\t%s\nS\t2\t%s\nS\t3\t%s\nS\t4\t%s\nS\t5\t%s\nS\t6\t%s\n" \
         "$BB1" "$UNIT" "$BB2" "$BB1" "$UNIT" "$BB2"
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t3\t+\t4\t+\t0M\n'
  printf 'L\t4\t+\t5\t+\t0M\nL\t5\t+\t5\t+\t0M\nL\t5\t+\t6\t+\t0M\n'
  printf 'L\t1\t+\t3\t+\t0M\nL\t4\t+\t6\t+\t0M\n'
  printf 'P\tpA\t1+,2+,2+,3+,4+,5+,5+,6+\t*\nP\tpB\t1+,2+,2+,3+,4+,5+,5+,6+\t*\n'
  printf 'P\tpC\t1+,3+,4+,6+\t*\n'; } > "$OUT/two.gfa"
"$BIN" bubble -i "$OUT/two.gfa" -r pA -o "$OUT/twob" --min-variant-bp 0 -q >/dev/null 2>&1
"$BIN" panphorte -i "$OUT/twob.sorted.gfa" -c "$OUT/twob.bubbles.csv" -o "$OUT/twop" \
       --min-array-prevalence 0.5 -q >/dev/null 2>&1
reps=$(pline "$OUT/twop.normalized.gfa" pA | tr ',' '\n' | tr -d '+-' | sort | uniq -c \
       | awk '$1>=2{print $2}' | wc -l | tr -d ' ')
[ "$reps" -ge 2 ] \
  && ok "two sites with an identical repeat unit get distinct REP nodes" \
  || bad "the two sites share a REP node, joining separate loci through one node"

# ---------------------------------------------------------------- rotation and strand resolve to one motif
# The same repeat written starting one base later, and on the other strand, is the same motif.
ROT="CGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAA"
{ printf 'H\tVN:Z:1.0\n'
  printf "S\t1\t%s\nS\t2\t%s\nS\t3\t%s\nS\t4\t%s\n" "$BB1" "$UNIT" "$ROT" "$BB2"
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t2\t+\t0M\nL\t2\t+\t4\t+\t0M\n'
  printf 'L\t1\t+\t3\t+\t0M\nL\t3\t+\t3\t+\t0M\nL\t3\t+\t4\t+\t0M\n'
  printf 'P\tfwd\t1+,2+,2+,4+\t*\nP\trot\t1+,3+,3+,4+\t*\n'; } > "$OUT/rot.gfa"
"$BIN" bubble -i "$OUT/rot.gfa" -r fwd -o "$OUT/rotb" --min-variant-bp 0 -q >/dev/null 2>&1
"$BIN" panphorte -i "$OUT/rotb.sorted.gfa" -c "$OUT/rotb.bubbles.csv" -o "$OUT/rotp" -q >/dev/null 2>&1
[ "$(rep "$OUT/rotp.panphorte.report.tsv" n_motifs)" = "1" ] \
  && ok "a rotated unit resolves to the same motif" \
  || bad "n_motifs is $(rep "$OUT/rotp.panphorte.report.tsv" n_motifs), expected 1 for a rotation"

# ---------------------------------------------------------------- exact sequence preservation
# Exact mode replaces identical copies with a node spelling the same bases, so every haplotype must
# still spell exactly what it spelled before.
"$BIN" panphorte -i "$OUT/cnb.sorted.gfa" -c "$OUT/cnb.bubbles.csv" -o "$OUT/keep" -q >/dev/null 2>&1
spell() {
  awk -v want="$2" '
    $1=="S"{seq[$2]=$3}
    $1=="P" && $2==want{
      n=split($3,a,","); out=""
      for(i=1;i<=n;i++){ s=a[i]; id=substr(s,1,length(s)-1); o=substr(s,length(s))
        t=seq[id]
        if(o=="-"){ r=""; for(j=length(t);j>0;j--){ c=substr(t,j,1)
            r=r (c=="A"?"T":c=="T"?"A":c=="C"?"G":c=="G"?"C":c) } t=r }
        out=out t }
      print out }' "$1"
}
allsame=1
for h in cn3 cn2 cn1 cn0; do
  [ "$(spell "$OUT/cnb.sorted.gfa" "$h")" = "$(spell "$OUT/keep.normalized.gfa" "$h")" ] || allsame=0
done
[ "$allsame" = "1" ] \
  && ok "exact folding preserves every haplotype's spelled sequence byte for byte" \
  || bad "exact folding changed a haplotype's sequence"

# ---------------------------------------------------------------- contract and inputs
{ printf 'H\tVN:Z:1.0\n'
  printf "S\t1\t%s\nS\t2\t%s\nS\t3\t%s\n" "$BB1" "$UNIT" "$BB2"
  printf 'L\t1\t+\t2\t+\t0M\nL\t1\t+\t3\t+\t0M\n'
  printf 'P\tarr\t1+,2+,3+\t*\nP\tdel\t1+,3+\t*\n'; } > "$OUT/bad.gfa"
"$BIN" panphorte -i "$OUT/bad.gfa" -c "$OUT/cnb.bubbles.csv" -o "$OUT/badp" -q >/dev/null 2>&1
[ "$?" -ne 0 ] && ok "a path step with no link behind it is refused" \
               || bad "a malformed graph was accepted and silently repaired"

printf 'bubble_id,source,sink,path_support,min_inside_bp,max_inside_bp,inside_nodes\n9,900,902,2,10,10,"901"\n' \
  > "$OUT/ghost.csv"
"$BIN" panphorte -i "$OUT/cnb.sorted.gfa" -c "$OUT/ghost.csv" -o "$OUT/ghostp" -q >/dev/null 2>&1
[ "$?" -ne 0 ] && ok "a bubbles CSV naming nodes absent from the GFA is refused" \
               || bad "a mismatched CSV ran to completion"

"$BIN" panphorte -i "$OUT/cnb.sorted.gfa" -c "$OUT/cnb.bubbles.csv" -o "$OUT/nope" --bubble-id 4242 \
       -q >/dev/null 2>&1
[ "$?" -ne 0 ] && ok "a --bubble-id that names nothing is refused" \
               || bad "--bubble-id 4242 ran over an empty selection"

"$BIN" panphorte -i "$OUT/cnb.sorted.gfa" -c "$OUT/cnb.bubbles.csv" -o "$OUT/cnb" -q >/dev/null 2>&1
[ "$?" -ne 0 ] && ok "an output that would overwrite an input is refused" \
               || bad "an output was allowed to overwrite an input"

# ---------------------------------------------------------------- reference spelling
"$BIN" panphorte -i "$OUT/cnb.sorted.gfa" -c "$OUT/cnb.bubbles.csv" -o "$OUT/rx" -r cn3 \
       --resnarl-min-variant-bp 0 -q >/dev/null 2>&1
"$BIN" panphorte -i "$OUT/cnb.sorted.gfa" -c "$OUT/cnb.bubbles.csv" -o "$OUT/rX" -r CN3 \
       --resnarl-min-variant-bp 0 -q >/dev/null 2>&1
if [ -s "$OUT/rx.bubbles.csv" ] && [ -s "$OUT/rX.bubbles.csv" ]; then
  cmp -s "$OUT/rx.bubbles.csv" "$OUT/rX.bubbles.csv" \
    && ok "a case-insensitive reference alias gives the same re-snarled CSV" \
    || bad "-r CN3 and -r cn3 produced different call-ready CSVs"
else
  bad "no re-snarled CSV was produced for the reference-spelling check"
fi

# ---------------------------------------------------------------- determinism
"$BIN" panphorte -i "$OUT/cnb.sorted.gfa" -c "$OUT/cnb.bubbles.csv" -o "$OUT/t1" \
       --min-similarity 0.90 --threads 1 -q >/dev/null 2>&1
"$BIN" panphorte -i "$OUT/cnb.sorted.gfa" -c "$OUT/cnb.bubbles.csv" -o "$OUT/t8" \
       --min-similarity 0.90 --threads 8 -q >/dev/null 2>&1
cmp -s "$OUT/t1.normalized.gfa" "$OUT/t8.normalized.gfa" \
  && ok "approximate output is byte-identical at 1 and 8 threads" \
  || bad "approximate output differs between 1 and 8 threads"

echo
if [ "$fails" -eq 0 ]; then echo "panphorte_stats: all assertions passed"; exit 0; fi
echo "panphorte_stats: $fails assertion(s) FAILED"; exit 1
