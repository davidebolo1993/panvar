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

# ---------------------------------------------------------------- rotation resolves to one motif
# The same repeat written starting one base later is the same MOTIF, and prevalence is asked of it
# once rather than of two half-supported candidates. It does NOT follow that both phases share one REP
# node: two phase-rotated linear sequences cannot both be spelled by one unsplit node while exact
# spelling is preserved, so each phase gets its own REP. Downstream copy number therefore still needs a
# REP -> (site, motif, phase) mapping to know they are the same site; that mapping is not implemented,
# and this test asserts the motif count only. (Strand equivalence is a separate claim and is not
# tested here -- this fixture is a rotation, not a reverse-oriented traversal.)
ROT="CGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAA"
{ printf 'H\tVN:Z:1.0\n'
  printf "S\t1\t%s\nS\t2\t%s\nS\t3\t%s\nS\t4\t%s\n" "$BB1" "$UNIT" "$ROT" "$BB2"
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t2\t+\t0M\nL\t2\t+\t4\t+\t0M\n'
  printf 'L\t1\t+\t3\t+\t0M\nL\t3\t+\t3\t+\t0M\nL\t3\t+\t4\t+\t0M\n'
  printf 'P\tfwd\t1+,2+,2+,4+\t*\nP\trot\t1+,3+,3+,4+\t*\n'; } > "$OUT/rot.gfa"
"$BIN" bubble -i "$OUT/rot.gfa" -r fwd -o "$OUT/rotb" --min-variant-bp 0 -q >/dev/null 2>&1
"$BIN" panphorte -i "$OUT/rotb.sorted.gfa" -c "$OUT/rotb.bubbles.csv" -o "$OUT/rotp" -q >/dev/null 2>&1
[ "$(rep "$OUT/rotp.panphorte.report.tsv" n_motifs)" = "1" ] \
  && ok "a rotated unit is counted as one motif, not two half-supported candidates" \
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

# ---------------------------------------------------------------- approximate mode: direct alleles
# The approximate branch kept the inside-node-only interval finder after the exact branch was fixed,
# so one array carrier among nine direct alleles still read prevalence 1/1 and folded through a 50%
# gate.
"$BIN" panphorte -i "$OUT/prevb.sorted.gfa" -c "$OUT/prevb.bubbles.csv" -o "$OUT/apx" \
       --min-similarity 0.90 -q >/dev/null 2>&1
[ "$(rep "$OUT/apx.panphorte.report.tsv" n_traversing)" = "10" ] \
  && ok "approximate mode counts direct source->sink crossings in the denominator" \
  || bad "approximate n_traversing is $(rep "$OUT/apx.panphorte.report.tsv" n_traversing), expected 10"
[ "$(rep "$OUT/apx.panphorte.report.tsv" normalized)" = "no" ] \
  && ok "approximate mode refuses prevalence 1/10 at the default gate" \
  || bad "approximate mode folded a site carried by one haplotype in ten"

# ---------------------------------------------------------------- approximate mode: interruptions
# Two copies of a 16 bp unit separated by a 64 bp gap: the gap is 64/96 = 0.667 of the span, so the
# pair is one array only when --max-interruption-frac allows that much. The option was applied when
# seeding the unit and nowhere afterwards, so the pair was accepted at every setting and
# interruptions_bp always read 0.
U16="ACGTTGCAACGTTGCA"
G64="TTTTTTTTGGGGCCCCTTTTTTTTGGGGCCCCTTTTTTTTGGGGCCCCTTTTTTTTGGGGCCCC"
{ printf 'H\tVN:Z:1.0\n'
  printf "S\t1\t%s\nS\t2\t%s\nS\t3\t%s\nS\t4\t%s\n" "$BB1" "$U16" "$G64" "$BB2"
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t2\t+\t0M\nL\t2\t+\t4\t+\t0M\n'
  printf 'L\t2\t+\t3\t+\t0M\nL\t3\t+\t2\t+\t0M\n'
  printf 'P\tclean\t1+,2+,2+,4+\t*\nP\tsep\t1+,2+,3+,2+,4+\t*\n'; } > "$OUT/gap.gfa"
"$BIN" bubble -i "$OUT/gap.gfa" -r clean -o "$OUT/gapb" --min-variant-bp 0 -q >/dev/null 2>&1
"$BIN" panphorte -i "$OUT/gapb.sorted.gfa" -c "$OUT/gapb.bubbles.csv" -o "$OUT/g0" \
       --min-similarity 0.90 --min-unit-bp 16 --max-interruption-frac 0 --min-array-prevalence 0.5 \
       -q >/dev/null 2>&1
"$BIN" panphorte -i "$OUT/gapb.sorted.gfa" -c "$OUT/gapb.bubbles.csv" -o "$OUT/g9" \
       --min-similarity 0.90 --min-unit-bp 16 --max-interruption-frac 0.9 --min-array-prevalence 0.5 \
       -q >/dev/null 2>&1
[ "$(rep "$OUT/g0.panphorte.report.tsv" n_motif_carriers)" = "1" ] \
  && ok "at --max-interruption-frac 0 a 0.667 gap splits the pair into two arrays of one" \
  || bad "frac 0 gave $(rep "$OUT/g0.panphorte.report.tsv" n_motif_carriers) carriers, expected 1"
[ "$(rep "$OUT/g9.panphorte.report.tsv" n_motif_carriers)" = "2" ] \
  && ok "at --max-interruption-frac 0.9 the same gap is tolerated and the pair is one array" \
  || bad "frac 0.9 gave $(rep "$OUT/g9.panphorte.report.tsv" n_motif_carriers) carriers, expected 2"
[ "$(rep "$OUT/g9.panphorte.report.tsv" interruptions_bp)" = "64" ] \
  && ok "the interrupting bases accepted into an array are reported (64 bp)" \
  || bad "interruptions_bp is $(rep "$OUT/g9.panphorte.report.tsv" interruptions_bp), expected 64"

# ---------------------------------------------------------------- approximate mode: partial boundaries
# A copy whose boundary falls inside a node cannot be folded without splitting that node, and rounding
# to the nearest boundary DELETES the bases between the node edge and the copy edge -- sequence outside
# the copy. Until splitting exists the copy is declined, so the haplotype keeps its original nodes and
# spells exactly what it spelled before.
{ printf 'H\tVN:Z:1.0\n'
  printf "S\t1\t%s\nS\t2\t%s\nS\t3\tTTT%sGGG\nS\t4\t%s\n" "$BB1" "$UNIT" "$UNIT" "$BB2"
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t2\t+\t0M\nL\t2\t+\t4\t+\t0M\n'
  printf 'L\t1\t+\t3\t+\t0M\nL\t3\t+\t4\t+\t0M\n'
  printf 'P\tcleanA\t1+,2+,2+,4+\t*\nP\tcleanB\t1+,2+,2+,4+\t*\nP\tflank\t1+,3+,4+\t*\n'; } \
  > "$OUT/flank.gfa"
"$BIN" bubble -i "$OUT/flank.gfa" -r cleanA -o "$OUT/flankb" --min-variant-bp 0 -q >/dev/null 2>&1
before_flank=$(spell "$OUT/flankb.sorted.gfa" flank)
"$BIN" panphorte -i "$OUT/flankb.sorted.gfa" -c "$OUT/flankb.bubbles.csv" -o "$OUT/flankp" \
       --min-similarity 0.90 --min-array-prevalence 0.5 -q >/dev/null 2>&1
[ "$(spell "$OUT/flankp.normalized.gfa" flank)" = "$before_flank" ] \
  && ok "a copy with a mid-node boundary is declined, not truncated to the node edge" \
  || bad "the flanking haplotype lost $(( ${#before_flank} - $(spell "$OUT/flankp.normalized.gfa" flank | wc -c) )) bases"
# Site-wide, not per copy: folding the rest would leave the site half REP and half literal, and `call`
# counts REP occurrences -- so the refused haplotype would be reported CN 0 while carrying one copy.
# The flank haplotype's ONLY copy is unfoldable, so it would reach the site literally while its
# neighbours use the REP -- that is the false-zero condition, and the site is refused.
[ "$(rep "$OUT/flankp.panphorte.report.tsv" normalized)" = "no" ] \
  && ok "a site is refused when a haplotype would be left with no foldable copy" \
  || bad "the site normalized around a haplotype left entirely literal"
[ "$(rep "$OUT/flankp.panphorte.report.tsv" status)" = "partial_boundary" ] \
  && ok "the refusal is reported as partial_boundary" \
  || bad "status is $(rep "$OUT/flankp.panphorte.report.tsv" status), expected partial_boundary"
for h in cleanA cleanB flank; do
  [ "$(spell "$OUT/flankb.sorted.gfa" $h)" = "$(spell "$OUT/flankp.normalized.gfa" $h)" ] \
    || bad "$h changed sequence at a refused site"
done
ok "every haplotype at a refused site keeps its sequence"

# ---------------------------------------------------------------- a refused site is never called CN 0
# The point of site-wide refusal: a half-REP, half-literal site makes `call` report CN 0 for a
# haplotype that carries a copy. End to end, through call --cn, no such record may appear.
{ printf 'H\tVN:Z:1.0\n'
  printf "S\t1\t%s\nS\t2\t%s\nS\t5\t%s\nS\t3\tTTT%sGGG\nS\t4\t%s\n" \
         "$BB1" "$UNIT" "$UNIT" "$UNIT" "$BB2"
  for e in "1 2" "2 5" "5 4" "1 3" "3 4"; do
    set -- $e; printf "L\t%s\t+\t%s\t+\t0M\n" "$1" "$2"
  done
  printf 'P\tcleanA\t1+,2+,5+,4+\t*\nP\tcleanB\t1+,2+,5+,4+\t*\nP\tflank\t1+,3+,4+\t*\n'; } > "$OUT/e2e.gfa"
"$BIN" bubble -i "$OUT/e2e.gfa" -r cleanA -o "$OUT/e2eb" --min-variant-bp 0 -q >/dev/null 2>&1
"$BIN" panphorte -i "$OUT/e2eb.sorted.gfa" -c "$OUT/e2eb.bubbles.csv" -o "$OUT/e2ep" \
       --min-similarity 0.90 --min-array-prevalence 0.5 -r cleanA --resnarl-min-variant-bp 0 \
       -q >/dev/null 2>&1
if [ -s "$OUT/e2ep.bubbles.csv" ]; then
  "$BIN" call -i "$OUT/e2ep.normalized.sorted.gfa" -b "$OUT/e2ep" -r cleanA -o "$OUT/e2ec" --cn \
         -q >/dev/null 2>&1
  if [ -f "$OUT/e2ec.region.vcf" ]; then
    [ "$(grep -c 'CN=0' "$OUT/e2ec.region.vcf")" = "0" ] \
      && ok "a site refused for a partial boundary yields no CN=0 call" \
      || bad "call reports CN=0 for a haplotype that carries a copy: $(grep -m1 'CN=0' "$OUT/e2ec.region.vcf" | cut -c1-80)"
  else
    bad "call produced no region VCF for the refused site"
  fi
else
  bad "the end-to-end fixture produced no call-ready CSV"
fi

# ---------------------------------------------------------------- the refusal is survivable, loudly
"$BIN" panphorte -i "$OUT/e2eb.sorted.gfa" -c "$OUT/e2eb.bubbles.csv" -o "$OUT/e2eo" \
       --min-similarity 0.90 --min-array-prevalence 0.5 --allow-partial-boundary \
       > "$OUT/e2eo.log" 2>&1
[ "$(rep "$OUT/e2eo.panphorte.report.tsv" normalized)" = "yes" ] \
  && ok "--allow-partial-boundary folds the site instead of refusing it" \
  || bad "--allow-partial-boundary did not restore folding"
grep -q "CN 0" "$OUT/e2eo.log" \
  && ok "and warns that the affected haplotypes will be called CN 0" \
  || bad "--allow-partial-boundary folded without warning about the CN consequence"

# ---------------------------------------------------------------- REP provenance is emitted
[ -s "$OUT/cnp.panphorte.rep_provenance.tsv" ] \
  && ok "a REP provenance table is written" \
  || bad "no rep_provenance.tsv was written"
head -1 "$OUT/cnp.panphorte.rep_provenance.tsv" | grep -q "canonical_motif" \
  && ok "provenance names the site motif each REP stands for" \
  || bad "provenance header lacks canonical_motif"

# ---------------------------------------------------------------- a declined copy beside folded ones
# The opposite case, and the common one at a real array: a haplotype has several copies and only one is
# unfoldable. It still reaches the site through the REP node, so its copy number is undercounted by
# that copy rather than zeroed -- refusing the whole site here would throw away the fold for every
# haplotype to prevent an error of one copy in several. Measured at LPA KIV-2: 466 of 466 haplotypes
# fold at least one copy (median 20), so the false-zero condition never arises there.
# `mixed` reaches the motif twice: once on node 2, which is exactly one unit and folds, and once inside
# node 3, whose copy is flanked by 3 bp on each side and cannot be folded without splitting the node.
{ printf 'H\tVN:Z:1.0\n'
  printf "S\t1\t%s\nS\t2\t%s\nS\t5\t%s\nS\t3\tTTT%sGGG\nS\t4\t%s\n" \
         "$BB1" "$UNIT" "$UNIT" "$UNIT" "$BB2"
  for e in "1 2" "2 5" "5 4" "2 3" "3 4"; do
    set -- $e; printf "L\t%s\t+\t%s\t+\t0M\n" "$1" "$2"
  done
  printf 'P\tcleanA\t1+,2+,5+,4+\t*\nP\tcleanB\t1+,2+,5+,4+\t*\nP\tmixed\t1+,2+,3+,4+\t*\n'; } \
  > "$OUT/mix.gfa"
"$BIN" bubble -i "$OUT/mix.gfa" -r cleanA -o "$OUT/mixb" --min-variant-bp 0 -q >/dev/null 2>&1
"$BIN" panphorte -i "$OUT/mixb.sorted.gfa" -c "$OUT/mixb.bubbles.csv" -o "$OUT/mixp" \
       --min-similarity 0.90 --min-array-prevalence 0.5 > "$OUT/mixp.log" 2>&1
if [ "$(rep "$OUT/mixp.panphorte.report.tsv" status)" = "partial_boundary" ]; then
  bad "a site was refused although every haplotype still folds at least one copy"
else
  ok "a declined copy beside folded ones does not refuse the site"
fi

echo
if [ "$fails" -eq 0 ]; then echo "panphorte_stats: all assertions passed"; exit 0; fi
echo "panphorte_stats: $fails assertion(s) FAILED"; exit 1
