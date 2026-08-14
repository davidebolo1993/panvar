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
# value of a FORMAT key for one sample in the first data record of a VCF. CN is a FORMAT field, so
# grepping for the literal "CN=0" can never find a zero-copy sample and the assertion passes vacuously.
fmt() { awk -F'\t' -v s="$2" -v k="$3" '
  /^#CHROM/ {for (i = 10; i <= NF; i++) if ($i == s) col = i; next}
  /^#/ {next}
  col && !done {split($9, f, ":"); for (i in f) if (f[i] == k) j = i;
                split($col, v, ":"); print v[j]; done = 1}' "$1"; }

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
# A copy boundary rarely lands on a node edge, and rounding to the nearest one DELETES the bases between
# the node edge and the copy edge -- sequence outside the copy. The copy is folded anyway: the whole
# containing step range is replaced and those outside bases come back as fragment nodes around the REP
# step, so the haplotype spells exactly what it spelled before AND its copy is counted.
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
  && ok "a copy with a mid-node boundary is not truncated to the node edge" \
  || bad "the flanking haplotype lost $(( ${#before_flank} - $(spell "$OUT/flankp.normalized.gfa" flank | wc -c) )) bases"
# The copy sits inside a 70 bp node with 3 bp of flank on each side. Neither declining it (an undercount
# of a whole repeat unit) nor rounding it away (deleting the flanks) is necessary.
[ "$(rep "$OUT/flankp.panphorte.report.tsv" normalized)" = "yes" ] \
  && ok "a copy with flanking bases in its node still folds" \
  || bad "the copy was not folded, so its haplotype is undercounted by a whole unit"
[ "$(rep "$OUT/flankp.panphorte.report.tsv" copies_declined_partial_boundary)" = "0" ] \
  && ok "no copy is declined for a mid-node boundary" \
  || bad "$(rep "$OUT/flankp.panphorte.report.tsv" copies_declined_partial_boundary) copies declined"
for h in cleanA cleanB flank; do
  [ "$(spell "$OUT/flankb.sorted.gfa" $h)" = "$(spell "$OUT/flankp.normalized.gfa" $h)" ] \
    || bad "$h lost or gained bases when its copy was folded"
done
ok "the flanking bases survive as fragments, byte for byte"

# ---------------------------------------------------------------- no carrier is ever called CN 0
# The condition the whole boundary policy exists to prevent: a half-REP, half-literal site makes `call`
# read CN 0 for a haplotype that carries a copy, because copy number is counted from REP occurrences.
# `flank` carries exactly one copy, inside a node with 3 bp on each side, so its CN must be 1.
# Also the topology check: the nodes the rewrite replaced must be GONE from the delivered graph, not
# left beside the normalized route as a dead branch for the re-snarl to pick up.
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
    e2e_cn=$(fmt "$OUT/e2ec.region.vcf" flank CN)
    [ -n "$e2e_cn" ] && [ "$e2e_cn" -ge 1 ] \
      && ok "a haplotype carrying one copy is called CN $e2e_cn, not 0" \
      || bad "call reports CN '$e2e_cn' for a haplotype that carries a copy"
  else
    bad "call produced no region VCF"
  fi
else
  bad "the end-to-end fixture produced no call-ready CSV"
fi

# Every node the delivered paths walk must exist, every adjacency must have a link, and no segment may
# be left that no path visits: the approximate rewrite used to mark nothing for removal, so the replaced
# nodes and their links survived beside the REP route and the re-snarl swallowed them back in.
awk -v out="$OUT" '
  $1 == "S" {seg[$2] = 1}
  $1 == "L" {link[$2 "" $3 ">" $4 "" $5] = 1}
  $1 == "P" {n = split($3, st, ","); for (i = 1; i <= n; i++) {
               node = substr(st[i], 1, length(st[i]) - 1); o = substr(st[i], length(st[i]))
               used[node] = 1
               if (!(node in seg)) absent++
               if (i > 1) {
                 fo = (po == "+") ? "-" : "+"; ro = (o == "+") ? "-" : "+"
                 if (!((pn "" po ">" node "" o) in link) && !((node "" ro ">" pn "" fo) in link)) dangling++
               }
               pn = node; po = o }}
  END {for (s in seg) if (!(s in used)) orphan++
       printf "%d %d %d\n", absent + 0, dangling + 0, orphan + 0}' \
  "$OUT/e2ep.normalized.sorted.gfa" > "$OUT/e2e.topo"
read -r e2e_absent e2e_dangling e2e_orphan < "$OUT/e2e.topo"
[ "$e2e_absent" = "0" ] && [ "$e2e_dangling" = "0" ] \
  && ok "the normalized graph has no absent node and no dangling adjacency" \
  || bad "normalized graph: $e2e_absent absent nodes, $e2e_dangling adjacencies with no link"
[ "$e2e_orphan" = "0" ] \
  && ok "no replaced node survives in the delivered graph" \
  || bad "$e2e_orphan replaced node(s) left behind, so the old branch survives beside the REP route"
[ "$(rep "$OUT/e2ep.panphorte.report.tsv" nodes_collapsed)" != "0" ] \
  && ok "nodes_collapsed counts the nodes the rewrite replaced" \
  || bad "nodes_collapsed is 0 although the rewrite replaced nodes"
# The re-snarled site must be exactly the normalized route -- fragment, REP, fragment -- with no extra
# interior node standing for the branch the rewrite was supposed to remove.
e2e_inside=$(awk -F',' 'NR==2{print $6}' "$OUT/e2ep.bubbles.csv")
[ "$e2e_inside" = "3" ] \
  && ok "the re-snarled bubble holds only the normalized route (3 interior nodes)" \
  || bad "the re-snarled bubble has $e2e_inside interior nodes, so a dead branch was picked up"

# ---------------------------------------------------------------- overlapping sites are refused
# Two bubbles claiming the same interior describe one piece of sequence twice, and folding both rewrites
# one span inside the other -- the same array folded at two scales. The acceptance check catches it, but
# only after every haplotype has been aligned, so the input is rejected up front instead. This is not
# hypothetical: `bubble` emits such a pair at ANKRD36C, where a site enclosing all ten others normalized
# with the same 5616 bp unit as the site inside it.
{ printf 'bubble_id,source,source_orient,sink,sink_orient,inside_node_count,total_node_count,'
  printf 'path_support,distinct_alleles,ref_allele_support,alt_allele_support_max,'
  printf 'alt_allele_support_min,min_inside_bp,max_inside_bp,inside_nodes\n'
  printf '1,1,+,4,+,2,4,3,2,2,1,1,64,70,"2;5"\n'
  printf '2,1,+,4,+,2,4,3,2,2,1,1,64,70,"5;3"\n'; } > "$OUT/nested.bubbles.csv"
"$BIN" panphorte -i "$OUT/e2e.gfa" -c "$OUT/nested.bubbles.csv" -o "$OUT/nestedp" \
       --min-similarity 0.90 -q > "$OUT/nested.log" 2>&1
nested_rc=$?
[ "$nested_rc" != "0" ] && grep -q "both claim interior node" "$OUT/nested.log" \
  && ok "two bubbles claiming the same interior node are refused, naming both" \
  || bad "overlapping sites were accepted (exit $nested_rc): $(head -1 "$OUT/nested.log")"
[ -z "$(ls "$OUT"/nestedp.* 2>/dev/null)" ] \
  && ok "the refused run writes nothing" \
  || bad "the refused run left output behind"

# ---------------------------------------------------------------- a replaced node that survives
# A node the rewrite replaced HERE may still be walked somewhere else, so it stays -- but its arcs at
# this site do not. Removing the node is not enough on its own: the link is what keeps the old route
# walkable, and the re-snarl reads links, not paths. `flank2` walks node 2 downstream of the bubble, so
# 2 survives while 1+>2+, 2+>2+ and 2+>4+ must all go.
{ printf 'H\tVN:Z:1.0\n'
  printf "S\t1\t%s\nS\t2\t%s\nS\t3\tTTT%sGGG\nS\t4\t%s\nS\t5\t%s\nS\t6\t%s\n" \
         "$BB1" "$UNIT" "$UNIT" "$BB2" "$BB1" "$BB2"
  for e in "1 2" "2 2" "2 4" "1 3" "3 4" "4 5" "5 2" "2 6" "5 6"; do
    set -- $e; printf "L\t%s\t+\t%s\t+\t0M\n" "$1" "$2"
  done
  printf 'P\tcleanA\t1+,2+,2+,4+,5+,6+\t*\nP\tcleanB\t1+,2+,2+,4+,5+,6+\t*\n'
  printf 'P\tflank2\t1+,3+,4+,5+,2+,6+\t*\n'; } > "$OUT/reuse.gfa"
# Hand-written CSV: run through `bubble` the reused node drags the whole graph into one snarl, and the
# point here is a site whose replaced node is walked OUTSIDE it.
{ printf 'bubble_id,source,source_orient,sink,sink_orient,inside_node_count,total_node_count,'
  printf 'path_support,distinct_alleles,ref_allele_support,alt_allele_support_max,'
  printf 'alt_allele_support_min,min_inside_bp,max_inside_bp,inside_nodes\n'
  printf '1,1,+,4,+,2,4,3,2,2,1,1,64,70,"2;3"\n'; } > "$OUT/reuse.bubbles.csv"
"$BIN" panphorte -i "$OUT/reuse.gfa" -c "$OUT/reuse.bubbles.csv" -o "$OUT/reusep" \
       --min-similarity 0.90 --min-array-prevalence 0.5 -q >/dev/null 2>&1
if [ -s "$OUT/reusep.normalized.gfa" ]; then
  grep -q '^S	2	' "$OUT/reusep.normalized.gfa" \
    && ok "a replaced node still walked elsewhere is kept" \
    || bad "a node another path still walks was removed"
  stale=0
  for l in '1	+	2	+' '2	+	2	+' '2	+	4	+'; do
    grep -q "^L	$l" "$OUT/reusep.normalized.gfa" && stale=$((stale + 1))
  done
  [ "$stale" = "0" ] \
    && ok "the arcs only the pre-rewrite paths walked are pruned" \
    || bad "$stale stale link(s) keep the replaced branch walkable"
  live=0
  for l in '5	+	2	+' '2	+	6	+'; do
    grep -q "^L	$l" "$OUT/reusep.normalized.gfa" && live=$((live + 1))
  done
  [ "$live" = "2" ] \
    && ok "the arcs a final path still walks are kept" \
    || bad "pruning removed a link a path still traverses ($live of 2 kept)"
else
  bad "the reuse fixture produced no normalized graph"
fi

# ---------------------------------------------------------------- the same case, reverse-traversed
# A GFA stores each link once, in whichever direction it happens to store it, and a+ -> b+ IS b- -> a-.
# Here the carriers walk the site backwards, so every arc they use is recorded under the dual of the
# stored L. Keyed on the stored direction alone the pruner found none of them and the entire replaced
# branch stayed walkable -- measured: 0 links removed and 1+>2+, 2+>2+, 2+>4+ all still present.
{ printf 'H\tVN:Z:1.0\n'
  printf "S\t1\t%s\nS\t2\t%s\nS\t3\tTTT%sGGG\nS\t4\t%s\nS\t5\t%s\nS\t6\t%s\n" \
         "$BB1" "$UNIT" "$UNIT" "$BB2" "$BB1" "$BB2"
  for e in "1 2" "2 2" "2 4" "1 3" "3 4" "5 2" "2 6"; do
    set -- $e; printf "L\t%s\t+\t%s\t+\t0M\n" "$1" "$2"
  done
  printf 'P\trevA\t4-,2-,2-,1-\t*\nP\trevB\t4-,2-,2-,1-\t*\nP\trevC\t4-,3-,1-\t*\n'
  printf 'P\tkeep\t5+,2+,6+\t*\n'; } > "$OUT/rev.gfa"
"$BIN" panphorte -i "$OUT/rev.gfa" -c "$OUT/reuse.bubbles.csv" -o "$OUT/revp" \
       --min-similarity 0.90 --min-array-prevalence 0.5 -q >/dev/null 2>&1
if [ -s "$OUT/revp.normalized.gfa" ]; then
  grep -q '^S	2	' "$OUT/revp.normalized.gfa" \
    && ok "reverse: a replaced node still walked elsewhere is kept" \
    || bad "reverse: a node another path still walks was removed"
  rstale=0
  for l in '1	+	2	+' '2	+	2	+' '2	+	4	+'; do
    grep -q "^L	$l" "$OUT/revp.normalized.gfa" && rstale=$((rstale + 1))
  done
  [ "$rstale" = "0" ] \
    && ok "reverse: an arc stored one way and walked the other is still recognised as obsolete" \
    || bad "reverse: $rstale stale link(s) survived because the key was not canonical"
  rlive=0
  for l in '5	+	2	+' '2	+	6	+'; do
    grep -q "^L	$l" "$OUT/revp.normalized.gfa" && rlive=$((rlive + 1))
  done
  [ "$rlive" = "2" ] \
    && ok "reverse: the arcs a final path still walks are kept" \
    || bad "reverse: pruning removed a link a path still traverses ($rlive of 2 kept)"
else
  bad "the reverse-traversal fixture produced no normalized graph"
fi

# Same canonical key, second consequence: ensure() adds an L only when the arc is absent in EITHER
# direction. Keyed on the stored direction it re-added every adjacency the graph already held as a
# dual, so a run that folded NOTHING still doubled the link set -- C4 6532 -> 12374, CYP2D6 8171 ->
# 15606. No link may appear twice, once in each direction.
dupe_arcs() {  # links present in both directions, counted once per pair
  awk -F'\t' '$1=="L"{
    fo=($3=="+")?"-":"+"; to=($5=="+")?"-":"+"
    a=$2 $3 ">" $4 $5; b=$4 to ">" $2 fo
    k=(a<b)?a:b; n[k]++
  } END {d=0; for (x in n) if (n[x] > 1) d++; print d}' "$1"
}
[ "$(dupe_arcs "$OUT/revp.normalized.gfa")" = "0" ] \
  && ok "no arc is emitted twice, once in each direction" \
  || bad "$(dupe_arcs "$OUT/revp.normalized.gfa") arc(s) emitted in both directions"

# ---------------------------------------------------------------- the replaced route must not survive
# A link kept because another path needs it is not on its own a defect: at a real array it is the normal
# case -- LPA keeps 2009, because an arc inside the repeat unit is crossed once per copy and folding
# replaces all but the crossing that falls outside the folded span. What matters is whether the
# surviving nodes and links still spell the OLD allele end to end -- then the site is represented twice,
# folded and unfolded, and the re-snarl reports two alleles for one haplotype. Here `stubX` keeps 1+>2+ alive and
# `stubY` keeps 2+>2+ and 2+>4+, neither crossing the site, so cleanA's original route is fully intact
# after its rewrite.
{ printf 'H\tVN:Z:1.0\n'
  printf "S\t1\t%s\nS\t2\t%s\nS\t3\tTTT%sGGG\nS\t4\t%s\nS\t7\t%s\nS\t8\t%s\n" \
         "$BB1" "$UNIT" "$UNIT" "$BB2" "$(head -c 300 /dev/zero | tr '\0' 'T')" \
         "$(head -c 300 /dev/zero | tr '\0' 'C')"
  for e in "1 2" "2 2" "2 4" "1 3" "3 4" "2 7" "8 2"; do
    set -- $e; printf "L\t%s\t+\t%s\t+\t0M\n" "$1" "$2"
  done
  printf 'P\tcleanA\t1+,2+,2+,4+\t*\nP\tcleanB\t1+,2+,2+,4+\t*\nP\tflank\t1+,3+,4+\t*\n'
  printf 'P\tstubX\t1+,2+,7+\t*\nP\tstubY\t8+,2+,2+,4+\t*\n'; } > "$OUT/twin.gfa"
"$BIN" panphorte -i "$OUT/twin.gfa" -c "$OUT/reuse.bubbles.csv" -o "$OUT/twinp" \
       --min-similarity 0.90 -q > "$OUT/twin.log" 2>&1
twin_rc=$?
[ "$twin_rc" != "0" ] && grep -q "ORIGINAL route" "$OUT/twin.log" \
  && ok "a rewrite whose replaced route stays walkable is refused, naming the paths" \
  || bad "a site folded twice was accepted (exit $twin_rc): $(head -c 120 "$OUT/twin.log")"
[ -z "$(ls "$OUT"/twinp.* 2>/dev/null)" ] \
  && ok "and it writes nothing -- there is no override for it" \
  || bad "the refused run left output behind"
# The other side of the same rule: the reuse fixture keeps node 2 and one of its arcs alive, so the
# obsolete-and-live count is non-zero, but the entry and exit arcs go and the old route is broken. It
# must run, not be refused -- refusing on the arc count alone would block LPA.
[ -s "$OUT/reusep.normalized.gfa" ] \
  && ok "a shared link whose route no longer reconnects does not refuse the run" \
  || bad "the reuse fixture was refused although its replaced route is broken"

# ...and the check has to be asked PER EDIT. A haplotype normalized at two sites has two independent old
# routes; testing the whole original walk conflates them, so the first site's route being properly gone
# makes the walk unwalkable and hides the second site's route surviving intact. Here site A's nodes are
# removed outright while `stubX` and `stubY` keep every node and arc of site B alive, so main, main2 and
# flankA are each broken at A and intact at B. Measured: the whole-walk form accepts this fixture and
# writes a graph; the per-edit form rejects all three.
UB="TTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAA"
MID="$(head -c 2000 /dev/zero | tr '\0' 'T')"
S8="$(head -c 300 /dev/zero | tr '\0' 'C')"
S9="$(head -c 300 /dev/zero | tr '\0' 'A')"
{ printf 'H\tVN:Z:1.0\n'
  printf "S\t1\t%s\nS\t2\t%s\nS\t3\tTTT%sGGG\nS\t4\t%s\nS\t5\t%s\nS\t7\t%s\nS\t6\t%s\nS\t8\t%s\nS\t9\t%s\n" \
         "$BB1" "$UNIT" "$UNIT" "$MID" "$UB" "$S9" "$BB2" "$S8" "$S9"
  for e in "1 2" "2 2" "2 4" "1 3" "3 4" "4 5" "5 5" "5 6" "4 7" "7 6" "5 8" "9 5"; do
    set -- $e; printf "L\t%s\t+\t%s\t+\t0M\n" "$1" "$2"
  done
  printf 'P\tmain\t1+,2+,2+,4+,5+,5+,6+\t*\nP\tmain2\t1+,2+,2+,4+,5+,5+,6+\t*\n'
  printf 'P\tflankA\t1+,3+,4+,5+,5+,6+\t*\nP\tflankB\t1+,2+,2+,4+,7+,6+\t*\n'
  printf 'P\tstubX\t4+,5+,8+\t*\nP\tstubY\t9+,5+,5+,6+\t*\n'; } > "$OUT/twosite.gfa"
{ printf 'bubble_id,source,source_orient,sink,sink_orient,inside_node_count,total_node_count,'
  printf 'path_support,distinct_alleles,ref_allele_support,alt_allele_support_max,'
  printf 'alt_allele_support_min,min_inside_bp,max_inside_bp,inside_nodes\n'
  printf '1,1,+,4,+,2,4,4,2,2,1,1,64,300,"2;3"\n'
  printf '2,4,+,6,+,2,4,4,2,2,1,1,64,300,"5;7"\n'; } > "$OUT/twosite.bubbles.csv"
"$BIN" panphorte -i "$OUT/twosite.gfa" -c "$OUT/twosite.bubbles.csv" -o "$OUT/twositep" \
       --min-similarity 0.90 -q > "$OUT/twosite.log" 2>&1
ts_rc=$?
[ "$ts_rc" != "0" ] && grep -q "bubble 2" "$OUT/twosite.log" \
  && ok "a surviving route at the second of two sites is not masked by the first being removed" \
  || bad "the second site's surviving route was missed (exit $ts_rc): $(head -c 140 "$OUT/twosite.log")"

# ---------------------------------------------------------------- overlap refusal respects --bubble-id
# Only the SELECTED sites can overlap: refusing on a pair the run was never going to touch blocked the
# one safe way to work at a locus whose CSV carries an overlap, which is to name one member of it.
for want in 1 2; do
  "$BIN" panphorte -i "$OUT/e2e.gfa" -c "$OUT/nested.bubbles.csv" -o "$OUT/sel$want" \
         --min-similarity 0.90 --bubble-id "$want" -q > "$OUT/sel$want.log" 2>&1
  [ "$?" = "0" ] \
    && ok "selecting overlapping bubble $want on its own runs" \
    || bad "--bubble-id $want was refused although nothing it selects overlaps: $(head -1 "$OUT/sel$want.log")"
done
"$BIN" panphorte -i "$OUT/e2e.gfa" -c "$OUT/nested.bubbles.csv" -o "$OUT/selboth" \
       --min-similarity 0.90 --bubble-id 1 --bubble-id 2 -q > "$OUT/selboth.log" 2>&1
[ "$?" != "0" ] && grep -q "both claim interior node" "$OUT/selboth.log" \
  && ok "selecting both members of an overlapping pair is still refused" \
  || bad "--bubble-id 1 --bubble-id 2 was accepted"

# ---------------------------------------------------------------- two copies inside one node
# Two copies of a 16 bp unit sit in a single 96 bp node, separated by 64 bp. There is no step boundary
# between them, so mapping each copy to its own containing step range declined the second -- at an
# array that is a whole repeat unit missing from the copy number. Copies sharing a range are emitted as
# one block: fragment, REP, fragment, REP, fragment.
U16="ACGTTGCAACGTTGCA"
G64="TTTTTTTTGGGGCCCCTTTTTTTTGGGGCCCCTTTTTTTTGGGGCCCCTTTTTTTTGGGGCCCC"
{ printf 'H\tVN:Z:1.0\n'
  printf "S\t1\t%s\nS\t2\t%s\nS\t3\t%s%s%s\nS\t4\t%s\n" "$BB1" "$U16" "$U16" "$G64" "$U16" "$BB2"
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t2\t+\t0M\nL\t2\t+\t4\t+\t0M\n'
  printf 'L\t1\t+\t3\t+\t0M\nL\t3\t+\t4\t+\t0M\n'
  printf 'P\tclean\t1+,2+,2+,4+\t*\nP\tsep\t1+,3+,4+\t*\n'; } > "$OUT/tio.gfa"
"$BIN" bubble -i "$OUT/tio.gfa" -r clean -o "$OUT/tiob" --min-variant-bp 0 -q >/dev/null 2>&1
sep_before=$(spell "$OUT/tiob.sorted.gfa" sep)
"$BIN" panphorte -i "$OUT/tiob.sorted.gfa" -c "$OUT/tiob.bubbles.csv" -o "$OUT/tiop" \
       --min-similarity 0.90 --min-unit-bp 16 --max-interruption-frac 0.9 --min-array-prevalence 0.5 \
       -q >/dev/null 2>&1
[ "$(rep "$OUT/tiop.panphorte.report.tsv" copies_declined_partial_boundary)" = "0" ] \
  && ok "two copies sharing one node are both folded, neither declined" \
  || bad "$(rep "$OUT/tiop.panphorte.report.tsv" copies_declined_partial_boundary) copies declined for sharing a node"
[ "$(rep "$OUT/tiop.panphorte.report.tsv" max_copies)" = "2" ] \
  && ok "the haplotype's copy number is 2, not 1" \
  || bad "max_copies is $(rep "$OUT/tiop.panphorte.report.tsv" max_copies), expected 2"
[ "$(spell "$OUT/tiop.normalized.gfa" sep)" = "$sep_before" ] \
  && ok "the 64 bp between the two copies survives as a fragment, byte for byte" \
  || bad "the intervening sequence was lost when both copies were folded"
# sep must traverse the REP twice with the gap fragment between: REP, frag, REP.
n=$(awk '$1=="P" && $2=="sep"{print $3}' "$OUT/tiop.normalized.gfa" | tr ',' '\n' | wc -l | tr -d ' ')
[ "$n" = "5" ] \
  && ok "sep is anchor, REP, fragment, REP, anchor (5 steps)" \
  || bad "sep has $n steps, expected 5"


# ---------------------------------------------------------------- REP provenance is emitted
[ -s "$OUT/cnp.panphorte.rep_provenance.tsv" ] \
  && ok "a REP provenance table is written" \
  || bad "no rep_provenance.tsv was written"
head -1 "$OUT/cnp.panphorte.rep_provenance.tsv" | grep -q "canonical_motif" \
  && ok "provenance names the site motif each REP stands for" \
  || bad "provenance header lacks canonical_motif"
# The join key has to name a node in the graph delivered BESIDE it. --reference-path renumbers, so the
# id a REP was created with is not the id it ships under, and the table named nodes that do not exist.
prov_out=$(awk -F'\t' 'NR==2{print $2}' "$OUT/e2ep.panphorte.rep_provenance.tsv")
prov_in=$(awk -F'\t' 'NR==2{print $1}' "$OUT/e2ep.panphorte.rep_provenance.tsv")
[ -n "$prov_out" ] && grep -q "^S	$prov_out	" "$OUT/e2ep.normalized.sorted.gfa" \
  && ok "output_rep_node names a node in the delivered graph" \
  || bad "provenance output_rep_node '$prov_out' is not a segment of the sorted graph"
grep -q "^L	$prov_out	+	$prov_out	+" "$OUT/e2ep.normalized.sorted.gfa" \
  && ok "the delivered REP node carries the self-loop copy number is counted from" \
  || bad "node $prov_out has no self-loop, so it is not the REP"
[ "$prov_in" != "$prov_out" ] \
  && ok "the created id and the delivered id are both reported, and here they differ" \
  || bad "created_rep_node and output_rep_node are equal; the sort renumbering was not applied"

# Sorting can FLIP a REP, and then the sequence panphorte built the node with is not the sequence the
# delivered node spells. Here the reference walks the array on the minus strand, so --flip reverses it.
{ printf 'H\tVN:Z:1.0\n'
  printf "S\t1\t%s\nS\t2\t%s\nS\t3\tTTT%sGGG\nS\t4\t%s\n" "$BB1" "$UNIT" "$UNIT" "$BB2"
  for e in "1 2" "2 2" "2 4" "1 3" "3 4"; do
    set -- $e; printf "L\t%s\t+\t%s\t+\t0M\n" "$1" "$2"
  done
  printf 'P\trefrev\t4-,2-,2-,1-\t*\nP\tcarrier\t1+,2+,2+,4+\t*\nP\tflank\t1+,3+,4+\t*\n'; } \
  > "$OUT/flip.gfa"
"$BIN" panphorte -i "$OUT/flip.gfa" -c "$OUT/reuse.bubbles.csv" -o "$OUT/flipp" \
       --min-similarity 0.90 -r refrev -q >/dev/null 2>&1
fl_out=$(awk -F'\t' 'NR==2{print $2}' "$OUT/flipp.panphorte.rep_provenance.tsv" 2>/dev/null)
fl_made=$(awk -F'\t' 'NR==2{print $5}' "$OUT/flipp.panphorte.rep_provenance.tsv" 2>/dev/null)
fl_deliv=$(awk -F'\t' 'NR==2{print $6}' "$OUT/flipp.panphorte.rep_provenance.tsv" 2>/dev/null)
fl_seq=$(awk -v n="$fl_out" '$1=="S" && $2==n {print $3}' "$OUT/flipp.normalized.sorted.gfa")
[ -n "$fl_seq" ] && [ "$fl_deliv" = "$fl_seq" ] \
  && ok "output_phase_unit is the sequence the delivered REP actually spells" \
  || bad "output_phase_unit does not match segment $fl_out in the sorted graph"
[ -n "$fl_made" ] && [ "$fl_made" != "$fl_deliv" ] \
  && ok "and it differs from created_phase_unit when the sort flips the REP" \
  || bad "the REP was not flipped, so this fixture proves nothing"

# ---------------------------------------------------------------- the GTF is an input like any other
# Checked only against the gene CSV, a GTF naming some other output was read and then overwritten by the
# commit. And clearing it when there is no --reference-path removed it from the preflight too, so the
# gap was widest exactly where the annotation is skipped.
"$BIN" panphorte -i "$OUT/e2e.gfa" -c "$OUT/e2eb.bubbles.csv" -o "$OUT/gawith" \
       --gtf "$OUT/gawith.panphorte.report.tsv" -r cleanA -q > "$OUT/gawith.log" 2>&1
ga_rc=$?
[ "$ga_rc" != "0" ] && grep -q "same file as input" "$OUT/gawith.log" \
  && ok "a GTF that names an output is refused (with --reference-path)" \
  || bad "GTF aliasing an output was accepted with --reference-path (exit $ga_rc)"
"$BIN" panphorte -i "$OUT/e2e.gfa" -c "$OUT/e2eb.bubbles.csv" -o "$OUT/gawithout" \
       --gtf "$OUT/gawithout.panphorte.report.tsv" -q > "$OUT/gawithout.log" 2>&1
ga_rc=$?
[ "$ga_rc" != "0" ] && grep -q "same file as input" "$OUT/gawithout.log" \
  && ok "a GTF that names an output is refused with no --reference-path, where it is not even used" \
  || bad "GTF aliasing an output was accepted without --reference-path (exit $ga_rc)"

# ---------------------------------------------------------------- copies on unequal node boundaries
# `mixed` reaches the motif twice: once on a node that is exactly one unit, and once inside a node whose
# copy is flanked by 3 bp on each side. Both fold, so its copy number is 2 -- the assertion that matters.
# The old check here asked only whether the site was refused, which is now structurally unreachable, and
# the fixture did not even seed (status=no_seed), so it passed whatever the fold did. `alt` bypasses the
# array so that the first unit node is interior rather than the site's own source anchor; without it the
# bubble starts AT that node and each haplotype has a single copy inside, which is what stopped it
# seeding.
{ printf 'H\tVN:Z:1.0\n'
  printf "S\t1\t%s\nS\t2\t%s\nS\t5\t%s\nS\t3\tTTT%sGGG\nS\t4\t%s\nS\t6\t%s\n" \
         "$BB1" "$UNIT" "$UNIT" "$UNIT" "$BB2" "$(head -c 200 /dev/zero | tr '\0' 'C')"
  for e in "1 2" "2 5" "5 4" "2 3" "3 4" "1 6" "6 4"; do
    set -- $e; printf "L\t%s\t+\t%s\t+\t0M\n" "$1" "$2"
  done
  printf 'P\tcleanA\t1+,2+,5+,4+\t*\nP\tcleanB\t1+,2+,5+,4+\t*\nP\tmixed\t1+,2+,3+,4+\t*\n'
  printf 'P\talt\t1+,6+,4+\t*\n'; } > "$OUT/mix.gfa"
"$BIN" bubble -i "$OUT/mix.gfa" -r cleanA -o "$OUT/mixb" --min-variant-bp 0 -q >/dev/null 2>&1
mix_before=$(spell "$OUT/mixb.sorted.gfa" mixed)
"$BIN" panphorte -i "$OUT/mixb.sorted.gfa" -c "$OUT/mixb.bubbles.csv" -o "$OUT/mixp" \
       --min-similarity 0.90 --min-array-prevalence 0.5 > "$OUT/mixp.log" 2>&1
[ "$(rep "$OUT/mixp.panphorte.report.tsv" status)" = "normalized" ] \
  && ok "the mixed-boundary site seeds and folds" \
  || bad "the mixed-boundary site reports $(rep "$OUT/mixp.panphorte.report.tsv" status)"
mix_copies=$(awk -F'\t' '$1=="mixed"{print $4}' "$OUT/mixp.panphorte.copies.tsv" 2>/dev/null)
[ "$mix_copies" = "2" ] \
  && ok "a copy on a node boundary and a copy inside a node both fold (CN 2)" \
  || bad "mixed folds '$mix_copies' copies, expected 2"
[ "$(spell "$OUT/mixp.normalized.gfa" mixed)" = "$mix_before" ] \
  && ok "and the 3 bp flanking each in-node copy survive" \
  || bad "mixed lost bases when its two unequal copies were folded"

echo
if [ "$fails" -eq 0 ]; then echo "panphorte_stats: all assertions passed"; exit 0; fi
echo "panphorte_stats: $fails assertion(s) FAILED"; exit 1
