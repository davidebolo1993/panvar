#!/usr/bin/env bash
# panphorte_stats.sh - contract assertions for `panvar panphorte`.
#
# Every fixture here is small enough that the expected outcome is derivable by hand.
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
# The checks below are for the single-node fold (one unit node a haplotype visits once per copy),
# which --lossy-units keeps; approximate runs here use it. The unit graph, the default, is checked against
# them further down (anchor passes = unit visits, every path spelled exactly).
panph() {
  local a prev="" lossy=""
  for a in "$@"; do
    if [ "$prev" = "--min-similarity" ] && awk -v s="$a" 'BEGIN{exit !((s < 1) || (s > 1 && s < 100))}'; then
      lossy="--lossy-units"
    fi
    prev=$a
  done
  if [ -n "$lossy" ]; then "$BIN" panphorte "$@" "$lossy"; else "$BIN" panphorte "$@"; fi
}
# value of <column> in the first data row of the report
rep() { awk -F'\t' -v w="$2" 'NR==1{for(i=1;i<=NF;i++)c[$i]=i;next} NR==2{print $(c[w])}' "$1"; }
# The one-node folds of a unit graph table (exact mode's units, or --lossy-units ones): their nodes, one per line.
onenode() { awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)c[$i]=i;next} $(c["n_nodes"])==1{print $(c["nodes"])}' "$1"; }
pline() { awk -v p="$2" '$1=="P" && $2==p {print $3}' "$1"; }
# value of a FORMAT key for one sample in the first data record of a VCF. CN is a FORMAT field, so
# grepping for the literal "CN=0" can never find a zero-copy sample and the assertion passes vacuously.
fmt() { awk -F'\t' -v s="$2" -v k="$3" '
  /^#CHROM/ {for (i = 10; i <= NF; i++) if ($i == s) col = i; next}
  /^#/ {next}
  col && !done {split($9, f, ":"); for (i in f) if (f[i] == k) j = i;
                split($col, v, ":"); print v[j]; done = 1}' "$1"; }

UNIT="TATACGCCATTCAATAACAAGCCAATTCGCTGGGTCAACTCCCAGCCAGACGCCCCGTTGGCGC"   # 64 bp, primitive: a unit of repeated short motifs is a short tandem repeat, which panphorte leaves unfolded
# Backbones: random, not homopolymer -- a 2 kb poly-A is a tandem repeat, which the re-snarl extends
# sites to. A Park-Miller generator, exact in any awk; $2 seeds it.
rnd() { awk -v n="$1" -v seed="$2" 'BEGIN{ b="ACGT"; x=seed; for(i=0;i<n;i++){ x=(x*16807)%2147483647
          printf "%s", substr(b, int(x*4/2147483647)+1, 1) } }'; }
BB1="$(rnd 2000 101)"
BB2="$(rnd 2000 202)"

# ---------------------------------------------------------------- prevalence counts every haplotype
# One two-copy carrier and nine haplotypes that cross the site directly. True prevalence is 1/10, so a
# 50% gate must refuse: the denominator is every crossing haplotype, not only those with an interior node.
{ printf 'H\tVN:Z:1.0\n'
  printf "S\t1\t%s\nS\t2\t%s\nS\t3\t%s\n" "$BB1" "$UNIT" "$BB2"
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t1\t+\t3\t+\t0M\n'
  printf 'P\tarr\t1+,2+,2+,3+\t*\n'
  for i in 1 2 3 4 5 6 7 8 9; do printf "P\td%s\t1+,3+\t*\n" "$i"; done; } > "$OUT/prev.gfa"
"$BIN" bubble -i "$OUT/prev.gfa" -r arr -o "$OUT/prevb" --min-variant-bp 0 --merge-nearby-bp 0 -q >/dev/null 2>&1
panph -i "$OUT/prevb.sorted.gfa" -c "$OUT/prevb.bubbles.csv" -o "$OUT/prevp" -q >/dev/null 2>&1
[ "$(rep "$OUT/prevp.panphorte.report.tsv" normalized)" = "no" ] \
  && ok "one carrier among nine deletion alleles does not clear a 50% gate" \
  || bad "prevalence 1/10 was normalized at the default gate"
[ "$(rep "$OUT/prevp.panphorte.report.tsv" n_traversing)" = "10" ] \
  && ok "the prevalence denominator counts direct source->sink crossings" \
  || bad "n_traversing is $(rep "$OUT/prevp.panphorte.report.tsv" n_traversing), expected 10"
panph -i "$OUT/prevb.sorted.gfa" -c "$OUT/prevb.bubbles.csv" -o "$OUT/prevp2" \
       --min-array-prevalence 0.05 -q >/dev/null 2>&1
[ "$(rep "$OUT/prevp2.panphorte.report.tsv" normalized)" = "yes" ] \
  && ok "the same site folds once the gate is lowered below its true prevalence" \
  || bad "prevalence 1/10 did not fold at a 0.05 gate"

# ---------------------------------------------------------------- CN 0/1/2/3 through one site unit
# --min-copies confirms the SITE; it is not a per-haplotype threshold. Every haplotype carrying at
# least one copy must fold through the same node, so copy number read from unit multiplicity is 3/2/1/0.
{ printf 'H\tVN:Z:1.0\n'
  printf "S\t1\t%s\nS\t2\t%s\nS\t3\t%s\n" "$BB1" "$UNIT" "$BB2"
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t1\t+\t3\t+\t0M\n'
  printf 'P\tcn3\t1+,2+,2+,2+,3+\t*\nP\tcn2\t1+,2+,2+,3+\t*\n'
  printf 'P\tcn1\t1+,2+,3+\t*\nP\tcn0\t1+,3+\t*\n'; } > "$OUT/cn.gfa"
"$BIN" bubble -i "$OUT/cn.gfa" -r cn3 -o "$OUT/cnb" --min-variant-bp 0 --merge-nearby-bp 0 -q >/dev/null 2>&1
panph -i "$OUT/cnb.sorted.gfa" -c "$OUT/cnb.bubbles.csv" -o "$OUT/cnp" -q >/dev/null 2>&1
if [ "$(rep "$OUT/cnp.panphorte.report.tsv" normalized)" = "yes" ]; then
  ok "a site carried at 3/2/1/0 copies is folded"
  rep_node=$(pline "$OUT/cnp.normalized.gfa" cn3 | tr ',' '\n' | sed -n '2p' | tr -d '+-')
  for want in "cn3 3" "cn2 2" "cn1 1"; do
    set -- $want
    n=$(pline "$OUT/cnp.normalized.gfa" "$1" | tr ',' '\n' | grep -c "^${rep_node}[+-]$")
    [ "$n" = "$2" ] && ok "$1 traverses the site unit $2x" \
                    || bad "$1 traverses the unit ${n}x, expected $2"
  done
  n0=$(pline "$OUT/cnp.normalized.gfa" cn0 | tr ',' '\n' | grep -c "^${rep_node}[+-]$")
  [ "$n0" = "0" ] && ok "cn0 does not traverse the unit at all" \
                  || bad "cn0 traverses the unit ${n0}x, expected 0"
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
"$BIN" bubble -i "$OUT/motif.gfa" -r pC -o "$OUT/motifb" --min-variant-bp 0 --merge-nearby-bp 0 -q >/dev/null 2>&1
panph -i "$OUT/motifb.sorted.gfa" -c "$OUT/motifb.bubbles.csv" -o "$OUT/motifp" \
       --min-array-prevalence 0.75 -q >/dev/null 2>&1
[ "$(rep "$OUT/motifp.panphorte.report.tsv" n_motifs)" = "2" ] \
  && ok "two unrelated repeats in one bubble are counted as two motifs" \
  || bad "n_motifs is $(rep "$OUT/motifp.panphorte.report.tsv" n_motifs), expected 2"
[ "$(rep "$OUT/motifp.panphorte.report.tsv" normalized)" = "no" ] \
  && ok "neither motif clears a 0.75 gate on its own half of the cohort" \
  || bad "a 0.75 gate was cleared by pooling two unrelated motifs"

# ---------------------------------------------------------------- unit nodes are site-local
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
"$BIN" bubble -i "$OUT/two.gfa" -r pA -o "$OUT/twob" --min-variant-bp 0 --merge-nearby-bp 0 -q >/dev/null 2>&1
panph -i "$OUT/twob.sorted.gfa" -c "$OUT/twob.bubbles.csv" -o "$OUT/twop" \
       --min-array-prevalence 0.5 -q >/dev/null 2>&1
reps=$(pline "$OUT/twop.normalized.gfa" pA | tr ',' '\n' | tr -d '+-' | sort | uniq -c \
       | awk '$1>=2{print $2}' | wc -l | tr -d ' ')
[ "$reps" -ge 2 ] \
  && ok "two sites with an identical repeat unit get distinct unit nodes" \
  || bad "the two sites share a unit node, joining separate loci through one node"

# ---------------------------------------------------------------- rotation resolves to one motif
# The same repeat written starting one base later is the same MOTIF, and prevalence is asked of it
# once rather than of two half-supported candidates. It does NOT follow that both phases share one unit
# node: two phase-rotated linear sequences cannot both be spelled by one unsplit node while exact
# spelling is preserved, so each phase gets its own unit. unit_graph.tsv gives both the same
# input_bubble_id and canonical_motif, which is how call joins them; this test asserts the motif count.
# (Strand equivalence is a separate claim and is not tested here -- this fixture is a rotation, not a
# reverse-oriented traversal.)
ROT="ATACGCCATTCAATAACAAGCCAATTCGCTGGGTCAACTCCCAGCCAGACGCCCCGTTGGCGCT"
{ printf 'H\tVN:Z:1.0\n'
  printf "S\t1\t%s\nS\t2\t%s\nS\t3\t%s\nS\t4\t%s\n" "$BB1" "$UNIT" "$ROT" "$BB2"
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t2\t+\t0M\nL\t2\t+\t4\t+\t0M\n'
  printf 'L\t1\t+\t3\t+\t0M\nL\t3\t+\t3\t+\t0M\nL\t3\t+\t4\t+\t0M\n'
  printf 'P\tfwd\t1+,2+,2+,4+\t*\nP\trot\t1+,3+,3+,4+\t*\n'; } > "$OUT/rot.gfa"
"$BIN" bubble -i "$OUT/rot.gfa" -r fwd -o "$OUT/rotb" --min-variant-bp 0 --merge-nearby-bp 0 -q >/dev/null 2>&1
panph -i "$OUT/rotb.sorted.gfa" -c "$OUT/rotb.bubbles.csv" -o "$OUT/rotp" -q >/dev/null 2>&1
[ "$(rep "$OUT/rotp.panphorte.report.tsv" n_motifs)" = "1" ] \
  && ok "a rotated unit is counted as one motif, not two half-supported candidates" \
  || bad "n_motifs is $(rep "$OUT/rotp.panphorte.report.tsv" n_motifs), expected 1 for a rotation"

# ---------------------------------------------------------------- exact sequence preservation
# Exact mode replaces identical copies with a node spelling the same bases, so every haplotype must
# still spell exactly what it spelled before.
panph -i "$OUT/cnb.sorted.gfa" -c "$OUT/cnb.bubbles.csv" -o "$OUT/keep" -q >/dev/null 2>&1
spell() {
  awk -v want="$2" '
    $1=="S"{seq[$2]=$3}
    $1=="P" && $2==want{
      n=split($3,a,","); out=""
      for(i=1;i<=n;i++){ s=a[i]; id=substr(s,1,length(s)-1); o=substr(s,length(s))
        t=seq[id]
        if(o=="-"){ r=""; for(j=length(t);j>0;j--){ c=substr(t,j,1)
            r=r (c=="A"?"T":c=="T"?"A":c=="C"?"G":c=="G"?"C":c=="a"?"t":c=="t"?"a":c=="c"?"g":c=="g"?"c":c) } t=r }
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
panph -i "$OUT/bad.gfa" -c "$OUT/cnb.bubbles.csv" -o "$OUT/badp" -q >/dev/null 2>&1
[ "$?" -ne 0 ] && ok "a path step with no link behind it is refused" \
               || bad "a malformed graph was accepted and silently repaired"

printf 'bubble_id,source,sink,path_support,min_inside_bp,max_inside_bp,inside_nodes\n9,900,902,2,10,10,"901"\n' \
  > "$OUT/ghost.csv"
panph -i "$OUT/cnb.sorted.gfa" -c "$OUT/ghost.csv" -o "$OUT/ghostp" -q >/dev/null 2>&1
[ "$?" -ne 0 ] && ok "a bubbles CSV naming nodes absent from the GFA is refused" \
               || bad "a mismatched CSV ran to completion"

panph -i "$OUT/cnb.sorted.gfa" -c "$OUT/cnb.bubbles.csv" -o "$OUT/nope" --bubble-id 4242 \
       -q >/dev/null 2>&1
[ "$?" -ne 0 ] && ok "a --bubble-id that names nothing is refused" \
               || bad "--bubble-id 4242 ran over an empty selection"

# Node ids alone do not tie a table to a graph (sorted graphs number nodes in order): a site whose interior
# length no haplotype spells is a table from another graph.
awk -F, -v OFS=, 'NR==1{for(i=1;i<=NF;i++)c[$i]=i; print; next} {$(c["min_inside_bp"])=99999; $(c["max_inside_bp"])=99999; print}' \
  "$OUT/cnb.bubbles.csv" > "$OUT/other.csv"
panph -i "$OUT/cnb.sorted.gfa" -c "$OUT/other.csv" -o "$OUT/otherp" -q > "$OUT/other.log" 2>&1
[ "$?" -ne 0 ] && grep -q "does not describe this graph" "$OUT/other.log" \
  && ok "a site table whose interiors no haplotype spells is refused" \
  || bad "a table from another graph was accepted: $(cat "$OUT/other.log")"
# Passes back to back are one span, as bubble measures them: a site every haplotype passes twice in a row,
# whose table gives the two passes' length together, describes this graph.
{ printf 'H\tVN:Z:1.0\n'
  printf "S\t1\t%s\nS\t2\t%s\nS\t3\t%s\nS\t4\t%s\nS\t5\t%s\n" "$BB1" "ACGTACGTAC" "${UNIT:0:60}" "TTGCATTGCA" "$BB2"
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t3\t+\t4\t+\t0M\nL\t4\t+\t2\t+\t0M\nL\t4\t+\t5\t+\t0M\n'
  printf 'P\th1\t1+,2+,3+,4+,2+,3+,4+,5+\t*\nP\th2\t1+,2+,3+,4+,2+,3+,4+,5+\t*\n'; } > "$OUT/b2b.gfa"
printf 'bubble_id,source,sink,path_support,min_inside_bp,max_inside_bp,inside_nodes\n1,2,4,2,140,140,"3"\n' > "$OUT/b2b.csv"
panph -i "$OUT/b2b.gfa" -c "$OUT/b2b.csv" -o "$OUT/b2bp" -q > "$OUT/b2b.log" 2>&1
[ "$?" -eq 0 ] && ok "a site passed twice back to back is measured as one span, as bubble does" \
  || bad "back-to-back passes refused: $(cat "$OUT/b2b.log")"
# An inversion site is in the table but not read: --bubble-id naming it says so.
awk -F, -v OFS=, 'NR==1{for(i=1;i<=NF;i++)c[$i]=i; print; next} {print} NR==2{$(c["bubble_id"])=77; $(c["kind"])="inversion"; print}' \
  "$OUT/cnb.bubbles.csv" > "$OUT/inv.csv"
panph -i "$OUT/cnb.sorted.gfa" -c "$OUT/inv.csv" -o "$OUT/invp" --bubble-id 77 -q > "$OUT/inv.log" 2>&1
[ "$?" -ne 0 ] && grep -q "names an inversion site" "$OUT/inv.log" \
  && ok "a --bubble-id naming an inversion site says so" \
  || bad "--bubble-id on an inversion site: $(cat "$OUT/inv.log")"
panph -i "$OUT/cnb.sorted.gfa" -c "$OUT/cnb.bubbles.csv" -o "$OUT/mu0" --min-unit-bp 0 -q > "$OUT/mu0.log" 2>&1
[ "$?" -ne 0 ] && grep -q "at least 1" "$OUT/mu0.log" \
  && ok "--min-unit-bp 0 is refused" || bad "--min-unit-bp 0: $(cat "$OUT/mu0.log")"
# A missing GTF is found before the fold, not after it.
panph -i "$OUT/cnb.sorted.gfa" -c "$OUT/cnb.bubbles.csv" -o "$OUT/nogtf" -r arr --gtf "$OUT/none.gtf" > "$OUT/nogtf.log" 2>&1
[ "$?" -ne 0 ] && ! grep -q "folded:" "$OUT/nogtf.log" \
  && ok "a missing --gtf is refused before any folding" || bad "missing GTF: $(cat "$OUT/nogtf.log")"

panph -i "$OUT/cnb.sorted.gfa" -c "$OUT/cnb.bubbles.csv" -o "$OUT/cnb" -q >/dev/null 2>&1
[ "$?" -ne 0 ] && ok "an output that would overwrite an input is refused" \
               || bad "an output was allowed to overwrite an input"

# ---------------------------------------------------------------- reference spelling
panph -i "$OUT/cnb.sorted.gfa" -c "$OUT/cnb.bubbles.csv" -o "$OUT/rx" -r cn3 \
       --resnarl-min-variant-bp 0 -q >/dev/null 2>&1
panph -i "$OUT/cnb.sorted.gfa" -c "$OUT/cnb.bubbles.csv" -o "$OUT/rX" -r CN3 \
       --resnarl-min-variant-bp 0 -q >/dev/null 2>&1
if [ -s "$OUT/rx.bubbles.csv" ] && [ -s "$OUT/rX.bubbles.csv" ]; then
  cmp -s "$OUT/rx.bubbles.csv" "$OUT/rX.bubbles.csv" \
    && ok "a case-insensitive reference alias gives the same re-snarled CSV" \
    || bad "-r CN3 and -r cn3 produced different call-ready CSVs"
else
  bad "no re-snarled CSV was produced for the reference-spelling check"
fi

# ---------------------------------------------------------------- determinism
panph -i "$OUT/cnb.sorted.gfa" -c "$OUT/cnb.bubbles.csv" -o "$OUT/t1" \
       --min-similarity 0.90 --threads 1 -q >/dev/null 2>&1
panph -i "$OUT/cnb.sorted.gfa" -c "$OUT/cnb.bubbles.csv" -o "$OUT/t8" \
       --min-similarity 0.90 --threads 8 -q >/dev/null 2>&1
cmp -s "$OUT/t1.normalized.gfa" "$OUT/t8.normalized.gfa" \
  && ok "approximate output is byte-identical at 1 and 8 threads" \
  || bad "approximate output differs between 1 and 8 threads"

# ---------------------------------------------------------------- approximate mode: direct alleles
# Approximate mode counts every crossing, direct ones included: counting only crossings with an interior
# node, one array carrier among nine direct alleles would read prevalence 1/1 and fold through a 50% gate.
panph -i "$OUT/prevb.sorted.gfa" -c "$OUT/prevb.bubbles.csv" -o "$OUT/apx" \
       --min-similarity 0.90 -q >/dev/null 2>&1
[ "$(rep "$OUT/apx.panphorte.report.tsv" n_traversing)" = "10" ] \
  && ok "approximate mode counts direct source->sink crossings in the denominator" \
  || bad "approximate n_traversing is $(rep "$OUT/apx.panphorte.report.tsv" n_traversing), expected 10"
[ "$(rep "$OUT/apx.panphorte.report.tsv" normalized)" = "no" ] \
  && ok "approximate mode refuses prevalence 1/10 at the default gate" \
  || bad "approximate mode folded a site carried by one haplotype in ten"

# ---------------------------------------------------------------- approximate mode: interruptions
# Two copies of a 16 bp unit separated by a 64 bp gap: the gap is 64/96 = 0.667 of the span, so the
# pair is one array only when --max-interruption-frac allows that much. The option applies to the copies
# found, not only to seeding the unit, and interruptions_bp reports the gap.
U16="GTTGTCTATGCCAGGG"
G64="TTTTTTTTGGGGCCCCTTTTTTTTGGGGCCCCTTTTTTTTGGGGCCCCTTTTTTTTGGGGCCCC"
{ printf 'H\tVN:Z:1.0\n'
  printf "S\t1\t%s\nS\t2\t%s\nS\t3\t%s\nS\t4\t%s\n" "$BB1" "$U16" "$G64" "$BB2"
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t2\t+\t0M\nL\t2\t+\t4\t+\t0M\n'
  printf 'L\t2\t+\t3\t+\t0M\nL\t3\t+\t2\t+\t0M\n'
  printf 'P\tclean\t1+,2+,2+,4+\t*\nP\tsep\t1+,2+,3+,2+,4+\t*\n'; } > "$OUT/gap.gfa"
"$BIN" bubble -i "$OUT/gap.gfa" -r clean -o "$OUT/gapb" --min-variant-bp 0 --merge-nearby-bp 0 -q >/dev/null 2>&1
panph -i "$OUT/gapb.sorted.gfa" -c "$OUT/gapb.bubbles.csv" -o "$OUT/g0" \
       --min-similarity 0.90 --min-unit-bp 16 --max-interruption-frac 0 --min-array-prevalence 0.5 \
       -q >/dev/null 2>&1
panph -i "$OUT/gapb.sorted.gfa" -c "$OUT/gapb.bubbles.csv" -o "$OUT/g9" \
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
# containing step range is replaced and those outside bases come back as fragment nodes around the unit
# step, so the haplotype spells exactly what it spelled before AND its copy is counted.
{ printf 'H\tVN:Z:1.0\n'
  printf "S\t1\t%s\nS\t2\t%s\nS\t3\tTTT%sGGG\nS\t4\t%s\n" "$BB1" "$UNIT" "$UNIT" "$BB2"
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t2\t+\t0M\nL\t2\t+\t4\t+\t0M\n'
  printf 'L\t1\t+\t3\t+\t0M\nL\t3\t+\t4\t+\t0M\n'
  printf 'P\tcleanA\t1+,2+,2+,4+\t*\nP\tcleanB\t1+,2+,2+,4+\t*\nP\tflank\t1+,3+,4+\t*\n'; } \
  > "$OUT/flank.gfa"
"$BIN" bubble -i "$OUT/flank.gfa" -r cleanA -o "$OUT/flankb" --min-variant-bp 0 --merge-nearby-bp 0 -q >/dev/null 2>&1
before_flank=$(spell "$OUT/flankb.sorted.gfa" flank)
panph -i "$OUT/flankb.sorted.gfa" -c "$OUT/flankb.bubbles.csv" -o "$OUT/flankp" \
       --min-similarity 0.90 --min-array-prevalence 0.5 -q >/dev/null 2>&1
[ "$(spell "$OUT/flankp.normalized.gfa" flank)" = "$before_flank" ] \
  && ok "a copy with a mid-node boundary is not truncated to the node edge" \
  || bad "the flanking haplotype lost $(( ${#before_flank} - $(spell "$OUT/flankp.normalized.gfa" flank | wc -c) )) bases"
# The copy sits inside a 70 bp node with 3 bp of flank on each side. Neither declining it (an undercount
# of a whole repeat unit) nor rounding it away (deleting the flanks) is necessary.
[ "$(rep "$OUT/flankp.panphorte.report.tsv" normalized)" = "yes" ] \
  && ok "a copy with flanking bases in its node still folds" \
  || bad "the copy was not folded, so its haplotype is undercounted by a whole unit"
for h in cleanA cleanB flank; do
  [ "$(spell "$OUT/flankb.sorted.gfa" $h)" = "$(spell "$OUT/flankp.normalized.gfa" $h)" ] \
    || bad "$h lost or gained bases when its copy was folded"
done
ok "the flanking bases survive as fragments, byte for byte"

# ---------------------------------------------------------------- approximate mode: shared fragments
# Fragments are interned by site, role, the step they follow and their bases, never minted per path.
# Minted per occurrence, two haplotypes that walk exactly the same route would each get a private copy of
# every flank fragment, so identical haplotypes would leave panphorte walking different routes and refine
# would read the site as a skeleton mismatch.
{ printf 'H\tVN:Z:1.0\n'
  printf "S\t1\t%s\nS\t2\t%s\nS\t3\tTTT%sGGG\nS\t4\t%s\n" "$BB1" "$UNIT" "$UNIT" "$BB2"
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t2\t+\t0M\nL\t2\t+\t4\t+\t0M\n'
  printf 'L\t1\t+\t3\t+\t0M\nL\t3\t+\t4\t+\t0M\n'
  printf 'P\tcleanA\t1+,2+,2+,4+\t*\nP\tcleanB\t1+,2+,2+,4+\t*\n'
  printf 'P\tflankA\t1+,3+,4+\t*\nP\tflankB\t1+,3+,4+\t*\n'; } > "$OUT/share.gfa"
"$BIN" bubble -i "$OUT/share.gfa" -r cleanA -o "$OUT/shareb" --min-variant-bp 0 --merge-nearby-bp 0 -q >/dev/null 2>&1
panph -i "$OUT/shareb.sorted.gfa" -c "$OUT/shareb.bubbles.csv" -o "$OUT/sharep" \
       --min-similarity 0.90 --min-array-prevalence 0.5 -q >/dev/null 2>&1
[ "$(rep "$OUT/sharep.panphorte.report.tsv" normalized)" = "yes" ] \
  && ok "the shared-fragment fixture folds" \
  || bad "the shared-fragment fixture did not fold, so the assertions below test nothing"
[ -n "$(pline "$OUT/sharep.normalized.gfa" flankA)" ] \
  && [ "$(pline "$OUT/sharep.normalized.gfa" flankA)" = "$(pline "$OUT/sharep.normalized.gfa" flankB)" ] \
  && ok "two identical haplotypes still walk one route after folding" \
  || bad "identical haplotypes walk different routes: $(pline "$OUT/sharep.normalized.gfa" flankA) vs $(pline "$OUT/sharep.normalized.gfa" flankB)"
share_frag=$(awk '$1 == "S" && ($3 == "TTT" || $3 == "GGG")' "$OUT/sharep.normalized.gfa" | wc -l | tr -d ' ')
[ "$share_frag" = "2" ] \
  && ok "each flank fragment is one node shared by both carriers" \
  || bad "$share_frag TTT/GGG fragment nodes, expected 2 (one per flank)"
for h in cleanA cleanB flankA flankB; do
  [ "$(spell "$OUT/shareb.sorted.gfa" $h)" = "$(spell "$OUT/sharep.normalized.gfa" $h)" ] \
    || bad "$h lost or gained bases when its fragments were shared"
done
ok "sharing fragments keeps every haplotype's bases"

# ---------------------------------------------------------------- no carrier is ever called CN 0
# The condition the whole boundary policy exists to prevent: a half-unit, half-literal site makes `call`
# read CN 0 for a haplotype that carries a copy, because copy number is counted from unit occurrences.
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
"$BIN" bubble -i "$OUT/e2e.gfa" -r cleanA -o "$OUT/e2eb" --min-variant-bp 0 --merge-nearby-bp 0 -q >/dev/null 2>&1
panph -i "$OUT/e2eb.sorted.gfa" -c "$OUT/e2eb.bubbles.csv" -o "$OUT/e2ep" \
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
# be left that no path visits: replaced nodes and their links left beside the unit route would be taken
# back into the site by the re-snarl.
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
  || bad "$e2e_orphan replaced node(s) left behind, so the old branch survives beside the unit route"
[ "$(rep "$OUT/e2ep.panphorte.report.tsv" nodes_collapsed)" != "0" ] \
  && ok "nodes_collapsed counts the nodes the rewrite replaced" \
  || bad "nodes_collapsed is 0 although the rewrite replaced nodes"
# The re-snarled site must be exactly the normalized route -- fragment, unit, fragment -- with no extra
# interior node standing for the branch the rewrite was supposed to remove.
e2e_inside=$(awk -F',' 'NR==2{print $6}' "$OUT/e2ep.bubbles.csv")
[ "$e2e_inside" = "3" ] \
  && ok "the re-snarled bubble holds only the normalized route (3 interior nodes)" \
  || bad "the re-snarled bubble has $e2e_inside interior nodes, so a dead branch was picked up"

# ---------------------------------------------------------------- overlapping sites are refused
# Two bubbles claiming the same interior describe one piece of sequence twice, and folding both rewrites
# one span inside the other -- the same array folded at two scales. The acceptance check catches it, but
# only after every haplotype has been aligned, so the input is rejected up front instead. `bubble` can
# emit such a pair: a site enclosing others whose array has the same unit as a site inside it.
{ printf 'bubble_id,source,source_orient,sink,sink_orient,inside_node_count,total_node_count,'
  printf 'path_support,distinct_alleles,ref_allele_support,alt_allele_support_max,'
  printf 'alt_allele_support_min,min_inside_bp,max_inside_bp,inside_nodes\n'
  printf '1,1,+,4,+,2,4,3,2,2,1,1,70,128,"2;5"\n'
  printf '2,1,+,4,+,2,4,3,2,2,1,1,64,70,"5;3"\n'; } > "$OUT/nested.bubbles.csv"
panph -i "$OUT/e2e.gfa" -c "$OUT/nested.bubbles.csv" -o "$OUT/nestedp" \
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
panph -i "$OUT/reuse.gfa" -c "$OUT/reuse.bubbles.csv" -o "$OUT/reusep" \
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
panph -i "$OUT/rev.gfa" -c "$OUT/reuse.bubbles.csv" -o "$OUT/revp" \
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
# direction. Keyed on the stored direction it would add again every adjacency the graph already holds
# as a dual, so a run that folds NOTHING would double the link set. No link may appear twice, once in
# each direction.
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
# case, because an arc inside the repeat unit is crossed once per copy and folding replaces all but the
# crossing that falls outside the folded span. What matters is whether the
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
panph -i "$OUT/twin.gfa" -c "$OUT/reuse.bubbles.csv" -o "$OUT/twinp" \
       --min-similarity 0.90 -q > "$OUT/twin.log" 2>&1
twin_rc=$?
# status of bubble <id> in a panphorte report
pstatus() { awk -F'\t' -v id="$2" 'NR==1{for(i=1;i<=NF;i++)c[$i]=i;next} $1==id{print $(c["status"])}' "$1"; }
# Such a site is declined, not the run refused: it is left unfolded and the others are folded.
[ "$twin_rc" = "0" ] && [ "$(pstatus "$OUT/twinp.panphorte.report.tsv" 1)" = "declined_shared_route" ] \
  && ok "a rewrite whose replaced route stays walkable declines that site (declined_shared_route)" \
  || bad "a site whose replaced route stays walkable: exit $twin_rc, status '$(pstatus "$OUT/twinp.panphorte.report.tsv" 1)'"
[ "$(awk -F'\t' '$1=="P" && $2=="cleanA"{print $3}' "$OUT/twinp.normalized.gfa")" = "1+,2+,2+,4+" ] \
  && ok "and the declined site is left exactly as it was" \
  || bad "the declined site was rewritten: cleanA is $(awk -F'\t' '$1=="P" && $2=="cleanA"{print $3}' "$OUT/twinp.normalized.gfa")"
# The other side of the same rule: the reuse fixture keeps node 2 and one of its arcs alive, so the
# obsolete-and-live count is non-zero, but the entry and exit arcs go and the old route is broken. It
# must run, not be refused -- refusing on the arc count alone would block every real array.
[ -s "$OUT/reusep.normalized.gfa" ] \
  && ok "a shared link whose route no longer reconnects does not refuse the run" \
  || bad "the reuse fixture was refused although its replaced route is broken"

# ...and the check has to be asked PER EDIT. A haplotype normalized at two sites has two independent old
# routes; testing the whole original walk conflates them, so the first site's route being properly gone
# makes the walk unwalkable and hides the second site's route surviving intact. Here site A's nodes are
# removed outright while `stubX` and `stubY` keep every node and arc of site B alive, so main, main2 and
# flankA are each broken at A and intact at B. The whole-walk form would accept this fixture and write a
# graph; the per-edit form rejects all three.
UB="TGGCTTTGTCATGCTAGACCACGGGCCTTTACAGGGGACTAGACCTTCTGCCGGCGAGAGGTCC"
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
panph -i "$OUT/twosite.gfa" -c "$OUT/twosite.bubbles.csv" -o "$OUT/twositep" \
       --min-similarity 0.90 -q > "$OUT/twosite.log" 2>&1
ts_rc=$?
[ "$ts_rc" = "0" ] && [ "$(pstatus "$OUT/twositep.panphorte.report.tsv" 2)" = "declined_shared_route" ] \
  && [ "$(pstatus "$OUT/twositep.panphorte.report.tsv" 1)" = "normalized" ] \
  && ok "a surviving route at the second of two sites declines that site alone, not masked by the first" \
  || bad "two sites: exit $ts_rc, statuses '$(pstatus "$OUT/twositep.panphorte.report.tsv" 1)' / '$(pstatus "$OUT/twositep.panphorte.report.tsv" 2)', expected normalized / declined_shared_route"

# ---------------------------------------------------------------- overlap refusal respects --bubble-id
# Only the SELECTED sites can overlap: refusing on a pair the run does not touch would block the one safe
# way to work at a locus whose CSV carries an overlap, which is to name one member of it.
for want in 1 2; do
  panph -i "$OUT/e2e.gfa" -c "$OUT/nested.bubbles.csv" -o "$OUT/sel$want" \
         --min-similarity 0.90 --bubble-id "$want" -q > "$OUT/sel$want.log" 2>&1
  [ "$?" = "0" ] \
    && ok "selecting overlapping bubble $want on its own runs" \
    || bad "--bubble-id $want was refused although nothing it selects overlaps: $(head -1 "$OUT/sel$want.log")"
done
panph -i "$OUT/e2e.gfa" -c "$OUT/nested.bubbles.csv" -o "$OUT/selboth" \
       --min-similarity 0.90 --bubble-id 1 --bubble-id 2 -q > "$OUT/selboth.log" 2>&1
[ "$?" != "0" ] && grep -q "both claim interior node" "$OUT/selboth.log" \
  && ok "selecting both members of an overlapping pair is still refused" \
  || bad "--bubble-id 1 --bubble-id 2 was accepted"

# ---------------------------------------------------------------- two copies inside one node
# Two copies of a 16 bp unit sit in a single 96 bp node, separated by 64 bp. There is no step boundary
# between them, so mapping each copy to its own containing step range declined the second -- at an
# array that is a whole repeat unit missing from the copy number. Copies sharing a range are emitted as
# one block: fragment, unit, fragment, unit, fragment.
U16="GTTGTCTATGCCAGGG"
G64="TTTTTTTTGGGGCCCCTTTTTTTTGGGGCCCCTTTTTTTTGGGGCCCCTTTTTTTTGGGGCCCC"
{ printf 'H\tVN:Z:1.0\n'
  printf "S\t1\t%s\nS\t2\t%s\nS\t3\t%s%s%s\nS\t4\t%s\n" "$BB1" "$U16" "$U16" "$G64" "$U16" "$BB2"
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t2\t+\t0M\nL\t2\t+\t4\t+\t0M\n'
  printf 'L\t1\t+\t3\t+\t0M\nL\t3\t+\t4\t+\t0M\n'
  printf 'P\tclean\t1+,2+,2+,4+\t*\nP\tsep\t1+,3+,4+\t*\n'; } > "$OUT/tio.gfa"
"$BIN" bubble -i "$OUT/tio.gfa" -r clean -o "$OUT/tiob" --min-variant-bp 0 --merge-nearby-bp 0 -q >/dev/null 2>&1
sep_before=$(spell "$OUT/tiob.sorted.gfa" sep)
panph -i "$OUT/tiob.sorted.gfa" -c "$OUT/tiob.bubbles.csv" -o "$OUT/tiop" \
       --min-similarity 0.90 --min-unit-bp 16 --max-interruption-frac 0.9 --min-array-prevalence 0.5 \
       -q >/dev/null 2>&1
[ "$(rep "$OUT/tiop.panphorte.report.tsv" normalized)" = "yes" ] \
  && ok "two copies sharing one node are both folded, neither declined" \
  || bad "the site with two copies in one node was not folded"
[ "$(rep "$OUT/tiop.panphorte.report.tsv" max_copies)" = "2" ] \
  && ok "the haplotype's copy number is 2, not 1" \
  || bad "max_copies is $(rep "$OUT/tiop.panphorte.report.tsv" max_copies), expected 2"
[ "$(spell "$OUT/tiop.normalized.gfa" sep)" = "$sep_before" ] \
  && ok "the 64 bp between the two copies survives as a fragment, byte for byte" \
  || bad "the intervening sequence was lost when both copies were folded"
# sep must traverse the unit twice with the gap fragment between: unit, frag, unit.
n=$(awk '$1=="P" && $2=="sep"{print $3}' "$OUT/tiop.normalized.gfa" | tr ',' '\n' | wc -l | tr -d ' ')
[ "$n" = "5" ] \
  && ok "sep is anchor, unit, fragment, unit, anchor (5 steps)" \
  || bad "sep has $n steps, expected 5"


# ---------------------------------------------------------------- the unit graph table is emitted
[ -s "$OUT/cnp.panphorte.unit_graph.tsv" ] \
  && ok "a unit graph table is written (exact mode too)" \
  || bad "no unit_graph.tsv was written"
head -1 "$OUT/cnp.panphorte.unit_graph.tsv" | grep -q "canonical_motif" \
  && ok "the table names the site motif each fold stands for" \
  || bad "unit graph header lacks canonical_motif"
# The join key has to name a node in the graph delivered BESIDE it. --reference-path renumbers, so the
# id a fold is created with is not the id it ships under, and the table would name nodes that do not exist.
prov_out=$(onenode "$OUT/e2ep.panphorte.unit_graph.tsv" | head -1)
[ -n "$prov_out" ] && grep -q "^S	$prov_out	" "$OUT/e2ep.normalized.sorted.gfa" \
  && ok "a one-node fold's node is a node of the delivered graph" \
  || bad "unit node '$prov_out' is not a segment of the sorted graph"
grep -q "^L	$prov_out	+	$prov_out	+" "$OUT/e2ep.normalized.sorted.gfa" \
  && ok "the delivered one-node fold carries the self-loop its copies walk" \
  || bad "node $prov_out has no self-loop, so it is not the fold"

# Sorting can FLIP a one-node fold, and then the sequence panphorte built the node with is not the sequence
# the delivered node spells; the table's unit is the unit as folded, read on either strand. Here the reference walks the site on the minus strand through the flanked single
# copy, so the unit is seeded from the carriers' forward copies (two of the four paths, meeting the
# prevalence floor) and --flip reverses the node. A reference walking the copies itself has the node built
# in its own orientation, and nothing is flipped.
{ printf 'H\tVN:Z:1.0\n'
  printf "S\t1\t%s\nS\t2\t%s\nS\t3\tTTT%sGGG\nS\t4\t%s\n" "$BB1" "$UNIT" "$UNIT" "$BB2"
  for e in "1 2" "2 2" "2 4" "1 3" "3 4"; do
    set -- $e; printf "L\t%s\t+\t%s\t+\t0M\n" "$1" "$2"
  done
  printf 'P\trefrev\t4-,3-,1-\t*\nP\tcarrier\t1+,2+,2+,4+\t*\nP\tcarrier2\t1+,2+,2+,2+,4+\t*\n'
  printf 'P\tflank\t1+,3+,4+\t*\n'; } \
  > "$OUT/flip.gfa"
panph -i "$OUT/flip.gfa" -c "$OUT/reuse.bubbles.csv" -o "$OUT/flipp" \
       --min-similarity 0.90 -r refrev -q >/dev/null 2>&1
fl_out=$(onenode "$OUT/flipp.panphorte.unit_graph.tsv" 2>/dev/null | head -1)
fl_made=$(rep "$OUT/flipp.panphorte.unit_graph.tsv" unit 2>/dev/null)
fl_seq=$(awk -v n="$fl_out" '$1=="S" && $2==n {print $3}' "$OUT/flipp.normalized.sorted.gfa")
fl_rc=$(printf '%s' "$fl_made" | rev | tr ACGT TGCA)
[ -n "$fl_seq" ] && [ "$fl_seq" = "$fl_rc" ] \
  && ok "the sort flips the fold, and the table's unit is the delivered node read on the other strand" \
  || bad "the delivered node $fl_out is not the table's unit reverse-complemented (fixture proves nothing)"

# ---------------------------------------------------------------- the GTF is an input like any other
# A GTF naming any output, not only the gene CSV, would be read and then overwritten by the commit. It is
# checked with and without --reference-path, since without it the annotation is skipped but the file is
# still an input.
# The output is made to exist, so the GTF is readable and the alias check is what refuses it.
: > "$OUT/gawith.panphorte.report.tsv"
panph -i "$OUT/e2e.gfa" -c "$OUT/e2eb.bubbles.csv" -o "$OUT/gawith" \
       --gtf "$OUT/gawith.panphorte.report.tsv" -r cleanA -q > "$OUT/gawith.log" 2>&1
ga_rc=$?
[ "$ga_rc" != "0" ] && grep -q "same file as input" "$OUT/gawith.log" \
  && ok "a GTF that names an output is refused (with --reference-path)" \
  || bad "GTF aliasing an output was accepted with --reference-path (exit $ga_rc)"
: > "$OUT/gawithout.panphorte.report.tsv"
panph -i "$OUT/e2e.gfa" -c "$OUT/e2eb.bubbles.csv" -o "$OUT/gawithout" \
       --gtf "$OUT/gawithout.panphorte.report.tsv" -q > "$OUT/gawithout.log" 2>&1
ga_rc=$?
[ "$ga_rc" != "0" ] && grep -q "same file as input" "$OUT/gawithout.log" \
  && ok "a GTF that names an output is refused with no --reference-path, where it is not even used" \
  || bad "GTF aliasing an output was accepted without --reference-path (exit $ga_rc)"

# ---------------------------------------------------------------- copies on unequal node boundaries
# `mixed` reaches the motif twice: once on a node that is exactly one unit, and once inside a node whose
# copy is flanked by 3 bp on each side. Both fold, so its copy number is 2. `alt` bypasses the array so
# that the first unit node is interior rather than the site's own source anchor; without it the bubble
# starts AT that node and each haplotype has a single copy inside, and the site does not seed.
{ printf 'H\tVN:Z:1.0\n'
  printf "S\t1\t%s\nS\t2\t%s\nS\t5\t%s\nS\t3\tTTT%sGGG\nS\t4\t%s\nS\t6\t%s\n" \
         "$BB1" "$UNIT" "$UNIT" "$UNIT" "$BB2" "$(head -c 200 /dev/zero | tr '\0' 'C')"
  for e in "1 2" "2 5" "5 4" "2 3" "3 4" "1 6" "6 4"; do
    set -- $e; printf "L\t%s\t+\t%s\t+\t0M\n" "$1" "$2"
  done
  printf 'P\tcleanA\t1+,2+,5+,4+\t*\nP\tcleanB\t1+,2+,5+,4+\t*\nP\tmixed\t1+,2+,3+,4+\t*\n'
  printf 'P\talt\t1+,6+,4+\t*\n'; } > "$OUT/mix.gfa"
"$BIN" bubble -i "$OUT/mix.gfa" -r cleanA -o "$OUT/mixb" --min-variant-bp 0 --merge-nearby-bp 0 -q >/dev/null 2>&1
mix_before=$(spell "$OUT/mixb.sorted.gfa" mixed)
panph -i "$OUT/mixb.sorted.gfa" -c "$OUT/mixb.bubbles.csv" -o "$OUT/mixp" \
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

# ---------------------------------------------------------------- --refine: flanks realigned around the array
# After folding, the bases around the array keep the input graph's representation: the same flank can
# be spelled through different node chains by different haplotypes. --refine realigns the flanks before
# and after each carrier's array block (one alignment per side), keeps the block and every non-carrier
# as folded, and accepts the rebuild only when it has fewer nodes and steps.
#
# Carriers A, A2, C and the reverse-crossing R spell the left flank ACGTAC as 5+6; B spells it 8+9.
# The right flank GGTT is 7 for A, A2, B and 10+11 for C and R. D crosses without copies.
rf_mk() {   # <name> <path lines "NAME STEPS"...>: the shared refine fixture graph
  local n=$1; shift; local paths=("$@")
  { printf 'H\tVN:Z:1.0\n'
    printf "S\t1\t%s\nS\t2\t%s\nS\t4\t%s\n" "$BB1" "$UNIT" "$BB2"
    printf "S\t5\tACG\nS\t6\tTAC\nS\t7\tGGTT\nS\t8\tAC\nS\t9\tGTAC\nS\t10\tGG\nS\t11\tTT\nS\t12\tTTTT\nS\t13\tCA\n"
    local e
    for e in "1 5" "5 6" "6 2" "2 2" "2 7" "7 4" "1 8" "8 9" "9 2" "2 10" "10 11" "11 4" "1 12" "12 4" \
             "2 13" "13 2" "5 2"; do
      printf "L\t%s\t+\t%s\t+\t0M\n" "${e% *}" "${e#* }"
    done
    for e in "${paths[@]}"; do printf 'P\t%s\t%s\t*\n' "${e%% *}" "${e#* }"; done; } > "$OUT/$n.gfa"
  "$BIN" bubble -i "$OUT/$n.gfa" -r A -o "$OUT/${n}b" --min-variant-bp 0 --merge-nearby-bp 0 -q >/dev/null 2>&1
}
rf_run() {  # <name> <out> [panphorte args...]
  local n=$1 o=$2; shift 2
  panph -i "$OUT/${n}b.sorted.gfa" -c "$OUT/${n}b.bubbles.csv" -r A --resnarl-min-variant-bp 0 \
         -o "$OUT/$o" -q "$@" >/dev/null 2>&1
}
rcol() { awk -F'\t' -v w="$2" 'NR==1{for(i=1;i<=NF;i++)c[$i]=i;next} NR==2{print $(c[w])}' "$1"; }
# unit occurrences per path: steps on nodes that carry a self-link
repcounts() { awk '$1=="L" && $2==$4 {self[$2]=1}
                   $1=="P" {n=split($3,a,","); c=0; for(i=1;i<=n;i++) if (substr(a[i],1,length(a[i])-1) in self) c++
                            print $2, c}' "$1" "$1" | sort; }
# interior node ids of a path, without the two anchors
inner() { pline "$1" "$2" | tr ',' '\n' | sed '1d;$d'; }
SPLIT=("A 1+,5+,6+,2+,2+,7+,4+" "A2 1+,5+,6+,2+,2+,7+,4+" "B 1+,8+,9+,2+,2+,2+,7+,4+"
       "C 1+,5+,6+,2+,2+,10+,11+,4+" "R 4-,11-,10-,2-,2-,6-,5-,1-" "D 1+,12+,4+")

rf_mk rf1 "${SPLIT[@]}"
rf_run rf1 rf1n
rf_run rf1 rf1r --refine
[ ! -e "$OUT/rf1n.panphorte.refine.tsv" ] \
  && ok "without --refine no refine table is written" \
  || bad "a run without --refine wrote the refine table"
[ "$(rcol "$OUT/rf1r.panphorte.refine.tsv" refine_status)" = "refined" ] \
  && ok "--refine rebuilds the flanks of a folded site whose flanks are split differently" \
  || bad "--refine status '$(rcol "$OUT/rf1r.panphorte.refine.tsv" refine_status)', expected refined"
rf_counts="$(rcol "$OUT/rf1r.panphorte.refine.tsv" nodes_before)>$(rcol "$OUT/rf1r.panphorte.refine.tsv" nodes_after) $(rcol "$OUT/rf1r.panphorte.refine.tsv" steps_before)>$(rcol "$OUT/rf1r.panphorte.refine.tsv" steps_after)"
[ "$rf_counts" = "8>3 18>11" ] \
  && ok "flank nodes 8 -> 3 and steps 18 -> 11 (one node per flank, the non-carrier untouched)" \
  || bad "flank nodes/steps before>after '$rf_counts', expected '8>3 18>11'"
rf_ok=1
for h in A A2 B C R D; do
  [ "$(spell "$OUT/rf1n.normalized.sorted.gfa" $h)" = "$(spell "$OUT/rf1r.normalized.sorted.gfa" $h)" ] || rf_ok=0
done
[ "$rf_ok" = 1 ] && ok "every haplotype spells exactly what it spelled after folding" \
                 || bad "--refine changed a haplotype's sequence"
[ "$(repcounts "$OUT/rf1n.normalized.sorted.gfa")" = "$(repcounts "$OUT/rf1r.normalized.sorted.gfa")" ] \
  && ok "every haplotype walks the unit node the same number of times" \
  || bad "--refine changed a unit copy count"
# B spelled the left flank through other nodes than A; after --refine both walk one node for it.
rf_l_a=$(inner "$OUT/rf1r.normalized.sorted.gfa" A | head -1)
rf_l_b=$(inner "$OUT/rf1r.normalized.sorted.gfa" B | head -1)
rf_r_a=$(inner "$OUT/rf1r.normalized.sorted.gfa" A | tail -1)
rf_r_c=$(inner "$OUT/rf1r.normalized.sorted.gfa" C | tail -1)
[ "$rf_l_a" = "$rf_l_b" ] && [ "$rf_r_a" = "$rf_r_c" ] \
  && ok "carriers that spelled a flank through different chains now share its node" \
  || bad "flanks still scattered: left A=$rf_l_a B=$rf_l_b, right A=$rf_r_a C=$rf_r_c"
# R crosses the site in reverse: its route must be C's route reversed, not a reverse-complemented copy.
rf_c=$(pline "$OUT/rf1r.normalized.sorted.gfa" C)
rf_r=$(pline "$OUT/rf1r.normalized.sorted.gfa" R | tr ',' '\n' | tail -r 2>/dev/null || pline "$OUT/rf1r.normalized.sorted.gfa" R | tr ',' '\n' | tac)
rf_r=$(echo "$rf_r" | sed 's/+$/#/;s/-$/+/;s/#$/-/' | paste -sd, -)
[ "$rf_c" = "$rf_r" ] \
  && ok "a carrier crossing in reverse walks the same route as its forward twin" \
  || bad "reverse carrier route '$rf_r' differs from C '$rf_c'"
# The fixture's CA node is used only by the interruption fixture below, so it is an orphan in this graph
# with or without --refine; what --refine must not do is add one.
rf_topo() { awk '$1 == "S" {seg[$2] = 1}
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
          printf "%d %d %d\n", absent + 0, dangling + 0, orphan + 0}' "$1"; }
rf_topo_r=$(rf_topo "$OUT/rf1r.normalized.sorted.gfa")
rf_orph_n=$(rf_topo "$OUT/rf1n.normalized.sorted.gfa" | cut -d' ' -f3)
[ "$rf_topo_r" = "0 0 $rf_orph_n" ] \
  && ok "the refined graph has no absent node or dangling adjacency, and no new orphan segment" \
  || bad "refined graph topology (absent dangling orphan): $rf_topo_r; unrefined orphans: $rf_orph_n"

# A folded array may hold fragments between its copies, repeated: unit CA unit CA unit. They belong to the
# array, which --refine keeps verbatim; only an unfolded repeat in a FLANK is a duplication to refuse.
rf_mk rf3 "${SPLIT[@]}" "I 1+,8+,9+,2+,13+,2+,13+,2+,10+,11+,4+"
rf_run rf3 rf3n --min-similarity 0.90
rf_run rf3 rf3r --min-similarity 0.90 --refine
[ "$(rcol "$OUT/rf3r.panphorte.refine.tsv" refine_status)" = "refined" ] \
  && ok "a fragment repeated between folded copies does not block --refine" \
  || bad "status '$(rcol "$OUT/rf3r.panphorte.refine.tsv" refine_status)' with a repeated between-copy fragment"
[ "$(spell "$OUT/rf3n.normalized.sorted.gfa" I)" = "$(spell "$OUT/rf3r.normalized.sorted.gfa" I)" ] \
  && [ "$(repcounts "$OUT/rf3n.normalized.sorted.gfa")" = "$(repcounts "$OUT/rf3r.normalized.sorted.gfa")" ] \
  && ok "and the array block, fragments included, is kept as folded" \
  || bad "the array block with its between-copy fragments changed"

# A flank duplication: I walks the flank 8,9 twice before its copies. That is not the site's motif, so
# nothing folds it, and it is the duplication --refine must not linearize: the site is left exactly as
# folded. (This used I's interrupted copies 2,13,2,13,2 until exact mode learnt to top those up as
# single copies -- they are folded now, and read CN 3 instead of 0.)
rf_mk rf6 "${SPLIT[@]}" "I 1+,8+,9+,8+,9+,2+,2+,7+,4+"
printf 'L\t9\t+\t8\t+\t0M\n' >> "$OUT/rf6.gfa"
"$BIN" bubble -i "$OUT/rf6.gfa" -r A -o "$OUT/rf6b" --min-variant-bp 0 --merge-nearby-bp 0 -q >/dev/null 2>&1
rf_run rf6 rf2n
rf_run rf6 rf2r --refine
[ "$(rcol "$OUT/rf2r.panphorte.refine.tsv" refine_status)" = "skipped:unfolded_dup" ] \
  && ok "a site with an unfolded duplication in a flank is skipped" \
  || bad "status '$(rcol "$OUT/rf2r.panphorte.refine.tsv" refine_status)', expected skipped:unfolded_dup"
# The skip is found mid-scan; the row still counts the 7 crossing haplotypes and the 6 carrying copies.
[ "$(rcol "$OUT/rf2r.panphorte.refine.tsv" paths):$(rcol "$OUT/rf2r.panphorte.refine.tsv" carriers)" = "7:6" ] \
  && ok "and the skipped row reports its crossing haplotypes and carriers (7, 6)" \
  || bad "skipped row paths:carriers is $(rcol "$OUT/rf2r.panphorte.refine.tsv" paths):$(rcol "$OUT/rf2r.panphorte.refine.tsv" carriers), expected 7:6"
rf_same=1
for f in normalized.sorted.gfa bubbles.csv bandage_nodes.csv panphorte.report.tsv panphorte.unit_graph.tsv; do
  cmp -s "$OUT/rf2n.$f" "$OUT/rf2r.$f" || rf_same=0
done
[ "$rf_same" = 1 ] && ok "and every output is byte-identical to the run without --refine" \
                   || bad "a skipped site still changed the output"

# Flanks already shared: realignment cannot make the site simpler, so it is kept as folded.
rf_mk rf4 "A 1+,5+,2+,2+,7+,4+" "B 1+,5+,2+,2+,2+,7+,4+" "D 1+,12+,4+"
rf_run rf4 rf4n
rf_run rf4 rf4r --refine
[ "$(rcol "$OUT/rf4r.panphorte.refine.tsv" refine_status)" = "kept_original:not_simpler" ] \
  && cmp -s "$OUT/rf4n.normalized.sorted.gfa" "$OUT/rf4r.normalized.sorted.gfa" \
  && ok "a rebuild that is not simpler is not kept, and the graph is unchanged" \
  || bad "status '$(rcol "$OUT/rf4r.panphorte.refine.tsv" refine_status)' or the graph changed for a site already minimal"

# Sites are aligned in reference orientation and the reference must cross each one.
panph -i "$OUT/rf1b.sorted.gfa" -c "$OUT/rf1b.bubbles.csv" -o "$OUT/rf5" --refine -q \
       >/dev/null 2>&1
rf5_rc=$?
[ "$rf5_rc" -ne 0 ] && [ -z "$(ls "$OUT" | grep '^rf5\.')" ] \
  && ok "--refine without --reference-path is refused before anything is written" \
  || bad "--refine without --reference-path exited $rf5_rc or wrote: $(ls "$OUT" | grep '^rf5\.' | tr '\n' ' ')"

# ---------------------------------------------------------------- approximate unit phased at the array edge
# A tandem repeat has no natural start, and a base-level seed is cut wherever its probe lands: here 15 bases into a copy.
# Folded in that phase, k copies would become a head fragment, k-1 unit steps and a tail fragment, so every haplotype
# would lose one copy and a one-copy haplotype would get no unit at all. The unit is rotated to the array's edge when
# that folds more copies. Each interior is ONE node ('A' + k copies), so no node-step period exists and
# the base-level seed is what decides the phase.
PU="GAGCCCACCTCAGAGCCCGCCCCCAGCCCGACCACCCCG"
{ printf 'H\tVN:Z:1.0\n'
  printf "S\t1\t%s\nS\t2\tA%s\nS\t3\tA%s%s\nS\t4\tA%s%s%s\nS\t5\t%s\n" "$BB1" "$PU" "$PU" "$PU" "$PU" "$PU" "$PU" "$BB2"
  for e in "1 2" "1 3" "1 4" "2 5" "3 5" "4 5"; do printf "L\t%s\t+\t%s\t+\t0M\n" "${e% *}" "${e#* }"; done
  printf 'P\tk2a\t1+,3+,5+\t*\nP\tk2b\t1+,3+,5+\t*\nP\tk3a\t1+,4+,5+\t*\nP\tk3b\t1+,4+,5+\t*\n'
  printf 'P\tk1a\t1+,2+,5+\t*\nP\tk1b\t1+,2+,5+\t*\n'; } > "$OUT/phase.gfa"
"$BIN" bubble -i "$OUT/phase.gfa" -r k2a -o "$OUT/phaseb" --min-variant-bp 0 --merge-nearby-bp 0 -q >/dev/null 2>&1
panph -i "$OUT/phaseb.sorted.gfa" -c "$OUT/phaseb.bubbles.csv" -o "$OUT/phasep" --min-similarity 0.90 \
       --min-unit-bp 20 --min-array-prevalence 0.5 > "$OUT/phasep.log" 2>&1
phase_rep=$(awk '$1=="L" && $2==$4 {print $2}' "$OUT/phasep.normalized.gfa" | head -1)
if [ -z "$phase_rep" ]; then
  bad "the phase fixture folded nothing ($(tail -1 "$OUT/phasep.log"))"
else
  [ "$(awk -v r="$phase_rep" '$1=="S" && $2==r {print $3}' "$OUT/phasep.normalized.gfa")" = "$PU" ] \
    && ok "the folded unit is in the copies' own phase (rotated at the array edge)" \
    || bad "the folded unit is a rotation of the copies: $(awk -v r="$phase_rep" '$1=="S" && $2==r {print $3}' "$OUT/phasep.normalized.gfa")"
  phase_ok=1
  for h in k1a:1 k1b:1 k2a:2 k2b:2 k3a:3 k3b:3; do
    n=$(pline "$OUT/phasep.normalized.gfa" "${h%%:*}" | tr ',' '\n' | grep -c "^${phase_rep}[+-]$")
    [ "$n" = "${h##*:}" ] || { phase_ok=0; bad "${h%%:*} walks the unit node $n times, expected ${h##*:}"; }
  done
  [ "$phase_ok" = 1 ] && ok "every copy folds: one, two and three copies walk the unit node 1, 2 and 3 times"
  phase_seq=1
  for h in k1a k2a k3a; do
    [ "$(spell "$OUT/phaseb.sorted.gfa" $h)" = "$(spell "$OUT/phasep.normalized.gfa" $h)" ] || phase_seq=0
  done
  [ "$phase_seq" = 1 ] && ok "and exact copies keep every haplotype's sequence" || bad "re-phasing changed a haplotype's sequence"
  grep -q "rotated to start at the array edge" "$OUT/phasep.log" \
    && ok "the re-phasing is logged" || bad "the re-phasing was not logged"
fi

# ---------------------------------------------------------------- family commit is a transaction
# Panphorte writes a family of files (normalized GFA, report, unit graph table, copies, bubbles, Bandage).
# The shared StagedOutputs commit makes that atomic; this checks it on panphorte's own family.
panph -i "$OUT/prevb.sorted.gfa" -c "$OUT/prevb.bubbles.csv" -o "$OUT/txnp" -q >/dev/null 2>&1
if [ ! -s "$OUT/txnp.normalized.gfa" ]; then
  bad "the panphorte rollback fixture produced no family to protect"
else
  tp_n=$(ls "$OUT" | grep -c '^txnp\.')
  tp_gfa=$(cksum < "$OUT/txnp.normalized.gfa")
  tp_rep=$(cksum < "$OUT/txnp.panphorte.report.tsv")
  PANVAR_TEST_FAIL_COMMIT_AT=3 panph -i "$OUT/prevb.sorted.gfa" \
      -c "$OUT/prevb.bubbles.csv" -o "$OUT/txnp" -q >/dev/null 2>&1
  [ "$?" -ne 0 ] && ok "an injected failure in the panphorte family exits non-zero" \
                 || bad "an injected failure in the panphorte family exited 0"
  [ "$(ls "$OUT" | grep -c '^txnp\.')" = "$tp_n" ] \
    && ok "every file of the previous panphorte family is still present ($tp_n)" \
    || bad "the panphorte family lost a file: $(ls "$OUT" | grep '^txnp\.' | tr '\n' ' ')"
  [ "$(cksum < "$OUT/txnp.normalized.gfa")" = "$tp_gfa" ] \
    && [ "$(cksum < "$OUT/txnp.panphorte.report.tsv")" = "$tp_rep" ] \
    && ok "and the graph and report are byte-unchanged" \
    || bad "a panphorte output was replaced by a run that did not complete"
  # The set-aside window, which the bubble fixture reaches but this family had never exercised.
  PANVAR_TEST_FAIL_COMMIT_AFTER_SETASIDE=3 panph -i "$OUT/prevb.sorted.gfa" \
      -c "$OUT/prevb.bubbles.csv" -o "$OUT/txnp" -q >/dev/null 2>&1
  [ "$(cksum < "$OUT/txnp.normalized.gfa" 2>/dev/null)" = "$tp_gfa" ] \
    && ok "and also when the failure lands between set-aside and install" \
    || bad "the set-aside window lost a panphorte output"
  [ "$(ls "$OUT" | grep -c -- '-tmp\.\|-prev\.')" = "0" ] \
    && ok "panphorte leaves no staging or reserve residue" \
    || bad "panphorte left staging/reserve residue"
fi

# ---------------------------------------------------------------- exact mode: P1-P4
# Synthetic sites for exact folding (--min-similarity 1.0). U is a 60 bp unit (the minimum is 50).
"${PYTHON:-python3}" - "$OUT" <<'PYEOF'
import random, sys
out = sys.argv[1]
r = random.Random(17)
rnd = lambda n: ''.join(r.choice('ACGT') for _ in range(n))
rc = lambda s: s.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
L, R, U, U2 = rnd(80), rnd(80), rnd(60), rnd(60)
def write(name, nodes, edges, paths):
    with open(f'{out}/{name}.gfa', 'w') as o:
        o.write('H\tVN:Z:1.0\n')
        for k, v in nodes.items(): o.write(f'S\t{k}\t{v}\n')
        for a, b in edges: o.write(f'L\t{a[:-1]}\t{a[-1]}\t{b[:-1]}\t{b[-1]}\t0M\n')
        for k, v in paths.items(): o.write(f'P\t{k}\t{v}\t*\n')
# P1: the sink boundary spells the unit too: a copy every haplotype carries, so the site is widened over it
write('p1', {'1': L, '2': U, '3': U, '4': R}, [('1+','2+'),('2+','2+'),('2+','3+'),('1+','3+'),('3+','4+')],
      {'A': '1+,2+,2+,3+,4+', 'B': '1+,2+,2+,2+,3+,4+', 'C': '1+,3+,4+'})
# P2: a one-copy-each haplotype carrying the unit then its reverse complement
write('p2', {'1': L, '2': U, '3': R, '4': rc(U)}, [('1+','2+'),('2+','2+'),('2+','3+'),('2+','4+'),('4+','3+'),('1+','3+')],
      {'A': '1+,2+,2+,3+', 'B': '1+,2+,2+,3+', 'C': '1+,2+,2+,3+', 'T': '1+,2+,4+,3+', 'D': '1+,3+'})
# P3: one copy of the site's motif next to an array of an unrelated motif carried by nobody else
write('p3', {'1': L, '2': U, '3': R, '5': U2}, [('1+','2+'),('2+','2+'),('2+','3+'),('2+','5+'),('5+','5+'),('5+','3+'),('1+','3+')],
      {'A': '1+,2+,2+,3+', 'B': '1+,2+,2+,3+', 'C': '1+,2+,2+,3+', 'Q': '1+,2+,3+', 'P': '1+,2+,5+,5+,3+', 'D': '1+,3+'})
# P4: the folded unit is two copies of the motif (node 4 = UU)
write('p4', {'1': L, '4': U + U, '3': R}, [('1+','4+'),('4+','4+'),('4+','3+'),('1+','3+')],
      {'A': '1+,4+,4+,3+', 'B': '1+,4+,4+,4+,3+', 'D': '1+,3+'})
PYEOF
# Unit node visits per path, over the one-node folds panphorte reports (output ids)
repvisits() { awk -F'\t' 'FNR==NR{rep[$1]=1; next} $1=="P"{n=split($3,a,","); c=0
                for(i=1;i<=n;i++) if (substr(a[i],1,length(a[i])-1) in rep) c++; printf "%s %d ", $2, c}' \
                <(onenode "$1.panphorte.unit_graph.tsv") "$1.normalized.sorted.gfa"; }
for n in p1 p2 p3 p4; do
  "$BIN" bubble -i "$OUT/$n.gfa" -r A -o "$OUT/${n}b" --min-variant-bp 0 --merge-nearby-bp 0 -q >/dev/null 2>&1
  panph -i "$OUT/${n}b.sorted.gfa" -c "$OUT/${n}b.bubbles.csv" -r A --min-similarity 1.0 \
     --resnarl-min-variant-bp 0 -o "$OUT/${n}p" -q > "$OUT/${n}p.log" 2>&1
  echo $? > "$OUT/${n}p.rc"
done
[ "$(cat "$OUT/p1p.rc")" = 0 ] && [ "$(repvisits "$OUT/p1p")" = "A 3 B 4 C 1 " ] \
  && ok "exact mode counts a sink spelling the unit as the copy every haplotype shares (A 3, B 4, C 1)" \
  || bad "exact mode with a unit-spelling sink: rc $(cat "$OUT/p1p.rc"), unit visits '$(repvisits "$OUT/p1p" 2>/dev/null)'"
[ "$(cat "$OUT/p2p.rc")" = 0 ] && repvisits "$OUT/p2p" | grep -q "T 2 " \
  && ok "a unit and its reverse complement are folded as two copies, losslessly (T 2)" \
  || bad "unit + reverse complement: rc $(cat "$OUT/p2p.rc") ($(tail -1 "$OUT/p2p.log")), unit visits '$(repvisits "$OUT/p2p" 2>/dev/null)'"
p3c=$(repvisits "$OUT/p3p" | tr ' ' '\n' | awk 'prev=="P"{print; exit} {prev=$0}')
p3n=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)c[$i]=i;next} {print $(c["paths_normalized"])}' "$OUT/p3p.panphorte.report.tsv" | head -1)
[ "$p3c" = "1" ] && [ "$p3n" = "5" ] \
  && ok "a lone copy beside an unrelated array is folded (P reads 1, not 0); 5 paths normalized" \
  || bad "lone copy beside an unrelated array: P count '$p3c', paths_normalized '$p3n' (expected 1, 5)"
q4=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)c[$i]=i;next} {print $(c["unit_bp"])":"$(c["copy_quantum"])}' "$OUT/p4p.panphorte.unit_graph.tsv" | head -1)
# A node holding two whole 60 bp copies is folded as the 60 bp unit walked twice, not as a 120 bp unit
# walked once (which read CN 1 for two copies): A carries 2 such nodes, 4 copies; B 3 nodes, 6.
[ "$q4" = "60:1" ] && [ "$(repvisits "$OUT/p4p")" = "A 4 B 6 D 0 " ] \
  && ok "a node of two whole copies is the unit walked twice (60 bp, A 4, B 6)" \
  || bad "unit_bp:copy_quantum is '$q4' (expected 60:1), unit visits '$(repvisits "$OUT/p4p" 2>/dev/null)' (expected A 4 B 6 D 0)"

# ---------------------------------------------------------------- crossings, units and counts: Q1-Q6
# Each fixture failed on the binary before the per-crossing rework (bcc134e): Q1 aborted at 0.95 and
# lost a copy at 1.0, Q2 read the short carrier 0, Q3 folded the sub-repeat, Q4 dropped the lone copy,
# Q5 reported 0 interruption bases and Q6 reported max_copies 2.
"${PYTHON:-python3}" - "$OUT" <<'PYEOF'
import random, sys
out = sys.argv[1]
def rnd(n, seed):
    r = random.Random(seed)
    return ''.join(r.choice('ACGT') for _ in range(n))
def write(name, nodes, paths):
    edges = set()
    def key(a, ao, b, bo):
        f = (a, ao, b, bo); r = (b, '-' if bo == '+' else '+', a, '-' if ao == '+' else '+')
        return min(f, r)
    for w in paths.values():
        st = [(x[:-1], x[-1]) for x in w.split(',')]
        for (a, ao), (b, bo) in zip(st, st[1:]): edges.add(key(a, ao, b, bo))
    with open(f'{out}/{name}.gfa', 'w') as o:
        o.write('H\tVN:Z:1.0\n')
        for k, v in nodes.items(): o.write(f'S\t{k}\t{v}\n')
        for a, ao, b, bo in sorted(edges): o.write(f'L\t{a}\t{ao}\t{b}\t{bo}\t0M\n')
        for k, v in paths.items(): o.write(f'P\t{k}\t{v}\t*\n')
# Q1: 'dup' crosses site 1 (anchors 1/3) twice, once with two copies and once with one
write('q1', {'0': rnd(3000, 1), '1': rnd(500, 2), '2': rnd(100, 3), '3': rnd(500, 4), '5': rnd(500, 5),
             '6': rnd(120, 6), '7': rnd(500, 7), '9': rnd(3000, 8)},
      {'a1': '0+,1+,2+,2+,3+,5+,6+,6+,7+,9+', 'a2': '0+,1+,2+,2+,2+,3+,5+,6+,6+,6+,7+,9+',
       'a3': '0+,1+,2+,2+,3+,5+,6+,6+,7+,9+', 'dup': '0+,1+,2+,2+,3+,5+,6+,6+,7+,1+,2+,3+,9+'})
with open(f'{out}/q1.csv', 'w') as o:
    o.write('bubble_id,source,source_orient,sink,sink_orient,inside_node_count,total_node_count,path_support,'
            'distinct_alleles,ref_allele_support,alt_allele_support_max,alt_allele_support_min,min_inside_bp,'
            'max_inside_bp,inside_nodes\n1,1,+,3,+,1,3,4,2,1,1,1,0,0,"2"\n2,5,+,7,+,1,3,4,2,1,1,1,0,0,"6"\n')
# Q2: a one-copy carrier whose copy has a 3 bp deletion (99% identity), so it is shorter than the unit
U = rnd(300, 71)
write('q2', {'1': rnd(3000, 1), '2': U, '3': U[:150] + U[153:], '4': rnd(3000, 2)},
      {'h3': '1+,2+,2+,2+,4+', 'h2': '1+,2+,2+,4+', 'h2b': '1+,2+,2+,4+', 'h1d': '1+,3+,4+', 'h2d': '1+,3+,2+,4+',
       'h0': '1+,4+'})
# Q3: an 860 bp module (m1 s s s m2) in 2/1/0 copies; the 60 bp s repeat inside it is in every module copy
M = 'm1+,s+,s+,s+,m2+'
q3 = {f'two{i}': f'1+,{M},{M},9+' for i in range(3)}
q3.update({f'one{i}': f'1+,{M},9+' for i in range(3)})
q3['zero'] = '1+,9+'
write('q3', {'1': rnd(40000, 1), 'm1': rnd(400, 2), 's': rnd(60, 3), 'm2': rnd(400, 4), '9': rnd(40000, 5)}, q3)
# Q4/Q6: a 64 bp unit with a 200 bp spacer; a lone copy after the spacer is still a copy
Uj, Sp = rnd(64, 61), rnd(200, 62)
q4n = {'1': rnd(3000, 1), '2': Uj, '3': Sp, '4': rnd(3000, 2)}
write('q4', q4n, {'h3': '1+,2+,2+,3+,2+,4+', 'h3b': '1+,2+,2+,3+,2+,4+', 'h2': '1+,2+,2+,3+,4+', 'h0': '1+,3+,4+'})
write('q6', q4n, {'h4': '1+,2+,2+,3+,2+,2+,4+', 'h4b': '1+,2+,2+,3+,2+,2+,4+', 'h2': '1+,2+,2+,3+,4+',
                  'h0': '1+,3+,4+'})
# Q5: a 40 bp spacer inside the node that also holds three 100 bp copies (U U Sp U)
Uc, Sc = rnd(100, 31), rnd(40, 32)
write('q5', {'1': rnd(2000, 1), '2': Uc + Uc + Sc + Uc, '3': Uc + Uc, '5': Uc, '4': rnd(2000, 2)},
      {'h3': '1+,2+,4+', 'h3b': '1+,2+,4+', 'h2': '1+,3+,4+', 'h1': '1+,5+,4+', 'h0': '1+,4+'})
PYEOF
# Unit node visits per path in a run without --reference-path (the graph is written unsorted)
qvisits() { awk -F'\t' 'FNR==NR{rep[$1]=1; next} $1=="P"{n=split($3,a,","); c=0
              for(i=1;i<=n;i++) if (substr(a[i],1,length(a[i])-1) in rep) c++; printf "%s %d ", $2, c}' \
              <(onenode "$1.panphorte.unit_graph.tsv") "$1.normalized.gfa"; }
qfold() { local n=$1 tag=$2; shift 2
          panph -i "$OUT/${n}b.sorted.gfa" -c "$OUT/${n}b.bubbles.csv" -o "$OUT/$n$tag" -q "$@" \
              >/dev/null 2>"$OUT/$n$tag.log"; echo $? > "$OUT/$n$tag.rc"; }
qbub() { "$BIN" bubble -i "$OUT/$1.gfa" -r "$2" -o "$OUT/${1}b" --min-variant-bp 0 --merge-nearby-bp 0 -q >/dev/null 2>&1; }

# Q1: a haplotype crossing a site twice has each crossing folded on its own: 2 + 1 copies at site 1
# and 2 at site 2. Joining the crossings made one array spanning both, which overlapped site 2's edit.
for s in 0.95 1.0; do
  panph -i "$OUT/q1.gfa" -c "$OUT/q1.csv" -o "$OUT/q1_$s" --min-similarity "$s" -q \
      >/dev/null 2>"$OUT/q1_$s.log"; echo $? > "$OUT/q1_$s.rc"
  [ "$(cat "$OUT/q1_$s.rc")" = 0 ] && [ "$(qvisits "$OUT/q1_$s")" = "a1 4 a2 6 a3 4 dup 5 " ] \
    && ok "a double crossing folds per crossing at --min-similarity $s (dup 2+1+2 = 5)" \
    || bad "double crossing at $s: rc $(cat "$OUT/q1_$s.rc") ($(tail -1 "$OUT/q1_$s.log")), unit visits '$(qvisits "$OUT/q1_$s" 2>/dev/null)'"
done
# Q2: one copy carrying a 3 bp deletion is shorter than the unit, but still within the alignment band
qbub q2 h3; qfold q2 p --min-similarity 0.95
[ "$(qvisits "$OUT/q2p")" = "h3 3 h2 2 h2b 2 h1d 1 h2d 2 h0 0 " ] \
  && ok "a copy a few bases shorter than the unit is counted (h1d 1)" \
  || bad "short copy: unit visits '$(qvisits "$OUT/q2p" 2>/dev/null)'"
# Q3: a module carried by 6 of 7 haplotypes contains a 60 bp repeat present in every module copy. The
# module, not its sub-repeat, is the unit whose count varies: 2/1/0, not 6/3/0. At the default 0.5 gate
# the module's 3/7 carriers of two copies fall short and the site is left alone.
qbub q3 two0; qfold q3 p4 --min-similarity 0.95 --min-array-prevalence 0.4
[ "$(rep "$OUT/q3p4.panphorte.report.tsv" unit_bp)" = 980 ] \
  && [ "$(qvisits "$OUT/q3p4")" = "two0 2 two1 2 two2 2 one0 1 one1 1 one2 1 zero 0 " ] \
  && ok "a module containing a sub-repeat is folded as the module (980 bp, 2/1/0)" \
  || bad "module vs sub-repeat: unit_bp $(rep "$OUT/q3p4.panphorte.report.tsv" unit_bp), unit visits '$(qvisits "$OUT/q3p4" 2>/dev/null)'"
qfold q3 p5 --min-similarity 0.95
[ "$(rep "$OUT/q3p5.panphorte.report.tsv" status)" = below_prevalence ] \
  && ok "and at the default gate the module site is below prevalence, not folded as its sub-repeat" \
  || bad "module site at the default gate: status $(rep "$OUT/q3p5.panphorte.report.tsv" status)"
# Q4: in exact mode a lone copy after a 200 bp spacer is still a copy of the site's unit (h3 reads 3)
qbub q4 h3; qfold q4 ex --min-similarity 1.0
[ "$(qvisits "$OUT/q4ex")" = "h3 3 h3b 3 h2 2 h0 0 " ] \
  && ok "exact mode folds a lone copy after a spacer (h3 3)" \
  || bad "exact lone copy after a spacer: unit visits '$(qvisits "$OUT/q4ex" 2>/dev/null)'"
# Q5: a 40 bp spacer inside the node that also holds the copies is an interruption: 2 carriers x 40 bp
qbub q5 h3; qfold q5 p --min-similarity 0.95 --max-interruption-frac 0.25
[ "$(rep "$OUT/q5p.panphorte.report.tsv" interruptions_bp)" = 80 ] \
  && [ "$(qvisits "$OUT/q5p")" = "h3 3 h3b 3 h2 2 h1 1 h0 0 " ] \
  && ok "a spacer inside a copy-bearing node counts as interruption bases (80)" \
  || bad "in-node spacer: interruptions_bp $(rep "$OUT/q5p.panphorte.report.tsv" interruptions_bp), unit visits '$(qvisits "$OUT/q5p" 2>/dev/null)'"
# Q6: max_copies is one haplotype's total at the site (2 + 2 across a spacer), in both modes
qbub q6 h4
for s in 1.0 0.95; do
  qfold q6 "_$s" --min-similarity "$s"
  [ "$(rep "$OUT/q6_$s.panphorte.report.tsv" max_copies)" = 4 ] \
    && ok "max_copies is a haplotype's total across a spacer at --min-similarity $s (4)" \
    || bad "max_copies at $s is $(rep "$OUT/q6_$s.panphorte.report.tsv" max_copies), expected 4"
done

# ---------------------------------------------------------------- a unit made of a short motif
# A 64 bp unit that is eight copies of an 8 bp motif is a short tandem repeat that the 50 bp floor made
# look like a copy unit; counting 64 bp copies of it would read a haplotype one 8 bp copy longer as a whole
# 64 bp copy up. It is reported short_motif and not folded. The same arrays of a primitive 64 bp unit fold.
STR8="ACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCA"
for kind in str8 prim; do
  u=$STR8; [ $kind = prim ] && u=$UNIT
  { printf 'H\tVN:Z:1.0\n'
    printf "S\t1\t%s\nS\t2\t%s\nS\t4\t%s\n" "$BB1" "$u" "$BB2"
    for e in "1 2" "2 2" "2 4" "1 4"; do set -- $e; printf "L\t%s\t+\t%s\t+\t0M\n" "$1" "$2"; done
    printf 'P\tref\t1+,2+,2+,4+\t*\nP\tc3\t1+,2+,2+,2+,4+\t*\nP\tc4\t1+,2+,2+,2+,2+,4+\t*\nP\tz\t1+,4+\t*\n'; } > "$OUT/sm_$kind.gfa"
  "$BIN" bubble -i "$OUT/sm_$kind.gfa" -r ref -o "$OUT/sm_${kind}b" --min-variant-bp 0 --merge-nearby-bp 0 -q >/dev/null 2>&1
  panph -i "$OUT/sm_${kind}b.sorted.gfa" --bubble-prefix "$OUT/sm_${kind}b" -o "$OUT/sm_${kind}p" \
         --min-similarity 0.95 -r ref -q >/dev/null 2>&1
done
[ "$(rep "$OUT/sm_str8p.panphorte.report.tsv" status)" = "short_motif" ] \
  && [ "$(rep "$OUT/sm_str8p.panphorte.report.tsv" normalized)" = "no" ] \
  && ok "a unit of eight 8 bp motifs is reported short_motif and not folded" \
  || bad "a short-motif unit: normalized $(rep "$OUT/sm_str8p.panphorte.report.tsv" normalized), status $(rep "$OUT/sm_str8p.panphorte.report.tsv" status)"
[ "$(rep "$OUT/sm_primp.panphorte.report.tsv" normalized)" = "yes" ] \
  && ok "  ... while the same arrays of a primitive 64 bp unit fold" \
  || bad "the primitive-unit control did not fold (status $(rep "$OUT/sm_primp.panphorte.report.tsv" status))"

# ---------------------------------------------------------------- the re-snarl keeps repeat extension
# A 72 bp (AC)n microsatellite holding a 2 bp deletion and an expansion: bubble extends the expansion's
# site over the repeat (see bubble_stats). Nothing here folds, and the call-ready CSV panphorte writes
# under --reference-path, re-snarled from the graph, must keep that site whole rather than split it
# back to the expansion's own boundaries -- unless --no-fuse-ambiguous says so.
ac() { local o="" i; for ((i = 0; i < $1; i++)); do o="${o}AC"; done; printf '%s' "$o"; }
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\t1\t%s\nS\t9\t%s\n' "$(rnd 300 303)" "$(rnd 300 304)"
  printf 'S\t2\t%s\nS\t3\tAC\nS\t5\t%s\nS\t6\t%s\nS\t7\t%s\nS\t8\t%s\n' "$(ac 10)" "$(ac 10)" "$(ac 5)" "$(ac 40)" "$(ac 10)"
  for e in "1 2" "2 3" "3 5" "2 5" "5 6" "5 7" "6 8" "7 8" "8 9"; do
    set -- $e; printf 'L\t%s\t+\t%s\t+\t0M\n' "$1" "$2"; done
  printf 'P\tref\t1+,2+,3+,5+,6+,8+,9+\t*\nP\thdel\t1+,2+,5+,6+,8+,9+\t*\nP\thexp\t1+,2+,3+,5+,7+,8+,9+\t*\n'; } > "$OUT/tr.gfa"
"$BIN" bubble -i "$OUT/tr.gfa" -r ref -o "$OUT/trb" --merge-nearby-bp 0 -q >/dev/null 2>&1
panph -i "$OUT/trb.sorted.gfa" -c "$OUT/trb.bubbles.csv" -o "$OUT/trp" --reference-path ref -q >/dev/null 2>&1
panph -i "$OUT/trb.sorted.gfa" -c "$OUT/trb.bubbles.csv" -o "$OUT/trpn" --reference-path ref \
       --no-fuse-ambiguous -q >/dev/null 2>&1
csv() { awk -F, 'NR==1{for(i=1;i<=NF;i++)c[$i]=i;next} {printf "%s%s:%s", (n++?",":""), $(c["source"]), $(c["sink"])}' "$1"; }
[ "$(csv "$OUT/trb.bubbles.csv")" = "2:7" ] && [ "$(csv "$OUT/trp.bubbles.csv")" = "2:7" ] \
  && ok "the re-snarled CSV keeps the site bubble extended over the repeat (2:7)" \
  || bad "repeat extension after the re-snarl: bubble $(csv "$OUT/trb.bubbles.csv"), panphorte $(csv "$OUT/trp.bubbles.csv"), expected 2:7"
[ "$(csv "$OUT/trpn.bubbles.csv")" = "4:7" ] \
  && ok "  ... and --no-fuse-ambiguous re-snarls it to its own boundaries (4:7)" \
  || bad "--no-fuse-ambiguous gave $(csv "$OUT/trpn.bubbles.csv"), expected 4:7"

# ---------------------------------------------------------------- a dimer of the unit is not the unit
# A VNTR of a 60 bp unit U whose copies come in two variants, U and V (two bases apart, so V is a copy
# of U at --min-similarity 0.95). u4_* carry U x4, u3 U x3, u1 one U, and uv_* alternate U V U V: an
# exact period of two copies, so they propose the 120 bp unit UV. Folded on UV, every haplotype read
# half its copies and u1 none. UV is only copies of U, so its votes go to U: unit 60, CN 4/3/4/1.
# hor2 gives UV the majority outright (5 uv against 4 U haplotypes); the same holds.
dimer_gfa() {   # dimer_gfa <out> <number of uv haplotypes>
  "${PYTHON:-python3}" - "$1" "$2" <<'PYEOF'
import random, sys
r = random.Random(41)
rnd = lambda n: ''.join(r.choice('ACGT') for _ in range(n))
U = rnd(60)
V = U[:20] + ('A' if U[20] != 'A' else 'C') + U[21:40] + ('G' if U[40] != 'G' else 'T') + U[41:]
P = {'ref': '1+,2+,2+,2+,2+,4+', 'u3': '1+,2+,2+,2+,4+', 'u1': '1+,2+,4+'}
for i in range(1, 5): P[f'u4_{i}'] = '1+,2+,2+,2+,2+,4+'
for i in range(1, int(sys.argv[2]) + 1): P[f'uv_{i}'] = '1+,2+,3+,2+,3+,4+'
with open(sys.argv[1], 'w') as o:
    o.write('H\tVN:Z:1.0\n')
    for k, v in (('1', rnd(2000)), ('2', U), ('3', V), ('4', rnd(2000))): o.write(f'S\t{k}\t{v}\n')
    for a, b in (('1', '2'), ('2', '2'), ('2', '3'), ('3', '2'), ('3', '4'), ('2', '4')): o.write(f'L\t{a}\t+\t{b}\t+\t0M\n')
    for k, v in P.items(): o.write(f'P\t{k}\t{v}\t*\n')
PYEOF
}
for d in "hor 3" "hor2 5"; do
  set -- $d
  dimer_gfa "$OUT/$1.gfa" "$2"
  "$BIN" bubble -i "$OUT/$1.gfa" -r ref -o "$OUT/$1b" --min-variant-bp 0 --merge-nearby-bp 0 -q >/dev/null 2>&1
  panph -i "$OUT/$1b.sorted.gfa" -c "$OUT/$1b.bubbles.csv" -o "$OUT/$1p" --min-similarity 0.95 \
         --reference-path ref -q >/dev/null 2>&1
  "$BIN" call -i "$OUT/$1p.normalized.sorted.gfa" --bubble-prefix "$OUT/$1p" -r ref -o "$OUT/$1c" --cn -q >/dev/null 2>&1
  got="$(rep "$OUT/$1p.panphorte.report.tsv" unit_bp):$(grep -v '^#' "$OUT/$1c.region.vcf" | grep -o 'RU_LEN=[0-9]*;REF_CN=[0-9]*\|REF_CN=[0-9]*;' | head -1):$(fmt "$OUT/$1c.region.vcf" u1 CN):$(grep -vc '^#' "$OUT/$1c.region.vcf")"
  case "$got" in
    60:*REF_CN=4*:1:1) ok "$1: alternating copy variants fold on the 60 bp unit, not its dimer (unit:REF_CN:u1 CN:records = $got)";;
    *) bad "$1: unit:REF_CN:u1 CN:records = $got, expected 60, REF_CN=4, u1 CN 1, one record";;
  esac
done

# What counts as copies of the unit is what folding on it reads, in any phase or orientation, and only
# for a whole multiple. Each case has the 60 bp U array of the test above (ref and u4_* x4, u3 x3, u1 x1)
# and eight haplotypes walking:
#   seed  V W V W, V and W changed at bases 15 and 35/40, inside several of the unit's seed windows
#   tri   V W U V W U, a trimer whose rotation starts on two variant copies
#   rc    U rc(U) U rc(U), copies alternating orientation
#   frac  U H U H U H, H the first half of U: 1.5 units, not a multiple, so it stays its own 90 bp unit
#   del   U D U D, D = U less 2 bp at base 3, before the copy's first seed: it is anchored inside the
#         copy before it, and is aligned from that copy's end
#   ins   U I U I, I = U with 2 bp inserted at base 5: anchored 2 bp late, likewise
#   tridel U E S U E S, E = U less base 5, S = U with base 30 changed
nm_gfa() {   # nm_gfa <out> <case>
  "${PYTHON:-python3}" - "$1" "$2" <<'PYEOF'
import random, sys
r = random.Random(41)
rnd = lambda n: ''.join(r.choice('ACGT') for _ in range(n))
def sub(s, ps):
    s = list(s)
    for p in ps: s[p] = 'A' if s[p] != 'A' else 'C'
    return ''.join(s)
F, G = rnd(2000), rnd(2000)
U = rnd(60)
S = {'1': F, '2': U, '3': sub(U, [15, 35]), '4': sub(U, [15, 40]), '5': U[:30], '6': U[:3] + U[5:],
     '7': U[:5] + 'GT' + U[5:], '8': U[:5] + U[6:], '10': sub(U, [30]), '9': G}
walk = {'seed': '3+,4+,3+,4+', 'tri': '3+,4+,2+,3+,4+,2+', 'rc': '2+,2-,2+,2-', 'frac': '2+,5+,2+,5+,2+,5+',
        'del': '2+,6+,2+,6+', 'ins': '2+,7+,2+,7+', 'tridel': '2+,8+,10+,2+,8+,10+'}[sys.argv[2]]
P = {'ref': '2+,2+,2+,2+', 'u3': '2+,2+,2+', 'u1': '2+'}
for i in range(1, 5): P[f'u4_{i}'] = '2+,2+,2+,2+'
for i in range(1, 9): P[f'uv_{i}'] = walk
L = set()
flip = lambda x: x[:-1] + ('-' if x[-1] == '+' else '+')
for k in P:
    P[k] = '1+,' + P[k] + ',9+'
    st = P[k].split(',')
    for a, b in zip(st, st[1:]): L.add(min((a, b), (flip(b), flip(a))))
with open(sys.argv[1], 'w') as o:
    o.write('H\tVN:Z:1.0\n')
    used = {x[:-1] for v in P.values() for x in v.split(',')}
    for k, v in S.items():
        if k in used: o.write(f'S\t{k}\t{v}\n')
    for a, b in sorted(L): o.write(f'L\t{a[:-1]}\t{a[-1]}\t{b[:-1]}\t{b[-1]}\t0M\n')
    for k, v in P.items(): o.write(f'P\t{k}\t{v}\t*\n')
PYEOF
}
copies_of() { awk -F'\t' -v p="$2" 'NR==1{for(i=1;i<=NF;i++)c[$i]=i;next} $1==p{print $c["copies"]}' "$1"; }
for d in "seed 60 4:4:1" "tri 60 4:6:1" "rc 60 4:4:1" "frac 90 2:3:" "del 60 4:4:1" "ins 60 4:4:1" "tridel 60 4:6:1"; do
  set -- $d
  nm_gfa "$OUT/nm_$1.gfa" "$1"
  "$BIN" bubble -i "$OUT/nm_$1.gfa" -r ref -o "$OUT/nm_$1b" --min-variant-bp 0 --merge-nearby-bp 0 -q >/dev/null 2>&1
  panph -i "$OUT/nm_$1b.sorted.gfa" -c "$OUT/nm_$1b.bubbles.csv" -o "$OUT/nm_$1p" --min-similarity 0.95 \
         --reference-path ref -q >/dev/null 2>&1
  t="$OUT/nm_$1p.panphorte.copies.tsv"
  got="$(rep "$OUT/nm_$1p.panphorte.report.tsv" unit_bp) $(copies_of "$t" ref):$(copies_of "$t" uv_1):$(copies_of "$t" u1)"
  [ "$got" = "$2 $3" ] && ok "$1: unit and copies of ref:uv:u1 are $got" || bad "$1: unit and copies of ref:uv:u1 are $got, expected $2 $3"
done

# A module of five variant sub-copies of U (one base changed in each) and 10 bp of its own, 2/1/0 copies.
# One haplotype also repeats the first sub-copy, so it alone proposes that 60 bp sub-copy as a unit. The
# module is within --min-similarity of five sub-copies, but one haplotype's proposal does not move the
# population's votes: the unit stays the 310 bp module and the copies 2/1.
"${PYTHON:-python3}" - "$OUT/hmod.gfa" <<'PYEOF'
import random, sys
r = random.Random(41)
rnd = lambda n: ''.join(r.choice('ACGT') for _ in range(n))
def sub(s, ps):
    s = list(s)
    for p in ps: s[p] = 'A' if s[p] != 'A' else 'C'
    return ''.join(s)
F, G = rnd(2000), rnd(2000)
U = rnd(60)
S = {'1': F, '2': U, '3': sub(U, [7]), '4': sub(U, [19]), '5': sub(U, [31]), '6': sub(U, [43]), '7': rnd(10), '9': G}
M = '2+,3+,4+,5+,6+,7+'
P = {'ref': M + ',' + M, 'm2a': M + ',' + M, 'm2b': M + ',' + M, 'm2c': M + ',' + M, 'm1a': M, 'm1b': M,
     'out': '2+,' + M, 'm0': ''}
L = set()
for k in P:
    P[k] = '1+,' + (P[k] + ',' if P[k] else '') + '9+'
    st = P[k].split(',')
    for a, b in zip(st, st[1:]): L.add((a, b))
with open(sys.argv[1], 'w') as o:
    o.write('H\tVN:Z:1.0\n')
    for k, v in S.items(): o.write(f'S\t{k}\t{v}\n')
    for a, b in sorted(L): o.write(f'L\t{a[:-1]}\t{a[-1]}\t{b[:-1]}\t{b[-1]}\t0M\n')
    for k, v in P.items(): o.write(f'P\t{k}\t{v}\t*\n')
PYEOF
"$BIN" bubble -i "$OUT/hmod.gfa" -r ref -o "$OUT/hmodb" --min-variant-bp 0 --merge-nearby-bp 0 -q >/dev/null 2>&1
panph -i "$OUT/hmodb.sorted.gfa" -c "$OUT/hmodb.bubbles.csv" -o "$OUT/hmodp" --min-similarity 0.95 \
       --reference-path ref -q >/dev/null 2>&1
t="$OUT/hmodp.panphorte.copies.tsv"
got="$(rep "$OUT/hmodp.panphorte.report.tsv" unit_bp) $(copies_of "$t" ref):$(copies_of "$t" m1a):$(copies_of "$t" out)"
[ "$got" = "310 2:1:1" ] && ok "one haplotype repeating a module's sub-copy does not move the unit (unit copies ref:m1a:out $got)" \
                         || bad "module with an outlier: unit copies ref:m1a:out $got, expected 310 2:1:1"

# A 69 bp unit that is an (AAAG)n repeat running into (AAGG)n through a gained G: past the frameshift
# its bases fall in other phases, so no per-phase consensus reaches 90%, but it fits a run of AAAG at
# 93% with the indel. It is a short tandem repeat, not a copy unit: short_motif, not folded.
W=AAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGGAAAGAAGGAAGGAAGGAAGG
{ printf 'H\tVN:Z:1.0\nS\t1\t%s\nS\t2\t%s\nS\t3\t%s\n' "$(rnd 600 71)" "$W" "$(rnd 600 72)"
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\n'
  printf 'P\tref\t1+,2+,2+,3+\t*\nP\th1\t1+,2+,2+,2+,3+\t*\nP\th2\t1+,2+,2+,2+,2+,3+\t*\nP\th3\t1+,2+,3+\t*\n'; } > "$OUT/fs.gfa"
"$BIN" bubble -i "$OUT/fs.gfa" -r ref -o "$OUT/fsb" --min-variant-bp 0 --merge-nearby-bp 0 -q >/dev/null 2>&1
panph -i "$OUT/fsb.sorted.gfa" -c "$OUT/fsb.bubbles.csv" -o "$OUT/fsp" --min-similarity 0.95 \
       --reference-path ref -q >/dev/null 2>&1
[ "$(rep "$OUT/fsp.panphorte.report.tsv" status):$(rep "$OUT/fsp.panphorte.report.tsv" normalized)" = "short_motif:no" ] \
  && ok "a short tandem repeat with a frameshift is short_motif and not folded" \
  || bad "frameshifted STR: $(rep "$OUT/fsp.panphorte.report.tsv" status):$(rep "$OUT/fsp.panphorte.report.tsv" normalized), expected short_motif:no"

# A 64 bp unit of A and T at 75% A, with no period: about two copies of its own 27 bp half fit at 90%
# with indels, as a skewed alphabet allows; two copies are left to the per-phase test, which it fails.
# It is a copy unit and folds.
W=AATAAAAATAAAAAAATTAAAAAAATAAATAAAAAAAAAAATTATAAAATTAATAAAAAAATAA
{ printf 'H\tVN:Z:1.0\nS\t1\t%s\nS\t2\t%s\nS\t3\t%s\n' "$(rnd 600 73)" "$W" "$(rnd 600 74)"
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\n'
  printf 'P\tref\t1+,2+,2+,3+\t*\nP\th1\t1+,2+,2+,2+,3+\t*\nP\th2\t1+,2+,2+,2+,2+,3+\t*\nP\th3\t1+,2+,3+\t*\n'; } > "$OUT/at.gfa"
"$BIN" bubble -i "$OUT/at.gfa" -r ref -o "$OUT/atb" --min-variant-bp 0 --merge-nearby-bp 0 -q >/dev/null 2>&1
panph -i "$OUT/atb.sorted.gfa" -c "$OUT/atb.bubbles.csv" -o "$OUT/atp" --min-similarity 0.95 \
       --reference-path ref -q >/dev/null 2>&1
[ "$(rep "$OUT/atp.panphorte.report.tsv" status)" = "normalized" ] \
  && ok "an A-rich unit with no period is not taken for a short motif" \
  || bad "A-rich unit: $(rep "$OUT/atp.panphorte.report.tsv" status), expected normalized"

# ---------------------------------------------------------------- letter case does not change copies
# The same graph with the U node in lower case (soft-masked), its variant V in upper case, folds and counts
# as the all-upper-case one: bases are compared without regard to case. Graphs built from soft-masked
# assemblies mix the two.
awk -F'\t' 'BEGIN{OFS="\t"} $1=="S" && $2=="2"{$3=tolower($3)} {print}' "$OUT/hor.gfa" > "$OUT/horlc.gfa"
"$BIN" bubble -i "$OUT/horlc.gfa" -r ref -o "$OUT/horlcb" --min-variant-bp 0 --merge-nearby-bp 0 -q >/dev/null 2>&1
panph -i "$OUT/horlcb.sorted.gfa" -c "$OUT/horlcb.bubbles.csv" -o "$OUT/horlcp" --min-similarity 0.95 \
       --reference-path ref -q >/dev/null 2>&1
[ "$(cut -f2- "$OUT/horlcp.panphorte.copies.tsv")" = "$(cut -f2- "$OUT/horp.panphorte.copies.tsv")" ] \
  && ok "a lower-case (soft-masked) graph gives the same copies as the upper-case one" \
  || bad "letter case changed the copies: $(cut -f2- "$OUT/horlcp.panphorte.copies.tsv" | tr '\n' ' ')"

# ---- fixtures in synthetic_data/panphorte_cases ------------------------------------------------------
CASES="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/synthetic_data/panphorte_cases"
# <report> <bubble id> <column>
pcol() { awk -F'\t' -v id="$2" -v want="$3" 'NR==1{for(i=1;i<=NF;i++)c[$i]=i;next} $1==id{print $(c[want])}' "$1"; }
pcase() {  # <fixture> <ref> <tag> [panphorte args...]: bubble then panphorte
  local f=$1 r=$2 t=$3; shift 3
  "$BIN" bubble -i "$CASES/$f.gfa" -r "$r" -o "$OUT/${t}_b" --min-variant-bp 0 --merge-nearby-bp 0 -q >/dev/null 2>&1
  "$BIN" panphorte -i "$OUT/${t}_b.sorted.gfa" -b "$OUT/${t}_b" -o "$OUT/${t}_p" -q "$@" > "$OUT/$t.log" 2>&1
}

# Contigs that start or end inside a site keep its original route walkable once the haplotypes are
# folded. The site is declined and the run completes.
pcase partial_walks ref pw
[ "$(pcol "$OUT/pw_p.panphorte.report.tsv" 1 status)" = "declined_shared_route" ] \
  && ok "a site whose route partial walks keep is declined, the run completes" \
  || bad "partial walks: status '$(pcol "$OUT/pw_p.panphorte.report.tsv" 1 status)' ($(head -c 100 "$OUT/pw.log"))"

# --refine takes the unit nodes folding created, not "nodes with a self-loop": a unit no haplotype walks
# twice in a row has none, and realigned into a flank node hX would lose a copy. hX's second copy, a
# rotation the reference lacks, is walked forward.
pcase refine_rep_no_loop ref rr -r ref --refine
[ "$(awk -F'\t' '$1=="P" && $2=="hX"{print $3}' "$OUT/rr_p.normalized.sorted.gfa")" = "1+,2+,3+,4+,5+,6+" ] \
  && ok "--refine keeps a unit that has no self-loop (hX still walks two unit copies)" \
  || bad "--refine: hX is $(awk -F'\t' '$1=="P" && $2=="hX"{print $3}' "$OUT/rr_p.normalized.sorted.gfa")"

# A rewritten P line's overlap list has one entry per link, or is `*`: the input's 0M list (five
# entries) stayed on a four-step walk.
pcase p_overlaps arr1 po
[ "$(awk -F'\t' '$1=="P" && $2=="arr1"{print $4}' "$OUT/po_p.normalized.gfa")" = "*" ] \
  && ok "a rewritten P line does not keep its old overlap list" \
  || bad "rewritten P line overlaps: $(awk -F'\t' '$1=="P" && $2=="arr1"{print $3, $4}' "$OUT/po_p.normalized.gfa")"

# Letter case: the period seed compares bases as the aligner does (otherwise a soft-masked second copy in one
# node would read as no seed), and exact mode keys motifs by upper case (otherwise soft-masked and upper-case
# carriers of one repeat would be two motifs, prevalence 0.5 for 3 of 4 carriers).
pcase case_period arr cp --min-similarity 0.95
[ "$(pcol "$OUT/cp_p.panphorte.report.tsv" 1 status)" = "normalized" ] \
  && ok "a soft-masked copy seeds the unit by period" \
  || bad "soft-masked period seed: status '$(pcol "$OUT/cp_p.panphorte.report.tsv" 1 status)'"
pcase case_motif ref cm --min-array-prevalence 0.6
[ "$(pcol "$OUT/cm_p.panphorte.report.tsv" 1 n_motifs)" = "1" ] && [ "$(pcol "$OUT/cm_p.panphorte.report.tsv" 1 status)" = "normalized" ] \
  && ok "soft-masked and upper-case carriers of one repeat are one motif" \
  || bad "case motif: n_motifs $(pcol "$OUT/cm_p.panphorte.report.tsv" 1 n_motifs), status $(pcol "$OUT/cm_p.panphorte.report.tsv" 1 status)"

# A unit under 16 bp cannot seed the approximate copy detector: said so, not "below_prevalence".
pcase short_unit ref su --min-unit-bp 10 --min-similarity 0.9
[ "$(pcol "$OUT/su_p.panphorte.report.tsv" 1 status)" = "unit_under_16bp" ] \
  && ok "an approximate unit under 16 bp is reported as such" \
  || bad "short unit: status '$(pcol "$OUT/su_p.panphorte.report.tsv" 1 status)'"

# A run at a prefix an earlier run used removes that run's files it does not write itself: here the
# sorted graph, site table, Bandage colors and copies of an approximate -r run, left beside an exact
# unsorted one.
"$BIN" panphorte -i "$OUT/cp_b.sorted.gfa" -b "$OUT/cp_b" -o "$OUT/stl" -r arr --min-similarity 0.95 -q >/dev/null 2>&1
"$BIN" panphorte -i "$OUT/cp_b.sorted.gfa" -b "$OUT/cp_b" -o "$OUT/stl" -q >/dev/null 2>&1
[ "$(ls "$OUT"/stl.* | xargs -n1 basename | tr '\n' ' ')" = "stl.normalized.gfa stl.panphorte.report.tsv stl.panphorte.unit_graph.tsv " ] \
  && ok "outputs an earlier run left at the prefix are removed" \
  || bad "files at the prefix after the second run: $(ls "$OUT"/stl.* | xargs -n1 basename | tr '\n' ' ')"

# A haplotype stored reversed folds as the same haplotype stored forward would: every crossing is read
# source -> sink. Read in its own direction, h1rev would cut its copies from the other end, keep the partial
# copy at the other end and get fragment nodes of its own, and the site would read one allele more than it has.
mirror() { awk -F'\t' -v a="$2" -v b="$3" '$1=="P" && $2==a{x=$3} $1=="P" && $2==b{y=$3}
  END{n=split(x,s,","); r=""; for(i=n;i>=1;i--){st=s[i]; o=substr(st,length(st),1); r=r (r==""?"":",") substr(st,1,length(st)-1) (o=="+"?"-":"+")}
      print (r==y) ? "yes" : "no"}' "$1"; }
pcase reversed_copy_in_node ref rci -r ref --min-similarity 0.95
[ "$(mirror "$OUT/rci_p.normalized.sorted.gfa" h1 h1rev)" = "yes" ] \
  && [ "$(cut -d, -f9 "$OUT/rci_p.bubbles.csv" | tail -n +2 | tr -d '\n')" = "3" ] \
  && ok "approximate: a reversed haplotype walks the mirror of its forward twin (3 alleles, not 4)" \
  || bad "approximate reversed twin: mirror $(mirror "$OUT/rci_p.normalized.sorted.gfa" h1 h1rev), alleles $(cut -d, -f9 "$OUT/rci_p.bubbles.csv" | tail -n +2 | tr '\n' ' ')"
pcase reversed_partial ref rpa -r ref
[ "$(mirror "$OUT/rpa_p.normalized.sorted.gfa" h1 h1rev)" = "yes" ] \
  && [ "$(cut -d, -f9 "$OUT/rpa_p.bubbles.csv" | tail -n +2 | tr -d '\n')" = "2" ] \
  && ok "exact: a reversed haplotype keeps its partial copy at the same end as its forward twin (2 alleles)" \
  || bad "exact reversed twin: mirror $(mirror "$OUT/rpa_p.normalized.sorted.gfa" h1 h1rev), alleles $(cut -d, -f9 "$OUT/rpa_p.bubbles.csv" | tail -n +2 | tr '\n' ' ')"

# Exact mode: a node holding two whole copies (UU) is the unit walked twice, not a 128 bp unit walked once.
pcase two_copies_one_node ref tco -r ref
[ "$(pcol "$OUT/tco_p.panphorte.report.tsv" 1 min_copies)" = "2" ] \
  && [ "$(awk -F'\t' '$1=="P" && $2=="h2"{print $3}' "$OUT/tco_p.normalized.sorted.gfa")" = "1+,2+,2+,3+" ] \
  && ok "exact: a node of two whole copies is walked as the unit twice (h2 CN 2)" \
  || bad "two copies in one node: min_copies $(pcol "$OUT/tco_p.panphorte.report.tsv" 1 min_copies), h2 $(awk -F'\t' '$1=="P" && $2=="h2"{print $3}' "$OUT/tco_p.normalized.sorted.gfa")"

# Exact mode applies the short-motif rule too: eight 8 bp motifs are a short tandem repeat, not a 64 bp
# unit (the approximate run of the same fixture is checked above).
panph -i "$OUT/sm_str8b.sorted.gfa" --bubble-prefix "$OUT/sm_str8b" -o "$OUT/sm_str8e" -r ref -q >/dev/null 2>&1
[ "$(rep "$OUT/sm_str8e.panphorte.report.tsv" status)" = "short_motif" ] \
  && ok "exact mode reports a unit of 8 bp motifs short_motif and leaves it" \
  || bad "exact short motif: status $(rep "$OUT/sm_str8e.panphorte.report.tsv" status)"
# Its report row names the unit tried and its carriers, as the approximate run's does.
sm_e="$(rep "$OUT/sm_str8e.panphorte.report.tsv" unit_bp) $(rep "$OUT/sm_str8e.panphorte.report.tsv" n_motif_carriers) $(rep "$OUT/sm_str8e.panphorte.report.tsv" prevalence)"
sm_a="$(rep "$OUT/sm_str8p.panphorte.report.tsv" unit_bp) $(rep "$OUT/sm_str8p.panphorte.report.tsv" n_motif_carriers) $(rep "$OUT/sm_str8p.panphorte.report.tsv" prevalence)"
[ "$sm_e" = "64 3 0.750000" ] && [ "$sm_e" = "$sm_a" ] \
  && ok "exact mode reports the short-motif unit tried (64 bp, 3 carriers, prevalence 0.75)" \
  || bad "exact short-motif row: unit_bp carriers prevalence '$sm_e' (approximate: '$sm_a')"

# ---- the unit graph (default with --min-similarity < 1) ------------------------------------------------
# Copies are detected approximately as above, but folded onto a graph of the unit: segments every copy
# shares, each copy's own differences as nodes of their own. Every path keeps its sequence, and a
# haplotype's copies are counted at the anchor (the nodes covering a unit position no copy deletes, one
# pass per copy):
# the same numbers the single-node fold gives as unit visits above.
anchors() {  # anchor passes per path, in a run without --reference-path (unsorted graph)
  awk -F'\t' 'FNR==1 && FNR==NR{for(i=1;i<=NF;i++)col[$i]=i} FNR==NR{v=$(col["anchor_nodes"]); if(FNR>1 && v!=".") {k=split(v,x,";"); for(j=1;j<=k;j++) an[x[j]]=1}; next} $1=="P"{n=split($3,a,","); c=0
              for(i=1;i<=n;i++) if (substr(a[i],1,length(a[i])-1) in an) c++; printf "%s %d ", $2, c}' \
              "$1.panphorte.unit_graph.tsv" "$1.normalized.gfa"; }
respelled() {  # paths of <out.gfa> that spell differently from <in.gfa>
  python3 -I - "$1" "$2" <<'PY'
import sys, re
def load(f):
    L, P = {}, {}
    for l in open(f):
        x = l.rstrip('\n').split('\t')
        if x[0] == 'S': L[x[1]] = x[2]
        elif x[0] == 'P': P[x[1]] = [(s[:-1], s[-1] == '-') for s in x[2].split(',') if s]
    return L, P
C = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
def spell(L, st): return ''.join(L[n].translate(C)[::-1] if r else L[n] for n, r in st)
Li, Pi = load(sys.argv[1]); Lo, Po = load(sys.argv[2])
print(sum(1 for n in Pi if n not in Po or spell(Li, Pi[n]) != spell(Lo, Po[n])))
PY
}
ug_check() {  # <tag> <input gfa> <prefix> <expected counts>
  [ "$(anchors "$3")" = "$4" ] && [ "$(respelled "$2" "$3.normalized.gfa")" = "0" ] \
    && ok "unit graph, $1: copies counted at the anchor as the single-node fold counts them, every path exact" \
    || bad "unit graph, $1: anchor passes '$(anchors "$3" 2>/dev/null)' (expected '$4'), $(respelled "$2" "$3.normalized.gfa" 2>/dev/null) path(s) respelled"
}
"$BIN" panphorte -i "$OUT/q1.gfa" -c "$OUT/q1.csv" -o "$OUT/ug_q1" --min-similarity 0.95 -q >/dev/null 2>&1
ug_check "double crossing" "$OUT/q1.gfa" "$OUT/ug_q1" "a1 4 a2 6 a3 4 dup 5 "
"$BIN" panphorte -i "$OUT/q2b.sorted.gfa" -c "$OUT/q2b.bubbles.csv" -o "$OUT/ug_q2" --min-similarity 0.95 -q >/dev/null 2>&1
ug_check "a copy with a 3 bp deletion" "$OUT/q2b.sorted.gfa" "$OUT/ug_q2" "h3 3 h2 2 h2b 2 h1d 1 h2d 2 h0 0 "
"$BIN" panphorte -i "$OUT/q3b.sorted.gfa" -c "$OUT/q3b.bubbles.csv" -o "$OUT/ug_q3" --min-similarity 0.95 \
  --min-array-prevalence 0.4 -q >/dev/null 2>&1
ug_check "module with a sub-repeat" "$OUT/q3b.sorted.gfa" "$OUT/ug_q3" "two0 2 two1 2 two2 2 one0 1 one1 1 one2 1 zero 0 "
"$BIN" panphorte -i "$OUT/q5b.sorted.gfa" -c "$OUT/q5b.bubbles.csv" -o "$OUT/ug_q5" --min-similarity 0.95 \
  --max-interruption-frac 0.25 -q >/dev/null 2>&1
ug_check "spacer inside a copy-bearing node" "$OUT/q5b.sorted.gfa" "$OUT/ug_q5" "h3 3 h3b 3 h2 2 h1 1 h0 0 "
# case_period (above): the two copies differ only in letter case, so no segment is shared: the lower-case copy is one
# replacement of the whole unit. The anchor is a unit position, counted on either node covering it.
[ "$(awk -F'\t' 'NR==2{print $6, $7}' "$OUT/cp_p.panphorte.unit_graph.tsv" | tr -d '[0-9]')" = " ;" ] \
  && [ "$(anchors "$OUT/cp_p")" = "arr 2 del 0 " ] \
  && ok "unit graph without a shared segment: copies counted at a position (anchor of 2 nodes)" \
  || bad "case-only copies: anchor '$(awk -F'\t' 'NR==2{print $6, $7}' "$OUT/cp_p.panphorte.unit_graph.tsv")', passes '$(anchors "$OUT/cp_p")'"
# Copies that differ: the reference's third copy has a SNV at 30, h2's second copy one at 40.
#   default            every copy keeps its sequence: cuts at 30 and 40, 5 segments + 2 SNV nodes
#   --lossy-units -r   the reference keeps its sequence (cut at 30 only: 3 segments + 1 SNV node) and h2
#                      walks the plain unit, its SNV dropped
#   --lossy-units      no reference given: every copy walks the unit, one node looped per copy
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\t1\t%s\nS\t2\t%s\nS\t3\t%s\nS\t5\t%s\nS\t4\t%s\n' "$BB1" "$UNIT" "${UNIT:0:30}A${UNIT:31}" \
    "${UNIT:0:40}$([ "${UNIT:40:1}" = C ] && echo G || echo C)${UNIT:41}" "$BB2"
  for e in "1 2" "2 2" "2 3" "3 2" "3 4" "2 4" "1 4" "2 5" "5 2"; do set -- $e; printf 'L\t%s\t+\t%s\t+\t0M\n' "$1" "$2"; done
  printf 'P\tref\t1+,2+,3+,4+\t*\nP\th1\t1+,2+,2+,4+\t*\nP\th2\t1+,2+,5+,2+,4+\t*\nP\th0\t1+,4+\t*\n'; } > "$OUT/ugd.gfa"
"$BIN" bubble -i "$OUT/ugd.gfa" -r ref -o "$OUT/ugdb" --min-variant-bp 0 --merge-nearby-bp 0 -q >/dev/null 2>&1
"$BIN" panphorte -i "$OUT/ugdb.sorted.gfa" -b "$OUT/ugdb" -o "$OUT/ugd_u" --min-similarity 0.95 -q >/dev/null 2>&1
"$BIN" panphorte -i "$OUT/ugdb.sorted.gfa" -b "$OUT/ugdb" -o "$OUT/ugd_r" --min-similarity 0.95 -r ref \
  --lossy-units -q >/dev/null 2>&1
"$BIN" panphorte -i "$OUT/ugdb.sorted.gfa" -b "$OUT/ugdb" -o "$OUT/ugd_l" --min-similarity 0.95 --lossy-units -q >/dev/null 2>&1
ugn() { rep "$1.panphorte.unit_graph.tsv" n_nodes; }
[ "$(respelled "$OUT/ugdb.sorted.gfa" "$OUT/ugd_u.normalized.gfa")" = "0" ] && [ "$(ugn "$OUT/ugd_u")" = "7" ] \
  && [ "$(anchors "$OUT/ugd_u")" = "ref 2 h1 2 h2 3 h0 0 " ] \
  && ok "unit graph: two SNV copies keep their sequence (7 nodes), copies counted 2/2/3/0" \
  || bad "default SNV fixture: respelled $(respelled "$OUT/ugdb.sorted.gfa" "$OUT/ugd_u.normalized.gfa"), nodes $(ugn "$OUT/ugd_u"), anchors '$(anchors "$OUT/ugd_u")'"
ref_kept() { python3 -I - "$1" "$2" <<'PY'
import sys
def load(f):
    L, P = {}, {}
    for l in open(f):
        x = l.rstrip('\n').split('\t')
        if x[0] == 'S': L[x[1]] = x[2]
        elif x[0] == 'P': P[x[1]] = [(s[:-1], s[-1] == '-') for s in x[2].split(',') if s]
    return L, P
C = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
def spell(L, st): return ''.join(L[n].translate(C)[::-1] if r else L[n] for n, r in st)
(Li, Pi), (Lo, Po) = load(sys.argv[1]), load(sys.argv[2])
print(' '.join(n for n in Pi if spell(Li, Pi[n]) != spell(Lo, Po[n])) or '-')
PY
}
[ "$(ref_kept "$OUT/ugdb.sorted.gfa" "$OUT/ugd_r.normalized.sorted.gfa")" = "h2" ] && [ "$(ugn "$OUT/ugd_r")" = "4" ] \
  && ok "--lossy-units -r: the reference keeps its SNV copy (4 nodes), only h2's own SNV is dropped" \
  || bad "--lossy-units -r: respelled '$(ref_kept "$OUT/ugdb.sorted.gfa" "$OUT/ugd_r.normalized.sorted.gfa")', nodes $(ugn "$OUT/ugd_r")"
[ "$(ref_kept "$OUT/ugdb.sorted.gfa" "$OUT/ugd_l.normalized.gfa")" = "ref h2" ] && [ "$(ugn "$OUT/ugd_l")" = "1" ] \
  && [ "$(onenode "$OUT/ugd_l.panphorte.unit_graph.tsv" | wc -l | tr -d ' ')" = "1" ] \
  && ok "--lossy-units without a reference: one node looped per copy, a one-node unit graph" \
  || bad "--lossy-units without -r: respelled '$(ref_kept "$OUT/ugdb.sorted.gfa" "$OUT/ugd_l.normalized.gfa")', nodes $(ugn "$OUT/ugd_l")"
"$BIN" panphorte -i "$OUT/ugdb.sorted.gfa" -b "$OUT/ugdb" -o "$OUT/ugd_x" --min-similarity 1.0 --lossy-units -q >/dev/null 2>&1
[ "$?" -ne 0 ] && ok "--lossy-units is refused in exact mode" || bad "--lossy-units was accepted in exact mode"

# ---- fixtures with their own site tables -------------------------------------------------------------
pcsv() {  # <fixture> <tag> [panphorte args...]: panphorte on the fixture's own site table
  local f=$1 t=$2; shift 2
  "$BIN" panphorte -i "$CASES/$f.gfa" -c "$CASES/$f.csv" -o "$OUT/${t}_p" --threads 1 -q "$@" > "$OUT/$t.log" 2>&1
}
cps() { awk -F'\t' 'NR>1{printf "%s %s ", $1, $4}' "$1.panphorte.copies.tsv"; }
# The seed and its phase read every crossing source -> sink, as detection does: read as each path is stored,
# a unit seeded on the reverse-stored haplotypes would be out of phase for the forward ones, and every
# haplotype would read one copy short.
pcase reverse_majority ref rvm --min-similarity 0.9
[ "$(awk -F'\t' 'NR>1 && $1 !~ /^S6#/{printf "%s %s ", $1, $4}' "$OUT/rvm_p.panphorte.copies.tsv")" = \
  "ref 3 S1#1#ctg1 1 S2#1#ctg2 3 S3#1#ctg3 1 S4#1#ctg4 4 S5#1#ctg5 2 S7#1#ctg7 4 S8#1#ctg8 4 S9#1#ctg9 4 S10#1#ctg10 3 S11#1#ctg11 3 S12#1#ctg12 5 S13#1#ctg13 1 " ] \
  && ok "approximate: haplotypes stored either way get their true copy numbers" \
  || bad "reverse-stored majority: copies '$(cps "$OUT/rvm_p")'"
# A copy missing the unit's first or last 10 bases: the first copy's implied start lies before the interior, and
# its free end must not run into the second copy, or two copies read as one.
pcsv edge_deletions ed --min-similarity 0.95
[ "$(cps "$OUT/ed_p")" = "ref 2 exact 2 del_start 2 del_end 2 del_mid 2 " ] \
  && ok "approximate: copies with a deletion at the unit's start or end are each counted" \
  || bad "terminal deletions: copies '$(cps "$OUT/ed_p")'"
# A node holding two copies (UU), or copies alternating two variants (U V U V), has a node period of two copies:
# the unit is the shortest period, not the dimer, or copy numbers halve and a one-copy carrier reads 0.
pcase unit_two_copies_per_node ref u2n --min-similarity 0.95
pcase unit_alternating ref ual --min-similarity 0.95
[ "$(cps "$OUT/u2n_p")" = "ref 4 h1 5 h2 2 h3 1 h4 6 " ] && [ "$(cps "$OUT/ual_p")" = "ref 4 h1 3 h2 2 h3 1 h4 6 h5 5 " ] \
  && ok "approximate: the unit is the shortest period (two copies per node, alternating variants)" \
  || bad "unit period: UU '$(cps "$OUT/u2n_p")', UV '$(cps "$OUT/ual_p")'"
# Copies drawn from six variants of a 120 bp unit, either strand: as separate motifs each variant would get one
# or two haplotypes' votes, a tie would decide the seed, and the winner could fall below prevalence. Near-identical
# motifs vote as one, and the seed is the one the plain vote picks unless a member of its group detects strictly
# more copies.
pcase variant_votes ref vv --min-similarity 0.9
[ "$(cps "$OUT/vv_p")" = "ref 3 S1#1#ctg1 1 S3#1#ctg3:1-4200 1 S4#1#ctg4:1-4443 3 S5#1#ctg5:1-4444 3 S6#1#ctg6 3 S7#1#ctg7:1-4200 1 S8#1#ctg8 3 S9#1#ctg9 2 S10#1#ctg10 3 S11#1#ctg11:1-4562 4 S12#1#ctg12 3 S13#1#ctg13 4 " ] \
  && ok "approximate: variants of one unit pool their votes, every haplotype at its true copy number" \
  || bad "variant votes: copies '$(cps "$OUT/vv_p")'"
# A contig ending inside the array keeps the first copies' nodes and links; joined to the spacer the rewritten
# haplotypes keep literal, they made a hybrid route (literal copies, then folded ones). Declined in both modes.
pcsv partial_subroute psr1
pcsv partial_subroute psr2 --min-similarity 0.95
[ "$(pcol "$OUT/psr1_p.panphorte.report.tsv" 1 status) $(pcol "$OUT/psr2_p.panphorte.report.tsv" 1 status)" = \
  "declined_shared_route declined_shared_route" ] \
  && ok "a site where an old prefix joins the folded route is declined (exact and approximate)" \
  || bad "hybrid route: status '$(pcol "$OUT/psr1_p.panphorte.report.tsv" 1 status)' / '$(pcol "$OUT/psr2_p.panphorte.report.tsv" 1 status)'"
# The unit-graph table is an output like the others: naming the input CSV, it is refused before anything is written.
cp "$CASES/partial_subroute.csv" "$OUT/al.panphorte.unit_graph.tsv"
"$BIN" panphorte -i "$CASES/partial_subroute.gfa" -c "$OUT/al.panphorte.unit_graph.tsv" -o "$OUT/al" \
  --min-similarity 0.95 -q >/dev/null 2>&1; rc=$?
[ "$rc" -ne 0 ] && cmp -s "$CASES/partial_subroute.csv" "$OUT/al.panphorte.unit_graph.tsv" \
  && ok "an input named like the unit-graph output is refused and left intact" \
  || bad "unit-graph alias: exit $rc, input $(head -c 20 "$OUT/al.panphorte.unit_graph.tsv")"
# --refine: a path crossing the site twice would give one window from its first anchor to its last, and the
# backbone between the crossings would be rebuilt as flank sequence. Such a site is skipped.
pcsv refine_two_crossings rtc -r h1 --refine
[ "$(awk -F'\t' 'NR==2{print $4}' "$OUT/rtc_p.panphorte.refine.tsv")" = "skipped:multiple_crossings" ] \
  && ok "--refine skips a site a path crosses twice" \
  || bad "two crossings: refine status '$(awk -F'\t' 'NR==2{print $4}' "$OUT/rtc_p.panphorte.refine.tsv")'"
# The phase sample holds interiors that carry the unit: with the commonest interiors all non-carriers, the seeded
# rotation stood and every carrier lost a copy. Path names end in their true copy number (_k<N>).
pcsv phase_noncarriers pnc --min-similarity 0.95
[ "$(pcol "$OUT/pnc_p.panphorte.report.tsv" 1 status)" = "normalized" ] && \
  [ "$(awk -F'\t' 'NR>1{split($1,a,"_k"); if (a[2] != $4) bad++} END{print NR-1, bad+0}' "$OUT/pnc_p.panphorte.copies.tsv")" = "12 0" ] \
  && ok "approximate: the phase is set on carriers, every carrier at its true copy number" \
  || bad "phase on non-carriers: $(pcol "$OUT/pnc_p.panphorte.report.tsv" 1 status), copies '$(cps "$OUT/pnc_p")'"
# A repeat outside the winning pooled group does not lead the candidates: five haplotypes with six copies of A
# do not outvote eight carrying four close variants of B.
pcsv pool_unrelated pun --min-similarity 0.95
[ "$(pcol "$OUT/pun_p.panphorte.report.tsv" 1 status) $(pcol "$OUT/pun_p.panphorte.report.tsv" 1 n_motif_carriers)" = "normalized 8" ] \
  && ok "approximate: the pooled vote is not overridden by an unrelated repeat's copy count" \
  || bad "pooled vote: $(pcol "$OUT/pun_p.panphorte.report.tsv" 1 status) $(pcol "$OUT/pun_p.panphorte.report.tsv" 1 n_motif_carriers)"
# Exact mode: n_motif_carriers is the prevalence's numerator, so short-motif carriers are not counted.
pcsv short_motif_carriers smc
[ "$(pcol "$OUT/smc_p.panphorte.report.tsv" 1 n_motif_carriers) $(pcol "$OUT/smc_p.panphorte.report.tsv" 1 prevalence)" = "4 0.333333" ] \
  && ok "exact: carriers and prevalence count the same motifs" \
  || bad "exact carriers: $(pcol "$OUT/smc_p.panphorte.report.tsv" 1 n_motif_carriers) $(pcol "$OUT/smc_p.panphorte.report.tsv" 1 prevalence)"

# ---------------------------------------------------------------- exact units read along the reference
# A unit node is stored the way its copies read along the reference's direction, also where the reference
# carries no copy and the sort does not touch the node. Here the unit's reverse complement sorts first, so
# a node stored by the lesser strand would have every copy walk it reversed, read by call as inverted.
UG="G$(rnd 62 4242)T"
{ printf 'H\tVN:Z:1.0\n'
  printf "S\t1\t%s\nS\t2\t%s\nS\t3\t%s\n" "$BB1" "$UG" "$BB2"
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t1\t+\t3\t+\t0M\n'
  printf 'P\tref\t1+,3+\t*\n'
  printf 'P\th0\t1+,2+,2+,3+\t*\nP\th1\t1+,2+,2+,2+,3+\t*\nP\th2\t1+,2+,2+,3+\t*\nP\th3\t1+,2+,2+,2+,3+\t*\n'; } > "$OUT/noref.gfa"
"$BIN" bubble -i "$OUT/noref.gfa" -r ref -o "$OUT/norefb" --min-variant-bp 0 --merge-nearby-bp 0 -q >/dev/null 2>&1
panph -i "$OUT/norefb.sorted.gfa" -b "$OUT/norefb" -o "$OUT/norefp" -r ref -q >/dev/null 2>&1
nr_unit=$(onenode "$OUT/norefp.panphorte.unit_graph.tsv" | head -1)
nr_walks=$(awk -F'\t' -v u="$nr_unit" '$1=="P" && $2 ~ /^h/ {n=split($3,a,","); f=0; r=0
             for(i=1;i<=n;i++){ if (a[i]==u"+") f++; if (a[i]==u"-") r++ } printf "%s:%d/%d ", $2, f, r}' \
             "$OUT/norefp.normalized.sorted.gfa")
[ -n "$nr_unit" ] && [ "$nr_walks" = "h0:2/0 h1:3/0 h2:2/0 h3:3/0 " ] \
  && ok "exact: with no reference copy the carriers walk the unit forward (forward/reverse passes)" \
  || bad "exact, no reference copy: unit '$nr_unit', passes forward/reverse '$nr_walks'"

# ---------------------------------------------------------------- a soft-masked copy inside a node
# A node holding two copies, the second soft-masked, is two copies: the period ignores letter case, as the
# motif does. Each copy keeps its own bases, so the fold stays exact.
{ printf 'H\tVN:Z:1.0\n'
  printf "S\t1\t%s\nS\t2\t%s%s\nS\t3\t%s\n" "$BB1" "$UNIT" "$(echo "$UNIT" | tr 'ACGT' 'acgt')" "$BB2"
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t1\t+\t3\t+\t0M\n'
  printf 'P\tref\t1+,2+,3+\t*\nP\th0\t1+,2+,2+,3+\t*\nP\th1\t1+,2+,2+,2+,3+\t*\nP\th2\t1+,2+,2+,3+\t*\nP\th3\t1+,3+\t*\n'; } > "$OUT/mixcase.gfa"
"$BIN" bubble -i "$OUT/mixcase.gfa" -r ref -o "$OUT/mixcaseb" --min-variant-bp 0 --merge-nearby-bp 0 -q >/dev/null 2>&1
panph -i "$OUT/mixcaseb.sorted.gfa" -b "$OUT/mixcaseb" -o "$OUT/mixcasep" -r ref -q >/dev/null 2>&1
# passes per path through the anchors of a unit graph table: <prefix>
anchor_passes() { awk -F'\t' 'NR==FNR {if (FNR==1) {for(i=1;i<=NF;i++) c[$i]=i; next}
                     n=split($(c["anchor_nodes"]),a,";"); for(i=1;i<=n;i++) anc[a[i]]=1; next}
                   $1=="P" {n=split($3,s,","); k=0; for(i=1;i<=n;i++) if (substr(s[i],1,length(s[i])-1) in anc) k++
                            printf "%s:%d ", $2, k}' "$1.panphorte.unit_graph.tsv" "$1.normalized.sorted.gfa"; }
mc="$(rep "$OUT/mixcasep.panphorte.report.tsv" unit_bp) $(rep "$OUT/mixcasep.panphorte.report.tsv" max_copies)"
[ "$mc" = "64 6" ] && [ "$(anchor_passes "$OUT/mixcasep")" = "ref:2 h0:4 h1:6 h2:4 h3:0 " ] \
  && ok "exact: a node of a copy and a soft-masked copy counts two copies (unit 64 bp)" \
  || bad "mixed-case node: unit_bp max_copies '$mc', passes '$(anchor_passes "$OUT/mixcasep")'"

# ---------------------------------------------------------------- --refine with soft-masked flanks
# The refine fixture with distinct left flanks (B's ends GTAA) and the flank nodes in lower case: the
# alignment keeps each flank's own bases, so the rebuild reproduces them and the site is refined.
awk -F'\t' 'BEGIN{OFS="\t"} $1=="S" && $2=="9" {$3="GTAA"} $1=="S" && $2>=5 && $2<=11 {$3=tolower($3)} {print}' \
  "$OUT/rf1.gfa" > "$OUT/rflc.gfa"
"$BIN" bubble -i "$OUT/rflc.gfa" -r A -o "$OUT/rflcb" --min-variant-bp 0 --merge-nearby-bp 0 -q >/dev/null 2>&1
rf_run rflc rflcn
rf_run rflc rflcr --refine
rflc_ok=1
for h in A A2 B C R D; do
  [ "$(spell "$OUT/rflcn.normalized.sorted.gfa" $h)" = "$(spell "$OUT/rflcr.normalized.sorted.gfa" $h)" ] || rflc_ok=0
done
[ "$(rcol "$OUT/rflcr.panphorte.refine.tsv" refine_status)" = "refined" ] && [ "$rflc_ok" = 1 ] \
  && ok "--refine realigns soft-masked flanks and every haplotype keeps its bases" \
  || bad "--refine, soft-masked flanks: status '$(rcol "$OUT/rflcr.panphorte.refine.tsv" refine_status)', spelled kept $rflc_ok"

# ---------------------------------------------------------------- copies every haplotype shares beside the site
# A graph that lays every copy position out as its own nodes leaves the first copy, which every haplotype
# carries, outside the site as shared backbone. The site is widened over it, so each haplotype's copies are
# all counted: 1, 3, 6 and 3 copies of a 60 bp unit cut 24/30/6, in exact and approximate mode.
U60="$(rnd 60 6060)"
P1=${U60:0:24}; P2=${U60:24:30}; P3=${U60:54:6}
{ printf 'H\tVN:Z:1.0\n'
  printf "S\t1\t%s\nS\t2\t%s\nS\t3\t%s\nS\t4\t%s\nS\t99\t%s\n" "$BB1" "$P1" "$P2" "$P3" "$BB2"
  for c in 1 2 3 4 5; do b=$((10 * c))
    printf "S\t%d\t%s\nS\t%d\t%s\nS\t%d\t%s\n" $b "$P1" $((b + 1)) "$P2" $((b + 2)) "$P3"
    printf "L\t%d\t+\t%d\t+\t0M\nL\t%d\t+\t%d\t+\t0M\nL\t%d\t+\t99\t+\t0M\n" $b $((b + 1)) $((b + 1)) $((b + 2)) $((b + 2))
    [ $c -lt 5 ] && printf "L\t%d\t+\t%d\t+\t0M\n" $((b + 2)) $((b + 10))
  done
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t3\t+\t4\t+\t0M\nL\t4\t+\t99\t+\t0M\nL\t4\t+\t10\t+\t0M\n'
  walk() { local s="1+,2+,3+,4+" c; for ((c = 1; c < $1; c++)); do b=$((10 * c)); s="$s,$b+,$((b + 1))+,$((b + 2))+"; done; echo "$s,99+"; }
  printf 'P\tref\t%s\t*\nP\ts1\t%s\t*\nP\ts2\t%s\t*\nP\ts3\t%s\t*\n' "$(walk 1)" "$(walk 3)" "$(walk 6)" "$(walk 3)"; } > "$OUT/unrolled.gfa"
"$BIN" bubble -i "$OUT/unrolled.gfa" -r ref -o "$OUT/unrolledb" --min-variant-bp 0 --merge-nearby-bp 0 -q >/dev/null 2>&1
panph -i "$OUT/unrolledb.sorted.gfa" -b "$OUT/unrolledb" -o "$OUT/unrolledx" -r ref -q >/dev/null 2>&1
panph -i "$OUT/unrolledb.sorted.gfa" -b "$OUT/unrolledb" -o "$OUT/unrolleda" -r ref --min-similarity 0.9 -q >/dev/null 2>&1
"$BIN" panphorte -i "$OUT/unrolledb.sorted.gfa" -b "$OUT/unrolledb" -o "$OUT/unrolledg" -r ref --min-similarity 0.9 -q >/dev/null 2>&1
[ "$(anchor_passes "$OUT/unrolledx")" = "ref:1 s1:3 s2:6 s3:3 " ] \
  && ok "exact: a copy every haplotype shares beside the site is counted (1/3/6/3)" \
  || bad "exact, shared boundary copy: passes '$(anchor_passes "$OUT/unrolledx")', expected 1/3/6/3"
[ "$(anchor_passes "$OUT/unrolleda")" = "ref:1 s1:3 s2:6 s3:3 " ] && [ "$(anchor_passes "$OUT/unrolledg")" = "ref:1 s1:3 s2:6 s3:3 " ] \
  && ok "approximate: the shared boundary copy is counted too (lossy and unit graph)" \
  || bad "approximate, shared boundary copy: passes '$(anchor_passes "$OUT/unrolleda")' / '$(anchor_passes "$OUT/unrolledg")'"
ur_ok=1
for h in ref s1 s2 s3; do [ "$(spell "$OUT/unrolledb.sorted.gfa" $h)" = "$(spell "$OUT/unrolledx.normalized.sorted.gfa" $h)" ] || ur_ok=0; done
[ "$ur_ok" = 1 ] && ok "and the widened fold keeps every haplotype's sequence" || bad "the widened fold changed a haplotype's sequence"

# ---------------------------------------------------------------- a lone copy cut into other nodes
# One site whose 60 bp unit is node 3 in the array and nodes 5,6,7 (24/30/6) for the copy after it. The
# array 3,3 is found first; the copy cut into 5,6,7 still folds as its own copy, on either strand, so
# every haplotype counts 2/3/4/2/3 copies.
U60="$(rnd 60 6161)"; P1=${U60:0:24}; P2=${U60:24:30}; P3=${U60:54:6}
spl_walk() { local s="1+,2+" c; for ((c = 0; c < $1; c++)); do s="$s,3+"; done; [ "$2" = 1 ] && s="$s,5+,6+,7+"; echo "$s,8+,9+"; }
{ printf 'H\tVN:Z:1.0\n'
  printf "S\t1\t%s\nS\t2\t%s\nS\t3\t%s\nS\t5\t%s\nS\t6\t%s\nS\t7\t%s\nS\t8\t%s\nS\t9\t%s\n" \
    "$BB1" "$(rnd 100 6262)" "$U60" "$P1" "$P2" "$P3" "$(rnd 100 6363)" "$BB2"
  for e in "1 2" "2 3" "3 3" "3 5" "2 5" "5 6" "6 7" "7 8" "2 8" "8 9"; do set -- $e; printf 'L\t%s\t+\t%s\t+\t0M\n' "$1" "$2"; done
  printf 'P\tref\t%s\t*\nP\th1\t%s\t*\nP\th2\t%s\t*\nP\th3\t%s\t*\nP\th4\t%s\t*\nP\th5\t%s\t*\n' \
    "$(spl_walk 1 1)" "$(spl_walk 2 1)" "$(spl_walk 3 1)" "$(spl_walk 1 1)" "$(spl_walk 2 1)" "$(spl_walk 0 0)"; } > "$OUT/spl.gfa"
SITE_HDR='bubble_id,source,source_orient,sink,sink_orient,path_support,min_inside_bp,max_inside_bp,inside_nodes'
printf '%s\n1,2,+,8,+,6,0,240,"3;5;6;7"\n' "$SITE_HDR" > "$OUT/spl.csv"
# h1 and h2 stored on the other strand
awk -F'\t' 'BEGIN{OFS="\t"} $1=="P" && ($2=="h1" || $2=="h2") {n=split($3,a,","); r=""
  for(i=n;i>=1;i--){x=a[i]; o=substr(x,length(x)); r=r (r==""?"":",") substr(x,1,length(x)-1) (o=="+"?"-":"+")}; $3=r} {print}' \
  "$OUT/spl.gfa" > "$OUT/splrev.gfa"
panph -i "$OUT/spl.gfa" -c "$OUT/spl.csv" -o "$OUT/splx" -r ref -q >/dev/null 2>&1
panph -i "$OUT/splrev.gfa" -c "$OUT/spl.csv" -o "$OUT/splrevx" -r ref -q >/dev/null 2>&1
spl_ok=1
for h in ref h1 h2 h3 h4 h5; do
  [ "$(spell "$OUT/spl.gfa" $h)" = "$(spell "$OUT/splx.normalized.sorted.gfa" $h)" ] || spl_ok=0
  [ "$(spell "$OUT/splrev.gfa" $h)" = "$(spell "$OUT/splrevx.normalized.sorted.gfa" $h)" ] || spl_ok=0
done
[ "$(anchor_passes "$OUT/splx")" = "ref:2 h1:3 h2:4 h3:2 h4:3 h5:0 " ] \
  && [ "$(anchor_passes "$OUT/splrevx")" = "ref:2 h1:3 h2:4 h3:2 h4:3 h5:0 " ] && [ "$spl_ok" = 1 ] \
  && ok "exact: a copy cut into other nodes beside an array folds too (2/3/4/2/3, either strand)" \
  || bad "split lone copy: passes '$(anchor_passes "$OUT/splx")' / reversed '$(anchor_passes "$OUT/splrevx")', spelled kept $spl_ok"

# ---------------------------------------------------------------- where a site is not widened
# Unit passes of one site: steps on the anchor nodes of its unit graph rows.
site_passes() { awk -F'\t' -v id="$2" 'NR==FNR {if (FNR==1) {for(i=1;i<=NF;i++) c[$i]=i; next}
                     if ($1==id) {n=split($(c["anchor_nodes"]),a,";"); for(i=1;i<=n;i++) anc[a[i]]=1}; next}
                   $1=="P" {n=split($3,s,","); k=0; for(i=1;i<=n;i++) if (substr(s[i],1,length(s[i])-1) in anc) k++
                            printf "%s:%d ", $2, k}' "$1.panphorte.unit_graph.tsv" "$1.normalized.sorted.gfa"; }
# Two sites, A = 2..5 and B = 8..12, with one copy of a 60 bp unit (5,6,7 cut 24/30/6) every haplotype walks
# between them. A reads it from its sink as P1 P2 P3, B from its source as P2 P3 P1: both could take 6 and 7.
# Sites are widened in id order and a node one site took is no other's, so A takes the copy and B keeps
# its own: each copy is counted once, at A.
U60="$(rnd 60 7070)"; P1=${U60:0:24}; P2=${U60:24:30}; P3=${U60:54:6}
nb_walk() { local s="1+,2+" c; for ((c = 0; c < $1; c++)); do s="$s,3+,4+,20+"; done; s="$s,5+,6+,7+,8+"
            for ((c = 0; c < $2; c++)); do s="$s,9+,10+,11+"; done; echo "$s,12+,13+"; }
{ printf 'H\tVN:Z:1.0\n'
  printf "S\t1\t%s\nS\t2\t%s\nS\t12\t%s\nS\t13\t%s\n" "$BB1" "$(rnd 100 7171)" "$(rnd 100 7272)" "$BB2"
  printf "S\t3\t%s\nS\t4\t%s\nS\t20\t%s\n" "$P1" "$P2" "$P3"
  printf "S\t5\t%s\nS\t6\t%s\nS\t7\t%s\nS\t8\t%s\nS\t9\t%s\nS\t10\t%s\nS\t11\t%s\n" "$P1" "$P2" "$P3" "$P1" "$P2" "$P3" "$P1"
  for e in "1 2" "2 3" "3 4" "4 20" "20 3" "20 5" "2 5" "5 6" "6 7" "7 8" "8 9" "9 10" "10 11" "11 9" "11 12" "8 12" "12 13"; do
    set -- $e; printf 'L\t%s\t+\t%s\t+\t0M\n' "$1" "$2"; done
  printf 'P\tref\t%s\t*\nP\th1\t%s\t*\nP\th2\t%s\t*\nP\th3\t%s\t*\nP\th4\t%s\t*\nP\th5\t%s\t*\n' \
    "$(nb_walk 1 1)" "$(nb_walk 2 2)" "$(nb_walk 0 3)" "$(nb_walk 3 0)" "$(nb_walk 2 1)" "$(nb_walk 1 2)"; } > "$OUT/nb.gfa"
printf '%s\n1,2,+,5,+,6,0,180,"3;4;20"\n2,8,+,12,+,6,0,180,"9;10;11"\n' "$SITE_HDR" > "$OUT/nb.csv"
panph -i "$OUT/nb.gfa" -c "$OUT/nb.csv" -o "$OUT/nbx" -r ref -q >/dev/null 2>&1; nb_rc=$?
"$BIN" panphorte -i "$OUT/nb.gfa" -c "$OUT/nb.csv" -o "$OUT/nba" -r ref --min-similarity 0.95 -q >/dev/null 2>&1; nba_rc=$?
nb_ok=1
for h in ref h1 h2 h3 h4 h5; do
  [ "$(spell "$OUT/nb.gfa" $h)" = "$(spell "$OUT/nbx.normalized.sorted.gfa" $h)" ] || nb_ok=0
  [ "$(spell "$OUT/nb.gfa" $h)" = "$(spell "$OUT/nba.normalized.sorted.gfa" $h)" ] || nb_ok=0
done
nb_want="ref:2 h1:3 h2:1 h3:4 h4:3 h5:2 | ref:1 h1:2 h2:3 h3:0 h4:1 h5:2 "
[ "$nb_rc:$nba_rc" = "0:0" ] && [ "$(site_passes "$OUT/nbx" 1)| $(site_passes "$OUT/nbx" 2)" = "$nb_want" ] \
  && [ "$(site_passes "$OUT/nba" 1)| $(site_passes "$OUT/nba" 2)" = "$nb_want" ] && [ "$nb_ok" = 1 ] \
  && ok "two neighbouring sites never widen over the same copy: it is counted once, at the first" \
  || bad "neighbouring widenings: exit $nb_rc/$nba_rc, exact '$(site_passes "$OUT/nbx" 1)| $(site_passes "$OUT/nbx" 2)', approximate '$(site_passes "$OUT/nba" 1)| $(site_passes "$OUT/nba" 2)', spelled kept $nb_ok"
# Runs <name>.gfa with <name>.csv in exact (x) and approximate (a) mode; prints each mode's exit code, the
# passes of site <id>, the status of site <other> and whether every haplotype is spelled as before.
nb_both() {
  local name=$1 id=$2 other=$3 m h kept
  for m in x a; do
    if [ $m = x ]; then panph -i "$OUT/$name.gfa" -c "$OUT/$name.csv" -o "$OUT/$name$m" -r ref -q >/dev/null 2>&1
    else "$BIN" panphorte -i "$OUT/$name.gfa" -c "$OUT/$name.csv" -o "$OUT/$name$m" -r ref --min-similarity 0.95 \
           -q >/dev/null 2>&1; fi
    printf '%s %s| %s ' "$?" "$(site_passes "$OUT/$name$m" "$id")" \
      "$(awk -F'\t' -v id="$other" '$1==id {print $NF}' "$OUT/$name$m.panphorte.report.tsv" 2>/dev/null)"
    kept=1
    for h in ref h1 h2 h3 h4 h5; do
      [ "$(spell "$OUT/$name.gfa" $h)" = "$(spell "$OUT/$name$m.normalized.sorted.gfa" $h)" ] || kept=0; done
    printf '%s; ' "$kept"
  done; }
# The nb layout with other walks: <name> then "a,m" (A's and B's own copies) for ref h1 .. h5.
nb_variant() {
  local name=$1 h; shift
  { grep -v '^P' "$OUT/nb.gfa"
    for h in ref h1 h2 h3 h4 h5; do printf 'P\t%s\t%s\t*\n' $h "$(nb_walk "${1%,*}" "${1#*,}")"; shift; done; } > "$OUT/$name.gfa"; }
# A site that does not fold edits nothing and takes nothing: A is widened over the shared copy, but too few
# haplotypes carry its array, so B takes the copy and reads 1 + its own copies. Also with the ids swapped,
# B first and left unfolded, A taking the copy.
nb_variant nbu 0,1 0,2 0,3 2,0 0,1 0,2
cp "$OUT/nb.csv" "$OUT/nbu.csv"
nb_variant nbv 1,0 2,0 3,0 0,2 1,0 2,0
printf '%s\n1,8,+,12,+,6,0,180,"9;10;11"\n2,2,+,5,+,6,0,180,"3;4;20"\n' "$SITE_HDR" > "$OUT/nbv.csv"
nbu_want="0 ref:2 h1:3 h2:4 h3:1 h4:2 h5:3 | below_prevalence 1; "
[ "$(nb_both nbu 2 1)" = "$nbu_want$nbu_want" ] && [ "$(nb_both nbv 2 1)" = "$nbu_want$nbu_want" ] \
  && ok "a widened site left unfolded leaves the shared copy to its neighbour (2/3/4/1/2/3, either id order)" \
  || bad "unfolded neighbour: '$(nb_both nbu 2 1)' / ids swapped '$(nb_both nbv 2 1)', expected '$nbu_want' per mode"
# Two sites with a 60 bp copy on single nodes on each facing side: A = 2..5 takes node 5 and B = 6..8 takes
# node 6, each with the other's old boundary as its new one. The edits lie side by side, so both are widened
# and every copy is counted once: A reads a + 1, B m + 1. Also with the ids swapped.
U60="$(rnd 60 7373)"
sd_walk() { local s="1+,2+" c; for ((c = 0; c < $1; c++)); do s="$s,3+"; done; s="$s,5+,6+"
            for ((c = 0; c < $2; c++)); do s="$s,7+"; done; echo "$s,8+,9+"; }
{ printf 'H\tVN:Z:1.0\n'
  printf "S\t1\t%s\nS\t2\t%s\nS\t8\t%s\nS\t9\t%s\n" "$BB1" "$(rnd 100 7474)" "$(rnd 100 7575)" "$BB2"
  printf "S\t3\t%s\nS\t5\t%s\nS\t6\t%s\nS\t7\t%s\n" "$U60" "$U60" "$U60" "$U60"
  for e in "1 2" "2 3" "3 3" "3 5" "2 5" "5 6" "6 7" "7 7" "7 8" "6 8" "8 9"; do
    set -- $e; printf 'L\t%s\t+\t%s\t+\t0M\n' "$1" "$2"; done
  printf 'P\tref\t%s\t*\nP\th1\t%s\t*\nP\th2\t%s\t*\nP\th3\t%s\t*\nP\th4\t%s\t*\nP\th5\t%s\t*\n' \
    "$(sd_walk 1 2)" "$(sd_walk 2 3)" "$(sd_walk 1 0)" "$(sd_walk 2 4)" "$(sd_walk 3 2)" "$(sd_walk 0 1)"; } > "$OUT/sd.gfa"
cp "$OUT/sd.gfa" "$OUT/sdv.gfa"
printf '%s\n1,2,+,5,+,6,0,180,"3"\n2,6,+,8,+,6,0,240,"7"\n' "$SITE_HDR" > "$OUT/sd.csv"
printf '%s\n1,6,+,8,+,6,0,240,"7"\n2,2,+,5,+,6,0,180,"3"\n' "$SITE_HDR" > "$OUT/sdv.csv"
sd_a="ref:2 h1:3 h2:2 h3:3 h4:4 h5:1 "; sd_b="ref:3 h1:4 h2:1 h3:5 h4:3 h5:2 "
sd_want() { printf '0 %s| normalized 1; 0 %s| normalized 1; ' "$1" "$1"; }
[ "$(nb_both sd 1 2)" = "$(sd_want "$sd_a")" ] && [ "$(nb_both sd 2 1)" = "$(sd_want "$sd_b")" ] \
  && [ "$(nb_both sdv 1 2)" = "$(sd_want "$sd_b")" ] && [ "$(nb_both sdv 2 1)" = "$(sd_want "$sd_a")" ] \
  && ok "neighbours each widened over their own copy on the facing side both fold, each copy counted once" \
  || bad "side-by-side widenings: A '$(nb_both sd 1 2)' B '$(nb_both sd 2 1)', ids swapped B '$(nb_both sdv 1 2)' A '$(nb_both sdv 2 1)'; expected A $sd_a, B $sd_b"
# Every haplotype with copies on both sides, so the copies every crossing shares run two deep past each
# site: each side stops before the other site's node, and A still reads a + 1, B m + 1.
{ grep -v '^P' "$OUT/sd.gfa"
  printf 'P\tref\t%s\t*\nP\th1\t%s\t*\nP\th2\t%s\t*\nP\th3\t%s\t*\nP\th4\t%s\t*\nP\th5\t%s\t*\n' \
    "$(sd_walk 1 1)" "$(sd_walk 2 2)" "$(sd_walk 1 3)" "$(sd_walk 3 1)" "$(sd_walk 2 1)" "$(sd_walk 1 2)"; } > "$OUT/sd2.gfa"
cp "$OUT/sd2.gfa" "$OUT/sd2v.gfa"; cp "$OUT/sd.csv" "$OUT/sd2.csv"; cp "$OUT/sdv.csv" "$OUT/sd2v.csv"
sd_a="ref:2 h1:3 h2:2 h3:4 h4:3 h5:2 "; sd_b="ref:2 h1:3 h2:4 h3:2 h4:2 h5:3 "
[ "$(nb_both sd2 1 2)" = "$(sd_want "$sd_a")" ] && [ "$(nb_both sd2 2 1)" = "$(sd_want "$sd_b")" ] \
  && [ "$(nb_both sd2v 1 2)" = "$(sd_want "$sd_b")" ] && [ "$(nb_both sd2v 2 1)" = "$(sd_want "$sd_a")" ] \
  && ok "neighbours whose shared copies run two deep each take the copy on their own side" \
  || bad "two-deep side-by-side widenings: A '$(nb_both sd2 1 2)' B '$(nb_both sd2 2 1)', ids swapped B '$(nb_both sd2v 1 2)' A '$(nb_both sd2v 2 1)'; expected A $sd_a, B $sd_b"
# The shared copy's nodes 6 and 7 are another site's interior: A is not widened over them, in either mode.
printf '%s\n1,2,+,5,+,6,0,180,"3;4;20"\n2,5,+,8,+,6,36,36,"6;7"\n' "$SITE_HDR" > "$OUT/nbo.csv"
panph -i "$OUT/nb.gfa" -c "$OUT/nbo.csv" -o "$OUT/nbox" -r ref --bubble-id 1 -q >/dev/null 2>&1
"$BIN" panphorte -i "$OUT/nb.gfa" -c "$OUT/nbo.csv" -o "$OUT/nboa" -r ref --bubble-id 1 --min-similarity 0.95 -q >/dev/null 2>&1
[ "$(site_passes "$OUT/nbox" 1)" = "ref:1 h1:2 h2:0 h3:3 h4:2 h5:1 " ] \
  && [ "$(site_passes "$OUT/nboa" 1)" = "ref:1 h1:2 h2:0 h3:3 h4:2 h5:1 " ] \
  && ok "a site is not widened over another site's nodes" \
  || bad "widened over another site's nodes: passes '$(site_passes "$OUT/nbox" 1)' / '$(site_passes "$OUT/nboa" 1)', expected 1/2/0/3/2/1"
# A table whose other site holds node 8 inside: A's new sink would lie inside it, so A is not widened.
printf '%s\n1,2,+,5,+,6,0,180,"3;4;20"\n2,12,+,13,+,6,0,0,"8"\n' "$SITE_HDR" > "$OUT/nbi.csv"
"$BIN" panphorte -i "$OUT/nb.gfa" -c "$OUT/nbi.csv" -o "$OUT/nbia" -r ref --bubble-id 1 --min-similarity 0.95 -q >/dev/null 2>&1
[ "$(site_passes "$OUT/nbia" 1)" = "ref:1 h1:2 h2:0 h3:3 h4:2 h5:1 " ] \
  && ok "a site is not widened to a boundary inside another site" \
  || bad "widened to a boundary inside another site: passes '$(site_passes "$OUT/nbia" 1)', expected 1/2/0/3/2/1"
# A contig walks the shared copy without crossing A: widened, its walk would keep the copy's old route.
{ cat "$OUT/nb.gfa"; printf 'P\tfrag\t6+,7+,8+,9+,10+,11+,12+,13+\t*\n'; } > "$OUT/nbf.gfa"
printf '%s\n1,2,+,5,+,6,0,180,"3;4;20"\n' "$SITE_HDR" > "$OUT/nbf.csv"
panph -i "$OUT/nbf.gfa" -c "$OUT/nbf.csv" -o "$OUT/nbfx" -r ref -q >/dev/null 2>&1
"$BIN" panphorte -i "$OUT/nbf.gfa" -c "$OUT/nbf.csv" -o "$OUT/nbfa" -r ref --min-similarity 0.95 -q >/dev/null 2>&1
[ "$(site_passes "$OUT/nbfx" 1)" = "ref:1 h1:2 h2:0 h3:3 h4:2 h5:1 frag:0 " ] \
  && [ "$(site_passes "$OUT/nbfa" 1)" = "ref:1 h1:2 h2:0 h3:3 h4:2 h5:1 frag:0 " ] \
  && ok "a site is not widened over a copy a walk passes outside its crossings" \
  || bad "widened over a copy walked elsewhere: passes '$(site_passes "$OUT/nbfx" 1)' / '$(site_passes "$OUT/nbfa" 1)'"
# The shared copies are a whole-unit node and a copy cut 24/30/6, unlike the interior's: the site is
# widened over both and both fold in exact mode, so every haplotype reads its true count.
U60="$(rnd 60 8080)"; P1=${U60:0:24}; P2=${U60:24:30}; P3=${U60:54:6}
cut_walk() { local s="1+,2+" c; for ((c = 0; c < $1; c++)); do s="$s,3+"; done; echo "$s,4+,5+,6+,7+,8+,12+,13+"; }
{ printf 'H\tVN:Z:1.0\n'
  printf "S\t1\t%s\nS\t2\t%s\nS\t3\t%s\nS\t4\t%s\nS\t5\t%s\nS\t6\t%s\nS\t7\t%s\nS\t8\t%s\nS\t12\t%s\nS\t13\t%s\n" \
    "$BB1" "$(rnd 100 8181)" "$U60" "$U60" "$P1" "$P2" "$P3" "$(rnd 100 8282)" "$(rnd 100 8383)" "$BB2"
  for e in "1 2" "2 3" "3 3" "3 4" "2 4" "4 5" "5 6" "6 7" "7 8" "8 12" "12 13"; do set -- $e; printf 'L\t%s\t+\t%s\t+\t0M\n' "$1" "$2"; done
  printf 'P\tref\t%s\t*\nP\th1\t%s\t*\nP\th2\t%s\t*\nP\th3\t%s\t*\n' "$(cut_walk 1)" "$(cut_walk 2)" "$(cut_walk 3)" "$(cut_walk 0)"; } > "$OUT/cut.gfa"
printf '%s\n1,2,+,4,+,4,0,180,"3"\n' "$SITE_HDR" > "$OUT/cut.csv"
panph -i "$OUT/cut.gfa" -c "$OUT/cut.csv" -o "$OUT/cutx" -r ref -q >/dev/null 2>&1
cut_ok=1
for h in ref h1 h2 h3; do [ "$(spell "$OUT/cut.gfa" $h)" = "$(spell "$OUT/cutx.normalized.sorted.gfa" $h)" ] || cut_ok=0; done
[ "$(anchor_passes "$OUT/cutx")" = "ref:3 h1:4 h2:5 h3:2 " ] && [ "$cut_ok" = 1 ] \
  && ok "exact: shared copies cut unlike the interior's are widened over and fold (3/4/5/2)" \
  || bad "shared copies cut differently: passes '$(anchor_passes "$OUT/cutx")', expected 3/4/5/2, spelled kept $cut_ok"

echo
if [ "$fails" -eq 0 ]; then echo "panphorte_stats: all assertions passed"; exit 0; fi
echo "panphorte_stats: $fails assertion(s) FAILED"; exit 1
