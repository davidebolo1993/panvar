#!/usr/bin/env bash
# genotype_frag_stats.sh - contract assertions for `panvar genotype-frag`.
#
# The fragment prototype shipped with no registered test at all, and three of the defects it has had
# so far would have been caught by one:
#
#   * a shortlist-wide anchor cap made the ANSWER depend on --max-haplotypes -- the same reads put the
#     truth at rank 2 with 48 haplotypes and rank 1 with 96;
#   * an unbounded Gaussian insert term let one mis-anchored pair cost 35,000 nats, so a sample's own
#     haplotype scored 15x worse than an unrelated one;
#   * the per-block projection sorted allele1/allele2 by index at every block, which is invisible to
#     an unordered allele-pair comparison and turns the spelled output into a chimera of the two
#     homologues.
#
# Each of those is an assertion below. Every expected value is derivable from the fixture rather than
# recorded from a previous run: reads are exact substrings of two known haplotypes, so the correct
# pair, the correct per-block alleles and the correct spelled sequences are all known by construction.
#
#   genotype_frag_stats.sh <panvar-binary> <out-dir>
set -uo pipefail

BIN="${1:?usage: genotype_frag_stats.sh <panvar> <outdir>}"
OUT="${2:?}"
mkdir -p "$OUT"; OUT="$OUT/run.$$.$(date +%s)"; rm -rf "$OUT"; mkdir -p "$OUT"
fails=0
ok()  { printf "  ok   %s\n" "$1"; }
bad_early() { printf "  FAIL %s\n" "$1"; exit 1; }
bad() { printf "  FAIL %s\n" "$1"; fails=$((fails + 1)); }

if ! "$BIN" genotype-frag --help >/dev/null 2>&1; then
  bad_early "this binary has no genotype-frag (configure with -DPANVAR_ENABLE_EXPERIMENTAL_GENOTYPE=ON); every assertion below would be vacuous"
fi

# Deterministic pseudo-random ACGT, same seeded LCG as genotype_stats.sh. The sequences must be fixed
# across runs or the syncmers, and therefore every recruitment and every alignment, would move.
seq_of() { awk -v n="$1" -v s="$2" 'BEGIN{x=s; b="ACGT";
           for(i=0;i<n;i++){x=(1103515245*x+12345)%2147483648; printf "%s", substr(b,(int(x/65536)%4)+1,1)} }'; }

# M is SHORT on purpose. The two bubbles must be close enough that a 350 bp fragment spans both, or
# hapAD/hapCB and hapAB/hapCD carry the identical diploid content {X1,X2,Y1,Y2} and are exactly tied
# -- no caller could separate them and the phase assertion below would be testing an unanswerable
# question. Verified: at M=900 the correct code fails, and it is right to. With M=100 the pairs are
# separated by physical read linkage across the junction, which is precisely the evidence production
# discards in count_reads and this architecture keeps.
L=$(seq_of 900 11); M=$(seq_of 100 12); N=$(seq_of 900 13)
X1=$(seq_of 400 21); X2=$(seq_of 400 22)
Y1=$(seq_of 400 31); Y2=$(seq_of 400 32)

# Two bubbles between three backbones. Four haplotypes spanning both combinations, so a projection
# that sorts allele indices per block can silently return a DIFFERENT panel haplotype pair that is
# also present -- which is exactly what makes the phase assertion below sharp rather than academic.
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\t1\t%s\n' "$L";  printf 'S\t2\t%s\n' "$X1"; printf 'S\t3\t%s\n' "$X2"
  printf 'S\t4\t%s\n' "$M";  printf 'S\t5\t%s\n' "$Y1"; printf 'S\t6\t%s\n' "$Y2"
  printf 'S\t7\t%s\n' "$N"
  for a in 2 3; do printf 'L\t1\t+\t%s\t+\t0M\nL\t%s\t+\t4\t+\t0M\n' "$a" "$a"; done
  for a in 5 6; do printf 'L\t4\t+\t%s\t+\t0M\nL\t%s\t+\t7\t+\t0M\n' "$a" "$a"; done
  printf 'P\tref\t1+,2+,4+,5+,7+\t*\n'
  for i in 1 2 3; do printf 'P\thapAB%d\t1+,2+,4+,5+,7+\t*\n' "$i"; done
  for i in 1 2 3; do printf 'P\thapCD%d\t1+,3+,4+,6+,7+\t*\n' "$i"; done
  for i in 1 2 3; do printf 'P\thapAD%d\t1+,2+,4+,6+,7+\t*\n' "$i"; done
  for i in 1 2 3; do printf 'P\thapCB%d\t1+,3+,4+,5+,7+\t*\n' "$i"; done
} > "$OUT/g.gfa"
"$BIN" bubble -i "$OUT/g.gfa" -r ref -o "$OUT/bub" --min-variant-bp 0 -q >/dev/null 2>&1

# The sample is hapAD / hapCB, NOT hapAB / hapCD, and the choice is the whole point of the phase
# assertion. Allele indices are assigned in first-seen path order, so bubble 1 is X1=0, X2=1 and
# bubble 2 is Y1=0, Y2=1. A hapAB/hapCD sample is (0,0) and (1,1): sorting the two alleles per block
# is a no-op on it and the assertion would pass against the defect it exists to catch -- verified,
# it did. A hapAD/hapCB sample is (0,1) and (1,0), so sorting yields (0,0) and (1,1), which spells
# hapAB/hapCD -- a real panel pair, wrong sample, and invisible to any unordered comparison.
H1="${L}${X1}${M}${Y2}${N}"     # hapAD*: allele 0 at bubble 1, allele 1 at bubble 2
H2="${L}${X2}${M}${Y1}${N}"     # hapCB*: allele 1 at bubble 1, allele 0 at bubble 2
printf '>truth1\n%s\n' "$H1" > "$OUT/truth1.fa"
printf '>truth2\n%s\n' "$H2" > "$OUT/truth2.fa"

# Proper FR pairs: mate 1 forward at i, mate 2 reverse-complemented from the far end of a 350 bp
# fragment. Both mates carry a name that differs only in the /1 and /2 suffix, which is what the
# fragment loader joins on.
emit_pairs() { awk -v s="$1" -v tag="$2" -v step="$3" 'BEGIN{
    rl=150; ins=350; n=length(s);
    for(i=1;i+ins-1<=n;i+=step){
      r1=substr(s,i,rl); r2=substr(s,i+ins-rl,rl);
      rc=""; for(j=length(r2);j>0;j--){c=substr(r2,j,1);
        rc=rc (c=="A"?"T":c=="C"?"G":c=="G"?"C":"A")}
      printf ">%s_%d/1\n%s\n>%s_%d/2\n%s\n", tag,i,r1,tag,i,rc } }'; }

emit_pairs "$H1" hA 25 >  "$OUT/reads.fa"
emit_pairs "$H2" hB 25 >> "$OUT/reads.fa"
NREADS=$(grep -c '>' "$OUT/reads.fa")
NFRAG=$((NREADS / 2))

run() { "$BIN" genotype-frag -i "$OUT/g.gfa" -b "$OUT/bub" -o "$1" -R "$OUT/reads.fa" \
          --haplotype-mode --fragment-len 350 --fragment-sd 50 "${@:2}" -q 2>/dev/null; }

# ------------------------------------------------------------------ mates are joined by read name
run "$OUT/base" -t 2
got=$("$BIN" genotype-frag -i "$OUT/g.gfa" -b "$OUT/bub" -o "$OUT/cnt" -R "$OUT/reads.fa" \
      --haplotype-mode -t 2 2>&1 | sed -nE 's/.*-> ([0-9]+) fragments \(([0-9]+) paired.*/\1 \2/p')
set -- $got
[ "${1:-0}" = "$NFRAG" ] && [ "${2:-0}" = "$NFRAG" ] \
  && ok "$NREADS reads join into $NFRAG fragments, all paired (the /1 and /2 suffix is stripped)" \
  || bad "expected $NFRAG paired fragments, got '${1:-}' fragments / '${2:-}' paired"

# Split R1/R2 files must give exactly what interleaved gives. Neither has to be declared, so a
# regression here would silently halve the evidence for every real cohort run.
awk '/\/1$/{p=1;print;next} /\/2$/{p=0;next} p' "$OUT/reads.fa" > "$OUT/r1.fa"
awk '/\/2$/{p=1;print;next} /\/1$/{p=0;next} p' "$OUT/reads.fa" > "$OUT/r2.fa"
split_frag=$("$BIN" genotype-frag -i "$OUT/g.gfa" -b "$OUT/bub" -o "$OUT/split" \
             -R "$OUT/r1.fa" -R "$OUT/r2.fa" --haplotype-mode -t 2 2>&1 \
             | sed -nE 's/.*-> ([0-9]+) fragments.*/\1/p')
[ "$split_frag" = "$NFRAG" ] \
  && ok "split R1/R2 files give the same $NFRAG fragments as interleaved" \
  || bad "split input gave $split_frag fragments, interleaved gave $NFRAG"

# ------------------------------------------------------------------------- the right pair is found
# Asserted on SEQUENCE, not on haplotype names. `ref` spells the same bases as hapAB* and the three
# hapAB copies spell the same bases as each other, so a name comparison would fail on a correct call
# for choosing a different label for identical sequence.
called_is_truth() {  # <prefix> -> "yes"/"no"
  "$BIN" genotype-frag -i "$OUT/g.gfa" -b "$OUT/bub" -o "$1.sp" \
    --spell-calls "$1.hap_blocks.tsv" -q >/dev/null 2>&1
  [ -s "$1.sp.called.fa" ] || { echo no; return; }
  local a b
  a=$(awk '/^>called_1/{f=1;next} /^>/{f=0} f' "$1.sp.called.fa" | tr -d '\n')
  b=$(awk '/^>called_2/{f=1;next} /^>/{f=0} f' "$1.sp.called.fa" | tr -d '\n')
  if { [ "$a" = "$H1" ] && [ "$b" = "$H2" ]; } || { [ "$a" = "$H2" ] && [ "$b" = "$H1" ]; }
  then echo yes; else echo no; fi
}
[ "$(called_is_truth "$OUT/base")" = "yes" ] \
  && ok "reads from hapAD and hapCB recover a pair spelling exactly those two haplotypes" \
  || bad "the called pair does not spell the two truth haplotypes"

# ----------------------------------------------------------------- PHASE: allele1 is one homologue
# The defect: sorting allele1/allele2 by index at every block. Both alternatives are IN THE PANEL
# here (hapAD, hapCB), so a sorted projection spells a real but wrong pair and nothing else would
# notice. Spelling the call and comparing to the truth catches it.
"$BIN" genotype-frag -i "$OUT/g.gfa" -b "$OUT/bub" -o "$OUT/spell" \
  --spell-calls "$OUT/base.hap_blocks.tsv" -q >/dev/null 2>&1
c1=$(awk '/^>called_1/{f=1;next} /^>/{f=0} f' "$OUT/spell.called.fa" 2>/dev/null | tr -d '\n')
if [ -n "$c1" ] && { [ "$c1" = "$H1" ] || [ "$c1" = "$H2" ]; }; then
  ok "called_1 is ONE homologue end to end (phase preserved across both bubbles)"
else
  bad "called_1 matches neither truth haplotype: allele1/allele2 are sorted per block, so concatenating them switches homologue mid-locus"
fi

# --------------------------------------------------- the answer does not depend on shortlist size
# The measured defect: a shortlist-WIDE anchor occurrence cap, which made anchoring a function of how
# many haplotypes were shortlisted. The cap is per haplotype now, and this is what holds it there.
# Sizes are all comfortably above the answer's needs. Smaller ones are NOT invariance failures: the
# panel holds four sequence-identical copies of each haplotype, so a shortlist of 4 can be filled
# entirely by duplicates of one of them and legitimately never see the other homologue. That is a
# candidate-generation property worth knowing -- it is a live lead for the gstm1 donors whose answer
# never reaches scoring -- but it is not what this assertion is about.
prev=""
for mh in 8 10 13; do
  run "$OUT/mh$mh" --max-haplotypes "$mh" -t 2
  cur=$(sed -n 2p "$OUT/mh$mh.hap_pairs.tsv" | cut -f2,3)
  if [ -n "$prev" ] && [ "$cur" != "$prev" ]; then
    bad "top pair changed with --max-haplotypes $mh: '$cur' vs '$prev' (anchoring must not depend on shortlist size)"
  fi
  prev="$cur"
done
[ -n "$prev" ] && ok "top pair is identical at --max-haplotypes 8, 10 and 13"

# The two mechanisms that make that invariance hold, asserted SEPARATELY -- with only four distinct
# sequences here, deduplication alone satisfies the test above, so it would pass with tie extension
# removed. It did.
run "$OUT/dedup" --max-haplotypes 13 -t 2
nuniq=$(( $(wc -l < "$OUT/dedup.hap_scores.tsv") - 1 ))
[ "$nuniq" = "4" ] \
  && ok "13 panel paths collapse to the 4 DISTINCT sequences; duplicates do not consume slots" \
  || bad "shortlist holds $nuniq entries, expected the 4 distinct sequences"

# Every haplotype here contains every allele the reads carry, so containment is 1.0 for all four and
# the cut is a tie. Asking for 2 must return all 4: cutting through a tie picks by path order.
run "$OUT/tie" --max-haplotypes 2 -t 2
ntie=$(( $(wc -l < "$OUT/tie.hap_scores.tsv") - 1 ))
[ "$ntie" = "4" ] \
  && ok "--max-haplotypes 2 against a 4-way containment tie keeps all 4, rather than cutting arbitrarily" \
  || bad "shortlist kept $ntie of a 4-way tie; the cut is being decided by path order"

# ---------------------------------------------------------------------------- thread invariance
run "$OUT/t1" -t 1
run "$OUT/t4" -t 4
if cmp -s "$OUT/t1.hap_blocks.tsv" "$OUT/t4.hap_blocks.tsv"; then
  ok "-t 1 and -t 4 produce byte-identical block calls"
else
  bad "threading changes the output; the fragment loop has a race or an order dependence"
fi

# ------------------------------------------------------- a wildly wrong insert is BOUNDED, not fatal
# Declaring a 20 bp library against a real 350 bp one makes every pair grossly discordant; the call
# must survive it.
#
# WHAT THIS DOES AND DOES NOT PROVE. Mutation-tested: replacing the bounded insert mixture with a raw
# Gaussian does NOT fail this assertion, and that is a fact about the model rather than a weak test.
# The per-fragment outlier component floors every fragment at read_ll(bg_divergence * len, len),
# which is at worst len*log(error_rate) ~ -1382 over a 300 bp fragment, and log-sum-exp then absorbs
# any insert penalty larger than that. So the outlier mixture ALREADY bounds the insert channel's
# worst case, and the insert mixture is a second layer whose independent effect this fixture cannot
# show. The insert term can still change RANKING among candidates that both sit above the floor,
# which is not tested here.
#
# Consequence for the record: the NA18939 repair reported in FRAGMENT_PROTOTYPE_RESULT.md changed two
# things at once, and on this evidence it is JOINT MATE PLACEMENT that was load-bearing -- mates no
# longer settling on different copies of a duplication -- not the insert bounding.
run "$OUT/badins" --fragment-len 20 --fragment-sd 5 -t 2
sc=$(sed -n 2p "$OUT/badins.hap_pairs.tsv" | cut -f4)
# Quantitative, because "finite" is not a bound: an unbounded Gaussian is finite too. Declaring a
# 20 bp library against a real 350 bp one puts every pair at z=(350-20)/5=66, so an unbounded term
# costs z^2/2 = 2178 nats PER FRAGMENT -- about -550,000 over this fixture. The bounded mixture
# floors each fragment at log(discordant_rate/discordant_span) = -13.8, so a few thousand in total.
# The threshold sits two orders of magnitude between them.
bounded=$(awk -v s="$sc" 'BEGIN{print (s==s+0 && s>-50000) ? "yes" : "no"}')
[ "$bounded" = "yes" ] \
  && ok "a grossly misdeclared insert costs a BOUNDED amount ($sc > -50000)" \
  || bad "score is $sc under a misdeclared insert; an unbounded Gaussian would land near -550000"
# THE OLD ASSERTION HERE WAS "and the correct pair still wins on sequence alone". It was wrong, and
# it passed only because of a defect. Declaring 20 +/- 5 against a real 350 bp library collapses the
# insert prior's support to the single length min_len, so NO fragment has a valid FR placement and
# every fragment falls to the background mixture. Verified against the exhaustive reference on the
# differential fixture: correctly declared it separates hA/hB -1797.0, hA/hC -9476.3, hE/hF -6269.0;
# misdeclared it returns -21754.822399868175 for all three, an exact tie. The MODEL has no
# discrimination left under a misdeclared insert.
#
# The accelerated scorer used to keep some, because it paired a placed mate with a synthetic partner
# at the band boundary -- a state the reference does not have. So this assertion was requiring the
# two scorers to DIFFER. What is true, and what is now asserted, is that the caller degrades into a
# large equivalence class rather than reporting a confident wrong pair.
# EVERY pair, not merely more than one. The reference ties all of them, so "> 1" would accept two
# tied pairs while the documented result is a complete tie -- a test claiming one thing and checking
# a much weaker one.
NTOP=$(awk -F'\t' 'NR>1 && $5==0' "$OUT/badins.hap_pairs.tsv" | wc -l | tr -d ' ')
NALL=$(( $(wc -l < "$OUT/badins.hap_pairs.tsv") - 1 ))
[ -n "$NTOP" ] && [ "$NTOP" = "$NALL" ] && [ "${NALL:-0}" -gt 1 ] \
  && ok "a misdeclared insert ties EVERY pair ($NTOP of $NALL), exactly as the exhaustive reference does" \
  || bad "a misdeclared insert tied only $NTOP of $NALL pairs; the reference ties all of them"
# ------------------------------------------------------ a third read under one name is not silent
{ cat "$OUT/reads.fa"; printf '>hA_1/1\n%s\n' "$(printf '%s' "$H1" | cut -c1-150)"; } > "$OUT/dup.fa"
dup=$("$BIN" genotype-frag -i "$OUT/g.gfa" -b "$OUT/bub" -o "$OUT/dup" -R "$OUT/dup.fa" \
      --haplotype-mode -t 2 2>&1 | grep -c "WARNING.*shared a name")
[ "$dup" -ge 1 ] \
  && ok "a third read under one fragment name is reported, not silently dropped" \
  || bad "a duplicate read name produced no warning; fragment counts could be wrong with no trace"

# ============================================================ joint depth: the duplication fixture
# The decisive case for one-fragment-one-placement. A haplotype carrying a segment TWICE, against one
# carrying it once, with reads generated from the single-copy sample. Every read from that segment has
# two valid placements on the two-copy candidate.
#
#   best-placement scoring pins all of them to one arbitrary copy, leaves the other falsely empty, and
#   then charges the candidate for absence its own heuristic invented;
#   complete assignment distributes them across both copies without counting any read twice.
#
# The assertion is not that the two-copy candidate loses -- it should, the sample has one copy -- but
# that the model refuses to run without an external depth rate, that the invariant holds, and that the
# answer does not depend on thread count.
D=$(seq_of 700 41)      # the duplicated segment
U=$(seq_of 500 42)      # unique spacer
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\t1\t%s\n' "$L"; printf 'S\t2\t%s\n' "$D"; printf 'S\t3\t%s\n' "$U"
  printf 'S\t4\t%s\n' "$D"; printf 'S\t5\t%s\n' "$N"
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t3\t+\t4\t+\t0M\nL\t4\t+\t5\t+\t0M\n'
  printf 'L\t3\t+\t5\t+\t0M\n'
  printf 'P\tref\t1+,2+,3+,5+\t*\n'
  for i in 1 2 3; do printf 'P\tone%d\t1+,2+,3+,5+\t*\n' "$i"; done
  for i in 1 2 3; do printf 'P\ttwo%d\t1+,2+,3+,4+,5+\t*\n' "$i"; done
} > "$OUT/dup.gfa"
"$BIN" bubble -i "$OUT/dup.gfa" -r ref -o "$OUT/dupbub" --min-variant-bp 0 -q >/dev/null 2>&1
DH="${L}${D}${U}${N}"
emit_pairs "$DH" dupA 25 >  "$OUT/dupreads.fa"
emit_pairs "$DH" dupB 25 >> "$OUT/dupreads.fa"

if "$BIN" genotype-frag -i "$OUT/dup.gfa" -b "$OUT/dupbub" -o "$OUT/dupj" -R "$OUT/dupreads.fa" \
     --haplotype-mode --joint-depth -q >/dev/null 2>&1; then
  bad "--joint-depth ran without --haploid-depth; the rate must come from outside the candidate set"
else
  ok "--joint-depth refuses to run without an externally supplied --haploid-depth"
fi

jrun() { "$BIN" genotype-frag -i "$OUT/dup.gfa" -b "$OUT/dupbub" -o "$1" -R "$OUT/dupreads.fa" \
           --haplotype-mode --joint-depth --haploid-depth 0.05 --joint-top-pairs 0 \
           --fragment-len 350 --fragment-sd 50 "${@:2}" -q 2>&1; }
if jerr=$(jrun "$OUT/dj1" -t 1) && [ -s "$OUT/dj1.hap_pairs.tsv" ]; then
  ok "joint depth completes on the duplication fixture (window/assignment invariant held)"
  jrun "$OUT/dj4" -t 4 >/dev/null 2>&1
  if cmp -s "$OUT/dj1.hap_blocks.tsv" "$OUT/dj4.hap_blocks.tsv"; then
    ok "joint depth is identical at -t 1 and -t 4"
  else bad "joint depth differs across thread counts"; fi
  # The single-copy sample must not be called as the two-copy haplotype. NOTE, measured: alignment
  # alone already gets this right on this fixture, so this assertion does not by itself demonstrate
  # the assignment mechanism -- it guards against the joint model BREAKING a case alignment handles.
  # What the fixture does test of the mechanism is the invariant above: with every read having two
  # valid placements on the two-copy candidate, window counts still sum to exactly the number of
  # assigned fragments, so no read was counted twice.
  top=$(sed -n 2p "$OUT/dj1.hap_pairs.tsv" | cut -f2,3)
  case "$top" in
    *two*) bad "a one-copy sample was called as the two-copy haplotype: $top" ;;
    *) ok "a one-copy sample is not called as the two-copy haplotype" ;;
  esac
else
  bad "joint depth failed on the duplication fixture: $jerr"
fi

# ------------------------------------------- homozygous intensity: the log(2) per placed fragment
# A homozygous pair is two copies of one sequence, so its exposure doubles AND its event intensity
# doubles. Doubling only the exposure charges the pair for coverage it is not credited with -- log(2)
# per placed fragment, thousands of nats -- and silently penalises every homozygous call.
#
# The fixture is a genuinely homozygous sample: reads from ONE haplotype at full depth. The correct
# call is that haplotype against itself, and under the bug it loses to a heterozygous pair.
emit_pairs "$H1" homo 12 > "$OUT/homo.fa"
if "$BIN" genotype-frag -i "$OUT/g.gfa" -b "$OUT/bub" -o "$OUT/homo" -R "$OUT/homo.fa" \
     --haplotype-mode --joint-marginal --haploid-depth 0.05 --joint-top-pairs 0 \
     --fragment-len 350 --fragment-sd 50 -t 2 -q >/dev/null 2>&1; then
  hs=$(awk -F'\t' 'NR==2{print ($2==$3) ? "hom" : "het"}' "$OUT/homo.hap_pairs.tsv")
  [ "$hs" = "hom" ] \
    && ok "a homozygous sample is called homozygous under --joint-marginal (intensity doubled with exposure)" \
    || bad "a homozygous sample was called heterozygous: the homozygous intensity is missing its factor of two"
else
  bad "--joint-marginal failed on the homozygous fixture"
fi

# ---------------------------------------------- same-strand mates must not form a fragment
# A real correctness defect, not test scaffolding: mate placements were combined without checking
# orientation, so two same-strand reads could be joined into a "fragment" the library cannot produce.
#
# The fixture emits the SAME pairs twice: once correctly (mate 2 reverse-complemented) and once with
# mate 2 left forward. A valid FR pair must score better than a same-strand one, because the
# same-strand pair should earn no concordant-fragment evidence at all.
emit_same_strand() { awk -v s="$1" -v tag="$2" -v step="$3" 'BEGIN{
    rl=150; ins=350; n=length(s);
    for(i=1;i+ins-1<=n;i+=step){
      printf ">%s_%d/1\n%s\n>%s_%d/2\n%s\n", tag,i,substr(s,i,rl),tag,i,substr(s,i+ins-rl,rl) } }'; }
emit_same_strand "$H1" ss 25 > "$OUT/ss_bad.fa"
emit_pairs      "$H1" ss 25 > "$OUT/ss_good.fa"
# Scored through the ACCELERATED path, because that is where the FR guard lives. Checking this with
# --reference-score would pass even if the fast guard were deleted, which is the defect the fixture
# exists for.
fastscore() {   # <reads> -> the top pair's score under the accelerated scorer
  "$BIN" genotype-frag -i "$OUT/g.gfa" -b "$OUT/bub" -o "$OUT/fs" -R "$1" \
    --haplotype-mode --rung-zero --hamming-emission --top-pairs 1 \
    --haploid-depth 0.05 --fragment-len 350 --fragment-sd 50 -t 2 -q >/dev/null 2>&1
  sed -n 2p "$OUT/fs.hap_pairs.tsv" | cut -f4
}
SG=$(fastscore "$OUT/ss_good.fa"); SB=$(fastscore "$OUT/ss_bad.fa")
awk -v a="${SG:-0}" -v b="${SB:-0}" 'BEGIN{ exit !(a > b + 1) }' \
  && ok "ACCELERATED: same-strand mates earn less than a valid FR pair ($SG vs $SB)" \
  || bad "accelerated scorer credits same-strand mates like FR pairs: $SG vs $SB"

# ------------------------------------------- strand prior: a haplotype and its revcomp score alike
# The 1/2-per-strand factor and the both-orientations search must agree, and the scorer must be
# strand-symmetric. Not cosmetic: real panel sequence contains haplotypes spelled reverse-complemented
# by block concatenation, which is what made an earlier recall measurement read 48% instead of 100%.
# Strand symmetry through the ACCELERATED path too: build a second graph whose node sequences are
# reverse-complemented, so every path spells the reverse complement of its counterpart. The same reads
# must score identically against it.
rcof() { printf '%s' "$1" | tr 'ACGTacgt' 'TGCAtgca' | rev; }
# The node ORDER reverses as well as each node's sequence: rc(A B C) is rc(C) rc(B) rc(A), not
# rc(A) rc(B) rc(C). Reverse-complementing nodes in place spells a scrambled sequence, which is what
# the first version of this fixture did -- and it looked exactly like an asymmetric scorer.
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\t1\t%s\n' "$(rcof "$N")"
  printf 'S\t2\t%s\n' "$(rcof "$Y1")"; printf 'S\t3\t%s\n' "$(rcof "$Y2")"
  printf 'S\t4\t%s\n' "$(rcof "$M")"
  printf 'S\t5\t%s\n' "$(rcof "$X1")"; printf 'S\t6\t%s\n' "$(rcof "$X2")"
  printf 'S\t7\t%s\n' "$(rcof "$L")"
  for a in 2 3; do printf 'L\t1\t+\t%s\t+\t0M\nL\t%s\t+\t4\t+\t0M\n' "$a" "$a"; done
  for a in 5 6; do printf 'L\t4\t+\t%s\t+\t0M\nL\t%s\t+\t7\t+\t0M\n' "$a" "$a"; done
  # rc(L X1 M Y1 N) = rc(N) rc(Y1) rc(M) rc(X1) rc(L) = 1,2,4,5,7
  printf 'P\tref\t1+,2+,4+,5+,7+\t*\n'
  for i in 1 2 3; do printf 'P\thapAB%d\t1+,2+,4+,5+,7+\t*\n' "$i"; done
  for i in 1 2 3; do printf 'P\thapCD%d\t1+,3+,4+,6+,7+\t*\n' "$i"; done
  for i in 1 2 3; do printf 'P\thapAD%d\t1+,3+,4+,5+,7+\t*\n' "$i"; done
  for i in 1 2 3; do printf 'P\thapCB%d\t1+,2+,4+,6+,7+\t*\n' "$i"; done
} > "$OUT/grc.gfa"
"$BIN" bubble -i "$OUT/grc.gfa" -r ref -o "$OUT/rcbub" --min-variant-bp 0 -q >/dev/null 2>&1
"$BIN" genotype-frag -i "$OUT/grc.gfa" -b "$OUT/rcbub" -o "$OUT/frc" -R "$OUT/reads.fa" \
  --haplotype-mode --rung-zero --hamming-emission --top-pairs 1 \
  --haploid-depth 0.05 --fragment-len 350 --fragment-sd 50 -t 2 -q >/dev/null 2>&1
"$BIN" genotype-frag -i "$OUT/g.gfa" -b "$OUT/bub" -o "$OUT/ffw" -R "$OUT/reads.fa" \
  --haplotype-mode --rung-zero --hamming-emission --top-pairs 1 \
  --haploid-depth 0.05 --fragment-len 350 --fragment-sd 50 -t 2 -q >/dev/null 2>&1
FS=$(sed -n 2p "$OUT/ffw.hap_pairs.tsv" | cut -f4)
RS=$(sed -n 2p "$OUT/frc.hap_pairs.tsv" | cut -f4)
awk -v a="${FS:-0}" -v b="${RS:-1}" 'BEGIN{ d=a-b; if(d<0) d=-d; exit !(d < 0.5) }' \
  && ok "ACCELERATED: a reverse-complemented panel scores identically ($FS vs $RS)" \
  || bad "accelerated strand handling is asymmetric: $FS vs $RS"

# ------------------------------------------------------- the post-exclusion sequence dump
# The LPA benchmark spent a whole branch charging a representation difference to the genotyper: the
# caller scored haplotypes spelled from a FOLDED graph while the truth had been spelled from the
# unfolded stage, ~4000 edits apart on a 300 kb locus, and every accuracy number carried that drift.
# The round-trip invariant could not see it -- it checks the block spelling against the GFA path and
# never the GFA path against the assembly. --dump-scored-sequences is what makes the comparison
# possible, so what it emits has to be the bytes the scorer uses and not a second spelling of them.
md5_of() { if command -v md5sum >/dev/null 2>&1; then md5sum | cut -d' ' -f1;
           else md5 -q; fi; }
PYTHON_BIN="${PYTHON:-python3}"
command -v "$PYTHON_BIN" >/dev/null 2>&1 || PYTHON_BIN=""

# Standalone: no --reads at all, because the audit case is a graph and a panel and nothing else.
"$BIN" genotype-frag -i "$OUT/g.gfa" -b "$OUT/bub" -o "$OUT/dmp" \
  --dump-scored-sequences "$OUT/dmp" -q >/dev/null 2>&1
if [ -s "$OUT/dmp.scored_sequences.tsv" ] && [ -s "$OUT/dmp.scored_sequences.fa" ]; then
  ok "--dump-scored-sequences runs with no --reads"
else
  bad "--dump-scored-sequences wrote nothing without --reads"
fi

# 13 paths: ref + three each of hapAB/hapCD/hapAD/hapCB.
NROW=$(( $(wc -l < "$OUT/dmp.scored_sequences.tsv") - 1 ))
[ "$NROW" = 13 ] && ok "dump covers all 13 panel paths" \
                 || bad "dump covers $NROW paths, expected 13"

# Every path must round-trip against its own GFA spelling, and the column must SAY so -- a dump that
# reported the same md5 in both columns by construction would assert nothing.
NRT=$(awk -F'\t' 'NR>1 && $8=="yes"' "$OUT/dmp.scored_sequences.tsv" | wc -l | tr -d ' ')
[ "$NRT" = 13 ] && ok "all 13 round-trip: block spelling == GFA path spelling" \
                || bad "only $NRT of 13 round-trip against the GFA path"

# The bytes are the fixture's, checked against the system md5 rather than against the caller's own
# other column. hapAD* spells H1 by construction; this is the assertion that would have caught LPA.
WANT=$(printf '%s' "$H1" | md5_of)
GOT=$(awk -F'\t' '$2=="hapAD1"{print $4}' "$OUT/dmp.scored_sequences.tsv")
[ "$GOT" = "$WANT" ] && ok "dumped hapAD1 md5 matches the known haplotype ($WANT)" \
                     || bad "dumped hapAD1 md5 $GOT, the fixture haplotype is $WANT"
LEN=$(awk -F'\t' '$2=="hapAD1"{print $3}' "$OUT/dmp.scored_sequences.tsv")
[ "$LEN" = "${#H1}" ] && ok "dumped hapAD1 length matches (${#H1} bp)" \
                      || bad "dumped hapAD1 is $LEN bp, the fixture haplotype is ${#H1} bp"

# POST-exclusion, which is the word in the flag's name. An excluded path must leave the panel group
# and reappear under held_out still spelling its own sequence -- if the dump ignored --exclude-
# haplotypes it would show the panel the caller does not score, which is the failure mode it exists
# to rule out.
"$BIN" genotype-frag -i "$OUT/g.gfa" -b "$OUT/bub" -o "$OUT/dmpx" \
  --dump-scored-sequences "$OUT/dmpx" --exclude-haplotypes 'hapAD1,hapCB1' -q >/dev/null 2>&1
NPAN=$(awk -F'\t' 'NR>1 && $1=="panel"' "$OUT/dmpx.scored_sequences.tsv" | wc -l | tr -d ' ')
NHLD=$(awk -F'\t' 'NR>1 && $1=="held_out"' "$OUT/dmpx.scored_sequences.tsv" | wc -l | tr -d ' ')
{ [ "$NPAN" = 11 ] && [ "$NHLD" = 2 ]; } \
  && ok "exclusion is reflected in the dump: 11 panel, 2 held out" \
  || bad "after excluding two paths the dump shows $NPAN panel / $NHLD held out, expected 11 / 2"
grep -q "^panel"$'\t'"hapAD1"$'\t' "$OUT/dmpx.scored_sequences.tsv" \
  && bad "an excluded path is still listed in the panel group" \
  || ok "an excluded path is gone from the panel group"
XGOT=$(awk -F'\t' '$1=="held_out" && $2=="hapAD1"{print $4}' "$OUT/dmpx.scored_sequences.tsv")
[ "$XGOT" = "$WANT" ] && ok "a held-out path still spells its own sequence from its own blocks" \
                      || bad "held-out hapAD1 md5 $XGOT, expected $WANT"

# The FASTA and the TSV must agree, or the hash travels while the sequence does not.
FGOT=""
[ -n "$PYTHON_BIN" ] && FGOT=$("$PYTHON_BIN" - "$OUT/dmp.scored_sequences.fa" <<'PYEOF' 2>/dev/null || true
import sys, hashlib
keep, buf = False, []
for line in open(sys.argv[1]):
    if line[0] == '>':
        if keep: break
        keep = line[1:].split()[0] == 'hapAD1'
    elif keep: buf.append(line.strip())
print(hashlib.md5("".join(buf).encode()).hexdigest())
PYEOF
)
if [ -n "$FGOT" ]; then
  [ "$FGOT" = "$WANT" ] && ok "the dumped FASTA record matches the md5 in the TSV" \
                        || bad "FASTA record hashes to $FGOT, the TSV says $WANT"
fi

echo
if [ "$fails" -eq 0 ]; then echo "genotype-frag stats: all assertions passed"; else
  echo "genotype-frag stats: $fails assertion(s) failed"; fi
exit "$fails"
