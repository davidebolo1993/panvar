#!/usr/bin/env bash
# genotype_bounded_search.sh - acceptance gates for the BOUNDED-COMPLETE placement search.
#
#   genotype_bounded_search.sh <panvar-binary> <out-dir>
#
# WRITTEN BEFORE THE IMPLEMENTATION, deliberately. Every acceleration on this branch that was
# measured only after it existed turned out to be measuring something else: recruitment labels were
# not valid factor scopes, a projector assumed source->sink order, a "projection" of a named pair
# returned the called pair. The conditions below are frozen first so the search is built against
# them rather than described by them afterwards.
#
# WHY IT IS NEEDED, concretely: --reference-score enumerates every fragment start on both
# haplotypes and does not finish in ten minutes on c4's 226kb/252kb pair even for 76 fragments. The
# exact audit of the 76 fragments that decide c4/NA19240 is blocked on this, as is any trustworthy
# fragment linkage factor.
#
# THE CONTRACT:
#   1. COMPLETENESS. Every mate placement within the declared edit band is found. No occurrence cap,
#      no top-k. A cap makes completeness a function of a tuning parameter.
#   2. STATE FIDELITY. Strand, start, end, insert length and haplotype identity are preserved. A
#      placement is a state, not a position: the fragment model's state is (start, insert,
#      orientation) and collapsing any component changes the likelihood.
#   3. VALID FR ONLY, which is three conditions and not one: forward/reverse orientation, the
#      reverse mate DOWNSTREAM of the forward one, and an insert length inside the prior's support.
#      Orientation alone admits pairs the model assigns zero probability.
#   4. MULTIPLICITY PRESERVED MATHEMATICALLY, not as a record count.
#        * the same placement discovered through several anchors is ONE state;
#        * different repeat-copy coordinates are DISTINCT origins;
#        * identical-likelihood origins MAY be compressed to one representative plus
#          log(multiplicity), provided placement mass and exposure both stay exact.
#      Requiring a particular number of stored records would forbid precisely the compression LPA
#      needs. The gate compares LIKELIHOOD and EXPOSURE against exhaustive enumeration, never
#      storage. Collapsing distinct origins destroys copy-number information in a marginal placement
#      model and must change neither placement mass nor exposure.
#   5. EXHAUSTIVE AGREEMENT on small fixtures. "Unseeded" means no usable CURRENT syncmer seed --
#      NOT "no exact seed anywhere": under the pigeonhole construction a read within d edits must
#      contain an exact piece, ambiguous bases aside, so the latter would test an impossibility.
#
#      TWO STAGES, because the emissions differ and must not be conflated:
#
#      STAGE 1 -- HAMMING-COMPLETE. Pigeonhole recruitment plus exhaustive verification of candidate
#      starts, certified against the current reference. reference_emission is FIXED-POSITION
#      HAMMING: it counts mismatches at one offset and has no gap model (src/genotype_fragments.cpp,
#      reference_emission). It can therefore certify substitution, reverse strand, repeat, junction
#      and unseeded fixtures, and nothing else. This stage alone unblocks the c4 76-fragment audit,
#      which was run with --hamming-emission.
#
#      STAGE 2 -- INDEL-AWARE, and NOT gated by this file until it has its own oracle. With indels,
#      start, reference consumption, end coordinate and insert length all need a defined alignment
#      contract, and an exhaustive edit-distance oracle must be written against that same
#      definition. Activating an indel fixture against the Hamming reference would compare two
#      DIFFERENT EMISSIONS and call the difference an acceleration defect -- the exact
#      reference-versus-accelerated confusion this branch has already paid for twice.
#   6. REPORTED WORK: candidate starts, verified starts, placements kept, runtime, and the
#      omitted-mass bound. A search that cannot say what it skipped cannot be trusted not to skip.
#
# THE ACCEPTANCE CONDITIONS. An earlier draft demanded
#
#     bounded likelihood mass == exhaustive likelihood mass
#
# which is IMPOSSIBLE and would have made this gate unsatisfiable. The exhaustive reference
# integrates every start and every insert length and does NOT truncate at max_divergence --
# reference_emission is finite at every position, which is exactly why a union-of-spans scope came
# out locus-wide in the oracle work and had to be replaced by a counterfactual test. So out-of-band
# states carry small but NONZERO probability, and a search complete only within the band cannot
# equal the total. Requiring it would force the bounded search to evaluate the outside-band tail
# exactly, defeating its purpose, or would need a separate hard-band reference -- a DIFFERENT MODEL
# from the current exhaustive one, not a faster implementation of it.
#
# What must hold instead:
#
#   A. bounded IN-BAND mass == exhaustive IN-BAND mass        (exact, every in-band placement found)
#   B. bounded exposure     == exhaustive exposure            (exact, analytic)
#   C. the omitted mass is BOUNDED, and the bound is valid:
#
#          M_found <= M_exhaustive <= M_found + M_omitted_bound
#
#      in log space, log M_exhaustive lies in
#          [ log M_found , log_add(log M_found, log M_omitted_bound) ]
#
#   D. per fragment, per candidate AND per diploid pair -- a per-fragment bound does not imply the
#      pair bound, the same reason the scope oracle's guarantee had to be made diploid;
#   E. the exhaustive reference PAIR SCORE lies inside the reported interval;
#   F. TWO TOLERANCES, recorded separately, because they are different kinds of quantity:
#        * NUMERICAL, for the exact in-band equalities A and B -- floating-point agreement;
#        * APPROXIMATION, for the omitted mass in C -- expressed as NATS PER FRAGMENT or as an
#          omitted-mass fraction, never as one raw-likelihood threshold. A shared absolute threshold
#          is meaningless across loci whose totals differ by orders of magnitude: c4 pair scores run
#          near -1.1e6 while the deciding margin is ~500 nats.
#
#   G. CERTIFICATION USES THE INTERVALS, not the width alone:
#
#          certified(g)  iff  L(g) > max over other genotypes of U(other)
#
#      A narrow interval does not certify a winner whose interval still overlaps a competitor's.
#      Genotypes whose intervals overlap are reported as an equivalence set. At c4/NA19240 the
#      winner already changes between the pessimistic and optimistic tail treatments, so it fails
#      even the weaker same-winner test and "unresolved" is the correct output there.
#
# This is not a weaker gate. It stops the test both from demanding an impossible equality and from
# passing merely because omitted likelihood rounds to zero on fixtures small enough to hide it.
# Exact multiplicity-aware compression stays permitted internally: it must change neither the
# in-band mass nor the exposure.
#
# BUILD ORDER for stage 1, frozen with the contract so the sequence is not re-litigated mid-build:
#   1. complete SINGLE-MATE Hamming search -- split into d+1 pieces, index EVERY occurrence with no
#      cap, infer the fixed start from each occurrence, deduplicate starts, verify Hamming <= d;
#   2. both strands, and reads with no usable current syncmer seed;
#   3. combine mate placements into valid FR fragment states (orientation, reverse mate downstream,
#      insert inside the prior's support);
#   4. preserve repeat-origin multiplicity while deduplicating rediscoveries of the SAME state --
#      the same placement reached through several pieces is one state, two repeat copies are two;
#   5. exact in-band mass, analytic exposure, and the outside-band bound;
#   6. propagate those bounds through haplotypes, diploid pairs, and the background mixture -- the
#      bound must survive the mixture, since log[(1-eta)*lambda*M + eta*P_bg] is where it is used;
#   7. apply certification (G) and emit the equivalence set.
# Mutation-test every fixture BEFORE running c4. A fixture that cannot fail is not evidence, and
# five vacuous checks have already been caught on this branch by asking that question late.
#
# TWO DIAGNOSTIC QUANTITIES, kept in separate columns and never merged:
#   * BEST in-band placement -- what the current max-placement caller actually uses, and the
#     quantity that reclassifies c4's 23 competitor-only and 53 truth-only fragments;
#   * SUMMED in-band placement mass -- what the reference model requires and what any future
#     linkage factor consumes.
# They answer different questions. Sharing one column would make a max-placement result look like
# evidence about the marginal model, which is precisely the confusion that produced the retracted
# "collapsed multiplicity caused the c4 failures" claim.
#
# FIRST SCIENTIFIC RESULT, in this order: the c4 table of which of the 23 competitor-only and 53
# truth-only fragments remain EXCLUSIVE under complete in-band search, then the resulting lower and
# upper genotype intervals. Not "the search runs".
#
# STATUS: the search is NOT IMPLEMENTED. This file exits 77 (skip) until the entry point exists, so
# it can be registered now and start gating the moment there is something to gate. It must never be
# made to pass by weakening a condition.
set -uo pipefail
BIN="${1:?usage: genotype_bounded_search.sh <panvar> <outdir>}"
OUT="${2:?}"; mkdir -p "$OUT"

PY="${PYTHON:-python3}"
fails=0
ok()  { printf '  ok   %s\n' "$*"; }
bad() { printf '  FAIL %s\n' "$*"; fails=$((fails+1)); }

if ! "$BIN" genotype-frag --help 2>&1 | grep -q -- "--bounded-search"; then
  echo "SKIP: --bounded-search is not implemented yet; gates are frozen and waiting"
  echo "      acceptance: in-band mass and exposure EXACT; omitted mass validly BOUNDED, with"
  echo "                  log M_exhaustive in [log M_found, logadd(log M_found, log M_bound)],"
  echo "                  per fragment/candidate/pair; reference score inside the interval;"
  echo "                  interval width within tolerance before certifying. 6 fixture classes"
  echo "                  (substitution, indel, reverse strand, repeat, junction, unseeded);"
  echo "                  completeness within band, no cap/top-k, valid-FR join, exact"
  echo "                  multiplicity-aware compression allowed, work and omitted-mass reported"
  exit 77
fi

# ---- FIXTURES ---------------------------------------------------------------------------------
# SMALL by necessity, not convenience: the exhaustive reference is O(|hap| x |read|) per cell, and on
# real c4 (131 paths x 226kb) a ten-read comparison is ~1e11 character operations and does not
# finish. Small fixtures are what makes exhaustive agreement checkable at all.
OUT="$OUT/run.$$"; rm -rf "$OUT"; mkdir -p "$OUT"
seq_of() { awk -v n="$1" -v s="$2" 'BEGIN{x=s; b="ACGT";
           for(i=0;i<n;i++){x=(1103515245*x+12345)%2147483648; printf "%s", substr(b,(int(x/65536)%4)+1,1)} }'; }

# RPT is 500 bp, LONGER than the 350 bp insert, so a whole fragment fits inside one copy and hapA
# offers TWO fragment origins. At 200 bp no fragment fits and pairR_ had a single origin on both
# paths -- the multiplicity assertion was then about single mates only.
U=$(seq_of 400 11); RPT=$(seq_of 500 12); V=$(seq_of 400 13); W=$(seq_of 400 14)
# hapA: unique - repeat - unique - repeat - unique. Two identical copies, so a read inside the
# repeat has TWO distinct origins that must both survive; collapsing them is the copy-number loss.
HA="${U}${RPT}${V}${RPT}${W}"
HB="${U}${RPT}${V}"
# hapN carries AMBIGUOUS BASES. An N in the read meeting an N in the haplotype costs NO Hamming
# mismatch, so a read can sit at zero edits while every one of its d+1 pieces is spoiled by
# ambiguity -- the pigeonhole then proposes nothing. Ns in the READ ALONE cannot produce this: each
# would be a mismatch against an ACGT haplotype, so spoiling d+1 pieces needs more than d
# mismatches and the read is out of band anyway. Measured: that first version placed nowhere and the
# non-vacuity assertion caught it.
# NINE and SIXTEEN adjacent identical copies. A fragment inside one copy therefore has 9 and 16
# origins -- above any plausible top-2 or top-8 cap. The 1-vs-2 case cannot distinguish "keeps every
# origin" from "keeps the first eight", which is exactly the production failure this work targets
# (--placement-topk 8 is a real flag in this codebase).
HK9=""; for i in 1 2 3 4 5 6 7 8 9; do HK9="${HK9}${RPT}"; done
HK16=""; for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16; do HK16="${HK16}${RPT}"; done
HN=$("$PY" - "$(seq_of 500 15)" <<'PYEOF'
import sys
h=list(sys.argv[1].upper())
for x in range(10, len(h), 13): h[x]='N'
print("".join(h))
PYEOF
)
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\t1\t%s\n' "$HA"; printf 'S\t2\t%s\n' "$HB"; printf 'S\t3\t%s\n' "$HN"
  printf 'L\t1\t+\t2\t+\t0M\n'; printf 'L\t2\t+\t3\t+\t0M\n'
  printf 'S\t4\t%s\n' "$HK9"; printf 'S\t5\t%s\n' "$HK16"
  printf 'L\t3\t+\t4\t+\t0M\nL\t4\t+\t5\t+\t0M\n'
  printf 'P\tref\t1+\t*\n'; printf 'P\thapB\t2+\t*\n'; printf 'P\thapN\t3+\t*\n'
  printf 'P\thapK9\t4+\t*\n'; printf 'P\thapK16\t5+\t*\n'; } > "$OUT/g.gfa"
"$BIN" bubble -i "$OUT/g.gfa" -r ref -o "$OUT/b" --min-variant-bp 0 -q >/dev/null 2>&1

# A REAL PAIR: mate 1 forward at the fragment start, mate 2 reverse-complemented from the fragment
# END, so the implied insert is INSERT and lands inside the prior's support. The first version wrote
# both mates from the SAME 120 bp window, making every implied insert 120 against a support of
# [200,500] -- so no valid-FR state could form and the state comparison silently had nothing to
# compare. Measured: 0 cells, 0 states.
INSERT=350
emit_pair() {  # emit_pair <name> <hapseq> <1-based fragment start>
  local nm="$1" h="$2" st="$3"
  local m1 m2
  m1=$(printf '%s' "$h" | cut -c$((st))-$((st+119)))
  m2=$(printf '%s' "$h" | cut -c$((st+INSERT-120))-$((st+INSERT-1)) | rev | tr ACGTacgt TGCAtgca)
  printf '>%s/1\n%s\n>%s/2\n%s\n' "$nm" "$m1" "$nm" "$m2" >> "$OUT/reads.fa"
}
emit_pair150() {  # a 150+150 pair, so max(r1+r2) is 300 while min stays 240
  local nm="$1" h="$2" st="$3"
  local m1 m2
  m1=$(printf '%s' "$h" | cut -c$((st))-$((st+149)))
  m2=$(printf '%s' "$h" | cut -c$((st+INSERT-150))-$((st+INSERT-1)) | rev | tr ACGTacgt TGCAtgca)
  printf '>%s/1\n%s\n>%s/2\n%s\n' "$nm" "$m1" "$nm" "$m2" >> "$OUT/reads.fa"
}
emit_pairB() {  # the SAME fragment with mate roles swapped: mate 2 forward, mate 1 reverse
  local nm="$1" h="$2" st="$3"
  local f r
  f=$(printf '%s' "$h" | cut -c$((st))-$((st+119)))
  r=$(printf '%s' "$h" | cut -c$((st+INSERT-120))-$((st+INSERT-1)) | rev | tr ACGTacgt TGCAtgca)
  printf '>%s/1\n%s\n>%s/2\n%s\n' "$nm" "$r" "$nm" "$f" >> "$OUT/reads.fa"
}
emit() {  # single-mate fixtures: both records from one window, for the MATE-level comparison only
  printf '>%s/1\n%s\n' "$1" "$2" >> "$OUT/reads.fa"
  printf '>%s/2\n%s\n' "$1" "$(printf '%s' "$2" | rev | tr ACGTacgt TGCAtgca)" >> "$OUT/reads.fa"
}
: > "$OUT/reads.fa"
EXACT=$(printf '%s' "$HA" | cut -c101-220)                 # unique flank, exact
SUBST=$(printf '%s' "$HA" | cut -c301-420 | sed 's/^\(.\{10\}\)./\1N/;s/N/A/')   # one substitution
INREP=$(printf '%s' "$HA" | cut -c451-570)                 # inside the repeat: two origins in hapA
JUNCT=$(printf '%s' "$HA" | cut -c580-699)                 # spans repeat->unique junction
# PIGEONHOLE-CRITICAL read. The other fixtures are near-exact and are found even with too FEW
# pieces, so they cannot detect a broken pigeonhole -- measured: cutting d+1 to d left every
# assertion passing. This read carries exactly d mismatches positioned to spoil all d pieces of the
# mutant split while leaving one piece of the correct d+1 split clean. read=120, d=6:
#   correct: 7 pieces of 17 -> mismatches at 5,25,45,65,85,105 hit pieces 0,1,2,3,5,6; piece 4 CLEAN
#   mutant : 6 pieces of 20 -> the same positions hit all six; nothing is proposed and the
#            placement is LOST, which is exactly what the exhaustive comparison must catch.
# GENERATED from the production band, never hard-coded. mate_band_edits is floor(div*len)+1, so a
# 120 bp read at 5% gives d=7, not the 6 an open-coded floor() produces -- and a fixture built for
# the wrong d tests a narrower band than production accepts.
PIG=$("$PY" - "$HA" 120 0.05 <<'PYEOF'
import sys
h=sys.argv[1]; L=int(sys.argv[2]); div=float(sys.argv[3])
d=int(div*L)+1                      # mate_band_edits
npc, npm = d+1, d                   # correct split vs the mutant one
Pc, Pm = L//npc, L//npm
sub={'A':'C','C':'G','G':'T','T':'A','a':'C','c':'G','g':'T','t':'A'}
# one mismatch inside each MUTANT piece, chosen to leave at least one CORRECT piece clean
pos=[]
for i in range(npm):
    lo,hi=i*Pm,min((i+1)*Pm,L)-1
    # prefer a position whose correct-piece index is already used, so a correct piece stays clean
    best=None
    for x in range(lo,hi+1):
        c=min(x//Pc, npc-1)
        if best is None or (c in [min(y//Pc,npc-1) for y in pos]): best=x; break
    pos.append(best if best is not None else lo)
clean=set(range(npc))-{min(x//Pc,npc-1) for x in pos}
if not clean or len(pos)!=d:
    sys.stderr.write("fixture construction failed: d=%d clean=%s\n"%(d,sorted(clean))); sys.exit(3)
b=list(h[700:700+L].upper())
for x in pos: b[x]=sub[b[x]]
sys.stderr.write("d=%d correct=%dx%d mutant=%dx%d clean_correct_pieces=%s\n"%(d,npc,Pc,npm,Pm,sorted(clean)))
print("".join(b))
PYEOF
)
emit pigeon_ "$PIG"
# AMBIGUITY. An N in the read meeting an N in the haplotype costs no mismatch, so a read can be
# inside the band while every piece is spoiled by ambiguity. Such a read MUST take the exhaustive
# fallback; skipping the pieces without falling back loses the placement silently.
AMB=$(printf '%s' "$HN" | cut -c101-220)   # exact substring of hapN, Ns included: zero edits
emit amb_ "$AMB"
emit sub_    "$SUBST"
emit uniq_   "$EXACT"
emit rep_    "$INREP"
emit junc_   "$JUNCT"
# PAIRED fixtures, for the fragment-state comparison. pairU sits in unique sequence; pairR starts
# inside the first repeat copy so hapA offers TWO fragment origins and hapB one.
emit_pair  pairU_ "$HA" 101
emit_pair  pairR_ "$HA" 451     # wholly inside repeat copy 1 (401..900): two origins on hapA
emit_pairB pairB_ "$HA" 151     # mate roles swapped: exercises library orientation B
# MIXED MATE LENGTHS. With uniform 120+120 pairs the minimum and maximum combined length are both
# 240, so taking the min instead of the max is invisible. This 150+150 pair makes the maximum 300,
# and the shared lower support bound must follow the LONGEST fragment.
emit_pair150 pair150_ "$HA" 1301
# TAIL-ONLY. Both mates carry MORE than d mismatches at their true position, so no in-band state
# exists anywhere -- lower = -inf -- while the untruncated reference still has finite mass there.
# This is the case that tests whether the omitted bound bounds ANYTHING: with an empty in-band set
# the interval is [-inf, bound] and the reference must sit inside it.
TAILP=$("$PY" - "$HA" 120 0.05 <<'PYEOF'
import sys
h=sys.argv[1]; L=int(sys.argv[2]); div=float(sys.argv[3])
d=int(div*L)+1
sub={'A':'C','C':'G','G':'T','T':'A'}
def mut(seg, k):
    b=list(seg.upper())
    for i in range(k): b[(i*7+3) % len(b)]=sub[b[(i*7+3) % len(b)]]
    return "".join(b)
st=1001   # HA is 2200 bp; 2001+350 runs past the end and truncated mate 2, dropping the pair
m1=mut(h[st-1:st-1+L], d+4)                      # d+4 mismatches: comfortably out of band
tail=h[st-1+350-L:st-1+350]
rc=lambda s: s[::-1].translate(str.maketrans("ACGT","TGCA"))
m2=rc(mut(tail, d+4))
print(m1); print(m2)
PYEOF
)
printf '>tail_/1\n%s\n>tail_/2\n%s\n' "$(printf '%s' "$TAILP" | sed -n 1p)" \
  "$(printf '%s' "$TAILP" | sed -n 2p)" >> "$OUT/reads.fa"

# ---- A/B: exact agreement with the exhaustive scan ---------------------------------------------
if ! "$BIN" genotype-frag -i "$OUT/g.gfa" -b "$OUT/b" -o "$OUT/o" -R "$OUT/reads.fa" \
      --max-divergence 0.05 --fragment-len 350 --fragment-sd 50 \
      --bounded-search "$OUT/bs.tsv" -q >/dev/null 2>&1; then
  bad "the bounded search exited nonzero -- it disagreed with the exhaustive reference"
elif [ ! -s "$OUT/bs.tsv" ]; then
  bad "no bounded-search output"
else
  N=$(awk -F'\t' 'NR>1' "$OUT/bs.tsv" | wc -l | tr -d ' ')
  DIS=$(awk -F'\t' 'NR>1 && $7!="yes"' "$OUT/bs.tsv" | wc -l | tr -d ' ')
  [ "$N" -gt 0 ] && ok "compared $N (mate,strand,haplotype) cells against exhaustive" \
                 || bad "no cells compared"
  [ "$DIS" = 0 ] && ok "bounded == exhaustive on every cell (condition A)" \
                 || bad "$DIS cell(s) disagree with the exhaustive scan"
  # NON-VACUITY: placements must actually be found, or agreement is 0 == 0 everywhere.
  PL=$(awk -F'\t' 'NR>1{s+=$5} END{print s+0}' "$OUT/bs.tsv")
  [ "$PL" -gt 0 ] && ok "and the comparison is not vacuous ($PL placements found)" \
                  || bad "zero placements anywhere; agreement is 0==0 and proves nothing"
  # REPEAT MULTIPLICITY: a read inside the duplicated unit must have TWO origins on hapA and ONE on
  # hapB. This is the condition that distinguishes state dedup from multiplicity collapse.
  RA=$(awk -F'\t' 'NR>1 && index($1,"rep_")==1 && $4=="ref" && $5>0{print $5; exit}' "$OUT/bs.tsv")
  RB=$(awk -F'\t' 'NR>1 && index($1,"rep_")==1 && $4=="hapB" && $5>0{print $5; exit}' "$OUT/bs.tsv")
  [ "${RA:-0}" -ge 2 ] && ok "a read in the duplicated unit keeps BOTH origins on the 2-copy path ($RA)" \
                       || bad "the duplicated unit yielded ${RA:-0} origin(s); multiplicity collapsed"
  [ "${RB:-0}" = 1 ] && ok "and exactly one on the 1-copy path ($RB), so copy number is visible" \
                     || bad "the 1-copy path yielded ${RB:-0} origins, expected 1"
  # BOTH STRANDS exercised.
  # The pigeonhole-critical read MUST place, or the mutation test above is vacuous.
  PG=$(awk -F'\t' 'NR>1 && index($1,"pigeon_")==1 && $5>0' "$OUT/bs.tsv" | wc -l | tr -d ' ')
  PGE=$(awk -F'\t' 'NR>1 && index($1,"pigeon_")==1 && $8>=0{print $8; exit}' "$OUT/bs.tsv")
  [ "${PG:-0}" -gt 0 ] \
    && ok "the pigeonhole-critical read (d mismatches, one clean piece) places, at $PGE edits" \
    || bad "the pigeonhole-critical read did not place; a broken d+1 split would go undetected"
  # THE BAND MUST BE PRODUCTION'"'"'S. floor(div*len) is one edit narrower than mate_band_edits and
  # would certify a search that loses exactly the boundary placements.
  MB=$(awk -F'\t' 'NR>1 && index($1,"pigeon_")==1{print $9; exit}' "$OUT/bs.tsv")
  EXP=$("$PY" -c "print(int(0.05*120)+1)")
  [ "${MB:-0}" = "$EXP" ] \
    && ok "the gate uses the production band ($MB = mate_band_edits(0.05,120))" \
    || bad "the gate ran at max_edits=${MB:-?}, production uses $EXP; a narrower band certifies less"
  # AMBIGUOUS READS take the exhaustive fallback and still place.
  AF=$(awk -F'\t' 'NR>1 && index($1,"amb_")==1 && $14==1' "$OUT/bs.tsv" | wc -l | tr -d ' ')
  AP=$(awk -F'\t' 'NR>1 && index($1,"amb_")==1 && $4=="hapN" && $5>0' "$OUT/bs.tsv" | wc -l | tr -d ' ')
  [ "${AF:-0}" -gt 0 ] \
    && ok "an ambiguous read takes the exhaustive fallback ($AF cells)" \
    || bad "no ambiguous read fell back; N-spoiled pieces would lose placements silently"
  [ "${AP:-0}" -gt 0 ] \
    && ok "and still places ($AP cells), so the fallback is not merely skipping it" \
    || bad "the ambiguous read places nowhere; the fallback assertion above is vacuous"
  # THE VALID-FR RULE, on synthetic coordinates. Three conditions, not one: opposite strands, the
  # reverse mate downstream, and the insert inside the prior's support. Checked directly because a
  # state-set comparison can agree while both sides share the same wrong rule.
  if [ -s "$OUT/bs.tsv.fr.tsv" ]; then
    FRW=$(awk -F'\t' 'NR>1 && $7!=$8' "$OUT/bs.tsv.fr.tsv" | wc -l | tr -d ' ')
    FRN=$(awk -F'\t' 'NR>1' "$OUT/bs.tsv.fr.tsv" | wc -l | tr -d ' ')
    FRA=$(awk -F'\t' 'NR>1 && $8==1' "$OUT/bs.tsv.fr.tsv" | wc -l | tr -d ' ')
    FRR=$(awk -F'\t' 'NR>1 && $8==0' "$OUT/bs.tsv.fr.tsv" | wc -l | tr -d ' ')
    [ "${FRW:-1}" = 0 ] && ok "the valid-FR rule is right on all $FRN synthetic cases" \
                        || bad "$FRW valid-FR case(s) disagree with the expected verdict"
    { [ "${FRA:-0}" -gt 0 ] && [ "${FRR:-0}" -gt 0 ]; } \
      && ok "and covers both verdicts ($FRA accepted, $FRR rejected: below/above support, reverse upstream)" \
      || bad "the FR cases are one-sided ($FRA accept, $FRR reject); a constant rule would pass"
  else
    bad "no valid-FR case table was written"
  fi
  # FRAGMENT STATES: built independently from the bounded and exhaustive placement vectors and
  # compared as complete sets. Single-mate agreement cannot catch a lost library orientation, a
  # wrong combination of the four vectors, collapsed states, or lost haplotype identity.
  if [ ! -s "$OUT/bs.tsv.states.tsv" ]; then
    bad "no fragment-state table; enumerate_fragment_states never ran"
  else
    SN=$(awk -F'\t' 'NR>1{s+=$3} END{print s+0}' "$OUT/bs.tsv.states.tsv")
    SD=$(awk -F'\t' 'NR>1 && $5!="yes"' "$OUT/bs.tsv.states.tsv" | wc -l | tr -d ' ')
    SA=$(awk -F'\t' 'NR>1{s+=$6} END{print s+0}' "$OUT/bs.tsv.states.tsv")
    SB=$(awk -F'\t' 'NR>1{s+=$7} END{print s+0}' "$OUT/bs.tsv.states.tsv")
    # THE PRODUCTION INSERT SUPPORT. make_insert_prior uses 4 sigmas and floors lo at the mates'
    # combined length: 120 bp mates give [240,550]. A hard-coded mean +/- 3sd gives [200,500], which
    # both admits states the model rejects and omits valid ones between 501 and 550.
    ILO=$(awk -F'\t' 'NR==2{print $8}' "$OUT/bs.tsv.states.tsv")
    IHI=$(awk -F'\t' 'NR==2{print $9}' "$OUT/bs.tsv.states.tsv")
    # The reads mix 120+120 and 150+150 pairs, so max(r1+r2)=300 and min=240. Production floors lo
    # at the MAXIMUM, so the shared support is [300,550]; taking the minimum would give [240,550]
    # and this assertion is what tells them apart.
    { [ "${ILO:-0}" = 300 ] && [ "${IHI:-0}" = 550 ]; } \
      && ok "the shared insert support follows the LONGEST fragment ([$ILO,$IHI])" \
      || bad "insert support is [${ILO:-?},${IHI:-?}]; production floors lo at max(r1+r2)=300, not min=240"
    # HAPLOTYPE IDENTITY as a key, tested directly: two haplotypes' identical states must not
    # collapse. Per-haplotype vectors cannot show this -- each holds one haplotype value.
    if [ -s "$OUT/bs.tsv.hapkey.tsv" ]; then
      HU=$(awk -F'\t' 'NR==2{print $3}' "$OUT/bs.tsv.hapkey.tsv")
      H0=$(awk -F'\t' 'NR==2{print $1}' "$OUT/bs.tsv.hapkey.tsv")
      { [ "${H0:-0}" -gt 0 ] && [ "${HU:-0}" = $(( H0 * 2 )) ]; } \
        && ok "haplotype is part of the state key (2 x $H0 states stay $HU after unique)" \
        || bad "unioning two haplotypes identical states gave ${HU:-?} from 2 x ${H0:-?}; hap is not in the key"
    else
      bad "no haplotype-key table was written"
    fi
    # IN-BAND MASS, in the reference's own terms. Two relations, and they are different claims:
    #   * bounded == exhaustive: same states, same formula, so only FP associativity may separate;
    #   * bounded <= reference: the reference integrates EVERY start and insert and does not
    #     truncate at the divergence band, so the in-band sum cannot exceed the integral containing
    #     it. This caught a real parity bug -- the state mass omitted log(0.5) per orientation while
    #     the reference averages over them, so it sat exactly log(2) = 0.6931 ABOVE the reference on
    #     every cell.
    MC=$(awk -F'\t' 'NR>1{n++} END{print n+0}' "$OUT/bs.tsv.states.tsv")
    MD=$(awk -F'\t' 'NR>1 && $11!=$12{n++} END{print n+0}' "$OUT/bs.tsv.states.tsv")
    MV=$(awk -F'\t' 'NR>1 && $11>$13+1e-9{n++} END{print n+0}' "$OUT/bs.tsv.states.tsv")
    { [ "${MC:-0}" -gt 0 ] && [ "${MD:-1}" = 0 ]; } \
      && ok "bounded and exhaustive in-band MASS agree on all $MC cells" \
      || bad "$MD of ${MC:-0} cells differ in in-band mass"
    [ "${MV:-1}" = 0 ] \
      && ok "and no cell exceeds the untruncated reference (M_in_band <= M_reference)" \
      || bad "$MV cell(s) have in-band mass ABOVE the reference integral that contains it"
    # CONDITION C, the interval. Containment on every cell, and at least one cell with a POSITIVE
    # width -- otherwise the bound is being "satisfied" by a quantity that is always zero.
    IC=$(awk -F'\t' 'NR>1 && $16==1{n++} END{print n+0}' "$OUT/bs.tsv.states.tsv")
    IN=$(awk -F'\t' 'NR>1{n++} END{print n+0}' "$OUT/bs.tsv.states.tsv")
    IW=$(awk -F'\t' 'NR>1 && $17>0{n++} END{print n+0}' "$OUT/bs.tsv.states.tsv")
    { [ "${IN:-0}" -gt 0 ] && [ "${IC:-0}" = "${IN:-1}" ]; } \
      && ok "the exact reference lies inside [lower, upper] on all $IN cells" \
      || bad "the reference falls outside the interval on $(( IN - IC )) of $IN cells"
    [ "${IW:-0}" -gt 0 ] \
      && ok "and $IW cell(s) have a POSITIVE interval width, so the bound is not always zero" \
      || bad "every interval has zero width; the containment assertion above is vacuous"
    # TAIL-ONLY: every placement out of band, so the in-band set is EMPTY and lower = -inf, while
    # the untruncated reference still has finite mass. This is the case the bound exists for, and
    # emitting only cells with states hid it completely -- the fixture produced no rows at all.
    TN=$(awk -F'\t' 'NR>1 && index($1,"tail_")==1{n++} END{print n+0}' "$OUT/bs.tsv.states.tsv")
    TZ=$(awk -F'\t' 'NR>1 && index($1,"tail_")==1 && $3==0 && $11=="-inf" && $16==1{n++} END{print n+0}' "$OUT/bs.tsv.states.tsv")
    TF=$(awk -F'\t' 'NR>1 && index($1,"tail_")==1 && $13!="-inf" && $13<0{n++} END{print n+0}' "$OUT/bs.tsv.states.tsv")
    { [ "${TN:-0}" -gt 0 ] && [ "${TZ:-0}" = "${TN:-1}" ] && [ "${TF:-0}" = "${TN:-1}" ]; } \
      && ok "tail-only: $TN cells with an EMPTY in-band set, finite reference mass, all contained" \
      || bad "tail-only case is not exercised (rows=$TN empty+contained=$TZ finite_ref=$TF)"
    # CONDITION F, the adaptive tail. The flat bound charges "one mate at d+1" to every omitted
    # state and was 93 nats loose on the tail-only fragment -- an interval of [-inf,-40.45] against
    # a true reference of -133.03, which certifies nothing. Deepening the exact band to D and
    # bounding only beyond D recovers it exactly.
    AC=$(awk -F'\t' 'NR>1 && $24==1{n++} END{print n+0}' "$OUT/bs.tsv.states.tsv")
    AN=$(awk -F'\t' 'NR>1{n++} END{print n+0}' "$OUT/bs.tsv.states.tsv")
    [ "${AC:-0}" = "${AN:-1}" ] \
      && ok "the adaptive-tail interval contains the exact reference on all $AN cells" \
      || bad "$(( AN - AC )) cell(s) have the reference outside the adaptive interval"
    # The tightening must be demonstrated where the flat bound failed, not merely on easy cells.
    TW=$(awk -F'\t' 'NR>1 && index($1,"tail_")==1 && $2=="ref"{print $22; exit}' "$OUT/bs.tsv.states.tsv")
    TL=$(awk -F'\t' 'NR>1 && index($1,"tail_")==1 && $2=="ref"{print $11; exit}' "$OUT/bs.tsv.states.tsv")
    TT=$(awk -F'\t' 'NR>1 && index($1,"tail_")==1 && $2=="ref"{print $23; exit}' "$OUT/bs.tsv.states.tsv")
    { [ "${TL:-x}" = "-inf" ] && [ "${TT:-0}" = 1 ]; } \
      && ok "on the tail-only fragment the flat lower bound is -inf but the adaptive one is tight (width $TW)" \
      || bad "the tail-only tightening is not demonstrated (flat_lower=${TL:-?} tol_ok=${TT:-?})"
    # And a cell it CANNOT tighten must be reported uncertifiable, never certified on an unchecked
    # width. hapB/hapN hold the read far outside any reachable band.
    UN=$(awk -F'\t' 'NR>1 && index($1,"tail_")==1 && $2!="ref" && $23==0{n++} END{print n+0}' "$OUT/bs.tsv.states.tsv")
    [ "${UN:-0}" -gt 0 ] \
      && ok "and $UN cell(s) it cannot tighten are flagged WIDE-RAW-MASS, not silently accepted" \
      || bad "no cell is flagged wide; the tolerance check cannot be failing anywhere"
    # NAMING, deliberately. A wide RAW-MASS bound is not an uncertifiable GENOTYPE. Certification is
    # a statement about the final per-fragment diploid contribution,
    #     log((1-eta)*lambda*(M_a + M_b) + eta*P_bg),
    # and the two diverge in BOTH directions: a candidate with a hugely uncertain but tiny mass can
    # still give a narrow contribution once the background dominates, while two individually
    # acceptable candidate bounds can accumulate into a consequential genotype interval. Condition F
    # belongs on the contribution; this column stays a diagnostic until that layer exists.
    # MIXED MASS: substantial mass BOTH inside and outside the band, so upper = logadd(lower, bound)
  # is a genuinely non-trivial addition rather than lower + epsilon. At eps=0.001 the tail is
  # suppressed by ~8 nats per extra mismatch and can never be substantial; a HIGHER error rate makes
  # it so. This tests the bound's arithmetic, not a typical sequencing regime, and skipping it
  # because the usual regime is quiet would leave logadd untested. eps=0.15 was not enough -- the
  # bound still sat ~10 nats below the in-band mass and logadd moved the upper bound by 1e-4.
  "$BIN" genotype-frag -i "$OUT/g.gfa" -b "$OUT/b" -o "$OUT/mx" -R "$OUT/reads.fa" \
    --max-divergence 0.05 --fragment-len 350 --fragment-sd 50 --error-rate 0.30 \
    --bounded-search "$OUT/mx.tsv" -q >/dev/null 2>&1
  if [ ! -s "$OUT/mx.tsv.states.tsv" ]; then
    bad "the mixed-mass arm produced no states"
  else
    MXC=$(awk -F'\t' 'NR>1 && $16==1{n++} END{print n+0}' "$OUT/mx.tsv.states.tsv")
    MXN=$(awk -F'\t' 'NR>1{n++} END{print n+0}' "$OUT/mx.tsv.states.tsv")
    # a cell with BOTH finite in-band mass AND an omitted bound within 20 nats of it: the logadd
    # genuinely moves the upper bound rather than being swamped by one term.
    MXB=$(awk -F'\t' 'NR>1 && $11!="-inf" && $14!="-inf" && $17>0.05{n++} END{print n+0}' "$OUT/mx.tsv.states.tsv")
    [ "${MXC:-0}" = "${MXN:-1}" ] \
      && ok "mixed mass (eps=0.30): the reference stays inside the interval on all $MXN cells" \
      || bad "mixed mass: $(( MXN - MXC )) cell(s) fall outside the interval"
    [ "${MXB:-0}" -gt 0 ] \
      && ok "and $MXB cell(s) have in-band and omitted mass comparable, so logadd does real work" \
      || bad "no cell has comparable in-band and omitted mass; logadd is untested"
  fi

  # REFINEMENT INVARIANTS. Deepening D must only reallocate TAIL accounting.
    LM=$(awk -F'\t' 'NR>1 && $25==0{n++} END{print n+0}' "$OUT/bs.tsv.states.tsv")
    UM=$(awk -F'\t' 'NR>1 && $26==0{n++} END{print n+0}' "$OUT/bs.tsv.states.tsv")
    IS=$(awk -F'\t' 'NR>1 && $27==0{n++} END{print n+0}' "$OUT/bs.tsv.states.tsv")
    [ "${LM:-1}" = 0 ] && ok "exact mass never DECREASES as D grows" \
                       || bad "$LM cell(s) lost exact mass at a deeper D"
    [ "${UM:-1}" = 0 ] && ok "the upper bound never INCREASES as D grows" \
                       || bad "$UM cell(s) had the upper bound rise at a deeper D"
    # THE ONE THAT MATTERS MOST: the production-band mass must be IDENTICAL at every D. If it moves,
    # D is changing the in-band model rather than resolving its tail, and in-band numbers computed
    # at different depths are different quantities.
    [ "${IS:-1}" = 0 ] \
      && ok "the production-band mass is unchanged at every D (deepening only reallocates the tail)" \
      || bad "$IS cell(s) changed their production-band mass as D grew; D is altering the model"
    # OMITTED MULTIPLICITY, isolated to ONE edit class. reference-minus-in-band is the wrong
    # quantity: extra copies also add unrelated starts, junctions, insert lengths and orientations,
    # all with finite emission, so that difference measures everything the extra sequence brought --
    # and log-space subtraction is fragile near cancellation besides. So K states at the SAME
    # likelihood are summed directly: class_mass - log(count) must be IDENTICAL across different K,
    # which is exactly log M_K - log M_1 = log K.
    if [ ! -s "$OUT/bs.tsv.editclass.tsv" ]; then
      bad "no edit-class table; multiplicity is not isolated"
    else
      KM=$(awk -F'\t' 'NR>1 && $5>1{n++} END{print n+0}' "$OUT/bs.tsv.editclass.tsv")
      KU=$(awk -F'\t' 'NR>1{printf "%.6f\n",$7}' "$OUT/bs.tsv.editclass.tsv" | sort -u | wc -l | tr -d ' ')
      KV=$(awk -F'\t' 'NR>1{printf "%.6f\n",$7}' "$OUT/bs.tsv.editclass.tsv" | sort | uniq -c | awk '$1>1{n++} END{print n+0}')
      # ABOVE ANY CAP, and EXACT rather than ">=". The unit is 500 bp and each copy contributes
      # exactly one state at this class, so 9 copies span 8*500 = 4000 and 16 span 15*500 = 7500.
      # A count reached by extra starts inside the first eight copies would have the right TOTAL and
      # the wrong SPAN, which is why the span is asserted and not just the count.
      K9=$(awk -F'\t' 'NR>1 && $2=="hapK9" && $5>1{print $5; exit}' "$OUT/bs.tsv.editclass.tsv")
      K16=$(awk -F'\t' 'NR>1 && $2=="hapK16" && $5>1{print $5; exit}' "$OUT/bs.tsv.editclass.tsv")
      S9=$(awk -F'\t' 'NR>1 && $2=="hapK9" && $5>1{print $10; exit}' "$OUT/bs.tsv.editclass.tsv")
      S16=$(awk -F'\t' 'NR>1 && $2=="hapK16" && $5>1{print $10; exit}' "$OUT/bs.tsv.editclass.tsv")
      { [ "${K9:-0}" = 9 ] && [ "${S9:-0}" = 4000 ]; } \
        && ok "hapK9: exactly 9 states spanning all 9 copies (span $S9 = 8 x 500)" \
        || bad "hapK9 has ${K9:-0} states spanning ${S9:-0}; expected 9 spanning 4000"
      { [ "${K16:-0}" = 16 ] && [ "${S16:-0}" = 7500 ]; } \
        && ok "hapK16: exactly 16 states spanning all copies (span $S16 = 15 x 500), above any top-8 cap" \
        || bad "hapK16 has ${K16:-0} states spanning ${S16:-0}; expected 16 spanning 7500"
      [ "${KM:-0}" -gt 0 ] \
        && ok "multiplicity is present: $KM edit class(es) hold more than one distinct state" \
        || bad "no edit class has count > 1; the log K property cannot be tested"
      [ "${KV:-0}" -gt 0 ] \
        && ok "and per-state mass is IDENTICAL across different K, so class mass differs by log K" \
        || bad "no per-state mass is shared between a K=1 and a K>1 class; log K is unverified"
    fi
    # D-BOUNDARY, the terminal case: at D = read length every state is in band, so the residual is
    # empty, the bound is -inf, and lower == reference exactly.
    if [ ! -s "$OUT/bs.tsv.dboundary.tsv" ]; then
      bad "no D-boundary table"
    else
      DN=$(awk -F'\t' 'NR>1{n++} END{print n+0}' "$OUT/bs.tsv.dboundary.tsv")
      DA=$(awk -F'\t' 'NR>1 && $7==1{n++} END{print n+0}' "$OUT/bs.tsv.dboundary.tsv")
      DI=$(awk -F'\t' 'NR>1 && $3=="-inf"{n++} END{print n+0}' "$OUT/bs.tsv.dboundary.tsv")
      { [ "${DN:-0}" -gt 0 ] && [ "${DA:-0}" = "${DN:-1}" ]; } \
        && ok "at D = read length, lower == upper == exhaustive on all $DN cells" \
        || bad "$(( DN - DA )) cell(s) disagree at the D boundary"
      [ "${DI:-0}" = "${DN:-1}" ] \
        && ok "and the residual set is empty everywhere (bound = -inf)" \
        || bad "$(( DN - DI )) cell(s) still have a non-empty residual at D = read length"
    fi
    # EXPOSURE is analytic and depends on haplotype length, so two different haplotypes must differ.
    NE=$(awk -F'\t' 'NR>1{print $18}' "$OUT/bs.tsv.states.tsv" | sort -u | wc -l | tr -d ' ')
    [ "${NE:-0}" -ge 2 ] \
      && ok "analytic exposure varies with haplotype length ($NE distinct values)" \
      || bad "exposure has ${NE:-0} distinct value(s); it is not tracking haplotype length"
    [ "$SN" -gt 0 ] && ok "fragment states built ($SN) from bounded and exhaustive placements" \
                    || bad "zero fragment states; the comparison below is vacuous"
    [ "${SD:-1}" = 0 ] && ok "bounded and exhaustive fragment-state SETS are identical" \
                       || bad "$SD cell(s) have differing fragment-state sets"
    { [ "$SA" -gt 0 ] && [ "$SB" -gt 0 ]; } \
      && ok "both library orientations occur (A=$SA, B=$SB)" \
      || bad "only one orientation occurs (A=$SA, B=$SB); dropping the other would go undetected"
    # multiplicity at FRAGMENT level: the repeat is longer than the insert, so a fragment wholly
    # inside copy 1 has a twin inside copy 2 on the two-copy path and none on the one-copy path.
    MR=$(awk -F'\t' 'NR>1 && index($1,"pairR_")==1 && $2=="ref"{print $3; exit}' "$OUT/bs.tsv.states.tsv")
    MB=$(awk -F'\t' 'NR>1 && index($1,"pairR_")==1 && $2=="hapB"{print $3; exit}' "$OUT/bs.tsv.states.tsv")
    [ "${MR:-0}" -ge 2 ] \
      && ok "a fragment inside the duplicated unit keeps BOTH origins on the 2-copy path ($MR)" \
      || bad "the 2-copy path gave ${MR:-0} fragment origin(s); fragment-level multiplicity collapsed"
    [ "${MB:-0}" = 1 ] \
      && ok "and one on the 1-copy path ($MB), so copy number is visible at fragment level" \
      || bad "the 1-copy path gave ${MB:-0} fragment origins, expected 1"
  fi
  FS=$(awk -F'\t' 'NR>1 && $3=="+" && $5>0' "$OUT/bs.tsv" | wc -l | tr -d ' ')
  RS=$(awk -F'\t' 'NR>1 && $3=="-" && $5>0' "$OUT/bs.tsv" | wc -l | tr -d ' ')
  { [ "$FS" -gt 0 ] && [ "$RS" -gt 0 ]; } \
    && ok "both strands place ($FS forward, $RS reverse cells)" \
    || bad "only one strand places (fwd=$FS rev=$RS); the reverse case is untested"
  # WORK REPORTED, and the pigeonhole actually filtering rather than silently falling back.
  FB=$(awk -F'\t' 'NR>1 && $13==1' "$OUT/bs.tsv" | wc -l | tr -d ' ')
  VER=$(awk -F'\t' 'NR>1{s+=$12} END{print s+0}' "$OUT/bs.tsv")
  [ "$VER" -gt 0 ] && ok "verified starts reported ($VER), $FB cell(s) on the exhaustive fallback" \
                   || bad "no verified-start counts reported; the search cannot say what it skipped"
fi

echo
if [ "$fails" -eq 0 ]; then echo "bounded search: stage-1 A/B assertions passed"; else
  echo "bounded search: $fails assertion(s) failed"; fi
echo "ASSERTED: A/B (state sets, in-band mass), C (omitted-mass bound, tail-only case),"
echo "          adaptive tail tightening on RAW CANDIDATE MASS, and analytic exposure."
echo "NOT YET CONDITION F: F is a tolerance on the final per-fragment DIPLOID CONTRIBUTION"
echo "          log((1-eta)*lambda*(M_a+M_b) + eta*P_bg), not on raw log M. What is measured here"
echo "          is a raw-mass diagnostic; it neither implies nor is implied by a certifiable"
echo "          genotype."
echo "          Mixed mass (eps=0.30, logadd moves the upper bound by 0.066 nats), and the three"
echo "          refinement invariants: exact mass never decreases with D, the upper bound never"
echo "          increases, and the PRODUCTION-BAND mass is unchanged at every D."
echo "          Omitted multiplicity, isolated to one edit class (log K), and the D-boundary"
echo "          terminal case (residual empty, lower == upper == exhaustive)."
echo "NOT ASSERTED: diploid/background propagation, and G (certification by"
echo "              L(g) > max U(other) with overlapping intervals)."
exit "$fails"
