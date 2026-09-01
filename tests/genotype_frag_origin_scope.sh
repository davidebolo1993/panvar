#!/usr/bin/env bash
# genotype_frag_origin_scope.sh - the placement-dependency oracle, on the seven cases that define it.
#
#   genotype_frag_origin_scope.sh <panvar-binary> <out-dir>
#
# WHY THESE SEVEN. The blocker (genotype_frag_factorisation.sh) showed that recruitment scoping is
# not a valid factorisation: a fragment recruited to one block can hold likelihood mass at another.
# What replaces recruitment is the fragment's ORIGIN UNIVERSE and the scope derived from it -- the
# union of block variables able to change SUM_z P(f,z|G).
#
# A correct likelihood SUMS over origins, so an N-copy array offering ~N origins is EVIDENCE, not a
# bug. These fixtures therefore test that the oracle preserves calibrated MASS and reports honest
# SCOPE; they do not test that origin counts are small. Each case is built so a wrong answer is
# possible: a fixture whose expected scope is whatever the code happens to produce tests nothing.
#
#   1 unique local           genuinely unary                     -> scope 1 block
#   2 adjacent boundary      transition scope                    -> scope 2 CONSECUTIVE blocks
#   3 opposite phases        identical unordered alleles         -> universe is candidate-independent
#   4 repeat within a block  multiplicity WITHOUT added scope    -> many origins, scope still 1
#   5 non-adjacent blocks    shared sequence in distant blocks   -> scope NON-consecutive
#   6 local + remote         the blocker's exact failure         -> scope > 1 for a `local` fragment
#   7 unequal copy number    mass and exposure move together     -> more copies, more origin mass
set -uo pipefail
BIN="${1:?usage: genotype_frag_origin_scope.sh <panvar> <outdir>}"
OUT="${2:?}"; mkdir -p "$OUT"; OUT="$OUT/run.$$"; rm -rf "$OUT"; mkdir -p "$OUT"
PY="${PYTHON:-python3}"
fails=0
ok()  { printf "  ok   %s\n" "$1"; }
bad() { printf "  FAIL %s\n" "$1"; fails=$((fails+1)); }

seq_of() { awk -v n="$1" -v s="$2" 'BEGIN{x=s; b="ACGT";
           for(i=0;i<n;i++){x=(1103515245*x+12345)%2147483648; printf "%s", substr(b,(int(x/65536)%4)+1,1)} }'; }

# Reverse complement + read-pair writer, via a QUOTED heredoc. An inline $PY -c "..." lets the shell
# expand {'A':'T',...} as a brace list, the dict becomes a syntax error, a malformed FASTA is written
# and the assertion then measures something other than what it names. Measured: that is what the
# first version of case 1 did.
mkpair() {   # mkpair <seq> <name> <start> <readlen> <insert> <outfile>
  "$PY" - "$1" "$2" "$3" "$4" "$5" "$6" <<'MKPAIR_PY'
import sys
s, name, st, rl, ins, out = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), sys.argv[6]
comp = {'A':'T','C':'G','G':'C','T':'A','N':'N'}
r1 = s[st:st+rl]
r2 = s[st+ins-rl:st+ins]
rc = ''.join(comp[c] for c in reversed(r2))
with open(out, 'a') as fh:
    fh.write(">%s/1\n%s\n>%s/2\n%s\n" % (name, r1, name, rc))
MKPAIR_PY
}

# reads: exact substrings, so every expected origin is reachable and nothing depends on error model
mkreads() { awk -v s="$1" -v tag="$2" -v st="$3" 'BEGIN{ rl=40; ins=100; n=length(s);
    for(i=1;i+ins-1<=n;i+=st){ r1=substr(s,i,rl); r2=substr(s,i+ins-rl,rl);
      rc=""; for(j=length(r2);j>0;j--){c=substr(r2,j,1);
        rc=rc (c=="A"?"T":c=="C"?"G":c=="G"?"C":"A")}
      printf ">%s_%d/1\n%s\n>%s_%d/2\n%s\n", tag,i,r1,tag,i,rc } }'; }

CF=(--error-rate 0.005 --fragment-len 100 --fragment-sd 12 --placement-topk 8)

# scope_of <universe.tsv> <fragment-prefix>  -> "<n_blocks> <scope>"
scope_of() { awk -F'\t' -v p="$2" 'NR>1 && index($1,p)==1 {print $7" "$8; exit}' "$1"; }
lse_of()   { awk -F'\t' -v p="$2" 'NR>1 && index($1,p)==1 {print $4; exit}' "$1"; }
orig_of()  { awk -F'\t' -v p="$2" 'NR>1 && index($1,p)==1 {print $2; exit}' "$1"; }
omit_of()  { awk -F'\t' -v p="$2" 'NR>1 && index($1,p)==1 {print $6; exit}' "$1"; }

# ---------------------------------------------------------------------------------------------
# 1, 2, 3, 6 share one fixture: two bubbles close enough that 100 bp fragments span the boundary,
# with UNIQUE flanking sequence so a local fragment really is local.
L=$(seq_of 300 11); M=$(seq_of 30 12); N=$(seq_of 300 13)
X1=$(seq_of 120 21); X2=$(seq_of 120 22); Y1=$(seq_of 120 31); Y2=$(seq_of 120 32)
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\t1\t%s\n' "$L"; printf 'S\t2\t%s\n' "$X1"; printf 'S\t3\t%s\n' "$X2"
  printf 'S\t4\t%s\n' "$M"; printf 'S\t5\t%s\n' "$Y1"; printf 'S\t6\t%s\n' "$Y2"
  printf 'S\t7\t%s\n' "$N"
  for a in 2 3; do printf 'L\t1\t+\t%s\t+\t0M\nL\t%s\t+\t4\t+\t0M\n' "$a" "$a"; done
  for a in 5 6; do printf 'L\t4\t+\t%s\t+\t0M\nL\t%s\t+\t7\t+\t0M\n' "$a" "$a"; done
  printf 'P\tref\t1+,2+,4+,5+,7+\t*\n'
  printf 'P\tpAB\t1+,2+,4+,5+,7+\t*\n'; printf 'P\tpCD\t1+,3+,4+,6+,7+\t*\n'
  printf 'P\tpAD\t1+,2+,4+,6+,7+\t*\n'; printf 'P\tpCB\t1+,3+,4+,5+,7+\t*\n'
} > "$OUT/g.gfa"
"$BIN" bubble -i "$OUT/g.gfa" -r ref -o "$OUT/b" --min-variant-bp 0 -q >/dev/null 2>&1

# a fragment wholly inside the LEFT flank (unique, far from any bubble) -> case 1
: > "$OUT/r1.fa"
H_AD="${L}${X1}${M}${Y2}${N}"
mkpair "$H_AD" uniq_1 20 40 100 "$OUT/r1.fa"
# a fragment spanning the two bubbles -> case 2 / 6. L is 300 and X1 is 120, so a start at 400 sits
# 20 bp before the end of bubble 1 and a 100 bp insert reaches across M into bubble 2.
mkpair "$H_AD" span_1 400 40 100 "$OUT/r1.fa"

"$BIN" genotype-frag -i "$OUT/g.gfa" -b "$OUT/b" -o "$OUT/u1" -R "$OUT/r1.fa" \
  --origin-universe "$OUT/u1.tsv" "${CF[@]}" -q >/dev/null 2>&1
if [ ! -s "$OUT/u1.tsv" ]; then
  bad "the oracle produced no universe for fixture 1"
else
  read -r NB SC <<<"$(scope_of "$OUT/u1.tsv" uniq_)"
  [ "${NB:-0}" = 1 ] && ok "1 unique local fragment: scope is ONE block ($SC)" \
                     || bad "1 unique local fragment: scope $NB blocks ($SC), expected 1"
  read -r NB2 SC2 <<<"$(scope_of "$OUT/u1.tsv" span_)"
  CONSEC=$($PY -c "
v=[int(x) for x in '${SC2:-}'.split(',') if x!='' and x!='.']
print('yes' if len(v)==2 and v[1]-v[0]==1 else 'no')" 2>/dev/null)
  { [ "${NB2:-0}" = 2 ] && [ "$CONSEC" = yes ]; } \
    && ok "2 boundary-spanning fragment: scope is TWO CONSECUTIVE blocks ($SC2)" \
    || bad "2 boundary-spanning fragment: scope $NB2 ($SC2), expected two consecutive"
  # 6: the span fragment is recruited to one block by the recruiter but depends on two
  [ "${NB2:-0}" -gt 1 ] \
    && ok "6 a fragment whose recruitment label is narrower than its scope is detected ($SC2)" \
    || bad "6 no fragment with scope wider than one block; the blocker's case is unreachable here"
fi

# 3 -- candidate independence. The universe is enumerated over the whole candidate set BEFORE any
# pair is chosen, so two runs must agree exactly. pAD/pCB and pAB/pCD carry identical unordered
# alleles and differ only in pairing: if the scope moved between them, the factor topology would be
# genotype-dependent, which is the thing this design exists to prevent.
"$BIN" genotype-frag -i "$OUT/g.gfa" -b "$OUT/b" -o "$OUT/u3a" -R "$OUT/r1.fa" \
  --origin-universe "$OUT/u3a.tsv" "${CF[@]}" -q >/dev/null 2>&1
"$BIN" genotype-frag -i "$OUT/g.gfa" -b "$OUT/b" -o "$OUT/u3b" -R "$OUT/r1.fa" \
  --origin-universe "$OUT/u3b.tsv" --truth-haplotypes pAD,pCB "${CF[@]}" -q >/dev/null 2>&1
if [ -s "$OUT/u3a.tsv" ] && [ -s "$OUT/u3b.tsv" ]; then
  cmp -s "$OUT/u3a.tsv" "$OUT/u3b.tsv" \
    && ok "3 the origin universe is identical regardless of which pair is named" \
    || bad "3 naming a candidate pair CHANGED the universe: the scope is genotype-dependent"
else
  bad "3 could not produce both universes"
fi

# ---------------------------------------------------------------------------------------------
# 4 -- a tandem repeat INSIDE one block. Multiplicity must raise the origin count without widening
# scope: copies of the same unit in one block are alternatives for the same block variable.
U=$(seq_of 60 41); FL=$(seq_of 300 42); FR=$(seq_of 300 43); ALT=$(seq_of 60 44)
REP4="${U}${U}${U}${U}"
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\t1\t%s\n' "$FL"; printf 'S\t2\t%s\n' "$REP4"; printf 'S\t3\t%s\n' "$ALT"
  printf 'S\t4\t%s\n' "$FR"
  for a in 2 3; do printf 'L\t1\t+\t%s\t+\t0M\nL\t%s\t+\t4\t+\t0M\n' "$a" "$a"; done
  printf 'P\tref\t1+,2+,4+\t*\n'; printf 'P\trep\t1+,2+,4+\t*\n'; printf 'P\talt\t1+,3+,4+\t*\n'
} > "$OUT/g4.gfa"
"$BIN" bubble -i "$OUT/g4.gfa" -r ref -o "$OUT/b4" --min-variant-bp 0 -q >/dev/null 2>&1
: > "$OUT/r4.fa"
mkpair "$U" inrep_1 0 40 60 "$OUT/r4.fa"

"$BIN" genotype-frag -i "$OUT/g4.gfa" -b "$OUT/b4" -o "$OUT/u4" -R "$OUT/r4.fa" \
  --origin-universe "$OUT/u4.tsv" "${CF[@]}" -q >/dev/null 2>&1
if [ -s "$OUT/u4.tsv" ]; then
  read -r NB4 SC4 <<<"$(scope_of "$OUT/u4.tsv" inrep_)"
  NO4=$(orig_of "$OUT/u4.tsv" inrep_)
  [ "${NB4:-0}" = 1 ] && ok "4 repeats within a block: $NO4 origins, scope still ONE block ($SC4)" \
                      || bad "4 repeats within a block widened scope to $NB4 ($SC4), expected 1"
  [ "${NO4:-0}" -gt 1 ] && ok "4 multiplicity is present ($NO4 origins), so the check is not vacuous" \
                        || bad "4 only $NO4 origin: the repeat produced no multiplicity"
else
  bad "4 could not produce the repeat universe"
fi

# ---------------------------------------------------------------------------------------------
# 5 -- the SAME sequence in two NON-adjacent blocks. Scope must be non-consecutive: this is the
# shared/higher-order case that cannot be written as a unary or transition factor at all.
#
# THREE bubbles, shared unit in the FIRST and the THIRD. Two bubbles put it in blocks 1 and 2 --
# consecutive -- so the assertion passed as a boundary case while testing nothing about non-local
# scope. Measured: that is exactly what the first version of this case did.
SH=$(seq_of 80 51); A1=$(seq_of 200 52); A2=$(seq_of 200 53); A3=$(seq_of 200 54); A4=$(seq_of 200 57)
V1=$(seq_of 80 55); V2=$(seq_of 80 56); V3=$(seq_of 80 58); W2=$(seq_of 80 59)
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\t1\t%s\n' "$A1";  printf 'S\t2\t%s\n' "$SH"; printf 'S\t3\t%s\n' "$V1"
  printf 'S\t4\t%s\n' "$A2";  printf 'S\t5\t%s\n' "$W2"; printf 'S\t6\t%s\n' "$V3"
  printf 'S\t7\t%s\n' "$A3";  printf 'S\t8\t%s\n' "$SH"; printf 'S\t9\t%s\n' "$V2"
  printf 'S\t10\t%s\n' "$A4"
  for a in 2 3; do printf 'L\t1\t+\t%s\t+\t0M\nL\t%s\t+\t4\t+\t0M\n' "$a" "$a"; done
  for a in 5 6; do printf 'L\t4\t+\t%s\t+\t0M\nL\t%s\t+\t7\t+\t0M\n' "$a" "$a"; done
  for a in 8 9; do printf 'L\t7\t+\t%s\t+\t0M\nL\t%s\t+\t10\t+\t0M\n' "$a" "$a"; done
  printf 'P\tref\t1+,2+,4+,5+,7+,8+,10+\t*\n'
  printf 'P\tsh\t1+,2+,4+,5+,7+,8+,10+\t*\n'
  printf 'P\tvv\t1+,3+,4+,6+,7+,9+,10+\t*\n'
} > "$OUT/g5.gfa"
"$BIN" bubble -i "$OUT/g5.gfa" -r ref -o "$OUT/b5" --min-variant-bp 0 -q >/dev/null 2>&1
: > "$OUT/r5.fa"
mkpair "$SH" shared_1 0 40 80 "$OUT/r5.fa"
"$BIN" genotype-frag -i "$OUT/g5.gfa" -b "$OUT/b5" -o "$OUT/u5" -R "$OUT/r5.fa" \
  --origin-universe "$OUT/u5.tsv" "${CF[@]}" -q >/dev/null 2>&1
if [ -s "$OUT/u5.tsv" ]; then
  read -r NB5 SC5 <<<"$(scope_of "$OUT/u5.tsv" shared_)"
  NONADJ=$($PY -c "
v=[int(x) for x in '${SC5:-}'.split(',') if x not in ('','.')]
print('yes' if len(v)>=2 and (v[-1]-v[0]+1)!=len(v) else 'no')" 2>/dev/null)
  [ "$NONADJ" = yes ] \
    && ok "5 sequence shared between distant blocks: scope is NON-consecutive ($SC5)" \
    || bad "5 expected a non-consecutive scope, got $NB5 blocks ($SC5)"
else
  bad "5 could not produce the shared-sequence universe"
fi

# ---------------------------------------------------------------------------------------------
# 7 -- unequal copy number. More copies must raise the summed origin MASS: that is what makes
# multiplicity evidence rather than noise. If the two scored identically, dosage would be invisible
# to the exact model too, and no emission built on it could ever recover copy number.
REP2="${U}${U}"; REP6="${U}${U}${U}${U}${U}${U}"
for CN in 2 6; do
  eval "R=\$REP$CN"
  { printf 'H\tVN:Z:1.0\n'
    printf 'S\t1\t%s\n' "$FL"; printf 'S\t2\t%s\n' "$R"; printf 'S\t3\t%s\n' "$ALT"
    printf 'S\t4\t%s\n' "$FR"
    for a in 2 3; do printf 'L\t1\t+\t%s\t+\t0M\nL\t%s\t+\t4\t+\t0M\n' "$a" "$a"; done
    printf 'P\tref\t1+,2+,4+\t*\n'; printf 'P\tcn\t1+,2+,4+\t*\n'; printf 'P\talt\t1+,3+,4+\t*\n'
  } > "$OUT/g7_$CN.gfa"
  "$BIN" bubble -i "$OUT/g7_$CN.gfa" -r ref -o "$OUT/b7_$CN" --min-variant-bp 0 -q >/dev/null 2>&1
  "$BIN" genotype-frag -i "$OUT/g7_$CN.gfa" -b "$OUT/b7_$CN" -o "$OUT/u7_$CN" -R "$OUT/r4.fa" \
    --origin-universe "$OUT/u7_$CN.tsv" "${CF[@]}" -q >/dev/null 2>&1
done
if [ -s "$OUT/u7_2.tsv" ] && [ -s "$OUT/u7_6.tsv" ]; then
  M2=$(lse_of "$OUT/u7_2.tsv" inrep_); M6=$(lse_of "$OUT/u7_6.tsv" inrep_)
  O2=$(orig_of "$OUT/u7_2.tsv" inrep_); O6=$(orig_of "$OUT/u7_6.tsv" inrep_)
  GREW=$($PY -c "print('yes' if float('$M6') > float('$M2') + 1e-9 else 'no')" 2>/dev/null)
  [ "$GREW" = yes ] \
    && ok "7 unequal copy number: 6 copies carry more origin mass than 2 ($M6 > $M2)" \
    || bad "7 copy number did not change origin mass ($M2 vs $M6): dosage is invisible to the oracle"
  [ "${O6:-0}" -gt "${O2:-0}" ] \
    && ok "7 and more origins with more copies ($O2 -> $O6), so the mass change is the repeat" \
    || bad "7 origin count did not grow with copy number ($O2 -> $O6)"
else
  bad "7 could not produce both copy-number universes"
fi

# ---------------------------------------------------------------------------------------------
# 9 -- PANEL COMPOSITION MUST NOT CHANGE SCOPE.
#
# Scope is decided PER CANDIDATE and then unioned. Deciding it on the aggregate over all candidates
# would make it depend on which haplotypes the panel happens to hold: duplicating one path would add
# mass and could flip a block in or out. Duplicating a path is the sharpest form of that -- it adds
# no information whatsoever, so any change in scope is definitionally an artifact.
{ cat "$OUT/g5.gfa"; printf 'P\tsh_dup1\t1+,2+,4+,5+,7+,8+,10+\t*\n'
                     printf 'P\tsh_dup2\t1+,2+,4+,5+,7+,8+,10+\t*\n'
                     printf 'P\tsh_dup3\t1+,2+,4+,5+,7+,8+,10+\t*\n'; } > "$OUT/g5dup.gfa"
"$BIN" bubble -i "$OUT/g5dup.gfa" -r ref -o "$OUT/b5dup" --min-variant-bp 0 -q >/dev/null 2>&1
"$BIN" genotype-frag -i "$OUT/g5dup.gfa" -b "$OUT/b5dup" -o "$OUT/u5dup" -R "$OUT/r5.fa" \
  --origin-universe "$OUT/u5dup.tsv" "${CF[@]}" -q >/dev/null 2>&1
if [ -s "$OUT/u5dup.tsv" ]; then
  read -r NBD SCD <<<"$(scope_of "$OUT/u5dup.tsv" shared_)"
  read -r NB5b SC5b <<<"$(scope_of "$OUT/u5.tsv" shared_)"
  [ "$SCD" = "$SC5b" ] \
    && ok "9 duplicating a panel path leaves the scope unchanged ($SC5b)" \
    || bad "9 duplicating a path changed the scope: $SC5b -> $SCD; scope depends on panel composition"
else
  bad "9 could not produce the duplicated-panel universe"
fi

# 10 -- A BLOCK THAT MATTERS TO ONE LOW-MASS CANDIDATE MUST SURVIVE.
#
# The aggregate rule drowned exactly this case: a candidate carrying little total mass can still have
# a block whose allele changes ITS likelihood, and unioning per-candidate decisions is what keeps it.
# Here `vv` holds none of the shared unit, so it is the low-mass candidate; `sh` holds it twice.
# Block 2 must NOT enter scope through vv (the background absorbs its changes), while blocks 1 and 3
# must remain through sh.
read -r NB10 SC10 <<<"$(scope_of "$OUT/u5.tsv" shared_)"
HAS13=$($PY -c "
v=set(x for x in '${SC10:-}'.split(',') if x not in ('','.'))
print('yes' if {'1','3'} <= v and '2' not in v else 'no')" 2>/dev/null)
[ "$HAS13" = yes ] \
  && ok "10 blocks that matter are kept and a background-dominated one is not ($SC10)" \
  || bad "10 scope $SC10 -- expected the two shared blocks and not the middle one"

# ---------------------------------------------------------------------------------------------
# 11 -- ANTIPARALLEL AND UNPROJECTABLE PATHS USE THE AUTHORITATIVE WALK.
#
# The oracle must take walk bytes and a VERIFIED walk-to-block map, never rebuild a candidate by
# concatenating block alleles: the concatenation is reverse-complemented for an antiparallel path and
# short for an unprojectable one. A candidate with no verifiable map must be skipped, not guessed --
# guessing puts arbitrary blocks into some fragment's dependency set.
{ cat "$OUT/g.gfa"
  printf 'P\trevX\t7-,6-,4-,3-,1-\t*\n'          # antiparallel: block spelling is the RC
  printf 'P\ttruncX\t4+,6+,7+\t*\n'; } > "$OUT/gaf.gfa"   # starts mid-locus: unprojectable
"$BIN" bubble -i "$OUT/gaf.gfa" -r ref -o "$OUT/baf" --min-variant-bp 0 -q >/dev/null 2>&1
AFLOG=$("$BIN" genotype-frag -i "$OUT/gaf.gfa" -b "$OUT/baf" -o "$OUT/uaf" -R "$OUT/r1.fa" \
         --origin-universe "$OUT/uaf.tsv" "${CF[@]}" 2>&1 >/dev/null | grep -i "candidate frames")
NANTI=$(printf '%s' "$AFLOG" | sed -nE 's/.*\(([0-9]+) antiparallel.*/\1/p')
NSKIP=$(printf '%s' "$AFLOG" | sed -nE 's/.*, ([0-9]+) skipped.*/\1/p')
[ "${NANTI:-0}" -ge 1 ] \
  && ok "11 an antiparallel path is mapped by mirroring, not rebuilt ($NANTI)" \
  || bad "11 no antiparallel candidate was recognised; the mirroring path is untested ($AFLOG)"
[ "${NSKIP:-0}" -ge 1 ] \
  && ok "11 a path with no verified walk-to-block map is SKIPPED, not guessed ($NSKIP)" \
  || bad "11 the unprojectable path was not skipped ($AFLOG)"
# and the scope must still be right with those candidates present
if [ -s "$OUT/uaf.tsv" ]; then
  read -r NBA SCA <<<"$(scope_of "$OUT/uaf.tsv" uniq_)"
  [ "${NBA:-0}" = 1 ] \
    && ok "11 scope is unchanged by the presence of antiparallel/unprojectable paths ($SCA)" \
    || bad "11 scope became $NBA blocks ($SCA) once those paths were present"
fi

# ---------------------------------------------------------------------------------------------
# 12 -- ANTIPARALLEL ALONE. Gate 11 puts a reverse path alongside its forward twin, so the FORWARD
# path can preserve the expected scope while the reverse origin is broken and simply never decides
# anything. Here the reverse candidate is the only one carrying the sequence.
#
# The specific defect this catches: for an antiparallel candidate block_of decreases along walk
# coordinates, so an origin's (start,end) come back REVERSED. A reversed interval makes the
# touched-block loop visit nothing and makes the in-scope test accept the origin vacuously -- the
# origin is then exempt from every scope decision, silently.
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\t1\t%s\n' "$L"; printf 'S\t2\t%s\n' "$X1"; printf 'S\t3\t%s\n' "$X2"
  printf 'S\t4\t%s\n' "$M"; printf 'S\t5\t%s\n' "$Y1"; printf 'S\t6\t%s\n' "$Y2"
  printf 'S\t7\t%s\n' "$N"
  for a in 2 3; do printf 'L\t1\t+\t%s\t+\t0M\nL\t%s\t+\t4\t+\t0M\n' "$a" "$a"; done
  for a in 5 6; do printf 'L\t4\t+\t%s\t+\t0M\nL\t%s\t+\t7\t+\t0M\n' "$a" "$a"; done
  printf 'P\tref\t1+,2+,4+,5+,7+\t*\n'
  printf 'P\tfwdOnly\t1+,2+,4+,6+,7+\t*\n'
} > "$OUT/gfwd.gfa"
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\t1\t%s\n' "$L"; printf 'S\t2\t%s\n' "$X1"; printf 'S\t3\t%s\n' "$X2"
  printf 'S\t4\t%s\n' "$M"; printf 'S\t5\t%s\n' "$Y1"; printf 'S\t6\t%s\n' "$Y2"
  printf 'S\t7\t%s\n' "$N"
  for a in 2 3; do printf 'L\t1\t+\t%s\t+\t0M\nL\t%s\t+\t4\t+\t0M\n' "$a" "$a"; done
  for a in 5 6; do printf 'L\t4\t+\t%s\t+\t0M\nL\t%s\t+\t7\t+\t0M\n' "$a" "$a"; done
  printf 'P\tref\t1+,2+,4+,5+,7+\t*\n'
  printf 'P\trevOnly\t7-,6-,4-,2-,1-\t*\n'
} > "$OUT/grev.gfa"
# fwdOnly/revOnly deliberately spell H_AD (1,2,4,6,7) and its reverse complement. The earlier
# version used 1,3,4,6,7, which no fragment here matches: span_ then had ZERO origins, its scope was
# empty on BOTH panels, and "same scope" compared . to . -- a gate that passes whatever the frame
# code does. The whole point of the case is an antiparallel origin that genuinely crosses a block
# boundary, so the fixture has to contain one.
"$BIN" bubble -i "$OUT/gfwd.gfa" -r ref -o "$OUT/bfwd" --min-variant-bp 0 -q >/dev/null 2>&1
"$BIN" bubble -i "$OUT/grev.gfa" -r ref -o "$OUT/brev" --min-variant-bp 0 -q >/dev/null 2>&1
"$BIN" genotype-frag -i "$OUT/gfwd.gfa" -b "$OUT/bfwd" -o "$OUT/ufwd" -R "$OUT/r1.fa" \
  --origin-universe "$OUT/ufwd.tsv" "${CF[@]}" -q >/dev/null 2>&1
"$BIN" genotype-frag -i "$OUT/grev.gfa" -b "$OUT/brev" -o "$OUT/urev" -R "$OUT/r1.fa" \
  --origin-universe "$OUT/urev.tsv" "${CF[@]}" -q >/dev/null 2>&1
if [ -s "$OUT/ufwd.tsv" ] && [ -s "$OUT/urev.tsv" ]; then
  # NOT compared across the two panels by block index. Reversing the walk changes what the bubble
  # decomposition emits -- here bfwd has bubbles (1+,4+){2;3} and (4+,7+){5;6} while brev has
  # (2+,4+){1;3} and (4+,6+){5;7} -- so "same scope" across panels would be comparing indices drawn
  # from different block sets, and neither the indices nor the inside-node sets correspond. The
  # earlier version of this gate passed only because span_ had no origins on either side. The
  # orientation invariant that IS meaningful is same-graph, and case 11 asserts it.
  if diff -q "$OUT/bfwd.bubbles.csv" "$OUT/brev.bubbles.csv" >/dev/null 2>&1; then
    bad "12 the two panels now decompose identically; use the direct index comparison instead"
  else
    ok "12 the reverse panel decomposes into DIFFERENT blocks, so only same-graph scopes compare"
  fi
  # NON-VACUITY. An empty scope on both sides would satisfy the equality above without exercising
  # any frame logic. The spanning fragment must resolve to two consecutive blocks on BOTH panels --
  # on the reverse panel that interval can only come out right if the block span is ordered.
  read -r NRS SRS <<<"$(scope_of "$OUT/urev.tsv" span_)"
  read -r NFS SFS <<<"$(scope_of "$OUT/ufwd.tsv" span_)"
  read -r NRU SRU <<<"$(scope_of "$OUT/urev.tsv" uniq_)"
  [ "${NRU:-0}" = 1 ] \
    && ok "12 a local fragment on the antiparallel panel stays in ONE block ($SRU)" \
    || bad "12 local fragment on the antiparallel panel spans $NRU blocks ($SRU), expected 1"
  CONS12=$($PY -c "
v=[int(x) for x in '${SRS:-}'.split(',') if x!='' and x!='.']
print('yes' if len(v)==2 and v[1]-v[0]==1 else 'no')" 2>/dev/null)
  { [ "${NRS:-0}" = 2 ] && [ "$CONS12" = yes ]; } \
    && ok "12 the ANTIPARALLEL spanning origin crosses two consecutive blocks ($SRS), not nothing" \
    || bad "12 antiparallel span scope is $NRS ($SRS); a vacuous . cannot detect a reversed interval"
  [ "${NFS:-0}" = 2 ] \
    && ok "12 and the forward panel agrees it spans two blocks ($SFS)" \
    || bad "12 forward panel span scope is $NFS ($SFS), expected 2"
  # the reverse origin must also be REJECTED when a required block is removed, i.e. it is genuinely
  # subject to the scope test rather than vacuously exempt
  BH=$(awk -F'\t' 'NR>1 && index($1,"span_")==1 {print $10; exit}' "$OUT/urev.tsv")
  [ "$BH" = yes ] \
    && ok "12 the reverse-only candidate satisfies the joint per-candidate bound" \
    || bad "12 the reverse-only candidate does not satisfy its bound ($BH)"
else
  bad "12 could not produce forward-only and reverse-only universes"
fi

# 13 -- the JOINT bound must be reported and must hold. Testing blocks one at a time bounds nothing
# about removing them together.
NB13=$(awk -F'\t' 'NR>1 && $10=="NO"' "$OUT/u1.tsv" | wc -l | tr -d ' ')
[ "${NB13:-0}" = 0 ] \
  && ok "13 every fragment's scope satisfies the joint per-candidate bound" \
  || bad "13 $NB13 fragments report bound_holds=NO"
WB=$(awk -F'\t' 'NR>1{if($9+0>m) m=$9+0} END{printf "%.9f", m+0}' "$OUT/u1.tsv")
ok "13 worst achieved joint bound across fragments: $WB"

# ---------------------------------------------------------------------------------------------
# THE RECONCILIATION GATE. Scope-restricted scoring must EQUAL the whole-locus reference for every
# candidate -- not merely differ by a constant. A scope that holds all the mass loses nothing, so
# there is no constant left to absorb, and exposure is charged once for the locus rather than per
# factor.
#
# This is the property recruitment scoping FAILED: cropping by recruitment label dropped origins
# whose availability depended on phase, and the whole-minus-factors difference moved with the
# candidate (spread 126.98, isolated to two fragments at 63.4905 each). See
# tests/genotype_frag_factorisation.sh.
RESID=""
while read -r P Q; do
  [ -z "$P" ] && continue
  LINE=$("$BIN" genotype-frag -i "$OUT/g.gfa" -b "$OUT/b" -o "$OUT/rc" -R "$OUT/r1.fa" \
          --reconcile-scope "$P" "$Q" "${CF[@]}" -q 2>/dev/null | grep -E '^-?[0-9]' | tail -1)
  D=$($PY -c "
p='$LINE'.split()
print('%.9f' % abs(float(p[2])) if len(p)==3 else 'NA')" 2>/dev/null)
  RESID="$RESID $D"
done <<'PAIRS'
pAD pCB
pAB pCD
pAB pAB
pAD pAD
pCD pCB
PAIRS
WORST=$($PY -c "
v=[x for x in '$RESID'.split() if x!='NA']
print('NA' if len(v)<4 else '%.9f' % max(map(float,v)))" 2>/dev/null)
case "$WORST" in
  NA) bad "8 could not reconcile scope-restricted against whole-locus for every candidate" ;;
  *) $PY -c "import sys; sys.exit(0 if float('$WORST') < 1e-6 else 1)" 2>/dev/null \
       && ok "8 scope-restricted scoring EQUALS the whole-locus reference for every candidate (worst |residual| $WORST)" \
       || bad "8 scope-restricted scoring differs from the reference (worst |residual| $WORST)" ;;
esac

# and it must be checked on BOTH phases, or it cannot detect a scope blind to linkage
ok "8 the candidate set above includes both phases (pAD/pCB and pAB/pCD, identical unordered alleles)"

# ---------------------------------------------------------------------------------------------
# 14 -- the GREEDY REPAIR path. Case 13 only ever observed initial_bound == 0: every block that
# mattered was already in scope, so the add-back loop never ran and the code that repairs a failed
# joint bound was never executed by any test. This fixture forces it. Six bubbles carry the SAME
# 200 bp repeat, so a fragment inside it places six times with equal mass. Dropping ONE copy costs
# log(6/5) = 0.182 nats, which is under a 0.25 tolerance, so the per-block test excludes all six --
# and their combined removal takes the mass to zero. The bound must fail, the loop must add copies
# back, and it must stop at five, where log(6/5) <= 0.25 holds again.
R14=$(seq_of 200 91)
# Segment ids are the topological numbering 1..19 (u_{i-1}=3i-2, a_i=3i-1, b_i=3i). Names like
# "u0"/"a1" get renumbered by the sorter, the walk-to-block map then fails to verify, and BOTH
# candidates are dropped with "no verified walk-to-block map" -- the universe comes back empty and
# every threshold below is compared against a zero that means "nothing ran".
u_() { echo $((3 * $1 + 1)); }; a_() { echo $((3 * $1 - 1)); }; b_() { echo $((3 * $1)); }
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\t1\t%s\n' "$(seq_of 300 90)"
  for i in 1 2 3 4 5 6; do
    printf 'S\t%s\t%s\n' "$(a_ $i)" "$R14"
    printf 'S\t%s\t%s\n' "$(b_ $i)" "$(seq_of 200 $((100+i)))"
    printf 'S\t%s\t%s\n' "$(u_ $i)" "$(seq_of 300 $((110+i)))"
  done
  for i in 1 2 3 4 5 6; do
    for x in "$(a_ $i)" "$(b_ $i)"; do
      printf 'L\t%s\t+\t%s\t+\t0M\n' "$(u_ $((i-1)))" "$x"
      printf 'L\t%s\t+\t%s\t+\t0M\n' "$x" "$(u_ $i)"
    done
  done
  PR="1+"; PA="1+"
  for i in 1 2 3 4 5 6; do
    PR="$PR,$(a_ $i)+,$(u_ $i)+"; PA="$PA,$(b_ $i)+,$(u_ $i)+"
  done
  printf 'P\tref\t%s\t*\n' "$PR"; printf 'P\taltB\t%s\t*\n' "$PA"
} > "$OUT/g14.gfa"
"$BIN" bubble -i "$OUT/g14.gfa" -r ref -o "$OUT/b14" --min-variant-bp 0 -q >/dev/null 2>&1
H14="$(seq_of 300 90)"; for i in 1 2 3 4 5 6; do H14="${H14}${R14}$(seq_of 300 $((110+i)))"; done
: > "$OUT/r14.fa"
mkpair "$H14" rep14 320 40 100 "$OUT/r14.fa"   # 320..420 lies inside the first copy (300..500)
"$BIN" genotype-frag -i "$OUT/g14.gfa" -b "$OUT/b14" -o "$OUT/u14" -R "$OUT/r14.fa" \
  --origin-universe "$OUT/u14.tsv" --scope-tol 0.25 "${CF[@]}" -q >/dev/null 2>&1
if [ ! -s "$OUT/u14.tsv" ]; then
  bad "14 the oracle produced no universe for the repeat fixture"
else
  NO14=$(awk -F'\t' 'NR>1 && index($1,"rep14")==1 {print $2; exit}' "$OUT/u14.tsv")
  [ "${NO14:-0}" -gt 0 ] 2>/dev/null \
    && ok "14 the repeat fragment has origins ($NO14), so the thresholds below are not vacuous" \
    || bad "14 the repeat fragment has NO origins; see $OUT/u14.tsv.skipped"
  read -r IB AB BA WP NB14 <<<"$(awk -F'\t' 'NR>1 && index($1,"rep14")==1 {
      print $11" "$9" "$12" "$13" "$7; exit}' "$OUT/u14.tsv")"
  $PY -c "import sys; sys.exit(0 if float('${IB:-0}') > 0.25 else 1)" 2>/dev/null \
    && ok "14 the per-block test alone leaves a FAILING joint bound (initial $IB > 0.25)" \
    || bad "14 initial bound is ${IB:-?}, so the repair loop is never entered and stays untested"
  [ "${BA:-0}" -gt 0 ] 2>/dev/null \
    && ok "14 the greedy loop adds blocks back ($BA added)" \
    || bad "14 blocks_added is ${BA:-?}; the add-back code did not run"
  $PY -c "import sys; sys.exit(0 if float('${AB:-1}') <= 0.25 else 1)" 2>/dev/null \
    && ok "14 and it repairs the bound (achieved $AB <= 0.25 over $NB14 blocks)" \
    || bad "14 the repaired bound is still ${AB:-?}, above the 0.25 tolerance"
  # the worst case must be a genotype, not a haplotype: with one candidate carrying all the mass
  # the worst PAIR is its homozygote, which only appears if a == b is enumerated
  [ "${WP:-}" = "0:0" ] || [ "${WP:-}" = "1:1" ] \
    && ok "14 the worst case is a HOMOZYGOTE pair ($WP), so a==b is in the enumeration" \
    || bad "14 worst pair is ${WP:-?}; expected a homozygote for a single-mass panel"
fi

# ---------------------------------------------------------------------------------------------
# 15 -- the instrument must fail loudly rather than emit plausible zeros. Each of these three was a
# way to get a well-formed table out of a run that measured nothing.
# 15a: a panel whose candidates cannot all be mapped yields origins=0, scope=., bound=0 for every
#      fragment -- indistinguishable from "the scope is trivially small" unless the run refuses.
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\tzz1\t%s\n' "$(seq_of 300 71)"; printf 'S\tzz2\t%s\n' "$(seq_of 200 72)"
  printf 'S\tzz3\t%s\n' "$(seq_of 200 73)"; printf 'S\tzz4\t%s\n' "$(seq_of 300 74)"
  printf 'L\tzz1\t+\tzz2\t+\t0M\nL\tzz2\t+\tzz4\t+\t0M\n'
  printf 'L\tzz1\t+\tzz3\t+\t0M\nL\tzz3\t+\tzz4\t+\t0M\n'
  printf 'P\tref\tzz1+,zz2+,zz4+\t*\n'; printf 'P\taltZ\tzz1+,zz3+,zz4+\t*\n'
} > "$OUT/g15.gfa"
"$BIN" bubble -i "$OUT/g15.gfa" -r ref -o "$OUT/b15" --min-variant-bp 0 -q >/dev/null 2>&1
if "$BIN" genotype-frag -i "$OUT/g15.gfa" -b "$OUT/b15" -o "$OUT/u15" -R "$OUT/r14.fa" \
     --origin-universe "$OUT/u15.tsv" "${CF[@]}" -q >/dev/null 2>&1; then
  NU=$(awk -F'\t' 'NR>1' "$OUT/u15.tsv" 2>/dev/null | wc -l | tr -d ' ')
  bad "15a an all-skipped panel produced a table of $NU rows instead of failing"
else
  ok "15a an all-skipped panel is an instrument failure, not a table of zeros"
fi
[ -s "$OUT/u15.tsv.skipped" ] \
  && ok "15a and the sidecar still names the skipped candidates" \
  || bad "15a the .skipped sidecar was not written, so the cause is unrecoverable"

# 15b: --scope-tol is a threshold the joint bound is compared against. A negative value makes the
#      bound pass against something no residual can meet; a NaN makes every comparison false.
for BADTOL in -1 nan -0.5; do
  if "$BIN" genotype-frag -i "$OUT/g14.gfa" -b "$OUT/b14" -o "$OUT/u15b" -R "$OUT/r14.fa" \
       --origin-universe "$OUT/u15b.tsv" --scope-tol "$BADTOL" "${CF[@]}" -q >/dev/null 2>&1; then
    bad "15b --scope-tol $BADTOL was accepted; the joint bound is then compared against nonsense"
  else
    ok "15b --scope-tol $BADTOL is rejected"
  fi
done
"$BIN" genotype-frag -i "$OUT/g14.gfa" -b "$OUT/b14" -o "$OUT/u15z" -R "$OUT/r14.fa" \
  --origin-universe "$OUT/u15z.tsv" --scope-tol 0 "${CF[@]}" -q >/dev/null 2>&1
[ -s "$OUT/u15z.tsv" ] \
  && ok "15b --scope-tol 0 is allowed and means the exact structural scope" \
  || bad "15b --scope-tol 0 was rejected, but zero is a meaningful request"

# 15c: --scope-tol must actually reach --reconcile-scope. It was parsed into a local that only the
#      --origin-universe branch read, so reconciliation silently ran at the 1e-6 default whatever
#      was asked for -- a flag that appears to be honoured and is not.
R15A=$("$BIN" genotype-frag -i "$OUT/g14.gfa" -b "$OUT/b14" -o "$OUT/u15c" -R "$OUT/r14.fa" \
        --reconcile-scope ref altB --scope-tol 1e-6 "${CF[@]}" -q 2>/dev/null \
        | grep -E '^-?[0-9]' | tail -1)
R15B=$("$BIN" genotype-frag -i "$OUT/g14.gfa" -b "$OUT/b14" -o "$OUT/u15d" -R "$OUT/r14.fa" \
        --reconcile-scope ref altB --scope-tol 12 "${CF[@]}" -q 2>/dev/null \
        | grep -E '^-?[0-9]' | tail -1)
if [ -z "$R15A" ] || [ -z "$R15B" ]; then
  bad "15c reconciliation produced no output at either tolerance"
elif [ "$R15A" = "$R15B" ]; then
  bad "15c --scope-tol does not reach --reconcile-scope: 1e-6 and 12 both give $R15A"
else
  ok "15c --scope-tol reaches --reconcile-scope (a 12-nat tolerance changes the residual)"
fi

echo
if [ "$fails" -eq 0 ]; then echo "origin scope: all assertions passed"; else
  echo "origin scope: $fails assertion(s) failed"; fi
exit "$fails"
