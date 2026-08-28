#!/usr/bin/env bash
# genotype_frag_differential.sh - accelerated scorer against the exact reference.
#
#   genotype_frag_differential.sh <panvar-binary> <out-dir>
#
# The reference implements the contract; the fast path is supposed to approximate it. This measures
# by how much, on a small synthetic panel where the reference is affordable.
#
# IT DOES NOT TUNE ANYTHING. The fast path is not adjusted to make the numbers agree -- the reference
# defines the model, and the fast path either approximates it within a stated bound or it does not.
# The output is that bound.
#
# What is compared, and why not the absolute score: the two carry different additive constants, so
# only DIFFERENCES between pairs are meaningful. A caller's behaviour is determined by which pair wins
# and by the margins, so the comparison is (a) top-pair agreement and (b) the error in pairwise score
# differences relative to the reference's own spread.
set -uo pipefail
BIN="${1:?usage: genotype_frag_differential.sh <panvar> <outdir>}"
OUT="${2:?}"; mkdir -p "$OUT"; OUT="$OUT/run.$$"; rm -rf "$OUT"; mkdir -p "$OUT"
PY="${PYTHON:-python3}"
fails=0
ok()  { printf "  ok   %s\n" "$1"; }
bad() { printf "  FAIL %s\n" "$1"; fails=$((fails+1)); }
gt() { "$PY" -c "import sys; a,b,t=float(sys.argv[1]),float(sys.argv[2]),float(sys.argv[3]); sys.exit(0 if a-b>t else 1)" "$1" "$2" "${3:-0.5}"; }

seq_of() { awk -v n="$1" -v s="$2" 'BEGIN{x=s; b="ACGT";
           for(i=0;i<n;i++){x=(1103515245*x+12345)%2147483648; printf "%s", substr(b,(int(x/65536)%4)+1,1)} }'; }

# A panel of six haplotypes over two bubbles: L + X{1,2,3} + M + Y{1,2} + N. Small enough that the
# reference can enumerate every start, varied enough that pair ranking is non-trivial.
L=$(seq_of 500 3); M=$(seq_of 260 4); N=$(seq_of 500 5)
X1=$(seq_of 300 11); X2=$(seq_of 300 12); X3=$(seq_of 300 13); X4=$(seq_of 300 14)
Y1=$(seq_of 300 21); Y2=$(seq_of 300 22)
{ printf 'H\tVN:Z:1.0\n'
  # Nodes MUST be numbered in reference order: the bubble caller renumbers into reference order and a
  # node that sits early in the path but late in the numbering produces a decomposition that spells
  # haplotypes wrongly -- silently. The first version of this fixture put the extra allele after the
  # sink and spelled one haplotype to 0 bp, which then looked like a scorer disagreement.
  printf 'S\t1\t%s\n' "$L"
  printf 'S\t2\t%s\n' "$X1"; printf 'S\t3\t%s\n' "$X2"
  printf 'S\t4\t%s\n' "$X3"; printf 'S\t5\t%s\n' "$X4"
  printf 'S\t6\t%s\n' "$M"
  printf 'S\t7\t%s\n' "$Y1"; printf 'S\t8\t%s\n' "$Y2"
  printf 'S\t9\t%s\n' "$N"
  for a in 2 3 4 5; do printf 'L\t1\t+\t%s\t+\t0M\nL\t%s\t+\t6\t+\t0M\n' "$a" "$a"; done
  for a in 7 8; do printf 'L\t6\t+\t%s\t+\t0M\nL\t%s\t+\t9\t+\t0M\n' "$a" "$a"; done
  printf 'P\tref\t1+,5+,6+,7+,9+\t*\n'
  printf 'P\thA\t1+,2+,6+,7+,9+\t*\n'; printf 'P\thB\t1+,3+,6+,8+,9+\t*\n'
  printf 'P\thC\t1+,4+,6+,7+,9+\t*\n'; printf 'P\thD\t1+,2+,6+,8+,9+\t*\n'
  printf 'P\thE\t1+,3+,6+,7+,9+\t*\n'; printf 'P\thF\t1+,4+,6+,8+,9+\t*\n'
} > "$OUT/g.gfa"
"$BIN" bubble -i "$OUT/g.gfa" -r ref -o "$OUT/bub" --min-variant-bp 0 -q >/dev/null 2>&1

# bash 3.2 (macOS default) has no associative arrays, so these are plain variables.
S_hA="${L}${X1}${M}${Y1}${N}"; S_hB="${L}${X2}${M}${Y2}${N}"
S_hC="${L}${X3}${M}${Y1}${N}"; S_hD="${L}${X1}${M}${Y2}${N}"
S_hE="${L}${X2}${M}${Y1}${N}"; S_hF="${L}${X3}${M}${Y2}${N}"
# ref is a CANDIDATE like any other -- path 1+,5+,6+,7+,9+ -- and the shortlist scores it. It had no
# FASTA, so the reference loop produced seven empty scores that silently dropped out of the
# comparison and reappeared as seven "extra" accelerated pairs.
S_ref="${L}${X4}${M}${Y1}${N}"
seq_for() { case "$1" in hA) printf '%s' "$S_hA";; hB) printf '%s' "$S_hB";; hC) printf '%s' "$S_hC";;
                        hD) printf '%s' "$S_hD";; hE) printf '%s' "$S_hE";; hF) printf '%s' "$S_hF";;
                        ref) printf '%s' "$S_ref";; esac; }
for h in hA hB hC hD hE hF ref; do printf '>%s\n%s\n' "$h" "$(seq_for "$h")" > "$OUT/$h.fa"; done

# a diploid sample: hA / hB
emit_pairs() { awk -v s="$1" -v tag="$2" -v step="$3" -v off="${4:-0}" 'BEGIN{
    rl=120; ins=350; n=length(s);
    for(i=off;i+ins-1<=n;i+=step){
      r1=substr(s,i+1,rl); r2=substr(s,i+ins-rl+1,rl);
      rc=""; for(j=length(r2);j>0;j--){c=substr(r2,j,1);
        rc=rc (c=="A"?"T":c=="C"?"G":c=="G"?"C":"A")}
      printf ">%s_%d/1\n%s\n>%s_%d/2\n%s\n", tag,i,r1,tag,i,rc } }'; }
{ emit_pairs "$S_hA" fA 20; emit_pairs "$S_hB" fB 20; } > "$OUT/reads.fa"
NF=$(grep -c '/1$' "$OUT/reads.fa")
echo "  panel 6 haplotypes + ref, sample hA/hB, $NF fragments"

# FIXTURE SELF-CHECK. Every haplotype must spell to the length the fixture intends, and every one must
# reach the shortlist. A fixture that decomposes wrongly produces a disagreement that looks like a
# scorer defect -- which is what the first version of this test did.
"$BIN" genotype-frag -i "$OUT/g.gfa" -b "$OUT/bub" -o "$OUT/chk" -R "$OUT/reads.fa" \
  --haplotype-mode --top-pairs 1 -q >/dev/null 2>&1
WANT=$(( ${#L} + 300 + ${#M} + 300 + ${#N} ))
badbp=0
while IFS=$'\t' read -r nm bp rest; do
  [ "$nm" = "haplotype" ] && continue
  [ "$bp" = "$WANT" ] || { printf "  FIXTURE %s spells %s bp, expected %s\n" "$nm" "$bp" "$WANT"; badbp=1; }
done < "$OUT/chk.hap_scores.tsv"
NSL=$(( $(wc -l < "$OUT/chk.hap_scores.tsv") - 1 ))
if [ "$badbp" -ne 0 ] || [ "$NSL" -lt 7 ]; then
  bad "fixture is malformed: $NSL of 7 haplotypes shortlisted, spelled-length check $([ $badbp -eq 0 ] && echo ok || echo FAILED)"
  echo; echo "differential: fixture invalid, no comparison attempted"; exit 1
fi
ok "fixture self-check: all 7 haplotypes spell $WANT bp and reach the shortlist"

P="--haploid-depth 0.05 --fragment-len 350 --fragment-sd 50 --error-rate 0.01"

# reference: every unordered pair, over the haplotypes THE SHORTLIST ACTUALLY CONTAINS rather than a
# hardcoded list. Hardcoded, the enumeration silently stopped covering the shortlist the moment a
# scoring change altered which haplotypes reach it: the key-set assertion then reported 7 "extra"
# accelerated pairs that were really 7 reference pairs never computed. Deriving the list keeps the
# two sides comparing the same candidates by construction.
HAPS=$(awk -F'\t' 'NR>1{print $1}' "$OUT/chk.hap_scores.tsv" | sort -u)
: > "$OUT/ref.tsv"
for i in $HAPS; do
  for j in $HAPS; do
    [[ "$i" > "$j" ]] && continue
    v=$("$BIN" genotype-frag --reference-score "$OUT/$i.fa" "$OUT/$j.fa" -R "$OUT/reads.fa" $P 2>/dev/null)
    printf '%s/%s\t%s\n' "$i" "$j" "$v" >> "$OUT/ref.tsv"
  done
done

# accelerated: one run, all pairs, marginal form (the closest to the contract)
"$BIN" genotype-frag -i "$OUT/g.gfa" -b "$OUT/bub" -o "$OUT/fast" -R "$OUT/reads.fa" \
  --haplotype-mode --joint-marginal --joint-top-pairs 0 --placement-topk 64 \
  --top-pairs 100 $P -q >/dev/null 2>&1
awk -F'\t' 'NR>1{
    a=$2; b=$3; if(a>b){t=a;a=b;b=t}
    print a"/"b"\t"$4 }' "$OUT/fast.hap_pairs.tsv" > "$OUT/fast.tsv"

"$PY" - "$OUT/ref.tsv" "$OUT/fast.tsv" <<'PYEOF'
import sys
def load(p):
    d={}
    for l in open(p):
        k,v=l.rstrip("\n").split("\t")
        try: d[k]=float(v)
        except ValueError: pass
    return d
ref, fast = load(sys.argv[1]), load(sys.argv[2])
common=[k for k in ref if k in fast]
if len(common) < 3:
    print(f"  .... only {len(common)} pairs in common; cannot compare"); sys.exit(1)
missing = [k for k in ref if k not in fast]
if missing:
    print(f"  .... {len(missing)} reference pair(s) absent from the accelerated output: {missing[:4]}")
    print("  .... refusing to compare a subset -- a dropped pair is a finding, not a filter")
    sys.exit(1)
rb=max(ref,key=ref.get); fb=max((k for k in common),key=lambda k: fast[k])
print(f"  .... reference best {rb}   accelerated best {fb}   ({len(common)} pairs compared)")
# differences relative to each scorer's own best, which removes the additive constant
rd={k: ref[k]-ref[rb] for k in common}
fd={k: fast[k]-fast[fb] for k in common}
err=[abs(rd[k]-fd[k]) for k in common]
spread=max(rd.values())-min(rd.values())
print(f"  .... reference spread {spread:.1f} nats; max error in pairwise differences {max(err):.1f}"
      f" ({100*max(err)/spread:.1f}% of spread)")
# A percentage of the TOTAL spread can hide a large error next to the winner, which is the only place
# it can change a call. What matters is the error among the top pairs against the winner's margin:
# the acceleration is decision-safe only while the former is smaller than the latter.
order_ref_all=sorted(common,key=lambda k:-ref[k])
margin = ref[order_ref_all[0]] - ref[order_ref_all[1]] if len(order_ref_all) > 1 else float("inf")
top = order_ref_all[:5]
top_err = max(abs(rd[k]-fd[k]) for k in top)
print(f"  .... winner's margin {margin:.1f} nats; max error among the top 5 pairs {top_err:.1f}"
      f"  -> {'DECISION-SAFE' if top_err < margin else 'NOT decision-safe'}")
order_ref=sorted(common,key=lambda k:-ref[k])
order_fast=sorted(common,key=lambda k:-fast[k])
top3=order_ref[:3]==order_fast[:3]
print(f"  .... reference order {order_ref[:3]}")
print(f"  .... accelerated     {order_fast[:3]}")
ok = (rb==fb)
sys.exit(0 if ok else 1)
PYEOF
[ $? -eq 0 ] && ok "accelerated scorer picks the same pair as the exact reference" \
             || bad "accelerated scorer picks a DIFFERENT pair from the exact reference"

# Per-pair agreement, not merely the winner: with the models identical the two scorers should differ
# by ONE candidate-independent constant. A constant that varies by pair is a model difference hiding
# behind a correct ranking.
"$BIN" genotype-frag -i "$OUT/g.gfa" -b "$OUT/bub" -o "$OUT/rz0" -R "$OUT/reads.fa" \
  --haplotype-mode --rung-zero --hamming-emission --top-pairs 100 $P -q >/dev/null 2>&1
awk -F'\t' 'NR>1{a=$2;b=$3; if(a>b){t=a;a=b;b=t} print a"/"b"\t"$4}' "$OUT/rz0.hap_pairs.tsv" > "$OUT/rz0.tsv"
"$PY" - "$OUT/ref.tsv" "$OUT/rz0.tsv" <<'PYEOF'
import sys, statistics
def load(p):
    d={}
    for l in open(p):
        k,v=l.rstrip("\n").split("\t")
        try: d[k]=float(v)
        except ValueError: pass
    return d
ref, rz = load(sys.argv[1]), load(sys.argv[2])
# EXACT key-set equality. Accepting "any common subset of at least three" would pass an accelerated
# scorer that silently dropped most pairs -- including the truth -- which is precisely the failure
# this suite already had once, when a duplicated `ref` sequence removed the truth pair from the
# comparison and the remaining pairs agreed nicely.
missing = sorted(set(ref) - set(rz))
extra   = sorted(set(rz) - set(ref))
if missing or extra:
    print(f"  .... key sets differ: {len(missing)} missing {missing[:4]}, {len(extra)} extra {extra[:4]}")
    sys.exit(1)
common=sorted(ref)
off=[rz[k]-ref[k] for k in common]
c=statistics.median(off)
resid=max(abs(o-c) for o in off)
print(f"  .... per-pair offset: median {c:.2f} nats, max deviation from it {resid:.2f} over "
      f"{len(common)} pairs")
sys.exit(0 if resid < 1.0 else 1)
PYEOF
[ $? -eq 0 ] && ok "every pair agrees with the reference to within 1 nat of one shared constant" \
             || bad "the offset between the scorers varies by pair -- a model difference remains"

# ============================================================== the case truncation actually bites
# A panel where one haplotype carries a segment TWICE. Every fragment from that segment has two valid
# placements, which is exactly the mass placement_topk discards, so this is where the accelerated
# scorer should diverge from the contract if it is going to. Reported as a bound at several topk
# values rather than asserted at one -- the point is to measure the approximation, not to pick a knob.
echo
echo "  repeat panel: one haplotype carries a 300 bp segment twice"
SEG=$(seq_of 300 31); SPC=$(seq_of 200 32); SPCALT=$(seq_of 200 33)
ONE="${L}${SEG}${SPC}${N}"
TWO="${L}${SEG}${SEG}${SPC}${N}"
# ref carries its own spacer. Sharing a sequence with `one` would deduplicate the two out of the
# shortlist and the comparison would then silently run without the pair it is about -- which is what
# the first version of this block did.
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\t1\t%s\n' "$L"; printf 'S\t2\t%s\n' "$SEG"; printf 'S\t3\t%s\n' "$SEG"
  printf 'S\t4\t%s\n' "$SPC"; printf 'S\t5\t%s\n' "$SPCALT"; printf 'S\t6\t%s\n' "$N"
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t3\t+\t4\t+\t0M\n'
  printf 'L\t2\t+\t4\t+\t0M\nL\t2\t+\t5\t+\t0M\n'
  printf 'L\t4\t+\t6\t+\t0M\nL\t5\t+\t6\t+\t0M\n'
  printf 'P\tref\t1+,2+,5+,6+\t*\n'
  printf 'P\tone\t1+,2+,4+,6+\t*\n'; printf 'P\ttwo\t1+,2+,3+,4+,6+\t*\n'
} > "$OUT/r.gfa"
"$BIN" bubble -i "$OUT/r.gfa" -r ref -o "$OUT/rbub" --min-variant-bp 0 -q >/dev/null 2>&1
printf '>one\n%s\n' "$ONE" > "$OUT/one.fa"; printf '>two\n%s\n' "$TWO" > "$OUT/two.fa"
# TWO INDEPENDENT STREAMS. The sample is homozygous for TWO, so at lambda per HOMOLOGUE it carries
# two homologues' worth of fragments. Emitting one stream gives the sample half the coverage both
# candidate models assume, which is itself a factor of two and would be mistaken for one in the
# scorer. Offset starts so the two streams are not identical fragments.
{ emit_pairs "$TWO" tA 20; emit_pairs "$TWO" tB 20 10; } > "$OUT/rreads.fa"
# The observed fragment count must match what lambda and the pair's exposure predict, or the fixture
# and the model disagree before any scorer is compared.
NRF=$(grep -c '/1$' "$OUT/rreads.fa")
"$PY" - "$NRF" "${#TWO}" <<'PYEOF'
import sys, math
n, N = int(sys.argv[1]), int(sys.argv[2])
mu, sd, disc, lam = 350.0, 50.0, 0.01, 0.05
lo, hi = int(mu-4*sd), int(mu+4*sd); span = hi-lo+1
w = [ (1-disc)*math.exp(-0.5*((L-mu)/sd)**2)/(sd*math.sqrt(2*math.pi)) + disc/span
      for L in range(lo, hi+1) ]
tot = sum(w)
E = sum((wi/tot)*max(0, N-L+1) for wi, L in zip(w, range(lo, hi+1)))
expected = lam * 2 * E                       # homozygous: two homologues
print(f"  .... {n} fragments observed; model expects lambda*E(two,two) = {expected:.0f}")
sys.exit(0 if abs(n-expected) < 0.35*expected else 1)
PYEOF
[ $? -eq 0 ] && ok "observed fragment count agrees with lambda x exposure for the homozygous pair" \
             || bad "fragment count disagrees with the model's expectation; the fixture is miscalibrated"

R2=$("$BIN" genotype-frag --reference-score "$OUT/two.fa" "$OUT/two.fa" -R "$OUT/rreads.fa" $P 2>/dev/null)
R1=$("$BIN" genotype-frag --reference-score "$OUT/one.fa" "$OUT/one.fa" -R "$OUT/rreads.fa" $P 2>/dev/null)
if gt "$R2" "$R1"; then
  ok "reference prefers the 2-copy haplotype on 2-copy reads ($R2 > $R1)"
else
  bad "reference does not prefer the 2-copy haplotype: $R2 vs $R1"
fi
REFDIFF=$("$PY" -c "import sys;print(f'{float(sys.argv[1])-float(sys.argv[2]):.1f}')" "$R2" "$R1")
# UNBOUNDED RECRUITED PLACEMENTS first -- deliberately not called "no approximations", because it is
# not. It removes the anchor cap, the start binning and top-k, but it still reaches placements through
# syncmer recruitment and banded local alignment rather than enumerating every (start, L, strand)
# state as the reference does. Those two remain, and until this matches the reference exactly they are
# candidate explanations for any difference. If it does not match, the rungs below it are
# uninterpretable.
# With --hamming-emission the last known model difference is removed too: the reference scores a
# placement at a FIXED position, and banded local alignment can slide a read to a better offset,
# creating placements the reference rejects. Measured on this fixture: 96% without it, 100% with.
"$BIN" genotype-frag -i "$OUT/r.gfa" -b "$OUT/rbub" -o "$OUT/rz" -R "$OUT/rreads.fa" \
  --haplotype-mode --rung-zero --hamming-emission --top-pairs 20 $P -q >/dev/null 2>&1
rzt=$(awk -F'\t' 'NR>1 && $2=="two" && $3=="two"{print $4}' "$OUT/rz.hap_pairs.tsv")
rzo=$(awk -F'\t' 'NR>1 && $2=="one" && $3=="one"{print $4}' "$OUT/rz.hap_pairs.tsv")
if [ -n "$rzt" ] && [ -n "$rzo" ]; then
  "$PY" - "$R2" "$R1" "$rzt" "$rzo" <<'PYEOF'
import sys
r2, r1, t, o = (float(x) for x in sys.argv[1:5])
ref_sep, rz_sep = r2 - r1, t - o
print(f"  .... unbounded recruited placements: two/two {t:.1f} (reference {r2:.1f}), separation "
      f"{rz_sep:.1f} of {ref_sep:.1f} = {100*rz_sep/ref_sep:.0f}%")
# TWO-SIDED. A one-sided ">= 0.9x" passes a result at 150%, which is just as much a disagreement --
# an accelerated scorer that over-separates has invented evidence.
sys.exit(0 if abs(rz_sep / ref_sep - 1.0) <= 0.10 else 1)
PYEOF
  [ $? -eq 0 ] && ok "unbounded recruited placements reproduce the reference separation within +/-10%" \
               || bad "unbounded recruited placements do NOT match the reference -- the models still differ"

  # ABSOLUTE agreement within a PREDECLARED 1-nat tolerance -- not "floating-point equality". The
  # observed residual is about 0.3 nats, which is orders of magnitude above numerical rounding and is
  # presumably the remaining recruitment approximation; calling it floating point would hide the one
  # quantity the ladder is about to measure.
  #
  # It has to be an absolute check, not a relative one: only an absolute comparison catches a term
  # that shifts every candidate by the same amount. Mutation-tested -- deleting the 1/2-per-strand
  # factor passes every ranking-based and offset-based assertion in this suite and is caught only
  # here, shifting the score by 146*log(2) = 101.2 nats.
  "$PY" - "$R2" "$rzt" <<'PYEOF'
import sys
ref, fast = float(sys.argv[1]), float(sys.argv[2])
print(f"  .... absolute agreement on two/two: reference {ref:.3f}, accelerated {fast:.3f}, "
      f"difference {abs(ref-fast):.3f}")
sys.exit(0 if abs(ref - fast) < 1.0 else 1)
PYEOF
  [ $? -eq 0 ] && ok "accelerated score agrees with the reference absolutely, within the 1-nat tolerance" \
               || bad "absolute disagreement beyond 1 nat -- some per-fragment term is not shared"
else
  bad "rung zero did not score both pairs"
fi

printf "  .... reference separates them by %s nats; then the approximations, one knob at a time:\n" "$REFDIFF"
for tk in 2 8 32; do
  "$BIN" genotype-frag -i "$OUT/r.gfa" -b "$OUT/rbub" -o "$OUT/rk$tk" -R "$OUT/rreads.fa" \
    --haplotype-mode --joint-marginal --joint-top-pairs 0 --placement-topk $tk \
    --top-pairs 20 $P -q >/dev/null 2>&1
  t=$(awk -F'\t' 'NR>1 && $2=="two" && $3=="two"{print $4}' "$OUT/rk$tk.hap_pairs.tsv")
  o=$(awk -F'\t' 'NR>1 && $2=="one" && $3=="one"{print $4}' "$OUT/rk$tk.hap_pairs.tsv")
  if [ -z "$t" ] || [ -z "$o" ]; then
    printf "  .... topk %-3s FIXTURE: one/one or two/two absent from the shortlist output\n" "$tk"
    continue
  fi
  if [ -n "$t" ] && [ -n "$o" ]; then
    d=$("$PY" -c "import sys;print(f'{float(sys.argv[1])-float(sys.argv[2]):.1f}')" "$t" "$o")
    rel=$("$PY" -c "import sys;a,b=float(sys.argv[1]),float(sys.argv[2]);print(f'{100*abs(a-b)/abs(b):.0f}')" "$d" "$REFDIFF")
    printf "  .... topk %-3s separation %-10s (%s%% of the reference's)\n" "$tk" "$d" "$rel"
  else
    printf "  .... topk %-3s one of the pairs was not scored\n" "$tk"
  fi
done

echo
if [ "$fails" -eq 0 ]; then echo "differential: agreement on the winning pair"; else
  echo "differential: $fails disagreement(s) -- the accelerated path does not implement the contract"; fi
exit "$fails"
