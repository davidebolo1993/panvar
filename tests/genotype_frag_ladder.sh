#!/usr/bin/env bash
# genotype_frag_ladder.sh - which approximations are decision-safe?
#
#   genotype_frag_ladder.sh <panvar-binary> <out-dir>
#
# Runs only AFTER the accelerated scorer reproduces the exact reference with every model difference
# removed (tests/genotype_frag_differential.sh). Each rung adds ONE approximation to the previous one,
# in a frozen order, so that a loss can be attributed:
#
#   rung 1  unrestricted syncmer recruitment          (against exhaustive starts)
#   rung 2  + start binning / deduplication
#   rung 3  + anchor-occurrence cap
#   rung 4  + placement top-k
#
# Local alignment is deliberately NOT a rung. It is a different emission model, not an acceleration,
# and needs its own exact fixed-start edit-likelihood oracle before any number comparing it with
# Hamming is an accuracy measurement rather than a model comparison.
#
# Reported at every rung: retained placement probability MASS (not count), the error in the
# certified-versus-competitor score difference, the reference winner's margin, whether the winner
# changes, runtime and peak memory. A rung is decision-safe only while its score error is smaller than
# the margin it has to preserve.
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/lib/experiment.sh"
BIN="${1:?usage: genotype_frag_ladder.sh <panvar> <outdir>}"
OUT="${2:?}"; mkdir -p "$OUT"; OUT="$OUT/run.$$"; rm -rf "$OUT"; mkdir -p "$OUT"
PY="${PYTHON:-python3}"
exp_init "frag-ladder" "$BIN" "$OUT/provenance.txt"

seq_of() { awk -v n="$1" -v s="$2" 'BEGIN{x=s; b="ACGT";
           for(i=0;i<n;i++){x=(1103515245*x+12345)%2147483648; printf "%s", substr(b,(int(x/65536)%4)+1,1)} }'; }
emit_pairs() { awk -v s="$1" -v tag="$2" -v step="$3" -v off="${4:-0}" 'BEGIN{
    rl=60; ins=150; n=length(s);
    for(i=off;i+ins-1<=n;i+=step){
      r1=substr(s,i+1,rl); r2=substr(s,i+ins-rl+1,rl);
      rc=""; for(j=length(r2);j>0;j--){c=substr(r2,j,1);
        rc=rc (c=="A"?"T":c=="C"?"G":c=="G"?"C":"A")}
      printf ">%s_%d/1\n%s\n>%s_%d/2\n%s\n", tag,i,r1,tag,i,rc } }'; }

# A HIGH-COPY panel, because the approximations must actually bind or the ladder measures nothing.
# With a 2-copy repeat, top-k 2 keeps both copies and an occurrence cap of 8 never fires -- measured,
# all four rungs then report identical scores and the ladder is vacuous.
#
# Here a 100 bp segment is repeated 10 times against 9, so a fragment inside the array has ten valid
# placements: the anchor cap of 8 removes its syncmers outright, and top-k 2 keeps two placements of
# ten. Scaled down (100 bp segments, 60 bp mates, 150 bp inserts) so the exhaustive reference is still
# affordable.
L=$(seq_of 200 3); N=$(seq_of 200 5)
SEG=$(seq_of 100 31); SPC=$(seq_of 100 32); SPCALT=$(seq_of 100 33)
rep() { local out="" i; for ((i=0;i<$1;i++)); do out="$out$SEG"; done; printf '%s' "$out"; }
MANY="${L}$(rep 10)${SPC}${N}"      # 10 copies -- the truth
FEWER="${L}$(rep 9)${SPC}${N}"      # 9 copies  -- one unit short
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\t1\t%s\n' "$L"
  for i in $(seq 2 11); do printf 'S\t%d\t%s\n' "$i" "$SEG"; done
  printf 'S\t12\t%s\n' "$SPC"; printf 'S\t13\t%s\n' "$SPCALT"; printf 'S\t14\t%s\n' "$N"
  printf 'L\t1\t+\t2\t+\t0M\n'
  for i in $(seq 2 10); do printf 'L\t%d\t+\t%d\t+\t0M\n' "$i" "$((i+1))"; done
  printf 'L\t11\t+\t12\t+\t0M\n'          # 10 copies -> spacer
  printf 'L\t10\t+\t12\t+\t0M\n'          # 9 copies  -> spacer
  printf 'L\t2\t+\t13\t+\t0M\n'           # ref: one copy then its own spacer
  printf 'L\t12\t+\t14\t+\t0M\nL\t13\t+\t14\t+\t0M\n'
  printf 'P\tref\t1+,2+,13+,14+\t*\n'
  printf 'P\tmany\t1+,2+,3+,4+,5+,6+,7+,8+,9+,10+,11+,12+,14+\t*\n'
  printf 'P\tfewer\t1+,2+,3+,4+,5+,6+,7+,8+,9+,10+,12+,14+\t*\n'
} > "$OUT/g.gfa"
"$BIN" bubble -i "$OUT/g.gfa" -r ref -o "$OUT/bub" --min-variant-bp 0 -q >/dev/null 2>&1
printf '>fewer\n%s\n' "$FEWER" > "$OUT/one.fa"; printf '>many\n%s\n' "$MANY" > "$OUT/two.fa"
{ emit_pairs "$MANY" tA 20; emit_pairs "$MANY" tB 20 10; } > "$OUT/reads.fa"
echo "  panel: 10-copy array against 9-copy, sample homozygous for 10; $(grep -c '/1$' "$OUT/reads.fa") fragments"

P="--haploid-depth 0.05 --fragment-len 150 --fragment-sd 20 --error-rate 0.01"
R2=$("$BIN" genotype-frag --reference-score "$OUT/two.fa" "$OUT/two.fa" -R "$OUT/reads.fa" $P \
       --dump-fragment-mass "$OUT/ref.mass" 2>/dev/null)
R1=$("$BIN" genotype-frag --reference-score "$OUT/one.fa" "$OUT/one.fa" -R "$OUT/reads.fa" $P 2>/dev/null)
MARGIN=$("$PY" -c "print(f'{float('$R2')-float('$R1'):.1f}')")
echo "  reference: many/many $R2, fewer/fewer $R1, margin $MARGIN nats (winner: many)"
echo
printf "  %-34s %9s %9s %9s %8s %8s %6s\n" rung many/many sep err_vs_ref mass% sec winner

rung() {   # label, extra flags
  local lbl="$1"; shift
  local t0 t1 t o sep err mass rss
  t0=$(date +%s)
  /usr/bin/time -l "$BIN" genotype-frag -i "$OUT/g.gfa" -b "$OUT/bub" -o "$OUT/r" -R "$OUT/reads.fa" \
    --haplotype-mode --joint-marginal --hamming-emission --joint-top-pairs 0 --top-pairs 20 \
    --dump-fragment-mass "$OUT/r.mass" $P "$@" -q > /dev/null 2> "$OUT/time.txt" || true
  t1=$(date +%s)
  rss=$(awk '/maximum resident set size/{printf "%.0f", $1/1048576}' "$OUT/time.txt")
  t=$(awk -F'\t' 'NR>1 && $2=="many" && $3=="many"{print $4}' "$OUT/r.hap_pairs.tsv")
  o=$(awk -F'\t' 'NR>1 && $2=="fewer" && $3=="fewer"{print $4}' "$OUT/r.hap_pairs.tsv")
  local win; win=$(sed -n 2p "$OUT/r.hap_pairs.tsv" | awk -F'\t' '{print ($2==$3)?$2:$2"/"$3}')
  "$PY" - "$R2" "$t" "$o" "$MARGIN" "$OUT/ref.mass" "$OUT/r.mass" "$lbl" "$((t1-t0))" "$win" "${rss:-0}" <<'PYEOF'
import sys, math
ref2, t, o, margin, refmass, fastmass, lbl, secs, win, rss = sys.argv[1:11]
ref2, t, o, margin = float(ref2), float(t or 0), float(o or 0), float(margin)
def load(p):
    d={}
    for l in open(p):
        if l.startswith("#") or l.startswith("fragment"): continue
        k,v=l.rstrip("\n").split("\t")
        d[k]=float(v)
    return d
rm, fm = load(refmass), load(fastmass)
common=[k for k in rm if k in fm and rm[k] > -1e300]
# retained placement MASS, weighted by the reference's own mass so that dominant placements dominate
num=sum(math.exp(fm[k]-rm[k]) * math.exp(rm[k]-max(rm.values())) for k in common)
den=sum(math.exp(rm[k]-max(rm.values())) for k in common)
retained = 100.0 * num / den if den else 0.0
sep = t - o
err = abs(t - ref2)
safe = "yes" if err < margin else "NO"
print(f"  {lbl:<34} {t:9.1f} {sep:9.1f} {err:9.2f} {retained:7.1f}% {secs:>7}s {win:>6}"
      f"   decision-safe: {safe}, peak {rss} MB")
PYEOF
}

rung "1 unrestricted recruitment"        --rung-zero
rung "2 + start binning (64)"            --rung-zero --placement-bin 64
rung "3 + anchor cap (8)"                --rung-zero --placement-bin 64 --max-anchor-occ 8
rung "4 + placement top-k (2)"           --rung-zero --placement-bin 64 --max-anchor-occ 8 --placement-topk 2
echo
echo "  a rung is decision-safe only while its score error is smaller than the margin it must preserve"
echo
echo "  CAVEAT on this fixture: the reference's own margin between 10 copies and 9 is only ~0.6 nats,"
echo "  because placement multiplicity favours the longer array by n*log(10/9) = 14.3 nats while its"
echo "  extra exposure penalises it by lambda*2*100 = 10.0. The two nearly cancel, so a one-unit copy"
echo "  difference is near-degenerate at this depth FOR THE EXACT MODEL TOO. Every approximation"
echo "  therefore flips the winner, and the decision-safety column carries no information here --"
echo "  the retained-mass column is what to read. A fixture with a larger copy difference is needed"
echo "  before decision-safety at an array can be assessed at all."
