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
# FASTA basenames must match the PANEL path names, or a named-pair mass dump cannot be requested and
# the reference winner cannot be looked up in the accelerated output.
printf '>fewer\n%s\n' "$FEWER" > "$OUT/fewer.fa"; printf '>many\n%s\n' "$MANY" > "$OUT/many.fa"
{ emit_pairs "$MANY" tA 20; emit_pairs "$MANY" tB 20 10; } > "$OUT/reads.fa"
echo "  panel: 10-copy array against 9-copy, sample homozygous for 10; $(grep -c '/1$' "$OUT/reads.fa") fragments"

P="--haploid-depth 0.05 --fragment-len 150 --fragment-sd 20 --error-rate 0.01"

# EVERY diplotype, so the margin is the reference winner's own margin. Scoring only many/many against
# fewer/fewer gives a two-copy diploid contrast, not the nearest alternative, and cannot say what the
# exact model's global winner is -- which is the quantity "the approximation flipped the winner" is
# about.
printf '>ref\n%s\n' "${L}${SEG}${SPCALT}${N}" > "$OUT/ref.fa"
: > "$OUT/refscores.tsv"
for i in ref fewer many; do
  for j in ref fewer many; do
    [[ "$i" > "$j" ]] && continue
    v=$("$BIN" genotype-frag --reference-score "$OUT/$i.fa" "$OUT/$j.fa" -R "$OUT/reads.fa" $P 2>/dev/null)
    printf '%s/%s\t%s\n' "$i" "$j" "$v" >> "$OUT/refscores.tsv"
  done
done
read -r REFWIN R2 MARGIN < <("$PY" - "$OUT/refscores.tsv" <<'PYEOF'
import sys
rows=[l.split("\t") for l in open(sys.argv[1]) if l.strip()]
sc=sorted(((float(v), k) for k,v in rows), reverse=True)
print(sc[0][1], f"{sc[0][0]:.6f}", f"{sc[0][0]-sc[1][0]:.2f}")
PYEOF
)
echo "  reference winner: $REFWIN at $R2; margin over the runner-up: $MARGIN nats"
sed 's/^/    /' "$OUT/refscores.tsv"
echo
# reference mass for the winner, so the accelerated dump is compared with the SAME pair
RW1=${REFWIN%%/*}; RW2=${REFWIN##*/}
"$BIN" genotype-frag --reference-score "$OUT/$RW1.fa" "$OUT/$RW2.fa" -R "$OUT/reads.fa" $P \
  --dump-fragment-mass "$OUT/ref.mass" >/dev/null 2>&1
# MANDATORY FIXTURE SELF-CHECK. The accelerated path scores the panel's SPELLED haplotypes while the
# reference scores the FASTA files, and nothing else in this harness notices if they differ. Measured:
# a 10-copy tandem array built as a chain of identical nodes spelled `many` at 1100 bp against a
# 1500 bp FASTA and `fewer` at 0 bp, and the ladder still produced a full table of plausible
# percentages. Refuse to run rather than report numbers from two different sequences.
"$BIN" genotype-frag -i "$OUT/g.gfa" -b "$OUT/bub" -o "$OUT/chk" -R "$OUT/reads.fa" \
  --haplotype-mode --top-pairs 1 $P -q >/dev/null 2>&1
badfix=0
while IFS=$'\t' read -r nm bp rest; do
  [ "$nm" = "haplotype" ] && continue
  want=$(awk '!/^>/{n+=length($0)} END{print n+0}' "$OUT/$nm.fa" 2>/dev/null)
  if [ -z "$want" ] || [ "$want" = "0" ] || [ "$bp" != "$want" ]; then
    printf "  FIXTURE MISMATCH: panel spells %s as %s bp, reference FASTA is %s bp\n" "$nm" "$bp" "${want:-missing}"
    badfix=1
  fi
done < "$OUT/chk.hap_scores.tsv"
if [ "$badfix" -ne 0 ]; then
  echo
  echo "  The panel and the reference are not scoring the same sequences. No ladder number from this"
  echo "  fixture means anything, so none is produced. A tandem array built as a chain of identical"
  echo "  nodes does not decompose into the alleles this harness assumes; the fixture needs rebuilding"
  echo "  before the ladder can run on a high-copy case."
  exit 1
fi
ok_fixture="all panel haplotypes spell exactly what the reference scores"
printf "  ok   %s\n\n" "$ok_fixture"
printf "  %-38s %10s %10s %8s %7s %6s\n" rung "winner_sc" "err_vs_ref" "mass%" sec winner

rung() {   # label, tag, extra flags
  local lbl="$1" tag="$2"; shift 2
  local pfx="$OUT/rung_$tag"
  rm -f "$pfx".* 2>/dev/null || true
  local t0 t1 rc=0
  t0=$(date +%s)
  /usr/bin/time -l "$BIN" genotype-frag -i "$OUT/g.gfa" -b "$OUT/bub" -o "$pfx" -R "$OUT/reads.fa" \
    --haplotype-mode --joint-marginal --hamming-emission --joint-top-pairs 0 --top-pairs 20 \
    --dump-fragment-mass "$pfx.mass" --dump-mass-pair "$RW1,$RW2" $P "$@" -q \
    > /dev/null 2> "$pfx.time" || rc=$?
  t1=$(date +%s)
  # No `|| true`: a failed rung must not be reported using the previous rung's files.
  if [ $rc -ne 0 ] || [ ! -s "$pfx.hap_pairs.tsv" ] || [ ! -s "$pfx.mass" ]; then
    printf "  %-38s FAILED (rc=%d) -- not reported\n" "$lbl" "$rc"; return 1
  fi
  local rss; rss=$(awk '/maximum resident set size/{printf "%.0f", $1/1048576}' "$pfx.time")
  local sc; sc=$(awk -F'\t' -v a="$RW1" -v b="$RW2" 'NR>1 && (($2==a&&$3==b)||($2==b&&$3==a)){print $4}' "$pfx.hap_pairs.tsv")
  local win; win=$(sed -n 2p "$pfx.hap_pairs.tsv" | awk -F'\t' '{print ($2==$3)?$2:$2"/"$3}')
  "$PY" - "$R2" "${sc:-nan}" "$MARGIN" "$OUT/ref.mass" "$pfx.mass" "$lbl" "$((t1-t0))" "$win" "${rss:-0}" "$RW1" "$RW2" <<'PYEOF'
import sys, math
ref2, sc, margin, refmass, fastmass, lbl, secs, win, rss, w1, w2 = sys.argv[1:12]
ref2, margin = float(ref2), float(margin)
try: sc = float(sc)
except ValueError: print(f"  {lbl:<38} winner pair absent from the accelerated output"); sys.exit(1)
def load(p):
    d, hdr = {}, ""
    for l in open(p):
        if l.startswith("# pair"): hdr = l.strip()
        if l.startswith("#") or l.startswith("fragment"): continue
        k, v = l.rstrip("\n").split("\t"); d[k] = float(v)
    return d, hdr
rm, _ = load(refmass); fm, hdr = load(fastmass)
# The accelerated dump must describe the SAME pair as the reference dump.
if hdr and not (w1 in hdr and w2 in hdr):
    print(f"  {lbl:<38} mass dump describes a different pair: {hdr}"); sys.exit(1)
common = [k for k in rm if k in fm and rm[k] > -1e300]
mx = max(rm[k] for k in common)
num = sum(math.exp(fm[k]-rm[k]) * math.exp(rm[k]-mx) for k in common)
den = sum(math.exp(rm[k]-mx) for k in common)
retained = 100.0*num/den if den else 0.0
err = abs(sc - ref2)
print(f"  {lbl:<38} {sc:10.1f} {err:10.2f} {retained:7.1f}% {secs:>6}s {win:>6}"
      f"   safe: {'yes' if err < margin else 'NO'}, {rss} MB")
PYEOF
}

rung "1 recruitment, exact (start,end) states" r1 --rung-zero
rung "2 + midpoint dedup (16 bp)"              r2 --rung-zero --placement-dedup 16
rung "3 + start binning (64)"                  r3 --rung-zero --placement-dedup 16 --placement-bin 64
rung "4 + anchor cap (8)"                      r4 --rung-zero --placement-dedup 16 --placement-bin 64 --max-anchor-occ 8
rung "5 + placement top-k (2)"                 r5 --rung-zero --placement-dedup 16 --placement-bin 64 --max-anchor-occ 8 --placement-topk 2
echo
echo "  A rung is decision-safe only while its score error is smaller than the reference winner's"
echo "  margin. Midpoint deduplication is now its own rung: it was previously hard-coded inside what"
echo "  was called unrestricted recruitment, so any mass loss attributed to recruitment may have been"
echo "  this instead."
