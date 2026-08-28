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

# A HIGH-COPY panel as ONE THREE-ALLELE BUBBLE, not a chain of identical nodes.
#
#   source --+-- ref   : SEG + alternate spacer --+-- sink
#            +-- fewer : SEG x9  + spacer         |
#            +-- many  : SEG x10 + spacer         |
#
# The fragment scorer works on spelled sequence, so the array only has to exist INSIDE an allele --
# it does not need to be a path through repeated nodes, and representing it that way does not
# decompose: measured, a chain of ten identical nodes spelled one haplotype at 1100 bp against its
# true 1500 and another at 0 bp, and the harness still produced a full table of numbers.
#
# This still exercises everything intended: repeated syncmers occur 9-10 times within the spelled
# allele, so --max-anchor-occ 8 binds; top-k sees many implied starts; the exact scorer enumerates
# every repeat placement; and the bubble decomposition is unambiguous.
L=$(seq_of 200 3); N=$(seq_of 200 5)
SEG=$(seq_of 100 31); SPC=$(seq_of 100 32); SPCALT=$(seq_of 100 33)
rep() { local out="" i; for ((i=0;i<$1;i++)); do out="$out$SEG"; done; printf '%s' "$out"; }
A_REF="${SEG}${SPCALT}"; A_FEW="$(rep 9)${SPC}"; A_MANY="$(rep 10)${SPC}"
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\t1\t%s\n' "$L"
  printf 'S\t2\t%s\n' "$A_REF"; printf 'S\t3\t%s\n' "$A_FEW"; printf 'S\t4\t%s\n' "$A_MANY"
  printf 'S\t5\t%s\n' "$N"
  for a in 2 3 4; do printf 'L\t1\t+\t%s\t+\t0M\nL\t%s\t+\t5\t+\t0M\n' "$a" "$a"; done
  printf 'P\tref\t1+,2+,5+\t*\n'
  printf 'P\tfewer\t1+,3+,5+\t*\n'
  printf 'P\tmany\t1+,4+,5+\t*\n'
} > "$OUT/g.gfa"
"$BIN" bubble -i "$OUT/g.gfa" -r ref -o "$OUT/bub" --min-variant-bp 0 -q >/dev/null 2>&1
printf '>ref\n%s\n'   "${L}${A_REF}${N}"  > "$OUT/ref.fa"
printf '>fewer\n%s\n' "${L}${A_FEW}${N}"  > "$OUT/fewer.fa"
printf '>many\n%s\n'  "${L}${A_MANY}${N}" > "$OUT/many.fa"
{ emit_pairs "${L}${A_MANY}${N}" tA 20; emit_pairs "${L}${A_MANY}${N}" tB 20 10; } > "$OUT/reads.fa"
echo "  panel: 10-copy array against 9-copy, one bubble; sample homozygous for 10; $(grep -c '/1$' "$OUT/reads.fa") fragments"

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
        # first two fields only: the dump grew a mates_seeded and a contrib column, and unpacking
        # into exactly two names turned that into a crash rather than a reading.
        f = l.rstrip("\n").split("\t")
        if len(f) < 2: continue
        d[f[0]] = float(f[1])
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

rung "1 exact recruited states (incl. strand)"  r1 --rung-zero
rung "2 + midpoint dedup (16 bp)"              r2 --rung-zero --placement-dedup 16
rung "3 + start binning (64)"                  r3 --rung-zero --placement-dedup 16 --placement-bin 64
rung "4 + anchor cap (8)"                      r4 --rung-zero --placement-dedup 16 --placement-bin 64 --max-anchor-occ 8
rung "5 + placement top-k (2)"                 r5 --rung-zero --placement-dedup 16 --placement-bin 64 --max-anchor-occ 8 --placement-topk 2
# The cumulative order hides top-k's OWN cost: by rung 5 the anchor cap has already deleted the
# array's syncmers, so there are no placements left for top-k to truncate and its increment is zero
# for a reason that has nothing to do with top-k. Measured separately, without the cap.
rung "  top-k (2) alone, no anchor cap"        r6 --rung-zero --placement-dedup 16 --placement-bin 64 --placement-topk 2
rung "  anchor cap (8) alone, no top-k"        r7 --rung-zero --placement-dedup 16 --placement-bin 64 --max-anchor-occ 8
# The proposed replacement: no anchor discarded for being frequent, no top-k, equal-likelihood
# placements grouped with their multiplicity, and pruning by omitted MASS instead of by count.
# The ISOLATED arm: --rung-zero so exact contract exposure is used, exactly as rungs 1-3 use it.
# Invoking it without --rung-zero changes grouping AND the exposure model at once, and the 0.19-nat
# difference that produced was entirely the exposure -- not the grouping, which is exact.
rung "  multiplicity-aware compression"        r8 --rung-zero --multiplicity-aware
echo
echo "  A rung is decision-safe only while its score error is smaller than the reference winner's"
echo "  margin. Midpoint deduplication is now its own rung: it was previously hard-coded inside what"
echo "  was called unrestricted recruitment, so any mass loss attributed to recruitment may have been"
echo "  this instead."
