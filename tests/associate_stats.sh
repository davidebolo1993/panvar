#!/usr/bin/env bash
# associate_stats.sh - deterministic statistical assertions for `panvar associate`.
#
# Every expected number here is either an external reference (R's lm()) or a structural invariant, so a
# regression shows up as a failed assertion rather than a plausible-looking p-value. No RNG, no network,
# no R at run time: R's answers are hard-coded with the command that produced them.
#
#   associate_stats.sh <panvar-binary> <out-dir>
set -uo pipefail

BIN="${1:?usage: associate_stats.sh <panvar> <outdir>}"
OUT="${2:?}"
mkdir -p "$OUT"
fails=0
ok()   { printf "  ok   %s\n" "$1"; }
bad()  { printf "  FAIL %s\n" "$1"; fails=$((fails + 1)); }
# numeric comparison within a relative tolerance
close() { awk -v a="$1" -v b="$2" -v t="${3:-1e-6}" 'BEGIN{d=(a>b?a-b:b-a); r=(b<0?-b:b); exit !(d <= t*(r>1?r:1))}'; }
col()  { awk -F'\t' -v want="$2" -v id="$3" 'NR==1{for(i=1;i<=NF;i++)c[$i]=i;next} $1==id{print $(c[want])}' "$1"; }

# ---------------------------------------------------------------- linear: Student-t, not the normal
# R: g <- c(0,1,2,0,1,2,0,1,2,0,1,2); cv <- 1:12
#    y <- c(2.1,3.4,5.2,2.8,3.9,5.9,3.2,4.6,6.1,3.7,4.9,6.8)
#    summary(lm(y ~ g + cv))  ->  beta 1.353888889, se 0.06632396484, p 7.58033045421e-09 on 9 df
# The normal tail would give 1.27e-92 -- 83 orders of magnitude out, which is what makes this a test
# rather than a formality.
printf 'f1, A, B, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2\n' > "$OUT/lin.bimbam"
: > "$OUT/lin.samples"; for i in $(seq 1 12); do echo "s$i" >> "$OUT/lin.samples"; done
{ printf 'sample\tphenotype\tcv\n'
  ys=(2.1 3.4 5.2 2.8 3.9 5.9 3.2 4.6 6.1 3.7 4.9 6.8)
  for i in $(seq 1 12); do printf 's%d\t%s\t%d\n' "$i" "${ys[$((i-1))]}" "$i"; done; } > "$OUT/lin.pheno"
gzip -cf "$OUT/lin.bimbam" > "$OUT/lin.bimbam.gz"; gzip -cf "$OUT/lin.samples" > "$OUT/lin.samples.gz"
"$BIN" associate --genotypes "$OUT/lin.bimbam.gz" --samples "$OUT/lin.samples.gz" \
  --phenotype "$OUT/lin.pheno" --min-maf 0 -o "$OUT/lin" --quiet >/dev/null 2>&1
p=$(col "$OUT/lin.assoc.tsv" p f1); b=$(col "$OUT/lin.assoc.tsv" beta f1); se=$(col "$OUT/lin.assoc.tsv" se f1)
close "$b"  1.353888889    1e-6 && ok "linear beta matches R lm() ($b)"        || bad "linear beta $b != 1.353888889"
close "$se" 0.06632396484  1e-6 && ok "linear se matches R lm() ($se)"         || bad "linear se $se != 0.06632396484"
close "$p"  7.58033045421e-09 1e-4 && ok "linear p is the Student-t tail ($p)" || bad "linear p $p != 7.58033045421e-09 (normal tail would be 1.27e-92)"

# ---------------------------------------------------------------- logistic: complete separation
# The genotype predicts the outcome perfectly, so the maximum-likelihood fit diverges. The contract is:
# a valid score p, no effect estimate, and effect_status saying why -- not a silent Wald fallback.
printf 'sep, A, B, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1\ncom, A, B, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1\n' > "$OUT/sep.bimbam"
{ printf 'sample\tphenotype\n'
  for i in $(seq 1 6);  do printf 's%d\t0\n' "$i"; done
  for i in $(seq 7 12); do printf 's%d\t1\n' "$i"; done; } > "$OUT/sep.pheno"
gzip -cf "$OUT/sep.bimbam" > "$OUT/sep.bimbam.gz"
"$BIN" associate --genotypes "$OUT/sep.bimbam.gz" --samples "$OUT/lin.samples.gz" \
  --phenotype "$OUT/sep.pheno" --min-maf 0 -o "$OUT/sep" --quiet >/dev/null 2>&1
st=$(col "$OUT/sep.assoc.tsv" effect_status sep); sp=$(col "$OUT/sep.assoc.tsv" p sep)
lo=$(col "$OUT/sep.assoc.tsv" log_or sep)
[ "$st" = "separation" ]                  && ok "separated feature reports effect_status=separation" || bad "effect_status=$st, expected separation"
[ -n "$sp" ] && [ "$sp" != "NA" ]         && ok "separated feature still has a score p ($sp)"        || bad "separated feature has no p"
[ "$lo" = "NA" ]                          && ok "separated feature has no effect estimate"           || bad "log_or=$lo, expected NA under separation"
mc=$(col "$OUT/sep.assoc.tsv" mac_case com); mo=$(col "$OUT/sep.assoc.tsv" mac_ctrl com)
[ "$mc" != "." ] && [ "$mo" != "." ]      && ok "case/control carrier counts reported ($mc/$mo)"     || bad "mac_case/mac_ctrl not reported"

# ---------------------------------------------------------------- kinship validation
n=12
seq 1 $n | awk -v n=$n '{r="";for(j=1;j<=n;j++)r=r (j==$1?"1":"0") (j<n?"\t":"");print r}'  > "$OUT/k_ok.txt"
seq 1 $n | awk -v n=$n '{r="";for(j=1;j<=n;j++)r=r (j==$1?"-1":"0") (j<n?"\t":"");print r}' > "$OUT/k_neg.txt"
awk 'NR==1{print "1\t0.7" ; next} {print}' "$OUT/k_ok.txt" > "$OUT/k_ragged.txt"
awk -v n=$n 'NR==2{r="0.7";for(j=2;j<=n;j++)r=r"\t"(j==2?"1":"0");print r;next}{print}' "$OUT/k_ok.txt" > "$OUT/k_asym.txt"
for k in k_ok k_neg k_ragged k_asym; do
  err=$("$BIN" associate --genotypes "$OUT/lin.bimbam.gz" --samples "$OUT/lin.samples.gz" \
        --phenotype "$OUT/lin.pheno" --min-maf 0 --kinship "$OUT/$k.txt" --pca 1 \
        -o "$OUT/kt" --quiet 2>&1 | grep -c "Error")
  if [ "$k" = "k_ok" ]; then
    [ "$err" -eq 0 ] && ok "valid GRM accepted" || bad "valid GRM rejected"
  else
    [ "$err" -gt 0 ] && ok "invalid GRM rejected ($k)" || bad "$k accepted but is not a valid GRM"
  fi
done

echo
if [ "$fails" -eq 0 ]; then echo "associate_stats: all assertions passed"; exit 0; fi
echo "associate_stats: $fails assertion(s) FAILED"; exit 1
