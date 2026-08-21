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
# Firth's penalised likelihood gives a FINITE estimate where maximum likelihood diverges, so the
# contract is a real effect size here, not NA.
awk -v v="$lo" 'BEGIN{exit !(v+0 > 1 && v+0 < 50)}' \
  && ok "separated feature has a finite Firth effect estimate ($lo)" || bad "log_or=$lo, expected a finite Firth estimate"
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

# ------------------------------------------------- SPA and the exact support boundary
# The score statistic S = sum_i Gt_i (y_i - mu_i) is a sum of independent Bernoulli terms, so its exact
# null distribution can be convolved. These reference values come from that enumeration (R, 2^18 outcomes
# folded by convolution), and they are why SPA exists: the normal tail is 4.7x anti-conservative here.
#
#   g <- c(1,1,1,1, rep(0,14)); y <- c(1,1,1, rep(0,15))
#   -> exact two-sided p = 0.00182437858954, normal = 0.000385746755682, z = 3.5496
{ printf 'rare, A, B'; for i in $(seq 1 18); do
    if [ "$i" -le 4 ]; then printf ', 1'; else printf ', 0'; fi; done; printf '\n'; } > "$OUT/spa.bimbam"
: > "$OUT/spa.samples"; for i in $(seq 1 18); do echo "s$i" >> "$OUT/spa.samples"; done
{ printf 'sample\tphenotype\n'
  for i in $(seq 1 3);  do printf 's%d\t1\n' "$i"; done
  for i in $(seq 4 18); do printf 's%d\t0\n' "$i"; done; } > "$OUT/spa.pheno"
gzip -cf "$OUT/spa.bimbam" > "$OUT/spa.bimbam.gz"; gzip -cf "$OUT/spa.samples" > "$OUT/spa.samples.gz"
"$BIN" associate --genotypes "$OUT/spa.bimbam.gz" --samples "$OUT/spa.samples.gz" \
  --phenotype "$OUT/spa.pheno" --min-maf 0 -o "$OUT/spa" --quiet >/dev/null 2>&1
sp=$(col "$OUT/spa.assoc.tsv" p rare); sm=$(col "$OUT/spa.assoc.tsv" p_method rare)
[ "$sm" = "score_spa" ] && ok "saddlepoint used past the |z| cutoff (p_method=$sm)" \
                        || bad "p_method=$sm, expected score_spa"
# must be materially closer to the exact tail than the normal approximation is
# Pinned against an independent Lugannani-Rice implementation in R on the same quantities, so this
# catches a numerical regression rather than only a gross one. (Exact 0.00182437858954; the normal
# would give 0.000386, so here the saddlepoint is much the better of the two.)
close "$sp" 0.00132055 1e-4 && ok "SPA p matches the reference implementation ($sp)" \
                            || bad "SPA p $sp != 0.00132055"

# ---- asymmetric support: the two tails are in different regimes and must be evaluated separately.
# Covariate-dependent mu makes smax (6.0749) and -smin (6.1295) differ, so a rule that treats "at or
# beyond the boundary" as one case gets this wrong. Both thresholds here are INTERIOR, so both tails are
# saddlepoint. R reference on the same gt/mu: upper 0.000116838 + lower 2.77864e-07 = 0.000117116.
# Note the exact tail is 0.000184813 and the NORMAL would give 0.000139 -- at n=16 the saddlepoint
# expansion is not yet accurate and is the further of the two. That is a property of SPA at tiny n, not
# a defect here, which is why this asserts the implementation rather than any superiority claim.
{ printf 'asym, A, B, 2, 0, 0, 0, 2, 0, 0, 1, 2, 0, 0, 0, 2, 0, 0, 0\n'; } > "$OUT/asym.bimbam"
: > "$OUT/asym.samples"; for i in $(seq 1 16); do echo "s$i" >> "$OUT/asym.samples"; done
{ printf 'sample\tphenotype\tcv\n'
  ys=(1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0); cs=(3 1 4 1 5 9 2 6 5 3 5 8 9 7 9 3)
  for i in $(seq 1 16); do printf 's%d\t%s\t%s\n' "$i" "${ys[$((i-1))]}" "${cs[$((i-1))]}"; done; } > "$OUT/asym.pheno"
gzip -cf "$OUT/asym.bimbam" > "$OUT/asym.bimbam.gz"; gzip -cf "$OUT/asym.samples" > "$OUT/asym.samples.gz"
"$BIN" associate --genotypes "$OUT/asym.bimbam.gz" --samples "$OUT/asym.samples.gz" \
  --phenotype "$OUT/asym.pheno" --min-maf 0 -o "$OUT/asym" --quiet >/dev/null 2>&1
ap=$(col "$OUT/asym.assoc.tsv" p asym); am=$(col "$OUT/asym.assoc.tsv" p_method asym)
[ "$am" = "score_spa" ] && ok "asymmetric support, both tails interior -> saddlepoint ($am)" \
                        || bad "p_method=$am on asymmetric support, expected score_spa"
close "$ap" 0.000117116 1e-4 && ok "asymmetric two-sided p matches the reference ($ap)" \
                             || bad "asymmetric p $ap != 0.000117116"

# At the support boundary the saddlepoint is at infinity, but the probability is a single Bernoulli
# configuration and is exact: 2 * 0.5^12 = 0.00048828125 for the separation fixture above.
bp=$(col "$OUT/sep.assoc.tsv" p sep); bm=$(col "$OUT/sep.assoc.tsv" p_method sep)
[ "$bm" = "score_exact" ] && ok "support boundary uses the exact tail (p_method=$bm)" \
                          || bad "p_method=$bm at the support boundary, expected score_exact"
close "$bp" 0.00048828125 1e-6 && ok "boundary p is exact ($bp)" \
                               || bad "boundary p $bp != 0.00048828125"

# ------------------------------------------------- Firth against an independent reference
# Firth by its definition in R (penalised score U* = X'(y - mu + h(0.5 - mu)), h the leverages):
#   y <- g <- c(rep(0,6), rep(1,6))  ->  beta = 5.129898715, se = 2.241794153
close "$lo" 5.129898715 1e-5 && ok "Firth beta matches the reference ($lo)" \
                             || bad "Firth beta $lo != 5.129898715"
fse=$(col "$OUT/sep.assoc.tsv" se sep)
close "$fse" 2.241794153 1e-5 && ok "Firth se matches the reference ($fse)" \
                              || bad "Firth se $fse != 2.241794153"

# ------------------------------------------------- Li-Ji Meff on spectra with known answers
# Four 0/1 patterns repeated 4x give exactly orthogonal centred columns, so the correlation matrix is
# the identity (Meff = k); identical columns give an all-ones matrix (Meff = 1); two orthogonal blocks
# of identical columns give eigenvalues (2,2,0,0) (Meff = 2).
mk_meff() { # <name> <rows...>
  local nm="$1"; shift
  : > "$OUT/$nm.bimbam"; for r in "$@"; do echo "$r" >> "$OUT/$nm.bimbam"; done
  gzip -cf "$OUT/$nm.bimbam" > "$OUT/$nm.bimbam.gz"
  printf 'feature_id\tlayer\tencoding\tbubbles\tnodes\tsvtype\tgene\tAF\tAN\n' > "$OUT/$nm.annot.tsv"
  awk -F', ' '{print $1"\tvariant\tdosage\t.\t.\t.\t.\t.\t."}' "$OUT/$nm.bimbam" >> "$OUT/$nm.annot.tsv"
  gzip -cf "$OUT/$nm.annot.tsv" > "$OUT/$nm.annot.tsv.gz"
}
A="0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1"
B="0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1"
C="0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0"
mk_meff indep  "a, A, B, $A" "b, A, B, $B" "c, A, B, $C"
mk_meff same   "a, A, B, $A" "b, A, B, $A" "c, A, B, $A"
mk_meff blocks "a1, A, B, $A" "a2, A, B, $A" "b1, A, B, $B" "b2, A, B, $B"
: > "$OUT/m.samples"; for i in $(seq 1 16); do echo "s$i" >> "$OUT/m.samples"; done
gzip -cf "$OUT/m.samples" > "$OUT/m.samples.gz"
{ printf 'sample\tphenotype\n'
  vs=(0.3 -1.2 0.7 1.9 -0.4 0.2 -1.7 0.9 1.1 -0.6 0.4 -0.2 1.4 -1.1 0.8 -0.9)
  for i in $(seq 1 16); do printf 's%d\t%s\n' "$i" "${vs[$((i-1))]}"; done; } > "$OUT/m.pheno"
for cfg in "indep:3" "same:1" "blocks:2"; do
  nm=${cfg%%:*}; want=${cfg##*:}
  "$BIN" associate --genotypes "$OUT/$nm.bimbam.gz" --samples "$OUT/m.samples.gz" \
    --feature-annot "$OUT/$nm.annot.tsv.gz" --phenotype "$OUT/m.pheno" --min-maf 0 \
    -o "$OUT/m.$nm" --quiet >/dev/null 2>&1
  got=$(awk -F'\t' '$1=="meff_eigen"{print $2}' "$OUT/m.$nm.summary.tsv")
  [ "$got" = "$want" ] && ok "Li-Ji Meff on '$nm' = $got (expected $want)" \
                       || bad "Li-Ji Meff on '$nm' = $got, expected $want"
done

# ------------------------------------------------- input validation
dup_s="$OUT/dup.samples"; gunzip -c "$OUT/lin.samples.gz" > "$dup_s"; echo "s1" >> "$dup_s"
gzip -cf "$dup_s" > "$dup_s.gz"
e=$("$BIN" associate --genotypes "$OUT/lin.bimbam.gz" --samples "$dup_s.gz" \
     --phenotype "$OUT/lin.pheno" --min-maf 0 -o "$OUT/dv" --quiet 2>&1 | grep -c "duplicate sample id")
[ "$e" -gt 0 ] && ok "duplicate genotype sample id rejected" || bad "duplicate genotype sample id accepted"
{ head -2 "$OUT/lin.pheno"; sed -n '2p' "$OUT/lin.pheno"; } > "$OUT/dup.pheno"
e=$("$BIN" associate --genotypes "$OUT/lin.bimbam.gz" --samples "$OUT/lin.samples.gz" \
     --phenotype "$OUT/dup.pheno" --min-maf 0 -o "$OUT/dv" --quiet 2>&1 | grep -c "duplicate sample id")
[ "$e" -gt 0 ] && ok "duplicate phenotype row rejected" || bad "duplicate phenotype row accepted"
for opt in "--min-maf 2" "--ld-r2 -1" "--min-ac -5" "--cojo-p 0" "--cojo-p 1.5" "--pca -2"; do
  set -- $opt
  e=$("$BIN" associate --genotypes "$OUT/lin.bimbam.gz" --samples "$OUT/lin.samples.gz" \
       --phenotype "$OUT/lin.pheno" "$1" "$2" -o "$OUT/dv" --quiet 2>&1 | grep -c "Error")
  [ "$e" -gt 0 ] && ok "out-of-range option rejected ($opt)" || bad "$opt accepted"
done


# ---------------------------------------------------------------- output contract
# associate wrote its two files straight to their destinations: a failure between them left a new
# assoc.tsv beside a stale summary.tsv, and nothing stopped -o from naming an input.
run_assoc() { "$BIN" associate --genotypes "$OUT/lin.bimbam.gz" --samples "$OUT/lin.samples.gz" \
                --phenotype "$OUT/lin.pheno" --min-maf 0 "$@" --quiet >/dev/null 2>&1; }

# -o naming the phenotype input must be refused BEFORE the input is touched. The phenotype IS the
# file the run would write, so without the preflight the input is consumed and then replaced.
cp "$OUT/lin.pheno" "$OUT/alias.assoc.tsv"
alias_before=$(cksum < "$OUT/alias.assoc.tsv")
"$BIN" associate --genotypes "$OUT/lin.bimbam.gz" --samples "$OUT/lin.samples.gz" \
  --phenotype "$OUT/alias.assoc.tsv" --min-maf 0 -o "$OUT/alias" --quiet >/dev/null 2>&1
[ "$?" -ne 0 ] && [ "$(cksum < "$OUT/alias.assoc.tsv")" = "$alias_before" ] \
  && ok "an output aliasing the phenotype input is refused, and the input survives" \
  || bad "associate overwrote, or accepted, an output aliasing its own input"

# The family replaces atomically: a failure part-way leaves the previous pair untouched.
run_assoc -o "$OUT/txn"
a_before=$(cksum < "$OUT/txn.assoc.tsv"); s_before=$(cksum < "$OUT/txn.summary.tsv")
PANVAR_TEST_FAIL_COMMIT_AT=2 run_assoc -o "$OUT/txn" --model linear
[ "$?" -ne 0 ] && ok "an injected commit failure exits non-zero" \
               || bad "an injected commit failure exited 0"
[ "$(cksum < "$OUT/txn.assoc.tsv")" = "$a_before" ] \
  && [ "$(cksum < "$OUT/txn.summary.tsv")" = "$s_before" ] \
  && ok "both previous associate outputs survive a failed commit unchanged" \
  || bad "a failed commit replaced part of the associate family"
PANVAR_TEST_FAIL_COMMIT_AFTER_SETASIDE=2 run_assoc -o "$OUT/txn" --model linear
[ "$(cksum < "$OUT/txn.assoc.tsv" 2>/dev/null)" = "$a_before" ] \
  && [ "$(cksum < "$OUT/txn.summary.tsv" 2>/dev/null)" = "$s_before" ] \
  && ok "and also when the failure lands between set-aside and install" \
  || bad "the set-aside window lost an associate output"
res=$(ls "$OUT" | grep -c -- '-tmp\.\|-prev\.')
[ "$res" = "0" ] && ok "associate leaves no staging or reserve residue" \
                 || bad "$res residue file(s) left by associate"

echo
if [ "$fails" -eq 0 ]; then echo "associate_stats: all assertions passed"; exit 0; fi
echo "associate_stats: $fails assertion(s) FAILED"; exit 1
