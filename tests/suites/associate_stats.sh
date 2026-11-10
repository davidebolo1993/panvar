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
# numeric comparison within a relative tolerance, |a - b| <= t |b| (absolute only for b = 0), so a p-value of
# 1e-9 is checked to its own digits. Both sides must BE numbers: awk compares "NA" with a number as strings,
# and "NA" > "1.2" made d negative, so an NA result passed every comparison.
close() { awk -v a="$1" -v b="$2" -v t="${3:-1e-6}" 'BEGIN{
  num = "^[-+]?([0-9]+[.]?[0-9]*|[.][0-9]+)([eE][-+]?[0-9]+)?$"
  if (a !~ num || b !~ num) exit 1
  a += 0; b += 0; d = (a > b ? a - b : b - a); r = (b < 0 ? -b : b); exit !(d <= t * (r > 0 ? r : 1))}'; }
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
        --phenotype "$OUT/lin.pheno" --min-maf 0 --kinship "$OUT/$k.txt" --model lmm \
        -o "$OUT/kt" --quiet 2>&1 | grep -c "Error")
  if [ "$k" = "k_ok" ]; then
    [ "$err" -eq 0 ] && ok "valid GRM accepted" || bad "valid GRM rejected"
  else
    [ "$err" -gt 0 ] && ok "invalid GRM rejected ($k)" || bad "$k accepted but is not a valid GRM"
  fi
done

# ------------------------------------------------- SPA and the exact support boundary
# Few carriers make the score a coarse lattice, where the normal tail fails. Without covariates the carriers
# are tested given the number of cases, which is Fisher's exact test:
#
#   g <- c(1,1,1,1, rep(0,14)); y <- c(1,1,1, rep(0,15))
#   -> fisher.test, one tail (the other cannot be as extreme): choose(15,1) / choose(18,4) = 0.00490196
#      (normal tail 0.000386, ten times too small)
{ printf 'rare, A, B'; for i in $(seq 1 18); do
    if [ "$i" -le 4 ]; then printf ', 1'; else printf ', 0'; fi; done; printf '\n'; } > "$OUT/spa.bimbam"
: > "$OUT/spa.samples"; for i in $(seq 1 18); do echo "s$i" >> "$OUT/spa.samples"; done
{ printf 'sample\tphenotype\n'
  for i in $(seq 1 3);  do printf 's%d\t1\n' "$i"; done
  for i in $(seq 4 18); do printf 's%d\t0\n' "$i"; done; } > "$OUT/spa.pheno"
gzip -cf "$OUT/spa.bimbam" > "$OUT/spa.bimbam.gz"; gzip -cf "$OUT/spa.samples" > "$OUT/spa.samples.gz"
# 4 carriers: the tail is enumerated, and is Fisher's p itself.
"$BIN" associate --genotypes "$OUT/spa.bimbam.gz" --samples "$OUT/spa.samples.gz" \
  --phenotype "$OUT/spa.pheno" --min-maf 0 -o "$OUT/spa_x" --quiet >/dev/null 2>&1
xp=$(col "$OUT/spa_x.assoc.tsv" p rare); xm=$(col "$OUT/spa_x.assoc.tsv" p_method rare)
[ "$xm" = "score_exact_carriers" ] && ok "4 carriers: the score's tail is enumerated (p_method=$xm)" \
                                   || bad "p_method=$xm with 4 carriers, expected score_exact_carriers"
close "$xp" 0.00490196078431 1e-6 && ok "and it is Fisher's exact p ($xp)" || bad "enumerated p $xp != 0.00490196078431"
# The saddlepoint, with enumeration off.
"$BIN" associate --genotypes "$OUT/spa.bimbam.gz" --samples "$OUT/spa.samples.gz" \
  --phenotype "$OUT/spa.pheno" --min-maf 0 --exact-carriers 0 -o "$OUT/spa" --quiet >/dev/null 2>&1
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
  --phenotype "$OUT/asym.pheno" --min-maf 0 -o "$OUT/asym_x" --quiet >/dev/null 2>&1
close "$(col "$OUT/asym_x.assoc.tsv" p asym)" 0.000184813 1e-5 \
  && ok "16 samples, enumerated in full: the exact tail 0.000184813" \
  || bad "asymmetric enumerated p $(col "$OUT/asym_x.assoc.tsv" p asym) != 0.000184813"
"$BIN" associate --genotypes "$OUT/asym.bimbam.gz" --samples "$OUT/asym.samples.gz" \
  --phenotype "$OUT/asym.pheno" --min-maf 0 --exact-carriers 0 -o "$OUT/asym" --quiet >/dev/null 2>&1
ap=$(col "$OUT/asym.assoc.tsv" p asym); am=$(col "$OUT/asym.assoc.tsv" p_method asym)
[ "$am" = "score_spa" ] && ok "asymmetric support, both tails interior -> saddlepoint ($am)" \
                        || bad "p_method=$am on asymmetric support, expected score_spa"
close "$ap" 0.000117116 1e-4 && ok "asymmetric two-sided p matches the reference ($ap)" \
                             || bad "asymmetric p $ap != 0.000117116"

# The separation fixture above: its 6 carriers are its 6 cases. Without covariates the exact tail is
# Fisher's, two-sided: 2 / choose(12, 6) = 0.00216450216. With enumeration off, the saddlepoint's
# boundary rule takes the single Bernoulli configuration instead: 2 * 0.5^12 = 0.00048828125.
close "$(col "$OUT/sep.assoc.tsv" p sep)" 0.00216450216 1e-6 && ok "enumerated boundary p is Fisher's" \
  || bad "enumerated boundary p $(col "$OUT/sep.assoc.tsv" p sep) != 0.00216450216"
"$BIN" associate --genotypes "$OUT/sep.bimbam.gz" --samples "$OUT/lin.samples.gz" \
  --phenotype "$OUT/sep.pheno" --min-maf 0 --exact-carriers 0 -o "$OUT/sep_b" --quiet >/dev/null 2>&1
bp=$(col "$OUT/sep_b.assoc.tsv" p sep); bm=$(col "$OUT/sep_b.assoc.tsv" p_method sep)
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

# Dosage columns pair with --samples names by position, so a row of any other width is unusable. The
# realistic cause is a samples.txt from the other describe folder (haplotype/ vs sample/): fewer
# names than columns would pair every dosage with the wrong sample, more names would skip every row and
# write an empty table, both with exit 0, so either must be refused.
bw_run() { "$BIN" associate --genotypes "$1" --samples "$2" --phenotype "$OUT/lin.pheno" \
             --min-maf 0 -o "$OUT/bw" --quiet 2>&1; }
gunzip -c "$OUT/lin.samples.gz" | head -11 | gzip -c > "$OUT/bw11.samples.gz"
{ gunzip -c "$OUT/lin.samples.gz"; echo s13; } | gzip -c > "$OUT/bw13.samples.gz"
{ cat "$OUT/lin.bimbam"; printf 'f,2, A, B, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2\n'; } | gzip -c > "$OUT/bwcomma.bimbam.gz"
: | gzip -c > "$OUT/bwempty.bimbam.gz"
rm -f "$OUT/bw.assoc.tsv"
e=$(bw_run "$OUT/lin.bimbam.gz" "$OUT/bw11.samples.gz")
echo "$e" | grep -q "carries 12 dosage(s) but --samples .* lists 11 name(s)" \
  && ok "fewer --samples names than dosage columns rejected" \
  || bad "fewer --samples names than dosage columns accepted: $e"
e=$(bw_run "$OUT/lin.bimbam.gz" "$OUT/bw13.samples.gz")
echo "$e" | grep -q "carries 12 dosage(s) but --samples .* lists 13 name(s)" \
  && ok "more --samples names than dosage columns rejected" \
  || bad "more --samples names than dosage columns accepted: $e"
e=$(bw_run "$OUT/bwcomma.bimbam.gz" "$OUT/lin.samples.gz")
echo "$e" | grep -q "line 2 has 16 comma-separated field(s), expected 15" \
  && ok "a row with an extra field is rejected and located" \
  || bad "a row with an extra field was accepted: $e"
e=$(bw_run "$OUT/bwempty.bimbam.gz" "$OUT/lin.samples.gz")
echo "$e" | grep -q "contains no genotype rows" \
  && ok "an empty genotype file is rejected" || bad "an empty genotype file was accepted: $e"
[ ! -e "$OUT/bw.assoc.tsv" ] && ok "and none of the rejected runs wrote a result table" \
                             || bad "a rejected run still wrote $OUT/bw.assoc.tsv"


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

# ------------------------------------------------- rare carriers with a covariate
# n = 200, 20 cases, one covariate; the 3 carriers are all cases. The normal tail gives 1.6e-7 here. With
# a covariate the carriers' outcomes are enumerated at their fitted chances: all three cases is the most
# extreme outcome, so p is the product of the carriers' chances from the covariates-only fit:
#   R: mu <- fitted(glm(y ~ cv, family = binomial)); prod(mu[1:3]) = 0.00100104779
# (Fisher's exact test, ignoring the covariate: 8.7e-4).
awk 'BEGIN{printf "rare, A, B"; for(i=1;i<=200;i++) printf ", %d", (i<=3); printf "\n"}' > "$OUT/rc.bimbam"
awk 'BEGIN{for(i=1;i<=200;i++) print "s" i}' > "$OUT/rc.samples"
awk 'BEGIN{print "sample\tphenotype\tcv"; for(i=1;i<=200;i++) printf "s%d\t%d\t%.1f\n", i, (i<=20), ((i*37)%101)/10}' > "$OUT/rc.pheno"
"$BIN" associate --genotypes "$OUT/rc.bimbam" --samples "$OUT/rc.samples" --phenotype "$OUT/rc.pheno" \
  --min-maf 0 -o "$OUT/rc" --quiet >/dev/null 2>&1
close "$(col "$OUT/rc.assoc.tsv" p rare)" 0.00100104779 1e-4 \
  && ok "3 rare carriers with a covariate: enumerated p $(col "$OUT/rc.assoc.tsv" p rare) (glm: 1.001e-3)" \
  || bad "rare-carrier p $(col "$OUT/rc.assoc.tsv" p rare) != 0.00100104779"
"$BIN" associate --genotypes "$OUT/rc.bimbam" --samples "$OUT/rc.samples" --phenotype "$OUT/rc.pheno" \
  --min-maf 0 --exact-carriers 0 -o "$OUT/rc0" --quiet >/dev/null 2>&1
rm0=$(col "$OUT/rc0.assoc.tsv" p_method rare); rp0=$(col "$OUT/rc0.assoc.tsv" p rare)
[ "$rm0" = "score_spa" ] && awk -v p="$rp0" 'BEGIN{exit !(p > 1e-4)}' \
  && ok "without enumeration the saddlepoint bounds the unreachable tail, conservatively (p=$rp0)" \
  || bad "saddlepoint without enumeration: $rm0 p=$rp0, expected score_spa and p > 1e-4 (normal: 1.6e-7)"

# ------------------------------------------------- Firth by step-halving
# One carrier, a case, among 27 cases and 33 controls, two covariates. A Newton step clamped to +-5 with no
# line search cycles here (3.32, -1.68, 3.32, ...) and never converges. Reference: the penalised
# likelihood l(b) + 0.5 log|X'WX| maximised directly (scipy BFGS then Nelder-Mead):
#   beta 1.214688, se 2.341763
awk 'BEGIN{printf "firth, A, B"; for(i=1;i<=60;i++) printf ", %d", (i==5); printf "\n"}' > "$OUT/fi.bimbam"
awk 'BEGIN{for(i=1;i<=60;i++) print "s" i}' > "$OUT/fi.samples"
awk 'BEGIN{print "sample\tphenotype\tage\tsex"; for(i=1;i<=60;i++) printf "s%d\t%d\t%d\t%d\n", i, (i<=27), 30+(i*13)%40, i%2}' > "$OUT/fi.pheno"
"$BIN" associate --genotypes "$OUT/fi.bimbam" --samples "$OUT/fi.samples" --phenotype "$OUT/fi.pheno" \
  --min-maf 0 -o "$OUT/fi" --quiet >/dev/null 2>&1
close "$(col "$OUT/fi.assoc.tsv" log_or firth)" 1.214688 1e-5 && close "$(col "$OUT/fi.assoc.tsv" se firth)" 2.341763 1e-5 \
  && ok "Firth converges where a clamped step cycles (log_or $(col "$OUT/fi.assoc.tsv" log_or firth))" \
  || bad "Firth log_or/se $(col "$OUT/fi.assoc.tsv" log_or firth)/$(col "$OUT/fi.assoc.tsv" se firth) != 1.214688/2.341763"

# ------------------------------------------------- --pca: the exact eigenvectors
# K = sum_j lambda_j v_j v_j' over the 16 Walsh vectors, lambda = 5, 4.95, 4.9, ... (and 0.1 for the
# constant one). The top two are v1 and v2 exactly, so --pca 2 must equal passing them as covariates.
# Subspace iteration from a random start does not separate eigenvalues this close (p 0.130, not 0.0749).
cat > "$OUT/walsh.awk" <<'AWKEOF'
function pc(x,  c){c=0; while(x){c+=x%2; x=int(x/2)} return c}
function band(a,b,  r,m){r=0;m=1; while(a>0&&b>0){ if(a%2&&b%2) r+=m; a=int(a/2); b=int(b/2); m*=2} return r}
function v(i,j){ return ((pc(band(i,j))%2)?-1:1)/4 }
BEGIN{ n=16; lam[0]=0.1; for(j=1;j<n;j++) lam[j]=5-0.05*(j-1)
  if (mode=="K") for(a=0;a<n;a++){ line=""; for(b=0;b<n;b++){ s=0; for(j=0;j<n;j++) s+=lam[j]*v(a,j)*v(b,j); line=line (b?"\t":"") sprintf("%.12f",s)} print line }
  if (mode=="P") { print "sample\tphenotype\tpc1\tpc2"; for(a=0;a<n;a++) printf "s%d\t%.3f\t%.6f\t%.6f\n", a+1, ((a*7)%11)/3+((a*5)%3), v(a,1), v(a,2) }
  if (mode=="Y") { print "sample\tphenotype"; for(a=0;a<n;a++) printf "s%d\t%.3f\n", a+1, ((a*7)%11)/3+((a*5)%3) }
  if (mode=="G") { printf "g, A, B"; for(a=0;a<n;a++) printf ", %s", (a==3||a==11)?"NA":(a*3)%5%3; printf "\n" } }
AWKEOF
for m in K P Y G; do awk -v mode=$m -f "$OUT/walsh.awk" > "$OUT/w.$m"; done
awk 'BEGIN{for(i=1;i<=16;i++) print "s" i}' > "$OUT/w.samples"
"$BIN" associate --genotypes "$OUT/w.G" --samples "$OUT/w.samples" --phenotype "$OUT/w.Y" --kinship "$OUT/w.K" \
  --pca 2 --min-maf 0 -o "$OUT/wpca" --quiet >/dev/null 2>&1
"$BIN" associate --genotypes "$OUT/w.G" --samples "$OUT/w.samples" --phenotype "$OUT/w.P" \
  --min-maf 0 -o "$OUT/wexp" --quiet >/dev/null 2>&1
close "$(col "$OUT/wpca.assoc.tsv" p g)" "$(col "$OUT/wexp.assoc.tsv" p g)" 1e-5 \
  && ok "--pca 2 equals the top two eigenvectors as covariates (p $(col "$OUT/wpca.assoc.tsv" p g))" \
  || bad "--pca 2 p $(col "$OUT/wpca.assoc.tsv" p g) != explicit-PC p $(col "$OUT/wexp.assoc.tsv" p g)"
# The LMM fits every sample, its two missing dosages mean-imputed, and counts the 14 it observed.
"$BIN" associate --genotypes "$OUT/w.G" --samples "$OUT/w.samples" --phenotype "$OUT/w.Y" --kinship "$OUT/w.K" \
  --model lmm --min-maf 0 -o "$OUT/wlmm" --quiet >/dev/null 2>&1
[ "$(col "$OUT/wlmm.assoc.tsv" n g)" = "14" ] && ok "LMM n counts the observed dosages (14 of 16)" \
  || bad "LMM n = $(col "$OUT/wlmm.assoc.tsv" n g), expected 14"
{ cat "$OUT/w.K"; awk 'NR==1' "$OUT/w.K"; } > "$OUT/w.K17"
"$BIN" associate --genotypes "$OUT/w.G" --samples "$OUT/w.samples" --phenotype "$OUT/w.Y" --kinship "$OUT/w.K17" \
  --pca 1 --min-maf 0 -o "$OUT/wk17" --quiet >/dev/null 2>&1
[ "$?" -ne 0 ] && ok "a kinship matrix with a row more than the samples is refused" \
               || bad "a 17-row kinship matrix was accepted for 16 samples"

# ------------------------------------------------- inputs that cannot be analysed
refuse() { # refuse <label> <phenotype file>
  "$BIN" associate --genotypes "$OUT/lin.bimbam.gz" --samples "$OUT/lin.samples.gz" --phenotype "$2" \
    --min-maf 0 -o "$OUT/rf" --quiet >/dev/null 2>&1
  [ "$?" -ne 0 ] && ok "$1 is refused" || bad "$1 was accepted"
}
mkph() { awk -v hdr="$1" -v expr="$2" 'BEGIN{print hdr; for(i=1;i<=12;i++){ y=(i*7)%5+0.5*i; c=i; print "s" i "\t" eval_(expr,i,y,c)}}
  function eval_(e,i,y,c){ if(e=="const") return "5"; if(e=="one") return "1"; if(e=="dup") return y "\t" y;
    if(e=="ccov") return y "\t3"; if(e=="coll") return y "\t" c "\t" (3*c+1); if(e=="junk") return (i==4?"1.5abc":y) }' > "$OUT/ph.$2"; }
mkph $'sample\tphenotype' const;          refuse "a constant phenotype" "$OUT/ph.const"
mkph $'sample\tphenotype' one;            refuse "a binary phenotype with one class" "$OUT/ph.one"
mkph $'sample\tphenotype\tcv' dup;        refuse "the phenotype given again as a covariate" "$OUT/ph.dup"
mkph $'sample\tphenotype\tcv' ccov;       refuse "a constant covariate" "$OUT/ph.ccov"
mkph $'sample\tphenotype\tc0\tc1' coll;   refuse "two collinear covariates" "$OUT/ph.coll"
mkph $'sample\tphenotype' junk;           refuse "a phenotype cell that is not a number (1.5abc)" "$OUT/ph.junk"
# A feature whose observed samples all share one phenotype value has nothing to fit: residuals of
# rounding noise must not give it a beta and a p.
printf 'fc, A, B, 0, 1, 2, 1, 0, 2, 1, 3, NA, NA, NA, NA\n' > "$OUT/fc.bimbam"
{ printf 'sample\tphenotype\n'; for i in $(seq 1 8); do printf 's%d\t1\n' $i; done
  for i in $(seq 9 12); do printf 's%d\t%d\n' $i $((i % 3)); done; } > "$OUT/fc.pheno"
"$BIN" associate --genotypes "$OUT/fc.bimbam" --samples "$OUT/lin.samples" --phenotype "$OUT/fc.pheno" \
  --min-maf 0 --model linear -o "$OUT/fc" --quiet >/dev/null 2>&1
[ "$(awk -F'\t' '$1=="features_tested"{print $2}' "$OUT/fc.summary.tsv")" = "0" ] \
  && ok "a feature whose samples share one phenotype value is not tested" \
  || bad "a feature over a constant phenotype was tested"
# The singularity test is on the covariates' correlation, not their units: lm()'s numbers again with
# the covariate at 1e-9 of its scale (an absolute pivot floor of 1e-12 dropped every feature).
{ printf 'sample\tphenotype\tcv\n'; ys=(2.1 3.4 5.2 2.8 3.9 5.9 3.2 4.6 6.1 3.7 4.9 6.8)
  for i in $(seq 1 12); do printf 's%d\t%s\t%s\n' "$i" "${ys[$((i-1))]}" "${i}e-9"; done; } > "$OUT/lin7.pheno"
"$BIN" associate --genotypes "$OUT/lin.bimbam.gz" --samples "$OUT/lin.samples.gz" --phenotype "$OUT/lin7.pheno" \
  --min-maf 0 -o "$OUT/lin7" --quiet >/dev/null 2>&1
close "$(col "$OUT/lin7.assoc.tsv" p f1)" 7.58033045421e-09 1e-4 \
  && ok "a covariate at 1e-9 of its scale gives lm()'s p unchanged" \
  || bad "covariate at scale 1e-9: p $(col "$OUT/lin7.assoc.tsv" p f1), expected 7.58033045421e-09"
# The PLINK/GEMMA "-9" is missing, as is any spelling of the number -9 ("-9.0" below).
{ printf 'sample\tphenotype\tcv\n'; ys=(-9 3.4 5.2 2.8 3.9 5.9 3.2 4.6 6.1 3.7 4.9 6.8)
  for i in $(seq 1 12); do printf 's%d\t%s\t%d\n' "$i" "${ys[$((i-1))]}" "$i"; done; } > "$OUT/na9.pheno"
"$BIN" associate --genotypes "$OUT/lin.bimbam.gz" --samples "$OUT/lin.samples.gz" --phenotype "$OUT/na9.pheno" \
  --min-maf 0 -o "$OUT/na9" --quiet >/dev/null 2>&1
[ "$(awk -F'\t' '$1=="samples_used"{print $2}' "$OUT/na9.summary.tsv")" = "11" ] \
  && ok "a phenotype of -9 is missing (11 samples used)" || bad "-9 was not read as missing"

# ------------------------------------------------- low_af and Meff
# Exactly --min-ac (3) carriers of 30 is not low: 1 - 27/30 is 0.0999..., which times 30 fell under 3.
awk 'BEGIN{printf "v1, A, B"; for(i=1;i<=30;i++) printf ", %d", (i<=3); printf "\n"}' > "$OUT/ac.bimbam"
awk 'BEGIN{for(i=1;i<=30;i++) print "s" i}' > "$OUT/ac.samples"
awk 'BEGIN{print "sample\tphenotype"; for(i=1;i<=30;i++) printf "s%d\t%.2f\n", i, (i*7)%11/3}' > "$OUT/ac.pheno"
"$BIN" associate --genotypes "$OUT/ac.bimbam" --samples "$OUT/ac.samples" --phenotype "$OUT/ac.pheno" \
  --unit variant --min-maf 0 -o "$OUT/acx" --quiet >/dev/null 2>&1
[ "$(col "$OUT/acx.assoc.tsv" low_af v1)" = "0" ] && ok "exactly --min-ac carriers is not low_af" \
  || bad "3 carriers of 30 flagged low_af=$(col "$OUT/acx.assoc.tsv" low_af v1) at --min-ac 3"
# The feature unit has no effective-test count: its family-wise threshold is --permutations.
printf 'feature_id\tlayer\tencoding\tbubbles\tnodes\nf1\tkmer\tsyncmer\t5;6;8\t1;2;3\n' | gzip -c > "$OUT/km.annot.gz"
"$BIN" associate --genotypes "$OUT/lin.bimbam.gz" --samples "$OUT/lin.samples.gz" --phenotype "$OUT/lin.pheno" \
  --feature-annot "$OUT/km.annot.gz" --min-maf 0 -o "$OUT/km" --quiet >/dev/null 2>&1
[ "$(awk -F'\t' '$1=="meff"{print $2}' "$OUT/km.summary.tsv"):$(col "$OUT/km.assoc.tsv" p_bonf_meff f1)" = "NA:." ] \
  && ok "feature unit: no Meff, p_bonf_meff '.'" \
  || bad "feature unit gave Meff $(awk -F'\t' '$1=="meff"{print $2}' "$OUT/km.summary.tsv")"

# ------------------------------------------------- conditioning on partly missing records; scale and parsing
sv() { awk -F'\t' -v k="$2" '$1==k{print $2}' "$1"; }
# A record nested in a deletion is missing ('.') exactly for the deletion's carriers. Both carry an
# effect, and each is an independent signal: N given D is N's own test (D is 0 on N's samples), D given
# N enters N with its mean-filled dosage. Requiring every conditioning variant observed would make both
# fits fail, leaving one signal and a shadow with p_conditional NA. Python statsmodels OLS:
#   N on its 28 observed samples p 0.0006694058699; y ~ D + N(mean-filled) p_D 0.0001558621697
awk 'BEGIN{ printf "D, A, B"; for(i=1;i<=40;i++) printf ", %d", (i<=12); printf "\n"
  printf "N, A, B"; for(i=1;i<=40;i++) printf ", %s", (i<=12 ? "NA" : (i<=23 && i%2) ? 1 : 0); printf "\n" }' > "$OUT/nest.bimbam"
awk 'BEGIN{for(i=1;i<=40;i++) print "s" i}' > "$OUT/nest.samples"
awk 'BEGIN{print "sample\tphenotype"; for(i=1;i<=40;i++){ n=(i>12 && i<=23 && i%2); printf "s%d\t%.10g\n", i, 1.5*(i<=12)+1.5*n+((i*7)%11)/4-1.2 }}' > "$OUT/nest.pheno"
"$BIN" associate --genotypes "$OUT/nest.bimbam" --samples "$OUT/nest.samples" --phenotype "$OUT/nest.pheno" \
  --unit variant --min-maf 0 -o "$OUT/nest" --quiet >/dev/null 2>&1
[ "$(col "$OUT/nest.assoc.tsv" cond_role D)" = "signal" ] && [ "$(col "$OUT/nest.assoc.tsv" cond_role N)" = "signal" ] \
  && close "$(col "$OUT/nest.assoc.tsv" p_conditional N)" 0.0006694058699 1e-6 \
  && close "$(col "$OUT/nest.assoc.tsv" p_conditional D)" 0.0001558621697 1e-5 \
  && [ "$(col "$OUT/nest.assoc.tsv" n_conditional D)" = "40" ] \
  && ok "a deletion and the record nested in it are two signals (conditional p as statsmodels)" \
  || bad "nested record: roles $(col "$OUT/nest.assoc.tsv" cond_role D)/$(col "$OUT/nest.assoc.tsv" cond_role N), p_conditional $(col "$OUT/nest.assoc.tsv" p_conditional D)/$(col "$OUT/nest.assoc.tsv" p_conditional N)"
# Two rows that test the same given the signals (max = sum - min, so given the sum, min and max tie): the one
# stronger on its own enters, not the one first in the file (max is written before min). The phenotype follows
# both haplotype values, the smaller more: the sum enters, then min (stronger than max on its own), and max is
# collinear.
# Haplotype values 0-3 from a fixed linear congruential sequence, so the sum leaves min and max free.
awk -v PH="$OUT/tie.pheno" 'BEGIN{ x = 12345
  for (i = 1; i <= 60; i++) { x = (x * 1103 + 12345) % 65536; a[i] = int(x / 4096) % 4
                              x = (x * 1103 + 12345) % 65536; b[i] = int(x / 4096) % 4 }
  printf "T:max, A, B"; for (i = 1; i <= 60; i++) printf ", %d", (a[i] > b[i] ? a[i] : b[i]); printf "\n"
  printf "T, A, B"; for (i = 1; i <= 60; i++) printf ", %d", a[i] + b[i]; printf "\n"
  printf "T:min, A, B"; for (i = 1; i <= 60; i++) printf ", %d", (a[i] < b[i] ? a[i] : b[i]); printf "\n"
  print "sample\tphenotype" > PH
  for (i = 1; i <= 60; i++) printf "s%d\t%.10g\n", i, 1.2 * (a[i] < b[i] ? a[i] : b[i]) + 0.6 * (a[i] > b[i] ? a[i] : b[i]) + ((i * 13) % 17) / 20 > PH }' \
  > "$OUT/tie.bimbam"
awk 'BEGIN{for(i=1;i<=60;i++) print "s" i}' > "$OUT/tie.samples"
"$BIN" associate --genotypes "$OUT/tie.bimbam" --samples "$OUT/tie.samples" --phenotype "$OUT/tie.pheno" \
  --unit variant --min-maf 0 --cojo-p 0.05 -o "$OUT/tie" --quiet >/dev/null 2>&1
[ "$(col "$OUT/tie.assoc.tsv" cond_role T:min)" = "signal" ] && [ "$(col "$OUT/tie.assoc.tsv" cond_role T)" = "signal" ] \
  && [ "$(col "$OUT/tie.assoc.tsv" cond_role T:max)" = "collinear" ] \
  && ok "of two rows tied given the signals (min and max given the sum), the one stronger on its own enters" \
  || bad "min/max tie: T:min $(col "$OUT/tie.assoc.tsv" cond_role T:min), T:max $(col "$OUT/tie.assoc.tsv" cond_role T:max), T $(col "$OUT/tie.assoc.tsv" cond_role T)"
# No signal selected: nothing was conditioned on, so no row is a shadow.
"$BIN" associate --genotypes "$OUT/nest.bimbam" --samples "$OUT/nest.samples" --phenotype "$OUT/nest.pheno" \
  --unit variant --min-maf 0 --cojo-p 1e-12 -o "$OUT/nest0" --quiet >/dev/null 2>&1
[ "$(col "$OUT/nest0.assoc.tsv" cond_role D)" = "." ] && [ "$(sv "$OUT/nest0.summary.tsv" cojo_independent_signals)" = "0" ] \
  && ok "with no signal selected the roles are '.'" || bad "no signal selected: D is '$(col "$OUT/nest0.assoc.tsv" cond_role D)'"
# lambda_gc is read back from p: lm()'s p 7.58033045421e-09 is chi-square 33.3792 (scipy chi2.isf),
# lambda 73.3727; the t statistic squared (20.41^2) would give 916.
close "$(sv "$OUT/lin.summary.tsv" lambda_gc)" 73.37267376 1e-6 && ok "lambda_gc comes from p, not from t^2" \
  || bad "lambda_gc $(sv "$OUT/lin.summary.tsv" lambda_gc), expected 73.37267376"
# A covariate far from 0 relative to its spread (1e4 + i/200) is i again after centring: lm()'s p.
{ printf 'sample\tphenotype\tcv\n'; ys=(2.1 3.4 5.2 2.8 3.9 5.9 3.2 4.6 6.1 3.7 4.9 6.8)
  for i in $(seq 1 12); do printf 's%d\t%s\t%s\n' "$i" "${ys[$((i-1))]}" "$(awk -v i=$i 'BEGIN{printf "%.3f", 10000+i/200}')"; done; } > "$OUT/lin8.pheno"
"$BIN" associate --genotypes "$OUT/lin.bimbam.gz" --samples "$OUT/lin.samples.gz" --phenotype "$OUT/lin8.pheno" \
  --min-maf 0 -o "$OUT/lin8" --quiet >/dev/null 2>&1
close "$(col "$OUT/lin8.assoc.tsv" p f1)" 7.58033045421e-09 1e-4 \
  && ok "a covariate of 1e4 + small spread is not the intercept again (lm()'s p)" \
  || bad "covariate 1e4 + i/200: p '$(col "$OUT/lin8.assoc.tsv" p f1)', expected 7.58033045421e-09"
# The LMM's variance ratio is in K's units: 1e6 K and 1e-6 K give the same test. Here delta is near 0
# (3.5e-5 K units), so a fixed [1e-5, 1e5] grid would floor it for 1e-6 K (p 0.19561, not 0.195025).
for sc in 1e6 1e-6; do
  awk -v s=$sc '{for(i=1;i<=NF;i++) $i=sprintf("%.12g", $i*s)} 1' OFS='\t' "$OUT/w.K" > "$OUT/w.K$sc"
  "$BIN" associate --genotypes "$OUT/w.G" --samples "$OUT/w.samples" --phenotype "$OUT/w.Y" --kinship "$OUT/w.K$sc" \
    --model lmm --min-maf 0 -o "$OUT/wlmm$sc" --quiet >/dev/null 2>&1
  close "$(col "$OUT/wlmm$sc.assoc.tsv" p g)" "$(col "$OUT/wlmm.assoc.tsv" p g)" 1e-6 \
    && ok "the LMM gives the same p for K and $sc K" \
    || bad "LMM p $(col "$OUT/wlmm$sc.assoc.tsv" p g) for $sc K, $(col "$OUT/wlmm.assoc.tsv" p g) for K"
done
# PLINK case/control coding 1/2 is binary: the logistic model, as for 0/1.
awk -F'\t' 'NR==1{print;next}{$2=$2+1; print}' OFS='\t' "$OUT/fi.pheno" > "$OUT/fi12.pheno"
"$BIN" associate --genotypes "$OUT/fi.bimbam" --samples "$OUT/fi.samples" --phenotype "$OUT/fi12.pheno" \
  --min-maf 0 -o "$OUT/fi12" --quiet >/dev/null 2>&1
[ "$(sv "$OUT/fi12.summary.tsv" model)" = "logistic" ] && [ "$(col "$OUT/fi12.assoc.tsv" p firth)" = "$(col "$OUT/fi.assoc.tsv" p firth)" ] \
  && ok "a 1/2-coded phenotype is PLINK case/control" || bad "a 1/2 phenotype ran model $(sv "$OUT/fi12.summary.tsv" model)"
# A sidecar's extra columns are read by name, in any order: the place on the reference reaches assoc.tsv.
printf 'feature_id\tlayer\tencoding\tbubbles\tnodes\tend\tchrom\tsvtype\tpos\nD\tvariant\tdosage\t1\t1\t2500\tchr9\tDEL\t1001\nN\tvariant\tdosage\t1\t2\t.\t.\tSNV\t.\n' > "$OUT/pos.annot"
"$BIN" associate --genotypes "$OUT/nest.bimbam" --samples "$OUT/nest.samples" --phenotype "$OUT/nest.pheno" \
  --feature-annot "$OUT/pos.annot" --unit variant --min-maf 0 -o "$OUT/pos" --quiet >/dev/null 2>&1
check_pos="$(col "$OUT/pos.assoc.tsv" chrom D):$(col "$OUT/pos.assoc.tsv" pos D)-$(col "$OUT/pos.assoc.tsv" end D):$(col "$OUT/pos.assoc.tsv" pos N)"
[ "$check_pos" = "chr9:1001-2500:." ] && ok "chrom, pos and end come from the sidecar by name; '.' where it has none" \
  || bad "sidecar place on the reference: $check_pos"
# Refusals: a sidecar that annotates none of the genotype rows, a duplicated feature id, and an
# --exact-carriers whose 2^N table would not fit.
printf 'feature_id\tlayer\tencoding\tbubbles\tnodes\tsvtype\tgene\taf\tan\nX1\tvariant\t.\t1\t1\tDEL\t.\t0.5\t40\n' > "$OUT/x.annot"
"$BIN" associate --genotypes "$OUT/nest.bimbam" --samples "$OUT/nest.samples" --phenotype "$OUT/nest.pheno" \
  --feature-annot "$OUT/x.annot" --min-maf 0 -o "$OUT/rf" --quiet >/dev/null 2>&1
[ "$?" -ne 0 ] && ok "a --feature-annot that annotates no genotype row is refused" || bad "a non-matching --feature-annot was accepted"
{ cat "$OUT/nest.bimbam"; head -1 "$OUT/nest.bimbam"; } > "$OUT/dup.bimbam"
"$BIN" associate --genotypes "$OUT/dup.bimbam" --samples "$OUT/nest.samples" --phenotype "$OUT/nest.pheno" \
  --min-maf 0 -o "$OUT/rf" --quiet >/dev/null 2>&1
[ "$?" -ne 0 ] && ok "a duplicated feature id is refused" || bad "a duplicated feature id was accepted"
"$BIN" associate --genotypes "$OUT/fi.bimbam" --samples "$OUT/fi.samples" --phenotype "$OUT/fi.pheno" \
  --exact-carriers 21 -o "$OUT/rf" --quiet >/dev/null 2>&1
[ "$?" -ne 0 ] && ok "--exact-carriers 21 is refused" || bad "--exact-carriers 21 was accepted"

# --genotypes is read more than once: a pipe would give the conditional pass an empty input (all NA, exit 0).
"$BIN" associate --genotypes <(cat "$OUT/nest.bimbam") --samples "$OUT/nest.samples" --phenotype "$OUT/nest.pheno" \
  --unit variant --min-maf 0 -o "$OUT/rf" --quiet >/dev/null 2>&1
[ "$?" -ne 0 ] && ok "--genotypes from a pipe is refused" || bad "--genotypes from a pipe was accepted"
# Numeric options are parsed strictly: "0.05%" is not 0.05, nor "4.9" 4.
for bad_opt in "--min-maf 0.05%" "--exact-carriers 4.9" "--cojo-p 1e-3junk"; do
  # shellcheck disable=SC2086
  "$BIN" associate --genotypes "$OUT/nest.bimbam" --samples "$OUT/nest.samples" --phenotype "$OUT/nest.pheno" \
    $bad_opt -o "$OUT/rf" --quiet >/dev/null 2>&1
  [ "$?" -ne 0 ] && ok "$bad_opt is refused" || bad "$bad_opt was accepted"
done

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

# ---------------------------------------------------------------- numerical edge cases, inputs, permutations
L="--genotypes $OUT/lin.bimbam.gz --samples $OUT/lin.samples.gz --min-maf 0 --quiet"
# A near-perfect fit is a result, not rounding noise: y = 2 g + 1e-6 * (pattern) gives beta 2.
{ printf 'sample\tphenotype\n'; i=0; for gv in 0 1 2 0 1 2 0 1 2 0 1 2; do i=$((i+1));
  printf 's%d\t%.10f\n' "$i" "$(awk -v g=$gv -v i=$i 'BEGIN{print 2*g + 1e-6*((i*7)%5)}')"; done; } > "$OUT/near.pheno"
"$BIN" associate $L --phenotype "$OUT/near.pheno" -o "$OUT/near" >/dev/null 2>&1
close "$(col "$OUT/near.assoc.tsv" beta f1)" 2 1e-5 && ok "a near-perfect fit is tested (beta 2)" \
  || bad "near-perfect fit: beta '$(col "$OUT/near.assoc.tsv" beta f1)'"
# The LMM does not depend on where the phenotype sits: + 1e9 gives the same p.
awk 'NR==1{print;next}{printf "%s\t%.6f\t%s\n", $1, $2 + 1e9, $3}' "$OUT/lin.pheno" > "$OUT/shift.pheno"
"$BIN" associate $L --phenotype "$OUT/lin.pheno" --model lmm --kinship "$OUT/k_ok.txt" -o "$OUT/lm0" >/dev/null 2>&1
"$BIN" associate $L --phenotype "$OUT/shift.pheno" --model lmm --kinship "$OUT/k_ok.txt" -o "$OUT/lm9" >/dev/null 2>&1
close "$(col "$OUT/lm9.assoc.tsv" p f1)" "$(col "$OUT/lm0.assoc.tsv" p f1)" 1e-4 \
  && ok "LMM p unchanged by a phenotype offset" || bad "LMM offset: $(col "$OUT/lm9.assoc.tsv" p f1) vs $(col "$OUT/lm0.assoc.tsv" p f1)"
# Tied kinship eigenvalues at the --pca cut leave the PCs undefined.
"$BIN" associate $L --phenotype "$OUT/lin.pheno" --kinship "$OUT/k_ok.txt" --pca 1 -o "$OUT/tie" > "$OUT/tie.err" 2>&1
grep -q 'not unique' "$OUT/tie.err" && ok "tied PCs are refused" || bad "tied PCs accepted"
# A phenotype of -9 is missing; a short row is refused; a trailing tab is no column.
awk 'NR==2{$2="-9.0"} {print}' OFS='\t' "$OUT/lin.pheno" > "$OUT/m9.pheno"
"$BIN" associate $L --phenotype "$OUT/m9.pheno" -o "$OUT/m9" >/dev/null 2>&1
[ "$(col "$OUT/m9.assoc.tsv" n f1)" = "11" ] && ok "-9 phenotype read as missing (n 11)" || bad "-9: n $(col "$OUT/m9.assoc.tsv" n f1)"
awk 'NR==3{print $1"\t"$2; next} {print}' "$OUT/lin.pheno" > "$OUT/short.pheno"
"$BIN" associate $L --phenotype "$OUT/short.pheno" -o "$OUT/short" >/dev/null 2>&1
[ "$?" -ne 0 ] && ok "a phenotype row shorter than the header is refused" || bad "short phenotype row accepted"
awk '{print $0"\t"}' "$OUT/lin.pheno" > "$OUT/tab.pheno"
"$BIN" associate $L --phenotype "$OUT/tab.pheno" -o "$OUT/tab" >/dev/null 2>&1
close "$(col "$OUT/tab.assoc.tsv" p f1)" 7.58033045421e-09 1e-5 && ok "a trailing tab is no column" || bad "trailing tab: p '$(col "$OUT/tab.assoc.tsv" p f1)'"
# The statistic column is `stat`; permutations give p_perm, reproducibly, and are refused with the LMM.
head -1 "$OUT/lin.assoc.tsv" | tr '\t' '\n' | grep -qx stat && ok "the statistic column is stat" || bad "no stat column"
"$BIN" associate $L --phenotype "$OUT/lin.pheno" --permutations 99 --seed 3 -o "$OUT/pm1" >/dev/null 2>&1
"$BIN" associate $L --phenotype "$OUT/lin.pheno" --permutations 99 --seed 3 -o "$OUT/pm2" >/dev/null 2>&1
pp=$(col "$OUT/pm1.assoc.tsv" p_perm f1)
close "$pp" 0.01 1e-9 && ok "a p beyond every shuffle has p_perm 1/(N+1)" || bad "p_perm $pp"
cmp -s "$OUT/pm1.assoc.tsv" "$OUT/pm2.assoc.tsv" && ok "permutations are reproducible with --seed" || bad "permutation output differs between runs"
[ "$(col "$OUT/lin.assoc.tsv" p_perm f1)" = "." ] && ok "p_perm is '.' without --permutations" || bad "p_perm without permutations"
# The LMM's shuffles run in the rotated coordinates where the samples are independent.
"$BIN" associate $L --phenotype "$OUT/lin.pheno" --model lmm --kinship "$OUT/k_ok.txt" --permutations 20 -o "$OUT/pml" >/dev/null 2>&1
[ "$?" -eq 0 ] && [ -n "$(awk -F'\t' '$1=="permutation_threshold" && $2!="NA"{print $2}' "$OUT/pml.summary.tsv")" ] \
  && ok "--permutations with the LMM gives a threshold" || bad "--permutations with the LMM failed"
# Non-numeric node ids reach their genes, through an edge key too.
printf 'feature_id\tlayer\tencoding\tbubbles\tnodes\nf1\tgraph\tedge\t1\tv1+>v2-\n' | gzip -c > "$OUT/vg.annot.gz"
printf 'node_id\tgenes\nv1\tGENEA\nv2\tGENEB\n1\tWRONG\n' > "$OUT/vg.node_genes.tsv"
"$BIN" associate $L --phenotype "$OUT/lin.pheno" --feature-annot "$OUT/vg.annot.gz" --node-genes "$OUT/vg.node_genes.tsv" -o "$OUT/vg" >/dev/null 2>&1
[ "$(col "$OUT/vg.assoc.tsv" gene f1)" = "GENEA;GENEB" ] && ok "node ids are read whole (v1, not 1)" || bad "genes '$(col "$OUT/vg.assoc.tsv" gene f1)'"
# Help: a bare command and -h anywhere exit 0.
"$BIN" associate >/dev/null 2>&1; r1=$?; "$BIN" associate --min-maf 0 -h >/dev/null 2>&1; r2=$?
[ "$r1:$r2" = "0:0" ] && ok "bare command and -h anywhere exit 0" || bad "help exit codes $r1:$r2"
# Logistic estimates do not depend on a covariate's units.
awk 'BEGIN{print "sample\tphenotype\tc"; for(i=1;i<=40;i++) printf "s%d\t%d\t%.6f\n", i, ((i*7)%10<4), (i*3)%11}' > "$OUT/lg1.pheno"
awk 'NR==1{print;next}{printf "%s\t%s\t%.15g\n", $1, $2, $3*1e-9}' OFS='\t' "$OUT/lg1.pheno" > "$OUT/lg9.pheno"
awk 'BEGIN{printf "g1, A, B"; for(i=1;i<=40;i++) printf ", %d", (i*5)%3; printf "\n"}' > "$OUT/lg.bimbam"
awk 'BEGIN{for(i=1;i<=40;i++) print "s" i}' > "$OUT/lg.samples"
for t in 1 9; do "$BIN" associate --genotypes "$OUT/lg.bimbam" --samples "$OUT/lg.samples" --phenotype "$OUT/lg$t.pheno" \
  --min-maf 0 -o "$OUT/lg$t" --quiet >/dev/null 2>&1; done
close "$(col "$OUT/lg9.assoc.tsv" log_or g1)" "$(col "$OUT/lg1.assoc.tsv" log_or g1)" 1e-6 \
  && ok "logistic log_or independent of covariate units" || bad "log_or $(col "$OUT/lg1.assoc.tsv" log_or g1) vs $(col "$OUT/lg9.assoc.tsv" log_or g1)"
# Feature-unit conditioning keeps the target's samples; a missing lead dosage is mean-imputed.
awk 'BEGIN{printf "a, A, B"; for(i=1;i<=20;i++) printf ", %s", (i==1?"NA":(i%2)); printf "\n";
           printf "b, A, B"; for(i=1;i<=20;i++) printf ", %d", (i%3==0); printf "\n"}' > "$OUT/fc.bimbam"
awk 'BEGIN{for(i=1;i<=20;i++) print "s" i}' > "$OUT/fc.samples"
awk 'BEGIN{print "sample\tphenotype"; for(i=1;i<=20;i++) printf "s%d\t%.3f\n", i, 3*(i%2) + 0.7*(i%3==0) + ((i*7)%5)/10}' > "$OUT/fc.pheno"
"$BIN" associate --genotypes "$OUT/fc.bimbam" --samples "$OUT/fc.samples" --phenotype "$OUT/fc.pheno" \
  --unit feature --min-maf 0 -o "$OUT/fc" --quiet >/dev/null 2>&1
[ "$(col "$OUT/fc.assoc.tsv" n_conditional b)" = "$(col "$OUT/fc.assoc.tsv" n b)" ] \
  && ok "feature conditional n equals marginal n ($(col "$OUT/fc.assoc.tsv" n b))" \
  || bad "feature conditional n $(col "$OUT/fc.assoc.tsv" n_conditional b) vs $(col "$OUT/fc.assoc.tsv" n b)"
# The modal dosage is found on the values themselves, so a rescaled matrix has the same minor frequency.
awk 'BEGIN{printf "r, A, B"; for(i=1;i<=12;i++) printf ", %d", 2+(i%5); printf "\n";
           printf "s, A, B"; for(i=1;i<=12;i++) printf ", %.6f", (2+(i%5))*2/6; printf "\n"}' > "$OUT/sc.bimbam"
"$BIN" associate --genotypes "$OUT/sc.bimbam" --samples "$OUT/lin.samples" --phenotype "$OUT/lin.pheno" --min-maf 0 -o "$OUT/sc" --quiet >/dev/null 2>&1
[ "$(col "$OUT/sc.assoc.tsv" minor_freq r)" = "$(col "$OUT/sc.assoc.tsv" minor_freq s)" ] \
  && ok "minor_freq unchanged by rescaling ($(col "$OUT/sc.assoc.tsv" minor_freq r))" \
  || bad "minor_freq $(col "$OUT/sc.assoc.tsv" minor_freq r) vs $(col "$OUT/sc.assoc.tsv" minor_freq s)"
# Li-Ji with more records than samples: 20 records, copies of two uncorrelated vectors, Meff 2.
awk 'BEGIN{for(k=1;k<=20;k++){printf "v%d, A, B", k; for(i=1;i<=12;i++) printf ", %d", (k<=10 ? (i%2) : (int((i-1)/2)%2)); printf "\n"}}' > "$OUT/lj.bimbam"
"$BIN" associate --genotypes "$OUT/lj.bimbam" --samples "$OUT/lin.samples" --phenotype "$OUT/lin.pheno" --unit variant \
  --min-maf 0 -o "$OUT/lj" --quiet >/dev/null 2>&1
[ "$(awk -F'\t' '$1=="meff_eigen"{print $2}' "$OUT/lj.summary.tsv")" = "2" ] && ok "Li-Ji through the samples' Gram matrix: Meff 2" \
  || bad "Li-Ji Meff $(awk -F'\t' '$1=="meff_eigen"{print $2}' "$OUT/lj.summary.tsv")"

# Several matrices are joined by sample name: a second matrix in another column order, with a sample the
# first lacks and without one it has, tests like one matrix over the shared samples.
awk 'BEGIN{for(k=1;k<=3;k++){printf "m%d, A, B", k; for(i=1;i<=12;i++) printf ", %d", (k==1 ? i%3 : k==2 ? (i%4>1) : (i%5)%3); printf "\n"}}' > "$OUT/mj.all"
head -n 2 "$OUT/mj.all" > "$OUT/mj.a.bimbam"
awk -F', ' 'NR==3{printf "%s, %s, %s, 1", $1, $2, $3; for(i=14;i>=4;i--) printf ", %s", $i; printf "\n"}' "$OUT/mj.all" > "$OUT/mj.b.bimbam"
{ echo s99; for i in $(seq 11 -1 1); do echo "s$i"; done; } > "$OUT/mj.b.samples"
awk -F', ' '{printf "%s, %s, %s", $1, $2, $3; for(i=4;i<=14;i++) printf ", %s", $i; printf "\n"}' "$OUT/mj.all" > "$OUT/mj.one.bimbam"
head -n 11 "$OUT/lin.samples" > "$OUT/mj.one.samples"
printf 'feature_id\tlayer\tencoding\tbubbles\tnodes\nm1\tvariant\tdosage\t1\t1\nm2\tvariant\tdosage\t1\t2\nm3\tsnv\tdosage\t.\t.\n' > "$OUT/mj.annot"
"$BIN" associate --genotypes "$OUT/mj.a.bimbam" --samples "$OUT/lin.samples" --genotypes "$OUT/mj.b.bimbam" \
  --samples "$OUT/mj.b.samples" --feature-annot "$OUT/mj.annot" --phenotype "$OUT/lin.pheno" --min-maf 0 \
  -o "$OUT/mj" --quiet >/dev/null 2>&1
"$BIN" associate --genotypes "$OUT/mj.one.bimbam" --samples "$OUT/mj.one.samples" --feature-annot "$OUT/mj.annot" \
  --phenotype "$OUT/lin.pheno" --min-maf 0 -o "$OUT/mj1" --quiet >/dev/null 2>&1
cmp -s "$OUT/mj.assoc.tsv" "$OUT/mj1.assoc.tsv" && [ -s "$OUT/mj.assoc.tsv" ] \
  && ok "two matrices joined by sample name test as one over the shared samples" \
  || bad "joined matrices differ from the single matrix over the shared samples"
[ "$(cut -f1 "$OUT/mj.layers.tsv" | sed 1d | sort | paste -sd, -)" = "snv,variant" ] \
  && ok "layers.tsv has one row per layer" || bad "layers.tsv rows: $(cut -f1 "$OUT/mj.layers.tsv" | paste -sd, -)"
lr=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)c[$i]=i;next} {print $(c["best_feature"]) ":" $(c["r2_with_top"])}' "$OUT/mj.layers.tsv" | head -1)
[ "${lr#*:}" = "1" ] && ok "layers.tsv: the top layer's r2_with_top is 1 ($lr)" || bad "layers.tsv top row $lr"
"$BIN" associate --genotypes "$OUT/mj.a.bimbam" --samples "$OUT/lin.samples" --genotypes "$OUT/mj.b.bimbam" \
  --phenotype "$OUT/lin.pheno" -o "$OUT/mjx" --quiet >/dev/null 2>&1 \
  && bad "a --genotypes without its --samples was accepted" || ok "one --samples per --genotypes is required"
"$BIN" associate --genotypes "$OUT/mj.a.bimbam" --samples "$OUT/lin.samples" --genotypes "$OUT/mj.a.bimbam" \
  --samples "$OUT/lin.samples" --phenotype "$OUT/lin.pheno" -o "$OUT/mjx" --quiet >/dev/null 2>&1 \
  && bad "an id in two matrices was accepted" || ok "an id in two matrices is refused"

# Two records both past the p floor (1e-300): the stronger statistic enters the selection, whichever
# comes first in the file.
awk 'BEGIN{n=400; for(i=1;i<=n;i++){g[i]=i%3; w[i]=(i<=2 ? (g[i]+1)%3 : g[i])}
  printf "weak, A, B"; for(i=1;i<=n;i++) printf ", %d", w[i]; printf "\n"
  printf "strong, A, B"; for(i=1;i<=n;i++) printf ", %d", g[i]; printf "\n"
  for(i=1;i<=n;i++) print "s" i > "'"$OUT"'/fl.samples"
  print "sample\tphenotype" > "'"$OUT"'/fl.pheno"
  for(i=1;i<=n;i++) printf "s%d\t%.4f\n", i, g[i] + ((i*37)%11-5)/100 > "'"$OUT"'/fl.pheno"}' > "$OUT/fl.bimbam"
"$BIN" associate --genotypes "$OUT/fl.bimbam" --samples "$OUT/fl.samples" --phenotype "$OUT/fl.pheno" --unit variant \
  --min-maf 0 -o "$OUT/fl" --quiet >/dev/null 2>&1
[ "$(col "$OUT/fl.assoc.tsv" p weak)" = "1e-300" ] && [ "$(col "$OUT/fl.assoc.tsv" cond_role strong)" = "signal" ] \
  && [ "$(col "$OUT/fl.assoc.tsv" cond_role weak)" = "shadow" ] \
  && ok "a tie at the p floor goes to the larger statistic" \
  || bad "floor tie: weak p $(col "$OUT/fl.assoc.tsv" p weak) role $(col "$OUT/fl.assoc.tsv" cond_role weak), strong role $(col "$OUT/fl.assoc.tsv" cond_role strong)"

[ "$(awk -F'\t' 'NR==2{print $1}' "$OUT/fl.assoc.tsv")" = "strong" ] && [ "$(col "$OUT/fl.assoc.tsv" is_lead strong)" = "1" ] \
  && ok "the same tie orders assoc.tsv and leads the clump" \
  || bad "floor tie: first row $(awk -F'\t' 'NR==2{print $1}' "$OUT/fl.assoc.tsv"), strong is_lead $(col "$OUT/fl.assoc.tsv" is_lead strong)"
# Past the floor, rows rank on their unfloored p, which depends on the degrees of freedom: `short`, observed
# on 400 samples with less noise, has the larger t, but `broad`, on all 1000, the smaller p.
awk 'BEGIN{n=1000; for(i=1;i<=n;i++){g[i]=i%3; e[i]=((i*37)%11-5)/5*(i<=400 ? 0.001 : 0.02)}
  printf "short, A, B"; for(i=1;i<=n;i++) printf ", %s", (i<=400 ? g[i] : "NA"); printf "\n"
  printf "broad, A, B"; for(i=1;i<=n;i++) printf ", %d", g[i]; printf "\n"
  for(i=1;i<=n;i++) print "s" i > "'"$OUT"'/df.samples"
  print "sample\tphenotype" > "'"$OUT"'/df.pheno"
  for(i=1;i<=n;i++) printf "s%d\t%.6f\n", i, g[i] + e[i] > "'"$OUT"'/df.pheno"}' > "$OUT/df.bimbam"
"$BIN" associate --genotypes "$OUT/df.bimbam" --samples "$OUT/df.samples" --phenotype "$OUT/df.pheno" --unit variant \
  --min-maf 0 -o "$OUT/df" --quiet >/dev/null 2>&1
ts=$(col "$OUT/df.assoc.tsv" stat short); tb=$(col "$OUT/df.assoc.tsv" stat broad)
awk -v a="$ts" -v b="$tb" 'BEGIN{exit !(a > b)}' && [ "$(awk -F'\t' 'NR==2{print $1}' "$OUT/df.assoc.tsv")" = "broad" ] \
  && [ "$(col "$OUT/df.assoc.tsv" cond_role broad)" = "signal" ] && [ "$(col "$OUT/df.layers.tsv" best_feature .)" = "broad" ] \
  && ok "past the floor, the smaller unfloored p ranks first (t $tb on 1000 over t $ts on 400)" \
  || bad "df ranking: first $(awk -F'\t' 'NR==2{print $1}' "$OUT/df.assoc.tsv"), t short $ts broad $tb"
# An id repeated across matrices is refused even when one copy would be dropped by --min-maf.
awk 'BEGIN{printf "dupx, A, B"; for(i=1;i<=12;i++) printf ", 0"; printf "\n"}' > "$OUT/dupa.bimbam"
awk 'BEGIN{printf "dupx, A, B"; for(i=1;i<=12;i++) printf ", %d", i%3; printf "\n"}' > "$OUT/dupb.bimbam"
"$BIN" associate --genotypes "$OUT/dupa.bimbam" --samples "$OUT/lin.samples" --genotypes "$OUT/dupb.bimbam" \
  --samples "$OUT/lin.samples" --phenotype "$OUT/lin.pheno" --unit feature -o "$OUT/dupx" --quiet > "$OUT/dupx.err" 2>&1
grep -q "has the feature id 'dupx'" "$OUT/dupx.err" \
  && ok "an id repeated in a filtered row is refused" || bad "an id repeated in a filtered row was accepted"
# The same id twice across --feature-annot files is refused.
"$BIN" associate --genotypes "$OUT/mj.a.bimbam" --samples "$OUT/lin.samples" --feature-annot "$OUT/mj.annot" \
  --feature-annot "$OUT/mj.annot" --phenotype "$OUT/lin.pheno" -o "$OUT/mjx" --quiet > "$OUT/mjx.err" 2>&1
grep -q "lists the id 'm1' more than once" "$OUT/mjx.err" \
  && ok "an id in two --feature-annot files is refused" || bad "an id in two --feature-annot files was accepted"
# --unit auto reads the annotations of the rows tested: sidecar rows of k-mers not in the matrix do not count.
{ cat "$OUT/mj.annot"; for k in $(seq 1 10); do printf 'K%d\tkmer\tsyncmer\t1\t1\n' "$k"; done; } > "$OUT/mjk.annot"
"$BIN" associate --genotypes "$OUT/mj.one.bimbam" --samples "$OUT/mj.one.samples" --feature-annot "$OUT/mjk.annot" \
  --phenotype "$OUT/lin.pheno" --min-maf 0 -o "$OUT/mjk" --quiet >/dev/null 2>&1
[ "$(awk -F'\t' '$1=="unit"{print $2}' "$OUT/mjk.summary.tsv")" = "variant" ] \
  && ok "--unit auto ignores annotations of rows not in the matrix" \
  || bad "--unit auto chose $(awk -F'\t' '$1=="unit"{print $2}' "$OUT/mjk.summary.tsv")"
# The COJO entry is Bonferroni by default and the permutation threshold with --permutations.
src=$(awk -F'\t' '$1=="cojo_entry_source"{print $2}' "$OUT/mj1.summary.tsv")
"$BIN" associate --genotypes "$OUT/mj.one.bimbam" --samples "$OUT/mj.one.samples" --feature-annot "$OUT/mj.annot" \
  --phenotype "$OUT/lin.pheno" --min-maf 0 --permutations 99 -o "$OUT/mjp" --quiet >/dev/null 2>&1
srcp=$(awk -F'\t' '$1=="cojo_entry_source"{print $2}' "$OUT/mjp.summary.tsv")
thr=$(awk -F'\t' '$1=="permutation_threshold"{print $2}' "$OUT/mjp.summary.tsv"); ent=$(awk -F'\t' '$1=="cojo_entry_threshold"{print $2}' "$OUT/mjp.summary.tsv")
[ "$src" = "bonferroni" ] && [ "$srcp" = "permutation" ] && [ "$thr" = "$ent" ] \
  && ok "COJO entry: Bonferroni by default, the permutation threshold with --permutations ($ent)" \
  || bad "COJO entry sources $src / $srcp, threshold $thr vs entry $ent"
# p_perm = (1 + k) / (B + 1): under 20 shuffles no p_perm can fall below 0.05, so the selection falls back to
# Bonferroni. With 19, p_perm >= 1/20, so a floating-point boundary must not admit the smallest shuffle
# minimum either; with 20, the smallest is 1/21.
for B in 10 19 20; do
  "$BIN" associate --genotypes "$OUT/fl.bimbam" --samples "$OUT/fl.samples" --phenotype "$OUT/fl.pheno" --unit variant \
    --min-maf 0 --permutations "$B" -o "$OUT/fp$B" > /dev/null 2> "$OUT/fp$B.err"
done
[ "$(sv "$OUT/fp10.summary.tsv" cojo_entry_source)" = "bonferroni" ] && [ "$(sv "$OUT/fp19.summary.tsv" cojo_entry_source)" = "bonferroni" ] \
  && [ "$(sv "$OUT/fp19.summary.tsv" permutation_threshold)" = "0" ] && [ "$(col "$OUT/fp19.assoc.tsv" cond_role strong)" = "signal" ] \
  && ok "under 20 shuffles the selection uses Bonferroni" \
  || bad "few shuffles: sources $(sv "$OUT/fp10.summary.tsv" cojo_entry_source) / $(sv "$OUT/fp19.summary.tsv" cojo_entry_source), threshold(19) $(sv "$OUT/fp19.summary.tsv" permutation_threshold)"
[ "$(sv "$OUT/fp20.summary.tsv" cojo_entry_source)" = "permutation" ] && close "$(col "$OUT/fp20.assoc.tsv" p_perm strong)" 0.0476190476 1e-6 \
  && ok "with 20 shuffles the permutation threshold is used, and p_perm is at least 1/21" \
  || bad "20 shuffles: source $(sv "$OUT/fp20.summary.tsv" cojo_entry_source), p_perm $(col "$OUT/fp20.assoc.tsv" p_perm strong)"
grep -q "19 permutations cannot give p_perm below 0.05 (20 or more needed)" "$OUT/fp19.err" && ! grep -q "cannot give" "$OUT/fp20.err" \
  && ok "the fallback warning asks for 20 shuffles, and is silent at 20" || bad "fallback warning: $(grep 'cannot give' "$OUT/fp19.err")"
"$BIN" associate --genotypes "$OUT/lin.bimbam.gz" --samples "$OUT/lin.samples.gz" --phenotype "$OUT/lin.pheno" --unit variant \
  --min-maf 0 --model lmm --kinship "$OUT/k_ok.txt" --permutations 19 -o "$OUT/fpl" > /dev/null 2> "$OUT/fpl.err"
[ "$(sv "$OUT/fpl.summary.tsv" cojo_entry_source)" = "bonferroni" ] && grep -q "cannot give p_perm below 0.05" "$OUT/fpl.err" \
  && ok "the LMM's selection falls back to Bonferroni the same way, and says so" \
  || bad "LMM, 19 shuffles: source $(sv "$OUT/fpl.summary.tsv" cojo_entry_source), warning '$(grep 'cannot give' "$OUT/fpl.err")'"
# Logistic tails past the floor: each row separates its observed cases and controls exactly, `weak` on 350
# cases and 800 controls (log10 p -306.9, larger score), `strong` on 520 and 520 (log10 p -312.8).
awk 'BEGIN{n=1600
  printf "weak, A, B"; for(i=0;i<n;i++) printf ", %s", (i<350 ? 1 : (i>=600 && i<1400 ? 0 : "NA")); printf "\n"
  printf "strong, A, B"; for(i=0;i<n;i++) printf ", %s", (i<520 ? 1 : (i>=600 && i<1120 ? 0 : "NA")); printf "\n"
  for(i=0;i<n;i++) print "s" i > "'"$OUT"'/lg.samples"
  print "sample\tphenotype" > "'"$OUT"'/lg.pheno"
  for(i=0;i<n;i++) printf "s%d\t%d\n", i, (i<600) > "'"$OUT"'/lg.pheno"}' > "$OUT/lg.bimbam"
"$BIN" associate --genotypes "$OUT/lg.bimbam" --samples "$OUT/lg.samples" --phenotype "$OUT/lg.pheno" --unit variant \
  --min-maf 0 -o "$OUT/lg" --quiet >/dev/null 2>&1
[ "$(awk -F'\t' 'NR==2{print $1}' "$OUT/lg.assoc.tsv")" = "strong" ] && [ "$(col "$OUT/lg.assoc.tsv" cond_role strong)" = "signal" ] \
  && ok "logistic corrected tails rank on their unfloored p" \
  || bad "logistic floor: first $(awk -F'\t' 'NR==2{print $1}' "$OUT/lg.assoc.tsv"), strong role $(col "$OUT/lg.assoc.tsv" cond_role strong)"
# --unit auto: rows without an annotation count as features, so one record among 30 unannotated rows is
# tested per feature.
{ cat "$OUT/mj.one.bimbam"; for k in $(seq 1 30); do awk -v k="$k" 'BEGIN{printf "U%d, A, B", k; for(i=1;i<=11;i++) printf ", %d", (i*k)%3; printf "\n"}'; done; } > "$OUT/mju.bimbam"
"$BIN" associate --genotypes "$OUT/mju.bimbam" --samples "$OUT/mj.one.samples" --feature-annot "$OUT/mj.annot" \
  --phenotype "$OUT/lin.pheno" --min-maf 0 -o "$OUT/mju" --quiet >/dev/null 2>&1
[ "$(sv "$OUT/mju.summary.tsv" unit)" = "feature" ] && ok "--unit auto counts unannotated rows as features" \
  || bad "--unit auto with unannotated rows chose $(sv "$OUT/mju.summary.tsv" unit)"
# Under the LMM, t lives in rotated coordinates: no best_r2.
awk 'BEGIN{for(i=1;i<=12;i++){for(j=1;j<=12;j++) printf "%s%s", (j>1?" ":""), (i==j?1:(int((i-1)/4)==int((j-1)/4)?0.25:0)); printf "\n"}}' > "$OUT/lmm.K"
"$BIN" associate --genotypes "$OUT/lin.bimbam.gz" --samples "$OUT/lin.samples.gz" --phenotype "$OUT/lin.pheno" \
  --model lmm --kinship "$OUT/lmm.K" --min-maf 0 -o "$OUT/lmmr" --quiet >/dev/null 2>&1
[ "$(col "$OUT/lmmr.layers.tsv" best_r2 .)" = "NA" ] && ok "no best_r2 under the LMM" \
  || bad "LMM best_r2 $(col "$OUT/lmmr.layers.tsv" best_r2 .)"
# The LMM runs the conditional selection in its rotated coordinates; without shuffles, at Bonferroni.
"$BIN" associate --genotypes "$OUT/lin.bimbam.gz" --samples "$OUT/lin.samples.gz" --phenotype "$OUT/lin.pheno" \
  --model lmm --kinship "$OUT/lmm.K" --unit variant --min-maf 0 -o "$OUT/lmmv" --quiet >/dev/null 2>&1
[ "$(sv "$OUT/lmmv.summary.tsv" cojo_entry_source)" = "bonferroni" ] && [ "$(sv "$OUT/lmmv.summary.tsv" cojo_independent_signals)" = "1" ] \
  && ok "the LMM's conditional selection keeps the one signal" \
  || bad "LMM selection: $(sv "$OUT/lmmv.summary.tsv" cojo_entry_source) $(sv "$OUT/lmmv.summary.tsv" cojo_independent_signals)"

# --haplotypes: per-haplotype columns, in an order unrelated to the samples', combined per sample. Sample i
# carries f1's dosage g as haplotypes (0,0), (1,0) or (0,1), (1,1); s12's second haplotype is two contigs.
# In r1 one of them is observed, so s12 is known; in r2 both are, so s12 is missing (n 11).
awk 'BEGIN{
  for (i = 12; i >= 1; i--) { c[++n] = "s" i "a"; if (i == 12) { c[++n] = "s12b_1"; c[++n] = "s12b_2" } else c[++n] = "s" i "b" }
  for (k = 1; k <= n; k++) print c[k] > "'"$OUT"'/hx.columns"
  print "sample\thap1\thap2" > "'"$OUT"'/hx.table"
  for (i = 1; i <= 12; i++) print "s" i "\ts" i "a\t" (i == 12 ? "s12b_1,s12b_2" : "s" i "b") > "'"$OUT"'/hx.table"
  for (r = 1; r <= 2; r++) {
    printf "r%d, A, B", r
    for (k = 1; k <= n; k++) {
      split(c[k], q, /[ab]/); i = substr(q[1], 2) + 0; g = (i - 1) % 3; side = (c[k] ~ /a$/) ? 1 : 2
      if (c[k] == "s12b_2") { printf ", %s", (r == 1 ? "NA" : "0"); continue }
      v = (g == 2) ? 1 : (g == 0 ? 0 : ((i % 2 == 1) == (side == 1) ? 1 : 0))
      printf ", %d", v }
    printf "\n" } }' > "$OUT/hx.bimbam"
# The same rows written per sample: the sum is f1, min is g == 2, max is g >= 1.
awk 'BEGIN{ split("sum min max", e, " ")
  for (k = 1; k <= 3; k++) { printf "r1%s, A, B", (k == 1 ? "" : ":" e[k])
    for (i = 1; i <= 12; i++) { g = (i - 1) % 3; printf ", %d", (k == 1 ? g : (k == 2 ? g == 2 : g >= 1)) }
    printf "\n" } }' > "$OUT/hs.bimbam"
"$BIN" associate --genotypes "$OUT/hx.bimbam" --samples "$OUT/hx.columns" --haplotypes "$OUT/hx.table" \
  --encoding sum,min,max --phenotype "$OUT/lin.pheno" --min-maf 0 -o "$OUT/hx" --quiet >/dev/null 2>&1
"$BIN" associate --genotypes "$OUT/hs.bimbam" --samples "$OUT/lin.samples" --phenotype "$OUT/lin.pheno" \
  --min-maf 0 -o "$OUT/hs" --quiet >/dev/null 2>&1
hx_ok=1
for id in r1 r1:min r1:max; do
  close "$(col "$OUT/hx.assoc.tsv" p "$id")" "$(col "$OUT/hs.assoc.tsv" p "$id")" 1e-9 || hx_ok=0
  [ "$(col "$OUT/hx.assoc.tsv" n "$id")" = 12 ] || hx_ok=0
done
[ $hx_ok = 1 ] && close "$(col "$OUT/hx.assoc.tsv" p r1)" 7.58033045421e-09 1e-4 \
  && ok "--haplotypes: sum, min and max per sample, as the same rows written per sample" \
  || bad "--haplotypes encodings: p r1 $(col "$OUT/hx.assoc.tsv" p r1) / $(col "$OUT/hs.assoc.tsv" p r1), min $(col "$OUT/hx.assoc.tsv" p r1:min) / $(col "$OUT/hs.assoc.tsv" p r1:min)"
[ "$(col "$OUT/hx.assoc.tsv" n r2)" = 11 ] && ok "--haplotypes: a haplotype with two observed contigs is missing" \
  || bad "--haplotypes: r2 n $(col "$OUT/hx.assoc.tsv" n r2), expected 11"
"$BIN" associate --genotypes "$OUT/lin.bimbam.gz" --samples "$OUT/lin.samples.gz" --encoding min \
  --phenotype "$OUT/lin.pheno" --min-maf 0 -o "$OUT/hy" --quiet > "$OUT/hy.err" 2>&1
grep -q "needs --haplotypes" "$OUT/hy.err" && ok "--encoding min without --haplotypes is refused" \
  || bad "--encoding min without --haplotypes: $(tail -n 1 "$OUT/hy.err")"
"$BIN" associate --genotypes "$OUT/hx.bimbam" --samples "$OUT/lin.samples" --haplotypes "$OUT/hx.table" \
  --phenotype "$OUT/lin.pheno" --min-maf 0 -o "$OUT/hz" --quiet > "$OUT/hz.err" 2>&1
grep -q "which --samples .* does not list" "$OUT/hz.err" && ok "--haplotypes naming a column the matrix lacks is refused" \
  || bad "--haplotypes with a sample-level matrix: $(tail -n 1 "$OUT/hz.err")"

# Fine-mapping. 400 samples: three groups of four correlated rows (A, B, C) and eight independent ones;
# the phenotype follows A1 and B1. Reference: a numpy port of susieR's IBSS with refinement on the same
# rows (L 10, estimated prior and residual variances): PIP A1 0.999409, B1 0.999318, B2 0.0707531, and
# credible sets {A1} and {B1}; on a phenotype of noise alone, no credible set.
awk -v F="$OUT" '
# Deterministic fixture: 400 samples, 3 clusters of 4 correlated 0/1/2 rows + 8 independent rows.
# The phenotype follows A1 and B1 (clusters A, B); cluster C and the rest are decoys. LCG noise.
function rnd() { seed = (seed * 1103515245 + 12345) % 2147483648; return seed / 2147483648 }
function gauss() { return sqrt(-2 * log(rnd() + 1e-12)) * cos(6.283185307179586 * rnd()) }
function geno(f) { return (rnd() < f) + (rnd() < f) }
BEGIN {
  seed = 7; n = 400
  for (i = 1; i <= n; i++) {
    for (c = 1; c <= 3; c++) { base = geno(0.3); for (k = 1; k <= 4; k++) g[c, k, i] = (rnd() < 0.12 ? geno(0.3) : base) }
    for (k = 1; k <= 8; k++) h[k, i] = geno(0.25)
    y[i] = 0.5 * g[1, 1, i] + 0.4 * g[2, 1, i] + gauss()
    z[i] = gauss()
  }
  split("A B C", nm, " ")
  for (c = 1; c <= 3; c++) for (k = 1; k <= 4; k++) { printf "%s%d, A, B", nm[c], k > (F "/fm.bimbam"); for (i = 1; i <= n; i++) printf ", %d", g[c, k, i] > (F "/fm.bimbam"); printf "\n" > (F "/fm.bimbam") }
  for (k = 1; k <= 8; k++) { printf "N%d, A, B", k > (F "/fm.bimbam"); for (i = 1; i <= n; i++) printf ", %d", h[k, i] > (F "/fm.bimbam"); printf "\n" > (F "/fm.bimbam") }
  for (i = 1; i <= n; i++) print "s" i > (F "/fm.samples")
  print "sample\tphenotype" > (F "/fm.pheno"); for (i = 1; i <= n; i++) printf "s%d\t%.6f\n", i, y[i] > (F "/fm.pheno")
  print "sample\tphenotype" > (F "/fm0.pheno"); for (i = 1; i <= n; i++) printf "s%d\t%.6f\n", i, z[i] > (F "/fm0.pheno")
}
'
fm_run() { "$BIN" associate --genotypes "$OUT/fm.bimbam" --samples "$OUT/fm.samples" --unit variant --min-maf 0 --quiet "$@" >/dev/null 2>&1; }
fm_run --phenotype "$OUT/fm.pheno" -o "$OUT/fm"
fm_run --phenotype "$OUT/fm0.pheno" -o "$OUT/fm0"
fm_run --phenotype "$OUT/fm.pheno" --finemap-l 0 -o "$OUT/fmoff"
close "$(col "$OUT/fm.assoc.tsv" pip A1)" 0.999409 1e-4 && close "$(col "$OUT/fm.assoc.tsv" pip B1)" 0.999318 1e-4 \
  && close "$(col "$OUT/fm.assoc.tsv" pip B2)" 0.0707531 1e-3 \
  && ok "fine-mapping PIPs match the susieR port" \
  || bad "PIPs A1 $(col "$OUT/fm.assoc.tsv" pip A1), B1 $(col "$OUT/fm.assoc.tsv" pip B1), B2 $(col "$OUT/fm.assoc.tsv" pip B2)"
[ "$(awk -F'\t' 'NR > 1 { printf "%s;", $8 }' "$OUT/fm.cs.tsv")" = "A1;B1;" ] && [ "$(col "$OUT/fm.assoc.tsv" cs B1)" = 2 ] \
  && ok "one credible set per causal row" || bad "credible sets: $(awk -F'\t' 'NR > 1 { printf "%s;", $8 }' "$OUT/fm.cs.tsv")"
[ "$(sv "$OUT/fm0.summary.tsv" credible_sets)" = 0 ] && ok "no credible set on noise" \
  || bad "noise phenotype: $(sv "$OUT/fm0.summary.tsv" credible_sets) credible sets"
[ "$(col "$OUT/fmoff.assoc.tsv" pip A1)" = "." ] && [ "$(sv "$OUT/fmoff.summary.tsv" credible_sets)" = "NA" ] \
  && ok "--finemap-l 0 turns fine-mapping off" || bad "--finemap-l 0: pip $(col "$OUT/fmoff.assoc.tsv" pip A1)"

# Rare-variant groups. 1000 samples; gene GA has four rare records, three raising the phenotype; GB four
# null ones; GC one record, where burden and SKAT are that record's score test; and common rows.
awk -v F="$OUT" '
function rnd() { seed = (seed * 1103515245 + 12345) % 2147483648; return seed / 2147483648 }
function gauss() { return sqrt(-2 * log(rnd() + 1e-12)) * cos(6.283185307179586 * rnd()) }
BEGIN {
  seed = 5; n = 1000; split("GA GB GC", gn, " "); split("4 4 1", gm, " ")
  print "feature_id\tlayer\tencoding\tbubbles\tnodes\tsvtype\tgene\tAF\tAN" > (F "/gr.annot")
  for (g = 1; g <= 3; g++) for (k = 1; k <= gm[g]; k++) {
    id = gn[g] "_" k; printf "%s, A, B", id > (F "/gr.bimbam")
    for (i = 1; i <= n; i++) { x = (rnd() < 0.004); D[g, k, i] = x; printf ", %d", x > (F "/gr.bimbam") }
    printf "\n" > (F "/gr.bimbam"); print id "\tvariant\tdosage\t" g "\t.\tSNV\t" gn[g] "\t.\t." > (F "/gr.annot") }
  for (k = 1; k <= 5; k++) { printf "C%d, A, B", k > (F "/gr.bimbam"); for (i = 1; i <= n; i++) printf ", %d", (rnd() < 0.3) + (rnd() < 0.3) > (F "/gr.bimbam"); printf "\n" > (F "/gr.bimbam")
    print "C" k "\tvariant\tdosage\t9\t.\tSNV\t.\t.\t." > (F "/gr.annot") }
  print "sample\tphenotype" > (F "/gr.pheno")
  for (i = 1; i <= n; i++) { print "s" i > (F "/gr.samples"); printf "s%d\t%.6f\n", i, 1.2 * (D[1, 1, i] + D[1, 2, i] + D[1, 3, i]) + gauss() > (F "/gr.pheno") } }'
"$BIN" associate --genotypes "$OUT/gr.bimbam" --samples "$OUT/gr.samples" --feature-annot "$OUT/gr.annot" --phenotype "$OUT/gr.pheno" \
  --unit variant --min-maf 0 --group gene -o "$OUT/gr" --quiet >/dev/null 2>&1
gcol() { awk -F'\t' -v want="$2" -v id="$3" 'NR==1{for(i=1;i<=NF;i++)c[$i]=i;next} $1==id{print $(c[want])}' "$1"; }
close "$(gcol "$OUT/gr.groups.tsv" burden_p GC)" "$(gcol "$OUT/gr.groups.tsv" skat_p GC)" 1e-9 \
  && [ "$(gcol "$OUT/gr.groups.tsv" skato_p GC)" = "$(gcol "$OUT/gr.groups.tsv" skat_p GC)" ] \
  && close "$(gcol "$OUT/gr.groups.tsv" skat_p GC)" "$(col "$OUT/gr.assoc.tsv" p GC_1)" 2e-2 \
  && ok "a one-record group: burden, SKAT and SKAT-O are the record's own test" \
  || bad "one-record group: burden $(gcol "$OUT/gr.groups.tsv" burden_p GC), SKAT $(gcol "$OUT/gr.groups.tsv" skat_p GC), row p $(col "$OUT/gr.assoc.tsv" p GC_1)"
awk -v a="$(gcol "$OUT/gr.groups.tsv" skato_p GA)" -v b="$(gcol "$OUT/gr.groups.tsv" skato_p GB)" 'BEGIN{exit !(a < 1e-4 && b > 0.01)}' \
  && ok "SKAT-O finds the gene with the effect, not the null one" \
  || bad "SKAT-O GA $(gcol "$OUT/gr.groups.tsv" skato_p GA), GB $(gcol "$OUT/gr.groups.tsv" skato_p GB)"
# Without shuffles the selection keeps GA_1 and GA_3 (Bonferroni over 14 rows); given them, GA adds little.
[ "$(sv "$OUT/gr.summary.tsv" groups_conditioned_on)" = 2 ] \
  && awk -v c="$(gcol "$OUT/gr.groups.tsv" skato_p_conditional GA)" 'BEGIN{exit !(c > 0.01)}' \
  && ok "groups tested again given the selected signals" \
  || bad "conditional groups: on $(sv "$OUT/gr.summary.tsv" groups_conditioned_on), GA $(gcol "$OUT/gr.groups.tsv" skato_p_conditional GA)"
# The same groups under the logistic model (phenotypes above 0.524 as cases, 313 of 1000) and the LMM (families of four): the gene
# with the effect still leads, every p is a probability.
awk -F'\t' 'NR == 1 { print; next } { print $1 "\t" ($2 > 0.524 ? 1 : 0) }' "$OUT/gr.pheno" > "$OUT/grb.pheno"
awk 'BEGIN{n=1000; for(i=1;i<=n;i++){for(j=1;j<=n;j++) printf "%s%s", (j>1?" ":""), (i==j?1:(int((i-1)/4)==int((j-1)/4)?0.5:0)); printf "\n"}}' > "$OUT/gr.K"
"$BIN" associate --genotypes "$OUT/gr.bimbam" --samples "$OUT/gr.samples" --feature-annot "$OUT/gr.annot" --phenotype "$OUT/grb.pheno" \
  --unit variant --min-maf 0 --group gene --finemap-l 0 -o "$OUT/grb" --quiet >/dev/null 2>&1
"$BIN" associate --genotypes "$OUT/gr.bimbam" --samples "$OUT/gr.samples" --feature-annot "$OUT/gr.annot" --phenotype "$OUT/gr.pheno" \
  --unit variant --min-maf 0 --group gene --finemap-l 0 --model lmm --kinship "$OUT/gr.K" -o "$OUT/grl" --quiet >/dev/null 2>&1
grp_ok() { awk -F'\t' 'NR == 1 { next } { for (c = 4; c <= 6; c++) if (!($c >= 0 && $c <= 1)) bad = 1; if (NR == 2 && $1 != "GA") bad = 1 } END { exit bad || NR < 2 }' "$1"; }
grp_ok "$OUT/grb.groups.tsv" && grp_ok "$OUT/grl.groups.tsv" && ok "group tests under the logistic model and the LMM" \
  || bad "logistic/LMM groups: $(awk -F'\t' 'NR == 2 { print $1, $6 }' "$OUT/grb.groups.tsv") / $(awk -F'\t' 'NR == 2 { print $1, $6 }' "$OUT/grl.groups.tsv")"
"$BIN" associate --genotypes "$OUT/hx.bimbam" --samples "$OUT/hx.columns" --haplotypes "$OUT/hx.table" --encoding min,max \
  --group all --phenotype "$OUT/lin.pheno" --min-maf 0 -o "$OUT/hg" --quiet > "$OUT/hg.err" 2>&1
grep -q "must include sum" "$OUT/hg.err" && ok "--group without the sum encoding is refused" || bad "--group with --encoding min,max: $(tail -n 1 "$OUT/hg.err")"
# A selected signal that separates cases from controls leaves no logistic null model for the conditional test.
"$BIN" associate --genotypes "$OUT/sep.bimbam.gz" --samples "$OUT/lin.samples.gz" --phenotype "$OUT/sep.pheno" --unit variant \
  --group all --rare-maf 0.5 --min-maf 0 --finemap-l 0 -o "$OUT/sepg" --quiet > "$OUT/sepg.err" 2>&1
[ "$(gcol "$OUT/sepg.groups.tsv" skato_p_conditional all)" = "NA" ] && grep -q "admit no null model" "$OUT/sepg.err" \
  && ok "a separating signal: no conditional group test, and a warning" \
  || bad "separating signal: conditional $(gcol "$OUT/sepg.groups.tsv" skato_p_conditional all)"
"$BIN" associate --genotypes "$OUT/gr.bimbam" --samples "$OUT/gr.samples" --phenotype "$OUT/gr.pheno" --unit feature \
  --group gene -o "$OUT/grf" --quiet > "$OUT/grf.err" 2>&1
grep -q "needs the variant unit" "$OUT/grf.err" && ok "--group in the feature unit is refused" || bad "--group, feature unit: $(tail -n 1 "$OUT/grf.err")"

# ------------------------------------------------- exact tails, classes and one-record groups
# Ten carriers, three of them among 100 cases in 1,000: without covariates, Fisher's exact test.
#   R: phyper(2, 100, 900, 10, lower.tail = FALSE) = 0.0692371494 (the score with the non-carriers' terms
#   gives 0.038, half-counting the observed outcome)
awk -v F="$OUT" 'BEGIN{ printf "k3, A, B" > (F "/fx.bimbam"); for(i=1;i<=1000;i++) printf ", %d", (i<=3 || (i>100 && i<=107)) > (F "/fx.bimbam"); printf "\n" > (F "/fx.bimbam")
  for(i=1;i<=1000;i++) print "s" i > (F "/fx.samples"); print "sample\tphenotype" > (F "/fx.pheno"); for(i=1;i<=1000;i++) printf "s%d\t%d\n", i, (i<=100) > (F "/fx.pheno") }'
"$BIN" associate --genotypes "$OUT/fx.bimbam" --samples "$OUT/fx.samples" --phenotype "$OUT/fx.pheno" --min-maf 0 --finemap-l 0 -o "$OUT/fx" --quiet >/dev/null 2>&1
close "$(col "$OUT/fx.assoc.tsv" p k3)" 0.0692371494 1e-6 && ok "rare carriers without covariates: Fisher's exact p" \
  || bad "rare carriers: p $(col "$OUT/fx.assoc.tsv" p k3), Fisher 0.0692371494"
# The same carriers in other units (0.01 per copy): the same p, whichever tail is used.
awk -v F="$OUT" 'BEGIN{ for(r=1;r<=2;r++){ printf "u%d, A, B", r > (F "/fu.bimbam"); for(i=1;i<=1000;i++) printf ", %s", ((i<=3 || (i>20 && i<=34)) ? (r==1 ? "1" : "0.01") : "0") > (F "/fu.bimbam"); printf "\n" > (F "/fu.bimbam") }
  print "sample\tphenotype" > (F "/fu.pheno"); for(i=1;i<=1000;i++) printf "s%d\t%d\n", i, (i<=20) > (F "/fu.pheno") }'
"$BIN" associate --genotypes "$OUT/fu.bimbam" --samples "$OUT/fx.samples" --phenotype "$OUT/fu.pheno" --min-maf 0 --exact-carriers 0 --finemap-l 0 -o "$OUT/fu" --quiet >/dev/null 2>&1
[ -n "$(col "$OUT/fu.assoc.tsv" p u1)" ] && [ "$(col "$OUT/fu.assoc.tsv" p u1)" = "$(col "$OUT/fu.assoc.tsv" p u2)" ] \
  && [ "$(col "$OUT/fu.assoc.tsv" p_method u2)" != "score_normal_fallback" ] && ok "the saddlepoint does not depend on dosage units" \
  || bad "units: p $(col "$OUT/fu.assoc.tsv" p u1) vs $(col "$OUT/fu.assoc.tsv" p u2) ($(col "$OUT/fu.assoc.tsv" p_method u2))"
# One record, strongly associated: it is the signal.
awk -v F="$OUT" 'BEGIN{ printf "only, A, B" > (F "/f1.bimbam"); for(i=1;i<=30;i++) printf ", %d", i%3 > (F "/f1.bimbam"); printf "\n" > (F "/f1.bimbam")
  for(i=1;i<=30;i++) print "s" i > (F "/f1.samples"); print "sample\tphenotype" > (F "/f1.pheno"); for(i=1;i<=30;i++) printf "s%d\t%.3f\n", i, 2*(i%3) + ((i*7)%5)/10 > (F "/f1.pheno") }'
"$BIN" associate --genotypes "$OUT/f1.bimbam" --samples "$OUT/f1.samples" --phenotype "$OUT/f1.pheno" --unit variant --min-maf 0 --finemap-l 0 -o "$OUT/f1" --quiet >/dev/null 2>&1
[ "$(col "$OUT/f1.assoc.tsv" cond_role only)" = "signal" ] && [ "$(sv "$OUT/f1.summary.tsv" cojo_independent_signals)" = 1 ] \
  && ok "a single significant record is the signal" || bad "single record: role $(col "$OUT/f1.assoc.tsv" cond_role only)"
# A covariate of -9 is a value; an empty last cell is a missing value.
awk 'BEGIN{print "sample\tphenotype\tcv"; for(i=1;i<=12;i++) printf "s%d\t%s\t%s\n", i, (i*7)%5+0.5*i, (i==4 ? "-9" : i)}' > "$OUT/m9.pheno"
awk 'BEGIN{print "sample\tphenotype\tcv"; for(i=1;i<=12;i++) printf "s%d\t%s\t%s\n", i, (i*7)%5+0.5*i, (i==4 ? "-9.0" : i)}' > "$OUT/m9f.pheno"
awk 'BEGIN{print "sample\tphenotype\tcv"; for(i=1;i<=12;i++) printf "s%d\t%s\t%s\n", i, (i*7)%5+0.5*i, (i==4 ? "" : i)}' > "$OUT/me.pheno"
for t in m9 m9f me; do "$BIN" associate --genotypes "$OUT/lin.bimbam.gz" --samples "$OUT/lin.samples.gz" --phenotype "$OUT/$t.pheno" --min-maf 0 -o "$OUT/$t" --quiet >/dev/null 2>&1; done
[ "$(sv "$OUT/m9.summary.tsv" samples_used)" = 12 ] && [ "$(col "$OUT/m9.assoc.tsv" p f1)" = "$(col "$OUT/m9f.assoc.tsv" p f1)" ] \
  && [ "$(sv "$OUT/me.summary.tsv" samples_used)" = 11 ] && ok "a covariate of -9 is kept; an empty last cell is missing" \
  || bad "-9 covariate: samples $(sv "$OUT/m9.summary.tsv" samples_used); empty cell: $(sv "$OUT/me.summary.tsv" samples_used)"
# A dosage whose middle class carries the effect: the additive line misses it, the classes do not.
#   R: anova(lm(y ~ g), lm(y ~ factor(g))) p 9.18663143e-16; against lm(y ~ 1): 7.25815563e-15
awk -v F="$OUT" 'BEGIN{ printf "u, A, B" > (F "/gt.bimbam"); for(i=1;i<=60;i++) printf ", %d", int((i-1)/20) > (F "/gt.bimbam"); printf "\n" > (F "/gt.bimbam")
  for(i=1;i<=60;i++) print "s" i > (F "/gt.samples"); print "sample\tphenotype" > (F "/gt.pheno")
  for(i=1;i<=60;i++){ g=int((i-1)/20); printf "s%d\t%.10g\n", i, (g==1) + ((i*37)%11)/10 - 0.5 > (F "/gt.pheno") } }'
"$BIN" associate --genotypes "$OUT/gt.bimbam" --samples "$OUT/gt.samples" --phenotype "$OUT/gt.pheno" --min-maf 0 --finemap-l 0 -o "$OUT/gt" --quiet >/dev/null 2>&1
close "$(col "$OUT/gt.assoc.tsv" p_nonadditive u)" 9.18663143e-16 1e-4 && close "$(col "$OUT/gt.assoc.tsv" p_genotypic u)" 7.25815563e-15 1e-4 \
  && ok "the genotypic and non-additive tests match R's anova" \
  || bad "genotypic: $(col "$OUT/gt.assoc.tsv" p_genotypic u) / $(col "$OUT/gt.assoc.tsv" p_nonadditive u)"
# --rint: R's qnorm((rank(y) - 0.5) / n), then lm: p 9.849320197e-09.
awk -v F="$OUT" 'BEGIN{ printf "g, A, B" > (F "/ri.bimbam"); for(i=1;i<=40;i++) printf ", %d", i%3 > (F "/ri.bimbam"); printf "\n" > (F "/ri.bimbam")
  for(i=1;i<=40;i++) print "s" i > (F "/ri.samples"); print "sample\tphenotype" > (F "/ri.pheno")
  for(i=1;i<=40;i++) printf "s%d\t%.12g\n", i, exp((i%3 + ((i*13)%7)/3)/2) > (F "/ri.pheno") }'
"$BIN" associate --genotypes "$OUT/ri.bimbam" --samples "$OUT/ri.samples" --phenotype "$OUT/ri.pheno" --min-maf 0 --rint --finemap-l 0 -o "$OUT/ri" --quiet >/dev/null 2>&1
close "$(col "$OUT/ri.assoc.tsv" p g)" 9.849320197e-09 1e-5 && ok "--rint matches R's rank-normal transform" || bad "--rint p $(col "$OUT/ri.assoc.tsv" p g)"
# A rare record selected as a signal and alone in its group: nothing left to test given it (NA).
awk -v F="$OUT" 'BEGIN{ printf "rare, A, B" > (F "/gs.bimbam"); for(i=1;i<=100;i++) printf ", %d", (i<=5) > (F "/gs.bimbam"); printf "\n" > (F "/gs.bimbam")
  printf "dec, A, B" > (F "/gs.bimbam"); for(i=1;i<=100;i++) printf ", %d", i%2 > (F "/gs.bimbam"); printf "\n" > (F "/gs.bimbam")
  for(i=1;i<=100;i++) print "s" i > (F "/gs.samples"); print "sample\tphenotype" > (F "/gs.pheno")
  for(i=1;i<=100;i++) printf "s%d\t%.4f\n", i, 3*(i<=5) + ((i*37)%11)/10 > (F "/gs.pheno") }'
"$BIN" associate --genotypes "$OUT/gs.bimbam" --samples "$OUT/gs.samples" --phenotype "$OUT/gs.pheno" --unit variant --min-maf 0 \
  --group all --rare-maf 0.1 --finemap-l 0 -o "$OUT/gs" --quiet >/dev/null 2>&1
[ "$(col "$OUT/gs.assoc.tsv" cond_role rare)" = "signal" ] && [ "$(gcol "$OUT/gs.groups.tsv" skato_p_conditional all)" = "NA" ] \
  && ok "a group made only of a selected signal is untestable given it (NA)" \
  || bad "group given its own signal: $(gcol "$OUT/gs.groups.tsv" skato_p_conditional all)"
# Case/control, a one-record group without covariates: its exact tests are Fisher's, as the record's own.
"$BIN" associate --genotypes "$OUT/fx.bimbam" --samples "$OUT/fx.samples" --phenotype "$OUT/fx.pheno" --unit variant --min-maf 0 \
  --group all --rare-maf 0.1 --finemap-l 0 -o "$OUT/fxg" --quiet >/dev/null 2>&1
close "$(gcol "$OUT/fxg.groups.tsv" skat_p all)" 0.0692371494 1e-6 && close "$(gcol "$OUT/fxg.groups.tsv" skato_p all)" 0.0692371494 1e-6 \
  && ok "few carriers, case/control: exact group tests (Fisher for one record)" \
  || bad "binary group: SKAT $(gcol "$OUT/fxg.groups.tsv" skat_p all), SKAT-O $(gcol "$OUT/fxg.groups.tsv" skato_p all)"

# Fourteen samples (every one enumerable), five carriers, one of them among eight cases: still Fisher's test.
#   R: fisher.test(table(g, y))$p.value = 0.09090909091 (the whole-sample enumeration gives 0.042)
awk -v F="$OUT" 'BEGIN{ printf "c5, A, B" > (F "/f14.bimbam"); for(i=1;i<=14;i++) printf ", %d", (i<=5) > (F "/f14.bimbam"); printf "\n" > (F "/f14.bimbam")
  for(i=1;i<=14;i++) print "s" i > (F "/f14.samples"); print "sample\tphenotype" > (F "/f14.pheno")
  for(i=1;i<=14;i++) printf "s%d\t%d\n", i, (i==1 || (i>=6 && i<=12)) > (F "/f14.pheno") }'
"$BIN" associate --genotypes "$OUT/f14.bimbam" --samples "$OUT/f14.samples" --phenotype "$OUT/f14.pheno" --min-maf 0 --finemap-l 0 -o "$OUT/f14" --quiet >/dev/null 2>&1
close "$(col "$OUT/f14.assoc.tsv" p c5)" 0.09090909091 1e-6 && ok "a few samples without covariates: Fisher's exact p" \
  || bad "14 samples: p $(col "$OUT/f14.assoc.tsv" p c5), Fisher 0.0909"
# A class merged into its neighbour keeps the models nested: the dosage joins the class model.
#   R: cls <- pmin(g, 2); anova(lm(y ~ g), lm(y ~ g + factor(cls))) p 0.4659560631; against lm(y ~ 1) 2.537310399e-75
awk -v F="$OUT" 'BEGIN{ n=61; printf "m, A, B" > (F "/gm.bimbam"); for(i=1;i<=n;i++){ g=(i<=20?0:(i<=40?1:(i<=60?2:3))); G[i]=g; printf ", %d", g > (F "/gm.bimbam") } printf "\n" > (F "/gm.bimbam")
  for(i=1;i<=n;i++) print "s" i > (F "/gm.samples"); print "sample\tphenotype" > (F "/gm.pheno")
  for(i=1;i<=n;i++){ e=((i-1)%3==0?-0.1:((i-1)%3==1?0:0.1)); printf "s%d\t%.10g\n", i, 2*G[i]+e > (F "/gm.pheno") } }'
"$BIN" associate --genotypes "$OUT/gm.bimbam" --samples "$OUT/gm.samples" --phenotype "$OUT/gm.pheno" --min-maf 0 --finemap-l 0 -o "$OUT/gm" --quiet >/dev/null 2>&1
close "$(col "$OUT/gm.assoc.tsv" p_nonadditive m)" 0.4659560631 1e-4 && close "$(col "$OUT/gm.assoc.tsv" p_genotypic m)" 2.537310399e-75 1e-3 \
  && ok "merged dose classes: nested tests as R's anova" || bad "merged classes: $(col "$OUT/gm.assoc.tsv" p_genotypic m) / $(col "$OUT/gm.assoc.tsv" p_nonadditive m)"
# The exact group test does not depend on dosage units: the Fisher fixture at 1e-8 per copy.
sed 's/, 1/, 0.00000001/g' "$OUT/fx.bimbam" > "$OUT/fxs.bimbam"
"$BIN" associate --genotypes "$OUT/fxs.bimbam" --samples "$OUT/fx.samples" --phenotype "$OUT/fx.pheno" --unit variant --min-maf 0 \
  --group all --rare-maf 0.1 --finemap-l 0 -o "$OUT/fxs" --quiet >/dev/null 2>&1
close "$(gcol "$OUT/fxs.groups.tsv" skat_p all)" 0.0692371494 1e-6 && ok "exact group tests do not depend on dosage units" \
  || bad "scaled group: SKAT $(gcol "$OUT/fxs.groups.tsv" skat_p all)"
# Feature unit: a lead missing exactly where the target is carried; its indicator would be the target itself.
awk -v F="$OUT" 'BEGIN{ printf "lead, A, B" > (F "/fi.bimbam"); for(i=0;i<60;i++) printf ", %s", (i<20 ? "NA" : i%2) > (F "/fi.bimbam"); printf "\n" > (F "/fi.bimbam")
  printf "tgt, A, B" > (F "/fi.bimbam"); for(i=0;i<60;i++) printf ", %d", (i<20) > (F "/fi.bimbam"); printf "\n" > (F "/fi.bimbam")
  for(i=0;i<60;i++) print "s" i > (F "/fi.samples"); print "sample\tphenotype" > (F "/fi.pheno")
  for(i=0;i<60;i++){ l=(i<20?0.5:i%2); printf "s%d\t%.10g\n", i, 2*l+4*(i<20)+((7*i)%11)/100 > (F "/fi.pheno") } }'
"$BIN" associate --genotypes "$OUT/fi.bimbam" --samples "$OUT/fi.samples" --phenotype "$OUT/fi.pheno" --unit feature --min-maf 0 --finemap-l 0 -o "$OUT/fi" --quiet >/dev/null 2>&1
[ "$(col "$OUT/fi.assoc.tsv" cond_role tgt)" = "conditioned" ] && ok "feature unit: an indicator copying the target is left out" \
  || bad "feature indicator: role $(col "$OUT/fi.assoc.tsv" cond_role tgt)"

# With covariates, counting over the carriers holds the non-carriers fixed, which is conservative once the
# carriers are a sizeable share: 6 carriers of 100 (all cases, 4 other cases, one 0/1 covariate) take the
# saddlepoint even past |z| = 2; the same 6 carriers among 200 (3%) are counted.
for nn in 100 200; do
  awk -v F="$OUT" -v n="$nn" 'BEGIN{ printf "c6, A, B" > (F "/cv" n ".bimbam"); for(i=1;i<=n;i++) printf ", %d", (i<=6) > (F "/cv" n ".bimbam"); printf "\n" > (F "/cv" n ".bimbam")
    for(i=1;i<=n;i++) print "s" i > (F "/cv" n ".samples"); print "sample\tphenotype\tc" > (F "/cv" n ".pheno")
    for(i=1;i<=n;i++) printf "s%d\t%d\t%d\n", i, (i<=10), (i*7%3==0) > (F "/cv" n ".pheno") }'
  "$BIN" associate --genotypes "$OUT/cv$nn.bimbam" --samples "$OUT/cv$nn.samples" --phenotype "$OUT/cv$nn.pheno" --min-maf 0 \
    --finemap-l 0 -o "$OUT/cv$nn" --quiet >/dev/null 2>&1
done
[ "$(col "$OUT/cv100.assoc.tsv" p_method c6)" = "score_spa" ] && [ "$(col "$OUT/cv200.assoc.tsv" p_method c6)" = "score_exact_carriers" ] \
  && ok "with covariates, carriers above 5% of the samples take the saddlepoint, below it the count" \
  || bad "covariates, 6 carriers: $(col "$OUT/cv100.assoc.tsv" p_method c6) of 100, $(col "$OUT/cv200.assoc.tsv" p_method c6) of 200"

# One residual degree of freedom is enough for a t test: three samples, g (0,1,2), y (0,1,3).
#   R: summary(lm(c(0,1,3) ~ c(0,1,2)))  ->  beta 1.5, p 0.1210377
printf 'g, A, B, 0, 1, 2\n' > "$OUT/df1.bimbam"; printf 's1\ns2\ns3\n' > "$OUT/df1.samples"
printf 'sample\tphenotype\ns1\t0\ns2\t1\ns3\t3\n' > "$OUT/df1.pheno"
"$BIN" associate --genotypes "$OUT/df1.bimbam" --samples "$OUT/df1.samples" --phenotype "$OUT/df1.pheno" --min-maf 0 --finemap-l 0 -o "$OUT/df1" --quiet >/dev/null 2>&1
close "$(col "$OUT/df1.assoc.tsv" p g)" 0.1210377 1e-5 && ok "a regression with one residual degree of freedom is tested" \
  || bad "one residual df: p $(col "$OUT/df1.assoc.tsv" p g)"

# Feature unit: a lead constant on the target's samples changes nothing, so the target keeps its own p.
#   D = 1 for samples 1-12; N missing there, else 1 for odd samples up to 23; y = 3D + 1.5N + noise.
#   R, N on its 28 samples: p 0.000669406
awk -v F="$OUT" 'BEGIN{ printf "D, A, B" > (F "/fl.bimbam"); for(i=1;i<=40;i++) printf ", %d", (i<=12) > (F "/fl.bimbam")
  printf "\nN, A, B" > (F "/fl.bimbam"); for(i=1;i<=40;i++) printf ", %s", (i<=12 ? "NA" : ((i<=23 && i%2==1) ? 1 : 0)) > (F "/fl.bimbam"); printf "\n" > (F "/fl.bimbam")
  for(i=1;i<=40;i++) print "s" i > (F "/fl.samples"); print "sample\tphenotype" > (F "/fl.pheno")
  for(i=1;i<=40;i++){ D=(i<=12); N=(i<=12?0:((i<=23 && i%2==1)?1:0)); printf "s%d\t%.10g\n", i, 3*D+1.5*N+((7*i)%11)/4-1.2 > (F "/fl.pheno") } }'
"$BIN" associate --genotypes "$OUT/fl.bimbam" --samples "$OUT/fl.samples" --phenotype "$OUT/fl.pheno" --unit feature --min-maf 0 -o "$OUT/fl" --quiet >/dev/null 2>&1
[ "$(col "$OUT/fl.assoc.tsv" cond_role N)" = conditioned ] && close "$(col "$OUT/fl.assoc.tsv" p_conditional N)" 0.000669406 1e-5 \
  && ok "feature unit: a lead constant on the target's samples is left out of its conditional test" \
  || bad "constant lead: N $(col "$OUT/fl.assoc.tsv" cond_role N) $(col "$OUT/fl.assoc.tsv" p_conditional N)"

# A lead observed only where it is constant (1 in samples 0-59, missing in 60-99, 0 after) varies on the
# target's samples (60-199) only through its mean-filled values, which copy its missingness indicator: it is
# left out and the indicator kept, in both units.
#   R: k <- i >= 60; summary(lm(y[k] ~ T[k] + miss[k]))  ->  p 4.57314764e-06
awk -v F="$OUT" 'BEGIN{ printf "L, A, B" > (F "/lc.bimbam"); for(i=0;i<200;i++) printf ", %s", (i<60 ? 1 : (i<100 ? "NA" : 0)) > (F "/lc.bimbam")
  printf "\nT, A, B" > (F "/lc.bimbam"); for(i=0;i<200;i++) printf ", %s", (i<60 ? "NA" : ((i*7)%5<2)) > (F "/lc.bimbam"); printf "\n" > (F "/lc.bimbam")
  for(i=0;i<200;i++) print "s" i > (F "/lc.samples"); print "sample\tphenotype" > (F "/lc.pheno")
  for(i=0;i<200;i++){ t = (i<60 ? 0 : ((i*7)%5<2)); printf "s%d\t%.10g\n", i, 2*(i<60) + 0.8*(i>=60 && i<100) + 0.5*t + ((i*37)%11)/5 - 1 > (F "/lc.pheno") } }'
for u in feature variant; do
  "$BIN" associate --genotypes "$OUT/lc.bimbam" --samples "$OUT/lc.samples" --phenotype "$OUT/lc.pheno" --unit $u --min-maf 0 \
    --cojo-p 1e-3 --finemap-l 0 -o "$OUT/lc_$u" --quiet >/dev/null 2>&1
done
close "$(col "$OUT/lc_feature.assoc.tsv" p_conditional T)" 4.57314764e-06 1e-5 && close "$(col "$OUT/lc_variant.assoc.tsv" p_conditional T)" 4.57314764e-06 1e-5 \
  && ok "a lead varying only where it is missing: its indicator alone is conditioned on" \
  || bad "lead copying its indicator: feature $(col "$OUT/lc_feature.assoc.tsv" p_conditional T), variant $(col "$OUT/lc_variant.assoc.tsv" p_conditional T)"

# Case/control groups past the exact limit are calibrated by the saddlepoint: 2,000 samples, one case, a
# record of 17 carriers including it. The normal tail gives 3e-27; the group must give the record's own
# saddlepoint p (the exact value is 17/2000 = 0.0085).
awk -v F="$OUT" 'BEGIN{ printf "r1, A, B" > (F "/g1.bimbam"); for(i=1;i<=2000;i++) printf ", %d", (i<=17) > (F "/g1.bimbam"); printf "\n" > (F "/g1.bimbam")
  for(i=1;i<=2000;i++) print "s" i > (F "/g1.samples"); print "sample\tphenotype" > (F "/g1.pheno"); for(i=1;i<=2000;i++) printf "s%d\t%d\n", i, (i==1) > (F "/g1.pheno") }'
"$BIN" associate --genotypes "$OUT/g1.bimbam" --samples "$OUT/g1.samples" --phenotype "$OUT/g1.pheno" --unit variant --group all \
  --min-maf 0 --finemap-l 0 -o "$OUT/g1" --quiet >/dev/null 2>&1
gp=$(gcol "$OUT/g1.groups.tsv" burden_p all)
close "$gp" "$(col "$OUT/g1.assoc.tsv" p r1)" 1e-4 && awk -v p="$gp" 'BEGIN{exit !(p > 1e-4)}' \
  && ok "case/control group past the exact limit: the saddlepoint-calibrated p ($gp)" || bad "one-case group: burden p $gp"

# The layers' delta_r2 interval re-picks each layer's best row and scores it out of sample: on a phenotype of
# pure noise, 300 noise records in one layer against 1 in another must not be ranked (the interval spans 0).
awk -v F="$OUT" 'BEGIN{ x = 12345; n = 40
  print "feature_id\tlayer\tencoding\tbubbles\tnodes\tsvtype\tgene\tAF\tAN" > (F "/dl.annot")
  for (j = 0; j <= 300; j++) { id = (j < 300 ? "A" j : "B0"); printf "%s, X, Y", id > (F "/dl.bimbam")
    for (i = 0; i < n; i++) { g = 0; for (k = 0; k < 2; k++) { x = (x * 16807) % 2147483647; if (x / 2147483647 < 0.3) g++ } printf ", %d", g > (F "/dl.bimbam") }
    printf "\n" > (F "/dl.bimbam"); printf "%s\t%s\tdosage\t.\t.\tSNV\t.\t.\t.\n", id, (j < 300 ? "A" : "B") > (F "/dl.annot") }
  for (i = 0; i < n; i++) print "s" i > (F "/dl.samples"); print "sample\tphenotype" > (F "/dl.pheno")
  for (i = 0; i < n; i++) { x = (x * 16807) % 2147483647; u1 = x / 2147483647; x = (x * 16807) % 2147483647; u2 = x / 2147483647
    printf "s%d\t%.6f\n", i, sqrt(-2 * log(u1)) * cos(6.283185307 * u2) > (F "/dl.pheno") } }'
"$BIN" associate --genotypes "$OUT/dl.bimbam" --samples "$OUT/dl.samples" --feature-annot "$OUT/dl.annot" --phenotype "$OUT/dl.pheno" \
  --unit variant --min-maf 0 --finemap-l 0 -o "$OUT/dl" --quiet >/dev/null 2>&1
dl=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)c[$i]=i;next} NR==3{print $(c["delta_r2_low"]), $(c["delta_r2_high"])}' "$OUT/dl.layers.tsv")
awk -v a="${dl% *}" -v b="${dl#* }" 'BEGIN{exit !(a + 0 <= 0 && b + 0 >= 0 && a != "NA")}' \
  && ok "delta_r2 on noise: the interval spans 0 ($dl)" || bad "delta_r2 on noise: interval $dl"

# Conditioning columns that repeat one another. DEL is 0/1; SNV1 and SNV2 are missing exactly in DEL's
# carriers, so their missingness indicators equal DEL's dosage and each other; C is independent; R1 is a
# rare record of 3 carriers. The repeated columns are left out, so C, whose own p clears Bonferroni, is a
# signal and not collinear, and every conditional p is the ordinary least squares one on the design with
# the duplicates removed. Python numpy/scipy OLS, the conditioning records mean-filled:
#   C | DEL, SNV1, SNV2: 3.452139115e-09; DEL | SNV1, SNV2, C: 6.84871296e-30
#   SNV1 | SNV2, C on its 267 observed samples: 2.700428126e-21; R1 | DEL, SNV1, SNV2, C: 0.03066710752
awk -v F="$OUT" 'function u() { x = (x * 16807) % 2147483647; return x / 2147483647 }
BEGIN{ x = 4242; n = 400
  for (i = 0; i < n; i++) {
    del[i] = (u() < 0.3); c[i] = (u() < 0.4) + (u() < 0.4); r1[i] = (i < 3)
    if (del[i]) { s1[i] = "NA"; s2[i] = "NA"; v1 = 0.5; v2 = 0.5 }
    else { s1[i] = (u() < 0.5) + (u() < 0.5); s2[i] = (u() < 0.3) + (u() < 0.3); v1 = s1[i]; v2 = s2[i] }
    a = u(); b = u(); e = sqrt(-2 * log(a)) * cos(6.283185307 * b)
    y[i] = 2 * del[i] + 0.9 * v1 + 0.7 * v2 + 0.5 * c[i] + 1.5 * r1[i] + e }
  split("DEL SNV1 SNV2 C R1", ids, " ")
  for (k = 1; k <= 5; k++) { printf "%s, A, B", ids[k] > (F "/cc.bimbam")
    for (i = 0; i < n; i++) printf ", %s", (k == 1 ? del[i] : k == 2 ? s1[i] : k == 3 ? s2[i] : k == 4 ? c[i] : r1[i]) > (F "/cc.bimbam")
    printf "\n" > (F "/cc.bimbam") }
  print "sample\tphenotype" > (F "/cc.pheno")
  for (i = 0; i < n; i++) { print "s" i > (F "/cc.samples"); printf "s%d\t%.6f\n", i, y[i] > (F "/cc.pheno") } }'
"$BIN" associate --genotypes "$OUT/cc.bimbam" --samples "$OUT/cc.samples" --phenotype "$OUT/cc.pheno" --unit variant --group all \
  --min-maf 0 --finemap-l 0 -o "$OUT/cc" > /dev/null 2> "$OUT/cc.err"
[ "$(col "$OUT/cc.assoc.tsv" cond_role C)" = "signal" ] && [ "$(sv "$OUT/cc.summary.tsv" cojo_independent_signals)" = "4" ] \
  && ok "redundant conditioning columns: the independent record is a signal, 4 selected" \
  || bad "redundant conditioning columns: C is $(col "$OUT/cc.assoc.tsv" cond_role C), $(sv "$OUT/cc.summary.tsv" cojo_independent_signals) signals"
close "$(col "$OUT/cc.assoc.tsv" p_conditional C)" 3.452139115e-09 1e-4 && close "$(col "$OUT/cc.assoc.tsv" p_conditional DEL)" 6.84871296e-30 1e-4 \
  && close "$(col "$OUT/cc.assoc.tsv" p_conditional SNV1)" 2.700428126e-21 1e-4 && close "$(col "$OUT/cc.assoc.tsv" p_conditional R1)" 0.03066710752 1e-4 \
  && [ "$(col "$OUT/cc.assoc.tsv" cond_role R1)" = "shadow" ] \
  && ok "and the conditional p-values are OLS on the design without the repeated columns" \
  || bad "conditional p: C $(col "$OUT/cc.assoc.tsv" p_conditional C), DEL $(col "$OUT/cc.assoc.tsv" p_conditional DEL), SNV1 $(col "$OUT/cc.assoc.tsv" p_conditional SNV1), R1 $(col "$OUT/cc.assoc.tsv" p_conditional R1) ($(col "$OUT/cc.assoc.tsv" cond_role R1))"
[ "$(gcol "$OUT/cc.groups.tsv" skato_p_conditional all)" != "NA" ] && ! grep -q "admit no null model" "$OUT/cc.err" \
  && ok "the group test given the same signals has a null model" \
  || bad "group given redundant signals: skato_p_conditional $(gcol "$OUT/cc.groups.tsv" skato_p_conditional all)"
# The same with 10000 added to DEL's dosage: its column still repeats SNV1's and SNV2's indicators and is
# left out where they are in, and every conditional p is unchanged.
awk -F', ' -v OFS=', ' '$1 == "DEL" { for (i = 4; i <= NF; i++) $i += 10000 } 1' "$OUT/cc.bimbam" > "$OUT/cco.bimbam"
"$BIN" associate --genotypes "$OUT/cco.bimbam" --samples "$OUT/cc.samples" --phenotype "$OUT/cc.pheno" --unit variant \
  --min-maf 0 --finemap-l 0 -o "$OUT/cco" --quiet >/dev/null 2>&1
cco_ok=1
for id in DEL SNV1 SNV2 C R1; do
  close "$(col "$OUT/cco.assoc.tsv" p_conditional $id)" "$(col "$OUT/cc.assoc.tsv" p_conditional $id)" 1e-4 || cco_ok=0
  [ "$(col "$OUT/cco.assoc.tsv" cond_role $id)" = "$(col "$OUT/cc.assoc.tsv" cond_role $id)" ] || cco_ok=0
done
[ "$cco_ok" = 1 ] && ok "a repeated conditioning column with a constant offset is still left out (same conditional p)" \
  || bad "offset DEL: $(for id in DEL SNV1 SNV2 C R1; do printf '%s %s/%s ' $id "$(col "$OUT/cco.assoc.tsv" p_conditional $id)" "$(col "$OUT/cc.assoc.tsv" p_conditional $id)"; done)"

# ------------------------------------------------- a conditioning record that separates cases from controls
# A is common, with an effect; R has 12 carriers, all cases (Rc: 30 carriers, all controls); N0 and N1
# are null. Given R, the logistic null model has no maximum: its limit leaves the test on R's
# non-carriers, so A given R must be A's own test on those samples, and A, R both signals.
sep_mk() { # <name> <mode: case|ctrl> <with covariate: 0|1>
  awk -v F="$OUT/$1" -v mode="$2" -v cv="$3" 'function u() { x = (x * 16807) % 2147483647; return x / 2147483647 }
  BEGIN{ x = 9173; n = 800; k = 0
    for (i = 0; i < n; i++) {
      a[i] = (u() < 0.4) + (u() < 0.4); n0[i] = (u() < 0.3) + (u() < 0.3); n1[i] = (u() < 0.2); c[i] = (u() < 0.5)
      y[i] = (u() < 1 / (1 + exp(-(-1.6 + 1.1 * a[i] + (cv ? 0.4 * c[i] : 0))))) }
    for (i = 0; i < n; i++) { r[i] = 0; if (y[i] == (mode == "case") && k < (mode == "case" ? 12 : 30)) { r[i] = 1; k++ } }
    split("A R N0 N1", ids, " ")
    for (j = 1; j <= 4; j++) { printf "%s, A, B", ids[j] > (F ".bimbam")
      for (i = 0; i < n; i++) printf ", %d", (j == 1 ? a[i] : j == 2 ? r[i] : j == 3 ? n0[i] : n1[i]) > (F ".bimbam")
      printf "\n" > (F ".bimbam") }
    printf "sample\tphenotype%s\n", (cv ? "\tcv" : "") > (F ".pheno")
    for (i = 0; i < n; i++) { print "s" i > (F ".samples"); if (cv) printf "s%d\t%d\t%d\n", i, y[i], c[i] > (F ".pheno"); else printf "s%d\t%d\n", i, y[i] > (F ".pheno")
      if (!r[i]) print "s" i > (F ".rest") } }'
  # R's non-carriers alone, without R: the reference for every conditional test given R alone
  awk -F', ' -v F="$OUT/$1" 'BEGIN{ while ((getline l < (F ".bimbam")) > 0) if (l ~ /^R,/) { m = split(l, f, ", "); for (i = 4; i <= m; i++) keep[i] = (f[i] == 0) } }
    $1 != "R" { o = $1 ", A, B"; for (i = 4; i <= NF; i++) if (keep[i]) o = o ", " $i; print o }' "$OUT/$1.bimbam" > "$OUT/$1_rest.bimbam"; }
sep_mk sepc case 0; sep_mk sepk ctrl 0; sep_mk sepv case 1
for nm in sepc sepk sepv; do
  "$BIN" associate --genotypes "$OUT/$nm.bimbam" --samples "$OUT/$nm.samples" --phenotype "$OUT/$nm.pheno" --unit variant \
    --min-maf 0 --finemap-l 0 -o "$OUT/$nm" --quiet >/dev/null 2>&1
  "$BIN" associate --genotypes "$OUT/${nm}_rest.bimbam" --samples "$OUT/$nm.rest" --phenotype "$OUT/$nm.pheno" --unit variant \
    --min-maf 0 --finemap-l 0 --cojo-p 1e-300 -o "$OUT/${nm}_rest" --quiet >/dev/null 2>&1
done
for nm in sepc sepk sepv; do
  [ "$(col "$OUT/$nm.assoc.tsv" cond_role A)" = signal ] && [ "$(col "$OUT/$nm.assoc.tsv" cond_role R)" = signal ] \
    && [ "$(col "$OUT/$nm.assoc.tsv" cond_role N0)" = shadow ] && [ "$(col "$OUT/$nm.assoc.tsv" cond_role N1)" = shadow ] \
    && [ "$(col "$OUT/$nm.assoc.tsv" n_conditional A)" = "$(wc -l < "$OUT/$nm.rest" | tr -d ' ')" ] \
    && close "$(col "$OUT/$nm.assoc.tsv" p_conditional A)" "$(col "$OUT/${nm}_rest.assoc.tsv" p A)" 1e-9 \
    && ok "a separating conditioning record ($nm): A given R is A's test on R's non-carriers, both signals" \
    || bad "separating record ($nm): A $(col "$OUT/$nm.assoc.tsv" cond_role A) p_cond $(col "$OUT/$nm.assoc.tsv" p_conditional A) (on R's non-carriers $(col "$OUT/${nm}_rest.assoc.tsv" p A)), R $(col "$OUT/$nm.assoc.tsv" cond_role R), N0 $(col "$OUT/$nm.assoc.tsv" cond_role N0)"
done
# R missing in two cases and three controls: its mean-filled dosage no longer separates alone, but beside
# its missingness indicator it does; the missing samples stay (they could take any value), the carriers go.
awk -F'\t' 'NR > 1 { print $2 }' "$OUT/sepc.pheno" > "$OUT/sepc.y"
awk -F', ' -v OFS=', ' 'NR == FNR { y[FNR - 1] = $1; next }
  $1 == "R" { for (i = 104; i <= NF; i++) if ($i == 0) { if (y[i - 4] == 0 && k0 < 3) { $i = "NA"; k0++ } else if (y[i - 4] == 1 && k1 < 2) { $i = "NA"; k1++ } } } 1' \
  "$OUT/sepc.y" "$OUT/sepc.bimbam" > "$OUT/sepm.bimbam"
"$BIN" associate --genotypes "$OUT/sepm.bimbam" --samples "$OUT/sepc.samples" --phenotype "$OUT/sepc.pheno" --unit variant \
  --min-maf 0 --finemap-l 0 -o "$OUT/sepm" --quiet >/dev/null 2>&1
[ "$(col "$OUT/sepm.assoc.tsv" cond_role A)" = signal ] && [ "$(col "$OUT/sepm.assoc.tsv" n_conditional A)" = 788 ] \
  && [ "$(col "$OUT/sepm.assoc.tsv" cond_role N0)" = shadow ] \
  && ok "a separating record with missing samples: tested on all but its carriers" \
  || bad "separating record with NA: A $(col "$OUT/sepm.assoc.tsv" cond_role A) on $(col "$OUT/sepm.assoc.tsv" n_conditional A), N0 $(col "$OUT/sepm.assoc.tsv" cond_role N0)"
# The same in the feature unit, R the lead: every other row is conditioned on it, on R's non-carriers.
awk -F', ' '$1 != "A"' "$OUT/sepc.bimbam" > "$OUT/sepf.bimbam"
"$BIN" associate --genotypes "$OUT/sepf.bimbam" --samples "$OUT/sepc.samples" --phenotype "$OUT/sepc.pheno" --unit feature \
  --min-maf 0 -o "$OUT/sepf" --quiet >/dev/null 2>&1
"$BIN" associate --genotypes "$OUT/sepc_rest.bimbam" --samples "$OUT/sepc.rest" --phenotype "$OUT/sepc.pheno" --unit feature \
  --min-maf 0 -o "$OUT/sepf_rest" --quiet >/dev/null 2>&1
[ "$(col "$OUT/sepf.assoc.tsv" cond_role R)" = lead ] && [ "$(col "$OUT/sepf.assoc.tsv" cond_role N0)" = conditioned ] \
  && close "$(col "$OUT/sepf.assoc.tsv" p_conditional N0)" "$(col "$OUT/sepf_rest.assoc.tsv" p N0)" 1e-9 \
  && ok "a separating lead (feature unit): the others are conditioned, on its non-carriers" \
  || bad "separating lead: R $(col "$OUT/sepf.assoc.tsv" cond_role R), N0 $(col "$OUT/sepf.assoc.tsv" cond_role N0) p_cond $(col "$OUT/sepf.assoc.tsv" p_conditional N0), reference $(col "$OUT/sepf_rest.assoc.tsv" p N0)"
# The lead missing in five of its non-carriers: beside its indicator it separates as before, and on the samples
# left its mean-filled dosage is the indicator scaled, one column. N0 given R is then N0's test on R's
# non-carriers with the indicator as a covariate.
awk -F', ' '$1 != "A"' "$OUT/sepm.bimbam" > "$OUT/sepmf.bimbam"
"$BIN" associate --genotypes "$OUT/sepmf.bimbam" --samples "$OUT/sepc.samples" --phenotype "$OUT/sepc.pheno" --unit feature \
  --min-maf 0 -o "$OUT/sepmf" --quiet >/dev/null 2>&1
awk -F', ' '$1 == "R" { for (i = 4; i <= NF; i++) if ($i == "NA") print "s" (i - 4) }' "$OUT/sepm.bimbam" > "$OUT/sepmf.na"
awk -F'\t' -v OFS='\t' 'NR == FNR { na[$1] = 1; next } FNR == 1 { print $0, "miss"; next } { print $0, ($1 in na) ? 1 : 0 }' \
  "$OUT/sepmf.na" "$OUT/sepc.pheno" > "$OUT/sepmf_rest.pheno"
"$BIN" associate --genotypes "$OUT/sepc_rest.bimbam" --samples "$OUT/sepc.rest" --phenotype "$OUT/sepmf_rest.pheno" --unit variant \
  --min-maf 0 --finemap-l 0 --cojo-p 1e-300 -o "$OUT/sepmf_rest" --quiet >/dev/null 2>&1
[ "$(wc -l < "$OUT/sepmf.na" | tr -d ' ')" = 5 ] && [ "$(col "$OUT/sepmf.assoc.tsv" cond_role N0)" = conditioned ] \
  && [ "$(col "$OUT/sepmf.assoc.tsv" n_conditional N0)" = 788 ] \
  && close "$(col "$OUT/sepmf.assoc.tsv" p_conditional N0)" "$(col "$OUT/sepmf_rest.assoc.tsv" p N0)" 1e-6 \
  && ok "a separating lead with missing dosages (feature unit): the others are conditioned, on its non-carriers" \
  || bad "separating lead with NA: N0 $(col "$OUT/sepmf.assoc.tsv" cond_role N0) on $(col "$OUT/sepmf.assoc.tsv" n_conditional N0) p_cond $(col "$OUT/sepmf.assoc.tsv" p_conditional N0), reference $(col "$OUT/sepmf_rest.assoc.tsv" p N0)"
# Covariates on the samples a separating record leaves: R is 0 in 40 controls, 2 in 40 cases and 1 in the
# 120 others, so A given R is A's test on those 120. C is 1 on the 80 removed, constant on the rest; C2 is
# C3 rescaled there and alternates on the removed; C3 varies and has an effect. Given R, C and one of C2, C3
# add nothing: the test is A's on the 120 with C3 alone, not collinear for want of a null fit.
awk -v F="$OUT/sepcv" 'function u() { x = (x * 16807) % 2147483647; return x / 2147483647 }
  BEGIN { x = 4242; n = 200
    for (i = 0; i < n; i++) {
      if (i < 40) { r = 0; y = 0 } else if (i < 80) { r = 2; y = 1 } else r = 1
      a[i] = (u() < 0.4) + (u() < 0.4); c3[i] = (u() < 0.5)
      if (r == 1) y = (u() < 1 / (1 + exp(-(-0.8 + 1.0 * a[i] + 0.8 * c3[i]))))
      R[i] = r; Y[i] = y; c[i] = (r != 1); c2[i] = (r != 1) ? 5 * (i % 2) : 2 * c3[i] + 1 }
    printf "R, A, B" > (F ".bimbam"); for (i = 0; i < n; i++) printf ", %d", R[i] > (F ".bimbam")
    printf "\nA, A, B" > (F ".bimbam"); for (i = 0; i < n; i++) printf ", %d", a[i] > (F ".bimbam"); printf "\n" > (F ".bimbam")
    printf "A, A, B" > (F "_rest.bimbam"); for (i = 0; i < n; i++) if (R[i] == 1) printf ", %d", a[i] > (F "_rest.bimbam")
    printf "\n" > (F "_rest.bimbam")
    print "sample\tphenotype\tC\tC2\tC3" > (F ".pheno"); print "sample\tphenotype\tC3" > (F "_rest.pheno")
    for (i = 0; i < n; i++) { print "s" i > (F ".samples"); printf "s%d\t%d\t%d\t%d\t%d\n", i, Y[i], c[i], c2[i], c3[i] > (F ".pheno")
      if (R[i] == 1) { print "s" i > (F ".rest"); printf "s%d\t%d\t%d\n", i, Y[i], c3[i] > (F "_rest.pheno") } } }'
"$BIN" associate --genotypes "$OUT/sepcv.bimbam" --samples "$OUT/sepcv.samples" --phenotype "$OUT/sepcv.pheno" --unit variant \
  --min-maf 0 --finemap-l 0 -o "$OUT/sepcv" --quiet >/dev/null 2>&1
"$BIN" associate --genotypes "$OUT/sepcv_rest.bimbam" --samples "$OUT/sepcv.rest" --phenotype "$OUT/sepcv_rest.pheno" --unit variant \
  --min-maf 0 --finemap-l 0 --cojo-p 1e-300 -o "$OUT/sepcv_rest" --quiet >/dev/null 2>&1
[ "$(col "$OUT/sepcv.assoc.tsv" cond_role R)" = signal ] && [ "$(col "$OUT/sepcv.assoc.tsv" n_conditional A)" = 120 ] \
  && close "$(col "$OUT/sepcv.assoc.tsv" p_conditional A)" "$(col "$OUT/sepcv_rest.assoc.tsv" p A)" 1e-6 \
  && ok "covariates constant or collinear on the samples a separating record leaves are left out of its test" \
  || bad "covariates after separation: R $(col "$OUT/sepcv.assoc.tsv" cond_role R), A $(col "$OUT/sepcv.assoc.tsv" cond_role A) on $(col "$OUT/sepcv.assoc.tsv" n_conditional A) p_cond $(col "$OUT/sepcv.assoc.tsv" p_conditional A), reference $(col "$OUT/sepcv_rest.assoc.tsv" p A)"
# A record that separates them completely leaves nothing to test given it: the other is collinear.
"$BIN" associate --genotypes "$OUT/sep.bimbam.gz" --samples "$OUT/lin.samples.gz" --phenotype "$OUT/sep.pheno" --unit variant \
  --min-maf 0 --finemap-l 0 -o "$OUT/sepx" --quiet >/dev/null 2>&1
[ "$(col "$OUT/sepx.assoc.tsv" cond_role sep)" = signal ] && [ "$(col "$OUT/sepx.assoc.tsv" cond_role com)" = collinear ] \
  && ok "a completely separating signal: the other record is collinear" \
  || bad "complete separation: sep $(col "$OUT/sepx.assoc.tsv" cond_role sep), com $(col "$OUT/sepx.assoc.tsv" cond_role com)"

# ------------------------------------------------- conditioning columns with a large constant offset
# 40 samples; S and C agree in 36 (r^2 0.64), the phenotype is S plus noise: given S, C adds nothing
# (p_conditional 1), and S is the one signal. The same with 10000 or -10000 added to S, or S in units of
# 1e-6: conditioning depends neither on a constant added to a column nor on its scale.
off_mk() { awk -v off="$2" -v sc="$3" -v F="$OUT/$1" 'BEGIN{ split("18 2 2 18", cnt, " "); split("0 0 1 1", gg, " "); split("0 1 0 1", cc, " "); n = 0
  for (k = 1; k <= 4; k++) for (j = 0; j < cnt[k]; j++) { g[n] = gg[k]; c[n] = cc[k]; y[n] = gg[k] + (j % 2 ? -0.1 : 0.1); n++ }
  printf "S, A, B" > (F ".bimbam"); for (i = 0; i < n; i++) printf ", %.10g", off + sc * g[i] > (F ".bimbam"); printf "\n" > (F ".bimbam")
  printf "C, A, B" > (F ".bimbam"); for (i = 0; i < n; i++) printf ", %d", c[i] > (F ".bimbam"); printf "\n" > (F ".bimbam")
  print "sample\tphenotype" > (F ".pheno"); for (i = 0; i < n; i++) { print "s" i > (F ".samples"); printf "s%d\t%.1f\n", i, y[i] > (F ".pheno") } }'
  "$BIN" associate --genotypes "$OUT/$1.bimbam" --samples "$OUT/$1.samples" --phenotype "$OUT/$1.pheno" --unit variant \
    --model linear --min-maf 0 --finemap-l 0 -o "$OUT/$1" --quiet >/dev/null 2>&1; }
for v in "of0 0 1" "ofp 10000 1" "ofn -10000 1" "ofs 0 1e-6"; do
  set -- $v; off_mk "$1" "$2" "$3"
  [ "$(col "$OUT/$1.assoc.tsv" cond_role S)" = signal ] && [ "$(col "$OUT/$1.assoc.tsv" cond_role C)" = shadow ] \
    && close "$(col "$OUT/$1.assoc.tsv" p_conditional C)" 1 1e-9 && [ "$(sv "$OUT/$1.summary.tsv" cojo_independent_signals)" = 1 ] \
    && ok "S = $2 + $3 g: C given S has p 1, one signal" \
    || bad "S = $2 + $3 g: C $(col "$OUT/$1.assoc.tsv" cond_role C) p_cond $(col "$OUT/$1.assoc.tsv" p_conditional C), $(sv "$OUT/$1.summary.tsv" cojo_independent_signals) signals"
done
# Near-collinear but informative: T is S but in about 1 in 7 samples (r^2 about 0.8), and both act on the
# phenotype. Each given the other keeps its p when 10000 is added to T, or T is in units of 1e-3.
nc_mk() { awk -v off="$2" -v sc="$3" -v F="$OUT/$1" 'function u() { x = (x * 16807) % 2147483647; return x / 2147483647 }
  BEGIN{ x = 5150; n = 300
    for (i = 0; i < n; i++) { s[i] = (u() < 0.4) + (u() < 0.4); t[i] = (u() < 0.15) ? (u() < 0.4) + (u() < 0.4) : s[i]
      a = u(); b = u(); y[i] = 0.6 * s[i] + 0.6 * t[i] + sqrt(-2 * log(a)) * cos(6.283185307 * b) }
    printf "S, A, B" > (F ".bimbam"); for (i = 0; i < n; i++) printf ", %d", s[i] > (F ".bimbam"); printf "\n" > (F ".bimbam")
    printf "T, A, B" > (F ".bimbam"); for (i = 0; i < n; i++) printf ", %.10g", off + sc * t[i] > (F ".bimbam"); printf "\n" > (F ".bimbam")
    print "sample\tphenotype" > (F ".pheno"); for (i = 0; i < n; i++) { print "s" i > (F ".samples"); printf "s%d\t%.6f\n", i, y[i] > (F ".pheno") } }'
  "$BIN" associate --genotypes "$OUT/$1.bimbam" --samples "$OUT/$1.samples" --phenotype "$OUT/$1.pheno" --unit variant \
    --model linear --min-maf 0 --finemap-l 0 --cojo-p 0.05 -o "$OUT/$1" --quiet >/dev/null 2>&1; }
nc_mk nc0 0 1; nc_mk ncp 10000 1; nc_mk ncs 0 1e-3
for nm in ncp ncs; do
  close "$(col "$OUT/$nm.assoc.tsv" p_conditional S)" "$(col "$OUT/nc0.assoc.tsv" p_conditional S)" 1e-6 \
    && close "$(col "$OUT/$nm.assoc.tsv" p_conditional T)" "$(col "$OUT/nc0.assoc.tsv" p_conditional T)" 1e-6 \
    && [ "$(col "$OUT/$nm.assoc.tsv" cond_role T)" = signal ] && [ "$(col "$OUT/$nm.assoc.tsv" cond_role S)" = signal ] \
    && ok "near-collinear informative columns ($nm): the conditional p-values do not depend on T's offset or units" \
    || bad "near-collinear ($nm): S $(col "$OUT/$nm.assoc.tsv" p_conditional S) vs $(col "$OUT/nc0.assoc.tsv" p_conditional S), T $(col "$OUT/$nm.assoc.tsv" p_conditional T) vs $(col "$OUT/nc0.assoc.tsv" p_conditional T) ($(col "$OUT/$nm.assoc.tsv" cond_role T))"
done

# ------------------------------------------------- tails that must not depend on the dosage's units
# The 17-carrier record of 2,000 samples with one case, in units of 0.001 per copy: the group's saddlepoint
# calibration gives the same p as in units of 1.
awk -F', ' -v OFS=', ' '{ for (i = 4; i <= NF; i++) if ($i != 0) $i = "0.001"; print }' "$OUT/g1.bimbam" > "$OUT/g1u.bimbam"
for nm in g1 g1u; do
  "$BIN" associate --genotypes "$OUT/$nm.bimbam" --samples "$OUT/g1.samples" --phenotype "$OUT/g1.pheno" --unit variant --group all \
    --min-maf 0 --finemap-l 0 -o "$OUT/${nm}s" --quiet >/dev/null 2>&1
done
close "$(gcol "$OUT/g1us.groups.tsv" burden_p all)" "$(gcol "$OUT/g1s.groups.tsv" burden_p all)" 1e-6 \
  && close "$(gcol "$OUT/g1us.groups.tsv" acat_p all)" "$(gcol "$OUT/g1s.groups.tsv" acat_p all)" 1e-6 \
  && ok "the group's saddlepoint calibration does not depend on dosage units" \
  || bad "group units: burden $(gcol "$OUT/g1us.groups.tsv" burden_p all) vs $(gcol "$OUT/g1s.groups.tsv" burden_p all), acat $(gcol "$OUT/g1us.groups.tsv" acat_p all) vs $(gcol "$OUT/g1s.groups.tsv" acat_p all)"
# Three carriers among 200 samples with a covariate, every fitted chance 0.1: the carriers' exact tail is
# 0.1^3 = 0.001 whether a carrier's dosage is 1 or 1e-10.
awk -v F="$OUT" 'BEGIN{ for (r = 1; r <= 2; r++) { printf "t%d, A, B", r > (F "/tw.bimbam"); for (i = 0; i < 200; i++) printf ", %s", (i < 3 ? (r == 1 ? "1" : "1e-10") : "0") > (F "/tw.bimbam"); printf "\n" > (F "/tw.bimbam") }
  print "sample\tphenotype\tcv" > (F "/tw.pheno"); for (i = 0; i < 200; i++) { print "s" i > (F "/tw.samples"); printf "s%d\t%d\t%d\n", i, (i < 20), i % 2 > (F "/tw.pheno") } }'
"$BIN" associate --genotypes "$OUT/tw.bimbam" --samples "$OUT/tw.samples" --phenotype "$OUT/tw.pheno" --unit variant --min-maf 0 \
  --finemap-l 0 -o "$OUT/tw" --quiet >/dev/null 2>&1
close "$(col "$OUT/tw.assoc.tsv" p t1)" 0.001 1e-6 && close "$(col "$OUT/tw.assoc.tsv" p t2)" 0.001 1e-6 \
  && [ "$(col "$OUT/tw.assoc.tsv" p_method t2)" = score_exact_carriers ] \
  && ok "carriers of dosage 1e-10 are counted exactly, as those of dosage 1 (p 0.001)" \
  || bad "tiny carrier dosages: p $(col "$OUT/tw.assoc.tsv" p t1) and $(col "$OUT/tw.assoc.tsv" p t2) ($(col "$OUT/tw.assoc.tsv" p_method t2)), expected 0.001"

# ------------------------------------------------- one residual degree of freedom in the selection
# Three samples, dosages 0, 1, 2, phenotype 0, 1, 2.01: the regression has one residual df, its p 0.00182861
# is below Bonferroni's 0.05, and the conditional selection takes it as the signal, as the marginal test does.
printf 'v1, A, B, 0, 1, 2\n' > "$OUT/df3.bimbam"; printf 's1\ns2\ns3\n' > "$OUT/df3.samples"
printf 'sample\tphenotype\ns1\t0\ns2\t1\ns3\t2.01\n' > "$OUT/df3.pheno"
"$BIN" associate --genotypes "$OUT/df3.bimbam" --samples "$OUT/df3.samples" --phenotype "$OUT/df3.pheno" --model linear --unit variant \
  --min-maf 0 --finemap-l 0 -o "$OUT/df3" --quiet >/dev/null 2>&1
[ "$(col "$OUT/df3.assoc.tsv" cond_role v1)" = signal ] && [ "$(sv "$OUT/df3.summary.tsv" cojo_independent_signals)" = 1 ] \
  && ok "a record with one residual df can be selected" \
  || bad "one residual df: v1 is '$(col "$OUT/df3.assoc.tsv" cond_role v1)', $(sv "$OUT/df3.summary.tsv" cojo_independent_signals) signals"

# ------------------------------------------------- the groups' log counts only grouped records
# Three rare records, one with a gene: one group of one record.
awk -v F="$OUT" 'BEGIN{ split("R0 R1 R2", ids, " "); for (r = 1; r <= 3; r++) { printf "%s, A, B", ids[r] > (F "/rg.bimbam"); for (i = 0; i < 400; i++) printf ", %d", (i == r || i == r + 10) > (F "/rg.bimbam"); printf "\n" > (F "/rg.bimbam") }
  printf "feature_id\tlayer\tencoding\tbubbles\tnodes\tsvtype\tgene\tAF\tAN\n" > (F "/rg.annot")
  for (r = 1; r <= 3; r++) printf "%s\tvariant\tdosage\t%d\t.\tDEL\t%s\t.\t.\n", ids[r], r, (r == 1 ? "GENEA" : ".") > (F "/rg.annot")
  print "sample\tphenotype" > (F "/rg.pheno"); for (i = 0; i < 400; i++) { print "s" i > (F "/rg.samples"); printf "s%d\t%.2f\n", i, (i * 37) % 101 / 10 > (F "/rg.pheno") } }'
"$BIN" associate --genotypes "$OUT/rg.bimbam" --samples "$OUT/rg.samples" --phenotype "$OUT/rg.pheno" --feature-annot "$OUT/rg.annot" \
  --unit variant --group gene --min-maf 0 --finemap-l 0 -o "$OUT/rg" > /dev/null 2> "$OUT/rg.err"
grep -q "groups by gene: 1; rare records in them: 1$" "$OUT/rg.err" && ok "the groups' log counts the records in groups" \
  || bad "groups log: $(grep 'groups by' "$OUT/rg.err")"

echo
if [ "$fails" -eq 0 ]; then echo "associate_stats: all assertions passed"; exit 0; fi
echo "associate_stats: $fails assertion(s) FAILED"; exit 1
