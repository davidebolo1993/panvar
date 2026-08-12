#!/usr/bin/env bash
# associate_null.sh - type-I error calibration for `panvar associate`.
#
# A p-value is only interpretable if it is uniform when nothing is going on. Two null designs:
#
#   MODE=permute (default)  shuffle the phenotype table's sample labels. Severs every genotype-phenotype
#                           link and preserves the phenotype/covariate joint distribution -- but it also
#                           severs the genotype-COVARIATE relationship, so it does not test covariate
#                           adjustment under confounding.
#   MODE=parametric         simulate the phenotype from the fitted covariate-only model
#                           (y* ~ Bernoulli(p-hat) for a binary trait, y* ~ N(Xb-hat, sigma-hat) else),
#                           leaving the real genotype-covariate structure intact. This is the design
#                           that actually tests adjustment, and the one to use for tail work.
#
# Reported:
#   type-I error at 0.05 / 0.01 / 0.001, with the uncertainty taken ACROSS replicates -- features within
#     one replicate are correlated, so a pooled binomial interval understates the spread
#   lambda_GC               -- median chi-square / 0.4549, should be ~1
#   per-feature uniformity  -- across replicates a single feature's p-values ARE independent, so a KS
#                              test per feature is valid where a pooled one is not
#   min-p 5% quantile       -- the regional maxT threshold, assuming nothing about independence
#
#   associate_null.sh <genotypes.bimbam.gz> <samples.txt.gz> <phenotype.tsv> [n_reps] [extra assoc args]
#
# Env: PANVAR_BIN, PYTHON, RSCRIPT, SEED, OUT, MODE.
set -uo pipefail

GENO="${1:?usage: associate_null.sh <genotypes> <samples> <phenotype> [n_reps] [extra...]}"
SAMPLES="${2:?}"
PHENO="${3:?}"
NREP="${4:-200}"
shift 4 2>/dev/null || shift 3
EXTRA=("$@"); [[ ${#EXTRA[@]} -eq 0 ]] && EXTRA=()

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/.." && pwd)"
BIN="${PANVAR_BIN:-$REPO/build/panvar}"
PY="${PYTHON:-python3}"
RS="${RSCRIPT:-$HOME/miniconda3/bin/Rscript}"
SEED="${SEED:-42}"
MODE="${MODE:-permute}"
OUT="${OUT:-${TMPDIR:-/tmp}/panvar_assoc_null}"
mkdir -p "$OUT"
: > "$OUT/pvals.txt"

cat > "$OUT/_permute.py" <<'PYEOF'
import random, sys
src, dst, seed = sys.argv[1], sys.argv[2], int(sys.argv[3])
rows = [l.rstrip('\n').split('\t') for l in open(src) if l.strip()]
hdr, body = rows[0], rows[1:]
# Permute the sample column only: every phenotype/covariate row keeps its internal structure and its
# NA pattern but is reassigned to a random individual, so genotype is independent of all of it.
labels = [r[0] for r in body]
random.Random(seed).shuffle(labels)
with open(dst, 'w') as o:
    o.write('\t'.join(hdr) + '\n')
    for lab, r in zip(labels, body):
        o.write('\t'.join([lab] + r[1:]) + '\n')
PYEOF

cat > "$OUT/_parametric.R" <<'REOF'
a <- commandArgs(trailingOnly = TRUE)
d <- read.delim(a[1], check.names = FALSE)
set.seed(as.integer(a[3]))
yname <- names(d)[2]
cov <- names(d)[-(1:2)]
keep <- stats::complete.cases(d)
f <- stats::as.formula(paste0("`", yname, "` ~ ",
       if (length(cov)) paste(sprintf("`%s`", cov), collapse = " + ") else "1"))
binary <- all(stats::na.omit(d[[yname]]) %in% c(0, 1))
if (binary) {
  m <- stats::glm(f, data = d[keep, , drop = FALSE], family = stats::binomial())
  d[[yname]][keep] <- stats::rbinom(sum(keep), 1, stats::fitted(m))
} else {
  m <- stats::lm(f, data = d[keep, , drop = FALSE])
  d[[yname]][keep] <- stats::fitted(m) + stats::rnorm(sum(keep), 0, stats::sigma(m))
}
write.table(d, a[2], sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
REOF

echo "== $MODE null, $NREP replicates: $(basename "$PHENO") =="
for ((r = 0; r < NREP; r++)); do
  if [ "$MODE" = "parametric" ]; then
    "$RS" "$OUT/_parametric.R" "$PHENO" "$OUT/pheno.$r.tsv" "$((SEED + r))" >/dev/null 2>&1 || continue
  else
    "$PY" "$OUT/_permute.py" "$PHENO" "$OUT/pheno.$r.tsv" "$((SEED + r))" || continue
  fi
  "$BIN" associate --genotypes "$GENO" --samples "$SAMPLES" --phenotype "$OUT/pheno.$r.tsv" \
    -o "$OUT/rep.$r" --quiet ${EXTRA[@]+"${EXTRA[@]}"} >/dev/null 2>&1 || continue
  awk -F'\t' -v R="$r" 'NR==1{for(i=1;i<=NF;i++)c[$i]=i;next} $(c["p"])!="NA" && $(c["p"])!="" \
    {print R"\t"$(c["feature_id"])"\t"$(c["p"])}' "$OUT/rep.$r.assoc.tsv" >> "$OUT/pvals.txt"
  rm -f "$OUT/rep.$r.assoc.tsv" "$OUT/rep.$r.summary.tsv" "$OUT/pheno.$r.tsv"
done

"$RS" - "$OUT/pvals.txt" <<'RAGG'
a <- commandArgs(trailingOnly = TRUE)
d <- read.table(a[1], sep = "\t", quote = "", col.names = c("rep", "feature", "p"))
d$p <- suppressWarnings(as.numeric(d$p))
d <- d[is.finite(d$p) & d$p >= 0 & d$p <= 1, ]
if (!nrow(d)) { cat("  no p-values collected\n"); quit(status = 1) }
p <- d$p; m <- length(p); nrep <- length(unique(d$rep)); nf <- length(unique(d$feature))
chi <- qchisq(p, df = 1, lower.tail = FALSE)
lam <- median(chi) / qchisq(0.5, df = 1, lower.tail = FALSE)
cat(sprintf("  %d p-values = %d replicates x %d features\n", m, nrep, nf))
# The pooled rate is unbiased, but features within one replicate are CORRELATED, so a pooled binomial
# interval understates the spread. The variation of the per-replicate rate ACROSS replicates assumes
# nothing about independence within one, and is the honest uncertainty.
for (al in c(0.05, 0.01, 0.001)) {
  per <- tapply(d$p, d$rep, function(v) mean(v < al))
  obs <- mean(per); se <- stats::sd(per) / sqrt(length(per))
  flag <- if (is.finite(se) && se > 0 && abs(obs - al) > 1.96 * se) "  <-- MISCALIBRATED" else ""
  cat(sprintf("  type-I at %-6g : %.5f  (nominal %g, across-replicate +/- %.5f, %.0f expected events)%s\n",
              al, obs, al, 1.96 * se, al * m, flag))
}
cat(sprintf("  lambda_GC        : %.4f  (expected ~1)%s\n", lam,
            if (abs(lam - 1) > 0.15) "  <-- INFLATED/DEFLATED" else ""))
ks <- sapply(split(d$p, d$feature), function(v)
  if (length(v) >= 30) suppressWarnings(ks.test(v, "punif")$p.value) else NA_real_)
ks <- ks[is.finite(ks)]
if (length(ks)) {
  bad <- sum(ks < 0.05)
  cat(sprintf("  per-feature uniformity: %d/%d reject uniform at 0.05 (expect ~%.1f)%s\n",
              bad, length(ks), 0.05 * length(ks),
              if (bad > qbinom(0.95, length(ks), 0.05)) "  <-- NOT UNIFORM" else ""))
  cat(sprintf("  worst feature KS p : %.4g\n", min(ks)))
}
# Regional family-wise threshold with no independence assumption: the empirical 5th percentile of the
# per-replicate minimum p. Directly comparable to the Bonferroni / Meff thresholds the run reports.
minp <- tapply(d$p, d$rep, min)
cat(sprintf("  min-p 5%% quantile  : %.3g  (maxT regional threshold; Bonferroni 0.05/n = %.3g)\n",
            stats::quantile(minp, 0.05), 0.05 / nf))
RAGG
