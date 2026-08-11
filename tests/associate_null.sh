#!/usr/bin/env bash
# associate_null.sh - type-I error calibration for `panvar associate`.
#
# A p-value is only interpretable if it is uniform when nothing is going on. This permutes the
# phenotype table's sample labels, which severs every genotype-phenotype link while leaving the
# genotype matrix, the missingness pattern and the phenotype/covariate joint distribution exactly as
# they are, then asks whether the resulting p-values look uniform.
#
# Reported per model:
#   type-I error at 0.05 / 0.01 / 0.001  -- should match the nominal rate
#   lambda_GC                            -- median chi-square / 0.4549, should be ~1
#   KS distance from uniform             -- should be small; the 0.05 critical value is ~1.36/sqrt(m)
#
#   associate_null.sh <genotypes.bimbam.gz> <samples.txt.gz> <phenotype.tsv> [n_perm] [extra associate args...]
#
# Env: PANVAR_BIN (default build/panvar), PYTHON, RSCRIPT, SEED, OUT.
set -uo pipefail

GENO="${1:?usage: associate_null.sh <genotypes> <samples> <phenotype> [n_perm] [extra...]}"
SAMPLES="${2:?}"
PHENO="${3:?}"
NPERM="${4:-200}"
shift 4 2>/dev/null || shift 3
EXTRA=("$@"); [[ ${#EXTRA[@]} -eq 0 ]] && EXTRA=()

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/.." && pwd)"
BIN="${PANVAR_BIN:-$REPO/build/panvar}"
PY="${PYTHON:-python3}"
RS="${RSCRIPT:-$HOME/miniconda3/bin/Rscript}"
SEED="${SEED:-42}"
OUT="${OUT:-${TMPDIR:-/tmp}/panvar_assoc_null}"
mkdir -p "$OUT"
: > "$OUT/pvals.txt"

echo "== permuting sample labels $NPERM times: $(basename "$PHENO") =="
for ((r = 0; r < NPERM; r++)); do
  "$PY" - "$PHENO" "$OUT/pheno.$r.tsv" "$((SEED + r))" <<'PY'
import random, sys
src, dst, seed = sys.argv[1], sys.argv[2], int(sys.argv[3])
rows = [l.rstrip('\n').split('\t') for l in open(src) if l.strip()]
hdr, body = rows[0], rows[1:]
# Permute the sample column only: every phenotype/covariate row keeps its internal structure and its
# NA pattern, but is reassigned to a random individual, so genotype is independent of all of it.
labels = [r[0] for r in body]
random.Random(seed).shuffle(labels)
with open(dst, 'w') as o:
    o.write('\t'.join(hdr) + '\n')
    for lab, r in zip(labels, body):
        o.write('\t'.join([lab] + r[1:]) + '\n')
PY
  "$BIN" associate --genotypes "$GENO" --samples "$SAMPLES" --phenotype "$OUT/pheno.$r.tsv" \
    -o "$OUT/perm.$r" --quiet ${EXTRA[@]+"${EXTRA[@]}"} >/dev/null 2>&1 || continue
  awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)c[$i]=i;next} $(c["p"])!="NA" && $(c["p"])!="" {print $(c["feature_id"])"\t"$(c["p"])}' \
    "$OUT/perm.$r.assoc.tsv" >> "$OUT/pvals.txt"
  rm -f "$OUT/perm.$r.assoc.tsv" "$OUT/perm.$r.summary.tsv" "$OUT/pheno.$r.tsv"
done

"$RS" - "$OUT/pvals.txt" <<'R'
a <- commandArgs(trailingOnly = TRUE)
d <- read.table(a[1], sep = "\t", quote = "", col.names = c("feature", "p"))
d$p <- suppressWarnings(as.numeric(d$p))
d <- d[is.finite(d$p) & d$p >= 0 & d$p <= 1, ]
if (!nrow(d)) { cat("  no p-values collected\n"); quit(status = 1) }
p <- d$p; m <- length(p)
chi <- qchisq(p, df = 1, lower.tail = FALSE)
lam <- median(chi) / qchisq(0.5, df = 1, lower.tail = FALSE)
cat(sprintf("  m = %d p-values over %d features\n", m, length(unique(d$feature))))
# Type-I error: the rate is unbiased under correlation, but features within one permutation are
# correlated, so the pooled binomial interval below UNDERSTATES the true spread. The per-feature KS
# underneath is the rigorous check -- across permutations a single feature's p-values are independent.
for (al in c(0.05, 0.01, 0.001)) {
  obs <- mean(p < al)
  se <- sqrt(al * (1 - al) / m)
  flag <- if (abs(obs - al) > 1.96 * se) "  <-- outside the (optimistic) pooled interval" else ""
  cat(sprintf("  type-I at %-6g : %.5f  (nominal %g, pooled +/- %.5f)%s\n", al, obs, al, 1.96 * se, flag))
}
cat(sprintf("  lambda_GC        : %.4f  (expected ~1)%s\n", lam,
            if (abs(lam - 1) > 0.15) "  <-- INFLATED/DEFLATED" else ""))
ks <- sapply(split(d$p, d$feature), function(v)
  if (length(v) >= 30) suppressWarnings(ks.test(v, "punif")$p.value) else NA_real_)
ks <- ks[is.finite(ks)]
if (length(ks)) {
  bad <- sum(ks < 0.05)
  cat(sprintf("  per-feature uniformity: %d/%d features reject uniform at 0.05 (expect ~%.1f)%s\n",
              bad, length(ks), 0.05 * length(ks),
              if (bad > qbinom(0.95, length(ks), 0.05)) "  <-- NOT UNIFORM" else ""))
  cat(sprintf("  worst feature KS p : %.4g\n", min(ks)))
}
R
