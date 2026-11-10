# reference.R <out dir>: recompute panvar associate's results on validate.sh's runs with the reference tools,
# print one line per comparison (the largest difference over the rows compared, against its tolerance) and
# write them to <out dir>/comparison.tsv. Exits non-zero if any comparison is outside its tolerance.
suppressMessages({ library(SKAT); library(susieR); library(SPAtest) })
# The SKAT package computes its tails with Davies' method at an absolute accuracy of 1e-6, a few percent of a p
# of 1e-5; the reference uses the same package with Davies asked for 1e-12.
local({
  davies <- SKAT:::SKAT_davies
  assignInNamespace("SKAT_davies", function(q, lambda, h = rep(1, length(lambda)), delta = rep(0, length(lambda)),
                                            sigma = 0, lim = 10000, acc = 1e-4)
    davies(q, lambda, h, delta, sigma, lim = 1e6, acc = 1e-12), "SKAT")
})
out <- commandArgs(TRUE)[1]
U <- file.path(out, "unrelated"); FAM <- file.path(out, "families")

read_bimbam <- function(p) {
  f <- strsplit(readLines(p), ",\\s*")
  M <- t(vapply(f, function(x) as.numeric(x[-(1:3)]), numeric(length(f[[1]]) - 3)))
  rownames(M) <- vapply(f, `[`, "", 1); M
}
rd <- function(p) read.delim(p, na.strings = c("NA", "."), stringsAsFactors = FALSE, check.names = FALSE)
nlp <- function(p) -log10(pmax(p, 1e-300))
res <- list()
add <- function(check, reference, rows, metric, value, tol) {
  ok <- rows > 0 && is.finite(value) && value <= tol
  res[[length(res) + 1]] <<- data.frame(check = check, reference = reference, rows = rows, metric = metric,
                                        value = signif(value, 3), tolerance = tol, pass = ok)
}
dlp <- function(a, b) max(abs(nlp(a) - nlp(b)))

G <- read_bimbam(file.path(U, "g.bimbam"))
ann <- rd(file.path(U, "annot.tsv"))
PQ <- rd(file.path(U, "pheno.q.tsv")); PB <- rd(file.path(U, "pheno.b.tsv"))
X <- as.matrix(PQ[, c("sex", "pc1", "pc2")]); n <- nrow(X)
yq <- PQ$phenotype; yb <- PB$phenotype
q <- rd(file.path(U, "q.assoc.tsv")); rownames(q) <- q$feature_id
modal <- function(g) { t <- table(g); as.numeric(names(t)[which.max(t)]) }

# ---------------------------------------------------------------- linear model: R lm
ids <- q$feature_id
fits <- t(vapply(ids, function(id) summary(lm(yq ~ G[id, ] + X))$coefficients[2, c(1, 2, 4)], numeric(3)))
add("linear: beta", "R lm", length(ids), "max |diff| / se", max(abs(q$beta - fits[, 1]) / fits[, 2]), 1e-4)
add("linear: se", "R lm", length(ids), "max relative diff", max(abs(q$se - fits[, 2]) / fits[, 2]), 1e-4)
add("linear: p", "R lm", length(ids), "max |diff| of -log10 p", dlp(q$p, fits[, 3]), 1e-4)

# ---------------------------------------------------------------- dose classes: R anova (rows with no class under 5)
gen <- c(); nad <- c(); pg <- c(); pn <- c()
for (id in ids) {
  g <- G[id, ]; tb <- table(g)
  if (length(tb) < 2 || min(tb) < 5) next
  f0 <- lm(yq ~ X); fa <- lm(yq ~ g + X); fc <- lm(yq ~ factor(g) + X)
  gen <- c(gen, q[id, "p_genotypic"]); pg <- c(pg, anova(f0, fc)$`Pr(>F)`[2])
  if (length(tb) >= 3) { nad <- c(nad, q[id, "p_nonadditive"]); pn <- c(pn, anova(fa, fc)$`Pr(>F)`[2]) }
}
add("dose classes: p_genotypic", "R anova", length(gen), "max |diff| of -log10 p", dlp(gen, pg), 1e-4)
add("dose classes: p_nonadditive", "R anova", length(nad), "max |diff| of -log10 p", dlp(nad, pn), 1e-4)

# ---------------------------------------------------------------- rank-normal transform: R qnorm(rank)
qr <- rd(file.path(U, "qr.assoc.tsv")); rownames(qr) <- qr$feature_id
z <- qnorm((rank(yq) - 0.5) / n)
pr <- vapply(ids, function(id) summary(lm(z ~ G[id, ] + X))$coefficients[2, 4], 0)
add("--rint: p", "R qnorm((rank - 0.5) / n) + lm", length(ids), "max |diff| of -log10 p", dlp(qr[ids, "p"], pr), 1e-4)

# ---------------------------------------------------------------- multiple testing: R p.adjust, Li-Ji
add("Benjamini-Hochberg q", "R p.adjust(\"BH\")", nrow(q), "max relative diff",
    max(abs(q$q_bh - p.adjust(q$p, "BH")) / p.adjust(q$p, "BH")), 1e-4)
add("Bonferroni", "R p.adjust(\"bonferroni\")", nrow(q), "max relative diff",
    max(abs(q$p_bonf - p.adjust(q$p, "bonferroni")) / p.adjust(q$p, "bonferroni")), 1e-4)
lam <- eigen(cor(t(G[ids, ])), symmetric = TRUE, only.values = TRUE)$values
lam[lam < 1e-8 * max(lam)] <- 0
liji <- ceiling(sum((lam >= 1) + (lam - floor(lam))) - 1e-9)
sq <- rd(file.path(U, "q.summary.tsv")); meff <- as.numeric(sq$value[sq$key == "meff_eigen"])
add("Meff (Li-Ji)", "Li and Ji 2005, eigenvalues in R", length(ids), "|diff|", abs(meff - liji), 0)

# ---------------------------------------------------------------- logistic score test: R glm (Rao), SPAtest, exact counts
b <- rd(file.path(U, "b.assoc.tsv")); rownames(b) <- b$feature_id
g0 <- glm(yb ~ X, family = binomial); mu <- fitted(g0)
exact_carriers <- function(g, y, mu, intercept_only) {
  m <- modal(g); car <- which(g != m); w <- g[car] - m
  if (intercept_only) {   # given the number of cases: the carriers' cases, by weight class
    K <- sum(y); nn <- length(y); cls <- split(seq_along(car), w)
    obs <- sum(w * y[car]); E <- sum(w) * K / nn
    grid <- expand.grid(lapply(cls, function(k) 0:length(k)))
    val <- as.matrix(grid) %*% as.numeric(names(cls)); kk <- rowSums(grid)
    lw <- rowSums(sapply(seq_along(cls), function(j) lchoose(length(cls[[j]]), grid[[j]]))) +
          lchoose(nn - length(car), K - kk) - lchoose(nn, K)
    a <- abs(obs - E); tol <- 1e-9 * max(sum(abs(w)), 1)
    return(min(1, sum(exp(lw[abs(val - E) >= a - tol]))))
  }
  D <- sum(w * (y[car] - mu[car])); a <- abs(D)
  conf <- as.matrix(expand.grid(rep(list(0:1), length(car))))
  vals <- conf %*% w - sum(w * mu[car])
  pr <- apply(conf, 1, function(r) prod(ifelse(r == 1, mu[car], 1 - mu[car])))
  tol <- 1e-9 * max(sum(abs(w)), 1)
  min(1, sum(pr[vals >= a - tol | vals <= -a + tol]))
}
logistic_checks <- function(tab, y, cov, label) {
  nul <- if (is.null(cov)) glm(y ~ 1, family = binomial) else glm(y ~ cov, family = binomial)
  mu <- fitted(nul)
  for (meth in c("score", "score_spa", "score_exact_carriers")) {
    rows <- tab$feature_id[tab$p_method == meth]
    if (!length(rows)) next
    ref <- vapply(rows, function(id) {
      g <- G[id, ]
      if (meth == "score") {
        f1 <- if (is.null(cov)) glm(y ~ g, family = binomial) else glm(y ~ g + cov, family = binomial)
        anova(nul, f1, test = "Rao")$`Pr(>Chi)`[2]
      } else if (meth == "score_spa") {
        # SPAtest reads genotypes 0-2; the saddlepoint p does not depend on the dosage's scale.
        g2 <- if (max(g) > 2) 2 * g / max(g) else g
        ScoreTest_SPA(genos = matrix(g2, nrow = 1), pheno = y, cov = cov, method = "SPA", minmac = 1, Cutoff = 2)$p.value
      } else exact_carriers(g, y, mu, is.null(cov))
    }, 0)
    tool <- c(score = "R glm, anova(test = \"Rao\")", score_spa = "SPAtest ScoreTest_SPA(method = \"SPA\")",
              score_exact_carriers = if (is.null(cov)) "hypergeometric counts in R" else "carrier enumeration in R")[meth]
    add(sprintf("logistic%s: %s", label, meth), tool, length(rows), "max |diff| of -log10 p", dlp(tab[rows, "p"], ref),
        if (meth == "score_spa") 1e-3 else 1e-4)
  }
}
logistic_checks(b, yb, X, "")
b0 <- rd(file.path(U, "b0.assoc.tsv")); rownames(b0) <- b0$feature_id
logistic_checks(b0, yb, NULL, ", no covariates")
ok <- b$feature_id[b$effect_status == "ok" & b$p_method == "score"]
lor <- t(vapply(ok, function(id) summary(glm(yb ~ G[id, ] + X, family = binomial))$coefficients[2, 1:2], numeric(2)))
add("logistic: log_or", "R glm", length(ok), "max |diff| / se", max(abs(b[ok, "log_or"] - lor[, 1]) / lor[, 2]), 1e-3)

# ---------------------------------------------------------------- rare-record groups: SKAT package, ACAT-O in R
acat <- function(p) {
  p <- pmin(pmax(p[is.finite(p)], 1e-300), 0.999)
  t <- mean(ifelse(p < 1e-15, 1 / (p * pi), tan((0.5 - p) * pi)))
  if (t > 1e15) 1 / (t * pi) else 0.5 - atan(t) / pi
}
skat_checks <- function(prefix, y, binary) {
  gr <- rd(file.path(U, paste0(prefix, ".groups.tsv")))
  obj <- if (binary) SKAT_Null_Model(y ~ X, out_type = "D", Adjustment = FALSE) else SKAT_Null_Model(y ~ X, out_type = "C")
  nul <- glm(y ~ X, family = if (binary) binomial else gaussian); mu <- fitted(nul)
  pb <- ps <- po <- pa <- rb <- rs <- ro <- ra <- c()
  for (k in seq_len(nrow(gr))) {
    feats <- strsplit(gr$features[k], ",")[[1]]
    Z <- sapply(feats, function(id) abs(G[id, ] - modal(G[id, ])))
    mf <- colMeans(Z > 0); w <- dbeta(mf, 1, 25)
    if (binary) {   # compared where panvar uses the plain asymptotic test: beyond the exact counts, no |z| > 2
      if (gr$carriers[k] <= 16) next
      W <- mu * (1 - mu); Xc <- cbind(1, X)
      P <- function(v) v - Xc %*% solve(t(Xc) %*% (W * Xc), t(Xc) %*% (W * v))
      zs <- apply(cbind(Z %*% diag(w, ncol(Z)), Z %*% w), 2, function(v) {
        a <- P(v); sum(v * (y - mu)) / sqrt(sum(W * a^2)) })
      if (any(abs(zs) > 2)) next
    }
    sk <- function(ww, ...) SKAT(Z, obj, weights = ww, is_check_genotype = FALSE, ...)$p.value
    b1 <- sk(w, r.corr = 1); s1 <- sk(w, method = "davies"); o1 <- sk(w, method = "optimal.adj")
    one <- rep(1, length(w))
    a1 <- acat(c(b1, s1, o1, sk(one, r.corr = 1), sk(one, method = "davies"), sk(one, method = "optimal.adj")))
    pb <- c(pb, gr$burden_p[k]); ps <- c(ps, gr$skat_p[k]); po <- c(po, gr$skato_p[k]); pa <- c(pa, gr$acat_p[k])
    rb <- c(rb, b1); rs <- c(rs, s1); ro <- c(ro, o1); ra <- c(ra, a1)
  }
  lab <- if (binary) "groups, case/control" else "groups, measurement"
  add(paste0(lab, ": burden"), "SKAT package, r.corr = 1", length(pb), "max |diff| of -log10 p", if (length(pb)) dlp(pb, rb) else NA, 1e-3)
  add(paste0(lab, ": SKAT"), "SKAT package, method = \"davies\" (accuracy 1e-12)", length(ps), "max |diff| of -log10 p", if (length(ps)) dlp(ps, rs) else NA, 1e-3)
  add(paste0(lab, ": SKAT-O"), "SKAT package, method = \"optimal.adj\" (accuracy 1e-12)", length(po), "max |diff| of -log10 p", if (length(po)) dlp(po, ro) else NA, 0.01)
  add(paste0(lab, ": ACAT-O"), "ACAT (Liu et al. 2019) over the SKAT package's six p", length(pa), "max |diff| of -log10 p", if (length(pa)) dlp(pa, ra) else NA, 0.01)
}
skat_checks("q", yq, FALSE)
skat_checks("b", yb, TRUE)

# ---------------------------------------------------------------- fine-mapping: susieR on covariate-adjusted data
C <- cbind(1, X)
res_on <- function(v) qr.resid(qr(C), v)
Gr <- apply(t(G[ids, ]), 2, res_on); yr <- res_on(yq)
fit <- susie(Gr, yr, L = 10, refine = TRUE)
pip <- susie_get_pip(fit); names(pip) <- ids
# In the credible sets the two must agree. Elsewhere a weak effect can differ: panvar searches each effect's prior
# variance over its whole range and keeps the effect when any variance raises the evidence, while susieR's optimiser
# can drop it; that spreads a little probability over rows outside the sets (reported, with a loose bound).
in_cs <- ids[!is.na(q[ids, "cs"])]
add("fine-mapping: PIP in credible sets", "susieR susie(L = 10, refine = TRUE)", length(in_cs), "max |diff|",
    max(abs(q[in_cs, "pip"] - pip[in_cs])), 1e-3)
add("fine-mapping: PIP elsewhere", "susieR susie(L = 10, refine = TRUE)", length(ids) - length(in_cs), "max |diff|",
    max(abs(q[setdiff(ids, in_cs), "pip"] - pip[setdiff(ids, in_cs)])), 0.1)
cs_r <- vapply(susie_get_cs(fit, X = Gr)$cs, function(i) paste(sort(ids[i]), collapse = ","), "")
cs_p <- tapply(q$feature_id[!is.na(q$cs)], q$cs[!is.na(q$cs)], function(v) paste(sort(v), collapse = ","))
add("fine-mapping: credible sets", "susieR susie_get_cs", length(cs_r), "sets differing",
    length(setdiff(cs_r, cs_p)) + length(setdiff(cs_p, cs_r)), 0)

# ---------------------------------------------------------------- conditional selection: forward and backward, R lm
bar <- 1e-4
cond_p <- function(id, given) {
  if (!length(given)) return(summary(lm(yq ~ G[id, ] + X))$coefficients[2, 4])
  summary(lm(yq ~ G[id, ] + X + t(G[given, , drop = FALSE])))$coefficients[2, 4]
}
sel <- c(); dropped <- c()
repeat {
  cand <- setdiff(ids, c(sel, dropped))
  pc <- vapply(cand, cond_p, 0, given = sel)
  if (!length(pc) || min(pc) >= bar) break
  sel <- c(sel, names(which.min(pc)))
  repeat {   # a selected record the later ones make redundant is dropped for good
    pk <- vapply(sel, function(s) cond_p(s, setdiff(sel, s)), 0)
    if (length(sel) < 2 || max(pk) < bar) break
    worst <- names(which.max(pk)); sel <- setdiff(sel, worst); dropped <- c(dropped, worst)
  }
}
sig_p <- sort(q$feature_id[q$cond_role %in% "signal"])
add("conditional selection: signals", "forward-backward selection in R lm (GCTA-COJO procedure)", length(sel),
    "signals differing", length(setdiff(sel, sig_p)) + length(setdiff(sig_p, sel)), 0)
shadows <- q$feature_id[q$cond_role %in% "shadow"]
pcs <- vapply(shadows, cond_p, 0, given = sel)
add("conditional selection: p given the signals", "R lm with the signals as covariates", length(shadows),
    "max |diff| of -log10 p", dlp(q[shadows, "p_conditional"], pcs), 1e-4)

# ---------------------------------------------------------------- mixed model: GEMMA; null-variance model (EMMAX) in R
GF <- read_bimbam(file.path(FAM, "g.bimbam"))
PF <- rd(file.path(FAM, "pheno.tsv")); yk <- PF$phenotype; nk <- length(yk)
K <- as.matrix(read.table(file.path(FAM, "kinship.txt")))
k <- rd(file.path(FAM, "k.assoc.tsv")); rownames(k) <- k$feature_id
gm <- read.delim(file.path(FAM, "gemma.assoc.txt")); rownames(gm) <- gm$rs
both <- intersect(k$feature_id, gm$rs)
add("mixed model: beta", "GEMMA -lmm 1", length(both), "max |diff| / se", max(abs(k[both, "beta"] - gm[both, "beta"]) / gm[both, "se"]), 1e-3)
add("mixed model: p", "GEMMA -lmm 1 (p_wald)", length(both), "max |diff| of -log10 p", dlp(k[both, "p"], gm[both, "p_wald"]), 1e-3)
e <- eigen(K, symmetric = TRUE); Ut <- t(e$vectors); d <- e$values
Xk <- cbind(1, PF$sex); Xr <- Ut %*% Xk; yr <- as.vector(Ut %*% yk)
reml <- function(ld) {
  w <- 1 / (d + exp(ld)); XtWX <- t(Xr) %*% (w * Xr)
  r <- yr - Xr %*% solve(XtWX, t(Xr) %*% (w * yr)); s2 <- sum(w * r^2) / (nk - ncol(Xk))
  -0.5 * ((nk - ncol(Xk)) * log(s2) + sum(log(d + exp(ld))) + as.numeric(determinant(XtWX)$modulus))
}
lo <- log(mean(d) * 1e-5); hi <- log(mean(d) * 1e5)
grid <- seq(lo, hi, length.out = 101); best <- grid[which.max(vapply(grid, reml, 0))]
ld <- optimize(reml, c(max(lo, best - 0.25), min(hi, best + 0.25)), maximum = TRUE, tol = 1e-10)$maximum
w <- 1 / (d + exp(ld))
emmax <- vapply(both, function(id) {
  Xg <- cbind(Xr, as.vector(Ut %*% GF[id, ])); A <- t(Xg) %*% (w * Xg)
  beta <- solve(A, t(Xg) %*% (w * yr)); r <- yr - Xg %*% beta
  s2 <- sum(w * r^2) / (nk - ncol(Xg)); tt <- beta[ncol(Xg)] / sqrt(s2 * solve(A)[ncol(Xg), ncol(Xg)])
  2 * pt(-abs(tt), nk - ncol(Xg))
}, 0)
kf <- rd(file.path(FAM, "kf.assoc.tsv")); rownames(kf) <- kf$feature_id
add("mixed model, --lmm-fast: p", "EMMAX procedure in R (null REML variance, then GLS t)", length(both),
    "max |diff| of -log10 p", dlp(kf[both, "p"], emmax), 1e-3)
sk <- rd(file.path(FAM, "k.summary.tsv"))
pve <- as.numeric(sub(".*=\\s*", "", grep("pve estimate in the null model", readLines(file.path(FAM, "gemma.log.txt")), value = TRUE)))
add("mixed model: background share", "GEMMA pve in the null model", 1, "|diff|", abs(as.numeric(sk$value[sk$key == "lmm_h2"]) - pve), 1e-3)

tab <- do.call(rbind, res)
write.table(tab, file.path(out, "comparison.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
w <- c(42, 7, 26, 10, 10)
cat(sprintf("%-46s %5s  %-24s %10s %10s  %s\n", "check", "rows", "metric", "value", "tolerance", "reference"))
for (i in seq_len(nrow(tab)))
  cat(sprintf("%-46s %5d  %-24s %10s %10s  %s%s\n", tab$check[i], tab$rows[i], tab$metric[i], format(tab$value[i]),
              format(tab$tolerance[i]), tab$reference[i], if (tab$pass[i]) "" else "   << FAIL"))
nf <- sum(!tab$pass)
cat(if (nf) sprintf("\nvalidate_associate: %d comparison(s) outside tolerance\n", nf) else "\nvalidate_associate: every comparison within tolerance\n")
quit(status = if (nf) 1 else 0)
