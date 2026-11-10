# simulate.R <out dir>: the datasets dev/validate_associate/validate.sh compares panvar associate on.
#
#   unrelated/   600 samples; covariates sex, pc1, pc2; 90 records in six genes: common 0/1/2 records, copy
#                numbers 0-6, records tagging a copy number, one dominant record, rare records (4-12 carriers)
#                and groups of them. A measured phenotype (pheno.q.tsv), a case/control one with about 30% cases
#                (pheno.b.tsv), and the case/control one without covariates (pheno.b0.tsv).
#   families/    150 families of four (two parents, two children; genotypes passed on by Mendel's rules), a
#                kinship matrix of the pedigree, and a measured phenotype with a family background.
# Every file is written with a fixed seed, so the comparison is reproducible.
args <- commandArgs(TRUE)
out <- if (length(args)) args[1] else "validate_associate"
set.seed(20261107)

write_bimbam <- function(G, path) {
  lines <- vapply(seq_len(nrow(G)), function(i) paste(c(rownames(G)[i], "A", "B", format(G[i, ], trim = TRUE)), collapse = ", "), "")
  writeLines(lines, path)
}
write_tsv <- function(df, path) write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE)

# ---------------------------------------------------------------- unrelated samples
dir.create(file.path(out, "unrelated"), recursive = TRUE, showWarnings = FALSE)
n <- 600
ids <- sprintf("u%03d", seq_len(n))
sex <- rbinom(n, 1, 0.5); pc1 <- rnorm(n); pc2 <- rnorm(n)
rows <- list(); gene <- character(0)
add <- function(id, g, gn) { rows[[id]] <<- g; gene[id] <<- gn }
for (j in 1:30) add(sprintf("snv%02d", j), rbinom(n, 2, runif(1, 0.05, 0.45)), sprintf("G%d", (j - 1) %% 6 + 1))
for (j in 1:6) add(sprintf("cn%d", j), pmin(6, rpois(n, 2)), sprintf("G%d", j))
cn1 <- rows[["cn1"]]
for (j in 1:6) {   # records tagging cn1: correlated but not copies of it
  t <- cn1; flip <- sample(n, round(n * (0.05 + 0.05 * j))); t[flip] <- pmin(6, rpois(length(flip), 2))
  add(sprintf("tag%d", j), t, "G1")
}
add("dom1", rbinom(n, 2, 0.3), "G2")
for (j in 1:47) {   # rare records, 4 to 12 carriers each, in the six genes
  g <- integer(n); g[sample(n, sample(4:12, 1))] <- 1L
  add(sprintf("rare%02d", j), g, sprintf("G%d", (j - 1) %% 6 + 1))
}
G <- do.call(rbind, rows); colnames(G) <- ids
rare_g3 <- colSums(G[grepl("^rare", rownames(G)) & gene[rownames(G)] == "G3", , drop = FALSE])
lin <- 0.45 * G["cn1", ] + 0.35 * G["snv05", ] + 0.6 * (G["dom1", ] > 0) + 0.8 * rare_g3 +
       0.3 * sex + 0.5 * pc1 - 0.2 * pc2
yq <- lin + rnorm(n)
liab <- -1.6 + 0.6 * (lin - mean(lin))
yb <- rbinom(n, 1, plogis(liab))
write_bimbam(G, file.path(out, "unrelated", "g.bimbam"))
writeLines(ids, file.path(out, "unrelated", "samples.txt"))
write_tsv(data.frame(feature_id = rownames(G), layer = ifelse(grepl("^cn", rownames(G)), "cn",
                     ifelse(grepl("^tag", rownames(G)), "tag", "snv")), encoding = "dosage", bubbles = ".", nodes = ".",
                     svtype = ".", gene = gene[rownames(G)], AF = ".", AN = "."),
          file.path(out, "unrelated", "annot.tsv"))
write_tsv(data.frame(sample = ids, phenotype = round(yq, 6), sex = sex, pc1 = round(pc1, 6), pc2 = round(pc2, 6)),
          file.path(out, "unrelated", "pheno.q.tsv"))
write_tsv(data.frame(sample = ids, phenotype = yb, sex = sex, pc1 = round(pc1, 6), pc2 = round(pc2, 6)),
          file.path(out, "unrelated", "pheno.b.tsv"))
write_tsv(data.frame(sample = ids, phenotype = yb), file.path(out, "unrelated", "pheno.b0.tsv"))

# ---------------------------------------------------------------- families
dir.create(file.path(out, "families"), recursive = TRUE, showWarnings = FALSE)
nf <- 150; m <- 40
fid <- rep(seq_len(nf), each = 4); role <- rep(c("p1", "p2", "c1", "c2"), nf)
fids <- sprintf("f%03d_%s", fid, role)
H <- array(0L, c(m, 4 * nf, 2))   # two haplotypes per person
for (r in seq_len(m)) {
  af <- runif(1, 0.05, 0.5)
  for (f in seq_len(nf)) {
    p <- 4 * (f - 1)
    H[r, p + 1, ] <- rbinom(2, 1, af); H[r, p + 2, ] <- rbinom(2, 1, af)
    for (c in 3:4) H[r, p + c, ] <- c(H[r, p + 1, sample(2, 1)], H[r, p + 2, sample(2, 1)])
  }
}
GF <- H[, , 1] + H[, , 2]; rownames(GF) <- sprintf("fsnv%02d", seq_len(m)); colnames(GF) <- fids
K <- diag(4 * nf)
for (f in seq_len(nf)) {
  p <- 4 * (f - 1)
  for (a in 1:4) for (b in 1:4) if (a != b && !(a <= 2 && b <= 2)) K[p + a, p + b] <- 0.5
}
sexf <- rbinom(4 * nf, 1, 0.5)
background <- as.vector(t(chol(K + 1e-9 * diag(4 * nf))) %*% rnorm(4 * nf))
yk <- 0.4 * GF["fsnv03", ] + 0.25 * GF["fsnv11", ] + 0.3 * sexf + 0.9 * background + rnorm(4 * nf)
write_bimbam(GF, file.path(out, "families", "g.bimbam"))
writeLines(fids, file.path(out, "families", "samples.txt"))
write.table(K, file.path(out, "families", "kinship.txt"), sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)
write_tsv(data.frame(sample = fids, phenotype = round(yk, 6), sex = sexf), file.path(out, "families", "pheno.tsv"))
cat("simulated:", nrow(G), "records x", n, "unrelated samples;", nrow(GF), "records x", 4 * nf, "family members\n")
