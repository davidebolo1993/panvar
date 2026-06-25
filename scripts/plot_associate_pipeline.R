#!/usr/bin/env Rscript
# Pipeline Manhattan for `panvar associate`: one facet per processing stage, colouring which features
# survive each stage, so the funnel TEST -> FILTER MAF -> [CLUMP] -> CORRECT -> CONDITION is explicit.
# The variant tier has a CLUMP stage (genotype-r^2 LD-clumping); the graph/k-mer feature tiers do not.
# See docs/gwas/example.md and docs/algorithms/associate.md.
#   plot_associate_pipeline.R --assoc <real.assoc.tsv> --unfiltered <minmaf0.assoc.tsv> \
#       --summary <real.summary.tsv> --min-maf <X> --out <prefix> [--title T]
# --assoc       : the real run (correct Meff/clump/COJO; post-MAF features)        [required]
# --unfiltered  : a `--min-maf 0` run of the SAME data (to show TEST + FILTER MAF) [optional]
suppressWarnings(suppressMessages(library(ggplot2)))

args <- commandArgs(trailingOnly = TRUE)
get <- function(f, d = NULL) { i <- match(f, args); if (is.na(i) || i == length(args)) return(d); args[i + 1] }
assoc <- get("--assoc"); unfilt <- get("--unfiltered"); summ <- get("--summary")
out <- get("--out"); title <- get("--title", "panvar associate")
min_maf <- as.numeric(get("--min-maf", "0.01"))
W <- as.numeric(get("--width", "9")); H <- as.numeric(get("--height", "11.5")); dpi <- as.numeric(get("--dpi", "150"))
if (is.null(assoc) || is.null(out)) stop("usage: plot_associate_pipeline.R --assoc <a> --out <p> [--unfiltered u --summary s --min-maf X --title T]")

B <- read.delim(assoc, sep = "\t", header = TRUE, check.names = FALSE)
B <- B[is.finite(B$p), ]
if (nrow(B) == 0) stop("no finite p in ", assoc)

meff <- nrow(B); unit <- "feature"
if (!is.null(summ) && file.exists(summ)) {
  s <- read.delim(summ, sep = "\t", header = TRUE, check.names = FALSE)
  v <- s$value[s$key == "meff"]; if (length(v) == 1) meff <- as.numeric(v)
  u <- s$value[s$key == "unit"]; if (length(u) == 1) unit <- as.character(u)
}
is_variant <- unit == "variant"
bonf_meff <- 0.05 / meff
fdr_p <- suppressWarnings(max(B$p[is.finite(B$q_bh) & B$q_bh < 0.05])); if (!is.finite(fdr_p)) fdr_p <- NA_real_

first_int <- function(x) suppressWarnings(as.numeric(sub("^[^0-9]*([0-9]+).*$", "\\1", x)))
xof <- function(df) { v <- first_int(df$nodes); if (all(!is.finite(v))) v <- seq_len(nrow(df)); v[!is.finite(v)] <- min(v[is.finite(v)], na.rm = TRUE); v }
B$x <- xof(B); B$logp <- -log10(pmax(B$p, 1e-300))

A <- NULL
if (!is.null(unfilt) && file.exists(unfilt)) {
  A <- read.delim(unfilt, sep = "\t", header = TRUE, check.names = FALSE)
  A <- A[is.finite(A$p), ]; A$x <- xof(A); A$logp <- -log10(pmax(A$p, 1e-300))
}

role_of <- function(df) if ("cond_role" %in% names(df)) as.character(df$cond_role) else rep(".", nrow(df))
# conditional y: anchor (signal/lead) shown at its marginal height, others at p_conditional (NA dropped)
cond_y <- function(df) {
  anchor <- role_of(df) %in% c("signal", "lead")
  ifelse(anchor, df$logp, ifelse(is.finite(df$p_conditional), -log10(pmax(df$p_conditional, 1e-300)), NA_real_))
}

mk <- function(x, y, cat, stage) data.frame(x = x, y = y, cat = cat, stage = stage, stringsAsFactors = FALSE)
parts <- list()
src_test <- if (!is.null(A)) A else B
parts[["TEST"]] <- mk(src_test$x, src_test$logp, "tested", "TEST")
if (!is.null(A)) {
  pass <- A$minor_freq >= min_maf
  parts[["FILTER MAF"]] <- mk(A$x, A$logp, ifelse(pass, "kept", "dropped (MAF)"), "FILTER MAF")
}
if (is_variant) {
  lead <- if ("is_lead" %in% names(B)) B$is_lead else rep(0, nrow(B))
  parts[["CLUMP"]] <- mk(B$x, B$logp, ifelse(lead == 1, "clump lead", "shadow (LD)"), "CLUMP")
}
corr_cat <- ifelse(B$p < bonf_meff, "Bonferroni (Meff)",
            ifelse(is.finite(B$q_bh) & B$q_bh < 0.05, "BH-FDR", "ns"))
parts[["CORRECT"]] <- mk(B$x, B$logp, corr_cat, "CORRECT")
cy <- cond_y(B); role <- role_of(B)
cond_cat <- ifelse(role %in% c("signal", "lead"), "signal",
            ifelse(role == "collinear", "collinear", "conditioned (collapsed)"))
keepc <- is.finite(cy)
parts[["CONDITION"]] <- mk(B$x[keepc], cy[keepc], cond_cat[keepc], "CONDITION")

stage_levels <- c("TEST", if (!is.null(A)) "FILTER MAF", if (is_variant) "CLUMP", "CORRECT", "CONDITION")
facet_labels <- c(
  "TEST"        = "TEST - g ~ phenotype + covariates -> p",
  "FILTER MAF"  = "FILTER MAF - keep minor freq >= min-maf",
  "CLUMP"       = "CLUMP - LD r2 > 0.8 -> leads vs shadows (sets Meff)",
  "CORRECT"     = "CORRECT - BH-FDR (blue) + Bonferroni x Meff (orange)",
  "CONDITION"   = "CONDITION - COJO conditional p (shadows collapse)")
dat <- do.call(rbind, parts)
dat$stage <- factor(dat$stage, levels = stage_levels, labels = facet_labels[stage_levels])

thr <- data.frame()
addln <- function(stage, y, col, lty) data.frame(stage = facet_labels[[stage]], y = y, col = col, lty = lty, stringsAsFactors = FALSE)
thr <- rbind(thr, addln("CORRECT", -log10(bonf_meff), "#d95f02", "dashed"))
if (is.finite(fdr_p)) thr <- rbind(thr, addln("CORRECT", -log10(fdr_p), "#2c7fb8", "dashed"))
thr <- rbind(thr, addln("CONDITION", -log10(bonf_meff), "#d95f02", "dashed"))
thr$stage <- factor(thr$stage, levels = facet_labels[stage_levels])

cols <- c("tested" = "grey55", "kept" = "#2c7fb8", "dropped (MAF)" = "grey80",
          "clump lead" = "#1b9e77", "shadow (LD)" = "grey70",
          "Bonferroni (Meff)" = "#d95f02", "BH-FDR" = "#2c7fb8", "ns" = "grey78",
          "signal" = "#e7298a", "conditioned (collapsed)" = "grey70", "collinear" = "grey85")

lab <- NULL
if ("gene" %in% names(B)) {
  gs <- B[role %in% c("signal", "lead") & !is.na(B$gene) & B$gene != ".", ]
  if (nrow(gs) > 0) {
    gs <- gs[order(gs$p), ]; gs <- gs[!duplicated(gs$gene), ]
    lab <- data.frame(x = gs$x, y = -log10(pmax(gs$p, 1e-300)), gene = gs$gene,
                      stage = facet_labels[["CONDITION"]], stringsAsFactors = FALSE)
    lab$stage <- factor(lab$stage, levels = facet_labels[stage_levels])
  }
}

p <- ggplot(dat, aes(x, y, colour = cat)) +
  geom_hline(data = thr, aes(yintercept = y), colour = thr$col, linetype = thr$lty) +
  geom_point(size = 1.2, alpha = 0.85) +
  facet_wrap(~stage, ncol = 1) +
  scale_colour_manual(values = cols, name = NULL, drop = TRUE) +
  labs(title = title,
       subtitle = sprintf("%s tier - %d tested - Meff = %d - Bonferroni x Meff = %.2g", unit, nrow(B), meff, bonf_meff),
       x = "graph order (node id)", y = expression(-log[10](p)~~"/"~~-log[10](p[conditional]))) +
  theme_bw(base_size = 11) + theme(legend.position = "top")
if (!is.null(lab)) {
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(data = lab, aes(x, y, label = gene), inherit.aes = FALSE,
      size = 3, colour = "black", min.segment.length = 0, box.padding = 0.4, max.overlaps = 20)
  } else {
    p <- p + geom_text(data = lab, aes(x, y, label = gene), inherit.aes = FALSE, size = 3, vjust = -0.6, colour = "black")
  }
}
ggsave(paste0(out, ".pipeline.png"), p, width = W, height = H, dpi = dpi)
ggsave(paste0(out, ".pipeline.pdf"), p, width = W, height = H)
cat("Wrote:", paste0(out, ".pipeline.png"), "\n")
