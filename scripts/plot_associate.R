#!/usr/bin/env Rscript
# Manhattan (before/after multiple-testing correction) + QQ from `panvar associate` output.
#
#   Rscript plot_associate.R --assoc out.assoc.tsv [--summary out.summary.tsv] --out prefix [--title T]
#
# assoc.tsv columns: feature_id, layer, bubbles, nodes, n, minor_freq, beta|log_or, se, z, p, p_bonf, q_bh.
# Writes <prefix>.manhattan.png/pdf (raw -log10 p with the nominal 0.05 and region-wide Bonferroni
# threshold lines; points passing FDR<0.05 highlighted) and <prefix>.qq.png/pdf.
suppressWarnings(suppressMessages(library(ggplot2)))

args <- commandArgs(trailingOnly = TRUE)
get <- function(flag, default = NULL) {
  i <- match(flag, args); if (is.na(i) || i == length(args)) return(default); args[i + 1]
}
assoc <- get("--assoc"); out <- get("--out"); summary_path <- get("--summary"); title <- get("--title", "panvar associate")
if (is.null(assoc) || is.null(out)) stop("usage: plot_associate.R --assoc <assoc.tsv> --out <prefix> [--summary s] [--title T]")

d <- read.delim(assoc, sep = "\t", header = TRUE, check.names = FALSE)
d <- d[is.finite(d$p), ]
if (nrow(d) == 0) stop("no finite p-values in ", assoc)
n_tests <- nrow(d)

# region-wide Bonferroni threshold (0.05 / n_tests); read n_tests from summary if given (authoritative)
if (!is.null(summary_path) && file.exists(summary_path)) {
  s <- read.delim(summary_path, sep = "\t", header = TRUE, check.names = FALSE)
  v <- s$value[s$key == "features_tested"]; if (length(v) == 1) n_tests <- as.numeric(v)
}
bonf <- 0.05 / n_tests
# FDR line: the largest raw p whose q < 0.05 (so points above it pass BH); NA if none
fdr_p <- suppressWarnings(max(d$p[is.finite(d$q_bh) & d$q_bh < 0.05]))
if (!is.finite(fdr_p)) fdr_p <- NA_real_

# x = graph order (first integer in the nodes field); k-mer features (nodes ".") parked at the left
first_int <- function(x) suppressWarnings(as.numeric(sub("^[^0-9]*([0-9]+).*$", "\\1", x)))
d$x <- first_int(d$nodes)
has_x <- any(is.finite(d$x))
if (has_x) {
  rng <- range(d$x[is.finite(d$x)]); span <- max(1, diff(rng))
  d$x[!is.finite(d$x)] <- rng[1] - 0.05 * span
} else {
  d$x <- seq_len(nrow(d))
}
d$logp <- -log10(pmax(d$p, 1e-300))
d$sig <- ifelse(d$p < bonf, "Bonferroni",
                ifelse(is.finite(d$q_bh) & d$q_bh < 0.05, "FDR<0.05", "ns"))
d$sig <- factor(d$sig, levels = c("ns", "FDR<0.05", "Bonferroni"))

cols <- c("ns" = "grey70", "FDR<0.05" = "#2c7fb8", "Bonferroni" = "#d95f02")
p_man <- ggplot(d, aes(x, logp, colour = sig)) +
  geom_point(size = 1.4, alpha = 0.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", colour = "grey50") +
  geom_hline(yintercept = -log10(bonf), linetype = "dashed", colour = "#d95f02") +
  { if (is.finite(fdr_p)) geom_hline(yintercept = -log10(fdr_p), linetype = "dashed", colour = "#2c7fb8") } +
  scale_colour_manual(values = cols, name = NULL) +
  labs(title = title, subtitle = sprintf("n=%d tests; Bonferroni 0.05/n = %.2g (dashed orange); nominal 0.05 (dotted)", n_tests, bonf),
       x = "graph order (node id; k-mers parked left)", y = expression(-log[10](p))) +
  theme_bw(base_size = 12) + theme(legend.position = "top")
ggsave(paste0(out, ".manhattan.png"), p_man, width = 10, height = 4.5, dpi = 150)
ggsave(paste0(out, ".manhattan.pdf"), p_man, width = 10, height = 4.5)

# QQ with genomic-inflation lambda
po <- sort(pmax(d$p, 1e-300)); m <- length(po)
qq <- data.frame(exp = -log10(ppoints(m)), obs = -log10(po))
lambda <- median(qchisq(1 - po, df = 1)) / qchisq(0.5, df = 1)
p_qq <- ggplot(qq, aes(exp, obs)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey60") +
  geom_point(size = 1.2, alpha = 0.8, colour = "#2c7fb8") +
  labs(title = paste0(title, " - QQ"), subtitle = sprintf("genomic inflation lambda = %.3f", lambda),
       x = expression(expected~-log[10](p)), y = expression(observed~-log[10](p))) +
  theme_bw(base_size = 12)
ggsave(paste0(out, ".qq.png"), p_qq, width = 5, height = 5, dpi = 150)
ggsave(paste0(out, ".qq.pdf"), p_qq, width = 5, height = 5)
cat("Wrote:", paste0(out, ".manhattan.png"), "and", paste0(out, ".qq.png"), "\n")
