#!/usr/bin/env Rscript
# Manhattan (before/after multiple-testing correction) + QQ from `panvar associate` output.
#
#   Rscript plot_associate.R --assoc out.assoc.tsv [--summary out.summary.tsv] --out prefix [--title T]
#
# assoc.tsv columns: feature_id, layer, bubbles, nodes, n, minor_freq, beta|log_or, se, z, p, p_bonf, q_bh, gene.
# Writes <prefix>.manhattan.png/pdf -- two stacked panels: BEFORE correction (raw -log10 p, nominal +
# Bonferroni lines) and AFTER correction (BH -log10 q, q=0.05 line). x = node id (graph) or per-k-mer
# index ordered by node id (k-mers); genes flagged with ggrepel labels from the `gene` column.
suppressWarnings(suppressMessages(library(ggplot2)))

args <- commandArgs(trailingOnly = TRUE)
get <- function(flag, default = NULL) {
  i <- match(flag, args); if (is.na(i) || i == length(args)) return(default); args[i + 1]
}
assoc <- get("--assoc"); out <- get("--out"); summary_path <- get("--summary"); title <- get("--title", "panvar associate")
man_w <- as.numeric(get("--width", "10")); man_h <- as.numeric(get("--height", "7")); dpi <- as.numeric(get("--dpi", "150"))
if (is.null(assoc) || is.null(out)) stop("usage: plot_associate.R --assoc <assoc.tsv> --out <prefix> [--summary s] [--title T] [--width W --height H --dpi D]")

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

# x-axis. Graph features ARE node/edge ids -> use the node id (genomic order). k-mers are many-per-node
# (and several can share a node), so use a per-k-mer INDEX ordered by node id: each k-mer is its own
# column, grouped by locus, so you can see how many k-mers pop up. The kmer<->node(<->gene) link now
# lives in feature_annot / the assoc `nodes` & `gene` columns (so this stays traceable).
first_int <- function(x) suppressWarnings(as.numeric(sub("^[^0-9]*([0-9]+).*$", "\\1", x)))
node1 <- first_int(d$nodes)
layer <- if ("layer" %in% names(d)) d$layer else rep(".", nrow(d))
is_kmer <- sum(layer == "kmer") > sum(layer == "graph")
if (is_kmer) {
  ord <- order(ifelse(is.finite(node1), node1, max(node1[is.finite(node1)], 0) + 1), d$p)
  d$x <- NA_real_; d$x[ord] <- seq_len(nrow(d))   # per-k-mer index, ordered by node id then p
  xlab <- "k-mer index (ordered by node id)"
} else {
  d$x <- node1
  if (any(is.finite(d$x))) {
    rng <- range(d$x[is.finite(d$x)]); d$x[!is.finite(d$x)] <- rng[1] - 0.05 * max(1, diff(rng))
  } else d$x <- seq_len(nrow(d))
  xlab <- "graph order (node id)"
}
d$logp <- -log10(pmax(d$p, 1e-300))
d$logq <- -log10(pmax(ifelse(is.finite(d$q_bh), d$q_bh, 1), 1e-300))
d$sig <- ifelse(d$p < bonf, "Bonferroni",
                ifelse(is.finite(d$q_bh) & d$q_bh < 0.05, "FDR<0.05", "ns"))
d$sig <- factor(d$sig, levels = c("ns", "FDR<0.05", "Bonferroni"))

# Two stacked panels on the SAME run: BEFORE correction (raw -log10 p, with the nominal 0.05 and
# region-wide Bonferroni lines -- every nominally-significant node shows here, LD/structure noise
# included) and AFTER correction (Benjamini-Hochberg -log10 q, with the q=0.05 line -- only the
# features that survive multiple-testing stay up; the rest collapse toward 0). Point colour is the
# overall verdict in both panels, so you can see which 'before' peaks are real vs noise.
lv <- c("before correction: raw -log10(p)", "after correction: BH -log10(q)")
long <- rbind(
  data.frame(x = d$x, y = d$logp, sig = d$sig, panel = lv[1]),
  data.frame(x = d$x, y = d$logq, sig = d$sig, panel = lv[2]))
long$panel <- factor(long$panel, levels = lv)
thr <- data.frame(
  panel = factor(c(lv[1], lv[1], lv[2]), levels = lv),
  yint  = c(-log10(0.05), -log10(bonf), -log10(0.05)),
  col   = c("grey50", "#d95f02", "#2c7fb8"),
  lty   = c("dotted", "dashed", "dashed"))

# Gene flags: one label per gene, placed at its most-significant feature, in BOTH panels (at its raw
# -log10 p height in the "before" panel and its -log10 q height in the "after" panel). Only genes whose
# top feature SURVIVES correction (Bonferroni or BH FDR<0.05) are labelled, so near-nominal neighbours
# aren't flagged as noise. Needs the `gene` column (`associate --node-genes call.node_genes.tsv`).
lab <- NULL
if ("gene" %in% names(d)) {
  g <- d[!is.na(d$gene) & d$gene != "." & (d$p < bonf | (is.finite(d$q_bh) & d$q_bh < 0.05)), ]
  if (nrow(g) > 0) {
    g <- g[order(g$p), ]; g <- g[!duplicated(g$gene), ]   # best feature per gene
    lab <- rbind(
      data.frame(x = g$x, y = -log10(pmax(g$p, 1e-300)), gene = g$gene, panel = lv[1]),
      data.frame(x = g$x, y = -log10(pmax(ifelse(is.finite(g$q_bh), g$q_bh, 1), 1e-300)),
                 gene = g$gene, panel = lv[2]))
    lab$panel <- factor(lab$panel, levels = lv)
  }
}

cols <- c("ns" = "grey70", "FDR<0.05" = "#2c7fb8", "Bonferroni" = "#d95f02")
p_man <- ggplot(long, aes(x, y, colour = sig)) +
  geom_hline(data = thr, aes(yintercept = yint), colour = thr$col, linetype = thr$lty) +
  geom_point(size = 1.3, alpha = 0.8) +
  facet_wrap(~panel, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = cols, name = NULL, drop = FALSE) +
  labs(title = title,
       subtitle = sprintf("n=%d tests; Bonferroni 0.05/n = %.2g (dashed orange); nominal/FDR 0.05 (dotted grey / dashed blue)", n_tests, bonf),
       x = xlab, y = NULL) +
  theme_bw(base_size = 12) + theme(legend.position = "top")
if (!is.null(lab)) {
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p_man <- p_man + ggrepel::geom_text_repel(data = lab, aes(x, y, label = gene),
      inherit.aes = FALSE, size = 3, min.segment.length = 0, max.overlaps = 20,
      box.padding = 0.4, colour = "black")
  } else {
    p_man <- p_man + geom_text(data = lab, aes(x, y, label = gene), inherit.aes = FALSE,
      size = 3, vjust = -0.6, colour = "black")
  }
}
ggsave(paste0(out, ".manhattan.png"), p_man, width = man_w, height = man_h, dpi = dpi)
ggsave(paste0(out, ".manhattan.pdf"), p_man, width = man_w, height = man_h)

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
ggsave(paste0(out, ".qq.png"), p_qq, width = 5, height = 5, dpi = dpi)
ggsave(paste0(out, ".qq.pdf"), p_qq, width = 5, height = 5)
cat("Wrote:", paste0(out, ".manhattan.png"), "and", paste0(out, ".qq.png"), "\n")
