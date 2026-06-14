#!/usr/bin/env Rscript

# Manhattan + QQ plots for the k-mer GWAS demo, comparing the presence/absence and the
# count (multiplicity) tests side by side. Consumes scripts/gwas_demo.py's <out>.assoc.tsv
# (columns: kmer, bubble_id, pos, nodes, variant, n_carriers, max_count, pa_p, count_p,
# pa_bonf, count_bonf, pa_q, count_q).
#
#   --assoc <path>   the .assoc.tsv                                  (required)
#   --out <prefix>   writes <prefix>.manhattan.{png,pdf} + <prefix>.qq.{png,pdf}
#   --title <str>    plot title (optional)

if (!requireNamespace("ggplot2", quietly = TRUE)) stop("needs ggplot2 (conda install -c conda-forge r-ggplot2)")
args <- commandArgs(trailingOnly = TRUE)
opt <- list(assoc = NULL, out = NULL, title = "k-mer GWAS", x = "ref")
i <- 1
while (i <= length(args)) {
  a <- args[[i]]; v <- function() args[[i + 1]]
  if (a == "--assoc") { opt$assoc <- v(); i <- i + 2 }
  else if (a == "--out") { opt$out <- v(); i <- i + 2 }
  else if (a == "--title") { opt$title <- v(); i <- i + 2 }
  else if (a == "--x") { opt$x <- v(); i <- i + 2 }
  else i <- i + 1
}
if (is.null(opt$assoc) || is.null(opt$out)) stop("usage: plot_gwas.R --assoc <assoc.tsv> --out <prefix> [--title T] [--x ref|nodes]")
if (!opt$x %in% c("ref", "nodes")) stop("--x must be 'ref' (reference bp) or 'nodes' (graph order)")
d <- read.delim(opt$assoc, sep = "\t", header = TRUE, check.names = FALSE)
outdir <- dirname(opt$out); if (!identical(outdir, ".") && !dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# x coordinate: reference bp (--x ref) or graph node order (--x nodes, places every variant
# incl. reference-disjoint ones). Testing is reference-free; this only sets the plot axis.
xcol <- if (opt$x == "nodes") d$node_min else d$pos
xlab <- if (opt$x == "nodes") "graph order (node id)" else "approx. reference position (bp)"

# long form: one row per (k-mer, method)
mk <- function(method, p, q) data.frame(method = method, kmer = d$kmer, bubble_id = d$bubble_id,
                                        pos = suppressWarnings(as.numeric(xcol)), variant = d$variant,
                                        p = suppressWarnings(as.numeric(p)), q = suppressWarnings(as.numeric(q)))
L <- rbind(mk("presence/absence", d$pa_p, d$pa_q), mk("count (multiplicity)", d$count_p, d$count_q))
L$method <- factor(L$method, levels = c("presence/absence", "count (multiplicity)"))

# ---- Manhattan: x = chosen coordinate (park position-less k-mers at the left), y = -log10 q
has_pos <- any(is.finite(L$pos))
if (has_pos) {
  rng <- range(L$pos[is.finite(L$pos)]); span <- max(1, diff(rng))
  L$x <- ifelse(is.finite(L$pos), L$pos, rng[1] - 0.05 * span)
} else {
  L$x <- ave(seq_len(nrow(L)), L$method, FUN = seq_along); xlab <- "k-mer (index)"
}
L$neglogq <- -log10(pmax(L$q, 1e-300))
L$sig <- L$q < 0.05
man <- ggplot2::ggplot(L, ggplot2::aes(x = x, y = neglogq)) +
  ggplot2::geom_point(ggplot2::aes(color = sig), size = 1.4, alpha = 0.8) +
  ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999") +
  ggplot2::scale_color_manual(values = c("FALSE" = "#9ecae1", "TRUE" = "#cb181d"),
                              labels = c("FALSE" = "q >= 0.05", "TRUE" = "q < 0.05"), name = NULL) +
  ggplot2::facet_wrap(~method, ncol = 1, scales = "free_y") +
  ggplot2::labs(title = paste0(opt$title, " - Manhattan"),
                subtitle = "BH q per k-mer; dashed = q 0.05. Causal CNV peak appears only in the count panel.",
                x = xlab, y = expression(-log[10](q))) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                 strip.text = ggplot2::element_text(face = "bold"), legend.position = "bottom")

# ---- QQ: observed vs expected -log10 p, per method, with genomic-inflation lambda
qqd <- do.call(rbind, lapply(levels(L$method), function(mn) {
  s <- L[L$method == mn & is.finite(L$p), ]
  p <- sort(pmax(s$p, 1e-300)); n <- length(p)
  if (n == 0) return(NULL)
  exp_p <- (seq_len(n) - 0.5) / n
  lam <- median(qchisq(1 - p, df = 1)) / qchisq(0.5, df = 1)
  data.frame(method = sprintf("%s (lambda=%.2f)", mn, lam),
             exp = -log10(exp_p), obs = -log10(p))
}))
qq <- ggplot2::ggplot(qqd, ggplot2::aes(x = exp, y = obs)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, color = "#999999", linetype = "dashed") +
  ggplot2::geom_point(size = 1.3, alpha = 0.8, color = "#2171b5") +
  ggplot2::facet_wrap(~method, ncol = 2, scales = "free") +
  ggplot2::labs(title = paste0(opt$title, " - QQ"),
                x = expression(expected~-log[10](p)), y = expression(observed~-log[10](p))) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"), strip.text = ggplot2::element_text(face = "bold"))

ggplot2::ggsave(paste0(opt$out, ".manhattan.png"), man, width = 9, height = 6, dpi = 170)
ggplot2::ggsave(paste0(opt$out, ".manhattan.pdf"), man, width = 9, height = 6)
ggplot2::ggsave(paste0(opt$out, ".qq.png"), qq, width = 9, height = 4, dpi = 170)
ggplot2::ggsave(paste0(opt$out, ".qq.pdf"), qq, width = 9, height = 4)
cat("Wrote:", paste0(opt$out, ".manhattan.png"), "and", paste0(opt$out, ".qq.png"), "\n")
