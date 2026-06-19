#!/usr/bin/env Rscript

# Faceted GWAS comparison for ONE substrate (e.g. all k-mer runs, or all graph node/edge runs),
# so the trait that shows the stronger association is obvious at a glance. Consumes one or more of
# scripts/gwas_demo.py's <out>.assoc.tsv (columns: kmer, bubble_id, pos, node_min, nodes, variant,
# n_carriers, max_count, pa_p, count_p, pa_bonf, count_bonf, pa_q, count_q). Run it once per substrate
# (the typical use is two figures: one k-mer, one graph).
#
#   --assoc <label>=<path>   a labelled assoc.tsv (label = trait, e.g. continuous / binary); repeatable
#   --out <prefix>           writes <prefix>.png and <prefix>.pdf
#   --title <str>            plot title (optional)
#   --x ref|nodes            x axis: reference bp (ref) or graph node order (nodes, default)

if (!requireNamespace("ggplot2", quietly = TRUE)) stop("needs ggplot2 (conda install -c conda-forge r-ggplot2)")
args <- commandArgs(trailingOnly = TRUE)
opt <- list(out = NULL, title = "GWAS comparison", x = "nodes")
assocs <- list()  # label -> path, in given order
i <- 1
while (i <= length(args)) {
  a <- args[[i]]; v <- function() { if (i + 1 > length(args)) stop(paste("missing value after", a)); args[[i + 1]] }
  if (a == "--assoc") { kv <- v(); eq <- regexpr("=", kv, fixed = TRUE)
    if (eq < 1) stop("--assoc expects label=path")
    assocs[[substr(kv, 1, eq - 1)]] <- substr(kv, eq + 1, nchar(kv)); i <- i + 2 }
  else if (a == "--out") { opt$out <- v(); i <- i + 2 }
  else if (a == "--title") { opt$title <- v(); i <- i + 2 }
  else if (a == "--x") { opt$x <- v(); i <- i + 2 }
  else if (a %in% c("-h", "--help")) {
    cat("Usage: plot_gwas_compare.R --assoc <label>=<assoc.tsv> [--assoc ...] --out <prefix> [--title T] [--x ref|nodes]\n"); quit(status = 0)
  } else i <- i + 1
}
if (length(assocs) == 0 || is.null(opt$out)) stop("usage: plot_gwas_compare.R --assoc <label>=<path> [...] --out <prefix>")
if (!opt$x %in% c("ref", "nodes")) stop("--x must be 'ref' or 'nodes'")
outdir <- dirname(opt$out); if (!identical(outdir, ".") && !dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# long form: one row per (trait, test, feature). test = the two association models run on the same file.
trait_levels <- names(assocs)
L <- do.call(rbind, lapply(trait_levels, function(lab) {
  d <- read.delim(assocs[[lab]], sep = "\t", header = TRUE, check.names = FALSE)
  if (nrow(d) == 0) return(NULL)
  xcol <- if (opt$x == "nodes") d$node_min else d$pos
  mk <- function(test, p, q) data.frame(trait = lab, test = test,
                                        pos = suppressWarnings(as.numeric(xcol)),
                                        p = suppressWarnings(as.numeric(p)), q = suppressWarnings(as.numeric(q)))
  rbind(mk("count (multiplicity)", d$count_p, d$count_q),
        mk("presence/absence", d$pa_p, d$pa_q))
}))
if (is.null(L) || nrow(L) == 0) stop("no association rows read")
L$trait <- factor(L$trait, levels = trait_levels)
L$test <- factor(L$test, levels = c("count (multiplicity)", "presence/absence"))

# x: chosen coordinate; park position-less features at the left of each (trait,test) panel.
xlab <- if (opt$x == "nodes") "graph order (node id)" else "approx. reference position (bp)"
if (any(is.finite(L$pos))) {
  rng <- range(L$pos[is.finite(L$pos)]); span <- max(1, diff(rng))
  L$x <- ifelse(is.finite(L$pos), L$pos, rng[1] - 0.05 * span)
} else { L$x <- ave(seq_len(nrow(L)), L$trait, L$test, FUN = seq_along); xlab <- "feature (index)" }
L$neglogq <- -log10(pmax(L$q, 1e-300))
L$sig <- !is.na(L$q) & L$q < 0.05

man <- ggplot2::ggplot(L, ggplot2::aes(x = x, y = neglogq)) +
  ggplot2::geom_point(ggplot2::aes(color = sig), size = 1.3, alpha = 0.8) +
  ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999") +
  ggplot2::scale_color_manual(values = c("FALSE" = "#9ecae1", "TRUE" = "#cb181d"),
                              labels = c("FALSE" = "q >= 0.05", "TRUE" = "q < 0.05"), name = NULL) +
  ggplot2::facet_grid(test ~ trait, scales = "free_y") +
  ggplot2::labs(title = paste0(opt$title, " - association comparison"),
                subtitle = "BH q per feature; dashed = q 0.05. Taller peak = stronger association.",
                x = xlab, y = expression(-log[10](q))) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                 strip.text = ggplot2::element_text(face = "bold"), legend.position = "bottom")

nc <- length(trait_levels)
ggplot2::ggsave(paste0(opt$out, ".png"), man, width = max(7, 3.2 * nc + 1), height = 6, dpi = 170)
ggplot2::ggsave(paste0(opt$out, ".pdf"), man, width = max(7, 3.2 * nc + 1), height = 6)
cat("Wrote:", paste0(opt$out, ".png"), "\n")
