#!/usr/bin/env Rscript
# Copy-number concordance (called vs ground-truth BED), one facet per gene, from
# compare_copy_number.py --dump-tsv (cols: gene sample truth_cn called_cn baseline_offset is_reference).
#   Rscript plot_cn_correlation.R --table cn_table.tsv --out results/cn_correlation
# Writes <out>.loci.png (lower-case locus totals, 2x2) and <out>.genes.png (upper-case CYP2D6/CYP2D7 split).
suppressWarnings(suppressMessages({
  library(ggplot2)
}))

args <- commandArgs(trailingOnly = TRUE)
get <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) return(default)
  args[i + 1]
}
table_path <- get("--table")
out_prefix <- get("--out", "cn_correlation")
if (is.null(table_path)) stop("need --table <tsv>")
# --dpi applies to both plots; --width/--height override the per-plot defaults when given.
dpi <- as.numeric(get("--dpi", "150"))
w_arg <- get("--width"); h_arg <- get("--height")
wnum <- function(x, d) if (is.null(x)) d else as.numeric(x)

d <- read.delim(table_path, stringsAsFactors = FALSE)
d <- d[!is.na(d$truth_cn) & !is.na(d$called_cn), ]
if (nrow(d) == 0) stop("no rows in table")
if (is.null(d$baseline_offset)) d$baseline_offset <- 0
d$ref <- ifelse(d$is_reference == 1, "reference", "haplotype")

# gene CN = absolute VCF CN minus the constant paralog baseline (0 LPA/C4, +2 GSTM1, +1 CYP2D6)
d$gene_cn <- d$called_cn - d$baseline_offset

# integer-only ticks: step 1 for narrow panels, coarser for wide ones (lpa spans ~24) to avoid crowding.
int_breaks <- function(lim) {
  hi <- floor(lim[2] + 1e-9)
  by <- if (hi <= 8) 1 else ceiling(hi / 8)
  seq(0, hi, by = by)
}

# Build a faceted concordance plot from a subset of genes (per-facet r, square axes, y=x line).
make_plot <- function(sub, ncol) {
  genes <- sort(unique(sub$gene))
  labs <- do.call(rbind, lapply(genes, function(g) {
    s <- sub[sub$gene == g, ]
    r <- suppressWarnings(tryCatch(cor(s$truth_cn, s$gene_cn), error = function(e) NA))
    # exact-match %: stays meaningful when truth is constant (CYP2D7=1), where Pearson r is undefined
    match <- mean(round(s$gene_cn) == s$truth_cn) * 100
    rlab <- if (is.na(r)) "r = n/a (truth constant)" else sprintf("r = %.3f", r)
    data.frame(gene = g,
               label = sprintf("%s\nmatch = %.0f%%\nn = %d", rlab, match, nrow(s)),
               stringsAsFactors = FALSE)
  }))
  # per-gene shared range so each panel is exactly square (invisible anchors define both free axes).
  lims <- do.call(rbind, lapply(genes, function(g) {
    s <- sub[sub$gene == g, ]
    hi <- max(c(s$truth_cn, s$gene_cn)) + 1   # one unit of headroom above the data
    data.frame(gene = g, v = c(-0.4, hi), stringsAsFactors = FALSE)
  }))
  ggplot(sub, aes(truth_cn, gene_cn)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
    geom_jitter(aes(colour = ref, size = ref), width = 0.12, height = 0.12, alpha = 0.75) +
    geom_blank(data = lims, aes(x = v, y = v)) +
    geom_text(data = labs, aes(x = -Inf, y = Inf, label = label),
              hjust = -0.1, vjust = 1.2, size = 3.2, colour = "grey20") +
    facet_wrap(~ gene, scales = "free", ncol = ncol) +
    scale_x_continuous(breaks = int_breaks, expand = expansion(mult = 0.02)) +
    scale_y_continuous(breaks = int_breaks, expand = expansion(mult = 0.02)) +
    scale_colour_manual(values = c(haplotype = "#2c7fb8", reference = "#d95f02"), name = NULL) +
    scale_size_manual(values = c(haplotype = 1.6, reference = 3.2), guide = "none") +
    labs(x = "ground-truth copy number",
         y = "panvar copy number (paralog baseline removed)") +
    theme_bw(base_size = 12) +
    theme(legend.position = "top", panel.grid.minor = element_blank(), aspect.ratio = 1)
}

# Split by label case: per-paralog split rows are the two upper-case genes; everything else is a locus total.
paralogs <- c("CYP2D6", "CYP2D7")
loci <- d[!(d$gene %in% paralogs), ]
genes <- d[d$gene %in% paralogs, ]

# Plot 1: the four loci as total counts (c4, cyp2d6, gstm1, lpa) in a 2x2 grid.
if (nrow(loci) > 0) {
  p_loci <- make_plot(loci, ncol = 2)
  ggsave(paste0(out_prefix, ".loci.png"), p_loci, width = wnum(w_arg, 8), height = wnum(h_arg, 7), dpi = dpi)
  cat("Wrote:", paste0(out_prefix, ".loci.png"), "\n")
} else {
  cat("(no locus-total rows; skipping .loci.png)\n")
}

# Plot 2: resolved CYP2D6 / CYP2D7 per-paralog split, one column of two stacked panels.
if (nrow(genes) > 0) {
  p_genes <- make_plot(genes, ncol = 1)
  ggsave(paste0(out_prefix, ".genes.png"), p_genes, width = wnum(w_arg, 4.6), height = wnum(h_arg, 8), dpi = dpi)
  cat("Wrote:", paste0(out_prefix, ".genes.png"), "\n")
} else {
  cat("(no CYP2D6/CYP2D7 per-gene rows; skipping .genes.png)\n")
}
