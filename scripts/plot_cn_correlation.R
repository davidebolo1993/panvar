#!/usr/bin/env Rscript
# Copy-number concordance: panvar calls (absolute CN from the final region VCF) vs the ground-truth
# BED, one facet per gene. Each point is a haplotype; the reference haplotype is highlighted. The
# Pearson correlation (and the y = x line) summarize agreement per gene.
#
#   Rscript plot_cn_correlation.R --table cn_table.tsv --out results/cn_correlation
#
# Input TSV (produced by scripts/compare_copy_number.py --dump-tsv):
#   gene  sample  truth_cn  called_cn  is_reference
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

d <- read.delim(table_path, stringsAsFactors = FALSE)
d <- d[!is.na(d$truth_cn) & !is.na(d$called_cn), ]
if (nrow(d) == 0) stop("no rows in table")
if (is.null(d$baseline_offset)) d$baseline_offset <- 0
d$ref <- ifelse(d$is_reference == 1, "reference", "haplotype")


# Recovered GENE copy number: absolute VCF CN minus the constant paralog baseline (0 for LPA/C4, where
# the call already counts the gene; +2 for GSTM1 and +1 for CYP2D6, the always-present folded paralogs).
d$gene_cn <- d$called_cn - d$baseline_offset

# Per-gene Pearson r + n + the subtracted baseline, placed in each facet.
genes <- sort(unique(d$gene))
labs <- do.call(rbind, lapply(genes, function(g) {
  s <- d[d$gene == g, ]
  r <- tryCatch(cor(s$truth_cn, s$gene_cn), error = function(e) NA)
  off <- s$baseline_offset[1]
  data.frame(gene = g,
             label = sprintf("r = %.3f\nn = %d\nbaseline %+d", r, nrow(s), off),
             stringsAsFactors = FALSE)
}))

# per-gene shared range so each panel is exactly square (x and y span the same integer range).
# invisible anchor points at -0.4 and (sharedmax + 1 headroom) fully define both free axes, so the
# range is deterministic and identical on x and y (not perturbed by the random jitter of real points).
lims <- do.call(rbind, lapply(genes, function(g) {
  s <- d[d$gene == g, ]
  hi <- max(c(s$truth_cn, s$gene_cn)) + 1   # one unit of headroom above the data
  data.frame(gene = g, v = c(-0.4, hi), stringsAsFactors = FALSE)
}))

# integer-only ticks: step 1 for narrow panels, coarser for wide ones (lpa spans ~24) to avoid crowding.
int_breaks <- function(lim) {
  hi <- floor(lim[2] + 1e-9)
  by <- if (hi <= 8) 1 else ceiling(hi / 8)
  seq(0, hi, by = by)
}

p <- ggplot(d, aes(truth_cn, gene_cn)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  geom_jitter(aes(colour = ref, size = ref), width = 0.12, height = 0.12, alpha = 0.75) +
  geom_blank(data = lims, aes(x = v, y = v)) +
  # pin the annotation to the upper-left corner of each panel (away from the points)
  geom_text(data = labs, aes(x = -Inf, y = Inf, label = label),
            hjust = -0.1, vjust = 1.2, size = 3.2, colour = "grey20") +
  facet_wrap(~ gene, scales = "free") +
  scale_x_continuous(breaks = int_breaks, expand = expansion(mult = 0.02)) +
  scale_y_continuous(breaks = int_breaks, expand = expansion(mult = 0.02)) +
  scale_colour_manual(values = c(haplotype = "#2c7fb8", reference = "#d95f02"), name = NULL) +
  scale_size_manual(values = c(haplotype = 1.6, reference = 3.2), guide = "none") +
  labs(x = "ground-truth copy number",
       y = "panvar copy number (paralog baseline removed)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "top", panel.grid.minor = element_blank(), aspect.ratio = 1)

png_path <- paste0(out_prefix, ".png")
ggsave(png_path, p, width = 9, height = 7, dpi = 150)
cat("Wrote:", png_path, "\n")
