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
             x = min(s$truth_cn), y = max(c(s$gene_cn, s$truth_cn)),
             label = sprintf("r = %.3f\nn = %d\nbaseline %+d", r, nrow(s), off),
             stringsAsFactors = FALSE)
}))

p <- ggplot(d, aes(truth_cn, gene_cn)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  geom_jitter(aes(colour = ref, size = ref), width = 0.12, height = 0.12, alpha = 0.75) +
  geom_text(data = labs, aes(x = x, y = y, label = label),
            hjust = 0, vjust = 1, size = 3.2, colour = "grey20") +
  facet_wrap(~ gene, scales = "free") +
  scale_colour_manual(values = c(haplotype = "#2c7fb8", reference = "#d95f02"), name = NULL) +
  scale_size_manual(values = c(haplotype = 1.6, reference = 3.2), guide = "none") +
  labs(x = "ground-truth copy number",
       y = "panvar copy number (paralog baseline removed)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "top", panel.grid.minor = element_blank())

png_path <- paste0(out_prefix, ".png")
ggsave(png_path, p, width = 9, height = 7, dpi = 150)
cat("Wrote:", png_path, "\n")
