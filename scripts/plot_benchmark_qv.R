#!/usr/bin/env Rscript
# Round-trip QV summary (cosigt-style): one stacked bar per locus, the fraction of haplotypes whose
# reconstruction falls in each QV category. Reads the combined benchmark table
# (locus sample sum_delta sum_aln_len qv band qv_ratio quintile identity), one row per haplotype.
#   Rscript plot_benchmark_qv.R --table results/benchmark_qv.tsv --out results/benchmark_qv [--category quintile|band|identity]
# Default category is `quintile` (qv/qv_max, length-fair): a perfect reconstruction of any locus size
# lands in the top bin. `band` is the raw cosigt QV band; `identity` is 1 - Σδ/ΣS, binned linearly so a
# few uncalled SNPs over a long region do not drop the bar.
suppressWarnings(suppressMessages(library(ggplot2)))

args <- commandArgs(trailingOnly = TRUE)
get <- function(flag, default = NULL) { i <- match(flag, args); if (is.na(i) || i == length(args)) return(default); args[i + 1] }
table_path <- get("--table"); out_prefix <- get("--out", "benchmark_qv")
category <- get("--category", "quintile")
dpi <- as.numeric(get("--dpi", "150"))
if (is.null(table_path)) stop("need --table <tsv>")

d <- read.delim(table_path, stringsAsFactors = FALSE)
if (nrow(d) == 0) stop("no rows in table")

# Continuous "how close to 100%" view: the per-haplotype reconstruction identity (1 - Σδ/ΣS) as a
# distribution per locus, on a "nines" axis (-log10(1 - identity)) so the near-1.0 spread is visible.
# Bar/quintile bins collapse when everything is near-perfect; this does not.
if (category == "dist") {
  d$nines <- pmin(6, -log10(pmax(1e-6, 1 - d$identity)))  # 3 = 99.9%, 4 = 99.99%, ...
  med <- sort(tapply(d$identity, d$locus, median), decreasing = TRUE)
  d$locus <- factor(d$locus, levels = names(med))
  p <- ggplot(d, aes(locus, nines)) +
    geom_violin(fill = "#a6cee3", colour = NA, alpha = 0.5, scale = "width") +
    geom_boxplot(width = 0.18, outlier.size = 0.5, fill = "white") +
    scale_y_continuous(breaks = 1:6, labels = c("90%", "99%", "99.9%", "99.99%", "99.999%", "≥99.9999%")) +
    labs(x = NULL, y = "reconstruction identity  (1 - Σδ/ΣS)",
         title = "Round-trip reconstruction identity per locus") +
    theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank())
  out <- paste0(out_prefix, ".dist.png")
  ggsave(out, p, width = as.numeric(get("--width", "7")), height = as.numeric(get("--height", "5")), dpi = dpi)
  cat("Wrote:", out, "\n"); quit(status = 0)
}

if (category == "band") {
  d$cat <- factor(d$band, levels = c("<17", "17-23", "23-33", ">33"))
  pal <- c("<17" = "#d7191c", "17-23" = "#fdae61", "23-33" = "#a6d96a", ">33" = "#1a9641")
  ylab <- "haplotypes (%)"; title <- "Round-trip QV — cosigt bands"
} else if (category == "identity") {
  b <- cut(d$identity, breaks = c(-Inf, 0.9, 0.99, 0.999, 0.9999, Inf),
           labels = c("<0.9", "0.9-0.99", "0.99-0.999", "0.999-0.9999", "≥0.9999"))
  d$cat <- factor(b, levels = c("<0.9", "0.9-0.99", "0.99-0.999", "0.999-0.9999", "≥0.9999"))
  pal <- c("<0.9" = "#d7191c", "0.9-0.99" = "#fdae61", "0.99-0.999" = "#ffffbf",
           "0.999-0.9999" = "#a6d96a", "≥0.9999" = "#1a9641")
  ylab <- "haplotypes (%)"; title <- "Round-trip reconstruction identity (1 - Σδ/ΣS)"
} else {
  lv <- c("0.0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0")
  d$cat <- factor(d$quintile, levels = lv)
  pal <- c("0.0-0.2" = "#d7191c", "0.2-0.4" = "#fdae61", "0.4-0.6" = "#ffffbf",
           "0.6-0.8" = "#a6d96a", "0.8-1.0" = "#1a9641")
  ylab <- "haplotypes (%)"; title <- "Round-trip QV / QV_max (length-fair, top = best)"
}

# Per-locus percentages, ordered by best category share so the plot reads left-to-right.
tab <- as.data.frame(table(locus = d$locus, cat = d$cat))
tot <- tapply(tab$Freq, tab$locus, sum)
tab$pct <- 100 * tab$Freq / tot[as.character(tab$locus)]
best <- levels(d$cat)[nlevels(d$cat)]
ord <- names(sort(tapply(tab$pct[tab$cat == best], tab$locus[tab$cat == best], sum), decreasing = TRUE))
tab$locus <- factor(tab$locus, levels = ord)

p <- ggplot(tab, aes(locus, pct, fill = cat)) +
  geom_col(width = 0.7, colour = "grey30", linewidth = 0.2) +
  scale_fill_manual(values = pal, name = NULL, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = ylab, title = title) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right", panel.grid.major.x = element_blank())

out <- paste0(out_prefix, ".", category, ".png")
ggsave(out, p, width = as.numeric(get("--width", "7")), height = as.numeric(get("--height", "5")), dpi = dpi)
cat("Wrote:", out, "\n")
