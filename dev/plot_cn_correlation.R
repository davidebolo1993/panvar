#!/usr/bin/env Rscript
# Called against expected copy number, one panel per gene or locus, from compare_copy_number.py tables.
# See --help for usage.
suppressWarnings(suppressMessages(library(ggplot2)))
set.seed(1)   # the jitter: the same tables give the same figure

args <- commandArgs(trailingOnly = TRUE)
usage <- function(status = 0) {
  cat(paste(c(
    "plot_cn_correlation.R - called against expected copy number, one panel per gene or locus",
    "",
    "Usage:",
    "  Rscript plot_cn_correlation.R --table <cn.tsv> [--table <cn.tsv> ...] --out <prefix>",
    "",
    "  --table <path>   a compare_copy_number.py -o table; repeat for several comparisons",
    "  --out <prefix>   writes <prefix>.loci.png (module totals), <prefix>.genes.png (per-gene copy numbers)",
    "                   and <prefix>.summary.tsv (per panel: haplotypes, exact matches, r, hybrid copies)",
    "",
    "Haplotypes call found a hybrid copy in (the gene table's hybrid_copies) are drawn apart."), collapse = "\n"), "\n")
  quit(status = status)
}
if (length(args) == 0 || any(args %in% c("-h", "--help"))) usage(0)
tables <- args[which(args == "--table") + 1]
i <- match("--out", args)
out_prefix <- if (is.na(i) || i == length(args)) NULL else args[i + 1]
if (length(tables) == 0 || is.null(out_prefix)) usage(1)

d <- do.call(rbind, lapply(tables, function(p) {
  if (!file.exists(p)) stop("not found: ", p)
  t <- read.delim(p, stringsAsFactors = FALSE)
  if (is.null(t$hybrid)) t$hybrid <- 0
  t[, c("label", "sample", "truth_cn", "called_cn", "is_reference", "level", "hybrid")]
}))
if (nrow(d) == 0) stop("no rows in the tables")
d$kind <- factor(ifelse(d$is_reference == 1, "reference", ifelse(d$hybrid == 1, "hybrid copy", "haplotype")),
                 c("haplotype", "hybrid copy", "reference"))

summary <- do.call(rbind, lapply(split(d, list(d$level, d$label), drop = TRUE), function(s) {
  # no correlation from one haplotype (sd is NA) or a constant column
  r <- if (nrow(s) > 1 && sd(s$truth_cn) > 0 && sd(s$called_cn) > 0) round(cor(s$truth_cn, s$called_cn), 3) else NA
  data.frame(label = s$label[1], level = s$level[1], n = nrow(s), exact = sum(s$truth_cn == s$called_cn),
             match_pct = round(100 * mean(s$truth_cn == s$called_cn), 1), r = r,
             hybrids = sum(s$hybrid == 1), hybrids_exact = sum(s$hybrid == 1 & s$truth_cn == s$called_cn))
}))
summary <- summary[order(summary$level, summary$label), ]
write.table(summary, paste0(out_prefix, ".summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("wrote", paste0(out_prefix, ".summary.tsv"), "\n")

int_breaks <- function(lim) { hi <- floor(lim[2] + 1e-9); seq(0, hi, by = if (hi <= 8) 1 else ceiling(hi / 8)) }

draw <- function(level, suffix) {
  sub <- d[d$level == level, ]
  if (nrow(sub) == 0) return(invisible())
  sub <- sub[order(sub$kind), ]   # hybrids and the reference drawn over the other haplotypes
  sm <- summary[summary$level == level, ]
  labs <- data.frame(label = sm$label, text = paste0(sprintf("%s\nexact = %d/%d",
                     ifelse(is.na(sm$r), "r = n/a (constant)", sprintf("r = %.3f", sm$r)), sm$exact, sm$n),
                     ifelse(sm$hybrids > 0, sprintf("\nhybrid copy: %d (%d exact)", sm$hybrids, sm$hybrids_exact), "")))
  # square panels: both axes span the same range
  lims <- do.call(rbind, lapply(unique(sub$label), function(g) {
    s <- sub[sub$label == g, ]
    data.frame(label = g, v = c(-0.4, max(c(s$truth_cn, s$called_cn)) + 1))
  }))
  n <- length(unique(sub$label)); ncol <- min(4, n)
  p <- ggplot(sub, aes(truth_cn, called_cn)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
    geom_jitter(aes(colour = kind, size = kind, shape = kind), width = 0.12, height = 0.12, alpha = 0.75) +
    geom_blank(data = lims, aes(x = v, y = v)) +
    geom_text(data = labs, aes(x = -Inf, y = Inf, label = gsub("(^|\n)", "\\1  ", text)), hjust = 0, vjust = 1.2, size = 3.2,
              colour = "grey20") +
    facet_wrap(~ label, scales = "free", ncol = ncol) +
    scale_x_continuous(breaks = int_breaks, expand = expansion(mult = 0.02)) +
    scale_y_continuous(breaks = int_breaks, expand = expansion(mult = 0.02)) +
    scale_colour_manual(values = c(haplotype = "#2c7fb8", "hybrid copy" = "#c51b7d", reference = "#d95f02"),
                        name = NULL, drop = TRUE) +
    scale_shape_manual(values = c(haplotype = 16, "hybrid copy" = 17, reference = 16), name = NULL, drop = TRUE) +
    scale_size_manual(values = c(haplotype = 1.6, "hybrid copy" = 2.4, reference = 3.2), guide = "none") +
    labs(x = "expected copy number (truth set)", y = "panvar copy number") +
    theme_bw(base_size = 12) +
    theme(legend.position = "top", panel.grid.minor = element_blank(), aspect.ratio = 1)
  out <- paste0(out_prefix, suffix)
  ggsave(out, p, width = 3.5 * ncol + 1, height = 3.5 * ceiling(n / ncol) + 1, dpi = 150)
  cat("wrote", out, "\n")
}
draw("locus", ".loci.png")
draw("gene", ".genes.png")
