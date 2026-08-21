#!/usr/bin/env Rscript
# plot_genotype_calibration.R - is a genotype quality score believable?
#
# A GQ of 30 claims a 1-in-1000 error rate. This plots what actually happened against what GQ
# promised, so overconfidence is visible rather than assumed. Points below the diagonal are
# overconfident, which is the dangerous direction: the caller is surer than it should be.
#
# Input is the per-call table written by tests/genotype_sim.sh:
#   locus depth error loo pair block_kind n_alleles n_markers gq explained correct filter
#
# Usage:
#   plot_genotype_calibration.R --calls <calls.tsv> --out <prefix> [--title T] [--dpi N]

suppressWarnings(suppressMessages({
  ok <- require(ggplot2, quietly = TRUE)
}))

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, default = NA) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1]
}
if ("--help" %in% args || "-h" %in% args || is.na(getarg("--calls"))) {
  cat(paste(readLines(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))[2:12],
            collapse = "\n"), "\n")
  quit(status = 0)
}

calls_path <- getarg("--calls")
out_prefix <- getarg("--out", "genotype_calibration")
title <- getarg("--title", "Genotype quality calibration")
dpi <- as.numeric(getarg("--dpi", "150"))

d <- read.delim(calls_path, stringsAsFactors = FALSE)
if (nrow(d) == 0) stop("no calls in ", calls_path)

# Bin on the GQ scale rather than uniformly: GQ is logarithmic, so equal-width bins would put almost
# everything in the top bucket and hide exactly the range where miscalibration matters.
breaks <- c(0, 5, 15, 25, 40, 60, 90, 100)
d$bin <- cut(d$gq, breaks = breaks, include.lowest = TRUE, right = FALSE)
agg <- aggregate(cbind(correct, n = rep(1, nrow(d))) ~ bin, data = d, FUN = sum)
agg$observed <- agg$correct / agg$n
mids <- (head(breaks, -1) + tail(breaks, -1)) / 2
agg$implied <- 1 - 10^(-mids[as.integer(agg$bin)] / 10)
# Wilson interval: n per bin is small, and a bare point estimate would overstate what we know.
z <- 1.96
p <- agg$observed; n <- agg$n
den <- 1 + z^2 / n
agg$lo <- pmax(0, (p + z^2 / (2 * n) - z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / den)
agg$hi <- pmin(1, (p + z^2 / (2 * n) + z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / den)

write.table(agg[, c("bin", "n", "observed", "implied", "lo", "hi")],
            paste0(out_prefix, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

if (!ok) {
  cat("ggplot2 not available; wrote ", out_prefix, ".tsv only\n", sep = "")
  quit(status = 0)
}

pl <- ggplot(agg, aes(x = implied, y = observed)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.01, colour = "grey40") +
  geom_point(aes(size = n), colour = "#377EB8") +
  scale_x_continuous("accuracy implied by GQ", limits = c(0, 1), labels = scales::percent) +
  scale_y_continuous("observed accuracy", limits = c(0, 1), labels = scales::percent) +
  scale_size_continuous("calls", range = c(2, 7)) +
  labs(title = title,
       subtitle = "dashed = perfect calibration; below it = overconfident") +
  theme_bw(base_size = 11)

ggsave(paste0(out_prefix, ".png"), pl, width = 6, height = 5, dpi = dpi)
cat("wrote ", out_prefix, ".png and ", out_prefix, ".tsv\n", sep = "")
