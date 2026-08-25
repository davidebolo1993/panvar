#!/usr/bin/env Rscript
# Per-block genotyping error, one plot per locus: two rows (leave-zero-out above, leave-one-out
# below), one violin pair per block (real reads against simulated), blocks in chain order.
#
# The metric is NORMALIZED EDIT DISTANCE to truth, 1 - identity, so blocks of wildly different length
# sit on one axis. The same metric is used in both rows, which is the point: LZO is the control where
# the answer is in the panel and the error should be ~0, and any distance from 0 there is read
# acquisition or model error rather than an unrepresentable truth.
#
# LOO carries an irreducible part -- the individual's own haplotypes are removed, so even a perfect
# model can only reach the nearest remaining panel pair. That floor (1 - best_identity) is drawn as a
# black point per block, so the violin ABOVE the point is the model's own error and the point itself
# is what the panel costs. Reporting the violin alone would blame the model for a missing haplotype.
#
#   Rscript plot_genotype_blocks.R --table all_blocks.tsv --out DIR [--loci c4,lpa] [--dpi 150]
suppressWarnings(suppressMessages({library(ggplot2)}))

args <- commandArgs(trailingOnly = TRUE)
usage <- function(status = 0) {
  cat(paste(c(
    "plot_genotype_blocks.R - per-block genotyping error, one plot per locus.",
    "",
    "  Two rows per locus: leave-zero-out (truth in the panel) above, leave-one-out below.",
    "  One violin pair per block: real reads against simulated, on normalized edit distance",
    "  to truth (1 - identity). The black point marks the panel floor (1 - best_identity),",
    "  which is what remains unreachable once the individual's haplotypes are removed.",
    "",
    "Usage:",
    "  Rscript plot_genotype_blocks.R --table <all_blocks.tsv> --out <dir> [options]",
    "",
    "Required:",
    "  --table <f>     per-block cohort table with locus, read_source, regime, block_index,",
    "                  block_kind, sample, identity, best_identity, exact",
    "  --out <dir>     output directory; writes <dir>/<locus>.png and <dir>/summary.tsv",
    "",
    "Options:",
    "  --loci <a,b>    only these loci (default: all present)",
    "  --dpi <n>       raster resolution (default 150)",
    "  --width <in>    plot width; default scales with the block count",
    "  -h, --help      this message",
    ""), collapse = "\n"))
  quit(status = status)
}
if (length(args) == 0 || "-h" %in% args || "--help" %in% args) usage()

getarg <- function(flag, default = NA) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i == length(args)) stop(paste(flag, "needs a value"))
  args[i + 1]
}
table_path <- getarg("--table")
out_dir    <- getarg("--out")
if (is.na(table_path) || is.na(out_dir)) usage(1)
if (!file.exists(table_path)) stop(paste("no such table:", table_path))
dpi   <- as.numeric(getarg("--dpi", 150))
width_override <- suppressWarnings(as.numeric(getarg("--width", NA)))
only_loci <- getarg("--loci", NA)

d <- read.delim(table_path, stringsAsFactors = FALSE, check.names = FALSE)
need <- c("locus", "read_source", "regime", "block_index", "block_kind",
          "sample", "identity", "best_identity", "filter")
missing <- setdiff(need, names(d))
if (length(missing)) stop(paste("table lacks columns:", paste(missing, collapse = ", ")))

d$identity      <- suppressWarnings(as.numeric(d$identity))
d$best_identity <- suppressWarnings(as.numeric(d$best_identity))
d <- d[is.finite(d$identity), ]
if (!nrow(d)) stop("no rows with a numeric identity")

# Distance to truth, and the part of it the panel cannot avoid once the donor is held out.
d$err   <- 1 - d$identity
d$floor <- ifelse(is.finite(d$best_identity), 1 - d$best_identity, NA_real_)

if (!is.na(only_loci)) {
  keep <- strsplit(only_loci, ",")[[1]]
  d <- d[d$locus %in% keep, ]
  if (!nrow(d)) stop("no rows for the requested loci")
}

d$regime <- factor(ifelse(d$regime %in% c("LZO", "0"), "leave-zero-out (truth in panel)",
                          "leave-one-out (donor removed)"),
                   levels = c("leave-zero-out (truth in panel)", "leave-one-out (donor removed)"))
d$reads <- factor(d$read_source, levels = c("real", "simulated"))
d <- d[!is.na(d$reads), ]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
pal <- c(real = "#B2182B", simulated = "#2166AC")

# Per (locus, regime, block, reads): the exact rate is what "genotyped properly" means -- the called
# allele pair IS the truth pair. Kept next to the error so a block that is usually exact but
# occasionally catastrophic is not confused with one that is uniformly mediocre.
if ("exact" %in% names(d)) {
  d$exact <- suppressWarnings(as.numeric(d$exact))
} else {
  d$exact <- NA_real_
}
agg <- aggregate(cbind(err, exact, floor) ~ locus + regime + block_index + block_kind + reads,
                 data = d, FUN = function(x) mean(x, na.rm = TRUE), na.action = na.pass)
names(agg)[names(agg) == "err"] <- "mean_err"
n_obs <- aggregate(err ~ locus + regime + block_index + reads, data = d, FUN = length)
names(n_obs)[names(n_obs) == "err"] <- "n"
agg <- merge(agg, n_obs, by = c("locus", "regime", "block_index", "reads"))
agg <- agg[order(agg$locus, agg$regime, agg$block_index, agg$reads), ]
write.table(agg, file.path(out_dir, "summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

for (loc in sort(unique(d$locus))) {
  sub <- d[d$locus == loc, ]
  blocks <- sort(unique(sub$block_index))
  # Label each block with its kind, since "which block is the bubble" is the first thing asked.
  kind_of <- sapply(blocks, function(b) {
    k <- unique(sub$block_kind[sub$block_index == b]); k[1]
  })
  lab <- paste0(blocks, "\n", substr(kind_of, 1, 4))
  sub$blk <- factor(sub$block_index, levels = blocks, labels = lab)

  fl <- unique(sub[, c("regime", "blk", "floor", "block_index")])
  fl <- fl[is.finite(fl$floor), ]
  fl <- aggregate(floor ~ regime + blk, data = fl, FUN = mean)

  # Which of this error the module already declines. A reader has to be able to tell "we call this
  # wrong" from "we decline to call this", because only the first is a defect -- the second is the
  # module correctly saying the block cannot be genotyped on the evidence it has.
  sub$kept <- ifelse(sub$filter == "PASS", "reported (PASS)", "declined")
  sub$kept <- factor(sub$kept, levels = c("reported (PASS)", "declined"))

  p <- ggplot(sub, aes(x = blk, y = err, fill = reads)) +
    # Violins collapse to a line when every sample is exact, which is the correct picture for a
    # solved block; the jittered points keep the sample count visible when that happens.
    geom_violin(position = position_dodge(width = 0.8), width = 0.75,
                scale = "width", alpha = 0.55, colour = NA, na.rm = TRUE) +
    geom_point(position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.8,
                                               seed = 42),
               aes(colour = reads, shape = kept, size = kept, alpha = kept), na.rm = TRUE) +
    scale_shape_manual(values = c("reported (PASS)" = 16, "declined" = 1), name = NULL) +
    scale_size_manual(values = c("reported (PASS)" = 0.75, "declined" = 0.9), guide = "none") +
    scale_alpha_manual(values = c("reported (PASS)" = 0.75, "declined" = 0.45), guide = "none") +
    geom_point(data = fl, aes(x = blk, y = floor), inherit.aes = FALSE,
               shape = 18, size = 1.9, colour = "black", na.rm = TRUE) +
    facet_wrap(~ regime, ncol = 1, scales = "free_y") +
    scale_fill_manual(values = pal, name = "reads") +
    scale_colour_manual(values = pal, guide = "none") +
    labs(title = paste0(loc, " - per-block genotyping error"),
         subtitle = paste0("normalized edit distance to truth (1 - identity). Filled = reported, ",
                           "hollow = declined by the module. Black diamond = panel floor ",
                           "(1 - best_identity), unreachable once the donor is held out."),
         x = "block (chain order, kind abbreviated)", y = "1 - identity") +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey92", colour = NA),
          strip.text = element_text(face = "bold"),
          legend.position = "top",
          plot.subtitle = element_text(size = 7, colour = "grey30"))

  w <- if (is.finite(width_override)) width_override else max(6, 0.42 * length(blocks) + 2.2)
  ggsave(file.path(out_dir, paste0(loc, ".png")), p,
         width = w, height = 6.0, dpi = dpi, limitsize = FALSE)
  cat(sprintf("wrote %s (%d blocks, %d observations)\n",
              file.path(out_dir, paste0(loc, ".png")), length(blocks), nrow(sub)))
}
