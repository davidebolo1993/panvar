#!/usr/bin/env Rscript
# Unified round-trip benchmark plot: per-gene reconstruction anatomy as two independent plots arranged
# side by side (each with its own legend, no facet strips).
#   LEFT  ("Reconstruction"): stacked bar = Reconstructed (identity) + Residual, as % of the aligned
#         haplotype sequence. Use --left-ymin to zoom the y-axis (e.g. 90) since bars are near 100%.
#   RIGHT ("Residual composition"): the residual (= 100% - Reconstructed) split into Not-callable
#         (sub-threshold) vs Mis-called (>= threshold), also as % of aligned sequence, with the y-axis
#         auto-scaled to the largest residual so the (usually tiny) shortfall is legible.
# Needs a benchmark table with per-haplotype rows carrying locus, sum_aln_len, sum_delta,
# sub_threshold_bp, over_threshold_bp (emitted by `panvar benchmark`). Genes are worst-first, up to
# --per-row genes per row (each plot wraps the same way).
#   Rscript plot_benchmark.R --table combined.tsv --out benchmark [--per-row 30 --left-ymin 0 --top 0 --dpi 150]
suppressWarnings(suppressMessages({library(ggplot2); library(grid)}))

args <- commandArgs(trailingOnly = TRUE)
usage <- function(status = 0) {
  cat(paste(c(
    "plot_benchmark.R - round-trip reconstruction anatomy, one bar per locus.",
    "  LEFT  (Reconstruction): stacked bar = Reconstructed (identity) + Residual, as % of the",
    "        aligned haplotype sequence.",
    "  RIGHT (Residual composition): that residual split into Not-callable (< threshold) and",
    "        Mis-called (>= threshold), y-axis auto-scaled so the usually tiny shortfall is legible.",
    "  Loci are ordered worst-first by mis-called bp.",
    "",
    "Usage:",
    "  Rscript plot_benchmark.R --table <combined.tsv> --out <prefix> [options]",
    "",
    "Required:",
    "  --table <path>           per-haplotype benchmark rows; needs the columns locus, sum_aln_len,",
    "                           sum_delta, sub_threshold_bp, over_threshold_bp (from `panvar benchmark`)",
    "",
    "Optional:",
    "  --out <prefix>           output prefix; writes <prefix>.png (default benchmark)",
    "  --per-row <N>            loci per row; both plots wrap the same way (default 30)",
    "  --top <N>                keep only the N worst loci (default 0 = all)",
    "  --left-ymin <y>          zoom the left y-axis, e.g. 90, since bars sit near 100% (default 0)",
    "  --dpi <n>                PNG resolution (default 150)",
    "  -h, --help               show this help",
    "",
    "Note: `benchmark` reads the post-filter variant_nodes.tsv, so calls dropped by a `call` filter",
    "such as --min-maf are counted as Mis-called rather than as filtered."), collapse = "\n"), "\n")
  quit(status = status)
}
if (length(args) == 0 || any(args %in% c("-h", "--help"))) usage(0)
get <- function(f, d = NULL) {
  i <- match(f, args)
  if (is.na(i) || i == length(args)) return(d)
  v <- args[i + 1]
  if (grepl("^--", v)) return(d)            # next token is another flag -> this flag got no value
  v
}
numarg <- function(f, d) { v <- suppressWarnings(as.numeric(get(f, d))); if (is.na(v)) as.numeric(d) else v }
tp <- get("--table"); out <- get("--out", "benchmark")
per_row <- as.integer(numarg("--per-row", "30")); if (per_row < 1) per_row <- 30L
top <- as.integer(numarg("--top", "0"))
left_ymin <- numarg("--left-ymin", "0"); dpi <- numarg("--dpi", "150")
if (is.null(tp)) usage(1)

d <- read.delim(tp, check.names = FALSE, stringsAsFactors = FALSE)
need <- c("locus", "sum_aln_len", "sum_delta", "sub_threshold_bp", "over_threshold_bp")
miss <- setdiff(need, names(d))
if (length(miss)) stop("table missing columns: ", paste(miss, collapse = ", "),
                       "  (re-run `panvar benchmark` to add the residual-split columns)")

a <- aggregate(cbind(sum_aln_len, sum_delta, sub_threshold_bp, over_threshold_bp) ~ locus, d, sum)
a <- a[a$sum_aln_len > 0, ]
a$recon <- 100 * (a$sum_aln_len - a$sum_delta) / a$sum_aln_len   # identity, % of aligned
a$resid <- 100 * a$sum_delta / a$sum_aln_len                     # residual,  % of aligned
a$ncall <- 100 * a$sub_threshold_bp / a$sum_aln_len             # residual split, absolute % of aligned
a$mis   <- 100 * a$over_threshold_bp / a$sum_aln_len

a <- a[order(-a$over_threshold_bp, -a$sum_delta, a$locus), ]     # worst-first
if (top > 0) a <- head(a, top)
if (nrow(a) == 0) stop("no loci to plot")
a$rank <- seq_len(nrow(a)); a$row <- ((a$rank - 1) %/% per_row) + 1
lv <- a$locus; nr <- max(a$row)

mkdf <- function(comps) do.call(rbind, lapply(names(comps), function(cn)
  data.frame(locus = a$locus, row = a$row, component = cn, value = a[[comps[[cn]]]], stringsAsFactors = FALSE)))
left  <- mkdf(c("Reconstructed" = "recon", "Residual" = "resid"))
right <- mkdf(c("Not-callable" = "ncall", "Mis-called" = "mis"))
left$locus  <- factor(left$locus,  levels = lv); left$component  <- factor(left$component,  levels = c("Reconstructed", "Residual"))
right$locus <- factor(right$locus, levels = lv); right$component <- factor(right$component, levels = c("Not-callable", "Mis-called"))

panel <- function(df, cols, ymin, ymax, title) {
  ggplot(df, aes(locus, value, fill = component)) +
    geom_col(width = 0.85, position = position_stack(reverse = TRUE)) +
    facet_wrap(~row, ncol = 1, scales = "free_x") +
    scale_fill_manual(values = cols, name = NULL) +
    coord_cartesian(ylim = c(ymin, ymax)) +
    labs(x = "gene", y = "% of aligned sequence", title = title) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, size = 11),
          strip.background = element_blank(), strip.text = element_blank(),
          panel.grid.major.x = element_blank(), axis.text.x = element_text(angle = 60, hjust = 1, size = 7))
}
p_left  <- panel(left,  c("Reconstructed" = "#3a9679", "Residual" = "#dcdcdc"), left_ymin, 100, "Reconstruction")
p_right <- panel(right, c("Not-callable" = "#74add1", "Mis-called" = "#d73027"), 0, max(a$resid) * 1.05, "Residual composition")

w <- max(9, min(22, 2 * (0.30 * per_row + 1.2)))
h <- 1.2 + 2.3 * nr
png(paste0(out, ".png"), width = w, height = h, units = "in", res = dpi)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
print(p_left,  vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p_right, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
invisible(dev.off())
cat(sprintf("wrote %s.png  (%d genes, %d row(s), %d per row; left y %g-100)\n", out, nrow(a), nr, per_row, left_ymin))
