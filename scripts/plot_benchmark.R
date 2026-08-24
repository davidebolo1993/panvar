#!/usr/bin/env Rscript
# Unified round-trip benchmark plot: per-gene reconstruction anatomy as two independent plots arranged
# side by side (each with its own legend, no facet strips).
#   LEFT  ("Reconstruction"): stacked bar = Reconstructed (identity) + Residual, as % of the aligned
#         haplotype sequence. Use --left-ymin to zoom the y-axis (e.g. 90) since bars are near 100%.
#   RIGHT ("Variation found"): the TRUTH EVENT LEDGER, as % of aligned sequence -- bases in events
#         at or above the call threshold that no record covers (Missed) against bases in events below
#         it (Below threshold). This is sized from the walks, NOT from the alignment: the previous
#         version split the residual by contiguous edit-run length, which classified a clean 60 bp
#         deletion as a dozen sub-threshold runs and read Mis-called = 0 at every real locus.
# Needs a benchmark table with per-haplotype rows carrying locus, sum_aln_len, sum_delta,
# truth_missed_bp, truth_below_bp (emitted by `panvar benchmark`). Genes are worst-first, up to
# --per-row genes per row (each plot wraps the same way).
#   Rscript plot_benchmark.R --table combined.tsv --out benchmark [--per-row 30 --left-ymin 0 --top 0 --dpi 150]
suppressWarnings(suppressMessages({library(ggplot2); library(grid)}))

args <- commandArgs(trailingOnly = TRUE)
usage <- function(status = 0) {
  cat(paste(c(
    "plot_benchmark.R - round-trip reconstruction anatomy, one bar per locus.",
    "  LEFT  (Reconstruction): stacked bar = Reconstructed (identity) + Residual, as % of the",
    "        aligned haplotype sequence.",
    "  RIGHT (Variation found): the truth-event ledger -- bases in above-threshold events no record",
    "        covers (Missed) against bases in below-threshold events, y-axis auto-scaled.",
    "  Loci are ordered worst-first by missed bp.",
    "",
    "Usage:",
    "  Rscript plot_benchmark.R --table <combined.tsv> --out <prefix> [options]",
    "",
    "Required:",
    "  --loss <path>            optional per-locus loss partition; adds the 'Where the loss lives'",
    "                           panel. NOT written by `panvar benchmark` -- it is assembled from each",
    "                           locus's <prefix>.qv_summary.tsv (the loss_bp rows) by",
    "                           scripts/regen_results.sh, as results/benchmark_loss.tsv",
    "  --locus <name>           label a single run's <prefix>.qv_by_haplotype.tsv, which has no",
    "                           `locus` column of its own",
    "  --table <path>           per-haplotype benchmark rows; needs the columns locus, sum_aln_len,",
    "                           sum_delta, truth_missed_bp, truth_below_bp (from `panvar benchmark`)",
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
    "such as --min-maf leave their events attributed to no record, so they appear as Missed",
    "rather than as filtered."), collapse = "\n"), "\n")
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
lp <- get("--loss", "")
per_row <- as.integer(numarg("--per-row", "30")); if (per_row < 1) per_row <- 30L
top <- as.integer(numarg("--top", "0"))
left_ymin <- numarg("--left-ymin", "0"); dpi <- numarg("--dpi", "150")
if (is.null(tp)) usage(1)

d <- read.delim(tp, check.names = FALSE, stringsAsFactors = FALSE)
# A single `panvar benchmark` run writes <prefix>.qv_by_haplotype.tsv with NO locus column -- that is
# added by scripts/regen_results.sh when it concatenates loci into results/benchmark_qv.tsv. Accept a
# single run directly by naming it, rather than telling the user to re-run benchmark, which can never
# add the column.
locus_name <- get("--locus", "")
if (!("locus" %in% names(d)) && nzchar(locus_name)) d$locus <- locus_name
need <- c("locus", "sum_aln_len", "sum_delta", "truth_missed_bp", "truth_below_bp")
# The region-VCF level is optional so a table without it still plots, but the left panel then shows
# only the graph CEILING -- which reads ~100% everywhere and is not what a consumer of the VCF gets.
has_gt <- all(c("gt_sum_delta", "gt_sum_aln_len") %in% names(d)) &&
          any(!is.na(suppressWarnings(as.numeric(d$gt_sum_delta))))
miss <- setdiff(need, names(d))
if (length(miss)) {
  hint <- if (identical(miss, "locus"))
    paste0("this looks like one run's <prefix>.qv_by_haplotype.tsv, which carries no `locus` column.",
           " Pass --locus <name> to plot it on its own, or use results/benchmark_qv.tsv written by",
           " scripts/regen_results.sh, which prepends it.")
  else "expected the per-haplotype table written by `panvar benchmark`."
  stop("table missing columns: ", paste(miss, collapse = ", "), "\n  ", hint)
}
if (!has_gt)
  message("note: no region-VCF columns (gt_*) -- benchmark was run without --vcf, so only the graph",
          " ceiling is shown. Re-run `panvar benchmark --vcf <call>.region.vcf` for the level a",
          " consumer of the VCF actually gets.")

agg_cols <- c("sum_aln_len", "sum_delta", "truth_missed_bp", "truth_below_bp")
if (has_gt) agg_cols <- c(agg_cols, "gt_sum_delta", "gt_sum_aln_len")
a <- aggregate(d[agg_cols], by = list(locus = d$locus), FUN = sum)
a <- a[a$sum_aln_len > 0, ]
a$recon <- 100 * (a$sum_aln_len - a$sum_delta) / a$sum_aln_len   # GRAPH ceiling, % of aligned
if (has_gt) {
  a$gtrecon <- 100 * (a$gt_sum_aln_len - a$gt_sum_delta) / a$gt_sum_aln_len   # what the VCF gives
  a$gtresid <- 100 - a$gtrecon
}
a$resid <- 100 * a$sum_delta / a$sum_aln_len                     # residual,  % of aligned
a$ncall <- 100 * a$truth_below_bp / a$sum_aln_len               # ledger, absolute % of aligned
a$mis   <- 100 * a$truth_missed_bp / a$sum_aln_len

a <- a[order(-a$truth_missed_bp, -a$sum_delta, a$locus), ]       # worst-first
if (top > 0) a <- head(a, top)
if (nrow(a) == 0) stop("no loci to plot")
a$rank <- seq_len(nrow(a)); a$row <- ((a$rank - 1) %/% per_row) + 1
lv <- a$locus; nr <- max(a$row)

mkdf <- function(comps) do.call(rbind, lapply(names(comps), function(cn)
  data.frame(locus = a$locus, row = a$row, component = cn, value = a[[comps[[cn]]]], stringsAsFactors = FALSE)))
left  <- mkdf(c("Reconstructed" = "recon", "Residual" = "resid"))
right <- mkdf(c("Below threshold" = "ncall", "Missed" = "mis"))
left$locus  <- factor(left$locus,  levels = lv); left$component  <- factor(left$component,  levels = c("Reconstructed", "Residual"))
right$locus <- factor(right$locus, levels = lv); right$component <- factor(right$component, levels = c("Below threshold", "Missed"))

panel <- function(df, cols, ymin, ymax, title, legend_rows = 1) {
  ggplot(df, aes(locus, value, fill = component)) +
    geom_col(width = 0.85, position = position_stack(reverse = TRUE)) +
    facet_wrap(~row, ncol = 1, scales = "free_x") +
    # A five-term legend does not fit on one row beside three other panels, and the entry that gets
    # clipped is the last one -- which here is the SIGNED term, the one most in need of its label.
    scale_fill_manual(values = cols, name = NULL,
                      guide = guide_legend(nrow = legend_rows, byrow = TRUE)) +
    coord_cartesian(ylim = c(ymin, ymax)) +
    labs(x = "gene", y = "% of aligned sequence", title = title) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, size = 11),
          strip.background = element_blank(), strip.text = element_blank(),
          panel.grid.major.x = element_blank(), axis.text.x = element_text(angle = 60, hjust = 1, size = 7))
}
# The panels NAME their level. The left one is the graph ceiling: it implants the haplotype's own true
# block wherever a call shares a node with it, so it reads near 100% at every locus and is not what a
# consumer of the VCF gets. Publishing it under the bare word "Reconstruction" is how "every haplotype
# reconstructs above 99.9%" ended up in the walkthrough describing a number nobody could reproduce
# from the emitted records.
p_left  <- panel(left,  c("Reconstructed" = "#3a9679", "Residual" = "#dcdcdc"), left_ymin, 100,
                 "Graph ceiling (implants the true block)")
p_right <- panel(right, c("Below threshold" = "#74add1", "Missed" = "#d73027"), 0,
                 max(c(a$ncall, a$mis, a$resid)) * 1.05, "Variation found")
if (has_gt) {
  gt <- mkdf(c("Reconstructed" = "gtrecon", "Residual" = "gtresid"))
  gt$locus <- factor(gt$locus, levels = lv)
  gt$component <- factor(gt$component, levels = c("Reconstructed", "Residual"))
  p_gt <- panel(gt, c("Reconstructed" = "#2c7fb8", "Residual" = "#dcdcdc"),
                max(0, min(a$gtrecon) - 5), 100, "From the VCF alone (what a consumer gets)")
}

# The loss partition: five consecutive terms that sum EXACTLY to the genotype residual, so the panel
# answers WHY reconstruction falls short rather than by how much. Read left to right it is one
# question per step -- did we even try (out_of_scope, below the size threshold and never asked for),
# did we find it (discovery), did we put it on the right haplotype (carrier_missed), did we describe
# it correctly (representation), did we add something that is not there (false_positive_damage).
# The last term is SIGNED: an edit applied where there is no eligible truth event occasionally helps.
has_loss <- FALSE
if (nzchar(lp) && file.exists(lp)) {
  L <- read.delim(lp, check.names = FALSE, stringsAsFactors = FALSE)
  L <- L[L$locus %in% a$locus & L$baseline_bp > 0, , drop = FALSE]
  if (nrow(L) > 0) {
    has_loss <- TRUE
    terms <- c("Not attempted (sub-threshold)" = "out_of_scope",
               "Not found"                      = "discovery_or_attribution",
               "Wrong haplotype"                = "carrier_missed",
               "Wrongly represented"            = "representation",
               "Spurious edits"                 = "false_positive_damage")
    loss <- do.call(rbind, lapply(names(terms), function(cn)
      data.frame(locus = L$locus,
                 row = a$row[match(L$locus, a$locus)],
                 component = cn,
                 value = 100 * L[[terms[[cn]]]] / L$baseline_bp,
                 stringsAsFactors = FALSE)))
    loss$locus <- factor(loss$locus, levels = lv)
    loss$component <- factor(loss$component, levels = names(terms))
    p_loss <- panel(loss,
                    c("Not attempted (sub-threshold)" = "#bdbdbd",
                      "Not found"                     = "#d73027",
                      "Wrong haplotype"               = "#fc8d59",
                      "Wrongly represented"           = "#4575b4",
                      "Spurious edits"                = "#984ea3"),
                    min(0, min(loss$value)), max(loss$value) * 1.1,
                    "Where the loss lives (% of baseline)", legend_rows = 3)
  }
}

ncol <- if (has_gt) 3 else 2
if (has_loss) ncol <- ncol + 1
w <- max(9, min(30, ncol * (0.30 * per_row + 1.2)))
h <- 1.2 + 2.3 * nr
png(paste0(out, ".png"), width = w, height = h, units = "in", res = dpi)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, ncol)))
col <- 1
print(p_left, vp = viewport(layout.pos.row = 1, layout.pos.col = col)); col <- col + 1
if (has_gt)   { print(p_gt,   vp = viewport(layout.pos.row = 1, layout.pos.col = col)); col <- col + 1 }
if (has_loss) { print(p_loss, vp = viewport(layout.pos.row = 1, layout.pos.col = col)); col <- col + 1 }
print(p_right, vp = viewport(layout.pos.row = 1, layout.pos.col = col))
invisible(dev.off())
cat(sprintf("wrote %s.png  (%d genes, %d row(s), %d per row; left y %g-100)\n", out, nrow(a), nr, per_row, left_ymin))
