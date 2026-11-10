#!/usr/bin/env Rscript
# Precision, recall and F1 of panvar call and svim-asm (compare_sv_callers.py's tables) as one figure.
# See --help for usage.
suppressWarnings(suppressMessages(library(ggplot2)))

args <- commandArgs(trailingOnly = TRUE)
usage <- function(status = 0) {
  cat(paste(c(
    "plot_sv_callers.R - compare_sv_callers.py's tables as one figure",
    "",
    "Usage:",
    "  Rscript plot_sv_callers.R [options] <table.tsv> [<table.tsv> ...]",
    "",
    "One table per locus (compare_sv_callers.py -o). Rows: precision, recall, F1; columns: the callers' own",
    "per-haplotype calls and their merged multi-sample output; the loci along the x axis.",
    "",
    "Optional:",
    "  --by-size          columns by event size instead (50-999 bp, 1-10 kb, >10 kb), for --level",
    "  --level <own|merged>  the level drawn with --by-size (default merged)",
    "  --out <path>       output image (default callers.png beside the first table)",
    "  --dpi <N>          resolution (default 150)",
    "  --width <in>, --height <in>",
    "                     size in inches (default 1.2 per locus + 4 per column, 7)",
    "  -h, --help         show this help"), collapse = "\n"), "\n")
  quit(status = status)
}
if (length(args) == 0 || any(args %in% c("-h", "--help"))) usage(0)
valued <- c("--out", "--dpi", "--width", "--height", "--level")
get <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i == length(args)) stop("missing value after ", flag)
  args[i + 1]
}
by_size <- "--by-size" %in% args
taken <- unlist(lapply(valued, function(f) { i <- match(f, args); if (is.na(i)) integer(0) else c(i, i + 1) }))
taken <- c(taken, which(args == "--by-size"))
tables <- if (length(taken)) args[-taken] else args
if (length(tables) == 0) usage(1)
level <- get("--level", "merged")
out <- get("--out", file.path(dirname(tables[1]), if (by_size) "callers.by_size.png" else "callers.png"))
dpi <- as.numeric(get("--dpi", "150"))

d <- do.call(rbind, lapply(tables, function(p) read.delim(p, colClasses = "character", na.strings = "NA", check.names = FALSE)))
loci <- unique(d$locus)
ncol <- if (by_size) 3 else 2
width <- as.numeric(get("--width", as.character(1.2 * length(loci) * ncol / 2 + 4 * ncol / 2 + 2)))
height <- as.numeric(get("--height", "7"))

sizes <- c("50-999 bp", "1-10 kb", ">10 kb")
if (by_size) {
  d <- d[d$level == level & d$stratum %in% sizes, ]
  d$column <- factor(d$stratum, sizes)
} else {
  d <- d[d$stratum == "all", ]
  d$column <- factor(ifelse(d$level == "own", "each haplotype's own calls", "merged output (VCF)"),
                     c("each haplotype's own calls", "merged output (VCF)"))
}
long <- do.call(rbind, lapply(c("precision", "recall", "f1"), function(m)
  data.frame(locus = d$locus, caller = d$caller, column = d$column,
             metric = c(precision = "precision", recall = "recall", f1 = "F1")[[m]],
             value = suppressWarnings(as.numeric(d[[m]])),
             n = if (m == "precision") d$calls else d$truth_events)))
long$locus <- factor(long$locus, loci)
long$metric <- factor(long$metric, c("precision", "recall", "F1"))
long$caller <- factor(long$caller, c("panvar call", "svim-asm"))
long$label <- ifelse(is.na(long$value), "", sprintf("%.2f", long$value))

dodge <- position_dodge(width = 0.8)
fig <- ggplot(long, aes(locus, value, fill = caller)) +
  geom_col(width = 0.75, position = dodge) +
  geom_text(aes(label = label), position = dodge, vjust = -0.35, size = 2.6) +
  facet_grid(metric ~ column, switch = "y") +
  scale_y_continuous(limits = c(0, 1.12), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  scale_fill_manual(values = c("panvar call" = "#2e75b6", "svim-asm" = "#e69138")) +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "top", panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        strip.placement = "outside", strip.background = element_rect(fill = "grey95"),
        axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(out, fig, width = width, height = height, dpi = dpi)
cat("wrote", out, "\n")
