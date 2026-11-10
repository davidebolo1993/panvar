#!/usr/bin/env Rscript
# Graph structure side by side, one panel per measure. See --help for usage.
suppressWarnings(suppressMessages(library(ggplot2)))

args <- commandArgs(trailingOnly = TRUE)
usage <- function(status = 0) {
  cat(paste(c(
    "plot_compare_graphs.R - compare_graphs.py's tables as one figure",
    "",
    "Usage:",
    "  Rscript plot_compare_graphs.R --prefix <prefix> [options]",
    "",
    "Required:",
    "  --prefix <prefix>  compare_graphs.py's output prefix (<prefix>.graphs.tsv, <prefix>.paths.tsv)",
    "",
    "Optional:",
    "  --out <path>       output image (default <prefix>.png)",
    "  --dpi <N>          resolution (default 150)",
    "  --width <in>, --height <in>",
    "                     size in inches (default 10 x 9)",
    "  -h, --help         show this help"), collapse = "\n"), "\n")
  quit(status = status)
}
if (length(args) == 0 || any(args %in% c("-h", "--help"))) usage(0)
get <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i == length(args)) stop("missing value after ", flag)
  args[i + 1]
}
prefix <- get("--prefix")
if (is.null(prefix)) usage(1)
out <- get("--out", paste0(prefix, ".png"))
dpi <- as.numeric(get("--dpi", "150"))
width <- as.numeric(get("--width", "10"))
height <- as.numeric(get("--height", "9"))

read <- function(p) read.delim(p, colClasses = "character", na.strings = "NA", check.names = FALSE)
g <- read(paste0(prefix, ".graphs.tsv"))
p <- read(paste0(prefix, ".paths.tsv"))
if (is.null(g$locus)) g$locus <- "locus"
if (is.null(p$locus)) p$locus <- "locus"
stages <- unique(g$graph)
loci <- unique(g$locus)
several <- length(loci) > 1   # then the loci are the x axis and the stages sit side by side

# One value per graph: drawn as bars.
scalar <- function(col, name, scale = 1) {
  data.frame(locus = g$locus, graph = g$graph, measure = name, value = as.numeric(g[[col]]) * scale)
}
bars <- rbind(
  scalar("nodes", "nodes"),
  scalar("edges", "edges"),
  scalar("sequence_bp", "sequence (Mb)", 1e-6),
  scalar("single_base_nodes", "1 bp nodes (%)", 100),
  scalar("hubs", "hubs"),
  scalar("max_degree", "max degree"),
  scalar("nodes_per_kb", "nodes/kb"),
  scalar("graph_over_median_haplotype", "graph bp / median haplotype bp"),
  scalar("self_loops", "self-loops"),
  data.frame(locus = g$locus, graph = g$graph, measure = "paths as in first (%)",
             value = 100 * as.numeric(g$paths_same_sequence) / as.numeric(g$paths_in_first)))

# One value per path: drawn as boxes.
path_measure <- function(values, name) data.frame(locus = p$locus, graph = p$graph, measure = name, value = values)
boxes <- rbind(
  path_measure(as.numeric(p$revisits), "revisits per path"),
  path_measure(as.numeric(p$steps) / (as.numeric(p$length_bp) / 1000), "steps per kb"),
  path_measure(as.numeric(p$length_diff_bp), "length vs first (bp)"))
unmatched <- sum(is.na(suppressWarnings(as.numeric(p$length_diff_bp))))
if (unmatched > 0) message(unmatched, " paths have no counterpart in their locus's first graph; ",
                           "they are left out of 'length vs first (bp)'")
boxes <- boxes[is.finite(boxes$value), ]

levels_measure <- c(unique(bars$measure), unique(boxes$measure))
for (v in c("bars", "boxes")) {
  d <- base::get(v)
  d$measure <- factor(d$measure, levels_measure)
  d$graph <- factor(d$graph, stages)
  d$locus <- factor(d$locus, loci)
  d$x <- if (several) d$locus else d$graph
  assign(v, d)
}
# A measure that is zero everywhere still gets an axis from 0 to 1, not one centred on 0.
zero <- tapply(bars$value, as.character(bars$measure), function(v) all(is.na(v) | v == 0))
floor_rows <- bars[0, ]
if (any(zero)) {
  floor_rows <- data.frame(locus = factor(loci[1], loci), graph = factor(stages[1], stages),
                           measure = factor(names(zero)[zero], levels_measure), value = 1)
  floor_rows$x <- if (several) floor_rows$locus else floor_rows$graph
}

dodge <- position_dodge(width = 0.75)
fig <- ggplot(mapping = aes(x, value, fill = graph)) +
  geom_col(data = bars, width = if (several) 0.7 else 0.6, position = dodge) +
  geom_blank(data = floor_rows) +
  expand_limits(y = 0) +
  geom_boxplot(data = boxes, width = if (several) 0.7 else 0.6, outlier.size = 0.6, position = dodge) +
  facet_wrap(~measure, scales = "free_y", ncol = 4) +
  scale_fill_brewer(palette = "Set2") +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "top", panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = if (several) element_text(angle = 45, hjust = 1) else element_text())

ggsave(out, fig, width = width, height = height, dpi = dpi)
cat("wrote", out, "\n")
