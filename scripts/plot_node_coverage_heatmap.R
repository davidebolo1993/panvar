#!/usr/bin/env Rscript

usage <- function(status = 0) {
  msg <- c(
    "Usage:",
    "  plot_node_coverage_heatmap.R --table <node_counts.tsv[.gz]> --out <output_prefix> [options]",
    "  plot_node_coverage_heatmap.R <node_counts.tsv[.gz]> <output_prefix>",
    "",
    "Options:",
    "  --table <path>       panvar inspect node-count table",
    "  --out <prefix>       Output prefix; writes <prefix>.png and <prefix>.pdf",
    "  --value <mode>       total, forward, or reverse (default: total)",
    "  --transform <mode>   raw or log1p (default: raw)",
    "  --cluster-rows       Cluster paths by coverage profile",
    "  --cluster-cols       Cluster nodes by coverage profile",
    "  --max-paths <N>      Keep at most N paths, selected by total coverage (default: all)",
    "  --max-nodes <N>      Keep at most N nodes, selected by total coverage (default: all)",
    "  --width <inches>     Figure width (default: auto)",
    "  --height <inches>    Figure height (default: auto)",
    "  -h, --help           Show this help"
  )
  cat(paste(msg, collapse = "\n"), "\n", sep = "")
  quit(status = status)
}

require_ggplot2 <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "This plotting helper requires ggplot2. Install it with: conda install -y -c conda-forge r-ggplot2",
      call. = FALSE
    )
  }
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0 || any(args %in% c("-h", "--help"))) {
  usage(0)
}

opts <- list(
  table = NULL,
  out = NULL,
  value = "total",
  transform = "raw",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  max_paths = 0,
  max_nodes = 0,
  width = NA_real_,
  height = NA_real_
)
positional <- character()

i <- 1
while (i <= length(args)) {
  arg <- args[[i]]
  read_value <- function(flag) {
    if (i + 1 > length(args)) {
      stop(paste("Missing value after", flag), call. = FALSE)
    }
    args[[i + 1]]
  }
  if (arg == "--table") {
    opts$table <- read_value(arg)
    i <- i + 2
  } else if (arg == "--out") {
    opts$out <- read_value(arg)
    i <- i + 2
  } else if (arg == "--value") {
    opts$value <- read_value(arg)
    i <- i + 2
  } else if (arg == "--transform") {
    opts$transform <- read_value(arg)
    i <- i + 2
  } else if (arg == "--cluster-rows") {
    opts$cluster_rows <- TRUE
    i <- i + 1
  } else if (arg == "--cluster-cols") {
    opts$cluster_cols <- TRUE
    i <- i + 1
  } else if (arg == "--max-paths") {
    opts$max_paths <- as.integer(read_value(arg))
    i <- i + 2
  } else if (arg == "--max-nodes") {
    opts$max_nodes <- as.integer(read_value(arg))
    i <- i + 2
  } else if (arg == "--width") {
    opts$width <- as.numeric(read_value(arg))
    i <- i + 2
  } else if (arg == "--height") {
    opts$height <- as.numeric(read_value(arg))
    i <- i + 2
  } else if (startsWith(arg, "-")) {
    stop(paste("Unknown option:", arg), call. = FALSE)
  } else {
    positional <- c(positional, arg)
    i <- i + 1
  }
}

if (is.null(opts$table) && length(positional) >= 1) {
  opts$table <- positional[[1]]
}
if (is.null(opts$out) && length(positional) >= 2) {
  opts$out <- positional[[2]]
}
if (is.null(opts$table) || is.null(opts$out)) {
  usage(1)
}
if (!opts$value %in% c("total", "forward", "reverse")) {
  stop("--value must be one of: total, forward, reverse", call. = FALSE)
}
if (!opts$transform %in% c("raw", "log1p")) {
  stop("--transform must be one of: raw, log1p", call. = FALSE)
}
if (is.na(opts$max_paths) || opts$max_paths < 0 || is.na(opts$max_nodes) || opts$max_nodes < 0) {
  stop("--max-paths/--max-nodes must be non-negative integers", call. = FALSE)
}

require_ggplot2()

open_input <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    gzfile(path, open = "rt")
  } else {
    file(path, open = "rt")
  }
}

con <- open_input(opts$table)
on.exit(close(con), add = TRUE)
tab <- read.delim(con, sep = "\t", header = TRUE, check.names = FALSE, quote = "", comment.char = "")
if (!all(c("path_name", "path_length_bp") %in% names(tab))) {
  stop("Input table must contain path_name and path_length_bp columns", call. = FALSE)
}

node_cols <- grep("^node\\.", names(tab), value = TRUE)
if (length(node_cols) == 0) {
  stop("Input table has no node.* columns", call. = FALSE)
}
if (nrow(tab) == 0) {
  stop("Input table has no path rows", call. = FALSE)
}

cells <- as.character(unlist(tab[node_cols], use.names = FALSE))
if (opts$value == "total") {
  values <- suppressWarnings(as.numeric(sub(":.*$", "", cells)))
} else if (opts$value == "forward") {
  values <- suppressWarnings(as.numeric(sub("^[^:]*:([^:]*):.*$", "\\1", cells)))
} else {
  values <- suppressWarnings(as.numeric(sub("^.*:", "", cells)))
}
values[is.na(values)] <- 0
mat <- matrix(values, nrow = nrow(tab), ncol = length(node_cols), byrow = FALSE)
rownames(mat) <- make.unique(as.character(tab$path_name), sep = ".")
colnames(mat) <- sub("^node\\.", "", node_cols)

if (opts$max_paths > 0 && nrow(mat) > opts$max_paths) {
  keep <- order(rowSums(mat), decreasing = TRUE)[seq_len(opts$max_paths)]
  keep <- sort(keep)
  mat <- mat[keep, , drop = FALSE]
}
if (opts$max_nodes > 0 && ncol(mat) > opts$max_nodes) {
  keep <- order(colSums(mat), decreasing = TRUE)[seq_len(opts$max_nodes)]
  keep <- sort(keep)
  mat <- mat[, keep, drop = FALSE]
}

plot_mat <- mat
if (opts$transform == "log1p") {
  plot_mat <- log1p(plot_mat)
}

if (opts$cluster_rows && nrow(plot_mat) > 1) {
  row_order <- hclust(dist(plot_mat))$order
  mat <- mat[row_order, , drop = FALSE]
  plot_mat <- plot_mat[row_order, , drop = FALSE]
}
if (opts$cluster_cols && ncol(plot_mat) > 1) {
  col_order <- hclust(dist(t(plot_mat)))$order
  mat <- mat[, col_order, drop = FALSE]
  plot_mat <- plot_mat[, col_order, drop = FALSE]
}

out_dir <- dirname(opts$out)
if (!identical(out_dir, ".") && !dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

if (is.na(opts$width)) {
  opts$width <- max(8, min(24, 4 + ncol(plot_mat) / 85))
}
if (is.na(opts$height)) {
  opts$height <- max(7, min(24, 4 + nrow(plot_mat) / 32))
}

sparse_ticks <- function(n, max_ticks) {
  if (n <= 0) {
    integer()
  } else {
    unique(pmax(1, pmin(n, round(seq(1, n, length.out = min(n, max_ticks))))))
  }
}

path_pos <- seq_len(nrow(plot_mat))
node_pos <- seq_len(ncol(plot_mat))
plot_df <- data.frame(
  node_pos = rep(node_pos, each = nrow(plot_mat)),
  path_pos = rep(path_pos, times = ncol(plot_mat)),
  value = as.vector(plot_mat)
)
plot_df$path_y <- nrow(plot_mat) - plot_df$path_pos + 1

x_breaks <- sparse_ticks(ncol(plot_mat), 30)
y_breaks_original <- sparse_ticks(nrow(plot_mat), 35)
y_breaks <- nrow(plot_mat) - y_breaks_original + 1
fill_name <- if (opts$transform == "log1p") "log1p(count)" else "count"

p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = node_pos, y = path_y, fill = value)) +
  ggplot2::geom_raster() +
  ggplot2::scale_fill_gradientn(
    colours = c("#ffffff", "#fff1ec", "#fcbba1", "#fb6a4a", "#cb181d", "#67000d"),
    name = fill_name
  ) +
  ggplot2::scale_x_continuous(
    breaks = x_breaks,
    labels = colnames(plot_mat)[x_breaks],
    expand = c(0, 0)
  ) +
  ggplot2::scale_y_continuous(
    breaks = y_breaks,
    labels = rownames(plot_mat)[y_breaks_original],
    expand = c(0, 0)
  ) +
  ggplot2::labs(
    title = sprintf("Node %s coverage", opts$value),
    subtitle = sprintf("%d paths x %d bubble-internal nodes", nrow(plot_mat), ncol(plot_mat)),
    x = "Bubble-internal nodes",
    y = "Paths"
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(face = "bold"),
    axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
    axis.text.y = ggplot2::element_text(size = 5),
    legend.position = "right"
  )

png_path <- paste0(opts$out, ".png")
pdf_path <- paste0(opts$out, ".pdf")

ggplot2::ggsave(png_path, p, width = opts$width, height = opts$height, units = "in", dpi = 180, limitsize = FALSE)
ggplot2::ggsave(pdf_path, p, width = opts$width, height = opts$height, units = "in", limitsize = FALSE)

cat("Wrote:", png_path, "\n")
cat("Wrote:", pdf_path, "\n")
