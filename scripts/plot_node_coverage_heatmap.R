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
    "  --node-lengths <path>  panvar inspect node_lengths.tsv; scales x by node bp length",
    "  --length-transform <mode>  raw, sqrt, or log1p tile-width scaling (default: sqrt)",
    "  --clusters <path>    panvar inspect clusters.tsv; keep only representative paths",
    "  --cluster-by <path>  panvar inspect clusters.tsv; keep all paths but group/order rows by",
    "                       cluster (representative first), with a separator between clusters",
    "  --max-paths <N>      Keep at most N paths, selected by total coverage (default: all)",
    "  --max-nodes <N>      Keep at most N nodes, selected by total coverage (default: all)",
    "  --width <inches>     Figure width (default: auto)",
    "  --height <inches>    Figure height (default: auto)",
    "  --dpi <int>          PNG resolution (default: 300)",
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
  node_lengths = NULL,
  length_transform = "sqrt",
  clusters = NULL,
  cluster_by = NULL,
  max_paths = 0,
  max_nodes = 0,
  width = NA_real_,
  height = NA_real_,
  dpi = 300
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
  } else if (arg == "--node-lengths") {
    opts$node_lengths <- read_value(arg)
    i <- i + 2
  } else if (arg == "--length-transform") {
    opts$length_transform <- read_value(arg)
    i <- i + 2
  } else if (arg == "--clusters") {
    opts$clusters <- read_value(arg)
    i <- i + 2
  } else if (arg == "--cluster-by") {
    opts$cluster_by <- read_value(arg)
    i <- i + 2
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
  } else if (arg == "--dpi") {
    opts$dpi <- as.numeric(read_value(arg))
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
if (!opts$length_transform %in% c("raw", "sqrt", "log1p")) {
  stop("--length-transform must be one of: raw, sqrt, log1p", call. = FALSE)
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

# Read a TSV and close its connection (avoids leaked-connection warnings).
read_tsv <- function(path) {
  con <- open_input(path)
  on.exit(close(con))
  read.delim(con, sep = "\t", header = TRUE, check.names = FALSE, quote = "", comment.char = "")
}

# Order rows by the inspect clusters.tsv (cluster_id, rep first, then name; unlisted rows last).
# Returns the reordered matrix + the row indices where a new cluster starts (for separators).
order_rows_by_clusters <- function(mat, clusters_path) {
  cl <- read_tsv(clusters_path)
  if (!all(c("cluster_id", "members") %in% names(cl))) {
    stop("--cluster-by table must contain cluster_id and members columns", call. = FALSE)
  }
  has_rep <- "representative_path" %in% names(cl)
  path_cluster <- new.env(parent = emptyenv())
  path_isrep <- new.env(parent = emptyenv())
  for (r in seq_len(nrow(cl))) {
    cid <- suppressWarnings(as.integer(cl$cluster_id[r]))
    mem <- strsplit(as.character(cl$members[r]), ";", fixed = TRUE)[[1]]
    mem <- mem[nzchar(mem)]
    repp <- if (has_rep) as.character(cl$representative_path[r]) else NA_character_
    for (m in mem) {
      assign(m, cid, envir = path_cluster)
      assign(m, identical(m, repp), envir = path_isrep)
    }
  }
  rn <- rownames(mat)
  cid_vec <- vapply(rn, function(x)
    if (exists(x, envir = path_cluster, inherits = FALSE)) get(x, envir = path_cluster) else NA_integer_,
    integer(1))
  isrep_vec <- vapply(rn, function(x)
    exists(x, envir = path_isrep, inherits = FALSE) && isTRUE(get(x, envir = path_isrep)),
    logical(1))
  cid_sort <- ifelse(is.na(cid_vec), .Machine$integer.max, cid_vec)
  ord <- order(cid_sort, !isrep_vec, rn)
  cid_ordered <- cid_sort[ord]
  boundaries <- if (length(cid_ordered) > 1) {
    which(cid_ordered[-1] != cid_ordered[-length(cid_ordered)]) + 1
  } else {
    integer(0)
  }
  list(mat = mat[ord, , drop = FALSE], boundaries = boundaries)
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

if (!is.null(opts$clusters)) {
  cl <- read_tsv(opts$clusters)
  if (!"representative_path" %in% names(cl)) {
    stop("--clusters table must contain a representative_path column", call. = FALSE)
  }
  reps <- unique(as.character(cl$representative_path))
  tab <- tab[as.character(tab$path_name) %in% reps, , drop = FALSE]
  if (nrow(tab) == 0) {
    stop("No path rows match the cluster representatives", call. = FALSE)
  }
}

len_by_node <- NULL
if (!is.null(opts$node_lengths)) {
  nl <- read_tsv(opts$node_lengths)
  if (!all(c("node_id", "length_bp") %in% names(nl))) {
    stop("--node-lengths table must contain node_id and length_bp columns", call. = FALSE)
  }
  len_by_node <- stats::setNames(as.numeric(nl$length_bp), as.character(nl$node_id))
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

# Row order: grouped by the external inspect clusters file when given (no auto-clustering).
cluster_boundaries <- integer(0)
if (!is.null(opts$cluster_by)) {
  grouped <- order_rows_by_clusters(mat, opts$cluster_by)
  ord <- match(rownames(grouped$mat), rownames(mat))
  mat <- grouped$mat
  plot_mat <- plot_mat[ord, , drop = FALSE]
  cluster_boundaries <- grouped$boundaries
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

n_paths <- nrow(plot_mat)
n_nodes <- ncol(plot_mat)
use_length <- !is.null(len_by_node)

# Per-column x extents. Default: unit-width columns (equal-spaced, like geom_raster).
# With --node-lengths: tile widths = (transformed) node bp, preserving column order.
if (use_length) {
  raw_len <- as.numeric(len_by_node[colnames(plot_mat)])
  w <- switch(opts$length_transform,
              raw = raw_len,
              sqrt = sqrt(pmax(raw_len, 0)),
              log1p = log1p(pmax(raw_len, 0)))
  pos_min <- suppressWarnings(min(w[is.finite(w) & w > 0]))
  if (!is.finite(pos_min)) pos_min <- 1
  w[!is.finite(w) | w <= 0] <- pos_min * 0.25
  x1 <- cumsum(w)
  x0 <- x1 - w
  xmid <- (x0 + x1) / 2
} else {
  x0 <- seq_len(n_nodes) - 1
  x1 <- seq_len(n_nodes)
  xmid <- seq_len(n_nodes) - 0.5
}

path_pos <- seq_len(n_paths)
plot_df <- data.frame(
  node_idx = rep(seq_len(n_nodes), each = n_paths),
  path_pos = rep(path_pos, times = n_nodes),
  value = as.vector(plot_mat)
)
plot_df$path_y <- n_paths - plot_df$path_pos + 1
plot_df$xmin <- x0[plot_df$node_idx]
plot_df$xmax <- x1[plot_df$node_idx]
plot_df$ymin <- plot_df$path_y - 0.5
plot_df$ymax <- plot_df$path_y + 0.5

x_break_idx <- sparse_ticks(n_nodes, 30)
x_breaks <- if (use_length) xmid[x_break_idx] else x_break_idx
x_labels <- colnames(plot_mat)[x_break_idx]
y_breaks_original <- sparse_ticks(n_paths, 35)
y_breaks <- n_paths - y_breaks_original + 1
fill_name <- if (opts$transform == "log1p") "log1p(count)" else "count"
pal <- c("#ffffff", "#fff1ec", "#fcbba1", "#fb6a4a", "#cb181d", "#67000d")

if (use_length) {
  p <- ggplot2::ggplot(plot_df) +
    ggplot2::geom_rect(ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = value))
} else {
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = node_idx, y = path_y, fill = value)) +
    ggplot2::geom_raster()
}

p <- p +
  ggplot2::scale_fill_gradientn(colours = pal, name = fill_name) +
  ggplot2::scale_x_continuous(
    breaks = x_breaks,
    labels = x_labels,
    expand = c(0, 0)
  ) +
  ggplot2::scale_y_continuous(
    breaks = y_breaks,
    labels = rownames(plot_mat)[y_breaks_original],
    expand = c(0, 0)
  ) +
  ggplot2::labs(
    title = sprintf("Node %s coverage", opts$value),
    subtitle = sprintf("%d paths x %d bubble-internal nodes%s", n_paths, n_nodes,
                       if (use_length) " (x scaled by node bp)" else ""),
    x = if (use_length) "Bubble-internal nodes (width = bp)" else "Bubble-internal nodes",
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

# Separator lines between adjacent inspect clusters (only with --cluster-by).
if (length(cluster_boundaries) > 0) {
  p <- p + ggplot2::geom_hline(
    yintercept = n_paths - cluster_boundaries + 1.5,
    colour = "grey30", linewidth = 0.3
  )
}

png_path <- paste0(opts$out, ".png")
pdf_path <- paste0(opts$out, ".pdf")

ggplot2::ggsave(png_path, p, width = opts$width, height = opts$height, units = "in", dpi = opts$dpi, limitsize = FALSE)
ggplot2::ggsave(pdf_path, p, width = opts$width, height = opts$height, units = "in", limitsize = FALSE)

cat("Wrote:", png_path, "\n")
cat("Wrote:", pdf_path, "\n")
