#!/usr/bin/env Rscript

usage <- function(status = 0) {
  msg <- c(
    "Usage:",
    "  plot_edge_coverage_heatmap.R --table <edge_counts.tsv[.gz]> --out <output_prefix> [options]",
    "  plot_edge_coverage_heatmap.R <edge_counts.tsv[.gz]> <output_prefix>",
    "",
    "Options:",
    "  --table <path>       panvar inspect edge-count table",
    "  --out <prefix>       Output prefix; writes <prefix>.png and <prefix>.pdf",
    "  --transform <mode>   raw or log1p (default: raw)",
    "  --clusters <path>    panvar inspect clusters.tsv; keep only representative paths",
    "  --cluster-by <path>  panvar inspect clusters.tsv; keep all paths but group/order rows by",
    "                       cluster (representative first), with a separator between clusters",
    "  --max-paths <N>      Keep at most N paths, selected by total traversals (default: all)",
    "  --max-edges <N>      Keep at most N edges, selected by total traversals (default: all)",
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
  transform = "raw",
  clusters = NULL,
  cluster_by = NULL,
  max_paths = 0,
  max_edges = 0,
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
  } else if (arg == "--transform") {
    opts$transform <- read_value(arg)
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
  } else if (arg == "--max-edges") {
    opts$max_edges <- as.integer(read_value(arg))
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
if (!opts$transform %in% c("raw", "log1p")) {
  stop("--transform must be one of: raw, log1p", call. = FALSE)
}
if (is.na(opts$max_paths) || opts$max_paths < 0 || is.na(opts$max_edges) || opts$max_edges < 0) {
  stop("--max-paths/--max-edges must be non-negative integers", call. = FALSE)
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

# Reorder a coverage matrix so cluster-mates (from a panvar inspect clusters.tsv) are
# adjacent: clusters ascending by cluster_id, representative first within each cluster,
# then path name; any rows absent from the cluster file sort last. Returns the reordered
# matrix and the mat-row indices (>1) where a new cluster begins (for separator lines).
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

edge_cols <- grep("^edge\\.", names(tab), value = TRUE)
if (length(edge_cols) == 0) {
  stop("Input table has no edge.* columns", call. = FALSE)
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

values <- suppressWarnings(as.numeric(as.character(unlist(tab[edge_cols], use.names = FALSE))))
values[is.na(values)] <- 0
mat <- matrix(values, nrow = nrow(tab), ncol = length(edge_cols), byrow = FALSE)
rownames(mat) <- make.unique(as.character(tab$path_name), sep = ".")
colnames(mat) <- sub("^edge\\.", "", edge_cols)

if (opts$max_paths > 0 && nrow(mat) > opts$max_paths) {
  keep <- order(rowSums(mat), decreasing = TRUE)[seq_len(opts$max_paths)]
  keep <- sort(keep)
  mat <- mat[keep, , drop = FALSE]
}
if (opts$max_edges > 0 && ncol(mat) > opts$max_edges) {
  keep <- order(colSums(mat), decreasing = TRUE)[seq_len(opts$max_edges)]
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

path_pos <- seq_len(nrow(plot_mat))
edge_pos <- seq_len(ncol(plot_mat))
plot_df <- data.frame(
  edge_pos = rep(edge_pos, each = nrow(plot_mat)),
  path_pos = rep(path_pos, times = ncol(plot_mat)),
  value = as.vector(plot_mat)
)
plot_df$path_y <- nrow(plot_mat) - plot_df$path_pos + 1

x_breaks <- sparse_ticks(ncol(plot_mat), 30)
y_breaks_original <- sparse_ticks(nrow(plot_mat), 35)
y_breaks <- nrow(plot_mat) - y_breaks_original + 1
fill_name <- if (opts$transform == "log1p") "log1p(traversals)" else "traversals"

p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = edge_pos, y = path_y, fill = value)) +
  ggplot2::geom_raster() +
  ggplot2::scale_fill_gradientn(
    colours = c("#ffffff", "#eef3fb", "#c6dbef", "#6baed6", "#2171b5", "#08306b"),
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
    title = "Edge traversal coverage",
    subtitle = sprintf("%d paths x %d bubble-internal edges", nrow(plot_mat), ncol(plot_mat)),
    x = "Bubble-internal edges (from>to, orientation-aware)",
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
    yintercept = nrow(plot_mat) - cluster_boundaries + 1.5,
    colour = "grey30", linewidth = 0.3
  )
}

png_path <- paste0(opts$out, ".png")
pdf_path <- paste0(opts$out, ".pdf")

ggplot2::ggsave(png_path, p, width = opts$width, height = opts$height, units = "in", dpi = opts$dpi, limitsize = FALSE)
ggplot2::ggsave(pdf_path, p, width = opts$width, height = opts$height, units = "in", limitsize = FALSE)

cat("Wrote:", png_path, "\n")
cat("Wrote:", pdf_path, "\n")
