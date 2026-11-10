#!/usr/bin/env Rscript
# One inspect site as a heatmap: a row per path (per crossing), a column per node or link. See --help for usage.

usage <- function(status = 0) {
  cat(paste(c(
    "plot_inspect.R - one site of panvar inspect as a heatmap",
    "",
    "Usage:",
    "  Rscript plot_inspect.R --table <site folder or file> [options]",
    "",
    "Required:",
    "  --table <path>          a site's folder (<prefix>/bubble_<N>/) or any of its files",
    "",
    "Optional:",
    "  --features <kind>       nodes or edges (default nodes)",
    "  --out <prefix>          output prefix, writes <prefix>.png and <prefix>.pdf",
    "                          (default: in the site's folder, <name>.bubble_<N>.<features>)",
    "  --value <mode>          nodes only: total, forward or reverse passes (default total)",
    "  --transform <mode>      raw or log1p (default raw)",
    "  --length-scale <mode>   nodes only: column width from node length, sqrt, raw, log1p or none",
    "                          (default sqrt)",
    "  --representatives-only  with clusters: one row per cluster, its representative path",
    "  --no-clusters           do not group the rows by cluster",
    "  --max-paths <N>, --max-features <N>",
    "                          keep the N rows or columns with the most passes (default all)",
    "  --width <in>, --height <in>",
    "                          figure size (default from the table's size)",
    "  --dpi <N>               PNG resolution (default 300)",
    "  -h, --help              show this help"), collapse = "\n"), "\n")
  quit(status = status)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0 || any(args %in% c("-h", "--help"))) usage(0)
flags <- c("--representatives-only", "--no-clusters")
valued <- c("--table", "--features", "--out", "--value", "--transform", "--length-scale", "--max-paths",
            "--max-features", "--width", "--height", "--dpi")
opt <- list()
i <- 1
while (i <= length(args)) {
  a <- args[[i]]
  if (a %in% flags) { opt[[a]] <- TRUE; i <- i + 1; next }
  if (!a %in% valued) stop("unknown option: ", a, " (see --help)", call. = FALSE)
  if (i == length(args)) stop("missing value after ", a, call. = FALSE)
  opt[[a]] <- args[[i + 1]]
  i <- i + 2
}
get_opt <- function(flag, default = NULL) if (is.null(opt[[flag]])) default else opt[[flag]]
one_of <- function(flag, default, allowed) {
  v <- get_opt(flag, default)
  if (!v %in% allowed) stop(flag, " must be one of: ", paste(allowed, collapse = ", "), call. = FALSE)
  v
}
count_opt <- function(flag) {
  v <- suppressWarnings(as.integer(get_opt(flag, "0")))
  if (is.na(v) || v < 0) stop(flag, " must be a non-negative integer", call. = FALSE)
  v
}
table_arg <- get_opt("--table")
if (is.null(table_arg)) usage(1)
features <- one_of("--features", "nodes", c("nodes", "edges"))
value <- one_of("--value", "total", c("total", "forward", "reverse"))
transform <- one_of("--transform", "raw", c("raw", "log1p"))
length_scale <- one_of("--length-scale", "sqrt", c("sqrt", "raw", "log1p", "none"))
max_paths <- count_opt("--max-paths")
max_features <- count_opt("--max-features")
if (features == "edges" && !is.null(opt[["--value"]]))
  stop("--value applies to nodes: a link's passes are not split by direction", call. = FALSE)

if (!requireNamespace("ggplot2", quietly = TRUE))
  stop("plot_inspect.R needs ggplot2 (conda install -c conda-forge r-ggplot2)", call. = FALSE)

# The site's prefix, from any of its files: everything before inspect's suffix.
suffixes <- c(".node_counts.tsv", ".edge_counts.tsv", ".node_lengths.tsv", ".clusters.tsv", ".paths.fa.gz",
              ".walks.fa.gz", ".alleles.fa.gz", ".alleles.walks.fa.gz", ".nodes.json.gz", ".genotypes.tsv.gz")
# A site folder (<prefix>/bubble_<N>/) holds <name>.bubble_<N>.<suffix>: the stem is read off its files.
table_arg <- sub("/+$", "", table_arg)
if (!file.exists(table_arg)) stop("not found: ", table_arg, call. = FALSE)
stem <- table_arg
if (dir.exists(table_arg)) {
  found <- list.files(table_arg, pattern = "\\.bubble_[0-9]+\\.(node|edge)_counts\\.tsv(\\.gz)?$", full.names = TRUE)
  if (length(found) == 0) stop("no inspect count table in ", table_arg, call. = FALSE)
  stem <- sub("\\.(node|edge)_counts\\.tsv(\\.gz)?$", "", found[[1]])
}
for (s in suffixes) {
  for (z in c("", ".gz")) {
    if (endsWith(table_arg, paste0(s, z))) stem <- substr(table_arg, 1, nchar(table_arg) - nchar(paste0(s, z)))
  }
}
find_file <- function(suffix) {
  for (p in c(paste0(stem, suffix), paste0(stem, suffix, ".gz"))) if (file.exists(p)) return(p)
  NULL
}
counts_path <- find_file(if (features == "nodes") ".node_counts.tsv" else ".edge_counts.tsv")
if (is.null(counts_path))
  stop("no ", if (features == "nodes") "node" else "edge", " count table beside ", table_arg, " (looked for ",
       stem, if (features == "nodes") ".node_counts.tsv" else ".edge_counts.tsv", ")", call. = FALSE)
lengths_path <- if (features == "nodes" && length_scale != "none") find_file(".node_lengths.tsv") else NULL
clusters_path <- if (is.null(opt[["--no-clusters"]])) find_file(".clusters.tsv") else NULL
if (!is.null(opt[["--representatives-only"]]) && is.null(clusters_path))
  stop("--representatives-only needs ", stem, ".clusters.tsv (inspect --cluster)", call. = FALSE)
out <- get_opt("--out", paste0(stem, ".", features))

# Every column as text, so a path or node named 001 keeps its name; numbers are converted where used.
read_tsv <- function(path) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con))
  read.delim(con, header = TRUE, check.names = FALSE, quote = "", comment.char = "", colClasses = "character",
             na.strings = character(0))
}

tab <- read_tsv(counts_path)
prefix_col <- if (features == "nodes") "^node\\." else "^edge\\."
cols <- grep(prefix_col, names(tab), value = TRUE)
if (!"path_name" %in% names(tab) || length(cols) == 0 || nrow(tab) == 0)
  stop(counts_path, " has no rows or no ", if (features == "nodes") "node" else "edge", " columns", call. = FALSE)
cells <- as.character(unlist(tab[cols], use.names = FALSE))
field <- if (features == "edges") 1 else switch(value, total = 1, forward = 2, reverse = 3)
values <- suppressWarnings(as.numeric(vapply(strsplit(cells, ":", fixed = TRUE),
                                             function(x) if (length(x) >= field) x[[field]] else NA_character_, "")))
values[is.na(values)] <- 0
mat <- matrix(values, nrow = nrow(tab), dimnames = list(as.character(tab$path_name), sub(prefix_col, "", cols)))

# Rows grouped by cluster, each crossing once, the representative first; rows in no cluster last.
boundaries <- integer(0)
if (!is.null(clusters_path)) {
  cl <- read_tsv(clusters_path)
  cluster_of <- setNames(rep(NA_integer_, nrow(mat)), rownames(mat))
  representative <- setNames(rep(FALSE, nrow(mat)), rownames(mat))
  for (r in seq_len(nrow(cl))) {
    members <- strsplit(cl$members[r], ";", fixed = TRUE)[[1]]
    members <- members[members %in% names(cluster_of)]
    cluster_of[members] <- as.integer(cl$cluster_id[r])
    if (cl$representative_path[r] %in% names(representative)) representative[cl$representative_path[r]] <- TRUE
  }
  if (!is.null(opt[["--representatives-only"]])) {
    mat <- mat[representative, , drop = FALSE]
    cluster_of <- cluster_of[representative]
    representative <- representative[representative]
  }
  key <- ifelse(is.na(cluster_of), .Machine$integer.max, cluster_of)
  ord <- order(key, !representative, rownames(mat))
  mat <- mat[ord, , drop = FALSE]
  key <- key[ord]
  if (length(key) > 1) boundaries <- which(key[-1] != key[-length(key)]) + 1
}
if (max_paths > 0 && nrow(mat) > max_paths) {
  keep <- sort(order(rowSums(mat), decreasing = TRUE)[seq_len(max_paths)])
  mat <- mat[keep, , drop = FALSE]
  boundaries <- integer(0)   # the groups no longer hold every member
}
if (max_features > 0 && ncol(mat) > max_features)
  mat <- mat[, sort(order(colSums(mat), decreasing = TRUE)[seq_len(max_features)]), drop = FALSE]
plot_mat <- if (transform == "log1p") log1p(mat) else mat

# Column extents: unit width, or with node lengths a width from each node's length.
n_rows <- nrow(plot_mat)
n_cols <- ncol(plot_mat)
widths <- rep(1, n_cols)
if (!is.null(lengths_path)) {
  nl <- read_tsv(lengths_path)
  len <- as.numeric(setNames(nl$length_bp, nl$node_id)[colnames(plot_mat)])
  widths <- switch(length_scale, raw = len, sqrt = sqrt(pmax(len, 0)), log1p = log1p(pmax(len, 0)))
  smallest <- suppressWarnings(min(widths[is.finite(widths) & widths > 0]))
  widths[!is.finite(widths) | widths <= 0] <- if (is.finite(smallest)) smallest / 4 else 1
}
x1 <- cumsum(widths)
x0 <- x1 - widths
df <- data.frame(col = rep(seq_len(n_cols), each = n_rows), row = rep(seq_len(n_rows), times = n_cols),
                 value = as.vector(plot_mat))
df$y <- n_rows - df$row + 1
ticks <- function(n, most) unique(pmax(1, pmin(n, round(seq(1, n, length.out = min(n, most))))))
xt <- ticks(n_cols, 30)
yt <- ticks(n_rows, 35)
what <- if (features == "nodes") "node" else "link"

p <- ggplot2::ggplot(df) +
  ggplot2::geom_rect(ggplot2::aes(xmin = x0[col], xmax = x1[col], ymin = y - 0.5, ymax = y + 0.5, fill = value)) +
  ggplot2::scale_fill_gradientn(colours = c("#ffffff", "#fff1ec", "#fcbba1", "#fb6a4a", "#cb181d", "#67000d"),
                                name = if (transform == "log1p") "log1p(passes)" else "passes") +
  ggplot2::scale_x_continuous(breaks = (x0 + x1)[xt] / 2, labels = colnames(plot_mat)[xt], expand = c(0, 0)) +
  ggplot2::scale_y_continuous(breaks = n_rows - yt + 1, labels = rownames(plot_mat)[yt], expand = c(0, 0)) +
  ggplot2::labs(title = basename(stem), x = paste0(what, "s, in graph order"), y = "paths") +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(panel.grid = ggplot2::element_blank(), plot.title = ggplot2::element_text(face = "bold"),
                 axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
                 axis.text.y = ggplot2::element_text(size = 5))
if (length(boundaries) > 0)
  p <- p + ggplot2::geom_hline(yintercept = n_rows - boundaries + 1.5, colour = "grey30", linewidth = 0.3)

width <- as.numeric(get_opt("--width", max(8, min(24, 4 + n_cols / 85))))
height <- as.numeric(get_opt("--height", max(7, min(24, 4 + n_rows / 32))))
dpi <- as.numeric(get_opt("--dpi", "300"))
if (!dir.exists(dirname(out))) dir.create(dirname(out), recursive = TRUE)
ggplot2::ggsave(paste0(out, ".png"), p, width = width, height = height, units = "in", dpi = dpi, limitsize = FALSE)
ggplot2::ggsave(paste0(out, ".pdf"), p, width = width, height = height, units = "in", limitsize = FALSE)
cat("wrote", paste0(out, ".png"), "and", paste0(out, ".pdf"), "\n")
