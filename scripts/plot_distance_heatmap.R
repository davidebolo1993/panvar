#!/usr/bin/env Rscript

usage <- function(exit_code = 0) {
  msg <- paste(
    "Usage:",
    "  plot_distance_heatmap.R --matrix <distance_matrix_norm.tsv> [--alleles <alleles.tsv>] --out <output_prefix>",
    "",
    "Positional shorthand:",
    "  plot_distance_heatmap.R <distance_matrix_norm.tsv> [alleles.tsv] <output_prefix>",
    "",
    "If --alleles is omitted, the script looks for alleles.tsv in the same directory as the matrix.",
    "Outputs: <output_prefix>.png and <output_prefix>.pdf",
    sep = "\n"
  )
  writeLines(msg, con = if (exit_code == 0) stdout() else stderr())
  quit(save = "no", status = exit_code)
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

matrix_path <- NULL
alleles_path <- NULL
out_prefix <- NULL

if (any(grepl("^--", args))) {
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (key %in% c("--matrix", "-m")) {
      i <- i + 1
      if (i > length(args)) stop("Missing value after ", key, call. = FALSE)
      matrix_path <- args[[i]]
    } else if (key %in% c("--alleles", "-a")) {
      i <- i + 1
      if (i > length(args)) stop("Missing value after ", key, call. = FALSE)
      alleles_path <- args[[i]]
    } else if (key %in% c("--out", "-o")) {
      i <- i + 1
      if (i > length(args)) stop("Missing value after ", key, call. = FALSE)
      out_prefix <- args[[i]]
    } else {
      stop("Unknown option: ", key, call. = FALSE)
    }
    i <- i + 1
  }
} else if (length(args) == 2) {
  matrix_path <- args[[1]]
  out_prefix <- args[[2]]
} else if (length(args) == 3) {
  matrix_path <- args[[1]]
  alleles_path <- args[[2]]
  out_prefix <- args[[3]]
} else {
  usage(1)
}

if (is.null(matrix_path) || is.null(out_prefix)) {
  usage(1)
}
if (!file.exists(matrix_path)) {
  stop("Distance matrix not found: ", matrix_path, call. = FALSE)
}

if (is.null(alleles_path)) {
  candidate <- file.path(dirname(matrix_path), "alleles.tsv")
  if (file.exists(candidate)) {
    alleles_path <- candidate
  }
}
if (!is.null(alleles_path) && !file.exists(alleles_path)) {
  stop("Alleles table not found: ", alleles_path, call. = FALSE)
}

require_ggplot2()

mat <- read.table(matrix_path, header = TRUE, sep = "\t", row.names = 1,
                  check.names = FALSE, quote = "", comment.char = "")
dist_mat <- as.matrix(mat)
storage.mode(dist_mat) <- "numeric"
dist_mat <- (dist_mat + t(dist_mat)) / 2
diag(dist_mat) <- 0

n <- nrow(dist_mat)
if (n < 2) {
  stop("Need at least two alleles for a heatmap", call. = FALSE)
}

ids <- rownames(dist_mat)
labels_by_id <- setNames(paste0("A", ids), ids)
cluster_by_id <- NULL
cluster_colors <- NULL

if (!is.null(alleles_path)) {
  alleles <- read.table(alleles_path, header = TRUE, sep = "\t",
                        check.names = FALSE, quote = "", comment.char = "")
  alleles$allele_id <- as.character(alleles$allele_id)
  alleles <- alleles[match(ids, alleles$allele_id), ]
  if (any(is.na(alleles$cluster_id))) {
    stop("Could not match all matrix allele IDs to alleles.tsv", call. = FALSE)
  }

  labels_by_id <- setNames(
    paste0("A", alleles$allele_id, " C", alleles$cluster_id, " n=", alleles$path_support),
    alleles$allele_id
  )
  clusters <- sort(unique(alleles$cluster_id))
  palette <- c(
    "#4E79A7", "#F28E2B", "#59A14F", "#E15759", "#76B7B2", "#EDC948",
    "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC", "#2F4B7C", "#A05195"
  )
  cluster_colors <- setNames(palette[((seq_along(clusters) - 1) %% length(palette)) + 1],
                             as.character(clusters))
  cluster_by_id <- setNames(as.character(alleles$cluster_id), alleles$allele_id)
}

hc <- hclust(as.dist(dist_mat), method = "average")
ordered_ids <- ids[hc$order]
ordered_mat <- dist_mat[ordered_ids, ordered_ids, drop = FALSE]

axis_labels <- labels_by_id[ordered_ids]
if (n > 80) {
  axis_labels <- rep("", n)
}

row_pos <- seq_len(n)
col_pos <- seq_len(n)
plot_df <- data.frame(
  col_pos = rep(col_pos, each = n),
  row_pos = rep(row_pos, times = n),
  distance = as.vector(ordered_mat)
)
plot_df$row_y <- n - plot_df$row_pos + 1

sparse_ticks <- function(n, max_ticks) {
  if (n <= 0) {
    integer()
  } else {
    unique(pmax(1, pmin(n, round(seq(1, n, length.out = min(n, max_ticks))))))
  }
}

max_axis_ticks <- if (n <= 40) n else if (n <= 100) 40 else 25
x_breaks <- sparse_ticks(n, max_axis_ticks)
y_breaks_original <- sparse_ticks(n, max_axis_ticks)
y_breaks <- n - y_breaks_original + 1

p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = col_pos, y = row_y, fill = distance)) +
  ggplot2::geom_raster() +
  ggplot2::scale_fill_gradientn(
    colours = c("#102A43", "#2F80ED", "#F7F7F2", "#F2994A", "#8C1D18"),
    limits = c(0, max(ordered_mat, na.rm = TRUE)),
    name = "normalized\ndistance"
  ) +
  ggplot2::scale_x_continuous(
    breaks = x_breaks,
    labels = axis_labels[x_breaks],
    expand = c(0, 0)
  ) +
  ggplot2::scale_y_continuous(
    breaks = y_breaks,
    labels = axis_labels[y_breaks_original],
    expand = c(0, 0)
  ) +
  ggplot2::labs(
    title = "Allele distance heatmap",
    subtitle = sprintf("%d alleles, average-linkage order", n),
    x = "Alleles",
    y = "Alleles"
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(face = "bold"),
    axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
    axis.text.y = ggplot2::element_text(size = 6),
    legend.position = "right"
  )

if (!is.null(cluster_by_id)) {
  cluster_df <- data.frame(
    pos = seq_len(n),
    cluster = factor(cluster_by_id[ordered_ids], levels = names(cluster_colors))
  )
  p <- p +
    ggplot2::geom_point(
      data = cluster_df,
      ggplot2::aes(x = pos, y = n + 0.65, colour = cluster),
      inherit.aes = FALSE,
      shape = 15,
      size = 2.1
    ) +
    ggplot2::geom_point(
      data = cluster_df,
      ggplot2::aes(x = 0.35, y = n - pos + 1, colour = cluster),
      inherit.aes = FALSE,
      shape = 15,
      size = 2.1
    ) +
    ggplot2::scale_colour_manual(values = cluster_colors, name = "cluster") +
    ggplot2::coord_cartesian(xlim = c(0.25, n + 0.5), ylim = c(0.5, n + 0.85), clip = "off") +
    ggplot2::theme(plot.margin = ggplot2::margin(8, 14, 8, 14))
}

out_dir <- dirname(out_prefix)
if (!identical(out_dir, ".") && !dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

png_path <- paste0(out_prefix, ".png")
pdf_path <- paste0(out_prefix, ".pdf")
size_in <- max(7, min(18, 4.5 + n * 0.10))

ggplot2::ggsave(png_path, p, width = size_in, height = size_in, units = "in", dpi = 180, limitsize = FALSE)
ggplot2::ggsave(pdf_path, p, width = size_in, height = size_in, units = "in", limitsize = FALSE)

message("Wrote: ", png_path)
message("Wrote: ", pdf_path)
