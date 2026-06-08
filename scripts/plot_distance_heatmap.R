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
labels <- paste0("A", ids)
side_colors <- NULL
cluster_legend <- NULL
cluster_colors <- NULL

if (!is.null(alleles_path)) {
  alleles <- read.table(alleles_path, header = TRUE, sep = "\t",
                        check.names = FALSE, quote = "", comment.char = "")
  alleles$allele_id <- as.character(alleles$allele_id)
  alleles <- alleles[match(ids, alleles$allele_id), ]
  if (any(is.na(alleles$cluster_id))) {
    stop("Could not match all matrix allele IDs to alleles.tsv", call. = FALSE)
  }

  clusters <- sort(unique(alleles$cluster_id))
  palette <- c(
    "#4E79A7", "#F28E2B", "#59A14F", "#E15759", "#76B7B2", "#EDC948",
    "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC", "#2F4B7C", "#A05195"
  )
  cluster_colors <- setNames(palette[((seq_along(clusters) - 1) %% length(palette)) + 1],
                             as.character(clusters))
  side_colors <- cluster_colors[as.character(alleles$cluster_id)]
  cluster_legend <- clusters
  labels <- paste0("A", alleles$allele_id, " C", alleles$cluster_id, " n=", alleles$path_support)
}

if (n > 80) {
  labels <- rep("", n)
}

hc <- hclust(as.dist(dist_mat), method = "average")
heat_cols <- grDevices::colorRampPalette(
  c("#102A43", "#2F80ED", "#F7F7F2", "#F2994A", "#8C1D18")
)(256)

draw_heatmap <- function() {
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  par(bg = "white", cex.main = 1.1)

  heatmap_args <- list(
    x = dist_mat,
    Rowv = as.dendrogram(hc),
    Colv = as.dendrogram(hc),
    symm = TRUE,
    scale = "none",
    revC = TRUE,
    col = heat_cols,
    labRow = labels,
    labCol = labels,
    margins = c(8, 8),
    main = "Allele distance heatmap",
    xlab = "Alleles",
    ylab = "Alleles"
  )
  if (!is.null(side_colors)) {
    heatmap_args$RowSideColors <- side_colors
    heatmap_args$ColSideColors <- side_colors
  }
  do.call(heatmap, heatmap_args)

  if (!is.null(cluster_legend)) {
    legend("topright",
           legend = paste0("Cluster ", cluster_legend),
           fill = cluster_colors[as.character(cluster_legend)],
           border = NA, bty = "n", cex = 0.8)
  }
}

png_path <- paste0(out_prefix, ".png")
pdf_path <- paste0(out_prefix, ".pdf")
px <- max(1400, min(4200, 45 * n))

grDevices::png(png_path, width = px, height = px, res = 180)
draw_heatmap()
grDevices::dev.off()

grDevices::pdf(pdf_path, width = max(7, min(18, n * 0.16)), height = max(7, min(18, n * 0.16)))
draw_heatmap()
grDevices::dev.off()

message("Wrote: ", png_path)
message("Wrote: ", pdf_path)
