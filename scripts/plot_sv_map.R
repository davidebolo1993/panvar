#!/usr/bin/env Rscript

# Node-level structural-variant map (odgi-viz style): rows = haplotypes, columns = the
# bubble's internal nodes in REFERENCE-SORTED order (the same node order and look as
# plot_node_coverage_heatmap.R), each cell colored by the CALLED EVENT the haplotype
# carries on that node. Because every event is painted only on its own EVENT_NODES, an
# event occupies just its node columns (in genomic order) instead of spanning the row,
# and the figure lines up 1:1 with the node/edge-coverage heatmap of the same bubble.
#
# Cell colors:
#   white  = haplotype does not traverse the node (count 0, no event)
#   grey   = traverses the node, reference-like (no event here)
#   DEL red / INS-NOVEL green / INS-DUP purple / INV orange
#   DUP blue (shaded by per-sample CN) / multiallelic teal (shaded by allele index)
#
# Inputs (one bubble):
#   --node-counts <path>  inspect <prefix>.bubble_<N>.node_counts.tsv   (required; substrate + order)
#   --vcf <path>          panvar call VCF                               (required; events + carriers + CN)
#   --out <prefix>        writes <prefix>.png and <prefix>.pdf          (required)
#   --bubble-id <N>       keep only VCF records with INFO BUBBLE_ID=N (match the node_counts file)
#   --node-lengths <path> inspect <...>.node_lengths.tsv: scale column widths by node bp
#   --length-transform    raw | sqrt | log1p column-width scaling (default sqrt)
#   --clusters <path>     inspect <...>.clusters.tsv: order + annotate rows by walk cluster
#   --reference-path <s>  pin this haplotype (substring match) as the top row
#   --max-nodes <N>       keep at most N node columns (by total coverage)
#   --max-paths <N>       keep at most N haplotype rows (by total coverage)
#   --width / --height    figure size (inches; default auto)

args <- commandArgs(trailingOnly = TRUE)
usage <- function(status = 0) {
  cat(paste(c(
    "Usage:",
    "  plot_sv_map.R --node-counts <bubble_N.node_counts.tsv> --vcf <call.vcf> --out <prefix>",
    "                [--bubble-id N] [--node-lengths nl.tsv] [--length-transform sqrt]",
    "                [--clusters clusters.tsv] [--reference-path NAME] [--max-nodes N]",
    "                [--max-paths N] [--width in] [--height in]"), collapse = "\n"), "\n")
  quit(status = status)
}
if (length(args) == 0 || any(args %in% c("-h", "--help"))) usage(0)
opts <- list(node_counts = NULL, vcf = NULL, out = NULL, bubble_id = NULL,
             node_lengths = NULL, length_transform = "sqrt", clusters = NULL,
             reference_path = NULL, max_nodes = 0, max_paths = 0,
             width = NA_real_, height = NA_real_)
i <- 1
while (i <= length(args)) {
  a <- args[[i]]; val <- function() { if (i + 1 > length(args)) stop(paste("missing value after", a)); args[[i + 1]] }
  if (a == "--node-counts") { opts$node_counts <- val(); i <- i + 2 }
  else if (a == "--vcf") { opts$vcf <- val(); i <- i + 2 }
  else if (a == "--out") { opts$out <- val(); i <- i + 2 }
  else if (a == "--bubble-id") { opts$bubble_id <- as.integer(val()); i <- i + 2 }
  else if (a == "--node-lengths") { opts$node_lengths <- val(); i <- i + 2 }
  else if (a == "--length-transform") { opts$length_transform <- val(); i <- i + 2 }
  else if (a == "--clusters") { opts$clusters <- val(); i <- i + 2 }
  else if (a == "--reference-path") { opts$reference_path <- val(); i <- i + 2 }
  else if (a == "--max-nodes") { opts$max_nodes <- as.integer(val()); i <- i + 2 }
  else if (a == "--max-paths") { opts$max_paths <- as.integer(val()); i <- i + 2 }
  else if (a == "--width") { opts$width <- as.numeric(val()); i <- i + 2 }
  else if (a == "--height") { opts$height <- as.numeric(val()); i <- i + 2 }
  else if (startsWith(a, "-")) stop(paste("unknown option:", a)) else i <- i + 1
}
for (r in c("node_counts", "vcf", "out")) if (is.null(opts[[r]])) usage(1)
if (!opts$length_transform %in% c("raw", "sqrt", "log1p")) stop("--length-transform must be raw|sqrt|log1p")
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("needs ggplot2 (conda install -c conda-forge r-ggplot2)")

open_input <- function(p) if (grepl("\\.gz$", p, ignore.case = TRUE)) gzfile(p, open = "rt") else file(p, open = "rt")
read_tsv <- function(p) { con <- open_input(p); on.exit(close(con)); read.delim(con, sep = "\t", header = TRUE, check.names = FALSE, quote = "", comment.char = "") }
info_get <- function(info, key) {
  m <- regmatches(info, regexpr(paste0("(^|;)", key, "=[^;]*"), info))
  if (length(m) == 0) return(NA_character_)
  sub(paste0("^;?", key, "="), "", m)
}

# ---- substrate: node_counts.tsv -> matrix (paths x nodes), node order = column order ----
tab <- read_tsv(opts$node_counts)
if (!all(c("path_name", "path_length_bp") %in% names(tab))) stop("node_counts table needs path_name + path_length_bp columns")
node_cols <- grep("^node\\.", names(tab), value = TRUE)
if (length(node_cols) == 0) stop("node_counts table has no node.* columns")
if (nrow(tab) == 0) stop("node_counts table has no path rows")
cells <- as.character(unlist(tab[node_cols], use.names = FALSE))
counts <- suppressWarnings(as.numeric(sub(":.*$", "", cells)))   # total = first field of total:fwd:rev
counts[is.na(counts)] <- 0
mat <- matrix(counts, nrow = nrow(tab), ncol = length(node_cols), byrow = FALSE)
rownames(mat) <- make.unique(as.character(tab$path_name), sep = ".")
colnames(mat) <- sub("^node\\.", "", node_cols)

# ---- parse VCF: records of this bubble, with EVENT_NODES + per-sample GT/CN ----
vcon <- open_input(opts$vcf); vlines <- readLines(vcon); close(vcon)
hdr <- vlines[grepl("^#CHROM", vlines)][1]
vcf_samples <- strsplit(hdr, "\t")[[1]][-(1:9)]
recs <- vlines[!grepl("^#", vlines)]
# node_counts is per-bubble (one bubble's nodes); infer the bubble id from the filename if
# not given so records from other bubbles are not mis-painted onto this bubble's columns.
if (is.null(opts$bubble_id)) {
  m <- regmatches(basename(opts$node_counts), regexpr("bubble_([0-9]+)", basename(opts$node_counts)))
  if (length(m) == 1) opts$bubble_id <- as.integer(sub("bubble_", "", m))
}
V <- list()
for (ln in recs) {
  f <- strsplit(ln, "\t")[[1]]
  info <- f[8]
  bid <- info_get(info, "BUBBLE_ID")
  if (!is.null(opts$bubble_id) && (is.na(bid) || as.integer(bid) != opts$bubble_id)) next
  svt <- info_get(info, "SVTYPE")
  if (is.na(svt)) svt <- if (!is.na(info_get(info, "NALLELES"))) "MULTI" else "."
  enodes <- info_get(info, "EVENT_NODES")
  enodes <- if (is.na(enodes)) character() else strsplit(enodes, ",")[[1]]
  gt <- sub(":.*$", "", f[-(1:9)])
  cn <- suppressWarnings(as.integer(sub("^[^:]*:", "", f[-(1:9)])))
  V[[length(V) + 1]] <- list(id = f[3], svt = svt, sub = info_get(info, "INS_SUBTYPE"),
                             nodes = enodes, gt = setNames(gt, vcf_samples), cn = setNames(cn, vcf_samples))
}
if (length(V) == 0) stop("no VCF records for this bubble (check --bubble-id)")

# ---- restrict columns/rows (by coverage) before building the overlay ----
if (opts$max_nodes > 0 && ncol(mat) > opts$max_nodes) {
  keepc <- sort(order(colSums(mat), decreasing = TRUE)[seq_len(opts$max_nodes)])
  mat <- mat[, keepc, drop = FALSE]
}
if (opts$max_paths > 0 && nrow(mat) > opts$max_paths) {
  keepr <- sort(order(rowSums(mat), decreasing = TRUE)[seq_len(opts$max_paths)])
  mat <- mat[keepr, , drop = FALSE]
}
nodes <- colnames(mat)
paths <- rownames(mat)
node_idx <- setNames(seq_along(nodes), nodes)

# ---- category + value per cell: start from coverage (none/ref), overlay events on carriers ----
# A DUP is encoded with a single representative node (the peak / module source), so painting
# only EVENT_NODES would render it as one invisible column. Instead a DUP carrier's WHOLE
# traversed module is shaded by per-node multiplicity (the local copy count from the coverage
# substrate) -- the "CN palette over the coverage nodes". DEL/INS/INV/MULTI keep their specific
# EVENT_NODES and are painted AFTER, on top of the DUP background.
cat_mat <- ifelse(mat > 0, "ref", "none")            # paths x nodes
val_mat <- matrix(NA_real_, nrow(mat), ncol(mat), dimnames = dimnames(mat))  # multiplicity / allele index
dup_max <- 1; multi_max <- 1
carriers_of <- function(r) intersect(names(r$gt)[!is.na(r$gt) & r$gt != "." & r$gt != "0" & r$gt != "0/0"], paths)

# pass 1: DUP background (all traversed module nodes, shade = per-node multiplicity).
# Guard: only a DUP whose representative node belongs to this bubble paints here.
for (r in V) {
  if (r$svt != "DUP") next
  if (length(intersect(r$nodes, nodes)) == 0) next
  for (s in carriers_of(r)) {
    ri <- match(s, paths)
    for (ci in seq_len(ncol(mat))) {
      if (mat[ri, ci] > 0) { cat_mat[ri, ci] <- "DUP"; val_mat[ri, ci] <- mat[ri, ci]; dup_max <- max(dup_max, mat[ri, ci]) }
    }
  }
}
# pass 2: DEL/INS/INV/MULTI on their specific EVENT_NODES (override the DUP background)
for (r in V) {
  if (r$svt == "DUP") next
  cols_here <- node_idx[intersect(r$nodes, nodes)]
  if (length(cols_here) == 0) next
  cat_s <- if (r$svt == "DEL") "DEL"
    else if (r$svt == "INV") "INV"
    else if (r$svt == "INS") (if (!is.na(r$sub) && r$sub == "DUP") "INS_DUP" else "INS_NOVEL")
    else "MULTI"
  for (s in carriers_of(r)) {
    ri <- match(s, paths)
    for (ci in cols_here) {
      cat_mat[ri, ci] <- cat_s
      if (cat_s == "MULTI") { v <- suppressWarnings(as.integer(r$gt[[s]])); if (is.na(v) || v < 1) v <- 1; val_mat[ri, ci] <- v; multi_max <- max(multi_max, v) }
    }
  }
}

# ---- row ordering: reference pinned top, then clusters, then name ----
clust <- setNames(rep(NA_integer_, length(paths)), paths)
if (!is.null(opts$clusters)) {
  ct <- read_tsv(opts$clusters)
  if ("members" %in% names(ct) && "cluster_id" %in% names(ct)) {
    for (k in seq_len(nrow(ct))) {
      mem <- strsplit(as.character(ct$members[k]), ";")[[1]]
      clust[intersect(paths, mem)] <- as.integer(ct$cluster_id[k])
    }
  }
}
ref_row <- NA_character_
if (!is.null(opts$reference_path)) {
  hit <- paths[grepl(opts$reference_path, paths, ignore.case = TRUE, fixed = FALSE)]
  if (length(hit) >= 1) ref_row <- hit[1]
}
ord <- order(ifelse(is.na(clust), .Machine$integer.max, clust), paths)
paths <- paths[ord]
if (!is.na(ref_row)) paths <- c(ref_row, setdiff(paths, ref_row))
mat <- mat[paths, , drop = FALSE]; cat_mat <- cat_mat[paths, , drop = FALSE]; val_mat <- val_mat[paths, , drop = FALSE]

# ---- colors ----
COL <- c(none = "#ffffff", ref = "#e8e8e8", DEL = "#cb181d", INV = "#f16913",
         INS_NOVEL = "#238b45", INS_DUP = "#6a51a3")
dup_ramp <- grDevices::colorRampPalette(c("#c6dbef", "#08306b"))
multi_ramp <- grDevices::colorRampPalette(c("#a6e3d7", "#00665a"))
cell_color <- function(cat, v) {
  if (cat == "DUP") { vv <- if (is.na(v) || v < 1) 1 else v; return(dup_ramp(max(2, dup_max))[min(vv, dup_max)]) }
  if (cat == "MULTI") { vv <- if (is.na(v) || v < 1) 1 else v; return(multi_ramp(max(2, multi_max))[min(vv, multi_max)]) }
  COL[[cat]]
}

# ---- geometry: column widths optionally scaled by node bp (like the coverage heatmap) ----
n_nodes <- length(nodes); n_paths <- length(paths)
use_length <- !is.null(opts$node_lengths)
if (use_length) {
  nl <- read_tsv(opts$node_lengths)
  len_by_node <- setNames(as.numeric(nl$length_bp), as.character(nl$node_id))
  raw_len <- as.numeric(len_by_node[nodes])
  w <- switch(opts$length_transform, raw = raw_len, sqrt = sqrt(pmax(raw_len, 0)), log1p = log1p(pmax(raw_len, 0)))
  pos_min <- suppressWarnings(min(w[is.finite(w) & w > 0])); if (!is.finite(pos_min)) pos_min <- 1
  w[!is.finite(w) | w <= 0] <- pos_min * 0.25
  x1 <- cumsum(w); x0 <- x1 - w; xmid <- (x0 + x1) / 2
} else { x0 <- seq_len(n_nodes) - 1; x1 <- seq_len(n_nodes); xmid <- seq_len(n_nodes) - 0.5 }

df <- expand.grid(ci = seq_len(n_nodes), ri = seq_len(n_paths))
df$fill <- mapply(function(ci, ri) cell_color(cat_mat[ri, ci], val_mat[ri, ci]), df$ci, df$ri)
df$path_y <- n_paths - df$ri + 1
df$xmin <- x0[df$ci]; df$xmax <- x1[df$ci]; df$ymin <- df$path_y - 0.5; df$ymax <- df$path_y + 0.5

# ---- render ----
sparse <- function(n, k) if (n <= 0) integer() else unique(pmax(1, pmin(n, round(seq(1, n, length.out = min(n, k))))))
xt <- sparse(n_nodes, 40); yt <- sparse(n_paths, 50)
x_breaks <- if (use_length) xmid[xt] else xt
out_dir <- dirname(opts$out); if (!identical(out_dir, ".") && !dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (is.na(opts$width)) opts$width <- max(7, min(30, 3 + n_nodes / 12))
if (is.na(opts$height)) opts$height <- max(4, min(40, 1.5 + n_paths / 9))

p <- ggplot2::ggplot(df) +
  ggplot2::geom_rect(ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill)) +
  ggplot2::scale_fill_identity() +
  ggplot2::scale_x_continuous(breaks = x_breaks, labels = nodes[xt], expand = c(0, 0)) +
  ggplot2::scale_y_continuous(breaks = n_paths - yt + 1, labels = paths[yt], expand = c(0, 0)) +
  ggplot2::labs(title = "Structural-variant map (haplotypes x bubble nodes)",
                subtitle = sprintf("%d haplotypes x %d nodes%s%s", n_paths, n_nodes,
                                   if (use_length) " (x scaled by node bp)" else "",
                                   if (!is.na(ref_row)) sprintf("; reference '%s' on top", ref_row) else ""),
                x = if (use_length) "Bubble nodes (reference-sorted; width = bp)" else "Bubble nodes (reference-sorted)",
                y = "Haplotypes") +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(panel.grid = ggplot2::element_blank(),
                 plot.title = ggplot2::element_text(face = "bold"),
                 axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5),
                 axis.text.y = ggplot2::element_text(size = 4))

# optional cluster color bar on the left
if (!is.null(opts$clusters)) {
  bar_x <- if (use_length) -max(x1) * 0.02 else -1
  cb <- data.frame(path_y = n_paths - seq_len(n_paths) + 1, cl = factor(clust[paths]))
  p <- p + ggplot2::geom_tile(data = cb, ggplot2::aes(x = bar_x, y = path_y, color = cl), fill = NA,
                              width = if (use_length) max(x1) * 0.015 else 0.8, height = 0.9, na.rm = TRUE) +
    ggplot2::scale_color_discrete(name = "cluster", na.value = "#cccccc") +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(shape = 15, size = 4)))
}

# manual SV color legend
legend_levels <- c("not traversed" = COL[["none"]], "reference-like" = COL[["ref"]],
                   "DEL" = COL[["DEL"]], "INV" = COL[["INV"]], "INS/NOVEL" = COL[["INS_NOVEL"]],
                   "INS/DUP" = COL[["INS_DUP"]], "DUP (shade=copies)" = dup_ramp(3)[3],
                   "multiallelic" = multi_ramp(3)[3])
leg <- data.frame(lab = factor(names(legend_levels), levels = names(legend_levels)), x = NA_real_, y = NA_real_)
p <- p + ggplot2::geom_point(data = leg, ggplot2::aes(x = x, y = y, shape = lab), na.rm = TRUE) +
  ggplot2::scale_shape_manual(name = "SV", values = rep(15, length(legend_levels))) +
  ggplot2::guides(shape = ggplot2::guide_legend(override.aes = list(color = unname(legend_levels), size = 5)))

png_path <- paste0(opts$out, ".png"); pdf_path <- paste0(opts$out, ".pdf")
ggplot2::ggsave(png_path, p, width = opts$width, height = opts$height, units = "in", dpi = 180, limitsize = FALSE)
ggplot2::ggsave(pdf_path, p, width = opts$width, height = opts$height, units = "in", limitsize = FALSE)
cat("Wrote:", png_path, "\n"); cat("Wrote:", pdf_path, "\n")
