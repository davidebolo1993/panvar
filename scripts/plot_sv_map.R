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
#   DUP blue, shaded by the haplotype's ABSOLUTE copy number (VCF FORMAT:CN) -- painted on the REP node
#       for a self-loop DUP, or across the whole traversed module for a coverage DUP; losses included
#   multiallelic teal (allele index)
#
# Inputs (one bubble):
#   --node-counts <path>  inspect <prefix>.bubble_<N>.node_counts.tsv   (required; substrate + order)
#   --vcf <path>          panvar call VCF                               (required; events + carriers + CN)
#   --out <prefix>        writes <prefix>.png and <prefix>.pdf          (required)
#   --bubble-id <N>       keep only VCF records with INFO BUBBLE_ID=N (match the node_counts file)
#   --node-lengths <path> inspect <...>.node_lengths.tsv: scale column widths by node bp
#   --length-transform    raw | sqrt | log1p column-width scaling (default sqrt)
#   --clusters <path>     inspect <...>.clusters.tsv: keep only cluster representative rows
#   --cluster-by <path>   inspect <...>.clusters.tsv: keep all rows but group/order by cluster
#   --reference-path <s>  pin this haplotype (substring match) as the top row
#   --max-nodes <N>       keep at most N node columns (by total coverage)
#   --max-paths <N>       keep at most N haplotype rows (by total coverage)
#   --width / --height    figure size (inches; default auto)

args <- commandArgs(trailingOnly = TRUE)
usage <- function(status = 0) {
  cat(paste(c(
    "plot_sv_map.R - node-level structural-variant map (odgi-viz style).",
    "  rows = haplotypes, columns = a bubble's internal nodes in reference-sorted order;",
    "  each cell is colored by the called event the haplotype carries on that node, so the",
    "  figure lines up 1:1 with the inspect node/edge-coverage heatmaps of the same bubble.",
    "",
    "Usage:",
    "  Rscript plot_sv_map.R --node-counts <bubble_N.node_counts.tsv> --vcf <call.vcf> --out <prefix> [options]",
    "",
    "Required:",
    "  --node-counts <path>     inspect <prefix>.bubble_<N>.node_counts.tsv (substrate + node order)",
    "  --vcf <path>             panvar call VCF (events, carriers, GT:CN)",
    "  --out <prefix>           output prefix; writes <prefix>.png and <prefix>.pdf",
    "",
    "Optional:",
    "  --bubble-id <N>          keep only VCF records with INFO BUBBLE_ID=N (else inferred from filename)",
    "  --node-lengths <path>    inspect <...>.node_lengths.tsv: scale column widths by node bp",
    "  --length-transform <t>   raw | sqrt | log1p column-width scaling (default sqrt)",
    "  --clusters <path>        inspect <...>.clusters.tsv: keep only the cluster representative rows",
    "  --cluster-by <path>      inspect <...>.clusters.tsv: keep all rows but group/order them by",
    "                           cluster (representative first), with a separator between clusters",
    "  --reference-path <s>     pin this haplotype (substring match) as the top row",
    "  --max-nodes <N>          keep at most N node columns (by total coverage)",
    "  --max-paths <N>          keep at most N haplotype rows (by total coverage)",
    "  --width / --height <in>  figure size in inches (default auto)",
    "  -h, --help               show this help"), collapse = "\n"), "\n")
  quit(status = status)
}
if (length(args) == 0 || any(args %in% c("-h", "--help"))) usage(0)
opts <- list(node_counts = NULL, vcf = NULL, out = NULL, bubble_id = NULL,
             node_lengths = NULL, length_transform = "sqrt", clusters = NULL,
             cluster_by = NULL, reference_path = NULL, max_nodes = 0, max_paths = 0,
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
  else if (a == "--cluster-by") { opts$cluster_by <- val(); i <- i + 2 }
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
  ref_cn <- suppressWarnings(as.integer(info_get(info, "REF_CN")))
  V[[length(V) + 1]] <- list(id = f[3], svt = svt, sub = info_get(info, "INS_SUBTYPE"),
                             nodes = enodes, ref_cn = ref_cn,
                             gt = setNames(gt, vcf_samples), cn = setNames(cn, vcf_samples))
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
# --clusters: keep only cluster representative rows (like the coverage heatmaps)
if (!is.null(opts$clusters)) {
  ct0 <- read_tsv(opts$clusters)
  if (!"representative_path" %in% names(ct0)) stop("--clusters table needs a representative_path column")
  reps <- unique(as.character(ct0$representative_path))
  keep_rep <- rownames(mat) %in% reps
  if (!any(keep_rep)) stop("no path rows match the cluster representatives")
  mat <- mat[keep_rep, , drop = FALSE]
}
nodes <- colnames(mat)
paths <- rownames(mat)
node_idx <- setNames(seq_along(nodes), nodes)

# ---- category + value per cell: start from coverage (none/ref), overlay the called CN / events ----
# A DUP is shaded by the haplotype's ABSOLUTE copy number (the VCF FORMAT:CN), for EVERY haplotype that
# traverses the locus -- not just gain/loss carriers -- so losses (CN < REF_CN), reference (CN = REF_CN)
# and gains (CN > REF_CN) all read as a sequential shade. Two cases, by how the DUP is encoded:
#   * REP / self-loop DUP: EVENT_NODES is a real inside node (the repeat unit), so the CN is painted on
#     that column -- one clean copy-number strip (e.g. LPA KIV-2 at the REP node).
#   * whole-module coverage DUP: the event node is the bubble boundary (not an inside column), so the CN
#     is painted across every inside node the haplotype traverses -- the whole module row takes one CN
#     shade, so C4/CYP2D6 copy number reads as horizontal bands.
# DEL/INS/INV/MULTI keep their specific EVENT_NODES and are painted AFTER, on top of the DUP background.
cat_mat <- ifelse(mat > 0, "ref", "none")            # paths x nodes
val_mat <- matrix(NA_real_, nrow(mat), ncol(mat), dimnames = dimnames(mat))  # called CN / allele index
dup_max <- 1; multi_max <- 1
ref_cn_seen <- NA_integer_
# carriers (gain/loss, GT=1) -- used only for the DEL/INS/INV layer below
carriers_of <- function(r) intersect(names(r$gt)[!is.na(r$gt) & r$gt != "." & r$gt != "0" & r$gt != "0/0"], paths)
# every haplotype that has a CN at this record (carriers explicit, reference-like = REF_CN)
cn_haps_of <- function(r) intersect(names(r$cn)[!is.na(r$cn)], paths)

# pass 1: DUP background -- shade by absolute CN (VCF), every CN-bearing haplotype.
for (r in V) {
  if (r$svt != "DUP") next
  cols_event <- node_idx[intersect(r$nodes, nodes)]   # REP / self-loop node columns, if any
  for (s in cn_haps_of(r)) {
    ri <- match(s, paths)
    cn <- r$cn[[s]]; if (is.na(cn)) next
    if (length(cols_event) > 0) {
      cols_here <- cols_event                          # REP DUP: the repeat-unit column(s)
    } else {
      cols_here <- which(mat[ri, ] > 0)                # coverage DUP: the whole traversed module
    }
    for (ci in cols_here) { cat_mat[ri, ci] <- "DUP"; val_mat[ri, ci] <- cn }
    if (cn > dup_max) dup_max <- cn
  }
  if (is.na(ref_cn_seen) && !is.na(r$ref_cn)) ref_cn_seen <- r$ref_cn
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

# ---- row ordering: optional cluster grouping (--cluster-by), then reference pinned on top ----
# Mirrors the inspect coverage heatmaps: --cluster-by groups cluster-mates adjacently
# (representative first), with a thin separator line between clusters and NO color sidebar.
cid_of <- setNames(rep(NA_integer_, length(paths)), paths)
cluster_boundaries <- integer(0)
if (!is.null(opts$cluster_by)) {
  ct <- read_tsv(opts$cluster_by)
  if (!all(c("cluster_id", "members") %in% names(ct))) stop("--cluster-by table needs cluster_id + members columns")
  has_rep <- "representative_path" %in% names(ct)
  isrep_of <- setNames(rep(FALSE, length(paths)), paths)
  for (k in seq_len(nrow(ct))) {
    mem <- strsplit(as.character(ct$members[k]), ";", fixed = TRUE)[[1]]; mem <- mem[nzchar(mem)]
    cid_of[intersect(paths, mem)] <- as.integer(ct$cluster_id[k])
    if (has_rep) isrep_of[intersect(paths, as.character(ct$representative_path[k]))] <- TRUE
  }
  cid_sort <- ifelse(is.na(cid_of[paths]), .Machine$integer.max, cid_of[paths])
  paths <- paths[order(cid_sort, !isrep_of[paths], paths)]
} else {
  paths <- paths[order(paths)]
}
ref_row <- NA_character_
if (!is.null(opts$reference_path)) {
  hit <- paths[grepl(opts$reference_path, paths, ignore.case = TRUE, fixed = FALSE)]
  if (length(hit) >= 1) ref_row <- hit[1]
}
if (!is.na(ref_row)) paths <- c(ref_row, setdiff(paths, ref_row))
# separator rows where the cluster id changes, computed on the FINAL order (--cluster-by only)
if (!is.null(opts$cluster_by)) {
  fc <- ifelse(is.na(cid_of[paths]), .Machine$integer.max, cid_of[paths])
  if (length(fc) > 1) cluster_boundaries <- which(fc[-1] != fc[-length(fc)]) + 1
}
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
                subtitle = sprintf("%d haplotypes x %d nodes%s%s%s", n_paths, n_nodes,
                                   if (use_length) " (x scaled by node bp)" else "",
                                   if (!is.na(ref_cn_seen)) sprintf("; DUP shade = absolute CN (REF_CN=%d)", ref_cn_seen) else "",
                                   if (!is.na(ref_row)) sprintf("; reference '%s' on top", ref_row) else ""),
                x = if (use_length) "Bubble nodes (reference-sorted; width = bp)" else "Bubble nodes (reference-sorted)",
                y = "Haplotypes") +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(panel.grid = ggplot2::element_blank(),
                 plot.title = ggplot2::element_text(face = "bold"),
                 axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5),
                 axis.text.y = ggplot2::element_text(size = 4))

# separator lines between adjacent inspect clusters (only with --cluster-by); no color sidebar
if (length(cluster_boundaries) > 0) {
  p <- p + ggplot2::geom_hline(yintercept = n_paths - cluster_boundaries + 1.5,
                               colour = "grey30", linewidth = 0.3)
}

# manual SV color legend
legend_levels <- c("not traversed" = COL[["none"]], "reference-like" = COL[["ref"]],
                   "DEL" = COL[["DEL"]], "INV" = COL[["INV"]], "INS/NOVEL" = COL[["INS_NOVEL"]],
                   "INS/DUP" = COL[["INS_DUP"]], "DUP (shade=absolute CN)" = dup_ramp(3)[3],
                   "multiallelic" = multi_ramp(3)[3])
leg <- data.frame(lab = factor(names(legend_levels), levels = names(legend_levels)), x = NA_real_, y = NA_real_)
p <- p + ggplot2::geom_point(data = leg, ggplot2::aes(x = x, y = y, shape = lab), na.rm = TRUE) +
  ggplot2::scale_shape_manual(name = "SV", values = rep(15, length(legend_levels))) +
  ggplot2::guides(shape = ggplot2::guide_legend(override.aes = list(color = unname(legend_levels), size = 5)))

png_path <- paste0(opts$out, ".png"); pdf_path <- paste0(opts$out, ".pdf")
ggplot2::ggsave(png_path, p, width = opts$width, height = opts$height, units = "in", dpi = 180, limitsize = FALSE)
ggplot2::ggsave(pdf_path, p, width = opts$width, height = opts$height, units = "in", limitsize = FALSE)
cat("Wrote:", png_path, "\n"); cat("Wrote:", pdf_path, "\n")
