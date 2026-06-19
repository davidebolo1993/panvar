#!/usr/bin/env Rscript

# Whole-VCF variant map (oncoprint / genotype-matrix style): ONE figure for an entire panvar
# call VCF. Rows = haplotypes, columns = the called variants, GROUPED BY BUBBLE (one facet per
# bubble), each column labelled by its variant ID. Every cell is the called event the haplotype
# carries at that variant, so the figure reads like the VCF itself -- not like the graph. (For the
# node-order / graph-structure view of a single bubble, use plot_sv_map.R + the inspect coverage
# heatmaps; this map is the at-a-glance whole-VCF summary.)
#
# Cell colors:
#   grey   = reference-like (haplotype does not carry this variant)
#   DEL red / INS-NOVEL green / INS-DUP purple / INV orange / multiallelic teal (allele index)
#   DUP blue, shaded by the haplotype's ABSOLUTE copy number (VCF FORMAT:CN) for EVERY haplotype
#       (loss = light, reference = mid, gain = dark), so the DUP column reads as a CN gradient
#
# Inputs:
#   --vcf <path>          panvar call VCF                                  (required; events + GT:CN)
#   --out <prefix>        writes <prefix>.png and <prefix>.pdf            (required)
#   --clusters <path>     inspect <...>.clusters.tsv: keep only cluster representative rows
#   --cluster-by <path>   inspect <...>.clusters.tsv: keep all rows but group/order them by cluster
#   --reference-path <s>  pin this haplotype (substring match) as the top (or leftmost, if --flip) row
#   --max-paths <N>       keep at most N haplotype rows (those carrying the most events)
#   --title <str>         plot title (optional)
#   --flip                transpose: variants on Y, haplotypes on X (legend moves to the bottom)
#   --scale               scale each variant's rectangle along the variant axis by its size
#                         (|SVLEN| for DEL/INS/INV/multiallelic; RU_LEN repeat-unit bp for DUP)
#   --scale-transform <t> raw | sqrt | log1p size->extent scaling for --scale (default sqrt)
#   --width / --height    figure size (inches; default auto)

args <- commandArgs(trailingOnly = TRUE)
usage <- function(status = 0) {
  cat(paste(c(
    "plot_vcf_map.R - whole-VCF variant map (oncoprint style).",
    "  rows = haplotypes, columns = called variants grouped by bubble (one facet per bubble),",
    "  each cell colored by the called event the haplotype carries; one figure for the whole VCF.",
    "",
    "Usage:",
    "  Rscript plot_vcf_map.R --vcf <call.vcf> --out <prefix> [options]",
    "",
    "Required:",
    "  --vcf <path>             panvar call VCF (events, carriers, GT:CN)",
    "  --out <prefix>           output prefix; writes <prefix>.png and <prefix>.pdf",
    "",
    "Optional:",
    "  --clusters <path>        inspect <...>.clusters.tsv: keep only the cluster representative rows",
    "  --cluster-by <path>      inspect <...>.clusters.tsv: keep all rows but group/order them by cluster",
    "  --reference-path <s>     pin this haplotype (substring match) as the top/leftmost row",
    "  --max-paths <N>          keep at most N haplotype rows (by number of events carried)",
    "  --title <str>            plot title",
    "  --flip                   transpose: variants on Y, haplotypes on X (legend at bottom)",
    "  --scale                  scale each variant rectangle by its size (SVLEN; RU_LEN for DUP)",
    "  --scale-transform <t>    raw | sqrt | log1p scaling for --scale (default sqrt)",
    "  --width / --height <in>  figure size in inches (default auto)",
    "  -h, --help               show this help"), collapse = "\n"), "\n")
  quit(status = status)
}
if (length(args) == 0 || any(args %in% c("-h", "--help"))) usage(0)
opts <- list(vcf = NULL, out = NULL, clusters = NULL, cluster_by = NULL,
             reference_path = NULL, max_paths = 0, title = NULL,
             flip = FALSE, scale = FALSE, scale_transform = "sqrt",
             width = NA_real_, height = NA_real_)
i <- 1
while (i <= length(args)) {
  a <- args[[i]]; val <- function() { if (i + 1 > length(args)) stop(paste("missing value after", a)); args[[i + 1]] }
  if (a == "--vcf") { opts$vcf <- val(); i <- i + 2 }
  else if (a == "--out") { opts$out <- val(); i <- i + 2 }
  else if (a == "--clusters") { opts$clusters <- val(); i <- i + 2 }
  else if (a == "--cluster-by") { opts$cluster_by <- val(); i <- i + 2 }
  else if (a == "--reference-path") { opts$reference_path <- val(); i <- i + 2 }
  else if (a == "--max-paths") { opts$max_paths <- as.integer(val()); i <- i + 2 }
  else if (a == "--title") { opts$title <- val(); i <- i + 2 }
  else if (a == "--flip") { opts$flip <- TRUE; i <- i + 1 }
  else if (a == "--scale") { opts$scale <- TRUE; i <- i + 1 }
  else if (a == "--scale-transform") { opts$scale_transform <- val(); i <- i + 2 }
  else if (a == "--width") { opts$width <- as.numeric(val()); i <- i + 2 }
  else if (a == "--height") { opts$height <- as.numeric(val()); i <- i + 2 }
  else if (startsWith(a, "-")) stop(paste("unknown option:", a)) else i <- i + 1
}
for (r in c("vcf", "out")) if (is.null(opts[[r]])) usage(1)
if (!opts$scale_transform %in% c("raw", "sqrt", "log1p")) stop("--scale-transform must be raw|sqrt|log1p")
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("needs ggplot2 (conda install -c conda-forge r-ggplot2)")

open_input <- function(p) if (grepl("\\.gz$", p, ignore.case = TRUE)) gzfile(p, open = "rt") else file(p, open = "rt")
read_tsv <- function(p) { con <- open_input(p); on.exit(close(con)); read.delim(con, sep = "\t", header = TRUE, check.names = FALSE, quote = "", comment.char = "") }
info_get <- function(info, key) {
  m <- regmatches(info, regexpr(paste0("(^|;)", key, "=[^;]*"), info))
  if (length(m) == 0) return(NA_character_)
  sub(paste0("^;?", key, "="), "", m)
}

# ---- parse VCF: every record with its bubble, type, and per-sample GT/CN ----
vcon <- open_input(opts$vcf); vlines <- readLines(vcon); close(vcon)
hdr <- vlines[grepl("^#CHROM", vlines)][1]
if (is.na(hdr)) stop("no #CHROM header in VCF")
samples <- strsplit(hdr, "\t")[[1]][-(1:9)]
recs <- vlines[!grepl("^#", vlines) & nzchar(vlines)]
if (length(recs) == 0) stop("no variant records in VCF")
V <- list()
for (ln in recs) {
  f <- strsplit(ln, "\t")[[1]]
  info <- f[8]
  svt <- info_get(info, "SVTYPE")
  if (is.na(svt)) svt <- if (!is.na(info_get(info, "NALLELES"))) "MULTI" else "."
  bid <- suppressWarnings(as.integer(info_get(info, "BUBBLE_ID")))
  gt <- sub(":.*$", "", f[-(1:9)])
  cn <- suppressWarnings(as.integer(sub("^[^:]*:", "", f[-(1:9)])))
  # Variant size for --scale: the representative ALT-allele size = |SVLEN|; for a DUP the REPEAT-UNIT
  # length (RU_LEN, one copy) instead, so a high-copy DUP is not drawn huge just because of its CN.
  svlen <- suppressWarnings(abs(as.numeric(info_get(info, "SVLEN"))))
  ru <- suppressWarnings(as.numeric(info_get(info, "RU_LEN")))
  size_bp <- if (svt == "DUP" && !is.na(ru)) ru else svlen
  V[[length(V) + 1]] <- list(id = f[3], svt = svt, sub = info_get(info, "INS_SUBTYPE"),
                             bubble = if (is.na(bid)) 0L else bid,
                             ref_cn = suppressWarnings(as.integer(info_get(info, "REF_CN"))),
                             size_bp = if (is.na(size_bp)) NA_real_ else size_bp,
                             gt = setNames(gt, samples), cn = setNames(cn, samples))
}
# columns: order by bubble, then by appearance in the VCF
ord <- order(vapply(V, function(r) r$bubble, integer(1)), seq_along(V))
V <- V[ord]
var_ids <- vapply(V, function(r) r$id, character(1))
var_bub <- vapply(V, function(r) r$bubble, integer(1))
var_size <- vapply(V, function(r) r$size_bp, numeric(1))
n_var <- length(V)

# ---- restrict rows: keep all samples, optionally collapse to cluster reps / cap row count ----
paths <- samples
if (!is.null(opts$clusters)) {
  ct0 <- read_tsv(opts$clusters)
  if (!"representative_path" %in% names(ct0)) stop("--clusters table needs a representative_path column")
  reps <- unique(as.character(ct0$representative_path))
  keep <- paths[paths %in% reps]
  if (length(keep) == 0) stop("no samples match the cluster representatives")
  paths <- keep
}

# ---- per-(path,variant) category + value ----
# DUP: every haplotype is shaded by its absolute CN (loss/ref/gain), not just gain/loss carriers.
# DEL/INS/INV/MULTI: carriers only (GT calls an alt allele); everyone else is reference-like.
is_carrier <- function(g) !is.na(g) & g != "." & g != "0" & g != "0/0" & g != "0|0"
cat_mat <- matrix("ref", nrow = length(paths), ncol = n_var, dimnames = list(paths, var_ids))
val_mat <- matrix(NA_real_, nrow = length(paths), ncol = n_var, dimnames = list(paths, var_ids))
dup_max <- 1; multi_max <- 1; ref_cn_seen <- NA_integer_
for (j in seq_len(n_var)) {
  r <- V[[j]]
  if (r$svt == "DUP") {
    cn <- r$cn[paths]
    have <- which(!is.na(cn))
    cat_mat[have, j] <- "DUP"; val_mat[have, j] <- cn[have]
    if (length(have)) dup_max <- max(dup_max, max(cn[have]))
    if (is.na(ref_cn_seen) && !is.na(r$ref_cn)) ref_cn_seen <- r$ref_cn
    next
  }
  g <- r$gt[paths]
  carr <- which(is_carrier(g))
  if (length(carr) == 0) next
  cat_s <- if (r$svt == "DEL") "DEL"
    else if (r$svt == "INV") "INV"
    else if (r$svt == "INS") (if (!is.na(r$sub) && r$sub == "DUP") "INS_DUP" else "INS_NOVEL")
    else "MULTI"
  cat_mat[carr, j] <- cat_s
  if (cat_s == "MULTI") {
    v <- suppressWarnings(as.integer(g[carr])); v[is.na(v) | v < 1] <- 1
    val_mat[carr, j] <- v; multi_max <- max(multi_max, max(v))
  }
}

# drop haplotypes with no event at all (all-reference rows add height but no signal), unless that
# would empty the plot; keep the reference even if event-free so it can be pinned.
ref_row <- NA_character_
if (!is.null(opts$reference_path)) {
  hit <- paths[grepl(opts$reference_path, paths, ignore.case = TRUE)]
  if (length(hit) >= 1) ref_row <- hit[1]
}
has_event <- apply(cat_mat, 1, function(z) any(z != "ref"))
keep <- has_event | (rownames(cat_mat) == ref_row)
if (any(keep)) { paths <- paths[keep]; cat_mat <- cat_mat[keep, , drop = FALSE]; val_mat <- val_mat[keep, , drop = FALSE] }

# ---- row ordering: cluster grouping (optional), then oncoprint memo-sort, reference pinned ----
cid_of <- setNames(rep(NA_integer_, length(paths)), paths)
cluster_boundaries <- integer(0)
isrep_of <- setNames(rep(FALSE, length(paths)), paths)
if (!is.null(opts$cluster_by)) {
  ct <- read_tsv(opts$cluster_by)
  if (!all(c("cluster_id", "members") %in% names(ct))) stop("--cluster-by table needs cluster_id + members columns")
  has_rep <- "representative_path" %in% names(ct)
  for (k in seq_len(nrow(ct))) {
    mem <- strsplit(as.character(ct$members[k]), ";", fixed = TRUE)[[1]]; mem <- mem[nzchar(mem)]
    cid_of[intersect(paths, mem)] <- as.integer(ct$cluster_id[k])
    if (has_rep) isrep_of[intersect(paths, as.character(ct$representative_path[k]))] <- TRUE
  }
  cid_sort <- ifelse(is.na(cid_of[paths]), .Machine$integer.max, cid_of[paths])
  paths <- paths[order(cid_sort, !isrep_of[paths], paths)]
} else {
  # oncoprint memo-sort: group haplotypes with the same event pattern adjacently. Build a per-row
  # string key from the column categories (left-to-right, so high-frequency columns dominate).
  rank <- c(none = 0, ref = 0, DEL = 1, INV = 2, INS_NOVEL = 3, INS_DUP = 4, MULTI = 5, DUP = 6)
  keymat <- matrix(sprintf("%d", rank[cat_mat]), nrow = nrow(cat_mat))
  key <- apply(keymat, 1, paste, collapse = "")
  paths <- paths[order(key, decreasing = TRUE)]
}
if (!is.na(ref_row)) paths <- c(ref_row, setdiff(paths, ref_row))
if (opts$max_paths > 0 && length(paths) > opts$max_paths) {
  pin <- if (!is.na(ref_row)) ref_row else character()
  paths <- c(pin, setdiff(paths, pin)[seq_len(opts$max_paths - length(pin))])
}
cat_mat <- cat_mat[paths, , drop = FALSE]; val_mat <- val_mat[paths, , drop = FALSE]
if (!is.null(opts$cluster_by)) {
  fc <- ifelse(is.na(cid_of[paths]), .Machine$integer.max, cid_of[paths])
  if (length(fc) > 1) cluster_boundaries <- which(fc[-1] != fc[-length(fc)]) + 1
}

# ---- colors ----
COL <- c(ref = "#e8e8e8", DEL = "#cb181d", INV = "#f16913",
         INS_NOVEL = "#238b45", INS_DUP = "#6a51a3")
dup_ramp <- grDevices::colorRampPalette(c("#c6dbef", "#08306b"))
multi_ramp <- grDevices::colorRampPalette(c("#a6e3d7", "#00665a"))
cell_color <- function(cat, v) {
  if (cat == "DUP") { vv <- if (is.na(v) || v < 1) 1 else v; return(dup_ramp(max(2, dup_max))[min(vv, dup_max)]) }
  if (cat == "MULTI") { vv <- if (is.na(v) || v < 1) 1 else v; return(multi_ramp(max(2, multi_max))[min(vv, multi_max)]) }
  COL[[cat]]
}

# ---- variant-axis layout: each variant gets an extent along its axis (equal, or scaled by size) ----
# --scale draws each variant's rectangle proportional to its size (|SVLEN| for DEL/INS/INV, RU_LEN for
# DUP). Positions are GLOBAL cumulative (variants are bubble-ordered), so under facet free-scales each
# bubble panel clips cleanly to its own contiguous range and space="free" sizes panels by content bp.
n_paths <- length(paths)
if (opts$scale) {
  sz <- var_size; sz[!is.finite(sz) | sz <= 0] <- NA
  w <- switch(opts$scale_transform, raw = sz, sqrt = sqrt(sz), log1p = log1p(sz))
  posmin <- suppressWarnings(min(w[is.finite(w) & w > 0])); if (!is.finite(posmin)) posmin <- 1
  w[!is.finite(w) | w <= 0] <- posmin * 0.3   # unknown / zero-size variant -> a thin slab
  w <- w / mean(w)                            # normalize so mean width ~ 1 (stable figure proportions)
} else {
  w <- rep(1, n_var)
}
v1 <- cumsum(w); v0 <- v1 - w; vmid <- (v0 + v1) / 2

# ---- build one rectangle per (path, variant) ----
G <- expand.grid(j = seq_len(n_var), ri = seq_len(n_paths))
G$fill <- mapply(function(j, ri) cell_color(cat_mat[ri, j], val_mat[ri, j]), G$j, G$ri)
G$bubble <- factor(sprintf("bubble %d", var_bub[G$j]), levels = sprintf("bubble %d", sort(unique(var_bub))))
G$vmin <- v0[G$j]; G$vmax <- v1[G$j]                 # variant-axis extent
# hap axis is categorical (equal width 1); orientation decides which axis is which.
if (!opts$flip) {
  hp <- n_paths - G$ri + 1                           # haplotypes on Y, first path on top
  G$xmin <- G$vmin; G$xmax <- G$vmax; G$ymin <- hp - 0.5; G$ymax <- hp + 0.5
} else {
  hp <- G$ri                                         # haplotypes on X, first path on left
  G$ymin <- G$vmin; G$ymax <- G$vmax; G$xmin <- hp - 0.5; G$xmax <- hp + 0.5
}

# ---- render ----
sparse <- function(n, k) if (n <= 0) integer() else unique(pmax(1, pmin(n, round(seq(1, n, length.out = min(n, k))))))
ht <- sparse(n_paths, 50)                            # sparse haplotype ticks
hap_pos <- if (!opts$flip) n_paths - ht + 1 else ht
out_dir <- dirname(opts$out); if (!identical(out_dir, ".") && !dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
var_axis_lab <- sprintf("Called variant (grouped by bubble%s)", if (opts$scale) sprintf("; %s width = bp", opts$scale_transform) else "")
title <- if (!is.null(opts$title)) opts$title else "Variant map (haplotypes x called variants)"

p <- ggplot2::ggplot(G) +
  ggplot2::geom_rect(ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
                     colour = "white", linewidth = 0.1) +
  ggplot2::scale_fill_identity() +
  ggplot2::labs(title = title) +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(panel.grid = ggplot2::element_blank(),
                 plot.title = ggplot2::element_text(face = "bold"),
                 strip.text = ggplot2::element_text(size = 7, face = "bold"))

if (!opts$flip) {
  p <- p +
    ggplot2::facet_grid(. ~ bubble, scales = "free_x", space = "free_x") +
    ggplot2::scale_x_continuous(breaks = vmid, labels = var_ids, expand = c(0, 0)) +
    ggplot2::scale_y_continuous(breaks = hap_pos, labels = paths[ht], expand = c(0, 0)) +
    ggplot2::labs(x = var_axis_lab, y = "Haplotypes") +
    ggplot2::theme(panel.spacing.x = ggplot2::unit(2, "pt"),
                   axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
                   axis.text.y = ggplot2::element_text(size = 4))
  if (length(cluster_boundaries) > 0)
    p <- p + ggplot2::geom_hline(yintercept = n_paths - cluster_boundaries + 1.5, colour = "grey30", linewidth = 0.3)
} else {
  p <- p +
    ggplot2::facet_grid(bubble ~ ., scales = "free_y", space = "free_y") +
    ggplot2::scale_y_continuous(breaks = vmid, labels = var_ids, expand = c(0, 0)) +
    ggplot2::scale_x_continuous(breaks = hap_pos, labels = paths[ht], expand = c(0, 0)) +
    ggplot2::labs(x = "Haplotypes", y = var_axis_lab) +
    ggplot2::theme(panel.spacing.y = ggplot2::unit(2, "pt"),
                   strip.text.y = ggplot2::element_text(angle = 0),
                   axis.text.y = ggplot2::element_text(size = 6),
                   axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 4),
                   legend.position = "bottom")
  if (length(cluster_boundaries) > 0)
    p <- p + ggplot2::geom_vline(xintercept = cluster_boundaries - 0.5, colour = "grey30", linewidth = 0.3)
}

# manual SV color legend (square swatches); bottom when flipped, else right (ggplot default)
legend_levels <- c("reference-like" = COL[["ref"]], "DEL" = COL[["DEL"]], "INV" = COL[["INV"]],
                   "INS/NOVEL" = COL[["INS_NOVEL"]], "INS/DUP" = COL[["INS_DUP"]],
                   "DUP (shade=absolute CN)" = dup_ramp(3)[3], "multiallelic" = multi_ramp(3)[3])
leg <- data.frame(lab = factor(names(legend_levels), levels = names(legend_levels)), x = NA_real_, y = NA_real_)
p <- p + ggplot2::geom_point(data = leg, ggplot2::aes(x = x, y = y, shape = lab), na.rm = TRUE) +
  ggplot2::scale_shape_manual(name = "call", values = rep(15, length(legend_levels))) +
  ggplot2::guides(shape = ggplot2::guide_legend(nrow = if (opts$flip) 1 else NULL,
                                                override.aes = list(color = unname(legend_levels), size = 5)))

# figure size (unless given): the variant axis grows with variant count, the haplotype axis with
# haplotype count; default orientation puts variants on X, flipped puts them on Y.
var_in <- max(7, min(28, 2.5 + n_var * 0.55)); hap_in <- max(4, min(40, 1.8 + n_paths / 11))
if (is.na(opts$width))  opts$width  <- if (!opts$flip) var_in else hap_in
if (is.na(opts$height)) opts$height <- if (!opts$flip) hap_in else var_in

png_path <- paste0(opts$out, ".png"); pdf_path <- paste0(opts$out, ".pdf")
ggplot2::ggsave(png_path, p, width = opts$width, height = opts$height, units = "in", dpi = 180, limitsize = FALSE)
ggplot2::ggsave(pdf_path, p, width = opts$width, height = opts$height, units = "in", limitsize = FALSE)
cat("Wrote:", png_path, "\n"); cat("Wrote:", pdf_path, "\n")
