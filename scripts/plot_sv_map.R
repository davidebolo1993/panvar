#!/usr/bin/env Rscript

# Per-bubble structural-variant colormap, in the style of plot_node_coverage_heatmap.R.
# Rows = haplotypes, x = ALL of the bubble's walked nodes (length-scaled tiles) ordered by
# the reference. Each cell is colored by the SV the haplotype carries at that node relative
# to the reference; nodes a haplotype traverses without a call are light grey; nodes it does
# not traverse are white.
#
# Inputs (all standard panvar outputs):
#   --vcf          panvar call region (or per-bubble) VCF        -> SV type/subtype/CN/carriers
#   --node-counts  panvar inspect <...>.bubble_<id>.node_counts.tsv  -> per-(path,node) traversal matrix
#   --node-track   panvar call <prefix>.node_track.tsv           -> reference-ordered, length-scaled x-axis
#   --bubble-id    bubble to plot
#   --out          output prefix (writes <out>.png and <out>.pdf)

usage <- function(status = 0) {
  msg <- c(
    "Usage:",
    "  plot_sv_map.R --vcf <region.vcf> --node-counts <bubble_N.node_counts.tsv[.gz]>",
    "                --node-track <call.node_track.tsv> --bubble-id <N> --out <prefix> [options]",
    "",
    "Options:",
    "  --vcf <path>            panvar call VCF (region or per-bubble)",
    "  --node-counts <path>    panvar inspect per-bubble node_counts.tsv[.gz]",
    "  --node-track <path>     panvar call node_track.tsv",
    "  --bubble-id <N>         bubble id to plot (required)",
    "  --out <prefix>          output prefix; writes <prefix>.png and <prefix>.pdf",
    "  --length-transform <m>  node-width scaling: raw, sqrt, or log1p (default: sqrt)",
    "  --max-paths <N>         keep at most N haplotypes that carry a call (plus the reference)",
    "  --width <inches>        figure width (default: auto)",
    "  --height <inches>       figure height (default: auto)",
    "  -h, --help              show this help"
  )
  cat(paste(msg, collapse = "\n"), "\n", sep = "")
  quit(status = status)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0 || any(args %in% c("-h", "--help"))) usage(0)

opts <- list(vcf = NULL, node_counts = NULL, node_track = NULL, bubble_id = NULL,
             out = NULL, length_transform = "sqrt", max_paths = 0,
             width = NA_real_, height = NA_real_)
i <- 1
while (i <= length(args)) {
  a <- args[[i]]
  val <- function(flag) { if (i + 1 > length(args)) stop(paste("Missing value after", flag), call. = FALSE); args[[i + 1]] }
  if (a == "--vcf") { opts$vcf <- val(a); i <- i + 2 }
  else if (a == "--node-counts") { opts$node_counts <- val(a); i <- i + 2 }
  else if (a == "--node-track") { opts$node_track <- val(a); i <- i + 2 }
  else if (a == "--bubble-id") { opts$bubble_id <- val(a); i <- i + 2 }
  else if (a == "--out") { opts$out <- val(a); i <- i + 2 }
  else if (a == "--length-transform") { opts$length_transform <- val(a); i <- i + 2 }
  else if (a == "--max-paths") { opts$max_paths <- as.integer(val(a)); i <- i + 2 }
  else if (a == "--width") { opts$width <- as.numeric(val(a)); i <- i + 2 }
  else if (a == "--height") { opts$height <- as.numeric(val(a)); i <- i + 2 }
  else if (startsWith(a, "-")) stop(paste("Unknown option:", a), call. = FALSE)
  else i <- i + 1
}
for (req in c("vcf", "node_counts", "node_track", "bubble_id", "out")) {
  if (is.null(opts[[req]])) usage(1)
}
if (!opts$length_transform %in% c("raw", "sqrt", "log1p")) {
  stop("--length-transform must be one of: raw, sqrt, log1p", call. = FALSE)
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("This plotting helper requires ggplot2. Install it with: conda install -y -c conda-forge r-ggplot2",
       call. = FALSE)
}
bubble_id <- as.integer(opts$bubble_id)

# Return a source read.delim/readLines can open AND close themselves (avoids leaked connections).
open_input <- function(path) if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path) else path

# ---- node_track: the reference-ordered, length-scaled x-axis for this bubble ----
nt <- read.delim(open_input(opts$node_track), sep = "\t", header = TRUE, check.names = FALSE,
                 quote = "", colClasses = "character")
nt <- nt[as.integer(nt$bubble_id) == bubble_id, , drop = FALSE]
if (nrow(nt) == 0) stop(sprintf("bubble %d not found in node_track", bubble_id), call. = FALSE)
nt$order <- as.integer(nt$order)
nt$length_bp <- as.numeric(nt$length_bp)
nt$is_cn <- as.integer(nt$is_cn)
nt <- nt[order(nt$order), , drop = FALSE]
w <- switch(opts$length_transform, raw = nt$length_bp, sqrt = sqrt(pmax(nt$length_bp, 0)),
            log1p = log1p(pmax(nt$length_bp, 0)))
w[!is.finite(w) | w <= 0] <- min(w[w > 0], na.rm = TRUE) * 0.25  # keep zero/odd nodes visible
nt$w <- w
nt$x1 <- cumsum(nt$w)
nt$x0 <- nt$x1 - nt$w
nt$xmid <- (nt$x0 + nt$x1) / 2
node_ids <- nt$node_id
is_cn <- setNames(nt$is_cn == 1, node_ids)

# ---- node_counts: per-(path,node) traversal (total:forward:reverse) ----
nc <- read.delim(open_input(opts$node_counts), sep = "\t", header = TRUE, check.names = FALSE,
                 quote = "", comment.char = "")
nc_node_cols <- grep("^node\\.", names(nc), value = TRUE)
nc_ids <- sub("^node\\.", "", nc_node_cols)
total_at <- function(cells) suppressWarnings(as.numeric(sub(":.*$", "", cells)))
M <- matrix(0, nrow = nrow(nc), ncol = length(node_ids), dimnames = list(nc$path_name, node_ids))
common <- intersect(node_ids, nc_ids)
for (nid in common) M[, nid] <- total_at(nc[[paste0("node.", nid)]])
M[is.na(M)] <- 0

# ---- VCF: per-record SV type/subtype/CN and per-sample GT for this bubble ----
vlines <- readLines(open_input(opts$vcf))
hdr <- vlines[grepl("^#CHROM", vlines)][1]
samples <- strsplit(hdr, "\t")[[1]][-(1:9)]
recs <- vlines[!grepl("^#", vlines)]
info_get <- function(info, key) {
  m <- regmatches(info, regexpr(paste0("(^|;)", key, "=[^;]*"), info))
  if (length(m) == 0) return(NA_character_)
  sub(paste0("^;?", key, "="), "", m)
}
records <- list()
for (ln in recs) {
  f <- strsplit(ln, "\t")[[1]]
  info <- f[8]
  if (is.na(info_get(info, "BUBBLE_ID")) || as.integer(info_get(info, "BUBBLE_ID")) != bubble_id) next
  gtcn <- f[-(1:9)]
  gt <- sub(":.*$", "", gtcn)
  cn <- suppressWarnings(as.integer(sub("^[^:]*:", "", gtcn)))
  records[[length(records) + 1]] <- list(
    svtype = info_get(info, "SVTYPE"),
    subtype = info_get(info, "INS_SUBTYPE"),
    nodes = strsplit(info_get(info, "EVENT_NODES"), ",")[[1]],
    gt = setNames(gt, samples), cn = setNames(cn, samples))
}
if (length(records) == 0) stop(sprintf("no VCF records for bubble %d", bubble_id), call. = FALSE)

# Samples that traverse the bubble (GT != '.') and those carrying any call.
traverses <- rep(FALSE, length(samples)); names(traverses) <- samples
carriers_any <- rep(FALSE, length(samples)); names(carriers_any) <- samples
for (r in records) {
  traverses[r$gt != "."] <- TRUE
  carriers_any[r$gt == "1"] <- TRUE
}
keep <- names(traverses)[traverses]
if (opts$max_paths > 0 && length(keep) > opts$max_paths) {
  carry <- keep[carriers_any[keep]]
  if (length(carry) > opts$max_paths) carry <- carry[seq_len(opts$max_paths)]
  keep <- carry
}
keep <- keep[keep %in% rownames(M)]
if (length(keep) == 0) stop("no traversing haplotypes with node-count rows", call. = FALSE)

# ---- category + per-cell color ----
COL <- c(none_white = "#ffffff", invariant = "#e0e0e0", DEL = "#cb181d", INV = "#f16913",
         INS_NOVEL = "#238b45", INS_DUP = "#6a51a3")
dup_ramp <- grDevices::colorRampPalette(c("#c6dbef", "#08306b"))
prio <- c(DEL = 1, INS_NOVEL = 2, INS_DUP = 2, INV = 3, DUP = 4)  # higher wins

ns <- length(keep); nn <- length(node_ids)
cat_mat <- matrix("", ns, nn, dimnames = list(keep, node_ids))
pri_mat <- matrix(0L, ns, nn, dimnames = list(keep, node_ids))
cn_mat <- matrix(NA_integer_, ns, nn, dimnames = list(keep, node_ids))
set_cat <- function(rows, cols, label, cnvals = NULL) {
  if (length(rows) == 0 || length(cols) == 0) return(invisible())
  p <- prio[[label]]
  for (cl in cols) {
    sel <- rows[pri_mat[rows, cl] < p]
    if (length(sel) == 0) next
    cat_mat[sel, cl] <<- label; pri_mat[sel, cl] <<- p
    if (!is.null(cnvals)) cn_mat[sel, cl] <<- cnvals[sel]
  }
}
for (r in records) {
  carr <- intersect(keep, names(r$gt)[r$gt == "1"])
  if (length(carr) == 0) next
  evn <- intersect(r$nodes, node_ids)
  if (r$svtype == "DEL") {
    for (nid in evn) { absent <- carr[M[carr, nid] == 0]; set_cat(absent, nid, "DEL") }
  } else if (r$svtype == "INS") {
    lab <- if (!is.na(r$subtype) && r$subtype == "DUP") "INS_DUP" else "INS_NOVEL"
    for (nid in evn) { present <- carr[M[carr, nid] > 0]; set_cat(present, nid, lab) }
  } else if (r$svtype == "INV") {
    for (nid in evn) { present <- carr[M[carr, nid] > 0]; set_cat(present, nid, "INV") }
  } else if (r$svtype == "DUP") {
    cncols <- intersect(node_ids, c(r$nodes, names(is_cn)[is_cn]))
    cncols <- intersect(cncols, r$nodes)
    set_cat(carr, cncols, "DUP", cnvals = r$cn)
  }
}

# Per-cell fill: overlay if categorized, else grey (traversed) / white (absent).
fill_for <- function(s, nid) {
  cat <- cat_mat[s, nid]
  if (cat == "DUP") {
    cnv <- cn_mat[s, nid]; if (is.na(cnv) || cnv < 1) cnv <- 1
    pal <- dup_ramp(max(2, dup_max)); return(pal[min(cnv, dup_max)])
  }
  if (cat != "") return(COL[[cat]])
  if (M[s, nid] > 0) return(COL[["invariant"]])
  COL[["none_white"]]
}
dup_max <- suppressWarnings(max(cn_mat[!is.na(cn_mat)], 1))
df <- expand.grid(s = keep, col = seq_len(nn), stringsAsFactors = FALSE)
df$node_id <- node_ids[df$col]
df$x0 <- nt$x0[df$col]; df$x1 <- nt$x1[df$col]; df$xmid <- nt$xmid[df$col]
df$y <- match(df$s, keep)
df$fill <- mapply(fill_for, df$s, df$node_id)

# ---- render ----
out_dir <- dirname(opts$out)
if (!identical(out_dir, ".") && !dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (is.na(opts$width)) opts$width <- max(8, min(28, 4 + max(nt$x1) / 120))
if (is.na(opts$height)) opts$height <- max(5, min(30, 2 + ns / 6))

# Legend via a dummy categorical layer (identity fill for the heatmap itself).
legend_levels <- c("not traversed" = COL[["none_white"]], "invariant" = COL[["invariant"]],
                   "DEL" = COL[["DEL"]], "INV" = COL[["INV"]], "INS/NOVEL" = COL[["INS_NOVEL"]],
                   "INS/DUP" = COL[["INS_DUP"]], "DUP (shade=CN)" = dup_ramp(3)[3])
sparse_ticks <- function(n, k) if (n <= 0) integer() else unique(pmax(1, pmin(n, round(seq(1, n, length.out = min(n, k))))))
xt <- sparse_ticks(nn, 30)
yt <- sparse_ticks(ns, 40)

p <- ggplot2::ggplot(df) +
  ggplot2::geom_rect(ggplot2::aes(xmin = x0, xmax = x1, ymin = y - 0.5, ymax = y + 0.5, fill = fill)) +
  ggplot2::scale_fill_identity() +
  ggplot2::scale_x_continuous(breaks = nt$xmid[xt], labels = node_ids[xt], expand = c(0, 0)) +
  ggplot2::scale_y_continuous(breaks = yt, labels = keep[yt], expand = c(0, 0),
                              trans = "reverse") +
  ggplot2::labs(
    title = sprintf("Bubble %d structural-variant map", bubble_id),
    subtitle = sprintf("%d haplotypes x %d nodes (x = reference order, width = %s bp)", ns, nn, opts$length_transform),
    x = "Bubble nodes (reference order, length-scaled)", y = "Haplotypes") +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(panel.grid = ggplot2::element_blank(),
                 plot.title = ggplot2::element_text(face = "bold"),
                 axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5),
                 axis.text.y = ggplot2::element_text(size = 4))

# Manual color legend.
leg <- data.frame(lab = factor(names(legend_levels), levels = names(legend_levels)),
                  x = NA_real_, y = NA_real_)
p <- p + ggplot2::geom_point(data = leg, ggplot2::aes(x = x, y = y, color = lab), na.rm = TRUE) +
  ggplot2::scale_color_manual(name = "SV", values = legend_levels, drop = FALSE) +
  ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(shape = 15, size = 5)))

png_path <- paste0(opts$out, ".png"); pdf_path <- paste0(opts$out, ".pdf")
ggplot2::ggsave(png_path, p, width = opts$width, height = opts$height, units = "in", dpi = 180, limitsize = FALSE)
ggplot2::ggsave(pdf_path, p, width = opts$width, height = opts$height, units = "in", limitsize = FALSE)
cat("Wrote:", png_path, "\n"); cat("Wrote:", pdf_path, "\n")
