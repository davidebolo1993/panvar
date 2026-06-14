#!/usr/bin/env Rscript

# Tube-map schematics for slides. A reference "spine" of node boxes along x with a few
# example haplotypes drawn as ribbons threading the boxes. Two modes, both built only
# from panvar outputs (no extra tools):
#
#   --mode panphorte   one figure per normalized bubble: BEFORE = the repeat unit drawn
#                      `copies` times; AFTER = a single REP box with a self-loop x copies.
#                      Input: <prefix>.panphorte.copies.tsv (approximate-mode panphorte).
#
#   --mode sv          how each called SV type departs from the reference in the graph.
#                      Reference spine (node_track.tsv, in_reference nodes) + one carrier
#                      ribbon per SV type (variant_paths.tsv sub_walk), faceted by type.
#                      Inputs: --vcf, --variant-paths, --node-track.
#
# Common: --out <prefix> (writes .png + .pdf), --bubble-id N, --width/--height.

suppressWarnings(suppressMessages({ ok <- requireNamespace("ggplot2", quietly = TRUE) }))
if (!ok) stop("needs ggplot2 (conda install -c conda-forge r-ggplot2)")

args <- commandArgs(trailingOnly = TRUE)
usage <- function(status = 0) {
  cat(paste(c(
    "Usage:",
    "  plot_graph_schematic.R --mode panphorte --copies <copies.tsv> --out <prefix> [--bubble-id N]",
    "  plot_graph_schematic.R --mode sv --vcf <call.vcf> --variant-paths <variant_paths.tsv> \\",
    "                         --node-track <node_track.tsv> --out <prefix> [--bubble-id N] [--variant-id ID]",
    "  [--max-haplotypes 4] [--width in] [--height in]"), collapse = "\n"), "\n")
  quit(status = status)
}
if (length(args) == 0 || any(args %in% c("-h", "--help"))) usage(0)
opt <- list(mode = NULL, copies = NULL, vcf = NULL, variant_paths = NULL, node_track = NULL,
            out = NULL, bubble_id = NULL, variant_id = NULL, max_haplotypes = 4,
            width = NA_real_, height = NA_real_)
i <- 1
while (i <= length(args)) {
  a <- args[[i]]; v <- function() { if (i + 1 > length(args)) stop(paste("missing value after", a)); args[[i + 1]] }
  if (a == "--mode") { opt$mode <- v(); i <- i + 2 }
  else if (a == "--copies") { opt$copies <- v(); i <- i + 2 }
  else if (a == "--vcf") { opt$vcf <- v(); i <- i + 2 }
  else if (a == "--variant-paths") { opt$variant_paths <- v(); i <- i + 2 }
  else if (a == "--node-track") { opt$node_track <- v(); i <- i + 2 }
  else if (a == "--out") { opt$out <- v(); i <- i + 2 }
  else if (a == "--bubble-id") { opt$bubble_id <- as.integer(v()); i <- i + 2 }
  else if (a == "--variant-id") { opt$variant_id <- v(); i <- i + 2 }
  else if (a == "--max-haplotypes") { opt$max_haplotypes <- as.integer(v()); i <- i + 2 }
  else if (a == "--width") { opt$width <- as.numeric(v()); i <- i + 2 }
  else if (a == "--height") { opt$height <- as.numeric(v()); i <- i + 2 }
  else if (startsWith(a, "-")) stop(paste("unknown option:", a)) else i <- i + 1
}
if (is.null(opt$mode) || is.null(opt$out)) usage(1)
open_input <- function(p) if (grepl("\\.gz$", p, ignore.case = TRUE)) gzfile(p, open = "rt") else file(p, open = "rt")
read_tsv <- function(p) { con <- open_input(p); on.exit(close(con)); read.delim(con, sep = "\t", header = TRUE, check.names = FALSE, quote = "", colClasses = "character") }
short_label <- function(s) sub("[:#].*$", "", s)  # trim "#hap..:coords" tails for legibility
sqrtw <- function(bp) { w <- sqrt(pmax(as.numeric(bp), 1)); w / max(w, na.rm = TRUE) }  # 0..1 relative widths
save_plot <- function(p, w, h) {
  d <- dirname(opt$out); if (!identical(d, ".") && !dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(paste0(opt$out, ".png"), p, width = w, height = h, units = "in", dpi = 180, limitsize = FALSE)
  ggplot2::ggsave(paste0(opt$out, ".pdf"), p, width = w, height = h, units = "in", limitsize = FALSE)
  cat("Wrote:", paste0(opt$out, ".png"), "\n"); cat("Wrote:", paste0(opt$out, ".pdf"), "\n")
}

# ---------------------------------------------------------------- panphorte mode
if (opt$mode == "panphorte") {
  if (is.null(opt$copies)) stop("--mode panphorte needs --copies <prefix>.panphorte.copies.tsv")
  cp <- read_tsv(opt$copies)
  need <- c("path_name", "bubble_id", "copies", "unit_bp", "orientations", "n_long", "n_short")
  if (!all(need %in% names(cp))) stop(paste("copies.tsv missing columns:", paste(setdiff(need, names(cp)), collapse = ",")))
  cp$copies <- as.integer(cp$copies); cp$unit_bp <- as.numeric(cp$unit_bp)
  cp$n_long <- as.integer(cp$n_long); cp$n_short <- as.integer(cp$n_short)
  bid <- if (!is.null(opt$bubble_id)) opt$bubble_id else as.integer(cp$bubble_id[which.max(cp$copies)])
  cp <- cp[as.integer(cp$bubble_id) == bid, , drop = FALSE]
  if (nrow(cp) == 0) stop("no copies rows for bubble ", bid)

  # pick up to N example haplotypes spanning the copy-number range (distinct copy counts)
  ord <- order(cp$copies)
  pick <- unique(round(seq(1, nrow(cp), length.out = min(opt$max_haplotypes, nrow(cp)))))
  ex <- cp[ord[pick], , drop = FALSE]
  ex <- ex[!duplicated(ex$path_name), , drop = FALSE]

  unit_w <- 0.9                       # box width for a full ("long") unit copy
  flank_w <- 0.6; gap <- 0.15
  boxes <- list(); loops <- list(); labs <- list()
  rows <- nrow(ex)
  for (ri in seq_len(rows)) {
    e <- ex[ri, ]; y <- rows - ri + 1
    # orientations are a per-copy +/- string (e.g. "+-+"), occasionally comma-joined
    ors <- strsplit(gsub(",", "", e$orientations), "")[[1]]; if (length(ors) != e$copies) ors <- rep("+", e$copies)
    lab <- sprintf("%s (%d copies)", short_label(e$path_name), e$copies)
    labs[[length(labs) + 1]] <- data.frame(x = -0.2, y = y, text = lab, hjust = 1)
    for (phase in c("before", "after")) {
      x <- 0
      # left flank
      boxes[[length(boxes) + 1]] <- data.frame(xmin = x, xmax = x + flank_w, y = y, phase = phase, fill = "flank", glab = "")
      x <- x + flank_w + gap
      if (phase == "before") {
        n_short <- ifelse(is.na(e$n_short), 0, e$n_short); n_long <- e$copies - n_short
        for (k in seq_len(e$copies)) {
          w <- if (k > n_long) unit_w * 0.5 else unit_w     # short copies drawn narrower
          fl <- if (ors[k] == "-") "unit_rev" else "unit"
          boxes[[length(boxes) + 1]] <- data.frame(xmin = x, xmax = x + w, y = y, phase = phase, fill = fl,
                                                   glab = if (ors[k] == "-") "u-" else "u")
          x <- x + w + gap
        }
      } else {
        # single REP box with a self-loop annotated x copies
        boxes[[length(boxes) + 1]] <- data.frame(xmin = x, xmax = x + unit_w, y = y, phase = phase, fill = "rep", glab = "U")
        loops[[length(loops) + 1]] <- data.frame(x = x + unit_w / 2, y = y, phase = phase,
                                                 text = sprintf("×%d", e$copies))
        x <- x + unit_w + gap
      }
      # right flank
      boxes[[length(boxes) + 1]] <- data.frame(xmin = x, xmax = x + flank_w, y = y, phase = phase, fill = "flank", glab = "")
    }
  }
  B <- do.call(rbind, boxes); L <- do.call(rbind, loops); LB <- do.call(rbind, labs)
  B$phase <- factor(B$phase, levels = c("before", "after"), labels = c("before normalization", "after normalization"))
  L$phase <- factor(L$phase, levels = c("before", "after"), labels = c("before normalization", "after normalization"))
  pal <- c(flank = "#d9d9d9", unit = "#4292c6", unit_rev = "#f16913", rep = "#fec44f")
  p <- ggplot2::ggplot(B) +
    ggplot2::geom_rect(ggplot2::aes(xmin = xmin, xmax = xmax, ymin = y - 0.32, ymax = y + 0.32, fill = fill),
                       color = "#404040", linewidth = 0.3) +
    ggplot2::geom_text(ggplot2::aes(x = (xmin + xmax) / 2, y = y, label = glab), size = 2.6) +
    ggplot2::geom_text(data = L, ggplot2::aes(x = x, y = y + 0.42, label = text), size = 3, color = "#cc4c02") +
    ggplot2::geom_curve(data = L, ggplot2::aes(x = x - 0.3, y = y + 0.34, xend = x + 0.3, yend = y + 0.34),
                        curvature = -1.2, linewidth = 0.4, color = "#cc4c02",
                        arrow = grid::arrow(length = grid::unit(0.06, "in"))) +
    ggplot2::geom_text(data = LB, ggplot2::aes(x = x, y = y, label = text), hjust = 1, size = 3) +
    ggplot2::scale_fill_manual(values = pal, name = NULL,
                               labels = c(flank = "flank", unit = "unit copy (+)", unit_rev = "unit copy (-)", rep = "REP node")) +
    ggplot2::facet_wrap(~phase, ncol = 1) +
    ggplot2::coord_cartesian(xlim = c(-3.2, max(B$xmax) + 0.3)) +
    ggplot2::labs(title = sprintf("panphorte normalization - bubble %d", bid),
                  subtitle = "tandem unit copies (before) collapse to one REP self-loop (after)") +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                   strip.text = ggplot2::element_text(face = "bold", hjust = 0),
                   legend.position = "bottom")
  w <- if (is.na(opt$width)) max(8, min(24, 4 + max(B$xmax) / 2)) else opt$width
  h <- if (is.na(opt$height)) max(4, 1.2 + rows * 1.1) else opt$height
  save_plot(p, w, h)
  quit(status = 0)
}

# ---------------------------------------------------------------------- sv mode
if (opt$mode == "sv") {
  for (r in c("vcf", "variant_paths", "node_track")) if (is.null(opt[[r]])) usage(1)
  info_get <- function(info, key) {
    m <- regmatches(info, regexpr(paste0("(^|;)", key, "=[^;]*"), info))
    if (length(m) == 0) return(NA_character_); sub(paste0("^;?", key, "="), "", m)
  }
  parse_walk <- function(s) {
    if (is.na(s) || s == "") return(list(node = character(), rev = logical()))
    toks <- regmatches(s, gregexpr("[<>][^<>]+", s))[[1]]
    list(node = sub("^[<>]", "", toks), rev = substr(toks, 1, 1) == "<")
  }
  nt <- read_tsv(opt$node_track)
  nt$bubble_id <- as.integer(nt$bubble_id); nt$order <- as.integer(nt$order); nt$length_bp <- as.numeric(nt$length_bp)
  vp <- read_tsv(opt$variant_paths); vp$bubble_id <- as.integer(vp$bubble_id)
  vcon <- open_input(opt$vcf); vlines <- readLines(vcon); close(vcon)
  recs <- vlines[!grepl("^#", vlines)]
  # collect bubble records: id, bubble, svtype, and (for DUP) reference + max carrier CN
  V <- lapply(recs, function(ln) {
    f <- strsplit(ln, "\t")[[1]]
    cn <- suppressWarnings(as.integer(sub("^[^:]*:", "", f[-(1:9)])))
    list(id = f[3], info = f[8], bubble = suppressWarnings(as.integer(info_get(f[8], "BUBBLE_ID"))),
         svt = info_get(f[8], "SVTYPE"),
         ref_cn = suppressWarnings(as.integer(info_get(f[8], "REF_CN"))),
         cn_max = suppressWarnings(max(cn[is.finite(cn)], na.rm = TRUE)))
  })
  V <- Filter(function(r) !is.na(r$bubble) && !is.na(r$svt), V)
  bid <- if (!is.null(opt$bubble_id)) opt$bubble_id else {
    tb <- sort(table(sapply(V, function(r) r$bubble)), decreasing = TRUE); as.integer(names(tb)[1]) }
  V <- Filter(function(r) r$bubble == bid, V)
  if (length(V) == 0) stop("no SV records for bubble ", bid)
  # choose one representative variant per SV type (or a single --variant-id)
  if (!is.null(opt$variant_id)) V <- Filter(function(r) r$id == opt$variant_id, V)
  chosen <- list(); seen <- character()
  for (r in V) if (!(r$svt %in% seen)) { chosen[[length(chosen) + 1]] <- r; seen <- c(seen, r$svt) }

  # per-node bp for ALL inside nodes of the bubble (spine + off-spine event nodes)
  nb <- nt[nt$bubble_id == bid, , drop = FALSE]
  len_bp <- setNames(nb$length_bp, nb$node_id)
  kb <- function(bp) ifelse(bp >= 1000, sprintf("%.1f kb", bp / 1000), sprintf("%d bp", as.integer(round(bp))))
  W <- function(bp) sqrt(pmax(bp, 1))                         # sqrt-compressed segment width

  # Build a local, type-aware figure per SV: flanks for context, then the event region.
  # reference row (y=1) and carrier row (y=0) start at the same left flank so the carrier
  # being shorter (DEL) or longer (INS/DUP) is the visible signal.
  FLW <- 6                                                    # flank box width (in sqrt-bp units)
  boxes <- list(); labs <- list(); loops <- list()
  add <- function(x0, w, y, facet, fill, lab = "") boxes[[length(boxes) + 1]] <<-
    data.frame(x = x0 + w / 2, w = w, y = y, facet = facet, fill = fill, lab = lab, stringsAsFactors = FALSE)
  for (r in chosen) {
    facet <- sprintf("%s  (%s)", r$svt, r$id)
    ev <- info_get(r$info, "EVENT_NODES"); ev <- if (is.na(ev)) character() else strsplit(ev, ",")[[1]]
    ev_bp <- sum(as.numeric(len_bp[ev]), na.rm = TRUE); if (!is.finite(ev_bp) || ev_bp <= 0) ev_bp <- 1
    ew <- W(ev_bp)
    # reference row
    add(0, FLW, 1, facet, "flank")
    if (r$svt == "DEL") { add(FLW, ew, 1, facet, "ev_del", kb(ev_bp)); add(FLW + ew, FLW, 1, facet, "flank") }
    else if (r$svt == "INV") { add(FLW, ew, 1, facet, "kept"); add(FLW + ew, FLW, 1, facet, "flank") }
    else if (r$svt == "DUP") { add(FLW, ew, 1, facet, "kept", kb(ev_bp)); add(FLW + ew, FLW, 1, facet, "flank") }
    else { add(FLW, FLW, 1, facet, "flank") }                 # INS: no reference event region
    # carrier row
    add(0, FLW, 0, facet, "flank")
    if (r$svt == "DEL") { add(FLW, FLW, 0, facet, "flank")    # deletion -> straight to right flank (shorter)
    } else if (r$svt == "INS") { add(FLW, ew, 0, facet, "ins", kb(ev_bp)); add(FLW + ew, FLW, 0, facet, "flank")
    } else if (r$svt == "INV") { add(FLW, ew, 0, facet, "inv", paste0(kb(ev_bp), " (rev)")); add(FLW + ew, FLW, 0, facet, "flank")
    } else if (r$svt == "DUP") {
      cn <- if (is.finite(r$cn_max) && !is.na(r$cn_max)) r$cn_max else (if (is.finite(r$ref_cn)) r$ref_cn + 1 else 2)
      cn <- max(cn, 2); x <- FLW
      for (k in seq_len(cn)) { add(x, ew, 0, facet, "dup"); x <- x + ew }
      loops[[length(loops) + 1]] <- data.frame(x = FLW + ew * cn / 2, y = 0, facet = facet, text = sprintf("x%d", cn))
      add(x, FLW, 0, facet, "flank")
    }
  }
  B <- do.call(rbind, boxes); Lp <- if (length(loops)) do.call(rbind, loops) else NULL
  lab_seg <- B[nzchar(B$lab), , drop = FALSE]
  pal <- c(flank = "#d9d9d9", kept = "#bdbdbd", ins = "#238b45", inv = "#f16913", dup = "#2171b5", ev_del = "#cb181d")
  p <- ggplot2::ggplot(B) +
    ggplot2::geom_tile(ggplot2::aes(x = x, y = y, width = w, height = 0.5, fill = fill),
                       color = "#404040", linewidth = 0.25) +
    ggplot2::geom_text(data = lab_seg, ggplot2::aes(x = x, y = y, label = lab), size = 2.4) +
    ggplot2::scale_fill_manual(values = pal, name = NULL,
                               breaks = c("flank", "ev_del", "ins", "inv", "dup"),
                               labels = c(flank = "flank / reference", ev_del = "deleted region (DEL)",
                                          ins = "inserted (INS)", inv = "inverted (INV)",
                                          dup = "duplicated copy (DUP)")) +
    ggplot2::scale_y_continuous(breaks = c(0, 1), labels = c("carrier", "reference"), limits = c(-0.6, 1.6)) +
    ggplot2::facet_wrap(~facet, ncol = 1, scales = "free_x") +
    ggplot2::labs(title = sprintf("Called SV types in the graph - bubble %d", bid),
                  subtitle = "reference (top) vs one carrier haplotype (bottom), zoomed to the event; width = sqrt(bp)",
                  x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                   strip.text = ggplot2::element_text(face = "bold", hjust = 0),
                   panel.grid = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank(),
                   legend.position = "bottom")
  if (!is.null(Lp)) p <- p +
    ggplot2::geom_text(data = Lp, ggplot2::aes(x = x, y = y + 0.42, label = text), size = 2.8, color = "#08519c")
  w <- if (is.na(opt$width)) max(8, min(26, 3 + max(B$x + B$w) / 3)) else opt$width
  h <- if (is.na(opt$height)) max(4, 1.2 + length(chosen) * 1.8) else opt$height
  save_plot(p, w, h)
  quit(status = 0)
}

stop("--mode must be 'panphorte' or 'sv'")
