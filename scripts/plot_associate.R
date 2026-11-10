#!/usr/bin/env Rscript
# associate runs along the reference: -log10 p, and -log10 p given the signals. See --help for usage.
suppressWarnings(suppressMessages({ library(ggplot2); library(ggrepel) }))

args <- commandArgs(trailingOnly = TRUE)
usage <- function(status = 0) {
  cat(paste(c(
    "plot_associate.R - associate runs along the reference",
    "",
    "Usage: Rscript plot_associate.R [options] [<label>=]<associate prefix> ...",
    "",
    "Top: -log10 p over each row's reference span, with the run's thresholds (selection in dark red).",
    "Bottom (variant unit): -log10 p given the signals. Several runs: a block each, on its own y scale.",
    "",
    "Options:",
    "  --top <N>          rows named per panel and run (default 5; 1 with several runs)",
    "  --layers <a,b>     only these sources (default all)",
    "  --out <path>       .png or .pdf (default <prefix>.region.png; associate.region.png for several)",
    "  --dpi <N>          resolution (default 150)",
    "  --width, --height  size in inches (default from the runs and rows)",
    "  -h, --help         show this help"), collapse = "\n"), "\n")
  quit(status = status)
}
if (length(args) == 0 || any(args %in% c("-h", "--help"))) usage(0)
valued <- c("--out", "--layers", "--dpi", "--width", "--height", "--top")
get <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i == length(args)) stop("missing value after ", flag)
  args[i + 1]
}
taken <- unlist(lapply(valued, function(f) { i <- match(f, args); if (is.na(i)) integer(0) else c(i, i + 1) }))
runs <- if (length(taken)) args[-taken] else args
if (any(startsWith(runs, "--"))) stop("unknown option: ", runs[startsWith(runs, "--")][1])
if (length(runs) == 0) usage(1)
top_n_given <- get("--top")
dpi <- as.numeric(get("--dpi", "150"))
keep_layers <- if (is.null(get("--layers"))) NULL else strsplit(get("--layers"), ",")[[1]]

read <- function(p) read.delim(p, colClasses = "character", check.names = FALSE, na.strings = c("NA", "."))
# Break marks: a white cut across the panel at the gap, and a stroke over each cut end of the frame.
GeomBreakMarks <- ggproto("GeomBreakMarks", Geom, required_aes = "y", default_aes = aes(colour = "grey20"),
  draw_key = draw_key_blank,
  draw_panel = function(data, panel_params, coord) {
    d <- coord$transform(transform(data, x = 0), panel_params)
    lo <- min(d$y); hi <- max(d$y)
    cut <- grid::rectGrob(x = grid::unit(0, "npc") - grid::unit(4, "pt"), y = grid::unit(lo, "npc"),
                          width = grid::unit(1, "npc") + grid::unit(8, "pt"), height = grid::unit(hi - lo, "npc"),
                          just = c("left", "bottom"), gp = grid::gpar(fill = "white", col = NA))
    side <- rep(c(0, 1), each = 2); yy <- rep(c(lo, hi), 2)
    strokes <- grid::segmentsGrob(x0 = grid::unit(side, "npc") - grid::unit(5, "pt"), x1 = grid::unit(side, "npc") + grid::unit(5, "pt"),
                                  y0 = grid::unit(yy, "npc") - grid::unit(3, "pt"), y1 = grid::unit(yy, "npc") + grid::unit(3, "pt"),
                                  gp = grid::gpar(col = d$colour[1], lwd = 1.4))
    grid::gTree(children = grid::gList(cut, strokes))
  })
geom_break_marks <- function(mapping = NULL, data = NULL, ...)
  layer(geom = GeomBreakMarks, mapping = mapping, data = data, stat = "identity", position = "identity",
        show.legend = FALSE, inherit.aes = FALSE, params = list(...))
source_order <- c("variant", "snv", "svim", "kmer", "graph")
threshold_keys <- c(Bonferroni = "bonferroni_threshold", permutation = "permutation_threshold")
rows <- list(); lines <- list(); labels <- character(0); prefixes <- character(0); info <- list()
for (r in runs) {
  label <- NA
  if (grepl("^[^/=]+=", r)) { label <- sub("=.*", "", r); r <- sub("^[^=]*=", "", r) }
  if (is.na(label)) label <- basename(r)
  label <- toupper(label)
  if (label %in% labels) stop("two runs labelled ", label, ": give each a <label>=")
  labels <- c(labels, label); prefixes <- c(prefixes, r)
  a <- read(paste0(r, ".assoc.tsv")); s <- read(paste0(r, ".summary.tsv"))
  key <- function(k) { v <- s$value[s$key == k]; if (length(v) == 0 || is.na(v[1])) NA_character_ else v[1] }
  for (c in c("chrom", "pos", "end"))
    if (!c %in% names(a)) stop(r, ".assoc.tsv has no ", c, " column: re-run associate on describe -r sidecars")
  if (!is.null(keep_layers)) a <- a[a$layer %in% keep_layers, ]
  for (c in c("pos", "end", "p", "p_conditional", "stat", "n")) a[[c]] <- as.numeric(a[[c]])
  placed <- !is.na(a$pos) & !is.na(a$end) & !is.na(a$p)
  if (sum(!placed) > 0) message(label, ": ", sum(!placed), " row(s) with no place on the reference or no p left out")
  a <- a[placed, ]
  if (nrow(a) == 0) stop(label, ": no row has a place on the reference")
  chrom <- unique(a$chrom)
  if (length(chrom) > 1) stop(label, ": rows lie on several chromosomes (", paste(chrom, collapse = ", "), "); a run is drawn along one")
  # At the 1e-300 floor a t test's p is rebuilt from its statistic (df from the row's samples, or every sample for
  # lmm), so such rows keep their order; a score test's tail cannot be rebuilt and stays at the floor.
  ncov <- suppressWarnings(as.numeric(key("covariates"))); if (is.na(ncov)) ncov <- 0
  n_all <- suppressWarnings(as.numeric(key("samples_used")))
  a$y <- -log10(pmax(a$p, 1e-300))
  deep <- a$p <= 1e-300 & is.finite(a$stat) & a$p_method %in% c("t", "lmm")
  if (any(deep)) {
    df <- ifelse(a$p_method[deep] == "lmm" & !is.na(n_all), n_all, a$n[deep]) - ncov - 2
    a$y[deep] <- -(log(2) + pt(-abs(a$stat[deep]), df = pmax(1, df), log.p = TRUE)) / log(10)
  }
  a$yc <- -log10(pmax(a$p_conditional, 1e-300))
  a$run <- label
  a$x0 <- a$pos / 1e6; a$x1 <- a$end / 1e6
  rows[[label]] <- a
  # The run's thresholds; the one the selection used is marked.
  l <- data.frame(name = character(0), p = numeric(0))
  for (nm in names(threshold_keys)) {
    v <- suppressWarnings(as.numeric(key(threshold_keys[[nm]])))
    if (!is.na(v) && v > 0) l[nrow(l) + 1, ] <- list(nm, v)
  }
  entry <- suppressWarnings(as.numeric(key("cojo_entry_threshold")))
  l$selection <- !is.na(entry) & abs(l$p - entry) <= 1e-9 * entry
  if (!is.na(entry) && !any(l$selection)) l[nrow(l) + 1, ] <- list("selection", entry, TRUE)
  l$run <- rep(label, nrow(l))
  lines[[label]] <- l
  info[[label]] <- list(chrom = chrom, model = key("model"), unit = key("unit"), n = nrow(a),
                        variant = identical(key("unit"), "variant") && any(!is.na(a$yc)))
}
a <- do.call(rbind, rows)
lt <- do.call(rbind, lines)
a$run <- factor(a$run, labels)
layers <- c(intersect(source_order, unique(a$layer)), sort(setdiff(unique(a$layer), source_order)))
palette <- setNames(c("#b45f06", "#2e75b6", "#6aa84f", "#8e7cc3", "#999999", "#cc4125", "#45818e")[seq_along(layers)], layers)
a$layer <- factor(a$layer, layers)
a$signal <- a$cond_role %in% c("signal", "lead")
a$in_cs <- !is.na(a$cs)

panels <- c(top = "-log10 p", bottom = "-log10 p given the signals")
two_rows <- any(vapply(info, function(i) i$variant, logical(1)))
long <- transform(a, panel = "top", yy = y)
if (two_rows) {
  # variant-unit runs only: a feature run's p given its lead is another test
  variant_runs <- labels[vapply(info, function(i) i$variant, logical(1))]
  b <- a[!a$signal & !is.na(a$yc) & a$run %in% variant_runs, ]
  if (nrow(b)) long <- rbind(long, transform(b, panel = "bottom", yy = yc))
}
long$panel <- factor(long$panel, names(panels))
long$wide <- long$x1 > long$x0
# Several runs: one row, a block [k - 1, k] of the x axis per run, as chromosomes in a genome-wide plot.
many <- length(labels) > 1
top_n <- as.integer(if (is.null(top_n_given)) (if (many) 1 else 5) else top_n_given)
if (many) {
  block <- function(x, r) {
    k <- match(as.character(r), labels)
    lo <- tapply(c(a$x0, a$x1), rep(as.character(a$run), 2), min)[as.character(r)]
    hi <- tapply(c(a$x0, a$x1), rep(as.character(a$run), 2), max)[as.character(r)]
    (k - 1) + 0.04 + 0.92 * ifelse(hi > lo, (x - lo) / (hi - lo), 0.5)
  }
  long$x0 <- block(long$x0, long$run); long$x1 <- block(long$x1, long$run)
}
long$xm <- (long$x0 + long$x1) / 2
# Several runs: each block maps its -log10 p to [0, 1] on its own scale, linear up to 0.95; when a few rows stand
# above three times its 99.5th percentile, the bulk and the lines take 0.62, a break 0.62-0.68, those few the rest.
if (many) {
  scale_of <- list()
  for (r in labels) for (pn in levels(long$panel)) {
    sel <- long$run == r & long$panel == pn
    if (!any(sel)) next
    v <- long$yy[sel]
    th <- -log10(lt$p[lt$run == r & (pn == "top" | lt$selection)])
    cut <- 3 * max(as.numeric(quantile(v, 0.995)), 1)
    hi <- v[v > cut]
    if (pn == "top" && length(hi) && any(v <= cut)) {
      b1 <- max(c(v[v <= cut], th)) * 1.15; b2 <- min(hi) * 0.95; tp <- max(hi) * 1.04
      f <- local({ b1 <- b1; b2 <- b2; tp <- tp
        function(y) ifelse(y <= b1, 0.62 * y / b1,
                           ifelse(y >= b2, 0.68 + 0.27 * (y - b2) / (tp - b2), 0.62 + 0.06 * (y - b1) / (b2 - b1))) })
      lo_t <- pretty(c(0, b1), n = 4); hi_t <- pretty(c(b2, tp), n = 2)
      ticks <- c(lo_t[lo_t <= b1 & 0.62 * lo_t / b1 <= 0.57], hi_t[hi_t >= b2 & hi_t <= tp])   # none just under the break
      gap <- TRUE
    } else {
      tp <- max(c(v, th)) * 1.08
      f <- local({ tp <- tp; function(y) 0.95 * y / tp })
      ticks <- pretty(c(0, tp), n = 4); ticks <- ticks[ticks <= tp]
      gap <- FALSE
    }
    scale_of[[paste(r, pn)]] <- list(f = f, ticks = ticks, gap = gap)
  }
  for (k in names(scale_of)) {
    sel <- paste(long$run, long$panel) == k
    long$yy[sel] <- scale_of[[k]]$f(long$yy[sel])
  }
  block_ticks <- do.call(rbind, lapply(names(scale_of), function(k) {
    r <- sub(" (top|bottom)$", "", k); pn <- sub(".* ", "", k)
    data.frame(k = match(r, labels), panel = factor(pn, names(panels)), value = scale_of[[k]]$ticks,
               d = scale_of[[k]]$f(scale_of[[k]]$ticks))
  }))
  block_gaps <- do.call(rbind, lapply(names(scale_of)[vapply(scale_of, function(s) s$gap, logical(1))], function(k) {
    r <- sub(" (top|bottom)$", "", k)
    data.frame(k = match(r, labels), panel = factor("top", names(panels)), lo = 0.62, hi = 0.68)
  }))
}

# Named rows: the strongest per run and panel; the top panel's first with its genes (bold for one run).
named <- do.call(rbind, lapply(split(long, list(long$run, long$panel), drop = TRUE), function(d) {
  d <- d[order(-d$yy), ][seq_len(min(top_n, nrow(d))), ]
  d$first <- !many & d$panel == "top" & seq_len(nrow(d)) == 1
  lead <- d$panel == "top" & seq_len(nrow(d)) == 1
  d$text <- ifelse(lead & !is.na(d$gene), paste0(d$feature_id, "\n", d$gene), d$feature_id)
  d
}))

# Threshold lines in the top row, the selection's also in the bottom row; their values go in the run's title.
lt$run <- factor(lt$run, labels)
lt$y <- -log10(lt$p)
lt$colour <- ifelse(lt$selection, "#a61c00", "grey35")
strip_text <- vapply(labels, function(r) {
  d <- lt[lt$run == r, ]
  d <- d[order(d$name), ]
  values <- paste0(d$name, " ", formatC(d$p, format = "e", digits = 1), ifelse(d$selection, " (selection)", ""),
                   collapse = " · ")
  paste0(r, " · ", info[[r]]$chrom, " (Mb)", if (nrow(d)) paste0("\n", values) else "")
}, "")
lines_long <- rbind(transform(lt, panel = "top"), if (two_rows) transform(lt[lt$selection, ], panel = "bottom"))
lines_long$panel <- factor(lines_long$panel, names(panels))
if (many) {
  # on the run's own scale; no line where the run has no rows
  has_scale <- paste(lines_long$run, lines_long$panel) %in% names(scale_of)
  lines_long <- lines_long[has_scale, ]
  for (i in seq_len(nrow(lines_long)))
    lines_long$y[i] <- scale_of[[paste(lines_long$run[i], lines_long$panel[i])]]$f(lines_long$y[i])
} else if (two_rows) {
  lines_long <- lines_long[lines_long$panel == "top" | lines_long$run %in% unique(long$run[long$panel == "bottom"]), ]
}
# Several runs: lines span their run's block, named with its selection threshold.
k_of <- match(as.character(lines_long$run), labels)
lines_long$x <- k_of - 1 + 0.01; lines_long$xend <- k_of - 0.01
block_text <- vapply(labels, function(r) {
  d <- lt[lt$run == r & lt$selection, ]
  paste0(r, " \u00b7 ", info[[r]]$chrom, if (nrow(d)) paste0("\nselection ", formatC(d$p[1], format = "e", digits = 1)) else "")
}, "")
shade <- if (many) expand.grid(k = seq(2, length(labels), by = 2), panel = factor(levels(long$panel)[levels(long$panel) %in% long$panel], names(panels))) else NULL
# One run: a broken axis when a few rows stand above three times the 99.5th percentile; linear below the break
# for the bulk and the lines, a short linear part above it for those few, each with its own ticks.
ytop <- long$yy[long$panel == "top"]
q <- as.numeric(quantile(ytop, 0.995))
high <- ytop[ytop > 3 * max(q, 1)]
broken <- !many && length(high) > 0 && any(ytop <= 3 * max(q, 1))
if (broken) {
  low_max <- max(c(ytop[ytop <= 3 * max(q, 1)], lt$y))
  b1 <- low_max * 1.15                       # top of the lower part
  b2 <- min(high) * 0.95                     # bottom of the upper part
  top_y <- max(high) * 1.04
  gap <- 0.08 * b1                           # the break, in lower-part units
  up <- 0.45 * b1 / (top_y - b2)             # the upper part takes about a third of the height
  fwd <- function(y) ifelse(y <= b1, y, ifelse(y >= b2, b1 + gap + (y - b2) * up, b1 + (y - b1) / (b2 - b1) * gap))
  inv <- function(t) ifelse(t <= b1, t, ifelse(t >= b1 + gap, b2 + (t - b1 - gap) / up, b1 + (t - b1) / gap * (b2 - b1)))
  ticks <- function(lim) {
    if (lim[2] <= b1) return(pretty(lim, n = 5))   # the bottom row
    lo <- pretty(c(0, b1), n = 5); lo <- lo[lo <= b1]
    hi <- pretty(c(b2, top_y), n = 3); hi <- hi[hi >= b2 & hi <= top_y]
    v <- c(lo, hi); v[v >= lim[1] & v <= lim[2]]
  }
  y_scale <- scale_y_continuous(trans = scales::trans_new("broken", fwd, inv), breaks = ticks,
                                expand = expansion(mult = c(0.02, 0.06)))
  marks <- data.frame(run = factor(labels, labels), panel = factor("top", names(panels)), y = c(b1, b2))
} else if (many) {
  y_scale <- scale_y_continuous(limits = c(-0.01, 1.02), breaks = NULL, expand = expansion(add = 0))
} else {
  y_scale <- scale_y_continuous(expand = expansion(mult = c(0.02, 0.15)))
}

fig <- ggplot(long) +
  # every other block shaded
  { if (many && nrow(shade)) geom_rect(data = shade, aes(xmin = k - 1, xmax = k, ymin = -Inf, ymax = Inf), inherit.aes = FALSE,
                                       fill = "grey95", colour = NA) } +
  # the break over the frame, under the rows
  { if (broken) geom_break_marks(data = marks, aes(y = y)) } +
  # several runs: each block's gridlines, ticks and break
  { if (many) geom_segment(data = block_ticks, aes(x = k - 1, xend = k, y = d, yend = d), inherit.aes = FALSE,
                           colour = "grey90", linewidth = 0.3) } +
  { if (many) geom_text(data = block_ticks, aes(x = k - 1 + 0.008, y = d, label = value), inherit.aes = FALSE,
                        hjust = 0, vjust = -0.3, size = 2.3, colour = "grey45") } +
  # an outer block's gap reaches past the panel edge, cutting the frame
  { if (many && !is.null(block_gaps)) geom_rect(data = block_gaps, aes(xmin = ifelse(k == 1, -0.008, k - 1),
                                                                       xmax = ifelse(k == length(labels), k + 0.008, k),
                                                                       ymin = lo, ymax = hi),
                                                inherit.aes = FALSE, fill = "white", colour = NA) } +
  { if (many && !is.null(block_gaps)) geom_segment(data = rbind(transform(block_gaps, x = k - 1, yb = lo), transform(block_gaps, x = k - 1, yb = hi),
                                                                transform(block_gaps, x = k, yb = lo), transform(block_gaps, x = k, yb = hi)),
                                                   aes(x = x - 0.012, xend = x + 0.012, y = yb - 0.012, yend = yb + 0.012),
                                                   inherit.aes = FALSE, colour = "grey20", linewidth = 0.6) } +
  { if (many) geom_segment(data = lines_long, aes(x = x, xend = xend, y = y, yend = y, linetype = name, colour = I(colour)),
                           linewidth = 0.5)
    else geom_hline(data = lines_long, aes(yintercept = y, linetype = name, colour = I(colour)), linewidth = 0.5) } +
  geom_segment(data = long[long$wide, ], aes(x = x0, xend = x1, y = yy, yend = yy, colour = layer), linewidth = 1.1) +
  geom_point(data = long[!long$wide, ], aes(x = x0, y = yy, colour = layer), size = 1.1, alpha = 0.7) +
  geom_point(data = long[long$in_cs & long$panel == "top", ], aes(x = xm, y = yy), shape = 21, size = 3,
             colour = "grey15", fill = NA, stroke = 0.6) +
  geom_label_repel(data = named, aes(x = xm, y = yy, label = text, fontface = ifelse(first, "bold", "plain")),
                   size = 2.4, colour = "grey10", fill = "white", label.size = 0.2, box.padding = 0.4,
                   min.segment.length = 0, segment.colour = "grey40", segment.size = 0.3,
                   arrow = arrow(length = unit(0.012, "npc"), type = "closed"), max.overlaps = Inf, seed = 1) +
  scale_colour_manual(values = palette, name = NULL, breaks = layers) +
  scale_linetype_manual(values = c(Bonferroni = "dotted", permutation = "dashed", selection = "solid"), name = NULL,
                        guide = guide_legend(override.aes = list(colour = "grey35"))) +
  { if (many) facet_grid(panel ~ ., switch = "y", labeller = labeller(panel = panels))
    else facet_grid(panel ~ run, scales = "free", switch = "y", labeller = labeller(panel = panels, run = setNames(strip_text, labels))) } +
  { if (many) scale_x_continuous(breaks = seq_along(labels) - 0.5, labels = block_text, expand = expansion(add = 0)) } +
  labs(x = NULL, y = NULL,
       title = if (length(labels) == 1) paste0(labels, " · ", info[[1]]$model, " · ", info[[1]]$unit, " unit · ",
                                              info[[1]]$n, " rows") else NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position = "top", panel.grid.minor = element_blank(), strip.placement = "outside",
        strip.background = element_blank(), strip.text = element_text(size = 9), plot.title = element_text(size = 10))
fig <- fig + y_scale
# The panel is exactly [0, number of runs], so the outer break strokes sit on the frame.
if (many) fig <- fig + theme(panel.grid.major.x = element_blank(), axis.ticks.x = element_blank()) +
  coord_cartesian(xlim = c(0, length(labels)), ylim = c(-0.01, 1.02), expand = FALSE, clip = "off")
# Several runs: the frame is a layer over the shading, its side edges cut at an outer block's break.
if (many) {
  K <- length(labels); y0 <- -0.01; y1 <- 1.02
  gap_at <- function(k) !is.null(block_gaps) && k %in% block_gaps$k
  frame <- do.call(rbind, lapply(levels(droplevels(long$panel)), function(pn) {
    edge <- function(x, k) if (pn == "top" && gap_at(k)) data.frame(x = x, y = c(y0, 0.68), yend = c(0.62, y1)) else
      data.frame(x = x, y = y0, yend = y1)
    rbind(data.frame(x = c(0, 0), xend = c(K, K), y = c(y0, y1), yend = c(y0, y1)),
          transform(edge(0, 1), xend = x), transform(edge(K, K), xend = x))[, c("x", "xend", "y", "yend")] |>
      transform(panel = factor(pn, names(panels)))
  }))
  fig$layers <- append(fig$layers, list(geom_segment(data = frame, aes(x = x, xend = xend, y = y, yend = yend),
                                                     inherit.aes = FALSE, colour = "grey20", linewidth = 0.6)),
                       after = if (nrow(shade)) 1 else 0)
  fig <- fig + theme(panel.border = element_blank(), panel.background = element_rect(fill = "white", colour = NA),
                     legend.key = element_rect(fill = "white", colour = NA))
}
if (broken)
  fig <- fig + coord_cartesian(clip = "off") +
    # the frame as the panel background, so the break's white cut interrupts it
    theme(panel.border = element_blank(), panel.background = element_rect(fill = "white", colour = "grey20", linewidth = 0.6),
          legend.key = element_rect(fill = "white", colour = NA))
out <- get("--out", if (length(prefixes) > 1) file.path(dirname(prefixes[1]), "associate.region.png") else
  paste0(prefixes[1], ".region.png"))
if (!grepl("\\.(png|pdf)$", out)) stop("--out must end in .png or .pdf")
width <- as.numeric(get("--width", as.character(if (many) max(10, 3 + 2.5 * length(labels)) else 10)))
height <- as.numeric(get("--height", if (two_rows) "7.5" else "5"))
ggsave(out, fig, width = width, height = height, dpi = dpi)
cat("wrote", out, "\n")
