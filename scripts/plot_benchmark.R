#!/usr/bin/env Rscript
# What call's output leaves of the haplotypes' variation, by cause, for one benchmark run or several.
# See --help for usage.
suppressWarnings(suppressMessages(library(ggplot2)))

args <- commandArgs(trailingOnly = TRUE)
usage <- function(status = 0) {
  cat(paste(c(
    "plot_benchmark.R - benchmark runs as one figure",
    "",
    "Usage: Rscript plot_benchmark.R [options] [<label>=]<prefix>[:<call min bp>] ...",
    "",
    "Bars: share of the differences present still left, by cause, for call's events and the VCF; all",
    "differences and those at call's threshold or more. Above --max-grouped runs, one bar per run.",
    "",
    "Options:",
    "  --call-min-bp <N>  call's --min-sv-bp for runs without :<N> (default 50)",
    "  --only-events, --only-bases, --only-haplotype, --only-sites",
    "                     keep one row or column of panels",
    "  --per-row <N>      runs per row (default 50; 25 with one bar per run)",
    "  --max-grouped <N>  above N runs, one bar per run (default 10)",
    "  --stage <s>        one bar per run: call or vcf (default vcf)",
    "  --size <s>         one bar per run: all or above (default above)",
    "  --out <path>       .png or .pdf; the numbers go beside it as .tsv (default <prefix>.png;",
    "                     benchmark.png for several)",
    "  --dpi <N>          resolution (default 150)",
    "  --width, --height  size in inches (default from the runs and rows)",
    "  -h, --help         show this help"), collapse = "\n"), "\n")
  quit(status = status)
}
if (length(args) == 0 || any(args %in% c("-h", "--help"))) usage(0)
valued <- c("--out", "--call-min-bp", "--per-row", "--dpi", "--width", "--height", "--max-grouped", "--stage", "--size")
switches <- c("--only-events", "--only-bases", "--only-haplotype", "--only-sites")
get <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i == length(args)) stop("missing value after ", flag)
  args[i + 1]
}
taken <- unlist(lapply(valued, function(f) { i <- match(f, args); if (is.na(i)) integer(0) else c(i, i + 1) }))
taken <- c(taken, which(args %in% switches))
runs <- if (length(taken)) args[-taken] else args
bad <- runs[startsWith(runs, "--")]
if (length(bad)) stop("unknown option: ", bad[1])
if (length(runs) == 0) usage(1)
call_min_default <- as.numeric(get("--call-min-bp", "50"))
per_row <- as.integer(get("--per-row", "50"))
dpi <- as.numeric(get("--dpi", "150"))
measures <- c(events = "share of the differences", bases = "share of their bases")
if ("--only-events" %in% args) measures <- measures["events"]
if ("--only-bases" %in% args) measures <- measures["bases"]
scopes <- c(haplotype = "whole haplotype", sites = "inside the sites")
if ("--only-haplotype" %in% args) scopes <- scopes["haplotype"]
if ("--only-sites" %in% args) scopes <- scopes["sites"]
if ("--only-events" %in% args && "--only-bases" %in% args) stop("--only-events and --only-bases exclude each other")
if ("--only-haplotype" %in% args && "--only-sites" %in% args) stop("--only-haplotype and --only-sites exclude each other")

causes <- c("incomplete", "unresolved", "not_carrier", "no_record", "outside_bubbles")
reasons <- c("copy_count", "no_sequence", "imprecise", "representative", "other")
part_levels <- c(paste("missed:", causes), paste("lost in VCF:", reasons))
palette <- setNames(c("#1f4e79", "#2e75b6", "#6fa8dc", "#9fc5e8", "#cfe2f3",
                      "#b45f06", "#e69138", "#f6b26b", "#cc4125", "#f9cb9c"), part_levels)
read <- function(p) read.delim(p, colClasses = "character", check.names = FALSE)

# One row per run, scope, size class, stage and part: what is left and what was present.
rows <- list(); labels <- character(0); prefixes <- character(0)
for (r in runs) {
  label <- NA
  if (grepl("^[^/=]+=", r)) { label <- sub("=.*", "", r); r <- sub("^[^=]*=", "", r) }
  t_call <- call_min_default
  if (grepl(":[0-9]+$", r)) { t_call <- as.numeric(sub(".*:", "", r)); r <- sub(":[0-9]+$", "", r) }
  if (is.na(label)) label <- basename(r)
  # names are shown in upper case, so two that differ only in case would be one
  if (toupper(label) %in% toupper(labels)) stop("two runs labelled ", label, ": give each a <label>=")
  labels <- c(labels, label); prefixes <- c(prefixes, r)
  s <- read(paste0(r, ".summary.tsv"))
  floor_bp <- as.numeric(s[[2]][s[[1]] == "min_sv_bp"])
  if (!"left_by_vcf_events" %in% s[[1]]) stop(r, ".summary.tsv is from an older benchmark: re-run benchmark")
  if (floor_bp >= t_call) message(label, ": benchmarked at ", floor_bp, " bp, not below call's ", t_call,
                                   " bp: both sizes show the same differences")
  present <- read(paste0(r, ".present.tsv"))
  if (!"bubbles" %in% names(present)) stop(r, ".present.tsv has no bubbles column: re-run benchmark")
  events <- read(paste0(r, ".events.tsv"))
  part <- ifelse(events$reason %in% causes, paste("missed:", events$reason), paste("lost in VCF:", events$reason))
  for (sc in names(scopes)) for (sz in c("all", "above")) {
    keep_p <- (sc == "haplotype" | present$bubbles != ".") & (sz == "all" | as.numeric(present$size_bp) >= t_call)
    keep_e <- (sc == "haplotype" | events$bubbles != ".") & (sz == "all" | as.numeric(events$size_bp) >= t_call)
    n_p <- sum(keep_p); bp_p <- sum(as.numeric(present$size_bp[keep_p]))
    for (st in c("call", "VCF")) {
      k <- keep_e & events$stage == (if (st == "call") "missed_by_call" else "left_by_vcf")
      if (!any(k)) {
        rows[[length(rows) + 1]] <- data.frame(run = label, call_min_bp = t_call, benchmark_min_bp = floor_bp,
          scope = sc, size = sz, stage = st, part = NA, left = 0, left_bp = 0, present = n_p, present_bp = bp_p)
        next
      }
      a <- aggregate(list(left = rep(1, sum(k)), left_bp = as.numeric(events$size_bp[k])), list(part = part[k]), sum)
      rows[[length(rows) + 1]] <- data.frame(run = label, call_min_bp = t_call, benchmark_min_bp = floor_bp,
        scope = sc, size = sz, stage = st, part = a$part, left = a$left, left_bp = a$left_bp,
        present = n_p, present_bp = bp_p)
    }
  }
}
d <- do.call(rbind, rows)
unknown <- setdiff(na.omit(d$part), part_levels)
if (length(unknown)) stop("unknown cause in events.tsv: ", unknown[1])

out <- get("--out", if (length(prefixes) > 1) file.path(dirname(dirname(prefixes[1])), "benchmark.png")
                    else paste0(prefixes[1], ".png"))
if (!grepl("\\.(png|pdf)$", out)) stop("--out must end in .png or .pdf")
table_out <- sub("\\.[a-z]+$", ".tsv", out)
write.table(d, table_out, sep = "\t", quote = FALSE, row.names = FALSE, na = ".")

# Many runs: one panel, one bar per run, for the measure, scope, rebuild and size class asked for.
max_grouped <- as.integer(get("--max-grouped", "10"))
if (length(labels) > max_grouped) {
  m <- if ("--only-bases" %in% args) "bases" else "events"
  sc <- if ("--only-sites" %in% args) "sites" else "haplotype"
  st <- get("--stage", "VCF")
  if (tolower(st) == "vcf") st <- "VCF"
  if (!st %in% c("call", "VCF")) stop("--stage must be call or vcf")
  sz <- get("--size", "above")
  if (!sz %in% c("all", "above")) stop("--size must be all or above")
  one <- d[d$scope == sc & d$size == sz & d$stage == st, ]
  one$n_left <- if (m == "events") one$left else one$left_bp
  one$of <- if (m == "events") one$present else one$present_bp
  one$value <- ifelse(one$of > 0, 100 * one$n_left / one$of, 0)
  one$part <- factor(one$part, part_levels)
  one$run <- factor(toupper(one$run), toupper(labels))
  # runs per row: --per-row, else 25
  row_n <- as.integer(get("--per-row", "25"))
  one$row <- (match(one$run, toupper(labels)) - 1) %/% row_n
  tot <- aggregate(list(n_left = one$n_left, value = one$value), list(run = one$run, of = one$of, row = one$row), sum)
  unit <- if (m == "events") "" else if (max(tot$of) < 1e4) "bp" else if (max(tot$of) < 1e8) "kb" else "Mb"
  div <- c(bp = 1, kb = 1e3, Mb = 1e6)
  fmt <- function(v) {
    if (m == "events") return(as.character(v))
    x <- v / div[unit]
    ifelse(x == 0, "0", ifelse(x >= 100 | unit == "bp", format(round(x), big.mark = ",", trim = TRUE), sprintf("%.1f", x)))
  }
  tot$text <- paste0(fmt(tot$n_left), " / ", fmt(tot$of))
  t_names <- unique(d$call_min_bp)
  size_text <- if (sz == "all") "all differences" else if (length(t_names) == 1) paste0("differences of ", t_names, " bp or more") else "differences of call's threshold or more"
  title <- paste0("share of the ", if (m == "events") "differences" else paste0("bases (", unit, ")"), " still left \u00b7 ",
                  c(haplotype = "whole haplotype", sites = "inside the sites")[sc], " \u00b7 ",
                  if (st == "call") "call's own events" else "the VCF", " \u00b7 ", size_text)
  segs <- one[!is.na(one$part) & one$value > 0, ]
  fig <- ggplot() +
    # every run first, in order, so a run with nothing left keeps its place
    geom_blank(data = tot, aes(x = run, y = 0)) +
    geom_col(data = segs, aes(x = run, y = value, fill = part), width = 0.7, colour = "white", linewidth = 0.25) +
    geom_text(data = tot, aes(x = run, y = value, label = text), angle = 90, hjust = -0.08, size = 2.3, colour = "grey25") +
    scale_fill_manual(values = palette, breaks = part_levels, drop = TRUE, name = NULL) +
    scale_x_discrete(drop = TRUE) +
    facet_wrap(~row, ncol = 1, scales = "free_x") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = "% still left", title = title) +
    theme_bw(base_size = 10) +
    theme(legend.position = "top", panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), plot.title = element_text(size = 10),
          strip.text = element_blank(), strip.background = element_blank(), panel.spacing.y = unit(0.35, "in")) +
    guides(fill = guide_legend(ncol = 5))
  if (any(tot$value > 100)) fig <- fig + geom_hline(yintercept = 100, linetype = "dashed", colour = "grey50")
  rows <- max(one$row) + 1
  # wide enough for the title and the legend however few runs a row holds
  width <- as.numeric(get("--width", as.character(max(9, 2 + 0.25 * min(row_n, length(labels))))))
  height <- as.numeric(get("--height", as.character(1.5 + 3 * rows)))
  ggsave(out, fig, width = width, height = height, dpi = dpi, limitsize = FALSE)
  cat("wrote", out, "and", table_out, "\n")
  quit(status = 0)
}

# Long form, one value per measure. Per run, call's events then the VCF, each with a bar per size class (all, and
# call's threshold or more; one when benchmark ran at call's threshold).
share <- function(left, of) ifelse(of > 0, 100 * left / of, NA)
long <- rbind(transform(d, measure = "events", value = share(left, present), of = present),
              transform(d, measure = "bases", value = share(left_bp, present_bp), of = present_bp))
long <- long[long$measure %in% names(measures) & long$scope %in% names(scopes), ]
sizes_differ <- any(d$benchmark_min_bp < d$call_min_bp)
if (!sizes_differ) long <- long[long$size == "all", ]
group_at <- c(call = -0.22, VCF = 0.22)
size_at <- if (sizes_differ) c(all = -0.08, above = 0.08) else c(all = 0)
half <- if (sizes_differ) 0.07 else 0.1
long$bar <- paste(long$stage, long$size, sep = ".")
long$i <- match(long$run, labels)
long$chunk <- (long$i - 1) %/% per_row
long$x <- (long$i - 1) %% per_row + 1 + group_at[long$stage] + size_at[long$size]
long$part <- factor(long$part, part_levels)
long <- long[order(long$measure, long$scope, long$run, long$bar, long$part), ]
key <- paste(long$measure, long$scope, long$run, long$bar)
long$v0 <- ifelse(is.na(long$value), 0, long$value)
long$ymax <- ave(long$v0, key, FUN = cumsum)
long$ymin <- long$ymax - long$v0
chunks <- max(long$chunk) + 1
# Bases are written in one unit per panel, named in its title: kb once what was present reaches 10 kb.
unit_key <- function(m, s, ch) paste(m, s, ch)
bp_max <- tapply(long$of[long$measure == "bases"], unit_key("bases", long$scope, long$chunk)[long$measure == "bases"], max)
unit_of <- function(m, s, ch) {
  if (m != "bases") return("")
  v <- bp_max[unit_key(m, s, ch)]
  if (is.na(v) || v < 1e4) "bp" else if (v < 1e8) "kb" else "Mb"
}
unit_div <- c(bp = 1, kb = 1e3, Mb = 1e6)
panel_of <- function(m, s, ch) {
  u <- unit_of(m, s, ch)
  lab <- paste0(measures[m], if (nzchar(u)) paste0(" (", u, ")") else "", " · ", scopes[s])
  if (chunks > 1) lab <- paste0(lab, " · runs ", ch * per_row + 1, "-", min((ch + 1) * per_row, length(labels)))
  lab
}
grid <- expand.grid(s = names(scopes), ch = seq_len(chunks) - 1, m = names(measures), stringsAsFactors = FALSE)
panel_levels <- mapply(panel_of, grid$m, grid$s, grid$ch)
long$panel <- factor(mapply(panel_of, long$measure, long$scope, long$chunk), panel_levels)
long$n_left <- ifelse(long$measure == "events", long$left, long$left_bp)
one_t <- length(unique(d$call_min_bp)) == 1
size_name <- function(sz, t) ifelse(sz == "all", "all", if (one_t) paste0("≥", t) else "≥T")
long$size_label <- factor(size_name(long$size, long$call_min_bp), unique(size_name(c("all", "above"), d$call_min_bp[1])))
# Per bar: its top, and what is left of what was present, written above it.
bars <- unique(long[, c("panel", "run", "bar", "stage", "size_label", "x", "of", "measure", "scope", "chunk")])
k_bar <- paste(long$panel, long$run, long$bar)
k_bars <- paste(bars$panel, bars$run, bars$bar)
bars$n_left <- ave(long$n_left, k_bar, FUN = sum)[match(k_bars, k_bar)]
bars$top <- ave(long$ymax, k_bar, FUN = max)[match(k_bars, k_bar)]
# in the panel's unit: whole numbers from 100 up, one decimal below
in_unit <- function(v, u) {
  x <- v / unit_div[u]
  ifelse(x == 0, "0", ifelse(x >= 100 | u == "bp", format(round(x), big.mark = ",", trim = TRUE), sprintf("%.1f", x)))
}
bars$unit <- mapply(function(m, sc, ch) unit_of(m, sc, ch), bars$measure, bars$scope, bars$chunk)
bars$text <- ifelse(bars$measure == "events", paste0(bars$n_left, " / ", bars$of),
                    paste0(mapply(in_unit, bars$n_left, bars$unit), " / ", mapply(in_unit, bars$of, bars$unit)))
groups <- unique(long[, c("panel", "run", "i", "stage")])
groups$x <- (groups$i - 1) %% per_row + 1 + group_at[groups$stage]
names_at <- unique(long[, c("panel", "run", "i", "call_min_bp")])
names_at$x <- (names_at$i - 1) %% per_row + 1
names_at$text <- if (one_t) toupper(names_at$run) else paste0(toupper(names_at$run), " (T = ", names_at$call_min_bp, ")")
over <- unique(long[!is.na(long$value) & long$ymax > 100, "panel", drop = FALSE])
segs <- long[!is.na(long$part) & long$v0 > 0, ]
row_runs <- min(per_row, length(labels))
# every panel spans a full row of runs and at least 0-1%, so bars keep their width
frame <- expand.grid(panel = factor(panel_levels, panel_levels), x = c(0.55, row_runs + 0.45), y = c(0, 1))

# Bars stacked by cause; under them the size class (when there are two), call or VCF, and the run.
fig <- ggplot() +
  geom_blank(data = frame, aes(x = x, y = y)) +
  geom_rect(data = segs, aes(xmin = x - half, xmax = x + half, ymin = ymin, ymax = ymax, fill = part),
            colour = "white", linewidth = 0.25) +
  geom_text(data = bars, aes(x = x, y = top, label = text), angle = 90, hjust = -0.08, size = 2.3, colour = "grey25") +
  geom_text(data = groups, aes(x = x, y = -Inf, label = stage), vjust = if (sizes_differ) 3.2 else 1.6, size = 3) +
  geom_text(data = names_at, aes(x = x, y = -Inf, label = text), vjust = if (sizes_differ) 4.8 else 3.4, size = 3.2) +
  scale_fill_manual(values = palette, breaks = part_levels, drop = TRUE, name = NULL) +
  scale_x_continuous(breaks = NULL, expand = expansion(add = 0.05)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.55))) +
  facet_wrap(~panel, ncol = length(scopes), scales = "free", drop = FALSE) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = "% still left") +
  theme_bw(base_size = 10) +
  theme(legend.position = "top", legend.box = "vertical", legend.margin = margin(0, 0, 0, 0), panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        strip.background = element_rect(fill = "grey95"), panel.spacing.y = unit(if (sizes_differ) 0.75 else 0.6, "in"),
        plot.margin = margin(5.5, 5.5, if (sizes_differ) 42 else 30, 5.5)) +
  guides(fill = guide_legend(ncol = min(5, max(2, row_runs * length(scopes)))))
# with two size classes, each bar is named under it
if (sizes_differ) fig <- fig + geom_text(data = bars, aes(x = x, y = -Inf, label = size_label), vjust = 1.5, size = 2.3,
                                         colour = "grey30")
if (nrow(over) > 0) fig <- fig + geom_hline(data = over, aes(yintercept = 100), linetype = "dashed", colour = "grey50")
width <- as.numeric(get("--width", as.character(2 + length(scopes) * (1 + 1.1 * row_runs))))
height <- as.numeric(get("--height", as.character(1.5 + length(measures) * chunks * 2.8)))
ggsave(out, fig, width = width, height = height, dpi = dpi, limitsize = FALSE)
cat("wrote", out, "and", table_out, "\n")
