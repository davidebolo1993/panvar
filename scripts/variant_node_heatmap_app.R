#!/usr/bin/env Rscript
# Interactive viewer of call's records on the graph (Shiny + plotly), from a bundle made by
# build_variant_node_data.R. See --help for usage.
if (any(commandArgs(trailingOnly = TRUE) %in% c("-h", "--help"))) {
  cat(paste(c(
    "variant_node_heatmap_app.R - call's records on the graph, interactively (Shiny + plotly)",
    "",
    "Usage:",
    "  VN_RDS=<bundle.rds> Rscript variant_node_heatmap_app.R",
    "",
    "  VN_RDS   a bundle from build_variant_node_data.R",
    "",
    "Panels for the chosen site, or the whole locus (each can be hidden):",
    "  passes     a row per haplotype, the nodes in bins of equal bases, coloured by the most passes through",
    "             a node in the bin (white none, grey once, red twice or more)",
    "  records    a row per record, dark where its nodes are; zooming here or in passes zooms both",
    "  genotypes  a row per haplotype, a column per record, coloured by the event carried; on a CNV, by copies",
    "             gained or lost against the reference",
    "Click a record to show its carriers."),
    collapse = "\n"), "\n")
  quit(status = 0)
}
suppressWarnings(suppressMessages({library(shiny); library(plotly); library(data.table)}))

MAX_BINS <- 800L   # columns of the passes panel: bounds the data sent to the browser
GAP_BINS <- 6L     # blank columns between the regions shown
VAR_PX <- 18L
PASSES_MAX_PX <- 480L   # the passes panel scrolls past this height, so the records below it stay in view
MARGIN_L <- 170L        # same left margin on both panels: their columns line up
hap_px <- function(n) if (n > 100L) 10L else 16L   # haplotype rows: shorter when there are many

load_bundle <- function() {
  p <- Sys.getenv("VN_RDS", "")
  if (p == "" || !file.exists(p)) stop("set VN_RDS to a bundle from build_variant_node_data.R")
  b <- readRDS(p)
  for (k in c("nodes", "coverage", "variants", "genotypes", "bubbles")) b[[k]] <- as.data.table(b[[k]])
  b$coverage <- b$coverage[!is.na(count) & count > 0]
  setkey(b$nodes, order)
  # sample#hap, made unique where two contigs of one haplotype share it (rows are drawn by label)
  short <- sub("([^#]+#[^#]+)#.*", "\\1", b$haplotypes)
  dup <- duplicated(short) | duplicated(short, fromLast = TRUE)
  short[dup] <- paste0(short[dup], " (", ave(seq_along(short), short, FUN = seq_along)[dup], ")")
  b$short <- setNames(short, b$haplotypes)
  b
}

# The node-order ranges shown: the site +/- flank, or every site (and the locus ends) for "all"; ranges
# closer than `gap` nodes are merged.
regions_for <- function(b, bubble, flank = 20L, ends = 20L, gap = 20L) {
  N <- nrow(b$nodes)
  bs <- if (identical(bubble, "all")) b$bubbles else b$bubbles[bubble_id == bubble]
  r <- data.table(lo = pmax(1L, bs$lo - flank), hi = pmin(N, bs$hi + flank))
  if (identical(bubble, "all"))
    r <- rbind(r, data.table(lo = 1L, hi = min(ends, N)), data.table(lo = max(1L, N - ends + 1L), hi = N))
  setorder(r, lo)
  out <- r[1]
  for (i in seq_len(nrow(r))[-1]) {
    if (r$lo[i] <= out$hi[nrow(out)] + gap) out$hi[nrow(out)] <- max(out$hi[nrow(out)], r$hi[i])
    else out <- rbind(out, r[i])
  }
  out
}

# Columns: each region cut into bins of equal bases (at least one per region, MAX_BINS in all), with a few
# blank columns between regions. Returns the bins (label, region) and which bins each node covers.
build_bins <- function(b, reg) {
  nd <- b$nodes[order %in% unlist(Map(seq, reg$lo, reg$hi))]
  nd[, region := findInterval(order, reg$lo)]
  bp <- nd[, .(bp = sum(length)), by = region]
  bp[, nbins := pmax(1L, as.integer(floor(bp / sum(bp) * (MAX_BINS - GAP_BINS * (.N - 1)))))]
  nd[, `:=`(s = cumsum(length) - length, e = cumsum(length)), by = region]
  nd[bp, `:=`(w = i.bp / i.nbins, nbins = i.nbins), on = "region"]
  nd[, `:=`(b0 = pmin(nbins - 1L, as.integer(s %/% w)), b1 = pmin(nbins - 1L, as.integer((e - 1) %/% w)))]
  offs <- cumsum(c(0L, head(bp$nbins + GAP_BINS, -1)))
  nd[, base := offs[region]]
  # each bin a node reaches, with the node's bases inside it
  node_bin <- nd[, {l <- seq(b0, b1); .(bin = base + l, ov = pmin(e, (l + 1) * w) - pmax(s, l * w))},
                 by = .(node_id, order, gene, length)]
  total <- sum(bp$nbins) + GAP_BINS * (nrow(bp) - 1L)
  lab <- node_bin[, .(first = min(order), last = max(order), n = uniqueN(node_id), bp = round(sum(ov)),
                      genes = paste(setdiff(unique(unlist(strsplit(gene, ";", fixed = TRUE))), "."), collapse = ";")),
                  by = bin]
  labels <- sprintf("gap: %d", seq_len(total))   # unique, as the column labels must be
  first_id <- b$nodes$node_id[lab$first]; last_id <- b$nodes$node_id[lab$last]
  labels[lab$bin + 1L] <- sprintf("column: %d<br>nodes: %s<br>node count: %d<br>bases: %d%s", lab$bin + 1L,
                                  ifelse(first_id == last_id, first_id, paste0(first_id, "..", last_id)),
                                  lab$n, lab$bp, ifelse(lab$genes == "", "", paste0("<br>genes: ", gsub(";", ", ", lab$genes))))
  list(node_bin = node_bin, labels = labels, total = total, used = sort(unique(node_bin$bin)) + 1L)
}

# 0 to 8 passes: white none, grey once, red from twice, darker with more
covscale <- list(c(0, "#ffffff"), c(0.06, "#d9d9d9"), c(0.125, "#9e9e9e"), c(0.13, "#fc9272"),
                 c(0.25, "#ef3b2c"), c(1, "#67000d"))

passes_panel <- function(b, B, hp) {
  cov <- b$coverage[haplotype %in% hp][B$node_bin, on = "node_id", nomatch = NULL, allow.cartesian = TRUE]
  cov <- cov[, .(count = max(count)), by = .(haplotype, bin)]
  Z <- matrix(NA_real_, nrow = length(hp), ncol = B$total)
  Z[, B$used] <- 0
  Z[cbind(match(cov$haplotype, hp), cov$bin + 1L)] <- cov$count
  plot_ly(z = Z, x = B$labels, y = unname(b$short[hp]), type = "heatmap", zmin = 0, zmax = 8,
          colorscale = covscale, showscale = FALSE, source = "vn", ygap = 1,
          hovertemplate = "haplotype: %{y}<br>%{x}<br>passes: %{z}<extra></extra>") |>
    layout(yaxis = list(title = "", autorange = "reversed", type = "category", tickmode = "linear", dtick = 1,
                        tickfont = list(size = 9), automargin = FALSE),
           xaxis = list(title = "", showticklabels = FALSE, ticks = "", showgrid = FALSE, zeroline = FALSE, type = "category"),
           margin = list(l = MARGIN_L, r = 20, t = 10, b = 10))
}

shown_variants <- function(b, bubble) {
  if (is.null(bubble) || bubble == "") return(character(0))
  if (identical(bubble, "all")) b$variants$variant_id else b$variants[bubble_id == bubble, variant_id]
}

describe_variant <- function(v) {
  sprintf("record: %s<br>type: %s%s<br>position: %s<br>carriers: %s of %s%s%s", v$variant_id, v$svtype,
          ifelse(v$svlen == ".", "", sprintf("<br>size: %s bp", v$svlen)), v$pos, v$ac, v$an,
          ifelse(v$cn_method == ".", "", paste0("<br>copies counted: ", tolower(v$cn_method))),
          ifelse(v$genes == ".", "", paste0("<br>genes: ", gsub(",", ", ", v$genes))))
}

records_panel <- function(b, B, vids) {
  vm <- b$variants[match(vids, variant_id)]
  hits <- rbindlist(lapply(seq_along(vids), function(i) {
    ns <- strsplit(vm$nodes[i], ",", fixed = TRUE)[[1]]
    B$node_bin[node_id %in% ns, .(row = i, bin)]
  }))
  Z <- matrix(NA_real_, nrow = length(vids), ncol = B$total)
  Z[, B$used] <- 0
  # the description along the whole row; past 100 records only on its own columns, to keep the page small
  desc <- describe_variant(vm)
  TX <- if (length(vids) <= 100L) matrix(desc, nrow = length(vids), ncol = B$total)
        else matrix("", nrow = length(vids), ncol = B$total)
  if (nrow(hits)) {
    Z[cbind(hits$row, hits$bin + 1L)] <- 1
    TX[cbind(hits$row, hits$bin + 1L)] <- desc[hits$row]
  }
  plot_ly(z = Z, x = B$labels, y = vids, type = "heatmap", zmin = 0, zmax = 1, showscale = FALSE,
          colorscale = list(c(0, "#ffffff"), c(1, "#8b0000")), source = "rec", ygap = 2,
          text = TX, hoverinfo = "text") |>
    layout(yaxis = list(title = "", autorange = "reversed", type = "category", tickmode = "linear", dtick = 1,
                        tickfont = list(size = 11), automargin = FALSE),
           xaxis = list(title = "", showticklabels = FALSE, ticks = "", showgrid = FALSE, zeroline = FALSE, type = "category"),
           margin = list(l = MARGIN_L, r = 20, t = 10, b = 10))
}

# Genotype classes, one colour each: no allele, reference-like, then the event carried.
GT_CLASSES <- c("no allele", "reference", "DEL", "INS", "INV", "SUB", "CN gain", "CN loss", "CN other")
GT_COLOURS <- c("#ffffff", "#ececec", "#e41a1c", "#4daf4a", "#ff7f00", "#a65628", "#377eb8", "#984ea3", "#9ecae1")

genotype_panel <- function(b, hp, vids) {
  vm <- b$variants[match(vids, variant_id)]
  g <- b$genotypes[variant_id %in% vids & haplotype %in% hp]
  g[vm, `:=`(svtype = i.svtype, ref_cn = i.ref_cn), on = "variant_id"]
  g[, `:=`(cnn = suppressWarnings(as.numeric(cn)), refn = suppressWarnings(as.numeric(ref_cn)))]
  g[, cls := fifelse(gt == "." | gt == "", "no allele",
             fifelse(svtype == "CNV" & !is.na(cnn) & !is.na(refn),
                     fifelse(cnn > refn, "CN gain", fifelse(cnn < refn, "CN loss",
                                                            fifelse(startsWith(gt, "1"), "CN other", "reference"))),
                     fifelse(startsWith(gt, "1"), svtype, "reference")))]
  g[!cls %in% GT_CLASSES, cls := "CN other"]
  Z <- matrix(0, nrow = length(hp), ncol = length(vids))
  TX <- matrix("", nrow = length(hp), ncol = length(vids))
  r <- match(g$haplotype, hp); k <- match(g$variant_id, vids)
  Z[cbind(r, k)] <- match(g$cls, GT_CLASSES) - 1
  TX[cbind(r, k)] <- sprintf("haplotype: %s<br>record: %s<br>GT: %s%s", b$short[g$haplotype], g$variant_id, g$gt,
                             ifelse(g$cn == ".", "", paste0("<br>CN: ", g$cn)))
  n <- length(GT_CLASSES)
  scale <- do.call(c, lapply(seq_len(n), function(i)
    list(list((i - 1) / n, GT_COLOURS[i]), list(i / n, GT_COLOURS[i]))))
  plot_ly(z = Z, x = vids, y = unname(b$short[hp]), text = TX, type = "heatmap", zmin = -0.5, zmax = n - 0.5,
          colorscale = scale, source = "gt", xgap = 1, ygap = 1, hoverinfo = "text",
          colorbar = list(tickvals = seq_len(n) - 1, ticktext = GT_CLASSES, title = "", lenmode = "pixels",
                          len = 260, yanchor = "top", y = 1)) |>
    layout(yaxis = list(title = "", autorange = "reversed", type = "category", tickmode = "linear", dtick = 1,
                        tickfont = list(size = 9)),
           xaxis = list(title = "records", tickangle = -45, type = "category"))
}

# Passes and records are two plots on the same bins; passes scrolls on its own so the records stay in view.
passes_height <- function(nh) as.integer(20L + nh * hap_px(nh))
records_height <- function(nv) as.integer(20L + nv * VAR_PX)

# The opening selection: the reference, the first carrier of each record and, on a CNV, the first haplotype of
# each copy number.
rep_haps <- function(b) {
  sel <- b$reference
  for (v in unique(b$variants$variant_id)) {
    if (b$variants[variant_id == v, svtype][1] == "CNV") {
      g <- b$genotypes[variant_id == v & cn != "."]
      for (cv in unique(g$cn)) sel <- c(sel, g[cn == cv, haplotype][1])
    } else {
      sel <- c(sel, b$genotypes[variant_id == v & startsWith(gt, "1"), haplotype][1])
    }
  }
  unique(sel[!is.na(sel)])
}

if (!isTRUE(getOption("vn.testing"))) {
  B0 <- load_bundle()
  choices <- setNames(B0$haplotypes, B0$short[B0$haplotypes])
  ui <- fluidPage(
    tags$style(HTML(".selectize-control{max-width:230px;} .selectize-input{max-height:260px; overflow-y:auto;}")),
    titlePanel("panvar calls", windowTitle = "panvar calls"),
    sidebarLayout(
      sidebarPanel(width = 3,
        selectInput("bubble", "Site", c("all (whole locus)" = "all", B0$bubbles$bubble_id), "all"),
        selectizeInput("haps", "Haplotypes (clear = all)", choices = NULL, multiple = TRUE),
        checkboxGroupInput("panels", "Panels", c("haplotype passes" = "passes", "records" = "records",
                                                 "genotypes" = "genotypes"), c("passes", "records", "genotypes")),
        helpText("Opens on the reference and one carrier per record. Clear the box to show every haplotype.",
                 "Click a record to show its carriers.")),
      mainPanel(width = 9, uiOutput("passes_ui"), uiOutput("rec_ui"), uiOutput("gt_ui"))))
  server <- function(input, output, session) {
    updateSelectizeInput(session, "haps", choices = choices, selected = rep_haps(B0), server = TRUE)
    shows <- function(panel) panel %in% input$panels
    site <- reactive(if (is.null(input$bubble) || input$bubble == "") NULL else input$bubble)
    hp <- reactive(if (length(input$haps)) input$haps else B0$haplotypes)
    vids <- reactive(shown_variants(B0, site()))
    bins <- reactive({ req(site()); build_bins(B0, regions_for(B0, site())) })
    output$passes <- renderPlotly({ req(shows("passes")); event_register(passes_panel(B0, bins(), hp()), "plotly_relayout") })
    output$passes_ui <- renderUI({
      if (!shows("passes") || is.null(site())) return(NULL)
      div(style = paste0("max-height:", PASSES_MAX_PX, "px; overflow-y:auto;"),
          plotlyOutput("passes", height = paste0(passes_height(length(hp())), "px")))
    })
    output$rec <- renderPlotly({ req(shows("records"), length(vids()) > 0)
      event_register(records_panel(B0, bins(), vids()), "plotly_relayout") })
    output$rec_ui <- renderUI({
      if (!shows("records") || !length(vids())) return(NULL)
      plotlyOutput("rec", height = paste0(records_height(length(vids())), "px"))
    })
    # Zoom along x in one panel, the other follows.
    follow <- function(from, to) observeEvent(event_data("plotly_relayout", source = from), {
      e <- event_data("plotly_relayout", source = from)
      x <- if (!is.null(e[["xaxis.range[0]"]])) list(xaxis = list(range = c(e[["xaxis.range[0]"]], e[["xaxis.range[1]"]])))
           else if (isTRUE(e[["xaxis.autorange"]])) list(xaxis = list(autorange = TRUE))
      if (!is.null(x)) plotlyProxy(to, session) |> plotlyProxyInvoke("relayout", x)
    })
    follow("vn", "rec"); follow("rec", "passes")
    output$gt <- renderPlotly({
      if (!length(vids()) || !shows("genotypes")) return(NULL)
      genotype_panel(B0, hp(), vids())
    })
    output$gt_ui <- renderUI({
      if (!shows("genotypes")) return(NULL)
      nh <- length(hp())
      plotlyOutput("gt", height = paste0(max(420L, 160L + nh * hap_px(nh)), "px"))   # room for the legend
    })
    # A click on a record selects its carriers and the reference; on a CNV every haplotype with a copy number,
    # as GT 1 marks only the larger changes there.
    show_carriers <- function(v) {
      if (!is.null(v) && v %in% B0$variants$variant_id) {
        cnv <- B0$variants[variant_id == v, svtype] == "CNV"
        car <- if (cnv) B0$genotypes[variant_id == v & cn != ".", haplotype]
               else B0$genotypes[variant_id == v & startsWith(gt, "1"), haplotype]
        updateSelectizeInput(session, "haps", selected = unique(c(B0$reference, car)))
      }
    }
    observeEvent(event_data("plotly_click", source = "rec"), show_carriers(event_data("plotly_click", source = "rec")$y))
    observeEvent(event_data("plotly_click", source = "gt"), show_carriers(event_data("plotly_click", source = "gt")$x))
  }
  shinyApp(ui, server)
}
