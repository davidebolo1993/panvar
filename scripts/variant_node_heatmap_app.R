#!/usr/bin/env Rscript
# Interactive node-coverage + variant-track viewer (Shiny + plotly). Bundle from build_variant_node_data.R:
#   VN_RDS=results/real_data/gstm1/call/variant_nodes.rds Rscript scripts/variant_node_heatmap_app.R
# Top: one row per selected haplotype, X = variant regions (each bubble + flanks; long invariant stretches
# collapsed to a fixed gap that shrinks when bubbles are close), node width proportional to length. Colour =
# per-walk coverage: white = not traversed, grey = x1, red = x2+ (so a DUP reads as a copy-number gradient).
# Bottom (toggle): one row per VCF variant_id, red at its union nodes. Hover shows node/gene/coverage and,
# on a variant node, that haplotype's GT (+ CN/CNBP for DUPs).
if (any(commandArgs(trailingOnly = TRUE) %in% c("-h", "--help"))) {
  cat(paste(c(
    "variant_node_heatmap_app.R - interactive node-coverage + variant-track viewer (Shiny + plotly).",
    "",
    "Usage:",
    "  VN_RDS=<bundle.rds> Rscript variant_node_heatmap_app.R",
    "",
    "  VN_RDS   path to a bundle from build_variant_node_data.R (required, passed as an env var)",
    "",
    "Top panel: one row per haplotype, X = variant regions, node colour = per-walk coverage",
    "  (white none, grey x1, red x2+). Bottom panel (toggle): one row per VCF variant at its nodes."),
    collapse = "\n"), "\n")
  quit(status = 0)
}
suppressWarnings(suppressMessages({library(shiny); library(plotly); library(data.table)}))

load_bundle <- function() {
  p <- Sys.getenv("VN_RDS", "")
  if (p == "" || !file.exists(p)) stop("set VN_RDS to a bundle from build_variant_node_data.R")
  b <- readRDS(p)
  for (k in c("nodes","coverage","variants","genotypes","bubbles")) b[[k]] <- as.data.table(b[[k]])
  b$coverage <- b$coverage[!is.na(count) & count > 0]
  setkey(b$nodes, order); setkey(b$coverage, haplotype, node_id)
  b$short <- setNames(sub("([^#]+#[^#]+)#.*", "\\1", b$haplotypes), b$haplotypes)  # sample#hap label
  b
}

# Regions (node-order ranges) to display: each bubble +/- flank (+ locus ends when showing all), merged when
# within `gap` of each other. Returns a data.table(lo, hi) sorted, non-overlapping.
regions_for <- function(b, bubble, flank = 20L, ends = 20L, gap = 20L) {
  N <- nrow(b$nodes)
  bs <- if (identical(bubble, "all")) b$bubbles else b$bubbles[bubble_id == bubble]
  r <- data.table(lo = pmax(1L, bs$lo - flank), hi = pmin(N, bs$hi + flank))
  if (identical(bubble, "all")) r <- rbind(r, data.table(lo = 1L, hi = min(ends, N)), data.table(lo = max(1L, N - ends + 1L), hi = N))
  setorder(r, lo)
  merged <- r[1];
  for (i in seq_len(nrow(r))[-1]) {
    if (r$lo[i] <= merged$hi[nrow(merged)] + gap) merged$hi[nrow(merged)] <- max(merged$hi[nrow(merged)], r$hi[i])
    else merged <- rbind(merged, r[i])
  }
  merged
}

# Column layout over the displayed regions: each node -> a run of columns ~ proportional to bp, a fixed gap
# strip between regions. Returns nd (displayed nodes, ordered) and col->node index (NA on gap columns).
build_display <- function(b, reg, target = 1600L, gapcols = 24L) {
  ord <- unlist(Map(function(lo, hi) lo:hi, reg$lo, reg$hi))
  nd <- b$nodes[order %in% ord]; setorder(nd, order)
  tot <- sum(nd$length)
  nd[, w := pmax(1L, round(length / tot * max(target - nrow(nd), nrow(nd))))]
  # region boundary after which to insert a gap
  regend <- cumsum(vapply(seq_len(nrow(reg)), function(i) sum(nd$order >= reg$lo[i] & nd$order <= reg$hi[i]), integer(1)))
  coln <- integer(0); k <- 0L
  for (i in seq_len(nrow(nd))) {
    coln <- c(coln, rep.int(i, nd$w[i]))
    k <- k + 1L
    if (k %in% regend && k < nrow(nd)) coln <- c(coln, rep.int(NA_integer_, gapcols))
  }
  list(nd = nd, coln = coln)
}

covscale <- list(c(0, "#ffffff"), c(0.0001, "#e8e8e8"), c(0.18, "#9e9e9e"),
                 c(0.36, "#fcae91"), c(0.6, "#fb6a4a"), c(1, "#67000d"))

render_view <- function(b, bubble = "all", haps = NULL, show_variants = TRUE) {
  if (is.null(bubble) || bubble == "")            # deferred: nothing rendered until a bubble is picked
    return(plot_ly() |> layout(xaxis = list(visible = FALSE), yaxis = list(visible = FALSE),
             annotations = list(text = "Pick a bubble to render  (or 'all' for the whole locus — slower).",
                                showarrow = FALSE, font = list(size = 16, color = "#666"))))
  reg <- regions_for(b, bubble)
  D <- build_display(b, reg); nd <- D$nd; coln <- D$coln; keep <- !is.na(coln)
  hp <- if (length(haps)) haps else b$haplotypes
  if (nrow(nd) == 0 || length(hp) == 0) return(plotly_empty())

  # coverage matrix (hap x node) over displayed nodes
  cov <- b$coverage[node_id %in% nd$node_id & haplotype %in% hp]
  Cm <- matrix(0L, nrow = length(hp), ncol = nrow(nd), dimnames = list(NULL, nd$node_id))
  if (nrow(cov)) Cm[cbind(match(cov$haplotype, hp), match(cov$node_id, nd$node_id))] <- cov$count
  Zc <- matrix(NA_real_, nrow = length(hp), ncol = length(coln)); Zc[, keep] <- Cm[, coln[keep], drop = FALSE]

  # every variant a node belongs to (a node can be in the DUP union and an INS/DEL), for the genotype hover
  rowof <- setNames(seq_len(nrow(nd)), nd$node_id)
  n2vs <- vector("list", nrow(nd))
  for (i in seq_len(nrow(b$variants))) {
    rr <- rowof[strsplit(b$variants$nodes[i], ",", fixed = TRUE)[[1]]]; rr <- rr[!is.na(rr)]
    for (r in rr) n2vs[[r]] <- c(n2vs[[r]], b$variants$variant_id[i])
  }
  # fast genotype suffix: "variant_id|haplotype" -> " | id GT=.. [CN=.. CNBP=..]"
  g <- b$genotypes
  gtmap <- setNames(ifelse(g$cn != ".",
                    sprintf(" | %s GT=%s CN=%s CNBP=%s", g$variant_id, g$gt, g$cn, g$cnbp),
                    sprintf(" | %s GT=%s", g$variant_id, g$gt)),
                    paste(g$variant_id, g$haplotype, sep = "|"))
  base <- sprintf("node %s | gene %s | %d bp", nd$node_id, nd$gene, nd$length)
  Cell <- matrix(base[rep(seq_len(nrow(nd)), each = length(hp))], nrow = length(hp))  # broadcast base per node
  for (j in which(lengths(n2vs) > 0)) {                                # variant nodes: append every variant's genotype
    suf <- rep("", length(hp))
    for (v in n2vs[[j]]) { s <- gtmap[paste(v, hp, sep = "|")]; s[is.na(s)] <- ""; suf <- paste0(suf, s) }
    Cell[, j] <- paste0(base[j], suf)
  }
  Cell <- matrix(paste0(b$short[hp], "<br>", Cell, " | cov ", Cm), nrow = length(hp))  # bake haplotype + coverage
  TXc <- matrix("", nrow = length(hp), ncol = length(coln)); TXc[, keep] <- Cell[, coln[keep], drop = FALSE]

  cov_hm <- plot_ly(z = Zc, y = b$short[hp], text = TXc, hoverinfo = "text", type = "heatmap", zmin = 0, zmax = 8,
                    ygap = 2, colorscale = covscale, showscale = FALSE, source = "vn") |>
    layout(yaxis = list(title = "", showticklabels = length(hp) <= 80, autorange = "reversed"),
           xaxis = list(showticklabels = FALSE, zeroline = FALSE))

  if (!show_variants) return(cov_hm |> layout(xaxis = list(title = "variant regions (width ∝ node length)")))

  # variant tracks: one row per variant. When a bubble is selected, show only that bubble's variants;
  # in the "all" view, show the variants that fall in the displayed node regions.
  vids <- shown_variants(b, bubble)
  if (!length(vids)) return(cov_hm |> layout(xaxis = list(title = "variant regions (width ∝ node length)")))
  Vm <- matrix(0L, nrow = length(vids), ncol = nrow(nd), dimnames = list(vids, nd$node_id))
  for (i in seq_along(vids)) {
    ns <- strsplit(b$variants[variant_id == vids[i], nodes], ",", fixed = TRUE)[[1]]
    Vm[i, colnames(Vm) %in% ns] <- 1L
  }
  Zv <- matrix(NA_real_, nrow = length(vids), ncol = length(coln)); Zv[, keep] <- Vm[, coln[keep], drop = FALSE]
  vmeta <- b$variants[match(vids, variant_id)]
  svl <- if ("svlen" %in% names(vmeta)) vmeta$svlen else rep(".", nrow(vmeta))   # older bundles lack svlen
  TXv <- matrix(sprintf("%s | %s | SVLEN %s | POS %s | gene %s", vmeta$variant_id, vmeta$svtype, svl, vmeta$pos, vmeta$gene),
                nrow = length(vids), ncol = length(coln))
  var_hm <- plot_ly(z = Zv, y = vids, text = TXv, hoverinfo = "text", type = "heatmap", zmin = 0, zmax = 1,
                    ygap = 2, colorscale = list(c(0, "#ffffff"), c(1, "#8b0000")), showscale = FALSE, source = "vn") |>
    layout(yaxis = list(title = "", showticklabels = TRUE, autorange = "reversed"),
           xaxis = list(title = "variant regions (width ∝ node length)", showticklabels = FALSE, zeroline = FALSE))

  hf <- (length(hp) * HAP_PX) / (length(hp) * HAP_PX + length(vids) * VAR_PX)  # split by each panel's row count
  subplot(cov_hm, var_hm, nrows = 2, shareX = TRUE, heights = c(hf, 1 - hf), titleY = TRUE)
}

# Representative default selection: the calling reference + one carrier per variant (one per distinct CN
# for DUPs), so the opening view is small but shows every called allele; the user adds/removes from there.
rep_haps <- function(b) {
  sel <- b$reference
  for (v in unique(b$variants$variant_id)) {
    g <- b$genotypes[variant_id == v & startsWith(gt, "1")]
    if (!nrow(g)) next
    if (any(g$cn != ".")) for (cv in unique(g$cn)) sel <- c(sel, g[cn == cv, haplotype][1])
    else sel <- c(sel, g$haplotype[1])
  }
  unique(sel[!is.na(sel)])
}

# Per-row pixel budget so haplotype cells stay short while every variant row keeps a readable label; the plot
# height then follows the row count (fixed heights thinned the 16 variant labels down to ~6 in the whole-locus view).
HAP_PX <- 22L; VAR_PX <- 18L
shown_variants <- function(b, bubble) {
  if (is.null(bubble) || bubble == "") return(character(0))
  if (identical(bubble, "all")) b$variants$variant_id else b$variants[bubble_id == bubble, variant_id]
}
view_height <- function(b, bubble, haps, show_vars) {
  if (is.null(bubble) || bubble == "") return(320L)
  nh <- if (length(haps)) length(haps) else length(b$haplotypes)
  nv <- if (isTRUE(show_vars)) length(shown_variants(b, bubble)) else 0L
  as.integer(90L + nh * HAP_PX + nv * VAR_PX)
}

if (!isTRUE(getOption("vn.testing"))) {
  B <- load_bundle()
  choices <- setNames(B$haplotypes, B$short[B$haplotypes])   # display short sample#hap, value = full name
  ui <- fluidPage(
    tags$style(HTML("#hm{min-width:1100px;} .selectize-control{max-width:230px;}
                     .selectize-input{max-height:260px; overflow-y:auto;}")),
    titlePanel("panvar — node coverage & variants per haplotype"),
    sidebarLayout(
      sidebarPanel(width = 3,
        selectInput("bubble", "Bubble", c("all (whole locus)" = "all", B$bubbles$bubble_id), "all"),
        selectizeInput("haps", "Haplotypes (clear = all)", choices = NULL, multiple = TRUE),
        checkboxInput("vars", "show VCF calls", TRUE),
        helpText("Opens on representative haplotypes (reference + one per allele/CN). Clear the box to show",
                 "all. White = not traversed, grey = traversed x1, red = x2+ (copy number).",
                 "Hover a node for coverage/genotype/gene information.")),
      mainPanel(width = 9, uiOutput("plot_ui"))))
  server <- function(input, output, session) {
    updateSelectizeInput(session, "haps", choices = choices, selected = rep_haps(B), server = TRUE)
    output$hm <- renderPlotly(render_view(B, input$bubble, input$haps, input$vars))
    output$plot_ui <- renderUI(plotlyOutput("hm",
      height = paste0(view_height(B, input$bubble, input$haps, input$vars), "px")))
    # changing the bubble resets the haplotypes back to the representative subset (drops a prior variant-click selection)
    observeEvent(input$bubble, updateSelectizeInput(session, "haps", selected = rep_haps(B)), ignoreInit = TRUE)
    # click a variant track row -> show every haplotype carrying that variant (+ the reference)
    observeEvent(event_data("plotly_click", source = "vn"), {
      cd <- event_data("plotly_click", source = "vn")
      v <- cd$y
      if (!is.null(v) && v %in% B$variants$variant_id) {
        car <- B$genotypes[variant_id == v & startsWith(gt, "1"), haplotype]
        updateSelectizeInput(session, "haps", selected = unique(c(B$reference, car)))
      }
    })
  }
  shinyApp(ui, server)
}
