#!/usr/bin/env Rscript
# Assemble the .rds bundle (nodes, coverage, variants, genotypes, bubbles) for variant_node_heatmap_app.R.
# See --help for usage.
suppressWarnings(suppressMessages(library(data.table)))
args <- commandArgs(trailingOnly = TRUE)
usage <- function(status = 0) {
  cat(paste(c(
    "build_variant_node_data.R - assemble the .rds bundle for variant_node_heatmap_app.R",
    "  (nodes, per-haplotype coverage, variants, genotypes, sites).",
    "",
    "Usage:",
    "  Rscript build_variant_node_data.R --gfa <graph.gfa> --prefix <call prefix> --bubbles <bubbles.csv> \\",
    "      [--out <bundle.rds>]",
    "",
    "Required:",
    "  --gfa <path>        the graph call was run on",
    "  --prefix <prefix>   call's output prefix: reads <prefix>.region.vcf and <prefix>.variant_nodes.tsv, and",
    "                      <prefix>.node_genes.tsv when it is there",
    "  --bubbles <path>    the site table call read",
    "",
    "Optional:",
    "  --out <path>        output .rds (default <prefix>.viewer.rds)",
    "  -h, --help          show this help"), collapse = "\n"), "\n")
  quit(status = status)
}
if (length(args) == 0 || any(args %in% c("-h", "--help"))) usage(0)
get <- function(f, d = NULL) {
  i <- match(f, args)
  if (is.na(i)) return(d)
  if (i == length(args)) stop("missing value after ", f)
  args[i + 1]
}
gfa <- get("--gfa"); prefix <- get("--prefix"); bcsv <- get("--bubbles")
if (is.null(gfa) || is.null(prefix) || is.null(bcsv)) usage(1)
outp <- get("--out", paste0(prefix, ".viewer.rds"))
vcfp <- paste0(prefix, ".region.vcf"); vnp <- paste0(prefix, ".variant_nodes.tsv")
ngp <- paste0(prefix, ".node_genes.tsv")
for (p in c(gfa, vcfp, vnp, bcsv)) if (!file.exists(p)) stop("not found: ", p)

# --- graph: node lengths, and each path's nodes in order ---
L <- readLines(gfa)
sf <- tstrsplit(sub("^S\t", "", L[startsWith(L, "S\t")]), "\t", fixed = TRUE)
nodes <- data.table(node_id = sf[[1]], length = nchar(sf[[2]]))
lenmap <- setNames(nodes$length, nodes$node_id)
# P lines, and W lines named as panvar names them (sample#hap#seqid:start+1-end).
P <- L[startsWith(L, "P\t")]
W <- L[startsWith(L, "W\t")]
paths_ids <- list(); haplotypes <- character(0)
if (length(P)) {
  pf <- tstrsplit(P, "\t", fixed = TRUE)
  paths_ids <- lapply(strsplit(pf[[3]], ",", fixed = TRUE), function(s) sub("[+-]$", "", s))
  haplotypes <- pf[[2]]
}
if (length(W)) {
  wf <- tstrsplit(W, "\t", fixed = TRUE)
  numeric_start <- grepl("^[0-9]+$", wf[[5]])
  wnames <- ifelse(numeric_start,
                   sprintf("%s#%s#%s:%s-%s", wf[[2]], wf[[3]], wf[[4]],
                           format(suppressWarnings(as.numeric(wf[[5]])) + 1, scientific = FALSE, trim = TRUE),
                           ifelse(grepl("^[0-9]+$", wf[[6]]), wf[[6]], "*")),
                   sprintf("%s#%s#%s", wf[[2]], wf[[3]], wf[[4]]))
  paths_ids <- c(paths_ids, lapply(regmatches(wf[[7]], gregexpr("[^<>]+", wf[[7]])), identity))
  haplotypes <- c(haplotypes, wnames)
}
if (!length(haplotypes)) stop(gfa, " has no P or W paths")
names(paths_ids) <- haplotypes
rm(L, P, W, sf)
coverage <- rbindlist(lapply(haplotypes, function(h) {
  tb <- table(paths_ids[[h]])
  data.table(haplotype = h, node_id = names(tb), count = as.integer(tb))
}))

# --- position of each node along the reference; a node off it takes the median position, over the walks
#     visiting it, of the last reference node before it ---
vl <- readLines(vcfp)
refname <- sub(".*##reference=", "", grep("^##reference=", vl, value = TRUE)[1])
refkey <- if (refname %in% haplotypes) refname else grep(refname, haplotypes, fixed = TRUE, value = TRUE)[1]
if (is.na(refkey)) stop("the VCF's reference ", refname, " is not a path of ", gfa)
refpos <- setNames(rep(NA_real_, nrow(nodes)), nodes$node_id)
ref_ids <- paths_ids[[refkey]]
first <- !duplicated(ref_ids)
refpos[ref_ids[first]] <- (cumsum(lenmap[ref_ids]) - lenmap[ref_ids])[first]
proj <- rbindlist(lapply(haplotypes, function(h) {
  ids <- paths_ids[[h]]
  data.table(node_id = ids, pos = nafill(unname(refpos[ids]), "locf"))[is.na(refpos[node_id]) & !is.na(pos)]
}))
gpos <- refpos
if (nrow(proj)) { pm <- proj[, .(pos = median(pos)), by = node_id]; gpos[pm$node_id] <- pm$pos }
gpos[is.na(gpos)] <- max(gpos, na.rm = TRUE) + 1   # placed nowhere: at the end
nodes[, gpos := gpos[node_id]]
setorder(nodes, gpos, node_id)
nodes[, order := .I]
order_of <- setNames(nodes$order, nodes$node_id)

# --- genes per node ---
nodes[, gene := "."]
if (file.exists(ngp)) {
  ng <- fread(ngp, sep = "\t", header = TRUE, colClasses = "character")
  nodes[ng, gene := i.genes, on = "node_id"]
}

# --- the VCF: records with the INFO fields the viewer shows, and each haplotype's genotype ---
body <- vl[!startsWith(vl, "#")]
samples <- strsplit(vl[startsWith(vl, "#CHROM")][1], "\t", fixed = TRUE)[[1]][-(1:9)]
info_of <- function(info, key) {
  m <- regmatches(info, regexpr(paste0("(^|;)", key, "=[^;]*"), info))
  if (length(m)) sub(paste0(".*", key, "="), "", m) else "."
}
recs <- lapply(body, function(bl) {
  f <- strsplit(bl, "\t", fixed = TRUE)[[1]]
  info <- f[8]
  rec <- data.table(variant_id = f[3], pos = as.integer(f[2]), end = info_of(info, "END"),
                    svtype = info_of(info, "SVTYPE"), svlen = info_of(info, "SVLEN"),
                    ac = info_of(info, "AC"), an = info_of(info, "AN"), af = info_of(info, "AF"),
                    cn_method = info_of(info, "CN_METHOD"), ref_cn = info_of(info, "REF_CN"),
                    genes = utils::URLdecode(info_of(info, "GENES")),
                    imprecise = grepl("(^|;)IMPRECISE(;|$)", info))
  keys <- strsplit(f[9], ":", fixed = TRUE)[[1]]
  vals <- strsplit(f[-(1:9)], ":", fixed = TRUE)
  pick <- function(k) {
    i <- match(k, keys)
    if (is.na(i)) rep(".", length(vals)) else vapply(vals, function(v) if (length(v) >= i) v[i] else ".", "")
  }
  gt <- data.table(variant_id = f[3], haplotype = samples, gt = pick("GT"), cn = pick("CN"),
                   cnbp = pick("CNBP"), cnresid = pick("CNRESID"))
  list(rec, gt)
})
records <- rbindlist(lapply(recs, `[[`, 1))
genotypes <- rbindlist(lapply(recs, `[[`, 2))
if (!nrow(records)) {   # a run with no records: empty panels
  records <- data.table(variant_id = character(0), pos = integer(0), end = character(0), svtype = character(0),
                        svlen = character(0), ac = character(0), an = character(0), af = character(0),
                        cn_method = character(0), ref_cn = character(0), genes = character(0), imprecise = logical(0))
  genotypes <- data.table(variant_id = character(0), haplotype = character(0), gt = character(0), cn = character(0),
                          cnbp = character(0), cnresid = character(0))
}

# A --diploid VCF has a column per sample: split back into its haplotype paths in call's slot order (numbered
# haplotypes ascending, then each unnumbered path). GT and HCN give a value per slot; CNBP and CNRESID are not
# split. The reference, not a column, reads as the reference allele.
if (any(startsWith(vl, "##panvar_columns=diploid"))) {
  others <- setdiff(haplotypes, refkey)
  sample_of <- sub("#.*", "", others)
  hapf <- ifelse(grepl("#", others), sub("^[^#]*#([^#]*).*", "\\1", others), "0")
  slots <- lapply(split(seq_along(others), sample_of), function(i) {
    num <- grepl("^[0-9]+$", hapf[i])
    numbered <- split(others[i][num], as.numeric(hapf[i][num]))
    c(unname(numbered[order(as.numeric(names(numbered)))]), as.list(others[i][!num]))
  })
  rows <- list()
  keys <- lapply(body, function(bl) strsplit(strsplit(bl, "\t", fixed = TRUE)[[1]][9], ":", fixed = TRUE)[[1]])
  for (r in seq_along(body)) {
    f <- strsplit(body[r], "\t", fixed = TRUE)[[1]]
    vals <- strsplit(f[-(1:9)], ":", fixed = TRUE)
    gi <- match("GT", keys[[r]]); hi <- match("HCN", keys[[r]])
    for (k in seq_along(samples)) {
      sl <- slots[[samples[k]]]
      if (is.null(sl)) next
      gts <- strsplit(vals[[k]][gi], "[|/]")[[1]]
      hcn <- if (is.na(hi) || length(vals[[k]]) < hi) rep(".", length(gts)) else strsplit(vals[[k]][hi], "[|/]")[[1]]
      for (s in seq_along(sl))
        rows[[length(rows) + 1]] <- data.table(variant_id = f[3], haplotype = sl[[s]],
                                               gt = if (s <= length(gts)) gts[s] else ".",
                                               cn = if (s <= length(hcn)) hcn[s] else ".", cnbp = ".", cnresid = ".")
    }
    rows[[length(rows) + 1]] <- data.table(variant_id = f[3], haplotype = refkey, gt = "0",
                                           cn = records[variant_id == f[3], fifelse(svtype == "CNV", ref_cn, ".")],
                                           cnbp = ".", cnresid = ".")
  }
  genotypes <- if (length(rows)) rbindlist(rows) else genotypes[0]
}


# --- sites, as node-order spans ---
bub <- fread(bcsv, sep = ",", header = TRUE, colClasses = "character")
if ("kind" %in% names(bub)) bub <- bub[kind != "inversion"]
bubbles <- rbindlist(lapply(seq_len(nrow(bub)), function(i) {
  ids <- c(bub$source[i], bub$sink[i], strsplit(bub$inside_nodes[i], ";", fixed = TRUE)[[1]])
  o <- order_of[ids]; o <- o[!is.na(o)]
  if (!length(o)) return(NULL)
  data.table(bubble_id = bub$bubble_id[i], lo = min(o), hi = max(o))
}))

vn <- fread(vnp, sep = "\t", header = TRUE, colClasses = "character")
variants <- merge(vn[, .(variant_id, bubble_id, nodes = node_ids)], records, by = "variant_id", sort = FALSE)
setorder(variants, pos)
setkey(nodes, node_id)

saveRDS(list(nodes = nodes, coverage = coverage, variants = variants, genotypes = genotypes,
             bubbles = bubbles, haplotypes = haplotypes, reference = refkey), outp)
cat(sprintf("wrote %s: %d nodes, %d variants, %d sites, %d haplotypes (reference %s)\n",
            outp, nrow(nodes), nrow(variants), nrow(bubbles), length(haplotypes), refkey))
