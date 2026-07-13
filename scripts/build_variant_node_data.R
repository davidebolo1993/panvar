#!/usr/bin/env Rscript
# Assemble the data bundle for the interactive node-coverage + variant-track viewer
# (scripts/variant_node_heatmap_app.R). Outputs an .rds with:
#   nodes     : node_id, length, gpos (reference genomic coordinate; inserted nodes projected to their
#               insertion locus), order (rank of gpos = genome order), gene, bubble_id
#   coverage  : haplotype, node_id, count   (per-walk traversal multiplicity; count>=1 only)
#   variants  : variant_id, bubble_id, svtype, pos, gene, nodes (union set from variant_nodes.tsv)
#   genotypes : variant_id, haplotype, gt, cn, cnbp
#   bubbles   : bubble_id, lo, hi   (gpos-order span of the bubble's nodes)
#   haplotypes: character; reference: the calling reference path name
#   Rscript build_variant_node_data.R --gfa <call.gfa> --variant-nodes <vn.tsv> --vcf <region.vcf> \
#       --bubbles <panphorte.bubbles.csv> [--node-genes <ng.tsv>] --out <bundle.rds>
suppressWarnings(suppressMessages(library(data.table)))
args <- commandArgs(trailingOnly = TRUE)
get <- function(f, d = NULL) { i <- match(f, args); if (is.na(i) || i == length(args)) d else args[i + 1] }
gfa <- get("--gfa"); vnp <- get("--variant-nodes"); vcfp <- get("--vcf")
bcsv <- get("--bubbles"); ngp <- get("--node-genes"); outp <- get("--out", "variant_node_data.rds")
stopifnot(!is.null(gfa), !is.null(vnp), !is.null(vcfp))

# --- GFA: node lengths (S) + ordered per-path node sequences (P) ---
L <- readLines(gfa)
S <- L[startsWith(L, "S\t")]
sf <- tstrsplit(sub("^S\t", "", S), "\t", fixed = TRUE)
nodes <- data.table(node_id = sf[[1]], length = nchar(sf[[2]]))
lenmap <- setNames(nodes$length, nodes$node_id)
paths_ids <- lapply(L[startsWith(L, "P\t")], function(p) {
  f <- strsplit(p, "\t", fixed = TRUE)[[1]]; sub("[+-]$", "", strsplit(f[3], ",", fixed = TRUE)[[1]]) })
names(paths_ids) <- vapply(L[startsWith(L, "P\t")], function(p) strsplit(p, "\t", fixed = TRUE)[[1]][2], character(1))
haplotypes <- names(paths_ids)
coverage <- rbindlist(lapply(haplotypes, function(h) {
  tb <- table(paths_ids[[h]]); data.table(haplotype = h, node_id = names(tb), count = as.integer(tb)) }))

# --- reference path (##reference) -> per-node genomic position; inserted nodes projected to the last
#     reference node seen before them along a walk (median across the haplotypes carrying them) ---
vl <- readLines(vcfp)
refname <- sub(".*##reference=", "", grep("^##reference=", vl, value = TRUE)[1])
refkey <- if (refname %in% haplotypes) refname else grep(refname, haplotypes, fixed = TRUE, value = TRUE)[1]
refpos <- rep(NA_real_, nrow(nodes)); names(refpos) <- nodes$node_id
{ cum <- 0; for (n in paths_ids[[refkey]]) { if (is.na(refpos[n])) refpos[n] <- cum; cum <- cum + lenmap[[n]] } }
proj <- rbindlist(lapply(haplotypes, function(h) {
  ids <- paths_ids[[h]]; rp <- refpos[ids]
  filled <- nafill(rp, "locf")
  data.table(node_id = ids, pos = filled)[is.na(refpos[node_id]) & !is.na(pos)] }))
gpos <- refpos
if (nrow(proj)) { pm <- proj[, .(pos = median(pos)), by = node_id]; gpos[pm$node_id] <- pm$pos }
gpos[is.na(gpos)] <- max(gpos, na.rm = TRUE) + 1  # unplaceable -> end
nodes[, gpos := gpos[node_id]]
setorder(nodes, gpos, node_id)
nodes[, order := .I]
order_of <- setNames(nodes$order, nodes$node_id)
setkey(nodes, node_id)

# --- genes ---
if (!is.null(ngp) && file.exists(ngp)) {
  ng <- fread(ngp, sep = "\t", header = TRUE, colClasses = "character"); nodes[ng, gene := i.genes, on = "node_id"]
}
if (!("gene" %in% names(nodes))) nodes[, gene := NA_character_]
nodes[is.na(gene), gene := "."]

# --- variant_nodes.tsv + POS/genotypes from the VCF ---
vn <- fread(vnp, sep = "\t", header = TRUE, colClasses = "character")
samples <- strsplit(vl[startsWith(vl, "#CHROM")][1], "\t", fixed = TRUE)[[1]][-(1:9)]
gt_rows <- vector("list", 0); vpos <- list(); vsvlen <- list()
for (bl in vl[!startsWith(vl, "#")]) {
  f <- strsplit(bl, "\t", fixed = TRUE)[[1]]; vid <- f[3]; vpos[[vid]] <- as.integer(f[2])
  m <- regmatches(f[8], regexpr("(^|;)SVLEN=[^;]+", f[8])); vsvlen[[vid]] <- if (length(m)) sub(".*SVLEN=", "", m) else "."
  keys <- strsplit(f[9], ":", fixed = TRUE)[[1]]; gi <- match("GT", keys); ci <- match("CN", keys); bi <- match("CNBP", keys)
  vals <- strsplit(f[-(1:9)], ":", fixed = TRUE)
  pick <- function(k) if (is.na(k)) rep(".", length(vals)) else vapply(vals, function(v) if (length(v) >= k) v[k] else ".", character(1))
  gt_rows[[length(gt_rows) + 1]] <- data.table(variant_id = vid, haplotype = samples, gt = pick(gi), cn = pick(ci), cnbp = pick(bi))
}
genotypes <- rbindlist(gt_rows)

# --- bubble node-order spans (from the bubbles' node sets, now in genomic order) ---
bub <- fread(bcsv, sep = ",", header = TRUE, colClasses = "character")
bubbles <- rbindlist(lapply(seq_len(nrow(bub)), function(i) {
  ids <- c(bub$source[i], bub$sink[i], strsplit(bub$inside_nodes[i], ";", fixed = TRUE)[[1]])
  o <- order_of[ids]; o <- o[!is.na(o)]; if (!length(o)) return(NULL)
  data.table(bubble_id = bub$bubble_id[i], lo = min(o), hi = max(o)) }))
nodes[, bubble_id := NA_character_]
for (i in seq_len(nrow(bubbles))) nodes[order >= bubbles$lo[i] & order <= bubbles$hi[i] & is.na(bubble_id), bubble_id := bubbles$bubble_id[i]]

gene_of <- setNames(nodes$gene, nodes$node_id)
variants <- vn[, .(variant_id, bubble_id, svtype, nodes = node_ids)]
variants[, pos := unlist(vpos[variant_id])]
variants[, svlen := unlist(vsvlen[variant_id])]
variants[, gene := vapply(strsplit(nodes, ",", fixed = TRUE), function(ns) {
  g <- unique(unlist(strsplit(gene_of[ns], ";", fixed = TRUE))); g <- g[!is.na(g) & g != "."]
  if (length(g)) paste(head(g, 3), collapse = ";") else "." }, character(1))]

saveRDS(list(nodes = nodes, coverage = coverage, variants = variants, genotypes = genotypes,
             bubbles = bubbles, haplotypes = haplotypes, reference = refkey), outp)
cat(sprintf("Wrote %s: %d nodes, %d coverage rows, %d variants, %d bubbles, %d haplotypes (ref %s)\n",
            outp, nrow(nodes), nrow(coverage), nrow(variants), nrow(bubbles), length(haplotypes), refkey))
