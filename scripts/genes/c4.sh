#!/usr/bin/env bash
# C4 (RCCX) per-gene driver -- DATA ONLY (no plots; plot commands are commented at the bottom).
# Topology: PGGB-collapsed paralog cluster -> --cn coverage route (total-module CN = C4A+C4B).
# Called on the panphorte graph (universal substrate; panphorte leaves this cluster untouched at
# min-similarity 0.97, so it equals the bubble graph).
#   scripts/genes/c4.sh
if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  echo "C4 (RCCX) driver: PGGB-collapsed paralog cluster; module CN + C4A/C4B per-gene split."
  echo "Usage: [PANVAR_BIN=..] [PYTHON=..] [RSCRIPT=..] [THREADS=..] scripts/genes/c4.sh"
  echo "Writes results/real_data/c4/. See scripts/genes/_common.sh for the env vars."
  exit 0
fi
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; source "$HERE/_common.sh"
region=c4; d="$OUT/$region"

run_gene_data "$region" "$DATA/c4.gfa.gz" 0.95 panphorte "--cn" 0.97 || exit 1

# ---- copy-number validation vs ground truth (optional; uncomment) -------------------------------
# "$PY" "$REPO/scripts/compare_copy_number.py" --vcf "$d/call/call.region.vcf" --label c4 \
#   --mode gene-count --truth "$DATA/c4.bed" --genes C4A,C4B

# ---- plots (uncomment what you need) -----------------------------------------------------------
# ref="$(ref_of "$DATA/c4.gfa.gz")"; cbub="$(main_dup_bubble "$d/call/call.region.vcf")"
# ip="$d/inspect/inspect.bubble_${cbub}"
# "$RS" "$REPO/scripts/plot_vcf_map.R"            --vcf "$d/call/call.region.vcf" --reference-path "$ref" --flip --dpi 600 --title "c4 variant map" --out "$d/c4_vcf_map_flipped"
# "$RS" "$REPO/scripts/plot_vcf_map.R"            --vcf "$d/call/call.region.vcf" --reference-path "$ref" --clusters "$ip.clusters.tsv" --title "c4 (cluster reps)" --out "$d/c4_vcf_map_clustered"
# "$RS" "$REPO/scripts/plot_node_coverage_heatmap.R" --table "$ip.node_counts.tsv" --node-lengths "$ip.node_lengths.tsv" --cluster-by "$ip.clusters.tsv" --out "$d/c4_node_heatmap"
# "$RS" "$REPO/scripts/plot_edge_coverage_heatmap.R" --table "$ip.edge_counts.tsv" --cluster-by "$ip.clusters.tsv" --out "$d/c4_edge_heatmap"
