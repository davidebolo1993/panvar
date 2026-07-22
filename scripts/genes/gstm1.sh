#!/usr/bin/env bash
# GSTM1 per-gene driver -- DATA ONLY (no plots; plot commands are commented at the bottom).
# Topology: PGGB-collapsed paralog cluster (GSTM1/2/5, reference folds 3x) -> --cn coverage route
# (total-module CN, validated against the GSTM1 count after the paralog baseline). Called on the
# panphorte graph (universal substrate; panphorte leaves this cluster untouched, so it equals bubble).
#   scripts/genes/gstm1.sh
if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  echo "GSTM1 driver: PGGB-collapsed paralog cluster; total-module coverage CN + per-gene split."
  echo "Usage: [PANVAR_BIN=..] [PYTHON=..] [RSCRIPT=..] [THREADS=..] scripts/genes/gstm1.sh"
  echo "Writes results/real_data/gstm1/. See scripts/genes/_common.sh for the env vars."
  exit 0
fi
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; source "$HERE/_common.sh"
region=gstm1; d="$OUT/$region"

run_gene_data "$region" "$DATA/gstm1.gfa.gz" 0.95 panphorte "--cn" 0.97 || exit 1

# ---- copy-number validation vs ground truth (optional; uncomment) -------------------------------
# "$PY" "$REPO/scripts/compare_copy_number.py" --vcf "$d/call/call.region.vcf" --label gstm1 \
#   --mode direct --truth "$DATA/gstm1.bed"

# ---- plots (uncomment what you need) -----------------------------------------------------------
# ref="$(ref_of "$DATA/gstm1.gfa.gz")"; cbub="$(main_dup_bubble "$d/call/call.region.vcf")"
# ip="$d/inspect/inspect.bubble_${cbub}"
# "$RS" "$REPO/scripts/plot_vcf_map.R"            --vcf "$d/call/call.region.vcf" --reference-path "$ref" --flip --dpi 600 --title "gstm1 variant map" --out "$d/gstm1_vcf_map_flipped"
# "$RS" "$REPO/scripts/plot_vcf_map.R"            --vcf "$d/call/call.region.vcf" --reference-path "$ref" --clusters "$ip.clusters.tsv" --title "gstm1 (cluster reps)" --out "$d/gstm1_vcf_map_clustered"
# "$RS" "$REPO/scripts/plot_node_coverage_heatmap.R" --table "$ip.node_counts.tsv" --node-lengths "$ip.node_lengths.tsv" --cluster-by "$ip.clusters.tsv" --out "$d/gstm1_node_heatmap"
# "$RS" "$REPO/scripts/plot_edge_coverage_heatmap.R" --table "$ip.edge_counts.tsv" --cluster-by "$ip.clusters.tsv" --out "$d/gstm1_edge_heatmap"
