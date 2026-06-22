#!/usr/bin/env bash
# CYP2D6 per-gene driver -- DATA ONLY (no plots; plot commands are commented at the bottom).
# Topology: PGGB-collapsed paralog cluster (CYP2D6/2D7/2D8P) -> call on the BUBBLE graph with
# --cn-from-coverage (total-module copy number). panphorte runs at min-similarity 0.97.
#   scripts/genes/cyp2d6.sh
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; source "$HERE/_common.sh"
region=cyp2d6; d="$OUT/$region"

run_gene_data "$region" "$DATA/cyp2d6.gfa.gz" 0.97 bubble "--cn-from-coverage" 0.97 || exit 1

# ---- copy-number validation vs ground truth (optional; uncomment) -------------------------------
# "$PY" "$REPO/scripts/compare_copy_number.py" --vcf "$d/call/call.region.vcf" --label cyp2d6 \
#   --mode gene-count --truth "$DATA/cyp2d6.bed" --genes CYP2D6,CYP2D7
# per-gene split (reliable genes from --gtf realignment):
# for g in CYP2D6 CYP2D7; do "$PY" "$REPO/scripts/compare_copy_number.py" --mode per-gene \
#   --truth "$DATA/cyp2d6.bed" --genes "$g" --dup-gene-cn "$d/call/call.dup_gene_cn.tsv" --label "$g" --offset 0; done

# ---- plots (uncomment what you need) -----------------------------------------------------------
# ref="$(ref_of "$DATA/cyp2d6.gfa.gz")"; cbub="$(main_dup_bubble "$d/call/call.region.vcf")"
# ip="$d/inspect/inspect.bubble_${cbub}"
# "$RS" "$REPO/scripts/plot_vcf_map.R"            --vcf "$d/call/call.region.vcf" --reference-path "$ref" --flip --dpi 600 --title "cyp2d6 variant map" --out "$d/cyp2d6_vcf_map_flipped"
# "$RS" "$REPO/scripts/plot_vcf_map.R"            --vcf "$d/call/call.region.vcf" --reference-path "$ref" --clusters "$ip.clusters.tsv" --title "cyp2d6 (cluster reps)" --out "$d/cyp2d6_vcf_map_clustered"
# "$RS" "$REPO/scripts/plot_node_coverage_heatmap.R" --table "$ip.node_counts.tsv" --node-lengths "$ip.node_lengths.tsv" --cluster-by "$ip.clusters.tsv" --out "$d/cyp2d6_node_heatmap"
# "$RS" "$REPO/scripts/plot_edge_coverage_heatmap.R" --table "$ip.edge_counts.tsv" --cluster-by "$ip.clusters.tsv" --out "$d/cyp2d6_edge_heatmap"
