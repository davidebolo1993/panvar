#!/usr/bin/env bash
# publish_plots.sh - copy the plots the documentation embeds from results/ into docs/img/.
#
# The docs reference docs/img/*.png; the pipeline writes results/**. Keeping the two in step by hand
# is how they drift: docs/img/benchmark_qv.png sat two weeks behind results/benchmark_qv.png, so the
# published figure described a run nobody could reproduce. Every copy goes through this script, and it
# refuses rather than publishing a figure older than the run that is supposed to have produced it.
#
# Run it AFTER a successful scripts/regen_results.sh, from a clean worktree.
#
#   scripts/publish_plots.sh [--check]
#
#   --check   report what is stale or missing and exit non-zero; copy nothing. For a release gate.
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/.." && pwd)"
RES="$REPO/results"
IMG="$REPO/docs/img"
CHECK=0
[[ "${1:-}" == "--check" ]] && CHECK=1

# <source under results/>  <name under docs/img/>
MAP=(
  "benchmark_qv.png|benchmark_qv.png"
  "cn_correlation.genes.png|cn_correlation.genes.png"
  "cn_correlation.loci.png|cn_correlation.loci.png"
  "real_data/lpa/plots/lpa_vcf_map.png|lpa_vcf_map.png"
  "real_data/lpa/plots/lpa_node_heatmap.png|lpa_inspect_node_heatmap.png"
  "real_data/lpa/gwas/associate/assoc_graph_quant.manhattan.png|assoc_graph_quant.manhattan.png"
  "real_data/lpa/gwas/associate/assoc_graph_quant.pipeline.png|assoc_graph_quant.pipeline.png"
  "real_data/lpa/gwas/associate/assoc_kmers_quant.manhattan.png|assoc_kmers_quant.manhattan.png"
  "real_data/lpa/gwas/associate/assoc_kmers_quant.pipeline.png|assoc_kmers_quant.pipeline.png"
  "real_data/lpa/gwas/associate/assoc_variant_quant.manhattan.png|assoc_variant_quant.manhattan.png"
  "real_data/lpa/gwas/associate/assoc_variant_quant.pipeline.png|assoc_variant_quant.pipeline.png"
)

# Bandage screenshots (lpa_bubble_bandage, lpa_panphorte_bandage), the panphorte-stage heatmap, the
# node-coverage figure and the MYOM2 rebuild pair are produced by hand or by a tool this pipeline does
# not drive, so they are not published from here and are not checked for staleness.
MANUAL="lpa_bubble_bandage.png lpa_panphorte_bandage.png lpa_inspect_node_heatmap_panphorte.png
        lpa_node_coverage.png myom2.original.png myom2.rebuild.png"

mkdir -p "$IMG"
missing=0; stale=0; copied=0
for entry in "${MAP[@]}"; do
  src="$RES/${entry%%|*}"; dst="$IMG/${entry##*|}"
  if [[ ! -f "$src" ]]; then
    echo "MISSING in results/: ${entry%%|*}  (run scripts/regen_results.sh first)"
    missing=$((missing + 1)); continue
  fi
  if [[ "$CHECK" == "1" ]]; then
    if [[ ! -f "$dst" ]] || [[ "$src" -nt "$dst" ]] || ! cmp -s "$src" "$dst"; then
      echo "STALE: docs/img/${entry##*|} does not match results/${entry%%|*}"
      stale=$((stale + 1))
    fi
  else
    cmp -s "$src" "$dst" || { cp "$src" "$dst"; copied=$((copied + 1)); echo "published ${entry##*|}"; }
  fi
done

for m in $MANUAL; do
  [[ -f "$IMG/$m" ]] || { echo "MISSING hand-made figure: docs/img/$m"; missing=$((missing + 1)); }
done

if [[ "$CHECK" == "1" ]]; then
  if [[ $((missing + stale)) -eq 0 ]]; then echo "docs/img is in step with results/"; exit 0; fi
  echo "$stale stale, $missing missing"; exit 1
fi
echo "published $copied file(s); $missing missing"
[[ "$missing" -eq 0 ]]
