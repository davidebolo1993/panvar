#!/usr/bin/env bash
# LPA (KIV-2) driver: tandem repeat, called on the panphorte graph with --cn-from-multiplicity.
# Runs the data pipeline then the GWAS (cohort sim -> describe --samples -> associate, region scan +
# structure demo). Needs a numpy-capable python.
#   PYTHON=~/miniconda3/bin/python scripts/genes/lpa.sh        # env: N (cohort size), SIM (null markers)
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; source "$HERE/_common.sh"
region=lpa; d="$OUT/$region"

run_gene_data "$region" "$DATA/lpa.gfa.gz" 0.97 panphorte "--cn-from-multiplicity" 0.95 || exit 1

# ---- copy-number validation vs ground truth (optional; uncomment) -------------------------------
# "$PY" "$REPO/scripts/compare_copy_number.py" --vcf "$d/call/call.region.vcf" --label lpa \
#   --mode repeats --truth "$DATA/lpa.repeats.tsv" --name-col hapl --count-col copies

# ---- GWAS: structured cohort -> describe --samples -> associate ---------------------------------
N="${N:-300}"; SIM="${SIM:-3000}"; LMM_MAX_N="${LMM_MAX_N:-4000}"
ref="$(ref_of "$DATA/lpa.gfa.gz")"
pgfa="$d/panphorte/panphorte.normalized.sorted.gfa"; pfx="$d/panphorte/panphorte"
gw="$d/gwas"; real="$gw/real"; mkdir -p "$real"
KOPT=(); [ "$N" -le "$LMM_MAX_N" ] && KOPT=(--kinship-out "$real/kinship.tsv")
"$PY" "$HERE/../../tests/gwas/make_lpa_phenotype.py" "$d/panphorte/panphorte.panphorte.copies.tsv" \
  "$real" --n "$N" --sim-markers "$SIM" "${KOPT[@]}" || exit 1
NGOPT=(); [ -f "$d/call/call.node_genes.tsv" ] && NGOPT=(--node-genes "$d/call/call.node_genes.tsv")
"$BIN" describe -i "$pgfa" --bubble-prefix-in "$pfx" --out-dir "$gw/desc" --kmer-size 31 \
  --no-wide-matrix --variant-nodes "$d/call/call.variant_nodes.tsv" --samples "$real/samples.tsv" --quiet || exit 1

DESC="$gw/desc"; SAMP="$DESC/bimbam.samples.samples.txt.gz"; ANNOT="$DESC/feature_annot.samples.tsv.gz"
echo "[lpa] region scan (associate, PC-adjusted)"
for sub in graph kmers; do for mode in quant binary; do
  "$BIN" associate --genotypes "$DESC/bimbam_${sub}.samples.bimbam.gz" --samples "$SAMP" \
    --feature-annot "$ANNOT" "${NGOPT[@]}" --phenotype "$real/pheno.${mode}.tsv" --min-maf 0.02 \
    -o "$gw/assoc_${sub}_${mode}"
done; done
echo "[lpa] structure-correction demo (naive / PC / LMM) on the synthetic panel"
SG="$real/geno.sim.bimbam.gz"; SS="$real/sim.samples.txt"; SA="$real/feature_annot.sim.tsv.gz"
if [ -f "$SG" ]; then
  "$BIN" associate --genotypes "$SG" --samples "$SS" --feature-annot "$SA" --phenotype "$real/pheno.quant.nopc.tsv" --min-maf 0.02 -o "$gw/sim_naive"
  "$BIN" associate --genotypes "$SG" --samples "$SS" --feature-annot "$SA" --phenotype "$real/pheno.quant.tsv"      --min-maf 0.02 -o "$gw/sim_pc"
  [ "$N" -le "$LMM_MAX_N" ] && "$BIN" associate --genotypes "$SG" --samples "$SS" --feature-annot "$SA" \
    --phenotype "$real/pheno.quant.nopc.tsv" --model lmm --make-kinship --min-maf 0.02 -o "$gw/sim_lmm"
fi
echo "[lpa] GWAS done -> $gw"

# ---- plots (uncomment what you need) -----------------------------------------------------------
# cbub="$(main_dup_bubble "$d/call/call.region.vcf")"; ip="$d/inspect/inspect.bubble_${cbub}"
# "$RS" "$REPO/scripts/plot_vcf_map.R"            --vcf "$d/call/call.region.vcf" --reference-path "$ref" --flip --dpi 600 --title "lpa variant map" --out "$d/lpa_vcf_map_flipped"
# "$RS" "$REPO/scripts/plot_node_coverage_heatmap.R" --table "$ip.node_counts.tsv" --node-lengths "$ip.node_lengths.tsv" --cluster-by "$ip.clusters.tsv" --out "$d/lpa_node_heatmap"
# "$RS" "$REPO/scripts/plot_edge_coverage_heatmap.R" --table "$ip.edge_counts.tsv" --cluster-by "$ip.clusters.tsv" --out "$d/lpa_edge_heatmap"
# for s in assoc_graph_quant assoc_kmers_quant sim_naive sim_pc sim_lmm; do \
#   "$RS" "$REPO/scripts/plot_associate.R" --assoc "$gw/$s.assoc.tsv" --summary "$gw/$s.summary.tsv" --out "$gw/$s" --title "lpa $s"; done
