#!/usr/bin/env bash
# Regenerate everything under results/ (gitignored): the full pipeline per region + copy-number
# validation vs the committed ground truth, plus the synthetic smoke. Needs Rscript for plots and a
# numpy/scipy python for the LPA GWAS (set RSCRIPT / PYTHON, e.g. conda base).
#   scripts/regen_results.sh [region ...]      # default: all regions + synthetic
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/.." && pwd)"
BIN="${PANVAR_BIN:-$REPO/build/panvar}"
PY="${PYTHON:-python3}"
RS="${RSCRIPT:-Rscript}"
DATA="$REPO/tests/real_data"
OUT="$REPO/results/real_data"
THREADS="${THREADS:-0}"
# optional GTF for gene annotation; skipped if absent
GTF="${GTF:-$DATA/Homo_sapiens.GRCh38.116.gtf.gz}"
GTFOPT=(); [[ -f "$GTF" ]] && GTFOPT=(--gtf "$GTF")

have_r() { command -v "$RS" >/dev/null 2>&1; }
ref_of() {  # pick a reference path: prefer GRCh38, then CHM13, else the first path
  gzcat "$1" 2>/dev/null | awk -F'\t' '($1=="P"||$1=="W"){n=$2; if(n~/[Gg][Rr][Cc]h38/){print n;exit} if(!first)first=n} END{if(first&&!done)print first}' | head -1
}

# run_region <region> <gfa.gz> <panphorte_extra> <call_graph: bubble|panphorte> <call_extra>
# Which graph `call` reads depends on topology; see docs/algorithms/call.md (CN-topology table).
run_region() {
  local region="$1" gfa="$2" pan_extra="$3" call_graph="$4" call_extra="$5" finer_clusters="${6:-0}"
  local d="$OUT/$region"
  echo "############ $region ############"
  mkdir -p "$d"/{bubble,panphorte,call,describe,inspect,plots}
  local ref; ref="$(ref_of "$gfa")"
  echo "[$region] reference: $ref ; call graph: $call_graph"

  # 1) bubble: internal sort+flip along the reference + cactus snarl finder. Writes
  #    bubble.sorted.gfa (the reference-oriented graph) which the rest of the pipeline consumes.
  "$BIN" bubble -i "$gfa" -o "$d/bubble/bubble" -r "$ref" "${GTFOPT[@]}" --quiet || return 1
  local sgfa="$d/bubble/bubble.sorted.gfa"
  # 2) panphorte on the sorted graph (normalized + re-sorted along the reference)
  "$BIN" panphorte -i "$sgfa" --bubble-prefix-in "$d/bubble/bubble" -o "$d/panphorte/panphorte" \
    --reference-path "$ref" --threads "$THREADS" "${GTFOPT[@]}" --quiet $pan_extra || return 1
  local pgfa="$d/panphorte/panphorte.normalized.sorted.gfa"
  # Pick the graph + bubble prefix `call`/`describe`/`inspect` read, per topology.
  local cgfa cpfx
  if [[ "$call_graph" == "bubble" ]]; then cgfa="$sgfa"; cpfx="$d/bubble/bubble";
  else cgfa="$pgfa"; cpfx="$d/panphorte/panphorte"; fi
  # 3) call
  "$BIN" call -i "$cgfa" --bubble-prefix-in "$cpfx" --reference-path "$ref" \
    -o "$d/call/call" --threads "$THREADS" "${GTFOPT[@]}" --quiet $call_extra || return 1
  # 4) describe (variant-restricted markers)
  "$BIN" describe -i "$cgfa" --bubble-prefix-in "$cpfx" --out-dir "$d/describe" \
    --variant-nodes "$d/call/call.variant_nodes.tsv" --no-wide-matrix --threads "$THREADS" --quiet || return 1
  # 4b) benchmark: round-trip QV of the calls (cosigt-style bands), if variant nodes were written
  if [[ -s "$d/call/call.variant_nodes.tsv" ]]; then
    "$BIN" benchmark -i "$cgfa" --bubble-prefix-in "$cpfx" --reference-path "$ref" \
      --variant-nodes "$d/call/call.variant_nodes.tsv" -o "$d/call/benchmark" \
      --threads "$THREADS" --quiet || echo "  (benchmark skipped)"
  fi
  # 5) inspect the largest-DUP bubble (the CN module) on the call graph, with walk clustering
  local bub; bub="$(awk -F'\t' '/SVTYPE=DUP/{
      sv=$8; sub(/.*SVLEN=/,"",sv); sub(/;.*/,"",sv); if(sv<0)sv=-sv;
      b=$8;  sub(/.*BUBBLE_ID=/,"",b);  sub(/;.*/,"",b);
      if(sv+0>m){m=sv+0;bb=b}} END{print bb}' "$d/call/call.region.vcf")"
  if [[ -n "${bub:-}" ]]; then
    "$BIN" inspect -i "$cgfa" --bubble-prefix-in "$cpfx" --bubble-id "$bub" \
      --cluster --cluster-similarity 0.97 -o "$d/inspect/inspect" --quiet || true
  fi
  # 6) plots (skipped if Rscript missing): whole-VCF variant map + per-bubble coverage heatmaps
  if have_r; then
    "$RS" "$HERE/plot_vcf_map.R" --vcf "$d/call/call.region.vcf" --reference-path "$ref" \
      --out "$d/plots/${region}_vcf_map" >/dev/null 2>&1 || echo "  (plot_vcf_map skipped)"
    "$RS" "$HERE/plot_vcf_map.R" --vcf "$d/call/call.region.vcf" --reference-path "$ref" --flip --dpi 600 \
      --out "$d/plots/${region}_vcf_map_flipped" >/dev/null 2>&1 || echo "  (plot_vcf_map flipped skipped)"
    if [[ -n "${bub:-}" ]]; then
      local ip="$d/inspect/inspect.bubble_${bub}"
      local nc="$ip.node_counts.tsv" ec="$ip.edge_counts.tsv" nl="$ip.node_lengths.tsv" cl="$ip.clusters.tsv"
      # representative-map clusters: call-substrate by default; finer_clusters=1 uses the pre-panphorte
      # bubble graph instead (LPA's collapsed graph has too few walk clusters to show)
      local cl_map="$cl"
      if [[ "$finer_clusters" == "1" ]]; then
        local bbub; bbub="$(awk -F',' 'NR>1 && $8+0>m{m=$8+0; b=$1} END{print b}' "$d/bubble/bubble.bubbles.csv")"
        if [[ -n "$bbub" ]]; then
          "$BIN" inspect -i "$sgfa" --bubble-prefix-in "$d/bubble/bubble" --bubble-id "$bbub" \
            --cluster --cluster-similarity 0.97 -o "$d/inspect/inspect_bubblestage" --quiet || true
          local clb="$d/inspect/inspect_bubblestage.bubble_${bbub}.clusters.tsv"
          [[ -f "$clb" ]] && cl_map="$clb"
        fi
      fi
      # cluster-representative-only variant maps (both orientations): one row per walk cluster.
      if [[ -f "$cl_map" ]]; then
        "$RS" "$HERE/plot_vcf_map.R" --vcf "$d/call/call.region.vcf" --reference-path "$ref" \
          --clusters "$cl_map" \
          --out "$d/plots/${region}_vcf_map_clustered" >/dev/null 2>&1 || echo "  (plot_vcf_map clustered skipped)"
        "$RS" "$HERE/plot_vcf_map.R" --vcf "$d/call/call.region.vcf" --reference-path "$ref" --flip \
          --clusters "$cl_map" \
          --out "$d/plots/${region}_vcf_map_clustered_flipped" >/dev/null 2>&1 || echo "  (plot_vcf_map clustered flipped skipped)"
      fi
      if [[ -f "$nc" ]]; then
        "$RS" "$HERE/plot_node_coverage_heatmap.R" --table "$nc" --node-lengths "$nl" \
          --cluster-by "$cl" --out "$d/plots/${region}_node_heatmap" >/dev/null 2>&1 || echo "  (node heatmap skipped)"
      fi
      [[ -f "$ec" ]] && { "$RS" "$HERE/plot_edge_coverage_heatmap.R" --table "$ec" --cluster-by "$cl" \
        --out "$d/plots/${region}_edge_heatmap" >/dev/null 2>&1 || echo "  (edge heatmap skipped)"; }
    fi
  else
    echo "  (R not found; skipping plots)"
  fi
  echo "[$region] done -> $d"
}

CN_TABLE="$REPO/results/cn_table.tsv"   # combined per-haplotype table for the concordance plot

validate_cn() {  # <region> <compare args...>
  local region="$1"; shift
  echo "---- $region copy-number vs ground truth ----"
  "$PY" "$HERE/compare_copy_number.py" --vcf "$OUT/$region/call/call.region.vcf" --label "$region" \
    --dump-tsv "$CN_TABLE" "$@" \
    || echo "  (validation could not run)"
}

REGIONS=("$@"); [[ ${#REGIONS[@]} -eq 0 ]] && REGIONS=(lpa c4 gstm1 cyp2d6 acot synthetic)
mkdir -p "$REPO/results"; rm -f "$CN_TABLE"

for region in "${REGIONS[@]}"; do
  case "$region" in
    lpa)
      # Tandem repeat (KIV-2): panphorte folds the variable tandem; call counts on the REP node.
      run_region lpa "$DATA/lpa.gfa.gz" "--min-similarity 0.70" panphorte "--cn" 1
      validate_cn lpa --mode repeats --truth "$DATA/lpa.repeats.tsv" --name-col hapl --count-col copies
      # LPA GWAS demo (needs python numpy/scipy); reuses the dedicated driver
      bash "$HERE/../tests/gwas/run_lpa_real.sh" "$BIN" "$OUT/lpa/gwas" "$PY" "$RS" || echo "  (gwas demo skipped)"
      ;;
    c4)
      # PGGB-collapsed paralog cluster (each RCCX module = one C4 gene). Total CN from folded-node
      # coverage -> 131/131 exact against the C4A+C4B gene count. Called on the panphorte graph (the
      # universal substrate; panphorte leaves this cluster untouched, so it equals the bubble graph).
      run_region c4 "$DATA/c4.gfa.gz" "--min-similarity 0.70" panphorte "--cn"
      validate_cn c4 --mode gene-count --truth "$DATA/c4.bed" --genes C4A,C4B
      # Resolved per-gene split (C4A/C4B from the k-mer per-site resolver, --gtf) vs per-gene truth,
      # absolute (offset 0). C4A/C4B are near-identical, so the split rides the per-site consensus.
      if [[ -s "$OUT/c4/call/call.dup_gene_cn.tsv" ]]; then
        for g in C4A C4B; do
          echo "---- c4 $g per-gene split vs ground truth ----"
          "$PY" "$HERE/compare_copy_number.py" --mode per-gene --truth "$DATA/c4.bed" --genes "$g" \
            --dup-gene-cn "$OUT/c4/call/call.dup_gene_cn.tsv" --label "$g" --offset 0 \
            --dump-tsv "$CN_TABLE" || echo "  (per-gene $g compare skipped)"
        done
      fi
      ;;
    gstm1)
      # PGGB-collapsed paralog cluster (GSTM1/2/5, reference folds 3x): total-module coverage,
      # validated against the GSTM1 count after the paralog baseline. Called on the panphorte graph.
      run_region gstm1 "$DATA/gstm1.gfa.gz" "--min-similarity 0.70" panphorte "--cn"
      validate_cn gstm1 --mode direct --truth "$DATA/gstm1.bed"
      # Per-gene split (--gtf) vs the full pangene truth (gstm.bed, per-molecule gene counts), absolute
      # (offset 0). We resolve GSTM1/2/4/5; GSTM3 is stable single-copy outside the module, so not scored.
      if [[ -s "$OUT/gstm1/call/call.dup_gene_cn.tsv" ]]; then
        for g in GSTM1 GSTM2 GSTM4 GSTM5; do
          echo "---- gstm1 $g per-gene split vs ground truth ----"
          "$PY" "$HERE/compare_copy_number.py" --mode per-gene --truth "$DATA/gstm.bed" --genes "$g" \
            --dup-gene-cn "$OUT/gstm1/call/call.dup_gene_cn.tsv" --label "$g" --offset 0 \
            --dump-tsv "$CN_TABLE" || echo "  (per-gene $g compare skipped)"
        done
      fi
      ;;
    cyp2d6)
      # PGGB-collapsed paralog cluster (CYP2D6/2D7/2D8P). Folded-node coverage = total module copies;
      # validated against the collapsed D6+D7 truth. The few residual misses carry an extra
      # (unannotated) CYP2D8P / 2D7-hybrid module that the gene BED does not count -- real biology.
      # Called on the panphorte graph (leaves this cluster untouched, so it equals the bubble graph).
      run_region cyp2d6 "$DATA/cyp2d6.gfa.gz" "" panphorte "--cn"
      validate_cn cyp2d6 --mode gene-count --truth "$DATA/cyp2d6.bed" --genes CYP2D6,CYP2D7
      # Resolved per-gene split (divergent paralogs, --gtf private-k-mer dosage) vs per-gene truth: CYP2D6
      # and CYP2D7 separately, absolute (offset 0). Adds CYP2D6/CYP2D7 facets to the concordance plot.
      if [[ -s "$OUT/cyp2d6/call/call.dup_gene_cn.tsv" ]]; then
        for g in CYP2D6 CYP2D7; do
          echo "---- cyp2d6 $g per-gene split vs ground truth ----"
          "$PY" "$HERE/compare_copy_number.py" --mode per-gene --truth "$DATA/cyp2d6.bed" --genes "$g" \
            --dup-gene-cn "$OUT/cyp2d6/call/call.dup_gene_cn.tsv" --label "$g" --offset 0 \
            --dump-tsv "$CN_TABLE" || echo "  (per-gene $g compare skipped)"
        done
      fi
      ;;
    acot)
      # chr14 ACOT paralog cluster: panphorte folds ACOT1/ACOT2 into one ~14.5 kb RU. Total module CN vs
      # the ACOT1+ACOT2 gene count is exact, and the private-k-mer split (--gtf) resolves the two variable
      # paralogs -- ACOT1 (a common presence/null, 0-3 copies) and ACOT2 (mostly single). ACOT4/ACOT6 are
      # stable single-copy outside the folded module (like GSTM3), so they are not resolved.
      run_region acot "$DATA/acot.gfa.gz" "--min-similarity 0.70" panphorte "--cn"
      validate_cn acot --mode gene-count --truth "$DATA/acot.bed" --genes ACOT1,ACOT2
      if [[ -s "$OUT/acot/call/call.dup_gene_cn.tsv" ]]; then
        for g in ACOT1 ACOT2; do
          echo "---- acot $g per-gene split vs ground truth ----"
          "$PY" "$HERE/compare_copy_number.py" --mode per-gene --truth "$DATA/acot.bed" --genes "$g" \
            --dup-gene-cn "$OUT/acot/call/call.dup_gene_cn.tsv" --label "$g" --offset 0 \
            --dump-tsv "$CN_TABLE" || echo "  (per-gene $g compare skipped)"
        done
      fi
      ;;
    synthetic)
      echo "############ synthetic ############"
      mkdir -p "$REPO/results/synthetic_data"
      bash "$REPO/tests/synthetic_smoke.sh" "$BIN" "$REPO/tests/synthetic_data" \
        "$REPO/results/synthetic_data" || echo "  (synthetic smoke failed)"
      ;;
    *) echo "unknown region: $region" ;;
  esac
done

# Copy-number concordance plot (called-vs-truth, faceted by gene) from the combined table.
if [[ -s "$CN_TABLE" ]] && have_r; then
  "$RS" "$HERE/plot_cn_correlation.R" --table "$CN_TABLE" --out "$REPO/results/cn_correlation" \
    || echo "  (cn_correlation plot skipped)"
elif [[ -s "$CN_TABLE" ]]; then
  echo "  (R not found; skipping cn_correlation plot; table at $CN_TABLE)"
fi

# Round-trip QV summary: combine each locus's per-haplotype benchmark, then a cosigt-style stacked bar.
QV_TABLE="$REPO/results/benchmark_qv.tsv"
printf 'locus\tsample\tsum_delta\tsum_aln_len\tqv\tband\tqv_ratio\tquintile\tidentity\n' > "$QV_TABLE"
for region in "${REGIONS[@]}"; do
  bh="$OUT/$region/call/benchmark.qv_by_haplotype.tsv"
  [[ -s "$bh" ]] && awk -F'\t' -v L="$region" 'NR>1{print L"\t"$1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$8"\t"$9"\t"$10}' "$bh" >> "$QV_TABLE"
done
if [[ $(wc -l < "$QV_TABLE") -gt 1 ]] && have_r; then
  for cat in dist quintile band identity; do
    "$RS" "$HERE/plot_benchmark_qv.R" --table "$QV_TABLE" --out "$REPO/results/benchmark_qv" --category "$cat" \
      || echo "  (benchmark_qv $cat plot skipped)"
  done
fi
echo "ALL DONE -> $REPO/results"
