#!/usr/bin/env bash
# Regenerate all panvar results into results/ (kept local; gitignored).
#
# For each real gene region (lpa, c4, gstm1, cyp2d6) runs the full pipeline
#   bubble -> panphorte -> call -> describe -> inspect -> plots -> copy-number validation
# into results/real_data/<region>/<module>/<module>.*, and the synthetic smoke into
# results/synthetic_data/. Copy number is checked against the committed ground truth
# (tests/real_data/{c4,cyp2d6,gstm1}.bed, lpa.repeats.tsv).
#
# R plots use Rscript; if you use conda, run `conda activate base` first (or set RSCRIPT).
# Python needs numpy+scipy for the LPA GWAS (set PYTHON to an env that has them).
#
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

have_r() { command -v "$RS" >/dev/null 2>&1; }
ref_of() {  # pick a reference path: prefer GRCh38, then CHM13, else the first path
  gzcat "$1" 2>/dev/null | awk -F'\t' '($1=="P"||$1=="W"){n=$2; if(n~/[Gg][Rr][Cc]h38/){print n;exit} if(!first)first=n} END{if(first&&!done)print first}' | head -1
}

# run_region <region> <gfa.gz> <panphorte_extra> <call_graph: bubble|panphorte> <call_extra>
#
# The graph `call` consumes depends on the region's topology (see docs/references.md / the CN story):
#   * panphorte : tandem-repeat regions (LPA KIV-2) - panphorte folds the variable tandem into a REP
#                 node, then `call --cn-from-multiplicity` counts copies on that graph.
#   * bubble    : PGGB-collapsed paralog clusters (C4, CYP2D6) - the copies live as node multiplicity
#                 on shared nodes, so `call --cn-from-coverage` counts them over the FULL bubble walk
#                 of the *unfolded* (bubble) graph; panphorte folding would hide that multiplicity.
# panphorte is still run for every region (its copies.tsv carries C4's long/short composition), but
# for the bubble-call regions its graph is not what `call` reads.
run_region() {
  local region="$1" gfa="$2" pan_extra="$3" call_graph="$4" call_extra="$5"
  local d="$OUT/$region"
  echo "############ $region ############"
  mkdir -p "$d"/{bubble,panphorte,call,describe,inspect,plots}
  local ref; ref="$(ref_of "$gfa")"
  echo "[$region] reference: $ref ; call graph: $call_graph"

  # 1) bubble: internal sort+flip along the reference + cactus snarl finder. Writes
  #    bubble.sorted.gfa (the reference-oriented graph) which the rest of the pipeline consumes.
  "$BIN" bubble -i "$gfa" -o "$d/bubble/bubble" -r "$ref" --quiet || return 1
  local sgfa="$d/bubble/bubble.sorted.gfa"
  # 2) panphorte on the sorted graph (normalized + re-sorted along the reference)
  "$BIN" panphorte -i "$sgfa" --bubble-prefix-in "$d/bubble/bubble" -o "$d/panphorte/panphorte" \
    --reference-path "$ref" --threads "$THREADS" --quiet $pan_extra || return 1
  local pgfa="$d/panphorte/panphorte.normalized.sorted.gfa"
  # Pick the graph + bubble prefix `call`/`describe`/`inspect` read, per topology.
  local cgfa cpfx
  if [[ "$call_graph" == "bubble" ]]; then cgfa="$sgfa"; cpfx="$d/bubble/bubble";
  else cgfa="$pgfa"; cpfx="$d/panphorte/panphorte"; fi
  # 3) call
  "$BIN" call -i "$cgfa" --bubble-prefix-in "$cpfx" --reference-path "$ref" \
    -o "$d/call/call" --threads "$THREADS" --quiet $call_extra || return 1
  # 4) describe (variant-restricted markers)
  "$BIN" describe -i "$cgfa" --bubble-prefix-in "$cpfx" --out-dir "$d/describe" \
    --variant-nodes "$d/call/call.variant_nodes.tsv" --no-wide-matrix --threads "$THREADS" --quiet || return 1
  # 5) inspect the main event bubble on the call substrate, WITH walk clustering, so the node ids line
  #    up 1:1 with the call VCF and rows can be ordered/representative-reduced by cluster. The "main"
  #    bubble is the one with the LARGEST DUP (|SVLEN|) — the copy-number module the comparator scores —
  #    not merely the first DUP in coordinate order.
  local bub; bub="$(awk -F'\t' '/SVTYPE=DUP/{
      sv=$8; sub(/.*SVLEN=/,"",sv); sub(/;.*/,"",sv); if(sv<0)sv=-sv;
      b=$8;  sub(/.*BUBBLE_ID=/,"",b);  sub(/;.*/,"",b);
      if(sv+0>m){m=sv+0;bb=b}} END{print bb}' "$d/call/call.region.vcf")"
  if [[ -n "${bub:-}" ]]; then
    "$BIN" inspect -i "$cgfa" --bubble-prefix-in "$cpfx" --bubble-id "$bub" \
      --cluster --cluster-similarity 0.97 -o "$d/inspect/inspect" --quiet || true
  fi
  # 6) plots (R; skipped if Rscript missing). The headline is the whole-VCF variant map (oncoprint
  #    style: haplotypes x called variants, grouped by bubble), which needs only the VCF. The per-bubble
  #    node/edge coverage heatmaps + node-order SV map (odgi-viz style, ordered by walk cluster) give the
  #    graph-structure view of the main event bubble.
  if have_r; then
    # whole-VCF variant map in BOTH orientations: default (haplotypes on Y) and flipped (variants on Y).
    "$RS" "$HERE/plot_vcf_map.R" --vcf "$d/call/call.region.vcf" --reference-path "$ref" \
      --title "$region variant map" --out "$d/plots/${region}_vcf_map" >/dev/null 2>&1 || echo "  (plot_vcf_map skipped)"
    "$RS" "$HERE/plot_vcf_map.R" --vcf "$d/call/call.region.vcf" --reference-path "$ref" --flip \
      --title "$region variant map" --out "$d/plots/${region}_vcf_map_flipped" >/dev/null 2>&1 || echo "  (plot_vcf_map flipped skipped)"
    if [[ -n "${bub:-}" ]]; then
      local ip="$d/inspect/inspect.bubble_${bub}"
      local nc="$ip.node_counts.tsv" ec="$ip.edge_counts.tsv" nl="$ip.node_lengths.tsv" cl="$ip.clusters.tsv"
      if [[ -f "$nc" ]]; then
        "$RS" "$HERE/plot_node_coverage_heatmap.R" --table "$nc" --node-lengths "$nl" \
          --cluster-by "$cl" --out "$d/plots/${region}_node_heatmap" >/dev/null 2>&1 || echo "  (node heatmap skipped)"
        "$RS" "$HERE/plot_sv_map.R" --node-counts "$nc" --vcf "$d/call/call.region.vcf" \
          --node-lengths "$nl" --cluster-by "$cl" --reference-path "$ref" \
          --out "$d/plots/${region}_sv_map" >/dev/null 2>&1 || echo "  (plot_sv_map skipped)"
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

REGIONS=("$@"); [[ ${#REGIONS[@]} -eq 0 ]] && REGIONS=(lpa c4 gstm1 cyp2d6 synthetic)
mkdir -p "$REPO/results"; rm -f "$CN_TABLE"

for region in "${REGIONS[@]}"; do
  case "$region" in
    lpa)
      # Tandem repeat (KIV-2): panphorte folds the variable tandem; call counts on the REP node.
      run_region lpa "$DATA/lpa.gfa.gz" "--min-similarity 0.70" panphorte "--cn-from-multiplicity"
      validate_cn lpa --mode repeats --truth "$DATA/lpa.repeats.tsv" --name-col hapl --count-col copies
      # LPA GWAS demo (needs python numpy/scipy); reuses the dedicated driver
      bash "$HERE/../tests/gwas/run_lpa_real.sh" "$BIN" "$OUT/lpa/gwas" "$PY" "$RS" || echo "  (gwas demo skipped)"
      ;;
    c4)
      # PGGB-collapsed paralog cluster (each RCCX module = one C4 gene). Total CN from full-walk
      # coverage on the bubble graph -> 131/131 exact against the C4A+C4B gene count.
      run_region c4 "$DATA/c4.gfa.gz" "--min-similarity 0.70" bubble "--cn-from-coverage"
      validate_cn c4 --mode gene-count --truth "$DATA/c4.bed" --genes C4A,C4B
      ;;
    gstm1)
      # Deletion/CNV (ref CN=1): peak node multiplicity over the segdup cluster (bubble graph).
      run_region gstm1 "$DATA/gstm1.gfa.gz" "--min-similarity 0.70" bubble "--cn-from-multiplicity"
      validate_cn gstm1 --mode direct --truth "$DATA/gstm1.bed"
      ;;
    cyp2d6)
      # PGGB-collapsed paralog cluster (CYP2D6/2D7/2D8P). Full-walk coverage = total module copies;
      # validated against the collapsed D6+D7 truth. The few residual misses carry an extra
      # (unannotated) CYP2D8P / 2D7-hybrid module that the gene BED does not count -- real biology.
      run_region cyp2d6 "$DATA/cyp2d6.gfa.gz" "" bubble "--cn-from-coverage"
      validate_cn cyp2d6 --mode gene-count --truth "$DATA/cyp2d6.bed" --genes CYP2D6,CYP2D7
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
echo "ALL DONE -> $REPO/results"
