#!/usr/bin/env bash
# Regenerate everything under results/ (gitignored): the full pipeline per region + copy-number
# validation vs the committed ground truth, plus the synthetic smoke.
set -uo pipefail

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  cat <<'EOF'
regen_results.sh - regenerate everything under results/ (gitignored).

Usage:
  scripts/regen_results.sh [region ...]

  region ...   subset to run (default: lpa c4 gstm1 cyp2d6 acot ankrd36c synthetic)

For each region: bubble -> panphorte -> refine -> call -> describe -> benchmark, plus copy-number
validation against the committed ground truth, plus the cn-correlation and benchmark-QV plots.

Environment:
  PANVAR_BIN   panvar binary        (default build/panvar)
  PYTHON       python for the plots + LPA GWAS; needs numpy/scipy (default python3)
  RSCRIPT      Rscript for the plots (default Rscript)
  THREADS      worker threads       (default 0 = auto)
  REFINE       1 = call on the refined graph, 0 = on the panphorte graph (default 1)
  GTF          gene GTF for annotation (default tests/real_data/Homo_sapiens.GRCh38.116.gtf.gz)
  STRICT       1 = release mode (default): every required stage must succeed, the GTF and R must be
               present, and a locus is published only after all of its stages pass. 0 = exploratory:
               failures are reported and the run continues, which is what "skipped" used to mean.
EOF
  exit 0
fi

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/.." && pwd)"
BIN="${PANVAR_BIN:-$REPO/build/panvar}"
PY="${PYTHON:-python3}"
RS="${RSCRIPT:-Rscript}"
DATA="$REPO/tests/real_data"
OUT="$REPO/results/real_data"
THREADS="${THREADS:-0}"
# refine (POA-realign bubble interiors to remove pggb artifacts) runs in-pipeline by default; the
# calls come from the refined graph. Set REFINE=0 to fall back to calling on the panphorte graph.
REFINE="${REFINE:-1}"
# optional GTF for gene annotation; skipped if absent
GTF="${GTF:-$DATA/Homo_sapiens.GRCh38.116.gtf.gz}"
GTFOPT=(); [[ -f "$GTF" ]] && GTFOPT=(--gtf "$GTF")
STRICT="${STRICT:-1}"

have_r() { command -v "$RS" >/dev/null 2>&1; }

# In strict mode a failure is a failure. The previous script turned several into the word "skipped"
# and could still print ALL DONE, so a published tree could be missing whole stages while reading as a
# complete run -- and the stale files below made that indistinguishable from success.
FAILURES=0
die_or_warn() {   # die_or_warn <message>
  if [[ "$STRICT" == "1" ]]; then echo "regen_results: FAILED: $1" >&2; exit 1; fi
  echo "  ($1 -- continuing, STRICT=0)"; FAILURES=$((FAILURES + 1))
}

if [[ "$STRICT" == "1" ]]; then
  [[ -f "$GTF" ]] || { echo "regen_results: STRICT needs the GTF at $GTF (set GTF= or STRICT=0)" >&2; exit 1; }
  have_r || { echo "regen_results: STRICT needs Rscript ($RS) for the plots (set RSCRIPT= or STRICT=0)" >&2; exit 1; }
  command -v "$PY" >/dev/null 2>&1 || { echo "regen_results: STRICT needs python ($PY)" >&2; exit 1; }
  [[ -x "$BIN" ]] || { echo "regen_results: STRICT needs the panvar binary at $BIN" >&2; exit 1; }
fi

# What produced this tree. Without it a results directory cannot be tied to a commit or a binary, and
# an A/B against it is only as trustworthy as someone's memory of how it was made.
write_manifest() {
  local mf="$REPO/results/MANIFEST.txt"
  {
    echo "generated    $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
    echo "commit       $(git -C "$REPO" rev-parse HEAD 2>/dev/null || echo unknown)"
    echo "worktree     $(git -C "$REPO" diff --quiet 2>/dev/null && echo clean || echo DIRTY)"
    echo "binary       $BIN"
    echo "binary_sha   $( (shasum -a 256 "$BIN" 2>/dev/null || sha256sum "$BIN" 2>/dev/null) | awk '{print $1}')"
    echo "version      $("$BIN" --help 2>&1 | head -1)"
    echo "regions      ${REGIONS[*]}"
    echo "threads      $THREADS"
    echo "refine       $REFINE"
    echo "strict       $STRICT"
    echo "gtf          $GTF"
    echo "python       $PY ($("$PY" --version 2>&1))"
    echo "rscript      $RS ($("$RS" --version 2>&1 | head -1))"
  } > "$mf"
  echo "wrote $mf"
}
ref_of() {  # pick a reference path: prefer GRCh38, then CHM13, else the first path
  gunzip -c "$1" 2>/dev/null | awk -F'\t' '($1=="P"||$1=="W"){n=$2; if(n~/[Gg][Rr][Cc]h38/){print n;exit} if(!first)first=n} END{if(first&&!done)print first}' | head -1
}

# run_region <region> <gfa.gz> <panphorte_extra> <call_graph: bubble|panphorte> <call_extra>
run_region() {
  local region="$1" gfa="$2" pan_extra="$3" call_graph="$4" call_extra="$5" finer_clusters="${6:-0}"
  local d="$OUT/$region"
  echo "############ $region ############"
  # Start from an EMPTY locus directory. Reusing it let output from an older run survive alongside the
  # new one with nothing marking which was which -- measured: ANKRD36C kept an inspect.bubble_11 family
  # from a run whose bubble numbering no longer exists, and C4 kept a bubble_4 family from July. An A/B
  # against such a tree silently compares the current build with a months-old one.
  [[ -n "$region" && "$d" == "$OUT/"* ]] || { echo "refusing to clear '$d'" >&2; return 1; }
  rm -rf "$d"
  mkdir -p "$d"/{bubble,panphorte,refine,call,describe,inspect,benchmark,plots}
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
  local ppfx="$d/panphorte/panphorte"
  # 2b) refine (default on; REFINE=0 to skip): POA-realign bubble interiors to remove pggb alignment
  #     artifacts. DUP-safe (folded REP nodes held fixed); emits its own sorted GFA + bubbles the rest consume.
  if [[ "$REFINE" == "1" && "$call_graph" == "panphorte" ]]; then
    "$BIN" refine -i "$pgfa" --bubble-prefix-in "$ppfx" --reference-path "$ref" \
      -o "$d/refine/refine" "${GTFOPT[@]}" --quiet || return 1
    pgfa="$d/refine/refine.normalized.sorted.gfa"; ppfx="$d/refine/refine"
  fi
  # Pick the graph + bubble prefix `call`/`describe`/`inspect` read, per topology.
  local cgfa cpfx
  if [[ "$call_graph" == "bubble" ]]; then cgfa="$sgfa"; cpfx="$d/bubble/bubble";
  else cgfa="$pgfa"; cpfx="$ppfx"; fi
  # 3) call
  "$BIN" call -i "$cgfa" --bubble-prefix-in "$cpfx" --reference-path "$ref" \
    -o "$d/call/call" --threads "$THREADS" "${GTFOPT[@]}" --allele-vcf --quiet $call_extra || return 1
  # 4) describe (variant-restricted markers)
  "$BIN" describe -i "$cgfa" --bubble-prefix-in "$cpfx" --out-dir "$d/describe" \
    --variant-nodes "$d/call/call.variant_nodes.tsv" --no-wide-matrix --threads "$THREADS" --quiet || return 1
  # 4b) benchmark: round-trip QV of the calls - if variant nodes were written
  if [[ -s "$d/call/call.variant_nodes.tsv" ]]; then
    "$BIN" benchmark -i "$cgfa" --bubble-prefix-in "$cpfx" --reference-path "$ref" \
      --variant-nodes "$d/call/call.variant_nodes.tsv" --vcf "$d/call/call.region.vcf" \
      -o "$d/benchmark/benchmark" --threads "$THREADS" --quiet || return 1
    # The SAME benchmark against the allele VCF. These answer different questions and the project
    # quotes both, so regenerating only the first left the headline claim -- the allele VCF is a
    # lossless serialization, 0 bp residual -- resting on numbers nothing in the tree reproduced.
    # The region VCF measures the interpreted output; this measures the representation.
    if [[ -s "$d/call/call.alleles.vcf" ]]; then
      "$BIN" benchmark -i "$cgfa" --bubble-prefix-in "$cpfx" --reference-path "$ref" \
        --variant-nodes "$d/call/call.variant_nodes.tsv" --vcf "$d/call/call.alleles.vcf" \
        -o "$d/benchmark/allele" --threads "$THREADS" --quiet || return 1
    fi
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

REGIONS=("$@"); [[ ${#REGIONS[@]} -eq 0 ]] && REGIONS=(lpa c4 gstm1 cyp2d6 acot ankrd36c synthetic)
mkdir -p "$REPO/results"; rm -f "$CN_TABLE"

for region in "${REGIONS[@]}"; do
  case "$region" in
    lpa)
      run_region lpa "$DATA/lpa.gfa.gz" "--min-similarity 0.95" panphorte "--cn" 1
      validate_cn lpa --mode repeats --truth "$DATA/lpa.repeats.tsv" --name-col hapl --count-col copies
      # LPA GWAS demo (needs python numpy/scipy); reuses the dedicated driver
      # reuse the lpa stages already built above (no rebuild into gwas/): call ran on the refined graph
      lgfa="$OUT/lpa/refine/refine.normalized.sorted.gfa"; lpfx="$OUT/lpa/refine/refine"
      [[ -s "$lgfa" ]] || { lgfa="$OUT/lpa/panphorte/panphorte.normalized.sorted.gfa"; lpfx="$OUT/lpa/panphorte/panphorte"; }
      REUSE_GRAPH="$lgfa" REUSE_BUBBLE_PREFIX="$lpfx" REUSE_CALL_PREFIX="$OUT/lpa/call/call" \
        REUSE_COPIES="$OUT/lpa/panphorte/panphorte.panphorte.copies.tsv" \
        bash "$HERE/../tests/gwas/run_lpa_real.sh" "$BIN" "$OUT/lpa/gwas" "$PY" "$RS" || echo "  (gwas demo skipped)"
      ;;
    c4)
      run_region c4 "$DATA/c4.gfa.gz" "--min-similarity 0.95" panphorte "--cn"
      validate_cn c4 --mode gene-count --truth "$DATA/c4.bed" --genes C4A,C4B
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
      run_region gstm1 "$DATA/gstm1.gfa.gz" "--min-similarity 0.95" panphorte "--cn"
      # TWO claims, and they need two comparisons. The region VCF's DUP carries the TOTAL copy number of
      # the collapsed module -- the paralogs together -- so it is scored against the summed gene count,
      # as c4/cyp2d6/acot already are. Scoring it against GSTM1 alone (the old --mode direct check) read
      # 0% and looked like a failure, when it was only asking the module for one gene's dosage. The
      # per-gene split below is what answers that, from the sidecar.
      validate_cn gstm1 --mode gene-count --truth "$DATA/gstm.bed" --genes GSTM1,GSTM2,GSTM4,GSTM5
      # Kept as a record of the distinction: GSTM1 alone against the module CN, which it is not.
      echo "---- gstm1 GSTM1-only vs the module CN (expected to disagree; see the per-gene split) ----"
      "$PY" "$HERE/compare_copy_number.py" --mode direct --truth "$DATA/gstm1.bed" \
        --vcf "$OUT/gstm1/call/call.region.vcf" --label gstm1_gene_only || true
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
      run_region cyp2d6 "$DATA/cyp2d6.gfa.gz" "--min-similarity 0.95" panphorte "--cn"
      validate_cn cyp2d6 --mode gene-count --truth "$DATA/cyp2d6.bed" --genes CYP2D6,CYP2D7
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
      run_region acot "$DATA/acot.gfa.gz" "--min-similarity 0.95" panphorte "--cn"
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
    ankrd36c)
      run_region ankrd36c "$DATA/ankrd36c.gfa.gz" "--min-similarity 0.95" panphorte "--cn"
      ;;
    synthetic)
      echo "############ synthetic ############"
      mkdir -p "$REPO/results/synthetic_data"
      bash "$REPO/tests/synthetic_smoke.sh" "$BIN" "$REPO/tests/synthetic_data" \
        "$REPO/results/synthetic_data" || echo "  (synthetic smoke failed)"
      ;;
    *) echo "unknown region: $region" ;;
  esac
  [[ $? -eq 0 ]] || die_or_warn "region $region did not complete"
done

# A locus-level reconstruction summary, so the walkthrough and the module pages can be checked against
# a file rather than against a recollection. Two levels, never one: `region_vcf` is what a consumer of
# the compact interpreted output actually reconstructs, and `allele` is the serialization ceiling,
# which reaches 0 residual because it spells every allele out. Neither involves reads.
RECON="$REPO/results/reconstruction.tsv"
printf 'locus\tvcf\tbaseline_bp\tgraph_bp\tcalled_bp\tcarrier_bp\tregion_vcf_bp\trecovered_pct\n' > "$RECON"
for region in "${REGIONS[@]}"; do
  for kind in benchmark allele; do
    f="$OUT/$region/benchmark/$kind.qv.tsv"
    [[ -s "$f" ]] || continue
    awk -F'\t' -v L="$region" -v K="$([[ $kind == benchmark ]] && echo region || echo allele)" '
      NR==1 { for (i=1;i<=NF;i++) c[$i]=i; next }
      { r+=$(c["ref_delta"]); g+=$(c["delta"]); cl+=$(c["called_delta"]); ca+=$(c["carrier_delta"]); gt+=$(c["gt_delta"]) }
      END { printf "%s\t%s\t%d\t%d\t%d\t%d\t%d\t%.4f\n", L, K, r, g, cl, ca, gt, (r>0 ? 100*(1-gt/r) : 0) }' "$f" >> "$RECON"
  done
done
echo "wrote $RECON"

# Copy-number concordance plot (called-vs-truth, faceted by gene) from the combined table.
if [[ -s "$CN_TABLE" ]] && have_r; then
  "$RS" "$HERE/plot_cn_correlation.R" --table "$CN_TABLE" --out "$REPO/results/cn_correlation" \
    || echo "  (cn_correlation plot skipped)"
elif [[ -s "$CN_TABLE" ]]; then
  echo "  (R not found; skipping cn_correlation plot; table at $CN_TABLE)"
fi

# Round-trip QV summary: combine each locus's per-haplotype benchmark, then the reconstruction-anatomy plot.
QV_TABLE="$REPO/results/benchmark_qv.tsv"
# Columns are looked up BY NAME. They have moved once already -- the truth-event ledger was inserted
# ahead of the residual columns -- and a positional read of a table that gained a column silently
# aggregates the wrong quantity rather than failing.
QV_COLS='sample sum_delta sum_aln_len qv identity truth_missed_bp truth_below_bp called_sum_delta carrier_sum_delta gt_sum_delta gt_sum_aln_len ref_sum_delta'
printf 'locus\t%s\n' "$(printf '%s\n' $QV_COLS | paste -sd$'\t' -)" > "$QV_TABLE"
for region in "${REGIONS[@]}"; do
  bh="$OUT/$region/benchmark/benchmark.qv_by_haplotype.tsv"
  [[ -s "$bh" ]] || continue
  awk -F'\t' -v L="$region" -v WANT="$QV_COLS" '
    NR==1 { for (i = 1; i <= NF; i++) at[$i] = i
            n = split(WANT, w, " ")
            for (k = 1; k <= n; k++) if (!(w[k] in at)) { print "regen_results: benchmark table has no column " w[k] > "/dev/stderr"; exit 1 }
            next }
    { line = L; for (k = 1; k <= n; k++) line = line "\t" $at[w[k]]; print line }' "$bh" >> "$QV_TABLE"
done
# The loss partition, per locus: five consecutive terms that sum EXACTLY to the region-VCF residual.
# It lives in each locus's qv_summary.tsv rather than the per-haplotype table, so it is collected here
# for the plot. It is the only view that says WHY reconstruction falls short rather than by how much.
LOSS_TABLE="$REPO/results/benchmark_loss.tsv"
printf 'locus\tbaseline_bp\tout_of_scope\tdiscovery_or_attribution\tcarrier_missed\trepresentation\tfalse_positive_damage\n' > "$LOSS_TABLE"
for region in "${REGIONS[@]}"; do
  qs="$OUT/$region/benchmark/benchmark.qv_summary.tsv"
  [[ -s "$qs" ]] || continue
  awk -F'\t' -v L="$region" '
    $1=="gt_gap"  && $3=="baseline_delta"           { base=$4 }
    $1=="loss_bp" && $3=="out_of_scope"             { a=$4 }
    $1=="loss_bp" && $3=="discovery_or_attribution" { b=$4 }
    $1=="loss_bp" && $3=="carrier_missed"           { c=$4 }
    $1=="loss_bp" && $3=="representation"           { d=$4 }
    $1=="loss_bp" && $3=="false_positive_damage"    { e=$4 }
    END { if (base != "") printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n", L, base, a+0, b+0, c+0, d+0, e+0 }' \
    "$qs" >> "$LOSS_TABLE"
done
echo "wrote $LOSS_TABLE"

if [[ $(wc -l < "$QV_TABLE") -gt 1 ]] && have_r; then
  "$RS" "$HERE/plot_benchmark.R" --table "$QV_TABLE" --loss "$LOSS_TABLE" \
    --out "$REPO/results/benchmark_qv" --per-row 12 || die_or_warn "benchmark plot"
fi
write_manifest
if [[ "$FAILURES" -gt 0 ]]; then
  echo "regen_results: completed with $FAILURES failed stage(s) (STRICT=0). This tree is NOT release-grade." >&2
  exit 1
fi
echo "ALL DONE -> $REPO/results"
