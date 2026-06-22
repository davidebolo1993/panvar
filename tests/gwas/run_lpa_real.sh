#!/usr/bin/env bash
# LPA pangenome-association demo on the real LPA graph (tests/real_data/lpa.gfa.gz), end to end with
# the native modules: bubble -> panphorte -> call -> describe --samples -> ASSOCIATE.
#
# KIV-2 copy number is read per real haplotype from panphorte; real haplotypes are paired into diploid
# samples (cosigt-style) with a literature-plausible INVERSE Lp(a) phenotype (synthetic values, real
# topology) plus synthetic covariates (Age, Sex, PC1-3) and ~5% NA phenotypes. `describe --samples`
# emits per-sample BIMBAM dosage; `panvar associate` runs the GWAS (GLM + covariates, MAF filter on the
# final genotypes, region-wide Bonferroni/BH correction). The KIV-2 locus is recovered as the top hit
# with a NEGATIVE effect (more KIV-2 -> lower Lp(a)).
#
#   run_lpa_real.sh <panvar_bin> <out_dir> [python] [Rscript]
set -euo pipefail

PANVAR_BIN="${1:?usage: run_lpa_real.sh <panvar_bin> <out_dir> [python] [Rscript]}"
OUT_DIR="${2:?need out_dir}"
PY="${3:-python3}"
RS="${4:-Rscript}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
GFA="$REPO/tests/real_data/lpa.gfa.gz"
REAL="$OUT_DIR/real"
mkdir -p "$OUT_DIR" "$REAL"
REF="$(gzcat "$GFA" | awk -F'\t' '($1=="P"||$1=="W"){n=$2; if(g==""&&n~/[Gg][Rr][Cc]h38/)g=n; if(f=="")f=n} END{print (g!=""?g:f)}')"
echo "reference path: $REF"

echo "== bubble -> panphorte (KIV-2 copies) =="
"$PANVAR_BIN" bubble    -i "$GFA" -o "$OUT_DIR/bub" -r "$REF" --quiet >/dev/null
SGFA="$OUT_DIR/bub.sorted.gfa"
"$PANVAR_BIN" panphorte -i "$SGFA" --bubble-prefix-in "$OUT_DIR/bub" -o "$OUT_DIR/pan" \
  --reference-path "$REF" --min-similarity 0.90 --quiet >/dev/null
PGFA="$OUT_DIR/pan.normalized.sorted.gfa"
# synthetic cohort: cosigt samples.tsv + pheno.{quant,binary}.tsv (covariates + NA)
"$PY" "$HERE/make_lpa_phenotype.py" "$OUT_DIR/pan.panphorte.copies.tsv" "$REAL" 200

echo "== call -> describe --samples (per-sample BIMBAM dosage) =="
GTF="${GTF:-$REPO/tests/real_data/Homo_sapiens.GRCh38.116.gtf.gz}"
GTFOPT=(); [ -f "$GTF" ] && GTFOPT=(--gtf "$GTF")
"$PANVAR_BIN" call -i "$PGFA" --bubble-prefix-in "$OUT_DIR/pan" --reference-path "$REF" \
  -o "$OUT_DIR/call" --cn-from-multiplicity ${GTFOPT[@]+"${GTFOPT[@]}"} --quiet >/dev/null
"$PANVAR_BIN" describe -i "$PGFA" --bubble-prefix-in "$OUT_DIR/pan" --out-dir "$OUT_DIR/desc" \
  --kmer-size 31 --no-wide-matrix --variant-nodes "$OUT_DIR/call.variant_nodes.tsv" \
  --samples "$REAL/samples.tsv" --quiet >/dev/null

DESC="$OUT_DIR/desc"
SAMP="$DESC/bimbam.samples.samples.txt.gz"
ANNOT="$DESC/feature_annot.samples.tsv.gz"
# associate: both substrates x both trait codings. graph = node/edge dosage (carries the KIV-2 self-loop
# REP node); kmer = k-mer multiplicity. quant = Lp(a) level (linear); binary = case/control (logistic).
for sub in graph kmers; do
  GENO="$DESC/bimbam_${sub}.samples.bimbam.gz"
  for mode in quant binary; do
    echo "== associate ($sub / $mode) =="
    "$PANVAR_BIN" associate --genotypes "$GENO" --samples "$SAMP" --feature-annot "$ANNOT" \
      --phenotype "$REAL/pheno.${mode}.tsv" --min-maf 0.02 -o "$OUT_DIR/assoc_${sub}_${mode}"
    "$RS" "$REPO/scripts/plot_associate.R" --assoc "$OUT_DIR/assoc_${sub}_${mode}.assoc.tsv" \
      --summary "$OUT_DIR/assoc_${sub}_${mode}.summary.tsv" --out "$OUT_DIR/assoc_${sub}_${mode}" \
      --title "LPA Lp(a) ($sub, $mode)" >/dev/null 2>&1 || echo "  (plot skipped: ggplot2?)"
  done
done

echo "== sanity: associate recovers the KIV-2 locus (bubble 7) with a negative effect =="
"$PY" - "$OUT_DIR/assoc_graph_quant.assoc.tsv" <<'PY'
import sys
rows = [l.split("\t") for l in open(sys.argv[1]).read().splitlines()[1:]]
# cols: 0 feature 1 layer 2 bubbles 3 nodes 4 n 5 minor_freq 6 beta 7 se 8 z 9 p 10 p_bonf 11 q_bh
top = rows[0]
beta, p_bonf = float(top[6]), float(top[10])
ok = ("7" in top[2].split(";")) and beta < 0 and p_bonf < 0.05
print(f"  top feature={top[0]} bubble={top[2]} beta={beta:.3f} p_bonf={p_bonf:.1e}  -> {'PASS' if ok else 'CHECK'}")
PY
echo "== DONE: $OUT_DIR (assoc_{graph,kmer}_{quant,binary}.{assoc.tsv,summary.tsv,manhattan.png,qq.png}) =="
