#!/usr/bin/env bash
# LPA association demo on the real graph: bubble -> panphorte -> call -> describe --samples -> associate.
# Produces (1) a region scan that recovers KIV-2 as the top hit, and (2) a structure-correction demo on a
# synthetic genome-wide panel (naive lambda>>1 -> ~1 with PCs/LMM). See docs/gwas/example.md. Needs a
# numpy-capable python.
#   run_lpa_real.sh <panvar_bin> <out_dir> [python] [Rscript] [--big]
# Env: N (cohort size; --big sets 10000), SIM (null markers, default 3000), LMM_MAX_N (LMM size cap).
set -euo pipefail

PANVAR_BIN="${1:?usage: run_lpa_real.sh <panvar_bin> <out_dir> [python] [Rscript] [--big]}"
OUT_DIR="${2:?need out_dir}"
PY="${3:-python3}"
RS="${4:-Rscript}"
[ "${5:-}" = "--big" ] && N="${N:-10000}"
N="${N:-300}"
SIM="${SIM:-3000}"
LMM_MAX_N="${LMM_MAX_N:-4000}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
GFA="$REPO/tests/real_data/lpa.gfa.gz"
REAL="$OUT_DIR/real"
mkdir -p "$OUT_DIR" "$REAL"
REF="$(gzcat "$GFA" | awk -F'\t' '($1=="P"||$1=="W"){n=$2; if(g==""&&n~/[Gg][Rr][Cc]h38/)g=n; if(f=="")f=n} END{print (g!=""?g:f)}')"
echo "reference path: $REF ; cohort N=$N ; sim null markers=$SIM"

echo "== bubble -> panphorte (KIV-2 copies) =="
"$PANVAR_BIN" bubble    -i "$GFA" -o "$OUT_DIR/bub" -r "$REF" --quiet >/dev/null
SGFA="$OUT_DIR/bub.sorted.gfa"
"$PANVAR_BIN" panphorte -i "$SGFA" --bubble-prefix-in "$OUT_DIR/bub" -o "$OUT_DIR/pan" \
  --reference-path "$REF" --min-similarity 0.97 --quiet >/dev/null
PGFA="$OUT_DIR/pan.normalized.sorted.gfa"

# structured cohort: cosigt samples.tsv + pheno.{quant,binary}{,.nopc}.tsv (subpops, PCs, NA) +
# synthetic structure-demo panel; kinship GRM only when N is small enough for a dense matrix file.
KOPT=(); [ "$N" -le "$LMM_MAX_N" ] && KOPT=(--kinship-out "$REAL/kinship.tsv")
"$PY" "$HERE/make_lpa_phenotype.py" "$OUT_DIR/pan.panphorte.copies.tsv" "$REAL" \
  --n "$N" --sim-markers "$SIM" "${KOPT[@]}"

echo "== call -> describe --samples (per-sample BIMBAM dosage) =="
GTF="${GTF:-$REPO/tests/real_data/Homo_sapiens.GRCh38.116.gtf.gz}"
GTFOPT=(); [ -f "$GTF" ] && GTFOPT=(--gtf "$GTF")
"$PANVAR_BIN" call -i "$PGFA" --bubble-prefix-in "$OUT_DIR/pan" --reference-path "$REF" \
  -o "$OUT_DIR/call" --cn-from-multiplicity ${GTFOPT[@]+"${GTFOPT[@]}"} --quiet >/dev/null
NGOPT=(); [ -f "$OUT_DIR/call.node_genes.tsv" ] && NGOPT=(--node-genes "$OUT_DIR/call.node_genes.tsv")
"$PANVAR_BIN" describe -i "$PGFA" --bubble-prefix-in "$OUT_DIR/pan" --out-dir "$OUT_DIR/desc" \
  --kmer-size 31 --no-wide-matrix --variant-nodes "$OUT_DIR/call.variant_nodes.tsv" \
  --samples "$REAL/samples.tsv" --quiet >/dev/null

DESC="$OUT_DIR/desc"
SAMP="$DESC/bimbam.samples.samples.txt.gz"
ANNOT="$DESC/feature_annot.samples.tsv.gz"
plot() {  # plot <assoc_prefix> <title>
  "$RS" "$REPO/scripts/plot_associate.R" --assoc "$1.assoc.tsv" --summary "$1.summary.tsv" \
    --out "$1" --title "$2" >/dev/null 2>&1 || echo "  (plot skipped: ggplot2?)"
}

echo "== 1) REGION SCAN: associate on the real KIV-2 features (PC-adjusted) =="
# graph = node/edge dosage (carries the KIV-2 self-loop REP node); kmer = k-mer multiplicity.
# quant = log10 Lp(a) (linear); binary = high-risk case/control (logistic). pheno tables carry PCs.
for sub in graph kmers; do
  GENO="$DESC/bimbam_${sub}.samples.bimbam.gz"
  for mode in quant binary; do
    echo "   -- associate ($sub / $mode) --"
    "$PANVAR_BIN" associate --genotypes "$GENO" --samples "$SAMP" --feature-annot "$ANNOT" \
      "${NGOPT[@]}" --phenotype "$REAL/pheno.${mode}.tsv" --min-maf 0.02 -o "$OUT_DIR/assoc_${sub}_${mode}"
    plot "$OUT_DIR/assoc_${sub}_${mode}" "LPA Lp(a) ($sub, $mode)"
  done
done

echo "== 2) STRUCTURE-CORRECTION demo on the synthetic genome-wide-like panel =="
SGENO="$REAL/geno.sim.bimbam.gz"; SSAMP="$REAL/sim.samples.txt"; SANNOT="$REAL/feature_annot.sim.tsv.gz"
if [ -f "$SGENO" ]; then
  echo "   -- naive (Age+Sex, no PCs): expect inflated lambda --"
  "$PANVAR_BIN" associate --genotypes "$SGENO" --samples "$SSAMP" --feature-annot "$SANNOT" \
    --phenotype "$REAL/pheno.quant.nopc.tsv" --min-maf 0.02 -o "$OUT_DIR/sim_naive"
  plot "$OUT_DIR/sim_naive" "LPA structure demo (naive)"
  echo "   -- + ancestry PC covariates: lambda back toward 1 --"
  "$PANVAR_BIN" associate --genotypes "$SGENO" --samples "$SSAMP" --feature-annot "$SANNOT" \
    --phenotype "$REAL/pheno.quant.tsv" --min-maf 0.02 -o "$OUT_DIR/sim_pc"
  plot "$OUT_DIR/sim_pc" "LPA structure demo (PC-adjusted)"
  if [ "$N" -le "$LMM_MAX_N" ]; then
    echo "   -- LMM with panel-derived GRM (--make-kinship): lambda ~ 1 --"
    "$PANVAR_BIN" associate --genotypes "$SGENO" --samples "$SSAMP" --feature-annot "$SANNOT" \
      --phenotype "$REAL/pheno.quant.nopc.tsv" --model lmm --make-kinship --min-maf 0.02 -o "$OUT_DIR/sim_lmm"
    plot "$OUT_DIR/sim_lmm" "LPA structure demo (LMM)"
  else
    echo "   (LMM skipped: N=$N > LMM_MAX_N=$LMM_MAX_N; use --pca / PC covariates at this scale)"
  fi
fi

echo "== sanity: region scan recovers KIV-2 (bubble 7, negative effect); structure correction lowers lambda =="
"$PY" - "$OUT_DIR/assoc_graph_quant.assoc.tsv" "$OUT_DIR/sim_naive.summary.tsv" "$OUT_DIR/sim_pc.summary.tsv" <<'PY'
import sys
rows = [l.split("\t") for l in open(sys.argv[1]).read().splitlines()[1:]]
top = rows[0]
beta, p_bonf = float(top[6]), float(top[10])
ok = ("7" in top[2].split(";")) and beta < 0 and p_bonf < 0.05
print(f"  region top feature={top[0]} bubble={top[2]} beta={beta:.3f} p_bonf={p_bonf:.1e} -> {'PASS' if ok else 'CHECK'}")
def lam(p):
    try:
        d = dict(l.split("\t") for l in open(p).read().splitlines()[1:]); return float(d["lambda_gc"])
    except Exception:
        return float("nan")
ln, lp = lam(sys.argv[2]), lam(sys.argv[3])
print(f"  structure-demo lambda: naive={ln:.2f} -> PC-adjusted={lp:.2f} -> {'PASS' if lp < ln else 'CHECK'}")
PY
echo "== DONE: $OUT_DIR (assoc_*, sim_{naive,pc,lmm}.{assoc.tsv,summary.tsv,manhattan.png,qq.png}) =="
