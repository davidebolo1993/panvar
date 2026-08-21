#!/usr/bin/env bash
# LPA association demo on the real graph: (reuse or build) -> describe --samples -> associate.
# Produces (1) a region scan that recovers KIV-2 as the top hit, and (2) a structure-correction demo on a
# synthetic structure panel (naive lambda>>1 -> ~1 with the supplied PCs; the LMM is a negative
# control here because this fixture's GRM does not capture the simulated structure). See docs/gwas.md. Needs a
# numpy-capable python.
#   run_lpa_real.sh <panvar_bin> <out_dir> [python] [Rscript] [--big]
# Env: N (cohort size; --big sets 6000, the Moli-sani WGS scale), SIM (null markers, default 3000),
#      LMM_MAX_N (LMM size cap).
# Reuse env (regen sets these so the GWAS does not rebuild the pipeline; unset = self-contained build
# under <out_dir>/prep): REUSE_GRAPH, REUSE_BUBBLE_PREFIX, REUSE_CALL_PREFIX, REUSE_COPIES.
set -euo pipefail

PANVAR_BIN="${1:?usage: run_lpa_real.sh <panvar_bin> <out_dir> [python] [Rscript] [--big]}"
OUT_DIR="${2:?need out_dir}"
PY="${3:-python3}"
RS="${4:-Rscript}"
[ "${5:-}" = "--big" ] && N="${N:-6000}"   # ~Moli-sani WGS cohort scale
N="${N:-300}"
SIM="${SIM:-3000}"
LMM_MAX_N="${LMM_MAX_N:-4000}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
GFA="$REPO/tests/real_data/lpa.gfa.gz"
REAL="$OUT_DIR/real"
ASSOC="$OUT_DIR/associate"       # all association outputs (region scan + structure demo) live here
mkdir -p "$OUT_DIR" "$REAL" "$ASSOC"
REF="$(gunzip -c "$GFA" | awk -F'\t' '($1=="P"||$1=="W"){n=$2; if(g==""&&n~/[Gg][Rr][Cc]h38/)g=n; if(f=="")f=n} END{print (g!=""?g:f)}')"
echo "reference path: $REF ; cohort N=$N ; sim null markers=$SIM"

# --- (reuse or build) the call substrate: graph, bubble prefix, call prefix, copies table -------------
if [ -n "${REUSE_GRAPH:-}" ]; then
  echo "== reusing pre-built pipeline (no rebuild): $REUSE_GRAPH =="
  PGFA="$REUSE_GRAPH"; BPFX="${REUSE_BUBBLE_PREFIX:?need REUSE_BUBBLE_PREFIX}"
  CALL="${REUSE_CALL_PREFIX:?need REUSE_CALL_PREFIX}"; COPIES="${REUSE_COPIES:?need REUSE_COPIES}"
else
  PREP="$OUT_DIR/prep"; mkdir -p "$PREP"
  echo "== bubble -> panphorte (KIV-2 copies) =="
  "$PANVAR_BIN" bubble    -i "$GFA" -o "$PREP/bub" -r "$REF" --quiet >/dev/null
  "$PANVAR_BIN" panphorte -i "$PREP/bub.sorted.gfa" --bubble-prefix-in "$PREP/bub" -o "$PREP/pan" \
    --reference-path "$REF" --min-similarity 0.97 --quiet >/dev/null
  GTF="${GTF:-$REPO/tests/real_data/Homo_sapiens.GRCh38.116.gtf.gz}"
  GTFOPT=(); [ -f "$GTF" ] && GTFOPT=(--gtf "$GTF")
  "$PANVAR_BIN" call -i "$PREP/pan.normalized.sorted.gfa" --bubble-prefix-in "$PREP/pan" \
    --reference-path "$REF" -o "$PREP/call" --cn ${GTFOPT[@]+"${GTFOPT[@]}"} --quiet >/dev/null
  PGFA="$PREP/pan.normalized.sorted.gfa"; BPFX="$PREP/pan"; CALL="$PREP/call"
  COPIES="$PREP/pan.panphorte.copies.tsv"
fi
NGOPT=(); [ -f "${CALL}.node_genes.tsv" ] && NGOPT=(--node-genes "${CALL}.node_genes.tsv")

# structured cohort: cosigt samples.tsv + pheno.{quant,binary}{,.nopc}.tsv (subpops, PCs, NA) +
# synthetic structure-demo panel; kinship GRM only when N is small enough for a dense matrix file.
KOPT=(); [ "$N" -le "$LMM_MAX_N" ] && KOPT=(--kinship-out "$REAL/kinship.tsv")
"$PY" "$HERE/make_lpa_phenotype.py" "$COPIES" "$REAL" --n "$N" --sim-markers "$SIM" ${KOPT[@]+"${KOPT[@]}"}

echo "== describe --samples (per-sample BIMBAM dosage) =="
"$PANVAR_BIN" describe -i "$PGFA" --bubble-prefix-in "$BPFX" --out-dir "$OUT_DIR/desc" \
  --kmer-size 31 --no-wide-matrix --variant-nodes "${CALL}.variant_nodes.tsv" \
  --variant-vcf "${CALL}.region.vcf" --samples "$REAL/samples.tsv" --quiet >/dev/null

DESC="$OUT_DIR/desc"   # per-substrate sample-level inputs live under $DESC/sample/<substrate>/
plot() {  # plot <assoc_prefix> <title>
  # Report what actually failed. This used to guess "(plot skipped: ggplot2?)", which sent me looking
  # for a missing R package that was installed -- the real cause was $RS not being on PATH.
  local err
  err="$("$RS" "$REPO/scripts/plot_associate.R" --assoc "$1.assoc.tsv" --summary "$1.summary.tsv" \
    --out "$1" --title "$2" 2>&1 >/dev/null)" \
    || echo "  (plot skipped: ${err:-$RS failed with no message}$([ -x "$(command -v "$RS")" ] || echo " -- '$RS' is not executable; set RSCRIPT=/path/to/Rscript)")"
}
# pipeline <real_prefix> <unfiltered_prefix> <min_maf> <title>: faceted per-stage Manhattan
# (TEST -> FILTER MAF -> [CLUMP] -> CORRECT -> CONDITION). Needs an extra --min-maf 0 run for the
# TEST/FILTER stages (so the MAF-dropped features are visible).
pipeline() {
  local err
  err="$("$RS" "$REPO/scripts/plot_associate_pipeline.R" --assoc "$1.assoc.tsv" --unfiltered "$2.assoc.tsv" \
    --summary "$1.summary.tsv" --min-maf "$3" --out "$1" --title "$4" 2>&1 >/dev/null)" \
    || echo "  (pipeline plot skipped: ${err:-$RS failed with no message}$([ -x "$(command -v "$RS")" ] || echo " -- '$RS' is not executable; set RSCRIPT=/path/to/Rscript)")"
}

echo "== 1) REGION SCAN: associate on the real KIV-2 features (PC-adjusted) =="
# graph = node/edge dosage (carries the KIV-2 self-loop REP node); kmer = k-mer multiplicity.
# quant = log10 Lp(a) (linear); binary = high-risk case/control (logistic). pheno tables carry PCs.
for sub in graph kmers; do
  SD="$DESC/sample/$sub"; GENO="$SD/bimbam_${sub}.bimbam.gz"
  SAMP="$SD/samples.txt.gz"; ANNOT="$SD/feature_annot.${sub}.tsv.gz"
  for mode in quant binary; do
    echo "   -- associate ($sub / $mode) --"
    "$PANVAR_BIN" associate --genotypes "$GENO" --samples "$SAMP" --feature-annot "$ANNOT" \
      ${NGOPT[@]+"${NGOPT[@]}"} --phenotype "$REAL/pheno.${mode}.tsv" --min-maf 0.02 -o "$ASSOC/assoc_${sub}_${mode}"
    plot "$ASSOC/assoc_${sub}_${mode}" "LPA Lp(a) ($sub, $mode)"
    if [ "$mode" = quant ]; then  # extra unfiltered run + per-stage pipeline plot
      "$PANVAR_BIN" associate --genotypes "$GENO" --samples "$SAMP" --feature-annot "$ANNOT" \
        ${NGOPT[@]+"${NGOPT[@]}"} --phenotype "$REAL/pheno.quant.tsv" --min-maf 0 -o "$ASSOC/assoc_${sub}_quant.unfiltered" >/dev/null
      pipeline "$ASSOC/assoc_${sub}_quant" "$ASSOC/assoc_${sub}_quant.unfiltered" 0.02 "LPA Lp(a) ($sub)"
    fi
  done
done

echo "== 1b) VARIANT-LEVEL scan: one test per SV call (Li-Ji Meff + LD labels) =="
# This coarser unit tests the SV calls directly, so n_tests = #variant records rather than the larger
# set of correlated k-mers/nodes. Records can still be correlated or encode linked parts of one event.
VD="$DESC/sample/variant"
for mode in quant binary; do
  echo "   -- associate (variant / $mode) --"
  "$PANVAR_BIN" associate --genotypes "$VD/bimbam_variant.bimbam.gz" \
    --samples "$VD/samples.txt.gz" \
    --feature-annot "$VD/feature_annot.variant.tsv.gz" \
    ${NGOPT[@]+"${NGOPT[@]}"} --phenotype "$REAL/pheno.${mode}.tsv" -o "$ASSOC/assoc_variant_${mode}"
  plot "$ASSOC/assoc_variant_${mode}" "LPA Lp(a) (variant, $mode)"
  if [ "$mode" = quant ]; then  # extra unfiltered run + per-stage pipeline plot (variant tier has CLUMP)
    "$PANVAR_BIN" associate --genotypes "$VD/bimbam_variant.bimbam.gz" \
      --samples "$VD/samples.txt.gz" --feature-annot "$VD/feature_annot.variant.tsv.gz" \
      ${NGOPT[@]+"${NGOPT[@]}"} --phenotype "$REAL/pheno.quant.tsv" --min-maf 0 -o "$ASSOC/assoc_variant_quant.unfiltered" >/dev/null
    pipeline "$ASSOC/assoc_variant_quant" "$ASSOC/assoc_variant_quant.unfiltered" 0.01 "LPA Lp(a) (variants)"
  fi
done

echo "== 2) STRUCTURE-CORRECTION demo on a synthetic panel with many null markers =="
SGENO="$REAL/geno.sim.bimbam.gz"; SSAMP="$REAL/sim.samples.txt"; SANNOT="$REAL/feature_annot.sim.tsv.gz"
if [ -f "$SGENO" ]; then
  echo "   -- naive (Age+Sex, no PCs): expect inflated lambda --"
  "$PANVAR_BIN" associate --genotypes "$SGENO" --samples "$SSAMP" --feature-annot "$SANNOT" \
    --phenotype "$REAL/pheno.quant.nopc.tsv" --min-maf 0.02 -o "$ASSOC/sim_naive"
  plot "$ASSOC/sim_naive" "LPA structure demo (naive)"
  echo "   -- + ancestry PC covariates: lambda back toward 1 --"
  "$PANVAR_BIN" associate --genotypes "$SGENO" --samples "$SSAMP" --feature-annot "$SANNOT" \
    --phenotype "$REAL/pheno.quant.tsv" --min-maf 0.02 -o "$ASSOC/sim_pc"
  plot "$ASSOC/sim_pc" "LPA structure demo (PC-adjusted)"
  if [ "$N" -le "$LMM_MAX_N" ] && [ -f "$REAL/kinship.tsv" ]; then
    echo "   -- LMM with the panel GRM (--kinship); compare its lambda_gc against sim_pc --"
    "$PANVAR_BIN" associate --genotypes "$SGENO" --samples "$SSAMP" --feature-annot "$SANNOT" \
      --phenotype "$REAL/pheno.quant.nopc.tsv" --model lmm --kinship "$REAL/kinship.tsv" --min-maf 0.02 -o "$ASSOC/sim_lmm"
    plot "$ASSOC/sim_lmm" "LPA structure demo (LMM)"
  else
    echo "   (LMM skipped: N=$N > LMM_MAX_N=$LMM_MAX_N or no kinship.tsv; use PC covariates at this scale)"
  fi
fi

# KIV-2 = the largest-DUP bubble (id is graph-dependent, so detect it rather than hardcode)
KIV2_BUB="$(awk -F'\t' '/SVTYPE=DUP/{sv=$8;sub(/.*SVLEN=/,"",sv);sub(/;.*/,"",sv);if(sv<0)sv=-sv;
  b=$8;sub(/.*BUBBLE_ID=/,"",b);sub(/;.*/,"",b); if(sv+0>m){m=sv+0;bb=b}} END{print bb}' "${CALL}.region.vcf")"
echo "== sanity: region scan recovers KIV-2 (bubble $KIV2_BUB, negative effect); structure correction lowers lambda =="
"$PY" - "$ASSOC/assoc_graph_quant.assoc.tsv" "$ASSOC/sim_naive.summary.tsv" "$ASSOC/sim_pc.summary.tsv" "$KIV2_BUB" <<'PY'
import sys
lines = open(sys.argv[1]).read().splitlines()
# BY NAME. These were positional, and `associate` has since inserted n_conditional at index 5 and
# added p_method/effect_status/mac_case/mac_ctrl -- so top[6] had become minor_freq and top[10] p.
# The check was asking whether a FREQUENCY is negative, which it never is, so it could not pass: it
# printed CHECK on a result that was in fact correct (bubble 7, beta -0.019, p_bonf 3e-08).
hdr = lines[0].split("\t")
ci = {c: i for i, c in enumerate(hdr)}
for need in ("feature_id", "bubbles", "beta", "p_bonf"):
    if need not in ci:
        sys.exit(f"{sys.argv[1]}: no '{need}' column; header is {hdr}")
top = lines[1].split("\t")
beta, p_bonf = float(top[ci["beta"]]), float(top[ci["p_bonf"]])
bubbles = top[ci["bubbles"]]
ok = (sys.argv[4] in bubbles.split(";")) and beta < 0 and p_bonf < 0.05
print(f"  region top feature={top[ci['feature_id']]} bubble={bubbles} beta={beta:.4f} "
      f"p_bonf={p_bonf:.1e} -> {'PASS' if ok else 'CHECK'}")
def lam(p):
    try:
        d = dict(l.split("\t") for l in open(p).read().splitlines()[1:]); return float(d["lambda_gc"])
    except Exception:
        return float("nan")
ln, lp = lam(sys.argv[2]), lam(sys.argv[3])
print(f"  structure-demo lambda: naive={ln:.2f} -> PC-adjusted={lp:.2f} -> {'PASS' if lp < ln else 'CHECK'}")
PY
echo "== DONE: $ASSOC (assoc_*, sim_{naive,pc,lmm}.{assoc.tsv,summary.tsv,manhattan.png,qq.png}) =="
