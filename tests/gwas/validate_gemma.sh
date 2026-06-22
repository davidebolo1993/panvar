#!/usr/bin/env bash
# Validate `panvar associate` against GEMMA on identical BIMBAM inputs.
#
# We feed the SAME synthetic genome-wide-like panel (make_lpa_phenotype.py --sim-markers: real KIV-2 dosage
# + subpop-stratified null SNPs) and the SAME phenotype/covariates to both engines and compare the genotype
# effect (beta) and significance (-log10 p) per feature:
#   * linear model : panvar `associate --model linear` (PC covariates)  vs  GEMMA `-lm 4` (Wald)
#   * mixed model  : panvar `associate --model lmm --kinship K`         vs  GEMMA `-lmm 4 -k K`
# Concordance is the Pearson r of beta and of -log10 p over the shared features, plus top-hit agreement.
# BIMBAM is GEMMA's native mean-genotype format, so no reformatting of the genotypes is needed.
#
#   validate_gemma.sh <panvar_bin> <out_dir> [python] [gemma_bin] [copies.tsv]
set -uo pipefail

BIN="${1:?usage: validate_gemma.sh <panvar_bin> <out_dir> [python] [gemma_bin] [copies.tsv]}"
OUT="${2:?need out_dir}"
PY="${3:-python3}"
GEMMA="${4:-gemma}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
COPIES="${5:-$REPO/results/real_data/lpa/panphorte/panphorte.panphorte.copies.tsv}"
N="${N:-800}"; SIM="${SIM:-2000}"
mkdir -p "$OUT"

if ! command -v "$GEMMA" >/dev/null 2>&1 && [ ! -x "$GEMMA" ]; then
  echo "GEMMA not available ($GEMMA) -- skipping live comparison."
  echo "panvar's BIMBAM export is GEMMA's native mean-genotype format, so it loads unchanged; install"
  echo "gemma (bioconda) to run this check."
  exit 0
fi
if [ ! -f "$COPIES" ]; then
  echo "copies.tsv not found ($COPIES); run the LPA pipeline (scripts/genes/lpa.sh) first." ; exit 0
fi

echo "== build a structured cohort + synthetic panel (n=$N, sim=$SIM) =="
"$PY" "$HERE/make_lpa_phenotype.py" "$COPIES" "$OUT" --n "$N" --sim-markers "$SIM" --kinship-out "$OUT/kinship.tsv"

GENO_GZ="$OUT/geno.sim.bimbam.gz"; GENO="$OUT/geno.sim.bimbam"
gzcat "$GENO_GZ" > "$GENO"                                   # GEMMA wants a plain BIMBAM
SAMP="$OUT/sim.samples.txt"; ANNOT="$OUT/feature_annot.sim.tsv.gz"

# GEMMA phenotype (one value/line, sample order, NA kept) + covariates (intercept + Age Sex PC1-3, no NA).
# pheno.quant.tsv rows are already in sim.samples order (make_lpa writes them together).
"$PY" - "$OUT/pheno.quant.tsv" "$OUT/gemma.pheno.txt" "$OUT/gemma.covar.txt" <<'PY'
import sys
rows = [l.rstrip("\n").split("\t") for l in open(sys.argv[1])][1:]  # sample phenotype Age Sex PC1 PC2 PC3
with open(sys.argv[2], "w") as ph, open(sys.argv[3], "w") as cv:
    for r in rows:
        ph.write((r[1] if r[1] not in ("", "NA") else "NA") + "\n")
        cv.write("1\t" + "\t".join(r[2:7]) + "\n")        # intercept + covariates
PY

run_assoc() { # <model-args...> -o <prefix>
  "$BIN" associate --genotypes "$GENO_GZ" --samples "$SAMP" --feature-annot "$ANNOT" --min-maf 0 "$@"
}

echo "== LINEAR: panvar associate (PC covariates) vs GEMMA -lm 4 =="
run_assoc --phenotype "$OUT/pheno.quant.tsv" --model linear -o "$OUT/pv_lm" -q
"$GEMMA" -g "$GENO" -p "$OUT/gemma.pheno.txt" -c "$OUT/gemma.covar.txt" -lm 4 -maf 0 -miss 1 \
  -outdir "$OUT" -o gemma_lm >/dev/null 2>&1 || echo "  (gemma -lm failed)"

echo "== MIXED: panvar associate --model lmm --kinship vs GEMMA -lmm 4 -k (our GRM) =="
run_assoc --phenotype "$OUT/pheno.quant.nopc.tsv" --model lmm --kinship "$OUT/kinship.tsv" -o "$OUT/pv_lmm" -q
"$GEMMA" -g "$GENO" -p "$OUT/gemma.pheno.txt" -k "$OUT/kinship.tsv" -lmm 4 -maf 0 -miss 1 \
  -outdir "$OUT" -o gemma_lmm >/dev/null 2>&1 || echo "  (gemma -lmm failed)"

echo "== concordance =="
"$PY" - "$OUT" <<'PY'
import sys, math
out = sys.argv[1]
def read_pv(p):  # panvar assoc.tsv -> {feature: (beta, p)}
    d = {}
    L = open(p).read().splitlines()
    h = L[0].split("\t"); bi = h.index("beta") if "beta" in h else 6; pi = h.index("p")
    for l in L[1:]:
        f = l.split("\t")
        try: d[f[0]] = (float(f[bi]), float(f[pi]))
        except: pass
    return d
def read_gemma(p):  # gemma .assoc.txt -> {rs: (beta, p_wald)}
    d = {}
    L = open(p).read().splitlines()
    h = L[0].split("\t"); ri=h.index("rs"); bi=h.index("beta"); pi=h.index("p_wald")
    for l in L[1:]:
        f = l.split("\t")
        try: d[f[ri]] = (float(f[bi]), float(f[pi]))
        except: pass
    return d
def pearson(xs, ys):
    n=len(xs);
    if n<3: return float("nan")
    mx=sum(xs)/n; my=sum(ys)/n
    sxy=sum((x-mx)*(y-my) for x,y in zip(xs,ys))
    sx=math.sqrt(sum((x-mx)**2 for x in xs)); sy=math.sqrt(sum((y-my)**2 for y in ys))
    return sxy/(sx*sy) if sx>0 and sy>0 else float("nan")
def compare(tag, pv_path, gm_path):
    import os
    if not os.path.exists(gm_path): print(f"  {tag}: GEMMA output missing -> skip"); return
    pv=read_pv(pv_path); gm=read_gemma(gm_path)
    keys=[k for k in pv if k in gm]
    nlp=lambda p: -math.log10(max(p,1e-300))
    rb=pearson([pv[k][0] for k in keys],[gm[k][0] for k in keys])
    rp=pearson([nlp(pv[k][1]) for k in keys],[nlp(gm[k][1]) for k in keys])
    pv_top=min(pv,key=lambda k:pv[k][1]); gm_top=min(gm,key=lambda k:gm[k][1])
    print(f"  {tag}: features={len(keys)}  r(beta)={rb:.4f}  r(-log10p)={rp:.4f}  "
          f"top: panvar={pv_top} gemma={gm_top} {'MATCH' if pv_top==gm_top else 'differ'}")
compare("linear", f"{out}/pv_lm.assoc.tsv",  f"{out}/gemma_lm.assoc.txt")
compare("lmm",    f"{out}/pv_lmm.assoc.tsv", f"{out}/gemma_lmm.assoc.txt")
PY
echo "== DONE: $OUT =="
