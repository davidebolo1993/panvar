#!/usr/bin/env bash
# Check associate against the tools its methods follow, on data from simulate.R: R's lm and glm, SPAtest, SKAT,
# susieR and GEMMA, and R versions of the exact counts, the conditional selection and EMMAX. reference.R prints
# one line per comparison (also <out dir>/comparison.tsv) and fails if any is outside its tolerance.
# --calibrate also runs calibrate.sh: null simulations for what no tool computes.
#   dev/validate_associate/validate.sh [--calibrate] [panvar] [out dir] [Rscript] [gemma]
# Needs R with SKAT, susieR and SPAtest, and GEMMA.
set -euo pipefail
here=$PWD
cd "$(dirname "$0")/../.."
calibrate=0
if [ "${1:-}" = "--calibrate" ]; then calibrate=1; shift; fi
arg() { if [ -z "$1" ]; then echo "$2"; elif [ "${1#/}" != "$1" ]; then echo "$1"; else echo "$here/$1"; fi; }
PANVAR=$(arg "${1:-}" build/panvar)
OUT=$(arg "${2:-}" tests/examples/validate_associate)
RSCRIPT=$(arg "${3:-}" "$(command -v Rscript || true)")
GEMMA=$(arg "${4:-}" "$(command -v gemma || true)")
[ -x "$PANVAR" ] || { echo "no panvar binary at $PANVAR"; exit 1; }
[ -x "$RSCRIPT" ] || { echo "no Rscript: give its path (R with SKAT, susieR, SPAtest)"; exit 1; }
[ -x "$GEMMA" ] || { echo "no gemma: give its path (bioconda gemma)"; exit 1; }
mkdir -p "$OUT"

"$RSCRIPT" dev/validate_associate/simulate.R "$OUT"
U=$OUT/unrelated F=$OUT/families
run() { local log=$1; shift; "$@" > "$log" 2>&1 || { echo "failed: $* (log: $log)"; tail -n 3 "$log"; exit 1; }; }
common=(--genotypes "$U/g.bimbam" --samples "$U/samples.txt" --feature-annot "$U/annot.tsv" --unit variant --min-maf 0
        --rare-maf 0.025 --group gene --cojo-p 1e-4)
run "$U/q.log"  "$PANVAR" associate "${common[@]}" --phenotype "$U/pheno.q.tsv" -o "$U/q"
run "$U/qr.log" "$PANVAR" associate "${common[@]}" --phenotype "$U/pheno.q.tsv" --rint -o "$U/qr"
run "$U/b.log"  "$PANVAR" associate "${common[@]}" --phenotype "$U/pheno.b.tsv" -o "$U/b"
run "$U/b0.log" "$PANVAR" associate "${common[@]}" --phenotype "$U/pheno.b0.tsv" -o "$U/b0"
fam=(--genotypes "$F/g.bimbam" --samples "$F/samples.txt" --phenotype "$F/pheno.tsv" --min-maf 0 --model lmm
     --kinship "$F/kinship.txt" --finemap-l 0)
run "$F/k.log"  "$PANVAR" associate "${fam[@]}" -o "$F/k"
run "$F/kf.log" "$PANVAR" associate "${fam[@]}" --lmm-fast -o "$F/kf"
awk -F'\t' 'NR > 1 { print $2 }' "$F/pheno.tsv" > "$F/gemma.pheno.txt"
awk -F'\t' 'NR > 1 { print 1, $3 }' "$F/pheno.tsv" > "$F/gemma.covar.txt"
run "$F/gemma.log" "$GEMMA" -g "$F/g.bimbam" -p "$F/gemma.pheno.txt" -c "$F/gemma.covar.txt" -k "$F/kinship.txt" \
  -lmm 1 -maf 0 -outdir "$F" -o gemma

status=0
"$RSCRIPT" dev/validate_associate/reference.R "$OUT" || status=$?
if [ $calibrate = 1 ]; then bash dev/validate_associate/calibrate.sh "$PANVAR" "$OUT" || status=1; fi
exit $status
