#!/usr/bin/env bash
# Compare associate with GEMMA, linear (-lm 1) and mixed (-lmm 1, same kinship), on families built from the
# associate example cohorts (make_families.py); prints how far effects, standard errors and p-values differ.
#   dev/gemma/validate_gemma.sh [panvar binary] [output dir] [gemma binary] [associate data dir]
set -euo pipefail
here=$PWD
cd "$(dirname "$0")/../.."
arg() { if [ -z "$1" ]; then echo "$2"; elif [ "${1#/}" != "$1" ]; then echo "$1"; else echo "$here/$1"; fi; }
PANVAR=$(arg "${1:-}" build/panvar)
OUT=$(arg "${2:-}" tests/examples/gemma)
GEMMA=$(arg "${3:-}" "$(command -v gemma || true)")
EXAMPLES=$(arg "${4:-}" tests/synthetic_data/associate)
[ -x "$GEMMA" ] || { echo "no gemma: give its path (bioconda gemma)"; exit 1; }
command -v python3 > /dev/null || { echo "python3 not found on PATH"; exit 1; }

for locus in gstm lpa; do
  [ -d "$EXAMPLES/$locus/haplotype" ] || { echo "no $EXAMPLES/$locus: run dev/associate_inputs.sh first"; exit 1; }
  d=$OUT/$locus
  mkdir -p "$d"
  echo "== $locus: $(python3 dev/gemma/make_families.py "$EXAMPLES/$locus" "$d")"
  for m in lm lmm; do
    k=(); [ $m = lmm ] && k=(--model lmm --kinship "$d/kinship.tsv")
    "$PANVAR" associate --genotypes "$d/geno.bimbam.gz" --samples "$d/samples.txt" --phenotype "$d/phenotype.tsv" \
      --min-maf 0 ${k[@]+"${k[@]}"} -o "$d/panvar_$m" -q
    k=(); [ $m = lmm ] && k=(-k "$d/kinship.tsv")
    "$GEMMA" -g "$d/geno.bimbam.gz" -p "$d/gemma.pheno.txt" -c "$d/gemma.covar.txt" ${k[@]+"${k[@]}"} -$m 1 -maf 0 -miss 1 \
      -outdir "$d" -o gemma_$m > "$d/gemma_$m.out" 2>&1 || { echo "gemma -$m failed: $d/gemma_$m.out"; exit 1; }
  done
  "$GEMMA" -g "$d/geno.raw.bimbam.gz" -p "$d/gemma.pheno.txt" -c "$d/gemma.covar.txt" -lm 1 -maf 0 -miss 1 \
    -outdir "$d" -o gemma_raw > "$d/gemma_raw.out" 2>&1 || { echo "gemma on the raw matrix failed"; exit 1; }

  python3 - "$d" <<'PY'
import gzip, math, sys
d = sys.argv[1]
def table(path, key, cols):
    with open(path) as fh:
        h = next(fh).rstrip("\n").split("\t")
        out = {}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            try:
                out[f[h.index(key)]] = tuple(float(f[h.index(c)]) for c in cols)
            except ValueError:
                pass
        return out
def nlp(p):
    return -math.log10(max(p, 1e-300))
def pearson(x, y):
    mx, my = sum(x) / len(x), sum(y) / len(y)
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    return sxy / math.sqrt(sum((a - mx) ** 2 for a in x) * sum((b - my) ** 2 for b in y))
rows, na = 0, set()
for line in gzip.open(f"{d}/geno.bimbam.gz", "rt"):
    f = [x.strip() for x in line.split(",")]
    rows += 1
    if "NA" in f[3:]:
        na.add(f[0])
for m in ("lm", "lmm"):
    pv = table(f"{d}/panvar_{m}.assoc.tsv", "feature_id", ("beta", "se", "p"))
    gm = table(f"{d}/gemma_{m}.assoc.txt", "rs", ("beta", "se", "p_wald"))
    both = [k for k in pv if k in gm]
    # rows with a missing dosage are left out: GEMMA fills it with the row's mean, panvar's linear model drops the sample
    keys = [k for k in both if k not in na]
    dl = {k: abs(nlp(pv[k][2]) - nlp(gm[k][2])) for k in keys}
    worst = max(keys, key=dl.get)
    top_pv, top_gm = min(pv, key=lambda k: pv[k][2]), min(gm, key=lambda k: gm[k][2])
    print(f"  {m:4s} {len(both)} of {rows} rows in both ({len(pv) - len(both)} panvar only, {len(gm) - len(both)} gemma only); "
          f"{len(keys)} without missing dosages: beta differs by up to {max(abs(pv[k][0] - gm[k][0]) / gm[k][1] for k in keys):.1e} se, "
          f"se by {max(abs(pv[k][1] - gm[k][1]) / gm[k][1] for k in keys):.1e} (relative), -log10 p by {dl[worst]:.2g} "
          f"(at p {gm[worst][2]:.2g}), correlation of -log10 p {pearson([nlp(pv[k][2]) for k in keys], [nlp(gm[k][2]) for k in keys]):.5f}; "
          f"top row {'the same' if top_pv == top_gm else top_pv + ' vs ' + top_gm} (p {pv[top_pv][2]:.2g} vs {gm[top_gm][2]:.2g})")
h2 = next(l.split("\t")[1].strip() for l in open(f"{d}/panvar_lmm.summary.tsv") if l.startswith("lmm_h2\t"))
pve = next(l.split("=")[1].strip() for l in open(f"{d}/gemma_lmm.log.txt") if "pve estimate in the null model" in l)
raw = sum(1 for _ in open(f"{d}/gemma_raw.assoc.txt")) - 1
print(f"  share of variance the kinship explains: panvar {h2}, gemma {pve}")
print(f"  {len(na)} rows with missing dosages, not compared (GEMMA fills them with the mean, as panvar's LMM does; panvar's linear model leaves the sample out)")
print(f"  unscaled matrix: gemma tests {raw} of {rows} rows (it skips rows with a dosage above 2)")
PY
done
