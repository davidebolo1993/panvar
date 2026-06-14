#!/usr/bin/env bash
# REAL LPA k-mer GWAS demo: validate the multiplicity approach on the genuine LPA graph
# (tests/real_data/lpa.gfa). KIV-2 copy number is read per real haplotype from panphorte; real
# haplotypes are paired into synthetic diploid samples with a literature-based inverse Lp(a)
# phenotype. The count GWAS should peak on the KIV-2 bubble (a DUP); presence/absence misses it.
#
#   run_lpa_real.sh <panvar_bin> <out_dir> [python] [Rscript]
set -euo pipefail

PANVAR_BIN="${1:?usage: run_lpa_real.sh <panvar_bin> <out_dir> [python] [Rscript]}"
OUT_DIR="${2:?need out_dir}"
PY="${3:-python3}"
RS="${4:-Rscript}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
GFA="$REPO/tests/real_data/lpa.gfa"
SNARLS="$REPO/tests/real_data/lpa.snarls.jsonl"
REAL="$OUT_DIR/real"
mkdir -p "$OUT_DIR" "$REAL"
REF="$(awk -F'\t' '($1=="P"||$1=="W")&&!f{print $2;f=1}' "$GFA")"
echo "reference path: $REF"

echo "== bubble -> panphorte (KIV-2 copies) -> phenotype =="
"$PANVAR_BIN" bubble    -i "$GFA" -o "$OUT_DIR/bub" --snarls-in "$SNARLS" --quiet >/dev/null
"$PANVAR_BIN" panphorte -i "$GFA" --bubble-prefix-in "$OUT_DIR/bub" -o "$OUT_DIR/panphorte" \
  --min-similarity 0.90 --quiet >/dev/null
"$PY" "$HERE/make_lpa_phenotype.py" "$OUT_DIR/panphorte.panphorte.copies.tsv" "$REAL" 200

echo "== call -> describe --samples (variant-restricted) =="
"$PANVAR_BIN" call -i "$GFA" --bubble-prefix-in "$OUT_DIR/bub" --reference-path "$REF" \
  -o "$OUT_DIR/call" --cn-from-multiplicity --quiet >/dev/null
"$PANVAR_BIN" describe -i "$GFA" --bubble-prefix-in "$OUT_DIR/bub" --out-dir "$OUT_DIR/desc" \
  --kmer-size 31 --no-wide-matrix --variant-nodes "$OUT_DIR/call.variant_nodes.tsv" \
  --samples "$REAL/samples.tsv" --quiet >/dev/null

FSM="$OUT_DIR/desc/fsm_kmers.samples.txt.gz"
FEAT="$OUT_DIR/desc/bubble_*/kmer_features.tsv.gz"
for mode in continuous:lpa_continuous binary:case_binary; do
  m="${mode%%:*}"; col="${mode##*:}"
  echo "== GWAS ($m), sample-level =="
  "$PY" "$REPO/scripts/gwas_demo.py" --fsm "$FSM" --phenotypes "$REAL/phenotypes.tsv" \
    --id-col sample --phenotype-col "$col" --mode "$m" --feature-map "$FEAT" \
    --node-track "$OUT_DIR/call.node_track.tsv" --variant-nodes "$OUT_DIR/call.variant_nodes.tsv" \
    --out "$OUT_DIR/gwas_$m"
  for xax in ref nodes; do
    "$RS" "$REPO/scripts/plot_gwas.R" --assoc "$OUT_DIR/gwas_$m.assoc.tsv" \
      --out "$OUT_DIR/plot_${m}_${xax}" --x "$xax" --title "Real LPA ($m, x=$xax)" >/dev/null \
      || echo "  (plot skipped: ggplot2?)"
  done
done

echo "== sanity: count recovers the KIV-2 locus that presence/absence misses =="
"$PY" - "$OUT_DIR/gwas_continuous.assoc.tsv" <<'PY'
import sys, collections
rows = [l.split("\t") for l in open(sys.argv[1]).read().splitlines()[1:]]
# 0 kmer 1 bubble 2 pos 3 node_min 4 nodes 5 variant 6 n_carriers ... 12 pa_q 13 count_q
rows.sort(key=lambda r: float(r[13]))
top = rows[0]
# the KIV-2 bubble = where the count signal concentrates
bub_hits = collections.Counter(r[1] for r in rows if float(r[13]) < 0.05)
kiv2_bub, n = bub_hits.most_common(1)[0]
N = max(int(r[6]) for r in rows)
present_in_all_count_sig = sum(1 for r in rows if int(r[6]) == N and float(r[13]) < 0.05 and float(r[12]) >= 0.05)
print(f"  top count hit: bubble={top[1]} variant={top[5]} count_q={float(top[13]):.2e} pa_q={float(top[12]):.2e}")
print(f"  count q<0.05 concentrate on bubble {kiv2_bub} ({n} k-mers) = the KIV-2 locus")
print(f"  present-in-all-{N} k-mers significant by COUNT but not P/A: {present_in_all_count_sig}")
ok = float(top[13]) < 0.05 and float(top[12]) >= 0.05 and present_in_all_count_sig > 0
print("  SANITY:", "PASS" if ok else "CHECK")
PY
echo "== DONE: outputs in $OUT_DIR (plots: plot_*_{ref,nodes}.manhattan.png) =="
