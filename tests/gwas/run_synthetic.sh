#!/usr/bin/env bash
# Synthetic LPA-CNV SAMPLE-LEVEL k-mer GWAS demo: prove that k-mer MULTIPLICITY recovers a
# copy-number trait (LPA KIV-2) that presence/absence is blind to. Models a haplotype panel
# (graph paths) + a diploid cohort (cosigt-style samples.tsv); describe --samples aggregates
# per-haplotype counts into per-sample dosage; GWAS tests samples vs phenotype. Self-checks the
# designed-in ground truth (tests/gwas/synthetic/truth.tsv).
#
#   run_synthetic.sh <panvar_bin> <out_dir> [python] [Rscript]
set -euo pipefail

PANVAR_BIN="${1:?usage: run_synthetic.sh <panvar_bin> <out_dir> [python] [Rscript]}"
OUT_DIR="${2:?need out_dir}"
PY="${3:-python3}"
RS="${4:-Rscript}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
DATA="$HERE/synthetic"

mkdir -p "$OUT_DIR"
echo "== regenerate synthetic cohort (haplotype panel + diploid samples) =="
"$PY" "$HERE/make_lpa_cnv.py"

echo "== bubble -> call -> describe --samples =="
"$PANVAR_BIN" bubble  -i "$DATA/lpa.gfa" -o "$OUT_DIR/bub" --snarls-in "$DATA/lpa.snarls.jsonl" --quiet >/dev/null
"$PANVAR_BIN" call    -i "$DATA/lpa.gfa" --bubble-prefix-in "$OUT_DIR/bub" --reference-path "hap00" \
  -o "$OUT_DIR/call" --cn-from-multiplicity --quiet >/dev/null
"$PANVAR_BIN" describe -i "$DATA/lpa.gfa" --bubble-prefix-in "$OUT_DIR/bub" \
  --out-dir "$OUT_DIR/desc" --kmer-size 31 --no-wide-matrix --samples "$DATA/samples.tsv" --quiet >/dev/null

FSM="$OUT_DIR/desc/fsm_kmers.samples.txt.gz"      # SAMPLE-level dosage file
FEAT="$OUT_DIR/desc/bubble_*/kmer_features.tsv.gz"
[[ -f "$FSM" ]] || { echo "error: sample fsm not written (describe --samples failed)"; exit 1; }

for mode in continuous:lpa_continuous binary:case_binary; do
  m="${mode%%:*}"; col="${mode##*:}"
  echo "== GWAS ($m / $col), sample-level =="
  "$PY" "$REPO/scripts/gwas_demo.py" --fsm "$FSM" --phenotypes "$DATA/phenotypes.tsv" \
    --id-col sample --phenotype-col "$col" --mode "$m" --feature-map "$FEAT" \
    --node-track "$OUT_DIR/call.node_track.tsv" --variant-nodes "$OUT_DIR/call.variant_nodes.tsv" \
    --out "$OUT_DIR/gwas_$m"
  for xax in ref nodes; do
    "$RS" "$REPO/scripts/plot_gwas.R" --assoc "$OUT_DIR/gwas_$m.assoc.tsv" \
      --out "$OUT_DIR/plot_${m}_${xax}" --x "$xax" --title "Synthetic LPA-CNV ($m, x=$xax)" >/dev/null \
      || echo "  (plot skipped: ggplot2?)"
  done
done

echo "== self-check against ground truth =="
"$PY" - "$DATA/truth.tsv" "$OUT_DIR/gwas_continuous.assoc.tsv" <<'PY'
import sys
truth, assoc = sys.argv[1], sys.argv[2]
tnodes = {l.split("\t")[0]: set(l.split("\t")[1].split(",")) for l in open(truth).read().splitlines()[1:]}
rows = [l.split("\t") for l in open(assoc).read().splitlines()[1:]]
# columns: 0 kmer 1 bubble 2 pos 3 node_min 4 nodes 5 variant 6 n_carriers 7 max_count
#          8 pa_p 9 count_p 10 pa_bonf 11 count_bonf 12 pa_q 13 count_q
def nodes_of(r): return set(r[4].replace(",", ";").split(";"))
def has(ns, r): return bool(ns & nodes_of(r))
N = max(int(r[6]) for r in rows)
ok = True
# KIV-2 pure multiplicity: unit k-mers present in EVERY sample -> count-sig, P/A-not
sel = [r for r in rows if has(tnodes["KIV2_VNTR"], r)]
pure = [r for r in sel if int(r[6]) == N and nodes_of(r) <= tnodes["KIV2_VNTR"]]
cmin = min(float(r[13]) for r in pure); pamax = max(float(r[12]) for r in pure)
print(f"  KIV-2 present-in-all unit k-mers: {len(pure)} | min count_q={cmin:.1e} | max pa_q={pamax:.1e}")
if not pure or cmin >= 0.05: ok = False; print("  FAIL: KIV-2 not count-significant")
if pure and pamax < 0.05: ok = False; print("  FAIL: KIV-2 should NOT be P/A-significant")
# decoy A: neither ; decoy B: both
for loc, want in (("decoy_A", False), ("decoy_B", True)):
    s = [r for r in rows if has(tnodes[loc], r)]
    cm = min(float(r[13]) for r in s); pm = min(float(r[12]) for r in s)
    print(f"  {loc}: {len(s)} k-mers | min count_q={cm:.1e} | min pa_q={pm:.1e}")
    sigboth = cm < 0.05 and pm < 0.05
    if want and not sigboth: ok = False; print(f"  FAIL: {loc} should be significant in both")
    if not want and (cm < 0.05 or pm < 0.05): ok = False; print(f"  FAIL: {loc} should be non-significant")
print("  SELF-CHECK:", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
PY

echo "== pyseer cross-check (presence/absence, optional) =="
if command -v pyseer >/dev/null 2>&1; then
  awk -F'\t' 'NR==1{print "samples\tbinary"} NR>1{print $1"\t"$5}' "$DATA/phenotypes.tsv" > "$OUT_DIR/pheno.pyseer.tsv"
  pyseer --phenotypes "$OUT_DIR/pheno.pyseer.tsv" --kmers "$FSM" --no-distances --min-af 0 --max-af 1 \
    > "$OUT_DIR/pyseer.binary.tsv" 2> "$OUT_DIR/pyseer.binary.log" || true
  echo "  pyseer rows: $(($(wc -l < "$OUT_DIR/pyseer.binary.tsv")-1)); KIV-2 (AF=1) k-mers are pre-filtered -> see $OUT_DIR/pyseer.binary.log"
else
  echo "  pyseer not installed; install with: pip install pyseer   (then re-run)"
fi
echo "== DONE: outputs in $OUT_DIR =="
