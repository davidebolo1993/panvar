#!/usr/bin/env bash
# LPA pangenome-association demo on the real LPA graph (tests/real_data/lpa.gfa.gz).
#
# Pipeline: bubble -> panphorte -> call -> describe --samples -> association. KIV-2 copy number is
# read per real haplotype from panphorte; real haplotypes are paired into diploid samples with a
# literature-plausible inverse Lp(a) phenotype (synthetic values, real topology). The association is
# run on BOTH feature substrates and BOTH trait codings:
#
#   substrate kmer  : k-mer multiplicity   (describe fsm_kmers.samples.txt.gz)  -- primary, tested
#   substrate graph : node/edge dosage     (describe fsm_graph.samples.txt.gz)  -- complementary
#   mode continuous : Lp(a) level          mode binary : case/control (median split)
#
# In every case the COUNT model recovers the KIV-2 copy-number locus that a PRESENCE/ABSENCE model
# is blind to. call and describe consume the panphorte-normalized/sorted graph + panphorte prefix.
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
# Prefer a GRCh38 path as the reference (else the first path). awk reads all input (no early exit, so
# no SIGPIPE under `set -o pipefail`).
REF="$(gzcat "$GFA" | awk -F'\t' '($1=="P"||$1=="W"){n=$2; if(g==""&&n~/[Gg][Rr][Cc]h38/)g=n; if(f=="")f=n} END{print (g!=""?g:f)}')"
echo "reference path: $REF"

echo "== bubble -> panphorte (KIV-2 copies, normalized+sorted graph) =="
"$PANVAR_BIN" bubble    -i "$GFA" -o "$OUT_DIR/bub" -r "$REF" --quiet >/dev/null
SGFA="$OUT_DIR/bub.sorted.gfa"
"$PANVAR_BIN" panphorte -i "$SGFA" --bubble-prefix-in "$OUT_DIR/bub" -o "$OUT_DIR/pan" \
  --reference-path "$REF" --min-similarity 0.90 --quiet >/dev/null
PGFA="$OUT_DIR/pan.normalized.sorted.gfa"   # panphorte graph: the input for call + describe
"$PY" "$HERE/make_lpa_phenotype.py" "$OUT_DIR/pan.panphorte.copies.tsv" "$REAL" 200

echo "== call -> describe --samples (variant-restricted, both substrates) =="
"$PANVAR_BIN" call -i "$PGFA" --bubble-prefix-in "$OUT_DIR/pan" --reference-path "$REF" \
  -o "$OUT_DIR/call" --cn-from-multiplicity --quiet >/dev/null
"$PANVAR_BIN" describe -i "$PGFA" --bubble-prefix-in "$OUT_DIR/pan" --out-dir "$OUT_DIR/desc" \
  --kmer-size 31 --no-wide-matrix --variant-nodes "$OUT_DIR/call.variant_nodes.tsv" \
  --samples "$REAL/samples.tsv" --quiet >/dev/null

NT="$OUT_DIR/call.node_track.tsv"; VN="$OUT_DIR/call.variant_nodes.tsv"
KFSM="$OUT_DIR/desc/fsm_kmers.samples.txt.gz"; KFEAT="$OUT_DIR/desc/bubble_*/kmer_features.tsv.gz"
GFSM="$OUT_DIR/desc/fsm_graph.samples.txt.gz"
for sub in kmer graph; do
  fsm="$KFSM"; extra=(--feature-map "$KFEAT")
  [ "$sub" = graph ] && { fsm="$GFSM"; extra=(); }
  for mode in continuous:lpa_continuous binary:case_binary; do
    m="${mode%%:*}"; col="${mode##*:}"
    echo "== association ($sub / $m) =="
    "$PY" "$REPO/scripts/gwas_demo.py" --fsm "$fsm" --phenotypes "$REAL/phenotypes.tsv" \
      --id-col sample --phenotype-col "$col" --mode "$m" --substrate "$sub" ${extra[@]+"${extra[@]}"} \
      --node-track "$NT" --variant-nodes "$VN" --out "$OUT_DIR/gwas_${sub}_${m}"
    "$RS" "$REPO/scripts/plot_gwas.R" --assoc "$OUT_DIR/gwas_${sub}_${m}.assoc.tsv" \
      --out "$OUT_DIR/plot_${sub}_${m}" --x nodes --title "LPA ($sub, $m)" >/dev/null \
      || echo "  (plot skipped: ggplot2?)"
  done
  # one faceted figure per substrate: continuous vs binary x count vs presence/absence, so the
  # stronger-association trait/test reads off the peak heights at a glance.
  label="k-mer"; [ "$sub" = graph ] && label="graph (node/edge)"
  "$RS" "$REPO/scripts/plot_gwas_compare.R" --x nodes --title "LPA $label" \
    --assoc "continuous=$OUT_DIR/gwas_${sub}_continuous.assoc.tsv" \
    --assoc "binary=$OUT_DIR/gwas_${sub}_binary.assoc.tsv" \
    --out "$OUT_DIR/gwas_${sub}_compare" >/dev/null || echo "  (compare plot skipped: ggplot2?)"
done

echo "== sanity: COUNT recovers the KIV-2 locus that PRESENCE/ABSENCE misses =="
"$PY" - "$OUT_DIR/gwas_kmer_continuous.assoc.tsv" "$OUT_DIR/gwas_graph_continuous.assoc.tsv" <<'PY'
import sys
for path in sys.argv[1:]:
    rows = [l.split("\t") for l in open(path).read().splitlines()[1:]]
    # cols: 0 feature 1 bubble 5 variant 6 n_carriers 12 pa_q 13 count_q
    rows.sort(key=lambda r: float(r[13]))
    top = rows[0]
    N = max(int(r[6]) for r in rows)
    all_count_only = sum(1 for r in rows if int(r[6]) == N and float(r[13]) < 0.05 and float(r[12]) >= 0.05)
    ok = float(top[13]) < 0.05 and float(top[12]) >= 0.05 and all_count_only > 0
    print(f"  {path.split('/')[-1]:<32} top variant={top[5]} count_q={float(top[13]):.1e} "
          f"pa_q={float(top[12]):.1e}  present-in-all & count-only-sig={all_count_only}  -> "
          f"{'PASS' if ok else 'CHECK'}")
PY
echo "== DONE: outputs in $OUT_DIR (plots: plot_{kmer,graph}_{continuous,binary}.manhattan.png) =="
