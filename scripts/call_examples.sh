#!/usr/bin/env bash
# Call variants and copy number (per gene with a GTF) on the graphs from refine_examples.sh. Each graph also gets
# a bundle for the viewer: VN_RDS=<output dir>/<graph>/<graph>.viewer.rds Rscript scripts/variant_node_heatmap_app.R
#   scripts/call_examples.sh [panvar binary] [output dir] [refine_examples.sh output dir] [GTF]
set -euo pipefail
here=$PWD
cd "$(dirname "$0")/.."
# Relative paths given as arguments are read from the calling directory.
arg() { if [ -z "$1" ]; then echo "$2"; elif [ "${1#/}" != "$1" ]; then echo "$1"; else echo "$here/$1"; fi; }
# Run a step with stderr to a log; on failure show the log's end and stop.
run() { local log=$1; shift; "$@" 2> "$log" || { echo "failed: $* (log: $log)"; tail -n 3 "$log"; exit 1; }; }
command -v Rscript > /dev/null || { echo "Rscript not found on PATH"; exit 1; }

PANVAR=$(arg "${1:-}" build/panvar)
OUT=$(arg "${2:-}" tests/examples/call)
REFINED=$(arg "${3:-}" tests/examples/refine)
GTF=$(arg "${4:-}" tests/real_data/Homo_sapiens.GRCh38.116.loci.gtf.gz)
genes=()
if [ -f "$GTF" ]; then genes=(--gtf "$GTF"); else echo "no GTF at $GTF: no per-gene copy numbers"; fi

call() {  # <name>
  local name=$1
  [ -f "$REFINED/$name/$name.bubbles.csv" ] ||
    { echo "no $REFINED/$name/$name.bubbles.csv: run scripts/refine_examples.sh first"; exit 1; }
  mkdir -p "$OUT/$name"
  run "$OUT/$name/$name.call.log" "$PANVAR" call -i "$REFINED/$name/$name.normalized.sorted.gfa" \
    -b "$REFINED/$name/$name" -r GRCh38 -o "$OUT/$name/$name" --cn --threads 2 ${genes[@]+"${genes[@]}"}
  echo "$name: $(grep -E "\] records: " "$OUT/$name/$name.call.log" | sed 's/^\[call [0-9:]*\] //')"
  Rscript scripts/build_variant_node_data.R --gfa "$REFINED/$name/$name.normalized.sorted.gfa" \
    --prefix "$OUT/$name/$name" --bubbles "$REFINED/$name/$name.bubbles.csv" > /dev/null
}

for name in lpa lpa_lossy c4 gstm cyp2d acot myom2_rebuilt; do
  call "$name"
done
