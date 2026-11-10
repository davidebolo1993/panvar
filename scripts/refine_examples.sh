#!/usr/bin/env bash
# Refine the graphs from panphorte_examples.sh and compare each with its input.
#   scripts/refine_examples.sh [panvar binary] [output dir] [panphorte_examples.sh output dir]
set -euo pipefail
here=$PWD
cd "$(dirname "$0")/.."
# Relative paths given as arguments are read from the calling directory.
arg() { if [ -z "$1" ]; then echo "$2"; elif [ "${1#/}" != "$1" ]; then echo "$1"; else echo "$here/$1"; fi; }
# Run a step with stderr to a log; on failure show the log's end and stop.
run() { local log=$1; shift; "$@" 2> "$log" || { echo "failed: $* (log: $log)"; tail -n 3 "$log"; exit 1; }; }
for tool in python3 Rscript; do command -v "$tool" > /dev/null || { echo "$tool not found on PATH"; exit 1; }; done

PANVAR=$(arg "${1:-}" build/panvar)
OUT=$(arg "${2:-}" tests/examples/refine)
FOLDED=$(arg "${3:-}" tests/examples/panphorte)
mkdir -p "$OUT"
printf 'locus\tstage\tpath\n' > "$OUT/graphs.tsv"

refine() {  # <name> <graph> <site table prefix>
  local name=$1 gfa=$2 sites=$3
  mkdir -p "$OUT/$name"
  run "$OUT/$name/$name.refine.log" "$PANVAR" refine -i "$gfa" -b "$sites" -r GRCh38 -o "$OUT/$name/$name"
  printf '%s\tinput\t%s\n%s\trefined\t%s\n' "$name" "$gfa" "$name" "$OUT/$name/$name.normalized.sorted.gfa" \
    >> "$OUT/graphs.tsv"
  echo "$name: $(awk -F'\t' 'NR>1 {n[$5 ":" $6]++} END {for (k in n) printf "%s %d; ", k, n[k]}' \
    "$OUT/$name/$name.refine.report.tsv")"
}

for name in lpa lpa_lossy c4 gstm cyp2d acot myom2_rebuilt; do
  [ -f "$FOLDED/$name/$name.bubbles.csv" ] || { echo "no $FOLDED/$name/$name.bubbles.csv: run scripts/panphorte_examples.sh first"; exit 1; }
  refine "$name" "$FOLDED/$name/$name.normalized.sorted.gfa" "$FOLDED/$name/$name"
done

python3 scripts/compare_graphs.py --list "$OUT/graphs.tsv" -o "$OUT/compare"
Rscript scripts/plot_compare_graphs.R --prefix "$OUT/compare"
