#!/usr/bin/env bash
# Benchmark the calls of call_examples.sh at 10 bp, so differences under call's 50 bp count too, and plot what
# call's events and its VCF leave, by cause (benchmark.png).
#   scripts/benchmark_examples.sh [panvar binary] [output dir] [refine_examples.sh output dir]
#                                 [call_examples.sh output dir]
set -euo pipefail
here=$PWD
cd "$(dirname "$0")/.."
arg() { if [ -z "$1" ]; then echo "$2"; elif [ "${1#/}" != "$1" ]; then echo "$1"; else echo "$here/$1"; fi; }
# bm_result <benchmark log>: its result lines, on one line.
bm_result() { grep -E '\] (haplotypes scored|present|missed by call|lost in the VCF|carriers disagreeing)' "$1" |
  sed 's/^\[benchmark [0-9:]*\] //' | paste -sd ';' - | sed 's/; */; /g'; }
run() { local log=$1; shift; "$@" 2> "$log" || { echo "failed: $* (log: $log)"; tail -n 3 "$log"; exit 1; }; }
command -v Rscript > /dev/null || { echo "Rscript not found on PATH"; exit 1; }

PANVAR=$(arg "${1:-}" build/panvar)
OUT=$(arg "${2:-}" tests/examples/benchmark)
REFINED=$(arg "${3:-}" tests/examples/refine)
CALLED=$(arg "${4:-}" tests/examples/call)

prefixes=()
for name in lpa lpa_lossy c4 gstm cyp2d acot myom2_rebuilt; do
  gfa=$REFINED/$name/$name.normalized.sorted.gfa call=$CALLED/$name/$name d=$OUT/$name
  [ -f "$call.region.vcf" ] || { echo "no $call.region.vcf: run scripts/call_examples.sh first"; exit 1; }
  mkdir -p "$d"
  run "$d/$name.benchmark.log" "$PANVAR" benchmark -i "$gfa" -b "$REFINED/$name/$name" -r GRCh38 --calls "$call" \
    -o "$d/$name" --min-sv-bp 10 --threads 2
  echo "$name: $(bm_result "$d/$name.benchmark.log")"
  prefixes+=("$d/$name")
done
Rscript scripts/plot_benchmark.R --call-min-bp 50 --out "$OUT/benchmark.png" "${prefixes[@]}"
