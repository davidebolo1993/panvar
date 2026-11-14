#!/usr/bin/env bash
# Project the loci GTF onto every haplotype of the paralog graphs, and print each gene's copies per haplotype.
#   scripts/transcripts_examples.sh [panvar binary] [output dir] [GTF]
set -euo pipefail
here=$PWD
cd "$(dirname "$0")/.."
# Relative paths given as arguments are read from the calling directory.
arg() { if [ -z "$1" ]; then echo "$2"; elif [ "${1#/}" != "$1" ]; then echo "$1"; else echo "$here/$1"; fi; }
# Run a step with stderr to a log; on failure show the log's end and stop.
run() { local log=$1; shift; "$@" 2> "$log" || { echo "failed: $* (log: $log)"; tail -n 3 "$log"; exit 1; }; }

PANVAR=$(arg "${1:-}" build/panvar)
OUT=$(arg "${2:-}" tests/examples/transcripts)
GTF=$(arg "${3:-}" tests/real_data/Homo_sapiens.GRCh38.116.loci.gtf.gz)
[ -f "$GTF" ] || { echo "no GTF at $GTF"; exit 1; }
mkdir -p "$OUT"

for name in c4 gstm cyp2d acot; do
  mkdir -p "$OUT/$name"
  run "$OUT/$name/$name.transcripts.log" "$PANVAR" transcripts -i "tests/real_data/paralogs/$name.gfa.gz" -r GRCh38 \
    --gtf "$GTF" -o "$OUT/$name/$name" -t 2
  echo "$name: $(grep -o 'haplotype-specific transcripts: [0-9]*' "$OUT/$name/$name.transcripts.log")"
  awk -F'\t' 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next} {printf "  %-16s %s\n", $2, $c["haplotypes_by_copies"]}' \
    "$OUT/$name/$name.hst.summary.tsv"
done
