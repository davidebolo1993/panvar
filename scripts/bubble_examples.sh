#!/usr/bin/env bash
# Find the sites of the paralog and tandem-repeat graphs, and of the tangled locus rebuilt by rebuild_examples.sh.
#   scripts/bubble_examples.sh [panvar binary] [output dir] [rebuild_examples.sh output dir]
set -euo pipefail
here=$PWD
cd "$(dirname "$0")/.."
# Relative paths given as arguments are read from the calling directory.
arg() { if [ -z "$1" ]; then echo "$2"; elif [ "${1#/}" != "$1" ]; then echo "$1"; else echo "$here/$1"; fi; }
# Run a step with stderr to a log; on failure show the log's end and stop.
run() { local log=$1; shift; "$@" 2> "$log" || { echo "failed: $* (log: $log)"; tail -n 3 "$log"; exit 1; }; }

PANVAR=$(arg "${1:-}" build/panvar)
OUT=$(arg "${2:-}" tests/examples/bubble)
REBUILT=$(arg "${3:-}" tests/examples/rebuild)
DATA=tests/real_data

bubble() {  # <name> <graph>
  mkdir -p "$OUT/$1"
  run "$OUT/$1/$1.bubble.log" "$PANVAR" bubble -i "$2" -r GRCh38 -o "$OUT/$1/$1"
  echo "$1: $(grep -E '\] sites: ' "$OUT/$1/$1.bubble.log" | sed 's/^\[bubble [0-9:]*\] //')"
}

for locus in c4 gstm cyp2d acot; do
  bubble "$locus" "$DATA/paralogs/$locus.gfa.gz"
done
bubble lpa "$DATA/tandem_repeat/lpa.gfa.gz"

[ -f "$REBUILT/myom2/myom2.reference.gfa" ] ||
  { echo "no $REBUILT/myom2/myom2.reference.gfa: run scripts/rebuild_examples.sh first"; exit 1; }
bubble myom2_rebuilt "$REBUILT/myom2/myom2.reference.gfa"
