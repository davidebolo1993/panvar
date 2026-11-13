#!/usr/bin/env bash
# Count a synthetic genome set (a reference with a locus and a diverged paralog of it, and two more genomes
# holding the locus) two ways: all genomes at once, and grown step by step (the reference, then a genome
# added, then one counted in a separate job and attached), merged and compacted; check that both hold the
# same counts, and locate the locus in every genome.
#   scripts/background_examples.sh [panvar binary] [output dir] [threads]
set -euo pipefail
here=$PWD
cd "$(dirname "$0")/.."
# Relative paths given as arguments are read from the calling directory.
arg() { if [ -z "$1" ]; then echo "$2"; elif [ "${1#/}" != "$1" ]; then echo "$1"; else echo "$here/$1"; fi; }
# Run a step with stderr to a log; on failure show the log's end and stop.
run() { local log=$1; shift; "$@" 2> "$log" || { echo "failed: $* (log: $log)"; tail -n 3 "$log"; exit 1; }; }
command -v python3 > /dev/null || { echo "python3 not found on PATH"; exit 1; }

PANVAR=$(arg "${1:-}" build/panvar)
OUT=$(arg "${2:-}" tests/examples/background)
THREADS=${3:-2}
rm -rf "$OUT"
mkdir -p "$OUT/tmp"
python3 tests/tools/locus_fixture.py "$OUT/data" > /dev/null
D=$OUT/data
BG=("$PANVAR" background -t "$THREADS" --partitions 64 --memory 512M --tmp-dir "$OUT/tmp")
genomes=(--genome ref="$D/ref.fa" --genome ga="$D/ga.fa" --genome gb="$D/gb.fa")

# all at once
run "$OUT/once.log" "${BG[@]}" -o "$OUT/once.bgk" "${genomes[@]}" --coordinate-reference ref
# grown: the reference with its anchors, a genome added, one counted apart and attached, merged, compacted
run "$OUT/grown.1.log" "${BG[@]}" -o "$OUT/grown.bgk" --genome ref="$D/ref.fa" --coordinate-reference ref
run "$OUT/grown.2.log" "${BG[@]}" -o "$OUT/grown.bgk" --add --genome ga="$D/ga.fa"
run "$OUT/grown.3.log" "${BG[@]}" -o "$OUT/grown.bgk" --segment-only "$OUT/seg_gb" --genome gb="$D/gb.fa"
run "$OUT/grown.4.log" "${BG[@]}" -o "$OUT/grown.bgk" --attach "$OUT/seg_gb"
run "$OUT/grown.5.log" "${BG[@]}" -o "$OUT/grown.bgk" --merge
run "$OUT/grown.6.log" "${BG[@]}" -o "$OUT/grown.bgk" --compact
echo "once:   $(grep -h 'segment: ' "$OUT/once.log" | sed 's/^\[background [0-9:]*\] //')"
echo "grown:  $(grep -h 'merged: ' "$OUT/grown.5.log" | sed 's/^\[background [0-9:]*\] //')"

# the same counts either way, in the store and in its dictionary
for other in "$OUT/grown.bgk" "$OUT/grown.bgk.sshash"; do
  run "$OUT/verify.$(basename "$other").log" "${BG[@]}" -o "$OUT/once.bgk" --verify "$other"
  echo "verify $(basename "$other"): $(grep -h 'checked: ' "$OUT/verify.$(basename "$other").log" | sed 's/^\[background [0-9:]*\] //')"
done
# the locus in every genome, from the reference's anchors
run "$OUT/locate.log" "${BG[@]}" -o "$OUT/grown.bgk" --locate chr1:120000-160000 > "$OUT/locate.tsv"
column -t "$OUT/locate.tsv"
rm -rf "$OUT/tmp"
