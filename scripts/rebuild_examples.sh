#!/usr/bin/env bash
# Rebuild a few loci, seeded by the reference (-r) and by the richest haplotype, and compare each with its input.
#   scripts/rebuild_examples.sh [panvar binary] [output dir] [threads]
set -euo pipefail
here=$PWD
cd "$(dirname "$0")/.."
# Relative paths given as arguments are read from the calling directory.
arg() { if [ -z "$1" ]; then echo "$2"; elif [ "${1#/}" != "$1" ]; then echo "$1"; else echo "$here/$1"; fi; }
# Run a step with stderr to a log; on failure show the log's end and stop.
run() { local log=$1; shift; "$@" 2> "$log" || { echo "failed: $* (log: $log)"; tail -n 3 "$log"; exit 1; }; }
for tool in python3 Rscript; do command -v "$tool" > /dev/null || { echo "$tool not found on PATH"; exit 1; }; done

PANVAR=$(arg "${1:-}" build/panvar)
OUT=$(arg "${2:-}" tests/examples/rebuild)
THREADS=${3:-2}
DATA=tests/real_data
mkdir -p "$OUT"

rebuild() {  # <locus> <stage> <input graph> [rebuild options]
  local locus=$1 stage=$2 gfa=$3
  shift 3
  local out=$OUT/$locus/$locus.$stage
  mkdir -p "$OUT/$locus"
  run "$out.log" "$PANVAR" rebuild -i "$gfa" -o "$out.gfa" -t "$THREADS" "$@"
  # a rejected rebuild writes its input unchanged: label it so
  local label=$stage
  grep -q "warning: rejected:" "$out.log" && label=${stage}_rejected
  printf '%s\t%s\t%s\n' "$locus" "$label" "$out.gfa" >> "$OUT/graphs.tsv"
  echo "$locus $stage: $(grep -iE "accepted|rejected|input is healthy" "$out.log" | tail -n 1 |
    sed 's/^\[rebuild [0-9:]*\] //')"
}

printf 'locus\tstage\tpath\n' > "$OUT/graphs.tsv"

myom2=$DATA/rebuild/myom2.gfa.gz
printf 'myom2\tinput\t%s\n' "$myom2" >> "$OUT/graphs.tsv"
rebuild myom2 reference "$myom2" -r GRCh38 --base-align --clip-to-reference
rebuild myom2 seed "$myom2"

# The paralog graphs pass the gate: --force rebuilds them.
for locus in c4 gstm cyp2d acot; do
  printf '%s\tinput\t%s\n' "$locus" "$DATA/paralogs/$locus.gfa.gz" >> "$OUT/graphs.tsv"
  rebuild "$locus" reference "$DATA/paralogs/$locus.gfa.gz" -r GRCh38 --force
  rebuild "$locus" seed "$DATA/paralogs/$locus.gfa.gz" --force
done

python3 scripts/compare_graphs.py --list "$OUT/graphs.tsv" -o "$OUT/compare"
Rscript scripts/plot_compare_graphs.R --prefix "$OUT/compare"
