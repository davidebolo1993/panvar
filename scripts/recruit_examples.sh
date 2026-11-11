#!/usr/bin/env bash
# Index the synthetic locus of background_examples.sh against its grown store, then recruit pairs drawn from
# the locus, from its diverged paralog elsewhere in the reference, and from unrelated sequence: the locus
# pairs should be certain, the paralog's never.
#   scripts/recruit_examples.sh [panvar binary] [output dir] [background_examples.sh output dir] [threads]
set -euo pipefail
here=$PWD
cd "$(dirname "$0")/.."
# Relative paths given as arguments are read from the calling directory.
arg() { if [ -z "$1" ]; then echo "$2"; elif [ "${1#/}" != "$1" ]; then echo "$1"; else echo "$here/$1"; fi; }
# Run a step with stderr to a log; on failure show the log's end and stop.
run() { local log=$1; shift; "$@" 2> "$log" || { echo "failed: $* (log: $log)"; tail -n 3 "$log"; exit 1; }; }
command -v python3 > /dev/null || { echo "python3 not found on PATH"; exit 1; }

PANVAR=$(arg "${1:-}" build/panvar)
OUT=$(arg "${2:-}" tests/examples/recruit)
BGDIR=$(arg "${3:-}" tests/examples/background)
THREADS=${4:-2}
[ -d "$BGDIR/grown.bgk" ] || { echo "no $BGDIR/grown.bgk: run scripts/background_examples.sh first"; exit 1; }
D=$BGDIR/data
rm -rf "$OUT"
mkdir -p "$OUT"

# the index, with its calibration (kept in $OUT/index; a rerun with the same graph and store reuses it)
run "$OUT/index.log" "$PANVAR" recruit --locus "$D/locus.gfa" --background "$BGDIR/grown.bgk" -r ref --save-index "$OUT/index" \
  --index-only -t "$THREADS"
grep -h "locus: locus; private\|calibration: short" "$OUT/index.log" | sed 's/^\[recruit [0-9:]*\] //'

# a sample: 2,000 pairs from each source, error-free, the source in each read's name
for src in "loc 120000 160000 1" "par 231000 259000 2" "far 0 100000 3"; do
  read -r name s e seed <<< "$src"
  python3 tests/tools/locus_reads.py "$D/ref.fa" chr1 "$s" "$e" 2000 "$name" "$seed" "$OUT/$name"
done
cat "$OUT"/loc_1.fq "$OUT"/par_1.fq "$OUT"/far_1.fq > "$OUT/sample_1.fq"
cat "$OUT"/loc_2.fq "$OUT"/par_2.fq "$OUT"/far_2.fq > "$OUT/sample_2.fq"
rm -f "$OUT"/{loc,par,far}_[12].fq
run "$OUT/sample.log" "$PANVAR" recruit --locus "$OUT/index/locus.plx" -1 "$OUT/sample_1.fq" -2 "$OUT/sample_2.fq" \
  -o "$OUT/sample.tsv" --reads-out "$OUT/reads" -t "$THREADS"
grep -h "site cut" "$OUT/sample.log" | sed 's/^\[recruit [0-9:]*\] //'

# fates by source (fragments the screen passed over are not in the table)
{ printf 'source\tpairs\tscreened\tcertain\tambiguous\tunrecruited\n'
  for src in loc par far; do
    awk -F'\t' -v s="$src" 'NR > 1 && substr($1, 1, 3) == s {n++; f[$3]++}
      END {printf "%s\t2000\t%d\t%d\t%d\t%d\n", s, n, f["certain"], f["ambiguous"], f["unrecruited"]}' "$OUT/sample.tsv"
  done
} > "$OUT/fates.tsv"
column -t "$OUT/fates.tsv"
