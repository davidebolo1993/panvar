#!/usr/bin/env bash
# Compare call with svim-asm on the haplotypes of each benchmark_examples.sh graph: svim-asm's calls
# (assembly_variants.py) and call's are scored against benchmark's present.tsv (compare_sv_callers.py) and
# plotted per graph (callers.png) and by event size (callers.by_size.png).
#   dev/benchmark_callers.sh [benchmark_examples.sh output dir] [refine_examples.sh output dir]
#                            [call_examples.sh output dir] [svim-asm tools dir] [truvari tools dir]
# Tools dirs: minimap2, k8, paftools.js, samtools and svim-asm; bcftools, tabix and truvari (default: PATH).
set -euo pipefail
here=$PWD
cd "$(dirname "$0")/.."
arg() { if [ -z "$1" ]; then echo "$2"; elif [ "${1#/}" != "$1" ]; then echo "$1"; else echo "$here/$1"; fi; }
run() { local log=$1; shift; "$@" 2> "$log" || { echo "failed: $* (log: $log)"; tail -n 3 "$log"; exit 1; }; }
for tool in python3 Rscript; do command -v "$tool" > /dev/null || { echo "$tool not found on PATH"; exit 1; }; done

OUT=$(arg "${1:-}" tests/examples/benchmark)
REFINED=$(arg "${2:-}" tests/examples/refine)
CALLED=$(arg "${3:-}" tests/examples/call)
svim=$(command -v svim-asm || true)
TOOLS=$(arg "${4:-}" "${svim%/svim-asm}")
[ -x "$TOOLS/svim-asm" ] || { echo "no svim-asm: give the directory holding minimap2, k8, paftools.js, samtools and svim-asm"; exit 1; }
truvari=$(command -v truvari || true)
MERGE=$(arg "${5:-}" "${truvari%/truvari}")
[ -x "$MERGE/truvari" ] || { echo "no truvari: give the directory holding bcftools, tabix and truvari"; exit 1; }

tables=()
for name in lpa lpa_lossy c4 gstm cyp2d acot myom2_rebuilt; do
  gfa=$REFINED/$name/$name.normalized.sorted.gfa call=$CALLED/$name/$name d=$OUT/$name
  [ -f "$d/$name.present.tsv" ] || { echo "no $d/$name.present.tsv: run scripts/benchmark_examples.sh first"; exit 1; }
  mkdir -p "$d/svim"
  run "$d/svim/$name.log" python3 dev/assembly_variants.py --gfa "$gfa" --reference GRCh38 --tools "$TOOLS" \
    --merge-tools "$MERGE" --out-dir "$d/svim" --threads 2 --sv-tables "$d/svim/$name"
  python3 dev/compare_sv_callers.py --label "$name" --benchmark "$d/$name" --calls "$call" --svim "$d/svim/$name" \
    -o "$d/$name.callers.tsv" --details "$d/$name.callers.events.tsv"
  tables+=("$d/$name.callers.tsv")
done
Rscript dev/plot_sv_callers.R --out "$OUT/callers.png" "${tables[@]}"
Rscript dev/plot_sv_callers.R --by-size --out "$OUT/callers.by_size.png" "${tables[@]}"
