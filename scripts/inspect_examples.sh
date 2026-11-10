#!/usr/bin/env bash
# Inspect the longest site of a few loci from bubble_examples.sh (clusters, allele tables) and draw heatmaps of
# node and link passes: every path, and one row per cluster.
#   scripts/inspect_examples.sh [panvar binary] [output dir] [bubble_examples.sh output dir]
set -euo pipefail
here=$PWD
cd "$(dirname "$0")/.."
# Relative paths given as arguments are read from the calling directory.
arg() { if [ -z "$1" ]; then echo "$2"; elif [ "${1#/}" != "$1" ]; then echo "$1"; else echo "$here/$1"; fi; }
# Run a step with stderr to a log; on failure show the log's end and stop.
run() { local log=$1; shift; "$@" 2> "$log" || { echo "failed: $* (log: $log)"; tail -n 3 "$log"; exit 1; }; }
command -v Rscript > /dev/null || { echo "Rscript not found on PATH"; exit 1; }

PANVAR=$(arg "${1:-}" build/panvar)
OUT=$(arg "${2:-}" tests/examples/inspect)
BUBBLES=$(arg "${3:-}" tests/examples/bubble)

longest_site() {  # <bubbles.csv>
  awk -F, 'NR==1{for(i=1;i<=NF;i++)c[$i]=i; next}
           $(c["kind"]) != "inversion" && $(c["max_inside_bp"])+0 > m {m=$(c["max_inside_bp"])+0; id=$(c["bubble_id"])}
           END{print id}' "$1"
}

for locus in c4 acot gstm lpa; do
  in=$BUBBLES/$locus/$locus
  [ -f "$in.bubbles.csv" ] || { echo "no $in.bubbles.csv: run scripts/bubble_examples.sh first"; exit 1; }
  site=$(longest_site "$in.bubbles.csv")
  mkdir -p "$OUT/$locus"
  run "$OUT/$locus/$locus.inspect.log" "$PANVAR" inspect -i "$in.sorted.gfa" -b "$in" -o "$OUT/$locus" --bubble-id "$site" -r GRCh38 \
    --cluster --cluster-similarity 0.97 --genotype-table
  dir=$OUT/$locus/bubble_$site
  for features in nodes edges; do
    Rscript scripts/plot_inspect.R --table "$dir" --features "$features" --dpi 150
    Rscript scripts/plot_inspect.R --table "$dir" --features "$features" --dpi 150 --representatives-only \
      --out "$dir/$locus.bubble_$site.$features.representatives"
  done
  echo "$locus: site $site; clusters: $(( $(wc -l < "$dir/$locus.bubble_$site.clusters.tsv") - 1 ))"
done
