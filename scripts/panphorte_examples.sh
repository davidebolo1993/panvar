#!/usr/bin/env bash
# Fold the graphs from bubble_examples.sh: lpa twice (lpa keeps every copy's differences, lpa_lossy only the
# reference's) and the other loci. Then inspect lpa's longest folded repeat before and after each fold, and
# compare the three lpa graphs.
#   scripts/panphorte_examples.sh [panvar binary] [output dir] [bubble_examples.sh output dir]
set -euo pipefail
here=$PWD
cd "$(dirname "$0")/.."
# Relative paths given as arguments are read from the calling directory.
arg() { if [ -z "$1" ]; then echo "$2"; elif [ "${1#/}" != "$1" ]; then echo "$1"; else echo "$here/$1"; fi; }
# Run a step with stderr to a log; on failure show the log's end and stop.
run() { local log=$1; shift; "$@" 2> "$log" || { echo "failed: $* (log: $log)"; tail -n 3 "$log"; exit 1; }; }
for tool in python3 Rscript; do command -v "$tool" > /dev/null || { echo "$tool not found on PATH"; exit 1; }; done

PANVAR=$(arg "${1:-}" build/panvar)
OUT=$(arg "${2:-}" tests/examples/panphorte)
BUBBLES=$(arg "${3:-}" tests/examples/bubble)
locus=lpa

fold() {  # <name> <locus> [panphorte options]
  local name=$1 in=$BUBBLES/$2/$2
  shift 2
  [ -f "$in.bubbles.csv" ] || { echo "no $in.bubbles.csv: run scripts/bubble_examples.sh first"; exit 1; }
  mkdir -p "$OUT/$name"
  run "$OUT/$name/$name.panphorte.log" "$PANVAR" panphorte -i "$in.sorted.gfa" -b "$in" -o "$OUT/$name/$name" \
    -r GRCh38 --min-similarity 0.95 --threads 2 "$@"
  echo "$name: $(grep -E "\] folded:" "$OUT/$name/$name.panphorte.log" | sed 's/^\[panphorte [0-9:]*\] //')"
}

folded_site() {  # <name> <site before the fold>: the site holding its unit after the fold
  local pp=$OUT/$1/$1
  local unit_node
  unit_node=$(awk -F'\t' -v b="$2" 'NR>1 && $1==b {split($9, n, ";"); print n[1]; exit}' "$pp.panphorte.unit_graph.tsv")
  awk -F, -v n="$unit_node" 'NR==1{for(i=1;i<=NF;i++)c[$i]=i; next}
      {gsub(/"/, "", $(c["inside_nodes"])); k=split($(c["inside_nodes"]), a, ";")
       for(i=1;i<=k;i++) if (a[i]==n) {print $(c["bubble_id"]); exit}}' "$pp.bubbles.csv"
}

inspect_site() {  # <name> <graph> <site table prefix> <site>
  mkdir -p "$OUT/inspect"
  run "$OUT/inspect/$1.inspect.log" "$PANVAR" inspect -i "$2" -b "$3" -o "$OUT/inspect/$1" --bubble-id "$4" \
    -r GRCh38 --genotype-table
  Rscript scripts/plot_inspect.R --table "$OUT/inspect/$1/bubble_$4" --dpi 150
}

fold "$locus" "$locus"
fold "${locus}_lossy" "$locus" --lossy-units
for other in c4 gstm cyp2d acot myom2_rebuilt; do
  fold "$other" "$other"
done

# lpa's fold with the longest unit, inspected before and after each fold.
in=$BUBBLES/$locus/$locus
before=$(awk -F'\t' 'NR>1 && $2=="yes" && $3+0 > m {m=$3+0; id=$1} END{print id}' "$OUT/$locus/$locus.panphorte.report.tsv")
[ -n "$before" ] || { echo "$locus: nothing was folded, so there is no folded site to inspect"; exit 1; }
inspect_site "${locus}_unfolded" "$in.sorted.gfa" "$in" "$before"
for name in "$locus" "${locus}_lossy"; do
  after=$(folded_site "$name" "$before")
  echo "$name: site $before before the fold is site $after after it"
  inspect_site "${name}_folded" "$OUT/$name/$name.normalized.sorted.gfa" "$OUT/$name/$name" "$after"
done

# The three lpa graphs compared.
printf 'locus\tstage\tpath\n%s\tunfolded\t%s\n%s\tfolded\t%s\n%s\tlossy\t%s\n' \
  "$locus" "$in.sorted.gfa" "$locus" "$OUT/$locus/$locus.normalized.sorted.gfa" \
  "$locus" "$OUT/${locus}_lossy/${locus}_lossy.normalized.sorted.gfa" > "$OUT/graphs.tsv"
python3 scripts/compare_graphs.py --list "$OUT/graphs.tsv" -o "$OUT/compare"
Rscript scripts/plot_compare_graphs.R --prefix "$OUT/compare"
