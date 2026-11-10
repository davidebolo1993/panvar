#!/usr/bin/env bash
# What a lower call threshold recovers: each example graph's sites found once at 10 bp (bubble --min-variant-bp,
# panphorte and refine --resnarl-min-variant-bp), then call alone at 50, 40, 30, 20 and 10 bp (--min-sv-bp), each
# benchmarked at one size (default 10 bp). One plot per threshold, every graph in it: thresholds.<bp>.png.
#   scripts/benchmark_thresholds.sh [panvar binary] [output dir] [rebuild_examples.sh output dir] [bench min bp]
# REUSE=1 keeps the sites when their graph is newer than the binary, and a call or benchmark newer than what it
# reads (a benchmark also at the same bench min).
{   # read whole first, so editing the file during a run does not change it
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
OUT=$(arg "${2:-}" tests/examples/benchmark/thresholds)
REBUILT=$(arg "${3:-}" tests/examples/rebuild)
BENCH_MIN=${4:-10}
DATA=tests/real_data
THRESHOLDS=(50 40 30 20 10)
[ -f "$REBUILT/myom2/myom2.reference.gfa" ] ||
  { echo "no $REBUILT/myom2/myom2.reference.gfa: run scripts/rebuild_examples.sh first"; exit 1; }

newer() {  # <file> <inputs...>: REUSE=1 and the file is newer than every input
  local f
  [ "${REUSE:-0}" = 1 ] && [ -f "$1" ] || return 1
  for f in "${@:2}"; do [ -f "$f" ] && [ "$1" -nt "$f" ] || return 1; done
}
input() {  # <name>: the graph bubble_examples.sh starts from
  case $1 in
    lpa|lpa_lossy) echo "$DATA/tandem_repeat/lpa.gfa.gz" ;;
    myom2_rebuilt) echo "$REBUILT/myom2/myom2.reference.gfa" ;;
    *) echo "$DATA/paralogs/$1.gfa.gz" ;;
  esac
}

names=(lpa lpa_lossy c4 gstm cyp2d acot myom2_rebuilt)
for name in "${names[@]}"; do
  # the sites, once, at 10 bp; only the refined graph is kept
  d=$OUT/$name/sites10 s=$OUT/$name/sites10/$name
  mkdir -p "$d"
  if ! newer "$s.normalized.sorted.gfa" "$PANVAR"; then
    fold_opts=(--min-similarity 0.95)
    [ "$name" = lpa_lossy ] && fold_opts+=(--lossy-units)
    run "$d/bubble.log" "$PANVAR" bubble -i "$(input "$name")" -r GRCh38 -o "$s.bubble" --min-variant-bp 10
    run "$d/panphorte.log" "$PANVAR" panphorte -i "$s.bubble.sorted.gfa" -b "$s.bubble" -o "$s.folded" -r GRCh38 \
      --threads 2 --resnarl-min-variant-bp 10 "${fold_opts[@]}"
    run "$d/refine.log" "$PANVAR" refine -i "$s.folded.normalized.sorted.gfa" -b "$s.folded" -r GRCh38 -o "$s" \
      --resnarl-min-variant-bp 10
    rm -f "$s.bubble.sorted.gfa" "$s.folded.normalized.sorted.gfa"
  fi
  # call alone at each threshold on those sites
  for t in "${THRESHOLDS[@]}"; do
    c=$OUT/$name/call$t p=$OUT/$name/call$t/$name
    mkdir -p "$c"
    newer "$p.region.vcf" "$s.normalized.sorted.gfa" "$PANVAR" ||
      run "$c/call.log" "$PANVAR" call -i "$s.normalized.sorted.gfa" -b "$s" -r GRCh38 -o "$p" --cn --min-sv-bp "$t" \
        --no-per-bubble-vcf --threads 2
    if ! { newer "$p.summary.tsv" "$p.region.vcf" "$PANVAR" && grep -q "^min_sv_bp	$BENCH_MIN\$" "$p.summary.tsv"; }; then
      run "$c/benchmark.log" "$PANVAR" benchmark -i "$s.normalized.sorted.gfa" -b "$s" -r GRCh38 --calls "$p" \
        -o "$p" --min-sv-bp "$BENCH_MIN" --threads 2
    fi
    echo "$name, sites at 10 bp, call at $t bp: $(bm_result "$c/benchmark.log")"
  done
done

# One plot per threshold, every graph in it.
for t in "${THRESHOLDS[@]}"; do
  runs=()
  for name in "${names[@]}"; do runs+=("$name=$OUT/$name/call$t/$name:$t"); done
  Rscript scripts/plot_benchmark.R --out "$OUT/thresholds.$t.png" "${runs[@]}"
done
exit
}
