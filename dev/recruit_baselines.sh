#!/usr/bin/env bash
# Baselines beside panvar recruit on the benchmark's read sets (dev/recruit_simulate.sh):
#   bbduk     BBDuk with every 31-mer of the locus haplotypes as bait (a pair is taken when a mate holds one)
#   mirabait  MIRA's mirabait, the same bait at k 31 (when mirabait is on PATH; MIRA has no macOS build)
# Locus, near, decoy and boundary sets are recruited to their own locus; the CHM13 genome set to the five
# loci together. Then dev/recruit_baselines_score.py tables them beside panvar's certain and recruited
# fragments (and Locityper's, from a dev/locityper_compare.sh run). A step that completed is marked
# (<output>.done) and skipped on a rerun. Seconds and peak memory per run in times.tsv.
#   recruit_baselines.sh <out dir> <targets dir with <locus>.fa> [sim dir] [locityper_compare out dir]
# BBDUK, MIRABAIT: the executables (default: on PATH; BBDuk needs Java).
set -euo pipefail

abs() { echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"; }
OUT="$1"; TARGETS="$(abs "$2")"; SIM="$(abs "${3:-$HOME/panvar-data/sim}")"; LT_RUN="${4:-}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BBDUK="${BBDUK:-bbduk.sh}"; MB="${MIRABAIT:-mirabait}"
mkdir -p "$OUT"
cd "$OUT"
LOCI=(c4 cyp2d gstm acot lpa)
paired() { [[ $1 == short || $1 == illumina ]]; }
[[ -s times.tsv ]] || echo -e "tool\tset\ttech\tseconds\tmax_rss_kb" > times.tsv
if /usr/bin/time -v true 2> /dev/null; then TIME=(/usr/bin/time -v); else TIME=(/usr/bin/time -l); fi
timed() {   # timed <tool> <set> <tech> <command...>
  local tool=$1 set=$2 tech=$3 t0 t1 rss; shift 3
  t0=$(python3 -c 'import time; print(time.time())')
  "${TIME[@]}" "$@" > /dev/null 2> time.tmp || { cat time.tmp >&2; exit 1; }
  t1=$(python3 -c 'import time; print(time.time())')
  rss=$(awk '/Maximum resident set size/ {print $NF} /maximum resident set size/ {print int($1 / 1024)}' time.tmp)
  echo -e "$tool\t$set\t$tech\t$(python3 -c "print(round($t1 - $t0, 2))")\t$rss" >> times.tsv
}
done_or_run() { local marker=$1; shift; [[ -f $marker ]] && return 0; "$@" && touch "$marker"; }
have_mirabait=0; command -v "$MB" > /dev/null && have_mirabait=1

all_targets() { [[ -s all.targets.fa ]] || cat $(for L in "${LOCI[@]}"; do echo "$TARGETS/$L.fa"; done) > all.targets.fa; }
run_set() {   # run_set <locus, or all> <tech> <reads prefix> <label>
  local L=$1 tech=$2 pre=$3 label=$4 in=() mbin=()
  if paired "$tech"; then
    in=(in="${pre}_1.fq.gz" in2="${pre}_2.fq.gz"); mbin=("${pre}_1.fq.gz" "${pre}_2.fq.gz")
  else
    in=(in="${pre}.fq.gz"); mbin=("${pre}.fq.gz")
  fi
  local ref=$TARGETS/$L.fa; [[ $L == all ]] && { all_targets; ref=all.targets.fa; }
  done_or_run bbduk.$label.done timed bbduk "$label" "$tech" "$BBDUK" "${in[@]}" ref="$ref" k=31 hdist=0 outm=bbduk.$label.fq.gz \
    threads=2 -Xmx6g overwrite=t
  if [[ $have_mirabait == 1 ]]; then
    # paired files are baited together (-p): a pair is kept when either mate has a bait k-mer
    local p=(); paired "$tech" && p=(-p)
    done_or_run mirabait.$label.done timed mirabait "$label" "$tech" "$MB" -b "$ref" -k 31 -n 1 "${p[@]}" -o mirabait.$label "${mbin[@]}"
  fi
}
for L in "${LOCI[@]}"; do
  for tech in short illumina ont hifi; do
    for src in locus near decoy boundary; do
      run_set "$L" "$tech" "$SIM/$L/${tech}_$src" "$L.$tech.$src"
    done
  done
done
for tech in short illumina ont hifi; do
  run_set all "$tech" "$SIM/new/${tech}_new" "all.$tech.new"
done
rm -f time.tmp
python3 "$ROOT/dev/recruit_baselines_score.py" "$OUT" "$SIM" ${LT_RUN:+"$LT_RUN"} > results.tsv
column -t results.tsv
