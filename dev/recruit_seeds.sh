#!/usr/bin/env bash
# The seed comparison: the benchmark's read sets (dev/recruit_simulate.sh) recruited with each screen
# (closed syncmers, minimizers at the same density, every k-mer) against the GRCh38-only indexes; then
# dev/recruit_seeds_score.py tables certain and recruited fragments per screen, and seconds per screen.
# A run that completed is marked (<output>.done) and skipped on a rerun.
#   recruit_seeds.sh <out dir> [sim dir]
set -euo pipefail

mkdir -p "$1"
OUT="$(cd "$1" && pwd)"; SIM="$(cd "${2:-$HOME/panvar-data/sim}" && pwd)"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
P="$ROOT/build/panvar"
cd "$OUT"
LOCI=(c4 cyp2d gstm acot lpa)
err_of() { case $1 in short|illumina) echo 0.002;; ont) echo 0.03;; hifi) echo 0.001;; esac; }
paired() { [[ $1 == short || $1 == illumina ]]; }
done_or_run() { local marker=$1; shift; [[ -f $marker ]] && return 0; "$@"; touch "$marker"; }
[[ -s times.tsv ]] || echo -e "screen\tset\ttech\tseconds" > times.tsv
run() {   # run <screen> <set> <tech> <reads prefix> <index>...
  local screen=$1 set=$2 tech=$3 pre=$4 in t0 t1; shift 4
  if paired "$tech"; then in=(-1 "${pre}_1.fq.gz" -2 "${pre}_2.fq.gz"); else in=(-1 "$pre.fq.gz"); fi
  local idx=(); for x in "$@"; do idx+=(--locus "$x"); done
  t0=$(python3 -c 'import time; print(time.time())')
  "$P" recruit "${idx[@]}" "${in[@]}" --screen "$screen" --min-sites 2 --error-rate "$(err_of "$tech")" -o "$screen.$set.$tech.tsv" -t 2 -q
  t1=$(python3 -c 'import time; print(time.time())')
  echo -e "$screen\t$set\t$tech\t$(python3 -c "print(round($t1 - $t0, 2))")" >> times.tsv
}
for screen in syncmer minimizer all; do
  for L in "${LOCI[@]}"; do
    for tech in short illumina ont hifi; do
      for src in locus near decoy boundary; do
        done_or_run "$screen.$L.$src.$tech.done" run "$screen" "$L.$src" "$tech" "$SIM/$L/${tech}_$src" "$SIM/$L/g38/$L.plx"
      done
    done
  done
  for tech in short illumina ont hifi; do
    all=(); for L in "${LOCI[@]}"; do all+=("$SIM/$L/g38/$L.plx"); done
    done_or_run "$screen.new.$tech.done" run "$screen" new "$tech" "$SIM/new/${tech}_new" "${all[@]}"
  done
done
python3 "$ROOT/dev/recruit_seeds_score.py" "$OUT" "$SIM" > results.tmp
mv results.tmp results.tsv
column -t results.tsv
