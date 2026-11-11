#!/usr/bin/env bash
# The recruitment benchmark on the read sets of dev/recruit_simulate.sh: every set is recruited once and
# dev/recruit_score.py scores every decision rule on the same evidence.
#   main arm  each locus's sets against its index on GRCh38 + CHM13 (short locus reads also mate 1 alone),
#             and the GRCh38 genome set against the five loci together
#   g38 arm   the locus sets against the index on GRCh38 only, and the CHM13 genome set, which that
#             background never saw, against the five loci together
# Seconds per set are in times.tsv.
#   recruit_benchmark.sh [sim dir]
set -euo pipefail

SIM="${1:-$HOME/panvar-data/sim}"
[[ -z "$SIM" ]] && SIM="$HOME/panvar-data/sim"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
P="$ROOT/build/panvar"
OUT="$SIM/benchmark"
mkdir -p "$OUT"
cd "$SIM"
LOCI=(c4 cyp2d gstm acot lpa)

count() { echo $(( $(zcat < "$1" | wc -l) / 4 )); }
err_of() { case $1 in short) echo 0.002;; ont) echo 0.03;; hifi) echo 0.001;; esac; }
sets=()
echo -e "arm\tlocus\ttech\tsource\tfragments\tseconds" > "$OUT/times.tsv"
recruit() {   # recruit <arm> <label> <tech> <source> <table> <reads1> [reads2] -- <index>...
  local arm=$1 label=$2 tech=$3 source=$4 table=$5 r1=$6 r2=$7; shift 7
  local idx=() n t0
  for x in "$@"; do idx+=(--locus "$x"); done
  n=$(count "$r1")
  t0=$(date +%s)
  if [[ -n $r2 ]]; then "$P" recruit "${idx[@]}" -1 "$r1" -2 "$r2" -o "$table" --error-rate "$(err_of "$tech")" -t 2 -q
  else "$P" recruit "${idx[@]}" -1 "$r1" -o "$table" --error-rate "$(err_of "$tech")" -t 2 -q; fi
  echo -e "$arm\t$label\t$tech\t$source\t$n\t$(( $(date +%s) - t0 ))" >> "$OUT/times.tsv"
  LAST_N=$n
}
for L in "${LOCI[@]}"; do
  echo "== $L $(date +%H:%M:%S)"
  for src in locus near decoy boundary; do
    for tech in short ont hifi; do
      t="$OUT/main.$L.$tech.$src.tsv"
      if [[ $tech == short ]]; then recruit main "$L" short "$src" "$t" "$L/short_${src}_1.fq.gz" "$L/short_${src}_2.fq.gz" "$L/$L.plx"
      else recruit main "$L" "$tech" "$src" "$t" "$L/${tech}_$src.fq.gz" "" "$L/$L.plx"; fi
      sets+=(main "$L" "$tech" "$src" "$LAST_N" "$t")
      if [[ $src == locus ]]; then
        t="$OUT/g38.$L.$tech.locus.tsv"
        if [[ $tech == short ]]; then recruit g38 "$L" short locus "$t" "$L/short_locus_1.fq.gz" "$L/short_locus_2.fq.gz" "$L/g38/$L.plx"
        else recruit g38 "$L" "$tech" locus "$t" "$L/${tech}_locus.fq.gz" "" "$L/g38/$L.plx"; fi
        sets+=(g38 "$L" "$tech" locus "$LAST_N" "$t")
      fi
    done
  done
  t="$OUT/main.$L.short.locus1.tsv"
  recruit main "$L" short locus1 "$t" "$L/short_locus_1.fq.gz" "" "$L/$L.plx"
  sets+=(main "$L" short locus1 "$LAST_N" "$t")
done
for arm in main g38; do
  g=genome; [[ $arm == g38 ]] && g=new
  idx=(); for L in "${LOCI[@]}"; do if [[ $arm == main ]]; then idx+=("$L/$L.plx"); else idx+=("$L/g38/$L.plx"); fi; done
  for tech in short ont hifi; do
    echo "== $arm $g $tech $(date +%H:%M:%S)"
    t="$OUT/$arm.$g.$tech.tsv"
    if [[ $tech == short ]]; then recruit $arm all short $g "$t" "$g/short_${g}_1.fq.gz" "$g/short_${g}_2.fq.gz" "${idx[@]}"
    else recruit $arm all "$tech" $g "$t" "$g/${tech}_$g.fq.gz" "" "${idx[@]}"; fi
    for L in "${LOCI[@]}"; do sets+=($arm "$L" "$tech" $g "$LAST_N" "$t"); done
  done
done
python3 "$ROOT/dev/recruit_score.py" --out "$OUT/summary.tsv" "${sets[@]}"
echo "== done $(date +%H:%M:%S)"
