#!/usr/bin/env bash
# Locityper's recruitment beside panvar recruit on the benchmark's read sets (dev/recruit_simulate.sh),
# with the same background content: GRCh38 (the proof-of-concept chromosomes) as the reference whose
# k-mer counts Locityper uses (jellyfish, k 25, as its documentation recommends), and as panvar's
# GRCh38-only background (the g38 indexes). Locityper's targets are each locus's haplotypes, spelled from
# the same held-out graphs. Locus, near, decoy and boundary sets are recruited to their own locus; the
# CHM13 genome set (a genome neither background holds) to the five loci together. Then
# dev/locityper_score.py compares them. Sets already recruited are kept (a rerun resumes).
#   locityper_compare.sh <out dir> <jellyfish counts> <targets dir with <locus>.fa> [sim dir]
# LOCITYPER=<locityper executable>; panvar from build/.
set -euo pipefail

OUT="$1"; JF="$2"; TARGETS="$3"; SIM="${4:-$HOME/panvar-data/sim}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
P="$ROOT/build/panvar"
LT="${LOCITYPER:-locityper}"
mkdir -p "$OUT"
LOCI=(c4 cyp2d gstm acot lpa)
preset() { case $1 in short|illumina) echo illumina;; ont) echo ont;; hifi) echo hifi;; esac; }
err_of() { case $1 in short|illumina) echo 0.002;; ont) echo 0.03;; hifi) echo 0.001;; esac; }
paired() { [[ $1 == short || $1 == illumina ]]; }
[[ -s "$OUT/times.tsv" ]] || echo -e "tool\tset\ttech\tseconds\tmax_rss_kb" > "$OUT/times.tsv"
timed() {   # timed <tool> <set> <tech> <command...>
  local tool=$1 set=$2 tech=$3; shift 3
  /usr/bin/time -l "$@" > /dev/null 2> "$OUT/time.tmp" || { cat "$OUT/time.tmp" >&2; exit 1; }
  echo -e "$tool\t$set\t$tech\t$(awk '/ real/ {print $1}' "$OUT/time.tmp")\t$(awk '/maximum resident/ {print int($1 / 1024)}' "$OUT/time.tmp")" >> "$OUT/times.tsv"
}
reads_of() {   # reads_of <tech> <prefix>: the read arguments, for locityper and for panvar
  if paired "$1"; then LT_IN=(-i "$2_1.fq.gz" "$2_2.fq.gz"); PV_IN=(-1 "$2_1.fq.gz" -2 "$2_2.fq.gz")
  else LT_IN=(-i "$2.fq.gz"); PV_IN=(-1 "$2.fq.gz"); fi
}
for L in "${LOCI[@]}"; do
  for tech in short illumina ont hifi; do
    for src in locus near decoy boundary; do
      reads_of "$tech" "$SIM/$L/${tech}_$src"
      [[ -s "$OUT/lt.$L.$tech.$src.fq" || -f "$OUT/lt.$L.$tech.$src.done" ]] ||
      timed locityper "$L.$src" "$tech" "$LT" recruit "${LT_IN[@]}" -s "$TARGETS/$L.fa" -j "$JF" -x "$(preset "$tech")" \
        -o "$OUT/lt.$L.$tech.$src.fq" -@ 4
      touch "$OUT/lt.$L.$tech.$src.done"
      [[ $src == locus ]] && continue   # the g38 arm's locus tables are the benchmark's
      [[ -s "$OUT/pv.$L.$tech.$src.tsv" ]] ||
      timed panvar "$L.$src" "$tech" "$P" recruit --locus "$SIM/$L/g38/$L.plx" "${PV_IN[@]}" --error-rate "$(err_of "$tech")" \
        -o "$OUT/pv.$L.$tech.$src.tsv" -t 4 -q
    done
  done
done
for tech in short illumina ont hifi; do
  reads_of "$tech" "$SIM/new/${tech}_new"
  timed locityper new "$tech" "$LT" recruit "${LT_IN[@]}" -s "$TARGETS/{}.fa" -j "$JF" -x "$(preset "$tech")" -o "$OUT/lt.{}.$tech.new.fq" -@ 4
  idx=(); for L in "${LOCI[@]}"; do idx+=(--locus "$SIM/$L/g38/$L.plx"); done
  timed panvar new "$tech" "$P" recruit "${idx[@]}" "${PV_IN[@]}" --error-rate "$(err_of "$tech")" -o "$OUT/pv.all.$tech.new.tsv" -t 4 -q
done
rm -f "$OUT/time.tmp"
python3 "$ROOT/dev/locityper_score.py" "$OUT" "$SIM" > "$OUT/comparison.tsv"
column -t "$OUT/comparison.tsv"
