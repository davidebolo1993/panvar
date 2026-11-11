#!/usr/bin/env bash
# The recruitment benchmark on the read sets of dev/recruit_simulate.sh: every set is recruited once and
# dev/recruit_score.py scores it against its truth.
#   main arm  each locus's sets against its index on GRCh38 + CHM13 (short-read locus pairs also mate 1
#             alone), and the GRCh38 genome set against the five loci together
#   g38 arm   the locus sets against the index on GRCh38 only, and the CHM13 genome set, which that
#             background never saw, against the five loci together
#   wg arm    when whole-genome indexes exist (in <locus>/$WG_DIR, default wg): every locus-side set
#             against them
#   wgnew     when the whole-CHM13 sets exist: against the five loci together, in every arm
# Seconds and CPU seconds per set are in times.tsv.
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
TECHS=(short illumina ont hifi)
WG_DIR="${WG_DIR:-wg}"
WG=1
for L in "${LOCI[@]}"; do [[ -f "$L/$WG_DIR/$L.plx" ]] || WG=0; done

err_of() { case $1 in short|illumina) echo 0.002;; ont) echo 0.03;; hifi) echo 0.001;; esac; }
paired() { [[ $1 == short || $1 == illumina ]]; }
sets=()
echo -e "arm\tlocus\ttech\tsource\tseconds\tcpu_seconds" > "$OUT/times.tsv"
recruit() {   # recruit <arm> <label> <tech> <source> <table> <reads1> [reads2] -- <index>...
  local arm=$1 label=$2 tech=$3 source=$4 table=$5 r1=$6 r2=$7; shift 7
  local idx=() t0 cpu
  for x in "$@"; do idx+=(--locus "$x"); done
  local args=(recruit "${idx[@]}" -1 "$r1" -o "$table" --error-rate "$(err_of "$tech")" -t 2 -q)
  [[ -n $r2 ]] && args+=(-2 "$r2")
  t0=$(date +%s)
  cpu=$( { /usr/bin/time -p "$P" "${args[@]}" > /dev/null; } 2>&1 | awk '/^user|^sys/ {s += $2} END {print s}')
  echo -e "$arm\t$label\t$tech\t$source\t$(( $(date +%s) - t0 ))\t$cpu" >> "$OUT/times.tsv"
}
one() {   # one <arm> <locus> <tech> <source> <reads prefix> <index...>: recruit a set and queue its scoring
  local arm=$1 L=$2 tech=$3 src=$4 pre=$5; shift 5
  local t="$OUT/$arm.$L.$tech.$src.tsv" r1 r2=""
  if paired "$tech"; then r1="${pre}_1.fq.gz"; r2="${pre}_2.fq.gz"; else r1="$pre.fq.gz"; fi
  [[ $src == locus1 ]] && r2=""
  recruit "$arm" "$L" "$tech" "$src" "$t" "$r1" "$r2" "$@"
  sets+=("$arm" "$L" "$tech" "$src" "$r1" "$t")
}
arms=(main g38); [[ $WG == 1 ]] && arms+=(wg)
for L in "${LOCI[@]}"; do
  echo "== $L $(date +%H:%M:%S)"
  for tech in "${TECHS[@]}"; do
    for src in locus near decoy boundary; do
      one main "$L" "$tech" "$src" "$L/${tech}_$src" "$L/$L.plx"
      [[ $src == locus ]] && one g38 "$L" "$tech" locus "$L/${tech}_locus" "$L/g38/$L.plx"
      [[ $WG == 1 ]] && one wg "$L" "$tech" "$src" "$L/${tech}_$src" "$L/$WG_DIR/$L.plx"
    done
    paired "$tech" && one main "$L" "$tech" locus1 "$L/${tech}_locus" "$L/$L.plx"
  done
done
negatives=()
for arm in "${arms[@]}"; do
  case $arm in main) negatives+=("main genome");; g38) negatives+=("g38 new");; esac
  [[ -s wgnew/illumina_wgnew_1.fq.gz ]] && negatives+=("$arm wgnew")
done
for pair in "${negatives[@]}"; do
  read -r arm g <<< "$pair"
  case $arm in main) dir=.;; g38) dir=g38;; wg) dir=$WG_DIR;; esac
  idx=(); for L in "${LOCI[@]}"; do if [[ $arm == main ]]; then idx+=("$L/$L.plx"); else idx+=("$L/$dir/$L.plx"); fi; done
  for tech in "${TECHS[@]}"; do
    [[ $g == wgnew && $tech == short ]] && continue
    echo "== $arm $g $tech $(date +%H:%M:%S)"
    t="$OUT/$arm.$g.$tech.tsv"
    if paired "$tech"; then recruit $arm all "$tech" $g "$t" "$g/${tech}_${g}_1.fq.gz" "$g/${tech}_${g}_2.fq.gz" "${idx[@]}"; r1="$g/${tech}_${g}_1.fq.gz"
    else recruit $arm all "$tech" $g "$t" "$g/${tech}_$g.fq.gz" "" "${idx[@]}"; r1="$g/${tech}_$g.fq.gz"; fi
    for L in "${LOCI[@]}"; do sets+=($arm "$L" "$tech" $g "$r1" "$t"); done
  done
done
python3 "$ROOT/dev/recruit_score.py" --sim "$SIM" --out "$OUT/summary.tsv" "${sets[@]}"
echo "== done $(date +%H:%M:%S)"
