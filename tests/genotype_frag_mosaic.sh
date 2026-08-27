#!/usr/bin/env bash
# genotype_frag_mosaic.sh - how much could a mosaic model recover, at most?
#
#   genotype_frag_mosaic.sh <locus> [n_donors]
#
# H_mosaic = complete-pair floor - free block-mosaic floor, per donor. Needs no reads: it is a
# property of the panel and the truth. Reports the switch count the free optimum uses, because a
# gain bought with one switch per boundary is not a gain any evidence could support.
set -uo pipefail
LOCUS="${1:?usage: genotype_frag_mosaic.sh <locus> [n_donors]}"
DONORS_N="${2:-10}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; REPO="$(cd "$HERE/.." && pwd)"
BIN="${PANVAR_BIN:-$REPO/build_genotype_review/panvar}"; PY="${PYTHON:-python3}"
OUT="${OUT:-${TMPDIR:-/tmp}/panvar_mosaic}/$LOCUS"; SEED="${SEED:-42}"
G="$REPO/results/real_data/$LOCUS/bubble/bubble.sorted.gfa"
PFX="$REPO/results/real_data/$LOCUS/bubble/bubble"
mkdir -p "$OUT"
NAMES=(); while IFS= read -r l; do NAMES+=("$l"); done < <("$PY" "$REPO/scripts/spell_paths.py" -i "$G" --list | cut -f1)
DONORS=(); while IFS= read -r l; do DONORS+=("$l"); done < <(
  printf '%s\n' "${NAMES[@]}" | awk -F'#' '{c[$1]++; if(c[$1]==1) first[$1]=$0; else if(c[$1]==2) second[$1]=$0}
    END{for(s in c) if(c[s]>=2) print first[s]"\t"second[s]}' | sort)
ND=${#DONORS[@]}
RES="$OUT/mosaic.tsv"
printf 'locus\tdonor\tcomplete\tfree\th_mosaic\tswitches\tpen10\tpen100\n' > "$RES"
for ((p=0; p<DONORS_N && p<ND; p++)); do
  d=$(( (SEED + p * 7919) % ND ))
  H1="${DONORS[$d]%%$'\t'*}"; H2="${DONORS[$d]##*$'\t'}"; DONOR="${H1%%#*}"
  "$BIN" genotype-frag -i "$G" -b "$PFX" -o "$OUT/m" --mosaic-floor \
    --truth-haplotypes "$H1,$H2" --exclude-haplotypes "$H1,$H2" \
    --switch-penalties 10,100 -t 8 -q >/dev/null 2>&1
  [[ -s "$OUT/m.mosaic_floor.tsv" ]] || { echo "  $DONOR skipped"; continue; }
  read -r C F SW P10 P100 < <(awk -F'\t' 'NR>1{
      if($2=="complete") c+=$4;
      else if($2=="free"){f+=$4; sw+=$5}
      else if($2=="penalised" && $3==10) p10+=$4;
      else if($2=="penalised" && $3==100) p100+=$4}
    END{print c+0, f+0, sw+0, p10+0, p100+0}' "$OUT/m.mosaic_floor.tsv")
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$LOCUS" "$DONOR" "$C" "$F" "$((C-F))" "$SW" "$P10" "$P100" >> "$RES"
  printf "  %-10s complete=%-8s free=%-8s H_mosaic=%-8s switches=%-4s pen10=%-8s\n" "$DONOR" "$C" "$F" "$((C-F))" "$SW" "$P10"
done
echo
awk -F'\t' 'NR>1{n++; c+=$3; f+=$4; h+=$5; s+=$6; p+=$7; if($5>100)big++}
  END{if(n){printf "n=%d  complete=%d  free=%d  H_mosaic=%d (%.0f%% of complete)  switches=%d  pen10=%d\n",
      n,c,f,h,100*h/(c?c:1),s,p; printf "donors where a mosaic could recover >100 edits: %d/%d\n", big+0, n}}' "$RES"
echo "rows: $RES"
