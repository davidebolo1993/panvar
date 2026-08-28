#!/usr/bin/env bash
# genotype_frag_exact_floor.sh - validate a locus's panel floor with true global alignment.
#
#   genotype_frag_exact_floor.sh <locus> [n_donors] [topk]
#
# The panel floor used everywhere else is alignment-derived: minimap2 chooses what to align and every
# base it declines is charged in full. That is a defensible bound but it is not a distance, and at a
# tandem array -- where split and partial alignments are the norm -- it is exactly the number a joint
# model would be judged against. So: shortlist the nearest few candidates with minimap2, compute a
# TRUE Needleman-Wunsch distance for those, and report both. If they differ materially the exact one
# governs, and every excess figure for the locus has to be restated against it.
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/lib/experiment.sh"
REPO="$(cd "$HERE/.." && pwd)"
LOCUS="${1:?usage: genotype_frag_exact_floor.sh <locus> [n_donors] [topk]}"
DONORS_N="${2:-10}"; TOPK="${3:-5}"
BIN="${PANVAR_BIN:-$REPO/build_genotype_review/panvar}"; PY="${PYTHON:-python3}"
OUT="${OUT:-${TMPDIR:-/tmp}/panvar_exactfloor}/$LOCUS"; SEED="${SEED:-42}"
G="$REPO/results/real_data/$LOCUS/bubble/bubble.sorted.gfa"
mkdir -p "$OUT"
exp_init "exact-floor-$LOCUS" "$BIN" "$OUT/provenance.txt"
[[ -d "$OUT/panel_all" ]] || "$PY" "$REPO/scripts/spell_paths.py" -i "$G" -o "$OUT/panel_all" >/dev/null
RES="$OUT/exact_floor.tsv"
printf 'locus\tdonor\tapprox\texact\tapprox_h1\tapprox_h2\texact_h1\texact_h2\tapprox_name1\tapprox_name2\texact_name1\texact_name2\n' > "$RES"

NAMES=(); while IFS= read -r l; do NAMES+=("$l"); done < <("$PY" "$REPO/scripts/spell_paths.py" -i "$G" --list | cut -f1)
DONORS=(); while IFS= read -r l; do DONORS+=("$l"); done < <(
  printf '%s\n' "${NAMES[@]}" | awk -F'#' '{c[$1]++; if(c[$1]==1) first[$1]=$0; else if(c[$1]==2) second[$1]=$0}
    END{for(s in c) if(c[s]>=2) print first[s]"\t"second[s]}' | sort)
ND=${#DONORS[@]}

for ((p=0; p<DONORS_N && p<ND; p++)); do
  d=$(( (SEED + p * 7919) % ND ))
  H1="${DONORS[$d]%%$'\t'*}"; H2="${DONORS[$d]##*$'\t'}"; DONOR="${H1%%#*}"
  rm -rf "$OUT/fa"; "$PY" "$REPO/scripts/spell_paths.py" -i "$G" -o "$OUT/fa" --paths "$H1,$H2" >/dev/null
  T1=$(ls "$OUT"/fa/*.fa | head -1); T2=$(ls "$OUT"/fa/*.fa | tail -1)
  cat $(ls "$OUT"/panel_all/*.fa | grep -v "$DONOR") > "$OUT/panel_loo.fa"
  # Written to a file and rc checked explicitly rather than captured in $(...): under `set -e` a
  # command substitution that fails takes the whole run down before the failure can be reported, and
  # a partially-written line is indistinguishable from a short one.
  rc=0
  EXACT_BINARY="$BIN" EXACT_TOPK="$TOPK" EXACT_TMP="$OUT" \
    "$PY" "$REPO/scripts/genotype_panel_floor.py" "$T1" "$T2" "$OUT/panel_loo.fa" "$DONOR" \
    > "$OUT/line.txt" 2>>"$OUT/provenance.txt" || rc=$?
  if [[ $rc -ne 0 || ! -s "$OUT/line.txt" ]]; then
    echo "  $DONOR: exact floor failed (rc=$rc)" >&2; continue
  fi
  line=$(cat "$OUT/line.txt")
  nf=$(awk -F'\t' 'NR==1{print NF}' "$OUT/line.txt")
  if [[ "${nf:-0}" -lt 11 ]]; then echo "  $DONOR: exact mode produced no result" >&2; continue; fi
  printf '%s\t%s\n' "$LOCUS" "$line" >> "$RES"
  awk -F'\t' -v D="$DONOR" '{printf "  %-10s approx=%-8s exact=%-8s  diff=%+d   (per hap approx %s/%s exact %s/%s)\n",
      D,$2,$3,$3-$2,$4,$5,$6,$7}' <<<"$line"
done

echo
awk -F'\t' 'NR>1{n++; a+=$3; e+=$4} END{if(n){
  printf "n=%d  approximate floor total=%d  EXACT floor total=%d  (%+.1f%%)\n", n,a,e,100.0*(e-a)/(a?a:1)}}' "$RES"
echo "rows: $RES"
