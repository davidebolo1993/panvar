#!/usr/bin/env bash
# genotype_frag_candidate_loss.sh - how much achievable sequence accuracy did the SHORTLIST lose?
#
#   genotype_frag_candidate_loss.sh <locus> [n_donors]
#
# The earlier measure asked whether one named floor pair survived shortlisting. That repeats a mistake
# this project has already had to retract: floor haplotypes can be tied or sequence-equivalent, so
# probing one representative understates the shortlist. What matters is not whether a particular pair
# survived but whether the shortlist still contains SOME pair as good:
#
#   E*_shortlist = min over pairs in the shortlist of E((a,b), T)
#   E_candidate_loss = E*_shortlist - E*_panel
#
# Only a positive difference proves candidate generation lost achievable accuracy. Both floors use the
# same independent-homologue decomposition, so each is one pass rather than a pair search.
set -uo pipefail
LOCUS="${1:?usage: genotype_frag_candidate_loss.sh <locus> [n_donors]}"
DONORS_N="${2:-10}"; DEPTH="${3:-30}"; ERR="${4:-0.001}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; REPO="$(cd "$HERE/.." && pwd)"
BIN="${PANVAR_BIN:-$REPO/build_genotype_review/panvar}"; PY="${PYTHON:-python3}"
OUT="${OUT:-${TMPDIR:-/tmp}/panvar_candloss}/$LOCUS"; SEED="${SEED:-42}"; MAXHAP="${MAXHAP:-48}"
G="$REPO/results/real_data/$LOCUS/bubble/bubble.sorted.gfa"
PFX="$REPO/results/real_data/$LOCUS/bubble/bubble"
mkdir -p "$OUT"
NAMES=(); while IFS= read -r l; do NAMES+=("$l"); done < <("$PY" "$REPO/scripts/spell_paths.py" -i "$G" --list | cut -f1)
DONORS=(); while IFS= read -r l; do DONORS+=("$l"); done < <(
  printf '%s\n' "${NAMES[@]}" | awk -F'#' '{c[$1]++; if(c[$1]==1) first[$1]=$0; else if(c[$1]==2) second[$1]=$0}
    END{for(s in c) if(c[s]>=2) print first[s]"\t"second[s]}' | sort)
ND=${#DONORS[@]}
[[ ! -d "$OUT/panel_all" ]] && "$PY" "$REPO/scripts/spell_paths.py" -i "$G" -o "$OUT/panel_all" >/dev/null
RES="$OUT/candloss.tsv"
printf 'locus\tdonor\tpanel_floor\tshortlist_floor\tcandidate_loss\tcalled\ttotal_excess\tshortlist_n\n' > "$RES"

for ((p=0; p<DONORS_N && p<ND; p++)); do
  d=$(( (SEED + p * 7919) % ND ))
  H1="${DONORS[$d]%%$'\t'*}"; H2="${DONORS[$d]##*$'\t'}"; DONOR="${H1%%#*}"
  rm -rf "$OUT/fa"; "$PY" "$REPO/scripts/spell_paths.py" -i "$G" -o "$OUT/fa" --paths "$H1,$H2" >/dev/null
  T1=$(ls "$OUT"/fa/*.fa | head -1); T2=$(ls "$OUT"/fa/*.fa | tail -1)
  rm -f "$OUT/r_1.fq" "$OUT/r_2.fq"; h=0
  for f in "$OUT"/fa/*.fa; do
    L=$(awk '!/^>/{n+=length($0)} END{print n}' "$f")
    wgsim -N $(( DEPTH * L / 2 / 300 )) -1 150 -2 150 -d 350 -s 50 -e "$ERR" -r 0 -R 0 -X 0 \
      -S $((SEED + p*97 + h)) "$f" "$OUT/a_1.fq" "$OUT/a_2.fq" >/dev/null 2>&1
    cat "$OUT/a_1.fq" >> "$OUT/r_1.fq"; cat "$OUT/a_2.fq" >> "$OUT/r_2.fq"; h=$((h+1))
  done
  gzip -f "$OUT/r_1.fq" "$OUT/r_2.fq"
  cat $(ls "$OUT"/panel_all/*.fa | grep -v "$DONOR") > "$OUT/panel_loo.fa"
  read -r _ PFLOOR _ _ _ _ < <("$PY" "$REPO/scripts/genotype_panel_floor.py" "$T1" "$T2" "$OUT/panel_loo.fa" f)

  "$BIN" genotype-frag -i "$G" -b "$PFX" -o "$OUT/frag" -R "$OUT/r_1.fq.gz" -R "$OUT/r_2.fq.gz" \
    --haplotype-mode --exclude-haplotypes "$H1,$H2" --max-haplotypes "$MAXHAP" -q >/dev/null 2>&1
  # exactly the haplotypes the caller could choose among
  awk -F'\t' 'NR>1{print $1}' "$OUT/frag.hap_scores.tsv" > "$OUT/sl.txt"
  SLN=$(wc -l < "$OUT/sl.txt" | tr -d ' ')
  : > "$OUT/shortlist.fa"
  while IFS= read -r n; do
    f=$(grep -l "^>${n}$" "$OUT"/panel_all/*.fa 2>/dev/null | head -1)
    [[ -n "$f" ]] && cat "$f" >> "$OUT/shortlist.fa"
  done < "$OUT/sl.txt"
  if [[ -s "$OUT/shortlist.fa" ]]; then
    read -r _ SFLOOR _ _ _ _ < <("$PY" "$REPO/scripts/genotype_panel_floor.py" "$T1" "$T2" "$OUT/shortlist.fa" s)
  else SFLOOR="NA"; fi

  "$BIN" genotype-frag -i "$G" -b "$PFX" -o "$OUT/sp" --exclude-haplotypes "$H1,$H2" \
    --spell-calls "$OUT/frag.hap_blocks.tsv" -q >/dev/null 2>&1
  awk '/^>called_1/{f=1;next} /^>called_2/{f=0} f' "$OUT/sp.called.fa" | sed '1i\
>c1' > "$OUT/c1.fa"
  awk '/^>called_2/{f=1;next} f' "$OUT/sp.called.fa" | sed '1i\
>c2' > "$OUT/c2.fa"
  read -r _ CALLED _ _ _ _ _ _ < <("$PY" "$REPO/scripts/genotype_pair_sequence_distance.py" \
      "$T1" "$T2" "$OUT/c1.fa" "$OUT/c2.fa" c)
  CL=$([[ "$SFLOOR" == "NA" ]] && echo NA || echo $((SFLOOR - PFLOOR)))
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$LOCUS" "$DONOR" "$PFLOOR" "$SFLOOR" "$CL" "$CALLED" \
    "$((CALLED - PFLOOR))" "$SLN" >> "$RES"
  printf "  %-10s panel_floor=%-7s shortlist_floor=%-7s candidate_loss=%-8s total_excess=%-8s\n" \
    "$DONOR" "$PFLOOR" "$SFLOOR" "$CL" "$((CALLED - PFLOOR))"
done
echo
awk -F'\t' 'NR>1 && $5!="NA"{n++; cl+=$5; ex+=$7; if($5>100)b++}
  END{if(n) printf "n=%d  total candidate loss=%d  total excess=%d  candidate loss is %.0f%% of excess; donors losing >100: %d\n",
      n, cl, ex, 100*cl/(ex?ex:1), b+0}' "$RES"
echo "rows: $RES"
