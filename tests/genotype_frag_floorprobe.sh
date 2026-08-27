#!/usr/bin/env bash
# genotype_frag_floorprobe.sh - for the donors that miss the panel floor, WHY?
#
#   genotype_frag_floorprobe.sh <locus> <donor-index-list>
#
# The per-donor excess distribution is bimodal: about half the donors sit at the sequence floor and a
# minority miss it by tens of thousands of edits. A median cannot see that and reported "solved". This
# splits the misses the only way that decides what to fix:
#
#   neither floor haplotype shortlisted   -> candidate generation
#   one shortlisted                       -> candidate generation, partial
#   both shortlisted, floor pair ranked 1 -> projection or spelling, not selection
#   both shortlisted, another pair won    -> the likelihood
set -uo pipefail
LOCUS="${1:?usage: genotype_frag_floorprobe.sh <locus> [n_donors]}"
DONORS_N="${2:-10}"; DEPTH="${3:-30}"; ERR="${4:-0.001}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; REPO="$(cd "$HERE/.." && pwd)"
BIN="${PANVAR_BIN:-$REPO/build_genotype_review/panvar}"; PY="${PYTHON:-python3}"
OUT="${OUT:-${TMPDIR:-/tmp}/panvar_floorprobe}/$LOCUS"; SEED="${SEED:-42}"; MAXHAP="${MAXHAP:-48}"
G="$REPO/results/real_data/$LOCUS/bubble/bubble.sorted.gfa"
PFX="$REPO/results/real_data/$LOCUS/bubble/bubble"
mkdir -p "$OUT"
NAMES=(); while IFS= read -r l; do NAMES+=("$l"); done < <("$PY" "$REPO/scripts/spell_paths.py" -i "$G" --list | cut -f1)
DONORS=(); while IFS= read -r l; do DONORS+=("$l"); done < <(
  printf '%s\n' "${NAMES[@]}" | awk -F'#' '{c[$1]++; if(c[$1]==1) first[$1]=$0; else if(c[$1]==2) second[$1]=$0}
    END{for(s in c) if(c[s]>=2) print first[s]"\t"second[s]}' | sort)
ND=${#DONORS[@]}
[[ ! -d "$OUT/panel_all" ]] && "$PY" "$REPO/scripts/spell_paths.py" -i "$G" -o "$OUT/panel_all" >/dev/null

RES="$OUT/floorprobe.tsv"
printf 'locus\tdonor\tfloor\tcalled\texcess\tfloor_hap1_shortlisted\tfloor_hap2_shortlisted\tfloor_pair_rank\tfloor_pair_delta\tverdict\n' > "$RES"

for ((p=0; p<DONORS_N && p<ND; p++)); do
  d=$(( (SEED + p * 7919) % ND ))
  H1="${DONORS[$d]%%$'\t'*}"; H2="${DONORS[$d]##*$'\t'}"; DONOR="${H1%%#*}"
  rm -rf "$OUT/fa"; "$PY" "$REPO/scripts/spell_paths.py" -i "$G" -o "$OUT/fa" --paths "$H1,$H2" >/dev/null
  T1=$(ls "$OUT"/fa/*.fa | head -1); T2=$(ls "$OUT"/fa/*.fa | tail -1)
  rm -f "$OUT/r_1.fq" "$OUT/r_2.fq"; h=0
  for f in "$OUT"/fa/*.fa; do
    L=$(awk 'NR>1{n+=length($0)} END{print n}' "$f")
    wgsim -N $(( DEPTH * L / 2 / 300 )) -1 150 -2 150 -d 350 -s 50 -e "$ERR" -r 0 -R 0 -X 0 \
      -S $((SEED + p*97 + h)) "$f" "$OUT/a_1.fq" "$OUT/a_2.fq" >/dev/null 2>&1
    cat "$OUT/a_1.fq" >> "$OUT/r_1.fq"; cat "$OUT/a_2.fq" >> "$OUT/r_2.fq"; h=$((h+1))
  done
  gzip -f "$OUT/r_1.fq" "$OUT/r_2.fq"
  cat $(ls "$OUT"/panel_all/*.fa | grep -v "$DONOR") > "$OUT/panel_loo.fa"
  read -r _ FLOOR _ _ B1 B2 < <("$PY" "$REPO/scripts/genotype_panel_floor.py" "$T1" "$T2" "$OUT/panel_loo.fa" f)

  "$BIN" genotype-frag -i "$G" -b "$PFX" -o "$OUT/frag" -R "$OUT/r_1.fq.gz" -R "$OUT/r_2.fq.gz" \
    --haplotype-mode --exclude-haplotypes "$H1,$H2" --max-haplotypes "$MAXHAP" \
    --probe-haplotypes "$B1,$B2" ${FRAG_EXTRA:-} -q >/dev/null 2>&1
  "$BIN" genotype-frag -i "$G" -b "$PFX" -o "$OUT/sp" --exclude-haplotypes "$H1,$H2" \
    --spell-calls "$OUT/frag.hap_blocks.tsv" -q >/dev/null 2>&1
  "$PY" - "$OUT/sp.called.fa" "$OUT" <<'PYS'
import sys
s={}; n=None
for line in open(sys.argv[1]):
    if line[0]=='>': n=line[1:].strip(); s[n]=[]
    else: s[n].append(line.strip())
for i,(k,v) in enumerate(s.items(),1):
    open(f"{sys.argv[2]}/called_{i}.fa","w").write(f">{k}\n{''.join(v)}\n")
PYS
  read -r _ CALLED _ _ _ _ _ < <("$PY" "$REPO/scripts/genotype_pair_sequence_distance.py" \
      "$T1" "$T2" "$OUT/called_1.fa" "$OUT/called_2.fa" c)
  read -r _ _ SL1 RANK _ DELTA _ _ < <(tail -1 "$OUT/frag.hap_probes.tsv")
  # shortlist membership of each floor haplotype separately: "one of two" is a different failure
  s1=$(awk -F'\t' -v n="$B1" 'NR>1 && $1==n{print 1; exit}' "$OUT/frag.hap_scores.tsv"); s1=${s1:-0}
  s2=$(awk -F'\t' -v n="$B2" 'NR>1 && $1==n{print 1; exit}' "$OUT/frag.hap_scores.tsv"); s2=${s2:-0}
  EX=$((CALLED - FLOOR))
  if   [ "$s1$s2" = "00" ]; then V="candidates:neither"
  elif [ "$s1$s2" != "11" ]; then V="candidates:one"
  elif [ "$RANK" = "1" ]; then V="not-selection"
  else V="likelihood"; fi
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$LOCUS" "$DONOR" "$FLOOR" "$CALLED" "$EX" "$s1" "$s2" "$RANK" "$DELTA" "$V" >> "$RES"
  printf "  %-10s floor=%-7s called=%-7s excess=%-7s shortlisted=%s%s rank=%-5s %s\n" \
    "$DONOR" "$FLOOR" "$CALLED" "$EX" "$s1" "$s2" "$RANK" "$V"
done
echo
awk -F'\t' 'NR>1{v[$10]++; if($5>100) big[$10]++} END{
  print "verdict for donors missing the floor by >100 edits:"
  for(k in big) printf "  %-22s %d\n", k, big[k]
  print "all donors:"; for(k in v) printf "  %-22s %d\n", k, v[k]}' "$RES"
echo "rows: $RES"
