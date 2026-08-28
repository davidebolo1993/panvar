#!/usr/bin/env bash
# genotype_frag_multiseed.sh - is each mechanism systematic, or one simulation draw?
#
#   genotype_frag_multiseed.sh <locus> <donor>[,<donor>...] [n_seeds]
#
# Every mechanism in the record so far rests on one seed. The frozen diagnostic donors carry the
# argument for the joint model, so each needs to be shown stable before anything is built on it.
# Reports, per donor per seed: panel floor, shortlist floor, called distance, excess, and the length
# decomposition -- so a donor whose excess swings between seeds is visible as such rather than
# averaged into a mechanism.
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/lib/experiment.sh"
REPO="$(cd "$HERE/.." && pwd)"
LOCUS="${1:?usage: genotype_frag_multiseed.sh <locus> <donor[,donor...]> [n_seeds]}"
DONOR_LIST="${2:?}"; NSEED="${3:-5}"; DEPTH="${DEPTH:-30}"; ERR="${ERR:-0.001}"
BIN="${PANVAR_BIN:-$REPO/build_genotype_review/panvar}"; PY="${PYTHON:-python3}"
OUT="${OUT:-${TMPDIR:-/tmp}/panvar_multiseed}/$LOCUS"; MAXHAP="${MAXHAP:-48}"
G="$REPO/results/real_data/$LOCUS/bubble/bubble.sorted.gfa"
PFX="$REPO/results/real_data/$LOCUS/bubble/bubble"
mkdir -p "$OUT"
exp_init "multiseed-$LOCUS" "$BIN" "$OUT/provenance.txt"
[[ -d "$OUT/panel_all" ]] || "$PY" "$REPO/scripts/spell_paths.py" -i "$G" -o "$OUT/panel_all" >/dev/null
bp() { awk '!/^>/{n+=length($0)} END{print n+0}' "$@"; }
RES="$OUT/multiseed.tsv"
printf 'locus\tdonor\tseed\tpanel_floor\tshortlist_floor\tcalled\texcess\tcand_loss\ttruth_bp\tcalled_bp\tlen_err\n' > "$RES"

IFS=',' read -ra DONORS <<< "$DONOR_LIST"
for DONOR in "${DONORS[@]}"; do
  H1=$("$PY" "$REPO/scripts/spell_paths.py" -i "$G" --list | cut -f1 | grep "^${DONOR}#1" | head -1)
  H2=$("$PY" "$REPO/scripts/spell_paths.py" -i "$G" --list | cut -f1 | grep "^${DONOR}#2" | head -1)
  [[ -n "$H1" && -n "$H2" ]] || { echo "  $DONOR: not two haplotypes in $LOCUS, skipped" >&2; continue; }
  rm -rf "$OUT/fa"; "$PY" "$REPO/scripts/spell_paths.py" -i "$G" -o "$OUT/fa" --paths "$H1,$H2" >/dev/null
  T1=$(ls "$OUT"/fa/*.fa | head -1); T2=$(ls "$OUT"/fa/*.fa | tail -1)
  TBP=$(( $(bp "$T1") + $(bp "$T2") ))
  cat $(ls "$OUT"/panel_all/*.fa | grep -v "$DONOR") > "$OUT/panel_loo.fa"
  read -r _ PFLOOR _ _ _ _ < <("$PY" "$REPO/scripts/genotype_panel_floor.py" "$T1" "$T2" "$OUT/panel_loo.fa" f)

  for ((s=1; s<=NSEED; s++)); do
    SEED=$((1000 + s * 7919))
    rm -f "$OUT/r_1.fq" "$OUT/r_2.fq"; h=0
    for f in "$OUT"/fa/*.fa; do
      L=$(bp "$f")
      wgsim -N $(( DEPTH * L / 2 / 300 )) -1 150 -2 150 -d 350 -s 50 -e "$ERR" -r 0 -R 0 -X 0 \
        -S $((SEED + h)) "$f" "$OUT/a_1.fq" "$OUT/a_2.fq" >/dev/null 2>&1
      cat "$OUT/a_1.fq" >> "$OUT/r_1.fq"; cat "$OUT/a_2.fq" >> "$OUT/r_2.fq"; h=$((h+1))
    done
    gzip -f "$OUT/r_1.fq" "$OUT/r_2.fq"

    exp_run "$OUT/f${s}" "$OUT/f${s}.hap_blocks.tsv" "$OUT/f${s}.hap_scores.tsv" -- \
      "$BIN" genotype-frag -i "$G" -b "$PFX" -o "$OUT/f${s}" \
      -R "$OUT/r_1.fq.gz" -R "$OUT/r_2.fq.gz" --haplotype-mode \
      --exclude-haplotypes "$H1,$H2" --max-haplotypes "$MAXHAP" -q || { echo "  $DONOR seed $s FAILED"; continue; }

    : > "$OUT/shortlist.fa"
    while IFS= read -r n; do
      src=$(grep -l "^>${n}$" "$OUT"/panel_all/*.fa 2>/dev/null | head -1)
      [[ -n "$src" ]] && cat "$src" >> "$OUT/shortlist.fa"
    done < <(awk -F'\t' 'NR>1{print $1}' "$OUT/f${s}.hap_scores.tsv")
    read -r _ SFLOOR _ _ _ _ < <("$PY" "$REPO/scripts/genotype_panel_floor.py" "$T1" "$T2" "$OUT/shortlist.fa" s)

    exp_run "$OUT/sp${s}" "$OUT/sp${s}.called.fa" -- \
      "$BIN" genotype-frag -i "$G" -b "$PFX" -o "$OUT/sp${s}" \
      --exclude-haplotypes "$H1,$H2" --spell-calls "$OUT/f${s}.hap_blocks.tsv" -q || continue
    awk '/^>called_1/{f=1;next} /^>called_2/{f=0} f' "$OUT/sp${s}.called.fa" | sed '1i\
>c1' > "$OUT/c1.fa"
    awk '/^>called_2/{f=1;next} f' "$OUT/sp${s}.called.fa" | sed '1i\
>c2' > "$OUT/c2.fa"
    read -r _ CALLED _ _ _ _ _ _ < <("$PY" "$REPO/scripts/genotype_pair_sequence_distance.py" \
        "$T1" "$T2" "$OUT/c1.fa" "$OUT/c2.fa" c)
    CBP=$(( $(bp "$OUT/c1.fa") + $(bp "$OUT/c2.fa") ))
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$LOCUS" "$DONOR" "$SEED" "$PFLOOR" \
      "$SFLOOR" "$CALLED" "$((CALLED-PFLOOR))" "$((SFLOOR-PFLOOR))" "$TBP" "$CBP" \
      "$(( CBP > TBP ? CBP-TBP : TBP-CBP ))" >> "$RES"
    printf "  %-10s seed %d  floor=%-7s excess=%-8s cand_loss=%-7s len_err=%s\n" \
      "$DONOR" "$s" "$PFLOOR" "$((CALLED-PFLOOR))" "$((SFLOOR-PFLOOR))" "$(( CBP > TBP ? CBP-TBP : TBP-CBP ))"
  done
done

echo
"$PY" - "$RES" <<'PYEOF'
import sys, collections, statistics
rows=[l.rstrip('\n').split('\t') for l in open(sys.argv[1])][1:]
by=collections.defaultdict(list)
for r in rows: by[(r[0],r[1])].append(int(r[6]))
print(f"{'locus':<9} {'donor':<10} {'n':>2} {'min':>9} {'median':>9} {'max':>9}   verdict")
for (loc,d),v in sorted(by.items()):
    if not v: continue
    med=statistics.median(v)
    spread = (max(v)-min(v))
    verdict = "systematic" if min(v) > 100 else ("never fails" if max(v) <= 100 else "SEED-DEPENDENT")
    print(f"{loc:<9} {d:<10} {len(v):>2} {min(v):>9} {int(med):>9} {max(v):>9}   {verdict}")
PYEOF
echo "rows: $RES"
