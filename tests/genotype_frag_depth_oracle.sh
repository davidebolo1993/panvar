#!/usr/bin/env bash
# genotype_frag_depth_oracle.sh - can TOTAL DOSAGE rescue the catastrophic tail, before any EM?
#
#   genotype_frag_depth_oracle.sh <locus> [n_donors]
#
# Pre-registered, three arms, run before implementing a windowed read-assignment model:
#
#   A  alignment only                      the current default
#   B  + total-fragment Poisson            lambda fitted ONCE outside any candidate
#   C  + truth total length                DIAGNOSTIC ONLY: no caller can know this
#
# The question C answers is the one that decides whether to build the joint model at all: if even
# PERFECT total-length knowledge does not close the catastrophic donors, their failure is not total
# dosage, and a depth model aimed at total dosage will not fix them. B answers whether an externally
# calibrated estimate is enough on its own.
#
# Declared before running: A < B < C in total excess if the tail is a dosage problem. If C does not
# beat A materially, the tail is homologue allocation or composition and this line closes.
set -uo pipefail
LOCUS="${1:?usage: genotype_frag_depth_oracle.sh <locus> [n_donors]}"
DONORS_N="${2:-10}"; DEPTH="${3:-30}"; ERR="${4:-0.001}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; REPO="$(cd "$HERE/.." && pwd)"
BIN="${PANVAR_BIN:-$REPO/build_genotype_review/panvar}"; PY="${PYTHON:-python3}"
OUT="${OUT:-${TMPDIR:-/tmp}/panvar_depth_oracle}/$LOCUS"; SEED="${SEED:-42}"; MAXHAP="${MAXHAP:-48}"
G="$REPO/results/real_data/$LOCUS/bubble/bubble.sorted.gfa"
PFX="$REPO/results/real_data/$LOCUS/bubble/bubble"
mkdir -p "$OUT"
NAMES=(); while IFS= read -r l; do NAMES+=("$l"); done < <("$PY" "$REPO/scripts/spell_paths.py" -i "$G" --list | cut -f1)
DONORS=(); while IFS= read -r l; do DONORS+=("$l"); done < <(
  printf '%s\n' "${NAMES[@]}" | awk -F'#' '{c[$1]++; if(c[$1]==1) first[$1]=$0; else if(c[$1]==2) second[$1]=$0}
    END{for(s in c) if(c[s]>=2) print first[s]"\t"second[s]}' | sort)
ND=${#DONORS[@]}
[[ ! -d "$OUT/panel_all" ]] && "$PY" "$REPO/scripts/spell_paths.py" -i "$G" -o "$OUT/panel_all" >/dev/null
RES="$OUT/oracle.tsv"
printf 'locus\tdonor\tarm\tfloor\tcalled\texcess\n' > "$RES"
bp() { awk '!/^>/{n+=length($0)} END{print n+0}' "$@"; }

for ((p=0; p<DONORS_N && p<ND; p++)); do
  d=$(( (SEED + p * 7919) % ND ))
  H1="${DONORS[$d]%%$'\t'*}"; H2="${DONORS[$d]##*$'\t'}"; DONOR="${H1%%#*}"
  rm -rf "$OUT/fa"; "$PY" "$REPO/scripts/spell_paths.py" -i "$G" -o "$OUT/fa" --paths "$H1,$H2" >/dev/null
  T1=$(ls "$OUT"/fa/*.fa | head -1); T2=$(ls "$OUT"/fa/*.fa | tail -1)
  TBP=$(( $(bp "$T1") + $(bp "$T2") ))
  rm -f "$OUT/r_1.fq" "$OUT/r_2.fq"; h=0
  for f in "$OUT"/fa/*.fa; do
    L=$(bp "$f")
    wgsim -N $(( DEPTH * L / 2 / 300 )) -1 150 -2 150 -d 350 -s 50 -e "$ERR" -r 0 -R 0 -X 0 \
      -S $((SEED + p*97 + h)) "$f" "$OUT/a_1.fq" "$OUT/a_2.fq" >/dev/null 2>&1
    cat "$OUT/a_1.fq" >> "$OUT/r_1.fq"; cat "$OUT/a_2.fq" >> "$OUT/r_2.fq"; h=$((h+1))
  done
  gzip -f "$OUT/r_1.fq" "$OUT/r_2.fq"
  cat $(ls "$OUT"/panel_all/*.fa | grep -v "$DONOR") > "$OUT/panel_loo.fa"
  read -r _ FLOOR _ _ _ _ < <("$PY" "$REPO/scripts/genotype_panel_floor.py" "$T1" "$T2" "$OUT/panel_loo.fa" f)

  arm() {  # name, extra flags
    "$BIN" genotype-frag -i "$G" -b "$PFX" -o "$OUT/a" -R "$OUT/r_1.fq.gz" -R "$OUT/r_2.fq.gz" \
      --haplotype-mode --exclude-haplotypes "$H1,$H2" --max-haplotypes "$MAXHAP" ${2:-} -q >/dev/null 2>&1
    "$BIN" genotype-frag -i "$G" -b "$PFX" -o "$OUT/s" --exclude-haplotypes "$H1,$H2" \
      --spell-calls "$OUT/a.hap_blocks.tsv" -q >/dev/null 2>&1
    awk '/^>called_1/{f=1;next} /^>called_2/{f=0} f' "$OUT/s.called.fa" > "$OUT/c1.fa"
    awk '/^>called_2/{f=1;next} f' "$OUT/s.called.fa" > "$OUT/c2.fa"
    sed -i.bak '1i\
>c
' "$OUT/c1.fa" 2>/dev/null; sed -i.bak '1i\
>c
' "$OUT/c2.fa" 2>/dev/null
    read -r _ CALLED _ _ _ _ _ _ < <("$PY" "$REPO/scripts/genotype_pair_sequence_distance.py" \
        "$T1" "$T2" "$OUT/c1.fa" "$OUT/c2.fa" c)
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$LOCUS" "$DONOR" "$1" "$FLOOR" "$CALLED" "$((CALLED-FLOOR))" >> "$RES"
    printf "%s=%s " "$1" "$((CALLED-FLOOR))"
  }
  printf "  %-10s floor=%-7s " "$DONOR" "$FLOOR"
  arm A ""
  arm B "--total-depth"
  arm C "--truth-total-bp $TBP"
  echo
done

echo
"$PY" - "$RES" <<'PYEOF'
import sys, collections
rows=[l.rstrip('\n').split('\t') for l in open(sys.argv[1])][1:]
by=collections.defaultdict(dict)
for r in rows: by[r[1]][r[2]]=int(r[5])
arms=['A','B','C']
ds=[d for d,v in by.items() if all(a in v for a in arms)]
print(f"{'donor':<10}" + "".join(f"{a:>12}" for a in arms))
for d in sorted(ds, key=lambda x:-by[x]['A']):
    print(f"{d:<10}" + "".join(f"{by[d][a]:>12}" for a in arms))
print(f"{'TOTAL':<10}" + "".join(f"{sum(by[d][a] for d in ds):>12}" for a in arms))
print(f"{'MEDIAN':<10}" + "".join(f"{sorted(by[d][a] for d in ds)[len(ds)//2]:>12}" for a in arms))
big=[d for d in ds if by[d]['A']>100]
print(f"\ndonors with excess >100 under alignment-only: {len(big)}")
for d in big:
    print(f"  {d:<10} A={by[d]['A']:<9} B={by[d]['B']:<9} C={by[d]['C']:<9}"
          f"  {'C rescues' if by[d]['C'] < by[d]['A']/2 else 'C does NOT rescue'}")
PYEOF
echo "rows: $RES"
