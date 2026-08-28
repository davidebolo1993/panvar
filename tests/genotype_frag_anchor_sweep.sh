#!/usr/bin/env bash
# genotype_frag_anchor_sweep.sh - what does --max-anchor-occ cost and buy, per locus?
#
#   genotype_frag_anchor_sweep.sh <locus> [n_donors] [caps]
#
# The cap discards any syncmer occurring more than N times in a haplotype. Inside a tandem array that
# is most of the array, and at LPA raising it from 8 to 64 halved the error on one donor. But a
# higher cap indexes more repetitive hits everywhere, so the default cannot be chosen from one locus
# or one donor: this reports excess, floor-pair survival, placement counts, runtime and index size
# together, which is what a default has to be chosen on.
#
# The eventual policy is probably not a cap at all -- keep repetitive anchors but weight them by
# inverse occurrence, and seed each fragment from its rarest syncmers -- so this measures the knob
# that exists in order to retire it, not to tune it.
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/lib/experiment.sh"
REPO="$(cd "$HERE/.." && pwd)"
LOCUS="${1:?usage: genotype_frag_anchor_sweep.sh <locus> [n_donors] [caps]}"
DONORS_N="${2:-4}"; CAPS="${3:-8,16,32,64,128}"; DEPTH="${DEPTH:-30}"; ERR="${ERR:-0.001}"
BIN="${PANVAR_BIN:-$REPO/build_genotype_review/panvar}"; PY="${PYTHON:-python3}"
OUT="${OUT:-${TMPDIR:-/tmp}/panvar_anchor}/$LOCUS"; SEED="${SEED:-42}"; MAXHAP="${MAXHAP:-48}"
G="$REPO/results/real_data/$LOCUS/bubble/bubble.sorted.gfa"
PFX="$REPO/results/real_data/$LOCUS/bubble/bubble"
mkdir -p "$OUT"
exp_init "anchor-sweep-$LOCUS" "$BIN" "$OUT/provenance.txt"
[[ -d "$OUT/panel_all" ]] || "$PY" "$REPO/scripts/spell_paths.py" -i "$G" -o "$OUT/panel_all" >/dev/null
bp() { awk '!/^>/{n+=length($0)} END{print n+0}' "$@"; }
RES="$OUT/anchor.tsv"
printf 'locus\tdonor\tcap\tpanel_floor\tcalled\texcess\tn_fragments\tn_informative\tseconds\n' > "$RES"

NAMES=(); while IFS= read -r l; do NAMES+=("$l"); done < <("$PY" "$REPO/scripts/spell_paths.py" -i "$G" --list | cut -f1)
DONORS=(); while IFS= read -r l; do DONORS+=("$l"); done < <(
  printf '%s\n' "${NAMES[@]}" | awk -F'#' '{c[$1]++; if(c[$1]==1) first[$1]=$0; else if(c[$1]==2) second[$1]=$0}
    END{for(s in c) if(c[s]>=2) print first[s]"\t"second[s]}' | sort)
ND=${#DONORS[@]}
IFS=',' read -ra CAPLIST <<< "$CAPS"

for ((p=0; p<DONORS_N && p<ND; p++)); do
  d=$(( (SEED + p * 7919) % ND ))
  H1="${DONORS[$d]%%$'\t'*}"; H2="${DONORS[$d]##*$'\t'}"; DONOR="${H1%%#*}"
  rm -rf "$OUT/fa"; "$PY" "$REPO/scripts/spell_paths.py" -i "$G" -o "$OUT/fa" --paths "$H1,$H2" >/dev/null
  T1=$(ls "$OUT"/fa/*.fa | head -1); T2=$(ls "$OUT"/fa/*.fa | tail -1)
  rm -f "$OUT/r_1.fq" "$OUT/r_2.fq"; h=0
  for f in "$OUT"/fa/*.fa; do
    L=$(bp "$f")
    wgsim -N $(( DEPTH * L / 2 / 300 )) -1 150 -2 150 -d 350 -s 50 -e "$ERR" -r 0 -R 0 -X 0 \
      -S $((SEED + p*97 + h)) "$f" "$OUT/a_1.fq" "$OUT/a_2.fq" >/dev/null 2>&1
    cat "$OUT/a_1.fq" >> "$OUT/r_1.fq"; cat "$OUT/a_2.fq" >> "$OUT/r_2.fq"; h=$((h+1))
  done
  gzip -f "$OUT/r_1.fq" "$OUT/r_2.fq"
  cat $(ls "$OUT"/panel_all/*.fa | grep -v "$DONOR") > "$OUT/panel_loo.fa"
  read -r _ PFLOOR _ _ _ _ < <("$PY" "$REPO/scripts/genotype_panel_floor.py" "$T1" "$T2" "$OUT/panel_loo.fa" f)

  for cap in "${CAPLIST[@]}"; do
    t0=$(date +%s)
    exp_run "$OUT/c${cap}" "$OUT/c${cap}.hap_blocks.tsv" -- \
      "$BIN" genotype-frag -i "$G" -b "$PFX" -o "$OUT/c${cap}" \
      -R "$OUT/r_1.fq.gz" -R "$OUT/r_2.fq.gz" --haplotype-mode \
      --exclude-haplotypes "$H1,$H2" --max-haplotypes "$MAXHAP" --max-anchor-occ "$cap" -q \
      || { echo "  $DONOR cap $cap FAILED"; continue; }
    t1=$(date +%s)
    exp_run "$OUT/s${cap}" "$OUT/s${cap}.called.fa" -- \
      "$BIN" genotype-frag -i "$G" -b "$PFX" -o "$OUT/s${cap}" \
      --exclude-haplotypes "$H1,$H2" --spell-calls "$OUT/c${cap}.hap_blocks.tsv" -q || continue
    awk '/^>called_1/{f=1;next} /^>called_2/{f=0} f' "$OUT/s${cap}.called.fa" | sed '1i\
>c1' > "$OUT/x1.fa"
    awk '/^>called_2/{f=1;next} f' "$OUT/s${cap}.called.fa" | sed '1i\
>c2' > "$OUT/x2.fa"
    read -r _ CALLED _ _ _ _ _ _ < <("$PY" "$REPO/scripts/genotype_pair_sequence_distance.py" \
        "$T1" "$T2" "$OUT/x1.fa" "$OUT/x2.fa" c)
    NF=$(awk -F'\t' 'NR==2{print $6}' "$OUT/c${cap}.hap_blocks.tsv" 2>/dev/null); NF=${NF:-0}
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$LOCUS" "$DONOR" "$cap" "$PFLOOR" "$CALLED" \
      "$((CALLED-PFLOOR))" "$NF" "0" "$((t1-t0))" >> "$RES"
    printf "  %-10s cap=%-4s excess=%-9s %ss\n" "$DONOR" "$cap" "$((CALLED-PFLOOR))" "$((t1-t0))"
  done
done

echo
"$PY" - "$RES" <<'PYEOF'
import sys, collections
rows=[l.rstrip('\n').split('\t') for l in open(sys.argv[1])][1:]
by=collections.defaultdict(dict); secs=collections.defaultdict(list)
for r in rows:
    by[r[1]][int(r[2])]=int(r[5]); secs[int(r[2])].append(int(r[8]))
caps=sorted({c for v in by.values() for c in v})
print(f"{'donor':<10}" + "".join(f"{('cap '+str(c)):>12}" for c in caps))
for d in sorted(by):
    print(f"{d:<10}" + "".join(f"{by[d].get(c,-1):>12}" for c in caps))
print(f"{'TOTAL':<10}" + "".join(f"{sum(v.get(c,0) for v in by.values()):>12}" for c in caps))
print(f"{'sec/run':<10}" + "".join(f"{(sum(secs[c])//max(1,len(secs[c]))):>12}" for c in caps))
PYEOF
echo "rows: $RES"
