#!/usr/bin/env bash
# genotype_frag_gap.sh - decompose the prototype's distance from its own non-mosaic ceiling.
#
#   genotype_frag_gap.sh <locus> [n_donors] [depth] [error]
#
# The prototype reaches 63% at cyp2d6 where the best fixed pair of complete panel haplotypes reaches
# 85%. That 22-point gap is a MIXTURE of three failures and no algorithm change is worth making until
# they are separated:
#
#   1. candidate generation - the ceiling pair was never shortlisted, so it was never scored;
#   2. placement           - it was shortlisted but the reads did not anchor onto it;
#   3. likelihood          - it was shortlisted, reads placed, and it still lost.
#
# Only the third is an argument for replacing best-placement scoring with a marginalised one. The
# ceiling searches the whole panel while the prototype shortlists 48, so (1) is a live confound and
# has to be excluded before (3) can be claimed.
#
# Env: PANVAR_BIN, SEED, OUT, MAXHAP (shortlist size, default 48; set to a large number for the
#      full-panel arm of the shortlist ladder)
set -uo pipefail

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  sed -n '2,20p' "$0" | sed 's/^# \{0,1\}//'
  exit 0
fi

LOCUS="${1:?usage: genotype_frag_gap.sh <locus> [n_donors] [depth] [error]}"
DONORS_N="${2:-10}"
DEPTH="${3:-30}"
ERR="${4:-0.001}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/.." && pwd)"
BIN="${PANVAR_BIN:-$REPO/build_genotype_review/panvar}"
PY="${PYTHON:-python3}"
OUT="${OUT:-${TMPDIR:-/tmp}/panvar_frag_gap}/$LOCUS"
SEED="${SEED:-42}"
MAXHAP="${MAXHAP:-48}"

G="$REPO/results/real_data/$LOCUS/bubble/bubble.sorted.gfa"
PFX="$REPO/results/real_data/$LOCUS/bubble/bubble"
[[ -s "$G" ]] || { echo "no bubble-stage graph for $LOCUS"; exit 1; }
mkdir -p "$OUT"
REF="$("$PY" "$REPO/scripts/spell_paths.py" -i "$G" --list | cut -f1 | grep -i grch38 | head -1)"
[[ -n "$REF" ]] || REF="$("$PY" "$REPO/scripts/spell_paths.py" -i "$G" --list | cut -f1 | head -1)"

NAMES=(); while IFS= read -r l; do NAMES+=("$l"); done < <("$PY" "$REPO/scripts/spell_paths.py" -i "$G" --list | cut -f1)
DONORS=(); while IFS= read -r l; do DONORS+=("$l"); done < <(
  printf '%s\n' "${NAMES[@]}" | awk -F'#' '{c[$1]++; if(c[$1]==1) first[$1]=$0; else if(c[$1]==2) second[$1]=$0}
    END{for(s in c) if(c[s]>=2) print first[s]"\t"second[s]}' | sort)
ND=${#DONORS[@]}

GAP="$OUT/gap.tsv"
printf 'locus\tdonor\tmaxhap\tn_optima\tceiling_exact\tn_representable\tprototype_exact\tin_shortlist\tceiling_rank\tceiling_delta\tplaced1\tplaced2\tbest_placed\n' > "$GAP"

echo "locus $LOCUS: $DONORS_N donors at ${DEPTH}x, leave-one-out, shortlist $MAXHAP"

for ((p=0; p<DONORS_N && p<ND; p++)); do
  d=$(( (SEED + p * 7919) % ND ))
  H1="${DONORS[$d]%%$'\t'*}"; H2="${DONORS[$d]##*$'\t'}"
  DONOR="${H1%%#*}"
  rm -rf "$OUT/fa"; "$PY" "$REPO/scripts/spell_paths.py" -i "$G" -o "$OUT/fa" --paths "$H1,$H2" >/dev/null
  rm -f "$OUT/r_1.fq" "$OUT/r_2.fq"; hidx=0
  for f in "$OUT"/fa/*.fa; do
    L=$(awk 'NR>1{n+=length($0)} END{print n}' "$f")
    wgsim -N $(( DEPTH * L / 2 / 300 )) -1 150 -2 150 -d 350 -s 50 -e "$ERR" -r 0 -R 0 -X 0 \
      -S $((SEED + p*97 + hidx)) "$f" "$OUT/a_1.fq" "$OUT/a_2.fq" >/dev/null 2>&1
    cat "$OUT/a_1.fq" >> "$OUT/r_1.fq"; cat "$OUT/a_2.fq" >> "$OUT/r_2.fq"; hidx=$((hidx+1))
  done
  gzip -f "$OUT/r_1.fq" "$OUT/r_2.fq"

  "$BIN" genotype -i "$G" --bubble-prefix-in "$PFX" -r "$REF" -o "$OUT/prod" \
    -R "$OUT/r_1.fq.gz" -R "$OUT/r_2.fq.gz" --truth-haplotypes "$H1,$H2" \
    --exclude-haplotypes "$H1,$H2" --dump-haplotype-alleles "$OUT/alleles.tsv" -q >/dev/null 2>&1
  [[ -s "$OUT/alleles.tsv" ]] || { echo "  $DONOR: production produced no allele dump, skipped"; continue; }
  "$PY" "$REPO/scripts/genotype_pair_ceiling.py" "$OUT/alleles.tsv" "$OUT/prod.genotypes.tsv" \
    "$LOCUS" 1 "$DONOR" "$OUT/optimum.tsv" > /dev/null
  read -r _ C1 C2 NOPT CEX CREP _ < <(sed -n 2p "$OUT/optimum.tsv")
  # Probe EVERY tied optimum and keep the best-ranked. Probing one of several understates the
  # likelihood whenever another optimum ranks better, and the ceiling here has a median of four.
  PROBES=(); while IFS=$'\t' read -r _ u v _; do PROBES+=(--probe-haplotypes "$u,$v"); done < <(tail -n +2 "$OUT/optimum.tsv")

  "$BIN" genotype-frag -i "$G" -b "$PFX" -o "$OUT/frag" -R "$OUT/r_1.fq.gz" -R "$OUT/r_2.fq.gz" \
    --haplotype-mode --truth-haplotypes "$H1,$H2" --exclude-haplotypes "$H1,$H2" \
    --max-haplotypes "$MAXHAP" "${PROBES[@]}" -q >/dev/null 2>&1
  [[ -s "$OUT/frag.hap_probes.tsv" ]] || { echo "  $DONOR: no probe output, skipped"; continue; }
  read -r INSL RANK DELTA PL1 PL2 < <(awk -F'\t' 'NR>1{
      if ($3==1) { s=1; if (b=="" || ($4>0 && $4<b)) { b=$4; d=$6; p1=$7; p2=$8 } } }
    END{ printf "%d\t%s\t%s\t%s\t%s\n", s+0, (b==""?-2:b), (d==""?0:d), p1+0, p2+0 }' "$OUT/frag.hap_probes.tsv")
  PEX=$(awk -F'\t' 'NR>1&&$10==1{r++; if($11==1)e++} END{print e+0}' "$OUT/frag.hap_blocks.tsv")
  BESTPL=$(awk -F'\t' 'NR>1{if($4>m)m=$4} END{print m+0}' "$OUT/frag.hap_scores.tsv")
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$LOCUS" "$DONOR" "$MAXHAP" "$NOPT" "$CEX" "$CREP" "$PEX" "$INSL" "$RANK" "$DELTA" "$PL1" "$PL2" "$BESTPL" >> "$GAP"
  echo "  $DONOR: ceiling $CEX/$CREP ($NOPT tied optima), prototype $PEX; ceiling pair shortlisted=$INSL rank=$RANK delta=$DELTA"
done

echo
"$PY" - "$GAP" <<'PYEOF'
import sys
rows=[l.rstrip('\n').split('\t') for l in open(sys.argv[1])]
h={n:i for i,n in enumerate(rows[0])}; rows=rows[1:]
if not rows: sys.exit("no rows")
n=len(rows)
not_short=[r for r in rows if r[h['in_shortlist']]=='0']
short=[r for r in rows if r[h['in_shortlist']]=='1']
won=[r for r in short if r[h['ceiling_rank']]=='1']
lost=[r for r in short if r[h['ceiling_rank']] not in ('1','-2')]
print(f"donors: {n}")
print(f"  ceiling pair NEVER SHORTLISTED  : {len(not_short):2d}   -> candidate generation")
print(f"  shortlisted and ranked 1        : {len(won):2d}   -> likelihood already agrees with the ceiling")
print(f"  shortlisted and LOST            : {len(lost):2d}   -> likelihood, the only case the placement rewrite targets")
if lost:
    ds=sorted(float(r[h['ceiling_delta']]) for r in lost)
    print(f"    their score deltas: median {ds[len(ds)//2]:.0f}, range {ds[0]:.0f} .. {ds[-1]:.0f}")
    print(f"    placement is NOT the failure where placed fragments match the best haplotype's:")
    for r in lost:
        pl=max(int(r[h['placed1']]),int(r[h['placed2']])); bp=int(r[h['best_placed']])
        print(f"      {r[h['donor']]:<10} ceiling {r[h['ceiling_exact']]}/{r[h['n_representable']]} vs prototype {r[h['prototype_exact']]},"
              f" rank {r[h['ceiling_rank']]}, placed {pl} of best {bp} ({100*pl//max(1,bp)}%)")
ties=[int(r[h['n_optima']]) for r in rows]
print(f"  tied optima per donor: median {sorted(ties)[len(ties)//2]}, max {max(ties)}"
      f"  (a ceiling with many tied optima is a weaker target than it looks)")
PYEOF
echo
echo "per-donor rows: $GAP"
