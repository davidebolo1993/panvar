#!/usr/bin/env bash
# synthetic_bench.sh - end-to-end validation on a locus whose answers are known by construction.
#
# Real loci cannot validate this pipeline: their truth is whatever the graph encodes, so a wrong call
# and an unrepresentable one look identical. Here every variant is planted, and every design is emitted
# twice as an exact twin -- so holding out a sample's two haplotypes leaves its twins in the panel and
# the expected leave-one-out accuracy is 100%. Anything short of that is a bug in panvar.
#
#   synthetic_bench.sh [out_dir] [n_pairs] [depth] [error]
#
# Env: PANVAR_BIN (default build/panvar), PYTHON, SEED, GEN_EXTRA (flags for make_synthetic_locus.py),
#      GENOTYPE_EXTRA (flags for `panvar genotype`), KEEP=1 to reuse an existing locus.
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/.." && pwd)"
BIN="${PANVAR_BIN:-$REPO/build/panvar}"
PY="${PYTHON:-python3}"
OUT="${1:-${TMPDIR:-/tmp}/panvar_synthetic}"
PAIRS="${2:-8}"
DEPTH="${3:-30}"
ERR="${4:-0.001}"
SEED="${SEED:-7}"

mkdir -p "$OUT"
if [[ "${KEEP:-0}" != "1" || ! -s "$OUT/graph.gfa" ]]; then
  $PY "$REPO/scripts/make_synthetic_locus.py" -o "$OUT" --seed "$SEED" ${GEN_EXTRA:-} || exit 1
fi

echo
echo "== bubble =="
"$BIN" bubble -i "$OUT/graph.gfa" -r ref -o "$OUT/bub" 2>&1 | grep -E "found|Error" || true

# Every planted SV must appear as exactly one bubble, and no bubble may contain another.
$PY - "$OUT" << 'PYEOF'
import csv, sys
out = sys.argv[1]
want = {r["site"]: r for r in csv.DictReader(open(f"{out}/truth.bubbles.tsv"), delimiter="\t")}
got = list(csv.DictReader(open(f"{out}/bub.bubbles.csv")))
print(f"   bubbles found {len(got)}, planted {len(want)}")
ins = {g["bubble_id"]: set(g["inside_nodes"].split(";")) for g in got}
nested = [a for a in ins for b in ins if a != b and ins[a] and ins[a] <= ins[b]]
print(f"   nested bubbles: {len(nested)}" + ("  <-- BUG" if nested else "  ok"))
if len(got) != len(want):
    print("   <-- BUG: bubble count does not match the planted sites")
PYEOF

echo
echo "== call =="
# --classify-ins so a tandem expansion is subtyped DUP rather than left as a bare INS. The other
# route to the same answer is panphorte first, which folds the array into a looped REP node and lets
# call emit a native DUP record; both are checked because the pipeline supports both.
"$BIN" call -i "$OUT/bub.sorted.gfa" -b "$OUT/bub" -r ref -o "$OUT/call" --classify-ins 2>&1 \
  | grep -E "Error" | head -3 || true
if [[ -s "$OUT/call.region.vcf" ]]; then
  echo "   variants called (planted: DEL 400, INS 600, DUP 500-unit ladder, INV 700):"
  grep -v '^#' "$OUT/call.region.vcf" | awk -F'\t' '{
      t=""; if (match($8,/SVTYPE=[A-Z]+/)) t=substr($8,RSTART+7,RLENGTH-7);
      l=""; if (match($8,/SVLEN=-?[0-9]+/)) l=substr($8,RSTART+6,RLENGTH-6);
      s=""; if (match($8,/INS_SUBTYPE=[A-Z]+/)) s=" "substr($8,RSTART+12,RLENGTH-12);
      printf "     POS=%-8s %-4s SVLEN=%-7s%s\n", $2,t,l,s }'
fi

echo
echo "== genotype: leave-zero-out then leave-one-out =="
printf '%-22s %-28s %-28s\n' "sample" "leave-zero-out" "leave-one-out (twins kept)"
NAMES=()
while IFS= read -r l; do NAMES+=("$l"); done < <(awk -F'\t' 'NR>1{print $1}' "$OUT/truth.haplotypes.tsv")
N=${#NAMES[@]}
lz_ok=0; lz_tot=0; lo_ok=0; lo_tot=0; lzb_ok=0; lzb_tot=0; lob_ok=0; lob_tot=0
for ((p=0; p<PAIRS; p++)); do
  i=$(( (p * 2) % N ))
  j=$(( (p * 2 + 3) % N ))
  [[ "$i" == "$j" ]] && j=$(( (j + 1) % N ))
  H1="${NAMES[$i]}"; H2="${NAMES[$j]}"
  F1="$OUT/hap_${H1/\#/_}.fa"; F2="$OUT/hap_${H2/\#/_}.fa"
  rm -f "$OUT/r_1.fq" "$OUT/r_2.fq"
  for f in "$F1" "$F2"; do
    L=$(awk 'NR>1{n+=length($0)} END{print n}' "$f")
    wgsim -N $(( DEPTH * L / 2 / 300 )) -1 150 -2 150 -d 350 -s 50 -e "$ERR" -r 0 -R 0 -X 0 \
      -S $((SEED + p)) "$f" "$OUT/a_1.fq" "$OUT/a_2.fq" >/dev/null 2>&1
    cat "$OUT/a_1.fq" >> "$OUT/r_1.fq"; cat "$OUT/a_2.fq" >> "$OUT/r_2.fq"
  done
  gzip -f "$OUT/r_1.fq" "$OUT/r_2.fq"

  run() {   # $1 = extra flags, $2 = out prefix
    "$BIN" genotype -i "$OUT/bub.sorted.gfa" -b "$OUT/bub" -r ref -o "$2" \
      -R "$OUT/r_1.fq.gz" -R "$OUT/r_2.fq.gz" --truth-haplotypes "$H1,$H2" \
      $1 ${GENOTYPE_EXTRA:-} 2>&1 | grep "truth check"
  }
  parse() { sed -E 's/.*check: ([0-9]+)\/([0-9]+).*bubble blocks ([0-9]+)\/([0-9]+).*/\1 \2 \3 \4/' <<<"$1"; }

  LZ=$(parse "$(run "" "$OUT/gt_lz")")
  LO=$(parse "$(run "--exclude-haplotypes $H1,$H2" "$OUT/gt_lo")")
  read -r a b c d <<<"$LZ"; read -r e f g h <<<"$LO"
  lz_ok=$((lz_ok+a)); lz_tot=$((lz_tot+b)); lzb_ok=$((lzb_ok+c)); lzb_tot=$((lzb_tot+d))
  lo_ok=$((lo_ok+e)); lo_tot=$((lo_tot+f)); lob_ok=$((lob_ok+g)); lob_tot=$((lob_tot+h))
  printf '%-22s %-28s %-28s\n' "$(cut -d'#' -f1 <<<"$H1")+$(cut -d'#' -f1 <<<"$H2")" \
    "blocks $a/$b bubbles $c/$d" "blocks $e/$f bubbles $g/$h"
done

echo
echo "TOTAL leave-zero-out: blocks $lz_ok/$lz_tot, bubbles $lzb_ok/$lzb_tot"
echo "TOTAL leave-one-out : blocks $lo_ok/$lo_tot, bubbles $lob_ok/$lob_tot"
echo "(with exact twins both should be 100% -- the panel can represent the sample exactly)"

# With --twin-divergence the sample's allele may be absent from the panel, so exact match is the wrong
# test: the right answer is the most similar allele available. The identity oracle answers that
# directly -- rank 1 means we chose it.
if [[ -s "$OUT/gt_lo.accuracy.tsv" ]]; then
  awk -F'\t' 'NR>1 && $14 > 0 && $15 > 0 {
      n += 2; if ($14 == 1) k++; if ($15 == 1) k++
      if ($10 > 0) { gap += $10 - ($12 + $13) / 2; g++ }
    } END {
      if (n) printf "ORACLE (last pair): picked the most similar available allele in %d/%d haplotype-blocks (%.0f%%)\n", k, n, 100*k/n
      if (g) printf "                    mean identity shortfall against the best available: %.6f\n", gap/g
    }' "$OUT/gt_lo.accuracy.tsv"
fi
