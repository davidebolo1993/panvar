#!/usr/bin/env bash
# synthetic_bench.sh - end-to-end validation on a locus whose answers are known by construction.
#
# Real loci cannot validate this pipeline: their truth is whatever the graph encodes, so a wrong call
# and an unrepresentable one look identical. Here every variant is planted.
#
# TWIN_DIVERGENCE defaults to 40 SNPs, NOT to 0. With exact twins, holding out a sample's two
# haplotypes leaves an identical haplotype in the panel, so leave-one-out has a perfectly representable
# answer and is not an off-panel test at all. Measured: with exact twins every case scores 36/36 blocks
# and 16/16 bubbles at 30x, 10x AND 5x -- 2.5x per haplotype -- because picking an identical panel
# haplotype barely needs coverage. The ladder was not saturated by generous depth; the task was
# trivial, and no reduction in depth could make it discriminate.
#
# So leave-ZERO-out is the implementation check and should be 100%: the sample's own haplotypes are in
# the panel and anything short of exact is a bug. Leave-ONE-out is the off-panel measurement and should
# NOT be 100%; with 40 SNPs of divergence it runs 11-16 of 16 bubbles and degrades with depth, which is
# the behaviour a real off-panel test has. Set TWIN_DIVERGENCE=0 to reproduce figures measured under
# the old default, and read them knowing what they were.
#
#   synthetic_bench.sh [out_dir] [n_pairs] [depth] [error]
#
# Env: TWIN_DIVERGENCE (default 40; 0 restores exact twins and makes leave-one-out trivial),
#      PANVAR_BIN (default build/panvar), PYTHON, SEED, GEN_EXTRA (flags for make_synthetic_locus.py),
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
  # Applied unless the caller already pinned it, so an explicit GEN_EXTRA still wins.
  TWIN_FLAG=""
  case " ${GEN_EXTRA:-} " in
    *" --twin-divergence "*) ;;
    *) TWIN_FLAG="--twin-divergence ${TWIN_DIVERGENCE:-40}" ;;
  esac
  $PY "$REPO/scripts/make_synthetic_locus.py" -o "$OUT" --seed "$SEED" $TWIN_FLAG ${GEN_EXTRA:-} || exit 1
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
printf '%-22s %-28s %-28s\n' "sample" "leave-zero-out" "leave-one-out (off-panel)"
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
  hidx=0
  for f in "$F1" "$F2"; do
    L=$(awk 'NR>1{n+=length($0)} END{print n}' "$f")
    # A distinct seed per homologue. Sharing one makes their coverage fluctuations identical, which is
    # the one thing a diploid simulation must not do: it removes exactly the per-haplotype depth noise
    # that decides heterozygous from homozygous. genotype_sim.sh was fixed for this long ago; this
    # harness was not (defect V14), so every absolute ladder figure predates the fix.
    hseed=$((SEED + p*97 + hidx)); hidx=$((hidx+1))
    wgsim -N $(( DEPTH * L / 2 / 300 )) -1 150 -2 150 -d 350 -s 50 -e "$ERR" -r 0 -R 0 -X 0 \
      -S "$hseed" "$f" "$OUT/a_1.fq" "$OUT/a_2.fq" >/dev/null 2>&1
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
  # Accumulate before the next pair overwrites it; the summary below reads every pair, not the last.
  if [[ -s "$OUT/gt_lo.accuracy.tsv" ]]; then
    [[ -s "$OUT/gt_lo.accuracy.all.tsv" ]] || head -1 "$OUT/gt_lo.accuracy.tsv" > "$OUT/gt_lo.accuracy.all.tsv"
    tail -n +2 "$OUT/gt_lo.accuracy.tsv" >> "$OUT/gt_lo.accuracy.all.tsv"
  fi
  read -r a b c d <<<"$LZ"; read -r e f g h <<<"$LO"
  lz_ok=$((lz_ok+a)); lz_tot=$((lz_tot+b)); lzb_ok=$((lzb_ok+c)); lzb_tot=$((lzb_tot+d))
  lo_ok=$((lo_ok+e)); lo_tot=$((lo_tot+f)); lob_ok=$((lob_ok+g)); lob_tot=$((lob_tot+h))
  printf '%-22s %-28s %-28s\n' "$(cut -d'#' -f1 <<<"$H1")+$(cut -d'#' -f1 <<<"$H2")" \
    "blocks $a/$b bubbles $c/$d" "blocks $e/$f bubbles $g/$h"
done

echo
echo "TOTAL leave-zero-out: blocks $lz_ok/$lz_tot, bubbles $lzb_ok/$lzb_tot"
echo "TOTAL leave-one-out : blocks $lo_ok/$lo_tot, bubbles $lob_ok/$lob_tot"
echo "(leave-zero-out should be 100%: the sample's haplotypes are in the panel, so anything short is a"
echo " bug. leave-one-out is NOT guaranteed to be: with diverged twins it is an off-panel measurement,"
echo " and a better genotyper may legitimately reach 100% there.)"
echo "(--twin-divergence perturbs SNP sites only, so twins still SHARE every structural bubble allele."
echo " The exact score therefore measures known bubble alleles on an unseen background; the ORACLE line"
echo " below is the nearest-available-allele measurement over every block.)"

# With --twin-divergence the sample's allele may be absent from the panel, so exact match is the wrong
# test: the right answer is the most similar allele available. The identity oracle answers that
# directly -- rank 1 means we chose it.
# A fully mosaic sample: a different donor in every region of the chain, so the backbones that carry a
# marker-poor bubble are themselves from different donors and nothing along the chain is stable. This is
# the case a linkage-based model should struggle with, so it is worth checking explicitly rather than
# assuming. Both mosaics are held out, so the panel cannot contain either.
if [[ -s "$OUT/hap_mosaic0.fa" && -s "$OUT/hap_mosaic1.fa" ]]; then
  echo
  echo "== fully mosaic sample (donor switches at every region; both held out) =="
  rm -f "$OUT/m_1.fq" "$OUT/m_2.fq"
  for h in mosaic0 mosaic1; do
    L=$(awk 'NR>1{n+=length($0)} END{print n}' "$OUT/hap_$h.fa")
    wgsim -N $(( DEPTH * L / 2 / 300 )) -1 150 -2 150 -d 350 -s 50 -e "$ERR" -r 0 -R 0 -X 0 \
      -S "$SEED" "$OUT/hap_$h.fa" "$OUT/a_1.fq" "$OUT/a_2.fq" >/dev/null 2>&1
    cat "$OUT/a_1.fq" >> "$OUT/m_1.fq"; cat "$OUT/a_2.fq" >> "$OUT/m_2.fq"
  done
  gzip -f "$OUT/m_1.fq" "$OUT/m_2.fq"
  "$BIN" genotype -i "$OUT/bub.sorted.gfa" -b "$OUT/bub" -r ref -o "$OUT/gt_mo" \
    -R "$OUT/m_1.fq.gz" -R "$OUT/m_2.fq.gz" --truth-haplotypes "mosaic0,mosaic1" \
    --exclude-haplotypes "mosaic0,mosaic1" ${GENOTYPE_EXTRA:-} 2>&1 | grep -E "truth check" | sed 's/^/  /'
  awk -F'\t' 'NR>1 && $2!="bubble"{n++; if ($10!=prev) sw++; prev=$10}
              END{printf "  haplotype1 switched at %d of %d non-bubble blocks (a mosaic should switch often)\n", sw-1, n}' \
    "$OUT/gt_mo.genotypes.tsv"
fi

# Every held-out pair, not the last one. `gt_lo.accuracy.tsv` is rewritten at the same prefix on each
# iteration, so reading it after the loop reported ONE pair under a header that said otherwise.
#
# The two numbers answer different questions and both matter. The exact-allele score above counts only
# blocks whose truth allele is in the panel, so it measures genotyping of representable alleles. This
# one asks whether, over EVERY block, the call is the most similar allele available -- which is the
# behaviour the real-data problem is about, and which can be much worse than the exact score suggests.
if [[ -s "$OUT/gt_lo.accuracy.all.tsv" ]]; then
  awk -F'\t' -v PAIRS="$PAIRS" 'FNR>1 && $14 > 0 && $15 > 0 {
      n += 2; if ($14 == 1) k++; if ($15 == 1) k++
      if ($10 > 0) { gap += $10 - ($12 + $13) / 2; g++ }
    } END {
      if (n) printf "ORACLE (all %d pairs): picked the most similar available allele in %d/%d haplotype-blocks (%.0f%%)\n", PAIRS, k, n, 100*k/n
      if (g) printf "                       mean identity shortfall against the best available: %.6f\n", gap/g
    }' "$OUT/gt_lo.accuracy.all.tsv"
fi
