#!/usr/bin/env bash
# Per-gene copy number from certain reads (dev/recruit_gene_cn.py) on the hard cases, against the truth
# tables of tests/real_data:
#   gstm1_null  a sample without GSTM1 on either haplotype, held out of the graph: GSTM1 must read 0
#   rccx        the C4 benchmark sample (HG00733): C4A, C4B, CYP21A2, TNXB, all inside the locus
#   cyp2d       the CYP2D benchmark sample (HG00733), whose first haplotype carries an extra CYP2D unit
#               the truth table names no gene (a hybrid): CYP2D6 reads above 2
#   lpa         KIV-2 copies of the LPA benchmark sample (HG002)
#   every gene of the five benchmark loci
# Short pairs (wgsim, 15x per haplotype for the held-out sample; the benchmark's own for the others),
# recruited against the GRCh38-only background (the g38 arm). Writes <out>/cn.tsv.
#   recruit_hard_cases.sh <out dir> [sim dir] [GRCh38-only store] [null sample]
set -euo pipefail

abs() { echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"; }
OUT="$1"; SIM="$(abs "${2:-$HOME/panvar-data/sim}")"; BG="$(abs "${3:-$SIM/bg_g38.bgk}")"; NULL="${4:-HG00097}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
P="$ROOT/build/panvar"
RD="$ROOT/tests/real_data"
CN=(python3 "$ROOT/dev/recruit_gene_cn.py")
mkdir -p "$OUT"
cd "$OUT"
# a step's outputs are trusted only once its marker is written, after it succeeded
done_or_run() { local marker=$1; shift; [[ -f $marker ]] && return 0; "$@"; touch "$marker"; }
src() { python3 "$ROOT/dev/recruit_sources.py" "$1" --gfa "$2" --sample "$3" --out "$4.tmp"; mv "$4.tmp" "$4"; }

# GSTM1-null sample: held out, its index, its reads
mkdir -p gstm1_null
[[ -s gstm1_null/gstm.gfa ]] || src held-out "$RD/paralogs/gstm.gfa.gz" "$NULL" gstm1_null/gstm.gfa
[[ -s gstm1_null/locus.fa ]] || src locus "$RD/paralogs/gstm.gfa.gz" "$NULL" gstm1_null/locus.fa
[[ -f gstm1_null/gstm.plx ]] || "$P" recruit --locus gstm1_null/gstm.gfa --background "$BG" -r grch38 --index-only --save-index gstm1_null -t 2 -q
simulate() {
  local len
  len=$(awk '!/^>/ {n += length($0)} END {print n}' gstm1_null/locus.fa)
  wgsim -N $((len * 15 / 300)) -1 150 -2 150 -d 400 -s 50 -e 0.002 -r 0 -R 0 -X 0 -S 11 gstm1_null/locus.fa gstm1_null/reads_1.fq gstm1_null/reads_2.fq > /dev/null
  gzip -f gstm1_null/reads_1.fq gstm1_null/reads_2.fq
}
done_or_run gstm1_null/reads.done simulate

recruit() {   # recruit <name> <plx> <reads_1> <reads_2>: certain reads in <name>.reads/
  rm -rf "$1.reads"
  "$P" recruit --locus "$2" -1 "$3" -2 "$4" --error-rate 0.002 -o "$1.tsv" --reads-out "$1.reads" --certain-only -t 2 -q
}
done_or_run gstm1_null.done recruit gstm1_null gstm1_null/gstm.plx gstm1_null/reads_1.fq.gz gstm1_null/reads_2.fq.gz
for entry in "c4 HG00733 paralogs/c4.bed" "cyp2d HG00733 paralogs/cyp2d.bed" "gstm HG002 paralogs/gstm.bed" "acot HG002 paralogs/acot.bed" "lpa HG002 -"; do
  read -r L _ _ <<< "$entry"
  done_or_run "$L.done" recruit "$L" "$SIM/$L/g38/$L.plx" "$SIM/$L/short_locus_1.fq.gz" "$SIM/$L/short_locus_2.fq.gz"
done

{
  echo -e "case\tsample\tunit\tkmers\tmedian_depth\tcopy_number\tcalled\ttruth\tmatch"
  "${CN[@]}" --gfa gstm1_null/gstm.gfa --sample "$NULL" --genes "$RD/paralogs/gstm.bed" --reads gstm1_null.reads/*.fq.gz | tail -n +2 | sed 's/^/gstm1_null\t/'
  for entry in "c4 HG00733 paralogs/c4.bed" "cyp2d HG00733 paralogs/cyp2d.bed" "gstm HG002 paralogs/gstm.bed" "acot HG002 paralogs/acot.bed"; do
    read -r L S B <<< "$entry"
    "${CN[@]}" --gfa "$SIM/$L/$L.gfa" --sample "$S" --genes "$RD/$B" --reads "$L".reads/*.fq.gz | tail -n +2 | sed "s/^/$L\t/"
  done
  "${CN[@]}" --gfa "$SIM/lpa/lpa.gfa" --sample HG002 --repeats "$RD/tandem_repeat/lpa.repeats.tsv" --reads lpa.reads/*.fq.gz | tail -n +2 | sed 's/^/lpa\t/'
} > cn.tmp 2> cn.log
mv cn.tmp cn.tsv
column -t cn.tsv
# fates of the CYP2D sample's fragments, by reason (the hybrid shows as conflict, if anywhere)
echo; echo "cyp2d fragments by fate and reason:"
awk -F'\t' 'NR > 1 {n[$3 " " $4]++} END {for (k in n) print n[k], k}' cyp2d.tsv | sort -k1nr
