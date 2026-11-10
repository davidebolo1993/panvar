#!/usr/bin/env bash
# Compare the copy numbers of scripts/call_examples.sh with the truth sets in tests/real_data (module totals and,
# where call ran with a GTF, per-gene splits) and plot the comparison, in <call output dir>/copy_number.
#   dev/call_truth.sh [call_examples.sh output dir]
set -euo pipefail
here=$PWD
cd "$(dirname "$0")/.."
arg() { if [ -z "$1" ]; then echo "$2"; elif [ "${1#/}" != "$1" ]; then echo "$1"; else echo "$here/$1"; fi; }
for tool in python3 Rscript; do command -v "$tool" > /dev/null || { echo "$tool not found on PATH"; exit 1; }; done

OUT=$(arg "${1:-}" tests/examples/call)
DATA=tests/real_data

mkdir -p "$OUT/copy_number"
tables=()
compare() {  # <level> <label> <compare_copy_number.py options>
  # Tables are named by level too: a gene and a locus can share a name, and some file systems ignore case.
  local stem=$OUT/copy_number/$1.$2 label=$2
  shift 2
  python3 dev/compare_copy_number.py --label "$label" -o "$stem.tsv" "$@" > "$stem.txt"
  tables+=(--table "$stem.tsv")
}
locus() {  # <name> <truth> <genes in the module's copies>
  local vcf=$OUT/$1/$1.region.vcf
  [ -f "$vcf" ] || { echo "no $vcf: run scripts/call_examples.sh first"; exit 1; }
  compare locus "$1" --vcf "$vcf" --truth "$2" --genes "$3"
}
per_gene() {  # <name> <truth BED>: every gene in both call's per-gene table and the truth, per site
  local name=$1 truth=$2 table=$OUT/$1/$1.dup_gene_cn.tsv g site sites
  [ -f "$table" ] || return 0
  for g in $(comm -12 <(awk -F'\t' 'NR > 1 { print $4 }' "$table" | sort -u) \
                      <(awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) if ($i == "gene") c = i; next } { print $c }' "$truth" | sort -u)); do
    sites=$(awk -F'\t' -v g="$g" 'NR > 1 && $4 == g { print $1 }' "$table" | sort -u)
    for site in $sites; do
      local label=$g; [ "$(echo "$sites" | wc -l)" -gt 1 ] && label=$g@$site
      compare gene "$label" --gene-table "$table" --truth "$truth" --genes "$g" --site "$site"
    done
  done
}

# Truth: KIV-2 copies for the tandem repeat; for the paralogs, pangene's count of the genes in the copies.
P=$DATA/paralogs
locus lpa "$DATA/tandem_repeat/lpa.repeats.tsv" ""
locus lpa_lossy "$DATA/tandem_repeat/lpa.repeats.tsv" ""
locus c4 "$P/c4.bed" C4A,C4B
locus gstm "$P/gstm.bed" GSTM1,GSTM2
locus cyp2d "$P/cyp2d.bed" CYP2D6,CYP2D7
locus acot "$P/acot.bed" ACOT1,ACOT2
per_gene c4 "$P/c4.bed"
per_gene gstm "$P/gstm.bed"
per_gene cyp2d "$P/cyp2d.bed"
per_gene acot "$P/acot.bed"

Rscript dev/plot_cn_correlation.R "${tables[@]}" --out "$OUT/copy_number/cn_correlation"
column -t "$OUT/copy_number/cn_correlation.summary.tsv"
