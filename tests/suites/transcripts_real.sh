#!/usr/bin/env bash
# panvar transcripts on real hyperpolymorphic loci: a single-copy gene must get one copy on nearly every
# haplotype, by nodes or, where its exons lie on allele nodes, by sequence.
#   transcripts_real.sh <panvar_bin> <tests/real_data dir> <out_dir>
set -uo pipefail

if [[ $# -lt 3 ]]; then
  echo "Usage: $0 <panvar_bin> <real_data_dir> <out_dir>" >&2
  exit 2
fi
BIN="$1"; DATA="$2"; OUT="$3"
rm -rf "$OUT"; mkdir -p "$OUT"; cd "$OUT" || exit 2
GTF="$DATA/Homo_sapiens.GRCh38.116.loci.gtf.gz"

pass=0; fail=0
ok() { echo "  ok   $1"; pass=$((pass + 1)); }
bad() { echo "  FAIL $1"; fail=$((fail + 1)); }
check() { if eval "$2"; then ok "$1"; else bad "$1"; fi; }
# haplotypes with exactly one copy of a gene, and all haplotypes, from the summary's haplotypes_by_copies
one() { awk -F'\t' -v g="$2" '$2==g {n=split($6, a, ","); for (i=1; i<=n; i++) {split(a[i], c, ":"); all+=c[2]; if (c[1]==1) one=c[2]}} END {print one+0, all+0}' "$1.hst.summary.tsv"; }

for pair in hla_a:HLA-A hla_b:HLA-B hla_c:HLA-C; do
  locus=${pair%%:*}; gene=${pair#*:}
  "$BIN" transcripts -i "$DATA/hla/$locus.gfa.gz" -r GRCh38 --gtf "$GTF" -o "$locus" -t 2 -q || { bad "$locus: transcripts did not run"; continue; }
  read -r n all <<< "$(one "$locus" "$gene")"
  check "$gene has one copy on at least 95% of the haplotypes ($n of $all)" '[[ $all -gt 0 && $((n * 100)) -ge $((all * 95)) ]]'
done
check "copies of HLA-C are placed by sequence where its nodes fall short ($(grep -c placed_by_sequence hla_c.hst.tsv) rows)" \
  '[[ $(awk -F"\t" "\$4==\"HLA-C\" && \$25 ~ /placed_by_sequence/ {print \$7}" hla_c.hst.tsv | sort -u | wc -l) -gt 100 ]]'

echo "transcripts_real: $pass passed, $fail failed"
[[ $fail -eq 0 ]]
