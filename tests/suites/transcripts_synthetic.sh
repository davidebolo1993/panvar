#!/usr/bin/env bash
# panvar transcripts on a synthetic locus: an exon SNP, an exon deletion, a splice-site SNP, a deletion across
# an exon boundary, a duplicated and a deleted gene, a minus-strand gene and a walk written reversed. Every
# row's status, notes and sequence against transcripts spliced from each haplotype independently; the same
# output at 1 and 2 threads; input errors stop the run.
#   transcripts_synthetic.sh <panvar_bin> <tests/tools dir> <out_dir>
set -uo pipefail

if [[ $# -lt 3 ]]; then
  echo "Usage: $0 <panvar_bin> <tools_dir> <out_dir>" >&2
  exit 2
fi
BIN="$1"; TOOLS="$2"; OUT="$3"
rm -rf "$OUT"; mkdir -p "$OUT"; cd "$OUT" || exit 2
python3 "$TOOLS/transcript_fixture.py" write .

pass=0; fail=0
ok() { echo "  ok   $1"; pass=$((pass + 1)); }
bad() { echo "  FAIL $1"; fail=$((fail + 1)); }
check() { if eval "$2"; then ok "$1"; else bad "$1"; fi; }
col() { awk -F'\t' -v t="$1" -v h="$2" -v c="$3" 'NR==1{for(i=1;i<=NF;i++) ix[$i]=i; next}
         $1==t && $7 ~ "^"h"#" {print $ix[c]}' t1.hst.tsv; }

"$BIN" transcripts -i locus.gfa -r ref --gtf genes.gtf -o t1 -t 1 -q || { echo "FAIL transcripts did not run"; exit 1; }
"$BIN" transcripts -i locus.gfa -r ref --gtf genes.gtf -o t2 -t 2 -q
check "every row matches the haplotypes' own transcripts" 'python3 "$TOOLS/transcript_fixture.py" check . t1 > check.log; cat check.log; grep -q "wrong: 0" check.log'
check "2 threads write the same files as 1" 'cmp -s t1.hst.tsv t2.hst.tsv && cmp -s t1.hst.fa t2.hst.fa && cmp -s t1.hst.summary.tsv t2.hst.summary.tsv'
check "the reference walk gives the reference transcripts" '[[ "$(awk -F"\t" "\$7 ~ /^ref#/ {print \$19}" t1.hst.tsv | sort -u)" == yes ]]'
check "a duplicated gene gives one transcript per copy" '[[ "$(col TC1 h5 copy | tr "\n" ,)" == "1,2," && "$(col TC1 h5 copies | sort -u)" == 2 ]]'
check "copies lie apart on the walk" '[[ "$(col TC1 h5 path_start | tr "\n" ,)" == "8000,8900," ]]'
check "a deleted gene has no copy" '[[ "$(col TC1 h6 status)" == no_copy && "$(col TC1 h6 copy)" == 0 ]]'
check "an exon deletion drops only the transcripts using it" '[[ "$(col TA1 h2 status)" == missing_exon && "$(col TA1 h2 detail)" == e2:missing && "$(col TA2 h2 status)" == ok ]]'
check "a boundary inside a deletion moves to the nearer shared base" '[[ "$(col TA1 h4 status)" == boundary_moved && "$(col TA1 h4 moved_boundaries)" == 1 ]]'
check "a reversed walk is read in the reference orientation" '[[ "$(col TA1 h7 orientation)" == - && "$(col TA1 h7 path_start)" == 8700 && "$(col TA1 h7 path_end)" == 11000 ]]'
check "the minus-strand transcript lists its exons 5 to 3" '[[ "$(col TB1 ref exon_starts)" == "5600,5000" && "$(col TB1 ref exon_lengths)" == "200,100" ]]'
check "FASTA names are transcript|gene|haplotype|copy" 'grep -q "^>TC1|GC|h5#1#ctg5:0-12000|2 source=projected status=ok$" t1.hst.fa'
check "the summary counts copies per haplotype" '[[ "$(awk -F"\t" "\$2==\"GC\"{print \$5}" t1.hst.summary.tsv)" == "0:1,1:6,2:2" ]]'

# input errors
"$BIN" transcripts -i locus.gfa -r ref -o e1 -q 2> e1.log
check "a missing --gtf is refused" 'grep -q "Missing required input: --gtf" e1.log'
"$BIN" transcripts -i locus.gfa -r nosuchpath --gtf genes.gtf -o e2 -q 2> e2.log
check "an unknown reference path is refused" 'grep -q "reference path not found" e2.log && [[ ! -e e2.hst.tsv ]]'
sed 's/^P\tref#0#chr1:1000-13000/P\trefpath/' locus.gfa > noint.gfa
"$BIN" transcripts -i noint.gfa -r refpath --gtf genes.gtf -o e3 -q 2> e3.log
check "a reference path without contig:start-end is refused" 'grep -q "carries no contig:start-end" e3.log'

echo "transcripts_synthetic: $pass passed, $fail failed"
[[ $fail -eq 0 ]]
