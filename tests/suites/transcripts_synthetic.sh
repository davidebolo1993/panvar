#!/usr/bin/env bash
# panvar transcripts on synthetic loci (transcript_fixture.py cases): exon SNPs and deletions, splice-site
# variants, a deletion across an exon boundary, duplicated and deleted genes, minus-strand genes and walks
# written reversed, the other cases its CASES list, and 100 seeded random cases. Every row's status, notes
# and sequence against transcripts spliced from each haplotype independently; the same output at 1 and 2
# threads; input errors stop the run.
#   transcripts_synthetic.sh <panvar_bin> <tests/tools dir> <out_dir>
set -uo pipefail

if [[ $# -lt 3 ]]; then
  echo "Usage: $0 <panvar_bin> <tools_dir> <out_dir>" >&2
  exit 2
fi
BIN="$1"; TOOLS="$2"; OUT="$3"
rm -rf "$OUT"; mkdir -p "$OUT"; cd "$OUT" || exit 2
python3 "$TOOLS/transcript_fixture.py" write main .

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
check "the reference walk gives the reference transcripts" '[[ "$(awk -F"\t" "\$7 ~ /^ref#/ {print \$20}" t1.hst.tsv | sort -u)" == yes ]]'
check "a duplicated gene gives one transcript per copy" '[[ "$(col TC1 h5 copy | tr "\n" ,)" == "1,2," && "$(col TC1 h5 copies | sort -u)" == 2 ]]'
check "copies lie apart on the walk" '[[ "$(col TC1 h5 path_start | tr "\n" ,)" == "8000,8900," ]]'
check "a deleted gene has no copy" '[[ "$(col TC1 h6 status)" == no_copy && "$(col TC1 h6 copy)" == 0 ]]'
check "an exon deletion drops only the transcripts using it" '[[ "$(col TA1 h2 status)" == missing_exon && "$(col TA1 h2 detail)" == e2:missing && "$(col TA2 h2 status)" == ok ]]'
check "a boundary inside a deletion moves to the nearer shared base" '[[ "$(col TA1 h4 status)" == boundary_moved && "$(col TA1 h4 moved_boundaries)" == 1 ]]'
check "a reversed walk is read in the reference orientation" '[[ "$(col TA1 h7 orientation)" == - && "$(col TA1 h7 path_start)" == 8700 && "$(col TA1 h7 path_end)" == 11000 ]]'
check "the minus-strand transcript lists its exons 5 to 3" '[[ "$(col TB1 ref exon_starts)" == "5600,5000" && "$(col TB1 ref exon_lengths)" == "200,100" ]]'
check "FASTA names are transcript|gene|haplotype|copy" 'grep -q "^>TC1|GC|h5#1#ctg5:0-12000|2 source=projected status=ok$" t1.hst.fa'
check "the summary counts copies per haplotype" '[[ "$(awk -F"\t" "\$2==\"GC\"{print \$6}" t1.hst.summary.tsv)" == "0:1,1:6,2:2" ]]'

# a transcript whose exons overlap on the reference: the reference walk does not give it back, and its rows,
# all missing an exon, keep that status and are noted
cp genes.gtf naa.gtf
a='gene_id "GA"; gene_name "GA"; gene_biotype "protein_coding"; transcript_id "TA9"; transcript_name "TA9";'
printf 'chr1\ttest\texon\t2001\t2301\t.\t+\t.\t%s\nchr1\ttest\texon\t2251\t2401\t.\t+\t.\t%s\n' "$a" "$a" >> naa.gtf
"$BIN" transcripts -i locus.gfa -r ref --gtf naa.gtf -o naa -t 2 -q 2> naa.log
check "a transcript the reference walk does not give back is noted on every row that keeps a more severe status" \
  'grep -q "does not give back 1 transcript(s) as annotated: TA9" naa.log && [[ "$(awk -F"\t" "\$1==\"TA9\" {print \$12, (\$25 ~ /not_as_annotated/)}" naa.hst.tsv | sort -u)" == "missing_exon 1" ]]'

# the checker refuses a table or FASTA naming one transcript, haplotype and copy twice
{ cat t1.hst.tsv; sed -n 2p t1.hst.tsv; } > dup_row.hst.tsv; cp t1.hst.fa dup_row.hst.fa
cp t1.hst.tsv dup_fa.hst.tsv; { cat t1.hst.fa; head -2 t1.hst.fa; } > dup_fa.hst.fa
check "the checker refuses a repeated row" '! python3 "$TOOLS/transcript_fixture.py" check . dup_row > /dev/null 2>&1'
check "the checker refuses a repeated FASTA name" '! python3 "$TOOLS/transcript_fixture.py" check . dup_fa > /dev/null 2>&1'

# every case, each row against the haplotypes' own transcripts
for c in $(python3 "$TOOLS/transcript_fixture.py" cases); do
  [[ $c == main ]] && continue
  mkdir -p "case_$c"
  python3 "$TOOLS/transcript_fixture.py" write "$c" "case_$c"
  extra=()
  [[ -e "case_$c/hap_annotations.tsv" ]] && extra=(--hap-annotations "case_$c/hap_annotations.tsv" --contig-alias "case_$c/alias.tsv")
  "$BIN" transcripts -i "case_$c/locus.gfa" -r ref --gtf "case_$c/genes.gtf" -o "case_$c/out" -t 2 -q ${extra[@]+"${extra[@]}"} \
    2> "case_$c/log"
  check "case $c: the reference walk gives back its annotation" '! grep -q "does not give back" "case_$c/log"'
  check "case $c: every row matches the haplotypes' own transcripts" 'python3 "$TOOLS/transcript_fixture.py" check "case_$c" "case_$c/out" > "case_$c/check.log"; cat "case_$c/check.log"; grep -q "wrong: 0" "case_$c/check.log"'
done

# seeded random cases (random_case): genes, SNPs, deletions, substitution blocks and whole-gene tandem
# duplications, every row against the haplotypes' own transcripts
wrong=""
for seed in $(seq 1 100); do
  d="random_$seed"; mkdir -p "$d"
  python3 "$TOOLS/transcript_fixture.py" write "random:$seed" "$d"
  "$BIN" transcripts -i "$d/locus.gfa" -r ref --gtf "$d/genes.gtf" -o "$d/out" -t 2 -q 2> "$d/log"
  if python3 "$TOOLS/transcript_fixture.py" check "$d" "$d/out" > "$d/check.log" && ! grep -q "does not give back" "$d/log"; then
    rm -rf "$d"
  else
    wrong+=" $seed"
  fi
done
check "100 random cases match the haplotypes' own transcripts${wrong:+ (seeds failing:$wrong)}" '[[ -z "$wrong" ]]'

# two transcripts over one tandem stretch whose last units end apart: each places its own units
mkdir -p pair && python3 "$TOOLS/transcript_fixture.py" write tandem_pair pair
"$BIN" transcripts -i pair/locus.gfa -r ref --gtf pair/genes.gtf -o pair/both -t 1 -q 2> pair/both.log
"$BIN" transcripts -i pair/locus.gfa -r ref --gtf pair/alone.gtf -o pair/alone -t 1 -q 2> pair/alone.log
check "a tandem transcript is placed the same beside another one with a shorter last unit" \
  'cmp -s <(grep "^TT2" pair/both.hst.tsv) <(grep "^TT2" pair/alone.hst.tsv) && ! grep -q "does not give back" pair/both.log'

# a gene on allele nodes, and one cut below the share
dv() { awk -F'\t' -v t="$1" -v h="$2" '$1==t && $7 ~ "^"h"#" {print $12":"$8}' "${3:-case_divergent/out}.hst.tsv"; }
check "a gene whose exons lie on allele nodes is placed by sequence" \
  '[[ "$(dv TH1 h1)" == ok:1 && "$(dv TH1 h4)" == ok:1 ]] && [[ $(grep -c "placed_by_sequence;best_share:0.750;identity:0.950" case_divergent/out.hst.tsv) -eq 4 ]]'
check "a gene cut below the share is a partial copy, with no sequence and no copy counted" \
  '[[ "$(dv TP1 h2)" == partial_copy:0 ]] && ! grep -q "^>TP1|GP|h2#" case_divergent/out.hst.fa && [[ "$(awk -F"\t" "\$2==\"GP\"{print \$6}" case_divergent/out.hst.summary.tsv)" == "0:2,1:4" ]]'
check "a partial copy on a walk written reversed has orientation -" \
  '[[ "$(awk -F"\t" "\$1==\"TP1\" && \$7 ~ /^h[25]#/ {print \$13}" case_divergent/out.hst.tsv | tr "\n" ,)" == "+,-," ]]'
check "an exon the stretch's alignment misses is found alone" '[[ $(grep -c "e4:aligned_alone" case_intron_shrunk/out.hst.tsv) -eq 2 ]]'
"$BIN" transcripts -i case_divergent/locus.gfa -r ref --gtf case_divergent/genes.gtf -o strict --min-identity 0.99 -q
check "below --min-identity the gene has no copy" '[[ "$(dv TH1 h1 strict)" == no_copy:0 ]]'

# haplotypes' own annotations: copies they disagree on, and annotated transcripts alone
printf 'molecule\tgene\tstart\tend\tstrand\nh1#1#ctg1:0-12000\tGB\t0\t1\t0\nh1#1#ctg1:0-12000\tGB\t0\t1\t0\n' > gb.bed
check "a copy the projection lacks and the annotation has off other genes' copies counts" \
  'python3 "$TOOLS/../../dev/transcripts_copies.py" case_annotated/out.hst.tsv gb.bed | grep -q "^GB	1	1	2	2	"'
grep -e haplotype -e '^h3#' case_annotated/hap_annotations.tsv > h3_only.tsv
"$BIN" transcripts -i case_annotated/locus.gfa -r ref --gtf case_annotated/genes.gtf --hap-annotations h3_only.tsv \
  --contig-alias case_annotated/alias.tsv -o h3_only -q
check "copies are compared only for genes an annotation has" \
  '[[ "$(awk -F"\t" "NR > 1 && \$2 ~ /^G[ABC]\$/ {print \$NF}" h3_only.hst.summary.tsv | sort -u)" == . ]]'
check "an annotated gene under another name counts by its source gene" \
  '[[ "$(awk -F"\t" "\$2==\"GC\" {print \$NF}" case_annotated/out.hst.summary.tsv)" != *h3#1* ]]'
check "gene records sharing a source gene are copies of one gene" \
  '[[ "$(awk -F"\t" "\$2==\"GC\" {print \$NF}" case_annotated/out.hst.summary.tsv)" != *h5#1* ]]'
check "an exon without its transcript record is reported" 'grep -q "1 exon row(s) whose transcript record is missing" case_annotated/log'
check "the summary notes a gene's copies where the annotation differs" 'awk -F"\t" "\$2==\"GC\"" case_annotation_rules/out.hst.summary.tsv | grep -q "h6#1:1/2"'
"$BIN" transcripts -i case_annotated/locus.gfa --hap-annotations case_annotated/hap_annotations.tsv \
  --contig-alias case_annotated/alias.tsv -o ann_only -q
check "without --gtf, only the annotated transcripts" '[[ "$(tail -n +2 ann_only.hst.tsv | cut -f10 | sort | uniq -c | tr -s " ")" == " 13 annotated" ]]'

# the rules between the projection and the annotation (case annotation_rules): rows checked above; the
# annotation preferred; the summary's counts
R=case_annotation_rules
"$BIN" transcripts -i $R/locus.gfa -r ref --gtf $R/genes.gtf --hap-annotations $R/hap_annotations.tsv --contig-alias $R/alias.tsv \
  --prefer-annotation -o $R/pref -t 2 -q
check "with --prefer-annotation, annotated records win and none is kept out" \
  'python3 "$TOOLS/transcript_fixture.py" check $R $R/pref expected_prefer.tsv > $R/pref.log; cat $R/pref.log; grep -q "wrong: 0" $R/pref.log'
sm() { awk -F'\t' -v g="$2" -v c="$3" 'NR==1{for(i=1;i<=NF;i++) ix[$i]=i; next} $1==g {print $ix[c]}' "$1.hst.summary.tsv"; }   # by gene id
check "the summary counts the records each rule keeps out" \
  '[[ "$(sm $R/out GA discordant) $(sm $R/out GC paralog_lift) $(sm $R/out o4 paralog_lift) $(sm $R/out ENSGX9 paralog_lift)" == "3 2 1 0" ]]'
check "a paralog lift is no copy, a filled copy is one" '[[ "$(sm $R/out GC haplotypes_by_copies)" == "1:7,2:2" && "$(sm $R/pref GC haplotypes_by_copies)" == "1:6,2:3" ]]'
check "without annotations the summary has no columns for their statuses" '! head -1 t1.hst.summary.tsv | grep -q -e discordant -e paralog_lift'
"$BIN" transcripts -i locus.gfa -r ref --gtf genes.gtf --prefer-annotation -o e0 -q 2> e0.log
check "--prefer-annotation without --hap-annotations is refused" 'grep -q "needs --hap-annotations" e0.log'

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
