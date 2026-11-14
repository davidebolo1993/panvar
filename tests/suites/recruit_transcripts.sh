#!/usr/bin/env bash
# panvar recruit with a locus's transcripts in its index: a gene whose first three exons a paralog elsewhere
# shares (their k-mers shared), so across them only the transcript's junctions are private. The junction
# k-mers join the index flagged transcript-only; they turn shared when a genome holds a processed pseudogene,
# when the background's transcriptome holds the paralog's transcript (not for the gene's own record there),
# and when the paralog's locus is indexed in the same run. Genomic k-mers keep their labels throughout, and
# an index without transcripts its bytes.
#   recruit_transcripts.sh <panvar_bin> <tests/tools dir> <out_dir>
set -uo pipefail

if [[ $# -lt 3 ]]; then
  echo "Usage: $0 <panvar_bin> <tools_dir> <out_dir>" >&2
  exit 2
fi
BIN="$1"; TOOLS="$2"; OUT="$3"
rm -rf "$OUT"; mkdir -p "$OUT/tmp"; cd "$OUT" || exit 2
python3 "$TOOLS/junction_fixture.py" data
read -r J J12 JS < data/junction_kmers.txt   # the transcript's junction k-mers: all, those the paralog's has too,
                                             # and those of the record with a changed splice site
JO=$((J - J12))

pass=0; fail=0
ok() { echo "  ok   $1"; pass=$((pass + 1)); }
bad() { echo "  FAIL $1"; fail=$((fail + 1)); }
check() { if eval "$2"; then ok "$1"; else bad "$1"; fi; }
BG=("$BIN" background -t 2 --partitions 64 --memory 64M -q --tmp-dir "$OUT/tmp")
RC=("$BIN" recruit -t 2 --background s1 -r ref --index-only --no-calibration)
TX=(--locus data/locus.gfa --transcripts data/locus.hst.fa)
# transcript-only k-mers of a dump by class; fragments certain in a table; a dump's genomic k-mers
tx() { awk -F'\t' 'NR>1 && $2!="D" && $8==1 {n[$2]++} END {printf "P%d S%d R%d", n["P"], n["S"], n["R"]}' "$1"; }
certain() { awk -F'\t' 'NR>1 && $3=="certain"' "$1" | wc -l | tr -d ' '; }
genomic() { awk -F'\t' 'NR>1 && $8==0' "$1" | cut -f1-6 | sort; }
dump() { "$BIN" recruit --dump-index "$1/${2:-locus}.plx=$1.${2:-locus}.tsv"; }
score() { "$BIN" recruit --locus "$1/locus.plx" -i "data/$2" -o "$3" --min-sites 2 -q; }

"$BIN" transcripts -i data/locus.gfa -r ref --gtf data/genes.gtf -o data/locus -q
"$BIN" transcripts -i data/locus_b.gfa -r ref --gtf data/genes_b.gtf -o data/locus_b -q
"${BG[@]}" -o s1 --genome ref=data/ref.fa --genome ga=data/ga.fa --coordinate-reference ref
"${BG[@]}" -o s2 --genome ref=data/ref.fa --genome ga=data/ga.fa --genome gp=data/gp.fa --coordinate-reference ref
"${BG[@]}" -o s3 --genome ref=data/ref.fa --genome ga=data/ga.fa --coordinate-reference ref

echo "== junction k-mers against the genomes"
"${RC[@]}" --save-index idx0 --locus data/locus.gfa 2> log0.txt
"${RC[@]}" --save-index idx1 "${TX[@]}" 2> log1.txt; rc=$?
"$BIN" recruit -t 2 --background s2 -r ref --index-only --no-calibration --save-index idx2 --locus data/locus.gfa \
  --transcripts locus=data/locus.hst.fa 2> log2.txt
for i in 0 1 2; do dump idx$i; done
check "the index builds with the transcripts" '[[ $rc -eq 0 && -f idx1/locus.plx ]] && grep -q "k-mers only in transcripts: $J" log1.txt'
check "the four junctions add their $J transcript-only k-mers, all private ($(tx idx1.locus.tsv))" '[[ "$(tx idx1.locus.tsv)" == "P$J S0 R0" ]]'
check "the shared exons' k-mers are shared with the paralog" '[[ $(awk -F"\t" "NR>1 && \$8==0 && \$2==\"S\"" idx1.locus.tsv | wc -l) -gt 500 ]]'
check "each junction k-mer is on a node" '[[ $(awk -F"\t" "NR>1 && \$8==1 && \$7==\".\"" idx1.locus.tsv | wc -l) -eq 0 ]]'
check "a processed pseudogene elsewhere makes them shared ($(tx idx2.locus.tsv))" '[[ "$(tx idx2.locus.tsv)" == "P0 S$J R0" ]] && [[ $(awk -F"\t" "NR>1 && \$8==1 && \$3!=1" idx2.locus.tsv | wc -l) -eq 0 ]]'
check "the report counts them" 'grep -q "^locus	transcript_only	$J (private $J, shared 0, repeat 0)" idx1/locus.index.tsv && grep -q "^locus	transcriptome	none in the background" idx1/locus.index.tsv'
check "an index without transcripts has no transcript-only k-mers" '[[ "$(tx idx0.locus.tsv)" == "P0 S0 R0" ]] && ! grep -q "^locus	transcript" idx0/locus.index.tsv'
check "the genomic k-mers keep their labels" 'cmp -s <(genomic idx0.locus.tsv) <(genomic idx1.locus.tsv)'

for i in 0 1 2; do score idx$i reads.fq r$i.tsv; done
check "reads across the first two junctions are certain only with the transcripts ($(certain r0.tsv), $(certain r1.tsv), $(certain r2.tsv) of 200)" \
  '[[ $(certain r0.tsv) -eq 0 && $(certain r1.tsv) -eq 200 && $(certain r2.tsv) -eq 0 ]]'
score idx2 reads_other.fq other2.tsv; score idx1 reads_other.fq other1.tsv
check "reads of the processed pseudogene into its flank and tail are not certain" '[[ $(grep -c "^pp.*	certain	" other2.tsv) -eq 0 ]]'
check "reads across an exon skip no transcript has are not certain" '[[ $(grep -c "^skip.*	certain	" other1.tsv) -eq 0 ]]'

echo "== reuse"
cp -r idx1 idxr
"${RC[@]}" --save-index idxr "${TX[@]}" 2> reuse.txt
check "an index with the same transcripts is reused" 'grep -q "current, reused" reuse.txt'
"${RC[@]}" --save-index idxr --locus data/locus.gfa 2> drop.txt
check "the same graph without the transcripts is rebuilt, and the log says what is lost" 'grep -q "built with transcripts, none given now: rebuilding" drop.txt'
"${RC[@]}" --save-index idxr "${TX[@]}" --transcripts-ok-only 2> okonly.txt
check "given transcripts again, it is rebuilt with them, and the log says so" 'grep -q "built without transcripts, given now: rebuilding" okonly.txt'
check "taking status ok only is another index" 'grep -q "status ok only" idxr/locus.index.tsv'
"${RC[@]}" --save-index idxr "${TX[@]}" 2> other.txt
check "other settings are a rebuild without that note" 'grep -q "built from another graph, background or settings: rebuilding" other.txt'
"$BIN" recruit -t 2 --background s1 -r ref --index-only --save-index idx4 "${TX[@]}" 2> cal.txt; rc=$?
check "a calibrated index with transcripts builds" '[[ $rc -eq 0 ]] && grep -q "calibration: short reads" cal.txt'
score idx4 reads.fq r4.tsv; dump idx4
"$BIN" recruit -t 2 --background s1 -r ref --index-only --save-index idx4 "${TX[@]}" 2> cal2.txt
check "a calibrated index with transcripts loads again: reused, the same reads certain, the same k-mers" \
  'grep -q "current, reused" cal2.txt && [[ $(certain r4.tsv) -eq 200 ]] && cmp -s <(cut -f1-5,7,8 idx1.locus.tsv | sort) <(cut -f1-5,7,8 idx4.locus.tsv | sort)'
cat data/locus.hst.fa data/splice_changed.fa > data/sc.hst.fa   # one more record, its splice site changed
"${RC[@]}" --save-index idxc --locus data/locus.gfa --transcripts data/sc.hst.fa 2> sc.txt; dump idxc
check "a transcript with a changed splice site is left out, and the log counts it by status ($(tx idxc.locus.tsv))" \
  '[[ "$(tx idxc.locus.tsv)" == "P$J S0 R0" ]] && grep -q "by status: ok 1/1, splice_site_changed 0/1" sc.txt && grep -q "^locus	transcripts_by_status	ok 1/1, splice_site_changed 0/1 (taken/all)" idxc/locus.index.tsv'
"${RC[@]}" --save-index idxc --locus data/locus.gfa --transcripts data/sc.hst.fa --transcripts-splice-changed 2> sc2.txt; dump idxc
check "--transcripts-splice-changed takes it: its junction joins the index ($(tx idxc.locus.tsv))" \
  '[[ "$(tx idxc.locus.tsv)" == "P$((J + JS)) S0 R0" ]] && grep -q "rebuilding" sc2.txt && grep -q "splice_site_changed 1/1" sc2.txt && grep -q "with changed splice sites" idxc/locus.index.tsv'
"${RC[@]}" --save-index idxc --locus data/locus.gfa --transcripts data/sc.hst.fa --transcripts-splice-changed --transcripts-ok-only 2> sc3.txt
check "with --transcripts-ok-only it is refused" 'grep -q "not both" sc3.txt'
cp data/locus.hst.fa data/short.hst.fa   # a table whose first row was cut short
{ sed -n 1p data/locus.hst.tsv; sed -n 2p data/locus.hst.tsv | cut -f1-10; sed -n '3,$p' data/locus.hst.tsv; } > data/short.hst.tsv
"${RC[@]}" --save-index idxt --locus data/locus.gfa --transcripts data/short.hst.fa 2> short.txt; rc=$?
check "a table row with fewer fields than the header is skipped with a warning" \
  '[[ $rc -eq 0 ]] && grep -q "warning: locus: locus; .*short.hst.tsv: 1 row(s) with fewer fields than the header" short.txt'
"${RC[@]}" --save-index idx3 --locus data/locus.gfa --transcripts other=data/locus.hst.fa 2> bad.txt
check "transcripts for a locus not given are refused" 'grep -q "no --locus graph is named other" bad.txt'

echo "== the transcriptome"
"${BG[@]}" -o s1 --add --transcriptome tx=data/tx_own.fa; rc=$?
check "a transcriptome is added to the store" '[[ $rc -eq 0 ]] && grep -q "\"transcriptomes\"" s1/manifest.json && [[ -d s1/tx-0 ]]'
check "its manifest is format version 3, which a panvar without transcriptomes refuses; a store without one stays version 2" \
  'grep -q "\"version\": 3," s1/manifest.json && grep -q "\"version\": 2," s3/manifest.json'
"${RC[@]}" --save-index idx0b --locus data/locus.gfa 2> log0b.txt
check "an index without transcripts is the same file with a transcriptome in the store" 'cmp -s idx0/locus.plx idx0b/locus.plx'
cp -r idx1 idxs
"${RC[@]}" --save-index idxs "${TX[@]}" 2> stale.txt
check "an index with transcripts built before the transcriptome is rebuilt" 'grep -q "rebuilding" stale.txt'
"${RC[@]}" --save-index idx5 "${TX[@]}" 2> log5.txt; dump idx5
check "the gene's own record there is subtracted: its junctions stay private ($(tx idx5.locus.tsv))" \
  '[[ "$(tx idx5.locus.tsv)" == "P$J S0 R0" ]] && grep -q "records of the locus.s genes subtracted: 1;" log5.txt'
"${BG[@]}" -o s1 --add --transcriptome tx=data/tx_paralog.fa 2> replace.txt; rc=$?
check "a transcriptome of the same name replaces it" '[[ $rc -eq 0 && $(grep -c "\"dir\": \"tx-" s1/manifest.json) -eq 1 ]]'
"${RC[@]}" --save-index idx6 "${TX[@]}" 2> log6.txt; dump idx6
check "the paralog's transcript makes the junctions it shares shared, the others stay private ($(tx idx6.locus.tsv))" \
  '[[ "$(tx idx6.locus.tsv)" == "P$JO S$J12 R0" ]] && [[ $(awk -F"\t" "NR>1 && \$8==1 && \$2==\"S\" && \$3!=1" idx6.locus.tsv | wc -l) -eq 0 ]]'
check "the report names the transcriptome and counts them" \
  'grep -q "^locus	transcriptome	tx (records of the locus.s genes subtracted: 1)" idx6/locus.index.tsv && grep -q "^locus	transcript_only_elsewhere	in other transcripts: $J12; in other loci of the run: 0" idx6/locus.index.tsv'
check "the genomic k-mers keep their labels with a transcriptome" 'cmp -s <(genomic idx1.locus.tsv) <(genomic idx6.locus.tsv)'
score idx5 reads.fq r5.tsv; score idx6 reads.fq r6.tsv; score idx6 reads_own.fq o6.tsv; score idx6 reads_other.fq other6.tsv
check "reads of the paralog's mRNA are not certain once the transcriptome holds its transcript ($(certain r5.tsv), $(certain r6.tsv) of 200)" \
  '[[ $(certain r5.tsv) -eq 200 && $(certain r6.tsv) -eq 0 ]]'
check "the gene's own reads, across its own junctions, stay certain ($(certain o6.tsv) of 200)" '[[ $(certain o6.tsv) -eq 200 ]]'
check "reads across an exon skip no transcript has are not certain with a transcriptome" '[[ $(grep -c "^skip.*	certain	" other6.tsv) -eq 0 ]]'
"${BG[@]}" -o s1 --clean
check "cleaning keeps the listed transcriptome and removes the replaced one" '[[ $(ls -d s1/tx-* | wc -l) -eq 1 ]]'
"${BG[@]}" -o s1 --add --genome gn=data/gn.fa && "${BG[@]}" -o s1 --merge; rc=$?
check "merging the genomes' segments leaves the transcriptome listed" \
  '[[ $rc -eq 0 && $(grep -c "\"dir\": \"tx-" s1/manifest.json) -eq 1 && $(grep -c "\"dir\": \"seg-" s1/manifest.json) -eq 1 ]]'
"${BG[@]}" -o s1 --compact; rc=$?
"${RC[@]}" --save-index idx7 "${TX[@]}" 2> log7.txt; dump idx7
check "a merged and compacted background reads the transcriptome beside it: the same labels" \
  '[[ $rc -eq 0 ]] && grep -q "background: compacted" log7.txt && cmp -s idx6.locus.tsv idx7.locus.tsv'
"${BG[@]}" -o s1 --add --genome x=data/ref.fa --transcriptome tx2=data/tx_own.fa 2> mixed.txt
check "a transcriptome is not added in a run with genomes" 'grep -q "transcriptome goes with --add alone" mixed.txt'

echo "== the paralog's locus in the same run"
RC3=("$BIN" recruit -t 2 --background s3 -r ref --index-only --no-calibration)
"${RC3[@]}" --save-index idx8 --locus data/locus.gfa --locus data/locus_b.gfa --transcripts locus=data/locus.hst.fa \
  --transcripts locus_b=data/locus_b.hst.fa 2> log8.txt; rc=$?
dump idx8; dump idx8 locus_b
check "junction k-mers the other locus's transcripts hold are shared, the others private ($(tx idx8.locus.tsv))" \
  '[[ $rc -eq 0 && "$(tx idx8.locus.tsv)" == "P$JO S$J12 R0" ]] && grep -q "in other loci of the run (locus_b): $J12" idx8/locus.index.tsv'
check "the genomic k-mers keep their labels beside another locus" 'cmp -s <(genomic idx1.locus.tsv) <(genomic idx8.locus.tsv)'
"${RC3[@]}" --save-index idx9 --locus data/locus_b.gfa --transcripts data/locus_b.hst.fa 2> log9.txt; dump idx9 locus_b
check "indexed alone, the paralog's junction k-mers are all private ($(tx idx9.locus_b.tsv))" \
  '[[ "$(tx idx9.locus_b.tsv)" == "P$(wc -l < data/junction_nodes_b.tsv | tr -d " ") S0 R0" ]]'
check "over several nodes, a reversed walk and a minus-strand transcript, each junction k-mer is on the node of its first base" \
  'cmp -s <(awk -F"\t" "NR>1 && \$8==1 {print \$1\"\t\"\$7}" idx9.locus_b.tsv | sort) data/junction_nodes_b.tsv'

echo "recruit_transcripts: $pass passed, $fail failed"
[[ $fail -eq 0 ]]
