#!/usr/bin/env bash
# panvar recruit's input layer on a synthetic locus: the same pairs as FASTQ pairs, interleaved FASTQ, a
# name-sorted SAM and a coordinate-sorted SAM (mates far apart, on other contigs, unmapped; secondary,
# supplementary and QC-failed records to skip) give the same table at 1 and 8 threads; --dry-run counts
# the records; a sample sheet writes one table per sample; broken inputs stop the run with the file named.
#   recruit_input.sh <panvar_bin> <tests/tools dir> <out_dir>
set -uo pipefail

if [[ $# -lt 3 ]]; then
  echo "Usage: $0 <panvar_bin> <tools_dir> <out_dir>" >&2
  exit 2
fi
BIN="$1"; TOOLS="$2"; OUT="$3"
rm -rf "$OUT"; mkdir -p "$OUT/tmp"; cd "$OUT" || exit 2
python3 "$TOOLS/locus_fixture.py" data

pass=0; fail=0
ok() { echo "  ok   $1"; pass=$((pass + 1)); }
bad() { echo "  FAIL $1"; fail=$((fail + 1)); }
check() { if eval "$2"; then ok "$1"; else bad "$1"; fi; }

"$BIN" background -t 2 --partitions 64 --memory 64M -q --tmp-dir "$OUT/tmp" -o s \
  --genome ref=data/ref.fa --genome ga=data/ga.fa --genome gb=data/gb.fa --coordinate-reference ref
"$BIN" recruit -t 2 -q --background s -r ref --index-only --save-index idx --locus data/locus.gfa
IDX=idx/locus.plx

# pairs from the locus, its paralog and elsewhere, as FASTQ and as SAM
python3 "$TOOLS/locus_reads.py" data/ref.fa chr1 120000 160000 1500 loc 1 a
python3 "$TOOLS/locus_reads.py" data/ref.fa chr1 231000 259000 1500 par 2 b
python3 "$TOOLS/locus_reads.py" data/ref.fa chr1 0 100000 1500 far 3 c
cat a_1.fq b_1.fq c_1.fq > r_1.fq
cat a_2.fq b_2.fq c_2.fq > r_2.fq
paste -d '\n' <(paste - - - - < r_1.fq) <(paste - - - - < r_2.fq) | tr '\t' '\n' > inter.fq
python3 "$TOOLS/reads_sam.py" r_1.fq r_2.fq name.sam coord.sam 11

sorted() { { head -1 "$1"; tail -n +2 "$1" | sort; } > "$1.sorted"; }
"$BIN" recruit --locus "$IDX" -1 r_1.fq -2 r_2.fq -o fq.tsv -t 1 -q; sorted fq.tsv
check "FASTQ pairs are recruited ($(awk -F'\t' 'NR>1 && $3=="certain"' fq.tsv | wc -l | tr -d ' ') certain)" '[[ $(awk -F"\t" '"'"'NR>1 && $3=="certain"'"'"' fq.tsv | wc -l) -gt 1000 ]]'
for t in 1 8; do
  "$BIN" recruit --locus "$IDX" -i inter.fq --interleaved -o inter.$t.tsv -t $t -q; sorted inter.$t.tsv
  "$BIN" recruit --locus "$IDX" -i name.sam -o name.$t.tsv -t $t -q; sorted name.$t.tsv
  "$BIN" recruit --locus "$IDX" -i coord.sam -o coord.$t.tsv -t $t 2> coord.$t.log; sorted coord.$t.tsv
  check "interleaved FASTQ gives the same table, $t thread(s)" 'cmp -s inter.$t.tsv.sorted fq.tsv.sorted'
  check "a name-sorted SAM gives the same table, $t thread(s)" 'cmp -s name.$t.tsv.sorted fq.tsv.sorted'
  check "a coordinate-sorted SAM gives the same table, $t thread(s)" 'cmp -s coord.$t.tsv.sorted fq.tsv.sorted'
done
check "far mates are fetched by a second pass ($(grep -o '[0-9]* mates fetched' coord.1.log))" 'grep -q "second pass" coord.1.log && grep -q "mates fetched" coord.1.log'
check "the coordinate-sorted SAM counts every fragment" 'grep -q "4500 fragments (4500 pairs)" coord.1.log'
check "QC-failed reads never reach the table" '! grep -q qcfail coord.1.tsv'
check "8 threads write the same rows in the same order as 1" 'cmp -s coord.1.tsv coord.8.tsv && cmp -s name.1.tsv name.8.tsv'

# mates that overlap wholly (150 bp fragments): the pair holds exactly the sites of one mate
python3 "$TOOLS/locus_reads.py" data/ref.fa chr1 120000 160000 1000 ovl 9 ovl 150
"$BIN" recruit --locus "$IDX" -1 ovl_1.fq -2 ovl_2.fq -o ovl_pair.tsv -q
"$BIN" recruit --locus "$IDX" -i ovl_1.fq -o ovl_one.tsv -q
check "wholly overlapping mates count each site once" 'diff <(cut -f1,8 ovl_pair.tsv | sort) <(cut -f1,8 ovl_one.tsv | sort) > /dev/null'

# mates lost to the QC filter: their partners are scored alone and reported
awk -F'\t' 'BEGIN {OFS = "\t"} !/^@/ && int($2 / 128) % 2 == 1 && int($2 / 512) % 2 == 0 && ++n % 25 == 0 {$2 += 512} {print}' coord.sam > orphans.sam
"$BIN" recruit --locus "$IDX" -i orphans.sam -o orphans.tsv -t 2 2> orphans.log; rc=$?
check "reads whose mate was filtered are scored alone ($(grep -o '[0-9]* reads without their mate' orphans.log))" '[[ $rc -eq 0 ]] && grep -q "reads without their mate" orphans.log'

# an unsorted file with every read 1 before every read 2: past the hold limit, only reads hitting a screen wait
{ grep '^@' name.sam | sed 's/SO:queryname/SO:unsorted/'; grep -v '^@' name.sam | awk -F'\t' 'int($2 / 64) % 2 == 1'; grep -v '^@' name.sam | awk -F'\t' 'int($2 / 64) % 2 == 0'; } > split.sam
PANVAR_TEST_HOLD_MAX=100 "$BIN" recruit --locus "$IDX" -i split.sam -o split.tsv -t 2 2> split.log; sorted split.tsv
check "an unsorted file of all read 1 then all read 2 gives the same table with a small hold limit" 'cmp -s split.tsv.sorted fq.tsv.sorted && grep -q "keeping only reads that hit" split.log'

sed 's/SO:queryname/SO:coordinate/' name.sam > liar.sam
"$BIN" recruit --locus "$IDX" -i liar.sam -o x.tsv -q 2> err3.txt; rc=$?
check "a file claiming coordinate order without it stops the run" '[[ $rc -ne 0 ]] && grep -q "not sorted by coordinate" err3.txt'

if command -v samtools > /dev/null; then
  samtools view -b -o coord.bam coord.sam && samtools index coord.bam
  "$BIN" recruit --locus "$IDX" -i coord.bam -o bam.tsv -t 4 2> bam.log; sorted bam.tsv
  check "a coordinate-sorted BAM with index gives the same table, mates fetched through the index" 'cmp -s bam.tsv.sorted fq.tsv.sorted && grep -q "mates fetched" bam.log && ! grep -q "second pass" bam.log'
  # mate fields pointing to the wrong place: the index misses, the second pass finds the mates
  awk -F'\t' 'BEGIN {OFS = "\t"} !/^@/ && $8 > 150000 {$8 += 7} {print}' coord.sam | samtools view -b -o stale.bam - && samtools index stale.bam
  "$BIN" recruit --locus "$IDX" -i stale.bam -o stale.tsv -t 4 2> stale.log; sorted stale.tsv
  check "stale mate positions: the mates are found by a second pass, same table" 'cmp -s stale.tsv.sorted fq.tsv.sorted && grep -q "not where the index put them" stale.log'
  for c in chrA chrB; do printf '>%s\n' $c; head -c 20000000 /dev/zero | tr '\0' 'A' | fold -w 60; echo; done > zero.fa
  samtools faidx zero.fa && samtools view -C -T zero.fa -o coord.cram coord.bam 2> /dev/null && samtools index coord.cram
  "$BIN" recruit --locus "$IDX" -i coord.cram --cram-reference zero.fa -o cram.tsv -t 4 -q; sorted cram.tsv
  check "a CRAM gives the same table" 'cmp -s cram.tsv.sorted fq.tsv.sorted'
  python3 -c "import sys; d=open('coord.bam','rb').read(); import struct; o=0; ends=[]
while o < len(d):
    bsize=struct.unpack('<H', d[o+16:o+18])[0]+1; o+=bsize; ends.append(o)
open('cut.bam','wb').write(d[:ends[len(ends)//2]])"
  "$BIN" recruit --locus "$IDX" -i cut.bam -o x.tsv -q 2> err4.txt; rc=$?
  check "a BAM cut at a block boundary is refused as truncated" '[[ $rc -ne 0 ]] && grep -q "truncated" err4.txt'
else
  echo "  skip BAM, CRAM and truncation checks (samtools not installed)"
fi

# recruited reads as FASTQ, the summary, a preset
"$BIN" recruit --locus "$IDX" -1 r_1.fq -2 r_2.fq -o out.tsv --reads-out reads_out --preset illumina -q
check "certain pairs are written as FASTQ, mates in step" '[[ $(gzip -dc reads_out/locus.certain_1.fq.gz | awk "NR % 4 == 1" | wc -l) -eq $(awk -F"\t" '"'"'$3=="certain"'"'"' out.tsv | wc -l) && $(gzip -dc reads_out/locus.certain_2.fq.gz | wc -l) -eq $(gzip -dc reads_out/locus.certain_1.fq.gz | wc -l) ]]'
check "the summary gives the error rate and the backbone depth" 'grep -q "^run.error_rate.0.00200.given" out.summary.tsv && grep -q "backbone_depth" out.summary.tsv'
"$BIN" recruit --locus "$IDX" -1 r_1.fq -2 r_2.fq -o est.tsv 2> est.log
check "without a rate, it is estimated from the base qualities" 'grep -q "error rate 0.0010 (estimated from" est.log'

# more loci than 64: seventy copies of the locus, each read scored for every one
mkdir -p many && for i in $(seq -w 1 70); do cp data/locus.gfa many/l$i.gfa; done
"$BIN" recruit -t 4 -q --background s -r ref --index-only --save-index many_idx $(for i in $(seq -w 1 70); do printf -- "--locus many/l%s.gfa " $i; done)
"$BIN" recruit -t 4 -q $(for f in many_idx/*.plx; do printf -- "--locus %s " $f; done) -1 a_1.fq -2 a_2.fq -o many.tsv
check "seventy loci in one run: a certain pair is certain in all seventy, the others listed" 'awk -F"\t" '"'"'NR > 1 && $3 == "certain" && $24 == 70 {n++} NR > 1 && $3 == "certain" && split($23, o, ",") != 69 {bad++} END {exit !(n > 0 && bad == 0)}'"'"' many.tsv'

"$BIN" recruit -i name.sam -1 r_1.fq -2 r_2.fq --dry-run > dry.tsv
check "--dry-run counts records by kind" 'awk -F"\t" '"'"'$2=="name.sam" && $3 > 9000 && $4 > 0 && $5 > 0 && $6 == 50 && $10 == 9000'"'"' dry.tsv | grep -q . && awk -F"\t" '"'"'$2=="r_1.fq" && $3 == 4500'"'"' dry.tsv | grep -q .'

printf 's1\tr_1.fq\tr_2.fq\ns2\tname.sam\ns2\tinter.fq\n' > sheet.tsv
"$BIN" recruit --locus "$IDX" --samples sheet.tsv -o per_sample -q
check "a sample sheet writes one table per sample, files of a sample pooled" '[[ -f per_sample/s1.tsv && $(tail -n +2 per_sample/s2.tsv | wc -l) -gt $(tail -n +2 per_sample/s1.tsv | wc -l) ]]'

sed -n '5,$p' r_2.fq > shifted_2.fq
"$BIN" recruit --locus "$IDX" -1 r_1.fq -2 shifted_2.fq -o x.tsv -q 2> err1.txt; rc=$?
check "mates out of step stop the run" '[[ $rc -ne 0 ]] && grep -q "out of step\|fewer reads" err1.txt'
head -c 300000 name.sam > cut.sam
"$BIN" recruit --locus "$IDX" -i cut.sam -o x.tsv -q 2> err2.txt; rc=$?
check "a damaged SAM stops the run, the file named" '[[ $rc -ne 0 ]] && grep -q "cut.sam" err2.txt'

echo
echo "recruit_input: $pass passed, $fail failed"
[[ $fail -eq 0 ]]
