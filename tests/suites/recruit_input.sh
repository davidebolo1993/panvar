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

"$BIN" recruit -i name.sam -1 r_1.fq -2 r_2.fq --dry-run > dry.tsv
check "--dry-run counts records by kind" 'awk -F"\t" '"'"'$2=="name.sam" && $3 > 9000 && $4 > 0 && $5 > 0 && $6 == 50 && $10 == 4500'"'"' dry.tsv | grep -q . && awk -F"\t" '"'"'$2=="r_1.fq" && $3 == 4500'"'"' dry.tsv | grep -q .'

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
