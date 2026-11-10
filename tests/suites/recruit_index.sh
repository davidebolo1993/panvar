#!/usr/bin/env bash
# panvar recruit --index-only on a synthetic locus whose reference carries a diverged paralog: outside
# counts against brute force over the genomes, classes from them, decoys from the paralog, the screen,
# reuse of a current index and rebuilding when the graph, background or settings change.
#   recruit_index.sh <panvar_bin> <tests/tools dir> <out_dir>
set -uo pipefail

if [[ $# -lt 3 ]]; then
  echo "Usage: $0 <panvar_bin> <tools_dir> <out_dir>" >&2
  exit 2
fi
BIN="$1"; TOOLS="$2"; OUT="$3"
rm -rf "$OUT"; mkdir -p "$OUT/tmp"; cd "$OUT" || exit 2
python3 "$TOOLS/locus_fixture.py" data
BG=("$BIN" background -t 2 --partitions 64 --memory 64M -q --tmp-dir "$OUT/tmp")
RC=("$BIN" recruit -t 2 --background s -r ref --index-only --save-index idx)

pass=0; fail=0
ok() { echo "  ok   $1"; pass=$((pass + 1)); }
bad() { echo "  FAIL $1"; fail=$((fail + 1)); }
check() { if eval "$2"; then ok "$1"; else bad "$1"; fi; }
field() { awk -F'\t' -v k="$2" '$1=="locus" && $2==k {print $3}' "$1"; }

"${BG[@]}" -o s --genome ref=data/ref.fa --genome ga=data/ga.fa --genome gb=data/gb.fa --coordinate-reference ref
"${RC[@]}" --locus data/locus.gfa 2> log1.txt; rc=$?
check "the index builds" '[[ $rc -eq 0 && -f idx/locus.plx && -f idx/locus.index.tsv ]]'
"$BIN" recruit --dump-index idx/locus.plx=dump.tsv
python3 "$TOOLS/locus_outside.py" idx/locus.index.tsv dump.tsv 8 data/ref.fa data/ga.fa data/gb.fa > brute.txt; rc=$?
check "outside counts equal brute force, classes follow, decoys are not locus k-mers ($(tail -1 brute.txt))" '[[ $rc -eq 0 ]]'
check "the interval is read from the reference path" '[[ "$(field idx/locus.index.tsv interval)" == "chr1:120000-160000" ]]'
check "the locus is found in every genome" '[[ $(awk -F"\t" '"'"'$1=="span" && $2!="genome"'"'"' idx/locus.index.tsv | cut -f2 | sort -u | wc -l | tr -d " ") == 3 ]]'
check "the paralog makes shared k-mers and decoys" '[[ $(field idx/locus.index.tsv shared | cut -d" " -f1) -gt 1000 && $(field idx/locus.index.tsv decoys) -gt 1000 ]]'
check "the screen holds only private and shared k-mers, about one in ten" 'awk -F"\t" '"'"'NR>1 && $2!="D" {n++; s+=$6} END {exit !(s > n/20 && s < n/5)}'"'"' dump.tsv'
check "ceilings per read length" 'awk -F"\t" '"'"'$1=="ceiling" && $2==150 && $4>0.9 && $4<=1'"'"' idx/locus.index.tsv | grep -q .'

"${RC[@]}" --locus data/locus.gfa 2> log2.txt
check "a current index is reused" 'grep -q "is current, reused" log2.txt'
"${RC[@]}" --locus data/locus.gfa --informative-max 2 2> log3.txt
check "another setting rebuilds it" 'grep -q "rebuilding" log3.txt && [[ $(field idx/locus.index.tsv repeat | cut -d" " -f1) -ge 0 ]]'
{ sed "s/>chr1/>x1/" data/ga.fa; sed "s/>chr1/>x2/" data/gb.fa; } > data/gc.fa
"${BG[@]}" -o s --add --genome gc=data/gc.fa
"${RC[@]}" --locus data/locus.gfa --informative-max 2 2> log4.txt
check "a background that counts more genomes rebuilds it" 'grep -q "rebuilding" log4.txt'
"$BIN" recruit --dump-index idx/locus.plx=dump4.tsv
python3 "$TOOLS/locus_outside.py" idx/locus.index.tsv dump4.tsv 2 data/ref.fa data/ga.fa data/gb.fa data/gc.fa > brute4.txt; rc=$?
check "after it, outside counts still equal brute force ($(tail -1 brute4.txt))" '[[ $rc -eq 0 ]]'
"${RC[@]}" --locus data/locus.gfa -r nosuchpath > /dev/null 2>&1; rc=$?
check "a reference path the graph lacks is refused" '[[ $rc -ne 0 ]]'

# the same genomes indexed (.fai), and bgzipped with .gzi, read by interval: the same index
if command -v samtools > /dev/null && command -v bgzip > /dev/null; then
  mkdir -p z
  for g in ref ga gb; do
    cp data/$g.fa z/$g.fa && samtools faidx z/$g.fa
    bgzip -c data/$g.fa > z/$g.fa.gz && samtools faidx z/$g.fa.gz
  done
  for kind in fa fa.gz; do
    "${BG[@]}" -o s_$kind --genome ref=z/ref.$kind --genome ga=z/ga.$kind --genome gb=z/gb.$kind --coordinate-reference ref
    "$BIN" recruit -t 2 -q --background s_$kind -r ref --index-only --save-index idx_$kind --locus data/locus.gfa
    "$BIN" recruit --dump-index idx_$kind/locus.plx=dump_$kind.tsv
  done
  "${BG[@]}" -o s0 --genome ref=data/ref.fa --genome ga=data/ga.fa --genome gb=data/gb.fa --coordinate-reference ref
  "$BIN" recruit -t 2 -q --background s0 -r ref --index-only --save-index idx0 --locus data/locus.gfa
  "$BIN" recruit --dump-index idx0/locus.plx=dump0.tsv
  check "genomes read through .fai, or bgzip and .gzi, give the same index" 'cmp -s dump0.tsv dump_fa.tsv && cmp -s dump0.tsv dump_fa.gz.tsv'
else
  echo "  skip indexed and bgzip genomes (samtools or bgzip not installed)"
fi

echo
echo "recruit_index: $pass passed, $fail failed"
[[ $fail -eq 0 ]]
