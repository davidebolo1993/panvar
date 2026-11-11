#!/usr/bin/env bash
# panvar recruit --index-only on a synthetic locus whose reference carries a diverged paralog: outside
# counts against brute force over the genomes, classes from them, decoys from the paralog, the screen,
# calibration and the cut it gives recruit, reuse of a current index and rebuilding when the graph,
# background or settings change.
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

# reads: from the locus, from its paralog, from elsewhere in the reference
python3 "$TOOLS/locus_reads.py" data/ref.fa chr1 120000 160000 2000 loc 1 reads_locus
python3 "$TOOLS/locus_reads.py" data/ref.fa chr1 231000 259000 2000 par 2 reads_paralog
python3 "$TOOLS/locus_reads.py" data/ref.fa chr1 0 100000 2000 far 3 reads_far
cat reads_locus_1.fq reads_paralog_1.fq reads_far_1.fq > reads_1.fq
cat reads_locus_2.fq reads_paralog_2.fq reads_far_2.fq > reads_2.fq
"$BIN" recruit --locus idx/locus.plx -1 reads_1.fq -2 reads_2.fq -o recruited.tsv -t 2 -q; rc=$?
fate() { awk -F'\t' -v p="$1" -v f="$2" 'NR>1 && index($1, p)==1 && $3==f' recruited.tsv | wc -l | tr -d ' '; }
check "reads from the locus are certain ($(fate loc certain) of 2000)" '[[ $rc -eq 0 && $(fate loc certain) -ge 1900 ]]'
check "reads from the paralog are never certain, and mostly unrecruited as homolog ($(fate par unrecruited) of 2000)" '[[ $(fate par certain) == 0 && $(fate par unrecruited) -ge 1000 ]]'
check "reads from elsewhere are not even screened" '[[ $(awk -F"\t" '"'"'NR>1 && index($1, "far")==1'"'"' recruited.tsv | wc -l | tr -d " ") == 0 ]]'
# across the locus's end (its last few hundred bases are shared: the copies' ends are uncertain): pairs of
# 1400 bp fragments with one mate well inside and the other past the end; single 2 kb reads past the end
python3 "$TOOLS/locus_reads.py" data/ref.fa chr1 158700 160400 300 edge 4 reads_edge 1400
"$BIN" recruit --locus idx/locus.plx -1 reads_edge_1.fq -2 reads_edge_2.fq -o edge.tsv -t 2 -q
python3 "$TOOLS/locus_reads.py" data/ref.fa chr1 158600 161000 300 long 5 reads_long 2000 2000
"$BIN" recruit --locus idx/locus.plx -1 reads_long_1.fq -o long.tsv -t 2 -q
reason() { awk -F'\t' -v r="$2" 'NR>1 && $3=="certain" && $4==r' "$1" | wc -l | tr -d ' '; }
check "pairs with a mate past the end are certain, rescued ($(reason edge.tsv rescued) of 300)" '[[ $(reason edge.tsv rescued) -ge 270 ]]'
python3 "$TOOLS/locus_reads.py" data/ref.fa chr1 115000 165000 30 span 6 reads_span 48000 48000
"$BIN" recruit --locus idx/locus.plx -1 reads_span_1.fq -o span.tsv -t 2 -q
check "reads holding the whole locus between two flanks are certain, partial ($(reason span.tsv partial) of 30)" '[[ $(reason span.tsv partial) == 30 ]]'
check "reads past the end are certain, partial ($(reason long.tsv partial) of 300), segments within the read" '[[ $(reason long.tsv partial) -ge 270 ]] && awk -F"\t" '"'"'NR>1 && $4=="partial" {split($6, g, "-"); if (g[2] - g[1] >= 2000 || g[1] != 0 && g[2] != 2000) bad = 1} END {exit bad}'"'"' long.tsv'

# a graph genome counts one haplotype per sample#haplotype of its paths
"$BIN" background -t 2 --partitions 64 --memory 64M -q --tmp-dir "$OUT/tmp" -o s_graph --genome ref=data/ref.fa --genome lg=data/locus.gfa --coordinate-reference ref
"$BIN" recruit -t 2 -q --background s_graph -r ref --index-only --save-index idx_graph --locus data/locus.gfa 2> graph.txt; rc=$?
check "a graph genome's PanSN paths count as haplotypes (1 FASTA + 3 paths)" '[[ $rc -eq 0 ]] && grep -q "haplotypes.4 in the background" idx_graph/locus.index.tsv'
cp data/ga.fa data/ga_moved.fa
"$BIN" background -t 2 --partitions 64 --memory 64M -q --tmp-dir "$OUT/tmp" -o s_moved --genome ref=data/ref.fa --genome ga=data/ga_moved.fa --coordinate-reference ref
printf '>chr1\nACGT\n' >> data/ga_moved.fa
"$BIN" recruit -t 2 -q --background s_moved -r ref --index-only --save-index idx_moved --locus data/locus.gfa 2> moved.txt; rc=$?
check "a genome file changed since it was counted is refused" '[[ $rc -ne 0 ]] && grep -q "was replaced" moved.txt'

# the decoy scan split across threads gives the same index
"$BIN" recruit -t 4 -q --background s -r ref --index-only --save-index idx4 --locus data/locus.gfa
"$BIN" recruit --dump-index idx4/locus.plx=dump4t.tsv
check "an index built on 4 threads equals the 1-thread one" 'cmp -s dump.tsv dump4t.tsv'
cal() { awk -F'\t' '$1=="calibration" && $2!="class"' "$1"; }
check "calibration: three read classes at cuts 1 to 8" '[[ $(cal idx/locus.index.tsv | wc -l | tr -d " ") == 24 && $(cal idx/locus.index.tsv | cut -f2 | sort -u | tr "\n" " ") == "hifi ont short " ]]'
check "calibration: certain counts within their totals, falling as the cut rises" 'cal idx/locus.index.tsv | awk -F"\t" '"'"'{if ($8 > $7 || $11 > $10 || $13 > $12) bad = 1; if ($6 > 1 && ($8 > l8 || $11 > l11 || $13 > l13)) bad = 1; l8 = $8; l11 = $11; l13 = $13} END {exit bad}'"'"''
held=$(cal idx/locus.index.tsv | awk -F'\t' '$2=="short" && $6==2 {print $8 " of " $7}')
check "calibration: the paralog gives off-target reads, the held-out sample locus reads ($held certain at 2 sites)" 'cal idx/locus.index.tsv | awk -F"\t" '"'"'$2=="short" && $6==2 && $12 > 1000 && $7 > 500 && $8 / $7 > 0.9 {ok = 1} END {exit !ok}'"'"''
check "calibration on 4 threads equals the 1-thread one" '[[ "$(cal idx/locus.index.tsv)" == "$(cal idx4/locus.index.tsv)" ]]'
check "recruit takes each locus's calibrated cut and reports it" 'grep -q "^locus.locus.site_cut.2.calibrated (short reads)" recruited.summary.tsv'
"$BIN" recruit --locus idx/locus.plx -1 reads_1.fq -2 reads_2.fq -o fixed.tsv -t 2 -q --min-sites 3
check "--min-sites gives every locus that cut" 'grep -q "^locus.locus.site_cut.3.given" fixed.summary.tsv'
"$BIN" recruit --locus idx/locus.plx -1 reads_1.fq -2 reads_2.fq -o both.tsv -q --min-sites 3 --false-certain-rate 1e-3 2> both.txt; rc=$?
check "--min-sites and --false-certain-rate together are refused" '[[ $rc -ne 0 ]] && grep -q "not both" both.txt'
"$BIN" recruit -t 2 -q --background s -r ref --index-only --save-index idx_nocal --locus data/locus.gfa --no-calibration
"$BIN" recruit --locus idx_nocal/locus.plx -1 reads_1.fq -2 reads_2.fq -o nocal.tsv -t 2 -q
check "an index without calibration: no calibration rows, the cut 2, the same table" '[[ -z "$(cal idx_nocal/locus.index.tsv)" ]] && grep -q "site_cut.2.the index carries no calibration" nocal.summary.tsv && cmp -s recruited.tsv nocal.tsv'
mkdir -p twin/a twin/b && cp data/locus.gfa twin/a/locus.gfa && cp data/locus.gfa twin/b/locus.gfa
"$BIN" recruit -t 2 -q --background s -r ref --index-only --save-index idx_twin --locus twin/a/locus.gfa --locus twin/b/locus.gfa 2> twin.txt; rc=$?
check "two graphs of one name are refused" '[[ $rc -ne 0 ]] && grep -q "rename one graph" twin.txt'

"${RC[@]}" --locus data/locus.gfa 2> log2.txt
check "a current index is reused" 'grep -q "is current, reused" log2.txt'
"${RC[@]}" --locus data/locus.gfa --no-calibration 2> log2b.txt
check "an index asked without calibration is rebuilt" 'grep -q "rebuilding" log2b.txt'
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
