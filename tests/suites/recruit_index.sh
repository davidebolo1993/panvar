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

# a reference path whose label is not where its sequence is, or not its length, is refused
sed 's/^P\tref#1#chr1:120000-160000/P\tref#1#chr1:125000-165000/' data/locus.gfa > shifted.gfa
"${RC[@]}" --save-index idx_shift --locus shifted.gfa 2> shift.txt; rc=$?
check "a reference label off its sequence is refused" '[[ $rc -ne 0 ]] && grep -q "label does not match the reference" shift.txt'
# a label on another path's sequence: the interval holds the locus's k-mers, but not the reference path's
python3 - data/ref.fa > swapped.fa <<'PY'
import sys
s = []
for line in open(sys.argv[1]):
    if line.startswith(">"):
        if s: break
        continue
    s.append(line.strip())
r = "".join(s)
print(">ref#1#chr1:120000-160000"); print(r[121000:161000])
print(">gx#1#chr1:120000-160000"); print(r[120000:160000])
PY
"${RC[@]}" --save-index idx_swap --locus swapped.fa 2> swap.txt; rc=$?
check "a reference label on another path's sequence is refused" '[[ $rc -ne 0 ]] && grep -q "does not spell the path" swap.txt'
sed 's/^P\tref#1#chr1:120000-160000/P\tref#1#chr1:120000-150000/' data/locus.gfa > short_label.gfa
"${RC[@]}" --save-index idx_len --locus short_label.gfa 2> len.txt; rc=$?
check "a reference label of another length is refused" '[[ $rc -ne 0 ]] && grep -q "length is not the path" len.txt'
# a link that overlaps would join segments into k-mers no genome has: refused
awk -F'\t' 'BEGIN {OFS = "\t"} $1 == "S" && $2 == 1 {print "S", "1", substr($3, 1, 20000); print "S", "1b", substr($3, 19971); next}
  $1 == "P" && $2 ~ /^ref/ {print "P", $2, "1+,1b+", "*"; next} {print}
  END {print "L", "1", "+", "1b", "+", "30M"}' data/locus.gfa > overlap.gfa
"${RC[@]}" --save-index idx_ovl --locus overlap.gfa 2> ovl.txt; rc=$?
check "a link with an overlap is refused" '[[ $rc -ne 0 ]] && grep -q "overlaps by 30 bases" ovl.txt'
# a short locus (fewer anchors inside than a chain needs): copies only between flank anchors, and only the
# interval projected from them taken as certain
python3 - data/ref.fa data/ga.fa > short.gfa <<'PY'
import sys
def chrom(path):
    s = []
    for line in open(path):
        if line.startswith(">"):
            if s: break
            continue
        s.append(line.strip())
    return "".join(s)
r = chrom(sys.argv[1])[130000:130600]
print("H\tVN:Z:1.0"); print(f"S\t1\t{r}"); print("P\tref#1#chr1:130000-130600\t1+\t*")
PY
"${RC[@]}" --save-index idx_short --locus short.gfa 2> short.txt; rc=$?
check "a short locus builds and is found in the other genomes" '[[ $rc -eq 0 && $(awk -F"\t" '"'"'$1=="span" && $2!="genome" && $2!="ref"'"'"' idx_short/short.index.tsv | wc -l) -ge 2 ]]'
check "a short locus's copies are certain only around the projected interval" 'awk -F"\t" '"'"'$1=="span" && $2!="genome" && $2!="ref" {n++; if ($9 - $8 > 0.2 * ($5 - $4)) bad++} END {exit !(n > 0 && !bad)}'"'"' idx_short/short.index.tsv'
# one locus failing late: the others' indexes are still saved, the run fails naming it
PANVAR_TEST_FAIL_LOCUS=short "${RC[@]}" --save-index idx_two --locus data/locus.gfa --locus short.gfa 2> two.txt; rc=$?
check "a locus failing late leaves the others' indexes saved" '[[ $rc -ne 0 && -f idx_two/locus.plx && ! -f idx_two/short.plx ]] && grep -q "the other loci" two.txt'
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
# a FASTA locus, the graph's paths spelled out one record each (the fixture's steps are all forward): the graph's
# index but for the nodes, and the same calls
mkdir -p fa
awk -F'\t' '$1 == "S" {s[$2] = $3} $1 == "P" {n = split($3, st, ","); q = ""; for (i = 1; i <= n; i++) q = q s[substr(st[i], 1, length(st[i]) - 1)]
  print ">" $2; print q}' data/locus.gfa > fa/locus.fa
"${RC[@]}" --save-index idx_fa --locus fa/locus.fa; rc=$?
"$BIN" recruit --dump-index idx_fa/locus.plx=dump_fa.tsv
check "a FASTA locus gives its graph's k-mers, classes and counts, without nodes" '[[ $rc -eq 0 ]] && cmp -s <(cut -f1-6 dump.tsv) <(cut -f1-6 dump_fa.tsv) && [[ $(tail -n +2 dump_fa.tsv | cut -f7 | sort -u) == "." ]]'
check "and its report and calibration but for the graph and its node sets" 'diff <(grep -vE "^locus	(graph|node_sets)	" idx/locus.index.tsv) <(grep -vE "^locus	(graph|node_sets)	" idx_fa/locus.index.tsv) > /dev/null'
"$BIN" recruit --locus idx_fa/locus.plx -1 reads_1.fq -2 reads_2.fq -o recruited_fa.tsv -t 2 -q
nodes_col=$(head -1 recruited.tsv | tr '\t' '\n' | grep -n '^nodes$' | cut -d: -f1)
blank() { awk -v c="$nodes_col" 'BEGIN {FS = OFS = "\t"} {$c = ""; print}' "$1"; }
check "and the same calls, nodes \".\"" 'cmp -s <(blank recruited.tsv) <(blank recruited_fa.tsv) && [[ $(tail -n +2 recruited_fa.tsv | cut -f"$nodes_col" | sort -u) == "." ]]'
{ cat fa/locus.fa; printf ">ga#1#chr1:120000-160000\n"; sed -n 2p fa/locus.fa; } > fa/dup.fa
"${RC[@]}" --save-index idx_dup --locus fa/dup.fa 2> dup.txt; rc=$?
check "a FASTA locus with two records of one name is refused" '[[ $rc -ne 0 ]] && grep -q "two records named ga#1#chr1:120000-160000" dup.txt'
printf 'not a locus\n' > fa/text.fa
"${RC[@]}" --save-index idx_text --locus fa/text.fa 2> text.txt; rc=$?
check "a file neither GFA nor FASTA is refused" '[[ $rc -ne 0 ]] && grep -q "neither a GFA graph nor a FASTA file" text.txt'
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
# a long read through the locus with the paralog (a homologous flank, its decoys) beyond one end and one stray
# paralog 40-mer beyond the other: the stray k-mers are no flank, so the locus stretch alone is judged
python3 - data/ref.fa > reads_flank.fq <<'PY'
import sys
seq = "".join(l.strip() for l in open(sys.argv[1]) if not l.startswith(">"))
locus = seq[120000:160000]
stray = next(seq[a:a + 40] for a in range(232000, 258000, 7) if seq[a:a + 40] not in locus)
read = seq[100000:110000] + stray + seq[110000:260000]
print("@flank\n" + read + "\n+\n" + "I" * len(read))
PY
"$BIN" recruit --locus idx/locus.plx -i reads_flank.fq -o flank.tsv -t 2 -q
check "a stray decoy k-mer beyond one end does not make a long read judged whole ($(awk -F'\t' 'NR > 1 {print $3 "/" $4}' flank.tsv))" '[[ "$(awk -F"\t" '"'"'NR > 1 {print $3}'"'"' flank.tsv)" == certain ]]'

# a locus copy broken inside a tandem array: the array's middle sits on a contig of its own, out of the
# anchors' reach; it is found as a piece and subtracted, and the copies' uncertain ends that are array count
# as the locus. A paralogous array in another genome is not a piece.
"$BIN" background -t 2 --partitions 64 --memory 64M -q --tmp-dir "$OUT/tmp" -o s_split --genome ref=data/split/ref.fa --genome gs=data/split/gs.fa --coordinate-reference ref
"$BIN" recruit -t 2 -q --background s_split -r ref --index-only --save-index idx_split --locus data/split/locus.gfa; rc=$?
"$BIN" recruit --dump-index idx_split/locus.plx=dump_split.tsv
pieces() { awk -F'\t' '$1=="piece" && $2!="genome" {print $2 ":" $3 ":" $4 "-" $5}' "$1" | tr '\n' ' '; }
check "a broken copy's orphan piece is found ($(pieces idx_split/locus.index.tsv))" '[[ $rc -eq 0 && "$(pieces idx_split/locus.index.tsv)" == "gs:x:0-17700 " ]]'
check "outside counts equal brute force with the piece subtracted" 'python3 "$TOOLS/locus_outside.py" idx_split/locus.index.tsv dump_split.tsv 8 data/split/ref.fa data/split/gs.fa > brute_split.txt'
check "the array's k-mers stay private, none demoted at the copies' ends" '[[ $(field idx_split/locus.index.tsv span_demoted) == 0 ]] && awk -F"\t" '"'"'NR>1 && $2!="D" && $4>=3 {n++; p+=($2=="P")} END {exit !(n>5000 && p==n)}'"'"' dump_split.tsv'
"$BIN" background -t 2 --partitions 64 --memory 64M -q --tmp-dir "$OUT/tmp" -o s_split_p --genome ref=data/split/ref.fa --genome gs=data/split/gs.fa --genome gp=data/split/gp.fa --coordinate-reference ref
"$BIN" recruit -t 2 -q --background s_split_p -r ref --index-only --save-index idx_split_p --locus data/split/locus.gfa; rc=$?
"$BIN" recruit --dump-index idx_split_p/locus.plx=dump_split_p.tsv
check "a paralogous array is not a piece; its k-mers count outside" '[[ $rc -eq 0 && "$(pieces idx_split_p/locus.index.tsv)" == "gs:x:0-17700 " ]] && python3 "$TOOLS/locus_outside.py" idx_split_p/locus.index.tsv dump_split_p.tsv 8 data/split/ref.fa data/split/gs.fa data/split/gp.fa > brute_split_p.txt && awk -F"\t" '"'"'NR>1 && $2=="S" && $4>=3 {n++} END {exit !(n>5000)}'"'"' dump_split_p.tsv'
# a genome holding a piece named as a decoy genome: the piece is the locus there, as its copies are, so it
# gives no decoy window and calibration reads none of it as off-target
"$BIN" recruit -t 2 -q --background s_split_p -r ref --index-only --save-index idx_split_d --locus data/split/locus.gfa --decoy-genomes ref,gs,gp; rc=$?
offcut2() { awk -F'\t' '$1=="calibration" && $2=="short" && $6==2 {print $13 " of " $12}' "$1"; }
check "a piece is no decoy window, nor off-target in calibration (short off-target certain at 2 sites: $(offcut2 idx_split_d/locus.index.tsv))" '[[ $rc -eq 0 && "$(pieces idx_split_d/locus.index.tsv)" == "gs:x:0-17700 " ]] && ! awk -F"\t" '"'"'$1=="decoy_window" && $2=="gs" && $3=="x"'"'"' idx_split_d/locus.index.tsv | grep -q . && [[ "$(offcut2 idx_split_d/locus.index.tsv)" == 0\ of* ]]'
# a genome without a copy: a paralog dense in the locus's k-mers but with other sequence on both sides, on a
# contig without a copy, is not a piece; a contig of the array alone is
"$BIN" background -t 2 --partitions 64 --memory 64M -q --tmp-dir "$OUT/tmp" -o s_split_c --genome ref=data/split/ref.fa --genome gc=data/split_c/gc.fa --coordinate-reference ref
"$BIN" recruit -t 2 -q --background s_split_c -r ref --index-only --save-index idx_split_c --locus data/split/locus.gfa; rc=$?
"$BIN" recruit --dump-index idx_split_c/locus.plx=dump_split_c.tsv
check "a paralog mid-contig is not a piece, an array contig is ($(pieces idx_split_c/locus.index.tsv))" '[[ $rc -eq 0 && "$(pieces idx_split_c/locus.index.tsv)" == "gc:o:0-16500 " ]] && python3 "$TOOLS/locus_outside.py" idx_split_c/locus.index.tsv dump_split_c.tsv 8 data/split/ref.fa data/split_c/gc.fa > brute_split_c.txt && awk -F"\t" '"'"'NR>1 && $2=="S" && $4>=3 {n++} END {exit !(n>5000)}'"'"' dump_split_c.tsv'
# a contig shorter than a window, of the array alone, is a piece whole
"$BIN" background -t 2 --partitions 64 --memory 64M -q --tmp-dir "$OUT/tmp" -o s_split_h --genome ref=data/split/ref.fa --genome gh=data/split_h/gh.fa --coordinate-reference ref
"$BIN" recruit -t 2 -q --background s_split_h -r ref --index-only --save-index idx_split_h --locus data/split/locus.gfa; rc=$?
"$BIN" recruit --dump-index idx_split_h/locus.plx=dump_split_h.tsv
check "a short array contig is a piece ($(pieces idx_split_h/locus.index.tsv))" '[[ $rc -eq 0 && "$(pieces idx_split_h/locus.index.tsv)" == "gh:s:0-4500 " ]] && python3 "$TOOLS/locus_outside.py" idx_split_h/locus.index.tsv dump_split_h.tsv 8 data/split/ref.fa data/split_h/gh.fa > brute_split_h.txt && awk -F"\t" '"'"'NR>1 && $2!="D" && $4>=3 {n++; p+=($2=="P")} END {exit !(n>5000 && p==n)}'"'"' dump_split_h.tsv'
# a much longer array on one contig: the copy's end is estimated and stops inside it; the contig is scanned
"$BIN" background -t 2 --partitions 64 --memory 64M -q --tmp-dir "$OUT/tmp" -o s_tail --genome ref=data/split_t/ref.fa --genome gl=data/split_t/gl.fa --coordinate-reference ref
"$BIN" recruit -t 2 -q --background s_tail -r ref --index-only --save-index idx_tail --locus data/split_t/locus.gfa; rc=$?
"$BIN" recruit --dump-index idx_tail/locus.plx=dump_tail.tsv
check "an array past a copy's estimated end is found on its contig ($(pieces idx_tail/locus.index.tsv))" '[[ $rc -eq 0 && "$(pieces idx_tail/locus.index.tsv)" == "gl:chr1:161970-228000 " ]] && python3 "$TOOLS/locus_outside.py" idx_tail/locus.index.tsv dump_tail.tsv 8 data/split_t/ref.fa data/split_t/gl.fa > brute_tail.txt && awk -F"\t" '"'"'NR>1 && $2!="D" && $4>=3 {n++; p+=($2=="P")} END {exit !(n>5000 && p==n)}'"'"' dump_tail.tsv'
# a whole copy, and more of its array on a contig of its own: the tandem k-mers left with a few copies
# outside send the index looking in every genome; a locus without an array never does
"$BIN" background -t 2 --partitions 64 --memory 64M -q --tmp-dir "$OUT/tmp" -o s_orphan --genome ref=data/split_o/ref.fa --genome go=data/split_o/go.fa --coordinate-reference ref
"$BIN" recruit -t 2 --background s_orphan -r ref --index-only --save-index idx_orphan --locus data/split_o/locus.gfa 2> orphan.txt; rc=$?
"$BIN" recruit --dump-index idx_orphan/locus.plx=dump_orphan.tsv
check "an array contig next to a whole copy is found ($(pieces idx_orphan/locus.index.tsv))" '[[ $rc -eq 0 && "$(pieces idx_orphan/locus.index.tsv)" == "go:orphan:0-22000 " ]] && grep -q "every genome scanned for pieces" orphan.txt && python3 "$TOOLS/locus_outside.py" idx_orphan/locus.index.tsv dump_orphan.tsv 8 data/split_o/ref.fa data/split_o/go.fa > brute_orphan.txt && awk -F"\t" '"'"'NR>1 && $2!="D" && $4>=3 {n++; p+=($2=="P")} END {exit !(n>5000 && p==n)}'"'"' dump_orphan.tsv'
check "a locus without an array is not looked for in every genome" '! grep -q "every genome scanned" log1.txt'
# arrays running into contig ends through stubs with errors: a piece and a copy's uncertain end reach the end
"$BIN" background -t 2 --partitions 64 --memory 64M -q --tmp-dir "$OUT/tmp" -o s_ends --genome ref=data/split_e/ref.fa --genome ge=data/split_e/ge.fa --genome gt=data/split_e/gt.fa --coordinate-reference ref
"$BIN" recruit -t 2 -q --background s_ends -r ref --index-only --save-index idx_ends --locus data/split_e/locus.gfa; rc=$?
"$BIN" recruit --dump-index idx_ends/locus.plx=dump_ends.tsv
check "array stubs at contig ends count as the locus ($(pieces idx_ends/locus.index.tsv))" '[[ $rc -eq 0 && "$(pieces idx_ends/locus.index.tsv)" == "ge:o:0-23700 " && $(field idx_ends/locus.index.tsv span_demoted) == 0 ]] && python3 "$TOOLS/locus_outside.py" idx_ends/locus.index.tsv dump_ends.tsv 8 data/split_e/ref.fa data/split_e/ge.fa data/split_e/gt.fa > brute_ends.txt && awk -F"\t" '"'"'NR>1 && $2!="D" {n++; p+=($2=="P")} END {exit !(n>50000 && p==n)}'"'"' dump_ends.tsv'
# more of the array mid-way along the copy's contig, far from the copy: no piece, its k-mers count outside
"$BIN" background -t 2 --partitions 64 --memory 64M -q --tmp-dir "$OUT/tmp" -o s_far --genome ref=data/split_f/ref.fa --genome gf=data/split_f/gf.fa --coordinate-reference ref
"$BIN" recruit -t 2 --background s_far -r ref --index-only --save-index idx_far --locus data/split_f/locus.gfa 2> far.txt; rc=$?
"$BIN" recruit --dump-index idx_far/locus.plx=dump_far.tsv
check "an array far from the copy on its contig is no piece ($(pieces idx_far/locus.index.tsv))" '[[ $rc -eq 0 && -z "$(pieces idx_far/locus.index.tsv)" ]] && grep -q "every genome scanned for pieces" far.txt && awk -F"\t" '"'"'NR>1 && $2=="S" && $4>=3 {n++} END {exit !(n>5000)}'"'"' dump_far.tsv'
# copies whose ends differ from the reference: anchored in the middle only, walked out to their ends through
# the locus's own k-mers, which stay private; the far copy of a stretch of one end stays outside (shared)
"$BIN" background -t 2 --partitions 64 --memory 64M -q --tmp-dir "$OUT/tmp" -o s_walk --genome ref=data/ends/ref.fa --genome gd=data/ends/gd.fa --coordinate-reference ref
"$BIN" recruit -t 2 -q --background s_walk -r ref --index-only --save-index idx_walk --locus data/ends/locus.gfa; rc=$?
"$BIN" recruit --dump-index idx_walk/locus.plx=dump_walk.tsv
walk() { awk -F'\t' '$1=="span" && $2=="gd" {print $8 "-" $9 ":" $10 "-" $11}' "$1"; }
check "a copy anchored in its middle is walked out to its ends ($(walk idx_walk/locus.index.tsv))" '[[ $rc -eq 0 && "$(walk idx_walk/locus.index.tsv)" =~ ^1[23][0-9]{4}-1[45][0-9]{4}:120000-160000$ ]]'
check "its ends' k-mers stay private, none demoted; the far copy's are shared" '[[ $(field idx_walk/locus.index.tsv span_demoted) == 0 ]] && python3 "$TOOLS/locus_outside.py" idx_walk/locus.index.tsv dump_walk.tsv 8 data/ends/ref.fa data/ends/gd.fa > brute_walk.txt && awk -F"\t" '"'"'NR>1 && $2=="S" {s++} NR>1 && $2=="P" {p++} END {exit !(s>1900 && s<2100 && p>55000)}'"'"' dump_walk.tsv'
# a repeat element between the anchored middle and a differing end is crossed; the same element in the flank is not
"$BIN" background -t 2 --partitions 64 --memory 64M -q --tmp-dir "$OUT/tmp" -o s_rep --genome ref=data/rep/ref.fa --genome gr=data/rep/gr.fa --coordinate-reference ref
"$BIN" recruit -t 2 -q --background s_rep -r ref --index-only --save-index idx_rep --locus data/rep/locus.gfa; rc=$?
"$BIN" recruit --dump-index idx_rep/locus.plx=dump_rep.tsv
rwalk() { awk -F'\t' '$1=="span" && $2=="gr" {print $8 "-" $9 ":" $10 "-" $11}' "$1"; }
check "a copy is walked across a repeat inside it, not into one past its end ($(rwalk idx_rep/locus.index.tsv))" '[[ $rc -eq 0 && "$(rwalk idx_rep/locus.index.tsv)" =~ :120000-160000$ ]]'
check "the end past the repeat stays private" '[[ $(field idx_rep/locus.index.tsv span_demoted) == 0 ]] && python3 "$TOOLS/locus_outside.py" idx_rep/locus.index.tsv dump_rep.tsv 8 data/rep/ref.fa data/rep/gr.fa > brute_rep.txt'
"$BIN" recruit -t 4 -q --background s_split -r ref --index-only --save-index idx_split4 --locus data/split/locus.gfa
"$BIN" recruit --dump-index idx_split4/locus.plx=dump_split4.tsv
check "the piece scan on 4 threads gives the same index" 'cmp -s dump_split.tsv dump_split4.tsv'
"$BIN" recruit -t 1 -q --background s_split_p -r ref --index-only --save-index idx_split_p1 --locus data/split/locus.gfa
"$BIN" recruit -t 4 -q --background s_split_p -r ref --index-only --save-index idx_split_p4 --locus data/split/locus.gfa
"$BIN" recruit --dump-index idx_split_p1/locus.plx=dump_split_p1.tsv
"$BIN" recruit --dump-index idx_split_p4/locus.plx=dump_split_p4.tsv
check "genomes scanned in parallel give the index of one at a time" 'cmp -s dump_split_p1.tsv dump_split_p4.tsv && cmp -s <(grep "^piece" idx_split_p1/locus.index.tsv) <(grep "^piece" idx_split_p4/locus.index.tsv)'

# a graph genome counts one haplotype per sample#haplotype of its paths
"$BIN" background -t 2 --partitions 64 --memory 64M -q --tmp-dir "$OUT/tmp" -o s_graph --genome ref=data/ref.fa --genome lg=data/locus.gfa --coordinate-reference ref
"$BIN" recruit -t 2 -q --background s_graph -r ref --index-only --save-index idx_graph --locus data/locus.gfa 2> graph.txt; rc=$?
check "a graph genome's PanSN paths count as haplotypes (1 FASTA + 3 paths)" '[[ $rc -eq 0 ]] && grep -q "haplotypes.4 in the background" idx_graph/locus.index.tsv'
# a haplotype named twice over the same range is refused (a span could be cut from the wrong one)
awk -F'\t' 'BEGIN {OFS = "\t"} {print} $1 == "P" && $2 ~ /^ga#/ {print}' data/locus.gfa > dup.gfa
"$BIN" background -t 2 --partitions 64 --memory 64M -q --tmp-dir "$OUT/tmp" -o s_dup --genome ref=data/ref.fa --genome lg=dup.gfa --coordinate-reference ref 2> dup.txt; rc=$?
check "a graph genome naming a haplotype twice is refused" '[[ $rc -ne 0 ]] && grep -q "is given twice" dup.txt'
# a damaged index is refused when loaded, not scored with
cp idx/locus.plx damaged.plx
python3 - damaged.plx <<'PY'
import sys
p = sys.argv[1]; b = bytearray(open(p, "rb").read()); b[len(b) // 2] ^= 0x10; open(p, "wb").write(bytes(b))
PY
"$BIN" recruit --locus damaged.plx -1 reads_1.fq -2 reads_2.fq -o damaged.tsv -q 2> damaged.txt; rc=$?
check "a damaged index is refused" '[[ $rc -ne 0 ]] && grep -q "corrupt" damaged.txt'

# a package: the compacted counts, the anchors joined and packed, the manifest; recruit builds the same index from it
"${BG[@]}" -o s --compact
"${BG[@]}" -o s --package pkg; rc=$?
check "a store packages: counts.sshash, the manifest, one packed anchor file" '[[ $rc -eq 0 && -f pkg/counts.sshash && -f pkg/manifest.json && -f pkg/seg-0/anchors-ref.bin && ! -e pkg/seg-0/parts-0.bin ]]'
packed=$(wc -c < pkg/seg-0/anchors-ref.bin | tr -d ' '); unpacked=$(cat s/seg-*/anchors-ref.bin | wc -c | tr -d ' ')
check "packed anchors are well under half the size ($packed of $unpacked bytes)" '[[ $((packed * 2)) -lt $unpacked ]]'
"$BIN" recruit -t 2 -q --background pkg -r ref --index-only --save-index idx_pkg --locus data/locus.gfa; rc=$?
"$BIN" recruit --dump-index idx_pkg/locus.plx=dump_pkg.tsv
check "an index built from the package equals the store's" '[[ $rc -eq 0 ]] && cmp -s dump.tsv dump_pkg.tsv'
"${BG[@]}" -o s --locate chr1:120000-160000 > loc_store.tsv; "${BG[@]}" -o pkg --locate chr1:120000-160000 > loc_pkg.tsv
check "--locate reads the packed anchors alike" '[[ -s loc_store.tsv ]] && cmp -s loc_store.tsv loc_pkg.tsv'
"${BG[@]}" -o pkg --add --genome gx=data/ga.fa 2> pkg_add.txt; rc=$?
check "a package cannot be changed" '[[ $rc -ne 0 ]] && grep -q "packaged background" pkg_add.txt'
# genome files moved, or recompressed: their new paths are checked against what was counted
mkdir -p moved && gzip -c data/ga.fa > moved/ga.fa.gz && cp data/gb.fa moved/gb.fa
printf 'ga\tmoved/ga.fa.gz\ngb\tmoved/gb.fa\n' > moved.tsv
"$BIN" background -o pkg --relocate moved.tsv 2> reloc.txt; rc=$?
check "relocated genomes, one recompressed: each checked by its sequence" '[[ $rc -eq 0 ]] && grep -q "genomes: 2; each checked against the sequence counted" reloc.txt'
"$BIN" recruit -t 2 -q --background pkg -r ref --index-only --save-index idx_moved --locus data/locus.gfa; rc=$?
"$BIN" recruit --dump-index idx_moved/locus.plx=dump_moved.tsv
check "the index from the relocated genomes is the same" '[[ $rc -eq 0 ]] && cmp -s dump.tsv dump_moved.tsv'
printf 'ga\tdata/gb.fa\n' > wrong.tsv
"${BG[@]}" -o pkg --relocate wrong.tsv 2> wrong.txt; rc=$?
check "a genome file other than the one counted is refused" '[[ $rc -ne 0 ]] && grep -q "is not the genome counted" wrong.txt'
# a file of the size counted but another sequence (one base changed) is refused too
python3 - data/gb.fa moved/gb_mut.fa <<'PY'
import sys
lines = open(sys.argv[1]).read().split("\n")
i = next(k for k, l in enumerate(lines) if l and not l.startswith(">"))
c = lines[i][0]; lines[i] = ("A" if c != "A" else "C") + lines[i][1:]
open(sys.argv[2], "w").write("\n".join(lines))
PY
printf 'gb\tmoved/gb_mut.fa\n' > same_size.tsv
"${BG[@]}" -o pkg --relocate same_size.tsv 2> same_size.txt; rc=$?
check "a file of the size counted but another sequence is refused" '[[ $rc -ne 0 && $(wc -c < moved/gb_mut.fa) == $(wc -c < data/gb.fa) ]] && grep -q "is not the genome counted" same_size.txt'
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
check "recruit takes each locus's score thresholds and reports them" 'grep -q "^locus.locus.score.short reads: r_in" recruited.summary.tsv && grep -q "^locus.locus.certain_threshold.[0-9.-]*: held-out recall" recruited.summary.tsv && grep -q "^locus.locus.likely_threshold" recruited.summary.tsv'
check "the table appends score and tier" '[[ "$(head -1 recruited.tsv | cut -f25-)" == "score	tier" ]] && awk -F"\t" '"'"'NR > 1 && ($25 !~ /^-?[0-9]+\.[0-9][0-9]$/ || ($3 == "certain") != ($26 == "certain")) {bad = 1} END {exit bad}'"'"' recruited.tsv'
"$BIN" recruit --locus idx/locus.plx -1 reads_1.fq -2 reads_2.fq -o rules1.tsv -t 2 -q --rules v1
check "the first rules: the calibrated cut, reported" 'grep -q "^locus.locus.site_cut.2.calibrated (short reads)" rules1.summary.tsv'
"$BIN" recruit --locus idx/locus.plx -1 reads_1.fq -2 reads_2.fq -o fixed.tsv -t 2 -q --min-sites 3
check "--min-sites gives every locus that cut" 'grep -q "^locus.locus.site_cut.3.given" fixed.summary.tsv'
"$BIN" recruit --locus idx/locus.plx -1 reads_1.fq -2 reads_2.fq -o ontcut.tsv -t 2 -q --preset ont --rules v1
check "the first rules, ont reads: the calibrated cut starts at 3" 'grep -q -E "^locus.locus.site_cut.[3-8] calibrated \(ont reads\).*ont reads start at 3" ontcut.summary.tsv'
check "summary rows keep four fields" '[[ $(awk -F"\t" '"'"'NF != 4'"'"' recruited.summary.tsv | wc -l | tr -d " ") == 0 ]]'
"$BIN" recruit --locus idx/locus.plx -1 reads_1.fq -2 reads_2.fq -o both.tsv -q --min-sites 3 --false-certain-rate 1e-3 2> both.txt; rc=$?
check "--min-sites and --false-certain-rate together are refused" '[[ $rc -ne 0 ]] && grep -q "not both" both.txt'
"$BIN" recruit -t 2 -q --background s -r ref --index-only --save-index idx_nocal --locus data/locus.gfa --no-calibration
"$BIN" recruit --locus idx_nocal/locus.plx -1 reads_1.fq -2 reads_2.fq -o nocal.tsv -t 2 -q
check "an index without calibration: no calibration rows, the cut 2, the first rules' calls, no score" '[[ -z "$(cal idx_nocal/locus.index.tsv)" ]] && grep -q "site_cut.2 (the index carries no calibration)" nocal.summary.tsv && cmp -s <(cut -f1-24 rules1.tsv) <(cut -f1-24 nocal.tsv) && [[ "$(tail -n +2 nocal.tsv | cut -f25-26 | sort -u)" == ".	." ]]'
"$BIN" recruit --locus idx_nocal/locus.plx -1 reads_1.fq -2 reads_2.fq -o nocal_ont.tsv -t 2 -q --preset ont
check "an index without calibration, ont reads: the cut 3" 'grep -q "site_cut.3 (the index carries no calibration)" nocal_ont.summary.tsv'
mkdir -p twin/a twin/b && cp data/locus.gfa twin/a/locus.gfa && cp data/locus.gfa twin/b/locus.gfa
"$BIN" recruit -t 2 -q --background s -r ref --index-only --save-index idx_twin --locus twin/a/locus.gfa --locus twin/b/locus.gfa 2> twin.txt; rc=$?
check "two graphs of one name are refused" '[[ $rc -ne 0 ]] && grep -q "rename one graph" twin.txt'

"${RC[@]}" --locus data/locus.gfa 2> log2.txt
check "a current index is reused" 'grep -q "; current, reused" log2.txt'
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
  # the piece fixtures read by interval (around the copies and at contig ends only): the same pieces and index
  same=0; cases=0
  for spec in "split split/locus.gfa split/ref:ref split/gs:gs" "split_p split/locus.gfa split/ref:ref split/gs:gs split/gp:gp" \
              "split_c split/locus.gfa split/ref:ref split_c/gc:gc" "split_h split/locus.gfa split/ref:ref split_h/gh:gh" \
              "tail split_t/locus.gfa split_t/ref:ref split_t/gl:gl" "orphan split_o/locus.gfa split_o/ref:ref split_o/go:go" \
              "ends split_e/locus.gfa split_e/ref:ref split_e/ge:ge split_e/gt:gt" "far split_f/locus.gfa split_f/ref:ref split_f/gf:gf"; do
    read -r c gfa rest <<< "$spec"; args=()
    for x in $rest; do f=${x%%:*}; n=${x##*:}; mkdir -p zp/$c
      bgzip -c data/$f.fa > zp/$c/$n.fa.gz && samtools faidx zp/$c/$n.fa.gz; args+=(--genome $n=zp/$c/$n.fa.gz); done
    "$BIN" background -t 2 --partitions 64 --memory 64M -q --tmp-dir "$OUT/tmp" -o s_z$c "${args[@]}" --coordinate-reference ref
    "$BIN" recruit -t 2 -q --background s_z$c -r ref --index-only --save-index idx_z$c --locus data/$gfa
    "$BIN" recruit --dump-index idx_z$c/locus.plx=dump_z$c.tsv
    cases=$((cases + 1))
    cmp -s dump_$c.tsv dump_z$c.tsv && cmp -s <(grep "^piece" idx_$c/locus.index.tsv) <(grep "^piece" idx_z$c/locus.index.tsv) && same=$((same + 1))
  done
  check "genomes read by interval give the pieces and index of a full read ($same of $cases)" '[[ $same == $cases ]]'
else
  echo "  skip indexed and bgzip genomes (samtools or bgzip not installed)"
fi

# the score. A locus whose homologs are held by other genomes only: each calibration genome held out of the background
# in turn, the locus k-mers only it shares read as private, so reads from its homolog carry chance private sites
"$BIN" background -t 2 --partitions 64 --memory 64M -q --tmp-dir "$OUT/tmp" -o s_logo --genome ref=data/logo/ref.fa --genome gu=data/logo/gu.fa --genome gd=data/logo/gd.fa --coordinate-reference ref
LG=("$BIN" recruit -q --background s_logo -r ref --index-only --locus data/logo/locus.gfa)
"${LG[@]}" -t 2 --save-index idx_logo --decoy-genomes ref,gd,gu; rc=$?
"$BIN" recruit --dump-index idx_logo/locus.plx=dump_logo.tsv
cg() { awk -F'\t' -v g="$2" -v c="$3" '$1=="calibration_genome" && $2==g {print $c}' "$1"; }
shared_logo=$(awk -F'\t' 'NR > 1 && $2 == "S"' dump_logo.tsv | wc -l | tr -d ' ')
check "the decoy genomes other than the reference are held out in calibration, the one holding the most homolog k-mers first ($(field idx_logo/locus.index.tsv calibration_genomes))" '[[ $rc -eq 0 && "$(field idx_logo/locus.index.tsv calibration_genomes)" == "gu,gd (up to 8 of the decoy genomes, by the homolog k-mers each adds)" ]]'
check "each relabels as private the locus k-mers only it shares (gu $(cg idx_logo/locus.index.tsv gu 6) + gd $(cg idx_logo/locus.index.tsv gd 6) of $shared_logo shared)" '[[ $(( $(cg idx_logo/locus.index.tsv gu 6) + $(cg idx_logo/locus.index.tsv gd 6) )) == $shared_logo && $(cg idx_logo/locus.index.tsv gu 6) -gt 1000 && $(cg idx_logo/locus.index.tsv gd 6) -gt 1000 ]]'
genome_private() { awk -F'\t' -v g="$2" '$1=="score_genome" && $2=="short" && $3==g {print $6 " of " $5}' "$1"; }
check "its homolog's reads are scored with chance private sites (gu: $(genome_private idx_logo/locus.index.tsv gu) screened)" 'awk -F"\t" '"'"'$1=="score_genome" && $2=="short" && $3=="gu" && $6 > 0 && $5 >= $6 {ok = 1} END {exit !ok}'"'"' idx_logo/locus.index.tsv'
"${LG[@]}" -t 2 --save-index idx_logo1 --decoy-genomes ref,gd,gu --calibration-max 1
check "--calibration-max 1: the genome adding the most homolog k-mers is taken, not the first listed ($(field idx_logo1/locus.index.tsv calibration_genomes))" '[[ "$(field idx_logo1/locus.index.tsv calibration_genomes)" == "gu (up to 1 of the decoy genomes, by the homolog k-mers each adds)" ]]'
"${LG[@]}" -t 2 --save-index idx_logo_c --decoy-genomes ref,gd,gu --calibration-genomes gd,gu --calibration-max 1
check "--calibration-genomes names the candidates: of them the locus takes the one adding the most ($(field idx_logo_c/locus.index.tsv calibration_genomes))" '[[ "$(field idx_logo_c/locus.index.tsv calibration_genomes)" == "gu (up to 1 of the 2 given, by the homolog k-mers each adds)" ]]'
"${LG[@]}" -t 2 --save-index idx_logo_m --decoy-genomes ref,gd,gu --calibration-max 0 2> max0.txt; rc=$?
check "--calibration-max 0 is refused" '[[ $rc -ne 0 ]] && grep -q "at least 1" max0.txt'
# candidates that are not decoy genomes are not scanned whole: their windows are where the anchors place the
# coordinate reference's (the paralog of the first fixture, which ga and gb hold too)
"$BIN" recruit -t 2 --background s -r ref --index-only --save-index idx_cand --locus data/locus.gfa --calibration-genomes ga,gb 2> cand.txt; rc=$?
cwin=$(awk -F'\t' '$1=="calibration_window" && $2=="ga" {print $3 ":" $4 "-" $5}' idx_cand/locus.index.tsv | head -1)
check "candidates that are not decoy genomes: their windows placed by the anchors, the genomes not scanned whole (ga: $cwin)" '[[ $rc -eq 0 ]] && awk -F"\t" '"'"'$1=="calibration_window" && $2=="ga" && $3=="chr1" && $4 > 225000 && $4 < 235000 && $5 > 255000 && $5 < 265000 {ok = 1} END {exit !ok}'"'"' idx_cand/locus.index.tsv && ! grep -q "genome scanned: ga" cand.txt && grep -q "calibration genomes: scanned: ga" cand.txt && [[ $(field idx_cand/locus.index.tsv decoy_windows | cut -d" " -f1) == $(field idx/locus.index.tsv decoy_windows | cut -d" " -f1) ]]'
check "a candidate adding no homolog k-mer to those taken is not held out ($(field idx_cand/locus.index.tsv calibration_genomes))" '[[ "$(field idx_cand/locus.index.tsv calibration_genomes)" == "ga (up to 8 of the 2 given, by the homolog k-mers each adds)" ]]'
"$BIN" recruit -t 2 -q --background s -r ref --index-only --save-index idx_cand_ref --locus data/locus.gfa --calibration-genomes ref,ga 2> candref.txt; rc=$?
check "the coordinate reference as a calibration genome is refused" '[[ $rc -ne 0 ]] && grep -q "coordinate reference" candref.txt'
"${LG[@]}" -t 2 --save-index idx_logo_n
check "no genome to hold out besides the reference: said so" '[[ "$(field idx_logo_n/locus.index.tsv calibration_genomes)" == none\ \(* ]]'
"${LG[@]}" -t 1 --save-index idx_logo_t1 --decoy-genomes ref,gd,gu
check "calibration with genomes held out on 1 thread equals 2 threads" 'grep -q "^score_genome.short.gu" idx_logo_t1/locus.index.tsv && cmp -s idx_logo/locus.index.tsv idx_logo_t1/locus.index.tsv && cmp -s idx_logo/locus.plx idx_logo_t1/locus.plx'
# a long read of the stretch gd shares, with two private sites spliced in 150 bases apart, and a 350-base read of it around
# them: the same sites, but the long read holds hundreds of k-mers that would carry more sites if it were the locus
python3 - data/logo/ref.fa > reads_len.fq <<'PY'
import sys
seq = "".join(l.strip() for l in open(sys.argv[1]) if not l.startswith(">"))
read = seq[146000:149000] + seq[125000:125040] + seq[149040:149150] + seq[126000:126040] + seq[149190:152000]
short = read[2950:3300]
for name, r in (("long", read), ("short", short)):
    print(f"@{name}\n{r}\n+\n{'I' * len(r)}")
PY
"$BIN" recruit --locus idx_logo/locus.plx -i reads_len.fq -o len.tsv --preset hifi -t 2 -q
lenrow() { awk -F'\t' -v f="$1" '$1==f {print "sites " $8 ", score " $25}' len.tsv; }
check "a long stretch with few sites scores below a short dense one (long: $(lenrow long); short: $(lenrow short))" 'awk -F"\t" '"'"'$1=="long" {l = $25; lp = $8} $1=="short" {s = $25; sp = $8} END {exit !(lp == 2 && sp == 2 && l + 0 < s + 0)}'"'"' len.tsv'
# one private site is never certain: with every threshold as low as it goes, a read holding one spliced-in site is
# likely, the same read with a second site certain
python3 - data/logo/ref.fa > reads_one.fq <<'PY'
import sys
seq = "".join(l.strip() for l in open(sys.argv[1]) if not l.startswith(">"))
one = seq[149000:149050] + seq[125000:125040] + seq[149090:149200]
two = seq[149000:149050] + seq[125000:125040] + seq[149090:149150] + seq[126000:126040] + seq[149190:149260]
for name, r in (("one", one), ("two", two)):
    print(f"@{name}\n{r}\n+\n{'I' * len(r)}")
PY
"$BIN" recruit --locus idx_logo/locus.plx -i reads_one.fq -o one.tsv --preset illumina --certain-rate 0.9 --likely-rate 0.9 -t 2 -q
onerow() { awk -F'\t' -v f="$1" '$1==f {print $3 "/" $4 "/" $8 " sites/tier " $26}' one.tsv; }
check "one private site is never certain, two can be (one: $(onerow one); two: $(onerow two))" '[[ "$(onerow one)" == "likely/single_site/1 sites/tier likely" && "$(onerow two)" == certain/*/2\ sites/tier\ certain ]]'
# thresholds follow the asked rate; the table is the same on any thread count. Pairs from the locus and from each homolog
python3 "$TOOLS/locus_reads.py" data/logo/ref.fa chr1 120000 160000 2000 loc 7 reads_lg_locus
python3 "$TOOLS/locus_reads.py" data/logo/gu.fa u 20000 32500 1000 hu 8 reads_lg_u
python3 "$TOOLS/locus_reads.py" data/logo/gd.fa d 20000 30000 1000 hd 9 reads_lg_d
cat reads_lg_locus_1.fq reads_lg_u_1.fq reads_lg_d_1.fq > reads_lg_1.fq
cat reads_lg_locus_2.fq reads_lg_u_2.fq reads_lg_d_2.fq > reads_lg_2.fq
prev=""; mono=1; counts=""
for rate in 0.1 0.01 0.001 0.0001; do
  "$BIN" recruit --locus idx_logo/locus.plx -1 reads_lg_1.fq -2 reads_lg_2.fq -o rate_$rate.tsv --certain-rate $rate --likely-rate 0.1 -t 2 -q
  th=$(awk -F'\t' '$3=="certain_threshold" {split($4, a, ":"); print a[1]}' rate_$rate.summary.tsv 2> /dev/null); [[ $th == any ]] && th=-1e9
  [[ -z $th ]] && mono=0
  n=$(awk -F'\t' 'NR > 1 && $3=="certain"' rate_$rate.tsv | wc -l | tr -d ' ')
  if [[ -n $prev ]] && ! awk -v a="$prev" -v b="$th" -v m="$prev_n" -v n="$n" 'BEGIN {exit !(b + 0 >= a + 0 && n <= m)}'; then mono=0; fi
  prev=$th; prev_n=$n; counts="$counts $rate:$th:$n"
done
check "a lower --certain-rate never lowers a threshold nor adds certain fragments ($counts)" '[[ $mono == 1 ]]'
"$BIN" recruit --locus idx_logo/locus.plx -1 reads_lg_1.fq -2 reads_lg_2.fq -o t1.tsv -t 1 -q --likely-rate 0.1
"$BIN" recruit --locus idx_logo/locus.plx -1 reads_lg_1.fq -2 reads_lg_2.fq -o t2.tsv -t 2 -q --likely-rate 0.1 --reads-out t2_reads
check "scores and tiers on 1 thread equal 2 threads" 'cmp -s t1.tsv t2.tsv && [[ $(awk -F"\t" '"'"'NR > 1 && $25 != "."'"'"' t1.tsv | wc -l) -gt 0 ]]'
nlikely=$(awk -F'\t' 'NR > 1 && $3=="likely"' t2.tsv | wc -l | tr -d ' ')
check "likely fragments go to their own files ($nlikely)" '[[ $nlikely -gt 0 && $(gzip -dc t2_reads/locus.likely_1.fq.gz | awk "END {print NR / 4}") == $nlikely ]]'
"$BIN" recruit --locus idx_logo/locus.plx -1 reads_lg_1.fq -2 reads_lg_2.fq -o bad.tsv --certain-rate 0.01 --likely-rate 0.001 -q 2> bad.txt; rc=$?
check "--likely-rate below --certain-rate is refused" '[[ $rc -ne 0 ]] && grep -q "likely-rate must be at least" bad.txt'
# a genome holding a copy of part of the locus that no anchors place (its pieces in reverse order): held out, its
# reads there are the locus's, not off-target, and certain reads stay certain
"$BIN" background -t 2 --partitions 64 --memory 64M -q --tmp-dir "$OUT/tmp" -o s_copy --genome ref=data/logo/ref.fa --genome gc=data/logo/gc.fa --coordinate-reference ref
"$BIN" recruit -q -t 2 --background s_copy -r ref --index-only --locus data/logo/locus.gfa --save-index idx_copy --decoy-genomes ref,gc
copy_bases=$(cg idx_copy/locus.index.tsv gc 8); copy_frag=$(awk -F'\t' '$1=="score_genome" && $2=="short" && $3=="gc" {print $7}' idx_copy/locus.index.tsv)
copy_recall=$(awk -F'\t' '$1=="score_tier" && $2=="short" && $3=="0.0001" {print $6}' idx_copy/locus.index.tsv)
check "a held-out genome's stretch that reads as a copy of the locus is not off-target (gc: ${copy_bases:-none} bases, ${copy_frag:-no} fragments left out; certain recall ${copy_recall:-none})" '[[ ${copy_bases:-0} -ge 5000 && ${copy_frag:-0} -gt 0 ]] && awk -v r="${copy_recall:-0}" '"'"'BEGIN {exit !(r >= 0.8)}'"'"' && [[ $(awk -F"\t" '"'"'$1=="span" && $2=="gc"'"'"' idx_copy/locus.index.tsv | wc -l | tr -d " ") == 1 ]]'
# the off-target rates are the median genome's: gd, which holds much of the locus in pieces, does not set them
"$BIN" background -t 2 --partitions 64 --memory 64M -q --tmp-dir "$OUT/tmp" -o s_logo3 --genome ref=data/logo/ref.fa --genome gu=data/logo/gu.fa --genome gd=data/logo/gd.fa --genome ge=data/logo/ge.fa --coordinate-reference ref
"$BIN" recruit -q -t 2 --background s_logo3 -r ref --index-only --locus data/logo/locus.gfa --save-index idx_logo3 --decoy-genomes ref,gu,gd,ge
median_ok() {   # the short class's r_out is the middle one of its three genomes' own
  awk -F'\t' '$1=="score_rates" && $2=="short" {r = $4} $1=="score_genome" && $2=="short" {v[++n] = $8} END {
    if (n != 3) exit 1
    for (i = 1; i <= 3; i++) for (j = i + 1; j <= 3; j++) if (v[j] + 0 < v[i] + 0) {t = v[i]; v[i] = v[j]; v[j] = t}
    exit !(v[2] == r && v[3] > v[2] && v[2] > v[1]) }' "$1"
}
rates=$(awk -F'\t' '$1=="score_genome" && $2=="short" {print $3 " " $8}' idx_logo3/locus.index.tsv | tr '\n' ',')
check "the off-target site rate is the median of the held-out genomes' own (r_out $(awk -F'\t' '$1=="score_rates" && $2=="short" {print $4}' idx_logo3/locus.index.tsv); $rates)" 'median_ok idx_logo3/locus.index.tsv'
# a rate the off-target fragments are too few to show: the threshold is read off their fitted tail, never below the
# one they give, and the run says so; a rate they can show is the one observed
tier() { awk -F'\t' -v r="$2" -v c="$3" '$1=="score_tier" && $2=="short" && $3==r {print $c}' "$1"; }
check "a rate too low to show: the threshold extrapolated from the off-target tail ($(tier idx/locus.index.tsv 0.0001 4), $(tier idx/locus.index.tsv 0.0001 12)), said in the summary; a rate shown: observed ($(tier idx/locus.index.tsv 0.01 4))" '[[ "$(tier idx/locus.index.tsv 0.0001 11)" == extrapolated && "$(tier idx/locus.index.tsv 0.01 11)" == observed && "$(tier idx/locus.index.tsv 0.01 12)" == "." ]] && grep -q "^locus.locus.certain_threshold.*too few to show the rate: the threshold is extrapolated from the tail of the off-target scores (scale .*; the fragments seen give" recruited.summary.tsv && ! grep -q "^locus.locus.likely_threshold.*extrapolated" recruited.summary.tsv'
# an index of an earlier format is refused, saying so
cp idx_logo/locus.plx old.plx
python3 - old.plx <<'PY'
import struct, sys
p = sys.argv[1]; b = bytearray(open(p, "rb").read()); b[8:16] = struct.pack("<Q", 15); open(p, "wb").write(bytes(b))
PY
"$BIN" recruit --locus old.plx -1 reads_1.fq -2 reads_2.fq -o old.tsv -q 2> old.txt; rc=$?
check "an index of format 15 is refused: rebuild it" '[[ $rc -ne 0 ]] && grep -q "of format 15; this panvar reads format 17" old.txt && grep -q "rebuild" old.txt'
# the first rules reproduce a reference build's tables (RECRUIT_V1_BIN: a panvar of the first rules)
if [[ -n "${RECRUIT_V1_BIN:-}" ]]; then
  same=0; cases=0
  "$RECRUIT_V1_BIN" recruit -t 2 -q --background s -r ref --index-only --save-index idx_v1ref --locus data/locus.gfa 2> /dev/null
  "$BIN" recruit -t 2 -q --background s -r ref --index-only --save-index idx_v1new --locus data/locus.gfa
  for x in "reads_1.fq -2 reads_2.fq" "reads_edge_1.fq -2 reads_edge_2.fq" "reads_long_1.fq" "reads_span_1.fq" "reads_flank.fq:i"; do
    args=(-1 ${x%%:*}); [[ $x == *:i ]] && args=(-i ${x%%:*})
    "$RECRUIT_V1_BIN" recruit --locus idx_v1ref/locus.plx ${args[@]} -o v1ref_$cases.tsv -t 2 -q
    "$BIN" recruit --locus idx_v1new/locus.plx ${args[@]} -o v1new_$cases.tsv -t 2 -q --rules v1
    cmp -s <(cut -f1-24 v1new_$cases.tsv) v1ref_$cases.tsv && cmp -s <(grep -v "^input" v1new_$cases.summary.tsv) <(grep -v "^input" v1ref_$cases.summary.tsv) && same=$((same + 1))
    cases=$((cases + 1))
  done
  check "--rules v1 gives the reference build's tables and summaries ($same of $cases inputs)" '[[ $same == $cases ]]'
fi

# a full store relocated: its segments' records follow, so later steps read the new files
cp -r s s_rel && mkdir -p moved_s && cp data/ga.fa moved_s/ga.fa
printf 'ga\tmoved_s/ga.fa\n' > moved_s.tsv
"${BG[@]}" -o s_rel --relocate moved_s.tsv; rc=$?
check "a store's relocated genomes are recorded in its segments too" '[[ $rc -eq 0 ]] && grep -q "moved_s/ga.fa" s_rel/seg-*/segment.tsv && ! grep -q "data/ga.fa" s_rel/seg-*/segment.tsv'

echo
echo "recruit_index: $pass passed, $fail failed"
[[ $fail -eq 0 ]]
