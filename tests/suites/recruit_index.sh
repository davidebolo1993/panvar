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
check "relocated genomes: one by size, one recompressed by its sequence checksum" '[[ $rc -eq 0 ]] && grep -q "1 the size counted, 1 another size but the sequence counted" reloc.txt'
"$BIN" recruit -t 2 -q --background pkg -r ref --index-only --save-index idx_moved --locus data/locus.gfa; rc=$?
"$BIN" recruit --dump-index idx_moved/locus.plx=dump_moved.tsv
check "the index from the relocated genomes is the same" '[[ $rc -eq 0 ]] && cmp -s dump.tsv dump_moved.tsv'
printf 'ga\tdata/gb.fa\n' > wrong.tsv
"${BG[@]}" -o pkg --relocate wrong.tsv 2> wrong.txt; rc=$?
check "a genome file other than the one counted is refused" '[[ $rc -ne 0 ]] && grep -q "is not the genome counted" wrong.txt'
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
"$BIN" recruit --locus idx/locus.plx -1 reads_1.fq -2 reads_2.fq -o ontcut.tsv -t 2 -q --preset ont
check "ont reads: the calibrated cut starts at 3" 'grep -q -E "^locus.locus.site_cut.[3-8] calibrated \(ont reads\).*ont reads start at 3" ontcut.summary.tsv'
check "summary rows keep four fields" '[[ $(awk -F"\t" '"'"'NF != 4'"'"' recruited.summary.tsv | wc -l | tr -d " ") == 0 ]]'
"$BIN" recruit --locus idx/locus.plx -1 reads_1.fq -2 reads_2.fq -o both.tsv -q --min-sites 3 --false-certain-rate 1e-3 2> both.txt; rc=$?
check "--min-sites and --false-certain-rate together are refused" '[[ $rc -ne 0 ]] && grep -q "not both" both.txt'
"$BIN" recruit -t 2 -q --background s -r ref --index-only --save-index idx_nocal --locus data/locus.gfa --no-calibration
"$BIN" recruit --locus idx_nocal/locus.plx -1 reads_1.fq -2 reads_2.fq -o nocal.tsv -t 2 -q
check "an index without calibration: no calibration rows, the cut 2, the same table" '[[ -z "$(cal idx_nocal/locus.index.tsv)" ]] && grep -q "site_cut.2 (the index carries no calibration)" nocal.summary.tsv && cmp -s recruited.tsv nocal.tsv'
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
else
  echo "  skip indexed and bgzip genomes (samtools or bgzip not installed)"
fi

echo
echo "recruit_index: $pass passed, $fail failed"
[[ $fail -eq 0 ]]
