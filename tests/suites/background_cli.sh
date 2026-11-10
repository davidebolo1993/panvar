#!/usr/bin/env bash
# panvar background end to end: stores built, added to, killed at every commit step, added to by two jobs
# at once, built as standalone segments and attached, merged in ranges, compacted, anchored after the fact,
# and run out of disk. After each step the store's dump must equal an independent count of exactly the
# genomes committed (tests/tools/kmer_count.py), or the run must fail without changing it.
#   background_cli.sh <panvar_bin> <tests/tools dir> <out_dir>
set -uo pipefail

if [[ $# -lt 3 ]]; then
  echo "Usage: $0 <panvar_bin> <tools_dir> <out_dir>" >&2
  exit 2
fi
BIN="$1"; TOOLS="$2"; OUT="$3"
rm -rf "$OUT"; mkdir -p "$OUT"; cd "$OUT" || exit 2
python3 "$TOOLS/background_fixture.py" data
BG=("$BIN" background -t 2 --partitions 64 --memory 64M -q --tmp-dir "$OUT/tmp")
mkdir -p tmp

pass=0; fail=0
ok() { echo "  ok   $1"; pass=$((pass + 1)); }
bad() { echo "  FAIL $1"; fail=$((fail + 1)); }
check() { if eval "$2"; then ok "$1"; else bad "$1"; fi; }

# the expected counts of a set of genomes, and a store's dump, both sorted
truth() { local out="$1"; shift; python3 "$TOOLS/kmer_count.py" "$@" > "$out"; }
dump() { "${BG[@]}" -o "$1" --dump "$2.raw" && sort "$2.raw" > "$2"; }
same() { dump "$1" "$1.dump" && cmp -s "$1.dump" "$2"; }
listed() { grep -o '"dir": "seg-[0-9]*"' "$1/manifest.json" | wc -l | tr -d ' '; }
on_disk() { ls -d "$1"/seg-* 2>/dev/null | wc -l | tr -d ' '; }

truth t_ref_g1 data/ref.fa data/g1.fa
truth t_12 data/ref.fa data/g1.fa data/g2.fa
truth t_123 data/ref.fa data/g1.fa data/g2.fa data/g3.fa
truth t_all data/ref.fa data/g1.fa data/g2.fa data/g3.fa data/g4.fa
# more genomes, each a sequence no other file has
two() { { sed "s/^>/>a_/" "$1"; sed "s/^>/>b_/" "$2"; } > "$3"; }   # two genomes as one, contig names kept apart
two data/g4.fa data/g1.fa data/g4b.fa
two data/g1.fa data/g3.fa data/g6.fa
two data/g2.fa data/g3.fa data/g7.fa
two data/g2.fa data/g4.fa data/g5.fa
cp data/g2.fa data/g2copy.fa

echo "== failed and interrupted runs leave the store as it was"
"${BG[@]}" -o s --genome ref=data/ref.fa --genome g1=data/g1.fa --coordinate-reference reff 2> err.txt; rc=$?
check "a mistyped coordinate reference fails before counting" '[[ $rc -ne 0 ]] && grep -q "not among the genomes" err.txt && [[ $(listed s) == 0 ]]'
"${BG[@]}" -o s --genome ref=data/ref.fa --genome g1=data/g1.fa --coordinate-reference ref
check "run again correctly, every count once" 'same s t_ref_g1'
PANVAR_TEST_FAIL_AT=count-before-commit "${BG[@]}" -o s --add --genome g2=data/g2.fa; rc=$?
check "killed after counting, before the commit: exit 137, store unchanged" '[[ $rc -eq 137 ]] && same s t_ref_g1'
"${BG[@]}" -o s --clean
check "--clean removes the dead job's staging" '[[ -z "$(ls -A s/.staging 2>/dev/null)" ]]'
"${BG[@]}" -o s --add --genome g2=data/g2.fa
check "the same --add again commits it" 'same s t_12'
PANVAR_TEST_FAIL_AT=commit-before-manifest "${BG[@]}" -o s --add --genome g3=data/g3.fa; rc=$?
check "killed between moving the segment in and writing the manifest: unlisted, ignored" '[[ $rc -eq 137 ]] && same s t_12 && [[ $(on_disk s) -gt $(listed s) ]]'
"${BG[@]}" -o s --add --genome g3=data/g3.fa
check "the same --add again commits it under a free name" 'same s t_123'
"${BG[@]}" -o s --clean
check "--clean removes the unlisted segment" '[[ $(on_disk s) == $(listed s) ]]'
"${BG[@]}" -o s --add --genome g1b=data/g1.fa --genome g1=data/g4.fa 2> err.txt; rc=$?
check "a genome name already in the store is refused" '[[ $rc -ne 0 ]] && grep -q "already in the store" err.txt && same s t_123'

echo "== two jobs adding at once both commit"
cp -r s s_conc
"${BG[@]}" -o s_conc --add --genome g4=data/g4.fa & p1=$!
"${BG[@]}" -o s_conc --add --genome g4b=data/g4b.fa & p2=$!
wait $p1; r1=$?; wait $p2; r2=$?
truth t_conc data/ref.fa data/g1.fa data/g2.fa data/g3.fa data/g4.fa data/g4b.fa
check "both exit 0 and both genomes are counted" '[[ $r1 -eq 0 && $r2 -eq 0 ]] && same s_conc t_conc'

echo "== standalone segments, attached"
"${BG[@]}" -o s --segment-only seg_g4 --genome g4=data/g4.fa
check "a segment built for the store is not part of it yet" 'same s t_123 && [[ -f seg_g4/segment.tsv && -f seg_g4/anchors-ref.bin ]]'
"${BG[@]}" -o s --attach seg_g4
check "attached, its genome counts" 'same s t_all && [[ ! -d seg_g4 ]]'

echo "== segment folders and failed attaches leave everything where it was"
cp -r s s_att
mkdir -p precious && echo keep > precious/results.txt
"${BG[@]}" -o s_att --segment-only precious --genome g6=data/g6.fa 2> err.txt; rc=$?
check "--segment-only refuses a folder that exists, and leaves it alone" '[[ $rc -ne 0 && -f precious/results.txt ]]'
PANVAR_TEST_WAIT_AT=segment-built:$OUT/go1 "${BG[@]}" -o s_att --segment-only seg_x --genome g6=data/g6.fa & p1=$!
for _ in $(seq 1 200); do [[ -d seg_x ]] && break; sleep 0.1; done
"${BG[@]}" -o s_att --segment-only seg_x --genome g7=data/g7.fa 2> err.txt; rc=$?
touch go1; wait $p1; r1=$?
check "a second job given the same --segment-only folder is refused; the first completes" '[[ $rc -ne 0 && $r1 -eq 0 ]] && grep -q "g6" seg_x/segment.tsv && ! grep -q "g7" seg_x/segment.tsv'
"${BG[@]}" -o s_att --segment-only seg_y --genome g2copy=data/g2copy.fa
"${BG[@]}" -o s_att --attach seg_x seg_y 2> err.txt; rc=$?
check "attaching a segment whose genome is already counted fails; the other source stays" '[[ $rc -ne 0 && -f seg_x/segment.tsv && -f seg_y/segment.tsv ]] && grep -q "same sequence" err.txt && same s_att t_all'
PANVAR_TEST_FAIL_AT=attach-before-manifest "${BG[@]}" -o s_att --attach seg_x; rc=$?
check "an attach killed before its commit: store unchanged, source still there" '[[ $rc -eq 137 && -f seg_x/segment.tsv ]] && same s_att t_all'
"${BG[@]}" -o s_att --clean && "${BG[@]}" -o s_att --attach seg_x
truth t_all6 data/ref.fa data/g1.fa data/g2.fa data/g3.fa data/g4.fa data/g6.fa
check "attached after --clean, its genome counts once" 'same s_att t_all6 && [[ $(on_disk s_att) == $(listed s_att) && ! -d seg_x ]]'

echo "== merging in ranges"
PANVAR_TEST_FAIL_AT=merge-range "${BG[@]}" -o s --merge --range 0/3; rc=$?
check "a range job killed after writing its part: store unchanged" '[[ $rc -eq 137 ]] && same s t_all && [[ $(listed s) == 4 ]]'
"${BG[@]}" -o s --merge --range 1/3 && "${BG[@]}" -o s --merge --range 0/3
check "two of three ranges done: not committed yet" 'same s t_all && [[ $(listed s) == 4 ]]'
"${BG[@]}" -o s --merge --range 2/3
check "the last range commits one segment with the same counts" 'same s t_all && [[ $(listed s) == 1 ]]'
check "the merged-away segments stay for jobs still reading them" '[[ $(on_disk s) == 5 ]]'
"${BG[@]}" -o s --clean
check "--clean removes them" '[[ $(on_disk s) == 1 ]] && same s t_all'
cp -r s_att s_mb
PANVAR_TEST_FAIL_AT=merge-before-manifest "${BG[@]}" -o s_mb --merge; rc=$?
check "a merge killed just before its commit: store unchanged" '[[ $rc -eq 137 ]] && same s_mb t_all6 && [[ $(listed s_mb) == 5 ]]'
"${BG[@]}" -o s_mb --merge --range 0/2 && "${BG[@]}" -o s_mb --add --genome g7=data/g7.fa && "${BG[@]}" -o s_mb --merge --range 1/2
truth t_all67 data/ref.fa data/g1.fa data/g2.fa data/g3.fa data/g4.fa data/g6.fa data/g7.fa
check "an --add between range jobs: the merge commits the segments it planned, the new one stays" 'same s_mb t_all67 && [[ $(listed s_mb) == 2 ]]'
cp -r s s_m2
"${BG[@]}" -o s_m2 --add --genome g4b=data/g4b.fa
PANVAR_TEST_FAIL_AT=merge-after-commit "${BG[@]}" -o s_m2 --merge; rc=$?
truth t_m2 data/ref.fa data/g1.fa data/g2.fa data/g3.fa data/g4.fa data/g4b.fa
check "killed after the merge commit: counts right" '[[ $rc -eq 137 ]] && same s_m2 t_m2 && [[ $(listed s_m2) == 1 && $(on_disk s_m2) == 3 ]]'
"${BG[@]}" -o s_m2 --clean
check "--clean removes the merged-away segments" '[[ $(on_disk s_m2) == 1 ]]'

echo "== compaction"
"${BG[@]}" -o s --compact
"${BG[@]}" -o s --verify s.sshash; rc=$?
check "the compacted background holds exactly the store's k-mers and counts" '[[ $rc -eq 0 ]]'
"${BG[@]}" -o s --verify s_m2 > /dev/null 2>&1; rc=$?
check "--verify fails against a store with more k-mers" '[[ $rc -ne 0 ]]'
"${BG[@]}" -o s --compact --local-strings
"${BG[@]}" -o s --verify s.sshash; rc=$?
check "the partition-local string builder gives the same background" '[[ $rc -eq 0 ]]'
PANVAR_TEST_WAIT_AT=compact-locked:$OUT/go3 "${BG[@]}" -o s --compact & p3=$!
for _ in $(seq 1 200); do [[ -d s/.compact ]] && break; sleep 0.1; done
"${BG[@]}" -o s --compact 2> err.txt; rc=$?
touch go3; wait $p3; r3=$?
check "a second --compact while one runs is refused; the first completes" '[[ $rc -ne 0 && $r3 -eq 0 ]] && grep -q "in use" err.txt'
"${BG[@]}" -o s_m2 --compact && cp -r s_m2 s_m3 && cp s_m2.sshash s_m3.sshash
"${BG[@]}" -o s_m3 --add --genome g5=data/g5.fa
"${BG[@]}" -o s_m3 --verify s_m3.sshash 2> err.txt; rc=$?
check "a dictionary compacted before an --add is refused as stale" '[[ $rc -ne 0 ]] && grep -q stale err.txt'
"${BG[@]}" -o s_m3 --add --genome g5again=data/g5.fa 2> err.txt; rc=$?
check "a genome under a new name but already counted is refused" '[[ $rc -ne 0 ]] && grep -qE "already counted|same sequence" err.txt'

echo "== locating, and anchoring an existing store on a new reference"
"${BG[@]}" -o s --locate chr1:120000-160000 > loc.tsv
check "every genome is listed" '[[ $(cut -f1 loc.tsv | sort -u | grep -vc genome) == 5 ]]'
has_copy() {  # <file> <genome> <contig> <strand> <min start> <max start> <min end> <max end>
  awk -F'\t' -v g="$2" -v c="$3" -v s="$4" -v a="$5" -v b="$6" -v x="$7" -v y="$8" \
    '$1==g && $3==c && $6==s && $4>=a && $4<=b && $5>=x && $5<=y {f=1} END {exit !f}' "$1"
}
no_overlap() {  # <file> <genome>
  awk -F'\t' -v g="$2" '$1==g && $2>0 {print $3, $4, $5}' "$1" | sort -k1,1 -k2,2n |
    awk '$1==c && $2<e {bad=1} {c=$1; e=$3} END {exit bad}'
}
check "g1: the locus 20 kb further on" 'has_copy loc.tsv g1 chr1 + 134000 140000 180000 186000'
check "g2: a reversed copy on a, and the reference copy on b" 'has_copy loc.tsv g2 a - 62000 65500 104500 108000 && has_copy loc.tsv g2 b + 114000 120000 160000 166000'
check "g3: two tandem copies, not overlapping" 'no_overlap loc.tsv g3 && [[ $(awk -F"\t" '"'"'$1=="g3" && $2>0'"'"' loc.tsv | wc -l | tr -d " ") == 2 ]]'
"${BG[@]}" -o s --locate chr1:160000-120000 > /dev/null 2>&1; rc=$?
check "a reversed interval is refused" '[[ $rc -ne 0 ]]'
"${BG[@]}" -o s --locate chr1:290000-400000 > /dev/null 2>&1; rc=$?
check "an interval past the contig is refused" '[[ $rc -ne 0 ]]'
"${BG[@]}" -o r --genome ref=data/ref.fa --genome g1=data/g1.fa
"${BG[@]}" -o r --add --genome g2=data/g2.fa
"${BG[@]}" -o r --coordinate-reference ref
"${BG[@]}" -o r --locate chr1:120000-160000 > loc2.tsv
check "a reference added later anchors every committed genome" '[[ $(awk -F"\t" '"'"'NR>1 && $2>0 {print $1}'"'"' loc2.tsv | sort -u | wc -l | tr -d " ") == 3 ]]'
"${BG[@]}" -o r2 --genome ref=data/ref.fa --genome g1=data/g1.fa
PANVAR_TEST_WAIT_AT=count-before-commit:$OUT/go2 "${BG[@]}" -o r2 --add --genome g2=data/g2.fa & p2=$!
sleep 1
PANVAR_TEST_FAIL_AT=reference-before-commit "${BG[@]}" -o r2 --coordinate-reference ref; rc=$?
check "adding a reference killed before its commit: no reference" '[[ $rc -eq 137 ]] && ! grep -q fingerprint r2/manifest.json'
"${BG[@]}" -o r2 --coordinate-reference ref
touch go2; wait $p2; r2rc=$?
"${BG[@]}" -o r2 --locate chr1:120000-160000 > loc3.tsv
check "an --add counting while a reference was added maps itself on it and commits" '[[ $r2rc -eq 0 && $(awk -F"\t" '"'"'NR>1 && $2>0 {print $1}'"'"' loc3.tsv | sort -u | wc -l | tr -d " ") == 3 ]]'

echo "== inputs and limits"
"${BG[@]}" -o bad1 --genome a/b=data/g1.fa > /dev/null 2>&1; rc=$?
check "a genome name with a slash is refused" '[[ $rc -ne 0 ]]'
"${BG[@]}" -o bad2 --genome r=data/reads.fq 2> err.txt; rc=$?
check "FASTQ is refused as a genome" '[[ $rc -ne 0 ]] && grep -q FASTQ err.txt'
cat data/g1.fa data/g4.fa > data/twice.fa
"${BG[@]}" -o bad3 --genome t=data/twice.fa 2> err.txt; rc=$?
check "a genome naming a contig twice is refused" '[[ $rc -ne 0 ]] && grep -q "twice" err.txt'
printf ">e\nNNNN\n" > data/empty.fa
"${BG[@]}" -o e --genome e=data/empty.fa && "${BG[@]}" -o e --compact && "${BG[@]}" -o e --verify e.sshash; rc=$?
check "a store without k-mers builds, compacts and verifies" '[[ $rc -eq 0 ]]'
( ulimit -f 400; trap "" XFSZ; exec "${BG[@]}" -o full --genome ref=data/ref.fa --genome g1=data/g1.fa --no-dedup ) > /dev/null 2>&1 &
pid=$!
for _ in $(seq 1 120); do kill -0 $pid 2>/dev/null || break; sleep 0.5; done
if kill -0 $pid 2>/dev/null; then kill -9 $pid; hung=1; else hung=0; fi
wait $pid; rc=$?
check "out of disk space: the run fails instead of hanging" '[[ $hung -eq 0 && $rc -ne 0 ]]'

echo
echo "background_cli: $pass passed, $fail failed"
[[ $fail -eq 0 ]]
