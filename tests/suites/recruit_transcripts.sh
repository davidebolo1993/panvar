#!/usr/bin/env bash
# panvar recruit with a locus's transcripts in its index: a gene whose exons a paralog elsewhere shares (its
# exon k-mers shared), so only the transcript's junctions are private. The junction k-mers join the index
# flagged transcript-only and private, turn shared when another genome holds a processed pseudogene, and
# make reads across both junctions certain only with the transcripts; an index without them is unchanged.
#   recruit_transcripts.sh <panvar_bin> <tests/tools dir> <out_dir>
set -uo pipefail

if [[ $# -lt 3 ]]; then
  echo "Usage: $0 <panvar_bin> <tools_dir> <out_dir>" >&2
  exit 2
fi
BIN="$1"; TOOLS="$2"; OUT="$3"
rm -rf "$OUT"; mkdir -p "$OUT/tmp"; cd "$OUT" || exit 2
python3 "$TOOLS/junction_fixture.py" data
J=$(cat data/junction_kmers.txt)

pass=0; fail=0
ok() { echo "  ok   $1"; pass=$((pass + 1)); }
bad() { echo "  FAIL $1"; fail=$((fail + 1)); }
check() { if eval "$2"; then ok "$1"; else bad "$1"; fi; }
BG=("$BIN" background -t 2 --partitions 64 --memory 64M -q --tmp-dir "$OUT/tmp")
RC=("$BIN" recruit -t 2 --background s1 -r ref --index-only --no-calibration)
# transcript-only k-mers of a dump by class, and those of the locus's exons (not transcript-only)
tx() { awk -F'\t' 'NR>1 && $2!="D" && $8==1 {n[$2]++} END {printf "P%d S%d R%d", n["P"], n["S"], n["R"]}' "$1"; }
certain() { awk -F'\t' 'NR>1 && $3=="certain"' "$1" | wc -l | tr -d ' '; }

"$BIN" transcripts -i data/locus.gfa -r ref --gtf data/genes.gtf -o data/locus -q
"${BG[@]}" -o s1 --genome ref=data/ref.fa --genome ga=data/ga.fa --coordinate-reference ref
"${BG[@]}" -o s2 --genome ref=data/ref.fa --genome ga=data/ga.fa --genome gp=data/gp.fa --coordinate-reference ref

"${RC[@]}" --save-index idx0 --locus data/locus.gfa 2> log0.txt
"${RC[@]}" --save-index idx1 --locus data/locus.gfa --transcripts data/locus.hst.fa 2> log1.txt; rc=$?
"$BIN" recruit -t 2 --background s2 -r ref --index-only --no-calibration --save-index idx2 --locus data/locus.gfa \
  --transcripts locus=data/locus.hst.fa 2> log2.txt
for i in 0 1 2; do "$BIN" recruit --dump-index idx$i/locus.plx=dump$i.tsv; done
check "the index builds with the transcripts" '[[ $rc -eq 0 && -f idx1/locus.plx ]] && grep -q "k-mers only in transcripts: $J" log1.txt'
check "the two junctions add their $J transcript-only k-mers, all private ($(tx dump1.tsv))" '[[ "$(tx dump1.tsv)" == "P$J S0 R0" ]]'
check "the exons' k-mers are shared with the paralog" '[[ $(awk -F"\t" "NR>1 && \$8==0 && \$2==\"S\"" dump1.tsv | wc -l) -gt 1000 ]]'
check "each junction k-mer is on the node of its first base" '[[ $(awk -F"\t" "NR>1 && \$8==1 && \$7==\".\"" dump1.tsv | wc -l) -eq 0 ]]'
check "a processed pseudogene elsewhere makes them shared ($(tx dump2.tsv))" '[[ "$(tx dump2.tsv)" == "P0 S$J R0" ]] && [[ $(awk -F"\t" "NR>1 && \$8==1 && \$3!=1" dump2.tsv | wc -l) -eq 0 ]]'
check "the report counts them" 'grep -q "^locus	transcript_only	$J (private $J, shared 0, repeat 0)" idx1/locus.index.tsv'
check "an index without transcripts has no transcript-only k-mers" '[[ "$(tx dump0.tsv)" == "P0 S0 R0" ]] && ! grep -q transcript_only idx0/locus.index.tsv'
check "the genomic k-mers keep their labels" 'cmp -s <(awk -F"\t" "NR>1 && \$8==0" dump0.tsv | cut -f1-6 | sort) <(awk -F"\t" "NR>1 && \$8==0" dump1.tsv | cut -f1-6 | sort)'

"$BIN" recruit --locus idx0/locus.plx -i data/reads.fq -o r0.tsv --min-sites 2 -q
"$BIN" recruit --locus idx1/locus.plx -i data/reads.fq -o r1.tsv --min-sites 2 -q
"$BIN" recruit --locus idx2/locus.plx -i data/reads.fq -o r2.tsv --min-sites 2 -q
check "reads across both junctions are certain only with the transcripts ($(certain r0.tsv), $(certain r1.tsv), $(certain r2.tsv) of 200)" \
  '[[ $(certain r0.tsv) -eq 0 && $(certain r1.tsv) -eq 200 && $(certain r2.tsv) -eq 0 ]]'

"${RC[@]}" --save-index idx1 --locus data/locus.gfa --transcripts data/locus.hst.fa 2> reuse.txt
check "an index with the same transcripts is reused" 'grep -q "current, reused" reuse.txt'
"${RC[@]}" --save-index idx1 --locus data/locus.gfa 2> drop.txt
check "the same graph without the transcripts is rebuilt" 'grep -q "rebuilding" drop.txt'
"${RC[@]}" --save-index idx1 --locus data/locus.gfa --transcripts data/locus.hst.fa --transcripts-ok-only 2> okonly.txt
check "taking status ok only is another index" 'grep -q "rebuilding" okonly.txt && grep -q "status ok only" idx1/locus.index.tsv'
"$BIN" recruit -t 2 --background s1 -r ref --index-only --save-index idx4 --locus data/locus.gfa --transcripts data/locus.hst.fa 2> cal.txt; rc=$?
check "a calibrated index with transcripts builds" '[[ $rc -eq 0 ]] && grep -q "calibration: short reads" cal.txt'
"${RC[@]}" --save-index idx3 --locus data/locus.gfa --transcripts other=data/locus.hst.fa 2> bad.txt
check "transcripts for a locus not given are refused" 'grep -q "no --locus graph is named other" bad.txt'

echo "recruit_transcripts: $pass passed, $fail failed"
[[ $fail -eq 0 ]]
