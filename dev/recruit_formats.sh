#!/usr/bin/env bash
# The input layer on real alignments: a FASTQ pair is aligned (minimap2 -ax sr) and written as a
# coordinate-sorted BAM with index, the same without index, a name-sorted BAM, a CRAM with index and an
# unaligned BAM; each is recruited at 1 and 8 threads and its table, sorted, compared with the FASTQ's.
# Seconds per run in times.tsv.
#   recruit_formats.sh <reads_1.fq.gz> <reads_2.fq.gz> <reference.fa(.gz), with .fai> <out dir> <index.plx>...
# minimap2 and samtools on PATH (MINIMAP2=... to use another minimap2, MM2_INDEX=... a prebuilt index of
# the reference).
set -euo pipefail

abs() { echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"; }
R1="$(abs "$1")"; R2="$(abs "$2")"; REF="$(abs "$3")"; OUT="$4"; shift 4
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
P="$ROOT/build/panvar"
MM2="${MINIMAP2:-minimap2}"
idx=(); for x in "$@"; do idx+=(--locus "$(abs "$x")"); done
mkdir -p "$OUT"
cd "$OUT"

if [[ ! -s coord.bam.bai ]]; then
  "$MM2" -ax sr -t 4 -R '@RG\tID:r\tSM:sample' "${MM2_INDEX:-$REF}" "$R1" "$R2" 2> map.log | samtools sort -@ 2 -m 1G -o coord.bam -
  samtools index coord.bam
  cp coord.bam noindex.bam
  samtools sort -n -@ 2 -o name.bam coord.bam
  samtools view -C -T "$REF" -@ 2 -o coord.cram coord.bam
  samtools index coord.cram
  samtools import -1 "$R1" -2 "$R2" -o unaligned.bam 2> /dev/null
fi

sorted() { { head -1 "$1"; tail -n +2 "$1" | sort; } > "$1.sorted"; }
echo -e "input\tthreads\tseconds\tsame_as_fastq" > times.tsv
run() {   # run <label> <threads> <args...>
  local label=$1 t=$2 t0 t1 same; shift 2
  t0=$(date +%s.%N 2>/dev/null || python3 -c 'import time; print(time.time())')
  "$P" recruit "${idx[@]}" "$@" -o "$label.$t.tsv" -t "$t" 2> "$label.$t.log"
  t1=$(date +%s.%N 2>/dev/null || python3 -c 'import time; print(time.time())')
  sorted "$label.$t.tsv"
  same=yes; cmp -s "$label.$t.tsv.sorted" fastq.1.tsv.sorted || same=no
  echo -e "$label\t$t\t$(python3 -c "print(round($t1 - $t0, 2))")\t$same" >> times.tsv
}
run fastq 1 -1 "$R1" -2 "$R2"
for t in 1 8; do
  [[ $t == 8 ]] && run fastq 8 -1 "$R1" -2 "$R2"
  run coord $t -i coord.bam
  run noindex $t -i noindex.bam
  run name $t -i name.bam
  run cram $t -i coord.cram --cram-reference "$REF"
  run unaligned $t -i unaligned.bam
done
column -t times.tsv
grep -h "mates fetched\|second pass" coord.8.log noindex.8.log | cut -c20- || true
