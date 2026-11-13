#!/usr/bin/env bash
# Recruitment of a real sample scored against the sample's own assembly: every read is aligned to each
# of the sample's two haplotypes, where each locus lies is found by aligning the locus's GRCh38 interval
# to them, and each fragment recruit screened (and each fragment aligned into a locus) is classed by
# where its mates align (dev/recruit_assembly_score.py). The indexes should be built without the sample.
#   recruit_assembly_truth.sh <out dir> <reads> <hap1.fa[.gz]> <hap2.fa[.gz]> <GRCh38.fa with .fai> <index.plx>...
# <reads>: one CRAM, BAM or FASTQ (single reads, or pairs with INTERLEAVED=1), or R1,R2 for a FASTQ pair.
# THREADS (default 16), PRESET (sr, map-hifi or map-ont; default sr), CRAM_REFERENCE (default the GRCh38
# FASTA), MINIMAP2 and SAMTOOLS (executables), TMPDIR (temporary FASTQ of a CRAM or BAM: about a third of
# its size). A step that completed is marked (<output>.done) and skipped on a rerun; a failing step stops
# the run unmarked.
set -euo pipefail

abs() { echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"; }
mkdir -p "$1"
OUT="$(cd "$1" && pwd)"; READS="$2"; H1="$(abs "$3")"; H2="$(abs "$4")"; REF="$(abs "$5")"; shift 5
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
P="$ROOT/build/panvar"
T="${THREADS:-16}"; PRESET="${PRESET:-sr}"
MM2="${MINIMAP2:-minimap2}"; ST="${SAMTOOLS:-samtools}"
CRAM_REF="${CRAM_REFERENCE:-$REF}"
mkdir -p "${TMPDIR:-$OUT}"
TMP="$(mktemp -d "$(cd "${TMPDIR:-$OUT}" && pwd)/truth.XXXXXX")"
trap 'rm -rf "$TMP"' EXIT
plx=(); for x in "$@"; do plx+=("$(abs "$x")"); done
if [[ $READS == *,* ]]; then R1="$(abs "${READS%%,*}")"; R2="$(abs "${READS#*,}")"; IN=(-1 "$R1" -2 "$R2"); else READS="$(abs "$READS")"; IN=(-i "$READS"); fi
[[ ${INTERLEAVED:-0} == 1 ]] && IN+=(--interleaved)
cd "$OUT"
# the step runs as a plain command, so set -e stops the run inside it; the marker is written only after it
done_or_run() { local marker=$1; shift; [[ -f $marker ]] && return 0; "$@"; touch "$marker"; }
stamp() { echo "== $* $(date +%H:%M:%S)"; }

# 1. recruit
idx=(); for x in "${plx[@]}"; do idx+=(--locus "$x"); done
stamp recruit
run_recruit() { "$P" recruit "${idx[@]}" "${IN[@]}" --cram-reference "$CRAM_REF" -o recruit.tsv -t "$T" 2> recruit.log; }
done_or_run recruit.done run_recruit

# 2. each locus on the two haplotypes: its GRCh38 interval aligned to each (minimap2 asm5)
stamp loci on the assembly
map_loci() {
  : > loci.fa
  for x in "${plx[@]}"; do
    tsv="${x%.plx}.index.tsv"
    name=$(awk -F'\t' '$1=="locus" && $2=="name" {print $3}' "$tsv")
    iv=$(awk -F'\t' '$1=="locus" && $2=="interval" {print $3}' "$tsv")
    c=${iv%%:*}; s=${iv#*:}; s0=${s%-*}; e0=${s#*-}
    "$ST" faidx "$REF" "$c:$((s0 + 1))-$e0" | sed "1s/.*/>$name/" >> loci.fa
  done
  # each haplotype on its own: in one search the other haplotype's copy is a secondary and may be dropped
  : > loci.paf
  for h in 1 2; do
    hf=$H1; [[ $h == 2 ]] && hf=$H2
    "$MM2" -x asm5 -t "$T" "$hf" loci.fa 2>> loci_map.log | awk -v OFS='\t' -v p="h${h}_" '{$6 = p $6; print}' >> loci.paf
  done
  python3 "$ROOT/dev/recruit_assembly_score.py" intervals loci.fa loci.paf > loci.bed.tmp
  mv loci.bed.tmp loci.bed
}
done_or_run loci.done map_loci
cat loci.bed

# 3. every read aligned to each haplotype on its own (aligned to both at once, a read the two share would
# have MAPQ 0); kept: primary and supplementary alignments in a locus, and those of screened fragments
stamp align
align() {
  cut -f1 recruit.tsv | tail -n +2 | sort -u > screened.txt
  # sets of reads aligned on their own: a pair's two files; or single reads (with pairs, the reads whose
  # mate is missing)
  local sets=()
  if [[ ${IN[0]} == -1 ]]; then
    sets=("$R1"$'\t'"$R2")
  elif [[ $READS == *.cram || $READS == *.bam || $READS == *.sam ]]; then
    if [[ $PRESET == sr ]]; then
      "$ST" collate -@ "$T" -u -O --reference "$CRAM_REF" "$READS" "$TMP/collate" |
        "$ST" fastq -@ "$T" -n -1 "$TMP/r1.fq.gz" -2 "$TMP/r2.fq.gz" -0 "$TMP/r0.fq.gz" -s "$TMP/rs.fq.gz" -c 1 -
      sets=("$TMP/r1.fq.gz"$'\t'"$TMP/r2.fq.gz")
      local f
      for f in "$TMP/r0.fq.gz" "$TMP/rs.fq.gz"; do
        if [[ $( (gzip -cd "$f" | head -c 1 | wc -c) || true) -gt 0 ]]; then sets+=("$f"); fi
      done
    else
      "$ST" fastq -@ "$T" --reference "$CRAM_REF" -0 "$TMP/r.fq.gz" -c 1 "$READS" > /dev/null
      sets=("$TMP/r.fq.gz")
    fi
  else
    sets=("$READS")   # interleaved pairs are paired by minimap2's fragment mode (-x sr)
  fi
  : > aligned.tmp
  local h set files
  for h in 1 2; do
    local hf=$H1; [[ $h == 2 ]] && hf=$H2
    for set in "${sets[@]}"; do
    IFS=$'\t' read -r -a files <<< "$set"
    "$MM2" -ax "$PRESET" -t "$T" "$hf" "${files[@]}" 2>> align.log |
    "$ST" view -F 0x104 - |
    awk -F'\t' -v OFS='\t' -v p="h${h}_" 'FNR == NR && FILENAME == ARGV[1] {keep[$1] = 1; next}
      FILENAME == ARGV[2] {c[++n] = $1; s[n] = $2; e[n] = $3; next}
      {
        $3 = p $3; len = 0; cig = $6
        while (match(cig, /[0-9]+[MDN=X]/)) { len += substr(cig, RSTART, RLENGTH - 1); cig = substr(cig, RSTART + RLENGTH) }
        inside = 0
        for (i = 1; i <= n; i++) if ($3 == c[i] && $4 - 1 < e[i] && $4 - 1 + len > s[i]) inside = 1
        if (inside || ($1 in keep)) print $1, (int($2 / 128) % 2 ? 2 : 1), $3, $4 - 1, $4 - 1 + len, $5
      }' screened.txt loci.bed - >> aligned.tmp
    done
  done
  mv aligned.tmp aligned.tsv
}
done_or_run aligned.done align

# 4. scores
stamp score
python3 "$ROOT/dev/recruit_assembly_score.py" score loci.bed aligned.tsv recruit.tsv > truth.tsv.tmp
mv truth.tsv.tmp truth.tsv
column -t truth.tsv
