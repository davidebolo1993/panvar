#!/usr/bin/env bash
# Recruitment of a real sample scored against the sample's own assembly: every read is aligned to each
# of the sample's two haplotypes, where each locus lies is found by aligning the locus's GRCh38 interval
# to them, and each fragment recruit screened (and each fragment aligned into a locus) is classed by
# where its mates align (dev/recruit_assembly_score.py). The indexes should be built without the sample.
#   recruit_assembly_truth.sh <out dir> <reads> <hap1.fa[.gz]> <hap2.fa[.gz]> <GRCh38.fa with .fai> <index.plx>...
# <reads>: one CRAM, BAM or FASTQ (single or interleaved), or R1,R2 for a FASTQ pair.
# THREADS (default 16), PRESET (sr, map-hifi or map-ont; default sr), CRAM_REFERENCE (default the GRCh38
# FASTA), MINIMAP2 and SAMTOOLS (executables), TMPDIR (temporary FASTQ of a CRAM or BAM: about a third of
# its size). A step that completed is marked (<output>.done) and skipped on a rerun.
set -euo pipefail

abs() { echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"; }
OUT="$1"; READS="$2"; H1="$(abs "$3")"; H2="$(abs "$4")"; REF="$(abs "$5")"; shift 5
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
P="$ROOT/build/panvar"
T="${THREADS:-16}"; PRESET="${PRESET:-sr}"
MM2="${MINIMAP2:-minimap2}"; ST="${SAMTOOLS:-samtools}"
CRAM_REF="${CRAM_REFERENCE:-$REF}"
TMP="${TMPDIR:-$OUT/tmp}"
plx=(); for x in "$@"; do plx+=("$(abs "$x")"); done
if [[ $READS == *,* ]]; then R1="$(abs "${READS%%,*}")"; R2="$(abs "${READS#*,}")"; IN=(-1 "$R1" -2 "$R2"); else READS="$(abs "$READS")"; IN=(-i "$READS"); fi
mkdir -p "$OUT" "$TMP"
cd "$OUT"
done_or_run() { local marker=$1; shift; [[ -f $marker ]] && return 0; "$@" && touch "$marker"; }
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
    "$MM2" -x asm5 -t "$T" --secondary=yes -N 5 "$hf" loci.fa 2>> loci_map.log | awk -v OFS='\t' -v p="h${h}_" '{$6 = p $6; print}' >> loci.paf
  done
  python3 "$ROOT/dev/recruit_assembly_score.py" intervals loci.fa loci.paf > loci.bed
}
done_or_run loci.done map_loci
cat loci.bed

# 3. every read aligned to each haplotype on its own (aligned to both at once, a read the two share would
# have MAPQ 0); kept: primary alignments in a locus, and those of screened fragments
stamp align
align() {
  cut -f1 recruit.tsv | tail -n +2 | sort -u > screened.txt
  local src
  if [[ ${IN[0]} == -1 ]]; then
    src=("$R1" "$R2")
  elif [[ $READS == *.cram || $READS == *.bam || $READS == *.sam ]]; then
    if [[ $PRESET == sr ]]; then
      "$ST" collate -@ "$T" -u -O --reference "$CRAM_REF" "$READS" "$TMP/collate" |
        "$ST" fastq -@ "$T" -n -1 "$TMP/r1.fq.gz" -2 "$TMP/r2.fq.gz" -0 /dev/null -s /dev/null -c 1 -
      src=("$TMP/r1.fq.gz" "$TMP/r2.fq.gz")
    else
      "$ST" fastq -@ "$T" --reference "$CRAM_REF" -0 "$TMP/r.fq.gz" -c 1 "$READS" > /dev/null
      src=("$TMP/r.fq.gz")
    fi
  else
    src=("$READS")
  fi
  : > aligned.tmp
  for h in 1 2; do
    local hf=$H1; [[ $h == 2 ]] && hf=$H2
    "$MM2" -ax "$PRESET" -t "$T" "$hf" "${src[@]}" 2>> align.log |
    "$ST" view -F 0x904 - |
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
  mv aligned.tmp aligned.tsv
  rm -f "$TMP/r1.fq.gz" "$TMP/r2.fq.gz" "$TMP/r.fq.gz"
}
done_or_run aligned.done align

# 4. scores
stamp score
python3 "$ROOT/dev/recruit_assembly_score.py" score loci.bed aligned.tsv recruit.tsv > truth.tsv
column -t truth.tsv
