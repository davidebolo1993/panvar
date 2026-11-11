#!/usr/bin/env bash
# Recruitment on real reads: a sample's pairs aligned to GRCh38 (a BAM, local or remote with its index) are
# fetched around each test locus (the locus and 200 kb on each side) and in the locus index's GRCh38 decoy
# windows, recruited against an index built without the sample, and compared with where the aligner put
# them (dev/recruit_real_score.py).
#   recruit_real.sh <bam> <sample> <out dir> [sim dir] [background store]
set -euo pipefail

BAM="$1"; SAMPLE="$2"; OUT="$3"
SIM="${4:-$HOME/panvar-data/sim}"
BG="${5:-$HOME/panvar-data/poc/poc.bgk}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
P="$ROOT/build/panvar"
RD="$ROOT/tests/real_data"
mkdir -p "$OUT"
cd "$OUT"
LOCI=("c4 $RD/paralogs/c4.gfa.gz" "cyp2d $RD/paralogs/cyp2d.gfa.gz" "gstm $RD/paralogs/gstm.gfa.gz"
      "acot $RD/paralogs/acot.gfa.gz" "lpa $RD/tandem_repeat/lpa.gfa.gz")
sets=()
for entry in "${LOCI[@]}"; do
  read -r L gfa <<< "$entry"
  echo "== $L $(date +%H:%M:%S)"
  mkdir -p "$L"
  # an index on the benchmark's graph, with the sample held out too when that graph has it
  if [[ ! -s "$L/$L.gfa" ]]; then
    if awk -v s="$SAMPLE" '($1 == "P" || $1 == "W") && index($2, s "#") == 1 {found = 1; exit} END {exit !found}' "$SIM/$L/$L.gfa"; then
      python3 "$ROOT/dev/recruit_sources.py" held-out --gfa "$SIM/$L/$L.gfa" --sample "$SAMPLE" --out "$L/$L.gfa"
    else
      cp "$SIM/$L/$L.gfa" "$L/$L.gfa"
    fi
  fi
  [[ -f "$L/$L.plx" ]] || "$P" recruit --locus "$L/$L.gfa" --background "$BG" -r grch38 --index-only --save-index "$L" -t 2 -q
  interval=$(awk -F'\t' '$1=="locus" && $2=="interval" {print $3}' "$L/$L.index.tsv")
  chrom=${interval%%:*}; span=${interval#*:}; s=${span%-*}; e=${span#*-}
  awk -F'\t' -v c="$chrom" -v s="$s" -v e="$e" 'BEGIN {OFS = "\t"; print c, (s > 200000 ? s - 200000 : 0), e + 200000}
    $1 == "decoy_window" && $2 == "grch38" {print $3, $4, $5}' "$L/$L.index.tsv" > "$L/regions.bed"
  if [[ ! -s "$L/reads_2.fq.gz" ]]; then
    if [[ ! -s "$L/fetched.bam" ]]; then
      samtools view -b -M -F 0xD04 -o "$L/fetching.bam" "$BAM" $(awk '{printf "%s:%d-%d ", $1, $2 + 1, $3}' "$L/regions.bed")
      mv "$L/fetching.bam" "$L/fetched.bam"
    fi
    samtools view -F 0x900 "$L/fetched.bam" | awk -F'\t' 'BEGIN {OFS = "\t"} {m = int($2 / 64) % 2 ? 1 : 2; c = $6; n = 0
      while (match(c, /[0-9]+[MDN=X]/)) {n += substr(c, RSTART, RLENGTH - 1); c = substr(c, RSTART + RLENGTH)}
      print $1, m, $3, $4, $5, $4 + n}' > "$L/aligned.tsv"
    samtools collate -u -O "$L/fetched.bam" "$L/collate" | samtools fastq -F 0x900 -1 "$L/reads_1.fq.gz" -2 "$L/reads_2.fq.gz" -s /dev/null -0 /dev/null -n - 2> /dev/null
    rm -f "$L/fetched.bam"
  fi
  "$P" recruit --locus "$L/$L.plx" -1 "$L/reads_1.fq.gz" -2 "$L/reads_2.fq.gz" -o "$L/recruited.tsv" -t 2 -q
  sets+=("$L" "$interval" "$L/regions.bed" "$L/aligned.tsv" "$L/reads_1.fq.gz" "$L/recruited.tsv")
done
python3 "$ROOT/dev/recruit_real_score.py" --out summary.tsv "${sets[@]}"
column -t summary.tsv
