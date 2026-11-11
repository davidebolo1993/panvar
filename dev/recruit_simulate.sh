#!/usr/bin/env bash
# Simulated reads for the recruitment benchmark. For each test locus a sample's paths, and its relatives',
# are held out of the graph; its index is built without them, against two backgrounds:
#   main     GRCh38 + CHM13 (the proof-of-concept store)
#   g38      GRCh38 only, so CHM13 is a genome the index has never seen (its reads are new negatives)
# Reads come from these sources, each to its own file so the file is the truth:
#   locus     the held-out sample's two haplotypes (15x each short, 10x each long)
#   near      GRCh38 2 Mb on each side of the locus (30x short, 20x long)
#   decoy     the main index's decoy windows in GRCh38 (30x short); for long reads widened by 10 kb of
#             context on each side, never into the locus (20x)
#   boundary  15 kb on each side of each locus edge, in GRCh38 (30x short, 20x long): partly locus
#   genome    the GRCh38 chromosomes of the background with the five loci masked (1x), shared by all loci
#   new       the CHM13 chromosomes with the five loci's CHM13 copies masked (1x): negatives for g38
# Short reads: wgsim, 2x150, inserts 400 +- 50, 0.2% errors, no added variants. Long reads: pbsim3 ONT
# (ERRHMM-ONT-HQ; about 96% accuracy comes out), and HiFi-like reads from dev/simulate_hifi.py (0.1%
# errors; pbsim3 makes HiFi only through PacBio's ccs and caps one pass at 99%), 15 kb mean; long reads
# that are mostly N (masked sequence) are dropped.
#   recruit_simulate.sh [out_dir] [background store] [GRCh38 fasta] [CHM13 fasta] [pbsim3 dir]
set -euo pipefail

OUT="${1:-$HOME/panvar-data/sim}"
BG="${2:-$HOME/panvar-data/poc/poc.bgk}"
REF="${3:-$HOME/panvar-data/poc/grch38_poc.fa.gz}"
T2T="${4:-$HOME/panvar-data/poc/chm13_poc.fa.gz}"
PB="${5:-$HOME/tools/pbsim3}"
[[ -z "$OUT" ]] && OUT="$HOME/panvar-data/sim"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
P="$ROOT/build/panvar"
SRC="python3 $ROOT/dev/recruit_sources.py"
RD="$ROOT/tests/real_data"
mkdir -p "$OUT"
cd "$OUT"

# locus, graph, held-out sample, everyone held out (the sample and its relatives in the graph)
LOCI=("c4 $RD/paralogs/c4.gfa.gz HG00733 HG00733,HG00731,HG00732"
      "cyp2d $RD/paralogs/cyp2d.gfa.gz HG00733 HG00733,HG00731,HG00732"
      "gstm $RD/paralogs/gstm.gfa.gz HG002 HG002" "acot $RD/paralogs/acot.gfa.gz HG002 HG002"
      "lpa $RD/tandem_repeat/lpa.gfa.gz HG002 HG002")

short() {   # short <fasta> <coverage> <prefix> <seed>
  local len n
  len=$(awk '!/^>/ {n += length($0)} END {print n}' "$1")
  n=$(( len * $2 / 300 ))
  wgsim -N "$n" -1 150 -2 150 -d 400 -s 50 -e 0.002 -r 0 -R 0 -X 0 -S "$4" "$1" "$3_1.fq" "$3_2.fq" > /dev/null
  gzip -f "$3_1.fq" "$3_2.fq"
}
long() {   # long <ont|hifi> <fasta> <depth> <prefix> <seed>
  if [[ $1 == hifi ]]; then
    python3 "$ROOT/dev/simulate_hifi.py" "$2" "$3" "$4.fq.gz" --seed "$5" --prefix "$(basename "$4")_"
    return
  fi
  local model=(--method errhmm --errhmm "$PB/data/ERRHMM-ONT-HQ.model" --accuracy-mean 0.98 --length-mean 15000 --length-sd 10000)
  # pbsim3 stops with SIGPIPE on its second sequence here, so each sequence is simulated on its own
  local parts="$4.parts" i=0 f
  mkdir -p "$parts"
  awk -v d="$parts" '/^>/ {n++} {print > (d "/" n ".fa")}' "$2"
  : > "$4.log"
  for f in "$parts"/*.fa; do
    i=$((i + 1))
    "$PB/src/pbsim" --strategy wgs --genome "$f" --depth "$3" --seed $(( $5 * 1000 + i )) --prefix "$parts/out$i" \
      --id-prefix "$(basename "$4")_$i." "${model[@]}" >> "$4.log" 2>&1
  done
  for f in "$parts"/out*_0001.fq.gz; do zcat < "$f"; done |
    awk 'NR % 4 == 1 {h = $0} NR % 4 == 2 {s = $0} NR % 4 == 3 {p = $0} NR % 4 == 0 {n = gsub(/[Nn]/, "N", s); if (n * 2 <= length(s)) print h "\n" s "\n" p "\n" $0}' |
    gzip > "$4.fq.gz"
  for f in "$parts"/*; do rm -f "$f"; done
  rmdir "$parts"
}

# the GRCh38-only background (the g38 arm)
BG38="$OUT/bg_g38.bgk"
[[ -f "$BG38/manifest.json" ]] || "$P" background -o "$BG38" --genome grch38="$REF" --coordinate-reference grch38 -t 2 --memory 4G -q

masks=(); masks_t2t=()
seed=1
for entry in "${LOCI[@]}"; do
  read -r name gfa sample out_samples <<< "$entry"
  echo "== $name (held out: $out_samples) $(date +%H:%M:%S)"
  mkdir -p "$name/g38"
  $SRC held-out --gfa "$gfa" --sample "$out_samples" --out "$name/$name.gfa"
  "$P" recruit --locus "$name/$name.gfa" --background "$BG" -r grch38 --index-only --save-index "$name" -t 2 -q
  "$P" recruit --locus "$name/$name.gfa" --background "$BG38" -r grch38 --index-only --save-index "$name/g38" -t 2 -q
  interval=$(awk -F'\t' '$1=="locus" && $2=="interval" {print $3}' "$name/$name.index.tsv")
  masks+=(--mask "$interval")
  while read -r c s e; do masks_t2t+=(--mask "$c:$s-$e"); done < <(awk -F'\t' '$1=="span" && $2=="chm13" {print $3, $4, $5}' "$name/$name.index.tsv")
  $SRC locus --gfa "$gfa" --sample "$sample" --out "$name/locus.fa"
  $SRC near --fasta "$REF" --region "$interval" --out "$name/near.fa"
  $SRC decoy --fasta "$REF" --genome-name grch38 --report "$name/$name.index.tsv" --out "$name/decoy.fa"
  $SRC decoy --fasta "$REF" --genome-name grch38 --report "$name/$name.index.tsv" --context 10000 --out "$name/decoy_long.fa"
  $SRC boundary --fasta "$REF" --region "$interval" --flank 15000 --out "$name/boundary.fa"
  for src in locus near decoy boundary; do
    cov=30; dep=20; [[ $src == locus ]] && { cov=15; dep=10; }
    lsrc=$src; [[ $src == decoy ]] && lsrc=decoy_long
    short "$name/$src.fa" $cov "$name/short_$src" $((seed++))
    long ont "$name/$lsrc.fa" $dep "$name/ont_$src" $((seed++))
    long hifi "$name/$lsrc.fa" $dep "$name/hifi_$src" $((seed++))
  done
done
echo "== genome and new genome $(date +%H:%M:%S)"
mkdir -p genome new
$SRC genome --fasta "$REF" "${masks[@]}" --out genome/genome.fa
$SRC genome --fasta "$T2T" "${masks_t2t[@]}" --out new/genome.fa
for g in genome new; do
  short $g/genome.fa 1 $g/short_$g $((seed++))
  long ont $g/genome.fa 1 $g/ont_$g $((seed++))
  long hifi $g/genome.fa 1 $g/hifi_$g $((seed++))
done
echo "== done $(date +%H:%M:%S)"
du -sh "$OUT"/* | sed 's/^/   /'
