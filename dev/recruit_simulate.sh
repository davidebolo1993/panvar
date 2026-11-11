#!/usr/bin/env bash
# Simulated reads for the recruitment benchmark. For each test locus a sample's paths, and its relatives',
# are held out of the graph; its index is built without them, against these backgrounds:
#   main     GRCh38 + CHM13 (the proof-of-concept store)
#   g38      GRCh38 only, so CHM13 is a genome the index has never seen (its reads are new negatives)
#   wg       with WG=1 and WG_REF: the whole GRCh38 (the other two hold only some chromosomes)
# Reads come from these sources, each to its own file so the file is the truth:
#   locus     the held-out sample's two haplotypes (15x each short, 10x each long)
#   near      GRCh38 2 Mb on each side of the locus (30x short, 20x long)
#   decoy     the main index's decoy windows in GRCh38 (30x short); for long reads widened by 10 kb of
#             context on each side, never into the locus (20x)
#   boundary  15 kb on each side of each locus edge, in GRCh38 (30x short, 20x long): partly locus
#   genome    the GRCh38 chromosomes of the background with the five loci masked (1x), shared by all loci
#   new       the CHM13 chromosomes with the five loci's CHM13 copies masked (1x): negatives for g38
#   wgnew     with WG=1, the whole CHM13 with the loci masked (0.5x ART pairs, 0.2x long): negatives with
#             every chromosome, most of them in no background but wg's
# Technologies:
#   short     wgsim, 2x150, inserts 400 +- 50, 0.2% uniform substitutions
#   illumina  ART's HiSeq 2500 profile, 2x150, inserts 400 +- 50: errors rising along the read, indels
#   ont       pbsim3 ERRHMM-ONT-HQ (about 96% accuracy comes out), 15 kb mean
#   hifi      dev/simulate_hifi.py, 0.1% errors (pbsim3 makes HiFi only through PacBio's ccs and caps one
#             pass at 99%), 15 kb mean
# ART pairs and ONT reads carry their source span in their names (dev/sim_truth.py), wgsim and HiFi reads
# their own way; long reads that are mostly N (masked sequence) are dropped. Sets already made are kept
# (boundary sets made before spans were in the names must be removed to be remade).
#   recruit_simulate.sh [out dir] [background store] <GRCh38 fasta> <CHM13 fasta> [pbsim3 dir]
# The out dir defaults to tests/examples/recruit_sim, the store to <out dir>/bg.bgk (counted from the two
# genomes when missing). pbsim3: a source tree (src/pbsim, data/), else PBSIM_DIR, else pbsim on PATH with
# PBSIM_MODEL the ERRHMM-ONT-HQ model file.
#   WG=1 WG_T2T=<whole CHM13 fasta> [WG_REF=<whole GRCh38 fasta>] recruit_simulate.sh ...
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
abs() { echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"; }
OUT="${1:-$ROOT/tests/examples/recruit_sim}"
BG="${2:-}"
REF="${3:-}"; T2T="${4:-}"
[[ -n "$REF" && -n "$T2T" ]] || { echo "recruit_simulate.sh: give the GRCh38 and CHM13 FASTAs (third and fourth arguments)" >&2; exit 1; }
REF="$(abs "$REF")"; T2T="$(abs "$T2T")"
PB="${5:-${PBSIM_DIR:-}}"
PBSIM="${PB:+$PB/src/pbsim}"; PBSIM="${PBSIM:-pbsim}"
PBSIM_MODEL="${PBSIM_MODEL:-${PB:+$PB/data/ERRHMM-ONT-HQ.model}}"
ART="${ART:-art_illumina}"
WG="${WG:-0}"
WG_REF="${WG_REF:-}"
WG_T2T="${WG_T2T:-}"
mkdir -p "$OUT"
OUT="$(cd "$OUT" && pwd)"
BG="${BG:+$(abs "$BG")}"; BG="${BG:-$OUT/bg.bgk}"
P="$ROOT/build/panvar"
SRC="python3 $ROOT/dev/recruit_sources.py"
TRUTH="python3 $ROOT/dev/sim_truth.py"
RD="$ROOT/tests/real_data"
cd "$OUT"
# the GRCh38 + CHM13 background of the main arm
[[ -f "$BG/manifest.json" ]] || "$P" background -o "$BG" --genome grch38="$REF" --genome chm13="$T2T" --coordinate-reference grch38 -t 2 --memory 4G -q

# locus, graph, held-out sample, everyone held out (the sample and its relatives in the graph)
LOCI=("c4 $RD/paralogs/c4.gfa.gz HG00733 HG00733,HG00731,HG00732"
      "cyp2d $RD/paralogs/cyp2d.gfa.gz HG00733 HG00733,HG00731,HG00732"
      "gstm $RD/paralogs/gstm.gfa.gz HG002 HG002" "acot $RD/paralogs/acot.gfa.gz HG002 HG002"
      "lpa $RD/tandem_repeat/lpa.gfa.gz HG002 HG002")

plain() {   # plain <fasta(.gz)> <scratch>: a path to the sequence uncompressed
  if [[ $1 == *.gz ]]; then zcat < "$1" > "$2"; echo "$2"; else echo "$1"; fi
}
short() {   # short <fasta> <coverage> <prefix> <seed>
  [[ -s "$3_2.fq.gz" ]] && return
  local fa len n
  fa=$(plain "$1" "$3.tmp.fa")
  len=$(awk '!/^>/ {n += length($0)} END {print n}' "$fa")
  n=$(awk -v l="$len" -v c="$2" 'BEGIN {printf "%d", l * c / 300}')
  wgsim -N "$n" -1 150 -2 150 -d 400 -s 50 -e 0.002 -r 0 -R 0 -X 0 -S "$4" "$fa" "$3_1.fq" "$3_2.fq" > /dev/null
  gzip -f "$3_1.fq" "$3_2.fq"
  rm -f "$3.tmp.fa"
}
illumina() {   # illumina <fasta> <coverage> <prefix> <seed>
  [[ -s "$3_2.fq.gz" ]] && return
  local fa
  fa=$(plain "$1" "$3.tmp.fa")
  "$ART" -ss HS25 -i "$fa" -p -l 150 -f "$2" -m 400 -s 50 -rs "$4" -o "$3.art_" -sam -na -q > "$3.art.log" 2>&1
  $TRUTH art "$3.art_" "$3_1.fq.gz" "$3_2.fq.gz"
  rm -f "$3.art_1.fq" "$3.art_2.fq" "$3.art_.sam" "$3.tmp.fa"
}
long() {   # long <ont|hifi> <fasta> <depth> <prefix> <seed>
  [[ -s "$4.fq.gz" ]] && return
  if [[ $1 == hifi ]]; then
    python3 "$ROOT/dev/simulate_hifi.py" "$2" "$3" "$4.fq.gz" --seed "$5" --prefix "$(basename "$4")_"
    return
  fi
  [[ -n "$PBSIM_MODEL" ]] || { echo "recruit_simulate.sh: set PBSIM_MODEL to pbsim3's ERRHMM-ONT-HQ.model (or give its source tree)" >&2; exit 1; }
  local model=(--method errhmm --errhmm "$PBSIM_MODEL" --accuracy-mean 0.98 --length-mean 15000 --length-sd 10000)
  # pbsim3 stops with SIGPIPE on its second sequence here, so each sequence is simulated on its own
  local parts="$4.parts" i=0 f name
  mkdir -p "$parts"
  if [[ $2 == *.gz ]]; then zcat < "$2"; else cat "$2"; fi | awk -v d="$parts" '/^>/ {n++} {print > (d "/" n ".fa")}'
  : > "$4.log"
  : > "$4.fq"
  for f in "$parts"/*.fa; do
    i=$((i + 1))
    name=$(head -1 "$f" | cut -c2- | cut -d' ' -f1)
    "$PBSIM" --strategy wgs --genome "$f" --depth "$3" --seed $(( $5 * 1000 + i )) --prefix "$parts/out$i" \
      --id-prefix "$(basename "$4")_$i." "${model[@]}" >> "$4.log" 2>&1
    if [[ -s "$parts/out${i}_0001.fq.gz" ]]; then $TRUTH maf "$parts/out${i}_0001.fq.gz" "$parts/out${i}_0001.maf.gz" "$name" >> "$4.fq"; fi
    rm -f "$parts/out${i}"_* "$f"
  done
  gzip -f "$4.fq"
  rmdir "$parts"
}

# the GRCh38-only background (the g38 arm), only while an index still needs it
BG38="$OUT/bg_g38.bgk"
need38=0
for entry in "${LOCI[@]}"; do read -r name _ <<< "$entry"; [[ -f "$name/g38/$name.plx" ]] || need38=1; done
if [[ $need38 == 1 && ! -f "$BG38/manifest.json" ]]; then
  "$P" background -o "$BG38" --genome grch38="$REF" --coordinate-reference grch38 -t 2 --memory 4G -q
fi

masks=(); masks_t2t=()
seed=1
for entry in "${LOCI[@]}"; do
  read -r name gfa sample out_samples <<< "$entry"
  echo "== $name (held out: $out_samples) $(date +%H:%M:%S)"
  mkdir -p "$name/g38"
  [[ -s "$name/$name.gfa" ]] || $SRC held-out --gfa "$gfa" --sample "$out_samples" --out "$name/$name.gfa"
  [[ -f "$name/$name.plx" ]] || "$P" recruit --locus "$name/$name.gfa" --background "$BG" -r grch38 --index-only --save-index "$name" -t 2 -q
  [[ -f "$name/g38/$name.plx" ]] || "$P" recruit --locus "$name/$name.gfa" --background "$BG38" -r grch38 --index-only --save-index "$name/g38" -t 2 -q
  interval=$(awk -F'\t' '$1=="locus" && $2=="interval" {print $3}' "$name/$name.index.tsv")
  masks+=(--mask "$interval")
  while read -r c s e; do masks_t2t+=(--mask "$c:$s-$e"); done < <(awk -F'\t' '$1=="span" && $2=="chm13" {print $3, $4, $5}' "$name/$name.index.tsv")
  [[ -s "$name/locus.fa" ]] || $SRC locus --gfa "$gfa" --sample "$sample" --out "$name/locus.fa"
  [[ -s "$name/near.fa" ]] || $SRC near --fasta "$REF" --region "$interval" --out "$name/near.fa"
  [[ -s "$name/decoy.fa" ]] || $SRC decoy --fasta "$REF" --genome-name grch38 --report "$name/$name.index.tsv" --out "$name/decoy.fa"
  [[ -s "$name/decoy_long.fa" ]] || $SRC decoy --fasta "$REF" --genome-name grch38 --report "$name/$name.index.tsv" --context 10000 --out "$name/decoy_long.fa"
  [[ -s "$name/boundary.fa" ]] || $SRC boundary --fasta "$REF" --region "$interval" --flank 15000 --out "$name/boundary.fa"
  for src in locus near decoy boundary; do
    cov=30; dep=20; [[ $src == locus ]] && { cov=15; dep=10; }
    lsrc=$src; [[ $src == decoy ]] && lsrc=decoy_long
    short "$name/$src.fa" $cov "$name/short_$src" $((seed++))
    long ont "$name/$lsrc.fa" $dep "$name/ont_$src" $((seed++))
    long hifi "$name/$lsrc.fa" $dep "$name/hifi_$src" $((seed++))
    illumina "$name/$src.fa" $cov "$name/illumina_$src" $((seed + 500))
  done
done
echo "== genome and new genome $(date +%H:%M:%S)"
mkdir -p genome new
[[ -s genome/genome.fa.gz ]] || { $SRC genome --fasta "$REF" "${masks[@]}" --out genome/genome.fa; gzip genome/genome.fa; }
[[ -s new/genome.fa.gz ]] || { $SRC genome --fasta "$T2T" "${masks_t2t[@]}" --out new/genome.fa; gzip new/genome.fa; }
for g in genome new; do
  short $g/genome.fa.gz 1 $g/short_$g $((seed++))
  long ont $g/genome.fa.gz 1 $g/ont_$g $((seed++))
  long hifi $g/genome.fa.gz 1 $g/hifi_$g $((seed++))
  illumina $g/genome.fa.gz 1 $g/illumina_$g $((seed + 500))
done

if [[ $WG == 1 ]]; then
  echo "== whole genome $(date +%H:%M:%S)"
  [[ -n $WG_T2T ]] || { echo "WG=1 needs WG_T2T" >&2; exit 1; }
  if [[ -n $WG_REF ]]; then   # the whole-GRCh38 background (about 15 GB, and 3 GB more while counting)
    BGW="$OUT/bg_wg.bgk"
    needw=0
    for entry in "${LOCI[@]}"; do read -r name _ <<< "$entry"; [[ -f "$name/wg/$name.plx" ]] || needw=1; done
    if [[ $needw == 1 && ! -f "$BGW/manifest.json" ]]; then
      "$P" background -o "$BGW" --genome grch38="$WG_REF" --coordinate-reference grch38 -t 2 --memory 6G -q
    fi
    for entry in "${LOCI[@]}"; do
      read -r name _ <<< "$entry"
      mkdir -p "$name/wg"
      [[ -f "$name/wg/$name.plx" ]] || "$P" recruit --locus "$name/$name.gfa" --background "$BGW" -r grch38 --index-only --save-index "$name/wg" -t 2 -q
    done
  fi
  mkdir -p wgnew
  [[ -s wgnew/genome.fa.gz ]] || { $SRC genome --fasta "$WG_T2T" "${masks_t2t[@]}" --out wgnew/genome.fa; gzip -1 wgnew/genome.fa; }
  illumina wgnew/genome.fa.gz 0.5 wgnew/illumina_wgnew $((seed + 700))
  long ont wgnew/genome.fa.gz 0.2 wgnew/ont_wgnew $((seed + 701))
  long hifi wgnew/genome.fa.gz 0.2 wgnew/hifi_wgnew $((seed + 702))
fi
echo "== done $(date +%H:%M:%S)"
du -sh "$OUT"/* | sed 's/^/   /'
