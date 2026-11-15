#!/usr/bin/env bash
# RNA recruitment benchmark. Per locus, groups of held-out samples (and their relatives) are removed from
# the graph; reads come from their haplotype-specific transcripts (projected onto their own walks), from
# molecules the annotation lacks, and from negatives elsewhere in the transcriptome (dev/rna_sim.py).
# The same reads are recruited in each arm:
#   a  genomic walks only
#   b  + the held-out graph's transcripts (panvar transcripts)
#   c  + a transcriptome segment in the background
#   d  as c, with the haplotypes' own annotations merged (loci given CAT_<locus>)
# Arms a and b are indexed before the transcriptome joins the store, c and d after; reads are recruited
# with each arm's .plx as built. Steps already done are kept, so a stopped run resumes.
#   rna_simulate.sh <out dir> <store> <transcriptome.fa.gz> [pbsim3 dir]
# The store holds the genomes only (no transcriptome yet, unless every a/b index exists); the
# transcriptome comes from `rna_sim.py transcriptome` and is added to the store as `ens`.
# Environment:
#   LOCI      loci to run, from the table below (default: all)
#   NGROUPS   held-out groups per locus (4), GROUP_SIZE samples per group (2)
#   N_SHORT N_LONG N_ONT   locus fragments per group and class (3000, 1000, 1000)
#   N_PARALOG N_PSEUDO N_TX  negative fragments per group and class (1000, 500, 3000)
#   CHROMS    negatives only from these chromosomes, for a background of some chromosomes
#   CAT_<locus>=<tsv>  haplotypes' annotations (transcripts --hap-annotations) for arm d
#   CONTIG_ALIAS  --contig-alias for those annotations
#   THREADS   (2)
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
abs() { echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"; }
[[ $# -ge 3 ]] || { echo "usage: rna_simulate.sh <out dir> <store> <transcriptome.fa.gz> [pbsim3 dir]" >&2; exit 1; }
mkdir -p "$1"
OUT="$(cd "$1" && pwd)"; BG="$(abs "$2")"; TX="$(abs "$3")"
PB="${4:-${PBSIM_DIR:-}}"
PBSIM="${PB:+$PB/src/pbsim}"; PBSIM="${PBSIM:-pbsim}"
PBSIM_MODEL="${PBSIM_MODEL:-${PB:+$PB/data/ERRHMM-ONT-HQ.model}}"
P="$ROOT/build/panvar"
SIM="python3 $ROOT/dev/rna_sim.py"
RD="$ROOT/tests/real_data"
GTF="$RD/Homo_sapiens.GRCh38.116.loci.gtf.gz"
T="${THREADS:-2}"
NGROUPS="${NGROUPS:-4}"; GROUP_SIZE="${GROUP_SIZE:-2}"
TIMES="$OUT/times.tsv"
[[ -f $TIMES ]] || printf 'step\tlocus\tgroup\tarm\tseconds\n' > "$TIMES"

# locus, graph, the graph's reference path for transcripts, symbols of the locus genes' families elsewhere
# (paralog negatives; genes inside the locus are never negatives)
LOCUS_TABLE=(
  "gstm paralogs/gstm.gfa.gz GRCh38 ^(GST|MGST|AMPD[0-9]|EPS8|GNA[IOTZ]|GNAT|RPL7$)"
  "acot paralogs/acot.gfa.gz GRCh38 ^(ACOT|BAAT|THEM[0-9]|NT5C|NDUFB8|RIOX|HEATR|DNAL|NUMB)"
  "c4 paralogs/c4.gfa.gz grch38 ^(C3|C5|C2|CFB|A2M|A2ML1|PZP|CD109|CPAMD8|CYP21A2|TNXB|EHMT1|ATF6|MTREX|NELFE|DXO)$"
  "cyp2d paralogs/cyp2d.gfa.gz grch38 ^(CYP2|NDUFA6|NAGA$|SMDT1|OLA1|TCF20|WBP2)"
  "lpa tandem_repeat/lpa.gfa.gz GRCh38 ^(PLG|LPA|HGF|HGFAC|MST1|F12|PLAT|PLAU|KREMEN)"
  "hla_a hla/hla_a.gfa.gz GRCh38 ^(HLA-|MIC[A-E]$|HCG4)"
  "hla_b hla/hla_b.gfa.gz GRCh38 ^(HLA-|MIC[A-E]$|HCG4)"
)
# held-out samples are taken in this order, then the graph's others; relatives leave the graph with them
PREFERRED="HG002 HG00733 NA19240 HG005 HG02257 HG01891 HG03516 NA18906 HG00514 NA12878"
relatives() {
  case $1 in
    HG002) echo HG003,HG004 ;; HG00733) echo HG00731,HG00732 ;; NA19240) echo NA19238,NA19239 ;;
    HG00514) echo HG00512,HG00513 ;; NA12878) echo NA12891,NA12892 ;; *) echo "" ;;
  esac
}
timed() {   # timed <step> <locus> <group> <arm> <command...>
  local s; s=$(date +%s)
  "${@:5}"
  printf '%s\t%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "$4" $(( $(date +%s) - s )) >> "$TIMES"
}
has_tx() { grep -q '"transcriptomes"' "$BG/manifest.json"; }
free_gb() { df -Pk "$OUT" | awk 'NR==2 {print int($4 / 1048576)}'; }
check_disk() { [[ $(free_gb) -ge 5 ]] || { echo "rna_simulate.sh: under 5 GB free, stopping" >&2; exit 1; }; }

want() { [[ -z ${LOCI:-} ]] || [[ " $LOCI " == *" $1 "* ]]; }
cat_of() { local v="CAT_$1"; echo "${!v:-}"; }

# held-out groups per locus, written once so a resumed run keeps them
groups_of() {   # groups_of <locus> <gfa>: lines "g<i> <sample,...> <sample+relatives,...>"
  local f="$OUT/$1/groups.txt"
  if [[ ! -s $f ]]; then
    local have; have=$(zcat < "$2" | awk '$1=="P" {split($2, a, "#"); if (a[1] !~ /^(GRCh38|grch38|CHM13|chm13|chm13v2)$/) print a[1]}' | sort -u)
    local picked=() s rel
    for s in $PREFERRED $have; do
      [[ ${#picked[@]} -ge $((NGROUPS * GROUP_SIZE)) ]] && break
      grep -qx "$s" <<< "$have" || continue
      [[ " ${picked[*]:-} " == *" $s "* ]] && continue
      # a relative of a sample already taken is not a sample of its own
      rel=0; for x in "${picked[@]:-}"; do if [[ ",$(relatives "$x")," == *",$s,"* ]]; then rel=1; fi; done
      [[ $rel == 1 ]] && continue
      picked+=("$s")
    done
    local i=0 g out
    : > "$f"
    while [[ $((i * GROUP_SIZE)) -lt ${#picked[@]} ]]; do
      g=("${picked[@]:$((i * GROUP_SIZE)):$GROUP_SIZE}")
      out=$(IFS=,; echo "${g[*]}")
      for s in "${g[@]}"; do rel=$(relatives "$s"); if [[ -n $rel ]]; then out="$out,$rel"; fi; done
      echo "g$((i + 1)) $(IFS=,; echo "${g[*]}") $out" >> "$f"
      i=$((i + 1))
    done
  fi
  cat "$f"
}

index() {   # index <locus> <dir> <arm> [--transcripts ...]
  local L=$1 D=$2 arm=$3; shift 3
  [[ -f $D/idx_$arm/$L.plx ]] && return
  mkdir -p "$D/idx_$arm"
  timed index "$L" "$(basename "$D")" "$arm" nice -n 15 "$P" recruit --locus "$D/$L.gfa" --background "$BG" -r grch38 \
    "$@" --index-only --save-index "$D/idx_$arm" -t "$T" -q
}

arms_of() { local a="a b c"; [[ -n $(cat_of "$1") ]] && a="$a d"; echo $a; }

# 1. graphs, truth, projections, the indexes of arms a and b
for entry in "${LOCUS_TABLE[@]}"; do
  read -r L gfa ref regex <<< "$entry"
  want "$L" || continue
  gfa="$RD/$gfa"
  mkdir -p "$OUT/$L"
  while read -r g samples out_samples; do
    D="$OUT/$L/$g"; mkdir -p "$D"; check_disk
    echo "== $L $g: $samples (out of the graph: $out_samples) $(date +%H:%M:%S)"
    [[ -s $D/$L.gfa ]] || python3 "$ROOT/dev/recruit_sources.py" held-out --gfa "$gfa" --sample "$out_samples" --out "$D/$L.gfa"
    [[ -s $D/truth_graph.gfa ]] || $SIM only --gfa "$gfa" --sample "$samples" --keep "$ref" --out "$D/truth_graph.gfa"
    [[ -s $D/truth.hst.tsv ]] || timed project_truth "$L" "$g" . nice -n 15 "$P" transcripts -i "$D/truth_graph.gfa" -r "$ref" --gtf "$GTF" -o "$D/truth" -t "$T" -q
    needs_b=0; for arm in b c; do [[ -f $D/idx_$arm/$L.plx ]] || needs_b=1; done
    if [[ $needs_b == 1 && ! -s $D/b.hst.fa ]]; then
      rm -f "$D/b.hst.tsv.gz"
      timed project "$L" "$g" b nice -n 15 "$P" transcripts -i "$D/$L.gfa" -r "$ref" --gtf "$GTF" -o "$D/b" -t "$T" -q
    fi
    cat=$(cat_of "$L")
    if [[ -n $cat && ! -f $D/idx_d/$L.plx && ! -s $D/d.hst.fa ]]; then
      rm -f "$D/d.hst.tsv.gz"
      # the held-out samples' own annotations leave with them
      grep -v -E "^($(echo "$out_samples" | tr ',' '|'))#" "$cat" > "$D/cat.tsv"
      timed project "$L" "$g" d nice -n 15 "$P" transcripts -i "$D/$L.gfa" -r "$ref" --gtf "$GTF" --hap-annotations "$D/cat.tsv" \
        ${CONTIG_ALIAS:+--contig-alias "$CONTIG_ALIAS"} -o "$D/d" -t "$T" -q
    fi
    if ! [[ -f $D/idx_a/$L.plx && -f $D/idx_b/$L.plx ]] && has_tx; then
      echo "rna_simulate.sh: $BG holds a transcriptome, so arm b of $L $g cannot be indexed; give a store without one" >&2; exit 1
    fi
    index "$L" "$D" a
    index "$L" "$D" b --transcripts "$L=$D/b.hst.fa"
  done < <(groups_of "$L" "$gfa")
done

# 2. the transcriptome joins the store
if ! has_tx; then
  echo "== transcriptome $(date +%H:%M:%S)"
  timed transcriptome . . . nice -n 15 "$P" background -o "$BG" --add --transcriptome ens="$TX" -t "$T" --memory 4G -q
fi

# 3. arms c and d; 4. reads, recruitment and scores
for entry in "${LOCUS_TABLE[@]}"; do
  read -r L gfa ref regex <<< "$entry"
  want "$L" || continue
  sets=()
  while read -r g samples out_samples; do
    D="$OUT/$L/$g"; check_disk
    index "$L" "$D" c --transcripts "$L=$D/b.hst.fa"
    if [[ -n $(cat_of "$L") ]]; then index "$L" "$D" d --transcripts "$L=$D/d.hst.fa"; fi
    rm -f "$D"/b.hst.fa "$D"/d.hst.fa   # large, and remade in seconds if wanted
    for f in "$D"/b.hst.tsv "$D"/d.hst.tsv; do if [[ -f $f ]]; then gzip -1f "$f"; fi; done
    interval=$(awk -F'\t' '$1=="locus" && $2=="interval" {print $3}' "$D/idx_a/$L.index.tsv")
    if [[ ! -s $D/sim/truth.tsv ]]; then
      seed=$(( $(cksum <<< "$L$g" | cut -d' ' -f1) % 100000 ))
      timed simulate "$L" "$g" . nice -n 15 $SIM simulate --gfa "$D/truth_graph.gfa" --hst "$D/truth" --samples "$samples" --tx "$TX" \
        --interval "$interval" ${CHROMS:+--chroms "$CHROMS"} --paralog-regex "$regex" --out "$D/sim" --seed "$seed" --prefix "${g}_" \
        --n-short "${N_SHORT:-3000}" --n-long "${N_LONG:-1000}" --n-ont "${N_ONT:-1000}" \
        --n-paralog "${N_PARALOG:-1000}" --n-pseudo "${N_PSEUDO:-500}" --n-tx "${N_TX:-3000}" \
        --pbsim "$PBSIM" --pbsim-model "$PBSIM_MODEL"
    fi
    mkdir -p "$D/rec"
    for arm in $(arms_of "$L"); do
      for cls in short kinnex ont; do
        case $cls in
          short) in=(-1 "$D/sim/short_1.fq.gz" -2 "$D/sim/short_2.fq.gz"); pre=illumina ;;
          kinnex) in=(-i "$D/sim/kinnex.fq.gz"); pre=hifi ;;
          ont) in=(-i "$D/sim/ont.fq.gz"); pre=ont ;;
        esac
        out="$D/rec/$arm.$cls.tsv"
        [[ -s $out ]] || timed recruit "$L" "$g" "$arm" nice -n 15 "$P" recruit --locus "$D/idx_$arm/$L.plx" "${in[@]}" --preset $pre -o "$out" -t "$T" -q
        sets+=("${arm}:${cls}:${out}")
      done
    done
  done < <(groups_of "$L" "$RD/$gfa")
  # every group's reads scored together (read names carry the group)
  awk 'FNR > 1 || NR == 1' "$OUT/$L"/g*/sim/truth.tsv > "$OUT/$L/truth.tsv"
  for arm in $(arms_of "$L"); do
    for cls in short kinnex ont; do
      awk 'FNR > 1 || NR == 1' "$OUT/$L"/g*/rec/$arm.$cls.tsv > "$OUT/$L/$arm.$cls.tsv"
    done
  done
  sets=()
  for arm in $(arms_of "$L"); do for cls in short kinnex ont; do sets+=("${arm}:${cls}:$OUT/$L/$arm.$cls.tsv"); done; done
  echo "== $L scores"
  timed score "$L" . . python3 "$ROOT/dev/rna_score.py" --locus "$L" --truth "$OUT/$L/truth.tsv" --out "$OUT/$L/scores.tsv" "${sets[@]}" \
    | tee "$OUT/$L/scores.txt" | column -t
done
echo "== done $(date +%H:%M:%S)"
du -sh "$OUT"/* 2>/dev/null | sed 's/^/   /'
