#!/usr/bin/env bash
# RNA recruitment benchmark. Per run (one locus, or loci indexed together), groups of held-out samples
# (and their relatives) are removed from the graphs; reads come from their haplotype-specific transcripts
# (projected onto their own walks), from molecules the annotation lacks, from their own annotations' (CAT)
# transcripts where given, and from negatives elsewhere in the transcriptome (dev/rna_sim.py). The same
# reads are recruited in each arm:
#   a  genomic walks only
#   b  + the held-out graphs' transcripts (panvar transcripts)
#   c  + a transcriptome segment in the background
#   d  as c, with the haplotypes' own annotations merged (runs whose loci are given CAT_<locus>)
# Arms a and b are indexed before the transcriptome joins the store, c and d after; reads are recruited
# with each arm's .plx as built. A run's loci are indexed in one call, so each is the others' companion,
# and each locus's reads are recruited against all of them. Steps already done are kept, so a stopped run
# resumes; what a finished group no longer needs is removed.
#   rna_simulate.sh <out dir> <store> <transcriptome.fa.gz> [pbsim3 dir]
# The store holds the genomes only (no transcriptome yet, unless every a/b index exists); the
# transcriptome comes from `rna_sim.py transcriptome` and is added to the store as `ens`.
# Environment:
#   RUNS      runs to do, from the table below (default: all)
#   NGROUPS   held-out groups per run (4), GROUP_SIZE samples per group (2), SEEDS expression draws (3)
#   N_SHORT N_LONG N_ONT   locus fragments per group, seed and class (3000, 1000, 1000)
#   N_PARALOG N_PSEUDO N_TX  negative fragments per group, seed and class (1000, 500, 3000)
#   CHROMS    negatives only from these chromosomes, for a background of some chromosomes
#   CAT_<locus>=<tsv>     haplotypes' annotations (transcripts --hap-annotations) for arm d
#   CATTRUTH_<locus>=<dir>  the held-out haplotypes' own annotations, one GFF3 each, as a second truth
#   CONTIG_ALIAS  --contig-alias for the arm d annotations
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
SCORE="python3 $ROOT/dev/rna_score.py"
RD="$ROOT/tests/real_data"
GTF="$RD/Homo_sapiens.GRCh38.116.loci.gtf.gz"
T="${THREADS:-2}"
NGROUPS="${NGROUPS:-4}"; GROUP_SIZE="${GROUP_SIZE:-2}"; SEEDS="${SEEDS:-3}"
TIMES="$OUT/times.tsv"
[[ -f $TIMES ]] || printf 'step\trun\tgroup\tarm\tseconds\n' > "$TIMES"

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
# runs: name and loci indexed together
RUN_TABLE=("gstm gstm" "acot acot" "c4 c4" "cyp2d cyp2d" "lpa lpa" "hla hla_a,hla_b")
# graphs of related genes outside a run: their held-out haplotypes' transcripts are paralog negatives
HAPNEG_hla="hla/hla_c.gfa.gz hla/hla_e.gfa.gz hla/hla_f.gfa.gz hla/hla_g.gfa.gz hla/hla_j.gfa.gz hla/hla_l.gfa.gz hla/hla_n.gfa.gz hla/hla_s.gfa.gz"
# held-out samples are taken in this order, then the graph's others; relatives leave the graph with them
PREFERRED="HG002 NA24385 HG00733 NA19240 HG005 HG02257 HG01891 HG03516 NA18906 HG00514 NA12878"
relatives() {
  case $1 in
    HG002) echo HG003,HG004 ;; HG00733) echo HG00731,HG00732 ;; NA19240) echo NA19238,NA19239 ;;
    HG00514) echo HG00512,HG00513 ;; NA12878) echo NA12891,NA12892 ;; *) echo "" ;;
  esac
}
timed() {   # timed <step> <run> <group> <arm> <command...>
  local s; s=$(date +%s)
  "${@:5}"
  printf '%s\t%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "$4" $(( $(date +%s) - s )) >> "$TIMES"
}
has_tx() { grep -q '"transcriptomes"' "$BG/manifest.json"; }
free_gb() { df -Pk "$OUT" | awk 'NR==2 {print int($4 / 1048576)}'; }
check_disk() { [[ $(free_gb) -ge 5 ]] || { echo "rna_simulate.sh: under 5 GB free, stopping" >&2; exit 1; }; }
want() { [[ -z ${RUNS:-} ]] || [[ " $RUNS " == *" $1 "* ]]; }
var() { local v="$1"; echo "${!v:-}"; }
locus_field() {   # locus_field <locus> <1-based field of the table>
  local e; for e in "${LOCUS_TABLE[@]}"; do read -r -a f <<< "$e"; [[ ${f[0]} == "$1" ]] && { echo "${f[$(( $2 - 1 ))]}"; return; }; done
  echo "rna_simulate.sh: no locus $1" >&2; exit 1
}

# held-out groups of a run, written once so a resumed run keeps them
groups_of() {   # groups_of <run> <gfa>: lines "g<i> <sample,...> <sample+relatives,...>"
  local f="$OUT/$1/groups.txt"
  if [[ ! -s $f ]]; then
    local have; have=$(zcat < "$2" | awk '$1=="P" {split($2, a, "#"); if (a[1] !~ /^(GRCh38|grch38|CHM13|chm13|chm13v2)$/) print a[1]}' | sort -u)
    local picked=() s rel x
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

project() {   # project <run> <dir> <locus> <arm b|d>: the held-out graph's transcripts, <dir>/<locus>.<arm>.hst.*
  local R=$1 D=$2 L=$3 arm=$4 ref cat
  [[ -s $D/$L.$arm.hst.fa ]] && return
  ref=$(locus_field "$L" 3)
  if [[ $arm == b ]]; then
    timed project "$R" "$(basename "$D")" b nice -n 15 "$P" transcripts -i "$D/$L.gfa" -r "$ref" --gtf "$GTF" -o "$D/$L.b" -t "$T" -q
  else
    cat=$(var "CAT_$L")
    # the held-out samples' own annotations leave with them
    grep -v -E "^($(tr ',' '|' <<< "$OUT_SAMPLES"))#" "$cat" > "$D/$L.cat.tsv"
    timed project "$R" "$(basename "$D")" d nice -n 15 "$P" transcripts -i "$D/$L.gfa" -r "$ref" --gtf "$GTF" \
      --hap-annotations "$D/$L.cat.tsv" ${CONTIG_ALIAS:+--contig-alias "$CONTIG_ALIAS"} -o "$D/$L.d" -t "$T" -q
  fi
}

index() {   # index <run> <dir> <arm> <loci,...>: one call for the run's loci, transcripts from b or d
  local R=$1 D=$2 arm=$3 loci=$4 L args=() src
  local done_all=1
  for L in ${loci//,/ }; do [[ -f $D/idx_$arm/$L.plx ]] || done_all=0; done
  [[ $done_all == 1 ]] && return
  mkdir -p "$D/idx_$arm"
  src=$arm; [[ $arm == c ]] && src=b
  for L in ${loci//,/ }; do
    args+=(--locus "$D/$L.gfa")
    if [[ $arm != a ]]; then
      if [[ $src == d && -z $(var "CAT_$L") ]]; then src=b; fi
      project "$R" "$D" "$L" "$src"
      args+=(--transcripts "$L=$D/$L.$src.hst.fa")
    fi
  done
  timed index "$R" "$(basename "$D")" "$arm" nice -n 15 "$P" recruit "${args[@]}" --background "$BG" -r grch38 \
    --index-only --save-index "$D/idx_$arm" -t "$T" -q
  # the FASTAs are large and remade in seconds; the tables are kept for what each source adds
  for L in ${loci//,/ }; do
    rm -f "$D/$L.$src.hst.fa"
    if [[ -f $D/$L.$src.hst.tsv ]]; then gzip -1f "$D/$L.$src.hst.tsv"; fi
  done
}

has_cat() { local L; for L in ${1//,/ }; do [[ -n $(var "CAT_$L") ]] && return 0; done; return 1; }
arms_of() { if has_cat "$1"; then echo a b c d; else echo a b c; fi; }

# 1. graphs, truth, the indexes of arms a and b
for entry in "${RUN_TABLE[@]}"; do
  read -r R loci <<< "$entry"
  want "$R" || continue
  first=${loci%%,*}
  mkdir -p "$OUT/$R"
  while read -r g samples OUT_SAMPLES; do
    D="$OUT/$R/$g"; mkdir -p "$D"; check_disk
    [[ -f $D/done ]] && continue
    echo "== $R $g: $samples (out of the graphs: $OUT_SAMPLES) $(date +%H:%M:%S)"
    for L in ${loci//,/ }; do
      gfa="$RD/$(locus_field "$L" 2)"; ref=$(locus_field "$L" 3)
      [[ -s $D/$L.gfa ]] || python3 "$ROOT/dev/recruit_sources.py" held-out --gfa "$gfa" --sample "$OUT_SAMPLES" --out "$D/$L.gfa"
      [[ -s $D/$L.truth_graph.gfa ]] || $SIM only --gfa "$gfa" --sample "$samples" --keep "$ref" --out "$D/$L.truth_graph.gfa"
      [[ -s $D/$L.truth.hst.fa ]] || timed project_truth "$R" "$g" . nice -n 15 "$P" transcripts -i "$D/$L.truth_graph.gfa" -r "$ref" --gtf "$GTF" -o "$D/$L.truth" -t "$T" -q
    done
    ab=1; for L in ${loci//,/ }; do [[ -f $D/idx_a/$L.plx && -f $D/idx_b/$L.plx ]] || ab=0; done
    if [[ $ab == 0 ]] && has_tx; then
      echo "rna_simulate.sh: $BG holds a transcriptome, so arm b of $R $g cannot be indexed; give a store without one" >&2; exit 1
    fi
    index "$R" "$D" a "$loci"
    index "$R" "$D" b "$loci"
  done < <(groups_of "$R" "$RD/$(locus_field "$first" 2)")
done

# 2. the transcriptome joins the store
if ! has_tx; then
  echo "== transcriptome $(date +%H:%M:%S)"
  timed transcriptome . . . nice -n 15 "$P" background -o "$BG" --add --transcriptome ens="$TX" -t "$T" --memory 4G -q
fi

# 3. arms c and d, reads, recruitment; what the group no longer needs is removed
for entry in "${RUN_TABLE[@]}"; do
  read -r R loci <<< "$entry"
  want "$R" || continue
  first=${loci%%,*}
  while read -r g samples OUT_SAMPLES; do
    D="$OUT/$R/$g"; check_disk
    [[ -f $D/done ]] && continue
    index "$R" "$D" c "$loci"
    if has_cat "$loci"; then index "$R" "$D" d "$loci"; fi
    # other graphs' transcripts of these samples, and the companions' own, as haplotype paralog negatives
    hapneg=$(var "HAPNEG_$R")
    pick='BEGIN {split(s, a, ","); for (i in a) keep[a[i]] = 1} /^>/ {split(substr($1, 2), f, "|"); split(f[3], h, "#"); p = (h[1] in keep)} p'
    if [[ -n $hapneg && ! -s $D/hapneg.fa ]]; then
      for x in $hapneg; do
        n=$(basename "$x" .gfa.gz)
        $SIM only --gfa "$RD/$x" --sample "$samples" --keep GRCh38 --out "$D/neg_$n.gfa" 2> /dev/null
        timed project_negatives "$R" "$g" . nice -n 15 "$P" transcripts -i "$D/neg_$n.gfa" -r GRCh38 --gtf "$GTF" -o "$D/neg_$n" -t "$T" -q
        awk -v s="$samples" "$pick" "$D/neg_$n.hst.fa" >> "$D/hapneg.fa.tmp"
        rm -f "$D/neg_$n".*
      done
      mv "$D/hapneg.fa.tmp" "$D/hapneg.fa"
    fi
    for L in ${loci//,/ }; do
      interval=$(awk -F'\t' '$1=="locus" && $2=="interval" {print $3}' "$D/idx_a/$L.index.tsv")
      negs=()
      if [[ -s $D/hapneg.fa ]]; then negs+=(--hap-negatives "$D/hapneg.fa"); fi
      for C in ${loci//,/ }; do
        if [[ $C != "$L" ]]; then awk -v s="$samples" "$pick" "$D/$C.truth.hst.fa" > "$D/companion_$C.fa"; negs+=(--hap-negatives "$D/companion_$C.fa"); fi
      done
      cattruth=$(var "CATTRUTH_$L")
      if [[ ! -s $D/sim_$L/s$SEEDS/truth.tsv ]]; then
        seed=$(( $(cksum <<< "$R$L$g" | cut -d' ' -f1) % 100000 ))
        regex=$(locus_field "$L" 4)
        timed simulate "$R" "$g" "$L" nice -n 15 $SIM simulate --gfa "$D/$L.truth_graph.gfa" --hst "$D/$L.truth" --samples "$samples" \
          --tx "$TX" --interval "$interval" ${CHROMS:+--chroms "$CHROMS"} --paralog-regex "$regex" --out "$D/sim_$L" \
          --seed "$seed" --seeds "$SEEDS" --prefix "${g}_" ${cattruth:+--cat-dir "$cattruth"} ${negs[@]+"${negs[@]}"} \
          --n-short "${N_SHORT:-3000}" --n-long "${N_LONG:-1000}" --n-ont "${N_ONT:-1000}" \
          --n-paralog "${N_PARALOG:-1000}" --n-pseudo "${N_PSEUDO:-500}" --n-tx "${N_TX:-3000}" \
          --pbsim "$PBSIM" --pbsim-model "$PBSIM_MODEL"
      fi
      rm -f "$D"/companion_*.fa
      mkdir -p "$D/rec"
      for arm in $(arms_of "$loci"); do
        idx=(); for C in ${loci//,/ }; do idx+=(--locus "$D/idx_$arm/$C.plx"); done
        for k in $(seq 1 "$SEEDS"); do
          S="$D/sim_$L/s$k"
          for cls in short kinnex ont; do
            case $cls in
              short) in=(-1 "$S/short_1.fq.gz" -2 "$S/short_2.fq.gz"); pre=illumina ;;
              kinnex) in=(-i "$S/kinnex.fq.gz"); pre=hifi ;;
              ont) in=(-i "$S/ont.fq.gz"); pre=ont ;;
            esac
            out="$D/rec/$arm.$L.$cls.s$k.tsv"
            if [[ ! -s $out.gz ]]; then
              timed recruit "$R" "$g" "$arm" nice -n 15 "$P" recruit "${idx[@]}" "${in[@]}" --preset $pre -o "$out" -t "$T" -q
              cut -f1-4 "$out" | gzip -1 > "$out.gz"; rm -f "$out" "${out%.tsv}.summary.tsv"
            fi
          done
        done
      done
    done
    # the group is done: keep the truth, the tables and the index reports
    for L in ${loci//,/ }; do
      for k in $(seq 1 "$SEEDS"); do rm -f "$D/sim_$L/s$k"/*.fq.gz; done
      rm -f "$D/$L.gfa" "$D/$L.truth_graph.gfa" "$D/$L.truth.hst.fa" "$D/$L.cat.tsv"
      if [[ -f $D/$L.truth.hst.tsv ]]; then gzip -1f "$D/$L.truth.hst.tsv"; fi
    done
    rm -f "$D"/idx_*/*.plx "$D/hapneg.fa"
    touch "$D/done"
    echo "== $R $g done $(date +%H:%M:%S), $(free_gb) GB free"
  done < <(groups_of "$R" "$RD/$(locus_field "$first" 2)")

  # 4. scores per locus and seed, every group together (read names carry the group), then over seeds
  for L in ${loci//,/ }; do
    seeds_out=()
    for k in $(seq 1 "$SEEDS"); do
      awk 'FNR > 1 || NR == 1' "$OUT/$R"/g*/sim_$L/s$k/truth.tsv > "$OUT/$R/$L.truth.s$k.tsv"
      sets=()
      for arm in $(arms_of "$loci"); do
        for cls in short kinnex ont; do
          t="$OUT/$R/$L.$arm.$cls.s$k.tsv"
          { echo -e "fragment\tlocus\tfate\treason"; for f in "$OUT/$R"/g*/rec/$arm.$L.$cls.s$k.tsv.gz; do gzip -dc "$f" | tail -n +2; done; } > "$t"
          sets+=("${arm}:${cls}:$t")
        done
      done
      timed score "$R" . "$L" $SCORE score --locus "$L" --truth "$OUT/$R/$L.truth.s$k.tsv" --out "$OUT/$R/$L.scores.s$k.tsv" "${sets[@]}" > "$OUT/$R/$L.scores.s$k.txt"
      for s in "${sets[@]}"; do rm -f "${s##*:}"; done
      gzip -1f "$OUT/$R/$L.truth.s$k.tsv"
      seeds_out+=("$OUT/$R/$L.scores.s$k.tsv")
    done
    $SCORE merge --out "$OUT/$R/$L.scores.tsv" "${seeds_out[@]}"
    echo "== $R $L scored"
  done
done
echo "== done $(date +%H:%M:%S), $(free_gb) GB free"
du -sh "$OUT"/* 2>/dev/null | sed 's/^/   /'
